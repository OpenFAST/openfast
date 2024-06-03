!**********************************************************************************************************************************
! The HydroDyn and HydroDyn_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. HydroDyns_Types will be auto-generated based on a description of the variables for the module.
!
! "HydroDyn" should be replaced with the name of your module. Example: HydroDyn
! "HydroDyn" (in HydroDyn_*) should be replaced with the module name or an abbreviation of it. Example: HD
!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
!
!    This file is part of HydroDyn.
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
MODULE HydroDyn

   USE HydroDyn_Types   
   USE NWTC_Library
   use Morison
   USE WAMIT
   USE WAMIT2
   USE SeaState
   USE HydroDyn_Input
   USE HydroDyn_Output

#ifdef USE_FIT
   USE FIT_MODULES
   USE FIT_Types
#endif      
   IMPLICIT NONE
   
   PRIVATE

  
   TYPE(ProgDesc), PARAMETER            :: HydroDyn_ProgDesc = ProgDesc( 'HydroDyn', '', '' )

    
   
   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: HydroDyn_Init                           ! Initialization routine
   PUBLIC :: HydroDyn_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: HydroDyn_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: HydroDyn_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: HydroDyn_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: HydroDyn_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   !PUBLIC :: HydroDyn_UpdateDiscState                ! Tight coupling routine for updating discrete states
      
   PUBLIC :: HD_JacobianPInput                 ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the inputs(u)
   PUBLIC :: HD_JacobianPContState             ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the continuous
                                               !   states(x)
   PUBLIC :: HD_JacobianPDiscState             ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the discrete
                                               !   states(xd)
   PUBLIC :: HD_JacobianPConstrState           ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the constraint
                                               !   states(z)
   PUBLIC :: HD_GetOP                          !< Routine to pack the operating point values (for linearization) into arrays
   
   CONTAINS
   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE HydroDyn_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(HydroDyn_InitInputType),       INTENT(INOUT)  :: InitInp     !< Input data for initialization routine. [INOUT because of a move_alloc() statement]
      TYPE(HydroDyn_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
      TYPE(HydroDyn_ParameterType),       INTENT(  OUT)  :: p           !< Parameters      
      TYPE(HydroDyn_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
      TYPE(HydroDyn_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
      TYPE(HydroDyn_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
      TYPE(HydroDyn_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states            
      TYPE(HydroDyn_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated; 
                                                                        !!   only the output mesh is initialized)
      TYPE(HydroDyn_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables           
      REAL(DbKi),                         INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                        !!   (1) HydroDyn_UpdateStates() is called in loose coupling &
                                                                        !!   (2) HydroDyn_UpdateDiscState() is called in tight coupling.
                                                                        !!   Input is the suggested time from the glue code; 
                                                                        !!   Output is the actual coupling interval that will be used 
                                                                        !!   by the glue code.
      TYPE(HydroDyn_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      
         ! Local variables
         
      CHARACTER(1024)                        :: SummaryName                         ! name of the HydroDyn summary file   
      TYPE(HydroDyn_InputFile)               :: InputFileData                       !< Data from input file
      TYPE(FileInfoType)                     :: InFileInfo                          !< The derived type for holding the full input file for parsing -- we may pass this in the future
!      LOGICAL                                :: hasWAMITOuts                        ! Are there any WAMIT-related outputs
!      LOGICAL                                :: hasMorisonOuts                      ! Are there any Morison-related outputs
!      INTEGER                                :: numHydroOuts                        ! total number of WAMIT and Morison outputs
      INTEGER                                :: I, J, k, iBody                                ! Generic counters
         ! These are dummy variables to satisfy the framework, but are not used 
         
 
#ifdef USE_FIT
         ! FIT - related data
      TYPE(FIT_InitInputType)                :: FITInitData 
      TYPE(FIT_InputType)                    :: FIT_u                             ! FIT module initial guess for the input; the input mesh is not defined because it is not used by the waves module
      TYPE(FIT_ParameterType)                :: FIT_p                             ! FIT module parameters
      TYPE(FIT_ContinuousStateType)          :: FIT_x                             ! FIT module initial continuous states
      TYPE(FIT_DiscreteStateType)            :: FIT_xd                            ! FIT module discrete states
      TYPE(FIT_ConstraintStateType)          :: FIT_z                             ! FIT module initial guess of the constraint states
      TYPE(FIT_OtherStateType)               :: FIT_OtherState                    ! FIT module other/optimization states 
      TYPE(FIT_OutputType)                   :: FIT_y                             ! FIT module outputs  
      TYPE(FIT_InitOutputType)               :: FIT_InitOut                       ! Initialization Outputs from the FIT module initialization
#endif

   
         ! WAMIT Mesh
      real(R8Ki)                             :: theta(3), orientation(3,3)        
      
      INTEGER(IntKi)                         :: ErrStat2                            ! local error status
      CHARACTER(ErrMsgLen)                   :: ErrMsg2                             ! local error message
      CHARACTER(*), PARAMETER                :: RoutineName = 'HydroDyn_Init'
   

         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      p%UnOutFile = -1
      
      p%WaveField    =>  InitInp%WaveField

      
         ! Initialize the NWTC Subroutine Library
         
      CALL NWTC_Init(  )
     
        
         ! Display the module information

      CALL DispNVD( HydroDyn_ProgDesc )        
      

      IF ( InitInp%UseInputFile ) THEN
         CALL ProcessComFile( InitInp%InputFile, InFileInfo, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      ELSE
         CALL NWTC_Library_CopyFileInfoType( InitInp%PassedFileData, InFileInfo, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF          
      ENDIF

      ! For diagnostic purposes, the following can be used to display the contents
      ! of the InFileInfo data structure.
      ! call Print_FileInfo_Struct( CU, InFileInfo ) ! CU is the screen -- different number on different systems.


      ! Parse all HydroDyn-related input and populate the InputFileData structure 
      CALL HydroDyn_ParseInput( InitInp%InputFile, InitInp%OutRootName, InFileInfo, InputFileData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
      InputFileData%Morison%WaveField => InitInp%WaveField
      InputFileData%WAMIT%WaveField   => InitInp%WaveField
      InputFileData%WAMIT2%WaveField  => InitInp%WaveField

         
      
         ! Verify all the necessary initialization data. Do this at the HydroDynInput module-level 
         !   because the HydroDynInput module is also responsible for parsing all this 
         !   initialization data from a file

      CALL HydroDynInput_ProcessInitData( InitInp, Interval, InputFileData, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
      
      
        ! Since the Convolution Radiation module is currently the only module which requires knowledge of the time step size, 
        !  we will set Hydrodyn's time step to be that of the Convolution radiation module if it is being used.  Otherwise, we
        !  will set it to be equal to the glue-codes
      IF ((InputFileData%PotMod == 1) .AND. (InputFileData%WAMIT%RdtnMod == 1) ) THEN
         
         
         p%DT = InputFileData%WAMIT%Conv_Rdtn%RdtnDT
 
#ifdef USE_FIT
      ELSE IF (Initlocal%PotMod == 2) THEN
         ! This is the FIT potential flow model and the time step needs to be >= the driver timestep, and and integer multiple if larger
         ! We example WaveDT for this timestep size because FIT is tied to WaveDT
         IF ( ( .NOT. EqualRealNos(mod(real(Initlocal%WaveDT,ReKi), real(Interval,ReKi)) , 0.0_ReKi) ) .OR. Initlocal%WaveDT <= 0.0_DbKi ) THEn
            CALL SetErrStat(ErrID_Fatal,'The value of WaveDT is not greater than zero and an integer multiple of the glue code timestep.',ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF 
         ELSE
            p%DT = Interval  !  If the above check is ok, then we can still set the module DT to the glue-code dt
            
         END IF
#endif

      ELSE
         
         p%DT = Interval
      END IF  
      
         ! Open a summary of the HydroDyn Initialization. Note: OutRootName must be set by the caller because there may not be an input file to obtain this rootname from.
         
      IF ( InputFileData%HDSum ) THEN 
         
         SummaryName = TRIM(InitInp%OutRootName)//'.sum'
         CALL HDOut_OpenSum( InputFileData%UnSum, SummaryName, HydroDyn_ProgDesc, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      
      ELSE
         
         InputFileData%UnSum = -1
         
      END IF
      
         ! Copy Additional preload, stiffness, and damping to the parameters
      p%AddF0        = InputFileData%AddF0
      p%AddCLin      = InputFileData%AddCLin
      p%AddBLin      = InputFileData%AddBLin
      p%AddBQuad     = InputFileData%AddBQuad
      
      
         ! Set summary unit number in Morison initialization input data
      InputFileData%Morison%UnSum         = InputFileData%UnSum
    
         ! Were visualization meshes requested?
      p%VisMeshes = InitInp%VisMeshes


         ! Now call each sub-module's *_Init subroutine
         ! to fully initialize each sub-module based on the necessary initialization data
   
            ! Is there a WAMIT body? 
         
         IF ( InputFileData%PotMod == 1 ) THEN
            InputFileData%WAMIT%WaveField  => InitInp%WaveField
            
            p%nWAMITObj              = InputFileData%nWAMITObj      ! All the data for the various WAMIT bodies are stored in a single WAMIT file
            p%vecMultiplier          = InputFileData%vecMultiplier  ! Multiply all vectors and matrices row/column lengths by NBody
            InputFileData%WAMIT%NBodyMod = InputFileData%NBodyMod
            InputFileData%WAMIT%Gravity  = InitInp%Gravity
            p%NBody                  = InputFileData%NBody
            p%NBodyMod               = InputFileData%NBodyMod
            call AllocAry( m%F_PtfmAdd, 6*InputFileData%NBody, "m%F_PtfmAdd", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            call AllocAry( m%F_Waves  , 6*InputFileData%NBody, "m%F_Waves"  , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
               ! Determine how many WAMIT modules we need based on NBody and NBodyMod
            if (p%NBodyMod == 1) then
               InputFileData%WAMIT%NBody    = InputFileData%NBody   ! The WAMIT object will contain all NBody WAMIT bodies
               
                  ! Allocate WAMIT InitInp arrays based on NBodyMod and copy the inputfile data into the WAMIT init data (entire arrays' worth for NBodyMod=1
               call AllocAry( InputFileData%WAMIT%PtfmVol0    , InputFileData%NBody, "PtfmVol0"    , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefxt   , InputFileData%NBody, "PtfmRefxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefyt   , InputFileData%NBody, "PtfmRefyt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefzt   , InputFileData%NBody, "PtfmRefzt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefztRot, InputFileData%NBody, "PtfmRefztRot", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmCOBxt   , InputFileData%NBody, "PtfmCOBxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmCOByt   , InputFileData%NBody, "PtfmCOByt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               allocate( p%WAMIT(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( x%WAMIT(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array x%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( xd%WAMIT(        1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array xd%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( OtherState%WAMIT(1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array OtherState%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( y%WAMIT(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array y%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( m%WAMIT(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( m%u_WAMIT(       1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%u_WAMIT.', ErrStat, ErrMsg, RoutineName )
                     
               InputFileData%WAMIT%PtfmVol0       = InputFileData%PtfmVol0    
               InputFileData%WAMIT%WAMITULEN      = InputFileData%WAMITULEN(1)   
               InputFileData%WAMIT%PtfmRefxt      = InputFileData%PtfmRefxt   
               InputFileData%WAMIT%PtfmRefyt      = InputFileData%PtfmRefyt   
               InputFileData%WAMIT%PtfmRefzt      = InputFileData%PtfmRefzt   
               InputFileData%WAMIT%PtfmRefztRot   = InputFileData%PtfmRefztRot
               InputFileData%WAMIT%PtfmCOBxt      = InputFileData%PtfmCOBxt   
               InputFileData%WAMIT%PtfmCOByt      = InputFileData%PtfmCOByt   
            else
               InputFileData%WAMIT%NBody    = 1    ! Each WAMIT object will only contain one of the NBody WAMIT bodies

                  ! Allocate WAMIT InitInp arrays based on NBodyMod and copy the inputfile data into the 1st WAMIT body init data for NBodyMod > 1
               call AllocAry( InputFileData%WAMIT%PtfmVol0    , 1, "PtfmVol0"    , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefxt   , 1, "PtfmRefxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefyt   , 1, "PtfmRefyt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefzt   , 1, "PtfmRefzt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmRefztRot, 1, "PtfmRefztRot", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmCOBxt   , 1, "PtfmCOBxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               call AllocAry( InputFileData%WAMIT%PtfmCOByt   , 1, "PtfmCOByt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               allocate( p%WAMIT(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( x%WAMIT(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array x%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( xd%WAMIT(        InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array xd%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( OtherState%WAMIT(InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array OtherState%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( y%WAMIT(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array y%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( m%WAMIT(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT.', ErrStat, ErrMsg, RoutineName )
               allocate( m%u_WAMIT(       InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%u_WAMIT.', ErrStat, ErrMsg, RoutineName )
               InputFileData%WAMIT%PtfmVol0    (1)  = InputFileData%PtfmVol0    (1)
               InputFileData%WAMIT%WAMITULEN        = InputFileData%WAMITULEN   (1)
               InputFileData%WAMIT%PtfmRefxt   (1)  = InputFileData%PtfmRefxt   (1)
               InputFileData%WAMIT%PtfmRefyt   (1)  = InputFileData%PtfmRefyt   (1)
               InputFileData%WAMIT%PtfmRefzt   (1)  = InputFileData%PtfmRefzt   (1)
               InputFileData%WAMIT%PtfmRefztRot(1)  = InputFileData%PtfmRefztRot(1)
               InputFileData%WAMIT%PtfmCOBxt   (1)  = InputFileData%PtfmCOBxt   (1)
               InputFileData%WAMIT%PtfmCOByt   (1)  = InputFileData%PtfmCOByt   (1)
  
            end if
            
            if ( ErrStat >= AbortErrLev ) then
               call CleanUp()
               return
            end if
            

            CALL WAMIT_Init(InputFileData%WAMIT, m%u_WAMIT(1), p%WAMIT(1), x%WAMIT(1), xd%WAMIT(1), z%WAMIT, OtherState%WAMIT(1), &
                                    y%WAMIT(1), m%WAMIT(1), Interval, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
                CALL CleanUp()
                RETURN
            END IF
               
               
                  ! For NBodyMod > 1 and NBody > 1, set the body info and init the WAMIT body
            do i = 2, p%nWAMITObj
                !-----------------------------------------
                ! Initialize the WAMIT Calculations 
                !-----------------------------------------
                InputFileData%WAMIT%WAMITFile        = InputFileData%PotFile     (i)
                InputFileData%WAMIT%PtfmVol0    (1)  = InputFileData%PtfmVol0    (i)
                InputFileData%WAMIT%WAMITULEN        = InputFileData%WAMITULEN   (i)
                InputFileData%WAMIT%PtfmRefxt   (1)  = InputFileData%PtfmRefxt   (i)
                InputFileData%WAMIT%PtfmRefyt   (1)  = InputFileData%PtfmRefyt   (i)
                InputFileData%WAMIT%PtfmRefzt   (1)  = InputFileData%PtfmRefzt   (i)
                InputFileData%WAMIT%PtfmRefztRot(1)  = InputFileData%PtfmRefztRot(i)
                InputFileData%WAMIT%PtfmCOBxt   (1)  = InputFileData%PtfmCOBxt   (i)
                InputFileData%WAMIT%PtfmCOByt   (1)  = InputFileData%PtfmCOByt   (i)
 
                CALL WAMIT_Init(InputFileData%WAMIT, m%u_WAMIT(i), p%WAMIT(i), x%WAMIT(i), xd%WAMIT(i), z%WAMIT, OtherState%WAMIT(i), &
                                        y%WAMIT(i), m%WAMIT(i), Interval, ErrStat2, ErrMsg2 )
                CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                IF ( ErrStat >= AbortErrLev ) THEN
                    CALL CleanUp()
                    RETURN
                END IF
            end do
               
               ! Generate Summary file information for WAMIT module
                   ! Compute the load contribution from hydrostatics:
            IF ( InputFileData%UnSum > 0 ) THEN
                do iBody = 1, InputFileData%NBody
                    WRITE( InputFileData%UnSum, '(A18,I5)')          'WAMIT Model - Body',iBody
                    WRITE( InputFileData%UnSum, '(A18)')             '------------------'
                    WRITE( InputFileData%UnSum, '(A42,2X,ES15.6)') 'Displaced volume (m^3)                 :', InputFileData%PtfmVol0(iBody)
                    WRITE( InputFileData%UnSum, '(A42,2X,ES15.6)') 'X-offset of the center of buoyancy (m) :', InputFileData%PtfmCOBxt(iBody)
                    WRITE( InputFileData%UnSum, '(A42,2X,ES15.6)') 'Y-offset of the center of buoyancy (m) :', InputFileData%PtfmCOByt(iBody)
                    WRITE( InputFileData%UnSum,  '(/)' ) 
                    WRITE( InputFileData%UnSum, '(A81)' ) 'Buoyancy loads from members modelled with WAMIT, summed about ( 0.0, 0.0, 0.0 )'
                    WRITE( InputFileData%UnSum, '(18x,6(2X,A20))' ) ' BuoyFxi ', ' BuoyFyi ', ' BuoyFzi ', ' BuoyMxi ', ' BuoyMyi ', ' BuoyMzi '
                    WRITE( InputFileData%UnSum, '(18x,6(2X,A20))' ) '   (N)   ', '   (N)   ', '   (N)   ', '  (N-m)  ', '  (N-m)  ', '  (N-m)  '
                    WRITE( InputFileData%UnSum, '(A18,6(2X,ES20.6))') '  External:       ',0.0,0.0,p%WaveField%RhoXg*InputFileData%PtfmVol0(iBody),p%WaveField%RhoXg*InputFileData%PtfmVol0(iBody)*InputFileData%PtfmCOByt(iBody), -p%WaveField%RhoXg*InputFileData%PtfmVol0(iBody)*InputFileData%PtfmCOBxt(iBody), 0.0   ! and the moment about Y due to the COB being offset from the WAMIT reference point
                end do
            END IF


               !-----------------------------------------
               ! Initialize the WAMIT2 Calculations
               !-----------------------------------------

               ! Only call the WAMIT2_Init if one of the flags is set for a calculation
            IF ( InputFileData%WAMIT2%MnDriftF .OR. InputFileData%WAMIT2%NewmanAppF .OR. InputFileData%WAMIT2%DiffQTFF .OR. InputFileData%WAMIT2%SumQTFF ) THEN

                  ! Flag required for indicating when to try using arrays that are allocated
               p%WAMIT2used   = .TRUE.

                  ! Copy Waves initialization output into the initialization input type for the WAMIT module
               InputFileData%WAMIT2%Gravity     = InitInp%Gravity

                  ! Set values for all NBodyMods
               InputFileData%WAMIT2%NBodyMod    = InputFileData%NBodyMod        ! There are restrictions in WAMIT2 on which files may be used for MnDriftF or NewmanAppF for BodyMod > 1
               InputFileData%WAMIT2%WAMITULEN   = InputFileData%WAMITULEN(1)

                  ! Determine how many WAMIT2 modules we need based on NBody and NBodyMod
               if (p%NBodyMod == 1) then
                  InputFileData%WAMIT2%NBody     = InputFileData%NBody    ! The WAMIT2 object will contain all NBody WAMIT2 bodies

                     ! Allocate WAMIT2 InitInp arrays based on NBodyMod and copy the inputfile data into the WAMIT2 init data (entire arrays' worth for NBodyMod=1
                  call AllocAry( InputFileData%WAMIT2%PtfmRefxt   , InputFileData%NBody, "PtfmRefxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefyt   , InputFileData%NBody, "PtfmRefyt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefzt   , InputFileData%NBody, "PtfmRefzt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefztRot, InputFileData%NBody, "PtfmRefztRot", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  allocate( p%WAMIT2(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( y%WAMIT2(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array y%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( m%WAMIT2(         1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  InputFileData%WAMIT2%PtfmRefxt     = InputFileData%PtfmRefxt
                  InputFileData%WAMIT2%PtfmRefyt     = InputFileData%PtfmRefyt
                  InputFileData%WAMIT2%PtfmRefzt     = InputFileData%PtfmRefzt
                  InputFileData%WAMIT2%PtfmRefztRot  = InputFileData%PtfmRefztRot

               else
                  InputFileData%WAMIT2%NBody  = 1_IntKi               ! The WAMIT2 object will contain all NBody WAMIT2 bodies

                     ! Allocate WAMIT2 InitInp arrays based on NBodyMod and copy the inputfile data into the 1st WAMIT body init data for NBodyMod > 1
                  call AllocAry( InputFileData%WAMIT2%PtfmRefxt   , 1, "PtfmRefxt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefyt   , 1, "PtfmRefyt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefzt   , 1, "PtfmRefzt"   , ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  call AllocAry( InputFileData%WAMIT2%PtfmRefztRot, 1, "PtfmRefztRot", ErrStat2, ErrMsg2 ); call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  allocate( p%WAMIT2(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( y%WAMIT2(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array y%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  allocate( m%WAMIT2(         InputFileData%NBody), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT2.', ErrStat, ErrMsg, RoutineName )
                  InputFileData%WAMIT2%PtfmRefxt   (1)  = InputFileData%PtfmRefxt   (1)
                  InputFileData%WAMIT2%PtfmRefyt   (1)  = InputFileData%PtfmRefyt   (1)
                  InputFileData%WAMIT2%PtfmRefzt   (1)  = InputFileData%PtfmRefzt   (1)
                  InputFileData%WAMIT2%PtfmRefztRot(1)  = InputFileData%PtfmRefztRot(1)

               endif

                if ( ErrStat >= AbortErrLev ) then
                   call CleanUp()
                   return
                end if

               CALL WAMIT2_Init(InputFileData%WAMIT2, p%WAMIT2(1), y%WAMIT2(1), m%WAMIT2(1), ErrStat2, ErrMsg2 )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL CleanUp()
                  RETURN
               END IF

                     ! For NBodyMod > 1 and NBody > 1, set the body info and init the WAMIT2 body
               do i = 2, p%nWAMITObj
                   InputFileData%WAMIT2%WAMITFile        = InputFileData%PotFile     (i)
                   InputFileData%WAMIT2%WAMITULEN        = InputFileData%WAMITULEN   (i)
                   InputFileData%WAMIT2%PtfmRefxt   (1)  = InputFileData%PtfmRefxt   (i)
                   InputFileData%WAMIT2%PtfmRefyt   (1)  = InputFileData%PtfmRefyt   (i)
                   InputFileData%WAMIT2%PtfmRefzt   (1)  = InputFileData%PtfmRefzt   (i)
                   InputFileData%WAMIT2%PtfmRefztRot(1)  = InputFileData%PtfmRefztRot(i)

                   CALL WAMIT2_Init(InputFileData%WAMIT2, p%WAMIT2(i), y%WAMIT2(i), m%WAMIT2(i), ErrStat2, ErrMsg2 )
                   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
                   IF ( ErrStat >= AbortErrLev ) THEN
                       CALL CleanUp()
                       RETURN
                   END IF
               end do


            ELSE
                  ! Flag used in output handling to indicate when to ignore WAMIT2 outputs.
               p%WAMIT2used   = .FALSE.
               allocate( p%WAMIT2(1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array p%WAMIT2.', ErrStat, ErrMsg, RoutineName )
               allocate( m%WAMIT2(1), stat = ErrStat2 ); if (ErrStat2 /=0) call SetErrStat( ErrID_Fatal, 'Failed to allocate array m%WAMIT2.', ErrStat, ErrMsg, RoutineName )
               p%WAMIT2%MnDriftF    = .FALSE.
               p%WAMIT2%NewmanAppF  = .FALSE.
               p%WAMIT2%DiffQTFF    = .FALSE.
               p%WAMIT2%SumQTFF     = .FALSE.
            ENDIF

#ifdef USE_FIT 
         ELSE IF ( InputFileData%PotMod == 2  ) THEN  ! FIT 
            ! Set up the Initialization data for FIT
               ! General
            FITInitData%InputFile      = InputFileData%PotFile
            FITInitData%Gravity        = InputFileData%Gravity
            FITInitData%Rho            = p%WaveField%WtrDens
            FITInitData%time_end       = InitInp%TMax
            FITInitData%dtime          = InitInp%WaveDT  ! Set the FIT module's timestep equal to the WaveDT timestep, this was checked earlier to make sure it is an integer muliple of the glue-code timestep!
               ! Waves
               ! Need to pre-process the incoming wave data to be compatible with FIT
            
            FITInitData%N_omega        = p%WaveField%NStepWave2
            FITInitData%Wave_angle     = p%WaveField%WaveDir
            
               ! allocate waves data arrays for FIT
            CALL AllocAry( FITInitData%Wave_amp, FITInitData%N_omega, "Wave_amp", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL AllocAry( FITInitData%Wave_omega, FITInitData%N_omega, "Wave_omega", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL AllocAry( FITInitData%Wave_number, FITInitData%N_omega, "Wave_number", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL AllocAry( FITInitData%Wave_phase, FITInitData%N_omega, "Wave_phase", ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL Cleanup()
                  RETURN
               END IF
               
               ! Populate wave arrays (Need to double chech this part. It doesn't look right!)
            Np = 2*(p%WaveField%WaveDOmega + 1)
            DO I = 1 , p%WaveField%NStepWave2
               
               dftreal        = p%WaveField%WaveElevC0( 1, ABS(I ) )
               dftimag        = p%WaveField%WaveElevC0( 2, ABS(I ) )*SIGN(1,I)
               FITInitData%Wave_amp   (I) = sqrt( dftreal**2 + dftimag**2 )  * 2.0 / Np
               FITInitData%Wave_omega (I) = I*p%WaveField%WaveDOmega
               FITInitData%Wave_number(I) = I*p%WaveField%WaveDOmega**2. / InputFileData%Gravity
               FITInitData%Wave_phase (I) = atan2( dftimag, dftreal ) 
              
            END DO         
         
  
              ! Output
            FITInitData%RootName       = trim(InputFileData%OutRootName)//'.FIT'
                              
      
            CALL FIT_Init(FITInitData, u%FIT, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, y%FIT, Interval, FIT_InitOut, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
#endif

         END IF
   
      ! Are there Morison elements?
      IF ( InputFileData%Morison%NMembers > 0 ) THEN
     
            ! Were visualization meshes requested?
         InputFileData%Morison%VisMeshes = p%VisMeshes
        
            ! Initialize the Morison Element Calculations 
         CALL Morison_Init(InputFileData%Morison, u%Morison, p%Morison, x%Morison, xd%Morison, z%Morison, OtherState%Morison, &
                               y%Morison, m%Morison, Interval, InitOut%Morison, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         
      END IF  ! Has Morison elements
    
!===============================================
      p%PotMod = InputFileData%Potmod      
      IF ( InputFileData%UnSum > 0 ) THEN
      
  
         IF ( InputFileData%PotMod == 1 .AND.  InputFileData%WAMIT%RdtnMod == 1) THEN
            ! Write the header for this section:  Note: When NBodyMod = 1 the kernel is now 6*NBody by 6*Nbody in size,
            !   and we have NBody 6 by 6 kernels for NBodyMod=2 or 3
            if (p%NBodyMod == 1) then
               ! NBodyMod=1 kernel printout which is 6*NBody x 6*NBody long
               WRITE( InputFileData%UnSum,  '(//)' ) 
               WRITE( InputFileData%UnSum,  '(A)' ) 'Radiation memory effect kernel'
               WRITE( InputFileData%UnSum,  '(//)' ) 
               
               WRITE( InputFileData%UnSum, '(1X,A10,2X,A10)',ADVANCE='no' )    '    n    ' , '     t    '
               do i = 1,6*p%NBody
                  do j = 1,6*p%NBody
                     WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) 'K'//trim(num2lstr(i))//trim(num2lstr(j))
                  end do
               end do
               write(InputFileData%UnSum,'()')  ! end of line character
               
             
               WRITE( InputFileData%UnSum, '(1X,A10,2X,A10)',ADVANCE='no' )    '   (-)   ' , '    (s)   ' 
               do i = 1,6*p%NBody
                  do j = 1,6*p%NBody
                     if ( mod(i-1,6)+1 < 4 ) then
                        if ( mod(j-1,6)+1 < 4  ) then  
                           WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) ' (kg/s^2) '
                        else
                           WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) ' (kgm/s^2) '
                        end if  
                     else
                        if  ( mod(j-1,6)+1 < 4 )  then
                           WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) ' (kgm/s^2) '
                        else
                           WRITE( InputFileData%UnSum, '(2X,A16)',ADVANCE='no' ) '(kgm^2/s^2)'
                        end if                      
                     end if
                  end do
               end do
               write(InputFileData%UnSum,'()')  ! end of line character
                 
               do k= 0,p%WAMIT(1)%Conv_Rdtn%NStepRdtn-1
                  WRITE( InputFileData%UnSum, '(1X,I10,2X,E12.5)',ADVANCE='no' ) K, K*p%WAMIT(1)%Conv_Rdtn%RdtnDT
                  do i = 1,6*p%NBody
                     do j = 1,6*p%NBody
                        WRITE( InputFileData%UnSum, '(2X,ES16.5)',ADVANCE='no' ) p%WAMIT(1)%Conv_Rdtn%RdtnKrnl(k,i,j)
                     end do
                  end do
                  write(InputFileData%UnSum,'()')  ! end of line character
               end do
               
            else
               do j = 1,p%nWAMITObj
                  WRITE( InputFileData%UnSum,  '(//)' ) 
                  WRITE( InputFileData%UnSum,  '(A)' ) 'Radiation memory effect kernel'
                  WRITE( InputFileData%UnSum,  '(//)' ) 
                  WRITE( InputFileData%UnSum, '(1X,A10,2X,A10,21(2X,A16))' )    '    n    ' , '     t    ', '   K11    ', '   K12    ', '    K13   ', '    K14    ', '    K15    ', '    K16    ', '    K22   ', '    K23   ', '    K24    ', '    K25    ', '    K26    ', '    K33    ', '    K34    ', '    K35    ',     'K36    ', '    K44    ', '    K45    ', '    K46    ', '    K55    ', '    K56    ', '    K66    '
                  WRITE( InputFileData%UnSum, '(1X,A10,2X,A10,21(2X,A16))' )    '   (-)   ' , '    (s)   ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kg/s^2) ', ' (kg/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kg/s^2)  ', ' (kgm/s^2) ', ' (kgm/s^2) ', ' (kgm/s^2) ', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)', '(kgm^2/s^2)'

               ! Write the data
                  DO I = 0,p%WAMIT(j)%Conv_Rdtn%NStepRdtn-1
                     WRITE( InputFileData%UnSum, '(1X,I10,2X,E12.5,21(2X,ES16.5))' ) I, I*p%WAMIT(j)%Conv_Rdtn%RdtnDT, &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,1), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,2), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,3), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,4), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,1,6), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,2), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,3), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,4), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,2,6), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,3,3), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,3,4), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,3,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,3,6), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,4,4), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,4,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,4,6), &
                        p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,5,5), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,5,6), p%WAMIT(j)%Conv_Rdtn%RdtnKrnl(I,6,6)
                  END DO
               end do
            end if
         END IF
         
      END IF

!==========================================
      
         ! Close the summary file
      IF ( InputFileData%HDSum ) THEN
         CALL HDOut_CloseSum( InputFileData%UnSum, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL CleanUp()
               RETURN
            END IF
      END IF
      
   ! Create the input mesh associated with kinematics of the platform reference point       
      CALL MeshCreate( BlankMesh        = u%PRPMesh         &
                     ,IOS               = COMPONENT_INPUT   &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = ErrStat2          &
                     ,ErrMess           = ErrMsg2           &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.            )    
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
          
      CALL MeshPositionNode (u%PRPMesh                             &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , ErrStat2                           &
                              , ErrMsg2                            )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)    

      CALL MeshConstructElement (  u%PRPMesh             &
                                    , ELEMENT_POINT      &                         
                                    , ErrStat2           &
                                    , ErrMsg2            &
                                    , 1                  )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   
      CALL MeshCommit ( u%PRPMesh             &
                        , ErrStat2            &
                        , ErrMsg2             )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
      
      u%PRPMesh%RemapFlag  = .TRUE.


      ! Create the input mesh associated with kinematics of the various WAMIT bodies
      IF (p%PotMod >= 1) THEN
         CALL MeshCreate( BlankMesh        = u%WAMITMesh       &
                        ,IOS               = COMPONENT_INPUT   &
                        ,Nnodes            = p%NBody           &
                        ,ErrStat           = ErrStat2          &
                        ,ErrMess           = ErrMsg2           &
                        ,TranslationDisp   = .TRUE.            &
                        ,Orientation       = .TRUE.            &
                        ,TranslationVel    = .TRUE.            &
                        ,RotationVel       = .TRUE.            &
                        ,TranslationAcc    = .TRUE.            &
                        ,RotationAcc       = .TRUE.            )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         do iBody = 1, p%NBody
            theta = (/ 0.0_R8Ki, 0.0_R8Ki, 0.0_R8Ki /)
            orientation = EulerConstruct(theta)

            CALL MeshPositionNode (u%WAMITMesh                                                                              &
                                 , iBody                                                                                    &
                                 , (/InputFileData%PtfmRefxt(iBody), InputFileData%PtfmRefyt(iBody), InputFileData%PtfmRefzt(iBody)/)   &
                                 , ErrStat2                                                                                 &
                                 , ErrMsg2                                                                                  &
                                 , orientation )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

            CALL MeshConstructElement (  u%WAMITMesh        &
                                       , ELEMENT_POINT      &
                                       , ErrStat2           &
                                       , ErrMsg2            &
                                       , iBody              )
               CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         end do

         CALL MeshCommit ( u%WAMITMesh           &
                           , ErrStat2            &
                           , ErrMsg2             )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      ! Output mesh for loads at each WAMIT body
         CALL MeshCopy (   SrcMesh    = u%WAMITMesh            &
                        ,DestMesh     = y%WAMITMesh            &
                        ,CtrlCode     = MESH_SIBLING           &
                        ,IOS          = COMPONENT_OUTPUT       &
                        ,ErrStat      = ErrStat2               &
                        ,ErrMess      = ErrMsg2                &
                        ,Force        = .TRUE.                 &
                        ,Moment       = .TRUE.                 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:y%WAMITMesh')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF
         u%WAMITMesh%RemapFlag  = .TRUE.
         y%WAMITMesh%RemapFlag  = .TRUE.
      ENDIF    ! PotMod > 1


   ! Create helper mesh to map all Hydrodynamics loads to the platform reference point to (0,0,0)
      CALL MeshCreate (  BlankMesh      = m%AllHdroOrigin   &
                     ,IOS               = COMPONENT_OUTPUT  &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = ErrStat2          &
                     ,ErrMess           = ErrMsg2           &
                     ,Force             = .TRUE.            &
                     ,Moment            = .TRUE.            )

         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'HydroDyn_Init:m%AllHdroOrigin')

      CALL MeshPositionNode (m%AllHdroOrigin                       &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , ErrStat2                           &
                              , ErrMsg2                            )
      
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL MeshConstructElement (  m%AllHdroOrigin       &
                                    , ELEMENT_POINT      &                         
                                    , ErrStat2           &
                                    , ErrMsg2            &
                                    , 1                  )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      CALL MeshCommit ( m%AllHdroOrigin     &
                        , ErrStat2            &
                        , ErrMsg2             )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF
      m%AllHdroOrigin%RemapFlag  = .TRUE.

         ! Set some more parameters
      p%OutSwtch      = InputFileData%OutSwtch 
      p%Delim         = ''
      p%OutFmt        = InputFileData%OutFmt
      p%OutSFmt       = InputFileData%OutSFmt
      p%NumOuts       = InputFileData%NumOuts


      CALL HDOUT_Init( HydroDyn_ProgDesc, InitInp%OutRootName, InputFileData, y,  p, m, InitOut, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF

      if ( y%WAMITMesh%Committed ) then 
         call MeshMapCreate( y%WAMITMesh,    m%AllHdroOrigin, m%HD_MeshMap%W_P_2_PRP_P, ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      end if
      
      IF ( y%Morison%Mesh%Committed ) THEN 
         CALL MeshMapCreate( y%Morison%Mesh, m%AllHdroOrigin, m%HD_MeshMap%M_P_2_PRP_P,  ErrStat2, ErrMsg2  );CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      ENDIF
      
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF         
         
      ! Define initialization-routine output here:
      InitOut%Ver = HydroDyn_ProgDesc         
         
      
      !............................................................................................
      ! Module Variables:
      !............................................................................................

      call HydroDyn_InitVars(u, p, x, y, m, InitOut, InputFileData, InitInp%Linearize, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      !............................................................................................
      ! Initialize Jacobian:
      !............................................................................................
      ! if (InitInp%Linearize) then      
      !    call HD_Init_Jacobian( p, u, y, InitOut, ErrStat2, ErrMsg2)
      !    call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      ! end if

      IF ( p%OutSwtch == 1 ) THEN ! Only HD-level output writing
         ! HACK  WE can tell FAST not to write any HD outputs by simply deallocating the WriteOutputHdr array!
         DEALLOCATE ( InitOut%WriteOutputHdr )
      END IF
      
         ! Destroy the local initialization data
      CALL CleanUp()
         
CONTAINS
!................................
   SUBROUTINE CleanUp()
      ! Use DEALLOCATEpointers = .false.
      ! NOTE: All of the pointer data originated in SeaState, and SeaState is responsible for deallocating the data
      !        all other modules are responsible for nullifying their versions of the pointers when they are done with the data

      CALL HydroDyn_DestroyInputFile( InputFileData,   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      CALL NWTC_Library_DestroyFileInfoType(InFileInfo,ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)  

   END SUBROUTINE CleanUp
!................................
END SUBROUTINE HydroDyn_Init


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE HydroDyn_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )

      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(HydroDyn_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(HydroDyn_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(HydroDyn_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other/optimization states            
      TYPE(HydroDyn_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      

      INTEGER(IntKi)                         :: ErrStat2                            ! local error status
      CHARACTER(ErrMsgLen)                   :: ErrMsg2                             ! local error message
      CHARACTER(*), PARAMETER                :: RoutineName = 'HydroDyn_End'

         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               

      
            
         ! Write the HydroDyn-level output file data FROM THE LAST COMPLETED TIME STEP if the user requested module-level output
         ! and the current time has advanced since the last stored time step.
         
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3) THEN               
         CALL HDOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      END IF          
      
         ! Close files here:  
      CALL HDOut_CloseOutput( p, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
          

         ! Destroy the input data: (ignore errors)
      CALL HydroDyn_DestroyInput( u, ErrStat2, ErrMsg2 )


         ! Destroy the parameter data: (ignore errors)
      CALL HydroDyn_DestroyParam( p, ErrStat2, ErrMsg2 )


         ! Destroy the state data: (ignore errors)
         
      CALL HydroDyn_DestroyContState(   x,           ErrStat2, ErrMsg2 )
      CALL HydroDyn_DestroyDiscState(   xd,          ErrStat2, ErrMsg2 )
      CALL HydroDyn_DestroyConstrState( z,           ErrStat2, ErrMsg2 )
      CALL HydroDyn_DestroyOtherState(  OtherState,  ErrStat2, ErrMsg2 )

         ! Destroy misc variables: (ignore errors)
      CALL HydroDyn_DestroyMisc( m, ErrStat2, ErrMsg2 )

         ! Destroy the output data: (ignore errors)
      CALL HydroDyn_DestroyOutput( y, ErrStat2, ErrMsg2 )
      

END SUBROUTINE HydroDyn_End

subroutine HydroDyn_InitVars(u, p, x, y, m, InitOut, InputFileData, Linearize, ErrStat, ErrMsg)
   type(HydroDyn_InputType),           intent(inout)  :: u              !< An initial guess for the input; input mesh must be defined
   type(HydroDyn_ParameterType),       intent(inout)  :: p              !< Parameters
   type(HydroDyn_ContinuousStateType), intent(inout)  :: x              !< Continuous state
   type(HydroDyn_OutputType),          intent(inout)  :: y              !< Initial system outputs (outputs are not calculated;
   type(HydroDyn_MiscVarType),         intent(inout)  :: m              !< Misc variables for optimization (not copied in glue code)
   type(HydroDyn_InitOutputType),      intent(inout)  :: InitOut        !< Output for initialization routine
   type(HydroDyn_InputFile),           intent(in)     :: InputFileData  !< Input file data
   logical,                      intent(in)     :: Linearize      !< Flag to initialize linearization variables
   integer(IntKi),               intent(out)    :: ErrStat        !< Error status of the operation
   character(*),                 intent(out)    :: ErrMsg         !< Error message if ErrStat /= ErrID_None

   character(*), parameter       :: RoutineName = 'HydroDyn_InitVars'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2

   integer(IntKi)                :: i, j, k
   real(R8Ki)                    :: PerturbTrans, PerturbRot, Perturbs(6)
   character(10)                 :: BodyDesc
   character(10), parameter      :: dofLabels(6) = &
      ['PtfmSg', 'PtfmSw', 'PtfmHv', 'PtfmR ', 'PtfmP ', 'PtfmY ']

   ! Allocate space for variables (deallocate if already allocated)
   if (associated(p%Vars)) deallocate(p%Vars)
   allocate(p%Vars, stat=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, "Error allocating vars", ErrStat, ErrMsg, RoutineName)
      return
   end if

   ErrStat = ErrID_None
   ErrMsg = ""

   ! Associate pointer in init output
   InitOut%Vars => p%Vars

   !----------------------------------------------------------------------------
   ! Continuous State Variables
   !----------------------------------------------------------------------------

   ! Need to determine how many wamit body objects there are
   p%totalExctnStates = 0
   p%totalRdtnStates = 0
   do j = 1, p%nWAMITObj
      p%totalExctnStates = p%totalExctnStates + p%WAMIT(j)%SS_Exctn%numStates  ! numStates defaults to zero in the case where ExctnMod = 0 instead of 2
      p%totalRdtnStates  = p%totalRdtnStates  + p%WAMIT(j)%SS_Rdtn%numStates   ! numStates defaults to zero in the case where RdtnMod  = 0 instead of 2
   end do
   p%totalStates = p%totalExctnStates + p%totalRdtnStates
   
   ! Initialize body description to empty
   BodyDesc = ""

   ! Get excitation
   do k = 1, p%nWAMITObj
      if (p%WAMIT(k)%SS_Exctn%numStates == 0) cycle
      if (p%NBody > 1) BodyDesc = 'B'//trim(Num2LStr(k))
      call MV_AddVar(p%Vars%x, "WAMIT("//trim(Num2LStr(k))//")%SS_Exctn", FieldScalar, &
                     Flags=VF_DerivOrder1, &
                     Num=p%WAMIT(k)%SS_Exctn%numStates, &
                     Perturb=20000.0_R8Ki * D2R_D, &
                     LinNames=[((trim(BodyDesc)//'Exctn'//trim(dofLabels(j))//Num2LStr(i), i = 1, p%WAMIT(k)%SS_Exctn%spDOF(j)), j = 1, 6)])
   end do

   do k = 1, p%nWAMITObj
      if (p%WAMIT(k)%SS_Rdtn%numStates == 0) cycle
      if (p%NBody > 1) BodyDesc = 'B'//trim(Num2LStr(k))
      call MV_AddVar(p%Vars%x, "WAMIT("//trim(Num2LStr(k))//")%SS_Rdtn", FieldScalar, &
                     Flags=VF_DerivOrder1, &
                     Num=p%WAMIT(k)%SS_Rdtn%numStates, &
                     Perturb=2.0_R8Ki * D2R_D , &
                     LinNames=[((trim(BodyDesc)//'Rdtn'//trim(dofLabels(j))//Num2LStr(i), i = 1, p%WAMIT(k)%SS_Rdtn%spDOF(j)), j = 1, 6)])
   end do

   !----------------------------------------------------------------------------
   ! Input variables
   !----------------------------------------------------------------------------

   ! Translation and rotation perturbations
   PerturbTrans = 0.02_R8Ki*D2R * max(real(p%WaveField%EffWtrDpth, R8Ki), 1.0_R8Ki) 
   PerturbRot = 2*D2R

   ! Create perturbation array (order based on MotionFields)
   Perturbs = [PerturbTrans, &  ! FieldTransDisp
               PerturbRot, &    ! FieldOrientation
               PerturbTrans, &  ! FieldTransVel
               PerturbRot, &    ! FieldAngularVel
               PerturbTrans, &  ! FieldTransAcc
               PerturbRot]      ! FieldAngularAcc

   call MV_AddMeshVar(p%Vars%u, "Morison", MotionFields, u%Morison%Mesh, &
                      VarIdx=p%iVarMorisonMotionMesh, &
                      Perturbs=Perturbs)

   call MV_AddMeshVar(p%Vars%u, "WAMIT", MotionFields, u%WAMITMesh, &
                      VarIdx=p%iVarWAMITMotionMesh, &
                      Perturbs=Perturbs)

   call MV_AddMeshVar(p%Vars%u, "Platform-RefPt", MotionFields, u%PRPMesh, &
                      VarIdx=p%iVarPRPMotionMesh, &
                      Perturbs=Perturbs)

   call MV_AddVar(p%Vars%u, "WaveElev0", FieldScalar, &
                  VarIdx=p%iVarWaveElev0, &
                  Flags=VF_ExtLin + VF_Linearize, &
                  LinNames=['Extended input: wave elevation at platform ref point, m'])

   call MV_AddVar(p%Vars%u, "HWindSpeed", FieldScalar, &
                  VarIdx=p%iVarHWindSpeed, &
                  Flags=VF_ExtLin + VF_Linearize, &
                  LinNames=['Extended input: horizontal current speed (steady/uniform wind), m/s'])

   call MV_AddVar(p%Vars%u, "PLexp", FieldScalar, &
                  VarIdx=p%iVarPLexp, &
                  Flags=VF_ExtLin + VF_Linearize, &
                  LinNames=['Extended input: vertical power-law shear exponent, -'])

   call MV_AddVar(p%Vars%u, "PropagationDir", FieldScalar, &
                  VarIdx=p%iVarPropagationDir, &
                  Flags=VF_ExtLin + VF_Linearize, &
                  LinNames=['Extended input: propagation direction, rad'])

   !----------------------------------------------------------------------------
   ! Output variables
   !----------------------------------------------------------------------------

   call MV_AddMeshVar(p%Vars%y, "MorisonLoads", LoadFields, y%Morison%Mesh, VarIdx=p%iVarMorisonLoadMesh)

   call MV_AddMeshVar(p%Vars%y, "WAMITLoads", LoadFields, y%WAMITMesh, VarIdx=p%iVarWAMITLoadMesh)

   call MV_AddVar(p%Vars%y, "WriteOutput", FieldScalar, &
                  VarIdx=p%iVarWriteOut, &
                  Flags=VF_WriteOut, &
                  Num=p%NumTotalOuts, &
                  LinNames=[(WriteOutputLinName(i), i = 1, p%NumTotalOuts)])

   !----------------------------------------------------------------------------
   ! Initialize Variables and Jacobian data
   !----------------------------------------------------------------------------

   call MV_InitVarsJac(p%Vars, m%Jac, Linearize, ErrStat2, ErrMsg2); if (Failed()) return

   call HydroDyn_CopyContState(x, m%x_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return
   call HydroDyn_CopyContState(x, m%dxdt_lin, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return
   call HydroDyn_CopyInput(u, m%u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return
   call HydroDyn_CopyOutput(y, m%y_lin, MESH_NEWCOPY, ErrStat2, ErrMsg2); if (Failed()) return

contains
   character(LinChanLen) function WriteOutputLinName(idx)
      integer(IntKi), intent(in) :: idx
      WriteOutputLinName = trim(InitOut%WriteOutputHdr(idx))//', '//trim(InitOut%WriteOutputUnt(idx))
   end function
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
   end function Failed
end subroutine


!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE HydroDyn_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

      REAL(DbKi),                         INTENT(IN   )  :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n               !< Current step of the simulation: t = n*Interval
      TYPE(HydroDyn_InputType),           INTENT(INOUT ) :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   )  :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p               !< Parameters
      TYPE(HydroDyn_ContinuousStateType), INTENT(INOUT)  :: x               !< Input: Continuous states at t;
                                                                            !!   Output: Continuous states at t + Interval
      TYPE(HydroDyn_DiscreteStateType),   INTENT(INOUT)  :: xd              !< Input: Discrete states at t;
                                                                            !!   Output: Discrete states at t + Interval
      TYPE(HydroDyn_ConstraintStateType), INTENT(INOUT)  :: z               !< Input: Constraint states at t;
                                                                            !!   Output: Constraint states at t + Interval
      TYPE(HydroDyn_OtherStateType),      INTENT(INOUT)  :: OtherState      !< Other states: Other states at t;
                                                                            !!   Output: Other states at t + Interval
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m               !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

         ! Local variables
      INTEGER                                            :: I, iWAMIT       ! Generic loop counters
      INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
      INTEGER                                            :: nTime           ! number of inputs 

      TYPE(WAMIT_InputType), ALLOCATABLE                 :: Inputs_WAMIT(:)  
      TYPE(Morison_InputType), ALLOCATABLE               :: Inputs_Morison(:)  
      TYPE(Morison_InputType)                            :: u_Morison
      CHARACTER(*), PARAMETER                            :: RoutineName = 'HydroDyn_UpdateStates'
      
          ! Create dummy variables required by framework but which are not used by the module
            
#ifdef USE_FIT      
      TYPE(FIT_InputType), ALLOCATABLE                 :: Inputs_FIT(:) 
      TYPE(FIT_ConstraintStateType)      :: FIT_z              ! constraint states
      TYPE(FIT_ContinuousStateType)      :: FIT_x              ! Input: Continuous states at t;
#endif      
      
         ! Initialize variables

      ErrStat   = ErrID_None           ! no error has occurred
      ErrMsg    = ""
      nTime = size(Inputs)

      IF (INPUTS(1)%Morison%Mesh%Committed) THEN
      
         ALLOCATE( Inputs_Morison(nTime), STAT = ErrStat2 )
         IF (ErrStat2 /=0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Failed to allocate array Inputs_Morison.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
         
         DO i=1,nTime
            CALL Morison_CopyInput(Inputs(i)%Morison, Inputs_Morison(i), MESH_NEWCOPY, ErrStat2, ErrMsg2)
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END DO
         CALL Morison_CopyInput(Inputs(1)%Morison, u_Morison, MESH_NEWCOPY, ErrStat2, ErrMsg2)
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
         CALL Morison_Input_ExtrapInterp(Inputs_Morison, InputTimes, u_Morison, t, ErrStat2, ErrMsg2) ! get inputs at time t
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
         IF (ErrStat < AbortErrLev) THEN
            ! Update the discrete states of Morison - The state of the high-pass velocity filter
            CALL Morison_UpdateDiscState( t, u_Morison, p%Morison, x%Morison, xd%Morison, &
                                   z%Morison, OtherState%Morison, m%Morison, ErrStat2, ErrMsg2 )
               call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         END IF
            
         call Morison_DestroyInput(u_Morison, ErrStat2, ErrMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
         do i=1,size(Inputs_Morison)
            call Morison_DestroyInput(Inputs_Morison(i), ErrStat2, ErrMsg2)
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         end do
         deallocate(Inputs_Morison)
            
      END IF
         
         ! Return without doing any work if the we are not using a potential flow model
      IF ( p%PotMod == 0  ) THEN
         RETURN
      ELSEIF ( p%PotMod == 1 ) THEN
       
         ALLOCATE( Inputs_WAMIT(nTime), STAT = ErrStat2 )
         IF (ErrStat2 /=0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Failed to allocate array Inputs_WAMIT.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF

         if ( p%NBodyMod == 1 .or. p%NBody == 1 ) then
            ! For this NBodyMod or NBody=1, there is only one WAMIT object, so copy the necessary inputs and then call WAMIT_UpdateStates
            do I=1,nTime
                  ! Copy the inputs from the HD mesh into the WAMIT mesh         
               call MeshCopy( Inputs(I)%WAMITMesh, Inputs_WAMIT(I)%Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
               call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
            end do
         
            if (ErrStat < AbortErrLev) then    ! if there was an error copying the input meshes, we'll skip this step and then cleanup the temporary input meshes     
                  ! Update the WAMIT module states
      
               call WAMIT_UpdateStates( t, n, Inputs_WAMIT, InputTimes, p%WAMIT(1), x%WAMIT(1), xd%WAMIT(1), z%WAMIT, OtherState%WAMIT(1), m%WAMIT(1), ErrStat2, ErrMsg2 )
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         

            end if
         
         else
         
            ! We have multiple WAMIT objects

               ! Loop over number of inputs and copy them into an array of WAMIT inputs
            do iWAMIT = 1, p%nWAMITObj
            
               do I=1,nTime
                     ! We need to create to valid mesh data structures in our Inputs_WAMIT(I)%Mesh using the miscvar version as a template, but the actually data will be generated below      
                  call MeshCopy( m%u_WAMIT(iWAMIT)%Mesh, Inputs_WAMIT(I)%Mesh, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               end do
               if (ErrStat > AbortErrLev) exit
               
               do I=1,nTime
                     ! We need to copy the iWAMIT-th node data from the Inputs(I)%WAMITMesh onto the 1st node of the Inputs_WAMIT(I)%Mesh
                  Inputs_WAMIT(I)%Mesh%TranslationDisp(:,1)  = Inputs(I)%WAMITMesh%TranslationDisp(:,iWAMIT)
                  Inputs_WAMIT(I)%Mesh%Orientation    (:,:,1)= Inputs(I)%WAMITMesh%Orientation    (:,:,iWAMIT)
                  Inputs_WAMIT(I)%Mesh%TranslationVel (:,1)  = Inputs(I)%WAMITMesh%TranslationVel (:,iWAMIT)
                  Inputs_WAMIT(I)%Mesh%RotationVel    (:,1)  = Inputs(I)%WAMITMesh%RotationVel    (:,iWAMIT)
                  Inputs_WAMIT(I)%Mesh%TranslationAcc (:,1)  = Inputs(I)%WAMITMesh%TranslationAcc (:,iWAMIT)
                  Inputs_WAMIT(I)%Mesh%RotationAcc    (:,1)  = Inputs(I)%WAMITMesh%RotationAcc    (:,iWAMIT)               
               end do
            
                  ! UpdateStates for the iWAMIT-th body
               call WAMIT_UpdateStates( t, n, Inputs_WAMIT, InputTimes, p%WAMIT(iWAMIT), x%WAMIT(iWAMIT), xd%WAMIT(iWAMIT), z%WAMIT, OtherState%WAMIT(iWAMIT), m%WAMIT(iWAMIT), ErrStat2, ErrMsg2 )
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )     
                  if (ErrStat > AbortErrLev) exit
               
            end do
         
         end if
       
            ! deallocate temporary inputs
         do I=1,nTime
            call WAMIT_DestroyInput( Inputs_WAMIT(I), ErrStat2, ErrMsg2 )
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         end do
      
         deallocate(Inputs_WAMIT)

#ifdef USE_FIT      
   ELSE IF ( p%PotMod == 2 ) THEN  ! FIT
      
      ALLOCATE( Inputs_FIT(nTime), STAT = ErrStat2 )
      IF (ErrStat2 /=0) THEN
         CALL SetErrStat( ErrID_Fatal, 'Failed to allocate array Inputs_FIT.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

         
         ! Loop over number of inputs and copy them into an array of FIT inputs
      
      DO I=1,nTime
         
            ! Copy the inputs from the HD mesh into the FIT input variables
         
            ! Determine the rotational angles from the direction-cosine matrix
         rotdisp = GetSmllRotAngs ( Inputs(I)%WAMITMesh%Orientation(:,:,1), ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_UpdateStates' )   
         Inputs_FIT(I)%roll     = rotdisp(1)
         Inputs_FIT(I)%pitch    = rotdisp(2)
         Inputs_FIT(I)%yaw      = rotdisp(3)
         Inputs_FIT(I)%si_t(:)  = Inputs(I)%WAMITMesh%TranslationDisp(:,1)             
         Inputs_FIT(I)%vel_t(:) = Inputs(I)%WAMITMesh%TranslationVel (:,1)  
      END DO
      
         
         
         ! Update the FIT module states
     
      CALL FIT_UpdateStates( t, n, Inputs_FIT, InputTimes, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
     

         ! deallocate temporary inputs
      DO I=1,nTime
         CALL FIT_DestroyInput( Inputs_FIT(I), ErrStat2, ErrMsg2 )     
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      END DO
      
      DEALLOCATE(Inputs_FIT) 
#endif

   END IF
   
   !CALL Cleanup()
   
contains
   subroutine Cleanup()
   
      if (allocated(Inputs_Morison)) then
         do i=1,size(Inputs_Morison)
            call Morison_DestroyInput(Inputs_Morison(i), ErrStat2, ErrMsg2)
         end do
         deallocate(Inputs_Morison)
      end if
      call Morison_DestroyInput(u_Morison, ErrStat2, ErrMsg2)
      
   
      if (allocated(Inputs_WAMIT)) then
         do i=1,size(Inputs_WAMIT)
            call Wamit_DestroyInput(Inputs_WAMIT(i), ErrStat2, ErrMsg2)
         end do
         deallocate(Inputs_WAMIT)
      end if
      
   end subroutine Cleanup
      
END SUBROUTINE HydroDyn_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE HydroDyn_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )   
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (note that this is intent out because we're copying the u%WAMITMesh into m%u_wamit%mesh)
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(HydroDyn_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                        !!   nectivity information does not have to be recalculated) + for previous WriteOutput results
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !! Error message if ErrStat /= ErrID_None

      INTEGER                                            :: I, J        ! Generic counters
      
      INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None

#ifdef USE_FIT       
      TYPE(FIT_ContinuousStateType)        :: FIT_x             ! Initial continuous states
      TYPE(FIT_ConstraintStateType)        :: FIT_z             ! Initial guess of the constraint states 
      TYPE(FIT_InputType)                  :: Inputs_FIT
#endif      
     
      REAL(ReKi)                           :: q(6*p%NBody), qdot(6*p%NBody), qdotsq(6*p%NBody), qdotdot(6*p%NBody)
      REAL(ReKi)                           :: rotdisp(3)                              ! small angle rotational displacements
      integer(IntKi)                       :: iBody, indxStart, indxEnd  ! Counters
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""
       
      
      ! Write the Hydrodyn-level output file data FROM THE LAST COMPLETED TIME STEP if the user requested module-level output
      ! and the current time has advanced since the last stored time step. Note that this must be done before filling y%WriteOutput
      ! so that we don't get recent results. Also note that this may give strange results in the .HD.out files of linearization simulations.
         
      IF ( (p%OutSwtch == 1 .OR. p%OutSwtch == 3) .AND. ( Time > m%LastOutTime ) ) THEN               
         CALL HDOut_WriteOutputs( m%LastOutTime, y, p, m%Decimate, ErrStat2, ErrMsg2 )         
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
      END IF
      m%LastOutTime   = Time ! time associated with the next values of y%WriteOutput
      
         
      
         !-------------------------------------------------------------------
         ! Additional stiffness, damping forces.  These need to be placed on a point mesh which is located at the WAMIT reference point (WRP).
         ! This mesh will need to get mapped by the glue code for use by either ElastoDyn or SubDyn.
         !-------------------------------------------------------------------



      if ( p%PotMod == 1 ) then
         do iBody = 1, p%NBody
               ! Determine the rotational angles from the direction-cosine matrix
            rotdisp = GetSmllRotAngs ( u%WAMITMesh%Orientation(:,:,iBody), ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
            indxStart = (iBody-1)*6+1
            indxEnd   = indxStart+5
            q      (indxStart:indxEnd)   = reshape((/real(u%WAMITMesh%TranslationDisp(:,iBody),ReKi),rotdisp(:)/),(/6/))
            qdot   (indxStart:indxEnd)   = reshape((/u%WAMITMesh%TranslationVel(:,iBody),u%WAMITMesh%RotationVel(:,iBody)/),(/6/))
            qdotsq (indxStart:indxEnd)   = abs(qdot(indxStart:indxEnd))*qdot(indxStart:indxEnd)
            qdotdot(indxStart:indxEnd)   = reshape((/u%WAMITMesh%TranslationAcc(:,iBody),u%WAMITMesh%RotationAcc(:,iBody)/),(/6/))
         end do
   !FIXME: Error handling appears to be broken here.
         if ( p%NBodyMod == 1 ) then
              
                  ! Compute the load contirbution from user-supplied added stiffness and damping        
               m%F_PtfmAdd = p%AddF0(:,1) - matmul(p%AddCLin(:,:,1), q) - matmul(p%AddBLin(:,:,1), qdot) - matmul(p%AddBQuad(:,:,1), qdotsq)
            do iBody = 1, p%NBody
               indxStart = (iBody-1)*6+1
               indxEnd   = indxStart+5
                  ! Attach to the output point mesh
               y%WAMITMesh%Force (:,iBody) = m%F_PtfmAdd(indxStart:indxStart+2)
               y%WAMITMesh%Moment(:,iBody) = m%F_PtfmAdd(indxStart+3:indxEnd)
            end do
         
         else
            do iBody = 1, p%NBody
               indxStart = (iBody-1)*6+1
               indxEnd   = indxStart+5
  
               m%F_PtfmAdd(indxStart:indxEnd) = p%AddF0(:,iBody) - matmul(p%AddCLin(:,:,iBody), q(indxStart:indxEnd)) - matmul(p%AddBLin(:,:,iBody), qdot(indxStart:indxEnd)) - matmul(p%AddBQuad(:,:,iBody), qdotsq(indxStart:indxEnd))
      
                  ! Attach to the output point mesh
               y%WAMITMesh%Force (:,iBody) = m%F_PtfmAdd(indxStart:indxStart+2)
               y%WAMITMesh%Moment(:,iBody) = m%F_PtfmAdd(indxStart+3:indxEnd)
            end do
         
         end if
       
         m%F_Waves = 0.0_ReKi

         if (allocated(m%u_WAMIT)) then               ! Check that we allocated u_WAMIT, otherwise there is an error checking the mesh
            if ( m%u_WAMIT(1)%Mesh%Committed ) then  ! Make sure we are using WAMIT / there is a valid mesh

               if ( p%NBodyMod == 1 .or. p%NBody == 1 ) then
                     ! Copy the inputs from the HD mesh into the WAMIT mesh
                  call MeshCopy( u%WAMITMesh, m%u_WAMIT(1)%Mesh, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
                     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
                        if ( ErrStat >= AbortErrLev ) return

                  call WAMIT_CalcOutput( Time, m%u_WAMIT(1), p%WAMIT(1), x%WAMIT(1), xd%WAMIT(1),  &
                                          z%WAMIT, OtherState%WAMIT(1), y%WAMIT(1), m%WAMIT(1), ErrStat2, ErrMsg2 )
                     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
                  do iBody=1,p%NBody
                     y%WAMITMesh%Force (:,iBody) = y%WAMITMesh%Force (:,iBody) + y%WAMIT(1)%Mesh%Force (:,iBody)
                     y%WAMITMesh%Moment(:,iBody) = y%WAMITMesh%Moment(:,iBody) + y%WAMIT(1)%Mesh%Moment(:,iBody)

                  end do
                     ! Copy the F_Waves1 information to the HydroDyn level so we can combine it with the 2nd order
                  m%F_Waves   = m%F_Waves + m%WAMIT(1)%F_Waves1
               else
                  do iBody=1,p%NBody

                        ! We need to copy the iWAMIT-th node data from the Inputs(I)%WAMITMesh onto the 1st node of the Inputs_WAMIT(I)%Mesh
                     m%u_WAMIT(iBody)%Mesh%TranslationDisp(:,1)  = u%WAMITMesh%TranslationDisp(:,iBody)
                     m%u_WAMIT(iBody)%Mesh%Orientation    (:,:,1)= u%WAMITMesh%Orientation    (:,:,iBody)
                     m%u_WAMIT(iBody)%Mesh%TranslationVel (:,1)  = u%WAMITMesh%TranslationVel (:,iBody)
                     m%u_WAMIT(iBody)%Mesh%RotationVel    (:,1)  = u%WAMITMesh%RotationVel    (:,iBody)
                     m%u_WAMIT(iBody)%Mesh%TranslationAcc (:,1)  = u%WAMITMesh%TranslationAcc (:,iBody)
                     m%u_WAMIT(iBody)%Mesh%RotationAcc    (:,1)  = u%WAMITMesh%RotationAcc    (:,iBody)

                     call WAMIT_CalcOutput( Time, m%u_WAMIT(iBody), p%WAMIT(iBody), x%WAMIT(iBody), xd%WAMIT(iBody),  &
                                          z%WAMIT, OtherState%WAMIT(iBody), y%WAMIT(iBody), m%WAMIT(iBody), ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
                     y%WAMITMesh%Force (:,iBody) = y%WAMITMesh%Force (:,iBody) + y%WAMIT(iBody)%Mesh%Force (:,1)
                     y%WAMITMesh%Moment(:,iBody) = y%WAMITMesh%Moment(:,iBody) + y%WAMIT(iBody)%Mesh%Moment(:,1)

                        ! Copy the F_Waves1 information to the HydroDyn level so we can combine it with the 2nd order
                     indxStart = (iBody-1)*6+1
                     indxEnd   = indxStart+5
                     m%F_Waves(indxStart:indxEnd) = m%F_Waves(indxStart:indxEnd) + m%WAMIT(iBody)%F_Waves1
                  end do
               end if
            end if   ! m%u_WAMIT(1)%Mesh%Committed
         end if      ! m%u_WAMIT is allocated



            ! Second order
         if (p%WAMIT2used) then

            if ( p%NBodyMod == 1 .or. p%NBody == 1 ) then
               call WAMIT2_CalcOutput( Time, p%WaveField, p%WAMIT2(1), y%WAMIT2(1), m%WAMIT2(1), ErrStat2, ErrMsg2 )
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
               do iBody=1,p%NBody
                  y%WAMITMesh%Force (:,iBody) = y%WAMITMesh%Force (:,iBody) + y%WAMIT2(1)%Mesh%Force (:,iBody)
                  y%WAMITMesh%Moment(:,iBody) = y%WAMITMesh%Moment(:,iBody) + y%WAMIT2(1)%Mesh%Moment(:,iBody)
               end do
                  ! Add F_Waves2 to m%F_Waves
               m%F_Waves  = m%F_Waves + m%WAMIT2(1)%F_Waves2
            else
               do iBody=1,p%NBody

                  call WAMIT2_CalcOutput( Time, p%WaveField, p%WAMIT2(iBody), y%WAMIT2(iBody), m%WAMIT2(iBody), ErrStat2, ErrMsg2 )
                     call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )
                  y%WAMITMesh%Force (:,iBody) = y%WAMITMesh%Force (:,iBody) + y%WAMIT2(iBody)%Mesh%Force (:,1)
                  y%WAMITMesh%Moment(:,iBody) = y%WAMITMesh%Moment(:,iBody) + y%WAMIT2(iBody)%Mesh%Moment(:,1)

                     ! Copy the F_Waves1 information to the HydroDyn level so we can combine it with the 2nd order
                  indxStart = (iBody-1)*6+1
                  indxEnd   = indxStart+5
                  m%F_Waves(indxStart:indxEnd) = m%F_Waves(indxStart:indxEnd) + m%WAMIT2(iBody)%F_Waves2
               end do
            end if
         end if         !  p%WAMIT2used

#ifdef USE_FIT          
      ELSE IF ( p%PotMod ==2 ) THEN !FIT
         Inputs_FIT%roll     = rotdisp(1)
         Inputs_FIT%pitch    = rotdisp(2)
         Inputs_FIT%yaw      = rotdisp(3)
         Inputs_FIT%si_t(:)  = u%WAMITMesh%TranslationDisp(:,1)             
         Inputs_FIT%vel_t(:) = u%WAMITMesh%TranslationVel (:,1)  
         CALL FIT_CalcOutput( Time, Inputs_FIT, p%FIT, FIT_x, xd%FIT, FIT_z, OtherState%FIT, y%FIT, ErrStat2, ErrMsg2 ) 
         
            ! Add FIT forces to the HydroDyn output mesh
         y%WAMITMesh%Force (:,1) = y%WAMITMesh%Force (:,1) + y%FIT%F(:)
         y%WAMITMesh%Moment(:,1) = y%WAMITMesh%Moment(:,1) + y%FIT%M(:)
#endif  
         
      END IF
      


      IF ( u%Morison%Mesh%Committed ) THEN  ! Make sure we are using Morison / there is a valid mesh
         CALL Morison_CalcOutput( Time, u%Morison, p%Morison, x%Morison, xd%Morison,  &
                                z%Morison, OtherState%Morison, y%Morison, m%Morison, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      END IF
      
         ! Integrate all the mesh loads onto the platfrom reference Point (PRP) at (0,0,0)
      m%F_Hydro = CalcLoadsAtWRP( y, u, m%AllHdroOrigin, m%HD_MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
        
      
      
         ! Map calculated results into the first p%NumOuts values of the y%WriteOutput Array
      CALL HDOut_MapOutputs( p, y, m%WAMIT, m%WAMIT2, m%F_PtfmAdd, m%F_Waves, m%F_Hydro, u%PRPMesh, q, qdot, qdotdot, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'HydroDyn_CalcOutput' )                  
      
      
         ! Aggregate the sub-module outputs 
      IF (p%Morison%NumOuts > 0) THEN
         J = p%NumOuts + 1
         DO I=1, p%Morison%NumOuts
            y%WriteOutput(J) = y%Morison%WriteOutput(I)
            J = J + 1
         END DO
      END IF
      
      m%LastOutTime   = Time
      
END SUBROUTINE HydroDyn_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE HydroDyn_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )  
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)
      TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters                             
      TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states                    
      TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
      TYPE(HydroDyn_ContinuousStateType), INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      integer(IntKi)             :: iWAMIT        ! loop counter
      CHARACTER(*), PARAMETER    :: RoutineName = 'HydroDyn_CalcContStateDeriv'
               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute the first time derivatives of the continuous states here:
   if ( .not. m%u_WAMIT(1)%Mesh%Committed ) return  ! Make sure we are using WAMIT / there is a valid mesh 

   call HydroDyn_CopyContState( x, dxdt, MESH_NEWCOPY, ErrStat, ErrMsg)
      if ( ErrStat >= AbortErrLev ) return
   
   if ( p%NBodyMod == 1 .or. p%NBody == 1 ) then
         ! Copy the inputs from the HD mesh into the WAMIT mesh
      call MeshCopy( u%WAMITMesh, m%u_WAMIT(1)%Mesh, MESH_UPDATECOPY, ErrStat, ErrMsg )   
         if ( ErrStat >= AbortErrLev ) return
      
      call WAMIT_CalcContStateDeriv( Time, m%u_WAMIT(1), p%WAMIT(1), x%WAMIT(1), xd%WAMIT(1), z%WAMIT, OtherState%WAMIT(1), m%WAMIT(1), dxdt%WAMIT(1), ErrStat, ErrMsg ) 
   else
           
      ! We have multiple WAMIT objects
         
         ! Loop over number of inputs and copy them into an array of WAMIT inputs
      do iWAMIT = 1, p%nWAMITObj
    
            ! We need to copy the iWAMIT-th node data from the Inputs(I)%WAMITMesh onto the 1st node of the Inputs_WAMIT(I)%Mesh
         m%u_WAMIT(iWAMIT)%Mesh%TranslationDisp(:,1)  = u%WAMITMesh%TranslationDisp(:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%Orientation    (:,:,1)= u%WAMITMesh%Orientation    (:,:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%TranslationVel (:,1)  = u%WAMITMesh%TranslationVel (:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%RotationVel    (:,1)  = u%WAMITMesh%RotationVel    (:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%TranslationAcc (:,1)  = u%WAMITMesh%TranslationAcc (:,iWAMIT)
         m%u_WAMIT(iWAMIT)%Mesh%RotationAcc    (:,1)  = u%WAMITMesh%RotationAcc    (:,iWAMIT)               

            ! UpdateStates for the iWAMIT-th body
         call WAMIT_CalcContStateDeriv( Time, m%u_WAMIT(iWAMIT), p%WAMIT(iWAMIT), x%WAMIT(iWAMIT), xd%WAMIT(iWAMIT), z%WAMIT, OtherState%WAMIT(iWAMIT), m%WAMIT(iWAMIT), dxdt%WAMIT(iWAMIT), ErrStat, ErrMsg ) 
            if (ErrStat > AbortErrLev) exit
               
      end do
         
   end if
   
END SUBROUTINE HydroDyn_CalcContStateDeriv




!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE HydroDyn_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )   
   
   REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds   
   TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (intent OUT only because we're copying the input mesh)              
   TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters                           
   TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
   TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other/optimization states                    
   TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
   TYPE(HydroDyn_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using  
                                                                     !!     the input values described above      
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""               
      
   ! Nothing to do here since none of the sub-modules have contraint states
   z_residual = z  
    
         ! Solve for the constraint states here:


END SUBROUTINE HydroDyn_CalcConstrStateResidual


!----------------------------------------------------------------------------------------------------------------------------------
function CalcLoadsAtWRP( y, u, AllHdroOrigin, MeshMapData, ErrStat, ErrMsg )

   type(HydroDyn_OutputType),  intent(inout)  :: y                        ! Hydrodyn outputs
   type(HydroDyn_InputType),   intent(in   )  :: u                        ! Hydrodyn inputs
   type(MeshType),             intent(inout)  :: AllHdroOrigin            ! This is the mesh which data is mapped onto.  We pass it in to avoid allocating it at each call
   type(HD_ModuleMapType),     intent(inout)  :: MeshMapData              ! Mesh mapping data structures 
   integer(IntKi),             intent(  out)  :: ErrStat                  ! Error status of the operation
   character(*),               intent(  out)  :: ErrMsg                   ! Error message if ErrStat /= ErrID_None                                                         
   real(ReKi)                                 :: CalcLoadsAtWRP(6)

      ! local variables
   integer(IntKi)                                 :: ErrStat2             ! temporary Error status of the operation
   character(ErrMsgLen)                           :: ErrMsg2              ! temporary Error message if ErrStat /= ErrID_None
   
   CalcLoadsAtWRP = 0.0_ReKi
   
   if ( y%WAMITMesh%Committed  ) then

      ! Just transfer the loads because the meshes are at the same location (0,0,0)
      call Transfer_Point_to_Point( y%WAMITMesh, AllHdroOrigin, MeshMapData%W_P_2_PRP_P, ErrStat2, ErrMsg2, u%WAMITMesh, u%PRPMesh )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CalcLoadsAtWRP')
            if (ErrStat >= AbortErrLev) return
            
      CalcLoadsAtWRP(1:3)  = CalcLoadsAtWRP(1:3)  + AllHdroOrigin%Force(:,1)
      CalcLoadsAtWRP(4:6)  = CalcLoadsAtWRP(4:6)  + AllHdroOrigin%Moment(:,1)

   end if      
      
   
   if ( y%Morison%Mesh%Committed ) then 

      call Transfer_Point_to_Point( y%Morison%Mesh, AllHdroOrigin, MeshMapData%M_P_2_PRP_P, ErrStat2, ErrMsg2,  u%Morison%Mesh, u%PRPMesh )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CalcLoadsAtWRP')
            if (ErrStat >= AbortErrLev) return
 
      CalcLoadsAtWRP(1:3)  = CalcLoadsAtWRP(1:3) + AllHdroOrigin%Force(:,1)
      CalcLoadsAtWRP(4:6)  = CalcLoadsAtWRP(4:6) + AllHdroOrigin%Moment(:,1)
         
   end if
   
end function CalcLoadsAtWRP
 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in HD_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and dZ/du are returned.
SUBROUTINE HD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu, FlagFilter )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions (Y) with respect 
                                                                               !!   to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) with 
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) with 
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z) with 
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   integer(IntKi), OPTIONAL,             INTENT(IN   )           :: FlagFilter !< Flag filter for variable calculation

   CHARACTER(*), PARAMETER       :: RoutineName = 'HD_JacobianPInput'
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   logical                       :: IsFullLin
   integer(IntKi)                :: FlagFilterLoc
   INTEGER(IntKi)                :: i, j, k, col
   INTEGER(IntKi)                :: startingI, startingJ, bOffset, offsetI
   
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Set full linearization flag and local filter flag
   if (present(FlagFilter)) then
      IsFullLin = FlagFilter == VF_None
      FlagFilterLoc = FlagFilter
   else
      IsFullLin = .true.
      FlagFilterLoc = VF_None
   end if
   
   ! make a copy of the inputs to perturb
   call HydroDyn_CopyInput(u, m%u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2)
   if (Failed()) return
   
   ! Pack inputs into array
   call HD_PackInputValues(p, u, m%Jac%u)
   if (Failed()) return

   ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
   IF ( PRESENT( dYdu ) ) THEN

      ! allocate dYdu if necessary
      if (.not. allocated(dYdu)) then
         call AllocAry(dYdu, p%Vars%Ny, p%Vars%Nu, 'dYdu', ErrStat2, ErrMsg2)
         if (Failed()) return
      end if
      
      ! Loop through input variables
      do i = 1, size(p%Vars%u)

         ! If variable flag not in flag filter, skip
         if (.not. MV_HasFlags(p%Vars%u(i), FlagFilterLoc)) cycle

         ! If variable is extended input, skip
         if (MV_HasFlags(p%Vars%u(i), VF_ExtLin)) cycle

         ! Loop through number of linearization perturbations in variable
         do j = 1, p%Vars%u(i)%Num

            ! Calculate positive perturbation
            call MV_Perturb(p%Vars%u(i), j, 1, m%Jac%u, m%Jac%u_perturb)
            call HD_UnpackInputValues(p, m%Jac%u_perturb, m%u_perturb)
            call HydroDyn_CalcOutput(t, m%u_perturb, p, x, xd, z, OtherState, m%y_lin, m, ErrStat2, ErrMsg2); if (Failed()) return
            call HD_PackOutputValues(p, m%y_lin, m%Jac%y_pos, IsFullLin)

            ! Calculate negative perturbation
            call MV_Perturb(p%Vars%u(i), j, -1, m%Jac%u, m%Jac%u_perturb)
            call HD_UnpackInputValues(p, m%Jac%u_perturb, m%u_perturb)
            call HydroDyn_CalcOutput(t, m%u_perturb, p, x, xd, z, OtherState, m%y_lin, m, ErrStat2, ErrMsg2); if (Failed()) return
            call HD_PackOutputValues(p, m%y_lin, m%Jac%y_neg, IsFullLin)

            ! Calculate column index
            col = p%Vars%u(i)%iLoc(1) + j - 1

            ! Get partial derivative via central difference and store in full linearization array
            call MV_ComputeCentralDiff(p%Vars%y, p%Vars%u(i)%Perturb, m%Jac%y_pos, m%Jac%y_neg, dYdu(:,col))
         end do
      end do

      ! Set extended inputs
      dYdu(:, p%Vars%u(p%iVarWaveElev0)%iLoc(1)) = 0.0_R8Ki
      dYdu(:, p%Vars%u(p%iVarHWindSpeed)%iLoc(1)) = 0.0_R8Ki
      dYdu(:, p%Vars%u(p%iVarPLexp)%iLoc(1)) = 0.0_R8Ki
      dYdu(:, p%Vars%u(p%iVarPropagationDir)%iLoc(1)) = 0.0_R8Ki
      
   END IF
   
   ! Calculate the partial derivative of the continuous state functions (X) with respect to the inputs (u) here:
   IF ( PRESENT( dXdu ) ) THEN

      ! For the case where either RdtnMod=0 and ExtcnMod=0 and hence %SS_Rdtn data or %SS_Exctn data is not valid then we do not have states, so simply return
      ! The key here is to never allocate the dXdu and related state Jacobian arrays because then the glue-code will behave properly
      
      ! allocate dXdu if necessary
      if (.not. allocated(dXdu)) then
         call AllocAry(dXdu, p%Vars%Nx, p%Vars%Nu, 'dXdu', ErrStat2, ErrMsg2)
         if (Failed()) return
      end if
      
      offsetI = 0
      dXdu = 0.0_R8Ki
      
      do j = 1,p%nWAMITObj
         do i = 1,p%WAMIT(j)%SS_Exctn%numStates
            dXdu(offsetI+i,p%Vars%Nu) = p%WAMIT(j)%SS_Exctn%B(i) ! B is numStates by 1
         end do
         offsetI = offsetI + p%WAMIT(j)%SS_Exctn%numStates
      end do

      startingI = p%totalStates - p%totalRdtnStates
      startingJ = p%Vars%Nu - 4 - 18 - 4*3*p%NBody !  subtract 4 for extended inputs and 4*3*NBody to place us at the beginning of the velocity inputs
      ! B is numStates by 6*NBody where NBody =1 if NBodyMod=2 or 3, but could be >1 for NBodyMod=1
      if ( p%NBodyMod == 1 ) then
         ! Example for NBodyMod=1 and NBody = 2,
         ! dXdu(:,startingIndx + 1) = p%WAMIT(1)%SS_Rdtn%B(:,1)
         ! dXdu(:,startingIndx + 2) = p%WAMIT(1)%SS_Rdtn%B(:,2)
         ! dXdu(:,startingIndx + 3) = p%WAMIT(1)%SS_Rdtn%B(:,3)
         ! dXdu(:,startingIndx + 4) = p%WAMIT(1)%SS_Rdtn%B(:,7)
         ! dXdu(:,startingIndx + 5) = p%WAMIT(1)%SS_Rdtn%B(:,8)
         ! dXdu(:,startingIndx + 6) = p%WAMIT(1)%SS_Rdtn%B(:,9)
         ! dXdu(:,startingIndx + 7) = p%WAMIT(1)%SS_Rdtn%B(:,4)   ! start of rotationalVel
         ! dXdu(:,startingIndx + 8) = p%WAMIT(1)%SS_Rdtn%B(:,5)
         ! dXdu(:,startingIndx + 9) = p%WAMIT(1)%SS_Rdtn%B(:,6)
         ! dXdu(:,startingIndx +10) = p%WAMIT(1)%SS_Rdtn%B(:,10)
         ! dXdu(:,startingIndx +11) = p%WAMIT(1)%SS_Rdtn%B(:,11)
         ! dXdu(:,startingIndx +12) = p%WAMIT(1)%SS_Rdtn%B(:,12)

         do i = 1,p%WAMIT(1)%SS_Rdtn%numStates
            k=0
            ! First set all translationalVel components
            do j = 1, p%WAMIT(1)%NBody
               bOffset = (j-1)*6
               dXdu(startingI+i,startingJ+k+1) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+1) 
               dXdu(startingI+i,startingJ+k+2) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+2)
               dXdu(startingI+i,startingJ+k+3) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+3)   
               k = k + 3
            end do
            ! Now set all rotationalVel components
            do j = 1, p%WAMIT(1)%NBody
               bOffset = (j-1)*6
               dXdu(startingI+i,startingJ+k+1) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+4)
               dXdu(startingI+i,startingJ+k+2) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+5)
               dXdu(startingI+i,startingJ+k+3) = p%WAMIT(1)%SS_Rdtn%B(i,bOffset+6)   
               k = k + 3
            end do
         end do         
      else
         ! Example NBodyMod=2or3 and NBody = 2,
         ! dXdu(:,startingIndx + 1) = p%WAMIT(1)%SS_Rdtn%B(:,1)
         ! dXdu(:,startingIndx + 2) = p%WAMIT(1)%SS_Rdtn%B(:,2)
         ! dXdu(:,startingIndx + 3) = p%WAMIT(1)%SS_Rdtn%B(:,3)
         ! dXdu(:,startingIndx + 4) = p%WAMIT(2)%SS_Rdtn%B(:,1)
         ! dXdu(:,startingIndx + 5) = p%WAMIT(2)%SS_Rdtn%B(:,2)
         ! dXdu(:,startingIndx + 6) = p%WAMIT(2)%SS_Rdtn%B(:,3)
         ! dXdu(:,startingIndx + 7) = p%WAMIT(1)%SS_Rdtn%B(:,4)   ! start of rotationalVel
         ! dXdu(:,startingIndx + 8) = p%WAMIT(1)%SS_Rdtn%B(:,5)
         ! dXdu(:,startingIndx + 9) = p%WAMIT(1)%SS_Rdtn%B(:,6)
         ! dXdu(:,startingIndx +10) = p%WAMIT(2)%SS_Rdtn%B(:,4)
         ! dXdu(:,startingIndx +11) = p%WAMIT(2)%SS_Rdtn%B(:,5)
         ! dXdu(:,startingIndx +12) = p%WAMIT(2)%SS_Rdtn%B(:,6)

         k=0
         offsetI=0
         ! First set all translationalVel components
         do j = 1, p%nWAMITObj
            do i = 1,p%WAMIT(j)%SS_Rdtn%numStates
               dXdu(startingI+offsetI+i,startingJ+k+1) = p%WAMIT(j)%SS_Rdtn%B(i,1) 
               dXdu(startingI+offsetI+i,startingJ+k+2) = p%WAMIT(j)%SS_Rdtn%B(i,2)
               dXdu(startingI+offsetI+i,startingJ+k+3) = p%WAMIT(j)%SS_Rdtn%B(i,3) 
            end do
            k = k + 3
            offsetI = offsetI + p%WAMIT(j)%SS_Rdtn%numStates
         end do
         ! Now set all rotationalVel components
         offsetI=0
         do j = 1, p%nWAMITObj
            do i = 1,p%WAMIT(j)%SS_Rdtn%numStates
               dXdu(startingI+offsetI+i,startingJ+k+1) = p%WAMIT(1)%SS_Rdtn%B(i,4)
               dXdu(startingI+offsetI+i,startingJ+k+2) = p%WAMIT(1)%SS_Rdtn%B(i,5)
               dXdu(startingI+offsetI+i,startingJ+k+3) = p%WAMIT(1)%SS_Rdtn%B(i,6)   
            end do
            k = k + 3
            offsetI = offsetI + p%WAMIT(j)%SS_Rdtn%numStates
         end do
      end if   
      
   END IF

   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu)
   END IF

   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu)
   END IF
   
contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
   end function Failed  
END SUBROUTINE HD_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
SUBROUTINE HD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx, FlagFilter )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdx(:,:)  !< Partial derivatives of output functions (Y) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdx(:,:)  !< Partial derivatives of continuous state functions (X) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddx(:,:) !< Partial derivatives of discrete state functions (Xd) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdx(:,:)  !< Partial derivatives of constraint state functions (Z) with respect 
                                                                               !!   to the continuous states (x) [intent in to avoid deallocation]
   integer(IntKi), OPTIONAL,             INTENT(IN   )           :: FlagFilter !< Flag filter for variable calculation
   
   CHARACTER(*), PARAMETER                           :: RoutineName = 'HD_JacobianPContState'
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   logical                                           :: IsFullLin
   integer(IntKi)                                    :: FlagFilterLoc
   INTEGER(IntKi)                                    :: i, j, k, col, sOffset

   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Set full linearization flag and local filter flag
   if (present(FlagFilter)) then
      IsFullLin = FlagFilter == VF_None
      FlagFilterLoc = FlagFilter
   else
      IsFullLin = .true.
      FlagFilterLoc = VF_None
   end if
   
   ! Copy State values to perturb
   call HydroDyn_CopyContState(x, m%x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2); if (Failed()) return
   call HD_PackStateValues(p, x, m%Jac%x)
   
   ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:
   if (present(dYdx)) then

      ! allocate dYdx if necessary
      if (.not. allocated(dYdx)) then
         call AllocAry(dYdx, p%Vars%Ny, p%Vars%Nx, 'dYdx', ErrStat2, ErrMsg2); if (Failed()) return
      end if
      
      ! Loop through state variables
      do i = 1, size(p%Vars%x)

         ! If variable flag not in flag filter, skip
         if (.not. MV_HasFlags(p%Vars%x(i), FlagFilterLoc)) cycle

         ! Loop through number of linearization perturbations in variable
         do j = 1, p%Vars%x(i)%Num

            ! Calculate positive perturbation
            call MV_Perturb(p%Vars%x(i), j, 1, m%Jac%x, m%Jac%x_perturb)
            call HD_UnpackStateValues(p, m%Jac%x_perturb, m%x_perturb)
            call HydroDyn_CalcOutput(t, u, p, m%x_perturb, xd, z, OtherState, m%y_lin, m, ErrStat2, ErrMsg2); if (Failed()) return
            call HD_PackOutputValues(p, m%y_lin, m%Jac%y_pos, IsFullLin)

            ! Calculate negative perturbation
            call MV_Perturb(p%Vars%x(i), j, -1, m%Jac%x, m%Jac%x_perturb)
            call HD_UnpackStateValues(p, m%Jac%x_perturb, m%x_perturb)
            call HydroDyn_CalcOutput(t, u, p, m%x_perturb, xd, z, OtherState, m%y_lin, m, ErrStat2, ErrMsg2); if (Failed()) return
            call HD_PackOutputValues(p, m%y_lin, m%Jac%y_neg, IsFullLin)

            ! Calculate column index
            col = p%Vars%x(i)%iLoc(1) + j - 1

            ! Get partial derivative via central difference and store in full linearization array
            call MV_ComputeCentralDiff(p%Vars%y, p%Vars%x(i)%Perturb, m%Jac%y_pos, m%Jac%y_neg, dYdx(:,col))
         end do
      end do
            
   end if

   ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:
   IF (present(dXdx)) then

      ! allocate dXdu if necessary
      if (.not. allocated(dXdx)) then
         call AllocAry(dXdx, p%Vars%Nx, p%Vars%Nx, 'dXdx', ErrStat2, ErrMsg2); if (Failed()) return
      end if

      dXdx = 0.0_R8Ki
      
      ! Analytical Jacobians from State-space models
      if ( p%totalExctnStates > 0 ) then
         sOffset = 0
         do k=1,p%nWAMITObj
            do j=1,p%WAMIT(k)%SS_Exctn%numStates   
               do i=1,p%WAMIT(k)%SS_Exctn%numStates ! Loop through all active (enabled) DOFs
                  dXdx(i+sOffset, j+sOffset) = p%WAMIT(k)%SS_Exctn%A(i,j)
               end do
            end do
            sOffset = sOffset + p%WAMIT(k)%SS_Exctn%numStates
         end do
      end if
      if ( p%totalRdtnStates > 0 ) then
         sOffset = 0
         do k=1,p%nWAMITObj
            do j=1,p%WAMIT(k)%SS_Rdtn%numStates   
               do i=1,p%WAMIT(k)%SS_Rdtn%numStates ! Loop through all active (enabled) DOFs
                  dXdx(i+p%totalExctnStates+sOffset, j+p%totalExctnStates+sOffset) = p%WAMIT(k)%SS_Rdtn%A(i,j)
               end do
            end do
            sOffset = sOffset + p%WAMIT(k)%SS_Rdtn%numStates
         end do
      end if
      
   END IF

   IF ( PRESENT( dXddx ) ) THEN
      if (allocated(dXddx)) deallocate(dXddx)
   END IF

   IF ( PRESENT( dZdx ) ) THEN
      if (allocated(dZdx)) deallocate(dZdx)
   END IF
   
contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
   end function Failed  
END SUBROUTINE HD_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and dZ/dxd are returned.
SUBROUTINE HD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdxd(:,:) !< Partial derivatives of output functions
                                                                               !!  (Y) with respect to the discrete
                                                                               !!  states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdxd(:,:) !< Partial derivatives of continuous state
                                                                               !!   functions (X) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddxd(:,:)!< Partial derivatives of discrete state
                                                                               !!   functions (Xd) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdxd(:,:) !< Partial derivatives of constraint state
                                                                               !!   functions (Z) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Calculate the partial derivative of the output functions (Y) with respect to the discrete states (xd) here:
   IF ( PRESENT( dYdxd ) ) THEN
      ! allocate and set dYdxd
   END IF

   ! Calculate the partial derivative of the continuous state functions (X) with respect to the discrete states (xd) here:
   IF ( PRESENT( dXdxd ) ) THEN
      ! allocate and set dXdxd
   END IF

   ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the discrete states (xd) here:
   IF ( PRESENT( dXddxd ) ) THEN
      ! allocate and set dXddxd
   END IF

   ! Calculate the partial derivative of the constraint state functions (Z) with respect to the discrete states (xd) here:
   IF ( PRESENT( dZdxd ) ) THEN
      ! allocate and set dZdxd
   END IF

END SUBROUTINE HD_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and dZ/dz are returned.
SUBROUTINE HD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is   
                                                                               !!   available here so that mesh parameter information (i.e.,  
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdz(:,:)  !< Partial derivatives of output functions (Y) with respect 
                                                                               !!  to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdz(:,:)  !< Partial derivatives of continuous state functions (X) with respect 
                                                                               !!  to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddz(:,:) !< Partial derivatives of discrete state functions (Xd) with respect 
                                                                               !!  to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdz(:,:)  !< Partial derivatives of constraint state functions (Z) with respect 
                                                                               !! to the constraint states (z) [intent in to avoid deallocation]
   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Calculate the partial derivative of the output functions (Y) with respect to the constraint states (z) here:
   IF ( PRESENT( dYdz ) ) THEN
      ! allocate and set dYdz
   END IF

   ! Calculate the partial derivative of the continuous state functions (X) with respect to the constraint states (z) here:
   IF ( PRESENT( dXdz ) ) THEN
      ! allocate and set dXdz
   END IF

   ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the constraint states (z) here:
   IF ( PRESENT( dXddz ) ) THEN
      ! allocate and set dXddz
   END IF

   ! Calculate the partial derivative of the constraint state functions (Z) with respect to the constraint states (z) here:
   IF ( PRESENT( dZdz ) ) THEN
      ! allocate and set dZdz
   END IF


END SUBROUTINE HD_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE HD_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(HydroDyn_OutputType),                  INTENT(IN   )           :: y          !< Output at operating point
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: u_op(:)    !< values of linearized inputs
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: y_op(:)    !< values of linearized outputs
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: x_op(:)    !< values of linearized continuous states
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: xd_op(:)   !< values of linearized discrete states
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: z_op(:)    !< values of linearized constraint states

   CHARACTER(*), PARAMETER                           :: RoutineName = 'HD_GetOP'
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   INTEGER(IntKi)                                    :: i, j, index

   ErrStat = ErrID_None
   ErrMsg  = ''

   !..................................
   IF ( PRESENT( u_op ) ) THEN
      if (.not. allocated(u_op)) then
         call AllocAry(u_op, p%Vars%Nu, 'u_op', ErrStat2, ErrMsg2)
         if (Failed()) return
      end if
      call HD_PackInputValues(p, u, u_op)
   END IF

   !..................................
   if ( PRESENT( y_op ) ) then
      if (.not. allocated(y_op)) then
         call AllocAry(y_op, p%Vars%Ny, 'y_op', ErrStat2, ErrMsg2)
         if (Failed()) return
      end if
      call HD_PackOutputValues(p, y, y_op, .true.)
   end if   

   !..................................
   IF ( PRESENT( x_op ) ) THEN
      if (p%Vars%Nx == 0) return
      if ( y%WAMITMesh%Committed ) then
         if (.not. allocated(x_op)) then 
            call AllocAry(x_op, p%Vars%Nx, 'x_op', ErrStat2, ErrMsg2)
            if (Failed()) return
         end if
         call HD_PackStateValues(p, x, x_op)
      end if
   END IF

   !..................................
   IF ( PRESENT( dx_op ) ) THEN
      if (p%Vars%Nx == 0) return
      if ( y%WAMITMesh%Committed ) then
         if (.not. allocated(dx_op)) then 
            call AllocAry(dx_op, p%Vars%Nx, 'dx_op', ErrStat2, ErrMsg2)
            if (Failed()) return
         end if
         call HydroDyn_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, m%dxdt_lin, ErrStat2, ErrMsg2 ) 
         if (Failed()) return
         call HD_PackStateValues(p, m%dxdt_lin, dx_op)
      end if    
   END IF

   !..................................
   IF ( PRESENT( xd_op ) ) THEN
   END IF
   
   !..................................
   IF ( PRESENT( z_op ) ) THEN
   END IF

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
   end function Failed
END SUBROUTINE HD_GetOP

subroutine HD_PackStateValues(p, x, ary)
   type(HydroDyn_ParameterType), intent(in)        :: p
   type(HydroDyn_ContinuousStateType), intent(in)  :: x
   real(R8Ki), intent(out)                         :: ary(:)
   integer(IntKi)                                  :: i, j, k
   k = 1
   do j = 1, p%nWAMITObj
      do i = 1,p%WAMIT(j)%SS_Exctn%numStates ! Loop through all DOFs
         ary(k) = x%WAMIT(j)%SS_Exctn%x(i)
         k = k + 1
      end do
   end do
   do j = 1, p%nWAMITObj
      do i = 1,p%WAMIT(j)%SS_Rdtn%numStates ! Loop through all DOFs
         ary(k) = x%WAMIT(j)%SS_Rdtn%x(i)
         k = k + 1
      end do
   end do
end subroutine

subroutine HD_UnpackStateValues(p, ary, x)
   type(HydroDyn_ParameterType), intent(in)           :: p
   real(R8Ki), intent(in)                             :: ary(:)
   type(HydroDyn_ContinuousStateType), intent(inout)  :: x
   integer(IntKi)                                     :: i, j, k
   k = 1
   do j = 1, p%nWAMITObj
      do i = 1,p%WAMIT(j)%SS_Exctn%numStates ! Loop through all DOFs
         x%WAMIT(j)%SS_Exctn%x(i) = ary(k)
         k = k + 1
      end do
   end do
   do j = 1, p%nWAMITObj
      do i = 1,p%WAMIT(j)%SS_Rdtn%numStates ! Loop through all DOFs
         x%WAMIT(j)%SS_Rdtn%x(i) = ary(k)
         k = k + 1
      end do
   end do
end subroutine

subroutine HD_PackInputValues(p, u, Ary)
   type(HydroDyn_ParameterType), intent(in)  :: p
   type(HydroDyn_InputType), intent(in)      :: u
   real(R8Ki), intent(out)                   :: Ary(:)
   integer(IntKi)                            :: i
   call MV_Pack(p%Vars%u, p%iVarMorisonMotionMesh, u%Morison%Mesh, Ary)
   call MV_Pack(p%Vars%u, p%iVarWAMITMotionMesh, u%WAMITMesh, Ary)
   call MV_Pack(p%Vars%u, p%iVarPRPMotionMesh, u%PRPMesh, Ary)
   call MV_Pack(p%Vars%u, p%iVarWaveElev0, 0.0_R8Ki, Ary)      ! Extended input
   call MV_Pack(p%Vars%u, p%iVarHWindSpeed, 0.0_R8Ki, Ary)     ! Extended input
   call MV_Pack(p%Vars%u, p%iVarPLexp, 0.0_R8Ki, Ary)          ! Extended input
   call MV_Pack(p%Vars%u, p%iVarPropagationDir, 0.0_R8Ki, Ary) ! Extended input

!FIXME: when sea current from IfW/FlowField is enabled, this code must be updated and enabled
!      !------------------------------
!      ! Extended inputs -- Linearization is only possible with Steady or Uniform Wind, so take advantage of that here
!      !     Module/Mesh/Field:  HWindSpeed      = 37
!      !     Module/Mesh/Field:  PLexp           = 38
!      !     Module/Mesh/Field:  PropagationDir  = 39
!      call IfW_UniformWind_GetOP(p_AD%FlowField%Uniform, t, .false. , OP_out)
!      ! HWindSpeed
!      u_op(index) = OP_out(1);   index = index + 1
!      ! PLexp
!      u_op(index) = OP_out(2);   index = index + 1
!      ! PropagationDir (include AngleH in calculation if any)
!      u_op(index) = OP_out(3) + p_AD%FlowField%PropagationDir;   index = index + 1
end subroutine

subroutine HD_UnpackInputValues(p, Ary, u)
   type(HydroDyn_ParameterType), intent(in)  :: p
   real(R8Ki), intent(in)                    :: Ary(:)
   type(HydroDyn_InputType), intent(inout)   :: u
   integer(IntKi)                            :: i
   call MV_Unpack(p%Vars%u, p%iVarMorisonMotionMesh, Ary, u%Morison%Mesh)
   call MV_Unpack(p%Vars%u, p%iVarWAMITMotionMesh, Ary, u%WAMITMesh)
   call MV_Unpack(p%Vars%u, p%iVarPRPMotionMesh, Ary, u%PRPMesh)
   ! call MV_Unpack(p%Vars%u, p%iVarWaveElev0, Ary, )   ! Extended input
   ! u_op(index) = 0.0_R8Ki; index=index+1   ! WaveElev0 -- linearization not allowed for non-zero
   ! u_op(index) = 0.0_R8Ki; index=index+1   ! HWindSpeed
   ! u_op(index) = 0.0_R8Ki; index=index+1   ! PLexp
   ! u_op(index) = 0.0_R8Ki; index=index+1   ! PropagationDir

!FIXME: when sea current from IfW/FlowField is enabled, this code must be updated and enabled
!      !------------------------------
!      ! Extended inputs -- Linearization is only possible with Steady or Uniform Wind, so take advantage of that here
!      !     Module/Mesh/Field:  HWindSpeed      = 37
!      !     Module/Mesh/Field:  PLexp           = 38
!      !     Module/Mesh/Field:  PropagationDir  = 39
!      call IfW_UniformWind_GetOP(p_AD%FlowField%Uniform, t, .false. , OP_out)
!      ! HWindSpeed
!      u_op(index) = OP_out(1);   index = index + 1
!      ! PLexp
!      u_op(index) = OP_out(2);   index = index + 1
!      ! PropagationDir (include AngleH in calculation if any)
!      u_op(index) = OP_out(3) + p_AD%FlowField%PropagationDir;   index = index + 1
end subroutine

subroutine HD_PackOutputValues(p, y, Ary, PackWriteOutput)
   type(HydroDyn_ParameterType), intent(in)  :: p
   type(HydroDyn_OutputType), intent(in)     :: y
   real(R8Ki), intent(out)                   :: Ary(:)
   logical, intent(in)                       :: PackWriteOutput
   integer(IntKi)                            :: i
   call MV_Pack(p%Vars%y, p%iVarMorisonLoadMesh, y%Morison%Mesh, Ary)
   call MV_Pack(p%Vars%y, p%iVarWAMITLoadMesh, y%WAMITMesh, Ary)
   call MV_Pack(p%Vars%y, p%iVarWriteOut, y%WriteOutput, Ary)
end subroutine

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE HydroDyn
!**********************************************************************************************************************************
