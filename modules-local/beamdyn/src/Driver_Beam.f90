!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD       
!   
!    This file is part of the NWTC Subroutine Library.
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
PROGRAM BeamDyn_Driver_Program

   USE BeamDyn_driver_subs  ! all other modules inherited through this one

   IMPLICIT NONE

   ! global glue-code-specific variables

   INTEGER(IntKi)                   :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                  :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   REAL(DbKi)                       :: dt_global        ! fixed/constant global time step
   REAL(DbKi)                       :: t_initial        ! time at initialization
   REAL(DbKi)                       :: t_final          ! time at simulation end 
   REAL(DbKi)                       :: t_global         ! global-loop time marker
   INTEGER(IntKi)                   :: n_t_final        ! total number of time steps
   INTEGER(IntKi)                   :: n_t_global       ! global-loop time counter
   INTEGER(IntKi)                   :: BD_interp_order  ! order of interpolation/extrapolation

   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BD_InitInputType)           :: BD_InitInput
   TYPE(BD_ParameterType)           :: BD_Parameter
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousState
   TYPE(BD_InitOutputType)          :: BD_InitOutput
   TYPE(BD_DiscreteStateType)       :: BD_DiscreteState
   TYPE(BD_ConstraintStateType)     :: BD_ConstraintState
   TYPE(BD_OtherStateType)          :: BD_OtherState
   TYPE(BD_MiscVarType)             :: BD_MiscVar
   TYPE(BD_InputType) ,ALLOCATABLE  :: BD_Input(:)
   REAL(DbKi),         ALLOCATABLE  :: BD_InputTimes(:)
   TYPE(BD_OutputType),ALLOCATABLE  :: BD_Output(:)
   REAL(DbKi),ALLOCATABLE           :: BD_OutputTimes(:)
   INTEGER(IntKi)                   :: DvrOut 
   REAL(BDKi)                       :: temp_POS(3)
   REAL(BDKi)                       :: temp_CRV(3)
   REAL(BDKi)                       :: temp_CRV2(3)
   REAL(R8Ki)                       :: DCM(3,3)          ! must be same type as mesh orientation fields
   REAL(ReKi)                       :: Pos(3)            ! must be same type as mesh position fields
   REAL(BDKi)                       :: TmpDCM(3,3)
   
   TYPE(BD_DriverInternalType)      :: DvrData

   CHARACTER(256)                   :: DvrInputFile
   CHARACTER(256)                   :: RootName


   ! local variables
   Integer(IntKi)                          :: i               ! counter for various loops
   REAL(R8Ki)                              :: start, finish
   REAL(BDKi) ,        ALLOCATABLE  :: IniVelo(:,:)        ! Initial Position Vector between origins of Global and blade frames [-]

   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'BeamDynDriverProgram'
   
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'BeamDyn Driver', 'v2.00.00', '9-May-2017' )  ! The version number of this program.
   

   ! -------------------------------------------------------------------------
   ! Initialization of library (especially for screen output)
   ! -------------------------------------------------------------------------      
   CALL NWTC_Init()
      ! Display the copyright notice
   CALL DispCopyrightLicense( version )   
      ! Tell our users what they're running
   CALL WrScr( ' Running '//GetNVD( version )//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )
      ! Give Git Hash info
#ifdef GIT_VERSION_INFO
   CALL WrScr( " Git version info: "//GIT_VERSION_INFO//NewLine )
#endif
   
   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------   
   
   CALL GET_COMMAND_ARGUMENT(1,DvrInputFile)
   CALL GetRoot(DvrInputFile,RootName)
   CALL BD_ReadDvrFile(DvrInputFile,t_initial,t_final,dt_global,BD_InitInput,DvrData,ErrStat,ErrMsg)
      CALL CheckError()
   BD_InitInput%RootName         = TRIM(BD_Initinput%InputFile)
   BD_InitInput%RootDisp(:)      = 0.0_R8Ki
   BD_InitInput%RootOri(:,:)     = 0.0_R8Ki
   BD_InitInput%RootOri(1:3,1:3) = BD_InitInput%GlbRot(1:3,1:3)
   t_global = t_initial
   n_t_final = ((t_final - t_initial) / dt_global )

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   BD_interp_order = 1

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BD_Input(BD_interp_order + 1)) 
   ALLOCATE(BD_InputTimes(BD_interp_order + 1)) 
   ALLOCATE(BD_Output(BD_interp_order + 1)) 
   ALLOCATE(BD_OutputTimes(BD_interp_order + 1)) 

   CALL BD_Init(BD_InitInput             &
                   , BD_Input(1)         &
                   , BD_Parameter        &
                   , BD_ContinuousState  &
                   , BD_DiscreteState    &
                   , BD_ConstraintState  &
                   , BD_OtherState       &
                   , BD_Output(1)        &
                   , BD_MiscVar          &
                   , dt_global           &
                   , BD_InitOutput       &
                   , ErrStat             &
                   , ErrMsg )
      CALL CheckError()
   
   call CreateMultiPointMeshes()   
   call Transfer_MultipointLoads()
   
!bjj: this is the driver's hack to get initial velocities for the input-output solve      
   CALL AllocAry(IniVelo,BD_Parameter%dof_node,BD_Parameter%node_total,'IniVelo',ErrStat,ErrMsg); 
      CALL CheckError()
   IniVelo = BD_ContinuousState%dqdt

   
   CALL Dvr_InitializeOutputFile(DvrOut,BD_InitOutput,RootName,ErrStat,ErrMsg)
      CALL CheckError()

   BD_InputTimes(1)  = t_initial
   BD_InputTimes(2)  = t_initial 
   BD_OutputTimes(1) = t_initial
   BD_OutputTimes(2) = t_initial

   !CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), BD_Parameter, BD_InitInput,IniVelo,ErrStat, ErrMsg)
   CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), BD_Parameter, DvrData,IniVelo,ErrStat, ErrMsg)
   CALL BD_CopyInput (BD_Input(1) , BD_Input(2) , MESH_NEWCOPY, ErrStat, ErrMsg)
   CALL BD_CopyOutput(BD_Output(1), BD_Output(2), MESH_NEWCOPY, ErrStat, ErrMsg)
      CALL CheckError()

      ! calcoutput
   CALL CPU_TIME(start)

   DO n_t_global = 0, n_t_final
     BD_InputTimes(2)  = BD_InputTimes(1) 
     BD_InputTimes(1)  = t_global + dt_global
     BD_OutputTimes(2) = BD_OutputTimes(1) 
     BD_OutputTimes(1) = t_global + dt_global
     CALL BD_InputSolve( BD_InputTimes(1), BD_Input(1), BD_Parameter, DvrData, IniVelo, ErrStat, ErrMsg)
     CALL BD_InputSolve( BD_InputTimes(2), BD_Input(2), BD_Parameter, DvrData, IniVelo, ErrStat, ErrMsg)
        CALL CheckError()
     CALL BD_CalcOutput( t_global, BD_Input(2), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                             BD_ConstraintState, &
                             BD_OtherState,  BD_Output(2), BD_MiscVar, ErrStat, ErrMsg)
        CALL CheckError()

     CALL Dvr_WriteOutputLine(t_global,DvrOut,BD_Parameter%OutFmt,BD_Output(2),ErrStat,ErrMsg)
        CALL CheckError()

     IF(BD_Parameter%analysis_type .EQ. BD_STATIC_ANALYSIS .AND. n_t_global > 8) EXIT 

     CALL BD_UpdateStates( t_global, n_t_global, BD_Input, BD_InputTimes, BD_Parameter, &
                               BD_ContinuousState, &
                               BD_DiscreteState, BD_ConstraintState, &
                               BD_OtherState, BD_MiscVar, ErrStat, ErrMsg )
        CALL CheckError()

      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial

   ENDDO
      
   CALL CPU_TIME(finish)
   
   WRITE(*,*) 'Start: ' , start
   WRITE(*,*) 'Finish: ', finish
   WRITE(*,*) 'Time: '  , finish-start

   CALL Dvr_End()

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE Dvr_End()

      character(ErrMsgLen)                          :: errMsg2                 ! temporary Error message if ErrStat /=
      integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
      character(*), parameter                       :: RoutineName = 'Dvr_End'

      IF(DvrOut >0) CLOSE(DvrOut)

      DO i=1,BD_interp_order + 1
          CALL BD_End( BD_Input(i), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                           BD_ConstraintState, BD_OtherState, BD_Output(i), BD_MiscVar, ErrStat2, ErrMsg2 )
      ENDDO 
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      IF(ALLOCATED(BD_InputTimes )) DEALLOCATE(BD_InputTimes )
      IF(ALLOCATED(BD_OutputTimes)) DEALLOCATE(BD_OutputTimes)
      if(allocated(DvrData%MultiPointLoad)) deallocate(DvrData%MultiPointLoad)

      if (ErrStat >= AbortErrLev) then      
         CALL ProgAbort( 'BeamDyn Driver encountered simulation error level: '&
             //TRIM(GetErrStr(ErrStat)), TrapErrors=.FALSE., TimeWait=3._ReKi )  ! wait 3 seconds (in case they double-clicked and got an error)
      else
         call NormStop()
      end if
   END SUBROUTINE Dvr_End
!----------------------------------------------------------------------------------------------------------------------------------
   subroutine CheckError()
   
      if (ErrStat /= ErrID_None) then
         call WrScr(TRIM(ErrMsg))
         
         if (ErrStat >= AbortErrLev) then
            call Dvr_End()
         end if
      end if
         
   end subroutine CheckError
!----------------------------------------------------------------------------------------------------------------------------------
   subroutine CreateMultiPointMeshes()
   
   ! DvrData%NumPointLoads is at least 1
   
   !.......................
   ! Mesh for multi-point loading on blades
   !.......................
   CALL MeshCreate( BlankMesh        = DvrData%mplMotion  &
                   ,IOS              = COMPONENT_INPUT    &
                   ,NNodes           = DvrData%NumPointLoads      &
                   ,TranslationDisp  = .TRUE.             &
                   ,Orientation      = .TRUE.             &
                   ,TranslationVel   = .TRUE.             &
                   ,RotationVel      = .TRUE.             &
                   ,TranslationAcc   = .TRUE.             &
                   ,RotationAcc      = .TRUE.             &
                   ,ErrStat          = ErrStat2           &
                   ,ErrMess          = ErrMsg2             )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! these nodes are placed along the key point line (as are the GLL nodes)
   DO i = 1,DvrData%NumPointLoads

       call Find_IniNode(BD_InitOutput%kp_coordinate, BD_Parameter, 1, BD_InitOutput%kp_total, DvrData%MultiPointLoad(i,1), temp_POS, temp_CRV)
       
       Pos = BD_Parameter%GlbPos + MATMUL(BD_Parameter%GlbRot,temp_POS)
       
       temp_CRV2 = MATMUL(BD_Parameter%GlbRot,temp_CRV)
       CALL BD_CrvCompose(temp_CRV,BD_Parameter%Glb_crv,temp_CRV2,FLAG_R1R2) !temp_CRV = p%Glb_crv composed with temp_CRV2

       CALL BD_CrvMatrixR(temp_CRV,TmpDCM) ! returns TmpDCM (the transpose of the DCM orientation matrix)

       ! possible type conversions here:
       DCM = TRANSPOSE(TmpDCM)

       ! set the reference position and orientation for each node.
       CALL MeshPositionNode ( Mesh    = DvrData%mplMotion     &
                              ,INode   = i             &
                              ,Pos     = Pos           &
                              ,ErrStat = ErrStat2      &
                              ,ErrMess = ErrMsg2       &
                              ,Orient  = DCM           )
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ENDDO
   
   DO i = 1,DvrData%NumPointLoads
       CALL MeshConstructElement( Mesh     = DvrData%mplMotion      &
                                 ,Xelement = ELEMENT_POINT    &
                                 ,P1       = i                &
                                 ,ErrStat  = ErrStat2         &
                                 ,ErrMess  = ErrMsg2          )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ENDDO
   CALL MeshCommit ( Mesh    = DvrData%mplMotion       &
                    ,ErrStat = ErrStat2        &
                    ,ErrMess = ErrMsg2         )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL CheckError()

   !.......................
   ! sibling mesh for motions, which are needed to transfer loads
   !.......................

   CALL MeshCopy ( SrcMesh  = DvrData%mplMotion    &
                 , DestMesh = DvrData%mplLoads     &
                 , CtrlCode = MESH_SIBLING         &
                 , IOS      = COMPONENT_INPUT      &
                 , Force    = .TRUE.               &
                 , Moment   = .TRUE.               &
                 , ErrStat  = ErrStat2             &
                 , ErrMess  = ErrMsg2              )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL CheckError()

      
   !.......................
   ! "sibling" mesh for BD motions, which are needed to transfer loads
   ! BD's input and output meshes may not be at the same location, so we need another copy
   !.......................
      
   CALL MeshCopy ( SrcMesh  = BD_Input(1)%PointLoad        &
                 , DestMesh = DvrData%y_BldMotion_at_u_point       &
                 , CtrlCode = MESH_COUSIN                  &  ! Like a sibling, except using new memory for position/refOrientation and elements
                 , IOS      = COMPONENT_OUTPUT             &
                 , Orientation     = .TRUE.                &
                 , TranslationDisp = .TRUE.                &
                 , ErrStat         = ErrStat2              &
                 , ErrMess         = ErrMsg2               )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL CheckError()
      
      
   !.......................
   ! initialize the mapping between the BD output motions and driver mpl motion mesh:
   !.......................
      
   CALL MeshMapCreate( BD_Output(1)%BldMotion, DvrData%mplMotion, DvrData%Map_BldMotion_to_mplMotion, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !.......................
   ! initialize the mapping between the driver mpl loads and BD input point loads mesh:
   !.......................
   CALL MeshMapCreate( DvrData%mplLoads, BD_Input(1)%PointLoad, DvrData%Map_mplLoads_to_PointLoad, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL CheckError()
      
   !.......................
   ! initialize the mapping between the BD output motions and BD motions at the nodes on the y%BldMotion mesh:
   !.......................
      
   CALL MeshMapCreate( BD_Output(1)%BldMotion, DvrData%y_BldMotion_at_u_point, DvrData%Map_y_BldMotion_to_u_point, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL CheckError()
      
      
   DvrData%mplMotion%remapFlag = .false.
   DvrData%mplLoads%remapFlag = .false.
   DvrData%y_BldMotion_at_u_point%remapFlag = .false.
      
   end subroutine CreateMultiPointMeshes
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Transfer_MultipointLoads()

      ! the mpl loads are constant over the course of the simulation, but they could be set up to vary
   DO i = 1,DvrData%NumPointLoads
      DvrData%mplLoads%Force(1:3,i)  = DvrData%MultiPointLoad(i,2:4)
      DvrData%mplLoads%Moment(1:3,i) = DvrData%MultiPointLoad(i,5:7)
   ENDDO
      
      ! get the motions from BD to the two meshes needed in this loads transfer:
   CALL Transfer_Line2_to_Point( BD_Output(1)%BldMotion, DvrData%y_BldMotion_at_u_point, DvrData%Map_y_BldMotion_to_u_point, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   CALL Transfer_Line2_to_Point( BD_Output(1)%BldMotion, DvrData%mplMotion, DvrData%Map_BldMotion_to_mplMotion, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      
      ! transfer the mpl loads to BD's input point mesh at the GLL nodes:
   CALL Transfer_Point_to_Point( DvrData%mplLoads, BD_Input(1)%PointLoad, DvrData%Map_mplLoads_to_PointLoad, ErrStat2, ErrMsg2, DvrData%mplMotion, DvrData%y_BldMotion_at_u_point)  
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL CheckError()

end subroutine Transfer_MultipointLoads
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_InputSolve( t, u,  p, DvrData, IniVelo, ErrStat, ErrMsg)
 
   REAL(DbKi),                  INTENT(IN   ) :: t
   TYPE(BD_InputType),          INTENT(INOUT) :: u
   TYPE(BD_ParameterType),      INTENT(IN   ) :: p
   TYPE(BD_DriverInternalType), INTENT(INOUT) :: DvrData
   REAL(BDKi),                  INTENT(IN   ) :: IniVelo(:,:)
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat          ! Error status of the operation
   CHARACTER(*),                INTENT(  OUT) :: ErrMsg           ! Error message if ErrStat /= ErrID_None
                                         
   ! local variables                     
   INTEGER(IntKi)                        :: i                ! do-loop counter
   REAL(BDKi)                            :: temp_vec(3)
   REAL(BDKi)                            :: temp_rr(3)
   REAL(BDKi)                            :: temp_r0(3)
   REAL(BDKi)                            :: temp_theta(3)
   REAL(BDKi)                            :: DCM(3,3)


   ! ----------------------------------------------------------------------
   
   ErrStat = ErrID_None
   ErrMsg  = ''

   temp_r0 = p%GlbPos

   temp_theta = MATMUL(p%GlbRot,IniVelo(4:6,1))
   temp_theta = temp_theta*t

   temp_vec(:) = 4.0_BDKi*TAN(temp_theta(:)/4.0_BDKi)
   CALL BD_CrvMatrixR(temp_vec,DCM)
   
   ! gather point forces and line forces

   ! Point mesh: RootMotion 
   ! Calculate root displacements and rotations
   u%RootMotion%Orientation(:,:,1) = MATMUL(p%GlbRot,DCM)
   !u%RootMotion%Orientation(:,:,1) = DCM
   
   !temp_rr(:) = MATMUL(u%RootMotion%Orientation(:,:,1),temp_r0)
   temp_rr(:) = MATMUL(DCM,temp_r0)
   !temp_rr(:) = u%RootMotion%TranslationDisp(:,:,1)
   u%RootMotion%Orientation(:,:,1)   = TRANSPOSE(u%RootMotion%Orientation(:,:,1))
   u%RootMotion%TranslationDisp(:,:) = 0.0_BDKi
   u%RootMotion%TranslationDisp(:,1) = temp_rr(:) - temp_r0(:)
   ! END Calculate root displacements and rotations
   
   ! Calculate root translational and angular velocities
   u%RootMotion%RotationVel(:,1) = MATMUL(p%GlbRot,IniVelo(4:6,1))
   
   u%RootMotion%TranslationVel(:,1) = MATMUL(SkewSymMat(real(u%RootMotion%RotationVel(:,1),BDKi)),temp_rr)
   ! END Calculate root translational and angular velocities


   ! Calculate root translational and angular accelerations
   u%RootMotion%RotationAcc   (:,:) = 0.0_BDKi
   u%RootMotion%TranslationAcc(:,1) = MATMUL(SkewSymMat(real(u%RootMotion%RotationVel(:,1),BDKi)), &
                                             u%RootMotion%TranslationVel(:,1))
   ! END Calculate root translational and angular accelerations

   
   
   ! set up the point load input:
   ! @VA: if we want to apply these at different positions, we should call Transfer_MultipointLoads(); alternatively, we could store the result of u%PointLoad the first
   ! time we call this and just use that instead of doing another transfer here.
   CALL Transfer_Point_to_Point( DvrData%mplLoads, u%PointLoad, DvrData%Map_mplLoads_to_PointLoad, ErrStat2, ErrMsg2, DvrData%mplMotion, DvrData%y_BldMotion_at_u_point)  
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL CheckError()
   
   u%PointLoad%Force(1:3,u%PointLoad%NNodes)  = u%PointLoad%Force(1:3,u%PointLoad%NNodes)  + DvrData%TipLoad(1:3)
   u%PointLoad%Moment(1:3,u%PointLoad%NNodes) = u%PointLoad%Moment(1:3,u%PointLoad%NNodes) + DvrData%TipLoad(4:6)
   
   ! LINE2 mesh: DistrLoad
   DO i=1,u%DistrLoad%NNodes
      u%DistrLoad%Force(:,i) =  DvrData%DistrLoad(1:3)
      u%DistrLoad%Moment(:,i)=  DvrData%DistrLoad(4:6)
   ENDDO

END SUBROUTINE BD_InputSolve

END PROGRAM BeamDyn_Driver_Program
