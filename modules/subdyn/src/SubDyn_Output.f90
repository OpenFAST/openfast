!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.
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
MODULE SubDyn_Output
   USE NWTC_Library
   USE SubDyn_Types
   USE SD_FEM
   USE SubDyn_Output_Params, only: MNfmKe, MNfmMe, MNTDss, MNRDe, MNTRAe, IntfSS, IntfTRss, IntfTRAss, ReactSS
   USE SubDyn_Output_Params, only: ParamIndxAry, ParamUnitsAry, ValidParamAry, SSqm01, SSqmd01, SSqmdd01

   IMPLICIT NONE

   ! The maximum number of output channels which can be output by the code.
   INTEGER(IntKi),PUBLIC, PARAMETER      :: MaxOutPts = 2265

   PRIVATE
      ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: SDOut_CloseSum
   PUBLIC :: SDOut_OpenSum
   PUBLIC :: SDOut_MapOutputs
   PUBLIC :: SDOut_OpenOutput
   PUBLIC :: SDOut_CloseOutput
   PUBLIC :: SDOut_WriteOutputNames
   PUBLIC :: SDOut_WriteOutputUnits
   PUBLIC :: SDOut_WriteOutputs
   PUBLIC :: SDOut_Init
   PUBLIC :: SD_Init_Jacobian
   PUBLIC :: SD_Perturb_u
   PUBLIC :: SD_Perturb_x
   PUBLIC :: SD_Compute_dY
   PUBLIC :: SD_Compute_dX

CONTAINS


!> This subroutine initializes the output module, checking if the output parameter list (OutList)
! contains valid names, and opening the output file if there are any requested outputs
SUBROUTINE SDOut_Init( Init, y,  p, misc, InitOut, WtrDpth, ErrStat, ErrMsg )
   TYPE(SD_InitType),               INTENT( INOUT ) :: Init                 ! data needed to initialize the output module
   TYPE(SD_OutputType),             INTENT( INOUT ) :: y                    ! SubDyn module's output data
   TYPE(SD_ParameterType), target,  INTENT( INOUT ) :: p                    ! SubDyn module paramters
   TYPE(SD_MiscVarType),            INTENT( INOUT ) :: misc                 ! SubDyn misc/optimization variables
   TYPE(SD_InitOutputType ),        INTENT( INOUT ) :: InitOut              ! SubDyn module initialization output data
   REAL(ReKi),                      INTENT( IN    ) :: WtrDpth              ! water depth from initialization routine
   INTEGER,                         INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred
   CHARACTER(*),                    INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER(IntKi)                 :: ErrStat2      ! Error status of the operation
   CHARACTER(ErrMsgLen)           :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                 :: I,J,K2 !Counters
   INTEGER(IntKi)                 :: iMember  ! Member index (not member ID)
   INTEGER(IntKi)                 :: iElem  ! Index of element in Element List
   INTEGER(IntKi)                 :: iNode  ! Index of node in Node list
   INTEGER(IntKi)                 :: iiElem ! Loop counter on element index
   INTEGER(IntKi)                 :: nElemPerNode, nNodesPerElem ! Number of elements connecting to a node, Number of nodes per elem 
   type(MeshAuxDataType), pointer :: pLst                                                   !< Alias to shorten notation and highlight code similarities
   real(ReKi), allocatable :: T_TIreact(:,:) ! Transpose of TIreact, temporary
   ErrStat = 0      
   ErrMsg=""

   p%OutAllDims=12*p%NMembers*2    !size of AllOut Member Joint forces

   ! Check that the variables in OutList are valid      
   CALL SDOut_ChkOutLst( Init%SSOutList, p,  ErrStat2, ErrMsg2 ); if(Failed()) return

   ! --- Allocation (size 0 if not outputs)
   !IF ( ALLOCATED( p%OutParam ) .AND. p%NumOuts > 0 ) THEN           ! Output has been requested           
   ! Allocate SDWrOuput which is used to store a time step's worth of output channels, prior to writing to a file.
   CALL AllocAry(misc%SDWrOutput       , p%NumOuts + p%OutAllInt*p%OutAllDims, 'SDWrOutupt' , ErrStat2, ErrMsg2) ; if(Failed()) return
   ! Allocate WriteOuput  
   CALL AllocAry(y%WriteOutput         , p%NumOuts + p%OutAllInt*p%OutAllDims, 'WriteOutput', ErrStat2, ErrMsg2); if(Failed()) return
   allocate(misc%AllOuts(0:MaxOutPts + p%OutAllInt*p%OutAllDims)) ! Need to start at 0... 
   ! Header, and Units, copy of data already available in the OutParam data structure ! TODO TODO TODO remove copy
   CALL AllocAry(InitOut%WriteOutputHdr, p%NumOuts + p%OutAllint*p%OutAllDims, 'WriteOutputHdr', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(InitOut%WriteOutputUnt, p%NumOuts + p%OutAllint*p%OutAllDims, 'WriteOutputUnt', ErrStat2, ErrMsg2); if(Failed()) return
   misc%SDWrOutput  = 0.0_ReKi
   misc%LastOutTime = 0.0_DbKi
   misc%Decimat     = 0
   y%WriteOutput = 0
   DO I = 1,p%NumOuts+p%OutAllint*p%OutAllDims
      InitOut%WriteOutputHdr(I) = TRIM( p%OutParam(I)%Name  )
      InitOut%WriteOutputUnt(I) = TRIM( p%OutParam(I)%Units )      
   END DO  
     
   !_________________________________ OUTPUT FOR REQUESTED MEMBERS _______________________________
   DO I=1,p%NMOutputs
      pLst => p%MOutLst(I) ! Alias to shorten notations
      CALL AllocAry(pLst%NodeIDs,    pLst%NoutCnt   , 'MOutLst(I)%NodeIDs', ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%ElmIDs,     pLst%NoutCnt, 2, 'MOutLst(I)%ElmIDs' , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%ElmNds,     pLst%NoutCnt, 2, 'MOutLst(I)%ElmNds' , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%Me, 12, 12, pLst%NoutCnt, 2, 'MOutLst(I)%Me'     , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%Ke, 12, 12, pLst%NoutCnt, 2, 'MOutLst(I)%Ke'     , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%Fg,     12, pLst%NoutCnt, 2, 'MOutLst(I)%Fg'     , ErrStat2, ErrMsg2); if(Failed()) return

      ! NOTE: len(MemberNodes) >2 if nDiv>1
      iMember = FINDLOCI(Init%Members(:,1), pLst%MemberID) ! Reindexing from MemberID to 1:nMembers
      pLst%NodeIDs(1:pLst%NoutCnt)=Init%MemberNodes(iMember, pLst%NodeCnt)  ! We are storing the actual node numbers corresponding to what the user ordinal number is requesting
      pLst%ElmIDs=0  !Initialize to 0
      pLst%ElmNds=0  !Initialize to 0

      DO J=1,pLst%NoutCnt ! loop on requested nodes for that member
         iNode        = pLst%NodeIDs(J)           ! Index of requested node in node list
         nElemPerNode = Init%NodesConnE(iNode, 1) ! Number of elements connecting to the j-th node
         ! Finding 1 or max 2 elements that belong to the member and connect to the node
         K2=0 ! Counter so that max 2 elements are included: NOTE: I belive more than 2 should be an error
         DO iiElem = 1, nElemPerNode
            iElem = Init%NodesConnE(iNode, iiElem+1) ! iiElem-th Element Number
            IF (ThisElementIsAlongMember(iElem, iNode, iMember)) THEN
               IF (K2 == 2) EXIT ! we found both elements already, error...
               K2=K2+1
               call ConfigOutputNode_MKF_ID(pLst, iElem, iiNode=J, iStore=K2, NodeID2=iNode)
            END IF    
         ENDDO  ! iiElem, nElemPerNode
      ENDDO !J, Noutcnt
   ENDDO  !I, NMOutputs
 
   !_________________________________ OUTPUT FOR ALL MEMBERS __________________________________
   IF (p%OutAll) THEN  !I need to store all member end forces and moments 

      ! MOutLst2: nodal output info by members, for all members, First and Last Node
      ALLOCATE ( p%MOutLst2(p%NMembers), STAT = ErrStat2 ); ErrMsg2 = 'Error allocating p%MOutLst2 array in SDOut_Init'; if(Failed()) return

      DO iMember=1,p%NMembers
         pLst => p%MOutLst2(iMember) ! Alias
         pLst%MemberID = Init%Members(iMember,1)
         nNodesPerElem = count(Init%MemberNodes(iMember,:) >0 ) 
         CALL AllocAry(pLst%NodeIDs, nNodesPerElem, 'MOutLst2(I)%NodeIDs', ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%ElmIDs,     2, 1, 'MOutLst2(I)%ElmIDs'     , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%ElmNds,     2, 1, 'MOutLst2(I)%ElmNds'     , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Me, 12, 12, 2, 1, 'MOutLst2(I)%Me'         , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Ke, 12, 12, 2, 1, 'MOutLst2(I)%Ke'         , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Fg,     12, 2, 1, 'MOutLst2(I)%Fg'         , ErrStat2, ErrMsg2); if(Failed()) return
         pLst%NodeIDs(1:nNodesPerElem) = Init%MemberNodes(iMember,1:nNodesPerElem) ! We are storing  the actual node numbers in the member
         !ElmIDs could contain the same element twice if Ndiv=1
         pLst%ElmIDs=0  !Initialize to 0
         DO J=1,nNodesPerElem,nNodesPerElem-1 ! loop on first and last node of member
            iNode        = pLst%NodeIDs(J)           ! Index of requested node in node list
            nElemPerNode = Init%NodesConnE(iNode, 1) ! Number of elements connecting to the 1st or last node of the member
            K2= J/(nNodesPerElem)+1  ! 1 (first node) or 2 (last node) depending on J
            DO iiElem=1, nElemPerNode
               iElem = Init%NodesConnE(iNode,iiElem+1) ! iiElem-th Element Number in the set of elements attached to the selected node
               IF (ThisElementIsAlongMember(iElem, iNode, iMember)) THEN
                  call ConfigOutputNode_MKF_ID(pLst, iElem, iiNode=K2, iStore=1, NodeID2=iNode)
                  EXIT   !We found the element for that node, exit loop on elements
               ENDIF
            ENDDO
         ENDDO ! Loop on divisions
      ENDDO ! Loop on members
   ENDIF ! OutAll
   !_____________________________________REACTIONS_____________________________________________
   ! --- Check if reaction requested by user
   p%OutReact = .FALSE.
   DO I=1,p%NumOuts
      if ( ANY( p%OutParam(I)%Indx == ReactSS) ) THEN ! bjj: removed check of first 5 characters being "React" because (1) cases matter and (2) we can also ask for "-React*" or "mREACT"
         p%OutReact   =.TRUE.  
         EXIT
      ENDIF
   ENDDO
   IF (p%OutReact) THEN  !I need to store all constrained forces and moments; WE do not allow more than one member to be connected at a constrained joint for the time being
      ! MOutLst3: nodal output info by members, for the members involved in reaction
      ALLOCATE(p%MOutLst3(p%nNodes_C), STAT = ErrStat2); ErrMsg2 = 'Error allocating p%MOutLst3 array in SDOut_Init'; if(Failed()) return

      DO I=1,p%nNodes_C  !For all constrained node
         pLst => p%MOutLst3(I)
         iNode        = p%Nodes_C(I,1)           ! Note: Nodes_C has been reindexed
         nElemPerNode = Init%NodesConnE(iNode,1) ! Number of elements connecting to the joint
         CALL AllocAry(pLst%ElmIDs,      1, nElemPerNode, ' p%MOutLst3(I)%ElmIds', ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%ElmNds,      1, nElemPerNode, ' p%MOutLst3(I)%ElmNds', ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Me, 12, 12 , 1, nElemPerNode, ' p%MOutLst3(I)%Me'    , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Ke, 12, 12 , 1, nElemPerNode, ' p%MOutLst3(I)%Ke'    , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Fg,     12 , 1, nElemPerNode, ' p%MOutLst3(I)%Fg'    , ErrStat2, ErrMsg2); if(Failed()) return
         DO iiElem = 1, nElemPerNode
            iElem = Init%NodesConnE(iNode, iiElem+1) ! iiElem-th Element Number in the set of elements attached to the selected node 
            call ConfigOutputNode_MKF_ID(pLst, iElem, iiNode=1, iStore=iiElem, NodeID2=iNode) 
         ENDDO
      ENDDO
      ! Compute p%TIreact, rigid transf. matrix from reaction DOFs to base structure point (0,0,-WD)
      CALL AllocAry(p%TIreact, 6, p%nDOFC__, 'TIReact  ', ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(T_TIreact, p%nDOFC__, 6, 'TIReact_T', ErrStat2, ErrMsg2); if(Failed()) return
      call RigidTrnsf(Init, p, (/0.0_Reki, 0.0_ReKi, -WtrDpth /), p%IDC__, p%nDOFC__, T_TIreact, ErrStat2, ErrMsg2); if(Failed()) return
      p%TIreact=transpose(T_TIreact)
      deallocate(T_TIreact)
   ENDIF
   RETURN

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SDOut_Init') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed

   !> Returns true if an element is connected to node iNode, and along member iMember
   LOGICAL FUNCTION ThisElementIsAlongMember(iElem, iNode, iMember)
      integer(IntKi), intent(in) :: iElem   !< Element index 
      integer(IntKi), intent(in) :: iNode   !< Node index
      integer(IntKi), intent(in) :: iMember !< Member index
      integer(IntKi), dimension(2) :: ElemNodes  ! Node IDs for element under consideration (may not be consecutive numbers)
      integer(IntKi)               :: iOtherNode ! Other node than iNode for element iElem
      ElemNodes = p%Elems(iElem,2:3) ! 1st and 2nd node of the element
      ! Check that the other node belongs to the member
      IF      (ElemNodes(1) == iNode) then
         iOtherNode=ElemNodes(2)
      else if (ElemNodes(2) == iNode) then
         iOtherNode=ElemNodes(1)
      else
         ThisElementIsAlongMember=.false. ! Not along member since nodes don't match
         return 
      endif
      ! Being along the member means the second node of the element is in the node list of the member
      ThisElementIsAlongMember= ANY(Init%MemberNodes(iMember,:) == iOtherNode)
   END FUNCTION

   !> Set different "data" for a given output node, and possibly store more than one "data" per node:
   !! The "data" is:
   !!   - Mass, stiffness matrices and constant element force (gravity and cable)
   !!   - A flag whether the node is the 1st or second node of an element 
   !! The "data" is stored at the index (iiNode,iStore):
   !!   - iiNode: node index within the list of nodes that are to be used for output for this member
   !!   - iStore: index over the number of "data" stored per node. E.g. Member1 and 2 connecting to a node  
   SUBROUTINE ConfigOutputNode_MKF_ID(pLst, iElem, iiNode, iStore, NodeID2)
      type(MeshAuxDataType), intent(inout)       :: pLst   !< Info for one member output
      integer(IntKi)       , intent(in)          :: iElem  !< Element index to which the node belong
      integer(IntKi)       , intent(in)          :: iiNode !< Index over the nodes of a given member (>2 if nDIV>1)
      integer(IntKi)       , intent(in)          :: iStore !< Storage index, used several informations are stored per node
      integer(IntKi)       , intent(in)          :: NodeID2 !< If ElemNode(2) == NodeID2, then it's the second node
      integer(IntKi), dimension(2) :: ElemNodes  ! Node IDs for element under consideration (may not be consecutive numbers)
      REAL(FEKi)                   :: FCe(12) ! Pretension force from cable element
      pLst%ElmIDs(iiNode,iStore) = iElem              ! This array has for each joint requested  the elements' ID to get results for
      ElemNodes = p%Elems(iElem,2:3) ! 1st and 2nd node of the k-th element
      if (ElemNodes(2) == NodeID2) then 
         pLst%ElmNds(iiNode,iStore) = 2 ! store whether first or second node of element  
      else
         pLst%ElmNds(iiNode,iStore) = 1 ! store whether first or second node of element
      endif
      ! --- Element Me, Ke, Fg, Fce
      CALL ElemM(p%ElemProps(iElem),         pLst%Me(:,:,iiNode,iStore))
      CALL ElemK(p%ElemProps(iElem),         pLst%Ke(:,:,iiNode,iStore))
      CALL ElemF(p%ElemProps(iElem), Init%g, pLst%Fg(:,iiNode,iStore), FCe)
      ! NOTE: Removing this force contribution for now 
      ! The output of subdyn will just be the "Kx" part for now
      !pLst%Fg(:,iiNode,iStore) = pLst%Fg(:,iiNode,iStore) + FCe(1:12) ! Adding cable element force 
      pLst%Fg(:,iiNode,iStore) = FCe(1:12) ! Adding cable element force 
   END SUBROUTINE ConfigOutputNode_MKF_ID


END SUBROUTINE SDOut_Init
!------------------------------------------------------------------------------------------------------
!> Writes the data stored in the y variable to the correct indexed postions in WriteOutput
!! This is called by SD_CalcOutput() at each time step.
!! This routine does fill Allouts
!! note that this routine assumes m%u_TP and m%udotdot_TP have been set before calling 
!!     this routine (which is done in SD_CalcOutput() and SD CalcContStateDeriv)
SUBROUTINE SDOut_MapOutputs(u,p,x, y, m, AllOuts, ErrStat, ErrMsg )
   type(SD_InputType),            intent( in )     :: u                    ! SubDyn module's input data
   type(SD_ContinuousStateType),  intent( in )     :: x                    ! SubDyn module's states data
   type(SD_OutputType),           intent( inout )  :: y                    ! SubDyn module's output data
   type(SD_ParameterType), target,intent( in    )  :: p                    ! SubDyn module's parameter data
   type(SD_MiscVarType),          intent( inout )  :: m                    ! Misc/optimization variables
   real(ReKi),                    intent(   out )  :: AllOuts(0:MaxOutPts+p%OutAllInt*p%OutAllDims) ! Array of output data for all possible outputs
   integer(IntKi),                intent(   out )  :: ErrStat              ! Error status of the operation
   character(*),                  intent(   out )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   !locals
   integer(IntKi)                 :: iMemberOutput, iiNode, iSDNode, iMeshNode, I, J, L, L2      ! Counters
   integer(IntKi)                 :: maxOutModes  ! maximum modes to output, the minimum of 99 or p%nDOFM
   real(ReKi), dimension (6)      :: FM_elm, FK_elm, Fext  !output static and dynamic forces and moments
   real(ReKi), dimension (6)      :: FM_elm2, FK_elm2      !output static and dynamic forces and moments
   real(FEKi), dimension (3,3)    :: DIRCOS    !direction cosice matrix (global to local) (3x3)
   real(ReKi), allocatable        :: ReactNs(:)    !6*Nreact reactions
   integer(IntKi)                 :: sgn !+1/-1 for node force calculations
   type(MeshAuxDataType), pointer :: pLst       !< Info for a given member-output (Alias to shorten notation)
   integer(IntKi), pointer        :: DOFList(:) !< List of DOF indices for a given Nodes (Alias to shorten notation)
   ErrStat = ErrID_None   
   ErrMsg  = ""
   
   AllOuts = 0.0_ReKi  ! initialize for those outputs that aren't valid (and thus aren't set in this routine)
         
   ! --------------------------------------------------------------------------------
   ! --- Requested member-outputs (Node kinematics and loads)
   ! --------------------------------------------------------------------------------
   ! p%MOutLst has the mapping for the member, node, elements per node, to be used
   ! MXNYZZZ   will need to connects to p%MOutLst(X)%ElmIDs(Y,1:2) if it is a force or accel; else to u%UFL(p%MOutLst(X)%NodeIDs(Y)) 
   if (p%NumOuts > 0) then  !bjj: some of these fields aren't allocated when NumOuts==0
      ! Loop over member-outputs requested
      DO iMemberOutput=1,p%NMOutputs
         pLst=>p%MOutLst(iMemberOutput) ! List for a given member-output 
         DO iiNode=1,pLst%NOutCnt !Iterate on requested nodes for that member 
            ! --- Forces (potentially averaged on 2 elements) 
            call ElementForce(pLst, iiNode, 1, FM_elm, FK_elm, sgn, DIRCOS, .false.)
            FM_elm2=sgn*FM_elm
            FK_elm2=sgn*FK_elm
            IF (pLst%ElmIDs(iiNode,2) .NE. 0) THEN  ! Second element exist
               ! NOTE: forces are computed in the coordinate system of the first element for averaging
               call ElementForce(pLst, iiNode, 2, FM_elm, FK_elm, sgn, DIRCOS, .true.) ! True= we use DIRCOS from element above
               FM_elm2=0.5*( FM_elm2 + sgn*FM_elm ) ! Now Average
               FK_elm2=0.5*( FK_elm2 + sgn*FK_elm)  ! Now Average
            ENDIF
            ! Static (elastic) component of reaction forces and moments at MαNβ along local member coordinate system
            !    "MαNβFKxe, MαNβFKye, MαNβFKze, MαNβMKxe, MαNβMKye, MαNβMKze"
            AllOuts(MNfmKe  (:,iiNode,iMemberOutput)) = FK_elm2  !static forces and moments (6) Local Ref
            ! Dynamic (inertial) component of reaction forces and moments at MαNβ along local member coordinate system
            !    "MαNβFMxe, MαNβFMye, MαNβFMze, MαNβMMxe, MαNβMMye, MαNβMMze"
            AllOuts(MNfmMe  (:,iiNode,iMemberOutput)) = FM_elm2  !dynamic forces and moments (6) Local Ref

            ! --- Displacements and acceleration
            DOFList => p%NodesDOF(pLst%NodeIDs(iiNode))%List
            ! Displacement- Translational -no need for averaging since it is a node translation - In global reference SS
            !     "MαNβTDxss, MαNβTDyss, MαNβTDzss"
            AllOuts(MNTDss (:,iiNode,iMemberOutput))       = m%U_full(DOFList(1:3))
            ! Displacement- Rotational - need direction cosine matrix to tranform rotations  - In Local reference Element Ref Sys
            !     "MαNβRDxss, MαNβRDye, MαNβRDze"
            AllOuts(MNRDe (:,iiNode,iMemberOutput))        = matmul(DIRCOS,m%U_full(DOFList(4:6))) !local ref
            ! Accelerations- I need to get the direction cosine matrix to tranform displacement and rotations
            !     "MαNβTAxe, MαNβTAye, MαNβTAze"
            !     "MαNβRAxe, MαNβRAye, MαNβRAze"
            AllOuts(MNTRAe (1:3,iiNode,iMemberOutput))     = matmul(DIRCOS,m%U_full_dotdot(DOFList(1:3))) ! translational accel local ref
            AllOuts(MNTRAe (4:6,iiNode,iMemberOutput))     = matmul(DIRCOS,m%U_full_dotdot(DOFList(4:6))) ! rotational accel  local ref
        ENDDO  ! iiNode, Loop on requested nodes for that member
     ENDDO ! iMemberOutput, Loop on member outputs
   END IF
  
   ! --------------------------------------------------------------------------------
   ! --- All nodal loads from stiffness and mass matrix 
   ! --------------------------------------------------------------------------------
   ! "MaaaJbFKxe, MaaaJbMKxe MaaaJbFMxe, MaaaJbMMxe for member aaa and node b."
   IF (p%OutAll) THEN 
      DO iMemberOutput=1,p%NMembers    !Cycle on all members
         pLst=>p%MOutLst2(iMemberOutput)
         DO iiNode=1,2 !Iterate on requested nodes for that member (first and last)  
            call ElementForce(pLst, iiNode, 1, FM_elm, FK_elm, sgn, DIRCOS, .false.)
            ! Store in All Outs
            L  = MaxOutPts+(iMemberOutput-1)*24+(iiNode-1)*12+1
            L2 = L+11
            AllOuts( L:L2 ) =sgn* (/FK_elm,FM_elm/)
         ENDDO !iiNode, nodes 1 and 2
      ENDDO ! iMemberOutput, Loop on members
   ENDIF
  
   ! --------------------------------------------------------------------------------
   ! --- Interface kinematics and loads (TP/platform reference point)
   ! --------------------------------------------------------------------------------
   ! Total interface reaction forces and moments in SS coordinate system
   !    "IntfFXss, IntfFYss, IntfFZss, IntfMXss, IntfMYss, IntfMZss,"
   AllOuts(IntfSS(1:nDOFL_TP))= - (/y%Y1Mesh%Force (:,1), y%Y1Mesh%Moment(:,1)/) !-y%Y1  !Note this is the force that the TP applies to the Jacket, opposite to what the GLue Code needs thus "-" sign
   ! Interface translations and rotations in SS coordinate system 
   !    "IntfTDXss, IntfTDYss, IntfTDZss, IntfRDXss, IntfRDYss IntfRDZss"
   AllOuts(IntfTRss(1:nDOFL_TP))=m%u_TP 
   ! Interface Translational and rotational accelerations in SS coordinate system
   !    "IntfTAXss, IntfTAYss, IntfTAZss, IntfRAXss, IntfRAYss IntfRAZss"
   AllOuts(IntfTRAss(1:nDOFL_TP))= m%udotdot_TP 

   ! --------------------------------------------------------------------------------
   ! --- Modal parameters "SSqmXX, SSqmdotXX, SSqmddXX" amplitude, speed and acceleration
   ! --------------------------------------------------------------------------------
   maxOutModes = min(p%nDOFM,99) ! We only have space for the first 99 values
   IF ( maxOutModes > 0 ) THEN 
      !BJJ: TODO: is there a check to see if we requested these channels but didn't request the modes? (i.e., retain 2 modes but asked for 75th mode?)
      AllOuts(SSqm01  :SSqm01  +maxOutModes-1) = x%qm      (1:maxOutModes)
      AllOuts(SSqmd01 :SSqmd01 +maxOutModes-1) = x%qmdot   (1:maxOutModes)
      AllOuts(SSqmdd01:SSqmdd01+maxOutModes-1) = m%qmdotdot(1:maxOutModes)
   END IF
   
   ! --------------------------------------------------------------------------------}
   ! --- Base reaction loads
   ! --------------------------------------------------------------------------------{
   ! Total base reaction forces and moments at the (0.,0.,-WtrDpth) location in SS coordinate system
   !    "ReactFXss, ReactFYss, ReactFZss, ReactMXss, ReactMYss, ReactMZss"
   IF (p%OutReact) THEN 
      ALLOCATE ( ReactNs(6*p%nNodes_C), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for ReactNs array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      ReactNs = 0.0_ReKi !Initialize
      DO I=1,p%nNodes_C   !Do for each constrained node, they are ordered as given in the input file and so as in the order of y2mesh
         FK_elm2=0._ReKi !Initialize for cumulative force
         FM_elm2=0._ReKi !Initialize
         pLst => p%MOutLst3(I)
         !Find the joint forces
         DO J=1,SIZE(pLst%ElmIDs(1,:))  !for all the elements connected (normally 1)
            iiNode = 1
            call ElementForce(pLst, iiNode, J, FM_elm, FK_elm, sgn, DIRCOS, .false.)
            !transform back to global, need to do 3 at a time since cosine matrix is 3x3
            DO L=1,2  
               FM_elm2((L-1)*3+1:L*3) = FM_elm2((L-1)*3+1:L*3) + matmul(transpose(DIRCOS),FM_elm((L-1)*3+1:L*3))  !sum forces at joint in GLOBAL REF
               FK_elm2((L-1)*3+1:L*3) = FK_elm2((L-1)*3+1:L*3) + matmul(transpose(DIRCOS),FK_elm((L-1)*3+1:L*3))  !signs may be wrong, we will fix that later;  
               ! I believe this is all fixed in terms of signs now ,RRD 5/20/13
            ENDDO           
         ENDDO
         ! FK_elm2 ! + FM_elm2  !removed the inertial component 12/13 !Not sure why I need an intermediate step here, but the sum would not work otherwise
         ! NEED TO ADD HYDRODYNAMIC FORCES AT THE RESTRAINT NODES
         iSDNode   = p%Nodes_C(I,1) 
         iMeshNode = iSDNode ! input and Y2 mesh nodes are the same as subdyn
         Fext =  (/ u%LMesh%Force(:,iMeshNode), u%LMesh%Moment(:,iMeshNode) /)
         ReactNs((I-1)*6+1:6*I) = FK_elm2 - Fext  !Accumulate reactions from all nodes in GLOBAL COORDINATES
      ENDDO
      ! Store into AllOuts
      AllOuts( ReactSS(1:nDOFL_TP) ) = matmul(p%TIreact,ReactNs)
   ENDIF
   if (allocated(ReactNs)) deallocate(ReactNs)
contains

   subroutine ElementForce(pLst, iiNode, JJ, FM_elm, FK_elm, sgn, DIRCOS, bUseInputDirCos)
      type(MeshAuxDataType), intent(in)          :: pLst   !< Info for one member output
      integer(IntKi)       , intent(in)          :: iiNode !< Index over the nodes of a given member (>2 if nDIV>1)
      integer(IntKi)       , intent(in)          :: JJ     !< TODO: interpretation: index over other member connected to the current member (for averaging)
      real(FEKi), dimension (3,3), intent(inout) :: DIRCOS  !direction cosice matrix (global to local) (3x3)
      real(ReKi), dimension (6), intent(out)     :: FM_elm, FK_elm  !output static and dynamic forces and moments
      integer(IntKi), intent(out)                :: sgn !+1/-1 for node force calculations
      logical, intent(in)                        :: bUseInputDirCos !< If True, use DIRCOS from input, otherwise, use element DirCos
      ! Local
      integer(IntKi)                          :: iElem !< Element index/number
      integer(IntKi)                          :: FirstOrSecond !< 1 or 2  if first node or second node
      integer(IntKi), dimension(2)            :: ElemNodes  ! Node IDs for element under consideration (may not be consecutive numbers)
      real(ReKi)    , dimension(12)           :: X_e, Xdd_e ! Displacement and acceleration for an element
      integer(IntKi), dimension(2), parameter :: NodeNumber_To_Sign = (/-1, +1/)

      iElem         = pLst%ElmIDs(iiNode,JJ)             ! element number
      FirstOrSecond = pLst%ElmNds(iiNode,JJ)             ! first or second node of the element to be considered
      sgn           = NodeNumber_To_Sign(FirstOrSecond) ! Assign sign depending if it's the 1st or second node
      ElemNodes     = p%Elems(iElem,2:3)                ! first and second node ID associated with element iElem
      X_e(1:6)      = m%U_full_elast (p%NodesDOF(ElemNodes(1))%List(1:6)) 
      X_e(7:12)     = m%U_full_elast (p%NodesDOF(ElemNodes(2))%List(1:6)) 
      Xdd_e(1:6)    = m%U_full_dotdot(p%NodesDOF(ElemNodes(1))%List(1:6)) 
      Xdd_e(7:12)   = m%U_full_dotdot(p%NodesDOF(ElemNodes(2))%List(1:6)) 
      if (.not. bUseInputDirCos) then
         DIRCOS=transpose(p%ElemProps(iElem)%DirCos)! global to local
      endif
      CALL CALC_NODE_FORCES( DIRCOS, pLst%Me(:,:,iiNode,JJ),pLst%Ke(:,:,iiNode,JJ), Xdd_e, X_e, pLst%Fg(:,iiNode,JJ), FirstOrSecond, FM_elm, FK_elm) 
   end subroutine ElementForce

   !====================================================================================================
   !> Calculates static and dynamic forces for a given element, using K and M of the element, and 
   !output quantities Udotdot and Y2 containing the 
   !and K2 indicating wheter the 1st (1) or 2nd (2) node is to be picked
   !----------------------------------------------------------------------------------------------------
   SUBROUTINE CALC_NODE_FORCES(DIRCOS,Me,Ke,Udotdot,Y2 ,Fg, FirstOrSecond, FM_nod, FK_nod)
      Real(FEKi), DIMENSION (3,3),   INTENT(IN)  :: DIRCOS    !direction cosice matrix (global to local) (3x3)
      Real(FEKi), DIMENSION (12,12), INTENT(IN)  :: Me,Ke    !element M and K matrices (12x12) in GLOBAL REFERENCE (DIRCOS^T K DIRCOS)
      Real(ReKi), DIMENSION (12),    INTENT(IN)  :: Udotdot, Y2     !acceleration and velocities, gravity forces
      Real(FEKi), DIMENSION (12),    INTENT(IN)  :: Fg     !acceleration and velocities, gravity forces
      Integer(IntKi),                INTENT(IN)  :: FirstOrSecond !1 or 2 depending on node of interest
      REAL(ReKi), DIMENSION (6),    INTENT(OUT)  :: FM_nod, FK_nod  !output static and dynamic forces and moments
      !Locals
      INTEGER(IntKi) :: L !counter
      REAL(DbKi), DIMENSION(12)                    :: FM_glb, FF_glb, FM_elm, FF_elm  ! temporary storage 

      FM_glb = matmul(Me,Udotdot)   ! GLOBAL REFERENCE
      FF_glb = matmul(Ke,Y2)        ! GLOBAL REFERENCE
      FF_glb = FF_glb - Fg          ! GLOBAL REFERENCE ! NOTE: Fg is now 0, only the "Kx" part in Fk
      DO L=1,4 ! Transforming coordinates 3 at a time
         FM_elm((L-1)*3+1:L*3) =  matmul(DIRCOS, FM_glb( (L-1)*3+1:L*3 ) )
         FF_elm((L-1)*3+1:L*3) =  matmul(DIRCOS, FF_glb( (L-1)*3+1:L*3 ) ) 
      ENDDO
      FM_nod = FM_elm(6*(FirstOrSecond-1)+1:FirstOrSecond*6) ! k2=1, 1:6,  k2=2  7:12 
      FK_nod = FF_elm(6*(FirstOrSecond-1)+1:FirstOrSecond*6) 

   END SUBROUTINE CALC_NODE_FORCES 
END SUBROUTINE SDOut_MapOutputs


!====================================================================================================
SUBROUTINE SDOut_CloseSum( UnSum, ErrStat, ErrMsg )
   INTEGER,                 INTENT( IN    )   :: UnSum                ! the unit number for the SubDyn summary file          
   INTEGER,                 INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),            INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER                                     :: Stat                 ! status from I/) operation 
   ErrStat = ErrID_None         
   ErrMsg  = ""
   ! Write any closing information in the summary file
   IF ( UnSum > 0 ) THEN
      WRITE (UnSum,'(/,A/)', IOSTAT=Stat)  '#This summary file was closed on '//CurDate()//' at '//CurTime()//'.'
      IF (Stat /= 0) THEN
         ErrStat = ErrID_FATAL
         ErrMsg  = ' Problem writing to summary file.'
      END IF
      ! Close the file
      CLOSE( UnSum, IOSTAT=Stat )
      IF (Stat /= 0) THEN
         ErrStat = ErrID_FATAL
         ErrMsg  = TRIM(ErrMsg)//' Problem closing summary file.'
      END IF
      IF ( ErrStat /= ErrID_None ) ErrMsg = 'SDOut_CloseSum'//TRIM(ErrMsg)
   END IF                      
END SUBROUTINE SDOut_CloseSum            

!====================================================================================================
SUBROUTINE SDOut_OpenSum( UnSum, SummaryName, SD_Prog, ErrStat, ErrMsg )
   INTEGER,                 INTENT(   OUT )   :: UnSum                ! the unit number for the SubDyn summary file          
   CHARACTER(*),            INTENT( IN    )   :: SummaryName          ! the name of the SubDyn summary file
   TYPE(ProgDesc),          INTENT( IN    )   :: SD_Prog              ! the name/version/date of the  program
   INTEGER,                 INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),            INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   integer                                    :: ErrStat2
   ErrStat = ErrID_None         
   ErrMsg  = ""       

   CALL GetNewUnit( UnSum )
   CALL OpenFOutFile ( UnSum, SummaryName, ErrStat, ErrMsg ) 
   IF ( ErrStat >= AbortErrLev ) THEN
      ErrMsg  = 'Failed to open SubDyn summary file: '//TRIM(ErrMsg)
      RETURN
   END IF
      
   ! Write the summary file header
   WRITE (UnSum,'(/,A/)', IOSTAT=ErrStat2)  '#This summary file was generated by '//TRIM( SD_Prog%Name )//&
                     ' '//TRIM( SD_Prog%Ver )//' on '//CurDate()//' at '//CurTime()//'.'
END SUBROUTINE SDOut_OpenSum 

!====================================================================================================
SUBROUTINE SDOut_OpenOutput( ProgVer, OutRootName,  p, InitOut, ErrStat, ErrMsg )
! This subroutine initialized the output module, checking if the output parameter list (OutList)
! contains valid names, and opening the output file if there are any requested outputs
!----------------------------------------------------------------------------------------------------
   ! Passed variables
   TYPE(ProgDesc),                INTENT( IN    ) :: ProgVer
   CHARACTER(*),                  INTENT( IN    ) :: OutRootName          ! Root name for the output file
   TYPE(SD_ParameterType),        INTENT( INOUT ) :: p   
   TYPE(SD_InitOutPutType ),      INTENT( IN    ) :: InitOut              !
   INTEGER,                       INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred           
   CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER                                        :: I                    ! Generic loop counter      
   CHARACTER(1024)                                :: OutFileName          ! The name of the output file  including the full path.
   CHARACTER(200)                                 :: Frmt                 ! a string to hold a format statement
   INTEGER                                        :: ErrStat2              
   ErrStat = ErrID_None  
   ErrMsg  = ""
   ! Open the output file, if necessary, and write the header
   IF ( ALLOCATED( p%OutParam ) .AND. p%NumOuts > 0 ) THEN           ! Output has been requested so let's open an output file            
      ! Open the file for output
      OutFileName = TRIM(OutRootName)//'.out'
      CALL GetNewUnit( p%UnJckF )
   
      CALL OpenFOutFile ( p%UnJckF, OutFileName, ErrStat, ErrMsg ) 
      IF ( ErrStat >= AbortErrLev ) THEN
         ErrMsg = ' Error opening SubDyn-level output file: '//TRIM(ErrMsg)
         RETURN
      END IF
       
      ! Write the output file header
      WRITE (p%UnJckF,'(/,A/)', IOSTAT=ErrStat2)  'These predictions were generated by '//TRIM(GETNVD(ProgVer))//&
                      ' on '//CurDate()//' at '//CurTime()//'.'
      
      WRITE(p%UnJckF, '(//)') ! add 3 lines to make file format consistant with FAST v8 (headers on line 7; units on line 8) [this allows easier post-processing]
      
      ! Write the names of the output parameters:
      Frmt = '(A8,'//TRIM(Int2LStr(p%NumOuts+p%OutAllInt*p%OutAllDims))//'(:,A,'//TRIM( p%OutSFmt )//'))'
      WRITE(p%UnJckF,Frmt, IOSTAT=ErrStat2)  TRIM( 'Time' ), ( p%Delim, TRIM( InitOut%WriteOutputHdr(I) ), I=1,p%NumOuts+p%OutAllInt*p%OutAllDims )
      
      ! Write the units of the output parameters:                 
      WRITE(p%UnJckF,Frmt, IOSTAT=ErrStat2)  TRIM( 's'), ( p%Delim, TRIM( InitOut%WriteOutputUnt(I) ), I=1,p%NumOuts+p%OutAllInt*p%OutAllDims )
   END IF   ! there are any requested outputs   
END SUBROUTINE SDOut_OpenOutput

!====================================================================================================


!====================================================================================================
SUBROUTINE SDOut_CloseOutput ( p, ErrStat, ErrMsg )
! This function cleans up after running the SubDyn output module. It closes the output file,
! releases memory, and resets the number of outputs requested to 0.
!----------------------------------------------------------------------------------------------------
   TYPE(SD_ParameterType),  INTENT( INOUT )       :: p                    ! data for this instance of the floating platform module
   INTEGER,                       INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred
   CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   LOGICAL                               :: Err

   ErrStat = 0
   ErrMsg  = ""
   Err     = .FALSE.

   ! Close our output file
   CLOSE( p%UnJckF, IOSTAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
 
   ! Make sure ErrStat is non-zero if an error occurred
   IF ( Err ) ErrStat = ErrID_Fatal
   RETURN

END SUBROUTINE SDOut_CloseOutput
!====================================================================================================

SUBROUTINE SDOut_WriteOutputNames( UnJckF, p, ErrStat, ErrMsg )

   INTEGER,                      INTENT( IN    ) :: UnJckF            ! file unit for the output file
   TYPE(SD_ParameterType),  INTENT( IN    ) :: p                    ! SubDyn module's parameter data
   INTEGER,                      INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                 INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   CHARACTER(200)                         :: Frmt                        ! a string to hold a format statement
   INTEGER                                :: I                           ! Generic loop counter
   
   ErrStat = ErrID_None   
   ErrMsg  = ""
   
   Frmt = '(A8,'//TRIM(Int2LStr(p%NumOuts+p%OutAllInt*p%OutAllDims))//'(:,A,'//TRIM( p%OutSFmt )//'))'

   WRITE(UnJckF,Frmt)  TRIM( p%OutParam(0)%Name ), ( p%Delim, TRIM( p%OutParam(I)%Name ), I=1,p%NumOuts+p%OutAllInt*p%OutAllDims )
      
END SUBROUTINE SDOut_WriteOutputNames

!====================================================================================================

SUBROUTINE SDOut_WriteOutputUnits( UnJckF, p, ErrStat, ErrMsg )
   INTEGER,                      INTENT( IN    ) :: UnJckF            ! file unit for the output file
   TYPE(SD_ParameterType),  INTENT( IN    ) :: p                    ! SubDyn module's parameter data
   INTEGER,                      INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                 INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   CHARACTER(200)                         :: Frmt                        ! a string to hold a format statement
   INTEGER                                :: I                           ! Generic loop counter
   ErrStat = ErrID_None   
   ErrMsg  = ""
   
   Frmt = '(A8,'//TRIM(Int2LStr(p%NumOuts+p%OutAllInt*p%OutAllDims))//'(:,A,'//TRIM( p%OutSFmt )//'))'

   WRITE(UnJckF,Frmt)  TRIM( p%OutParam(0)%Units ), ( p%Delim, TRIM( p%OutParam(I)%Units ), I=1,p%NumOuts+p%OutAllInt*p%OutAllDims )
      
END SUBROUTINE SDOut_WriteOutputUnits

!====================================================================================================
SUBROUTINE SDOut_WriteOutputs( UnJckF, Time, SDWrOutput, p, ErrStat, ErrMsg )
! This subroutine writes the data stored in WriteOutputs (and indexed in OutParam) to the file
! opened in SDOut_Init()
!---------------------------------------------------------------------------------------------------- 
   INTEGER,                      INTENT( IN    ) :: UnJckF               ! file unit for the output file
   REAL(DbKi),                   INTENT( IN    ) :: Time                 ! Time for this output
   REAL(ReKi),                   INTENT( IN    ) :: SDWrOutput(:)        ! SubDyn module's output data
   TYPE(SD_ParameterType),       INTENT( IN    ) :: p                    ! SubDyn module's parameter data
   INTEGER,                      INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                 INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER                                :: I                           ! Generic loop counter
   CHARACTER(200)                         :: Frmt                        ! a string to hold a format statement
   ErrStat = ErrID_None   
   ErrMsg  = ""
   
      ! Initialize ErrStat and determine if it makes any sense to write output
   IF ( .NOT. ALLOCATED( p%OutParam ) .OR. UnJckF < 0 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = ' To write outputs for SubDyn there must be a valid file ID and OutParam must be allocated.'
      RETURN
   ELSE
      ErrStat = ErrID_None
   END IF

      ! Write the output parameters to the file
   Frmt = '(F10.4,'//TRIM(Int2LStr(p%NumOuts+p%OutAllInt*p%OutAllDims))//'(:,A,'//TRIM( p%OutFmt )//'))'

   WRITE(UnJckF,Frmt)  Time, ( p%Delim, SDWrOutput(I), I=1,p%NumOuts+p%OutAllInt*p%OutAllDims )

END SUBROUTINE SDOut_WriteOutputs

!====================================================================================================


!====================================================================================================
SUBROUTINE SDOut_ChkOutLst( OutList, p, ErrStat, ErrMsg )
! This routine checks the names of inputted output channels, checks to see if any of them are ill-
! conditioned (returning an error if so), and assigns the OutputDataType settings (i.e, the index,  
! name, and units of the output channels). 
! NOTE OutParam is populated here
!----------------------------------------------------------------------------------------------------    
   TYPE(SD_ParameterType),   INTENT( INOUT ) :: p                    ! SubDyn module parameter data
   CHARACTER(ChanLen),       INTENT( IN    ) :: OutList (:)          ! An array holding the names of the requested output channels.         
   INTEGER,                  INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred           
   CHARACTER(*),             INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables.
   INTEGER                                   :: I,J,K                                         ! Generic loop-counting index.
   INTEGER                                   :: INDX                                      ! Index for valid arrays
   CHARACTER(ChanLen)                        :: OutListTmp                                ! A string to temporarily hold OutList(I).
   !CHARACTER(28), PARAMETER               :: OutPFmt    = "( I4, 3X,A 10,1 X, A10 )"   ! Output format parameter output list.
   CHARACTER(ChanLen), DIMENSION(24)         :: ToTUnits,ToTNames,ToTNames0
   LOGICAL                  :: InvalidOutput(0:MaxOutPts)                        ! This array determines if the output channel is valid for this configuration
   LOGICAL                  :: CheckOutListAgain
   ErrStat = ErrID_None   
   ErrMsg  = ""
   
   InvalidOutput            = .FALSE.

      ! mark invalid output channels:
   DO k=p%nDOFM+1,99
      InvalidOutput(SSqm01  +k-1) = .true.
      InvalidOutput(SSqmd01 +k-1) = .true.
      InvalidOutput(SSqmdd01+k-1) = .true.
   END DO
         
   DO I=1,9
          !I know el # and whether it is 1st node or second node
      if (I <= p%NMOutputs) then
         INDX=p%MOutLst(I)%NOutCnt+1
      else
         INDX = 1
      end if
            
      DO J=INDX,9 !Iterate on requested nodes for that member 
         !Forces and moments
         InvalidOutput(MNfmKe  (:,J,I)) = .true.  !static forces and moments (6) Local Ref
         InvalidOutput(MNfmMe  (:,J,I)) = .true.  !dynamic forces and moments (6) Local Ref
         !Displacement
         InvalidOutput(MNTDss  (:,J,I)) = .true.  !Translational
         InvalidOutput(MNRDe   (:,J,I)) = .true.  !Rotational
         !Accelerations
         InvalidOutput(MNTRAe  (:,J,I)) = .true.  !translational accel local ref
      END DO
   END DO
  
   !-------------------------------------------------------------------------------------------------
   ! ALLOCATE the OutParam array
   !-------------------------------------------------------------------------------------------------    
   ALLOCATE ( p%OutParam(1:p%NumOuts+p%OutAllInt*p%OutAllDims) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      ErrMsg  = ' Error allocating memory for the OutParam array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
     
   
   !-------------------------------------------------------------------------------------------------
   ! Set index, name, and units for the output channels
   ! If a selected output channel is not available in this module, set error flag and return.
   !-------------------------------------------------------------------------------------------------   
   !!!p%OutParam(0)%Name  = 'Time'    ! OutData(0) is the time channel by default.
   !!!p%OutParam(0)%Units = '(sec)'   !
   !!!p%OutParam(0)%Indx  = Time
   !!!p%OutParam(0)%SignM = 1
     
   DO I = 1,p%NumOuts
   
      p%OutParam(I)%Name = OutList(I)   
      OutListTmp         = OutList(I)
   
   
      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !   channel name with a '-', '_', 'm', or 'M' character indicating "minus".
      
      CheckOutListAgain = .FALSE.
      
      IF      ( INDEX( '-_', OutListTmp(1:1) ) > 0 ) THEN
         p%OutParam(I)%SignM = -1     ! ex, '-TipDxc1' causes the sign of TipDxc1 to be switched.
         OutListTmp                   = OutListTmp(2:)
      ELSE IF ( INDEX( 'mM', OutListTmp(1:1) ) > 0 ) THEN ! We'll assume this is a variable name for now, (if not, we will check later if OutListTmp(2:) is also a variable name)
         CheckOutListAgain            = .TRUE.
         p%OutParam(I)%SignM = 1
      ELSE
         p%OutParam(I)%SignM = 1
      END IF
      
      CALL Conv2UC( OutListTmp )    ! Convert OutListTmp to upper case
   
   
      Indx =  IndexCharAry( OutListTmp(1:9), ValidParamAry )
      
      IF ( CheckOutListAgain .AND. Indx < 1 ) THEN    ! Let's assume that "M" really meant "minus" and then test again         
         p%OutParam(I)%SignM = -1            ! ex, 'MTipDxc1' causes the sign of TipDxc1 to be switched.
         OutListTmp                   = OutListTmp(2:)
         
         Indx = IndexCharAry( OutListTmp(1:9), ValidParamAry )         
      END IF
      
      IF ( Indx > 0 ) THEN
         p%OutParam(I)%Indx = ParamIndxAry(Indx)
         IF ( InvalidOutput( ParamIndxAry(Indx) ) ) THEN
            p%OutParam(I)%Units = 'INVALID' 
            p%OutParam(I)%SignM =  0           
         ELSE
            p%OutParam(I)%Units = ParamUnitsAry(Indx)
         END IF
      ELSE
         ErrMsg  = p%OutParam(I)%Name//' is not an available output channel.'
         ErrStat = ErrID_Fatal
         p%OutParam(I)%Units = 'INVALID'  
         p%OutParam(I)%Indx  =  0
         p%OutParam(I)%SignM =  0                              ! this will print all zeros
      END IF
      
   END DO
   
   IF (p%OutAll) THEN   !Finish populating the OutParam with all the joint forces and moments
       ToTNames0=RESHAPE(SPREAD( (/"FKxe", "FKye", "FKze", "MKxe", "MKye", "MKze", "FMxe", "FMye", "FMze", "MMxe", "MMye", "MMze"/), 2, 2), (/24/) )
       ToTUnits=RESHAPE(SPREAD( (/"(N)  ","(N)  ","(N)  ", "(N*m)","(N*m)","(N*m)", "(N)  ","(N)  ","(N)  ", "(N*m)","(N*m)","(N*m)"/), 2, 2), (/24/) )
       DO I=1,p%NMembers
           DO K=1,2
            DO J=1,12  
             TotNames(J+(K-1)*12)=TRIM("M"//Int2Lstr(I))//TRIM("J"//Int2Lstr(K))//TRIM(TotNames0(J))
            ENDDO  
           ENDDO
           p%OutParam(p%NumOuts+(I-1)*12*2+1:p%NumOuts+I*12*2)%Name  = ToTNames
           p%OutParam(p%NumOuts+(I-1)*12*2+1:p%NumOuts+I*12*2)%Units = ToTUnits
       ENDDO
       p%OutParam(p%NumOuts+1:p%NumOuts+p%OutAllDims)%SignM = 1
       p%OutParam(p%NumOuts+1:p%NumOuts+p%OutAllDims)%Indx= MaxOutPts+(/(J, J=1, p%OutAllDims)/) 
   ENDIF

END SUBROUTINE SDOut_ChkOutLst
!====================================================================================================
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine initializes the array that maps rows/columns of the Jacobian to specific mesh fields.
!! Do not change the order of this packing without changing subroutine !
SUBROUTINE SD_Init_Jacobian(Init, p, u, y, InitOut, ErrStat, ErrMsg)
   TYPE(SD_InitType)                 , INTENT(IN   ) :: Init                  !< Init
   TYPE(SD_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(SD_InputType)                , INTENT(IN   ) :: u                     !< inputs
   TYPE(SD_OutputType)               , INTENT(IN   ) :: y                     !< outputs
   TYPE(SD_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Initialization output data (for Jacobian row/column names)
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'SD_Init_Jacobian'
   real(ReKi) :: dx, dy, dz, maxDim
   ! local variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! --- System dimension
   dx = maxval(Init%Nodes(:,2))- minval(Init%Nodes(:,2))
   dy = maxval(Init%Nodes(:,3))- minval(Init%Nodes(:,3))
   dz = maxval(Init%Nodes(:,4))- minval(Init%Nodes(:,4))
   maxDim = max(dx, dy, dz)
   
   ! --- System dimension
   call Init_Jacobian_y(); if (Failed()) return
   call Init_Jacobian_x(); if (Failed()) return
   call Init_Jacobian_u(); if (Failed()) return

contains
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Init_Jacobian') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed
   !> This routine initializes the Jacobian parameters and initialization outputs for the linearized outputs.

   SUBROUTINE Init_Jacobian_y()
      INTEGER(IntKi) :: index_next, i
      ! Number of outputs
      p%Jac_ny = y%Y1Mesh%nNodes * 6     & ! 3 forces + 3 moments at each node
               + y%Y2Mesh%nNodes * 18    & ! 6 displacements + 6 velocities + 6 accelerations at each node
               + y%Y3Mesh%nNodes * 18    & ! 6 displacements + 6 velocities + 6 accelerations at each node
               + p%NumOuts                 ! WriteOutput values 
      ! Storage info for each output (names, rotframe)
      call AllocAry(InitOut%LinNames_y, p%Jac_ny, 'LinNames_y',ErrStat2,ErrMsg2); if(ErrStat2/=ErrID_None) return
      call AllocAry(InitOut%RotFrame_y, p%Jac_ny, 'RotFrame_y',ErrStat2,ErrMsg2); if(ErrStat2/=ErrID_None) return
      ! Names
      index_next = 1
      call PackLoadMesh_Names(  y%Y1Mesh, 'Interface displacement', InitOut%LinNames_y, index_next)
      call PackMotionMesh_Names(y%Y2Mesh, 'Nodes motion mixed'    , InitOut%LinNames_y, index_next)
      call PackMotionMesh_Names(y%Y3Mesh, 'Nodes motion full'     , InitOut%LinNames_y, index_next)
      do i=1,p%NumOuts
         InitOut%LinNames_y(i+index_next-1) = trim(InitOut%WriteOutputHdr(i))//', '//trim(InitOut%WriteOutputUnt(i))
      end do
      ! RotFrame
      InitOut%RotFrame_y(:) = .false.
   END SUBROUTINE Init_Jacobian_y

   !> This routine initializes the Jacobian parameters and initialization outputs for the linearized continuous states.
   SUBROUTINE Init_Jacobian_x()
      INTEGER(IntKi) :: i
      p%Jac_nx = p%nDOFM ! qm 
      ! allocate space for the row/column names and for perturbation sizes
      CALL AllocAry(InitOut%LinNames_x  , 2*p%Jac_nx, 'LinNames_x'  , ErrStat2, ErrMsg2); if(ErrStat/=ErrID_None) return
      CALL AllocAry(InitOut%RotFrame_x  , 2*p%Jac_nx, 'RotFrame_x'  , ErrStat2, ErrMsg2); if(ErrStat/=ErrID_None) return
      CALL AllocAry(InitOut%DerivOrder_x, 2*p%Jac_nx, 'DerivOrder_x', ErrStat2, ErrMsg2); if(ErrStat/=ErrID_None) return
      ! default perturbations, p%dx:
      p%dx(1) = 2.0_ReKi*D2R_D   ! deflection states in rad and rad/s
      p%dx(2) = 2.0_ReKi*D2R_D   ! deflection states in rad and rad/s
      InitOut%RotFrame_x   = .false.
      InitOut%DerivOrder_x = 2
      ! set linearization output names:
      do i=1,p%Jac_nx
         InitOut%LinNames_x(i) = 'Craig-Bampton mode '//trim(num2lstr(i))//' amplitude, -'; 
      end do
      do i=1,p%Jac_nx
         InitOut%LinNames_x(i+p%Jac_nx) = 'First time derivative of '//trim(InitOut%LinNames_x(i))//'/s'
         InitOut%RotFrame_x(i+p%Jac_nx) = InitOut%RotFrame_x(i)
      end do
   END SUBROUTINE Init_Jacobian_x

   SUBROUTINE Init_Jacobian_u()
      REAL(R8Ki)     :: perturb
      INTEGER(IntKi) :: i, j, idx, nu, i_meshField
      ! Number of inputs
      nu = u%TPMesh%nNodes * 18 & ! 3 Translation Displacements + 3 orientations + 6 velocities + 6 accelerations at each node
         + u%LMesh%nNodes  * 6    ! 3 forces + 3 moments at each node
      ! --- Info of linearized inputs (Names, RotFrame, IsLoad)
      call AllocAry(InitOut%LinNames_u, nu, 'LinNames_u', ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return
      call AllocAry(InitOut%RotFrame_u, nu, 'RotFrame_u', ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return
      call AllocAry(InitOut%IsLoad_u  , nu, 'IsLoad_u'  , ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return
      InitOut%RotFrame_u = .false. ! every input is on a mesh, which stores values in the global (not rotating) frame
      idx = 1
      call PackMotionMesh_Names(u%TPMesh, 'TPMesh', InitOut%LinNames_u, idx) ! all 6 motion fields
      InitOut%IsLoad_u(1:idx-1) = .false. ! the TPMesh inputs are not loads
      InitOut%IsLoad_u(idx:)    = .true.  ! the remaining inputs are loads
      call PackLoadMesh_Names(  u%LMesh,   'LMesh', InitOut%LinNames_u, idx)

      ! --- Jac_u_indx:  matrix to store index to help us figure out what the ith value of the u vector really means
      ! (see perturb_u ... these MUST match )
      ! column 1 indicates module's mesh and field
      ! column 2 indicates the first index (x-y-z component) of the field
      ! column 3 is the node
      call allocAry( p%Jac_u_indx, nu, 3, 'p%Jac_u_indx', ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return
      idx = 1
      !Module/Mesh/Field: u%TPMesh%TranslationDisp  = 1;
      !Module/Mesh/Field: u%TPMesh%Orientation      = 2;
      !Module/Mesh/Field: u%TPMesh%TranslationVel   = 3;
      !Module/Mesh/Field: u%TPMesh%RotationVel      = 4;
      !Module/Mesh/Field: u%TPMesh%TranslationAcc   = 5;
      !Module/Mesh/Field: u%TPMesh%RotationAcc      = 6;
      do i_meshField = 1,6
         do i=1,u%TPMesh%nNodes
            do j=1,3
               p%Jac_u_indx(idx,1) =  i_meshField
               p%Jac_u_indx(idx,2) =  j !component idx:  j
               p%Jac_u_indx(idx,3) =  i !Node:   i
               idx = idx + 1
            end do !j
         end do !i
      end do
      !Module/Mesh/Field: u%LMesh%Force   = 7;
      !Module/Mesh/Field: u%LMesh%Moment  = 8;
      do i_meshField = 7,8
         do i=1,u%LMesh%nNodes
            do j=1,3
               p%Jac_u_indx(idx,1) =  i_meshField
               p%Jac_u_indx(idx,2) =  j !component idx:  j
               p%Jac_u_indx(idx,3) =  i !Node:   i
               idx = idx + 1
            end do !j
         end do !i
      end do

      ! --- Default perturbations, p%du:
      call allocAry( p%du, 8, 'p%du', ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return ! 8 = number of unique values in p%Jac_u_indx(:,1)
      perturb   = 2.0_R8Ki*D2R_D
      p%du( 1) = perturb       ! u%TPMesh%TranslationDisp  = 1;
      p%du( 2) = perturb       ! u%TPMesh%Orientation      = 2;
      p%du( 3) = perturb       ! u%TPMesh%TranslationVel   = 3;
      p%du( 4) = perturb       ! u%TPMesh%RotationVel      = 4;
      p%du( 5) = perturb       ! u%TPMesh%TranslationAcc   = 5;
      p%du( 6) = perturb       ! u%TPMesh%RotationAcc      = 6;
      p%du( 7) = 170*maxDim**2 ! u%LMesh%Force             = 7;
      p%du( 8) =  14*maxDim**3 ! u%LMesh%Moment            = 8;
   END SUBROUTINE Init_Jacobian_u

END SUBROUTINE SD_Init_Jacobian
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the u array (and mesh/field it corresponds to)
!! Do not change this without making sure subroutine beamdyn::init_jacobian is consistant with this routine!
SUBROUTINE SD_Perturb_u( p, n, perturb_sign, u, du )
   TYPE(SD_ParameterType)              , INTENT(IN   ) :: p            !< parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: n            !< number of array element to use
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(SD_InputType)                  , INTENT(INOUT) :: u            !< perturbed SD inputs
   REAL( R8Ki )                        , INTENT(  OUT) :: du           !< amount that specific input was perturbed
   ! local variables
   INTEGER :: fieldIndx
   INTEGER :: node
   fieldIndx = p%Jac_u_indx(n,2)
   node      = p%Jac_u_indx(n,3)
   du = p%du(  p%Jac_u_indx(n,1) )
   ! determine which mesh we're trying to perturb and perturb the input:
   SELECT CASE( p%Jac_u_indx(n,1) )
   CASE ( 1) !Module/Mesh/Field: u%TPMesh%TranslationDisp = 1;
      u%TPMesh%TranslationDisp( fieldIndx,node) = u%TPMesh%TranslationDisp( fieldIndx,node) + du * perturb_sign
   CASE ( 2) !Module/Mesh/Field: u%TPMesh%Orientation = 2;
      CALL PerturbOrientationMatrix( u%TPMesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
   CASE ( 3) !Module/Mesh/Field: u%TPMesh%TranslationVel = 3;
      u%TPMesh%TranslationVel( fieldIndx,node) = u%TPMesh%TranslationVel( fieldIndx,node) + du * perturb_sign
   CASE ( 4) !Module/Mesh/Field: u%TPMesh%RotationVel = 4;
      u%TPMesh%RotationVel(fieldIndx,node) = u%TPMesh%RotationVel(fieldIndx,node) + du * perturb_sign
   CASE ( 5) !Module/Mesh/Field: u%TPMesh%TranslationAcc = 5;
      u%TPMesh%TranslationAcc( fieldIndx,node) = u%TPMesh%TranslationAcc( fieldIndx,node) + du * perturb_sign
   CASE ( 6) !Module/Mesh/Field: u%TPMesh%RotationAcc = 6;
      u%TPMesh%RotationAcc(fieldIndx,node) = u%TPMesh%RotationAcc(fieldIndx,node) + du * perturb_sign
   CASE ( 7) !Module/Mesh/Field: u%LMesh%Force = 7;
      u%LMesh%Force(fieldIndx,node) = u%LMesh%Force(fieldIndx,node) + du * perturb_sign 
   CASE ( 8) !Module/Mesh/Field: u%LMesh%Moment  = 8;
      u%LMesh%Moment(fieldIndx,node) = u%LMesh%Moment(fieldIndx,node) + du * perturb_sign
   END SELECT
END SUBROUTINE SD_Perturb_u
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine beamdyn::init_jacobian is consistant with this routine!
SUBROUTINE SD_Compute_dY(p, y_p, y_m, delta, dY)
   TYPE(SD_ParameterType)            , INTENT(IN   ) :: p     !< parameters
   TYPE(SD_OutputType)               , INTENT(IN   ) :: y_p   !< SD outputs at \f$ u + \Delta_p u \f$ or \f$ z + \Delta_p z \f$ (p=plus)
   TYPE(SD_OutputType)               , INTENT(IN   ) :: y_m   !< SD outputs at \f$ u - \Delta_m u \f$ or \f$ z - \Delta_m z \f$ (m=minus)
   REAL(R8Ki)                        , INTENT(IN   ) :: delta !< difference in inputs or states \f$ delta_p = \Delta_p u \f$ or \f$ delta_p = \Delta_p x \f$
   REAL(R8Ki)                        , INTENT(INOUT) :: dY(:) !< column of dYdu or dYdx: \f$ \frac{\partial Y}{\partial u_i} = \frac{y_p - y_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial Y}{\partial z_i} = \frac{y_p - y_m}{2 \, \Delta x}\f$
   ! local variables:
   INTEGER(IntKi) :: i              ! loop over outputs
   INTEGER(IntKi) :: indx_first     ! index indicating next value of dY to be filled
   indx_first = 1
   call PackLoadMesh_dY(  y_p%Y1Mesh, y_m%Y1Mesh, dY, indx_first)
   call PackMotionMesh_dY(y_p%Y2Mesh, y_m%Y2Mesh, dY, indx_first, UseSmlAngle=.true.) ! all 6 motion fields
   call PackMotionMesh_dY(y_p%Y3Mesh, y_m%Y3Mesh, dY, indx_first, UseSmlAngle=.true.) ! all 6 motion fields
   do i=1,p%NumOuts
      dY(i+indx_first-1) = y_p%WriteOutput(i) - y_m%WriteOutput(i)
   end do
   dY = dY / (2.0_R8Ki*delta)
END SUBROUTINE SD_Compute_dY
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the x array (and mesh/field it corresponds to)
!! Do not change this without making sure subroutine sd_init_jacobian is consistant with this routine!
SUBROUTINE SD_Perturb_x( p, fieldIndx, mode, perturb_sign, x, dx )
   TYPE(SD_ParameterType)      , INTENT(IN   ) :: p            !< parameters
   INTEGER( IntKi )            , INTENT(IN   ) :: fieldIndx    !< field in the state type: 1=displacements; 2=velocities
   INTEGER( IntKi )            , INTENT(IN   ) :: mode         !< node number
   INTEGER( IntKi )            , INTENT(IN   ) :: perturb_sign !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(SD_ContinuousStateType), INTENT(INOUT) :: x            !< perturbed SD states
   REAL( R8Ki )                , INTENT(  OUT) :: dx           !< amount that specific state was perturbed
   if (fieldIndx==1) then
      dx=p%dx(1)
      x%qm(mode)    = x%qm(mode)    + dx * perturb_sign
   else
      dx=p%dx(2)
      x%qmdot(mode) = x%qmdot(mode) + dx * perturb_sign
   end if
END SUBROUTINE SD_Perturb_x
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine sd_init_jacobian is consistant with this routine!
SUBROUTINE SD_Compute_dX(p, x_p, x_m, delta, dX)
   TYPE(SD_ParameterType)      , INTENT(IN   ) :: p                                !< parameters
   TYPE(SD_ContinuousStateType), INTENT(IN   ) :: x_p                              !< SD continuous states at \f$ u + \Delta_p u \f$ or \f$ x + \Delta_p x \f$ (p=plus)
   TYPE(SD_ContinuousStateType), INTENT(IN   ) :: x_m                              !< SD continuous states at \f$ u - \Delta_m u \f$ or \f$ x - \Delta_m x \f$ (m=minus)
   REAL(R8Ki)                  , INTENT(IN   ) :: delta                            !< difference in inputs or states \f$ delta_p = \Delta_p u \f$ or \f$ delta_p = \Delta_p x \f$
   REAL(R8Ki)                  , INTENT(INOUT) :: dX(:)                            !< column of dXdu or dXdx: \f$ \frac{\partial X}{\partial u_i} = \frac{x_p - x_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial X}{\partial x_i} = \frac{x_p - x_m}{2 \, \Delta x}\f$
   INTEGER(IntKi) :: i ! loop over modes
   do i=1,p%Jac_nx
      dX(i) = x_p%qm(i) - x_m%qm(i)
   end do
   do i=1,p%Jac_nx
      dX(p%Jac_nx+i) = x_p%qmdot(i) - x_m%qmdot(i)
   end do
   dX = dX / (2.0_R8Ki*delta)
END SUBROUTINE SD_Compute_dX

END MODULE SubDyn_Output
