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
   USE SubDyn_Output_Params, only: MNfmKe, MNTDss, MNRDe, MNTRAe, IntfSS, IntfTRss, IntfTRAss, IntfTRe, ReactSS, RBTRDss, RBTRVss, RBTRAss
   USE SubDyn_Output_Params, only: ParamIndxAry, ParamUnitsAry, ValidParamAry, SSqm01, SSqmd01, SSqmdd01, OutStrLenM1

   IMPLICIT NONE

   ! The maximum number of output channels which can be output by the code.
   INTEGER(IntKi),PUBLIC, PARAMETER      :: MaxOutPts = 16575

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

CONTAINS


!> This subroutine initializes the output module, checking if the output parameter list (OutList)
! contains valid names, and opening the output file if there are any requested outputs
SUBROUTINE SDOut_Init( Init, y,  p, misc, InitOut, WtrDpth, ErrStat, ErrMsg )
   TYPE(SD_InitType),               INTENT( INOUT ) :: Init                 ! data needed to initialize the output module
   TYPE(SD_OutputType),             INTENT( INOUT ) :: y                    ! SubDyn module's output data
   TYPE(SD_ParameterType), target,  INTENT( INOUT ) :: p                    ! SubDyn module parameters
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
   INTEGER(IntKi)                 :: nElemPerNode, nNodesPerMember ! Number of elements connecting to a node, number of nodes per member
   type(MeshAuxDataType), pointer :: pLst                                                   !< Alias to shorten notation and highlight code similarities
   real(ReKi), allocatable :: T_TIreact(:,:) ! Transpose of TIreact, temporary
   ErrStat = 0
   ErrMsg=""

   p%OutAllDims=6*p%NMembers*2    !size of AllOut Member Joint forces

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
      CALL AllocAry(pLst%Ke, 12, 12, pLst%NoutCnt, 2, 'MOutLst(I)%Ke'     , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%Fg,     12, pLst%NoutCnt, 2, 'MOutLst(I)%Fg'     , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%extrap,     pLst%NoutCnt   , 'MOutLst(I)%extrap' , ErrStat2, ErrMsg2); if(Failed()) return

      ! NOTE: len(MemberNodes) >2 if nDiv>1
      iMember = FINDLOCI(Init%Members(:,1), pLst%MemberID) ! Reindexing from MemberID to 1:nMembers
      nNodesPerMember = count(Init%MemberNodes(iMember,:)>0_IntKi)
      pLst%NodeIDs(1:pLst%NoutCnt)=Init%MemberNodes(iMember, pLst%NodeCnt)  ! We are storing the actual node numbers corresponding to what the user ordinal number is requesting
      pLst%ElmIDs=0  !Initialize to 0
      pLst%ElmNds=0  !Initialize to 0
      pLst%extrap=.false.

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
         if ( (K2==2_IntKi).or.(nNodesPerMember==2_IntKi) ) cycle ! No need to proceed further if we have an interior node or only one element
         ! Save neighboring element info for force extrapolation to an end node if more than 1 element per member
         if (iNode == Init%MemberNodes(iMember,1)) then                ! First node of the member
            iNode = Init%MemberNodes(iMember,2)                        ! Index of the second node
         else if (iNode == Init%MemberNodes(iMember,nNodesPerMember)) then ! Last node of the member
            iNode = Init%MemberNodes(iMember,nNodesPerMember-1_IntKi)      ! Index of the second to last node
         end if
         do iiElem = 1, 2 ! Should have exactly two elements connecting to an interior node
            iElem = Init%NodesConnE(iNode, iiElem+1) ! iiElem-th element Number; no need to call ThisElementIsAlongMember since interior node but kept for safety for now
            if ( ThisElementIsAlongMember(iElem, iNode, iMember) .and. iElem /= pLst%ElmIDs(J,1_IntKi)) then
                call ConfigOutputNode_MKF_ID(pLst, iElem, iiNode=J, iStore=2, NodeID2=iNode)
                pLst%extrap(J) = .true.
            end if
         end do  ! iiElem, nElemPerNode
      ENDDO !J, Noutcnt
   ENDDO  !I, NMOutputs

   !_________________________________ OUTPUT FOR ALL MEMBERS __________________________________
   IF (p%OutAll) THEN  !I need to store all member end forces and moments

      ! MOutLst2: nodal output info by members, for all members, First and Last Node
      ALLOCATE ( p%MOutLst2(p%NMembers), STAT = ErrStat2 ); ErrMsg2 = 'Error allocating p%MOutLst2 array in SDOut_Init'; if(Failed()) return

      DO iMember=1,p%NMembers
         pLst => p%MOutLst2(iMember) ! Alias to shorten notations
         CALL AllocAry(pLst%NodeIDs,    2   , 'MOutLst(I)%NodeIDs', ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%ElmIDs,     2, 2, 'MOutLst(I)%ElmIDs' , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%ElmNds,     2, 2, 'MOutLst(I)%ElmNds' , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Ke, 12, 12, 2, 2, 'MOutLst(I)%Ke'     , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Fg,     12, 2, 2, 'MOutLst(I)%Fg'     , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%extrap,     2   , 'MOutLst(I)%extrap' , ErrStat2, ErrMsg2); if(Failed()) return
         pLst%MemberID = Init%Members(iMember,1)
         nNodesPerMember = count(Init%MemberNodes(iMember,:)>0_IntKi)
         pLst%NodeIDs(1) = Init%MemberNodes(iMember,1)             ! First node of the member
         pLst%NodeIDs(2) = Init%MemberNodes(iMember,nNodesPerMember) ! Last node of the member
         pLst%ElmIDs=0  !Initialize to 0
         pLst%ElmNds=0  !Initialize to 0
         pLst%extrap=.false.
         DO J=1,2 ! loop on requested nodes for that member
            iNode        = pLst%NodeIDs(J)           ! Index of requested node in node list
            nElemPerNode = Init%NodesConnE(iNode, 1) ! Number of elements connecting to the j-th node
            ! Finding the element that belongs to the member and connect to the node
            DO iiElem = 1, nElemPerNode
               iElem = Init%NodesConnE(iNode, iiElem+1) ! iiElem-th Element Number
               IF (ThisElementIsAlongMember(iElem, iNode, iMember)) THEN
                  call ConfigOutputNode_MKF_ID(pLst, iElem, iiNode=J, iStore=1, NodeID2=iNode)
                  exit ! End nodes can only have one connected element on the member
               END IF
            ENDDO  ! iiElem, nElemPerNode
            if ( nNodesPerMember==2_IntKi ) cycle ! No need to proceed further if we only have one element
            ! Save neighboring element info for force extrapolation to an end node if more than 1 element per member
            if (iNode == Init%MemberNodes(iMember,1)) then                ! First node of the member
               iNode = Init%MemberNodes(iMember,2)                        ! Index of the second node
            else if (iNode == Init%MemberNodes(iMember,nNodesPerMember)) then ! Last node of the member
               iNode = Init%MemberNodes(iMember,nNodesPerMember-1_IntKi)      ! Index of the second to last node
            end if
            do iiElem = 1,2 ! Should have exactly two elements connecting to an interior node
               iElem = Init%NodesConnE(iNode, iiElem+1) ! iiElem-th element Number; no need to call ThisElementIsAlongMember since interior node but kept for safety for now
               if ( ThisElementIsAlongMember(iElem, iNode, iMember) .and. iElem /= pLst%ElmIDs(J,1_IntKi)) then
                   call ConfigOutputNode_MKF_ID(pLst, iElem, iiNode=J, iStore=2, NodeID2=iNode)
                   pLst%extrap(J) = .true.
               end if
            end do  ! iiElem, nElemPerNode
         ENDDO !J, Noutcnt
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
      call RigidTrnsf(Init, p, (/0.0_Reki, 0.0_ReKi, -WtrDpth /), p%IDC__, p%nDOFC__, 1_IntKi, T_TIreact, ErrStat2, ErrMsg2); if(Failed()) return
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
   !!   - Mass, stiffness matrices and constant element force vector Fg:
   !!       - Beam elements: gravity fixed-end bending moments only (force DOFs zeroed because self-weight forces
   !!         are corrected elsewhere via force extrapolation that addresses hydro loads)
   !!       - Cable elements: initial pretension nodal force vector (to recover total tension T_pretension + k*delta)
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
      ! CALL ElemM(p%ElemProps(iElem),         pLst%Me(:,:,iiNode,iStore))
      CALL ElemK(p%ElemProps(iElem),         pLst%Ke(:,:,iiNode,iStore))
      CALL ElemF(p%ElemProps(iElem), Init%g, pLst%Fg(:,iiNode,iStore), FCe)
      ! Fg is set differently depending on element type:
      !   Beam: keep only the gravity fixed-end bending moments; zero the force DOFs (1:3, 7:9) because
      !         self-weight forces are corrected elsewhere via force extrapolation that addresses hydro loads.
      !   Cable: replace Fg with FCe (pretension), so CALC_NODE_FORCES recovers total tension: T_pretension + k*delta.
      ! Note: for floating systems, pLst%Fg bending moments for beams are recomputed in ElementForce
      !       using the current element orientation. Translational force DOFs are zeroed and not used.
      if (p%ElemProps(iElem)%eType == idMemberBeamCirc .or. &
          p%ElemProps(iElem)%eType == idMemberBeamRect .or. &
          p%ElemProps(iElem)%eType == idMemberBeamArb) then
         pLst%Fg(1:3,iiNode,iStore) = 0.0_FeKi
         pLst%Fg(7:9,iiNode,iStore) = 0.0_FeKi
      else if (p%ElemProps(iElem)%eType == idMemberCable) then
         pLst%Fg(:,iiNode,iStore) = FCe(1:12)
      endif
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
   real(ReKi), dimension (6)      :: FK_elm, FK_elm2   ! output elastic forces and moments
   real(ReKi), dimension (6)      :: Fext      ! external forces and moments
  real(FEKi), dimension (3,3)    :: DIRCOS    ! direction cosine matrix (global to local) (3x3)
   real(ReKi), allocatable        :: ReactNs(:)    ! 6*Nreact reactions
   integer(IntKi)                 :: sgn ! +1/-1 for node force calculations
   type(MeshAuxDataType), pointer :: pLst       !< Info for a given member-output (Alias to shorten notation)
   integer(IntKi), pointer        :: DOFList(:) !< List of DOF indices for a given Nodes (Alias to shorten notation)
   real(R8Ki), dimension(3,3)     :: Rg2b  ! Rotation matrix global 2 body (Guyan) coordinates
   real(R8Ki), dimension(6,6)     :: RRg2b ! Rotation matrix global 2 body (Guyan) coordinates, acts on a 6-vector
   real(R8Ki), dimension(6,6)     :: RRb2g ! Rotation matrix global 2 body (Guyan) coordinates, acts on a 6-vector
   integer(IntKi)                 :: iTP, nTP
   INTEGER(IntKi)                 :: ErrStat2      ! Error status of the operation
   CHARACTER(ErrMsgLen)           :: ErrMsg2       ! Error message if ErrStat /= ErrID_None

   ErrStat = ErrID_None
   ErrMsg  = ""

   if ( p%Floating ) then
      ! For floating, m%U_full_dotdot is currently in the earth-fixed frame.
      ! Need to transform back to the Guyan/rigid-body frame.
      if ( p%TP1IsRBRefPt ) then
         Rg2b = EulerConstructZYX(x%qR(4:6))
      else
         Rg2b = u%TPMesh(1)%Orientation(:,:,1)
      endif
   else
      call Eye(Rg2b, ErrStat2, ErrMsg2)
   end if
   RRg2b = 0.0_R8Ki
   RRg2b(1:3,1:3) = Rg2b
   RRg2b(4:6,4:6) = Rg2b
   RRb2g = transpose(RRg2b)

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
            ! --- Forces (potentially averaged across or extrapolated from 2 elements)
            call ElementForce(pLst, iiNode, 1, FK_elm, sgn, DIRCOS, .false.)
            FK_elm=sgn*FK_elm
            IF (pLst%ElmIDs(iiNode,2) .NE. 0) THEN  ! Second element exist
               if (pLst%extrap(iiNode)) then ! Linearly extrapolate forces to end nodes
                  ! NOTE: forces are computed in the coordinate system of the first element for extrapolating
                  call ElementForce(pLst, iiNode, 2, FK_elm2, sgn, DIRCOS, .true.) ! True= we use DIRCOS from element above
                  FK_elm(1:3) = 1.5_ReKi * FK_elm(1:3) - 0.5_ReKi * sgn*FK_elm2(1:3) ! Now extrapolate
               else ! Average/interpolate forces and moments at internal nodes
                  ! NOTE: forces are computed in the coordinate system of the first element for averaging
                  call ElementForce(pLst, iiNode, 2, FK_elm2, sgn, DIRCOS, .true.) ! True= we use DIRCOS from element above
                  FK_elm = 0.5_ReKi * ( FK_elm + sgn*FK_elm2 ) ! Now Average
               end if
            ENDIF
            ! Elastic component of reaction forces and moments at MαNβ along local member coordinate system
            !    "MαNβFKxe, MαNβFKye, MαNβFKze, MαNβMKxe, MαNβMKye, MαNβMKze"
            AllOuts(MNfmKe  (:,iiNode,iMemberOutput)) = FK_elm  !elastic forces and moments (6) Local Ref

            ! --- Displacements and acceleration
            ! Revolute joints can have more than 6 DOF. Need element-based DOF look-up instead.
            ! DOFList => p%NodesDOF(pLst%NodeIDs(iiNode))%List
            select case (pLst%ElmNds(iiNode,1_IntKi))
            case (1_IntKi) ! First node of the element
               DOFList => p%ElemsDOF(1:6 ,pLst%ElmIDs(iiNode,1_IntKi))
            case (2_IntKi) ! Second node of the element
               DOFList => p%ElemsDOF(7:12,pLst%ElmIDs(iiNode,1_IntKi))
            end select
            ! Displacement- Translational -no need for averaging since it is a node translation - In global reference SS
            !     "MαNβTDxss, MαNβTDyss, MαNβTDzss"
            AllOuts(MNTDss (:,iiNode,iMemberOutput))       = m%U_full(DOFList(1:3))
            ! Displacement- Rotational - need direction cosine matrix to tranform rotations  - In Local reference Element Ref Sys <- Need to rethink this for large platform rotation
            !     "MαNβRDxe, MαNβRDye, MαNβRDze"
            AllOuts(MNRDe (:,iiNode,iMemberOutput))        = matmul(DIRCOS,m%U_full_elast(DOFList(4:6))) ! Element elastic rotation only in Guyan frame for floating. Full motion for fixed-bottom.
            ! Accelerations- I need to get the direction cosine matrix to tranform displacement and rotations
            !     "MαNβTAxe, MαNβTAye, MαNβTAze"
            !     "MαNβRAxe, MαNβRAye, MαNβRAze"
            AllOuts(MNTRAe (1:3,iiNode,iMemberOutput))     = matmul(DIRCOS,matmul(Rg2b,m%U_full_dotdot(DOFList(1:3)))) ! translational accel local ref
            AllOuts(MNTRAe (4:6,iiNode,iMemberOutput))     = matmul(DIRCOS,matmul(Rg2b,m%U_full_dotdot(DOFList(4:6)))) ! rotational    accel local ref
        ENDDO  ! iiNode, Loop on requested nodes for that member
     ENDDO ! iMemberOutput, Loop on member outputs
   END IF

   ! --------------------------------------------------------------------------------
   ! --- All nodal loads from stiffness and mass matrix
   ! --------------------------------------------------------------------------------
   ! "MaaaJbFKxe, MaaaJbMKxe for member aaa and node b."
   IF (p%OutAll) THEN
      DO iMemberOutput=1,p%NMembers    !Cycle on all members
         pLst=>p%MOutLst2(iMemberOutput)
         DO iiNode=1,2 !Iterate on requested nodes for that member (first and last)
            ! --- Forces (potentially extrapolated from 2 elements)
            call ElementForce(pLst, iiNode, 1, FK_elm, sgn, DIRCOS, .false.)
            FK_elm=sgn*FK_elm
            if ( (pLst%ElmIDs(iiNode,2)/=0) .and. pLst%extrap(iiNode) ) then  ! Extrapolate forces to end nodes
               ! NOTE: forces are computed in the coordinate system of the first element for extrapolating
               call ElementForce(pLst, iiNode, 2, FK_elm2, sgn, DIRCOS, .true.) ! True= we use DIRCOS from element above
               FK_elm(1:3) = 1.5_ReKi * FK_elm(1:3) - 0.5_ReKi * sgn*FK_elm2(1:3) ! Now extrapolate
            end if
            ! Store in All Outs
            L  = MaxOutPts+(iMemberOutput-1)*12+(iiNode-1)*6+1
            L2 = L+5
            AllOuts( L:L2 ) = FK_elm
         ENDDO !iiNode, nodes 1 and 2
      ENDDO ! iMemberOutput, Loop on members
   ENDIF

   ! --------------------------------------------------------------------------------
   ! --- Interface kinematics and loads (TP/platform reference point)
   ! --------------------------------------------------------------------------------
   if (p%TP1IsRBRefPt) then
      nTP = p%nTP-1
   else
      nTP = p%nTP
   end if

   ! Total interface reaction forces and moments in SS coordinate system
   !    "IntfFXss, IntfFYss, IntfFZss, IntfMXss, IntfMYss, IntfMZss,"
   do iTP = 1,nTP
      AllOuts(IntfSS(1:6,iTP)) = - (/y%Y1Mesh(iTP)%Force(:,1), y%Y1Mesh(iTP)%Moment(:,1)/) !-y%Y1  !Note this is the force that the TP applies to the Jacket, opposite to what the GLue Code needs thus "-" sign
   end do

   ! Interface translations and rotations in SS coordinate system
   !    "IntfTDXss, IntfTDYss, IntfTDZss, IntfRDXss, IntfRDYss IntfRDZss"
   do iTP = 1,nTP
      AllOuts(IntfTRss(1:3,iTP)) = u%TPMesh(iTP)%TranslationDisp(:,1)
      AllOuts(IntfTRss(4:6,iTP)) = EulerExtractZYX(u%TPMesh(iTP)%Orientation(:,:,1))
   end do

   ! Interface Translational and rotational accelerations in SS coordinate system
   !    "IntfTAXss, IntfTAYss, IntfTAZss, IntfRAXss, IntfRAYss IntfRAZss"
   do iTP = 1,nTP
      AllOuts(IntfTRAss(1:3,iTP)) = u%TPMesh(iTP)%TranslationAcc(:,1)
      AllOuts(IntfTRAss(4:6,iTP)) = u%TPMesh(iTP)%RotationAcc(:,1)
   end do

   ! Interface elastic translational and rotational deflection in rigid-body coordinate system relative to rigid-body configuration
   !    "IntfTDXe, IntfTDYe, IntfTDZe, IntfRDXe, IntfRDYe IntfRDZe"
   if (p%floating) then
      if (p%TP1IsRBRefPt) then
         do iTP = 1,nTP
            AllOuts(IntfTRe(1:6,iTP)) = m%u_TP( iTP*6+1:iTP*6+6 )
         end do
      else
         AllOuts(IntfTRe(1:6,1  )) = 0.0
      endif
   else
      do iTP = 1,nTP
         AllOuts(IntfTRe(1:6,iTP)) = m%u_TP( (iTP-1)*6+1:(iTP-1)*6+6 )
      end do
   end if

   ! --------------------------------------------------------------------------------
   ! --- Interface kinematics and loads (TP/platform reference point)
   ! --------------------------------------------------------------------------------
   if (p%floating) then
      AllOuts(RBTRDss) = m%u_TP(1:6)
      AllOuts(RBTRVss) = matmul( RRb2g, m%udot_TP(1:6)    )
      AllOuts(RBTRAss) = matmul( RRb2g, m%udotdot_TP(1:6) )
   else
      AllOuts(RBTRDss) = 0.0
      AllOuts(RBTRVss) = 0.0
      AllOuts(RBTRAss) = 0.0
   end if

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
      IF (p%SlDNonLinear) THEN
         ! When SoilDyn nonlinear loads are active (e.g., SoilDyn CalcOption = 3), SubDyn reaction loads are incomplete (only the linear part is included)
         ! The total reaction at each base reaction joint is available in the SoilDyn output sensors (e.g., "Sld1Fxg Sld1Fyg Sld1Fzg Sld1Mxg Sld1Myg Sld1Mzg")
         AllOuts( ReactSS(1:nDOFL_TP) ) = NaN
      ELSE
         ALLOCATE ( ReactNs(6*p%nNodes_C), STAT = ErrStat )
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Error allocating space for ReactNs array.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         ReactNs = 0.0_ReKi !Initialize
         DO I=1,p%nNodes_C   !Do for each constrained node, they are ordered as given in the input file and so as in the order of y2mesh
            FK_elm2=0._ReKi !Initialize for cumulative force
            pLst => p%MOutLst3(I)
            !Find the joint forces
            DO J=1,SIZE(pLst%ElmIDs(1,:))  !for all the elements connected (normally 1)
               iiNode = 1
               call ElementForce(pLst, iiNode, J, FK_elm, sgn, DIRCOS, .false.)
               !transform back to global, need to do 3 at a time since cosine matrix is 3x3
               DO L=1,2
                  FK_elm2((L-1)*3+1:L*3) = FK_elm2((L-1)*3+1:L*3) + matmul(transpose(DIRCOS),FK_elm((L-1)*3+1:L*3))  !signs may be wrong, we will fix that later;
                  ! I believe this is all fixed in terms of signs now ,RRD 5/20/13
               ENDDO
            ENDDO
            ! NEED TO ADD HYDRODYNAMIC FORCES AT THE RESTRAINT NODES
            iSDNode   = p%Nodes_C(I,1)
            iMeshNode = iSDNode ! input and Y2 mesh nodes are the same as subdyn
            Fext =  (/ u%LMesh%Force(:,iMeshNode), u%LMesh%Moment(:,iMeshNode) /) + p%FG(p%NodesDOF(iMeshNode)%List(1:6))
            Fext(1:3) = Fext(1:3) + p%FC(p%NodesDOF(iMeshNode)%List(1:3))
            ReactNs((I-1)*6+1:6*I) = FK_elm2 - Fext  !Accumulate reactions from all nodes in GLOBAL COORDINATES
         ENDDO
         ! Store into AllOuts
         AllOuts( ReactSS(1:nDOFL_TP) ) = matmul(p%TIreact,ReactNs)
      ENDIF
   ENDIF
   if (allocated(ReactNs)) deallocate(ReactNs)
contains

   subroutine ElementForce(pLst, iiNode, JJ, FK_elm, sgn, DIRCOS, bUseInputDirCos)
      type(MeshAuxDataType),       intent(in)    :: pLst   !< Info for one member output
      integer(IntKi),              intent(in)    :: iiNode !< Index over the nodes of a given member (>2 if nDIV>1)
      integer(IntKi),              intent(in)    :: JJ     !< TODO: interpretation: index over other member connected to the current member (for averaging)
      real(FEKi), dimension (3,3), intent(inout) :: DIRCOS  !direction cosice matrix (global to local) (3x3)
      real(ReKi), dimension (6),   intent(out)   :: FK_elm  !output elastic forces and moments
      integer(IntKi),              intent(out)   :: sgn !+1/-1 for node force calculations
      logical,                     intent(in)    :: bUseInputDirCos !< If True, use DIRCOS from input, otherwise, use element DirCos
      ! Local
      integer(IntKi)                          :: iElem !< Element index/number
      integer(IntKi)                          :: FirstOrSecond !< 1 or 2  if first node or second node
      integer(IntKi), dimension(2)            :: ElemNodes  ! Node IDs for element under consideration (may not be consecutive numbers)
      real(ReKi)    , dimension(12)           :: X_e        ! Deflection of an element
      real(FEKi)    , dimension(12)           :: Fg_e ! Gravity force (beam elements) or initial pretension (cable elements), re-oriented for floating systems
      real(FEKi)    , dimension(3,3)          :: CurDirCos ! Current element direction cosine matrix in the floating body frame
      integer(IntKi), dimension(2), parameter :: NodeNumber_To_Sign = (/-1, +1/)

      iElem         = pLst%ElmIDs(iiNode,JJ)             ! element number
      FirstOrSecond = pLst%ElmNds(iiNode,JJ)             ! first or second node of the element to be considered
      sgn           = NodeNumber_To_Sign(FirstOrSecond) ! Assign sign depending if it's the 1st or second node
      ElemNodes     = p%Elems(iElem,2:3)                ! first and second node ID associated with element iElem
      ! Note that a node can have more than 6DOF if it is a revolute joint; must use element-based DOF lookup instead
      ! X_e(1:6)      = m%U_full_elast (p%NodesDOF(ElemNodes(1))%List(1:6))   ! For floating, m%U_full_elast is the CB+SIM elastic deformation only in the Guyan (rigid-body) frame
      ! X_e(7:12)     = m%U_full_elast (p%NodesDOF(ElemNodes(2))%List(1:6))   ! No additional transformation required
      X_e(1:6)      = m%U_full_elast (p%ElemsDOF(1:6 ,iElem))   ! For floating, m%U_full_elast is the CB+SIM elastic deformation only in the Guyan (rigid-body) frame
      X_e(7:12)     = m%U_full_elast (p%ElemsDOF(7:12,iElem))   ! No additional transformation required

      ! Load Fg: gravity force (beam elements self-weight) or initial pretension (cable elements), computed at initialization.
      ! For floating systems:
      !   - Beam self-weight force components are rotated from the initial to the current body/Guyan frame.
      !   - Beam self-weight bending moment components are recomputed using the current element orientation.
      !   - Cable pretension force is already expressed along the local z-axis. No need to rotate.
      Fg_e = real(pLst%Fg(:,iiNode,JJ), R8Ki)
      if (p%Floating) then
         if (p%ElemProps(iElem)%eType == idMemberBeamCirc .or. &
             p%ElemProps(iElem)%eType == idMemberBeamRect .or. &
             p%ElemProps(iElem)%eType == idMemberBeamArb) then
            ! Beam elements self-weight forces were zeroed out. Otherwise, they would require this rotation as well
            ! Fg_e(1:3) = matmul(Rg2b, Fg_e(1:3))
            ! Fg_e(7:9) = matmul(Rg2b, Fg_e(7:9))

            ! Recompute beam self-weight bending moments using the current element orientation
            ! CurDirCos = Rb2g * DirCos0 (DirCos0 is the element direction cosine matrix at initialization)
            CurDirCos = matmul(transpose(Rg2b), p%ElemProps(iElem)%DirCos)
            Fg_e(4)  = -p%ElemProps(iElem)%Length**2 * p%ElemProps(iElem)%Rho * p%ElemProps(iElem)%Area * p%g / 12.0_FEKi * CurDirCos(2,3)
            Fg_e(5)  =  p%ElemProps(iElem)%Length**2 * p%ElemProps(iElem)%Rho * p%ElemProps(iElem)%Area * p%g / 12.0_FEKi * CurDirCos(1,3)
            Fg_e(6)  = 0.0_FEKi   ! no torsional self-weight moment
            Fg_e(10) = -Fg_e(4)
            Fg_e(11) = -Fg_e(5)
            Fg_e(12) = 0.0_FEKi
         endif
      endif
      if (.not. bUseInputDirCos) then
         DIRCOS=transpose(p%ElemProps(iElem)%DirCos)! global to local
      endif
      CALL CALC_NODE_FORCES( DIRCOS, pLst%Ke(:,:,iiNode,JJ), X_e, Fg_e, FirstOrSecond, FK_elm)
   end subroutine ElementForce

   !====================================================================================================
   !> Calculates elastic forces for a given element, using K of the element
   !  Fg is the beam element gravity and the initial cable pretension load vector.
   !  FirstOrSecond selects whether the node of interest is the first (1) or second (2) node of the element.
   !----------------------------------------------------------------------------------------------------
   SUBROUTINE CALC_NODE_FORCES(DIRCOS, Ke, Y2, Fg, FirstOrSecond, FK_nod)
      Real(FEKi), DIMENSION (3,3),   INTENT(IN)  :: DIRCOS          ! direction cosice matrix (global to local) (3x3)
      Real(FEKi), DIMENSION (12,12), INTENT(IN)  :: Ke              ! element K matrices (12x12) in GLOBAL REFERENCE (DIRCOS^T K DIRCOS)
      Real(ReKi), DIMENSION (12),    INTENT(IN)  :: Y2              ! element elastic deflection
      Real(FEKi), DIMENSION (12),    INTENT(IN)  :: Fg              ! element load vector from gravity and initial cable pretension (orientation dependent for floating, constant for fixed-bottom)
      Integer(IntKi),                INTENT(IN)  :: FirstOrSecond   ! 1 or 2 depending on node of interest
      REAL(ReKi), DIMENSION (6),     INTENT(OUT) :: FK_nod          ! output elastic forces and moments
      !Locals
      INTEGER(IntKi)                             :: L               !counter
      REAL(DbKi), DIMENSION(12)                  :: FF_glb, FF_elm  ! temporary storage

      FF_glb = matmul(Ke,Y2) - Fg ! GLOBAL REFERENCE (Guyan/rigid-body frame if floating)
      DO L=1,4 ! Transforming coordinates 3 at a time
         FF_elm((L-1)*3+1:L*3) =  matmul(DIRCOS, FF_glb( (L-1)*3+1:L*3 ) ) 
      ENDDO
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

   ! Initialize to -1 to indicate that the output file unit is not valid
   p%UnJckF = -1

   ! No outputs requested, so just return
   if ((.not. allocated(p%OutParam)) .or. (p%NumOuts == 0)) then
      call WrScr('SubDyn: no outputs were requested, so separate output file will not be generated.')
      return
   end if

   ! Open the file for output
   OutFileName = TRIM(OutRootName)//'.out'
   call GetNewUnit( p%UnJckF )

   call OpenFOutFile ( p%UnJckF, OutFileName, ErrStat, ErrMsg )
   if (ErrStat >= AbortErrLev) then
      ErrMsg = ' Error opening SubDyn-level output file: '//TRIM(ErrMsg)
      return
   end if

   ! Write the output file header
   write(p%UnJckF,'(/,A/)', IOSTAT=ErrStat2)  'These predictions were generated by '//TRIM(GETNVD(ProgVer))//&
                  ' on '//CurDate()//' at '//CurTime()//'.'

   write(p%UnJckF, '(//)') ! add 3 lines to make file format consistant with FAST v8 (headers on line 7; units on line 8) [this allows easier post-processing]

   ! Write the names of the output parameters:
   Frmt = '(A8,'//TRIM(Int2LStr(p%NumOuts+p%OutAllInt*p%OutAllDims))//'(:,A,'//TRIM( p%OutSFmt )//'))'
   write(p%UnJckF,Frmt, IOSTAT=ErrStat2)  TRIM( 'Time' ), ( p%Delim, TRIM( InitOut%WriteOutputHdr(I) ), I=1,p%NumOuts+p%OutAllInt*p%OutAllDims )

   ! Write the units of the output parameters:
   write(p%UnJckF,Frmt, IOSTAT=ErrStat2)  TRIM( 's'), ( p%Delim, TRIM( InitOut%WriteOutputUnt(I) ), I=1,p%NumOuts+p%OutAllInt*p%OutAllDims )
END SUBROUTINE SDOut_OpenOutput

!====================================================================================================


!====================================================================================================
! SDOut_CloseOutput closes the output file, if open, after running the SubDyn output module.
!----------------------------------------------------------------------------------------------------
SUBROUTINE SDOut_CloseOutput ( p, ErrStat, ErrMsg )
   type(SD_ParameterType),  INTENT( INOUT ) :: p        ! data for this instance of the floating platform module
   integer(IntKi),          INTENT(   OUT ) :: ErrStat  ! a non-zero value indicates an error occurred
   character(*),            INTENT(   OUT ) :: ErrMsg   ! Error message if ErrStat /= ErrID_None
   integer(IntKi)                           :: Stat

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! If file is not open, return
   if (p%UnJckF == -1) return

   ! Close our output file
   close(p%UnJckF, iostat=Stat)
   if (Stat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Problem closing SubDyn output file.'
   end if

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

   ! If output file is not open, return
   if (p%UnJckF == -1) return

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
   CHARACTER(ChanLen), DIMENSION(12)         :: ToTUnits,ToTNames,ToTNames0
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

   DO I=1,99
          !I know el # and whether it is 1st node or second node
      if (I <= p%NMOutputs) then
         INDX=p%MOutLst(I)%NOutCnt+1
      else
         INDX = 1
      end if

      DO J=INDX,9 !Iterate on requested nodes for that member
         !Forces and moments
         InvalidOutput(MNfmKe  (:,J,I)) = .true.  !elastic forces and moments (6) Local Ref
         !Displacement
         InvalidOutput(MNTDss  (:,J,I)) = .true.  !Translational
         InvalidOutput(MNRDe   (:,J,I)) = .true.  !Rotational
         !Accelerations
         InvalidOutput(MNTRAe  (:,J,I)) = .true.  !translational accel local ref
      END DO
   END DO

   IF (p%TP1IsRBRefPt) THEN
      DO I=p%nTP,9
         InvalidOutput(IntfSS(:,I))    = .true.
         InvalidOutput(IntfTRe(:,I))   = .true.
         InvalidOutput(IntfTRss(:,I))  = .true.
         InvalidOutput(IntfTRAss(:,I)) = .true.
      END DO
   ELSE
      DO I=p%nTP+1,9
         InvalidOutput(IntfSS(:,I))    = .true.
         InvalidOutput(IntfTRe(:,I))   = .true.
         InvalidOutput(IntfTRss(:,I))  = .true.
         InvalidOutput(IntfTRAss(:,I)) = .true.
      END DO
   END IF

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

      CALL Conv2UC( OutListTmp )    ! Convert OutListTmp to upper case

      ! Interface output backward compatibility
      k = INDEX( OutListTmp, 'INTF' )
      IF ( k>0 .and. INDEX( '0123456789', OutListTmp(k+4:k+4) ) <= 0 ) THEN
         OutListTmp = OutListTmp(1:k+3)//'1'//OutListTmp(k+4:)
      END IF

      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !   channel name with a '-', '_', 'm', or 'M' character indicating "minus".

      CheckOutListAgain = .FALSE.

      IF      ( INDEX( '-_', OutListTmp(1:1) ) > 0 ) THEN
         p%OutParam(I)%SignM = -1     ! ex, '-TipDxc1' causes the sign of TipDxc1 to be switched.
         OutListTmp                   = OutListTmp(2:)
      ELSE IF ( INDEX( 'mM', OutListTmp(1:1) ) > 0 ) THEN ! We'll assume this is a variable name for now, (if not, we will check later if OutListTmp(2:) is also a variable name)
         CheckOutListAgain  = .TRUE.
         p%OutParam(I)%SignM = 1
      ELSE
         p%OutParam(I)%SignM = 1
      END IF

      if ( INDEX( 'mM', OutListTmp(1:1) ) > 0 .and. INDEX( '0123456789', OutListTmp(2:2) ) > 0 .and. INDEX( 'nN', OutListTmp(3:3) ) > 0 ) then ! an old-style output without the leading zero on the member number
         OutListTmp = OutListTmp(1:1)//'0'//OutListTmp(2:)
         CheckOutListAgain  = .FALSE.
      end if

      Indx =  IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )

      IF ( CheckOutListAgain .AND. Indx < 1 ) THEN    ! Let's assume that "M" really meant "minus" and then test again
         p%OutParam(I)%SignM = -1            ! ex, 'MTipDxc1' causes the sign of TipDxc1 to be switched.
         OutListTmp                   = OutListTmp(2:)

         Indx = IndexCharAry( OutListTmp(1:10), ValidParamAry )
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
       ToTNames0=RESHAPE(SPREAD( (/"FKxe" ,"FKye" ,"FKze" ,"MKxe" ,"MKye" ,"MKze" /), 2, 2), (/12/) )
       ToTUnits =RESHAPE(SPREAD( (/"(N)  ","(N)  ","(N)  ","(N*m)","(N*m)","(N*m)"/), 2, 2), (/12/) )
       DO I=1,p%NMembers
           DO K=1,2
            DO J=1,6
             TotNames(J+(K-1)*6)=TRIM("M"//Int2Lstr(I))//TRIM("J"//Int2Lstr(K))//TRIM(ToTNames0(J))
            ENDDO  
           ENDDO
           p%OutParam(p%NumOuts+(I-1)*6*2+1:p%NumOuts+I*6*2)%Name  = ToTNames
           p%OutParam(p%NumOuts+(I-1)*6*2+1:p%NumOuts+I*6*2)%Units = ToTUnits
       ENDDO
       p%OutParam(p%NumOuts+1:p%NumOuts+p%OutAllDims)%SignM = 1
       p%OutParam(p%NumOuts+1:p%NumOuts+p%OutAllDims)%Indx= MaxOutPts+(/(J, J=1, p%OutAllDims)/)
   ENDIF

END SUBROUTINE SDOut_ChkOutLst
!====================================================================================================

END MODULE SubDyn_Output
