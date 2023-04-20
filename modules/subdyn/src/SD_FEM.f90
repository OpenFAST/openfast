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
!**********************************************************************************************************************************
MODULE SD_FEM
  USE NWTC_Library
  USE SubDyn_Types
  USE FEM
  IMPLICIT NONE

 
  INTEGER(IntKi),   PARAMETER  :: MaxMemJnt       = 20                    ! Maximum number of members at one joint
  INTEGER(IntKi),   PARAMETER  :: MaxOutChs       = 2000                  ! Max number of Output Channels to be read in
  INTEGER(IntKi),   PARAMETER  :: nDOFL_TP        = 6  !TODO rename me    ! 6 degrees of freedom (length of u subarray [UTP])
   
  ! values of these parameters are ordered by their place in SubDyn input file:
  INTEGER(IntKi),   PARAMETER  :: JointsCol       = 9                    ! Number of columns in Joints (JointID, JointXss, JointYss, JointZss, JointType, JointDirX JointDirY JointDirZ JointStiff)
  INTEGER(IntKi),   PARAMETER  :: InterfCol       = 7                     ! Number of columns in interf matrix (JointID,ItfTDxss,ItfTDYss,ItfTDZss,ItfRDXss,ItfRDYss,ItfRDZss)
  INTEGER(IntKi),   PARAMETER  :: ReactCol        = 7                     ! Number of columns in reaction matrix (JointID,ItfTDxss,ItfTDYss,ItfTDZss,ItfRDXss,ItfRDYss,ItfRDZss)
  INTEGER(IntKi),   PARAMETER  :: MaxNodesPerElem = 2                     ! Maximum number of nodes per element (currently 2)
  INTEGER(IntKi),   PARAMETER  :: MembersCol      = MaxNodesPerElem + 3+1+1 ! Number of columns in Members (MemberID,MJointID1,MJointID2,MPropSetID1,MPropSetID2,COSMID) 
  INTEGER(IntKi),   PARAMETER  :: PropSetsBCol    = 6                     ! Number of columns in PropSets  (PropSetID,YoungE,ShearG,MatDens,XsecD,XsecT)  !bjj: this really doesn't need to store k, does it? or is this supposed to be an ID, in which case we shouldn't be storing k (except new property sets), we should be storing IDs
  INTEGER(IntKi),   PARAMETER  :: PropSetsXCol    = 10                    ! Number of columns in XPropSets (PropSetID,YoungE,ShearG,MatDens,XsecA,XsecAsx,XsecAsy,XsecJxx,XsecJyy,XsecJ0)
  INTEGER(IntKi),   PARAMETER  :: PropSetsCCol    = 5                     ! Number of columns in CablePropSet (PropSetID, EA, MatDens, T0)
  INTEGER(IntKi),   PARAMETER  :: PropSetsRCol    = 2                     ! Number of columns in RigidPropSet (PropSetID, MatDens)
  INTEGER(IntKi),   PARAMETER  :: COSMsCol        = 10                    ! Number of columns in (cosine matrices) COSMs (COSMID,COSM11,COSM12,COSM13,COSM21,COSM22,COSM23,COSM31,COSM32,COSM33)
  INTEGER(IntKi),   PARAMETER  :: CMassCol        = 11                    ! Number of columns in Concentrated Mass (CMJointID,JMass,JMXX,JMYY,JMZZ, Optional:JMXY,JMXZ,JMYZ,CGX,CGY,CGZ)
  ! Indices in Members table
  INTEGER(IntKi),   PARAMETER  :: iMType= 6 ! Index in Members table where the type is stored
  INTEGER(IntKi),   PARAMETER  :: iMDirCosID = 7 ! Index in Members table where the type is stored
  INTEGER(IntKi),   PARAMETER  :: iMProp= 4 ! Index in Members table where the PropSet1 and 2 are stored

  ! Indices in Joints table
  INTEGER(IntKi),   PARAMETER  :: iJointType= 5  ! Index in Joints where the joint type is stored
  INTEGER(IntKi),   PARAMETER  :: iJointDir= 6 ! Index in Joints where the joint-direction are stored
  INTEGER(IntKi),   PARAMETER  :: iJointStiff= 9 ! Index in Joints where the joint-stiffness is stored

  ! ID for joint types
  INTEGER(IntKi),   PARAMETER  :: idJointCantilever = 1
  INTEGER(IntKi),   PARAMETER  :: idJointUniversal  = 2
  INTEGER(IntKi),   PARAMETER  :: idJointPin        = 3
  INTEGER(IntKi),   PARAMETER  :: idJointBall       = 4

  ! ID for member types
  INTEGER(IntKi),   PARAMETER  :: idMemberBeamCirc   = 1
  INTEGER(IntKi),   PARAMETER  :: idMemberCable      = 2
  INTEGER(IntKi),   PARAMETER  :: idMemberRigid      = 3
  INTEGER(IntKi),   PARAMETER  :: idMemberBeamArb    = 4

  ! Types of Boundary Conditions
  INTEGER(IntKi),   PARAMETER  :: idBC_Fixed    = 11 ! Fixed BC
  INTEGER(IntKi),   PARAMETER  :: idBC_Internal = 12 ! Free BC
  INTEGER(IntKi),   PARAMETER  :: idBC_Leader   = 13 ! TODO, and maybe "BC" not appropriate here

  ! Types of Static Improvement Methods
  INTEGER(IntKi),   PARAMETER  :: idSIM_None     = 0
  INTEGER(IntKi),   PARAMETER  :: idSIM_Full     = 1
  INTEGER(IntKi)               :: idSIM_Valid(2)  = (/idSIM_None, idSIM_Full/)

  ! Types of Guyan Damping
  INTEGER(IntKi),   PARAMETER  :: idGuyanDamp_None     = 0
  INTEGER(IntKi),   PARAMETER  :: idGuyanDamp_Rayleigh = 1
  INTEGER(IntKi),   PARAMETER  :: idGuyanDamp_66       = 2 
  INTEGER(IntKi)               :: idGuyanDamp_Valid(3) = (/idGuyanDamp_None, idGuyanDamp_Rayleigh, idGuyanDamp_66 /)
  
  INTEGER(IntKi),   PARAMETER  :: SDMaxInpCols    = MAX(JointsCol,InterfCol,MembersCol,PropSetsBCol,PropSetsXCol,COSMsCol,CMassCol)

  ! Output Formats
  INTEGER(IntKi),   PARAMETER  :: idOutputFormatNone  = 0
  INTEGER(IntKi),   PARAMETER  :: idOutputFormatJSON  = 1


  ! Implementation Flags
  LOGICAL, PARAMETER :: DEV_VERSION    = .false.
  LOGICAL, PARAMETER :: BC_Before_CB   = .true.
  LOGICAL, PARAMETER :: ANALYTICAL_LIN = .true.
  LOGICAL, PARAMETER :: GUYAN_RIGID_FLOATING = .true.


CONTAINS
!------------------------------------------------------------------------------------------------------
! --- Helper functions
!------------------------------------------------------------------------------------------------------
!> Maps nodes to elements 
!! allocate NodesConnE and NodesConnN                                                                               
SUBROUTINE NodeCon(Init,p, ErrStat, ErrMsg)
   TYPE(SD_InitType),              INTENT( INOUT ) :: Init
   TYPE(SD_ParameterType),         INTENT( IN    ) :: p
   INTEGER(IntKi),                 INTENT(   OUT ) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(   OUT ) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER(IntKi) :: I,J,K  !counter

   ! The row index is the number of the real node, i.e. ID, 1st col has number of elements attached to node, and 2nd col has element numbers (up to 10)                                    
   CALL AllocAry(Init%NodesConnE, p%nNodes, MaxMemJnt+1,'NodesConnE', ErrStat, ErrMsg); if (ErrStat/=0) return;
   CALL AllocAry(Init%NodesConnN, p%nNodes, MaxMemJnt+2,'NodesConnN', ErrStat, ErrMsg); if (ErrStat/=0) return;
   Init%NodesConnE = 0                                                                                                    
   Init%NodesConnN = -99999 ! Not Used

   ! find the node connectivity, nodes/elements that connect to a common node                                             
   DO I = 1, p%nNodes                                                                                                   
      !Init%NodesConnN(I, 1) = NINT( Init%Nodes(I, 1) )      !This should not be needed, could remove the extra 1st column like for the other array                                                                      
      k = 0                                                                                                               
      DO J = 1, Init%NElem                          !This should be vectorized                                                                      
         IF ( ( NINT(Init%Nodes(I, 1))==p%Elems(J, 2)) .OR. (NINT(Init%Nodes(I, 1))==p%Elems(J, 3) ) ) THEN   !If i-th nodeID matches 1st node or 2nd of j-th element                                                                   
            k = k + 1                                                                                                     
            if (k+1 > MaxMemJnt+1) then 
               CALL SetErrStat(ErrID_Fatal, 'Maximum number of members reached on node number '//trim(Num2LStr(NINT(Init%Nodes(I,1))))//&
                  &' (index in Joint list, not JointID)). The maximum number of member per node is hardcoded to MaxMemJnt='//trim(num2lstr(MaxMemJnt))//&
                  &'. Recompile the code by changing `MaxMemJnt` in SD_FEM.f90.', ErrStat, ErrMsg, 'NodeCon');
               return
            endif
            Init%NodesConnE(I, k + 1) = p%Elems(J, 1)                                                                  
         ENDIF                                                                                                            
      ENDDO                                                                                                               
      Init%NodesConnE(I, 1) = k    !Store how many elements connect i-th node in 2nd column                                                                                       
   ENDDO                            

END SUBROUTINE NodeCon

!----------------------------------------------------------------------------
!> Check if two elements are connected
!! returns true if they are, and return which node (1 or 2) of each element is involved
LOGICAL FUNCTION ElementsConnected(p, ie1, ie2, iWhichNode_e1, iWhichNode_e2)
   TYPE(SD_ParameterType),       INTENT(IN)  :: p
   INTEGER(IntKi),               INTENT(IN)  :: ie1, ie2 ! Indices of elements
   INTEGER(IntKi),               INTENT(OUT) :: iWhichNode_e1, iWhichNode_e2 ! 1 or 2 if node 1 or node 2
   if      ((p%Elems(ie1, 2) == p%Elems(ie2, 2))) then ! node 1 connected to node 1
      iWhichNode_e1=1
      iWhichNode_e2=1
      ElementsConnected=.True.
   else if((p%Elems(ie1, 2) == p%Elems(ie2, 3))) then  ! node 1 connected to node 2
      iWhichNode_e1=1
      iWhichNode_e2=2
      ElementsConnected=.True.
   else if((p%Elems(ie1, 3) == p%Elems(ie2, 2))) then  ! node 2 connected to node 1
      iWhichNode_e1=2
      iWhichNode_e2=1
      ElementsConnected=.True.
   else if((p%Elems(ie1, 3) == p%Elems(ie2, 3))) then  ! node 2 connected to node 2
      iWhichNode_e1=2
      iWhichNode_e2=2
      ElementsConnected=.True.
   else
      ElementsConnected=.False.
      iWhichNode_e1=-1
      iWhichNode_e2=-1
   endif
END FUNCTION ElementsConnected

!> Loop through a list of elements and returns a list of unique joints
TYPE(IList) FUNCTION NodesList(p, Elements)
   use IntegerList, only: init_list, append, find, sort
   use IntegerList, only: print_list
   TYPE(SD_ParameterType),       INTENT(IN)  :: p
   integer(IntKi), dimension(:), INTENT(IN)  :: Elements
   integer(IntKi)  :: ie, ei, j1, j2
   INTEGER(IntKi)  :: ErrStat2
   CHARACTER(ErrMsgLen) :: ErrMsg2

   call init_list(NodesList, 0, 0, ErrStat2, ErrMsg2)
   do ie = 1, size(Elements)
      ei = Elements(ie)  ! Element index
      j1 = p%Elems(ei,2) ! Joint 1 
      j2 = p%Elems(ei,3) ! Joint 2
      ! Append joints indices if not in list already
      if (find(NodesList, j1, ErrStat2, ErrMsg2)<=0) call append(NodesList, j1, ErrStat2, ErrMsg2)
      if (find(NodesList, j2, ErrStat2, ErrMsg2)<=0) call append(NodesList, j2, ErrStat2, ErrMsg2)
      ! Sorting required by find function
      call sort(NodesList, ErrStat2, ErrMsg2)
   enddo
   if (DEV_VERSION) then
      call print_list(NodesList, 'Joint list')
   endif
END FUNCTION NodesList
!------------------------------------------------------------------------------------------------------
!> Returns list of rigid link elements (Er) 
TYPE(IList) FUNCTION RigidLinkElements(Init, p, ErrStat, ErrMsg)
   use IntegerList, only: init_list, append
   use IntegerList, only: print_list
   TYPE(SD_InitType),            INTENT(INOUT) :: Init
   TYPE(SD_ParameterType),       INTENT(INOUT) :: p
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)  :: ie       !< Index on elements
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! --- Establish a list of rigid link elements
   call init_list(RigidLinkElements, 0, 0, ErrStat, ErrMsg);

   do ie = 1, Init%NElem
      if (p%ElemProps(ie)%eType == idMemberRigid) then
         call append(RigidLinkElements, ie, ErrStat, ErrMsg);
      endif
   end do
   if (DEV_VERSION) then
      call print_list(RigidLinkElements,'Rigid element list')
   endif
END FUNCTION RigidLinkElements

!------------------------------------------------------------------------------------------------------
!> Returns true if one of the element connected to the node is a rigid link
LOGICAL FUNCTION NodeHasRigidElem(iJoint, Init, p, ei)
   integer(IntKi),               intent(in)    :: iJoint
   type(SD_InitType),            intent(in)    :: Init
   type(SD_ParameterType),       intent(in)    :: p
   integer(IntKi),               intent(  out) :: ei !< Element index that connects do iJoint rigidly
   ! Local variables
   integer(IntKi) :: ie       !< Loop index on elements

   NodeHasRigidElem = .False. ! default return value
   ! Loop through elements connected to node J 
   do ie = 1, Init%NodesConnE(iJoint, 1)
      ei = Init%NodesConnE(iJoint, ie+1)
      if (p%ElemProps(ei)%eType == idMemberRigid) then
         NodeHasRigidElem = .True.
         return  ! we exit as soon as one rigid member is found
      endif
   enddo
   ei=-1
END FUNCTION NodeHasRigidElem
!------------------------------------------------------------------------------------------------------
!> Returns a rigid body transformation matrix from nDOF to 6 reference DOF: T_ref (6 x nDOF), such that Uref = T_ref.U_subset
!! Typically called to get: 
!!    - the transformation from the interface points to the TP point
!!    - the transformation from the bottom nodes to SubDyn origin (0,0,)
SUBROUTINE RigidTrnsf(Init, p, RefPoint, DOF, nDOF, T_ref, ErrStat, ErrMsg)
   TYPE(SD_InitType),      INTENT(IN   )  :: Init        ! Input data for initialization routine
   TYPE(SD_ParameterType), INTENT(IN   )  :: p        
   REAL(ReKi),             INTENT(IN   )  :: RefPoint(3) ! Coordinate of the reference point 
   INTEGER(IntKi),         INTENT(IN   )  :: nDOF        ! Number of DOFS 
   INTEGER(IntKi),         INTENT(IN   )  :: DOF(nDOF)  ! DOF indices that are used to create the transformation matrix
   REAL(ReKi),             INTENT(  OUT)  :: T_ref(nDOF,6)  ! matrix that relates the subset of DOFs to the reference point
   INTEGER(IntKi),         INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER                                :: I, iDOF, iiDOF, iNode, nDOFPerNode
   REAL(ReKi)                             :: dx, dy, dz
   REAL(ReKi), dimension(6)               :: Line
   ErrStat = ErrID_None
   ErrMsg  = ""
   T_ref(:,:)=0
   DO I = 1, nDOF
      iDOF        = DOF(I) ! DOF index in constrained system
      iNode       = p%DOFred2Nodes(iDOF,1) ! First column is node 
      nDOFPerNode = p%DOFred2Nodes(iDOF,2) ! Second column is number of DOF per node
      iiDOF       = p%DOFred2Nodes(iDOF,3) ! Third column is dof index for this joint (1-6 for cantilever)

      if ((iiDOF<1) .or. (iiDOF>6)) then
         ErrMsg  = 'RigidTrnsf, node DOF number is not valid. DOF:'//trim(Num2LStr(iDOF))//' Node:'//trim(Num2LStr(iNode))//' iiDOF:'//trim(Num2LStr(iiDOF)); ErrStat = ErrID_Fatal
         return
      endif
      if (nDOFPerNode/=6) then
         ErrMsg  = 'RigidTrnsf, node doesnt have 6 DOFs. DOF:'//trim(Num2LStr(iDOF))//' Node:'//trim(Num2LStr(iNode))//' nDOF:'//trim(Num2LStr(nDOFPerNode)); ErrStat = ErrID_Fatal
         return
      endif
      
      dx = Init%Nodes(iNode, 2) - RefPoint(1)
      dy = Init%Nodes(iNode, 3) - RefPoint(2)
      dz = Init%Nodes(iNode, 4) - RefPoint(3)

      CALL RigidTransformationLine(dx,dy,dz,iiDOF,Line) !returns Line
      T_ref(I, 1:6) = Line
   ENDDO
END SUBROUTINE RigidTrnsf

!------------------------------------------------------------------------------------------------------
! --- Main routines, more or less listed in order in which they are called
!------------------------------------------------------------------------------------------------------
!> Reindexing 
! - Removes the notion of "ID" and use Index instead
! - Creates Nodes (use indices instead of ID), similar to Joints array
! - Creates Elems (use indices instead of ID)  similar to Members array
! - Updates Reacts (use indices instead of ID)
! - Updates Interf (use indices instead of ID)
SUBROUTINE SD_ReIndex_CreateNodesAndElems(Init,p, ErrStat, ErrMsg)
   TYPE(SD_InitType),            INTENT(INOUT)  ::Init
   TYPE(SD_ParameterType),       INTENT(INOUT)  ::p
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variable
   INTEGER                       :: I, n, iMem, iNode, JointID, mID, jType, iJoint, iInterf
   INTEGER(IntKi)                :: mType !< Member Type
   CHARACTER(1255)               :: sType !< String for element type
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! TODO See if Elems is actually used elsewhere
   CALL AllocAry(p%Elems         ,Init%NElem,MembersCol        ,'p%Elems'         ,ErrStat2,ErrMsg2); if(Failed())return
   CALL AllocAry(Init%Nodes      ,p%nNodes  ,JointsCol         ,'Init%Nodes'      ,ErrStat2,ErrMsg2); if(Failed())return
   CALL AllocAry(p%NodeID2JointID,p%nNodes                     ,'p%NodeID2JointID',ErrStat2,ErrMsg2); if(Failed())return

   ! --- Initialize Nodes
   Init%Nodes       = -999999 ! Init to unphysical values
   p%NodeID2JointID = -999999
   do I = 1,Init%NJoints
      Init%Nodes(I, 1) = I                                     ! JointID replaced by index I
      Init%Nodes(I, 2:JointsCol) = Init%Joints(I, 2:JointsCol) ! All the rest is copied
      p%NodeID2JointID(I) = Init%Joints(I,1)                   ! JointID
   enddo

   ! --- Re-Initialize Reactions, pointing to index instead of JointID
   do I = 1, p%nNodes_C
      JointID=p%Nodes_C(I,1)
      p%Nodes_C(I,1) = FINDLOCI(Init%Joints(:,1), JointID ) ! Replace JointID with Index
      if (p%Nodes_C(I,1)<=0) then
         CALL Fatal('Reaction joint table: line '//TRIM(Num2LStr(I))//' refers to JointID '//trim(Num2LStr(JointID))//' which is not in the joint list!')
         return
      endif
   enddo

   ! --- Re-Initialize interface joints, pointing to index instead of JointID
   do I = 1, p%nNodes_I
      JointID=p%Nodes_I(I,1)
      p%Nodes_I(I,1) = FINDLOCI(Init%Joints(:,1), JointID )
      if (p%Nodes_I(I,1)<=0) then
         CALL Fatal('Interface joint table: line '//TRIM(Num2LStr(I))//' refers to JointID '//trim(Num2LStr(JointID))//' which is not in the joint list!')
         return
      endif
   enddo

   ! Change numbering in concentrated mass matrix
   do I = 1, Init%NCMass
      JointID = Init%CMass(I,1)
      Init%CMass(I,1) = FINDLOCI(Init%Joints(:,1), JointID )
      if (Init%CMass(I,1)<=0) then
         CALL Fatal('Concentrated mass table: line '//TRIM(Num2LStr(I))//' refers to JointID '//trim(Num2LStr(JointID))//' which is not in the joint list!')
         return
      endif
   enddo

   ! --- Initialize Elems, starting with each member as an element (we'll take NDiv into account later)
   p%Elems = 0
   ! --- Replacing "MemberID"  "JointID", and "PropSetID" by simple index in this tables
   DO iMem = 1, p%NMembers
      mID = Init%Members(iMem, 1)
      ! Column 1  : member index (instead of MemberID)
      p%Elems(iMem,     1)  = iMem
      mType =  Init%Members(iMem, iMType) ! 
      ! Column 2-3: Joint index (instead of JointIDs)
      p%Elems(iMem,     1)  = iMem  ! NOTE: element/member number (not MemberID)
      do iNode=2,3
         JointID = Init%Members(iMem, iNode)
         iJoint = FINDLOCI(Init%Joints(:,1), JointID ) 
         p%Elems(iMem,iNode) = iJoint
         if (p%Elems(iMem,iNode)<=0) then
            CALL Fatal(' MemberID '//TRIM(Num2LStr(mID))//' has JointID'//TRIM(Num2LStr(iNode-1))//' = '// TRIM(Num2LStr(JointID))//' which is not in the joint list!')
            return
         endif
         if (mType==idMemberRigid) then
            ! Check that rigid link are not connected to ball/pin/universal joints
            jType = int(Init%Nodes(iJoint, iJointType))
            if (jType /= idJointCantilever) then
               CALL Fatal('All joints of a rigid link should be cantilever (not ball/pin/universal). The problematic member is MemberID='//TRIM(Num2LStr(mID))//' (which is a rigid link) involving joint JointID='// TRIM(Num2LStr(JointID))// ' (which is not a cantilever joint).')
               return
            endif
            ! Check that rigid links are not connected to the interface
            iInterf = FINDLOCI(p%Nodes_I(:,1), iJoint )
            if (iInterf>=1) then
               CALL WrScr('[WARNING] There might be a bug when rigid links are connected to the interface nodes (mostly if cables are involved). The problematic member is MemberID='//TRIM(Num2LStr(mID))//' (which is a rigid link) involving joint JointID='// TRIM(Num2LStr(JointID))// ' (which is in an interface joint).')
            endif
         endif
      enddo
      ! Column 4-5: PropIndex 1-2 (instead of PropSetID1&2)
      ! NOTE: this index has different meaning depending on the member type !
      DO n=iMProp,iMProp+1

         if (mType==idMemberBeamCirc) then
            sType='Member circular cross-section property'
            p%Elems(iMem,n) = FINDLOCI(Init%PropSetsB(:,1), Init%Members(iMem, n) ) 
         else if (mType==idMemberCable) then
            sType='Cable property'
            p%Elems(iMem,n) = FINDLOCI(Init%PropSetsC(:,1), Init%Members(iMem, n) ) 
         else if (mType==idMemberRigid) then
            sType='Rigid property'
            p%Elems(iMem,n) = FINDLOCI(Init%PropSetsR(:,1), Init%Members(iMem, n) ) 
         else if (mType==idMemberBeamArb) then
            sType='Member arbitrary cross-section property'
            p%Elems(iMem,n) = FINDLOCI(Init%PropSetsX(:,1), Init%Members(iMem, n) )
         else
            ! Should not happen
            print*,'Element type unknown',mType
            STOP
         end if
         ! Test that the two properties match for non-beam 
         if (mType/=idMemberBeamCirc) then
             if (Init%Members(iMem, iMProp)/=Init%Members(iMem, iMProp+1)) then
                ! NOTE: for non circular beams, we could just check that E, rho, G are the same for both properties
                call Fatal('Property IDs should be the same at both joints for arbitrary beams, rigid links, and cables. Check member with ID: '//TRIM(Num2LStr(Init%Members(iMem,1))))
                return
             endif
         endif
         if (p%Elems(iMem,n)<=0) then
            CALL Fatal('For MemberID '//TRIM(Num2LStr(Init%Members(iMem,1)))//', the PropSetID'//TRIM(Num2LStr(n-3))//' is not in the'//trim(sType)//' table!')
            return
         endif
      END DO !n, loop through property ids         
      ! Column 6: member type
      p%Elems(iMem, iMType) = Init%Members(iMem, iMType) ! 
      ! Column 7: member type

      if (p%Elems(iMem,  iMDirCosID)/=-1) then
         p%Elems(iMem,  iMDirCosID) = FINDLOCI(Init%COSMs(:,1), Init%Members(iMem,  iMDirCosID) )
      endif

   END DO !iMem, loop through members
    
   ! TODO in theory, we shouldn't need these anymore
   ! deallocate(Init%Members)
   ! deallocate(Init%Joints)
CONTAINS
   LOGICAL FUNCTION Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_ReIndex_CreateNodesAndElems') 
      Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed
   SUBROUTINE Fatal(ErrMsg_in)
      CHARACTER(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'SD_ReIndex_CreateNodesAndElems');
   END SUBROUTINE Fatal
END SUBROUTINE SD_ReIndex_CreateNodesAndElems

!----------------------------------------------------------------------------
!> Divide (split) members into nDIV element. Only beams are split.
SUBROUTINE SD_Discrt(Init,p, ErrStat, ErrMsg)
   TYPE(SD_InitType),            INTENT(INOUT)  ::Init
   TYPE(SD_ParameterType),       INTENT(INOUT)  ::p
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variable
   INTEGER                       :: I, J, Node1, Node2, Prop1, Prop2   
   INTEGER                       :: NNE      ! number of nodes per element
   INTEGER                       :: MaxNProp
   REAL(ReKi), ALLOCATABLE       :: TempProps(:, :)
   REAL(ReKi), ALLOCATABLE       :: TempPropsX(:, :)
   INTEGER, ALLOCATABLE          :: TempMembers(:, :)
   INTEGER                       :: knode, kelem, kprop, nprop
   INTEGER                       :: iDirCos
   REAL(ReKi)                    :: x1, y1, z1, x2, y2, z2, dx, dy, dz, dd, dt, d1, d2, t1, t2
   LOGICAL                       :: CreateNewProp
   INTEGER(IntKi)                :: nMemberCable, nMemberRigid, nMemberBeamCirc, nMemberBeamArb !< Number of memebers per type
   INTEGER(IntKi)                :: eType !< Element Type
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! number of nodes per element
   IF( ( Init%FEMMod >= 0 ) .and. (Init%FEMMod <= 3) ) THEN
      NNE = 2 
   ELSE
      CALL Fatal('FEMMod '//TRIM(Num2LStr(Init%FEMMod))//' not implemented.'); return
   ENDIF
   
   ! --- Total number of element   
   nMemberBeamCirc = count(Init%Members(:,iMType) == idMemberBeamCirc)
   nMemberCable    = count(Init%Members(:,iMType) == idMemberCable)
   nMemberRigid    = count(Init%Members(:,iMType) == idMemberRigid)
   nMemberBeamArb  = count(Init%Members(:,iMType) == idMemberBeamArb)
   Init%NElem = (nMemberBeamCirc + nMemberBeamArb)*Init%NDiv + nMemberCable + nMemberRigid  ! NOTE: only Beams are divided
   IF ( (nMemberBeamCirc+nMemberRigid+nMemberCable+nMemberBeamArb) /= size(Init%Members,1)) then
      CALL Fatal(' Member list contains an element which is not a beam, a cable or a rigid link'); return
   ENDIF

   ! Total number of nodes - Depends on division and number of nodes per element
   p%nNodes = Init%NJoints + ( Init%NDiv - 1 )*(nMemberBeamCirc) ! TODO add nMemberBeamArb when support for division provided
   
   ! check the number of interior modes
   IF ( p%nDOFM > 6*(p%nNodes - p%nNodes_I - p%nNodes_C) ) THEN
      CALL Fatal(' NModes must be less than or equal to '//TRIM(Num2LStr( 6*(p%nNodes - p%nNodes_I - p%nNodes_C) ))); return
   ENDIF
   
   ! TODO replace this with an integer list!
   CALL AllocAry(Init%MemberNodes,p%NMembers,    Init%NDiv+1,'Init%MemberNodes',ErrStat2, ErrMsg2); if(Failed()) return ! for two-node element only, otherwise the number of nodes in one element is different

   ! --- Reindexing JointsID and MembersID into Nodes and Elems arrays
   ! NOTE: need NNode and NElem 
   CALL SD_ReIndex_CreateNodesAndElems(Init, p, ErrStat2, ErrMsg2);  if(Failed()) return
   
   ! --- Perform some sanity checks (irrespectively of NDiv) Would be better to do that before Reindexing...
    do I = 1, p%NMembers !the first p%NMembers rows of p%Elems contain the element information
       ! Member data
       Node1 = p%Elems(I, 2)
       Node2 = p%Elems(I, 3)
       Prop1 = p%Elems(I, iMProp  )
       Prop2 = p%Elems(I, iMProp+1)
       eType = p%Elems(I, iMType  )

       if ( Node1==Node2 ) THEN
          CALL Fatal(' Same starting and ending node in the member. (See member at position '//trim(num2lstr(I))//' in member list)')
          return
       endif

       if (eType==idMemberBeamCirc) then
          if  ( ( .not. EqualRealNos(Init%PropSetsB(Prop1, 2),Init%PropSetsB(Prop2, 2) ) ) &
           .or. ( .not. EqualRealNos(Init%PropSetsB(Prop1, 3),Init%PropSetsB(Prop2, 3) ) ) &
           .or. ( .not. EqualRealNos(Init%PropSetsB(Prop1, 4),Init%PropSetsB(Prop2, 4) ) ) ) then
             call Fatal(' Material E, G and rho in a member must be the same (See member at position '//trim(num2lstr(I))//' in member list)')
             return
          endif
       else if (eType==idMemberBeamArb) then
          if  (Prop1 /= Prop2 ) then
             call Fatal(' Members using arbitrary cross section properties must have the same properties on both ends. See member at position '//trim(num2lstr(I))//' in member list)')
             return
          endif
          !if  ( ( .not. EqualRealNos(Init%PropSetsX(Prop1, 2),Init%PropSetsX(Prop2, 2) ) ) &
          ! .or. ( .not. EqualRealNos(Init%PropSetsX(Prop1, 3),Init%PropSetsX(Prop2, 3) ) ) &
          ! .or. ( .not. EqualRealNos(Init%PropSetsX(Prop1, 4),Init%PropSetsX(Prop2, 4) ) ) ) then
          !   call Fatal(' Material E, G and rho in a member must be the same (See member at position '//trim(num2lstr(I))//' in member list)')
          !   return
          !endif
       endif ! is beam
    enddo

  
    Init%MemberNodes = 0
    ! --- Setting up MemberNodes (And Elems, Props, Nodes if divisions)
    if (Init%NDiv==1) then
       ! NDiv = 1
       Init%MemberNodes(1:p%NMembers, 1:2) = p%Elems(1:Init%NElem, 2:3) 
       Init%NPropB = Init%NPropSetsB

    else if (Init%NDiv > 1) then

       ! Discretize structure according to NDiv 
       ! - Elems is fully reinitialized, connectivity needs to be done again using SetNewElem
       ! - Nodes are not  reinitialized, but appended to NNodes
       ! 

       ! Initialize Temp arrays that will contain user inputs + input from the subdivided members
       !  We don't know how many properties will be needed, so allocated to size MaxNProp
       ! TODO add Init%NPropSetsX and use PropSetXCol in the future or allocate a new TempProps
       MaxNProp   = Init%NPropSetsB  + Init%NElem*NNE ! Maximum possible number of property sets (temp): This is property set per element node, for all elements (bjj, added Init%NPropSets to account for possibility of entering many unused prop sets)
       CALL AllocAry(TempMembers, p%NMembers,    MembersCol , 'TempMembers', ErrStat2, ErrMsg2); if(Failed()) return
       CALL AllocAry(TempProps,  MaxNProp,      PropSetsBCol,'TempProps',  ErrStat2, ErrMsg2); if(Failed()) return
       TempProps = -9999.
       TempProps(1:Init%NPropSetsB, :) = Init%PropSetsB  
       TempMembers                      = p%Elems(1:p%NMembers,:)
        
       p%Elems(:,:) = -9999. ! Reinitialized. Elements will be ordered by member subdivisions (see setNewElem)

       kelem = 0
       knode = Init%NJoints

       kprop = Init%NPropSetsB

       DO I = 1, p%NMembers !the first p%NMembers rows of p%Elems contain the element information
          ! Member data
          Node1 = TempMembers(I, 2)
          Node2 = TempMembers(I, 3)
          Prop1 = TempMembers(I, iMProp  )
          Prop2 = TempMembers(I, iMProp+1)
          eType = TempMembers(I, iMType  )
          iDirCos = TempMembers(I,  iMDirCosID)
          
          if (eType==idMemberRigid .OR. eType==idMemberCable) then
             ! --- Cables and rigid links are not subdivided and have same prop at nodes
             ! No need to create new properties or new nodes
             Init%MemberNodes(I, 1) = Node1
             Init%MemberNodes(I, 2) = Node2
             kelem = kelem + 1
             CALL SetNewElem(kelem, Node1, Node2, eType, Prop1, Prop1, p, iDirCos)                
             cycle
          endif

          ! --- Subdivision of beams
          Init%MemberNodes(I,           1) = Node1
          Init%MemberNodes(I, Init%NDiv+1) = Node2

          x1 = Init%Nodes(Node1, 2)
          y1 = Init%Nodes(Node1, 3)
          z1 = Init%Nodes(Node1, 4)

          x2 = Init%Nodes(Node2, 2)
          y2 = Init%Nodes(Node2, 3)
          z2 = Init%Nodes(Node2, 4)
          
          dx = ( x2 - x1 )/Init%NDiv
          dy = ( y2 - y1 )/Init%NDiv
          dz = ( z2 - z1 )/Init%NDiv

          if (eType == idMemberBeamCirc) then

            d1 = TempProps(Prop1, 5)
            t1 = TempProps(Prop1, 6)

            d2 = TempProps(Prop2, 5)
            t2 = TempProps(Prop2, 6)
            
            dd = ( d2 - d1 )/Init%NDiv
            dt = ( t2 - t1 )/Init%NDiv

             ! If both dd and dt are 0, no interpolation is needed, and we can use the same property set for new nodes/elements. otherwise we'll have to create new properties for each new node
           
            CreateNewProp = .NOT. ( EqualRealNos( dd , 0.0_ReKi ) .AND.  EqualRealNos( dt , 0.0_ReKi ) ) 

          elseif (eType == idMemberBeamArb) then
    
            CreateNewProp = .FALSE.
            CALL WrScr('[WARNING] Members with non-circular cross-sections are currently not divided (member at position '//TRIM(Num2LStr(I))//' ).')
         
          endif

          ! node connect to Node1
          knode = knode + 1
          Init%MemberNodes(I, 2) = knode
          CALL SetNewNode(knode, x1+dx, y1+dy, z1+dz, Init); if (ErrStat>ErrID_None) return;
          
          IF ( CreateNewProp ) THEN   
               ! create a new property set 
               kprop = kprop + 1
               !                  k,  E1,                  G1,                  rho1,                d,     t,    
               CALL SetNewProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3), TempProps(Prop1, 4), d1+dd, t1+dt, TempProps)           
               kelem = kelem + 1
               CALL SetNewElem(kelem, Node1, knode, eType, Prop1, kprop, p, iDirCos); if (ErrStat>ErrID_None) return;
               nprop = kprop
          ELSE
               kelem = kelem + 1
               CALL SetNewElem(kelem, Node1, knode, eType, Prop1, Prop1, p, iDirCos); if (ErrStat>ErrID_None) return;             
               nprop = Prop1 
          ENDIF
          
          ! interior nodes
          DO J = 2, (Init%NDiv-1)
             knode = knode + 1
             Init%MemberNodes(I, J+1) = knode

             CALL SetNewNode(knode, x1 + J*dx, y1 + J*dy, z1 + J*dz, Init) ! Set Init%Nodes(knode,:)
             
             IF ( CreateNewProp ) THEN   
                  ! create a new property set 
                  kprop = kprop + 1
                  !                  k,  E1,                  G1,                  rho1,                     d,          t
                  CALL SetNewProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3), Init%PropSetsB(Prop1, 4), d1 + J*dd, t1 + J*dt,  TempProps)           
                  kelem = kelem + 1
                  CALL SetNewElem(kelem, knode-1, knode, eType, nprop, kprop, p, iDirCos); if (ErrStat>ErrID_None) return;
                  nprop = kprop
             ELSE
                  kelem = kelem + 1
                  CALL SetNewElem(kelem, knode-1, knode, eType, nprop, nprop, p, iDirCos); if (ErrStat>ErrID_None) return;
             ENDIF
          ENDDO
          
          ! the element connect to Node2
          kelem = kelem + 1
          CALL SetNewElem(kelem, knode, Node2, eType, nprop, Prop2, p, iDirCos); if (ErrStat>ErrID_None) return;
       ENDDO ! loop over all members
       !
       Init%NPropB = kprop
       if(knode/=size(Init%Nodes,1)) then
          call Fatal('Implementation error. Number of nodes wrongly estimated.');return
       endif
       if(kelem/=size(p%Elems,1)) then
          call Fatal('Implementation error. Number of elements wrongly estimated.');return
       endif

    ENDIF ! if NDiv is greater than 1

    ! set the props in Init
    CALL AllocAry(Init%PropsB, Init%NPropB, PropSetsBCol, 'Init%PropsBeams', ErrStat2, ErrMsg2); if(Failed()) return

    if (Init%NDiv==1) then
       Init%PropsB(1:Init%NPropB, 1:PropSetsBCol) = Init%PropSetsB(1:Init%NPropB, 1:PropSetsBCol)
    else if (Init%NDiv>1) then
       Init%PropsB(1:Init%NPropB, 1:PropSetsBCol) = TempProps(1:Init%NPropB, 1:PropSetsBCol)
    endif

    ! --- Cables and rigid link properties (these cannot be subdivided, so direct copy of inputs)
    Init%NPropC = Init%NPropSetsC
    Init%NPropR = Init%NPropSetsR
    CALL AllocAry(Init%PropsC, Init%NPropC, PropSetsCCol, 'Init%PropsCable', ErrStat2, ErrMsg2); if(Failed()) return
    CALL AllocAry(Init%PropsR, Init%NPropR, PropSetsRCol, 'Init%PropsRigid', ErrStat2, ErrMsg2); if(Failed()) return
    Init%PropsC(1:Init%NPropC, 1:PropSetsCCol) = Init%PropSetsC(1:Init%NPropC, 1:PropSetsCCol)
    Init%PropsR(1:Init%NPropR, 1:PropSetsRCol) = Init%PropSetsR(1:Init%NPropR, 1:PropSetsRCol)

    CALL CleanUp_Discrt()

CONTAINS
   LOGICAL FUNCTION Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Discrt') 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) CALL CleanUp_Discrt()
   END FUNCTION Failed

   SUBROUTINE Fatal(ErrMsg_in)
      CHARACTER(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'SD_Discrt');
      CALL CleanUp_Discrt()
   END SUBROUTINE Fatal

   SUBROUTINE CleanUp_Discrt()
      ! deallocate temp matrices
      IF (ALLOCATED(TempProps))   DEALLOCATE(TempProps)
      IF (ALLOCATED(TempMembers)) DEALLOCATE(TempMembers)
   END SUBROUTINE CleanUp_Discrt

   !> Set properties of node k
   SUBROUTINE SetNewNode(k, x, y, z, Init)
      TYPE(SD_InitType),      INTENT(INOUT) :: Init
      INTEGER,                INTENT(IN)    :: k
      REAL(ReKi),             INTENT(IN)    :: x, y, z
      if (k>size(Init%Nodes,1)) then
         call Fatal('Implementation Error. Attempt to add more node than space allocated.');
         return
      endif
      Init%Nodes(k, 1)                     = k
      Init%Nodes(k, 2)                     = x
      Init%Nodes(k, 3)                     = y
      Init%Nodes(k, 4)                     = z
      Init%Nodes(k, iJointType)            = idJointCantilever ! Note: all added nodes are Cantilever
      ! Properties below are for non-cantilever joints
      Init%Nodes(k, iJointDir:iJointDir+2) = 0.0_ReKi ! NOTE: irrelevant for cantilever nodes
      Init%Nodes(k, iJointStiff)           = 0.0_ReKi ! NOTE: irrelevant for cantilever nodes
   END SUBROUTINE SetNewNode
   
   !> Set properties of element k
   SUBROUTINE SetNewElem(k, n1, n2, etype, p1, p2, p, iDirCos)
      INTEGER,                INTENT(IN   )   :: k
      INTEGER,                INTENT(IN   )   :: n1
      INTEGER,                INTENT(IN   )   :: n2
      INTEGER,                INTENT(IN   )   :: eType
      INTEGER,                INTENT(IN   )   :: p1
      INTEGER,                INTENT(IN   )   :: p2
      INTEGER,                INTENT(IN   )   :: iDirCos
      TYPE(SD_ParameterType), INTENT(INOUT)   :: p
      if (k>size(p%Elems,1)) then
         call Fatal('Implementation Error. Attempt to add more element than space allocated.');
         return
      endif
      p%Elems(k, 1)        = k
      p%Elems(k, 2)        = n1
      p%Elems(k, 3)        = n2
      p%Elems(k, iMProp  ) = p1
      p%Elems(k, iMProp+1) = p2
      p%Elems(k, iMType)   = eType
      p%Elems(k,  iMDirCosID)   = iDirCos
   END SUBROUTINE SetNewElem

   !> Set material properties of element k,  NOTE: this is only for a beam
   SUBROUTINE SetNewProp(k, E, G, rho, d, t, TempProps)
      INTEGER   , INTENT(IN)   :: k
      REAL(ReKi), INTENT(IN)   :: E, G, rho, d, t
      REAL(ReKi), INTENT(INOUT):: TempProps(:, :)
      if (k>size(TempProps,1)) then
         call Fatal('Implementation Error. Attempt to add more properties than space allocated.');
         return
      endif
      TempProps(k, 1) = k
      TempProps(k, 2) = E
      TempProps(k, 3) = G
      TempProps(k, 4) = rho
      TempProps(k, 5) = d
      TempProps(k, 6) = t
   END SUBROUTINE SetNewProp

END SUBROUTINE SD_Discrt


!> Store relative vector between nodes and TP point, to later compute Guyan rigid body motion
subroutine StoreNodesRelPos(Init, p, ErrStat, ErrMsg)
   type(SD_InitType),      intent(in   ) :: Init
   type(SD_ParameterType), intent(inout) :: p
   integer(IntKi),         intent(out)   :: ErrStat     ! Error status of the operation
   character(*),           intent(out)   :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   integer(Intki)       :: i
   integer(IntKi)       :: ErrStat2
   character(ErrMsgLen) :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! NOTE: using efficient memory order
   call AllocAry(p%DP0, 3, size(Init%Nodes,1), 'DP0', ErrStat2, ErrMsg2); if(Failed()) return

   do i = 1, size(Init%Nodes,1)
      p%DP0(1, i) = Init%Nodes(i, 2) - Init%TP_RefPoint(1)
      p%DP0(2, i) = Init%Nodes(i, 3) - Init%TP_RefPoint(2)
      p%DP0(3, i) = Init%Nodes(i, 4) - Init%TP_RefPoint(3)
   enddo

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, Errstat, ErrMsg, 'StoreNodesRelPos') 
      failed =  ErrStat >= AbortErrLev
   end function Failed
end subroutine StoreNodesRelPos



!------------------------------------------------------------------------------------------------------
!> Set Element properties p%ElemProps, different properties are set depening on element type..
SUBROUTINE SetElementProperties(Init, p, ErrStat, ErrMsg)
   TYPE(SD_InitType),            INTENT(IN   ) :: Init
   TYPE(SD_ParameterType),       INTENT(INOUT) :: p
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER                  :: I
   INTEGER                  :: N1, N2     ! starting node and ending node in the element
   INTEGER                  :: P1, P2     ! property set numbers for starting and ending nodes
   INTEGER                  :: iDirCos
   REAL(ReKi)               :: D1, D2, t1, t2, E, G, rho ! properties of a section
   REAL(FEKi)               :: DirCos(3, 3)              ! direction cosine matrices
   REAL(ReKi)               :: L                         ! length of the element
   REAL(ReKi)               :: r1, r2, t, Iyy, Jzz, Ixx, A, kappa, kappa_x, kappa_y, nu, ratioSq, D_inner, D_outer
   LOGICAL                  :: shear
   INTEGER(IntKi)           :: eType !< Member type
   REAL(ReKi)               :: Point1(3), Point2(3) ! (x,y,z) positions of two nodes making up an element
   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(ErrMsgLen)     :: ErrMsg2
   ErrMsg  = ""
   ErrStat = ErrID_None
   
   ALLOCATE( p%ElemProps(Init%NElem), STAT=ErrStat2); ErrMsg2='Error allocating p%ElemProps'
   if(Failed()) return
   
   ! Loop over all elements and set ElementProperties
   do I = 1, Init%NElem
      N1 = p%Elems(I, 2)
      N2 = p%Elems(I, 3)
      
      P1    = p%Elems(I, iMProp  )
      P2    = p%Elems(I, iMProp+1)
      eType = p%Elems(I, iMType)
      iDirCos = p%Elems(I,  iMDirCosID)

      ! --- Properties common to all element types: L, DirCos (and Area and rho)
      Point1 = Init%Nodes(N1,2:4)
      Point2 = Init%Nodes(N2,2:4)

      if (iDirCos/=-1) then
         CALL GetDirCos(Point1, Point2, DirCos, L, ErrStat2, ErrMsg2); if(Failed()) return ! sets L
         
         ! overwrites direction cosines
         DirCos(1, 1) =  Init%COSMs(iDirCos, 2)
         DirCos(2, 1) =  Init%COSMs(iDirCos, 3)
         DirCos(3, 1) =  Init%COSMs(iDirCos, 4)
         DirCos(1, 2) =  Init%COSMs(iDirCos, 5)
         DirCos(2, 2) =  Init%COSMs(iDirCos, 6)
         DirCos(3, 2) =  Init%COSMs(iDirCos, 7)
         DirCos(1, 3) =  Init%COSMs(iDirCos, 8)
         DirCos(2, 3) =  Init%COSMs(iDirCos, 9)
         DirCos(3, 3) =  Init%COSMs(iDirCos, 10)

      else
         CALL GetDirCos(Point1, Point2, DirCos, L, ErrStat2, ErrMsg2); if(Failed()) return ! L and DirCos
      endif



      p%ElemProps(i)%eType  = eType
      p%ElemProps(i)%Length = L
      p%ElemProps(i)%DirCos = DirCos

      ! Init to excessive values to detect any issue
      p%ElemProps(i)%Ixx     = -9.99e+36
      p%ElemProps(i)%Iyy     = -9.99e+36
      p%ElemProps(i)%Jzz     = -9.99e+36
      p%ElemProps(i)%Kappa_x   = -9.99e+36
      p%ElemProps(i)%Kappa_y   = -9.99e+36
      p%ElemProps(i)%YoungE  = -9.99e+36
      p%ElemProps(i)%ShearG  = -9.99e+36
      p%ElemProps(i)%Area    = -9.99e+36
      p%ElemProps(i)%Rho     = -9.99e+36
      p%ElemProps(i)%T0      = -9.99e+36

      ! --- Properties that are specific to some elements
      if (eType==idMemberBeamCirc) then
         E   = Init%PropsB(P1, 2) ! TODO E2 
         G   = Init%PropsB(P1, 3) ! TODO G2
         rho = Init%PropsB(P1, 4) ! TODO rho2
         D1  = Init%PropsB(P1, 5)
         t1  = Init%PropsB(P1, 6)
         D2  = Init%PropsB(P2, 5)
         t2  = Init%PropsB(P2, 6)
         r1 = 0.25*(D1 + D2)
         t  = 0.5*(t1+t2)
         if ( EqualRealNos(t, 0.0_ReKi) ) then
            r2 = 0
         else
            r2 = r1 - t
         endif
         A = Pi_D*(r1*r1-r2*r2)
         Ixx = 0.25*Pi_D*(r1**4-r2**4)
         Iyy = Ixx
         Jzz = 2.0*Ixx
         
         if( Init%FEMMod == 1 ) then ! uniform Euler-Bernoulli
            Shear = .false.
            kappa = 0
         elseif( Init%FEMMod == 3 ) then ! uniform Timoshenko
            Shear = .true.
          ! kappa = 0.53            
            ! equation 13 (Steinboeck et al) in SubDyn Theory Manual 
            nu = E / (2.0_ReKi*G) - 1.0_ReKi
            D_outer = 2.0_ReKi * r1  ! average (outer) diameter
            D_inner = D_outer - 2*t  ! remove 2x thickness to get inner diameter
            ratioSq = ( D_inner / D_outer)**2
            kappa =   ( 6.0 * (1.0 + nu) **2 * (1.0 + ratioSq)**2 ) &
                    / ( ( 1.0 + ratioSq )**2 * ( 7.0 + 14.0*nu + 8.0*nu**2 ) + 4.0 * ratioSq * ( 5.0 + 10.0*nu + 4.0 *nu**2 ) )
         endif
         ! Storing Beam specific properties
         p%ElemProps(i)%Ixx    = Ixx
         p%ElemProps(i)%Iyy    = Iyy
         p%ElemProps(i)%Jzz    = Jzz
         p%ElemProps(i)%Shear  = Shear
         p%ElemProps(i)%Kappa_x  = kappa
         p%ElemProps(i)%Kappa_y  = kappa
         p%ElemProps(i)%YoungE = E
         p%ElemProps(i)%ShearG = G
         p%ElemProps(i)%Area   = A
         p%ElemProps(i)%Rho    = rho
         p%ElemProps(i)%D      = (/D1, D2/)

      else if (eType==idMemberBeamArb) then

         p%ElemProps(i)%eType  = 1
         if( Init%FEMMod == 1 ) then ! uniform Euler-Bernoulli
            Shear = .false.
         elseif( Init%FEMMod == 3 ) then ! uniform Timoshenko
            Shear = .true.
         endif
         ! Storing Beam specific properties
         ! Here we are averaging the values at both extremities which is different from what is done for regular beams.
         ! The averaging should have no effect for the material properties because the beam is assumed to be isotropic (E, G, rho constant).
         E   = (Init%PropSetsX(P1, 2) + Init%PropSetsX(P2, 2)) / 2
         G   = (Init%PropSetsX(P1, 3) + Init%PropSetsX(P2, 3)) / 2
         rho = (Init%PropSetsX(P1, 4) + Init%PropSetsX(P2, 4)) / 2
         ! Averaging will have an impact on geometry, shear and inertia
         ! but we are currently forcing the property ID to be the same, so no effect.
         A   = (Init%PropSetsX(P1, 5) + Init%PropSetsX(P2, 5)) / 2
         Kappa_x   = (Init%PropSetsX(P1, 6) + Init%PropSetsX(P2, 6)) / 2 / A
         Kappa_y   = (Init%PropSetsX(P1, 7) + Init%PropSetsX(P2, 7)) / 2 / A
         Ixx   = (Init%PropSetsX(P1, 8) + Init%PropSetsX(P2, 8)) / 2
         Iyy   = (Init%PropSetsX(P1, 9) + Init%PropSetsX(P2, 9)) / 2
         Jzz   = (Init%PropSetsX(P1, 10) + Init%PropSetsX(P2, 10)) / 2
         D1 = 2._ReKi*(A/PI)**0.5 !Approximation, this value should not be used
         D2 = D1

         p%ElemProps(i)%Ixx    = Ixx
         p%ElemProps(i)%Iyy    = Iyy
         p%ElemProps(i)%Jzz    = Jzz
         p%ElemProps(i)%Shear  = Shear
         p%ElemProps(i)%Kappa_x  = kappa_x
         p%ElemProps(i)%Kappa_y  = kappa_y
         p%ElemProps(i)%YoungE = E
         p%ElemProps(i)%ShearG = G
         p%ElemProps(i)%Area   = A
         p%ElemProps(i)%Rho    = rho
         p%ElemProps(i)%D      = (/D1, D2/)

      else if (eType==idMemberCable) then
         if (DEV_VERSION) then
            print*,'Member',I,'is a cable'
         endif
         p%ElemProps(i)%Area   = 1                       ! Arbitrary set to 1
         p%ElemProps(i)%YoungE = Init%PropsC(P1, 2)/1    ! Young's modulus, E=EA/A  [N/m^2]
         p%ElemProps(i)%Rho    = Init%PropsC(P1, 3)      ! Material density [kg/m3]
         p%ElemProps(i)%T0     = Init%PropsC(P1, 4)      ! Pretension force [N]
         p%ElemProps(i)%D      = min(sqrt(1/Pi)*4, L*0.05) ! For plotting only

      else if (eType==idMemberRigid) then
         if (DEV_VERSION) then
            print*,'Member',I,'is a rigid link'
         endif
         p%ElemProps(i)%Area   = 1                  ! Arbitrary set to 1
         p%ElemProps(i)%Rho    = Init%PropsR(P1, 2)
         p%ElemProps(i)%D      = min(sqrt(1/Pi)*4, L*0.05) ! For plotting only

      else
         ! Should not happen
         print*,'Element type unknown',eType
         STOP
      end if
   enddo ! I end loop over elements
CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SetElementProperties') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed
END SUBROUTINE SetElementProperties 


!> Distribute global DOF indices corresponding to Nodes, Elements, BCs, Reactions
!! For Cantilever Joint -> Condensation into 3 translational and 3 rotational DOFs
!! For other joint type -> Condensation of the 3 translational DOF
!!                      -> Keeping 3 rotational DOF for each memeber connected to the joint
SUBROUTINE DistributeDOF(Init, p, ErrStat, ErrMsg)
   use IntegerList, only: init_list, len
   TYPE(SD_InitType),            INTENT(INOUT) :: Init
   TYPE(SD_ParameterType),       INTENT(INOUT) :: p
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   integer(IntKi) :: iNode, k
   integer(IntKi) :: iPrev ! Cumulative counter over the global DOF
   integer(IntKi) :: iElem ! 
   integer(IntKi) :: idElem
   integer(IntKi) :: nRot ! Number of rotational DOFs (multiple of 3) to be used at the joint
   integer(IntKi) :: iOff ! Offset, 0 or 6, depending if node 1 or node 2
   integer(IntKi), dimension(6) :: DOFNode_Old
   integer(IntKi)           :: ErrStat2
   character(ErrMsgLen)     :: ErrMsg2
   ErrMsg  = ""
   ErrStat = ErrID_None

   allocate(p%NodesDOF(1:p%nNodes), stat=ErrStat2)
   ErrMsg2="Error allocating NodesDOF"
   if(Failed()) return

   call AllocAry(p%ElemsDOF, 12, Init%NElem, 'ElemsDOF', ErrStat2, ErrMsg2); if(Failed()) return;
   p%ElemsDOF=-9999

   iPrev =0
   do iNode = 1, p%nNodes
      ! --- Distribute to joints iPrev + 1:6, or, iPrev + 1:(3+3m)
      if (int(Init%Nodes(iNode,iJointType)) == idJointCantilever ) then
         nRot=3
      else
         nRot= 3*Init%NodesConnE(iNode,1) ! Col1: number of elements connected to this joint
      endif
      call init_list(p%NodesDOF(iNode), 3+nRot, iPrev, ErrStat2, ErrMsg2)
      p%NodesDOF(iNode)%List(1:(3+nRot)) = (/ ((iElem+iPrev), iElem=1,3+nRot) /)

      ! --- Distribute to members
      do iElem = 1, Init%NodesConnE(iNode,1) ! members connected to joint iJ
         idElem = Init%NodesConnE(iNode,iElem+1)
         if (iNode == p%Elems(idElem, 2)) then ! Current joint is Elem node 1
            iOff = 0
         else                              ! Current joint is Elem node 2
            iOff = 6
         endif
         p%ElemsDOF(iOff+1:iOff+3, idElem) =  p%NodesDOF(iNode)%List(1:3)
         if (int(Init%Nodes(iNode,iJointType)) == idJointCantilever ) then
            p%ElemsDOF(iOff+4:iOff+6, idElem) = p%NodesDOF(iNode)%List(4:6)
         else
            p%ElemsDOF(iOff+4:iOff+6, idElem) = p%NodesDOF(iNode)%List(3*iElem+1:3*iElem+3)   
         endif
      enddo ! iElem, loop on members connect to joint
      iPrev = iPrev + len(p%NodesDOF(iNode))
   enddo ! iNode, loop on joints

   ! --- Safety check
   if (any(p%ElemsDOF<0)) then
      ErrStat=ErrID_Fatal
      ErrMsg ="Implementation error in Distribute DOF, some member DOF were not allocated"
   endif

   ! --- Safety check (backward compatibility, only valid if all joints are Cantilever)
   if (p%nNodes == count( Init%Nodes(:, iJointType) == idJointCantilever)) then
      do idElem = 1, Init%NElem
         iNode = p%Elems(idElem, 2)
         DOFNode_Old= (/ ((iNode*6-5+k), k=0,5) /)
         if ( any( (p%ElemsDOF(1:6, idElem) /= DOFNode_Old)) ) then
            ErrStat=ErrID_Fatal
            ErrMsg ="Implementation error in Distribute DOF, DOF indices have changed for iElem="//trim(Num2LStr(idElem))
            return
         endif
      enddo
   else
      ! Safety check does not apply if some joints are non-cantilever
   endif

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SetElementProperties') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed

END SUBROUTINE DistributeDOF


!> Checks reaction BC, adn remap 0s and 1s 
SUBROUTINE CheckBCs(p, ErrStat, ErrMsg)
   TYPE(SD_ParameterType),INTENT(INOUT) :: p
   INTEGER(IntKi),        INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),          INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi) :: I, J, iNode
   ErrMsg  = ""
   ErrStat = ErrID_None
   DO I = 1, p%nNodes_C
      iNode = p%Nodes_C(I,1) ! Node index
      DO J = 1, 6
         if (p%Nodes_C(I,J+1)==1) then ! User input 1=Constrained/Fixed (should be eliminated)
            p%Nodes_C(I, J+1)       = idBC_Fixed
         else if (p%Nodes_C(I,J+1)==0) then ! User input 0=Free, fill be part of Internal DOF
            p%Nodes_C(I, J+1)       = idBC_Internal
         else if (p%Nodes_C(I,J+1)==2) then ! User input 2=Leader DOF
            p%Nodes_C(I, J+1)       = idBC_Leader
            ErrStat=ErrID_Fatal
            ErrMsg='BC 2 not allowed for now, node '//trim(Num2LStr(iNode))
         else
            ErrStat=ErrID_Fatal
            ErrMsg='Wrong boundary condition input for reaction node '//trim(Num2LStr(iNode))
         endif
      ENDDO
   ENDDO
END SUBROUTINE CheckBCs

!> Check interface inputs, and remap 0s and 1s 
SUBROUTINE CheckIntf(p, ErrStat, ErrMsg)
   TYPE(SD_ParameterType),INTENT(INOUT) :: p
   INTEGER(IntKi),        INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),          INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi) :: I, J, iNode
   ErrMsg  = ""
   ErrStat = ErrID_None
   DO I = 1, p%nNodes_I
      iNode = p%Nodes_I(I,1) ! Node index
      DO J = 1, 6 ! ItfTDXss    ItfTDYss    ItfTDZss    ItfRDXss    ItfRDYss    ItfRDZss
         if     (p%Nodes_I(I,J+1)==1) then ! User input 1=Leader DOF
            p%Nodes_I(I,J+1)          = idBC_Leader
         elseif (p%Nodes_I(I,J+1)==0) then ! User input 0=Fixed DOF
            p%Nodes_I(I,J+1)          = idBC_Fixed
            ErrStat = ErrID_Fatal
            ErrMsg  = 'Fixed boundary condition not yet supported for interface nodes, node:'//trim(Num2LStr(iNode))
         else
            ErrStat = ErrID_Fatal
            ErrMsg  = 'Wrong boundary condition input for interface node'//trim(Num2LStr(iNode))
         endif
      ENDDO
   ENDDO
END SUBROUTINE CheckIntf


!------------------------------------------------------------------------------------------------------
!> Assemble stiffness and mass matrix, and gravity force vector
SUBROUTINE AssembleKM(Init, p, ErrStat, ErrMsg)
   TYPE(SD_InitType),            INTENT(INOUT) :: Init
   TYPE(SD_ParameterType),       INTENT(INOUT) :: p
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER                  :: I, J, K
   INTEGER                  :: iGlob
   REAL(FEKi)               :: Ke(12,12), Me(12, 12), FGe(12) ! element stiffness and mass matrices gravity force vector
   REAL(FEKi)               :: FCe(12) ! Pretension force from cable element
   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(ErrMsgLen)     :: ErrMsg2
   INTEGER(IntKi)           :: iNode !< Node index
   integer(IntKi), dimension(12) :: IDOF !  12 DOF indices in global unconstrained system
   real(ReKi), dimension(6,6) :: M66  ! Mass matrix of an element node
   real(ReKi) :: m, x, y, z, Jxx, Jyy, Jzz, Jxy, Jxz, Jyz
   INTEGER    :: jGlob, kGlob
   ErrMsg  = ""
   ErrStat = ErrID_None
   
   ! total unconstrained degrees of freedom of the system 
   p%nDOF = nDOF_Unconstrained()
   if (DEV_VERSION) then
      print*,'nDOF_unconstrained:',p%nDOF, ' (if all Cantilever, it would be: ',6*p%nNodes,')'
   endif

   CALL AllocAry( Init%K, p%nDOF, p%nDOF , 'Init%K',  ErrStat2, ErrMsg2); if(Failed()) return; ! system stiffness matrix 
   CALL AllocAry( Init%M, p%nDOF, p%nDOF , 'Init%M',  ErrStat2, ErrMsg2); if(Failed()) return; ! system mass matrix 
   CALL AllocAry( p%FG,   p%nDOF,          'p%FG'  ,  ErrStat2, ErrMsg2); if(Failed()) return; ! system gravity force vector 
   Init%K  = 0.0_FEKi
   Init%M  = 0.0_FEKi
   p%FG    = 0.0_FEKi

   ! loop over all elements, compute element matrices and assemble into global matrices
   DO i = 1, Init%NElem
      ! --- Element Me,Ke,Fg, Fce
      CALL ElemM(p%ElemProps(i), Me)
      CALL ElemK(p%ElemProps(i), Ke)
      CALL ElemF(p%ElemProps(i), Init%g, FGe, FCe)

      ! --- Assembly in global unconstrained system
      IDOF = p%ElemsDOF(1:12, i)
      p%FG     ( IDOF )  = p%FG( IDOF )   + FGe(1:12)+ FCe(1:12) ! Note: gravity and pretension cable forces
      Init%K(IDOF, IDOF) = Init%K( IDOF, IDOF) + Ke(1:12,1:12)
      Init%M(IDOF, IDOF) = Init%M( IDOF, IDOF) + Me(1:12,1:12)
   ENDDO
      
   ! Add concentrated mass to mass matrix
   DO I = 1, Init%nCMass
      iNode = NINT(Init%CMass(I, 1)) ! Note index where concentrated mass is to be added
      ! Safety check (otherwise we might have more than 6 DOF)
      if (Init%Nodes(iNode,iJointType) /= idJointCantilever) then
         ErrMsg2='Concentrated mass is only for cantilever joints. Problematic node: '//trim(Num2LStr(iNode)); ErrStat2=ErrID_Fatal;
         if(Failed()) return
      endif
      ! Mass matrix of a rigid body
      M66 = 0.0_ReKi
      m   = Init%CMass(I,2)
      Jxx = Init%CMass(I,3 ); Jxy = Init%CMass(I,6 ); x = Init%CMass(I,9 );
      Jyy = Init%CMass(I,4 ); Jxz = Init%CMass(I,7 ); y = Init%CMass(I,10);
      Jzz = Init%CMass(I,5 ); Jyz = Init%CMass(I,8 ); z = Init%CMass(I,11);
      call rigidBodyMassMatrix(m, Jxx, Jyy, Jzz, Jxy, Jxz, Jyz, x, y, z, M66)
      ! Adding
      DO J = 1, 6
         jGlob = p%NodesDOF(iNode)%List(J)
         DO K = 1, 6
            kGlob = p%NodesDOF(iNode)%List(K)
            Init%M(jGlob, kGlob) = Init%M(jGlob, kGlob) + M66(J,K)
         ENDDO
      ENDDO
   ENDDO ! Loop on concentrated mass

   ! Add concentrated mass induced gravity force
   DO I = 1, Init%nCMass
       iNode = NINT(Init%CMass(I, 1)) ! Note index where concentrated mass is to be added
       iGlob = p%NodesDOF(iNode)%List(3) ! uz
       p%FG(iGlob) = p%FG(iGlob) - Init%CMass(I, 2)*Init%g 
   ENDDO

   CALL CleanUp_AssembleKM()
   
CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call Cleanup_AssembleKM()
   END FUNCTION Failed
   
   SUBROUTINE Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'AssembleKM');
      CALL CleanUp_AssembleKM()
   END SUBROUTINE Fatal

   SUBROUTINE CleanUp_AssembleKM()
      !pass
   END SUBROUTINE CleanUp_AssembleKM

   INTEGER(IntKi) FUNCTION nDOF_Unconstrained()
      integer(IntKi) :: i
      integer(IntKi) :: m
      nDOF_Unconstrained=0
      do i = 1,p%nNodes
         if (int(Init%Nodes(i,iJointType)) == idJointCantilever ) then
            nDOF_Unconstrained = nDOF_Unconstrained + 6
         else
            m = Init%NodesConnE(i,1) ! Col1: number of elements connected to this joint
            nDOF_Unconstrained = nDOF_Unconstrained + 3 + 3*m
         endif
      end do
   END FUNCTION
   
END SUBROUTINE AssembleKM

!> Map control cable index to control channel index
!! Also set the InitOut%CableCChanRqst logical array to indicate which channels were requested
!!    The array element is set to true for the corresponding requested channel (SD does not need
!!    to request a contiguous block of channels)
subroutine ControlCableMapping(Init, uInit, p, InitOut, ErrStat, ErrMsg)
   type(SD_InitType),            intent(in   ) :: Init        !< init
   type(SD_InputType),           intent(inout) :: uInit       !< init input guess
   type(SD_ParameterType),       intent(inout) :: p           !< param
   type(SD_InitOutputType),      intent(inout) :: InitOut     !< Output for initialization routine
   integer(IntKi),               intent(  out) :: ErrStat     !< Error status of the operation
   character(*),                 intent(  out) :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)           :: i, nCC, idCProp, iElem !< index, number of controlable cables, id of Cable Prop
   integer(IntKi)           ::  maxCC                 !< max control chan number
   integer(IntKi)           :: ErrStat2
   character(ErrMsgLen)     :: ErrMsg2
   ErrMsg  = ""
   ErrStat = ErrID_None

   ! --- Count number of Controllable cables
   nCC = 0
   maxCC = 0
   do i = 1, size(p%ElemProps)
      if (p%ElemProps(i)%eType==idMemberCable) then
         idCProp= p%Elems(i,iMProp)
         if (NINT(Init%PropsC(idCProp, 5 ))>0) then
            !print*,'Cable Element',i,'controllable with channel',Init%PropsC(idCProp, 5 )
            nCC=nCC+1
            maxCC = max( maxCC, NINT(Init%PropsC(idCProp,5)) )
         endif
      endif
   enddo
   if (nCC>0) then
      call WrScr('Number of controllable cables: '//trim(num2lstr(nCC)))
   endif
   call AllocAry( p%CtrlElem2Channel, nCC, 2, 'p%CtrlElem2Channel', ErrStat2, ErrMsg2); if(Failed()) return; ! Constant cable force

   ! --- create array for telling calling code which channels are requested  -- leave unallocated if none found!!!
   if (nCC>0) then
      if (allocated(InitOut%CableCChanRqst)) deallocate(InitOut%CableCChanRqst)
      call AllocAry(InitOut%CableCChanRqst, maxCC, 'InitOut%CableCChanRqst', ErrStat2, ErrMsg2); if(Failed()) return;
      InitOut%CableCChanRqst = .FALSE.    ! Initialize to false
   endif

   ! --- Store mapping 
   nCC = 0
   do i = 1, size(p%ElemProps)
      if (p%ElemProps(i)%eType==idMemberCable) then
         idCProp= p%Elems(i,iMProp)
         if (NINT(Init%PropsC(idCProp, 5 ))>0) then
            nCC=nCC+1
            p%CtrlElem2Channel(nCC, 1) = i ! Element index (in p%Elems and p%ElemProps)
            p%CtrlElem2Channel(nCC, 2) = NINT(Init%PropsC(idCProp,5),IntKi) ! Control channel
            InitOut%CableCChanRqst(NINT(Init%PropsC(idCProp, 5 ),IntKi)) = .TRUE.
         endif
      endif
   enddo

   ! --- DeltaL Guess for inputs
   if (allocated(uInit%CableDeltaL)) deallocate(uInit%CableDeltaL)
   call AllocAry(uInit%CableDeltaL, nCC, 'uInit%CableDeltaL', ErrStat2, ErrMsg2); if(Failed()) return;
   do i = 1, nCC
       iElem    = p%CtrlElem2Channel(i,1)
       ! DeltaL 0 = - Le T0 / (EA + T0) = - Le eps0 / (1+eps0)
       !uInit%CableDeltaL(i) = - p%ElemProps(iElem)%Length * p%ElemProps(iElem)%T0  / (p%ElemProps(iElem)%YoungE*p%ElemProps(iElem)%Area   +  p%ElemProps(iElem)%T0)
       uInit%CableDeltaL(i) = 0.0_ReKi
   enddo

contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ControlCableMapping') 
        Failed =  ErrStat >= AbortErrLev
   end function Failed
end subroutine ControlCableMapping

!> Init for control Cable force
!! The change of cable forces due to the control is linear, so we just store a "unit" force vector
!! We will just scale this vector at each time step based on the control input (Tcontrol):
!!   Fcontrol =  (Tcontrol-T0) * Funit
!! We store it in "non-reduced" system since it will added to the external forces
SUBROUTINE ControlCableForceInit(p, m, ErrStat, ErrMsg)
   TYPE(SD_ParameterType),       INTENT(IN   ) :: p  !< Parameters 
   TYPE(SD_MiscVarType),         INTENT(INOUT) :: m  !< Misc
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER                  :: iCC, iElem
   REAL(FEKi)               :: FCe(12) ! Pretension force from cable element
   integer(IntKi), dimension(12) :: IDOF !  12 DOF indices in global unconstrained system
   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(ErrMsgLen)     :: ErrMsg2
   ErrMsg  = ""
   ErrStat = ErrID_None
   
   ! Allocating necessary arrays
   CALL AllocAry( m%FC_unit , p%nDOF, 'm%FC0' , ErrStat2, ErrMsg2); if(Failed()) return; ! Control cable force
   m%FC_unit = 0.0_ReKi

   ! loop over all elements, compute element matrices and assemble into global matrices
   DO iCC = 1, size(p%CtrlElem2Channel,1) 
      iElem = p%CtrlElem2Channel(iCC,1)
      CALL ElemF_Cable(1.0_ReKi, p%ElemProps(iElem)%DirCos, FCe) !< NOTE: using unitary load T0=1.0_ReKi
      ! --- Assembly in global unconstrained system
      IDOF = p%ElemsDOF(1:12, iElem)
      m%FC_unit( IDOF )    = m%FC_unit( IDOF ) + FCe(1:12) 
   ENDDO
   ! Transforming the vector into reduced, direct elimination system:
   !FC_red = matmul(transpose(p%T_red), FC)
   !if(allocated(FC)) deallocate(FC)

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ControlCableForceInit') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed
END SUBROUTINE ControlCableForceInit

!> Add soil stiffness and mass to global system matrices
!! Soil stiffness can come from two sources: 
!!   - "SSI" matrices (specified at reaction nodes)
!!   - "Soil" matrices (specified at Initalization)
SUBROUTINE InsertSoilMatrices(M, K, NodesDOF, Init, p, ErrStat, ErrMsg, Substract)
   real(FEKi), dimension(:,:),   intent(inout) :: M
   real(FEKi), dimension(:,:),   intent(inout) :: K
   type(IList),dimension(:),     intent(in   ) :: NodesDOF !< Map from Node Index to DOF lists
   type(SD_InitType),            intent(inout) :: Init ! TODO look for closest indices elsewhere
   type(SD_ParameterType),       intent(in   ) :: p
   integer(IntKi),               intent(  out) :: ErrStat     ! Error status of the operation
   character(*),                 intent(  out) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   logical, optional,            intent(in   ) :: SubStract   ! If present, and if true, substract instead of adding
   integer                    :: I, J, iiNode, nDOF
   integer                    :: iDOF, jDOF, iNode  !< DOF and node indices
   real(FEKi), dimension(6,6) :: K_soil, M_soil ! Auxiliary matrices for soil
   real(ReKi)                 :: Dist
   ErrMsg  = ""
   ErrStat = ErrID_None
   ! --- SSI matrices
   ! TODO consider doing the 21 -> 6x6 conversion while reading
   ! 6x6 matrix goes to one node of one element only
   do iiNode = 1, p%nNodes_C ! loop on constrained nodes
      iNode = p%Nodes_C(iiNode,1)
      nDOF=size(NodesDOF(iNode)%List)
      if (nDOF/=6) then
         ErrMsg='SSI soil matrix is to be inserted at SubDyn node '//Num2LStr(iNode)//', but this node has '//num2lstr(nDOF)//' DOFs';
         ErrStat=ErrID_Fatal; return
      endif
      call Array21_to_6by6(Init%SSIK(:,iiNode), K_soil)
      call Array21_to_6by6(Init%SSIM(:,iiNode), M_soil)
      if (present(Substract)) then
         if (Substract) then
            K_soil = - K_soil
            M_soil = - M_soil
         endif
      endif
      do I = 1, 6
         iDOF = NodesDOF(iNode)%List(I)   ! DOF index
         do J = 1, 6
            jDOF = NodesDOF(iNode)%List(J)   ! DOF index
            K(iDOF, jDOF) = K(iDOF, jDOF) + K_soil(I,J)
            M(iDOF, jDOF) = M(iDOF, jDOF) + M_soil(I,J)
         enddo
      enddo
   enddo
   ! --- "Soil" matrices
   if (allocated(Init%Soil_K)) then
      do iiNode = 1,size(Init%Soil_Points,2)
         ! --- Find closest node
         call FindClosestNodes(Init%Soil_Points(1:3,iiNode), Init%Nodes, iNode, Dist);
         if (Dist>0.1_ReKi) then
            ErrMsg='Closest SubDyn Node is node '//Num2LStr(iNode)//', which is more than 0.1m away from soildyn point '//num2lstr(iiNode);
            ErrStat=ErrID_Fatal; return
         endif
         Init%Soil_Nodes(iiNode) = iNode
         ! --- Insert/remove from matrices
         nDOF=size(NodesDOF(iNode)%List)
         if (nDOF/=6) then
            ErrMsg='Soil matrix is to be inserted at SubDyn node '//Num2LStr(iNode)//', but this node has '//num2lstr(nDOF)//' DOFs';
            ErrStat=ErrID_Fatal; return
         endif
         K_soil = Init%Soil_K(1:6,1:6,iiNode)
         if (present(Substract)) then
            if (Substract) then
               K_soil = - K_soil
            endif
         endif
         do I = 1, 6
            iDOF = NodesDOF(iNode)%List(I)   ! DOF index
            do J = 1, 6
               jDOF = NodesDOF(iNode)%List(J)   ! DOF index
               K(iDOF, jDOF) = K(iDOF, jDOF) + K_soil(I,J)
            enddo
         enddo
         if (.not.present(Substract)) then
            CALL WrScr('   Soil stiffness inserted at SubDyn node '//trim(Num2LStr(iNode)))
            print*,'    ',K_Soil(1,1:6)
            print*,'    ',K_Soil(2,1:6)
            print*,'    ',K_Soil(3,1:6)
            print*,'    ',K_Soil(4,1:6)
            print*,'    ',K_Soil(5,1:6)
            print*,'    ',K_Soil(6,1:6)
         endif
      enddo
   endif
contains
   !> Convert a flatten array of 21 values into a symmetric  6x6 matrix
   SUBROUTINE Array21_to_6by6(A21, M66)
      use NWTC_LAPACK, only: LAPACK_TPTTR 
      real(FEKi), dimension(21) , intent(in)  :: A21
      real(FEKi), dimension(6,6), intent(out) :: M66
      integer :: j
      M66 = 0.0_ReKi
      ! Reconstruct from sparse elements
      CALL LAPACK_TPTTR('U',6,A21,M66,6, ErrStat, ErrMsg)
      ! Ensuring symmetry
      do j=1,6
         M66(j,j) = M66(j,j)/2
      enddo  
      M66=M66+TRANSPOSE(M66) 
   END SUBROUTINE Array21_to_6by6
END SUBROUTINE InsertSoilMatrices

!------------------------------------------------------------------------------------------------------
!> Find closest node index to a point, returns distance as well
SUBROUTINE FindClosestNodes(Point, Nodes, iNode, Dist)
   real(ReKi), dimension(3),     intent(IN   ) :: Point  !< Point coordinates
   real(ReKi), dimension(:,:),   intent(IN   ) :: Nodes  !< List of nodes, Positions are in columns 2-4...
   integer(IntKi),               intent(  OUT) :: iNode  !< Index of closest node
   real(ReKi),                   intent(  OUT) :: Dist   !< Distance from Point to node iNode
   integer(IntKi) :: I
   real(ReKi) :: min_dist, loc_dist
   ! 
   min_dist=999999._ReKi
   iNode=-1
   do i = 1, size(Nodes,1)
      loc_dist = sqrt((Point(1) - Nodes(i,2))**2 + (Point(2) - Nodes(i,3))**2+ (Point(3) - Nodes(i,4))**2) 
      if (loc_dist<min_dist) then
         iNode=i
         min_dist = loc_dist
      endif
   enddo
   Dist=min_dist
END SUBROUTINE FindClosestNodes

!------------------------------------------------------------------------------------------------------
!> Build transformation matrix T, such that x= T.x~ where x~ is the reduced vector of DOF
!! Variables set by this routine
!! - p%NodesDOFred(iNode)=[list of DOF]: Created for each node, the list of DOF of this node in the 
!!         reduced system. 
!!         NOTE: follower nodes in rigid assembly have no DOFred (convention)
!! - p%nDOF_red: number of DOF in reduced system (<= nDOF)
!! - p%reduced: true if a reduction is needed, i.e. a T matrix is needed, and nDOF_red<nDOF
!!
!! Variables returned:
!! - T_red: retuction matrix such that x= T_red.x~ where x~ is the reduced vector of DOF
SUBROUTINE BuildTMatrix(Init, p, RA, RAm1, T_red, ErrStat, ErrMsg)
   use IntegerList, only: init_list, find, pop, destroy_list, len
   use IntegerList, only: print_list
   TYPE(SD_InitType),            INTENT(IN   ) :: Init
   TYPE(SD_ParameterType),target,INTENT(INOUT) :: p
   type(IList), dimension(:),    INTENT(IN   ) :: RA   !< RA(a) = [e1,..,en]  list of elements forming a rigid link assembly
   integer(IntKi), dimension(:), INTENT(IN   ) :: RAm1 !< RA^-1(e) = a , for a given element give the index of a rigid assembly
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   real(FEKi), dimension(:,:), allocatable :: T_red !< Transformation matrix for DOF elimination
   ! Local  
   real(ReKi), dimension(:,:), allocatable   :: Tc
   integer(IntKi), dimension(:), allocatable :: INodesID !< List of unique nodes involved in Elements
   integer(IntKi), dimension(:), allocatable :: IDOFOld !< 
   integer(IntKi), dimension(:), pointer :: IDOFNew !< 
   real(ReKi), dimension(6,6) :: I6       !< Identity matrix of size 6
   integer(IntKi) :: iPrev
   type(IList) :: IRA !< list of rigid assembly indices to process
   integer(IntKi) :: aID, ia ! assembly ID, and index in IRA
   integer(IntKi) :: iNode, iNodeSel, iNodeRemaining, iiNodeRemaining
   integer(IntKi) :: er !< Index of one rigid element belong to a rigid assembly
   integer(IntKi) :: JType
   integer(IntKi) :: I
   integer(IntKi) :: nc !< Number of DOF after constraints applied
   integer(IntKi) :: nj
   real(ReKi)  :: phat(3) !< Directional vector of the joint
   type(IList), dimension(:), allocatable :: RA_DOFred ! DOF indices for each rigid assembly, in reduced system
   INTEGER(IntKi)       :: ErrStat2
   CHARACTER(ErrMsgLen) :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! --- Misc inits
   nullify(IDOFNew)
   I6(1:6,1:6)=0; do i = 1,6 ; I6(i,i)=1_ReKi; enddo ! I6 =  eye(6)
   allocate(p%NodesDOFred(1:p%nNodes), stat=ErrStat2); if(Failed()) return; ! Indices of DOF for each joint, in reduced system
   allocate(RA_DOFred(1:size(RA)), stat=ErrStat2); if(Failed()) return; ! Indices of DOF for each rigid assmbly, in reduced system

   p%nDOF_red = nDOF_ConstraintReduced()
   p%reduced  = reductionNeeded()      ! True if reduction needed, allow for optimization if not needed

   if (DEV_VERSION) then
      print*,'nDOF constraint elim', p%nDOF_red , '/' , p%nDOF
   endif
   CALL AllocAry( T_red, p%nDOF, p%nDOF_red, 'p%T_red',  ErrStat2, ErrMsg2); if(Failed()) return; ! system stiffness matrix 
   T_red=0.0_FeKi
   call init_list(IRA, size(RA), 0, ErrStat2, ErrMsg2); if(Failed()) return;
   IRA%List(1:size(RA)) = (/(ia , ia = 1,size(RA))/)
   if (DEV_VERSION) then
      call print_list(IRA, 'List of RA indices')
   endif

   ! --- For each node:
   !  - create list of indices I      in the assembled vector of DOF
   !  - create list of indices Itilde in the reduced vector of DOF
   !  - increment iPrev by the number of DOF of Itilde
   iPrev =0 
   do iNode = 1, p%nNodes
      iNodeSel = iNode ! Unless changed by Rigid assembly, using this index
      if (allocated(Tc))      deallocate(Tc)
      if (allocated(IDOFOld)) deallocate(IDOFOld)
      JType = int(Init%Nodes(iNodeSel,iJointType))
      if(JType == idJointCantilever ) then
         if ( NodeHasRigidElem(iNodeSel, Init, p, er)) then ! return True and element index "er" if element is rigid
            ! --- The joint is involved in a rigid link assembly
            aID = RAm1(er) ! ID of rigid assembly for element er
            if (aID<0) then
               call Fatal('No rigid assembly attributed to node'//trim(Num2LStr(iNodeSel))//'. RAm1 wrong'); return
            endif
            ia  = find(IRA, aID, ErrStat2, ErrMsg2); if(Failed()) return ! We "pop" IRA, so index and ID are different
            if (DEV_VERSION) then
               print'(4X,A,I5,A,I5,A,I5)','Node',iNodeSel, ' is involved in RA:', aID, '. Current index in list of RA', ia
            endif
            if ( ia <= 0) then
               ! --- This rigid assembly has already been processed, simple triggers below
               ! OLD: The DOF list is taken from the stored RA DOF list
               ! call init_list(p%NodesDOFred(iNodeSel), RA_DOFred(aID)%List, ErrStat2, ErrMsg2)
               ! NEW: this node has no DOFs, so we set an empty list of DOFred for this node
               !call init_list(p%NodesDOFred(iNodeSel), 0, 0, ErrStat2, ErrMsg2)
               if (DEV_VERSION) then
                  print*,'   The RA',aID,', has already been processed!'! The following node has no reduced DOF'
                  !print*,'   but based on its RA, we can list its Itilde DOF:'
                  !print*,'   N',iNodeSel,'I ',p%NodesDOF(iNodeSel)%List(1:6)
                  !print*,'   N',iNodeSel,'It',RA_DOFred(aID)%List
               endif
               cycle ! We pass to the next joint, important so that:
               !     - we don't increase iPrev
               !     - we don't set Tc
               !     - p%NodesDOFred is not set (assuming it has already been done)
            else
               ! --- Proceeding the rigid assembly
               ! Returns TC and INodesID, do not change other variables
               call RAElimination( RA(aID)%List, Tc, INodesID, Init, p, ErrStat2, ErrMsg2); if(Failed()) return;
               aID = pop(IRA, ia, ErrStat2, ErrMsg2) ! this assembly has been processed, remove it from IRA list
               nj = size(INodesID) ! Number of nodes in this rigid assembly
               allocate(IDOFOld(1:6*nj))
               do I=1, nj
                  IDOFOld( (I-1)*6+1 : I*6 ) = p%NodesDOF(INodesID(I))%List(1:6)
               enddo

               ! Storing DOF list for this RA (Note: same as NodesDOFred below, only for debug)
               nc=size(Tc,2) ! Should be 6 
               call init_list(RA_DOFred(aID), (/ (iprev + i, i=1,nc) /), ErrStat2, ErrMsg2);

               ! --- Processing trigger for leader/follower Nodes
               iNodeSel = INodesID(1)  ! The first index returned is the leader of the assembly, we use this from now on
               do iiNodeRemaining=2,size(INodesID) ! start at 2 because 1 is always the leader
                  iNodeRemaining = INodesID(iiNodeRemaining)
                  ! OLD: The DOF list is taken from the stored RA DOF list
                  ! call init_list(p%NodesDOFred(iNode), RA_DOFred(aID)%List, ErrStat2, ErrMsg2)
                  ! NEW: this node has no DOFs, so we set an empty list of DOFred for this node
                  call init_list(p%NodesDOFred(iNodeRemaining), 0, 0, ErrStat2, ErrMsg2)
                  if (DEV_VERSION) then
                     print'(4X,A,I5,A,I5,I5)','Node',iNodeRemaining,' has no reduced DOF since its the follower of leader node ',INodesID(1),iNodeSel
                  endif
               enddo
            endif
         else
            ! --- Regular cantilever joint
            ! TODO/NOTE: We could apply fixed constraint/BC here, returning Tc as a 6xn matrix with n<6
            !            Extreme case would be Tc: 6*0, in which case NodesDOFred would be empty ([])
            allocate(Tc(1:6,1:6))
            allocate(IDOFOld(1:6))
            Tc=I6
            IDOFOld = p%NodesDOF(iNodeSel)%List(1:6)
         endif
      else
         ! --- Ball/Pin/Universal joint
         allocate(IDOFOld(1:len(p%NodesDOF(iNodeSel))))
         IDOFOld(:) = p%NodesDOF(iNodeSel)%List(:)
         phat = Init%Nodes(iNodeSel, iJointDir:iJointDir+2)
         ! Return Tc, do not change other variable
         call JointElimination(Init%NodesConnE(iNodeSel,:), JType, phat, p, Tc, ErrStat2, ErrMsg2); if(Failed()) return
      endif ! Cantilever or Special Joint
      nc=size(Tc,2) 
      call init_list(p%NodesDOFred(iNodeSel), nc, 0, ErrStat2, ErrMsg2)
      p%NodesDOFred(iNodeSel)%List(1:nc) = (/ (iprev + i, i=1,nc) /)
      IDOFNew => p%NodesDOFred(iNodeSel)%List(1:nc) ! alias to shorten notations
      if (DEV_VERSION) then
         ! KEEP ME, VERY USEFUL
         print*,'N',iNodeSel,'I ',IDOFOld
         print*,'N',iNodeSel,'It',IDOFNew
      endif
      T_red(IDOFOld, IDOFNew) = Tc
      iPrev = iPrev + nc
   enddo
   if (DEV_VERSION) then
      print'(A)','--- End of BuildTMatrix'
      print*,'   - T_red set'
      print*,'   - p%nDOF_red', p%nDOF_red
      print*,'   - p%reduced ', p%reduced
      print*,'   - p%NodesDOFred: (list of reduced DOF indices per node) '
      do iNode = 1, p%nNodes
         print*,'N',iNode, 'It', p%NodesDOFred(iNode)%List(:)
      enddo
   endif

   ! --- Safety checks
   if (len(IRA)>0) then 
      call Fatal('Not all rigid assemblies were processed'); return
   endif
   if (iPrev /= p%nDOF_red) then 
      call Fatal('Inconsistency in number of reduced DOF'); return
   endif
   call CleanUp_BuildTMatrix()
contains
   LOGICAL FUNCTION Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'BuildTMatrix') 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call CleanUp_BuildTMatrix()
   END FUNCTION Failed

   SUBROUTINE Fatal(ErrMsg_in)
      CHARACTER(len=*), intent(in) :: ErrMsg_in
      CALL SetErrStat(ErrID_Fatal, ErrMsg_in, ErrStat, ErrMsg, 'BuildTMatrix');
   END SUBROUTINE Fatal

   SUBROUTINE CleanUp_BuildTMatrix()
      nullify(IDOFNew)
      call destroy_list(IRA, ErrStat2, ErrMsg2)
      if (allocated(Tc)     ) deallocate(Tc)
      if (allocated(IDOFOld)) deallocate(IDOFOld)
      if (allocated(INodesID)) deallocate(INodesID)
      if (allocated(RA_DOFred)) deallocate(RA_DOFred)
   END SUBROUTINE CleanUp_BuildTMatrix

   !> Returns number of DOF after constraint reduction (via the matrix T)
   INTEGER(IntKi) FUNCTION nDOF_ConstraintReduced()
      integer(IntKi) :: iNode
      integer(IntKi) :: ia ! Index on rigid link assembly
      integer(IntKi) :: m  ! Number of elements connected to a joint
      integer(IntKi) :: NodeType
      nDOF_ConstraintReduced = 0

      ! Rigid assemblies contribution
      nDOF_ConstraintReduced = nDOF_ConstraintReduced + 6*size(RA)

      ! Contribution from all the other joints
      do iNode = 1, p%nNodes
         m = Init%NodesConnE(iNode,1) ! Col1: number of elements connected to this joint
         NodeType = Init%Nodes(iNode,iJointType)

         if    (NodeType == idJointPin ) then
            nDOF_ConstraintReduced = nDOF_ConstraintReduced + 5 + 1*m
            print'(4X,A,I5,A,I5)','Node',iNode, ' is a pin joint, number of members involved: ',m

         elseif(NodeType == idJointUniversal ) then
            nDOF_ConstraintReduced = nDOF_ConstraintReduced + 4 + 2*m
            print'(4X,A,I5,A,I5)','Node',iNode, ' is an universal joint, number of members involved: ',m

         elseif(NodeType == idJointBall ) then
            nDOF_ConstraintReduced = nDOF_ConstraintReduced + 3 + 3*m
            print'(4X,A,I5,A,I5)','Node',iNode, ' is a ball joint, number of members involved: ',m

         elseif(NodeType == idJointCantilever ) then
            if ( NodeHasRigidElem(iNode, Init, p, er)) then
               ! This joint is involved in a rigid link assembly, we skip it (accounted for above)
               print'(4X,A,I5,A,I5)','Node',iNode, ' is involved in a Rigid assembly'
            else
               ! That's a regular Cantilever joint
               nDOF_ConstraintReduced = nDOF_ConstraintReduced + 6
               !print*,'Node',iNode, 'is a regular cantilever'
            endif
         else
            ErrMsg='Wrong joint type'; ErrStat=ErrID_Fatal
         endif
      end do
   END FUNCTION nDOF_ConstraintReduced

   !> return true if reduction needed (i.e. special joints, special elements)
   logical FUNCTION reductionNeeded()
      integer(IntKi) :: i
      integer(IntKi) :: myType
      reductionNeeded=.false.
      ! Rigid or cable links
      do i =1,size(p%Elems,1)
         myType = p%Elems(i, iMType)
         if (any((/idMemberCable, idMemberRigid/)==myType)) then
            reductionNeeded=.true.
            return
         endif
      enddo
      ! Special joints
      do i = 1, p%nNodes
         myType = Init%Nodes(i,iJointType)
         if (any((/idJointPin, idJointUniversal, idJointBall/)==myType)) then
            reductionNeeded=.true.
            return
         endif
      enddo
   end FUNCTION reductionNeeded

END SUBROUTINE BuildTMatrix
!------------------------------------------------------------------------------------------------------
!> Assemble stiffness and mass matrix, and gravity force vector
SUBROUTINE DirectElimination(Init, p, ErrStat, ErrMsg)
   use NWTC_LAPACK, only: LAPACK_GEMM
   use IntegerList, only: len
   TYPE(SD_InitType),            INTENT(INOUT) :: Init
   TYPE(SD_ParameterType),target,INTENT(INOUT) :: p
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER(IntKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   ! Varaibles for rigid assembly
   type(IList), dimension(:), allocatable    :: RA       !< RA(a) = [e1,..,en]  list of elements forming a rigid link assembly
   integer(IntKi), dimension(:), allocatable :: RAm1 !< RA^-1(e) = a , for a given element give the index of a rigid assembly
   real(FEKi), dimension(:,:), allocatable :: MM, KK
   real(FEKi), dimension(:,:), allocatable :: Temp
   integer(IntKi) :: nDOF, iDOF, nDOFPerNode, iNode, iiDOF, i,j
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Setup list of rigid link assemblies (RA) and the inverse function RA^{-1}
   call RigidLinkAssemblies(Init, p, RA, RAm1, ErrStat2, ErrMsg2); if(Failed()) return
   call BuildTMatrix(Init, p, RA, RAm1, p%T_red, ErrStat2, ErrMsg2); if (Failed()) return
   if (allocated(RAm1)) deallocate(RAm1)
   if (allocated(RA  )) deallocate(RA  )

   ! --- DOF elimination for system matrices and RHS vector
   nDOF = p%nDOF_red
   if (p%reduced) then
      ! Temporary backup of M and K of full system (Flang compiler failed when move_alloc was used here, so arrays are allocated, moved, and deallocated manually)
      CALL AllocAry(KK, size(Init%K,1), size(Init%K,2), 'KK',  ErrStat2, ErrMsg2); if(Failed()) return; ! system stiffness matrix 
      CALL AllocAry(MM, size(Init%M,1), size(Init%M,2), 'MM',  ErrStat2, ErrMsg2); if(Failed()) return; ! system mass matrix 
      KK = Init%K
      MM = Init%M
      deallocate(Init%K, Init%M)
      !  Reallocating
      CALL AllocAry( Init%K,      nDOF, nDOF,       'Init%K'   ,  ErrStat2, ErrMsg2); if(Failed()) return; ! system stiffness matrix 
      CALL AllocAry( Init%M,      nDOF, nDOF,       'Init%M'   ,  ErrStat2, ErrMsg2); if(Failed()) return; ! system mass matrix 
      CALL AllocAry( Temp   ,size(MM,1), nDOF,      'Temp'     ,  ErrStat2, ErrMsg2); if(Failed()) return; 
      CALL AllocAry( p%T_red_T,nDOF   , size(MM,1), 'T_red_T' ,  ErrStat2, ErrMsg2); if(Failed()) return; 
      ! --- Elimination (stack expensive)
      !Init%M  = matmul(transpose(p%T_red), matmul(MM, p%T_red))
      !Init%K  = matmul(transpose(p%T_red), matmul(KK, p%T_red))
      !p%T_red_T = transpose(p%T_red)
      do i = 1, size(p%T_red,1)
         do j = 1, size(p%T_red,2)
            p%T_red_T(j,i) = p%T_red(i,j)
         enddo
      enddo
      !Temp    = matmul(MM, p%T_red)
      CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, MM     , p%T_red, 0.0_FeKi, Temp  , ErrStat2, ErrMsg2); if(Failed()) return
      !Init%M  = matmul(p%T_red_T, Temp)
      CALL LAPACK_gemm( 'T', 'N', 1.0_FeKi, p%T_red, Temp   , 0.0_FeKi, Init%M, ErrStat2, ErrMsg2); if(Failed()) return
      !Temp    = matmul(KK, p%T_red)
      CALL LAPACK_gemm( 'N', 'N', 1.0_FeKi, KK     , p%T_red, 0.0_FeKi, Temp  , ErrStat2, ErrMsg2); if(Failed()) return
      !Init%K  = matmul(p%T_red_T, Temp)
      CALL LAPACK_gemm( 'T', 'N', 1.0_FeKi, p%T_red, Temp   , 0.0_FeKi, Init%K, ErrStat2, ErrMsg2); if(Failed()) return
      if (allocated(Temp))    deallocate(Temp)
   endif
   !CALL AllocAry( Init%D,      nDOF, nDOF,  'Init%D'   ,  ErrStat2, ErrMsg2); if(Failed()) return; ! system damping matrix 
   !Init%D = 0 !< Used for additional damping 

   ! --- Creating a convenient Map from DOF to Nodes
   call AllocAry(p%DOFred2Nodes, p%nDOF_red, 3, 'DOFred2Nodes', ErrStat2, ErrMsg2); if(Failed()) return;
   p%DOFred2Nodes=-999
   do iNode=1,p%nNodes
      nDOFPerNode = len(p%NodesDOFred(iNode))
      do iiDOF = 1, nDOFPerNode
         iDOF = p%NodesDOFred(iNode)%List(iiDOF)
         p%DOFred2Nodes(iDOF,1) = iNode       ! First column is Node index
         p%DOFred2Nodes(iDOF,2) = nDOFPerNode ! Second column is number of DOF per node
         p%DOFred2Nodes(iDOF,3) = iiDOF       ! Third column is number of DOF per node
      enddo
   enddo

   call CleanUp_DirectElimination()

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'DirectElimination') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp_DirectElimination()
   END FUNCTION Failed
   SUBROUTINE CleanUp_DirectElimination()
      ! Cleaning up memory
      if (allocated(MM  )) deallocate(MM  )
      if (allocated(KK  )) deallocate(KK  )
      if (allocated(RA  )) deallocate(RA  )
      if (allocated(RAm1)) deallocate(RAm1)
      if (allocated(Temp)) deallocate(Temp)
   END SUBROUTINE CleanUp_DirectElimination
END SUBROUTINE DirectElimination

!------------------------------------------------------------------------------------------------------
!> Returns constraint matrix Tc for a rigid assembly (RA) formed by a set of elements. 
!!   x_c = Tc.x_c_tilde  
!! where x_c are all the DOF of the rigid assembly, and x_c_tilde are the 6 reduced DOF (leader DOF)
SUBROUTINE RAElimination(Elements, Tc, INodesID, Init, p, ErrStat, ErrMsg)
   use IntegerList, only: init_list, len, append, print_list, pop, destroy_list, get, unique, find
   integer(IntKi), dimension(:), INTENT(IN   ) :: Elements !< List of elements
   real(ReKi), dimension(:,:), allocatable     :: Tc
   integer(IntKi), dimension(:), allocatable   :: INodesID !< List of unique nodes involved in Elements
   TYPE(SD_InitType),            INTENT(IN   ) :: Init
   TYPE(SD_ParameterType),       INTENT(IN   ) :: p
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat  !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg   !< Error message if ErrStat /= ErrID_None
   ! Local variables
   type(IList)          :: LNodesID     !< List of nodes id involved in element
   type(IList)          :: LNodesInterf !< List of nodes id involved in interface
   integer(IntKi)       :: NodeID   !< NodeID
   integer(IntKi)       :: iTmp     !< Temporary index
   integer(IntKi)       :: iNodeID  !< Loop index on node ID list
   integer(IntKi)       :: iiMainNode !< Index of main node selected for rigid assembly within INodesID list
   integer(IntKi)       :: iMainNode !< Main node index
   integer(IntKi)       :: nNodes  !< Number of Nodes involved in RA
   integer(IntKi)       :: iFound  !< Loop index on node ID list
   integer(IntKi)       :: i       !< Loop index 
   real(ReKi)           :: TRigid(6,6) ! Transformation matrix such that xi = T.x1
   real(ReKi)           :: P1(3), Pi(3) ! Nodal points
   INTEGER(IntKi)       :: ErrStat2
   CHARACTER(ErrMsgLen) :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! --- List of nodes stored first in LINodes than moved to INodes
   LNodesID = NodesList(p, Elements)
   if (DEV_VERSION) then
      print*,'   --- RAElimination, Processing a rigid assembly'
      print*,'   Nodes involved in assembly (before any manip) ',LNodesID%List
   endif
   call unique(LNodesID, ErrStat2, ErrMsg2);
   if (DEV_VERSION) then
      print*,'   Nodes involved in assembly (selecting unique) ',LNodesID%List
   endif

   !--- Look for potential interface node
   call init_list(LNodesInterf, 0, 0, ErrStat2, ErrMsg2);
   do iNodeID = 1, len(LNodesID)
      NodeID = LNodesID%List(iNodeID)
      iFound =  FINDLOCI( p%Nodes_I(:,1), NodeID)
      if (iFound>0) then
         call append(LNodesInterf, NodeID, ErrStat2, ErrMsg2)
         ! This node is an interface node
         print'(4X,A,I5,A)','Node',NodeID, ' is an interface node, selecting it for the rigid assembly'
      endif
   enddo

   ! --- Decide which node will be the main node of the rigid assembly
   if      (len(LNodesInterf)==0) then
      iiMainNode = 1 ! By default we select the first node
   else if (len(LNodesInterf)==1) then
      ! Finding the index of the interface node
      iMainNode  = pop(LNodesInterf, ErrStat2, ErrMsg2)
      iiMainNode = find(LNodesID, iMainNode, ErrStat2, ErrMsg2);
   else
      ErrStat=ErrID_Fatal
      ErrMsg='Cannot have several interface nodes linked within a same rigid assembly'
      return
   endif
   call destroy_list(LNodesInterf, ErrStat2, ErrMsg2)
   if (DEV_VERSION) then
      print'(4X,A,I5)','We will select the node at position ',iiMainNode
      print*,'   in the following list of nodes:               ',LNodesID%List
   endif

   ! --- Extracting index array from list
   if (allocated(INodesID)) deallocate(INodesID)
   call move_alloc(LNodesID%List, INodesID)
   call destroy_list(LNodesID, ErrStat2, ErrMsg2)

   ! --- Order list of joints with main node first (swapping iMainNode with INodes(1))
   iTmp                 = INodesID(1)
   INodesID(1)          = INodesID(iiMainNode)
   INodesID(iiMainNode) = iTmp
   print*,'   Nodes involved in assembly:',INodesID

   ! --- Building Transformation matrix
   nNodes =size(INodesID)
   allocate(Tc(6*nNodes,6)) ! NOTE: do not deallocate, this is an ouput of this function
   Tc(:,:)=0
   ! I6 for first node since it's the "leader"
   do i = 1,6 ; Tc(i,i)=1_ReKi; enddo ! I6 =  eye(6)
   ! Rigid transformation matrix for the other nodes 
   P1 = Init%Nodes(INodesID(1), 2:4) ! reference node coordinates
   do i = 2, nNodes
      Pi = Init%Nodes(INodesID(i), 2:4) ! follower node coordinates
      call GetRigidTransformation(P1, Pi, TRigid, ErrStat2, ErrMsg2)
      Tc( ((i-1)*6)+1:6*i, 1:6) = TRigid(1:6,1:6)
      if (DEV_VERSION) then
         print'(4X,A,3(F6.1),A,3(F6.1))','Rigid transformation from ref point',P1,' to ',Pi
      endif
   enddo
END SUBROUTINE RAElimination
!------------------------------------------------------------------------------------------------------
!> Returns constraint matrix Tc for a joint involving several Elements
!!   x_c = Tc.x_c_tilde  
!! where
!    x_c       are all the DOF of the joint (3 translation + 3*m, m the number of elements) 
!    x_c_tilde are the nc reduced DOF 
SUBROUTINE JointElimination(Elements, JType, phat, p, Tc, ErrStat, ErrMsg)
   use IntegerList, only: init_list, len, append, print_list, pop, destroy_list, get
   integer(IntKi), dimension(:), INTENT(IN   ) :: Elements !< List of elements involved at a joint
   integer(IntKi),               INTENT(IN   ) :: JType !< Joint type
   real(ReKi),                   INTENT(IN   ) :: phat(3) !< Directional vector of the joint
   TYPE(SD_ParameterType),       INTENT(IN   ) :: p
   real(ReKi), dimension(:,:), allocatable     :: Tc  !< Transformation matrix from eliminated to full
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat  !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg   !< Error message if ErrStat /= ErrID_None
   ! Local variables
   !type(IList)          :: I !< List of indices for Nodes involved in interface
   integer(IntKi)       :: i, j, ie, ne       !< Loop index 
   integer(IntKi)       :: nDOFr     !< Number of reduced DOF
   integer(IntKi)       :: nDOFt     !< Number of total DOF *nreduced)
   real(ReKi)           :: e1(3), e2(3), e3(3) ! forming orthonormal basis with phat 
   integer(IntKi)       :: ErrStat2
   character(ErrMsgLen) :: ErrMsg2
   real(FEKi), dimension(:,:), allocatable :: Tc_rot !< Part of Tc just for rotational DOF
   real(FEKi), dimension(:,:), allocatable :: Tc_rot_m1 !< Inverse of Tc_rot
   real(ReKi) :: ColMean
   ErrStat = ErrID_None
   ErrMsg  = ""

   ne = Elements(1) ! TODO TODO
   nDOFt = 3 + 3*ne

   ! The elements already share the same translational DOF

   if    (JType == idJointPin ) then
      nDOFr = 5 + 1*ne
      allocate(Tc  (nDOFt, nDOFr)); 
      allocate(Tc_rot_m1(nDOFr-3, nDOFt-3)); 
      Tc(:,:)=0
      Tc_rot_m1(:,:)=0

      ! Normalizing 
      e3= phat/sqrt(phat(1)**2 + phat(2)**2 + phat(3)**2)
      call GetOrthVectors(e3, e1, e2, ErrStat2, ErrMsg2);
      ! Forming Tcm1, inverse of Tc
      do ie=1,ne
         Tc_rot_m1(1   , (ie-1)*3+1:ie*3 ) = e1(1:3)/ne
         Tc_rot_m1(2   , (ie-1)*3+1:ie*3 ) = e2(1:3)/ne
         Tc_rot_m1(ie+2, (ie-1)*3+1:ie*3 ) = e3(1:3)
      enddo
      ! Pseudo inverse:
      call PseudoInverse(Tc_rot_m1, Tc_rot, ErrStat2, ErrMsg2)
      ! --- Forming Tc
      do i = 1,3    ; Tc(i,i)=1_ReKi; enddo !  I3 for translational DOF
      Tc(4:nDOFt,4:nDOFr)=Tc_rot(1:nDOFt-3, 1:nDOFr-3)
      do i = 1,size(Tc,1); do ie = 1,size(Tc,2)
         if (abs(Tc(i,ie))<1e-13) then
            Tc(i,ie)=0.0_ReKi
         endif; enddo;
      enddo;
      deallocate(Tc_rot)
      deallocate(Tc_rot_m1)

   elseif(JType == idJointUniversal ) then
      if (ne/=2) then
         ErrMsg='JointElimination: universal joints should only connect two elements.'; ErrStat=ErrID_Fatal
         return
      endif
      nDOFr = 4 + 2*ne
      allocate(Tc(nDOFt, nDOFr)); 
      allocate(Tc_rot_m1(nDOFr-3, nDOFt-3)); 
      Tc(:,:)=0
      Tc_rot_m1(:,:)=0 ! Important init
      ! Forming the inverse of Tc_rot
      Tc_rot_m1(1,1:3) = p%ElemProps(Elements(1))%DirCos(:,3)/2._ReKi
      Tc_rot_m1(1,4:6) = p%ElemProps(Elements(2))%DirCos(:,3)/2._ReKi
      Tc_rot_m1(2,1:3) = p%ElemProps(Elements(1))%DirCos(:,1)
      Tc_rot_m1(3,1:3) = p%ElemProps(Elements(1))%DirCos(:,2)
      Tc_rot_m1(4,4:6) = p%ElemProps(Elements(2))%DirCos(:,1)
      Tc_rot_m1(5,4:6) = p%ElemProps(Elements(2))%DirCos(:,2)
      ! Pseudo inverse
      call PseudoInverse(Tc_rot_m1, Tc_rot, ErrStat2, ErrMsg2)
      ! --- Forming Tc
      do i = 1,3    ; Tc(i,i)=1_ReKi; enddo !  I3 for translational DOF
      Tc(4:nDOFt,4:nDOFr)=Tc_rot(1:nDOFt-3, 1:nDOFr-3)
      deallocate(Tc_rot)
      deallocate(Tc_rot_m1)

   elseif(JType == idJointBall      ) then
      nDOFr = 3 + 3*ne
      allocate(Tc(nDOFt, nDOFr)); 
      Tc(:,:)=0
      do i = 1,3    ; Tc(i,i)=1_ReKi; enddo !  I3 for translational DOF
      do i = 3,nDOFr; Tc(i,i)=1_ReKi; enddo ! Identity for other DOF as well

   else
      ErrMsg='JointElimination: Wrong joint type'; ErrStat=ErrID_Fatal
   endif
   !do i=1,nDOFt
   !   print*,'Tc',Tc(i,:)
   !enddo
   ! --- Safety check
   do j =1, size(Tc,2)
      ColMean=0; do i=1,size(Tc,1) ; ColMean = ColMean + abs(Tc(i,j)); enddo
      ColMean = ColMean/size(Tc,1)
      if (ColMean<1e-6) then
         ErrMsg='JointElimination: a reduced degree of freedom has a singular mapping.'; ErrStat=ErrID_Fatal
         return
      endif
   enddo

END SUBROUTINE JointElimination

!------------------------------------------------------------------------------------------------------
!> Setup a list of rigid link assemblies (RA)
!! Variables created by this routine:
!! - RA(ia)= [e1,..,en]  list of elements forming each rigid link assembly "ia".
!!                       Needed for BuildTMatrix
!! - RAm1(e)=(RA^-1(e)= a) : for a given element give the index of a rigid assembly. 
!!                       Needed for BuildTMatrix
!!
SUBROUTINE RigidLinkAssemblies(Init, p, RA, RAm1, ErrStat, ErrMsg)
   use IntegerList, only: init_list, len, append, print_list, pop, destroy_list, get
   TYPE(SD_InitType),            INTENT(INOUT) :: Init
   TYPE(SD_ParameterType),       INTENT(INOUT) :: p
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   type(IList), dimension(:), allocatable    :: RA   !< RA(a) = [e1,..,en]  list of elements forming a rigid link assembly
   integer(IntKi), dimension(:), allocatable :: RAm1 !< RA^-1(e) = a , for a given element give the index of a rigid assembly
   ! Local variables
   type(IList)                               :: Er    !< List of rigid elements
   type(IList)                               :: Ea    !< List of elements in a rigid assembly
   integer(IntKi)                            :: nRA  !< Number of rigid assemblies
   integer(IntKi)                            :: ie  !< Index on elements
   integer(IntKi)                            :: ia  !< Index on assemblies
   integer(IntKi)                            :: e0  !< Index of an element
   INTEGER(IntKi)       :: ErrStat2
   CHARACTER(ErrMsgLen) :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""
   allocate(RAm1(1:Init%NElem)) ! NOTE: do not deallocate, this is an "output" of this function
   RAm1(1:Init%NElem) = -1

   ! --- Establish a list of rigid link elements
   Er = RigidLinkElements(Init, p, ErrStat2, ErrMsg2)
   nRA=0
   do while (len(Er)>0)
      nRA=nRA+1
      ! Creating List Ea of elements of a given assembly
      call init_list(Ea, 0, 0, ErrStat2, ErrMsg2);
      e0 = pop(Er, ErrStat2, ErrMsg2);
      call append(Ea, e0, ErrStat2, ErrMsg2);
      call AddNeighbors(e0, Er, Ea)
      if (DEV_VERSION) then
         call print_list(Ea,'Rigid assembly (loop 1) element list')
      endif
      do ie = 1, len(Ea)
         e0 = get(Ea, ie, ErrStat2, ErrMsg2)
         RAm1(e0) = nRA ! Index of rigid assembly that this element belongs to
      enddo
      call destroy_list(Ea, ErrStat2, ErrMsg2)
   enddo
   call destroy_list(Er, ErrStat2, ErrMsg2)

   ! --- Creating RA, array of lists of assembly elements.
   ! Note: exactly the same as all the Ea created above, but we didn't know the total number of RA
   allocate(RA(1:nRA)) ! NOTE: do not deallocate, this is an "output" of this function
   do ia = 1, nRA
      call init_list(RA(ia), 0, 0, ErrStat2, ErrMsg2)
   enddo
   do ie = 1, Init%NElem
      ia = RAm1(ie) ! Index of the assembly the element belongs to: RA^{-1}(ie) = ia
      if (ia>0) then
         call append(RA(ia), ie, ErrStat2, ErrMsg2)
      endif
   enddo
   if (DEV_VERSION) then
      do ia = 1, nRA
         call print_list(RA(ia),'Rigid assembly (loop 2) element list')
      enddo
   endif
CONTAINS
   !> The neighbor-elements of element e0 (that are found within the list Er) are added to the list Ea  
   RECURSIVE SUBROUTINE AddNeighbors(e0, Er, Ea) 
      integer(IntKi), intent(in) :: e0  !< Index of an element
      type(IList), intent(inout) :: Er  !< List of rigid elements
      type(IList), intent(inout) :: Ea  !< List of elements in a rigid assembly
      type(IList)     :: En             !< List of neighbors of e0
      integer (IntKi) :: ik
      integer (IntKi) :: ek, ek2
      integer (IntKi) :: iWhichNode_e0, iWhichNode_ek
      call init_list(En, 0, 0, ErrStat2, ErrMsg2)
      ! Loop through all elements, setup list of e0-neighbors, add them to Ea, remove them from Er
      ik=0
      do while (ik< len(Er))
         ik=ik+1
         ek = Er%List(ik)
         if (ElementsConnected(p, e0, ek, iWhichNode_e0, iWhichNode_ek)) then
            if (DEV_VERSION) then
               print*,'Element ',ek,'is connected to element',e0,'via its node',iWhichNode_ek
            endif
            ! Remove element from Er (a rigid element can belong to only one assembly)
            ek2 =  pop(Er, ik,  ErrStat2, ErrMsg2) ! same as ek before
            ik=ik-1
            if (ek/=ek2) then
               print*,'Problem in popping',ek,ek2
               STOP
            endif
            call append(En, ek, ErrStat2, ErrMsg2)
            call append(Ea, ek, ErrStat2, ErrMsg2)
         endif
      enddo
      ! Loop through neighbors and recursively add neighbors of neighbors
      do ik = 1, len(En)
         ek = En%List(ik)
         call AddNeighbors(ek, Er, Ea)
      enddo
      call destroy_list(En, ErrStat2, ErrMsg2)
   END SUBROUTINE AddNeighbors

END SUBROUTINE RigidLinkAssemblies


!------------------------------------------------------------------------------------------------------
!> Add stiffness and damping to some joints
!! NOTE: damping was removed around 13/07/2020
SUBROUTINE InsertJointStiffDamp(p, Init, ErrStat, ErrMsg)
   TYPE(SD_ParameterType),target,INTENT(IN   ) :: p
   TYPE(SD_InitType),            INTENT(INOUT) :: Init
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi) :: iNode, JType, iStart, i
   integer(IntKi) :: nFreeRot ! Number of free rot DOF
   integer(IntKi) :: nMembers ! Number of members attached to this node
   integer(IntKi) :: nSpace   ! Number of spaces between diagonal "bands" (0:pin, 1:univ, 2:ball)
   real(ReKi) :: StifAdd
   real(ReKi), dimension(:,:), allocatable :: K_Add ! Stiffness matrix added to global system
   integer(IntKi), dimension(:), pointer :: Ifreerot
   ErrStat = ErrID_None
   ErrMsg  = ""
   do iNode = 1, p%nNodes
      JType   = int(Init%Nodes(iNode,iJointType))
      StifAdd = Init%Nodes(iNode, iJointStiff)
      if(JType == idJointCantilever ) then
         ! Cantilever joints should not have damping or stiffness
         if(StifAdd>0) then 
            ErrMsg='InsertJointStiffDamp: Additional stiffness should be 0 for cantilever joints. Index of problematic node: '//trim(Num2LStr(iNode)); ErrStat=ErrID_Fatal;
            return
         endif
      else
         ! Ball/Univ/Pin joints have damping/stiffness inserted at indices of "free rotation"
         nMembers = Init%NodesConnE(iNode,1) ! Col1: number of elements connected to this joint
         if      ( JType == idJointBall      ) then; iStart=4; nSpace=2;
         else if ( JType == idJointUniversal ) then; iStart=5; nSpace=1;
         else if ( JType == idJointPin       ) then; iStart=6; nSpace=0;
         endif
         Ifreerot=>p%NodesDOFred(iNode)%List(iStart:)
         nFreeRot = size(Ifreerot)
         ! Creating matrices of 0, and -K and nK on diagonals
         allocate(K_Add(1:nFreeRot,1:nFreeRot)); 
         call ChessBoard(K_Add, -StifAdd, 0._ReKi, nSpace=nSpace, diagVal=(nMembers-1)*StifAdd)
         ! Ball/Pin/Universal joints
         if(StifAdd>0) then 
            !print*,'Stiffness Add, Node:',iNode,'DOF:', Ifreerot
            !do i=1,nFreeRot
            !   print*,'K Add',K_Add(i,:)
            !enddo
            Init%K(Ifreerot,Ifreerot) = Init%K(Ifreerot,Ifreerot) + K_Add
         endif
         if(allocated(K_Add)) deallocate(K_Add)
      endif
   enddo
END SUBROUTINE InsertJointStiffDamp

!> Returns true if the substructure can be considered "fixed bottom"
LOGICAL FUNCTION isFixedBottom(Init, p)
   TYPE(SD_InitType),  INTENT(IN   ) :: Init
   TYPE(SD_ParameterType),INTENT(IN   ) :: p
   isFixedBottom=.not.isFloating(Init,p)
   !INTEGER(IntKi) :: i, nFixed
   !nFixed=0
   !do i =1,size(p%Nodes_C,1)
   !   if (ALL(p%Nodes_C(I,2:7)==idBC_Fixed)) then
   !      nFixed=nFixed+1
   !   elseif (Init%SSIfile(I)/='') then
   !      nFixed=nFixed+1
   !   endif
   !enddo
   !bFixed = nFixed >=1
END FUNCTION isFixedBottom

!> True if a structure is floating, no fixed BC at the bottom
logical function isFloating(Init, p) 
   type(SD_InitType),     intent(in   ):: Init
   type(SD_ParameterType),intent(in   ) :: p
   integer(IntKi) :: i
   !isFloating=size(p%Nodes_C)>0
   isFloating=.True.
   do i =1,size(p%Nodes_C,1)
      if ((all(p%Nodes_C(I,2:7)==idBC_Internal)) .and. (Init%SSIfile(i)=='')) then
         continue
      else
         isFloating=.False.
         return
      endif
   enddo
end function isFloating

SUBROUTINE ElemM(ep, Me)
   TYPE(ElemPropType), INTENT(IN) :: eP        !< Element Property
   REAL(FEKi), INTENT(OUT)        :: Me(12, 12)
   REAL(FEKi) :: L0, Eps0
   if (ep%eType==idMemberBeamCirc) then
      !Calculate Ke, Me to be used for output
      CALL ElemM_Beam(eP%Area, eP%Length, eP%Ixx, eP%Iyy, eP%Jzz,  eP%rho, eP%DirCos, Me)

   else if (ep%eType==idMemberCable) then
      Eps0 = ep%T0/(ep%YoungE*ep%Area)
      L0   = ep%Length/(1+Eps0)  ! "rest length" for which pretension would be 0
      CALL ElemM_Cable(ep%Area, L0, ep%rho, ep%DirCos, Me)

   else if (ep%eType==idMemberRigid) then
      if ( EqualRealNos(eP%rho, 0.0_ReKi) ) then
         Me=0.0_FEKi
      else
         CALL ElemM_Cable(ep%Area, real(ep%Length,FEKi), ep%rho, ep%DirCos, Me)
         !CALL ElemM_(A, L, rho, DirCos, Me)
      endif
   endif
END SUBROUTINE ElemM

SUBROUTINE ElemK(ep, Ke)
   TYPE(ElemPropType), INTENT(IN) :: eP        !< Element Property
   REAL(FEKi), INTENT(OUT)        :: Ke(12, 12)

   if (ep%eType==idMemberBeamCirc) then
      CALL ElemK_Beam( eP%Area, eP%Length, eP%Ixx, eP%Iyy, eP%Jzz, eP%Shear, eP%Kappa_x, eP%Kappa_y, eP%YoungE, eP%ShearG, eP%DirCos, Ke)

   else if (ep%eType==idMemberCable) then
      CALL ElemK_Cable(ep%Area, ep%Length, ep%YoungE, ep%T0, eP%DirCos, Ke)

   else if (ep%eType==idMemberRigid) then
      Ke = 0.0_FEKi
   endif
END SUBROUTINE ElemK

SUBROUTINE ElemF(ep, gravity, Fg, Fo)
   TYPE(ElemPropType), INTENT(IN) :: eP        !< Element Property
   REAL(ReKi), INTENT(IN)     :: gravity       !< acceleration of gravity
   REAL(FEKi), INTENT(OUT)    :: Fg(12)
   REAL(FEKi), INTENT(OUT)    :: Fo(12)
   if (ep%eType==idMemberBeamCirc) then
      Fo(1:12)=0.0_FEKi
   else if (ep%eType==idMemberCable) then
      CALL ElemF_Cable(ep%T0, ep%DirCos, Fo)
   else if (ep%eType==idMemberRigid) then
      Fo(1:12)=0.0_FEKi
   endif
   CALL ElemG( eP%Area, eP%Length, eP%rho, eP%DirCos, Fg, gravity )
END SUBROUTINE ElemF

!> Return skew symmetric matrix
SUBROUTINE skew(x,M33)
   real(ReKi), intent(in   ) :: x(3)
   real(ReKi), intent(  out) :: M33(3,3)
   M33(1 , :)=(/0.0_ReKi , -x(3)        , x (2)   /)
   M33(2 , :)=(/  x(3 )  , 0.0_ReKi     , -x(1)   /)
   M33(3 , :)=(/ -x(2 )  , x(1)        , 0.0_ReKi /)
END SUBROUTINE

!>Transform inertia matrix with respect to point P to the inertia matrix with respect to the COG
!!NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system.
SUBROUTINE translateInertiaMatrixToCOG(I_P, Mass, r_PG, I_G)
   real(ReKi), intent(in   ) :: I_P(3,3) !< Inertia matrix 3x3 with respect to point P
   real(ReKi), intent(in   ) :: Mass     !< Mass of the body
   real(ReKi), intent(in   ) :: r_PG(3)  !< vector from P to COG 
   real(ReKi), intent(  out) :: I_G(3,3) !< Inertia matrix (3x3) with respect to COG
   real(ReKi) :: S1(3,3) 
   call skew(r_PG, S1) 
   I_G = I_P + Mass * MATMUL(S1, S1)
END SUBROUTINE

!>Transform mass matrix with respect to point P to the mass matrix with respect to the COG
SUBROUTINE translateMassMatrixToCOG(MM, MM_G)
   real(ReKi), intent(in   ) :: MM(6,6)   !< Mass matrix (6x6) with respect to point P
   real(ReKi), intent(  out) :: MM_G(6,6) !< Mass matrix with respect to COG
   real(ReKi) :: m        ! Mass of the body
   real(ReKi) :: r_PG(3)  ! Vector from point P to G
   real(ReKi) :: J_P(3,3),  J_G(3,3) 
   ! Distance from refpoint to COG
   call rigidBodyMassMatrixCOG(MM, r_PG)
   ! Inertia at ref point
   J_P = MM(4:6,4:6)
   ! Inertia at COG
   call translateInertiaMatrixToCOG(J_P, MM(1,1), r_PG, J_G) 
   ! Rigid body mass matrix at COG
   call rigidBodyMassMatrix(MM(1,1), J_G(1,1), J_G(2,2), J_G(3,3), J_G(1,2), J_G(1,3), J_G(2,3), 0.0_ReKi, 0.0_ReKi, 0.0_ReKi, MM_G)
END SUBROUTINE

!>Transform mass matrix with respect to point P1 to the mass matrix with respect to point P2
SUBROUTINE translateMassMatrixToP(MM1, r_P1P2, MM2)
   real(ReKi), intent(in   ) :: MM1(6,6) !< Mass matrix (6x6) with respect to point P1
   real(ReKi), intent(in   ) :: r_P1P2(3)!< vector from P1 to P2
   real(ReKi), intent(  out) :: MM2(6,6) !< Mass matrix with respect to point P2
   real(ReKi) :: MM_G(6,6) !< Mass matrix with respect to COG
   real(ReKi) :: m        ! Mass of the body
   real(ReKi) :: r_P1G(3), r_P2G(3)  ! vector from P to COG 
   real(ReKi) :: J_G(3,3) 
   ! Rigid body mass matrix at COG to get inertia at COG
   call translateMassMatrixToCOG(MM1, MM_G)
   J_G = MM_G(4:6,4:6)
   ! Distance from refpoint to COG
   call rigidBodyMassMatrixCOG(MM1, r_P1G)
   r_P2G=-r_P1P2+r_P1G
   ! Rigid body mass matrix at Point P2
   call rigidBodyMassMatrix(MM1(1,1), J_G(1,1), J_G(2,2), J_G(3,3), J_G(1,2), J_G(1,3), J_G(2,3), r_P2G(1), r_P2G(2), r_P2G(3), MM2)
END SUBROUTINE

!> Return Center of gravity location from a 6x6 mass matrix
SUBROUTINE rigidBodyMassMatrixCOG(MM, r_PG)
   real(ReKi), intent(in   ) :: MM(6,6) !< Mass matrix (6x6) with respect to point P
   real(ReKi), intent(  out) :: r_PG(3) !< vector from P to G (center of mass)
   r_PG = (/ 0.5_ReKi*( MM(2,6)-MM(3,5)), & ! Using average of Coeffs
             0.5_ReKi*(-MM(1,6)+MM(3,4)), &
             0.5_ReKi*( MM(1,5)-MM(2,4)) /)
   r_PG = r_PG/MM(1,1)
END SUBROUTINE

!> Rigid body mass matrix (6x6) at a given reference point P
SUBROUTINE rigidBodyMassMatrix(m, Jxx, Jyy, Jzz, Jxy, Jxz, Jyz, x, y, z, M66)
   real(ReKi), intent(in   ) :: m             !< Mass of body
   real(ReKi), intent(in   ) :: Jxx, Jyy, Jzz, Jxy, Jxz, Jyz !< Inertia of body at COG
   real(ReKi), intent(in   ) :: x, y, z       !< x,y,z position of center of gravity (COG) with respect to the reference point
   real(ReKi), intent(  out) :: M66(6,6)      !< Mass matrix (6x6) with respect to point P
   M66(1 , :)=(/ m       , 0._ReKi , 0._ReKi , 0._ReKi             ,  z*m                , -y*m                 /)
   M66(2 , :)=(/ 0._ReKi , m       , 0._ReKi , -z*m                , 0._ReKi             ,  x*m                 /)
   M66(3 , :)=(/ 0._ReKi , 0._ReKi , m       ,  y*m                , -x*m                , 0._ReKi              /)
   M66(4 , :)=(/ 0._ReKi , -z*m    ,  y*m    , Jxx + m*(y**2+z**2) , Jxy - m*x*y         , Jxz  - m*x*z         /)
   M66(5 , :)=(/  z*m    , 0._ReKi , -x*m    , Jxy - m*x*y         , Jyy + m*(x**2+z**2) , Jyz  - m*y*z         /)
   M66(6 , :)=(/ -y*m    , x*m     , 0._ReKi , Jxz - m*x*z         , Jyz - m*y*z         , Jzz  + m*(x**2+y**2) /)
END SUBROUTINE

END MODULE SD_FEM
