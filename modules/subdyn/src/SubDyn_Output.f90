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
   INTEGER(IntKi)       :: ErrStat2      ! Error status of the operation
   CHARACTER(ErrMsgLen) :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                                   :: I,J,K,K2,L,NconEls   !Counters
   INTEGER(IntKi)                                   :: Junk  !Temporary Holders
   INTEGER(IntKi)                                   :: iMember  ! Member index (not member ID)
   INTEGER(IntKi), Dimension(2)                     :: M   !counter for two nodes at a time
   type(MeshAuxDataType), pointer :: pLst !< Alias to shorten notation and highlight code similarities
   REAL(ReKi) :: FCe(12) ! Pretension force from cable element
   ErrStat = 0      
   ErrMsg=""

   p%OutAllDims=12*p%NMembers*2    !size of AllOut Member Joint forces

   ! Check that the variables in OutList are valid      
   CALL SDOut_ChkOutLst( Init%SSOutList, p,  ErrStat2, ErrMsg2 ); if(Failed()) return

   IF ( ALLOCATED( p%OutParam ) .AND. p%NumOuts > 0 ) THEN           ! Output has been requested           
   ! Allocate SDWrOuput which is used to store a time step's worth of output channels, prior to writing to a file.
   CALL AllocAry(misc%SDWrOutput, p%NumOuts +p%OutAllInt*p%OutAllDims, 'SDWrOutupt', ErrStat2, ErrMsg2) ; if(Failed()) return
   misc%SDWrOutput  = 0.0_ReKi
   misc%LastOutTime = 0.0_DbKi
   misc%Decimat     = 0
   
   !Allocate WriteOuput  
   CALL AllocAry(y%WriteOutput, p%NumOuts +p%OutAllInt*p%OutAllDims, 'WriteOutput', ErrStat2, ErrMsg2); if(Failed()) return
   y%WriteOutput = 0

  
  DO I=1,p%NMOutputs
   pLst => p%MOutLst(I)
   CALL AllocAry(pLst%NodeIDs,    pLst%NoutCnt           , 'MOutLst(I)%NodeIDs', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(pLst%ElmIDs,     pLst%NoutCnt, p%NAvgEls, 'MOutLst(I)%ElmIDs' , ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(pLst%ElmNds,     pLst%NoutCnt, p%NAvgEls, 'MOutLst(I)%ElmNds' , ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(pLst%Me, 12, 12, pLst%NoutCnt, p%NAvgEls, 'MOutLst(I)%Me'     , ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(pLst%Ke, 12, 12, pLst%NoutCnt, p%NAvgEls, 'MOutLst(I)%Ke'     , ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(pLst%Fg,     12, pLst%NoutCnt, p%NAvgEls, 'MOutLst(I)%Fg'     , ErrStat2, ErrMsg2); if(Failed()) return

   ! NOTE: MemberNodes >2 if nDiv>1
   iMember = FINDLOCI(Init%Members(:,1), pLst%MemberID) ! Reindexing from MemberID to 1:nMembers
   pLst%NodeIDs=Init%MemberNodes(iMember,pLst%NodeCnt)  ! We are storing the actual node numbers corresponding to what the user ordinal number is requesting
   pLst%ElmIDs=0  !Initialize to 0
   pLst%ElmNds=0  !Initialize to 0

   DO J=1,pLst%NoutCnt !Iterate on requested nodes for that member
      !I need to get at most 2 elements that belong to the same MOutLst(I)%MemberID
      !make use of MemberNodes and NodesConnE
      NconEls=Init%NodesConnE(pLst%NodeIDs(J),1)!Number of elements connecting to the j-th node

      K2=0    !Initialize counter
      DO K=1, NconEls 
         L=Init%NodesConnE(pLst%NodeIDs(J),k+1)  !k-th Element Number 
         M     = p%Elems(L,2:3) !1st and 2nd node of the k-th element
         !Select only the other node, not the one where elements connect to
         IF (M(1) .EQ. pLst%NodeIDs(J)) then
            Junk=M(2)
         else
            Junk=M(1)
         endif
         IF (ANY(Init%MemberNodes(iMember,:) .EQ. Junk)) THEN  !This means we are in the selected member
            IF (K2 .EQ. 2) EXIT
            K2=K2+1
            pLst%ElmIDs(J,K2)=L        !This array has for each node requested NODEID(J), for each memberMOutLst(I)%MemberID, the 2 elements to average from, it may have 1 if one of the numbers is 0 
            IF (M(2) .EQ. pLst%NodeIDs(J) )then 
               pLst%ElmNds(J,K2)=2 !store whether first or second node of element  
            else
               pLst%ElmNds(J,K2)=1 !store whether first or second node of element  
            endif
            ! --- Element Me, Ke, Fg, Fce
            CALL ElemM(p%ElemProps(L),         pLst%Me(:,:,J,K2))
            CALL ElemK(p%ElemProps(L),         pLst%Ke(:,:,J,K2))
            CALL ElemF(p%ElemProps(L), Init%g, pLst%Fg(:,J,K2), FCe)
            pLst%Fg(:,J,K2) = pLst%Fg(:,J,K2) + FCe(1:12) ! Adding cable element force 
         END IF    
      ENDDO  ! K, NconEls
     ENDDO !J, Noutcnt
   ENDDO  !I, NMOutputs

   END IF   ! there are any requested outputs   
 
   IF (p%OutAll) THEN  !I need to store all member end forces and moments 
     
    ALLOCATE ( p%MOutLst2(p%NMembers), STAT = ErrStat2 )     !this list contains different arrays for each of its elements
    ErrMsg2 = 'Error allocating p%MOutLst2 array in SDOut_Init'
     
    DO I=1,p%NMembers
      pLst => p%MOutLst2(I)
      CALL AllocAry(pLst%NodeIDs, Init%Ndiv+1, 'MOutLst2(I)%NodeIDs', ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%ElmIDs,     2, 1, 'MOutLst2(I)%ElmIDs'     , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%ElmNds,     2, 1, 'MOutLst2(I)%ElmNds'     , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%Me, 12, 12, 2, 1, 'MOutLst2(I)%Me'         , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%Ke, 12, 12, 2, 1, 'MOutLst2(I)%Ke'         , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(pLst%Fg,     12, 2, 1, 'MOutLst2(I)%Fg'         , ErrStat2, ErrMsg2); if(Failed()) return
      pLst%MemberID=Init%Members(I,1)
      pLst%NodeIDs=Init%MemberNodes(I,1:Init%Ndiv+1)  !We are storing  the actual node numbers in the member
      !Now I need to find out which elements are attached to those nodes and still belong to the member I
      !ElmIDs could contain the same element twice if Ndiv=1
      pLst%ElmIDs=0  !Initialize to 0
      DO J=1,Init%Ndiv+1,Init%Ndiv !Iterate on requested nodes for that member (first and last)
          !I need to get at most 2 elements that belong to the same I Member
          !make use of MemberNodes and NodesConnE
          NconEls=Init%NodesConnE(pLst%NodeIDs(J),1) !Number of elements connecting to the 1st or last node of the member
          K2= J/(Init%Ndiv+1)+1  !store this quantity used later, basically 1 or 2 depending on J
          DO K=1, NconEls 
             L=Init%NodesConnE(pLst%NodeIDs(J),k+1)  !k-th Element Number in the set of elements attached to the selected node 
             M=p%Elems(L,2:3) !1st and 2nd node of the k-th element
             !Select only the other node, not the one where elements connect to
             Junk=M(1)
             IF (M(1) .EQ. pLst%NodeIDs(J)) Junk=M(2)
             IF (ANY(Init%MemberNodes(I,:) .EQ. Junk)) THEN  !This means we are in the selected member
                  pLst%ElmIDs(K2,1)=L     !This array has for each node requested NODEID(J), for each member I, the element to get results for 
                  pLst%ElmNds(K2,1)=1                        !store whether first or second node of element  
                  IF (M(2) .EQ. pLst%NodeIDs(J) ) pLst%ElmNds(K2,1)=2 !store whether first or second node of element  
                  ! --- Element Me, Ke, Fg, Fce
                  CALL ElemM(p%ElemProps(L),         pLst%Me(:,:,K2, 1))
                  CALL ElemK(p%ElemProps(L),         pLst%Ke(:,:,K2, 1))
                  CALL ElemF(p%ElemProps(L), Init%g, pLst%Fg(:,K2,1), FCe)
                  pLst%Fg(:,K2,1) = pLst%Fg(:,K2,1) + FCe(1:12) ! Adding cable element force 
                  EXIT   !We found the element for that node, exit loop on elements
              ENDIF
          ENDDO
      ENDDO
    ENDDO    
   ENDIF
   !_____________________________________REACTIONS_____________________________________________
   p%OutReact = .FALSE.
   DO I=1,p%NumOuts
      if ( ANY( p%OutParam(I)%Indx == ReactSS) ) THEN ! bjj: removed check of first 5 characters being "React" because (1) cases matter and (2) we can also ask for "-React*" or "mREACT"
         p%OutReact   =.TRUE.  
         EXIT
      ENDIF
   ENDDO
 
   IF (p%OutReact) THEN  !I need to store all constrained forces and moments; WE do not allow more than one member to be connected at a constrained joint for the time being

      ALLOCATE ( p%MOutLst3(p%nNodes_C), STAT = ErrStat2)     !this list contains different arrays for each of its elements
      ErrMsg2 = 'Error allocating p%MOutLst3 array in SDOut_Init'
      if(Failed()) return

      DO I=1,p%nNodes_C  !For all constrained node
         pLst => p%MOutLst3(I)
         pLst%Noutcnt=p%Nodes_C(I,1) !Assign nodeID for list I, I am using Noutcnt as a temporary holder for it, since nodeID is n array
         NconEls=Init%NodesConnE(pLst%Noutcnt,1) !Number of elements connecting to the joint
         ! ElmIDs: element IDs connecting to the joint; (1,NconEls) and not (NconEls) as the same meshauxtype is used with other MOutLst
         ! Me: has for each selected joint, and for each element attached to that node, a 12x12 matrix (extra dimension redundant)
         ! Ke: has for each selected joint, and for each element attached to that node  a 12x12 matrix
         CALL AllocAry(pLst%ElmIDs,      1, NconEls, ' p%MOutLst3(I)%ElmIds', ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%ElmNds,      1, NconEls, ' p%MOutLst3(I)%ElmNds', ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Me, 12, 12 , 1, NconEls, ' p%MOutLst3(I)%Me'    , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Ke, 12, 12 , 1, NconEls, ' p%MOutLst3(I)%Ke'    , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(pLst%Fg,     12 , 1, NconEls, ' p%MOutLst3(I)%Fg'    , ErrStat2, ErrMsg2); if(Failed()) return

         DO K=1, NconEls 
            L=Init%NodesConnE(pLst%Noutcnt,k+1)  !k-th Element Number in the set of elements attached to the selected node 
            pLst%ElmIDs(1,K)=L     !This array has for each joint requested  the elements' ID to get results for  
            M=p%Elems(L,2:3) !1st and 2nd node of the k-th element
            !Select whether the joint is the 1st or second node of the element
            pLst%ElmNds(1,K)=1                        !store whether first or second node of element  
            IF (M(2) .EQ. pLst%Noutcnt ) pLst%ElmNds(1,K)=2 !store whether first or second node of element  
            ! --- Element Me, Ke, Fg, Fce
            CALL ElemM(p%ElemProps(L),         pLst%Me(:,:,1,K))
            CALL ElemK(p%ElemProps(L),         pLst%Ke(:,:,1,K))
            CALL ElemF(p%ElemProps(L), Init%g, pLst%Fg(:,1,K), FCe)
            pLst%Fg(:,1,K) = pLst%Fg(:,1,K) + FCe(1:12) ! Adding cable element force 
         ENDDO
      ENDDO
      ! Compute p%TIreact, matrix to calculate single point reaction at the base of structure
      CALL ReactMatx(Init, WtrDpth, p, ErrStat, ErrMsg)
   ENDIF
 
   ! These variables are to help follow the framework template, but the data in them is simply a copy of data
   ! already available in the OutParam data structure
   CALL AllocAry(InitOut%WriteOutputHdr, p%NumOuts+p%OutAllint*p%OutAllDims, 'WriteOutputHdr', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(InitOut%WriteOutputUnt, p%NumOuts+p%OutAllint*p%OutAllDims, 'WriteOutputUnt', ErrStat2, ErrMsg2); if(Failed()) return
   DO I = 1,p%NumOuts+p%OutAllint*p%OutAllDims
      InitOut%WriteOutputHdr(I) = TRIM( p%OutParam(I)%Name  )
      InitOut%WriteOutputUnt(I) = TRIM( p%OutParam(I)%Units )      
   END DO  
   
   RETURN

CONTAINS
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SDOut_Init') 
        Failed =  ErrStat >= AbortErrLev
!         if (Failed) call CleanUp()
   END FUNCTION Failed
END SUBROUTINE SDOut_Init

!------------------------------------------------------------------------------------------------------
!>This subroutine allocates and calculated TIreact, Matrix to go from local reactions at constrained nodes to single point reactions
SUBROUTINE ReactMatx(Init, WtrDpth, p, ErrStat, ErrMsg)
   TYPE(SD_InitType),      INTENT(IN   ) :: Init    !< Input data for initialization routine
   REAL(ReKi),             INTENT(IN   ) :: WtrDpth !< Water depth
   TYPE(SD_ParameterType), INTENT(INOUT) :: p       !< Parameter data
   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER                  :: I !counter
   INTEGER(IntKi)           :: iDOF, iiDOF, iNode, nDOFPerNode !
   REAL(ReKi)               :: dx, dy, dz ! distances from reaction points to subdyn origin (mudline)
   REAL(ReKi), dimension(6) :: Line
   ErrStat=ErrID_None
   ErrMsg=""
   
   CALL AllocAry(p%TIreact, 6, p%nDOFC, 'TIReact', ErrStat, ErrMsg); if ( ErrStat /= ErrID_None ) return
   
   ! --- TI: Transformation matrix from interface points to ref point
   p%TIreact(1:6,:)=0 !Initialize
   DO I = 1, p%nDOFC
      iDOF = p%IDC(I) ! DOF index in constrained system
      iNode       = p%DOFtilde2Nodes(iDOF,1) ! First column is node 
      nDOFPerNode = p%DOFtilde2Nodes(iDOF,2) ! Second column is number of DOF per node
      iiDOF       = p%DOFtilde2Nodes(iDOF,3) ! Third column is dof index for this joint (1-6 for cantilever)
      if ((iiDOF<1) .or. (iiDOF>6)) then
         ErrMsg  = 'ReactMatx, interface node DOF number is not valid. DOF:'//trim(Num2LStr(iDOF))//' Node:'//trim(Num2LStr(iNode))//' iiDOF:'//trim(Num2LStr(iiDOF)); ErrStat = ErrID_Fatal
         return
      endif
      if (nDOFPerNode/=6) then
         ErrMsg  = 'ReactMatx, interface node doesnt have 6 DOFs. DOF:'//trim(Num2LStr(iDOF))//' Node:'//trim(Num2LStr(iNode))//' nDOF:'//trim(Num2LStr(nDOFPerNode)); ErrStat = ErrID_Fatal
         return
      endif
      
      dx = Init%Nodes(iNode, 2)          
      dy = Init%Nodes(iNode, 3)          
      dz = Init%Nodes(iNode, 4) + WtrDpth

      CALL RigidTransformationLine(dx,dy,dz,iiDOF,Line) !returns Line
      p%TIreact(1:6, I) = Line
   enddo

END SUBROUTINE ReactMatx

!====================================================================================================
!> Writes the data stored in the y variable to the correct indexed postions in WriteOutput
!! This is called by SD_CalcOutput() at each time step.
!! This routine does fill Allouts
!! note that this routine assumes m%u_TP and m%udotdot_TP have been set before calling 
!!     this routine (which is done in SD_CalcOutput() and SD CalcContStateDeriv)
!---------------------------------------------------------------------------------------------------- 
SUBROUTINE SDOut_MapOutputs( CurrentTime, u,p,x, y, m, AllOuts, ErrStat, ErrMsg )
   real(DbKi),                    intent( in    )  :: CurrentTime          ! Current simulation time in seconds
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
   real(ReKi), dimension (3,3)    :: DIRCOS    !direction cosice matrix (global to local) (3x3)
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
            IF (pLst%ElmIDs(iiNode,p%NavgEls) .NE. 0) THEN  ! Second element exist
               ! NOTE: forces are computed in the coordinate system of the first element for averaging
               call ElementForce(pLst, iiNode, p%NavgEls, FM_elm, FK_elm, sgn, DIRCOS, .true.) ! True= we use DIRCOS from element above
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
      ReactNs = 0.0 !Initialize
      DO I=1,p%nNodes_C   !Do for each constrained node, they are ordered as given in the input file and so as in the order of y2mesh
         FK_elm2=0. !Initialize for cumulative force
         FM_elm2=0. !Initialize
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
         !junk= FK_elm2 ! + FM_elm2  !removed the inertial component 12/13 !Not sure why I need an intermediate step here, but the sum would not work otherwise
         !NEED TO ADD HYDRODYNAMIC FORCES AT THE RESTRAINT NODES
         !   The joind iD of the reaction, i.e. thre reaction node ID is within p%MOutLst3(I)%Noutcnt
         !Since constrained nodes are ordered as given in the input file and so as in the order of y2mesh, i Can do:
         iSDNode   = p%Nodes_C(I,1)
         iMeshNode = p%INodes_SD_to_Mesh(iSDNode)
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
      integer(IntKi)       , intent(in)          :: iiNode !< Index over the nodes of a given number
      integer(IntKi)       , intent(in)          :: JJ     !< TODO: interpretation: index over other member connected to the current member (for averaging)
      real(ReKi), dimension (3,3), intent(inout) :: DIRCOS  !direction cosice matrix (global to local) (3x3)
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
      X_e(1:6)      = m%U_full       (p%NodesDOF(ElemNodes(1))%List(1:6)) 
      X_e(7:12)     = m%U_full       (p%NodesDOF(ElemNodes(2))%List(1:6)) 
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
      Real(ReKi), DIMENSION (3,3),   INTENT(IN)  :: DIRCOS    !direction cosice matrix (global to local) (3x3)
      Real(ReKi), DIMENSION (12,12), INTENT(IN)  :: Me,Ke    !element M and K matrices (12x12) in GLOBAL REFERENCE (DIRCOS^T K DIRCOS)
      Real(ReKi), DIMENSION (12),    INTENT(IN)  :: Udotdot, Y2, Fg     !acceleration and velocities, gravity forces
      Integer(IntKi),                INTENT(IN)  :: FirstOrSecond !1 or 2 depending on node of interest
      REAL(ReKi), DIMENSION (6),    INTENT(OUT)  :: FM_nod, FK_nod  !output static and dynamic forces and moments
      !Locals
      INTEGER(IntKi) :: L !counter
      REAL(DbKi), DIMENSION(12)                    :: FM_glb, FF_glb, FM_elm, FF_elm  ! temporary storage 

      FM_glb = matmul(Me,Udotdot)   ! GLOBAL REFERENCE
      FF_glb = matmul(Ke,Y2)        ! GLOBAL REFERENCE
      FF_glb = FF_glb - Fg          ! GLOBAL REFERENCE
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
      WRITE (UnSum,'(/,A/)', IOSTAT=Stat)  'This summary file was closed on '//CurDate()//' at '//CurTime()//'.'
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
   WRITE (UnSum,'(/,A/)', IOSTAT=ErrStat2)  'This summary file was generated by '//TRIM( SD_Prog%Name )//&
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


END MODULE SubDyn_Output
