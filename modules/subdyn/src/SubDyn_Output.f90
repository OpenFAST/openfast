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
   TYPE(SD_InitType),        INTENT( INOUT ) :: Init                 ! data needed to initialize the output module     
   TYPE(SD_OutputType),      INTENT( INOUT ) :: y                    ! SubDyn module's output data
   TYPE(SD_ParameterType),   INTENT( INOUT ) :: p                    ! SubDyn module paramters
   TYPE(SD_MiscVarType),     INTENT( INOUT ) :: misc                 ! SubDyn misc/optimization variables
   TYPE(SD_InitOutputType ), INTENT( INOUT ) :: InitOut              ! SubDyn module initialization output data
   REAL(ReKi),               INTENT( IN    ) :: WtrDpth              ! water depth from initialization routine  
   INTEGER,                       INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred           
   CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   ! Local variables
   INTEGER(IntKi)       :: ErrStat2      ! Error status of the operation
   CHARACTER(ErrMsgLen) :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                                   :: I,J,K,K2,L,NconEls   !Counters
   INTEGER(IntKi)                                   :: Junk  !Temporary Holders
   INTEGER(IntKi), Dimension(2)                     :: M   !counter for two nodes at a time
   INTEGER(IntKi)                                   :: eType
   REAL(ReKi) :: FCe(12) ! Pretension force from cable element
   ErrStat = 0      
   ErrMsg=""

   p%OutAllDims=12*p%Nmembers*2    !size of AllOut Member Joint forces

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

  !Store mapping between nodes and elements      
  CALL NodeCon(Init,p,ErrStat2, ErrMsg2); if(Failed()) return
  
  DO I=1,p%NMOutputs
   CALL AllocAry(p%MOutLst(I)%NodeIDs,    p%MOutLst(I)%NoutCnt           , 'MOutLst(I)%NodeIDs', ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(p%MOutLst(I)%ElmIDs,     p%MOutLst(I)%NoutCnt, p%NAvgEls, 'MOutLst(I)%ElmIDs' , ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(p%MOutLst(I)%ElmNds,     p%MOutLst(I)%NoutCnt, p%NAvgEls, 'MOutLst(I)%ElmNds' , ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(p%MOutLst(I)%Me, 12, 12, p%MOutLst(I)%NoutCnt, p%NAvgEls, 'MOutLst(I)%Me'     , ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(p%MOutLst(I)%Ke, 12, 12, p%MOutLst(I)%NoutCnt, p%NAvgEls, 'MOutLst(I)%Ke'     , ErrStat2, ErrMsg2); if(Failed()) return
   CALL AllocAry(p%MOutLst(I)%Fg,     12, p%MOutLst(I)%NoutCnt, p%NAvgEls, 'MOutLst(I)%Fg'     , ErrStat2, ErrMsg2); if(Failed()) return

   p%MOutLst(I)%NodeIDs=Init%MemberNodes(p%MOutLst(I)%MemberID,p%MOutLst(I)%NodeCnt)  !We are storing the actual node numbers corresponding to what the user ordinal number is requesting
   p%MOutLst(I)%ElmIDs=0  !Initialize to 0
   p%MOutLst(I)%ElmNds=0  !Initialize to 0

   DO J=1,p%MOutLst(I)%NoutCnt !Iterate on requested nodes for that member
      !I need to get at most 2 elements that belong to the same MOutLst(I)%MemberID
      !make use of MemberNodes and NodesConnE
      NconEls=Init%NodesConnE(p%MOutLst(I)%NodeIDs(J),1)!Number of elements connecting to the j-th node

      K2=0    !Initialize counter
      DO K=1, NconEls 
         L=Init%NodesConnE(p%MOutLst(I)%NodeIDs(J),k+1)  !k-th Element Number 
         M     = p%Elems(L,2:3) !1st and 2nd node of the k-th element
         eType = p%Elems(L, iMType)
            
         !Select only the other node, not the one where elements connect to
          IF (M(1) .EQ. p%MOutLst(I)%NodeIDs(J)) then
            Junk=M(2)
         else
            Junk=M(1)
         endif
                        
         IF (ANY(Init%MemberNodes(p%MOutLst(I)%MemberID,:) .EQ. Junk)) THEN  !This means we are in the selected member
            IF (K2 .EQ. 2) EXIT
            K2=K2+1
            p%MOutLst(I)%ElmIDs(J,K2)=L        !This array has for each node requested NODEID(J), for each memberMOutLst(I)%MemberID, the 2 elements to average from, it may have 1 if one of the numbers is 0 
            IF (M(2) .EQ. p%MOutLst(I)%NodeIDs(J) )then 
               p%MOutLst(I)%ElmNds(J,K2)=2 !store whether first or second node of element  
            else
               p%MOutLst(I)%ElmNds(J,K2)=1 !store whether first or second node of element  
            endif
            ! --- Element Me, Ke, Fg, Fce
            CALL ElemM(p%ElemProps(L),         p%MOutLst(I)%Me(:,:,J,K2))
            CALL ElemK(p%ElemProps(L),         p%MOutLst(I)%Ke(:,:,J,K2))
            CALL ElemF(p%ElemProps(L), Init%g, p%MOutLst(I)%Fg(:,J,K2), FCe)
            p%MOutLst(I)%Fg(:,J,K2) = p%MOutLst(I)%Fg(:,J,K2) + FCe(1:12) ! Adding cable element force 
         END IF    
      ENDDO  ! K, NconEls
     ENDDO !J, Noutcnt
   ENDDO  !I, NMOutputs

   END IF   ! there are any requested outputs   
 
   IF (p%OutAll) THEN  !I need to store all member end forces and moments 
     
    ALLOCATE ( p%MOutLst2(p%NMembers), STAT = ErrStat2 )     !this list contains different arrays for each of its elements
    ErrMsg2 = 'Error allocating p%MOutLst2 array in SDOut_Init'
     
    DO I=1,p%NMembers
      CALL AllocAry(p%MOutLst2(I)%NodeIDs, Init%Ndiv+1, 'MOutLst2(I)%NodeIDs', ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(p%MOutLst2(I)%ElmIDs,     2, 1, 'MOutLst2(I)%ElmIDs'     , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(p%MOutLst2(I)%ElmNds,     2, 1, 'MOutLst2(I)%ElmNds'     , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(p%MOutLst2(I)%Me, 12, 12, 2, 1, 'MOutLst2(I)%Me'         , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(p%MOutLst2(I)%Ke, 12, 12, 2, 1, 'MOutLst2(I)%Ke'         , ErrStat2, ErrMsg2); if(Failed()) return
      CALL AllocAry(p%MOutLst2(I)%Fg,     12, 2, 1, 'MOutLst2(I)%Fg'         , ErrStat2, ErrMsg2); if(Failed()) return

      p%MOutLst2(I)%MemberID=I !Assign memberID for all members
      p%MOutLst2(I)%NodeIDs=Init%MemberNodes(I,1:Init%Ndiv+1)  !We are storing  the actual node numbers in the member
      
      !Now I need to find out which elements are attached to those nodes and still belong to the member I
      !ElmIDs could contain the same element twice if Ndiv=1
      p%MOutLst2(I)%ElmIDs=0  !Initialize to 0

      DO J=1,Init%Ndiv+1,Init%Ndiv !Iterate on requested nodes for that member (first and last)
          !I need to get at most 2 elements that belong to the same I Member
          !make use of MemberNodes and NodesConnE
          
          NconEls=Init%NodesConnE(p%MOutLst2(I)%NodeIDs(J),1) !Number of elements connecting to the 1st or last node of the member
          
          K2= J/(Init%Ndiv+1)+1  !store this quantity used later, basically 1 or 2 depending on J
         
          DO K=1, NconEls 
              L=Init%NodesConnE(p%MOutLst2(I)%NodeIDs(J),k+1)  !k-th Element Number in the set of elements attached to the selected node 
              M=p%Elems(L,2:3) !1st and 2nd node of the k-th element
             !Select only the other node, not the one where elements connect to
             Junk=M(1)
             IF (M(1) .EQ. p%MOutLst2(I)%NodeIDs(J)) Junk=M(2)
             
             IF (ANY(Init%MemberNodes(p%MOutLst2(I)%MemberID,:) .EQ. Junk)) THEN  !This means we are in the selected member
                  p%MOutLst2(I)%ElmIDs(K2,1)=L     !This array has for each node requested NODEID(J), for each member I, the element to get results for 
                  p%MOutLst2(I)%ElmNds(K2,1)=1                        !store whether first or second node of element  
                  IF (M(2) .EQ. p%MOutLst2(I)%NodeIDs(J) ) p%MOutLst2(I)%ElmNds(K2,1)=2 !store whether first or second node of element  
                  ! --- Element Me, Ke, Fg, Fce
                  CALL ElemM(p%ElemProps(L),         p%MOutLst2(I)%Me(:,:,K2, 1))
                  CALL ElemK(p%ElemProps(L),         p%MOutLst2(I)%Ke(:,:,K2, 1))
                  CALL ElemF(p%ElemProps(L), Init%g, p%MOutLst2(I)%Fg(:,K2,1), FCe)
                  p%MOutLst2(I)%Fg(:,K2,1) = p%MOutLst2(I)%Fg(:,K2,1) + FCe(1:12) ! Adding cable element force 
                   
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

      ALLOCATE ( p%MOutLst3(p%NReact), STAT = ErrStat2)     !this list contains different arrays for each of its elements
      ErrMsg2 = 'Error allocating p%MOutLst3 array in SDOut_Init'
      if(Failed()) return

      DO I=1,p%NReact  !For all constrained node
         p%MOutLst3(I)%Noutcnt=p%Reacts(I,1) !Assign nodeID for list I, I am using Noutcnt as a temporary holder for it, since nodeID is n array
         NconEls=Init%NodesConnE(p%MOutLst3(I)%Noutcnt,1) !Number of elements connecting to the joint
         ! ElmIDs: element IDs connecting to the joint; (1,NconEls) and not (NconEls) as the same meshauxtype is used with other MOutLst
         ! Me: has for each selected joint, and for each element attached to that node, a 12x12 matrix (extra dimension redundant)
         ! Ke: has for each selected joint, and for each element attached to that node  a 12x12 matrix
         CALL AllocAry(p%MOutLst3(I)%ElmIDs,      1, NconEls, ' p%MOutLst3(I)%ElmIds', ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(p%MOutLst3(I)%ElmNds,      1, NconEls, ' p%MOutLst3(I)%ElmNds', ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(p%MOutLst3(I)%Me, 12, 12 , 1, NconEls, ' p%MOutLst3(I)%Me'    , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(p%MOutLst3(I)%Ke, 12, 12 , 1, NconEls, ' p%MOutLst3(I)%Ke'    , ErrStat2, ErrMsg2); if(Failed()) return
         CALL AllocAry(p%MOutLst3(I)%Fg,     12 , 1, NconEls, ' p%MOutLst3(I)%Fg'    , ErrStat2, ErrMsg2); if(Failed()) return

         DO K=1, NconEls 
            L=Init%NodesConnE(p%MOutLst3(I)%Noutcnt,k+1)  !k-th Element Number in the set of elements attached to the selected node 
            p%MOutLst3(I)%ElmIDs(1,K)=L     !This array has for each joint requested  the elements' ID to get results for  
            M=p%Elems(L,2:3) !1st and 2nd node of the k-th element
            !Select whether the joint is the 1st or second node of the element
            p%MOutLst3(I)%ElmNds(1,K)=1                        !store whether first or second node of element  
            IF (M(2) .EQ. p%MOutLst3(I)%Noutcnt ) p%MOutLst3(I)%ElmNds(1,K)=2 !store whether first or second node of element  
            ! --- Element Me, Ke, Fg, Fce
            CALL ElemM(p%ElemProps(L),         p%MOutLst3(I)%Me(:,:,1,K))
            CALL ElemK(p%ElemProps(L),         p%MOutLst3(I)%Ke(:,:,1,K))
            CALL ElemF(p%ElemProps(L), Init%g, p%MOutLst3(I)%Fg(:,1,K), FCe)
            p%MOutLst3(I)%Fg(:,1,K) = p%MOutLst3(I)%Fg(:,1,K) + FCe(1:12) ! Adding cable element force 
         ENDDO
      ENDDO
      !Store the matrix that will let me calculate single point reaction at the base of structure
      CALL ReactMatx(Init, p, WtrDpth, ErrStat, ErrMsg)
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
SUBROUTINE ReactMatx(Init, p, WtrDpth, ErrStat, ErrMsg)
   TYPE(SD_InitType),      INTENT(  IN)  :: Init         ! Input data for initialization routine
   TYPE(SD_ParameterType), INTENT(  INOUT)  :: p         ! Parameter data
   REAL(ReKi),                   INTENT(IN)     :: WtrDpth
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   ! local variables
   INTEGER                             :: I !counter
   INTEGER                             :: rmndr !type-column index
   INTEGER                             ::  n  !node ID
   INTEGER(IntKi)                      :: DOFC !  DOFC = Init%NReact*6
   REAL(ReKi)                          :: x, y, z !coordinates
   ErrStat=ErrID_None
   ErrMsg=""
   
   DOFC = p%NReact*6 ! bjj, this is p%DOFC    !Total DOFs at the base of structure 
   
   CALL AllocAry(p%TIreact, 6,DOFC, 'p%TIReact', ErrStat, ErrMsg )
   if ( ErrStat /= ErrID_None ) return
   
   p%TIreact=0 !Initialize
   
   DO I=1,3  !Take care of first three rows
      p%TIreact(I,I:DOFC:6)=1
   ENDDO 

    !Other rows done per column actually  
   DO I = 1, DOFC
      
      n = p%Reacts(ceiling(I/6.0),1)  !Constrained Node ID (this works in the reordered/renumbered p%Reacts) ! TODO different DOF ordering
      
      x = Init%Nodes(n, 2)
      y = Init%Nodes(n, 3)
      z = Init%Nodes(n, 4) + WtrDpth
      
      rmndr = MOD(I, 6)  !It gives me the column index among the 6 different kinds
      SELECT CASE (rmndr)
         CASE (1); p%TIreact(4:6, I) = (/0.0_ReKi , z        , -y/)
         CASE (2); p%TIreact(4:6, I) = (/-z       , 0.0_ReKi , x/)
         CASE (3); p%TIreact(4:6, I) = (/y        , -x       , 0.0_ReKi/)
         CASE (4); p%TIreact(4:6, I) = (/1.0_ReKi , 0.0_ReKi , 0.0_ReKi/)
         CASE (5); p%TIreact(4:6, I) = (/0.0_ReKi , 1.0_ReKi , 0.0_ReKi/)
         CASE (0); p%TIreact(4:6, I) = (/0.0_ReKi , 0.0_ReKi , 1.0_ReKi/)
            
         CASE DEFAULT
            ErrStat = ErrID_Fatal
            ErrMsg  = 'Error calculating transformation matrix TIreact, wrong column index '
            RETURN
         END SELECT
   ENDDO
END SUBROUTINE ReactMatx

!====================================================================================================
SUBROUTINE SDOut_MapOutputs( CurrentTime, u,p,x, y, m, AllOuts, ErrStat, ErrMsg )
! This subroutine writes the data stored in the y variable to the correct indexed postions in WriteOutput
! This is called by SD_CalcOutput() at each time step.
! This routine does fill Allouts
! note that this routine assumes m%u_TP and m%udotdot_TP have been set before calling 
!     this routine (which is done in SD_CalcOutput() and SD CalcContStateDeriv)
!---------------------------------------------------------------------------------------------------- 
   REAL(DbKi),                    INTENT( IN    )  :: CurrentTime          ! Current simulation time in seconds
   TYPE(SD_InputType),            INTENT( IN )     :: u                    ! SubDyn module's input data
   TYPE(SD_ContinuousStateType),  INTENT( IN )     :: x                    ! SubDyn module's states data
   TYPE(SD_OutputType),           INTENT( INOUT )  :: y                    ! SubDyn module's output data
   TYPE(SD_ParameterType), target,INTENT( IN    )  :: p                    ! SubDyn module's parameter data
   TYPE(SD_MiscVarType),          INTENT( INOUT )  :: m                    ! Misc/optimization variables
   REAL(ReKi),                    INTENT(   OUT )  :: AllOuts(0:MaxOutPts+p%OutAllInt*p%OutAllDims) ! Array of output data for all possible outputs
   INTEGER(IntKi),                INTENT(   OUT )  :: ErrStat              ! Error status of the operation
   CHARACTER(*),                  INTENT(   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   !locals
   INTEGER(IntKi)                           ::I,J,K,K2,L,L2      ! Counters
   INTEGER(IntKi), DIMENSION(2)             ::K3    ! It stores Node IDs for element under consideration (may not be consecutive numbers)
   INTEGER(IntKi)                           :: maxOutModes  ! maximum modes to output, the minimum of 99 or p%Nmodes
   REAL(ReKi), DIMENSION (6)                :: FM_elm, FK_elm, junk  !output static and dynamic forces and moments
   REAL(ReKi), DIMENSION (6)                :: FM_elm2, FK_elm2  !output static and dynamic forces and moments
   Real(ReKi), DIMENSION (3,3)              :: DIRCOS    !direction cosice matrix (global to local) (3x3)
   Real(ReKi), ALLOCATABLE                  :: ReactNs(:)    !6*Nreact reactions
   REAL(ReKi)                               :: Tmp_Udotdot(12), Tmp_y2(12) !temporary storage for calls to CALC_LOCAL
   
   Real(reKi), DIMENSION( p%URbarL+p%DOFL+6*p%Nreact)      :: yout            ! modifications to Y2 and Udotdot to include constrained node DOFs
   Real(ReKi),  DIMENSION(p%URbarL+p%DOFL+6*p%Nreact)      ::uddout           ! modifications to Y2 and Udotdot to include constrained node DOFs
   Integer(IntKi)                              ::sgn !+1/-1 for node force calculations
   type(MeshAuxDataType), pointer :: pLst
   ErrStat = ErrID_None   
   ErrMsg  = ""
   
   AllOuts = 0.0_ReKi  ! initialize for those outputs that aren't valid (and thus aren't set in this routine)
   
   !Create a variable that lists Y2 and adds removed constrained nodes' dofs; we will be using it to carry out other calculations with a special indexing array
   yout =0 !Initialize and populate with Y2 data  
   yout(1:         p%UrbarL       ) = m%UR_bar
   yout(p%URbarL+1:p%URbarL+p%DOFL) = m%UL
  
   !Same for a variable that deals with Udotdot
   uddout =0 !Initialize and populate with Udotdot data
   uddout(1          : p%URbarL         ) = m%UR_bar_dotdot
   uddout(p%URbarL+1 : p%URbarL+p%DOFL  ) = m%UL_dotdot

   ! TODO TODO TODO, there is a lot of similarity between the three outputs sections with some code redundency
         
      ! Only generate member-based outputs for the number of user-requested member outputs
      !Now store and identify needed output as requested by user
      !p%MOutLst has the mapping for the member, node, elements per node, to be used
      !MXNYZZZ   will need to connects to p%MOutLst(X)%ElmIDs(Y,1:2) if it is a force or accel; else to u%UFL(p%MOutLst(X)%NodeIDs(Y)) 
      !Inertial Load for the elements that are needed
  if (p%NumOuts > 0) then  !bjj: some of these fields aren't allocated when NumOuts==0
     DO I=1,p%NMOutputs
          !I know el # and whether it is 1st node or second node
        pLst=>p%MOutLst(I)
        DO J=1,pLst%NOutCnt !Iterate on requested nodes for that member 
             !I need to average across potentially up to 2 elements
             !Calculate forces on 1st stored element, and if 2nd exists do averaging with the second
             K = pLst%ElmIDs(J,1)  !element number
             K2= pLst%ElmNds(J,1)  !first or second node of the element to be considered
             !Assign the sign depending on whether it is the 1st or second node
             sgn=-1
             IF (K2 .EQ. 2) sgn= +1
             K3=p%Elems(K,2:3)  !first and second node ID associated with element K 
             L =p%IDY((K3(1)-1)*6+1)! starting index for node K3(1) within yout
             L2=p%IDY((K3(2)-1)*6+1)! starting index for node K3(2) within yout
             DIRCOS=transpose(p%elemprops(K)%DirCos)! global to local dir-cosine matrix
             !I need to find the Udotdot() for the two nodes of the element of interest 
             !I need to move the displacements to the local ref system
             ! bjj: added these temporary storage variables so that the CALC_NODE_FORCES call doesn't created
             !      new temporary *every time* it's called
             Tmp_Udotdot(1: 6) = uddout( L : L+5  )
             Tmp_Udotdot(7:12) = uddout( L2 : L2+5 )
             Tmp_y2(1: 6)      = yout( L : L+5 )
             Tmp_y2(7:12)      = yout( L2 : L2+5 )
             CALL CALC_NODE_FORCES( DIRCOS,pLst%Me(:,:,J,1),pLst%Ke(:,:,J,1),Tmp_Udotdot, Tmp_y2,pLst%Fg(:,J,1), K2,FM_elm,FK_elm) 
             
             FM_elm2=sgn*FM_elm
             FK_elm2=sgn*FK_elm
             
             IF (p%MOutLst(I)%ElmIDs(J,p%NavgEls) .NE. 0) THEN  !element number
                 K  = pLst%ElmIDs(J,p%NavgEls)  !element number
                 K2 = pLst%ElmNds(J,p%NavgEls)  !first or second node of the element to be considered
                 !Assign the sign depending on whether it is the 1st or second node
                 sgn=-1
                 IF (K2 .EQ. 2) sgn= +1
                 K3=p%Elems(K,2:3)  !first and second node ID associated with element K 
                 L  = p%IDY((K3(1)-1)*6+1)! starting index for node K3(1) within yout ! TODO different DOF order
                 L2 = p%IDY((K3(2)-1)*6+1)! starting index for node K3(2) within yout
                 CALL CALC_NODE_FORCES(DIRCOS,pLst%Me(:,:,J,p%NavgEls),pLst%Ke(:,:,J,p%NavgEls),(/uddout( L : L+5  ),uddout( L2 : L2+5 )/), &
                                 (/yout( L : L+5 ), yout( L2 : L2+5 )/), pLst%Fg(:,J,p%NavgEls), K2,FM_elm,FK_elm ) 
                                   
                 FM_elm2=0.5*( FM_elm2 + sgn*FM_elm ) !Now Average
                 FK_elm2=0.5*( FK_elm2 + sgn*FK_elm) !Now Average
             ENDIF
           
              ! Store in AllOuts
              !Forces and moments
              AllOuts(MNfmKe  (:,J,I))     = FK_elm2  !static forces and moments (6) Local Ref
              AllOuts(MNfmMe  (:,J,I))     = FM_elm2  !dynamic forces and moments (6) Local Ref
              !Displacement- Translational -no need for averaging since it is a node translation - In global reference SS
              L=p%IDY( (p%MOutLst(I)%NodeIDs(J)-1)*6 +1 )! starting index for nodeID(J) within yout
              AllOuts(MNTDss (:,J,I))      = yout(L:L+2)
              !Displacement- Rotational - I need to get the direction cosine matrix to tranform rotations  - In Local reference Element Ref Sys
              AllOuts(MNRDe (:,J,I))        = matmul(DIRCOS,yout(L+3:L+5)  ) !local ref
              !Accelerations- I need to get the direction cosine matrix to tranform displacement and rotations
              AllOuts(MNTRAe (1:3,J,I))     = matmul(DIRCOS,uddout(L:L+2)  )   !translational accel local ref
              AllOuts(MNTRAe (4:6,J,I))     = matmul(DIRCOS,uddout(L+3:L+5) )  !rotational accel  local ref
              
        ENDDO  ! J, Loop on requested nodes for that member
        
     ENDDO ! I, Loop on member outputs
   END IF
  
   IF (p%OutAll) THEN  !NEED TO CALCULATE TOTAL FORCES
      DO I=1,p%NMembers    !Cycle on all members
         pLst=>p%MOutLst2(I)
         DO J=1,2 !Iterate on requested nodes for that member (first and last)  
                K =pLst%ElmIDs(J,1)  !element number
                K2=pLst%ElmNds(J,1)  !first or second node of the element to be considered
                !Assign the sign depending on whether it is the 1st or second node
                sgn=-1
                IF (K2 .EQ. 2) sgn= +1
                K3=p%Elems(K,2:3)  !first and second node ID associated with element K 
                L =p%IDY((K3(1)-1)*6+1)! starting index for node K3(1) within yout ! TODO different DOF order
                L2=p%IDY((K3(2)-1)*6+1)! starting index for node K3(2) within yout
                DIRCOS=transpose(p%elemprops(K)%DirCos)! global to local
                CALL CALC_NODE_FORCES( DIRCOS, pLst%Me(:,:,J,1),pLst%Ke(:,:,J,1),(/uddout( L : L+5  ),uddout( L2 : L2+5 )/), &
                                 (/yout( L : L+5 ), yout( L2 : L2+5 )/), pLst%Fg(:,J,1),  K2,FM_elm,FK_elm) 
                 ! Store in All Outs
                 L=MaxOutPts+(I-1)*24+(J-1)*12+1!start index
                 L2=L+11
                 AllOuts( L:L2 ) =sgn* (/FK_elm,FM_elm/)
              ENDDO !J, nodes 1 and 2
         ENDDO ! I, Loop on members
  ENDIF
  
  !Assign interface forces and moments 
  AllOuts(IntfSS(1:TPdofL))= - (/y%Y1Mesh%Force (:,1), y%Y1Mesh%Moment(:,1)/) !-y%Y1  !Note this is the force that the TP applies to the Jacket, opposite to what the GLue Code needs thus "-" sign
  !Assign interface translations and rotations at the TP ref point  
  AllOuts(IntfTRss(1:TPdofL))=m%u_TP 
  !Assign interface translations and rotations accelerations
  AllOuts(IntfTRAss(1:TPdofL))= m%udotdot_TP 

  ! Assign all SSqm, SSqmdot, SSqmdotdot
  ! We only have space for the first 99 values
  maxOutModes = min(p%Nmodes,99)
  IF ( maxOutModes > 0 ) THEN 
     !BJJ: TODO: is there a check to see if we requested these channels but didn't request the modes? (i.e., retain 2 modes but asked for 75th mode?)
     Allouts(SSqm01  :SSqm01  +maxOutModes-1) = x%qm      (1:maxOutModes)
     Allouts(SSqmd01 :SSqmd01 +maxOutModes-1) = x%qmdot   (1:maxOutModes)
     Allouts(SSqmdd01:SSqmdd01+maxOutModes-1) = m%qmdotdot(1:maxOutModes)
  END IF
   
  !Need to Calculate Reaction Forces Now, but only if requested
  IF (p%OutReact) THEN 
       
       ALLOCATE ( ReactNs(6*p%NReact), STAT = ErrStat )
       IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for ReactNs array.'
         ErrStat = ErrID_Fatal
         RETURN
       END IF
       
       ReactNs = 0.0 !Initialize
       
       DO I=1,p%NReact   !Do for each constrained node, they are ordered as given in the input file and so as in the order of y2mesh
          FK_elm2=0. !Initialize for cumulative force
          FM_elm2=0. !Initialize
          pLst => p%MOutLst3(I)
          !Find the joint forces
             DO J=1,SIZE(pLst%ElmIDs(1,:))  !for all the elements connected (normally 1)
                K  = pLst%ElmIDs(1,J) ! element number
                K2 = pLst%ElmNds(1,J) ! 1=first, 2=second node of the element to be considered
                K3 = p%Elems(K,2:3)  !first and second node ID associated with element K 
                L  = p%IDY((K3(1)-1)*6+1)! starting index for node K3(1) within yout ! TODO different DOF order
                L2 = p%IDY((K3(2)-1)*6+1)! starting index for node K3(2) within yout
                DIRCOS=transpose(p%elemprops(K)%DirCos)! global to local
                CALL CALC_NODE_FORCES( DIRCOS,pLst%Me(:,:,1,J),pLst%Ke(:,:,1,J),(/uddout( L : L+5  ),uddout( L2 : L2+5 )/), &
                                 (/yout( L : L+5 ), yout( L2 : L2+5 )/),  pLst%Fg(:,1,J),   K2,FM_elm,FK_elm) 
                !transform back to global, need to do 3 at a time since cosine matrix is 3x3
                DO L=1,2  
                   FM_elm2((L-1)*3+1:L*3) = FM_elm2((L-1)*3+1:L*3) + matmul(p%elemprops(K)%DirCos,FM_elm((L-1)*3+1:L*3))  !sum forces at joint in GLOBAL REF
                   FK_elm2((L-1)*3+1:L*3) = FK_elm2((L-1)*3+1:L*3) + matmul(p%elemprops(K)%DirCos,FK_elm((L-1)*3+1:L*3))  !signs may be wrong, we will fix that later;  
                ! I believe this is all fixed in terms of signs now ,RRD 5/20/13
                ENDDO           
             ENDDO
             !junk= FK_elm2 ! + FM_elm2  !removed the inertial component 12/13 !Not sure why I need an intermediate step here, but the sum would not work otherwise
             !NEED TO ADD HYDRODYNAMIC FORCES AT THE RESTRAINT NODES
             !   The joind iD of the reaction, i.e. thre reaction node ID is within p%MOutLst3(I)%Noutcnt
             !The index in Y2mesh is? 
             !Since constrained nodes are ordered as given in the input file and so as in the order of y2mesh, i Can do:
             junk =  (/u%LMesh%Force(:,p%NNodes_I+p%NNodes_L+I),u%LMesh%Moment(:,p%NNodes_I+p%NNodes_L+I)/)
             ReactNs((I-1)*6+1:6*I)=FK_elm2 - junk  !Accumulate reactions from all nodes in GLOBAL COORDINATES
       ENDDO
       ! Store into AllOuts
       AllOuts( ReactSS(1:TPdofL) ) = matmul(p%TIreact,ReactNs)
  ENDIF
  if (allocated(ReactNs)) deallocate(ReactNs)

END SUBROUTINE SDOut_MapOutputs

!====================================================================================================
   SUBROUTINE CALC_NODE_FORCES(DIRCOS,Me,Ke,Udotdot,Y2 ,Fg, K2,FM_nod,FK_nod)
   !This function calculates for the given element the static and dynamic forces, given K and M of the element, and 
   !output quantities Udotdot and Y2 containing the 
   !and K2 indicating wheter the 1st (1) or 2nd (2) node is to be picked
!----------------------------------------------------------------------------------------------------
        Real(ReKi), DIMENSION (3,3),   INTENT(IN)  :: DIRCOS    !direction cosice matrix (global to local) (3x3)
        Real(ReKi), DIMENSION (12,12), INTENT(IN)  :: Me,Ke    !element M and K matrices (12x12) in GLOBAL REFERENCE (DIRCOS^T K DIRCOS)
        Real(ReKi), DIMENSION (12),    INTENT(IN)  :: Udotdot, Y2, Fg     !acceleration and velocities, gravity forces
        Integer(IntKi),                INTENT(IN)  :: K2   !1 or 2 depending on node of interest
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
        FM_nod=FM_elm(6*(k2-1)+1:k2*6) ! k2=1, 1:6,  k2=2  7:12 
        FK_nod=FF_elm(6*(k2-1)+1:k2*6) 
   
   END SUBROUTINE CALC_NODE_FORCES 

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
   DO k=p%Nmodes+1,99
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
!         RETURN
         p%OutParam(I)%Units = 'INVALID'  
         p%OutParam(I)%Indx  =  0
         p%OutParam(I)%SignM =  0                              ! this will print all zeros
      END IF
      
   END DO
   
   IF (p%OutAll) THEN   !Finish populating the OutParam with all the joint forces and moments
       ToTNames0=RESHAPE(SPREAD( (/"FKxe", "FKye", "FKze", "MKxe", "MKye", "MKze", "FMxe", "FMye", "FMze", "MMxe", "MMye", "MMze"/), 2, 2), (/24/) )
       ToTUnits=RESHAPE(SPREAD( (/"(N)  ","(N)  ","(N)  ", "(N*m)","(N*m)","(N*m)", "(N)  ","(N)  ","(N)  ", "(N*m)","(N*m)","(N*m)"/), 2, 2), (/24/) )
       DO I=1,p%Nmembers
           DO K=1,2
            DO J=1,12  !looks like I cnanot vectorize TRIM etc in Fortran
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
