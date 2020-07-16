!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
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
!**********************************************************************************************************************************
!> This code implements the spatial mapping algorithms described in the following papers
!! -  Sprague, Michael A.; Jonkman, Jason M.; and Jonkman, Bonnie J., "FAST Modular Framework for Wind Turbine Simulation: New 
!!    Algorithms and Numerical Examples." Proceedings of the 53rd Aerospace Sciences Meeting, 2015, also published in tech report 
!!    NREL/CP-2C00-63203, National Renewable Energy Laboratory, Golden, CO.  http://www.nrel.gov/docs/fy16osti/63203.pdf
!! -  Jonkman, Jason M. and Jonkman, Bonnie J., "FAST Modular Framework for Wind Turbine Simulation:  
!!    Full-System Linearization.", J. Phys.: Conf. Ser. 753 082010,  also published in tech report 
!!    NREL/CP-5000-67015, National Renewable Energy Laboratory, Golden, CO. http://www.nrel.gov/docs/fy17osti/67015.pdf 
   
!**********************************************************************************************************************************
MODULE ModMesh_Mapping

   USE ModMesh
   USE NWTC_LAPACK

   IMPLICIT NONE

   PRIVATE

   !bjj: these types require the use of ModMesh.f90, thus they cannot be part of NWTC_Library_Types.f90 (though they are auto_generated with that code):

      !> Type that describes characteristics of the mapping between two meshes
   TYPE, PUBLIC :: MapType
      INTEGER(IntKi) :: OtherMesh_Element    !< Node (for point meshes) or Element (for line2 meshes) number on other mesh; for loads, other mesh is Dest, for motions/scalars, other mesh is Src
      REAL(R8Ki)     :: distance             !< magnitude of couple_arm
      REAL(R8Ki)     :: couple_arm(3)        !< Vector between a point and node 1 of an element (p_ODR - p_OSR)
      REAL(R8Ki)     :: shape_fn(2)          !< shape functions: 1-D element-level location [0,1] based on closest-line projection of point
   END TYPE MapType

      !> data structures (for linearization) containing jacobians of mapping between fields on different meshes
   TYPE, PUBLIC :: MeshMapLinearizationType
         ! values for motions:
      REAL(R8Ki),     ALLOCATABLE :: mi(:,:)           !< block matrix of motions that reflects identity (i.e., solely the mapping of one quantity to itself on another mesh) [-]
      REAL(R8Ki),     ALLOCATABLE :: fx_p(:,:)         !< block matrix of motions that reflects skew-symmetric (cross-product) matrix [-]
      REAL(R8Ki),     ALLOCATABLE :: tv_uD(:,:)        !< block matrix of translational velocity that is multiplied by destination translational displacement [-]
      REAL(R8Ki),     ALLOCATABLE :: tv_uS(:,:)        !< block matrix of translational velocity that is multiplied by source translational displacement [-]
      REAL(R8Ki),     ALLOCATABLE :: ta_uD(:,:)        !< block matrix of translational acceleration that is multiplied by destination translational displacement [-]
      REAL(R8Ki),     ALLOCATABLE :: ta_uS(:,:)        !< block matrix of translational acceleration that is multiplied by source translational displacement [-]
      REAL(R8Ki),     ALLOCATABLE :: ta_rv(:,:)        !< block matrix of translational acceleration that is multiplied by omega (RotationVel) [-]      
         ! values for loads:
      REAL(R8Ki),     ALLOCATABLE :: li(:,:)           !< block matrix of loads that reflects identity (i.e., solely the mapping on one quantity to itself on another mesh) [-]
      REAL(R8Ki),     ALLOCATABLE :: M_uS(:,:)         !< block matrix of moment that is multiplied by Source u (translationDisp) [-]
      REAL(R8Ki),     ALLOCATABLE :: M_uD(:,:)         !< block matrix of moment that is multiplied by Destination u (translationDisp) [-]
      REAL(R8Ki),     ALLOCATABLE :: M_f(:,:)          !< block matrix of moment that is multiplied by force [-]
   END TYPE MeshMapLinearizationType
   
   
      !> data structures to determine full mapping between fields on different meshes
   TYPE, PUBLIC :: MeshMapType
      TYPE(MapType),  ALLOCATABLE :: MapLoads(:)               !< mapping data structure for loads on the mesh
      TYPE(MapType),  ALLOCATABLE :: MapMotions(:)             !< mapping data structure for motions and/or scalars on the mesh [-]
      TYPE(MapType),  ALLOCATABLE :: MapSrcToAugmt(:)          !< for source line2 loads, we map between source and an augmented source mesh, then between augmented source and destination
      TYPE(MeshType)              :: Augmented_Ln2_Src         !< the augmented source mesh needed for some mapping types
      TYPE(MeshType)              :: Lumped_Points_Src         !< a lumped mesh needed for some mapping types, stored here for efficiency
#ifdef MESH_DEBUG     
      TYPE(MeshType)              :: Lumped_Points_Dest        
#endif
      INTEGER,        ALLOCATABLE :: LoadLn2_A_Mat_Piv(:)      !< The pivot values for the factorization of LoadLn2_A_Mat
      REAL(R8Ki),     ALLOCATABLE :: DisplacedPosition(:,:,:)  !< couple_arm +Scr%Disp - Dest%Disp for each mapped node (stored here for efficiency)
      REAL(R8Ki),     ALLOCATABLE :: LoadLn2_A_Mat(:,:)        !< The n-by-n (n=3xNNodes) matrix that makes up the diagonal of the [A 0; B A] matrix in the point-to-line load mapping
      REAL(R8Ki),     ALLOCATABLE :: LoadLn2_F(:,:)            !< The 3-components of the forces for each node of an element in the point-to-line load mapping (for each element)
      REAL(R8Ki),     ALLOCATABLE :: LoadLn2_M(:,:)            !< The 3-components of the moments for each node of an element in the point-to-line load mapping (for each element)
      
      TYPE(MeshMapLinearizationType) :: dM                     !< type that contains information for linearization matrices, partial M partial u (or y)                  
   END TYPE MeshMapType
   
      ! note that these parameters must be negative (positive indicates the node/element it is mapped to)
   INTEGER(IntKi),  PARAMETER   :: NODE_NOT_MAPPED = -1        !< constant that indicates a node is not mapped

   PUBLIC :: MeshMapCreate
   PUBLIC :: MeshMapDestroy
   PUBLIC :: MeshMapWrBin
   PUBLIC :: Transfer_Point_to_Point
   PUBLIC :: Transfer_Line2_to_Point
   PUBLIC :: Transfer_Point_to_Line2
   PUBLIC :: Transfer_Line2_to_Line2
   
   PUBLIC :: Linearize_Point_to_Point
   PUBLIC :: Linearize_Line2_to_Point
   PUBLIC :: Linearize_Point_to_Line2
   PUBLIC :: Linearize_Line2_to_Line2
   !PUBLIC :: Lump_Line2_to_Point
   
   PUBLIC :: WriteMappingTransferToFile ! routine for mesh-mapping debugging
   
   ! auto-generated routines, necessary for the FAST Registry:
   PUBLIC :: NWTC_Library_DestroyMeshMapType, NWTC_Library_CopyMeshMapType, NWTC_Library_PackMeshMapType, NWTC_Library_UnpackMeshMapType
   PUBLIC :: NWTC_Library_DestroyMapType,     NWTC_Library_CopyMapType,     NWTC_Library_PackMapType,     NWTC_Library_UnpackMapType

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!bjj: maybe the MeshMapCreate routine shouldn't actually allocate arrays; allocate them in
!   the "IF (RemapFlag)" sections so that if people add nodes during the simulation, the structures get reallocated to correct 
!   size? MeshMapCreate should maybe be MeshMapping_Init() and only check that fields are compatible, etc.  
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine takes two meshes, determines the sizes required for the mapping data structure, and then
!! allocates the mappings (different for loads and motions/scalars).
SUBROUTINE MeshMapCreate( Src, Dest, MeshMap, ErrStat, ErrMsg )

! note that MeshMap%MapSrcToAugmt is allocated in Create_Augmented_Ln2_Src_Mesh() along with the Augmented_Ln2_Src Mesh

   TYPE(MeshType),           INTENT(IN)     ::  Src         !< source mesh
   TYPE(MeshType),           INTENT(IN)     ::  Dest        !< destination mesh
 
   TYPE(MeshMapType),        INTENT(INOUT)  :: MeshMap      !< mapping data structure

   INTEGER(IntKi),           INTENT(  OUT)  :: ErrStat      !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT)  :: ErrMsg       !< Error message if ErrStat /= ErrID_None

      ! local variables:

   INTEGER(IntKi)                           :: PointsInMap, PointsInTmpMap
   INTEGER(IntKi)                           :: ElementNodes
   LOGICAL                                  :: MapCreated
   INTEGER(IntKi)                           :: ErrStat2
   CHARACTER(ErrMsgLen)                     :: ErrMsg2
   CHARACTER(*), PARAMETER                  :: RoutineName = 'MeshMapCreate'
   

   ErrStat = ErrID_None
   ErrMsg  = ''

   MapCreated = .FALSE.


   IF ( .NOT. Dest%Committed .OR. .NOT. Src%Committed ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = " Both meshes must be committed before they can be mapped."
      RETURN
   END IF


   ElementNodes = 1
   PointsInTmpMap = 0
   
      !................................................
      ! Allocate the mapping for Motions and Scalars (if both meshes have some):
      !................................................
   IF ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) THEN

      IF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN !Line2-to-Point and Line2-to-Line2 motion mapping
         ElementNodes = 2      
      END IF
            
      
      ! for motion fields, every destination node is mapped to a source element or node 
      
      PointsInMap = Dest%Nnodes
      PointsInTmpMap = MAX(PointsInTmpMap,PointsInMap)

      IF ( PointsInMap < 1 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'MeshMap%MapMotions not allocated because no nodes were found to map.', ErrStat, ErrMsg, RoutineName)
      ELSE

            ! Allocate the mapping structure:
         ALLOCATE( MeshMap%MapMotions(PointsInMap), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error trying to allocate MeshMap%MapMotions.', ErrStat, ErrMsg, RoutineName)
         ELSE
            MapCreated = .TRUE.
            
               ! set up the initial mappings so that we don't necessarially have to do this multiple times on the first time step (if calculating Jacobians)
            IF ( Dest%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! point-to-Line2 or Line2-to-Line2
                  
               IF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! Line2-to-Line2         
                  CALL CreateMotionMap_L2_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               ELSEIF ( Src%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-Line2
                  CALL CreateMotionMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               END IF
         
            ELSEIF ( Dest%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-point or Line2-to-point
         
               IF ( Src%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-point            
                  CALL CreateMotionMap_P_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
               ELSEIF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! Line2-to-point         
                  CALL CreateMotionMap_L2_to_P(Src, Dest, MeshMap, ErrStat2, ErrMsg2)         
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
               END IF         
                  
            END IF ! create initial mapping based on mesh element type
                                    
         END IF ! MapMotions created

      END IF ! Dest has nodes to map

   END IF !HasMotionFields


      !................................................
      ! Allocate the mapping for Loads:
      !................................................
   IF ( HasLoadFields(Src) .AND. HasLoadFields(Dest) ) THEN

      ! check that the appropriate combinations of source/destination force/moments exist:
      IF ( Src%FieldMask(MASKID_Force) ) THEN
         IF (.NOT. Dest%FieldMask(MASKID_Force) ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Destination mesh does not contain force but source mesh does.', ErrStat, ErrMsg, RoutineName)
         END IF
         IF (.NOT. Dest%FieldMask(MASKID_Moment) ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Destination mesh must contain moment when source mesh contains force.', ErrStat, ErrMsg, RoutineName)
         END IF
      END IF
      IF ( Src%FieldMask(MASKID_Moment) ) THEN
         IF (.NOT. Dest%FieldMask(MASKID_Moment) ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Destination mesh does not contain moment but source mesh does.', ErrStat, ErrMsg, RoutineName)
         END IF
      END IF

      
      ! get size of mapping:
      PointsInMap = Src%Nnodes 

      IF ( PointsInMap < 1 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'MeshMap%MapLoads not allocated because no nodes were found to map.', ErrStat, ErrMsg, RoutineName)
      ELSE

            ! Allocate the mapping structure:
         ALLOCATE( MeshMap%MapLoads(PointsInMap), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error trying to allocate MeshMap%MapLoads.', ErrStat, ErrMsg, RoutineName)
         ELSE
            MapCreated = .TRUE.
            
               ! set up the initial mappings so that we don't necessarially have to do this multiple times on the first time step (if calculating Jacobians)
            IF ( Dest%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! point-to-Line2 or Line2-to-Line2
         
               ElementNodes = 2
         
               IF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! Line2-to-Line2         
                  CALL CreateLoadMap_L2_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               ELSEIF ( Src%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-Line2
                  CALL CreateLoadMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               END IF
         
            ELSEIF ( Dest%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-point or Line2-to-point
         
               IF ( Src%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-point            
                  CALL CreateLoadMap_P_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
               ELSEIF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! Line2-to-point         
                  CALL CreateLoadMap_L2_to_P(Src, Dest, MeshMap, ErrStat2, ErrMsg2)         
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
               END IF         
                  
            END IF ! create initial mapping based on mesh element type
                        
         END IF ! MapLoads allocated

      END IF ! Src has nodes to transfer
                  
   END IF ! HasLoadFields

   IF ( .NOT. MapCreated ) THEN
      CALL SetErrStat( ErrID_Fatal, 'Neither MapMotions or MapLoads was allocated. Meshes may not have compatible fields for mapping.', ErrStat, ErrMsg, RoutineName)
      RETURN
   END IF


      !................................................
      ! Allocate the DisplacedPosition field:
      !................................................

   IF (.NOT. ALLOCATED (MeshMap%DisplacedPosition)) THEN
      CALL AllocAry( MeshMap%DisplacedPosition, 3, PointsInTmpMap, ElementNodes, 'MeshMap%DisplacedPosition', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END IF


END SUBROUTINE MeshMapCreate
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine destroys the elements of the MeshMapType
SUBROUTINE MeshMapDestroy( MeshMap, ErrStat, ErrMsg )

   TYPE(MeshMapType),        INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),           INTENT(  OUT)  :: ErrStat      ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT)  :: ErrMsg       ! Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                           :: ErrStat2
   CHARACTER(ErrMsgLen)                     :: ErrMsg2


   ErrStat = ErrID_None
   ErrMsg  = ""

   CALL NWTC_Library_DestroyMeshMapType(MeshMap, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapDestroy')
         
#ifdef MESH_DEBUG     
  CALL MeshDestroy( MeshMap%Lumped_Points_Dest, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapDestroy')
#endif   

END SUBROUTINE MeshMapDestroy
!----------------------------------------------------------------------------------------------------------------------------------
!> This creates a matrix to write to a binary file (for debugging)
SUBROUTINE MeshMapWrBin( UnIn, Src, Dest, MeshMap, ErrStat, ErrMsg, FileName )

   INTEGER,                  INTENT(INOUT)  ::  UnIn     !< fortran output unit

   TYPE(MeshType),           INTENT(IN   )  :: Src       !< source mesh
   TYPE(MeshType),           INTENT(IN   )  :: Dest      !< destination mesh
   TYPE(MeshMapType),        INTENT(IN   )  :: MeshMap   !< mapping data structure

   INTEGER(IntKi),           INTENT(  OUT)  :: ErrStat   !< Error status of the operation
   CHARACTER(*),             INTENT(  OUT)  :: ErrMsg    !< Error message if ErrStat /= ErrID_None
   CHARACTER(*), OPTIONAL,   INTENT(IN   )  :: FileName  !< Name of the file to write the output in

      ! local variables
   REAL(ReKi),                  allocatable :: MapMat(:,:)
   INTEGER(IntKi)                           :: ErrStat2  ! Temporary storage for local errors
   INTEGER(IntKi)                           :: I , j , k, n       ! loop counter
   INTEGER(B4Ki)                            :: MeshesInFile(3), MatInFile(2)
   CHARACTER(ErrMsgLen)                     :: ErrMsg2

   ErrMsg = ""
   ErrStat = ErrID_None
   MeshesInFile = 0
   MatInFile= 0
   
   IF (UnIn < 0) THEN
      CALL GetNewUnit( UnIn, ErrStat, ErrMsg )

      CALL OpenBOutFile ( UnIn, TRIM(FileName), ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   
   IF ( MeshMap%Augmented_Ln2_Src%Committed ) THEN
      ALLOCATE( MapMat(MeshMap%Augmented_Ln2_Src%Nnodes, Dest%Nnodes),STAT=ErrStat2)
   ELSE
      ALLOCATE( MapMat(Src%Nnodes, Dest%Nnodes),STAT=ErrStat2)
   END IF
   IF ( ErrStat2 /=0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg="Error allocating MapMat."
      RETURN
   END IF
   
   MapMat = 0.0_ReKi
   IF ( HasMotionFields(Src) ) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(1,B4Ki)

      IF ( Src%ElemTable(ELEMENT_POINT)%Nelem > 0 ) THEN
         DO i=1,Dest%Nnodes
            j=MeshMap%MapMotions(i)%OtherMesh_Element
            if ( j < 1 )  CYCLE
            MapMat(j,i) = 1 
         END DO
      ELSEIF ( Src%ElemTable(ELEMENT_LINE2)%Nelem > 0 ) THEN
         DO i=1,Dest%Nnodes
            j=MeshMap%MapMotions(i)%OtherMesh_Element
            if ( j < 1 )  CYCLE
            DO k=1,size(Src%ElemTable(ELEMENT_LINE2)%Elements(j)%ElemNodes)
               n=Src%ElemTable(ELEMENT_LINE2)%Elements(j)%ElemNodes(k)
               MapMat(n,i) = MeshMap%MapMotions(i)%shape_fn(k)
            END DO
         END DO !i
      ENDIF
      
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(SIZE(MapMat,1),B4Ki)
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(SIZE(MapMat,2),B4Ki)
      WRITE (UnIn, IOSTAT=ErrStat2)   MapMat

   else
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(0,B4Ki)
   END IF
   
 
   IF ( HasLoadFields(Src) ) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(1,B4Ki) ! contains a load field
      MapMat = 0.0_ReKi
      IF ( Src%ElemTable(ELEMENT_POINT)%Nelem > 0 ) THEN
         DO i=1,Src%Nnodes
            j=MeshMap%MapLoads(i)%OtherMesh_Element
            if ( j < 1 )  CYCLE
            MapMat(i,j) = 1 
         END DO
      ELSEIF ( Src%ElemTable(ELEMENT_LINE2)%Nelem > 0 ) THEN
         
         IF ( Dest%ElemTable(ELEMENT_POINT)%Nelem > 0 ) THEN
            DO i=1,Src%Nnodes
               j=MeshMap%MapLoads(i)%OtherMesh_Element
               if ( j < 1 )  CYCLE
               DO k=1,size(Dest%ElemTable(ELEMENT_POINT)%Elements(j)%ElemNodes)
                  n=Dest%ElemTable(ELEMENT_POINT)%Elements(j)%ElemNodes(k)
                  MapMat(i,n) = 1 !this isn't really "1"... it's a factor of lumping, etc, but we'll leave as is for now.
               END DO
            END DO !i
         ELSEIF (Dest%ElemTable(ELEMENT_LINE2)%Nelem > 0 ) THEN   
            DO i=1,MeshMap%Augmented_Ln2_Src%Nnodes
               j=MeshMap%MapLoads(i)%OtherMesh_Element
               if ( j < 1 )  CYCLE
               DO k=1,size(Dest%ElemTable(ELEMENT_LINE2)%Elements(j)%ElemNodes)
                  n=Dest%ElemTable(ELEMENT_LINE2)%Elements(j)%ElemNodes(k)
                  MapMat(i,n) = MeshMap%MapLoads(i)%shape_fn(k)
               END DO
            END DO !i
         END IF
         
      ENDIF
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(SIZE(MapMat,1),B4Ki)
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(SIZE(MapMat,2),B4Ki)
      WRITE (UnIn, IOSTAT=ErrStat2)   MapMat
   else
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(0,B4Ki)
   END IF   
   
   
   IF ( ALLOCATED(MapMat) ) DEALLOCATE(MapMat)
   
   MeshesInFile = 0
#ifdef MESH_DEBUG     
  IF (MeshMap%Augmented_Ln2_Src%committed)  MeshesInFile(1)=1
  IF (MeshMap%Lumped_Points_Src%committed)  MeshesInFile(2)=1
  IF (MeshMap%Lumped_Points_Dest%committed)  MeshesInFile(3)=1
  WRITE (UnIn, IOSTAT=ErrStat2)   MeshesInFile

  IF (MeshMap%Augmented_Ln2_Src%committed)  CALL MeshWrBin(UnIn, MeshMap%Augmented_Ln2_Src,  ErrStat2, ErrMsg2)
  IF (MeshMap%Lumped_Points_Src%committed)  CALL MeshWrBin(UnIn, MeshMap%Lumped_Points_Src,  ErrStat2, ErrMsg2)
  IF (MeshMap%Lumped_Points_Dest%committed) CALL MeshWrBin(UnIn, MeshMap%Lumped_Points_Dest, ErrStat2, ErrMsg2)
#else
  WRITE (UnIn, IOSTAT=ErrStat2)   MeshesInFile
#endif

END SUBROUTINE MeshMapWrBin
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for transfering data from Line2 mesh to Point Mesh
SUBROUTINE Transfer_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),         INTENT(IN   ) ::  Src      !< source mesh
   TYPE(MeshType),         INTENT(INOUT) ::  Dest     !< destination mesh
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  SrcDisp  !< a "functional" sibling of the source mesh required for loads transfer; Src contains loads and SrcDisp contains TranslationDisp and Orientation
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  DestDisp !< a "functional" sibling of the destination mesh required for loads transfer; Dest contains loads and DestDisp contains TranslationDisp and Orientation

   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap   !< mapping data structure

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat   !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg    !< Error message if ErrStat /= ErrID_None


   REAL(ReKi)                            :: LoadsScaleFactor  ! bjj: added this scaling factor to get loads in a better numerical range 
   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*), PARAMETER               :: RoutineName = 'Transfer_Line2_to_Point'

   
   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

   !.................
   ! Check to ensure that the source mesh is composed of Line2 elements and destination mesh is composed of Point elements
   !.................   
   
   if (Src%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Source mesh must have one or more Line2 elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif

   if (Dest%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Destination mesh must have one or more Point elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif


   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------

   if ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) then
      
      ! This is the same algorithm as Transfer_Line2_to_Line2 (motions)
      
      !........................
      ! Start: Create Mapping data (if remap is true)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMotionMap_L2_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN

      endif !remapping


      !........................
      ! Start: Transfer data
      !........................

      CALL Transfer_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN

   endif !algorithm for motions/scalars


   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------

   ! ------------------------------------------------------------------------------------------------------
   ! Start: If mapping force/moment from Line2 to Point, need to created a temporary
   ! PointMesh that has nodes at the same locations as the Line2 mesh; distributed force/moment
   ! from the Line2 mesh can then be lumped into the temporary Point mesh and transfered to the Destination
   ! Point mesh via the same algorithm as subroutine Transfer_Point_to_Point
   ! ------------------------------------------------------------------------------------------------------

   if ( HasLoadFields(Src) ) then

      !........................
      ! Create mapping (including the temporary src mesh)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateLoadMap_L2_to_P(Src, Dest, MeshMap, ErrStat2, ErrMsg2)         
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
         

      ELSE ! Check that the temporary mesh has been set

         IF ( .NOT. MeshMap%Lumped_Points_Src%Initialized ) THEN
            CALL SetErrStat( ErrID_Fatal, 'MeshMap%Lumped_Points_Src not initialized (set RemapFlag = TRUE).', ErrStat, ErrMsg, RoutineName)
            RETURN
         END IF
         
      END IF

      IF (.not. PRESENT(SrcDisp) .OR. .NOT. PRESENT(DestDisp) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'SrcDisp and DestDisp arguments are required for load transfer.', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF
         
      LoadsScaleFactor = GetLoadsScaleFactor ( Src )
         
      ! first, we take the source fields and transfer them to fields on the augmented source mesh:
      !  (we're also taking the SrcDisp field and putting it on our augmented mesh)
      CALL Transfer_Src_To_Augmented_Ln2_Src( Src, MeshMap, ErrStat2, ErrMsg2, SrcDisp, LoadsScaleFactor ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  
      ! then we lump the loads from the augmented source mesh:
      CALL Lump_Line2_to_Point( MeshMap%Augmented_Ln2_Src,  MeshMap%Lumped_Points_Src,  ErrStat2, ErrMsg2, LoadsScaleFactor=LoadsScaleFactor  ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
         
         
      !........................
      ! Transfer data
      !........................
         
      ! we already checked if SrcDisp is present and transferred the displacements to MeshMap%Augmented_Ln2_Src
      CALL Transfer_Loads_Point_to_Point( MeshMap%Lumped_Points_Src, Dest, MeshMap, ErrStat2, ErrMsg2, MeshMap%Augmented_Ln2_Src, DestDisp, LoadsScaleFactor )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
         

   end if !algorithm for loads


END SUBROUTINE Transfer_Line2_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing linearization matrices for data transfer from Line2 mesh to Point Mesh
!! \copydetails modmesh_mapping::linearize_line2_to_line2 
SUBROUTINE Linearize_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),         INTENT(IN   ) ::  Src      !  source mesh
   TYPE(MeshType),         INTENT(IN   ) ::  Dest     !  destination mesh
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  SrcDisp  !  a "functional" sibling of the source mesh required for loads transfer; Src contains loads and SrcDisp contains TranslationDisp and Orientation
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  DestDisp !  a "functional" sibling of the destination mesh required for loads transfer; Dest contains loads and DestDisp contains TranslationDisp and Orientation

   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap   !  mapping data structure

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat   !  Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg    !  Error message if ErrStat /= ErrID_None


   real(R8Ki), allocatable               :: M_A(:,:)        ! linearization matrix for augmented source mesh
   real(R8Ki), allocatable               :: M_SL_fm(:,:)    ! linearization matrix for source-mesh lumped force component of moment
   real(R8Ki), allocatable               :: M_SL_uSm(:,:)   ! linearization matrix for source-mesh lumped translational displacement component of moment
   real(R8Ki), allocatable               :: M_SL_li(:,:)    ! linearization matrix for source-mesh lumped load "identity" component 
   
   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*), PARAMETER               :: RoutineName = 'Linearize_Line2_to_Point'

   
   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! ------------------------------------------------------------------------------------------------------------------------------
   !> ### Mapping and Linearization of Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------

   if ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) then
      
      ! This is the same algorithm as Transfer_Line2_to_Line2 (motions)
      
      !........................
      !> * Create Mapping data (if remap is true)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMotionMap_L2_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN

      endif !remapping


      !........................
      !> * Get linearization matrices of data transfer
      !........................

      CALL Linearize_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN

   endif !algorithm for motions/scalars


   ! ------------------------------------------------------------------------------------------------------------------------------
   !> ### Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------

   ! ------------------------------------------------------------------------------------------------------
   ! Start: If mapping force/moment from Line2 to Point, follow same algorithm as Line2-to-Line2, but skip final "unlumping" step
   ! ------------------------------------------------------------------------------------------------------

   if ( HasLoadFields(Src) ) then

      !........................
      !> * Create mapping (including the temporary src mesh)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateLoadMap_L2_to_P(Src, Dest, MeshMap, ErrStat2, ErrMsg2)         
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
         

      ELSE ! Check that the temporary mesh has been set

         IF ( .NOT. MeshMap%Lumped_Points_Src%Initialized ) THEN
            CALL SetErrStat( ErrID_Fatal, 'MeshMap%Lumped_Points_Src not initialized (set RemapFlag = TRUE).', ErrStat, ErrMsg, RoutineName)
            RETURN
         END IF
         
      END IF


      IF (.not. PRESENT(SrcDisp) .OR. .NOT. PRESENT(DestDisp) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'SrcDisp and DestDisp arguments are required for load transfer linearization.', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF
      
      !........................
      !> * Linearize data
      !........................
 
      !........................
      !>  + Get individual transformation matrices
      !........................
      
         
      !>   1. Get the matrix that transfers the source fields to the augmented source mesh.
      !! (We're also taking the force field and and translational displacement field and 
      !! putting them on our [intermediate] augmented mesh.)
      CALL Linearize_Src_To_Augmented_Ln2_Src( Src, MeshMap, ErrStat2, ErrMsg2, SrcDisp ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  
      call move_alloc( MeshMap%dM%li, M_A )    
! ^^^ size of M_A is 3*MeshMap%Augmented_Ln2_Src%NNodes X 3*Src%Nnodes            
            
      !>   2. Get the matrices that lump the loads on the augmented source mesh.
      !! (We're also taking the force field and putting it on our [intermediate] lumped mesh.)
      CALL Linearize_Lump_Line2_to_Point( MeshMap%Augmented_Ln2_Src,  MeshMap%Lumped_Points_Src, MeshMap%dM, ErrStat2, ErrMsg2  ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
       
      call move_alloc( MeshMap%dM%m_uD, M_SL_uSm )    
      call move_alloc( MeshMap%dM%m_f,  M_SL_fm )    
      call move_alloc( MeshMap%dM%li,   M_SL_li )    
         
! ^^^ size of M_SL_i (as well as M_SL_uSm and M_SL_fm) is 3*MeshMap%Augmented_Ln2_Src%NNodes X 3*MeshMap%Augmented_Ln2_Src%NNodes         

      !>   3. Get the matrices that transfer point meshes to other point meshes.
      ! we already checked if SrcDisp is present and transferred the displacements to MeshMap%Augmented_Ln2_Src
      CALL Linearize_Loads_Point_to_Point( MeshMap%Lumped_Points_Src, Dest, MeshMap, ErrStat2, ErrMsg2, MeshMap%Augmented_Ln2_Src, DestDisp )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) then
            call cleanup()
            return
         end if
         
!^^ this creates matrices of size 3*Dest%NNodes x 3*MeshMap%Lumped_Points_Src%NNodes (m_uD is square, of size 3*Dest%NNodes X 3*Dest%NNodes)
! need to return size 3*Dest%NNodes x 3*Src%NNodes            
            
      !........................
      !>  + Multiply individual transformation matrices to form full linearization matrices
      !........................
         
      CALL FormMatrix_FullLinearization( MeshMap%dM, M_A, M_SL_fm, M_SL_uSm, M_SL_li, ErrStat2, ErrMsg2 )    
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                                          
      call cleanup()               

   end if !algorithm for loads
   

contains
subroutine cleanup()

      if (allocated(M_A     )) deallocate(M_A    )
      if (allocated(M_SL_uSm)) deallocate(M_SL_uSm)   
      if (allocated(M_SL_fm )) deallocate(M_SL_fm)   
      if (allocated(M_SL_li )) deallocate(M_SL_li)   
               
end subroutine cleanup            


END SUBROUTINE Linearize_Line2_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that forms the final linearization matrices for of loads from line2 meshes to either line2 or point matrices
subroutine FormMatrix_FullLinearization( dM, M_A, M_SL_fm, M_SL_uSm, M_SL_li, ErrStat, ErrMsg )

   type(MeshMapLinearizationType), intent(inout) :: dM              !< linearization data type currently filled with values from point-to-point or point-to-line2 linearization
   real(R8Ki), allocatable,        intent(inout) :: M_A(:,:)        !< linearization matrix for augmented source mesh
   real(R8Ki), allocatable,        intent(inout) :: M_SL_fm(:,:)    !< linearization matrix for source-mesh lumped force component of moment
   real(R8Ki), allocatable,        intent(inout) :: M_SL_uSm(:,:)   !< linearization matrix for source-mesh lumped source translational displacement component of moment
   real(R8Ki), allocatable,        intent(inout) :: M_SL_li(:,:)    !< linearization matrix for source-mesh lumped load "identity" component 
      
   INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None                                   
                                   
   real(R8Ki), allocatable                       :: M(:,:)          ! temporary transfer matrix for linearization (to make sure everything is the correct size)
   
   INTEGER(IntKi)                                :: ErrStat2
   CHARACTER(ErrMsgLen)                          :: ErrMsg2
   CHARACTER(*), PARAMETER                       :: RoutineName = 'FormMatrix_FullLinearization'

   

   ErrStat = ErrID_None
   ErrMsg  = ''
                  
   
   !> Matrix \f$ M_{uSm} = \left[  M_{uSm}^D + M_{li}^D M_{uSm}^{SL} \right] M^A  \f$.  
   if (allocated(dM%m_uS)) then
                     
      dM%m_uS = dM%m_uS + matmul( dM%li, M_SL_uSm)  

      deallocate(M_SL_uSm)
         
      call AllocAry(M, size(dM%m_uS,1), size(M_A,2), 'M', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat>=AbortErrLev) return
            
      M = matmul( dM%m_uS, M_A )
      call move_alloc( M, dM%m_uS )
         
   end if
   
   !> Matrix \f$ M_{uDm} = M_{uDm}^D \f$     
   ! already stored in dM%m_uD
      
   !> Matrix \f$ M_{fm} = \left[  M_{fm}^D M_{li}^{SL} + M_{li}^D M_{fm}^{SL} \right] M^A  \f$      
   if (allocated(dM%m_f)) then

      dM%m_f =          matmul( dM%m_f, M_SL_li)         
      dM%m_f = dM%m_f + matmul( dM%li,  M_SL_fm) 
         
      deallocate(M_SL_fm)
         
      call AllocAry(M, size(dM%m_f,1), size(M_A,2), 'M', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat>=AbortErrLev) return
                              
      M = matmul( dM%m_f, M_A )
      call move_alloc( M, dM%m_f )
                  
   end if
      
      
   !> Matrix \f$ M_{li} = M_{li}^D M_{li}^{SL} M^A \f$      
   if (allocated(dM%li)) then
         
      dM%li = matmul( dM%li, M_SL_li)         
                           
      deallocate(M_SL_li)
         
      call AllocAry(M, size(dM%li,1), size(M_A,2), 'M', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat>=AbortErrLev) return
                                       
      M = matmul( dM%li, M_A )
      call move_alloc( M, dM%li )
                  
   end if
      
end subroutine FormMatrix_FullLinearization
!----------------------------------------------------------------------------------------------------------------------------------
#ifdef __NEW_LINE2_TO_POINT_MOTION_TRANSFER
!!!! 2-Jun-2016: bjj: this new Transfer_Motions_Line2_to_Point routine doesn't seem to give any better results than the old one;
!!!! it's more computationally expensive so we'll keep the old one.
!!!
!!!!> This subroutine returns an interpolated value of the source orientation, computed by
!!!!! \f$\ RotationMatrix =  
!!!!! \Lambda\left( \sum\limits_{i=1}^{2} \log\left( \left[\theta^{SR}_{eSn_i}\right]^T \right) \phi_i \right)
!!!!! \Lambda\left( \sum\limits_{i=1}^{2} \log\left( \theta^S_{eSn_i}                   \right) \phi_i \right)
!!!!! \f$
!!!!! where \f$\log()\f$ is nwtc_num::dcm_logmap and \f$\Lambda()\f$ is nwtc_num::dcm_exp
!!!subroutine InterpSrcOrientation(i, Src, MeshMap, RotationMatrixD, ErrStat, ErrMsg) 
!!!
!!!   INTEGER(IntKi),                 INTENT(IN   )  :: i         !< current destination node
!!!   TYPE(MeshType),                 INTENT(IN   )  :: Src       !< The source (Line2) mesh with motion fields allocated
!!!
!!!   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   !< The mapping data
!!!   REAL(DbKi),                     INTENT(  OUT)  :: RotationMatrixD(3,3) !< interpolated value of MATMUL( TRANSPOSE( Src%RefOrientation(:,:,n) ), Src%Orientation(:,:,n) )
!!!
!!!   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   !< Error status of the operation
!!!   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    !< Error message if ErrStat /= ErrID_None
!!!
!!!
!!!   REAL(DbKi)                                     :: FieldValue(3,2)                ! Temporary variable to store values for DCM interpolation
!!!   REAL(DbKi)                                     :: RotationMatrixDS(3,3)
!!!   REAL(DbKi)                                     :: tensor_interp(3)
!!!   INTEGER(IntKi)                                 :: n, n1, n2                      ! temporary space for node numbers
!!!
!!!   ErrStat = ErrID_None
!!!   ErrMsg  = ""
!!!
!!!   n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
!!!   n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)
!!!   
!!!   
!!!      ! bjj: because of numerical issues when the angle of rotation is pi, (where 
!!!      ! DCM_exp( DCM_logmap (x) ) isn't quite x
!!!   if ( EqualRealNos( MeshMap%MapMotions(i)%shape_fn(1), 1.0_ReKi ) ) then
!!!      RotationMatrixD = MATMUL( TRANSPOSE( Src%RefOrientation(:,:,n1) ), Src%Orientation(:,:,n1) )
!!!
!!!   elseif ( EqualRealNos( MeshMap%MapMotions(i)%shape_fn(2), 1.0_ReKi ) ) then
!!!      RotationMatrixD = MATMUL( TRANSPOSE( Src%RefOrientation(:,:,n2) ), Src%Orientation(:,:,n2) )
!!!
!!!   else
!!!         
!!!         !.........
!!!         ! interpolate reference orientation (bjj: could be optimized to not do this every step!)
!!!      do n1=1,NumNodes(ELEMENT_LINE2)
!!!         n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!               
!!!         RotationMatrixD = TRANSPOSE( Src%RefOrientation(:,:,n) )
!!!         CALL DCM_logmap( RotationMatrixD, FieldValue(:,n1), ErrStat, ErrMsg )                  
!!!            IF (ErrStat >= AbortErrLev) RETURN            
!!!      end do
!!!            
!!!      CALL DCM_SetLogMapForInterp( FieldValue )  ! make sure we don't cross a 2pi boundary         
!!!         
!!!         ! interpolate tensors: 
!!!      tensor_interp =   MeshMap%MapMotions(i)%shape_fn(1)*FieldValue(:,1)  &
!!!                      + MeshMap%MapMotions(i)%shape_fn(2)*FieldValue(:,2)                
!!!            
!!!         ! convert back to DCM for transpose (Src%RefOrientation):
!!!      RotationMatrixD = DCM_exp( tensor_interp )
!!!         
!!!      
!!!         !.........
!!!         ! interpolate Src%Orientation
!!!      do n1=1,NumNodes(ELEMENT_LINE2)
!!!         n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!               
!!!         RotationMatrixDS = Src%Orientation(:,:,n)    ! possible change of precision
!!!         CALL DCM_logmap( RotationMatrixDS, FieldValue(:,n1), ErrStat, ErrMsg )                  
!!!            IF (ErrStat >= AbortErrLev) RETURN
!!!      end do
!!!            
!!!      CALL DCM_SetLogMapForInterp( FieldValue )  ! make sure we don't cross a 2pi boundary         
!!!         
!!!         ! interpolate tensors: 
!!!      tensor_interp =   MeshMap%MapMotions(i)%shape_fn(1)*FieldValue(:,1)  &
!!!                      + MeshMap%MapMotions(i)%shape_fn(2)*FieldValue(:,2)                
!!!            
!!!         ! convert back to DCM for transpose (Src%RefOrientation):
!!!      RotationMatrixDS = DCM_exp( tensor_interp )            
!!!
!!!         ! multiply so we have matmul( transpose(Src%RefOrientation), transpose(Src%Orientation) )
!!!      RotationMatrixD = matmul( RotationMatrixD, RotationMatrixDS )                   
!!!      
!!!   end if  
!!!
!!!end subroutine InterpSrcOrientation
!!!!----------------------------------------------------------------------------------------------------------------------------------
!!!!> Given a mapping, this routine transfers the motions from nodes on Line2 elements to nodes on another mesh.
!!!SUBROUTINE Transfer_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )
!!!
!!!   TYPE(MeshType),                 INTENT(IN   )  :: Src       !< The source (Line2) mesh with motion fields allocated
!!!   TYPE(MeshType),                 INTENT(INOUT)  :: Dest      !< The destination mesh
!!!
!!!   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   !< The mapping data
!!!
!!!   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   !< Error status of the operation
!!!   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    !< Error message if ErrStat /= ErrID_None
!!!
!!!      ! local variables
!!!   INTEGER(IntKi)            :: i , j                          ! counter over the nodes
!!!   INTEGER(IntKi)            :: n, n1                          ! temporary space for node numbers
!!!   REAL(ReKi)                :: TmpVec(3), TmpVec2(3)
!!!   REAL(R8Ki)                :: RotationMatrix(3,3)
!!!
!!!   REAL(DbKi)                :: FieldValue(3,2)                ! Temporary variable to store values for DCM interpolation
!!!   REAL(DbKi)                :: RotationMatrixD(3,3)
!!!   
!!!
!!!   ErrStat = ErrID_None
!!!   ErrMsg  = ""
!!!
!!!  !> Define \f$ \phi_1 = 1-\bar{l}^S \f$  and 
!!!  !!        \f$ \phi_2 =   \bar{l}^S \f$.
!!!
!!!   
!!!   
!!!      ! ---------------------------- ORIENTATION/Direction Cosine Matrix   ----------------------
!!!      !> Orientation: \f$\theta^D = \theta^{DR} 
!!!      !! \Lambda\left( \sum\limits_{i=1}^{2} \log\left( \left[\theta^{SR}_{eSn_i}\right]^T \right) \phi_i \right)
!!!      !! \Lambda\left( \sum\limits_{i=1}^{2} \log\left( \theta^S_{eSn_i}                   \right) \phi_i \right)
!!!      !! \f$
!!!      !! where \f$\log()\f$ is nwtc_num::dcm_logmap and \f$\Lambda()\f$ is nwtc_num::dcm_exp
!!!      
!!!      ! transfer direction cosine matrix, aka orientation
!!!   if ( Src%FieldMask(MASKID_Orientation) .AND. Dest%FieldMask(MASKID_Orientation) ) then
!!!
!!!      do i=1, Dest%Nnodes
!!!         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
!!!
!!!         call InterpSrcOrientation(i, Src, MeshMap, RotationMatrixD, ErrStat, ErrMsg)
!!!            if (ErrStat >= AbortErrLev) return
!!!            
!!!      CALL DCM_logmap( RotationMatrixD, FieldValue(:,1), ErrStat, ErrMsg )                  
!!!            
!!!         RotationMatrix = REAL( RotationMatrixD, R8Ki )
!!!         Dest%Orientation(:,:,i) = MATMUL( Dest%RefOrientation(:,:,i), RotationMatrix  )
!!!
!!!         RotationMatrixD = Dest%RefOrientation(:,:,i)
!!!      CALL DCM_logmap( RotationMatrixD, FieldValue(:,2), ErrStat, ErrMsg )                  
!!!         
!!!      end do
!!!
!!!   endif   
!!!   
!!!             
!!!      ! ---------------------------- Translation ------------------------------------------------
!!!      !> Translational Displacement: \f$\vec{u}^D = 
!!!      !!   \sum\limits_{i=1}^{2}\left( \vec{u}^S_{eSn_i} \phi_i \right) +
!!!      !! \left[
!!!      !! \Lambda\left( \sum\limits_{i=1}^{2} \log\left( \left[\theta^{S}_{eSn_i}\right]^T \right) \phi_i \right)
!!!      !! \Lambda\left( \sum\limits_{i=1}^{2} \log\left(       \theta^{SR}_{eSn_i}         \right) \phi_i \right)
!!!      !!  - I \right]
!!!      !! \left\{ \vec{p}^{ODR}-
!!!      !!   \sum\limits_{i=1}^{2}\left( \vec{p}^{OSR}_{eSn_i} \phi_i \right) 
!!!      !! \right\}
!!!      !! \f$
!!!
!!!   if ( Src%FieldMask(MASKID_TranslationDisp) .AND. Dest%FieldMask(MASKID_TranslationDisp) ) then
!!!      do i=1, Dest%Nnodes
!!!         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
!!!                                                               
!!!         Dest%TranslationDisp(:,i) = 0.0_ReKi
!!!         do n1=1,NumNodes(ELEMENT_LINE2)
!!!            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!            Dest%TranslationDisp(:,i) = Dest%TranslationDisp(:,i) + MeshMap%MapMotions(i)%shape_fn(n1)*Src%TranslationDisp(:,n)
!!!         end do
!!!                           
!!!         
!!!            ! if Src mesh has orientation, superpose Dest displacement with translation due to rotation and couple arm
!!!         if ( Src%FieldMask(MASKID_Orientation) ) then
!!!
!!!               ! Calculate RotationMatrix as O_S^T*O_SR
!!!            call InterpSrcOrientation(i, Src, MeshMap, RotationMatrixD, ErrStat, ErrMsg)
!!!               if (ErrStat >= AbortErrLev) return
!!!            RotationMatrix = transpose( RotationMatrixD )               
!!!            
!!!               ! subtract I
!!!            do j=1,3
!!!               RotationMatrix(j,j)= RotationMatrix(j,j) - 1.0_ReKi
!!!            end do
!!!
!!!               ! multiply by p_DR - p_SR
!!!            
!!!            do n1=1,NumNodes(ELEMENT_LINE2)
!!!               n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!               FieldValue(:,n1) = MeshMap%MapMotions(i)%shape_fn(n1) * Src%Position(:,n)
!!!            end do
!!!            TmpVec = Dest%Position(:,i) - sum(FieldValue,2)
!!!            Dest%TranslationDisp(:,i) = Dest%TranslationDisp(:,i) + matmul( RotationMatrix, TmpVec )
!!!            
!!!         end if
!!!         
!!!      end do
!!!
!!!   end if
!!!
!!!
!!!
!!!      ! ---------------------------- Calculated total displaced positions  ---------------------
!!!      ! these values are used in both the translational velocity and translational acceleration
!!!      ! calculations. The calculations rely on the TranslationDisp fields, which are calculated
!!!      ! earlier in this routine.
!!!   IF ( Src%FieldMask(MASKID_TranslationVel) .OR. Src%FieldMask(MASKID_TranslationAcc) ) THEN
!!!      DO i = 1,Dest%Nnodes
!!!         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
!!!      
!!!         do n1=1,NumNodes(ELEMENT_LINE2)
!!!            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!            FieldValue(:,n1) = MeshMap%MapMotions(i)%shape_fn(n1) * ( Src%Position(:,n) + Src%TranslationDisp(:,n) )
!!!         end do         
!!!         TmpVec    = sum(FieldValue,2)
!!!         MeshMap%DisplacedPosition(:,i,1) = TmpVec - Dest%Position(:,i) - Dest%TranslationDisp(:,i)
!!!         
!!!      END DO   
!!!   END IF
!!!   
!!!      ! ---------------------------- TranslationVel  --------------------------------------------
!!!      !> Translational Velocity: \f$\vec{v}^D = 
!!!      !!   \sum\limits_{i=1}^{2}\left( \vec{v}^S_{eSn_i} \phi_i \right)
!!!      !! + \left\{
!!!      !! \sum\limits_{i=1}^{2}\left( \left\{ \vec{p}^{OSR}_{eSn_i} + \vec{u}^S_{eSn_i} \right\} \phi_i \right)
!!!      !!                           - \left\{ \vec{p}^{ODR}         + \vec{u}^D \right\} 
!!!      !! \right\} \times 
!!!      !! \sum\limits_{i=1}^{2} \left( \vec{\omega}^S_{eSn_i} \phi_i \right)
!!!      !! \f$
!!!   
!!!   
!!!   if ( Src%FieldMask(MASKID_TranslationVel) .AND. Dest%FieldMask(MASKID_TranslationVel) ) then
!!!      do i=1, Dest%Nnodes
!!!         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
!!!
!!!         Dest%TranslationVel(:,i) = 0.0_ReKi
!!!         do n1=1,NumNodes(ELEMENT_LINE2)
!!!            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!            Dest%TranslationVel(:,i) = Dest%TranslationVel(:,i) + MeshMap%MapMotions(i)%shape_fn(n1)*Src%TranslationVel(:,n)
!!!         end do
!!!                  
!!!                  
!!!         if ( Src%FieldMask(MASKID_RotationVel) ) then
!!!            
!!!            do n1=1,NumNodes(ELEMENT_LINE2)
!!!               n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!               FieldValue(:,n1) = MeshMap%MapMotions(i)%shape_fn(n1)*Src%RotationVel(:,n)
!!!            end do         
!!!            TmpVec = sum(FieldValue,2)
!!!            
!!!            Dest%TranslationVel(:,i) = Dest%TranslationVel(:,i) + cross_product( MeshMap%DisplacedPosition(:,i,1) , TmpVec)
!!!                        
!!!         endif
!!!
!!!      end do
!!!
!!!   endif
!!!
!!!      ! ---------------------------- RotationVel  -----------------------------------------------
!!!      !> Rotational Velocity: \f$\vec{\omega}^D = \sum\limits_{i=1}^{2} 
!!!      !!              \vec{\omega}^S_{eSn_i} 
!!!      !!              \phi_i\f$
!!!
!!!   
!!!   if ( Src%FieldMask(MASKID_RotationVel) .AND. Dest%FieldMask(MASKID_RotationVel) ) then
!!!      do i=1, Dest%Nnodes
!!!         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
!!!                                                                        
!!!         Dest%RotationVel(:,i) = 0.0_ReKi
!!!         do n1=1,NumNodes(ELEMENT_LINE2)
!!!            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!            Dest%RotationVel(:,i) = Dest%RotationVel(:,i) + MeshMap%MapMotions(i)%shape_fn(n1)*Src%RotationVel(:,n)
!!!         end do
!!!      end do
!!!
!!!   end if
!!!
!!!      ! ---------------------------- TranslationAcc -----------------------------------------------
!!!      !> Translational Acceleration: \f$\vec{a}^D = 
!!!      !!   \sum\limits_{i=1}^{2}\left( \vec{a}^S_{eSn_i} \phi_i \right)
!!!      !! + \left\{
!!!      !! \sum\limits_{i=1}^{2}\left( \left\{ \vec{p}^{OSR}_{eSn_i} + \vec{u}^S_{eSn_i} \right\} \phi_i \right)
!!!      !!                           - \left\{ \vec{p}^{ODR}         + \vec{u}^D \right\} 
!!!      !! \right\} \times 
!!!      !!    \left(  \sum\limits_{i=1}^{2}\left( \vec{\alpha}^S_{eSn_i} \phi_i \right) \right)
!!!      !! +  \left(  \sum\limits_{i=1}^{2}\left( \vec{\omega}^S_{eSn_i} \phi_i \right) \right) 
!!!      !! \times \left\{  \left\{ 
!!!      !! \sum\limits_{i=1}^{2}\left( \left\{ \vec{p}^{OSR}_{eSn_i} + \vec{u}^S_{eSn_i} \right\} \phi_i \right)
!!!      !!                           - \left\{ \vec{p}^{ODR}         + \vec{u}^D \right\} 
!!!      !! \right\} \times 
!!!      !! \left(  \sum\limits_{i=1}^{2}\left( \vec{\omega}^S_{eSn_i} \phi_i \right) \right) 
!!!      !! \right\}
!!!      !! \f$
!!!
!!!         
!!!   if ( Src%FieldMask(MASKID_TranslationAcc) .AND. Dest%FieldMask(MASKID_TranslationAcc) ) then
!!!      do i=1, Dest%Nnodes
!!!         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
!!!
!!!
!!!         Dest%TranslationAcc(:,i) = 0.0_ReKi
!!!         do n1=1,NumNodes(ELEMENT_LINE2)
!!!            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!            Dest%TranslationAcc(:,i) = Dest%TranslationAcc(:,i) + MeshMap%MapMotions(i)%shape_fn(n1)*Src%TranslationAcc(:,n)
!!!         end do
!!!         
!!!         if ( Src%FieldMask(MASKID_RotationAcc) )  then
!!!                        
!!!            do n1=1,NumNodes(ELEMENT_LINE2)
!!!               n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!               FieldValue(:,n1) = MeshMap%MapMotions(i)%shape_fn(n1)*Src%RotationAcc(:,n)
!!!            end do         
!!!            TmpVec = sum(FieldValue,2)
!!!            
!!!            Dest%TranslationAcc(:,i) = Dest%TranslationAcc(:,i) + cross_product( MeshMap%DisplacedPosition(:,i,1) , TmpVec)
!!!         end if
!!!         
!!!         if ( Src%FieldMask(MASKID_RotationVel) )  then
!!!            
!!!            do n1=1,NumNodes(ELEMENT_LINE2)
!!!               n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!               FieldValue(:,n1) = MeshMap%MapMotions(i)%shape_fn(n1)*Src%RotationVel(:,n)
!!!            end do         
!!!            TmpVec = sum(FieldValue,2) ! omega
!!!                                    
!!!            TmpVec2 = cross_product( MeshMap%DisplacedPosition(:,i,1), TmpVec )
!!!            
!!!            Dest%TranslationAcc(:,i) = Dest%TranslationAcc(:,i) + cross_product( TmpVec, TmpVec2)
!!!                                           
!!!         endif
!!!
!!!      end do
!!!   endif
!!!
!!!
!!!      ! ---------------------------- RotationAcc  -----------------------------------------------
!!!      !> Rotational Acceleration: \f$\vec{\alpha}^D = \sum\limits_{i=1}^{2} 
!!!      !!              \vec{\alpha}^S_{eSn_i} 
!!!      !!              \phi_i\f$
!!!
!!!   if (Src%FieldMask(MASKID_RotationAcc) .AND. Dest%FieldMask(MASKID_RotationAcc) ) then
!!!      do i=1, Dest%Nnodes
!!!         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
!!!
!!!         Dest%RotationAcc(:,i) = 0.0_ReKi
!!!         do n1=1,NumNodes(ELEMENT_LINE2)
!!!            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!            Dest%RotationAcc(:,i) = Dest%RotationAcc(:,i) + MeshMap%MapMotions(i)%shape_fn(n1)*Src%RotationAcc(:,n)
!!!         end do         
!!!         
!!!      end do
!!!   end if
!!!
!!!      ! ---------------------------- Scalars  -----------------------------------------------
!!!      !> Scalar: \f$S^D = \sum\limits_{i=1}^{2} 
!!!      !!              S^S_{eSn_i} 
!!!      !!              \phi_i\f$
!!!
!!!   if (Src%FieldMask(MASKID_SCALAR) .AND. Dest%FieldMask(MASKID_SCALAR) ) then
!!!      do i=1, Dest%Nnodes
!!!         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
!!!
!!!         Dest%Scalars(:,i) = 0.0_ReKi
!!!         do n1=1,NumNodes(ELEMENT_LINE2)
!!!            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(n1)
!!!         
!!!            Dest%Scalars(:,i) = Dest%Scalars(:,i) + MeshMap%MapMotions(i)%shape_fn(n1)*Src%Scalars(:,n)
!!!         end do
!!!                  
!!!      end do
!!!   end if
!!!
!!!
!!!END SUBROUTINE Transfer_Motions_Line2_to_Point
!!!!----------------------------------------------------------------------------------------------------------------------------------
#else
!> Given a mapping, this routine transfers the motions from nodes on Line2 elements to nodes on another mesh.
SUBROUTINE Transfer_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src       !< The source (Line2) mesh with motion fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest      !< The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   !< The mapping data

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)            :: i , j                          ! counter over the nodes
   INTEGER(IntKi)            :: k                              ! counter components
   INTEGER(IntKi)            :: nScalars                       ! number of scalars transferred
   INTEGER(IntKi)            :: n, n1, n2                      ! temporary space for node numbers
   REAL(R8Ki)                :: FieldValueN1(3)                ! Temporary variable to store field values on element nodes
   REAL(R8Ki)                :: FieldValueN2(3)                ! Temporary variable to store field values on element nodes
   REAL(ReKi)                :: TmpVec(3)
   REAL(R8Ki)                :: RotationMatrix(3,3)

   REAL(DbKi)                :: FieldValue(3,2)                ! Temporary variable to store values for DCM interpolation
   REAL(DbKi)                :: RotationMatrixD(3,3)
   REAL(DbKi)                :: tensor_interp(3)
   

   ErrStat = ErrID_None
   ErrMsg  = ""

  !> Define \f$ \phi_1 = 1-\bar{l}^S \f$  and 
  !!        \f$ \phi_2 =   \bar{l}^S \f$.

!bjj: FieldValueN1 and FieldValueN2 should really be one matrix of DIM (3,2) now that we've modified some of the other data structures....
             
      ! ---------------------------- Translation ------------------------------------------------
      !> Translational Displacement: \f$\vec{u}^D = \sum\limits_{i=1}^{2}\left( 
      !!              \vec{u}^S_{eSn_i} + \left[\left[\theta^S_{eSn_i}\right]^T \theta^{SR}_{eSn_i} - I\right]\left\{\vec{p}^{ODR}-\vec{p}^{OSR}_{eSn_i}\right\}
      !!              \right) \phi_i\f$

      ! u_Dest1 = u_Src + [Orientation_Src^T * RefOrientation_Src - I] * [p_Dest - p_Src] at Source Node n1
      ! u_Dest2 = u_Src + [Orientation_Src^T * RefOrientation_Src - I] * [p_Dest - p_Src] at Source Node n2
      ! u_Dest = (1.-elem_position)*u_Dest1 + elem_position*u_Dest2
   if ( Src%FieldMask(MASKID_TranslationDisp) .AND. Dest%FieldMask(MASKID_TranslationDisp) ) then
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

            ! add the translation displacement portion part
         do j=1,NumNodes(ELEMENT_LINE2) ! number of nodes per line2 element
            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
         
            FieldValue(:,j) = Src%TranslationDisp(:,n)
         end do
            

            ! if Src mesh has orientation, superpose Dest displacement with translation due to rotation and couple arm
         if ( Src%FieldMask(MASKID_Orientation) ) then

            do j=1,NumNodes(ELEMENT_LINE2) ! number of nodes per line2 element
               n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
         
                  !Calculate RotationMatrix as O_S^T*O_SR
               RotationMatrix = TRANSPOSE( Src%Orientation(:,:,n ) )
               RotationMatrix = MATMUL( RotationMatrix, Src%RefOrientation(:,:,n) )

                  ! subtract I
               do k=1,3
                  RotationMatrix(k,k)= RotationMatrix(k,k) - 1.0_R8Ki
               end do
               
               FieldValue(:,j) = FieldValue(:,j) + MATMUL(RotationMatrix,(Dest%Position(:,i)-Src%Position(:,n)))
                              
            end do
                        
         end if

            ! now form a weighted average of the two points:
         Dest%TranslationDisp(:,i) = MeshMap%MapMotions(i)%shape_fn(1)*FieldValue(:,1)  &
                                   + MeshMap%MapMotions(i)%shape_fn(2)*FieldValue(:,2)

      end do

   end if

      ! ---------------------------- ORIENTATION/Direction Cosine Matrix   ----------------------
      !> Orientation: \f$\theta^D = \Lambda\left( \sum\limits_{i=1}^{2} 
      !!              \log\left( \theta^{DR}\left[\theta^{SR}_{eSn_i}\right]^T\theta^S_{eSn_i} \right)
      !!              \phi_i \right)\f$
      !! where \f$\log()\f$ is nwtc_num::dcm_logmap and \f$\Lambda()\f$ is nwtc_num::dcm_exp
   
   
   
      ! transfer direction cosine matrix, aka orientation

   if ( Src%FieldMask(MASKID_Orientation) .AND. Dest%FieldMask(MASKID_Orientation) ) then

      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
                  
         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

            ! bjj: added this IF statement because of numerical issues when the angle of rotation is pi, 
            !      (where DCM_exp( DCM_logmap (x) ) isn't quite x
         if ( EqualRealNos( MeshMap%MapMotions(i)%shape_fn(1), 1.0_R8Ki ) ) then
            
            RotationMatrixD = MATMUL( TRANSPOSE( Src%RefOrientation(:,:,n1) ), Src%Orientation(:,:,n1) )
            RotationMatrixD = MATMUL( Dest%RefOrientation(:,:,i), RotationMatrixD )
      
         elseif ( EqualRealNos( MeshMap%MapMotions(i)%shape_fn(2), 1.0_R8Ki ) ) then
            
            RotationMatrixD = MATMUL( TRANSPOSE( Src%RefOrientation(:,:,n2) ), Src%Orientation(:,:,n2) )
            RotationMatrixD = MATMUL( Dest%RefOrientation(:,:,i), RotationMatrixD )
      
         else

               ! calculate Rotation matrix for FieldValueN1 and convert to tensor:
            RotationMatrixD = MATMUL( TRANSPOSE( Src%RefOrientation(:,:,n1) ), Src%Orientation(:,:,n1) )
            RotationMatrixD = MATMUL( Dest%RefOrientation(:,:,i), RotationMatrixD )

            CALL DCM_logmap( RotationMatrixD, FieldValue(:,1), ErrStat, ErrMsg )
            IF (ErrStat >= AbortErrLev) RETURN

               ! calculate Rotation matrix for FieldValueN2 and convert to tensor:
            RotationMatrixD = MATMUL( TRANSPOSE( Src%RefOrientation(:,:,n2) ), Src%Orientation(:,:,n2) )
            RotationMatrixD = MATMUL( Dest%RefOrientation(:,:,i), RotationMatrixD )
         
            CALL DCM_logmap( RotationMatrixD, FieldValue(:,2), ErrStat, ErrMsg )                  
            IF (ErrStat >= AbortErrLev) RETURN
         
            CALL DCM_SetLogMapForInterp( FieldValue )  ! make sure we don't cross a 2pi boundary
         
         
               ! interpolate tensors: 
            tensor_interp =   MeshMap%MapMotions(i)%shape_fn(1)*FieldValue(:,1)  &
                            + MeshMap%MapMotions(i)%shape_fn(2)*FieldValue(:,2)    
                  
               ! convert back to DCM:
            RotationMatrixD = DCM_exp( tensor_interp )
                        
         end if
         
         Dest%Orientation(:,:,i) = REAL( RotationMatrixD, R8Ki )
             
      end do

   endif

      ! ---------------------------- Calculated total displaced positions  ---------------------
      ! these values are used in both the translational velocity and translational acceleration
      ! calculations. The calculations rely on the TranslationDisp fields, which are calculated
      ! earlier in this routine.
   IF ( Src%FieldMask(MASKID_TranslationVel) .OR. Src%FieldMask(MASKID_TranslationAcc) ) THEN
      DO i = 1,Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
      
         DO j=1,NumNodes(ELEMENT_LINE2) ! number of nodes per line2 element
            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
         
            MeshMap%DisplacedPosition(:,i,j) =    Src%Position(:,n) +  Src%TranslationDisp(:,n)  &
                                               - Dest%Position(:,i) - Dest%TranslationDisp(:,i)  
         end do
      
      END DO   
   END IF
   
      ! ---------------------------- TranslationVel  --------------------------------------------
      !> Translational Velocity: \f$\vec{v}^D = \sum\limits_{i=1}^{2}\left( 
      !!              \vec{v}^S_{eSn_i} 
      !!              + \left\{ \left\{ \vec{p}^{OSR}_{eSn_i} + \vec{u}^S_{eSn_i} \right\} - \left\{ \vec{p}^{ODR} + \vec{u}^D \right\} \right\} \times \vec{\omega}^S_{eSn_i}
      !!              \right) \phi_i\f$
   
   
   if ( Src%FieldMask(MASKID_TranslationVel) .AND. Dest%FieldMask(MASKID_TranslationVel) ) then
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         FieldValueN1 = Src%TranslationVel(:,n1)
         FieldValueN2 = Src%TranslationVel(:,n2)

         if ( Src%FieldMask(MASKID_RotationVel) ) then
            FieldValueN1 = FieldValueN1 + cross_product ( MeshMap%DisplacedPosition(:,i,1), Src%RotationVel(:,n1) )
            FieldValueN2 = FieldValueN2 + cross_product ( MeshMap%DisplacedPosition(:,i,2), Src%RotationVel(:,n2) )
         endif

         Dest%TranslationVel(:,i) = MeshMap%MapMotions(i)%shape_fn(1)*FieldValueN1  &
                                  + MeshMap%MapMotions(i)%shape_fn(2)*FieldValueN2


      end do

   endif

      ! ---------------------------- RotationVel  -----------------------------------------------
      !> Rotational Velocity: \f$\vec{\omega}^D = \sum\limits_{i=1}^{2} 
      !!              \vec{\omega}^S_{eSn_i} 
      !!              \phi_i\f$

   
   if ( Src%FieldMask(MASKID_RotationVel) .AND. Dest%FieldMask(MASKID_RotationVel) ) then
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         Dest%RotationVel(:,i) = MeshMap%MapMotions(i)%shape_fn(1)*Src%RotationVel(:,n1)  &
                               + MeshMap%MapMotions(i)%shape_fn(2)*Src%RotationVel(:,n2)
      end do

   end if

      ! ---------------------------- TranslationAcc -----------------------------------------------
      !> Translational Acceleration: \f$\vec{a}^D = \sum\limits_{i=1}^{2}\left( 
      !!              \vec{a}^S_{eSn_i} 
      !!            + \left\{ \left\{ \vec{p}^{OSR}_{eSn_i} + \vec{u}^S_{eSn_i} \right\} - \left\{ \vec{p}^{ODR} + \vec{u}^D \right\} \right\} \times \vec{\alpha}^S_{eSn_i}
      !!            + \vec{\omega}^S_{eSn_i} \times \left\{
      !!              \left\{ \left\{ \vec{p}^{OSR}_{eSn_i} + \vec{u}^S_{eSn_i} \right\} - \left\{ \vec{p}^{ODR} + \vec{u}^D \right\} \right\} \times \vec{\omega}^S_{eSn_i}
      !!              \right\}
      !!              \right) \phi_i\f$

   
   
   
   if ( Src%FieldMask(MASKID_TranslationAcc) .AND. Dest%FieldMask(MASKID_TranslationAcc) ) then
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         FieldValueN1 = Src%TranslationAcc(:,n1)
         FieldValueN2 = Src%TranslationAcc(:,n2)


         if ( Src%FieldMask(MASKID_RotationAcc) )  then
            FieldValueN1 = FieldValueN1 + cross_product( MeshMap%DisplacedPosition(:,i,1), Src%RotationAcc(:,n1) )
            FieldValueN2 = FieldValueN2 + cross_product( MeshMap%DisplacedPosition(:,i,2), Src%RotationAcc(:,n2) )
         endif

         if ( Src%FieldMask(MASKID_RotationVel) )  then
            TmpVec = cross_product( MeshMap%DisplacedPosition(:,i,1), Src%RotationVel(:,n1) )
            FieldValueN1 =  FieldValueN1 + cross_product( Src%RotationVel(:,n1), TmpVec )
            
            TmpVec = cross_product( MeshMap%DisplacedPosition(:,i,2), Src%RotationVel(:,n2) )
            FieldValueN2 =  FieldValueN2 + cross_product( Src%RotationVel(:,n2), TmpVec )
                                           
         endif

         Dest%TranslationAcc(:,i) = MeshMap%MapMotions(i)%shape_fn(1)*FieldValueN1  &
                                  + MeshMap%MapMotions(i)%shape_fn(2)*FieldValueN2

      end do
   endif


      ! ---------------------------- RotationAcc  -----------------------------------------------
      !> Rotational Acceleration: \f$\vec{\alpha}^D = \sum\limits_{i=1}^{2} 
      !!              \vec{\alpha}^S_{eSn_i} 
      !!              \phi_i\f$

   if (Src%FieldMask(MASKID_RotationAcc) .AND. Dest%FieldMask(MASKID_RotationAcc) ) then
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         Dest%RotationAcc(:,i) = MeshMap%MapMotions(i)%shape_fn(1)*Src%RotationAcc(:,n1)  &
                               + MeshMap%MapMotions(i)%shape_fn(2)*Src%RotationAcc(:,n2)


      end do
   end if

      ! ---------------------------- Scalars  -----------------------------------------------
      !> Scalar: \f$S^D = \sum\limits_{i=1}^{2} 
      !!              S^S_{eSn_i} 
      !!              \phi_i\f$

   if (Src%FieldMask(MASKID_SCALAR) .AND. Dest%FieldMask(MASKID_SCALAR) ) then
      nScalars = min(Dest%nScalars, Src%nScalars)
      
      if (Dest%nScalars > nScalars) then
         call SetErrStat(ErrID_Severe, "Not all scalars could be computed from source mesh (insufficient data).", ErrStat, ErrMsg, 'Transfer_Motions_Line2_to_Point')
         Dest%Scalars(nScalars+1:,:) = 0.0_ReKi
      end if
         
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         Dest%Scalars(1:nScalars,i) = MeshMap%MapMotions(i)%shape_fn(1)*Src%Scalars(1:nScalars,n1)  &
                                    + MeshMap%MapMotions(i)%shape_fn(2)*Src%Scalars(1:nScalars,n2)

      end do
   end if


END SUBROUTINE Transfer_Motions_Line2_to_Point
#endif 
!----------------------------------------------------------------------------------------------------------------------------------
!> Given a mapping, this routine transfers the motions from nodes on Line2 elements to nodes on another mesh.
SUBROUTINE Linearize_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src       !< The source (Line2) mesh with motion fields allocated
   TYPE(MeshType),                 INTENT(IN   )  :: Dest      !< The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   !< The mapping data

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   integer(intKi)                          :: i,j, k, n, d_start, d_end, s_start, s_end
   real(R8Ki)                              :: tmp, tmpVec(3), SSMat(3,3)
   real(R8Ki)                              :: RotVel(3)
   
   character(*), parameter :: RoutineName = 'Linearize_Motions_Line2_to_Point'
   ErrStat = ErrID_None
   ErrMsg  = ""

   
   if (.not. allocated(MeshMap%dM%mi) ) then
      call AllocAry(MeshMap%dM%mi, Dest%Nnodes*3, Src%Nnodes*3, 'dM%mi', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
   end if
      
   MeshMap%dM%mi = 0.0_R8Ki
   do i=1, Dest%Nnodes
                  
      do j=1,NumNodes(ELEMENT_LINE2)
            
         n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
            
         do k=1,3
            MeshMap%dM%mi( (i-1)*3+k, (n-1)*3+k ) = MeshMap%MapMotions(i)%shape_fn(j)
         end do   
            
      end do
         
   end do
      

   if (      (Src%FieldMask(MASKID_TranslationDisp) .AND. Dest%FieldMask(MASKID_TranslationDisp)) &
        .or. (Src%FieldMask(MASKID_TranslationVel ) .AND. Dest%FieldMask(MASKID_TranslationVel )) &
        .or. (Src%FieldMask(MASKID_TranslationAcc ) .AND. Dest%FieldMask(MASKID_TranslationAcc )) ) then
               
      
         ! calculate displaced positions at operating point:                           
      DO i = 1,Dest%Nnodes
      
         DO j=1,NumNodes(ELEMENT_LINE2) ! number of nodes per line2 element
            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
         
            MeshMap%DisplacedPosition(:,i,j) =     Src%Position(:,n) +  Src%TranslationDisp(:,n)  &
                                                - Dest%Position(:,i) - Dest%TranslationDisp(:,i)  
         end do
      
      END DO   
         
         
         
         ! MeshMap%dM%fx_p required for all three transfers:
      if (.not. allocated(MeshMap%dM%fx_p) ) then
         call AllocAry(MeshMap%dM%fx_p, Dest%Nnodes*3, Src%Nnodes*3, 'dM%fx_p', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      end if
                                    
      MeshMap%dM%fx_p = 0.0_R8Ki      
      do i=1, Dest%Nnodes

         d_start = (i-1)*3+1
         d_end   = d_start+2
            
         do j=1,NumNodes(ELEMENT_LINE2) 
               
            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
               
            s_start = (n - 1)*3+1
            s_end   = s_start+2
               
            tmpVec = MeshMap%DisplacedPosition(:,i,j) * MeshMap%MapMotions(i)%shape_fn(j) 
            MeshMap%dM%fx_p( d_start:d_end, s_start:s_end ) = SkewSymMat( tmpVec ) 
               
         end do
            
      end do
      
                  
      ! MeshMap%dM%tv_uS and MeshMap%dM%tv_uD required for translational velocity:         
      if ( Src%FieldMask(MASKID_TranslationVel) .AND. Dest%FieldMask(MASKID_TranslationVel) ) then
         
         if ( Src%FieldMask(MASKID_RotationVel) ) then
            if (.not. allocated(MeshMap%dM%tv_uD) ) then
               call AllocAry(MeshMap%dM%tv_uD, Dest%Nnodes*3, Dest%Nnodes*3, 'dM%tv_uD', ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  IF (ErrStat >= AbortErrLev) RETURN
            end if
                     
            if (.not. allocated(MeshMap%dM%tv_uS) ) then
               call AllocAry(MeshMap%dM%tv_uS, Dest%Nnodes*3, Src%Nnodes*3, 'dM%tv_uS', ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  IF (ErrStat >= AbortErrLev) RETURN
            end if
         
            MeshMap%dM%tv_uD = 0.0_R8Ki
            MeshMap%dM%tv_uS = 0.0_R8Ki
            do i=1, Dest%Nnodes

               d_start = (i-1)*3+1
               d_end   = d_start+2
                  
               do j=1,NumNodes(ELEMENT_LINE2) 
               
                  n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
                  
                  s_start = (n - 1)*3+1
                  s_end   = s_start+2
                                    
                  TmpVec = Src%RotationVel(:,n) * MeshMap%MapMotions(i)%shape_fn(j)
                  SSMat = SkewSymMat( TmpVec )
                  
                  MeshMap%dM%tv_uD( d_start:d_end, d_start:d_end ) = MeshMap%dM%tv_uD( d_start:d_end, d_start:d_end ) +  SSMat
                  MeshMap%dM%tv_uS( d_start:d_end, s_start:s_end ) = -SSMat
               end do
                  
            end do
            
         else
            if (allocated(MeshMap%dM%tv_uD)) deallocate(MeshMap%dM%tv_uD)               
            if (allocated(MeshMap%dM%tv_uS)) deallocate(MeshMap%dM%tv_uS)               
         end if !MASKID_RotationVel
         
      else
         if (allocated(MeshMap%dM%tv_uD)) deallocate(MeshMap%dM%tv_uD)               
         if (allocated(MeshMap%dM%tv_uS)) deallocate(MeshMap%dM%tv_uS)               
      end if !MASKID_TranslationVel
            
         
      if ( Src%FieldMask(MASKID_TranslationAcc ) .AND. Dest%FieldMask(MASKID_TranslationAcc ) ) then
            
         !-------------- ta_uD and ta_uS -----------------------------
         if (.not. allocated(MeshMap%dM%ta_uD) ) then
            call AllocAry(MeshMap%dM%ta_uD, Dest%Nnodes*3, Dest%Nnodes*3, 'dM%ta_uD', ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) RETURN
         end if
         if (.not. allocated(MeshMap%dM%ta_uS) ) then
            call AllocAry(MeshMap%dM%ta_uS, Dest%Nnodes*3, Src%Nnodes*3, 'dM%ta_uS', ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) RETURN
         end if
         
                        
         MeshMap%dM%ta_uD = 0.0_R8Ki
         MeshMap%dM%ta_uS = 0.0_R8Ki
         if ( Src%FieldMask(MASKID_RotationAcc) ) then            
            do i=1, Dest%Nnodes
            
               d_start = (i-1)*3+1
               d_end   = d_start+2
                  
               do j=1,NumNodes(ELEMENT_LINE2) 
               
                  n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
                  
                  s_start = (n - 1)*3+1
                  s_end   = s_start+2
                                    
                  TmpVec = Src%RotationAcc(:,n) * MeshMap%MapMotions(i)%shape_fn(j) 
                  SSMat = SkewSymMat( TmpVec ) 
                  MeshMap%dM%ta_uD( d_start:d_end, d_start:d_end ) =  MeshMap%dM%ta_uD( d_start:d_end, d_start:d_end ) + SSMat
                  MeshMap%dM%ta_uS( d_start:d_end, s_start:s_end ) = -SSMat
               end do
            end do
         end if
         
         if ( Src%FieldMask(MASKID_RotationVel) ) then            
                
            do i=1, Dest%Nnodes
            
               d_start = (i-1)*3+1
               d_end   = d_start+2
                  
               do j=1,NumNodes(ELEMENT_LINE2) 
               
                  n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
                  s_start = (n - 1)*3+1
                  s_end   = s_start+2
                              
                  TmpVec = Src%RotationVel(:,n) * MeshMap%MapMotions(i)%shape_fn(j)
                  RotVel = Src%RotationVel(:,n)
                  SSMat = OuterProduct( RotVel, TmpVec ) 
                  tmp   =  dot_product( RotVel, TmpVec )
                  do k=1,3
                     SSMat(k,k) = SSMat(k,k) - tmp
                  end do                        
                                    
                  MeshMap%dM%ta_uD( d_start:d_end, d_start:d_end ) =  MeshMap%dM%ta_uD( d_start:d_end, d_start:d_end ) + SSMat
                  MeshMap%dM%ta_uS( d_start:d_end, s_start:s_end ) =  MeshMap%dM%ta_uS( d_start:d_end, s_start:s_end ) - SSMat
                                                         
               end do
                  
            end do
         end if
            

                               
            
         !-------------- ta_rv -----------------------------
         if ( Src%FieldMask(MASKID_RotationVel) ) then
            if (.not. allocated(MeshMap%dM%ta_rv) ) then
               call AllocAry(MeshMap%dM%ta_rv, Dest%Nnodes*3, Src%Nnodes*3, 'dM%ta_rv', ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  IF (ErrStat >= AbortErrLev) RETURN
            end if
                        
            MeshMap%dM%ta_rv = 0.0_R8Ki             

            do i=1, Dest%Nnodes
            
               d_start = (i-1)*3+1
               d_end   = d_start+2
                  
               do j=1,NumNodes(ELEMENT_LINE2) 
               
                  n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
                  
                  s_start = (n - 1)*3+1
                  s_end   = s_start+2
                                    
                  RotVel = Src%RotationVel(:,n)
                  tmpVec = cross_product( RotVel, MeshMap%DisplacedPosition(:,i,j)  )
                  
                  SSMat = SkewSymMat( tmpVec ) + OuterProduct( MeshMap%DisplacedPosition(:,i,j), RotVel )
                  MeshMap%dM%ta_rv( d_start:d_end, s_start:s_end ) =  SSMat * MeshMap%MapMotions(i)%shape_fn(j) 
                                       
                  tmp=dot_product( MeshMap%DisplacedPosition(:,i,1), RotVel ) * MeshMap%MapMotions(i)%shape_fn(j)
                  do k=0,2
                     MeshMap%dM%ta_rv( d_start+k, s_start+k ) = MeshMap%dM%ta_rv( d_start+k, s_start+k ) - tmp 
                  end do
                                    
               end do                  
            end do
                            
         else
            if (allocated(MeshMap%dM%ta_rv)) deallocate(MeshMap%dM%ta_rv)            
         end if ! MASKID_RotationVel            
      else
         if (allocated(MeshMap%dM%ta_uD)) deallocate(MeshMap%dM%ta_uD)
         if (allocated(MeshMap%dM%ta_uS)) deallocate(MeshMap%dM%ta_uS)            
         if (allocated(MeshMap%dM%ta_rv)) deallocate(MeshMap%dM%ta_rv)            
      end if ! MASKID_TranslationAcc
         
   end if ! MASKID_TranslationDisp, MASKID_RotationVel, or MASKID_TranslationAcc      
                  

END SUBROUTINE Linearize_Motions_Line2_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine projects Mesh1 onto a Line2 mesh (Mesh2) to find the element mappings between the two meshes.
SUBROUTINE CreateMapping_ProjectToLine2(Mesh1, Mesh2, NodeMap, Mesh1_TYPE, ErrStat, ErrMsg)

   TYPE(MeshType),                 INTENT(IN   )  :: Mesh1                          !< The mesh in the outer mapping loop (Dest for Motions/Scalars; Src for Loads)
   TYPE(MeshType),                 INTENT(IN   )  :: Mesh2                          !< The mesh in the inner mapping loop (Src for Motions/Scalars; Dest for Loads)

   TYPE(MapType),                  INTENT(INOUT)  :: NodeMap(:)                     !< The mapping from Src to Dest

   INTEGER(IntKi),                 INTENT(IN   )  :: Mesh1_TYPE                     !< Type of Mesh1 elements to map
   INTEGER(IntKi),   PARAMETER                    :: Mesh2_TYPE  = ELEMENT_LINE2    !< Type of Mesh2 elements on map (MUST BE ELEMENT_LINE2)

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                         !< Error message if ErrStat /= ErrID_None
   

      ! local variables

!   INTEGER(IntKi)                                 :: ErrStat2                       ! Error status of the operation
!   CHARACTER(ErrMsgLen)                           :: ErrMsg2                        ! Error message if ErrStat2 /= ErrID_None
#ifdef DEBUG_MESHMAPPING
   CHARACTER(200)                                 :: DebugFileName                  ! File name for debugging file
#endif   
   CHARACTER(*), PARAMETER                        :: RoutineName = 'CreateMapping_ProjectToLine2' 
   
   
   REAL(ReKi)      :: denom
   REAL(ReKi)      :: dist
   REAL(ReKi)      :: min_dist
   REAL(ReKi)      :: elem_position
   REAL(SiKi)      :: elem_position_SiKi

   REAL(ReKi)      :: Mesh1_xyz(3)

   REAL(ReKi)      :: n1_n2_vector(3)     ! vector going from node 1 to node 2 in Line2 element
   REAL(ReKi)      :: n1_Point_vector(3)  ! vector going from node 1 in Line 2 element to Destination Point
   REAL(ReKi)      :: tmp(3)              ! temporary vector for cross product calculation


   INTEGER(IntKi)  :: iElem, iNode, i  ! do-loop counter for elements on Mesh1, associated node(S)
   INTEGER(IntKi)  :: jElem            ! do-loop counter for elements on Mesh2, associated node

   INTEGER(IntKi)  :: n1, n2           ! nodes associated with an element

   LOGICAL         :: found
   LOGICAL         :: on_element
   REAL(ReKi)      :: closest_elem_position
   INTEGER(IntKi)  :: closest_elem
   REAL(ReKi)      :: closest_elem_diff
   REAL(ReKi)      :: closest_elem_distance
   
#ifdef DEBUG_MESHMAPPING
   INTEGER(IntKi)       :: Un               ! unit number for debugging
   INTEGER(IntKi)       :: ErrStat2
   CHARACTER(ErrMsgLen) :: ErrMsg2
#endif
   


      ! initialization
   ErrStat = ErrID_None
   ErrMsg  = ""


   ! Map the source nodes to destination nodes:
   do n1=1,size(NodeMap)
      NodeMap(n1)%OtherMesh_Element = NODE_NOT_MAPPED ! initialize this so we know if we've mapped this node already (done only because we may have different elements)
   end do !n1
      


   do iElem = 1, Mesh1%ElemTable(Mesh1_TYPE)%nelem   ! number of Mesh1_TYPE elements on Mesh1
      do iNode = 1, SIZE( Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes )
         i = Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes(iNode)  ! the nodes on element iElem
         IF ( NodeMap(i)%OtherMesh_Element > 0 ) CYCLE  ! we already mapped this node; let's move on to the next iNode (or iElem)

         ! destination point
         Mesh1_xyz = Mesh1%Position(:, i)

         found = .false.
         min_dist = HUGE(min_dist)
         
            ! some values for finding mapping if there are some numerical issues
         closest_elem_diff = HUGE(min_dist)
         closest_elem = 0

         do jElem = 1, Mesh2%ElemTable(Mesh2_TYPE)%nelem  ! ELEMENT_LINE2 = Mesh2_TYPE

               ! write(*,*) 'i,jElem = ', i,jElem, 'found = ', found

               ! grab node numbers associated with the jElem_th element
            n1 = Mesh2%ElemTable(Mesh2_TYPE)%Elements(jElem)%ElemNodes(1)
            n2 = Mesh2%ElemTable(Mesh2_TYPE)%Elements(jElem)%ElemNodes(2)

               ! Calculate vectors used in projection operation

            n1_n2_vector    = Mesh2%Position(:,n2) - Mesh2%Position(:,n1)
            n1_Point_vector = Mesh1_xyz - Mesh2%Position(:,n1)

            denom           = DOT_PRODUCT( n1_n2_vector, n1_n2_vector )
            IF ( EqualRealNos( denom, 0.0_ReKi ) ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Division by zero because Line2 element nodes are in same position.', ErrStat, ErrMsg, RoutineName)
               RETURN
            END IF

               ! project point onto line defined by n1 and n2

            elem_position = DOT_PRODUCT(n1_n2_vector,n1_Point_vector) / denom

                  ! note: i forumlated it this way because Fortran doesn't necessarially do shortcutting and I don't want to call EqualRealNos if we don't need it:
            if ( elem_position .ge. 0.0_ReKi .and. elem_position .le. 1.0_ReKi ) then !we're ON the element (between the two nodes)
               on_element = .true.
            else
               elem_position_SiKi = REAL( elem_position, SiKi )
               if (EqualRealNos( elem_position_SiKi, 1.0_SiKi )) then !we're ON the element (at a node)
                  on_element = .true.
                  elem_position = 1.0_ReKi
               elseif (EqualRealNos( elem_position_SiKi,  0.0_SiKi )) then !we're ON the element (at a node)
                  on_element = .true.
                  elem_position = 0.0_ReKi
               else !we're not on the element
                  on_element = .false.
                  
                  if (.not. found) then ! see if we have are very close to the end of an element (numerical roundoff?)
                     if ( elem_position_SiKi < 0.0_SiKi ) then
                        if ( -elem_position_SiKi < closest_elem_diff ) then
                           closest_elem_diff = -elem_position_SiKi
                           closest_elem = jElem
                           closest_elem_position = 0.0_ReKi
                           closest_elem_distance    = sqrt(denom) * closest_elem_diff ! distance from end of element, in meters
                        end if
                     else
                        if ( elem_position_SiKi-1.0_SiKi < closest_elem_diff ) then
                           closest_elem_diff = elem_position_SiKi-1.0_SiKi
                           closest_elem = jElem
                           closest_elem_position = 1.0_ReKi
                           closest_elem_distance    = sqrt(denom) * closest_elem_diff ! distance from end of element, in meters
                        end if
                     end if
                  end if
                  
               end if
            end if

            if (on_element) then

               ! calculate distance between point and line (note: this is actually the distance squared);
               ! will only store information once we have determined the closest element
               tmp  = cross_product( n1_n2_vector, n1_Point_vector )
               dist = DOT_PRODUCT(tmp,tmp) / denom

               if (dist .lt. min_dist) then
                  found = .true.
                  min_dist = dist

                  NodeMap(i)%OtherMesh_Element = jElem
                  NodeMap(i)%shape_fn(1)       = 1.0_ReKi - elem_position
                  NodeMap(i)%shape_fn(2)       = elem_position

                  !NodeMap(i)%couple_arm        = n1_Point_vector

               end if !the point is closest to this line2 element

            endif

         end do !jElem

            ! if failed to find an element that the Point projected into, throw an error
         if (.not. found) then
            if ( closest_elem_distance <= 7.5e-3 ) then ! if it is within 7.5mm of the end of an element, we'll accept it
               NodeMap(i)%OtherMesh_Element = closest_elem
               NodeMap(i)%shape_fn(1)       = 1.0_ReKi - closest_elem_position
               NodeMap(i)%shape_fn(2)       = closest_elem_position
               CALL SetErrStat( ErrID_Info, 'Found close value for node '//trim(num2Lstr(i))//'. ('//trim(num2lstr(closest_elem_distance))//' m)', ErrStat, ErrMsg, RoutineName)
            end if
         
            if (NodeMap(i)%OtherMesh_Element .lt. 1 )  then
               CALL SetErrStat( ErrID_Fatal, 'Node '//trim(num2Lstr(i))//' does not project onto any line2 element.' &
                           //' Closest distance is '//trim(num2lstr(closest_elem_distance))//' m.', ErrStat, ErrMsg, RoutineName)
               
#ifdef DEBUG_MESHMAPPING
                  ! output some mesh information for debugging
               CALL GetNewUnit(Un,ErrStat2,ErrMsg2)
               DebugFileName='FAST_Meshes.'//trim(num2Lstr(Un))//'.dbg'
               CALL OpenFOutFile(Un,DebugFileName,ErrStat2,ErrMsg2)
               IF (ErrStat2 >= AbortErrLev) RETURN
               
               CALL SetErrStat( ErrID_Info, 'See '//trim(DebugFileName)//' for mesh debug information.', ErrStat, ErrMsg, RoutineName)               
               WRITE( Un, '(A,I5,A,I5,A,ES15.5,A)' ) 'Element ', closest_elem, ' is closest to node ', i, &
                                                   '. It has a relative position of ', closest_elem_diff, '.'

               WRITE( Un, '(A)') '************************************************** Mesh1 ***************************************************'
               WRITE( Un, '(A)') 'Mesh1 is the destination mesh for transfer of motions/scalars; it is the source mesh for transfer of loads.'
               WRITE( Un, '(A)') '************************************************************************************************************'
               CALL MeshPrintInfo ( Un, Mesh1 )
               WRITE( Un, '(A)') '************************************************** Mesh2 ***************************************************'
               WRITE( Un, '(A)') 'Mesh2 is the source mesh for transfer of motions/scalars; it is the destination mesh for transfer of loads.'
               WRITE( Un, '(A)') '************************************************************************************************************'
               CALL MeshPrintInfo ( Un, Mesh2 )
               ! CLOSE(Un) ! by not closing this, I can ensure unique file names.
#endif               
               
               RETURN
            endif

         end if !not found on projection to element

      end do !iNode
   end do !iElem

END SUBROUTINE CreateMapping_ProjectToLine2
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine creates a new mesh with the same positions as the Src mesh, except all of the elements are points. It adds fields for
!! forces, moments, and/or TranslationDisp, if they are part of the Src mesh.
SUBROUTINE Create_PointMesh(Src, Temp_Point_Src, ErrStat, ErrMsg)

   TYPE(MeshType),                 INTENT(IN   )  :: Src                !< The source mesh
   TYPE(MeshType),                 INTENT(INOUT)  :: Temp_Point_Src     !< A blank mesh to be created

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat            !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg             !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                 :: i !loop over the nodes
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   CALL MeshDestroy( Temp_Point_Src, ErrStat2, ErrMsg2, .TRUE. )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Create_PointMesh')
      IF (ErrStat >= AbortErrLev) RETURN


   CALL MeshCreate(   BlankMesh       = Temp_Point_Src               &
                     ,IOS             = Src%IOS                      &
                     ,NNodes          = Src%nnodes                   &
                     ,Force           = Src%FieldMask(maskid_force)  &
                     ,Moment          = Src%FieldMask(maskid_moment) &
                     ,TranslationDisp = Src%FieldMask(maskid_TranslationDisp)  & !                     ,Orientation     = Src%FieldMask(maskid_Orientation) &
                     ,ErrStat         = ErrStat2                     &
                     ,ErrMess         = ErrMsg2                 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Create_PointMesh')
      IF (ErrStat >= AbortErrLev) RETURN

   do i = 1, src%nnodes

      CALL MeshConstructElement ( Mesh = Temp_Point_Src         &
                                 ,Xelement = ELEMENT_POINT      &
                                 ,P1       = I                  &
                                 ,ErrStat  = ErrStat2           &
                                 ,ErrMess  = ErrMsg2            )

         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Create_PointMesh')
         IF (ErrStat >= AbortErrLev) RETURN
      
      CALL MeshPositionNode ( Mesh = Temp_Point_Src               &
                              ,INode = i                          &
                              ,Pos = Src%Position(:,i)            &
                              ,Orient = Src%RefOrientation(:,:,i) &
                              ,ErrStat   = ErrStat2               &
                              ,ErrMess   = ErrMsg2                )
      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Create_PointMesh')
         IF (ErrStat >= AbortErrLev) RETURN

   enddo

   CALL MeshCommit ( Temp_Point_Src, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Create_PointMesh')
      IF (ErrStat >= AbortErrLev) RETURN

END SUBROUTINE Create_PointMesh
!----------------------------------------------------------------------------------------------------------------------------------
!> routine that creats a map of line2 loads to points
SUBROUTINE CreateLoadMap_L2_to_P( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                            !< The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                           !< The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap                        !< mapping structure

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                         !< Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   ! augment the source mesh with nodes formed by projections from the destination mesh nodes
   CALL Create_Augmented_Ln2_Src_Mesh(Src, Dest, MeshMap, ELEMENT_POINT, ErrStat2, ErrMsg2)         
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateLoadMap_L2_to_P')
      IF (ErrStat >= AbortErrLev) RETURN
                  
   ! Create a temporary mesh for lumped point elements of the line2 mesh
   CALL Create_PointMesh( MeshMap%Augmented_Ln2_Src, MeshMap%Lumped_Points_Src, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateLoadMap_L2_to_P')
      IF (ErrStat >= AbortErrLev) RETURN

   ! in following call, Src is mesh to loop over, finding a corresponding point for each point in Dest
   CALL CreateMapping_NearestNeighbor( MeshMap%Lumped_Points_Src, Dest, MeshMap%MapLoads, ELEMENT_POINT, ELEMENT_POINT, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateLoadMap_L2_to_P')
      IF (ErrStat >= AbortErrLev) RETURN
            
END SUBROUTINE CreateLoadMap_L2_to_P
!----------------------------------------------------------------------------------------------------------------------------------
!> routine that creats a map of line2 motions to points
SUBROUTINE CreateMotionMap_L2_to_P( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             !< The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            !< The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap                         !< mapping data structure

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          !< Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   CALL CreateMapping_ProjectToLine2(Dest,Src, MeshMap%MapMotions, ELEMENT_POINT, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateMotionMap_L2_to_P')
      !IF (ErrStat >= AbortErrLev) RETURN
         
END SUBROUTINE CreateMotionMap_L2_to_P        
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that transfers data from a point mesh to a line2 mesh.
SUBROUTINE Transfer_Point_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),         INTENT(IN   ) :: Src       !< source (point) mesh
   TYPE(MeshType),         INTENT(INOUT) :: Dest      !< destination (line2) mesh
   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap   !< mapping data structure

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat   !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg    !< Error message if ErrStat /= ErrID_None

   TYPE(MeshType),OPTIONAL,INTENT(IN   ) :: SrcDisp   !< a "functional" sibling of the source mesh for load mapping; Src contains loads and SrcDisp contains TranslationDisp and Orientation
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) :: DestDisp  !< a "functional" sibling of the destination mesh for load mapping; Dest contains loads and DestDisp contains TranslationDisp and Orientation

   ! local variables

   REAL(ReKi)                            :: LoadsScaleFactor  ! bjj: added this scaling factor to get loads in a better numerical range 
   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*), PARAMETER               :: RoutineName= 'Transfer_Point_to_Line2'

   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

   !.................
   ! Check to ensure that the source mesh is composed of Point elements and destination mesh is composed of Line2 elements
   !.................   
   
   if (Src%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Source mesh must have one or more Point elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif
   
   if (Dest%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Destination mesh must have one or more Line2 elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif

   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) then

      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMotionMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
         
      end if

      !........................
      ! Transfer data
      !........................

         ! This is the same algorithm as Transfer_Point_to_Point
      CALL Transfer_Motions_Point_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN


   end if ! algorithm for motions/scalars

   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if ( HasLoadFields(Src) ) then

      IF (.not. PRESENT(SrcDisp) .OR. .NOT. PRESENT(DestDisp) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'SrcDisp and DestDisp arguments are required for load transfer.', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF

      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then
                           
         CALL CreateLoadMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
                           
      end if


      !........................
      ! Transfer data
      !........................
      LoadsScaleFactor = GetLoadsScaleFactor ( Src ) 

      CALL Transfer_Loads_Point_to_Line2( Src, Dest, MeshMap, ErrStat2, ErrMsg2, SrcDisp, DestDisp, LoadsScaleFactor )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
      
     
   end if ! algorithm for loads


END SUBROUTINE Transfer_Point_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that creates linearization matricies for data tranfer from a point mesh to a line2 mesh.
!! \copydetails modmesh_mapping::linearize_point_to_point  
SUBROUTINE Linearize_Point_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),            INTENT(IN   ) :: Src       ! source (point) mesh
   TYPE(MeshType),            INTENT(IN   ) :: Dest      ! destination (line2) mesh
   TYPE(MeshMapType),         INTENT(INOUT) :: MeshMap   ! mapping data structure of Src to Dest
                              
   INTEGER(IntKi),            INTENT(  OUT) :: ErrStat   ! Error status of the operation
   CHARACTER(*),              INTENT(  OUT) :: ErrMsg    ! Error message if ErrStat /= ErrID_None
                              
   TYPE(MeshType),OPTIONAL,   INTENT(IN   ) :: SrcDisp   ! a "functional" sibling of the source mesh for load mapping; Src contains loads and SrcDisp contains TranslationDisp and Orientation
   TYPE(MeshType),OPTIONAL,   INTENT(IN   ) :: DestDisp  ! a "functional" sibling of the destination mesh for load mapping; Dest contains loads and DestDisp contains TranslationDisp and Orientation

   ! local variables

   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*), PARAMETER               :: RoutineName= 'Linearize_Point_to_Line2'

   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''


   ! ------------------------------------------------------------------------------------------------------------------------------
   !> ### Mapping and Linearization of Data Transfer for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) then

      !........................
      !> * Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMotionMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
                                 
      end if

      !........................
      !> * Get linearization matrices for data transfer
      !........................

         ! This is the same algorithm as Transfer_Point_to_Point
      CALL Linearize_Motions_Point_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN


   end if ! algorithm for motions/scalars

   ! ------------------------------------------------------------------------------------------------------------------------------
   !. ### Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if ( HasLoadFields(Src) ) then

      IF (.not. PRESENT(SrcDisp) .OR. .NOT. PRESENT(DestDisp) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'SrcDisp and DestDisp arguments are required for load transfer linearization.', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF

      !........................
      !> * Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then
                           
         CALL CreateLoadMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
                           
      end if
         
         !>  + We need a motion mapping structure from the SrcDisp to DestDisp meshes so that we can form the proper matrices for moment
         !!  linearizaiton. (This is because we have both the source translational displacement and rotational displacement on the right side
         !! of the equation; alternatively we'd have to have the source AND destination translational displacements on the right hand side.)
      if ( .not. allocated(MeshMap%MapMotions)) then
         
                  ! Allocate the mapping structure:
         ALLOCATE( MeshMap%MapMotions(DestDisp%NNodes), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error trying to allocate MeshMap%MapMotions.', ErrStat, ErrMsg, RoutineName)
            return
         ELSE            
               ! set up the initial mappings so that we don't necessarially have to do this multiple times on the first time step (if calculating Jacobians)                  
            !note this is for SrcDisp and DestDisp meshes, but they are "functional siblings" of Src and Dest (i.e., may have different element types than point and line2)
            CALL CreateMotionMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) RETURN
         END IF
         
      elseif (SrcDisp%RemapFlag .or. DestDisp%RemapFlag ) then
                           
         !note this is for SrcDisp and DestDisp meshes, but they are "functional siblings" of Src and Dest (i.e., may have different element types than line2)
         CALL CreateMotionMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
                           
      end if

      !........................
      !> * Get linearization matrices for data transfer
      !........................

      CALL Linearize_Loads_Point_to_Line2( Src, Dest, MeshMap, ErrStat2, ErrMsg2, SrcDisp, DestDisp )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
      
     
   end if ! algorithm for loads


END SUBROUTINE Linearize_Point_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CreateLoadMap_P_to_L2( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            ! The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

      ! LOCAL variables:
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !........................
   ! Create matrix used to "unlump" loads later (needs to have element connectivity information to create it)
   !........................
   CALL Create_InverseLumping_Matrix( Dest, MeshMap, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateLoadMap_P_to_L2')
         
         
   CALL CreateMapping_ProjectToLine2(Src, Dest, MeshMap%MapLoads, ELEMENT_POINT, ErrStat2, ErrMsg2)            
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateLoadMap_P_to_L2')
      IF (ErrStat >= AbortErrLev) RETURN
   
      
END SUBROUTINE CreateLoadMap_P_to_L2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CreateMotionMap_P_to_L2( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            ! The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                 :: i
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   ! Each destination node (on a LINE2 mesh) needs a source
   ! in following call, Dest is mesh to looped over, finding a corresponding point for each point in Src
   CALL CreateMapping_NearestNeighbor( Dest, Src, MeshMap%MapMotions, ELEMENT_LINE2, ELEMENT_POINT, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateMotionMap_P_to_L2')
      IF (ErrStat >= AbortErrLev) RETURN

   ! bjj: for consistant definition of couple_arm (i.e. p_ODR-p_OSR), let's multiply by -1
   do i=1,SIZE(MeshMap%MapMotions)
      MeshMap%MapMotions(i)%couple_arm = -1._ReKi*MeshMap%MapMotions(i)%couple_arm
   end do
   
         
END SUBROUTINE CreateMotionMap_P_to_L2        


!----------------------------------------------------------------------------------------------------------------------------------
!> Given two point meshes, this routine transfers the source mesh values to the destination mesh using appropriate math.
SUBROUTINE Transfer_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),          INTENT(IN   ) ::  Src      !< source point mesh
   TYPE(MeshType),          INTENT(INOUT) ::  Dest     !< destination point mesh
   TYPE(MeshType), OPTIONAL,INTENT(IN   ) ::  SrcDisp  !< an optional mesh, which contains the displacements associated with the source if the source contains load information
   TYPE(MeshType), OPTIONAL,INTENT(IN   ) ::  DestDisp !< an optional mesh, which contains the displacements associated with the destination if the destination contains load information

   TYPE(MeshMapType),       INTENT(INOUT) ::  MeshMap  !< Mapping(s) between Src and Dest meshes

   INTEGER(IntKi),          INTENT(  OUT) :: ErrStat   !< Error status of the operation
   CHARACTER(*),            INTENT(  OUT) :: ErrMsg    !< Error message if ErrStat /= ErrID_None


   ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'Transfer_Point_to_Point'
   
   REAL(ReKi)                              :: LoadsScaleFactor  ! bjj: added this scaling factor to get loads in a better numerical range 



   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

   !.................
   ! Check to ensure that both the source and destination meshes are composed of Point elements
   !.................

   if (Src%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Source mesh must have one or more Point elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif

   if (Dest%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Destination mesh must have one or more Point elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif


   !> ------------------------------------------------------------------------------------------------------------------------------
   !! Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   !! ------------------------------------------------------------------------------------------------------------------------------
   !! If Src is displacements/velocities then loop over the Destination Mesh;
   !!    each motion/scalar in the destination mesh needs to be interpolated from the source mesh.


   if ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) then

      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMotionMap_P_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
            
      end if

      !........................
      ! Transfer data
      !........................

      CALL Transfer_Motions_Point_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN

   end if ! algorithm for motions/scalars

   !> ------------------------------------------------------------------------------------------------------------------------------
   !! Mapping and Transfer of Data for Mesh Load Fields
   !! ------------------------------------------------------------------------------------------------------------------------------
   !! IF Src is forces and/or moments, loop over Src Mesh;
   !!    each load in the source mesh needs to be placed somewhere in the destination mesh.

   if ( HasLoadFields(Src) ) then
            
      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then
         
         CALL CreateLoadMap_P_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Point_to_Point:CreateLoadMap_P_to_P')
            IF (ErrStat >= AbortErrLev) RETURN                        

      end if

      !........................
      ! Transfer data
      !........................

      IF ( PRESENT( SrcDisp ) .AND. PRESENT( DestDisp ) ) THEN
                  
         LoadsScaleFactor = GetLoadsScaleFactor ( Src )
         
         CALL Transfer_Loads_Point_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2, SrcDisp, DestDisp, LoadsScaleFactor )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN                        
            
      ELSE !bjj: we could check if Src also had TranslationDisp fields and call with SrcDisp=Src
         CALL SetErrStat( ErrID_Fatal, 'Invalid arguments to Transfer_Point_to_Point for meshes with load fields.', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF

   end if ! algorithm for loads



END SUBROUTINE Transfer_Point_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine creates the node-to-node (nearest neighbor). We map FROM Mesh1 to Mesh2
SUBROUTINE CreateMapping_NearestNeighbor( Mesh1, Mesh2, NodeMap, Mesh1_TYPE, Mesh2_TYPE, ErrStat, ErrMsg )
!.......................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Mesh1      !< The mesh in the outer mapping loop (Dest for Motions/Scalars; Src for Loads)
   TYPE(MeshType),                 INTENT(IN   )  :: Mesh2      !< The mesh in the inner mapping loop (Src for Motions/Scalars; Dest for Loads)

   TYPE(MapType),                  INTENT(INOUT)  :: NodeMap(:) !< The mapping from Src to Dest

   INTEGER(IntKi),                 INTENT(IN   )  :: Mesh1_TYPE !< Type of Mesh1 elements to map
   INTEGER(IntKi),                 INTENT(IN   )  :: Mesh2_TYPE !< Type of Mesh2 elements on map

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat    !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg     !< Error message if ErrStat /= ErrID_None

      ! local variables

   REAL(ReKi)      :: dist
   REAL(ReKi)      :: min_dist

   REAL(ReKi)      :: Mesh1_xyz(3)
   REAL(ReKi)      :: Mesh2_xyz(3)

   INTEGER(IntKi)  :: point_with_min_dist
   INTEGER(IntKi)  :: iElem, iNode, i  ! do-loop counter for elements on Mesh1, associated node(S)
   INTEGER(IntKi)  :: jElem, jNode, j  ! do-loop counter for elements on Mesh2, associated node

   LOGICAL         :: UseMesh2Node(Mesh2%NNodes) ! determines if the node on the second mesh is part of the mapping (i.e., contained in an element of the appropriate type)
     
   
      ! initialization
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Determine which nodes on mesh2 are going to be in the mapping
   UseMesh2Node = .FALSE.
   do jElem = 1, Mesh2%ElemTable(Mesh2_TYPE)%nelem  ! number of point elements on Mesh2
      do jNode = 1, SIZE( Mesh2%ElemTable(Mesh2_TYPE)%Elements(jElem)%ElemNodes )
         UseMesh2Node( Mesh2%ElemTable(Mesh2_TYPE)%Elements(jElem)%ElemNodes(jNode) ) = .TRUE.
      end do
   end do

   ! Map the source nodes to destination nodes:
   do i=1,size(NodeMap)
      NodeMap(i)%OtherMesh_Element = NODE_NOT_MAPPED ! initialize this so we know if we've mapped this node already (done only because we may have different elements)
   end do !n1
   

   do iElem = 1, Mesh1%ElemTable(Mesh1_TYPE)%nelem   ! number of Mesh1_TYPE elements on Mesh1 = number of points on Mesh1
      do iNode = 1, SIZE( Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes )
         i = Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes(iNode)  ! the nodes on element iElem
         IF ( NodeMap(i)%OtherMesh_Element > 0 ) CYCLE  ! we already mapped this node; let's move on


         ! Find the nearest neighbor node for this particular node

         ! initialize minimum distance marker at some huge number
         min_dist = HUGE(min_dist)
         point_with_min_dist = 0

         Mesh1_xyz = Mesh1%Position(:, i)

         do j = 1, Mesh2%NNodes
            IF ( .NOT. UseMesh2Node(j) ) CYCLE !This node isn't part of the elements we're mapping

            ! destination point
            Mesh2_xyz = Mesh2%Position(:, j)

            ! calculate distance  between source and desination; will only store information once we have determined
            ! the closest point
            dist = sqrt(  (Mesh1_xyz(1) - Mesh2_xyz(1))**2 &
                        + (Mesh1_xyz(2) - Mesh2_xyz(2))**2 &
                        + (Mesh1_xyz(3) - Mesh2_xyz(3))**2 )

            if (dist .lt. min_dist) then

               min_dist = dist
               point_with_min_dist = j

               !if (EqualRealNos(dist), 0.0_ReKi)) EXIT !we have an exact match so let's just stop looking

            endif

         end do !j

         if (point_with_min_dist .lt. 1 )  then
            CALL SetErrStat( ErrID_Fatal, 'Failed to find destination point associated with source point.', ErrStat, ErrMsg, 'CreateMapping_NearestNeighbor')
            RETURN
         endif

         NodeMap(i)%OtherMesh_Element = point_with_min_dist !bjj: For consistency, I really wish we had used element numbers here instead....

         NodeMap(i)%distance = min_dist

         NodeMap(i)%couple_arm = Mesh2%Position(:, point_with_min_dist) - Mesh1_xyz
         !bjj: this is the negative of the case where it's Mesh2=src, so we'll have to multiply by -1 outside this routine if that's the case

      end do !iNode
   end do !iElem


END SUBROUTINE CreateMapping_NearestNeighbor
!----------------------------------------------------------------------------------------------------------------------------------
!> Given a nearest-neighbor mapping, this routine transfers motions between nodes on the mesh.
SUBROUTINE Transfer_Motions_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src       !< The source mesh with motion fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest      !< The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   !< data for the mesh mapping

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)  :: nScalars
   INTEGER(IntKi)  :: i, j                                     ! counter over the nodes
   REAL(R8Ki)      :: RotationMatrix(3,3)
   REAL(ReKi)      :: TmpVec(3)


   ErrStat = ErrID_None
   ErrMsg  = ""



      ! ---------------------------- Translation ------------------------------------------------
      !> Translational Displacement: \f$\vec{u}^D = \vec{u}^S + \left[\left[\theta^S\right]^T \theta^{SR} - I\right]\left\{\vec{p}^{ODR}-\vec{p}^{OSR}\right\}\f$

      ! u_Dest = u_Src + [Orientation_Src^T * RefOrientation_Src - I] * [p_Dest - p_Src]
   if ( Src%FieldMask(MASKID_TranslationDisp) .AND. Dest%FieldMask(MASKID_TranslationDisp) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%TranslationDisp(:,i) = Src%TranslationDisp(:,MeshMap%MapMotions(i)%OtherMesh_Element)

            ! if Src mesh has orientation, superpose Dest displacement with translation due to rotation and couple arm
         if ( Src%FieldMask(MASKID_Orientation) ) then

               !Calculate RotationMatrix as O_S^T*O_SR
            RotationMatrix = TRANSPOSE( Src%Orientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element) )
            RotationMatrix = MATMUL( RotationMatrix, Src%RefOrientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element) )

               ! subtract I
            do j=1,3
               RotationMatrix(j,j)= RotationMatrix(j,j) - 1.0_R8Ki
            end do


            Dest%TranslationDisp(:,i) = Dest%TranslationDisp(:,i) + MATMUL(RotationMatrix, MeshMap%MapMotions(i)%couple_arm)

         end if

      end do

   end if


      ! ---------------------------- ORIENTATION/Direction Cosine Matrix   ----------------------
      !> Orientation: \f$\theta^D = \theta^{DR}\left[\theta^{SR}\right]^T\theta^S\f$

      ! transfer direction cosine matrix, aka orientation

   if ( Src%FieldMask(MASKID_Orientation) .AND. Dest%FieldMask(MASKID_Orientation) ) then

      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
         
         RotationMatrix = TRANSPOSE( Src%RefOrientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element) )
         RotationMatrix = MATMUL( Dest%RefOrientation(:,:,i), RotationMatrix )
         Dest%Orientation(:,:,i) = MATMUL( RotationMatrix, Src%Orientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element) )
      end do

   endif

      ! ---------------------------- Calculated total displaced positions  ---------------------
      ! these values are used in both the translational velocity and translational acceleration
      ! calculations. The calculations rely on the TranslationDisp fields, which are calculated
      ! earlier in this routine.
   IF ( Src%FieldMask(MASKID_TranslationVel) .OR. Src%FieldMask(MASKID_TranslationAcc) ) THEN

      DO i = 1,Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
         MeshMap%DisplacedPosition(:,i,1) =    Src%TranslationDisp(:,MeshMap%MapMotions(i)%OtherMesh_Element) &
                                            - Dest%TranslationDisp(:,i) &
                                            - MeshMap%MapMotions(i)%couple_arm
      END DO
      
   END IF
   
      ! ---------------------------- TranslationVel  --------------------------------------------
      !> Translational Velocity: \f$\vec{v}^D = \vec{v}^S 
      !!              + \left\{ \left\{ \vec{p}^{OSR} + \vec{u}^S \right\} - \left\{ \vec{p}^{ODR} + \vec{u}^D \right\} \right\} \times \vec{\omega}^S\f$

   if ( Src%FieldMask(MASKID_TranslationVel) .AND. Dest%FieldMask(MASKID_TranslationVel) ) then
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%TranslationVel(:,i) = Src%TranslationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element)

         if ( Src%FieldMask(MASKID_RotationVel) ) then
            Dest%TranslationVel(:,i) = Dest%TranslationVel(:,i) + &
                                       cross_product ( MeshMap%DisplacedPosition(:,i,1), &
                                                       Src%RotationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element) )
         endif
      end do

   endif

      ! ---------------------------- RotationVel  -----------------------------------------------
      !> Rotational Velocity: \f$\vec{\omega}^D = \vec{\omega}^S\f$

   if ( Src%FieldMask(MASKID_RotationVel) .AND. Dest%FieldMask(MASKID_RotationVel) ) then
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%RotationVel(:,i) = Src%RotationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element)
      end do
   end if

      ! ---------------------------- TranslationAcc -----------------------------------------------
      !> Translational Acceleration: \f$\vec{a}^D = \vec{a}^S 
      !!                        + \left\{ \left\{ \vec{p}^{OSR} + \vec{u}^S \right\} - \left\{ \vec{p}^{ODR} + \vec{u}^D \right\} \right\} \times \vec{\alpha}^S
      !!                        + \vec{\omega}^S \times \left\{
      !!                          \left\{ \left\{ \vec{p}^{OSR} + \vec{u}^S \right\} - \left\{ \vec{p}^{ODR} + \vec{u}^D \right\} \right\} \times \vec{\omega}^S
      !!                          \right\}
      !!\f$

   if ( Src%FieldMask(MASKID_TranslationAcc) .AND. Dest%FieldMask(MASKID_TranslationAcc) ) then
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%TranslationAcc(:,i) = Src%TranslationAcc(:,MeshMap%MapMotions(i)%OtherMesh_Element)

         if ( Src%FieldMask(MASKID_RotationAcc) )  then
            Dest%TranslationAcc(:,i) = Dest%TranslationAcc(:,i) + &
                                       cross_product( MeshMap%DisplacedPosition(:,i,1), &
                                                      Src%RotationAcc(:,MeshMap%MapMotions(i)%OtherMesh_Element) )
         endif

         if ( Src%FieldMask(MASKID_RotationVel) )  then
            TmpVec = cross_product( MeshMap%DisplacedPosition(:,i,1), Src%RotationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element) )

            Dest%TranslationAcc(:,i) = Dest%TranslationAcc(:,i) + &
                                       cross_product( Src%RotationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element), TmpVec)
         endif
      end do
   endif

      ! ---------------------------- RotationAcc  -----------------------------------------------
      !> Rotational Acceleration: \f$\vec{\alpha}^D = \vec{\alpha}^S\f$

   if (Src%FieldMask(MASKID_RotationAcc) .AND. Dest%FieldMask(MASKID_RotationAcc) ) then
      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%RotationAcc(:,i) = Src%RotationAcc(:,MeshMap%MapMotions(i)%OtherMesh_Element)
      end do
   end if

      ! ---------------------------- Scalars  -----------------------------------------------
      !> Scalars: \f$S^D = S^S\f$

   if (Src%FieldMask(MASKID_SCALAR) .AND. Dest%FieldMask(MASKID_SCALAR) ) then
      nScalars = min(Dest%nScalars, Src%nScalars)
      
      if (Dest%nScalars > nScalars) then
         call SetErrStat(ErrID_Severe, "Not all scalars could be computed from source mesh (insufficient data).", ErrStat, ErrMsg, 'Transfer_Motions_Point_to_Point')
         Dest%Scalars(nScalars+1:,:) = 0.0_ReKi
      end if

      do i=1, Dest%Nnodes
         !if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%Scalars(1:nScalars,i) = Src%Scalars(1:nScalars,MeshMap%MapMotions(i)%OtherMesh_Element)
      end do
   end if


END SUBROUTINE Transfer_Motions_Point_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
!> Given a nearest-neighbor mapping, this routine forms the linearization matrices of motion transfer between nodes on the mesh.
SUBROUTINE Linearize_Motions_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src       !< The source mesh with motion fields allocated
   TYPE(MeshType),                 INTENT(IN   )  :: Dest      !< The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   !< data for the mesh mapping

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'Linearize_Motions_Point_to_Point'
   
   integer(intKi)                          :: i,j,k, n, d_start, d_end, s_start, s_end
   real(r8Ki)                              :: tmp, tmpVec(3), SSMat(3,3)
   real(r8Ki)                              :: RotVel(3), RotAcc(3)
   
   

   ErrStat = ErrID_None
   ErrMsg  = ""


   
      if (.not. allocated(MeshMap%dM%mi) ) then
         call AllocAry(MeshMap%dM%mi, Dest%Nnodes*3, Src%Nnodes*3, 'dM%mi', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      end if
      
      MeshMap%dM%mi = 0.0_R8Ki      
      do i=1, Dest%Nnodes

         n = MeshMap%MapMotions(i)%OtherMesh_Element
         do j=1,3
            MeshMap%dM%mi( (i-1)*3+j, (n-1)*3+j ) = 1.0_R8Ki
         end do         
      end do
      

      if (      (Src%FieldMask(MASKID_TranslationDisp) .AND. Dest%FieldMask(MASKID_TranslationDisp)) &
           .or. (Src%FieldMask(MASKID_TranslationVel ) .AND. Dest%FieldMask(MASKID_TranslationVel )) &
           .or. (Src%FieldMask(MASKID_TranslationAcc ) .AND. Dest%FieldMask(MASKID_TranslationAcc )) ) then
               
      
            ! calculate displaced positions at operating point:
         DO i = 1,Dest%Nnodes
            MeshMap%DisplacedPosition(:,i,1) =    Src%TranslationDisp(:,MeshMap%MapMotions(i)%OtherMesh_Element) &
                                               - Dest%TranslationDisp(:,i) &
                                               - MeshMap%MapMotions(i)%couple_arm
         END DO
         
            ! MeshMap%dM%fx_p required for all three transfers:
         if (.not. allocated(MeshMap%dM%fx_p) ) then
            call AllocAry(MeshMap%dM%fx_p, Dest%Nnodes*3, Src%Nnodes*3, 'dM%fx_p', ErrStat2, ErrMsg2 )
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) RETURN
         end if
                                    
         MeshMap%dM%fx_p = 0.0_R8Ki      
         do i=1, Dest%Nnodes

            d_start = (i-1)*3+1
            d_end   = d_start+2
            s_start = (MeshMap%MapMotions(i)%OtherMesh_Element - 1)*3+1
            s_end   = s_start+2
            MeshMap%dM%fx_p( d_start:d_end, s_start:s_end ) = SkewSymMat( MeshMap%DisplacedPosition(:,i,1) )
         end do
      
                  
            ! MeshMap%dM%tv_uS and MeshMap%dM%tv_uD required for translational velocity:         
         if ( Src%FieldMask(MASKID_TranslationVel) .AND. Dest%FieldMask(MASKID_TranslationVel) ) then
                        
            if ( Src%FieldMask(MASKID_RotationVel) ) then
               
               if (.not. allocated(MeshMap%dM%tv_uD) ) then
                  call AllocAry(MeshMap%dM%tv_uD, Dest%Nnodes*3, Dest%Nnodes*3, 'dM%tv_uD', ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                     IF (ErrStat >= AbortErrLev) RETURN
               end if
               if (.not. allocated(MeshMap%dM%tv_uS) ) then
                  call AllocAry(MeshMap%dM%tv_uS, Dest%Nnodes*3, Src%Nnodes*3, 'dM%tv_uS', ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                     IF (ErrStat >= AbortErrLev) RETURN
               end if
                     
               MeshMap%dM%tv_uD = 0.0_R8Ki      
               MeshMap%dM%tv_uS = 0.0_R8Ki      
               do i=1, Dest%Nnodes

                  n = MeshMap%MapMotions(i)%OtherMesh_Element
                  d_start = (i-1)*3+1
                  d_end   = d_start+2
                  
                  s_start = (n - 1)*3+1
                  s_end   = s_start+2
                  
                  RotVel = Src%RotationVel(:,n)
                  SSMat = SkewSymMat( RotVel )
                  MeshMap%dM%tv_uD( d_start:d_end, d_start:d_end ) = MeshMap%dM%tv_uD( d_start:d_end, d_start:d_end ) + SSMat
                  MeshMap%dM%tv_uS( d_start:d_end, s_start:s_end ) = -SSMat
                  
               end do
            else
               if (allocated(MeshMap%dM%tv_uD)) deallocate(MeshMap%dM%tv_uD)               
               if (allocated(MeshMap%dM%tv_uS)) deallocate(MeshMap%dM%tv_uS)               
            end if !MASKID_RotationVel
         else
            if (allocated(MeshMap%dM%tv_uD)) deallocate(MeshMap%dM%tv_uD)               
            if (allocated(MeshMap%dM%tv_uS)) deallocate(MeshMap%dM%tv_uS)               
         end if !MASKID_TranslationVel
            
         
         if ( Src%FieldMask(MASKID_TranslationAcc ) .AND. Dest%FieldMask(MASKID_TranslationAcc ) ) then
            
            !-------------- ta_uD and ta_uS -----------------------------
            ! these are the matrix that relates destination and source displacements to translational acceleration
            
            if (.not. allocated(MeshMap%dM%ta_uD) ) then
               call AllocAry(MeshMap%dM%ta_uD, Dest%Nnodes*3, Dest%Nnodes*3, 'dM%ta_uD', ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  IF (ErrStat >= AbortErrLev) RETURN
            end if
            if (.not. allocated(MeshMap%dM%ta_uS) ) then
               call AllocAry(MeshMap%dM%ta_uS, Dest%Nnodes*3, Src%Nnodes*3, 'dM%ta_uS', ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  IF (ErrStat >= AbortErrLev) RETURN
            end if                     
                        
            MeshMap%dM%ta_uD = 0.0_R8Ki      
            MeshMap%dM%ta_uS = 0.0_R8Ki      
            if ( Src%FieldMask(MASKID_RotationAcc) ) then            
               do i=1, Dest%Nnodes
            
                  n = MeshMap%MapMotions(i)%OtherMesh_Element
                  d_start = (i-1)*3+1
                  d_end   = d_start+2
                  s_start = (n-1)*3+1
                  s_end   = s_start+2
                                    
                  RotAcc = Src%RotationAcc(:,n)
                  SSMat = SkewSymMat( RotAcc )
                  MeshMap%dM%ta_uD( d_start:d_end, d_start:d_end ) =  MeshMap%dM%ta_uD( d_start:d_end, d_start:d_end ) + SSMat
                  MeshMap%dM%ta_uS( d_start:d_end, s_start:s_end ) = -SSMat
                  
               end do
            end if
            
            if ( Src%FieldMask(MASKID_RotationVel) ) then            
            
               do i=1, Dest%Nnodes
            
                  n = MeshMap%MapMotions(i)%OtherMesh_Element
                  d_start = (i-1)*3+1
                  d_end   = d_start+2
                  s_start = (n-1)*3+1
                  s_end   = s_start+2
                                    
                  RotVel = Src%RotationVel(:,n)
                  SSMat = OuterProduct( RotVel, RotVel )
                  tmp   =  dot_product( RotVel, RotVel )
                  do k=1,3
                     SSMat(k,k) = SSMat(k,k) - tmp
                  end do                        
                                    
                  MeshMap%dM%ta_uD( d_start:d_end, d_start:d_end ) =  MeshMap%dM%ta_uD( d_start:d_end, d_start:d_end ) + SSMat
                  MeshMap%dM%ta_uS( d_start:d_end, s_start:s_end ) =  MeshMap%dM%ta_uS( d_start:d_end, s_start:s_end ) - SSMat
               end do
            end if
                               
            
            !-------------- ta_rv -----------------------------
            ! this is the matrix that relates rotational velocity to translational acceleration
            if ( Src%FieldMask(MASKID_RotationVel) ) then
            
               if (.not. allocated(MeshMap%dM%ta_rv) ) then
                  call AllocAry(MeshMap%dM%ta_rv, Dest%Nnodes*3, Src%Nnodes*3, 'dM%ta_rv', ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                     IF (ErrStat >= AbortErrLev) RETURN
               end if
                        
               MeshMap%dM%ta_rv = 0.0_R8Ki
               do i=1, Dest%Nnodes
            
                  n = MeshMap%MapMotions(i)%OtherMesh_Element
                  d_start = (i-1)*3+1
                  d_end   = d_start+2
                  s_start = (n - 1)*3+1
                  s_end   = s_start+2
                                    
                  RotVel = Src%RotationVel(:,n)
                  tmpVec = cross_product( RotVel, MeshMap%DisplacedPosition(:,i,1)  )
                  
                  MeshMap%dM%ta_rv( d_start:d_end, s_start:s_end ) = SkewSymMat( tmpVec  ) + &
                                                                     OuterProduct( MeshMap%DisplacedPosition(:,i,1), RotVel )
                  
                  tmp=dot_product( MeshMap%DisplacedPosition(:,i,1), RotVel )
                  do j=0,2
                     MeshMap%dM%ta_rv( d_start+j, s_start+j ) = MeshMap%dM%ta_rv( d_start+j, s_start+j ) - tmp
                  end do                        
                  
               end do
            else
               if (allocated(MeshMap%dM%ta_rv)) deallocate(MeshMap%dM%ta_rv)
            end if ! MASKID_RotationVel  
            
         else
            if (allocated(MeshMap%dM%ta_uD)) deallocate(MeshMap%dM%ta_uD)
            if (allocated(MeshMap%dM%ta_uS)) deallocate(MeshMap%dM%ta_uS)            
            if (allocated(MeshMap%dM%ta_rv)) deallocate(MeshMap%dM%ta_rv)            
         end if ! MASKID_TranslationAcc
         
      end if ! MASKID_TranslationDisp, MASKID_RotationVel, or MASKID_TranslationAcc   
   

END SUBROUTINE Linearize_Motions_Point_to_Point
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
!> Given a nearest-neighbor mapping, this routine transfers loads between point nodes on the meshes.
SUBROUTINE Transfer_Loads_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp,LoadsScaleFactor )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                !< The source mesh with loads fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest               !< The destination mesh
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp            !< A mesh that contains the displacements associated with the source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp           !< A mesh that contains the displacements associated with the destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap            !< The mapping data structure (from Dest to Src)
   REAL(ReKi),                     INTENT(IN)     :: LoadsScaleFactor   !< Scaling factor for loads (to help with numerical issues)

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat            !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg             !< Error message if ErrStat /= ErrID_None

      ! local variables
!   REAL(R8Ki)                                     :: RotationMatrix(3,3)
   REAL(ReKi)                                     :: torque(3), DisplacedPosition(3)
   INTEGER(IntKi)                                 :: i         ! counter over the nodes


   ErrStat = ErrID_None
   ErrMsg  = ""
   

   !> Force: \f$\vec{F}^D = \sum\limits_{eS} \vec{F}^S \f$
   
   !> Moment: \f$\vec{M}^D = \sum\limits_{eS} \left\{ \vec{M}^S 
   !!         + \left\{ \left\{ \vec{p}^{OSR} + \vec{u}^S \right\} - \left\{ \vec{p}^{ODR} + \vec{u}^D \right\} \right\} \times \vec{F}^S \right\}\f$

   
!bjj note that we already checked that the following two conditions apply in this case:
!   if Src%FieldMask(MASKID_FORCE),  Dest%FieldMask(MASKID_FORCE) and Dest%FieldMask(MASKID_MOMENT)
!   if Src%FieldMask(MASKID_MOMENT), Dest%FieldMask(MASKID_MOMENT)

   if (Dest%FieldMask(MASKID_MOMENT) ) Dest%Moment = 0. ! whole array initialization; required to handle superposition of moments
   
   if (Src%FieldMask(MASKID_FORCE) ) THEN
      Dest%Force  = 0.     ! whole array initialization; required to handle superposition of forces
      do i = 1, Src%NNodes
         !if ( MeshMap%MapLoads(i)%OtherMesh_Element < 1 )  CYCLE ! would only happen if we had non-point elements (or nodes not contained in an element)
         
         ! F_d += F_s
         Dest%Force(:,MeshMap%MapLoads(i)%OtherMesh_Element) = Dest%Force(:,MeshMap%MapLoads(i)%OtherMesh_Element) + (Src%Force(:,i) / LoadsScaleFactor)
      end do
      Dest%Force  =  Dest%Force  * LoadsScaleFactor
      
      
      ! M_d += torque
      if ( Dest%FieldMask(MASKID_MOMENT) ) then
         
            ! if the distance (which can never be less than zero) is greater than "zero" and there is a
            ! force in the source mesh, then we need to add a moment to the destination mesh to account
            ! for the mismatch between points

         do i = 1, Src%NNodes
               DisplacedPosition =       SrcDisp%TranslationDisp(:,i) + SrcDisp%Position(:,i) &
                                       - ( DestDisp%TranslationDisp(:,MeshMap%MapLoads(i)%OtherMesh_Element) &
                                         + DestDisp%Position(       :,MeshMap%MapLoads(i)%OtherMesh_Element) )                              
               ! calculation torque vector based on offset force: torque = couple_arm X Force   
               torque = Src%Force(:,i) / LoadsScaleFactor !not torque yet, but we're doing this cross product in two step to avoid tempoary memory storage
               torque = CROSS_PRODUCT( DisplacedPosition, torque )
               Dest%Moment(:,MeshMap%MapLoads(i)%OtherMesh_Element) = Dest%Moment(:,MeshMap%MapLoads(i)%OtherMesh_Element) + torque                                             
         enddo
      endif      
      
   end if
   
      
   if (Src%FieldMask(MASKID_MOMENT) ) then

      ! M_d += M_s      
      do i = 1, Src%NNodes
         !if ( MeshMap%MapLoads(i)%OtherMesh_Element < 1 )  CYCLE ! would only happen if we had non-point elements (or nodes not contained in an element)

         Dest%Moment(:,MeshMap%MapLoads(i)%OtherMesh_Element) = Dest%Moment(:,MeshMap%MapLoads(i)%OtherMesh_Element) + (Src%Moment(:,i) / LoadsScaleFactor)
      end do
                  
   endif
      
   if (Dest%FieldMask(MASKID_MOMENT) )   Dest%Moment =  Dest%Moment * LoadsScaleFactor 

END SUBROUTINE Transfer_Loads_Point_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
!> Given a nearest-neighbor mapping, this routine generates the linearization matrices for loads transfer between point nodes on the meshes.
SUBROUTINE Linearize_Loads_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                !< The source mesh with loads fields allocated
   TYPE(MeshType),                 INTENT(IN   )  :: Dest               !< The destination mesh
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp            !< A mesh that contains the displacements associated with the source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp           !< A mesh that contains the displacements associated with the destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap            !< The mapping data structure (from Dest to Src)

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat            !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg             !< Error message if ErrStat /= ErrID_None

      ! local variables
   integer(intKi)                                 :: i,j,n, d_start, d_end, s_start, s_end
   real(r8Ki)                                     :: DisplacedPosition(3), SSmat(3,3)
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   character(*), parameter                        :: RoutineName = 'Linearize_Loads_Point_to_Point'
   

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !> Matrix \f$ M_{li}^D = M_{li} \f$, stored in modmesh_mapping::meshmaplinearizationtype::li,
   !! is allocated to be size Dest\%NNodes*3, Src\%NNodes*3.
   !! This is the matrix that maps each field of the source mesh to an augmented mesh.

   ! identity for forces and/or moments:
   if (.not. allocated(MeshMap%dM%li) ) then
      call AllocAry(MeshMap%dM%li, Dest%Nnodes*3, Src%Nnodes*3, 'dM%li', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
   elseif (size(MeshMap%dM%li,1) /= Dest%Nnodes*3 .or. size(MeshMap%dM%li,2) /= Src%Nnodes*3) then
      deallocate(MeshMap%dM%li)
      call AllocAry(MeshMap%dM%li, Dest%Nnodes*3, Src%Nnodes*3, 'dM%li', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
   end if
      
   MeshMap%dM%li = 0.0_R8Ki
   do i = 1, Src%NNodes
      !if ( MeshMap%MapLoads(i)%OtherMesh_Element < 1 )  CYCLE ! would only happen if we had non-point elements (or nodes not contained in an element)
         
      n = MeshMap%MapLoads(i)%OtherMesh_Element
      do j=1,3
         MeshMap%dM%li( (n-1)*3+j, (i-1)*3+j ) = 1.0_R8Ki
      end do
   end do
                              
         
      ! M_uD, M_uS, and m_f for moments:
         
   if (Dest%FieldMask(MASKID_MOMENT) .AND. Src%FieldMask(MASKID_FORCE) ) then
            
      !> Matrix \f$ M_{uD}^D \f$, stored in modmesh_mapping::meshmaplinearizationtype::m_ud,
      !! is allocated to be size Dest\%NNodes*3, Dest\%NNodes*3.
      !! This is the block matrix that maps the destination translation displacement field 
      !! of the source mesh to the moment field of the destination mesh.             
      
      if (.not. allocated(MeshMap%dM%M_uD) ) then
         call AllocAry(MeshMap%dM%M_uD, Dest%Nnodes*3, Dest%Nnodes*3, 'dM%M_uD', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      elseif (size(MeshMap%dM%M_uD,1) /= Dest%Nnodes*3 .or. size(MeshMap%dM%M_uD,2) /= Dest%Nnodes*3) then
         deallocate(MeshMap%dM%M_uD)
         call AllocAry(MeshMap%dM%M_uD, Dest%Nnodes*3, Dest%Nnodes*3, 'dM%M_uD', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      end if
            
      !> Matrix \f$ M_{uS}^D \f$, stored in modmesh_mapping::meshmaplinearizationtype::m_us,
      !! is allocated to be size Dest\%NNodes*3, Src\%NNodes*3.
      !! This is the block matrix that maps the source translation displacement field 
      !! of the source mesh to the moment field of the destination mesh.             
      
      if (.not. allocated(MeshMap%dM%M_uS) ) then
         call AllocAry(MeshMap%dM%M_uS, Dest%Nnodes*3, Src%Nnodes*3, 'dM%M_uS', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      elseif (size(MeshMap%dM%M_uS,1) /= Dest%Nnodes*3 .or. size(MeshMap%dM%M_uS,2) /= Src%Nnodes*3) then
         deallocate(MeshMap%dM%M_uS)
         call AllocAry(MeshMap%dM%M_uS, Dest%Nnodes*3, Src%Nnodes*3, 'dM%M_uS', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      end if
      
      
      MeshMap%dM%M_uD = 0.0_R8Ki
      MeshMap%dM%M_uS = 0.0_R8Ki
      
      do i = 1, Src%NNodes         
         n = MeshMap%MapLoads(i)%OtherMesh_Element

         s_start = (i-1)*3+1
         s_end   = s_start+2
         
         d_start = (n-1)*3+1
         d_end   = d_start+2
                             
         SSMat = SkewSymMat( Src%Force(:,i) )
         
         MeshMap%dM%M_uD( d_start:d_end, d_start:d_end ) = MeshMap%dM%M_uD( d_start:d_end, d_start:d_end ) + SSMat
         MeshMap%dM%M_uS( d_start:d_end, s_start:s_end ) = MeshMap%dM%M_uS( d_start:d_end, s_start:s_end ) - SSMat
               
      end do      
      
      !> > Matrix \f$ M_{fm}^D \f$, stored in modmesh_mapping::meshmaplinearizationtype::m_f,
      !! > is allocated to be size Dest\%NNodes*3, Src\%NNodes*3.
      !! > This is the block matrix that maps the force field of the source mesh to the moment field of the destination mesh.             
      
      if (.not. allocated(MeshMap%dM%m_f) ) then
         call AllocAry(MeshMap%dM%m_f, Dest%Nnodes*3, Src%Nnodes*3, 'dM%m_f', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      elseif (size(MeshMap%dM%m_f,1) /= Dest%Nnodes*3 .or. size(MeshMap%dM%m_f,2) /= Src%Nnodes*3) then
         deallocate(MeshMap%dM%m_f)
         call AllocAry(MeshMap%dM%m_f, Dest%Nnodes*3, Src%Nnodes*3, 'dM%m_f', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      end if
            
      
      MeshMap%dM%m_f = 0.0_R8Ki      
      
      do i = 1, Src%NNodes         
         n = MeshMap%MapLoads(i)%OtherMesh_Element
               
         s_start = (i-1)*3+1
         s_end   = s_start+2

         d_start = (n-1)*3+1
         d_end   = d_start+2
                                                              
         DisplacedPosition =       SrcDisp%TranslationDisp(:,i) +  SrcDisp%Position(:,i) &
                              - ( DestDisp%TranslationDisp(:,n) + DestDisp%Position(:,n) )                              
                                             
         MeshMap%dM%m_f( d_start:d_end, s_start:s_end ) = SkewSymMat( DisplacedPosition )
               
      end do
                           
   endif    
      

END SUBROUTINE Linearize_Loads_Point_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
! This routine computes a scaling factor for loads fields to reduce numerical problems in computions that involve accelerations 
! and forces/moments.
FUNCTION GetLoadsScaleFactor( Src )

   TYPE(MeshType),                 INTENT(IN   )  :: Src        !< The source mesh with loads fields allocated
   REAL(ReKi)                                     :: GetLoadsScaleFactor  !< scaling factor for loads fields
   
   
   ! LOCAL:
   INTEGER                                         :: I, j
   REAL(ReKi)                                      :: MaxLoad

   
   GetLoadsScaleFactor = 1.0
   MaxLoad             = 0.0
   
   IF ( Src%FIELDMASK( MASKID_FORCE ) ) THEN
      
      DO I=1,Src%Nnodes
         DO J=1,3
            MaxLoad = MAX(MaxLoad, ABS(Src%Force(j,I) ) )
         END DO                  
      END DO
      
   END IF
   

   IF ( Src%FIELDMASK( MASKID_MOMENT ) ) THEN
      
      DO I=1,Src%Nnodes
         DO J=1,3
            MaxLoad = MAX(MaxLoad, ABS(Src%Moment(j,I) ) )
         END DO                  
      END DO
      
   END IF
   
   IF ( MaxLoad > 10. ) THEN
      GetLoadsScaleFactor = 10**MIN( NINT(log10(MaxLoad)), 15 )  ! Let's not get carried away and cause overflow; 10E15 is as far as we'll go
   END IF
   

END FUNCTION GetLoadsScaleFactor
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CreateLoadMap_P_to_P( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            ! The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

   
      ! in following call, Src is mesh to loop over, finding a corresponding point for each point in Dest
   CALL CreateMapping_NearestNeighbor( Src, Dest, MeshMap%MapLoads, ELEMENT_POINT, ELEMENT_POINT, ErrStat, ErrMsg )
         
END SUBROUTINE CreateLoadMap_P_to_P
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CreateMotionMap_P_to_P( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            ! The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                 :: i
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! in following call, Dest is mesh to looped over, finding a corresponding point for each point in Src
      CALL CreateMapping_NearestNeighbor( Dest, Src, MeshMap%MapMotions, ELEMENT_POINT, ELEMENT_POINT, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateMotionMap_P_to_P')
         IF (ErrStat >= AbortErrLev) RETURN
            
      ! bjj: for consistant definition of couple_arm (i.e. p_ODR-p_OSR), let's multiply by -1
      do i=1,SIZE(MeshMap%MapMotions)
         MeshMap%MapMotions(i)%couple_arm = -1.0_ReKi*MeshMap%MapMotions(i)%couple_arm
      end do
                     
END SUBROUTINE CreateMotionMap_P_to_P        
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that transfers data from a line2 mesh to another line2 mesh.
SUBROUTINE Transfer_Line2_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),         INTENT(IN   ) ::  Src      !< source Line2 mesh
   TYPE(MeshType),         INTENT(INOUT) ::  Dest     !< destination Line2 mesh
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  SrcDisp  !< a "functional" sibling of the source mesh; Src contains loads and SrcDisp contains TranslationDisp and Orientation
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  DestDisp !< a "functional" sibling of the destination mesh; Dest contains loads and DestDisp contains TranslationDisp and Orientation

   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap   !< mapping between Src and Dest meshes


   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat   !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg    !< Error message if ErrStat /= ErrID_None

      ! local variables
   REAL(ReKi)                            :: LoadsScaleFactor  ! Scaling factor for loads (to help with numerical issues)
   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*), PARAMETER               :: RoutineName = 'Transfer_Line2_to_Line2'   
   
   ErrStat = ErrID_None
   ErrMsg  = ''

   !.................
   ! Check to ensure that both the source and destination meshes are composed of Line2 elements
   !.................   
   if (Src%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Source mesh must have one or more Line2 elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif

   if (Dest%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Destination mesh must have one or more Line2 elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif
   

   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------

   if ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) then

      !........................
      ! Start: Create Mapping data (if remap is true)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMotionMap_L2_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN

      endif !remapping

      !........................
      ! Start: Transfer data
      !........................
         
      CALL Transfer_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN

   endif !algorithm for motions/scalars


   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if ( HasLoadFields(Src) ) then

      IF (.not. PRESENT(SrcDisp) .OR. .NOT. PRESENT(DestDisp) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'SrcDisp and DestDisp arguments are required for load transfer.', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF
            
      !.................
      ! other checks for available mesh fields (now done in AllocMapping routine)
      !.................


      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then
                  
         CALL CreateLoadMap_L2_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN         
            
      end if

      !........................
      ! Transfer data
      !........................

      LoadsScaleFactor = GetLoadsScaleFactor ( Src )
      
      ! first, we take the source fields and transfer them to fields on the augmented source mesh:
      !  (we're also taking the SrcDisp field and putting it on our augmented mesh)
      CALL Transfer_Src_To_Augmented_Ln2_Src( Src, MeshMap, ErrStat2, ErrMsg2, SrcDisp, LoadsScaleFactor ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
      
      ! then we lump the loads from the augmented source mesh:
      CALL Lump_Line2_to_Point( MeshMap%Augmented_Ln2_Src,  MeshMap%Lumped_Points_Src,  ErrStat2, ErrMsg2, LoadsScaleFactor=LoadsScaleFactor  ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
      
      CALL Transfer_Loads_Point_to_Line2( MeshMap%Lumped_Points_Src, Dest, MeshMap, ErrStat2, ErrMsg2, MeshMap%Augmented_Ln2_Src, DestDisp, LoadsScaleFactor )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN


   end if ! algorithm for loads


END SUBROUTINE Transfer_Line2_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine that computes the linearization matrices for data transfers from a mesh to a line2 mesh.
!! Mapping equations are as follows: 
!!
!! Rotational Displacement: \f$ \frac{\partial M_\Lambda}{\partial x} = \begin{bmatrix} M_{mi} \end{bmatrix} \f$ for source fields \f$ \left\{ \vec{\theta}^S\right\}\f$
!!
!! Translational Displacement: \f$ \frac{\partial M_u}{\partial x} = \begin{bmatrix} M_{mi} & M_{f_{\times p}} \end{bmatrix} \f$ for source fields
!! \f$ \left\{ \begin{matrix}
!!      \vec{u}^S \\
!!      \vec{\theta}^S
!! \end{matrix} \right\} \f$
!!
!! Rotational Velocity: \f$ \frac{\partial M_\omega}{\partial x} = \begin{bmatrix} M_{mi} \end{bmatrix} \f$ for source fields \f$\left\{\vec{\omega}^S\right\}\f$
!!
!! Translational Velocity: \f$ \frac{\partial M_v}{\partial x} = \begin{bmatrix} M_{tv\_uD} & M_{tv\_uS} & M_{mi} & M_{f_{\times p}} \end{bmatrix} \f$ for source fields
!! \f$ \left\{ \begin{matrix}
!!      \vec{u}^D \\
!!      \vec{u}^S \\
!!      \vec{v}^S \\
!!      \vec{\omega}^S
!! \end{matrix} \right\} \f$
!!
!! Rotational Acceleration: \f$ \frac{\partial M_\alpha}{\partial x} = \begin{bmatrix} M_{mi} \end{bmatrix} \f$ for source fields \f$\left\{\vec{\alpha}^S\right\}\f$
!!
!! Translational Acceleration: \f$ \frac{\partial M_a}{\partial x} = \begin{bmatrix} M_{ta\_uD} & M_{ta\_uS} & M_{ta\_rv} & M_{mi} & M_{f_{\times p}} \end{bmatrix} \f$ for source fields
!! \f$ \left\{ \begin{matrix}
!!      \vec{u}^D \\
!!      \vec{u}^S \\
!!      \vec{\omega}^S \\
!!      \vec{a}^S \\
!!      \vec{\alpha}^S
!! \end{matrix} \right\} \f$
!!
!! Scalar quantities: \f$ \frac{\partial M_S}{\partial x} = \begin{bmatrix} M_{mi} \end{bmatrix} \f$ for source fields \f$\left\{\vec{S}^S\right\}\f$
!!
!! Forces: \f$ \frac{\partial M_f}{\partial x} = \begin{bmatrix} M_{li} \end{bmatrix} \f$ for source fields \f$\left\{\vec{f}^S\right\}\f$
!!
!! Moments: \f$ \frac{\partial M_m}{\partial x} = \begin{bmatrix} M_{uDm} & M_{uSm} & M_{fm} & M_{li} \end{bmatrix} \f$ for source fields
!! \f$ \left\{ \begin{matrix}
!!      \vec{u}^D \\
!!      \vec{u}^S \\
!!      \vec{f}^S \\
!!      \vec{m}^S
!! \end{matrix} \right\} \f$
!!
!! \f$M_{mi}\f$ is modmesh_mapping::meshmaplinearizationtype::mi \n
!! \f$M_{f_{\times p}}\f$ is modmesh_mapping::meshmaplinearizationtype::fx_p \n
!! \f$M_{tv\_uD}\f$ is modmesh_mapping::meshmaplinearizationtype::tv_uD \n
!! \f$M_{tv\_uS}\f$ is modmesh_mapping::meshmaplinearizationtype::tv_uS \n
!! \f$M_{ta\_uD}\f$ is modmesh_mapping::meshmaplinearizationtype::ta_uD \n
!! \f$M_{ta\_uS}\f$ is modmesh_mapping::meshmaplinearizationtype::ta_uS \n
!! \f$M_{ta\_rv}\f$ is modmesh_mapping::meshmaplinearizationtype::ta_rv \n
!! \f$M_{li}\f$ is modmesh_mapping::meshmaplinearizationtype::li \n
!! \f$M_{uSm}\f$ is modmesh_mapping::meshmaplinearizationtype::m_us \n
!! \f$M_{uDm}\f$ is modmesh_mapping::meshmaplinearizationtype::m_ud \n
!! \f$M_{fm}\f$ is modmesh_mapping::meshmaplinearizationtype::m_f 
SUBROUTINE Linearize_Line2_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),         INTENT(IN   ) :: Src       !< source Line2 mesh
   TYPE(MeshType),         INTENT(IN   ) :: Dest      !< destination mesh
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) :: SrcDisp   !< a "functional" sibling of the source mesh; Src contains loads and SrcDisp contains TranslationDisp and Orientation
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) :: DestDisp  !< a "functional" sibling of the destination mesh; Dest contains loads and DestDisp contains TranslationDisp and Orientation

   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap   !< mapping between Src and Dest meshes

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat   !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg    !< Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*), PARAMETER               :: RoutineName = 'Linearize_Line2_to_Line2'   

   real(R8Ki), allocatable               :: M_A(:,:)        ! linearization matrix for augmented source mesh
   real(R8Ki), allocatable               :: M_SL_fm(:,:)    ! linearization matrix for source-mesh lumped force component of moment
   real(R8Ki), allocatable               :: M_SL_uSm(:,:)   ! linearization matrix for source-mesh lumped translational displacement component of moment
   real(R8Ki), allocatable               :: M_SL_li(:,:)    ! linearization matrix for source-mesh lumped load "identity" component 
        
   ErrStat = ErrID_None
   ErrMsg  = ''


   ! ------------------------------------------------------------------------------------------------------------------------------
   !> ### Mapping and Linearization of Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------

   if ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) then

      !........................
      !> * Create mapping (if remap is true)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMotionMap_L2_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN

      endif !remapping

      !........................
      !> * Get linearization matrices of data transfer
      !........................
         
      CALL Linearize_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN

   endif !algorithm for motions/scalars


   ! ------------------------------------------------------------------------------------------------------------------------------
   !> ### Mapping and Linearization of Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if ( HasLoadFields(Src) ) then

      IF (.not. PRESENT(SrcDisp) .OR. .NOT. PRESENT(DestDisp) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'SrcDisp and DestDisp arguments are required for load transfer linearization.', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF
            
      !.................
      ! other checks for available mesh fields (now done in AllocMapping routine)
      !.................


      !........................
      !> * Create mapping (if remap is true)
      ! Note that we also create a motion mapping for the source and destination displacement meshes (SrcDisp and DestDisp)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then
                  
         CALL CreateLoadMap_L2_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN         
            
      end if
      
      !........................
      !> * Linearize data
      !........................            
      
      !>  + Get individual transformation matrices:
      
      !>   1. Get the matrix that transfers the source fields to the augmented source mesh.
      !! (We're also taking the force field and and (source) translational displacement field and 
      !! putting them on our [intermediate] augmented mesh.)
      CALL Linearize_Src_To_Augmented_Ln2_Src( Src, MeshMap, ErrStat2, ErrMsg2, SrcDisp ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
         
      call move_alloc( MeshMap%dM%li, M_A )
         
! ^^^ size of M_A is 3*MeshMap%Augmented_Ln2_Src%NNodes X 3*Src%Nnodes
                  
      !>   2. Get the matrices that lump the loads on the augmented source mesh.
      !! (We're also taking the force field and putting it on our lumped mesh.)
      CALL Linearize_Lump_Line2_to_Point( MeshMap%Augmented_Ln2_Src,  MeshMap%Lumped_Points_Src, MeshMap%dM, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
      call move_alloc( MeshMap%dM%m_uD, M_SL_uSm )    
      call move_alloc( MeshMap%dM%m_f,  M_SL_fm )    
      call move_alloc( MeshMap%dM%li,   M_SL_li )
         
! ^^^ size of M_SL_li (as well as M_SL_uSm and M_SL_fm) is 3*MeshMap%Augmented_Ln2_Src%NNodes X 3*MeshMap%Augmented_Ln2_Src%NNodes
                                    
         !>   3. (and 4) Get the matrices transfering the lumped (point) loads to line2 loads.
      CALL Linearize_Loads_Point_to_Line2( MeshMap%Lumped_Points_Src, Dest, MeshMap, ErrStat2, ErrMsg2, MeshMap%Augmented_Ln2_Src, DestDisp )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

! ^^^ size of dM%li is 3*Dest%NNodes X 3*MeshMap%Augmented_Ln2_Src%NNodes (m_uD is square, of size 3*Dest%NNodes X 3*Dest%NNodes)
! need to return size of dM%li as 3*Dest%NNodes X 3*Src%NNodes:
                                                      
         
      !>  + Multiply individual transformation matrices to get full linearization matrices         
      CALL FormMatrix_FullLinearization( MeshMap%dM, M_A, M_SL_fm, M_SL_uSm, M_SL_li, ErrStat2, ErrMsg2 )    
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                                          
      call cleanup()               

   end if ! algorithm for loads

   
contains
subroutine cleanup()

   if (allocated(M_A     )) deallocate(M_A    )
   if (allocated(M_SL_li )) deallocate(M_SL_li)
   if (allocated(M_SL_uSm)) deallocate(M_SL_uSm)
   if (allocated(M_SL_fm )) deallocate(M_SL_fm)

end subroutine cleanup

END SUBROUTINE Linearize_Line2_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
!> Given a mapping, this routine transfers the loads from nodes on a point-element mesh to nodes on another Line2 mesh.
SUBROUTINE Transfer_Loads_Point_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp, LoadsScaleFactor )

   TYPE(MeshType),                 INTENT(IN   )  :: Src       !< The source (Line2) mesh with loads fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest      !< The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   !< The mapping data

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    !< Error message if ErrStat /= ErrID_None
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp   !< The source mesh's cooresponding position mesh
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp  !< The destination mesh's cooresponding position mesh
   REAL(ReKi),                     INTENT(IN)     :: LoadsScaleFactor  !< Scaling factor for loads (to help with numerical issues)

      ! local variables
   REAL(ReKi)                                     :: torque(3), DisplacedPosition(3)

   INTEGER(IntKi)                                 :: i,j       ! loop counters
   INTEGER(IntKi)                                 :: jElem     ! element number
   INTEGER(IntKi)                                 :: jNode     ! node number

   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   character(*), parameter                        :: RoutineName = 'Transfer_Loads_Point_to_Line2'
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""

      
   !> Start with transferring the loads like point-to-point (except split between two nodes of the dest element).
   ! [because the OtherMesh_Element is line2, we can't just call the routine here]
         
   if ( Dest%FieldMask(MASKID_MOMENT) ) Dest%Moment = 0._ReKi     ! whole array initialization; required to handle superposition of loads
   
      ! ---------------------------- Force ------------------------------------------------
      ! Transfer this source force (lumped) to destination force and/or moment fields (currently lumped):

      ! split the point source loads across the line2 destination element (note that they are still point loads on line2 elements here)
   
   if ( Src%FieldMask(MASKID_Force) ) THEN
      Dest%Force  = 0._ReKi     ! whole array initialization; required to handle superposition of loads
      ! ( Dest%FieldMask(MASKID_Force) ) is .TRUE.

      DO i = 1, Src%Nnodes  ! the source nodes (containing lumped loads)
            
         jElem = MeshMap%MapLoads(i)%OtherMesh_Element
         DO j = 1,2 ! number of nodes on dest 
            jNode = Dest%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(j)
            Dest%Force(:,jNode) = Dest%Force(:,jNode) + (Src%Force(:,i)/LoadsScaleFactor)*MeshMap%MapLoads(i)%shape_fn(j)             
         END DO !j

      END DO !loop through source nodes
      
            
            ! calculate the moment due to force
      !if ( Dest%FieldMask(MASKID_MOMENT) ) then (now must always be true)
                  
         ! if the distance (which can never be less than zero) is greater than "zero" and there is a
         ! force in the source mesh, then we need to add a moment to the destination mesh to account
         ! for the mismatch between points

        ! M_d += torque

         DO i = 1, Src%Nnodes  ! the source nodes
            
            jElem = MeshMap%MapLoads(i)%OtherMesh_Element
            DO j = 1,2 ! number of nodes on dest 
               jNode = Dest%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(j)
                        
               DisplacedPosition =   Src%Position(:,i)     +  SrcDisp%TranslationDisp(:,i)     &
                                  - Dest%Position(:,jNode) - DestDisp%TranslationDisp(:,jNode)  
                              
               torque = Src%Force(:,i) / LoadsScaleFactor !not torque yet, but we're doing this cross product in two steps to avoid tempoary memory storage
               torque = CROSS_PRODUCT( DisplacedPosition, torque )               
               Dest%Moment(:,jNode) = Dest%Moment(:,jNode) + torque*MeshMap%MapLoads(i)%shape_fn(j)            
            END DO !j
            

         END DO !loop through source elements


      !end if ! Dest has moment, too      
                  
      
   end if !source has force
   
      
      ! ---------------------------- Moments ------------------------------------------------
      ! Transfer this source moment to the destination moment field

   if ( Src%FieldMask(MASKID_Moment) ) then !
      
      
      DO i = 1,  Src%Nnodes  ! the source nodes
         
         jElem = MeshMap%MapLoads(i)%OtherMesh_Element
         DO j = 1,2 ! number of nodes on dest 
            jNode = Dest%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(j)
            Dest%Moment(:,jNode) = Dest%Moment(:,jNode) + (Src%Moment(:,i)/LoadsScaleFactor)*MeshMap%MapLoads(i)%shape_fn(j)             
         END DO !j

      END DO !loop through source nodes
      
      
   end if !src has moment field
      
#ifdef MESH_DEBUG     
   CALL Create_PointMesh(Dest, MeshMap%Lumped_Points_Dest, ErrStat, ErrMsg)
   IF (Dest%FieldMask(maskid_force))  MeshMap%Lumped_Points_Dest%Force  = Dest%Force
   IF (Dest%FieldMask(maskid_moment)) MeshMap%Lumped_Points_Dest%Moment = Dest%Moment
#endif

   !> now we can convert the lumped loads on the dest (line2) mesh into distributed loads on the same mesh:   
   CALL Convert_Point_To_Line2_Loads(Dest, MeshMap, ErrStat2, ErrMsg2, DestDisp)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      !IF (ErrStat >= AbortErrLev) RETURN
   
   if (Dest%FieldMask(MASKID_FORCE) ) Dest%Force  = Dest%Force  * LoadsScaleFactor     ! whole array initialization
   if (Dest%FieldMask(MASKID_MOMENT)) Dest%Moment = Dest%Moment * LoadsScaleFactor     ! whole array initialization
   
      
END SUBROUTINE Transfer_Loads_Point_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
!> Given a mapping, this routine contains the linearization matrices of the load transfers from nodes on a point-element mesh to 
!! nodes on another Line2 mesh. At the end of this routine, the MeshMap data structure will contain the four matrices for loads
!! transfer from a point to line2 mesh: \f$ M_{um}\f$, \f$ M_{tm}\f$, \f$M_{fm}\f$, and \f$M_{li}\f$, where
!! \f{equation}{ 
!!   \begin{aligned}
!! \left\{ \begin{matrix} \Delta\vec{u}^D \\  \Delta\vec{\theta}^D \\  \Delta\vec{f}^D \\   \Delta\vec{m}^D \end{matrix} \right\} & = 
!! \begin{bmatrix} M_{mi} & M_{f_{\times p}} & 0      & 0      \\
!!                 0        & M_{mi}         & 0      & 0      \\
!!                 0        & 0              & M_{li} & 0      \\
!!                 M_{um}   & M_{tm}         & M_{fm} & M_{li} \\ \end{bmatrix}
!! \left\{ \begin{matrix} \Delta\vec{u}^S \\  \Delta\vec{\theta}^S \\  \Delta\vec{F}^S \\   \Delta\vec{M}^S \end{matrix} \right\} \\ & =   
!! \begin{bmatrix} I                                                         & 0 & 0                                                                                                        & 0        \\
!!                 0                                                         & I & 0                                                                                                        & 0        \\
!!                 0                                                         & 0 &  \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1}                                                          & 0        \\
!!                -\begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1}M_{um}^{DL} & 0 & -\begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1}M_{fm}^{DL}\begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} & \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} \\ \end{bmatrix}
!! \begin{bmatrix} M_{mi}   & M_{f_{\times p}} & 0        & 0        \\
!!                 0        & M_{mi}           & 0        & 0        \\
!!                 0        & 0                & M_{li}^D & 0        \\
!!                 0        & M_{tm}^D         & M_{fm}^D & M_{li}^D \\ \end{bmatrix}
!! \left\{ \begin{matrix} \Delta\vec{u}^S \\  \Delta\vec{\theta}^S \\  \Delta\vec{F}^S \\   \Delta\vec{M}^S \end{matrix} \right\}
!! \end{aligned}
!! \f}
!! "S" refers to the source (point) mesh;"D" refers to the destination (line2) mesh; "DL" refers to the destination mesh with lumped (point) loads. \n
!! \f$M_{li}\f$ is modmesh_mapping::meshmaplinearizationtype::li \n
!! \f$M_{M_u}\f$ is modmesh_mapping::meshmaplinearizationtype::m_us \n
!! \f$M_{M_t}\f$ is modmesh_mapping::meshmaplinearizationtype::m_ud \n
!! \f$M_{M_f}\f$ is modmesh_mapping::meshmaplinearizationtype::m_f
SUBROUTINE Linearize_Loads_Point_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),                 INTENT(IN   )  :: Src       !< The source (Line2) mesh with loads fields allocated
   TYPE(MeshType),                 INTENT(IN   )  :: Dest      !< The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   !< The mapping data from Src to Dest

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    !< Error message if ErrStat /= ErrID_None
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp   !< The source mesh's cooresponding position mesh
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp  !< The destination mesh's cooresponding position mesh

      ! local variables
   REAL(R8Ki)                                     :: DisplacedPosition(3), SSMat(3,3)

   REAL(R8Ki), ALLOCATABLE                        :: li_D(:,:)
   REAL(R8Ki), ALLOCATABLE                        :: muS_D(:,:)
   REAL(R8Ki), ALLOCATABLE                        :: muD_D(:,:)
   REAL(R8Ki), ALLOCATABLE                        :: mf_D(:,:)
   
   integer(intKi)                                 :: i, j, jElem, k, n, d_start, d_end, s_start, s_end

   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   character(*), parameter                        :: RoutineName = 'Linearize_Loads_Point_to_Line2'
   
      
   ErrStat = ErrID_None
   ErrMsg  = ""

   
   !> Start by determining how the point load is split based on its projected location in the mapped destination Line2 element.
   !! (This is like Linearize_Loads_Point_to_Point, except we need to loop on elements instead of nodes [because of nearest neighbor mapping routine])
   
   !> > Matrix \f$ M_{li}^D \f$, stored in li_D,
   !! > is allocated to be size Dest\%NNodes*3, Src\%NNodes*3.
   !! > This is the matrix that maps each field of the source point mesh to another point mesh.
      
      ! "identity" (portion of node) for forces and/or moments:
   if (.not. allocated(li_D) ) then
      call AllocAry(li_D, Dest%Nnodes*3, Src%Nnodes*3, 'li_D', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) then
            call cleanup()
            RETURN
         end if         
   end if
      
   li_D = 0.0_R8Ki
   do i = 1, Src%NNodes
         
      jElem = MeshMap%MapLoads(i)%OtherMesh_Element
      do j=1,NumNodes(ELEMENT_LINE2)
         n = Dest%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(j)
      
         do k=1,3
            li_D( (n-1)*3+k, (i-1)*3+k ) = MeshMap%MapLoads(i)%shape_fn(j)
         end do
      end do
      
   end do
       
      ! m_uD, m_uS, and m_f for moments:
         
   if (Dest%FieldMask(MASKID_MOMENT) .AND. Src%FieldMask(MASKID_FORCE) ) then
            
      !> > Matrix \f$ M_{uD}^D \f$, stored in muD_D,
      !! > is allocated to be size Dest\%NNodes*3, Dest\%NNodes*3.
      !! > This is the block matrix that maps the u (translational displacement) field 
      !! > of the destination mesh to the moment field of the destination mesh.             
      
      if (.not. allocated(muD_D) ) then
         call AllocAry(muD_D, Dest%Nnodes*3, Dest%Nnodes*3, 'muD_D', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) then
               call cleanup()
               RETURN
            end if
      end if   
      
      !> > Matrix \f$ M_{uS}^D \f$, stored in muS_D,
      !! > is allocated to be size Dest\%NNodes*3, Src\%NNodes*3.
      !! > This is the block matrix that maps the u (translational displacement) field 
      !! > of the source mesh to the moment field of the destination mesh.             
      
      if (.not. allocated(muS_D) ) then
         call AllocAry(muS_D, Dest%Nnodes*3, Src%Nnodes*3, 'muS_D', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) then
               call cleanup()
               RETURN
            end if
      end if
                 
      muD_D = 0.0_R8Ki      
      muS_D = 0.0_R8Ki         
      do i = 1, Src%NNodes         
         jElem = MeshMap%MapLoads(i)%OtherMesh_Element
               
         s_start = (i-1)*3+1
         s_end   = s_start+2         
         
         do j=1,NumNodes(ELEMENT_LINE2)
            n = Dest%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(j)
            
            d_start = (n-1)*3+1
            d_end   = d_start+2
         
            SSmat = SkewSymMat( Src%Force(:,i)*MeshMap%MapLoads(i)%shape_fn(j) )
            
            muD_D( d_start:d_end, d_start:d_end ) = muD_D( d_start:d_end, d_start:d_end ) + SSmat
            muS_D( d_start:d_end, s_start:s_end ) = muS_D( d_start:d_end, s_start:s_end ) - SSmat           
               
         end do !j
         
      end do      
      
      
      !> > Matrix \f$ M_{fm}^D \f$, stored in mf_D,
      !! > is allocated to be size Dest\%NNodes*3, Src\%NNodes*3.
      !! > This is the block matrix that maps the force field of the source mesh to the moment field of the destination mesh.             
      
      if (.not. allocated(mf_D) ) then
         call AllocAry(mf_D, Dest%Nnodes*3, Src%Nnodes*3, 'mf_D', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) then
               call cleanup()
               RETURN
            end if
      end if
            
      
      mf_D = 0.0_R8Ki         
      do i = 1, Src%NNodes         
         jElem = MeshMap%MapLoads(i)%OtherMesh_Element
               
         s_start = (i - 1)*3+1
         s_end   = s_start+2         
         do j=1,NumNodes(ELEMENT_LINE2)
            n = Dest%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(j)
            
            d_start = (n-1)*3+1
            d_end   = d_start+2
            
            DisplacedPosition =       SrcDisp%TranslationDisp(:,i) +  SrcDisp%Position(:,i) &
                                 - ( DestDisp%TranslationDisp(:,n) + DestDisp%Position(:,n) )

            DisplacedPosition = DisplacedPosition*MeshMap%MapLoads(i)%shape_fn(j)  ! multiply the fraction here so we don't need to apply to any of the terms below

            mf_D( d_start:d_end, s_start:s_end ) = SkewSymMat( DisplacedPosition )
                           
         end do !j
         
      end do
            
   endif    
                  
   !> Then we can obtain the matrices that convert the lumped loads on the dest (line2) mesh into distributed loads on the same mesh:  
   !! For linearization analysis, we assume this destination mesh already contains the distributed loads (i.e., values at the operating point).
   CALL Linearize_Convert_Point_To_Line2_Loads(Dest, MeshMap, ErrStat2, ErrMsg2, DestDisp)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         RETURN
      end if
   
   !> Finally, we multiply the matrices together and return the non-zero block matrices that remain.
   
      n=Dest%Nnodes*3
      
   !> > Matrix \f$ M_{li} = \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} M_{li}^{D} \f$, stored in modmesh_mapping::meshmaplinearizationtype::li,
   !! > is allocated to be size Dest\%NNodes*3, Src\%NNodes*3.
   !! > This is the matrix that maps each field of the source mesh to an augmented mesh.
   !! > \n Note that we solve the equation \f$ M_{li}^{DL} M_{li} = M_{li}^{D} \f$ for \f$ M_{li} \f$.
      
   ! solve this before M_{fm} so that we can use M_{li} in the equation for M_{fm}

      ! solve for M_{li}:
   CALL LAPACK_getrs(TRANS='N',N=n,A=MeshMap%LoadLn2_A_Mat,IPIV=MeshMap%LoadLn2_A_Mat_piv, B=li_D, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return      
                 
   call move_alloc(li_D, MeshMap%dM%li)
      
   
   if (Dest%FieldMask(MASKID_MOMENT) .AND. Src%FieldMask(MASKID_FORCE) ) then               
      
      !> > Matrix \f$ M_{uSm} = -\begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} M_{uSm}^{DL} \f$      
      !! > stored in modmesh_mapping::meshmaplinearizationtype::m_us, is allocated to be size Dest\%NNodes*3, Src\%NNodes*3.
      !! > This is the matrix that maps the source u (translational displacement) to the destination moment.
      !! > \n Note that we solve the equation \f$ M_{li}^{DL} M_{uSm} = M_{uSm}^{DL} \f$ for \f$ M_{uSm} \f$.

         ! solve for M_{uSm}      
      CALL LAPACK_getrs(TRANS='N',N=n,A=MeshMap%LoadLn2_A_Mat,IPIV=MeshMap%LoadLn2_A_Mat_piv, B=muS_D, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return      
         
      call move_alloc(muS_D, MeshMap%dM%M_uS)        
       
      !> > Matrix \f$ M_{uDm} = -\begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} M_{uDm}^{DL} + \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} M_{uDm}^D \f$      
      !! > stored in modmesh_mapping::meshmaplinearizationtype::m_ud, is allocated to be size Dest\%NNodes*3, Dest\%NNodes*3.
      !! > This is the matrix that maps the destination u (translational displacement) to the destination moment.
      !! > \n Note that we solve the equation \f$ M_{li}^{DL} M_{uDm} = M_{uDm}^D - M_{uDm}^{DL} \f$ for \f$ M_{uDm} \f$.

         ! get RHS of equation:
      muD_D = muD_D - MeshMap%dM%m_uD     

         ! solve for M_{uDm}      
      CALL LAPACK_getrs(TRANS='N',N=n,A=MeshMap%LoadLn2_A_Mat,IPIV=MeshMap%LoadLn2_A_Mat_piv, B=muD_D, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return      
         
      call move_alloc(muD_D, MeshMap%dM%M_uD)      
            
      !> > Matrix  \f$ M_{fm} \f$ is stored in modmesh_mapping::meshmaplinearizationtype::m_f and is allocated to be size Dest\%NNodes*3, Src\%NNodes*3.
      !! > This is the matrix that maps the force to the moment. \n
      !! > \f$ \begin{aligned} M_{fm} & = 
      !!   - \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} M_{fm}^{DL}\begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1}M_{li}^D 
      !!   + \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} M_{fm}^D \\ & =
      !!     \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} \left( - M_{fm}^{DL} M_{li} + M_{fm}^D \right) \end{aligned}\f$, 
      !! > \n Note that we solve the equation \f$ M_{li}^{DL} M_{fm} = M_{fm}^D - M_{fm}^{DL} M_{li}  \f$ for \f$ M_{fm} \f$.
            
         ! get RHS of equation
      mf_D = mf_D - matmul( MeshMap%dM%m_f, MeshMap%dM%li )
      
         ! solve for M_{fm}      
      CALL LAPACK_getrs(TRANS='N',N=n,A=MeshMap%LoadLn2_A_Mat,IPIV=MeshMap%LoadLn2_A_Mat_piv, B=mf_D, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return      
                  
      call move_alloc(mf_D, MeshMap%dM%m_f)
                              
   end if
   
   
   call cleanup()
      
contains

   subroutine cleanup()
      if (allocated(li_D)) deallocate(li_D)
      if (allocated(muD_D)) deallocate(muD_D)
      if (allocated(muS_D)) deallocate(muS_D)
      if (allocated(mf_D)) deallocate(mf_D)
   end subroutine cleanup
   
END SUBROUTINE Linearize_Loads_Point_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine takes the lumped loads on nodes of a (line2) mesh and converts them to loads distributed across the line2 elements.
SUBROUTINE Convert_Point_To_Line2_Loads(Dest, MeshMap, ErrStat, ErrMsg, DestDisp)

   TYPE(MeshType),                 INTENT(INOUT)  :: Dest                           !< The mesh (on input, the nodal loads values are lumped; on output they are distributed along the element)
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp                       !< The mesh that contains the translation displacement for these loads
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap                        !< The mapping data

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                         !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi) :: jElem, n, i,j, n1, n2
   REAL(ReKi)     :: a_vec(3), sum_f(3), crossProd(3)
   REAL(ReKi)     :: c
   
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   character(*), parameter                        :: RoutineName = 'Convert_Point_To_Line2_Loads'
   
   
   n=3*Dest%Nnodes !also SIZE(MeshMap%LoadLn2_F,1) and SIZE(MeshMap%LoadLn2_M,1)
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! Convert the forces on each element
   
   IF ( Dest%FieldMask(MASKID_Force) ) then  
      
      MeshMap%LoadLn2_F = RESHAPE( Dest%Force, (/ n, 1 /) )
      
      ! After following call, LoadLn2_F contains the distributed forces:
      
      CALL LAPACK_getrs(TRANS='N',N=n,A=MeshMap%LoadLn2_A_Mat,IPIV=MeshMap%LoadLn2_A_Mat_piv, &
                        B=MeshMap%LoadLn2_F, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
      
      ! Transfer forces to the mesh fields
   
      Dest%Force =  RESHAPE( MeshMap%LoadLn2_F, (/ 3, Dest%Nnodes /) )
      
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN      
                  
   END IF ! Force
   
   
      ! Convert the moments on each element
      
   IF ( Dest%FieldMask(MASKID_Moment) ) then  
      
      MeshMap%LoadLn2_M = RESHAPE( Dest%Moment, (/ n, 1 /) )       
      
      IF ( Dest%FieldMask(MASKID_Force) ) then
                                                               
         DO jElem = 1,Dest%ElemTable(ELEMENT_LINE2)%nelem

            n1=Dest%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(1)
            n2=Dest%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(2)
            
            c = Dest%ElemTable(ELEMENT_LINE2)%Elements(jElem)%det_jac/6.0_ReKi  ! this contains an extra factor of 1/2, which comes from omitting it from the a_vec term
                                                
            a_vec =   ( DestDisp%Position(:,n2) + DestDisp%TranslationDisp(:,n2) ) &
                    - ( DestDisp%Position(:,n1) + DestDisp%TranslationDisp(:,n1) )
            
            ! We solved LoadLn2_F for the distributed value, so we'll use it here: (the B matrix is just the cross product terms)
            sum_f = Dest%Force(:,n1) + Dest%Force(:,n2)             
            crossProd = c * cross_product( a_vec, sum_f)
                                 
            ! subtract the force (using the distributed values) from the lumped moment terms:
            i = (n1-1)*3+1;  j = i+2
            MeshMap%LoadLn2_M(i:j,1)=MeshMap%LoadLn2_M(i:j,1) - crossProd
            
            i = (n2-1)*3+1;  j = i+2            
            MeshMap%LoadLn2_M(i:j,1)=MeshMap%LoadLn2_M(i:j,1) + crossProd
         END DO        
                           
      END IF ! moment due to force
      
      ! After following call, LoadLn2_M contains the distributed moments:
      
      CALL LAPACK_getrs(TRANS='N',N=n,A=MeshMap%LoadLn2_A_Mat,IPIV=MeshMap%LoadLn2_A_Mat_piv, &
                        B=MeshMap%LoadLn2_M, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
            
      ! Transfer moments to the mesh fields
   
      Dest%Moment =  RESHAPE( MeshMap%LoadLn2_M, (/ 3, Dest%Nnodes /) )

      
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN        
      
   END IF ! Moment
   

END SUBROUTINE Convert_Point_To_Line2_Loads    
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the linearized matrices that convert the lumped (point) loads on nodes of a (line2) mesh 
!! to loads distributed across the line2 elements. (i.e., inverse lumping)
SUBROUTINE Linearize_Convert_Point_To_Line2_Loads(Dest, MeshMap, ErrStat, ErrMsg, DestDisp)

   TYPE(MeshType),                 INTENT(IN   )  :: Dest                           !< The mesh (on linearization input the loads are distributed along the line2 element)
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap                        !< The mapping data
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp                       !< The mesh that contains the translation displacement for these loads

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                        !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                         !< Error message if ErrStat /= ErrID_None

   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   character(*), parameter                        :: RoutineName = 'Linearize_Convert_Point_To_Line2_Loads'
      
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      !> Matrix \f$ \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} \f$, stored in modmesh_mapping::meshmaplinearizationtype::li,
      !! is allocated to be size Dest\%NNodes*3, Dest\%NNodes*3.
      !! This is the block matrix that maps the individual load fields of the source mesh to the corresponding field in the destination mesh.
      !! This matrix also forms part of the terms of modmesh_mapping::meshmaplinearizationtype::m_uS and modmesh_mapping::meshmaplinearizationtype::m_f.
      !! > However, since this is easier done with a solve, we are not actually going to form the matrix here.
      
      ! Now I need to get the cross-product terms for the moments as well:         
      ! M_uD and m_f for moments:
         
   if (Dest%FieldMask(MASKID_MOMENT) .AND. Dest%FieldMask(MASKID_FORCE) ) then               
            
      !> Matrix \f$ M_{uDm}^{DLinv} = - \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} M_{uDm}^{DL} \f$, stored in modmesh_mapping::meshmaplinearizationtype::m_us,
      !! is allocated to be size Dest\%NNodes*3, Dest\%NNodes*3.
      !! This is the block matrix that relates the destination translational displacement field to the moment field.
      !! > However, since this is easier done with a solve, we are not actually going to form the matrix here. Instead, I'll return \f$ M_{uDm}^{DL} \f$ 
            
      
      !> Matrix \f$ M_{fm}^{DLinv} = - \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1} M_{fm}^{DL} \begin{bmatrix} M_{li}^{DL} \end{bmatrix}^{-1}\f$, 
      !! stored in modmesh_mapping::meshmaplinearizationtype::m_f,
      !! is allocated to be size Dest\%NNodes*3, Dest\%NNodes*3.
      !! This is the block matrix that relates the force field to the moment field.
      !! > However, since this is easier done with a solve, we are not actually going to form the matrix here. Instead, I'll return \f$ M_{fm}^{DL} \f$ 
      
            
         ! this call sets MeshMap%dM%m_uD = M_{uDm}^{DL} and MeshMap%dM%m_f = M_{fm}^{DL}
      call FormMatrix_Lump_Line2_to_Point( Dest, Meshmap%dM, ErrStat2, ErrMsg2, DestDisp )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
                                                   
   end if
   

END SUBROUTINE Linearize_Convert_Point_To_Line2_Loads    
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Create_Augmented_Ln2_Src_Mesh(Src, Dest, MeshMap, Dest_TYPE, ErrStat, ErrMsg)
!This routine creates a new line2 mesh. This augmented line2-element source mesh is formed by splitting the original Line2-element 
! source mesh at each location where a destination mesh node (either a point or line2 element) projects orthogonally onto the 
! line2-element source mesh.
! It adds fields for forces, moments, and/or TranslationDisp, if they are part of the Src mesh

! creates MeshMap%Augmented_Ln2_Src and
! MeshMap%MapSrcToAugmt
! and (reallocates) MeshMap%MapLoads if necessary
! 
! bjj: maybe this should be in the AllocMapping routine? Or maybe the allocMapping routine shouldn't actually allocate arrays; do
!   it in the "IF (RemapFlag)" sections so that if people add nodes during the simulation, the structures get reallocated to correct 
!   size? allocMapping should maybe be MeshMapping_Init() and only check that fields are compatible, etc.
! 
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            ! The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap
   INTEGER(IntKi),                 INTENT(IN   )  :: Dest_TYPE                       ! Type of Dest elements to map

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                 :: iElem, i            ! do-loop counter for nodes/elements on source
   INTEGER(IntKi)                                 :: jElem, jNode, j     ! do-loop counter for nodes/elements on destination
   
   INTEGER(IntKi)                                 :: max_new_nodes, max_nodes, Aug_Nnodes
   INTEGER(IntKi)                                 :: Aug_NElem, curr_Aug_NElem
   INTEGER(IntKi)                                 :: n1, n2
   REAL(ReKi)                                     :: p_ED(3), p_ES(3), n1S_nD_vector(3), position(3)
   REAL(ReKi)                                     :: p_ED_orig(3), denom_orig
   REAL(R8Ki)                                     :: RefOrientation(3,3)
   REAL(DbKi)                                     :: TmpVec(3), RefOrientationD(3,3), FieldValue(3,2)   ! values for interpolating direction cosine matrices
   REAL(ReKi)                                     :: denom, elem_position
   REAL(ReKi), PARAMETER                          :: TOL = sqrt(epsilon(elem_position))  ! we're not using EqualRealNos here because we don't want elements of zero length (EqualRealNos produces elements of zero length)
   REAL(ReKi)                                     :: L         ! length of newly created element(s)
   
   TYPE(MeshType)                                 :: Temp_Ln2_Src                   ! A temporary mesh
   INTEGER, allocatable                           :: Original_Src_Element(:)
   REAL(ReKi), allocatable                        :: shape_fn2(:)
   
   
   
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Create_Augmented_Ln2_Src_Mesh'
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   ! We're going to find out the maximum number of new nodes that we will add to this mesh and create a new mesh with those nodes.
   ! We will likely not actually need all those extra nodes, so we'll have to transfer them to another mesh to commit it later.
         
   ! first, we need to know how many (additional) nodes we might need:
   !   each node of each destination element could potentially split each element of the source mesh
   max_new_nodes = (NumNodes( Dest_TYPE ) * dest%ElemTable(Dest_TYPE)%nelem) * src%ElemTable(ELEMENT_LINE2)%nelem     ! max number of new nodes
   max_nodes     = max_new_nodes + Src%nnodes                                                                         ! max total number of nodes in new mesh
      
   ! create a temporary mesh that we can work with (add nodes and split elements):
   ! note that we don't have any fields, and we will never commit this mesh, either.
   
   CALL MeshCreate(   BlankMesh       = Temp_Ln2_Src                 &
                     ,IOS             = Src%IOS                      &
                     ,NNodes          = max_nodes                    &
                     ,ErrStat         = ErrStat2                     &
                     ,ErrMess         = ErrMsg2                      )  
   
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN      
      END IF
   
   
   CALL AllocAry( Original_Src_Element, src%ElemTable(ELEMENT_LINE2)%nelem+max_new_nodes, 'Original_Src_Element', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN      
      END IF
   CALL AllocAry( shape_fn2,            max_nodes, 'shape_fn2', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN      
      END IF
   
   
   ! we'll first add the nodes and Line2 elements from the original source mesh   
   do i = 1, src%nnodes

      CALL MeshPositionNode ( Mesh = Temp_Ln2_Src                 &
                              ,INode = i                          &
                              ,Pos = Src%Position(:,i)            &
                              ,Orient = Src%RefOrientation(:,:,i) &
                              ,ErrStat   = ErrStat2               &
                              ,ErrMess   = ErrMsg2                 )      

         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN      
         END IF
      
   !   Original_Src_Element(i)=-1  !Part of the original mesh
   enddo
   DO i=1,src%ElemTable(ELEMENT_LINE2)%nelem
      CALL MeshConstructElement ( Mesh = Temp_Ln2_Src           &
                                 ,Xelement = ELEMENT_LINE2      &
                                 ,P1       = Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)  &
                                 ,P2       = Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)  &
                                 ,ErrStat  = ErrStat2            &
                                 ,ErrMess  = ErrMsg2             )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN      
         END IF
      
      Original_Src_Element(i)=i  !Part of the original mesh      
   end do
   
   ! now we'll augment the mesh:
   Aug_Nnodes = Src%nnodes  ! number of nodes in the augmented mesh 
   Aug_NElem  = Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%nelem

   ! loop through the destination elements (NOTE: for point elements, this is the same as looping over nodes):
   DO jElem = 1,dest%ElemTable(Dest_TYPE)%nelem
      
      IF ( Dest_TYPE == ELEMENT_LINE2 ) THEN
         p_eD =   dest%Position(:, dest%ElemTable(Dest_TYPE)%Elements(jElem)%ElemNodes(2)) &
                - dest%Position(:, dest%ElemTable(Dest_TYPE)%Elements(jElem)%ElemNodes(1))
         p_eD_orig = p_eD ! save for later calculations (to allow point elements, too)
      END IF
      
     
      iElem = 1         
      curr_Aug_NElem = Aug_NElem
      j = 1   
      
      Src_Elements: DO WHILE ( iElem <= curr_Aug_NElem )  ! bjj: we're adding elements and nodes to Temp_Ln2_Src, so this must be flexible
         
         p_eS =   Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2)) &
                - Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1))
         
         IF ( Dest_TYPE == ELEMENT_POINT ) p_eD = p_eS
         denom = DOT_PRODUCT( p_eD , p_eS )
            
  
         IF ( .NOT. EqualRealNos( denom, 0.0_ReKi) ) THEN ! we ignore source elements that are parallel to the destination element (i.e., denom == 0)
            DO jNode = j, NumNodes( Dest_TYPE ) ! check each node of the destination element
               n1S_nD_vector =            dest%Position(:, dest%ElemTable(Dest_TYPE)%Elements(jElem)%ElemNodes(jNode)) &
                                - Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1))
               
               elem_position = DOT_PRODUCT( p_eD, n1S_nD_vector ) / denom                 ! Eq 37 (AIAA 2015 paper)
               
!bjj: todo: we need to set this TOL based on actual distances, not relative values (0,1) on an element....                 
!  for now, I've calculated the element length inside this tolerance and reserve the right to reject new nodes that create 0-length elements.
               
               ! if 0 < elem_position < 1, we create a new node and split this element
               IF ( elem_position > TOL .AND. elem_position < (1.0_ReKi - TOL)  ) THEN
                  
                  
                     ! Add a node (and therefore an element):
                  Aug_Nnodes = Aug_Nnodes + 1

                     ! calculate the position and orientation relative to the *original* source element:
                  n1=Src%ElemTable(ELEMENT_LINE2)%Elements(Original_Src_Element(iElem))%ElemNodes(1)
                  n2=Src%ElemTable(ELEMENT_LINE2)%Elements(Original_Src_Element(iElem))%ElemNodes(2)
                                          
                                          
                  p_eS  = Src%Position(:, n2) - Src%Position(:, n1)
                  IF ( Dest_TYPE == ELEMENT_POINT ) p_eD_orig = p_eS
                  denom_orig = DOT_PRODUCT( p_eD_orig , p_eS )   ! we don't need to check that this is zero because it's just a shorter version of the temp Temp_Ln2_Src element
                  n1S_nD_vector =   dest%Position(:, dest%ElemTable(Dest_TYPE)%Elements(jElem)%ElemNodes(jNode)) &
                                   - Src%Position(:, n1 )
                  shape_fn2(Aug_Nnodes) = DOT_PRODUCT( p_eD_orig, n1S_nD_vector ) / denom_orig       ! save this for later, when we need to map the mesh fields...
                  
                     ! interpolate position on the original souce element:                     
                     
                  position = (1.0_ReKi - shape_fn2(Aug_Nnodes)) * Src%Position(:, n1) &
                                       + shape_fn2(Aug_Nnodes)  * Src%Position(:, n2) 
                  
                  ! let's just verify that this new node (n1) doesn't give us zero-length elements:
                  ! (note we use the NEW (not original) source element, which may have been split)
                  p_eS = position - Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1) )
                  L = SQRT(dot_product(p_eS,p_eS)) ! length of new element
                  
                  IF ( L < MIN_LINE2_ELEMENT_LENGTH ) THEN ! this element is basically zero length
                        ! for numerical reasons, we really didn't want this node....
                     Aug_Nnodes = Aug_Nnodes - 1
                  ELSE
                     
                     ! let's verify the other node (n2) of this element doesn't give zero-length:
                     p_eS = position - Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2))
                     L = SQRT(dot_product(p_eS,p_eS)) ! length of new element
                     
                     IF ( L < MIN_LINE2_ELEMENT_LENGTH ) THEN ! this element is basically zero length
                        ! for numerical reasons, we really didn't want this node....
                        Aug_Nnodes = Aug_Nnodes - 1
                     ELSE   
   
                           ! we can add the node (and an element)
                        Aug_NElem  = Aug_NElem + 1 
                        CALL MeshSplitElement_2PT( Temp_Ln2_Src, ELEMENT_LINE2, ErrStat2, ErrMsg2, iElem, Aug_Nnodes  )                      
                           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                           IF (ErrStat >= AbortErrLev) THEN
                              CALL CleanUp()
                              RETURN      
                           END IF
                       
                        Original_Src_Element( Aug_NElem ) = Original_Src_Element( iElem )  ! this node can now be mapped to original source mesh element               
                       
                        
                     ! get the Reference orientation for this new node
                        ! convert DCMs to tensors: 

                        RefOrientationD = Src%RefOrientation(:, :, n1)
                        CALL DCM_logmap( RefOrientationD, FieldValue(:,1), ErrStat, ErrMsg )
                        IF (ErrStat >= AbortErrLev) RETURN
                  
                        RefOrientationD = Src%RefOrientation(:, :, n2)
                        CALL DCM_logmap( RefOrientationD, FieldValue(:,2), ErrStat, ErrMsg )
                        IF (ErrStat >= AbortErrLev) RETURN
         
                        CALL DCM_SetLogMapForInterp( FieldValue )  ! make sure we don't cross a 2pi boundary
                  
                           ! interpolate tensors: 
                        TmpVec = (1.0_ReKi - shape_fn2(Aug_Nnodes)) * FieldValue(:, 1) &
                                             + shape_fn2(Aug_Nnodes)  * FieldValue(:, 2) 
                              
                           ! convert back to DCM:
                        RefOrientationD = DCM_exp( TmpVec )
                        RefOrientation  = REAL(RefOrientationD, R8Ki)

                        
                        CALL MeshPositionNode ( Mesh       = Temp_Ln2_Src                      &
                                                ,INode     = Aug_Nnodes                        &
                                                ,Pos       = position                          & 
                                                ,Orient    = RefOrientation                    &
                                                ,ErrStat   = ErrStat2                          &
                                                ,ErrMess   = ErrMsg2                           )
                           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                           IF (ErrStat >= AbortErrLev) THEN
                              CALL CleanUp()
                              RETURN
                           END IF
                     
                        ! if we have to check a second node, we need to first recalculate p_eS and denom on Temp_Ln2_Src:
                        IF ( jNode < NumNodes( Dest_TYPE )) THEN
                           j = jNode+1 ! start on the next node
                           CYCLE Src_Elements 
                        END IF
                        
                     END IF ! doesn't cause zero-length element
                     
                     
                  END IF
                  
               END IF ! New node/element
               
            END DO !jNode
            j = 1
               
         END IF ! elements aren't parallel
         
            ! move to the next source element (i.e., CYCLE Src_Elements)
         iElem = iElem + 1
            
      END DO Src_Elements ! while
                  
   END DO !jElem
   
            
   !..................      
   ! Now let's create the actual augmented mesh (should be smaller than the Temp_Ln2_Src mesh):   
   !.................
   
   CALL MeshDestroy( MeshMap%Augmented_Ln2_Src, ErrStat2, ErrMsg2, .TRUE. )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN      
      END IF
   
   
   CALL MeshCreate(   BlankMesh       = MeshMap%Augmented_Ln2_Src              &
                     ,IOS             = Src%IOS                                &
                     ,NNodes          = Aug_Nnodes                             &
                     ,Force           = Src%FieldMask(maskid_force)            &
                     ,Moment          = Src%FieldMask(maskid_moment)           & !                     ,TranslationDisp = Src%FieldMask(maskid_TranslationDisp)  &
                     ,TranslationDisp = .TRUE.                                 &
                     ,ErrStat         = ErrStat2                               &
                     ,ErrMess         = ErrMsg2                                )

      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN      
      END IF   
   
   do i = 1, Aug_Nnodes
            
      CALL MeshPositionNode ( Mesh = MeshMap%Augmented_Ln2_Src                &
                              ,INode = i                                      &
                              ,Pos = Temp_Ln2_Src%Position(:,i)               &
                              ,Orient = Temp_Ln2_Src%RefOrientation(:,:,i)    &
                              ,ErrStat   = ErrStat2                           &
                              ,ErrMess   = ErrMsg2                            )      
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN      
         END IF      

   enddo
   DO i=1,Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%nelem
      CALL MeshConstructElement ( Mesh = MeshMap%Augmented_Ln2_Src                                           &
                                 ,Xelement = ELEMENT_LINE2                                                   &
                                 ,P1       = Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)  &
                                 ,P2       = Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)  &
                                 ,ErrStat  = ErrStat2                                                        &
                                 ,ErrMess  = ErrMsg2                                                         )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN      
         END IF   
   end do   
   
   CALL MeshCommit ( MeshMap%Augmented_Ln2_Src, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN      
      END IF   
   
   
      ! we're going to need the mapping data from source to augmented source:
   IF (Aug_Nnodes > Src%Nnodes) THEN
      IF (ALLOCATED(MeshMap%MapSrcToAugmt)) THEN
         IF (UBOUND(MeshMap%MapSrcToAugmt,1) < Aug_NNodes) THEN
            DEALLOCATE (MeshMap%MapSrcToAugmt)
         END IF
      END IF
      
      IF (.NOT. ALLOCATED(MeshMap%MapSrcToAugmt)) THEN
         ALLOCATE( MeshMap%MapSrcToAugmt((Src%Nnodes+1):Aug_NNodes), STAT=ErrStat2 ) 
         
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate MeshMap%MapSrcToAugmt.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN      
         END IF
      END IF
         
      j=Src%ElemTable(ELEMENT_LINE2)%nelem 
      do i=LBOUND(MeshMap%MapSrcToAugmt,1),UBOUND(MeshMap%MapSrcToAugmt,1)
         j = j+1
         MeshMap%MapSrcToAugmt(i)%OtherMesh_Element = Original_Src_Element( j ) ! we added just as many nodes as elements...
         MeshMap%MapSrcToAugmt(i)%shape_fn(2)       = shape_fn2( i )
         MeshMap%MapSrcToAugmt(i)%shape_fn(1)       = 1.0_ReKi - shape_fn2( i )
      end do
            
      IF (ALLOCATED(MeshMap%DisplacedPosition)) THEN
         IF ( SIZE(MeshMap%DisplacedPosition,2) < Aug_Nnodes )  THEN
         
            DEALLOCATE (MeshMap%DisplacedPosition)
         END IF
      END IF
      
      IF (.NOT. ALLOCATED(MeshMap%DisplacedPosition)) THEN            
         CALL AllocAry( MeshMap%DisplacedPosition, 3, Aug_Nnodes,2, 'MeshMap%DisplacedPosition', ErrStat2, ErrMsg2)               
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN      
            END IF
      END IF

      
      IF ( ALLOCATED(MeshMap%MapLoads)) THEN
         IF (SIZE(MeshMap%MapLoads) < Aug_NNodes) THEN
            DEALLOCATE( MeshMap%MapLoads )
         END IF
      END IF
      
      IF (.NOT. ALLOCATED(MeshMap%MapLoads)) THEN            
         ALLOCATE( MeshMap%MapLoads(Aug_Nnodes), STAT=ErrStat2 )
         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate MeshMap%MapLoads.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN      
         END IF      
      END IF
      
      
   END IF
   
      
   CALL CleanUp()
   RETURN
   
      ! clean up
CONTAINS
   SUBROUTINE CleanUp()
   
      CALL MeshDestroy( Temp_Ln2_Src, ErrStat2, ErrMsg2 )
      IF ( ALLOCATED(Original_Src_Element) ) DEALLOCATE(Original_Src_Element)
      IF ( ALLOCATED(shape_fn2)            ) DEALLOCATE(shape_fn2)
      
   END SUBROUTINE CleanUp
   
END SUBROUTINE Create_Augmented_Ln2_Src_Mesh
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine takes the loads from a line2 mesh and transfers them to an augmented line2 mesh (the src mesh with extra/augmented nodes).
!! It also transfers the translation displacement field to this augmented mesh so that we have only one mesh to deal with
!! instead of two. The first NNodes values of the two meshes are equal, and the additional, augmented nodes are split between
!! nodes on a line2 element. 
SUBROUTINE Transfer_Src_To_Augmented_Ln2_Src( Src, MeshMap, ErrStat, ErrMsg, SrcDisp, LoadsScaleFactor )
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Src                         !< The source mesh containing line2 loads
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap                     !< mesh mapping data, including the augmented source mesh

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                      !< Error message if ErrStat /= ErrID_None
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp                     !< The displacements associated with the source mesh
   REAL(ReKi),                     INTENT(IN)     :: LoadsScaleFactor            !< Scaling factor for loads (to help with numerical issues)

      ! local variables
   INTEGER(IntKi)                                 :: iElem, i  ! do-loop counter for nodes/elements on source   
   INTEGER(IntKi)                                 :: n1, n2
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !.......................................
   ! TranslationDisp (from external mesh to augmented one)
   !.......................................
   
   IF ( SrcDisp%FIELDMASK(MASKID_TranslationDisp) ) THEN      
   
      DO i = 1,Src%NNodes  !in case SrcDisp has more nodes that Src (e.g. ElastoDyn blade and tower meshes), I'm using Src%NNodes here
         MeshMap%Augmented_Ln2_Src%TranslationDisp(:,i) = SrcDisp%TranslationDisp(:,i)
      END DO
      
      DO i = (Src%NNodes+1),MeshMap%Augmented_Ln2_Src%NNodes
         iElem = MeshMap%MapSrcToAugmt(i)%OtherMesh_Element
         
         n1=SrcDisp%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1)
         n2=SrcDisp%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2)
         
         MeshMap%Augmented_Ln2_Src%TranslationDisp(:,i) = MeshMap%MapSrcToAugmt(i)%shape_fn(1) * SrcDisp%TranslationDisp(:, n1) &
                                                        + MeshMap%MapSrcToAugmt(i)%shape_fn(2) * SrcDisp%TranslationDisp(:, n2)   
                                                                                                            
      END DO
      
   END IF
         
   !.......................................
   ! Force
   !.......................................
   
   IF ( Src%FIELDMASK(MASKID_Force) ) THEN      
   
      DO i = 1,Src%NNodes
         MeshMap%Augmented_Ln2_Src%Force(:,i) = Src%Force(:,i)
      END DO
      
      DO i = (Src%NNodes+1),MeshMap%Augmented_Ln2_Src%NNodes
         iElem = MeshMap%MapSrcToAugmt(i)%OtherMesh_Element
         
         n1=Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1)
         n2=Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2)
         
         MeshMap%Augmented_Ln2_Src%Force(:,i) = MeshMap%MapSrcToAugmt(i)%shape_fn(1) * (Src%Force(:, n1)/LoadsScaleFactor) &
                                              + MeshMap%MapSrcToAugmt(i)%shape_fn(2) * (Src%Force(:, n2)/LoadsScaleFactor)   
         
         MeshMap%Augmented_Ln2_Src%Force(:,i) = MeshMap%Augmented_Ln2_Src%Force(:,i)*LoadsScaleFactor
                                                                                                            
      END DO
      
   END IF
   
   
   !.......................................
   ! Moment
   !.......................................
   
   IF ( Src%FIELDMASK(MASKID_Moment) ) THEN      
   
      DO i = 1,Src%NNodes
         MeshMap%Augmented_Ln2_Src%Moment(:,i) = Src%Moment(:,i)
      END DO
      
      DO i = (Src%NNodes+1),MeshMap%Augmented_Ln2_Src%NNodes
         iElem = MeshMap%MapSrcToAugmt(i)%OtherMesh_Element
         
         n1=Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1)
         n2=Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2)
         
         MeshMap%Augmented_Ln2_Src%Moment(:,i) = MeshMap%MapSrcToAugmt(i)%shape_fn(1) * (Src%Moment(:, n1)/LoadsScaleFactor) &
                                               + MeshMap%MapSrcToAugmt(i)%shape_fn(2) * (Src%Moment(:, n2)/LoadsScaleFactor)   
         
         MeshMap%Augmented_Ln2_Src%Moment(:,i) = MeshMap%Augmented_Ln2_Src%Moment(:,i)*LoadsScaleFactor
         
                                                                                                            
      END DO
      
   END IF   
   
   
         
END SUBROUTINE Transfer_Src_To_Augmented_Ln2_Src
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine linearizes Transfer_Src_To_Augmented_Ln2_Src (modmesh_mapping::transfer_src_to_augmented_ln2_src). 
SUBROUTINE Linearize_Src_To_Augmented_Ln2_Src( Src, MeshMap, ErrStat, ErrMsg, SrcDisp )
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Src                         !< The source mesh containing line2 loads
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap                     !< mesh mapping data, including the augmented source mesh

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                      !< Error message if ErrStat /= ErrID_None
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp                     !< The displacements associated with the source mesh

      ! local variables
   INTEGER(IntKi)                                 :: iElem, i  ! do-loop counter for nodes/elements on source   
   INTEGER(IntKi)                                 :: n, n1, n2, j
   
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   CHARACTER(*), PARAMETER                        :: RoutineName = 'Linearize_Src_To_Augmented_Ln2_Src'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
    
   !> Matrix \f$ M_{li}^A = M_{mi}^A = M^A \f$, stored in modmesh_mapping::meshmaplinearizationtype::li,
   !! is allocated to be size MeshMap\%Augmented_Ln2_Src\%NNodes*3, Src\%NNodes*3.
   !> This is the matrix that maps each field of the source mesh to an augmented mesh.
         
   if (.not. allocated(MeshMap%dM%li) ) then
      call AllocAry(MeshMap%dM%li, MeshMap%Augmented_Ln2_Src%NNodes*3, Src%Nnodes*3, 'dM%li', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
   elseif (size(MeshMap%dM%li,1) /= MeshMap%Augmented_Ln2_Src%NNodes*3 .or. size(MeshMap%dM%li,2) /= Src%Nnodes*3 ) then
         ! this would happen if we've called the linearization routine multiple times
      deallocate(MeshMap%dM%li)
      call AllocAry(MeshMap%dM%li, MeshMap%Augmented_Ln2_Src%NNodes*3, Src%Nnodes*3, 'dM%li', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN       
   end if
      
   MeshMap%dM%li = 0.0_R8Ki      
   do i=1, Src%Nnodes
      do j=1,3
         MeshMap%dM%li( (i-1)*3+j, (i-1)*3+j ) = 1.0_R8Ki
      end do      
   end do
      
   do i = (Src%NNodes+1),MeshMap%Augmented_Ln2_Src%NNodes
         
      iElem = MeshMap%MapSrcToAugmt(i)%OtherMesh_Element                  
      do n1=1,NumNodes(ELEMENT_LINE2)
            
         n = SrcDisp%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(n1)
            
         do j=1,3
            MeshMap%dM%li( (i-1)*3+j, (n-1)*3+j ) = MeshMap%MapSrcToAugmt(i)%shape_fn(n1)
         end do   
            
      end do
         
   end do
   
      
   !.......................................
   ! TranslationDisp (from external mesh to augmented one)
   !.......................................
   
   IF ( SrcDisp%FIELDMASK(MASKID_TranslationDisp) ) THEN      
   
      DO i = 1,Src%NNodes  !in case SrcDisp has more nodes that Src (e.g. ElastoDyn blade and tower meshes), I'm using Src%NNodes here
         MeshMap%Augmented_Ln2_Src%TranslationDisp(:,i) = SrcDisp%TranslationDisp(:,i)
      END DO
      
      DO i = (Src%NNodes+1),MeshMap%Augmented_Ln2_Src%NNodes
         iElem = MeshMap%MapSrcToAugmt(i)%OtherMesh_Element
         
         n1=SrcDisp%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1)
         n2=SrcDisp%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2)
         
         MeshMap%Augmented_Ln2_Src%TranslationDisp(:,i) = MeshMap%MapSrcToAugmt(i)%shape_fn(1) * SrcDisp%TranslationDisp(:, n1) &
                                                        + MeshMap%MapSrcToAugmt(i)%shape_fn(2) * SrcDisp%TranslationDisp(:, n2)   
                                                                                                            
      END DO
      
   END IF   
   
   !.......................................
   ! Force
   !> We need the distributed forces on the augmented mesh, so we transfer them in the linearization, too.
   !! Note that this is done without the LoadsScaleFactor.
   !.......................................
   
   IF ( Src%FIELDMASK(MASKID_Force) ) THEN      
   
      DO i = 1,Src%NNodes
         MeshMap%Augmented_Ln2_Src%Force(:,i) = Src%Force(:,i)
      END DO
      
      DO i = (Src%NNodes+1),MeshMap%Augmented_Ln2_Src%NNodes
         iElem = MeshMap%MapSrcToAugmt(i)%OtherMesh_Element
         
         n1=Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1)
         n2=Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2)
         
         MeshMap%Augmented_Ln2_Src%Force(:,i) = MeshMap%MapSrcToAugmt(i)%shape_fn(1) * Src%Force(:, n1) &
                                              + MeshMap%MapSrcToAugmt(i)%shape_fn(2) * Src%Force(:, n2)   
                                                                                                                     
      END DO
      
   END IF   
                                 
END SUBROUTINE Linearize_Src_To_Augmented_Ln2_Src
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine creates a matrix used to "unlump" loads later (needs to have element connectivity information to create it)
SUBROUTINE Create_InverseLumping_Matrix( Dest, MeshMap, ErrStat, ErrMsg )
   TYPE(MeshType),         INTENT(IN   ) ::  Dest        !< destination mesh
   TYPE(MeshMapType),      INTENT(INOUT) ::  MeshMap     !< mesh mapping data


   INTEGER(IntKi),         INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(ReKi)                              :: c, TwoC

   INTEGER(IntKi)                          :: N, n1, n2  ! node numbers
   INTEGER(IntKi)                          :: j,k, iElem, iComp
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'Create_InverseLumping_Matrix'


   ErrStat = ErrID_None
   ErrMsg = ""

   n = Dest%Nnodes* 3

   IF ( ALLOCATED(MeshMap%LoadLn2_A_Mat) ) THEN
      IF ( SIZE(MeshMap%LoadLn2_A_Mat,1) < n ) THEN
         DEALLOCATE(MeshMap%LoadLn2_A_Mat)
         DEALLOCATE(MeshMap%LoadLn2_A_Mat_piv)
         DEALLOCATE(MeshMap%LoadLn2_F)
         DEALLOCATE(MeshMap%LoadLn2_M)
         CALL AllocInvLumpingArrays()
         IF (ErrStat >= AbortErrLev) RETURN 
      END IF
   ELSE
      CALL AllocInvLumpingArrays()  
      IF (ErrStat >= AbortErrLev) RETURN 
   END IF
        
   
      ! factor the "A" matrix for the solve later:
         
   MeshMap%LoadLn2_A_Mat = 0.0_ReKi
                           
      ! loop over source mesh, integrating over each element
   do iElem = 1, Dest%ElemTable(ELEMENT_LINE2)%nelem

      n1 = Dest%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1)
      n2 = Dest%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2)

      c    = Dest%ElemTable(ELEMENT_LINE2)%Elements(iElem)%det_jac / 3.0_ReKi  != TwoNorm( p(n2)-p(n1) )/6
      TwoC = 2.0_ReKi * c

      do iComp=1,3 !3 components of force
         j = (n1-1)*3 + iComp      ! this is the index of the node/component of the lumped value
         k = (n2-1)*3 + iComp      ! the is the portion of the "lesser" node
         
         MeshMap%LoadLn2_A_Mat(j,j) = MeshMap%LoadLn2_A_Mat(j,j) + TwoC
         MeshMap%LoadLn2_A_Mat(k,k) = MeshMap%LoadLn2_A_Mat(k,k) + TwoC
         MeshMap%LoadLn2_A_Mat(j,k) = MeshMap%LoadLn2_A_Mat(j,k) + c
         MeshMap%LoadLn2_A_Mat(k,j) = MeshMap%LoadLn2_A_Mat(k,j) + c
      enddo !icomp
   enddo !i
      
   CALL LAPACK_getrf(n,n,MeshMap%LoadLn2_A_Mat,MeshMap%LoadLn2_A_Mat_piv, ErrStat, ErrMsg)      
           
!........................................................................................................   
CONTAINS
   SUBROUTINE AllocInvLumpingArrays
   
         CALL AllocAry( MeshMap%LoadLn2_A_Mat,     n, n,  'MeshMap%LoadLn2_A_Mat',    ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         CALL AllocAry( MeshMap%LoadLn2_A_Mat_piv, n,     'MeshMap%LoadLn2_A_Mat_piv',ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         CALL AllocAry( MeshMap%LoadLn2_F,         n, 1,  'MeshMap%LoadLn2_F',        ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         CALL AllocAry( MeshMap%LoadLn2_M,         n, 1,  'MeshMap%LoadLn2_M',        ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   END SUBROUTINE AllocInvLumpingArrays
   
END SUBROUTINE Create_InverseLumping_Matrix
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine creates the mapping from a line2 mesh with loads to another line2 mesh.
SUBROUTINE CreateLoadMap_L2_to_L2( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             !< The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            !< The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap                         !< mapping data

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          !< Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   character(*), parameter                        :: RoutineName = 'CreateLoadMap_L2_to_L2'

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      !........................
      !> + Create a matrix used to "unlump" loads later (needs to have element connectivity information to create it)
      !........................
      CALL Create_InverseLumping_Matrix( Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
            
      !........................
      !> + An augmented Line2-element source mesh is formed
      !! by splitting the original Line2-element source mesh at each location where a
      !! destination-mesh Line2-element node projects orthogonally from the destination mesh
      !........................
         
      CALL Create_Augmented_Ln2_Src_Mesh(Src, Dest, MeshMap, ELEMENT_LINE2, ErrStat2, ErrMsg2) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
                     
      !........................
      !> + For each Line2-element node of the augmented source mesh, 
      !! a nearest-neighbor Line2 element of the destination mesh is found (else aborted) in the reference configuration, 
      !! for which the source Line2-element node projects orthogonally onto the destination Line2-element domain.
      !! A destination-mesh Line2 element may be associated with multiple source-mesh Line2-element nodes.
      !........................
            
      CALL CreateMapping_ProjectToLine2(MeshMap%Augmented_Ln2_Src, Dest, MeshMap%MapLoads, ELEMENT_LINE2, ErrStat2, ErrMsg2)            
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
         
         
      !........................
      !> + Create a temporary mesh for lumped point elements of the line2 source 
      !........................
      CALL Create_PointMesh( MeshMap%Augmented_Ln2_Src, MeshMap%Lumped_Points_Src, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
         
         
END SUBROUTINE CreateLoadMap_L2_to_L2        
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine creates the mapping of motions between two meshes
SUBROUTINE CreateMotionMap_L2_to_L2( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             !< The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            !< The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap                         !< structure that contains data necessary to map these two meshes

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          !< Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   CALL CreateMapping_ProjectToLine2(Dest,Src, MeshMap%MapMotions, ELEMENT_LINE2, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateMotionMap_L2_to_L2')
      IF (ErrStat >= AbortErrLev) RETURN
            
END SUBROUTINE CreateMotionMap_L2_to_L2        
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine takes a Line2 mesh (SRC) with distributed forces and/or moments and lumps those to the concentrated forces
!! and/or moments at the endpoints (nodes) of the Line2 mesh.   These are placed in a Point mesh (DEST) where the Points are
!! assumed to be colocated with the nodes of the Line2 Mesh.
!!
!! This routine is expected to be used in the transfer of HydroDyn Output (distributed force/moment) to the SubDyn Input
!! (concentrated force/moment); routine is to be called from the Input_Ouput routine for HydroDyn to SubDyn.
!!
!! Formulation:
!! Distrbuted quantities are integrated over the element with two-point nodal (Gauss-Lobatto) quadrature.
!!
!! F(x_j) = int_{-1}^{1} f(xi) phi_j(xi) J dxi (exact)
!!
!! F(x_1) =~ f(-1) phi_1(-1) J w_1 + f(+1) phi_1(+1) J w_2 (two-point quadrature)
!! F(x_2) =~ f(-1) phi_2(-1) J w_1 + f(+1) phi_2(+1) J w_2 (two-point quadrature)
!!
!! which can be simplified to
!!
!! F(x_1) =~ f(-1) J
!! F(x_2) =~ f(+1) J
!!
!! since w_1 = w_2 = 1, phi_j(xi_i) = Delta_{j,i}  (Kronecker Delta)
!!
!! where x_j is the jth node in 3D physical coordinates, xi_i is the ith node in 1D element coordinates, j in {1,2},
!! xi in [-1,1] is the element natural coordinate, f(xi) is the distributed quantity, phi_j(xi) is the basis function
!! for the jth node, J is the determinant of the Jacobian
!! in the mapping from [X_1,X_2] to [-1,1] (X_j is the 3D location of the jth point)
!!
!! For Line2 elements, J = DIST(X_1,X_2) / 2
!!
SUBROUTINE Lump_Line2_to_Point( Line2_Src, Point_Dest, ErrStat, ErrMsg, SrcDisp, LoadsScaleFactor )

   TYPE(MeshType),         INTENT(IN   ) ::  Line2_Src   !< line2 source mesh
   TYPE(MeshType),         INTENT(INOUT) ::  Point_Dest  !< point destination mesh
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  SrcDisp     !< Another mesh that may contain the source's TranslationDisp field

   REAL(ReKi),             INTENT(IN)    :: LoadsScaleFactor  !< Scaling factor for loads (to help with numerical issues)

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   
   ! local variables
   REAL(ReKi) :: det_jac_Ovr3
!  REAL(ReKi) :: n1_n2_vector(3) ! vector going from node 1 to node 2 in a Line2 element
   REAL(ReKi) :: pCrossf(3)      ! a temporary vector storing cross product of positions with forces
   REAL(ReKi) :: p1(3), p2(3)    ! temporary position vectors
   REAL(ReKi) :: dp(3)           ! a temporary vector storing p2-p1

   INTEGER(IntKi) :: i
   INTEGER(IntKi) :: nnodes

   INTEGER(IntKi) :: n1
   INTEGER(IntKi) :: n2
   character(*), parameter               :: RoutineName = 'Lump_Line2_to_Point'


   ErrStat = ErrID_None
   ErrMsg = ""

#ifdef MESH_DEBUG     
   ! bjj: we shouldn't have to check this in production mode; this routine is an internal one only.
   
   if (Line2_Src%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Source mesh must have one or more Line2 Elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif

   if (Point_Dest%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Destination mesh must have one or more Point Elements.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif

   if (Line2_Src%nnodes .ne. Point_Dest%nnodes) then
      CALL SetErrStat( ErrID_Fatal, 'Source and Destination meshes must have same number of Nodes.', ErrStat, ErrMsg, RoutineName)
      RETURN
   endif

   if (Point_Dest%FieldMask(MASKID_FORCE) ) then
      if (.not. Line2_Src%FieldMask(MASKID_FORCE) ) then
         CALL SetErrStat( ErrID_Fatal, 'Destination mesh contains Force, but Source does not.', ErrStat, ErrMsg, RoutineName)
         RETURN
      endif
   endif
!bjj: this may not be covered... 
   if (Point_Dest%FieldMask(MASKID_MOMENT) ) then
      if (.not. Line2_SRC%FieldMask(MASKID_MOMENT) ) then
         CALL SetErrStat( ErrID_Fatal, 'Destination mesh contains Moment, but Source does not.', ErrStat, ErrMsg, RoutineName)
         RETURN
      endif
   endif
   
#endif

   ! see equations 40-43 of Spraque, Jonkman, and Jonkman [2014]
   
   nnodes = Point_Dest%nnodes  ! also = Line2_Src%nnodes

   if ( Point_Dest%FieldMask(MASKID_TRANSLATIONDISP) ) then
      Point_Dest%TranslationDisp = Line2_Src%TranslationDisp
   end if
   
   
   ! initialize force/moment in Dest
   if (Point_Dest%FieldMask(MASKID_FORCE) )  Point_Dest%Force  = 0.
   if (Point_Dest%FieldMask(MASKID_MOMENT) ) Point_Dest%Moment = 0.

   ! loop over source mesh, integrating over each element
   do i = 1, Line2_Src%ElemTable(ELEMENT_LINE2)%nelem

      n1 = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
      n2 = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)

      !n1_n2_vector = Line2_Src%Position(:,n2) - Line2_Src%Position(:,n1)
      !
      !det_jac_Ovr3 = SQRT( DOT_PRODUCT(n1_n2_vector,n1_n2_vector) ) / 6._ReKi  = L / 6 = det_jac/3.

      det_jac_Ovr3 = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%det_jac / 3.0_ReKi


      if (Point_Dest%FieldMask(MASKID_FORCE) ) then

         Point_Dest%Force(:,n1) = Point_Dest%Force(:,n1) + det_jac_Ovr3 * ( 2.*(Line2_Src%Force(:,n1)/LoadsScaleFactor) +    (Line2_Src%Force(:,n2)/LoadsScaleFactor) )
         Point_Dest%Force(:,n2) = Point_Dest%Force(:,n2) + det_jac_Ovr3 * (    (Line2_Src%Force(:,n1)/LoadsScaleFactor) + 2.*(Line2_Src%Force(:,n2)/LoadsScaleFactor) )

      endif

      if (Point_Dest%FieldMask(MASKID_MOMENT) ) then

         Point_Dest%Moment(:,n1) = Point_Dest%Moment(:,n1) + det_jac_Ovr3 * ( 2.*(Line2_Src%Moment(:,n1)/LoadsScaleFactor) +    (Line2_Src%Moment(:,n2)/LoadsScaleFactor) )
         Point_Dest%Moment(:,n2) = Point_Dest%Moment(:,n2) + det_jac_Ovr3 * (    (Line2_Src%Moment(:,n1)/LoadsScaleFactor) + 2.*(Line2_Src%Moment(:,n2)/LoadsScaleFactor) )


         if ( Point_Dest%FieldMask(MASKID_FORCE) ) then

            p1 = Line2_Src%Position(:,n1)
            p2 = Line2_Src%Position(:,n2)

            IF ( PRESENT(SrcDisp) ) THEN
               p1 = p1 + SrcDisp%TranslationDisp(:,n1)
               p2 = p2 + SrcDisp%TranslationDisp(:,n2)
            ELSEIF (Line2_Src%FieldMask(MASKID_TRANSLATIONDISP)) THEN
               p1 = p1 + Line2_Src%TranslationDisp(:,n1)
               p2 = p2 + Line2_Src%TranslationDisp(:,n2)
            END IF

            dp = p2-p1
            pCrossf = (Line2_Src%Force(:,n1)/LoadsScaleFactor) + (Line2_Src%Force(:,n2)/LoadsScaleFactor) !temp storage of f to avoid array creating in cross_product
            pCrossf = 0.5*det_jac_Ovr3 *cross_product( dp, pCrossf)

            Point_Dest%Moment(:,n1) = Point_Dest%Moment(:,n1) + pCrossf
            Point_Dest%Moment(:,n2) = Point_Dest%Moment(:,n2) - pCrossf

         end if ! src  moment AND force

      endif

   enddo

   if (Point_Dest%FieldMask(MASKID_FORCE) )  Point_Dest%Force  = Point_Dest%Force * LoadsScaleFactor
   if (Point_Dest%FieldMask(MASKID_MOMENT) ) Point_Dest%Moment = Point_Dest%Moment* LoadsScaleFactor
   

END SUBROUTINE Lump_Line2_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine linearizes the lumping of distributed loads to point loads.
SUBROUTINE Linearize_Lump_Line2_to_Point( Line2_Src, Point_Dest, dM, ErrStat, ErrMsg )

   TYPE(MeshType),         INTENT(IN   ) :: Line2_Src   !< line2 source mesh
   TYPE(MeshType),         INTENT(INOUT) :: Point_Dest  !< point destination mesh

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     !< Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      !< Error message if ErrStat /= ErrID_None

   type(MeshMapLinearizationType), intent(inout) :: dM  !< data structure for linearization

   
   ! local variables
   REAL(R8Ki)                            :: c, TwoC
  
   INTEGER(IntKi)                        :: n1, n2  ! node numbers
   INTEGER(IntKi)                        :: i, j, k, iComp  
   Integer(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   character(*),   parameter             :: RoutineName = 'Linearize_Lump_Line2_to_Point'


   ErrStat = ErrID_None
   ErrMsg = ""
         
      ! identity for forces and/or moments:

   !> Matrix \f$ M_{li}^{SL} \f$, stored in modmesh_mapping::meshmaplinearizationtype::li,
   !! is allocated to be size Point_Dest\%NNodes*3, Line2_Src\%NNodes*3 (a square matrix).
   !! This is the block matrix that maps each field of the distributed source mesh to its identical field on the point (lumped) mesh. 
   !! (i.e., force component of force; moment component of moment).
   
   if (.not. allocated(dM%li) ) then
      call AllocAry(dM%li, Point_Dest%Nnodes*3, Line2_Src%Nnodes*3, 'dM%li', ErrStat2, ErrMsg2 ) !Line2_Src%Nnodes should equal Point_Dest%Nnodes
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
   elseif (size(dM%li,1) /= Point_Dest%Nnodes*3 .or. size(dM%li,2) /= Line2_Src%Nnodes*3) then
      deallocate(dM%li)
      call AllocAry(dM%li, Point_Dest%Nnodes*3, Line2_Src%Nnodes*3, 'dM%li', ErrStat2, ErrMsg2 ) !Line2_Src%Nnodes should equal Point_Dest%Nnodes
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
   end if
      
   dM%li = 0.0_R8Ki                           
      ! loop over source mesh, integrating over each element
   do i = 1, Line2_Src%ElemTable(ELEMENT_LINE2)%nelem

      n1 = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
      n2 = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)
         
      c    = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%det_jac / 3.0_R8Ki  != TwoNorm( p(n2)-p(n1) )/6 = det_jac/3
      TwoC = 2.0_R8Ki * c

      do iComp=1,3 !3 components of force or moment
         j = (n1-1)*3 + iComp      ! this is the index of the node/component of the lumped value
         k = (n2-1)*3 + iComp      ! the is the portion of the "lesser" node
         
         dM%li(j,j) =dM%li(j,j) + TwoC
         dM%li(k,k) =dM%li(k,k) + TwoC
         dM%li(j,k) =dM%li(j,k) + c
         dM%li(k,j) =dM%li(k,j) + c
      enddo !icomp
   enddo !i = iElem         
      
                              
         
      ! M_uD and m_f for moments:
         
   if ( Point_Dest%FieldMask(MASKID_MOMENT) .AND. Line2_Src%FieldMask(MASKID_FORCE) ) then
            
      call FormMatrix_Lump_Line2_to_Point( Line2_Src, dM, ErrStat2, ErrMsg2, Line2_Src ) !Line2_Src is assumed to contain the displaced positions as well as the loads
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
      
      !> Matrix \f$ M_{um}^{SL} \f$, stored in modmesh_mapping::meshmaplinearizationtype::m_ud,
      !! is allocated to be size Point_Dest\%NNodes*3, Line2_Src\%NNodes*3 (a square matrix).
      !! This is the block matrix that relates the translational displacement (u) to the lumped moment.
      !! > Note this could be considered either m_ud or m_us, depending on context.         
                  
      !> Matrix \f$ M_{fm}^{SL} \f$, stored in modmesh_mapping::meshmaplinearizationtype::m_f,
      !! is allocated to be size Point_Dest\%NNodes*3, Line2_Src\%NNodes*3 (a square matrix).
      !> This is the block matrix that relates the force (f) to the lumped moment.
         

      !> We also transfer the force to the destination mesh so that we can use these values for creating the appropriate matrices later
      if (Point_Dest%FieldMask(MASKID_FORCE) ) then
         Point_Dest%Force  = 0.
         ! loop over source mesh, integrating over each element
         do i = 1, Line2_Src%ElemTable(ELEMENT_LINE2)%nelem

            n1 = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
            n2 = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)
        
            c = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%det_jac / 3.0_ReKi  != TwoNorm( p(n2)-p(n1) )/6 = det_jac_Ovr3
                           
            Point_Dest%Force(:,n1) = Point_Dest%Force(:,n1) + c * ( 2.0_ReKi*Line2_Src%Force(:,n1) +    Line2_Src%Force(:,n2) )
            Point_Dest%Force(:,n2) = Point_Dest%Force(:,n2) + c * (    Line2_Src%Force(:,n1) + 2.0_ReKi*Line2_Src%Force(:,n2) )

         end do
      end if
   end if
   
            
END SUBROUTINE Linearize_Lump_Line2_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FormMatrix_Lump_Line2_to_Point( Mesh, dM, ErrStat, ErrMsg, DispMesh )

   TYPE(MeshType),                  INTENT(IN   ) :: Mesh         !< mesh that has loads being lumped
   TYPE(MeshType),                  INTENT(IN   ) :: DispMesh     !< mesh that has displaced position of Mesh
   type(MeshMapLinearizationType),  intent(inout) :: dM           !< data structure for linearization

   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat      !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg       !< Error message if ErrStat /= ErrID_None

   
   ! local variables
   REAL(ReKi)                            :: c, x(3), SSmat(3,3)

   INTEGER(IntKi)                        :: n1_start, n1_end  ! node numbers
   INTEGER(IntKi)                        :: n2_start, n2_end  ! node numbers
   
   INTEGER(IntKi)                        :: n1, n2  ! node numbers
   INTEGER(IntKi)                        :: i  
   Integer(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   character(*),   parameter             :: RoutineName = 'FormMatrix_Lump_Line2_to_Point'


   ErrStat = ErrID_None
   ErrMsg = ""

      ! M_uD and m_f for moments:
                           
   !> Matrix \f$ M_{uDm}^{L} \f$, stored in modmesh_mapping::meshmaplinearizationtype::m_ud,
   !! is allocated to be size Mesh\%NNodes*3, Mesh\%NNodes*3 (a square matrix).
   !> This is the block matrix that relates the destination translational displacement (u) to the lumped moment.
      
      if (.not. allocated(dM%m_uD) ) then
         call AllocAry(dM%m_uD, Mesh%Nnodes*3, Mesh%Nnodes*3, 'dM%m_uD', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      elseif (size(dM%m_uD,1) /= Mesh%Nnodes*3 .or. size(dM%m_uD,2) /= Mesh%Nnodes*3) then
         deallocate(dM%m_uD)
         call AllocAry(dM%m_uD, Mesh%Nnodes*3, Mesh%Nnodes*3, 'dM%m_uD', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      end if
            
   !> Matrix \f$ M_{fm}^{L} \f$, stored in modmesh_mapping::meshmaplinearizationtype::m_f,
   !! is allocated to be size Mesh\%NNodes*3, Mesh\%NNodes*3 (a square matrix).
   !> This is the block matrix that relates the force (f) to the lumped moment.
      
      if (.not. allocated(dM%m_f) ) then
         call AllocAry(dM%m_f, Mesh%Nnodes*3, Mesh%Nnodes*3, 'dM%m_f', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      elseif (size(dM%m_f,1) /= Mesh%Nnodes*3 .or. size(dM%m_f,2) /= Mesh%Nnodes*3) then
         deallocate(dM%m_f)
         call AllocAry(dM%m_f, Mesh%Nnodes*3, Mesh%Nnodes*3, 'dM%m_f', ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
      end if
            
      
      dM%m_uD = 0.0_ReKi                 
      dM%m_f = 0.0_ReKi   
      do i = 1, Mesh%ElemTable(ELEMENT_LINE2)%nelem

         n1 = Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
         n2 = Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)

         n1_start = (n1-1)*3+1
         n1_end   = n1_start+2
         n2_start = (n2 - 1)*3+1
         n2_end   = n2_start+2
                  
         c    = Mesh%ElemTable(ELEMENT_LINE2)%Elements(i)%det_jac / 6.0_ReKi  != TwoNorm( p(n2)-p(n1) )/12 

            ! dM%m_uS (portion that gets multiplied by u [TranslationDisp])
         
         x     = Mesh%Force(:,n1) + Mesh%Force(:,n2)
         x     = x*c         
         SSmat = SkewSymMat(x)         

         dM%m_uD( n1_start:n1_end, n1_start:n1_end ) = dM%m_uD( n1_start:n1_end, n1_start:n1_end ) + SSmat
         dM%m_uD( n2_start:n2_end, n1_start:n1_end ) = dM%m_uD( n2_start:n2_end, n1_start:n1_end ) - SSmat                  
         
         dM%m_uD( n1_start:n1_end, n2_start:n2_end ) = dM%m_uD( n1_start:n1_end, n2_start:n2_end ) - SSmat
         dM%m_uD( n2_start:n2_end, n2_start:n2_end ) = dM%m_uD( n2_start:n2_end, n2_start:n2_end ) + SSmat
         
         
            ! m_f (portion that gets multiplied by f [Force])
         x     =  Mesh%Position(:,n2) + DispMesh%TranslationDisp(:,n2) - &
                 (Mesh%Position(:,n1) + DispMesh%TranslationDisp(:,n1))
         x     = x*c         
         SSmat = SkewSymMat(x)
         
         dM%m_f( n1_start:n1_end, n1_start:n1_end ) = dM%m_f( n1_start:n1_end, n1_start:n1_end ) + SSmat
         dM%m_f( n2_start:n2_end, n1_start:n1_end ) = dM%m_f( n2_start:n2_end, n1_start:n1_end ) - SSmat                  
                                   
         dM%m_f( n1_start:n1_end, n2_start:n2_end ) = dM%m_f( n1_start:n1_end, n2_start:n2_end ) + SSmat
         dM%m_f( n2_start:n2_end, n2_start:n2_end ) = dM%m_f( n2_start:n2_end, n2_start:n2_end ) - SSmat
      enddo !i = iElem    
         
   
END SUBROUTINE FormMatrix_Lump_Line2_to_Point 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine creates the linearization matrices for transfer of fields between point meshes.
!! Mapping equations are as follows: 
!!
!> Rotational Displacement: \f$ \frac{\partial M_\Lambda}{\partial x} = \begin{bmatrix} M_{mi} \end{bmatrix} \f$ for source fields \f$ \left\{ \vec{\theta}^S\right\}\f$
!!
!> Translational Displacement: \f$ \frac{\partial M_u}{\partial x} = \begin{bmatrix} M_{mi} & M_{f_{\times p}} \end{bmatrix} \f$ for source fields
!! \f$ \left\{ \begin{matrix}
!!      \vec{u}^S \\
!!      \vec{\theta}^S
!! \end{matrix} \right\} \f$
!!
!> Rotational Velocity: \f$ \frac{\partial M_\omega}{\partial x} = \begin{bmatrix} M_{mi} \end{bmatrix} \f$ for source fields \f$\left\{\vec{\omega}^S\right\}\f$
!!
!! Translational Velocity: \f$ \frac{\partial M_v}{\partial x} = \begin{bmatrix}  M_{tv\_uD} & M_{tv\_uS} & M_{mi} & M_{f_{\times p}} \end{bmatrix} \f$ for source fields
!! \f$ \left\{ \begin{matrix}
!!      \vec{u}^D \\
!!      \vec{u}^S \\
!!      \vec{v}^S \\
!!      \vec{\omega}^S
!! \end{matrix} \right\} \f$
!!
!! Rotational Acceleration: \f$ \frac{\partial M_\alpha}{\partial x} = \begin{bmatrix} M_{mi} \end{bmatrix} \f$ for source fields \f$\left\{\vec{\alpha}^S\right\}\f$
!!
!! Translational Acceleration: \f$ \frac{\partial M_a}{\partial x} = \begin{bmatrix} M_{ta\_uD} & M_{ta\_uS} & M_{ta\_rv} & M_{mi} & M_{f_{\times p}} \end{bmatrix} \f$ for source fields
!! \f$ \left\{ \begin{matrix}
!!      \vec{u}^D \\
!!      \vec{u}^S \\
!!      \vec{\omega}^S \\
!!      \vec{a}^S \\
!!      \vec{\alpha}^S
!! \end{matrix} \right\} \f$
!!
!! Scalar quantities: \f$ \frac{\partial M_S}{\partial x} = \begin{bmatrix} M_{mi} \end{bmatrix} \f$ for source fields \f$\left\{\vec{S}^S\right\}\f$
!!
!! Forces: \f$ \frac{\partial M_F}{\partial x} = \begin{bmatrix} M_{li} \end{bmatrix} \f$ for source fields \f$\left\{\vec{F}^S\right\}\f$
!!
!! Moments: \f$ \frac{\partial M_M}{\partial x} = \begin{bmatrix} M_{uDm} & M_{uSm} & M_{fm} & M_{li} \end{bmatrix} \f$ for source fields
!! \f$ \left\{ \begin{matrix}
!!      \vec{u}^D \\
!!      \vec{u}^S \\
!!      \vec{F}^S \\
!!      \vec{M}^S
!! \end{matrix} \right\} \f$
!!
!! \f$M_{mi}\f$ is modmesh_mapping::meshmaplinearizationtype::mi \n
!! \f$M_{f_{\times p}}\f$ is modmesh_mapping::meshmaplinearizationtype::fx_p \n
!! \f$M_{tv\_uD}\f$ is modmesh_mapping::meshmaplinearizationtype::tv_uD \n
!! \f$M_{tv\_uS}\f$ is modmesh_mapping::meshmaplinearizationtype::tv_uS \n
!! \f$M_{ta\_uD}\f$ is modmesh_mapping::meshmaplinearizationtype::ta_uD \n
!! \f$M_{ta\_uS}\f$ is modmesh_mapping::meshmaplinearizationtype::ta_uS \n
!! \f$M_{ta\_rv}\f$ is modmesh_mapping::meshmaplinearizationtype::ta_rv \n
!! \f$M_{li}\f$ is modmesh_mapping::meshmaplinearizationtype::li \n
!! \f$M_{uSm}\f$ is modmesh_mapping::meshmaplinearizationtype::m_us \n
!! \f$M_{uDm}\f$ is modmesh_mapping::meshmaplinearizationtype::m_ud \n
!! \f$M_{fm}\f$ is modmesh_mapping::meshmaplinearizationtype::m_f 
!!
SUBROUTINE Linearize_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),          INTENT(IN   ) :: Src       !< source point mesh
   TYPE(MeshType),          INTENT(IN   ) :: Dest      !< destination mesh
   TYPE(MeshType), OPTIONAL,INTENT(IN   ) :: SrcDisp   !< an optional mesh, which contains the displacements associated with the source if the source contains load information
   TYPE(MeshType), OPTIONAL,INTENT(IN   ) :: DestDisp  !< an optional mesh, which contains the displacements associated with the destination if the destination contains load information

   TYPE(MeshMapType),       INTENT(INOUT) :: MeshMap   !< Mapping(s) between Src and Dest meshes; also contains jacobian of this mapping

   INTEGER(IntKi),          INTENT(  OUT) :: ErrStat   !< Error status of the operation
   CHARACTER(*),            INTENT(  OUT) :: ErrMsg    !< Error message if ErrStat /= ErrID_None


   ! local variables
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2
   CHARACTER(*), PARAMETER                 :: RoutineName = 'Linearize_Point_to_Point'
   
   
   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! ------------------------------------------------------------------------------------------------------------------------------
   !> ### Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   ! If Src is displacements/velocities then loop over the Destination Mesh;
   !    each motion/scalar in the destination mesh needs to be interpolated from the source mesh.


   if ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) then

      !........................
      !> * Create mapping (if remapFlag)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMotionMap_P_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN
            
      end if

      !........................
      !> * Create linearization matrices for motions
      !........................
      
      call Linearize_Motions_Point_to_Point(Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN
                                                                               
   end if ! HasMotionFields: algorithm for motions/scalars

   ! ------------------------------------------------------------------------------------------------------------------------------
   !> ### Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------

   if ( HasLoadFields(Src) ) then
            
      !........................
      !> * Create mapping (if RemapFlag)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then
         
         CALL CreateLoadMap_P_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN                        

      end if

      !........................
      !> * Create linearization matrices for loads
      !........................      

      IF ( PRESENT( SrcDisp ) .AND. PRESENT( DestDisp ) ) THEN
                                       
         call Linearize_Loads_Point_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2, SrcDisp, DestDisp )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN                        
             
      ELSE !bjj: we could check if Src also had TranslationDisp fields and call with SrcDisp=Src
         CALL SetErrStat( ErrID_Fatal, 'Meshes with load fields must also have sibling meshes with displacment field.', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF

   end if ! algorithm for loads



END SUBROUTINE Linearize_Point_to_Point
!----------------------------------------------------------------------------------------------------------------------------------


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine is used for debugging mesh mappings. It checks to see that the sum of the loads on each mesh is the same as when
!! all the loads from each mesh are lumped to a single point.
SUBROUTINE WriteMappingTransferToFile(Mesh1_I,Mesh1_O,Mesh2_I,Mesh2_O,Map_Mod1_Mod2,Map_Mod2_Mod1,BinOutputName)

   TYPE(meshtype),    intent(in) :: mesh1_I              !< mesh 1 inputs
   TYPE(meshtype),    intent(in) :: mesh1_O              !< mesh 1 outputs  
   TYPE(meshtype),    intent(in) :: mesh2_I              !< mesh 2 inputs
   TYPE(meshtype),    intent(in) :: mesh2_O              !< mesh 2 outputs
   
   TYPE(MeshMapType), intent(in) :: Map_Mod1_Mod2        !< Data for mapping meshes from mesh 1 (outputs) to mesh 2 (inputs)
   TYPE(MeshMapType), intent(in) :: Map_Mod2_Mod1        !< Data for mapping meshes from mesh 2 (outputs) to mesh 1 (inputs)

   CHARACTER(*),      INTENT(IN) :: BinOutputName        !< name of binary output file
   
   
   ! local variables:
   TYPE(meshtype)                         :: mesh_Motion_1PT, mesh1_I_1PT, mesh2_O_1PT
   TYPE(MeshMapType)                      :: Map_Mod2_O_1PT, Map_Mod1_I_1PT
      
   INTEGER(IntKi)                         :: i
   INTEGER(IntKi)                         :: un_out
   INTEGER(IntKi)                         :: ErrStat          ! Error status of the operation
   CHARACTER(ErrMsgLen)                   :: ErrMsg           ! Error message if ErrStat /= ErrID_None
   CHARACTER(256)                         :: PrintWarnF, PrintWarnM, TmpValues

#ifdef MESH_DEBUG     
   TYPE(meshtype)                         :: mesh2_O_1PT_augmented, mesh2_O_1PT_lumped, mesh1_I_1PT_lumped
   TYPE(MeshMapType)                      :: Map_Mod2_O_1PT_augmented, Map_Mod2_O_1PT_lumped, Map_Mod1_I_1PT_lumped
#endif      
   
   
   !------------------------------------------------------------------------
   ! Make sure the meshes are committed before checking them:
   !------------------------------------------------------------------------
   
   IF (.NOT. mesh1_I%Committed .OR. .NOT. mesh1_O%Committed ) RETURN
   IF (.NOT. mesh2_I%Committed .OR. .NOT. mesh2_O%Committed ) RETURN
      
   !------------------------------------------------------------------------
   ! lump the loads to one point and compare:
   !------------------------------------------------------------------------       
   
   ! create one loads mesh with one point:
   CALL MeshCreate( BlankMesh       = mesh1_I_1PT        &
                  , IOS              = COMPONENT_INPUT   &
                  , NNodes           = 1                 &
                  , Force            = .TRUE.            &
                  , Moment           = .TRUE.            &
                  , ErrStat          = ErrStat           &
                  , ErrMess          = ErrMsg            )
      
   IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
         
         
   CALL MeshPositionNode ( mesh1_I_1PT, 1, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/), ErrStat, ErrMsg ) ; IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             
   CALL MeshConstructElement ( mesh1_I_1PT, ELEMENT_POINT, ErrStat, ErrMsg, 1 );                 IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                                       
   CALL MeshCommit ( mesh1_I_1PT, ErrStat, ErrMsg )                                       ;      IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      
   !.....         
   ! create a corresponding motion mesh with one point:
   
   CALL MeshCopy( mesh1_I_1PT, mesh_Motion_1PT, MESH_SIBLING, ErrStat, ErrMsg &
                  , IOS              = COMPONENT_OUTPUT  &
                  , TranslationDisp  = .TRUE.            ) ;                                     IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                     
   !.....         
   ! create a second loads mesh with one point:
   CALL MeshCopy( mesh1_I_1PT, mesh2_O_1PT, MESH_NEWCOPY, ErrStat, ErrMsg )  ! This thinks it's for input, but really it's for output. I don't think it matters...       
       
   !.....         
   ! create the mapping data structures:       
   CALL MeshMapCreate( Mesh1_I, Mesh1_I_1PT, Map_Mod1_I_1PT,  ErrStat, ErrMsg );                 IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshMapCreate( Mesh2_O, Mesh2_O_1PT, Map_Mod2_O_1PT,  ErrStat, ErrMsg );                 IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

   !.....         
   ! transfer MESH1_I (loads) to single point:
   
   IF ( mesh1_I%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN
      CALL Transfer_Point_to_Point( Mesh1_I, Mesh1_I_1PT, Map_Mod1_I_1PT,ErrStat,ErrMsg,mesh1_O,mesh_Motion_1PT);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))         
   ELSEIF ( mesh1_I%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN
      CALL Transfer_Line2_to_Point( Mesh1_I, Mesh1_I_1PT, Map_Mod1_I_1PT,ErrStat,ErrMsg,mesh1_O,mesh_Motion_1PT);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))          
   END IF
   
   !.....         
   ! transfer Mesh2_O (loads) to single point:      
   IF ( Mesh2_O%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN
      CALL Transfer_Point_to_Point( Mesh2_O, Mesh2_O_1PT, Map_Mod2_O_1PT,ErrStat,ErrMsg,mesh2_I,mesh_Motion_1PT);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
   ELSEIF ( Mesh2_O%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN 
      CALL Transfer_Line2_to_Point( Mesh2_O, Mesh2_O_1PT, Map_Mod2_O_1PT,ErrStat,ErrMsg,mesh2_I,mesh_Motion_1PT);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))          
   END IF
   !............
            
   ! dsisplay a warning if the point loads are not equal:         
   PrintWarnF=""
   PrintWarnM=""
   do i=1,3
      if (.NOT. equalrealnos(mesh1_I_1PT%Force( i,1),mesh2_O_1PT%Force( i,1)) ) PrintWarnF=NewLine//"  <----------- WARNING: Forces are not equal ----------->  "//NewLine//NewLine
      if (.NOT. equalrealnos(mesh1_I_1PT%Moment(i,1),mesh2_O_1PT%Moment(i,1)) ) PrintWarnM=NewLine//"  <----------- WARNING: Moments are not equal ----------->  "//NewLine//NewLine
   end do

   
   call wrscr(TRIM(PrintWarnF)//'Total Force:' )
   write(TmpValues,*) mesh1_I_1PT%Force;   call wrscr('     Mesh 1: '//TRIM(TmpValues))
   write(TmpValues,*) mesh2_O_1PT%Force;   call wrscr('     Mesh 2: '//TRIM(TmpValues))
   call wrscr(TRIM(PrintWarnM)//'Total Moment:' )
   write(TmpValues,*) mesh1_I_1PT%Moment;  call wrscr('     Mesh 1: '//TRIM(TmpValues))
   write(TmpValues,*) mesh2_O_1PT%Moment;  call wrscr('     Mesh 2: '//TRIM(TmpValues))
   !............
   
   !------------------------------------------------------------------------
   ! now we'll write all the mesh info to a file for debugging:   
   !------------------------------------------------------------------------

   un_out = -1
   CALL MeshWrBin ( un_out, Mesh1_I,         ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshWrBin ( un_out, Mesh1_O,         ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshWrBin ( un_out, Mesh2_I,         ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshWrBin ( un_out, Mesh2_O,         ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      
   CALL MeshWrBin ( un_out, mesh1_I_1PT,     ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshWrBin ( un_out, mesh2_O_1PT,     ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))         
   CALL MeshWrBin ( un_out, mesh_Motion_1PT, ErrStat, ErrMsg, BinOutputName);  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))   
      
   CALL MeshMapWrBin( un_out, Mesh1_O, Mesh2_I, Map_Mod1_Mod2, ErrStat, ErrMsg, BinOutputName );  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
   CALL MeshMapWrBin( un_out, Mesh2_O, Mesh1_I, Map_Mod2_Mod1, ErrStat, ErrMsg, BinOutputName );  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
      
   CALL MeshMapWrBin( un_out, Mesh1_I, Mesh1_I_1PT, Map_Mod1_I_1PT, ErrStat, ErrMsg, BinOutputName );  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 
   CALL MeshMapWrBin( un_out, Mesh2_O, Mesh2_O_1PT, Map_Mod2_O_1PT, ErrStat, ErrMsg, BinOutputName );  IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg)) 

   close( un_out )
   
   
   
#ifdef MESH_DEBUG  
      
if (LEN_TRIM(PrintWarnM)>0 .OR. LEN_TRIM(PrintWarnF)>0 ) THEN
      call wrscr(NewLine)
      call wrmatrix( mesh1_O%TranslationDisp, CU, 'f10.5')
      call wrscr(NewLine)
     ! call wrmatrix( mesh2_I%TranslationDisp, CU, 'f10.5')
     ! call wrmatrix( mesh_Motion_1PT%TranslationDisp, CU, 'f10.5')
     ! call wrmatrix( mesh_Motion_1PT%Position, CU, 'f10.5')


   if (Map_Mod2_Mod1%Augmented_Ln2_Src%committed) THEN
      CALL MeshCopy( Mesh2_O_1PT, mesh2_O_1PT_augmented, MESH_NEWCOPY, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      mesh2_O_1PT_augmented%force = 0.0
      mesh2_O_1PT_augmented%moment = 0.0
     ! Map_Mod2_Mod1%Augmented_Ln2_Src%TranslationDisp=0.0
      CALL MeshMapCreate( Map_Mod2_Mod1%Augmented_Ln2_Src,  mesh2_O_1PT_augmented, Map_Mod2_O_1PT_augmented, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL Transfer_Line2_to_Point(Map_Mod2_Mod1%Augmented_Ln2_Src,  mesh2_O_1PT_augmented, Map_Mod2_O_1PT_augmented, ErrStat, ErrMsg, Map_Mod2_Mod1%Augmented_Ln2_Src, mesh_Motion_1PT );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

      !print *, '---------Augmented_Ln2_Src:---------------'
      !call meshprintinfo(CU,Map_Mod2_Mod1%Augmented_Ln2_Src)
   end if
      
   if (Map_Mod2_Mod1%Lumped_Points_Src%committed) THEN
      CALL MeshCopy( Mesh2_O_1PT, mesh2_O_1PT_lumped,    MESH_NEWCOPY, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshMapCreate( Map_Mod2_Mod1%Lumped_Points_Src,  mesh2_O_1PT_lumped,    Map_Mod2_O_1PT_lumped,    ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL Transfer_Point_to_Point(Map_Mod2_Mod1%Lumped_Points_Src,  mesh2_O_1PT_lumped,    Map_Mod2_O_1PT_lumped,    ErrStat, ErrMsg, Map_Mod2_Mod1%Lumped_Points_Src, mesh_Motion_1PT );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))            

      !print *, '-------------Lumped_Points_Src:----------------'
      !call meshprintinfo(CU,Map_Mod2_Mod1%Lumped_Points_Src)   
   end if

   if (Map_Mod2_Mod1%Lumped_Points_Dest%committed) THEN
      CALL MeshCopy( Mesh1_I_1PT, mesh1_I_1PT_lumped,    MESH_NEWCOPY, ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
      CALL MeshMapCreate( Map_Mod2_Mod1%Lumped_Points_Dest, mesh1_I_1PT_lumped,    Map_Mod1_I_1PT_lumped,    ErrStat, ErrMsg );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))      
      CALL Transfer_Point_to_Point(Map_Mod2_Mod1%Lumped_Points_Dest, mesh1_I_1PT_lumped,    Map_Mod1_I_1PT_lumped,    ErrStat, ErrMsg, mesh1_O, mesh_Motion_1PT );       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))     

      !print *, '-----------Lumped_Points_Dest:---------------'
      !call meshprintinfo(CU,Map_Mod2_Mod1%Lumped_Points_Dest)   
   end if
         
            
      call wrscr('Total Force:' )     
      if (mesh2_O_1PT_augmented%committed) then
         write(TmpValues,*) mesh2_O_1PT_augmented%Force
         call wrscr('     Mesh 2 (augmented):'//TRIM(TmpValues)//trim(num2lstr(Map_Mod2_Mod1%Augmented_Ln2_Src%Nnodes  )) )
      end if
      
      if (mesh2_O_1PT_lumped%committed) then
         write(TmpValues,*) mesh2_O_1PT_lumped%Force
         call wrscr('     Mesh 2 (lumped):   '//TRIM(TmpValues)//trim(num2lstr(Map_Mod2_Mod1%Lumped_Points_Src%Nnodes  )) )
      end if
      
      if (mesh1_I_1PT_lumped%committed) then
         write(TmpValues,*) mesh1_I_1PT_lumped%Force
         call wrscr('     Mesh 1 (lumped):   '//TRIM(TmpValues)//trim(num2lstr(Map_Mod2_Mod1%Lumped_Points_Dest%Nnodes )) )
      end if
      
      call wrscr('Total Moment:' )
      if (mesh2_O_1PT_augmented%committed) then
         write(TmpValues,*) mesh2_O_1PT_augmented%Moment
         call wrscr('     Mesh 2 (augmented):'//TRIM(TmpValues)//trim(num2lstr(Map_Mod2_Mod1%Augmented_Ln2_Src%Nnodes  )) )
      end if
      
      if (mesh2_O_1PT_lumped%committed) then
         write(TmpValues,*) mesh2_O_1PT_lumped%Moment
         call wrscr('     Mesh 2 (lumped):   '//TRIM(TmpValues)//trim(num2lstr(Map_Mod2_Mod1%Lumped_Points_Src%Nnodes  )) )
      end if
      
      if (mesh1_I_1PT_lumped%committed) then
         write(TmpValues,*) mesh1_I_1PT_lumped%Moment
         call wrscr('     Mesh 1 (lumped):   '//TRIM(TmpValues)//trim(num2lstr(Map_Mod2_Mod1%Lumped_Points_Dest%Nnodes )) )
      end if
      
END IF
#endif
   
   
   !------------------------------------------------------------------------
   ! destroy local copies:
   !------------------------------------------------------------------------
   
   CALL MeshDestroy( mesh_Motion_1PT, ErrStat, ErrMsg ); IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshDestroy( mesh1_I_1PT, ErrStat, ErrMsg );     IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshDestroy( mesh2_O_1PT, ErrStat, ErrMsg );     IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
   call MeshMapDestroy(Map_Mod1_I_1PT, ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   call MeshMapDestroy(Map_Mod2_O_1PT, ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   
   
#ifdef MESH_DEBUG  
   CALL MeshDestroy( mesh2_O_1PT_augmented,      ErrStat, ErrMsg );   IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshDestroy( mesh2_O_1PT_lumped,         ErrStat, ErrMsg );   IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   CALL MeshDestroy( mesh1_I_1PT_lumped,         ErrStat, ErrMsg );   IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))

   call MeshMapDestroy(Map_Mod2_O_1PT_augmented, ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   call MeshMapDestroy(Map_Mod2_O_1PT_lumped,    ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
   call MeshMapDestroy(Map_Mod1_I_1PT_lumped,    ErrStat, ErrMsg);IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))
#endif      
   

END SUBROUTINE WriteMappingTransferToFile 
!----------------------------------------------------------------------------------------------------------------------------------

!==================================================================================================================================
!bjj: these routines require the use of ModMesh.f90, thus they cannot be part of NWTC_Library_Types.f90:
!STARTOFREGISTRYGENERATEDFILE 'NWTC_Library_Types.f90'
!
! WARNING This file is generated automatically by the FAST registry
! Do not edit.  Your changes to this file will be lost.
!
! FAST Registry (v3.02.00, 23-Jul-2016)
!*********************************************************************************************************************************
 SUBROUTINE NWTC_Library_CopyMapType( SrcMapTypeData, DstMapTypeData, CtrlCode, ErrStat, ErrMsg )
   TYPE(MapType), INTENT(IN) :: SrcMapTypeData
   TYPE(MapType), INTENT(INOUT) :: DstMapTypeData
   INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
! Local 
   INTEGER(IntKi)                 :: i,j,k
   INTEGER(IntKi)                 :: i1, i1_l, i1_u  !  bounds (upper/lower) for an array dimension 1
   INTEGER(IntKi)                 :: i2, i2_l, i2_u  !  bounds (upper/lower) for an array dimension 2
   INTEGER(IntKi)                 :: i3, i3_l, i3_u  !  bounds (upper/lower) for an array dimension 3
   INTEGER(IntKi)                 :: ErrStat2
   CHARACTER(ErrMsgLen)           :: ErrMsg2
   CHARACTER(*), PARAMETER        :: RoutineName = 'NWTC_Library_CopyMapType'
! 
   ErrStat = ErrID_None
   ErrMsg  = ""
    DstMapTypeData%OtherMesh_Element = SrcMapTypeData%OtherMesh_Element
    DstMapTypeData%distance = SrcMapTypeData%distance
    DstMapTypeData%couple_arm = SrcMapTypeData%couple_arm
    DstMapTypeData%shape_fn = SrcMapTypeData%shape_fn
 END SUBROUTINE NWTC_Library_CopyMapType

 SUBROUTINE NWTC_Library_DestroyMapType( MapTypeData, ErrStat, ErrMsg )
  TYPE(MapType), INTENT(INOUT) :: MapTypeData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
  CHARACTER(*),    PARAMETER :: RoutineName = 'NWTC_Library_DestroyMapType'
  INTEGER(IntKi)                 :: i, i1, i2, i3, i4, i5 
! 
  ErrStat = ErrID_None
  ErrMsg  = ""
 END SUBROUTINE NWTC_Library_DestroyMapType

 SUBROUTINE NWTC_Library_PackMapType( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )
  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)
  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)
  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)
  TYPE(MapType),  INTENT(IN) :: InData
  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat
  CHARACTER(*),     INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly
    ! Local variables
  INTEGER(IntKi)                 :: Re_BufSz
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_BufSz
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_BufSz
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i,i1,i2,i3,i4,i5
  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'NWTC_Library_PackMapType'
 ! buffers to store subtypes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)

  OnlySize = .FALSE.
  IF ( PRESENT(SizeOnly) ) THEN
    OnlySize = SizeOnly
  ENDIF
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_BufSz  = 0
  Db_BufSz  = 0
  Int_BufSz  = 0
      Int_BufSz  = Int_BufSz  + 1  ! OtherMesh_Element
      Db_BufSz   = Db_BufSz   + 1  ! distance
      Db_BufSz   = Db_BufSz   + SIZE(InData%couple_arm)  ! couple_arm
      Db_BufSz   = Db_BufSz   + SIZE(InData%shape_fn)  ! shape_fn
  IF ( Re_BufSz  .GT. 0 ) THEN 
     ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Db_BufSz  .GT. 0 ) THEN 
     ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Int_BufSz  .GT. 0 ) THEN 
     ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)

  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred = 1

    IntKiBuf(Int_Xferred) = InData%OtherMesh_Element
    Int_Xferred = Int_Xferred + 1
    DbKiBuf(Db_Xferred) = InData%distance
    Db_Xferred = Db_Xferred + 1
    DO i1 = LBOUND(InData%couple_arm,1), UBOUND(InData%couple_arm,1)
      DbKiBuf(Db_Xferred) = InData%couple_arm(i1)
      Db_Xferred = Db_Xferred + 1
    END DO
    DO i1 = LBOUND(InData%shape_fn,1), UBOUND(InData%shape_fn,1)
      DbKiBuf(Db_Xferred) = InData%shape_fn(i1)
      Db_Xferred = Db_Xferred + 1
    END DO
 END SUBROUTINE NWTC_Library_PackMapType

 SUBROUTINE NWTC_Library_UnPackMapType( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )
  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)
  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)
  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)
  TYPE(MapType), INTENT(INOUT) :: OutData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
    ! Local variables
  INTEGER(IntKi)                 :: Buf_size
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i
  INTEGER(IntKi)                 :: i1, i1_l, i1_u  !  bounds (upper/lower) for an array dimension 1
  INTEGER(IntKi)                 :: i2, i2_l, i2_u  !  bounds (upper/lower) for an array dimension 2
  INTEGER(IntKi)                 :: i3, i3_l, i3_u  !  bounds (upper/lower) for an array dimension 3
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'NWTC_Library_UnPackMapType'
 ! buffers to store meshes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred  = 1
    OutData%OtherMesh_Element = IntKiBuf(Int_Xferred)
    Int_Xferred = Int_Xferred + 1
    OutData%distance = REAL(DbKiBuf(Db_Xferred), R8Ki)
    Db_Xferred = Db_Xferred + 1
    i1_l = LBOUND(OutData%couple_arm,1)
    i1_u = UBOUND(OutData%couple_arm,1)
    DO i1 = LBOUND(OutData%couple_arm,1), UBOUND(OutData%couple_arm,1)
      OutData%couple_arm(i1) = REAL(DbKiBuf(Db_Xferred), R8Ki)
      Db_Xferred = Db_Xferred + 1
    END DO
    i1_l = LBOUND(OutData%shape_fn,1)
    i1_u = UBOUND(OutData%shape_fn,1)
    DO i1 = LBOUND(OutData%shape_fn,1), UBOUND(OutData%shape_fn,1)
      OutData%shape_fn(i1) = REAL(DbKiBuf(Db_Xferred), R8Ki)
      Db_Xferred = Db_Xferred + 1
    END DO
 END SUBROUTINE NWTC_Library_UnPackMapType

 SUBROUTINE NWTC_Library_CopyMeshMapLinearizationType( SrcMeshMapLinearizationTypeData, DstMeshMapLinearizationTypeData, CtrlCode, ErrStat, ErrMsg )
   TYPE(MeshMapLinearizationType), INTENT(IN) :: SrcMeshMapLinearizationTypeData
   TYPE(MeshMapLinearizationType), INTENT(INOUT) :: DstMeshMapLinearizationTypeData
   INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
! Local 
   INTEGER(IntKi)                 :: i,j,k
   INTEGER(IntKi)                 :: i1, i1_l, i1_u  !  bounds (upper/lower) for an array dimension 1
   INTEGER(IntKi)                 :: i2, i2_l, i2_u  !  bounds (upper/lower) for an array dimension 2
   INTEGER(IntKi)                 :: ErrStat2
   CHARACTER(ErrMsgLen)           :: ErrMsg2
   CHARACTER(*), PARAMETER        :: RoutineName = 'NWTC_Library_CopyMeshMapLinearizationType'
! 
   ErrStat = ErrID_None
   ErrMsg  = ""
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%mi)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%mi,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%mi,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%mi,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%mi,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%mi)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%mi(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%mi.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%mi = SrcMeshMapLinearizationTypeData%mi
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%fx_p)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%fx_p,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%fx_p,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%fx_p,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%fx_p,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%fx_p)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%fx_p(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%fx_p.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%fx_p = SrcMeshMapLinearizationTypeData%fx_p
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%tv_uD)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%tv_uD,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%tv_uD,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%tv_uD,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%tv_uD,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%tv_uD)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%tv_uD(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%tv_uD.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%tv_uD = SrcMeshMapLinearizationTypeData%tv_uD
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%tv_uS)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%tv_uS,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%tv_uS,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%tv_uS,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%tv_uS,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%tv_uS)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%tv_uS(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%tv_uS.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%tv_uS = SrcMeshMapLinearizationTypeData%tv_uS
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%ta_uD)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%ta_uD,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%ta_uD,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%ta_uD,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%ta_uD,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%ta_uD)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%ta_uD(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%ta_uD.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%ta_uD = SrcMeshMapLinearizationTypeData%ta_uD
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%ta_uS)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%ta_uS,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%ta_uS,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%ta_uS,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%ta_uS,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%ta_uS)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%ta_uS(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%ta_uS.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%ta_uS = SrcMeshMapLinearizationTypeData%ta_uS
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%ta_rv)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%ta_rv,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%ta_rv,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%ta_rv,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%ta_rv,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%ta_rv)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%ta_rv(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%ta_rv.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%ta_rv = SrcMeshMapLinearizationTypeData%ta_rv
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%li)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%li,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%li,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%li,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%li,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%li)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%li(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%li.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%li = SrcMeshMapLinearizationTypeData%li
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%M_uS)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%M_uS,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%M_uS,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%M_uS,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%M_uS,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%M_uS)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%M_uS(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%M_uS.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%M_uS = SrcMeshMapLinearizationTypeData%M_uS
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%M_uD)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%M_uD,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%M_uD,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%M_uD,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%M_uD,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%M_uD)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%M_uD(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%M_uD.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%M_uD = SrcMeshMapLinearizationTypeData%M_uD
ENDIF
IF (ALLOCATED(SrcMeshMapLinearizationTypeData%M_f)) THEN
  i1_l = LBOUND(SrcMeshMapLinearizationTypeData%M_f,1)
  i1_u = UBOUND(SrcMeshMapLinearizationTypeData%M_f,1)
  i2_l = LBOUND(SrcMeshMapLinearizationTypeData%M_f,2)
  i2_u = UBOUND(SrcMeshMapLinearizationTypeData%M_f,2)
  IF (.NOT. ALLOCATED(DstMeshMapLinearizationTypeData%M_f)) THEN 
    ALLOCATE(DstMeshMapLinearizationTypeData%M_f(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%M_f.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapLinearizationTypeData%M_f = SrcMeshMapLinearizationTypeData%M_f
ENDIF
 END SUBROUTINE NWTC_Library_CopyMeshMapLinearizationType

 SUBROUTINE NWTC_Library_DestroyMeshMapLinearizationType( MeshMapLinearizationTypeData, ErrStat, ErrMsg )
  TYPE(MeshMapLinearizationType), INTENT(INOUT) :: MeshMapLinearizationTypeData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
  CHARACTER(*),    PARAMETER :: RoutineName = 'NWTC_Library_DestroyMeshMapLinearizationType'
  INTEGER(IntKi)                 :: i, i1, i2, i3, i4, i5 
! 
  ErrStat = ErrID_None
  ErrMsg  = ""
IF (ALLOCATED(MeshMapLinearizationTypeData%mi)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%mi)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%fx_p)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%fx_p)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%tv_uD)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%tv_uD)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%tv_uS)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%tv_uS)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%ta_uD)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%ta_uD)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%ta_uS)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%ta_uS)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%ta_rv)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%ta_rv)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%li)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%li)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%M_uS)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%M_uS)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%M_uD)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%M_uD)
ENDIF
IF (ALLOCATED(MeshMapLinearizationTypeData%M_f)) THEN
  DEALLOCATE(MeshMapLinearizationTypeData%M_f)
ENDIF
 END SUBROUTINE NWTC_Library_DestroyMeshMapLinearizationType

 SUBROUTINE NWTC_Library_PackMeshMapLinearizationType( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )
  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)
  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)
  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)
  TYPE(MeshMapLinearizationType),  INTENT(IN) :: InData
  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat
  CHARACTER(*),     INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly
    ! Local variables
  INTEGER(IntKi)                 :: Re_BufSz
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_BufSz
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_BufSz
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i,i1,i2,i3,i4,i5
  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'NWTC_Library_PackMeshMapLinearizationType'
 ! buffers to store subtypes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)

  OnlySize = .FALSE.
  IF ( PRESENT(SizeOnly) ) THEN
    OnlySize = SizeOnly
  ENDIF
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_BufSz  = 0
  Db_BufSz  = 0
  Int_BufSz  = 0
  Int_BufSz   = Int_BufSz   + 1     ! mi allocated yes/no
  IF ( ALLOCATED(InData%mi) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! mi upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%mi)  ! mi
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! fx_p allocated yes/no
  IF ( ALLOCATED(InData%fx_p) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! fx_p upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%fx_p)  ! fx_p
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! tv_uD allocated yes/no
  IF ( ALLOCATED(InData%tv_uD) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! tv_uD upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%tv_uD)  ! tv_uD
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! tv_uS allocated yes/no
  IF ( ALLOCATED(InData%tv_uS) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! tv_uS upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%tv_uS)  ! tv_uS
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! ta_uD allocated yes/no
  IF ( ALLOCATED(InData%ta_uD) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! ta_uD upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%ta_uD)  ! ta_uD
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! ta_uS allocated yes/no
  IF ( ALLOCATED(InData%ta_uS) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! ta_uS upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%ta_uS)  ! ta_uS
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! ta_rv allocated yes/no
  IF ( ALLOCATED(InData%ta_rv) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! ta_rv upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%ta_rv)  ! ta_rv
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! li allocated yes/no
  IF ( ALLOCATED(InData%li) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! li upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%li)  ! li
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! M_uS allocated yes/no
  IF ( ALLOCATED(InData%M_uS) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! M_uS upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%M_uS)  ! M_uS
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! M_uD allocated yes/no
  IF ( ALLOCATED(InData%M_uD) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! M_uD upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%M_uD)  ! M_uD
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! M_f allocated yes/no
  IF ( ALLOCATED(InData%M_f) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! M_f upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%M_f)  ! M_f
  END IF
  IF ( Re_BufSz  .GT. 0 ) THEN 
     ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Db_BufSz  .GT. 0 ) THEN 
     ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Int_BufSz  .GT. 0 ) THEN 
     ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)

  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred = 1

  IF ( .NOT. ALLOCATED(InData%mi) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%mi,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%mi,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%mi,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%mi,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%mi,2), UBOUND(InData%mi,2)
        DO i1 = LBOUND(InData%mi,1), UBOUND(InData%mi,1)
          DbKiBuf(Db_Xferred) = InData%mi(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%fx_p) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%fx_p,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%fx_p,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%fx_p,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%fx_p,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%fx_p,2), UBOUND(InData%fx_p,2)
        DO i1 = LBOUND(InData%fx_p,1), UBOUND(InData%fx_p,1)
          DbKiBuf(Db_Xferred) = InData%fx_p(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%tv_uD) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%tv_uD,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%tv_uD,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%tv_uD,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%tv_uD,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%tv_uD,2), UBOUND(InData%tv_uD,2)
        DO i1 = LBOUND(InData%tv_uD,1), UBOUND(InData%tv_uD,1)
          DbKiBuf(Db_Xferred) = InData%tv_uD(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%tv_uS) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%tv_uS,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%tv_uS,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%tv_uS,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%tv_uS,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%tv_uS,2), UBOUND(InData%tv_uS,2)
        DO i1 = LBOUND(InData%tv_uS,1), UBOUND(InData%tv_uS,1)
          DbKiBuf(Db_Xferred) = InData%tv_uS(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%ta_uD) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%ta_uD,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%ta_uD,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%ta_uD,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%ta_uD,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%ta_uD,2), UBOUND(InData%ta_uD,2)
        DO i1 = LBOUND(InData%ta_uD,1), UBOUND(InData%ta_uD,1)
          DbKiBuf(Db_Xferred) = InData%ta_uD(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%ta_uS) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%ta_uS,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%ta_uS,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%ta_uS,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%ta_uS,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%ta_uS,2), UBOUND(InData%ta_uS,2)
        DO i1 = LBOUND(InData%ta_uS,1), UBOUND(InData%ta_uS,1)
          DbKiBuf(Db_Xferred) = InData%ta_uS(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%ta_rv) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%ta_rv,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%ta_rv,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%ta_rv,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%ta_rv,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%ta_rv,2), UBOUND(InData%ta_rv,2)
        DO i1 = LBOUND(InData%ta_rv,1), UBOUND(InData%ta_rv,1)
          DbKiBuf(Db_Xferred) = InData%ta_rv(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%li) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%li,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%li,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%li,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%li,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%li,2), UBOUND(InData%li,2)
        DO i1 = LBOUND(InData%li,1), UBOUND(InData%li,1)
          DbKiBuf(Db_Xferred) = InData%li(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%M_uS) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%M_uS,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%M_uS,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%M_uS,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%M_uS,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%M_uS,2), UBOUND(InData%M_uS,2)
        DO i1 = LBOUND(InData%M_uS,1), UBOUND(InData%M_uS,1)
          DbKiBuf(Db_Xferred) = InData%M_uS(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%M_uD) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%M_uD,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%M_uD,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%M_uD,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%M_uD,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%M_uD,2), UBOUND(InData%M_uD,2)
        DO i1 = LBOUND(InData%M_uD,1), UBOUND(InData%M_uD,1)
          DbKiBuf(Db_Xferred) = InData%M_uD(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%M_f) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%M_f,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%M_f,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%M_f,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%M_f,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%M_f,2), UBOUND(InData%M_f,2)
        DO i1 = LBOUND(InData%M_f,1), UBOUND(InData%M_f,1)
          DbKiBuf(Db_Xferred) = InData%M_f(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
 END SUBROUTINE NWTC_Library_PackMeshMapLinearizationType

 SUBROUTINE NWTC_Library_UnPackMeshMapLinearizationType( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )
  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)
  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)
  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)
  TYPE(MeshMapLinearizationType), INTENT(INOUT) :: OutData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
    ! Local variables
  INTEGER(IntKi)                 :: Buf_size
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i
  INTEGER(IntKi)                 :: i1, i1_l, i1_u  !  bounds (upper/lower) for an array dimension 1
  INTEGER(IntKi)                 :: i2, i2_l, i2_u  !  bounds (upper/lower) for an array dimension 2
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'NWTC_Library_UnPackMeshMapLinearizationType'
 ! buffers to store meshes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred  = 1
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! mi not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%mi)) DEALLOCATE(OutData%mi)
    ALLOCATE(OutData%mi(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%mi.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%mi,2), UBOUND(OutData%mi,2)
        DO i1 = LBOUND(OutData%mi,1), UBOUND(OutData%mi,1)
          OutData%mi(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! fx_p not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%fx_p)) DEALLOCATE(OutData%fx_p)
    ALLOCATE(OutData%fx_p(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%fx_p.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%fx_p,2), UBOUND(OutData%fx_p,2)
        DO i1 = LBOUND(OutData%fx_p,1), UBOUND(OutData%fx_p,1)
          OutData%fx_p(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! tv_uD not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%tv_uD)) DEALLOCATE(OutData%tv_uD)
    ALLOCATE(OutData%tv_uD(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%tv_uD.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%tv_uD,2), UBOUND(OutData%tv_uD,2)
        DO i1 = LBOUND(OutData%tv_uD,1), UBOUND(OutData%tv_uD,1)
          OutData%tv_uD(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! tv_uS not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%tv_uS)) DEALLOCATE(OutData%tv_uS)
    ALLOCATE(OutData%tv_uS(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%tv_uS.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%tv_uS,2), UBOUND(OutData%tv_uS,2)
        DO i1 = LBOUND(OutData%tv_uS,1), UBOUND(OutData%tv_uS,1)
          OutData%tv_uS(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! ta_uD not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%ta_uD)) DEALLOCATE(OutData%ta_uD)
    ALLOCATE(OutData%ta_uD(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%ta_uD.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%ta_uD,2), UBOUND(OutData%ta_uD,2)
        DO i1 = LBOUND(OutData%ta_uD,1), UBOUND(OutData%ta_uD,1)
          OutData%ta_uD(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! ta_uS not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%ta_uS)) DEALLOCATE(OutData%ta_uS)
    ALLOCATE(OutData%ta_uS(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%ta_uS.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%ta_uS,2), UBOUND(OutData%ta_uS,2)
        DO i1 = LBOUND(OutData%ta_uS,1), UBOUND(OutData%ta_uS,1)
          OutData%ta_uS(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! ta_rv not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%ta_rv)) DEALLOCATE(OutData%ta_rv)
    ALLOCATE(OutData%ta_rv(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%ta_rv.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%ta_rv,2), UBOUND(OutData%ta_rv,2)
        DO i1 = LBOUND(OutData%ta_rv,1), UBOUND(OutData%ta_rv,1)
          OutData%ta_rv(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! li not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%li)) DEALLOCATE(OutData%li)
    ALLOCATE(OutData%li(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%li.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%li,2), UBOUND(OutData%li,2)
        DO i1 = LBOUND(OutData%li,1), UBOUND(OutData%li,1)
          OutData%li(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! M_uS not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%M_uS)) DEALLOCATE(OutData%M_uS)
    ALLOCATE(OutData%M_uS(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%M_uS.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%M_uS,2), UBOUND(OutData%M_uS,2)
        DO i1 = LBOUND(OutData%M_uS,1), UBOUND(OutData%M_uS,1)
          OutData%M_uS(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! M_uD not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%M_uD)) DEALLOCATE(OutData%M_uD)
    ALLOCATE(OutData%M_uD(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%M_uD.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%M_uD,2), UBOUND(OutData%M_uD,2)
        DO i1 = LBOUND(OutData%M_uD,1), UBOUND(OutData%M_uD,1)
          OutData%M_uD(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! M_f not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%M_f)) DEALLOCATE(OutData%M_f)
    ALLOCATE(OutData%M_f(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%M_f.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%M_f,2), UBOUND(OutData%M_f,2)
        DO i1 = LBOUND(OutData%M_f,1), UBOUND(OutData%M_f,1)
          OutData%M_f(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
 END SUBROUTINE NWTC_Library_UnPackMeshMapLinearizationType

 SUBROUTINE NWTC_Library_CopyMeshMapType( SrcMeshMapTypeData, DstMeshMapTypeData, CtrlCode, ErrStat, ErrMsg )
   TYPE(MeshMapType), INTENT(INOUT) :: SrcMeshMapTypeData
   TYPE(MeshMapType), INTENT(INOUT) :: DstMeshMapTypeData
   INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode
   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
   CHARACTER(*),    INTENT(  OUT) :: ErrMsg
! Local 
   INTEGER(IntKi)                 :: i,j,k
   INTEGER(IntKi)                 :: i1, i1_l, i1_u  !  bounds (upper/lower) for an array dimension 1
   INTEGER(IntKi)                 :: i2, i2_l, i2_u  !  bounds (upper/lower) for an array dimension 2
   INTEGER(IntKi)                 :: i3, i3_l, i3_u  !  bounds (upper/lower) for an array dimension 3
   INTEGER(IntKi)                 :: ErrStat2
   CHARACTER(ErrMsgLen)           :: ErrMsg2
   CHARACTER(*), PARAMETER        :: RoutineName = 'NWTC_Library_CopyMeshMapType'
! 
   ErrStat = ErrID_None
   ErrMsg  = ""
IF (ALLOCATED(SrcMeshMapTypeData%MapLoads)) THEN
  i1_l = LBOUND(SrcMeshMapTypeData%MapLoads,1)
  i1_u = UBOUND(SrcMeshMapTypeData%MapLoads,1)
  IF (.NOT. ALLOCATED(DstMeshMapTypeData%MapLoads)) THEN 
    ALLOCATE(DstMeshMapTypeData%MapLoads(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%MapLoads.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DO i1 = LBOUND(SrcMeshMapTypeData%MapLoads,1), UBOUND(SrcMeshMapTypeData%MapLoads,1)
      CALL NWTC_Library_Copymaptype( SrcMeshMapTypeData%MapLoads(i1), DstMeshMapTypeData%MapLoads(i1), CtrlCode, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)
         IF (ErrStat>=AbortErrLev) RETURN
    ENDDO
ENDIF
IF (ALLOCATED(SrcMeshMapTypeData%MapMotions)) THEN
  i1_l = LBOUND(SrcMeshMapTypeData%MapMotions,1)
  i1_u = UBOUND(SrcMeshMapTypeData%MapMotions,1)
  IF (.NOT. ALLOCATED(DstMeshMapTypeData%MapMotions)) THEN 
    ALLOCATE(DstMeshMapTypeData%MapMotions(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%MapMotions.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DO i1 = LBOUND(SrcMeshMapTypeData%MapMotions,1), UBOUND(SrcMeshMapTypeData%MapMotions,1)
      CALL NWTC_Library_Copymaptype( SrcMeshMapTypeData%MapMotions(i1), DstMeshMapTypeData%MapMotions(i1), CtrlCode, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)
         IF (ErrStat>=AbortErrLev) RETURN
    ENDDO
ENDIF
IF (ALLOCATED(SrcMeshMapTypeData%MapSrcToAugmt)) THEN
  i1_l = LBOUND(SrcMeshMapTypeData%MapSrcToAugmt,1)
  i1_u = UBOUND(SrcMeshMapTypeData%MapSrcToAugmt,1)
  IF (.NOT. ALLOCATED(DstMeshMapTypeData%MapSrcToAugmt)) THEN 
    ALLOCATE(DstMeshMapTypeData%MapSrcToAugmt(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%MapSrcToAugmt.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DO i1 = LBOUND(SrcMeshMapTypeData%MapSrcToAugmt,1), UBOUND(SrcMeshMapTypeData%MapSrcToAugmt,1)
      CALL NWTC_Library_Copymaptype( SrcMeshMapTypeData%MapSrcToAugmt(i1), DstMeshMapTypeData%MapSrcToAugmt(i1), CtrlCode, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)
         IF (ErrStat>=AbortErrLev) RETURN
    ENDDO
ENDIF
      CALL MeshCopy( SrcMeshMapTypeData%Augmented_Ln2_Src, DstMeshMapTypeData%Augmented_Ln2_Src, CtrlCode, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat>=AbortErrLev) RETURN
      CALL MeshCopy( SrcMeshMapTypeData%Lumped_Points_Src, DstMeshMapTypeData%Lumped_Points_Src, CtrlCode, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF (ErrStat>=AbortErrLev) RETURN
IF (ALLOCATED(SrcMeshMapTypeData%LoadLn2_A_Mat_Piv)) THEN
  i1_l = LBOUND(SrcMeshMapTypeData%LoadLn2_A_Mat_Piv,1)
  i1_u = UBOUND(SrcMeshMapTypeData%LoadLn2_A_Mat_Piv,1)
  IF (.NOT. ALLOCATED(DstMeshMapTypeData%LoadLn2_A_Mat_Piv)) THEN 
    ALLOCATE(DstMeshMapTypeData%LoadLn2_A_Mat_Piv(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%LoadLn2_A_Mat_Piv.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapTypeData%LoadLn2_A_Mat_Piv = SrcMeshMapTypeData%LoadLn2_A_Mat_Piv
ENDIF
IF (ALLOCATED(SrcMeshMapTypeData%DisplacedPosition)) THEN
  i1_l = LBOUND(SrcMeshMapTypeData%DisplacedPosition,1)
  i1_u = UBOUND(SrcMeshMapTypeData%DisplacedPosition,1)
  i2_l = LBOUND(SrcMeshMapTypeData%DisplacedPosition,2)
  i2_u = UBOUND(SrcMeshMapTypeData%DisplacedPosition,2)
  i3_l = LBOUND(SrcMeshMapTypeData%DisplacedPosition,3)
  i3_u = UBOUND(SrcMeshMapTypeData%DisplacedPosition,3)
  IF (.NOT. ALLOCATED(DstMeshMapTypeData%DisplacedPosition)) THEN 
    ALLOCATE(DstMeshMapTypeData%DisplacedPosition(i1_l:i1_u,i2_l:i2_u,i3_l:i3_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%DisplacedPosition.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapTypeData%DisplacedPosition = SrcMeshMapTypeData%DisplacedPosition
ENDIF
IF (ALLOCATED(SrcMeshMapTypeData%LoadLn2_A_Mat)) THEN
  i1_l = LBOUND(SrcMeshMapTypeData%LoadLn2_A_Mat,1)
  i1_u = UBOUND(SrcMeshMapTypeData%LoadLn2_A_Mat,1)
  i2_l = LBOUND(SrcMeshMapTypeData%LoadLn2_A_Mat,2)
  i2_u = UBOUND(SrcMeshMapTypeData%LoadLn2_A_Mat,2)
  IF (.NOT. ALLOCATED(DstMeshMapTypeData%LoadLn2_A_Mat)) THEN 
    ALLOCATE(DstMeshMapTypeData%LoadLn2_A_Mat(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%LoadLn2_A_Mat.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapTypeData%LoadLn2_A_Mat = SrcMeshMapTypeData%LoadLn2_A_Mat
ENDIF
IF (ALLOCATED(SrcMeshMapTypeData%LoadLn2_F)) THEN
  i1_l = LBOUND(SrcMeshMapTypeData%LoadLn2_F,1)
  i1_u = UBOUND(SrcMeshMapTypeData%LoadLn2_F,1)
  i2_l = LBOUND(SrcMeshMapTypeData%LoadLn2_F,2)
  i2_u = UBOUND(SrcMeshMapTypeData%LoadLn2_F,2)
  IF (.NOT. ALLOCATED(DstMeshMapTypeData%LoadLn2_F)) THEN 
    ALLOCATE(DstMeshMapTypeData%LoadLn2_F(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%LoadLn2_F.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapTypeData%LoadLn2_F = SrcMeshMapTypeData%LoadLn2_F
ENDIF
IF (ALLOCATED(SrcMeshMapTypeData%LoadLn2_M)) THEN
  i1_l = LBOUND(SrcMeshMapTypeData%LoadLn2_M,1)
  i1_u = UBOUND(SrcMeshMapTypeData%LoadLn2_M,1)
  i2_l = LBOUND(SrcMeshMapTypeData%LoadLn2_M,2)
  i2_u = UBOUND(SrcMeshMapTypeData%LoadLn2_M,2)
  IF (.NOT. ALLOCATED(DstMeshMapTypeData%LoadLn2_M)) THEN 
    ALLOCATE(DstMeshMapTypeData%LoadLn2_M(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
      CALL SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%LoadLn2_M.', ErrStat, ErrMsg,RoutineName)
      RETURN
    END IF
  END IF
    DstMeshMapTypeData%LoadLn2_M = SrcMeshMapTypeData%LoadLn2_M
ENDIF
      CALL NWTC_Library_Copymeshmaplinearizationtype( SrcMeshMapTypeData%dM, DstMeshMapTypeData%dM, CtrlCode, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)
         IF (ErrStat>=AbortErrLev) RETURN
 END SUBROUTINE NWTC_Library_CopyMeshMapType

 SUBROUTINE NWTC_Library_DestroyMeshMapType( MeshMapTypeData, ErrStat, ErrMsg )
  TYPE(MeshMapType), INTENT(INOUT) :: MeshMapTypeData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
  CHARACTER(*),    PARAMETER :: RoutineName = 'NWTC_Library_DestroyMeshMapType'
  INTEGER(IntKi)                 :: i, i1, i2, i3, i4, i5 
! 
  ErrStat = ErrID_None
  ErrMsg  = ""
IF (ALLOCATED(MeshMapTypeData%MapLoads)) THEN
DO i1 = LBOUND(MeshMapTypeData%MapLoads,1), UBOUND(MeshMapTypeData%MapLoads,1)
  CALL NWTC_Library_Destroymaptype( MeshMapTypeData%MapLoads(i1), ErrStat, ErrMsg )
ENDDO
  DEALLOCATE(MeshMapTypeData%MapLoads)
ENDIF
IF (ALLOCATED(MeshMapTypeData%MapMotions)) THEN
DO i1 = LBOUND(MeshMapTypeData%MapMotions,1), UBOUND(MeshMapTypeData%MapMotions,1)
  CALL NWTC_Library_Destroymaptype( MeshMapTypeData%MapMotions(i1), ErrStat, ErrMsg )
ENDDO
  DEALLOCATE(MeshMapTypeData%MapMotions)
ENDIF
IF (ALLOCATED(MeshMapTypeData%MapSrcToAugmt)) THEN
DO i1 = LBOUND(MeshMapTypeData%MapSrcToAugmt,1), UBOUND(MeshMapTypeData%MapSrcToAugmt,1)
  CALL NWTC_Library_Destroymaptype( MeshMapTypeData%MapSrcToAugmt(i1), ErrStat, ErrMsg )
ENDDO
  DEALLOCATE(MeshMapTypeData%MapSrcToAugmt)
ENDIF
  CALL MeshDestroy( MeshMapTypeData%Augmented_Ln2_Src, ErrStat, ErrMsg )
  CALL MeshDestroy( MeshMapTypeData%Lumped_Points_Src, ErrStat, ErrMsg )
IF (ALLOCATED(MeshMapTypeData%LoadLn2_A_Mat_Piv)) THEN
  DEALLOCATE(MeshMapTypeData%LoadLn2_A_Mat_Piv)
ENDIF
IF (ALLOCATED(MeshMapTypeData%DisplacedPosition)) THEN
  DEALLOCATE(MeshMapTypeData%DisplacedPosition)
ENDIF
IF (ALLOCATED(MeshMapTypeData%LoadLn2_A_Mat)) THEN
  DEALLOCATE(MeshMapTypeData%LoadLn2_A_Mat)
ENDIF
IF (ALLOCATED(MeshMapTypeData%LoadLn2_F)) THEN
  DEALLOCATE(MeshMapTypeData%LoadLn2_F)
ENDIF
IF (ALLOCATED(MeshMapTypeData%LoadLn2_M)) THEN
  DEALLOCATE(MeshMapTypeData%LoadLn2_M)
ENDIF
  CALL NWTC_Library_Destroymeshmaplinearizationtype( MeshMapTypeData%dM, ErrStat, ErrMsg )
 END SUBROUTINE NWTC_Library_DestroyMeshMapType

 SUBROUTINE NWTC_Library_PackMeshMapType( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )
  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)
  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)
  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)
  TYPE(MeshMapType),  INTENT(IN) :: InData
  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat
  CHARACTER(*),     INTENT(  OUT) :: ErrMsg
  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly
    ! Local variables
  INTEGER(IntKi)                 :: Re_BufSz
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_BufSz
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_BufSz
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i,i1,i2,i3,i4,i5
  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'NWTC_Library_PackMeshMapType'
 ! buffers to store subtypes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)

  OnlySize = .FALSE.
  IF ( PRESENT(SizeOnly) ) THEN
    OnlySize = SizeOnly
  ENDIF
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_BufSz  = 0
  Db_BufSz  = 0
  Int_BufSz  = 0
  Int_BufSz   = Int_BufSz   + 1     ! MapLoads allocated yes/no
  IF ( ALLOCATED(InData%MapLoads) ) THEN
    Int_BufSz   = Int_BufSz   + 2*1  ! MapLoads upper/lower bounds for each dimension
   ! Allocate buffers for subtypes, if any (we'll get sizes from these) 
    DO i1 = LBOUND(InData%MapLoads,1), UBOUND(InData%MapLoads,1)
      Int_BufSz   = Int_BufSz + 3  ! MapLoads: size of buffers for each call to pack subtype
      CALL NWTC_Library_Packmaptype( Re_Buf, Db_Buf, Int_Buf, InData%MapLoads(i1), ErrStat2, ErrMsg2, .TRUE. ) ! MapLoads 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! MapLoads
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! MapLoads
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! MapLoads
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    END DO
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! MapMotions allocated yes/no
  IF ( ALLOCATED(InData%MapMotions) ) THEN
    Int_BufSz   = Int_BufSz   + 2*1  ! MapMotions upper/lower bounds for each dimension
    DO i1 = LBOUND(InData%MapMotions,1), UBOUND(InData%MapMotions,1)
      Int_BufSz   = Int_BufSz + 3  ! MapMotions: size of buffers for each call to pack subtype
      CALL NWTC_Library_Packmaptype( Re_Buf, Db_Buf, Int_Buf, InData%MapMotions(i1), ErrStat2, ErrMsg2, .TRUE. ) ! MapMotions 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! MapMotions
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! MapMotions
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! MapMotions
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    END DO
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! MapSrcToAugmt allocated yes/no
  IF ( ALLOCATED(InData%MapSrcToAugmt) ) THEN
    Int_BufSz   = Int_BufSz   + 2*1  ! MapSrcToAugmt upper/lower bounds for each dimension
    DO i1 = LBOUND(InData%MapSrcToAugmt,1), UBOUND(InData%MapSrcToAugmt,1)
      Int_BufSz   = Int_BufSz + 3  ! MapSrcToAugmt: size of buffers for each call to pack subtype
      CALL NWTC_Library_Packmaptype( Re_Buf, Db_Buf, Int_Buf, InData%MapSrcToAugmt(i1), ErrStat2, ErrMsg2, .TRUE. ) ! MapSrcToAugmt 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! MapSrcToAugmt
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! MapSrcToAugmt
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! MapSrcToAugmt
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
    END DO
  END IF
      Int_BufSz   = Int_BufSz + 3  ! Augmented_Ln2_Src: size of buffers for each call to pack subtype
      CALL MeshPack( InData%Augmented_Ln2_Src, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2, .TRUE. ) ! Augmented_Ln2_Src 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! Augmented_Ln2_Src
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! Augmented_Ln2_Src
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! Augmented_Ln2_Src
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
      Int_BufSz   = Int_BufSz + 3  ! Lumped_Points_Src: size of buffers for each call to pack subtype
      CALL MeshPack( InData%Lumped_Points_Src, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2, .TRUE. ) ! Lumped_Points_Src 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! Lumped_Points_Src
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! Lumped_Points_Src
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! Lumped_Points_Src
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
  Int_BufSz   = Int_BufSz   + 1     ! LoadLn2_A_Mat_Piv allocated yes/no
  IF ( ALLOCATED(InData%LoadLn2_A_Mat_Piv) ) THEN
    Int_BufSz   = Int_BufSz   + 2*1  ! LoadLn2_A_Mat_Piv upper/lower bounds for each dimension
      Int_BufSz  = Int_BufSz  + SIZE(InData%LoadLn2_A_Mat_Piv)  ! LoadLn2_A_Mat_Piv
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! DisplacedPosition allocated yes/no
  IF ( ALLOCATED(InData%DisplacedPosition) ) THEN
    Int_BufSz   = Int_BufSz   + 2*3  ! DisplacedPosition upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%DisplacedPosition)  ! DisplacedPosition
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! LoadLn2_A_Mat allocated yes/no
  IF ( ALLOCATED(InData%LoadLn2_A_Mat) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! LoadLn2_A_Mat upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%LoadLn2_A_Mat)  ! LoadLn2_A_Mat
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! LoadLn2_F allocated yes/no
  IF ( ALLOCATED(InData%LoadLn2_F) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! LoadLn2_F upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%LoadLn2_F)  ! LoadLn2_F
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! LoadLn2_M allocated yes/no
  IF ( ALLOCATED(InData%LoadLn2_M) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! LoadLn2_M upper/lower bounds for each dimension
      Db_BufSz   = Db_BufSz   + SIZE(InData%LoadLn2_M)  ! LoadLn2_M
  END IF
      Int_BufSz   = Int_BufSz + 3  ! dM: size of buffers for each call to pack subtype
      CALL NWTC_Library_Packmeshmaplinearizationtype( Re_Buf, Db_Buf, Int_Buf, InData%dM, ErrStat2, ErrMsg2, .TRUE. ) ! dM 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN ! dM
         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )
         DEALLOCATE(Re_Buf)
      END IF
      IF(ALLOCATED(Db_Buf)) THEN ! dM
         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )
         DEALLOCATE(Db_Buf)
      END IF
      IF(ALLOCATED(Int_Buf)) THEN ! dM
         Int_BufSz = Int_BufSz + SIZE( Int_Buf )
         DEALLOCATE(Int_Buf)
      END IF
  IF ( Re_BufSz  .GT. 0 ) THEN 
     ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Db_BufSz  .GT. 0 ) THEN 
     ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF ( Int_BufSz  .GT. 0 ) THEN 
     ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )
     IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMsg,RoutineName)
       RETURN
     END IF
  END IF
  IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)

  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred = 1

  IF ( .NOT. ALLOCATED(InData%MapLoads) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%MapLoads,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%MapLoads,1)
    Int_Xferred = Int_Xferred + 2

    DO i1 = LBOUND(InData%MapLoads,1), UBOUND(InData%MapLoads,1)
      CALL NWTC_Library_Packmaptype( Re_Buf, Db_Buf, Int_Buf, InData%MapLoads(i1), ErrStat2, ErrMsg2, OnlySize ) ! MapLoads 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%MapMotions) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%MapMotions,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%MapMotions,1)
    Int_Xferred = Int_Xferred + 2

    DO i1 = LBOUND(InData%MapMotions,1), UBOUND(InData%MapMotions,1)
      CALL NWTC_Library_Packmaptype( Re_Buf, Db_Buf, Int_Buf, InData%MapMotions(i1), ErrStat2, ErrMsg2, OnlySize ) ! MapMotions 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%MapSrcToAugmt) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%MapSrcToAugmt,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%MapSrcToAugmt,1)
    Int_Xferred = Int_Xferred + 2

    DO i1 = LBOUND(InData%MapSrcToAugmt,1), UBOUND(InData%MapSrcToAugmt,1)
      CALL NWTC_Library_Packmaptype( Re_Buf, Db_Buf, Int_Buf, InData%MapSrcToAugmt(i1), ErrStat2, ErrMsg2, OnlySize ) ! MapSrcToAugmt 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
    END DO
  END IF
      CALL MeshPack( InData%Augmented_Ln2_Src, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2, OnlySize ) ! Augmented_Ln2_Src 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      CALL MeshPack( InData%Lumped_Points_Src, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2, OnlySize ) ! Lumped_Points_Src 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
  IF ( .NOT. ALLOCATED(InData%LoadLn2_A_Mat_Piv) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%LoadLn2_A_Mat_Piv,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%LoadLn2_A_Mat_Piv,1)
    Int_Xferred = Int_Xferred + 2

      DO i1 = LBOUND(InData%LoadLn2_A_Mat_Piv,1), UBOUND(InData%LoadLn2_A_Mat_Piv,1)
        IntKiBuf(Int_Xferred) = InData%LoadLn2_A_Mat_Piv(i1)
        Int_Xferred = Int_Xferred + 1
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%DisplacedPosition) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%DisplacedPosition,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%DisplacedPosition,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%DisplacedPosition,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%DisplacedPosition,2)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%DisplacedPosition,3)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%DisplacedPosition,3)
    Int_Xferred = Int_Xferred + 2

      DO i3 = LBOUND(InData%DisplacedPosition,3), UBOUND(InData%DisplacedPosition,3)
        DO i2 = LBOUND(InData%DisplacedPosition,2), UBOUND(InData%DisplacedPosition,2)
          DO i1 = LBOUND(InData%DisplacedPosition,1), UBOUND(InData%DisplacedPosition,1)
            DbKiBuf(Db_Xferred) = InData%DisplacedPosition(i1,i2,i3)
            Db_Xferred = Db_Xferred + 1
          END DO
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%LoadLn2_A_Mat) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%LoadLn2_A_Mat,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%LoadLn2_A_Mat,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%LoadLn2_A_Mat,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%LoadLn2_A_Mat,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%LoadLn2_A_Mat,2), UBOUND(InData%LoadLn2_A_Mat,2)
        DO i1 = LBOUND(InData%LoadLn2_A_Mat,1), UBOUND(InData%LoadLn2_A_Mat,1)
          DbKiBuf(Db_Xferred) = InData%LoadLn2_A_Mat(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%LoadLn2_F) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%LoadLn2_F,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%LoadLn2_F,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%LoadLn2_F,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%LoadLn2_F,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%LoadLn2_F,2), UBOUND(InData%LoadLn2_F,2)
        DO i1 = LBOUND(InData%LoadLn2_F,1), UBOUND(InData%LoadLn2_F,1)
          DbKiBuf(Db_Xferred) = InData%LoadLn2_F(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( .NOT. ALLOCATED(InData%LoadLn2_M) ) THEN
    IntKiBuf( Int_Xferred ) = 0
    Int_Xferred = Int_Xferred + 1
  ELSE
    IntKiBuf( Int_Xferred ) = 1
    Int_Xferred = Int_Xferred + 1
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%LoadLn2_M,1)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%LoadLn2_M,1)
    Int_Xferred = Int_Xferred + 2
    IntKiBuf( Int_Xferred    ) = LBOUND(InData%LoadLn2_M,2)
    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%LoadLn2_M,2)
    Int_Xferred = Int_Xferred + 2

      DO i2 = LBOUND(InData%LoadLn2_M,2), UBOUND(InData%LoadLn2_M,2)
        DO i1 = LBOUND(InData%LoadLn2_M,1), UBOUND(InData%LoadLn2_M,1)
          DbKiBuf(Db_Xferred) = InData%LoadLn2_M(i1,i2)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
      CALL NWTC_Library_Packmeshmaplinearizationtype( Re_Buf, Db_Buf, Int_Buf, InData%dM, ErrStat2, ErrMsg2, OnlySize ) ! dM 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf
        Re_Xferred = Re_Xferred + SIZE(Re_Buf)
        DEALLOCATE(Re_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Db_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf
        Db_Xferred = Db_Xferred + SIZE(Db_Buf)
        DEALLOCATE(Db_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
      IF(ALLOCATED(Int_Buf)) THEN
        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1
        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf
        Int_Xferred = Int_Xferred + SIZE(Int_Buf)
        DEALLOCATE(Int_Buf)
      ELSE
        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1
      ENDIF
 END SUBROUTINE NWTC_Library_PackMeshMapType

 SUBROUTINE NWTC_Library_UnPackMeshMapType( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )
  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)
  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)
  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)
  TYPE(MeshMapType), INTENT(INOUT) :: OutData
  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat
  CHARACTER(*),    INTENT(  OUT) :: ErrMsg
    ! Local variables
  INTEGER(IntKi)                 :: Buf_size
  INTEGER(IntKi)                 :: Re_Xferred
  INTEGER(IntKi)                 :: Db_Xferred
  INTEGER(IntKi)                 :: Int_Xferred
  INTEGER(IntKi)                 :: i
  INTEGER(IntKi)                 :: i1, i1_l, i1_u  !  bounds (upper/lower) for an array dimension 1
  INTEGER(IntKi)                 :: i2, i2_l, i2_u  !  bounds (upper/lower) for an array dimension 2
  INTEGER(IntKi)                 :: i3, i3_l, i3_u  !  bounds (upper/lower) for an array dimension 3
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(ErrMsgLen)           :: ErrMsg2
  CHARACTER(*), PARAMETER        :: RoutineName = 'NWTC_Library_UnPackMeshMapType'
 ! buffers to store meshes, if any
  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)
  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)
  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)
    !
  ErrStat = ErrID_None
  ErrMsg  = ""
  Re_Xferred  = 1
  Db_Xferred  = 1
  Int_Xferred  = 1
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! MapLoads not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%MapLoads)) DEALLOCATE(OutData%MapLoads)
    ALLOCATE(OutData%MapLoads(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%MapLoads.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    DO i1 = LBOUND(OutData%MapLoads,1), UBOUND(OutData%MapLoads,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL NWTC_Library_Unpackmaptype( Re_Buf, Db_Buf, Int_Buf, OutData%MapLoads(i1), ErrStat2, ErrMsg2 ) ! MapLoads 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! MapMotions not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%MapMotions)) DEALLOCATE(OutData%MapMotions)
    ALLOCATE(OutData%MapMotions(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%MapMotions.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    DO i1 = LBOUND(OutData%MapMotions,1), UBOUND(OutData%MapMotions,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL NWTC_Library_Unpackmaptype( Re_Buf, Db_Buf, Int_Buf, OutData%MapMotions(i1), ErrStat2, ErrMsg2 ) ! MapMotions 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! MapSrcToAugmt not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%MapSrcToAugmt)) DEALLOCATE(OutData%MapSrcToAugmt)
    ALLOCATE(OutData%MapSrcToAugmt(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%MapSrcToAugmt.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    DO i1 = LBOUND(OutData%MapSrcToAugmt,1), UBOUND(OutData%MapSrcToAugmt,1)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL NWTC_Library_Unpackmaptype( Re_Buf, Db_Buf, Int_Buf, OutData%MapSrcToAugmt(i1), ErrStat2, ErrMsg2 ) ! MapSrcToAugmt 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
    END DO
  END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL MeshUnpack( OutData%Augmented_Ln2_Src, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2 ) ! Augmented_Ln2_Src 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL MeshUnpack( OutData%Lumped_Points_Src, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2 ) ! Lumped_Points_Src 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! LoadLn2_A_Mat_Piv not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%LoadLn2_A_Mat_Piv)) DEALLOCATE(OutData%LoadLn2_A_Mat_Piv)
    ALLOCATE(OutData%LoadLn2_A_Mat_Piv(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%LoadLn2_A_Mat_Piv.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i1 = LBOUND(OutData%LoadLn2_A_Mat_Piv,1), UBOUND(OutData%LoadLn2_A_Mat_Piv,1)
        OutData%LoadLn2_A_Mat_Piv(i1) = IntKiBuf(Int_Xferred)
        Int_Xferred = Int_Xferred + 1
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! DisplacedPosition not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i3_l = IntKiBuf( Int_Xferred    )
    i3_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%DisplacedPosition)) DEALLOCATE(OutData%DisplacedPosition)
    ALLOCATE(OutData%DisplacedPosition(i1_l:i1_u,i2_l:i2_u,i3_l:i3_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%DisplacedPosition.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i3 = LBOUND(OutData%DisplacedPosition,3), UBOUND(OutData%DisplacedPosition,3)
        DO i2 = LBOUND(OutData%DisplacedPosition,2), UBOUND(OutData%DisplacedPosition,2)
          DO i1 = LBOUND(OutData%DisplacedPosition,1), UBOUND(OutData%DisplacedPosition,1)
            OutData%DisplacedPosition(i1,i2,i3) = REAL(DbKiBuf(Db_Xferred), R8Ki)
            Db_Xferred = Db_Xferred + 1
          END DO
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! LoadLn2_A_Mat not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%LoadLn2_A_Mat)) DEALLOCATE(OutData%LoadLn2_A_Mat)
    ALLOCATE(OutData%LoadLn2_A_Mat(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%LoadLn2_A_Mat.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%LoadLn2_A_Mat,2), UBOUND(OutData%LoadLn2_A_Mat,2)
        DO i1 = LBOUND(OutData%LoadLn2_A_Mat,1), UBOUND(OutData%LoadLn2_A_Mat,1)
          OutData%LoadLn2_A_Mat(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! LoadLn2_F not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%LoadLn2_F)) DEALLOCATE(OutData%LoadLn2_F)
    ALLOCATE(OutData%LoadLn2_F(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%LoadLn2_F.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%LoadLn2_F,2), UBOUND(OutData%LoadLn2_F,2)
        DO i1 = LBOUND(OutData%LoadLn2_F,1), UBOUND(OutData%LoadLn2_F,1)
          OutData%LoadLn2_F(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! LoadLn2_M not allocated
    Int_Xferred = Int_Xferred + 1
  ELSE
    Int_Xferred = Int_Xferred + 1
    i1_l = IntKiBuf( Int_Xferred    )
    i1_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    i2_l = IntKiBuf( Int_Xferred    )
    i2_u = IntKiBuf( Int_Xferred + 1)
    Int_Xferred = Int_Xferred + 2
    IF (ALLOCATED(OutData%LoadLn2_M)) DEALLOCATE(OutData%LoadLn2_M)
    ALLOCATE(OutData%LoadLn2_M(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%LoadLn2_M.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
      DO i2 = LBOUND(OutData%LoadLn2_M,2), UBOUND(OutData%LoadLn2_M,2)
        DO i1 = LBOUND(OutData%LoadLn2_M,1), UBOUND(OutData%LoadLn2_M,1)
          OutData%LoadLn2_M(i1,i2) = REAL(DbKiBuf(Db_Xferred), R8Ki)
          Db_Xferred = Db_Xferred + 1
        END DO
      END DO
  END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )
        Re_Xferred = Re_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )
        Db_Xferred = Db_Xferred + Buf_size
      END IF
      Buf_size=IntKiBuf( Int_Xferred )
      Int_Xferred = Int_Xferred + 1
      IF(Buf_size > 0) THEN
        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)
        IF (ErrStat2 /= 0) THEN 
           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)
           RETURN
        END IF
        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )
        Int_Xferred = Int_Xferred + Buf_size
      END IF
      CALL NWTC_Library_Unpackmeshmaplinearizationtype( Re_Buf, Db_Buf, Int_Buf, OutData%dM, ErrStat2, ErrMsg2 ) ! dM 
        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        IF (ErrStat >= AbortErrLev) RETURN

      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )
      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )
      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)
 END SUBROUTINE NWTC_Library_UnPackMeshMapType

!*********************************************************************************************************************************
!ENDOFREGISTRYGENERATEDFILE


!----------------------------------------------------------------------------------------------------------------------------------
END MODULE ModMesh_Mapping
!----------------------------------------------------------------------------------------------------------------------------------
