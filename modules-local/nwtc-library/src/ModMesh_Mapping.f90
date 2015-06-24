!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2015  National Renewable Energy Laboratory
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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
! This code implements the spatial mapping algorithms described in 
!    Sprague, Michael A.; Jonkman, Jason M.; and Jonkman, Bonnie J., "FAST Modular Wind Turbine CAE Tool: Nonmatching Spatial and  
!    Temporal Meshes." Proceedings of the 52nd Aerospace Sciences Meeting, 2014, also published in tech report NREL/CP-2C00-60742
!    National Renewable Energy Laboratory, Golden, CO. http://www.nrel.gov/docs/fy14osti/60742.pdf
!**********************************************************************************************************************************
MODULE ModMesh_Mapping

   USE ModMesh
   USE NWTC_LAPACK

   IMPLICIT NONE

   PRIVATE

   !bjj: these types require the use of ModMesh.f90, thus they cannot be part of NWTC_Library_Types.f90 (though they are auto_generated with that code):

   TYPE, PUBLIC :: MapType
      INTEGER(IntKi) :: OtherMesh_Element    ! Node (for point meshes) or Element (for line2 meshes) number on other mesh; for loads, other mesh is Dest, for motions/scalars, other mesh is Src
      REAL(ReKi)     :: distance             ! magnitude of couple_arm
      REAL(ReKi)     :: couple_arm(3)        ! Vector between a point and node 1 of an element (p_ODR - p_OSR)
      REAL(ReKi)     :: shape_fn(2)          ! shape functions: 1-D element-level location [0,1] based on closest-line projection of point
   END TYPE MapType

   TYPE, PUBLIC :: MeshMapType
      TYPE(MapType),  ALLOCATABLE :: MapLoads(:)
      TYPE(MapType),  ALLOCATABLE :: MapMotions(:)
      TYPE(MapType),  ALLOCATABLE :: MapSrcToAugmt(:)         ! for source line2 loads, we map between source and an augmented source mesh, then betweeN augmented source and destination
      TYPE(MeshType)              :: Augmented_Ln2_Src
      TYPE(MeshType)              :: Lumped_Points_Src        ! Stored here for efficiency
#ifdef MESH_DEBUG     
      TYPE(MeshType)              :: Lumped_Points_Dest        
#endif
      INTEGER,        ALLOCATABLE :: LoadLn2_A_Mat_Piv(:)      ! The pivot values for the factorizatioin of LoadLn2_A_Mat
      REAL(ReKi),     ALLOCATABLE :: DisplacedPosition(:,:,:)  ! couple_arm +Scr%Disp - Dest%Disp for each mapped node (stored here for efficiency.)
      REAL(ReKi),     ALLOCATABLE :: LoadLn2_A_Mat(:,:)        ! The 6-by-6 matrix that makes up the diagonal of the [A 0; B A] matrix in the point-to-line load mapping
      REAL(ReKi),     ALLOCATABLE :: LoadLn2_F(:,:)            ! The 3-components of the forces for each node of an element in the point-to-line load mapping (for each element)
      REAL(ReKi),     ALLOCATABLE :: LoadLn2_M(:,:)            ! The 3-components of the moments for each node of an element in the point-to-line load mapping (for each element)
   END TYPE

      ! note that these parameters must be negative (positive indicates the node/element it is mapped to)
   INTEGER(IntKi),  PARAMETER   :: NODE_NOT_MAPPED = -1

   PUBLIC :: MeshMapCreate
   PUBLIC :: MeshMapDestroy
   PUBLIC :: MeshMapWrBin
   PUBLIC :: Transfer_Point_to_Point
   PUBLIC :: Transfer_Line2_to_Point
   PUBLIC :: Transfer_Point_to_Line2
   PUBLIC :: Transfer_Line2_to_Line2
   !PUBLIC :: Lump_Line2_to_Point
   
   ! auto-generated routines, necessary for the FAST Registry:
   PUBLIC :: NWTC_Library_DestroyMeshMapType, NWTC_Library_CopyMeshMapType, NWTC_Library_PackMeshMapType, NWTC_Library_UnpackMeshMapType
   PUBLIC :: NWTC_Library_DestroyMapType,     NWTC_Library_CopyMapType,     NWTC_Library_PackMapType,     NWTC_Library_UnpackMapType

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!bjj: maybe the MeshMapCreate routine shouldn't actually allocate arrays; allocate them in
!   the "IF (RemapFlag)" sections so that if people add nodes during the simulation, the structures get reallocated to correct 
!   size? allocMapping should maybe be MeshMapping_Init() and only check that fields are compatible, etc.  
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE MeshMapCreate( Src, Dest, MeshMap, ErrStat, ErrMsg )
! This subroutine takes two meshes, determines the sizes required for the mapping data structure, and then
! allocates the mappings (different for loads and motions/scalars).
!..................................................................................................................................

! note that MeshMap%MapSrcToAugmt is allocated in Create_Augmented_Ln2_Src_Mesh() along with the Augmented_Ln2_Src Mesh

   TYPE(MeshType),           INTENT(IN)     ::  Src
   TYPE(MeshType),           INTENT(IN)     ::  Dest

   TYPE(MeshMapType),        INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),           INTENT(  OUT)  :: ErrStat      ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT)  :: ErrMsg       ! Error message if ErrStat /= ErrID_None


      ! local variables:

   INTEGER(IntKi)                           :: PointsInMap, PointsInTmpMap
   INTEGER(IntKi)                           :: ElementNodes
   LOGICAL                                  :: MapCreated
   INTEGER(IntKi)                           :: ErrStat2
   CHARACTER(ErrMsgLen)                     :: ErrMsg2
   

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
         CALL SetErrStat( ErrID_Fatal, 'MeshMap%MapMotions not allocated because no nodes were found to map.', ErrStat, ErrMsg, 'MeshMapCreate')
      ELSE

            ! Allocate the mapping structure:
         ALLOCATE( MeshMap%MapMotions(PointsInMap), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error trying to allocate MeshMap%MapMotions.', ErrStat, ErrMsg, 'MeshMapCreate')
         ELSE
            MapCreated = .TRUE.
            
               ! set up the initial mappings so that we don't necessarially have to do this multiple times on the first time step (if calculating Jacobians)
            IF ( Dest%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! point-to-Line2 or Line2-to-Line2
                  
               IF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! Line2-to-Line2         
                  CALL CreateMotionMap_L2_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapCreate')
               ELSEIF ( Src%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-Line2
                  CALL CreateMotionMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapCreate')
               END IF
         
            ELSEIF ( Dest%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-point or Line2-to-point
         
               IF ( Src%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-point            
                  CALL CreateMotionMap_P_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapCreate')            
               ELSEIF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! Line2-to-point         
                  CALL CreateMotionMap_L2_to_P(Src, Dest, MeshMap, ErrStat2, ErrMsg2)         
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapCreate')            
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
            CALL SetErrStat( ErrID_Fatal, 'Destination mesh does not contain force but source mesh does.', ErrStat, ErrMsg, 'MeshMapCreate')
         END IF
         IF (.NOT. Dest%FieldMask(MASKID_Moment) ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Destination mesh must contain moment when source mesh contains force.', ErrStat, ErrMsg, 'MeshMapCreate')
         END IF
      END IF
      IF ( Src%FieldMask(MASKID_Moment) ) THEN
         IF (.NOT. Dest%FieldMask(MASKID_Moment) ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Destination mesh does not contain moment but source mesh does.', ErrStat, ErrMsg, 'MeshMapCreate')
         END IF
      END IF

      
      ! get size of mapping:
      PointsInMap = Src%Nnodes 

      IF ( PointsInMap < 1 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'MeshMap%MapLoads not allocated because no nodes were found to map.', ErrStat, ErrMsg, 'MeshMapCreate')
      ELSE

            ! Allocate the mapping structure:
         ALLOCATE( MeshMap%MapLoads(PointsInMap), STAT=ErrStat2 )
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error trying to allocate MeshMap%MapLoads.', ErrStat, ErrMsg, 'MeshMapCreate')
         ELSE
            MapCreated = .TRUE.
            
               ! set up the initial mappings so that we don't necessarially have to do this multiple times on the first time step (if calculating Jacobians)
            IF ( Dest%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! point-to-Line2 or Line2-to-Line2
         
               ElementNodes = 2
         
               IF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! Line2-to-Line2         
                  CALL CreateLoadMap_L2_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapCreate')
               ELSEIF ( Src%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-Line2
                  CALL CreateLoadMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapCreate')
               END IF
         
            ELSEIF ( Dest%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-point or Line2-to-point
         
               IF ( Src%ElemTable(ELEMENT_POINT)%nelem > 0 ) THEN ! point-to-point            
                  CALL CreateLoadMap_P_to_P( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapCreate')            
               ELSEIF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN ! Line2-to-point         
                  CALL CreateLoadMap_L2_to_P(Src, Dest, MeshMap, ErrStat2, ErrMsg2)         
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapCreate')            
               END IF         
                  
            END IF ! create initial mapping based on mesh element type
                        
         END IF ! MapLoads allocated

      END IF ! Src has nodes to transfer
                  
   END IF ! HasLoadFields

   IF ( .NOT. MapCreated ) THEN
      CALL SetErrStat( ErrID_Fatal, 'Neither MapMotions or MapLoads was allocated. Meshes may not have compatible fields for mapping.', ErrStat, ErrMsg, 'MeshMapCreate')
      RETURN
   END IF


      !................................................
      ! Allocate the DisplacedPosition field:
      !................................................

   IF (.NOT. ALLOCATED (MeshMap%DisplacedPosition)) THEN
      CALL AllocAry( MeshMap%DisplacedPosition, 3, PointsInTmpMap, ElementNodes, 'MeshMap%DisplacedPosition', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MeshMapCreate')
   END IF


END SUBROUTINE MeshMapCreate
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE MeshMapDestroy( MeshMap, ErrStat, ErrMsg )
! This routine destroys the elements of the MeshMapType
!..................................................................................................................................

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
SUBROUTINE MeshMapWrBin( UnIn, Src, Dest, MeshMap, ErrStat, ErrMsg, FileName )
! This creates a matrix to write to a binary file (for debugging)
!..................................................................................................................................
   INTEGER,    INTENT(INOUT)                ::  UnIn     ! fortran output unit


   TYPE(MeshType),           INTENT(IN   )  :: Src
   TYPE(MeshType),           INTENT(IN   )  :: Dest
   TYPE(MeshMapType),        INTENT(IN   )  :: MeshMap

   INTEGER(IntKi),           INTENT(  OUT)  :: ErrStat      ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT)  :: ErrMsg       ! Error message if ErrStat /= ErrID_None
   CHARACTER(*),       INTENT(IN), OPTIONAL :: FileName  ! Name of the file to write the output in

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
SUBROUTINE Transfer_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

! for transfering displacement/velocity/accerlation data from Line2 mesh to Point Mesh
!
! Also works for transfering
!
   TYPE(MeshType),         INTENT(IN   ) ::  Src
   TYPE(MeshType),         INTENT(INOUT) ::  Dest
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  SrcDisp  ! a "functional" sibling of the source mesh; Src contains loads and SrcDisp contains TranslationDisp and Orientaiton
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  DestDisp ! a "functional" sibling of the destination mesh; Dest contains loads and DestDisp contains TranslationDisp and Orientaiton

   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   REAL(ReKi)                              :: LoadsScaleFactor  ! bjj: added this scaling factor to get loads in a better numerical range 
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2

   
   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

   !.................
   ! Check to ensure that the source mesh is composed of Line2 elements and destination mesh is composed of Point elements
   !.................   
   
   if (Src%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Source mesh must have one or more Line2 elements.', ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
      RETURN
   endif

   if (Dest%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Destination mesh must have one or more Point elements.', ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
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
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
            IF (ErrStat >= AbortErrLev) RETURN

      endif !remapping


      !........................
      ! Start: Transfer data
      !........................

      CALL Transfer_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
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
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
            IF (ErrStat >= AbortErrLev) RETURN
         

      ELSE ! Check that the temporary mesh has been set

         IF ( .NOT. MeshMap%Lumped_Points_Src%Initialized ) THEN
            CALL SetErrStat( ErrID_Fatal, 'MeshMap%Lumped_Points_Src not initialized (set RemapFlag = TRUE).', ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
            RETURN
         END IF
         
      END IF

      !........................
      ! Lump the values on line2 elements to nodes on a point mesh
      ! These point and line2 meshes have the same nodes so there
      ! is no mapping here:
      !........................

      IF ( PRESENT(SrcDisp) ) THEN
         
         LoadsScaleFactor = GetLoadsScaleFactor ( Src )
         
         ! first, we take the source fields and transfer them to fields on the augmented source mesh:
         !  (we're also taking the SrcDisp field and putting it on our augmented mesh)
         CALL Transfer_Src_To_Augmented_Ln2_Src( Src, MeshMap, ErrStat2, ErrMsg2, SrcDisp, LoadsScaleFactor ) 
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
                  
         ! then we lump the loads from the augmented source mesh:
         CALL Lump_Line2_to_Point( MeshMap%Augmented_Ln2_Src,  MeshMap%Lumped_Points_Src,  ErrStat2, ErrMsg2, LoadsScaleFactor=LoadsScaleFactor  ) 
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
            IF (ErrStat >= AbortErrLev) RETURN
         
      ELSE
         CALL SetErrStat( ErrID_Fatal, 'Invalid arguments to routine for transfer of loads.', ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
         RETURN
      END IF


      !........................
      ! Transfer data
      !........................

      IF ( PRESENT( DestDisp ) ) THEN ! note that we already checked if SrcDisp is present 
         
         ! and transferred the displacements to MeshMap%Augmented_Ln2_Src
         CALL Transfer_Loads_Point_to_Point( MeshMap%Lumped_Points_Src, Dest, MeshMap, ErrStat2, ErrMsg2, MeshMap%Augmented_Ln2_Src, DestDisp, LoadsScaleFactor )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
            IF (ErrStat >= AbortErrLev) RETURN
         
      ELSE
         CALL SetErrStat( ErrID_Fatal, 'Invalid arguments to routine for transfer of loads.', ErrStat, ErrMsg, 'Transfer_Line2_to_Point')
         RETURN
      END IF


   end if !algorithm for loads


END SUBROUTINE Transfer_Line2_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )
! Given a mapping, this routine transfers the motions from nodes on Line2 elements to nodes on another mesh.
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src       ! The source (Line2) mesh with motion fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest      ! The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   ! The mapping data

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)            :: i , j                          ! counter over the nodes
   INTEGER(IntKi)            :: n, n1, n2                      ! temporary space for node numbers
   REAL(ReKi)                :: FieldValueN1(3)                ! Temporary variable to store field values on element nodes
   REAL(ReKi)                :: FieldValueN2(3)                ! Temporary variable to store field values on element nodes
   REAL(ReKi)                :: FieldValue(3,2)                ! Temporary variable to store values for DCM interpolation
   REAL(ReKi)                :: TmpVec(3)
   REAL(ReKi)                :: RotationMatrix(3,3)


   ErrStat = ErrID_None
   ErrMsg  = ""



!bjj: FieldValueN1 and FieldValueN2 should really be one matrix of DIM (3,2) now that we've modified some of the other data structures....
             
      ! ---------------------------- Translation ------------------------------------------------

      ! u_Dest1 = u_Src + [Orientation_Src^T * RefOrientation_Src - I] * [p_Dest - p_Src] at Source Node n1
      ! u_Dest2 = u_Src + [Orientation_Src^T * RefOrientation_Src - I] * [p_Dest - p_Src] at Source Node n2
      ! u_Dest = (1.-elem_position)*u_Dest1 + elem_position*u_Dest2
   if ( Src%FieldMask(MASKID_TranslationDisp) .AND. Dest%FieldMask(MASKID_TranslationDisp) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

                                                               
         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         FieldValueN1 = Src%TranslationDisp(:,n1)
         FieldValueN2 = Src%TranslationDisp(:,n2)

            ! if Src mesh has orientation, superpose Dest displacement with translation due to rotation and couple arm
         IF ( Src%FieldMask(MASKID_Orientation) ) THEN

            ! augment translation on node n1:
            
               !Calculate RotationMatrix as O_S^T*O_SR
            RotationMatrix = TRANSPOSE( Src%Orientation(:,:,n1 ) )
            RotationMatrix = MATMUL( RotationMatrix, Src%RefOrientation(:,:,n1) )

               ! subtract I
            do j=1,3
               RotationMatrix(j,j)= RotationMatrix(j,j) - 1.0_ReKi
            end do

            FieldValueN1 = FieldValueN1 + MATMUL(RotationMatrix,(Dest%Position(:,i)-Src%Position(:,n1)))
            
            ! augment translation on node n2:
            
               !Calculate RotationMatrix as O_S^T*O_SR
            RotationMatrix = TRANSPOSE( Src%Orientation(:,:,n2 ) )
            RotationMatrix = MATMUL( RotationMatrix, Src%RefOrientation(:,:,n2) )

               ! subtract I
            do j=1,3
               RotationMatrix(j,j)= RotationMatrix(j,j) - 1.0_ReKi
            end do

            FieldValueN2 = FieldValueN2 + MATMUL(RotationMatrix,(Dest%Position(:,i)-Src%Position(:,n2)))
            
         end if

         Dest%TranslationDisp(:,i) = MeshMap%MapMotions(i)%shape_fn(1)*FieldValueN1  &
                                   + MeshMap%MapMotions(i)%shape_fn(2)*FieldValueN2

      end do

   end if

      ! ---------------------------- ORIENTATION/Direction Cosine Matrix   ----------------------

      ! transfer direction cosine matrix, aka orientation

   if ( Src%FieldMask(MASKID_Orientation) .AND. Dest%FieldMask(MASKID_Orientation) ) then

      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
                  
         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

#ifdef __NN_ORIENTATIONS         
         IF ( NINT( MeshMap%MapMotions(i)%shape_fn(2) ) == 0 ) THEN
            n = n1
         ELSE
            n = n2
         END IF
         
         Dest%Orientation(:,:,i) = MATMUL( MATMUL( Dest%RefOrientation(:,:,i), TRANSPOSE( Src%RefOrientation(:,:,n) ) )&
                                          , Src%Orientation(:,:,n) )
#else         
#ifdef __ORIGINAL_LOGMAP    
            ! calculate Rotation matrix for FieldValueN1 and convert to tensor:
         RotationMatrix = MATMUL( MATMUL( Dest%RefOrientation(:,:,i), TRANSPOSE( Src%RefOrientation(:,:,n1) ) )&
                                 , Src%Orientation(:,:,n1) )

         CALL DCM_logmap( RotationMatrix, FieldValue(:,1), ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN

            ! calculate Rotation matrix for FieldValueN2 and convert to tensor:
         RotationMatrix = MATMUL( MATMUL( Dest%RefOrientation(:,:,i), TRANSPOSE( Src%RefOrientation(:,:,n2) ) )&
                                 , Src%Orientation(:,:,n2) )
         
         CALL DCM_logmap( RotationMatrix, FieldValue(:,2), ErrStat, ErrMsg )                  
         IF (ErrStat >= AbortErrLev) RETURN
         
         CALL DCM_SetLogMapForInterp( FieldValue )  ! make sure we don't cross a 2pi boundary
         
         
            ! interpolate tensors: 
         TmpVec =   MeshMap%MapMotions(i)%shape_fn(1)*FieldValue(:,1)  &
                  + MeshMap%MapMotions(i)%shape_fn(2)*FieldValue(:,2)    
                  
            ! convert back to DCM:
         Dest%Orientation(:,:,i) = DCM_exp( TmpVec )               
      
         
#else         
!this should be equivalent, with one less matrix multiply
         
            ! calculate Rotation matrix for FieldValueN1 and convert to tensor:
         RotationMatrix = MATMUL( TRANSPOSE( Src%RefOrientation(:,:,n1) ), Src%Orientation(:,:,n1) )
         
         CALL DCM_logmap( RotationMatrix, FieldValue(:,1), ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN
         
            ! calculate Rotation matrix for FieldValueN2 and convert to tensor:
         RotationMatrix = MATMUL( TRANSPOSE( Src%RefOrientation(:,:,n2) ), Src%Orientation(:,:,n2) )
         
         CALL DCM_logmap( RotationMatrix, FieldValue(:,2), ErrStat, ErrMsg )                  
         IF (ErrStat >= AbortErrLev) RETURN
         
         CALL DCM_SetLogMapForInterp( FieldValue )  ! make sure we don't cross a 2pi boundary
         
         
            ! interpolate tensors: 
         TmpVec =   MeshMap%MapMotions(i)%shape_fn(1)*FieldValue(:,1)  &
                  + MeshMap%MapMotions(i)%shape_fn(2)*FieldValue(:,2)    
                  
            ! convert back to DCM:
         Dest%Orientation(:,:,i) = MATMUL( Dest%RefOrientation(:,:,i), DCM_exp( TmpVec )  )
         
#endif         
#endif
         
      end do

   endif

      ! ---------------------------- Calculated total displaced positions  ---------------------
      ! these values are used in both the translational velocity and translational acceleration
      ! calculations. The calculations rely on the TranslationDisp fields, which are calculated
      ! earlier in this routine.
   IF ( Src%FieldMask(MASKID_TranslationVel) .OR. Src%FieldMask(MASKID_TranslationAcc) ) THEN
      DO i = 1,Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
      
         DO j=1,NumNodes(ELEMENT_LINE2) ! number of nodes per line2 element
            n = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(j)
         
            MeshMap%DisplacedPosition(:,i,j) =    Src%Position(:,n) +  Src%TranslationDisp(:,n)  &
                                               - Dest%Position(:,i) - Dest%TranslationDisp(:,i)  
         end do
      
      END DO   
   END IF
   
      ! ---------------------------- TranslationVel  --------------------------------------------

   if ( Src%FieldMask(MASKID_TranslationVel) .AND. Dest%FieldMask(MASKID_TranslationVel) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

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

   if ( Src%FieldMask(MASKID_RotationVel) .AND. Dest%FieldMask(MASKID_RotationVel) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         Dest%RotationVel(:,i) = MeshMap%MapMotions(i)%shape_fn(1)*Src%RotationVel(:,n1)  &
                               + MeshMap%MapMotions(i)%shape_fn(2)*Src%RotationVel(:,n2)
      end do

   end if

      ! ---------------------------- TranslationAcc -----------------------------------------------

   if ( Src%FieldMask(MASKID_TranslationAcc) .AND. Dest%FieldMask(MASKID_TranslationAcc) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE


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

   if (Src%FieldMask(MASKID_RotationAcc) .AND. Dest%FieldMask(MASKID_RotationAcc) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         Dest%RotationAcc(:,i) = MeshMap%MapMotions(i)%shape_fn(1)*Src%RotationAcc(:,n1)  &
                               + MeshMap%MapMotions(i)%shape_fn(2)*Src%RotationAcc(:,n2)


      end do
   end if

      ! ---------------------------- Scalars  -----------------------------------------------

   if (Src%FieldMask(MASKID_SCALAR) .AND. Dest%FieldMask(MASKID_SCALAR) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         Dest%Scalars(:,i) = MeshMap%MapMotions(i)%shape_fn(1)*Src%Scalars(:,n1)  &
                           + MeshMap%MapMotions(i)%shape_fn(2)*Src%Scalars(:,n2)

      end do
   end if


END SUBROUTINE Transfer_Motions_Line2_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CreateMapping_ProjectToLine2(Mesh1, Mesh2, Map, Mesh1_TYPE, ErrStat, ErrMsg)
!This routine projects Mesh1 onto a Line2 mesh (Mesh2) to find the element mappings between the two meshes.
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Mesh1                           ! The mesh in the outer mapping loop (Dest for Motions/Scalars; Src for Loads)
   TYPE(MeshType),                 INTENT(IN   )  :: Mesh2                           ! The mesh in the inner mapping loop (Src for Motions/Scalars; Dest for Loads)

   TYPE(MapType),                  INTENT(INOUT)  :: Map(:)                          ! The mapping from Src to Dest

   INTEGER(IntKi),                 INTENT(IN   )  :: Mesh1_TYPE                      ! Type of Mesh1 elements to map
   INTEGER(IntKi),   PARAMETER                    :: Mesh2_TYPE  = ELEMENT_LINE2     ! Type of Mesh2 elements on map (MUST BE ELEMENT_LINE2)

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                        ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                         ! Error message if ErrStat /= ErrID_None
   

      ! local variables

   REAL(ReKi)      :: denom
   REAL(ReKi)      :: dist
   REAL(ReKi)      :: min_dist
   REAL(ReKi)      :: elem_position

   REAL(ReKi)      :: Mesh1_xyz(3)

   REAL(ReKi)      :: n1_n2_vector(3)     ! vector going from node 1 to node 2 in Line2 element
   REAL(ReKi)      :: n1_Point_vector(3)  ! vector going from node 1 in Line 2 element to Destination Point
   REAL(ReKi)      :: tmp(3)              ! temporary vector for cross product calculation


   INTEGER(IntKi)  :: iElem, iNode, i  ! do-loop counter for elements on Mesh1, associated node(S)
   INTEGER(IntKi)  :: jElem            ! do-loop counter for elements on Mesh2, associated node

   INTEGER(IntKi)  :: n1, n2           ! nodes associated with an element

   LOGICAL         :: found
   LOGICAL         :: on_element
   


      ! initialization
   ErrStat = ErrID_None
   ErrMsg  = ""


   ! Map the source nodes to destination nodes:
   do n1=1,size(Map)
      Map(n1)%OtherMesh_Element = NODE_NOT_MAPPED ! initialize this so we know if we've mapped this node already (done only because we may have different elements)
   end do !n1
      


   do iElem = 1, Mesh1%ElemTable(Mesh1_TYPE)%nelem   ! number of Mesh1_TYPE elements on Mesh1
      do iNode = 1, SIZE( Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes )
         i = Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes(iNode)  ! the nodes on element iElem
         IF ( Map(i)%OtherMesh_Element > 0 ) CYCLE  ! we already mapped this node; let's move on to the next iNode (or iElem)

         ! destination point
         Mesh1_xyz = Mesh1%Position(:, i)

         found = .false.
         min_dist = HUGE(min_dist)

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
               CALL SetErrStat( ErrID_Fatal, 'Division by zero because Line2 element nodes are in same position.', ErrStat, ErrMsg, 'CreateMapping_ProjectToLine2')
               RETURN
            END IF

               ! project point onto line defined by n1 and n2

            elem_position = DOT_PRODUCT(n1_n2_vector,n1_Point_vector) / denom

                  ! note: i forumlated it this way because Fortran doesn't necessarially do shortcutting and I don't want to call EqualRealNos if we don't need it:
            if ( elem_position .ge. 0.0_ReKi .and. elem_position .le. 1.0_ReKi ) then !we're ON the element (between the two nodes)
               on_element = .true.
            elseif (EqualRealNos( elem_position, 1.0_ReKi )) then !we're ON the element (at a node)
               on_element = .true.
            elseif (EqualRealNos( elem_position,  0.0_ReKi )) then !we're ON the element (at a node)
               on_element = .true.
            else !we're not on the element
               on_element = .false.
            end if

            if (on_element) then

               ! calculate distance between point and line (note: this is actually the distance squared);
               ! will only store information once we have determined the closest element
               tmp  = cross_product( n1_n2_vector, n1_Point_vector )
               dist = DOT_PRODUCT(tmp,tmp) / denom

               if (dist .lt. min_dist) then
                  found = .true.
                  min_dist = dist

                  Map(i)%OtherMesh_Element = jElem
                  Map(i)%shape_fn(1)       = 1.0_ReKi - elem_position
                  Map(i)%shape_fn(2)       = elem_position

                  !Map(i)%couple_arm        = n1_Point_vector

               end if !the point is closest to this line2 element

            endif

         end do !jElem

            ! if failed to find an element that the Point projected into, throw an error
         if (.not. found) then

            if (Map(i)%OtherMesh_Element .lt. 1 )  then
               CALL SetErrStat( ErrID_Fatal, 'node does not project onto any line2 element', ErrStat, ErrMsg, 'CreateMapping_ProjectToLine2')
               RETURN
            endif

         end if !not found on projection to element

      end do !iNode
   end do !iElem

END SUBROUTINE CreateMapping_ProjectToLine2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Create_PointMesh(Src, Temp_Point_Src, ErrStat, ErrMsg)
!This routine creates a new mesh with the same positions as the Src mesh, except all of the elements are points. It adds fields for
! forces, moments, and/or TranslationDisp, if they are part of the Src mesh
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(INOUT)  :: Temp_Point_Src                  ! A blank mesh to be

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

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
SUBROUTINE CreateLoadMap_L2_to_P( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            ! The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

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
SUBROUTINE CreateMotionMap_L2_to_P( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            ! The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

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
SUBROUTINE Transfer_Point_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

! for transfering displacement/velocity/accerlation data from Point Mesh to Line2 Mesh

! Example: SubDyn  Output: Point Mesh: Displacement/Orientation/TranslationVel, etc
!          HyDroDyn Input: Line2 Mesh: "

   TYPE(MeshType),         INTENT(IN   ) :: Src
   TYPE(MeshType),         INTENT(INOUT) :: Dest
   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   TYPE(MeshType),OPTIONAL,INTENT(IN   ) :: SrcDisp  ! a "functional" sibling of the source mesh; Src contains loads and SrcDisp contains TranslationDisp and Orientaiton
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) :: DestDisp ! a "functional" sibling of the destination mesh; Dest contains loads and DestDisp contains TranslationDisp and Orientaiton

   ! local variables

   REAL(ReKi)                            :: LoadsScaleFactor  ! bjj: added this scaling factor to get loads in a better numerical range 
   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2

   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

      IF (ErrStat >= AbortErrLev) RETURN
   !.................
   ! Check to ensure that the source mesh is composed of Point elements and destination mesh is composed of Line2 elements
   !.................   
   
   if (Src%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Source mesh must have one or more Point elements.', ErrStat, ErrMsg, 'Transfer_Point_to_Line2')
      RETURN
   endif
   
   if (Dest%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Destination mesh must have one or more Line2 elements.', ErrStat, ErrMsg, 'Transfer_Point_to_Line2')
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
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Point_to_Line2')
            IF (ErrStat >= AbortErrLev) RETURN
         
      end if

      !........................
      ! Transfer data
      !........................

         ! This is the same algorithm as Transfer_Point_to_Point
      CALL Transfer_Motions_Point_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Point_to_Line2')
         IF (ErrStat >= AbortErrLev) RETURN


   end if ! algorithm for motions/scalars

   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if ( HasLoadFields(Src) ) then

      IF (.not. PRESENT(SrcDisp) .OR. .NOT. PRESENT(DestDisp) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'SrcDisp and DestDisp arguments are required for load transfer.', ErrStat, ErrMsg, 'Transfer_Point_to_Line2')
         RETURN
      END IF

      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then
                           
         CALL CreateLoadMap_P_to_L2( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Point_to_Line2')
            IF (ErrStat >= AbortErrLev) RETURN
                           
      end if


      !........................
      ! Transfer data
      !........................
      LoadsScaleFactor = GetLoadsScaleFactor ( Src ) 

      CALL Transfer_Loads_Point_to_Line2( Src, Dest, MeshMap, ErrStat2, ErrMsg2, SrcDisp, DestDisp, LoadsScaleFactor )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Point_to_Line2')
         IF (ErrStat >= AbortErrLev) RETURN
      
     
   end if ! algorithm for loads


END SUBROUTINE Transfer_Point_to_Line2
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
SUBROUTINE Transfer_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

   TYPE(MeshType),          INTENT(IN   ) ::  Src
   TYPE(MeshType),          INTENT(INOUT) ::  Dest
   TYPE(MeshType), OPTIONAL,INTENT(IN   ) ::  SrcDisp  ! an optional mesh, which contains the displacements associated with the source if the source contains load information
   TYPE(MeshType), OPTIONAL,INTENT(IN   ) ::  DestDisp ! an optional mesh, which contains the displacements associated with the destination if the destination contains load information

   TYPE(MeshMapType),       INTENT(INOUT) ::  MeshMap     ! Mapping(s) between Src and Dest meshes

   INTEGER(IntKi),          INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),            INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None


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


   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   ! If Src is displacements/velocities then loop over the Destination Mesh;
   !    each motion/scalar in the destination mesh needs to be interpolated from the source mesh.


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

   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   ! IF Src is forces and/or moments, loop over Src Mesh;
   !    each load in the source mesh needs to be placed somewhere in the destination mesh.

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
SUBROUTINE CreateMapping_NearestNeighbor( Mesh1, Mesh2, Map, Mesh1_TYPE, Mesh2_TYPE, ErrStat, ErrMsg )
! This routine creates the node-to-node (nearest neighbor). We map FROM Mesh1 to Mesh2
!.......................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Mesh1      ! The mesh in the outer mapping loop (Dest for Motions/Scalars; Src for Loads)
   TYPE(MeshType),                 INTENT(IN   )  :: Mesh2      ! The mesh in the inner mapping loop (Src for Motions/Scalars; Dest for Loads)

   TYPE(MapType),                  INTENT(INOUT)  :: Map(:)      ! The mapping from Src to Dest

   INTEGER(IntKi),                 INTENT(IN   )  :: Mesh1_TYPE  ! Type of Mesh1 elements to map
   INTEGER(IntKi),                 INTENT(IN   )  :: Mesh2_TYPE  ! Type of Mesh2 elements on map

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

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
   do i=1,size(Map)
      Map(i)%OtherMesh_Element = NODE_NOT_MAPPED ! initialize this so we know if we've mapped this node already (done only because we may have different elements)
   end do !n1
   

   do iElem = 1, Mesh1%ElemTable(Mesh1_TYPE)%nelem   ! number of Mesh1_TYPE elements on Mesh1 = number of points on Mesh1
      do iNode = 1, SIZE( Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes )
         i = Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes(iNode)  ! the nodes on element iElem
         IF ( Map(i)%OtherMesh_Element > 0 ) CYCLE  ! we already mapped this node; let's move on


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

         Map(i)%OtherMesh_Element = point_with_min_dist

         Map(i)%distance = min_dist

         Map(i)%couple_arm = Mesh2%Position(:, point_with_min_dist) - Mesh1_xyz
         !bjj: this is the negative of the case where it's Mesh2=src, so we'll have to multiply by -1 outside this routine if that's the case

      end do !iNode
   end do !iElem


END SUBROUTINE CreateMapping_NearestNeighbor
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_Motions_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )
! Given a nearest-neighbor mapping, this routine transfers motions between nodes on the mesh.
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src       ! The source mesh with motion fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest      ! The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   ! data for the mesh mapping

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)  :: i, j                                     ! counter over the nodes
   REAL(ReKi)      :: RotationMatrix(3,3)
   REAL(ReKi)      :: TmpVec(3)


   ErrStat = ErrID_None
   ErrMsg  = ""



      ! ---------------------------- Translation ------------------------------------------------

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
               RotationMatrix(j,j)= RotationMatrix(j,j) - 1.0_ReKi
            end do


            Dest%TranslationDisp(:,i) = Dest%TranslationDisp(:,i) + MATMUL(RotationMatrix, MeshMap%MapMotions(i)%couple_arm)

         end if

      end do

   end if


      ! ---------------------------- ORIENTATION/Direction Cosine Matrix   ----------------------

      ! transfer direction cosine matrix, aka orientation

   if ( Src%FieldMask(MASKID_Orientation) .AND. Dest%FieldMask(MASKID_Orientation) ) then

      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
         
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
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE
         MeshMap%DisplacedPosition(:,i,1) =    Src%TranslationDisp(:,MeshMap%MapMotions(i)%OtherMesh_Element) &
                                            - Dest%TranslationDisp(:,i) &
                                            - MeshMap%MapMotions(i)%couple_arm
      END DO
      
   END IF
   
      ! ---------------------------- TranslationVel  --------------------------------------------

   if ( Src%FieldMask(MASKID_TranslationVel) .AND. Dest%FieldMask(MASKID_TranslationVel) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%TranslationVel(:,i) = Src%TranslationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element)

         if ( Src%FieldMask(MASKID_RotationVel) ) then
            Dest%TranslationVel(:,i) = Dest%TranslationVel(:,i) + &
                                       cross_product ( MeshMap%DisplacedPosition(:,i,1), &
                                                       Src%RotationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element) )
         endif
      end do

   endif

      ! ---------------------------- RotationVel  -----------------------------------------------

   if ( Src%FieldMask(MASKID_RotationVel) .AND. Dest%FieldMask(MASKID_RotationVel) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%RotationVel(:,i) = Src%RotationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element)
      end do
   end if

      ! ---------------------------- TranslationAcc -----------------------------------------------

   if ( Src%FieldMask(MASKID_TranslationAcc) .AND. Dest%FieldMask(MASKID_TranslationAcc) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

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

   if (Src%FieldMask(MASKID_RotationAcc) .AND. Dest%FieldMask(MASKID_RotationAcc) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%RotationAcc(:,i) = Src%RotationAcc(:,MeshMap%MapMotions(i)%OtherMesh_Element)
      end do
   end if

      ! ---------------------------- Scalars  -----------------------------------------------

   if (Src%FieldMask(MASKID_SCALAR) .AND. Dest%FieldMask(MASKID_SCALAR) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%Scalars(:,i) = Src%Scalars(:,MeshMap%MapMotions(i)%OtherMesh_Element)
      end do
   end if


END SUBROUTINE Transfer_Motions_Point_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_Loads_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp,LoadsScaleFactor )
! Given a nearest-neighbor mapping, this routine transfers loads between nodes on the mesh.
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src        ! The source mesh with loads fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest       ! The destination mesh
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp    ! A mesh that contains the displacements associated with the source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp   ! A mesh that contains the displacements associated with the destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap    ! The mapping data structure (from Dest to Src)
   REAL(ReKi),                     INTENT(IN)     :: LoadsScaleFactor  ! Scaling factor for loads (to help with numerical issues)

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat    ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      ! local variables
!   REAL(ReKi)                                     :: RotationMatrix(3,3)
   REAL(ReKi)                                     :: torque(3), DisplacedPosition(3)
   INTEGER(IntKi)                                 :: i         ! counter over the nodes


   ErrStat = ErrID_None
   ErrMsg  = ""
   


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
FUNCTION GetLoadsScaleFactor( Src )

   TYPE(MeshType),                 INTENT(IN   )  :: Src        ! The source mesh with loads fields allocated
   REAL(ReKi)                                     :: GetLoadsScaleFactor
   
   
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
SUBROUTINE Transfer_Line2_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )

! for transfering displacement-type or force-type data from Line2 mesh to Line2 Mesh
!
   TYPE(MeshType),         INTENT(IN   ) ::  Src
   TYPE(MeshType),         INTENT(INOUT) ::  Dest
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  SrcDisp  ! a "functional" sibling of the source mesh; Src contains loads and SrcDisp contains TranslationDisp and Orientaiton
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  DestDisp ! a "functional" sibling of the destination mesh; Dest contains loads and DestDisp contains TranslationDisp and Orientaiton

   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap


   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   REAL(ReKi)                            :: LoadsScaleFactor  ! Scaling factor for loads (to help with numerical issues)
   INTEGER(IntKi)                        :: ErrStat2
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
      
   
   ErrStat = ErrID_None
   ErrMsg  = ''

   !.................
   ! Check to ensure that both the source and destination meshes are composed of Line2 elements
   !.................   
   if (Src%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Source mesh must have one or more Line2 elements.', ErrStat, ErrMsg, 'Transfer_Line2_to_Line2')
      RETURN
   endif

   if (Dest%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Destination mesh must have one or more Line2 elements.', ErrStat, ErrMsg, 'Transfer_Line2_to_Line2')
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
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Line2')
            IF (ErrStat >= AbortErrLev) RETURN

      endif !remapping

      !........................
      ! Start: Transfer data
      !........................
         
      CALL Transfer_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Line2')
         IF (ErrStat >= AbortErrLev) RETURN

   endif !algorithm for motions/scalars


   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if ( HasLoadFields(Src) ) then

      IF (.not. PRESENT(SrcDisp) .OR. .NOT. PRESENT(DestDisp) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'SrcDisp and DestDisp arguments are required for load transfer.', ErrStat, ErrMsg, 'Transfer_Line2_to_Line2')
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
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Line2')
            IF (ErrStat >= AbortErrLev) RETURN         
            
      end if

      !........................
      ! Transfer data
      !........................

      LoadsScaleFactor = GetLoadsScaleFactor ( Src )
      
      ! first, we take the source fields and transfer them to fields on the augmented source mesh:
      !  (we're also taking the SrcDisp field and putting it on our augmented mesh)
      CALL Transfer_Src_To_Augmented_Ln2_Src( Src, MeshMap, ErrStat2, ErrMsg2, SrcDisp, LoadsScaleFactor ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Line2')
         IF (ErrStat >= AbortErrLev) RETURN
      
      ! then we lump the loads from the augmented source mesh:
      CALL Lump_Line2_to_Point( MeshMap%Augmented_Ln2_Src,  MeshMap%Lumped_Points_Src,  ErrStat2, ErrMsg2, LoadsScaleFactor=LoadsScaleFactor  ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Line2')
         IF (ErrStat >= AbortErrLev) RETURN
      
      CALL Transfer_Loads_Point_to_Line2( MeshMap%Lumped_Points_Src, Dest, MeshMap, ErrStat2, ErrMsg2, MeshMap%Augmented_Ln2_Src, DestDisp, LoadsScaleFactor )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Line2_to_Line2')
         IF (ErrStat >= AbortErrLev) RETURN


   end if ! algorithm for loads


END SUBROUTINE Transfer_Line2_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_Loads_Point_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp, LoadsScaleFactor )
! Given a mapping, this routine transfers the loads from nodes on a point-element mesh to nodes on another Line2 mesh.
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src       ! The source (Line2) mesh with loads fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest      ! The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   ! The mapping data

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp   ! The source mesh's cooresponding position mesh
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp  ! The destination mesh's cooresponding position mesh
   REAL(ReKi),                     INTENT(IN)     :: LoadsScaleFactor  ! Scaling factor for loads (to help with numerical issues)

      ! local variables
   REAL(ReKi)                                     :: torque(3), DisplacedPosition(3)

   INTEGER(IntKi)                                 :: i,j       ! loop counters
   INTEGER(IntKi)                                 :: jElem     ! element number
   INTEGER(IntKi)                                 :: jNode     ! node number

   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""

      
   ! Start with transferring the loads like point-to-point (except split between two nodes of the dest element).
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
                              
               torque = Src%Force(:,i) / LoadsScaleFactor !not torque yet, but we're doing this cross product in two step to avoid tempoary memory storage
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

   ! now we can convert the lumped loads on the dest (line2) mesh into distributed loads on the same mesh:   
   CALL Convert_Point_To_Line2_Loads(Dest, MeshMap, ErrStat2, ErrMsg2, DestDisp)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Transfer_Loads_Point_to_Line2')
      !IF (ErrStat >= AbortErrLev) RETURN
   
   if (Dest%FieldMask(MASKID_FORCE) ) Dest%Force  = Dest%Force  * LoadsScaleFactor     ! whole array initialization
   if (Dest%FieldMask(MASKID_MOMENT)) Dest%Moment = Dest%Moment * LoadsScaleFactor     ! whole array initialization
   
      
END SUBROUTINE Transfer_Loads_Point_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Convert_Point_To_Line2_Loads(Dest, MeshMap, ErrStat, ErrMsg, DestDisp)
! This routine takes the lumped loads on nodes of a (line2) mesh and converts them to loads distributed across the line2 elements
! 
! Note that the matrix B is really just taking the cross product of terms.
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest                            ! The mesh (on input, the nodal loads values are lumped; on output they are distributed along the element)
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp                        ! The mesh that contains the translation displacement for these loads
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

   INTEGER(IntKi) :: jElem, n, i,j, n1, n2
   REAL(ReKi)     :: a_vec(3), sum_f(3), crossProd(3)
   REAL(ReKi)     :: c
   
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2
   
   
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
      
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Convert_Point_To_Line2_Loads')
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

      
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Convert_Point_To_Line2_Loads')
      IF (ErrStat >= AbortErrLev) RETURN        
      
   END IF ! Moment
   

END SUBROUTINE Convert_Point_To_Line2_Loads    
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
   INTEGER(IntKi)                                 :: n, n1, n2
   REAL(ReKi)                                     :: p_ED(3), p_ES(3), n1S_nD_vector(3), position(3)
   REAL(ReKi)                                     :: TmpVec(3), RefOrientation(3,3), FieldValue(3,2)   ! values for interpolating direction cosine matrices
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

   ! loop through the destination elements:
   DO jElem = 1,dest%ElemTable(Dest_TYPE)%nelem
      
      IF ( Dest_TYPE == ELEMENT_LINE2 ) THEN
         p_eD =   dest%Position(:, dest%ElemTable(Dest_TYPE)%Elements(jElem)%ElemNodes(2)) &
                - dest%Position(:, dest%ElemTable(Dest_TYPE)%Elements(jElem)%ElemNodes(1))
      ELSE
         p_eD =   dest%Position(:, dest%ElemTable(Dest_TYPE)%Elements(jElem)%ElemNodes(1))
      END IF
      
     
      iElem = 1         
      curr_Aug_NElem = Aug_NElem
      j = 1   
      
      Src_Elements: DO WHILE ( iElem <= curr_Aug_NElem )  ! bjj: we're adding elements and nodes to Temp_Ln2_Src, so this must be flexible
         
         p_eS =   Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2)) &
                - Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1))
         
         denom = DOT_PRODUCT( p_eD , p_eS )
            
  
         IF ( .NOT. EqualRealNos( denom, 0.0_ReKi) ) THEN ! we ignore source elements that are parallel to the destination element (i.e., denom == 0)
            DO jNode = j, NumNodes( Dest_TYPE ) 
               n1S_nD_vector =            dest%Position(:, dest%ElemTable(Dest_TYPE)%Elements(jElem)%ElemNodes(jNode)) &
                                - Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1))
               
               elem_position = DOT_PRODUCT( p_eD, n1S_nD_vector ) / denom
               
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
                  denom = DOT_PRODUCT( p_eD , p_eS )   ! we don't need to check that this is zero because it's just a shorter version of the temp Temp_Ln2_Src element
                  n1S_nD_vector =   dest%Position(:, dest%ElemTable(Dest_TYPE)%Elements(jElem)%ElemNodes(jNode)) &
                                   - Src%Position(:, n1 )
                  shape_fn2(Aug_Nnodes) = DOT_PRODUCT( p_eD, n1S_nD_vector ) / denom       ! save this for later, when we need to map the mesh fields...
                  
                     ! interpolate position on the original souce element:                     
                     
                  position = (1.0_ReKi - shape_fn2(Aug_Nnodes)) * Src%Position(:, n1) &
                                       + shape_fn2(Aug_Nnodes)  * Src%Position(:, n2) 
                  
                  ! let's just verify that this new node (n1) doesn't give us zero-length elements:
                  ! (note we use the NEW (not original) source element, which may have been split)
                  p_eS = position - Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(1) )                                    
                  L = SQRT(dot_product(p_eS,p_eS)) ! length of new element
                                    
                  IF ( L < TOL ) THEN ! this element is basically zero length
                        ! for numerical reasons, we really didn't want this node....
                        Aug_Nnodes = Aug_Nnodes - 1
                  ELSE
                     
                     ! let's verify the other node (n2) of this element doesn't give zero-length:
                     p_eS = position - Temp_Ln2_Src%Position(:, Temp_Ln2_Src%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(2))                                    
                     L = SQRT(dot_product(p_eS,p_eS)) ! length of new element
                     
                     IF ( L < TOL ) THEN ! this element is basically zero length
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
#ifdef __NN_ORIENTATIONS                        
                        ! set the RefOrientation based on proximity to original element's nodes:
                        IF ( NINT( shape_fn2(Aug_Nnodes) ) .EQ. 0 ) THEN
                           n = n1
                        ELSE
                           n = n2                        
                        END IF
                        RefOrientation = Src%RefOrientation(:, :, n) 
#else
                        
                           ! convert DCMs to tensors: 
                        CALL DCM_logmap( Src%RefOrientation(:, :, n1), FieldValue(:,1), ErrStat, ErrMsg )
                        IF (ErrStat >= AbortErrLev) RETURN
                  
                        CALL DCM_logmap( Src%RefOrientation(:, :, n2), FieldValue(:,2), ErrStat, ErrMsg )                  
                        IF (ErrStat >= AbortErrLev) RETURN
         
                        CALL DCM_SetLogMapForInterp( FieldValue )  ! make sure we don't cross a 2pi boundary
                  
                           ! interpolate tensors: 
                        TmpVec = (1.0_ReKi - shape_fn2(Aug_Nnodes)) * FieldValue(:, 1) &
                                           + shape_fn2(Aug_Nnodes)  * FieldValue(:, 2) 
                              
                           ! convert back to DCM:
                        RefOrientation = DCM_exp( TmpVec )

#endif
                        
                        
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
SUBROUTINE Transfer_Src_To_Augmented_Ln2_Src( Src, MeshMap, ErrStat, ErrMsg, SrcDisp, LoadsScaleFactor )
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp                         ! The displacements associated with the source mesh
   REAL(ReKi),                     INTENT(IN)     :: LoadsScaleFactor  ! Scaling factor for loads (to help with numerical issues)

      ! local variables
   INTEGER(IntKi)                                 :: iElem, i  ! do-loop counter for nodes/elements on source   
   INTEGER(IntKi)                                 :: n1, n2
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !.......................................
   ! TranslationDisp (from external mesh to augmented one)
   !.......................................
   
   IF ( SrcDisp%FIELDMASK(MASKID_TranslationDisp) ) THEN      
   
      DO i = 1,SrcDisp%NNodes
         MeshMap%Augmented_Ln2_Src%TranslationDisp(:,i) = SrcDisp%TranslationDisp(:,i)
      END DO
      
      DO i = (SrcDisp%NNodes+1),MeshMap%Augmented_Ln2_Src%NNodes
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
SUBROUTINE Create_InverseLumping_Matrix( Dest, MeshMap, ErrStat, ErrMsg )
!
!..................................................................................................................................
   TYPE(MeshType),         INTENT(IN   ) ::  Dest
   TYPE(MeshMapType),      INTENT(INOUT) ::  MeshMap


   INTEGER(IntKi),         INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(ReKi)                              :: c, TwoC

   INTEGER(IntKi)                          :: N, n1, n2  ! node numbers
   INTEGER(IntKi)                          :: j,k, iElem, iComp
   INTEGER(IntKi)                          :: ErrStat2
   CHARACTER(ErrMsgLen)                    :: ErrMsg2


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

      c    = Dest%ElemTable(ELEMENT_LINE2)%Elements(iElem)%det_jac / 3.0_ReKi
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
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Create_InverseLumping_Matrix')
            
         CALL AllocAry( MeshMap%LoadLn2_A_Mat_piv, n,     'MeshMap%LoadLn2_A_Mat_piv',ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Create_InverseLumping_Matrix')
            
         CALL AllocAry( MeshMap%LoadLn2_F,         n, 1,  'MeshMap%LoadLn2_F',        ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Create_InverseLumping_Matrix')
            
         CALL AllocAry( MeshMap%LoadLn2_M,         n, 1,  'MeshMap%LoadLn2_M',        ErrStat2, ErrMsg2)
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Create_InverseLumping_Matrix')

   END SUBROUTINE AllocInvLumpingArrays
   
END SUBROUTINE Create_InverseLumping_Matrix
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CreateLoadMap_L2_to_L2( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            ! The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

      ! Local variables:
   INTEGER(IntKi)                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                           :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      !........................
      ! Create matrix used to "unlump" loads later (needs to have element connectivity information to create it)
      !........................
      CALL Create_InverseLumping_Matrix( Dest, MeshMap, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateLoadMap_L2_to_L2')
         IF (ErrStat >= AbortErrLev) RETURN
            
      !........................
      ! An augmented Line2-element source mesh is formed
      ! by splitting the original Line2-element source mesh at each location where a
      ! destination-mesh Line2-element node projects orthogonally from the destination mesh
      !........................
         
      CALL Create_Augmented_Ln2_Src_Mesh(Src, Dest, MeshMap, ELEMENT_LINE2, ErrStat2, ErrMsg2) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateLoadMap_L2_to_L2')
         IF (ErrStat >= AbortErrLev) RETURN
                     
      !........................
      ! For each Line2-element node of the augmented source mesh, a nearest-neighbor of the augmented source mesh,
      ! a nearest-neighbor line2 element of the destination mesh is found (else aborted) in the reference configuration, 
      ! for which the source Line2-element node projects orthogonally onto the destination Line2-element domain.
      ! A destination-mesh Line2 element may be associated with multiple source-mesh Line2-element nodes.
      !........................
            
      CALL CreateMapping_ProjectToLine2(MeshMap%Augmented_Ln2_Src, Dest, MeshMap%MapLoads, ELEMENT_LINE2, ErrStat2, ErrMsg2)            
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateLoadMap_L2_to_L2')
         IF (ErrStat >= AbortErrLev) RETURN
         
         
      !........................
      ! Create a temporary mesh for lumped point elements of the line2 source 
      !........................
      CALL Create_PointMesh( MeshMap%Augmented_Ln2_Src, MeshMap%Lumped_Points_Src, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CreateLoadMap_L2_to_L2')
         IF (ErrStat >= AbortErrLev) RETURN
         
         
END SUBROUTINE CreateLoadMap_L2_to_L2        
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CreateMotionMap_L2_to_L2( Src, Dest, MeshMap, ErrStat, ErrMsg )

   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                            ! The destination mesh
   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                         ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

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
SUBROUTINE Lump_Line2_to_Point( Line2_Src, Point_Dest, ErrStat, ErrMsg, SrcDisp, LoadsScaleFactor )
!
! This subroutine takes a Line2 mesh (SRC) with distributed forces and/or moments and lumps those to the concentrated forces
! and/or moments at the endpoints (nodes) of the Line2 mesh.   These are placed in a Point mesh (DEST) where the Points are
! assumed to be colocated with the nodes of the Line2 Mesh.
!
! This routine is expected to be used in the transfer of HydroDyn Output (distributed force/moment) to the SubDyn Input
! (concentrated force/moment); routine is to be called from the Input_Ouput routine for HydroDyn to SubDyn.
!
! Formulation:
! Distrbuted quantities are integrated over the element with two-point nodal (Gauss-Lobatto) quadrature.
!
! F(x_j) = int_{-1}^{1} f(xi) phi_j(xi) J dxi (exact)
!
! F(x_1) =~ f(-1) phi_1(-1) J w_1 + f(+1) phi_1(+1) J w_2 (two-point quadrature)
! F(x_2) =~ f(-1) phi_2(-1) J w_1 + f(+1) phi_2(+1) J w_2 (two-point quadrature)
!
! which can be simplified to
!
! F(x_1) =~ f(-1) J
! F(x_2) =~ f(+1) J
!
! since w_1 = w_2 = 1, phi_j(xi_i) = Delta_{j,i}  (Kronecker Delta)
!
! where x_j is the jth node in 3D physical coordinates, xi_i is the ith node in 1D element coordinates, j in {1,2},
! xi in [-1,1] is the element natural coordinate, f(xi) is the distributed quantity, phi_j(xi) is the basis function
! for the jth node, J is the determinant of the Jacobian
! in the mapping from [X_1,X_2] to [-1,1] (X_j is the 3D location of the jth point)
!
! For Line2 elements, J = DIST(X_1,X_2) / 2
!
!..................................................................................................................................


   TYPE(MeshType),         INTENT(IN   ) ::  Line2_Src
   TYPE(MeshType),         INTENT(INOUT) ::  Point_Dest
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  SrcDisp     ! Another mesh that may contain the source's TranslationDisp field

   REAL(ReKi),             INTENT(IN)     :: LoadsScaleFactor  ! Scaling factor for loads (to help with numerical issues)

   INTEGER(IntKi),         INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

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


   ErrStat = ErrID_None
   ErrMsg = ""

#ifdef MESH_DEBUG     
   ! bjj: we shouldn't have to check this in production mode; this routine is an internal one only.
   
   if (Line2_Src%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Source mesh must have one or more Line2 Elements.', ErrStat, ErrMsg, 'Lump_Line2_to_Point')
      RETURN
   endif

   if (Point_Dest%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      CALL SetErrStat( ErrID_Fatal, 'Destination mesh must have one or more Point Elements.', ErrStat, ErrMsg, 'Lump_Line2_to_Point')
      RETURN
   endif

   if (Line2_Src%nnodes .ne. Point_Dest%nnodes) then
      CALL SetErrStat( ErrID_Fatal, 'Source and Destination meshes must have same number of Nodes.', ErrStat, ErrMsg, 'Lump_Line2_to_Point')
      RETURN
   endif

   if (Point_Dest%FieldMask(MASKID_FORCE) ) then
      if (.not. Line2_Src%FieldMask(MASKID_FORCE) ) then
         CALL SetErrStat( ErrID_Fatal, 'Destination mesh contains Force, but Source does not.', ErrStat, ErrMsg, 'Lump_Line2_to_Point')
         RETURN
      endif
   endif
!bjj: this may not be covered... 
   if (Point_Dest%FieldMask(MASKID_MOMENT) ) then
      if (.not. Line2_SRC%FieldMask(MASKID_MOMENT) ) then
         CALL SetErrStat( ErrID_Fatal, 'Destination mesh contains Moment, but Source does not.', ErrStat, ErrMsg, 'Lump_Line2_to_Point')
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
!==================================================================================================================================
!bjj: these routines require the use of ModMesh.f90, thus they cannot be part of NWTC_Library_Types.f90:
!STARTOFREGISTRYGENERATEDFILE './NWTC_Library_Types.f90'
!
! WARNING This file is generated automatically by the FAST registry
! Do not edit.  Your changes to this file will be lost.
!
! FAST Registry (v2.06.00, 14-Apr-2015)
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
   CHARACTER(1024)                :: ErrMsg2
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
  CHARACTER(1024)                :: ErrMsg2
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
      Re_BufSz   = Re_BufSz   + 1  ! distance
      Re_BufSz   = Re_BufSz   + SIZE(InData%couple_arm)  ! couple_arm
      Re_BufSz   = Re_BufSz   + SIZE(InData%shape_fn)  ! shape_fn
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

       IntKiBuf ( Int_Xferred:Int_Xferred+(1)-1 ) = InData%OtherMesh_Element
      Int_Xferred   = Int_Xferred   + 1
       ReKiBuf ( Re_Xferred:Re_Xferred+(1)-1 ) = InData%distance
      Re_Xferred   = Re_Xferred   + 1
       ReKiBuf ( Re_Xferred:Re_Xferred+(SIZE(InData%couple_arm))-1 ) = PACK(InData%couple_arm,.TRUE.)
      Re_Xferred   = Re_Xferred   + SIZE(InData%couple_arm)
       ReKiBuf ( Re_Xferred:Re_Xferred+(SIZE(InData%shape_fn))-1 ) = PACK(InData%shape_fn,.TRUE.)
      Re_Xferred   = Re_Xferred   + SIZE(InData%shape_fn)
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
  LOGICAL                        :: mask0
  LOGICAL, ALLOCATABLE           :: mask1(:)
  LOGICAL, ALLOCATABLE           :: mask2(:,:)
  LOGICAL, ALLOCATABLE           :: mask3(:,:,:)
  LOGICAL, ALLOCATABLE           :: mask4(:,:,:,:)
  LOGICAL, ALLOCATABLE           :: mask5(:,:,:,:,:)
  INTEGER(IntKi)                 :: i1, i1_l, i1_u  !  bounds (upper/lower) for an array dimension 1
  INTEGER(IntKi)                 :: i2, i2_l, i2_u  !  bounds (upper/lower) for an array dimension 2
  INTEGER(IntKi)                 :: i3, i3_l, i3_u  !  bounds (upper/lower) for an array dimension 3
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(1024)                :: ErrMsg2
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
      OutData%OtherMesh_Element = IntKiBuf( Int_Xferred ) 
      Int_Xferred   = Int_Xferred + 1
      OutData%distance = ReKiBuf( Re_Xferred )
      Re_Xferred   = Re_Xferred + 1
    i1_l = LBOUND(OutData%couple_arm,1)
    i1_u = UBOUND(OutData%couple_arm,1)
    ALLOCATE(mask1(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating mask1.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    mask1 = .TRUE. 
       OutData%couple_arm = UNPACK(ReKiBuf( Re_Xferred:Re_Xferred+(SIZE(OutData%couple_arm))-1 ), mask1, 0.0_ReKi )
      Re_Xferred   = Re_Xferred   + SIZE(OutData%couple_arm)
    DEALLOCATE(mask1)
    i1_l = LBOUND(OutData%shape_fn,1)
    i1_u = UBOUND(OutData%shape_fn,1)
    ALLOCATE(mask1(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating mask1.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    mask1 = .TRUE. 
       OutData%shape_fn = UNPACK(ReKiBuf( Re_Xferred:Re_Xferred+(SIZE(OutData%shape_fn))-1 ), mask1, 0.0_ReKi )
      Re_Xferred   = Re_Xferred   + SIZE(OutData%shape_fn)
    DEALLOCATE(mask1)
 END SUBROUTINE NWTC_Library_UnPackMapType

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
   CHARACTER(1024)                :: ErrMsg2
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
IF (ALLOCATED(MeshMapTypeData%LoadLn2_F)) THEN
   DEALLOCATE(MeshMapTypeData%LoadLn2_F)
ENDIF
IF (ALLOCATED(MeshMapTypeData%LoadLn2_A_Mat)) THEN
   DEALLOCATE(MeshMapTypeData%LoadLn2_A_Mat)
ENDIF
IF (ALLOCATED(MeshMapTypeData%LoadLn2_M)) THEN
   DEALLOCATE(MeshMapTypeData%LoadLn2_M)
ENDIF
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
  CHARACTER(1024)                :: ErrMsg2
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
      Re_BufSz   = Re_BufSz   + SIZE(InData%DisplacedPosition)  ! DisplacedPosition
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! LoadLn2_F allocated yes/no
  IF ( ALLOCATED(InData%LoadLn2_F) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! LoadLn2_F upper/lower bounds for each dimension
      Re_BufSz   = Re_BufSz   + SIZE(InData%LoadLn2_F)  ! LoadLn2_F
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! LoadLn2_A_Mat allocated yes/no
  IF ( ALLOCATED(InData%LoadLn2_A_Mat) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! LoadLn2_A_Mat upper/lower bounds for each dimension
      Re_BufSz   = Re_BufSz   + SIZE(InData%LoadLn2_A_Mat)  ! LoadLn2_A_Mat
  END IF
  Int_BufSz   = Int_BufSz   + 1     ! LoadLn2_M allocated yes/no
  IF ( ALLOCATED(InData%LoadLn2_M) ) THEN
    Int_BufSz   = Int_BufSz   + 2*2  ! LoadLn2_M upper/lower bounds for each dimension
      Re_BufSz   = Re_BufSz   + SIZE(InData%LoadLn2_M)  ! LoadLn2_M
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

      IF (SIZE(InData%LoadLn2_A_Mat_Piv)>0) IntKiBuf ( Int_Xferred:Int_Xferred+(SIZE(InData%LoadLn2_A_Mat_Piv))-1 ) = PACK(InData%LoadLn2_A_Mat_Piv,.TRUE.)
      Int_Xferred   = Int_Xferred   + SIZE(InData%LoadLn2_A_Mat_Piv)
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

      IF (SIZE(InData%DisplacedPosition)>0) ReKiBuf ( Re_Xferred:Re_Xferred+(SIZE(InData%DisplacedPosition))-1 ) = PACK(InData%DisplacedPosition,.TRUE.)
      Re_Xferred   = Re_Xferred   + SIZE(InData%DisplacedPosition)
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

      IF (SIZE(InData%LoadLn2_F)>0) ReKiBuf ( Re_Xferred:Re_Xferred+(SIZE(InData%LoadLn2_F))-1 ) = PACK(InData%LoadLn2_F,.TRUE.)
      Re_Xferred   = Re_Xferred   + SIZE(InData%LoadLn2_F)
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

      IF (SIZE(InData%LoadLn2_A_Mat)>0) ReKiBuf ( Re_Xferred:Re_Xferred+(SIZE(InData%LoadLn2_A_Mat))-1 ) = PACK(InData%LoadLn2_A_Mat,.TRUE.)
      Re_Xferred   = Re_Xferred   + SIZE(InData%LoadLn2_A_Mat)
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

      IF (SIZE(InData%LoadLn2_M)>0) ReKiBuf ( Re_Xferred:Re_Xferred+(SIZE(InData%LoadLn2_M))-1 ) = PACK(InData%LoadLn2_M,.TRUE.)
      Re_Xferred   = Re_Xferred   + SIZE(InData%LoadLn2_M)
  END IF
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
  LOGICAL                        :: mask0
  LOGICAL, ALLOCATABLE           :: mask1(:)
  LOGICAL, ALLOCATABLE           :: mask2(:,:)
  LOGICAL, ALLOCATABLE           :: mask3(:,:,:)
  LOGICAL, ALLOCATABLE           :: mask4(:,:,:,:)
  LOGICAL, ALLOCATABLE           :: mask5(:,:,:,:,:)
  INTEGER(IntKi)                 :: i1, i1_l, i1_u  !  bounds (upper/lower) for an array dimension 1
  INTEGER(IntKi)                 :: i2, i2_l, i2_u  !  bounds (upper/lower) for an array dimension 2
  INTEGER(IntKi)                 :: i3, i3_l, i3_u  !  bounds (upper/lower) for an array dimension 3
  INTEGER(IntKi)                 :: ErrStat2
  CHARACTER(1024)                :: ErrMsg2
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
    ALLOCATE(mask1(i1_l:i1_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating mask1.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    mask1 = .TRUE. 
      IF (SIZE(OutData%LoadLn2_A_Mat_Piv)>0) OutData%LoadLn2_A_Mat_Piv = UNPACK( IntKiBuf ( Int_Xferred:Int_Xferred+(SIZE(OutData%LoadLn2_A_Mat_Piv))-1 ), mask1, 0_IntKi )
      Int_Xferred   = Int_Xferred   + SIZE(OutData%LoadLn2_A_Mat_Piv)
    DEALLOCATE(mask1)
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
    ALLOCATE(mask3(i1_l:i1_u,i2_l:i2_u,i3_l:i3_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating mask3.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    mask3 = .TRUE. 
      IF (SIZE(OutData%DisplacedPosition)>0) OutData%DisplacedPosition = UNPACK(ReKiBuf( Re_Xferred:Re_Xferred+(SIZE(OutData%DisplacedPosition))-1 ), mask3, 0.0_ReKi )
      Re_Xferred   = Re_Xferred   + SIZE(OutData%DisplacedPosition)
    DEALLOCATE(mask3)
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
    ALLOCATE(mask2(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating mask2.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    mask2 = .TRUE. 
      IF (SIZE(OutData%LoadLn2_F)>0) OutData%LoadLn2_F = UNPACK(ReKiBuf( Re_Xferred:Re_Xferred+(SIZE(OutData%LoadLn2_F))-1 ), mask2, 0.0_ReKi )
      Re_Xferred   = Re_Xferred   + SIZE(OutData%LoadLn2_F)
    DEALLOCATE(mask2)
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
    ALLOCATE(mask2(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating mask2.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    mask2 = .TRUE. 
      IF (SIZE(OutData%LoadLn2_A_Mat)>0) OutData%LoadLn2_A_Mat = UNPACK(ReKiBuf( Re_Xferred:Re_Xferred+(SIZE(OutData%LoadLn2_A_Mat))-1 ), mask2, 0.0_ReKi )
      Re_Xferred   = Re_Xferred   + SIZE(OutData%LoadLn2_A_Mat)
    DEALLOCATE(mask2)
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
    ALLOCATE(mask2(i1_l:i1_u,i2_l:i2_u),STAT=ErrStat2)
    IF (ErrStat2 /= 0) THEN 
       CALL SetErrStat(ErrID_Fatal, 'Error allocating mask2.', ErrStat, ErrMsg,RoutineName)
       RETURN
    END IF
    mask2 = .TRUE. 
      IF (SIZE(OutData%LoadLn2_M)>0) OutData%LoadLn2_M = UNPACK(ReKiBuf( Re_Xferred:Re_Xferred+(SIZE(OutData%LoadLn2_M))-1 ), mask2, 0.0_ReKi )
      Re_Xferred   = Re_Xferred   + SIZE(OutData%LoadLn2_M)
    DEALLOCATE(mask2)
  END IF
 END SUBROUTINE NWTC_Library_UnPackMeshMapType
!*********************************************************************************************************************************
!ENDOFREGISTRYGENERATEDFILE


!----------------------------------------------------------------------------------------------------------------------------------
END MODULE ModMesh_Mapping
!----------------------------------------------------------------------------------------------------------------------------------
