!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
!
!    NWTC Subroutine Library is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with NWTC Subroutine Library.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
MODULE ModMesh_Mapping

   USE ModMesh

   IMPLICIT NONE

   PRIVATE

   TYPE, PUBLIC :: MapType
      INTEGER(IntKi) :: OtherMesh_Element    ! Node (for point meshes) or Element (for line2 meshes) number on other mesh; for loads, other mesh is Dest, for motions/scalars, other mesh is Src
      REAL(ReKi)     :: distance             ! magnitude of couple_arm
      REAL(ReKi)     :: couple_arm(3)        ! Vector between a point and node 1 of an element
!      REAL(ReKi)     :: elem_position        ! 1-D element-level location [0,1] based on closest-line projection of point
      REAL(ReKi)     :: shape_fn(2)          ! shape functions: 1-D element-level location [0,1] based on closest-line projection of point
   END TYPE MapType

   TYPE, PUBLIC :: MeshMapType
      TYPE(MapType),ALLOCATABLE :: MapLoads(:)
      TYPE(MapType),ALLOCATABLE :: MapMotions(:)
      TYPE(MeshType)            :: Temp_Lumped_Points_Src  ! Stored here for efficiency
      TYPE(MeshType)            :: Temp_Lumped_Points_Dest ! Stored here for efficiency
      REAL(ReKi),ALLOCATABLE    :: RotatedPosition(:,:)    ! Orientation^T * RefOrientation * couple_arm for each mapped node (NOTE that changes with the orientation field. stored here for efficiency.)
      REAL(ReKi),ALLOCATABLE    :: RotatedPosition_n2(:,:) ! Orientation^T * RefOrientation * distance vector for each mapped node (NOTE that changes with the orientation field. stored here for efficiency; needed for linear combinations of 2 nodes)
   END TYPE

      ! note that these parameters must be negative (positive indicates the node/element it is mapped to)
   INTEGER(IntKi),  PARAMETER   :: NODE_NOT_MAPPED = -1
   INTEGER(IntKi),  PARAMETER   :: NODE_IGNORED    = -2

   PUBLIC :: AllocMapping
   PUBLIC :: MeshMapDestroy
   PUBLIC :: Transfer_Point_to_Point
   PUBLIC :: Transfer_Line2_to_Point
   PUBLIC :: Transfer_Point_to_Line2
   PUBLIC :: Transfer_Line2_to_Line2
   !PUBLIC :: Lump_Line2_to_Point



CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AllocMapping( Src, Dest, MeshMap, ErrStat, ErrMsg )
! This subroutine takes two meshes, determines the sizes required for the mapping data structure, and then
! allocates the mappings (different for loads and motions/scalars).
!..................................................................................................................................

   TYPE(MeshType),           INTENT(IN)     ::  Src
   TYPE(MeshType),           INTENT(IN)     ::  Dest

   TYPE(MeshMapType),        INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),           INTENT(  OUT)  :: ErrStat      ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT)  :: ErrMsg       ! Error message if ErrStat /= ErrID_None


      ! local variables:

   INTEGER(IntKi)                           :: PointsInMap
   LOGICAL                                  :: MapCreated

   ErrStat = ErrID_None
   ErrMsg  = ''

   MapCreated = .FALSE.


   IF ( .NOT. Dest%Committed .OR. .NOT. Src%Committed ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = " Both meshes must be committed before they can be mapped."
      RETURN
   END IF


      !................................................
      ! Allocate the mapping for Motions and Scalars (if both meshes have some):
      !................................................
   IF ( HasMotionFields(Src) .AND. HasMotionFields(Dest) ) THEN

      PointsInMap = Dest%Nnodes

      IF ( PointsInMap < 1 ) THEN
         ErrStat = ErrID_Fatal
         IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg  = TRIM(ErrMsg)//' AllocMapping: MeshMap%MapMotions not allocated because no nodes were found to map.'
      ELSE

            ! Allocate the mapping structure:
         ALLOCATE( MeshMap%MapMotions(PointsInMap), STAT=ErrStat )
         IF ( ErrStat /= 0 ) THEN
            ErrStat = ErrID_Fatal
            IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
            ErrMsg  = ' AllocMapping: Error trying to allocate MeshMap%MapMotions.'
         ELSE
            MapCreated = .TRUE.
         END IF

      END IF

   END IF !HasMotionFields


      !................................................
      ! Allocate the mapping for Loads:
      !................................................
   IF ( HasLoadFields(Src) .AND. HasLoadFields(Dest) ) THEN


      ! check that the appropriate combinations of source/destination force/moments exist:
      IF ( Src%FieldMask(MASKID_Force) ) THEN
         IF (.NOT. Dest%FieldMask(MASKID_Force) ) THEN
            ErrStat = ErrID_Fatal
            IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
            ErrMsg  = ' AllocMapping: Destination mesh does not contain force but source mesh does.'
         END IF
      END IF
      IF ( Src%FieldMask(MASKID_Moment) ) THEN
         IF (.NOT. Dest%FieldMask(MASKID_Moment) ) THEN
            ErrStat = ErrID_Fatal
            IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
            ErrMsg  = ' AllocMapping: Destination mesh does not contain moment but source mesh moment.'
         END IF
      END IF

      ! get size of mapping:
      IF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 .AND. Dest%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN
         PointsInMap = Src%ElemTable(ELEMENT_LINE2)%nelem  !we're mapping elements to element
      ELSE
         PointsInMap = Src%Nnodes !otherwise, we're mapping node to element
      END IF


      IF ( PointsInMap < 1 ) THEN
         ErrStat = ErrID_Fatal
         IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg  = ' AllocMapping: MeshMap%MapLoads not allocated because no nodes were found to map.'
      ELSE

            ! Allocate the mapping structure:
         ALLOCATE( MeshMap%MapLoads(PointsInMap), STAT=ErrStat )
         IF ( ErrStat /= 0 ) THEN
            ErrStat = ErrID_Fatal
            IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
            ErrMsg  = ' AllocMapping: Error trying to allocate MeshMap%MapLoads.'
         ELSE
            MapCreated = .TRUE.
         END IF

      END IF

   END IF ! HasLoadFields

   IF ( .NOT. MapCreated ) THEN
      ErrStat = ErrID_Fatal
      IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg  = ' AllocMapping: Neither MapMotions or MapLoads was allocated. Meshes may not have compatible fields for mapping.'
   END IF


      !................................................
      ! Allocate the RotatedPosition field:
      !................................................

   IF ( Src%ElemTable(ELEMENT_LINE2)%nelem > 0 ) THEN !Line2-to-Point and Line2-to-Line2
      !PointsInMap = MAX(Dest%nelemlist,Src%nelemlist)
      PointsInMap = MAX(Dest%Nnodes,Src%Nnodes)

      ALLOCATE( MeshMap%RotatedPosition(3, PointsInMap), STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg  = ' AllocMapping: Error trying to allocate MeshMap%RotatedPosition.'
      END IF

      ALLOCATE( MeshMap%RotatedPosition_n2(3, PointsInMap), STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg  = ' AllocMapping: Error trying to allocate MeshMap%RotatedPosition_n2.'
      END IF
   ELSE
      PointsInMap = MAX(Dest%Nnodes,Src%Nnodes)

      ALLOCATE( MeshMap%RotatedPosition(3, PointsInMap), STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrStat = ErrID_Fatal
         IF ( LEN_TRIM(ErrMsg) /= 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg  = ' AllocMapping: Error trying to allocate MeshMap%RotatedPosition.'
      END IF

   END IF


END SUBROUTINE AllocMapping
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE MeshMapDestroy( MeshMap, ErrStat, ErrMsg )
! This routine destroys the elements of the MeshMapType
!..................................................................................................................................

   TYPE(MeshMapType),        INTENT(INOUT)  :: MeshMap

   INTEGER(IntKi),           INTENT(  OUT)  :: ErrStat      ! Error status of the operation
   CHARACTER(*),             INTENT(  OUT)  :: ErrMsg       ! Error message if ErrStat /= ErrID_None


   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (ALLOCATED(MeshMap%MapLoads          )) DEALLOCATE(MeshMap%MapLoads          )
   IF (ALLOCATED(MeshMap%MapMotions        )) DEALLOCATE(MeshMap%MapMotions        )
   IF (ALLOCATED(MeshMap%RotatedPosition   )) DEALLOCATE(MeshMap%RotatedPosition   )
   IF (ALLOCATED(MeshMap%RotatedPosition_n2)) DEALLOCATE(MeshMap%RotatedPosition_n2)

   CALL MeshDestroy( MeshMap%Temp_Lumped_Points_Src,  ErrStat, ErrMsg )
   CALL MeshDestroy( MeshMap%Temp_Lumped_Points_Dest, ErrStat, ErrMsg )

END SUBROUTINE MeshMapDestroy
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcOrientDisp )

! for transfering displacement/velocity/accerlation data from Line2 mesh to Point Mesh
!
! Also works for transfering
!
   TYPE(MeshType),         INTENT(IN   ) ::  Src
   TYPE(MeshType),         INTENT(INOUT) ::  Dest
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  SrcOrientDisp ! Mesh with same nodes as Src, and containing orientation and translational displacement information

   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None



   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

   if (Src%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Transfer_Line2_to_Point: Src mesh must have one or more Line2 Elements '
      RETURN
   endif

   if (Dest%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Transfer_Line2_to_Point: Destination mesh must have one or more Point Elements '
      RETURN
   endif


   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------

   if (     Src%FieldMask(MASKID_TranslationDisp)    &
       .or. Src%FieldMask(MASKID_Orientation)        &
       .or. Src%FieldMask(MASKID_TranslationVel)     &
       .or. Src%FieldMask(MASKID_ROTATIONVEL)        &
       .or. Src%FieldMask(MASKID_TRANSLATIONACC)     &
       .or. Src%FieldMask(MASKID_ROTATIONACC)        &
       .or. Src%FieldMask(MASKID_SCALAR)             &
       ) then


      !........................
      ! Check that the mapping data structure is the correct size:
      !........................

      if (SIZE(MeshMap%MapMotions) < Dest%nnodes) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Transfer_Line2_to_Point: MeshMap%MapMotions(:) should be allocated to Dest%nnodes.'
         RETURN
      endif


      !........................
      ! Start: Create Mapping data (if remap is true)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMapping_ProjectToLine2(Dest,Src, MeshMap%MapMotions, ELEMENT_POINT, ErrStat, ErrMsg)
            IF (ErrStat >= AbortErrLev) RETURN

      endif !remapping


      !........................
      ! Start: Transfer data
      !........................

      CALL Transfer_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )
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

   if ( Src%FieldMask(maskid_force) .or. Src%FieldMask(maskid_moment) ) then

      !........................
      ! Check that the mapping data structure is the correct size:
      !........................

      if (SIZE(MeshMap%MapLoads) < Src%nnodes) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Transfer_Line2_to_Point: MeshMap%MapLoads(:) should be allocated to Src%nnodes.'
         RETURN
      endif

      !.................
      ! other checks for available mesh fields (now done in alloc mapping routine)
      !.................

      !........................
      ! Create mapping (including the temporary src mesh)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         !........................
         ! Create a temporary mesh for lumped point elements of the line2 mesh
         !........................
         CALL Create_PointMesh( Src, MeshMap%Temp_Lumped_Points_Src, ErrStat, ErrMsg )
            IF (ErrStat >= AbortErrLev) RETURN

         ! in following call, Src is mesh to loop over, finding a corresponding point for each point in Dest
         CALL CreateMapping_NearestNeighbor( MeshMap%Temp_Lumped_Points_Src, Dest, MeshMap%MapLoads, ELEMENT_POINT, ELEMENT_POINT, ErrStat, ErrMsg )
            IF (ErrStat >= AbortErrLev) RETURN

      ELSE ! Check that the temporary mesh has been set

         IF ( .NOT. MeshMap%Temp_Lumped_Points_Src%Initialized ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = ' Error in Transfer_Line2_to_Point: MeshMap%Temp_Lumped_Points_Src not initialized (set RemapFlag = .TRUE.).'
            RETURN
         END IF
      END IF

      !........................
      ! Lump the values on line2 elements to nodes on a point mesh
      ! These point and line2 meshes have the same nodes so there
      ! is no mapping here:
      !........................

      IF ( PRESENT(SrcOrientDisp) ) THEN
         CALL Lump_Line2_to_Point( Src, MeshMap%Temp_Lumped_Points_Src, ErrStat, ErrMsg, SrcOrientDisp )
      ELSE
         CALL Lump_Line2_to_Point( Src, MeshMap%Temp_Lumped_Points_Src, ErrStat, ErrMsg )
      END IF
         IF (ErrStat >= AbortErrLev) RETURN


      !........................
      ! Transfer data
      !........................

      IF ( PRESENT( SrcOrientDisp ) ) THEN
         CALL Transfer_Loads_Point_to_Point( MeshMap%Temp_Lumped_Points_Src, Dest, MeshMap%MapLoads, ErrStat, ErrMsg, SrcOrientDisp )
      ELSE
         CALL Transfer_Loads_Point_to_Point( MeshMap%Temp_Lumped_Points_Src, Dest, MeshMap%MapLoads, ErrStat, ErrMsg )
      END IF
      IF (ErrStat >= AbortErrLev) RETURN


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
   INTEGER(IntKi)  :: i                                        ! counter over the nodes
   INTEGER(IntKi)            :: n, n1, n2                      ! temporary space for node numbers
   REAL(ReKi)                :: FieldValueN1(3)                ! Temporary variable to store field values on element nodes
   REAL(ReKi)                :: FieldValueN2(3)                ! Temporary variable to store field values on element nodes


   ErrStat = ErrID_None
   ErrMsg  = ""


      ! calculate MeshMap%RotatedPosition and MeshMap%RotatedPosition_n2:
   CALL Calculate_RotatedPosition_N2( Src, Dest, MeshMap, ErrStat, ErrMsg )

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

               ! augment translation
            FieldValueN1 = FieldValueN1 + MeshMap%RotatedPosition(:,i)    - (Dest%Position(:,i)-Src%Position(:,n1))
            FieldValueN2 = FieldValueN2 + MeshMap%RotatedPosition_n2(:,i) - (Dest%Position(:,i)-Src%Position(:,n2))

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


            ! Get the closest node (node closest in the projection)
            ! We can't interpolate the DCM, so we'll use the nearest neighbor.
         ! MeshMap%MapMotions(i)%shape_fn(1)*FieldValueN1 + MeshMap%MapMotions(i)%shape_fn(2)*FieldValueN2

         IF ( NINT( MeshMap%MapMotions(i)%shape_fn(2) ) == 0 ) THEN
            n = n1
         ELSE
            n = n2
         END IF

         Dest%Orientation(:,:,i) = MATMUL( MATMUL( Dest%RefOrientation(:,:,i), TRANSPOSE( Src%RefOrientation(:,:,n) ) )&
                                          , Src%Orientation(:,:,n) )
      end do

   endif

      ! ---------------------------- TranslationVel  --------------------------------------------

   if ( Src%FieldMask(MASKID_TranslationVel) .AND. Dest%FieldMask(MASKID_TranslationVel) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         FieldValueN1 = Src%TranslationVel(:,n1)
         FieldValueN2 = Src%TranslationVel(:,n2)

         if ( Src%FieldMask(MASKID_RotationVel) ) then
            FieldValueN1 = FieldValueN1 + cross_product ( Src%RotationVel(:,n1), MeshMap%RotatedPosition(:,i) )
            FieldValueN2 = FieldValueN2 + cross_product ( Src%RotationVel(:,n2), MeshMap%RotatedPosition_n2(:,i) )
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
            FieldValueN1 = FieldValueN1 + cross_product( Src%RotationAcc(:,n1),  MeshMap%RotatedPosition(   :,i) )
            FieldValueN2 = FieldValueN2 + cross_product( Src%RotationAcc(:,n2),  MeshMap%RotatedPosition_n2(:,i) )
         endif

         if ( Src%FieldMask(MASKID_RotationVel) )  then
            FieldValueN1 =  FieldValueN1 + cross_product( Src%RotationVel(:,n1), &
                                           cross_product( Src%RotationVel(:,n1), MeshMap%RotatedPosition(   :,i) ))
            FieldValueN2 =  FieldValueN2 + cross_product( Src%RotationVel(:,n2), &
                                           cross_product( Src%RotationVel(:,n2), MeshMap%RotatedPosition_n2(:,i) ))
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
SUBROUTINE CreateMapping_ProjectToLine2(Mesh1, Mesh2, Map, Mesh1_TYPE, ErrStat, ErrMsg, UseClosestNode)
!This routine projects Mesh1 onto a Line2 mesh (Mesh2) to find the element mappings between the two meshes.
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Mesh1                           ! The mesh in the outer mapping loop (Dest for Motions/Scalars; Src for Loads)
   TYPE(MeshType),                 INTENT(IN   )  :: Mesh2                           ! The mesh in the inner mapping loop (Src for Motions/Scalars; Dest for Loads)

   TYPE(MapType),                  INTENT(INOUT)  :: Map(:)                          ! The mapping from Src to Dest

   INTEGER(IntKi),                 INTENT(IN   )  :: Mesh1_TYPE                      ! Type of Mesh1 elements to map
   INTEGER(IntKi),   PARAMETER                    :: Mesh2_TYPE  = ELEMENT_LINE2     ! Type of Mesh2 elements on map (MUST BE ELEMENT_LINE2)

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                        ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                         ! Error message if ErrStat /= ErrID_None
   LOGICAL,          OPTIONAL,     INTENT(IN   )  :: UseClosestNode                 ! Option to ignore closest node (if UseClosestNode=FALSE) if node doesn't prlject on any elements

      ! local variables

   REAL(ReKi)      :: denom
   REAL(ReKi)      :: dist
   REAL(ReKi)      :: min_dist
   REAL(ReKi)      :: elem_position

   REAL(ReKi)      :: Mesh1_xyz(3)
   REAL(ReKi)      :: Mesh2_xyz(3)

   REAL(ReKi)      :: n1_n2_vector(3)     ! vector going from node 1 to node 2 in Line2 element
   REAL(ReKi)      :: n1_Point_vector(3)  ! vector going from node 1 in Line 2 element to Destination Point
   REAL(ReKi)      :: tmp(3)              ! temporary vector for cross product calculation

   INTEGER(IntKi)  :: point_with_min_dist
   INTEGER(IntKi)  :: elem_with_min_dist

   INTEGER(IntKi)  :: iElem, iNode, i  ! do-loop counter for elements on Mesh1, associated node(S)
   INTEGER(IntKi)  :: jElem, jNode, j  ! do-loop counter for elements on Mesh2, associated node

   INTEGER(IntKi)  :: n1, n2           ! nodes associated with an element

   INTEGER(IntKi)  :: Mesh2NodeElem(Mesh2%NNodes,2) ! Gives the mesh node a corresponding element and element node (which of the nodes making up the element)
   LOGICAL         :: UseMesh2Node(Mesh2%NNodes)    ! determines if the node on the second mesh is part of the mapping (i.e., contained in an element of the appropriate type)

   LOGICAL         :: found
   LOGICAL         :: on_element



      ! initialization
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Determine which nodes on mesh2 are going to be in the mapping
   UseMesh2Node = .FALSE.
   do jElem = 1, Mesh2%ElemTable(Mesh2_TYPE)%nelem  ! number of point elements on Mesh2
      do jNode = 1, SIZE( Mesh2%ElemTable(Mesh2_TYPE)%Elements(jElem)%ElemNodes )
         UseMesh2Node( Mesh2%ElemTable(Mesh2_TYPE)%Elements(jElem)%ElemNodes(jNode)    ) = .TRUE. ! use this node:
         Mesh2NodeElem( Mesh2%ElemTable(Mesh2_TYPE)%Elements(jElem)%ElemNodes(jNode),1 ) = jElem  !   the node is part of this element
         Mesh2NodeElem( Mesh2%ElemTable(Mesh2_TYPE)%Elements(jElem)%ElemNodes(jNode),2 ) = jNode  !   and it's the jNode of this element
      end do
   end do

   ! Map the source nodes to destination nodes:
   Map(:)%OtherMesh_Element = NODE_NOT_MAPPED ! initialize this so we know if we've mapped this node already (done only because we may have different elements)


   do iElem = 1, Mesh1%ElemTable(Mesh1_TYPE)%nelem   ! number of Mesh1_TYPE elements on Mesh1
      do iNode = 1, SIZE( Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes )
         i = Mesh1%ElemTable(Mesh1_TYPE)%Elements(iElem)%ElemNodes(iNode)  ! the nodes on element iElem
         IF ( Map(i)%OtherMesh_Element > 0 ) CYCLE  ! we already mapped this node; let's move on to the next iNode (or iElem)

         ! destination point
         Mesh1_xyz = Mesh1%Position(:, i)

         found = .false.
         min_dist = HUGE(min_dist)

         do jElem = 1, Mesh2%ElemTable(ELEMENT_LINE2)%nelem  ! ELEMENT_LINE2 = Mesh2_TYPE

               ! write(*,*) 'i,jElem = ', i,jElem, 'found = ', found

               ! grab node numbers associated with the jElem_th element
            n1 = Mesh2%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(1)
            n2 = Mesh2%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(2)

               ! Calculate vectors used in projection operation

            n1_n2_vector    = Mesh2%Position(:,n2) - Mesh2%Position(:,n1)
            n1_Point_vector = Mesh1_xyz - Mesh2%Position(:,n1)

            denom           = DOT_PRODUCT( n1_n2_vector, n1_n2_vector )
            IF ( EqualRealNos( denom, 0.0_ReKi ) ) THEN
               ErrStat = ErrID_Fatal
               ErrMsg  = ' Error in CreateMapping_ProjectToLine2: Division by zero because Line2 element nodes are in same position.'
               RETURN
            END IF

               ! project point onto line defined by n1 and n2

            elem_position = DOT_PRODUCT(n1_n2_vector,n1_Point_vector) / denom

                  ! note: i forumlated it this way because Fortran doesn't do shortcutting and I don't want to call EqualRealNos if we don't need it:
            if ( elem_position .le. 1.0_ReKi .and. elem_position .ge. 0.0_ReKi ) then !we're ON the element (between the two nodes)
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

                  Map(i)%OtherMesh_Element = jElem
                  !Map(i)%elem_position     = elem_position
                  Map(i)%shape_fn(1)       = 1.0_ReKi - elem_position
                  Map(i)%shape_fn(2)       = elem_position

                  !Map(i)%couple_arm        = n1_Point_vector

               !write(*,*) 'found the element ', jElem
               !write(*,*) 'elem_position ',  elem_position
               end if !the point is closest to this line2 element

            endif

         end do !jElem

            ! if failed to find an element that the Point projected into, find the nearest point that is in the Line2 Mesh
            ! compare with CreateMapping_NearestNeighbor algorithm
         if (.not. found) then

            IF ( PRESENT(UseClosestNode) ) THEN
               IF ( .NOT. UseClosestNode ) THEN !Let's ignore this node
                  Map(i)%OtherMesh_Element = NODE_IGNORED
                  CYCLE ! go to the next iNode (or iElem)
               END IF
            END IF

            ! Find the nearest neighbor node for this particular node

            ! initialize minimum distance marker at some huge number
            min_dist = HUGE(min_dist)
            point_with_min_dist = 0
            elem_with_min_dist = 0

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
                  elem_with_min_dist = Mesh2NodeElem(j,1)         ! node j is on this element
                  elem_position      = Mesh2NodeElem(j,2) - 1     ! (if we're at the element's node 1, we want elem_position to be 0; if we're at node 2, elem_position is 1.)

                  !if (EqualRealNos(dist), 0.0_ReKi)) EXIT !we have an exact match so let's just stop looking
               endif

            end do !j

            Map(i)%OtherMesh_Element = elem_with_min_dist
            !Map(i)%elem_position     = elem_position
            Map(i)%shape_fn(1)       = 1.0_ReKi - elem_position
            Map(i)%shape_fn(2)       = elem_position

            if (elem_with_min_dist .lt. 1 )  then
               ErrStat = ErrID_Fatal
               ErrMsg  = ' Error in CreateMapping_ProjectToLine2: Failed to find destination point associated with source point.'
               RETURN
            endif

         end if !not found on projection to element

      end do !iNode
   end do !iElem

END SUBROUTINE CreateMapping_ProjectToLine2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Create_PointMesh(Src, Temp_Point_Src, ErrStat, ErrMsg)
!This routine creates new mesh with the same positions as the Src mesh, except all of the elements are points. It adds fields for
! forces, moments, TranslationDisp, and/or Orientation, if they are part of the Src mesh
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Src                             ! The source mesh
   TYPE(MeshType),                 INTENT(INOUT)  :: Temp_Point_Src                  ! A blank mesh to be

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
   INTEGER(IntKi)                                 :: i !loop over the nodes


   CALL MeshDestroy( Temp_Point_Src, ErrStat, ErrMsg, .TRUE. )


   CALL MeshCreate(   BlankMesh       = Temp_Point_Src               &
                     ,IOS             = Src%IOS                      &
                     ,NNodes          = Src%nnodes                   &
                     ,Force           = Src%FieldMask(maskid_force)  &
                     ,Moment          = Src%FieldMask(maskid_moment) &
                     ,TranslationDisp = Src%FieldMask(maskid_TranslationDisp)  &
                     ,Orientation     = Src%FieldMask(maskid_Orientation) &
                     ,ErrStat         = ErrStat                      &
                     ,ErrMess         = ErrMsg                 )

   do i = 1, src%nnodes

      CALL MeshConstructElement ( Mesh = Temp_Point_Src         &
                                 ,Xelement = ELEMENT_POINT      &
                                 ,P1       = I                  &
                                 ,ErrStat  = ErrStat            &
                                 ,ErrMess  = ErrMsg             )

      CALL MeshPositionNode ( Mesh = Temp_Point_Src               &
                              ,INode = i                          &
                              ,Pos = Src%Position(:,i)            &
                              ,Orient = Src%RefOrientation(:,:,i) &
                              ,ErrStat   = ErrStat                &
                              ,ErrMess   = ErrMsg                 )

   enddo

   CALL MeshCommit ( Temp_Point_Src, ErrStat, ErrMsg )

END SUBROUTINE Create_PointMesh
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Calculate_RotatedPosition_N2( Src, Dest, MeshMap, ErrStat, ErrMsg )
! This routine calculates the terms used in many of the transformation equations:
!
! MeshMap%RotatedPosition(:,:,i) =
! Src%Orientation(:,:,n)^T * Src%RefOrientation(:,:,n) * (Dest%Position(:,i)-Src%Position(:,n))
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src       ! The source mesh with motion fields allocated
   TYPE(MeshType),                 INTENT(IN   )  :: Dest      ! The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   ! data for the mesh mapping

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(ReKi)                                     :: RotationMatrixN1(3,3)
   REAL(ReKi)                                     :: RotationMatrixN2(3,3)
   INTEGER(IntKi)                                 :: n1, n2                  ! nodes associated with an element

   INTEGER(IntKi)                                 :: i                       ! loop counter

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF ( Src%FieldMask(MASKID_Orientation) ) THEN
      DO i=1, Dest%Nnodes
         IF ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         RotationMatrixN1 = MATMUL( TRANSPOSE( Src%Orientation(:,:,n1) ) , Src%RefOrientation(:,:,n1) )
         RotationMatrixN2 = MATMUL( TRANSPOSE( Src%Orientation(:,:,n2) ) , Src%RefOrientation(:,:,n2) )

         MeshMap%RotatedPosition(:,i)    = MATMUL( RotationMatrixN1, Dest%Position(:,i)-Src%Position(:,n1)  )
         MeshMap%RotatedPosition_n2(:,i) = MATMUL( RotationMatrixN2, Dest%Position(:,i)-Src%Position(:,n2)  )
      END DO
   ELSE

      DO i = 1,Dest%Nnodes
         IF ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(MeshMap%MapMotions(i)%OtherMesh_Element)%ElemNodes(2)

         MeshMap%RotatedPosition(:,i)    = Dest%Position(:,i)-Src%Position(:,n1)
         MeshMap%RotatedPosition_n2(:,i) = Dest%Position(:,i)-Src%Position(:,n2)
      END DO

   END IF


END SUBROUTINE Calculate_RotatedPosition_N2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_Point_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg )

! for transfering displacement/velocity/accerlation data from Point Mesh to Line2 Mesh

! Example: SubDyn  Output: Point Mesh: Displacement/Orientation/TranslationVel, etc
!          HyDroDyn Input: Line2 Mesh: "

   TYPE(MeshType),         INTENT(IN   ) ::  Src
   TYPE(MeshType),         INTENT(INOUT) ::  Dest

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   TYPE(MeshMapType),      INTENT(INOUT) :: MeshMap

   ! local variables

   INTEGER(IntKi)  ::   i ! do-loop counter

   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

   if (Dest%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Transfer_Point_to_Line2: DEST mesh must have one or more Line2 Elements.'
      RETURN
   endif

   if (Src%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Transfer_Point_to_Line2: SRC mesh must have one or more Point Elements.'
      RETURN
   endif

   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if (    Src%FieldMask(MASKID_TranslationDisp)    &
      .or. Src%FieldMask(MASKID_Orientation)        &
      .or. Src%FieldMask(MASKID_TranslationVel)     &
      .or. Src%FieldMask(MASKID_ROTATIONVEL)        &
      .or. Src%FieldMask(MASKID_TRANSLATIONACC)     &
      .or. Src%FieldMask(MASKID_ROTATIONACC)        &
      .or. Src%FieldMask(MASKID_SCALAR)             &
      ) then

      !........................
      ! Check that the mapping data structure is the correct size:
      !........................

      if (SIZE(MeshMap%MapMotions) < Dest%nnodes) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Transfer_Point_to_Line2: MeshMap%MapMotions(:) should be allocated to Dest%nnodes. '
         RETURN
      endif


      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         ! Each destination node (on a LINE2 mesh) needs a source
         ! in following call, Dest is mesh to looped over, finding a corresponding point for each point in Src
         CALL CreateMapping_NearestNeighbor( Dest, Src, MeshMap%MapMotions, ELEMENT_LINE2, ELEMENT_POINT, ErrStat, ErrMsg )

         ! bjj: for consistant definition of couple_arm (i.e. p_ODR-p_OSR), let's multiply by -1
         do i=1,SIZE(MeshMap%MapMotions)
            MeshMap%MapMotions(i)%couple_arm = -1._ReKi*MeshMap%MapMotions(i)%couple_arm
         end do

         IF (ErrStat >= AbortErrLev) RETURN

      end if

      !........................
      ! Transfer data
      !........................

         ! This is the same algorithm as Transfer_Point_to_Point
      CALL Transfer_Motions_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN


   end if ! algorithm for motions/scalars

   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if (Src%FieldMask(MASKID_FORCE) .or. Src%FieldMask(MASKID_MOMENT)) then

         ErrStat = ErrID_Fatal
         ErrMsg  = " Transfer_Point_to_Line2: Transfer of loads is not implemented."

      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then


      end if


      !........................
      ! Transfer data
      !........................


   end if ! algorithm for loads


END SUBROUTINE Transfer_Point_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcOrient )

   TYPE(MeshType),          INTENT(IN   ) ::  Src
   TYPE(MeshType),          INTENT(INOUT) ::  Dest
   TYPE(MeshType), OPTIONAL,INTENT(IN   ) ::  SrcOrient ! an optional mesh, which contains the orientations associated with the source if the source contains load information

   TYPE(MeshMapType),      INTENT(INOUT) ::  MeshMap     ! Mapping(s) between Src and Dest meshes

   INTEGER(IntKi),         INTENT(  OUT) :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT) :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   ! local variables

   INTEGER(IntKi)  ::   i ! do-loop counter



   ! logic

   ErrStat = ErrID_None
   ErrMsg  = ''

   !.................
   ! Check to ensure that both the source and destination meshes are composed of Point elements
   !.................

   if (Src%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Transfer_Point_to_Point: Src mesh must have one or more Point Elements.'
      RETURN
   endif

   if (Dest%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Transfer_Point_to_Point: Destination mesh must have one or more Point Elements.'
      RETURN
   endif


   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   ! If Src is displacements/velocities then loop over the Destination Mesh;
   !    each motion/scalar in the destination mesh needs to be interpolated from the source mesh.

   if (    Src%FieldMask(MASKID_TranslationDisp)    &
      .or. Src%FieldMask(MASKID_Orientation)        &
      .or. Src%FieldMask(MASKID_TranslationVel)     &
      .or. Src%FieldMask(MASKID_ROTATIONVEL)        &
      .or. Src%FieldMask(MASKID_TRANSLATIONACC)     &
      .or. Src%FieldMask(MASKID_ROTATIONACC)        &
      .or. Src%FieldMask(MASKID_SCALAR)             &
      ) then

      !........................
      ! Check that the mapping data structure is the correct size:
      !........................

      if (SIZE(MeshMap%MapMotions) < Dest%nnodes) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Transfer_Point_to_Point: MeshMap%MapMotions(:) should be allocated to Dest%nnodes.'
         RETURN
      endif


      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         ! in following call, Dest is mesh to looped over, finding a corresponding point for each point in Src
         CALL CreateMapping_NearestNeighbor( Dest, Src, MeshMap%MapMotions, ELEMENT_POINT, ELEMENT_POINT, ErrStat, ErrMsg )

         ! bjj: for consistant definition of couple_arm (i.e. p_ODR-p_OSR), let's multiply by -1
         do i=1,SIZE(MeshMap%MapMotions)
            MeshMap%MapMotions(i)%couple_arm = -1._ReKi*MeshMap%MapMotions(i)%couple_arm
         end do

         IF (ErrStat >= AbortErrLev) RETURN

      end if

      !........................
      ! Transfer data
      !........................

      CALL Transfer_Motions_Point_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN

   end if ! algorithm for motions/scalars

   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   ! IF Src is forces and/or moments, loop over Src Mesh;
   !    each load in the source mesh needs to be placed somewhere in the destination mesh.

   if (Src%FieldMask(MASKID_FORCE) .or. Src%FieldMask(MASKID_MOMENT)) then

      !........................
      ! Check that the mapping data structure is the correct size:
      !........................

      if (SIZE(MeshMap%MapLoads) < Src%nnodes) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Transfer_Point_to_Point: MeshMap%MapLoads(:) should be allocated to Src%nnodes.'
         RETURN
      endif


      !.................
      ! other checks for available mesh fields (now done in alloc mapping routine):
      !.................

      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         ! in following call, Src is mesh to loop over, finding a corresponding point for each point in Dest
         CALL CreateMapping_NearestNeighbor( Src, Dest, MeshMap%MapLoads, ELEMENT_POINT, ELEMENT_POINT, ErrStat, ErrMsg )

         IF (ErrStat >= AbortErrLev) RETURN

      end if

      !........................
      ! Transfer data
      !........................

      IF ( PRESENT( SrcOrient ) ) THEN
         CALL Transfer_Loads_Point_to_Point( Src, Dest, MeshMap%MapLoads, ErrStat, ErrMsg, SrcOrient )
      ELSE
         CALL Transfer_Loads_Point_to_Point( Src, Dest, MeshMap%MapLoads, ErrStat, ErrMsg )
      END IF
      IF (ErrStat >= AbortErrLev) RETURN


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

   INTEGER(IntKi)  :: point_with_min_dist
   INTEGER(IntKi)  :: iElem, iNode, i  ! do-loop counter for elements on Mesh1, associated node(S)
   INTEGER(IntKi)  :: jElem, jNode, j  ! do-loop counter for elements on Mesh2, associated node

   REAL(ReKi)      :: Mesh1_xyz(3)
   REAL(ReKi)      :: Mesh2_xyz(3)

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
   Map(:)%OtherMesh_Element = NODE_NOT_MAPPED ! initialize this so we know if we've mapped this node already (done only because we may have different elements)

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
            ErrStat = ErrID_Fatal
            ErrMsg  = ' Error in CreateMapping_NearestNeighbor: Failed to find destination point associated with source point.'
            RETURN
         endif

         Map(i)%OtherMesh_Element = point_with_min_dist

         Map(i)%distance = min_dist

         Map(i)%couple_arm = Mesh2_xyz - Mesh1_xyz
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
   INTEGER(IntKi)  :: i                                        ! counter over the nodes
   REAL(ReKi)      :: RotationMatrix(3,3)


   ErrStat = ErrID_None
   ErrMsg  = ""


      ! calculate the rotated positions used in several equations below:
   CALL Calculate_RotatedPosition( Src, Dest, MeshMap, ErrStat, ErrMsg )


      ! ---------------------------- Translation ------------------------------------------------

      ! u_Dest = u_Src + [Orientation_Src^T * RefOrientation_Src - I] * [p_Dest - p_Src]
   if ( Src%FieldMask(MASKID_TranslationDisp) .AND. Dest%FieldMask(MASKID_TranslationDisp) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%TranslationDisp(:,i) = Src%TranslationDisp(:,MeshMap%MapMotions(i)%OtherMesh_Element)

            ! if Src mesh has orientation, superpose Dest displacement with translation due to rotation and couple arm
         if ( Src%FieldMask(MASKID_Orientation) ) then

            Dest%TranslationDisp(:,i) = Dest%TranslationDisp(:,i) + MeshMap%RotatedPosition(:,i) - MeshMap%MapMotions(i)%couple_arm

         end if

      end do

   end if

      ! ---------------------------- ORIENTATION/Direction Cosine Matrix   ----------------------

      ! transfer direction cosine matrix, aka orientation

   if ( Src%FieldMask(MASKID_Orientation) .AND. Dest%FieldMask(MASKID_Orientation) ) then

      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         RotationMatrix = MATMUL( Dest%RefOrientation(:,:,i), &
                                  TRANSPOSE( Src%RefOrientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element) ) )
         Dest%Orientation(:,:,i) = MATMUL( RotationMatrix, Src%Orientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element) )
      end do

   endif

      ! ---------------------------- TranslationVel  --------------------------------------------

   if ( Src%FieldMask(MASKID_TranslationVel) .AND. Dest%FieldMask(MASKID_TranslationVel) ) then
      do i=1, Dest%Nnodes
         if ( MeshMap%MapMotions(i)%OtherMesh_Element < 1 )  CYCLE

         Dest%TranslationVel(:,i) = Src%TranslationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element)

         if ( Src%FieldMask(MASKID_RotationVel) ) then
            Dest%TranslationVel(:,i) = Dest%TranslationVel(:,i) + &
                                       cross_product ( Src%RotationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element), &
                                       MeshMap%RotatedPosition(:,i)  )
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
                                       cross_product( Src%RotationAcc(:,MeshMap%MapMotions(i)%OtherMesh_Element), &
                                                      MeshMap%RotatedPosition(:,i) )
         endif

         if ( Src%FieldMask(MASKID_RotationVel) )  then
            Dest%TranslationAcc(:,i) = Dest%TranslationAcc(:,i) + &
                                       cross_product( Src%RotationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element), &
                                       cross_product( Src%RotationVel(:,MeshMap%MapMotions(i)%OtherMesh_Element), &
                                                                        MeshMap%RotatedPosition(:,i) ) )
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
SUBROUTINE Transfer_Loads_Point_to_Point( Src, Dest, Map, ErrStat, ErrMsg, SrcOrient )
! Given a nearest-neighbor mapping, this routine transfers loads between nodes on the mesh.
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src       ! The source mesh with loads fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest      ! The destination mesh
   TYPE(MeshType), OPTIONAL,       INTENT(IN   )  :: SrcOrient ! An optional mesh that contains the orientations associated with the source mesh

   TYPE(MapType),                  INTENT(INOUT)  :: Map(:)    ! The mapping from Dest to Src

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None

      ! local variables
   REAL(ReKi)                                     :: RotationMatrix(3,3)
   REAL(ReKi)                                     :: torque(3)
   INTEGER(IntKi)                                 :: i         ! counter over the nodes


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! initialize the RotationMatrix to Identity
   RotationMatrix(:,1) = (/ 1._ReKi, 0._ReKi, 0._ReKi /)
   RotationMatrix(:,2) = (/ 0._ReKi, 1._ReKi, 0._ReKi /)
   RotationMatrix(:,3) = (/ 0._ReKi, 0._ReKi, 1._ReKi /)


!bjj note that we already checked that the following two conditions apply in this case:
!    if Src%FieldMask(MASKID_FORCE),  Dest%FieldMask(MASKID_FORCE)
!    if Src%FieldMask(MASKID_MOMENT), Dest%FieldMask(MASKID_MOMENT)

   ! initialization; required to handle superposition of forces
   if (Dest%FieldMask(MASKID_FORCE)  ) Dest%Force  = 0.     ! whole array initialization
   if (Dest%FieldMask(MASKID_Moment) ) Dest%Moment = 0.     ! whole array initialization

   do i = 1, Src%NNodes
      if ( Map(i)%OtherMesh_Element < 1 )  CYCLE ! would only happen if we had non-point elements (or nodes not contained in an element)


      if (Src%FieldMask(MASKID_FORCE) ) &
               Dest%Force(:,Map(i)%OtherMesh_Element) = Dest%Force(:,Map(i)%OtherMesh_Element) + Src%Force(:,i)

      if (Dest%FieldMask(MASKID_MOMENT) ) then

         if (Src%FieldMask(MASKID_MOMENT) ) then

            Dest%Moment(:,Map(i)%OtherMesh_Element) = Dest%Moment(:,Map(i)%OtherMesh_Element) + Src%Moment(:,i)

         endif

         ! if the distance (which can never be less than zero) is greater than "zero" and there is a
         ! force in the source mesh, then we need to add a moment to the destination mesh to account
         ! for the mismatch between points

         if (Map(i)%distance .gt. 0_ReKi) then

            if (Src%FieldMask(MASKID_FORCE) ) then

               ! calculation torque vector based on offset force: torque = couple_arm X Force  (depends on sign of couple)

               IF ( PRESENT( SrcOrient) ) THEN
                  RotationMatrix = MATMUL( TRANSPOSE( SrcOrient%Orientation(:,:,i) ), SrcOrient%RefOrientation(:,:,i) )
               ELSEIF ( Src%FieldMask(MASKID_ORIENTATION) ) THEN
                  RotationMatrix = MATMUL( TRANSPOSE( Src%Orientation(:,:,i) ), Src%RefOrientation(:,:,i) )
               END IF

               torque = CROSS_PRODUCT( Src%Force(:,i), MATMUL(RotationMatrix, Map(i)%couple_arm) )
               Dest%Moment(:,Map(i)%OtherMesh_Element) = Dest%Moment(:,Map(i)%OtherMesh_Element) + torque

            endif

         endif

      endif

   enddo

END SUBROUTINE Transfer_Loads_Point_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Calculate_RotatedPosition( Src, Dest, MeshMap, ErrStat, ErrMsg )
! This routine calculates a term used in many of the transformation equations:
!
! MeshMap%RotatedPosition(:,:,i) =
! Src%Orientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element)^T * Src%RefOrientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element) * MeshMap%MapMotions(i)%couple_arm
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src       ! The source mesh with motion fields allocated
   TYPE(MeshType),                 INTENT(IN   )  :: Dest      ! The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   ! data for the mesh mapping

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(ReKi)                                     :: RotationMatrix(3,3)
   INTEGER(IntKi)                                 :: i                       ! loop counter



   ErrStat = ErrID_None
   ErrMsg  = ""


   IF ( Src%FieldMask(MASKID_Orientation) ) THEN

      DO i = 1,Dest%Nnodes

         RotationMatrix = MATMUL( TRANSPOSE(    Src%Orientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element) ) , &
                                             Src%RefOrientation(:,:,MeshMap%MapMotions(i)%OtherMesh_Element) )

         MeshMap%RotatedPosition(:,i) = MATMUL( RotationMatrix, MeshMap%MapMotions(i)%couple_arm )

      END DO

   ELSE

      ! Rotation matrix is Identity (assume Orientation is same as RefOrientation)
      DO i = 1,Dest%Nnodes

         MeshMap%RotatedPosition(:,i) = MeshMap%MapMotions(i)%couple_arm

      END DO

   END IF


END SUBROUTINE Calculate_RotatedPosition
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


   ErrStat = ErrID_None
   ErrMsg  = ''

   if (Src%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Transfer_Line2_to_Line2: Src mesh must have one or more Line2 Elements '
      RETURN
   endif

   if (Dest%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Transfer_Line2_to_Line2: Destination mesh must have one or more Line2 Elements '
      RETURN
   endif


   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Motion and Scalar Fields
   ! ------------------------------------------------------------------------------------------------------------------------------

   if (     Src%FieldMask(MASKID_TranslationDisp)    &
       .or. Src%FieldMask(MASKID_Orientation)        &
       .or. Src%FieldMask(MASKID_TranslationVel)     &
       .or. Src%FieldMask(MASKID_ROTATIONVEL)        &
       .or. Src%FieldMask(MASKID_TRANSLATIONACC)     &
       .or. Src%FieldMask(MASKID_ROTATIONACC)        &
       .or. Src%FieldMask(MASKID_SCALAR)             &
       ) then

      !........................
      ! Check that the mapping data structure is the correct size:
      !........................

      if (SIZE(MeshMap%MapMotions) < Dest%nnodes) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Transfer_Line2_to_Line2: MeshMap%MapMotions(:) should be allocated to Dest%nnodes.'
         RETURN
      endif


      !........................
      ! Start: Create Mapping data (if remap is true)
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then

         CALL CreateMapping_ProjectToLine2(Dest,Src, MeshMap%MapMotions, ELEMENT_LINE2, ErrStat, ErrMsg)
            IF (ErrStat >= AbortErrLev) RETURN

      endif !remapping

      !........................
      ! Start: Transfer data
      !........................
         ! This is the same algorithm as Transfer_Line2_to_Line2
      CALL Transfer_Motions_Line2_to_Point( Src, Dest, MeshMap, ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN

   endif !algorithm for motions/scalars


   ! ------------------------------------------------------------------------------------------------------------------------------
   ! Mapping and Transfer of Data for Mesh Load Fields
   ! ------------------------------------------------------------------------------------------------------------------------------
   if (Src%FieldMask(MASKID_FORCE) .or. Src%FieldMask(MASKID_MOMENT)) then

      !........................
      ! Check that the mapping data structure is the correct size:
      !........................

      !BJJ: THIS IS SUPPOSED TO BE ON THE ELEMENT LEVEL NOW.
      !if (SIZE(MeshMap%MapLoads) < Src%nnodes) then
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = ' Error in Transfer_Line2_to_Line2: MeshMap%MapLoads(:) should be allocated to Src%nnodes.'
      !   RETURN
      !endif

      IF (.not. PRESENT(SrcDisp) .OR. .NOT. PRESENT(DestDisp) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Transfer_Line2_to_Line2: SrcDisp and DestDisp arguments are required.'
         RETURN
      END IF

      !.................
      ! other checks for available mesh fields (now done in AllocMapping routine)
      !.................


      !........................
      ! Create mapping
      !........................

      if (Src%RemapFlag .or. Dest%RemapFlag ) then
         !........................
         ! Create a temporary mesh for lumped point elements of the line2 meshes
         !........................
         CALL Create_PointMesh( Src, MeshMap%Temp_Lumped_Points_Src, ErrStat, ErrMsg )
            IF (ErrStat >= AbortErrLev) RETURN

         CALL Create_PointMesh( Dest, MeshMap%Temp_Lumped_Points_Dest, ErrStat, ErrMsg )
            IF (ErrStat >= AbortErrLev) RETURN


         !........................
         ! Do the mapping
         !........................

         CALL CreateMapping_ProjectLine2ToLine2(Src, Dest, MeshMap%MapLoads, ErrStat, ErrMsg)
            IF (ErrStat >= AbortErrLev) RETURN

      end if

      !........................
      ! Transfer data
      !........................

      CALL Transfer_Loads_Line2_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )


   end if ! algorithm for loads


END SUBROUTINE Transfer_Line2_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Transfer_Loads_Line2_to_Line2( Src, Dest, MeshMap, ErrStat, ErrMsg, SrcDisp, DestDisp )
! Given a mapping, this routine transfers the loads from nodes on Line2 elements to nodes on another Line2 mesh.
!..................................................................................................................................

   TYPE(MeshType),                 INTENT(IN   )  :: Src       ! The source (Line2) mesh with loads fields allocated
   TYPE(MeshType),                 INTENT(INOUT)  :: Dest      ! The destination mesh

   TYPE(MeshMapType),              INTENT(INOUT)  :: MeshMap   ! The mapping data

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat   ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None
   TYPE(MeshType),                 INTENT(IN   )  :: SrcDisp   ! The source mesh's cooresponding position mesh
   TYPE(MeshType),                 INTENT(IN   )  :: DestDisp  ! The destination mesh's cooresponding position mesh
!   TYPE(MeshType), OPTIONAL,       INTENT(IN   )  :: SrcDisp   ! The source mesh's cooresponding position mesh
!   TYPE(MeshType), OPTIONAL,       INTENT(IN   )  :: DestDisp  ! The destination mesh's cooresponding position mesh

      ! local variables
!   REAL(ReKi)                                     :: RotationMatrix(3,3)
!   REAL(ReKi)                                     :: torque(3)

   REAL(ReKi), DIMENSION(3)                       :: xyz1, xyz2, Dest_xyz, tmp_vec ! Position + TranslationDisp for source and destination nodes
   REAL(ReKi)                                     :: Moment_ForceComp(3,Dest%NNodes) ! a temporary holder for the force component of the moment

   REAL(ReKi)                                     :: Scale_F(3)    ! scale for the forces
   REAL(ReKi)                                     :: Scale_M(3)    ! scale for the moment component of the moments
   REAL(ReKi)                                     :: Scale_MF(3)   ! scale for the force component of the moments
   REAL(ReKi), DIMENSION(3)                       :: Scale_F_Denom, Scale_M_Denom, Scale_MF_Denom !temporary vectors to hold values for denominator in calculation of variables Scale_*.
!   REAL(ReKi), DIMENSION(3)                       :: Scale_F_Num, Scale_F_Den, Scale_M_Num, Scale_M_Den, Scale_MF_Num, Scale_MF_Den !temporary vectors to hold values for numerator and denominator

   INTEGER(IntKi)                                 :: i         ! counter over the elements
   INTEGER(IntKi)                                 :: n1, n2    ! temporary space for node numbers

   REAL(ReKi), PARAMETER                          :: TOL = 0.1 ! Error limit: ABS(1-Scale_*) <= TOL [0.1=10% error]

!Possible combinations (if fields exist):
! Src%Force  Src%Moment  Dest%Force Dest%Moment
!----------  ----------  ---------- -----------
!  T             T           T           T
!  T             F           T           T
!  T             F           T           F
!  F             T           F           T
!  F             T           T           T


   ErrStat = ErrID_None
   ErrMsg  = ""


   !   ! set this to the identity matrix for now:
   !RotationMatrix(:,1) = (/ 1._ReKi, 0._ReKi, 0._ReKi /)
   !RotationMatrix(:,2) = (/ 0._ReKi, 1._ReKi, 0._ReKi /)
   !RotationMatrix(:,3) = (/ 0._ReKi, 0._ReKi, 1._ReKi /)


   ! initialization; required to handle superposition of loads
   if (Dest%FieldMask(MASKID_FORCE) ) Dest%Force  = 0._ReKi     ! whole array initialization
   if (Dest%FieldMask(MASKID_MOMENT)) Dest%Moment = 0._ReKi     ! whole array initialization


      ! ---------------------------- Force ------------------------------------------------
      ! Transfer this source force to destination force and/or moment fields:

   Scale_F        = 0.0_ReKi
   Scale_F_Denom  = 0.0_ReKi
   Scale_M        = 0.0_ReKi
   Scale_M_Denom  = 0.0_ReKi
   Scale_MF       = 0.0_ReKi
   Scale_MF_Denom = 0.0_ReKi

   Moment_ForceComp = 0.0_ReKi

   if ( Src%FieldMask(MASKID_Force) ) THEN

      ! ( Dest%FieldMask(MASKID_Force) ) is .TRUE.

      DO i = 1, Src%ElemTable(ELEMENT_LINE2)%nelem  ! the source elements

         IF ( MeshMap%MapLoads(i)%OtherMesh_Element > 0 ) THEN       ! Otherwise, we're ignoring it
            n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
            n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)

            Dest%Force(:,MeshMap%MapLoads(i)%OtherMesh_Element) = Dest%Force(:,    MeshMap%MapLoads(i)%OtherMesh_Element ) &
                                                                 + Src%Force(:,n1)*MeshMap%MapLoads(i)%shape_fn(1)         &
                                                                 + Src%Force(:,n2)*MeshMap%MapLoads(i)%shape_fn(2)
         END IF


      END DO !loop through source elements

     !........................
      ! Lump the values on line2 elements to nodes on a point mesh
      ! These point and line2 meshes have the same nodes so there is no mapping here:
      !........................

      CALL Lump_Line2_to_Point( Src,  MeshMap%Temp_Lumped_Points_Src,  ErrStat, ErrMsg, SrcDisp )
         IF (ErrStat >= AbortErrLev) RETURN

      !ELSE
      !   CALL Lump_Line2_to_Point( Src, MeshMap%Temp_Lumped_Points_Src, ErrStat, ErrMsg )
      !END IF
      !IF (ErrStat >= AbortErrLev) RETURN


         ! Lump the calculated dest%force:
      CALL Lump_Line2_to_Point( Dest, MeshMap%Temp_Lumped_Points_Dest, ErrStat, ErrMsg, DestDisp )
         IF (ErrStat >= AbortErrLev) RETURN


         ! Calculate Scale_F

      ! Compute the numerator of Scale_F
      DO i=1,MeshMap%Temp_Lumped_Points_Src%Nnodes
         Scale_F = Scale_F + MeshMap%Temp_Lumped_Points_Src%Force(:,i)
      END DO

      ! Compute the denominator of Scale_F
      DO i=1,MeshMap%Temp_Lumped_Points_Dest%Nnodes
         Scale_F_Denom = Scale_F_Denom + MeshMap%Temp_Lumped_Points_Dest%Force(:,i)
      END DO

print *, 'Scale_F_Numer=',Scale_F
print *, 'Scale_F_Denom=',Scale_F_Denom
      ! Divide numerator by denominator for Scale_F
      DO i=1,3
         IF ( .NOT. EqualRealNos( Scale_F_Denom(i), 0.0_ReKi ) ) THEN

            Scale_F(i) = Scale_F(i) / Scale_F_Denom(i)

            IF ( ABS(1.0 - Scale_F(i)) > TOL ) THEN
               ErrStat = ErrID_Warn
               IF (LEN_TRIM(ErrMsg)/=0) ErrMsg = TRIM(ErrMsg)//NewLine
               ErrMsg  = TRIM(ErrMsg)//" Transfer_Loads_Line2_to_Line2: Mesh mapping contains large errors (Scale_F("//&
                         TRIM(Num2LStr(i))//") is not close to 1). Consider using a different mesh."
               IF (ErrStat >= AbortErrLev) RETURN
            END IF

         ELSEIF ( EqualRealNos( Scale_F(i), 0.0_ReKi ) ) THEN
            Scale_F(i) = 1.0_ReKi  ! 0/0
         ELSE
            ErrStat = ErrID_Fatal
            IF (LEN_TRIM(ErrMsg)/=0) ErrMsg = TRIM(ErrMsg)//NewLine
            ErrMsg  = TRIM(ErrMsg)//" Transfer_Loads_Line2_to_Line2: Meshes cannot be mapped (Scale_F("//&
                      TRIM(Num2LStr(i))//") denominator is zero)."
            RETURN
         END IF
      END DO

      ! Dest%Force--it exists (check in allocMapping routine)--now contains F^D'

      ! scale the forces:
      Dest%Force(1,:) = Dest%Force(1,:) * Scale_F(1)
      Dest%Force(2,:) = Dest%Force(2,:) * Scale_F(2)
      Dest%Force(3,:) = Dest%Force(3,:) * Scale_F(3)

print *, 'Scale_F=',Scale_F

!bjj: we need to replace SrcDisp%TranslationDisp and DestDisp%TranslationDisp with some checks about if they're present, etc

            ! calculate M^D'_F and Scale_MF  (the force component of the destination moment)
      if ( Dest%FieldMask(MASKID_Moment) ) then

         DO i = 1, Src%ElemTable(ELEMENT_LINE2)%nelem  ! the source elements
print *, 'calculating values: element=', i, ' node=', MeshMap%MapLoads(i)%OtherMesh_Element

            IF ( MeshMap%MapLoads(i)%OtherMesh_Element > 0 ) THEN       ! Otherwise, we're ignoring it

               n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
               n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)

               xyz1  = Src%Position(:,n1) + SrcDisp%TranslationDisp(:,n1)
               xyz2  = Src%Position(:,n2) + SrcDisp%TranslationDisp(:,n2)

               tmp_vec =  xyz1*MeshMap%MapLoads(i)%shape_fn(1)                   &
                        + xyz2*MeshMap%MapLoads(i)%shape_fn(2)                   &
                        - Dest%Position(:,MeshMap%MapLoads(i)%OtherMesh_Element) &
                        - DestDisp%TranslationDisp(:,MeshMap%MapLoads(i)%OtherMesh_Element)
!print *, '[Src%Position(n1) + SrcDisp%TranslationDisp(n1)]*S1=', xyz1*MeshMap%MapLoads(i)%shape_fn(1)
!print *, '[Src%Position(n2) + SrcDisp%TranslationDisp(n2)]*S2=', xyz2*MeshMap%MapLoads(i)%shape_fn(2)
!print *, 'Dest%Position(m) + DestDisp%TranslationDisp(m)=', Dest%Position(:,MeshMap%MapLoads(i)%OtherMesh_Element)  +  DestDisp%TranslationDisp(:,MeshMap%MapLoads(i)%OtherMesh_Element)
!print *, 'tmp_vec=',tmp_vec
!print *,'========'
!
               Moment_ForceComp(:,MeshMap%MapLoads(i)%OtherMesh_Element) =                        &
                                    Moment_ForceComp(:,MeshMap%MapLoads(i)%OtherMesh_Element )    &
                                 + cross_product( tmp_vec,                                        &
                                                  Src%Force(:,n1)*MeshMap%MapLoads(i)%shape_fn(1) &
                                                + Src%Force(:,n2)*MeshMap%MapLoads(i)%shape_fn(2)  )
!print *, Moment_ForceComp
            END IF

         END DO !loop through source elements


         ! Compute the numerator of Scale_MF
         DO i=1,MeshMap%Temp_Lumped_Points_Src%Nnodes
            xyz1     = Src%Position(:,i) + SrcDisp%TranslationDisp(:,i)
            Scale_MF = Scale_MF + cross_product(xyz1,MeshMap%Temp_Lumped_Points_Src%Force(:,i))
         END DO

         ! Compute the denominator of Scale_MF
         DO i=1,MeshMap%Temp_Lumped_Points_Dest%Nnodes
            xyz1           = Dest%Position(:,i) + DestDisp%TranslationDisp(:,i)
            Scale_MF_Denom = Scale_MF_Denom + cross_product(xyz1,MeshMap%Temp_Lumped_Points_Dest%Force(:,i))
         END DO


print *, 'Scale_MF_Numer=',Scale_MF
print *, 'Scale_MF_Denom=',Scale_MF_Denom
         ! calculate Scale_MF
         DO i=1,3
            IF ( .NOT. EqualRealNos( Scale_MF_Denom(i), 0.0_ReKi ) ) THEN

               Scale_MF(i) = Scale_MF(i) / Scale_MF_Denom(i)

               IF ( ABS(1.0 - Scale_MF(i)) > TOL ) THEN
                  ErrStat = ErrID_Warn
                  IF (LEN_TRIM(ErrMsg)/=0) ErrMsg = TRIM(ErrMsg)//NewLine
                  ErrMsg  = TRIM(ErrMsg)//" Transfer_Loads_Line2_to_Line2: Mesh mapping contains large errors (Scale_MF("//&
                            TRIM(Num2LStr(i))//") is not close to 1). Consider using a different mesh."
                  IF (ErrStat >= AbortErrLev) RETURN
               END IF

            ELSEIF ( EqualRealNos( Scale_MF(i), 0.0_ReKi ) ) THEN
               Scale_MF(i) = 1.0_ReKi  ! 0/0
            ELSE
               ErrStat = ErrID_Fatal
               IF (LEN_TRIM(ErrMsg)/=0) ErrMsg = TRIM(ErrMsg)//NewLine
               ErrMsg  = TRIM(ErrMsg)//" Transfer_Loads_Line2_to_Line2: Meshes cannot be mapped (Scale_MF("//&
                         TRIM(Num2LStr(i))//") denominator is zero)."
               RETURN
            END IF
         END DO   ! 3 components of Scale_MF


         Moment_ForceComp(1,:) = Moment_ForceComp(1,:)* Scale_MF(1)
         Moment_ForceComp(2,:) = Moment_ForceComp(2,:)* Scale_MF(2)
         Moment_ForceComp(3,:) = Moment_ForceComp(3,:)* Scale_MF(3)

      end if ! transfer to destination moment field
print *, 'Scale_MF=',Scale_MF


   end if !source has force

      ! ---------------------------- Moments ------------------------------------------------
      ! Transfer this source moment to the destination moment field

   if ( Src%FieldMask(MASKID_Moment) ) then !

      DO i = 1, Src%ElemTable(ELEMENT_LINE2)%nelem  ! the source elements

         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)

            ! Compute the numerator of Scale_M
         Scale_M = Scale_M + Src%ElemTable(ELEMENT_LINE2)%Elements(i)%det_jac * ( Src%Moment(:,n1) + Src%Moment(:,n2) )

         IF ( MeshMap%MapLoads(i)%OtherMesh_Element > 0 ) THEN       ! Otherwise, we're ignoring it

            Dest%Moment(:,MeshMap%MapLoads(i)%OtherMesh_Element) =  Dest%Moment(:,   MeshMap%MapLoads(i)%OtherMesh_Element ) &
                                                                   + Src%Moment(:,n1)*MeshMap%MapLoads(i)%shape_fn(1)        &
                                                                   + Src%Moment(:,n2)*MeshMap%MapLoads(i)%shape_fn(2)
         END IF


      END DO !loop through source elements

         ! Calculate Scale_M:
      DO i=1,Dest%ElemTable(ELEMENT_LINE2)%nelem
         n1 = Dest%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
         n2 = Dest%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)

            ! Compute the denominator of Scale_M
         Scale_M_Denom = Scale_M_Denom + Dest%ElemTable(ELEMENT_LINE2)%Elements(i)%det_jac * &
                           ( Dest%Moment(:,n1) + Dest%Moment(:,n2) )


      END DO

print *, 'Scale_M_Numer=',Scale_M
print *, 'Scale_M_Denom=',Scale_M_Denom
      DO i=1,3
         IF ( .NOT. EqualRealNos( Scale_M_Denom(i), 0.0_ReKi ) ) THEN

            Scale_M(i) = Scale_M(i) / Scale_M_Denom(i)

            IF ( ABS(1.0 - Scale_M(i)) > TOL ) THEN
               ErrStat = ErrID_Warn
               IF (LEN_TRIM(ErrMsg)/=0) ErrMsg = TRIM(ErrMsg)//NewLine
               ErrMsg  = TRIM(ErrMsg)//" Transfer_Loads_Line2_to_Line2: Mesh mapping contains large errors (Scale_M("//&
                         TRIM(Num2LStr(i))//") is not close to 1). Consider using a different mesh."
               IF (ErrStat >= AbortErrLev) RETURN
            END IF

         ELSEIF ( EqualRealNos( Scale_M(i), 0.0_ReKi ) ) THEN
            Scale_M(i) = 1.0_ReKi  ! 0/0
         ELSE
            ErrStat = ErrID_Fatal
            IF (LEN_TRIM(ErrMsg)/=0) ErrMsg = TRIM(ErrMsg)//NewLine
            ErrMsg  = TRIM(ErrMsg)//" Transfer_Loads_Line2_to_Line2: Meshes cannot be mapped (Scale_M("//&
                      TRIM(Num2LStr(i))//") denominator is zero)."
            RETURN
         END IF
      END DO


            ! scale the moments:

      Dest%Moment(1,:) = Dest%Moment(1,:) * Scale_M(1) + Moment_ForceComp(1,:)
      Dest%Moment(2,:) = Dest%Moment(2,:) * Scale_M(2) + Moment_ForceComp(2,:)
      Dest%Moment(3,:) = Dest%Moment(3,:) * Scale_M(3) + Moment_ForceComp(3,:)

   end if !src has moment field

print *, 'Scale_M=',Scale_M

END SUBROUTINE Transfer_Loads_Line2_to_Line2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CreateMapping_ProjectLine2ToLine2(Src, Dest, Map, ErrStat, ErrMsg)
!This routine projects Dest onto Src to find the element mappings between the two meshes. This is used in the mapping for
! Line2-to-Line2 Loads.
!..................................................................................................................................
   TYPE(MeshType),                 INTENT(IN   )  :: Dest                           ! The mesh in the outer mapping loop (Dest)
   TYPE(MeshType),                 INTENT(IN   )  :: Src                            ! The mesh in the inner mapping loop (Src )

   TYPE(MapType),                  INTENT(INOUT)  :: Map(:)                         ! The mapping from Src to Dest

   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat                        ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg                         ! Error message if ErrStat /= ErrID_None

      ! local variables

   REAL(ReKi)      :: denom
   REAL(ReKi)      :: dist
   REAL(ReKi)      :: min_dist
   REAL(ReKi)      :: elem_position

   REAL(ReKi)      :: Dest_xyz(3)
   REAL(ReKi)      :: Src_xyz(3)

   REAL(ReKi)      :: n1_n2_vector(3)     ! vector going from node 1 to node 2 in Line2 element
   REAL(ReKi)      :: n1_Point_vector(3)  ! vector going from node 1 in Line 2 element to Destination Point
   REAL(ReKi)      :: tmp(3)              ! temporary vector for cross product calculation

   INTEGER(IntKi)  :: point_with_min_dist
   INTEGER(IntKi)  :: elem_with_min_dist

   INTEGER(IntKi)  :: iElem, iNode, i  ! do-loop counter for elements on Dest, associated node(S)
   INTEGER(IntKi)  :: jElem, jNode, j  ! do-loop counter for elements on Src, associated node

   INTEGER(IntKi)  :: n1, n2           ! nodes associated with an element

   LOGICAL         :: UseDestNode(Dest%NNodes)    ! determines if the node on the destination mesh is part of the mapping (i.e., contained in an element of the appropriate type)

   LOGICAL         :: on_element



      ! initialization
   ErrStat = ErrID_None
   ErrMsg  = ""


   ! Map the source nodes to destination nodes:
   Map(:)%OtherMesh_Element = NODE_IGNORED ! initialize this so we know if we've mapped this line2 element or if we're going to ignore it

   ! the index, j, of Map(j) corresponds to the source element
   ! the value of Map(j) is the destination node

   ! Determine which nodes on Dest can be in the mapping (we don't want to do the projection and search more than than once per node!)
   UseDestNode = .FALSE.
   do iElem = 1, Dest%ElemTable(ELEMENT_LINE2)%nelem  ! number of point elements on Dest
      do iNode = 1, SIZE( Dest%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes )

         UseDestNode( Dest%ElemTable(ELEMENT_LINE2)%Elements(iElem)%ElemNodes(iNode)    ) = .TRUE. ! use this node (it's part of a line2 element)

      end do
   end do

   ! initialize closest destination node on source mesh:

   Map(:)%distance = HUGE( Map(1)%distance )

   ! project each of the nodes on the destination mesh onto the source elements

   do i = 1, Dest%NNodes   ! a particular node on destination mesh (on a line2 element)
      IF (.NOT. UseDestNode(i) ) CYCLE

      ! point on destination node
      Dest_xyz = Dest%Position(:, i)

      do jElem = 1, Src%ElemTable(ELEMENT_LINE2)%nelem  ! let's see which source elements the projection lies on

            ! grab node numbers associated with the jElem_th element
         n1 = Src%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(1)
         n2 = Src%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(2)

            ! Calculate vectors used in projection operation
         n1_n2_vector    = Src%Position(:,n2) - Src%Position(:,n1)
         n1_Point_vector = Dest_xyz - Src%Position(:,n1)

         denom           = DOT_PRODUCT( n1_n2_vector, n1_n2_vector )
         IF ( EqualRealNos( denom, 0.0_ReKi ) ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = ' Error in CreateMapping_ProjectLine2ToLine2: Division by zero because Line2 element nodes are in same position.'
            RETURN
         END IF

            ! project point onto line defined by n1 and n2
         elem_position = DOT_PRODUCT(n1_n2_vector,n1_Point_vector) / denom

            ! note: i forumlated it this way because Fortran doesn't do shortcutting and I don't want to call EqualRealNos if we don't need it:
         if ( elem_position .le. 1.0_ReKi .and. elem_position .ge. 0.0_ReKi ) then !we're ON the element (between the two nodes)
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
  print *,'dest node=',i, ' source element=',jElem, ' dist=',dist, ' elem_position=',elem_position

            if (dist .lt. Map(jElem)%distance) then

               Map(jElem)%OtherMesh_Element = i
               Map(jElem)%shape_fn(1)       = 1.0_ReKi - elem_position
               Map(jElem)%shape_fn(2)       = elem_position
               Map(jElem)%distance          = dist

               !Map(jElem)%couple_arm        = n1_Point_vector

            end if !the point is closest to this line2 element

         endif !found on the element

      end do !jElem

   end do ! i

   do jElem = 1, Dest%ElemTable(ELEMENT_LINE2)%nelem  ! number of point elements on Dest
 print *,'source element=',jElem, ' dest node=',Map(jElem)%OtherMesh_Element, ' dist=',Map(jElem)%distance, ' elem_position=',Map(jElem)%shape_fn(2)
   end do ! jElem


END SUBROUTINE CreateMapping_ProjectLine2ToLine2
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Lump_Line2_to_Point( Line2_Src, Point_Dest, ErrStat, ErrMsg, SrcTransDisp )
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
   TYPE(MeshType),OPTIONAL,INTENT(IN   ) ::  SrcTransDisp  ! Another mesh that may contain the source's TranslationDisp field


   INTEGER(IntKi),         INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),           INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   REAL(ReKi) :: det_jac
!  REAL(ReKi) :: n1_n2_vector(3) ! vector going from node 1 to node 2 in a Line2 element
   REAL(ReKi) :: pCrossf(3)      ! a temporary vector storing cross product of positions with forces
   REAL(ReKi) :: p1(3), p2(3)    ! temporary position vectors

   INTEGER(IntKi) :: i
   INTEGER(IntKi) :: nnodes

   INTEGER(IntKi) :: n1
   INTEGER(IntKi) :: n2


   ErrStat = ErrID_None
   ErrMsg = ""

   if (Line2_Src%ElemTable(ELEMENT_LINE2)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Lump_Line2_to_Point: Src mesh must have one or more Line2 Elements '
      RETURN
   endif

   if (Point_Dest%ElemTable(ELEMENT_POINT)%nelem .eq. 0) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Lump_Line2_to_Point: dest mesh must have one or more Point Elements '
      RETURN
   endif

   if (Line2_Src%nnodes .ne. Point_Dest%nnodes) then
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in Lump_Line2_to_Point: Src and Dest meshes must have same number of Nodes '
      RETURN
   endif

   if (Point_Dest%FieldMask(MASKID_FORCE) ) then
      if (.not. Line2_Src%FieldMask(MASKID_FORCE) ) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Lump_Line2_to_Point: DEST mesh has FORCE but SRC Mesh Does not '
         RETURN
      endif
   endif

   if (Point_Dest%FieldMask(MASKID_MOMENT) ) then
      if (.not. Line2_SRC%FieldMask(MASKID_MOMENT) ) then
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in Lump_Line2_to_Point: DEST mesh has MOMENT but SRC Mesh Does not '
         RETURN
      endif
   endif

   nnodes = Point_Dest%nnodes  ! also = Line2_Src%nnodes

   ! initialize force/moment in Dest
   if (Point_Dest%FieldMask(MASKID_FORCE) )  Point_Dest%Force  = 0.
   if (Point_Dest%FieldMask(MASKID_MOMENT) ) Point_Dest%Moment = 0.

   ! loop over source mesh, integrating over each element
   do i = 1, Line2_Src%ElemTable(ELEMENT_LINE2)%nelem

      n1 = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)
      n2 = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)

      !n1_n2_vector = Line2_Src%Position(:,n2) - Line2_Src%Position(:,n1)
      !
      !det_jac = SQRT( DOT_PRODUCT(n1_n2_vector,n1_n2_vector) ) / 6._ReKi  ! = L / 6 = det_jac/3.

      det_jac = Line2_Src%ElemTable(ELEMENT_LINE2)%Elements(i)%det_jac / 3.0_ReKi


      if (Point_Dest%FieldMask(MASKID_FORCE) ) then

         !mas
         !det_jac = 0.5 * SQRT( DOT_PRODUCT(n1_n2_vector,n1_n2_vector) )
         !Point_Dest%Force(:,n1) = Point_Dest%Force(:,n1) + det_jac * Line2_Src%Force(:,n1)
         !Point_Dest%Force(:,n2) = Point_Dest%Force(:,n2) + det_jac * Line2_Src%Force(:,n2)

         Point_Dest%Force(:,n1) = Point_Dest%Force(:,n1) + det_jac * ( 2.*Line2_Src%Force(:,n1) +    Line2_Src%Force(:,n2) )
         Point_Dest%Force(:,n2) = Point_Dest%Force(:,n2) + det_jac * (    Line2_Src%Force(:,n1) + 2.*Line2_Src%Force(:,n2) )

      endif

      if (Point_Dest%FieldMask(MASKID_MOMENT) ) then
         !mas
         !det_jac = 0.5 * SQRT( DOT_PRODUCT(n1_n2_vector,n1_n2_vector) )
         !Point_Dest%Moment(:,n1) = Point_Dest%Moment(:,n1) + det_jac * Line2_Src%Moment(:,n1)
         !Point_Dest%Moment(:,n2) = Point_Dest%Moment(:,n2) + det_jac * Line2_Src%Moment(:,n2)

         Point_Dest%Moment(:,n1) = Point_Dest%Moment(:,n1) + det_jac * ( 2.*Line2_Src%Moment(:,n1) +    Line2_Src%Moment(:,n2) )
         Point_Dest%Moment(:,n2) = Point_Dest%Moment(:,n2) + det_jac * (    Line2_Src%Moment(:,n1) + 2.*Line2_Src%Moment(:,n2) )


         if ( Point_Dest%FieldMask(MASKID_FORCE) ) then

            p1 = Line2_Src%Position(:,n1)
            p2 = Line2_Src%Position(:,n2)

            IF ( PRESENT(SrcTransDisp) ) THEN
               p1 = p1 + SrcTransDisp%TranslationDisp(:,n1)
               p2 = p1 + SrcTransDisp%TranslationDisp(:,n2)
            ELSEIF (Line2_Src%FieldMask(MASKID_TRANSLATIONDISP)) THEN
               p1 = p1 + Line2_Src%TranslationDisp(:,n1)
               p2 = p1 + Line2_Src%TranslationDisp(:,n2)
            END IF

            pCrossf = 0.5*det_jac *cross_product( p2-p1, Line2_Src%Force(:,n1) + Line2_Src%Force(:,n2))

            Point_Dest%Moment(:,n1) = Point_Dest%Moment(:,n1) + pCrossf
            Point_Dest%Moment(:,n2) = Point_Dest%Moment(:,n2) - pCrossf

         end if ! src  moment AND force

      endif

   enddo


END SUBROUTINE Lump_Line2_to_Point
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE ModMesh_Mapping
!----------------------------------------------------------------------------------------------------------------------------------
