!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2014  National Renewable Energy Laboratory
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
MODULE ModMesh_Types
   USE NWTC_Num
   IMPLICIT NONE

   INTEGER, PUBLIC, PARAMETER :: COMPONENT_INPUT = 1
   INTEGER, PUBLIC, PARAMETER :: COMPONENT_OUTPUT = 2
   INTEGER, PUBLIC, PARAMETER :: COMPONENT_STATE = 3

!WARNING... if you add fields here add them to the buffer size computation MeshPack too
   INTEGER, PUBLIC, PARAMETER :: MASKID_FORCE = 1
   INTEGER, PUBLIC, PARAMETER :: MASKID_MOMENT = 2
   INTEGER, PUBLIC, PARAMETER :: MASKID_ORIENTATION = 3
   INTEGER, PUBLIC, PARAMETER :: MASKID_TRANSLATIONDISP = 4
   INTEGER, PUBLIC, PARAMETER :: MASKID_TRANSLATIONVEL = 5
   INTEGER, PUBLIC, PARAMETER :: MASKID_ROTATIONVEL = 6
   INTEGER, PUBLIC, PARAMETER :: MASKID_TRANSLATIONACC = 7
   INTEGER, PUBLIC, PARAMETER :: MASKID_ROTATIONACC = 8
   INTEGER, PUBLIC, PARAMETER :: MASKID_SCALAR = 9
   INTEGER, PUBLIC, PARAMETER :: FIELDMASK_SIZE = 9

! Format of the Int buffer
   INTEGER, PUBLIC, PARAMETER :: HDR_INTBUFSIZE  = 1
   INTEGER, PUBLIC, PARAMETER :: HDR_REALBUFSIZE = 2
   INTEGER, PUBLIC, PARAMETER :: HDR_DBLBUFSIZE  = 3
   INTEGER, PUBLIC, PARAMETER :: HDR_IOS         = 4
   INTEGER, PUBLIC, PARAMETER :: HDR_NUMNODES    = 5
   INTEGER, PUBLIC, PARAMETER :: HDR_NUMELEMREC  = 6
   INTEGER, PUBLIC, PARAMETER :: HDR_FIELDMASK   = 7
   INTEGER, PUBLIC, PARAMETER :: HDR_FIXEDLEN    = HDR_FIELDMASK+FIELDMASK_SIZE-1
   INTEGER, PUBLIC, PARAMETER :: HDR_FIRSTELEM   = HDR_FIXEDLEN+1

   INTEGER, PUBLIC, PARAMETER :: ELEMENT_POINT   = 1
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_LINE2   = 2
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_LINE3   = 3
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_TRI3    = 4
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_TRI6    = 5
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_QUAD4   = 6
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_QUAD8   = 7
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_TET4    = 8
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_TET10   = 9
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_HEX8    = 10
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_HEX20   = 11
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_WEDGE6  = 12
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_WEDGE15 = 13
   INTEGER, PUBLIC, PARAMETER :: NELEMKINDS = 13

   CHARACTER(*),DIMENSION(NELEMKINDS),PARAMETER::  ElemNames = (/&
       "Point  ","Line2  ","Line3  ","Tri3   ","Tri6   ",&
       "Quad4  ","Quad8  ","Tet4   ","Tet10  ","Hex8   ",&
       "Hex20  ","Wedge6 ","Wedge15"                   /)


   INTEGER, PUBLIC, PARAMETER :: MESH_NEWCOPY         = 1
   INTEGER, PUBLIC, PARAMETER :: MESH_SIBLING         = 2
   INTEGER, PUBLIC, PARAMETER :: MESH_UPDATECOPY      = 3
   INTEGER, PUBLIC, PARAMETER :: MESH_UPDATEREFERENCE = 4

   INTEGER, PUBLIC, PARAMETER :: MESH_NEXT  = -2
   INTEGER, PUBLIC, PARAMETER :: MESH_NOMOREELEMS  = -3
   INTEGER, PUBLIC, PARAMETER :: MESH_NOMOREELEMENTS  = MESH_NOMOREELEMS  !synonym
   INTEGER, PUBLIC, PARAMETER :: MESH_NOMORE  = MESH_NOMOREELEMS          !synonym

   LOGICAL, PARAMETER :: mesh_debug = .FALSE.


   TYPE, PUBLIC :: ElemRecType ! a particular element
      ! note: any fields added to this type must be copied in Mesh_MoveAlloc_ElemRecType()
     INTEGER                        :: Xelement      ! which kind of element
     INTEGER                        :: Nneighbors    ! how many neighbors
     REAL(ReKi)                     :: det_jac       ! determinant of the Jacobian (e.g., 1/2 the length of a line-2 element)
     INTEGER,           ALLOCATABLE :: ElemNodes(:)  ! the list of nodes
!     TYPE(ElemRecType), POINTER     :: Neighbors(:)  ! neighbor list (cannot be allocatable because it's a subtype of itself)
   END TYPE ElemRecType

   TYPE, PUBLIC :: ElemTabType   ! table of all elements of a particular type
     INTEGER                     :: nelem, maxelem
     INTEGER                     :: Xelement               ! which kind of element
     TYPE(ElemRecType), POINTER  :: Elements(:) => NULL()
   END TYPE ElemTabType

   TYPE, PUBLIC :: ElemListType ! all elements (may be different types, but not spatial dimensions)
     TYPE(ElemRecType), POINTER :: Element => NULL()                 ! bjj added initialization
   END TYPE ElemListType

   TYPE, PUBLIC :: MeshType
      LOGICAL :: initialized = .FALSE.                       ! Indicate whether this mesh is initialized
      LOGICAL :: committed = .FALSE.                         ! Indicate whether this mesh is committed
      LOGICAL :: fieldmask(FIELDMASK_SIZE) = .FALSE.         ! Dimension as number of allocatable fields, below
      LOGICAL,POINTER :: RemapFlag  => NULL()                ! false=no action/ignore; true=remap required
      INTEGER :: ios                                         ! Mesh type: input (1), output(2), or state(3)
      INTEGER :: Nnodes = 0                                  ! Number of nodes (vertices) in mesh

     ! Mesh elements
      TYPE(ElemTabType), POINTER :: ElemTable(:) => NULL()   ! Elements in the mesh, by type

     ! Mesh traversal
      INTEGER :: nelemlist                                   ! Number of elements in the list (ElemList)
      INTEGER :: maxelemlist                                 ! Maximum number of elements in the list
      INTEGER :: nextelem                                    ! Next element in the list
      TYPE(ElemListType), POINTER :: ElemList(:) => NULL()   ! All of the elements in the mesh

     ! Node position and reference orientation, which are always allocated (and shared between siblings):
      REAL(ReKi), POINTER :: Position(:,:) => NULL()         ! XYZ coordinate of node (3,:)
      REAL(ReKi), POINTER :: RefOrientation(:,:,:) => NULL() ! Original/reference orientation [DCM] (3,3,:)

! Here are some built in derived data types that can represent values at the nodes
! the last dimension of each of these has range 1:nnodes for the mesh being represented
! and they are indexed by the element arrays above
! only some of these would be allocated, depending on what's being represented
! on the mesh.
! Whether or not these are allocated is indicated in the fieldmask, which can
! be interrogated by a routine using an instance of the type. If you add a field
! here, be sure to change the table of parameters used to size and index fieldmask above.

      REAL(ReKi), ALLOCATABLE :: Force(:,:)              ! Field: Force vectors (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: Moment(:,:)             ! Field: Moment vectors (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: Orientation(:,:,:)      ! Field: Direction Cosine Matrix (DCM) (3,3,NNodes)
      REAL(ReKi), ALLOCATABLE :: TranslationDisp(:,:)    ! Field: Translational displacements (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: RotationVel(:,:)        ! Field: Rotational velocities (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: TranslationVel(:,:)     ! Field: Translational velocities (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: RotationAcc(:,:)        ! Field: Rotational accelerations (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: TranslationAcc(:,:)     ! Field: Translational accelerations (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: Scalars(:,:)            ! Scalars (nScalars,NNodes)
      INTEGER                 :: nScalars                ! Stores value of nScalars when created
!bjj: to be added later:       REAL(ReKi), ALLOCATABLE :: ElementScalars(nElementScalars,nelemlist)       ! Scalars associated with elements
!bjj: to be added later:       INTEGER             :: nElementScalars               ! Stores value of nElementScalars when created

     ! Keeping track of siblings:
      TYPE(MeshType),POINTER :: SiblingMesh      => NULL() ! Pointer to mesh's (only) sibling

   END TYPE MeshType

CONTAINS
   INTEGER FUNCTION NumNodes( Xelement )
     INTEGER, INTENT(IN) :: Xelement

     SELECT CASE ( Xelement )
       CASE ( ELEMENT_POINT )
         NumNodes = 1
       CASE ( ELEMENT_LINE2 )
         NumNodes = 2
       CASE ( ELEMENT_LINE3 )
         NumNodes = 3
       CASE ( ELEMENT_TRI3 )
         NumNodes = 3
       CASE ( ELEMENT_TRI6 )
         NumNodes = 6
       CASE ( ELEMENT_QUAD4 )
         NumNodes = 4
       CASE ( ELEMENT_QUAD8 )
         NumNodes = 8
       CASE ( ELEMENT_TET4 )
         NumNodes = 4
       CASE ( ELEMENT_TET10 )
         NumNodes = 10
       CASE ( ELEMENT_HEX8 )
         NumNodes = 8
       CASE ( ELEMENT_HEX20 )
         NumNodes = 20
       CASE ( ELEMENT_WEDGE6 )
         NumNodes = 6
       CASE ( ELEMENT_WEDGE15 )
         NumNodes = 15
       CASE DEFAULT
          NumNodes = 0
         CALL ProgAbort(' NumNodes: invalid argument Xelement = '//TRIM(Num2LStr(Xelement)))
     END SELECT
   END FUNCTION NumNodes

   LOGICAL FUNCTION HasMotionFields(Mesh)

      TYPE(MeshType), INTENT(IN) :: Mesh

      IF (     Mesh%FieldMask(MASKID_TRANSLATIONDISP)    &
          .or. Mesh%FieldMask(MASKID_ORIENTATION)        &
          .or. Mesh%FieldMask(MASKID_TRANSLATIONVEL)     &
          .or. Mesh%FieldMask(MASKID_ROTATIONVEL)        &
          .or. Mesh%FieldMask(MASKID_TRANSLATIONACC)     &
          .or. Mesh%FieldMask(MASKID_ROTATIONACC)        &
          .or. Mesh%FieldMask(MASKID_SCALAR)             &
          ) THEN

         HasMotionFields = .TRUE.
      ELSE
         HasMotionFields = .FALSE.
      END IF

   END FUNCTION HasMotionFields

   LOGICAL FUNCTION HasLoadFields(Mesh)

      TYPE(MeshType), INTENT(IN) :: Mesh

      IF (     Mesh%FieldMask(MASKID_FORCE)   &
          .or. Mesh%FieldMask(MASKID_MOMENT)  &
          ) THEN

         HasLoadFields = .TRUE.
      ELSE
         HasLoadFields = .FALSE.
      END IF

   END FUNCTION HasLoadFields

   SUBROUTINE Mesh_MoveAlloc_ElemRecType( Src, Dest )
      TYPE(ElemRecType), INTENT(INOUT) :: Src
      TYPE(ElemRecType), INTENT(INOUT) :: Dest
   
      Dest%Xelement   = Src%Xelement
      Dest%Nneighbors = Src%Nneighbors
      Dest%det_jac    = Src%det_jac
      CALL Move_Alloc( Src%ElemNodes,  Dest%ElemNodes )
             
   END SUBROUTINE Mesh_MoveAlloc_ElemRecType
END MODULE ModMesh_Types
