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
!> This module contains the type definition of of ModMesh, the FAST spatial mesh structure.   
MODULE ModMesh_Types
   USE NWTC_Num
   IMPLICIT NONE

   INTEGER, PUBLIC, PARAMETER :: COMPONENT_INPUT        = 1  !< parameter for "input mesh"
   INTEGER, PUBLIC, PARAMETER :: COMPONENT_OUTPUT       = 2  !< parameter for "output mesh"
   INTEGER, PUBLIC, PARAMETER :: COMPONENT_STATE        = 3  !< parameter for "state mesh" (not recommended to use)

!WARNING... if you add fields here add them to the buffer size computation MeshPack too
   INTEGER, PUBLIC, PARAMETER :: MASKID_FORCE           = 1  !< parameter for fields holding force
   INTEGER, PUBLIC, PARAMETER :: MASKID_MOMENT          = 2  !< parameter for fields holding moment
   INTEGER, PUBLIC, PARAMETER :: MASKID_ORIENTATION     = 3  !< parameter for fields holding orientation
   INTEGER, PUBLIC, PARAMETER :: MASKID_TRANSLATIONDISP = 4  !< parameter for fields holding translational displacement
   INTEGER, PUBLIC, PARAMETER :: MASKID_TRANSLATIONVEL  = 5  !< parameter for fields holding translational velocity
   INTEGER, PUBLIC, PARAMETER :: MASKID_ROTATIONVEL     = 6  !< parameter for fields holding rotational velocity
   INTEGER, PUBLIC, PARAMETER :: MASKID_TRANSLATIONACC  = 7  !< parameter for fields holding translational acceleration
   INTEGER, PUBLIC, PARAMETER :: MASKID_ROTATIONACC     = 8  !< parameter for fields holding rotational acceleration
   INTEGER, PUBLIC, PARAMETER :: MASKID_SCALAR          = 9  !< parameter for fields holding scalars
   INTEGER, PUBLIC, PARAMETER :: FIELDMASK_SIZE         = 9  !< maximum number of fields in a mesh 

   INTEGER, PUBLIC, PARAMETER :: ELEMENT_POINT          = 1  !< parameter for elements of point
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_LINE2          = 2  !< parameter for elements of 2-point lines
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_LINE3          = 3  !< parameter for elements of 3-point lines (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_TRI3           = 4  !< parameter for elements (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_TRI6           = 5  !< parameter for elements (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_QUAD4          = 6  !< parameter for elements (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_QUAD8          = 7  !< parameter for elements (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_TET4           = 8  !< parameter for elements (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_TET10          = 9  !< parameter for elements (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_HEX8           = 10 !< parameter for elements (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_HEX20          = 11 !< parameter for elements (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_WEDGE6         = 12 !< parameter for elements (currently unused)
   INTEGER, PUBLIC, PARAMETER :: ELEMENT_WEDGE15        = 13 !< parameter for elements (currently unused)                    
   INTEGER, PUBLIC, PARAMETER :: NELEMKINDS             = 13 !< parameter for maximum number of element kinds

   CHARACTER(*),DIMENSION(NELEMKINDS),PARAMETER::  ElemNames = (/&
       "Point  ","Line2  ","Line3  ","Tri3   ","Tri6   ",&
       "Quad4  ","Quad8  ","Tet4   ","Tet10  ","Hex8   ",&
       "Hex20  ","Wedge6 ","Wedge15"                   /)   !< element names


   INTEGER, PUBLIC, PARAMETER :: MESH_NEWCOPY         = 1   !< parameter for type of mesh copy: new mesh instance
   INTEGER, PUBLIC, PARAMETER :: MESH_SIBLING         = 2   !< parameter for type of mesh copy: new sibling (shares element and reference data; fields separate)
   INTEGER, PUBLIC, PARAMETER :: MESH_UPDATECOPY      = 3   !< parameter for type of mesh copy: updates fields in existing mesh
   INTEGER, PUBLIC, PARAMETER :: MESH_UPDATEREFERENCE = 4   !< parameter for type of mesh copy: updates reference fields in existing mesh
   INTEGER, PUBLIC, PARAMETER :: MESH_COUSIN          = 5   !< parameter for type of mesh copy: like sibling, but allocates memory for all data

   INTEGER, PUBLIC, PARAMETER :: MESH_NEXT            = -2  !< parameter for next element in mesh
   INTEGER, PUBLIC, PARAMETER :: MESH_NOMOREELEMS     = -3  !< parameter indicating no more elements in mesh
   INTEGER, PUBLIC, PARAMETER :: MESH_NOMOREELEMENTS  = MESH_NOMOREELEMS  !< synonym
   INTEGER, PUBLIC, PARAMETER :: MESH_NOMORE          = MESH_NOMOREELEMS  !< synonym

   LOGICAL, PARAMETER :: mesh_debug = .FALSE.


!   REAL(ReKi), PARAMETER            :: MIN_LINE2_ELEMENT_LENGTH = 0.001 ! 1 millimeter
   REAL(ReKi), PARAMETER            :: MIN_LINE2_ELEMENT_LENGTH = sqrt(epsilon(1.0_ReKi)) ! old length
   
      !> element record type: fields for a particular element
   TYPE, PUBLIC :: ElemRecType 
      ! note: any fields added to this type must be copied in Mesh_MoveAlloc_ElemRecType (modmesh_types::mesh_movealloc_elemrectype)
     INTEGER                        :: Xelement      !< which kind of element
     INTEGER                        :: Nneighbors    !< how many neighbors
     REAL(ReKi)                     :: det_jac       !< determinant of the Jacobian (e.g., 1/2 the length of a line-2 element)
     INTEGER,           ALLOCATABLE :: ElemNodes(:)  !< the list of nodes
!     TYPE(ElemRecType), POINTER     :: Neighbors(:)  ! neighbor list (cannot be allocatable because it's a subtype of itself)
   END TYPE ElemRecType

      !> table of all elements of a particular type
   TYPE, PUBLIC :: ElemTabType   
     INTEGER                     :: nelem                  !< number of elements of this type
     INTEGER                     :: maxelem                !< maximum elements currently allocated for this type 
     INTEGER                     :: Xelement               !< which kind of element
     TYPE(ElemRecType), POINTER  :: Elements(:) => NULL()  !< pointer to list of all elements of this type
   END TYPE ElemTabType

      !> table/list of all elements (may be different types, but not spatial dimensions)
   TYPE, PUBLIC :: ElemListType 
     TYPE(ElemRecType), POINTER :: Element => NULL()        !< pointer to a particular mesh element
   END TYPE ElemListType

      !> mesh data structure
   TYPE, PUBLIC :: MeshType
      LOGICAL :: initialized = .FALSE.                       !< Indicate whether this mesh is initialized
      LOGICAL :: committed = .FALSE.                         !< Indicate whether this mesh is committed
      LOGICAL :: fieldmask(FIELDMASK_SIZE) = .FALSE.         !< Dimension as number of allocatable fields, below
      LOGICAL,POINTER :: RemapFlag  => NULL()                !< false=no action/ignore; true=remap required
      INTEGER :: ios                                         !< Mesh type: input (1), output(2), or state(3)
      INTEGER :: refNode = 0                                 !< optional reference node (informational only)
      INTEGER :: Nnodes = 0                                  !< Number of nodes (vertices) in mesh

     ! Mesh elements
      TYPE(ElemTabType), POINTER :: ElemTable(:) => NULL()   !< A table of all elements in the mesh, by type

     ! Mesh traversal
      INTEGER :: nelemlist                                   !< Number of elements in the list (ElemList)
      INTEGER :: maxelemlist                                 !< Maximum number of elements in the list
      INTEGER :: nextelem                                    !< Next element in the list
      TYPE(ElemListType), POINTER :: ElemList(:) => NULL()   !< All of the elements in the mesh

     ! Node position and reference orientation, which are always allocated (and shared between siblings):
      REAL(ReKi), POINTER :: Position(:,:) => NULL()         !< XYZ coordinate of node (3,:)
      REAL(R8Ki), POINTER :: RefOrientation(:,:,:) => NULL() !< Original/reference orientation [DCM] (3,3,:)

! Here are some built in derived data types that can represent values at the nodes
! the last dimension of each of these has range 1:nnodes for the mesh being represented
! and they are indexed by the element arrays above
! only some of these would be allocated, depending on what's being represented
! on the mesh.
! Whether or not these are allocated is indicated in the fieldmask, which can
! be interrogated by a routine using an instance of the type. If you add a field
! here, be sure to change the table of parameters used to size and index fieldmask above.

      REAL(ReKi), ALLOCATABLE :: Force(:,:)              !< Field: Force vectors (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: Moment(:,:)             !< Field: Moment vectors (3,NNodes)
      REAL(R8Ki), ALLOCATABLE :: Orientation(:,:,:)      !< Field: Direction Cosine Matrix (DCM) (3,3,NNodes)
      REAL(R8Ki), ALLOCATABLE :: TranslationDisp(:,:)    !< Field: Translational displacements (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: RotationVel(:,:)        !< Field: Rotational velocities (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: TranslationVel(:,:)     !< Field: Translational velocities (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: RotationAcc(:,:)        !< Field: Rotational accelerations (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: TranslationAcc(:,:)     !< Field: Translational accelerations (3,NNodes)
      REAL(ReKi), ALLOCATABLE :: Scalars(:,:)            !< Scalars (nScalars,NNodes)
      INTEGER                 :: nScalars                !< Stores value of nScalars when created
!bjj: to be added later:       REAL(ReKi), ALLOCATABLE :: ElementScalars(nElementScalars,nelemlist)       ! Scalars associated with elements
!bjj: to be added later:       INTEGER             :: nElementScalars               ! Stores value of nElementScalars when created

     ! Keeping track of siblings:
      TYPE(MeshType),POINTER :: SiblingMesh      => NULL() !< Pointer to mesh's (only) sibling

   END TYPE MeshType

CONTAINS
!> This function returns the number of nodes in a given type of element.
   INTEGER FUNCTION NumNodes( Xelement )
     INTEGER, INTENT(IN) :: Xelement !< type of element

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

!> This function determines if a mesh contains any motion field (translational/rotational positions, velocities, accelerations or scalars).
   LOGICAL FUNCTION HasMotionFields(Mesh)

      TYPE(MeshType), INTENT(IN) :: Mesh  !< mesh to query

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

!> This function determines if a mesh contains any load field (force or motion).
   LOGICAL FUNCTION HasLoadFields(Mesh)

      TYPE(MeshType), INTENT(IN) :: Mesh !< mesh to query

      IF (     Mesh%FieldMask(MASKID_FORCE)   &
          .or. Mesh%FieldMask(MASKID_MOMENT)  &
          ) THEN

         HasLoadFields = .TRUE.
      ELSE
         HasLoadFields = .FALSE.
      END IF

   END FUNCTION HasLoadFields

!> This subroutine copies the element record data from one ElemRecType data structure to another. It calls the Fortran 2003 
!! intrinsic MOVE_ALLOC routine to move the address of the Src\%ElemNodes array to the Dest\%ElemNodes array without physically
!! copying any of the array. On exist Src\%ElemNodes will be deallocated. 
   SUBROUTINE Mesh_MoveAlloc_ElemRecType( Src, Dest )
      TYPE(ElemRecType), INTENT(INOUT) :: Src   !< mesh containing ElemNodes to be moved
      TYPE(ElemRecType), INTENT(INOUT) :: Dest  !< mesh that will receive the ElemNodes array from Src 
   
      Dest%Xelement   = Src%Xelement
      Dest%Nneighbors = Src%Nneighbors
      Dest%det_jac    = Src%det_jac
      if (allocated(Src%ElemNodes)) &   ! bjj: 9/12/15 added this because of invalid memory address (harmless?) found with Inspector
      CALL Move_Alloc( Src%ElemNodes,  Dest%ElemNodes )
             
   END SUBROUTINE Mesh_MoveAlloc_ElemRecType
END MODULE ModMesh_Types
