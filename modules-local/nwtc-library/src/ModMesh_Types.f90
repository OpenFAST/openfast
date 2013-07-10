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
   INTEGER, PUBLIC, PARAMETER :: MASKID_ADDEDMASS = 9
   INTEGER, PUBLIC, PARAMETER :: MASKID_SCALAR = 10
   INTEGER, PUBLIC, PARAMETER :: FIELDMASK_SIZE = 11 !bjj should this be 10?

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


   INTEGER, PUBLIC, PARAMETER :: MESH_NEWCOPY  = 1
   INTEGER, PUBLIC, PARAMETER :: MESH_SIBLING  = 2
   INTEGER, PUBLIC, PARAMETER :: MESH_UPDATECOPY  = 3

   INTEGER, PUBLIC, PARAMETER :: MESH_NEXT  = -2
   INTEGER, PUBLIC, PARAMETER :: MESH_NOMOREELEMS  = -3
   INTEGER, PUBLIC, PARAMETER :: MESH_NOMOREELEMENTS  = MESH_NOMOREELEMS  !synonym
   INTEGER, PUBLIC, PARAMETER :: MESH_NOMORE  = MESH_NOMOREELEMS          !synonym

   LOGICAL, PARAMETER :: mesh_debug = .FALSE.


   TYPE, PUBLIC :: ElemRecType ! a particular element
     INTEGER                     :: Xelement               ! which kind of element
     INTEGER                     :: Nneighbors             ! how many neighbors
     INTEGER, POINTER            :: ElemNodes(:) => NULL() ! pointer to the list of nodes
     TYPE(ElemRecType), POINTER  :: Neighbors(:) => NULL() ! neighbor list
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
      LOGICAL :: initialized = .FALSE.                ! Indicate whether this mesh is initialized     (bjj: added default initialization here)
      LOGICAL :: fieldmask(FIELDMASK_SIZE) = .FALSE.  ! Dimension as number of allocateable fields, below
      LOGICAL,POINTER :: RemapFlag  => NULL()         ! false=no action/ignore; true=remap required
      INTEGER :: ios                                  ! Mesh type: input (1), output(2), or state(3)
      INTEGER :: Nnodes                               ! Number of nodes (vertices) in mesh

      TYPE(ElemTabType), POINTER :: ElemTable(:) => NULL() ! Elements in the mesh, by type

     ! Mesh traversal
      INTEGER :: nelemlist                ! Number of elements in the list (ElemList)
      INTEGER :: maxelemlist              ! Maximum number of elements in the list
      INTEGER :: nextelem                 ! Next element in the list
      INTEGER :: spatial !bjj: unused?
      TYPE(ElemListType), POINTER :: ElemList(:) => NULL() ! All of the elements in the mesh

! Here are some built in derived data types that can represent values at the nodes
! the last dimension of each of these has range 1:nnodes for the mesh being represented
! and they are indexed by the element arrays above
! only some of these would be allocated, depending on what's being represented
! on the mesh.
! Whether or not these are allocated is indicated in the fieldmask, which can
! be interrogated by a routine using an instance of the type. If you add a field
! here, be sure to change the table of parameters used to size and index fieldmask above.
      REAL(ReKi), POINTER     :: Position(:,:) => NULL() ! XYZ coordinate of node (always allocated) (3,:)
      REAL(ReKi), ALLOCATABLE :: Force(:,:)              ! Field: Force vectors (3,:)
      REAL(ReKi), ALLOCATABLE :: Moment(:,:)             ! Field: Moment vectors (3,:)
      REAL(ReKi), ALLOCATABLE :: Orientation(:,:,:)      ! Field: Direction Cosine Matrix (DCM) (3,3,:)
      REAL(ReKi), ALLOCATABLE :: TranslationDisp(:,:)    ! Field: Translational displacements (3,:)
      REAL(ReKi), ALLOCATABLE :: RotationVel(:,:)        ! Field: Rotational velocities (3,:)
      REAL(ReKi), ALLOCATABLE :: TranslationVel(:,:)     ! Field: Translational velocities (3,:)
      REAL(ReKi), ALLOCATABLE :: RotationAcc(:,:)        ! Field: Rotational accelerations (3,:)
      REAL(ReKi), ALLOCATABLE :: TranslationAcc(:,:)     ! Field: Translational accelerations (3,:)
      REAL(ReKi), ALLOCATABLE :: AddedMass(:,:,:)        ! Field: Added mass matrix (6,6,:)
      REAL(ReKi), ALLOCATABLE :: Scalars(:,:)            ! Scalars (1st Dim is over Scalars; 2nd is nodes)
      INTEGER                 :: nScalars                ! store value of nScalars when created
!bjj: to be added later:       REAL(ReKi), ALLOCATABLE :: ElementScalars(:,:)       ! Scalars associated with elements
!bjj: to be added later:       INTEGER             :: nElementScalars               ! store value of nElementScalars when created

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
         CALL ProgAbort('NumNodes: internal error: invalid argument '//TRIM(Num2LStr(Xelement)))
     END SELECT
   END FUNCTION NumNodes

END MODULE ModMesh_Types
