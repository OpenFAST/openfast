MODULE ModMesh_Types
   USE Precision                                            ! This should be unnecessary if NWTC_IO is used.    - MLB
!   USE NWTC_Library
   USE NWTC_IO
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
   INTEGER, PUBLIC, PARAMETER :: FIELDMASK_SIZE = 11

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

   LOGICAL :: mesh_debug = .FALSE.


   TYPE, PUBLIC :: ElemRecType
     INTEGER                     :: Xelement     ! which kind of element
     INTEGER                     :: Nneighbors   ! how many neighbors
     INTEGER, POINTER            :: ElemNodes(:) => NULL() ! pointer to the list of nodes
     TYPE(ElemRecType), POINTER  :: Neighbors(:) => NULL() ! neighbor list
   END TYPE ElemRecType

   TYPE, PUBLIC :: ElemTabType
     INTEGER                         :: nelem, maxelem
     INTEGER                         :: Xelement     ! which kind of element
     TYPE(ElemRecType), POINTER      :: Elements(:) => NULL()
   END TYPE ElemTabType

   TYPE, PUBLIC :: ElemListType
     TYPE(ElemRecType), POINTER :: Element
   END TYPE ElemListType

   TYPE, PUBLIC :: MeshType
      LOGICAL initialized              ! indicate whether this mesh is initialized
      LOGICAL fieldmask(FIELDMASK_SIZE)! dimension as number of allocateable fields, below
      LOGICAL RemapFlag                ! false=no action/ignore; true=remap required
      INTEGER ios                      ! input (1), output(2), or state(3)
      INTEGER Nnodes                   ! Number of nodes (vertices) in mesh

      TYPE(ElemTabType), POINTER :: ElemTable(:) => NULL()

     ! Mesh traversal
      INTEGER nelemlist, maxelemlist, nextelem
      INTEGER spatial
      TYPE(ElemListType), pointer :: ElemList(:) => NULL()

! Here are some built in derived data types that can represent values at the nodes
! the last dimension of each of these has range 1:nnodes for the mesh being represented
! and they are indexed by the element arrays above
! only some of these would be allocated, depending on what's being represented
! on the mesh.
! Whether or not these are allocated is indicted in the fieldmask, which can
! be interrogated by a routine using an instance of the type. If you add a field
! here, be sure to change the table of parameters used to size and index fieldmask above.
      REAL(ReKi), pointer :: Position(:,:)  => NULL()     ! XYZ coordinate of node (always allocated) (3,:)
      REAL(ReKi), pointer :: Force(:,:)     => NULL()     !  (3,:)
      REAL(ReKi), pointer :: Moment(:,:)    => NULL()     !  (3,:)
      REAL(ReKi), pointer :: Orientation(:,:,:) => NULL() ! Direction Cosine Matrix (DCM) (3,3,:)
      REAL(ReKi), pointer :: TranslationDisp(:,:)   => NULL() ! Translational Displacements  (3,:)
      REAL(ReKi), pointer :: RotationVel(:,:)      => NULL() ! Rotational Velocities  (3,:)
      REAL(ReKi), pointer :: TranslationVel(:,:)   => NULL() ! Translational Velocities  (3,:)
      REAL(ReKi), pointer :: RotationAcc(:,:)      => NULL() ! Rotational Accelerations  (3,:)
      REAL(ReKi), pointer :: TranslationAcc(:,:)   => NULL() ! Translational Accelerations  (3,:)
      REAL(ReKi), pointer :: AddedMass(:,:,:)   => NULL() ! AddedMass (6,6,:)
      REAL(ReKi), pointer :: Scalars(:,:)       => NULL() ! Scalars (1st Dim is over Scalars)
      INTEGER             :: nScalars                     ! store value of nScalars when created

      TYPE(MeshType),POINTER :: SiblingMesh  => NULL()

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
         write(0,*)'NumNodes: internal error: invalid argument ', Xelement
         call progabort('')
     END SELECT
   END FUNCTION NumNodes

END MODULE ModMesh_Types
