MODULE ModMesh
! (c) 2012 National Renewable Energy Laboratory
!
!BJJ: This is a place holder for a module John Michalaches and Ilene Carpenter are writing.
! This will be part of the NWTC Subroutine Library.


 USE PRECISION

   INTEGER(IntKi), PARAMETER :: MESH_NEWCOPY    = 1
   INTEGER(IntKi), PARAMETER :: MESH_SIBLING    = 2
   INTEGER(IntKi), PARAMETER :: MESH_UPDATECOPY = 3


TYPE, PUBLIC :: MeshType
   LOGICAL                 :: committed             ! Indicate whether this mesh is committed
   INTEGER(IntKi)          :: ios                   ! COMPONENT_INPUT/OUTPUT/STATE/PARAMETER
   INTEGER(IntKi)          :: RemapFlag             ! Remap flag: zero=noaction/ignore
                                                    !             nonzero=module_defined
   INTEGER(IntKi)          :: Nnodes                ! Number of nodes (vertices) in mesh
   INTEGER(IntKi)          :: Nelements             ! Number of elements in mesh
   INTEGER(IntKi)          :: Npoint                ! Number of point elements
   INTEGER(IntKi)          :: Nline2                ! Number of 2-node line elements
   INTEGER(IntKi)          :: Nline3                ! Number of 3-node line elements
   INTEGER(IntKi)          :: Ntri3                 ! Number of 3-node triangle elements
   INTEGER(IntKi)          :: Ntri6                 ! Number of 6-node triangle elements
   INTEGER(IntKi)          :: Nquad4                ! Number of 4-node quadrilateral elements
   INTEGER(IntKi)          :: Nquad8                ! Number of 8-node quadrilateral elements
   INTEGER(IntKi)          :: Ntet4                 ! Number of 4-node tet elements
   INTEGER(IntKi)          :: Ntet10                ! Number of 10-node tet elements
   INTEGER(IntKi)          :: Nhex8                 ! Number of 8-node hex elements
   INTEGER(IntKi)          :: Nhex20                ! Number of 20-node hex elements
   INTEGER(IntKi)          :: Nwedge6               ! Number of 6-node wedge elements
   INTEGER(IntKi)          :: Nwedge15              ! Number of 15-node wedgeelements
   INTEGER(IntKi), POINTER :: element_point(:)      ! Point connectivity
   INTEGER(IntKi), POINTER :: element_line2(:,:)    ! 2-node line connectivity
   INTEGER(IntKi), POINTER :: element_line3(:,:)    ! 3-node line connectivity
   INTEGER(IntKi), POINTER :: element_tri3(:,:)     ! 3-node triangle connectivity
   INTEGER(IntKi), POINTER :: element_tri6(:,:)     ! 6-node triangle connectivity
   INTEGER(IntKi), POINTER :: element_quad4(:,:)    ! 4-node quad connectivity
   INTEGER(IntKi), POINTER :: element_quad8(:,:)    ! 8-node quad connectivity
   INTEGER(IntKi), POINTER :: element_tet4(:,:)     ! 4-node tet connectivity
   INTEGER(IntKi), POINTER :: element_tet10(:,:)   ! 10-node tet connectivity
   INTEGER(IntKi), POINTER :: element_hex8(:,:)     ! 8-node hex connectivity
   INTEGER(IntKi), POINTER :: element_hex20(:,:)   ! 20-node hex connectivity
   INTEGER(IntKi), POINTER :: element_wedge6(:,:)   ! 6-node wedge connectivity
   INTEGER(IntKi), POINTER :: element_wedge15(:,:) ! 15-node wedge connectivity
   REAL(ReKi),     POINTER :: Position(:,:)         ! XYZ coordinate of node
   REAL(ReKi),     POINTER :: Force(:,:)            ! Force vectors
   REAL(ReKi),     POINTER :: Moment(:,:)           ! Moment vectors
   REAL(ReKi),     POINTER :: Orientation(:,:,:)    ! Direction Cosine Matrix (DCM)
   REAL(ReKi),     POINTER :: Rotation(:,:)         ! Rotational Velocities
   REAL(ReKi),     POINTER :: Translation(:,:)      ! Translational Velocities
   REAL(ReKi),     POINTER :: AddedMass(:,:,:)      ! Added mass matrix
   REAL(ReKi),     POINTER :: Scalars(:,:)          ! Scalars (2nd Dim is over Scalars)
   TYPE(MeshType), POINTER :: YoungerSibling        ! Pointer to next sibling in list
   TYPE(MeshType), POINTER :: ElderSibling          ! Pointer to prev sibling in list
END TYPE MeshType


END MODULE
