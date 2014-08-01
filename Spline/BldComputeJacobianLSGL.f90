   SUBROUTINE BldComputeJacobianLSGL(rr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,jacobian)
   !_------------------------------------------------------------------------------------------------
   ! This subroutine 1) computes the jacobian of a element;
   !                 2) adjusts derivative of shape functions.
   ! For details, see
   ! Bauchau, O.A., "Flexible Multibody Dynamics", Springer, pp. 643
   !-------------------------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN)::  rr            ! rrth Gauss point location 
   REAL(ReKi),    INTENT(IN)::  Nuu0(:)       ! Element nodal initial position
   REAL(ReKi),    INTENT(IN)::  gp(:)         ! Gauss point location
   REAL(ReKi),    INTENT(IN)::  GLL_temp(:)   ! Gauss-Lobatto-Legendre point location
   REAL(ReKi),    INTENT(OUT):: jacobian      ! Jacobian of element
   REAL(ReKi),    INTENT(OUT):: hhx(:)        ! Shape function
   REAL(ReKi),    INTENT(OUT):: hpx(:)        ! Derivative of shape function
   INTEGER(IntKi),INTENT(IN)::  node_elem     ! Number of node per element
   INTEGER(IntKi),INTENT(IN)::  dof_node      ! Number of DoF per node
   INTEGER(IntKi),INTENT(IN)::  ngp           ! Total number of Gauss point
   INTEGER(IntKi),INTENT(IN)::  igp           ! ith Gauss point


   ! Local variables
   REAL(ReKi)::Gup0(3)
   INTEGER(IntKi)::inode
   INTEGER(IntKi)::temp_id
   INTEGER(IntKi)::i

   hhx = 0.0D0
   hpx = 0.0D0
   CALL diffmtc(node_elem-1,ngp,gp,GLL_temp,igp,hhx,hpx)

   Gup0 = 0.0D0
   DO inode=1,node_elem
       temp_id = (inode-1)*dof_node
       DO i=1,3
           Gup0(i) = Gup0(i) + hpx(inode)*Nuu0(temp_id+i)
       ENDDO
   ENDDO

   jacobian = 0.0D0
   jacobian = SQRT(DOT_PRODUCT(Gup0,Gup0))
   
   DO inode=1,node_elem
       hpx(inode) = hpx(inode)/jacobian
   ENDDO

   END SUBROUTINE BldComputeJacobianLSGL
