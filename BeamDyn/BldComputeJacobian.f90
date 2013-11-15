   SUBROUTINE BldComputeJacobian(rr,Nuu0,node_elem,dof_node,hhx,hpx,jacobian)
   REAL(ReKi),INTENT(IN):: rr,Nuu0(:)
   REAL(ReKi),INTENT(OUT):: jacobian
   REAL(ReKi),INTENT(OUT):: hhx(:),hpx(:)
   INTEGER(IntKi),INTENT(IN):: node_elem,dof_node

   REAL(ReKi)::Gup0(3),
   INTEGER(IntKi)::inode,temp_id,i

   hhx = 0.0D0
   hpx = 0.0D0
   CALL ShapeFunction1D(node_elem,rr,hhx,hpx)

   Gup0 = 0.0D0
   DO inode=1,node_elem
       temp_id = (inode-1)*dof_node
       DO i=1,3
           Gup0(i) = Gup0(i) + hpx(inode)*Nuu0(temp_id+i)
       ENDDO
   ENDDO

   jacobian = SQRT(DOT_PRODUCT(Gup0,Gup0))

   DO inode=1,node_elem
       hpx(inode) = hpx(inode)/jacobian
   ENDDO

   END SUBROUTINE BldComputeJacobian
