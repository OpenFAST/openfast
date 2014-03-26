   SUBROUTINE BldComputeJacobianLSGL(rr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,jacobian)
   REAL(ReKi),INTENT(IN):: rr ! Counter, gpr in ElementMatrix function
   REAL(ReKi),INTENT(IN):: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),INTENT(IN):: gp(:)
   REAL(ReKi),INTENT(IN):: GLL_temp(:)
   REAL(ReKi),INTENT(OUT):: jacobian
   REAL(ReKi),INTENT(OUT):: hhx(:)
   REAL(ReKi),INTENT(OUT):: hpx(:)
   INTEGER(IntKi),INTENT(IN):: node_elem
   INTEGER(IntKi),INTENT(IN):: dof_node
   INTEGER(IntKi),INTENT(IN):: ngp
   INTEGER(IntKi),INTENT(IN):: igp

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
   
   !WRITE(*,*) "Jacobian = ", jacobian

   DO inode=1,node_elem
       hpx(inode) = hpx(inode)/jacobian
   ENDDO

   END SUBROUTINE BldComputeJacobianLSGL
