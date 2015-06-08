   SUBROUTINE BldGaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
   !----------------------------------------------------------------------------------------
   ! This subroutine computes initial Gauss point values: uu0, E10, and Stif
   !----------------------------------------------------------------------------------------
   REAL(ReKi),    INTENT(IN   ):: hhx(:)         ! Shape function
   REAL(ReKi),    INTENT(IN   ):: hpx(:)         ! Derivative of shape function
   REAL(ReKi),    INTENT(IN   ):: Nuu0(:)        ! Element initial nodal position array
   REAL(ReKi),    INTENT(IN   ):: Nrr0(:)        ! Element initial nodal relative rotation array
   INTEGER(IntKi),INTENT(IN   ):: node_elem      ! Number of node in one element
   INTEGER(IntKi),INTENT(IN   ):: dof_node       ! DoF per node (=6)
   REAL(ReKi),    INTENT(  OUT):: uu0(:)         ! Initial position array at Gauss point
   REAL(ReKi),    INTENT(  OUT):: E10(:)         ! E10 = x_0^\prime at Gauss point

   ! Local variables
   REAL(ReKi)                  :: hhi
   REAL(ReKi)                  :: hpi
   REAL(ReKi)                  :: rot0_temp(3)
   REAL(ReKi)                  :: rotu_temp(3)
   REAL(ReKi)                  :: rot_temp(3)
   INTEGER(IntKi)              :: inode
   INTEGER(IntKi)              :: temp_id
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   
   uu0 = 0.0D0
   E10 = 0.0D0
   DO inode=1,node_elem
       hhi = hhx(inode)
       hpi = hpx(inode)
       temp_id = (inode-1)*dof_node
       temp_id2 = (inode-1)*dof_node/2
       DO i=1,3
           uu0(i) = uu0(i) + hhi*Nuu0(temp_id+i)
           uu0(i+3) = uu0(i+3) + hhi*Nrr0(temp_id2+i)
           E10(i) = E10(i) + hpi*Nuu0(temp_id+i)
       ENDDO
   ENDDO   

   rot0_temp = 0.0D0
   rotu_temp = 0.0D0
   DO i=1,3
       rot0_temp(i) = Nuu0(i+3)
       rotu_temp(i) = uu0(i+3)
   ENDDO
   rot_temp = 0.0D0
   CALL CrvCompose_temp(rot_temp,rot0_temp,rotu_temp,0)
   DO i=1,3
       uu0(i+3) = rot_temp(i)
   ENDDO

   END SUBROUTINE BldGaussPointDataAt0
