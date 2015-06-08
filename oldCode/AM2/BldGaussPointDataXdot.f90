   SUBROUTINE BldGaussPointDataXdot(hhx,Nvvv,node_elem,dof_node,vvv)
   !----------------------------------------------------------------------------------------
   ! This subroutine computes mass-related Gauss point values: 1) vvv, 2) mmm, 3) mEta,
   ! and 4) rho
   !----------------------------------------------------------------------------------------


   REAL(ReKi),    INTENT(IN   ):: hhx(:)        ! Shape function
   REAL(ReKi),    INTENT(IN   ):: Nvvv(:)       ! Element nodal velocity array
   INTEGER(IntKi),INTENT(IN   ):: node_elem     ! Number of node in one element
   INTEGER(IntKi),INTENT(IN   ):: dof_node      ! Number of DoF per node (=6)
   REAL(ReKi),    INTENT(  OUT):: vvv(:)        ! Velocitis at Gauss point (including linear and angular velocities)
   
   ! Local variables
   REAL(ReKi):: hhi
   INTEGER(IntKi):: inode,temp_id,i,j

   vvv = 0.0D0
   
   DO inode=1,node_elem
       hhi = hhx(inode)
       temp_id = (inode-1)*dof_node
       DO i=1,dof_node
           vvv(i) = vvv(i) + hhi * Nvvv(temp_id+i)
       ENDDO
   ENDDO

   END SUBROUTINE BldGaussPointDataXdot
