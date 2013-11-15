   SUBROUTINE BldGaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)

   REAL(ReKi),INTENT(IN):: hhx(:),hpx(:),Nuu0(:),Nrr0(:)
   REAL(ReKi),INTENT(INOUT):: uu0(:),E10(:)
   INTEGER(IntKi),INTENT(IN):: node_elem,dof_node,  

   REAL(ReKi):: hhi,hpi,rot0_temp(3),rotu_temp(3),rot_temp(3)
   INTEGER(IntKi):: inode,temp_id,temp_id2,i
   
   uu0 = 0.0D0
   E10 = 0.0D0
   DO inode=1,node_elem
       hhi = hhx(inode)
       hpi = hpx(inode)
       temp_id = (inode-1)*dof_node
       temp_id2 = (inode-1)*dof_node/2
       DO i=1,3
           uu0(i) = uu0(+i) + hhi*Nuu0(temp_id+i)
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
   CALL CrvCompose(rot_temp,rot0_temp,rotu_temp,0)
   DO i=1,3
       uu0(i+3) = rot_temp(i)
   ENDDO

   END SUBROUTINE BldGaussPointDataAt0
