   SUBROUTINE BldGaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,mEta0,rho0,node_elem,dof_node,vvv,aaa,mEta,rho)

   REAL(ReKi),INTENT(IN):: hhx(:),hpx(:),Nvvv(:),Naaa(:)
   REAL(ReKi),INTENT(IN)::RR0(:,:),mEta0(:),rho0(:,:)
   INTEGER(IntKi),INTENT(IN):: node_elem,dof_node
   REAL(ReKi),INTENT(OUT)::mEta(:),rho(:,:)
   REAL(ReKi),INTENT(OUT):: vvv(:),aaa(:)
   
   REAL(ReKi):: hhi,hpi
   INTEGER(IntKi):: inode,temp_id,i

   mEta = 0.0D0
   rho = 0.0D0
   mEta = MATMUL(RR0,mEta0)
   rho = MATMUL(RR0,MATMUL(rho0,TRANSPOSE(RR0)))

   vvv = 0.0D0
   aaa = 0.0D0
   DO inode=1,node_elem
       hhi = hhx(inode)
       hpi = hpx(inode)
       temp_id = (inode-1)*dof_node
       DO i=1,dof_node
           vvv(i) = vvv(i) + hhi * Nvvv(temp_id+i)
           aaa(i) = aaa(i) + hhi * Naaa(temp_id+i)
       ENDDO
   ENDDO

   END SUBROUTINE BldGaussPointDataMass
