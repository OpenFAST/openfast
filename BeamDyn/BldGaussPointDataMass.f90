   SUBROUTINE BldGaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,Nm00,NmEta0,Nrho0,node_elem,dof_node,&
                                   &vvv,aaa,mmm,mEta,rho)

   REAL(ReKi),INTENT(IN):: hhx(:),hpx(:),Nvvv(:),Naaa(:)
   REAL(ReKi),INTENT(IN)::RR0(:,:)
   REAL(ReKi),INTENT(IN)::Nm00(:),NmEta0(:,:),Nrho0(:,:,:)
   INTEGER(IntKi),INTENT(IN):: node_elem,dof_node
   REAL(ReKi),INTENT(OUT)::mmm,mEta(:),rho(:,:)
   REAL(ReKi),INTENT(OUT):: vvv(:),aaa(:)
   
   REAL(ReKi):: hhi,hpi
   INTEGER(IntKi):: inode,temp_id,i,j

   mmm = 0.0D0
   mEta = 0.0D0
   rho = 0.0D0

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
       mmm = mmm + hhi * Nm00(inode)
       DO i=1,3
           mEta(i) = mEta(i) + hhi * NmEta0(i,inode)
           DO j=1,3
               rho(i,j) = rho(i,j) + hhi * Nrho0(i,j,inode)
           ENDDO
       ENDDO
   ENDDO

   mEta = MATMUL(RR0,mEta)
   rho = MATMUL(RR0,MATMUL(rho,TRANSPOSE(RR0)))

   END SUBROUTINE BldGaussPointDataMass
