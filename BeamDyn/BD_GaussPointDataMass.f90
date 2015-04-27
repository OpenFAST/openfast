   SUBROUTINE BD_GaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,node_elem,dof_node,&
                                    vvv,aaa,vvp,mmm,mEta,rho)

   REAL(ReKi),     INTENT(IN   ):: hhx(:)
   REAL(ReKi),     INTENT(IN   ):: hpx(:)
   REAL(ReKi),     INTENT(IN   ):: Nvvv(:)
   REAL(ReKi),     INTENT(IN   ):: Naaa(:)
   REAL(ReKi),     INTENT(IN   ):: RR0(:,:)
   INTEGER(IntKi), INTENT(IN   ):: node_elem
   INTEGER(IntKi), INTENT(IN   ):: dof_node
   REAL(ReKi),     INTENT(  OUT):: vvv(:)
   REAL(ReKi),     INTENT(  OUT):: vvp(:)
   REAL(ReKi),     INTENT(  OUT):: aaa(:)
   REAL(ReKi),     INTENT(INOUT):: mmm
   REAL(ReKi),     INTENT(INOUT):: mEta(:)
   REAL(ReKi),     INTENT(INOUT):: rho(:,:)
   
   REAL(ReKi)                   :: hhi
   REAL(ReKi)                   :: hpi
   INTEGER(IntKi)               :: inode
   INTEGER(IntKi)               :: temp_id
   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j

   vvv(:) = 0.0D0
   vvp(:) = 0.0D0
   aaa(:) = 0.0D0

   DO inode=1,node_elem
       hhi = hhx(inode)
       hpi = hpx(inode)
       temp_id = (inode-1)*dof_node
       DO i=1,dof_node
           vvv(i) = vvv(i) + hhi * Nvvv(temp_id+i)
           vvp(i) = vvp(i) + hpi * Nvvv(temp_id+i)
           aaa(i) = aaa(i) + hhi * Naaa(temp_id+i)
       ENDDO
   ENDDO

   mEta = MATMUL(RR0,mEta)
   rho = MATMUL(RR0,MATMUL(rho,TRANSPOSE(RR0)))

   END SUBROUTINE BD_GaussPointDataMass
