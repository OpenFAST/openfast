   SUBROUTINE ComputeStrainForceGP(Nuu0,Nuuu,Nrr0,Nrrr,NStif0,ngp,node_elem,dof_node,eeeGP,fffGP)

   REAL(ReKi),INTENT(IN):: Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:)
   REAL(ReKi),INTENT(IN):: NStif0(:,:,:)
   INTEGER(IntKi),INTENT(IN):: ngp,node_elem,dof_node

   REAL(ReKi),INTENT(OUT):: eeeGP(:),fffGP(:)      

   REAL(ReKi),ALLOCATABLE:: gp(:),gw(:),hhx(:),hpx(:)
   
   REAL(ReKi),ALLOCATABLE:: GLL_temp(:),w_temp(:)
   
   REAL(ReKi):: uu0(6),E10(3),RR0(3,3),kapa(3),E1(3),Stif(6,6),cet
   REAL(ReKi):: uuu(6),uup(3),Jacobian,gpr
   REAL(ReKi):: eee(6),fff(6),temp_eee(6)
   REAL(ReKi):: tempS(3),tempK(3)
   REAL(ReKi):: e1s,k1s

   INTEGER(IntKi)::igp,i,j,m,n,temp_id
   INTEGER(IntKi)::allo_stat


   ALLOCATE(gp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gp = 0.0D0

   ALLOCATE(gw(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gw = 0.0D0

   ALLOCATE(hhx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hhx = 0.0D0

   ALLOCATE(hpx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hpx = 0.0D0
   
   ALLOCATE(GLL_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   GLL_temp = 0.0D0
   
   ALLOCATE(w_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   w_temp = 0.0D0
   

   CALL BeamDyn_gen_gll_LSGL(node_elem-1,GLL_temp,w_temp)
!   CALL BldSet1DGaussPointScheme(ngp,gp,gw)
   CALL BldGaussPointWeight(ngp,gp,gw)
   DO igp=1,ngp
       gpr = gp(igp)
       CALL BldComputeJacobianLSGL(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)
       CALL BldGaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,NStif0,node_elem,dof_node,uu0,E10,Stif)
       Stif(:,:) = NStif0(:,:,igp)
       CALL BldGaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)
       eee = 0.0D0
       tempS = 0.0D0
       tempK = 0.0D0
       DO i=1,3
           eee(i) = E1(i) - RR0(i,1)
           eee(i+3) = kapa(i)

           tempS(i) = eee(i)
           tempK(i) = eee(i+3)
       ENDDO
       temp_id = (igp-1)*dof_node
       DO i=1,dof_node
          eeeGP(temp_id+i) = eee(i)
       ENDDO
      
       tempS = MATMUL(TRANSPOSE(RR0),tempS)
       tempK = MATMUL(TRANSPOSE(RR0),tempK)
       
       temp_eee = 0.0D0
       DO i=1,3
           temp_eee(i) = tempS(i)
           temp_eee(i+3) = tempK(i)
       ENDDO
       
       fff = 0.0D0
       fff = MATMUL(NStif0(:,:,igp),temp_eee)

       k1s = temp_eee(4)
       fff(1) = fff(1) + 0.5D0*cet*k1s*k1s
       e1s = temp_eee(1)
       fff(4) = fff(4) + cet*e1s*k1s
       
       DO i=1,3
           tempS(i) = fff(i)
           tempK(i) = fff(i+3)
       ENDDO       
       tempS = MATMUL(RR0,tempS)
       tempK = MATMUL(RR0,tempK)
       DO i=1,3
           fffGP(temp_id+i) = tempS(i)
           fffGP(temp_id+i+3) = tempK(i)
       ENDDO

   ENDDO

   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
        ENDIF

   END SUBROUTINE ComputeStrainForceGP
