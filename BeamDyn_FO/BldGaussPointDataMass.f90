   SUBROUTINE BldGaussPointDataMass(hhx,hpx,Nvvv,RR0,Nm00,NmEta0,Nrho0,node_elem,dof_node,&
                                   &vvv,mmm,mEta,rho)
   !----------------------------------------------------------------------------------------
   ! This subroutine computes mass-related Gauss point values: 1) vvv, 2) mmm, 3) mEta,
   ! and 4) rho
   !----------------------------------------------------------------------------------------


   REAL(ReKi),    INTENT(IN)::  hhx(:)        ! Shape function
   REAL(ReKi),    INTENT(IN)::  hpx(:)        ! Derivative of shape function
   REAL(ReKi),    INTENT(IN)::  Nvvv(:)       ! Element nodal velocity array
   REAL(ReKi),    INTENT(IN)::  RR0(:,:)      ! Rotation matrix at Gauss point
   REAL(ReKi),    INTENT(IN)::  Nm00(:)       ! Element initial nodal mass density array
   REAL(ReKi),    INTENT(IN)::  NmEta0(:,:)   ! Element initial nodal m\Eta array 
   REAL(ReKi),    INTENT(IN)::  Nrho0(:,:,:)  ! Element initial nodal sectional tensor of inertia per unit span
   INTEGER(IntKi),INTENT(IN)::  node_elem     ! Number of node in one element
   INTEGER(IntKi),INTENT(IN)::  dof_node      ! Number of DoF per node (=6)
   REAL(ReKi),    INTENT(OUT):: mmm           ! Mass density at Gauss point
   REAL(ReKi),    INTENT(OUT):: mEta(:)       ! m\Eta resolved in inertia frame at Gauss point
   REAL(ReKi),    INTENT(OUT):: rho(:,:)      ! Tensor of inertia resolved in inertia frame at Gauss point
   REAL(ReKi),    INTENT(OUT):: vvv(:)        ! Velocitis at Gauss point (including linear and angular velocities)
   
   ! Local variables
   REAL(ReKi):: hhi,hpi
   INTEGER(IntKi):: inode,temp_id,i,j

   mmm = 0.0D0
   mEta = 0.0D0
   rho = 0.0D0

   vvv = 0.0D0
   
   DO inode=1,node_elem
       hhi = hhx(inode)
       hpi = hpx(inode)
       temp_id = (inode-1)*dof_node
       DO i=1,dof_node
           vvv(i) = vvv(i) + hhi * Nvvv(temp_id+i)
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
