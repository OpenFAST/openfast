   SUBROUTINE NodalDataMass(RR0,mEta0,rho0,mEta,rho)

   REAL(ReKi),INTENT(IN)::RR0(:,:),mEta0(:),rho0(:,:)
   REAL(ReKi),INTENT(INOUT)::mEta(:),rho(:,:)

   mEta = MATMUL(RR0,mEta0)
   rho = MATMUL(RR0,MATMUL(rho0,TRANSPOSE(RR0)))   

   END SUBROUTINE NodalDataMass 
