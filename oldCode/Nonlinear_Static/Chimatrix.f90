subroutine Chimatrix(Chim,dof_total,dof_node,norder,hhp,uf,nelem,nnode,Jacobian)
   
   integer dof_total, dof_node
   integer norder,nelem,nnode 
   double precision uf(dof_total)
   double precision hhp(norder+1,norder+1),Jacobian   
   double precision Chim(3,1) 
   
   double precision uprime,vprime,theta_prime,temp_theta
   integer m,temp_id
   
   
   uprime = 0.0d0
   vprime = 0.0d0
   theta_prime = 0.0d0
   
   do m=1,norder+1
      temp_id = ((nelem-1)*norder+m-1)*dof_node
      uprime = uprime + hhp(m,nnode)*uf(temp_id+1)
      vprime = vprime + hhp(m,nnode)*uf(temp_id+2)
      theta_prime = theta_prime + hhp(m,nnode)*uf(temp_id+3)
   enddo
   
   uprime = uprime/Jacobian
   vprime = vprime/Jacobian
   theta_prime = theta_prime/Jacobian
   
   Chim=0.0d0
   
   temp_id = (((nelem-1)*norder+nnode)-1)*dof_node
   temp_theta = uf(temp_id+dof_node)
   
   Chim(1,1) = COS(temp_theta)*uprime+SIN(temp_theta)*vprime-1.0d0 + COS(temp_theta)
   Chim(2,1) = -SIN(temp_theta)*uprime + COS(temp_theta)*vprime - SIN(temp_theta)
   Chim(3,1) = theta_prime
   
   return
   
end subroutine

