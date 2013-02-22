subroutine Amatrix(Am,dof_total,dof_node,norder,hhp,uf,nelem,nnode)

   integer dof_total, dof_node
   integer norder,nelem,nnode      
   double precision uf(dof_total)
   double precision hhp(norder+1,norder+1)      
   double precision Am(3,4)
   
   double precision uprime,vprime,temp_theta
   integer m,temp_id


   uprime = 0.0d0
   vprime = 0.0d0
   do m=1,norder+1
      temp_id = ((nelem-1)*norder+m-1)*dof_node
      uprime = uprime + hhp(m,nnode)*uf(temp_id+1)
      vprime = vprime + hhp(m,nnode)*uf(temp_id+2)
   enddo

   Am=0.0d0

   temp_id = (((nelem-1)*norder+nnode)-1)*dof_node
   temp_theta = uf(temp_id+dof_node)

   Am(1,1) = COS(temp_theta)
   Am(1,2) = SIN(temp_theta)
   Am(1,4) = -SIN(temp_theta) * uprime + COS(temp_theta) * vprime

   Am(2,1) = -SIN(temp_theta)
   Am(2,2) = COS(temp_theta)
   Am(2,4) = -COS(temp_theta) * uprime - SIN(temp_theta) * vprime

   Am(3,3) = 1.0d0

   return

end subroutine
      
      
      
      
      
