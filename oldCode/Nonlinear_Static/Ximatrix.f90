subroutine Ximatrix(Xim,dof_total,dof_node,norder,hhp,uf,nelem,nnode,&
                  &node_total,dmat,Jacobian)
   
   integer dof_total, dof_node,node_total
   integer norder,nelem,nnode 
   double precision uf(dof_total)
   double precision hhp(norder+1,norder+1),dmat(node_total,3), Jacobian  
   double precision Xim(4,4)
   
   double precision uprime,vprime,temp_theta,tempN(3,1)
   double precision Dm(3,3),Chim(3,1)
   integer m,temp_id
   
   
   call Dmatrix(Dm,node_total,dmat,nelem,nnode)
   call Chimatrix(Chim,dof_total,dof_node,norder,hhp,uf,nelem,nnode,Jacobian)
   
   tempN = MATMUL(Dm,Chim) 
   
   uprime = 0.0d0
   vprime = 0.0d0
   do m=1,norder+1
      temp_id = ((nelem-1)*norder+m-1)*dof_node
      uprime = uprime + hhp(m,nnode)*uf(temp_id+1)
      vprime = vprime + hhp(m,nnode)*uf(temp_id+2)
   enddo
   
   uprime = uprime/Jacobian
   vprime = vprime/Jacobian
   
   Xim = 0.0d0
   
   temp_id = (((nelem-1)*norder+nnode)-1)*dof_node
   temp_theta = uf(temp_id+dof_node)
   
   Xim(1,4) = -COS(temp_theta)*tempN(2,1)-SIN(temp_theta)*tempN(1,1)
   Xim(2,4) = -SIN(temp_theta)*tempN(2,1) + COS(temp_theta)*tempN(1,1)
   Xim(4,1) = Xim(1,4)
   Xim(4,2) = Xim(2,4)
   Xim(4,4) = -(COS(temp_theta)*uprime+SIN(temp_theta)*vprime+COS(temp_theta))*tempN(1,1) -&
            &(-SIN(temp_theta)*uprime+COS(temp_theta)*vprime - SIN(temp_theta))*tempN(2,1)
        
!   write(*,*) "Xim(1,4)",Xim(4,4)
            
   return
end subroutine
