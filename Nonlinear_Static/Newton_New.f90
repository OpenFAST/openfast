subroutine Newton_New(dof_node,dof_total,norder,node_total,elem_total,&
                     &hhp,uf,dmat,wj,niter,Jacobian)
   
   integer dof_node,dof_total,norder,node_total,elem_total,niter
   double precision hhp(norder+1,norder+1),wj(norder+1),dmat(node_total,3)
   double precision uf(dof_total),Jacobian
   
   double precision RHS(dof_total), KT(dof_total,dof_total)
   double precision ui(dof_total), bc(dof_total)

   double precision ui_old(dof_total)!mas ????
   double precision rel_change(dof_total)

   parameter(TOLF  = 1.0d-6)  ! Xiao & Zhong use 1d-5
   integer i,j,m
   double precision temp1,temp2,FmL
   
   ui=0.0d0
   FmL=3.14d0 
     
   do i=1, niter
   
       write(*,*) "niter=",i
       
       call AssembleRHS(RHS, dof_node, dof_total, uf, &
                  & norder, hhp, wj, node_total, dmat, &
                  & elem_total)
       if(i==1) RHS(dof_total) = RHS(dof_total) + FmL
       call AssembleKT(KT, dof_node, dof_total, norder, node_total, elem_total,&
                    & hhp, uf, dmat, wj, Jacobian)
                    
       do j=1,dof_total
           write(*,*) KT(2,j)
       enddo
       
       if(i.gt.1) stop
                    
       ui_old = ui
       
      ! Solve the linear system
      !====================================
       bc = 1.0d0

       bc(1) = 0.0d0
       bc(2) = 0.0d0
       bc(3) = 0.0d0  
                    
       
       call CGSolver(RHS, KT, ui, bc, dof_total)
              
       uf = uf + ui
       
       
!       write(*,*) "RHS"
!       do j=1,dof_total
!           write(*,*) RHS(j)
!       enddo
!       stop
       
       rel_change = ABS(ui-ui_old)
       
       call Norm(dof_total, uf, temp1)
       call Norm(dof_total, rel_change, temp2) 
       
       errx = temp2 - tolx*temp1
       
!       write(*,*) "niter=",i
!       write(*,*) "Norm"
!       write(*,*) temp1
!       write(*,*) temp2
!       write(*,*) errx
       
       if(errx .le. 0) return
       if(i==niter) then
           write(*,*) "The solution does not converge after the maximum number of iterations" 
           stop
       endif
                   
       
   enddo
     
   return
end subroutine
