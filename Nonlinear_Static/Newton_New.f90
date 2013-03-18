subroutine Newton_New(dof_node,dof_total,norder,node_total,elem_total,&
                     &hhp,uf,dmat,wj,niter,Jacobian, xloc, F_ext)
   
   integer dof_node,dof_total,norder,node_total,elem_total,niter
   double precision hhp(norder+1,norder+1),wj(norder+1),dmat(node_total,3)
   double precision uf(dof_total),Jacobian
   
   double precision RHS(dof_total), KT(dof_total,dof_total)
   double precision ui(dof_total), bc(dof_total)

   double precision ui_old(dof_total)!mas ????
   double precision rel_change(dof_total)

   double precision xloc(node_total)
   double precision F_ext
   
   double precision position(dof_total)

   parameter(TOLF  = 1.0d-5)  ! Xiao & Zhong use 1d-5
   integer i,j,m
   double precision temp1,temp2,FmL,errf,errx

   character*4 filename

   external filename
   
   ui=0.0d0

   yloc = 0.

   position = 0.

!   FmL=3.14d-5 
   do i = 1, node_total
      position(3*(i-1)+1) = xloc(i)
   enddo
   
   call outputvtk(position + uf, node_total, dof_total, 'initial' )         
     
   do i=1, niter
   
       write(*,*) "starting Newton iteration ",i
              
       call AssembleRHS(RHS, dof_node, dof_total, uf, &
                  & norder, hhp, wj, node_total, dmat, &
                  & elem_total, Jacobian, F_ext)
                  
       
!       if(i==1) RHS(dof_total) = RHS(dof_total) + FmL
       call AssembleKT(KT, dof_node, dof_total, norder, node_total, elem_total,&
                    & hhp, uf, dmat, wj, Jacobian)
                    
!       write(*,*) "KT"
!       do j=1,dof_total
!           write(*,*) KT(1,j)
!       enddo
!       if(i.gt.1) stop
                    
       errf = 0.0d0
       
       call Norm(dof_total, RHS, errf) 
       
!       write(*,*) "errf", errf
       
       if(errf .le. TOLF) return
                    
       ui_old = ui
       
      ! Solve the linear system
      !====================================
       bc = 1.0d0

       bc(1) = 0.0d0
       bc(2) = 0.0d0
       bc(3) = 0.0d0  
                    
       
       call CGSolver(RHS, KT, ui, bc, dof_total)
              
       
       
       
!       write(*,*) "ui"
!       do j=1,dof_total
!           write(*,*) ui(j)
!       enddo
!       if(i.gt.4) stop
       
       rel_change = ABS(ui-ui_old)

!       do j=1,3
!           RHS(i)=0.0d0
!       enddo
       
       call Norm(dof_total, uf, temp1)
       call Norm(dof_total, rel_change, temp2) 
       
       errx = temp2 !- tolx*temp1

!        call Norm(dof_total, RHS, errx)

        write(*,*) "errx",errx
 
!        if (i .eq. 10) stop
   
        if(errx .lt. TOLF*temp1) return

!        if(errx .le. TOLF) return
        
        uf = uf + ui
        
        !position = position + uf
 
       call outputvtk(position + uf, node_total, dof_total, 'output'//filename(i) )         
               
!       write(*,*) "Norm"
!       write(*,*) temp1
!       write(*,*) temp2
!       write(*,*) errx
       
!       if(i.gt.1) stop
       
!       if(errx .le. 0) return
       if(i==niter) then
           write(*,*) "The solution does not converge after the maximum number of iterations" 
           stop
       endif
                   
       
   enddo
     
   return
end subroutine
