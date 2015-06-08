subroutine LineSearch(dof_total,xold,fold,gradient,rhs,dx,fmin,stpmax,check,&
                     &dof_node,uf,norder,hhp,wj,node_total,dmat,elem_total)
  
   implicit none
 
   integer dof_total,dof_node,norder,node_total,elem_total
   double precision uf(dof_total),hhp(norder+1,norder+1)
   double precision wj(norder+1),dmat(node_total,3)
   double precision xold(dof_total),fold,fmin,stpmax
   double precision gradient(dof_total),rhs(dof_total),dx(dof_total)
   logical check

   double precision alf, toly, tolerance
   parameter(alf=1.0D-4) ! A small number to indicate sufficient decrease of the function
   parameter(toly=1.0D-9) ! a small number to calculate the minimum step size

   double precision tmp, slope,alamin, alam,tmplam,rhs1,rhs2,f2,alam2,a,b,disc

   check=.FALSE.

   call Norm(dof_total,dx,tmp)
   
   IF(tmp>stpmax) dx=dx*stpmax/tmp  ! Scale if attempted step is too big

   slope=0.0d0
   slope=DOT_PRODUCT(gradient,dx)

   IF(slope>=0.0D0) THEN
      write(*,*) "roundoff problem in line search"
      RETURN
   ENDIF

   alamin=toly/MAXVAL(ABS(dx)/MAX(ABS(xold),1.0D0))  ! Compute lambda_min

   alam=1.0D0

   DO 

      !write(*,*) xold
      
      uf=xold+alam*dx
      
      CALL AssembleRHS(RHS,dof_node,dof_total,uf,&
                     &norder,hhp,wj,node_total,dmat,&
                     &elem_total)
                     

      fmin=0.5D0*DOT_PRODUCT(rhs,rhs)
 
      IF(alam<alamin) THEN  ! CONVERGENCE ON deltaX, the calling program should verify the convergence
         uf=xold
         check=.TRUE.
         RETURN
      ELSE IF(fmin<=fold+alf*alam*slope) THEN
            RETURN  ! SUFFICIENT function decrease
      ELSE 
         IF(ABS(alam-1.0D0)<TOLERANCE) THEN
               tmplam=-slope/(2.0D0*(fmin-fold-slope))
         ELSE  
               rhs1=fmin-fold-alam*slope
               rhs2=f2-fold- alam2*slope
               a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
               b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
       
               IF(ABS(a)<TOLERANCE)THEN
                  tmplam=-slope/(2.0D0*b)
               ELSE 
                  disc=b*b-3.0D0*a*slope
                  IF(disc<0.0d0) THEN
                     tmplam=0.5D0*alam  ! for imagine roots choose the max allowed
                  ELSE 
                     tmplam=(-b+sqrt(disc))/(3.0D0*a)
                  ENDIF
               ENDIF
               IF(tmplam>0.5d0*alam) tmplam=0.5d0*alam 
         ENDIF
      ENDIF
  
      alam2=alam
      f2=fmin
      alam=MAX(tmplam,0.1D0*alam)

   ENDDO
   
   return

END SUBROUTINE
