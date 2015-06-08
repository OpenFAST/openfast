program driver
   
   implicit double precision (a-h,o-z)
   integer norder, elem_total, node_total    
   integer dof_node,dof_total,nquar
   integer niter
   integer nstep
 
   parameter (norder = 20)
   parameter (elem_total = 1 )
   parameter (node_total = elem_total*norder + 1)
   parameter (dof_node = 3)
   parameter (dof_total = dof_node*node_total)
     
   double precision  xj(norder+1), wj(norder+1)
   double precision  dloc(node_total),dmat(node_total,3)
   double precision  hhp(norder+1,norder+1)
   double precision  uf(dof_total),Jacobian
   double precision  uf_temp(dof_total)
   
   double precision  h0,h1,b0,b1
   double precision  FmL,FmL_temp
   
   integer i,j,temp_id
 
   open( unit = 20, file = 'node_loc.dat', status = 'unknown')
   open( unit = 25, file = 'node_mat.dat', status = 'unknown')
   open( unit = 30, file = 'u.dat', status = 'unknown')
   open( unit = 35, file = 'v.dat', status = 'unknown')      
   open( unit = 40, file = 'rot.dat', status = 'unknown')      
   
   call gen_gll(norder, xj, wj)
   
    
!  Define material constants
   Young = 180.0d+9  !1.95d11       !30.0d+06          !      !Young Modulus
   shear = 5.0/6.0  !5.d0/6.d0    !Shear Correction Factor
   poisson = 0.3d0     !0.0d0      !Poisson's Ratio
   G1 = 0.5d0 * shear * Young / ( 1. + poisson )    !Shear Modulus with Correction Factor
   rho = 7700.d0        !Density
    
!   write(*,*) 'G1 = ', G1
            
!  Define geometric parameters
   xmin = 0.0d0
   xmax = 10.0d+00

   blength = xmax - xmin              !Beam length

   base = 1.
   h_over_l = 1.d-1

   elem_length = blength / float(elem_total)    !Element size

   Jacobian = elem_length / 2.d0                   !Determinant of jacobian
   pi = acos(-1.d0)
   
   FmL = (5.833d+06)  !0.3   !2*pi*10 !(5.833d+08)*2*pi   !  !1.75d+07

   FmL = 1.  
   
   niter = 50   !100000
   nstep = 1 
   

   call gen_deriv(xj, hhp, norder+1)
   

   call NodeLoc(dloc, nquar, xmin, elem_length,& 
                & xj, norder, elem_total, node_total, blength)
    
   write(*,*) "nquar",nquar
   write(*,*) "nquar location",dloc(nquar)
    
   do i=1, node_total
      write(20,*) i, dloc(i)
   enddo
   
!   do i=1, node_total
!      write(*,*) "xj"
!      write(*,*) i, xj(i)
!   enddo
      
            
! Obtain EI and kGA at each node, which is also quadrature point for GLL

   h0 = blength * h_over_l
   h1 = blength * h_over_l
   b0 = base
   b1 = base
   
   dmat = 0.0d0
    
   call NodeMat(dmat, dloc, h0, h1, b0, b1, node_total, blength,&
              & Young, G1) 
    
!   write(*,*) "dmat"
!   do i=1,node_total
!       write(*,*) dmat(i,1)
!       write(*,*) dmat(i,2)
!       write(*,*) dmat(i,3)
!   enddo
!   stop
    
   do i=1, node_total
      write(25,*) i, dmat(i,1),dmat (i,2),dmat(i,3)
   enddo
      
   uf = 0.0d0
   uf_temp = 0.0d0

   write(*,*) "entering NewtonRaphson"
   write(*,*) "dof_node = ", dof_node
   write(*,*) "dof_total = ", dof_total
   write(*,*) "node_total = ", node_total
   write(*,*) "elem_total = ", elem_total


!   call NewtonRaphson(dof_node, dof_total, norder, node_total, elem_total,&
!                    & hhp, uf, dmat, wj, niter, Jacobian)

    delta_FmL = FmL / float(nstep)

    do i=1, nstep

       FmL_temp = float(i) * delta_FmL

       write(*,*) "Incremental Load",i

       write(*,*) "FmL_temp",FmL_temp

       call Newton_New(dof_node, dof_total, norder, node_total, elem_total,&
                        & hhp, uf_temp, dmat, wj, niter, Jacobian, dloc, FmL_temp)

        uf = uf + uf_temp
        uf_temp = uf
       
       write(*,*) "Finished NR"
       write(*,*) uf_temp

    enddo
                     
         
   
     
   do i=1, node_total    
          
      temp_id = (i-1)*dof_node
          
      write(30,*) dloc(i), uf_temp(temp_id+1)
      write(35,*) dloc(i), uf_temp(temp_id+2)
      write(40,*) dloc(i), uf_temp(temp_id+3)
   enddo
      
end program
      
      
      
      
