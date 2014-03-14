   SUBROUTINE BDyn_RK4()
!
! This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!   Define constants k1, k2, k3, and k4 as 
!        k1 = dt * f(t        , x_t        )
!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!        k4 = dt * f(t + dt   , x_t + k3   ).
!   Then the continuous states at t = t + dt are
!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!
! For details, see:
! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!   Runge-Kutta." ß16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!   Cambridge University Press, pp. 704-716, 1992.
!
!..................................................................................................................................

      ! Initialize ErrStat


      ! interpolate u to find u_interp = u(t)

      ! find xdot at t
      CALL BDyn_CalcContStateDeriv(uuN0,uuN,vvN,Stif0,m00,mEta0,rho0,&
                                  &node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                                  &time,qdot,qddot)

      k1_qdot = deltaT * qdot
      k1_qddot = deltaT * qddot
  
      x_tmp_q    = uuN  + 0.5 * k1_qdot
      x_tmp_qdot = qdot + 0.5 * k1_qddot

      CALL QD2VEL(x_tmp_qdot,dof_node,node_total,vvN_tmp)

      ! find xdot at t + dt/2
      CALL BDyn_CalcContStateDeriv(uuN0,x_tmp_q,vvN,Stif0,m00,mEta0,rho0,&
                                  &node_elem,dof_node,elem_total,dof_total,node_total,ngp,&
                                  &time+deltaT/2.0D0,qdot,qddot)

      k2%q    = p%dt * xdot%q
      k2%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + 0.5 * k2%q
      x_tmp%dqdt = x%dqdt + 0.5 * k2%dqdt

      ! find xdot at t + dt/2
      CALL Mod1_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )
     
      k3%q    = p%dt * xdot%q
      k3%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + k3%q
      x_tmp%dqdt = x%dqdt + k3%dqdt

      ! interpolate u to find u_interp = u(t + dt)
      CALL Mod1_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL Mod1_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k4%q    = p%dt * xdot%q
      k4%dqdt = p%dt * xdot%dqdt

      x%q    = x%q    +  ( k1%q    + 2. * k2%q    + 2. * k3%q    + k4%q    ) / 6.      
      x%dqdt = x%dqdt +  ( k1%dqdt + 2. * k2%dqdt + 2. * k3%dqdt + k4%dqdt ) / 6.      

      CALL MeshDestroy ( u_interp%PointMesh       &
                       , ErrStat  = ErrStat         &
                       , ErrMess  = ErrMsg           )


END SUBROUTINE Mod1_RK4
