
   SUBROUTINE AppliedNodalLoad(F_ext,time,dof_total)
   !***************************************************************************************
   ! This is the subroutine where the external force is input.
   !***************************************************************************************

   REAL(ReKi),    INTENT(INOUT)::F_ext(:)    ! Applied nodal force
   REAL(DbKi),    INTENT(IN)::   time        ! Current time
   INTEGER(IntKi),INTENT(IN)::   dof_total   ! Total degrees of freedom

   
   F_ext = 0.0D0

!   F_ext(dof_total-3) = 1.0D+02 * SIN(10.0D+00 * time) ! Input external force here

   END SUBROUTINE AppliedNodalLoad
