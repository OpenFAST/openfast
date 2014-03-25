! This is the subroutine where the external force !
! is input.                                       !

  SUBROUTINE AppliedNodalLoad(F_ext,time,dof_total)

   REAL(ReKi),INTENT(INOUT)::F_ext(:) ! Applied nodal force
   REAL(DbKi),INTENT(IN)::time ! Current time
   INTEGER(IntKi),INTENT(IN)::dof_total ! Total degrees of freedom

!   INTEGER::i
!   WRITE(*,*) "time = ", time
   
   F_ext = 0.0D0

!   F_ext(dof_total-3) = 1.0D+02 * SIN(10.0D+00 * time) ! Input external force here

!   WRITE(*,*) F_ext(dof_total-3)
!   STOP

   END SUBROUTINE AppliedNodalLoad
