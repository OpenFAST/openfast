   SUBROUTINE AppliedNodalLoad(F_ext,time,dof_total)

   REAL(ReKi),INTENT(INOUT)::F_ext(:)
   REAL(ReKi),INTENT(IN)::time
   INTEGER(IntKi),INTENT(IN)::dof_total

!   INTEGER::i

   F_ext = 0.0D0

   F_ext(dof_total-3) = 1.0D+05 * SIN(20 * time)



   END SUBROUTINE AppliedNodalLoad
