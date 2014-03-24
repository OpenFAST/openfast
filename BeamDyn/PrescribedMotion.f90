   SUBROUTINE PrescribedMotion(uuNf,vvNf,aaNf,time)

   REAL(DbKi),INTENT(IN):: time
   REAL(ReKi),INTENT(INOUT):: uuNf(:),vvNf(:),aaNf(:)

   INTEGER(IntKi):: i

   DO i=1,6
       uuNf(i) = 0.0D0
       vvNf(i) = 0.0D0
       aaNf(i) = 0.0D0
   ENDDO

!   uuNf(1) = 5.0D-02*time/3.0D0
!   vvNf(1) = 5.0D-02/3.0D0
   uuNf(5) = -4.0D0*TAN((3.1415926D0*time/3.0D0)/4.0D0)
   vvNf(5) = -3.1415926D0/3.0D0


   END SUBROUTINE PrescribedMotion
