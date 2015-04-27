   SUBROUTINE BD_ComputeIniCoef(kp_member,kp_coord,Coef)

   REAL(ReKi),    INTENT(IN   ):: kp_coord(:,:)
   INTEGER(IntKi),INTENT(IN   ):: kp_member
   REAL(ReKi),    INTENT(  OUT):: Coef(:,:,:)

   REAL(ReKi)    :: K(4*(kp_member-1),4*(kp_member-1))
   REAL(ReKi)    :: RHS(4*(kp_member-1))
   REAL(ReKi)    :: sol_temp(4*(kp_member-1))
   REAL(ReKi)    :: d
   REAL(ReKi)    :: temp
   INTEGER(IntKi):: indx(4*(kp_member-1))
   INTEGER(IntKi):: i
   INTEGER(IntKi):: j
   INTEGER(IntKi):: m
   INTEGER(IntKi):: temp_id1
   
   DO i=1,4
       K(:,:) = 0.0D0
       RHS(:) = 0.0D0

       K(1,3) = 2.0D0
       RHS(1) = 0.0D0
       DO j=1,kp_member-2
           temp_id1 = (j-1)*4
           K(temp_id1+2,temp_id1+1) = 1.0D0
           K(temp_id1+3,temp_id1+1:temp_id1+4) = 1.0D0
           K(temp_id1+4,temp_id1+2) = 1.0D0
           K(temp_id1+4,temp_id1+3) = 2.0D0
           K(temp_id1+4,temp_id1+4) = 3.0D0
           K(temp_id1+4,temp_id1+6) = -1.0D0
           K(temp_id1+5,temp_id1+3) = 2.0D0
           K(temp_id1+5,temp_id1+4) = 6.0D0
           K(temp_id1+5,temp_id1+7) = -2.0D0
           RHS(temp_id1+2) = kp_coord(j,i)
           RHS(temp_id1+3) = kp_coord(j+1,i)
       ENDDO
       temp_id1 = (kp_member-2)*4
       K(temp_id1+2,temp_id1+1) = 1.0D0
       K(temp_id1+3,temp_id1+1:temp_id1+4) = 1.0D0
       K(temp_id1+4,temp_id1+3) = 2.0D0
       K(temp_id1+4,temp_id1+4) = 6.0D0
       RHS(temp_id1+2) = kp_coord(kp_member-1,i)
       RHS(temp_id1+3) = kp_coord(kp_member,i)
       RHS(temp_id1+4) = 0.0D0
       CALL ludcmp(K,4*(kp_member-1),indx,d)
       CALL lubksb(K,4*(kp_member-1),indx,RHS,sol_temp)
       DO j=1,kp_member-1
           DO m=1,4
               temp_id1 = (j-1)*4+m
               Coef(j,m,i) = sol_temp(temp_id1)
           ENDDO
       ENDDO
   ENDDO


   END SUBROUTINE BD_ComputeIniCoef
