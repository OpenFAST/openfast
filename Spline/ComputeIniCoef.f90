   SUBROUTINE ComputeIniCoef(kp_member,kp_coord,Coef)

   REAL(ReKi),    INTENT(IN   ):: kp_coord(:,:)
   INTEGER(IntKi),INTENT(IN   ):: kp_member
   REAL(ReKi),    INTENT(  OUT):: Coef(:,:,:)

   REAL(ReKi)    :: K(4*member_total,4*member_total)
   REAL(ReKi)    :: RHS(4*member_total)
   REAL(ReKi)    :: sol_temp(4*member_total)
   REAL(ReKi)    :: d
   REAL(ReKi)    :: temp
   INTEGER(IntKi):: indx(4*member_total)
   INTEGER(IntKi):: i
   INTEGER(IntKi):: j
   INTEGER(IntKi):: m
   INTEGER(IntKi):: temp_id1
   
   DO i=1,3
       K(:,:) = 0.0D0
       RHS(:) = 0.0D0

       K(1,3) = 2.0D0
       RHS(1) = 0.0D0!tang_vec(1,i) !*(kp_coord(2,i)-kp_coord(1,i))
       IF(member_total==1) THEN
           K(2,1) = 1.0D0
           K(3,1:4) = 1.0D0
!           K(4,2) = 1.0D0
!           K(4,3) = 2.0D0
!           K(4,4) = 3.0D0
           K(4,3) = 2.0D0
           K(4,4) = 6.0D0
           RHS(2) = kp_coord(1,i)
           RHS(3) = kp_coord(2,i)
           RHS(4) = 0.0D0 !tang_vec(2,i)!*(kp_coord(2,i)-kp_coord(1,i))
       ELSE
           DO j=1,member_total-1
               temp_id1 = (j-1)*4
               K(temp_id1+2,temp_id1+1) = 1.0D0
               K(temp_id1+3,temp_id1+1:temp_id1+4) = 1.0D0
               K(temp_id1+4,temp_id1+2) = 1.0D0
               K(temp_id1+4,temp_id1+3) = 2.0D0
               K(temp_id1+4,temp_id1+4) = 3.0D0
               temp = ((kp_coord(j+1,i)-kp_coord(j,i))/(kp_coord(j+2,i)-kp_coord(j+1,i)))
               K(temp_id1+4,temp_id1+6) = -temp
               K(temp_id1+5,temp_id1+3) = 2.0D0
               K(temp_id1+5,temp_id1+4) = 6.0D0
               K(temp_id1+5,temp_id1+7) = -2.0D0*temp*temp
               RHS(temp_id1+2) = kp_coord(j,i)
               RHS(temp_id1+3) = kp_coord(j+1,i)
           ENDDO
           temp_id1 = (member_total-1)*4
           K(temp_id1+2,temp_id1+1) = 1.0D0
           K(temp_id1+3,temp_id1+1:temp_id1+4) = 1.0D0
!           K(temp_id1+4,temp_id1+2) = 1.0D0
           K(temp_id1+4,temp_id1+3) = 2.0D0
           K(temp_id1+4,temp_id1+4) = 6.0D0
           RHS(temp_id1+2) = kp_coord(member_total,i)
           RHS(temp_id1+3) = kp_coord(member_total+1,i)
           RHS(temp_id1+4) = 0.0D0 ! tang_vec(2,i)*(kp_coord(member_total+1,i)-kp_coord(member_total,i))
       ENDIF
       CALL ludcmp(K,4*member_total,indx,d)
       CALL lubksb(K,4*member_total,indx,RHS,sol_temp)
       DO j=1,member_total
           DO m=1,4
               temp_id1 = (j-1)*4+m
               Coef(j,m,i) = sol_temp(temp_id1)
           ENDDO
       ENDDO
   ENDDO


   END SUBROUTINE ComputeIniCoef
