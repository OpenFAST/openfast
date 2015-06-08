   SUBROUTINE UpdateDynamic(ainc,uf,vf,af,xf,coef,node_total,dof_node)

   REAL(ReKi), INTENT(IN):: ainc(:)
   REAL(DbKi),INTENT(IN)::coef(:)
   INTEGER(IntKi), INTENT(IN):: node_total, dof_node
   REAL(ReKi), INTENT(INOUT):: uf(:),vf(:),af(:),xf(:)

   REAL(ReKi):: rotf_temp(3), roti_temp(3), rot_temp(3)

   INTEGER(IntKi):: i, j, temp_id

   DO i=2, node_total
       temp_id = (i - 1) * dof_node
       rotf_temp = 0.0D0
       roti_temp = 0.0D0
       rot_temp = 0.0D0
       DO j = 1, 3
           uf(temp_id+j) = uf(temp_id+j) + coef(8) * ainc(temp_id+j)
           rotf_temp(j) = uf(temp_id+3+j)
           roti_temp(j) = coef(8) * ainc(temp_id+3+j)
       ENDDO
!       CALL CrvCompose_temp(rot_temp,roti_temp,rotf_temp,0)
       CALL CrvCompose(rot_temp,roti_temp,rotf_temp,0)
       DO j = 1, 3
           uf(temp_id+3+j) = rot_temp(j)
       ENDDO

       DO j=1, 6
           vf(temp_id+j) = vf(temp_id+j) + coef(7) * ainc(temp_id+j)
           af(temp_id+j) = af(temp_id+j) + ainc(temp_id+j)
           xf(temp_id+j) = xf(temp_id+j) + coef(9) * ainc(temp_id+j)
       ENDDO
   ENDDO

   END SUBROUTINE UpdateDynamic
