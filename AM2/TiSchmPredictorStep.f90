   SUBROUTINE TiSchmPredictorStep(uuNi,vvNi,aaNi,xxNi,coef,deltat,uuNf,vvNf,aaNf,xxNf,uuNf_inc,node_total,dof_node)

   REAL(ReKi),INTENT(IN)::uuNi(:),vvNi(:),aaNi(:),xxNi(:)
   REAL(DbKi),INTENT(IN)::deltat,coef(:)
   REAL(ReKi),INTENT(INOUT)::uuNf(:),vvNf(:),aaNf(:),xxNf(:),uuNf_inc(:)

   INTEGER(IntKi),INTENT(IN)::node_total,dof_node

   REAL(ReKi)::vi,ai,xi,tr(6),tr_temp(3),uuNi_temp(3),rot_temp(3)
   INTEGER::i,j,temp_id

   DO i=1,node_total

       DO j=1,6
           temp_id = (i - 1) * dof_node + j
           vi = vvNi(temp_id)
           ai = aaNi(temp_id)
           xi = xxNi(temp_id)
           tr(j) = deltat * vi + coef(1) * ai + coef(2) * xi
           vvNf(temp_id) = vi + coef(3) * ai + coef(4) * xi
           aaNf(temp_id) = 0.0D0
           xxNf(temp_id) = coef(5) * ai + coef(6) * xi
       ENDDO

       tr_temp = 0.0D0
       uuNi_temp = 0.0D0
       DO j=1,3
           temp_id = (i - 1) * dof_node + j
           uuNf(temp_id) = uuNi(temp_id) + tr(j)
           uuNf_inc(temp_id) = tr(j)
           uuNf_inc(temp_id+3) = tr(j+3)
           tr_temp(j) = tr(j+3)
           uuNi_temp(j) = uuNi(temp_id + 3)
       ENDDO
       rot_temp = 0.0D0
!       CALL CrvCompose_temp(rot_temp,tr_temp,uuNi_temp,0)
       CALL CrvCompose(rot_temp,tr_temp,uuNi_temp,0)
       DO j=1,3
           temp_id = (i - 1) * dof_node +j
           uuNf(temp_id + 3) = rot_temp(j)
       ENDDO

   ENDDO

   END SUBROUTINE TiSchmPredictorStep
