   SUBROUTINE ComputeUDN(node_total,dof_node,vvN,uuN,udN)

   REAL(ReKi),INTENT(IN):: vvN(:)
   REAL(ReKi),INTENT(IN):: uuN(:)
   INTEGER(IntKi),INTENT(IN):: node_total
   INTEGER(IntKi),INTENT(IN):: dof_node

   REAL(ReKi),INTENT(OUT):: udN(:)

   REAL(ReKi):: ome(3)
   REAL(ReKi):: cc(3)
   REAL(ReKi):: Hinv_temp(3,3)
   INTEGER:: i
   INTEGER:: j
   INTEGER:: temp_id

   udN = 0.0D0
   DO i=1,node_total
       ome = 0.0D0
       cc = 0.0D0
       Hinv_temp = 0.0D0
       DO j=1,3
           temp_id = (i-1)*dof_node+j
           udN(temp_id) = vvN(temp_id)
           ome(j) = vvN(temp_id+3)
           cc(j) = uuN(temp_id+3)
       ENDDO
       CALL CrvMatrixHinv(cc,Hinv_temp)
       cc = MATMUL(Hinv_temp,ome)
       DO j=1,3
           temp_id = (i-1)*dof_node+j
           udN(temp_id+3) = cc(j)
       ENDDO
   ENDDO

   END SUBROUTINE ComputeUDN
