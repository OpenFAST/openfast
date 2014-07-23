   SUBROUTINE ComputeIniNodalTwist(member_length,station_eta,IniTwist0,&
                                  &node_elem,elem_total,station_total,&
                                  &IniTwist_Nodal,ErrStat,ErrMsg)

   REAL(ReKi),    INTENT(IN   ):: member_length(:,:)
   REAL(ReKi),    INTENT(IN   ):: station_eta(:)
   REAL(ReKi),    INTENT(IN   ):: IniTwist0(:)
   INTEGER(IntKi),INTENT(IN   ):: node_elem
   INTEGER(IntKi),INTENT(IN   ):: elem_total
   INTEGER(IntKi),INTENT(IN   ):: station_total
   REAL(ReKi),    INTENT(  OUT):: IniTwist_Nodal(:)
   INTEGER(IntKi),INTENT(  OUT):: ErrStat
   CHARACTER(*),  INTENT(  OUT):: ErrMsg

   REAL(ReKi),ALLOCATABLE      :: temp_ratio(:,:)
   REAL(ReKi),ALLOCATABLE      :: temp_GLL(:)
   REAL(ReKi),ALLOCATABLE      :: temp_w(:)
   INTEGER(IntKi)              :: i 
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: k 
   INTEGER(IntKi)              :: temp_id

   CALL AllocAry(temp_ratio,node_elem,elem_total,'temp_ratio',ErrStat,ErrMsg)
   temp_ratio(:,:) = 0.0D0
   CALL AllocAry(temp_GLL,node_elem,'temp_GLL',ErrStat,ErrMsg)
   temp_GLL(:) = 0.0D0
   CALL AllocAry(temp_w,node_elem,'temp_weight_GLL',ErrStat,ErrMsg)
   temp_w(:) = 0.0D0

   CALL BeamDyn_gen_gll_LSGL(node_elem-1,temp_GLL,temp_w)

   DO i=1,node_elem
       temp_GLL(i) = (temp_GLL(i) + 1.0D0)/2.0D0
   ENDDO

   DO i=1,elem_total
       IF(i .EQ. 1) THEN
           DO j=1,node_elem
               temp_ratio(j,i) = temp_GLL(j)*member_length(i,2)
           ENDDO
       ELSE
           DO j=1,i-1
               temp_ratio(:,i) = temp_ratio(:,i) + member_length(j,2)
           ENDDO
           DO j=1,node_elem
               temp_ratio(j,i) = temp_ratio(j,i) + temp_GLL(j)*member_length(i,2)
           ENDDO
       ENDIF
   ENDDO

   i=1
   DO j=1,node_elem
       temp_id = j
       DO k=1,station_total
           IF(temp_ratio(j,i) <= station_eta(k)) THEN
               IF(temp_ratio(j,i) == station_eta(k)) THEN
                   IniTwist_Nodal(temp_id) = IniTwist0(k)
               ELSE
                   IniTwist_Nodal(temp_id) = 0.5D0*(IniTwist0(k-1) + IniTwist0(k))
               ENDIF
               EXIT
           ENDIF
       ENDDO
   ENDDO
   DO i=2,elem_total
       DO j=2,node_elem
           temp_id = (i-1)*(node_elem-1)+j
           DO k=1,station_total
               IF(temp_ratio(j,i) <= station_eta(k)) THEN
                   IF(temp_ratio(j,i) == station_eta(k)) THEN
                       IniTwist_Nodal(temp_id) = IniTwist0(k)
                   ELSE
                       IniTwist_Nodal(temp_id) = 0.5D0*(IniTwist0(k-1) + IniTwist0(k))
                   ENDIF
                   EXIT
               ENDIF
           ENDDO
       ENDDO
   ENDDO

   DEALLOCATE(temp_ratio)
   DEALLOCATE(temp_GLL)
   DEALLOCATE(temp_w)

   END SUBROUTINE ComputeIniNodalTwist
