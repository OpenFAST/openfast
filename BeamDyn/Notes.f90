          Nrrr_temp = 0.0d0
          IF(i==1) THEN
              Nuuu_temp1(1) = Nuuu(temp_id+4)
              Nuuu_temp1(2) = Nuuu(temp_id+5)
              Nuuu_temp1(3) = Nuuu(temp_id+6)
          ENDIF
          Nuuu_temp(1) = Nuuu(temp_id+4)
          Nuuu_temp(2) = Nuuu(temp_id+5)
          Nuuu_temp(3) = Nuuu(temp_id+6)
          CALL CrvCompose(Nrrr_temp,Nuuu_temp1,Nuuu_temp,1)
          temp_id = ((nelem-1)*norder+i-1)*3
          IF(nelem==1) THEN
              Nrrr(temp_id+1) = Nuuu_temp(1)
              Nrrr(temp_id+2) = Nuuu_temp(2)
              Nrrr(temp_id+3) = Nuuu_temp(3)
          ELSEIF(i/=1) THEN
              Nrrr(temp_id+1) = Nuuu_temp(1)
              Nrrr(temp_id+2) = Nuuu_temp(2)
              Nrrr(temp_id+3) = Nuuu_temp(3)
          ENDIF
          
          
          
          
          
          Nrr0_temp = 0.0d0
          IF(i==1) THEN
              Nuu0_temp1(1) = Nuu0(temp_id+4)
              Nuu0_temp1(2) = Nuu0(temp_id+5)
              Nuu0_temp1(3) = Nuu0(temp_id+6)
          ENDIF
          Nuu0_temp(1) = Nuu0(temp_id+4)
          Nuu0_temp(2) = Nuu0(temp_id+5)
          Nuu0_temp(3) = Nuu0(temp_id+6)
          CALL CrvCompose(Nrr0_temp,Nuu0_temp1,Nuu0_temp,1)
          temp_id = ((nelem-1)*norder+i-1)*3
          IF(nelem==1) THEN
              Nrr0(temp_id+1) = Nrr0_temp(1)
              Nrr0(temp_id+2) = Nrr0_temp(2)
              Nrr0(temp_id+3) = Nrr0_temp(3)
          ELSEIF(i/=1) THEN
              Nrr0(temp_id+1) = Nrr0_temp(1)
              Nrr0(temp_id+2) = Nrr0_temp(2)
              Nrr0(temp_id+3) = Nrr0_temp(3)
          ENDIF 
          
          
      DO i=1,node_elem
          CALL InertialForce()
          CALL ElasticForce()
          DO j=1,dof_node
              temp_count1 = (i-1)*dof_node+j
              DO k=1,dof_node
                  temp_count2 = (i-1)*dof_node+k
                  elk(temp_count1,temp_count2) = w(i)* Ki(j,k)*Jac
                  elk(temp_count1,temp_count2) = elk(temp_count1,temp_count2)+w(i)* Qe(j,k)*Jac
              ENDDO
          ENDDO
          DO j=1,node_elem
              DO k=1,dof_node
                  temp_count1 = (i-1)*dof_node+k
                  DO m=1,dof_node
                      temp_count2 = (j-1)*dof_node+m
                      elk(temp_count1,temp_count2) = elk(temp_count1,temp_count2)+w(i)*Pe(k,m)*hhp(j,i)
                  ENDDO
              ENDDO              
          ENDDO
          DO j=1,node_elem
              DO k=1,dof_node
                  temp_count1 = (j-1)*dof_node+k
                  DO m=1,dof_node
                      temp_count2 = (i-1)*dof_node+m
                      elk(temp_count1,temp_count2) = elk(temp_count1,temp_count2)+w(i)*Oe(k,m)*hhp(j,i)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO
      DO i=1,node_elem
          DO j=1,node_elem
              DO k=1,dof_node
                  temp_count1 = (i-1)*dof_node+k
                  DO m=1,dof_node
                      temp_count2 = (j-1)*dof_node+m
                      DO n=1,node_elem
                          CALL ElasticForce()
                          elk(temp_count1,temp_count2) = elk(temp_count1,temp_count2)+w(n)*hhp(i,n)*hhp(j,n)*Se(k,m)
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO