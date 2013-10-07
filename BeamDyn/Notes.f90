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
      
         SUBROUTINE TiSchmComputeCoefficients(beta,gama,deltat,alfaM,alfaF,coef)

   REAL(ReKi),INTENT(IN)::beta, gama, deltat, alfaM, alfaF

   REAL(ReKi),INTENT(INOUT):: coef(:)

   REAL(ReKi)::deltat2, oalfaM, tr0, tr1, tr2

   deltat2 = deltat * deltat
   oalfaM = 1.0D0 - alfaM
   tr0 =  alfaf / oalfaM
   tr1 = alfam / oalfaM
   tr2 = (1.0D0 - alfaF) / oalfaM

   coef(1) = beta * tr0 * deltat2
   coef(2) = (0.5D0 - beta/oalfaM) * deltat2
   coef(3) = gama * tr0 * deltat
   coef(4) = (1.0D0 - gama / oalfaM) * deltat
   coef(5) = tr0
   coef(6) = -tr1
   coef(7) = gama * tr2 * deltat
   coef(8) = beta * tr2 * deltat2
   coef(9) = tr2 

   END SUBROUTINE TiSchmComputeCoefficients


   SUBROUTINE TiSchmComputeParameters()

   REAL(ReKi),INTENT(IN)::rhoinf
   REAL(ReKi),INTENT(OUT)::alfam, alfaf, gama, beta

   REAL(ReKi)::tr0


   tr0 = rhoinf + 1.0D0
   alfam = (2.0D0 * rhoinf - 1.0D0) / tr0
   alfaf = rhoinf / tr0
   gama = 0.5D0 - alfam + alfaf
   beta = 0.25 * (1.0D0 - alfam + alfaf) * (1.0D0 - alfam + alfaf)
   

   END SUBROUTINE TiSchmComputeParameters
