   SUBROUTINE ComputeIniNodalCrvLS(uN0,node_elem,inode,cc)
   !-----------------------------------------------------
   ! This subroutine computes initial CRV parameters
   ! given geometry information
   !-----------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: uN0(:)
   INTEGER(IntKi),INTENT(IN   ):: node_elem
   INTEGER(IntKi),INTENT(IN   ):: inode
   REAL(ReKi),    INTENT(  OUT):: cc(:)       ! Initial Crv Parameter

   REAL(ReKi)                  :: GLL(node_elem)                    ! Derivatives of second-order shape function 
   REAL(ReKi)                  :: w(node_elem)                    ! Derivatives of second-order shape function 
   REAL(ReKi)                  :: hhx(node_elem)                    ! Derivatives of second-order shape function 
   REAL(ReKi)                  :: hpx(node_elem)                    ! Derivatives of second-order shape function 
   REAL(ReKi)                  :: e1(3)                     ! Unit tangent vector
   REAL(ReKi)                  :: e2(3)                     ! Unit normal vector
   REAL(ReKi)                  :: e3(3)                     ! Unit e3 = e1 * e2, cross-product
   REAL(ReKi)                  :: Rr(3,3)                   ! Initial rotation matrix
   REAL(ReKi)                  :: temp
   REAL(ReKi)                  :: temp2
   REAL(ReKi)                  :: Delta
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id

   CALL BD_gen_gll_LSGL(node_elem-1,GLL,w)
   CALL diffmtc(node_elem-1,node_elem,GLL,GLL,inode,hhx,hpx)

   Rr = 0.0D0

   e1 = 0.0D0
   DO i=1,3
       DO j=1,node_elem
           temp_id = (j-1)*6
           e1(i) = e1(i) + hpx(j)*uN0(temp_id+i)
       ENDDO
   ENDDO
   e1 = e1/Norm(e1)
   DO i=1,3
       Rr(i,1) = e1(i)
   ENDDO

   e2 = 0.0D0
   temp_id = (inode-1) * 6 + 4
   temp = uN0(temp_id)
   temp2 = ((e1(2)*COS(temp) + e1(3)*SIN(temp))/e1(1))
   Delta = SQRT(1.0D0 + temp2*temp2)
   e2(1) = -(e1(2)*COS(temp)+e1(3)*SIN(temp))/e1(1)
   e2(2) = COS(temp)
   e2(3) = SIN(temp)
   e2 = e2/Delta
   DO i=1,3
       Rr(i,2) = e2(i)
   ENDDO

   e3 = 0.0D0
   e3 = CrossProduct(e1,e2)
   DO i=1,3
       Rr(i,3) = e3(i)
   ENDDO

   CALL CrvExtractCrv(Rr,cc)

   END SUBROUTINE ComputeIniNodalCrvLS
