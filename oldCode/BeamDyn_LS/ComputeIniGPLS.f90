   SUBROUTINE ComputeIniGPLS(uN0,node_elem,inode,cc)
   !-----------------------------------------------------
   ! This subroutine computes initial CRV parameters
   ! given geometry information
   !-----------------------------------------------------

   REAL(ReKi),    INTENT(IN   ):: uN0(:)
   INTEGER(IntKi),INTENT(IN   ):: node_elem
   INTEGER(IntKi),INTENT(IN   ):: inode
   REAL(ReKi),    INTENT(  OUT):: cc(:)       ! Initial Gauss point location

   REAL(ReKi)                  :: GLL(node_elem)                    ! Derivatives of second-order shape function 
   REAL(ReKi)                  :: w(node_elem)                    ! Derivatives of second-order shape function 
   REAL(ReKi)                  :: hhx(node_elem)                    ! Derivatives of second-order shape function 
   REAL(ReKi)                  :: hpx(node_elem)                    ! Derivatives of second-order shape function 
   REAL(ReKi)                  :: gp(node_elem-1)                     ! Unit tangent vector
   REAL(ReKi)                  :: gw(node_elem-1)                     ! Unit normal vector
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: temp_id

   CALL BldGaussPointWeight(node_elem-1,gp,gw)
   CALL BD_gen_gll_LSGL(node_elem-1,GLL,w)
   CALL diffmtc(node_elem-1,node_elem-1,gp,GLL,inode,hhx,hpx)

   cc = 0.0D0
   DO i=1,3
       DO j=1,node_elem
           temp_id = (j-1)*6
           cc(i) = cc(i) + hhx(j)*uN0(temp_id+i)
       ENDDO
   ENDDO

   END SUBROUTINE ComputeIniGPLS
