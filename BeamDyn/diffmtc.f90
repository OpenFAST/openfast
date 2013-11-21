      SUBROUTINE diffmtc(np,ns,spts,npts,igp,hhx,hpx)
!
! calculate Lagrangian interpolant tensor at ns points where basis
! functions are assumed to be associated with (np+1) GLL points on
! [-1,1]
!
!     INPUT:
!       np             : polynomial order of basis funcitons
!       ns             : number of points at which to eval shape/deriv
!       spts(ns)       : location of ns points at which to eval
!       npts(np+1)     : location of the (np+1) GLL points
!
!     OUTPUT
!       dPhis(np+1,ns) : derivative evaluated at ns points
!       ps(np+1,ns)    : (np+1) shape functions evaluated at ns points 

      INTEGER(IntKi),INTENT(IN):: np,ns,igp
      REAL(ReKi),INTENT(IN):: spts(:)
      REAL(ReKi),INTENT(IN):: npts(:)
      
      REAL(ReKi),INTENT(OUT):: hhx(:),hpx(:) 
      
      REAL(ReKi):: dPhis(np+1,ns),Ps(np+1,ns)
      
      REAL(ReKi):: dnum,den
      REAL(ReKi),PARAMETER:: eps = 1.0D-08
      INTEGER(IntKi):: l,j,i,k
      
      do l = 1,np+1
        do j = 1,ns
          dPhis(l,j) = 0.
          den = 1.
          if ((abs(spts(j)-1.).LE.eps).AND.(l.EQ.np+1)) then
            dPhis(l,j) = float((np+1)*np)/4.
          elseif ((abs(spts(j)+1.).LE.eps).AND.(l.EQ.1)) then
            dPhis(l,j) = -float((np+1)*np)/4.
          elseif (abs(spts(j)-npts(l)).LE.eps) then
            dPhis(l,j) = 0.
          else
            do i = 1,np+1
              if (i.NE.l) then
                den = den*(npts(l)-npts(i))
              endif
              dnum = 1.
              do k = 1,np+1
                if ((k.NE.l).AND.(k.NE.i).AND.(i.NE.l)) then
                  dnum = dnum*(spts(j)-npts(k))
                elseif (i.EQ.l) then
                  dnum = 0.
                endif
              enddo
              dPhis(l,j) = dPhis(l,j) + dnum
            enddo
            dPhis(l,j) = dPhis(l,j)/den
          endif
        enddo
      enddo

      do l = 1,np+1
        do j = 1,ns
          Ps(l,j) = 0.
          dnum = 1.
          den = 1.
          if(abs(spts(j)-npts(l)).LE.eps) then
            Ps(l,j) = 1.
          else
            do k = 1,np+1
              if (k.NE.l) then
                den = den*(npts(l) - npts(k))
                dnum = dnum*(spts(j) - npts(k))
              endif
            enddo
            Ps(l,j) = dnum/den
          endif
        enddo
      enddo
      
      DO i=1,np+1
         hhx(i) = Ps(i,igp)
         hpx(i) = dPhis(i,igp)
      ENDDO

      END SUBROUTINE diffmtc
