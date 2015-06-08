      program rms_error

      implicit double precision (a-h,o-z)

      parameter (imax = 5000000)
      
      double precision t(imax), u(imax)
      double precision t_temp(3), u_temp(3)  ! four-point interp
      double precision tb(imax), ub(imax)
      double precision ub_interp(imax)

      character tmp

      external polint
      external locate

      open (unit = 20, file = 'Bench.dat', status = 'old')
      open (unit = 21, file = 'Test.dat', status = 'old')

! read in benchmark data
! first line of Bench.dat should be # followed by number of lines of data, 
! e.g.,
! # 2484
      read(20,*) tmp, numpts_b
      do i = 1, numpts_b
        read(20,*) tb(i), ub(i)
      enddo
      write(*,*) 'numpts_b is', numpts_b

! first line of Test.dat should be # followed by number of lines of data, 
! e.g.,
! # 2484

      read(21,*) tmp, numpts
! read in benchmark data
      do i = 1, numpts
        read(21,*) t(i), u(i)
      enddo
      write(*,*) 'numpts is', numpts



C interpolate ub to same grid as u, store in ub_interp
      do i = 1, numpts
        time = t(i)
        call locate(tb, numpts_b, time, locat)
        if (locat.eq.0) locat= 1
        do j = 1,3
          t_temp(j) = tb(locat+j-1)
          u_temp(j) = ub(locat+j-1)
        enddo
        call polint(t_temp, u_temp,3,time,ub_out,dy)
        ub_interp(i) = ub_out
        !write(*,*) time, ub_out
      enddo

C calculate normalized rms error and the normalized inf error
C use trap rule integration

      rms = 0.d0
      dnorm = 0.d0

      error_max = 0.d0

      do i = 1, numpts-1

         dmid_u  = (u(i) + u(i+1)) / 2.d0

         dmid_ub = (ub_interp(i) + ub_interp(i+1)) / 2.d0

         error = abs(ub_interp(i) - u(i))
         if (error .gt. error_max) error_max = error
 
         rms = rms + (dmid_u - dmid_ub)**2

         dnorm = dnorm + dmid_ub**2

      enddo

      total_time = t(numpts)

      dnorm = sqrt(total_time * dnorm / dble(numpts-1))

      rms = sqrt(total_time * rms / dble(numpts-1))

      write(*,*) 'errornorm ', dnorm

      write(*,*) 'normalized rms error = ', rms/dnorm

      write(*,*) 'error_max = ', error_max

      end

c--------------------------------------------------------------------
      subroutine locate(xx,n,x,j)
      integer j,n
      double precision x,xx(n)
      integer jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      end
c  (c) copr. 1986-92 numerical recipes software +k$<,(5cl.
c--------------------------------------------------------------------
      subroutine polint(xa,ya,n,x,y,dy)
      integer n,nmax
      double precision dy,x,y,xa(n),ya(n)
      parameter (nmax=10)
      integer i,m,ns
      double precision den,dif,dift,ho,hp,w,c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) stop 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
c  (c) copr. 1986-92 numerical recipes software +k$<,(5cl.
c--------------------------------------------------------------------
