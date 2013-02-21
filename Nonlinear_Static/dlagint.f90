      function dlagint(i, j, order,xj)


      integer i, j, order
      double precision xj(order+1)
      double precision dlagint

      double precision dlegen
      external dlegen


      dlagint = 0.d0
 

      if (i.eq.j) then
       if (i.eq.1) then
        dlagint = -dble((order+1)*order)/4.d0
       else if (i.eq.(order+1)) then
        dlagint = +dble((order+1)*order)/4.d0
       else
        dlagint = 0.d0 
       end if
      else
       dlagint = dlegen(order,xj(j))/(dlegen(order,xj(i))*(xj(j)-xj(i)))
      end if

      return
      end function
