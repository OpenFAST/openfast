      program generate_test_data

      implicit double precision (a-h,o-z)

      nb = 10000

      n = 800

      tmin = 0.

      tmax = 10.

      eps = 0.1

      dtb = (tmax - tmin) / nb

      dt = (tmax - tmin) / n

      open (unit = 20, file = 'Bench.dat', status = 'unknown')
      open (unit = 21, file = 'Test.dat', status = 'unknown')


      write(20,*) '# ',nb+1
      do i = 0, nb
      
         t = i*dtb
 
         write(20,*) t, sin(t)

      enddo

      write(21,*) '# ',n+1
      do i = 0, n
      
         t = i*dt
 
         write(21,*) t, sin((1.+eps)*t)

      enddo

      write(*,*) 'test these date in rms_error.f; analytical normalized'
      write(*,*) 'error is 0.6068828437420661'
      write(*,*) 'see verification_mathematica.nb for details'

      end

