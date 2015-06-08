!***********************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory (NREL)
!
!    This is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the
!    License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License 
!    along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!***********************************************************************
!    This code was created at NREL by Michael A. Sprague and Ignas 
!    Satkauskas  and was meant for open-source distribution.
!
!    Software was created under funding from a Shared Research Grant 
!    from the Center for Research and Education in Wind (CREW), during
!    the period 01 October 2011 - 31 January 2013.
!
!    http://crew.colorado.edu/ 
!
!    Questions?  Please contact Michael Sprague:
!    email:  michael.a.sprague@nrel.gov
!
!***********************************************************************

      function filename(i)

c     takes integer i, and returns a character expression of length 4
c     
c     For example,
c     filename(5) should return a character '0005'
c
c     or 
c
c     filename(15) should return a character '0015'

      integer i

      character(1) ione
      character(2) itwo
      character(3) ithree
      character(4) ifour

      character(4) filename

      if (i .lt. 10) then

        write(ione,'(i1)' ) i
        filename = '000'//ione
 
      elseif (i.lt. 100) then
     
        write(itwo,'(i2)' ) i
        filename = '00'//itwo

      elseif (i.lt. 1000) then
        write(ithree,'(i3)' ) i
        filename = '0'//ithree

      elseif (i.lt. 10000) then
        write(ifour,'(i4)' ) i
        filename = ifour

      else

        stop 'integer too big for function filename'

      endif

!1001  format(I1)
!1002  format(I2)
!1003  format(I3)
!1004  format(I4)

      return
      end 

