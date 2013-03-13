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

      subroutine outputvtk(u, nmax, dofmax, filename)

      implicit double precision (a-h,o-z)
 
      integer dofmax

      double precision u(dofmax)

      integer nelem(2,nmax-1)
 
      integer ihat_array(4)

      double precision xyz(2,nmax)
      double precision u_vtk(nmax)

      character filename*(*)

      character command*80

      numel = nmax - 1

      ihat = 0

      do i = 1, nmax 

        ilocal = (i-1)*3 + 1

        ihat =  ihat + 1 

        xyz(1,ihat) = u(ilocal)
        xyz(2,ihat) = u(ilocal+1)

        u_vtk(ihat) = 0.

      enddo

      ihat_array(1) = 0

      ielem = 0
      ihat = 0
      do i = 1, numel
        ihat = ihat + 1
        ielem = ielem + 1 
        nelem(1,ielem) = ihat
        nelem(2,ielem) = ihat + 1
      enddo

      numnodes = nmax 

      open (unit=20, file=filename//'.vtk', status='unknown')

      write(20,*)'# vtk DataFile Version 2.0'

      write(20,*) 'Really cool data'

      write(20,*) 'ASCII'

      write(20,*) 'DATASET UNSTRUCTURED_GRID'

      write(20,*) 'POINTS', numnodes, 'float'

      do i = 1,numnodes
        write(20,991) xyz(1,i)
        write(20,991) xyz(2,i)
        write(20,991) 0.
        !write(20,*) 0.1 * temp(i) / dmax
        write(20,*) ' '
      enddo

      write(20,*) 'CELLS', numel, 3*numel

      do i = 1,numel
        write(20,*) '2', (nelem(l,i)-1,l = 1,2)
      enddo

      write(20,*) 'CELL_TYPES', numel
      write(20,*) (3, i = 1,numel)

      write(20,*) 'POINT_DATA', numnodes
      write(20,*) 'SCALARS ','zero',' float 1'
      write(20,*) 'LOOKUP_TABLE default'

      do i = 1,numnodes
         write(20,991) 0.
      enddo

      call flush(20)

      close(unit=20)

      ifinal = 1

      if (ifinal .eq. 1) then
        !command = '/usr/bin/asa '//filename//'.vtk > tmp.vtk'
        command = 'asa/asa '//filename//'.vtk > tmp.vtk'
        call system(command)
        command = '/bin/mv -f tmp.vtk '//filename//'.vtk'
        call system(command)
      endif

991   format(E16.8)
993   format(3E16.8)

      return
      end

