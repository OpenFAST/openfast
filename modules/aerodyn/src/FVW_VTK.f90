module FVW_VTK
    !use PrecisionMod, only: ReKi
    use NWTC_Library, only: ReKi, GetNewUnit
    implicit none
!     character(8), parameter :: RFMT='F14.5'
    !character(8), parameter :: RFMT='E24.15E3'
    character(8), parameter :: RFMT='E17.8E3'
    character(8), parameter :: IFMT='I7'

   TYPE, PUBLIC :: FVW_VTK_Misc
      integer :: vtk_unit
      logical :: bFileOpen=.false.

      integer :: nData=0;
      integer :: nPoints=0;

      logical :: bBinary = .false.
      character(len=255) :: buffer

      ! Reference Frame
      real(ReKi),dimension(3,3) :: T_g2b
      real(ReKi),dimension(3)   :: PO_g
   END TYPE FVW_VTK_Misc

    character(1), parameter :: NL = char(10) ! New Line character

    interface vtk_dataset_structured_grid; module procedure &
            vtk_dataset_structured_grid_flat, &
            vtk_dataset_structured_grid_grid
    end interface

    interface vtk_point_data_vector; module procedure &
            vtk_point_data_vector_flat, &
            vtk_point_data_vector_grid2D,&
            vtk_point_data_vector_grid
    end interface
    interface vtk_point_data_scalar; module procedure &
            vtk_point_data_scalar_flat, &
            vtk_point_data_scalar_grid2D, &
            vtk_point_data_scalar_grid
    end interface
    interface vtk_cell_data_scalar; module procedure &
            vtk_cell_data_scalar_1d,&
            vtk_cell_data_scalar_2d 
    end interface

    public

contains

   subroutine vtk_misc_init(mvtk)
      type(FVW_VTK_Misc),intent(inout) :: mvtk
      mvtk%vtk_unit = -1           !< VTK output unit [-]
      mvtk%bFileOpen = .false.     !< binary file is open [-]
      mvtk%bBinary = .false.       !< write binary files [-]
      mvtk%nData = 0               !< number of data lines [-]
      mvtk%nPoints = 0             !< number of points [-]
   end subroutine

    !>
    subroutine set_vtk_binary_format(bBin,mvtk)
        logical, intent(in)::bBin
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        mvtk%bBinary=bBin
    end subroutine

    
    !> Save a coordinate transform
    ! ALL VTK Will be exported in this coordinate system!
    subroutine set_vtk_coordinate_transform(T_g2b_in,PO_g_in,mvtk)
        real(ReKi),dimension(3,3), intent(in) :: T_g2b_in
        real(ReKi),dimension(3)  , intent(in) :: PO_g_in
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        mvtk%T_g2b=T_g2b_in
        mvtk%PO_g=PO_g_in
    end subroutine

    logical function vtk_new_ascii_file(filename,label,mvtk)
        !use MainIO,     only: get_free_unit ,check_io
        !use MainIOData, only: bSTOP_ALLOWED
        !use FileSystem, only: file_exists
        !use Logging,    only: log_warning,log_error,log_info
        !
        character(len=*),intent(in)      :: filename
        character(len=*),intent(in)      :: label
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        !
        integer :: iostatvar
        logical :: b

        if (.not. mvtk%bFileOpen) then
            CALL GetNewUnit( mvtk%vtk_unit )   
            if (mvtk%bBinary) then
                ! Fortran 2003 stream, otherwise intel fortran !
                !form='UNFORMATTED',access='SEQUENTIAL',action='WRITE',convert='BIG_ENDIAN',recordtype='STREAM',buffered='YES',
               !print*,'Not available for this compiler' !COMPAQ-COMPILER
               !STOP !COMPAQ-COMPILER
!bjj: CONVERT is non-standard, so maybe this should be part of Sys*.f90? Like OpenUnfInpBEFile()?
                open(unit = mvtk%vtk_unit,file= trim(adjustl(filename)),form='UNFORMATTED',access = 'stream',& !OTHER-COMPILER
                    action = 'WRITE',convert= 'BIG_ENDIAN',iostat=iostatvar,status='replace') !OTHER-COMPILER
            else
                open(mvtk%vtk_unit,file=trim(adjustl(filename)),iostat=iostatvar,action="write",status='replace')
            endif
            if (iostatvar == 0) then
                if (mvtk%bBinary) then
                    write(mvtk%vtk_unit)'# vtk DataFile Version 3.0'//NL
                    write(mvtk%vtk_unit)trim(label)//NL
                    write(mvtk%vtk_unit)'BINARY'//NL
                else
                    write(mvtk%vtk_unit,'(a)') '# vtk DataFile Version 2.0'
                    write(mvtk%vtk_unit,'(a)') trim(label)
                    write(mvtk%vtk_unit,'(a)') 'ASCII'
                    write(mvtk%vtk_unit,'(a)') ' '
                endif

                mvtk%bFileOpen=.true.
                mvtk%nData=-1;
            endif
        else
            b=.false.
            !call log_error('VTK: Cannot open two vtk files at the same time, call vtk_close first')
        endif
        if (iostatvar ==0) then
           vtk_new_ascii_file=.true.
        else
           vtk_new_ascii_file=.false.
        endif
    end function

    subroutine vtk_close_file(mvtk)
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            close(mvtk%vtk_unit)
            mvtk%bFileOpen=.false.
        endif
    endsubroutine


    ! ------------------------------------------------------------------------- 
    ! --- POLYDATA STUFF
    ! ------------------------------------------------------------------------- 
    subroutine vtk_dataset_polydata(Points,mvtk,bladeFrame)
        real(ReKi), dimension(:,:),intent(in) :: Points  !< 3 x n
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        logical, intent(in) :: bladeFrame
        integer :: i
        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=size(Points,2)
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'DATASET POLYDATA'//NL
                write(mvtk%buffer,'(A,I0,A)') 'POINTS ', mvtk%nPoints ,' double'
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
                if (bladeFrame)  then
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit)matmul(mvtk%T_g2b,Points(1:3,i)-mvtk%PO_g)
                    enddo
                else
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit)Points(1:3,i)
                    enddo
                endif
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET POLYDATA'
                write(mvtk%vtk_unit,'(A,I0,A)') 'POINTS ', mvtk%nPoints ,' double'
                if (bladeFrame)  then
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit,'(3'//RFMT//')') matmul(mvtk%T_g2b,Points(1:3,i)-mvtk%PO_g)
                    enddo
                else
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit,'(3'//RFMT//')') Points(1:3,i)
                    enddo
                endif
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine


    subroutine vtk_lines(L,mvtk)
        integer, dimension(:,:),intent(in) :: L    !< 2 x n
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        integer :: i

        if ( mvtk%bFileOpen ) then
            mvtk%nData=size(L,2)
            if (mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0,A,I0)')'LINES ',mvtk%nData,' ',3*mvtk%nData
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit)2,L(1:2,i)
                enddo
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,I0,A,I0)')'LINES ',mvtk%nData,' ',3*mvtk%nData
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit,'(3'//IFMT//')') 2, L(1:2,i)
                enddo
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine

    subroutine vtk_quad(Q,mvtk)
        integer, dimension(:,:),intent(in) :: Q    !< 4 x n
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        integer :: i
        if ( mvtk%bFileOpen ) then
            mvtk%nData=size(Q,2)
            if (mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0,A,I0)')'POLYGONS ',mvtk%nData,' ',5*mvtk%nData
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit)4,Q(1:4,i)
                enddo
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,I0,A,I0)') 'POLYGONS ', mvtk%nData,' ',5*mvtk%nData
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit,'(5'//IFMT//')') 4, Q(1:4,i)
                enddo
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine

    ! ------------------------------------------------------------------------- 
    ! --- RECTILINEAR
    ! ------------------------------------------------------------------------- 
    subroutine vtk_dataset_rectilinear(v1,v2,v3,mvtk)
        real(ReKi), dimension(:),intent(in) :: v1,v2,v3  !<  n
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=size(v1)*size(v2)*size(v3)
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET RECTILINEAR_GRID'//NL
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', size(v1),' ',size(v2),' ',size(v3)
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,I0,A)') 'X_COORDINATES ', size(v1), ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)v1
                write(mvtk%vtk_unit)NL
                write(mvtk%buffer,'(A,I0,A)') 'Y_COORDINATES ', size(v2), ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)v2
                write(mvtk%vtk_unit)NL
                write(mvtk%buffer,'(A,I0,A)') 'Z_COORDINATES ', size(v3), ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)v3
                !write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET RECTILINEAR_GRID'
                write(mvtk%vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', size(v1),' ',size(v2),' ',size(v3)
                write(mvtk%vtk_unit,'(A,I0,A)') 'X_COORDINATES ', size(v1), ' double'
                write(mvtk%vtk_unit,'('//RFMT//')') v1
                write(mvtk%vtk_unit,'(A,I0,A)') 'Y_COORDINATES ', size(v2), ' double'
                write(mvtk%vtk_unit,'('//RFMT//')') v2
                write(mvtk%vtk_unit,'(A,I0,A)') 'Z_COORDINATES ', size(v3), ' double'
                write(mvtk%vtk_unit,'('//RFMT//')') v3
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine

    subroutine vtk_dataset_structured_points(x0,dx,n,mvtk)
        real(ReKi),     dimension(3), intent(in) :: x0 !< origin
        real(ReKi),     dimension(3), intent(in) :: dx !< spacing
        integer,        dimension(3), intent(in) :: n  !< length
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=n(1)*n(2)*n(3)
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET STRUCTURED_POINTS'//NL
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ',n(1),' ',n(2),' ',n(3)
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,3F16.8)') 'ORIGIN ', x0
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,3F16.8)') 'SPACING ', dx
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET STRUCTURED_POINTS'
                write(mvtk%vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n(1),' ',n(2),' ',n(3)
                write(mvtk%vtk_unit,'(A,3F16.8,A)') 'ORIGIN  ',x0
                write(mvtk%vtk_unit,'(A,3F16.8,A)') 'SPACING ',dx
            endif
        endif
    end subroutine


    ! ------------------------------------------------------------------------- 
    ! --- STRUCTURED GRID (Points dumped without for loop since memory is in proper order)
    ! ------------------------------------------------------------------------- 
    !> Subroutine using flat data as input (not in natural order)
    subroutine vtk_dataset_structured_grid_flat(D,n1,n2,n3,mvtk)
        integer , intent(in) :: n1,n2,n3
        real(ReKi), dimension(:,:),intent(in)::D
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=n1*n2*n3
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET STRUCTURED_GRID'//NL
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET STRUCTURED_GRID'
                write(mvtk%vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine

    !> Using Grid data 4d as input
    subroutine vtk_dataset_structured_grid_grid(D,n1,n2,n3,mvtk)
        integer , intent(in) :: n1,n2,n3
        real(ReKi), dimension(:,:,:,:),intent(in)::D
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=n1*n2*n3
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET STRUCTURED_GRID'//NL
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%buffer,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A)') 'DATASET STRUCTURED_GRID'
                write(mvtk%vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
                write(mvtk%vtk_unit,*) ' '
            endif
        endif
    end subroutine
    
    
    
    ! ------------------------------------------------------------------------- 
    ! --- POINT DATA
    ! ------------------------------------------------------------------------- 
    subroutine vtk_point_data_init(mvtk)
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if(mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0)')'POINT_DATA ',mvtk%nPoints
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
            else
                write(mvtk%vtk_unit,'(A,I0)') 'POINT_DATA ', mvtk%nPoints
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_flat(D,sname,mvtk)
        real(ReKi), dimension(:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_grid(D,sname,mvtk)
        real(ReKi), dimension(:,:,:,:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_grid2D(D,sname,mvtk)
        real(ReKi), dimension(:,:,:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    !>
    subroutine vtk_point_data_vector_flat(D,sname,mvtk)
        real(ReKi), dimension(:,:),intent(in) :: D  !< 3 x n
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine
    !>
    subroutine vtk_point_data_vector_grid(D,sname,mvtk)
        real(ReKi), dimension(:,:,:,:),intent(in) :: D  !< 3 x n
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine
    !>
    subroutine vtk_point_data_vector_grid2D(D,sname,mvtk)
        real(ReKi), dimension(:,:,:),intent(in) :: D  !< 
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine


    ! ------------------------------------------------------------------------- 
    ! --- CELL DATA
    ! ------------------------------------------------------------------------- 
    subroutine vtk_cell_data_init(mvtk)
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0)')'CELL_DATA ',mvtk%nData
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NL
            else
                write(mvtk%vtk_unit,'(A,I0)') 'CELL_DATA ', mvtk%nData
            endif
        endif
    end subroutine

    subroutine vtk_cell_data_scalar_1d(D,sname,mvtk)
        real(ReKi), dimension(:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double 1'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,fmt='(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    subroutine vtk_cell_data_scalar_2d(D,sname,mvtk)
        real(ReKi), dimension(:,:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double 1'//NL
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,fmt='(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine
    
    
    subroutine vtk_cell_data_vector(D,sname,mvtk)
        real(ReKi), dimension(:,:),intent(in) :: D  !< 3 x n
        character(len=*),intent(in) ::sname
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NL
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(mvtk%vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine

    ! --------------------------------------------------------------------------------}
    ! --- VTK Tools 
    ! --------------------------------------------------------------------------------{
    !> Exports a Plane From a mesh
    subroutine export_plane_grid3d(fname,v1,v2,v3,Values,mvtk)
        character(len=*),intent(in)             :: fname
        real(ReKi),dimension(:), intent(in)       :: v1,v2,v3
        real(ReKi),dimension(:,:,:,:), intent(in) :: Values
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        !  Variables
        integer :: nD

        ! Writting
        if ( vtk_new_ascii_file(trim(fname),'grid',mvtk)) then
            nD=size(Values,1)
            call vtk_dataset_rectilinear(v1,v2,v3,mvtk)
            ! Output as a structured grid, No need to reorder
            call vtk_point_data_init(mvtk)
            ! Could be a function of nDim, be careful
            if(nD==3) then
                call vtk_point_data_vector(Values(1:3,:,:,:),'Velocity',mvtk) ! Label...
            endif

            call vtk_close_file(mvtk)
        endif ! file opening
    end subroutine
    
    !> Exports a Plane From a mesh
    subroutine export_plane_grid2d(fname,v1,v2,v3,Values,mvtk)
        character(len=*),intent(in)             :: fname
        real(ReKi),dimension(:), intent(in)       :: v1,v2,v3
        real(ReKi),dimension(:,:,:), intent(in) :: Values
        type(FVW_VTK_Misc),intent(inout) :: mvtk
        !  Variables
        integer :: nD

        ! Writting
        if ( vtk_new_ascii_file(trim(fname),'plane',mvtk) ) then
            nD=size(Values,1)
            call vtk_dataset_rectilinear(v1,v2,v3,mvtk)
            ! Output as a structured grid, No need to reorder
            call vtk_point_data_init(mvtk)
            ! Could be a function of nDim, be careful
            if(nD==3) then
                call vtk_point_data_vector(Values(1:3,:,:),'Velocity',mvtk) ! Label...
            endif

            call vtk_close_file(mvtk)
        endif ! file opening
    end subroutine
end module FVW_VTK
