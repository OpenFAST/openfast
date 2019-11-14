module VTK
    !use PrecisionMod, only: ReKi
    use NWTC_Library, only: ReKi
    implicit none
!     character(8), parameter :: RFMT='F14.5'
    !character(8), parameter :: RFMT='E24.15E3'
    character(8), parameter :: RFMT='E17.8E3'
    character(8), parameter :: IFMT='I7'

    integer, save :: vtk_unit
    logical, save :: bFileOpen=.false.

    integer, save :: nData=0;
    integer, save :: nPoints=0;

    logical, save :: bOverWritWarned = .true.


    logical, save :: bBinary = .false.
    character(len=255), save :: buffer
    character(1), parameter :: NL = char(10) ! New Line character


    ! Reference Frame
    logical,save                 :: bChangeFrame =.false.
    real(ReKi),dimension(3,3),save :: T_g2b
    real(ReKi),dimension(3)  ,save :: PO_g


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

    public
    private:: vtk_unit, bFileOpen, nData, nPoints,bBinary, buffer
    private :: bChangeFrame, T_g2b, PO_g

contains


    !>
    subroutine set_vtk_binary_format(bBin)
        logical, intent(in)::bBin
        bBinary=bBin
    end subroutine

    
    !> Save a coordinate transform
    ! ALL VTK Will be exported in this coordinate system!
    subroutine set_vtk_coordinate_transform(T_g2b_in,PO_g_in)
        real(ReKi),dimension(3,3), intent(in) :: T_g2b_in
        real(ReKi),dimension(3)  , intent(in) :: PO_g_in
        !
        bChangeFrame=.true.
        T_g2b=T_g2b_in
        PO_g=PO_g_in
    end subroutine



    logical function vtk_new_ascii_file(filename,label,bDEBUG)
        !use MainIO,     only: get_free_unit ,check_io
        !use MainIOData, only: bSTOP_ALLOWED
        !use FileSystem, only: file_exists
        !use Logging,    only: log_warning,log_error,log_info
        !
        character(len=*),intent(in) :: filename
        character(len=*),intent(in) :: label
        logical, intent(in),optional ::bDEBUG
        !
        integer :: iostatvar
        character(len=255) :: iomessage
        logical :: b

        if (.not.bFileOpen) then
            !if (present(bDEBUG)) then
            !    if (bDEBUG) call log_info('Opening VTK file:'//filename)
            !endif

            !if (file_exists(filename) .and. .not.bOverWritWarned) then 
                !call log_warning('Overwritting vtk file '//filename)
                !call log_warning('Further overwritting warnings will not be displayed')
            !    bOverWritWarned=.true.
            !endif

            !vtk_unit=get_free_unit();
            vtk_unit=123455
            if (bBinary) then
                ! Fortran 2003 stream, otherwise intel fortran !
                !form='UNFORMATTED',access='SEQUENTIAL',action='WRITE',convert='BIG_ENDIAN',recordtype='STREAM',buffered='YES',
               !print*,'Not available for this compiler' !COMPAQ-COMPILER
               !STOP !COMPAQ-COMPILER
                open(unit = vtk_unit,file= trim(adjustl(filename)),form='UNFORMATTED',access = 'stream',& !OTHER-COMPILER
                    action = 'WRITE',convert= 'BIG_ENDIAN',iostat=iostatvar,status='replace') !OTHER-COMPILER
            else
                open(vtk_unit,file=trim(adjustl(filename)),iostat=iostatvar,action="write",status='replace')
            endif
            !b=check_io(iostatvar,iomessage);
            !if(.not.b .and. bSTOP_ALLOWED) then
            !    STOP ! io errors are fatal
            !endif

            if (iostatvar == 0) then
                if (bBinary) then
                    write(vtk_unit)'# vtk DataFile Version 3.0'//NL
                    write(vtk_unit)trim(label)//NL
                    write(vtk_unit)'BINARY'//NL
                else
                    write(vtk_unit,'(a)') '# vtk DataFile Version 2.0'
                    write(vtk_unit,'(a)') label
                    write(vtk_unit,'(a)') 'ASCII'
                    write(vtk_unit,'(a)') ' '
                endif

                bFileOpen=.true.
                nData=-1;
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

    subroutine vtk_close_file()
        if ( bFileOpen ) then
            close(vtk_unit)
            bFileOpen=.false.
        endif
    endsubroutine


    ! ------------------------------------------------------------------------- 
    ! --- POLYDATA STUFF
    ! ------------------------------------------------------------------------- 
    subroutine vtk_dataset_polydata(Points)
        real(ReKi), dimension(:,:),intent(in) :: Points  !< 3 x n
        integer :: i
        if ( bFileOpen ) then
            nPoints=size(Points,2)
            if (bBinary) then
                write(vtk_unit)'DATASET POLYDATA'//NL
                write(buffer,'(A,I0,A)') 'POINTS ', nPoints ,' double'
                write(vtk_unit)trim(buffer)//NL
                if (bChangeFrame)  then
                    do i=1,nPoints
                        write(vtk_unit)matmul(T_g2b,Points(1:3,i)-PO_g)
                    enddo
                else
                    do i=1,nPoints
                        write(vtk_unit)Points(1:3,i)
                    enddo
                endif
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A)') 'DATASET POLYDATA'
                write(vtk_unit,'(A,I0,A)') 'POINTS ', nPoints ,' double'
                if (bChangeFrame)  then
                    do i=1,nPoints
                        write(vtk_unit,'(3'//RFMT//')') matmul(T_g2b,Points(1:3,i)-PO_g)
                    enddo
                else
                    do i=1,nPoints
                        write(vtk_unit,'(3'//RFMT//')') Points(1:3,i)
                    enddo
                endif
                write(vtk_unit,*) ' '
            endif
        endif
    end subroutine


    subroutine vtk_lines(L)
        integer, dimension(:,:),intent(in) :: L    !< 2 x n

        integer :: i

        if ( bFileOpen ) then
            nData=size(L,2)
            if (bBinary) then
                write(buffer,'(A,I0,A,I0)')'LINES ',nData,' ',3*nData
                write(vtk_unit)trim(buffer)//NL
                do i=1,nData
                    write(vtk_unit)2,L(1:2,i)
                enddo
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A,I0,A,I0)')'LINES ',nData,' ',3*nData
                do i=1,nData
                    write(vtk_unit,'(3'//IFMT//')') 2, L(1:2,i)
                enddo
                write(vtk_unit,*) ' '
            endif
        endif
    end subroutine

    subroutine vtk_quad(Q)
        integer, dimension(:,:),intent(in) :: Q    !< 4 x n
        integer :: i
        if ( bFileOpen ) then
            nData=size(Q,2)
            if (bBinary) then
                write(buffer,'(A,I0,A,I0)')'POLYGONS ',nData,' ',5*nData
                write(vtk_unit)trim(buffer)//NL
                do i=1,nData
                    write(vtk_unit)4,Q(1:4,i)
                enddo
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A,I0,A,I0)') 'POLYGONS ', nData,' ',5*nData
                do i=1,nData
                    write(vtk_unit,'(5'//IFMT//')') 4, Q(1:4,i)
                enddo
                write(vtk_unit,*) ' '
            endif
        endif
    end subroutine

    ! ------------------------------------------------------------------------- 
    ! --- RECTILINEAR
    ! ------------------------------------------------------------------------- 
    subroutine vtk_dataset_rectilinear(v1,v2,v3)
        real(ReKi), dimension(:),intent(in) :: v1,v2,v3  !<  n

        if ( bFileOpen ) then
            nPoints=size(v1)*size(v2)*size(v3)
            if (bBinary) then
                write(vtk_unit) 'DATASET RECTILINEAR_GRID'//NL
                write(buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', size(v1),' ',size(v2),' ',size(v3)
                write(vtk_unit) trim(buffer)//NL
                write(buffer,'(A,I0,A)') 'X_COORDINATES ', size(v1), ' double'
                write(vtk_unit) trim(buffer)//NL
                write(vtk_unit)v1
                write(vtk_unit)NL
                write(buffer,'(A,I0,A)') 'Y_COORDINATES ', size(v2), ' double'
                write(vtk_unit) trim(buffer)//NL
                write(vtk_unit)v2
                write(vtk_unit)NL
                write(buffer,'(A,I0,A)') 'Z_COORDINATES ', size(v3), ' double'
                write(vtk_unit) trim(buffer)//NL
                write(vtk_unit)v3
                !write(vtk_unit)NL
            else
                write(vtk_unit,'(A)') 'DATASET RECTILINEAR_GRID'
                write(vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', size(v1),' ',size(v2),' ',size(v3)
                write(vtk_unit,'(A,I0,A)') 'X_COORDINATES ', size(v1), ' double'
                write(vtk_unit,'('//RFMT//')') v1
                write(vtk_unit,'(A,I0,A)') 'Y_COORDINATES ', size(v2), ' double'
                write(vtk_unit,'('//RFMT//')') v2
                write(vtk_unit,'(A,I0,A)') 'Z_COORDINATES ', size(v3), ' double'
                write(vtk_unit,'('//RFMT//')') v3
                write(vtk_unit,*) ' '
            endif
        endif
    end subroutine

    ! ------------------------------------------------------------------------- 
    ! --- STRUCTURED GRID (Points dumped without for loop since memory is in proper order)
    ! ------------------------------------------------------------------------- 
    !> Subroutine using flat data as input (not in natural order)
    subroutine vtk_dataset_structured_grid_flat(D,n1,n2,n3)
        integer , intent(in) :: n1,n2,n3
        real(ReKi), dimension(:,:),intent(in)::D
        if ( bFileOpen ) then
            nPoints=n1*n2*n3
            if (bBinary) then
                write(vtk_unit) 'DATASET STRUCTURED_GRID'//NL
                write(buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(vtk_unit) trim(buffer)//NL
                write(buffer,'(A,I0,A)') 'POINTS ', nPoints, ' double'
                write(vtk_unit) trim(buffer)//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A)') 'DATASET STRUCTURED_GRID'
                write(vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(vtk_unit,'(A,I0,A)') 'POINTS ', nPoints, ' double'
                write(vtk_unit,'(3'//RFMT//')')D
                write(vtk_unit,*) ' '
            endif
        endif
    end subroutine

    !> Using Grid data 4d as input
    subroutine vtk_dataset_structured_grid_grid(D,n1,n2,n3)
        integer , intent(in) :: n1,n2,n3
        real(ReKi), dimension(:,:,:,:),intent(in)::D

        if ( bFileOpen ) then
            nPoints=n1*n2*n3
            if (bBinary) then
                write(vtk_unit) 'DATASET STRUCTURED_GRID'//NL
                write(buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(vtk_unit) trim(buffer)//NL
                write(buffer,'(A,I0,A)') 'POINTS ', nPoints, ' double'
                write(vtk_unit) trim(buffer)//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A)') 'DATASET STRUCTURED_GRID'
                write(vtk_unit,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(vtk_unit,'(A,I0,A)') 'POINTS ', nPoints, ' double'
                write(vtk_unit,'(3'//RFMT//')')D
                write(vtk_unit,*) ' '
            endif
        endif
    end subroutine
    
    
    
    ! ------------------------------------------------------------------------- 
    ! --- POINT DATA
    ! ------------------------------------------------------------------------- 
    subroutine vtk_point_data_init()
        if ( bFileOpen ) then
            if(bBinary) then
                write(buffer,'(A,I0)')'POINT_DATA ',nPoints
                write(vtk_unit)trim(buffer)//NL
            else
                write(vtk_unit,'(A,I0)') 'POINT_DATA ', nPoints
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_flat(D,sname)
        real(ReKi), dimension(:),intent(in)::D
        character(len=*),intent(in) ::sname

        if ( bFileOpen ) then
            if (bBinary) then
                write(vtk_unit)'SCALARS '//trim(sname)//' double'//NL
                write(vtk_unit)'LOOKUP_TABLE default'//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_grid(D,sname)
        real(ReKi), dimension(:,:,:,:),intent(in)::D
        character(len=*),intent(in) ::sname

        if ( bFileOpen ) then
            if (bBinary) then
                write(vtk_unit)'SCALARS '//trim(sname)//' double'//NL
                write(vtk_unit)'LOOKUP_TABLE default'//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_grid2D(D,sname)
        real(ReKi), dimension(:,:,:),intent(in)::D
        character(len=*),intent(in) ::sname

        if ( bFileOpen ) then
            if (bBinary) then
                write(vtk_unit)'SCALARS '//trim(sname)//' double'//NL
                write(vtk_unit)'LOOKUP_TABLE default'//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    !>
    subroutine vtk_point_data_vector_flat(D,sname)
        real(ReKi), dimension(:,:),intent(in) :: D  !< 3 x n
        character(len=*),intent(in) ::sname
        if ( bFileOpen ) then
            if (bBinary) then
                write(vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine
    !>
    subroutine vtk_point_data_vector_grid(D,sname)
        real(ReKi), dimension(:,:,:,:),intent(in) :: D  !< 3 x n
        character(len=*),intent(in) ::sname
        if ( bFileOpen ) then
            if (bBinary) then
                write(vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine
    !>
    subroutine vtk_point_data_vector_grid2D(D,sname)
        real(ReKi), dimension(:,:,:),intent(in) :: D  !< 
        character(len=*),intent(in) ::sname
        if ( bFileOpen ) then
            if (bBinary) then
                write(vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine


    ! ------------------------------------------------------------------------- 
    ! --- CELL DATA
    ! ------------------------------------------------------------------------- 
    subroutine vtk_cell_data_init()
        if ( bFileOpen ) then
            if (bBinary) then
                write(buffer,'(A,I0)')'CELL_DATA ',nData
                write(vtk_unit)trim(buffer)//NL
            else
                write(vtk_unit,'(A,I0)') 'CELL_DATA ', nData
            endif
        endif
    end subroutine

    subroutine vtk_cell_data_scalar(D,sname)
        real(ReKi), dimension(:),intent(in)::D
        character(len=*),intent(in) ::sname

        if ( bFileOpen ) then
            if (bBinary) then
                write(vtk_unit)'SCALARS '//trim(sname)//' double 1'//NL
                write(vtk_unit)'LOOKUP_TABLE default'//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,fmt='(A,A,A)') 'SCALARS ', sname, ' double'
                write(vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine
    
    
    subroutine vtk_cell_data_vector(D,sname)
        real(ReKi), dimension(:,:),intent(in) :: D  !< 3 x n
        character(len=*),intent(in) ::sname
        if ( bFileOpen ) then
            if (bBinary) then
                write(vtk_unit)'VECTORS '//trim(sname)//' double'//NL
                write(vtk_unit)D
                write(vtk_unit)NL
            else
                write(vtk_unit,'(A,A,A)') 'VECTORS ', sname, ' double'
                write(vtk_unit,'(3'//RFMT//')')D
            endif
        endif
    end subroutine

    ! --------------------------------------------------------------------------------}
    ! --- VTK Tools 
    ! --------------------------------------------------------------------------------{
    !> Exports a Plane From a mesh
    subroutine export_plane_grid3d(fname,v1,v2,v3,Values)
        character(len=*),intent(in)             :: fname
        real(ReKi),dimension(:), intent(in)       :: v1,v2,v3
        real(ReKi),dimension(:,:,:,:), intent(in) :: Values
        !  Variables
        integer :: nD

        ! Writting
        if ( vtk_new_ascii_file(trim(fname),'plane', .false.)) then
            nD=size(Values,1)
            call vtk_dataset_rectilinear(v1,v2,v3)
            ! Output as a structured grid, No need to reorder
            call vtk_point_data_init()
            ! Could be a function of nDim, be careful
            if(nD==3) then
                call vtk_point_data_vector(Values(1:3,:,:,:),'Velocity') ! Label...
            endif

            call vtk_close_file()
        endif ! file opening
    end subroutine
    
    !> Exports a Plane From a mesh
    subroutine export_plane_grid2d(fname,v1,v2,v3,Values)
        character(len=*),intent(in)             :: fname
        real(ReKi),dimension(:), intent(in)       :: v1,v2,v3
        real(ReKi),dimension(:,:,:), intent(in) :: Values
        !  Variables
        integer :: nD

        ! Writting
        if ( vtk_new_ascii_file(trim(fname),'plane', .false.)) then
            nD=size(Values,1)
            call vtk_dataset_rectilinear(v1,v2,v3)
            ! Output as a structured grid, No need to reorder
            call vtk_point_data_init()
            ! Could be a function of nDim, be careful
            if(nD==3) then
                call vtk_point_data_vector(Values(1:3,:,:),'Velocity') ! Label...
            endif

            call vtk_close_file()
        endif ! file opening
    end subroutine
end module VTK
