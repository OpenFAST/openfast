!> Tools to read/write VTK files, for both
! - VTK ASCII format (vtk_* routines, and the vtk_misc type)
! - VTK ASCII format for structure points (SP)
! - VTK XML format (header and footer, the rest is found in ModMesh for now)
module VTK

   use Precision, only: IntKi, SiKi, ReKi
   use NWTC_Base, only: ErrID_None, ErrID_Fatal, AbortErrLev, ErrMsgLen, SetErrStat
   use NWTC_IO, only: GetNewUnit, NewLine, WrScr, ReadStr, OpenFOutFile
   use NWTC_IO, only: OpenFinpFile, ReadCom, Conv2UC

   implicit none

   character(*), parameter :: RFMT='E17.8E3'
   character(*), parameter :: IFMT='I7'

   ! Internal type to ensure the same options are used in between calls for the functions vtk_*
   TYPE, PUBLIC :: VTK_Misc
      integer :: vtk_unit
      logical :: bFileOpen=.false.

      integer :: nData=0;
      integer :: nPoints=0;

      logical :: bBinary = .false.
      character(len=255) :: buffer

      ! Reference Frame
      real(ReKi),dimension(3,3) :: T_g2b
      real(ReKi),dimension(3)   :: PO_g
   END TYPE VTK_Misc

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
         vtk_point_data_scalar_2D, &
         vtk_point_data_scalar_grid2D, &
         vtk_point_data_scalar_grid
   end interface
   interface vtk_cell_data_scalar; module procedure &
         vtk_cell_data_scalar_1d,&
         vtk_cell_data_scalar_2d 
   end interface

   private

   ! --- VTK ASCII routines
   public :: vtk_misc_init
   public :: set_vtk_binary_format
   public :: set_vtk_coordinate_transform
   public :: vtk_new_ascii_file
   public :: vtk_close_file
   public :: vtk_dataset_polydata
   public :: vtk_lines
   public :: vtk_quad
   public :: vtk_dataset_rectilinear
   public :: vtk_dataset_structured_points
   public :: vtk_dataset_structured_grid
   public :: vtk_point_data_init
   public :: vtk_point_data_scalar
   public :: vtk_point_data_vector
   public :: vtk_cell_data_init
   public :: vtk_cell_data_scalar
   public :: vtk_cell_data_vector
   public :: WrVTK_SP_header
   public :: WrVTK_SP_vectors3D
   public :: ReadVTK_SP_info
   public :: ReadVTK_SP_vectors
   ! --- VTK XML routines
   public :: WrVTK_header
   public :: WrVTK_footer
contains

!=======================================================================
!> This routine writes out the heading for an vtk xml file (associated footer generated in
!! nwtc_io::wrvtk_footer). It tries to open a text file for writing and returns the Unit number of the opened file.
   SUBROUTINE WrVTK_header( FileName, NumberOfPoints, NumberOfLines, NumberOfPolys, Un, ErrStat, ErrMsg ) 
   
      CHARACTER(*)    , INTENT(IN   )        :: FileName             !< Name of output file
      INTEGER(IntKi)  , INTENT(IN   )        :: NumberOfPoints       !< Number of points in this VTK file
      INTEGER(IntKi)  , INTENT(IN   )        :: NumberOfLines        !< Number of lines in this VTK file
      INTEGER(IntKi)  , INTENT(IN   )        :: NumberOfPolys        !< Number of polygons in this VTK file
      INTEGER(IntKi)  , INTENT(  OUT)        :: Un                   !< unit number of opened file
      INTEGER(IntKi)  , INTENT(  OUT)        :: ErrStat              !< error level/status of OpenFOutFile operation
      CHARACTER(*)    , INTENT(  OUT)        :: ErrMsg               !< message when error occurs
   
      !$OMP critical(fileopen)
      CALL GetNewUnit( Un, ErrStat, ErrMsg )      
      CALL OpenFOutFile ( Un, TRIM(FileName), ErrStat, ErrMsg )
      !$OMP end critical(fileopen)
         if (ErrStat >= AbortErrLev) return
      
      ! Write a VTP mesh file (Polygonal VTK file) with positions and polygons (surfaces)
      ! (note alignment of WRITE statements to make sure spaces are lined up in XML file)
      WRITE(Un,'(A)')         '<?xml version="1.0"?>'
      WRITE(Un,'(A)')         '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">' ! bjj note: we don't have binary data in this file, so byte_order shouldn't matter, right?
      WRITE(Un,'(A)')         '  <PolyData>'
      WRITE(Un,'(2(A,i7),A)') '    <Piece NumberOfPoints="', NumberOfPoints, '" NumberOfVerts="  0" NumberOfLines="', NumberOfLines, '"'
      WRITE(Un,'(A,i7,A)')    '           NumberOfStrips="  0" NumberOfPolys="',  NumberOfPolys, '">'
   
      RETURN
   END SUBROUTINE WrVTK_header
!=======================================================================
!> This routine writes out the footer for an vtk xml file (associated header generated  
!! in nwtc_io::wrvtk_header). It closes the file Un.
   SUBROUTINE WrVTK_footer( Un ) 
   
      INTEGER(IntKi)  , INTENT(IN   )        :: Un                   !< unit number of opened file
         
      WRITE(Un,'(A)')         '    </Piece>'
      WRITE(Un,'(A)')         '  </PolyData>'
      WRITE(Un,'(A)')         '</VTKFile>'
      CLOSE(Un)         
   
      RETURN
   END SUBROUTINE WrVTK_footer                
   
!=======================================================================
!> This routine reads the header for a vtk, ascii, structured_points dataset file,
!! including all the information about the structured points.  It tries to open a 
!! text file for reading and returns the Unit number of the opened file.
!! The caller is responsible for closing the file unit unless caller set Un = -1!
   SUBROUTINE ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat, ErrMsg ) 
   
      CHARACTER(*)    , INTENT(IN   )        :: FileName             !< Name of output file     
      CHARACTER(1024) , INTENT(  OUT)        :: descr                !< Line describing the contents of the file
      INTEGER(IntKi)  , INTENT(  OUT)        :: dims(3)              !< dimension of the 3D grid (nX,nY,nZ)
      REAL(ReKi)      , INTENT(  OUT)        :: origin(3)            !< the lower-left corner of the 3D grid (X0,Y0,Z0)
      REAL(ReKi)      , INTENT(  OUT)        :: gridSpacing(3)       !< spacing between grid points in each of the 3 directions (dX,dY,dZ)
      CHARACTER(1024) , INTENT(  OUT)        :: vecLabel
      INTEGER(IntKi)  , INTENT(INOUT)        :: Un                   !< unit number of opened file
      INTEGER(IntKi)  , INTENT(  OUT)        :: ErrStat              !< error level/status of OpenFOutFile operation
      CHARACTER(*)    , INTENT(  OUT)        :: ErrMsg               !< message when error occurs
   
      INTEGER(IntKi)              :: ErrStat2              ! local error level/status of OpenFOutFile operation
      CHARACTER(ErrMsgLen)        :: ErrMsg2               ! local message when error occurs
      CHARACTER(1024)             :: Dummy1, Dummy2
      CHARACTER(1024)             :: Line                  ! one line of the file
      CHARACTER(1024)             :: formatLbl
      CHARACTER(*), PARAMETER     :: RoutineName = 'ReadVTK_SP_info'
      INTEGER(IntKi)              :: sz, nPts, nArr, nums(2)
      LOGICAL                     :: closeOnReturn
      
      ErrStat = ErrID_None
      ErrMsg  = ''
      
      IF (Un == -1 ) THEN
         closeOnReturn = .TRUE.
      ELSE
         closeOnReturn = .FALSE.
      END IF
      
      !$OMP critical(fileopen)
      CALL GetNewUnit( Un, ErrStat, ErrMsg )      
      CALL OpenFInpFile ( Un, TRIM(FileName), ErrStat, ErrMsg )
      !$OMP end critical(fileopen)
         if (ErrStat >= AbortErrLev) return
      
       CALL ReadCom( Un, FileName, 'File header: Module Version (line 1)', ErrStat2, ErrMsg2, 0 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
      CALL ReadStr( Un, FileName, descr, 'descr', 'File Description line', ErrStat2, ErrMsg2, 0 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      formatLbl = ""   
      CALL ReadStr( Un, FileName, formatLbl, 'formatLbl', 'ASCII label', ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      call Conv2UC(formatLbl)
      if (INDEX(formatLbl, "ASCII" ) /= 1 ) THEN ! If this line doesn't contain the word ASCII, we have a bad file header
         CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: did not find ASCII label', ErrStat, ErrMsg, RoutineName )
      end if  
      Line = ""
      CALL ReadStr( Un, FileName, Line, "dataset", "DATASET STRUCTURED_POINTS", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "DATASET" ) /= 1 ) THEN ! If this line doesn't contain the word dataset, we have a bad file header
        CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: did not find DATASET label', ErrStat, ErrMsg, RoutineName )
      END IF 
      IF ( INDEX(Line, "STRUCTURED_POINTS" ) == 0 ) THEN ! If this line doesn't also contain the word structured_points, we have a bad file header
         CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: did not find STRUCTURED_POINTS label', ErrStat, ErrMsg, RoutineName )
      end if
        
         ! Dimensions
      Line = ""
      CALL ReadStr( Un, FileName, Line, "Dimensions", "DIMENSIONS data", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      Line = trim(Line)
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "DIMENSIONS" ) /= 1 ) THEN ! If this line doesn't contain the word dataset, we have a bad file header
         CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: did not find DIMENSIONS label', ErrStat, ErrMsg, RoutineName )
      ELSE
         sz = len(Line)
         Line = Line(12:sz)
         READ(Line,*, IOSTAT=ErrStat2)  dims
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Error reading "dims".', ErrStat, ErrMsg, RoutineName )
         end if
      END IF 
      
         ! Origin
      Line = ""
      CALL ReadStr( Un, FileName, Line, "Origin", "ORIGIN data", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      Line = trim(Line)
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "ORIGIN" ) /= 1 ) THEN ! If this line doesn't contain the word dataset, we have a bad file header
         CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: did not find ORIGIN label', ErrStat, ErrMsg, RoutineName )
      ELSE
         sz = len(Line)
         Line = Line(8:sz)
         READ(Line,*, IOSTAT=ErrStat2)  origin
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Error reading "origin".', ErrStat, ErrMsg, RoutineName )
         end if

      END IF 
      
         ! Spacing      
      Line = ""
      CALL ReadStr( Un, FileName, Line, "gridSpacing", "SPACING data", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      Line = trim(Line)
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "SPACING" ) /= 1 ) THEN ! If this line doesn't contain the word dataset, we have a bad file header
         CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: did not find SPACING label', ErrStat, ErrMsg, RoutineName )
      ELSE
         sz = len(Line)
         Line = Line(9:sz)
         READ(Line,*,IOSTAT=ErrStat2)  gridSpacing
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Error reading "gridSpacing".', ErrStat, ErrMsg, RoutineName )
         end if
         
      END IF 
      
         ! Point Data
      Line = ""
      CALL ReadStr( Un, FileName, Line, "Point_Data", "POINT_DATA", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      Line = trim(Line)
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "POINT_DATA" ) /= 1 ) THEN ! If this line doesn't contain the word dataset, we have a bad file header
         CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: did not find POINT_DATA label', ErrStat, ErrMsg, RoutineName )
      ELSE
         sz = len(Line)
         Line = Line(12:sz)
         READ(Line,*,IOSTAT=ErrStat2)  nPts
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Error reading "nPts".', ErrStat, ErrMsg, RoutineName )
         end if
         IF ( nPts /= ( dims(1)*dims(2)*dims(3) ) ) THEN ! Abort if DIMENSIONS AND POINT_DATA don't agree
            CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: POINT_DATA does not match DIMENSIONS', ErrStat, ErrMsg, RoutineName )
         END IF
      END IF 
      
         ! VECTOR or FIELD Label
      Line = ""
      CALL ReadStr( Un, FileName, Line, "VECTORS or FIELD", "VECTORS or FIELD label", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      Line = trim(Line)
      CALL Conv2UC( Line )
      IF ( ( INDEX(Line, "VECTORS" ) /= 1 ) .AND. ( INDEX(Line, "FIELD" ) /= 1 ) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: did not find VECTORS or FIELD label', ErrStat, ErrMsg, RoutineName )
      ELSE
         IF ( INDEX(Line, "FIELD" ) == 1 ) THEN ! Must be FIELD
            READ(Line,*,IOSTAT=ErrStat2) Dummy1, Dummy2, nArr
            if (ErrStat2 /= 0) then
                CALL SetErrStat( ErrID_Fatal, 'Error reading "nArr".', ErrStat, ErrMsg, RoutineName )
            ELSE IF ( nArr /= 1_IntKi ) THEN
                CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: FIELD label must have only 1 array', ErrStat, ErrMsg, RoutineName )
            END IF
            
            Line = ""
            CALL ReadStr( Un, FileName, Line, "Array", "Array definition", ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
            Line = trim(Line)
            Call Conv2UC( Line )
            sz = INDEX(Line, "FLOAT" )
            IF ( sz == 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Invalid FIELD datatype.  Must be set to float.', ErrStat, ErrMsg, RoutineName )
            ELSE        
               READ(Line,*,IOSTAT=ErrStat2) Dummy1, nums
               if (ErrStat2 /= 0) then
                  CALL SetErrStat( ErrID_Fatal, 'Error reading "nums".', ErrStat, ErrMsg, RoutineName )
               ELSEIF ( nums(1) /= 3_IntKi ) THEN                         ! Abort if we don't have 3-element vectors
                  CALL SetErrStat( ErrID_Fatal, 'Invalid FIELD datatype.  FIELD array must have 3 elements.', ErrStat, ErrMsg, RoutineName )
               ELSEIF ( nums(2) /= ( dims(1)*dims(2)*dims(3) ) ) THEN ! Abort if DIMENSIONS AND FIELD data don't agree
                  CALL SetErrStat( ErrID_Fatal, 'Invalid vtk structured_points file: FIELD array does not match DIMENSIONS', ErrStat, ErrMsg, RoutineName )
               END IF
            END IF
         ELSE                                    ! Must be VECTORS
            sz = INDEX(Line, "FLOAT" )
            IF ( sz == 0 ) THEN
               CALL SetErrStat( ErrID_Fatal, 'Invalid VECTORS datatype.  Must be set to float.', ErrStat, ErrMsg, RoutineName )
            ELSE        
               vecLabel = Line(9:sz-2)
            END IF
         END IF
      END IF
      
      IF ( (ErrStat >= AbortErrLev) .or. closeOnReturn ) THEN        
         close(Un)
         Un = -1
         RETURN
      END IF
      
      RETURN
   END SUBROUTINE ReadVTK_SP_info
   
!=======================================================================
!> This routine reads the vector data for a vtk, ascii, structured_points dataset file,
!! The Unit number of the  file is already assumed to be valid via a previous call to
!! ReadVTK_SP_info.  
   SUBROUTINE ReadVTK_SP_vectors( FileName, Un, dims, gridVals, ErrStat, ErrMsg ) 
   
      CHARACTER(*)    , INTENT(IN   )        :: FileName             !< Name of output file    
      INTEGER(IntKi)  , INTENT(IN   )        :: Un                   !< unit number of opened file
      INTEGER(IntKi)  , INTENT(IN   )        :: dims(3)              !< dimension of the 3D grid (nX,nY,nZ)
      REAL(SiKi)      , INTENT(  OUT)        :: gridVals(:,:,:,:)    !< 4D array of data, size (3,nX,nY,nZ), must be pre-allocated
      INTEGER(IntKi)  , INTENT(  OUT)        :: ErrStat              !< error level/status of OpenFOutFile operation
      CHARACTER(*)    , INTENT(  OUT)        :: ErrMsg               !< message when error occurs
      
      INTEGER                                :: ErrStat2
      
      ErrStat = ErrID_None
      ErrMsg  = ''
      
      READ(Un,*, IOSTAT=ErrStat2)  gridVals(1:3,1:dims(1),1:dims(2),1:dims(3))
      
      close(Un)
      if (ErrStat2 /= 0) then
         CALL SetErrStat( ErrID_Fatal, 'Invalid vtk file: '//trim(FileName)//'.', ErrStat, ErrMsg, 'ReadVTK_SP_vectors' )
      end if
      
   END SUBROUTINE ReadVTK_SP_vectors
   
!=======================================================================
!> This routine writes out the heading for an vtk, ascii, structured_points dataset file .
!! It tries to open a text file for writing and returns the Unit number of the opened file.
   SUBROUTINE WrVTK_SP_header( FileName, descr, Un, ErrStat, ErrMsg ) 
   
      CHARACTER(*)    , INTENT(IN   )        :: FileName             !< Name of output file     
      CHARACTER(*)    , INTENT(IN   )        :: descr                !< Line describing the contents of the file
      INTEGER(IntKi)  , INTENT(  OUT)        :: Un                   !< unit number of opened file
      INTEGER(IntKi)  , INTENT(  OUT)        :: ErrStat              !< error level/status of OpenFOutFile operation
      CHARACTER(*)    , INTENT(  OUT)        :: ErrMsg               !< message when error occurs
   
      !$OMP critical(fileopen)
      CALL GetNewUnit( Un, ErrStat, ErrMsg )      
      CALL OpenFOutFile ( Un, TRIM(FileName), ErrStat, ErrMsg )
      !$OMP end critical(fileopen)
         if (ErrStat >= AbortErrLev) return
      
      WRITE(Un,'(A)')  '# vtk DataFile Version 3.0'
      WRITE(Un,'(A)')  trim(descr)
      WRITE(Un,'(A)')  'ASCII'
      WRITE(Un,'(A)')  'DATASET STRUCTURED_POINTS'
      
      RETURN
   END SUBROUTINE WrVTK_SP_header
   
   
   SUBROUTINE WrVTK_SP_vectors3D( Un, dataDescr, dims, origin, gridSpacing, gridVals, ErrStat, ErrMsg ) 
   
      INTEGER(IntKi)  , INTENT(IN   )        :: Un                   !< unit number of previously opened file (via call to WrVTK_SP_header)
      CHARACTER(*)    , INTENT(IN   )        :: dataDescr            !< Short label describing the vector data
      INTEGER(IntKi)  , INTENT(IN   )        :: dims(3)              !< dimension of the 3D grid (nX,nY,nZ)
      REAL(ReKi)      , INTENT(IN   )        :: origin(3)            !< the lower-left corner of the 3D grid (X0,Y0,Z0)
      REAL(ReKi)      , INTENT(IN   )        :: gridSpacing(3)       !< spacing between grid points in each of the 3 directions (dX,dY,dZ)
      REAL(SiKi)      , INTENT(IN   )        :: gridVals(:,:,:,:)      !< 3D array of data, size (nX,nY,nZ)
      INTEGER(IntKi)  , INTENT(  OUT)        :: ErrStat              !< error level/status of OpenFOutFile operation
      CHARACTER(*)    , INTENT(  OUT)        :: ErrMsg               !< message when error occurs
 
      INTEGER(IntKi)                         :: nPts                 ! Total number of grid points 
      
      if ( .not. (Un > 0) ) then
         ErrStat = ErrID_Fatal
         ErrMsg  = 'WrVTK_SP_points: Invalid file unit, be sure to call WrVTK_SP_header prior to calling WrVTK_SP_points.'
         return
      end if
   
      ErrStat = ErrID_None
      ErrMsg  = ''
      nPts    = dims(1)*dims(2)*dims(3)
      
      ! Note: gridVals must be stored such that the left-most dimension is X and the right-most dimension is Z
      WRITE(Un,'(A,3(i5,1X))')    'DIMENSIONS ',  dims
      WRITE(Un,'(A,3(f10.2,1X))') 'ORIGIN '    ,  origin
      WRITE(Un,'(A,3(f10.2,1X))') 'SPACING '   ,  gridSpacing
      WRITE(Un,'(A,i15)')         'POINT_DATA ',  nPts
      WRITE(Un,'(A)')            'VECTORS '//trim(dataDescr)//' float'
      WRITE(Un,'(3(f10.2,1X))')   gridVals
      close(Un)
      RETURN
      
   END SUBROUTINE WrVTK_SP_vectors3D

   subroutine vtk_misc_init(mvtk)
      type(VTK_Misc),intent(inout) :: mvtk
      mvtk%vtk_unit = -1           !< VTK output unit [-]
      mvtk%bFileOpen = .false.     !< binary file is open [-]
      mvtk%bBinary = .false.       !< write binary files [-]
      mvtk%nData = 0               !< number of data lines [-]
      mvtk%nPoints = 0             !< number of points [-]
   end subroutine

    !>
    subroutine set_vtk_binary_format(bBin,mvtk)
        logical, intent(in)::bBin
        type(VTK_Misc),intent(inout) :: mvtk
        mvtk%bBinary=bBin
    end subroutine

    
    !> Save a coordinate transform
    ! ALL VTK Will be exported in this coordinate system!
    subroutine set_vtk_coordinate_transform(T_g2b_in,PO_g_in,mvtk)
        real(ReKi),dimension(3,3), intent(in) :: T_g2b_in
        real(ReKi),dimension(3)  , intent(in) :: PO_g_in
        type(VTK_Misc),intent(inout) :: mvtk
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
        type(VTK_Misc),intent(inout) :: mvtk
        !
        integer :: iostatvar
        logical :: b

        if (.not. mvtk%bFileOpen) then
            !$OMP critical(fileopen)
            CALL GetNewUnit( mvtk%vtk_unit )   
            if (mvtk%bBinary) then
                ! Fortran 2003 stream, otherwise intel fortran !
                !form='UNFORMATTED',access='SEQUENTIAL',action='WRITE',convert='BIG_ENDIAN',recordtype='STREAM',buffered='YES',
               !print*,'Not available for this compiler' !COMPAQ-COMPILER
               !STOP !COMPAQ-COMPILER
!bjj: CONVERT is non-standard, so maybe this should be part of Sys*.f90? Like OpenUnfInpBEFile()?
!eb : Commented out for now since it doesnt work anymore anyway (worked 5 years ago..). Adding binary would be great though!
                !open(unit = mvtk%vtk_unit,file= trim(adjustl(filename)),form='UNFORMATTED',access = 'stream',& !OTHER-COMPILER
                !    action = 'WRITE',convert= 'BIG_ENDIAN',iostat=iostatvar,status='replace') !OTHER-COMPILER
                call WrScr('Binary VTK output not available at the moment')
                STOP ! Temporary, but easiest for now.
            else
                open(mvtk%vtk_unit,file=trim(adjustl(filename)),iostat=iostatvar,action="write",status='replace')
            endif
            !$OMP end critical(fileopen)
            if (iostatvar == 0) then
                if (mvtk%bBinary) then
                    write(mvtk%vtk_unit)'# vtk DataFile Version 3.0'//NewLine
                    write(mvtk%vtk_unit)trim(label)//NewLine
                    write(mvtk%vtk_unit)'BINARY'//NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
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
        type(VTK_Misc),intent(inout) :: mvtk
        logical, intent(in) :: bladeFrame
        integer :: i
        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=size(Points,2)
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'DATASET POLYDATA'//NewLine
                write(mvtk%buffer,'(A,I0,A)') 'POINTS ', mvtk%nPoints ,' double'
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NewLine
                if (bladeFrame)  then
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit)matmul(mvtk%T_g2b,Points(1:3,i)-mvtk%PO_g)
                    enddo
                else
                    do i=1,mvtk%nPoints
                        write(mvtk%vtk_unit)Points(1:3,i)
                    enddo
                endif
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk

        integer :: i

        if ( mvtk%bFileOpen ) then
            mvtk%nData=size(L,2)
            if (mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0,A,I0)')'LINES ',mvtk%nData,' ',3*mvtk%nData
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NewLine
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit)2,L(1:2,i)
                enddo
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
        integer :: i
        if ( mvtk%bFileOpen ) then
            mvtk%nData=size(Q,2)
            if (mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0,A,I0)')'POLYGONS ',mvtk%nData,' ',5*mvtk%nData
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NewLine
                do i=1,mvtk%nData
                    write(mvtk%vtk_unit)4,Q(1:4,i)
                enddo
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=size(v1)*size(v2)*size(v3)
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET RECTILINEAR_GRID'//NewLine
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', size(v1),' ',size(v2),' ',size(v3)
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%buffer,'(A,I0,A)') 'X_COORDINATES ', size(v1), ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%vtk_unit)v1
                write(mvtk%vtk_unit)NewLine
                write(mvtk%buffer,'(A,I0,A)') 'Y_COORDINATES ', size(v2), ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%vtk_unit)v2
                write(mvtk%vtk_unit)NewLine
                write(mvtk%buffer,'(A,I0,A)') 'Z_COORDINATES ', size(v3), ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%vtk_unit)v3
                !write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=n(1)*n(2)*n(3)
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET STRUCTURED_POINTS'//NewLine
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ',n(1),' ',n(2),' ',n(3)
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%buffer,'(A,3F16.8)') 'ORIGIN ', x0
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%buffer,'(A,3F16.8)') 'SPACING ', dx
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=n1*n2*n3
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET STRUCTURED_GRID'//NewLine
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%buffer,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            mvtk%nPoints=n1*n2*n3
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit) 'DATASET STRUCTURED_GRID'//NewLine
                write(mvtk%buffer,'(A,I0,A,I0,A,I0)') 'DIMENSIONS ', n1,' ',n2,' ',n3
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%buffer,'(A,I0,A)') 'POINTS ', mvtk%nPoints, ' double'
                write(mvtk%vtk_unit) trim(mvtk%buffer)//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if(mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0)')'POINT_DATA ',mvtk%nPoints
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NewLine
            else
                write(mvtk%vtk_unit,'(A,I0)') 'POINT_DATA ', mvtk%nPoints
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_flat(D,sname,mvtk)
        real(ReKi), dimension(:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NewLine
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
            else
                write(mvtk%vtk_unit,'(A,A,A)') 'SCALARS ', sname, ' double'
                write(mvtk%vtk_unit,'(A)') 'LOOKUP_TABLE default'
                write(mvtk%vtk_unit,'(1'//RFMT//')')D
            endif
        endif
    end subroutine

    subroutine vtk_point_data_scalar_2D(D,sname,mvtk)
        real(ReKi), dimension(:,:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NewLine
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NewLine
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double'//NewLine
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%buffer,'(A,I0)')'CELL_DATA ',mvtk%nData
                write(mvtk%vtk_unit)trim(mvtk%buffer)//NewLine
            else
                write(mvtk%vtk_unit,'(A,I0)') 'CELL_DATA ', mvtk%nData
            endif
        endif
    end subroutine

    subroutine vtk_cell_data_scalar_1d(D,sname,mvtk)
        real(ReKi), dimension(:),intent(in)::D
        character(len=*),intent(in) ::sname
        type(VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double 1'//NewLine
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk

        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'SCALARS '//trim(sname)//' double 1'//NewLine
                write(mvtk%vtk_unit)'LOOKUP_TABLE default'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
        if ( mvtk%bFileOpen ) then
            if (mvtk%bBinary) then
                write(mvtk%vtk_unit)'VECTORS '//trim(sname)//' double'//NewLine
                write(mvtk%vtk_unit)D
                write(mvtk%vtk_unit)NewLine
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
        type(VTK_Misc),intent(inout) :: mvtk
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
        type(VTK_Misc),intent(inout) :: mvtk
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
end module VTK
