!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of Ambient Wind and Array Effects model for FAST.Farm.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE AWAE_IO
 
   use NWTC_Library
   use AWAE_Types

   
   implicit none
   
   type(ProgDesc), parameter  :: AWAE_Ver = ProgDesc( 'AWAE', 'v00.01.00', '25-Jan-2017' )
   character(*),   parameter  :: AWAE_Nickname = 'AWAE'
    
   public :: AWAE_IO_InitGridInfo
   public :: ReadLowResWindFile
   
   contains

!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine ScanDir(dir, listingName)
   character(*), intent(in ) :: dir
   character(*), intent(in ) :: listingName
#if _WIN32
   call system('dir "'//trim(dir)//'" /B /A:-D-S-H > '//trim(listingName))
#else
   call system('ls '//trim(dir)//' > '//trim(listingName))
#endif
end subroutine ScanDir
   
subroutine ReadLowResWindFileHeaders(p, errStat, errMsg)
   type(AWAE_ParameterType),       intent(in   )  :: p            !< Parameters
   integer(IntKi),                 intent(  out)  :: errStat      !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg       !< Error message if errStat /= ErrID_None
   
   real :: r
   integer :: Un, i,reason,nFiles
   character(LEN=100), dimension(:), allocatable :: fileNames

   errStat = ErrID_None
   errMsg  = ""

   ! get the files
   call ScanDir(trim(p%WindFilePath)//'\Low\','dirContents.txt')
   CALL GetNewUnit( Un, ErrStat, ErrMsg )
   CALL OpenFOutFile ( Un, 'dirContents.txt', ErrStat, ErrMsg )

   !how many
   nFiles = 0
   do
      read(Un,FMT='(a)',iostat=reason) r
      if (reason/=0) EXIT
      nFiles = nFiles+1
   end do
  
   allocate(fileNames(nFiles))
   rewind(Un)
   do i = 1,nFiles
      read(Un,'(a)') fileNames(i)
      
   end do
   close(Un)
   
!==============================================================================
    
end subroutine ReadLowResWindFileHeaders  
   

subroutine WriteDisWindFiles( n, p, y, m, errStat, errMsg )
   integer(IntKi),             intent(in   ) :: n            !<  Low-resolution time step increment
   type(AWAE_ParameterType),   intent(in   ) :: p            !< Parameters
   type(AWAE_OutputType),   intent(in   ) :: y            !< Outputs
   type(AWAE_MiscVarType),     intent(inout) :: m            !< Misc/optimization variables
   integer(IntKi),             intent(  out) :: errStat      !< Error status of the operation
   character(*),               intent(  out) :: errMsg       !< Error message if errStat /= ErrID_None
   
   CHARACTER(1024)                         :: FileName
   INTEGER(IntKi)                          :: Un                   ! unit number of opened file   
   INTEGER(IntKi)                          :: ErrStat2                        ! Temporary Error status
   CHARACTER(ErrMsgLen)                    :: ErrMsg2                         ! Temporary Error message
   CHARACTER(*),   PARAMETER               :: RoutineName = 'WriteDisWindFiles'
   INTEGER(IntKi)                          :: nt, n_hl, n_high 
   
   
   FileName = trim(p%WindFilePath)//trim(PathSep)//"Low"//trim(PathSep)//"Dis.t"//trim(num2lstr(n))//".vtk"
   call WrVTK_SP_header( FileName, "Low resolution disturbed wind for low-resolution timestep "//trim(num2lstr(n)), Un, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   call WrVTK_SP_vectors3D( Un, "LowDis", (/p%nX_low,p%nY_low,p%nZ_low/), (/p%X0_low,p%Y0_low,p%Z0_low/), (/p%dX_low,p%dY_low,p%dZ_low/), m%Vdist_low, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
    
   do nt= 1,p%NumTurbines
      do n_hl = 1,p%n_high_low
         n_high = n_hl + p%n_high_low*n
      FileName = trim(p%WindFilePath)//trim(PathSep)//"High"//trim(PathSep)//"Dis.t"//trim(num2lstr(n_high))//".vtk"
      call WrVTK_SP_header( FileName, "High resolution disturbed wind for high-resolution timestep "//trim(num2lstr(n_high)), Un, errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
      call WrVTK_SP_vectors3D( Un, "HighDis", (/p%nX_high,p%nY_high,p%nZ_high/), (/p%X0_high,p%Y0_high,p%Z0_high/), (/p%dX_high,p%dY_high,p%dZ_high/), y%Vdist_high(:,:,:,:,n_hl,nt), errStat2, errMsg2 )
         call SetErrStat(ErrStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
      end do   
   end do
   
end subroutine WriteDisWindFiles


!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine ReadLowResWindFile(n, p, Vamb_Low, errStat, errMsg)
   integer(IntKi),                 intent(in   )  :: n            !< Current simulation timestep increment (zero-based)
   type(AWAE_ParameterType),       intent(in   )  :: p            !< Parameters
   real(ReKi),                     intent(inout)  :: Vamb_Low(:,0:,0:,0:)         !< Array which will contain the low resolution grid of ambient wind velocities
   integer(IntKi),                 intent(  out)  :: errStat      !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg       !< Error message if errStat /= ErrID_None
  
   integer(IntKi)           :: dims(3)              !  dimension of the 3D grid (nX,nY,nZ)
   real(ReKi)               :: origin(3)            !  the lower-left corner of the 3D grid (X0,Y0,Z0)
   real(ReKi)               :: gridSpacing(3)       !  spacing between grid points in each of the 3 directions (dX,dY,dZ)
   integer(IntKi)           :: Un                   !  unit number of opened file
   character(1024)          :: FileName             ! Name of output file     
   character(1024)          :: descr                ! Line describing the contents of the file
   character(1024)          :: vecLabel             ! descriptor of the vector data
   
   errStat = ErrID_None
   errMsg  = ""
  ! TODO: Try to skip this and just jump to the correct file location for the vector reads
 
   FileName = trim(p%WindFilePath)//trim(PathSep)//"Low"//trim(PathSep)//"Amb.t"//trim(num2lstr(n))//".vtk"
   call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat, ErrMsg ) 
      if (ErrStat >= AbortErrLev) return
   call ReadVTK_SP_vectors( FileName, Un, dims, Vamb_Low, ErrStat, ErrMsg )
      
   
!==============================================================================

   
end subroutine ReadLowResWindFile

!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine ReadHighResWindFile(nt, n, p, Vamb_high, errStat, errMsg)

   integer(IntKi),                 intent(in   )  :: nt
   integer(IntKi),                 intent(in   )  :: n                       !< high-res time increment
   type(AWAE_ParameterType),       intent(in   )  :: p            !< Parameters
   real(ReKi),                     intent(inout)  :: Vamb_high(:,0:,0:,0:)         !< Array which will contain the low resolution grid of ambient wind velocities
   integer(IntKi),                 intent(  out)  :: errStat      !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg       !< Error message if errStat /= ErrID_None
  
   
   integer(IntKi)           :: dims(3)              !  dimension of the 3D grid (nX,nY,nZ)
   real(ReKi)               :: origin(3)            !  the lower-left corner of the 3D grid (X0,Y0,Z0)
   real(ReKi)               :: gridSpacing(3)       !  spacing between grid points in each of the 3 directions (dX,dY,dZ)
   integer(IntKi)           :: Un                   !  unit number of opened file
   character(1024)          :: FileName             ! Name of output file     
   character(1024)          :: descr                ! Line describing the contents of the file
   character(1024)          :: vecLabel             ! descriptor of the vector data
   
   errStat = ErrID_None
   errMsg  = ""
   
! TODO: Try to skip this and just jump to the correct file location for the vector reads
   FileName = trim(p%WindFilePath)//trim(PathSep)//"HighT"//trim(num2lstr(nt))//trim(PathSep)//"Amb.t"//trim(num2lstr(n))//".vtk"
   call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat, ErrMsg ) 
      if (ErrStat >= AbortErrLev) return
   call ReadVTK_SP_vectors( FileName, Un, dims, Vamb_high, ErrStat, ErrMsg ) 
      
   
!==============================================================================
  
end subroutine ReadHighResWindFile

subroutine AWAE_IO_InitGridInfo(InitInp, p, InitOut, errStat, errMsg)

   type(AWAE_InitInputType),    intent(in   ) :: InitInp        !< Input data for initialization routine
   type(AWAE_ParameterType),    intent(inout) :: p              !< Parameters
   type(AWAE_InitOutputType),   intent(  out) :: InitOut        !< Output for initialization routine
   integer(IntKi),              intent(  out) :: errStat
   character(*),                intent(  out) :: errMsg

   integer(IntKi)                             :: NumGrid_low, NumGrid_high
   integer(IntKi)                             :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                       :: errMsg2       ! temporary error message 
   character(*), parameter                    :: RoutineName = 'AWAE_IO_InitGridInfo'
   real(ReKi)                                 :: X0_low, Y0_low, Z0_low, dX_low, dY_low, dZ_low, dt_low, dt_high
   integer(IntKi)                             :: nXYZ_low, nt, nx_low, ny_low, nz_low, nXYZ_high, nx_high, ny_high, nz_high
   integer(IntKi)                             :: dims(3)              ! dimension of the 3D grid (nX,nY,nZ)
   real(ReKi)                                 :: origin(3)            ! the lower-left corner of the 3D grid (X0,Y0,Z0)
   real(ReKi)                                 :: gridSpacing(3)       ! spacing between grid points in each of the 3 directions (dX,dY,dZ)
   character(1024)                            :: FileName             ! Name of output file     
   character(1024)                            :: descr                ! Line describing the contents of the file
   character(1024)                            :: vecLabel             ! descriptor of the vector data
   integer(IntKi)                             :: Un                   ! file unit
   errStat = ErrID_None
   errMsg  = ""

   
   

   FileName = trim(p%WindFilePath)//trim(PathSep)//"Low"//trim(PathSep)//"Amb.t0.vtk"
   Un = -1 ! Set to force closing of file on return
   call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat, ErrMsg )     
      if (ErrStat >= AbortErrLev) return     
   
   p%X0_Low           = origin(1)
   p%Y0_low           = origin(2)
   p%Z0_low           = origin(3) 
   p%nX_Low           = dims(1)
   p%nY_low           = dims(2)
   p%nZ_low           = dims(3) 
   p%dX_low           = gridSpacing(1)
   p%dY_low           = gridSpacing(2)
   p%dZ_low           = gridSpacing(3)
   
   InitOut%X0_Low     = origin(1)
   InitOut%Y0_low     = origin(2)
   InitOut%Z0_low     = origin(3) 
   InitOut%nX_Low     = dims(1)
   InitOut%nY_low     = dims(2)
   InitOut%nZ_low     = dims(3) 
   InitOut%dX_low     = gridSpacing(1)
   InitOut%dY_low     = gridSpacing(2)
   InitOut%dZ_low     = gridSpacing(3)
   
   NumGrid_low        = p%nX_Low*p%nY_Low*p%nZ_Low
   p%n_wind_min = 100
   
   p%n_wind_max = ceiling(30.0_ReKi*pi*(2.0_ReKi*p%r(p%NumRadii-1))**2*p%dt/(gridSpacing(1)*gridSpacing(2)*gridSpacing(3)))

      ! Grid runs from (X0_low, Y0_low, Z0_low) to (X0_low + (p%nX_Low-1)*dX_low, Y0_low+ (p%nY_Low-1)*dY_low, Z0_low+ (p%nZ_Low-1)*dZ_low)
      ! (0,0,0) to (180,180,180) 
     
  
   allocate( p%Grid_low(3,NumGrid_low),stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for Grid_low.', errStat, errMsg, RoutineName )
         return
      end if
      
   nXYZ_low = 0
   do nz_low=0, p%nZ_low-1 
      do ny_low=0, p%nY_low-1
         do nx_low=0, p%nX_low-1
            nXYZ_low = nXYZ_low + 1
            p%Grid_low(1,nXYZ_low) = origin(1) + nx_low*gridSpacing(1)
            p%Grid_low(2,nXYZ_low) = origin(2) + ny_low*gridSpacing(2)
            p%Grid_low(3,nXYZ_low) = origin(3) + nz_low*gridSpacing(3)            
         end do
      end do
   end do
   
    ! Parse a high res wind input file to gather the grid information
   
   
   FileName = trim(p%WindFilePath)//trim(PathSep)//"HighT1"//trim(PathSep)//"Amb.t0.vtk"  !TODO: Should the turbine numbers be padding with leading zero(es)?
    ! TODO: Error checking to see that all p%NumTurbines turbines use the same nX, nY, nZ for the high res grids
   Un = -1 ! Set to force closing of file on return
   call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat, ErrMsg ) 
      if (ErrStat >= AbortErrLev) return 
   
   p%nX_high          = dims(1)
   p%nY_high          = dims(2)
   p%nZ_high          = dims(3)
   NumGrid_high       = p%nX_high*p%nY_high*p%nZ_high
   
   allocate( p%Grid_high(3,NumGrid_high,p%NumTurbines ),stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for Grid_high.', errStat, errMsg, RoutineName )
         return
      end if
      
   allocate( InitOut%X0_high(p%NumTurbines), InitOut%Y0_high(p%NumTurbines), InitOut%Z0_high(p%NumTurbines), stat=errStat2)   
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for InitOut origin arrays.', errStat, errMsg, RoutineName )
   allocate( InitOut%dX_high(p%NumTurbines), InitOut%dY_high(p%NumTurbines), InitOut%dZ_high(p%NumTurbines), stat=errStat2)   
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for InitOut spatial increment arrays.', errStat, errMsg, RoutineName )
   allocate( p%X0_high(p%NumTurbines), p%Y0_high(p%NumTurbines), p%Z0_high(p%NumTurbines), stat=errStat2)   
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for p origin arrays.', errStat, errMsg, RoutineName )
   allocate( p%dX_high(p%NumTurbines), p%dY_high(p%NumTurbines), p%dZ_high(p%NumTurbines), stat=errStat2)   
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for p spatial increment arrays.', errStat, errMsg, RoutineName )
   if (ErrStat >= AbortErrLev) return
   
   do nt = 1, p%NumTurbines 
      FileName = trim(p%WindFilePath)//trim(PathSep)//"HighT"//trim(num2lstr(nt))//trim(PathSep)//"Amb.t0.vtk"
      Un = -1 ! Set to force closing of file on return
      call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat, ErrMsg ) 
         if (ErrStat >= AbortErrLev) return 
      InitOut%X0_high(nt) = origin(1)
      InitOut%Y0_high(nt) = origin(2)
      InitOut%Z0_high(nt) = origin(3)
      
      InitOut%dX_high(nt) = gridSpacing(1)
      InitOut%dY_high(nt) = gridSpacing(2)
      InitOut%dZ_high(nt) = gridSpacing(3)
      p%X0_high(nt) = origin(1)
      p%Y0_high(nt) = origin(2)
      p%Z0_high(nt) = origin(3)
      p%dX_high(nt) = gridSpacing(1)
      p%dY_high(nt) = gridSpacing(2)
      p%dZ_high(nt) = gridSpacing(3)
      
      nXYZ_high = 0
      do nz_high=0, p%nZ_high-1 
         do ny_high=0, p%nY_high-1
            do nx_high=0, p%nX_high-1
               nXYZ_high = nXYZ_high + 1
               p%Grid_high(1,nXYZ_high,nt) = InitOut%X0_high(nt) + nx_high*InitOut%dX_high(nt)
               p%Grid_high(2,nXYZ_high,nt) = InitOut%Y0_high(nt) + ny_high*InitOut%dY_high(nt)
               p%Grid_high(3,nXYZ_high,nt) = InitOut%Z0_high(nt) + nz_high*InitOut%dZ_high(nt)            
            end do
         end do
      end do
      
   end do
   
   InitOut%nx_high = p%nx_high
   InitOut%ny_high = p%ny_high
   InitOut%nz_high = p%nz_high
   
   
   !TODO:  Perform any error checking on InitOut and all wind input files here  : Review Plan.
   
! End simulated read of low and high res ambient wind files   
!==============================================================================
   
end subroutine AWAE_IO_InitGridInfo

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AWAE_PrintSum(  p, u, y, ErrStat, ErrMsg )
! This routine generates the summary file, which contains a summary of input file options.

      ! passed variables
   !TYPE(AWAE_InitInput),        INTENT(IN)  :: InputFileData                        ! Input-file data
   type(AWAE_ParameterType),    intent(in)  :: p                                    ! Parameters
   type(AWAE_InputType),        intent(in)  :: u                                    ! inputs 
   type(AWAE_OutputType),       intent(in)  :: y                                    ! outputs
   integer(IntKi),            intent(out) :: ErrStat
   character(*),              intent(out) :: ErrMsg


      ! Local variables.

   INTEGER(IntKi)               :: I                                               ! Index for the nodes.
   INTEGER(IntKi)               :: UnSu                                            ! I/O unit number for the summary output file

   CHARACTER(*), PARAMETER      :: FmtDat    = '(A,T35,1(:,F13.3))'                ! Format for outputting mass and modal data.
   CHARACTER(*), PARAMETER      :: FmtDatT   = '(A,T35,1(:,F13.8))'                ! Format for outputting time steps.

   CHARACTER(30)                :: OutPFmt                                         ! Format to print list of selected output channels to summary file
   CHARACTER(100)               :: Msg                                             ! temporary string for writing appropriate text to summary file

   ! Open the summary file and give it a heading.
      
   CALL GetNewUnit( UnSu, ErrStat, ErrMsg )
   CALL OpenFOutFile ( UnSu, TRIM( p%OutFileRoot )//'.sum', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN

 

   CLOSE(UnSu)

RETURN
END SUBROUTINE AWAE_PrintSum
!----------------------------------------------------------------------------------------------------------------------------------



END MODULE AWAE_IO

