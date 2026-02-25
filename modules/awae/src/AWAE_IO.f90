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
   use VTK
   use AWAE_Types
   use iso_c_binding, only: c_char, c_int, c_double, c_float, c_null_char
   use amrex_utils
   
   implicit none
   
   type(ProgDesc), parameter  :: AWAE_Ver = ProgDesc( 'AWAE', '', '' )
   character(*),   parameter  :: AWAE_Nickname = 'AWAE'
    
   public :: AWAE_IO_InitGridInfo
   public :: ReadLowResWindVTK, ReadWindAMReX

   interface
      subroutine ReadVTK_inflow_info(FileName, Desc, dims, origin, gridSpacing, vecLabel, values, read_values, err_stat, err_msg) BIND(C,name='ReadVTK_inflow_info')     
         use iso_c_binding, only: c_char, c_int, c_double, c_float, c_null_char        
         implicit none
         character(kind=c_char), intent(in)  :: FileName(*)
         character(kind=c_char), intent(out) :: Desc(1024)
         integer(c_int), intent(out)         :: dims(3)
         real(c_double), intent(out)         :: origin(3)
         real(c_double), intent(out)         :: gridSpacing(3)
         character(kind=c_char), intent(out) :: vecLabel(1024)
         real(c_float), intent(out)          :: values(*)
         integer(c_int), intent(in)          :: read_values
         integer(c_int), intent(out)         :: err_stat
         character(kind=c_char), intent(out) :: err_msg(1024)
      end subroutine
   end interface
   
   contains

subroutine WriteDisWindFiles( n, WrDisSkp1, p, y, m, errStat, errMsg )
   integer(IntKi),             intent(in   ) :: n            !< Low-resolution time step increment
   integer(IntKi),             intent(in   ) :: WrDisSkp1    !< Number of low resolution time step increments per one output increment
   type(AWAE_ParameterType),   intent(in   ) :: p            !< Parameters
   type(AWAE_OutputType),      intent(in   ) :: y            !< Outputs
   type(AWAE_MiscVarType),     intent(inout) :: m            !< Misc/optimization variables
   integer(IntKi),             intent(  out) :: errStat      !< Error status of the operation
   character(*),               intent(  out) :: errMsg       !< Error message if errStat /= ErrID_None
   
   CHARACTER(*),PARAMETER         :: RoutineName = 'WriteDisWindFiles'
   INTEGER(IntKi)                 :: ErrStat2                ! Temporary Error status
   CHARACTER(ErrMsgLen)           :: ErrMsg2                 ! Temporary Error message
   CHARACTER(1024)                :: FileName
   INTEGER(IntKi)                 :: Un                      ! unit number of opened file
   INTEGER(IntKi)                 :: nt, n_out 
   REAL(ReKi)                     :: t_out
   character(p%VTK_tWidth)        :: Tstr  ! string for current VTK write-out step (padded with zeros)
   
   n_out = n/WrDisSkp1
   t_out = n*p%DT_low

   ! TimeStamp
   write(Tstr, '(i' // trim(Num2LStr(p%VTK_tWidth)) //'.'// trim(Num2LStr(p%VTK_tWidth)) // ')') n_out ! TODO use n instead..

   FileName = trim(p%OutFileVTKRoot)//".Low.Dis."//trim(Tstr)//".vtk"
   call WrVTK_SP_header( FileName, "Low resolution disturbed wind for time = "//trim(num2lstr(t_out))//" seconds.", Un, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   call WrVTK_SP_vectors3D( Un, "Velocity", p%LowRes%nXYZ, p%LowRes%oXYZ, p%LowRes%dXYZ, m%Vdist_low, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
    
   do nt = 1, p%NumTurbines
       ! We are only writing out the first of the high res data for a given low res time step
      
      FileName = trim(p%OutFileVTKRoot)//".HighT"//trim(num2lstr(nt))//".Dis."//trim(Tstr)//".vtk"
      call WrVTK_SP_header( FileName, "High resolution disturbed wind for time = "//trim(num2lstr(t_out))//" seconds.", Un, errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
      call WrVTK_SP_vectors3D( Un, "Velocity", p%HighRes(nt)%nXYZ, p%HighRes(nt)%oXYZ, p%HighRes(nt)%dXYZ, y%Vdist_high(nt)%data(:,:,:,:,0), errStat2, errMsg2 )
         call SetErrStat(ErrStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
       
   end do


end subroutine WriteDisWindFiles

!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine read the low res wind file (VTK) at a given time step `n`
subroutine ReadLowResWindVTK(n, p, Vamb_Low, errStat, errMsg)
   integer(IntKi),                 intent(in   )  :: n            !< Current simulation timestep increment (zero-based)
   type(AWAE_ParameterType),       intent(in   )  :: p            !< Parameters
   real(SiKi), contiguous,         intent(inout)  :: Vamb_Low(:,0:,0:,0:)         !< Array which will contain the low resolution grid of ambient wind velocities
   integer(IntKi),                 intent(  out)  :: errStat      !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg       !< Error message if errStat /= ErrID_None
  
   integer(IntKi)           :: dims(3)              ! Dimension of the 3D grid (nX,nY,nZ)
   real(R8Ki)               :: origin(3)            ! The lower-left corner of the 3D grid (X0,Y0,Z0)
   real(R8Ki)               :: gridSpacing(3)       ! Spacing between grid points in each of the 3 directions (dX,dY,dZ)
   character(kind=c_char)   :: FileName(2048)       ! Name of output file     
   character(kind=c_char)   :: desc(1024)           ! Line describing the contents of the file
   character(kind=c_char)   :: vecLabel(1024)       ! descriptor of the vector data
   
   errStat = ErrID_None
   errMsg  = ""
  
   FileName = transfer(trim(p%WindFilePath)//trim(PathSep)//"Low"//trim(PathSep)//"Amb.t"//trim(Num2LStr(n))//".vtk"//c_null_char, FileName)
   call ReadVTK_inflow_info(FileName, desc, dims, origin, gridSpacing, vecLabel, Vamb_Low, 1, ErrStat, ErrMsg)
   if (ErrStat /= ErrID_None) ErrMsg = "ReadLowResWindVTK:"//trim(ErrMsg)

end subroutine ReadLowResWindVTK

!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine read the high res wind file (VTK) at a given time step `n`
subroutine ReadHighResWindVTK(nt, n, p, Vamb_high, errStat, errMsg)

   integer(IntKi),                 intent(in   )  :: nt
   integer(IntKi),                 intent(in   )  :: n                       !< high-res time increment
   type(AWAE_ParameterType),       intent(in   )  :: p            !< Parameters
   real(SiKi), contiguous,         intent(inout)  :: Vamb_high(:,0:,0:,0:)         !< Array which will contain the low resolution grid of ambient wind velocities
   integer(IntKi),                 intent(  out)  :: errStat      !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg       !< Error message if errStat /= ErrID_None
  
   integer(IntKi)           :: dims(3)              ! Dimension of the 3D grid (nX,nY,nZ)
   real(R8Ki)               :: origin(3)            ! The lower-left corner of the 3D grid (X0,Y0,Z0)
   real(R8Ki)               :: gridSpacing(3)       ! Spacing between grid points in each of the 3 directions (dX,dY,dZ)
   character(kind=c_char)   :: FileName(2048)       ! Name of output file
   character(kind=c_char)   :: desc(1024)           ! Line describing the contents of the file
   character(kind=c_char)   :: vecLabel(1024)       ! descriptor of the vector data
   
   errStat = ErrID_None
   errMsg  = ""
   
   FileName = transfer(trim(p%WindFilePath)//trim(PathSep)//"HighT"//trim(num2lstr(nt))//trim(PathSep)//"Amb.t"//trim(num2lstr(n))//".vtk"//c_null_char, FileName)
   call ReadVTK_inflow_info(FileName, desc, dims, origin, gridSpacing, vecLabel, Vamb_high, 1, ErrStat, ErrMsg)
   if (ErrStat /= ErrID_None) ErrMsg = "ReadHighResWindVTK:"//trim(ErrMsg)

end subroutine ReadHighResWindVTK


!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine read the AMReX at a given time step `n`
subroutine ReadWindAMReX(sv, n, p, Vamb, ErrStat, ErrMsg)
   use amrex_utils

   integer(IntKi),           intent(in   )  :: sv              ! Sub-volume {0: low-res, 1+: high-res turbine}
   integer(IntKi),           intent(in   )  :: n               !< time increment
   type(AWAE_ParameterType), intent(in   )  :: p               !< Parameters
   real(SiKi), contiguous,   intent(inout)  :: Vamb(:,:,:,:)   !< Array which will contain the grid of ambient wind velocities
   integer(IntKi),           intent(  out)  :: ErrStat         !< Error status of the operation
   character(*),             intent(  out)  :: ErrMsg          !< Error message if errStat /= ErrID_None
  
   character(len=2048) :: FileName       ! Name of output file
   character(len=12)   :: DirIndex       ! Directory index suffix 
   integer(IntKi)      :: i

   ! If sub-volume is 0 then this is low-resolution file
   if (sv == 0) then
      write(DirIndex,'(i'//trim(Num2LStr(p%DirIndexLen))//')') p%DirStartNum + p%DirIndexDeltaLow * n
   else
      write(DirIndex,'(i'//trim(Num2LStr(p%DirIndexLen))//')') p%DirStartNum + p%DirIndexDeltaHigh * n
   end if

   ! Prepend zeros in front of index number
   do i = 1, p%DirIndexLen
      if (DirIndex(i:i) /= " ") exit
      DirIndex(i:i) = '0'
   end do

   FileName = trim(p%WindFilePath)//"_"//trim(num2lstr(sv))//"_"//DirIndex(1:p%DirIndexLen)
   call amrex_read_data(FileName, Vamb, ErrStat, ErrMsg)

end subroutine

!----------------------------------------------------------------------------------------------------------------------------------   
!> Flat array of Cartesian point coordinates
!! Grid runs from (X0, Y0, Z0) to (X0 + (p%nX-1)*dX, Y0+ (p%nY-1)*dY, Z0+ (p%nZ-1)*dZ)
subroutine flatCartGridCoordinates(Origin, n, d, GridP)
   Real(ReKi)    , intent(in ) :: Origin(3)  !< 
   integer(IntKi), intent(in ) :: n(3)       !< dimension    nx, ny, nz
   Real(ReKi)    , intent(in ) :: d(3)       !< grid spacing dx, dy, dz
   real(ReKi),     intent(out) :: GridP(:,:) !< Grid points, flatten 3 x nTot
   integer(IntKi) :: iXYZ 
   integer :: ix, iy, iz
   iXYZ = 0
   do iz=0, n(3)-1 
      do iy=0, n(2)-1
         do ix=0, n(1)-1
            iXYZ = iXYZ + 1
            GridP(1,iXYZ) = origin(1) + ix*d(1)
            GridP(2,iXYZ) = origin(2) + iy*d(2)
            GridP(3,iXYZ) = origin(3) + iz*d(3)            
         end do
      end do
   end do
end subroutine flatCartGridCoordinates
!----------------------------------------------------------------------------------------------------------------------------------   
!> Initialize low and high res grid from VTK or InflowWind 
!! Set grid points, perform sanity checks
subroutine AWAE_IO_InitGridInfo(InitInp, p, InitOut, errStat, errMsg)

   type(AWAE_InitInputType),    intent(in   ) :: InitInp        !< Input data for initialization routine
   type(AWAE_ParameterType),    intent(inout) :: p              !< Parameters
   type(AWAE_InitOutputType),   intent(  out) :: InitOut        !< Output for initialization routine
   integer(IntKi),              intent(  out) :: errStat
   character(*),                intent(  out) :: errMsg

   character(*), parameter                    :: RoutineName = 'AWAE_IO_InitGridInfo'
   integer(IntKi)                             :: ErrStat2      ! temporary error status of the operation
   character(ErrMsgLen)                       :: ErrMsg2       ! temporary error message
   integer(IntKi)                             :: i, j, k
   integer(IntKi)                             :: nXYZ_low, nx_low, ny_low, nz_low 
   integer(IntKi)                             :: nXYZ_high, nx_high, ny_high, nz_high
   integer(IntKi)                             :: dims(3)              ! dimension of the 3D grid (nX,nY,nZ)
   real(ReKi)                                 :: origin(3)            ! the lower-left corner of the 3D grid (X0,Y0,Z0)
   real(ReKi)                                 :: gridSpacing(3)       ! spacing between grid points in each of the 3 directions (dX,dY,dZ)
   real(ReKi)                                 :: gridSpacingWAT(3)    ! 
   real(DbKi)                                 :: Time 
   character(1024)                            :: FileName             ! Name of output file     
   character(1024)                            :: desc                 ! Line describing the contents of the file
   character(1024)                            :: vecLabel             ! descriptor of the vector data
   integer(IntKi)                             :: Un                   ! file unit
   integer(IntKi)                             :: n, nt, nh, n_high_low, nhigh
   real(ReKi)                                 :: gridRatio             ! Temporary real for checking WAT resolution
   character(ErrMsgLen)                       :: TmpMsg                ! Temporary Error message text for WAT resolution checks
   real(ReKi), parameter                      :: fstretch = 2.0_ReKi   ! stretching factor for checking WAT resolution
   real(ReKi)                                 :: TargetChunkSize
   integer(IntKi)                             :: nChunksX, nChunksY
   integer(IntKi)                             :: nChunkPointsX, nChunkPointsY
   integer(IntKi), allocatable                :: ChunkIndicesX(:,:), ChunkIndicesY(:,:)
   integer(IntKi)                             :: StartIndexNum, IndexDelta
   
   errStat = ErrID_None
   errMsg  = ""
   
   !============================================================================
   ! Low-resolution grid
   !============================================================================

   ! Get grid origin, spacing, and dimensions depending on wind type
   select case (p%Mod_AmbWind)
   
   ! VTK-based inflow
   case (1)

      ! Parse time 0.0, low res wind input file to gather the grid information
      FileName = trim(p%WindFilePath)//trim(PathSep)//"Low"//trim(PathSep)//"Amb.t0.vtk"
      Un = -1 ! Set to force closing of file on return
      call ReadVTK_SP_info(FileName, desc, dims, origin, gridSpacing, vecLabel, Un, ErrStat2, ErrMsg2)     
      if (Failed()) return     
        
   ! InflowWind-based inflow
   case (2,3)

      ! Grid dimensions from initialization input
      dims(1) = InitInp%InputFileData%nX_low
      dims(2) = InitInp%InputFileData%nY_low
      dims(3) = InitInp%InputFileData%nZ_low

      ! Grid origin from initialization input
      origin(1) = InitInp%InputFileData%X0_low
      origin(2) = InitInp%InputFileData%Y0_low
      origin(3) = InitInp%InputFileData%Z0_low

      ! Grid spacing from initialization input
      gridSpacing(1) = InitInp%InputFileData%dX_low
      gridSpacing(2) = InitInp%InputFileData%dY_low
      gridSpacing(3) = InitInp%InputFileData%dZ_low
   
   ! AMReX-based inflow
   case (4)

      ! Read first low-res file
      FileName = trim(p%WindFilePath)//"_0_"//p%DirStartIndex
      call amrex_read_header(FileName, Time, dims, gridSpacing, origin, ErrStat2, ErrMsg2)
      if (Failed()) return

      ! Search directory for time slices of this sub-volume
      call amrex_find_subvols(p%WindFilePath, 0, p%dt_low, p%NumDT, p%DirStartIndex, &
                              StartIndexNum, p%DirIndexDeltaLow, ErrStat2, ErrMsg2)
      if (Failed()) return

   end select

   !----------------------------------------------------------------------------
   ! Check low-res grid properties
   !----------------------------------------------------------------------------

   ! If the grid spacing is less than or equal to zero, return error
   if (any(gridSpacing <= 0.0_ReKi)) then
      call SetErrStat(ErrID_Fatal, 'Low resolution spatial resolution must be greater than zero in each spatial direction. ', errStat, errMsg, RoutineName)
      return
   endif

   ! If any grid dimensions are less than 2, return error
   if (any(dims < 2)) then
      call SetErrStat(ErrID_Fatal, 'Low resolution grid dimensions must be >= 2 in all directions.', errStat, errMsg, RoutineName)
      return
   end if

   !----------------------------------------------------------------------------
   ! Save properties in low-res grid
   !----------------------------------------------------------------------------

   p%LowRes%oXYZ = origin
   p%LowRes%dXYZ = gridSpacing
   p%LowRes%nXYZ = dims
   p%LowRes%nPoints = product(dims)
   p%LowRes%Size = gridSpacing * real(dims - 1, ReKi)
   p%LowRes%Center = origin + 0.5_ReKi * p%LowRes%Size
   
   ! Polar data
   p%dPol = (gridSpacing(1)+gridSpacing(2)+gridSpacing(3))/3.0_ReKi
   p%n_rp_max = ceiling(pi*((p%C_Meander*((p%NumRadii-1)*InitInp%InputFileData%dr+p%dPol))/p%dPol)**2.0_ReKi)

   ! Set coordinates of points in a flat array
   call AllocAry(p%LowRes%GridPoints, 3, p%LowRes%nPoints, 'LowRes%GridPoints', errStat2, errMsg2); if(Failed()) return
   call flatCartGridCoordinates(origin, dims, gridSpacing, p%LowRes%GridPoints)

   ! Save low-res grid initialization output
   InitOut%oXYZ_Low   = origin
   InitOut%nXYZ_Low   = dims
   InitOut%dXYZ_Low   = gridSpacing

   !----------------------------------------------------------------------------
   ! Check that all low-res wind files have the same grid dimensions
   !----------------------------------------------------------------------------

   if (InitInp%InputFileData%ChkWndFiles) then

      select case (p%Mod_AmbWind)

      ! VTK-based inflow
      case (1)

         ! Loop through step indices
         ! Start at 1 because the first step was already checked
         do n = 1, p%NumDT-1

            ! Read VTK header for low-res wind data
            FileName = trim(p%WindFilePath)//trim(PathSep)//"Low"//trim(PathSep)//"Amb.t"//trim(Num2LStr(n))//".vtk"
            Un = -1 ! Set to force closing of file on return
            call ReadVTK_SP_info(FileName, desc, dims, origin, gridSpacing, vecLabel, Un, ErrStat, ErrMsg) 
               if (ErrStat >= AbortErrLev) return 

            ! Check that VTK properties match grid
            call CheckLowResGridProps(p%LowRes, n, origin, gridSpacing, dims, ErrStat, ErrMsg)
         end do
     
      end select
   end if

   !----------------------------------------------------------------------------
   ! Divide low-resolution grid into chunks for selectively updating the domain
   !----------------------------------------------------------------------------

   ! Approximate a chunk size based on the max wake diameter
   TargetChunkSize = p%y(p%NumRadii-1)

   ! Target number of points in each chunk
   nChunkPointsX = nint(TargetChunkSize/p%LowRes%dXYZ(1), IntKi)
   nChunkPointsY = nint(TargetChunkSize/p%LowRes%dXYZ(2), IntKi)

   nChunksX = int(ceiling(real(p%LowRes%nXYZ(1), ReKi) / real(nChunkPointsX, ReKi)), IntKi)
   nChunksY = int(ceiling(real(p%LowRes%nXYZ(2), ReKi) / real(nChunkPointsY, ReKi)), IntKi)

   ! Allocate arrays to store chunk start and end indices for X dimension, check that each chunk has at least 2 points
   call AllocAry(ChunkIndicesX, 2, nChunksX, 'ChunkIndicesX', errStat2, errMsg2); if(Failed()) return
   call BalancedChunks(p%LowRes%nXYZ(1), nChunksX, ChunkIndicesX, errStat2, errMsg2); if(Failed()) return
   do i = 1, nChunksX
      if (ChunkIndicesX(2,i) - ChunkIndicesX(1,i) + 1 < 2) then
         call SetErrStat(ErrID_Fatal, 'Chunk size in X direction is less than 2.', errStat, errMsg, RoutineName)
         return
      end if
   end do
   
   ! Allocate arrays to store chunk start and end indices for Y dimension, check that each chunk has at least 2 points
   call AllocAry(ChunkIndicesY, 2, nChunksY, 'ChunkIndicesY', errStat2, errMsg2); if(Failed()) return
   call BalancedChunks(p%LowRes%nXYZ(2), nChunksY, ChunkIndicesY, errStat2, errMsg2); if(Failed()) return
   do i = 1, nChunksY
      if (ChunkIndicesY(2,i) - ChunkIndicesY(1,i) + 1 < 2) then
         call SetErrStat(ErrID_Fatal, 'Chunk size in Y direction is less than 2.', errStat, errMsg, RoutineName)
         return
      end if
   end do

   ! Allocate array to store chunk information
   allocate(p%LowRes%WakeChunks(nChunksX*nChunksY), stat=errStat2)
   if (errStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, 'Allocation failure for p%LowRes%WakeChunks array.', errStat, errMsg, RoutineName)
      return
   end if

   ! Loop through chunks and populate data
   k = 0
   do i = 1, nChunksX
      do j = 1, nChunksY
         k = k + 1
         associate (chunk => p%LowRes%WakeChunks(k))
            chunk%iChunk = [i, j, 1]
            chunk%iSubGridX = ChunkIndicesX(:,i)
            chunk%iSubGridY = ChunkIndicesY(:,j)
            chunk%iSubGridZ = [1, p%LowRes%nXYZ(3)] ! All points in Z
            chunk%nPoints = (ChunkIndicesX(2,i) - ChunkIndicesX(1,i) + 1) * &
                            (ChunkIndicesY(2,j) - ChunkIndicesY(1,j) + 1) * &
                            p%LowRes%nXYZ(3)
            call AllocAry(chunk%iGridPoints, chunk%nPoints, "iGridPoints", ErrStat2, ErrMsg2); if (Failed()) return
            call ChunkPointIndices(p%LowRes%nXYZ, chunk%iSubGridX, chunk%iSubGridY, chunk%iSubGridZ, chunk%iGridPoints)
            chunk%oXYZ = [p%LowRes%oXYZ(1) + (chunk%iSubGridX(1) - 1) * p%LowRes%dXYZ(1), &
                          p%LowRes%oXYZ(2) + (chunk%iSubGridY(1) - 1) * p%LowRes%dXYZ(2), &
                          p%LowRes%oXYZ(3)]
            chunk%Size = [(chunk%iSubGridX(2) - chunk%iSubGridX(1)) * p%LowRes%dXYZ(1), &
                          (chunk%iSubGridY(2) - chunk%iSubGridY(1)) * p%LowRes%dXYZ(2), &
                          p%LowRes%Size(3)]
            chunk%Center = chunk%oXYZ + chunk%Size / 2.0_ReKi
            chunk%Radius = norm2(chunk%Size(1:2)) / 2.0_ReKi ! XY only for k-d tree
         end associate
      end do
   end do

   !============================================================================
   ! High-resolution grid
   !============================================================================

   ! Allocate array of high-res parameter data
   allocate(p%HighRes(p%NumTurbines), stat=errStat2); if (errStat2 /= 0) then
      call SetErrStat(ErrID_Fatal, 'Allocation failure for p%HighRes array.', errStat, errMsg, RoutineName)
      return
   end if

   ! Allocate initialization output arrays
   call AllocAry(InitOut%oXYZ_high, 3, p%NumTurbines, 'oXYZ_high', errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitOut%dXYZ_high, 3, p%NumTurbines, 'dXYZ_high', errStat2, errMsg2); if(Failed()) return
   call AllocAry(InitOut%nXYZ_high, 3, p%NumTurbines, 'nXYZ_high', errStat2, errMsg2); if(Failed()) return

   ! If wake-added turbulence is enabled, calculate spacing
   if (p%WAT_Enabled) then
      gridSpacingWAT = [0.0_ReKi, 1/p%WAT_FlowField%Grid3D%InvDY, 1/p%WAT_FlowField%Grid3D%InvDZ]
   endif

   !----------------------------------------------------------------------------
   ! Loop through high-resolution grids (one per turbine)
   !----------------------------------------------------------------------------
   do nt = 1, p%NumTurbines 

      ! Get high-res grid origin, dimensions, and spacing based on wind type
      select case (p%Mod_AmbWind)

      ! VTK-based wind
      case (1)

         Un = -1 ! Set to force closing of file on return
         FileName = trim(p%WindFilePath)//trim(PathSep)//"HighT"//trim(num2lstr(nt))//trim(PathSep)//"Amb.t0.vtk"
         call ReadVTK_SP_info(FileName, desc, dims, origin, gridSpacing, vecLabel, Un, ErrStat2, ErrMsg2) 
         if(Failed()) return
      
      ! InflowWind-based wind
      case (2,3)

         ! Grid dimensions from initialization input (same for all high-res grids)
         dims(1) = InitInp%InputFileData%nX_high
         dims(2) = InitInp%InputFileData%nY_high
         dims(3) = InitInp%InputFileData%nZ_high

         ! Grid origin from initialization input
         origin(1) = InitInp%InputFileData%X0_high(nt)
         origin(2) = InitInp%InputFileData%Y0_high(nt)
         origin(3) = InitInp%InputFileData%Z0_high(nt)

         ! Grid spacing from initialization input
         gridSpacing(1) = InitInp%InputFileData%dX_high(nt)
         gridSpacing(2) = InitInp%InputFileData%dY_high(nt)
         gridSpacing(3) = InitInp%InputFileData%dZ_high(nt)

      ! AMReX-based wind
      case (4)

         FileName = trim(p%WindFilePath)//"_"//trim(Num2LStr(nt))//"_"//p%DirStartIndex
         call amrex_read_header(FileName, Time, dims, gridSpacing, origin, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! Search directory for time slices of this sub-volume
         call amrex_find_subvols(p%WindFilePath, nt, p%dt_high, p%NumDT*p%n_high_low, p%DirStartIndex, &
                                 StartIndexNum, IndexDelta, ErrStat2, ErrMsg2)
         if (Failed()) return

         ! If first turbine, save index delta, otherwise ensure that it is the same
         if (nt == 1) then
            p%DirIndexDeltaHigh = IndexDelta
         else if (p%DirIndexDeltaHigh /= IndexDelta) then
            call SetErrStat(ErrID_Fatal, "got different index delta for sub-volume "//trim(Num2LStr(nt))//" than for sub-volume 1", &
                            ErrStat, ErrMsg, RoutineName)
            return
         end if

      end select

      !-------------------------------------------------------------------------
      ! Check grid properties
      !-------------------------------------------------------------------------

      ! If the grid spacing is less than or equal to zero, return error
      if (any(gridSpacing <= 0.0_ReKi)) then
         call SetErrStat(ErrID_Fatal, 'The high resolution spatial resolution for Turbine '//trim(num2lstr(nt))//' must be greater than zero in each spatial direction. ', errStat, errMsg, RoutineName)
         return
      endif

      ! If any grid dimensions are less than 2, return error
      if (any(dims < 2)) then
         call SetErrStat(ErrID_Fatal, 'High resolution grid dimensions must be >= 2 in all directions.', errStat, errMsg, RoutineName)
         return
      end if
      
      ! Checks for grid spacing
      if (p%WAT_Enabled) then
         do k= 2, 3
            gridRatio = gridSpacing(k)/gridSpacingWAT(k)
            if (gridRatio < 1.0_ReKi/fstretch .or. gridRatio > fstretch)  then
               call SetErrStat(ErrID_Fatal, &
                               'Ratio of high res domain resolution to wake added turbulence resolution should be between '// &
                               trim(Num2LStr(1.0_ReKi/fstretch))//' and '//trim(Num2LStr(fstretch))//', but is '// &
                               trim(Num2LStr(gridSpacing(k)))//' / '//trim(Num2LStr(gridSpacingWAT(k)))// ' = '// &
                               trim(Num2LStr(gridRatio))//' for turbine '//trim(Num2LStr(nt))//' in X.', &
                               errStat, errMsg, RoutineName)
               return
            endif
         end do
      endif

      !-------------------------------------------------------------------------
      ! Set properties in high-res grid
      !-------------------------------------------------------------------------

      p%HighRes(nt)%WT_Position = InitInp%InputFileData%WT_Position(:,nt)
      p%HighRes(nt)%oXYZ = origin
      p%HighRes(nt)%dXYZ = gridSpacing
      p%HighRes(nt)%nXYZ = dims
      p%HighRes(nt)%nPoints = product(dims)
      p%HighRes(nt)%Size = gridSpacing * real(dims - 1, ReKi)
      p%HighRes(nt)%Center = origin + 0.5_ReKi * p%HighRes(nt)%Size
      p%HighRes(nt)%Radius = norm2(p%HighRes(nt)%Size(1:2)) / 2.0_ReKi   ! XY only for k-d tree

      ! Allocate and set coordinates of grid points in a flat array
      call AllocAry(p%HighRes(nt)%GridPoints, 3, p%HighRes(nt)%nPoints, 'HighRes%GridPoints', errStat2, errMsg2); if (Failed()) return

      ! Calculate grid coordinates as a flattened array
      call flatCartGridCoordinates(origin, dims, gridSpacing, p%HighRes(nt)%GridPoints)

      ! Set high-res parameters in InitOut
      InitOut%oXYZ_high(:,nt) = p%HighRes(nt)%oXYZ
      InitOut%dXYZ_high(:,nt) = p%HighRes(nt)%dXYZ
      InitOut%nXYZ_high(:,nt) = p%HighRes(nt)%nXYZ

      !-------------------------------------------------------------------------
      ! Check that all high-res wind files have the same grid properties
      !-------------------------------------------------------------------------
         
      if (InitInp%InputFileData%ChkWndFiles) then

         select case (p%Mod_AmbWind)

         ! VTK-based inflow
         case (1)

            ! We have already checked the first high-res files associated with n=0, 
            ! but need to check the remaining, so for simplicity of code we will repeat the check on the first file
            do n = 0, p%NumDT-1  

               ! We only have one high res input for for the very last low res time step
               if (n == (p%NumDT-1)) then
                  n_high_low = 1
               else
                  n_high_low = p%n_high_low
               end if
               
               ! Loop through high-res time slices
               do nh = 1, n_high_low

                  ! Step index
                  nhigh = nh + n*p%n_high_low - 1

                  ! Read VTK info
                  FileName = trim(p%WindFilePath)//trim(PathSep)//"HighT"//trim(Num2LStr(nt))//trim(PathSep)//"Amb.t"//trim(Num2LStr(nhigh))//".vtk"
                  Un = -1 ! Set to force closing of file on return
                  call ReadVTK_SP_info(FileName, desc, dims, origin, gridSpacing, vecLabel, Un, ErrStat2, ErrMsg2) 
                  if (Failed()) return

                  ! Check that file properties match high-res grid
                  call CheckHighResGridProps(p%HighRes(nt), nt, nhigh, origin, gridSpacing, dims, ErrStat2, ErrMsg2)
                  if (Failed()) return
               end do
            end do

         end select
      end if
   end do   

contains

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
   
end subroutine AWAE_IO_InitGridInfo

pure subroutine ChunkPointIndices(dims, ix, iy, iz, indices)
   integer, intent(in)     :: dims(3), ix(2), iy(2), iz(2)
   integer, intent(inout)  :: indices(:)
   integer :: i, j, k, n
   n = 0
   do k = iz(1), iz(2)
      do j = iy(1), iy(2)
         do i = ix(1), ix(2)
            n = n + 1
            indices(n) = i + (j - 1) * dims(1) + (k - 1) * dims(1) * dims(2)
         end do
      end do
   end do
end subroutine

! Divide 1..N points into n_chunks of nearly equal size
subroutine BalancedChunks(N, n_chunks, chunk_indices, ErrStat, ErrMsg)
   integer, intent(in)           :: N, n_chunks
   integer, intent(out)          :: chunk_indices(:,:)
   integer(IntKi), intent(out)   :: ErrStat
   character(*), intent(out)     :: ErrMsg

   character(*), parameter       :: RoutineName = 'BalancedChunks'
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   integer(IntKi)                :: base, remainder, i, offset

   ErrStat = ErrID_None
   ErrMsg  = ""

   if (n_chunks <= 0) then
      call SetErrStat(ErrID_Fatal, "n_chunks must be positive", ErrStat, ErrMsg, RoutineName)
      return
   end if

   if (N < 0) then
      call SetErrStat(ErrID_Fatal, "N must be nonnegative", ErrStat, ErrMsg, RoutineName)
      return
   end if

   base = N / n_chunks
   remainder = mod(N, n_chunks)

   offset = 0
   do i = 1, n_chunks
      ! First 'remainder' chunks get one extra point
      if (i <= remainder) then
         chunk_indices(1,i) = offset + 1
         chunk_indices(2,i) = offset + base + 1
         offset = offset + base + 1
      else
         chunk_indices(1,i) = offset + 1
         chunk_indices(2,i) = offset + base
         offset = offset + base
      end if
   end do
end subroutine




subroutine CheckLowResGridProps(Grid, step, origin, spacing, dims, ErrStat, ErrMsg)
   type(LRGParamType), intent(in)      :: Grid
   integer(IntKi), intent(in)          :: step
   real(ReKi), intent(in)              :: origin(3), spacing(3)
   integer(IntKi), intent(in)          :: dims(3)
   integer(IntKi), intent(out)         :: ErrStat      ! temporary error status of the operation
   character(ErrMsgLen), intent(out)   :: ErrMsg       ! temporary error message

   character(*), parameter             :: RoutineName = 'CheckLowResGridProps'

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Verify that low-res wind file dimensions match input
   if (any(dims /= Grid%nXYZ)) then
      call SetErrStat(ErrID_Fatal, 'The low resolution grid dimensions for time step '// &
                      trim(Num2LStr(step))//' do not match time step 0.', errStat, errMsg, RoutineName)
      return
   end if 

   ! Verify that low-res wind file origin matches input
   if (any(origin /= Grid%oXYZ)) then
      call SetErrStat(ErrID_Fatal, 'The low resolution grid origins for time step '// &
                      trim(Num2LStr(step))//' do not match time step 0.', errStat, errMsg, RoutineName)
      return
   end if 

   ! Verify that low-res wind file grid spacing matches input
   if (any(spacing /= Grid%dXYZ)) then
      call SetErrStat(ErrID_Fatal, 'The low resolution grid spacing for time step '// &
                      trim(Num2LStr(step))//' do not match time step 0.', errStat, errMsg, RoutineName)
      return
   end if  
end subroutine

subroutine CheckHighResGridProps(Grid, iWT, step, origin, spacing, dims, ErrStat, ErrMsg)
   type(HRGParamType), intent(in)      :: Grid
   integer(IntKi), intent(in)          :: iWT, step
   real(ReKi), intent(in)              :: origin(3), spacing(3)
   integer(IntKi), intent(in)          :: dims(3)
   integer(IntKi), intent(out)         :: ErrStat      ! temporary error status of the operation
   character(ErrMsgLen), intent(out)   :: ErrMsg       ! temporary error message

   character(*), parameter             :: RoutineName = 'CheckHighResGridProps'

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! If inflow grid dimensions do not match grid, return error
   if (any(dims /= Grid%nXYZ)) then
      call SetErrStat(ErrID_Fatal, 'The high resolution grid dimensions for turbine #'//trim(Num2LStr(iWT))// &
                      ' and high-res time step '//trim(Num2LStr(step))//' do not match turbine #1 and time step 0.', &
                      errStat, errMsg, RoutineName)
      return
   end if
   
   ! If inflow grid spacing does not match grid, return error
   if (any(spacing /= Grid%dXYZ)) then
      call SetErrStat(ErrID_Fatal, 'The high resolution grid spacing for turbine #'//trim(Num2LStr(iWT))// &
                      ' and high-res time step '//trim(Num2LStr(step))//' do not match time step 0.', &
                      errStat, errMsg, RoutineName)
      return
   end if
   
   ! If inflow grid origin does not match grid, return error
   if (any(origin /= Grid%oXYZ)) then
      call SetErrStat(ErrID_Fatal, 'The high resolution grid origin for turbine #'//trim(Num2LStr(iWT))// &
                      ' and high-res time step '//trim(Num2LStr(step))//' do not match time step 0.', &
                      errStat, errMsg, RoutineName)
      return
   end if
end subroutine

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

   errStat = ErrID_None
   errMsg  = ""

   ! Open the summary file and give it a heading.
      
   !$OMP critical(fileopen_critical)
   CALL GetNewUnit( UnSu, ErrStat, ErrMsg )
   CALL OpenFOutFile ( UnSu, TRIM( p%OutFileRoot )//'.sum', ErrStat, ErrMsg )
   !$OMP end critical(fileopen_critical)
   IF ( ErrStat >= AbortErrLev ) RETURN

 

   !$OMP critical(fileopen_critical)
   CLOSE(UnSu)
   !$OMP end critical(fileopen_critical)

RETURN
END SUBROUTINE AWAE_PrintSum
!----------------------------------------------------------------------------------------------------------------------------------



END MODULE AWAE_IO

