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
!> VTK output routines for the AWAE module (wake-plane wireframes, structured
!! grid data files, and ParaView .vtk.series index files).
module AWAE_vtk

   use NWTC_Library
   use AWAE_Types
   use VTK

#ifdef _OPENMP
   use OMP_LIB
#endif

   implicit none

   private

   public :: PlaneAxes
   public :: PlaneCorners
   public :: Write_Planes_WireFrame
   public :: Write_Planes_Data
   public :: Write_WakePlane_Data_File
   public :: Write_WakePlane_Series
   public :: Write_WireFrame_Series
   public :: Write_NullPlane

contains

!> Construct the orthonormal in-plane basis (yhat, zhat) in the global
!! inertial frame for a wake plane with unit normal `xhat`. yhat is the
!! horizontal in-plane direction (no Z component); zhat = xhat x yhat.
!! If xhat is purely vertical, yhat falls back to the global Y axis.
subroutine PlaneAxes(xhat, yhat, zhat)
   real(ReKi), intent(in   ) :: xhat(3)
   real(ReKi), intent(  out) :: yhat(3)
   real(ReKi), intent(  out) :: zhat(3)
   real(ReKi)                :: ynorm

   yhat  = (/ -xhat(2), xhat(1), 0.0_ReKi /)
   ynorm = TwoNorm(yhat)
   if (ynorm > 0.0_ReKi) then
      yhat = yhat / ynorm
   else
      yhat = (/ 0.0_ReKi, 1.0_ReKi, 0.0_ReKi /)
   end if

   zhat(1) = xhat(2)*yhat(3) - xhat(3)*yhat(2)
   zhat(2) = xhat(3)*yhat(1) - xhat(1)*yhat(3)
   zhat(3) = xhat(1)*yhat(2) - xhat(2)*yhat(1)
end subroutine PlaneAxes


!> Compute the four corners of a single wake plane in the global inertial
!! frame from the plane normal `xhat`, the plane center `pc`, and the
!! in-plane half-extents `hy` (horizontal) and `hz` (vertical-ish).
!! Corners are returned counter-clockwise about +xhat.
subroutine PlaneCorners(xhat, pc, hy, hz, corners)
   real(ReKi), intent(in   ) :: xhat(3)         !< Plane normal (unit vector)
   real(ReKi), intent(in   ) :: pc(3)           !< Plane center, global frame
   real(ReKi), intent(in   ) :: hy              !< In-plane horizontal half-extent
   real(ReKi), intent(in   ) :: hz              !< In-plane vertical-ish half-extent
   real(ReKi), intent(  out) :: corners(3,4)    !< Four corner positions, global frame

   real(ReKi) :: yhat(3), zhat(3)

   call PlaneAxes(xhat, yhat, zhat)

   ! Four corners, ordered counter-clockwise about +xhat
   corners(:,1) = pc - hy*yhat - hz*zhat
   corners(:,2) = pc + hy*yhat - hz*zhat
   corners(:,3) = pc + hy*yhat + hz*zhat
   corners(:,4) = pc - hy*yhat + hz*zhat
end subroutine PlaneCorners


!> Write the four edges of every active wake plane for each turbine as a
!! wireframe to a per-turbine VTK polydata file.  One file is produced per
!! turbine per output step, named .T<nt>.WakePlanesWireFrame.<t>.vtk.
subroutine Write_Planes_WireFrame(p, u, t, Tstr)
   type(AWAE_ParameterType), intent(in) :: p      !< AWAE parameters
   type(AWAE_InputType),     intent(in) :: u      !< AWAE inputs
   real(DbKi),               intent(in) :: t      !< Current simulation time in seconds
   character(*),             intent(in) :: Tstr   !< Zero-padded time-step string for filenames

   integer(IntKi)              :: nt_wp, np_wp, nActive, ipt
   real(ReKi)                  :: pc(3), corners(3,4)
   real(ReKi)                  :: hy, hz
   real(ReKi),     allocatable :: WPPoints(:,:)
   integer(IntKi), allocatable :: WPLines(:,:)
   type(VTK_Misc)              :: mvtk
   character(1024)             :: WPFileName
   integer(IntKi)              :: ntWidth
   character(16)               :: TurbNum, FmtStrT

   ! Plane half-extents in the local Y and Z (in-plane) directions
   hy = p%y(p%NumRadii-1)
   hz = p%z(p%NumRadii-1)

   ! Zero-padded width sufficient to represent all turbine indices
   ntWidth = max(1, int(floor(log10(real(max(p%NumTurbines, 1), ReKi))) + 1, IntKi))
   write(FmtStrT, '(A,I0,A,I0,A)') '(I', ntWidth, '.', ntWidth, ')'

   do nt_wp = 1, p%NumTurbines
      nActive = NINT(u%NumPlanes(nt_wp))
      if (nActive <= 0) cycle

      allocate(WPPoints(3, 4*nActive))
      allocate(WPLines(2, 4*nActive))

      do np_wp = 0, nActive - 1
         pc = u%p_plane(:, np_wp, nt_wp)
         call PlaneCorners(u%xhat_plane(:, np_wp, nt_wp), pc, hy, hz, corners)

         ipt = 4*np_wp
         WPPoints(:, ipt+1) = corners(:,1)
         WPPoints(:, ipt+2) = corners(:,2)
         WPPoints(:, ipt+3) = corners(:,3)
         WPPoints(:, ipt+4) = corners(:,4)

         ! Four edges of the closed quad outline (0-based point indices for VTK)
         WPLines(:, ipt+1) = (/ ipt,   ipt+1 /)
         WPLines(:, ipt+2) = (/ ipt+1, ipt+2 /)
         WPLines(:, ipt+3) = (/ ipt+2, ipt+3 /)
         WPLines(:, ipt+4) = (/ ipt+3, ipt   /)
      end do

      write(TurbNum, FmtStrT) nt_wp
      WPFileName = trim(p%OutFileFFvtkWakeRoot)//".T"//trim(TurbNum)// &
                   ".WakePlanesWireFrame."//trim(Tstr)//".vtk"

      call vtk_misc_init(mvtk)
      if (vtk_new_ascii_file(WPFileName, &
          "Wake plane wireframes for turbine "//trim(TurbNum)// &
          " at time = "//trim(num2lstr(t))//" seconds.", mvtk)) then
         call vtk_dataset_polydata(WPPoints, mvtk, .false.)
         call vtk_lines(WPLines, mvtk)
         call vtk_close_file(mvtk)
      end if

      deallocate(WPPoints, WPLines)

   end do
end subroutine Write_Planes_WireFrame


!> Write one VTK STRUCTURED_GRID file per active wake plane containing the
!! wake velocity sampled on the plane's structured Y-Z grid. Point
!! coordinates and velocity vectors are written in the global inertial frame.
subroutine Write_Planes_Data(p, u, m, n, t, Tstr)
   type(AWAE_ParameterType), intent(in)    :: p      !< AWAE parameters
   type(AWAE_InputType),     intent(in)    :: u      !< AWAE inputs
   type(AWAE_MiscVarType),   intent(inout) :: m      !< AWAE misc variables
   integer(IntKi),           intent(in)    :: n      !< Current low-resolution time step index
   real(DbKi),               intent(in)    :: t      !< Current simulation time in seconds
   character(*),             intent(in)    :: Tstr   !< Zero-padded time-step string for filenames

   integer(IntKi)              :: nt_wp, np_wp, nActive, ntWidth
   integer(IntKi)              :: jy, kz, idx, nY, nZ
   real(ReKi)                  :: pc(3), xhat(3), yhat(3), zhat(3)
   real(ReKi),     allocatable :: Pts(:,:), Vel(:,:)
   character(1024)             :: WPFileName
   character(16)               :: FmtStrWk, FmtStrT
   character(p%VTK_tWidthPlanes) :: PlaneNum
   character(6)                :: TurbNumStr

   ! Plane grid dimensions (p%y and p%z run from -NumRadii+1 to NumRadii-1)
   nY = 2*p%NumRadii - 1
   nZ = 2*p%NumRadii - 1

   ! Field width for plane index, based on MaxPlanes (1-based count)
   write(FmtStrWk, '(A,I0,A,I0,A)') '(I', p%VTK_tWidthPlanes, '.', p%VTK_tWidthPlanes, ')'
   ! Zero-padded width sufficient to represent all turbine indices
   ntWidth = max(1, int(floor(log10(real(max(p%NumTurbines, 1), ReKi))) + 1, IntKi))
   write(FmtStrT, '(A,I0,A,I0,A)') '(A1,I', ntWidth, '.', ntWidth, ')'

   allocate(Pts(3, nY*nZ))
   allocate(Vel(3, nY*nZ))

   do nt_wp = 1, p%NumTurbines
      ! Turbine number
      write(TurbNumStr, FmtStrT) "T", nt_wp

      ! Number of planes active
      nActive = NINT(u%NumPlanes(nt_wp))

      do np_wp = 0, nActive - 1
         pc   = u%p_plane(:,    np_wp, nt_wp)
         xhat = u%xhat_plane(:, np_wp, nt_wp)
         call PlaneAxes(xhat, yhat, zhat)

         ! Fill points and velocity in VTK natural order (jy fastest, then kz)
         do kz = 1, nZ
            do jy = 1, nY
               idx = jy + (kz-1)*nY
               Pts(:, idx) = pc + p%y(jy-p%NumRadii) * yhat &
                                + p%z(kz-p%NumRadii) * zhat
               Vel(:, idx) = u%Vx_wake(jy-p%NumRadii, kz-p%NumRadii, np_wp, nt_wp) * xhat &
                           + u%Vy_wake(jy-p%NumRadii, kz-p%NumRadii, np_wp, nt_wp) * yhat &
                           + u%Vz_wake(jy-p%NumRadii, kz-p%NumRadii, np_wp, nt_wp) * zhat
            end do
         end do

         ! Per-plane index string with consistent zero-padded width
         write(PlaneNum, FmtStrWk) np_wp

         WPFileName = trim(p%OutFileFFvtkWakeRoot)//"."//trim(TurbNumStr)//".WakePlane_"//trim(PlaneNum)// &
                      "."//trim(Tstr)//".vtk"

         call Write_WakePlane_Data_File(WPFileName, &
              "Wake plane "//trim(PlaneNum)//" at time = "// &
              trim(num2lstr(t))//" seconds.", nY, nZ, Pts, Vel)

         ! track what plane this was first written out at (initialized to huge, so this will catch only the first)
         if (m%WakeVTK_StartN(np_wp,nt_wp) > n)    m%WakeVTK_StartN(np_wp,nt_wp) = n
      end do
   end do

   deallocate(Pts, Vel)
end subroutine Write_Planes_Data


!> Helper: write a single 2D VTK STRUCTURED_GRID file for one wake plane
!! containing point coordinates and a "WakeVelocity" point-data vector.
subroutine Write_WakePlane_Data_File(WPFileName, label, n1, n2, Pts, Vel)
   character(*),   intent(in) :: WPFileName
   character(*),   intent(in) :: label
   integer(IntKi), intent(in) :: n1, n2
   real(ReKi),     intent(in) :: Pts(:,:)    !< 3 x (n1*n2)
   real(ReKi),     intent(in) :: Vel(:,:)    !< 3 x (n1*n2)
   type(VTK_Misc)             :: mvtk

   call vtk_misc_init(mvtk)
   if (vtk_new_ascii_file(WPFileName, label, mvtk)) then
      call vtk_dataset_structured_grid(Pts, n1, n2, 1, mvtk)
      call vtk_point_data_init(mvtk)
      call vtk_point_data_vector(Vel, "WakeVelocityDeficit", mvtk)
      call vtk_close_file(mvtk)
   end if
end subroutine Write_WakePlane_Data_File


!> Write the final ParaView .vtk.series index files for all wake planes
!! that were output during the simulation.  Called once from AWAE_End.
!! Each series file lists every VTK output step, referencing either the
!! actual per-plane VTK file or a null-data placeholder for steps before
!! the plane first appeared.
subroutine Write_WakePlane_Series(p, m)
   type(AWAE_ParameterType), intent(in) :: p   !< AWAE parameters
   type(AWAE_MiscVarType),   intent(in) :: m   !< AWAE misc variables

   integer(IntKi)              :: nt_wp, np_wp, ntWidth, n_final, n_out
   integer(IntKi)              :: UnSer, out_idx, SerErrStat
   character(ErrMsgLen)        :: SerErrMsg
   character(16)               :: FmtStrWk, FmtStrT
   character(p%VTK_tWidthPlanes) :: PlaneNum
   character(6)                :: TurbNumStr
   character(1024)             :: SeriesFile, baseName, EntryName, VTKprefix
   character(p%VTK_tWidth)     :: TstrOut
   character(32)               :: TimeStr
   real(DbKi)                  :: t_out
   logical                     :: firstEntry

   if (.not. p%WrPlanes) return

   ! Final low-res time step and total number of VTK output steps
   n_final = p%NumDT - 1
   n_out   = n_final / p%WrDisSkp1

   ! Compute basename once (shared across all series files)
   call GetPath(p%OutFileFFvtkWakeRoot, EntryName, baseName)

   ! Format strings for zero-padded plane and turbine numbers
   write(FmtStrWk, '(A,I0,A,I0,A)') '(I', p%VTK_tWidthPlanes, '.', p%VTK_tWidthPlanes, ')'
   ntWidth = max(1, int(floor(log10(real(max(p%NumTurbines, 1), ReKi))) + 1, IntKi))
   write(FmtStrT, '(A,I0,A,I0,A)') '(A1,I', ntWidth, '.', ntWidth, ')'

   do nt_wp = 1, p%NumTurbines
      write(TurbNumStr, FmtStrT) "T", nt_wp
      do np_wp = 0, p%MaxPlanes - 1
         ! Skip planes that were never written during the simulation
         if (m%WakeVTK_StartN(np_wp, nt_wp) > n_final) cycle

         write(PlaneNum, FmtStrWk) np_wp
         VTKprefix  = trim(TurbNumStr)//".WakePlane_"//trim(PlaneNum)
         SeriesFile = trim(p%OutFileFFvtkWakeRoot)//"."//trim(VTKprefix)//".vtk.series"

         !$OMP critical(fileopen_critical)
         call GetNewUnit(UnSer, SerErrStat, SerErrMsg)
         call OpenFOutFile(UnSer, SeriesFile, SerErrStat, SerErrMsg)
         !$OMP end critical(fileopen_critical)
         if (SerErrStat >= AbortErrLev) cycle

         write(UnSer, '(A)') '{'
         write(UnSer, '(A)') '  "file-series-version" : "1.0",'
         write(UnSer, '(A)') '  "files" : ['

         firstEntry = .true.
         do out_idx = 0, n_out
            t_out = real(out_idx, DbKi) * real(p%WrDisSkp1, DbKi) * p%dt_low
            write(TstrOut, '(i'//trim(Num2LStr(p%VTK_tWidth))//'.'// &
                  trim(Num2LStr(p%VTK_tWidth))//')') out_idx

            if (m%WakeVTK_StartN(np_wp, nt_wp) > out_idx) then
               EntryName = trim(p%OutFileFFvtkWakeNullData)
            else
               EntryName = trim(baseName)//"."//trim(VTKprefix)//"."//trim(TstrOut)//".vtk"
            endif

            write(TimeStr, '(F14.5)') t_out
            if (firstEntry) then
               write(UnSer, '(A,A,A,A,A)') '    { "name" : "', trim(EntryName), &
                                           '", "time" : ', trim(TimeStr), ' }'
               firstEntry = .false.
            else
               write(UnSer, '(A,A,A,A,A)') '   ,{ "name" : "', trim(EntryName), &
                                           '", "time" : ', trim(TimeStr), ' }'
            end if
         end do

         write(UnSer, '(A)') '  ]'
         write(UnSer, '(A)') '}'

         !$OMP critical(fileopen_critical)
         close(UnSer)
         !$OMP end critical(fileopen_critical)
      end do
   end do
end subroutine Write_WakePlane_Series

!> Write a zeroed-out "null" wake-plane VTK file used as a placeholder
!! in the ParaView series for time steps before a plane first appears.
subroutine Write_NullPlane(OutFileVTKwakeDir, p)
   character(*),             intent(in) :: OutFileVTKwakeDir  !< Directory for VTK wake output files
   type(AWAE_ParameterType), intent(in) :: p                  !< AWAE parameters

   real(ReKi), allocatable :: PtsVel(:,:)
   integer(IntKi)          :: nY, nZ

   nY = 2*p%NumRadii - 1
   nZ = 2*p%NumRadii - 1
   allocate(PtsVel(3, nY*nZ))
   PtsVel = 0.0_ReKi
   call Write_WakePlane_Data_File(trim(OutFileVTKwakeDir)//PathSep//p%OutFileFFvtkWakeNullData, &
        "Wake plane NULL at time = NULL seconds.", nY, nZ, PtsVel, PtsVel)
   deallocate(PtsVel)
end subroutine Write_NullPlane


!> Write one ParaView .vtk.series index file per turbine for the wireframe
!! VTK files produced by Write_Planes_WireFrame.  Called once from AWAE_End.
subroutine Write_WireFrame_Series(p)
   type(AWAE_ParameterType), intent(in) :: p   !< AWAE parameters

   integer(IntKi)              :: nt_wp, ntWidth, n_out
   integer(IntKi)              :: UnSer, out_idx, SerErrStat
   character(ErrMsgLen)        :: SerErrMsg
   character(16)               :: FmtStrT, TurbNum
   character(1024)             :: SeriesFile, baseName, EntryName
   character(p%VTK_tWidth)     :: TstrOut
   character(32)               :: TimeStr
   real(DbKi)                  :: t_out
   logical                     :: firstEntry

   if (.not. p%WrPlanes) return

   ! Total number of VTK output steps
   n_out = (p%NumDT - 1) / p%WrDisSkp1

   ! Compute basename once (filename portion of OutFileFFvtkWakeRoot)
   call GetPath(p%OutFileFFvtkWakeRoot, EntryName, baseName)

   ! Format string for zero-padded turbine number
   ntWidth = max(1, int(floor(log10(real(max(p%NumTurbines, 1), ReKi))) + 1, IntKi))
   write(FmtStrT, '(A,I0,A,I0,A)') '(I', ntWidth, '.', ntWidth, ')'

   do nt_wp = 1, p%NumTurbines
      write(TurbNum, FmtStrT) nt_wp

      SeriesFile = trim(p%OutFileFFvtkWakeRoot)//".T"//trim(TurbNum)// &
                   ".WakePlanesWireFrame.vtk.series"

      !$OMP critical(fileopen_critical)
      call GetNewUnit(UnSer, SerErrStat, SerErrMsg)
      call OpenFOutFile(UnSer, SeriesFile, SerErrStat, SerErrMsg)
      !$OMP end critical(fileopen_critical)
      if (SerErrStat >= AbortErrLev) cycle

      write(UnSer, '(A)') '{'
      write(UnSer, '(A)') '  "file-series-version" : "1.0",'
      write(UnSer, '(A)') '  "files" : ['

      firstEntry = .true.
      do out_idx = 0, n_out
         t_out = real(out_idx, DbKi) * real(p%WrDisSkp1, DbKi) * p%dt_low
         write(TstrOut, '(i'//trim(Num2LStr(p%VTK_tWidth))//'.'// &
               trim(Num2LStr(p%VTK_tWidth))//')') out_idx

         EntryName = trim(baseName)//".T"//trim(TurbNum)// &
                     ".WakePlanesWireFrame."//trim(TstrOut)//".vtk"

         write(TimeStr, '(F14.5)') t_out
         if (firstEntry) then
            write(UnSer, '(A,A,A,A,A)') '    { "name" : "', trim(EntryName), &
                                        '", "time" : ', trim(TimeStr), ' }'
            firstEntry = .false.
         else
            write(UnSer, '(A,A,A,A,A)') '   ,{ "name" : "', trim(EntryName), &
                                        '", "time" : ', trim(TimeStr), ' }'
         end if
      end do

      write(UnSer, '(A)') '  ]'
      write(UnSer, '(A)') '}'

      !$OMP critical(fileopen_critical)
      close(UnSer)
      !$OMP end critical(fileopen_critical)
   end do
end subroutine Write_WireFrame_Series

end module AWAE_vtk
