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
   
   type(ProgDesc), parameter  :: AWAE_Ver = ProgDesc( 'AWAE', '', '' )
   character(*),   parameter  :: AWAE_Nickname = 'AWAE'
    
   public :: AWAE_IO_InitGridInfo
   public :: ReadLowResWindFile
   
   contains

subroutine HiResWindCheck(n, nt, nX, nY, nZ, dX, dY, dZ, X0, Y0, Z0, dims, gridSpacing, origin, callingRoutine, errMsg, errStat)
   integer(IntKi),             intent(in   ) :: n               !< high-resolution time step number (0-based)
   integer(IntKi),             intent(in   ) :: nt              !< turbine number
   integer(IntKi),             intent(in   ) :: nX              !< number of grid points in the X-direction for turbine 1 at high-res time step 0 
   integer(IntKi),             intent(in   ) :: nY              !< number of grid points in the Y-direction for turbine 1 at high-res time step 0 
   integer(IntKi),             intent(in   ) :: nZ              !< number of grid points in the Z-direction for turbine 1 at high-res time step 0 
   real(ReKi),                 intent(in   ) :: dX              !< space between grid points in the X-direction for turbine 1 at high-res time step 0 
   real(ReKi),                 intent(in   ) :: dY              !< space between grid points in the Y-direction for turbine 1 at high-res time step 0 
   real(ReKi),                 intent(in   ) :: dZ              !< space between grid points in the Z-direction for turbine 1 at high-res time step 0 
   real(ReKi),                 intent(in   ) :: X0              !< starting X-location of the grid for turbine 1 at high-res time step 0 (m)
   real(ReKi),                 intent(in   ) :: Y0              !< starting Y-location of the grid for turbine 1 at high-res time step 0 (m)
   real(ReKi),                 intent(in   ) :: Z0              !< starting Z-location of the grid for turbine 1 at high-res time step 0 (m)
   integer(IntKi),             intent(in   ) :: dims(3)         !< dimensions of the grid for turbine nt at high-res time step n (m)
   real(ReKi),                 intent(in   ) :: gridSpacing(3)  !< spacing between grid points for turbine nt at high-res time step n (m)
   real(ReKi),                 intent(in   ) :: origin(3)       !< starting coordinates of the grid for turbine nt at high-res time step n (m)
   character(*),               intent(in   ) :: callingRoutine  !< string containing the name of the calling routine.
   integer(IntKi),             intent(  out) :: errStat         !< Error status of the operation
   character(*),               intent(  out) :: errMsg          !< Error message if errStat /= ErrID_None

      ! grid must have two points in each direction
   if ( (dims(1) < 2) .or. (dims(2) < 2) .or. (dims(3) < 2) ) then
      call SetErrStat ( ErrID_Fatal, 'The high resolution grid dimensions must contain a minimum of 2 nodes in each spatial direction. Turbine #'//trim(num2lstr(nt))//', time step '//trim(num2lstr(n)), errStat, errMsg, callingRoutine )
      return
   end if
   
      ! All turbines and all time steps must have the same grid dimensions due to array allocation assumptions
   if ( ( dims(1) .ne. nX ) .or. ( dims(2) .ne. nY ) .or. ( dims(3) .ne. nZ ) ) then
      call SetErrStat ( ErrID_Fatal, 'The high resolution grid dimensions for turbine #'//trim(num2lstr(nt))//' and high-res time step '//trim(num2lstr(n))//' do not match turbine #1 and time step 0.', errStat, errMsg, callingRoutine )
      return
   end if
   
      ! spacing must be consistent for a given turbine across all time steps
   if ( ( gridSpacing(1) .ne. dX ) .or. ( gridSpacing(2) .ne. dY ) .or. ( gridSpacing(3) .ne. dZ ) ) then
      call SetErrStat ( ErrID_Fatal, 'The high resolution grid spacing for turbine #'//trim(num2lstr(nt))//' and high-res time step '//trim(num2lstr(n))//' do not match time step 0.', errStat, errMsg, callingRoutine )
      return
   end if
   
      ! verify origin of any given turbine is not changing with time step.          
   if ( ( origin(1) .ne. X0 ) .or. ( origin(2) .ne. Y0 ) .or. ( origin(3) .ne. Z0 ) ) then
      call SetErrStat ( ErrID_Fatal, 'The high resolution grid origin for turbine #'//trim(num2lstr(nt))//' and high-res time step '//trim(num2lstr(n))//' do not match time step 0.', errStat, errMsg, callingRoutine )
      return
   end if
   
end subroutine HiResWindCheck

subroutine WriteDisWindFiles( n, WrDisSkp1, p, y, m, errStat, errMsg )
   integer(IntKi),             intent(in   ) :: n            !<  Low-resolution time step increment
   integer(IntKi),             intent(in   ) :: WrDisSkp1    !<  Number of low resolution time step increments per one output increment
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
   INTEGER(IntKi)                          :: nt, n_out 
   REAL(ReKi)                              :: t_out
   character(p%VTK_tWidth)                 :: Tstr  ! string for current VTK write-out step (padded with zeros)
   
   n_out = n/WrDisSkp1
   t_out = n*p%DT_low

   ! TimeStamp
   write(Tstr, '(i' // trim(Num2LStr(p%VTK_tWidth)) //'.'// trim(Num2LStr(p%VTK_tWidth)) // ')') n_out ! TODO use n instead..

   FileName = trim(p%OutFileVTKRoot)//".Low.Dis."//trim(Tstr)//".vtk"
   call WrVTK_SP_header( FileName, "Low resolution disturbed wind for time = "//trim(num2lstr(t_out))//" seconds.", Un, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   call WrVTK_SP_vectors3D( Un, "Velocity", (/p%nX_low,p%nY_low,p%nZ_low/), (/p%X0_low,p%Y0_low,p%Z0_low/), (/p%dX_low,p%dY_low,p%dZ_low/), m%Vdist_low, errStat2, errMsg2 )
      call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
    
   do nt= 1,p%NumTurbines
       ! We are only writing out the first of the high res data for a given low res time step
      
      FileName = trim(p%OutFileVTKRoot)//".HighT"//trim(num2lstr(nt))//".Dis."//trim(Tstr)//".vtk"
      call WrVTK_SP_header( FileName, "High resolution disturbed wind for time = "//trim(num2lstr(t_out))//" seconds.", Un, errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
      call WrVTK_SP_vectors3D( Un, "Velocity", (/p%nX_high,p%nY_high,p%nZ_high/), (/p%X0_high(nt),p%Y0_high(nt),p%Z0_high(nt)/), (/p%dX_high(nt),p%dY_high(nt),p%dZ_high(nt)/), y%Vdist_high(nt)%data(:,:,:,:,0), errStat2, errMsg2 )
         call SetErrStat(ErrStat2, errMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
       
   end do


end subroutine WriteDisWindFiles


!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine ReadLowResWindFile(n, p, Vamb_Low, errStat, errMsg)
   integer(IntKi),                 intent(in   )  :: n            !< Current simulation timestep increment (zero-based)
   type(AWAE_ParameterType),       intent(in   )  :: p            !< Parameters
   real(SiKi),                     intent(inout)  :: Vamb_Low(:,0:,0:,0:)         !< Array which will contain the low resolution grid of ambient wind velocities
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
  
   FileName = trim(p%WindFilePath)//trim(PathSep)//"Low"//trim(PathSep)//"Amb.t"//trim(num2lstr(n))//".vtk"
   Un = 0; ! Initialization different from -1, important to prevent file closing 
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
   real(SiKi),                     intent(inout)  :: Vamb_high(:,0:,0:,0:)         !< Array which will contain the low resolution grid of ambient wind velocities
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
   
   FileName = trim(p%WindFilePath)//trim(PathSep)//"HighT"//trim(num2lstr(nt))//trim(PathSep)//"Amb.t"//trim(num2lstr(n))//".vtk"
   Un = 0; ! Initialization different from -1, important to prevent file closing 
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

   integer(IntKi)                             :: NumGrid_high
   integer(IntKi)                             :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                       :: errMsg2       ! temporary error message 
   character(*), parameter                    :: RoutineName = 'AWAE_IO_InitGridInfo'
   integer(IntKi)                             :: nXYZ_low, nx_low, ny_low, nz_low, nXYZ_high, nx_high, ny_high, nz_high
   integer(IntKi)                             :: dims(3)              ! dimension of the 3D grid (nX,nY,nZ)
   real(ReKi)                                 :: origin(3)            ! the lower-left corner of the 3D grid (X0,Y0,Z0)
   real(ReKi)                                 :: gridSpacing(3)       ! spacing between grid points in each of the 3 directions (dX,dY,dZ)
   character(1024)                            :: FileName             ! Name of output file     
   character(1024)                            :: descr                ! Line describing the contents of the file
   character(1024)                            :: vecLabel             ! descriptor of the vector data
   integer(IntKi)                             :: Un                   ! file unit
   integer(IntKi)                             :: n, nt, nh, n_high_low, nhigh
   
   
   errStat = ErrID_None
   errMsg  = ""
   
   
   !---------------------------------------------------------------------------
   ! Parse time 0.0, low res wind input file to gather the grid 
   !    information and set data associated with the low res grid
   !---------------------------------------------------------------------------

   if ( p%Mod_AmbWind == 1 ) then
      FileName = trim(p%WindFilePath)//trim(PathSep)//"Low"//trim(PathSep)//"Amb.t0.vtk"
      Un = -1 ! Set to force closing of file on return
      call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat, ErrMsg )     
         if (ErrStat >= AbortErrLev) return     
      
      if ( (dims(1) < 2) .or. (dims(2) < 2) .or. (dims(3) < 2) ) then
         call SetErrStat ( ErrID_Fatal, 'The low resolution grid dimensions most contain a minimum of 2 nodes in each spatial direction.', errStat, errMsg, RoutineName )
         return
      end if
   else
      ! Using InflowWind, so data has been passed in via the InitInp data structure
      origin(1) = InitInp%InputFileData%X0_low
      origin(2) = InitInp%InputFileData%Y0_low
      origin(3) = InitInp%InputFileData%Z0_low
      dims(1)   = InitInp%InputFileData%nX_low
      dims(2)   = InitInp%InputFileData%nY_low
      dims(3)   = InitInp%InputFileData%nZ_low
      gridSpacing(1) = InitInp%InputFileData%dX_low
      gridSpacing(2) = InitInp%InputFileData%dY_low
      gridSpacing(3) = InitInp%InputFileData%dZ_low
          
   end if
   
   if ( (gridSpacing(1) <= 0.0_ReKi) .or. (gridSpacing(2) <= 0.0_ReKi) .or. (gridSpacing(3) <= 0.0_ReKi) ) &
      call SetErrStat ( ErrID_Fatal, 'The low resolution spatial resolution for Turbine 1 must be greater than zero in each spatial direction. ', errStat, errMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   p%X0_low           = origin(1)
   p%Y0_low           = origin(2)
   p%Z0_low           = origin(3) 
   p%nX_low           = dims(1)
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
   
   p%NumGrid_low      = p%nX_Low*p%nY_Low*p%nZ_Low
   
   p%dXYZ_Low = gridSpacing
   p%dpol = (gridSpacing(1)+gridSpacing(2)+gridSpacing(3))/3.0_ReKi
   p%n_rp_max = ceiling(pi*((p%C_Meander*(p%r(p%NumRadii-1)+p%dpol))/p%dpol)**2.0_ReKi)
      ! Grid runs from (X0_low, Y0_low, Z0_low) to (X0_low + (p%nX_Low-1)*dX_low, Y0_low+ (p%nY_Low-1)*dY_low, Z0_low+ (p%nZ_Low-1)*dZ_low)
      ! (0,0,0) to (180,180,180) 
     
  
   allocate( p%Grid_low(3,p%NumGrid_low),stat=errStat2)
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
   
   if ( (InitInp%InputFileData%ChkWndFiles) .and. (p%Mod_AmbWind == 1) ) then
      do n=1,p%NumDT-1  ! We have already checked the first low res time step
         
         FileName = trim(p%WindFilePath)//trim(PathSep)//"Low"//trim(PathSep)//"Amb.t"//trim(num2lstr(n))//".vtk"
         Un = -1 ! Set to force closing of file on return
         call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat, ErrMsg ) 
            if (ErrStat >= AbortErrLev) return 
            
            ! verify dims, origin, gridSpacing match the first input file
         if ( ( dims(1) .ne. p%nX_low ) .or. ( dims(2) .ne. p%nY_low ) .or. ( dims(3) .ne. p%nZ_low ) ) then
            call SetErrStat ( ErrID_Fatal, 'The low resolution grid dimensions for time step '//trim(num2lstr(n))//' do not match time step 0.', errStat, errMsg, RoutineName )
            return
         end if 
         if ( ( origin(1) .ne. InitOut%X0_Low ) .or. ( origin(2) .ne. InitOut%Y0_Low ) .or. ( origin(3) .ne. InitOut%Z0_Low ) ) then
            call SetErrStat ( ErrID_Fatal, 'The low resolution grid origins for time step '//trim(num2lstr(n))//' do not match time step 0.', errStat, errMsg, RoutineName )
            return
         end if 
         if ( ( gridSpacing(1) .ne. p%dX_low ) .or. ( gridSpacing(2) .ne. p%dY_low ) .or. ( gridSpacing(3) .ne. p%dZ_low ) ) then
            call SetErrStat ( ErrID_Fatal, 'The low resolution grid spacing for time step '//trim(num2lstr(n))//' do not match time step 0.', errStat, errMsg, RoutineName )
            return
         end if  
         
      end do
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
   
   if ( p%Mod_AmbWind == 1 ) then

      !---------------------------------------------------------------------------
      ! Parse turbine 1, 1st timestep, high res wind input file to gather the grid 
      !    information and set data associated with turbine 1
      !---------------------------------------------------------------------------
   
      FileName = trim(p%WindFilePath)//trim(PathSep)//"HighT1"//trim(PathSep)//"Amb.t0.vtk"
      Un = -1 ! Set to force closing of file on return
      call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, errStat2, errMsg2 ) 
         call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
         if (errStat >= AbortErrLev) return 
   else
      
      ! Using InflowWind, so data has been passed in via the InitInp data structure
      origin(1) = InitInp%InputFileData%X0_high(1)
      origin(2) = InitInp%InputFileData%Y0_high(1)
      origin(3) = InitInp%InputFileData%Z0_high(1)
      dims(1)   = InitInp%InputFileData%nX_high
      dims(2)   = InitInp%InputFileData%nY_high
      dims(3)   = InitInp%InputFileData%nZ_high
      gridSpacing(1) = InitInp%InputFileData%dX_high(1)
      gridSpacing(2) = InitInp%InputFileData%dY_high(1)
      gridSpacing(3) = InitInp%InputFileData%dZ_high(1)
      p%dt_high = InitInp%InputFileData%dt_high
      
   end if
   
   if ( (gridSpacing(1) <= 0.0_ReKi) .or. (gridSpacing(2) <= 0.0_ReKi) .or. (gridSpacing(3) <= 0.0_ReKi) ) &
      call SetErrStat ( ErrID_Fatal, 'The high resolution spatial resolution for Turbine 1 must be greater than zero in each spatial direction. ', errStat, errMsg, RoutineName )
      if (errStat >= AbortErrLev ) return

   p%nX_high          = dims(1)
   p%nY_high          = dims(2)
   p%nZ_high          = dims(3)
   p%X0_high(1)       = origin(1)
   p%Y0_high(1)       = origin(2)
   p%Z0_high(1)       = origin(3)
   p%dX_high(1)       = gridSpacing(1)
   p%dY_high(1)       = gridSpacing(2)
   p%dZ_high(1)       = gridSpacing(3)
   NumGrid_high       = p%nX_high*p%nY_high*p%nZ_high
   
   InitOut%X0_high(1) = origin(1)
   InitOut%Y0_high(1) = origin(2)
   InitOut%Z0_high(1) = origin(3)  
   InitOut%dX_high(1) = gridSpacing(1)
   InitOut%dY_high(1) = gridSpacing(2)
   InitOut%dZ_high(1) = gridSpacing(3)
      
   if ( p%Mod_AmbWind == 1 ) then
         ! Just using this to make sure dims are >=2 points in each direction
      call HiResWindCheck(0, 1, p%nX_high, p%nY_high, p%nZ_high, p%dX_high(1), p%dY_high(1), p%dZ_high(1), p%X0_high(1), p%Y0_high(1), p%Z0_high(1), dims, gridSpacing, origin, RoutineName, errMsg2, errStat2)
         call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
         if (errStat >= AbortErrLev) return 
   end if
   
   allocate( p%Grid_high(3,NumGrid_high,p%NumTurbines ),stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for Grid_high.', errStat, errMsg, RoutineName )
         return
      end if
      
   nXYZ_high = 0
   do nz_high=0, p%nZ_high-1 
      do ny_high=0, p%nY_high-1
         do nx_high=0, p%nX_high-1
            nXYZ_high = nXYZ_high + 1
            p%Grid_high(1,nXYZ_high,1) = InitOut%X0_high(1) + nx_high*InitOut%dX_high(1)
            p%Grid_high(2,nXYZ_high,1) = InitOut%Y0_high(1) + ny_high*InitOut%dY_high(1)
            p%Grid_high(3,nXYZ_high,1) = InitOut%Z0_high(1) + nz_high*InitOut%dZ_high(1)            
         end do
      end do
   end do
 
   !---------------------------------------------------------------------------
   ! Parse the remaining turbine's 1st timestep, high res wind input files to 
   !    gather the grid information and set data associated with those turbines
   !---------------------------------------------------------------------------

   do nt = 2, p%NumTurbines 
      
      if ( p%Mod_AmbWind == 1 ) then
         FileName = trim(p%WindFilePath)//trim(PathSep)//"HighT"//trim(num2lstr(nt))//trim(PathSep)//"Amb.t0.vtk"
         Un = -1 ! Set to force closing of file on return
         call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat2, ErrMsg2 ) 
            call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
            if (ErrStat >= AbortErrLev) return 
      else
         ! Using InflowWind, so data has been passed in via the InitInp data structure
         origin(1) = InitInp%InputFileData%X0_high(nt)
         origin(2) = InitInp%InputFileData%Y0_high(nt)
         origin(3) = InitInp%InputFileData%Z0_high(nt)
         gridSpacing(1) = InitInp%InputFileData%dX_high(nt)
         gridSpacing(2) = InitInp%InputFileData%dY_high(nt)
         gridSpacing(3) = InitInp%InputFileData%dZ_high(nt)

      end if
      
      if ( (gridSpacing(1) <= 0.0_ReKi) .or. (gridSpacing(2) <= 0.0_ReKi) .or. (gridSpacing(3) <= 0.0_ReKi) ) &
         call SetErrStat ( ErrID_Fatal, 'The high resolution spatial resolution for Turbine '//trim(num2lstr(nt))//' must be greater than zero in each spatial direction. ', errStat, errMsg, RoutineName )

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
      
      if ( p%Mod_AmbWind == 1 ) then
            ! Using this to make sure dims are >=2 points in each direction, and number of grid points in each direction matches turbine 1
         call HiResWindCheck(0, nt, p%nX_high, p%nY_high, p%nZ_high, p%dX_high(nt), p%dY_high(nt), p%dZ_high(nt), p%X0_high(nt), p%Y0_high(nt), p%Z0_high(nt), dims, gridSpacing, origin, RoutineName, errMsg2, errStat2)
            call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
            if (ErrStat >= AbortErrLev) return 
      end if
      
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
   
   
   if ( (InitInp%InputFileData%ChkWndFiles) .and. (p%Mod_AmbWind == 1) ) then
      do nt=1,p%NumTurbines
         do n=0,p%NumDT-1  ! We have already checked the first high-res files associated with n=0, but need to check the remaining, so for simplicity of code we will repeat the check on the first file
  
               ! We only have one high res input for for the very last low res timestep
            if ( n == (p%NumDT-1) ) then
               n_high_low = 1
            else
               n_high_low = p%n_high_low
            end if
            
            do nh=1,n_high_low
               nhigh = nh+n*p%n_high_low-1
               FileName = trim(p%WindFilePath)//trim(PathSep)//"HighT"//trim(num2lstr(nt))//trim(PathSep)//"Amb.t"//trim(num2lstr(nhigh))//".vtk"  !TODO: Should the turbine numbers be padding with leading zero(es)? 
               Un = -1 ! Set to force closing of file on return
               call ReadVTK_SP_info( FileName, descr, dims, origin, gridSpacing, vecLabel, Un, ErrStat2, ErrMsg2 ) 
                  call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
                  if (ErrStat >= AbortErrLev) return 
               
               call HiResWindCheck(nhigh, nt, p%nX_high, p%nY_high, p%nZ_high, p%dX_high(nt), p%dY_high(nt), p%dZ_high(nt), p%X0_high(nt), p%Y0_high(nt), p%Z0_high(nt), dims, gridSpacing, origin, RoutineName, errMsg2, errStat2)
                  call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
                  if (errStat >= AbortErrLev ) return
            
            end do
         end do     
      end do      
   end if
   
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

   errStat = ErrID_None
   errMsg  = ""

   ! Open the summary file and give it a heading.
      
   CALL GetNewUnit( UnSu, ErrStat, ErrMsg )
   CALL OpenFOutFile ( UnSu, TRIM( p%OutFileRoot )//'.sum', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN

 

   CLOSE(UnSu)

RETURN
END SUBROUTINE AWAE_PrintSum
!----------------------------------------------------------------------------------------------------------------------------------



END MODULE AWAE_IO

