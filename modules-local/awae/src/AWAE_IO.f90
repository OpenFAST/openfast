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
   call ScanDir(trim(p%WindFileRoot)//'\low-res','dirContents.txt')
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
   
!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine ReadLowResWindFile(t, p, Vamb_Low, errStat, errMsg)
   real(DbKi),                     intent(in   )  :: t            !< Current simulation timestep (s)
   type(AWAE_ParameterType),       intent(in   )  :: p            !< Parameters
   real(ReKi),                     intent(inout)  :: Vamb_Low(:,0:,0:,0:)         !< Array which will contain the low resolution grid of ambient wind velocities
   integer(IntKi),                 intent(  out)  :: errStat      !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg       !< Error message if errStat /= ErrID_None
  
   errStat = ErrID_None
   errMsg  = ""

!==============================================================================
! Place holder until we have an actual file format for the ambient wind data
! TODO: Replace with actual file I/O code
   
   Vamb_Low(1,:,:,:) = 8.0_ReKi  ! m/s
   Vamb_Low(2,:,:,:) = 0.0_ReKi  ! m/s
   Vamb_Low(3,:,:,:) = 0.0_ReKi  ! m/s
   
!==============================================================================

   
end subroutine ReadLowResWindFile

!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine ReadHighResWindFile(nt, n_hl, t, p, Vamb_high, errStat, errMsg)

   integer(IntKi),                 intent(in   )  :: nt
   integer(IntKi),                 intent(in   )  :: n_hl
   real(DbKi),                     intent(in   )  :: t            !< Current simulation timestep (s)
   type(AWAE_ParameterType),       intent(in   )  :: p            !< Parameters
   real(ReKi),                     intent(inout)  :: Vamb_high(:,0:,0:,0:)         !< Array which will contain the low resolution grid of ambient wind velocities
   integer(IntKi),                 intent(  out)  :: errStat      !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg       !< Error message if errStat /= ErrID_None
  
   errStat = ErrID_None
   errMsg  = ""

!==============================================================================
! Place holder until we have an actual file format for the ambient wind data
! TODO: Replace with actual file I/O code
   
   !Vamb_high(:,:,:,:,n_hl,nt) = 10.0_ReKi  ! m/s
   Vamb_high(1,:,:,:) = 8.0_ReKi  ! m/s
   Vamb_high(2,:,:,:) = 0.0_ReKi  ! m/s
   Vamb_high(3,:,:,:) = 0.0_ReKi  ! m/s
   
!==============================================================================
  
end subroutine ReadHighResWindFile

subroutine AWAE_IO_InitGridInfo(p, InitOut, errStat, errMsg)

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
   errStat = ErrID_None
   errMsg  = ""

!==============================================================================
! Start simulated read of low and high res ambient wind files 
! TODO: Replace this block with code which actually parses ambient wind data files
   
   call ReadLowResWindFileHeaders(p, errStat, errMsg)
   
   X0_low = -75.0_ReKi
   Y0_low = -500.0_ReKi
   Z0_low = 0.0_ReKi
   dX_low = 10.0_ReKi
   dY_low = 10.0_ReKi
   dZ_low = 10.0_ReKi
   nXYZ_low = 0
   dt_low = 1.0
   
   p%n_wind_min = 100
   p%n_wind_max = ceiling(30.0_ReKi*pi*(2.0_ReKi*p%r(p%NumRadii-1))**2*dt_low/(dX_low*dY_low*dZ_low))
   
      ! Parse a low res wind input file to gather the grid information
   p%nX_Low           = 151      ! 10 ! 
   p%nY_low           = 101      ! 10 ! 
   p%nZ_low           = 51      ! 10 ! 
   NumGrid_low        = p%nX_Low*p%nY_Low*p%nZ_Low
      ! Grid runs from (X0_low, Y0_low, Z0_low) to (X0_low + (p%nX_Low-1)*dX_low, Y0_low+ (p%nY_Low-1)*dY_low, Z0_low+ (p%nZ_Low-1)*dZ_low)
      ! (0,0,0) to (180,180,180) 
      ! Parse a high res wind input file to gather the grid information
   
   ! TODO: Error checking to see that all turbines use the same nX, nY, nZ for the high res grids
   p%nX_high          = 16 ! 38      ! 10 ! 
   p%nY_high          = 16 ! 38      ! 10 ! 
   p%nZ_high          = 16 ! 34      ! 10 ! 
   NumGrid_high       = p%nX_high*p%nY_high*p%nZ_high
   
      ! Determine the number of high res timesteps per a single low res time step by parsing folder names?
   !p%n_high_low       = 5
     
   allocate( p%Grid_low(3,NumGrid_low),stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for Grid_low.', errStat, errMsg, RoutineName )
         return
      end if
   do nz_low=0, p%nZ_low-1 
      do ny_low=0, p%nY_low-1
         do nx_low=0, p%nX_low-1
            nXYZ_low = nXYZ_low + 1
            p%Grid_low(1,nXYZ_low) = X0_low + nx_low*dX_low
            p%Grid_low(2,nXYZ_low) = Y0_low + ny_low*dY_low
            p%Grid_low(3,nXYZ_low) = Z0_low + nz_low*dZ_low            
         end do
      end do
   end do
   
   allocate( p%Grid_high(3,NumGrid_high,p%NumTurbines ),stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for Grid_high.', errStat, errMsg, RoutineName )
         return
      end if
      
   allocate( InitOut%X0_high(p%NumTurbines), InitOut%Y0_high(p%NumTurbines), InitOut%Z0_high(p%NumTurbines), stat=errStat2)   
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for InitOut origin arrays.', errStat, errMsg, RoutineName )
   allocate( InitOut%dX_high(p%NumTurbines), InitOut%dY_high(p%NumTurbines), InitOut%dZ_high(p%NumTurbines), stat=errStat2)   
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for InitOut spatial increment arrays.', errStat, errMsg, RoutineName )

   if (ErrStat >= AbortErrLev) return
         
   do nt = 1, p%NumTurbines 
      InitOut%X0_high(nt) = -75.0_ReKi
      InitOut%Y0_high(nt) = -75.0_ReKi
      InitOut%Z0_high(nt) = 0.0_ReKi
      
      InitOut%dX_high(nt) = 10.0_ReKi
      InitOut%dY_high(nt) = 10.0_ReKi
      InitOut%dZ_high(nt) = 10.0_ReKi
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
   
   !TODO:  Set the corresponding InitOut data 
   !TODO:  Perform any error checking on InitOut here
   
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

