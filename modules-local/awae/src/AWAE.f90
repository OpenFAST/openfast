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
!>  AWAE is a time-domain module for modeling Ambient Wind and Array Effects of one or more horizontal-axis wind turbines.
module AWAE
    
   use NWTC_Library
   use AWAE_Types
   use AWAE_IO
#ifdef PARALLEL_CODE
   use OMP_LIB
#endif
   
   implicit none

   private
         

   ! ..... Public Subroutines ...................................................................................................

   public :: AWAE_Init                           ! Initialization routine
   public :: AWAE_End                            ! Ending routine (includes clean up)
   public :: AWAE_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   public :: AWAE_CalcOutput                     ! Routine for computing outputs
   public :: AWAE_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   
   
      ! Unit testing routines
   public :: AWAE_TEST_Init_BadData    
   public :: AWAE_TEST_Init_GoodData    
   public :: AWAE_TEST_CalcOutput
#ifdef PARALLEL_CODE
   public :: AWAE_TEST_ExtractSlice
   public :: AWAE_TEST_LowResGridCalcs
#endif

   contains  
   
#ifdef PARALLEL_CODE
subroutine ExtractSliceOMP( sliceType, s, s0, szs, sz1, sz2, ds,  V, slice)


   integer(IntKi),      intent(in   ) :: sliceType  !< Type of slice: XYSlice, YZSlice, XZSlice
   real(ReKi),          intent(in   ) :: s          !< data value in meters of the interpolatan
   real(ReKi),          intent(in   ) :: s0         !< origin value in meters of the interpolatan
   integer(IntKi),      intent(in   ) :: szs
   integer(IntKi),      intent(in   ) :: sz1        !< 1st dimension of slice
   integer(IntKi),      intent(in   ) :: sz2        !< 2nd dimension of slice
   real(ReKi),          intent(in   ) :: ds
   real(SiKi),          intent(in   ) :: V(:,0:,0:,0:)
   real(SiKi),          intent(inout) :: slice(:,0:,0:)
   
   integer(IntKi)   :: s_grid0,s_grid1,i,j,m, procs
   real(SiKi)       :: s_grid, sd
   
      
      ! s is in meters but all the s_grid variables are in the grid units so that we can index into the grid arrays properly
      ! NOTE: The grid coordinates run from 0 to sz-1 
      
   s_grid  = real((s-s0)/ds,SiKi)
   
      ! Lower bounds of grid cell in interpolation direction
   s_grid0 = floor(s_grid)
   
      ! Upper bounds of grid cell in interpolation direction
   s_grid1 = s_grid0 + 1   
   
      ! fractional distance of requested slice from lower cell bounds in the range [0-1]
   sd = (s_grid-real(s_grid0,SiKi)) 
   
   if (s_grid0 == (szs-1)) s_grid1 = s_grid0  ! Handle case where s0 is the last index in the grid, in this case sd = 0.0, so the 2nd term in the interpolation will not contribute
  
   procs = omp_get_num_procs()
   call omp_set_num_threads(procs)
   
   !$OMP PARALLEL DO PRIVATE(i,j) SHARED(slice,sz1,sz2,sliceType,V,sd,s_grid0,s_grid1) DEFAULT(NONE)   
   do m = 0,(sz2)*(sz1)-1
      j = m/(sz1)
      i = mod(m,sz1)
      select case (sliceType)
      case (XYSlice)
         slice(:,i,j) = V(:,i,j,s_grid0)*(1.0_SiKi-sd) + V(:,i,j,s_grid1)*sd
      case (YZSlice)
         slice(:,i,j) = V(:,s_grid0,i,j)*(1.0_SiKi-sd) + V(:,s_grid1,i,j)*sd
      case (XZSlice)
         slice(:,i,j) = V(:,i,s_grid0,j)*(1.0_SiKi-sd) + V(:,i,s_grid1,j)*sd
      end select
   end do
   !$OMP END PARALLEL DO
   
end subroutine ExtractSliceOMP
#endif

subroutine ExtractSlice( sliceType, s, s0, szs, sz1, sz2, ds,  V, slice)
   
   integer(IntKi),      intent(in   ) :: sliceType  !< Type of slice: XYSlice, YZSlice, XZSlice
   real(ReKi),          intent(in   ) :: s          !< data value in meters of the interpolatan
   real(ReKi),          intent(in   ) :: s0         !< origin value in meters of the interpolatan
   integer(IntKi),      intent(in   ) :: szs
   integer(IntKi),      intent(in   ) :: sz1        !< 1st dimension of slice
   integer(IntKi),      intent(in   ) :: sz2        !< 2nd dimension of slice
   real(ReKi),          intent(in   ) :: ds
   real(SiKi),          intent(in   ) :: V(:,0:,0:,0:)
   real(SiKi),          intent(inout) :: slice(:,0:,0:)
   
   integer(IntKi)   :: s_grid0,s_grid1,i,j
   real(SiKi)       :: s_grid, sd
   
      
      ! s is in meters but all the s_grid variables are in the grid units so that we can index into the grid arrays properly
      ! NOTE: The grid coordinates run from 0 to sz-1 
      
   s_grid  = real((s-s0)/ds,SiKi)
   
      ! Lower bounds of grid cell in interpolation direction
   s_grid0 = floor(s_grid)
   
      ! Upper bounds of grid cell in interpolation direction
   s_grid1 = s_grid0 + 1   
   
      ! fractional distance of requested slice from lower cell bounds in the range [0-1]
   sd = (s_grid-real(s_grid0,SiKi)) 
   
   if (s_grid0 == (szs-1)) s_grid1 = s_grid0  ! Handle case where s0 is the last index in the grid, in this case sd = 0.0, so the 2nd term in the interpolation will not contribute
  
   do j = 0,sz2-1
      do i = 0,sz1-1
         select case (sliceType)
         case (XYSlice)
            slice(:,i,j) = V(:,i,j,s_grid0)*(1.0_SiKi-sd) + V(:,i,j,s_grid1)*sd
         case (YZSlice)
            slice(:,i,j) = V(:,s_grid0,i,j)*(1.0_SiKi-sd) + V(:,s_grid1,i,j)*sd
         case (XZSlice)
            slice(:,i,j) = V(:,i,s_grid0,j)*(1.0_SiKi-sd) + V(:,i,s_grid1,j)*sd
         end select
      end do
   end do
   
end subroutine ExtractSlice
!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine ComputeLocals(n, u, p, y, m, errStat, errMsg)
   integer(IntKi),                 intent(in   )  :: n           !< Current simulation time increment (zero-based)
   type(AWAE_InputType),           intent(in   )  :: u           !< Inputs at Time t
   type(AWAE_ParameterType),       intent(in   )  :: p           !< Parameters
   type(AWAE_OutputType),          intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   type(AWAE_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                 intent(  out)  :: errStat     !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None

   integer(IntKi)      :: nt, np, maxPln
   real(ReKi)          :: cosTerm, sinTerm, dp(3), rmax
   
   errStat = 0
   errMsg  = ""
   maxPln =   min(n,p%NumPlanes-2)
   rmax = p%r(p%NumRadii-1)
   do nt = 1,p%NumTurbines
      do np = 0, maxPln
         cosTerm = dot_product(u%xhat_plane(:,np+1,nt),u%xhat_plane(:,np,nt))
         if (EqualRealNos(cosTerm, 1.0_ReKi)) then
            sinTerm = 0.0_ReKi
         else
            sinTerm = sqrt(1.0_ReKi - cosTerm**2)
         end if
         
         dp      = u%p_plane(:,np+1,nt) - u%p_plane(:,np,nt)
         m%r_e(np,nt) = dot_product( u%xhat_plane(:,np  ,nt), dp )
         m%r_s(np,nt) = dot_product( u%xhat_plane(:,np+1,nt), dp )
         
         if (   sinTerm > ( max( m%r_e(np,nt), m%r_s(np,nt) ) / ( 100.0_ReKi*rmax ) ) ) then
            m%parallelFlag(np,nt) = .false.
            if ( u%D_wake(np,nt) > 0.0_ReKi ) then
               if ( m%r_e(np,nt) < rmax ) then
                  call SetErrStat( ErrID_Fatal, 'Radius to the wake center in the ending wake plane from the line where the starting and ending wake planes intersect for a given wake volume (volume='//trim(num2lstr(np))//',turbine='//trim(num2lstr(nt))//') is smaller than rmax: '//trim(num2lstr(rmax))//'.', errStat, errMsg, 'ComputeLocals' )
                  return
               end if
               if ( m%r_s(np,nt) < rmax ) then
                  call SetErrStat( ErrID_Fatal, 'Radius to the wake center in the starting wake plane from the line where the starting and ending wake planes intersect for a given wake volume (volume='//trim(num2lstr(np))//',turbine='//trim(num2lstr(nt))//') is smaller than rmax: '//trim(num2lstr(rmax))//'.', errStat, errMsg, 'ComputeLocals' )
                  return
               end if
            end if 
            m%r_e(np,nt) = m%r_e(np,nt) / sinTerm
            m%r_s(np,nt) = m%r_s(np,nt) / sinTerm
            m%rhat_s(:,np,nt) = (u%xhat_plane(:,np,nt)*cosTerm - u%xhat_plane(:,np+1,nt)        ) / sinTerm
            m%rhat_e(:,np,nt) = (u%xhat_plane(:,np,nt)         - u%xhat_plane(:,np+1,nt)*cosTerm) / sinTerm
            m%pvec_cs(:,np,nt) = u%p_plane(:,np  ,nt) - m%r_s(np,nt)*m%rhat_s(:,np,nt)
            m%pvec_ce(:,np,nt) = u%p_plane(:,np+1,nt) - m%r_e(np,nt)*m%rhat_e(:,np,nt)
         else
            m%parallelFlag(np,nt) = .true.
         end if
         
      end do
      
   end do
   
   
end subroutine ComputeLocals
!----------------------------------------------------------------------------------------------------------------------------------   
!> This function calculates jinc(x) = J_1(2*Pi*x)/x
real(ReKi) function jinc ( x )

   real(ReKi),      intent(in   ) :: x
   
   if ( EqualRealNos(x,0.0_ReKi) ) then
      jinc = Pi
   else
      jinc = BESSEL_JN( 1, TwoPi*x )/x
   end if
   
end function jinc
!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine LowResGridCalcOutput(n, u, p, y, m, errStat, errMsg)
   integer(IntKi),                 intent(in   )  :: n           !< Current simulation time increment (zero-based)
   type(AWAE_InputType),           intent(in   )  :: u           !< Inputs at Time t
   type(AWAE_ParameterType),       intent(in   )  :: p           !< Parameters
   type(AWAE_OutputType),          intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   type(AWAE_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                 intent(  out)  :: errStat     !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
   
   integer(IntKi)      :: nx, ny, nz, nt, np, nw, nx_low, ny_low, nz_low !< loop counters
   integer(IntKi)      :: nXYZ_low, n_wake       !< accumulating counters
   real(ReKi)          :: xhatBar_plane(3)       !< 
   real(ReKi)          :: x_end_plane
   real(ReKi)          :: x_start_plane
   real(ReKi)          :: r_vec_plane(3)
   real(ReKi)          :: tmp_xhatBar_plane
   real(ReKi)          :: r_tmp_plane
   real(ReKi)          :: D_wake_tmp
   real(ReKi)          :: Vx_wake_tmp
   real(ReKi)          :: Vr_wake_tmp(3)
   real(ReKi)          :: Vr_term(3)
   real(ReKi)          :: Vx_term
   real(ReKi)          :: Vsum_low(3)
   real(ReKi)          :: p_tmp_plane(3)
   real(ReKi)          :: tmp_vec(3)
   real(ReKi)          :: Vave_amb_low_norm
   real(ReKi)          :: delta, deltad
   real(ReKi)          :: wsum_tmp
   integer(IntKi)      :: ILo
   integer(IntKi)      :: maxPln
   character(*), parameter   :: RoutineName = 'LowResGridCalcOutput'   
   errStat = ErrID_None
   errMsg  = ""
   maxPln =  min(n,p%NumPlanes-2)
   
   !    ! read from file the ambient flow for the current time step
   ! call ReadLowResWindFile(n, p, m%Vamb_Low, errStat, errMsg)
   !    if ( errStat >= AbortErrLev ) then
   !       return
   !    end if
   
            
   nXYZ_low = 0
   m%N_wind(:,:) = 0
   
      ! Loop over the entire grid of low resolution ambient wind data to compute:
      !    1) the disturbed flow at each point and 2) the averaged disturbed velocity of each wake plane
   
   do nz_low=0, p%nZ_low-1 
      do ny_low=0, p%nY_low-1
         do nx_low=0, p%nX_low-1
            
               ! set the disturbed flow equal to the ambient flow for this time step
            m%Vdist_low(:,nx_low,ny_low,nz_low) = m%Vamb_low(:,nx_low,ny_low,nz_low)
            
            nXYZ_low = nXYZ_low + 1
            n_wake = 0
            xhatBar_plane = 0.0_ReKi
            
            do nt = 1,p%NumTurbines
                  
               x_end_plane = dot_product(u%xhat_plane(:,0,nt), (p%Grid_Low(:,nXYZ_low) - u%p_plane(:,0,nt)) )
               
               do np = 0, maxPln  
                  
                     ! Reset interpolation counter
                  ILo = 0
                  
                     ! Construct the endcaps of the current wake plane volume
                  x_start_plane = x_end_plane
                  x_end_plane = dot_product(u%xhat_plane(:,np+1,nt), (p%Grid_Low(:,nXYZ_low) - u%p_plane(:,np+1,nt)) )
                  
                     ! test if the point is within the endcaps of the wake volume
                  if ( ( x_start_plane >= 0.0_ReKi ) .and. ( x_end_plane < 0.0_ReKi ) ) then
                     
                     delta = x_start_plane / ( x_start_plane - x_end_plane )
                     deltad = (1.0_ReKi - delta)
                     if ( m%parallelFlag(np,nt) ) then
                        p_tmp_plane = delta*u%p_plane(:,np+1,nt) + deltad*u%p_plane(:,np,nt)
                     else
                        tmp_vec = delta*m%rhat_e(:,np,nt) + deltad*m%rhat_s(:,np,nt)
                        p_tmp_plane = delta*m%pvec_ce(:,np,nt) + deltad*m%pvec_cs(:,np,nt) + ( delta*m%r_e(np,nt) + deltad*m%r_s(np,nt) )* tmp_vec / TwoNorm(tmp_vec)
                     end if
                     
                        
                        
                     r_vec_plane = p%Grid_Low(:,nXYZ_low) - p_tmp_plane
                     r_tmp_plane = TwoNorm( r_vec_plane )
                     
                        ! test if the point is within radial finite-difference grid
                     if ( r_tmp_plane <= p%r(p%numRadii-1) ) then
                        
                        n_wake = n_wake + 1
                       ! m%r_plane(n_wake) = r_tmp_plane   ! Why do we need this??  GJH
                        
                        if ( EqualRealNos(r_tmp_plane, 0.0_ReKi) ) then         
                           m%rhat_plane(:,n_wake) = 0.0_ReKi
                        else
                           m%rhat_plane(:,n_wake) = ( r_vec_plane  ) / r_tmp_plane
                        end if
                        

                           ! given r_tmp_plane and Vx_wake at p%dr increments, find value of m%Vx_wake(@r_tmp_plane) using interpolation 
                        m%Vx_wake(n_wake) = delta*InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np+1,nt), ILo, p%NumRadii ) + deltad*InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np,nt), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                        m%Vr_wake(n_wake) = delta*InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np+1,nt), ILo, p%NumRadii ) + deltad*InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np,nt), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                        
                        
                        m%xhat_plane(:,n_wake) = delta*u%xhat_plane(:,np+1,nt) + deltad*u%xhat_plane(:,np,nt)
                        m%xhat_plane(:,n_wake) = m%xhat_plane(:,n_wake) / TwoNorm(m%xhat_plane(:,n_wake))
                        xhatBar_plane = xhatBar_plane + abs(m%Vx_wake(n_wake))*m%xhat_plane(:,n_wake)
      
                     end if  ! if the point is within radial finite-difference grid
                     
                        ! test if the point is within the radius of the wake volume cylinder                   
                     D_wake_tmp = delta*u%D_wake(np+1,nt) + deltad*u%D_wake(np,nt)        
                     if ( r_tmp_plane <= p%C_ScaleDiam*D_wake_tmp ) then
                        m%N_wind(np,nt) = m%N_wind(np,nt) + 1
                        
                        
                        if ( m%N_wind(np,nt) > p%n_wind_max ) then
                           call SetErrStat( ErrID_Fatal, 'The wake plane volume (plane='//trim(num2lstr(np))//',turbine='//trim(num2lstr(nt))//') contains more points than the maximum predicted points: 30*pi*DT(2*r*[Nr-1])**2/(dx*dy*dz)', errStat, errMsg, RoutineName )
                           return  ! if m%N_wind(np,nt) > p%n_wind_max then we will be indexing beyond the allocated memory for nx_wind,ny_wind,nz_wind arrays
                        end if
                        
                        select case ( p%Mod_Meander )
                        case (MeanderMod_Uniform) 
                           m%w(   m%N_wind(np,nt),np,nt) = 1.0_ReKi
                        case (MeanderMod_TruncJinc)  
                           m%w(   m%N_wind(np,nt),np,nt) = jinc( r_tmp_plane/( p%C_Meander*D_wake_tmp ) )
                        case (MeanderMod_WndwdJinc) 
                           m%w(   m%N_wind(np,nt),np,nt) = jinc( r_tmp_plane/( p%C_Meander*D_wake_tmp ) )*jinc( r_tmp_plane/( 2.0_ReKi*p%C_Meander*D_wake_tmp ) )
                        end select
                        
                        m%nx_wind(m%N_wind(np,nt),np,nt) = nx_low
                        m%ny_wind(m%N_wind(np,nt),np,nt) = ny_low
                        m%nz_wind(m%N_wind(np,nt),np,nt) = nz_low   
                        
                        if ( np == 0 ) then
                           if ( r_tmp_plane <= 0.5_ReKi*p%C_Meander*D_wake_tmp ) then
                              m%w_Amb(m%N_wind(np,nt),nt) = 1.0_ReKi
                           else
                              m%w_Amb(m%N_wind(np,nt),nt) = 0.0_ReKi
                           end if
                        end if
                                                
                     end if
                     exit
                  end if  ! if the point is within the endcaps of the wake volume                 
               end do     ! do np = 0, p%NumPlanes-2
            end do        ! do nt = 1,p%NumTurbines
            if (n_wake > 0) then
               tmp_xhatBar_plane = TwoNorm(xhatBar_plane)
               if ( EqualRealNos(tmp_xhatBar_plane, 0.0_ReKi) ) then
                  xhatBar_plane = 0.0_ReKi
               else
                  xhatBar_plane = xhatBar_plane / tmp_xhatBar_plane
               end if
               
               Vx_wake_tmp   = 0.0_ReKi
               Vr_wake_tmp   = 0.0_ReKi
               do nw = 1,n_wake 
                  Vr_term     = m%Vx_wake(nw)*m%xhat_plane(:,nw) + m%Vr_wake(nw)*m%rhat_plane(:,nw)
                  Vx_term     = dot_product( xhatBar_plane, Vr_term )
                  Vx_wake_tmp = Vx_wake_tmp + Vx_term*Vx_term
                  Vr_wake_tmp = Vr_wake_tmp + Vr_term
               end do
                  ! [I - XX']V = V - (V dot X)X
               Vr_wake_tmp = Vr_wake_tmp - dot_product(Vr_wake_tmp,xhatBar_plane)*xhatBar_plane               
               m%Vdist_low(:,nx_low,ny_low,nz_low) = m%Vdist_low(:,nx_low,ny_low,nz_low) + real(Vr_wake_tmp - xhatBar_plane*sqrt(Vx_wake_tmp),SiKi)
            end if  ! (n_wake > 0)
            
         end do ! do nx_low=0, p%nX_low-1
      end do    ! do ny_low=0, p%nY_low-1
   end do       ! do nz_low=0, p%nZ_low-1
   
   do nt = 1,p%NumTurbines
      if ( m%N_wind(0,nt) > 0 ) then
         
         Vsum_low  = 0.0_ReKi
         m%wsum(0) = 0.0_ReKi
      
         do nw=1,m%N_wind(0,nt)   
            Vsum_low  = Vsum_low  + m%w_Amb(nw,nt)*real(m%Vamb_Low(:, m%nx_wind(nw,0,nt), m%ny_wind(nw,0,nt), m%nz_wind(nw,0,nt)),ReKi)
            m%wsum(0) = m%wsum(0) + m%w_Amb(nw,nt)
         end do
         
         if ( EqualRealNos(m%wsum(0),0.0_ReKi) ) then
            call SetErrStat( ErrID_Fatal, 'The sum of the weightings for ambient spatial-averaging in the low-resolution domain associated with the wake volume at the rotor disk for turbine '//trim(num2lstr(nt))//' is zero.', errStat, errMsg, RoutineName )
            return     
         end if
         
         Vsum_low       = Vsum_low / m%wsum(0)  ! if N_wind gets large ( ~= 100,000 ) then this may not give enough precision in Vave_amb_low
         Vave_amb_low_norm  = TwoNorm(Vsum_low)
         if ( EqualRealNos(Vave_amb_low_norm,0.0_ReKi) ) then    
            call SetErrStat( ErrID_Fatal, 'The magnitude of the spatial-averaged ambient wind speed in the low-resolution domain associated with the wake volume at the rotor disk for turbine '//trim(num2lstr(nt))//' is zero.', errStat, errMsg, RoutineName )
            return     
         end if
      
         y%Vx_wind_disk(nt) = dot_product( u%xhat_plane(:,0,nt), Vsum_low )
         y%TI_amb(nt)       = 0.0_ReKi
         do nw=1,m%N_wind(0,nt)
            y%TI_amb(nt) = y%TI_amb(nt) + m%w_Amb(nw,nt)*TwoNorm( real(m%Vamb_Low(:, m%nx_wind(nw,0,nt), m%ny_wind(nw,0,nt), m%nz_wind(nw,0,nt)),ReKi) - Vsum_low )**2
         end do
         y%TI_amb(nt) = sqrt(y%TI_amb(nt)/(3.0*m%wsum(0)))/Vave_amb_low_norm
      else
         y%Vx_wind_disk(nt) = 0.0_ReKi
         y%TI_amb(nt)       = 0.0_ReKi 
      end if
      
      
      do np = 0, maxPln !p%NumPlanes-2
         if ( (u%D_wake(np,nt) > 0.0_ReKi) .and.  (m%N_wind(np,nt) < p%n_wind_min) ) then
            call SetErrStat( ErrID_Fatal, 'The number of points in the wake volume #'//trim(num2lstr(np))//' for turbine '//trim(num2lstr(nt))//' is '//trim(num2lstr(m%N_wind(np,nt)))//', which is less than the minimum threshold, '//trim(num2lstr(p%n_wind_min))//'.', errStat, errMsg, RoutineName )
            return     
         else if ( m%N_wind(np,nt) > 0  ) then            
            Vsum_low   = 0.0_ReKi
            m%wsum(np) = 0.0_ReKi
            do nw=1,m%N_wind(np,nt)   
               Vsum_low   = Vsum_low   + m%w(nw,np,nt)*m%Vdist_low( :, m%nx_wind(nw,np,nt),m%ny_wind(nw,np,nt),m%nz_wind(nw,np,nt) )
               m%wsum(np) = m%wsum(np) + m%w(nw,np,nt)
            end do
            y%V_plane(:,np,nt) = Vsum_low
         else
            y%V_plane(:,np,nt) = 0.0_ReKi
            m%wsum   (  np   ) = 0.0_ReKi
         end if
         
      end do
      
      if (  m%wsum(maxPln) > 0.0_ReKi ) then
         y%V_plane(:,maxPln+1,nt) =   y%V_plane(:,maxPln,nt)                          / m%wsum(maxPln)
      else
         y%V_plane(:,maxPln+1,nt) = 0.0_ReKi
      end if
      do np = maxPln, 1, -1
         wsum_tmp = m%wsum(np) + m%wsum(np-1)
         if ( wsum_tmp     > 0.0_ReKi ) then
            y%V_plane(:,np   ,nt) = ( y%V_plane(:,np    ,nt) + y%V_plane(:,np-1,nt) ) /   wsum_tmp
         end if
      end do
      if (  m%wsum(0     ) > 0.0_ReKi ) then
         y%V_plane(:,0       ,nt) =   y%V_plane(:,0     ,nt)                          / m%wsum(0     )
      end if
      
   end do

end subroutine LowResGridCalcOutput

#ifdef PARALLEL_CODE
!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine LowResGridCalcOutputOMP(n, u, p, y, Vdist_low, Vamb_low, m, errStat, errMsg)
   integer(IntKi),                 intent(in   )  :: n           !< Current simulation time increment (zero-based)
   type(AWAE_InputType),           intent(in   )  :: u           !< Inputs at Time t
   type(AWAE_ParameterType),       intent(in   )  :: p           !< Parameters
   type(AWAE_OutputType),          intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   real(SiKi),                     intent(inout)  :: Vdist_low(:,0:,0:,0:)
   real(SiKi),                     intent(inout)  :: Vamb_low(:,0:,0:,0:)
   type(AWAE_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                 intent(  out)  :: errStat     !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
   
   integer(IntKi)      :: i,nx, ny, nz, nt, np, nw, nx_low, ny_low, nz_low !< loop counters
   integer(IntKi)      :: nXYZ_low, n_wake       !< accumulating counters
   real(ReKi)          :: xhatBar_plane(3)       !< 
   real(ReKi)          :: x_end_plane
   real(ReKi)          :: x_start_plane
   real(ReKi)          :: r_vec_plane(3)
   real(ReKi)          :: tmp_xhatBar_plane
   real(ReKi)          :: r_tmp_plane
   real(ReKi)          :: D_wake_tmp
   real(ReKi)          :: Vx_wake_tmp
   real(ReKi)          :: Vr_wake_tmp(3)
   real(ReKi)          :: Vr_term(3)
   real(ReKi)          :: Vx_term
   real(ReKi)          :: Vsum_low(3)
   real(ReKi)          :: p_tmp_plane(3)
   real(ReKi)          :: tmp_vec(3)
   real(ReKi)          :: Vave_amb_low_norm
   real(ReKi)          :: delta, deltad
   real(ReKi)          :: wsum_tmp
   integer(IntKi)      :: ILo, j
   integer(IntKi)      :: maxPln, procs
   
   character(*), parameter   :: RoutineName = 'LowResGridCalcOutput'   
   errStat = ErrID_None
   errMsg  = ""
   maxPln =  min(n,p%NumPlanes-2)

   nXYZ_low = 0
   m%N_wind(:,:) = 0
   
      ! Loop over the entire grid of low resolution ambient wind data to compute:
      !    1) the disturbed flow at each point and 2) the averaged disturbed velocity of each wake plane
   procs = omp_get_num_procs()
   call omp_set_num_threads(procs)
  
   !$OMP PARALLEL DO PRIVATE(j, nx_low,ny_low,nz_low, nXYZ_low, n_wake, xhatBar_plane, x_end_plane, nt,np,ILo,x_start_plane,delta,deltad, p_tmp_plane,tmp_vec,r_vec_plane,r_tmp_plane,tmp_xhatBar_plane,Vx_wake_tmp,Vr_wake_tmp,nw,Vr_term,Vx_term) SHARED(Vdist_low,Vamb_low,m,u,p,maxPln,errStat, errMsg) DEFAULT(NONE)    
   do i = 0, p%nX_low*p%nY_low*p%nZ_low - 1
      
            nx_low = mod(i,p%nX_low)
            ny_low = mod(i/(p%nX_low),p%nY_low)
            nz_low = i / (p%nX_low*p%nY_low)

               ! set the disturbed flow equal to the ambient flow for this time step
            Vdist_low(:,nx_low,ny_low,nz_low) = Vamb_low(:,nx_low,ny_low,nz_low)

            nXYZ_low = i + 1
            n_wake = 0
            xhatBar_plane = 0.0_ReKi
            
            do nt = 1,p%NumTurbines
                  
               x_end_plane = dot_product(u%xhat_plane(:,0,nt), (p%Grid_Low(:,nXYZ_low) - u%p_plane(:,0,nt)) )
               
               do np = 0, maxPln  
                  
                     ! Reset interpolation counter
                  ILo = 0
                  
                     ! Construct the endcaps of the current wake plane volume
                  x_start_plane = x_end_plane
                  x_end_plane = dot_product(u%xhat_plane(:,np+1,nt), (p%Grid_Low(:,nXYZ_low) - u%p_plane(:,np+1,nt)) )
                  
                     ! test if the point is within the endcaps of the wake volume
                  if ( ( x_start_plane >= 0.0_ReKi ) .and. ( x_end_plane < 0.0_ReKi ) ) then
                     
                     delta = x_start_plane / ( x_start_plane - x_end_plane )
                     deltad = (1.0_ReKi - delta)
                     if ( m%parallelFlag(np,nt) ) then
                        p_tmp_plane = delta*u%p_plane(:,np+1,nt) + deltad*u%p_plane(:,np,nt)
                     else
                        tmp_vec = delta*m%rhat_e(:,np,nt) + deltad*m%rhat_s(:,np,nt)
                        p_tmp_plane = delta*m%pvec_ce(:,np,nt) + deltad*m%pvec_cs(:,np,nt) + ( delta*m%r_e(np,nt) + deltad*m%r_s(np,nt) )* tmp_vec / TwoNorm(tmp_vec)
                     end if
                     
                        
                        
                     r_vec_plane = p%Grid_Low(:,nXYZ_low) - p_tmp_plane
                     r_tmp_plane = TwoNorm( r_vec_plane )
                     
                        ! test if the point is within radial finite-difference grid
                     if ( r_tmp_plane <= p%r(p%numRadii-1) ) then
                        
                        n_wake = n_wake + 1
                       ! m%r_plane(n_wake) = r_tmp_plane   ! Why do we need this??  GJH
                        
                        if ( EqualRealNos(r_tmp_plane, 0.0_ReKi) ) then         
                           m%rhat_plane(:,n_wake) = 0.0_ReKi
                        else
                           m%rhat_plane(:,n_wake) = ( r_vec_plane  ) / r_tmp_plane
                        end if
                        

                           ! given r_tmp_plane and Vx_wake at p%dr increments, find value of m%Vx_wake(@r_tmp_plane) using interpolation 
                        m%Vx_wake(n_wake) = delta*InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np+1,nt), ILo, p%NumRadii ) + deltad*InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np,nt), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                        m%Vr_wake(n_wake) = delta*InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np+1,nt), ILo, p%NumRadii ) + deltad*InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np,nt), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                        
                        
                        m%xhat_plane(:,n_wake) = delta*u%xhat_plane(:,np+1,nt) + deltad*u%xhat_plane(:,np,nt)
                        m%xhat_plane(:,n_wake) = m%xhat_plane(:,n_wake) / TwoNorm(m%xhat_plane(:,n_wake))
                        xhatBar_plane = xhatBar_plane + abs(m%Vx_wake(n_wake))*m%xhat_plane(:,n_wake)
      
                     end if  ! if the point is within radial finite-difference grid
                     
                        ! test if the point is within the radius of the wake volume cylinder
                     D_wake_tmp = delta*u%D_wake(np+1,nt) + deltad*u%D_wake(np,nt)        
                     if ( r_tmp_plane <= p%C_ScaleDiam*D_wake_tmp ) then
                        m%N_wind(np,nt) = m%N_wind(np,nt) + 1
                        
                           
                        if ( m%N_wind(np,nt) > p%n_wind_max ) then
                           call SetErrStat( ErrID_Fatal, 'The wake plane volume (plane='//trim(num2lstr(np))//',turbine='//trim(num2lstr(nt))//') contains more points than the maximum predicted points: 30*pi*DT(2*r*[Nr-1])**2/(dx*dy*dz)', errStat, errMsg, 'LowResGridCalcOutput'  )
                           return  ! if m%N_wind(np,nt) > p%n_wind_max then we will be indexing beyond the allocated memory for nx_wind,ny_wind,nz_wind arrays
                        end if
                        
                        select case ( p%Mod_Meander )
                        case (MeanderMod_Uniform) 
                           m%w(   m%N_wind(np,nt),np,nt) = 1.0_ReKi
                        case (MeanderMod_TruncJinc)  
                           m%w(   m%N_wind(np,nt),np,nt) = jinc( r_tmp_plane/( p%C_Meander*D_wake_tmp ) )
                        case (MeanderMod_WndwdJinc) 
                           m%w(   m%N_wind(np,nt),np,nt) = jinc( r_tmp_plane/( p%C_Meander*D_wake_tmp ) )*jinc( r_tmp_plane/( 2.0_ReKi*p%C_Meander*D_wake_tmp ) )
                        end select
                        
                        m%nx_wind(m%N_wind(np,nt),np,nt) = nx_low
                        m%ny_wind(m%N_wind(np,nt),np,nt) = ny_low
                        m%nz_wind(m%N_wind(np,nt),np,nt) = nz_low   
                        
                        if ( np == 0 ) then
                           if ( r_tmp_plane <= 0.5_ReKi*p%C_Meander*D_wake_tmp ) then
                              m%w_Amb(m%N_wind(np,nt),nt) = 1.0_ReKi
                           else
                              m%w_Amb(m%N_wind(np,nt),nt) = 0.0_ReKi
                           end if
                        end if
                        
                     end if
                     exit
                  end if  ! if the point is within the endcaps of the wake volume                 
               end do     ! do np = 0, p%NumPlanes-2
            end do        ! do nt = 1,p%NumTurbines
            if (n_wake > 0) then
               tmp_xhatBar_plane = TwoNorm(xhatBar_plane)
               if ( EqualRealNos(tmp_xhatBar_plane, 0.0_ReKi) ) then
                  xhatBar_plane = 0.0_ReKi
               else
                  xhatBar_plane = xhatBar_plane / tmp_xhatBar_plane
               end if
               
               Vx_wake_tmp   = 0.0_ReKi
               Vr_wake_tmp   = 0.0_ReKi
               do nw = 1,n_wake 
                  Vr_term     = m%Vx_wake(nw)*m%xhat_plane(:,nw) + m%Vr_wake(nw)*m%rhat_plane(:,nw)
                  Vx_term     = dot_product( xhatBar_plane, Vr_term )
                  Vx_wake_tmp = Vx_wake_tmp + Vx_term*Vx_term
                  Vr_wake_tmp = Vr_wake_tmp + Vr_term
               end do
                  ! [I - XX']V = V - (V dot X)X
               Vr_wake_tmp = Vr_wake_tmp - dot_product(Vr_wake_tmp,xhatBar_plane)*xhatBar_plane               
               m%Vdist_low(:,nx_low,ny_low,nz_low) = m%Vdist_low(:,nx_low,ny_low,nz_low) + real(Vr_wake_tmp - xhatBar_plane*sqrt(Vx_wake_tmp),SiKi)
            end if  ! (n_wake > 0)
            
   end do ! do i : loop over all grid points
   !$OMP END PARALLEL DO
   
   do nt = 1,p%NumTurbines
      if ( m%N_wind(0,nt) > 0 ) then
         
         Vsum_low  = 0.0_ReKi
         m%wsum(0) = 0.0_ReKi
      
         do nw=1,m%N_wind(0,nt)   
            Vsum_low  = Vsum_low  + m%w_Amb(nw,nt)*real(m%Vamb_Low(:, m%nx_wind(nw,0,nt), m%ny_wind(nw,0,nt), m%nz_wind(nw,0,nt)),ReKi)
            m%wsum(0) = m%wsum(0) + m%w_Amb(nw,nt)
         end do
         
         if ( EqualRealNos(m%wsum(0),0.0_ReKi) ) then
            call SetErrStat( ErrID_Fatal, 'The sum of the weightings for ambient spatial-averaging in the low-resolution domain associated with the wake volume at the rotor disk for turbine '//trim(num2lstr(nt))//' is zero.', errStat, errMsg, RoutineName )
            return     
         end if
         
         Vsum_low       = Vsum_low / m%wsum(0)  ! if N_wind gets large ( ~= 100,000 ) then this may not give enough precision in Vave_amb_low
         Vave_amb_low_norm  = TwoNorm(Vsum_low)
         if ( EqualRealNos(Vave_amb_low_norm,0.0_ReKi) ) then    
            call SetErrStat( ErrID_Fatal, 'The magnitude of the spatial-averaged ambient wind speed in the low-resolution domain associated with the wake volume at the rotor disk for turbine '//trim(num2lstr(nt))//' is zero.', errStat, errMsg, RoutineName )
            return     
         end if
      
         y%Vx_wind_disk(nt) = dot_product( u%xhat_plane(:,0,nt), Vsum_low )
         y%TI_amb(nt)       = 0.0_ReKi
         do nw=1,m%N_wind(0,nt)
            y%TI_amb(nt) = y%TI_amb(nt) + m%w_Amb(nw,nt)*TwoNorm( real(m%Vamb_Low(:, m%nx_wind(nw,0,nt), m%ny_wind(nw,0,nt), m%nz_wind(nw,0,nt)),ReKi) - Vsum_low )**2
         end do
         y%TI_amb(nt) = sqrt(y%TI_amb(nt)/(3.0*m%wsum(0)))/Vave_amb_low_norm
      else
         y%Vx_wind_disk(nt) = 0.0_ReKi
         y%TI_amb(nt)       = 0.0_ReKi 
      end if
      
      
      do np = 0, maxPln !p%NumPlanes-2
         if ( (u%D_wake(np,nt) > 0.0_ReKi) .and.  (m%N_wind(np,nt) < p%n_wind_min) ) then
            call SetErrStat( ErrID_Fatal, 'The number of points in the wake volume #'//trim(num2lstr(np))//' for turbine '//trim(num2lstr(nt))//' is '//trim(num2lstr(m%N_wind(np,nt)))//', which is less than the minimum threshold, '//trim(num2lstr(p%n_wind_min))//'.', errStat, errMsg, RoutineName )
            return     
         else if ( m%N_wind(np,nt) > 0  ) then            
            Vsum_low   = 0.0_ReKi
            m%wsum(np) = 0.0_ReKi
            do nw=1,m%N_wind(np,nt)   
               Vsum_low   = Vsum_low   + m%w(nw,np,nt)*m%Vdist_low( :, m%nx_wind(nw,np,nt),m%ny_wind(nw,np,nt),m%nz_wind(nw,np,nt) )
               m%wsum(np) = m%wsum(np) + m%w(nw,np,nt)
            end do
            y%V_plane(:,np,nt) = Vsum_low
         else
            y%V_plane(:,np,nt) = 0.0_ReKi
            m%wsum   (  np   ) = 0.0_ReKi
         end if
         
      end do
      
      if (  m%wsum(maxPln) > 0.0_ReKi ) then
         y%V_plane(:,maxPln+1,nt) =   y%V_plane(:,maxPln,nt)                          / m%wsum(maxPln)
      else
         y%V_plane(:,maxPln+1,nt) = 0.0_ReKi
      end if
      do np = maxPln, 1, -1
         wsum_tmp = m%wsum(np) + m%wsum(np-1)
         if ( wsum_tmp     > 0.0_ReKi ) then
            y%V_plane(:,np   ,nt) = ( y%V_plane(:,np    ,nt) + y%V_plane(:,np-1,nt) ) /   wsum_tmp
         end if
      end do
      if (  m%wsum(0     ) > 0.0_ReKi ) then
         y%V_plane(:,0       ,nt) =   y%V_plane(:,0     ,nt)                          / m%wsum(0     )
      end if
      
   end do
   
end subroutine LowResGridCalcOutputOMP
#endif
!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine HighResGridCalcOutput(n, u, p, y, m, errStat, errMsg)
   integer(IntKi),                 intent(in   )  :: n           !< Current high-res, simulation time increment (zero-based)
   type(AWAE_InputType),           intent(in   )  :: u           !< Inputs at Time t
   type(AWAE_ParameterType),       intent(in   )  :: p           !< Parameters
   type(AWAE_OutputType),          intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   type(AWAE_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                 intent(  out)  :: errStat     !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
   
   integer(IntKi)      :: nx, ny, nz, nt, nt2, np, nw, nx_high, ny_high, nz_high, n_hl !< loop counters
   integer(IntKi)      :: nXYZ_high, n_wake       !< accumulating counters
   real(ReKi)          :: xhatBar_plane(3)       !< 
   real(ReKi)          :: tmp_xhatBar_plane
   real(ReKi)          :: x_end_plane
   real(ReKi)          :: x_start_plane
   real(ReKi)          :: r_vec_plane(3)
   real(ReKi)          :: r_tmp_plane
   real(ReKi)          :: Vx_wake_tmp
   real(ReKi)          :: Vr_wake_tmp(3)
   real(ReKi)          :: Vr_term(3)
   real(ReKi)          :: Vx_term
   real(ReKi)          :: Vsum_low(3)
   real(ReKi)          :: p_tmp_plane(3)
   real(ReKi)          :: tmp_vec(3)
   real(ReKi)          :: delta, deltad
   integer(IntKi)      :: ILo
   integer(IntKi)      :: maxPln
   integer(IntKi)      :: n_high_low
   character(*), parameter   :: RoutineName = 'HighResGridCalcOutput'
   errStat = ErrID_None
   errMsg  = ""

   
   maxPln =  min(n,p%NumPlanes-2)
     
      ! We only need one high res file for that last simulation time
   if ( (n/p%n_high_low) == (p%NumDT-1) ) then
      n_high_low = 1
   else
      n_high_low = p%n_high_low
   end if
   
   
      ! Loop over the entire grid of low resolution ambient wind data to compute:
      !    1) the disturbed flow at each point and 2) the averaged disturbed velocity of each wake plane
   

   do nt = 1,p%NumTurbines
      nXYZ_high = 0
      
            ! set the disturbed flow equal to the ambient flow for this time step
      y%Vdist_high(nt)%data = m%Vamb_high(nt)%data
      
      do nz_high=0, p%nZ_high-1 
         do ny_high=0, p%nY_high-1
            do nx_high=0, p%nX_high-1
              
               nXYZ_high = nXYZ_high + 1
               n_wake = 0
               xhatBar_plane = 0.0_ReKi

               do nt2 = 1,p%NumTurbines
                  if (nt /= nt2) then  
                     
                     x_end_plane = dot_product(u%xhat_plane(:,0,nt2), (p%Grid_high(:,nXYZ_high,nt) - u%p_plane(:,0,nt2)) )
               
                     do np = 0, maxPln !p%NumPlanes-2
                  
                           ! Reset interpolation counter
                        ILo = 0
                  
                           ! Construct the endcaps of the current wake plane volume
                        x_start_plane = x_end_plane
                        x_end_plane = dot_product(u%xhat_plane(:,np+1,nt2), (p%Grid_high(:,nXYZ_high,nt) - u%p_plane(:,np+1,nt2)) )
                  
                           ! test if the point is within the endcaps of the wake volume
                        if ( ( x_start_plane >= 0.0_ReKi ) .and. ( x_end_plane < 0.0_ReKi ) ) then
                           
                           delta = x_start_plane / ( x_start_plane - x_end_plane )
                           deltad = (1.0_ReKi - delta)
                           if ( m%parallelFlag(np,nt2) ) then
                              p_tmp_plane = delta*u%p_plane(:,np+1,nt2) + deltad*u%p_plane(:,np,nt2)
                           else
                              tmp_vec  = delta*m%rhat_e(:,np,nt2) + deltad*m%rhat_s(:,np,nt2)
                              p_tmp_plane = delta*m%pvec_ce(:,np,nt2) + deltad*m%pvec_cs(:,np,nt2) + ( delta*m%r_e(np,nt2) + deltad*m%r_s(np,nt2) )* tmp_vec / TwoNorm(tmp_vec)
                           end if
                     
                           r_vec_plane = p%Grid_high(:,nXYZ_high,nt) - p_tmp_plane
                           r_tmp_plane = TwoNorm( r_vec_plane )
                     
                              ! test if the point is within radial finite-difference grid
                           if ( r_tmp_plane <= p%r(p%numRadii-1) ) then
                        
                              n_wake = n_wake + 1
                             ! m%r_plane(n_wake) = r_tmp_plane   ! Why do we need this??  GJH
                        
                              if ( EqualRealNos(r_tmp_plane, 0.0_ReKi) ) then         
                                 m%rhat_plane(:,n_wake) = 0.0_ReKi
                              else
                                 m%rhat_plane(:,n_wake) = ( r_vec_plane  ) / r_tmp_plane
                              end if
                        
                                                 

                           ! given r_tmp_plane and Vx_wake at p%dr increments, find value of m%Vx_wake(@r_tmp_plane) using interpolation 
                              m%Vx_wake(n_wake) = delta*InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np+1,nt2), ILo, p%NumRadii ) + deltad*InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np,nt2), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                              m%Vr_wake(n_wake) = delta*InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np+1,nt2), ILo, p%NumRadii ) + deltad*InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np,nt2), ILo, p%NumRadii )!( XVal, XAry, YAry, ILo, AryLen )
                                              
                              m%xhat_plane(:,n_wake) = delta*u%xhat_plane(:,np+1,nt2) + deltad*u%xhat_plane(:,np,nt2)
                              m%xhat_plane(:,n_wake) = m%xhat_plane(:,n_wake) / TwoNorm(m%xhat_plane(:,n_wake))
                              xhatBar_plane = xhatBar_plane + abs(m%Vx_wake(n_wake))*m%xhat_plane(:,n_wake)    
                              
                           end if  ! if the point is within radial finite-difference grid
 
                           exit
                        end if  ! if the point is within the endcaps of the wake volume                 
                     end do     ! np = 0, p%NumPlanes-2
                  end if    ! nt /= nt2          
               end do        ! nt2 = 1,p%NumTurbines
               if (n_wake > 0) then
                  
                  tmp_xhatBar_plane = TwoNorm(xhatBar_plane)
                  if ( EqualRealNos(tmp_xhatBar_plane, 0.0_ReKi) ) then
                     xhatBar_plane = 0.0_ReKi
                  else
                     xhatBar_plane = xhatBar_plane / tmp_xhatBar_plane
                  end if
                  
                  Vx_wake_tmp   = 0.0_ReKi
                  Vr_wake_tmp   = 0.0_ReKi
                  do nw = 1,n_wake 
                     Vr_term     = m%Vx_wake(nw)*m%xhat_plane(:,nw) + m%Vr_wake(nw)*m%rhat_plane(:,nw)
                     Vx_term     = dot_product( xhatBar_plane, Vr_term )
                     Vx_wake_tmp = Vx_wake_tmp + Vx_term*Vx_term
                     Vr_wake_tmp = Vr_wake_tmp + Vr_term
                  end do
                     ! [I - XX']V = V - (V dot X)X
                  Vr_wake_tmp = Vr_wake_tmp - dot_product(Vr_wake_tmp,xhatBar_plane)*xhatBar_plane 
                  do n_hl=0, n_high_low-1
                     y%Vdist_high(nt)%data(:,nx_high,ny_high,nz_high,n_hl) = y%Vdist_high(nt)%data(:,nx_high,ny_high,nz_high,n_hl) + real(Vr_wake_tmp - xhatBar_plane*sqrt(Vx_wake_tmp),SiKi)
                  end do
                     
               end if  ! (n_wake > 0)
            
            end do ! nx_high=0, p%nX_high-1
         end do    ! ny_high=0, p%nY_high-1
      end do       ! nz_high=0, p%nZ_high-1
   end do          ! nt = 1,p%NumTurbines
   
   
end subroutine HighResGridCalcOutput
  


!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine AWAE_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, errStat, errMsg )
!..................................................................................................................................

   type(AWAE_InitInputType),       intent(in   ) :: InitInp       !< Input data for initialization routine
   type(AWAE_InputType),           intent(  out) :: u             !< An initial guess for the input; input mesh must be defined
   type(AWAE_ParameterType),       intent(  out) :: p             !< Parameters
   type(AWAE_ContinuousStateType), intent(  out) :: x             !< Initial continuous states
   type(AWAE_DiscreteStateType),   intent(  out) :: xd            !< Initial discrete states
   type(AWAE_ConstraintStateType), intent(  out) :: z             !< Initial guess of the constraint states
   type(AWAE_OtherStateType),      intent(  out) :: OtherState    !< Initial other states
   type(AWAE_OutputType),          intent(  out) :: y             !< Initial system outputs (outputs are not calculated;
                                                                  !!   only the output mesh is initialized)
   type(AWAE_MiscVarType),         intent(  out) :: m             !< Initial misc/optimization variables
   real(DbKi),                     intent(in   ) :: interval      !< Low-resolution (FAST.Farm driver/glue code) time step, s 
   type(AWAE_InitOutputType),      intent(  out) :: InitOut       !< Output for initialization routine
   integer(IntKi),                 intent(  out) :: errStat       !< Error status of the operation
   character(*),                   intent(  out) :: errMsg        !< Error message if errStat /= ErrID_None
   

      ! Local variables
   integer(IntKi)                                :: i,j             ! loop counter
   real(ReKi)                                    :: gridLoc       ! Location of requested output slice in grid coordinates [0,sz-1]                                      
   integer(IntKi)                                :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                          :: errMsg2       ! temporary error message                                           
   character(*), parameter                       :: RoutineName = 'AWAE_Init'
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
      ! Initialize the NWTC Subroutine Library

   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   call DispNVD( AWAE_Ver )
  
   p%OutFileRoot  = TRIM(InitInp%OutFileRoot)
   
   
   
      ! Validate the initialization inputs
   call ValidateInitInputData( InitInp%InputFileData, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName ) 
      if (errStat >= AbortErrLev) then
         return
      end if
      
      !............................................................................................
      ! Define parameters
      !............................................................................................
      
   
      
      ! set the rest of the parameters         
   p%NumPlanes        = InitInp%InputFileData%NumPlanes   
   p%NumRadii         = InitInp%InputFileData%NumRadii    
   p%NumTurbines      = InitInp%InputFileData%NumTurbines 
   p%WindFilePath     = InitInp%InputFileData%WindFilePath ! TODO: Make sure this wasn't specified with the trailing folder separator. Note: on Windows a trailing / or \ causes no problem! GJH 
   p%n_high_low       = InitInp%n_high_low
   p%dt               = InitInp%InputFileData%dt
   p%NumDT            = InitInp%NumDT
   p%NOutDisWindXY    = InitInp%InputFileData%NOutDisWindXY
   p%NOutDisWindYZ    = InitInp%InputFileData%NOutDisWindYZ
   p%NOutDisWindXZ    = InitInp%InputFileData%NOutDisWindXZ
   p%WrDisWind        = InitInp%InputFileData%WrDisWind
   p%WrDisSkp1        = nint(InitInp%InputFileData%WrDisDT / p%dt)
   p%Mod_Meander      = InitInp%InputFileData%Mod_Meander
   p%C_Meander        = InitInp%InputFileData%C_Meander
   select case ( p%Mod_Meander )
   case (MeanderMod_Uniform) 
      p%C_ScaleDiam   = 0.5_ReKi*p%C_Meander
   case (MeanderMod_TruncJinc)  
      p%C_ScaleDiam   = 0.5_ReKi*p%C_Meander*1.21967_ReKi
   case (MeanderMod_WndwdJinc) 
      p%C_ScaleDiam   = 0.5_ReKi*p%C_Meander*2.23313_ReKi
   end select

      
   call allocAry( p%OutDisWindZ, p%NOutDisWindXY, "OutDisWindZ", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         RETURN        
      end if
      
   call allocAry( p%OutDisWindX, p%NOutDisWindYZ, "OutDisWindX", ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            RETURN        
         end if
         
   call allocAry( p%OutDisWindY, p%NOutDisWindXZ, "OutDisWindY", ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            RETURN        
         end if      
         
   p%OutDisWindZ = InitInp%InputFileData%OutDisWindZ      
   p%OutDisWindX = InitInp%InputFileData%OutDisWindX
   p%OutDisWindY = InitInp%InputFileData%OutDisWindY
   
   allocate( p%r(0:p%NumRadii-1),stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for p%r.', errStat, errMsg, RoutineName )
         return
      end if
    
   do i = 0,p%NumRadii-1
      p%r(i)       = InitInp%InputFileData%dr*i     
   end do
   
   
      ! Obtain the precursor grid information by parsing the necessary input files
      ! This will establish certain parameters as well as all of the initialization outputs
      ! Sets:
      ! Parameters: nX_low, nY_low, nZ_low, nX_high, nY_high, nZ_high, Grid_low, 
      !             Grid_high, n_high_low, n_wind_max, n_wind_min
      ! InitOutput: X0_high, Y0_high, Z0_high, dX_high, dY_high, dZ_high, nX_high, nY_high, nZ_high
   

   call AWAE_IO_InitGridInfo(InitInp, p, InitOut, errStat2, errMsg2)
      call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if (errStat2 >= AbortErrLev) then      
         return
   end if
   
   InitOut%Ver = AWAE_Ver
   
      ! Test the request output wind locations against grid information
   
         ! XY plane slices
      do i = 1,p%NOutDisWindXY        
         gridLoc = (p%OutDisWindZ(i) - p%Z0_low) / p%dZ_low
         if ( ( gridLoc < 0.0_ReKi ) .or. ( gridLoc > real(p%nZ_low-1, ReKi) ) ) then
            call SetErrStat(ErrID_Fatal, "The requested low-resolution XY output slice location, Z="//TRIM(Num2LStr(p%OutDisWindZ(i)))//", is outside of the low-resolution grid.", errStat, errMsg, RoutineName )              
         end if
      end do
      
         ! XZ plane slices
      do i = 1,p%NOutDisWindXZ        
         gridLoc = (p%OutDisWindY(i) - p%Y0_low) / p%dY_low
         if ( ( gridLoc < 0.0_ReKi ) .or. ( gridLoc > real(p%nY_low-1, ReKi) ) ) then
            call SetErrStat(ErrID_Fatal, "The requested low-resolution XZ output slice location, Y="//TRIM(Num2LStr(p%OutDisWindY(i)))//", is outside of the low-resolution grid.", errStat, errMsg, RoutineName )
         end if
      end do
      
         ! XZ plane slices
      do i = 1,p%NOutDisWindYZ        
         gridLoc = (p%OutDisWindX(i) - p%X0_low) / p%dX_low
         if ( ( gridLoc < 0.0_ReKi ) .or. ( gridLoc > real(p%nX_low-1, ReKi) ) ) then
            call SetErrStat(ErrID_Fatal, "The requested low-resolution YZ output slice location, X="//TRIM(Num2LStr(p%OutDisWindX(i)))//", is outside of the low-resolution grid.", errStat, errMsg, RoutineName )
         end if
      end do
      if (errStat2 >= AbortErrLev) then      
         return
      end if
  
   
   !interval = InitOut%dt
   
      !............................................................................................
      ! Define and initialize inputs here 
      !............................................................................................
   
   allocate ( u%xhat_plane(3,0:p%NumPlanes-1,1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for u%xhat_plane.', errStat, errMsg, RoutineName )     
   allocate ( u%p_plane   (3,0:p%NumPlanes-1,1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for u%p_plane.', errStat, errMsg, RoutineName )  
   allocate ( u%Vx_wake   (0:p%NumRadii-1,0:p%NumPlanes-1,1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for u%Vx_wake.', errStat, errMsg, RoutineName )  
   allocate ( u%Vr_wake   (0:p%NumRadii-1,0:p%NumPlanes-1,1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for u%Vr_wake.', errStat, errMsg, RoutineName )  
   allocate ( u%D_wake    (0:p%NumPlanes-1,1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for u%D_wake.', errStat, errMsg, RoutineName )  
   if (errStat /= ErrID_None) return
   

         
      
      !............................................................................................
      ! Define outputs here
      !............................................................................................

   allocate ( y%V_plane(3,0:p%NumPlanes-1,1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%V_plane.', errStat, errMsg, RoutineName )     
   allocate ( y%Vdist_High(1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%Vdist_High.', errStat, errMsg, RoutineName )  
   do i = 1, p%NumTurbines
      allocate ( y%Vdist_High(i)%data(3,0:p%nX_high-1,0:p%nY_high-1,0:p%nZ_high-1,0:p%n_high_low-1), STAT=ErrStat2 )
         if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%Vdist_High%data.', errStat, errMsg, RoutineName ) 
      y%Vdist_High(i)%data    = 0.0_Siki
   end do
   
   allocate ( y%Vx_wind_disk   (1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%Vx_rel_disk.', errStat, errMsg, RoutineName )  
   allocate ( y%TI_amb   (1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%TI_amb.', errStat, errMsg, RoutineName )  
   if (errStat /= ErrID_None) return
   
      ! This next step is not strictly necessary
   y%V_plane       = 0.0_Reki
   y%Vx_wind_disk  = 0.0_Reki
   y%TI_amb        = 0.0_Reki
   
   
   if ( p%NOutDisWindXY > 0 ) then
      ALLOCATE ( m%OutVizXYPlane(3,p%nX_low, p%nY_low,1) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error allocating memory for the Fast.Farm OutVizXYPlane arrays.'
         RETURN
      ENDIF   
   end if
   if ( p%NOutDisWindYZ > 0 ) then
      ALLOCATE ( m%OutVizYZPlane(3,p%nY_low, p%nZ_low,1) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error allocating memory for the Fast.Farm OutVizYZPlane arrays.'
         RETURN
      ENDIF   
   end if
   if ( p%NOutDisWindXZ > 0 ) then
      ALLOCATE ( m%OutVizXZPlane(3,p%nX_low, p%nZ_low,1) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error allocating memory for the Fast.Farm OutVizXZPlane arrays.'
         RETURN
      ENDIF   
   end if
      !............................................................................................
      ! Initialize misc vars : Note these are not the correct initializations because
      ! that would require valid input data, which we do not have here.  Instead we will check for
      ! an firstPass flag on the miscVars and if it is false we will properly initialize these state
      ! in CalcOutput or UpdateStates, as necessary.
      !............................................................................................
      
  
      
   
      ! miscvars to avoid the allocation per timestep
 
   allocate ( m%Vamb_low   ( 3, 0:p%nX_low-1 , 0:p%nY_low-1 , 0:p%nZ_low-1 )                  , STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vamb_low.', errStat, errMsg, RoutineName )   
   allocate ( m%Vdist_low  ( 3, 0:p%nX_low-1 , 0:p%nY_low-1 , 0:p%nZ_low-1 )                  , STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vdist_low.', errStat, errMsg, RoutineName ) 
      
      
   !allocate ( m%Vamb_high  ( 3, 0:p%nX_high-1, 0:p%nY_high-1, 0:p%nZ_high-1 ), STAT=errStat2 ) 
   !   if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vamb_High.', errStat, errMsg, RoutineName )  
      
   allocate ( m%Vamb_high(1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vamb_high.', errStat, errMsg, RoutineName )  
   do i = 1, p%NumTurbines
         allocate ( m%Vamb_high(i)%data(3,0:p%nX_high-1,0:p%nY_high-1,0:p%nZ_high-1,0:p%n_high_low-1), STAT=ErrStat2 )
            if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vamb_high%data.', errStat, errMsg, RoutineName ) 
   end do   
      
   allocate ( m%N_wind     ( 0:p%NumPlanes-2,1:p%NumTurbines ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%N_wind.', errStat, errMsg, RoutineName )  
   allocate ( m%xhat_plane ( 3, 1:p%NumTurbines ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%xhat_plane.', errStat, errMsg, RoutineName )  
   allocate ( m%rhat_plane ( 3, 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%rhat_plane.', errStat, errMsg, RoutineName )  
   allocate ( m%Vx_wake    ( 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vx_wake.', errStat, errMsg, RoutineName )  
   allocate ( m%Vr_wake    ( 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vr_wake.', errStat, errMsg, RoutineName )  
      
      
   allocate ( m%w          ( p%n_wind_max, 0:p%NumPlanes-2, 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%w.', errStat, errMsg, RoutineName )  
   allocate ( m%w_Amb      ( p%n_wind_max, 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%w_Amb.', errStat, errMsg, RoutineName )  
   allocate ( m%wsum       ( 0:p%NumPlanes-2 ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%wsum.', errStat, errMsg, RoutineName )  

   allocate ( m%nx_wind    ( p%n_wind_max, 0:p%NumPlanes-2, 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%nx_wind.', errStat, errMsg, RoutineName )  
   allocate ( m%ny_wind    ( p%n_wind_max, 0:p%NumPlanes-2, 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%ny_wind.', errStat, errMsg, RoutineName )  
   allocate ( m%nz_wind    ( p%n_wind_max, 0:p%NumPlanes-2, 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%nz_wind.', errStat, errMsg, RoutineName )  
   
   allocate ( m%parallelFlag( 0:p%NumPlanes-2,1:p%NumTurbines ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%parallelFlag.', errStat, errMsg, RoutineName )
   allocate ( m%r_s( 0:p%NumPlanes-2,1:p%NumTurbines ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%r_s.', errStat, errMsg, RoutineName )
   allocate ( m%r_e( 0:p%NumPlanes-2,1:p%NumTurbines ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%r_e.', errStat, errMsg, RoutineName )
   allocate ( m%rhat_s( 3,0:p%NumPlanes-2,1:p%NumTurbines ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%rhat_s.', errStat, errMsg, RoutineName )
   allocate ( m%rhat_e( 3,0:p%NumPlanes-2,1:p%NumTurbines ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%rhat_e.', errStat, errMsg, RoutineName )
   allocate ( m%pvec_cs( 3,0:p%NumPlanes-2,1:p%NumTurbines ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%pvec_cs.', errStat, errMsg, RoutineName )
   allocate ( m%pvec_ce( 3,0:p%NumPlanes-2,1:p%NumTurbines ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%pvec_ce.', errStat, errMsg, RoutineName )
   if (errStat /= ErrID_None) return

   
   ! Read-in the ambient wind data for the initial calculate output
   
   call AWAE_UpdateStates( 0.0_DbKi, -1, u, p, x, xd, z, OtherState, m, errStat, errMsg )
   
  
      


end subroutine AWAE_Init

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine AWAE_End( u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
!..................................................................................................................................

      type(AWAE_InputType),           intent(inout)  :: u           !< System inputs
      type(AWAE_ParameterType),       intent(inout)  :: p           !< Parameters
      type(AWAE_ContinuousStateType), intent(inout)  :: x           !< Continuous states
      type(AWAE_DiscreteStateType),   intent(inout)  :: xd          !< Discrete states
      type(AWAE_ConstraintStateType), intent(inout)  :: z           !< Constraint states
      type(AWAE_OtherStateType),      intent(inout)  :: OtherState  !< Other states
      type(AWAE_OutputType),          intent(inout)  :: y           !< System outputs
      type(AWAE_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
      integer(IntKi),               intent(  out)  :: errStat     !< Error status of the operation
      character(*),                 intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None



         ! Initialize errStat

      errStat = ErrID_None
      errMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:


         ! Destroy the input data:

      call AWAE_DestroyInput( u, errStat, errMsg )


         ! Destroy the parameter data:

      call AWAE_DestroyParam( p, errStat, errMsg )


         ! Destroy the state data:

      call AWAE_DestroyContState(   x,           errStat, errMsg )
      call AWAE_DestroyDiscState(   xd,          errStat, errMsg )
      call AWAE_DestroyConstrState( z,           errStat, errMsg )
      call AWAE_DestroyOtherState(  OtherState,  errStat, errMsg )
      call AWAE_DestroyMisc(        m,           errStat, errMsg ) 

         ! Destroy the output data:

      call AWAE_DestroyOutput( y, errStat, errMsg )




end subroutine AWAE_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine AWAE_UpdateStates( t, n, u, p, x, xd, z, OtherState, m, errStat, errMsg )
!..................................................................................................................................

   real(DbKi),                     intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),                 intent(in   ) :: n          !< Current simulation time step n = 0,1,...
   type(AWAE_InputType),             intent(inout) :: u          !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
  ! real(DbKi),                     intent(in   ) :: utimes   !< Times associated with u(:), in seconds
   type(AWAE_ParameterType),         intent(in   ) :: p          !< Parameters
   type(AWAE_ContinuousStateType),   intent(inout) :: x          !< Input: Continuous states at t;
                                                               !!   Output: Continuous states at t + Interval
   type(AWAE_DiscreteStateType),     intent(inout) :: xd         !< Input: Discrete states at t;
                                                               !!   Output: Discrete states at t  + Interval
   type(AWAE_ConstraintStateType),   intent(inout) :: z          !< Input: Constraint states at t;
                                                               !!   Output: Constraint states at t+dt
   type(AWAE_OtherStateType),        intent(inout) :: OtherState !< Input: Other states at t;
                                                               !!   Output: Other states at t+dt
   type(AWAE_MiscVarType),           intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),                 intent(  out) :: errStat    !< Error status of the operation
   character(*),                   intent(  out) :: errMsg     !< Error message if errStat /= ErrID_None

   ! local variables
   type(AWAE_InputType)                           :: uInterp           ! Interpolated/Extrapolated input
   integer(intKi)                               :: errStat2          ! temporary Error status
   character(ErrMsgLen)                         :: errMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'AWAE_UpdateStates'
   real(DbKi)          :: t1, tread
   integer(IntKi)                               :: n_high_low, nt, n_hl
   
   errStat = ErrID_None
   errMsg  = ""
   
   ! Read the ambient wind data that is needed for t+dt, i.e., n+1
   
   !t1 = omp_get_wtime()   
      ! read from file the ambient flow for the n+1 time step
   call ReadLowResWindFile(n+1, p, m%Vamb_Low, errStat, errMsg)
      if ( errStat >= AbortErrLev ) then
         return
      end if
  ! tread =  ( omp_get_wtime() - t1 )       
  ! write(*,*) 'Time spent reading data:  ',tread            
   
      
   if ( (n+1) == (p%NumDT-1) ) then
      n_high_low = 1
   else
      n_high_low = p%n_high_low
   end if
 
   do nt = 1,p%NumTurbines
      do n_hl=0, n_high_low-1
            ! read from file the ambient flow for the current time step
         call ReadHighResWindFile(nt, (n+1)*p%n_high_low + n_hl, p, m%Vamb_high(nt)%data(:,:,:,:,n_hl), errStat, errMsg)
            if ( errStat >= AbortErrLev ) then
               return
            end if 
      end do
   end do
   
end subroutine AWAE_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine AWAE_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
! NOTE: no matter how many channels are selected for output, all of the outputs are calcalated
! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   real(DbKi),                     intent(in   )  :: t           !< Current simulation time in seconds
   type(AWAE_InputType),           intent(in   )  :: u           !< Inputs at Time t
   type(AWAE_ParameterType),       intent(in   )  :: p           !< Parameters
   type(AWAE_ContinuousStateType), intent(in   )  :: x           !< Continuous states at t
   type(AWAE_DiscreteStateType),   intent(in   )  :: xd          !< Discrete states at t
   type(AWAE_ConstraintStateType), intent(in   )  :: z           !< Constraint states at t
   type(AWAE_OtherStateType),      intent(in   )  :: OtherState  !< Other states at t
   type(AWAE_OutputType),          intent(inout)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   type(AWAE_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
   integer(IntKi),                 intent(  out)  :: errStat     !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None
   

   integer, parameter                           :: indx = 1  
   integer(intKi)                               :: i, j, k
   integer(intKi)                               :: errStat2
   character(ErrMsgLen)                         :: errMsg2
   character(*), parameter                      :: RoutineName = 'AWAE_CalcOutput'
   integer(intKi)                               :: n, n_high
   CHARACTER(1024)                              :: FileName
   INTEGER(IntKi)                               :: Un                   ! unit number of opened file   
   
   
   errStat = ErrID_None
   errMsg  = ""
   n = nint(t / p%dt)
   call ComputeLocals(n, u, p, y, m, errStat2, errMsg2)
      call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat2 >= AbortErrLev) then 
            return
      end if
   call LowResGridCalcOutput(n, u, p, y, m, errStat2, errMsg2)
   !call LowResGridCalcOutputOMP(n, u, p, y, m%Vdist_low, m%Vamb_low, m, errStat2, errMsg2)
   
      call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat2 >= AbortErrLev) then 
            return
      end if
   
      ! starting index for the high-res files
   n_high =  n*p%n_high_low
   call HighResGridCalcOutput(n_high, u, p, y, m, errStat2, errMsg2)
      call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat2 >= AbortErrLev) then 
            return
      end if

   if (mod(n,p%WrDisSkp1) == 0) then
      if ( p%WrDisWind  ) then
         call WriteDisWindFiles( n, p%WrDisSkp1, p, y, m, ErrStat2, ErrMsg2 )
      end if
   
         ! XY plane slices
      do k = 1,p%NOutDisWindXY
      
         call ExtractSlice( XYSlice, p%OutDisWindZ(k), p%Z0_low, p%nZ_low, p%nX_low, p%nY_low, p%dZ_low, m%Vdist_low, m%outVizXYPlane(:,:,:,1))        
            ! Create the output vtk file with naming <WindFilePath>/Low/DisXY<k>.t<n/p%WrDisSkp1>.vtk
         FileName = trim(p%OutFileRoot)//".Low.DisXY"//trim(num2lstr(k))//".t"//trim(num2lstr(n/p%WrDisSkp1))//".vtk"
         call WrVTK_SP_header( FileName, "Low resolution, disturbed wind of XY Slice at time = "//trim(num2lstr(t))//" seconds.", Un, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
         call WrVTK_SP_vectors3D( Un, "DisXY", (/p%nX_low,p%nY_low,1_IntKi/), (/p%X0_low,p%Y0_low,p%OutDisWindZ(k)/), (/p%dX_low,p%dY_low,p%dZ_low/), m%outVizXYPlane, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end do
      
   
         ! YZ plane slices
      do k = 1,p%NOutDisWindYZ
         call ExtractSlice( YZSlice, p%OutDisWindX(k), p%X0_low, p%nX_low, p%nY_low, p%nZ_low, p%dX_low, m%Vdist_low, m%outVizYZPlane(:,:,:,1))        
            ! Create the output vtk file with naming <WindFilePath>/Low/DisYZ<k>.t<n/p%WrDisSkp1>.vtk
         FileName = trim(p%OutFileRoot)//".Low.DisYZ"//trim(num2lstr(k))//".t"//trim(num2lstr(n/p%WrDisSkp1))//".vtk"
         call WrVTK_SP_header( FileName, "Low resolution, disturbed wind of YZ Slice at time = "//trim(num2lstr(t))//" seconds.", Un, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
         call WrVTK_SP_vectors3D( Un, "DisYZ", (/1,p%nY_low,p%nZ_low/), (/p%OutDisWindX(k),p%Y0_low,p%Z0_low/), (/p%dX_low,p%dY_low,p%dZ_low/), m%outVizYZPlane, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end do
      
         ! XZ plane slices
      do k = 1,p%NOutDisWindXZ
         call ExtractSlice( XZSlice, p%OutDisWindY(k), p%Y0_low, p%nY_low, p%nX_low, p%nZ_low, p%dY_low, m%Vdist_low, m%outVizXZPlane(:,:,:,1))        
            ! Create the output vtk file with naming <WindFilePath>/Low/DisXZ<k>.t<n/p%WrDisSkp1>.vtk
         FileName = trim(p%OutFileRoot)//".Low.DisXZ"//trim(num2lstr(k))//".t"//trim(num2lstr(n/p%WrDisSkp1))//".vtk"
         call WrVTK_SP_header( FileName, "Low resolution, disturbed wind of XZ Slice at time = "//trim(num2lstr(t))//" seconds.", Un, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
         call WrVTK_SP_vectors3D( Un, "DisXZ", (/p%nX_low,1,p%nZ_low/), (/p%X0_low,p%OutDisWindY(k),p%Z0_low/), (/p%dX_low,p%dY_low,p%dZ_low/), m%outVizXZPlane, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end do
   end if

end subroutine AWAE_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
subroutine AWAE_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, errStat, errMsg )
!..................................................................................................................................

   real(DbKi),                     intent(in   )   :: Time        !< Current simulation time in seconds
   type(AWAE_InputType),           intent(in   )   :: u           !< Inputs at Time
   type(AWAE_ParameterType),       intent(in   )   :: p           !< Parameters
   type(AWAE_ContinuousStateType), intent(in   )   :: x           !< Continuous states at Time
   type(AWAE_DiscreteStateType),   intent(in   )   :: xd          !< Discrete states at Time
   type(AWAE_ConstraintStateType), intent(in   )   :: z           !< Constraint states at Time (possibly a guess)
   type(AWAE_OtherStateType),      intent(in   )   :: OtherState  !< Other states at Time
   type(AWAE_MiscVarType),         intent(inout)   :: m           !< Misc/optimization variables
   type(AWAE_ConstraintStateType), intent(inout)   :: Z_residual  !< Residual of the constraint state equations using
                                                                !!     the input values described above
   integer(IntKi),                 intent(  out)   :: errStat     !< Error status of the operation
   character(*),                   intent(  out)   :: errMsg      !< Error message if errStat /= ErrID_None
   
      ! Local variables   
   integer, parameter                            :: indx = 1  
   integer(intKi)                                :: ErrStat2
   character(ErrMsgLen)                          :: ErrMsg2
   character(*), parameter                       :: RoutineName = 'AWAE_CalcConstrStateResidual'
   
   
   
   errStat = ErrID_None
   errMsg  = ""
  
end subroutine AWAE_CalcConstrStateResidual


   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine validates the inputs from the Wind_AmbientAndArray input files.
subroutine ValidateInitInputData( InputFileData, errStat, errMsg )
!..................................................................................................................................
      
      ! Passed variables:
   type(AWAE_InputFileType), intent(in)     :: InputFileData                     !< All the data in the Wind_AmbientAndArray input file
   integer(IntKi),           intent(out)    :: errStat                           !< Error status
   character(*),             intent(out)    :: errMsg                            !< Error message

   
      ! local variables
   integer(IntKi)                           :: k                                 ! Blade number
   integer(IntKi)                           :: j                                 ! node number
   character(*), parameter                  :: RoutineName = 'ValidateInitInputData'
   
   errStat = ErrID_None
   errMsg  = ""
   
   
   if (len_trim(InputFileData%WindFilePath) == 0) call SetErrStat ( ErrID_Fatal, 'WindFilePath must contain at least one character.', errStat, errMsg, RoutineName )   
   if (  InputFileData%NumTurbines <   1  )  call SetErrStat ( ErrID_Fatal, 'Number of turbines must be greater than zero.', errStat, errMsg, RoutineName )
   if (  InputFileData%NumPlanes   <   2  )  call SetErrStat ( ErrID_Fatal, 'Number of wake planes must be greater than one.', errStat, errMsg, RoutineName )
   if (  InputFileData%NumRadii    <   2  )  call SetErrStat ( ErrID_Fatal, 'Number of radii in the radial finite-difference grid must be greater than one.', errStat, errMsg, RoutineName )
   if (  InputFileData%dr          <=  0.0)  call SetErrStat ( ErrID_Fatal, 'dr must be greater than zero.', errStat, errMsg, RoutineName ) 
   if (.not. ((InputFileData%Mod_Meander == 1) .or. (InputFileData%Mod_Meander == 2) .or. (InputFileData%Mod_Meander == 3)) ) call SetErrStat ( ErrID_Fatal, 'Mod_Meander must be equal to 1, 2, or 3.', errStat, errMsg, RoutineName ) 
   if (  InputFileData%C_Meander   <   1.0_ReKi ) call SetErrStat ( ErrID_Fatal, 'C_Meander must not be less than 1.', errStat, errMsg, RoutineName ) 
   
end subroutine ValidateInitInputData



!=======================================================================
! Unit Tests
!=======================================================================

subroutine AWAE_TEST_Init_BadData(errStat, errMsg)

   integer(IntKi),           intent(out)    :: errStat                           !< Error status
   character(*),             intent(out)    :: errMsg                            !< Error message


   type(AWAE_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(AWAE_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(AWAE_ParameterType)       :: p             !< Parameters
   type(AWAE_ContinuousStateType) :: x             !< Initial continuous states
   type(AWAE_DiscreteStateType)   :: xd            !< Initial discrete states
   type(AWAE_ConstraintStateType) :: z             !< Initial guess of the constraint states
   type(AWAE_OtherStateType)      :: OtherState    !< Initial other states
   type(AWAE_OutputType)          :: y             !< Initial system outputs (outputs are not calculated;
                                                 !!   only the output mesh is initialized)
   type(AWAE_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                   :: interval      !< Coupling interval in seconds: the rate that
   
   type(AWAE_InitOutputType)      :: initOut                         !< Input data for initialization routine
   
                                                                        
   

   
      ! Set up the initialization inputs
   
    
   interval               = 0.0_DbKi
   InitInp%InputFileData%WindFilePath   = '' 
   InitInp%InputFileData%NumTurbines    = 0
   InitInp%InputFileData%NumPlanes      = 0
   InitInp%InputFileData%NumRadii       = 0
   InitInp%InputFileData%dr             = 0.0_ReKi
   InitInp%InputFileData%Mod_Meander    = 0
   InitInp%InputFileData%C_Meander      = 0.0_ReKi

    
   call AWAE_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, errStat, errMsg )
   
   return
   
end subroutine AWAE_TEST_Init_BadData

subroutine AWAE_TEST_SetGoodInitInpData(interval, InitInp)
   real(DbKi)            , intent(out)       :: interval
   type(AWAE_InitInputType), intent(out)       :: InitInp       !< Input data for initialization routine

   ! Based on NREL 5MW
   interval               = 2.0_DbKi
   InitInp%InputFileData%WindFilePath   = 'C:\Dev\OpenFAST-farm\OpenFAST-test\fast-farm\steady' 
   InitInp%InputFileData%WindFilePath   = 'Y:\Wind\Public\Projects\Projects F\FAST.Farm\AmbWind\04'
   InitInp%InputFileData%NumTurbines    = 1
   InitInp%InputFileData%NumPlanes      = 140
   InitInp%InputFileData%NumRadii       = 40
   InitInp%InputFileData%dr             = 5.0_ReKi
   InitInp%n_high_low                   = 6
   InitInp%InputFileData%dt             = 2.0_DbKi
   InitInp%NumDT                        = 1
   InitInp%InputFileData%NOutDisWindXY  = 0
   InitInp%InputFileData%NOutDisWindYZ  = 0
   InitInp%InputFileData%NOutDisWindXZ  = 0
   InitInp%InputFileData%WrDisWind      = .false.
   InitInp%InputFileData%WrDisDT        = 0.0
   InitInp%InputFileData%OutDisWindY    = 0
   InitInp%InputFileData%OutDisWindZ    = 0
   InitInp%InputFileData%OutDisWindX    = 0
   InitInp%InputFileData%Mod_Meander    = 1
   InitInp%InputFileData%C_Meander      = 2.0_ReKi


end subroutine AWAE_TEST_SetGoodInitInpData


subroutine AWAE_TEST_Init_GoodData(errStat, errMsg)

   integer(IntKi),           intent(out)    :: errStat                           !< Error status
   character(*),             intent(out)    :: errMsg                            !< Error message


   type(AWAE_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(AWAE_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(AWAE_ParameterType)       :: p             !< Parameters
   type(AWAE_ContinuousStateType) :: x             !< Initial continuous states
   type(AWAE_DiscreteStateType)   :: xd            !< Initial discrete states
   type(AWAE_ConstraintStateType) :: z             !< Initial guess of the constraint states
   type(AWAE_OtherStateType)      :: OtherState    !< Initial other states
   type(AWAE_OutputType)          :: y             !< Initial system outputs (outputs are not calculated;
                                                 !!   only the output mesh is initialized)
   type(AWAE_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                   :: interval      !< Coupling interval in seconds: the rate that
   
   type(AWAE_InitOutputType)      :: initOut                         !< Input data for initialization routine
   
                                                                        
   

   
      ! Set up the initialization inputs
   call AWAE_TEST_SetGoodInitInpData(interval, InitInp)

   call AWAE_Init( InitInp, u, p, x, xd, z, OtherState, y, m, interval, InitOut, errStat, errMsg )
   
   return
   
end subroutine AWAE_TEST_Init_GoodData

#ifdef PARALLEL_CODE
subroutine AWAE_TEST_LowResGridCalcs(errStat, errMsg)

   integer(IntKi),           intent(out)    :: errStat                           !< Error status
   character(*),             intent(out)    :: errMsg                            !< Error message

   real(DbKi)                               :: t1, t2, t3
   type(AWAE_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(AWAE_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(AWAE_ParameterType)       :: p             !< Parameters
   type(AWAE_ContinuousStateType) :: x             !< Initial continuous states
   type(AWAE_DiscreteStateType)   :: xd            !< Initial discrete states
   type(AWAE_ConstraintStateType) :: z             !< Initial guess of the constraint states
   type(AWAE_OtherStateType)      :: OtherState    !< Initial other states
   type(AWAE_OutputType)          :: y             !< Initial system outputs (outputs are not calculated;
                                                 !!   only the output mesh is initialized)
   type(AWAE_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                   :: interval      !< Coupling interval in seconds: the rate that
   real(DbKi)                     :: t
   type(AWAE_InitOutputType)      :: initOut                         !< Input data for initialization routine
   integer(IntKi)                 :: n,ix,iy,iz
   real(SiKi),allocatable                     :: Vamb_lowOMP(:,:,:,:)
   real(SiKi),allocatable                     :: Vdist_lowOMP(:,:,:,:)
      ! Set up the initialization inputs
   call AWAE_TEST_SetGoodInitInpData(interval, InitInp)

   call AWAE_Init( InitInp, u, p, x, xd, z, OtherState, y, m, interval, InitOut, errStat, errMsg )
   
   
   allocate ( Vamb_lowOMP   ( 3, 0:p%nX_low-1 , 0:p%nY_low-1 , 0:p%nZ_low-1 )                  , STAT=errStat ) 
      if (errStat /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vamb_low.', errStat, errMsg, 'AWAE_TEST_LowResGridCalcs' )   
   allocate ( Vdist_lowOMP  ( 3, 0:p%nX_low-1 , 0:p%nY_low-1 , 0:p%nZ_low-1 )                  , STAT=errStat ) 
      if (errStat /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vdist_low.', errStat, errMsg, 'AWAE_TEST_LowResGridCalcs' )  
   if (errStat > ErrID_None) return
   
   n = 0 ! starting time
   t = n*p%dt
      ! Update states is really reading the wind data for n+1 time step
   call AWAE_UpdateStates( t, n, u, p, x, xd, z, OtherState, m, errStat, errMsg )
   
   n=n+1 ! increment time
   
   t1 = omp_get_wtime()
   call LowResGridCalcOutput(n, u, p, y, m, errStat, errMsg)
      if (errStat >= AbortErrLev) then 
            return
      end if   
   t2 = omp_get_wtime()
   call LowResGridCalcOutputOMP(n, u, p, y, Vdist_lowOMP, m%Vamb_low, m, errStat, errMsg)
      if (errStat >= AbortErrLev) then 
            return
      end if   
   t3 = omp_get_wtime()
   write(*,*) 'Serial Low Res Grid took '//trim(num2lstr(t2-t1))//' seconds.'
   write(*,*) 'Parallel Low Res Grid took '//trim(num2lstr(t3-t2))//' seconds.'
   
   do iz = 0,p%nZ_low-1
      do iy = 0,p%nY_low-1
         do ix = 0,p%nX_low-1
            if ( .not. (EqualRealNos(Vdist_lowOMP(1,ix,iy,iz),m%Vdist_low(1,ix,iy,iz)) .and. EqualRealNos(Vdist_lowOMP(2,ix,iy,iz),m%Vdist_low(2,ix,iy,iz)) .and. EqualRealNos(Vdist_lowOMP(3,ix,iy,iz),m%Vdist_low(3,ix,iy,iz)) ) ) then
               WRITE(*,*) 'Values do not match between serial and parallel results'
            end if        
         end do
      end do
   end do      
   
end subroutine AWAE_TEST_LowResGridCalcs 

subroutine AWAE_TEST_ExtractSlice(errStat, errMsg)

   integer(IntKi),           intent(out)    :: errStat                           !< Error status
   character(*),             intent(out)    :: errMsg                            !< Error message
   
   integer(IntKi)                           :: nx, ny, nz, ix, iy, iz
   real(SiKi), allocatable                  :: V(:,:,:,:),sliceXY(:,:,:),sliceXY_OMP(:,:,:),sliceYZ(:,:,:),sliceXZ(:,:,:)
   real(ReKi)                               :: dx,dy,dz,x0,y0,z0,x,y,z
   real(DbKi)                               :: t1, t2, t3
   
      ! Create a volume of random number data
   nx = 500
   ny = 500
   nz = 100
   dx = 10.0
   dy = 10.0
   dz = 10.0
   x0 = 500.0
   y0 = 300.0
   z0 = 10.0
   !allocate( p%r(0:p%NumRadii-1),stat=errStat2)
   allocate(V(3,0:nx-1,0:ny-1,0:nz-1))
   allocate(sliceXY(3,0:nx-1,0:ny-1))
   allocate(sliceXY_OMP(3,0:nx-1,0:ny-1))
   do iz = 0,nz-1
      do iy = 0,ny-1
         do ix = 0,nx-1
            CALL RANDOM_NUMBER(V(1,ix,iy,iz) )
            CALL RANDOM_NUMBER(V(2,ix,iy,iz) )
            CALL RANDOM_NUMBER(V(3,ix,iy,iz) )            
         end do
      end do      
   end do
   CALL RANDOM_NUMBER(z)
   z = z*nz*dz + z0
   t1 = omp_get_wtime()
   call ExtractSlice( XYSlice, z, z0, nz, nx, ny, dz,  V, sliceXY)
   t2 = omp_get_wtime()
   call ExtractSliceOMP( XYSlice, z, z0, nz, nx, ny, dz,  V, sliceXY_OMP)
   t3 = omp_get_wtime()
   write(*,*) 'Serial Extract slice took '//trim(num2lstr(t2-t1))//' seconds.'
   write(*,*) 'Parallel Extract slice took '//trim(num2lstr(t3-t2))//' seconds.'
      do iy = 0,ny-1
         do ix = 0,nx-1
            if ( .not. (EqualRealNos(sliceXY(1,ix,iy),sliceXY_OMP(1,ix,iy)) .and. EqualRealNos(sliceXY(2,ix,iy),sliceXY_OMP(2,ix,iy)) .and. EqualRealNos(sliceXY(3,ix,iy),sliceXY_OMP(3,ix,iy)) ) ) then
               WRITE(*,*) 'Values do not match between serial and parallel results'
            end if        
         end do
      end do      
end subroutine AWAE_TEST_ExtractSlice
#endif

subroutine AWAE_TEST_CalcOutput(errStat, errMsg)

   integer(IntKi),           intent(out)    :: errStat                           !< Error status
   character(*),             intent(out)    :: errMsg                            !< Error message


   type(AWAE_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(AWAE_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(AWAE_ParameterType)       :: p             !< Parameters
   type(AWAE_ContinuousStateType) :: x             !< Initial continuous states
   type(AWAE_DiscreteStateType)   :: xd            !< Initial discrete states
   type(AWAE_ConstraintStateType) :: z             !< Initial guess of the constraint states
   type(AWAE_OtherStateType)      :: OtherState    !< Initial other states
   type(AWAE_OutputType)          :: y             !< Initial system outputs (outputs are not calculated;
                                                 !!   only the output mesh is initialized)
   type(AWAE_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                   :: interval      !< Coupling interval in seconds: the rate that
   
   type(AWAE_InitOutputType)      :: initOut                         !< Input data for initialization routine
   
   integer(IntKi)  :: nt, nr, np
   real(DbKi) :: t
   
   ! This example creates turbine 1 at the global coordinate [0,0,0]
   ! The data is hardcoded in: AWAE_IO_InitGridInfo() as follows:
   ! X0_low = -750.0_ReKi
   ! Y0_low = -500.0_ReKi
   ! Z0_low = 0.0_ReKi
   ! dX_low = 10.0_ReKi
   ! dY_low = 10.0_ReKi
   ! dZ_low = 10.0_ReKi
   !    ! Parse a low res wind input file to gather the grid information
   ! p%nX_Low           = 151    
   ! p%nY_low           = 101    
   ! p%nZ_low           = 51     
   !    ! Grid runs from (X0_low, Y0_low, Z0_low) to (X0_low + (p%nX_Low-1)*dX_low, Y0_low+ (p%nY_Low-1)*dY_low, Z0_low+ (p%nZ_Low-1)*dZ_low)
   !    ! (0,0,0) to (180,180,180) 
   !    ! Parse a high res wind input file to gather the grid information
   ! p%nX_high          = 16 
   ! p%nY_high          = 16 
   ! p%nZ_high          = 16 
   ! The low resolution grid extends from [-750,-500,0] to [750,500,500]
   ! The first turbine's grid is located at [
   
      ! Based on NREL 5MW
    interval               = 1.0_DbKi
    InitInp%InputFileData%WindFilePath   = 'C:\Dev\NWTC Github\FAST.Farm\data' 
    InitInp%InputFileData%NumTurbines    = 3
    InitInp%InputFileData%NumPlanes      = 500
    InitInp%InputFileData%NumRadii       = 40
    InitInp%InputFileData%dr             = 5.0_ReKi
   
      ! Initialize the module
   call AWAE_Init( InitInp, u, p, x, xd, z, OtherState, y, m, interval, InitOut, errStat, errMsg )
   if (errStat > ErrID_None) then 
      return
   end if
   
   
      ! Set up the inputs
   do nt = 1,p%NumTurbines
      do np = 0,p%NumPlanes-1
         do nr = 0,p%NumRadii-1         
            u%Vx_wake(nr,np,nt) = -1.0_ReKi      
            u%Vr_wake(nr,np,nt) =  0.1_ReKi    
         end do
      end do
   end do
   
  
   u%xhat_plane(1,:,:) = 1.0_ReKi  
   u%xhat_plane(2,:,:) = 0.0_ReKi  
   u%xhat_plane(3,:,:) = 0.0_ReKi 
   
   do nt = 1,p%NumTurbines
      do np = 0,p%NumPlanes-1
         u%p_plane(1,np,nt)    = 0.0_ReKi + 8.0*np*interval + 250.0_ReKi*(nt-1)
         u%p_Plane(2,np,nt)    = 0.0_ReKi
         u%p_Plane(3,np,nt)    = 90.0_ReKi
         u%D_wake(np,nt)       = 126.0_ReKi
      end do
   end do

   t = 0.0_DbKi
   
   call AWAE_CalcOutput(t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
   if (errStat > ErrID_None) then 
      return
   end if
  ! call AWAE_UpdateStates(t, 0, u, p, x, xd, z, OtherState, m, errStat, errMsg )
   
   !t = t + interval
   !call AWAE_CalcOutput(t, u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
   !
   !   ! Verify that xd and y are the same
   !
   !if (errStat == ErrID_None) then
   !   call AWAE_UpdateStates(0.0_DbKi, 1, u, p, x, xd, z, OtherState, m, errStat, errMsg )
   !end if
   
   return


end subroutine AWAE_TEST_CalcOutput
end module AWAE
