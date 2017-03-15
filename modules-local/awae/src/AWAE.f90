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
  
   contains  
   


!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine LowResGridCalcOutput(t, u, p, y, m, errStat, errMsg)
   real(DbKi),                     intent(in   )  :: t           !< Current simulation time in seconds
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
   real(ReKi)          :: Vx_wake_tmp
   real(ReKi)          :: Vr_wake_tmp(3)
   real(ReKi)          :: Vr_term(3)
   real(ReKi)          :: Vx_term
   real(ReKi)          :: Vsum_low(3)
   real(ReKi)          :: Vave_amb_low_norm
   integer(IntKi)      :: ILo
   character(*), parameter   :: RoutineName = 'LowResGridCalcOutput'   
   errStat = ErrID_None
   errMsg  = ""

   
      ! read from file the ambient flow for the current time step
   call ReadLowResWindFile(t, p, m%Vamb_Low, errStat, errMsg)
      if ( errStat > AbortErrLev ) then
         return
      end if
   
            
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
               
               do np = 0, p%NumPlanes-2
                  
                     ! Reset interpolation counter
                  ILo = 0
                  
                     ! Construct the endcaps of the current wake plane volume
                  x_start_plane = x_end_plane
                  x_end_plane = dot_product(u%xhat_plane(:,np+1,nt), (p%Grid_Low(:,nXYZ_low) - u%p_plane(:,np+1,nt)) )
                  
                     ! test if the point is within the endcaps of the wake volume
                  if ( ( x_start_plane >= 0.0_ReKi ) .and. ( x_end_plane < 0.0_ReKi ) ) then
                     r_vec_plane = p%Grid_Low(:,nXYZ_low) - u%p_plane(:,np,nt) - x_start_plane*u%xhat_plane(:,np,nt)
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
                        m%Vx_wake(n_wake) = InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np,nt), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                        m%Vr_wake(n_wake) = InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np,nt), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                        
                        
                        m%xhat_plane(:,n_wake) = u%xhat_plane(:,np,nt)
                        xhatBar_plane = xhatBar_plane + abs(m%Vx_wake(n_wake))*m%xhat_plane(:,n_wake)
      
                     end if  ! if the point is within radial finite-difference grid
                     
                        ! test if the point is within the radius of the wake volume cylinder
                     if ( r_tmp_plane <= u%D_wake(np,nt) ) then
                        m%N_wind(np,nt) = m%N_wind(np,nt) + 1
                        
                        ! TODO: Verify that m%N_wind(np,nt) <= MAX_1ST_DIM_SIZE of nx_wind
                        if ( m%N_wind(np,nt) > p%n_wind_max ) then
                           call SetErrStat( ErrID_Fatal, 'The wake plane volume (plane='//trim(num2lstr(np))//',turbine='//trim(num2lstr(nt))//') contains more points than the maximum predicted points: 30*pi*DT(2*r*[Nr-1])**2/(dx*dy*dz)', errStat, errMsg, RoutineName )
                           return
                        end if
                        
                        m%nx_wind(m%N_wind(np,nt),np,nt) = nx_low
                        m%ny_wind(m%N_wind(np,nt),np,nt) = ny_low
                        m%nz_wind(m%N_wind(np,nt),np,nt) = nz_low   
                        
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
               m%Vdist_low(:,nx_low,ny_low,nz_low) = m%Vdist_low(:,nx_low,ny_low,nz_low) + Vr_wake_tmp - xhatBar_plane*sqrt(Vx_wake_tmp)
            end if  ! (n_wake > 0)
            
         end do ! do nx_low=0, p%nX_low-1
      end do    ! do ny_low=0, p%nY_low-1
   end do       ! do nz_low=0, p%nZ_low-1
   
   do nt = 1,p%NumTurbines
      if ( m%N_wind(0,nt) > 0 ) then
         ! TODO: This is modified in Rev 7 and the current code is for Rev 6. 13/Feb/2017
         Vsum_low = 0.0_ReKi
      
         do nw=1,m%N_wind(0,nt)   
            Vsum_low = Vsum_low + m%Vamb_Low(:, m%nx_wind(nw,0,nt), m%ny_wind(nw,0,nt), m%nz_wind(nw,0,nt))
         end do
         Vsum_low       = Vsum_low / m%N_wind(0,nt)  ! if N_wind gets large ( ~= 100,000 ) then this may not give enough precision in Vave_amb_low
         Vave_amb_low_norm  = TwoNorm(Vsum_low)
         if ( EqualRealNos(Vave_amb_low_norm,0.0_ReKi) ) then    
            call SetErrStat( ErrID_Fatal, 'The magnitude of the spatial-averaged ambient wind speed in the low-resolution domain associated with the wake volume at the rotor disk for turbine '//trim(num2lstr(nt))//' is zero.', errStat, errMsg, RoutineName )
            return     
         end if
      
         y%Vx_wind_disk(nt) = dot_product( u%xhat_plane(:,0,nt), Vsum_low )
         y%TI_amb(nt)       = 0.0_ReKi
         do nw=1,m%N_wind(np,nt)
            y%TI_amb(nt) = y%TI_amb(nt) + TwoNorm( m%Vamb_Low(:, m%nx_wind(nw,0,nt), m%ny_wind(nw,0,nt), m%nz_wind(nw,0,nt)) - Vsum_low )**2
         end do
         y%TI_amb(nt) = sqrt(y%TI_amb(nt)/(3.0*m%N_wind(0,nt)))/Vave_amb_low_norm
      else
         y%Vx_wind_disk(nt) = 0.0_ReKi
         y%TI_amb(nt)       = 0.0_ReKi
      end if
      
      
      do np = 0, p%NumPlanes-2
         if ( (u%D_wake(np,nt) > 0.0_ReKi) .and.  (m%N_wind(np,nt) < p%n_wind_min) ) then
            call SetErrStat( ErrID_Fatal, 'The number of points in the wake volume '//trim(num2lstr(np))//' for turbine '//trim(num2lstr(nt))//' is less than the minimum threshold, '//trim(num2lstr(p%n_wind_min))//'.', errStat, errMsg, RoutineName )
            return     
         else if ( m%N_wind(np,nt) > 0  ) then            
            Vsum_low = 0.0_ReKi
            do nw=1,m%N_wind(np,nt)   
               Vsum_low = Vsum_low + m%Vdist_low( :, m%nx_wind(nw,np,nt),m%ny_wind(nw,np,nt),m%nz_wind(nw,np,nt) )
            end do
            y%V_plane(:,np,nt) = Vsum_low /  m%N_wind(np,nt)
         else
            y%V_plane(:,np,nt) = 0.0_ReKi
         end if
         
      end do
      
   end do
   
end subroutine LowResGridCalcOutput


!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine 
!!
subroutine HighResGridCalcOutput(t, u, p, y, m, errStat, errMsg)
   real(DbKi),                     intent(in   )  :: t           !< Current simulation time in seconds
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
   integer(IntKi)      :: ILo
   character(*), parameter   :: RoutineName = 'HighResGridCalcOutput'
   errStat = ErrID_None
   errMsg  = ""

   
   
            
   
   
      ! Loop over the entire grid of low resolution ambient wind data to compute:
      !    1) the disturbed flow at each point and 2) the averaged disturbed velocity of each wake plane
   

   do nt = 1,p%NumTurbines
      nXYZ_high = 0
      do n_hl=0, p%n_high_low-1
            ! read from file the ambient flow for the current time step
         call ReadHighResWindFile(nt, n_hl, t, p, m%Vamb_high, errStat, errMsg)
            if ( errStat > AbortErrLev ) then
               return
            end if
            
            ! set the disturbed flow equal to the ambient flow for this time step
         y%Vdist_high(:,:,:,:,n_hl,nt) = m%Vamb_high(:,:,:,:)
      end do
      
      do nz_high=0, p%nZ_high-1 
         do ny_high=0, p%nY_high-1
            do nx_high=0, p%nX_high-1
               
               
               nXYZ_high = nXYZ_high + 1
               n_wake = 0
               xhatBar_plane = 0.0_ReKi
            
               do nt2 = 1,p%NumTurbines
                  if (nt /= nt2) then  
                     
                     x_end_plane = dot_product(u%xhat_plane(:,0,nt2), (p%Grid_high(:,nXYZ_high,nt2) - u%p_plane(:,0,nt2)) )
               
                     do np = 0, p%NumPlanes-2
                  
                           ! Reset interpolation counter
                        ILo = 0
                  
                           ! Construct the endcaps of the current wake plane volume
                        x_start_plane = x_end_plane
                        x_end_plane = dot_product(u%xhat_plane(:,np+1,nt2), (p%Grid_high(:,nXYZ_high,nt2) - u%p_plane(:,np+1,nt2)) )
                  
                           ! test if the point is within the endcaps of the wake volume
                        if ( ( x_start_plane >= 0.0_ReKi ) .and. ( x_end_plane < 0.0_ReKi ) ) then
                           r_vec_plane = p%Grid_high(:,nXYZ_high,nt2) - u%p_plane(:,np,nt2) - x_start_plane*u%xhat_plane(:,np,nt2)
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
                              m%Vx_wake(n_wake) = InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np,nt2), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                              m%Vr_wake(n_wake) = InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np,nt2), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                                              
                              m%xhat_plane(:,n_wake) = u%xhat_plane(:,np,nt2)
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
                  do n_hl=0, p%n_high_low-1
                     y%Vdist_high(:,nx_high,ny_high,nz_high,n_hl,nt) = y%Vdist_high(:,nx_high,ny_high,nz_high,n_hl,nt) + Vr_wake_tmp - xhatBar_plane*sqrt(Vx_wake_tmp)
                  end do
                     
               end if  ! (n_wake > 0)
            
            end do ! nx_high=0, p%nX_high-1
         end do    ! ny_high=0, p%nY_high-1
      end do       ! nz_high=0, p%nZ_high-1
   end do          ! nt = 1,p%NumTurbines
   
   
end subroutine HighResGridCalcOutput

subroutine HighResGridCalcOutput2(t, u, p, y, m, errStat, errMsg)
   real(DbKi),                     intent(in   )  :: t           !< Current simulation time in seconds
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
   integer(IntKi)      :: ILo
   character(*), parameter   :: RoutineName = 'HighResGridCalcOutput'
   errStat = ErrID_None
   errMsg  = ""

   
   
            
   
   
      ! Loop over the entire grid of low resolution ambient wind data to compute:
      !    1) the disturbed flow at each point and 2) the averaged disturbed velocity of each wake plane
   
   
   
   do nt = 1,p%NumTurbines
     ! nXYZ_high = 0
      do n_hl=0, p%n_high_low-1
            ! read from file the ambient flow for the current time step
         call ReadHighResWindFile(nt, n_hl, t, p, m%Vamb_high, errStat, errMsg)
            if ( errStat > AbortErrLev ) then
               return
            end if
            
            ! set the disturbed flow equal to the ambient flow for this time step
         y%Vdist_high(:,:,:,:,n_hl,nt) = m%Vamb_high(:,:,:,:)
      end do
      
      !$omp parallel &
      !$omp& private(nz_high,ny_high,nx_high,nXYZ_high,n_wake,xhatBar_plane,nt2), &
      !$omp& private(x_end_plane,np,ILo,x_start_plane,r_vec_plane,r_tmp_plane,tmp_xhatBar_plane), &
      !$omp& private(Vx_wake_tmp,Vr_wake_tmp,Vr_term,Vx_term,n_hl ), & 
      !$omp& shared(nt,p%nX_high,p%nY_high,p%nZ_high,u%xhat_plane,p%Grid_high,u%p_plane,m%rhat_plane), &
      !$omp& shared(p%NumRadii,m%Vx_wake,m%Vr_wake,p%r,u%Vx_wake,u%Vr_wake,m%xhat_plane,y%Vdist_high), &
      !$omp& collapse(3)
      do nz_high=0, p%nZ_high-1 
         do ny_high=0, p%nY_high-1
            do nx_high=0, p%nX_high-1
         
               nXYZ_high =  nz_high*p%nX_high*p%nZ_high + ny_high*p%nX_high + nx_high !  nXYZ_high + 1
               n_wake = 0
               xhatBar_plane = 0.0_ReKi
            
               do nt2 = 1,p%NumTurbines
                  if (nt /= nt2) then  
                     
                     x_end_plane = dot_product(u%xhat_plane(:,0,nt2), (p%Grid_high(:,nXYZ_high,nt2) - u%p_plane(:,0,nt2)) )
               
                     do np = 0, p%NumPlanes-2
                  
                           ! Reset interpolation counter
                        ILo = 0
                  
                           ! Construct the endcaps of the current wake plane volume
                        x_start_plane = x_end_plane
                        x_end_plane = dot_product(u%xhat_plane(:,np+1,nt2), (p%Grid_high(:,nXYZ_high,nt2) - u%p_plane(:,np+1,nt2)) )
                  
                           ! test if the point is within the endcaps of the wake volume
                        if ( ( x_start_plane >= 0.0_ReKi ) .and. ( x_end_plane < 0.0_ReKi ) ) then
                           r_vec_plane = p%Grid_high(:,nXYZ_high,nt2) - u%p_plane(:,np,nt2) - x_start_plane*u%xhat_plane(:,np,nt2)
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
                              m%Vx_wake(n_wake) = InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np,nt2), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                              m%Vr_wake(n_wake) = InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np,nt2), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                                              
                              m%xhat_plane(:,n_wake) = u%xhat_plane(:,np,nt2)
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
                  do n_hl=0, p%n_high_low-1
                     y%Vdist_high(:,nx_high,ny_high,nz_high,n_hl,nt) = y%Vdist_high(:,nx_high,ny_high,nz_high,n_hl,nt) + Vr_wake_tmp - xhatBar_plane*sqrt(Vx_wake_tmp)
                  end do
                     
               end if  ! (n_wake > 0)
            
            end do ! nx_high=0, p%nX_high-1
         end do    ! ny_high=0, p%nY_high-1
   end do       ! nz_high=0, p%nZ_high-1
   !$omp end parallel
   end do          ! nt = 1,p%NumTurbines

   
end subroutine HighResGridCalcOutput2

!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine sets the initialization output data structure, which contains data to be returned to the calling program (e.g.,
!! FAST or Wind_AmbientAndArray_Driver)   
subroutine AWAE_SetInitOut(p, InputInp, InitOut, errStat, errMsg)

   type(AWAE_InitOutputType),       intent(  out)  :: InitOut          !< Initialization output data
   type(AWAE_InitInputType),        intent(in   )  :: InputInp         !< Initialization input  data 
   type(AWAE_ParameterType),        intent(in   )  :: p                !< Parameters
   integer(IntKi),                 intent(  out)  :: errStat          !< Error status of the operation
   character(*),                   intent(  out)  :: errMsg           !< Error message if errStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'AWAE_SetInitOut'
   
   
   
   integer(IntKi)                               :: i, j, k, f
   integer(IntKi)                               :: NumCoords

      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""

      ! Set output data
   
   InitOut%Ver = AWAE_Ver
 
end subroutine AWAE_SetInitOut

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
   integer(IntKi)                                :: i             ! loop counter
                                              
   integer(IntKi)                                :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                          :: errMsg2       ! temporary error message 
                                              
   !type(AWAE_InitInputType)                       :: InputFileData ! Data stored in the module's input file
   integer(IntKi)                                :: UnEcho        ! Unit number for the echo file
                                              
   character(*), parameter                       :: RoutineName = 'AWAE_Init'
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   UnEcho  = -1  ! TODO: May not need file echoing in this module

      ! Initialize the NWTC Subroutine Library

   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   call DispNVD( AWAE_Ver )
   
   
   ! TODO: For output reporting in this module we need to have Rootname include the turbine number 
   
   !p%OutFileRoot  = TRIM(InitInp%RootName)//'.WAA'
   
   
         
      
      ! Validate the initialization inputs
   call ValidateInitInputData( InitInp%InputFileData, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName ) 
      if (errStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
      
      !............................................................................................
      ! Define parameters
      !............................................................................................
      
   
      
      ! set the rest of the parameters
   !p%DT               = interval       ! TODO: So we don't actually obtain the DT from the driver, but from the wind data files    
   p%NumPlanes        = InitInp%InputFileData%NumPlanes   
   p%NumRadii         = InitInp%InputFileData%NumRadii    
   p%NumTurbines      = InitInp%InputFileData%NumTurbines 
   p%WindFileRoot     = InitInp%InputFileData%WindFilePath
   p%n_high_low       = InitInp%n_high_low
   
   allocate( p%r(0:p%NumRadii-1),stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for p%r.', errStat, errMsg, RoutineName )
         call Cleanup()
         return
      end if
    
   do i = 0,p%NumRadii-1
      p%r(i)       = InitInp%InputFileData%dr*i     
   end do
   
   
      ! Obtain the precursor grid information by parsing the necessary input files
      ! This will establish certain parameters as well as all of the initialization outputs
   call AWAE_IO_InitGridInfo(p, InitOut, errStat2, errMsg2)
      call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if (errStat2 > AbortErrLev) then
         
         call Cleanup() 
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
   allocate ( y%Vdist_High   (3,0:p%nX_high-1,0:p%nY_high-1,0:p%nZ_high-1,0:p%n_high_low-1, 1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%Vdist_High.', errStat, errMsg, RoutineName )  
   allocate ( y%Vx_wind_disk   (1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%Vx_rel_disk.', errStat, errMsg, RoutineName )  
   allocate ( y%TI_amb   (1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%TI_amb.', errStat, errMsg, RoutineName )  
   if (errStat /= ErrID_None) return
   
      ! This next step is not strictly necessary
   y%V_plane       = 0.0_Reki
   y%Vdist_High    = 0.0_Reki
   y%Vx_wind_disk  = 0.0_Reki
   y%TI_amb        = 0.0_Reki
   
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
   allocate ( m%Vamb_high  ( 3, 0:p%nX_high-1, 0:p%nY_high-1, 0:p%nZ_high-1 ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vamb_High.', errStat, errMsg, RoutineName )  
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
      
      ! TODO: How do we really want to set this first dimension
   allocate ( m%nx_wind    ( MAX_WAKE_VOL_PTS, 0:p%NumPlanes-2, 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%nx_wind.', errStat, errMsg, RoutineName )  
   allocate ( m%ny_wind    ( MAX_WAKE_VOL_PTS, 0:p%NumPlanes-2, 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%ny_wind.', errStat, errMsg, RoutineName )  
   allocate ( m%nz_wind    ( MAX_WAKE_VOL_PTS, 0:p%NumPlanes-2, 1:p%NumTurbines ), STAT=errStat2 ) 
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%nz_wind.', errStat, errMsg, RoutineName )  
   
   if (errStat /= ErrID_None) return
   
   ! TODO: This step isn't really needed
   m%Vamb_low     = 0.0_ReKi
   m%Vdist_low    = 0.0_ReKi
   m%Vamb_High    = 0.0_ReKi 
   m%N_wind       = 0 
   m%xhat_plane   = 0.0_ReKi 
   m%rhat_plane   = 0.0_ReKi
   m%Vx_wake      = 0.0_ReKi   
   m%Vr_wake      = 0.0_ReKi   
   m%nx_wind      = 0.0_ReKi   
   m%ny_wind      = 0.0_ReKi  
   m%nz_wind      = 0.0_ReKi 

   call Cleanup() 
      
contains
   subroutine Cleanup()

      ! TODO: May not need file echoing in this module

      IF ( UnEcho > 0 ) CLOSE( UnEcho )
      
   end subroutine Cleanup

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
   real(ReKi)                                   :: lstar, dx, Vx_wake_min, r_wake, V_planeDT(3), a_interp, norm2_xhat_plane, EddyTermA, EddyTermB  
   real(ReKi)                                   :: dy_HWkDfl(3)
   integer(intKi)                               :: i,j, ILo
   
   errStat = ErrID_None
   errMsg  = ""
 
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
   integer(intKi)                               :: i, j
   integer(intKi)                               :: errStat2
   character(ErrMsgLen)                         :: errMsg2
   character(*), parameter                      :: RoutineName = 'AWAE_CalcOutput'
   real(ReKi)                                   :: correction(3)
   real(ReKi)                                   :: x_plane
   
   errStat = ErrID_None
   errMsg  = ""

   call LowResGridCalcOutput(t, u, p, y, m, errStat2, errMsg2)
      call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat2 > AbortErrLev) then 
            return
      end if
   call HighResGridCalcOutput(t, u, p, y, m, errStat2, errMsg2)
      call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat2 > AbortErrLev) then 
            return
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
    
   call AWAE_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, errStat, errMsg )
   
   return
   
end subroutine AWAE_TEST_Init_BadData

subroutine AWAE_TEST_SetGoodInitInpData(interval, InitInp)
   real(DbKi)            , intent(out)       :: interval
   type(AWAE_InitInputType), intent(out)       :: InitInp       !< Input data for initialization routine

      ! Based on NREL 5MW
    interval               = 1.0_DbKi
    InitInp%InputFileData%WindFilePath   = 'C:\Dev\NWTC Github\FAST.Farm\data' 
    InitInp%InputFileData%NumTurbines    = 1
    InitInp%InputFileData%NumPlanes      = 500
    InitInp%InputFileData%NumRadii       = 40
    InitInp%InputFileData%dr             = 5.0_ReKi

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
