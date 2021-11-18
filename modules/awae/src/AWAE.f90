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
   use InflowWind_Types
   use InflowWind

#ifdef _OPENMP
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


   contains


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
         m%r_e(np,nt) = abs( dot_product( u%xhat_plane(:,np  ,nt), dp ) )
         m%r_s(np,nt) = abs( dot_product( u%xhat_plane(:,np+1,nt), dp ) )

         if (   sinTerm > ( max( m%r_e(np,nt), m%r_s(np,nt) ) / ( 100.0_ReKi*rmax ) ) ) then
            m%parallelFlag(np,nt) = .false.
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

   integer(IntKi)      :: nx, ny, nz, nt, np, nw, nx_low, ny_low, nz_low, nr, npsi, wamb, iwsum !< loop counters
   integer(IntKi)      :: nXYZ_low, n_wake, n_r_polar, n_psi_polar       !< accumulating counters
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
   real(ReKi)          :: p_tmp_plane(3)
   real(ReKi)          :: tmp_vec(3)
   real(ReKi)          :: Vave_amb_low_norm, Vamb_lowpol_tmp(3), Vdist_lowpol_tmp(3), Vamb_low_tmp(3,8)
   real(ReKi)          :: delta, deltad
   real(ReKi)          :: wsum_tmp, w
   real(ReKi)          :: tmp_x,tmp_y,tmp_z !, tm1, tm2
   real(ReKi)          :: xxplane(3), xyplane(3), yyplane(3), yxplane(3), psi_polar, r_polar, p_polar(3)
   real(ReKi)          :: yzplane_Y(3), xyplane_norm
   real(ReKi)          :: xplane_sq, yplane_sq, xysq_Z(3), xzplane_X(3), tmp_yhat_plane(3), tmp_zhat_plane(3)
   real(ReKi), ALLOCATABLE :: tmp_rhat_plane(:,:), tmp_xhat_plane(:,:)
   real(ReKi), ALLOCATABLE :: tmp_Vx_wake(:), tmp_Vr_wake(:)
   integer(IntKi)      :: ILo
   integer(IntKi)      :: maxPln, tmpPln, maxN_wake
   integer(IntKi)      :: i,np1,errStat2
   character(*), parameter   :: RoutineName = 'LowResGridCalcOutput'
   logical             :: within

   errStat = ErrID_None
   errMsg  = ""

   maxPln =  min(n,p%NumPlanes-2)
   tmpPln =  min(p%NumPlanes-1, n+1)
   
   
   
!#ifdef _OPENMP  
!   tm1 =  omp_get_wtime() 
!#endif 

   maxN_wake = p%NumTurbines*( p%NumPlanes-1 )
   ! Temporary variables needed by OpenMP 
   allocate ( tmp_xhat_plane ( 3, 1:maxN_wake ), STAT=errStat2 )
       if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for tmp_xhat_plane.', errStat, errMsg, RoutineName )
   allocate ( tmp_rhat_plane ( 3, 1:maxN_wake ), STAT=errStat2 )
       if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for tmp_rhat_plane.', errStat, errMsg, RoutineName )
   allocate ( tmp_Vx_wake    ( 1:maxN_wake ), STAT=errStat2 )
       if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for tmp_Vx_wake.', errStat, errMsg, RoutineName )
   allocate ( tmp_Vr_wake    ( 1:maxN_wake ), STAT=errStat2 )
       if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for tmp_Vr_wake.', errStat, errMsg, RoutineName )
   if (ErrStat >= AbortErrLev) return


      ! Loop over the entire grid of low resolution ambient wind data to compute:
      !    1) the disturbed flow at each point and 2) the averaged disturbed velocity of each wake plane

   !$OMP PARALLEL DO PRIVATE(nx_low,ny_low,nz_low, nXYZ_low, n_wake, xhatBar_plane, x_end_plane,nt,np,ILo,x_start_plane,delta,deltad,p_tmp_plane,tmp_vec,r_vec_plane,r_tmp_plane,tmp_xhatBar_plane, Vx_wake_tmp,Vr_wake_tmp,nw,Vr_term,Vx_term,tmp_x,tmp_y,tmp_z,tmp_xhat_plane,tmp_rhat_plane,tmp_Vx_wake,tmp_Vr_wake,i,np1,errStat2) SHARED(m,u,p,maxPln,errStat,errMsg) DEFAULT(NONE)
   !do nz_low=0, p%nZ_low-1
   !   do ny_low=0, p%yZ_low-1
   !      do nx_low=0, p%nX_low-1
   !         nXYZ_low=0
   do i = 0 , p%NumGrid_low - 1

            nx_low = mod(     i                        ,p%nX_low)
            ny_low = mod(int( i / (p%nX_low         ) ),p%nY_low)
            nz_low =     int( i / (p%nX_low*p%nY_low) )

               ! set the disturbed flow equal to the ambient flow for this time step
            m%Vdist_low(:,nx_low,ny_low,nz_low) = m%Vamb_low(:,nx_low,ny_low,nz_low)

            !nXYZ_low = nXYZ_low + 1
            nXYZ_low = i + 1
            n_wake = 0
            xhatBar_plane = 0.0_ReKi

            do nt = 1,p%NumTurbines

               ! H Long: replace intrinsic dot_product with explicit do product can save as much as 10% of total calculation time!
               !x_end_plane = dot_product(u%xhat_plane(:,0,nt), (p%Grid_Low(:,nXYZ_low) - u%p_plane(:,0,nt)) )
               tmp_x = u%xhat_plane(1,0,nt) * (p%Grid_Low(1,nXYZ_low) - u%p_plane(1,0,nt))
               tmp_y = u%xhat_plane(2,0,nt) * (p%Grid_Low(2,nXYZ_low) - u%p_plane(2,0,nt))
               tmp_z = u%xhat_plane(3,0,nt) * (p%Grid_Low(3,nXYZ_low) - u%p_plane(3,0,nt))
               x_end_plane = tmp_x + tmp_y + tmp_z

               do np = 0, maxPln
                     ! Reset interpolation counter
                  ILo = 0
                  np1 = np + 1
                     ! Construct the endcaps of the current wake plane volume
                  x_start_plane = x_end_plane
                  ! H Long: again, replace intrinsic dot_product
                  !x_end_plane = dot_product(u%xhat_plane(:,np+1,nt), (p%Grid_Low(:,nXYZ_low) - u%p_plane(:,np+1,nt)) )
                  tmp_x = u%xhat_plane(1,np1,nt) * (p%Grid_Low(1,nXYZ_low) - u%p_plane(1,np1,nt))
                  tmp_y = u%xhat_plane(2,np1,nt) * (p%Grid_Low(2,nXYZ_low) - u%p_plane(2,np1,nt))
                  tmp_z = u%xhat_plane(3,np1,nt) * (p%Grid_Low(3,nXYZ_low) - u%p_plane(3,np1,nt))


                  x_end_plane = tmp_x + tmp_y + tmp_z

                     ! test if the point is within the endcaps of the wake volume

                  if ( ( ( x_start_plane >= 0.0_ReKi ) .and. ( x_end_plane < 0.0_ReKi ) ) .or. &
                       ( ( x_start_plane <= 0.0_ReKi ) .and. ( x_end_plane > 0.0_ReKi ) )        ) then

                     if ( EqualRealNos( x_start_plane, x_end_plane ) ) then
                        delta = 0.5_ReKi
                     else
                        delta = x_start_plane / ( x_start_plane - x_end_plane )
                     end if
                     deltad = (1.0_ReKi - delta)
                     if ( m%parallelFlag(np,nt) ) then
                        p_tmp_plane = delta*u%p_plane(:,np1,nt) + deltad*u%p_plane(:,np,nt)
                     else
                        tmp_vec = delta*m%rhat_e(:,np,nt) + deltad*m%rhat_s(:,np,nt)
                        p_tmp_plane = delta*m%pvec_ce(:,np,nt) + deltad*m%pvec_cs(:,np,nt) + ( delta*m%r_e(np,nt) + deltad*m%r_s(np,nt) )* tmp_vec / TwoNorm(tmp_vec)
                     end if



                     r_vec_plane = p%Grid_Low(:,nXYZ_low) - p_tmp_plane
                     r_tmp_plane = TwoNorm( r_vec_plane )

                        ! test if the point is within radial finite-difference grid
                     if ( r_tmp_plane <= p%r(p%numRadii-1) ) then

                        n_wake = n_wake + 1


                        if ( EqualRealNos(r_tmp_plane, 0.0_ReKi) ) then
                           tmp_rhat_plane(:,n_wake) = 0.0_ReKi
                        else
                           tmp_rhat_plane(:,n_wake) = ( r_vec_plane  ) / r_tmp_plane
                        end if


                           ! given r_tmp_plane and Vx_wake at p%dr increments, find value of m%Vx_wake(@r_tmp_plane) using interpolation
                        tmp_Vx_wake(n_wake) = delta*InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np1,nt), ILo, p%NumRadii ) + deltad*InterpBin( r_tmp_plane, p%r, u%Vx_wake(:,np,nt), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
                        tmp_Vr_wake(n_wake) = delta*InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np1,nt), ILo, p%NumRadii ) + deltad*InterpBin( r_tmp_plane, p%r, u%Vr_wake(:,np,nt), ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )


                        tmp_xhat_plane(:,n_wake) = delta*u%xhat_plane(:,np1,nt) + deltad*u%xhat_plane(:,np,nt)
                        tmp_xhat_plane(:,n_wake) = tmp_xhat_plane(:,n_wake) / TwoNorm(tmp_xhat_plane(:,n_wake))
                        xhatBar_plane = xhatBar_plane + abs(tmp_Vx_wake(n_wake))*tmp_xhat_plane(:,n_wake)
                     end if  ! if the point is within radial finite-difference grid
                end if
            end do  ! do np = 0, p%NumPlanes-2
        end do      ! do nt = 1,p%NumTurbines

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
              Vr_term     = tmp_Vx_wake(nw)*tmp_xhat_plane(:,nw) + tmp_Vr_wake(nw)*tmp_rhat_plane(:,nw)
              Vx_term     = dot_product( xhatBar_plane, Vr_term )
              Vx_wake_tmp = Vx_wake_tmp + Vx_term*Vx_term
              Vr_wake_tmp = Vr_wake_tmp + Vr_term
           end do
                  ! [I - XX']V = V - (V dot X)X
           Vr_wake_tmp = Vr_wake_tmp - dot_product(Vr_wake_tmp,xhatBar_plane)*xhatBar_plane
           m%Vdist_low(:,nx_low,ny_low,nz_low) = m%Vdist_low(:,nx_low,ny_low,nz_low) + real(Vr_wake_tmp - xhatBar_plane*sqrt(Vx_wake_tmp),SiKi)
        end if  ! (n_wake > 0)
   end do
   !      end do ! do nx_low=0, p%nX_low-1
   !   end do    ! do ny_low=0, p%nY_low-1
   !end do       ! do nz_low=0, p%nZ_low-1
   !$OMP END PARALLEL DO

   do nt = 1,p%NumTurbines
         
      do np = 0,tmpPln 
      
      !!Defining yhat and zhat
         xxplane = (/u%xhat_plane(1,np,nt), 0.0_ReKi, 0.0_ReKi/)
         xyplane = (/0.0_ReKi, u%xhat_plane(1,np,nt), 0.0_ReKi/)
         yyplane = (/0.0_ReKi, u%xhat_plane(2,np,nt), 0.0_ReKi/)
         yxplane = (/u%xhat_plane(2,np,nt), 0.0_ReKi, 0.0_ReKi/)
         xyplane_norm = TwoNorm(xxplane+yyplane)

         IF (EqualRealNos(xyplane_norm, 0.0_ReKi)) THEN ! This should only be true during the first call to AWAE_CalcOutput at model initialization

            y%Vx_wind_disk(nt) = 0.0_ReKi
            y%TI_amb(      nt) = 0.0_ReKi
            y%V_plane(:,np,nt) = 0.0_ReKi

         ELSE                                           ! All subsequent calls to AWAE_CalcOutput


      ! Warn our kind users if wake planes leave the low-resolution domain:
            if ( u%p_plane(1,np,nt) < p%Grid_Low(1,            1) ) call SetErrStat(ErrID_Warn, 'The center of wake plane #'//trim(num2lstr(np))//' for turbine #'//trim(num2lstr(nt))//' has passed the lowest-most X boundary of the low-resolution domain.', errStat, errMsg, RoutineName)
            if ( u%p_plane(1,np,nt) > p%Grid_Low(1,p%NumGrid_low) ) call SetErrStat(ErrID_Warn, 'The center of wake plane #'//trim(num2lstr(np))//' for turbine #'//trim(num2lstr(nt))//' has passed the upper-most X boundary of the low-resolution domain.' , errStat, errMsg, RoutineName)
            if ( u%p_plane(2,np,nt) < p%Grid_Low(2,            1) ) call SetErrStat(ErrID_Warn, 'The center of wake plane #'//trim(num2lstr(np))//' for turbine #'//trim(num2lstr(nt))//' has passed the lowest-most Y boundary of the low-resolution domain.', errStat, errMsg, RoutineName)
            if ( u%p_plane(2,np,nt) > p%Grid_Low(2,p%NumGrid_low) ) call SetErrStat(ErrID_Warn, 'The center of wake plane #'//trim(num2lstr(np))//' for turbine #'//trim(num2lstr(nt))//' has passed the upper-most Y boundary of the low-resolution domain.' , errStat, errMsg, RoutineName)
            if ( u%p_plane(3,np,nt) < p%Grid_Low(3,            1) ) call SetErrStat(ErrID_Warn, 'The center of wake plane #'//trim(num2lstr(np))//' for turbine #'//trim(num2lstr(nt))//' has passed the lowest-most Z boundary of the low-resolution domain.', errStat, errMsg, RoutineName)
            if ( u%p_plane(3,np,nt) > p%Grid_Low(3,p%NumGrid_low) ) call SetErrStat(ErrID_Warn, 'The center of wake plane #'//trim(num2lstr(np))//' for turbine #'//trim(num2lstr(nt))//' has passed the upper-most Z boundary of the low-resolution domain.' , errStat, errMsg, RoutineName)
         
         
             xplane_sq = u%xhat_plane(1,np,nt)**2.0_ReKi
             yplane_sq = u%xhat_plane(2,np,nt)**2.0_ReKi
             xysq_Z = (/0.0_ReKi, 0.0_ReKi, xplane_sq+yplane_sq/)
             xzplane_X = (/u%xhat_plane(1,np,nt)*u%xhat_plane(3,np,nt), 0.0_ReKi, 0.0_ReKi/)
             yzplane_Y = (/0.0_ReKi, u%xhat_plane(2,np,nt)*u%xhat_plane(3,np,nt), 0.0_ReKi/)
             tmp_yhat_plane = (xyplane-yxplane)/xyplane_norm
             tmp_zhat_plane = (xysq_Z-xzplane_X-yzplane_Y)/xyplane_norm


             ! Calculate y%Vx_wind_disk and y%TI_amb at the rotor disk

             if ( np == 0 ) then

                Vsum_low  = 0.0_ReKi
                iwsum = 0
                n_r_polar = FLOOR((p%C_Meander*u%D_wake(np,nt))/(2.0_ReKi*p%dpol))

                do nr = 0,n_r_polar

                   r_polar = REAL(nr,ReKi)*p%dpol
                   n_psi_polar = MAX(CEILING(TwoPi*REAL(nr,ReKi))-1,0)

                   do npsi = 0,n_psi_polar

                      psi_polar = (TwoPi*REAL(npsi,ReKi))/(REAL(n_psi_polar+1,ReKi))
                      p_polar = u%p_plane(:,np,nt) + r_polar*COS(psi_polar)*tmp_yhat_plane + r_polar*SIN(psi_polar)*tmp_zhat_plane
                      Vamb_lowpol_tmp = INTERP3D( p_polar, p%Grid_Low(:,1), p%dXYZ_Low, m%Vamb_low, within, p%nX_low, p%nY_low, p%nZ_low, Vbox=Vamb_low_tmp )
                      if ( within ) then
                         Vsum_low = Vsum_low + Vamb_lowpol_tmp
                         do i = 1,8
                            iwsum = iwsum + 1
                            m%Vamb_lowpol(:,iwsum) = Vamb_low_tmp(:,i)
                         end do
                      end if

                   end do

                end do

                if ( iwsum == 0 ) then

                   call SetErrStat( ErrID_Fatal, 'The rotor plane for turbine '//trim(num2lstr(nt))//' has left the low-resolution domain (i.e., there are no points in the polar grid that lie within the low-resolution domain).', errStat, errMsg, RoutineName )
                   return

                else

                   Vsum_low = Vsum_low/REAL(iwsum/8,ReKi)   ! iwsum is always a multiple of 8
                   Vave_amb_low_norm  = TwoNorm(Vsum_low)
                   if ( EqualRealNos(Vave_amb_low_norm, 0.0_ReKi ) )  then
                      call SetErrStat( ErrID_Fatal, 'The magnitude of the spatial-averaged ambient wind speed in the low-resolution domain associated with the wake plane at the rotor disk for turbine #'//trim(num2lstr(nt))//' is zero.', errStat, errMsg, RoutineName )
                      return
                   else
                      y%Vx_wind_disk(nt) = dot_product( u%xhat_plane(:,np,nt),Vsum_low )
                      y%TI_amb(nt) = 0.0_ReKi
                      do wamb = 1, iwsum
                         y%TI_amb(nt) = y%TI_amb(nt)+TwoNorm(m%Vamb_lowpol(:,wamb)-Vsum_low)**2.0_ReKi
                      end do  !wamb
                      y%TI_amb(nt) = sqrt(y%TI_amb(nt)/(3.0_ReKi*REAL(iwsum,ReKi)))/Vave_amb_low_norm
                   end if !Vave_amb_low_norm

                 end if

             end if


             ! Calculate y%V_plane

             y%V_plane(:,np,nt) = 0.0_ReKi
             wsum_tmp = 0.0_ReKi
             n_r_polar = FLOOR((p%C_ScaleDiam*u%D_wake(np,nt))/p%dpol)

             do nr = 0, n_r_polar

                r_polar = REAL(nr,ReKi)*p%dpol

                select case ( p%Mod_Meander )
                case (MeanderMod_Uniform)
                   w = 1.0_ReKi
                case (MeanderMod_TruncJinc)
                   w = jinc( r_polar/(p%C_Meander*u%D_wake(np,nt) ) )
                case (MeanderMod_WndwdJinc)
                   w = jinc( r_polar/(p%C_Meander*u%D_wake(np,nt) ) )*jinc( r_polar/(2.0_ReKi*p%C_Meander*u%D_wake(np,nt) ) )
                end select

                n_psi_polar = MAX(CEILING(TwoPi*REAL(nr,ReKi))-1,0)

                do npsi = 0,n_psi_polar

                   psi_polar = (TwoPi*REAL(npsi,ReKi))/(REAL(n_psi_polar+1,ReKi))
                   p_polar = u%p_plane(:,np,nt) + r_polar*COS(psi_polar)*tmp_yhat_plane + r_polar*SIN(psi_polar)*tmp_zhat_plane
                   Vdist_lowpol_tmp = INTERP3D( p_polar, p%Grid_Low(:,1), p%dXYZ_Low, m%Vdist_low, within, p%nX_low, p%nY_low, p%nZ_low )
                   if ( within ) then
                      y%V_plane(:,np,nt) = y%V_plane(:,np,nt) + w*Vdist_lowpol_tmp
                      wsum_tmp = wsum_tmp + w
                   end if

                end do !npsi

             end do!nr

             if ( EqualRealNos( wsum_tmp, 0.0_ReKi ) ) then
                y%V_plane(:,np,nt) = 0.0_ReKi
             else
                y%V_plane(:,np,nt) = y%V_plane(:,np,nt)/wsum_tmp
             end if

         end if

      end do
   end do

!#ifdef _OPENMP  
!   tm2 =  omp_get_wtime() 
!   write(*,*)  'Total AWAE:LowResGridCalcOutput using '//trim(num2lstr(tm2-tm1))//' seconds'

!#endif 

   if (allocated(tmp_xhat_plane)) deallocate(tmp_xhat_plane)
   if (allocated(tmp_rhat_plane)) deallocate(tmp_rhat_plane)
   if (allocated(tmp_Vx_wake)) deallocate(tmp_Vx_wake)
   if (allocated(tmp_Vr_wake)) deallocate(tmp_Vr_wake)

end subroutine LowResGridCalcOutput


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
      n_high_low = 0
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
                        if ( ( ( x_start_plane >= 0.0_ReKi ) .and. ( x_end_plane < 0.0_ReKi ) ) .or. &
                             ( ( x_start_plane <= 0.0_ReKi ) .and. ( x_end_plane > 0.0_ReKi ) )        ) then

                           if ( EqualRealNos( x_start_plane, x_end_plane ) ) then
                              delta = 0.5_ReKi
                           else
                              delta = x_start_plane / ( x_start_plane - x_end_plane )
                           end if
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
                  do n_hl=0, n_high_low
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
   character(1024)                               :: rootDir, baseName, OutFileVTKDir ! Simulation root dir, basename for outputs
   integer(IntKi)                                :: i,j,nt        ! loop counter
   integer(IntKi)                                :: maxN_wake
   real(ReKi)                                    :: gridLoc       ! Location of requested output slice in grid coordinates [0,sz-1]
   integer(IntKi)                                :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                          :: errMsg2       ! temporary error message
   character(*), parameter                       :: RoutineName = 'AWAE_Init'
   type(InflowWind_InitInputType)                :: IfW_InitInp
   type(InflowWind_InitOutputType)               :: IfW_InitOut
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
   p%Mod_AmbWind      = InitInp%InputFileData%Mod_AmbWind
   p%NumPlanes        = InitInp%InputFileData%NumPlanes
   p%NumRadii         = InitInp%InputFileData%NumRadii
   p%NumTurbines      = InitInp%InputFileData%NumTurbines
   p%WindFilePath     = InitInp%InputFileData%WindFilePath ! TODO: Make sure this wasn't specified with the trailing folder separator. Note: on Windows a trailing / or \ causes no problem! GJH
   p%n_high_low       = InitInp%n_high_low
   p%dt_low           = InitInp%InputFileData%dt_low
   p%NumDT            = InitInp%NumDT
   p%NOutDisWindXY    = InitInp%InputFileData%NOutDisWindXY
   p%NOutDisWindYZ    = InitInp%InputFileData%NOutDisWindYZ
   p%NOutDisWindXZ    = InitInp%InputFileData%NOutDisWindXZ
   p%WrDisWind        = InitInp%InputFileData%WrDisWind
   p%WrDisSkp1        = nint(InitInp%InputFileData%WrDisDT / p%dt_low)
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

   ! --- Vtk Outputs
   call GetPath( p%OutFileRoot, rootDir, baseName ) 
   OutFileVTKDir    = trim(rootDir) // 'vtk_ff'  ! Directory for VTK outputs
   p%OutFileVTKRoot = trim(rootDir) // 'vtk_ff' // PathSep // trim(baseName) ! Basename for VTK files
   p%VTK_tWidth = CEILING( log10( real(p%NumDT, ReKi)/real(p%WrDisSkp1, ReKi) ) + 1) ! Length for time stamp
   if (p%WrDisWind .or. p%NOutDisWindXY>0 .or. p%NOutDisWindYZ>0 .or. p%NOutDisWindXZ>0) then
      call MKDIR(OutFileVTKDir) ! creating output directory
   end if


   allocate( p%r(0:p%NumRadii-1),stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for p%r.', errStat, errMsg, RoutineName )
         return
      end if

   do i = 0,p%NumRadii-1
      p%r(i)       = InitInp%InputFileData%dr*i
   end do

   allocate( p%WT_Position(3,p%NumTurbines),stat=errStat2)
      if (errStat2 /= 0) then
         call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for p%WT_Position.', errStat, errMsg, RoutineName )
         return
      end if
   p%WT_Position = InitInp%InputFileData%WT_Position

      ! Obtain the precursor grid information by parsing the necessary input files
      ! This will establish certain parameters as well as all of the initialization outputs
      ! Sets:
      ! Parameters: nX_low, nY_low, nZ_low, nX_high, nY_high, nZ_high, Grid_low,
      !             Grid_high, n_high_low, n_rp_max
      ! InitOutput: X0_high, Y0_high, Z0_high, dX_high, dY_high, dZ_high, nX_high, nY_high, nZ_high


   call AWAE_IO_InitGridInfo(InitInp, p, InitOut, errStat2, errMsg2)
      call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if (errStat2 >= AbortErrLev) then
         return
   end if

   if ( p%Mod_AmbWind > 1 ) then
      ! Using InflowWind, so initialize that module now
      IfW_InitInp%Linearize         = .false.
      IfW_InitInp%RootName          = TRIM(p%OutFileRoot)//'.IfW'
      IfW_InitInp%UseInputFile      = .TRUE.
      IfW_InitInp%InputFileName     = InitInp%InputFileData%InflowFile
      IfW_InitInp%lidar%Tmax        = 0.0_ReKi
      IfW_InitInp%lidar%HubPosition = 0.0_ReKi
      IfW_InitInp%lidar%SensorType  = SensorType_None
      IfW_InitInp%Use4Dext          = .false.

      if (      p%Mod_AmbWind == 2 ) then ! one InflowWind module

         ALLOCATE(p%IfW(0:0),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind parameter data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(x%IfW(0:0),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind continuous states data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(xd%IfW(0:0),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind discrete states data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(z%IfW(0:0),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind constraint states data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(OtherState%IfW(0:0),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind other states data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(m%IfW(0:0),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind miscvar data', ErrStat, ErrMsg, RoutineName )
            return
         end if

         ! Initialize InflowWind
         IfW_InitInp%FixedWindFileRootName = .false.
         IfW_InitInp%NumWindPoints         = p%NumGrid_low
      
         call InflowWind_Init( IfW_InitInp, m%u_IfW_Low, p%IfW(0), x%IfW(0), xd%IfW(0), z%IfW(0), OtherState%IfW(0), m%y_IfW_Low, m%IfW(0), Interval, IfW_InitOut, ErrStat2, ErrMsg2 )
            call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
         if (errStat2 >= AbortErrLev) then
            return
         end if

      else if ( p%Mod_AmbWind == 3 ) then ! multiple InflowWind modules

         ALLOCATE(p%IfW(0:p%NumTurbines),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind parameter data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(x%IfW(0:p%NumTurbines),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind continuous states data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(xd%IfW(0:p%NumTurbines),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind discrete states data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(z%IfW(0:p%NumTurbines),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind constraint states data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(OtherState%IfW(0:p%NumTurbines),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind other states data', ErrStat, ErrMsg, RoutineName )
            return
         end if
         ALLOCATE(m%IfW(0:p%NumTurbines),STAT=ErrStat2)
         if (ErrStat2 /= 0) then
            CALL SetErrStat( ErrID_Fatal, 'Could not allocate memory for InflowWind miscvar data', ErrStat, ErrMsg, RoutineName )
            return
         end if

         ! Initialize InflowWind for the low-resolution domain
         IfW_InitInp%FixedWindFileRootName = .true.
         IfW_InitInp%NumWindPoints         = p%NumGrid_low
         IfW_InitInp%TurbineID             = 0
      
         call InflowWind_Init( IfW_InitInp, m%u_IfW_Low, p%IfW(0), x%IfW(0), xd%IfW(0), z%IfW(0), OtherState%IfW(0), m%y_IfW_Low, m%IfW(0), Interval, IfW_InitOut, ErrStat2, ErrMsg2 )
            call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
         if (errStat2 >= AbortErrLev) then
            return
         end if

         ! Initialize InflowWind for each high-resolution domain
         IfW_InitInp%NumWindPoints         = p%nX_high*p%nY_high*p%nZ_high
         
         do nt = 1,p%NumTurbines

            IfW_InitInp%TurbineID          = nt
      
            call InflowWind_Init( IfW_InitInp, m%u_IfW_High, p%IfW(nt), x%IfW(nt), xd%IfW(nt), z%IfW(nt), OtherState%IfW(nt), m%y_IfW_High, m%IfW(nt), Interval, IfW_InitOut, ErrStat2, ErrMsg2 )
               call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
            if (errStat2 >= AbortErrLev) then
               return
            end if

         end do

      end if

         ! Set the position inputs once for the low-resolution grid
      m%u_IfW_Low%PositionXYZ = p%Grid_low

         ! Initialize the high-resolution grid inputs and outputs
       IF ( .NOT. ALLOCATED( m%u_IfW_High%PositionXYZ ) ) THEN
         call AllocAry(m%u_IfW_High%PositionXYZ, 3, p%nX_high*p%nY_high*p%nZ_high, 'm%u_IfW_High%PositionXYZ', ErrStat2, ErrMsg2)
            call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
         call AllocAry(m%y_IfW_High%VelocityUVW, 3, p%nX_high*p%nY_high*p%nZ_high, 'm%y_IfW_High%VelocityUVW', ErrStat2, ErrMsg2)
            call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
         call AllocAry(m%y_IfW_High%WriteOutput, size(m%y_IfW_Low%WriteOutput), 'm%y_IfW_High%WriteOutput', ErrStat2, ErrMsg2)
            call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
      END IF
      if (errStat2 >= AbortErrLev) then
            return
      end if

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
      if (errStat >= AbortErrLev) return


   !interval = InitOut%dt_low

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
   if (errStat >= AbortErrLev) return




      !............................................................................................
      ! Define outputs here
      !............................................................................................

   allocate ( y%V_plane(3,0:p%NumPlanes-1,1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%V_plane.', errStat, errMsg, RoutineName )
   allocate ( y%Vdist_High(1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%Vdist_High.', errStat, errMsg, RoutineName )
   do i = 1, p%NumTurbines
      allocate ( y%Vdist_High(i)%data(3,0:p%nX_high-1,0:p%nY_high-1,0:p%nZ_high-1,0:p%n_high_low), STAT=ErrStat2 )
         if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%Vdist_High%data.', errStat, errMsg, RoutineName )
      y%Vdist_High(i)%data    = 0.0_Siki
   end do

   allocate ( y%Vx_wind_disk   (1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%Vx_rel_disk.', errStat, errMsg, RoutineName )
   allocate ( y%TI_amb   (1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for y%TI_amb.', errStat, errMsg, RoutineName )
   if (errStat >= AbortErrLev) return

      ! This next step is not strictly necessary
   y%V_plane       = 0.0_Reki
   y%Vx_wind_disk  = 0.0_Reki
   y%TI_amb        = 0.0_Reki


   if ( p%NOutDisWindXY > 0 ) then
      ALLOCATE ( m%OutVizXYPlane(3,p%nX_low, p%nY_low,1) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error allocating memory for the Fast.Farm OutVizXYPlane arrays.'
         RETURN
      ENDIF
   end if
   if ( p%NOutDisWindYZ > 0 ) then
      ALLOCATE ( m%OutVizYZPlane(3,p%nY_low, p%nZ_low,1) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error allocating memory for the Fast.Farm OutVizYZPlane arrays.'
         RETURN
      ENDIF
   end if
   if ( p%NOutDisWindXZ > 0 ) then
      ALLOCATE ( m%OutVizXZPlane(3,p%nX_low, p%nZ_low,1) , STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 )  THEN
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
   allocate ( m%Vamb_lowpol   ( 3, 0:p%n_rp_max*8 ) , STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vamb_lowpol.', errStat, errMsg, RoutineName )
   allocate ( m%Vdist_low  ( 3, 0:p%nX_low-1 , 0:p%nY_low-1 , 0:p%nZ_low-1 )                  , STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vdist_low.', errStat, errMsg, RoutineName )

   allocate ( m%Vamb_high(1:p%NumTurbines), STAT=ErrStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vamb_high.', errStat, errMsg, RoutineName )
   do i = 1, p%NumTurbines
         allocate ( m%Vamb_high(i)%data(3,0:p%nX_high-1,0:p%nY_high-1,0:p%nZ_high-1,0:p%n_high_low), STAT=ErrStat2 )
            if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vamb_high%data.', errStat, errMsg, RoutineName )
   end do
   maxN_wake = p%NumTurbines*( p%NumPlanes-1 )
   allocate ( m%xhat_plane ( 3, 1:maxN_wake ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%xhat_plane.', errStat, errMsg, RoutineName )
   allocate ( m%rhat_plane ( 3, 1:maxN_wake ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%rhat_plane.', errStat, errMsg, RoutineName )
   allocate ( m%Vx_wake    ( 1:maxN_wake ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vx_wake.', errStat, errMsg, RoutineName )
   allocate ( m%Vr_wake    ( 1:maxN_wake ), STAT=errStat2 )
      if (errStat2 /= 0) call SetErrStat ( ErrID_Fatal, 'Could not allocate memory for m%Vr_wake.', errStat, errMsg, RoutineName )

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
   if (errStat >= AbortErrLev) return


   ! Read-in the ambient wind data for the initial calculate output

   call AWAE_UpdateStates( 0.0_DbKi, -1, u, p, x, xd, z, OtherState, m, errStat2, errMsg2 )
   call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )





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
      integer(IntKi),                 intent(  out)  :: errStat     !< Error status of the operation
      character(*),                   intent(  out)  :: errMsg      !< Error message if errStat /= ErrID_None

         ! Local variables
      integer(IntKi)                                 :: nt          !< loop counter


         ! Initialize errStat

      errStat = ErrID_None
      errMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:
 
         ! End all instances of the InflowWind module
      if (      p%Mod_AmbWind == 2 ) then
         call    InflowWind_End( m%u_IfW_Low, p%IfW(0 ), x%IfW(0 ), xd%IfW(0 ), z%IfW(0 ), OtherState%IfW(0 ), m%y_IfW_Low, m%IfW(0 ), errStat, errMsg )
      else if ( p%Mod_AmbWind == 3 ) then
         call    InflowWind_End( m%u_IfW_Low, p%IfW(0 ), x%IfW(0 ), xd%IfW(0 ), z%IfW(0 ), OtherState%IfW(0 ), m%y_IfW_Low, m%IfW(0 ), errStat, errMsg )
         do nt = 1,p%NumTurbines
            call InflowWind_End( m%u_IfW_Low, p%IfW(nt), x%IfW(nt), xd%IfW(nt), z%IfW(nt), OtherState%IfW(nt), m%y_IfW_Low, m%IfW(nt), errStat, errMsg )
         end do
      end if

      

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
!   real(DbKi)          :: t1, t2
   integer(IntKi)                               :: n_high_low, nt, n_hl, i,j,k,c
   
   errStat = ErrID_None
   errMsg  = ""
   
   ! Read the ambient wind data that is needed for t+dt, i.e., n+1
!#ifdef _OPENMP
!   t1 = omp_get_wtime()  
!#endif 
   
   if ( (n+1) == (p%NumDT-1) ) then
      n_high_low = 0
   else
      n_high_low = p%n_high_low
   end if

   if ( p%Mod_AmbWind == 1 ) then
         ! read from file the ambient flow for the n+1 time step
      call ReadLowResWindFile(n+1, p, m%Vamb_Low, errStat2, errMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
         if (errStat >= AbortErrLev) return 
   !#ifdef _OPENMP
   !   t2 = omp_get_wtime()      
   !   write(*,*) '        AWAE_UpdateStates: Time spent reading Low Res data : '//trim(num2lstr(t2-t1))//' seconds'            
   !#endif   
      
      do nt = 1,p%NumTurbines
         do n_hl=0, n_high_low
               ! read from file the ambient flow for the current time step
            call ReadHighResWindFile(nt, (n+1)*p%n_high_low + n_hl, p, m%Vamb_high(nt)%data(:,:,:,:,n_hl), errStat2, errMsg2)
               call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
               if (errStat >= AbortErrLev) return 
         end do
      end do

   else ! p%Mod_AmbWind == 2 .or. 3

      ! Set low-resolution inflow wind velocities
      call InflowWind_CalcOutput(t+p%dt_low, m%u_IfW_Low, p%IfW(0), x%IfW(0), xd%IfW(0), z%IfW(0), OtherState%IfW(0), m%y_IfW_Low, m%IfW(0), errStat2, errMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
         if (errStat >= AbortErrLev) return 
      c = 1
      do k = 0,p%nZ_low-1
         do j = 0,p%nY_low-1
            do i = 0,p%nX_low-1
               m%Vamb_Low(:,i,j,k) = m%y_IfW_Low%VelocityUVW(:,c)
               c = c+1
            end do
         end do
      end do
      ! Set the high-resoultion inflow wind velocities for each turbine
      if (  p%Mod_AmbWind == 2 ) then
         do nt = 1,p%NumTurbines
            m%u_IfW_High%PositionXYZ = p%Grid_high(:,:,nt)
            do n_hl=0, n_high_low
               call InflowWind_CalcOutput(t+p%dt_low+n_hl*p%DT_high, m%u_IfW_High, p%IfW(0), x%IfW(0), xd%IfW(0), z%IfW(0), OtherState%IfW(0), m%y_IfW_High, m%IfW(0), errStat2, errMsg2)
                  call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
                  if (errStat >= AbortErrLev) return 
               c = 1
               do k = 0,p%nZ_high-1
                  do j = 0,p%nY_high-1
                     do i = 0,p%nX_high-1
                        m%Vamb_high(nt)%data(:,i,j,k,n_hl) = m%y_IfW_High%VelocityUVW(:,c)
                        c = c+1
                     end do
                  end do
               end do
            end do
         end do

      else !p%Mod_AmbWind == 3
         do nt = 1,p%NumTurbines
            c = 1
            do k = 0,p%nZ_high-1
               do j = 0,p%nY_high-1
                  do i = 0,p%nX_high-1
                     m%u_IfW_High%PositionXYZ(:,c) = p%Grid_high(:,c,nt) - p%WT_Position(:,nt)
                     c = c+1
                  end do
               end do
            end do
            do n_hl=0, n_high_low
               call InflowWind_CalcOutput(t+p%dt_low+n_hl*p%DT_high, m%u_IfW_High, p%IfW(nt), x%IfW(nt), xd%IfW(nt), z%IfW(nt), OtherState%IfW(nt), m%y_IfW_High, m%IfW(nt), errStat2, errMsg2)
                  call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, RoutineName )
                  if (errStat >= AbortErrLev) return 
               c = 1
               do k = 0,p%nZ_high-1
                  do j = 0,p%nY_high-1
                     do i = 0,p%nX_high-1
                        m%Vamb_high(nt)%data(:,i,j,k,n_hl) = m%y_IfW_High%VelocityUVW(:,c)
                        c = c+1
                     end do
                  end do
               end do
            end do
         end do
      end if

   end if

!#ifdef _OPENMP
!   t1 = omp_get_wtime()      
!   write(*,*) '        AWAE_UpdateStates: Time spent reading High Res data : '//trim(num2lstr(t1-t2))//' seconds'             
!#endif 
   
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
   character(p%VTK_tWidth)                      :: Tstr          ! string for current VTK write-out step (padded with zeros)
   integer(intKi)                               :: i, j, k
   integer(intKi)                               :: errStat2
   character(ErrMsgLen)                         :: errMsg2
   character(*), parameter                      :: RoutineName = 'AWAE_CalcOutput'
   integer(intKi)                               :: n, n_high
   CHARACTER(1024)                              :: FileName
   INTEGER(IntKi)                               :: Un                   ! unit number of opened file


   errStat = ErrID_None
   errMsg  = ""
   n = nint(t / p%dt_low)
   call ComputeLocals(n, u, p, y, m, errStat2, errMsg2)
      call SetErrStat ( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if (errStat2 >= AbortErrLev) then
            return
      end if
   call LowResGridCalcOutput(n, u, p, y, m, errStat2, errMsg2)
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

      ! TimeStamp
      write(Tstr, '(i' // trim(Num2LStr(p%VTK_tWidth)) //'.'// trim(Num2LStr(p%VTK_tWidth)) // ')') n/p%WrDisSkp1 ! TODO use n instead..

         ! XY plane slices
      do k = 1,p%NOutDisWindXY

         call ExtractSlice( XYSlice, p%OutDisWindZ(k), p%Z0_low, p%nZ_low, p%nX_low, p%nY_low, p%dZ_low, m%Vdist_low, m%outVizXYPlane(:,:,:,1))
            ! Create the output vtk file with naming <WindFilePath>/Low/DisXY<k>.t<n/p%WrDisSkp1>.vtk
         FileName = trim(p%OutFileVTKRoot)//".Low.DisXY"//trim(num2lstr(k))//"."//trim(Tstr)//".vtk"
         call WrVTK_SP_header( FileName, "Low resolution, disturbed wind of XY Slice at time = "//trim(num2lstr(t))//" seconds.", Un, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
         call WrVTK_SP_vectors3D( Un, "Velocity", (/p%nX_low,p%nY_low,1_IntKi/), (/p%X0_low,p%Y0_low,p%OutDisWindZ(k)/), (/p%dX_low,p%dY_low,p%dZ_low/), m%outVizXYPlane, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end do


         ! YZ plane slices
      do k = 1,p%NOutDisWindYZ
         call ExtractSlice( YZSlice, p%OutDisWindX(k), p%X0_low, p%nX_low, p%nY_low, p%nZ_low, p%dX_low, m%Vdist_low, m%outVizYZPlane(:,:,:,1))
            ! Create the output vtk file with naming <WindFilePath>/Low/DisYZ<k>.t<n/p%WrDisSkp1>.vtk
         FileName = trim(p%OutFileVTKRoot)//".Low.DisYZ"//trim(num2lstr(k))//"."//trim(Tstr)//".vtk"
         call WrVTK_SP_header( FileName, "Low resolution, disturbed wind of YZ Slice at time = "//trim(num2lstr(t))//" seconds.", Un, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
         call WrVTK_SP_vectors3D( Un, "Velocity", (/1,p%nY_low,p%nZ_low/), (/p%OutDisWindX(k),p%Y0_low,p%Z0_low/), (/p%dX_low,p%dY_low,p%dZ_low/), m%outVizYZPlane, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end do

         ! XZ plane slices
      do k = 1,p%NOutDisWindXZ
         call ExtractSlice( XZSlice, p%OutDisWindY(k), p%Y0_low, p%nY_low, p%nX_low, p%nZ_low, p%dY_low, m%Vdist_low, m%outVizXZPlane(:,:,:,1))
            ! Create the output vtk file with naming <WindFilePath>/Low/DisXZ<k>.t<n/p%WrDisSkp1>.vtk
         FileName = trim(p%OutFileVTKRoot)//".Low.DisXZ"//trim(num2lstr(k))//"."//trim(Tstr)//".vtk"
         call WrVTK_SP_header( FileName, "Low resolution, disturbed wind of XZ Slice at time = "//trim(num2lstr(t))//" seconds.", Un, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
         call WrVTK_SP_vectors3D( Un, "Velocity", (/p%nX_low,1,p%nZ_low/), (/p%X0_low,p%OutDisWindY(k),p%Z0_low/), (/p%dX_low,p%dY_low,p%dZ_low/), m%outVizXZPlane, ErrStat2, ErrMsg2 )
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

   if ( (InputFileData%Mod_AmbWind < 1) .or. (InputFileData%Mod_AmbWind > 3) ) call SetErrStat ( ErrID_Fatal, 'Mod_AmbWind must be 1: high-fidelity precursor in VTK format, 2: one instance of InflowWind module, or 3: multiple instances of InflowWind module.', errStat, errMsg, RoutineName )
   if ( InputFileData%Mod_AmbWind == 1 ) then
      if (len_trim(InputFileData%WindFilePath) == 0) call SetErrStat ( ErrID_Fatal, 'WindFilePath must contain at least one character.', errStat, errMsg, RoutineName )
   else
      if (len_trim(InputFileData%InflowFile) == 0) call SetErrStat ( ErrID_Fatal, 'InflowFile must contain at least one character.', errStat, errMsg, RoutineName )
      if ( (InputFileData%nX_low < 2) .or. (InputFileData%nY_low < 2) .or. (InputFileData%nZ_low < 2) ) &
         call SetErrStat ( ErrID_Fatal, 'The low resolution grid dimensions must contain a minimum of 2 nodes in each spatial direction. ', errStat, errMsg, RoutineName )
      if ( (InputFileData%nX_high < 2) .or. (InputFileData%nY_high < 2) .or. (InputFileData%nY_high < 2) ) &
         call SetErrStat ( ErrID_Fatal, 'The high resolution grid dimensions must contain a minimum of 2 nodes in each spatial direction. ', errStat, errMsg, RoutineName )
      if ( (InputFileData%dX_low <= 0.0_ReKi) .or. (InputFileData%dY_low <= 0.0_ReKi) .or. (InputFileData%dY_low <= 0.0_ReKi) ) &
         call SetErrStat ( ErrID_Fatal, 'The low resolution spatial resolution must be greater than zero in each spatial direction. ', errStat, errMsg, RoutineName )
   end if

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
   InitInp%InputFileData%dt_low         = 2.0_DbKi
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

FUNCTION INTERP3D(p,p0,del,V,within,nX,nY,nZ,Vbox)
      !  I/O variables
         Real(ReKi), INTENT( IN    ) :: p(3)            !< Position where the 3D velocity field will be interpreted (m)
         Real(ReKi), INTENT( IN    ) :: p0(3)           !< Origin of the spatial domain (m)
         Real(ReKi), INTENT( IN    ) :: del(3)          !< XYZ-components of the spatial increment of the domain (m)
         INTEGER(IntKi), INTENT( IN) :: nX, nY, nZ      !< Size of XYZ spatial dimensions
         Real(SiKi), INTENT( IN    ) :: V(3,0:nX-1,0:nY-1,0:nZ-1)        !< 3D velocity field to be interpolated

         Real(SiKi) :: INTERP3D(3)     !Vint(3)         !< Interpolated velocity (m/s)
         Logical,    INTENT(   OUT ) :: within          !< Logical flag indicating weather or not the input position lies within the domain (flag)
         REAL(ReKi), OPTIONAL, INTENT(OUT) :: Vbox(3,8) !< Wind velocities at the 8 points in the 3D spatial domain surrounding the input position

      !  Local variables
         INTEGER(IntKi)        :: i
         Real(ReKi)            :: f(3), N(8)
         Real(SiKi)            :: Vtmp(3,8)
         INTEGER(IntKi)        :: n_lo(3), n_hi(3)


     !!! CHECK BOUNDS
   within = .TRUE.
   do i = 1, 3
      f(i) = (p(i)-p0(i))/del(i)
      n_lo(i) = FLOOR(f(i))
      n_hi(i) = n_lo(i)+1_IntKi
      f(i) = 2.0_ReKi*( f(i)-REAL(n_lo(i),ReKi) )-1.0_ReKi  ! convert to value between -1 and 1
      if (( n_lo(i) < 0) .OR. (n_hi(i) > size(V,i+1)-1)) THEN
         within = .FALSE.
      END IF
   end do

     !!! INTERPOLATE
   INTERP3D = 0.0_SiKi
   if (within) then
      
      N(1) = ((1.0_ReKi-f(1))*(1.0_ReKi-f(2))*(1.0_ReKi-f(3)))/8.0_ReKi
      N(2) = ((1.0_ReKi+f(1))*(1.0_ReKi-f(2))*(1.0_ReKi-f(3)))/8.0_ReKi
      N(3) = ((1.0_ReKi-f(1))*(1.0_ReKi+f(2))*(1.0_ReKi-f(3)))/8.0_ReKi
      N(4) = ((1.0_ReKi+f(1))*(1.0_ReKi+f(2))*(1.0_ReKi-f(3)))/8.0_ReKi
      N(5) = ((1.0_ReKi-f(1))*(1.0_ReKi-f(2))*(1.0_ReKi+f(3)))/8.0_ReKi
      N(6) = ((1.0_ReKi+f(1))*(1.0_ReKi-f(2))*(1.0_ReKi+f(3)))/8.0_ReKi
      N(7) = ((1.0_ReKi-f(1))*(1.0_ReKi+f(2))*(1.0_ReKi+f(3)))/8.0_ReKi
      N(8) = ((1.0_ReKi+f(1))*(1.0_ReKi+f(2))*(1.0_ReKi+f(3)))/8.0_ReKi
      Vtmp(:,1) = V(:,n_lo(1),n_lo(2),n_lo(3))
      Vtmp(:,2) = V(:,n_hi(1),n_lo(2),n_lo(3))
      Vtmp(:,3) = V(:,n_lo(1),n_hi(2),n_lo(3))
      Vtmp(:,4) = V(:,n_hi(1),n_hi(2),n_lo(3))
      Vtmp(:,5) = V(:,n_lo(1),n_lo(2),n_hi(3))
      Vtmp(:,6) = V(:,n_hi(1),n_lo(2),n_hi(3))
      Vtmp(:,7) = V(:,n_lo(1),n_hi(2),n_hi(3))
      Vtmp(:,8) = V(:,n_hi(1),n_hi(2),n_hi(3))

      do i=1,8

         ! To support complex terrain, the wind data will have NaNs at any point in the domain below the ground; throw away these points.
         if ( Is_NaN( REAL(vtmp(1,i),DbKi) ) .OR. Is_NaN( REAL(vtmp(2,i),DbKi) ) .OR. Is_NaN( REAL(vtmp(3,i),DbKi) ) ) then
            within = .FALSE.
            INTERP3D(:) = 0.0_SiKi
            EXIT
         end if

         INTERP3D(:) = INTERP3D(:) + N(i)*Vtmp(:,i)

      end do

   else

      vtmp = 0.0_SiKi

   end if

   ! Output the wind velocities at the 8 points in the 3D spatial domain surrounding the input position (if necessary)
   IF ( PRESENT( Vbox ) ) Vbox = REAL( Vtmp, ReKi )

END FUNCTION

end module AWAE
