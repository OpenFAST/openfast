!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
!
!    This file is part of AeroDyn.
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
MODULE AeroDyn_IO
 
   use NWTC_Library
   use AeroDyn_IO_Params
   use AeroDyn_Types
   use BEMTUncoupled, only : VelocityIsZero
   use FVW_Subs,      only : FVW_AeroOuts

   USE AeroDyn_AllBldNdOuts_IO
   
   implicit none

   type(ProgDesc), parameter  :: AD_Ver = ProgDesc( 'AeroDyn', '', '' )
   character(*),   parameter  :: AD_Nickname = 'AD'
      

  
     
contains
   
!> Compute maximum radius over all blades (contains hub radius), in "projected rotor plane"
!! Solely based on AD inputs,  needed for FVW since rLocal is not stored
PURE REAL(ReKi) FUNCTION Calc_MaxRadius(p, u) result(rmax)
   implicit none
   TYPE(RotParameterType),    INTENT(IN   ) :: p    !< The module parameters
   TYPE(RotInputType),        INTENT(IN   ) :: u    !< Inputs
   real(ReKi)     :: y_hat_disk(3), z_hat_disk(3), dr_gl(3), rLocal
   integer(IntKi) :: iB, j
   y_hat_disk = u%HubMotion%Orientation(2,:,1)
   z_hat_disk = u%HubMotion%Orientation(3,:,1)
   rmax = 0.0_ReKi
   do iB=1,p%numBlades
      do j=1,p%NumBlNds
         dr_gl =  u%BladeMotion(iB)%Position(:,j) - u%HubMotion%Position(:,1) ! vector hub center to node j in global coord
         rLocal = sqrt( dot_product(dr_gl, y_hat_disk)**2 + dot_product(dr_gl, z_hat_disk)**2 )
         rmax   = max(rmax, rLocal)
      end do !j=nodes
   end do !iB=blades
END FUNCTION Calc_MaxRadius

!> Rotor speed
PURE REAL(ReKi) FUNCTION Calc_Omega(u)
   TYPE(RotInputType),        INTENT(IN   ) :: u    !< Inputs
   Calc_Omega = dot_product(u%HubMotion%RotationVel(:,1), u%HubMotion%Orientation(1,:,1))
END FUNCTION Calc_Omega

!> Mean skew angle
REAL(ReKi) FUNCTION Calc_Chi0(V_diskAvg, V_dot_x)
   implicit none
   REAL(ReKi), INTENT(IN  ) :: V_diskAvg(3)
   REAL(ReKi), INTENT(IN  ) :: V_dot_x
   REAL(ReKi) :: V_norm, sy 
   V_norm = TwoNorm( V_diskAvg )
   if ( EqualRealNos( V_norm, 0.0_ReKi ) ) then
      Calc_Chi0 = 0.0_ReKi
   else
      ! make sure we don't have numerical issues that make the ratio outside +/-1
      sy = min(  1.0_ReKi, V_dot_x / V_norm )
      sy = max( -1.0_ReKi, sy )
      Calc_Chi0 = acos( sy )
   end if
END FUNCTION Calc_Chi0


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Calc_WriteOutput( p, p_AD, u, RotInflow, x, m, m_AD, y, OtherState, xd, indx, iRot, ErrStat, ErrMsg )
   
   TYPE(RotParameterType),       INTENT(IN   )  :: p                                 ! The rotor parameters
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p_AD                              ! The module parameters
   TYPE(RotInputType),           INTENT(IN   )  :: u                                 ! inputs
   TYPE(RotInflowType),          INTENT(IN   )  :: RotInflow                         ! other states%RotInflow at t (for DBEMT and UA)
   TYPE(RotContinuousStateType), INTENT(IN   )  :: x                                 !< Continuous states at t
   TYPE(RotMiscVarType),         INTENT(INOUT)  :: m                                 ! misc variables
   TYPE(AD_MiscVarType),         INTENT(INOUT)  :: m_AD                              ! misc variables
   TYPE(RotOutputType),          INTENT(IN   )  :: y                                 ! outputs
   TYPE(RotOtherStateType),      INTENT(IN   )  :: OtherState                        ! other states at t (for DBEMT and UA)
   TYPE(RotDiscreteStateType),   INTENT(IN   )  :: xd                                ! Discrete states
   integer,                      intent(in   )  :: indx                              ! index into m%BEMT_u(indx) array; 1=t and 2=t+dt (but not checked here)
   integer,                      intent(in   )  :: iRot                              ! Rotor index, needed for FVW
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat                           ! The error status code
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg                            ! The error message, if an error occurred

      ! local variables
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Calc_WriteOutput'
   INTEGER(intKi)                               :: ErrStat2
   CHARACTER(ErrMsgLen)                         :: ErrMsg2
   
   INTEGER(IntKi)                               :: j,k,beta
   REAL(ReKi)                                   :: tmp(3)
   REAL(ReKi)                                   :: tmpHubFB(3)
   REAL(ReKi)                                   :: tmpHubMB(3)   
   REAL(ReKi)                                   :: force(3)
   REAL(ReKi)                                   :: moment(3)
   REAL(ReKi)                                   :: denom, rmax, omega
   REAL(ReKi)                                   :: ct, st ! cosine, sine of theta
   REAL(ReKi)                                   :: cp, sp ! cosine, sine of phi
      
   
      ! start routine:
   ErrStat = ErrID_None
   ErrMsg  = ""
   

   ! Compute max radius and rotor speed
   if (p%NumBlades == 0) then
      rmax  = 0.0_ReKi
      omega = 0.0_ReKi
   elseif (p_AD%Wake_Mod /= WakeMod_FVW) then
      rmax = 0.0_ReKi
      do k=1,p%NumBlades
         do j=1,p%NumBlNds
            rmax = max(rmax, m%BEMT_u(indx)%rLocal(j,k) )
         end do !j=nodes
      end do !k=blades
      
!     rmax  = p%BEMT%rTipFixMax
      omega = m%BEMT_u(indx)%omega
   else
      rmax  = Calc_MaxRadius(p, u)
      omega = Calc_Omega(u)
   endif

   
   ! Common outputs to all AeroDyn submodules
   call Calc_WriteOutput_AD() ! need to call this before calling the BEMT vs FVW versions of outputs so that the integrated output quantities are known
   
   if (p_AD%Wake_Mod /= WakeMod_FVW) then
      call Calc_WriteOutput_BEMT()
   else
      call Calc_WriteOutput_FVW()
   endif

      ! set these for debugging
!   m%AllOuts( Debug1 ) = 0.0_ReKi !TwoNorm( m%BEMT%u_SkewWake(1)%v_qsw )
!   m%AllOuts( Debug2 ) = 0.0_ReKi !TwoNorm( x%BEMT%v_w )
!   m%AllOuts( Debug3 ) = 0.0_ReKi
   
CONTAINS
   !..........................................................................................
   subroutine Calc_WriteOutput_AD()
   
         ! tower outputs
      do beta=1,p%NTwOuts
         j = p%TwOutNd(beta)
      
         tmp = matmul( u%TowerMotion%Orientation(:,:,j) , RotInflow%Tower%InflowVel(:,j) )
         m%AllOuts( TwNVUnd(:,beta) ) = tmp
      
         tmp = matmul( u%TowerMotion%Orientation(:,:,j) , u%TowerMotion%TranslationVel(:,j) )
         m%AllOuts( TwNSTV(:,beta) ) = tmp
      
         m%AllOuts( TwNVrel(beta) ) = m%W_Twr(j)                           ! relative velocity   
         m%AllOuts( TwNDynP(beta) ) = 0.5 * p%AirDens * m%W_Twr(j)**2      ! dynamic pressure
         m%AllOuts( TwNRe(  beta) ) = p%TwrDiam(j) * m%W_Twr(j) / p%KinVisc / 1.0E6 ! reynolds number (in millions)
         m%AllOuts( TwNM(   beta) ) = m%W_Twr(j) / p%SpdSound               ! Mach number
         m%AllOuts( TwNFdx( beta) ) = m%X_Twr(j)         
         m%AllOuts( TwNFdy( beta) ) = m%Y_Twr(j)         
      end do ! out nodes
      
      if ( p%Buoyancy ) then
         do beta=1,p%NTwOuts
            j = p%TwOutNd(beta)
            
            tmp = matmul( u%TowerMotion%Orientation(:,:,j) , m%TwrBuoyLoad%Force(:,j) )
            m%AllOuts( TwNFbx(beta) ) = tmp(1)
            m%AllOuts( TwNFby(beta) ) = tmp(2)
            m%AllOuts( TwNFbz(beta) ) = tmp(3)
   
            tmp = matmul( u%TowerMotion%Orientation(:,:,j) , m%TwrBuoyLoad%Moment(:,j) )
            m%AllOuts( TwNMbx(beta) ) = tmp(1)
            m%AllOuts( TwNMby(beta) ) = tmp(2)
            m%AllOuts( TwNMbz(beta) ) = tmp(3)
         end do
      end if
      
         ! hub outputs
      if ( p%Buoyancy ) then
         tmpHubFB = matmul( u%HubMotion%Orientation(:,:,1) , m%HubFB )
         m%AllOuts( HbFbx ) = tmpHubFB(1)
         m%AllOuts( HbFby ) = tmpHubFB(2)
         m%AllOuts( HbFbz ) = tmpHubFB(3)
   
         tmpHubMB = matmul( u%HubMotion%Orientation(:,:,1) , m%HubMB )
         m%AllOuts( HbMbx ) = tmpHubMB(1)
         m%AllOuts( HbMby ) = tmpHubMB(2)
         m%AllOuts( HbMbz ) = tmpHubMB(3)
      else
         tmpHubFB = 0.0_ReKi ! initialize for integration later
         tmpHubMB = 0.0_ReKi ! initialize for integration later
      end if
   
         ! nacelle buoyancy outputs
      if ( p%Buoyancy ) then
         tmp = matmul( u%NacelleMotion%Orientation(:,:,1) , m%NacFB )
         m%AllOuts( NcFbx ) = tmp(1)
         m%AllOuts( NcFby ) = tmp(2)
         m%AllOuts( NcFbz ) = tmp(3)
   
         tmp = matmul( u%NacelleMotion%Orientation(:,:,1) , m%NacMB )
         m%AllOuts( NcMbx ) = tmp(1)
         m%AllOuts( NcMby ) = tmp(2)
         m%AllOuts( NcMbz ) = tmp(3)
      end if

         ! nacelle drag outputs
      if ( p%NacelleDrag ) then

         tmp = matmul( u%NacelleMotion%Orientation(:,:,1) , m%NacDragF )
         m%AllOuts( NcFdx ) = tmp(1)
         m%AllOuts( NcFdy ) = tmp(2)
         m%AllOuts( NcFdz ) = tmp(3)
   
         tmp = matmul( u%NacelleMotion%Orientation(:,:,1) , m%NacDragM )
         m%AllOuts( NcMdx ) = tmp(1)
         m%AllOuts( NcMdy ) = tmp(2)
         m%AllOuts( NcMdz ) = tmp(3)
      end if

         ! nacelle total forces and moments
      if ( p%Buoyancy .OR. p%NacelleDrag) then

         tmp = m%NacFi
         m%AllOuts( NcFxi ) = tmp(1)
         m%AllOuts( NcFyi ) = tmp(2)
         m%AllOuts( NcFzi ) = tmp(3)

         tmp = m%NacMi
         m%AllOuts( NcMxi ) = tmp(1)
         m%AllOuts( NcMyi ) = tmp(2)
         m%AllOuts( NcMzi ) = tmp(3)
      end if

         ! blade outputs
      do k=1,min(p%numBlades,AD_MaxBl_Out)    ! limit this
         do beta=1,p%NBlOuts
            j=p%BlOutNd(beta)

            tmp = matmul( m%orientationAnnulus(:,:,j,k), RotInflow%Blade(k)%InflowVel(:,j) )
            m%AllOuts( BNVUndx(beta,k) ) = tmp(1)
            m%AllOuts( BNVUndy(beta,k) ) = tmp(2)
            m%AllOuts( BNVUndz(beta,k) ) = tmp(3)

            tmp = matmul( m%orientationAnnulus(:,:,j,k), m%DisturbedInflow(:,j,k) )
            m%AllOuts( BNVDisx(beta,k) ) = tmp(1)
            m%AllOuts( BNVDisy(beta,k) ) = tmp(2)
            m%AllOuts( BNVDisz(beta,k) ) = tmp(3)

            tmp = matmul( m%orientationAnnulus(:,:,j,k), u%BladeMotion(k)%TranslationVel(:,j) )
            m%AllOuts( BNSTVx( beta,k) ) = tmp(1)
            m%AllOuts( BNSTVy( beta,k) ) = tmp(2)
            m%AllOuts( BNSTVz( beta,k) ) = tmp(3)
         
            m%AllOuts( BNCurve(beta,k) ) = m%Cant(j,k)*R2D
                  
            m%AllOuts( BNSigCr(   beta,k) ) = m%SigmaCavitCrit(j,k)
            m%AllOuts( BNSgCav(   beta,k) ) = m%SigmaCavit(j,k)

         end do ! nodes
      end do ! blades

      if ( p%Buoyancy ) then
         do k=1,min(p%numBlades,AD_MaxBl_Out)    ! limit this
            do beta=1,p%NBlOuts
               j=p%BlOutNd(beta)

               tmp = matmul( u%BladeMotion(k)%Orientation(:,:,j), m%BladeBuoyLoad(k)%Force(:,j) )
               m%AllOuts( BNFbn(beta,k) ) = tmp(1)
               m%AllOuts( BNFbt(beta,k) ) = tmp(2)
               m%AllOuts( BNFbs(beta,k) ) = tmp(3)

               tmp = matmul( u%BladeMotion(k)%Orientation(:,:,j), m%BladeBuoyLoad(k)%Moment(:,j) )
               m%AllOuts( BNMbn(beta,k) ) = tmp(1)
               m%AllOuts( BNMbt(beta,k) ) = tmp(2)
               m%AllOuts( BNMbs(beta,k) ) = tmp(3)
            end do ! nodes
         end do ! blades
      end if
      
      ! blade node tower clearance (requires tower influence calculation):
      if (p%TwrPotent /= TwrPotent_none .or. p%TwrShadow /= TwrShadow_none) then
         do k=1,min(p%numBlades,AD_MaxBl_Out)
            do beta=1,p%NBlOuts
               j=p%BlOutNd(beta)
               m%AllOuts( BNClrnc( beta,k) ) = m%TwrClrnc(j,k)
            end do
         end do
      end if

   
   

      m%AllOuts( RtSpeed ) = omega*RPS2RPM
      m%AllOuts( RtArea  ) = pi * rmax**2
      
      tmp = matmul( u%HubMotion%Orientation(:,:,1), m%V_DiskAvg )
      m%AllOuts( RtVAvgxh ) = tmp(1)
      m%AllOuts( RtVAvgyh ) = tmp(2)
      m%AllOuts( RtVAvgzh ) = tmp(3)

      

         ! integrate force/moments over blades by performing mesh transfer to hub point:
      force  = tmpHubFB
      moment = tmpHubMB
      do k=1,p%NumBlades
         call Transfer_Line2_to_Point( y%BladeLoad(k), m%HubLoad, m%B_L_2_H_P(k), ErrStat2, ErrMsg2, u%BladeMotion(k), u%HubMotion )
         force  = force  + m%HubLoad%force( :,1)
         moment = moment + m%HubLoad%moment(:,1)
         
         if (k<=size(BAeroFxi)) then
            ! Power contribution of blade wrt hub
            tmp = matmul( u%HubMotion%Orientation(:,:,1), m%HubLoad%moment(:,1) )
            m%AllOuts( BAeroPwr(k) ) = omega * tmp(1)
            
            ! In global, wrt hub! 
            m%AllOuts( BAeroFxi(k) ) = m%HubLoad%force(1,1)
            m%AllOuts( BAeroFyi(k) ) = m%HubLoad%force(2,1)
            m%AllOuts( BAeroFzi(k) ) = m%HubLoad%force(3,1)
            m%AllOuts( BAeroMxi(k) ) = m%HubLoad%moment(1,1)
            m%AllOuts( BAeroMyi(k) ) = m%HubLoad%moment(2,1)
            m%AllOuts( BAeroMzi(k) ) = m%HubLoad%moment(3,1)
         end if
      end do

        ! In global
      m%AllOuts( RtAeroFxi ) = force(1)
      m%AllOuts( RtAeroFyi ) = force(2)
      m%AllOuts( RtAeroFzi ) = force(3)
      m%AllOuts( RtAeroMxi ) = moment(1)
      m%AllOuts( RtAeroMyi ) = moment(2)
      m%AllOuts( RtAeroMzi ) = moment(3)
      tmp = matmul( u%HubMotion%Orientation(:,:,1), force )
      m%AllOuts( RtAeroFxh ) = tmp(1)
      m%AllOuts( RtAeroFyh ) = tmp(2)
      m%AllOuts( RtAeroFzh ) = tmp(3)
   
      tmp = matmul( u%HubMotion%Orientation(:,:,1), moment )
      m%AllOuts( RtAeroMxh ) = tmp(1)
      m%AllOuts( RtAeroMyh ) = tmp(2)
      m%AllOuts( RtAeroMzh ) = tmp(3)
      
      m%AllOuts( RtAeroPwr ) = omega * m%AllOuts( RtAeroMxh )
      
     
   
      ! Integrate force/moments over blades by performing mesh transfer to blade root points:
      do k=1,p%NumBlades
         call Transfer_Line2_to_Point( y%BladeLoad(k), m%BladeRootLoad(k), m%B_L_2_R_P(k), ErrStat2, ErrMsg2, u%BladeMotion(k), u%BladeRootMotion(k) )
      end do
      do k=1,min(p%NumBlades,size(BAeroFx))
         ! Transform force vector to blade root coordinate system
         tmp = matmul( u%BladeRootMotion(k)%Orientation(:,:,1), m%BladeRootLoad(k)%force( :,1) )
         m%AllOuts( BAeroFx(k) ) = tmp(1)
         m%AllOuts( BAeroFy(k) ) = tmp(2)
         m%AllOuts( BAeroFz(k) ) = tmp(3)
      
         ! Transform moment vector to blade root coordinate system
         tmp = matmul( u%BladeRootMotion(k)%Orientation(:,:,1), m%BladeRootLoad(k)%moment( :,1) )
         m%AllOuts( BAeroMx(k) ) = tmp(1)
         m%AllOuts( BAeroMy(k) ) = tmp(2)
         m%AllOuts( BAeroMz(k) ) = tmp(3)
      end do  ! k=blades
   
         ! rotor outputs
      if ( abs( m%V_dot_x ) < 0.04_ReKi .or. p%NumBlades == 0 ) then ! < 0.018 is close to "v_dot_x**3 not equal to 0" (cubed because of Cp denominator)
         m%AllOuts( RtAeroCp ) = 0.0_ReKi
         m%AllOuts( RtAeroCq ) = 0.0_ReKi
         m%AllOuts( RtAeroCt ) = 0.0_ReKi
      else
         denom = 0.5*p%AirDens*m%AllOuts( RtArea )*m%V_dot_x**2
         !denom = 0.5 * p%AirDens * (pi * rmax**2) * m%V_dot_x**2
      
         m%AllOuts( RtAeroCp ) = m%AllOuts( RtAeroPwr ) / (denom * m%V_dot_x)
         m%AllOuts( RtAeroCq ) = m%AllOuts( RtAeroMxh ) / (denom * rmax )
         m%AllOuts( RtAeroCt ) = m%AllOuts( RtAeroFxh ) /  denom
      end if
      

      ! TailFin
      if (p%TFinAero) then
         m%AllOuts( TFAlpha )           = m%TFinAlpha*R2D ! [deg]
         m%AllOuts( TFRe )              = m%TFinRe
         !m%AllOuts( TFM )              = m%TFinVrel/p%SpdSound
         m%AllOuts( TFVrel )            = m%TFinVrel
         m%AllOuts( TFVundxi:TFVundzi)  = m%TFinVund_i
         m%AllOuts( TFVindxi:TFVindzi ) = m%TFinVind_i
         m%AllOuts( TFVrelxi:TFVrelzi ) = m%TFinVrel_i
         m%AllOuts( TFSTVxi:TFSTVzi )   = m%TFinSTV_i
         m%AllOuts( TFFxi:TFFzi )       = m%TFinF_i
         m%AllOuts( TFMxi:TFMzi )       = m%TFinM_i
      endif

            
   end subroutine Calc_WriteOutput_AD
   !..........................................................................................
   subroutine Calc_WriteOutput_BEMT()
      REAL(R8Ki)                                   :: orient(3,3)
      REAL(R8Ki)                                   :: theta(3)
      REAL(ReKi)                                   :: Vind_s(3)  ! Induced velocity in "w" or "p" system
      REAL(ReKi)                                   :: denom !, rmax
      REAL(ReKi)                                   :: ct, st ! cosine, sine of theta
      REAL(ReKi)                                   :: cp, sp ! cosine, sine of phi
 

      ! Induced velocity in Global
      do k=1,min(p%numBlades,3)
         do j=1,u%BladeMotion(k)%NNodes
            !if(p%BEM_Mod==BEMMod_2D) then
            ! NOTE: if BEMMod_2D:   x & y are in "w" system (WithoutSweepPitchTwist)
            !       if BEMMod_3D:   x & y are in "l" system (local-polar system)
            Vind_s = (/ -m%BEMT_u(indx)%Vx(j,k)*m%BEMT_y%axInduction(j,k), m%BEMT_u(indx)%Vy(j,k)*m%BEMT_y%tanInduction(j,k), 0.0_ReKi /)
            m%Vind_i(:,j,k) = matmul(Vind_s, m%orientationAnnulus(:,:,j,k)) ! TODO rename orientationAnnulus
         enddo
      enddo


   
         ! blade outputs
      do k=1,min(p%numBlades,AD_MaxBl_Out)    ! limit this
         m%AllOuts( BAzimuth(k) ) = MODULO( m%BEMT_u(indx)%psi_s(k)*R2D, 360.0_ReKi )
       ! m%AllOuts( BPitch(  k) ) = calculated in SetInputsForBEMT
      
         do beta=1,p%NBlOuts
         
            j=p%BlOutNd(beta)
                  
            m%AllOuts( BNVrel( beta,k) ) = m%BEMT_y%Vrel(j,k)
            m%AllOuts( BNDynP( beta,k) ) = 0.5 * p%airDens * m%BEMT_y%Vrel(j,k)**2
            m%AllOuts( BNRe(   beta,k) ) = p%BEMT%chord(j,k) * m%BEMT_y%Vrel(j,k) / p%KinVisc / 1.0E6
            m%AllOuts( BNM(    beta,k) ) = m%BEMT_y%Vrel(j,k) / p%SpdSound

            m%AllOuts( BNVIndx(beta,k) ) = - m%BEMT_u(indx)%Vx(j,k) * m%BEMT_y%axInduction( j,k)
            m%AllOuts( BNVIndy(beta,k) ) =   m%BEMT_u(indx)%Vy(j,k) * m%BEMT_y%tanInduction(j,k)

            m%AllOuts( BNAxInd(beta,k) ) = m%BEMT_y%axInduction(j,k)
            m%AllOuts( BNTnInd(beta,k) ) = m%BEMT_y%tanInduction(j,k)

            m%AllOuts( BNAlpha(beta,k) ) = Rad2M180to180Deg( m%BEMT_y%phi(j,k) - m%BEMT_u(indx)%theta(j,k) )
            m%AllOuts( BNTheta(beta,k) ) = m%BEMT_u(indx)%theta(j,k)*R2D
            m%AllOuts( BNPhi(  beta,k) ) = m%BEMT_y%phi(j,k)*R2D

            m%AllOuts( BNCpmin(   beta,k) ) = m%BEMT_y%Cpmin(j,k)
   !         m%AllOuts( BNSigCr(   beta,k) ) = m%SigmaCavitCrit(j,k)
   !         m%AllOuts( BNSgCav(   beta,k) ) = m%SigmaCavit(j,k)

            !m%AllOuts( BNCl(   beta,k) ) = m%BEMT_y%Cl(j,k)
            !m%AllOuts( BNCd(   beta,k) ) = m%BEMT_y%Cd(j,k)
            cp=cos(m%BEMT_y%phi(j,k))
            sp=sin(m%BEMT_y%phi(j,k))
            m%AllOuts( BNCl(   beta,k) ) = m%BEMT_y%Cx(j,k)*cp + m%BEMT_y%Cy(j,k)*sp
            m%AllOuts( BNCd(   beta,k) ) = m%BEMT_y%Cx(j,k)*sp - m%BEMT_y%Cy(j,k)*cp
            m%AllOuts( BNCm(   beta,k) ) = m%BEMT_y%Cm(j,k)
            m%AllOuts( BNCx(   beta,k) ) = m%BEMT_y%Cx(j,k)
            m%AllOuts( BNCy(   beta,k) ) = m%BEMT_y%Cy(j,k)

            ct=cos(m%BEMT_u(indx)%theta(j,k))
            st=sin(m%BEMT_u(indx)%theta(j,k))
            m%AllOuts( BNCn(   beta,k) ) = m%BEMT_y%Cx(j,k)*ct + m%BEMT_y%Cy(j,k)*st
            m%AllOuts( BNCt(   beta,k) ) =-m%BEMT_y%Cx(j,k)*st + m%BEMT_y%Cy(j,k)*ct

            m%AllOuts( BNFl(   beta,k) ) =  m%X(j,k)*cp - m%Y(j,k)*sp
            m%AllOuts( BNFd(   beta,k) ) =  m%X(j,k)*sp + m%Y(j,k)*cp
            m%AllOuts( BNMm(   beta,k) ) =  m%M(j,k)
            m%AllOuts( BNFx(   beta,k) ) =  m%X(j,k)
            m%AllOuts( BNFy(   beta,k) ) = -m%Y(j,k)
            m%AllOuts( BNFn(   beta,k) ) =  m%X(j,k)*ct - m%Y(j,k)*st
            m%AllOuts( BNFt(   beta,k) ) = -m%X(j,k)*st - m%Y(j,k)*ct

            m%AllOuts( BNGam(  beta,k) ) = 0.5_ReKi * p%BEMT%chord(j,k) * m%BEMT_y%Vrel(j,k) * m%BEMT_y%Cl(j,k) ! "Gam" [m^2/s]
            
         end do ! nodes
      end do ! blades
   

      ! rotor outputs:
      if (p%NumBlades > 0) then
         m%AllOuts( RtSkew   ) = m%BEMT_u(indx)%chi0*R2D
         m%AllOuts( RtTSR    ) = m%BEMT_u(indx)%TSR
         m%AllOuts( DBEMTau1 ) = OtherState%BEMT%DBEMT%tau1
      else
         m%AllOuts( RtSkew   ) = 0.0_ReKi
         m%AllOuts( RtTSR    ) = 0.0_ReKi
         m%AllOuts( DBEMTau1 ) = 0.0_ReKi
      end if
      
   end subroutine Calc_WriteOutput_BEMT

   !..........................................................................................
   !> Similar to Calc_WriteOutput_BEMT. TODO Merge me
   !! NOTE: relies on the prior calculation of m%V_dot_x, and m%V_diskAvg (done in DiskAvgValues)
   !!                                          m%DisturbedInflow (done in SetInputs)
   !!       Make sure these are set!
   subroutine Calc_WriteOutput_FVW
      integer    :: iW

      ! Induced velocity in global
      ! FVW already return this, we do a simple copy from Wings to Blades
      do k=1,min(p%numBlades,3)
         iW = p_AD%FVW%Bld2Wings(iRot, k)
         do j=1,u%BladeMotion(k)%NNodes
            m%Vind_i(:,j,k) = m_AD%FVW_y%W(iW)%Vind(1:3,j)
         enddo
      enddo

      ! TODO TODO TODO ALL THIS SHOULD BE COMPUTED IN THE SAME MEMORY FORMAT AS AERODYN

         ! blade outputs
      do k=1,min(p%numBlades,3)
         iW=p_AD%FVW%Bld2Wings(iRot, k)

         do beta=1,p%NBlOuts
            j=p%BlOutNd(beta)

            m%AllOuts( BNVrel( beta,k) ) = m_AD%FVW%W(iW)%BN_Vrel(j)
            m%AllOuts( BNDynP( beta,k) ) = 0.5 * p%airDens * m_AD%FVW%W(iW)%BN_Vrel(j)**2
            m%AllOuts( BNRe(   beta,k) ) = m_AD%FVW%W(iW)%BN_Re(j)  / 1.0E6
            m%AllOuts( BNM(    beta,k) ) = m_AD%FVW%W(iW)%BN_Vrel(j) / p%SpdSound

            m%AllOuts( BNVIndx(beta,k) ) = -m_AD%FVW%W(iW)%BN_UrelWind_s(1,j) * m_AD%FVW%W(iW)%BN_AxInd(j)
            m%AllOuts( BNVIndy(beta,k) ) =  m_AD%FVW%W(iW)%BN_UrelWind_s(2,j) * m_AD%FVW%W(iW)%BN_TanInd(j)

            m%AllOuts( BNAxInd(beta,k) ) = m_AD%FVW%W(iW)%BN_AxInd(j)
            m%AllOuts( BNTnInd(beta,k) ) = m_AD%FVW%W(iW)%BN_TanInd(j)

            m%AllOuts( BNAlpha(beta,k) ) = m_AD%FVW%W(iW)%BN_alpha(j)*R2D
            m%AllOuts( BNTheta(beta,k) ) = m_AD%FVW%W(iW)%PitchAndTwist(j)*R2D
            m%AllOuts( BNPhi(  beta,k) ) = m_AD%FVW%W(iW)%BN_phi(j)*R2D

            m%AllOuts( BNCpmin(beta,k) ) = m_AD%FVW%W(iW)%BN_Cpmin(j)
            m%AllOuts( BNCl(   beta,k) ) = m_AD%FVW%W(iW)%BN_Cl(j)
            m%AllOuts( BNCd(   beta,k) ) = m_AD%FVW%W(iW)%BN_Cd(j)
            m%AllOuts( BNCm(   beta,k) ) = m_AD%FVW%W(iW)%BN_Cm(j)
            m%AllOuts( BNCx(   beta,k) ) = m_AD%FVW%W(iW)%BN_Cx(j)
            m%AllOuts( BNCy(   beta,k) ) = m_AD%FVW%W(iW)%BN_Cy(j)

            ct=cos(m_AD%FVW%W(iW)%PitchAndTwist(j))    ! cos(theta)
            st=sin(m_AD%FVW%W(iW)%PitchAndTwist(j))    ! sin(theta)
            m%AllOuts( BNCn(   beta,k) ) = m_AD%FVW%W(iW)%BN_Cx(j)*ct + m_AD%FVW%W(iW)%BN_Cy(j)*st
            m%AllOuts( BNCt(   beta,k) ) =-m_AD%FVW%W(iW)%BN_Cx(j)*st + m_AD%FVW%W(iW)%BN_Cy(j)*ct

            cp=cos(m_AD%FVW%W(iW)%BN_phi(j))
            sp=sin(m_AD%FVW%W(iW)%BN_phi(j))
            m%AllOuts( BNFl(   beta,k) ) =  m%X(j,k)*cp - m%Y(j,k)*sp
            m%AllOuts( BNFd(   beta,k) ) =  m%X(j,k)*sp + m%Y(j,k)*cp
            m%AllOuts( BNMm(   beta,k) ) =  m%M(j,k)
            m%AllOuts( BNFx(   beta,k) ) =  m%X(j,k)
            m%AllOuts( BNFy(   beta,k) ) = -m%Y(j,k)
            m%AllOuts( BNFn(   beta,k) ) =  m%X(j,k)*ct - m%Y(j,k)*st
            m%AllOuts( BNFt(   beta,k) ) = -m%X(j,k)*st - m%Y(j,k)*ct

            m%AllOuts( BNGam(  beta,k) ) = 0.5_ReKi * p_AD%FVW%W(iW)%chord_LL(j) * m_AD%FVW%W(iW)%BN_Vrel(j) * m_AD%FVW%W(iW)%BN_Cl(j) ! "Gam" [m^2/s]
         end do ! nodes
      end do ! blades


!     m%AllOuts( RtArea  ) = pi*rmax**2     ! TODO vertical axis
      m%AllOuts( RtSkew  ) = Calc_Chi0(m%V_diskAvg, m%V_dot_x) * R2D 

      if ( EqualRealNos( REAL(m%V_dot_x, SiKi), 0.0_SiKi ) ) then
        m%AllOuts( RtTSR )    = 0.0_ReKi
      else
        m%AllOuts( RtTSR )    = omega * rmax / m%V_dot_x
      end if
      m%AllOuts( DBEMTau1 ) = 0.0_ReKi ! not valid with FVW

   end subroutine Calc_WriteOutput_FVW

   
END SUBROUTINE Calc_WriteOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadInputFiles( InputFileName, InputFileData, Default_DT, OutFileRoot, NumBlades, AeroProjMod, UnEcho, calcCrvAngle, ErrStat, ErrMsg )
! This subroutine reads the input file and stores all the data in the AD_InputFile structure.
! It does not perform data validation.
!..................................................................................................................................

      ! Passed variables
   REAL(DbKi),              INTENT(IN)    :: Default_DT      ! The default DT (from glue code)

   CHARACTER(*),            INTENT(IN)    :: InputFileName   ! Name of the input file
   CHARACTER(*),            INTENT(IN)    :: OutFileRoot     ! The rootname of all the output files written by this routine.

   TYPE(AD_InputFile),      INTENT(INOUT) :: InputFileData   ! Data stored in the module's input file
   INTEGER(IntKi),          INTENT(INOUT) :: UnEcho          ! Unit number for the echo file

   INTEGER(IntKi),          INTENT(IN)    :: NumBlades(:)    ! Number of blades for this model per rotor
   INTEGER(IntKi),          INTENT(IN)    :: AeroProjMod(:)  ! AeroProjMod per rotor
   LOGICAL,                 INTENT(INOUT) :: calcCrvAngle(:) ! Whether this blade definition should calculate BlCrvAng (each blades and each rotor)
   INTEGER(IntKi),          INTENT(OUT)   :: ErrStat         ! The error status code
   CHARACTER(*),            INTENT(OUT)   :: ErrMsg          ! The error message, if an error occurred

      ! local variables

   INTEGER(IntKi)                         :: I
   INTEGER(IntKi)                         :: iR              ! Loop on rotor
   integer(IntKi)                         :: iBld            ! counter on blades
   INTEGER(IntKi)                         :: ErrStat2        ! The error status code
   CHARACTER(ErrMsgLen)                   :: ErrMsg2         ! The error message, if an error occurred

   CHARACTER(*), PARAMETER                :: RoutineName = 'ReadInputFiles'
   
   
      ! initialize values:

   ErrStat = ErrID_None
   ErrMsg  = ''

   InputFileData%DTAero = Default_DT  ! the glue code's suggested DT for the module (may be overwritten in ReadPrimaryFile())
   calcCrvAngle = .false.             ! initialize in case of early return

      ! get the blade input-file data
   iBld=1
   do iR = 1, size(InputFileData%rotors)
      
      ALLOCATE( InputFileData%rotors(iR)%BladeProps( NumBlades(iR) ), STAT = ErrStat2 )
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating memory for BladeProps.", ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
         RETURN
      END IF
         
   !FIXME: add options for passing the blade files.  This routine will need restructuring to handle that.
      DO I=1,NumBlades(iR)
         CALL ReadBladeInputs ( InputFileData%ADBlFile(iBld), InputFileData%rotors(iR)%BladeProps(I), AeroProjMod(iR), UnEcho, calcCrvAngle(iBld), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName//TRIM(':Blade')//TRIM(Num2LStr(I)))
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL Cleanup()
               RETURN
            END IF
         iBld = iBld+1 ! Increment blade counter
      END DO
   
   end do ! loop on rotors

   ! Read TailFin
   do iR = 1, size(InputFileData%rotors)
      if (InputFileData%rotors(iR)%TFinAero) then 
         call ReadTailFinInputs(InputFileData%rotors(iR)%TFinFile, InputFileData%rotors(iR)%TFin, UnEcho, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call Cleanup()
            return
         end if
      endif
   enddo ! iR, rotors
      

   CALL Cleanup ( )


CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up before exiting this subroutine
   !...............................................................................................................................

      ! IF ( UnEcho > 0 ) CLOSE( UnEcho )

   END SUBROUTINE Cleanup

END SUBROUTINE ReadInputFiles
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine parses the input file data stored in FileInfo_In and places it in the InputFileData structure for validating.
SUBROUTINE ParsePrimaryFileInfo( PriPath, InitInp, InputFile, RootName, NumBlades, interval, FileInfo_In, InputFileData, UnEc, ErrStat, ErrMsg )

   implicit    none

      ! Passed variables
   character(*),                    intent(in   )  :: PriPath           !< primary path
   type(AD_InitInputType),          intent(in   )  :: InitInp           !< Input data for initialization routine
   CHARACTER(*),                    intent(in   )  :: InputFile         !< Name of the file containing the primary input data
   CHARACTER(*),                    intent(in   )  :: RootName          !< The rootname of the echo file, possibly opened in this routine
   integer(IntKi),                  intent(in   )  :: NumBlades(:)      !< Number of blades per rotor we expect -- from InitInp
   real(DBKi),                      intent(in   )  :: interval          !< timestep
   type(AD_InputFile),              intent(inout)  :: InputFileData     !< All the data in the AD15 primary input file
   type(FileInfoType),              intent(in   )  :: FileInfo_In       !< The derived type for holding the file information.
   integer(IntKi),                  intent(  out)  :: UnEc              !< The local unit number for this module's echo file
   integer(IntKi),                  intent(  out)  :: ErrStat           !< Error status
   CHARACTER(ErrMsgLen),            intent(  out)  :: ErrMsg            !< Error message

      ! Local variables:
   integer(IntKi)                                  :: i                 !< generic counter
   integer(IntKi)                                  :: iR                !< Loop on rotors
   integer(IntKi)                                  :: numBladesTot      !< total number of blades
   integer(IntKi)                                  :: ErrStat2          !< Temporary Error status
   character(ErrMsgLen)                            :: ErrMsg2           !< Temporary Error message
   character(ErrMsgLen)                            :: ErrMsg_NoAllBldNdOuts
   integer(IntKi)                                  :: CurLine           !< current entry in FileInfo_In%Lines array
   real(ReKi)                                      :: TmpRe5(5)         !< temporary 8 number array for reading values in
   logical                                         :: TwrAeroLogical    !< convert TwrAero from logical (input file) to integer (new)
   character(1024)                                 :: sDummy            !< temporary string
   character(1024)                                 :: tmpOutStr         !< temporary string for writing to screen
   logical :: wakeModProvided, frozenWakeProvided, skewModProvided, AFAeroModProvided, UAModProvided, isLegalComment, firstWarn !< Temporary for legacy purposes
   logical :: AoA34_Missing
   integer :: UAMod_Old
   integer :: WakeMod_Old
   integer :: AFAeroMod_Old
   integer :: SkewMod_Old
   logical :: FrozenWake_Old
   character(*), parameter                         :: RoutineName = 'ParsePrimaryFileInfo'
   UAMod_Old      = -1
   WakeMod_Old    = -1
   AFAeroMod_Old  = -1
   SkewMod_Old    = -1
   FrozenWake_Old = .False.

   InputFileData%UA_Init%UA_OUTS    = 0
   InputFileData%UA_Init%d_34_to_ac = 0.5_ReKi
   

   ! Initialization
   ErrStat  =  ErrId_None
   ErrMsg   =  ""
   UnEc   = -1     ! Echo file unit.  >0 when used
   firstWarn=.False.


   CALL AllocAry( InputFileData%OutList, MaxOutPts, "Outlist", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Allocate array for holding the list of node outputs
   CALL AllocAry( InputFileData%BldNd_OutList, 2*BldNd_MaxOutPts, "BldNd_Outlist", ErrStat2, ErrMsg2 ) ! allow users to enter twice the number of unique outputs
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   numBladesTot=sum(NumBlades)

   !-------------------------------------------------------------------------------------------------
   ! General settings
   !-------------------------------------------------------------------------------------------------

   CurLine = 4    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2 )
         if (Failed()) return;

   if ( InputFileData%Echo ) then
      CALL OpenEcho ( UnEc, TRIM(RootName)//'.ech', ErrStat2, ErrMsg2 )
         if (Failed()) return;
      WRITE(UnEc, '(A)') 'Echo file for AeroDyn 15 primary input file: '//trim(InputFile)
      ! Write the first three lines into the echo file
      WRITE(UnEc, '(A)') FileInfo_In%Lines(1)
      WRITE(UnEc, '(A)') FileInfo_In%Lines(2)
      WRITE(UnEc, '(A)') FileInfo_In%Lines(3)

      CurLine = 4
      call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
   endif


      ! DTAero - Time interval for aerodynamic calculations {or default} (s):
   call ParseVarWDefault ( FileInfo_In, CurLine, "DTAero", InputFileData%DTAero, interval, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   ! WakeMod - LEGACY 
   call ParseVar( FileInfo_In, CurLine, "WakeMod", WakeMod_Old, ErrStat2, ErrMsg2, UnEc)
   wakeModProvided = legacyInputPresent('WakeMod', CurLine, ErrStat2, ErrMsg2, 'Wake_Mod=0 (WakeMod=0), Wake_Mod=1 (WakeMod=1), DBEMT_Mod>0 (WakeMod=2), Wake_Mod=3 (WakeMod=3)')
   ! Wake_Mod- Type of wake/induction model (switch) {0=none, 1=BEMT, 2=TBD, 3=OLAF}
   call ParseVar( FileInfo_In, CurLine, "Wake_Mod", InputFileData%Wake_Mod, ErrStat2, ErrMsg2, UnEc )
   if (newInputMissing('Wake_Mod', CurLine, errStat2, errMsg2)) then
      call WrScr('         Setting Wake_Mod to 1 (BEM active) as the input is Missing (typical behavior).')
      InputFileData%Wake_Mod = WakeMod_BEMT
   else
      if (wakeModProvided) then
         call LegacyAbort('Cannot have both Wake_Mod and WakeMod in the input file'); return
      endif
   endif


   ! AFAeroMod - Type of blade airfoil aerodynamics model (switch) {1=steady model, 2=Beddoes-Leishman unsteady model} [AFAeroMod must be 1 when linearizing]
   call ParseVar( FileInfo_In, CurLine, "AFAeroMod", AFAeroMod_Old, ErrStat2, ErrMsg2, UnEc )
   AFAeroModProvided = legacyInputPresent('AFAeroMod', CurLine, ErrStat2, ErrMsg2, 'UA_Mod=0 (AFAeroMod=1) or UA_Mod>1 (AFAeroMod=2)')
      ! TwrPotent - Type of tower influence on wind based on potential flow around the tower (switch) {0=none, 1=baseline potential flow, 2=potential flow with Bak correction}
   call ParseVar( FileInfo_In, CurLine, "TwrPotent", InputFileData%TwrPotent, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! TwrShadow - Type of tower influence on wind based on downstream tower shadow {0=none, 1=Powles model, 2=Eames model}
   call ParseVar( FileInfo_In, CurLine, "TwrShadow", InputFileData%TwrShadow, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! TwrAero - Calculate tower aerodynamic loads? (flag)
   call ParseVar( FileInfo_In, CurLine, "TwrAero", TwrAeroLogical, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      if (TwrAeroLogical) then
         InputFileData%TwrAero = TwrAero_NoVIV
      else
         InputFileData%TwrAero = TwrAero_None
      end if
      
   ! FrozenWake - Assume frozen wake during linearization? (flag) [used only when WakeMod=1 and when linearizing]
   call ParseVar( FileInfo_In, CurLine, "FrozenWake", FrozenWake_Old, ErrStat2, ErrMsg2, UnEc )
   frozenWakeProvided = legacyInputPresent('FrozenWake', Curline, ErrStat2, ErrMsg2, 'DBEMTMod=-1 (FrozenWake=True) or DBEMTMod>-1 (FrozenWake=False)')
      ! CavitCheck - Perform cavitation check? (flag) [AFAeroMod must be 1 when CavitCheck=true]
   call ParseVar( FileInfo_In, CurLine, "CavitCheck", InputFileData%CavitCheck, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! Buoyancy - Include buoyancy effects? (flag)
   call ParseVar( FileInfo_In, CurLine, "Buoyancy", InputFileData%Buoyancy, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! NacelleDrag - Include Nacelle Drag effects? (flag)
   call ParseVar( FileInfo_In, CurLine, "NacelleDrag", InputFileData%NacelleDrag, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! CompAA - Flag to compute AeroAcoustics calculation [only used when WakeMod=1 or 2]
   call ParseVar( FileInfo_In, CurLine, "CompAA", InputFileData%CompAA, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! AA_InputFile - Aeroacoustics input file
   call ParseVar( FileInfo_In, CurLine, "AA_InputFile", InputFileData%AA_InputFile, ErrStat2, ErrMsg2, UnEc, IsPath=.true. )
      if (Failed()) return
      IF ( PathIsRelative( InputFileData%AA_InputFile ) ) InputFileData%AA_InputFile = TRIM(PriPath)//TRIM(InputFileData%AA_InputFile)

   !======  Environmental Conditions  ===================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! AirDens - Air density {or default} (kg/m^3)
   call ParseVarWDefault( FileInfo_In, CurLine, "AirDens", InputFileData%AirDens, InitInp%defFldDens, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! KinVisc - Kinematic air viscosity {or default} (m^2/s)
   call ParseVarWDefault( FileInfo_In, CurLine, "KinVisc", InputFileData%KinVisc, InitInp%defKinVisc, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! SpdSound - Speed of sound {or default} (m/s)
   call ParseVarWDefault( FileInfo_In, CurLine, "SpdSound", InputFileData%SpdSound, InitInp%defSpdSound, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      InputFileData%UA_Init%a_s        = InputFileData%SpdSound

      ! Patm - Atmospheric pressure {or default} (Pa) [used only when CavitCheck=True]
   call ParseVarWDefault( FileInfo_In, CurLine, "Patm", InputFileData%Patm, InitInp%defPatm, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! Pvap - Vapour pressure of fluid {or default} (Pa) [used only when CavitCheck=True]
   call ParseVarWDefault( FileInfo_In, CurLine, "Pvap", InputFileData%Pvap, InitInp%defPvap, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

   !======  Blade-Element/Momentum Theory Options  ====================================================== [unused when WakeMod=0 or 3]
   call ParseCom (FileInfo_in, CurLine, sDummy, errStat2, errMsg2, UnEc, isLegalComment); if (Failed()) return

   ! BEM_Mod
   call ParseVar( FileInfo_In, CurLine, "BEM_Mod", InputFileData%BEM_Mod, ErrStat2, ErrMsg2, UnEc )
   if (newInputMissing('BEM_Mod', CurLine, errStat2, errMsg2)) then
      call WrScr('         Setting BEM_Mod to 1 (NoPitchSweepPitch) as the input is Missing (legacy behavior).')
      InputFileData%BEM_Mod = BEMMod_2D
   else
      call ParseCom (FileInfo_in, CurLine, sDummy, errStat2, errMsg2, UnEc, isLegalComment); if (Failed()) return
   endif

   ! SkewMod Legacy
   call ParseVar( FileInfo_In, CurLine, "SkewMod", SkewMod_Old, ErrStat2, ErrMsg2, UnEc )
   skewModProvided = legacyInputPresent('SkewMod', CurLine, ErrStat2, ErrMsg2, 'Skew_Mod=-1 (SkewMod=0), Skew_Mod=0 (SkewMod=1), Skew_Mod=1 (SkewMod>=2)')
   ! Skew_Mod-  Select skew model {0: No skew model at all, -1:Throw away non-normal component for linearization, 1: Glauert skew model, }
   call ParseVar( FileInfo_In, CurLine, "Skew_Mod", InputFileData%Skew_Mod, ErrStat2, ErrMsg2, UnEc )
   if (newInputMissing('Skew_Mod', CurLine, errStat2, errMsg2)) then
      call WrScr('         Setting Skew_Mod to 1 (skew active) as the input is Missing (typical behavior).')
      InputFileData%Skew_Mod = Skew_Mod_Active
   else
      if (skewModProvided) then
         call LegacyAbort('Cannot have both Skew_Mod and SkewMod in the input file'); return
      endif
   endif


   ! SkewMomCorr - Turn the skew momentum correction on or off [used only when SkewMod=1]
   call ParseVar( FileInfo_In, CurLine, "SkewMomCorr", InputFileData%SkewMomCorr, ErrStat2, ErrMsg2, UnEc )
   if (newInputMissing('SkewMomCorr', CurLine, errStat2, errMsg2)) then
      call WrScr('         Setting SkewMomCorr to False as the input is Missing (legacy behavior).')
      InputFileData%SkewMomCorr = .False.
   endif

   ! SkewRedistr_Mod - Type of skewed-wake correction model (switch) {0: no redistribution, 1=Glauert/Pitt/Peters, 2=Vortex Cylinder} [unsed only when SkewMod=1]
   call ParseVarWDefault( FileInfo_In, CurLine, "SkewRedistr_Mod", InputFileData%SkewRedistr_Mod, 1, ErrStat2, ErrMsg2, UnEc )
   if (newInputMissing('SkewRedistr_Mod', CurLine, errStat2, errMsg2)) then
      call WrScr('         Setting SkewRedistr_Mod to 1 as the input is Missing (legacy behavior).')
      InputFileData%SkewRedistr_Mod = 1
   endif

      ! SkewModFactor - Constant used in Pitt/Peters skewed wake model {or "default" is 15/32*pi} (-) [used only when SkewMod=2; unused when WakeMod=0 or 3]
   call ParseVarWDefault( FileInfo_In, CurLine, "SkewModFactor", InputFileData%SkewModFactor, (15.0_ReKi * pi / 32.0_ReKi), ErrStat2, ErrMsg2, UnEc )
   if( legacyInputPresent('SkewModFactor', CurLine, ErrStat2, ErrMsg2, 'Rename this parameter to SkewRedistrFactor')) then
      ! pass
   else
      call ParseVarWDefault( FileInfo_In, CurLine, "SkewRedistrFactor", InputFileData%SkewModFactor, (15.0_ReKi * pi / 32.0_ReKi), ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
      call ParseCom (FileInfo_in, CurLine, sDummy, errStat2, errMsg2, UnEc, isLegalComment); if (Failed()) return
   endif

      ! TipLoss - Use the Prandtl tip-loss model? (flag) [unused when WakeMod=0 or 3]
   call ParseVar( FileInfo_In, CurLine, "TipLoss", InputFileData%TipLoss, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! HubLoss - Use the Prandtl hub-loss model? (flag) [unused when WakeMod=0 or 3]
   call ParseVar( FileInfo_In, CurLine, "HubLoss", InputFileData%HubLoss, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! TanInd - Include tangential induction in BEMT calculations? (flag) [unused when WakeMod=0 or 3]
   call ParseVar( FileInfo_In, CurLine, "TanInd", InputFileData%TanInd, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! AIDrag - Include the drag term in the axial-induction calculation? (flag) [unused when WakeMod=0 or 3]
   call ParseVar( FileInfo_In, CurLine, "AIDrag", InputFileData%AIDrag, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! TIDrag - Include the drag term in the tangential-induction calculation? (flag) [unused when WakeMod=0,3 or TanInd=FALSE]
   call ParseVar( FileInfo_In, CurLine, "TIDrag", InputFileData%TIDrag, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! IndToler - Convergence tolerance for BEMT nonlinear solve residual equation {or "default"} (-) [unused when WakeMod=0 or 3]
   if (ReKi==SiKi) then
      call ParseVarWDefault( FileInfo_In, CurLine, "IndToler", InputFileData%IndToler, real(5E-5,ReKi), ErrStat2, ErrMsg2, UnEc )
   else
      call ParseVarWDefault( FileInfo_In, CurLine, "IndToler", InputFileData%IndToler, real(5D-10,ReKi), ErrStat2, ErrMsg2, UnEc )
   end if
      if (Failed()) return
      ! MaxIter - Maximum number of iteration steps (-) [unused when WakeMod=0]
   call ParseVar( FileInfo_In, CurLine, "MaxIter", InputFileData%MaxIter, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   ! ---  Shear
   call ParseCom (FileInfo_in, CurLine, sDummy, errStat2, errMsg2, UnEc, isLegalComment); if (Failed()) return
   call ParseVar( FileInfo_In, CurLine, "SectAvg"         , InputFileData%SectAvg, ErrStat2, ErrMsg2, UnEc ); 
   if (newInputMissing('SectAvg', CurLine, errStat2, errMsg2)) then
      call WrScr('         Setting SectAvg to False as the input is Missing (legacy behavior).')
      InputFileData%SectAvg = .false.
   else
      call ParseVarWDefault( FileInfo_In, CurLine, "SectAvgWeighting", InputFileData%SA_Weighting, 1      , ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
      call ParseVarWDefault( FileInfo_In, CurLine, "SectAvgNPoints"  , InputFileData%SA_nPerSec,   5      , ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
      call ParseVarWDefault( FileInfo_In, CurLine, "SectAvgPsiBwd"   , InputFileData%SA_PsiBwd,  -60._ReKi, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
      call ParseVarWDefault( FileInfo_In, CurLine, "SectAvgPsiFwd"   , InputFileData%SA_PsiFwd,   60._ReKi, ErrStat2, ErrMsg2, UnEc ); if (Failed()) return
      call ParseCom (FileInfo_in, CurLine, sDummy, errStat2, errMsg2, UnEc, isLegalComment); if (Failed()) return
   endif

   !======  Dynamic Blade-Element/Momentum Theory Options  ============================================== [used only when WakeMod=2]
   ! DBEMT_Mod - Type of dynamic BEMT (DBEMT) model {1=constant tau1, 2=time-dependent tau1} (-) [used only when WakeMod=2]
   call ParseVar( FileInfo_In, CurLine, "DBEMT_Mod", InputFileData%DBEMT_Mod, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! tau1_const - Time constant for DBEMT (s) [used only when WakeMod=2 and DBEMT_Mod=1]
   call ParseVar( FileInfo_In, CurLine, "tau1_const", InputFileData%tau1_const, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

   !======  OLAF -- cOnvecting LAgrangian Filaments (Free Vortex Wake) Theory Options  ================== [used only when WakeMod=3]
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! OLAFInputFileName - Input file for OLAF [used only when WakeMod=3]
   call ParseVar( FileInfo_In, CurLine, "OLAFInputFileName", InputFileData%FVWFileName, ErrStat2, ErrMsg2, UnEc, IsPath=.true. )
      if (Failed()) return
      IF ( PathIsRelative( InputFileData%FVWFileName ) ) InputFileData%FVWFileName = TRIM(PriPath)//TRIM(InputFileData%FVWFileName)

   !======  Beddoes-Leishman Unsteady Airfoil Aerodynamics Options  ===================================== [used only when AFAeroMod=2]
   call ParseCom (FileInfo_in, CurLine, sDummy, errStat2, errMsg2, UnEc, isLegalComment); if (Failed()) return
   ! AoA34 Sample the angle of attack (AoA) at the 3/4 chord or the AC point {default=True} [always used]
   call ParseVar( FileInfo_In, CurLine, "AoA34", InputFileData%AoA34, ErrStat2, ErrMsg2, UnEc )
   AoA34_Missing = newInputMissing('AoA34', CurLine, errStat2, errMsg2)
   ! UAMod (Legacy)
   call ParseVar( FileInfo_In, CurLine, "UAMod", UAMod_Old, ErrStat2, ErrMsg2, UnEc )
   UAModProvided = legacyInputPresent('UAMod', CurLine, ErrStat2, ErrMsg2, 'UA_Mod=0 (AFAeroMod=1), UA_Mod>1 (AFAeroMod=2 and UA_Mod=UAMod')
   ! UA_Mod - Unsteady Aero Model Switch (switch) {0=Quasi-steady (no UA),  2=Gonzalez's variant (changes in Cn,Cc,Cm), 3=Minnema/Pierce variant (changes in Cc and Cm)} 
   call ParseVar( FileInfo_In, CurLine, "UA_Mod", InputFileData%UA_Init%UAMod, ErrStat2, ErrMsg2, UnEc )
   if (newInputMissing('UA_Mod', CurLine, errStat2, errMsg2)) then
      ! We'll deal with it when we deal with AFAeroMod
      InputFileData%UA_Init%UAMod = UAMod_Old
      if (.not. UAModProvided) then
         call LegacyAbort('Need to provide either UA_Mod or UAMod in the input file'); return
      endif
   else
      if (UAModProvided) then
         call LegacyAbort('Cannot have both UA_Mod and UAMod in the input file'); return
      endif
   endif


      ! FLookup - Flag to indicate whether a lookup for f' will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files (flag) [used only when AFAeroMod=2]
   call ParseVar( FileInfo_In, CurLine, "FLookup", InputFileData%UA_Init%FLookup, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! IntegrationMethod - Switch to indicate which integration method UA uses (1=RK4, 2=AB4, 3=ABM4, 4=BDF2) (switch):
   call ParseVar( FileInfo_In, CurLine, "IntegrationMethod", InputFileData%UA_Init%IntegrationMethod, ErrStat2, ErrMsg2, UnEc )
      if (ErrStat2>= AbortErrLev) InputFileData%UA_Init%IntegrationMethod = UA_Method_ABM4
   
      ! UAStartRad - Starting radius for dynamic stall (fraction of rotor radius) [used only when AFAeroMod=2]:
   call ParseVar( FileInfo_In, CurLine, "UAStartRad", InputFileData%UAStartRad, ErrStat2, ErrMsg2, UnEc )
      if (ErrStat2>= AbortErrLev) InputFileData%UAStartRad = 0.0_ReKi
      ! UAEndRad - Ending radius for dynamic stall (fraction of rotor radius) [used only when AFAeroMod=2]:
   call ParseVar( FileInfo_In, CurLine, "UAEndRad", InputFileData%UAEndRad, ErrStat2, ErrMsg2, UnEc )
      if (ErrStat2>= AbortErrLev) InputFileData%UAEndRad = 1.0_ReKi

   !======  Airfoil Information =========================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! AFTabMod - Interpolation method for multiple airfoil tables {1=1D interpolation on AoA (first table only); 2=2D interpolation on AoA and Re; 3=2D interpolation on AoA and UserProp} (-)
   call ParseVar( FileInfo_In, CurLine, "AFTabMod", InputFileData%AFTabMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! InCol_Alfa - The column in the airfoil tables that contains the angle of attack (-)
   call ParseVar( FileInfo_In, CurLine, "InCol_Alfa", InputFileData%InCol_Alfa, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! InCol_Cl - The column in the airfoil tables that contains the lift coefficient (-)
   call ParseVar( FileInfo_In, CurLine, "InCol_Cl", InputFileData%InCol_Cl, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! InCol_Cd - The column in the airfoil tables that contains the drag coefficient (-)
   call ParseVar( FileInfo_In, CurLine, "InCol_Cd", InputFileData%InCol_Cd, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! InCol_Cm - The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)
   call ParseVar( FileInfo_In, CurLine, "InCol_Cm", InputFileData%InCol_Cm, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! InCol_Cpmin - The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column (-)
   call ParseVar( FileInfo_In, CurLine, "InCol_Cpmin", InputFileData%InCol_Cpmin, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! NumAFfiles - Number of airfoil files used (-)
   call ParseVar( FileInfo_In, CurLine, "NumAFfiles", InputFileData%NumAFfiles, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
         ! Allocate space to hold AFNames
      ALLOCATE( InputFileData%AFNames(InputFileData%NumAFfiles), STAT=ErrStat2)
         IF (ErrStat2 /= 0 ) THEN
            ErrStat2=ErrID_Fatal
            ErrMsg2 = "Error allocating AFNames."
            if (Failed()) return
         END IF
      ! AFNames - Airfoil file names (NumAFfiles lines) (quoted strings): -- NOTE: this line may not have a keyname with it
   DO I = 1,InputFileData%NumAFfiles         ! ParseChVar allows empty keynames.
      call ParseVar( FileInfo_In, CurLine, "", InputFileData%AFNames(I), ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      IF ( PathIsRelative( InputFileData%AFNames(I) ) ) InputFileData%AFNames(I) = TRIM(PriPath)//TRIM(InputFileData%AFNames(I))
   END DO

   !======  Rotor/Blade Properties  =====================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! UseBlCm - Include aerodynamic pitching moment in calculations?  (flag)
   call ParseVar( FileInfo_In, CurLine, "UseBlCm", InputFileData%UseBlCm, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! Allocate space for AD blade file names -- MaxBl is usually set to 3, but if we specify more blades, this will work still.
   call AllocAry( InputFileData%ADBlFile, max(MaxBl,NumBladesTot), 'ADBlFile', ErrStat2, ErrMsg2)
      if (Failed()) return
   do I =1,size(InputFileData%ADBlFile)  ! We expect MaxBl blade file lines.  We may want to revisit this idea later if we allow more thn 3 blades
      call ParseVar( FileInfo_In, CurLine, "", InputFileData%ADBlFile(i), ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
      IF ( PathIsRelative( InputFileData%ADBlFile(I) ) ) InputFileData%ADBlFile(I) = TRIM(PriPath)//TRIM(InputFileData%ADBlFile(I))
   enddo

   !======  Hub Properties ============================================================================== [used only when Buoyancy=True]
   do iR = 1,size(NumBlades) ! Loop on rotors
      if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
      CurLine = CurLine + 1
         ! VolHub - Hub volume (m^3)
      call ParseVar( FileInfo_In, CurLine, "VolHub", InputFileData%rotors(iR)%VolHub, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return   
         ! HubCenBx - Hub center of buoyancy x direction offset (m)
      call ParseVar( FileInfo_In, CurLine, "HubCenBx", InputFileData%rotors(iR)%HubCenBx, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
   end do

   !======  Nacelle Properties ========================================================================== [used only when Buoyancy=True]
   do iR = 1,size(NumBlades) ! Loop on rotors
      if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
      CurLine = CurLine + 1
         ! VolNac - Nacelle volume (m^3)
      call ParseVar( FileInfo_In, CurLine, "VolNac", InputFileData%rotors(iR)%VolNac, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
         ! NacCenB - Nacelle center of buoyancy x,y,z direction offsets (m)
      call ParseAry( FileInfo_In, CurLine, 'NacCenB', InputFileData%rotors(iR)%NacCenB, 3 , ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return

         ! NacArea - Projected area of the nacelle in X, Y, Z in the nacelle coordinate system (m^2)
      call ParseAry( FileInfo_In, CurLine, "NacArea", InputFileData%rotors(iR)%NacArea, 3, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
         ! NacCd - Drag coefficient for the nacelle areas defined above (-)
         call ParseAry( FileInfo_In, CurLine, "NacCd", InputFileData%rotors(iR)%NacCd, 3, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
         ! NacDragAC - Position of aerodynamic center of nacelle drag in nacelle coordinates (m)
         call ParseAry( FileInfo_In, CurLine, "NacDragAC", InputFileData%rotors(iR)%NacDragAC, 3, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
   end do

   !======  Tail fin aerodynamics ========================================================================
   do iR = 1,size(NumBlades) ! Loop on rotors
      if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
      CurLine = CurLine + 1
      ! NOTE: being nice with legacy input file. Uncomment in next release
      call ParseVar(FileInfo_In, CurLine, "TFinAero", InputFileData%rotors(iR)%TFinAero, ErrStat2, ErrMsg2, UnEc); 
      if (ErrStat2==ErrID_None) then
         call ParseVar(FileInfo_In, CurLine, "TFinFile", InputFileData%rotors(iR)%TFinFile, ErrStat2, ErrMsg2, UnEc, IsPath=.true.); if (Failed()) return
         IF ( PathIsRelative( InputFileData%rotors(iR)%TFinFile ) ) InputFileData%rotors(iR)%TFinFile = trim(PriPath) // trim(InputFileData%rotors(iR)%TFinFile)
      else
         call LegacyWarning('Tail Fin section (TFinAero, TFinFile) is missing from input file.')
         CurLine = CurLine - 1
      endif
   enddo

   !======  Tower Influence and Aerodynamics ============================================================ [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or Buoyancy=True]
   do iR = 1,size(NumBlades) ! Loop on rotors
      if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
      CurLine = CurLine + 1
         ! NumTwrNds - Number of tower nodes used in the analysis  (-) [used only when TwrPotent/=0, TwrShadow/=0, TwrAero=True, or Buoyancy=True]
      call ParseVar( FileInfo_In, CurLine, "NumTwrNds", InputFileData%rotors(iR)%NumTwrNds, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
         !TwrElev        TwrDiam        TwrCd        TwrTI        TwrCb
      if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Tower Table Header: '//FileInfo_In%Lines(CurLine)    ! Write section break to echo
      CurLine = CurLine + 1
         !(m)            (m)            (-)          (-)          (-)
      if ( InputFileData%Echo )   WRITE(UnEc, '(A)') 'Tower Table Header: '//FileInfo_In%Lines(CurLine)    ! Write section break to echo
      CurLine = CurLine + 1
         ! Allocate space for tower table
      CALL AllocAry( InputFileData%rotors(iR)%TwrElev, InputFileData%rotors(iR)%NumTwrNds, 'TwrElev',  ErrStat2, ErrMsg2)
         if (Failed()) return
      CALL AllocAry( InputFileData%rotors(iR)%TwrDiam, InputFileData%rotors(iR)%NumTwrNds, 'TwrDiam', ErrStat2, ErrMsg2)
         if (Failed()) return
      CALL AllocAry( InputFileData%rotors(iR)%TwrCd, InputFileData%rotors(iR)%NumTwrNds, 'TwrCd', ErrStat2, ErrMsg2)
         if (Failed()) return
      CALL AllocAry( InputFileData%rotors(iR)%TwrTI, InputFileData%rotors(iR)%NumTwrNds, 'TwrTI', ErrStat2, ErrMsg2)
         if (Failed()) return
      CALL AllocAry( InputFileData%rotors(iR)%TwrCb, InputFileData%rotors(iR)%NumTwrNds, 'TwrCb', ErrStat2, ErrMsg2)
         if (Failed()) return 
      do I=1,InputFileData%rotors(iR)%NumTwrNds
         call ParseAry ( FileInfo_In, CurLine, 'Properties for tower node '//trim( Int2LStr( I ) )//'.', TmpRe5, 5, ErrStat2, ErrMsg2, UnEc )
            if (Failed()) return;
         InputFileData%rotors(iR)%TwrElev(I) = TmpRe5( 1)
         InputFileData%rotors(iR)%TwrDiam(I) = TmpRe5( 2)
         InputFileData%rotors(iR)%TwrCd(I)   = TmpRe5( 3)
         InputFileData%rotors(iR)%TwrTI(I)   = TmpRe5( 4)
         InputFileData%rotors(iR)%TwrCb(I)   = TmpRe5( 5)
      end do
   enddo

   !======  Outputs  ====================================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! SumPrint - Generate a summary file listing input options and interpolated properties to "<rootname>.AD.sum"?  (flag)
   call ParseVar( FileInfo_In, CurLine, "SumPrint", InputFileData%SumPrint, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   InputFileData%UA_Init%WrSum      = InputFileData%SumPrint

      ! NBlOuts - Number of blade node outputs [0 - 9] (-)
   call ParseVar( FileInfo_In, CurLine, "NBlOuts", InputFileData%NBlOuts, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! Make sure we don't try to read in more than will fit in the pre-allocated BlOutNd array
   if ( InputFileData%NBlOuts > SIZE(InputFileData%BlOutNd) ) THEN
      CALL SetErrStat( ErrID_Warn, ' Warning: number of blade output nodes exceeds '//&
                        TRIM(Num2LStr(SIZE(InputFileData%BlOutNd))) //'.', ErrStat, ErrMsg, RoutineName )
      InputFileData%NBlOuts = SIZE(InputFileData%BlOutNd)
   endif
      ! BlOutNd - Blade nodes whose values will be output (-):
   call ParseAry( FileInfo_In, CurLine, "BlOutNd", InputFileData%BlOutNd, InputFileData%NBlOuts, ErrStat2, ErrMsg2, UnEc)
      if (Failed()) return

      ! NTwOuts - Number of blade node outputs [0 - 9] (-)
   call ParseVar( FileInfo_In, CurLine, "NTwOuts", InputFileData%NTwOuts, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! Make sure we don't try to read in more than will fit in the pre-allocated TwOutNd array
   if ( InputFileData%NTwOuts > SIZE(InputFileData%TwOutNd) ) THEN
      CALL SetErrStat( ErrID_Warn, ' Warning: number of blade output nodes exceeds '//&
                        TRIM(Num2LStr(SIZE(InputFileData%TwOutNd))) //'.', ErrStat, ErrMsg, RoutineName )
      InputFileData%NTwOuts = SIZE(InputFileData%TwOutNd)
   endif
      ! TwOutNd - Tower nodes whose values will be output (-):
   call ParseAry( FileInfo_In, CurLine, "TwOutNd", InputFileData%TwOutNd, InputFileData%NTwOuts, ErrStat2, ErrMsg2, UnEc)
      if (Failed()) return

   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
   call ReadOutputListFromFileInfo( FileInfo_In, CurLine, InputFileData%OutList, InputFileData%NumOuts, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return;

   !====== Legacy logic to match old and new input files ================================================
   ! NOTE: remove me in future release
   if (frozenWakeProvided) then
      if (FrozenWake_Old) then
         call WrScr('> FrozenWake=True  -> Setting DBEMT_Mod=-1')
         InputFileData%DBEMT_Mod = DBEMT_frozen
      else
         call WrScr('> FrozenWake=False -> Not changing DBEMT_Mod')
      endif
   endif
   if (wakeModProvided) then
      InputFileData%Wake_Mod = WakeMod_Old
      if (WakeMod_Old==1) then
         call WrScr('> WakeMod=1        -> Setting DBEMT_Mod=0')
         ! Turn off DBEMT
         InputFileData%DBEMT_Mod=DBEMT_none
      else if (WakeMod_Old==2) then
         call WrScr('> WakeMod=2        -> Setting Wake_Mod=1 (BEMT) (DBEMT_Mod needs to be >0)')
         InputFileData%Wake_Mod = WakeMod_BEMT
         if (InputFileData%DBEMT_Mod < DBEMT_none) then
            call LegacyAbort('DBEMT should be >0 when using legacy input WakeMod=2'); return
         endif
      endif
   endif
   if (AFAeroModProvided) then
      if (AFAeroMod_Old==1) then
         call WrScr('> AFAeroMod=1      -> Setting UA_Mod=0')
         InputFileData%UA_Init%UAMod = UA_None
         if (AoA34_Missing) then
            call WrScr('> Setting AoA34 to False as the input is Missing and UA is turned off (legacy behavior).')
            InputFileData%AoA34=.false.
         endif
      else if (AFAeroMod_Old==2) then
         call WrScr('> AFAeroMod=2      -> Not changing DBEMT_Mod')
         if (InputFileData%UA_Init%UAMod==0) then
            call LegacyAbort('Cannot set UA_Mod=0 with legacy option AFAeroMod=2 (inconsistent behavior).'); return
         else if (AoA34_Missing) then
            call WrScr('> Setting AoA34 to True as the input is Missing and UA is turned on (legacy behavior).')
            InputFileData%AoA34=.true.
         endif
      else
         call LegacyAbort('AFAeroMod should be 1 or 2'); return
      endif
   endif
   if (skewModProvided) then
      if (SkewMod_Old==0) then
         InputFileData%Skew_Mod = Skew_Mod_Orthogonal
      else if (SkewMod_Old==1) then
         InputFileData%Skew_Mod = Skew_Mod_None
      else if (SkewMod_Old==2) then
         InputFileData%Skew_Mod = Skew_Mod_Active
      else
         call LegacyAbort('Legacy option SkewMod is not 0, 1,2  which is not supported.'); return
      endif
   endif

   !====== Print new and old inputs =====================================================================
   if (wakeModProvided .or. frozenWakeProvided .or. skewModProvided .or. AFAeroModProvided .or. UAModProvided) then
      call printNewOldInputs()
   endif

   !======  Nodal Outputs  ==============================================================================
      ! In case there is something ill-formed in the additional nodal outputs section, we will simply ignore it.
      ! Expecting at least 5 more lines in the input file for this section
   if (FileInfo_In%NumLines < CurLine + 5) then
      ErrStat2 = ErrID_Fatal
      ErrMsg2  = ''
      if (FailedNodal()) return;
   endif

   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      ! BldNd_BladesOut - Number of blades to output all node information at.  Up to number of blades on turbine. (-)
      ! TODO:  In a future release, allow this to be an array of N blade numbers (change BldNd_BladesOut to an array if we do that).
      !        Will likely require reading this line in as a string (BldNd_BladesOut_Str) and parsing it
   call ParseVar( FileInfo_In, CurLine, "BldNd_BladesOut", InputFileData%BldNd_BladesOut, ErrStat2, ErrMsg2, UnEc )
      if (FailedNodal()) return
      ! BldNd_BlOutNd - Allow selecting a portion of the nodes to output. (-)
   call ParseVar( FileInfo_In, CurLine, "BldNd_BlOutNd", InputFileData%BldNd_BlOutNd_Str, ErrStat2, ErrMsg2, UnEc )
   if (ErrStat2 /= ErrID_None) then
      ! ParseVar won't read a string of numbers in quotes since the quotes are a delimiter, so we'll just copy the whole line here and move on 
      InputFileData%BldNd_BlOutNd_Str = FileInfo_In%Lines(CurLine)
      if ( InputFileData%Echo )   WRITE(UnEc, '(A)') InputFileData%BldNd_BlOutNd_Str    ! Write BldNd_BlOutNd_Str to echo

      CurLine = CurLine + 1
      ErrStat2 = ErrID_None
   end if
      
      ! OutList - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   call ReadOutputListFromFileInfo( FileInfo_In, CurLine, InputFileData%BldNd_OutList, InputFileData%BldNd_NumOuts, ErrStat2, ErrMsg2, UnEc )
         if (FailedNodal()) return;

!FIXME: improve logic on the node outputs
   ! Prevent segfault when no blades specified.  All logic tests on BldNd_NumOuts at present.
   if (InputFileData%BldNd_BladesOut <= 0)   InputFileData%BldNd_NumOuts = 0


   !====== Advanced Options =============================================================================
   if ((CurLine) >= size(FileInfo_In%Lines)) RETURN

   call WrScr(' - Reading advanced options for AeroDyn')
   do CurLine= CurLine, size(FileInfo_In%Lines)
      sDummy = FileInfo_In%Lines(CurLine)
      call Conv2UC(sDummy)  ! to uppercase
      if (index(sDummy, '!') == 1 .or. index(sDummy, '=') == 1 .or. index(sDummy, '#') == 1 .or. index(sDummy, '---') == 1) then
         ! pass comment lines
         elseif (index(sDummy, 'SECTAVG')>1) then
            read(sDummy, '(L1)') InputFileData%SectAvg
            write(tmpOutStr,*) '   >>> SectAvg        ',InputFileData%SectAvg
         elseif (index(sDummy, 'SA_PSIBWD')>1) then
            read(sDummy, *) InputFileData%SA_PsiBwd
            write(tmpOutStr,*) '   >>> SA_PsiBwd      ',InputFileData%SA_PsiBwd
         elseif (index(sDummy, 'SA_PSIFWD')>1) then
            read(sDummy, *) InputFileData%SA_PsiFwd
            write(tmpOutStr,*) '   >>> SA_PsiFwd      ',InputFileData%SA_PsiFwd
         elseif (index(sDummy, 'SA_NPERSEC')>1) then
            read(sDummy, *) InputFileData%SA_nPerSec
            write(tmpOutStr,*) '   >>> SA_nPerSec     ',InputFileData%SA_nPerSec
         else
            write(tmpOutStr,*) '[WARN] AeroDyn Line ignored: '//trim(sDummy)
      endif
      call WrScr(trim(tmpOutStr))
   enddo


   !---------------------- END OF FILE -----------------------------------------
   

   RETURN


CONTAINS
   !-------------------------------------------------------------------------------------------------
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ParsePrimaryFileInfo' )
      Failed = ErrStat >= AbortErrLev
      !if (Failed) then
      !endif
   end function Failed
   logical function FailedNodal()
      ErrMsg_NoAllBldNdOuts='AD15 Nodal Outputs: Nodal output section of AeroDyn input file not found or improperly formatted. Skipping nodal outputs.'
      ! TODO Use and ErrID_Fatal here
      FailedNodal = ErrStat2 >= AbortErrLev 
      if ( FailedNodal ) then
         InputFileData%BldNd_BladesOut = 0
         InputFileData%BldNd_NumOuts = 0
         call wrscr( trim(ErrMsg_NoAllBldNdOuts) )
      endif
   end function FailedNodal
   subroutine LegacyWarning(Message)
      character(len=*), intent(in) :: Message
      if (.not.FirstWarn) then
         call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
         call WrScr('[WARN] The AeroDyn input file is not at the latest format!' )
         call WrScr('       Visit: https://openfast.readthedocs.io/en/dev/source/user/api_change.html')
         call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
         FirstWarn=.true.
      endif
      call WrScr('> Issue: '//trim(Message))
   end subroutine LegacyWarning
   !-------------------------------------------------------------------------------------------------
   subroutine LegacyAbort(Message)
      character(len=*), intent(in) :: Message
      call SetErrStat( ErrID_Fatal, Message, ErrStat, ErrMsg, 'ParsePrimaryFileInfo' )
   end subroutine LegacyAbort
   !-------------------------------------------------------------------------------------------------
   logical function legacyInputPresent(varName, iLine, errStat, errMsg, varNameSubs)
      character(len=*),           intent(in  )  :: varName     !< Variable being read
      integer(IntKi),             intent(in   ) :: iLine       !< Line number
      integer(IntKi),             intent(inout) :: errStat     !< Error status
      character(ErrMsgLen),       intent(inout) :: errMsg      !< Error message
      character(len=*), optional, intent(in  )  :: varNameSubs !< Substituted variable
      legacyInputPresent = errStat == ErrID_None
      if (legacyInputPresent) then
         if (present(varNameSubs)) then
            call LegacyWarning(trim(varName)//' has now been removed.'//NewLine//'    Use: '//trim(varNameSubs)//'.')
         else
            call LegacyWarning(trim(varName)//' has now been removed.')
         endif
      else
         ! We are actually happy, this input should indeed not be present.
      endif
      ! We erase the error no matter what
      errStat = ErrID_None
      errMsg  = ''
   end function legacyInputPresent
   !-------------------------------------------------------------------------------------------------
   logical function newInputMissing(varName, iLine, errStat, errMsg, varNameSubs)
      character(len=*),           intent(in  )  :: varName     !< Variable being read
      integer(IntKi),             intent(in   ) :: iLine       !< Line number
      integer(IntKi),             intent(inout) :: errStat     !< Error status
      character(ErrMsgLen),       intent(inout) :: errMsg      !< Error message
      character(len=*), optional, intent(in  )  :: varNameSubs !< Substituted variable
      newInputMissing = errStat == ErrID_Fatal
      if (newInputMissing) then
         call LegacyWarning(trim(varName)//' should be present on line '//trim(num2lstr(iLine))//'.')
      else
         ! We are happy
      endif
      ! We erase the error
      errStat = ErrID_None
      errMsg  = ''
   end function newInputMissing

   !-------------------------------------------------------------------------------------------------
   subroutine printNewOldInputs()
      character(1024)   :: tmpStr
      ! Temporary HACK, for WakeMod=10, 11 or 12 use AeroProjMod 2 (will trigger PolarBEM)
      if (InputFileData%Wake_Mod==10) then
         call WrScr('[WARN] Wake_Mod=10 is a temporary hack. Setting BEM_Mod to 0')
         InputFileData%BEM_Mod = 0
      elseif (InputFileData%Wake_Mod==11) then
         call WrScr('[WARN] Wake_Mod=11 is a temporary hack. Setting BEM_Mod to 2')
         InputFileData%BEM_Mod = 2
      elseif (InputFileData%Wake_Mod==12) then
         call WrScr('[WARN] Wake_Mod=12 is a temporary hack. Setting BEM_Mod to 2')
         InputFileData%BEM_Mod = 2
      endif
      !====== Summary of new AeroDyn options ===============================================================
      ! NOTE: remove me in future release
      call WrScr('-------------- New AeroDyn inputs (with new meaning):')
      write (tmpStr,'(A20,I0)') 'Wake_Mod: '         , InputFileData%Wake_Mod;         call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'BEM_Mod:  '         , InputFileData%BEM_Mod;          call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,L1)') 'SectAvg:  '         , InputFileData%SectAvg;          call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'SectAvgWeighting:  ', InputFileData%SA_Weighting;     call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'SectAvgNPoints:    ', InputFileData%SA_nPerSec;       call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'DBEMT_Mod:'         , InputFileData%DBEMT_Mod;        call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'Skew_Mod:  '        , InputFileData%Skew_Mod;         call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,L1)') 'SkewMomCorr:'       , InputFileData%SkewMomCorr;      call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'SkewRedistr_Mod:'   , InputFileData%SkewRedistr_Mod;  call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,L1)') 'AoA34:    '         , InputFileData%AoA34;            call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'UA_Mod:   '         , InputFileData%UA_Init%UAMod;    call WrScr(trim(tmpStr))
      call WrScr('-------------- Old AeroDyn inputs:')
      write (tmpStr,'(A20,I0)') 'WakeMod:  ',  WakeMod_Old;      call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'SkewMod:  ',  SkewMod_Old;      call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'AFAeroMod:',  AFAeroMod_Old;    call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,L1)') 'FrozenWake:', FrozenWake_Old;   call WrScr(trim(tmpStr))
      write (tmpStr,'(A20,I0)') 'UAMod:     ', UAMod_Old;        call WrScr(trim(tmpStr))
      call WrScr('------------------------------------------------------')
   end subroutine printNewOldInputs

END SUBROUTINE ParsePrimaryFileInfo
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadBladeInputs ( ADBlFile, BladeKInputFileData, AeroProjMod, UnEc, calcCrvAngle, ErrStat, ErrMsg )
! This routine reads a blade input file.
!..................................................................................................................................


      ! Passed variables:

   TYPE(AD_BladePropsType),  INTENT(INOUT)  :: BladeKInputFileData                 ! Data for Blade K stored in the module's input file
   CHARACTER(*),             INTENT(IN)     :: ADBlFile                            ! Name of the blade input file data
   INTEGER(IntKi),           INTENT(IN)     :: AeroProjMod                         ! AeroProjMod
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                ! I/O unit for echo file. If present and > 0, write to UnEc
   LOGICAL,                  INTENT(INOUT)  :: calcCrvAngle                        ! Whether this blade definition should calculate BlCrvAng

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message


      ! Local variables:

   INTEGER(IntKi)                            :: I                                               ! A generic DO index.
   INTEGER( IntKi )                          :: UnIn                                            ! Unit number for reading file
   INTEGER(IntKi)                            :: ErrStat2 , IOS                                  ! Temporary Error status
   CHARACTER(ErrMsgLen)                      :: ErrMsg2                                         ! Temporary Err msg
   INTEGER,         PARAMETER                :: MaxCols = 10
   CHARACTER(NWTC_SizeOfNumWord*(MaxCols+1)) :: Line
   INTEGER(IntKi)                            :: Indx(MaxCols)
   CHARACTER(8),   PARAMETER                 :: AvailableChanNames(MaxCols) = (/'BLSPN   ', 'BLCRVAC ','BLSWPAC ','BLCRVANG','BLTWIST ','BLCHORD ', 'BLAFID  ', 'BLCB    ', 'BLCENBN ','BLCENBT ' /) ! in upper case only
   LOGICAL,        PARAMETER                 :: RequiredChanNames( MaxCols) = (/.true.    , .true.    ,.true.    ,.false.   ,.true.    ,.true.    , .true.    , .false.   , .false.   ,.false.    /)

   CHARACTER(*), PARAMETER                   :: RoutineName = 'ReadBladeInputs'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   UnIn = -1
   
   ! Open the input file for blade K.
   !$OMP critical(fileopen_critical)
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
   CALL OpenFInpFile ( UnIn, ADBlFile, ErrStat2, ErrMsg2 )
   !$OMP end critical(fileopen_critical)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- HEADER -------------------------------------------------------

      ! Skip the header.

   CALL ReadCom ( UnIn, ADBlFile, 'unused blade file header line 1', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom ( UnIn, ADBlFile, 'unused blade file header line 2', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   !  -------------- Blade properties table ------------------------------------------                                    
   CALL ReadCom ( UnIn, ADBlFile, 'Section header: Blade Properties', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! NumBlNds - Number of blade nodes used in the analysis (-):
   CALL ReadVar( UnIn, ADBlFile, BladeKInputFileData%NumBlNds, "NumBlNds", "Number of blade nodes used in the analysis (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat>= AbortErrLev ) THEN 
         CALL Cleanup()
         RETURN
      END IF

   CALL ReadCom ( UnIn, ADBlFile, 'Table header: names', ErrStat2, ErrMsg2, UnEc, Comment=Line )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom ( UnIn, ADBlFile, 'Table header: units', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   IF ( ErrStat>= AbortErrLev ) THEN 
      CALL Cleanup()
      RETURN
   END IF
   
      
      ! allocate space for blade inputs:
   CALL AllocAry( BladeKInputFileData%BlSpn,   BladeKInputFileData%NumBlNds, 'BlSpn',   ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlCrvAC, BladeKInputFileData%NumBlNds, 'BlCrvAC', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlSwpAC, BladeKInputFileData%NumBlNds, 'BlSwpAC', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlCrvAng,BladeKInputFileData%NumBlNds, 'BlCrvAng',ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlTwist, BladeKInputFileData%NumBlNds, 'BlTwist', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlChord, BladeKInputFileData%NumBlNds, 'BlChord', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlAFID,  BladeKInputFileData%NumBlNds, 'BlAFID',  ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlCb, BladeKInputFileData%NumBlNds, 'BlCb', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlCenBn, BladeKInputFileData%NumBlNds, 'BlCenBn', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlCenBt, BladeKInputFileData%NumBlNds, 'BlCenBt', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Return on error if we didn't allocate space for the next inputs
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
   
   ! Initialize in case these columns are missing (e.g., no buoyancy, or cant angle)
   BladeKInputFileData%BlCrvAng = 0.0_ReKi
   BladeKInputFileData%BlCb     = 0.0_ReKi
   BladeKInputFileData%BlCenBn  = 0.0_ReKi
   BladeKInputFileData%BlCenBt  = 0.0_ReKi
   
   
   ! figure out what columns are specified in this file and in what order:
   CALL GetInputColumnIndex(MaxCols, AvailableChanNames, RequiredChanNames, Line, Indx, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


   DO I=1,BladeKInputFileData%NumBlNds
      CALL ReadCom( UnIn, ADBlFile, 'Blade properties row '//TRIM(Num2LStr(I)), ErrStat2, ErrMsg2, UnEc, Comment=Line ) ! this will get echoed as a comment instead of a table
      ! Return on error if we couldn't read this line
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
         
      CALL ConvertLineToCols(Line, i, Indx, BladeKInputFileData, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   END DO

   BladeKInputFileData%BlTwist  = BladeKInputFileData%BlTwist*D2R
   BladeKInputFileData%BlCrvAng = BladeKInputFileData%BlCrvAng*D2R
   
   ! note, we will compute BlCrvAng later for APM_BEM_Polar or in the case the BlCrvAng column is missing from the file
   calcCrvAngle = AeroProjMod==APM_BEM_Polar .or. Indx(4) < 1
   
   if (Indx(4) > 0 .and. calcCrvAngle) then
      CALL SetErrStat(ErrID_Warn,'BlCrvAng will be calculated and overwrite the values specified in blade file "'//trim(ADBlFile)//'".', ErrStat, ErrMsg, RoutineName)
   end if
   
   !bjj: do we still need this???
   if (all(BladeKInputFileData%BlCrvAC.eq.0.0_ReKi)) then
      BladeKInputFileData%BlCrvAng = 0.0_ReKi
   endif
   
   !  -------------- END OF FILE --------------------------------------------

   CALL Cleanup()
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up local variables and closes files
   !...............................................................................................................................

      IF (UnIn > 0) CLOSE(UnIn)

   END SUBROUTINE Cleanup
   !...............................................................................................................................
END SUBROUTINE ReadBladeInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ConvertLineToCols(Line, i, Indx, BladeKInputFileData, ErrStat, ErrMsg)
   CHARACTER(*),             INTENT(IN   )   :: Line                                ! text containing line we are reading/parsing
   INTEGER(IntKi),           INTENT(IN   )   :: i                                   ! row of input table we are reading
   INTEGER(IntKi),           INTENT(IN   )   :: Indx(:)                             ! order of table columns, determined from headers in subroutine GetInputColumnIndex()
   TYPE(AD_BladePropsType),  INTENT(INOUT)   :: BladeKInputFileData                 ! Data for Blade K stored in the module's input file
   INTEGER(IntKi),           INTENT(  OUT)   :: ErrStat
   CHARACTER(*),             INTENT(  OUT)   :: ErrMsg
      
   INTEGER(IntKi)                            :: NumCols                             ! number of columns in the file
   CHARACTER(NWTC_SizeOfNumWord)             :: Words(size(Indx))
   INTEGER                                   :: IOS(size(Indx))
   INTEGER                                   :: c                                   ! column index
   CHARACTER(*), PARAMETER                   :: RoutineName = 'ConvertLineToCols'

      
   ErrStat = ErrID_None
   ErrMsg  = ""
   
!   CALL GetWords ( Line, Words, size(Indx), NumCols )
   !IF (MAXVAL(Indx) > NumCols) THEN
   !   CALL SetErrStat(ErrID_Fatal, "Required column is not available in the table on row "//trim(num2lstr(i))//".", ErrStat, ErrMsg, RoutineName)
   !   RETURN
   !END IF
   
   ! Read "words" as character strings:
   NumCols = MAXVAL(Indx)
   READ(Line, *, IOStat=IOS(1)) Words(1:NumCols)
      
   IOS = 0 ! initialize in case we don't read all of the columns
      
   ! Note: See order of variable AvailableChanNames in subroutine ReadBladeInputs() for these variables indices
   ! Also, we have checked that Indx is non zero and less than MaxCols for each of the required words
   c=Indx( 1); READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlSpn(I)
   c=Indx( 2); READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlCrvAC(I)
   c=Indx( 3); READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlSwpAC(I)
   c=Indx( 5); READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlTwist(I)
   c=Indx( 6); READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlChord(I)
   c=Indx( 7); READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlAFID(I)

   c=Indx(4)
   IF (c > 0) THEN
      READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlCrvAng(I)
   END IF

   c=Indx(8)
   IF (c > 0) THEN
      READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlCb(I)
   END IF

   c=Indx(9)
   IF (c > 0) THEN
      READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlCenBn(I)
   END IF

   c=Indx(10)
   IF (c > 0) THEN
      READ( Words(c), *, IOStat=IOS(c) ) BladeKInputFileData%BlCenBt(I)
   END IF

   IF (ANY(IOS /= 0)) THEN
      CALL SetErrStat(ErrID_Fatal, "Unable to read numeric data from all columns in the table on row "//trim(num2lstr(i))//".", ErrStat, ErrMsg, RoutineName)
      RETURN
   END IF
      
END SUBROUTINE ConvertLineToCols

!----------------------------------------------------------------------------------------------------------------------------------
!> Read Tail Fin inputs
SUBROUTINE ReadTailFinInputs(FileName, TFData, UnEc, ErrStat, ErrMsg)
   character(*),                    intent(in   )  :: FileName          !< Name of the file containing the Tail Fin aero data
   type(TFinInputFileType),         intent(inout)  :: TFData            !< All the data in the Tail Fin input file
   integer(IntKi),                  intent(in   )  :: UnEc              !< Echo unit 
   integer(IntKi),                  intent(  out)  :: ErrStat           !< Error status
   character(ErrMsgLen),            intent(  out)  :: ErrMsg            !< Error message
   ! Local
   type(FileInfoType)   :: FileInfo_In ! < The derived type for holding the file information.
   integer(IntKi)       :: iLine       !< current entry in FileInfo_In%Lines array
   integer(IntKi)       :: ErrStat2    !< Temporary Error status
   character(ErrMsgLen) :: ErrMsg2     !< Temporary Error message
   character(len=1024 ) :: DummyLine

   ! --- Read Tail fin input file into array of strings
   call ProcessComFile( FileName, FileInfo_In, ErrStat2, ErrMsg2)

   ! --- Parse the array of strings
   ! Skip the first two lines as they are known to be header lines and separators
   do iLine = 1,2 
      if ( UnEc>0 )   WRITE(UnEc, '(A)') FileInfo_In%Lines(iLine)    ! Write header to echo
   enddo
   iLine = 3 
   !====== General inputs ============================================================
   call ParseCom(FileInfo_in, iLine, DummyLine                          , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinMod'   , TFData%TFinMod       , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinArea'  , TFData%TFinArea      , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_In, iLine, 'TFinRefP_n', TFData%TFinRefP_n, 3 , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_In, iLine, 'TFinAngles', TFData%TFinAngles, 3 , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinIndMod', TFData%TFinIndMod    , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   !====== Polar-based model ================================ [used only when TFinMod=1]
   call ParseCom(FileInfo_in, iLine, DummyLine                          , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinAFID'  , TFData%TFinAFID      , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinChord' , TFData%TFinChord     , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   !====== Unsteady slender body model ===================== [used only when TFinMod=2]
   call ParseCom(FileInfo_in, iLine, DummyLine                          , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinKp'  , TFData%TFinKp      , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_In, iLine, 'TFinSigma'  , TFData%TFinSigma, 3 , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_In, iLine, 'TFinAStar', TFData%TFinAStar, 3 , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinKv'    , TFData%TFinKv       , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinCDc'   , TFData%TFinCDc       , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   
   ! TODO

   ! --- Triggers
   TFData%TFinAngles = TFData%TFinAngles*D2R ! deg2rad

   ! --- Validation on the fly
   if (all((/TFinAero_none,TFinAero_polar,TFinAero_USB/) /= TFData%TFinMod)) then
      call Fatal('TFinMod needs to be 0, 1 or 2')
   endif
   !if (all((/TFinIndMod_none,TFinIndMod_rotavg/) /= TFData%TFinIndMod)) then
   if (all((/TFinIndMod_none/) /= TFData%TFinIndMod)) then
      call Fatal('TFinIndMod needs to be 0')
   endif

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadTailFinInputs' )
      Failed = ErrStat >= AbortErrLev
   end function Failed

   subroutine Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      call SetErrStat(ErrID_Fatal, 'File:'//trim(FileName)//':'//trim(ErrMsg_in), ErrStat, ErrMsg, 'ReadTailFinInputs')
   end subroutine Fatal
   
END SUBROUTINE ReadTailFinInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_PrintSum( InputFileData, p, p_AD, u, y, NumBlades, BladeInputFileData, ErrStat, ErrMsg )
! This routine generates the summary file, which contains a summary of input file options.
   use YAML, only: yaml_write_var
      ! passed variables
   TYPE(AD_InputFile),        INTENT(IN)  :: InputFileData                        ! Input-file data
   TYPE(RotParameterType),    INTENT(IN)  :: p                                    ! Parameters
   TYPE(AD_ParameterType),    INTENT(IN)  :: p_AD                                 ! Parameters
   TYPE(AD_InputType),        INTENT(IN)  :: u                                    ! inputs 
   TYPE(AD_OutputType),       INTENT(IN)  :: y                                    ! outputs
   INTEGER(IntKi),            INTENT(IN)  :: NumBlades                            ! Number of blades for this rotor
   TYPE(AD_BladePropsType),   INTENT(IN)  :: BladeInputFileData(:)                ! Data for Bladex stored in the module's input file
   
   INTEGER(IntKi),            INTENT(OUT) :: ErrStat
   CHARACTER(*),              INTENT(OUT) :: ErrMsg


      ! Local variables.

   INTEGER(IntKi)               :: I                                               ! Index for the nodes.
   INTEGER(IntKi)               :: K                                               ! Index for the blades
   INTEGER(IntKi)               :: UnSu                                            ! I/O unit number for the summary output file

   CHARACTER(*), PARAMETER      :: FmtDat    = '(A,T41,1(:,F13.3))'                ! Format for outputting mass and modal data.
   CHARACTER(*), PARAMETER      :: FmtDatT   = '(A,T35,1(:,F13.8))'                ! Format for outputting time steps.

   CHARACTER(30)                :: OutPFmt                                         ! Format to print list of selected output channels to summary file
   CHARACTER(30)                :: OutPFmtS                                        ! Format to print list of selected output channels to summary file
   CHARACTER(100)               :: Msg                                             ! temporary string for writing appropriate text to summary file

   CHARACTER(ChanLen),PARAMETER :: TitleStr(2) = (/ 'Parameter', 'Units    ' /)
   CHARACTER(ChanLen),PARAMETER :: TitleStrLines(2) = (/ '---------------', '---------------' /)

   ! Open the summary file and give it a heading.
      
   !$OMP critical(fileopen_critical)
   CALL GetNewUnit( UnSu, ErrStat, ErrMsg )
   CALL OpenFOutFile ( UnSu, TRIM( p%RootName )//'.sum', ErrStat, ErrMsg )
   !$OMP end critical(fileopen_critical)
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Heading:
   WRITE (UnSu,'(/,A)')  'This summary information was generated by '//TRIM( GetNVD(AD_Ver) )// &
                         ' on '//CurDate()//' at '//CurTime()//'.'

   WRITE (UnSu,'(/,A)') '======  General Options  ============================================================================'
   ! WakeMod
   select case (p_AD%Wake_Mod)
      case (WakeMod_BEMT)
         Msg = 'Blade-Element/Momentum Theory'
      case (WakeMod_FVW)
         Msg = 'Free Vortex Wake Theory'
      case (WakeMod_None)
         Msg = 'steady'
      case default      
         Msg = 'unknown'      
   end select   
   WRITE (UnSu,Ec_IntFrmt) p_AD%Wake_Mod, 'WakeMod', 'Type of wake/induction model: '//TRIM(Msg)
   
   ! TwrPotent
   select case (p%TwrPotent)
      case (TwrPotent_baseline)
         Msg = 'baseline potential flow'
      case (TwrPotent_Bak)
         Msg = 'potential flow with Bak correction'
      case (TwrPotent_none)
         Msg = 'none'      
      case default      
         Msg = 'unknown'      
   end select
   WRITE (UnSu,Ec_IntFrmt) p%TwrPotent, 'TwrPotent', 'Type of tower influence on wind based on potential flow around the tower: '//TRIM(Msg)
   
   
   ! TwrShadow
   select case (p%TwrShadow)
      case (TwrShadow_Powles)
         Msg = 'Powles tower shadow model'
      case (TwrShadow_Eames)
         Msg = 'Eames tower shadow model with TI values from the table' 
      case (TwrShadow_none)
         Msg = 'none'
      case default
         Msg = 'unknown'
   end select
   WRITE (UnSu,Ec_IntFrmt) p%TwrShadow, 'TwrShadow', 'Type of tower influence on wind based on downstream tower shadow: '//TRIM(Msg)
   
   
   ! TwrAero
   select case (p%TwrAero)
      case (TwrAero_none)
         Msg = "none"
      case (TwrAero_NoVIV)
         Msg = "Tower aero calculated without VIV"
!      case (TwrAero_VIV)
!         Msg = "Tower aero calculated with VIV"
      case default
         Msg = 'unknown'
   end select
   WRITE (UnSu,Ec_LgFrmt) p%TwrAero, 'TwrAero', 'Tower aerodynamic loads: '//TRIM(Msg)

      ! CavitCheck
   if (p%CavitCheck) then
      Msg = 'Yes'
   else
      Msg = 'No'
   end if   
   WRITE (UnSu,Ec_LgFrmt) p%CavitCheck, 'CavitCheck', 'Perform cavitation check? '//TRIM(Msg)

      ! Buoyancy
   if (p%Buoyancy) then
      Msg = 'Yes'
   else
      Msg = 'No'
   end if   
   WRITE (UnSu,Ec_LgFrmt) p%Buoyancy, 'Buoyancy', 'Include buoyancy effects? '//TRIM(Msg)

      ! Nacelle Drag
   if (p%NacelleDrag) then
      Msg = 'Yes'
   else
      Msg = 'No'
   end if   
   WRITE (UnSu,Ec_LgFrmt) p%NacelleDrag, 'NacelleDrag', 'Include NacelleDrag effects? '//TRIM(Msg)


   if (p_AD%Wake_Mod/=WakeMod_none) then
      WRITE (UnSu,'(A)') '======  Blade-Element/Momentum Theory Options  ======================================================'
      
      ! SkewMod 
      select case (InputFileData%Skew_Mod)
         case (Skew_Mod_Orthogonal)
            Msg = 'orthogonal'
         case (Skew_Mod_None)
            Msg = 'no correction'
         case (Skew_Mod_Active)
            Msg = 'active'
         case default      
            Msg = 'unknown'      
      end select
      WRITE (UnSu,Ec_IntFrmt) InputFileData%Skew_Mod, 'Skew_Mod', 'Skewed-wake correction model: '//TRIM(Msg)
      
      
      ! TipLoss
      if (InputFileData%TipLoss) then
         Msg = 'Yes'
      else
         Msg = 'No'
      end if   
      WRITE (UnSu,Ec_LgFrmt) InputFileData%TipLoss, 'TipLoss', "Use the Prandtl tip-loss model? "//TRIM(Msg)      
      
      
      ! HubLoss 
      if (InputFileData%HubLoss) then
         Msg = 'Yes'
      else
         Msg = 'No'
      end if   
      WRITE (UnSu,Ec_LgFrmt) InputFileData%HubLoss, 'HubLoss', "Use the Prandtl hub-loss model? "//TRIM(Msg)      

      
      ! TanInd  
      if (InputFileData%TanInd) then
         Msg = 'Yes'
      else
         Msg = 'No'
      end if   
      WRITE (UnSu,Ec_LgFrmt) InputFileData%TanInd, 'TanInd', "Include tangential induction in BEMT calculations? "//TRIM(Msg)      

      
      ! AIDrag 
      if (InputFileData%AIDrag) then
         Msg = 'Yes'
      else
         Msg = 'No'
      end if   
      WRITE (UnSu,Ec_LgFrmt) InputFileData%AIDrag, 'AIDrag', "Include the drag term in the axial-induction calculation? "//TRIM(Msg)      
      
      ! TIDrag  
      if (InputFileData%TIDrag .and. InputFileData%TanInd) then
         Msg = 'Yes'
      else
         Msg = 'No'
      end if   
      WRITE (UnSu,Ec_LgFrmt) InputFileData%TIDrag, 'TIDrag', "Include the drag term in the tangential-induction calculation? "//TRIM(Msg)      
      
      ! IndToler
      WRITE (UnSu,Ec_ReFrmt) InputFileData%IndToler, 'IndToler', "Convergence tolerance for BEM induction factors (radians)"     
      
      ! MaxIter 
      
      
      select case (InputFileData%DBEMT_Mod)
         case (DBEMT_frozen)
            Msg = 'frozen-wake'
         case (DBEMT_none)
            Msg = 'quasi-steady'
         case (DBEMT_tauConst)
            Msg = 'dynamic - constant tau1'
         case (DBEMT_tauVaries)
            Msg = 'dynamic - time-dependent tau1'
         case (DBEMT_cont_tauConst)
            Msg = 'dynamic - continuous formulation with constant tau1'
         case default
            Msg = 'unknown'
      end select   
      
      WRITE (UnSu,Ec_IntFrmt) InputFileData%DBEMT_Mod, 'DBEMT_Mod', 'Type of dynamic BEMT (DBEMT) model: '//TRIM(Msg)
         
      if (InputFileData%DBEMT_Mod==DBEMT_tauConst) &
      WRITE (UnSu,Ec_ReFrmt) InputFileData%tau1_const, 'tau1_const', 'Time constant for DBEMT (s)'
         
      
   end if
   
   WRITE (UnSu,'(A)') '======================== Unsteady Airfoil Aerodynamics Options  ====================================='
   
   ! UAMod
   select case (InputFileData%UA_Init%UAMod)
      case (UA_None)
         Msg = 'none (quasi-steady airfoil aerodynamics)'
      case (UA_Baseline)
         Msg = 'baseline model (original)'
      case (UA_Gonzalez)
         Msg = "Gonzalez's variant (changes in Cn, Cc, and Cm)"
      case (UA_MinnemaPierce)
         Msg = 'Minnema/Pierce variant (changes in Cc and Cm)'      
      !case (4)
      !   Msg = 'DYSTOOL'      
      case (UA_HGM)
         Msg = 'HGM (continuous state)'
      case (UA_HGMV)
         Msg = 'HGMV (continuous state + vortex)'
      case (UA_OYE)
         Msg = 'Stieg Oye dynamic stall model'
      case (UA_BV)
         Msg = 'Boeing-Vertol dynamic stall model (e.g. used in CACTUS)'
      case default
         Msg = 'unknown'
   end select
   WRITE (UnSu,Ec_IntFrmt) InputFileData%UA_Init%UAMod, 'UA_Mod', 'Unsteady Aero Model: '//TRIM(Msg)


   ! FLookup
   if (InputFileData%UA_Init%FLookup) then
      Msg = 'Yes'
   else
      Msg = 'No, use best-fit exponential equations instead'
   end if   
   WRITE (UnSu,Ec_LgFrmt) InputFileData%UA_Init%FLookup, 'FLookup', "Use a lookup for f'? "//TRIM(Msg)


   ! IntegrationMethod
   select case (InputFileData%UA_Init%IntegrationMethod)
      case (UA_Method_RK4)
         Msg = 'fourth-order Runge-Kutta Method (RK4)'
      case (UA_Method_AB4)
         Msg = 'fourth-order Adams-Bashforth Method (AB4)'
      case (UA_Method_ABM4)
         Msg = "fourth-order Adams-Bashforth-Moulton Method (ABM4)"
      case (UA_Method_BDF2)
         Msg = '2nd-order backward differentiation formula (BDF2)'
      case default
         Msg = 'unknown'
   end select
   WRITE (UnSu,Ec_IntFrmt) InputFileData%UA_Init%IntegrationMethod, 'IntegrationMethod', 'Integration method for continuous UA models: '//TRIM(Msg)


   ! UAStartRad, UAEndRad
   WRITE (UnSu,"( 2X, F11.5,2X,A,T30,' - ',A )") InputFileData%UAStartRad, 'UAStartRad', 'Starting blade radius fraction for UA models (-)'  ! compare with Ec_ReFrmt format statement
   WRITE (UnSu,"( 2X, F11.5,2X,A,T30,' - ',A )") InputFileData%UAEndRad,    'UAEndRad', 'Ending blade radius fraction for UA models (-)'

   
   WRITE (UnSu,'(A)') '======  Outputs  ===================================================================================='
   
   OutPFmt = '( 49X, I11, 2X, I13 )'
   
   WRITE(UnSu,Ec_IntFrmt) p%NBlOuts,'NBlOuts','Number of blade nodes selected for output'
   if (p%NBlOuts > 0) then
      WRITE(UnSu,Ec_IntFrmt) p%NumBlNds,'NumBlNds','Number of blade nodes in the analysis'
      
      WRITE (UnSu,"(15x,A)")  'Blade nodes selected for output:  Output node  Analysis node'
      WRITE (UnSu,"(15x,A)")  '                                  -----------  -------------'
      DO I = 1,p%NBlOuts
         WRITE (UnSu,OutPFmt)  I, p%BlOutNd(I)
      END DO  
   end if
   
   WRITE(UnSu,Ec_IntFrmt) p%NTwOuts,'NTwOuts','Number of tower nodes selected for output'
   if (p%NTwOuts > 0) then
      WRITE(UnSu,Ec_IntFrmt) p%NumTwrNds,'NumTwrNds','Number of tower nodes in the analysis'
      WRITE (UnSu,"(15x,A)")  'Tower nodes selected for output:  Output node  Analysis node'
      WRITE (UnSu,"(15x,A)")  '                                  -----------  -------------'
      DO I = 1,p%NTwOuts
         WRITE (UnSu,OutPFmt)  I, p%TwOutNd(I)
      END DO  
   end if
   
   
   OutPFmt  = '( 15x, I4, 3X,A '//TRIM(Num2LStr(ChanLen))//',1 X, A'//TRIM(Num2LStr(ChanLen))//' )'
   OutPFmtS = '( 15x, A4, 3X,A '//TRIM(Num2LStr(ChanLen))//',1 X, A'//TRIM(Num2LStr(ChanLen))//' )'
   WRITE (UnSu,'(15x,A)')
   WRITE (UnSu,'(15x,A)')  'Requested Output Channels:'
   WRITE (UnSu,OutPFmtS )  "Col", TitleStr
   WRITE (UnSu,OutPFmtS )  "---", TitleStrLines
   DO I = 0,p%NumOuts
      WRITE (UnSu,OutPFmt)  I, p%OutParam(I)%Name, p%OutParam(I)%Units
   END DO             

   WRITE (UnSu,'(15x,A)')
   WRITE (UnSu,'(15x,A)')
   WRITE (UnSu,'(15x,A)')  'Requested Output Channels at each blade station:'
   WRITE (UnSu,OutPFmtS )  "Col", TitleStr
   WRITE (UnSu,OutPFmtS )  "---", TitleStrLines
   DO I = 1,p%BldNd_NumOuts
      WRITE (UnSu,OutPFmt)  I, p%BldNd_OutParam(I)%Name, p%BldNd_OutParam(I)%Units
   END DO

   ! Buoyancy parameters
   WRITE (UnSu,'(//,A,/)')  'Buoyancy parameters:'
   call yaml_write_var ( UnSu , 'Hub volume (m^3)' , p%VolHub , 'F15.3' , ErrStat , ErrMsg ) ! Buoyancy volume of the hub
   call yaml_write_var ( UnSu , 'Nacelle volume (m^3)' , p%VolNac , 'F11.3' , ErrStat , ErrMsg ) ! Buoyancy volume of the nacelle
   call yaml_write_var ( UnSu , 'Total blade volume (m^3)' , p%VolBl  , 'F7.3' , ErrStat , ErrMsg ) ! Buoyancy volume of all blades
   call yaml_write_var ( UnSu , 'Tower volume (m^3)' , p%VolTwr , 'F13.3' , ErrStat , ErrMsg ) ! Buoyancy volume of the tower

   WRITE (UnSu,'(/,/,A)') '======  Blade definitions  =========================================================================='
   
   DO k=1,NumBlades
      WRITE (UnSu,'(15x,A)')
      WRITE (UnSu,'(3x,A,I2,A)') '-----  Blade ', k, ' -----------------------------------------------------------------------------------------------------------'
      WRITE (UnSu,'(6(1x,A20))') 'BlSpn',   'BlCrvAC', 'BlSwpAC','BlCrvAng','BlTwist', 'BlChord'
      WRITE (UnSu,'(6(1x,A20))') '(m)',     '(m)',     '(m)',    '(deg)',   '(deg)',   '(m)'
      WRITE (UnSu,'(6(1x,A20))') '--------','--------','-------','--------','--------','--------'
      DO I=1,size(BladeInputFileData(K)%BlSpn)
         WRITE( UnSu, '(3( 1X,F20.6), 2( 1X,F20.4 ), 1( 1X,F20.6))') &
                               BladeInputFileData(K)%BlSpn(I),        BladeInputFileData(K)%BlCrvAC(I),      BladeInputFileData(K)%BlSwpAC(I), &
                               BladeInputFileData(K)%BlCrvAng(I)*R2D, BladeInputFileData(K)%BlTwist(I)*R2D,  BladeInputFileData(K)%BlChord(I)
      END DO
   END DO
   
   
   
   CLOSE(UnSu)

RETURN
END SUBROUTINE AD_PrintSum
!----------------------------------------------------------------------------------------------------------------------------------


!**********************************************************************************************************************************
! NOTE: The following lines of code were generated by a Matlab script called "Write_ChckOutLst.m"
!      using the parameters listed in the "OutListParameters.xlsx" Excel file. Any changes to these 
!      lines should be modified in the Matlab script and/or Excel worksheet as necessary. 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine checks to see if any requested output channel names (stored in the OutList(:)) are invalid. It returns a 
!! warning if any of the channels are not available outputs from the module.
!!  It assigns the settings for OutParam(:) (i.e, the index, name, and units of the output channels, WriteOutput(:)).
!!  the sign is set to 0 if the channel is invalid.
!! It sets assumes the value p%NumOuts has been set before this routine has been called, and it sets the values of p%OutParam here.
!! 
!! This routine was generated by Write_ChckOutLst.m using the parameters listed in OutListParameters.xlsx at 27-Oct-2022 11:00:28.
SUBROUTINE SetOutParam(OutList, p, p_AD, ErrStat, ErrMsg )
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables
   
   CHARACTER(ChanLen),        INTENT(IN)     :: OutList(:)                        !< The list of user-requested outputs
   TYPE(RotParameterType),    INTENT(INOUT)  :: p                                 !< The module parameters
   TYPE(AD_ParameterType),    INTENT(INOUT)  :: p_AD                              !< The module parameters
   INTEGER(IntKi),            INTENT(OUT)    :: ErrStat                           !< The error status code
   CHARACTER(*),              INTENT(OUT)    :: ErrMsg                            !< The error message, if an error occurred
   
      ! Local variables
   
   INTEGER                      :: ErrStat2                                        ! temporary (local) error status
   INTEGER                      :: I                                               ! Generic loop-counting index
   INTEGER                      :: J                                               ! Generic loop-counting index
   INTEGER                      :: INDX                                            ! Index for valid arrays
   LOGICAL                      :: InvalidOutput(0:MaxOutPts)                      ! This array determines if the output channel is valid for this configuration
   CHARACTER(*), PARAMETER      :: RoutineName = "SetOutParam"
   

      ! Initialize values
   ErrStat = ErrID_None
   ErrMsg = ""
   InvalidOutput = .FALSE.


!   ..... Developer must add checking for invalid inputs here: .....

   !bjj: do we want to avoid outputting this if we haven't used tower aero?
   
   if ( p%TwrPotent == TwrPotent_none .and. p%TwrShadow == TwrShadow_none ) then
      
         ! BNClrnc is set only when we're computing the tower influence
      do i = 1,size(BNClrnc,2)  ! all blades (need to do this in a loop because we need the index of InvalidOutput to be an array of rank one)
         InvalidOutput( BNClrnc(:,i) ) = .true.
      end do
      
   end if
      
   if (p%DBEMT_Mod == DBEMT_none) then
      InvalidOutput( DBEMTau1 ) = .true.
   end if

   if (.not. p%Buoyancy) then  ! Invalid buoyant loads
      InvalidOutput( HbFbx ) = .true.
      InvalidOutput( HbFby ) = .true.
      InvalidOutput( HbFbz ) = .true.
      InvalidOutput( HbMbx ) = .true.
      InvalidOutput( HbMby ) = .true.
      InvalidOutput( HbMbz ) = .true.
      InvalidOutput( NcFbx ) = .true.
      InvalidOutput( NcFby ) = .true.
      InvalidOutput( NcFbz ) = .true.
      InvalidOutput( NcMbx ) = .true.
      InvalidOutput( NcMby ) = .true.
      InvalidOutput( NcMbz ) = .true.
      InvalidOutput( TwNFbx ) = .true.
      InvalidOutput( TwNFby ) = .true.
      InvalidOutput( TwNFbz ) = .true.
      InvalidOutput( TwNMbx ) = .true.
      InvalidOutput( TwNMby ) = .true.
      InvalidOutput( TwNMbz ) = .true.
      do i = 1,size(BNFbn,2)
         InvalidOutput( BNFbn(:,i) ) = .true.
      end do
      do i = 1,size(BNFbt,2)
         InvalidOutput( BNFbt(:,i) ) = .true.
      end do
      do i = 1,size(BNFbs,2)
         InvalidOutput( BNFbs(:,i) ) = .true.
      end do
      do i = 1,size(BNMbn,2)
         InvalidOutput( BNMbn(:,i) ) = .true.
      end do
      do i = 1,size(BNMbt,2) 
         InvalidOutput( BNMbt(:,i) ) = .true.
      end do
      do i = 1,size(BNMbs,2)
         InvalidOutput( BNMbs(:,i) ) = .true.
      end do
   end if

   if (.not. p%NacelleDrag) then  ! Invalid Nacelle Drag loads
      InvalidOutput( NcFdx ) = .true.
      InvalidOutput( NcFdy ) = .true.
      InvalidOutput( NcFdz ) = .true.
      InvalidOutput( NcMdx ) = .true.
      InvalidOutput( NcMdy ) = .true.
      InvalidOutput( NcMdz ) = .true.   
   end if
   

   if (.not. (p%NacelleDrag .OR. p%Buoyancy)) then  ! Invalid Nacelle Total loads
      InvalidOutput( NcFxi ) = .true.
      InvalidOutput( NcFyi ) = .true.
      InvalidOutput( NcFzi ) = .true.
      InvalidOutput( NcMxi ) = .true.
      InvalidOutput( NcMyi ) = .true.
      InvalidOutput( NcMzi ) = .true.   
   end if


   DO i = p%NTwOuts+1,9  ! Invalid tower nodes
   
      InvalidOutput( TwNVUnd(:,i) ) = .true.
      InvalidOutput( TwNSTV( :,i) ) = .true.
      InvalidOutput( TwNVRel(  i) ) = .true.
      InvalidOutput( TwNDynP(  i) ) = .true.
      InvalidOutput( TwNRe(    i) ) = .true.
      InvalidOutput( TwNM(     i) ) = .true.
      InvalidOutput( TwNFdx(   i) ) = .true.
      InvalidOutput( TwNFdy(   i) ) = .true.
      InvalidOutput( TwNFbx(   i) ) = .true.
      InvalidOutput( TwNFby(   i) ) = .true.
      InvalidOutput( TwNFbz(   i) ) = .true.
      InvalidOutput( TwNMbx(   i) ) = .true.
      InvalidOutput( TwNMby(   i) ) = .true.
      InvalidOutput( TwNMbz(   i) ) = .true.
      
   END DO
   
   DO I = p%NumBlades+1,size(BAzimuth,1)  ! Invalid blades (Note: size(BAzimuth) should be AD_MaxBl_Out)
      
      InvalidOutput( BAzimuth( i) ) = .true.
      InvalidOutput( BPitch(   i) ) = .true.
      InvalidOutput( BNVUndx(:,i) ) = .true.
      InvalidOutput( BNVUndy(:,i) ) = .true.
      InvalidOutput( BNVUndz(:,i) ) = .true.
      InvalidOutput( BNVDisx(:,i) ) = .true.
      InvalidOutput( BNVDisy(:,i) ) = .true.
      InvalidOutput( BNVDisz(:,i) ) = .true.
      InvalidOutput( BNSTVx( :,i) ) = .true.
      InvalidOutput( BNSTVy( :,i) ) = .true.
      InvalidOutput( BNSTVz( :,i) ) = .true.
      InvalidOutput( BNVRel( :,i) ) = .true.
      InvalidOutput( BNDynP( :,i) ) = .true.
      InvalidOutput( BNRe(   :,i) ) = .true.
      InvalidOutput( BNM(    :,i) ) = .true.
      InvalidOutput( BNVIndx(:,i) ) = .true.
      InvalidOutput( BNVIndy(:,i) ) = .true.
      InvalidOutput( BNAxInd(:,i) ) = .true.
      InvalidOutput( BNTnInd(:,i) ) = .true.
      InvalidOutput( BNAlpha(:,i) ) = .true.
      InvalidOutput( BNTheta(:,i) ) = .true.
      InvalidOutput( BNPhi(  :,i) ) = .true.
      InvalidOutput( BNCurve(:,i) ) = .true.
      InvalidOutput( BNCl(   :,i) ) = .true.
      InvalidOutput( BNCd(   :,i) ) = .true.
      InvalidOutput( BNCm(   :,i) ) = .true.
      InvalidOutput( BNCx(   :,i) ) = .true.
      InvalidOutput( BNCy(   :,i) ) = .true.
      InvalidOutput( BNCn(   :,i) ) = .true.
      InvalidOutput( BNCt(   :,i) ) = .true.
      InvalidOutput( BNFl(   :,i) ) = .true.
      InvalidOutput( BNFd(   :,i) ) = .true.
      InvalidOutput( BNMm(   :,i) ) = .true.
      InvalidOutput( BNFx(   :,i) ) = .true.
      InvalidOutput( BNFy(   :,i) ) = .true.
      InvalidOutput( BNFn(   :,i) ) = .true.
      InvalidOutput( BNFt(   :,i) ) = .true.
      InvalidOutput( BNClrnc(:,i) ) = .true.
      InvalidOutput( BNGam(  :,i) ) = .true.
      InvalidOutput( BNSgCav(:,i) ) = .true.
      InvalidOutput( BNSigCr(:,i) ) = .true.
      InvalidOutput( BNCpMin(:,i) ) = .true.
      InvalidOutput( BAeroFx(  i) ) = .true.
      InvalidOutput( BAeroFy(  i) ) = .true.
      InvalidOutput( BAeroFz(  i) ) = .true.
      InvalidOutput( BAeroMx(  i) ) = .true.
      InvalidOutput( BAeroMy(  i) ) = .true.
      InvalidOutput( BAeroMz(  i) ) = .true.
      InvalidOutput( BNFbn(   :,i) ) = .true.
      InvalidOutput( BNFbt(   :,i) ) = .true.
      InvalidOutput( BNFbs(   :,i) ) = .true.
      InvalidOutput( BNMbn(   :,i) ) = .true.
      InvalidOutput( BNMbt(   :,i) ) = .true.
      InvalidOutput( BNMbs(   :,i) ) = .true.
               
   END DO
      
   DO I = p%NBlOuts+1,9  ! Invalid blade nodes
      
      InvalidOutput( BNVUndx(i,:) ) = .true.
      InvalidOutput( BNVUndy(i,:) ) = .true.
      InvalidOutput( BNVUndz(i,:) ) = .true.
      InvalidOutput( BNVDisx(i,:) ) = .true.
      InvalidOutput( BNVDisy(i,:) ) = .true.
      InvalidOutput( BNVDisz(i,:) ) = .true.
      InvalidOutput( BNSTVx( i,:) ) = .true.
      InvalidOutput( BNSTVy( i,:) ) = .true.
      InvalidOutput( BNSTVz( i,:) ) = .true.
      InvalidOutput( BNVRel( i,:) ) = .true.
      InvalidOutput( BNDynP( i,:) ) = .true.
      InvalidOutput( BNRe(   i,:) ) = .true.
      InvalidOutput( BNM(    i,:) ) = .true.
      InvalidOutput( BNVIndx(i,:) ) = .true.
      InvalidOutput( BNVIndy(i,:) ) = .true.
      InvalidOutput( BNAxInd(i,:) ) = .true.
      InvalidOutput( BNTnInd(i,:) ) = .true.
      InvalidOutput( BNAlpha(i,:) ) = .true.
      InvalidOutput( BNTheta(i,:) ) = .true.
      InvalidOutput( BNPhi(  i,:) ) = .true.
      InvalidOutput( BNCurve(i,:) ) = .true.
      InvalidOutput( BNCl(   i,:) ) = .true.
      InvalidOutput( BNCd(   i,:) ) = .true.
      InvalidOutput( BNCm(   i,:) ) = .true.
      InvalidOutput( BNCx(   i,:) ) = .true.
      InvalidOutput( BNCy(   i,:) ) = .true.
      InvalidOutput( BNCn(   i,:) ) = .true.
      InvalidOutput( BNCt(   i,:) ) = .true.
      InvalidOutput( BNFl(   i,:) ) = .true.
      InvalidOutput( BNFd(   i,:) ) = .true.
      InvalidOutput( BNMm(   i,:) ) = .true.
      InvalidOutput( BNFx(   i,:) ) = .true.
      InvalidOutput( BNFy(   i,:) ) = .true.
      InvalidOutput( BNFn(   i,:) ) = .true.
      InvalidOutput( BNFt(   i,:) ) = .true.
      InvalidOutput( BNClrnc(i,:) ) = .true.
      InvalidOutput( BNGam(  i,:) ) = .true.
      InvalidOutput( BNSgCav(i,:) ) = .true.
      InvalidOutput( BNSigCr(i,:) ) = .true.
      InvalidOutput( BNCpMin(i,:) ) = .true.
      InvalidOutput( BNFbn(  i,:) ) = .true.
      InvalidOutput( BNFbt(  i,:) ) = .true.
      InvalidOutput( BNFbs(  i,:) ) = .true.
      InvalidOutput( BNMbn(  i,:) ) = .true.
      InvalidOutput( BNMbt(  i,:) ) = .true.
      InvalidOutput( BNMbs(  i,:) ) = .true.
         
   END DO

!   ................. End of validity checking .................


   !-------------------------------------------------------------------------------------------------
   ! Allocate and set index, name, and units for the output channels
   ! If a selected output channel is not available in this module, set error flag.
   !-------------------------------------------------------------------------------------------------

   ALLOCATE ( p%OutParam(0:p%NumOuts) , STAT=ErrStat2 )
   IF ( ErrStat2 /= 0_IntKi )  THEN
      CALL SetErrStat( ErrID_Fatal,"Error allocating memory for the AeroDyn OutParam array.", ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

      ! Set index, name, and units for the time output channel:

   p%OutParam(0)%Indx  = Time
   p%OutParam(0)%Name  = "Time"    ! OutParam(0) is the time channel by default.
   p%OutParam(0)%Units = "(s)"
   p%OutParam(0)%SignM = 1


      ! Set index, name, and units for all of the output channels.
      ! If a selected output channel is not available by this module set ErrStat = ErrID_Warn.

   DO I = 1,p%NumOuts

      p%OutParam(I)%Name  = OutList(I)

      Indx = FindValidChannelIndx(OutList(I), ValidParamAry, p%OutParam(I)%SignM)

      IF ( Indx > 0 ) THEN ! we found the channel name
         IF ( InvalidOutput( ParamIndxAry(Indx) ) ) THEN  ! but, it isn't valid for these settings
            p%OutParam(I)%Indx  = 0                 ! pick any valid channel (I just picked "Time=0" here because it's universal)
            p%OutParam(I)%Units = "INVALID"
            p%OutParam(I)%SignM = 0
         ELSE
            p%OutParam(I)%Indx  = ParamIndxAry(Indx)
            p%OutParam(I)%Units = ParamUnitsAry(Indx) ! it's a valid output
         END IF
      ELSE ! this channel isn't valid
         p%OutParam(I)%Indx  = 0                    ! pick any valid channel (I just picked "Time=0" here because it's universal)
         p%OutParam(I)%Units = "INVALID"
         p%OutParam(I)%SignM = 0                    ! multiply all results by zero

         CALL SetErrStat(ErrID_Fatal, TRIM(p%OutParam(I)%Name)//" is not an available output channel.",ErrStat,ErrMsg,RoutineName)
      END IF

   END DO

   RETURN
END SUBROUTINE SetOutParam
!----------------------------------------------------------------------------------------------------------------------------------
!End of code generated by Matlab script
!**********************************************************************************************************************************


!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up the information needed for plotting VTK surfaces.
subroutine AD_SetVTKSurface(InitOutData_AD, u_AD, VTK_Surface, errStat, errMsg)
   type(AD_InitOutputType),       intent(inout) :: InitOutData_AD   !< The initialization output from AeroDyn
   type(AD_InputType),  target,   intent(in   ) :: u_AD
   type(AD_VTK_RotSurfaceType), allocatable, intent(inout) :: VTK_Surface(:) ! VTK_Surface for each rotor
   integer(IntKi),               intent(  out) :: errStat          !< Error status of the operation
   character(*),                 intent(  out) :: errMsg           !< Error message if errStat /= ErrID_None
   real(SiKi)              :: TwrDiam_top, TwrDiam_base, TwrRatio, TwrLength
   integer(IntKi)          :: topNode, baseNode, cylNode, tipNode, rootNode
   integer(IntKi)          :: k
   character(1024)         :: vtkroot
   integer(IntKi)          :: iWT
   integer(IntKi)          :: errStat2
   character(ErrMsgLen)    :: errMsg2
   type(MeshType), pointer :: Mesh
   integer(IntKi)          :: nBlades
   integer(IntKi)          :: nWT
   errStat = ErrID_None
   errMsg  = ""
   
   if (allocated(VTK_Surface)) then
      return ! The surfaces were already computed (for combined cases)
   endif

   nWT = size(u_AD%rotors)
   allocate(VTK_Surface(nWT),stat=errStat2); errMsg2='Allocating VTK_Surface'; if(Failed()) return

   ! --- Create surfaces for Nacelle, Base, Tower, Blades
   do iWT=1,nWT
      !.......................
      ! tapered tower
      !.......................
      Mesh=>u_AD%rotors(iWT)%TowerMotion
      if(associated(Mesh)) then
         if (Mesh%NNodes>0) then
            CALL AllocAry(VTK_Surface(iWT)%TowerRad, Mesh%NNodes,'VTK_Surface(iWT)%TowerRad',errStat2,errMsg2)
            topNode   = Mesh%NNodes - 1
            !baseNode  = Mesh%refNode
            baseNode  = 1 ! TODO TODO
            TwrLength = TwoNorm( Mesh%position(:,topNode) - Mesh%position(:,baseNode) ) ! this is the assumed length of the tower
            TwrRatio  = TwrLength / 87.6_SiKi  ! use ratio of the tower length to the length of the 5MW tower
            TwrDiam_top  = 3.87*TwrRatio
            TwrDiam_base = 6.0*TwrRatio
            TwrRatio = 0.5 * (TwrDiam_top - TwrDiam_base) / TwrLength
            do k=1,Mesh%NNodes
               TwrLength = TwoNorm( Mesh%position(:,k) - Mesh%position(:,baseNode) ) 
               VTK_Surface(iWT)%TowerRad(k) = 0.5*TwrDiam_Base + TwrRatio*TwrLength
            end do
         else
            !print*,'>>>> TOWER HAS NO NODES'
            ! TODO create a fake tower
            CALL AllocAry(VTK_Surface(iWT)%TowerRad, 0,'VTK_Surface(iWT)%TowerRad',errStat2,errMsg2)
         endif
      endif

      !.......................
      ! blade surfaces
      !.......................
      nBlades = size(u_AD%rotors(iWT)%BladeMotion)
      allocate(VTK_Surface(iWT)%BladeShape(nBlades),stat=errStat2); errMsg2='Allocating BladeShape'; if(Failed()) return
      if (allocated(InitOutData_AD%rotors(iWT)%BladeShape)) THEN
         do k=1, nBlades
            call move_alloc( InitOutData_AD%rotors(iWT)%BladeShape(k)%AirfoilCoords, VTK_Surface(iWT)%BladeShape(k)%AirfoilCoords )
         end do
      else
         call WrScr('Profile coordinates missing, using dummy coordinates for blade surface VTK outputs')
         rootNode = 1
         do K=1, nBlades
            tipNode  = u_AD%rotors(iWT)%BladeMotion(K)%NNodes
            cylNode  = min(3,u_AD%rotors(iWT)%BladeMotion(K)%Nnodes)
            call AD_SetVTKDefaultBladeParams(u_AD%rotors(iWT)%BladeMotion(K), VTK_Surface(iWT)%BladeShape(K), tipNode, rootNode, cylNode, errStat2, errMsg2, BlChord=InitOutData_AD%rotors(iWT)%BladeProps(k)%BlChord); if (Failed()) return
         end do                           
      endif
   enddo ! iWT, turbines

contains

   logical function Failed()
        call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'AD_SetVTKSurface') 
        Failed =  errStat >= AbortErrLev
   end function Failed

end subroutine AD_SetVTKSurface
!----------------------------------------------------------------------------------------------------------------------------------
!> Writes AeroDyn VTK surfaces (tower, hub, blades)
subroutine AD_WrVTK_Surfaces(u_AD, y_AD, RefPoint, VTK_Surface, VTK_count, OutFileRoot, tWidth, numSectors, HubRad)
   type(AD_InputType),       intent(in   ) :: u_AD
   type(AD_OutputType),      intent(in   ) :: y_AD
   type(AD_VTK_RotSurfaceType), intent(in   ) :: VTK_Surface(:) ! VTK_Surface for each rotor
   real(SiKi),               intent(in   ) :: RefPoint(3)
   real(SiKi),               intent(in   ) :: HubRad
   integer(IntKi)          , intent(in   ) :: VTK_count
   character(len=*),         intent(in   ) :: OutFileRoot
   integer,                  intent(in   ) :: tWidth
   integer,                  intent(in   ) :: numSectors
   logical, parameter       :: OutputFields = .FALSE.          ! due to confusion about what fields mean on a surface, we are going to just output the basic meshes if people ask for fields
   integer(IntKi)           :: k
   integer(IntKi)           :: errStat2
   character(ErrMsgLen)     :: errMSg2
   integer(IntKi)           :: iWT
   integer(IntKi)           :: nBlades
   integer(IntKi)           :: nWT
   character(10)            :: sWT

   nWT = size(u_AD%rotors)
   do iWT = 1, nWT
      if (nWT==1) then
         sWT = ''
      else
         sWT = '.T'//trim(num2lstr(iWT))
      endif

      ! Tower motions
      if (u_AD%rotors(iWT)%TowerMotion%nNodes>0) then
         call MeshWrVTK_Ln2Surface (RefPoint, u_AD%rotors(iWT)%TowerMotion, trim(OutFileRoot)//trim(sWT)//'.TowerSurface', &
                                    VTK_count, OutputFields, errStat2, errMsg2, tWidth, numSectors, VTK_Surface(iWT)%TowerRad )
      endif

      nBlades = size(u_AD%rotors(iWT)%BladeMotion)
    
      if (nBlades>0) then
         ! Hub
         call MeshWrVTK_PointSurface (RefPoint, u_AD%rotors(iWT)%HubMotion, trim(OutFileRoot)//trim(sWT)//'.HubSurface', &
                                      VTK_count, OutputFields, errStat2, errMsg2, tWidth , &
                                      NumSegments=numSectors, radius=HubRad)
      endif

      ! Blades
      do K=1,nBlades
         call MeshWrVTK_Ln2Surface (RefPoint, u_AD%rotors(iWT)%BladeMotion(K), trim(OutFileRoot)//trim(sWT)//'.Blade'//trim(num2lstr(k))//'Surface', &
                                    VTK_count, OutputFields, errStat2, errMsg2, tWidth , verts=VTK_Surface(iWT)%BladeShape(K)%AirfoilCoords &
                                    ,Sib=y_AD%rotors(iWT)%BladeLoad(k) )
      end do                  
   enddo

end subroutine AD_WrVTK_Surfaces
!----------------------------------------------------------------------------------------------------------------------------------
!> Writes AeroDyn VTK Lines and points (tower, hub, blades)
subroutine AD_WrVTK_LinesPoints(u_AD, y_AD, RefPoint, VTK_count, OutFileRoot, tWidth)
   type(AD_InputType),       intent(in   ) :: u_AD
   type(AD_OutputType),      intent(in   ) :: y_AD
   real(SiKi),               intent(in   ) :: RefPoint(3)
   integer(IntKi),           intent(in   ) :: VTK_count
   character(len=*),         intent(in   ) :: OutFileRoot
   integer,                  intent(in   ) :: tWidth
   logical, parameter       :: OutputFields = .TRUE.
   integer(IntKi)           :: k
   integer(IntKi)           :: errStat2
   character(ErrMsgLen)     :: errMSg2
   integer(IntKi)           :: iWT
   integer(IntKi)           :: nBlades
   integer(IntKi)           :: nWT
   character(10)            :: sWT

   nWT = size(u_AD%rotors)
   do iWT = 1, nWT
      if (nWT==1) then
         sWT = ''
      else
         sWT = '.T'//trim(num2lstr(iWT))
      endif

      ! Tower motions
      if (u_AD%rotors(iWT)%TowerMotion%nNodes>0) then
         call MeshWrVTK(RefPoint, u_AD%rotors(iWT)%TowerMotion, trim(OutFileRoot)//trim(sWT)//'.Tower', &
                        VTK_count, OutputFields, errStat2, errMsg2, tWidth )
      endif

      nBlades = size(u_AD%rotors(iWT)%BladeMotion)
    
      if (nBlades>0) then
         ! Hub
         call MeshWrVTK(RefPoint, u_AD%rotors(iWT)%HubMotion, trim(OutFileRoot)//trim(sWT)//'.Hub', &
                        VTK_count, OutputFields, errStat2, errMsg2, tWidth )
      endif

      ! Blades
      do K=1,nBlades
         call MeshWrVTK(RefPoint, u_AD%rotors(iWT)%BladeMotion(K), trim(OutFileRoot)//trim(sWT)//'.Blade'//trim(num2lstr(k)), &
                        VTK_count, OutputFields, errStat2, errMsg2, tWidth, Sib=y_AD%rotors(iWT)%BladeLoad(k) )
      end do                  
   enddo

end subroutine AD_WrVTK_LinesPoints
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine comes up with some default airfoils for blade surfaces for a given blade mesh, M.
SUBROUTINE AD_SetVTKDefaultBladeParams(M, BladeShape, tipNode, rootNode, cylNode, errStat, errMsg, BlChord)
   TYPE(MeshType),               INTENT(IN   ) :: M                !< The Mesh the defaults should be calculated for
   TYPE(AD_VTK_BLSurfaceType), INTENT(INOUT) :: BladeShape       !< BladeShape to set to default values
   INTEGER(IntKi),               INTENT(IN   ) :: rootNode         !< Index of root node (innermost node) for this mesh
   INTEGER(IntKi),               INTENT(IN   ) :: tipNode          !< Index of tip node (outermost node) for this mesh
   INTEGER(IntKi),               INTENT(IN   ) :: cylNode          !< Index of last node to have a cylinder shape
   INTEGER(IntKi),               INTENT(  OUT) :: errStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT) :: errMsg           !< Error message if errStat /= ErrID_None
   REAL(ReKi),  OPTIONAL,        INTENT(IN   ) :: BlChord(:)
   REAL(SiKi)                                  :: bladeLength, chord, pitchAxis
   REAL(SiKi)                                  :: bladeLengthFract, bladeLengthFract2, ratio, posLength ! temporary quantities               
   REAL(SiKi)                                  :: cylinderLength, x, y, angle               
   INTEGER(IntKi)                              :: i, j
   INTEGER(IntKi)                              :: errStat2
   CHARACTER(ErrMsgLen)                        :: errMsg2
   CHARACTER(*), PARAMETER                     :: RoutineName = 'SetVTKDefaultBladeParams'
   integer, parameter :: N = 66
   ! default airfoil shape coordinates; uses S809 values from http://wind.nrel.gov/airfoils/Shapes/S809_Shape.html:   
   real, parameter, dimension(N) :: xc=(/ 1.0,0.996203,0.98519,0.967844,0.945073,0.917488,0.885293,0.848455,0.80747,0.763042,0.715952,0.667064,0.617331,0.56783,0.519832,0.474243,0.428461,0.382612,0.33726,0.29297,0.250247,0.209576,0.171409,0.136174,0.104263,0.076035,0.051823,0.03191,0.01659,0.006026,0.000658,0.000204,0.0,0.000213,0.001045,0.001208,0.002398,0.009313,0.02323,0.04232,0.065877,0.093426,0.124111,0.157653,0.193738,0.231914,0.271438,0.311968,0.35337,0.395329,0.438273,0.48192,0.527928,0.576211,0.626092,0.676744,0.727211,0.776432,0.823285,0.86663,0.905365,0.938474,0.965086,0.984478,0.996141,1.0 /)
   real, parameter, dimension(N) :: yc=(/ 0.0,0.000487,0.002373,0.00596,0.011024,0.017033,0.023458,0.03028,0.037766,0.045974,0.054872,0.064353,0.074214,0.084095,0.093268,0.099392,0.10176,0.10184,0.10007,0.096703,0.091908,0.085851,0.078687,0.07058,0.061697,0.052224,0.042352,0.032299,0.02229,0.012615,0.003723,0.001942,-0.00002,-0.001794,-0.003477,-0.003724,-0.005266,-0.011499,-0.020399,-0.030269,-0.040821,-0.051923,-0.063082,-0.07373,-0.083567,-0.092442,-0.099905,-0.105281,-0.108181,-0.108011,-0.104552,-0.097347,-0.086571,-0.073979,-0.060644,-0.047441,-0.0351,-0.024204,-0.015163,-0.008204,-0.003363,-0.000487,0.000743,0.000775,0.00029,0.0 /)
   call AllocAry(BladeShape%AirfoilCoords, 2, N, M%NNodes, 'BladeShape%AirfoilCoords', errStat2, errMsg2)
      CALL SetErrStat(errStat2,errMsg2,errStat,errMsg,RoutineName)
      IF (errStat >= AbortErrLev) RETURN
   ! Chord length and pitch axis location are given by scaling law
   bladeLength       = TwoNorm( M%position(:,tipNode) - M%Position(:,rootNode) )
   cylinderLength    = TwoNorm( M%Position(:,cylNode) - M%Position(:,rootNode) )
   bladeLengthFract  = 0.22*bladeLength
   bladeLengthFract2 = bladeLength-bladeLengthFract != 0.78*bladeLength
   DO i=1,M%Nnodes
      posLength = TwoNorm( M%Position(:,i) - M%Position(:,rootNode) )
      IF (posLength .LE. bladeLengthFract) THEN
         ratio     = posLength/bladeLengthFract
         chord     =  (0.06 + 0.02*ratio)*bladeLength
         pitchAxis =   0.25 + 0.125*ratio
      ELSE
         chord     = (0.08 - 0.06*(posLength-bladeLengthFract)/bladeLengthFract2)*bladeLength
         pitchAxis = 0.375
      END IF
      if(present(BlChord)) then
         chord = BlChord(i)
      endif
      IF (posLength .LE. cylinderLength) THEN 
         ! create a cylinder for this node
         chord = chord/2.0_SiKi
         DO j=1,N
            ! normalized x,y coordinates for airfoil
            x = yc(j)
            y = xc(j) - 0.5
            angle = ATAN2( y, x)
               ! x,y coordinates for cylinder
            BladeShape%AirfoilCoords(1,j,i) = chord*COS(angle) ! x (note that "chord" is really representing chord/2 here)
            BladeShape%AirfoilCoords(2,j,i) = chord*SIN(angle) ! y (note that "chord" is really representing chord/2 here)
         END DO                                                     
      ELSE
         ! create an airfoil for this node
         DO j=1,N                  
            ! normalized x,y coordinates for airfoil, assuming an upwind turbine
            x = yc(j)
            y = xc(j) - pitchAxis
               ! x,y coordinates for airfoil
            BladeShape%AirfoilCoords(1,j,i) =  chord*x
            BladeShape%AirfoilCoords(2,j,i) =  chord*y                        
         END DO
      END IF
   END DO ! nodes on mesh
         
END SUBROUTINE AD_SetVTKDefaultBladeParams


!----------------------------------------------------------------------------------------------------------------------------------
!> Set the cant angle from the mid-chord reference line
!! SETCANTANGLE() will update the BlCrvAng based upon the projection of
!! the mid-chord reference line onto the X-Z plane.
!!
!! NOTE: this assumes that the cant and toe angles are zero and only the
!! local twist and chord length contirbute to the mid-chord location; in
!! the future an iterative approach could be taken to include the cant
!! angle influence on the mid-code node locations
subroutine setCantAngle( BladeKInputFileData )
! Note: we have already checked that BladeKInputFileData%NumBlNds > 2 and that zMidChord is in increasing order (not constant)
   TYPE(AD_BladePropsType),   INTENT(INOUT)  :: BladeKInputFileData

   REAL(ReKi)                                :: dx(          BladeKInputFileData%NumBlNds)
   REAL(ReKi)                                :: xMidChord(   BladeKInputFileData%NumBlNds)
   REAL(ReKi)                                :: zMidChord(   BladeKInputFileData%NumBlNds)
   REAL(ReKi)                                :: prebendSlope(BladeKInputFileData%NumBlNds)
   INTEGER                                   :: NumBlNds
   INTEGER                                   :: ii          ! loop counter
   
   NumBlNds = BladeKInputFileData%NumBlNds
   
   ! Compute mid-chord location in global system
   dx = 0.25_ReKi * BladeKInputFileData%BlChord * sin( BladeKInputFileData%BlTwist )  ! note element-by-element multiplication (not matrix multiply here); twist is in radians
            
   xMidChord = BladeKInputFileData%BlCrvAC + dx
   zMidChord = BladeKInputFileData%BlSpn

   
   ! Compute prebend slope relative to z-span

   ! Root node
   prebendSlope(1) = dfdxOfLagrangeInterpolant( zMidChord(1),               &
                                                zMidChord(1), xMidChord(1), &
                                                zMidChord(2), xMidChord(2), &
                                                zMidChord(3), xMidChord(3) )

   ! Internal nodes
   do ii = 2, NumBlNds - 1
         prebendSlope(ii) = dfdxOfLagrangeInterpolant( zMidChord(ii),                     &
                                                       zMidChord(ii-1), xMidChord(ii-1),  &
                                                       zMidChord(ii),   xMidChord(ii),    &
                                                       zMidChord(ii+1), xMidChord(ii+1) )
   end do

   ! Tip node
   prebendSlope(NumBlNds) = dfdxOfLagrangeInterpolant( zMidChord(NumBlNds),                          &
                                                       zMidChord(NumBlNds-2), xMidChord(NumBlNds-2), &
                                                       zMidChord(NumBlNds-1), xMidChord(NumBlNds-1), &
                                                       zMidChord(NumBlNds),   xMidChord(NumBlNds) )

   ! Convert slope to cant angle
   BladeKInputFileData%BlCrvAng = atan2( prebendSlope, 1.0_ReKi ) ! return value in radians
        
end subroutine setCantAngle
!----------------------------------------------------------------------------------------------------------------------------------
!> df/dx approximation from Lagrange interpolating polynomial
!! See Eqn 5 at https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html
REAL(ReKi) function dfdxOfLagrangeInterpolant(x, x1, f1, x2, f2, x3, f3) RESULT(dfdx)

   REAL(ReKi) , INTENT(IN)       :: x
   REAL(ReKi) , INTENT(IN)       :: x1
   REAL(ReKi) , INTENT(IN)       :: x2
   REAL(ReKi) , INTENT(IN)       :: x3
   REAL(ReKi) , INTENT(IN)       :: f1
   REAL(ReKi) , INTENT(IN)       :: f2
   REAL(ReKi) , INTENT(IN)       :: f3

   dfdx = f1 * (2*x-x2-x3)/((x1-x2)*(x1-x3)) &
        + f2 * (2*x-x1-x3)/((x2-x1)*(x2-x3)) &
        + f3 * (2*x-x1-x2)/((x3-x1)*(x3-x2));

end function dfdxOfLagrangeInterpolant
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE AeroDyn_IO
