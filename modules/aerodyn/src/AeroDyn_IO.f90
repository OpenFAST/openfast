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
SUBROUTINE Calc_WriteOutput( p, p_AD, u, x, m, m_AD, y, OtherState, xd, indx, iRot, ErrStat, ErrMsg )
   
   TYPE(RotParameterType),       INTENT(IN   )  :: p                                 ! The rotor parameters
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p_AD                              ! The module parameters
   TYPE(RotInputType),           INTENT(IN   )  :: u                                 ! inputs
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
   
   tmpHubFB  = 0.0_ReKi
   tmpHubMB  = 0.0_ReKi

   ! Compute max radius and rotor speed
   if (p_AD%WakeMod /= WakeMod_FVW) then
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

   
   call Calc_WriteOutput_AD() ! need to call this before calling the BEMT vs FVW versions of outputs so that the integrated output quantities are known
   
   if (p_AD%WakeMod /= WakeMod_FVW) then
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
      
         tmp = matmul( u%TowerMotion%Orientation(:,:,j) , u%InflowOnTower(:,j) )
         m%AllOuts( TwNVUnd(:,beta) ) = tmp
      
         tmp = matmul( u%TowerMotion%Orientation(:,:,j) , u%TowerMotion%TranslationVel(:,j) )
         m%AllOuts( TwNSTV(:,beta) ) = tmp
      
         m%AllOuts( TwNVrel(beta) ) = m%W_Twr(j)                           ! relative velocity   
         m%AllOuts( TwNDynP(beta) ) = 0.5 * p%AirDens * m%W_Twr(j)**2      ! dynamic pressure
         m%AllOuts( TwNRe(  beta) ) = p%TwrDiam(j) * m%W_Twr(j) / p%KinVisc / 1.0E6 ! reynolds number (in millions)
         m%AllOuts( TwNM(   beta) ) = m%W_Twr(j) / p%SpdSound               ! Mach number
         m%AllOuts( TwNFdx( beta) ) = m%X_Twr(j)         
         m%AllOuts( TwNFdy( beta) ) = m%Y_Twr(j)         
      
         if ( p%Buoyancy ) then
            tmp = matmul( u%TowerMotion%Orientation(:,:,j) , m%TwrBuoyLoad%Force(:,j) )
            m%AllOuts( TwNFbx(beta) ) = tmp(1)
            m%AllOuts( TwNFby(beta) ) = tmp(2)
            m%AllOuts( TwNFbz(beta) ) = tmp(3)
   
            tmp = matmul( u%TowerMotion%Orientation(:,:,j) , m%TwrBuoyLoad%Moment(:,j) )
            m%AllOuts( TwNMbx(beta) ) = tmp(1)
            m%AllOuts( TwNMby(beta) ) = tmp(2)
            m%AllOuts( TwNMbz(beta) ) = tmp(3)
         end if
   
      end do ! out nodes
   
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
      end if
   
         ! nacelle outputs
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

         ! blade outputs
      do k=1,min(p%numBlades,AD_MaxBl_Out)    ! limit this
         do beta=1,p%NBlOuts
            j=p%BlOutNd(beta)

            tmp = matmul( m%orientationAnnulus(:,:,j,k), u%InflowOnBlade(:,j,k) )
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
         
            m%AllOuts( BNCurve(beta,k) ) = m%Curve(j,k)*R2D
                  
            m%AllOuts( BNSigCr(   beta,k) ) = m%SigmaCavitCrit(j,k)
            m%AllOuts( BNSgCav(   beta,k) ) = m%SigmaCavit(j,k)

            if ( p%Buoyancy ) then
               tmp = matmul( u%BladeMotion(k)%Orientation(:,:,j), m%BladeBuoyLoad(k)%Force(:,j) )
               m%AllOuts( BNFbn(beta,k) ) = tmp(1)
               m%AllOuts( BNFbt(beta,k) ) = tmp(2)
               m%AllOuts( BNFbs(beta,k) ) = tmp(3)

               tmp = matmul( u%BladeMotion(k)%Orientation(:,:,j), m%BladeBuoyLoad(k)%Moment(:,j) )
               m%AllOuts( BNMbn(beta,k) ) = tmp(1)
               m%AllOuts( BNMbt(beta,k) ) = tmp(2)
               m%AllOuts( BNMbs(beta,k) ) = tmp(3)
            end if

         end do ! nodes
      end do ! blades
   


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
         
         if (k<=size(BAeroFxg)) then
            ! Power contribution of blade wrt hub
            tmp = matmul( u%HubMotion%Orientation(:,:,1), m%HubLoad%moment(:,1) )
            m%AllOuts( BAeroPwr(k) ) = omega * tmp(1)
            
            ! In global, wrt hub! 
            m%AllOuts( BAeroFxg(k) ) = m%HubLoad%force(1,1)
            m%AllOuts( BAeroFyg(k) ) = m%HubLoad%force(2,1)
            m%AllOuts( BAeroFzg(k) ) = m%HubLoad%force(3,1)
            m%AllOuts( BAeroMxg(k) ) = m%HubLoad%moment(1,1)
            m%AllOuts( BAeroMyg(k) ) = m%HubLoad%moment(2,1)
            m%AllOuts( BAeroMzg(k) ) = m%HubLoad%moment(3,1)
         end if
      end do

        ! In global
      m%AllOuts( RtAeroFxg ) = force(1)
      m%AllOuts( RtAeroFyg ) = force(2)
      m%AllOuts( RtAeroFzg ) = force(3)
      m%AllOuts( RtAeroMxg ) = moment(1)
      m%AllOuts( RtAeroMyg ) = moment(2)
      m%AllOuts( RtAeroMzg ) = moment(3)
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
      if ( EqualRealNos( real(m%V_dot_x,SiKi), 0.0_SiKi ) ) then
         m%AllOuts( RtTSR )    = 0.0_ReKi
         m%AllOuts( RtAeroCp ) = 0.0_ReKi
         m%AllOuts( RtAeroCq ) = 0.0_ReKi
         m%AllOuts( RtAeroCt ) = 0.0_ReKi
      else
         m%AllOuts( RtTSR )    = omega * rmax / m%V_dot_x

         denom = 0.5*p%AirDens*m%AllOuts( RtArea )*m%V_dot_x**2
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
      REAL(ReKi)                                   :: denom !, rmax
      REAL(ReKi)                                   :: ct, st ! cosine, sine of theta
      REAL(ReKi)                                   :: cp, sp ! cosine, sine of phi
 



   
         ! blade outputs
      do k=1,min(p%numBlades,AD_MaxBl_Out)    ! limit this
         m%AllOuts( BAzimuth(k) ) = MODULO( m%BEMT_u(indx)%psi(k)*R2D, 360.0_ReKi )
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

   !        m%AllOuts( BNCurve(beta,k) ) = m%Curve(j,k)*R2D

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
   
      m%AllOuts( RtSkew   ) = m%BEMT_u(indx)%chi0*R2D
!     m%AllOuts( RtTSR    ) = m%BEMT_u(indx)%TSR
      m%AllOuts( DBEMTau1 ) = OtherState%BEMT%DBEMT%tau1

      
   end subroutine Calc_WriteOutput_BEMT

   !..........................................................................................
   !> Similar to Calc_WriteOutput_BEMT. TODO Merge me
   !! NOTE: relies on the prior calculation of m%V_dot_x, and m%V_diskAvg (done in DiskAvgValues)
   !!                                          m%DisturbedInflow (done in SetInputs)
   !!       Make sure these are set!
   subroutine Calc_WriteOutput_FVW
      integer    :: iW

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
!             m%AllOuts( BNCurve(beta,k) ) = m%Curve(j,k)*R2D ! TODO

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

!      m%AllOuts( DBEMTau1 ) = 0.0_ReKi ! not valid with FVW

   end subroutine Calc_WriteOutput_FVW

END SUBROUTINE Calc_WriteOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadInputFiles( InputFileName, InputFileData, Default_DT, OutFileRoot, NumBlades, AeroProjMod, UnEcho, ErrStat, ErrMsg )
! This subroutine reads the input file and stores all the data in the AD_InputFile structure.
! It does not perform data validation.
!..................................................................................................................................

      ! Passed variables
   REAL(DbKi),              INTENT(IN)    :: Default_DT      ! The default DT (from glue code)

   CHARACTER(*),            INTENT(IN)    :: InputFileName   ! Name of the input file
   CHARACTER(*),            INTENT(IN)    :: OutFileRoot     ! The rootname of all the output files written by this routine.

   TYPE(AD_InputFile),      INTENT(INOUT) :: InputFileData   ! Data stored in the module's input file
   INTEGER(IntKi),          INTENT(INOUT) :: UnEcho          ! Unit number for the echo file

   INTEGER(IntKi),          INTENT(IN)    :: NumBlades(:)    ! Number of blades per rotor 
   INTEGER(IntKi),          INTENT(IN)    :: AeroProjMod(:)  ! AeroProjMod per rotor
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
         CALL ReadBladeInputs ( InputFileData%ADBlFile(iBld), InputFileData%rotors(iR)%BladeProps(I), AeroProjMod(iR), UnEcho, ErrStat2, ErrMsg2 )
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

   character(*), parameter                         :: RoutineName = 'ParsePrimaryFileInfo'

   ! Initialization
   ErrStat  =  ErrId_None
   ErrMsg   =  ""
   UnEc   = -1     ! Echo file unit.  >0 when used


   CALL AllocAry( InputFileData%OutList, MaxOutPts, "Outlist", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Allocate array for holding the list of node outputs
   CALL AllocAry( InputFileData%BldNd_OutList, BldNd_MaxOutPts, "BldNd_Outlist", ErrStat2, ErrMsg2 )
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
      ! WakeMod - Type of wake/induction model (switch) {0=none, 1=BEMT, 2=DBEMT, 3=OLAF}  [WakeMod cannot be 2 or 3 when linearizing]
   call ParseVar( FileInfo_In, CurLine, "WakeMod", InputFileData%WakeMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! AFAeroMod - Type of blade airfoil aerodynamics model (switch) {1=steady model, 2=Beddoes-Leishman unsteady model} [AFAeroMod must be 1 when linearizing]
   call ParseVar( FileInfo_In, CurLine, "AFAeroMod", InputFileData%AFAeroMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! TwrPotent - Type of tower influence on wind based on potential flow around the tower (switch) {0=none, 1=baseline potential flow, 2=potential flow with Bak correction}
   call ParseVar( FileInfo_In, CurLine, "TwrPotent", InputFileData%TwrPotent, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! TwrShadow - Type of tower influence on wind based on downstream tower shadow {0=none, 1=Powles model, 2=Eames model}
   call ParseVar( FileInfo_In, CurLine, "TwrShadow", InputFileData%TwrShadow, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! TwrAero - Calculate tower aerodynamic loads? (flag)
   call ParseVar( FileInfo_In, CurLine, "TwrAero", InputFileData%TwrAero, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! FrozenWake - Assume frozen wake during linearization? (flag) [used only when WakeMod=1 and when linearizing]
   call ParseVar( FileInfo_In, CurLine, "FrozenWake", InputFileData%FrozenWake, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! CavitCheck - Perform cavitation check? (flag) [AFAeroMod must be 1 when CavitCheck=true]
   call ParseVar( FileInfo_In, CurLine, "CavitCheck", InputFileData%CavitCheck, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! Buoyancy - Include buoyancy effects? (flag)
   call ParseVar( FileInfo_In, CurLine, "Buoyancy", InputFileData%Buoyancy, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! CompAA - Flag to compute AeroAcoustics calculation [only used when WakeMod=1 or 2]
   call ParseVar( FileInfo_In, CurLine, "CompAA", InputFileData%CompAA, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! AA_InputFile - Aeroacoustics input file
   call ParseVar( FileInfo_In, CurLine, "AA_InputFile", InputFileData%AA_InputFile, ErrStat2, ErrMsg2, UnEc )
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
      ! Patm - Atmospheric pressure {or default} (Pa) [used only when CavitCheck=True]
   call ParseVarWDefault( FileInfo_In, CurLine, "Patm", InputFileData%Patm, InitInp%defPatm, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! Pvap - Vapour pressure of fluid {or default} (Pa) [used only when CavitCheck=True]
   call ParseVarWDefault( FileInfo_In, CurLine, "Pvap", InputFileData%Pvap, InitInp%defPvap, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

   !======  Blade-Element/Momentum Theory Options  ====================================================== [unused when WakeMod=0 or 3]
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! SkewMod - Type of skewed-wake correction model (switch) {1=uncoupled, 2=Pitt/Peters, 3=coupled} [unused when WakeMod=0 or 3]
   call ParseVar( FileInfo_In, CurLine, "SkewMod", InputFileData%SkewMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! SkewModFactor - Constant used in Pitt/Peters skewed wake model {or "default" is 15/32*pi} (-) [used only when SkewMod=2; unused when WakeMod=0 or 3]
   call ParseVarWDefault( FileInfo_In, CurLine, "SkewModFactor", InputFileData%SkewModFactor, (15.0_ReKi * pi / 32.0_ReKi), ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
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

   !======  Dynamic Blade-Element/Momentum Theory Options  ============================================== [used only when WakeMod=2]
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
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
   call ParseVar( FileInfo_In, CurLine, "OLAFInputFileName", InputFileData%FVWFileName, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      IF ( PathIsRelative( InputFileData%FVWFileName ) ) InputFileData%FVWFileName = TRIM(PriPath)//TRIM(InputFileData%FVWFileName)

   !======  Beddoes-Leishman Unsteady Airfoil Aerodynamics Options  ===================================== [used only when AFAeroMod=2]
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! UAMod - Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez's variant (changes in Cn,Cc,Cm), 3=Minnema/Pierce variant (changes in Cc and Cm)} [used only when AFAeroMod=2]
   call ParseVar( FileInfo_In, CurLine, "UAMod", InputFileData%UAMod, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
      ! FLookup - Flag to indicate whether a lookup for f' will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files (flag) [used only when AFAeroMod=2]
   call ParseVar( FileInfo_In, CurLine, "FLookup", InputFileData%FLookup, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
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
   end do

   !======  Tail fin aerodynamics ========================================================================
   do iR = 1,size(NumBlades) ! Loop on rotors
      if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
      CurLine = CurLine + 1
      ! NOTE: being nice with legacy input file. Uncomment in next release
      call ParseVar(FileInfo_In, CurLine, "TFinAero", InputFileData%rotors(iR)%TFinAero, ErrStat2, ErrMsg2, UnEc); 
      if (ErrStat2==ErrID_None) then
         call ParseVar(FileInfo_In, CurLine, "TFinFile", InputFileData%rotors(iR)%TFinFile, ErrStat2, ErrMsg2, UnEc); if (Failed()) return
         InputFileData%rotors(iR)%TFinFile = trim(PriPath) // trim(InputFileData%rotors(iR)%TFinFile)
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
      ! BldNd_BlOutNd - Future feature will allow selecting a portion of the nodes to output.  Not implemented yet. (-)
      ! TODO: Parse this string into an array of nodes to output at (one idea is to set an array of boolean to T/F for which nodes to output).  At present, we ignore it entirely.
   call ParseVar( FileInfo_In, CurLine, "BldNd_BlOutNd", InputFileData%BldNd_BlOutNd_Str, ErrStat2, ErrMsg2, UnEc )
      if (FailedNodal()) return
      ! OutList - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   call ReadOutputListFromFileInfo( FileInfo_In, CurLine, InputFileData%BldNd_OutList, InputFileData%BldNd_NumOuts, ErrStat2, ErrMsg2, UnEc )
         if (FailedNodal()) return;

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
      FailedNodal = ErrStat2 >= AbortErrLev
      if ( FailedNodal ) then
         InputFileData%BldNd_BladesOut = 0
         InputFileData%BldNd_NumOuts = 0
         call wrscr( trim(ErrMsg_NoAllBldNdOuts) )
      endif
   end function FailedNodal
   subroutine LegacyWarning(Message)
      character(len=*), intent(in) :: Message
      call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      call WrScr('Warning: the AeroDyn input file is not at the latest format!' )
      call WrScr('         Visit: https://openfast.readthedocs.io/en/dev/source/user/api_change.html')
      call WrScr('> Issue: '//trim(Message))
      call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
   end subroutine LegacyWarning
   !-------------------------------------------------------------------------------------------------
END SUBROUTINE ParsePrimaryFileInfo
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadBladeInputs ( ADBlFile, BladeKInputFileData, AeroProjMod, UnEc, ErrStat, ErrMsg )
! This routine reads a blade input file.
!..................................................................................................................................


      ! Passed variables:

   TYPE(AD_BladePropsType),  INTENT(INOUT)  :: BladeKInputFileData                 ! Data for Blade K stored in the module's input file
   CHARACTER(*),             INTENT(IN)     :: ADBlFile                            ! Name of the blade input file data
   INTEGER(IntKi),           INTENT(IN)     :: AeroProjMod                         ! AeroProjMod
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                ! I/O unit for echo file. If present and > 0, write to UnEc

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message


      ! Local variables:

   INTEGER(IntKi)               :: I                                               ! A generic DO index.
   INTEGER( IntKi )             :: UnIn                                            ! Unit number for reading file
   INTEGER(IntKi)               :: ErrStat2 , IOS                                  ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                         ! Temporary Err msg
   CHARACTER(*), PARAMETER      :: RoutineName = 'ReadBladeInputs'
   CHARACTER(len=1024)          :: Line
   CHARACTER(len=50)            :: HeaderCols(10)                                  ! Header columns in file
   LOGICAL                      :: hasBuoyancy                                     ! Does file contain Buoyancy columns


   ErrStat = ErrID_None
   ErrMsg  = ""
   UnIn = -1
   
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


      ! Open the input file for blade K.

   CALL OpenFInpFile ( UnIn, ADBlFile, ErrStat2, ErrMsg2 )
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
   ! Check if 10 columns are present
   READ (Line,*, IOSTAT=ErrStat2) ( HeaderCols(I), I=1,10 )
   hasBuoyancy = .true.
   IF ( ErrStat2 < 0 )  THEN ! end of line reached
      hasBuoyancy = .false.
      !call WrScr('Blade input file is missing buoyancy columns.')
   ELSE IF ( ErrStat2 > 0 )  THEN
      CALL SetErrStat(ErrID_Fatal, 'Unexpected error while trying to infer column headers in blade file.', ErrStat, ErrMsg, RoutineName)
   endif


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

   IF (.not. hasBuoyancy) THEN
      BladeKInputFileData%BlCb    = 0.0_ReKi
      BladeKInputFileData%BlCenBn = 0.0_ReKi
      BladeKInputFileData%BlCenBt = 0.0_ReKi
   ENDIF
      
      ! Return on error if we didn't allocate space for the next inputs
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
            
   DO I=1,BladeKInputFileData%NumBlNds
      IF (hasBuoyancy) THEN
         READ( UnIn, *, IOStat=IOS ) BladeKInputFileData%BlSpn(I), BladeKInputFileData%BlCrvAC(I), BladeKInputFileData%BlSwpAC(I), &
                                     BladeKInputFileData%BlCrvAng(I), BladeKInputFileData%BlTwist(I), BladeKInputFileData%BlChord(I), &
                                     BladeKInputFileData%BlAFID(I), BladeKInputFileData%BlCb(I), BladeKInputFileData%BlCenBn(I), BladeKInputFileData%BlCenBt(I) 
      ELSE
         READ( UnIn, *, IOStat=IOS ) BladeKInputFileData%BlSpn(I), BladeKInputFileData%BlCrvAC(I), BladeKInputFileData%BlSwpAC(I), &
                                     BladeKInputFileData%BlCrvAng(I), BladeKInputFileData%BlTwist(I), BladeKInputFileData%BlChord(I), &
                                     BladeKInputFileData%BlAFID(I)
      ENDIF
         CALL CheckIOS( IOS, ADBlFile, 'Blade properties row '//TRIM(Num2LStr(I)), NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               ! Return on error if we couldn't read this line
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL Cleanup()
               RETURN
            END IF
         
         IF (UnEc > 0) THEN
            WRITE( UnEc, "(6(F9.4,1x),I9,4(F9.4,1x))", IOStat=IOS) BladeKInputFileData%BlSpn(I), BladeKInputFileData%BlCrvAC(I), BladeKInputFileData%BlSwpAC(I), &
                                  BladeKInputFileData%BlCrvAng(I), BladeKInputFileData%BlTwist(I), BladeKInputFileData%BlChord(I), &
                                  BladeKInputFileData%BlAFID(I), BladeKInputFileData%BlCb(I), BladeKInputFileData%BlCenBn(I), BladeKInputFileData%BlCenBt(I)
         END IF         
   END DO


   if (all(BladeKInputFileData%BlCrvAC.eq.0.0_ReKi)) then
        BladeKInputFileData%BlCrvAng = 0.0_ReKi
   else
      if (AeroProjMod==APM_BEM_NoSweepPitchTwist .or. AeroProjMod==APM_LiftingLine) then
         !call WrScr('>>> ReadBladeInputs: Not computing cant angle (BlCrvAng), AeroProjMod='//trim(num2lstr(AeroProjMod)))
      else if (AeroProjMod==APM_BEM_Polar) then
         call WrScr('>>> ReadBladeInputs: Computing cant angle (BlCrvAng), AeroProjMod='//trim(num2lstr(AeroProjMod)))
         call calcCantAngle(BladeKInputFileData%BlCrvAC,BladeKInputFileData%BlSpn,3, size(BladeKInputFileData%BlSpn),BladeKInputFileData%BlCrvAng) 
      else
         call SetErrStat(ErrID_Fatal, 'Unsupported AeroProjMod='//trim(num2lstr(AeroProjMod)), ErrStat, ErrMsg, RoutineName)
         call Cleanup()
         return
      endif
   endif
   BladeKInputFileData%BlCrvAng = BladeKInputFileData%BlCrvAng*D2R
   BladeKInputFileData%BlTwist  = BladeKInputFileData%BlTwist*D2R
                  
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

END SUBROUTINE ReadBladeInputs      
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
   call ParseVar(FileInfo_In, iLine, 'TFinChord' , TFData%TFinChord     , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinArea'  , TFData%TFinArea      , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_In, iLine, 'TFinRefP_n', TFData%TFinRefP_n, 3 , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_In, iLine, 'TFinAngles', TFData%TFinAngles, 3 , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinIndMod', TFData%TFinIndMod    , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   !====== Polar-based model ================================ [used only when TFinMod=1]
   call ParseCom(FileInfo_in, iLine, DummyLine                          , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_In, iLine, 'TFinAFID'  , TFData%TFinAFID      , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   !====== Unsteady slender body model ===================== [used only when TFinMod=2]
   call ParseCom(FileInfo_in, iLine, DummyLine                          , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   ! TODO

   ! --- Triggers
   TFData%TFinAngles = TFData%TFinAngles*D2R ! deg2rad

   ! --- Validation on the fly
   !if (all((/TFinAero_none,TFinAero_polar, TFinAero_USB/) /= TFData%TFinMod)) then
   if (all((/TFinAero_none,TFinAero_polar/) /= TFData%TFinMod)) then
      call Fatal('TFinMod needs to be 0, or 1')
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
SUBROUTINE AD_PrintSum( InputFileData, p, p_AD, u, y, ErrStat, ErrMsg )
! This routine generates the summary file, which contains a summary of input file options.
   use YAML, only: yaml_write_var
      ! passed variables
   TYPE(AD_InputFile),        INTENT(IN)  :: InputFileData                        ! Input-file data
   TYPE(RotParameterType),    INTENT(IN)  :: p                                    ! Parameters
   TYPE(AD_ParameterType),    INTENT(IN)  :: p_AD                                 ! Parameters
   TYPE(AD_InputType),        INTENT(IN)  :: u                                    ! inputs 
   TYPE(AD_OutputType),       INTENT(IN)  :: y                                    ! outputs
   INTEGER(IntKi),            INTENT(OUT) :: ErrStat
   CHARACTER(*),              INTENT(OUT) :: ErrMsg


      ! Local variables.

   INTEGER(IntKi)               :: I                                               ! Index for the nodes.
   INTEGER(IntKi)               :: UnSu                                            ! I/O unit number for the summary output file

   CHARACTER(*), PARAMETER      :: FmtDat    = '(A,T41,1(:,F13.3))'                ! Format for outputting mass and modal data.
   CHARACTER(*), PARAMETER      :: FmtDatT   = '(A,T35,1(:,F13.8))'                ! Format for outputting time steps.

   CHARACTER(30)                :: OutPFmt                                         ! Format to print list of selected output channels to summary file
   CHARACTER(30)                :: OutPFmtS                                        ! Format to print list of selected output channels to summary file
   CHARACTER(100)               :: Msg                                             ! temporary string for writing appropriate text to summary file

   CHARACTER(ChanLen),PARAMETER :: TitleStr(2) = (/ 'Parameter', 'Units    ' /)
   CHARACTER(ChanLen),PARAMETER :: TitleStrLines(2) = (/ '---------------', '---------------' /)

   ! Open the summary file and give it a heading.
      
   CALL GetNewUnit( UnSu, ErrStat, ErrMsg )
   CALL OpenFOutFile ( UnSu, TRIM( p%RootName )//'.sum', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Heading:
   WRITE (UnSu,'(/,A)')  'This summary information was generated by '//TRIM( GetNVD(AD_Ver) )// &
                         ' on '//CurDate()//' at '//CurTime()//'.'

   WRITE (UnSu,'(/,A)') '======  General Options  ============================================================================'
   ! WakeMod
   select case (p_AD%WakeMod)
      case (WakeMod_BEMT)
         Msg = 'Blade-Element/Momentum Theory'
      case (WakeMod_DBEMT)
         Msg = 'Dynamic Blade-Element/Momentum Theory'
      case (WakeMod_FVW)
         Msg = 'Free Vortex Wake Theory'
      case (WakeMod_None)
         Msg = 'steady'
      case default      
         Msg = 'unknown'      
   end select   
   WRITE (UnSu,Ec_IntFrmt) p_AD%WakeMod, 'WakeMod', 'Type of wake/induction model: '//TRIM(Msg)

   
   ! AFAeroMod
   select case (InputFileData%AFAeroMod)
      case (AFAeroMod_BL_unsteady)
         Msg = 'Beddoes-Leishman unsteady model'
      case (AFAeroMod_steady)
         Msg = 'steady'
      case default      
         Msg = 'unknown'      
   end select   
   WRITE (UnSu,Ec_IntFrmt) InputFileData%AFAeroMod, 'AFAeroMod', 'Type of blade airfoil aerodynamics model: '//TRIM(Msg)
   
   
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
   if (p%TwrAero) then
      Msg = 'Yes'
   else
      Msg = 'No'
   end if   
   WRITE (UnSu,Ec_LgFrmt) p%TwrAero, 'TwrAero', 'Calculate tower aerodynamic loads? '//TRIM(Msg)

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


   if (p_AD%WakeMod/=WakeMod_none) then
      WRITE (UnSu,'(A)') '======  Blade-Element/Momentum Theory Options  ======================================================'
      
      ! SkewMod 
      select case (InputFileData%SkewMod)
         case (SkewMod_Orthogonal)
            Msg = 'orthogonal'
         case (SkewMod_Uncoupled)
            Msg = 'uncoupled'
         case (SkewMod_PittPeters)
            Msg = 'Pitt/Peters' 
         case default      
            Msg = 'unknown'      
      end select
      WRITE (UnSu,Ec_IntFrmt) InputFileData%SkewMod, 'SkewMod', 'Type of skewed-wake correction model: '//TRIM(Msg)
      
      
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
      
      
      if (p_AD%WakeMod == WakeMod_DBEMT) then
         select case (InputFileData%DBEMT_Mod)
            case (DBEMT_tauConst)
               Msg = 'constant tau1'
            case (DBEMT_tauVaries)
               Msg = 'time-dependent tau1'
            case (DBEMT_cont_tauConst)
               Msg = 'continuous formulation with constant tau1'
            case default
               Msg = 'unknown'
         end select   
         
         WRITE (UnSu,Ec_IntFrmt) InputFileData%DBEMT_Mod, 'DBEMT_Mod', 'Type of dynamic BEMT (DBEMT) model: '//TRIM(Msg)
         
         if (InputFileData%DBEMT_Mod==DBEMT_tauConst) &
         WRITE (UnSu,Ec_ReFrmt) InputFileData%tau1_const, 'tau1_const', 'Time constant for DBEMT (s)'
         
      end if      
      
   end if
   
   if (InputFileData%AFAeroMod==AFAeroMod_BL_unsteady) then
      WRITE (UnSu,'(A)') '======  Beddoes-Leishman Unsteady Airfoil Aerodynamics Options  ====================================='
      
      ! UAMod
      select case (InputFileData%UAMod)
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
      WRITE (UnSu,Ec_IntFrmt) InputFileData%UAMod, 'UAMod', 'Unsteady Aero Model: '//TRIM(Msg)
   
   
      ! FLookup
      if (InputFileData%FLookup) then
         Msg = 'Yes'
      else
         Msg = 'No, use best-fit exponential equations instead'
      end if   
      WRITE (UnSu,Ec_LgFrmt) InputFileData%FLookup, 'FLookup', "Use a lookup for f'? "//TRIM(Msg)      
      
   end if
   
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
      
   if (p_AD%WakeMod /= WakeMod_DBEMT) then
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
   
   DO I = p%NumBlades+1,size(BAzimuth,1)  ! Invalid blades
      
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


subroutine calcCantAngle(f, xi,stencilSize,n,cantAngle)
! This subroutine calculates implicit cant angle based on the blade reference line that includes prebend.
    implicit none
    integer(IntKi), intent(in)  :: stencilSize, n 
    integer(IntKi)              :: i, j
    integer(IntKi)              :: sortInd(n)
    integer(IntKi)              :: info
    real(ReKi),  intent(in)     :: f(n), xi(n)
    real(ReKi)                  :: cx(stencilSize), cf(stencilSize), xiIn(stencilSize)
    real(ReKi)                  :: fIn(stencilSize), cPrime(n), fPrime(n), xiAbs(n)
    real(ReKi), intent(inout)   :: cantAngle(n)
    real(ReKi), parameter       :: ThisTol = 1e-6
    
    do i = 1,size(xi)
        
        xiAbs = abs(xi-xi(i))
        call hpsort_eps_epw (n, xiAbs, sortInd, ThisTol)
     
        if (i.eq.1) then
            fIn = f(1:stencilSize)
            xiIn = xi(1:stencilSize)
            call differ_stencil ( xi(i), 1, 2, xiIn, cx, info )
            if (info /= 0) return ! use default cantAngle in this case
            call differ_stencil ( xi(i), 1, 2, fIn, cf, info )
            if (info /= 0) return ! use default cantAngle in this case
        elseif (i.eq.size(xi)) then
            fIn = f(size(xi)-stencilSize +1:size(xi))
            xiIn = xi(size(xi)-stencilSize+1:size(xi))
            call differ_stencil ( xi(i), 1, 2, xiIn, cx, info )
            if (info /= 0) return ! use default cantAngle in this case
            call differ_stencil ( xi(i), 1, 2, fIn, cf, info )
            if (info /= 0) return ! use default cantAngle in this case
        else
            fIn = f(i-1:i+1)
            xiIn = xi(i-1:i+1)
            call differ_stencil ( xi(i), 1, 2, xiIn, cx, info )
            if (info /= 0) return ! use default cantAngle in this case
            call differ_stencil ( xi(i), 1, 2, fIn, cf, info )
            if (info /= 0) return ! use default cantAngle in this case
        endif
    
        cPrime(i) = 0.0
        fPrime(i) = 0.0

        do j = 1,size(cx)
            cPrime(i) = cPrime(i) + cx(j)*xiIn(j)
            fPrime(i) = fPrime(i) + cx(j)*fIn(j)            
        end do
        cantAngle(i) = atan2(fPrime(i),cPrime(i))*180_ReKi/pi
    end do
    
end subroutine calcCantAngle



subroutine differ_stencil ( x0, o, p, x, c, info )

!*****************************************************************************80
!
!! DIFFER_STENCIL computes finite difference coefficients.
!
!  Discussion:
!
!    We determine coefficients C to approximate the derivative at X0
!    of order O and precision P, using finite differences, so that 
!
!      d^o f(x)/dx^o (x0) = sum ( 0 <= i <= o+p-1 ) c(i) f(x(i)) 
!        + O(h^(p))
!
!    where H is the maximum spacing between X0 and any X(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 November 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X0, the point where the derivative is to 
!    be approximated.
!
!    Input, integer ( kind = 4 ) O, the order of the derivative to be 
!    approximated.  1 <= O.
!
!    Input, integer ( kind = 4 ) P, the order of the error, as a power of H.
!
!    Input, real ( kind = 8 ) X(O+P), the evaluation points.
!
!    Output, real ( kind = 8 ) C(O+P), the coefficients.
!
  implicit none

  integer(IntKi), intent(in)   :: o
  integer(IntKi), intent(in)   :: p

  real(ReKi)                   :: b(o+p)
  real(ReKi), intent(out)      :: c(o+p)
  real(ReKi)                   :: dx(o+p)
  integer(IntKi)               :: i
  integer(IntKi), intent(out)  :: info
  integer(IntKi)               :: job
  integer(IntKi)               :: n
  real(R8Ki)                   :: r8_factorial
  real(ReKi), intent(in)       :: x(o+p)
  real(ReKi), intent(in)       :: x0

  n = o + p

  dx(1:n) = x(1:n) - x0

  b(1:o+p) = 0.0D+00
  b(o+1) = 1.0D+00

  job = 0
  call r8vm_sl ( n, dx, b, c, job, info )

  if ( info /= 0 ) then
    call WrScr('DIFFER_STENCIL: Vandermonde linear system is singular.')
    return
  end if
    r8_factorial = 1.0D+00
  do i = 1,o
    r8_factorial = r8_factorial*i
  end do
  c(1:n) = c(1:n) * r8_factorial

  return
  
end subroutine differ_stencil

subroutine r8vm_sl ( n, a, b, x, job, info )

!*****************************************************************************80
!
!! R8VM_SL solves an R8VM linear system.
!
!  Discussion:
!
!    The R8VM storage format is used for an M by N Vandermonde matrix.
!    An M by N Vandermonde matrix is defined by the values in its second
!    row, which will be written here as X(1:N).  The matrix has a first 
!    row of 1's, a second row equal to X(1:N), a third row whose entries
!    are the squares of the X values, up to the M-th row whose entries
!    are the (M-1)th powers of the X values.  The matrix can be stored
!    compactly by listing just the values X(1:N).
!
!    Vandermonde systems are very close to singularity.  The singularity
!    gets worse as N increases, and as any pair of values defining
!    the matrix get close.  Even a system as small as N = 10 will
!    involve the 9th power of the defining values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 September 2003
!
!  Author:
!
!    John Burkardt.
!
!  Reference:
!
!    Gene Golub, Charles Van Loan,
!    Matrix Computations,
!    Third Edition,
!    Johns Hopkins, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns of 
!    the matrix.
!
!    Input, real ( kind = 8 ) A(N), the R8VM matrix.
!
!    Input, real ( kind = 8 ) B(N), the right hand side.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
!    Input, integer ( kind = 4 ) JOB, specifies the system to solve.
!    0, solve A * x = b.
!    nonzero, solve A' * x = b.
!
!    Output, integer ( kind = 4 ) INFO.
!    0, no error.
!    nonzero, at least two of the values in A are equal.
!
  implicit none

  integer (IntKi ), intent(in) :: n

  real(ReKi),     intent(in)   :: a(n)
  real(ReKi),     intent(in)   :: b(n)
  integer(IntKi)               :: i
  integer(IntKi), intent(out)  :: info
  integer(IntKi)               :: j
  integer(IntKi), intent(in)   :: job
  real(ReKi),     intent(out)  :: x(n)
!
!  Check for explicit singularity.
!
  info = 0

  do j = 1, n - 1
    do i = j + 1, n
      if ( a(i) == a(j) ) then
        info = 1
        return
      end if
    end do
  end do

  x(1:n) = b(1:n)

  if ( job == 0 ) then

    do j = 1, n - 1
      do i = n, j + 1, -1
        x(i) = x(i) - a(j) * x(i-1)
      end do
    end do

    do j = n - 1, 1, -1

      do i = j + 1, n
        x(i) = x(i) / ( a(i) - a(i-j) )
      end do

      do i = j, n - 1
        x(i) = x(i) - x(i+1)
      end do

    end do

  else

    do j = 1, n - 1
      do i = n, j + 1, -1
        x(i) = ( x(i) - x(i-1) ) / ( a(i) - a(i-j) )
      end do
    end do

    do j = n - 1, 1, -1
      do i = j, n - 1
        x(i) = x(i) - x(i+1) * a(j)
      end do
    end do

  end if

  return
end subroutine r8vm_sl

!                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from flib/hpsort_eps
  !---------------------------------------------------------------------
  subroutine hpsort_eps_epw (n, ra, ind, eps)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  !use kinds, ONLY : DP
  implicit none  
  !-input/output variables
  integer(IntKi), intent(in)    :: n  
  real, intent(in)        :: eps
  integer(IntKi)                :: ind (n)  
  real(ReKi)                    :: ra (n)
  !-local variables
  integer(IntKi)                :: i, ir, j, l, iind  
  real(ReKi)                    :: rra  
!
  ! initialize index array
  IF (ind (1) .eq.0) then  
     DO i = 1, n  
        ind (i) = i  
     ENDDO
  ENDIF
  ! nothing to order
  IF (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
  
    ! still in hiring phase
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    ELSE  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       IF ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  

  END DO sorting    
contains 

  !  internal function 
  !  compare two real number and return the result

  logical function hslt( a, b )
    REAL(ReKi) :: a, b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt

  !
end subroutine hpsort_eps_epw

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE AeroDyn_IO
