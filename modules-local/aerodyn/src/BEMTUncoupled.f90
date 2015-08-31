!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
module BEMTUnCoupled
 
   use NWTC_Library
   use BladeElement


   implicit none
   
   
   
   
   private
   
   public :: ComputeAirfoilCoefs
   public :: UncoupledErrFn
   public :: BEMTU_InductionWithResidual
   public :: ApplySkewedWakeCorrection
   
   public :: BEMTU_Wind
   contains
   
   
   subroutine BEMTU_Wind(phi, axInduction, tanInduction, Vx, Vy,  chord, theta, airDens, mu, AOA,  W, Re)

    
    ! in
    real(ReKi), intent(in) :: phi, axInduction, tanInduction, Vx, Vy
    real(ReKi), intent(in) :: chord, theta, airDens, mu

    ! out
    real(ReKi), intent(out) :: AOA,  Re, W
    
    ! locals
    !real(ReKi)              :: W
    
    ! angle of attack
    AOA = phi - theta 

    ! avoid numerical errors when angle is close to 0 or 90 deg
    ! and other induction factor is at some ridiculous value
    ! this only occurs when iterating on Reynolds number
    ! during the phi sweep where a solution has not been found yet
    !if ( abs(axInduction) > 10 ) then
    !    W = Vy*(1+tanInduction)/cos(phi)
    !else if ( abs(tanInduction) > 10 ) then
    !    W = Vx*(1-axInduction)/sin(phi)
    !else
        W = sqrt((Vx*(1-axInduction))**2 + (Vy*(1+tanInduction))**2)
    !end if

    Re = airDens * W * chord / mu
    if ( EqualRealNos(Re, 0.0_ReKi) ) Re = 0.001  ! Do this to avoid a singularity when we take log(Re) in the airfoil lookup.

end subroutine BEMTU_Wind

!----------------------------------------------------------------------------------------------------------------------------------  
subroutine ComputeAirfoilCoefs( phi, axInduction, tanInduction, Vx, Vy, chord, theta, airDens, mu, useAIDrag, useTIDrag, AFInfo, &
                      UA_Flag, p_UA, xd_UA, OtherState_UA, &
                      AOA, Re, Cl, Cd, Cx, Cy, Cm, errStat, errMsg )
! This routine is called from BEMTU_InductionWithResidual and possibly BEMT_CalcOutput.
! Determine the Cl, Cd, Cx, Cy coeficients for a given set of induction factors and inflow angle
!..................................................................................................................................
   real(ReKi),             intent(in   ) :: phi
   real(ReKi),             intent(in   ) :: axInduction
   real(ReKi),             intent(in   ) :: tanInduction
   real(ReKi),             intent(in   ) :: Vx
   real(ReKi),             intent(in   ) :: Vy
   real(ReKi),             intent(in   ) :: chord 
   real(ReKi),             intent(in   ) :: theta  
   real(ReKi),             intent(in   ) :: airDens
   real(ReKi),             intent(in   ) :: mu
   logical,                intent(in   ) :: useAIDrag
   logical,                intent(in   ) :: useTIDrag       
   type(AFInfoType),       intent(in   ) :: AFInfo
   logical,                intent(in   ) :: UA_Flag
   type(UA_ParameterType),       intent(in   ) :: p_UA           ! Parameters
   type(UA_DiscreteStateType),   intent(in   ) :: xd_UA          ! Discrete states at Time
   type(UA_OtherStateType),      intent(in   ) :: OtherState_UA  ! Other/optimization states
   real(ReKi),             intent(  out) :: AOA, Re, Cl, Cd, Cx, Cy, Cm
   integer(IntKi),         intent(  out) :: errStat       ! Error status of the operation
   character(*),           intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None 
   
   real(ReKi)                            :: W
      ! Compute AOA, Re, W based on current values of axInduction, tanInduction
   call BEMTU_Wind(phi, axInduction, tanInduction, Vx, Vy, chord, theta, airDens, mu, AOA, W, Re)
      
   call  BE_CalcOutputs( AFInfo, UA_Flag, AOA, W, log(Re), p_UA, xd_UA, OtherState_UA, Cl, Cd, Cm, errStat, errMsg)  
   !call  BE_CalcOutputs(AFInfo, AOA*R2D, log(Re), Cl, Cd, errStat, errMsg) ! AOA is in degrees in this look up table and Re is in log(Re)
   if (errStat >= AbortErrLev) then
      call SetErrStat( errStat, errMsg, errStat, errMsg, 'ComputeAirfoilCoefs' ) 
      return
   end if   
      
         ! Determine Cx, Cy from Cl, Cd and phi
   call BE_CalcCxCyCoefs(phi, useAIDrag, useTIDrag, Cl, Cd, Cx, Cy)
   
end subroutine ComputeAirfoilCoefs




                           ! This is the residual calculation for the uncoupled BEM solve
real(ReKi) function BEMTU_InductionWithResidual(phi, psi, chi0, numReIterations, airDens, mu, numBlades, rlocal, rtip, chord, theta,  AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, hubLossConst, tipLossConst, SkewWakeMod, &
                              UA_Flag, p_UA, xd_UA, OtherState_UA, &
                              AOA, Re, Cl, Cd, Cx, Cy, Cm, axInduction, tanInduction, chi, ErrStat, ErrMsg)
      


   real(ReKi),             intent(in   ) :: phi
   real(ReKi),             intent(in   ) :: psi
   real(ReKi),             intent(in   ) :: chi0
   integer,                intent(in   ) :: numReIterations
   real(ReKi),             intent(in   ) :: airDens
   real(ReKi),             intent(in   ) :: mu
   integer,                intent(in   ) :: numBlades
   real(ReKi),             intent(in   ) :: rlocal   
   real(ReKi),             intent(in   ) :: rtip   
   real(ReKi),             intent(in   ) :: chord 
   real(ReKi),             intent(in   ) :: theta         
   type(AFInfoType),       intent(in   ) :: AFInfo
   real(ReKi),             intent(in   ) :: Vx
   real(ReKi),             intent(in   ) :: Vy
   logical,                intent(in   ) :: useTanInd 
   logical,                intent(in   ) :: useAIDrag
   logical,                intent(in   ) :: useTIDrag
   logical,                intent(in   ) :: useHubLoss
   logical,                intent(in   ) :: useTipLoss
   real(ReKi),             intent(in   ) :: hubLossConst
   real(ReKi),             intent(in   ) :: tipLossConst
   integer,                intent(in   ) :: SkewWakeMod   ! Skewed wake model
   logical,                intent(in   ) :: UA_Flag
   type(UA_ParameterType),       intent(in   ) :: p_UA           ! Parameters
   type(UA_DiscreteStateType),   intent(in   ) :: xd_UA          ! Discrete states at Time
   type(UA_OtherStateType),      intent(in   ) :: OtherState_UA  ! Other/optimization states
   real(ReKi),             intent(  out) :: AOA, Re, Cl, Cd, Cx, Cy, Cm, axInduction, tanInduction, chi
   integer(IntKi),         intent(  out) :: ErrStat       ! Error status of the operation
   character(*),           intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
  
   ! Local variables
   
   
   real(ReKi)                            :: fzero, degAOA
   integer                               :: I, kk
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   BEMTU_InductionWithResidual = 0.0_ReKi
    
      ! Set the local version of the induction factors
   axInduction  = 0.0_ReKi  ! axInductionIN
   tanInduction = 0.0_ReKi  ! tanInductionIN
   
   
      ! If we say Re is dependent on axInduction, tanInduction, then we would create an iteration loop around the residual calculation
    
   !do I = 1,numReIterations
      
      call ComputeAirfoilCoefs( phi, axInduction, tanInduction, Vx, Vy, chord, theta, airDens, mu, useAIDrag, useTIDrag, AFInfo, &
                                .FALSE. , p_UA, xd_UA, OtherState_UA, &    ! Never use unsteady aero for this version of the airfoil coefs
                                AOA, Re, Cl, Cd, Cx, Cy, Cm, errStat, errMsg )       
      if (errStat >= AbortErrLev) then
         call SetErrStat( errStat, errMsg, errStat, errMsg, 'BEMTU_InductionWithResidual' ) 
         return
      end if
      
      if ( ( EqualRealNos(Vx, 0.0_ReKi) ) .or. ( EqualRealNos(Vy, 0.0_ReKi) ) ) then
         
         axInduction  = 0.0_ReKi
         tanInduction = 0.0_ReKi
         fzero        = 0.0_ReKi
         
      else
         
            ! Determine axInduction, tanInduction for the current Cl, Cd, phi
         call inductionFactors( rlocal, rtip, chord, phi, psi, chi0, Cx, Cy, numBlades, &
                                 Vx, Vy, useTanInd, useHubLoss, useTipLoss,  hubLossConst, tipLossConst,  SkewWakeMod, &
                                 fzero, axInduction, tanInduction, chi, errStat, errMsg)
         if (errStat >= AbortErrLev) then
            call SetErrStat( errStat, errMsg, errStat, errMsg, 'BEMTU_InductionWithResidual' ) 
            return
         end if
         
      end if
      
      BEMTU_InductionWithResidual = fzero  ! the residual
      
  ! end do
   
end function BEMTU_InductionWithResidual

      ! This is the residual calculation for the uncoupled BEM solve

real(ReKi) function UncoupledErrFn(phi, psi, chi0, numReIterations, airDens, mu, numBlades, rlocal, rtip, chord, theta, AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, hubLossConst, tipLossConst, SkewWakeMod, &
                              UA_Flag, p_UA, xd_UA, OtherState_UA, &
                              ErrStat, ErrMsg)
      


   real(ReKi),             intent(in   ) :: phi
   real(ReKi),             intent(in   ) :: psi
   real(ReKi),             intent(in   ) :: chi0
   integer,                intent(in   ) :: numReIterations
   real(ReKi),             intent(in   ) :: airDens
   real(ReKi),             intent(in   ) :: mu
   integer,                intent(in   ) :: numBlades
   real(ReKi),             intent(in   ) :: rlocal   
   real(ReKi),             intent(in   ) :: rtip   
   real(ReKi),             intent(in   ) :: chord 
   real(ReKi),             intent(in   ) :: theta         
   type(AFInfoType),       intent(in   ) :: AFInfo
   real(ReKi),             intent(in   ) :: Vx
   real(ReKi),             intent(in   ) :: Vy
   logical,                intent(in   ) :: useTanInd 
   logical,                intent(in   ) :: useAIDrag
   logical,                intent(in   ) :: useTIDrag
   logical,                intent(in   ) :: useHubLoss
   logical,                intent(in   ) :: useTipLoss
   real(ReKi),             intent(in   ) :: hubLossConst
   real(ReKi),             intent(in   ) :: tipLossConst
   integer,                intent(in   ) :: SkewWakeMod   ! Skewed wake model
   logical,                intent(in   ) :: UA_Flag
   type(UA_ParameterType),       intent(in   ) :: p_UA           ! Parameters
   type(UA_DiscreteStateType),   intent(in   ) :: xd_UA          ! Discrete states at Time
   type(UA_OtherStateType),      intent(in   ) :: OtherState_UA  ! Other/optimization states
   integer(IntKi),         intent(  out) :: ErrStat       ! Error status of the operation
   character(*),           intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
  
   ! Local variables
   
   
   real(ReKi)                            :: fzero, AOA, Re, Cl, Cd, Cx, Cy, Cm, axInduction, tanInduction, chi
   integer                               :: I
   
   ErrStat = ErrID_None
   ErrMsg  = ""
    
   
      
      UncoupledErrFn = BEMTU_InductionWithResidual(phi, psi, chi0, numReIterations, airDens, mu, numBlades, rlocal, rtip, chord, theta,  AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, hubLossConst, tipLossConst, SkewWakeMod, &
                              UA_Flag, p_UA, xd_UA, OtherState_UA, &
                              AOA, Re, Cl, Cd, Cx, Cy, Cm, axInduction, tanInduction, chi, ErrStat, ErrMsg)
      
   
   
end function UncoupledErrFn

                              
subroutine ApplySkewedWakeCorrection( Vx, Vy, azimuth, chi0, a, ap, tipRatio, phi, chi, ErrStat, ErrMsg )
   
   real(ReKi),                intent(in   ) :: Vx
   real(ReKi),                intent(in   ) :: Vy
   real(ReKi),                intent(in   ) :: azimuth
   real(ReKi),                intent(in   ) :: chi0 
   real(ReKi),                intent(inout) :: a 
   real(ReKi),                intent(inout) :: ap 
   real(ReKi),                intent(in   ) :: tipRatio            ! r/Rtip 
   real(ReKi),                intent(  out) :: phi
   real(ReKi),                intent(  out) :: chi
   integer(IntKi),            intent(  out) :: ErrStat       ! Error status of the operation
   character(*),              intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None   
   
      ! Local variables      
   real(ReKi)                               :: yawCorr, saz, x, y
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! Skewed wake correction
   
   saz = sin(azimuth)
   chi = chi0
   
   if ( abs(saz) > 0.005_ReKi ) then
      chi = (0.6_ReKi*a + 1.0_ReKi)*chi0
      
   !if (chi0 < 40.0*d2r .and. chi > 0.0 ) then
      ! TODO: Add check on chi to make sure it is < pi/2 and (positive check should be outside solve)  GJH 5/20/2015
      !yawCorr = max(0.0,chi0-0.5236)
      !yawCorr = min(0.785,yawCorr)
      yawCorr = (15.0_ReKi*pi/64.0_ReKi*tan(chi/2.0_ReKi) * (tipRatio) * saz)
               
      a = a * (1.0 +  yawCorr) ! *(-yawCorr/0.785 + 1) )
      !if ((a > 1.0 .AND. ayaw < 1.0) .OR. (a < 1.0 .AND. ayaw > 1.0 )) then
      !   call WrScr('Yaw correction crossed over 1.0.')
      !   !a = max(1.0, ayaw)
      !else if ((a < -1.0 .AND. ayaw > -1.0) .OR. (a > -1.0 .AND. ayaw < -1.0 )) then
      !   call WrScr('Yaw correction crossed over -1.0.')
      !
      !end if
         
   else
      chi = chi0
       !call WrScr('Warning: high yaw angle.  Not applying Pitt-Peters correction.')
   end if
      
   y = (1-a )*Vx
   x = (1+ap)*Vy
   
   if ( EqualRealNos(y, 0.0_ReKi) .OR. EqualRealNos(x, 0.0_ReKi) ) then
      a     = 0.0_ReKi
      ap    = 0.0_ReKi
      phi   = 0.0_ReKi
   else
      phi   = atan2(y,x)
   end if
   
   
end subroutine ApplySkewedWakeCorrection
                              
recursive subroutine inductionFactors(r , Rtip, chord, phi, azimuth, chi0, cn, ct, B, &
                              Vx, Vy, wakerotation, hubLoss, tipLoss, hubLossConst, tipLossConst, skewWakeMod, &
                              fzero, a, ap, chi, ErrStat, ErrMsg)

   implicit none

   ! in
   real(ReKi), intent(in) :: r, chord, Rtip, phi, cn, ct
   integer, intent(in) :: B
   real(ReKi), intent(in) :: Vx, Vy
   real(ReKi), intent(in) :: chi0, azimuth, hubLossConst, tipLossConst
   logical, intent(in) ::  hubLoss, tipLoss,  wakerotation
   integer, intent(in) :: skewWakeMod  ! useCd,
    
    
    
    

   ! out
   real(ReKi), intent(out) :: fzero, a, ap
   REAL(ReKi),       INTENT(  OUT) :: chi
   INTEGER(IntKi),   INTENT(  OUT) :: ErrStat       ! Error status of the operation
   CHARACTER(*),     INTENT(  OUT) :: ErrMsg        ! Error message if ErrStat /= ErrID_None

   ! local
    
      ! constants
   REAL(ReKi), PARAMETER :: c1 = 2.6296e-3 
   REAL(ReKi), PARAMETER :: c2 = 1.6222e-3 
   REAL(ReKi), PARAMETER :: c3 = -1.1111e-5 
   REAL(ReKi), PARAMETER :: c4 = -3.70371e-6 
   
   real(ReKi)  :: yawCorr
    
   
    
   real(ReKi) ::  sigma_p, sphi, cphi, lambda_r, saz !, pi
   real(ReKi) :: factortip, Ftip, factorhub, Fhub
   real(ReKi) :: k, kp,  F , ayaw  !cn, ct,
   real(ReKi) :: g1, g2, g3
   real(ReKi) :: Fsphi, sigma_pcn
   real(ReKi) :: phitemp
   !real(ReKi) :: chi

   errStat = ErrID_None
   errMsg  = ""
   
    
   if ( EqualRealNos(phi, 0.0_ReKi) ) then
      fzero =  0.0_ReKi
      a     =  0.0_ReKi
      ap    =  0.0_ReKi
      chi   =  0.0_ReKi  ! TODO: eliminate legacy return value
      return
   end if
    
   sigma_p = B/2.0_ReKi/pi*chord/r
   sphi = sin(phi)
   cphi = cos(phi)

   chi = 0.0_ReKi
   saz = sin(azimuth)
    
   !if ( EqualRealNos(azimuth, 3.141593) .OR. EqualRealNos(azimuth, -3.141593) ) then
   !   saz = 0.0_ReKi
   !else if ( EqualRealNos(azimuth, 1.570796) ) then
   !   saz = 1.0_ReKi
   !else if ( EqualRealNos(azimuth, 3*1.570796 ) .OR. EqualRealNos(azimuth, -1.570796) ) then
   !   saz = -1.0_ReKi
   !else     
   !   saz = sin(azimuth)
   !end if
    
    
      ! resolve into normal and tangential forces
      !if ( .not. useCd ) then
      !    cn = cl*cphi
      !    ct = cl*sphi
      !else
      !    cn = cl*cphi + cd*sphi
      !    ct = cl*sphi - cd*cphi
      !end if

      ! Prandtl's tip and hub loss factor
   !Ftip = 1.0_ReKi
   !! NOTE: check below isn't good enough: at tip Ftip should = 0.0, not 1.0
   !! if ( tipLoss .AND. (.NOT.(EqualRealNos(sphi, 0.0_DbKi))) ) then
   !if ( tipLoss  ) then
   !   factortip = B/2.0_ReKi*(Rtip - r)/(r*abs(sphi))
   !   Ftip      = 2.0_ReKi/pi*acos(exp(-factortip))
   !end if
   !
   !Fhub = 1.0_ReKi
   !!if ( hubLoss .AND. (.NOT.(EqualRealNos(sphi, 0.0_DbKi))) ) then
   !if ( hubLoss ) then
   !   factorhub = B/2.0_ReKi*(r - Rhub)/(Rhub*abs(sphi))
   !   Fhub      = 2.0_ReKi/pi*acos(exp(-factorhub))
   !end if
   !
   
         ! Prandtl's tip and hub loss factor
   Ftip = 1.0
   if ( tipLoss ) then
      factortip = tipLossConst/abs(sphi)
      Ftip = (2.0/pi)*acos(exp(-factortip))
   end if

   Fhub = 1.0
   if ( hubLoss ) then
      factorhub = hubLossConst/abs(sphi)
      Fhub = (2.0/pi)*acos(exp(-factorhub))
   end if
      
   F = Ftip * Fhub

   if ( EqualRealNos(F, 0.0_ReKi) ) then
      fzero =  0.0_ReKi
      a     =  1.0_ReKi
      ap    =  -1.0_ReKi
      chi   =  0.0_ReKi
      return
   end if
    
      ! bem parameters
    
      !Fsphi     = 4.0_ReKi*F*sphi**2 
      !sigma_pcn = sigma_p*cn
    
   k = sigma_p*cn/4.0_ReKi/F/sphi/sphi
    
    ! compute axial induction factor
   if (phi > 0.0_ReKi) then  ! momentum/empirical

 
        ! update axial induction factor
      if (k <= 2.0_ReKi/3.0_ReKi) then  ! momentum state
         if ( EqualRealNos( k, -1.0_ReKi) ) then
            k = k - 0.1_ReKi   ! Need to bump k to avoid singularities
         end if
      
         a = k/(1.0_ReKi+k)

      else  ! Glauert(Buhl) correction

         g1 = 2.0_ReKi*F*k - (10.0_ReKi/9-F)
         g2 = 2.0_ReKi*F*k - (4.0_ReKi/3-F)*F
         g3 = 2.0_ReKi*F*k - (25.0_ReKi/9-2*F)

         if (abs(g3) < 1e-6_ReKi) then  ! avoid singularity
               a = 1.0_ReKi - 1.0_ReKi/2.0/sqrt(g2)
         else
               a = (g1 - sqrt(g2)) / g3
         end if

      end if

   else  ! propeller brake region (a and ap not directly used but update anyway)

      if (k > 1.0_ReKi) then
      !if (sigma_pcn > Fsphi) then
         a =   k/(k-1.0_ReKi) !sigma_pcn / (sigma_pcn - Fsphi )  !
      else
         a = 0.0_ReKi  ! dummy value
      end if

   end if

    ! apply yaw correction
    !if (skewWakeMod) then
    !    
    !    chi = (0.6*a + 1.0)*chi0
    !    a = a * (1.0 + 15.0*pi/32*tan(chi/2.0) * r/Rtip * saz)
    !    a = min(a, 0.999999)
    !end if

    
   
    ! compute tangential induction factor
   kp = sigma_p*ct/4.0_ReKi/F/sphi/cphi
      ! Per conversation with Rick, we should only trigger this if phi = 0 , so we will return predefined values as if phi=0.0
   if (EqualRealNos(kp, 1.0_ReKi)) then
      fzero =  0.0_ReKi
      a     =  0.0_ReKi
      ap    =  0.0_ReKi
      chi   =  0.0_ReKi
      return
   end if
   
   ap = kp/(1.0_ReKi-kp)
!bjj: 3-jun-2015: TODO: was able to trigger divide-by-zero here using ccBlade_UAE.dvr without tiploss or hubloss
    
   if (.not. wakerotation) then
      ap = 0.0_ReKi
      kp = 0.0_ReKi
   end if

    !if ( skewWakeMod > SkewMod_Uncoupled ) then  
    !  phitemp = InflowAngle(Vx_in, Vy_in, REAL(a, ReKi), REAL(ap))
    !  call inductionFactors(r_in     , Rtip_in, chord_in, Rhub_in,  lambda_in, phitemp, azimuth_in, yaw_in  , cn_in, ct_in, B, &
    !                          Vx_in, Vy_in, wakerotation,   hubLoss , tipLoss   , 0, &
    !                          fzero_out, a_out,           ap_out,           chi_out, ErrStat, ErrMsg)
    !  return
    !end if  
    
    ! error function
   lambda_r = Vy/Vx
   if (phi > 0) then  ! momentum/empirical
      if ( EqualRealNos(a, 1.0_ReKi) ) then
         fzero = - cphi/lambda_r*(1-kp)
      else       
         fzero = sphi/(1-a) - cphi/lambda_r*(1-kp)
      end if
      
   else  ! propeller brake region
      fzero = sphi*(1-k) - cphi/lambda_r*(1-kp)
   end if
    
end subroutine inductionFactors
    

                              

                              
                              
end module BEMTUncoupled