!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
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
   use AirfoilInfo_Types
   use UnsteadyAero
   use UnsteadyAero_Types


   implicit none
   
   integer(IntKi), public, parameter  :: SkewMod_Uncoupled  = 1      ! Uncoupled (no correction) [-]
   integer(IntKi), public, parameter  :: SkewMod_PittPeters = 2      ! Pitt/Peters [-]
   integer(IntKi), public, parameter  :: SkewMod_Coupled    = 3      ! Coupled [-]
   
   
   private
   
   public :: Compute_UA_AirfoilCoefs
   public :: ComputeSteadyAirfoilCoefs
   public :: UncoupledErrFn
   public :: BEMTU_InductionWithResidual
   public :: ApplySkewedWakeCorrection
   public :: Transform_ClCd_to_CxCy
   
   public :: BEMTU_Wind
   contains
   
   
   subroutine BEMTU_Wind( axInduction, tanInduction, Vx, Vy,  chord, airDens, mu, W, Re )

    
    ! in
    real(ReKi), intent(in) :: axInduction, tanInduction, Vx, Vy
    real(ReKi), intent(in) :: chord, airDens, mu

    ! out
    real(ReKi), intent(out) :: Re, W
    
    
    

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

subroutine Transform_ClCd_to_CxCy( phi, useAIDrag, useTIDrag, Cl, Cd, Cx, Cy )
   real(ReKi),             intent(in   ) :: phi
   logical,                intent(in   ) :: useAIDrag
   logical,                intent(in   ) :: useTIDrag
   real(ReKi),             intent(in   ) :: Cl
   real(ReKi),             intent(in   ) :: Cd
   real(ReKi),             intent(  out) :: Cx
   real(ReKi),             intent(  out) :: Cy

   real(ReKi)      cphi, sphi

   cphi = cos(phi)
   sphi = sin(phi)
   
      ! resolve into normal (x) and tangential (y) forces
   if (  useAIDrag ) then
      Cx = Cl*cphi + Cd*sphi
   else      
      Cx = Cl*cphi
   end if
    
   if (  useTIDrag ) then     
      Cy = Cl*sphi - Cd*cphi
   else     
      Cy = Cl*sphi
   end if
   
end subroutine Transform_ClCd_to_CxCy

!----------------------------------------------------------------------------------------------------------------------------------  
subroutine ComputeSteadyAirfoilCoefs( AOA, Re, AFInfo, &
                      Cl, Cd, Cm, errStat, errMsg )
! This routine is called from BEMTU_InductionWithResidual and possibly BEMT_CalcOutput.
! Determine the Cl, Cd, Cm, coeficients for a given angle of attack
!..................................................................................................................................
   real(ReKi),             intent(in   ) :: AOA
   real(ReKi),             intent(in   ) :: Re           ! Unused in the current version!     
   type(AFInfoType),       intent(in   ) :: AFInfo
   real(ReKi),             intent(  out) :: Cl, Cd, Cm
   integer(IntKi),         intent(  out) :: errStat       ! Error status of the operation
   character(*),           intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None 
   

   real                            :: IntAFCoefs(4)                ! The interpolated airfoil coefficients.
   integer                         :: s1      
      
   ErrStat = ErrID_None
   ErrMsg  = ''
   IntAFCoefs = 0.0_ReKi ! initialize in case we only don't have 4 columns in the airfoil data (i.e., so cm is zero if not in the file)
      
      
    
      
   
      ! NOTE: we use Table(1) because the right now we can only interpolate with AOA and not Re or other variables.  If we had multiple tables stored
      ! for changes in other variables (Re, Mach #, etc) then then we would need to interpolate across tables.
      !
   s1 = size(AFInfo%Table(1)%Coefs,2)
   
   IntAFCoefs(1:s1) = CubicSplineInterpM( 1.0_ReKi*real( AOA*R2D, ReKi ) &
                                          , AFInfo%Table(1)%Alpha &
                                          , AFInfo%Table(1)%Coefs &
                                          , AFInfo%Table(1)%SplineCoefs &
                                          , ErrStat, ErrMsg )
   
   Cl = IntAFCoefs(1)
   Cd = IntAFCoefs(2)
   Cm = IntAFCoefs(3)
     
   
       
end subroutine ComputeSteadyAirfoilCoefs



!----------------------------------------------------------------------------------------------------------------------------------  
subroutine Compute_UA_AirfoilCoefs( AOA, U, Re, AFInfo, &
                      p_UA, xd_UA, OtherState_UA, OtherState_y_UA, &
                      Cl, Cd, Cm, errStat, errMsg )
! This routine is called from BEMTU_InductionWithResidual and possibly BEMT_CalcOutput.
! Determine the Cl, Cd, Cm coeficients for a given angle of attack
!..................................................................................................................................
   real(ReKi),                   intent(in   ) :: AOA
   real(ReKi),                   intent(in   ) :: U
   real(ReKi),                   intent(in   ) :: Re                 ! Unused in the current version!
   type(AFInfoType),             intent(in   ) :: AFInfo
   type(UA_ParameterType),       intent(in   ) :: p_UA               ! Parameters
   type(UA_DiscreteStateType),   intent(in   ) :: xd_UA              ! Discrete states at Time
   type(UA_OtherStateType),      intent(inout) :: OtherState_UA      ! Other/optimization states
   type(UA_OutputType),          intent(inout) :: OtherState_y_UA    !
   real(ReKi),                   intent(  out) :: Cl, Cd, Cm
   integer(IntKi),               intent(  out) :: errStat            ! Error status of the operation
   character(*),                 intent(  out) :: errMsg             ! Error message if ErrStat /= ErrID_None 
   
   integer(intKi)                              :: ErrStat2           ! temporary Error status
   character(ErrMsgLen)                        :: ErrMsg2            ! temporary Error message
   character(*), parameter                     :: RoutineName = 'Compute_UA_AirfoilCoefs'
  
   
   type(UA_InputType)              :: u_UA
   type(UA_OutputType)             :: y_UA          !
      
   ErrStat = ErrID_None
   ErrMsg  = ''

   u_UA%alpha = AOA   
   u_UA%Re    = Re
   u_UA%U     = U
   
   !bjj: TODO: this gets called element-by-element (not all at once). Are OtherState%iBladeNode and OtherState%iBlade set properly?
#ifdef DEBUG_v14
   call UA_CalcOutput2(u_UA, p_UA, xd_UA, OtherState_UA, AFInfo, OtherState_y_UA, errStat2, errMsg2 )
#else
   call UA_CalcOutput(u_UA, p_UA, xd_UA, OtherState_UA, AFInfo, OtherState_y_UA, errStat2, errMsg2 )
#endif
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName ) 
      if (errStat >= AbortErrLev) return

   Cl         = OtherState_y_UA%Cl
   Cd         = OtherState_y_UA%Cd
   Cm         = OtherState_y_UA%Cm
                  
       
end subroutine Compute_UA_AirfoilCoefs

                           ! This is the residual calculation for the uncoupled BEM solve
real(ReKi) function BEMTU_InductionWithResidual(phi, AOA, Re, numBlades, rlocal, rtip, chord, AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, hubLossConst, tipLossConst,  &
                              axInduction, tanInduction,  ErrStat, ErrMsg)
      


   real(ReKi),             intent(in   ) :: phi
   real(ReKi),             intent(in   ) :: AOA
   real(ReKi),             intent(in   ) :: Re
   integer,                intent(in   ) :: numBlades
   real(ReKi),             intent(in   ) :: rlocal   
   real(ReKi),             intent(in   ) :: rtip   
   real(ReKi),             intent(in   ) :: chord         
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
   real(ReKi),             intent(  out) :: axInduction, tanInduction
   integer(IntKi),         intent(  out) :: ErrStat       ! Error status of the operation
   character(*),           intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
  
   ! Local variables
   
   integer(intKi)                        :: ErrStat2           ! temporary Error status
   character(ErrMsgLen)                  :: ErrMsg2            ! temporary Error message
   character(*), parameter               :: RoutineName = 'BEMTU_InductionWithResidual'
   
   real(ReKi)                            :: fzero

   real(ReKi)                            :: Cl, Cd, Cx, Cy, Cm
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   BEMTU_InductionWithResidual = 0.0_ReKi
    
      ! Set the local version of the induction factors
   axInduction  = 0.0_ReKi  ! axInductionIN
   tanInduction = 0.0_ReKi  ! tanInductionIN
   
         
   if (( .NOT. EqualRealNos(Vx, 0.0_ReKi) ) .AND. ( .NOT. EqualRealNos(Vy, 0.0_ReKi) ) ) then 
      
      call ComputeSteadyAirfoilCoefs( AOA, Re, AFInfo, Cl, Cd, Cm, errStat2, errMsg2 )       
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName ) 
         if (ErrStat >= AbortErrLev) return
      
         ! Compute Cx, Cy given Cl, Cd and phi, we honor the useAIDrag and useTIDrag flag because Cx,Cy are only used for the solution of inductions
      call Transform_ClCd_to_CxCy( phi, useAIDrag, useTIDrag, Cl, Cd, Cx, Cy )  
      
         ! Determine axInduction, tanInduction for the current Cl, Cd, phi
      call inductionFactors( rlocal, rtip, chord, phi, Cx, Cy, numBlades, &
                              Vx, Vy, useTanInd, useHubLoss, useTipLoss,  hubLossConst, tipLossConst, &
                              fzero, axInduction, tanInduction, errStat2, errMsg2)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName ) 
         if (ErrStat >= AbortErrLev) return
      BEMTU_InductionWithResidual = fzero  ! the residual
      
   end if
      
   
end function BEMTU_InductionWithResidual

      ! This is the residual calculation for the uncoupled BEM solve

real(ReKi) function UncoupledErrFn(phi, AOA, Re, numBlades, rlocal, rtip, chord, AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, hubLossConst, tipLossConst, &
                              ErrStat, ErrMsg)
      


   real(ReKi),             intent(in   ) :: phi
   real(ReKi),             intent(in   ) :: AOA
   real(ReKi),             intent(in   ) :: Re
   integer,                intent(in   ) :: numBlades
   real(ReKi),             intent(in   ) :: rlocal   
   real(ReKi),             intent(in   ) :: rtip   
   real(ReKi),             intent(in   ) :: chord         
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
   integer(IntKi),         intent(  out) :: ErrStat       ! Error status of the operation
   character(*),           intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
  
   ! Local variables
   
   
   real(ReKi)                            :: axInduction, tanInduction
   
   ErrStat = ErrID_None
   ErrMsg  = ""
    
   UncoupledErrFn = BEMTU_InductionWithResidual(phi, AOA, Re, numBlades, rlocal, rtip, chord, AFInfo, &
                           Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, hubLossConst, tipLossConst, &
                           axInduction, tanInduction, ErrStat, ErrMsg)

   
end function UncoupledErrFn

                              
subroutine ApplySkewedWakeCorrection( Vx, Vy, azimuth, chi0, tipRatio, a, ap, chi, ErrStat, ErrMsg )
   
   real(ReKi),                intent(in   ) :: Vx
   real(ReKi),                intent(in   ) :: Vy
   real(ReKi),                intent(in   ) :: azimuth
   real(ReKi),                intent(in   ) :: chi0 
   real(ReKi),                intent(in   ) :: tipRatio            ! r/Rtip 
   real(ReKi),                intent(inout) :: a 
   real(ReKi),                intent(inout) :: ap 
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
      !bjj: modified 22-Sep-2015: RRD recommends 32 instead of 64 in the denominator (like AD14)
      yawCorr = (15.0_ReKi*pi/32.0_ReKi*tan(chi/2.0_ReKi) * (tipRatio) * saz)
               
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
   end if
      
   y = (1-a )*Vx
   x = (1+ap)*Vy
   
   if ( EqualRealNos(y, 0.0_ReKi) .OR. EqualRealNos(x, 0.0_ReKi) ) then
      a     = 0.0_ReKi
      ap    = 0.0_ReKi

   end if
   
end subroutine ApplySkewedWakeCorrection
                              
recursive subroutine inductionFactors(r , Rtip, chord, phi, cn, ct, B, &
                              Vx, Vy, wakerotation, hubLoss, tipLoss, hubLossConst, tipLossConst, &
                              fzero, a, ap, ErrStat, ErrMsg)

   implicit none

   ! in
   real(ReKi), intent(in) :: r, chord, Rtip, phi, cn, ct
   integer, intent(in) :: B
   real(ReKi), intent(in) :: Vx, Vy
   real(ReKi), intent(in) :: hubLossConst, tipLossConst
   logical, intent(in) ::  hubLoss, tipLoss,  wakerotation
   
    
    
    
    

   ! out
   real(ReKi), intent(out) :: fzero, a, ap
   INTEGER(IntKi),   INTENT(  OUT) :: ErrStat       ! Error status of the operation
   CHARACTER(*),     INTENT(  OUT) :: ErrMsg        ! Error message if ErrStat /= ErrID_None

   ! local
    
      ! constants
   REAL(ReKi), PARAMETER :: c1 = 2.6296e-3 
   REAL(ReKi), PARAMETER :: c2 = 1.6222e-3 
   REAL(ReKi), PARAMETER :: c3 = -1.1111e-5 
   REAL(ReKi), PARAMETER :: c4 = -3.70371e-6 
   

    
   
    
   real(ReKi) ::  sigma_p, sphi, cphi, lambda_r, saz !, pi
   real(ReKi) :: factortip, Ftip, factorhub, Fhub
   real(ReKi) :: k, kp,  F 
   real(ReKi) :: g1, g2, g3

 

   errStat = ErrID_None
   errMsg  = ""
   
      ! We are simply going to bail if we are using tiploss and tipLossConst = 0 or using hubloss and hubLossConst=0, regardless of phi!
   if ( ( tiploss .and. EqualRealNos(tipLossConst,0.0_ReKi) ) .or. ( hubloss .and. EqualRealNos(hubLossConst,0.0_ReKi) ) ) then
      fzero =  0.0_ReKi
      a     =  1.0_ReKi
      ap    =  -1.0_ReKi
      return      
   end if
   
   
   if ( EqualRealNos(phi, 0.0_ReKi) ) then
      fzero =  0.0_ReKi
      a     =  0.0_ReKi
      ap    =  0.0_ReKi
      return
   end if
    
   sigma_p = B/2.0_ReKi/pi*chord/r
   sphi = sin(phi)
   cphi = cos(phi)
   
    
   
    
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
      Ftip = (2.0/pi)*acos(min(1.0_ReKi,exp(-factortip)))
   end if

   Fhub = 1.0
   if ( hubLoss ) then
      factorhub = hubLossConst/abs(sphi)
      Fhub = (2.0/pi)*acos(min(1.0_ReKi,exp(-factorhub)))
   end if
      
   F = Ftip * Fhub

   
    
      ! bem parameters
 
    
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

   else  ! propeller brake region (a and ap not directly used but update anyway) !bjj: huh? when k is slightly larger than 1, a is definitely getting used (and causing issues)...

      if (k > 1.0_ReKi .and. .not. EqualRealNos(k, 1.0_ReKi) ) then
      !if (sigma_pcn > Fsphi) then
         a =   k/(k-1.0_ReKi) !sigma_pcn / (sigma_pcn - Fsphi )  !

         ! axial induction is blowing up, so I'm putting a band-aid here. BJJ 25-Feb-2016
         a = min(a, 10.0_ReKi ) 
      
      else
         a = 0.0_ReKi  ! dummy value
      end if

   end if

    
   
    ! compute tangential induction factor
   if ( cphi==0.0_ReKi ) then ! We don't want NaN here
      kp = HUGE(kp)
   else
      kp = sigma_p*ct/4.0_ReKi/F/sphi/cphi
   end if
   
      ! Per conversation with Rick, we should only trigger this if phi = 0 , so we will return predefined values as if phi=0.0
   if (EqualRealNos(kp, 1.0_ReKi)) then
      fzero =  0.0_ReKi
      a     =  0.0_ReKi
      ap    =  0.0_ReKi
      return
   end if
   
   ap = kp/(1.0_ReKi-kp)
   ! tangential induction is blowing up, so we're putting a band-aid here. GJH, JMJ, BJJ 1-Sep-2015
   if ( abs(ap) > 10.0_ReKi ) then
      ap = sign( 10.0_ReKi, ap )
   end if
      
   
!bjj: 3-jun-2015: TODO: was able to trigger divide-by-zero here using ccBlade_UAE.dvr without tiploss or hubloss
    
   if (.not. wakerotation) then
      ap = 0.0_ReKi
      kp = 0.0_ReKi
   end if

    
    
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