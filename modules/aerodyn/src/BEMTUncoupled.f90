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
module BEMTUnCoupled
 
   use NWTC_Library
   use AirfoilInfo_Types
   use UnsteadyAero
   use UnsteadyAero_Types


   implicit none
   
   integer(IntKi), public, parameter  :: SkewMod_Uncoupled  = 1            ! Uncoupled (no correction) [-]
   integer(IntKi), public, parameter  :: SkewMod_PittPeters = 2            ! Pitt/Peters [-]
   integer(IntKi), public, parameter  :: SkewMod_Coupled    = 3            ! Coupled [-]
   
   real(ReKi),     public, parameter  :: BEMT_MaxInduction(2) = (/1.5_ReKi, 1.0_ReKi /)  ! largest magnitude of axial (1) and tangential (2) induction factors
   real(ReKi),     public, parameter  :: BEMT_MinInduction(2) = -1.0_ReKi

   
   !1e-6 works for double precision, but not single precision 
   real(ReKi),     public, parameter  :: BEMT_epsilon2 = 10.0_ReKi*sqrt(epsilon(1.0_ReKi)) !this is the tolerance in radians for values around singularities in phi (i.e., phi=0 and phi=pi/2); must be large enough so that EqualRealNos(BEMT_epsilon2, 0.0_ReKi) is false
   
   
   private
   
   public :: Compute_UA_AirfoilCoefs
   public :: ComputeSteadyAirfoilCoefs
   public :: UncoupledErrFn
   public :: BEMTU_InductionWithResidual
   public :: ApplySkewedWakeCorrection
   public :: Transform_ClCd_to_CxCy
   
   public :: BEMTU_Wind
   public :: VelocityIsZero
contains
   
!..................................................................................................................................   
   function VelocityIsZero ( v )

      ! passed variables

   REAL(ReKi), INTENT(IN )         :: v                                 !< the velocity that needs to be compared with zero

   LOGICAL                         :: VelocityIsZero                    !< .true. if and only if the velocity is (almost) equal to zero

   
      
      VelocityIsZero = abs(v) < 0.001_ReKi ! tolerance in m/s for what we consider zero velocity for BEM computations
   
   end function VelocityIsZero
!..................................................................................................................................   
   
   subroutine BEMTU_Wind( axInduction, tanInduction, Vx, Vy,  chord, nu, W, Re )

    
    ! in
    real(ReKi), intent(in) :: axInduction, tanInduction, Vx, Vy
    real(ReKi), intent(in) :: chord, nu

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

    Re = W * chord / nu
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
                      Cl, Cd, Cm, Cpmin, errStat, errMsg )
! This routine is called from BEMTU_InductionWithResidual and possibly BEMT_CalcOutput.
! Determine the Cl, Cd, Cm, coeficients for a given angle of attack
!..................................................................................................................................
   real(ReKi),             intent(in   ) :: AOA
   real(ReKi),             intent(in   ) :: Re           ! Unused in the current version!     
   type(AFInfoType),       intent(in   ) :: AFInfo
   real(ReKi),             intent(  out) :: Cl, Cd, Cm, Cpmin
   integer(IntKi),         intent(  out) :: errStat       ! Error status of the operation
   character(*),           intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None 
   

   real                            :: IntAFCoefs(4)                ! The interpolated airfoil coefficients.
   real(reki)                      :: Alpha
   integer                         :: s1      
      
   ErrStat = ErrID_None
   ErrMsg  = ''
   IntAFCoefs = 0.0_ReKi ! initialize in case we only don't have 4 columns in the airfoil data (i.e., so cm is zero if not in the file)
      
      
    
      
   
      ! NOTE: we use Table(1) because the right now we can only interpolate with AOA and not Re or other variables.  If we had multiple tables stored
      ! for changes in other variables (Re, Mach #, etc) then then we would need to interpolate across tables.
      !
   s1 = size(AFInfo%Table(1)%Coefs,2)
   
   Alpha = AOA
   call MPi2Pi ( Alpha ) ! change AOA into range of -pi to pi
   IntAFCoefs(1:s1) = CubicSplineInterpM( Alpha  &
                                          , AFInfo%Table(1)%Alpha &
                                          , AFInfo%Table(1)%Coefs &
                                          , AFInfo%Table(1)%SplineCoefs &
                                          , ErrStat, ErrMsg )
   
  
   Cl    = IntAFCoefs(1)
   Cd    = IntAFCoefs(2)
   Cm    = 0.0_Reki  !Set these to zero unless there is data to be read in
   Cpmin = 0.0_Reki
     
   IF ( AFInfo%ColCm > 0 ) Cm = IntAFCoefs(AFInfo%ColCm)
         
   IF ( AFInfo%ColCpmin > 0 ) Cpmin = IntAFCoefs(AFInfo%ColCpmin)
      
             
end subroutine ComputeSteadyAirfoilCoefs
   
!----------------------------------------------------------------------------------------------------------------------------------  
subroutine Compute_UA_AirfoilCoefs( AOA, U, Re, AFInfo, &
                      p_UA, xd_UA, OtherState_UA, y_UA, m_UA, &
                      Cl, Cd, Cm, errStat, errMsg )
! This routine is called from BEMTU_InductionWithResidual and possibly BEMT_CalcOutput.
! Determine the Cl, Cd, Cm coeficients for a given angle of attack
!..................................................................................................................................
   real(ReKi),                   intent(in   ) :: AOA                !< angle of attack, radians
   real(ReKi),                   intent(in   ) :: U                  !< Vrel, m/s
   real(ReKi),                   intent(in   ) :: Re                 ! Unused in the current version!
   type(AFInfoType),             intent(in   ) :: AFInfo
   type(UA_ParameterType),       intent(in   ) :: p_UA               ! Parameters
   type(UA_DiscreteStateType),   intent(in   ) :: xd_UA              ! Discrete states at Time
   type(UA_OtherStateType),      intent(in   ) :: OtherState_UA      ! Other states at Time
   type(UA_OutputType),          intent(inout) :: y_UA               !
   type(UA_MiscVarType),         intent(inout) :: m_UA               ! misc/optimization variables
   real(ReKi),                   intent(  out) :: Cl, Cd, Cm
   integer(IntKi),               intent(  out) :: errStat            ! Error status of the operation
   character(*),                 intent(  out) :: errMsg             ! Error message if ErrStat /= ErrID_None 
   
   integer(intKi)                              :: ErrStat2           ! temporary Error status
   character(ErrMsgLen)                        :: ErrMsg2            ! temporary Error message
   character(*), parameter                     :: RoutineName = 'Compute_UA_AirfoilCoefs'
  
   
   type(UA_InputType)              :: u_UA
      
   ErrStat = ErrID_None
   ErrMsg  = ''

   u_UA%alpha = AOA   
   u_UA%Re    = Re
   u_UA%U     = U
   
   call UA_CalcOutput(u_UA, p_UA, xd_UA, OtherState_UA, AFInfo, y_UA, m_UA, errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName ) 
      if (errStat >= AbortErrLev) return

   Cl         = y_UA%Cl
   Cd         = y_UA%Cd
   Cm         = y_UA%Cm
                  
       
end subroutine Compute_UA_AirfoilCoefs
!----------------------------------------------------------------------------------------------------------------------------------
!>This is the residual calculation for the uncoupled BEM solve
real(ReKi) function BEMTU_InductionWithResidual(phi, AOA, Re, numBlades, rlocal, chord, AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, hubLossConst, tipLossConst,  &
                              axInduction, tanInduction,  IsValidSolution, ErrStat, ErrMsg)
      


   real(ReKi),             intent(in   ) :: phi
   real(ReKi),             intent(in   ) :: AOA
   real(ReKi),             intent(in   ) :: Re
   integer,                intent(in   ) :: numBlades
   real(ReKi),             intent(in   ) :: rlocal      
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
   logical,                intent(  out) :: IsValidSolution !< this is set to false if k<=1 in the propeller brake region or k<-1 in the momentum region, indicating an invalid solution
   integer(IntKi),         intent(  out) :: ErrStat       ! Error status of the operation
   character(*),           intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
  
   ! Local variables
   
   integer(intKi)                        :: ErrStat2           ! temporary Error status
   character(ErrMsgLen)                  :: ErrMsg2            ! temporary Error message
   character(*), parameter               :: RoutineName = 'BEMTU_InductionWithResidual'
   
   real(ReKi)                            :: fzero

   real(ReKi)                            :: Cl, Cd, Cx, Cy, Cm, Cpmin
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   BEMTU_InductionWithResidual = 0.0_ReKi
   IsValidSolution = .true.
   
   ! make these return values consistent with what is returned in inductionFactors routine:
    
      ! Set the local version of the induction factors
   if ( ( useTiploss .and. EqualRealNos(tipLossConst,0.0_ReKi) ) .or. ( useHubloss .and. EqualRealNos(hubLossConst,0.0_ReKi) ) ) then
      ! We are simply going to bail if we are using tiploss and tipLossConst = 0 or using hubloss and hubLossConst=0, regardless of phi! [do this before checking if Vx or Vy is zero or you'll get jumps in the induction and loads]
      axInduction  =  1.0_ReKi
      tanInduction =  0.0_ReKi
   elseif ( EqualRealNos(phi, 0.0_ReKi) .or. VelocityIsZero(Vx) .OR. VelocityIsZero(Vy) ) then 
      axInduction  =  0.0_ReKi
      tanInduction =  0.0_ReKi
   else !if ( (.NOT. VelocityIsZero(Vx)) .AND. (.NOT. VelocityIsZero(Vy)) ) then 

      call ComputeSteadyAirfoilCoefs( AOA, Re, AFInfo, Cl, Cd, Cm, Cpmin, errStat2, errMsg2 )       !bjj: would be nice if this could be done outside this routine (so we don't copy AFInfo so much)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName ) 
         if (ErrStat >= AbortErrLev) return
      
         ! Compute Cx, Cy given Cl, Cd and phi, we honor the useAIDrag and useTIDrag flag because Cx,Cy are only used for the solution of inductions
      call Transform_ClCd_to_CxCy( phi, useAIDrag, useTIDrag, Cl, Cd, Cx, Cy )  
      
      
         ! Determine axInduction, tanInduction for the current Cl, Cd, phi
      call inductionFactors( rlocal, chord, phi, Cx, Cy, numBlades, &
                              Vx, Vy, useTanInd, useHubLoss, useTipLoss,  hubLossConst, tipLossConst, &
                              fzero, axInduction, tanInduction, IsValidSolution)
      BEMTU_InductionWithResidual = fzero  ! the residual
      
   end if
      
   
end function BEMTU_InductionWithResidual


      ! This is the residual calculation for the uncoupled BEM solve

real(ReKi) function UncoupledErrFn(phi, theta, Re, numBlades, rlocal, chord, AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, hubLossConst, tipLossConst, &
                              IsValidSolution, ErrStat, ErrMsg)
      


   real(ReKi),             intent(in   ) :: phi
   real(ReKi),             intent(in   ) :: theta
   real(ReKi),             intent(in   ) :: Re
   integer,                intent(in   ) :: numBlades
   real(ReKi),             intent(in   ) :: rlocal      
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
   logical,                intent(  out) :: IsValidSolution !< this is set to false if k<=1 in the propeller brake region or k<-1 in the momentum region, indicating an invalid solution
   character(*),           intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
  
   ! Local variables
   
   real(ReKi)                            :: axInduction, tanInduction, AoA
   
   ErrStat = ErrID_None
   ErrMsg  = ""
    
   AOA = phi - theta
   
                
   UncoupledErrFn = BEMTU_InductionWithResidual(phi, AOA, Re, numBlades, rlocal, chord, AFInfo, &
                           Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, hubLossConst, tipLossConst, &
                           axInduction, tanInduction, IsValidSolution, ErrStat, ErrMsg)

   
end function UncoupledErrFn

                              
subroutine ApplySkewedWakeCorrection( yawCorrFactor, Vx, Vy, azimuth, chi0, tipRatio, a, ap, chi, FirstWarn )
   
   real(ReKi),                intent(in   ) :: yawCorrFactor ! set to 15*pi/32 previously; now allowed to be input (to better match data) 
   real(ReKi),                intent(in   ) :: Vx
   real(ReKi),                intent(in   ) :: Vy
   real(ReKi),                intent(in   ) :: azimuth
   real(ReKi),                intent(in   ) :: chi0 
   real(ReKi),                intent(in   ) :: tipRatio            ! r/Rtip 
   real(ReKi),                intent(inout) :: a 
   real(ReKi),                intent(inout) :: ap 
   real(ReKi),                intent(  out) :: chi
   logical(IntKi),            intent(inout) :: FirstWarn       ! If this is the first warning about invalid skew
   
      ! Local variables      
   real(ReKi)                               :: yawCorr
   real(ReKi)                               :: yawCorr_tan ! magnitude of the tan(chi/2) correction term (with possible limits)
   
   
   ! Skewed wake correction
      
   chi = (0.6_ReKi*a + 1.0_ReKi)*chi0
      
   call MPi2Pi( chi ) ! make sure chi is in [-pi, pi] before testing if it's outside a valid range
      
   if (abs(chi) > piBy2) then
         
      if (FirstWarn) then
         call WrScr( 'Warning: SkewedWakeCorrection encountered a large value of chi ('//trim(num2lstr(chi*R2D))// &
            ' deg), so the yaw correction will be limited. This warning will not be repeated though the condition may persist. See the AD15 chi output channels, and'// &
            ' consider turning off the Pitt/Peters skew model (set SkewMod=1) if this condition persists.')
         FirstWarn = .false.
      end if
         
      yawCorr_tan = sign( 1.0_ReKi, chi ) ! set to +/- 1 = +/- tan( pi/4 )
   else
      yawCorr_tan = tan(chi/2.0_ReKi)
   end if
      
      !bjj: modified 22-Sep-2015: RRD recommends 32 instead of 64 in the denominator (like AD14)
   yawCorr = ( yawCorrFactor * yawCorr_tan * (tipRatio) * sin(azimuth) ) ! bjj: note that when chi gets close to +/-pi this blows up
      
   a = a * (1.0 +  yawCorr)
   
   
end subroutine ApplySkewedWakeCorrection
!-----------------------------------------------------------------------------------------
!> This subroutine computes the induction factors (a) and (ap) along with the residual (fzero)
subroutine inductionFactors(r, chord, phi, cn, ct, B, Vx, Vy, wakerotation, useHubLoss, useTipLoss, hubLossConst, tipLossConst, &
                              fzero, a, ap, IsValidSolution)

   implicit none

   ! in
   real(ReKi), intent(in) :: r              !< local radial position [u%rlocal]
   real(ReKi), intent(in) :: chord          !< chord [p%chord]
   real(ReKi), intent(in) :: phi            !< angle between the plane of rotation and the direction of the local wind [y%phi]; must be in range [-pi,pi]
   real(ReKi), intent(in) :: cn             !< normal force coefficient (normal to the plane, not chord) of the jth node in the kth blade; [y%cx]
   real(ReKi), intent(in) :: ct             !< tangential force coefficient (tangential to the plane, not chord) of the jth node in the kth blade; [y%cy]
   integer,    intent(in) :: B              !< number of blades [p%numBlades]
   real(ReKi), intent(in) :: Vx             !< velocity component [u%Vx]
   real(ReKi), intent(in) :: Vy             !< velocity component [u%Vy]
   real(ReKi), intent(in) :: hubLossConst   !< hub loss constant [p%hubLossConst]
   real(ReKi), intent(in) :: tipLossConst   !< tip loss constant [p%tipLossConst]
   logical,    intent(in) :: useHubLoss     !< hub-loss flag [p%useHubLoss]
   logical,    intent(in) :: useTipLoss     !< tip-loss flag [p%useTipLoss]
   logical,    intent(in) :: wakerotation   !< Include tangential induction in BEMT calculations [flag] [p%useTanInd]
                   
   ! out
   real(ReKi), intent(out) :: fzero         !< residual of BEM equations
   real(ReKi), intent(out) :: a             !< axial induction [y%axInduction]
   real(ReKi), intent(out) :: ap            !< tangential induction, i.e., a-prime [y%tanInduction]
   logical,    intent(out) :: IsValidSolution !< this is set to false if k<=1 in the propeller brake region or k<-1 in the momentum region, indicating an invalid solution
   
   ! local
        
   real(ReKi) :: sigma_p   ! local solidity (B*chord/(TwoPi*r))
   real(ReKi) :: sphi, cphi, lambda_r
   real(ReKi) :: k, kp ! non-dimensional parameters 
   real(ReKi) :: F ! hub/tip loss correction factor
   real(ReKi) :: g1, g2, g3
   real(ReKi) :: temp  ! temporary variable so we don't have to calculate 2.0_ReKi*F*k multiple times
   real(ReKi), parameter :: InductionLimit = 1000000.0_ReKi
   real(ReKi), parameter :: MaxTnInd = BEMT_MaxInduction(2)
   real(ReKi), parameter :: MaxAxInd = BEMT_MaxInduction(1)
   real(ReKi), parameter :: MinTnInd = BEMT_MinInduction(2)
   real(ReKi), parameter :: MinAxInd = BEMT_MinInduction(1)
   
   logical    :: momentumRegion

   
   IsValidSolution  = .true.
   
   !.....................................................
   ! Some special cases (bjj commented out because we have taken care of these in BEMTU_InductionWithResidual, the only routine that calls this function)
   !.....................................................
   !if ( ( useTiploss .and. EqualRealNos(tipLossConst,0.0_ReKi) ) .or. ( useHubloss .and. EqualRealNos(hubLossConst,0.0_ReKi) ) ) then
   !   ! We are simply going to bail if we are using tiploss and tipLossConst = 0 or using hubloss and hubLossConst=0, regardless of phi!
   !   fzero =  0.0_ReKi
   !   a     =  1.0_ReKi
   !   ap    = -1.0_ReKi
   !   return      
   !else if ( EqualRealNos(phi, 0.0_ReKi) ) then 
   !   fzero =  0.0_ReKi
   !   a     =  1.0_ReKi
   !   if (wakerotation) then 
   !      ap = -1.0_ReKi
   !   else
   !      ap = 0.0_ReKi
   !   end if
   !   
   !   return
   !end if
   
   !.....................................................
   ! Temporary variables:
   !.....................................................
   sphi = sin(phi)
   cphi = cos(phi)
   
   
   !.....................................................
   ! Prandtl's tip and hub loss factor:
   !.....................................................

   F = getHubTipLossCorrection(sphi, useHubLoss, useTipLoss, hubLossConst, tipLossConst) ! Prandtl's tip and hub loss factor
   
    
   !.....................................................
   ! compute axial induction factor:
   !.....................................................
   sigma_p = B*chord/(TwoPi*r)  ! local solidity   
   k = sigma_p*cn/4.0_ReKi/F/sphi/sphi
    
   
   momentumRegion = (phi > 0.0_ReKi .and. Vx >= 0.0_ReKi) .or. (phi < 0.0_ReKi .and. Vx < 0.0_ReKi)
   if (momentumRegion) then  ! momentum/empirical

 
        ! update axial induction factor
      if (k <= 2.0_ReKi/3.0_ReKi) then  ! momentum state for a < 0.4
         
         if ( EqualRealNos(k,-1.0_ReKi) ) then
            a = -sign(InductionLimit, 1.0_ReKi+k)
         else   
            a = k/(1.0_ReKi+k)
         end if
         
         if (k<-1.0_ReKi) then ! k < -1 cannot be a solution in momentum region (this is equivalent to a>1.0)
            IsValidSolution = .false.
         end if
         
         
         ! note that we'll put a max on the magnitude of 'a' later 
                  
      else  ! Glauert(Buhl) correction for 0.4<a<1

         temp = 2.0_ReKi*F*k
         g1 = temp - (10.0_ReKi/9.0_ReKi-F)
         g2 = temp - ( 4.0_ReKi/3.0_ReKi-F)*F
         g3 = temp - (25.0_ReKi/9.0_ReKi-2.0_ReKi*F)

         if (abs(g3) < 1e-6_ReKi) then  ! avoid singularity
            a = 1.0_ReKi - 0.5_ReKi/sqrt(g2)
         else
            a = (g1 - sqrt(abs(g2))) / g3 ! bjj: g2 should always be > 0, but just in case there are numerical issues, I will add the abs() here
         end if

      end if
      
   else  ! propeller brake
      
            
      if ( EqualRealNos(k,1.0_ReKi) ) then
         IsValidSolution = .false.
         a = InductionLimit
      else
         a = k/(k-1.0_ReKi)
      end if

      
      if (k<=1.0_ReKi) then ! k <= 1 cannot be a solution in propeller brake region (this is equivalent to a<1.0)
         IsValidSolution = .false.
      else if (a > MaxAxInd) then ! propeller brake region is for induction factors > 1, but not too large; 
         ! note that we use k in the residual equation instead of a in the propeller brake region, so we can put the limit here
         a = MaxAxInd
      end if         
      
   end if

   !.....................................................
   ! compute tangential induction factor:
   !.....................................................
    
   if (wakerotation) then 
   
         ! compute tangential induction factor
      if ( EqualRealNos(cphi,0.0_ReKi) ) then
         
         ap = -1.0_ReKi
         kp =  sign(InductionLimit, ct*sphi)*sign(1.0_ReKi,Vx)
         
      else
         
         kp = sigma_p*ct/4.0_ReKi/F/sphi/cphi
         if (Vx < 0.0_ReKi) then 
            kp = -kp
         end if
         
      
         if ( EqualRealNos(kp,1.0_ReKi) ) then
            ap = sign(InductionLimit, 1.0_ReKi-kp)
         else
            ap = kp/(1.0_ReKi-kp)
         end if
         
         ! bandaid so that this doesn't blow up. Note that we're not using ap in the residual calculation, so we can modify it here.
         ap = min( ap, MaxTnInd)
         ap = max( ap, MinTnInd)
         
      end if
         
            
   else 
      
      ! we're not computing tangential induction:       
      ap = 0.0_ReKi
      kp = 0.0_ReKi
      
   end if

    
   !.....................................................
   ! error function (residual)
   !.....................................................
   lambda_r = Vy/Vx

   if (momentumRegion) then  ! momentum/empirical
      if ( EqualRealNos(a, 1.0_ReKi) ) then
         fzero = - cphi/lambda_r*(1-kp)
      else       
         fzero = sphi/(1-a) - cphi/lambda_r*(1-kp)

         ! bandaid so that axial induction doesn't blow up
         a = max(a,MinAxInd)
      end if
      
   else  ! propeller brake region
      fzero = sphi*(1-k) - cphi/lambda_r*(1-kp)
   end if
    
   if (.not. IsValidSolution) then
      a = 0.0_ReKi
      ap = 0.0_ReKi
   end if
   
   
end subroutine inductionFactors
!-----------------------------------------------------------------------------------------
!> This function computes \f$F\f$, the hub/tip loss correction
real(reKi) function getHubTipLossCorrection(sphi, useHubLoss, useTipLoss, hubLossConst, tipLossConst) result(F)

   real(ReKi), intent(in) :: sphi           !< sine of local inflow angle, sin(phi)
   real(ReKi), intent(in) :: hubLossConst   !< hub loss constant [p%hubLossConst]
   real(ReKi), intent(in) :: tipLossConst   !< tip loss constant [p%tipLossConst]
   logical,    intent(in) :: useHubLoss     !< hub-loss flag [p%useHubLoss]
   logical,    intent(in) :: useTipLoss     !< tip-loss flag [p%useTipLoss]


   real(ReKi) :: factortip, Ftip, factorhub, Fhub

   !.....................................................
   ! Prandtl's tip and hub loss factor:
   !.....................................................

   Ftip = 1.0_ReKi     ! default tip loss value 
   Fhub = 1.0_ReKi     ! default hub loss value
   
   if (.not. EqualRealNos(sphi,0.0_ReKi)) then
      if ( useTipLoss ) then
         factortip = tipLossConst/abs(sphi)
         Ftip = TwoByPi*acos(min(1.0_ReKi,exp(-factortip)))
         ! else Ftip = 1.0_ReKi  ! TwoByPi*Pi/2
      end if

      if ( useHubLoss ) then
         factorhub = hubLossConst/abs(sphi)
         Fhub = TwoByPi*acos(min(1.0_ReKi,exp(-factorhub)))
         ! else Ftip = 1.0_ReKi  ! TwoByPi*Pi/2
      end if
   end if
      
   F = Ftip * Fhub
   
end function getHubTipLossCorrection
!-----------------------------------------------------------------------------------------

end module BEMTUncoupled
