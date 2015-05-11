module BEMTUnCoupled
 
   use NWTC_Library
   use BladeElement


   implicit none
   
   
   
   
   private
   
   public :: ComputeAirfoilCoefs
   public :: UncoupledErrFn
   public :: BEMTU_InductionWithResidual
   
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
        W = sqrt((Vx*(1-axInduction))**2 + (Vy)**2 ) !*(1+tanInduction))**2)
    !end if

    Re = airDens * W * chord / mu


end subroutine BEMTU_Wind

!----------------------------------------------------------------------------------------------------------------------------------  
subroutine ComputeAirfoilCoefs( phi, axInduction, tanInduction, Vx, Vy, chord, theta, airDens, mu, useAIDrag, useTIDrag, AFInfo, &
                      UA_Flag, p_UA, xd_UA, OtherState_UA, &
                      AOA, Re, Cl, Cd, Cx, Cy, errStat, errMsg )
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
   real(ReKi),             intent(  out) :: AOA, Re, Cl, Cd, Cx, Cy
   integer(IntKi),         intent(  out) :: errStat       ! Error status of the operation
   character(*),           intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None 
   
   real(ReKi)                            :: W
      ! Compute AOA, Re based on current values of axInduction, tanInduction
   call BEMTU_Wind(phi, axInduction, tanInduction, Vx, Vy, chord, theta, airDens, mu, AOA, W, Re)
      
   call  BE_CalcOutputs( AFInfo, UA_Flag, AOA*R2D, W, log(Re), p_UA, xd_UA, OtherState_UA, Cl, Cd,  errStat, errMsg)  
   !call  BE_CalcOutputs(AFInfo, AOA*R2D, log(Re), Cl, Cd, errStat, errMsg) ! AOA is in degrees in this look up table and Re is in log(Re)
   if (errStat >= AbortErrLev) then
      call SetErrStat( errStat, errMsg, errStat, errMsg, 'ComputeAirfoilCoefs' ) 
      return
   end if   
      
         ! Determine Cx, Cy from Cl, Cd and phi
   call BE_CalcCxCyCoefs(phi, useAIDrag, useTIDrag, Cl, Cd, Cx, Cy)
   
end subroutine ComputeAirfoilCoefs




                           ! This is the residual calculation for the uncoupled BEM solve
real(ReKi) function BEMTU_InductionWithResidual(phi, psi, chi0, numReIterations, airDens, mu, numBlades, rlocal, rtip, chord, theta, rHub, lambda, AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, SkewWakeMod, &
                              UA_Flag, p_UA, xd_UA, OtherState_UA, &
                              AOA, Re, Cl, Cd, Cx, Cy, axInduction, tanInduction, chi, ErrStat, ErrMsg)
      


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
   real(ReKi),             intent(in   ) :: rHub
   real(ReKi),             intent(in   ) :: lambda
   type(AFInfoType),       intent(in   ) :: AFInfo
   real(ReKi),             intent(in   ) :: Vx
   real(ReKi),             intent(in   ) :: Vy
   logical,                intent(in   ) :: useTanInd 
   logical,                intent(in   ) :: useAIDrag
   logical,                intent(in   ) :: useTIDrag
   logical,                intent(in   ) :: useHubLoss
   logical,                intent(in   ) :: useTipLoss
   integer,                intent(in   ) :: SkewWakeMod   ! Skewed wake model
   logical,                intent(in   ) :: UA_Flag
   type(UA_ParameterType),       intent(in   ) :: p_UA           ! Parameters
   type(UA_DiscreteStateType),   intent(in   ) :: xd_UA          ! Discrete states at Time
   type(UA_OtherStateType),      intent(in   ) :: OtherState_UA  ! Other/optimization states
   real(ReKi),             intent(  out) :: AOA, Re, Cl, Cd, Cx, Cy, axInduction, tanInduction, chi
   integer(IntKi),         intent(  out) :: ErrStat       ! Error status of the operation
   character(*),           intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
  
   ! Local variables
   
   
   real(ReKi)                            :: fzero, degAOA
   integer                               :: I, kk
   
   ErrStat = ErrID_None
   ErrMsg  = ""
    
      ! Set the local version of the induction factors
   axInduction  = 0.0_ReKi  ! axInductionIN
   tanInduction = 0.0_ReKi  ! tanInductionIN
   
      ! If we say Re is dependent on axInduction, tanInduction, then we would create an iteration loop around the residual calculation
    
   do I = 1,numReIterations
      
      call ComputeAirfoilCoefs( phi, axInduction, tanInduction, Vx, Vy, chord, theta, airDens, mu, useAIDrag, useTIDrag, AFInfo, &
                                UA_Flag, p_UA, xd_UA, OtherState_UA, &
                                AOA, Re, Cl, Cd, Cx, Cy, errStat, errMsg )       
      if (errStat >= AbortErrLev) then
         call SetErrStat( errStat, errMsg, errStat, errMsg, 'BEMTU_InductionWithResidual' ) 
         return
      end if
      
         ! Determine axInduction, tanInduction for the current Cl, Cd, phi
      call inductionFactors( rlocal, rtip, chord, rHub, lambda, phi, psi, chi0, Cx, Cy, numBlades, &
                              Vx, Vy, useTanInd, useHubLoss, useTipLoss,  SkewWakeMod, &
                              fzero, axInduction, tanInduction, chi, errStat, errMsg)
      if (errStat >= AbortErrLev) then
         call SetErrStat( errStat, errMsg, errStat, errMsg, 'BEMTU_InductionWithResidual' ) 
         return
      end if
      
      BEMTU_InductionWithResidual = fzero  ! the residual
      
   end do
   
end function BEMTU_InductionWithResidual

      ! This is the residual calculation for the uncoupled BEM solve

real(ReKi) function UncoupledErrFn(phi, psi, chi0, numReIterations, airDens, mu, numBlades, rlocal, rtip, chord, theta, rHub, lambda, AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, SkewWakeMod, &
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
   real(ReKi),             intent(in   ) :: rHub
   real(ReKi),             intent(in   ) :: lambda
   type(AFInfoType),       intent(in   ) :: AFInfo
   real(ReKi),             intent(in   ) :: Vx
   real(ReKi),             intent(in   ) :: Vy
   logical,                intent(in   ) :: useTanInd 
   logical,                intent(in   ) :: useAIDrag
   logical,                intent(in   ) :: useTIDrag
   logical,                intent(in   ) :: useHubLoss
   logical,                intent(in   ) :: useTipLoss
   integer,                intent(in   ) :: SkewWakeMod   ! Skewed wake model
   logical,                intent(in   ) :: UA_Flag
   type(UA_ParameterType),       intent(in   ) :: p_UA           ! Parameters
   type(UA_DiscreteStateType),   intent(in   ) :: xd_UA          ! Discrete states at Time
   type(UA_OtherStateType),      intent(in   ) :: OtherState_UA  ! Other/optimization states
   integer(IntKi),         intent(  out) :: ErrStat       ! Error status of the operation
   character(*),           intent(  out) :: ErrMsg        ! Error message if ErrStat /= ErrID_None
  
   ! Local variables
   
   
   real(ReKi)                            :: fzero, AOA, Re, Cl, Cd, Cx, Cy, axInduction, tanInduction, chi
   integer                               :: I
   
   ErrStat = ErrID_None
   ErrMsg  = ""
    
   
      
      UncoupledErrFn = BEMTU_InductionWithResidual(phi, psi, chi0, numReIterations, airDens, mu, numBlades, rlocal, rtip, chord, theta, rHub, lambda, AFInfo, &
                              Vx, Vy, useTanInd, useAIDrag, useTIDrag, useHubLoss, useTipLoss, SkewWakeMod, &
                              UA_Flag, p_UA, xd_UA, OtherState_UA, &
                              AOA, Re, Cl, Cd, Cx, Cy, axInduction, tanInduction, chi, ErrStat, ErrMsg)
      
   
   
end function UncoupledErrFn

                              
                              
                              
recursive subroutine inductionFactors(r , Rtip, chord, Rhub, lambda, phi, azimuth, chi0, cn, ct, B, &
                              Vx, Vy, wakerotation, hubLoss, tipLoss, skewWakeMod, &
                              fzero, a, ap, chi, ErrStat, ErrMsg)

    implicit none

    ! in
    real(ReKi), intent(in) :: r, chord, Rhub, Rtip, phi, cn, ct
    REAL(ReKi),       INTENT(IN   ) :: lambda        ! Tip speed ratio
    integer, intent(in) :: B
    real(ReKi), intent(in) :: Vx, Vy
    real(ReKi), intent(in) :: chi0, azimuth
    logical, intent(in) ::  hubLoss, tipLoss, wakerotation
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
    real(ReKi) :: k, kp,  F   !cn, ct,
    real(ReKi) :: g1, g2, g3
    real(ReKi) :: Fsphi, sigma_pcn
    real(ReKi) :: phitemp
    !real(ReKi) :: chi

    errStat = ErrID_None
    errMsg  = ""
   
    
    if ( EqualRealNos(phi, 0.0_ReKi) ) then
      fzero =  0.0
      a     =  0.0
      ap    =  0.0
      return
    end if
    
   
    
    
    
    sigma_p = B/2.0_ReKi/pi*chord/r
    sphi = sin(phi)
    cphi = cos(phi)

    chi = 0.0_ReKi
    !saz = sin(azimuth)
    
    if ( EqualRealNos(azimuth, 3.141593) ) then
       saz = 0.0_ReKi
    else if ( EqualRealNos(azimuth, 1.570796) ) then
       saz = 1.0_ReKi
    else if ( EqualRealNos(azimuth, 3*1.570796 )) then
       saz = -1.0_ReKi
    else     
       saz = sin(azimuth)
    end if
    
    
    ! resolve into normal and tangential forces
    !if ( .not. useCd ) then
    !    cn = cl*cphi
    !    ct = cl*sphi
    !else
    !    cn = cl*cphi + cd*sphi
    !    ct = cl*sphi - cd*cphi
    !end if

    ! Prandtl's tip and hub loss factor
    Ftip = 1.0_ReKi
    ! NOTE: check below isn't good enough: at tip Ftip should = 0.0, not 1.0
   ! if ( tipLoss .AND. (.NOT.(EqualRealNos(sphi, 0.0_DbKi))) ) then
    if ( tipLoss  ) then
        factortip = B/2.0_ReKi*(Rtip - r)/(r*abs(sphi))
        Ftip      = 2.0_ReKi/pi*acos(exp(-factortip))
    end if

    Fhub = 1.0_ReKi
    !if ( hubLoss .AND. (.NOT.(EqualRealNos(sphi, 0.0_DbKi))) ) then
    if ( hubLoss ) then
        factorhub = B/2.0_ReKi*(r - Rhub)/(Rhub*abs(sphi))
        Fhub      = 2.0_ReKi/pi*acos(exp(-factorhub))
    end if

    F = Ftip * Fhub

    ! bem parameters
    
    !Fsphi     = 4.0_ReKi*F*sphi**2 
    !sigma_pcn = sigma_p*cn
    
    k = sigma_p*cn/4.0_ReKi/F/sphi/sphi

    ! compute axial induction factor
    if (phi > 0.0_ReKi) then  ! momentum/empirical

        
        
    
        ! update axial induction factor
        if (k <= 2.0_ReKi/3.0_ReKi) then  ! momentum state
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

       ! Skewed wake correction
         
   if ( skewWakeMod == 1 ) then
      chi = (0.6_ReKi*a + 1.0_ReKi)*chi0
      !yawCorr = max(0.0,chi0-0.5236)
      !yawCorr = min(0.785,yawCorr)
      yawCorr = (15.0_ReKi*pi/64.0_ReKi*tan(chi/2.0_ReKi) * (r/Rtip) * saz)
      a = a * (1.0 +  yawCorr) ! *(-yawCorr/0.785 + 1) )
        !a = min(a, 0.9)
   else if ( skewWakeMod == 2 ) then
      chi = (0.6_ReKi*a + 1.0_ReKi)*chi0
      a = a * ( 1.0_ReKi + 0.63_ReKi*(r/Rtip)**2 )*(c1 + c2*lambda + c3*lambda**2 + c4*lambda**3 )*( chi*180.0_ReKi/pi )*saz   ! chi needs to be in degrees here.
   end if
   
    ! compute tangential induction factor
    kp = sigma_p*ct/4.0_ReKi/F/sphi/cphi
    ap = kp/(1.0_ReKi-kp)

    if (.not. wakerotation) then
        ap = 0.0_ReKi
        kp = 0.0_ReKi
    end if

    !if ( skewWakeMod > 0 ) then  
    !  phitemp = InflowAngle(Vx_in, Vy_in, REAL(a, ReKi), REAL(ap))
    !  call inductionFactors(r_in     , Rtip_in, chord_in, Rhub_in,  lambda_in, phitemp, azimuth_in, yaw_in  , cn_in, ct_in, B, &
    !                          Vx_in, Vy_in, wakerotation,   hubLoss , tipLoss   , 0, &
    !                          fzero_out, a_out,           ap_out,           chi_out, ErrStat, ErrMsg)
    !  return
    !
    !end if  
    
    ! error function
    lambda_r = Vy/Vx
    if (phi > 0) then  ! momentum/empirical
        fzero = sphi/(1-a) - cphi/lambda_r*(1-kp)
    else  ! propeller brake region
        fzero = sphi*(1-k) - cphi/lambda_r*(1-kp)
    end if
    
end subroutine inductionFactors
    
recursive subroutine inductionFactors2(r_in     , Rtip_in, chord_in, Rhub_in,                      lambda_in, phi_in, azimuth_in, yaw_in  , cn_in, ct_in, B, &
                              Vx_in, Vy_in, wakerotation,   hubLoss , tipLoss   , yawcorrection, &
                              fzero_out, a_out,           ap_out,           chi_out, ErrStat, ErrMsg)

    implicit none

    integer, parameter :: dp = SELECTED_REAL_KIND(  6,  30 )  !kind(0.d0) !!

    ! in
    real(ReKi), intent(in) :: r_in, chord_in, Rhub_in, Rtip_in, phi_in, cn_in, ct_in
    REAL(ReKi),       INTENT(IN   ) :: lambda_in        ! Tip speed ratio
    integer, intent(in) :: B
    real(ReKi), intent(in) :: Vx_in, Vy_in
    real(ReKi), intent(in) :: yaw_in, azimuth_in
    logical, intent(in) ::  hubLoss, tipLoss, wakerotation
    integer, intent(in) :: yawcorrection  ! useCd,
    
    
    
    

    ! out
    real(ReKi), intent(out) :: fzero_out, a_out, ap_out
    REAL(ReKi),       INTENT(  OUT) :: chi_out
    INTEGER(IntKi),   INTENT(  OUT) :: ErrStat       ! Error status of the operation
    CHARACTER(*),     INTENT(  OUT) :: ErrMsg        ! Error message if ErrStat /= ErrID_None

    ! local
    
       ! constants
   REAL(dp), PARAMETER :: c1 = 2.6296e-3 
   REAL(dp), PARAMETER :: c2 = 1.6222e-3 
   REAL(dp), PARAMETER :: c3 = -1.1111e-5 
   REAL(dp), PARAMETER :: c4 = -3.70371e-6 
   
    real(dp)  :: r, chord, Rhub, Rtip, phi, cn, ct, yawCorr
    REAL(dp)  :: lambda        ! Tip speed ratio
    
    real(dp)  :: Vx, Vy
    real(dp)  :: yaw, azimuth
    
    real(dp) ::  sigma_p, sphi, cphi, lambda_r, saz !, pi
    real(dp) :: factortip, Ftip, factorhub, Fhub
    real(dp) :: k, kp,  F   !cn, ct,
    real(dp) :: g1, g2, g3
    real(dp) :: fzero, a, ap, Fsphi, sigma_pcn
    REAL(dp) :: chi
    real(ReKi) :: phitemp
    !real(dp) :: chi

    errStat = ErrID_None
    errMsg  = ""
   
    phi = phi_in
    if ( EqualRealNos(phi, 0.0_ReKi) ) then
      fzero_out =  0.0
      a_out     =  0.0
      ap_out    =  0.0
      return
    end if
    
    r = r_in
    chord = chord_in
    Rhub = Rhub_in
    Rtip = Rtip_in
    lambda = lambda_in
    
    cn = cn_in
    ct = ct_in
    Vx = Vx_in
    Vy = Vy_in
    yaw = yaw_in
    azimuth = azimuth_in
    ! constants
    !pi = 3.1415926535897932_dp
    sigma_p = B/2.0_dp/pi*chord/r
    sphi = sin(phi_in)
    cphi = cos(phi_in)

    chi = 0.0_dp
    saz = sin(azimuth_in)
    
    if ( EqualRealNos(azimuth_in, 3.141593) ) then
       saz = 0.0_dp
    else if ( EqualRealNos(azimuth_in, 1.570796) ) then
       saz = 1.0_dp
    else if ( EqualRealNos(azimuth_in, 3*1.570796 )) then
       saz = -1.0_dp
    else     
       saz = sin(azimuth)
    end if
    
    
    ! resolve into normal and tangential forces
    !if ( .not. useCd ) then
    !    cn = cl*cphi
    !    ct = cl*sphi
    !else
    !    cn = cl*cphi + cd*sphi
    !    ct = cl*sphi - cd*cphi
    !end if

    ! Prandtl's tip and hub loss factor
    Ftip = 1.0_dp
    ! NOTE: check below isn't good enough: at tip Ftip should = 0.0, not 1.0
   ! if ( tipLoss .AND. (.NOT.(EqualRealNos(sphi, 0.0_DbKi))) ) then
    if ( tipLoss  ) then
        factortip = B/2.0_dp*(Rtip - r)/(r*abs(sphi))
        Ftip      = 2.0_dp/pi*acos(exp(-factortip))
    end if

    Fhub = 1.0_dp
    !if ( hubLoss .AND. (.NOT.(EqualRealNos(sphi, 0.0_DbKi))) ) then
    if ( hubLoss ) then
        factorhub = B/2.0_dp*(r - Rhub)/(Rhub*abs(sphi))
        Fhub      = 2.0_dp/pi*acos(exp(-factorhub))
    end if

    F = Ftip * Fhub

    ! bem parameters
    
    !Fsphi     = 4.0_dp*F*sphi**2 
    !sigma_pcn = sigma_p*cn
    
    k = sigma_p*cn/4.0_dp/F/sphi/sphi

    ! compute axial induction factor
    if (phi > 0.0) then  ! momentum/empirical

        
        
    
        ! update axial induction factor
        if (k <= 2.0_dp/3.0) then  ! momentum state
            a = k/(1+k)

        else  ! Glauert(Buhl) correction

            g1 = 2.0_dp*F*k - (10.0_dp/9-F)
            g2 = 2.0_dp*F*k - (4.0_dp/3-F)*F
            g3 = 2.0_dp*F*k - (25.0_dp/9-2*F)

            if (abs(g3) < 1e-6_dp) then  ! avoid singularity
                a = 1.0_dp - 1.0_dp/2.0/sqrt(g2)
            else
                a = (g1 - sqrt(g2)) / g3
            end if

        end if

    else  ! propeller brake region (a and ap not directly used but update anyway)

        if (k > 1.0) then
        !if (sigma_pcn > Fsphi) then
            a =   k/(k-1) !sigma_pcn / (sigma_pcn - Fsphi )  !
        else
            a = 0.0_dp  ! dummy value
        end if

    end if

    ! apply yaw correction
    !if (yawcorrection) then
    !    
    !    chi = (0.6*a + 1.0)*yaw
    !    a = a * (1.0 + 15.0*pi/32*tan(chi/2.0) * r/Rtip * saz)
    !    a = min(a, 0.999999)
    !end if

       ! Skewed wake correction
         
   if ( yawcorrection == 1 ) then
      chi = (0.6*a + 1.0)*yaw
      !yawCorr = max(0.0,yaw-0.5236)
      !yawCorr = min(0.785,yawCorr)
      yawCorr = (15.0*pi/64*tan(chi/2.0) * (r/Rtip) * saz)
      a = a * (1.0 +  yawCorr) ! *(-yawCorr/0.785 + 1) )
        !a = min(a, 0.9)
   else if ( yawcorrection == 2 ) then
      chi = (0.6*a + 1.0)*yaw
      a = a * ( 1 + 0.63*(r/Rtip)**2 )*(c1 + c2*lambda + c3*lambda**2 + c4*lambda**3 )*( chi*180.0/pi )*saz   ! chi needs to be in degrees here.
   end if
   
    ! compute tangential induction factor
    kp = sigma_p*ct/4.0_dp/F/sphi/cphi
    ap = kp/(1-kp)

    if (.not. wakerotation) then
        ap = 0.0_dp
        kp = 0.0_dp
    end if

    !if ( yawcorrection > 0 ) then  
    !  phitemp = InflowAngle(Vx_in, Vy_in, REAL(a, ReKi), REAL(ap))
    !  call inductionFactors(r_in     , Rtip_in, chord_in, Rhub_in,  lambda_in, phitemp, azimuth_in, yaw_in  , cn_in, ct_in, B, &
    !                          Vx_in, Vy_in, wakerotation,   hubLoss , tipLoss   , 0, &
    !                          fzero_out, a_out,           ap_out,           chi_out, ErrStat, ErrMsg)
    !  return
    !
    !end if  
    
    ! error function
    lambda_r = Vy/Vx
    if (phi > 0) then  ! momentum/empirical
        fzero = sphi/(1-a) - cphi/lambda_r*(1-kp)
    else  ! propeller brake region
        fzero = sphi*(1-k) - cphi/lambda_r*(1-kp)
    end if

    ap_out = ap
    a_out  = a
    fzero_out = fzero
    chi_out   = chi
end subroutine inductionFactors2
                              
subroutine BEMT_InductionFactors( rlocal, rtip, chord, tipLossConst, hubLossConst, lambda, phi, psi, chi0, Cx, Cy, numBlades, &
                              Vx, Vy, useTanInd, useHubLoss, useTipLoss, SkewWakeMod, &
                              fzero, axInduction, tanInduction, chi, ErrStat, ErrMsg)
    
      ! in
   REAL(ReKi),       INTENT(IN   ) :: rlocal        ! straight-line distance from center-of-rotation to blade element node (m)
   REAL(ReKi),       INTENT(IN   ) :: rtip          ! straight-line distance from center-of-rotation to blade tip (m)
   REAL(ReKi),       INTENT(IN   ) :: chord         ! chord length for the blade element cross-section (m)
   REAL(ReKi),       INTENT(IN   ) :: tipLossConst  !
   REAL(ReKi),       INTENT(IN   ) :: hubLossConst  !
   REAL(ReKi),       INTENT(IN   ) :: lambda        ! Tip speed ratio
   REAL(ReKi),       INTENT(IN   ) :: phi           ! angle between local chord line and local wind vector (rad)
   REAL(ReKi),       INTENT(IN   ) :: psi           ! asimuth angle   (rad) 
   real(ReKi),       intent(in   ) :: chi0         ! Yaw angle
   REAL(ReKi),       INTENT(IN   ) :: Cx            ! coefficient of  (-)
   REAL(ReKi),       INTENT(IN   ) :: Cy            ! coefficient of  (-)
   INTEGER,          INTENT(IN   ) :: numBlades     ! number of blades (-)
   REAL(ReKi),       INTENT(IN   ) :: Vx            ! local x wind velocity (m/s)
   REAL(ReKi),       INTENT(IN   ) :: Vy            ! local y wind velocity (m/s)
   LOGICAL,          INTENT(IN   ) :: useTanInd     ! compute and use tangential induction factor
   LOGICAL,          INTENT(IN   ) :: useHubLoss    ! include hub loss factor
   LOGICAL,          INTENT(IN   ) :: useTipLoss    ! include tip loss factor
   INTEGER,          INTENT(IN   ) :: SkewWakeMod   ! Skewed wake model
    
   REAL(ReKi),       INTENT(  OUT) :: fzero
   REAL(ReKi),       INTENT(  OUT) :: axInduction
   REAL(ReKi),       INTENT(  OUT) :: tanInduction
   REAL(ReKi),       INTENT(  OUT) :: chi
   INTEGER(IntKi),   INTENT(  OUT) :: ErrStat       ! Error status of the operation
   CHARACTER(*),     INTENT(  OUT) :: ErrMsg        ! Error message if ErrStat /= ErrID_None



      ! local
   REAL(DbKi) :: sigma_p, sphi, cphi, lambda_r
   REAL(DbKi) :: factortip, Ftip, factorhub, Fhub
   REAL(DbKi) :: k, kp, F
   REAL(DbKi) :: g1, g2sqrt, g3
    
      ! constants
   REAL(ReKi), PARAMETER :: c1 = 2.6296e-3 
   REAL(ReKi), PARAMETER :: c2 = 1.6222e-3 
   REAL(ReKi), PARAMETER :: c3 = -1.1111e-5 
   REAL(ReKi), PARAMETER :: c4 = -3.70371e-6 

      ! Initialize ErrStat

   errStat = ErrID_None
   errMsg  = ""
    
   axInduction  = 0.0_ReKi
   tanInduction = 0.0_ReKi
    
   sigma_p = (numBlades*chord)/(2.0*pi*rlocal)                   ! solidity factor  ! checked GJH
   sphi = sin(phi)
   cphi = cos(phi)

if ( .NOT. (EqualRealNos(phi,0.0_ReKi) ) ) then
      
                              
      ! Prandtl's tip and hub loss factor
   Ftip = 1.0
   if ( useTipLoss ) then
      factortip = tipLossConst/abs(sphi)
      Ftip = (2.0/pi)*acos(exp(-factortip))
   end if

   Fhub = 1.0
   if ( useHubLoss ) then
      factorhub = hubLossConst/abs(sphi)
      Fhub = (2.0/pi)*acos(exp(-factorhub))
   end if

   F = Ftip * Fhub

      ! bem parameters
   k  = sigma_p*Cx/4.0_DbKi/F/sphi/sphi  ! checked GJH
    

      ! compute axial induction factor
   if (phi > 0) then  ! momentum/empirical

      ! update axial induction factor
      if (k <= 2.0/3.0) then  ! momentum region
            
         axInduction = k/(1.0_DbKi+k)
         if ( EqualRealNos(axInduction, 1.0_ReKi) ) then
            axInduction = 0.0 !1.0 - 1e-6  ! TODO:  This is a test to avoid a singularity GJH
         end if
            
      else  ! Glauert(Buhl) correction
            
         g2sqrt = sqrt(2*F*k - (4.0_DbKi /3.0 -   F)*F)
         g3     =      2*F*k - (25.0_DbKi/9.0 - 2*F)

         if ( g3 < 1e-6)  then  ! TODO:  Make this really our 'epsilon'
               axInduction = 1.0_DbKi - 1.0_DbKi/(2.0*g2sqrt)
         else
          
         !IF ( abs(k - (25.0/18/F - 1.0)) < 1e-6) THEN
         !    k = k + 1.0e-5  ! avoid singularity
         !END IF

               g1 = 2*F*k - (10.0_DbKi/9.0 -   F)
               axInduction  = (g1 - g2sqrt) / g3
         end if 
      end if

   else  ! propeller brake region (axInduction and tanInduction not directly used but update anyway)

      if (k > 1.0) then
         axInduction = k/(k-1)
         if ( EqualRealNos(axInduction, 1.0_ReKi) ) then
            axInduction = 0.0_ReKi
         end if
         
      else
         axInduction = 0.0  ! dummy value
      end if

   end if

      
      ! Estimated wake skew angle from Burton
   chi = (0.6*axInduction + 1) * chi0 
   
   
      ! Skewed wake correction
       
   if ( SkewWakeMod == 1 ) then
      axInduction = axInduction * ( 1 + (15.0*pi/64.0)*tan(chi/2.0)*(rlocal/rtip)*sin(psi) )  !TODO: Verify this equation.  It doesn't match previous AD versions 
  ! else if ( SkewWakeMod == 2 ) then
  !    axInduction = axInduction * ( 1 + 0.63*(rlocal/rtip)**2 )*(c1 + c2*lambda + c3*lambda**2 + c4*lambda**3 )*( chi*180.0/pi )*sin(psi)   ! chi needs to be in degrees here.
   end if
        
      ! compute tangential induction factor
   if (useTanInd) then
      kp = sigma_p*Cy/(4*F*sphi*cphi)
      tanInduction = kp/(1-kp)
   else
      tanInduction = 0.0_ReKi
      kp = 0.0_ReKi
   end if
   
   
   
end if ! if phi /= 0.0
                              
   ! error function
!if ( EqualRealNos(Vx, 0.0_ReKi) ) then
!   fzero = sphi/(1-axInduction) 
!else
   fzero = sphi/(1-axInduction) - (cphi*Vx)/(Vy*(1+tanInduction))
   !lambda_r = Vy/Vx
   !if (phi > 0) then  ! momentum/empirical
   !   fzero = sphi/(1-axInduction) - cphi/lambda_r*(1-kp)
   !else  ! propeller brake region
   !   fzero = sphi*(1-k) - cphi/lambda_r*(1-kp)
   !end if 
   
!end if

      ! error function
    
   !
    

END SUBROUTINE BEMT_InductionFactors
                              
                              
end module BEMTUncoupled