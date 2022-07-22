!> Biot-Savart law functions
!! NOTE: these functions should be independent of the framework types
module FVW_BiotSavart 

   use NWTC_Library, only: ReKi, IntKi
   use OMP_LIB

   implicit none

   real(ReKi),parameter :: PRECISION_UI  = epsilon(1.0_ReKi)/100 !< NOTE assuming problem of size 1
   real(ReKi),parameter :: MIN_EXP_VALUE=-10.0_ReKi
   real(ReKi),parameter :: MINDENOM=0.0_ReKi
!    real(ReKi),parameter :: MINDENOM=1e-15_ReKi
   real(ReKi),parameter :: MINNORM=1e-4

   integer(IntKi), parameter :: idRegNone       = 0
   integer(IntKi), parameter :: idRegRankine    = 1
   integer(IntKi), parameter :: idRegLambOseen  = 2
   integer(IntKi), parameter :: idRegVatistas   = 3
   integer(IntKi), parameter :: idRegOffset     = 4
   integer(IntKi), parameter :: idRegExp        = 1
   integer(IntKi), parameter :: idRegCompact    = 2
   integer(IntKi), parameter, dimension(5) :: idRegVALID      = (/idRegNone,idRegRankine,idRegLambOseen,idRegVatistas,idRegOffset/)
   integer(IntKi), parameter, dimension(3) :: idRegPartVALID  = (/idRegNone,idRegExp,idRegCompact/)

   real(ReKi),parameter    :: fourpi_inv =  0.25_ReKi / ACOS(-1.0_Reki )

contains


!> Induced velocity from one segment at one control points
subroutine ui_seg_11(DeltaPa, DeltaPb, SegGamma, RegFunction, RegParam1, Uind)
   ! Input/output arguments 
   real(ReKi), dimension(3), intent(in) :: DeltaPa    !< 3 x 1   Pcp-P1  [m]
   real(ReKi), dimension(3), intent(in) :: DeltaPb    !< 3 x 1   Pcp-P2  [m]
   real(ReKi), intent(in)               :: SegGamma   !< Circulation  [m^2/s]
   integer, intent(in)                  :: RegFunction!< Regularization model
   real(ReKi), intent(in)               :: RegParam1  !< Regularization parameter (core radius) [m]
   real(ReKi), dimension(3),intent(out) :: Uind       !< Induced velocity (no side effects) [m/s]
   ! Variables declaration 
   real(ReKi),dimension(3) :: crossprod      !< 
   real(ReKi)              :: denominator    !< 
   real(ReKi)              :: r_bar2         !< (r/rc)^2
   real(ReKi)              :: Kv             !< 
   real(ReKi)              :: norm_a, norm_b !< 
   real(ReKi)              :: norm2_r0       !< Squared length of the segment
   real(ReKi)              :: norm2_orth     !< Squared distance orthogonal to the segment
   real(ReKi)              :: xa, ya, za, xb, yb, zb !< Coordinates of X-Xa and X-Xb
   real(ReKi)              :: exp_value      !< 
   ! 
   Uind(1:3)=0.0_ReKi
   xa=DeltaPa(1); ya=DeltaPa(2); za=DeltaPa(3)
   xb=DeltaPb(1); yb=DeltaPb(2); zb=DeltaPb(3)
   norm_a      = sqrt(xa*xa + ya*ya + za*za)
   norm_b      = sqrt(xb*xb + yb*yb + zb*zb)
   denominator  = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb) ! |r1|*|r2|*(|r1|*|r2| + r1.r2)
   if (denominator>PRECISION_UI) then
      crossprod(1) = ya*zb-za*yb; crossprod(2) = za*xb-xa*zb; crossprod(3) = xa*yb-ya*xb
      norm2_orth   = crossprod(1)**2 + crossprod(2)**2 + crossprod(3)**2
      if (norm2_orth>PRECISION_UI) then ! On the singularity, Uind(1:3)=0.0_ReKi
         norm2_r0     = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb) 
         if (norm2_r0>PRECISION_UI) then ! segment of zero length
            ! --- Far field TODO
            ! --- Regularization (close field)
            norm2_orth   = norm2_orth/norm2_r0 ! d = (r1xr2)/r0
            select case (RegFunction) !
            case ( idRegNone )         ! No vortex core model
               Kv=1.0_ReKi
            case ( idRegRankine )      ! Rankine
               r_bar2    = norm2_orth/ RegParam1**2
               if (r_bar2<1) then 
                  Kv=r_bar2
               else
                  Kv=1.0_ReKi 
               end if 
            case ( idRegLambOseen )    ! Lamb-Oseen
               r_bar2    = norm2_orth/ RegParam1**2
               exp_value = -1.25643_ReKi*r_bar2
               if(exp_value<MIN_EXP_VALUE) then ! Remove me when Far distance implemented
                  Kv = 1.0_ReKi
               else
                  Kv = 1.0_ReKi-exp(exp_value)
               endif
            case ( idRegVatistas )    ! Vatistas n=2
               r_bar2 = norm2_orth/ RegParam1**2
               Kv     = r_bar2/sqrt(1.0_ReKi+r_bar2**2)
            case ( idRegOffset )      ! Cut-off radius 
               Kv        = 1.0_ReKi
               denominator=denominator+RegParam1**2*norm2_r0
            case default
               Kv=1.0_ReKi !< Should be an error
            end select 
            Kv=SegGamma*fourpi_inv*Kv*(norm_a+norm_b)/denominator
            Uind(1:3) = Kv*crossprod(1:3)
         endif
      endif
   endif
end subroutine ui_seg_11


!> Induced velocity from a list of segments defined by Connectivity (SegConnct) and Points (SegPoints)
!! NOTE: this function has side effects and expects Uind_out to be initialized!
!! The function can compute the velocity on part of the segments and part of the control points.
!! This feature is useful if some parallelization is used, while common storage vectors are used.
subroutine ui_seg(iCPStart, iCPEnd, CPs, &
      iSegStart, iSegEnd, nSegTot, nSegPTot, SegPoints, SegConnct, SegGamma,  &
      RegFunction, RegParam, Uind_out)
   real(ReKi), dimension(:,:),     intent(in)    :: CPs         !< Control points (3 x nCPs++)
   integer(IntKi),                 intent(in)    :: iCPStart    !< Index where we start in Control points array
   integer(IntKi),                 intent(in)    :: iCPEnd      !< Index where we end in Control points array
   real(ReKi), dimension(:,:),     intent(in)    :: SegPoints   !< Segment points (3 x nSegPTot)
   integer(IntKi), dimension(:,:), intent(in)    :: SegConnct   !< Connectivity, indices of segments points iSeg1, iSeg2, iDepth, iSpan
   real(ReKi), dimension(:),       intent(in)    :: SegGamma    !< Segment circulation (nSegTot)
   integer(IntKi),                 intent(in)    :: iSegStart   !< Index in SegConnct, and SegGamma where we start
   integer(IntKi),                 intent(in)    :: iSegEnd     !< Index in SegConnct, and SegGamma where we end
   integer(IntKi),                 intent(in)    :: nSegTot     !< Total number of segments
   integer(IntKi),                 intent(in)    :: nSegPTot    !< Total number of segment points
   integer(IntKi),                 intent(in)    :: RegFunction !< Regularization model
   real(ReKi), dimension(:),       intent(in)    :: RegParam    !< Regularization parameter (nSegTot)
   real(ReKi), dimension(:,:)    , intent(inout) :: Uind_out    !< Induced velocity vector - Side effects!!! (3 x nCPs++)
   ! Variables
   integer(IntKi) :: icp, is
   real(ReKi), dimension(3) :: Uind           !< 
   real(ReKi), dimension(3) :: P1, P2         !< Extremities of a given segment
   ! Variables declaration 
   real(ReKi),dimension(3) :: crossprod       !< 
   real(ReKi)              :: denominator     !< 
   real(ReKi)              :: r_bar2          !< (r/rc)^2
   real(ReKi)              :: Kv              !< 
   real(ReKi)              :: norm_a, norm_b  !< 
   real(ReKi)              :: norm2_orth      !< 
   real(ReKi)              :: norm2_r0        !< Squared length of the segment d = (r1xr2)/r0
   real(ReKi)              :: xa, ya, za, xb, yb, zb !< Coordinates of X-Xa and X-Xb
   real(ReKi)              :: exp_value       !< 

   ! Branching based on regularization model
   ! NOTE: copy paste of code is done for optimization!
   !       The only thing changing is the part labelled "regularization"
   select case (RegFunction) 
   case ( idRegNone ) ! No vortex core 
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,Uind,P1,P2,crossprod,denominator,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         do is=iSegStart,iSegEnd ! loop on selected segments 
            Uind = 0.0_ReKi
            P1   = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2   = SegPoints(1:3, SegConnct(2,is))
            xa=CPs(1,icp)-P1(1); ya=CPs(2,icp)-P1(2); za=CPs(3,icp)-P1(3);
            xb=CPs(1,icp)-P2(1); yb=CPs(2,icp)-P2(2); zb=CPs(3,icp)-P2(3);
            norm_a      = sqrt(xa*xa + ya*ya + za*za)
            norm_b      = sqrt(xb*xb + yb*yb + zb*zb)
            denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
            ! --- Far field TODO
            if (denominator>PRECISION_UI) then
               crossprod(1) = ya*zb-za*yb; crossprod(2) = za*xb-xa*zb; crossprod(3) = xa*yb-ya*xb
               norm2_orth   = crossprod(1)**2 + crossprod(2)**2 + crossprod(3)**2
               if (norm2_orth>PRECISION_UI) then ! On the singularity, Uind(1:3)=0.0_ReKi
                  norm2_r0     = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb) 
                  if (norm2_r0>PRECISION_UI) then
                     ! --- Far field TODO
                     ! --- NO Regularization (close field)
                     Kv        = SegGamma(is)*fourpi_inv*(norm_a+norm_b)/(denominator + MINDENOM)
                     Uind(1:3) = Kv*crossprod(1:3)
                  end if
               end if
            end if
            Uind_out(1:3,icp) = Uind_out(1:3,icp)+Uind(1:3)
         end do ! Loop on segments
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL
      
   case ( idRegRankine )      ! Rankine
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,Uind,P1,P2,crossprod,denominator,r_bar2,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         do is=iSegStart,iSegEnd ! loop on selected segments 
            Uind = 0.0_ReKi
            P1   = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2   = SegPoints(1:3, SegConnct(2,is))
            xa=CPs(1,icp)-P1(1); ya=CPs(2,icp)-P1(2); za=CPs(3,icp)-P1(3);
            xb=CPs(1,icp)-P2(1); yb=CPs(2,icp)-P2(2); zb=CPs(3,icp)-P2(3);
            norm_a      = sqrt(xa*xa + ya*ya + za*za)
            norm_b      = sqrt(xb*xb + yb*yb + zb*zb)
            denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
            if (denominator>PRECISION_UI) then
               crossprod(1) = ya*zb-za*yb; crossprod(2) = za*xb-xa*zb; crossprod(3) = xa*yb-ya*xb
               norm2_orth   = crossprod(1)**2 + crossprod(2)**2 + crossprod(3)**2
               if (norm2_orth>PRECISION_UI) then ! On the singularity, Uind(1:3)=0.0_ReKi
                  norm2_r0     = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb) 
                  if (norm2_r0>PRECISION_UI) then
                     ! --- Far field TODO
                     ! --- Regularization (close field) --- Rankine
                     norm2_orth = norm2_orth/norm2_r0 ! d = (r1xr2)/r0
                     r_bar2     = norm2_orth/ RegParam(is)**2
                     if (r_bar2<1) then 
                        Kv=r_bar2
                     else
                        Kv=1.0_ReKi 
                     end if 
                     Kv        = SegGamma(is)*fourpi_inv*Kv*(norm_a+norm_b)/(denominator + MINDENOM)
                     Uind(1:3) = Kv*crossprod(1:3)
                  end if
               end if ! denominator size or distances too small
            end if ! 
            Uind_out(1:3,icp) = Uind_out(1:3,icp)+Uind(1:3)
         end do ! Loop on segments
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL

   case ( idRegLambOseen )      ! LambOseen
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,Uind,P1,P2,crossprod,denominator,r_bar2,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb,exp_value) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         do is=iSegStart,iSegEnd ! loop on selected segments 
            Uind = 0.0_ReKi
            P1   = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2   = SegPoints(1:3, SegConnct(2,is))
            xa=CPs(1,icp)-P1(1); ya=CPs(2,icp)-P1(2); za=CPs(3,icp)-P1(3);
            xb=CPs(1,icp)-P2(1); yb=CPs(2,icp)-P2(2); zb=CPs(3,icp)-P2(3);
            norm_a      = sqrt(xa*xa + ya*ya + za*za)
            norm_b      = sqrt(xb*xb + yb*yb + zb*zb)
            denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
            if (denominator>PRECISION_UI) then
               crossprod(1) = ya*zb-za*yb; crossprod(2) = za*xb-xa*zb; crossprod(3) = xa*yb-ya*xb
               norm2_orth   = crossprod(1)**2 + crossprod(2)**2 + crossprod(3)**2
               if (norm2_orth>PRECISION_UI) then ! On the singularity, Uind(1:3)=0.0_ReKi
                  norm2_r0     = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb) 
                  if (norm2_r0>PRECISION_UI) then
                     ! --- Far field TODO
                     ! --- Regularization (close field) --- Lamb Oseen
                     norm2_orth = norm2_orth/norm2_r0 ! d = (r1xr2)/r0
                     r_bar2     = norm2_orth/ RegParam(is)**2
                     exp_value  = -1.25643_ReKi*r_bar2
                     if(exp_value<MIN_EXP_VALUE) then ! Remove me when Far distance implemented
                        Kv = 1.0_ReKi
                     else
                        Kv = 1.0_ReKi-exp(exp_value)
                     endif
                     Kv        = SegGamma(is)*fourpi_inv*Kv*(norm_a+norm_b)/(denominator + MINDENOM)
                     Uind(1:3) = Kv*crossprod(1:3)
                  endif 
               end if ! denominator size or distances too small
            end if ! 
            Uind_out(1:3,icp) = Uind_out(1:3,icp)+Uind(1:3)
         end do ! Loop on segments
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL

   case ( idRegVatistas )      ! Vatistas n=2
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,Uind,P1,P2,crossprod,denominator,r_bar2,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         do is=iSegStart,iSegEnd ! loop on selected segments 
            Uind = 0.0_ReKi
            P1   = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2   = SegPoints(1:3, SegConnct(2,is))
            xa=CPs(1,icp)-P1(1); ya=CPs(2,icp)-P1(2); za=CPs(3,icp)-P1(3);
            xb=CPs(1,icp)-P2(1); yb=CPs(2,icp)-P2(2); zb=CPs(3,icp)-P2(3);
            norm_a      = sqrt(xa*xa + ya*ya + za*za)
            norm_b      = sqrt(xb*xb + yb*yb + zb*zb)
            denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
            if (denominator>PRECISION_UI) then
               crossprod(1) = ya*zb-za*yb; crossprod(2) = za*xb-xa*zb; crossprod(3) = xa*yb-ya*xb
               norm2_orth   = crossprod(1)**2 + crossprod(2)**2 + crossprod(3)**2
               if (norm2_orth>PRECISION_UI) then ! On the singularity, Uind(1:3)=0.0_ReKi
                  norm2_r0     = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb) 
                  if (norm2_r0>PRECISION_UI) then
                     ! --- Far field TODO
                     ! --- Regularization (close field) --- Vatistas
                     norm2_orth = norm2_orth/norm2_r0 ! d = (r1xr2)/r0
                     r_bar2     = norm2_orth/RegParam(is)**2
                     Kv         = r_bar2/sqrt(1+r_bar2**2)
                     Kv         = SegGamma(is)*fourpi_inv*Kv*(norm_a+norm_b)/(denominator + MINDENOM)
                     Uind(1:3)  = Kv*crossprod(1:3)
                  end if
               end if ! denominator size or distances too small
            end if ! 
            Uind_out(1:3,icp) = Uind_out(1:3,icp)+Uind(1:3)
         end do ! Loop on segments
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL

   case ( idRegOffset )      ! Denominator offset
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,Uind,P1,P2,crossprod,denominator,r_bar2,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         do is=iSegStart,iSegEnd ! loop on selected segments 
            Uind      = 0.0_ReKi
            P1        = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2        = SegPoints(1:3, SegConnct(2,is))
            xa=CPs(1,icp)-P1(1); ya=CPs(2,icp)-P1(2); za=CPs(3,icp)-P1(3);
            xb=CPs(1,icp)-P2(1); yb=CPs(2,icp)-P2(2); zb=CPs(3,icp)-P2(3);
            norm_a      = sqrt(xa*xa + ya*ya + za*za)
            norm_b      = sqrt(xb*xb + yb*yb + zb*zb)
            denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
            if (denominator>PRECISION_UI) then
               crossprod(1) = ya*zb-za*yb; crossprod(2) = za*xb-xa*zb; crossprod(3) = xa*yb-ya*xb
               norm2_orth   = crossprod(1)**2 + crossprod(2)**2 + crossprod(3)**2
               if (norm2_orth>PRECISION_UI) then ! On the singularity, Uind(1:3)=0.0_ReKi
                  norm2_r0     = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb) 
                  if (norm2_r0>PRECISION_UI) then
                     ! --- Far field TODO
                     ! --- Regularization (close field) -- Offset
                     denominator = denominator+RegParam(is)**2*norm2_r0
                     Kv          = SegGamma(is)*fourpi_inv*(norm_a+norm_b)/(denominator + MINDENOM)
                     Uind(1:3)   = Kv*crossprod(1:3)
                  end if
               end if ! denominator size or distances too small
            end if ! 
            Uind_out(1:3,icp) = Uind_out(1:3,icp)+Uind(1:3)
         end do ! Loop on segments
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL
   case default
      print*,'[ERROR] Unknown RegFunction for segment',RegFunction
      STOP
   end select
end subroutine ui_seg

!> Induced velocity from `nPart` particles at `nCPs` control points. The velocity gradient is not computed
subroutine ui_part_nograd(CPs, Part, Alpha, RegFunction, RegParam, UIout, nCPs, nPart)
   integer(IntKi),               intent(in)    :: nCPs
   integer(IntKi),               intent(in)    :: nPart
   real(ReKi), dimension(:,:),   intent(in)    :: CPs         !< Control points (3 x nCPs)
   real(ReKi), dimension(:,:),   intent(inout) :: UIout       !< Induced velocity, with side effects! (3 x nCPs)
   real(ReKi), dimension(:,:),   intent(in)    :: Part        !< Particle positions (3 x nPart)
   real(ReKi), dimension(:,:),   intent(in)    :: Alpha       !< Particle intensity [m^3/s] (3 x nPart) omega dV= alpha
   integer(IntKi),               intent(in)    :: RegFunction !< Regularization function 
   real(ReKi), dimension(:),     intent(in)    :: RegParam    !< Regularization parameter (nPart)
   real(ReKi), dimension(3) :: UItmp   !< 
   real(ReKi), dimension(3) :: DP      !< 
   integer :: icp,ip
   ! TODO: inlining of regularization
   !$OMP PARALLEL DEFAULT(SHARED)
   !$OMP DO PRIVATE(icp,ip, DP, UItmp) schedule(runtime)
   do icp=1,nCPs ! loop on CPs 
      do ip=1,nPart ! loop on particles
         UItmp(1:3) = 0.0_ReKi
         DP(1:3)    = CPs(1:3,icp)-Part(1:3,ip)
         call ui_part_nograd_11(DP, Alpha(1:3,ip), RegFunction , RegParam(ip), UItmp)
         UIout(1:3,icp)=UIout(1:3,icp)+UItmp(1:3)
      enddo! loop on particles
   enddo ! loop CPs
   !$OMP END DO 
   !$OMP END PARALLEL
end subroutine ui_part_nograd

!> Induced velocity from 1 particle at 1 control point. The velocity gradient is not computed
subroutine ui_part_nograd_11(DeltaP, Alpha, RegFunction, RegParam, Ui)
   real(ReKi), dimension(3), intent(out) :: Ui          !< no side effects
   real(ReKi), dimension(3), intent(in)  :: DeltaP      !< CP-PP "control point - particle point"
   real(ReKi), dimension(3), intent(in)  :: Alpha       !< Particle intensity [m^2/s] alpha=om.dV
   integer(IntKi),           intent(in)  :: RegFunction !< 
   real(ReKi),               intent(in)  :: RegParam    !< 
   real(ReKi),dimension(3) :: C          !< Cross product of Alpha and r
   real(ReKi)              :: E          !< Exponential poart for the mollifider
   real(ReKi)              :: r3_inv     !< 
   real(ReKi)              :: rDeltaP    !< norm , distance between point and particle
   real(ReKi)              :: ScalarPart !< the part containing the inverse of the distance, but not 4pi, Mollifier
   rDeltaP=sqrt(DeltaP(1)**2+ DeltaP(2)**2+ DeltaP(3)**2)! norm
   if (rDeltaP<MINNORM) then !--- Exactly on the Singularity 
      Ui(1:3)  = 0.0_ReKi
      return
   else !--- Normal Procedure 
      C(1) = Alpha(2) * DeltaP(3) - Alpha(3) * DeltaP(2)
      C(2) = Alpha(3) * DeltaP(1) - Alpha(1) * DeltaP(3)
      C(3) = Alpha(1) * DeltaP(2) - Alpha(2) * DeltaP(1)
      select case (RegFunction) !
      case (idRegNone) ! No mollification
         r3_inv     = 1._ReKi/(rDeltaP**3)
         ScalarPart = r3_inv*fourpi_inv
      case (idRegExp) ! Exponential mollifier
         r3_inv     = 1._ReKi/(rDeltaP**3)
         E          = exp(-rDeltaP**3/RegParam**3)
         ScalarPart = (1._ReKi-E)*r3_inv*fourpi_inv
      case (idRegCompact) ! Compact support
         r3_inv     = 1._ReKi/sqrt(RegParam**6+rDeltaP**6)
         ScalarPart = r3_inv*fourpi_inv
      case default 
         print*,'[ERROR] Wrong regularization function for particles',RegFunction
         STOP
      end select 
      Ui(1:3)=C*ScalarPart
   end if ! not on singularity
end subroutine ui_part_nograd_11

!> Velocity induced by one vortex quad on nCPs Control Points
subroutine ui_quad_n1(CPs, nCPs, P1, P2, P3, P4, Gamm, RegFunction, RegParam, Uind)
   integer,                    intent(in)    :: nCPs        !< 
   real(ReKi), dimension(:,:), intent(in)    :: CPs         !< 3 x "nCPs"++
   real(ReKi), dimension(3),   intent(in)    :: P1,P2,P3,P4 !< Coordinates of vortex quadrilateral
   real(ReKi),                 intent(in)    :: Gamm
   integer(IntKi) ,            intent(in)    :: RegFunction !< Regularization model (e.g. LambOseen)
   real(ReKi),                 intent(in)    :: RegParam    !< Regularization parameter [m]
   real(ReKi), dimension(:,:), intent(inout) :: Uind        !< side effects!!! 3 x "nCPs++"
   real(ReKi), dimension(3) :: CP      !< 
   real(ReKi), dimension(3) :: Uindtmp !< 
   real(ReKi), dimension(3) :: DP1     !< 
   real(ReKi), dimension(3) :: DP2     !< 
   integer                :: icp
   ! 
   !OMP PARALLEL DEFAULT(SHARED)
   !OMP DO PRIVATE(icp,CP,Uindtmp,DP1,DP2) schedule(runtime)
   do icp=1,nCPs
      CP(1:3)=CPs(1:3,icp)
      ! 1-2 segment
      DP1=CP-P1; DP2=CP-P2; 
      call ui_seg_11 ( DP1, DP2, Gamm, RegFunction, RegParam, Uindtmp)
      Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
      ! 3-4 segment
      DP1=CP-P3; DP2=CP-P4; 
      call ui_seg_11 ( DP1, DP2, Gamm, RegFunction, RegParam, Uindtmp)
      Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
      ! 2-3 segment
      DP1=CP-P2; DP2=CP-P3; 
      call ui_seg_11 ( DP1, DP2, Gamm, RegFunction, RegParam, Uindtmp)
      Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
      ! 4-1 segment
      DP1=CP-P4; DP2=CP-P1; 
      call ui_seg_11 ( DP1, DP2, Gamm, RegFunction, RegParam, Uindtmp)
      Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
   end do  ! loop on CPs
   !OMP END DO 
   !OMP END PARALLEL
end subroutine  ui_quad_n1

end module FVW_BiotSavart
