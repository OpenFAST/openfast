!> Biot-Savart law functions
!! NOTE: these functions should be independent of the framework types
module FVW_BiotSavart 

   use NWTC_Library, only: ReKi, IntKi, Pi, EqualRealNos
   use OMP_LIB

   implicit none

   real(ReKi),parameter :: PRECISION_UI  = epsilon(1.0_ReKi)/100 !< NOTE assuming problem of size 1
   real(ReKi),parameter :: PRECISION_EPS =  epsilon(1.0_ReKi) !< Machine Precision For the given ReKi for problems of scale 1!
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
   real(ReKi),parameter    :: fourpi     =  4.00_ReKi * ACOS(-1.0_Reki )

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
   Uind(1:3) = 0.0_ReKi
   xa=DeltaPa(1); ya=DeltaPa(2); za=DeltaPa(3)
   xb=DeltaPb(1); yb=DeltaPb(2); zb=DeltaPb(3)
   norm_a      = sqrt(xa*xa + ya*ya + za*za)
   norm_b      = sqrt(xb*xb + yb*yb + zb*zb)
   denominator  = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb) ! |r1|*|r2|*(|r1|*|r2| + r1.r2)
   
   if (denominator <= PRECISION_UI) return
   crossprod(1) = ya*zb-za*yb; crossprod(2) = za*xb-xa*zb; crossprod(3) = xa*yb-ya*xb
   norm2_orth = crossprod(1)**2 + crossprod(2)**2 + crossprod(3)**2

   ! On the singularity, Uind(1:3)=0.0_ReKi
   if (norm2_orth <= PRECISION_UI) return
   norm2_r0 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) + (za-zb)*(za-zb) 

   ! segment of zero length
   if (norm2_r0 <= PRECISION_UI) return 

   ! --- Far field TODO
   ! --- Regularization (close field)
   norm2_orth = norm2_orth/norm2_r0 ! d = (r1xr2)/r0

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
      print*,'Unknown SgmtReg', RegFunction
      STOP ! Will never happen
      Kv=1.0_ReKi !< Should be an error
   end select 

   Kv = SegGamma*fourpi_inv*Kv*(norm_a+norm_b)/denominator
   Uind(1:3) = Kv*crossprod(1:3)
   
end subroutine ui_seg_11


!> Induced velocity from a list of segments defined by Connectivity (SegConnct) and Points (SegPoints)
!! NOTE: this function has side effects and expects Uind_out to be initialized!
!! The function can compute the velocity on part of the segments and part of the control points.
!! This feature is useful if some parallelization is used, while common storage vectors are used.
subroutine ui_seg(iCPStart, iCPEnd, CPs, &
      iSegStart, iSegEnd, SegPoints, SegConnct, SegGamma,  &
      RegFunction, RegParam, Uind_out)
   real(ReKi), dimension(:,:),     intent(in)    :: CPs         !< Control points (3 x nCPs++)
   integer(IntKi),                 intent(in)    :: iCPStart    !< Index where we start in Control points array
   integer(IntKi),                 intent(in)    :: iCPEnd      !< Index where we end in Control points array
   real(ReKi), dimension(:,:),     intent(in)    :: SegPoints   !< Segment points (3 x nSegPTot)
   integer(IntKi), dimension(:,:), intent(in)    :: SegConnct   !< Connectivity, indices of segments points iSeg1, iSeg2, iDepth, iSpan
   real(ReKi), dimension(:),       intent(in)    :: SegGamma    !< Segment circulation (nSegTot)
   integer(IntKi),                 intent(in)    :: iSegStart   !< Index in SegConnct, and SegGamma where we start
   integer(IntKi),                 intent(in)    :: iSegEnd     !< Index in SegConnct, and SegGamma where we end
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
   real(ReKi)              :: CPs_icp(3)      !< 

   ! Branching based on regularization model
   ! NOTE: copy paste of code is done for optimization!
   !       The only thing changing is the part labelled "regularization"
   select case (RegFunction) 
   case ( idRegNone ) ! No vortex core 
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,CPs_icp,Uind,P1,P2,crossprod,denominator,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         Uind = 0.0_ReKi
         CPs_icp = CPs(:,icp)
         do is=iSegStart,iSegEnd ! loop on selected segments 
            P1   = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2   = SegPoints(1:3, SegConnct(2,is))
            xa=CPs_icp(1)-P1(1); ya=CPs_icp(2)-P1(2); za=CPs_icp(3)-P1(3);
            xb=CPs_icp(1)-P2(1); yb=CPs_icp(2)-P2(2); zb=CPs_icp(3)-P2(3);
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
                     Uind(1:3) = Uind(1:3) + Kv*crossprod(1:3)
                  end if
               end if
            end if
         end do ! Loop on segments
         Uind_out(1:3,icp) = Uind_out(1:3,icp)+Uind(1:3)
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL
      
   case ( idRegRankine )      ! Rankine
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,CPs_icp,Uind,P1,P2,crossprod,denominator,r_bar2,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         Uind = 0.0_ReKi
         CPs_icp = CPs(:,icp)
         do is=iSegStart,iSegEnd ! loop on selected segments 
            P1   = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2   = SegPoints(1:3, SegConnct(2,is))
            xa=CPs_icp(1)-P1(1); ya=CPs_icp(2)-P1(2); za=CPs_icp(3)-P1(3);
            xb=CPs_icp(1)-P2(1); yb=CPs_icp(2)-P2(2); zb=CPs_icp(3)-P2(3);
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
                     Uind(1:3) = Uind(1:3) + Kv*crossprod(1:3)
                  end if
               end if ! denominator size or distances too small
            end if ! 
         end do ! Loop on segments
         Uind_out(1:3,icp) = Uind_out(1:3,icp) + Uind(1:3)
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL

   case ( idRegLambOseen )      ! LambOseen
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,CPs_icp,Uind,P1,P2,crossprod,denominator,r_bar2,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb,exp_value) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         Uind = 0.0_ReKi
         CPs_icp = CPs(:,icp)
         do is=iSegStart,iSegEnd ! loop on selected segments 
            P1   = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2   = SegPoints(1:3, SegConnct(2,is))
            xa=CPs_icp(1)-P1(1); ya=CPs_icp(2)-P1(2); za=CPs_icp(3)-P1(3);
            xb=CPs_icp(1)-P2(1); yb=CPs_icp(2)-P2(2); zb=CPs_icp(3)-P2(3);
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
                     Uind(1:3) = Uind(1:3) + Kv*crossprod(1:3)
                  endif 
               end if ! denominator size or distances too small
            end if ! 
         end do ! Loop on segments
         Uind_out(1:3,icp) = Uind_out(1:3,icp) + Uind(1:3)
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL

   case ( idRegVatistas )      ! Vatistas n=2
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,CPs_icp,Uind,P1,P2,crossprod,denominator,r_bar2,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         Uind = 0.0_ReKi
         CPs_icp = CPs(:,icp)
         do is=iSegStart,iSegEnd ! loop on selected segments 
            P1   = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2   = SegPoints(1:3, SegConnct(2,is))
            xa=CPs_icp(1)-P1(1); ya=CPs_icp(2)-P1(2); za=CPs_icp(3)-P1(3);
            xb=CPs_icp(1)-P2(1); yb=CPs_icp(2)-P2(2); zb=CPs_icp(3)-P2(3);
            norm_a      = sqrt(xa*xa + ya*ya + za*za)
            norm_b      = sqrt(xb*xb + yb*yb + zb*zb)
            denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
            ! denominator size or distances too small
            if (denominator <= PRECISION_UI) cycle
            crossprod(1) = ya*zb-za*yb; crossprod(2) = za*xb-xa*zb; crossprod(3) = xa*yb-ya*xb
            norm2_orth   = crossprod(1)**2 + crossprod(2)**2 + crossprod(3)**2
            ! On the singularity, cycle
            if (norm2_orth <= PRECISION_UI) cycle
            norm2_r0     = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb) 
            ! segment of zero length
            if (norm2_r0 <= PRECISION_UI) cycle
            ! --- Far field TODO
            ! --- Regularization (close field) --- Vatistas
            norm2_orth = norm2_orth/norm2_r0 ! d = (r1xr2)/r0
            r_bar2     = norm2_orth/RegParam(is)**2
            Kv         = r_bar2/sqrt(1.0_ReKi+r_bar2**2)
            Kv         = SegGamma(is)*fourpi_inv*Kv*(norm_a+norm_b)/(denominator + MINDENOM)
            Uind(1:3)  = Uind(1:3) + Kv*crossprod(1:3)
         end do ! Loop on segments
         Uind_out(1:3,icp) = Uind_out(1:3,icp) + Uind(1:3)
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL

   case ( idRegOffset )      ! Denominator offset
      !$OMP PARALLEL default(shared)
      !$OMP do private(icp,is,CPs_icp,Uind,P1,P2,crossprod,denominator,r_bar2,Kv,norm_a,norm_b,norm2_r0,norm2_orth,xa,ya,za,xb,yb,zb) schedule(runtime)
      do icp=iCPStart,iCPEnd ! loop on CPs 
         Uind      = 0.0_ReKi
         CPs_icp = CPs(:,icp)
         do is=iSegStart,iSegEnd ! loop on selected segments 
            P1        = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
            P2        = SegPoints(1:3, SegConnct(2,is))
            xa=CPs_icp(1)-P1(1); ya=CPs_icp(2)-P1(2); za=CPs_icp(3)-P1(3);
            xb=CPs_icp(1)-P2(1); yb=CPs_icp(2)-P2(2); zb=CPs_icp(3)-P2(3);
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
                     Uind(1:3)   = Uind(1:3) + Kv*crossprod(1:3)
                  end if
               end if ! denominator size or distances too small
            end if ! 
         end do ! Loop on segments
         Uind_out(1:3,icp) = Uind_out(1:3,icp)+Uind(1:3)
      enddo ! Loop on control points
      !$OMP END DO 
      !$OMP END PARALLEL
   case default
      print*,'[ERROR] Unknown RegFunction for segment',RegFunction
      STOP
   end select
end subroutine ui_seg

!> Induced velocity from `nPart` particles at `nCPs` control points. The velocity gradient is not computed
subroutine ui_part_nograd(nCPS, CPs, nPart, Part, Alpha, RegFunction, RegParam, UIout)
   integer(IntKi),               intent(in)    :: nCPs        !< Number of control points to use (nCPs<=size(CPs,2))
   integer(IntKi),               intent(in)    :: nPart       !< Number of particles to use (nPart<=size(Part,2))
   real(ReKi), dimension(:,:),   intent(in)    :: CPs         !< Control points (3 x nCPs+)
   real(ReKi), dimension(:,:),   intent(inout) :: UIout       !< Induced velocity, with side effects! (3 x nCPs+)
   real(ReKi), dimension(:,:),   intent(in)    :: Part        !< Particle positions (3 x nPart+)
   real(ReKi), dimension(:,:),   intent(in)    :: Alpha       !< Particle intensity [m^3/s] (3 x nPart+) omega dV= alpha
   integer(IntKi),               intent(in)    :: RegFunction !< Regularization function 
   real(ReKi), dimension(:),     intent(in)    :: RegParam    !< Regularization parameter (nPart+)
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


subroutine ui_quad_src_11(CP, Sigma, xi, eta, RefPoint, R_g2p, UI)
   real(ReKi),                 intent(in)  :: Sigma      !< Source panel intensity
   real(ReKi), dimension(3),   intent(in)  :: CP         !< Control Point
   real(ReKi), dimension(3),   intent(out) :: UI         !< Induced velocity
   real(ReKi), dimension(3),   intent(in)  :: RefPoint   !< Coordinate of panel origin in ref coordinates
   real(ReKi), dimension(4),   intent(in)  :: xi         !< Panel points coordinates
   real(ReKi), dimension(4),   intent(in)  :: eta        !< Panel points  coordinates
   real(ReKi), dimension(3,3), intent(in)  :: R_g2p !< 3 x 3, global 2 panel
   real(ReKi),parameter       :: eps_quadsource=1e-6_ReKi !!!!!!!!!!!!!!!!!! !< Used if z coordinate close to zero
   real(ReKi), dimension(3,3) :: tA                         !< 
   real(ReKi)                 :: d12, d23, d34, d41         !< 
   real(ReKi)                 :: m12, m23, m34, m41         !< 
   real(ReKi)                 :: xi1,  xi2,  xi3,  xi4      !< 
   real(ReKi)                 :: eta1,  eta2,  eta3,  eta4  !< 
   real(ReKi)                 :: e1,  e2,  e3,  e4          !< 
   real(ReKi)                 :: h1,  h2,  h3,  h4          !< 
   real(ReKi)                 :: r1,  r2,  r3,  r4          !< 
   real(ReKi)                 :: RJ12, RJ23, RJ34, RJ41     !< 
   real(ReKi)                 :: TAN12, TAN23, TAN34, TAN41 !< 
   real(ReKi), dimension(3)   :: Vp                         !< 
   real(ReKi), dimension(3)   :: DP                         !< 
   real(ReKi), dimension(3)   :: DPp                        !< 
   xi1=xi(1)
   xi2=xi(2)
   xi3=xi(3)
   xi4=xi(4)
   eta1=eta(1)
   eta2=eta(2)
   eta3=eta(3)
   eta4=eta(4)
   !param that are constant for each panel, distances and slopes - The slopes can be if divided by zero NaN => security required 
   d12 = sqrt((xi2-xi1)**2+(eta2-eta1)**2)
   d23 = sqrt((xi3-xi2)**2+(eta3-eta2)**2)
   d34 = sqrt((xi4-xi3)**2+(eta4-eta3)**2)
   d41 = sqrt((xi1-xi4)**2+(eta1-eta4)**2)

   ! transform control points in panel coordinate system using matrix 
   DP(1:3) = CP(1:3)-RefPoint(1:3)
   DPp      = matmul(R_g2p, DP)           ! transfo in element coordinate system, noted x,y,z, but in fact xi eta zeta
   ! scalars
   r1 = sqrt((DPp(1)-xi1)**2 + (DPp(2)-eta1)**2 + DPp(3)**2)
   r2 = sqrt((DPp(1)-xi2)**2 + (DPp(2)-eta2)**2 + DPp(3)**2)
   r3 = sqrt((DPp(1)-xi3)**2 + (DPp(2)-eta3)**2 + DPp(3)**2)
   r4 = sqrt((DPp(1)-xi4)**2 + (DPp(2)-eta4)**2 + DPp(3)**2)
   !
   e1 = DPp(3)**2 + (DPp(1)-xi1)**2
   e2 = DPp(3)**2 + (DPp(1)-xi2)**2
   e3 = DPp(3)**2 + (DPp(1)-xi3)**2
   e4 = DPp(3)**2 + (DPp(1)-xi4)**2
   ! 
   h1 = (DPp(2)-eta1)*(DPp(1)-xi1)
   h2 = (DPp(2)-eta2)*(DPp(1)-xi2)
   h3 = (DPp(2)-eta3)*(DPp(1)-xi3)
   h4 = (DPp(2)-eta4)*(DPp(1)-xi4)
   ! Velocities in element frame 
   ! --- Log term 
   ! Security - Katz Plotkin page 608 appendix D Code 11 
   if ( r1+r2-d12<=0.0_ReKi .or. d12 <=0.0_ReKi ) then
      RJ12=0._ReKi
   else
      RJ12=1/d12 * log((r1+r2-d12)/(r1+r2+d12))
   endif

   if ( r2+r3-d23<=0.0_ReKi .or. d23 <=0.0_ReKi ) then
      RJ23=0._ReKi
   else
      RJ23=1/d23 * log((r2+r3-d23)/(r2+r3+d23))
   endif

   if ( r3+r4-d34<=0.0_ReKi .or. d34 <=0.0_ReKi ) then
      RJ34=0._ReKi
   else
      RJ34=1/d34 * log((r3+r4-d34)/(r3+r4+d34))
   endif

   if ( r4+r1-d41<=0.0_ReKi .or. d41 <=0.0_ReKi ) then
      RJ41=0._ReKi
   else
      RJ41=1/d41 * log((r4+r1-d41)/(r4+r1+d41))
   endif
   ! --- Tan term 
   ! 12
   if (EqualRealNos(xi2,xi1)) then ! Security - Hess 1962 - page 47 - bottom
      TAN12=0._ReKi
   else
      m12=(eta2-eta1)/(xi2-xi1)
      if( abs(DPp(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
         TAN12=pi*aint((signit(1.0_ReKi,(m12*e1-h1)) - signit(1.0_ReKi,(m12*e2-h2)))/2) ! Security-Hess1962-page47-top
      else
         TAN12= atan((m12*e1-h1)/(DPp(3)*r1)) - atan((m12*e2-h2)/(DPp(3)*r2))
      endif
   endif
   ! 23
   if (EqualRealNos(xi3,xi2)) then ! Security - Hess 1962 - page 47 - bottom
      TAN23=0._ReKi
   else
      m23=(eta3-eta2)/(xi3-xi2)
      if( abs(DPp(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
         TAN23=pi*aint((signit(1.0_ReKi,(m23*e2-h2)) - signit(1.0_ReKi,(m23*e3-h3)))/2) ! Security-Hess1962-page47-top
      else
         TAN23= atan((m23*e2-h2)/(DPp(3)*r2)) - atan((m23*e3-h3)/(DPp(3)*r3))
      endif
   endif 
   ! 34
   if (EqualRealNos(xi4,xi3)) then ! Security - Hess 1962 - page 47 - bottom
      TAN34=0._ReKi
   else
      m34=(eta4-eta3)/(xi4-xi3)
      if( abs(DPp(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
         TAN34=pi*aint((signit(1.0_ReKi,(m34*e3-h3)) - signit(1.0_ReKi,(m34*e4-h4)))/2) ! Security-Hess1962-page47-top
      else
         TAN34= atan((m34*e3-h3)/(DPp(3)*r3)) - atan((m34*e4-h4)/(DPp(3)*r4))
      endif
   endif
   ! 41
   if (EqualRealNos(xi1,xi4)) then ! Security - Hess 1962 - page 47 - bottom
      TAN41=0._ReKi
   else
      m41=(eta1-eta4)/(xi1-xi4)
      if( abs(DPp(3))<eps_quadsource ) then ! case where z is too small, jumps may occur 
         TAN41=pi*aint((signit(1.0_ReKi,(m41*e4-h4)) - signit(1.0_ReKi,(m41*e1-h1)))/2) ! Security-Hess1962-page47-top
      else
         TAN41= atan((m41*e4-h4)/(DPp(3)*r4)) - atan((m41*e1-h1)/(DPp(3)*r1))
      endif
   endif
   ! --- Velocity  in Panel frame
   Vp(1)= Sigma/(fourpi)*( (eta2-eta1)*RJ12 + (eta3-eta2)*RJ23 + (eta4-eta3)*RJ34 + (eta1-eta4)*RJ41 )
   Vp(2)= Sigma/(fourpi)*( (xi1-xi2)  *RJ12 + (xi2-xi3)  *RJ23 +  (xi3-xi4) *RJ34 +  (xi4-xi1) *RJ41 )
   Vp(3)= Sigma/(fourpi)*( ( TAN12 ) + ( TAN23 ) + ( TAN34 ) + ( TAN41 ) )
   ! --- Velocity in Reference frame 
   UI(1:3) = matmul(transpose(R_g2p), Vp(1:3))
end subroutine  ui_quad_src_11

!> Induced velocity by several flat quadrilateral source panels on multiple control points (CPs)
subroutine ui_quad_src_nn(CPs, Sigmas, xi, eta, RefPoint, R_g2p, UI, nCPs, nPanels)
   integer,                            intent(in)    :: nCPs
   integer,                            intent(in)    :: nPanels
   real(ReKi), dimension(3,nCPs),      intent(in)    :: CPs      !< 3 x nCPs, coordinates of the control points
   real(ReKi), dimension(3,nCPs),      intent(inout) :: UI       !< 3 x nCPs, induced velocity on control points (side effects)
   real(ReKi), dimension(nPanels),     intent(in)    :: Sigmas   !< np,       intensities of Panels
   real(ReKi), dimension(3,nPanels),   intent(in)    :: RefPoint !< 3 x np  , coordinate of panel origin in ref coordinates
   real(ReKi), dimension(4,nPanels),   intent(in)    :: xi       !< 4 x np
   real(ReKi), dimension(4,nPanels),   intent(in)    :: eta      !< 4 x np
   real(ReKi), dimension(3,3,nPanels), intent(in)    :: R_g2p    !< 3 x 3 x np, transformation matrix global to panel
   real(ReKi) :: Uind_tmp(3) !< 
   real(ReKi) :: Uind_cum(3) !< 
   integer    :: ip, icp     !< loop index
   !$OMP PARALLEL DEFAULT(SHARED)
   !$OMP DO PRIVATE(icp, Uind_cum, Uind_tmp, ip) schedule(runtime)
   do icp=1,nCPs ! loop on Control Points
      Uind_cum = 0.0_ReKi
      do ip=1,nPanels !loop on panels 
         call ui_quad_src_11(CPs(:,icp), Sigmas(ip), xi(:,ip), eta(:,ip), RefPoint(:,ip), R_g2p(:,:,ip), Uind_tmp)
         Uind_cum = Uind_cum + Uind_tmp
      enddo
      UI(1:3,icp) = UI(1:3,icp) + Uind_cum
   end do ! control points
   !$OMP END DO 
   !$OMP END PARALLEL
end subroutine ui_quad_src_nn

elemental real(ReKi) function signit(ref, val)
  real(ReKi),intent(in) ::ref
  real(ReKi),intent(in) ::val
  if ( abs(val)>PRECISION_EPS ) then
      signit = sign(ref, val)
  else
      signit = 1.0_ReKi
  endif
endfunction

end module FVW_BiotSavart
