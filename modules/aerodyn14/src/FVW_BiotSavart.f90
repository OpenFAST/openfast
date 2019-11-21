module FVW_BiotSavart 

   use NWTC_Library, only: ReKi, IntKi

   implicit none

   real(ReKi),parameter :: PRECISION_UI  = epsilon(1.0_ReKi)/100 !< NOTE assuming problem of size 1
   real(ReKi),parameter :: MIN_EXP_VALUE=-10.0_ReKi
   real(ReKi),parameter :: MINDENOM=1e-15_ReKi
   real(ReKi),parameter :: MINNORMSIMP=1e-6_ReKi

   integer(IntKi), parameter :: idSegSmoothNone       = 0
   integer(IntKi), parameter :: idSegSmoothRankine    = 1
   integer(IntKi), parameter :: idSegSmoothLambOseen  = 2
   integer(IntKi), parameter :: idSegSmoothVatistas   = 3
   integer(IntKi), parameter :: idSegSmoothCompact    = 4
   integer(IntKi), parameter :: idSegSmoothCompactDim = 5
contains


!> Induced velocity from one segment at one control points
subroutine ui_seg_11(DeltaPa, DeltaPb, Gam, SmoothModel , SmoothParam, Ui)
   implicit none 
   ! Input/output arguments 
   real(ReKi), dimension(3), intent(in) :: DeltaPa      !< 3 x 1   Pcp-P1 
   real(ReKi), dimension(3), intent(in) :: DeltaPb      !< 3 x 1   Pcp-P2
   real(ReKi), intent(in)               :: Gam    !< 
   integer, intent(in)                  :: SmoothModel   !< 
   real(ReKi), intent(in)               :: SmoothParam   !< 
   real(ReKi), dimension(3),intent(out) :: Ui           !<  No Side effects
   ! Variables declaration 
   real(ReKi),dimension(3) :: crossprod      !< 
   real(ReKi),dimension(3) :: D              !< 
   real(ReKi)              :: denom          !< 
   real(ReKi)              :: denominator    !< 
   real(ReKi)              :: h2             !< Square of h
   real(ReKi)              :: h              !< Only used by one model
   real(ReKi)              :: Kv             !< 
   real(ReKi)              :: norm_a         !< 
   real(ReKi)              :: norm_b         !< 
   real(ReKi)              :: norm2_r0        !< 
   real(ReKi)              :: norm2_crossprod !< 
   real(ReKi)              :: xa             !< 
   real(ReKi)              :: ya             !< 
   real(ReKi)              :: za             !< 
   real(ReKi)              :: xb             !< 
   real(ReKi)              :: yb             !< 
   real(ReKi)              :: zb             !< 
   real(ReKi)              :: SmoothParam2    !< Square of viscous param
   real(ReKi)              :: exp_value      !< 
   real(ReKi),parameter    :: fourpi_inv =  0.25_ReKi / ACOS(-1.0_Reki )
   ! 
   xa=DeltaPa(1)
   ya=DeltaPa(2)
   za=DeltaPa(3)
   xb=DeltaPb(1)
   yb=DeltaPb(2)
   zb=DeltaPb(3)
   !--- Simple test if on an extremity point
   if(abs(xa)+abs(ya)+abs(za)<MINNORMSIMP) then
      Ui(1:3)=0.0_ReKi
   elseif(abs(xb)+abs(yb)+abs(zb)<MINNORMSIMP) then
      Ui(1:3)=0.0_ReKi
   else
      norm2_r0 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb)
      norm_b = sqrt(xb*xb + yb*yb + zb*zb)
      norm_a = sqrt(xa**2 + ya**2 + za**2)
      !denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
      denom       =               (norm_a*norm_b + xa*xb+ya*yb+za*zb)
      denominator = norm_a*norm_b*denom
      crossprod(1) = ya*zb-za*yb
      crossprod(2) = za*xb-xa*zb
      crossprod(3) = xa*yb-ya*xb
      norm2_crossprod=crossprod(1)**2+ crossprod(2)**2+ crossprod(3)**2
      if (norm2_r0<PRECISION_UI .or. denominator<PRECISION_UI .or. norm2_crossprod<PRECISION_UI) then
         !--- Exactly on the Singularity, velocity is zero, but gradient may be constant (is is wanted?)
         if(norm2_r0<PRECISION_UI) then
            !print*,'Exit3, because of segment too small'
         else
            if(denominator<PRECISION_UI) then
               !print*,'Exit3, because of denom'
            endif
            if(norm2_crossprod<PRECISION_UI) then
               !print*,'Exit3, because of crossprod'
            endif
         endif
         Ui(1:3)=0.0_ReKi
         if( SmoothModel==0) then 
            return
         end if
      else
         !denominator=denominator+MINDENOM
         denominator=denominator
         ! TODO far distance
         !--- Normal Procedure 
         ! smooth model  
         select case (SmoothModel) !
         case ( idSegSmoothNone ) !! No vortex core model
            Kv=1.0_ReKi
         case ( idSegSmoothRankine ) !!Rankine - t<=>rc
            ! orthogonal distance r1xr2/r0 
            h2 = norm2_crossprod/ norm2_r0 
            SmoothParam2=SmoothParam**2
            if (h2<SmoothParam2) then 
               Kv=h2/SmoothParam2
            else
               Kv=1.0_ReKi 
            end if 
         case ( idSegSmoothLambOseen ) !!Lamb-Oseen - vsic_param<=>rc 
            ! orthogonal distance r1xr2/r0 
            h2 = norm2_crossprod/ norm2_r0 
            SmoothParam2=SmoothParam**2
            exp_value = -1.25643_ReKi*(h2)/(SmoothParam2)
            if(exp_value<MIN_EXP_VALUE) then
               Kv = 1.0_ReKi
            else
               Kv = 1.0_ReKi-exp(exp_value)
            endif
         case ( idSegSmoothVatistas ) !!Vatistas n=2 - SmoothParam<=>rc 
            ! orthogonal distance r1xr2/r0 
            h2 = norm2_crossprod/ norm2_r0 
            SmoothParam2=SmoothParam**2
            ! h = (norm_a+norm_b)/2; 
            Kv = h2/sqrt(SmoothParam2**2+h2**2)
         case ( idSegSmoothCompact ) !!Cut-off radius as function of length: delta^2*norm(r0)^2  
            h2 = norm2_crossprod/ norm2_r0 
            Kv=1.0_ReKi
            ! delta*norm(r0)^2  
            denominator=denominator+SmoothParam*norm2_r0
            SmoothParam2=SmoothParam**2 ! TODO
         case ( idSegSmoothCompactDim ) !!Cut-off radius dimension fixed: (delta l_0)^2  
            h2 = norm2_crossprod/ norm2_r0 
            Kv=1.0_ReKi
            ! (delta l_0)^2  
            denominator=denominator+SmoothParam
            SmoothParam2=SmoothParam**2 ! TODO
         case ( 33 ) !!Vatistas n=2 - SmoothParamt<=>rc  
            ! See Hoydonck 2012. and Script MainRingSensitivityN
            ! Not matture enough
            h = (norm_a+norm_b)/2.0_ReKi
            SmoothParam2=SmoothParam**2
            if(h<2*sqrt(norm2_r0)) then
               ! use Vatistas n=2, normal one
               h2 = norm2_crossprod/norm2_r0
            else
               h2 = 1._ReKi ! TODO
            endif
            Kv = (h2)/sqrt(SmoothParam2**2+h2**2)
         case DEFAULT
            Kv=1.0_ReKi
         end select 

         Kv=Gam*fourpi_inv*Kv*(norm_a+norm_b)/denominator
         Ui(1:3) = Kv*crossprod(1:3)
      end if ! denominator size or distances too small
   end if ! 
   ! printf("!4.3f !4.3f !4.3f !4.3f !4.3fNewLine",Uout(1),Uout(2),Uout(3),Kv,denominator);  
end subroutine ui_seg_11


!> Induced velocity from a list of segments defined by Connectivity (SegConnct) and Points (SegPoints)
!! NOTE: this funciton has side effects and expect Uind_out to be initialized!
!! The function can compute the velocity on part of the segments and part of the control points.
!! This feature is useful if some parallelization is used, while common storage vectors are used.
!!
subroutine ui_seg(iCPStart, iCPEnd, nCPsTot, CPs, &
      iSegStart, iSegEnd, nSegTot, nSegPTot, SegPoints, SegConnct, SegGamma,  &
      SmoothModel, SmoothParam, Uind_out)
   use OMP_LIB
   implicit none
   real(ReKi), dimension(3,nCPsTot), intent(in)     :: CPs            !< Control points
   integer(IntKi), intent(in)                       :: iCPStart       !< Index where we start in Control points array
   integer(IntKi), intent(in)                       :: iCPEnd         !< Index where we end in Control points array
   integer(IntKi), intent(in)                       :: nCPsTot        !< Total number of control points
   real(ReKi), dimension(3,nSegPTot), intent(in)    :: SegPoints      !< Segment points
   integer(IntKi), dimension(2,nSegTot), intent(in) :: SegConnct      !< Connectivity, indices of segments points iSeg1, iSeg2
   real(ReKi), dimension(nSegTot), intent(in)       :: SegGamma       !< Segment circulation
   integer(IntKi),intent(in)                        :: iSegStart      !< Index in SegConnct, and SegGamma where we start
   integer(IntKi),intent(in)                        :: iSegEnd        !< Index in SegConnct, and SegGamma where we end
   integer(IntKi), intent(in)                       :: nSegTot        !< Total number of segments
   integer(IntKi), intent(in)                       :: nSegPTot       !< Total number of segment points
   integer(IntKi), intent(in)                       :: SmoothModel !< Smooth model 
   real(ReKi), dimension(nSegTot), intent(in)       :: SmoothParam !< Smooth parameter
   real(ReKi), dimension(3,nCPsTot), intent(inout)  :: Uind_out       !< Induced velocity vector - Side effects!!!
   ! Variables
   integer(IntKi) :: icp,is
   ! Part of argument arrays
   real(ReKi), dimension(3) :: Uind           !< 
   real(ReKi), dimension(3) :: P1             !< Point1 of a given segment
   real(ReKi), dimension(3) :: P2             !< Point2 of a given segment
   real(ReKi), dimension(3) :: DP1, DP2       !<
   real(ReKi)               :: norm2_r0       !<
   real(ReKi)               :: Gam            !<
   real(ReKi)               :: SmoothParam1    !< 
   ! Variables declaration 
   real(ReKi),dimension(3) :: crossprod       !< 
   real(ReKi)              :: denom           !< 
   real(ReKi)              :: denominator     !< 
   real(ReKi)              :: h2              !< Square of h
   real(ReKi)              :: h               !< Only used by one model
   real(ReKi)              :: Kv              !< 
   real(ReKi)              :: norm_a          !< 
   real(ReKi)              :: norm_b          !< 
   real(ReKi)              :: norm2_crossprod !< 
   real(ReKi)              :: xa              !< 
   real(ReKi)              :: ya              !< 
   real(ReKi)              :: za              !< 
   real(ReKi)              :: xb              !< 
   real(ReKi)              :: yb              !< 
   real(ReKi)              :: zb              !< 
   real(ReKi)              :: SmoothParam2     !< Square of viscous param
   real(ReKi)              :: exp_value       !< 
   real(ReKi),parameter    :: fourpi_inv =  0.25_ReKi / ACOS(-1.0_Reki )

   !OMP PARALLEL default(shared)
   !OMP do private(&
   !OMP& icp,is,Uind,P1,P2,DP1,DP2,norm2_r0,Gam,SmoothParam,&
   !OMP& crossprod,denom,denominator,h2,h,Kv,norm_a,norm_b,norm2_crossprod,xa,ya,za,xb,yb,zb,SmoothParam2,exp_value&
   !OMP& ) schedule(runtime)
   ! loop on CPs 
   do icp=iCPStart,iCPEnd
      do is=iSegStart,iSegEnd ! loop on selected segments 
         ! Kind of arguments for Segment 11
         Uind=0.0_ReKi
         P1         = SegPoints(1:3, SegConnct(1,is)) ! Segment extremity points
         P2         = SegPoints(1:3, SegConnct(2,is))
         DP1        = CPs(:,icp)-P1;
         DP2        = CPs(:,icp)-P2
         Gam        = SegGamma(is)
         SmoothParam1= SmoothParam(is)
         ! --- inlining of ui_seg_11
         xa=DP1(1)
         ya=DP1(2)
         za=DP1(3)
         xb=DP2(1)
         yb=DP2(2)
         zb=DP2(3)
         !--- Simple test if on an extremity point
         if(abs(xa)+abs(ya)+abs(za)<MINNORMSIMP) then
            !Uind(1:3)=0.0_ReKi
         elseif(abs(xb)+abs(yb)+abs(zb)<MINNORMSIMP) then
            !Uind(1:3)=0.0_ReKi
         else
            norm2_r0 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb) +(za-zb)*(za-zb)
            norm_b = sqrt(xb*xb + yb*yb + zb*zb)
            norm_a = sqrt(xa**2 + ya**2 + za**2)
            !denominator = norm_a*norm_b*(norm_a*norm_b + xa*xb+ya*yb+za*zb)
            denom       =               (norm_a*norm_b + xa*xb+ya*yb+za*zb)
            denominator = norm_a*norm_b*denom
            crossprod(1) = ya*zb-za*yb
            crossprod(2) = za*xb-xa*zb
            crossprod(3) = xa*yb-ya*xb
            norm2_crossprod=crossprod(1)**2+ crossprod(2)**2+ crossprod(3)**2
            ! check for singularity  
            if (denominator<PRECISION_UI .or. norm2_crossprod<PRECISION_UI) then
               !--- Exactly on the Singularity, velocity is zero
               !Uind(1:3)=0.0_ReKi
            else
               denominator=denominator+MINDENOM
               ! TODO far distance
               !--- Normal Procedure 
               ! smooth model  
               h2 = 0.0_ReKi
               select case ((SmoothModel)) !
               case ( idSegSmoothNone ) !! No vortex core model
                  Kv=1.0_ReKi
               case ( idSegSmoothRankine ) !!Rankine - t<=>rc
                  ! orthogonal distance r1xr2/r0 
                  h2          = norm2_crossprod/ norm2_r0
                  SmoothParam2 = SmoothParam1**2
                  if (h2<SmoothParam2) then 
                     Kv=h2/SmoothParam2
                  else
                     Kv=1.0_ReKi 
                  end if 
               case ( idSegSmoothLambOseen ) !!Lamb-Oseen - vsic_param<=>rc 
                  ! orthogonal distance r1xr2/r0 
                  h2          = norm2_crossprod/ norm2_r0
                  SmoothParam2 = SmoothParam1**2
                  exp_value = -1.25643_ReKi*(h2)/(SmoothParam2)
                  if(exp_value<MIN_EXP_VALUE) then
                     Kv = 1.0_ReKi
                  else
                     Kv = 1.0_ReKi-exp(exp_value)
                  endif
               case ( idSegSmoothVatistas ) !!Vatistas n=2 - SmoothParam<=>rc 
                  ! orthogonal distance r1xr2/r0 
                  h2          = norm2_crossprod/ norm2_r0
                  SmoothParam2 = SmoothParam1**2
                  ! h = (norm_a+norm_b)/2; 
                  Kv = h2/sqrt(SmoothParam2**2+h2**2)
               case ( idSegSmoothCompact ) !!Cut-off radius as function of length: delta^2*norm(r0)^2  
                  Kv=1.0_ReKi
                  denominator=denominator+SmoothParam1*SmoothParam1*norm2_r0
               case ( idSegSmoothCompactDim ) !!Cut-off radius dimension fixed: (delta l_0)^2  
                  Kv=1.0_ReKi
                  denominator=denominator+SmoothParam1*SmoothParam1
               case ( 33 ) !!Vatistas n=2 - SmoothParamt<=>rc  
                  ! See Hoydonck 2012
                  h           = (norm_a+norm_b)/2.0_ReKi
                  SmoothParam2 = SmoothParam1**2
                  if(h<2*sqrt(norm2_r0)) then
                     ! use Vatistas n=2, normal one
                     h2 = norm2_crossprod/norm2_r0
                  end if
                  Kv = (h2)/sqrt(SmoothParam2**2+h2**2)
               case default
                  Kv=1.0_ReKi
               end select 

               Kv=Gam*fourpi_inv*Kv*(norm_a+norm_b)/denominator
               Uind(1:3) = Kv*crossprod(1:3)
               !print "(A,3F21.16)","ui   :",Uind(1:3)
            end if ! denominator size or distances too small
         end if ! 
         ! --- END inlining of Segment 11 
         Uind_out(1:3,icp) = Uind_out(1:3,icp)+Uind(1:3)
      end do ! Loop on segments
   enddo ! Loop on control points
   !OMP END DO 
   !OMP END PARALLEL
end subroutine  


!> Velocity induced by one vortex quad on nCPs Control Points
subroutine ui_quad_n1(CPs, nCPs, P1, P2, P3, P4, Gamm, SmoothModel, SmoothParam, Uind)
   ! Arguments declarations 
   integer, intent(in)                        :: nCPs        !< 
   real(ReKi), dimension(3,nCPs), intent(in)  :: CPs         !< 
   real(ReKi), dimension(3), intent(in)       :: P1,P2,P3,P4 !< Coordinates of vortex quadrilateral
   real(ReKi),                    intent(in)  :: Gamm
   integer(IntKi) , intent(in)                :: SmoothModel !< Most likely should be 0
   real(ReKi),                    intent(in)  :: SmoothParam !< 
   real(ReKi), dimension(3,nCPs), intent(out) :: Uind        !< no side effects!!!
   ! Variable declarations 
   real(ReKi), dimension(3) :: CP           !< 
   real(ReKi), dimension(3) :: Uindtmp        !< 
   real(ReKi), dimension(3) :: DP1          !< 
   real(ReKi), dimension(3) :: DP2          !< 
   integer                :: icp
   ! 
   !OMP PARALLEL DEFAULT(SHARED)
   !OMP DO PRIVATE(icp,CP,Uindtmp,DP1,DP2) schedule(runtime)
   do icp=1,nCPs
      CP(1:3)=CPs(1:3,icp)
      ! 1-2 segment
      DP1=CP-P1; DP2=CP-P2; 
      call ui_seg_11 ( DP1, DP2, Gamm, SmoothModel, SmoothParam, Uindtmp)
      Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
      ! 3-4 segment
      DP1=CP-P3; DP2=CP-P4; 
      call ui_seg_11 ( DP1, DP2, Gamm, SmoothModel, SmoothParam, Uindtmp)
      Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
      ! 2-3 segment
      DP1=CP-P2; DP2=CP-P3; 
      call ui_seg_11 ( DP1, DP2, Gamm, SmoothModel, SmoothParam, Uindtmp)
      Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
      ! 4-1 segment
      DP1=CP-P4; DP2=CP-P1; 
      call ui_seg_11 ( DP1, DP2, Gamm, SmoothModel, SmoothParam, Uindtmp)
      Uind(1:3,icp) = Uind(1:3,icp)+Uindtmp(1:3)
   end do  ! loop on CPs
   !OMP END DO 
   !OMP END PARALLEL
end subroutine  ui_quad_n1

end module FVW_BiotSavart
