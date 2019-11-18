MODULE FVW_VortexTools
   ! Contains Typical Tools  for vortex methods

   ! Should be *independent* of the Framework and any derived type 

   ! Only low level functions !

   use NWTC_LIBRARY

   IMPLICIT NONE

CONTAINS

   subroutine LatticeToSegments(LatticePoints, LatticeGamma, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
      real(Reki), dimension(:,:,:),    intent(in   )  :: LatticePoints  !< Points 3 x nSpan x nDepth
      real(Reki), dimension(:,:),      intent(in   )  :: LatticeGamma   !< GammaPanl  nSpan x nDepth
      real(ReKi), dimension(:,:),      intent(inout)  :: SegPoints      !< 
      integer(IntKi), dimension(:,:),  intent(inout)  :: SegConnct      !< 
      real(ReKi),     dimension(:),    intent(inout)  :: SegGamma       !< 
      integer(IntKi),                  intent(  out)  :: iHeadP         !< Index indicating where to start in SegPoints
      integer(IntKi),                  intent(  out)  :: iHeadC         !< Index indicating where to start in SegConnct
      ! Local
      integer(IntKi) :: nSpan, nDepth
      integer(IntKi) :: iSpan, iDepth
      real(ReKi) :: c1,c2
      integer(IntKi) :: iHeadP0, iseg1, iseg2, iseg3 ,iseg4  !< Index indicating where to start in SegPoints
      real(ReKi) :: Gamma12
      real(ReKi) :: Gamma41

      nSpan = size(LatticePoints,2)
      nDepth= size(LatticePoints,3)


      iHeadP0=iHeadP ! Storing

      ! We will need all the points, we flatten the point array
      ! Loop order is important
      ! Points are flattened as follows:
      !   
      !   3---6
      !   |   |
      !   2---5
      !   |   |
      !   1---4
      !
      do iDepth = 1, nDepth
         do iSpan = 1, nSpan
            SegPoints(1:3,iHeadP) = LatticePoints(1:3, iSpan, iDepth)
            iHeadP=iHeadP+1
         enddo
      enddo

      ! --- Creating segments
      ! Naming convention for point indices and segments of a panel:
      !     2---3
      !     |   |
      !     1---4
      ! We go "Panel per panel" , for a given Panel, we create
      !   - Segment 1-2
      !   - Segment 1-4 
      !   - Segment 4-3 if     the last Depth panel 
      !   - Segment 2-3 if     the last Span panel 
      !  Circulation is defined positive as follows (clockwise):
      !     2->-3
      !     ^   v
      !     1-<-4
      do iDepth = 1, nDepth-1
         do iSpan = 1, nSpan-1
            iseg1 = iHeadP0 + (iSpan-1) +(iDepth-1)*nSpan  ! Point 1
            iseg2 = iHeadP0 + (iSpan  ) +(iDepth-1)*nSpan  ! Point 2
            iseg3 = iHeadP0 + (iSpan  ) +(iDepth  )*nSpan  ! Point 3
            iseg4 = iHeadP0 + (iSpan-1) +(iDepth  )*nSpan  ! Point 4
            if (iDepth==1) then
               Gamma12 = LatticeGamma(iSpan,iDepth)
            else
               Gamma12 = LatticeGamma(iSpan,iDepth)-LatticeGamma(iSpan,iDepth-1)
            endif
            if (iSpan==1) then
               Gamma41 = LatticeGamma(iSpan,iDepth)
            else
               Gamma41 = LatticeGamma(iSpan,iDepth)-LatticeGamma(iSpan-1,iDepth)
            endif
            !print*,iseg1,iseg2,iseg3,iseg4
            ! Segment 1-2
            SegConnct(1,iHeadC) = iseg1
            SegConnct(2,iHeadC) = iseg2
            SegGamma (iHeadC  ) = Gamma12
            iHeadC=iHeadC+1
            ! Segment 1-4
            SegConnct(1,iHeadC) = iseg1
            SegConnct(2,iHeadC) = iseg4
            SegGamma (iHeadC  ) = -Gamma41
            iHeadC=iHeadC+1
            ! Segment 4-3
            if (iDepth==nDepth-1) then
               SegConnct(1,iHeadC) = iseg4
               SegConnct(2,iHeadC) = iseg3
               SegGamma (iHeadC  ) = - LatticeGamma(iSpan,iDepth)
               iHeadC=iHeadC+1
            endif
            ! Segment 2-3
            if (iSpan==nSpan-1) then
               SegConnct(1,iHeadC) = iseg2
               SegConnct(2,iHeadC) = iseg3
               SegGamma (iHeadC  ) = LatticeGamma(iSpan,iDepth)
               iHeadC=iHeadC+1
            endif
         enddo
      enddo

   end subroutine 



END MODULE FVW_VortexTools
