module FVW_VortexTools
   ! Contains Typical Tools  for vortex methods
   ! 
   ! Should be *independent* of the Framework and any derived type 
   !    apart from the types defined here: Tree, Seg and Part
   !
   ! Only low-to-mid level functions!

   use NWTC_LIBRARY, only: ReKi, IntKi, num2lstr, ErrID_Fatal, ErrID_None, EqualRealNos !< Not desired
   use NWTC_LIBRARY, only: InterpArray
   use FVW_BiotSavart, only: fourpi_inv

   implicit none

   logical,parameter :: DEV_VERSION_VT = .False.

   ! --- TIC TOC MODULE
   ! Parameters
   logical, parameter    :: bSilentTicToc =.false.   !< If set to true, will only show times bigger than SmartTicTocVal
   logical, parameter    :: bSmartTicToc =.true.     !< If set to true, will only show times bigger than SmartTicTocVal
   integer, parameter    :: nmax_store=40            !< maximum number of storage in tic toc stack
   integer, parameter    :: nmax_label=46            !< maximum length of label
   real(ReKi), parameter :: SmartTicTocVal=0.02_ReKi !< duration below this value will not be displayed (if bSmartTicToc is True)
   ! Reals
   real(ReKi), save :: start_time,finish_time
   ! Tic toc stack arrays
   integer, save                                          :: npos = 0    ! < position on the stack
   integer, dimension(8), save                            :: time_array
   integer, dimension(nmax_store,8), save                 :: start_arrays = 0
   character(len=nmax_label), dimension(nmax_store), save :: labels
   ! --- END TIC TOC MODULE

   ! Tree parameters
   integer, parameter :: IK1 = selected_int_kind(1) ! to store particle branch number (from 1 to 8)
   integer,parameter :: M0 = 1, M1_1=2, M1_2=3, M1_3=4, M2_11=5, M2_21=6, M2_22=7, M2_31=8, M2_32=9, M2_33=10  ! For moment coefficients
   integer,parameter :: M0_000 = 1
   integer,parameter :: M1_100 = 2
   integer,parameter :: M1_010 = 3
   integer,parameter :: M1_001 = 4

   !> 
   type T_VPart
      real(ReKi), dimension(:,:), pointer :: P           =>null() 
      real(ReKi), dimension(:,:), pointer :: Alpha       =>null() 
      real(ReKi), dimension(:),   pointer :: RegParam    =>null() 
      integer(IntKi)                      :: RegFunction =-1
      integer(IntKi)                      :: n =-1
   end type T_VPart

   !>
   type T_VSgmt
      real(ReKi), dimension(:,:), pointer :: SP           =>null()
      integer(IntKi), dimension(:,:), pointer :: SConnct      =>null()
      real(ReKi), dimension(:), pointer :: SGamma       =>null()
      real(ReKi), dimension(:),   pointer :: RegParam    =>null()
      integer(IntKi)                      :: RegFunction =-1
      integer(IntKi)                      :: n =-1
   end type T_VSgmt

   !> The node type is recursive and is used to make a chained-list of nodes for the tree
   type T_Node
      real(ReKi)                         :: radius !< Typical dimension of a cell (max of x,y,z extent)
      real(ReKi),dimension(3)            :: center !< Used to store first the geometric center, then the vorticity center
      real(ReKi),dimension(3,10)         :: Moments
      integer,dimension(:),pointer       :: iPart=>null()  !< indexes of particles stored in this node
      integer,dimension(:),pointer       :: leaves=>null()  ! NOTE: leaves are introduced to save memory
      type(T_Node),dimension(:), pointer :: branches =>null()
      integer                            :: nPart = -1  ! Number of particles in branches and leaves of this node
   end type T_Node

   !> The type tree contains some basic data, a chained-list of nodes, and a pointer to the Particle data that were used
   type T_Tree
      type(T_VPart)  :: Part            !< Storage for all particles
      type(T_VSgmt)   :: Seg             !< Storage for all segments
      integer       :: iStep =-1       !< Time step at which the tree was built
      logical       :: bGrown =.false. !< Is the tree build
      real(ReKi)    :: DistanceDirect
      type(T_Node)  :: Root            !< Contains the chained-list of nodes
   end type T_Tree

   interface cut_tree
      module procedure cut_tree_parallel ; ! to switch between parallel and rec easily
   end interface
   interface cut_tree_segment
      module procedure cut_tree_segment_parallel ; ! to switch between parallel and rec easily
   end interface

contains

   !> Returns the discrete trailed vorticity (\Delta Gamma) based on a "bound" circulation Gamma
   !!    \Delta\Gamma_i  = \int_{r i-1}^{r i+1} d\Gamma/dr *dr  [m^2/s]
   !! NOTE: this is not gamma = d\Gamma/dr  \approx \Delta\Gamma /\Delta r 
   subroutine GammaTrailed(n, Gamma_b, Gamma_t)
      integer(IntKi),             intent(in)    :: n       !< number of points along the span
      real(Reki), dimension(n),   intent(in   ) :: Gamma_b !< Bound circulation
      real(ReKi), dimension(n+1), intent(out)   :: Gamma_t !< Trailed circulation
      integer(IntKi) :: i
      Gamma_t(1)   =   Gamma_b(1)
      Gamma_t(n+1) = - Gamma_b(n)
      do i = 2, n
         Gamma_t(i) = Gamma_b(i)-Gamma_b(i-1)
      enddo
   endsubroutine 

   !> Curvilinear distance along a set of connected points.
   subroutine CurvilinearDist(n, P, s)
      integer(IntKi),             intent(in)    :: n !< number of points 
      real(Reki), dimension(3,n), intent(in)    :: P !< Point cordinates
      real(ReKi), dimension(n),   intent(out)   :: s !< Curvilinear coordinate
      integer(IntKi) :: i
      s(1) = 0
      do i = 2, n
         s(i) = s(i-1) + sqrt((P(1,i)-P(1,i-1))**2 + (P(2,i)-P(2,i-1))**2 + (P(3,i)-P(3,i-1))**2)
      enddo
   endsubroutine 

   !> Place Tip and Root vorticity according to Gamma distribution
   !! Attempts to preserve the first moment of vorticity in the 15% region of tip and root
   !! Estimtes the regularization parameters based on the length of the tip region..
   subroutine PlaceTipRoot(n, Gamma_b, xV, Eps, iRoot, iTip, Gamma_max, EpsTip, EpsRoot)
      integer(IntKi),               intent(in)  :: n           !< number of vortex points
      real(ReKi), dimension(n)    , intent(in)  :: Gamma_b     !< Wake panel circulation
      real(Reki), dimension(3,n+1), intent(in)  :: xV          !< Vortex point nodal coordinates
      real(Reki), dimension(3,n  ), intent(in)  :: Eps         !< Vortex panels epsilon
      integer(IntKi),               intent(inout) :: iRoot, iTip !< Index of tip and root vortex, if <0, they are computed
      real(ReKi),                   intent(out) :: Gamma_max
      real(Reki),                   intent(out) :: EpsTip      !< Regularization of tip
      real(Reki),                   intent(out) :: EpsRoot     !< Regularization of root
      real(ReKi), dimension(n+1) :: Gamma_t
      real(ReKi), dimension(n+1) :: s
      real(ReKi) :: rRoot, rTip, lTip, lRoot
      integer(IntKi) :: i10, i90, iTipPanel,iRootPanel

      if (n>2) then
         Gamma_max = maxval(Gamma_b,1)
         call CurvilinearDist(n+1, xV, s)
         if (iTip<0) then

            ! Find position of tip and root
            call GammaTrailed(n, Gamma_b, Gamma_t)

            ! If circulation is constant then use first and last
            if(sum(abs(Gamma_t(:)))/(n+1)<1e-6) then
               iTip =n+1
               iRoot=1
            else
               i10 = minloc(abs(s(:)-0.15*s(n+1)),1)
               i90 = minloc(abs(s(:)-0.85*s(n+1)),1)

               rTip  = sum(Gamma_t(i90:)  * (s(i90:)))/ sum(Gamma_t(i90:))
               rRoot = sum(Gamma_t(1:i10) * s(1:i10)) / sum(Gamma_t(1:i10))
               iTip  = minloc(abs(rTip - s), 1) ! NOTE:  not accurate since epsilon has one dimension less..
               iRoot = minloc(abs(rRoot - s), 1)
               iTip  = max(min(iTip,n+1), i90)
               iRoot = min(max(iRoot,1) , i10)
            endif
         endif
         rTip  = s(iTip)
         rRoot = s(iRoot)
         iTipPanel  = max(min(iTip,n), 1)
         iRootPanel = min(max(iRoot,1),n)
         ! Mean regularization at the tip and root
         EpsTip  = sum(Eps(1,iTipPanel:))  /(n-iTipPanel+1)
         EpsRoot = sum(Eps(1,1:iRootPanel))/(iRootPanel)
         ! Scaling based on the "length" of the vortex, this will need further tuning
         lTip  = (s(n+1)-rTip )/3.14 ! Approximate radius if the tip has done half a turn
         lRoot = (rRoot       )/3.14
         EpsTip  = 1.3*(lTip+EpsTip) ! Tuning factors
         EpsRoot = 1.7*(lRoot+EpsRoot)
      else
         iTip =n+1
         iRoot=1
         EpsTip   = Eps(1,iTip-1)
         EpsRoot  = Eps(1,iRoot)
      endif
   endsubroutine PlaceTipRoot


   !> Flatten/ravel a 3D grid of vectors (each of size n)
   subroutine FlattenValues(GridValues, FlatValues, iHeadP)
      real(Reki), dimension(:,:,:,:),  intent(in   )  :: GridValues  !< Grid values n x nx x ny x nz
      real(ReKi), dimension(:,:),      intent(  out)  :: FlatValues  !< Flat values n x (nx x ny x nz)
      integer(IntKi),                  intent(inout)  :: iHeadP      !< Index indicating where to start in Values
      integer(IntKi) :: i,j,k
      do k = 1, size(GridValues,4)
         do j = 1, size(GridValues,3)
            do i = 1, size(GridValues,2)
               FlatValues(:,iHeadP) = GridValues(:, i, j, k)
               iHeadP=iHeadP+1
            enddo
         enddo
      enddo
   endsubroutine FlattenValues

   !> Flatten a 3D grid of vectors (each of size n)
   subroutine DeflateValues(FlatValues, GridValues, iHeadP)
      real(ReKi), dimension(:,:),      intent(in   )  :: FlatValues  !< Flat values n x (nx x ny x nz)
      real(Reki), dimension(:,:,:,:),  intent(  out)  :: GridValues  !< Grid values n x nx x ny x nz
      integer(IntKi),                  intent(inout)  :: iHeadP      !< Index indicating where to start in Values
      integer(IntKi) :: i,j,k
      do k = 1, size(GridValues,4)
         do j = 1, size(GridValues,3)
            do i = 1, size(GridValues,2)
               GridValues(:, i, j, k) = FlatValues(:,iHeadP) 
               iHeadP=iHeadP+1
            enddo
         enddo
      enddo
   endsubroutine DeflateValues

   subroutine VecToLattice(PointVectors, iDepthStart, LatticeVectors, iHeadP)
      real(Reki), dimension(:,:),        intent(in   )  :: PointVectors   !< nVal x n
      integer(IntKi),                    intent(in   )  :: iDepthStart    !< Start index for depth dimension
      real(ReKi), dimension(:,:,:),      intent(inout)  :: LatticeVectors !< nVal x nSpan x nDepth
      integer(IntKi),                    intent(inout)  :: iHeadP         !< Index indicating where to start in PointVectors
      integer(IntKi) :: iSpan, iDepth
      do iDepth = iDepthStart,  size(LatticeVectors,3)
         do iSpan = 1, size(LatticeVectors,2)
            LatticeVectors(:, iSpan, iDepth) = PointVectors(:, iHeadP)
            iHeadP=iHeadP+1
         enddo
      enddo
   end subroutine

   subroutine VecToLattice2D(PointVectors, LatticeVals, iHeadP)
      real(Reki), dimension(:,:),        intent(in   )  :: PointVectors   !< nVal x n
      real(ReKi), dimension(:,:),        intent(inout)  :: LatticeVals    !< nVal x nSpan
      integer(IntKi),                    intent(inout)  :: iHeadP         !< Index indicating where to start in PointVectors
      integer(IntKi) :: iSpan
      do iSpan = 1, size(LatticeVals,2)
         LatticeVals(:, iSpan) = PointVectors(:, iHeadP)
         iHeadP=iHeadP+1
      enddo
   end subroutine


   subroutine LatticeToPoints(LatticePoints, iDepthStart, Points, iHeadP)
      real(Reki), dimension(:,:,:),    intent(in   )  :: LatticePoints  !< Points 3 x nSpan x nDepth
      integer(IntKi),                  intent(in   )  :: iDepthStart    !< Start index for depth dimension
      real(ReKi), dimension(:,:),      intent(inout)  :: Points         !< 
      integer(IntKi),                  intent(inout)  :: iHeadP         !< Index indicating where to start in Points
      ! Local
      integer(IntKi) :: iSpan, iDepth
      ! Points are flattened as follows: (Loop order is important)
      !   
      !   3---6
      !   |   |
      !   2---5
      !   |   |
      !   1---4
      !
      do iDepth = iDepthStart,  size(LatticePoints,3)
         do iSpan = 1, size(LatticePoints,2)
            Points(1:3,iHeadP) = LatticePoints(1:3, iSpan, iDepth)
            iHeadP=iHeadP+1
         enddo
      enddo

   endsubroutine LatticeToPoints

   subroutine LatticeToPoints2D(LatticePoints, Points, iHeadP)
      real(Reki), dimension(:,:),     intent(in   )  :: LatticePoints  !< Points 3 x nSpan 
      real(ReKi), dimension(:,:),      intent(inout)  :: Points         !< 
      integer(IntKi),                  intent(inout)  :: iHeadP         !< Index indicating where to start in Points
      ! Local
      integer(IntKi) :: iSpan
      do iSpan = 1, size(LatticePoints,2)
         Points(1:3,iHeadP) = LatticePoints(1:3, iSpan)
         iHeadP=iHeadP+1
      enddo
   endsubroutine LatticeToPoints2D


   subroutine LatticeToSegments(LatticePoints, LatticeGamma, LatticeEpsilon, iDepthStart, SegPoints, SegConnct, SegGamma, SegEpsilon, iHeadP, iHeadC, bShedVorticity, bShedLastVorticity, bHackEpsilon )
      real(Reki), dimension(:,:,:),    intent(in   )  :: LatticePoints  !< Points  3 x nSpan x nDepth
      real(Reki), dimension(:,:),      intent(in   )  :: LatticeGamma   !< GammaPanl   nSpan x nDepth
      real(Reki), dimension(:,:,:),    intent(in   )  :: LatticeEpsilon !< EpsPanl 3 x nSpan x nDepth (one per dimension)
      integer(IntKi),                  intent(in   )  :: iDepthStart    !< Start index for depth dimension
      real(ReKi), dimension(:,:),      intent(inout)  :: SegPoints      !< 
      integer(IntKi), dimension(:,:),  intent(inout)  :: SegConnct      !< 
      real(ReKi),     dimension(:),    intent(inout)  :: SegGamma       !< 
      real(ReKi),     dimension(:),    intent(inout)  :: SegEpsilon     !< 
      integer(IntKi),                  intent(inout)  :: iHeadP         !< Index indicating where to start in SegPoints
      integer(IntKi),                  intent(inout)  :: iHeadC         !< Index indicating where to start in SegConnct
      logical       ,                  intent(in   )  :: bShedVorticity !< Shed vorticity is included if true
      logical       ,                  intent(in   )  :: bShedLastVorticity !< Shed the last vorticity segment if true
      logical       ,                  intent(in   )  :: bHackEpsilon   !< Unfortunate fix so that tip and root vortex have different epsilon for FW
      ! Local
      integer(IntKi) :: nSpan, nDepth
      integer(IntKi) :: iSpan, iDepth
      integer(IntKi) :: iHeadP0, iseg1, iseg2, iseg3 ,iseg4  !< Index indicating where to start in SegPoints
      real(ReKi) :: Gamma12, Eps12
      real(ReKi) :: Gamma41, Eps41

      nSpan = size(LatticePoints,2)
      nDepth= size(LatticePoints,3)


      iHeadP0=iHeadP ! Storing
      ! --- Flattening LatticePoints into SegPoints array, and increment iHeadP
      ! We will need all the points, we flatten the point array
      call LatticeToPoints(LatticePoints, iDepthStart, SegPoints, iHeadP)

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
      do iDepth = iDepthStart, nDepth-1
         do iSpan = 1, nSpan-1
            iseg1 = iHeadP0 + (iSpan-1) +(iDepth-1-iDepthStart+1)*nSpan  ! Point 1
            iseg2 = iHeadP0 + (iSpan  ) +(iDepth-1-iDepthStart+1)*nSpan  ! Point 2
            iseg3 = iHeadP0 + (iSpan  ) +(iDepth  -iDepthStart+1)*nSpan  ! Point 3
            iseg4 = iHeadP0 + (iSpan-1) +(iDepth  -iDepthStart+1)*nSpan  ! Point 4
            if (iDepth==iDepthStart) then
               Gamma12 = LatticeGamma(iSpan,iDepth)
               Eps12   = LatticeEpsilon(1,iSpan,iDepth) ! Using epsilon x for seg12&43. TODO might change in the future
            else
               Gamma12 = LatticeGamma(iSpan,iDepth)-LatticeGamma(iSpan,iDepth-1)
               Eps12   = (LatticeEpsilon(1,iSpan,iDepth)+LatticeEpsilon(1,iSpan,iDepth-1))/2.0_ReKi
            endif
            if (iSpan==1) then
               Gamma41 = LatticeGamma(iSpan,iDepth)
               if (bHackEpsilon) then
                  Eps41   = LatticeEpsilon(2,iSpan,iDepth) ! Using epsilon y for seg41 hacked
               else
                  Eps41   = LatticeEpsilon(3,iSpan,iDepth) ! Using epsilon z for seg23&41. TODO might change in the future
               endif
            else
               Gamma41 = LatticeGamma(iSpan,iDepth)-LatticeGamma(iSpan-1,iDepth)
               Eps41   = (LatticeEpsilon(3,iSpan,iDepth)+LatticeEpsilon(3,iSpan-1,iDepth))/2.0_ReKi
            endif
            ! Segment 1-2
            if (bShedVorticity) then
               SegConnct(1,iHeadC) = iseg1
               SegConnct(2,iHeadC) = iseg2
               SegConnct(3,iHeadC) = iDepth
               SegConnct(4,iHeadC) = iSpan
               SegGamma  (iHeadC ) = Gamma12
               SegEpsilon(iHeadC ) = Eps12
               iHeadC=iHeadC+1
            endif
            ! Segment 1-4
            SegConnct(1,iHeadC) = iseg1
            SegConnct(2,iHeadC) = iseg4
            SegConnct(3,iHeadC) = iDepth
            SegConnct(4,iHeadC) = iSpan
            SegGamma (iHeadC  ) = -Gamma41
            SegEpsilon(iHeadC ) = Eps41
            iHeadC=iHeadC+1
            ! Segment 4-3
            if (iDepth==nDepth-1) then
               if ((bShedVorticity) .and. (bShedLastVorticity)) then
                  SegConnct(1,iHeadC) = iseg4
                  SegConnct(2,iHeadC) = iseg3
                  SegConnct(3,iHeadC) = iDepth
                  SegConnct(4,iHeadC) = iSpan
                  SegGamma (iHeadC  ) = - LatticeGamma(iSpan,iDepth)
                  SegEpsilon(iHeadC ) = LatticeEpsilon(1,iSpan,iDepth) ! Using epsilon x 
                  iHeadC=iHeadC+1
               endif
            endif
            ! Segment 2-3
            if (iSpan==nSpan-1) then
               SegConnct(1,iHeadC) = iseg2
               SegConnct(2,iHeadC) = iseg3
               SegConnct(3,iHeadC) = iDepth
               SegConnct(4,iHeadC) = iSpan
               SegGamma (iHeadC  ) = LatticeGamma(iSpan,iDepth)
               SegEpsilon(iHeadC ) = LatticeEpsilon(3,iSpan,iDepth) ! Using epsilon z
               iHeadC=iHeadC+1
            endif
         enddo
      enddo
   end subroutine LatticeToSegments

    !> Convert segments between index iSegStart and iSegEnd to particles. 
    subroutine SegmentsToPart(SegPoints, SegConnct, SegGamma, SegEpsilon, iSegStart, iSegEnd, nPartPerSeg, PartPoints, PartAlpha, PartEpsilon, iHeadPart)
       real(ReKi),     dimension(:,:), intent(in   ) :: SegPoints     !< 
       integer(IntKi), dimension(:,:), intent(in   ) :: SegConnct     !< 
       real(ReKi),     dimension(:),   intent(in   ) :: SegGamma      !< 
       real(ReKi),     dimension(:),   intent(in   ) :: SegEpsilon    !< 
       integer,                        intent(in   ) :: iSegStart     !< Index where to start in Seg* vectors
       integer,                        intent(in   ) :: iSegEnd       !<
       integer,                        intent(in   ) :: nPartPerSeg   !< Segments will be dividied into nPartPerSeg particles
       real(ReKi),     dimension(:,:), intent(inout) :: PartPoints    !< Particle points (3 x nPart) 
       real(ReKi),     dimension(:,:), intent(inout) :: PartAlpha     !< Particle intensities (3 x nPart)
       real(ReKi),     dimension(:),   intent(inout) :: PartEpsilon   !< Particle regularization parameter (nPart)
       integer,    optional,           intent(inout) :: iHeadPart     !< Index where to start in Part* vectors
       real(ReKi), dimension(3) :: P1, P2, DP  !< Segment extremities
       real(ReKi), dimension(3) :: SegDir  !< direction vector
       real(ReKi), dimension(3) :: PartInt
       real(ReKi)     :: PartLen         !< Initial "length" of the blob
       real(ReKi)     :: PartEps         !< Regularization of the blob
       real(ReKi)     :: SegLen
       integer(IntKi) :: iPart           !< index in particle vectors
       integer(IntKi) :: iSeg, iSubPart
       if (present(iHeadPart)) then
          iPart = iHeadPart
       else
          iPart = 1
       endif
       ! loop on selected segments
       do iSeg=iSegStart,iSegEnd
          P1 = SegPoints(1:3,SegConnct(1,iSeg)) ! Segment extremities
          P2 = SegPoints(1:3,SegConnct(2,iSeg))
          DP = P2-P1
          SegLen  = sqrt(DP(1)**2 + DP(2)**2 + DP(3)**2)
          if (SegLen>0) then
             SegDir  = DP/SegLen                     ! Unit vector along segment direction
          else
             ! For now we set direction to zero. Part. Intensity will be zero. Multiple zero part may be created at P1
             ! In the future we should skip the creation of the particles all together
             SegDir  = 0
             if (DEV_VERSION_VT) then
                ! This is prone to happening with turbulence and single precision
                print*,'OLAF: encountered a segment of zero length'
                STOP
             endif
          endif
          PartInt = DP*SegGamma(iSeg)/nPartPerSeg ! alpha = Gamma.L/n = omega.dV [m^3/s]
          PartEps = SegEpsilon(iSeg)                   ! TODO this might need tuning depending on RegFunction and n_new
          PartLen = SegLen/nPartPerSeg
          do iSubPart=0,nPartPerSeg-1
             PartPoints(1:3, iPart) = P1(1:3) + (0.5_ReKi+iSubPart)*PartLen*SegDir(1:3) ! ds/2:ds:L 
             PartAlpha (1:3, iPart) = PartInt(1:3)
             PartEpsilon(    iPart) = PartEps
             iPart = iPart +1 
          enddo
       enddo 
       if (present(iHeadPart)) then
          iHeadPart=iPart
       endif
    end subroutine SegmentsToPart

   subroutine print_mean_4d(M, Label)
      real(ReKi), dimension(:,:,:,:), intent(in) :: M
      character(len=*), intent(in)               :: Label
      integer(IntKi) :: i, j, k
      real(ReKi), dimension(3) :: U
      !
      U(1:3)=0
      if ((size(M,4)*size(M,3)*size(M,2) )>0) then
         do i=1,size(M,4); do j=1,size(M,3); do k=1,size(M,2); 
            U(1:3)= U(1:3)+ M(1:3, k, j, i)
         enddo; enddo; enddo; 
         U(1:3)=U(1:3)/ (size(M,4)*size(M,3)*size(M,2))
      endif
      print'(A25,3F12.4)',trim(Label),U
   end subroutine

   subroutine print_mean_3d(M, Label)
      real(ReKi), dimension(:,:,:), intent(in) :: M
      character(len=*), intent(in)               :: Label
      integer(IntKi) :: i, j
      real(ReKi), dimension(3) :: U
      !
      U(1:3)=0
      if ((size(M,3)*size(M,2))>0) then
         do i=1,size(M,3); do j=1,size(M,2)
            U(1:3)= U(1:3)+ M(1:3, j, i)
         enddo; enddo;
         U(1:3)=U(1:3)/ (size(M,3)*size(M,2))
      endif
      print'(A25,3F12.4)',trim(Label),U
   end subroutine

   subroutine find_equal_points(array, varname)
      real(ReKi), dimension(:,:,:), intent(in) :: array ! 3*n*m
      character(len=*), intent(in) :: varname
      logical :: found
      integer :: iSpan, jSpan, iDepth, jDepth, tot
      real(ReKi), dimension(3) :: DP, P1, P2
      real(ReKi) ::s1
      ! Check if two points are the same
      found=.false.
      tot=0
      do iDepth = 1, size(array,3)
         do iSpan = 1, size(array,2)
            P1 = array(1:3, iSpan  , iDepth)
            do jDepth = iDepth, size(array,3)
               do jSpan = iSpan, size(array,2)
                  if ((iDepth==jDepth) .and. (iSpan == jSpan)) then
                     ! pass
                  else
                     P2 = array(1:3, jSpan, jDepth)
                     DP = P2-P1; s1  = sqrt(DP(1)**2 + DP(2)**2 + DP(3)**2)
                     if (s1<=0) then
                        tot=tot+1; found=.true.
                        if (tot<5) then
                           print*,'Two points are the same on the lattice!'
                           print*,'P1', P1, 'i', iSpan, iDepth
                           print*,'P2', P2, 'j', jSpan, jDepth
                        endif
                     endif

                  endif
               enddo
            enddo
         enddo
      enddo
      if (found) then
         print*,'>>>>>>>>>>>>> Points the same ',trim(varname), tot, size(array,2)*size(array,3)
         print*,'>>>>>>>>>>>>> Lattice size ', size(array,2), size(array,3)
         STOP
      endif
   end subroutine find_equal_points

   !> Perform interpolation from control points to nodes assuming CP are between nodes
   subroutine interpextrap_cp2node(xin, yin, xnew, ynew)
      real(ReKi), intent(in   ) :: xin(:)
      real(ReKi), intent(in   ) :: yin(:)
      real(ReKi), intent(in   ) :: xnew(:)
      real(ReKi), intent(  out) :: ynew(:)
      integer(IntKi) :: n
      n=size(xin)
      call InterpArray(xin, yin, xnew(2:n), ynew(2:n))
      ! Boundaries
      if (n>1) then ! If more than 2 panels, use extrapolation
         ynew(1)   = lin_extrap(xnew(1)  , xin(1), yin(1), xin(2)  , yin(2))
         ynew(n+1) = lin_extrap(xnew(n+1), xin(n), yin(n), xin(n-1), yin(n-1))
      else ! If one panel, duplicate the unique point on both side
         ynew(1)   = yin(1)
         ynew(n+1) = yin(n) !n=1
      endif
   contains
      !> Perform linear extrapolation to get value of y(x0), using y(x1) and y(x2)
      real(ReKi) function lin_extrap(x0, x1, y1, x2, y2) result(y0)
         real(ReKi), intent(in)   :: x0, x1, y1, x2, y2
         real(ReKi) :: a
         a = (x0-x1)/(x0-x2)
         y0 = 1._ReKi/(1._ReKi-a) * (y1-a*y2)
      end function lin_extrap
   end subroutine interpextrap_cp2node

   ! --------------------------------------------------------------------------------}
   ! --- Tree -Grow 
   ! --------------------------------------------------------------------------------{
   subroutine grow_tree_part(Tree, nPart, PartP, PartAlpha, PartRegFunction, PartRegParam, iStep)
      type(T_Tree),               intent(inout), target :: Tree            !< 
      integer(IntKi),             intent(in   )         :: nPart           !< 
      real(ReKi), dimension(:,:), intent(in   ), target :: PartP           !< 
      real(ReKi), dimension(:,:), intent(in   ), target :: PartAlpha       !< 
      integer(IntKi),             intent(in   )         :: PartRegFunction !< 
      real(ReKi), dimension(:),   intent(in   ), target :: PartRegParam    !< 
      integer(IntKi),             intent(in   )         :: iStep           !< 
      type(T_Node), pointer :: node !< Alias
      type(T_VPart), pointer :: Part !< Alias
      real(ReKi) :: max_x,max_y,max_z !< for domain dimension
      real(ReKi) :: min_x,min_y,min_z !< for domain dimension
      integer(IntKi) :: i

      ! Cutting tree if it already has content 
      if (associated(Tree%root%branches).or.associated(Tree%root%leaves)) then
         call cut_tree_parallel(Tree)
      endif
      ! Linking tree particles to given part, no copy!
      nullify(Tree%Part%P)
      nullify(Tree%Part%Alpha)
      nullify(Tree%Part%RegParam)
      Tree%Part%P           => PartP(:,1:nPart)
      Tree%Part%Alpha       => PartAlpha(:,1:nPart)
      Tree%Part%RegParam    => PartRegParam(1:nPart)
      Tree%Part%RegFunction = PartRegFunction
      Tree%Part%n           = nPart

      ! --- Handle special case for root node
      node => Tree%Root
      Part => Tree%Part
      node%nPart = Part%n
      if (Part%n==0) then
         ! Do nothing
         node%radius = -9999.99_ReKi
         node%center = -9999.99_ReKi
         node%Moments= -9999.99_ReKi
      else if (Tree%Part%n==1) then
         node%radius=0
         node%center(1:3)=Part%P(1:3,1)
         node%Moments=0.0_ReKi
         nullify(node%iPart)
         nullify(node%branches)
         allocate(node%leaves(1:1))
         node%leaves(1) = 1 !< index
      else
         ! Domain dimensions
         max_x=maxval(Part%P(1,1:Part%n)); max_y=maxval(Part%P(2,1:Part%n)); max_z=maxval(Part%P(3,1:Part%n))
         min_x=minval(Part%P(1,1:Part%n)); min_y=minval(Part%P(2,1:Part%n)); min_z=minval(Part%P(3,1:Part%n))

         ! Init of trunc
         ! Radius taken slightly bigger than domain extent. This radius will be divided by 2 successively
         node%radius = max(abs(max_x-min_x),abs(max_y-min_y),abs(max_z-min_z))*1.001_ReKi
         node%center = (/ (max_x+min_x)/2._ReKi, (max_y+min_y)/2._ReKi, (max_z+min_z)/2._ReKi /)
         node%Moments=0.0_ReKi
         if (DEV_VERSION_VT) then
            if(node%radius>1e6) then
                print*,'[Error] Domain extent too large, particle points must be invalid';
                print*, min_x, max_x, min_y, max_y, min_z, max_z
                STOP
            endif
            if(associated(node%iPart)) then ; print*,'[Error] Node part allocated'; STOP; endif
         endif
         allocate(node%iPart(1:Part%n))
         do i=1,Part%n
            node%iPart(i) = i
         end do
         if (DEV_VERSION_VT) then
            if(associated(node%branches)) then;  print*,'node branches allocated'; STOP; endif
            if(associated(node%leaves)) then;  print*,'node leaves allocated'; STOP; endif
         endif
         node%branches=>null()
         node%leaves=>null()
         ! --- Calling grow function on subbranches
         call grow_tree_part_parallel(Tree%root, Tree%Part)
!          call grow_tree_rec(Tree%root, Tree%Part)
      endif
      Tree%iStep  = iStep
      Tree%bGrown = .true.
      Tree%DistanceDirect = 2*sum(PartRegParam)/size(PartRegParam) ! 2*mean(eps), below that distance eps has a strong effect
   end subroutine grow_tree_part

   !> Recursive function to grow/setup a tree. 
   !! Note, needed preliminary calc are done by grow_tree before
   recursive subroutine grow_tree_part_rec(node, Part)
      type(T_Node), target     :: node !<
      type(T_VPart), intent(in) :: Part !<
      integer :: i
      !  Sub Step:
      !   -  compute moments and center for the current node
      !   -  allocate branches and leaves
      call grow_tree_part_substep(node, Part)
      ! Call grow_tree on branches
      if(associated(node%branches)) then
         do i = 1,size(node%branches)
            call grow_tree_part_rec(node%branches(i), Part)
         end do
      endif
   end subroutine  grow_tree_part_rec

   !> Perform a substep of tree growth, growing sub branches from a given node/cell
   !! Parent has already setup node%iPart, indices of the particle in this cell
   !! Steps are:
   !!   - Compute node center (barycenter of vorticity)
   !!   - Compute node moments
   !!   - Distribute particles in each 8 octants. Branches are not created for empty octant
   !!   - Allocate branches and leaves and distribute particles to them
   subroutine grow_tree_part_substep(node, Part)
      type(T_Node), intent(inout) :: node !< Current node we are growing from
      type(T_VPart), intent(in)    :: Part !< All particles info
      integer(IK1) :: iPartOctant                                     !< Index corresponding to which octant the particle falls into
      integer      :: nLeaves, nBranches
      integer      :: iLeaf, iOctant, iBranch
      integer      :: nPerBranchAcc(8) !< Accumulated counter on number of particles per branch
      integer      :: i !,j,k
      real(ReKi)   :: wTot, wLoc ! Total and local vorticity strength
      real(ReKi)   :: halfSize ! TODO remove me
      real(ReKi),dimension(3) :: locCenter, DeltaP,PartPos,PartAlpha
      real(ReKi),dimension(3) :: GeomC !< Geometric center from division of the domain in powers of 2
      real(ReKi),dimension(3) :: VortC !< Vorticity weighted center
      integer(IK1),dimension(:),allocatable :: PartOctant !< Stores the octant (1-8) where each particle belongs
      integer,dimension(8) :: npart_per_octant !< Number of particle per octant
      integer,dimension(8) :: octant2branches  !< Mapping between 8 octants, to index of non empty branch
      integer,dimension(8) :: octant2leaves    !< Idem for singleton/leaves 
      real(ReKi) :: max_x,max_y,max_z !< for domain dimension
      real(ReKi) :: min_x,min_y,min_z !< for domain dimension
      GeomC = node%center ! NOTE: we rely on the fact that our parent has set this to the Geometric value 
      VortC = 0.0_ReKi
      wTot = 0.0_ReKi
      ! --- Barycenter of vorticity of the node 
      do i = 1,node%nPart
         PartPos   = Part%P(:,node%iPart(i))
         PartAlpha = Part%Alpha(:,node%iPart(i))
         wLoc      = (PartAlpha(1)**2 + PartAlpha(2)**2 + PartAlpha(3)**2)**0.5_ReKi ! Vorticity norm
         VortC     = VortC + wLoc*PartPos                                   ! Sum coordinates weighted by vorticity
         wTot      = wTot + wLoc                                                     ! Total vorticity
      end do
      ! There is no vorticity, we make it a empty node and we exit 
      if(EqualRealNos(abs(wTot),0.0_ReKi)) then
         node%nPart=0
         if (associated(node%iPart)) deallocate(node%iPart)
         return ! NOTE: we exit 
      endif
      VortC = VortC/wTot ! barycenter of vorticity
      node%center  = VortC ! updating center, it's now the vorticity center, geomCenter will be use to divide cell though

      ! --- Calculation of moments about VortC
      do i = 1,node%nPart
         PartPos   = Part%P    (:,node%iPart(i))
         PartAlpha = Part%Alpha(:,node%iPart(i))
         DeltaP    = PartPos-VortC
         ! Order 0
         node%Moments(1:3,M0_000) = node%Moments(1:3,M0_000) + PartAlpha
         ! 1st order
         node%Moments(1:3,M1_100) = node%Moments(1:3,M1_100) + PartAlpha*DeltaP(1) ! 100
         node%Moments(1:3,M1_010) = node%Moments(1:3,M1_010) + PartAlpha*DeltaP(2) ! 010
         node%Moments(1:3,M1_001) = node%Moments(1:3,M1_001) + PartAlpha*DeltaP(3) ! 001
         ! 2nd order
         ! KEEP ME:
         !   do j=1,3
         !      do k=1,j
         !         node%Moments(1:3,3+j+k+j/3) = node%Moments(1:3,3+j+k+j/3) + PartAlpha*DeltaP(j)*DeltaP(k)
         !      end do
         !   end do
         node%Moments(1:3, 5) = node%Moments(1:3, 5) + PartAlpha*DeltaP(1)*DeltaP(1) ! j=1,k=1
         node%Moments(1:3, 6) = node%Moments(1:3, 6) + PartAlpha*DeltaP(2)*DeltaP(1) ! j=2,k=1
         node%Moments(1:3, 7) = node%Moments(1:3, 7) + PartAlpha*DeltaP(2)*DeltaP(2) ! j=2,k=2
         node%Moments(1:3, 8) = node%Moments(1:3, 8) + PartAlpha*DeltaP(3)*DeltaP(1) ! j=3,k=1
         node%Moments(1:3, 9) = node%Moments(1:3, 9) + PartAlpha*DeltaP(3)*DeltaP(2) ! j=3,k=2
         node%Moments(1:3,10) = node%Moments(1:3,10) + PartAlpha*DeltaP(3)*DeltaP(3) ! j=3,k=3
      end do

      ! --- Distributing particles to the 8 octants (based on the geometric center!)
      allocate (PartOctant(1:node%nPart))
      npart_per_octant(1:8)=0
      do i = 1,node%nPart
         PartPos      = Part%P(:,node%iPart(i))
         ! index corresponding to which octant the particle falls into
         iPartOctant = int(1,IK1)
         if (PartPos(1) > GeomC(1)) iPartOctant = iPartOctant + int(1,IK1)
         if (PartPos(2) > GeomC(2)) iPartOctant = iPartOctant + int(2,IK1)
         if (PartPos(3) > GeomC(3)) iPartOctant = iPartOctant + int(4,IK1)
         npart_per_octant(iPartOctant) = npart_per_octant(iPartOctant) + 1 ! Counter of particles per octant
         PartOctant(i)=iPartOctant ! Store in which octant particle i is
      end do

      ! --- Leaves and branches
      ! A node contains a combination of child nodes and leaves (single particles)
      ! TODO: introduce a "minimum cell size", (e.g. cell radius is less than the Distance for direct evaluation, then all should be leaves)
      nLeaves          = 0
      nBranches        = 0
      octant2branches  = 0
      octant2leaves    = 0
      do iOctant = 1,8
         if(npart_per_octant(iOctant)==1) then 
            nLeaves                = nLeaves+1
            octant2leaves(iOctant) = nLeaves
         else if(npart_per_octant(iOctant)>1) then
            if (npart_per_octant(iOctant)==node%nPart) then
               ! All particle falls into the same octant, if they all have the same location, we would divide forever.
               ! Quick fix below
               max_x=maxval(Part%P(1,node%iPart(:))); max_y=maxval(Part%P(2,node%iPart(:))); max_z=maxval(Part%P(3,node%iPart(:)))
               min_x=minval(Part%P(1,node%iPart(:))); min_y=minval(Part%P(2,node%iPart(:))); min_z=minval(Part%P(3,node%iPart(:)))
               if (max(abs(max_x-min_x),abs(max_y-min_y),abs(max_z-min_z))< 1.0e-5) then
                  nLeaves=node%nPart
                  allocate (node%leaves(1:nLeaves))
                  do i = 1,node%nPart
                     node%leaves(i)=node%iPart(i)
                  enddo
                  ! Cleanup and exit!
                  if (associated(node%iPart)) deallocate(node%iPart) ! Freeing memory
                  if (allocated(PartOctant)) deallocate(PartOctant)
                  return
               endif
            endif
            nBranches                = nBranches+1
            octant2branches(iOctant) = nBranches
         endif
      enddo
      if (DEV_VERSION_VT) then
         if (associated(node%branches)) then
            print*,'Tree build: error, branches associated'
            STOP
         endif
         if (associated(node%leaves)) then
            print*,'Tree build: error, leaves associated'
            STOP
         end if
      endif

      if(nBranches>0) allocate (node%branches(1:nBranches))
      if(nLeaves>0)   allocate (node%leaves(1:nLeaves))

      ! --- Initializing the branches nodes and leaves
      halfSize = node%radius/2._ReKi
      do iOctant = 1,8 ! there is max 8 octant
         iBranch     = octant2branches(iOctant)
         if (iBranch>0) then ! this node has branches 
            allocate(node%branches(iBranch)%iPart(1:npart_per_octant(iOctant)))
            node%branches(iBranch)%nPart=npart_per_octant(iOctant)
            ! NOTE: this is geometric center not barycenter
            locCenter = GeomC + 0.5*halfSize*(/ (-1)**(iOctant), (-1)**floor(0.5*real(iOctant-1)+1), (-1)**floor(0.25*real(iOctant-1)+1) /)
            ! Init of branches
            node%branches(iBranch)%radius  = halfSize  ! 
            node%branches(iBranch)%center  = locCenter ! NOTE: this is the geometric center
            node%branches(iBranch)%Moments = 0.0_ReKi  ! 
            node%branches(iBranch)%branches=>null()
            node%branches(iBranch)%leaves=>null()
         endif
         ! other cases are leaves or dead branches
      end do

      ! Store indices of the particles the sub-branch contains
      nPerBranchAcc(:) = 0
      do i = 1,node%nPart
         iBranch = octant2branches(PartOctant(i))
         if(iBranch>0) then
            nPerBranchAcc(iBranch)=nPerBranchAcc(iBranch)+1; node%branches(iBranch)%iPart(nPerBranchAcc(iBranch)) = node%iPart(i)
         else
            iLeaf = octant2leaves(PartOctant(i))
            if(iLeaf>0) then
               node%leaves(iLeaf)=node%iPart(i)
            else
               print*,'This particle do not belong to anybody!!',i
               STOP
            endif
         endif
      end do
      if (DEV_VERSION_VT) then
         do iBranch=1,nBranches
            if (nPerBranchAcc(iBranch)/=node%branches(iBranch)%nPart) then
               print*,'Grow tree: Problem in repartition of particles per branches'
               STOP
            endif
         enddo
      endif
      if (associated(node%iPart)) deallocate(node%iPart) ! Freeing memory
      if (allocated(PartOctant)) deallocate(PartOctant)
   end subroutine grow_tree_part_substep

   !> Grow a tree in "parallel", since recursive calls cannot be parallized easily, we unroll the different layer calls
   !! Note, needed preliminary calc are done by grow_tree before!
   subroutine grow_tree_part_parallel(Root, Part)
      type(T_Node), intent(inout) :: Root
      type(T_VPart), intent(in) :: Part
      integer :: i, nBranches
      integer :: i1
      integer :: i2

      ! ---  Unrolled version of grow_tree for the first node
      !  Sub Step:
      !   -  compute moments and center for the current node
      !   -  allocate branches and leaves
      call grow_tree_part_substep(Root, Part)

      if(.not. associated(Root%branches)) then
         nBranches=0
      else
         nBranches=size(Root%branches)
         if (nBranches==0) then
            print*,'No branches' ! This should not happen
            STOP
         else
            ! Call "grow_tree" on branches

            !$OMP PARALLEL default(shared)

            ! ---  Unrolled version of grow_tree for the second levels
            !$OMP do private(i) schedule(runtime)
            do i = 1,nBranches ! maximum 8 branches
               if(Root%branches(i)%nPart>1) then ! I dont think this test is needed
                  call grow_tree_part_substep(Root%branches(i), Part)
               endif
            end do
            !$OMP end do
            !$OMP barrier ! we need to be sure that all the branches were built
            ! ---  Unrolled version of grow_tree for third node levels
            !$OMP do private(i,i1,i2) schedule(runtime)
            do i = 1,nBranches*8 ! maximum 64 sub branches
               i1=(i-1)/8+1;
               i2=mod(i-1,8)+1;
               if(associated(Root%branches(i1)%branches)) then
                  if (i2<=size(Root%branches(i1)%branches)) then
                     call grow_tree_part_rec(Root%branches(i1)%branches(i2), Part)
                  endif
               endif
            enddo
            !$OMP end do 
            !Note: We could add more levels 
            !$OMP END PARALLEL
         endif
      endif
   end subroutine grow_tree_part_parallel
   ! --------------------------------------------------------------------------------}
   ! --- Tree -Grow for vortex lines
   ! --------------------------------------------------------------------------------{
   subroutine grow_tree_segment(Tree_Seg, nSeg, SegPoints, SegConnct, SegGamma, SegRegFunction, SegRegParam, iStep)
      type(T_Tree),               intent(inout), target :: Tree_Seg           !<
      integer(IntKi),                 intent(in   )         :: nSeg           !<
      real(ReKi), dimension(:,:),     intent(in   ), target :: SegPoints      !<
      integer(IntKi), dimension(:,:), intent(in   ), target :: SegConnct  !<
      real(ReKi), dimension(:),       intent(in   ), target :: SegGamma       !<
      integer(IntKi),                 intent(in   )         :: SegRegFunction !<
      real(ReKi), dimension(:),       intent(in   ), target :: SegRegParam    !<
      integer(IntKi),                 intent(in   )         :: iStep          !<
      type(T_Node), pointer :: node !< Alias
      type(T_VSgmt), pointer :: Seg !< Alias
      real(ReKi) :: max_x,max_y,max_z !< for domain dimension
      real(ReKi) :: min_x,min_y,min_z !< for domain dimension
      integer(IntKi) :: i

      ! Cutting tree if it already has content
      if (associated(Tree_Seg%root%branches).or.associated(Tree_Seg%root%leaves)) then
         call cut_tree_segment_parallel(Tree_Seg)
      endif
      ! Linking tree particles to given part, no copy!
      nullify(Tree_Seg%Seg%SP)
      nullify(Tree_Seg%Seg%SConnct)
      nullify(Tree_Seg%Seg%SGamma)
      nullify(Tree_Seg%Seg%RegParam)
      Tree_Seg%Seg%SP          => SegPoints
      Tree_Seg%Seg%SConnct     => SegConnct(:,1:nSeg)
      Tree_Seg%Seg%SGamma      => SegGamma(1:nSeg)
      Tree_Seg%Seg%RegParam    => SegRegParam(1:nSeg)
      Tree_Seg%Seg%RegFunction = SegRegFunction
      Tree_Seg%Seg%n           = nSeg

      ! --- Handle special case for root node
      node => Tree_Seg%Root
      Seg => Tree_Seg%Seg
      if (Seg%n==0) then
         ! Do nothing
         node%radius = -9999.99_ReKi
         node%center = -9999.99_ReKi
         node%Moments= -9999.99_ReKi
      else if (Tree_Seg%Seg%n==1) then
         node%radius=0
         node%center(1:3)=0.5_ReKi*(Seg%SP(1:3,1)+Seg%SP(1:3,2))
         node%Moments=0.0_ReKi
         nullify(node%iPart)
         nullify(node%branches)
         allocate(node%leaves(1:1))
         node%leaves(1) = Seg%n !< index
         node%nPart    = 1
      else
         ! Domain dimensions
         max_x=max(maxval(Seg%SP(1,Seg%SConnct(1,:))),maxval(Seg%SP(1,Seg%SConnct(2,:))))
         max_y=max(maxval(Seg%SP(2,Seg%SConnct(1,:))),maxval(Seg%SP(2,Seg%SConnct(2,:))))
         max_z=max(maxval(Seg%SP(3,Seg%SConnct(1,:))),maxval(Seg%SP(3,Seg%SConnct(2,:))))
         min_x=min(minval(Seg%SP(1,Seg%SConnct(1,:))),minval(Seg%SP(1,Seg%SConnct(2,:))))
         min_y=min(minval(Seg%SP(2,Seg%SConnct(1,:))),minval(Seg%SP(2,Seg%SConnct(2,:))))
         min_z=min(minval(Seg%SP(3,Seg%SConnct(1,:))),minval(Seg%SP(3,Seg%SConnct(2,:))))

         ! Init of trunc
         ! Radius taken slightly bigger than domain extent. This radius will be divided by 2 successively
         node%radius = max(abs(max_x-min_x),abs(max_y-min_y),abs(max_z-min_z))*1.001_ReKi
         node%center = (/ (max_x+min_x)/2._ReKi, (max_y+min_y)/2._ReKi, (max_z+min_z)/2._ReKi /)
         node%Moments=0.0_ReKi
         if (DEV_VERSION_VT) then
            if (any(Seg%SP(1,:)<-999.99_ReKi)) then
              print*,'Error in segment transmission to grow tree segment'
              STOP
            endif
            if(node%radius>1e6) then
                print*,'[Error] Domain extent too large, segment points must be invalid';
                print*, min_x, max_x, min_y, max_y, min_z, max_z
                STOP
            endif
            if(associated(node%iPart)) then ; print*,'[Error] Node part allocated'; STOP; endif
         endif
         allocate(node%iPart(1:Seg%n))
         do i=1,Seg%n
            node%iPart(i) = i
         end do
         if (DEV_VERSION_VT) then
            if(associated(node%branches)) then;  print*,'node branches allocated'; STOP; endif
            if(associated(node%leaves)) then;  print*,'node leaves allocated'; STOP; endif
         endif
         node%branches=>null()
         node%leaves=>null()
         node%nPart=Seg%n
         ! --- Calling grow function on subbranches
         call grow_tree_segment_parallel(Tree_Seg%root, Tree_Seg%Seg)
         !call grow_tree_segment_rec(Tree_Seg%root, Tree_Seg%Seg)
      endif
      Tree_Seg%iStep  = iStep
      Tree_Seg%bGrown = .true.
      Tree_Seg%DistanceDirect = 2*sum(SegRegParam)/size(SegRegParam) ! 2*mean(eps), below that distance eps has a strong effect ! TODO REMOVE
   end subroutine grow_tree_segment

   !> Recursive function to grow/setup a tree.
   !! Note, needed preliminary calc are done by grow_tree before
   recursive subroutine grow_tree_segment_rec(node, Seg)
      type(T_Node), target     :: node !<
      type(T_VSgmt), intent(in) :: Seg !<
      integer :: i
      !  Sub Step:
      !   -  compute moments and center for the current node
      !   -  allocate branches and leaves
      call grow_tree_segment_substep(node, Seg)
      ! Call grow_tree on branches
      if(associated(node%branches)) then
         do i = 1,size(node%branches)
            call grow_tree_segment_rec(node%branches(i), Seg)
         end do
      endif
   end subroutine  grow_tree_segment_rec

   !> Perform a substep of tree growth, growing sub branches from a given node/cell
   !! Parent has already setup node%iPart, indices of the particle in this cell
   !! Steps are:
   !!   - Compute node center (barycenter of vorticity)
   !!   - Compute node moments
   !!   - Distribute particles in each 8 octants. Branches are not created for empty octant
   !!   - Allocate branches and leaves and distribute particles to them
   subroutine grow_tree_segment_substep(node, Seg)
      type(T_Node), intent(inout) :: node !< Current node we are growing from
      type(T_VSgmt), intent(in)    :: Seg !< All particles info
      integer(IK1) :: iPartOctant                                     !< Index corresponding to which octant the particle falls into
      integer      :: nLeaves, nBranches
      integer      :: iLeaf, iOctant, iBranch
      integer      :: nPerBranchAcc(8) !< Accumulated counter on number of particles per branch
      integer      :: i!,j,k
      real(ReKi)   :: wTot ! Total vorticity strength
      real(ReKi)   :: halfSize ! TODO remove me
      real(ReKi),dimension(3) :: locCenter, SegCenter,DP, SegGammaVec
      real(ReKi),dimension(3) :: P1,P2 !< Segment extremities
      real(ReKi),dimension(3) :: GeomC !< Geometric center from division of the domain in powers of 2
      real(ReKi),dimension(3) :: VortC !< Vorticity weighted center
      integer(IK1),dimension(:),allocatable :: PartOctant !< Stores the octant (1-8) where each particle belongs
      integer,dimension(8) :: npart_per_octant !< Number of particle per octant
      integer,dimension(8) :: octant2branches  !< Mapping between 8 octants, to index of non empty branch
      integer,dimension(8) :: octant2leaves    !< Idem for singleton/leaves
      real(ReKi) :: max_x,max_y,max_z !< for domain dimension
      real(ReKi) :: min_x,min_y,min_z !< for domain dimension
      GeomC = node%center ! NOTE: we rely on the fact that our parent has set this to the Geometric value
      VortC = 0.0_ReKi
      wTot = 0.0_ReKi
      ! --- Barycenter of vorticity of the node
      do i = 1,node%nPart
         P1 = Seg%SP(1:3,Seg%SConnct(1,node%iPart(i)))
         P2 = Seg%SP(1:3,Seg%SConnct(2,node%iPart(i)))
         SegCenter = 0.5_ReKi*(P1+P2)
         VortC = VortC + abs(Seg%SGamma(node%iPart(i)))*SegCenter  ! Sum coordinates weighted by vorticity
         wTot   = wTot + abs(Seg%SGamma(node%iPart(i)))            ! Total vorticity
      end do
      ! There is no vorticity, we make it a empty node and we exit
      if(EqualRealNos(abs(wTot),0.0_ReKi)) then
         node%nPart=0
         if (associated(node%iPart)) deallocate(node%iPart)
         return ! NOTE: we exit
      endif
      VortC = VortC/wTot ! barycenter of vorticity
      node%center  = VortC ! updating center, it's now the vorticity center, geomCenter will be use to divide cell though

      ! --- Calculation of moments about VortC
      do i = 1,node%nPart
         P1 = Seg%SP(1:3,Seg%SConnct(1,node%iPart(i)))
         P2 = Seg%SP(1:3,Seg%SConnct(2,node%iPart(i)))
         DP = P2-P1
         SegGammaVec = Seg%SGamma(node%iPart(i))*DP            !Vorticity vector
         ! Order 0
         node%Moments(1:3,M0_000) = node%Moments(1:3,M0_000) + SegGammaVec
         node%Moments(1:3,M1_100) = node%Moments(1:3,M1_100) + SegGammaVec*(0.5_ReKi*(P1(1)+P2(1))-VortC(1)) ! 100
         node%Moments(1:3,M1_010) = node%Moments(1:3,M1_010) + SegGammaVec*(0.5_ReKi*(P1(2)+P2(2))-VortC(2)) ! 010
         node%Moments(1:3,M1_001) = node%Moments(1:3,M1_001) + SegGammaVec*(0.5_ReKi*(P1(3)+P2(3))-VortC(3)) ! 001
         ! 2nd order
         ! KEEP ME
         !do j=1,3
         !   do k=1,j
         !      if (j==k) then
         !        node%Moments(1:3,3+j+k+j/3) = node%Moments(1:3,3+j+k+j/3) + SegGammaVec*1/3.0_ReKi*(3*VortC(j)**2 -3*VortC(j)*(P1(j)+P2(j))+P1(j)**2+P1(j)*P2(j)+P2(j)**2)

         !      else
         !        node%Moments(1:3,3+j+k+j/3) = node%Moments(1:3,3+j+k+j/3) + SegGammaVec*1/6.0_ReKi*( 6*VortC(j)*VortC(k) -3*VortC(j)*(P1(k)+P2(k)) +P1(j)*(-3*VortC(k)+2*P1(k)+P2(k)) +P2(j)*(-3*VortC(k)+P1(k)+2*P2(k)))
         !      end if
         !   end do
         !end do

         !node%Moments(1:3,3+j+k+j/3) = node%Moments(1:3,3+j+k+j/3) + SegGammaVec*1/3.0_ReKi*(3*VortC(j)**2 -3*VortC(j)*(P1(j)+P2(j))+P1(j)**2+P1(j)*P2(j)+P2(j)**2)
         node%Moments(1:3,   5     ) = node%Moments(1:3,   5     ) + SegGammaVec*1/3.0_ReKi*(3*VortC(1)**2 -3*VortC(1)*(P1(1)+P2(1))+P1(1)**2+P1(1)*P2(1)+P2(1)**2) !j=1,k=1
         node%Moments(1:3,   7     ) = node%Moments(1:3,   7     ) + SegGammaVec*1/3.0_ReKi*(3*VortC(2)**2 -3*VortC(2)*(P1(2)+P2(2))+P1(2)**2+P1(2)*P2(2)+P2(2)**2) !j=2,k=2
         node%Moments(1:3,   10    ) = node%Moments(1:3,   10    ) + SegGammaVec*1/3.0_ReKi*(3*VortC(3)**2 -3*VortC(3)*(P1(3)+P2(3))+P1(3)**2+P1(3)*P2(3)+P2(3)**2) !j=3,k=3
         !node%Moments(1:3,3+j+k+j/3) = node%Moments(1:3,3+j+k+j/3) + SegGammaVec*1/6.0_ReKi*( 6*VortC(j)*VortC(k)  -3*VortC(j)*(P1(k)+P2(k))  +P1(j)*(-3*VortC(k)+2*P1(k)+P2(k))  +P2(j)*(-3*VortC(k)+P1(k)+2*P2(k)))
         node%Moments(1:3,   6     ) = node%Moments(1:3,   6     ) + SegGammaVec*1/6.0_ReKi*( 6*VortC(2)*VortC(1)  -3*VortC(2)*(P1(1)+P2(1))  +P1(2)*(-3*VortC(1)+2*P1(1)+P2(1))  +P2(2)*(-3*VortC(1)+P1(1)+2*P2(1))) !j=2, k=1
         node%Moments(1:3,   8     ) = node%Moments(1:3,   8     ) + SegGammaVec*1/6.0_ReKi*( 6*VortC(3)*VortC(1)  -3*VortC(3)*(P1(1)+P2(1))  +P1(3)*(-3*VortC(1)+2*P1(1)+P2(1))  +P2(3)*(-3*VortC(1)+P1(1)+2*P2(1))) !j=3, k=1
         node%Moments(1:3,   9     ) = node%Moments(1:3,   9     ) + SegGammaVec*1/6.0_ReKi*( 6*VortC(3)*VortC(2)  -3*VortC(3)*(P1(2)+P2(2))  +P1(3)*(-3*VortC(2)+2*P1(2)+P2(2))  +P2(3)*(-3*VortC(2)+P1(2)+2*P2(2))) !j=3, k=2
      end do

      ! --- Distributing particles to the 8 octants (based on the geometric center!)
      allocate (PartOctant(1:node%nPart))
      npart_per_octant(1:8)=0
      do i = 1,node%nPart
         P1 = Seg%SP(1:3,Seg%SConnct(1,node%iPart(i)))
         P2 = Seg%SP(1:3,Seg%SConnct(2,node%iPart(i)))
         SegCenter = 0.5_ReKi*(P1+P2) ! We use the segment center. We could consider looking at point 1 and 2 separately 
         ! index corresponding to which octant the particle falls into
         iPartOctant = int(1,IK1)
         if (SegCenter(1) > GeomC(1)) iPartOctant = iPartOctant + int(1,IK1)
         if (SegCenter(2) > GeomC(2)) iPartOctant = iPartOctant + int(2,IK1)
         if (SegCenter(3) > GeomC(3)) iPartOctant = iPartOctant + int(4,IK1)
         npart_per_octant(iPartOctant) = npart_per_octant(iPartOctant) + 1 ! Counter of particles per octant
         PartOctant(i)=iPartOctant ! Store in which octant particle i is
      end do

      ! --- Leaves and branches
      ! A node contains a combination of child nodes and leaves (single particles)
      ! TODO: introduce a "minimum cell size", (e.g. cell radius is less than the Distance for direct evaluation, then all should be leaves)
      nLeaves          = 0
      nBranches        = 0
      octant2branches  = 0
      octant2leaves    = 0
      do iOctant = 1,8
         if(npart_per_octant(iOctant)==1) then
            nLeaves                = nLeaves+1
            octant2leaves(iOctant) = nLeaves
         else if(npart_per_octant(iOctant)>1) then
            if (npart_per_octant(iOctant)==node%nPart) then
               ! All particle falls into the same octant, if they all have the same location, we would divide forever.
               ! Quick fix below
               max_x=max(maxval(Seg%SP(1,Seg%SConnct(1,node%iPart(1:(node%nPart))))),maxval(Seg%SP(1,Seg%SConnct(2,node%iPart(1:(node%nPart))))))
               max_y=max(maxval(Seg%SP(2,Seg%SConnct(1,node%iPart(1:(node%nPart))))),maxval(Seg%SP(2,Seg%SConnct(2,node%iPart(1:(node%nPart))))))
               max_z=max(maxval(Seg%SP(3,Seg%SConnct(1,node%iPart(1:(node%nPart))))),maxval(Seg%SP(3,Seg%SConnct(2,node%iPart(1:(node%nPart))))))
               min_x=min(minval(Seg%SP(1,Seg%SConnct(1,node%iPart(1:(node%nPart))))),minval(Seg%SP(1,Seg%SConnct(2,node%iPart(1:(node%nPart))))))
               min_y=min(minval(Seg%SP(2,Seg%SConnct(1,node%iPart(1:(node%nPart))))),minval(Seg%SP(2,Seg%SConnct(2,node%iPart(1:(node%nPart))))))
               min_z=min(minval(Seg%SP(3,Seg%SConnct(1,node%iPart(1:(node%nPart))))),minval(Seg%SP(3,Seg%SConnct(2,node%iPart(1:(node%nPart))))))
               if (max(abs(max_x-min_x),abs(max_y-min_y),abs(max_z-min_z))< 1.0e-5) then
                  nLeaves=node%nPart
                  allocate (node%leaves(1:nLeaves))
                  do i = 1,node%nPart
                     node%leaves(i)=node%iPart(i)
                  enddo
                  ! Cleanup and exit!
                  if (associated(node%iPart)) deallocate(node%iPart) ! Freeing memory
                  if (allocated(PartOctant)) deallocate(PartOctant)
                  return
               endif
            endif
            nBranches                = nBranches+1
            octant2branches(iOctant) = nBranches
         endif
      enddo
      if (DEV_VERSION_VT) then
         if (associated(node%branches)) then
            print*,'Tree build: error, branches associated'
            STOP
         endif
         if (associated(node%leaves)) then
            print*,'Tree build: error, leaves associated'
            STOP
         end if
      endif

      if(nBranches>0) allocate (node%branches(1:nBranches))
      if(nLeaves>0)   allocate (node%leaves(1:nLeaves))

      ! --- Initializing the branches nodes and leaves
      halfSize = node%radius/2._ReKi
      do iOctant = 1,8 ! there is max 8 octant
         iBranch     = octant2branches(iOctant)
         if (iBranch>0) then ! this node has branches
            allocate(node%branches(iBranch)%iPart(1:npart_per_octant(iOctant)))
            node%branches(iBranch)%nPart=npart_per_octant(iOctant)
            ! NOTE: this is geometric center not barycenter
            locCenter = GeomC + 0.5*halfSize*(/ (-1)**(iOctant), (-1)**floor(0.5*real(iOctant-1)+1), (-1)**floor(0.25*real(iOctant-1)+1) /)
            ! Init of branches
            node%branches(iBranch)%radius  = halfSize  !
            node%branches(iBranch)%center  = locCenter ! NOTE: this is the geometric center
            node%branches(iBranch)%Moments = 0.0_ReKi  !
            node%branches(iBranch)%branches=>null()
            node%branches(iBranch)%leaves=>null()
         endif
         ! other cases are leaves or dead branches
      end do

      ! Store indices of the particles the sub-branch contains
      nPerBranchAcc(:) = 0
      do i = 1,node%nPart
         iBranch = octant2branches(PartOctant(i))
         if(iBranch>0) then
            nPerBranchAcc(iBranch)=nPerBranchAcc(iBranch)+1; node%branches(iBranch)%iPart(nPerBranchAcc(iBranch)) = node%iPart(i)
         else
            iLeaf = octant2leaves(PartOctant(i))
            if(iLeaf>0) then
               node%leaves(iLeaf)=node%iPart(i)
            else
               print*,'This particle do not belong to anybody!!',i
               STOP
            endif
         endif
      end do
      if (associated(node%iPart)) deallocate(node%iPart) ! Freeing memory
      if (allocated(PartOctant)) deallocate(PartOctant)
   end subroutine grow_tree_segment_substep

   !> Grow a tree in "parallel", since recursive calls cannot be parallized easily, we unroll the different layer calls
   !! Note, needed preliminary calc are done by grow_tree before!
   subroutine grow_tree_segment_parallel(Root, Seg)
      type(T_Node), intent(inout) :: Root
      type(T_VSgmt), intent(in) :: Seg
      integer :: i, nBranches
      integer :: i1
      integer :: i2
      !write(*,*) 'Entering grow_tree_segment parallel', Seg%n
      ! ---  Unrolled version of grow_tree for the first node
      !  Sub Step:
      !   -  compute moments and center for the current node
      !   -  allocate branches and leaves
      call grow_tree_segment_substep(Root, Seg)

      if(.not. associated(Root%branches)) then
         nBranches=0
      else
         nBranches=size(Root%branches)
         if (nBranches==0) then
            print*,'No branches' ! This should not happen
            STOP
         else
            ! Call "grow_tree" on branches

            !$OMP PARALLEL default(shared)

            ! ---  Unrolled version of grow_tree for the second levels
            !$OMP do private(i) schedule(runtime)
            do i = 1,nBranches ! maximum 8 branches
               if(Root%branches(i)%nPart>1) then ! I dont think this test is needed
                  call grow_tree_segment_substep(Root%branches(i), Seg)
               endif
            end do
            !$OMP end do
            !$OMP barrier ! we need to be sure that all the branches were built
            ! ---  Unrolled version of grow_tree for third node levels
            !$OMP do private(i,i1,i2) schedule(runtime)
            do i = 1,nBranches*8 ! maximum 64 sub branches
               i1=(i-1)/8+1;
               i2=mod(i-1,8)+1;
               if(associated(Root%branches(i1)%branches)) then
                  if (i2<=size(Root%branches(i1)%branches)) then
                     call grow_tree_segment_rec(Root%branches(i1)%branches(i2), Seg)
                  endif
               endif
            enddo
            !$OMP end do
            !Note: We could add more levels
            !$OMP END PARALLEL
         endif
      endif
   end subroutine grow_tree_segment_parallel



   ! --------------------------------------------------------------------------------
   ! --- Cut tree 
   ! --------------------------------------------------------------------------------
   !> Cut a tree and all its subbranches in a recursive manner
   recursive subroutine cut_tree_rec(node)
      integer :: i
      type(T_Node),intent(inout) :: node
      call cut_substep(node)
      if (associated(node%branches)) then
         do i=1,size(node%branches)
            call cut_tree_rec(node%branches(i))
         end do
         deallocate(node%branches)
         node%branches=> null()
      end if
   end subroutine cut_tree_rec

   !> Perform a substep of tree cutting (used by recursive and parallel calls)
   subroutine cut_substep(node)
      type(T_Node), intent(inout) :: node
      if (associated(node%leaves)) then
         deallocate(node%leaves)
      end if
      if (DEV_VERSION_VT) then
         if (associated(node%iPart)) then
            print*,'The tree particles were not properly cleaned'
            STOP
            deallocate(node%iPart)
         end if
      endif
   end subroutine cut_substep
   
   !> Cut a tree and all its sub-branches, unrolled to use parallelization for the first 3 levels
   subroutine cut_tree_parallel(Tree, deallocPart, deallocSgmt)
      type(T_Tree), intent(inout) :: Tree
      logical, optional, intent(in)  :: deallocPart
      logical, optional, intent(in)  :: deallocSgmt
      integer :: i,i1,i2,nBranches,istat
      ! --- Deallocating data we are pointing to, only if user requests it
      if (present(deallocPart)) then
         if (deallocPart) then
            deallocate(Tree%Part%P        , stat=istat)
            deallocate(Tree%Part%AlphA    , stat=istat)
            deallocate(Tree%Part%RegParam , stat=istat)
         endif
      endif
      if (present(deallocSgmt)) then
         if (deallocSgmt) then
            deallocate(Tree%Seg%SP, Tree%Seg%SConnct, Tree%Seg%SGamma, Tree%Seg%RegParam, stat=istat)
         endif
      endif
      ! --- Unlinking particles 
      nullify(Tree%Part%P)
      nullify(Tree%Part%Alpha)
      nullify(Tree%Part%RegParam)
      nullify(Tree%Seg%SP)
      nullify(Tree%Seg%SConnct)
      nullify(Tree%Seg%SGamma)
      nullify(Tree%Seg%RegParam)
      ! ---  Unrolled version of cut_tree for the first node
      call cut_substep(Tree%root)
      if(associated(Tree%Root%branches)) then
         nBranches=size(Tree%Root%branches)
         !$OMP PARALLEL default(shared)

         ! ---  Unrolled version for the second levels
         !$OMP do private(i) schedule(runtime)
         do i = 1,nBranches ! maximum 8 branches
            call cut_substep(Tree%Root%branches(i))
         end do
         !$OMP end do
         !$OMP barrier ! we need to be sure that all the branches were cut

         ! ---  Unrolled version for third node levels
         !$OMP do private(i,i1,i2) schedule(runtime)
         do i = 1,nBranches*8 ! maximum 64 sub branches
            i1=(i-1)/8+1;
            i2=mod(i-1,8)+1;
            if(associated(Tree%Root%branches(i1)%branches)) then
               if (i2<=size(Tree%Root%branches(i1)%branches)) then
                  call cut_tree_rec(Tree%Root%branches(i1)%branches(i2))
               endif
            endif
         enddo
         !$OMP end do 
         !$OMP END PARALLEL

         ! --- Cleanup second level
         do i = 1,nBranches ! maximum 8 branches
            if (associated(Tree%root%branches(i)%branches)) then
               deallocate(Tree%root%branches(i)%branches)
               nullify(Tree%root%branches(i)%branches)
            endif
         end do

         ! --- Cleanup First level
         deallocate(Tree%root%branches)
         nullify(Tree%root%branches)
      endif
      if (DEV_VERSION_VT) then
         if (associated(Tree%root%branches)) then
            print*,'Tree cut: branches are still allocated'
            STOP
         endif
      endif
      Tree%iStep=-1
      Tree%root%nPart=-1
      Tree%bGrown=.false.
   end subroutine cut_tree_parallel

   subroutine print_tree(Tree)
      type(T_Tree) :: Tree
      character(len=1024) :: preffix
      preffix='root'
      print '(A, L1)', trim(preffix)//':partP_assoc = ',associated(Tree%Part%P)
      print '(A, L1)', trim(preffix)//':bGrown      = ',Tree%bGrown
      print '(A, I0)', trim(preffix)//':iStep       = ',Tree%iStep
      call print_tree_rec(Tree%Root, preffix)
   contains
      recursive subroutine print_tree_rec(node, preffix)
         type(T_Node), target     :: node !<
         character(len=*), intent(in) :: preffix
         integer :: i
         ! Test if there are enough particles on the node to build new branchess
         ! The case of only one particle should be handled upstream by allocating one leaf to the parent node
         print'(A, I0)'   ,trim(preffix)//':nPart       = ', node%nPart
         print'(A,3F12.3)',trim(preffix)//':center      =', node%center
         print'(A,1F12.3)',trim(preffix)//':radius      =', node%radius
         if(associated(node%leaves)) then
            do i = 1,size(node%leaves)
               print'(A, I0, A, I0)',trim(preffix)//':leaf', i,'=', node%leaves(i)
            end do
         endif
         if(associated(node%branches)) then
            do i = 1,size(node%branches)
               call print_tree_rec(node%branches(i), trim(preffix)//':branch'//trim(Num2LStr(i)))
            end do
         endif
      end subroutine  print_tree_rec
   end subroutine  print_tree
   ! --------------------------------------------------------------------------------
   ! --- Cut tree
   ! --------------------------------------------------------------------------------
   !> Cut a tree and all its sub-branches, unrolled to use parallelization for the first 3 levels
   subroutine cut_tree_segment_parallel(Tree_Seg)
      type(T_Tree), intent(inout) :: Tree_Seg
      integer :: i,i1,i2,nBranches
      ! --- Unlinking particles
      nullify(Tree_Seg%Seg%SP)
      nullify(Tree_Seg%Seg%SConnct)
      nullify(Tree_Seg%Seg%SGamma)
      nullify(Tree_Seg%Seg%RegParam)
      ! ---  Unrolled version of cut_tree for the first node
      call cut_substep(Tree_Seg%root)
      if(associated(Tree_Seg%Root%branches)) then
         nBranches=size(Tree_Seg%Root%branches)
         !$OMP PARALLEL default(shared)

         ! ---  Unrolled version for the second levels
         !$OMP do private(i) schedule(runtime)
         do i = 1,nBranches ! maximum 8 branches
            call cut_substep(Tree_Seg%Root%branches(i))
         end do
         !$OMP end do
         !$OMP barrier ! we need to be sure that all the branches were cut

         ! ---  Unrolled version for third node levels
         !$OMP do private(i,i1,i2) schedule(runtime)
         do i = 1,nBranches*8 ! maximum 64 sub branches
            i1=(i-1)/8+1;
            i2=mod(i-1,8)+1;
            if(associated(Tree_Seg%Root%branches(i1)%branches)) then
               if (i2<=size(Tree_Seg%Root%branches(i1)%branches)) then
                  call cut_tree_rec(Tree_Seg%Root%branches(i1)%branches(i2))
               endif
            endif
         enddo
         !$OMP end do
         !$OMP END PARALLEL

         ! --- Cleanup second level
         do i = 1,nBranches ! maximum 8 branches
            if (associated(Tree_Seg%root%branches(i)%branches)) then
               deallocate(Tree_Seg%root%branches(i)%branches)
               nullify(Tree_Seg%root%branches(i)%branches)
            endif
         end do

         ! --- Cleanup First level
         deallocate(Tree_Seg%root%branches)
         nullify(Tree_Seg%root%branches)
      endif
      if (DEV_VERSION_VT) then
         if (associated(Tree_Seg%root%branches)) then
            print*,'Tree cut: branches are still allocated'
            STOP
         endif
      endif
      Tree_Seg%iStep=-1
      Tree_Seg%root%nPart=-1
      Tree_Seg%bGrown=.false.
   end subroutine cut_tree_segment_parallel

   ! --------------------------------------------------------------------------------
   ! --- Velocity computation 
   ! --------------------------------------------------------------------------------
   subroutine ui_tree_part(Tree, icp_end, CPs, BranchFactor, DistanceDirect, Uind, ErrStat, ErrMsg)
      use FVW_BiotSavart, only: ui_part_nograd_11
      type(T_Tree), target,          intent(inout) :: Tree            !< 
      integer,                       intent(in   ) :: icp_end         !< Number of CPs to use <size(CPs,2)
      real(ReKi),                    intent(in   ) :: BranchFactor    !<
      real(ReKi),                    intent(in   ) :: DistanceDirect  !< Distance under which direct evaluation should be done no matter what the tree cell size is
      real(ReKi), dimension(:,:),    intent(in   ) :: CPs             !< Control Points  (3 x nCPs)
      real(ReKi), dimension(:,:),    intent(inout) :: Uind            !< Induced velocity at CPs, with side effects (3 x nCPs)
      integer(IntKi),                intent(  out) :: ErrStat         !< Error status of the operation
      character(*),                  intent(  out) :: ErrMsg          !< Error message if ErrStat /= ErrID_None
      real(ReKi), dimension(3) :: Uind_tmp !< 
      real(ReKi), dimension(3) :: CP       !< Current CP
      integer :: icp
      type(T_VPart), pointer :: Part ! Alias
      ErrStat = ErrID_None
      ErrMsg = ''
      Part => Tree%Part
      if(.not. associated(Part%P)) then
         ErrMsg='Ui Part Tree called but tree particles not associated'; ErrStat=ErrID_Fatal; return
      endif
      !$OMP PARALLEL DEFAULT(SHARED)
      !$OMP DO PRIVATE(icp,CP,Uind_tmp) schedule(runtime)
      do icp=1,icp_end
         CP = CPs(1:3,icp)
         Uind_tmp(1:3) = 0.0_ReKi
         call ui_tree_part_11(Tree%root, CP, Uind_tmp) !< SIDE EFFECTS
         Uind(1:3, icp) = Uind(1:3, icp) + Uind_tmp(1:3)
      enddo
      !$OMP END DO 
      !$OMP END PARALLEL
   contains
      !> Velocity at one control point from the entire tree
      recursive subroutine ui_tree_part_11(node, CP, Uind)
         real(ReKi), dimension(3), intent(in   ) :: CP
         real(ReKi), dimension(3), intent(inout) :: Uind  !< Velocity at control point, with side effect
         type(T_Node), intent(inout) :: node
         real(ReKi) :: distDirect
         real(ReKi),dimension(3) :: DeltaP, Uloc
         real(ReKi) :: r
         real(ReKi) :: coeff, coeff3, coeff5, coeff7, coeff7ij
         real(ReKi) :: x, y, z
         real(ReKi),dimension(3) :: phi
         real(ReKi) :: ieqj
         integer :: i,j
         integer :: iPart
         if (node%nPart<=0) then
            ! We skip the dead leaf
         elseif (.not.associated(node%branches)) then
            ! Loop on leaves
            if(associated(node%leaves)) then
               do i =1,size(node%leaves) 
                  iPart=node%leaves(i)
                  DeltaP = CP(1:3) - Part%P(1:3,iPart)
                  call  ui_part_nograd_11(DeltaP, Part%Alpha(1:3,iPart), Part%RegFunction, Part%RegParam(iPart), Uloc)
                  Uind(1:3) = Uind(1:3) + Uloc
               enddo
            endif
         else
            distDirect = max(BranchFactor*node%radius, DistanceDirect) ! Under this distance-> Direct eval., Above it -> quadrupole calculation
            DeltaP  = - node%center + CP(1:3)                          ! Vector between the control point and the center of the branch
            r       = sqrt( DeltaP(1)**2 + DeltaP(2)**2 + DeltaP(3)**2)
            ! Test if the control point is too close from the branch node so that a direct evaluation is needed 
            if (r<distDirect) then
               ! We are too close, perform direct evaluation using children (leaves and branches)
               if(associated(node%leaves)) then
                  do i =1,size(node%leaves) 
                     iPart=node%leaves(i)
                     DeltaP = CP(1:3) - Part%P(1:3,iPart)
                     call  ui_part_nograd_11(DeltaP, Part%Alpha(1:3,iPart), Part%RegFunction, Part%RegParam(iPart), Uloc)
                     Uind(1:3) = Uind(1:3) + Uloc
                  enddo
               endif
               if(associated(node%branches)) then
                  ! TODO: consider implementing a recursive method for that: direct call on all children
                  do i =1,size(node%branches)
                     call ui_tree_part_11(node%branches(i), CP, Uind)
                  end do
               endif

            else 
               ! We are far enough, use branch node quadrupole
               !call ui_expansion_order2(r, DeltaP, node%node%Moments, Uloc) 
               !! NOTE all of this is common between particle and segment
               coeff3 =  fourpi_inv/r**3
               x=DeltaP(1)
               y=DeltaP(2)
               z=DeltaP(3)
               ! Speed, order 0
               phi = node%Moments(1:3,M0)*coeff3
               Uloc(1) = phi(2)*z - phi(3)*y
               Uloc(2) = phi(3)*x - phi(1)*z
               Uloc(3) = phi(1)*y - phi(2)*x
               Uind = Uind + Uloc

               ! Speed, order 1
               Uloc=0.0_ReKi
               coeff5 = 3.0_ReKi*fourpi_inv/r**5
               coeff = coeff5*x
               Uloc(1) = Uloc(1)- coeff*(node%Moments(3,M1_1)*y-node%Moments(2,M1_1)*z)
               Uloc(2) = Uloc(2)- coeff*(node%Moments(1,M1_1)*z-node%Moments(3,M1_1)*x)
               Uloc(3) = Uloc(3)- coeff*(node%Moments(2,M1_1)*x-node%Moments(1,M1_1)*y)

               coeff = coeff5*y
               Uloc(1) = Uloc(1) - coeff*(node%Moments(3,M1_2)*y-node%Moments(2,M1_2)*z)
               Uloc(2) = Uloc(2) - coeff*(node%Moments(1,M1_2)*z-node%Moments(3,M1_2)*x)
               Uloc(3) = Uloc(3) - coeff*(node%Moments(2,M1_2)*x-node%Moments(1,M1_2)*y)

               coeff = coeff5*z
               Uloc(1) = Uloc(1) - coeff*(node%Moments(3,M1_3)*y-node%Moments(2,M1_3)*z)
               Uloc(2) = Uloc(2) - coeff*(node%Moments(1,M1_3)*z-node%Moments(3,M1_3)*x)
               Uloc(3) = Uloc(3) - coeff*(node%Moments(2,M1_3)*x-node%Moments(1,M1_3)*y)

               Uloc(1) = Uloc(1) + coeff3*node%Moments(3,M1_2) - coeff3*node%Moments(2,M1_3)
               Uloc(2) = Uloc(2) + coeff3*node%Moments(1,M1_3) - coeff3*node%Moments(3,M1_1)
               Uloc(3) = Uloc(3) + coeff3*node%Moments(2,M1_1) - coeff3*node%Moments(1,M1_2)
               Uind=Uind+Uloc

               ! Speed, order 2
               Uloc =0.0_ReKi
               coeff = coeff5*0.5_ReKi ! 3/2 * 1/(4 pi r**5)
               coeff7= -7.5_ReKi*fourpi_inv/r**7
               ! TODO: BAR_OLAF test fails (very slightly) with unrolled version...
               !! NOTE: KEEP ME
               do i =1,3
                  Uloc(1) = Uloc(1) + coeff* (y*node%Moments(3,5+2*(i/2)+3*(i/3)) - z*node%Moments(2,5+2*(i/2)+3*(i/3)))
                  Uloc(2) = Uloc(2) + coeff* (z*node%Moments(1,5+2*(i/2)+3*(i/3)) - x*node%Moments(3,5+2*(i/2)+3*(i/3)))
                  Uloc(3) = Uloc(3) + coeff* (x*node%Moments(2,5+2*(i/2)+3*(i/3)) - y*node%Moments(1,5+2*(i/2)+3*(i/3)))
                  do j=1,i
                     if (i==j) then
                        ieqj = 1.0_ReKi
                     else
                        ieqj = 2.0_ReKi
                     end if
                     coeff7ij = DeltaP(i)*DeltaP(j)*coeff7
                     Uloc(1) = Uloc(1) + ieqj * coeff7ij * ( y*node%Moments(3,3+i+j+i/3) - z*node%Moments(2,3+i+j+i/3))
                     Uloc(2) = Uloc(2) + ieqj * coeff7ij * ( z*node%Moments(1,3+i+j+i/3) - x*node%Moments(3,3+i+j+i/3))
                     Uloc(3) = Uloc(3) + ieqj * coeff7ij * ( x*node%Moments(2,3+i+j+i/3) - y*node%Moments(1,3+i+j+i/3))
                  end do
               end do
               !! r5 terms, unrolled version
               !!    Uloc(1) = Uloc(1) + coeff * (y*node%Moments(3,5+2*(i/2)+3*(i/3)) - z*node%Moments(2,5+2*(i/2)+3*(i/3)))
               !!    Uloc(2) = Uloc(2) + coeff * (z*node%Moments(1,5+2*(i/2)+3*(i/3)) - x*node%Moments(3,5+2*(i/2)+3*(i/3)))
               !!    Uloc(3) = Uloc(3) + coeff * (x*node%Moments(2,5+2*(i/2)+3*(i/3)) - y*node%Moments(1,5+2*(i/2)+3*(i/3)))
               ! i=1
               !Uloc(1) = Uloc(1) + coeff    * (y*node%Moments(3,5) - z*node%Moments(2,5))
               !Uloc(2) = Uloc(2) + coeff    * (z*node%Moments(1,5) - x*node%Moments(3,5))
               !Uloc(3) = Uloc(3) + coeff    * (x*node%Moments(2,5) - y*node%Moments(1,5))
               !! i=2
               !Uloc(1) = Uloc(1) + coeff    * (y*node%Moments(3,7) - z*node%Moments(2,7))
               !Uloc(2) = Uloc(2) + coeff    * (z*node%Moments(1,7) - x*node%Moments(3,7))
               !Uloc(3) = Uloc(3) + coeff    * (x*node%Moments(2,7) - y*node%Moments(1,7))
               !! i=3
               !Uloc(1) = Uloc(1) + coeff    * (y*node%Moments(3,8) - z*node%Moments(2,8))
               !Uloc(2) = Uloc(2) + coeff    * (z*node%Moments(1,8) - x*node%Moments(3,8))
               !Uloc(3) = Uloc(3) + coeff    * (x*node%Moments(2,8) - y*node%Moments(1,8))
               !! r7 terms, unrolled version
               !!    coeff7ij = DeltaP(i)*DeltaP(j)*coeff7 * (1 or 2 depdening if i=j)
               !!    Uloc(1) = Uloc(1) + ieqj * coeff7ij * ( y*node%Moments(3,3+i+j+i/3) - z*node%Moments(2,3+i+j+i/3))
               !!    Uloc(2) = Uloc(2) + ieqj * coeff7ij * ( z*node%Moments(1,3+i+j+i/3) - x*node%Moments(3,3+i+j+i/3))
               !!    Uloc(3) = Uloc(3) + ieqj * coeff7ij * ( x*node%Moments(2,3+i+j+i/3) - y*node%Moments(1,3+i+j+i/3))
               !! i=1, j=1
               !coeff7ij = DeltaP(1)*DeltaP(1)*coeff7 !=-7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               !Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,5) - z*node%Moments(2,5))
               !Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,5) - x*node%Moments(3,5))
               !Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,5) - y*node%Moments(1,5))
               !! i=2, j=1
               !coeff7ij = 2.0_ReKi* DeltaP(2)*DeltaP(1)*coeff7 !=-15_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               !Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,6) - z*node%Moments(2,6))
               !Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,6) - x*node%Moments(3,6))
               !Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,6) - y*node%Moments(1,6))
               !! i=2, j=2
               !coeff7ij = DeltaP(2)*DeltaP(2)*coeff7 !=-7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               !Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,7) - z*node%Moments(2,7))
               !Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,7) - x*node%Moments(3,7))
               !Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,7) - y*node%Moments(1,7))
               !! i=3, j=1
               !coeff7ij = 2.0_ReKi*DeltaP(3)*DeltaP(1)*coeff7 !=-15_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               !Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,8) - z*node%Moments(2,8))
               !Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,8) - x*node%Moments(3,8))
               !Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,8) - y*node%Moments(1,8))
               !! i=3, j=2
               !coeff7ij = 2.0_ReKi*DeltaP(3)*DeltaP(2)*coeff7 !=-15_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               !Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,9) - z*node%Moments(2,9))
               !Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,9) - x*node%Moments(3,9))
               !Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,9) - y*node%Moments(1,9))
               !! i=3, j=3
               !coeff7ij = DeltaP(3)*DeltaP(3)*coeff7 !=-7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               !Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,10) - z*node%Moments(2,10))
               !Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,10) - x*node%Moments(3,10))
               !Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,10) - y*node%Moments(1,10))

               ! r5 terms 2
               Uloc(1) = Uloc(1) + coeff5 * (y*node%Moments(3,M2_22) - z*node%Moments(2,M2_33))
               Uloc(2) = Uloc(2) + coeff5 * (z*node%Moments(1,M2_33) - x*node%Moments(3,M2_11))
               Uloc(3) = Uloc(3) + coeff5 * (x*node%Moments(2,M2_11) - y*node%Moments(1,M2_22))
               Uloc(1) = Uloc(1) + coeff5 * (z*node%Moments(3,M2_32) + x*node%Moments(3,M2_21) - x*node%Moments(2,M2_31) - y*node%Moments(2,M2_32))
               Uloc(2) = Uloc(2) + coeff5 * (x*node%Moments(1,M2_31) + y*node%Moments(1,M2_32) - y*node%Moments(3,M2_21) - z*node%Moments(3,M2_31))
               Uloc(3) = Uloc(3) + coeff5 * (y*node%Moments(2,M2_21) + z*node%Moments(2,M2_31) - z*node%Moments(1,M2_32) - x*node%Moments(1,M2_21))
               Uind = Uind + Uloc
            end if ! Far enough
         end if ! had more than 1 particles
      end subroutine ui_tree_part_11
   end subroutine  ui_tree_part

   subroutine ui_tree_segment(Tree, CPs, icp_end, BranchFactor, DistanceDirect, Uind, ErrStat, ErrMsg)
      use FVW_BiotSavart, only: ui_seg_11
      type(T_Tree), target,          intent(inout) :: Tree            !<
      integer,                       intent(in   ) :: icp_end         !<
      real(ReKi),                    intent(in   ) :: BranchFactor    !<
      real(ReKi),                    intent(in   ) :: DistanceDirect  !< Distance under which direct evaluation should be done no matter what the tree cell size is
      real(ReKi), dimension(:,:),    intent(in   ) :: CPs             !< Control Points  (3 x nCPs)
      real(ReKi), dimension(:,:),    intent(inout) :: Uind            !< Induced velocity at CPs, with side effects (3 x nCPs)
      integer(IntKi),                intent(  out) :: ErrStat         !< Error status of the operation
      character(*),                  intent(  out) :: ErrMsg          !< Error message if ErrStat /= ErrID_None
      real(ReKi), dimension(3) :: Uind_tmp !<
      real(ReKi), dimension(3) :: CP       !< Current CP
      integer :: icp
      type(T_VSgmt), pointer :: Seg ! Alias
      Seg => Tree%Seg
      if(.not. associated(Seg%SP)) then
         ErrMsg='Ui Sgmt Tree called but tree segments not associated'; ErrStat=ErrID_Fatal; return
      endif
      !$OMP PARALLEL DEFAULT(SHARED)
      !$OMP DO PRIVATE(icp, CP, Uind_tmp) schedule(runtime)
      do icp=1,icp_end
         CP = CPs(1:3,icp)
         Uind_tmp(1:3) = 0.0_ReKi
         call ui_tree_segment_11(Tree%root, CP, Uind_tmp) !< SIDE EFFECTS
         Uind(1:3, icp) = Uind(1:3, icp) + Uind_tmp(1:3)
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
   contains
      !> Velocity at one control point from the entire tree
      recursive subroutine ui_tree_segment_11(node, CP, Uind)
         real(ReKi), dimension(3), intent(in   ) :: CP
         real(ReKi), dimension(3), intent(inout) :: Uind  !< Velocity at control point, with side effect
         type(T_Node), intent(inout) :: node
         real(ReKi) :: distDirect
         real(ReKi),dimension(3) :: DeltaP, DeltaPa, DeltaPb, Uloc
         real(ReKi),dimension(3) :: Uquad
         real(ReKi) :: r
         real(ReKi) :: coeff, coeff3, coeff5, coeff7, coeff7ij
         real(ReKi) :: x, y, z
         real(ReKi),dimension(3) :: phi
         integer :: i
         integer :: iPart
         if (node%nPart<=0) then
            ! We skip the dead leaf
         elseif (.not.associated(node%branches)) then
            ! Loop on leaves
            if(associated(node%leaves)) then
               do i =1,size(node%leaves)
                  iPart=node%leaves(i)
                  DeltaPa = CP(1:3) - Seg%SP(1:3,Seg%SConnct(1,iPart))
                  DeltaPb = CP(1:3) - Seg%SP(1:3,Seg%SConnct(2,iPart))
                  call ui_seg_11(DeltaPa, DeltaPb, Seg%SGamma(iPart), Seg%RegFunction, Seg%RegParam(iPart), Uloc)
                  Uind(1:3) = Uind(1:3) + Uloc
               enddo
            endif
         else
            distDirect = max(BranchFactor*node%radius, DistanceDirect) ! Under this distance-> Direct eval., Above it -> quadrupole calculation
            DeltaP  = - node%center + CP(1:3)                          ! Vector between the control point and the center of the branch
            r       = sqrt( DeltaP(1)**2 + DeltaP(2)**2 + DeltaP(3)**2)
            ! Test if the control point is too close from the branch node so that a direct evaluation is needed
            if (r<distDirect) then
               ! We are too close, perform direct evaluation using children (leaves and branches)
               if(associated(node%leaves)) then
                  do i =1,size(node%leaves)
                     iPart=node%leaves(i)
                     DeltaPa = CP(1:3) - Seg%SP(1:3,Seg%SConnct(1,iPart))
                     DeltaPb = CP(1:3) - Seg%SP(1:3,Seg%SConnct(2,iPart))
                     call ui_seg_11(DeltaPa, DeltaPb, Seg%SGamma(iPart), Seg%RegFunction, Seg%RegParam(iPart), Uloc)
                     Uind(1:3) = Uind(1:3) + Uloc
                  enddo
               endif
               if(associated(node%branches)) then
                  ! TODO: consider implementing a recursive method for that: direct call on all children
                  do i =1,size(node%branches)
                     call ui_tree_segment_11(node%branches(i), CP, Uind)
                  end do
               endif

            else
               ! We are far enough, use branch node quadrupole
               !call ui_expansion_order2(r, DeltaP, node%node%Moments, Uloc) 
               !! NOTE all of this is common between particle and segment
               coeff3 =  fourpi_inv/r**3
               x=DeltaP(1)
               y=DeltaP(2)
               z=DeltaP(3)
               ! Speed, order 0
               phi = node%Moments(1:3,M0)*coeff3
               Uloc(1) = phi(2)*z - phi(3)*y
               Uloc(2) = phi(3)*x - phi(1)*z
               Uloc(3) = phi(1)*y - phi(2)*x
               Uquad = Uloc

               ! Speed, order 1
               Uloc=0.0_ReKi
               coeff5 = 3.0_ReKi*fourpi_inv/r**5
               coeff = coeff5*x
               Uloc(1) = Uloc(1)- coeff*(node%Moments(3,M1_1)*y-node%Moments(2,M1_1)*z)
               Uloc(2) = Uloc(2)- coeff*(node%Moments(1,M1_1)*z-node%Moments(3,M1_1)*x)
               Uloc(3) = Uloc(3)- coeff*(node%Moments(2,M1_1)*x-node%Moments(1,M1_1)*y)

               coeff = coeff5*y
               Uloc(1) = Uloc(1) - coeff*(node%Moments(3,M1_2)*y-node%Moments(2,M1_2)*z)
               Uloc(2) = Uloc(2) - coeff*(node%Moments(1,M1_2)*z-node%Moments(3,M1_2)*x)
               Uloc(3) = Uloc(3) - coeff*(node%Moments(2,M1_2)*x-node%Moments(1,M1_2)*y)

               coeff = coeff5*z
               Uloc(1) = Uloc(1) - coeff*(node%Moments(3,M1_3)*y-node%Moments(2,M1_3)*z)
               Uloc(2) = Uloc(2) - coeff*(node%Moments(1,M1_3)*z-node%Moments(3,M1_3)*x)
               Uloc(3) = Uloc(3) - coeff*(node%Moments(2,M1_3)*x-node%Moments(1,M1_3)*y)

               Uloc(1) = Uloc(1) + coeff3*node%Moments(3,M1_2) - coeff3*node%Moments(2,M1_3)
               Uloc(2) = Uloc(2) + coeff3*node%Moments(1,M1_3) - coeff3*node%Moments(3,M1_1)
               Uloc(3) = Uloc(3) + coeff3*node%Moments(2,M1_1) - coeff3*node%Moments(1,M1_2)
               Uquad = Uquad + Uloc

               ! Speed, order 2
               Uloc =0.0_ReKi
               coeff = coeff5*0.5_ReKi ! 3/2 * 1/(4 pi r**5)
               coeff7= -7.5_ReKi*fourpi_inv/r**7
               !! NOTE: KEEP ME
               !do i =1,3
               !   Uloc(1) = Uloc(1) + coeff* (y*node%Moments(3,5+2*(i/2)+3*(i/3)) - z*node%Moments(2,5+2*(i/2)+3*(i/3)))
               !   Uloc(2) = Uloc(2) + coeff* (z*node%Moments(1,5+2*(i/2)+3*(i/3)) - x*node%Moments(3,5+2*(i/2)+3*(i/3)))
               !   Uloc(3) = Uloc(3) + coeff* (x*node%Moments(2,5+2*(i/2)+3*(i/3)) - y*node%Moments(1,5+2*(i/2)+3*(i/3)))
               !   do j=1,i
               !      if (i==j) then
               !         ieqj = 1.0_ReKi
               !      else
               !         ieqj = 2.0_ReKi
               !      end if
               !      coeff7ij = DeltaP(i)*DeltaP(j)*coeff7
               !      Uloc(1) = Uloc(1) + ieqj * coeff7ij * ( y*node%Moments(3,3+i+j+i/3) - z*node%Moments(2,3+i+j+i/3))
               !      Uloc(2) = Uloc(2) + ieqj * coeff7ij * ( z*node%Moments(1,3+i+j+i/3) - x*node%Moments(3,3+i+j+i/3))
               !      Uloc(3) = Uloc(3) + ieqj * coeff7ij * ( x*node%Moments(2,3+i+j+i/3) - y*node%Moments(1,3+i+j+i/3))
               !   end do
               !end do
               !! r5 terms, unrolled version
               !!    Uloc(1) = Uloc(1) + coeff * (y*node%Moments(3,5+2*(i/2)+3*(i/3)) - z*node%Moments(2,5+2*(i/2)+3*(i/3)))
               !!    Uloc(2) = Uloc(2) + coeff * (z*node%Moments(1,5+2*(i/2)+3*(i/3)) - x*node%Moments(3,5+2*(i/2)+3*(i/3)))
               !!    Uloc(3) = Uloc(3) + coeff * (x*node%Moments(2,5+2*(i/2)+3*(i/3)) - y*node%Moments(1,5+2*(i/2)+3*(i/3)))
               ! i=1
               Uloc(1) = Uloc(1) + coeff    * (y*node%Moments(3,5) - z*node%Moments(2,5))
               Uloc(2) = Uloc(2) + coeff    * (z*node%Moments(1,5) - x*node%Moments(3,5))
               Uloc(3) = Uloc(3) + coeff    * (x*node%Moments(2,5) - y*node%Moments(1,5))
               ! i=2
               Uloc(1) = Uloc(1) + coeff    * (y*node%Moments(3,7) - z*node%Moments(2,7))
               Uloc(2) = Uloc(2) + coeff    * (z*node%Moments(1,7) - x*node%Moments(3,7))
               Uloc(3) = Uloc(3) + coeff    * (x*node%Moments(2,7) - y*node%Moments(1,7))
               ! i=3
               Uloc(1) = Uloc(1) + coeff    * (y*node%Moments(3,8) - z*node%Moments(2,8))
               Uloc(2) = Uloc(2) + coeff    * (z*node%Moments(1,8) - x*node%Moments(3,8))
               Uloc(3) = Uloc(3) + coeff    * (x*node%Moments(2,8) - y*node%Moments(1,8))
               ! r7 terms, unrolled version
               !    coeff7ij = DeltaP(i)*DeltaP(j)*coeff7 * (1 or 2 depdening if i=j)
               !    Uloc(1) = Uloc(1) + ieqj * coeff7ij * ( y*node%Moments(3,3+i+j+i/3) - z*node%Moments(2,3+i+j+i/3))
               !    Uloc(2) = Uloc(2) + ieqj * coeff7ij * ( z*node%Moments(1,3+i+j+i/3) - x*node%Moments(3,3+i+j+i/3))
               !    Uloc(3) = Uloc(3) + ieqj * coeff7ij * ( x*node%Moments(2,3+i+j+i/3) - y*node%Moments(1,3+i+j+i/3))
               ! i=1, j=1
               coeff7ij = DeltaP(1)*DeltaP(1)*coeff7 !=-7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,5) - z*node%Moments(2,5))
               Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,5) - x*node%Moments(3,5))
               Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,5) - y*node%Moments(1,5))
               ! i=2, j=1
               coeff7ij = 2.0_ReKi* DeltaP(2)*DeltaP(1)*coeff7 !=-15_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,6) - z*node%Moments(2,6))
               Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,6) - x*node%Moments(3,6))
               Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,6) - y*node%Moments(1,6))
               ! i=2, j=2
               coeff7ij = DeltaP(2)*DeltaP(2)*coeff7 !=-7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,7) - z*node%Moments(2,7))
               Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,7) - x*node%Moments(3,7))
               Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,7) - y*node%Moments(1,7))
               ! i=3, j=1
               coeff7ij = 2.0_ReKi*DeltaP(3)*DeltaP(1)*coeff7 !=-15_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,8) - z*node%Moments(2,8))
               Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,8) - x*node%Moments(3,8))
               Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,8) - y*node%Moments(1,8))
               ! i=3, j=2
               coeff7ij = 2.0_ReKi*DeltaP(3)*DeltaP(2)*coeff7 !=-15_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,9) - z*node%Moments(2,9))
               Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,9) - x*node%Moments(3,9))
               Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,9) - y*node%Moments(1,9))
               ! i=3, j=3
               coeff7ij = DeltaP(3)*DeltaP(3)*coeff7 !=-7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
               Uloc(1) = Uloc(1) + coeff7ij * (y*node%Moments(3,10) - z*node%Moments(2,10))
               Uloc(2) = Uloc(2) + coeff7ij * (z*node%Moments(1,10) - x*node%Moments(3,10))
               Uloc(3) = Uloc(3) + coeff7ij * (x*node%Moments(2,10) - y*node%Moments(1,10))

               ! r5 terms 2
               Uloc(1) = Uloc(1) + coeff5 * (y*node%Moments(3,M2_22) - z*node%Moments(2,M2_33))
               Uloc(2) = Uloc(2) + coeff5 * (z*node%Moments(1,M2_33) - x*node%Moments(3,M2_11))
               Uloc(3) = Uloc(3) + coeff5 * (x*node%Moments(2,M2_11) - y*node%Moments(1,M2_22))
               Uloc(1) = Uloc(1) + coeff5 * (z*node%Moments(3,M2_32) + x*node%Moments(3,M2_21) - x*node%Moments(2,M2_31) - y*node%Moments(2,M2_32))
               Uloc(2) = Uloc(2) + coeff5 * (x*node%Moments(1,M2_31) + y*node%Moments(1,M2_32) - y*node%Moments(3,M2_21) - z*node%Moments(3,M2_31))
               Uloc(3) = Uloc(3) + coeff5 * (y*node%Moments(2,M2_21) + z*node%Moments(2,M2_31) - z*node%Moments(1,M2_32) - x*node%Moments(1,M2_21))
               Uquad = Uquad + Uloc

               Uind = Uind + Uquad
            end if ! Far enough
         end if ! had more than 1 particles
      end subroutine ui_tree_segment_11
   end subroutine  ui_tree_segment


   !> Compute second order expansion of velocity for a given point based on moment up to order 2 (quadrupole)
   !! TODO: might be sensitive to numerics since r,x,y,z might be large
   pure subroutine ui_expansion_order2(r, DeltaP, Moments, Uind) 
      real(ReKi),                  intent(in   ) :: r      !< Norm of DeltaP
      real(ReKi), dimension(3),    intent(in   ) :: DeltaP !< ControlPoint - Vorticity Center
      real(ReKi), dimension(3,10), intent(in   ) :: Moments!< Moments from Taylor expansion
      real(ReKi), dimension(3),    intent(  out) :: Uind   !< Velocity
      real(ReKi) :: coeff, coeff3, coeff5, coeff7, coeff7ij
      real(ReKi) :: x, y, z
      real(ReKi),dimension(3) :: phi
      real(ReKi),dimension(3) :: Uloc
      !integer :: i,j
      !real(ReKi) :: ieqj
      Uind=0.0_ReKi

      !! NOTE all of this is common between particle and segment
      coeff3 =  fourpi_inv/r**3
      x=DeltaP(1)
      y=DeltaP(2)
      z=DeltaP(3)
      ! Speed, order 0
      phi = Moments(1:3,M0)*coeff3
      Uloc(1) = phi(2)*z - phi(3)*y
      Uloc(2) = phi(3)*x - phi(1)*z
      Uloc(3) = phi(1)*y - phi(2)*x
      Uind = Uind + Uloc

      ! Speed, order 1
      Uloc=0.0_ReKi
      coeff5 = 3.0_ReKi*fourpi_inv/r**5
      coeff = coeff5*x
      Uloc(1) = Uloc(1)- coeff*(Moments(3,M1_1)*y-Moments(2,M1_1)*z)
      Uloc(2) = Uloc(2)- coeff*(Moments(1,M1_1)*z-Moments(3,M1_1)*x)
      Uloc(3) = Uloc(3)- coeff*(Moments(2,M1_1)*x-Moments(1,M1_1)*y)

      coeff = coeff5*y
      Uloc(1) = Uloc(1) - coeff*(Moments(3,M1_2)*y-Moments(2,M1_2)*z)
      Uloc(2) = Uloc(2) - coeff*(Moments(1,M1_2)*z-Moments(3,M1_2)*x)
      Uloc(3) = Uloc(3) - coeff*(Moments(2,M1_2)*x-Moments(1,M1_2)*y)

      coeff = coeff5*z
      Uloc(1) = Uloc(1) - coeff*(Moments(3,M1_3)*y-Moments(2,M1_3)*z)
      Uloc(2) = Uloc(2) - coeff*(Moments(1,M1_3)*z-Moments(3,M1_3)*x)
      Uloc(3) = Uloc(3) - coeff*(Moments(2,M1_3)*x-Moments(1,M1_3)*y)

      Uloc(1) = Uloc(1) + coeff3*Moments(3,M1_2) - coeff3*Moments(2,M1_3)
      Uloc(2) = Uloc(2) + coeff3*Moments(1,M1_3) - coeff3*Moments(3,M1_1)
      Uloc(3) = Uloc(3) + coeff3*Moments(2,M1_1) - coeff3*Moments(1,M1_2)
      Uind=Uind+Uloc

      ! Speed, order 2
      Uloc =0.0_ReKi
      coeff = coeff5*0.5_ReKi ! 3/2 * 1/(4 pi r**5)
      coeff7= -7.5_ReKi*fourpi_inv/r**7
      !! NOTE: KEEP ME
      !do i =1,3
      !   Uloc(1) = Uloc(1) + coeff* (y*Moments(3,5+2*(i/2)+3*(i/3)) - z*Moments(2,5+2*(i/2)+3*(i/3)))
      !   Uloc(2) = Uloc(2) + coeff* (z*Moments(1,5+2*(i/2)+3*(i/3)) - x*Moments(3,5+2*(i/2)+3*(i/3)))
      !   Uloc(3) = Uloc(3) + coeff* (x*Moments(2,5+2*(i/2)+3*(i/3)) - y*Moments(1,5+2*(i/2)+3*(i/3)))
      !   do j=1,i
      !      if (i==j) then
      !         ieqj = 1.0_ReKi
      !      else
      !         ieqj = 2.0_ReKi
      !      end if
      !      coeff7ij = DeltaP(i)*DeltaP(j)*coeff7
      !      Uloc(1) = Uloc(1) + ieqj * coeff7ij * ( y*Moments(3,3+i+j+i/3) - z*Moments(2,3+i+j+i/3))
      !      Uloc(2) = Uloc(2) + ieqj * coeff7ij * ( z*Moments(1,3+i+j+i/3) - x*Moments(3,3+i+j+i/3))
      !      Uloc(3) = Uloc(3) + ieqj * coeff7ij * ( x*Moments(2,3+i+j+i/3) - y*Moments(1,3+i+j+i/3))
      !   end do
      !end do
      !! r5 terms, unrolled version
      !!    Uloc(1) = Uloc(1) + coeff * (y*Moments(3,5+2*(i/2)+3*(i/3)) - z*Moments(2,5+2*(i/2)+3*(i/3)))
      !!    Uloc(2) = Uloc(2) + coeff * (z*Moments(1,5+2*(i/2)+3*(i/3)) - x*Moments(3,5+2*(i/2)+3*(i/3)))
      !!    Uloc(3) = Uloc(3) + coeff * (x*Moments(2,5+2*(i/2)+3*(i/3)) - y*Moments(1,5+2*(i/2)+3*(i/3)))
      ! i=1
      Uloc(1) = Uloc(1) + coeff    * (y*Moments(3,5) - z*Moments(2,5))
      Uloc(2) = Uloc(2) + coeff    * (z*Moments(1,5) - x*Moments(3,5))
      Uloc(3) = Uloc(3) + coeff    * (x*Moments(2,5) - y*Moments(1,5))
      ! i=2
      Uloc(1) = Uloc(1) + coeff    * (y*Moments(3,7) - z*Moments(2,7))
      Uloc(2) = Uloc(2) + coeff    * (z*Moments(1,7) - x*Moments(3,7))
      Uloc(3) = Uloc(3) + coeff    * (x*Moments(2,7) - y*Moments(1,7))
      ! i=3
      Uloc(1) = Uloc(1) + coeff    * (y*Moments(3,8) - z*Moments(2,8))
      Uloc(2) = Uloc(2) + coeff    * (z*Moments(1,8) - x*Moments(3,8))
      Uloc(3) = Uloc(3) + coeff    * (x*Moments(2,8) - y*Moments(1,8))
      ! r7 terms, unrolled version
      !    coeff7ij = DeltaP(i)*DeltaP(j)*coeff7 * (1 or 2 depdening if i=j)
      !    Uloc(1) = Uloc(1) + ieqj * coeff7ij * ( y*Moments(3,3+i+j+i/3) - z*Moments(2,3+i+j+i/3))
      !    Uloc(2) = Uloc(2) + ieqj * coeff7ij * ( z*Moments(1,3+i+j+i/3) - x*Moments(3,3+i+j+i/3))
      !    Uloc(3) = Uloc(3) + ieqj * coeff7ij * ( x*Moments(2,3+i+j+i/3) - y*Moments(1,3+i+j+i/3))
      ! i=1, j=1
      coeff7ij = DeltaP(1)*DeltaP(1)*coeff7 !=-7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
      Uloc(1) = Uloc(1) + coeff7ij * (y*Moments(3,5) - z*Moments(2,5))
      Uloc(2) = Uloc(2) + coeff7ij * (z*Moments(1,5) - x*Moments(3,5))
      Uloc(3) = Uloc(3) + coeff7ij * (x*Moments(2,5) - y*Moments(1,5))
      ! i=2, j=1
      coeff7ij = 2.0_ReKi* DeltaP(2)*DeltaP(1)*coeff7 !=-15_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
      Uloc(1) = Uloc(1) + coeff7ij * (y*Moments(3,6) - z*Moments(2,6))
      Uloc(2) = Uloc(2) + coeff7ij * (z*Moments(1,6) - x*Moments(3,6))
      Uloc(3) = Uloc(3) + coeff7ij * (x*Moments(2,6) - y*Moments(1,6))
      ! i=2, j=2
      coeff7ij = DeltaP(2)*DeltaP(2)*coeff7 !=-7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
      Uloc(1) = Uloc(1) + coeff7ij * (y*Moments(3,7) - z*Moments(2,7))
      Uloc(2) = Uloc(2) + coeff7ij * (z*Moments(1,7) - x*Moments(3,7))
      Uloc(3) = Uloc(3) + coeff7ij * (x*Moments(2,7) - y*Moments(1,7))
      ! i=3, j=1
      coeff7ij = 2.0_ReKi*DeltaP(3)*DeltaP(1)*coeff7 !=-15_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
      Uloc(1) = Uloc(1) + coeff7ij * (y*Moments(3,8) - z*Moments(2,8))
      Uloc(2) = Uloc(2) + coeff7ij * (z*Moments(1,8) - x*Moments(3,8))
      Uloc(3) = Uloc(3) + coeff7ij * (x*Moments(2,8) - y*Moments(1,8))
      ! i=3, j=2
      coeff7ij = 2.0_ReKi*DeltaP(3)*DeltaP(2)*coeff7 !=-15_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
      Uloc(1) = Uloc(1) + coeff7ij * (y*Moments(3,9) - z*Moments(2,9))
      Uloc(2) = Uloc(2) + coeff7ij * (z*Moments(1,9) - x*Moments(3,9))
      Uloc(3) = Uloc(3) + coeff7ij * (x*Moments(2,9) - y*Moments(1,9))
      ! i=3, j=3
      coeff7ij = DeltaP(3)*DeltaP(3)*coeff7 !=-7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
      Uloc(1) = Uloc(1) + coeff7ij * (y*Moments(3,10) - z*Moments(2,10))
      Uloc(2) = Uloc(2) + coeff7ij * (z*Moments(1,10) - x*Moments(3,10))
      Uloc(3) = Uloc(3) + coeff7ij * (x*Moments(2,10) - y*Moments(1,10))

      ! r5 terms 2
      Uloc(1) = Uloc(1) + coeff5 * (y*Moments(3,M2_22) - z*Moments(2,M2_33))
      Uloc(2) = Uloc(2) + coeff5 * (z*Moments(1,M2_33) - x*Moments(3,M2_11))
      Uloc(3) = Uloc(3) + coeff5 * (x*Moments(2,M2_11) - y*Moments(1,M2_22))
      Uloc(1) = Uloc(1) + coeff5 * (z*Moments(3,M2_32) + x*Moments(3,M2_21) - x*Moments(2,M2_31) - y*Moments(2,M2_32))
      Uloc(2) = Uloc(2) + coeff5 * (x*Moments(1,M2_31) + y*Moments(1,M2_32) - y*Moments(3,M2_21) - z*Moments(3,M2_31))
      Uloc(3) = Uloc(3) + coeff5 * (y*Moments(2,M2_21) + z*Moments(2,M2_31) - z*Moments(1,M2_32) - x*Moments(1,M2_21))
      Uind = Uind + Uloc
   end subroutine ui_expansion_order2


   ! --------------------------------------------------------------------------------
   ! --- Vector analysis tools 
   ! --------------------------------------------------------------------------------
   !> Compute the curl of a 3d vector field (stored on a regular grid) using 4th order finite differences
   subroutine curl_regular_grid(f,rotf,ix,iy,iz,nx,ny,nz,dx,dy,dz)
      ! Arguments
      integer, intent(in)  :: nx,ny,nz !< Dimensions of grid
      integer, intent(in)  :: ix,iy,iz !< Dimensions of grid
      real(ReKi), intent(in) :: dx,dy,dz !< grid spacing
      real(ReKi), dimension(:,:,:,:), intent(in)  :: f   !< vector field
      real(ReKi), dimension(:,:,:,:), intent(out) :: rotf !< curl of f
      ! Variables
      integer :: i,j,k,l
      real(ReKi),dimension(3) :: dfi_dx, dfi_dy, dfi_dz
      real(ReKi) dx2,dy2,dz2
      dx2=2.0_ReKi*dx
      dy2=2.0_ReKi*dy
      dz2=2.0_ReKi*dz

      ! left boundary: forward difference   (a=0, b=5) and (a=1, b=3)&
      ! right boundary: backward difference (a=3, b=1) and (a=0, b=5)&
      ! elsewhere: centered difference (a=2,b=2)&
      do k = iz,nz
         do j = iy,ny
            do i = ix,nx
               ! x-derivatives of fy and fz
               if ((nx-ix)<=3) then
                  dfi_dx=0.0_ReKi
               else
                  do l=2,3
                     if (i==ix) then
                        dfi_dx(l) = ( - 25.0_ReKi/6.0_ReKi * f(l,i  ,j,k)&
                            &         +           8.0_ReKi * f(l,i+1,j,k)&
                            &         -           6.0_ReKi * f(l,i+2,j,k)&
                            &         +  8.0_ReKi/3.0_ReKi * f(l,i+3,j,k)&
                            &         -  1.0_ReKi/2.0_ReKi * f(l,i+4,j,k))/dx2
                     elseif (i==ix+1) then
                        dfi_dx(l)  = ( - 1.0_ReKi/2.0_ReKi * f(l,i-1,j,k)&
                           &           - 5.0_ReKi/3.0_ReKi * f(l,i  ,j,k)&
                           &           +          3.0_ReKi * f(l,i+1,j,k)&
                           &           -          1.0_ReKi * f(l,i+2,j,k)&
                           &           + 1.0_ReKi/6.0_ReKi * f(l,i+3,j,k))/dx2
                     elseif (i==nx-1) then
                        dfi_dx(l)  = ( - 1.0_ReKi/6.0_ReKi * f(l,i-3,j,k)&
                           &           +          1.0_ReKi * f(l,i-2,j,k)&
                           &           -          3.0_ReKi * f(l,i-1,j,k)&
                           &           + 5.0_ReKi/3.0_ReKi * f(l,i  ,j,k)&
                           &           + 1.0_ReKi/2.0_ReKi * f(l,i+1,j,k))/dx2
                     elseif (i==nx) then
                        dfi_dx(l) = (    1.0_ReKi/2.0_ReKi * f(l,i-4,j,k)&
                           &          -  8.0_ReKi/3.0_ReKi * f(l,i-3,j,k)&
                           &          +           6.0_ReKi * f(l,i-2,j,k)&
                           &          -           8.0_ReKi * f(l,i-1,j,k)&
                           &          + 25.0_ReKi/6.0_ReKi * f(l,i  ,j,k))/dx2
                     else
                        dfi_dx(l)  = (   1.0_ReKi/6.0_ReKi * f(l,i-2,j,k)&
                           &           - 4.0_ReKi/3.0_ReKi * f(l,i-1,j,k)&
                           &           + 4.0_ReKi/3.0_ReKi * f(l,i+1,j,k)&
                           &           - 1.0_ReKi/6.0_ReKi * f(l,i+2,j,k))/dx2
                     end if
                  end do
               endif

               ! y-derivatives of fx and fz
               if ((ny-iy)<=3) then
                  dfi_dy=0.0_ReKi
               else
                  do l=1,3,2
                     if (j==iy) then
                        dfi_dy(l) = ( - 25.0_ReKi/6.0_ReKi * f(l,i,j  ,k)&
                           &          +           8.0_ReKi * f(l,i,j+1,k)&
                           &          -           6.0_ReKi * f(l,i,j+2,k)&
                           &          +  8.0_ReKi/3.0_ReKi * f(l,i,j+3,k)&
                           &          -  1.0_ReKi/2.0_ReKi * f(l,i,j+4,k))/dy2
                     elseif (j==iy+1) then
                        dfi_dy(l)  = ( - 1.0_ReKi/2.0_ReKi * f(l,i,j-1,k)&
                           &           - 5.0_ReKi/3.0_ReKi * f(l,i,j  ,k)&
                           &           +          3.0_ReKi * f(l,i,j+1,k)&
                           &           -          1.0_ReKi * f(l,i,j+2,k)&
                           &           + 1.0_ReKi/6.0_ReKi * f(l,i,j+3,k))/dy2
                     elseif (j==ny-1) then
                        dfi_dy(l)  = ( - 1.0_ReKi/6.0_ReKi * f(l,i,j-3,k)&
                           &           +          1.0_ReKi * f(l,i,j-2,k)&
                           &           -          3.0_ReKi * f(l,i,j-1,k)&
                           &           + 5.0_ReKi/3.0_ReKi * f(l,i,j  ,k)&
                           &           + 1.0_ReKi/2.0_ReKi * f(l,i,j+1,k))/dy2
                     elseif (j==ny) then
                        dfi_dy(l) = (    1.0_ReKi/2.0_ReKi * f(l,i,j-4,k)&
                           &          -  8.0_ReKi/3.0_ReKi * f(l,i,j-3,k)&
                           &          +           6.0_ReKi * f(l,i,j-2,k)&
                           &          -           8.0_ReKi * f(l,i,j-1,k)&
                           &          + 25.0_ReKi/6.0_ReKi * f(l,i,j  ,k))/dy2
                     else
                        dfi_dy(l)  = (   1.0_ReKi/6.0_ReKi * f(l,i,j-2,k)&
                           &           - 4.0_ReKi/3.0_ReKi * f(l,i,j-1,k)&
                           &           + 4.0_ReKi/3.0_ReKi * f(l,i,j+1,k)&
                           &           - 1.0_ReKi/6.0_ReKi * f(l,i,j+2,k))/dy2
                     end if
                  end do
               endif

               ! z-derivatives of fy and fx
               if ((nz-iz)<=3) then
                  dfi_dz=0.0_ReKi
               else
                  do l=1,2
                     if (k==iz) then
                        dfi_dz(l) = ( - 25.0_ReKi/6.0_ReKi * f(l,i,j,k  )&
                           &          +           8.0_ReKi * f(l,i,j,k+1)&
                           &          -           6.0_ReKi * f(l,i,j,k+2)&
                           &          +  8.0_ReKi/3.0_ReKi * f(l,i,j,k+3)&
                           &          -  1.0_ReKi/2.0_ReKi * f(l,i,j,k+4))/dz2
                     elseif (k==iz+1) then
                        dfi_dz(l)  = ( - 1.0_ReKi/2.0_ReKi * f(l,i,j,k-1)&
                           &           - 5.0_ReKi/3.0_ReKi * f(l,i,j,k  )&
                           &           +          3.0_ReKi * f(l,i,j,k+1)&
                           &           -          1.0_ReKi * f(l,i,j,k+2)&
                           &           + 1.0_ReKi/6.0_ReKi * f(l,i,j,k+3))/dz2
                     elseif (k==nz-1) then
                        dfi_dz(l)  = ( - 1.0_ReKi/6.0_ReKi * f(l,i,j,k-3)&
                           &           +          1.0_ReKi * f(l,i,j,k-2)&
                           &           -          3.0_ReKi * f(l,i,j,k-1)&
                           &           + 5.0_ReKi/3.0_ReKi * f(l,i,j,k  )&
                           &           + 1.0_ReKi/2.0_ReKi * f(l,i,j,k+1))/dz2

                     elseif (k==nz) then
                        dfi_dz(l)  = (    1.0_ReKi/2.0_ReKi * f(l,i,j,k-4)&
                           &           -  8.0_ReKi/3.0_ReKi * f(l,i,j,k-3)&
                           &           +           6.0_ReKi * f(l,i,j,k-2)&
                           &           -           8.0_ReKi * f(l,i,j,k-1)&
                           &           + 25.0_ReKi/6.0_ReKi * f(l,i,j,k  ))/dz2
                     else
                        dfi_dz(l) = (   1.0_ReKi/6.0_ReKi * f(l,i,j,k-2)&
                           &          - 4.0_ReKi/3.0_ReKi * f(l,i,j,k-1)&
                           &          + 4.0_ReKi/3.0_ReKi * f(l,i,j,k+1)&
                           &          - 1.0_ReKi/6.0_ReKi * f(l,i,j,k+2))/dz2
                     end if
                  end do
               endif
               ! Rotational
               rotf(1,i,j,k)= dfi_dy(3)-dfi_dz(2)
               rotf(2,i,j,k)= dfi_dz(1)-dfi_dx(3)
               rotf(3,i,j,k)= dfi_dx(2)-dfi_dy(1)
            end do !i
         end do ! j
      end do ! k
      ! If not enough data..
      !if ((nx-ix)<=3) then
      !   rotf(2,:,:,:)=0.0_ReKi
      !   rotf(3,:,:,:)=0.0_ReKi
      !endif
      !if ((ny-iy)<=3) then
      !   rotf(1,:,:,:)=0.0_ReKi
      !   rotf(3,:,:,:)=0.0_ReKi
      !endif
      !if ((nz-iz)<=3) then
      !   rotf(1,:,:,:)=0.0_ReKi
      !   rotf(2,:,:,:)=0.0_ReKi
      !endif
   end subroutine curl_regular_grid



   ! --- TIC TOC MODULE
   !> Simpler version of matlab tic
   subroutine tic(label)
      character(len=*),intent(in),optional ::label !< Optional label will be displayed when calling toc
      character(len=nmax_label)::lbl
      integer :: i
      if (npos<nmax_store ) then
         if (present(label)) then
            do i=1,nmax_label
               if(i<=npos) then
                  lbl(i:i)='.'
               else
                  if(i-npos<=len_trim(label).and. i-npos>0)  then
                     lbl(i:i)=label((i-npos):(i-npos))
                  else
                     lbl(i:i)=' '
                  endif
               endif
            enddo
         else
            lbl=''
         endif
         !
         npos=npos+1
         call date_and_time(values=time_array)
         start_arrays(npos,1:8) = time_array(1:8) !< we store the whole array
         labels(npos)=lbl
      else
         write(*,*) 'TicToc: stack full'
      endif

   end subroutine tic

   !> Simpler version of matlab toc. Computes elapsed wallclock after a call to tic()
   subroutine toc(delta_t_out)
      real(ReKi), intent(out), optional ::  delta_t_out
      real(ReKi) ::  delta_t
      integer, dimension(8) :: v_dt
      if (npos<=nmax_store.and.npos>0 ) then
         call date_and_time(values=time_array)
         v_dt=time_array-start_arrays(npos,1:8)
         if (v_dt(4)<0. ) then
            v_dt(6)=v_dt(6)+30 ! approximate month change
         endif
         delta_t= v_dt(4)*86400+v_dt(5)*3600+v_dt(6)*60 +v_dt(7)+0.001*v_dt(8)
         if (.not. bSilentTicToc) then 
            if((bSmartTicToc .and. delta_t>SmartTicTocVal) .or.(.not.bSmartTicToc)) then
               write (6, '(A,A,A,A)') 'Time:    ',labels(npos),'- Time: ', pretty_time(delta_t)
            endif
         endif
         npos=npos-1
      else
         write(*,*) 'TicToc: stack error'
         delta_t=0
      endif
      if(present(delta_t_out)) then
         delta_t_out=delta_t
      endif
   end subroutine toc

   function pretty_time(t) 
      character(len=6) :: pretty_time
      real(ReKi), intent(in) :: t
      integer :: d,m,h,s,c
      if(t<0) then
         write(pretty_time,'(A)') '------'
      elseif(t<1) then
         !c=int(t*1000)
         !write(pretty_time,'(A,I3.3)')  '00.',c
         c=int(t*100)
         write(pretty_time,'(A,I2.2,A)') ' 0.',c,'s'
      elseif(t<60) then
         s=int(t)
         c=int((t-s)*100)
         write(pretty_time,'(I2,A,I2.2,A)') s,'.',c,'s'
      elseif(t<3600) then
         m=int(t/60)
         s=mod(int(t),60)
         write(pretty_time,'(I2,A,I2.2,A)') m,'m',s,'s'
      elseif(t<86400) then
         h=int(t/3600);
         m=int(mod(int(t),3600)/60);
         write(pretty_time,'(I2,A,I2.2,A)') h,'h',m,'m'
      elseif(t<8553600) then
         d=int(t/86400);
         h=int(mod(int(t),86400)/3600);
         write(pretty_time,'(I2,A,I2.2,A)') d,'d',h,'h'
      else
         pretty_time='+3mon.'
      endif
   end function pretty_time
   ! --- END TIC TOC MODULE

end module FVW_VortexTools
