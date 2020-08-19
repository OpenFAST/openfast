module FVW_VortexTools
   ! Contains Typical Tools  for vortex methods

   ! Should be *independent* of the Framework and any derived type 

   ! Only low level functions !

   use NWTC_LIBRARY

   implicit none

   ! Tree parameters
   integer, parameter :: IK1 = selected_int_kind(1) ! to store particle branch number (from 1 to 8)
   integer,parameter :: M0 = 1, M1_1=2, M1_2=3, M1_3=4, M2_11=5, M2_21=6, M2_22=7, M2_31=8, M2_32=9, M2_33=10  ! For moment coefficients
   integer,parameter :: M0_000 = 1
   integer,parameter :: M1_100 = 2
   integer,parameter :: M1_010 = 3
   integer,parameter :: M1_001 = 4

   !> 
   type T_Part
      real(ReKi), dimension(:,:), pointer :: P           =>null() 
      real(ReKi), dimension(:,:), pointer :: Alpha       =>null() 
      real(ReKi), dimension(:),   pointer :: RegParam    =>null() 
      integer(IntKi)                      :: RegFunction =-1
      integer(IntKi)                      :: n =-1
   end type T_Part

   !> The node type is recursive and is used to make a chained-list of nodes for the tree
   type T_Node
      real(ReKi)                         :: radius !< Typical dimension of a cell (max of x,y,z extent)
      real(ReKi),dimension(3)            :: center
      real(ReKi),dimension(3,10)         :: Moments
      integer,dimension(:),pointer       :: iPart=>null()  !< indexes of particles stored in this node
      integer,dimension(:),pointer       :: leaves=>null()  ! NOTE: leaves are introduced to save memory
      type(T_Node),dimension(:), pointer :: branches =>null() 
      integer                            :: nPart = -1  ! Number of particles in branches and leaves of this node
   end type T_Node

   !> The type tree contains some basic data, a chained-list of nodes, and a pointer to the Particle data that were used
   type T_Tree
      type(T_Part)  :: Part            !< Storage for all particles
      integer       :: iStep =-1       !< Time step at which the tree was built
      logical       :: bGrown =.false. !< Is the tree build
      type(T_Node)  :: Root            !< Contains the chained-list of nodes
   end type T_Tree

   interface cut_tree
      module procedure cut_tree_parallel ; ! to switch between parallel and rec easily
   end interface 

contains

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

   subroutine LatticeToSegments(LatticePoints, LatticeGamma, iDepthStart, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC, bShedVorticity, bShedLastVorticity )
      real(Reki), dimension(:,:,:),    intent(in   )  :: LatticePoints  !< Points 3 x nSpan x nDepth
      real(Reki), dimension(:,:),      intent(in   )  :: LatticeGamma   !< GammaPanl  nSpan x nDepth
      integer(IntKi),                  intent(in   )  :: iDepthStart    !< Start index for depth dimension
      real(ReKi), dimension(:,:),      intent(inout)  :: SegPoints      !< 
      integer(IntKi), dimension(:,:),  intent(inout)  :: SegConnct      !< 
      real(ReKi),     dimension(:),    intent(inout)  :: SegGamma       !< 
      integer(IntKi),                  intent(inout)  :: iHeadP         !< Index indicating where to start in SegPoints
      integer(IntKi),                  intent(inout)  :: iHeadC         !< Index indicating where to start in SegConnct
      logical       ,                  intent(in   )  :: bShedVorticity !< Shed vorticity is included if true
      logical       ,                  intent(in   )  :: bShedLastVorticity !< Shed the last vorticity segment if true
      ! Local
      integer(IntKi) :: nSpan, nDepth
      integer(IntKi) :: iSpan, iDepth
      integer(IntKi) :: iHeadP0, iseg1, iseg2, iseg3 ,iseg4  !< Index indicating where to start in SegPoints
      real(ReKi) :: Gamma12
      real(ReKi) :: Gamma41

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
            else
               Gamma12 = LatticeGamma(iSpan,iDepth)-LatticeGamma(iSpan,iDepth-1)
            endif
            if (iSpan==1) then
               Gamma41 = LatticeGamma(iSpan,iDepth)
            else
               Gamma41 = LatticeGamma(iSpan,iDepth)-LatticeGamma(iSpan-1,iDepth)
            endif
            ! Segment 1-2
            if (bShedVorticity) then
               SegConnct(1,iHeadC) = iseg1
               SegConnct(2,iHeadC) = iseg2
               SegConnct(3,iHeadC) = iDepth
               SegConnct(4,iHeadC) = iSpan
               SegGamma (iHeadC  ) = Gamma12
               iHeadC=iHeadC+1
            endif
            ! Segment 1-4
            SegConnct(1,iHeadC) = iseg1
            SegConnct(2,iHeadC) = iseg4
            SegConnct(3,iHeadC) = iDepth
            SegConnct(4,iHeadC) = iSpan
            SegGamma (iHeadC  ) = -Gamma41
            iHeadC=iHeadC+1
            ! Segment 4-3
            if (iDepth==nDepth-1) then
               if ((bShedVorticity) .and. (bShedLastVorticity)) then
                  SegConnct(1,iHeadC) = iseg4
                  SegConnct(2,iHeadC) = iseg3
                  SegConnct(3,iHeadC) = iDepth
                  SegConnct(4,iHeadC) = iSpan
                  SegGamma (iHeadC  ) = - LatticeGamma(iSpan,iDepth)
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
          SegDir  = DP/SegLen                     ! Unit vector along segment direction
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
      print'(A24,3F12.4)',trim(Label),U
   end subroutine

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
   subroutine grow_tree(Tree, PartP, PartAlpha, PartRegFunction, PartRegParam, iStep)
      type(T_Tree),               intent(inout), target :: Tree            !< 
      real(ReKi), dimension(:,:), intent(in   ), target :: PartP           !< 
      real(ReKi), dimension(:,:), intent(in   ), target :: PartAlpha       !< 
      integer(IntKi),             intent(in   )         :: PartRegFunction !< 
      real(ReKi), dimension(:),   intent(in   ), target :: PartRegParam    !< 
      integer(IntKi),             intent(in   )         :: iStep           !< 
      type(T_Node), pointer :: node !< Alias
      type(T_Part), pointer :: Part !< Alias
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
      Tree%Part%P           => PartP
      Tree%Part%Alpha       => PartAlpha
      Tree%Part%RegParam    => PartRegParam
      Tree%Part%RegFunction = PartRegFunction
      Tree%Part%n           = size(PartP,2)

      ! --- Handle special case for root node
      node => Tree%Root
      Part => Tree%Part
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
         node%leaves(1) = Part%n !< index
         node%nPart    = 1
      else
         ! Domain dimensions
         max_x=maxval(Part%P(1,1:Part%n)); max_y=maxval(Part%P(2,1:Part%n)); max_z=maxval(Part%P(3,1:Part%n))
         min_x=minval(Part%P(1,1:Part%n)); min_y=minval(Part%P(2,1:Part%n)); min_z=minval(Part%P(3,1:Part%n))

         ! Init of trunc
         ! Radius taken slightly bigger than domain extent. This radius will be divided by 2 successively
         node%radius = max(abs(max_x-min_x),abs(max_y-min_y),abs(max_z-min_z))*1.001_ReKi
         if(node%radius>1e6) then
             print*,'[Error] Domain extent too large, particle points must be invalid';
             print*, min_x, max_x, min_y, max_y, min_z, max_z
             STOP
         endif
         node%center = (/ (max_x+min_x)/2._ReKi, (max_y+min_y)/2._ReKi, (max_z+min_z)/2._ReKi /)
         node%Moments=0.0_ReKi
         if(associated(node%iPart)) then ; print*,'[Error] Node part allocated'; STOP; endif
         allocate(node%iPart(1:Part%n))
         do i=1,Part%n
            node%iPart(i) = i
         end do
         if(associated(node%branches)) then;  print*,'node branches allocated'; STOP; endif
         if(associated(node%leaves)) then;  print*,'node leaves allocated'; STOP; endif
         node%branches=>null()
         node%leaves=>null()
         node%nPart=Part%n
         ! --- Calling grow function on subbrances
         call grow_tree_parallel(Tree%root, Tree%Part)
!          call grow_tree_rec(Tree%root, Tree%Part)
      endif
      Tree%iStep  = iStep
      Tree%bGrown = .true.
   end subroutine grow_tree

   !> Recursive function to grow/setup a tree. 
   !! Note, needed preliminary calc are done by grow_tree before
   recursive subroutine grow_tree_rec(node, Part)
      type(T_Node), target     :: node !<
      type(T_Part), intent(in) :: Part !<
      integer :: i
      !  Sub Step:
      !   -  compute moments and center for the current node
      !   -  allocate branches and leaves
      call grow_tree_substep(node, Part)
      ! Call grow_tree on branches
      if(associated(node%branches)) then
         do i = 1,size(node%branches)
            call grow_tree_rec(node%branches(i), Part)
         end do
      endif
   end subroutine  grow_tree_rec

   !> Perform a substep of tree growth, growing sub branches from a given node/cell
   !! Parent has already setup node%iPart, indices of the particle in this cell
   !! Steps are:
   !!   - Compute node center (barycenter of vorticity)
   !!   - Compute node moments
   !!   - Distribute particles in each 8 octants. Branches are not created for empty octant
   !!   - Allocate branches and leaves and distribute particles to them
   subroutine grow_tree_substep(node, Part)
      type(T_Node), intent(inout) :: node !< Current node we are growing from
      type(T_Part), intent(in)    :: Part !< All particles info
      integer(IK1) :: iPartOctant                                     !< Index corresponding to which octant the particle falls into
      integer      :: nLeaves, nBranches
      integer      :: iLeaf, iOctant, iBranch
      integer      :: i1,i2,i3,i4,i5,i6,i7,i8
      integer      :: i,j,k
      real(ReKi)   :: wTot, wLoc ! Total and local vorticity strength
      real(ReKi)   :: halfSize ! TODO remove me
      real(ReKi),dimension(3) :: locCenter, DeltaP,PartPos,PartAlpha
      real(ReKi),dimension(3) :: nodeGeomCenter !< Geometric center from division of the domain in powers of 2
      real(ReKi),dimension(3) :: nodeBaryCenter !< Vorticity weighted center
      integer(IK1),dimension(:),allocatable :: PartOctant !< Stores the octant (1-8) where each particle belongs
      integer,dimension(8) :: npart_per_octant !< Number of particle per octant
      integer,dimension(8) :: octant2branches  !< Mapping between 8 octants, to index of non empty branch
      integer,dimension(8) :: octant2leaves    !< Idem for singleton/leaves 
      real(ReKi) :: max_x,max_y,max_z !< for domain dimension
      real(ReKi) :: min_x,min_y,min_z !< for domain dimension
      nodeGeomCenter = node%center ! NOTE: we rely on the fact that our parent has set this to the Geometric value 
      nodeBaryCenter = 0.0_ReKi
      wTot = 0.0_ReKi
      ! --- Barycenter of vorticity of the node 
      do i = 1,node%nPart
         PartPos        = Part%P(:,node%iPart(i))
         PartAlpha      = Part%Alpha(:,node%iPart(i))
         wLoc           = (PartAlpha(1)**2 + PartAlpha(2)**2 + PartAlpha(3)**2)**0.5_ReKi ! Vorticity norm
         nodeBaryCenter = nodeBaryCenter + wLoc*PartPos                                   ! Sum coordinates weighted by vorticity
         wTot           = wTot + wLoc                                                     ! Total vorticity
      end do
      ! There is no vorticity, we make it a empty node and we exit 
      if(EqualRealNos(abs(wTot),0.0_ReKi)) then
         node%nPart=0
         if (associated(node%iPart)) deallocate(node%iPart)
         return ! NOTE: we exit 
      endif
      nodeBaryCenter = nodeBaryCenter/wTot ! barycenter of vorticity
      node%center   = nodeBaryCenter  ! updating

      ! --- Calculation of moments about nodeBaryCenter
      do i = 1,node%nPart
         PartPos   = Part%P    (:,node%iPart(i))
         PartAlpha = Part%Alpha(:,node%iPart(i))
         DeltaP    = PartPos-nodeBaryCenter
         ! Order 0
         node%Moments(1:3,M0_000) = node%Moments(1:3,M0_000) + PartAlpha
         ! 1st order
         node%Moments(1:3,M1_100) = node%Moments(1:3,M1_100) + PartAlpha*DeltaP(1) ! 100
         node%Moments(1:3,M1_010) = node%Moments(1:3,M1_010) + PartAlpha*DeltaP(2) ! 010
         node%Moments(1:3,M1_001) = node%Moments(1:3,M1_001) + PartAlpha*DeltaP(3) ! 001
         ! 2nd order
         do j=1,3
            do k=1,j
               node%Moments(1:3,3+j+k+j/3) = node%Moments(1:3,3+j+k+j/3) + PartAlpha*DeltaP(j)*DeltaP(k)
            end do
         end do
      end do

      ! --- Distributing particles to the 8 octants (based on the geometric center!)
      allocate (PartOctant(1:node%nPart))
      npart_per_octant(1:8)=0
      do i = 1,node%nPart
         PartPos      = Part%P(:,node%iPart(i))
         ! index corresponding to which octant the particle falls into
         iPartOctant = int(1,IK1)
         if (PartPos(1) > nodeGeomCenter(1)) iPartOctant = iPartOctant + int(1,IK1)
         if (PartPos(2) > nodeGeomCenter(2)) iPartOctant = iPartOctant + int(2,IK1)
         if (PartPos(3) > nodeGeomCenter(3)) iPartOctant = iPartOctant + int(4,IK1)
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
      if (associated(node%branches)) then
         print*,'Tree build: error, branches associated'
         STOP
      endif
      if (associated(node%leaves)) then
         print*,'Tree build: error, leaves associated'
         STOP
      end if

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
            locCenter = nodeGeomCenter + 0.5*halfSize*(/ (-1)**(iOctant), (-1)**floor(0.5*real(iOctant-1)+1), (-1)**floor(0.25*real(iOctant-1)+1) /)
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
      i1=0; i2=0; i3=0; i4=0; i5=0; i6=0; i7=0; i8=0;
      do i = 1,node%nPart
         iBranch = octant2branches(PartOctant(i))
         if(iBranch>0) then
            select case(iBranch)
            case(1);i1=i1+1; node%branches(1)%iPart(i1) = node%iPart(i)
            case(2);i2=i2+1; node%branches(2)%iPart(i2) = node%iPart(i)
            case(3);i3=i3+1; node%branches(3)%iPart(i3) = node%iPart(i)
            case(4);i4=i4+1; node%branches(4)%iPart(i4) = node%iPart(i)
            case(5);i5=i5+1; node%branches(5)%iPart(i5) = node%iPart(i)
            case(6);i6=i6+1; node%branches(6)%iPart(i6) = node%iPart(i)
            case(7);i7=i7+1; node%branches(7)%iPart(i7) = node%iPart(i)
            case(8);i8=i8+1; node%branches(8)%iPart(i8) = node%iPart(i)
            end select
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
   end subroutine grow_tree_substep

   !> Grow a tree in "parallel", since recursive calls cannot be parallized easily, we unroll the different layer calls
   !! Note, needed preliminary calc are done by grow_tree before!
   subroutine grow_tree_parallel(Root, Part)
      type(T_Node), intent(inout) :: Root
      type(T_Part), intent(in) :: Part
      integer :: i, nBranches
      integer :: i1
      integer :: i2

      ! ---  Unrolled version of grow_tree for the first node
      !  Sub Step:
      !   -  compute moments and center for the current node
      !   -  allocate branches and leaves
      call grow_tree_substep(Root, Part)

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
                  call grow_tree_substep(Root%branches(i), Part)
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
                     call grow_tree_rec(Root%branches(i1)%branches(i2), Part)
                  endif
               endif
            enddo
            !$OMP end do 
            !Note: We could add more levels 
            !$OMP END PARALLEL
         endif
      endif
   end subroutine grow_tree_parallel



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
      if (associated(node%iPart)) then
         print*,'The tree particles were not properly cleaned'
         STOP
         deallocate(node%iPart)
      end if
   end subroutine cut_substep
   
   !> Cut a tree and all its sub-branches, unrolled to use parallelization for the first 3 levels
   subroutine cut_tree_parallel(Tree)
      type(T_Tree), intent(inout) :: Tree
      integer :: i,i1,i2,nBranches
      ! --- Unlinking particles 
      nullify(Tree%Part%P)
      nullify(Tree%Part%Alpha)
      nullify(Tree%Part%RegParam)
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
      if (associated(Tree%root%branches)) then
         print*,'Tree cut: branches are still allocated'
         STOP
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
         print'(A)'      ,trim(preffix)//':nPart       = '//Num2LStr(node%nPart)
         print'(A,3F12.3)',trim(preffix)//':center      =',node%center
         print'(A,1F12.3)',trim(preffix)//':radius      =',node%radius
         if(associated(node%leaves)) then
            do i = 1,size(node%leaves)
               print'(A)',trim(preffix)//':leaf'//trim(Num2LStr(i))//'='//trim(Num2LStr(node%leaves(i)))
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
   ! --- Velocity computation 
   ! --------------------------------------------------------------------------------
   subroutine ui_tree(Tree, CPs, ioff, icp_beg, icp_end, BranchFactor, DistanceDirect, Uind, ErrStat, ErrMsg)
      use FVW_BiotSavart, only: fourpi_inv, ui_part_nograd_11
      type(T_Tree), target,          intent(inout) :: Tree            !< 
      integer,                       intent(in   ) :: ioff            !< 
      integer,                       intent(in   ) :: icp_beg         !< 
      integer,                       intent(in   ) :: icp_end         !< 
      real(ReKi),                    intent(in   ) :: BranchFactor    !<
      real(ReKi),                    intent(in   ) :: DistanceDirect  !< Distance under which direct evaluation should be done no matter what the tree cell size is
      real(ReKi), dimension(:,:),    intent(in   ) :: CPs             !< Control Points  (3 x nCPs)
      real(ReKi), dimension(:,:),    intent(inout) :: Uind            !< Induced velocity at CPs, with side effects (3 x nCPs)
      integer(IntKi),                intent(  out) :: ErrStat         !< Error status of the operation
      character(*),                  intent(  out) :: ErrMsg          !< Error message if ErrStat /= ErrID_None
      real(ReKi), dimension(3) :: Uind_tmp !< 
      real(ReKi), dimension(3) :: CP       !< Current CP
      integer :: icp, nDirect, nQuad
      type(T_Part), pointer :: Part ! Alias
      Part => Tree%Part
      if(.not. associated(Part%P)) then
         ErrMsg='Ui Part Tree called but tree particles not associated'; ErrStat=ErrID_Fatal; return
      endif
      !$OMP PARALLEL DEFAULT(SHARED)
      !$OMP DO PRIVATE(icp,CP,Uind_tmp,nDirect,nQuad) schedule(runtime)
      do icp=icp_beg,icp_end
         CP = CPs(1:3,icp)
         Uind_tmp(1:3) = 0.0_ReKi
         nDirect =0
         nQuad =0
         call ui_tree_11(Tree%root, CP, Uind_tmp, nDirect, nQuad) !< SIDE EFFECTS
         !print*,'Number of direct calls, and quad calls',nDirect, nQuad
         Uind(1:3,ioff+icp-icp_beg+1) = Uind(1:3,ioff+icp-icp_beg+1) + Uind_tmp(1:3)
      enddo
      !$OMP END DO 
      !$OMP END PARALLEL
   contains
      !> Velocity at one control point from the entire tree
      recursive subroutine ui_tree_11(node, CP, Uind, nDirect, nQuad)
         real(ReKi),dimension(3),intent(inout) :: CP, Uind  !< Velocity at control point, with side effect
         integer, intent(inout) :: nDirect,nQuad
         type(T_Node), intent(inout) :: node
         real(ReKi) :: distDirect, coeff
         real(ReKi),dimension(3) :: DeltaP, phi, Uloc
         real(ReKi) :: x,y,z,mx,my,mz,r
         integer :: i,j,ieqj
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
                  nDirect=nDirect+1
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
                     nDirect=nDirect+1
                     Uind(1:3) = Uind(1:3) + Uloc
                  enddo
               endif
               if(associated(node%branches)) then
                  ! TODO: consider implementing a recursive method for that: direct call on all children
                  do i =1,size(node%branches)
                     call ui_tree_11(node%branches(i), CP, Uind, nDirect, nQuad)
                  end do
               endif

            else 
               ! We are far enough, use branch node quadrupole 
               x=DeltaP(1)
               y=DeltaP(2)
               z=DeltaP(3)
               phi = node%Moments(1:3,M0)/r**3*fourpi_inv

               ! Speed, order 0
               Uloc(1) = phi(2)*z - phi(3)*y
               Uloc(2) = phi(3)*x - phi(1)*z
               Uloc(3) = phi(1)*y - phi(2)*x
               Uind = Uind+Uloc

               ! Speed, order 1
               Uloc=0.0_ReKi
               coeff = 3.0_ReKi*x*fourpi_inv/r**5
               Uloc(1) = Uloc(1)- coeff*(node%Moments(3,M1_1)*y-node%Moments(2,M1_1)*z)
               Uloc(2) = Uloc(2)- coeff*(node%Moments(1,M1_1)*z-node%Moments(3,M1_1)*x)
               Uloc(3) = Uloc(3)- coeff*(node%Moments(2,M1_1)*x-node%Moments(1,M1_1)*y)

               coeff = 3.0_ReKi*y*fourpi_inv/r**5
               Uloc(1) = Uloc(1) - coeff*(node%Moments(3,M1_2)*y-node%Moments(2,M1_2)*z)
               Uloc(2) = Uloc(2) - coeff*(node%Moments(1,M1_2)*z-node%Moments(3,M1_2)*x)
               Uloc(3) = Uloc(3) - coeff*(node%Moments(2,M1_2)*x-node%Moments(1,M1_2)*y)

               coeff = 3.0_ReKi*z*fourpi_inv/r**5
               Uloc(1) = Uloc(1) - coeff*(node%Moments(3,M1_3)*y-node%Moments(2,M1_3)*z)
               Uloc(2) = Uloc(2) - coeff*(node%Moments(1,M1_3)*z-node%Moments(3,M1_3)*x)
               Uloc(3) = Uloc(3) - coeff*(node%Moments(2,M1_3)*x-node%Moments(1,M1_3)*y)

               coeff =   fourpi_inv/r**3
               Uloc(1) = Uloc(1) + coeff*node%Moments(3,M1_2) - coeff*node%Moments(2,M1_3)
               Uloc(2) = Uloc(2) + coeff*node%Moments(1,M1_3) - coeff*node%Moments(3,M1_1)
               Uloc(3) = Uloc(3) + coeff*node%Moments(2,M1_1) - coeff*node%Moments(1,M1_2)
               Uind=Uind+Uloc
               Uloc =0.0_ReKi
               do i =1,3
                  coeff = 1.5_ReKi*fourpi_inv/r**5
                  Uloc(1) = Uloc(1) + coeff * (y*node%Moments(3,5+2*(i/2)+3*(i/3)) - z*node%Moments(2,5+2*(i/2)+3*(i/3)))
                  Uloc(2) = Uloc(2) + coeff * (z*node%Moments(1,5+2*(i/2)+3*(i/3)) - x*node%Moments(3,5+2*(i/2)+3*(i/3)))
                  Uloc(3) = Uloc(3) + coeff * (x*node%Moments(2,5+2*(i/2)+3*(i/3)) - y*node%Moments(1,5+2*(i/2)+3*(i/3)))
                  do j=1,i
                     if (i==j) then
                        ieqj = 1
                     else
                        ieqj = 2 
                     end if
                     coeff = -7.5_ReKi*DeltaP(i)*DeltaP(j)*fourpi_inv/r**7
                     Uloc(1) = Uloc(1) + ieqj * coeff * ( y*node%Moments(3,3+i+j+i/3) - z*node%Moments(2,3+i+j+i/3))
                     Uloc(2) = Uloc(2) + ieqj * coeff * ( z*node%Moments(1,3+i+j+i/3) - x*node%Moments(3,3+i+j+i/3))
                     Uloc(3) = Uloc(3) + ieqj * coeff * ( x*node%Moments(2,3+i+j+i/3) - y*node%Moments(1,3+i+j+i/3))
                  end do
               end do
               coeff = 3.0_ReKi*fourpi_inv/r**5
               Uloc(1) = Uloc(1) + coeff * ( y*node%Moments(3,M2_22) - z*node%Moments(2,M2_33) )
               Uloc(2) = Uloc(2) + coeff * ( z*node%Moments(1,M2_33) - x*node%Moments(3,M2_11) )
               Uloc(3) = Uloc(3) + coeff * ( x*node%Moments(2,M2_11) - y*node%Moments(1,M2_22) )

               coeff = 3.0_ReKi*fourpi_inv/r**5 
               Uloc(1) = Uloc(1) + coeff * (z*node%Moments(3,M2_32) + x*node%Moments(3,M2_21) - x*node%Moments(2,M2_31) - y*node%Moments(2,M2_32))
               Uloc(2) = Uloc(2) + coeff * (x*node%Moments(1,M2_31) + y*node%Moments(1,M2_32) - y*node%Moments(3,M2_21) - z*node%Moments(3,M2_31))
               Uloc(3) = Uloc(3) + coeff * (y*node%Moments(2,M2_21) + z*node%Moments(2,M2_31) - z*node%Moments(1,M2_32) - x*node%Moments(1,M2_21))
               Uind(1:3) = Uind(1:3) + Uloc
               nQuad=nQuad+1
            end if ! Far enough
         end if ! had more than 1 particles
      end subroutine ui_tree_11
   end subroutine  ui_tree

end module FVW_VortexTools
