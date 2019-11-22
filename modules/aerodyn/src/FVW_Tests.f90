module FVW_Tests 

   use NWTC_Library

   use FVW_Types
   use FVW_Subs
   use FVW_VortexTools
   use FVW_Wings
   use FVW_IO
   use FVW_BiotSavart

   implicit none

contains

   !>
   subroutine Test_LatticeToSegment(iStat)
      integer(IntKi), intent(  out)  :: iStat !< Status for test
      ! Local
      integer(IntKi),dimension(:,:), allocatable :: SegConnct !< Segment connectivity
      real(ReKi),    dimension(:,:), allocatable :: SegPoints !< Segment Points
      real(ReKi),    dimension(:)  , allocatable :: SegGamma  !< Segment Circulation
      real(ReKi),    dimension(:),   allocatable :: SegSmooth !< 
      !
      real(ReKi),    dimension(:,:,:), allocatable :: LatticePoints1 !< Lattice Points
      real(ReKi),    dimension(:,:,:), allocatable :: LatticePoints2 !< Lattice Points
      real(ReKi),    dimension(:,:),   allocatable :: LatticeGamma1  !< Lattice Circulation
      real(ReKi),    dimension(:,:),   allocatable :: LatticeGamma2  !< Lattice Circulation
      real(ReKi),    dimension(:,:),   allocatable :: CPs   !< ControlPoints
      real(ReKi),    dimension(:,:),   allocatable :: Uind  !< Induced velocity
      integer(IntKi) :: iHeadC
      integer(IntKi) :: iHeadP
      integer(IntKi) :: i,j,k
      integer(IntKi) :: nP
      integer(IntKi) :: nC
      integer(IntKi) :: nP1, nP2
      integer(IntKi) :: nC1, nC2
      integer(IntKi) :: nDepth, nSpan
      integer(IntKi) :: SmoothModel

      ! --- Creating two lattice
      allocate(LatticePoints1(3,2,2)) 
      allocate(LatticePoints2(3,3,4)) 
      allocate(LatticeGamma1(1,1)) ; 
      allocate(LatticeGamma2(2,3)) ; 
      LatticeGamma1=1
      ! Test shed vorticity
      LatticeGamma2(:,1)=1
      LatticeGamma2(:,2)=2
      LatticeGamma2(:,3)=3
      ! Test trailed vorticity
!       LatticeGamma2(1,:)=1
!       LatticeGamma2(2,:)=2
      CALL MeshMe(LatticePoints1,(/0.,0.,0./))
      CALL MeshMe(LatticePoints2,(/0.,0.,1./))

      CALL WrVTK_Lattice('Points1.vtk',LatticePoints1, LatticeGamma1)
      CALL WrVTK_Lattice('Points2.vtk',LatticePoints2, LatticeGamma2)

      ! --- Convert lattice 1 to segments
      nSpan  = size(LatticePoints1,2)
      nDepth = size(LatticePoints1,3)
      nP1 = nSpan*nDepth
      nC1 = 2*(nSpan*nDepth)-nSpan-nDepth
      allocate(SegConnct(1:2,1:nC1)); SegConnct=-1
      allocate(SegPoints(1:3,1:nP1)); SegPoints=-1
      allocate(SegGamma (1:nC1)    ); SegGamma=-999
      allocate(SegSmooth(1:nC1)    ); SegSmooth=0.0_ReKi

      iHeadP=1
      iHeadC=1
      CALL LatticeToSegments(LatticePoints1, LatticeGamma1, 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
      CALL printall()
      CALL WrVTK_Segments('Points1_seg.vtk', SegPoints, SegConnct, SegGamma) 

      allocate(Uind(1:3,1) ); Uind=0.0_ReKi
      allocate(CPs (1:3,1) ); 
      CPs(1:3,1)=(/1.5,1.5,0./)
      SegSmooth=100.0_ReKi
      SmoothModel=0 ! No smooth
      CALL ui_seg(1, 1, 1, CPs, &
      1, nC1, nC1, nP1, SegPoints, SegConnct, SegGamma,   &
      SmoothModel, SegSmooth, Uind)
   print*,'Uind',Uind

      ! --- Convert lattice 2 to segments
      nSpan  = size(LatticePoints2,2)
      nDepth = size(LatticePoints2,3)
      nP2 = nSpan*nDepth
      nC2 = 2*(nSpan*nDepth)-nSpan-nDepth
      deallocate(SegConnct)
      deallocate(SegPoints)
      deallocate(SegGamma)
      allocate(SegConnct(1:2,1:nC2)); SegConnct=-1
      allocate(SegPoints(1:3,1:nP2)); SegPoints=-1
      allocate(SegGamma (1:nC2)    ); SegGamma=-9999
      iHeadP=1
      iHeadC=1
      CALL LatticeToSegments(LatticePoints2, LatticeGamma2, 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
      CALL printall()
      CALL WrVTK_Segments('Points2_seg.vtk', SegPoints, SegConnct, SegGamma) 

      ! --- Concatenate both
      nP = nP1 + nP2
      nC = nC1 + nC2
      iHeadP=1
      iHeadC=1
      deallocate(SegConnct)
      deallocate(SegPoints)
      deallocate(SegGamma)
      allocate(SegConnct(1:2,1:nC)); SegConnct=-1
      allocate(SegPoints(1:3,1:nP)); SegPoints=-1
      allocate(SegGamma (1:nC)    ); SegGamma=-9999
      CALL LatticeToSegments(LatticePoints1, LatticeGamma1, 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
      CALL LatticeToSegments(LatticePoints2, LatticeGamma2, 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
      CALL printall()
      CALL WrVTK_Segments('PointsBoth_seg.vtk', SegPoints, SegConnct, SegGamma) 


   contains
      subroutine printall()
         print*,'Points'
         do i=1,size(SegPoints,2)
            print*,'i',i,'Coords:', SegPoints(1:3,i)
         enddo
         print*,'Connectivity'
         do i=1,size(SegConnct,2)
            print*,'i',i,'Conn:', SegConnct(1:2,i),'Gam:', SegGamma(i)
         enddo
         print*,'-----------------------------'
      endsubroutine

      subroutine MeshMe(M,offset)
         real(ReKi), dimension(:,:,:), intent(inout) :: M
         real(ReKi), dimension(3)    , intent(in   ):: offset
         do j=1,size(M,3)
            do i=1,size(M,2)
               M(1,i,j)=i + offset(1)
               M(2,i,j)=j + offset(2)
               M(3,i,j)=0 + offset(3)
            enddo
         enddo 
      endsubroutine 
   endsubroutine 

end module FVW_Tests
