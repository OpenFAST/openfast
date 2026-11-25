!*********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026  National Laboratory of the Rockies
!
!    This file is part of ExternalInflow.
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
!*********************************************************************************************************************************
module test_ExtInfw_VelNodes

use, intrinsic :: iso_c_binding, only: c_float
use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_Library
use AeroDyn_Types, only: RotInputType, AD_DestroyRotInputType
use ExternalInflow, only: ExtInfw_NumVelNodes, ExtInfw_SetVelNodePositions

implicit none

private
public :: test_ExtInfw_VelNodes_suite

integer(IntKi), parameter :: nBlades = 3
real(ReKi), parameter     :: HubPos(3)  = [ -0.14_ReKi, 0.0_ReKi, 2.19_ReKi ]    ! rotor apex
real(ReKi), parameter     :: NacPos(3)  = [ 0.05_ReKi, 0.0_ReKi, 2.12_ReKi ]     ! nacelle reference point
real(ReKi), parameter     :: TFinPos(3) = [ 0.90_ReKi, 0.0_ReKi, 2.30_ReKi ]     ! tail fin reference point
real(ReKi), parameter     :: TwrTopZ = 2.0_ReKi                                  ! the tower goes from the origin to this height
real(ReKi), parameter     :: HubRad = 0.09_ReKi
real(ReKi), parameter     :: TipRad = 1.19_ReKi
real(ReKi), parameter     :: tol = 1.0E-6_ReKi                                   ! the positions are passed in single precision

contains

!> Collect all exported unit tests
subroutine test_ExtInfw_VelNodes_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_tower_and_nacelle", test_tower_and_nacelle), &
               new_unittest("test_tail_fin", test_tail_fin), &
               new_unittest("test_no_tower", test_no_tower), &
               new_unittest("test_array_size", test_array_size) &
               ]
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> A rotor with a tower and a nacelle, as AD15 has for a land-based turbine: AD15 reads the wind of the hub, of the nodes of each
!! blade, of the tower nodes and of the nacelle, in this order (AD_CalcWind_Rotor), so the external driver has to get the positions
!! of all of them, displaced.  The tower nodes and the nacelle are where their meshes are, not at the origin.
subroutine test_tower_and_nacelle(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter  :: nBladeNodes = 12
   integer(IntKi), parameter  :: nTowerNodes = 4
   type(RotInputType)         :: u
   real(c_float), allocatable :: px(:), py(:), pz(:)
   real(ReKi)                 :: x(3)
   integer(IntKi)             :: n, iTwr, iNac
   integer(IntKi)             :: ErrStat
   character(ErrMsgLen)       :: ErrMsg

   call create_rotor(nBladeNodes, nTowerNodes, .false., u, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return

   n = ExtInfw_NumVelNodes(u)
   call check(error, n, 1 + nBlades*nBladeNodes + nTowerNodes + 1, 'number of velocity points with a tower and a nacelle'); if (allocated(error)) return

   allocate(px(n), py(n), pz(n))
   px = -999.0_c_float
   py = -999.0_c_float
   pz = -999.0_c_float
   call ExtInfw_SetVelNodePositions(u, px, py, pz, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   call check_points(error, nBladeNodes, nTowerNodes, .false., px, py, pz); if (allocated(error)) return

      ! The points that follow the blades, spelled out: the tower base, the tower top, the nacelle
   iTwr = 1 + nBlades*nBladeNodes + 1
   iNac = 1 + nBlades*nBladeNodes + nTowerNodes + 1
   x = displaced([ 0.0_ReKi, 0.0_ReKi, 0.0_ReKi ])
   call check(error, real(pz(iTwr),ReKi), x(3), more='height of the tower base', thr=tol); if (allocated(error)) return
   x = displaced([ 0.0_ReKi, 0.0_ReKi, TwrTopZ ])
   call check(error, real(pz(iNac-1),ReKi), x(3), more='height of the tower top', thr=tol); if (allocated(error)) return
   x = displaced(NacPos)
   call check(error, real(px(iNac),ReKi), x(1), more='x of the nacelle', thr=tol); if (allocated(error)) return
   call check(error, real(pz(iNac),ReKi), x(3), more='height of the nacelle', thr=tol); if (allocated(error)) return

   call AD_DestroyRotInputType(u, ErrStat, ErrMsg)

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> With a tail fin (TFinAero), AD15 reads one more point after the nacelle
subroutine test_tail_fin(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter  :: nBladeNodes = 7
   integer(IntKi), parameter  :: nTowerNodes = 3
   type(RotInputType)         :: u
   real(c_float), allocatable :: px(:), py(:), pz(:)
   real(ReKi)                 :: x(3)
   integer(IntKi)             :: n
   integer(IntKi)             :: ErrStat
   character(ErrMsgLen)       :: ErrMsg

   call create_rotor(nBladeNodes, nTowerNodes, .true., u, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return

   n = ExtInfw_NumVelNodes(u)
   call check(error, n, 1 + nBlades*nBladeNodes + nTowerNodes + 2, 'number of velocity points with a tower, a nacelle and a tail fin'); if (allocated(error)) return

   allocate(px(n), py(n), pz(n))
   px = -999.0_c_float
   py = -999.0_c_float
   pz = -999.0_c_float
   call ExtInfw_SetVelNodePositions(u, px, py, pz, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   call check_points(error, nBladeNodes, nTowerNodes, .true., px, py, pz); if (allocated(error)) return

   x = displaced(TFinPos)
   call check(error, real(px(n),ReKi), x(1), more='x of the tail fin, the last point', thr=tol); if (allocated(error)) return
   call check(error, real(pz(n),ReKi), x(3), more='height of the tail fin, the last point', thr=tol); if (allocated(error)) return

   call AD_DestroyRotInputType(u, ErrStat, ErrMsg)

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Without tower nodes in AD15 (NumTwrNds = 0) the nacelle follows the last blade
subroutine test_no_tower(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter  :: nBladeNodes = 9
   type(RotInputType)         :: u
   real(c_float), allocatable :: px(:), py(:), pz(:)
   real(ReKi)                 :: x(3)
   integer(IntKi)             :: n
   integer(IntKi)             :: ErrStat
   character(ErrMsgLen)       :: ErrMsg

   call create_rotor(nBladeNodes, 0, .false., u, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return

   n = ExtInfw_NumVelNodes(u)
   call check(error, n, 1 + nBlades*nBladeNodes + 1, 'number of velocity points without a tower'); if (allocated(error)) return

   allocate(px(n), py(n), pz(n))
   px = -999.0_c_float
   py = -999.0_c_float
   pz = -999.0_c_float
   call ExtInfw_SetVelNodePositions(u, px, py, pz, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   call check_points(error, nBladeNodes, 0, .false., px, py, pz); if (allocated(error)) return

   x = displaced(NacPos)
   call check(error, real(px(n),ReKi), x(1), more='x of the nacelle, the last point', thr=tol); if (allocated(error)) return
   call check(error, real(pz(n),ReKi), x(3), more='height of the nacelle, the last point', thr=tol); if (allocated(error)) return

   call AD_DestroyRotInputType(u, ErrStat, ErrMsg)

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Arrays that do not have one entry for each point AD15 reads are refused, and left as they are
subroutine test_array_size(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter  :: nBladeNodes = 5
   integer(IntKi), parameter  :: nTowerNodes = 2
   type(RotInputType)         :: u
   real(c_float), allocatable :: px(:), py(:), pz(:)
   integer(IntKi)             :: n
   integer(IntKi)             :: ErrStat
   character(ErrMsgLen)       :: ErrMsg

   call create_rotor(nBladeNodes, nTowerNodes, .false., u, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return

   n = ExtInfw_NumVelNodes(u)
   allocate(px(n-1), py(n-1), pz(n-1))
   px = -999.0_c_float
   py = -999.0_c_float
   pz = -999.0_c_float
   call ExtInfw_SetVelNodePositions(u, px, py, pz, ErrStat, ErrMsg)
   call check(error, ErrStat, ErrID_Fatal, 'Arrays with one entry too few must give a fatal error'); if (allocated(error)) return
   call check(error, all(px == -999.0_c_float) .and. all(py == -999.0_c_float) .and. all(pz == -999.0_c_float), 'Arrays of the wrong size must not be written to'); if (allocated(error)) return

   call AD_DestroyRotInputType(u, ErrStat, ErrMsg)

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Compares every point with the displaced position of the node AD15 reads the wind of at that index: the hub, the nodes of each
!! blade from root to tip, the tower nodes from the base up, the nacelle, the tail fin
subroutine check_points(error, nBladeNodes, nTowerNodes, with_tail_fin, px, py, pz)
   type(error_type), allocatable, intent(out) :: error
   integer(IntKi),                intent(in ) :: nBladeNodes
   integer(IntKi),                intent(in ) :: nTowerNodes
   logical,                       intent(in ) :: with_tail_fin
   real(c_float),                 intent(in ) :: px(:)
   real(c_float),                 intent(in ) :: py(:)
   real(c_float),                 intent(in ) :: pz(:)

   integer(IntKi)                             :: j, k, n

   n = 1
   call check_point(error, n, displaced(HubPos), 'hub', px, py, pz); if (allocated(error)) return
   do k = 1, nBlades
      do j = 1, nBladeNodes
         n = n + 1
         call check_point(error, n, displaced(blade_node(k, j, nBladeNodes)), 'node '//trim(Num2LStr(j))//' of blade '//trim(Num2LStr(k)), px, py, pz); if (allocated(error)) return
      end do
   end do
   do j = 1, nTowerNodes
      n = n + 1
      call check_point(error, n, displaced(tower_node(j, nTowerNodes)), 'tower node '//trim(Num2LStr(j)), px, py, pz); if (allocated(error)) return
   end do
   n = n + 1
   call check_point(error, n, displaced(NacPos), 'nacelle', px, py, pz); if (allocated(error)) return
   if (with_tail_fin) then
      n = n + 1
      call check_point(error, n, displaced(TFinPos), 'tail fin', px, py, pz); if (allocated(error)) return
   end if
   call check(error, size(px), n, 'number of points checked')
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Compares the point n with the position x
subroutine check_point(error, n, x, label, px, py, pz)
   type(error_type), allocatable, intent(out) :: error
   integer(IntKi),                intent(in ) :: n
   real(ReKi),                    intent(in ) :: x(3)
   character(*),                  intent(in ) :: label
   real(c_float),                 intent(in ) :: px(:)
   real(c_float),                 intent(in ) :: py(:)
   real(c_float),                 intent(in ) :: pz(:)

   call check(error, real(px(n),ReKi), x(1), more='x of point '//trim(Num2LStr(n))//', the '//label, thr=tol); if (allocated(error)) return
   call check(error, real(py(n),ReKi), x(2), more='y of point '//trim(Num2LStr(n))//', the '//label, thr=tol); if (allocated(error)) return
   call check(error, real(pz(n),ReKi), x(3), more='z of point '//trim(Num2LStr(n))//', the '//label, thr=tol); if (allocated(error)) return
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Reference position of node j of n of blade k: a straight blade from the hub radius to the tip, in the plane of the rotor
function blade_node(k, j, n) result(x)
   integer(IntKi), intent(in) :: k
   integer(IntKi), intent(in) :: j
   integer(IntKi), intent(in) :: n
   real(ReKi)                 :: x(3)
   real(ReKi)                 :: r, azimuth

   r = HubRad + (TipRad - HubRad)*real(j-1,ReKi)/real(n-1,ReKi)
   azimuth = TwoPi*real(k-1,ReKi)/real(nBlades,ReKi)
   x = HubPos + r*[ 0.0_ReKi, -sin(azimuth), cos(azimuth) ]
end function

!---------------------------------------------------------------------------------------------------------------------------------
!> Reference position of tower node j of n: a vertical tower from the origin
function tower_node(j, n) result(x)
   integer(IntKi), intent(in) :: j
   integer(IntKi), intent(in) :: n
   real(ReKi)                 :: x(3)

   x = [ 0.0_ReKi, 0.0_ReKi, TwrTopZ*real(j-1,ReKi)/real(n-1,ReKi) ]
end function

!---------------------------------------------------------------------------------------------------------------------------------
!> Displacement of the node with the reference position x: not zero at the origin, different for every node, three components
function displacement(x) result(d)
   real(ReKi), intent(in) :: x(3)
   real(R8Ki)             :: d(3)

   d = [ 0.011_R8Ki + 0.020_R8Ki*real(x(3),R8Ki), -0.007_R8Ki + 0.030_R8Ki*real(x(2),R8Ki), 0.013_R8Ki + 0.004_R8Ki*real(x(1),R8Ki) - 0.002_R8Ki*real(x(3),R8Ki)**2 ]
end function

!---------------------------------------------------------------------------------------------------------------------------------
!> Displaced position of the node with the reference position x
function displaced(x) result(xd)
   real(ReKi), intent(in) :: x(3)
   real(ReKi)             :: xd(3)

   xd = real(real(x,R8Ki) + displacement(x), ReKi)
end function

!---------------------------------------------------------------------------------------------------------------------------------
!> Input meshes of an AeroDyn rotor, as AD15 creates them: the hub, the blades, the tower when it has nodes, the nacelle, and the
!! tail fin when AD15 computes its loads.  Every node is displaced.
subroutine create_rotor(nBladeNodes, nTowerNodes, with_tail_fin, u, ErrStat, ErrMsg)
   integer(IntKi),     intent(in   ) :: nBladeNodes
   integer(IntKi),     intent(in   ) :: nTowerNodes
   logical,            intent(in   ) :: with_tail_fin
   type(RotInputType), intent(inout) :: u
   integer(IntKi),     intent(  out) :: ErrStat
   character(*),       intent(  out) :: ErrMsg

   integer(IntKi)                    :: j, k
   real(ReKi), allocatable           :: xBlade(:,:), xTower(:,:)
   integer(IntKi)                    :: ErrStat2
   character(ErrMsgLen)              :: ErrMsg2
   character(*), parameter           :: RoutineName = 'create_rotor'

   ErrStat = ErrID_None
   ErrMsg  = ''

   call create_point_mesh(HubPos, u%HubMotion, ErrStat2, ErrMsg2); if (Failed()) return

   allocate(u%BladeMotion(nBlades), xBlade(3,nBladeNodes))
   do k = 1, nBlades
      do j = 1, nBladeNodes
         xBlade(:,j) = blade_node(k, j, nBladeNodes)
      end do
      call create_line_mesh(xBlade, u%BladeMotion(k), ErrStat2, ErrMsg2); if (Failed()) return
   end do

   if (nTowerNodes > 0) then
      allocate(xTower(3,nTowerNodes))
      do j = 1, nTowerNodes
         xTower(:,j) = tower_node(j, nTowerNodes)
      end do
      call create_line_mesh(xTower, u%TowerMotion, ErrStat2, ErrMsg2); if (Failed()) return
   end if

   call create_point_mesh(NacPos, u%NacelleMotion, ErrStat2, ErrMsg2); if (Failed()) return

   if (with_tail_fin) then
      call create_point_mesh(TFinPos, u%TFinMotion, ErrStat2, ErrMsg2); if (Failed()) return
   end if

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> A mesh of one displaced point at x
subroutine create_point_mesh(x, M, ErrStat, ErrMsg)
   real(ReKi),     intent(in   ) :: x(3)
   type(MeshType), intent(inout) :: M
   integer(IntKi), intent(  out) :: ErrStat
   character(*),   intent(  out) :: ErrMsg

   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(*), parameter       :: RoutineName = 'create_point_mesh'

   ErrStat = ErrID_None
   ErrMsg  = ''

   call MeshCreate(BlankMesh=M, IOS=COMPONENT_INPUT, Nnodes=1, ErrStat=ErrStat2, ErrMess=ErrMsg2, Orientation=.true., TranslationDisp=.true.); if (Failed()) return
   call MeshPositionNode(M, 1, x, ErrStat2, ErrMsg2); if (Failed()) return
   call MeshConstructElement(M, ELEMENT_POINT, ErrStat2, ErrMsg2, p1=1); if (Failed()) return
   call MeshCommit(M, ErrStat2, ErrMsg2); if (Failed()) return
   M%TranslationDisp(:,1) = displacement(x)

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> A line mesh with displaced nodes at x
subroutine create_line_mesh(x, M, ErrStat, ErrMsg)
   real(ReKi),     intent(in   ) :: x(:,:)
   type(MeshType), intent(inout) :: M
   integer(IntKi), intent(  out) :: ErrStat
   character(*),   intent(  out) :: ErrMsg

   integer(IntKi)                :: j, n
   integer(IntKi)                :: ErrStat2
   character(ErrMsgLen)          :: ErrMsg2
   character(*), parameter       :: RoutineName = 'create_line_mesh'

   ErrStat = ErrID_None
   ErrMsg  = ''

   n = size(x,2)
   call MeshCreate(BlankMesh=M, IOS=COMPONENT_INPUT, Nnodes=n, ErrStat=ErrStat2, ErrMess=ErrMsg2, Orientation=.true., TranslationDisp=.true.); if (Failed()) return
   do j = 1, n
      call MeshPositionNode(M, j, x(:,j), ErrStat2, ErrMsg2); if (Failed()) return
   end do
   do j = 1, n-1
      call MeshConstructElement(M, ELEMENT_LINE2, ErrStat2, ErrMsg2, p1=j, p2=j+1); if (Failed()) return
   end do
   call MeshCommit(M, ErrStat2, ErrMsg2); if (Failed()) return
   do j = 1, n
      M%TranslationDisp(:,j) = displacement(x(:,j))
   end do

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

end module
