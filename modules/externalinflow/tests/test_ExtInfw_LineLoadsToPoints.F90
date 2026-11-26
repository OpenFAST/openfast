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
module test_ExtInfw_LineLoadsToPoints

use, intrinsic :: iso_fortran_env, only: error_unit
use testdrive, only: new_unittest, unittest_type, error_type, check
use NWTC_Library
use ExternalInflow, only: ExtInfw_LineLoadsToPoints

implicit none

private
public :: test_ExtInfw_LineLoadsToPoints_suite

!> Meshes of one blade, as Init_ExtInfw sets them up
type :: BladeMeshes
   type(MeshType)    :: LineMotion     !< AeroDyn line2 mesh with the motions
   type(MeshType)    :: LineLoad       !< its sibling with the loads per unit length
   type(MeshType)    :: PointMotion    !< actuator force points with their motions
   type(MeshType)    :: PointLoad      !< their sibling with the loads
   type(MeshMapType) :: MotionMap      !< mapping of the motions from the line to the force points
end type BladeMeshes

contains

!> Collect all exported unit tests
subroutine test_ExtInfw_LineLoadsToPoints_suite(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
               new_unittest("test_uniform_line_load", test_uniform_line_load), &
               new_unittest("test_total_force_and_moment", test_total_force_and_moment), &
               new_unittest("test_first_moment_of_point_forces", test_first_moment_of_point_forces), &
               new_unittest("test_coincident_points", test_coincident_points), &
               new_unittest("test_no_line_load_mesh", test_no_line_load_mesh) &
               ]
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> A uniform line load q on a straight blade: each force point has to get q times its share of the span (from the mid-point to its
!! neighbour on one side to the mid-point on the other side), whatever the AeroDyn nodes are.  The CFD code divides the point load
!! by that share of the span to get the load per unit length it spreads, so any other distribution is a jagged body force there.
!! The loads of Transfer_Line2_to_Point for the same meshes are printed next to the new ones.
subroutine test_uniform_line_load(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter :: nCases = 4
   integer(IntKi), parameter :: nLine(nCases) = [ 20, 82, 82, 75 ]     ! AeroDyn nodes, clustered toward root and tip
   integer(IntKi), parameter :: nPts(nCases)  = [ 82, 40, 120, 75 ]    ! force points, uniform
   real(ReKi), parameter     :: r0 = 0.1_ReKi
   real(ReKi), parameter     :: r1 = 1.3_ReKi
   real(ReKi), parameter     :: q = 30.0_ReKi
   type(BladeMeshes)         :: B
   real(ReKi), allocatable   :: rLine(:), rPts(:), dr(:), F_old(:,:), M_old(:,:)
   real(ReKi)                :: tol, h, dev_old, dev_new
   integer(IntKi)            :: iCase, nL, nP, i
   integer(IntKi)            :: ErrStat
   character(ErrMsgLen)      :: ErrMsg

   tol = 1.0E3_ReKi*epsilon(1.0_ReKi)*q*(r1 - r0)

   do iCase = 1, nCases
      nL = nLine(iCase)
      nP = nPts(iCase)
      allocate(rLine(nL), rPts(nP), dr(nP), F_old(3,nP), M_old(3,nP))
      call span_nodes(nL, r0, r1, .true., rLine)
      call span_nodes(nP, r0, r1, .false., rPts)
      call create_blade(rLine, rPts, .false., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
      call set_motions(rLine, .false., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
      B%LineLoad%Force(1,:) = q

      h = rPts(2) - rPts(1)
      dr = h
      dr(1) = 0.5_ReKi*h
      dr(nP) = 0.5_ReKi*h

      call mesh_mapping_loads(B, F_old, M_old, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
      call ExtInfw_LineLoadsToPoints(B%LineLoad, B%LineMotion, B%PointMotion, B%MotionMap, B%PointLoad, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return

      dev_old = maxval(abs(F_old(1,:)/(q*dr) - 1.0_ReKi))
      dev_new = maxval(abs(B%PointLoad%Force(1,:)/(q*dr) - 1.0_ReKi))
      write (error_unit, '("#   ",i0," AeroDyn nodes, ",i0," force points, largest deviation of a point load from q*dr: Transfer_Line2_to_Point ",es9.2,", ExtInfw_LineLoadsToPoints ",es9.2)') nL, nP, dev_old, dev_new

      do i = 1, nP
         call check(error, B%PointLoad%Force(1,i), q*dr(i), more='force of point '//trim(Num2LStr(i))//' of '//trim(Num2LStr(nP))//' with '//trim(Num2LStr(nL))//' AeroDyn nodes', thr=tol); if (allocated(error)) return
         call check(error, B%PointLoad%Force(2,i), 0.0_ReKi, more='force across the load, point '//trim(Num2LStr(i)), thr=tol); if (allocated(error)) return
         call check(error, B%PointLoad%Force(3,i), 0.0_ReKi, more='force along the blade, point '//trim(Num2LStr(i)), thr=tol); if (allocated(error)) return
      end do

         ! Moments about the force points: none in the interior (the share of the span is centred on the point); the two end points
         ! carry the load of the half spacing next to them, which acts at a third of the spacing from them
      call check(error, B%PointLoad%Moment(2,1), q*h**2/6.0_ReKi, more='moment of the first point', thr=tol); if (allocated(error)) return
      call check(error, B%PointLoad%Moment(2,nP), -q*h**2/6.0_ReKi, more='moment of the last point', thr=tol); if (allocated(error)) return
      do i = 2, nP-1
         call check(error, B%PointLoad%Moment(2,i), 0.0_ReKi, more='moment of point '//trim(Num2LStr(i)), thr=tol); if (allocated(error)) return
      end do

      call destroy_blade(B)
      deallocate(rLine, rPts, dr, F_old, M_old)
   end do

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> The loads of the force points add up to the integral of the line load, and their moments about a common point to its moment,
!! on a curved and deflected blade with a force and a moment per unit length that vary along the span; with force points from
!! root to tip (as ExtInfw places them) and with force points that cover a part of the blade only.
subroutine test_total_force_and_moment(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter :: nCases = 2
   integer(IntKi), parameter :: nL = 30
   integer(IntKi), parameter :: nPts(nCases) = [ 17, 11 ]
   real(ReKi), parameter     :: r0 = 0.1_ReKi
   real(ReKi), parameter     :: r1 = 1.3_ReKi
   real(ReKi), parameter     :: rFirst(nCases) = [ 0.1_ReKi, 0.34_ReKi ]
   real(ReKi), parameter     :: rLast(nCases)  = [ 1.3_ReKi, 1.12_ReKi ]
   type(BladeMeshes)         :: B
   real(ReKi)                :: rLine(nL)
   real(ReKi), allocatable   :: rPts(:)
   real(R8Ki)                :: Ftot(3), Mtot(3), Fsum(3), Msum(3), xP(3)
   real(ReKi)                :: tol
   integer(IntKi)            :: iCase, nP, i, k
   integer(IntKi)            :: ErrStat
   character(ErrMsgLen)      :: ErrMsg

   tol = 1.0E3_ReKi*epsilon(1.0_ReKi)*10.0_ReKi

   call span_nodes(nL, r0, r1, .true., rLine)

   do iCase = 1, nCases
      nP = nPts(iCase)
      allocate(rPts(nP))
      call span_nodes(nP, rFirst(iCase), rLast(iCase), .false., rPts)
      call create_blade(rLine, rPts, .true., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
      call set_motions(rLine, .true., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
      call set_line_loads(rLine, B)
      call line_load_totals(B, Ftot, Mtot)

      call ExtInfw_LineLoadsToPoints(B%LineLoad, B%LineMotion, B%PointMotion, B%MotionMap, B%PointLoad, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return

      Fsum = 0.0_R8Ki
      Msum = 0.0_R8Ki
      do i = 1, nP
         xP = real(B%PointMotion%Position(:,i),R8Ki) + B%PointMotion%TranslationDisp(:,i)
         Fsum = Fsum + real(B%PointLoad%Force(:,i),R8Ki)
         Msum = Msum + real(B%PointLoad%Moment(:,i),R8Ki) + cross_product(xP, real(B%PointLoad%Force(:,i),R8Ki))
      end do

      do k = 1, 3
         call check(error, real(Fsum(k),ReKi), real(Ftot(k),ReKi), more='total force, component '//trim(Num2LStr(k))//', '//trim(Num2LStr(nP))//' force points', thr=tol); if (allocated(error)) return
         call check(error, real(Msum(k),ReKi), real(Mtot(k),ReKi), more='total moment, component '//trim(Num2LStr(k))//', '//trim(Num2LStr(nP))//' force points', thr=tol); if (allocated(error)) return
      end do

      call destroy_blade(B)
      deallocate(rPts)
   end do

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> On a straight blade with force points from root to tip, the point forces alone (a CFD code that does not read the point moments)
!! give the bending moment of the line load at the root, and the point moments add up to zero.
subroutine test_first_moment_of_point_forces(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter :: nL = 82
   integer(IntKi), parameter :: nP = 40
   real(ReKi), parameter     :: r0 = 0.1_ReKi
   real(ReKi), parameter     :: r1 = 1.3_ReKi
   type(BladeMeshes)         :: B
   real(ReKi)                :: rLine(nL), rPts(nP)
   real(R8Ki)                :: f1, f2, s1, s2, bending_line, bending_points, Msum(3)
   real(ReKi)                :: tol
   integer(IntKi)            :: i, j, k
   integer(IntKi)            :: ErrStat
   character(ErrMsgLen)      :: ErrMsg

   tol = 1.0E3_ReKi*epsilon(1.0_ReKi)*100.0_ReKi

   call span_nodes(nL, r0, r1, .true., rLine)
   call span_nodes(nP, r0, r1, .false., rPts)
   call create_blade(rLine, rPts, .false., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   call set_motions(rLine, .false., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   do j = 1, nL
      B%LineLoad%Force(1,j) = 30.0_ReKi*(1.0_ReKi + 2.0_ReKi*rLine(j) - rLine(j)**3)
   end do

      ! Bending moment at the root of the line load, which is linear on each element (Simpson's rule is exact)
   bending_line = 0.0_R8Ki
   do j = 1, nL-1
      s1 = real(rLine(j),R8Ki) - real(r0,R8Ki)
      s2 = real(rLine(j+1),R8Ki) - real(r0,R8Ki)
      f1 = real(B%LineLoad%Force(1,j),R8Ki)
      f2 = real(B%LineLoad%Force(1,j+1),R8Ki)
      bending_line = bending_line + (s2 - s1)*( s1*f1 + 4.0_R8Ki*0.5_R8Ki*(s1 + s2)*0.5_R8Ki*(f1 + f2) + s2*f2 )/6.0_R8Ki
   end do

   call ExtInfw_LineLoadsToPoints(B%LineLoad, B%LineMotion, B%PointMotion, B%MotionMap, B%PointLoad, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return

   bending_points = 0.0_R8Ki
   Msum = 0.0_R8Ki
   do i = 1, nP
      bending_points = bending_points + (real(rPts(i),R8Ki) - real(r0,R8Ki))*real(B%PointLoad%Force(1,i),R8Ki)
      Msum = Msum + real(B%PointLoad%Moment(:,i),R8Ki)
   end do

   call check(error, real(bending_points,ReKi), real(bending_line,ReKi), more='bending moment at the root from the point forces', thr=tol); if (allocated(error)) return
   do k = 1, 3
      call check(error, real(Msum(k),ReKi), 0.0_ReKi, more='sum of the point moments, component '//trim(Num2LStr(k)), thr=tol); if (allocated(error)) return
   end do

   call destroy_blade(B)

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> With one force point on each AeroDyn node, the loads are those of Transfer_Line2_to_Point (curved and deflected blade, force
!! and moment per unit length that vary along the span).
subroutine test_coincident_points(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter :: nL = 25
   real(ReKi), parameter     :: r0 = 0.1_ReKi
   real(ReKi), parameter     :: r1 = 1.3_ReKi
   type(BladeMeshes)         :: B
   real(ReKi)                :: rLine(nL), F_old(3,nL), M_old(3,nL)
   real(ReKi)                :: tol
   integer(IntKi)            :: i, k
   integer(IntKi)            :: ErrStat
   character(ErrMsgLen)      :: ErrMsg

   tol = 1.0E3_ReKi*epsilon(1.0_ReKi)*10.0_ReKi

   call span_nodes(nL, r0, r1, .true., rLine)
   call create_blade(rLine, rLine, .true., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   call set_motions(rLine, .true., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   call set_line_loads(rLine, B)

   call mesh_mapping_loads(B, F_old, M_old, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   call ExtInfw_LineLoadsToPoints(B%LineLoad, B%LineMotion, B%PointMotion, B%MotionMap, B%PointLoad, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return

   do i = 1, nL
      do k = 1, 3
         call check(error, B%PointLoad%Force(k,i), F_old(k,i), more='force of point '//trim(Num2LStr(i))//', component '//trim(Num2LStr(k)), thr=tol); if (allocated(error)) return
         call check(error, B%PointLoad%Moment(k,i), M_old(k,i), more='moment of point '//trim(Num2LStr(i))//', component '//trim(Num2LStr(k)), thr=tol); if (allocated(error)) return
      end do
   end do

   call destroy_blade(B)

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> AeroDyn does not create the mesh of the tower loads when it computes none (TwrAero off), while the mesh of the tower motions
!! exists: the tower force points then get no load, and no error.
subroutine test_no_line_load_mesh(error)
   type(error_type), allocatable, intent(out) :: error

   integer(IntKi), parameter :: nL = 5
   integer(IntKi), parameter :: nP = 7
   type(BladeMeshes)         :: B
   type(MeshType)            :: NoLoad
   real(ReKi)                :: rLine(nL), rPts(nP)
   integer(IntKi)            :: ErrStat
   character(ErrMsgLen)      :: ErrMsg

   call span_nodes(nL, 0.0_ReKi, 1.5_ReKi, .false., rLine)
   call span_nodes(nP, 0.0_ReKi, 1.5_ReKi, .false., rPts)
   call create_blade(rLine, rPts, .false., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   call set_motions(rLine, .false., B, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return
   B%PointLoad%Force = 1.0_ReKi
   B%PointLoad%Moment = 1.0_ReKi

   call ExtInfw_LineLoadsToPoints(NoLoad, B%LineMotion, B%PointMotion, B%MotionMap, B%PointLoad, ErrStat, ErrMsg); call check(error, ErrStat, ErrID_None, trim(ErrMsg)); if (allocated(error)) return

   call check(error, all(B%PointLoad%Force == 0.0_ReKi), 'The force points must get no force when there is no mesh of line loads'); if (allocated(error)) return
   call check(error, all(B%PointLoad%Moment == 0.0_ReKi), 'The force points must get no moment when there is no mesh of line loads'); if (allocated(error)) return

   call destroy_blade(B)

end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Distances from the rotor axis of n nodes between r0 and r1: clustered toward both ends (cosine spacing) or uniform
subroutine span_nodes(n, r0, r1, clustered, r)
   integer(IntKi), intent(in   ) :: n
   real(ReKi),     intent(in   ) :: r0
   real(ReKi),     intent(in   ) :: r1
   logical,        intent(in   ) :: clustered
   real(ReKi),     intent(  out) :: r(n)
   integer(IntKi)                :: j

   do j = 1, n
      if (clustered) then
         r(j) = r0 + 0.5_ReKi*(r1 - r0)*(1.0_ReKi - cos(Pi*real(j-1,ReKi)/real(n-1,ReKi)))
      else
         r(j) = r0 + (r1 - r0)*real(j-1,ReKi)/real(n-1,ReKi)
      end if
   end do
   r(1) = r0
   r(n) = r1
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Position of the blade axis at the distance r from the rotor axis: with prebend and sweep when curved, along z otherwise
function blade_axis(r, curved) result(x)
   real(ReKi), intent(in) :: r
   logical,    intent(in) :: curved
   real(ReKi)             :: x(3)

   x = [ 0.0_ReKi, 0.0_ReKi, r ]
   if (curved) x(1:2) = [ 0.05_ReKi*r**2, 0.02_ReKi*sin(3.0_ReKi*r) ]
end function

!---------------------------------------------------------------------------------------------------------------------------------
!> Meshes of a blade with AeroDyn nodes at rLine and force points at rPts.  The force points lie on the line elements, as those
!! of ExtInfw (CalcForceActuatorPositionsBlade).
subroutine create_blade(rLine, rPts, curved, B, ErrStat, ErrMsg)
   real(ReKi),        intent(in   ) :: rLine(:)
   real(ReKi),        intent(in   ) :: rPts(:)
   logical,           intent(in   ) :: curved
   type(BladeMeshes), intent(inout) :: B
   integer(IntKi),    intent(  out) :: ErrStat
   character(*),      intent(  out) :: ErrMsg

   integer(IntKi)                   :: nL, nP, i, j, jLower
   real(ReKi), allocatable          :: xLine(:,:)
   real(ReKi)                       :: xP(3), w
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   character(*), parameter          :: RoutineName = 'create_blade'

   ErrStat = ErrID_None
   ErrMsg  = ''

   nL = size(rLine)
   nP = size(rPts)
   allocate(xLine(3,nL))
   do j = 1, nL
      xLine(:,j) = blade_axis(rLine(j), curved)
   end do

      ! AeroDyn line: motions, and the sibling with the loads
   call MeshCreate(BlankMesh=B%LineMotion, IOS=COMPONENT_INPUT, Nnodes=nL, ErrStat=ErrStat2, ErrMess=ErrMsg2, Orientation=.true., TranslationDisp=.true., TranslationVel=.true., RotationVel=.true.); if (Failed()) return
   do j = 1, nL
      call MeshPositionNode(B%LineMotion, j, xLine(:,j), ErrStat2, ErrMsg2); if (Failed()) return
   end do
   do j = 1, nL-1
      call MeshConstructElement(B%LineMotion, ELEMENT_LINE2, ErrStat2, ErrMsg2, p1=j, p2=j+1); if (Failed()) return
   end do
   call MeshCommit(B%LineMotion, ErrStat2, ErrMsg2); if (Failed()) return
   call MeshCopy(SrcMesh=B%LineMotion, DestMesh=B%LineLoad, CtrlCode=MESH_SIBLING, IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.true., Moment=.true.); if (Failed()) return

      ! Force points: motions, and the sibling with the loads
   call MeshCreate(BlankMesh=B%PointMotion, IOS=COMPONENT_INPUT, Nnodes=nP, ErrStat=ErrStat2, ErrMess=ErrMsg2, Orientation=.true., TranslationDisp=.true., TranslationVel=.true., RotationVel=.true.); if (Failed()) return
   do i = 1, nP
      jLower = 1
      do j = 1, nL-2
         if (rLine(j+1) < rPts(i)) jLower = j + 1
      end do
      w = (rPts(i) - rLine(jLower)) / (rLine(jLower+1) - rLine(jLower))
      xP = (1.0_ReKi - w)*xLine(:,jLower) + w*xLine(:,jLower+1)
      call MeshPositionNode(B%PointMotion, i, xP, ErrStat2, ErrMsg2); if (Failed()) return
      call MeshConstructElement(B%PointMotion, ELEMENT_POINT, ErrStat2, ErrMsg2, p1=i); if (Failed()) return
   end do
   call MeshCommit(B%PointMotion, ErrStat2, ErrMsg2); if (Failed()) return
   call MeshCopy(SrcMesh=B%PointMotion, DestMesh=B%PointLoad, CtrlCode=MESH_SIBLING, IOS=COMPONENT_OUTPUT, ErrStat=ErrStat2, ErrMess=ErrMsg2, Force=.true., Moment=.true.); if (Failed()) return

   call MeshMapCreate(B%LineMotion, B%PointMotion, B%MotionMap, ErrStat2, ErrMsg2); if (Failed()) return

   B%LineLoad%Force  = 0.0_ReKi
   B%LineLoad%Moment = 0.0_ReKi

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Displaces the AeroDyn nodes (a deflection that grows with the distance from the rotor axis, when deflected) and moves the
!! force points with them, as SetExtInfwPositions does
subroutine set_motions(rLine, deflected, B, ErrStat, ErrMsg)
   real(ReKi),        intent(in   ) :: rLine(:)
   logical,           intent(in   ) :: deflected
   type(BladeMeshes), intent(inout) :: B
   integer(IntKi),    intent(  out) :: ErrStat
   character(*),      intent(  out) :: ErrMsg

   integer(IntKi)                   :: j
   real(R8Ki)                       :: r

   B%LineMotion%TranslationDisp = 0.0_R8Ki
   if (deflected) then
      do j = 1, size(rLine)
         r = real(rLine(j),R8Ki)
         B%LineMotion%TranslationDisp(:,j) = [ 0.1_R8Ki*r**2, -0.03_R8Ki*r, -0.01_R8Ki*r**2 ]
      end do
   end if

   call Transfer_Line2_to_Point(B%LineMotion, B%PointMotion, B%MotionMap, ErrStat, ErrMsg)
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> A force and a moment per unit length that vary along the span, with all three components
subroutine set_line_loads(rLine, B)
   real(ReKi),        intent(in   ) :: rLine(:)
   type(BladeMeshes), intent(inout) :: B

   integer(IntKi)                   :: j

   do j = 1, size(rLine)
      B%LineLoad%Force(:,j)  = [ 5.0_ReKi*cos(2.0_ReKi*rLine(j)), 3.0_ReKi + 2.0_ReKi*rLine(j)**2, 0.3_ReKi*rLine(j) ]
      B%LineLoad%Moment(:,j) = [ 0.1_ReKi*rLine(j), -0.2_ReKi, 0.4_ReKi*sin(rLine(j)) ]
   end do
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Integral along the line of the line load, and of its moment about the origin.  The loads per unit length and the displaced
!! positions are linear on each element, so Simpson's rule is exact; the length of an element is that of the reference positions.
subroutine line_load_totals(B, Ftot, Mtot)
   type(BladeMeshes), intent(in   ) :: B
   real(R8Ki),        intent(  out) :: Ftot(3)
   real(R8Ki),        intent(  out) :: Mtot(3)

   integer(IntKi)                   :: j
   real(R8Ki)                       :: x1(3), x2(3), f1(3), f2(3), m1(3), m2(3), length

   Ftot = 0.0_R8Ki
   Mtot = 0.0_R8Ki
   do j = 1, B%LineLoad%NNodes-1
      x1 = real(B%LineMotion%Position(:,j),R8Ki)   + B%LineMotion%TranslationDisp(:,j)
      x2 = real(B%LineMotion%Position(:,j+1),R8Ki) + B%LineMotion%TranslationDisp(:,j+1)
      f1 = real(B%LineLoad%Force(:,j),R8Ki)
      f2 = real(B%LineLoad%Force(:,j+1),R8Ki)
      m1 = real(B%LineLoad%Moment(:,j),R8Ki)
      m2 = real(B%LineLoad%Moment(:,j+1),R8Ki)
      length = sqrt( sum( ( real(B%LineMotion%Position(:,j+1),R8Ki) - real(B%LineMotion%Position(:,j),R8Ki) )**2 ) )
      Ftot = Ftot + length*0.5_R8Ki*(f1 + f2)
      Mtot = Mtot + length*( 0.5_R8Ki*(m1 + m2) + ( cross_product(x1, f1) + 4.0_R8Ki*cross_product(0.5_R8Ki*(x1 + x2), 0.5_R8Ki*(f1 + f2)) + cross_product(x2, f2) )/6.0_R8Ki )
   end do
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Loads of the force points from the mesh mapping of the NWTC library (Transfer_Line2_to_Point), which ExtInfw used before
!! ExtInfw_LineLoadsToPoints
subroutine mesh_mapping_loads(B, F, M, ErrStat, ErrMsg)
   type(BladeMeshes), intent(inout) :: B
   real(ReKi),        intent(  out) :: F(:,:)
   real(ReKi),        intent(  out) :: M(:,:)
   integer(IntKi),    intent(  out) :: ErrStat
   character(*),      intent(  out) :: ErrMsg

   type(MeshMapType)                :: LoadMap
   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2
   character(*), parameter          :: RoutineName = 'mesh_mapping_loads'

   ErrStat = ErrID_None
   ErrMsg  = ''

   call MeshMapCreate(B%LineLoad, B%PointLoad, LoadMap, ErrStat2, ErrMsg2); if (Failed()) return
   call Transfer_Line2_to_Point(B%LineLoad, B%PointLoad, LoadMap, ErrStat2, ErrMsg2, B%LineMotion, B%PointMotion); if (Failed()) return
   F = B%PointLoad%Force
   M = B%PointLoad%Moment
   call MeshMapDestroy(LoadMap, ErrStat2, ErrMsg2); if (Failed()) return

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine

!---------------------------------------------------------------------------------------------------------------------------------
!> Frees the meshes of a blade (the siblings with the loads go with the meshes of the motions)
subroutine destroy_blade(B)
   type(BladeMeshes), intent(inout) :: B

   integer(IntKi)                   :: ErrStat2
   character(ErrMsgLen)             :: ErrMsg2

   call MeshMapDestroy(B%MotionMap, ErrStat2, ErrMsg2)
   call MeshDestroy(B%LineMotion, ErrStat2, ErrMsg2)
   call MeshDestroy(B%PointMotion, ErrStat2, ErrMsg2)
end subroutine

end module
