!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026 Alliance for Sustainable Energy, LLC
!
!    This file is part of MoorDyn.
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
!
!**********************************************************************************************************************************
! Slip supports: a line threaded through a hoop on a Point. The hoop acts on a virtual contact point on the
! line; contact forces are spread to nearby nodes with linear (2-node) or cubic B-spline (4-node) weights.
MODULE MoorDyn_SlipSupport

   USE MoorDyn_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: SlipSupport_Setup
   PUBLIC :: SlipSupport_InitLineShape
   PUBLIC :: SlipSupport_CalcForces
   PUBLIC :: SlipSupport_ArmRelease
   PUBLIC :: SlipSupport_CommitRelease

CONTAINS

   !--------------------------------------------------------------
   SUBROUTINE SlipSupport_Setup(Slip, Line, p)

      TYPE(MD_SlipSupport),   INTENT(INOUT) :: Slip
      TYPE(MD_Line),          INTENT(IN   ) :: Line
      TYPE(MD_ParameterType), INTENT(IN   ) :: p

      REAL(DbKi) :: area, w, wWet, mLine

      area = 0.25_DbKi*Pi*Line%d**2
      w    = Line%rho*area
      Slip%mNode = Line%UnstrLen/REAL(Line%N, DbKi) * (w + p%rhoW*Line%Can*area)

      IF (Slip%kPerp <= 0.0_DbKi) Slip%kPerp = Slip%mNode*(0.3_DbKi/p%dtM0)**2

      ! critical damping of whole-line sliding (restoring rate ~2*wWet), capped for explicit stability
      wWet  = ABS(w - p%rhoW*area)*p%g
      mLine = Line%UnstrLen*(w + p%rhoW*Line%Cat*area)
      Slip%cIC = MIN(2.0_DbKi*SQRT(2.0_DbKi*wWet*mLine), 0.5_DbKi*Slip%mNode/p%dtM0)

      Slip%status = 1
      Slip%releasePending = 0
      Slip%kSeg   = 0
      Slip%xi     = 0.0_DbKi
      Slip%s      = 0.0_DbKi
      Slip%vSlip  = 0.0_DbKi
      Slip%Fn     = 0.0_DbKi
      Slip%Ff     = 0.0_DbKi
      Slip%Fhoop  = 0.0_DbKi

   END SUBROUTINE SlipSupport_Setup
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   ! Lay out a line's internal nodes along the polyline A -> hoop(s) -> B, with hoops in input order.
   SUBROUTINE SlipSupport_InitLineShape(m, p, iLine)

      TYPE(MD_MiscVarType),   INTENT(INOUT) :: m
      TYPE(MD_ParameterType), INTENT(IN   ) :: p
      INTEGER(IntKi),         INTENT(IN   ) :: iLine

      INTEGER(IntKi)          :: is, nV, I, J, N
      REAL(DbKi), ALLOCATABLE :: rV(:,:), cV(:)
      REAL(DbKi)              :: sNode, dTarget, frac

      N  = m%LineList(iLine)%N

      nV = 0
      DO is = 1, p%nSlips
         IF (m%SlipList(is)%LineID == iLine) nV = nV + 1
      END DO
      IF (nV == 0) RETURN

      ALLOCATE(rV(3, 0:nV+1), cV(0:nV+1))

      rV(:,0) = m%LineList(iLine)%r(:,0)
      J = 0
      DO is = 1, p%nSlips
         IF (m%SlipList(is)%LineID == iLine) THEN
            J = J + 1
            rV(:,J) = m%PointList(m%SlipList(is)%PointID)%r
         END IF
      END DO
      rV(:,nV+1) = m%LineList(iLine)%r(:,N)

      cV(0) = 0.0_DbKi
      DO J = 1, nV+1
         cV(J) = cV(J-1) + SQRT(SUM((rV(:,J) - rV(:,J-1))**2))
      END DO

      IF (cV(nV+1) <= 0.0_DbKi) THEN
         DEALLOCATE(rV, cV)
         RETURN
      END IF

      ! node arc length is mapped proportionally onto the polyline chord length
      sNode = 0.0_DbKi
      DO I = 1, N-1
         sNode   = sNode + m%LineList(iLine)%l(I)
         dTarget = sNode/m%LineList(iLine)%UnstrLen*cV(nV+1)
         DO J = 1, nV+1
            IF ((dTarget <= cV(J)) .OR. (J == nV+1)) EXIT
         END DO
         frac = (dTarget - cV(J-1))/MAX(cV(J) - cV(J-1), 1.0E-12_DbKi)
         m%LineList(iLine)%r(:,I)  = rV(:,J-1) + frac*(rV(:,J) - rV(:,J-1))
         m%LineList(iLine)%rd(:,I) = 0.0_DbKi
      END DO

      J = 0
      DO is = 1, p%nSlips
         IF (m%SlipList(is)%LineID == iLine) THEN
            J = J + 1
            m%SlipList(is)%s = cV(J)/cV(nV+1)*m%LineList(iLine)%UnstrLen
         END IF
      END DO

      DEALLOCATE(rV, cV)

   END SUBROUTINE SlipSupport_InitLineShape
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   ! Compute contact forces of all slip supports, filling Line%Fslip and Slip%Fhoop.
   ! Must be called after line and point kinematics are set and before Line_GetStateDeriv.
   SUBROUTINE SlipSupport_CalcForces(m, p)

      TYPE(MD_MiscVarType),   INTENT(INOUT) :: m
      TYPE(MD_ParameterType), INTENT(IN   ) :: p

      INTEGER(IntKi) :: is, iL, N, k, kLo, kHi, kBest, j, nW
      REAL(DbKi)     :: xh(3), vh(3), u, d2, d2Best, uBest
      REAL(DbKi)     :: tC(3), tCmag, rStar(3), vStar(3), delta(3), vRel(3), vRelPerp(3)
      REAL(DbKi)     :: Fperp(3), F(3), Fsum(3), mC, cPerp, vTolEff, muEff, overshoot
      INTEGER(IntKi), ALLOCATABLE :: idx(:)
      REAL(DbKi),     ALLOCATABLE :: w(:), dw(:), sNode(:), qN(:,:), Fnode(:,:)

      DO is = 1, p%nSlips
         m%LineList(m%SlipList(is)%LineID)%Fslip = 0.0_DbKi
      END DO

      DO is = 1, p%nSlips

         IF (m%SlipList(is)%status /= 1) CYCLE

         iL = m%SlipList(is)%LineID
         N  = m%LineList(iL)%N
         xh = m%PointList(m%SlipList(is)%PointID)%r
         vh = m%PointList(m%SlipList(is)%PointID)%rd

         IF (ALLOCATED(sNode)) THEN
            IF (SIZE(sNode) /= N+1) DEALLOCATE(idx, w, dw, sNode, qN, Fnode)
         END IF
         IF (.NOT. ALLOCATED(sNode)) ALLOCATE(idx(N+1), w(N+1), dw(N+1), sNode(0:N), qN(3,N+1), Fnode(3,N+1))
         sNode(0) = 0.0_DbKi
         DO j = 1, N
            sNode(j) = sNode(j-1) + m%LineList(iL)%l(j)
         END DO

         ! search neighbors of the cached segment; full search when not yet located
         IF (m%SlipList(is)%kSeg < 1) THEN
            kLo = 1
            kHi = N
         ELSE
            kLo = MAX(1, m%SlipList(is)%kSeg - 2)
            kHi = MIN(N, m%SlipList(is)%kSeg + 2)
         END IF

         d2Best = HUGE(1.0_DbKi)
         kBest  = kLo
         uBest  = 0.0_DbKi
         DO k = kLo, kHi
            CALL ProjectOnSegment(m%LineList(iL), m%SlipList(is)%shape, m%SlipList(is)%width, sNode, k, xh, u, d2, idx, w, dw)
            IF (d2 < d2Best) THEN
               d2Best = d2
               kBest  = k
               uBest  = u
            END IF
         END DO

         k = kBest
         u = uBest
         m%SlipList(is)%kSeg = k
         m%SlipList(is)%xi   = u

         CALL GetWeights(m%LineList(iL), m%SlipList(is)%shape, m%SlipList(is)%width, sNode, k, u, nW, idx, w, dw)

         rStar = 0.0_DbKi
         vStar = 0.0_DbKi
         DO j = 1, nW
            rStar = rStar + w(j)*m%LineList(iL)%r( :,idx(j))
            vStar = vStar + w(j)*m%LineList(iL)%rd(:,idx(j))
         END DO

         IF ((m%SlipList(is)%shape == 2) .OR. (m%SlipList(is)%width > 0.0_DbKi)) THEN
            tC = 0.0_DbKi
            DO j = 1, nW
               tC = tC + dw(j)*m%LineList(iL)%r(:,idx(j))
            END DO
         ELSE
            ! blend bisector tangents of the bounding nodes so the normal bisects kinks
            tC = (1.0_DbKi - u)*NodeTangent(m%LineList(iL), k-1) + u*NodeTangent(m%LineList(iL), k)
         END IF
         tCmag = SQRT(DOT_PRODUCT(tC, tC))
         IF (tCmag > 1.0E-6_DbKi) THEN
            tC = tC/tCmag
         ELSE
            tC = UnitOrZero(m%LineList(iL)%r(:,k) - m%LineList(iL)%r(:,k-1))
            IF (DOT_PRODUCT(tC, tC) <= 0.0_DbKi) CYCLE
         END IF

         ! a release is only committed by SlipSupport_CommitRelease after an accepted time step
         overshoot = 0.0_DbKi
         IF ((k == 1) .AND. (u <= 0.0_DbKi)) overshoot = -DOT_PRODUCT(xh - rStar, tC)
         IF ((k == N) .AND. (u >= 1.0_DbKi)) overshoot =  DOT_PRODUCT(xh - rStar, tC)
         IF (overshoot <= 0.1_DbKi*m%LineList(iL)%l(k)) m%SlipList(is)%releasePending = 0

         delta    = rStar - xh
         delta    = delta - DOT_PRODUCT(delta, tC)*tC
         vRel     = vStar - vh
         vRelPerp = vRel - DOT_PRODUCT(vRel, tC)*tC

         ! inertia felt at the contact point when a force is shared among nodes of mass mNode
         mC    = m%SlipList(is)%mNode/SUM(w(1:nW)**2)
         cPerp = 2.0_DbKi*m%SlipList(is)%zeta*SQRT(m%SlipList(is)%kPerp*mC)
         Fperp = -m%SlipList(is)%kPerp*delta - cPerp*vRelPerp

         ! split the normal force per node using each node's own tangent, so a frictionless hoop
         ! never pushes nodes along the line (avoids tension trapping at sharp bends)
         Fsum = 0.0_DbKi
         DO j = 1, nW
            qN(:,j) = NodeTangent(m%LineList(iL), idx(j))
            IF (DOT_PRODUCT(qN(:,j), qN(:,j)) <= 0.0_DbKi) qN(:,j) = tC
            Fnode(:,j) = w(j)*(Fperp - DOT_PRODUCT(Fperp, qN(:,j))*qN(:,j))
            Fsum = Fsum + Fnode(:,j)
         END DO

         ! normal load actually delivered to the line sets the friction capacity
         m%SlipList(is)%Fn    = SQRT(DOT_PRODUCT(Fsum, Fsum))
         m%SlipList(is)%vSlip = DOT_PRODUCT(vRel, tC)

         IF (m%IC_gen) THEN
            m%SlipList(is)%Ff = -m%SlipList(is)%cIC*m%SlipList(is)%vSlip
         ELSE IF (m%SlipList(is)%muS > 0.0_DbKi) THEN
            ! floor vTol so the peak friction slope muS*Fn/vTol stays within explicit stability
            vTolEff = MAX(m%SlipList(is)%vTol, m%SlipList(is)%muS*m%SlipList(is)%Fn*p%dtM0/mC)
            IF (m%SlipList(is)%muS > m%SlipList(is)%muK) THEN
               muEff = m%SlipList(is)%muK + (m%SlipList(is)%muS - m%SlipList(is)%muK)*EXP(-(m%SlipList(is)%vSlip/m%SlipList(is)%vS)**2)
            ELSE
               muEff = m%SlipList(is)%muK
            END IF
            m%SlipList(is)%Ff = -muEff*m%SlipList(is)%Fn*TANH(m%SlipList(is)%vSlip/vTolEff)
         ELSE
            m%SlipList(is)%Ff = 0.0_DbKi
         END IF

         Fsum = 0.0_DbKi
         DO j = 1, nW
            F = Fnode(:,j) + w(j)*m%SlipList(is)%Ff*qN(:,j)
            m%LineList(iL)%Fslip(:,idx(j)) = m%LineList(iL)%Fslip(:,idx(j)) + F
            Fsum = Fsum + F
         END DO
         m%SlipList(is)%Fhoop = -Fsum

         m%SlipList(is)%s = sNode(k-1) + u*m%LineList(iL)%l(k)

      END DO

   END SUBROUTINE SlipSupport_CalcForces
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   ! Closest contact parameter u in [0,1] on segment k of the contact curve, and squared distance to xh.
   SUBROUTINE ProjectOnSegment(Line, shape, width, sNode, k, xh, u, d2, idx, w, dw)

      TYPE(MD_Line),  INTENT(IN   ) :: Line
      INTEGER(IntKi), INTENT(IN   ) :: shape
      REAL(DbKi),     INTENT(IN   ) :: width
      REAL(DbKi),     INTENT(IN   ) :: sNode(0:)
      INTEGER(IntKi), INTENT(IN   ) :: k
      REAL(DbKi),     INTENT(IN   ) :: xh(3)
      REAL(DbKi),     INTENT(  OUT) :: u
      REAL(DbKi),     INTENT(  OUT) :: d2
      INTEGER(IntKi), INTENT(INOUT) :: idx(:)
      REAL(DbKi),     INTENT(INOUT) :: w(:)
      REAL(DbKi),     INTENT(INOUT) :: dw(:)

      INTEGER(IntKi) :: it, j, nW
      REAL(DbKi)     :: a(3), lvec(3), L2, C(3), Cp(3), h, uNew

      a    = Line%r(:,k-1)
      lvec = Line%r(:,k) - a
      L2   = DOT_PRODUCT(lvec, lvec)
      IF (L2 > 0.0_DbKi) THEN
         u = MIN(1.0_DbKi, MAX(0.0_DbKi, DOT_PRODUCT(xh - a, lvec)/L2))
      ELSE
         u = 0.0_DbKi
      END IF

      IF ((shape /= 2) .AND. (width <= 0.0_DbKi)) THEN
         d2 = SUM((a + u*lvec - xh)**2)
         RETURN
      END IF

      ! Gauss-Newton on the smoothed contact curve, started from the chord projection
      DO it = 1, 10
         CALL GetWeights(Line, shape, width, sNode, k, u, nW, idx, w, dw)
         C  = 0.0_DbKi
         Cp = 0.0_DbKi
         DO j = 1, nW
            C  = C  + w(j) *Line%r(:,idx(j))
            Cp = Cp + dw(j)*Line%r(:,idx(j))
         END DO
         h = DOT_PRODUCT(Cp, Cp)
         IF (h <= 0.0_DbKi) EXIT
         uNew = MIN(1.0_DbKi, MAX(0.0_DbKi, u - DOT_PRODUCT(C - xh, Cp)/h))
         IF (ABS(uNew - u) < 1.0E-10_DbKi) THEN
            u = uNew
            EXIT
         END IF
         u = uNew
      END DO

      CALL GetWeights(Line, shape, width, sNode, k, u, nW, idx, w, dw)
      C = 0.0_DbKi
      DO j = 1, nW
         C = C + w(j)*Line%r(:,idx(j))
      END DO
      d2 = SUM((C - xh)**2)

   END SUBROUTINE ProjectOnSegment
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   SUBROUTINE GetWeights(Line, shape, width, sNode, k, u, nW, idx, w, dw)

      TYPE(MD_Line),  INTENT(IN   ) :: Line
      INTEGER(IntKi), INTENT(IN   ) :: shape
      REAL(DbKi),     INTENT(IN   ) :: width
      REAL(DbKi),     INTENT(IN   ) :: sNode(0:)
      INTEGER(IntKi), INTENT(IN   ) :: k
      REAL(DbKi),     INTENT(IN   ) :: u
      INTEGER(IntKi), INTENT(  OUT) :: nW
      INTEGER(IntKi), INTENT(  OUT) :: idx(:)
      REAL(DbKi),     INTENT(  OUT) :: w(:)
      REAL(DbKi),     INTENT(  OUT) :: dw(:)

      IF (width > 0.0_DbKi) THEN
         CALL KernelWeights(Line, shape, width, sNode, k, u, nW, idx, w, dw)
      ELSE
         CALL ContactWeights(shape, Line%N, k, u, nW, idx, w, dw)
      END IF

   END SUBROUTINE GetWeights
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   ! Weights from a B-spline kernel (shape 1: hat, shape 2: cubic) of fixed scale 'width' in unstretched arc
   ! length, normalized to sum to one. Nodes are mirrored about each line end (r -> 2*r_end - r) and folded
   ! back, so the contact curve passes through the end nodes. dw is d/du on segment k.
   SUBROUTINE KernelWeights(Line, shape, width, sNode, k, u, nW, idx, w, dw)

      TYPE(MD_Line),  INTENT(IN   ) :: Line
      INTEGER(IntKi), INTENT(IN   ) :: shape
      REAL(DbKi),     INTENT(IN   ) :: width
      REAL(DbKi),     INTENT(IN   ) :: sNode(0:)
      INTEGER(IntKi), INTENT(IN   ) :: k
      REAL(DbKi),     INTENT(IN   ) :: u
      INTEGER(IntKi), INTENT(  OUT) :: nW
      INTEGER(IntKi), INTENT(  OUT) :: idx(:)
      REAL(DbKi),     INTENT(  OUT) :: w(:)
      REAL(DbKi),     INTENT(  OUT) :: dw(:)

      INTEGER(IntKi) :: N, j, jLo, jHi
      REAL(DbKi)     :: s, R, Ltot, Asum, dAsum, aK, daK

      N    = Line%N
      Ltot = sNode(N)
      s    = sNode(k-1) + u*Line%l(k)
      IF (shape == 2) THEN
         R = 2.0_DbKi*width
      ELSE
         R = width
      END IF

      jLo = k-1
      DO WHILE (jLo > 0)
         IF (sNode(jLo-1) <= s - R) EXIT
         jLo = jLo - 1
      END DO
      jHi = k
      DO WHILE (jHi < N)
         IF (sNode(jHi+1) >= s + R) EXIT
         jHi = jHi + 1
      END DO

      nW = jHi - jLo + 1
      DO j = jLo, jHi
         idx(j-jLo+1) = j
      END DO
      w(1:nW)  = 0.0_DbKi
      dw(1:nW) = 0.0_DbKi
      Asum  = 0.0_DbKi
      dAsum = 0.0_DbKi

      DO j = jLo, jHi
         CALL Sample(sNode(j), j)
         CALL Add(j, aK, daK)
      END DO

      IF (s - R < 0.0_DbKi) THEN
         DO j = 1, N
            IF (-sNode(j) <= s - R) EXIT
            CALL Sample(-sNode(j), j)
            CALL Add(0,  2.0_DbKi*aK,  2.0_DbKi*daK)
            CALL Add(j, -aK, -daK)
         END DO
      END IF

      IF (s + R > Ltot) THEN
         DO j = 1, N
            IF (2.0_DbKi*Ltot - sNode(N-j) >= s + R) EXIT
            CALL Sample(2.0_DbKi*Ltot - sNode(N-j), N-j)
            CALL Add(N,    2.0_DbKi*aK,  2.0_DbKi*daK)
            CALL Add(N-j, -aK, -daK)
         END DO
      END IF

      ! normalize; Width >= longest segment guarantees Asum > 0
      w(1:nW)  = w(1:nW)/Asum
      dw(1:nW) = (dw(1:nW) - w(1:nW)*dAsum)/Asum * Line%l(k)

   CONTAINS

      ! raw kernel value (a) and its d/ds (da) for a node at arc position sj with tributary length of node jt
      SUBROUTINE Sample(sj, jt)
         REAL(DbKi),     INTENT(IN) :: sj
         INTEGER(IntKi), INTENT(IN) :: jt
         REAL(DbKi)                 :: t, x, Kx, dKx, ds

         t = (sj - s)/width
         x = ABS(t)
         IF (shape == 2) THEN
            IF (x < 1.0_DbKi) THEN
               Kx  = (4.0_DbKi - 6.0_DbKi*x**2 + 3.0_DbKi*x**3)/6.0_DbKi
               dKx = (-12.0_DbKi*x + 9.0_DbKi*x**2)/6.0_DbKi
            ELSE IF (x < 2.0_DbKi) THEN
               Kx  = (2.0_DbKi - x)**3/6.0_DbKi
               dKx = -(2.0_DbKi - x)**2/2.0_DbKi
            ELSE
               Kx  = 0.0_DbKi
               dKx = 0.0_DbKi
            END IF
         ELSE
            IF (x < 1.0_DbKi) THEN
               Kx  = 1.0_DbKi - x
               dKx = -1.0_DbKi
            ELSE
               Kx  = 0.0_DbKi
               dKx = 0.0_DbKi
            END IF
         END IF

         ! tributary length; end nodes count as interior because of mirroring
         IF (jt == 0) THEN
            ds = Line%l(1)
         ELSE IF (jt == N) THEN
            ds = Line%l(N)
         ELSE
            ds = 0.5_DbKi*(Line%l(jt) + Line%l(jt+1))
         END IF

         aK  = Kx*ds
         daK = -SIGN(1.0_DbKi, t)*dKx*ds/width
      END SUBROUTINE Sample

      SUBROUTINE Add(jn, ai, dai)
         INTEGER(IntKi), INTENT(IN) :: jn
         REAL(DbKi),     INTENT(IN) :: ai, dai
         w( jn-jLo+1) = w( jn-jLo+1) + ai
         dw(jn-jLo+1) = dw(jn-jLo+1) + dai
         Asum  = Asum  + ai
         dAsum = dAsum + dai
      END SUBROUTINE Add

   END SUBROUTINE KernelWeights
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   ! Node indices and shape-function weights (and d/du) for contact parameter u on segment k.
   ! shape 1: linear on nodes k-1,k.  shape 2: uniform cubic B-spline on nodes k-2..k+1, with
   ! mirrored phantom nodes at the line ends folded back so the curve passes through end nodes.
   SUBROUTINE ContactWeights(shape, N, k, u, nW, idx, w, dw)

      INTEGER(IntKi), INTENT(IN   ) :: shape
      INTEGER(IntKi), INTENT(IN   ) :: N
      INTEGER(IntKi), INTENT(IN   ) :: k
      REAL(DbKi),     INTENT(IN   ) :: u
      INTEGER(IntKi), INTENT(  OUT) :: nW
      INTEGER(IntKi), INTENT(  OUT) :: idx(:)
      REAL(DbKi),     INTENT(  OUT) :: w(:)
      REAL(DbKi),     INTENT(  OUT) :: dw(:)

      INTEGER(IntKi) :: j, jr
      REAL(DbKi)     :: bs(4), dbs(4)

      nW  = 0
      idx = 0
      w   = 0.0_DbKi
      dw  = 0.0_DbKi

      IF (shape == 2) THEN
         bs  = (/ (1.0_DbKi - u)**3, 3.0_DbKi*u**3 - 6.0_DbKi*u**2 + 4.0_DbKi, &
                  -3.0_DbKi*u**3 + 3.0_DbKi*u**2 + 3.0_DbKi*u + 1.0_DbKi, u**3 /) / 6.0_DbKi
         dbs = (/ -(1.0_DbKi - u)**2, 3.0_DbKi*u**2 - 4.0_DbKi*u, &
                  -3.0_DbKi*u**2 + 2.0_DbKi*u + 1.0_DbKi, u**2 /) / 2.0_DbKi
         DO j = 1, 4
            jr = k - 3 + j
            IF (jr < 0) THEN
               CALL AddWeight(0,    2.0_DbKi*bs(j),  2.0_DbKi*dbs(j))
               CALL AddWeight(1,   -bs(j),          -dbs(j))
            ELSE IF (jr > N) THEN
               CALL AddWeight(N,    2.0_DbKi*bs(j),  2.0_DbKi*dbs(j))
               CALL AddWeight(N-1, -bs(j),          -dbs(j))
            ELSE
               CALL AddWeight(jr,   bs(j),           dbs(j))
            END IF
         END DO
      ELSE
         CALL AddWeight(k-1, 1.0_DbKi - u, -1.0_DbKi)
         CALL AddWeight(k,   u,             1.0_DbKi)
      END IF

   CONTAINS

      SUBROUTINE AddWeight(iNode, wi, dwi)
         INTEGER(IntKi), INTENT(IN) :: iNode
         REAL(DbKi),     INTENT(IN) :: wi, dwi
         INTEGER(IntKi)             :: jj
         DO jj = 1, nW
            IF (idx(jj) == iNode) THEN
               w(jj)  = w(jj)  + wi
               dw(jj) = dw(jj) + dwi
               RETURN
            END IF
         END DO
         nW      = nW + 1
         idx(nW) = iNode
         w(nW)   = wi
         dw(nW)  = dwi
      END SUBROUTINE AddWeight

   END SUBROUTINE ContactWeights
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   ! Arm release checks before a time step; CalcForces clears the flag at any evaluation that is not past a line end.
   SUBROUTINE SlipSupport_ArmRelease(m, p)

      TYPE(MD_MiscVarType),   INTENT(INOUT) :: m
      TYPE(MD_ParameterType), INTENT(IN   ) :: p

      INTEGER(IntKi) :: is

      DO is = 1, p%nSlips
         m%SlipList(is)%releasePending = m%SlipList(is)%status
      END DO

   END SUBROUTINE SlipSupport_ArmRelease
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   ! Permanently release supports whose hoop was past a line end in every force evaluation of the last time step.
   SUBROUTINE SlipSupport_CommitRelease(m, p, t)

      TYPE(MD_MiscVarType),   INTENT(INOUT) :: m
      TYPE(MD_ParameterType), INTENT(IN   ) :: p
      REAL(DbKi),             INTENT(IN   ) :: t

      INTEGER(IntKi)                        :: is
      CHARACTER(200)                        :: Msg

      DO is = 1, p%nSlips
         IF ((m%SlipList(is)%status == 1) .AND. (m%SlipList(is)%releasePending == 1)) THEN
            m%SlipList(is)%status         = 0
            m%SlipList(is)%releasePending = 0
            m%SlipList(is)%vSlip          = 0.0_DbKi
            m%SlipList(is)%Fn             = 0.0_DbKi
            m%SlipList(is)%Ff             = 0.0_DbKi
            m%SlipList(is)%Fhoop          = 0.0_DbKi
            Msg = 'Slip support '//TRIM(Num2LStr(m%SlipList(is)%IdNum))//' released Line '// &
                  TRIM(Num2LStr(m%SlipList(is)%LineID))//' at t = '//TRIM(Num2LStr(t))//' s.'
            CALL WrScr('   '//TRIM(Msg))
            IF (p%writeLog > 0) WRITE(p%UnLog,'(A)') TRIM(Msg)
         END IF
      END DO

   END SUBROUTINE SlipSupport_CommitRelease
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   ! unit bisector of the adjacent segment directions at node i (end nodes use their single segment)
   FUNCTION NodeTangent(Line, i) RESULT(q)

      TYPE(MD_Line),  INTENT(IN) :: Line
      INTEGER(IntKi), INTENT(IN) :: i
      REAL(DbKi)                 :: q(3)

      q = 0.0_DbKi
      IF (i > 0)      q = q + UnitOrZero(Line%r(:,i)   - Line%r(:,i-1))
      IF (i < Line%N) q = q + UnitOrZero(Line%r(:,i+1) - Line%r(:,i))
      q = UnitOrZero(q)

   END FUNCTION NodeTangent
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   FUNCTION UnitOrZero(v) RESULT(u)

      REAL(DbKi), INTENT(IN) :: v(3)
      REAL(DbKi)             :: u(3)
      REAL(DbKi)             :: vMag

      vMag = SQRT(DOT_PRODUCT(v, v))
      IF (vMag > 0.0_DbKi) THEN
         u = v/vMag
      ELSE
         u = 0.0_DbKi
      END IF

   END FUNCTION UnitOrZero
   !--------------------------------------------------------------

END MODULE MoorDyn_SlipSupport
