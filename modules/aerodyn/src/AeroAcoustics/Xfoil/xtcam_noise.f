***********************************************************************
C    Module:  xtcam.f
C 
C    Copyright (C) 2000 Harold Youngren, Mark Drela
C 
C    This program is free software; you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation; either version 2 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program; if not, write to the Free Software
C    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
C***********************************************************************


      SUBROUTINE GETCAM (XCM,YCM,NCM,XTK,YTK,NTK,
     &                   X,XP,Y,YP,S,N )
C------------------------------------------------------
C     Finds camber and thickness 
C     distribution for input airfoil 
C------------------------------------------------------
      REAL XCM(*), YCM(*)
      REAL XTK(*), YTK(*)
      REAL X(*),XP(*),Y(*),YP(*),S(*)
C
      CALL XLFIND(SL,X,XP,Y,YP,S,N)
      XL = SEVAL(SL,X,XP,S,N)
      YL = SEVAL(SL,Y,YP,S,N)
C
C---- go over each point, finding opposite points, getting camber and thickness
      DO 10 I=1, N
C------ coordinates of point on the opposite side with the same x value
        CALL SOPPS(SOPP, S(I), X,XP,Y,YP,S,N,SL)
        XOPP = SEVAL(SOPP,X,XP,S,N)
        YOPP = SEVAL(SOPP,Y,YP,S,N)
C
C------ get camber and thickness
        XCM(I) = 0.5*(X(I)+XOPP)
        YCM(I) = 0.5*(Y(I)+YOPP)
        XTK(I) = 0.5*(X(I)+XOPP)
        YTK(I) = 0.5*(Y(I)-YOPP)
        YTK(I) = ABS(YTK(I))
c        if (XOPP.gt.0.9) then
c         write(*,*) 'cm i,x,y ',i,xcm(i),ycm(i)
c         write(*,*) 'tk i,x,y ',i,xtk(i),ytk(i)
c        endif
   10 CONTINUE
C
C---- Tolerance for nominally identical points
      TOL = 1.0E-3 * (S(N)-S(1))
C
C---- Sort the camber points
      NCM = N+1
      XCM(N+1) = XL
      YCM(N+1) = YL
      CALL SORTOL(TOL,NCM,XCM,YCM)
C
C--- Reorigin camber from LE so camberlines start at Y=0  4/24/01 HHY 
C    policy now to generate camber independent of Y-offsets 
      YOF = YCM(1)
      DO I = 1, NCM
        YCM(I) = YCM(I) - YOF
      END DO
C
C---- Sort the thickness points
      NTK = N+1
      XTK(N+1) = XL
      YTK(N+1) = 0.0
      CALL SORTOL(TOL,NTK,XTK,YTK)
C
      RETURN
      END ! GETCAM


      SUBROUTINE GETMAX(X,Y,YP,N,XMAX,YMAX)
      REAL X(*), Y(*), YP(*)
C------------------------------------------------
C     Calculates camber or thickness highpoint 
C     and x position
C------------------------------------------------
C
      XLEN = X(N) - X(1)
      XTOL = XLEN * 1.0E-5
C
      CALL SEGSPL(Y,YP,X,N)
C
C---- get approx max point and rough interval size
      YMAX0 = Y(1)
      XMAX0 = X(1)
      DO 5 I = 2, N
        IF (ABS(Y(I)).GT.ABS(YMAX0)) THEN
          YMAX0 = Y(I)
          XMAX0 = 0.5*(X(I-1) + X(I))
          DDX = 0.5*ABS(X(I+1) - X(I-1))
        ENDIF
 5    CONTINUE
      XMAX = XMAX0
C
C---- do a Newton loop to refine estimate
      DO 10 ITER=1, 10
        YMAX  = SEVAL(XMAX,Y,YP,X,N)
        RES   = DEVAL(XMAX,Y,YP,X,N)
        RESP  = D2VAL(XMAX,Y,YP,X,N)
        IF (ABS(XLEN*RESP) .LT. 1.0E-6) GO TO 20
          DX = -RES/RESP
          DX = SIGN( MIN(0.5*DDX,ABS(DX)) , DX)
          XMAX = XMAX + DX
          IF(ABS(DX) .LT. XTOL) GO TO 20
   10 CONTINUE
C      WRITE(*,*)
C     &  'GETMAX: Newton iteration for max camber/thickness failed.'
      YMAX = YMAX0
      XMAX = XMAX0
C
 20   RETURN
      END ! GETMAX



      SUBROUTINE CPCAM(N,X,Y,DYDX,P,DPDX)
      REAL X(*), Y(*), DYDX(*), P(*), DPDX(*)
C------------------------------------------------------------------
C     Generates y(x) camberline from specified DCp(x) distribution.
C
C     Input:  N       number of points
C             X(.)    x array
C             P(.)    DCp array
C             DPDX(.) dDCp/dx array
C
C     Output: Y(.)    y(x) array
C             DYDX(.) dy/dx array
C------------------------------------------------------------------
C---- 1 / 4 pi
      DATA QOPI / 7.9577471545948E-02 /
C
C---- singular part of camber y(x) due to finite loadings P0,P1 at LE and TE
C-    dYSING/dX has logarithmic singularity at x=X0,X1
      YSING(XT) = QOPI*P1*((XT-X1)*LOG(MAX((X1-XT)/(X1-X0),1.E-6)) - XT)
     &          - QOPI*P0*((XT-X0)*LOG(MAX((XT-X0)/(X1-X0),1.E-6)) - XT)
C
      P0 = P(1)
      P1 = P(N)
C
      X0 = X(1)
      X1 = X(N)
C         
C---- calculate Cauchy integral for y'(x) with removed singularity
      DO I=1, N
        DYDX(I) = 0.0
        J = 1
        IF(I.EQ.J) THEN
         YP1 = DPDX(J)
        ELSE
         YP1 = (P(J) - P(I)) / (X(J) - X(I))
        ENDIF
        DO J=2, N
          IF(I.EQ.J) THEN
           YP2 = DPDX(J)
          ELSE
           YP2 = (P(J) - P(I)) / (X(J) - X(I))
          ENDIF
          DYDX(I) = DYDX(I) + 0.5*(YP1+YP2)*(X(J)-X(J-1))
          YP1 = YP2
        END DO
        DYDX(I) = QOPI*DYDX(I)
C
C------ add on removed part of Cauchy integral, further leaving out the
C-      possible infinities at LE and TE so that y(x) can be safely splined. 
C-      The infinities are analytically integrated, and added on to y(x)
C-      with the statement function YSING.
        IF(I.NE.1) THEN
         DYDX(I) = DYDX(I)
     &           - QOPI*(P(I) - P0)*LOG(X(I) - X0)
        ENDIF
        IF(I.NE.N) THEN
         DYDX(I) = DYDX(I)
     &           + QOPI*(P(I) - P1)*LOG(X1 - X(I))
        ENDIF
      END DO
C
C---- integrate regular part of y'(x) from LE
      Y(1) = 0.
      DO I=2, N
        Y(I) = Y(I-1)
     &       + 0.5*(DYDX(I) + DYDX(I-1))*(X(I) - X(I-1))
      END DO
C
C---- add on singular part
      DO I=1, N
        Y(I) = Y(I) + YSING(X(I))
      END DO
C
C---- add offset and angle of attack to get y(0) = y(1) = 0
      Y0 = Y(1)
      Y1 = Y(N)
      DO I=1, N
        Y(I) = Y(I)
     &       - Y0*(X1  -X(I))/(X1-X0)
     &       - Y1*(X(I)-X0  )/(X1-X0)
      END DO
C
      RETURN
      END ! CPCAM
