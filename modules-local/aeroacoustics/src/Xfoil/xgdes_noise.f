C***********************************************************************
C    Module:  xgdes.f
C 
C    Copyright (C) 2000 Mark Drela 
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

      SUBROUTINE GDES_Noise(Chord_Scale)
      INCLUDE 'XFOIL.INC'
      CHARACTER*4 COMAND, COMOLD
      LOGICAL LRECALC, LMODPL
      DIMENSION XBOX(2), YBOX(2), XRF(2)
C
      CHARACTER*128 COMARG, ARGOLD
C
      DIMENSION IINPUT(20)
      DIMENSION RINPUT(20)
      LOGICAL ERROR
C

C
C--------------------------------------------------------
C      ELSEIF(COMAND.EQ.'SCAL' .OR.
C     &       COMAND.EQ.'S   '      ) THEN

        XXFAC = Chord_Scale
        YYFAC = Chord_Scale

C
       DO I=1, NB
         XB(I) = XB(I)*XXFAC
         YB(I) = YB(I)*YYFAC
       ENDDO
C
C----- re-order if necessary to maintain counterclockwise ordering 
       IF(XXFAC*YYFAC .LT. 0.0) THEN
         DO I=1, NB/2
           XTMP = XB(I)
           YTMP = YB(I)
           XB(I) = XB(NB-I+1)
           YB(I) = YB(NB-I+1)
           XB(NB-I+1) = XTMP
           YB(NB-I+1) = YTMP
         ENDDO
       ENDIF
C
C----- re-spline new geometry
       CALL SCALC(XB,YB,SB,NB)
       CALL SEGSPL(XB,XBP,SB,NB)
       CALL SEGSPL(YB,YBP,SB,NB)
C
       CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
     &             SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &             EI11BA,EI22BA,APX1BA,APX2BA,
     &             EI11BT,EI22BT,APX1BT,APX2BT,
     &             THICKB,CAMBRB )
C
C       CALL NEWPEN(2)
C       CALL PLTAIR(XB,XBP,YB,YBP,SB,NB, XOFF,XSF,YOFF,YSF,'magenta')
C       CALL PLNEWP('magenta')
       LGEOPL = .FALSE.
C

C--------------------------------------------------------
C      ELSEIF(COMAND.EQ.'EXEC' .OR.
C     &       COMAND.EQ.'X   '      ) THEN
       CALL ABCOPY(.TRUE.)


      END ! GDES_Noise



      SUBROUTINE ABCOPY(LCONF)
      INCLUDE 'XFOIL.INC'
      LOGICAL LCONF
C
      IF(NB.LE.1) THEN
       WRITE(*,*) 'ABCOPY: Buffer airfoil not available.'
       RETURN
      ELSEIF(NB.GT.IQX-5) THEN
       WRITE(*,*) 'Maximum number of panel nodes  : ',IQX-5
       WRITE(*,*) 'Number of buffer airfoil points: ',NB
       WRITE(*,*) 'Current airfoil cannot be set.'
       WRITE(*,*) 'Try executing PANE at Top Level instead.'
       RETURN
      ENDIF
      IF(N.NE.NB) LBLINI = .FALSE.
C
      N = NB
      DO 101 I=1, N
        X(I) = XB(I)
        Y(I) = YB(I)
  101 CONTINUE
      LGSAME = .TRUE.
C
      IF(LBFLAP) THEN
       XOF = XBF
       YOF = YBF
       LFLAP = .TRUE.
      ENDIF
C
C---- strip out doubled points
      I = 1
 102  CONTINUE
      I = I+1
      IF(X(I-1).EQ.X(I) .AND. Y(I-1).EQ.Y(I)) THEN
        DO 104 J=I, N-1
          X(J) = X(J+1)
          Y(J) = Y(J+1)
 104    CONTINUE
        N = N-1
      ENDIF
      IF(I.LT.N) GO TO 102
C
      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)
      CALL NCALC(X,Y,S,N,NX,NY)
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE = SEVAL(SLE,X,XP,S,N)
      YLE = SEVAL(SLE,Y,YP,S,N)
      XTE = 0.5*(X(1)+X(N))
      YTE = 0.5*(Y(1)+Y(N))
      CHORD  = SQRT( (XTE-XLE)**2 + (YTE-YLE)**2 )
      CALL TECALC
      CALL APCALC
C
      LGAMU = .FALSE.
      LQINU = .FALSE.
      LWAKE = .FALSE.
      LQAIJ = .FALSE.
      LADIJ = .FALSE.
      LWDIJ = .FALSE.
      LIPAN = .FALSE.
      LVCONV = .FALSE.
      LSCINI = .FALSE.
CCC      LBLINI = .FALSE.
C
!      IF(LCONF) WRITE(*,1200) N
! 1200 FORMAT(/' Current airfoil nodes set from buffer airfoil nodes (',
!     &        I4,' )')
C
      RETURN
      END ! ABCOPY




      SUBROUTINE SCLXY
C---------------------------------------------------
C     Scale airfoil about LE, TE, or selected point 
C---------------------------------------------------
      INCLUDE 'XFOIL.INC'
      CHARACTER*1 VAR
C
      CALL LEFIND(SBLE,XB,XBP,YB,YBP,SB,NB)
      XLE = SEVAL(SBLE,XB,XBP,SB,NB)
      YLE = SEVAL(SBLE,YB,YBP,SB,NB)
      XTE = 0.5*(XB(1) + XB(NB))
      YTE = 0.5*(YB(1) + YB(NB))
C
      WRITE(*,*) 'Enter origin for airfoil scaling:'
      WRITE(*,*) '  L  scales about LE'
      WRITE(*,*) '  T  scales about TE'
      WRITE(*,*) '  P  scales about input point'
C      
      CALL ASKS('Select origin for scaling^',VAR)
      IF (VAR.EQ.'L') THEN
        XORG = XLE
        YORG = YLE
       ELSE IF (VAR.EQ.'T') THEN
        XORG = XTE
        YORG = YTE
       ELSE 
        XORG = 0.25
        YORG = 0.0
        CALL ASKR('Enter X origin for scaling^',XORG)
        CALL ASKR('Enter Y origin for scaling^',YORG)
      ENDIF       
C      
      SCL = 1.0
      CALL ASKR('Enter scaling factor about selected point^',SCL)
C
      DO 10 I=1, NB
        XB(I) = SCL*(XB(I) - XORG) + XORG 
        YB(I) = SCL*(YB(I) - YORG) + YORG 
   10 CONTINUE
C
      CALL SCALC(XB,YB,SB,NB)
      CALL SEGSPL(XB,XBP,SB,NB)
      CALL SEGSPL(YB,YBP,SB,NB)
C
      CALL GEOPAR(XB,XBP,YB,YBP,SB,NB,W1,
     &            SBLE,CHORDB,AREAB,RADBLE,ANGBTE,
     &            EI11BA,EI22BA,APX1BA,APX2BA,
     &            EI11BT,EI22BT,APX1BT,APX2BT,
     &            THICKB,CAMBRB )
C
      RETURN
      END ! SCLXY





      LOGICAL FUNCTION INSIDE(X,Y,N, XF,YF)
      DIMENSION X(N),Y(N)
C-------------------------------------
C     Returns .TRUE. if point XF,YF 
C     is inside contour X(i),Y(i).
C-------------------------------------
C
C---- integrate subtended angle around airfoil perimeter
      ANGLE = 0.0
      DO 10 I=1, N
        IP = I+1
        IF(I.EQ.N) IP = 1
        XB1 = X(I)  - XF
        YB1 = Y(I)  - YF
        XB2 = X(IP) - XF
        YB2 = Y(IP) - YF
        ANGLE = ANGLE + (XB1*YB2 - YB1*XB2)
     &                   / SQRT((XB1**2 + YB1**2)*(XB2**2 + YB2**2))
 10   CONTINUE
C
C---- angle = 0 if XF,YF is outside, angle = +/- 2 pi  if XF,YF is inside
      INSIDE = ABS(ANGLE) .GT. 1.0
C
      RETURN
      END ! INSIDE



      SUBROUTINE GETXYF(X,XP,Y,YP,S,N, TOPS,BOTS,XF,YF)
      DIMENSION X(N),XP(N),Y(N),YP(N),S(N)
C
      IF(XF .EQ. -999.0)
     &  CALL ASKR('Enter flap hinge x location^',XF)
C
C---- find top and bottom y at hinge x location
      TOPS = S(1) + (X(1) - XF)
      BOTS = S(N) - (X(N) - XF)
      CALL SINVRT(TOPS,XF,X,XP,S,N)      
      CALL SINVRT(BOTS,XF,X,XP,S,N)      
      TOPY = SEVAL(TOPS,Y,YP,S,N)
      BOTY = SEVAL(BOTS,Y,YP,S,N)
C
      WRITE(*,1000) TOPY, BOTY
 1000 FORMAT(/'  Top    surface:  y =', F8.4,'     y/t = 1.0'
     &       /'  Bottom surface:  y =', F8.4,'     y/t = 0.0')
C
      IF(YF .EQ. -999.0)
     & CALL ASKR(
     &  'Enter flap hinge y location (or 999 to specify y/t)^',YF)
C
      IF(YF .EQ. 999.0) THEN
        CALL ASKR('Enter flap hinge relative y/t location^',YREL)
        YF = TOPY*YREL + BOTY*(1.0-YREL)
      ENDIF
C
      RETURN
      END ! GETXYF





