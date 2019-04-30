C***********************************************************************
C    Module:  plutil.f
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
      SUBROUTINE SCALIT(II,Y,YOFF,YSF)
      DIMENSION Y(II)
C-------------------------------------------------------------
C     Y(1:II)  array whose scaling factor is to be determined
C     YOFF     offset of Y array  (Y-YOFF is actually scaled)
C     YSF      Y scaling factor
C-------------------------------------------------------------
C
      AG2 = LOG10(2.0)
      AG5 = LOG10(5.0)
C
      YMAX = ABS(Y(1) - YOFF)
      DO 10 I=2, II
        YMAX = MAX( YMAX , ABS(Y(I)-YOFF) )
   10 CONTINUE
C
      IF(YMAX .EQ. 0.0) YMAX = 1.0E-8
      YLOG = LOG10(YMAX)
C
C---- find log of nearest power of 10 above YMAX
      YLOG1 = AINT(YLOG+100.0) - 99.0
 
C---- find log of nearest 2x(power of 10) above YMAX
      YLOG2 = YLOG1 + AG2
      IF(YLOG2-1.0.GT.YLOG) YLOG2 = YLOG2 - 1.0
C
C---- find log of nearest 5x(power of 10) above YMAX
      YLOG5 = YLOG1 + AG5
      IF(YLOG5-1.0.GT.YLOG) YLOG5 = YLOG5 - 1.0
C
C---- find log of smallest upper bound
      GMIN = MIN( YLOG1 , YLOG2 , YLOG5 )
C
C---- set scaling factor
      YSF = 10.0**(-GMIN)
C
      RETURN
      END