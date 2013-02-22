* *******************************************************************
* COPYRIGHT (c) 1980 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Original date 1 December 1993
C       Toolpack tool decs employed.
C       Arg dimensions changed to *.
C 1/4/99 Size of MARK increased to 100.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE KB07AD(COUNT,N,INDEX)
      INTEGER N
      DOUBLE PRECISION COUNT(*)
      INTEGER INDEX(*)
      DOUBLE PRECISION AV,X
      INTEGER I,IF,IFK,IFKA,INT,INTEST,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M,
     +        MLOOP
      INTEGER MARK(100)
      DO 10 I = 1,N
        INDEX(I) = I
   10 CONTINUE
      IF (N.EQ.1) GO TO 200
      IF (N.GE.1) GO TO 30
      WRITE (6,FMT=20)
   20 FORMAT (/,/,/,20X,' ***KB07AD***NO NUMBERS TO BE SORTED ** ',
     + 'RETURN TO CALLING PROGRAM')
      GO TO 200
   30 M = 12
      LA = 2
      IS = 1
      IF = N
      DO 190 MLOOP = 1,N
        IFKA = IF - IS
        IF ((IFKA+1).GT.M) GO TO 70
C********* FINAL SORTING ***
        IS1 = IS + 1
        DO 60 J = IS1,IF
          I = J
   40     IF (COUNT(I-1).LT.COUNT(I)) GO TO 60
          IF (COUNT(I-1).GT.COUNT(I)) GO TO 50
          IF (INDEX(I-1).LT.INDEX(I)) GO TO 60
   50     AV = COUNT(I-1)
          COUNT(I-1) = COUNT(I)
          COUNT(I) = AV
          INT = INDEX(I-1)
          INDEX(I-1) = INDEX(I)
          INDEX(I) = INT
          I = I - 1
          IF (I.GT.IS) GO TO 40
   60   CONTINUE
        LA = LA - 2
        GO TO 170
   70   IY = (IS+IF)/2
        X = COUNT(IY)
        INTEST = INDEX(IY)
        COUNT(IY) = COUNT(IF)
        INDEX(IY) = INDEX(IF)
        K = 1
        IFK = IF
        DO 110 I = IS,IF
          IF (X.GT.COUNT(I)) GO TO 110
          IF (X.LT.COUNT(I)) GO TO 80
          IF (INTEST.GT.INDEX(I)) GO TO 110
   80     IF (I.GE.IFK) GO TO 120
          COUNT(IFK) = COUNT(I)
          INDEX(IFK) = INDEX(I)
          K1 = K
          DO 100 K = K1,IFKA
            IFK = IF - K
            IF (COUNT(IFK).GT.X) GO TO 100
            IF (COUNT(IFK).LT.X) GO TO 90
            IF (INTEST.LE.INDEX(IFK)) GO TO 100
   90       IF (I.GE.IFK) GO TO 130
            COUNT(I) = COUNT(IFK)
            INDEX(I) = INDEX(IFK)
            GO TO 110
  100     CONTINUE
          GO TO 120
  110   CONTINUE
  120   COUNT(IFK) = X
        INDEX(IFK) = INTEST
        IP = IFK
        GO TO 140
  130   COUNT(I) = X
        INDEX(I) = INTEST
        IP = I
  140   IF ((IP-IS).GT. (IF-IP)) GO TO 150
        MARK(LA) = IF
        MARK(LA-1) = IP + 1
        IF = IP - 1
        GO TO 160
  150   MARK(LA) = IP - 1
        MARK(LA-1) = IS
        IS = IP + 1
  160   LNGTH = IF - IS
        IF (LNGTH.LE.0) GO TO 180
        LA = LA + 2
        GO TO 190
  170   IF (LA.LE.0) GO TO 200
  180   IF = MARK(LA)
        IS = MARK(LA-1)
  190 CONTINUE
  200 RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1980 Hyprotech UK
* All rights reserved.
*
* None of the comments in this Copyright notice between the lines
* of asterisks shall be removed or altered in any way.
*
* This Package is intended for compilation without modification,
* so most of the embedded comments have been removed.
*
* ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
* SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
*
* Please note that for an ACADEMIC Licence:
*
* 1. The Packages may only be used for academic research or teaching
*    purposes by the Licensee, and must not be copied by the Licensee for
*    use by any other persons. Use of the Packages in any commercial
*    application shall be subject to prior written agreement between
*    Hyprotech UK Limited and the Licensee on suitable terms and
*    conditions, which will include financial conditions.
* 2. All information on the Package is provided to the Licensee on the
*    understanding that the details thereof are confidential.
* 3. All publications issued by the Licensee that include results obtained
*    with the help of one or more of the Packages shall acknowledge the
*    use of the Packages. The Licensee will notify the Numerical Analysis
*    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
* 4. The Packages may be modified by or on behalf of the Licensee
*    for such use in research applications but at no time shall such
*    Packages or modifications thereof become the property of the
*    Licensee. The Licensee shall make available free of charge to the
*    copyright holder for any purpose all information relating to
*    any modification.
* 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
*    direct or consequential loss or damage whatsoever arising out of
*    the use of Packages by the Licensee.
* *******************************************************************
*
* Original date 1 December 1993
C       Toolpack tool decs employed.
C       Arg dimensions changed to *.
C 1/4/99 Size of MARK increased to 100.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE KB08AD(COUNT,N,INDEX)
      INTEGER N
      DOUBLE PRECISION COUNT(*)
      INTEGER INDEX(*)
      DOUBLE PRECISION AV,X
      INTEGER I,IF,IFK,IFKA,INT,INTEST,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M,
     +        MLOOP
      INTEGER MARK(100)
      DO 10 I = 1,N
        INDEX(I) = I
   10 CONTINUE
      IF (N.EQ.1) GO TO 200
      IF (N.GE.1) GO TO 30
      WRITE (6,FMT=20)
   20 FORMAT (/,/,/,20X,' ***KB08AD***NO NUMBERS TO BE SORTED ** ',
     +       'RETURN TO CALLING PROGRAM')
      GO TO 200
   30 M = 12
      LA = 2
      IS = 1
      IF = N
      DO 190 MLOOP = 1,N
        IFKA = IF - IS
        IF ((IFKA+1).GT.M) GO TO 70
C********* FINAL SORTING ***
        IS1 = IS + 1
        DO 60 J = IS1,IF
          I = J
   40     IF (COUNT(I-1).GT.COUNT(I)) GO TO 60
          IF (COUNT(I-1).LT.COUNT(I)) GO TO 50
          IF (INDEX(I-1).LT.INDEX(I)) GO TO 60
   50     AV = COUNT(I-1)
          COUNT(I-1) = COUNT(I)
          COUNT(I) = AV
          INT = INDEX(I-1)
          INDEX(I-1) = INDEX(I)
          INDEX(I) = INT
          I = I - 1
          IF (I.GT.IS) GO TO 40
   60   CONTINUE
        LA = LA - 2
        GO TO 170
   70   IY = (IS+IF)/2
        X = COUNT(IY)
        INTEST = INDEX(IY)
        COUNT(IY) = COUNT(IF)
        INDEX(IY) = INDEX(IF)
        K = 1
        IFK = IF
        DO 110 I = IS,IF
          IF (X.LT.COUNT(I)) GO TO 110
          IF (X.GT.COUNT(I)) GO TO 80
          IF (INTEST.GT.INDEX(I)) GO TO 110
   80     IF (I.GE.IFK) GO TO 120
          COUNT(IFK) = COUNT(I)
          INDEX(IFK) = INDEX(I)
          K1 = K
          DO 100 K = K1,IFKA
            IFK = IF - K
            IF (COUNT(IFK).LT.X) GO TO 100
            IF (COUNT(IFK).GT.X) GO TO 90
            IF (INTEST.LE.INDEX(IFK)) GO TO 100
   90       IF (I.GE.IFK) GO TO 130
            COUNT(I) = COUNT(IFK)
            INDEX(I) = INDEX(IFK)
            GO TO 110
  100     CONTINUE
          GO TO 120
  110   CONTINUE
  120   COUNT(IFK) = X
        INDEX(IFK) = INTEST
        IP = IFK
        GO TO 140
  130   COUNT(I) = X
        INDEX(I) = INTEST
        IP = I
  140   IF ((IP-IS).GT. (IF-IP)) GO TO 150
        MARK(LA) = IF
        MARK(LA-1) = IP + 1
        IF = IP - 1
        GO TO 160
  150   MARK(LA) = IP - 1
        MARK(LA-1) = IS
        IS = IP + 1
  160   LNGTH = IF - IS
        IF (LNGTH.LE.0) GO TO 180
        LA = LA + 2
        GO TO 190
  170   IF (LA.LE.0) GO TO 200
  180   IF = MARK(LA)
        IS = MARK(LA-1)
  190 CONTINUE
  200 RETURN
      END
* COPYRIGHT (c) 1989 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date 28th feb 2005

C 28th February 2005 Version 1.0.0. New routine to replace za02.

      DOUBLE PRECISION FUNCTION ZA12AD(X)
      DOUBLE PRECISION X

C  Returns the current CPU time in seconds when compiled
C  with any Fortran 95 compiler.

      CALL CPU_TIME(ZA12AD)

      END

