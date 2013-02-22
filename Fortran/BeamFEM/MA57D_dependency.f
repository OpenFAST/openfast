C *******************************************************************
C COPYRIGHT (c) 1995 Timothy A. Davis, Patrick Amestoy and
C           Council for the Central Laboratory of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
C SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
C
C Please note that for an ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date 30 November 1995
C  April 2001: call to MC49 changed to MC59 to make routine threadsafe
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.
C 23 May 2007 Version 1.1.0. Absolute value of hash taken to cover the
C            case of integer overflow.
C            Comments with character in column 2 corrected.
C 2 August 2007 Version 2.0.0 Dense row handling added, error & warning
C            messages added, iovflo added, interface changed. MC47I/ID
C            added.
C 31 October 2007 Version 2.1.0 Corrected tree formation when handling
C            full variables
C

      SUBROUTINE MC47ID(ICNTL)
      INTEGER ICNTL(10)
      INTEGER I
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      ICNTL(4) = 1
      ICNTL(5) = 2139062143
      DO 100 I=6,10
        ICNTL(I) = 0
 100  CONTINUE
      RETURN
      END
      SUBROUTINE MC47AD(N, NE, PE, IW, IWLEN,
     *      ICNTL,INFO, RINFO)
      INTEGER N, NE, PE(N+1), IWLEN, IW(IWLEN), INFO(10)
      INTEGER ICNTL(10)
      DOUBLE PRECISION RINFO(10)
      INTEGER DEGREE
      DOUBLE PRECISION DUMMY(1)
      INTEGER ELEN,HEAD,I,II,I1,I2,J,LAST,LEN,LENIW,LP,MP,
     *        NEXT,NV,PFREE,W,WP
      INTEGER ICT59(10),INFO59(10),IOUT,JOUT,IDUP,JNFO(10)
      EXTERNAL MC59AD,MC34AD,MC47BD
      DO 5 J = 1,10
         INFO(J) = 0
 5    CONTINUE
      LP = ICNTL(1)
      WP = ICNTL(2)
      MP = ICNTL(3)
      IF (N.LT.1) THEN
        INFO(1) = -1
        IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
     +       '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +       'N has value ',N
        GO TO 1000
      ENDIF
      IF (PE(1).LT.1) THEN
        IF (2*NE+N.GT.IWLEN) THEN
          INFO(1) = -2
         IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
     +        '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +        'IWLEN has value ',IWLEN
          GO TO 1000
        ENDIF
      ELSE
        IF (NE+N.GT.IWLEN) THEN
          INFO(1) = -2
         IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10)')
     +        '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +        'IWLEN has value ',IWLEN
          GO TO 1000
        ENDIF
      ENDIF
      IF (MP.GE.0) THEN
        WRITE(MP,'(/A)') 'Entry to MC47A/AD'
        WRITE(MP,'(A,I10,A,I10,A)') 'Matrix of order',N,' with',NE,
     *                            ' entries'
        IF (PE(1).LT.0)  THEN
          WRITE(MP,'(A)') 'Matrix input in coordinate form'
          WRITE(MP,'(A/(4(I8,I8)))') 'Row and column indices',
     *          (IW(I),IW(NE+I),I=1,NE)
        ELSE
          WRITE(MP,'(A)') 'Matrix input by columns'
          DO 10 J=1,N
            WRITE(MP,'(A,I4/(10I8))') 'Column',J,
     *                                (IW(I),I=PE(J),PE(J+1)-1)
   10     CONTINUE
        ENDIF
      ENDIF
      LAST   = IWLEN  - N + 1
      ELEN   = LAST   - N
      NV     = ELEN   - N
      W      = NV     - N
      DEGREE = W      - N
      HEAD   = DEGREE - N
      NEXT   = HEAD   - N
      LEN    = NEXT   - N
      LENIW = LEN-1
      INFO(6) = 0
      INFO(7) = 0
      IF (PE(1).LT.0) THEN
        DO 20 I=1,NE
          IF (IW(I).LE.IW(NE+I)) THEN
            IF (IW(I).EQ.IW(NE+I) .AND. IW(I).NE.0) THEN
              INFO(7) = INFO(7) + 1
            ELSE
              IF (IW(I).GT.0) INFO(6) = INFO(6) + 1
            ENDIF
            IW(I)=0
          ENDIF
   20   CONTINUE
        ICT59(1) = 0
        ICT59(2) = 1
        ICT59(3) = 1
        ICT59(4) = LP
        ICT59(5) = -1
        ICT59(6) = 0
        CALL MC59AD(ICT59,N,N,NE,IW,NE,IW(NE+1),1,DUMMY,
     *              N+1,PE,N+1,IW(2*NE+1),INFO59)
        IDUP  = INFO59(3)
        IOUT  = INFO59(4)
        JOUT  = INFO59(5)
      ELSE
        IDUP = 0
        IOUT = 0
        JOUT = 0
        DO 30 I = 1,N
          IW(NE+I) = 0
   30   CONTINUE
        DO 50 J=1,N
          I1 = PE(J)
          PE(J) = I1-(IOUT+IDUP)
          I2 = PE(J+1)-1
          IF (I2.LT.I1-1) THEN
            INFO(1) = -3
            GO TO 1000
          ENDIF
          DO 40 II = I1,I2
            I = IW(II)
            IF (I.LE.J .OR. I.GT.N) THEN
              IF (I.EQ.J) INFO(7) = INFO(7) + 1
              IF (I.GT.0 .AND. I.LT.J) INFO(6) = INFO(6) + 1
              IOUT = IOUT + 1
            ELSE
              IF (IW(NE+I).EQ.J) THEN
                IDUP = IDUP + 1
              ELSE
                IW(NE+I)=J
                IW(II-(IOUT+IDUP)) = I
              ENDIF
            ENDIF
   40     CONTINUE
   50   CONTINUE
        PE(N+1) = NE - (IOUT+IDUP) + 1
      ENDIF
      IF (IDUP.GT.0) THEN
        INFO(1) = 1
        INFO(4) = IDUP
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of duplicates found: ',INFO(4)
      ELSE
        INFO(4) = 0
      ENDIF
      IF (IOUT+ JOUT - INFO(7) .GT.0 ) THEN
        INFO(1) = 1
        INFO(5) = IOUT + JOUT - INFO(7)
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of out of range entries found and ignored: ',
     +       INFO(5)
      ELSE
        INFO(5) = 0
      ENDIF
      IF (INFO(6).GT.0) THEN
         INFO(1) = 1
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of entries in upper triangle found and ignored: ',
     +        INFO(6)
      ENDIF
      IF (INFO(7).GT.0) THEN
         INFO(1) = 1
        IF (WP.GE.0) WRITE(WP,'(/A,I3/A,I10)')
     +       '**** Warning from MC47AD **** INFO(1) =',INFO(1),
     +       'Number of entries in diagonals found and ignored: ',
     +        INFO(7)
      ENDIF
      IF (NE-(IOUT+IDUP).EQ.0) THEN
        INFO(1) = -4
        IF (LP.GE.0) WRITE(LP,'(/A,I3/A)')
     +       '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +       'Matrix is null'
        GO TO 1000
      ENDIF
      IF (LENIW.LT.2*(PE(N+1)-1)) THEN
        INFO(1) = -2
         IF (LP.GE.0) WRITE(LP,'(/A,I3/A,I10/A,I10)')
     +        '**** Error return from MC47AD **** INFO(1) =',INFO(1),
     +        'IWLEN has value ',IWLEN,
     +        'Should be at least', 2*(PE(N+1)-1)+8*N
        GO TO 1000
      ENDIF
      CALL MC34AD(N,IW,PE,.FALSE.,DUMMY,IW(W))
      PFREE = PE(N+1)
      DO 60 I=1,N
        IW(LEN+I-1) = PE(I+1) - PE(I)
   60 CONTINUE
      CALL MC47BD(N,LENIW,PE,PFREE,IW(LEN),IW,IW(NV),
     *            IW(ELEN),IW(LAST),IW(DEGREE),
     *            IW(HEAD),IW(NEXT),IW(W), ICNTL,JNFO, RINFO)
      INFO(2) = JNFO(1)
      INFO(3) = PFREE+8*N
      INFO(8) = JNFO(2)
      IF (MP.GE.0) THEN
        WRITE(MP,'(/A)') 'Exit from MC47A/AD'
        WRITE(MP,'(A/(7I10))') 'INFO(1-10):',(INFO(I),I=1,10)
        WRITE(MP,'(A/(8I10))') 'Parent array',(PE(I),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Permutation',(IW(ELEN+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Inverse permutation',
     *                         (IW(LAST+I-1),I=1,N)
        WRITE(MP,'(A/(8I10))') 'Degree array',(IW(NV+I-1),I=1,N)
      ENDIF
 1000 RETURN
      END
      SUBROUTINE MC47BD (N, IWLEN, PE, PFREE, LEN, IW, NV,
     $                   ELEN, LAST, DEGREE,
     $                   HEAD, DENXT, W, ICNTL, JNFO, RJNFO)
      INTEGER N, IWLEN, PE(N), PFREE, LEN(N), IW(IWLEN), NV(N),
     $        ELEN(N), LAST(N),  DEGREE(N),
     $         HEAD(N), DENXT(N), W(N), ICNTL(10), JNFO(10)
      DOUBLE PRECISION RJNFO(10)
      INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, HASH, HMOD, I,
     $     IDUMMY, ILAST, INEXT, IOVFLO,J, JDUMMY, JLAST, JNEXT, K,
     $     KNT1, KNT2, KNT3, LASTD,  LENJ, LN, MAXMEM, ME,
     $     MEM, MINDEG, NBD, NCMPA, NDME, NEL, NELME, NEWMEM,
     $     NFULL, NLEFT, NRLADU, NVI, NVJ, NVPIV, P, P1, P2, P3, PDST,
     $     PEE, PEE1, PEND, PJ, PME, PME1, PME2, PN, PSRC, RSTRT,
     $     SLENME, THRESH, THRESM, WE, WFLG, WNVI,X
     $
      DOUBLE PRECISION RELDEN, SM, STD, OPS
      LOGICAL IDENSE
      INTRINSIC MAX, MIN, MOD
      DO 2 I = 1,10
         RJNFO(I) = 0.0
         JNFO(I) = 0
 2    CONTINUE
      DMAX = 0
      HMOD = MAX (1, N-1)
      IOVFLO = ICNTL(5)
      LASTD = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 1
      NBD   = 0
      NCMPA = 0
      NEL = 0
      NFULL  = 0
      NRLADU = 0
      RSTRT = 0
      OPS = 0.00
      THRESH = ICNTL(4)
      WFLG = 2
      IF (THRESH.GT.0) THEN
         THRESM  = 0
         RELDEN = 0.0
         SM = 0
         DO 5 I=1,N
            THRESM = MAX(THRESM, LEN(I))
            IF (LEN(I).GT.0) THEN
               RELDEN = RELDEN + LEN(I)
               SM = SM + (LEN(I) * LEN(I))
            END IF
            LAST (I) = 0
            HEAD (I) = 0
            NV (I) = 1
            DEGREE (I) = LEN (I)
            IF (DEGREE(I) .EQ. 0) THEN
               NEL = NEL + 1
               ELEN (I) = -NEL
               PE (I) = 0
               W (I) = 0
               NRLADU = NRLADU + 1
               OPS = OPS + 1
            ELSE
               W (I) = 1
               ELEN (I) = 0
            ENDIF
 5       CONTINUE
         IF (N .EQ. NEL) GOTO 265
         RELDEN = RELDEN/(N-NEL)
         SM = SM/(N-NEL-NFULL) - RELDEN*RELDEN
         STD = SQRT(ABS(SM))
         IF (STD .LE. RELDEN) THEN
            THRESM = -1
         ELSE
            THRESM = INT(9*RELDEN + 0.5*STD*((STD/(RELDEN+0.01))**1.5)+
     *           2*RELDEN*RELDEN/(STD+0.01) +1)
         END IF
      ELSE
         THRESM = THRESH
         DO 10 I = 1, N
            LAST (I) = 0
            HEAD (I) = 0
            NV (I) = 1
            DEGREE (I) = LEN (I)
            IF (DEGREE(I) .EQ. 0) THEN
               NEL = NEL + 1
               ELEN (I) = -NEL
               PE (I) = 0
               W (I) = 0
               NRLADU = NRLADU + 1
               OPS = OPS + 1
            ELSE
               W (I) = 1
               ELEN (I) = 0
            ENDIF
 10      CONTINUE
      ENDIF
      IF (THRESM.GE.0) THEN
         IF (THRESM.GE.N) THEN
            THRESM = -1
         ELSE IF (THRESM.EQ.0) THEN
            THRESM = N
         ENDIF
      ENDIF
      DO 20 I = 1, N
         DEG = DEGREE (I)
         IF (DEG .GT. 0) THEN
            IF ( (THRESM.GE.0) .AND.
     &           (DEG+1.GE.THRESM.OR.DEG+1.GE.N-NEL )) THEN
               NBD = NBD+1
               IF (DEG+1.NE.N-NEL) THEN
                  DEGREE(I) = DEGREE(I)+N+1
                  DEG = N
                  INEXT = HEAD (DEG)
                  IF (INEXT .NE. 0) LAST (INEXT) = I
                  DENXT (I) = INEXT
                  HEAD (DEG) = I
                  LAST(I)  = 0
                  IF (LASTD.EQ.0) THEN
                     LASTD=I
                  END IF
               ELSE
                  NFULL = NFULL+1
                  DEGREE(I) = N+1
                  DEG = N
                  IF (LASTD.EQ.0) THEN
                     LASTD     = I
                     HEAD(DEG) = I
                     DENXT(I)   = 0
                     LAST(I)   = 0
                  ELSE
                        DENXT(LASTD) = I
                        LAST(I)     = LASTD
                        LASTD       = I
                        DENXT(I)     = 0
                  ENDIF
               ENDIF
            ELSE
               INEXT = HEAD (DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               DENXT (I) = INEXT
               HEAD (DEG) = I
            ENDIF
         ENDIF
 20   CONTINUE
      IF (NBD.EQ.0 .AND. THRESH.GT.0) THEN
         THRESM = -1
      END IF
 30   IF (NEL .LT. N) THEN
         DO 40 DEG = MINDEG, N
            ME = HEAD (DEG)
            IF (ME .GT. 0) GO TO 50
 40      CONTINUE
 50      MINDEG = DEG
         IF (DEG.LT.N)  THEN
            INEXT = DENXT (ME)
            IF (INEXT .NE. 0) LAST (INEXT) = 0
            HEAD (DEG) = INEXT
         ELSE
            IF (DEGREE(ME).EQ.N+1) GO TO 263
            RSTRT = RSTRT + 1
            RELDEN = 0.0
            SM = 0
            IF (WFLG .GT. IOVFLO-NBD-1) THEN
               DO  51 X = 1, N
                  IF (W (X) .NE. 0) W (X) = 1
 51            CONTINUE
               WFLG = 2
            END IF
            WFLG = WFLG + 1
            DO 57 IDUMMY = 1,N
               INEXT = DENXT (ME)
               IF (INEXT .NE. 0) THEN
                  LAST (INEXT) = 0
               ELSE
                  LASTD = 0
               ENDIF
               DENXT(ME) = 0
               W(ME)      = WFLG
               P1 = PE(ME)
               P2 = P1 + LEN(ME) -1
               LN       = P1
               ELN      = P1
               DO 55 P=P1,P2
                  E= IW(P)
                  IF (W(E).EQ.WFLG) GO TO 55
                  W(E) = WFLG
                  DO 52 JDUMMY = 1,N
                     IF ( PE(E) .GE. 0 ) GOTO 53
                     E = -PE(E)
                     IF (W(E) .EQ.WFLG) GOTO 55
                     W(E) = WFLG
 52               CONTINUE
 53               IF (ELEN(E).LT.0) THEN
                     DENXT(E) = DENXT(E) - NV(ME)
                     IW(LN) = IW(ELN)
                     IW(ELN) = E
                     LN  = LN+1
                     ELN = ELN + 1
                     PEE1 = PE(E)
                     DO 54 PEE = PEE1, PEE1+LEN(E)-1
                        X = IW(PEE)
                        IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                           DENXT(ME) = DENXT(ME) + NV(X)
                           W(X) = WFLG
                        ENDIF
 54                  CONTINUE
                  ELSE
                     DENXT(ME) = DENXT(ME) + NV(E)
                     IW(LN)=E
                     LN = LN+1
                  ENDIF
 55            CONTINUE
               WFLG     = WFLG + 1
               LEN(ME)  = LN-P1
               ELEN(ME) = ELN- P1
               NDME = DENXT(ME)+NV(ME)
               IF (DENXT(ME).EQ.0) DENXT(ME) =1
               IF (NDME .LT. NBD) THEN
                  RELDEN = RELDEN + NV(ME)*NDME
                  SM = SM + NV(ME)*NDME*NDME
                  DEGREE(ME) = DENXT(ME)
                  DEG = DEGREE(ME)
                  MINDEG = MIN(DEG,MINDEG)
                  JNEXT = HEAD(DEG)
                  IF (JNEXT.NE. 0) LAST (JNEXT) = ME
                  DENXT(ME) = JNEXT
                  HEAD(DEG) = ME
               ELSE
                  DEGREE(ME) = N+1
                  DEG = DENXT(ME)
                  MINDEG = MIN(DEG,MINDEG)
                  DEG = N
                  P1 = PE(ME)
                  P2 = P1 + ELEN(ME) - 1
                  DO 56 PJ=P1,P2
                     E= IW(PJ)
                     DENXT (E) = DENXT(E) + NV(ME)
 56               CONTINUE
                  DEG = N
                  NFULL = NFULL +NV(ME)
                  IF (LASTD.EQ.0) THEN
                     LASTD     = ME
                     HEAD(N) = ME
                     DENXT(ME)   = 0
                     LAST(ME)   = 0
                     IF (INEXT.EQ.0) INEXT = LASTD
                  ELSE
                        DENXT(LASTD) = ME
                        LAST(ME)     = LASTD
                        LASTD        = ME
                        DENXT(ME)     = 0
                        IF (INEXT.EQ.0) INEXT = LASTD
                  ENDIF
               END IF
               ME    = INEXT
               IF (ME.EQ.0) GO TO 58
               IF (DEGREE(ME).LE.(N+1) ) GOTO 58
 57         CONTINUE
 58         HEAD (N) = ME
            IF (NBD.EQ.NFULL) THEN
               RELDEN = 0
               SM = 0
            ELSE
               RELDEN = (RELDEN + NFULL*NBD)/(NBD)
               SM = (SM + NFULL*NBD*NBD)/(NBD) - RELDEN*RELDEN
            END IF
            STD = SQRT(ABS(SM))
            THRESM = INT(9*RELDEN+0.5*STD*((STD/(RELDEN + 0.01))**1.5)
     *           + 2*RELDEN*RELDEN/(STD+0.01) +1)
            THRESM = MIN(THRESM,NBD)
            IF (THRESM.GE.NBD) THEN
               THRESM = N
            END IF
            NBD = NFULL
            GOTO 30
         ENDIF
         ELENME = ELEN (ME)
         ELEN (ME) = - (NEL + 1)
         NVPIV = NV (ME)
         NEL = NEL + NVPIV
         DENXT(ME) = 0
         NV (ME) = -NVPIV
         DEGME = 0
         IF (ELENME .EQ. 0) THEN
            PME1 = PE (ME)
            PME2 = PME1 - 1
            DO 60 P = PME1, PME1 + LEN (ME) - 1
               I = IW (P)
               NVI = NV (I)
               IF (NVI .GT. 0) THEN
                  DEGME = DEGME + NVI
                  NV (I) = -NVI
                  PME2 = PME2 + 1
                  IW (PME2) = I
                  IF (DEGREE(I).LE.N) THEN
                     ILAST = LAST (I)
                     INEXT = DENXT (I)
                     IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                     IF (ILAST .NE. 0) THEN
                        DENXT (ILAST) = INEXT
                     ELSE
                        HEAD (DEGREE (I)) = INEXT
                     ENDIF
                  ELSE
                     DENXT(ME) = DENXT(ME) + NVI
                  ENDIF
               ENDIF
 60         CONTINUE
            NEWMEM = 0
         ELSE
            P  = PE (ME)
            PME1 = PFREE
            SLENME = LEN (ME) - ELENME
            DO 120 KNT1 = 1, ELENME
               E = IW (P)
               P = P + 1
               PJ = PE (E)
               LN = LEN (E)
               DO 110 KNT2 = 1, LN
                  I = IW (PJ)
                  PJ = PJ + 1
                  NVI = NV (I)
                  IF (NVI .GT. 0) THEN
                     IF (PFREE .GT. IWLEN) THEN
                        PE (ME) = P
                        LEN (ME) = LEN (ME) - KNT1
                        IF (LEN (ME) .EQ. 0) PE (ME) = 0
                        PE (E) = PJ
                        LEN (E) = LN - KNT2
                        IF (LEN (E) .EQ. 0) PE (E) = 0
                        NCMPA = NCMPA + 1
                        DO 70 J = 1, N
                           PN = PE (J)
                           IF (PN .GT. 0) THEN
                              PE (J) = IW (PN)
                              IW (PN) = -J
                           ENDIF
 70                     CONTINUE
                        PDST = 1
                        PSRC = 1
                        PEND = PME1 - 1
                        DO 91 IDUMMY = 1, IWLEN
                           IF (PSRC .GT. PEND) GO TO 95
                           J = -IW (PSRC)
                           PSRC = PSRC + 1
                           IF (J .GT. 0) THEN
                              IW (PDST) = PE (J)
                              PE (J) = PDST
                              PDST = PDST + 1
                              LENJ = LEN (J)
                              DO 90 KNT3 = 0, LENJ - 2
                                 IW (PDST + KNT3) = IW (PSRC + KNT3)
 90                           CONTINUE
                              PDST = PDST + LENJ - 1
                              PSRC = PSRC + LENJ - 1
                           ENDIF
 91                     END DO
 95                     P1 = PDST
                        DO 100 PSRC = PME1, PFREE - 1
                           IW (PDST) = IW (PSRC)
                           PDST = PDST + 1
 100                    CONTINUE
                        PME1 = P1
                        PFREE = PDST
                        PJ = PE (E)
                        P = PE (ME)
                     ENDIF
                     DEGME = DEGME + NVI
                     NV (I) = -NVI
                     IW (PFREE) = I
                     PFREE = PFREE + 1
                     IF (DEGREE(I).LE.N) THEN
                        ILAST = LAST (I)
                        INEXT = DENXT (I)
                        IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                        IF (ILAST .NE. 0) THEN
                           DENXT (ILAST) = INEXT
                        ELSE
                           HEAD (DEGREE (I)) = INEXT
                        ENDIF
                     ELSE
                        DENXT(ME) = DENXT(ME) + NVI
                     ENDIF
                  ENDIF
 110           CONTINUE
                  PE (E) = -ME
                  W (E) = 0
 120        CONTINUE
            KNT1 = ELENME + 1
            E = ME
            PJ = P
            LN = SLENME
            DO 126 KNT2 = 1, LN
               I = IW (PJ)
               PJ = PJ + 1
               NVI = NV (I)
               IF (NVI .GT. 0) THEN
                  IF (PFREE .GT. IWLEN) THEN
                     PE (ME) = P
                     LEN (ME) = LEN (ME) - KNT1
                     IF (LEN (ME) .EQ. 0) PE (ME) = 0
                     PE (E) = PJ
                     LEN (E) = LN - KNT2
                     IF (LEN (E) .EQ. 0) PE (E) = 0
                     NCMPA = NCMPA + 1
                     DO 121 J = 1, N
                        PN = PE (J)
                        IF (PN .GT. 0) THEN
                           PE (J) = IW (PN)
                           IW (PN) = -J
                        ENDIF
 121                 CONTINUE
                     PDST = 1
                     PSRC = 1
                     PEND = PME1 - 1
                     DO 123 IDUMMY = 1,IWLEN
                        IF (PSRC .GT. PEND) GO TO 124
                        J = -IW (PSRC)
                        PSRC = PSRC + 1
                        IF (J .GT. 0) THEN
                           IW (PDST) = PE (J)
                           PE (J) = PDST
                           PDST = PDST + 1
                           LENJ = LEN (J)
                           DO 122 KNT3 = 0, LENJ - 2
                              IW (PDST + KNT3) = IW (PSRC + KNT3)
 122                       CONTINUE
                           PDST = PDST + LENJ - 1
                           PSRC = PSRC + LENJ - 1
                        ENDIF
 123                 END DO
 124                 P1 = PDST
                     DO 125 PSRC = PME1, PFREE - 1
                        IW (PDST) = IW (PSRC)
                        PDST = PDST + 1
 125                 CONTINUE
                     PME1 = P1
                     PFREE = PDST
                     PJ = PE (E)
                     P = PE (ME)
                  END IF
                  DEGME = DEGME + NVI
                  NV (I) = -NVI
                  IW (PFREE) = I
                  PFREE = PFREE + 1
                  IF (DEGREE(I).LE.N) THEN
                     ILAST = LAST (I)
                     INEXT = DENXT (I)
                     IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                     IF (ILAST .NE. 0) THEN
                        DENXT (ILAST) = INEXT
                     ELSE
                        HEAD (DEGREE (I)) = INEXT
                     ENDIF
                  ELSE
                     DENXT(ME) = DENXT(ME) + NVI
                  ENDIF
               ENDIF
 126           CONTINUE
            PME2 = PFREE - 1
            NEWMEM = PFREE - PME1
            MEM = MEM + NEWMEM
            MAXMEM = MAX (MAXMEM, MEM)
         ENDIF
         DEGREE (ME) = DEGME
         PE (ME) = PME1
         LEN (ME) = PME2 - PME1 + 1
         IF (WFLG .GT. IOVFLO-N) THEN
            DO 130 X = 1, N
               IF (W (X) .NE. 0) W (X) = 1
 130        CONTINUE
            WFLG = 2
         ENDIF
         IF (NBD.GT.0) THEN
            DO 150 PME = PME1, PME2
               I = IW (PME)
               IF (DEGREE(I).GT.N) GOTO 150
               ELN = ELEN (I)
               IF (ELN .GT. 0) THEN
                  NVI = -NV (I)
                  WNVI = WFLG - NVI
                  DO 140 P = PE (I), PE (I) + ELN - 1
                     E = IW (P)
                     WE = W (E)
                     IF (WE .GE. WFLG) THEN
                        WE = WE - NVI
                     ELSE IF (WE .NE. 0) THEN
                        WE = DEGREE (E) + WNVI - DENXT(E)
                     ENDIF
                     W (E) = WE
 140              CONTINUE
               ENDIF
 150        CONTINUE
         ELSE
            DO 152 PME = PME1, PME2
               I = IW (PME)
               ELN = ELEN (I)
               IF (ELN .GT. 0) THEN
                  NVI = -NV (I)
                  WNVI = WFLG - NVI
                  DO 151 P = PE (I), PE (I) + ELN - 1
                     E = IW (P)
                     WE = W (E)
                     IF (WE .GE. WFLG) THEN
                        WE = WE - NVI
                     ELSE IF (WE .NE. 0) THEN
                        WE = DEGREE (E) + WNVI
                     ENDIF
                     W (E) = WE
 151              CONTINUE
               ENDIF
 152        CONTINUE
         END IF
         IF (NBD.GT.0) THEN
            DO 180 PME = PME1, PME2
               I = IW (PME)
               IF (DEGREE(I).GT.N) GOTO 180
               P1 = PE (I)
               P2 = P1 + ELEN (I) - 1
               PN = P1
               HASH = 0
               DEG = 0
               DO 160 P = P1, P2
                  E = IW (P)
                  DEXT = W (E) - WFLG
                  IF (DEXT .GT. 0) THEN
                     DEG = DEG + DEXT
                     IW (PN) = E
                     PN = PN + 1
                     HASH = HASH+E
                  ELSE IF ((DEXT .EQ. 0) .AND.
     &                    (DENXT(ME).EQ.NBD)) THEN
                     PE (E) = -ME
                     W (E)  = 0
                  ELSE IF (DEXT.EQ.0) THEN
                     IW(PN) = E
                     PN     = PN+1
                     HASH = HASH + E
                  ENDIF
 160           CONTINUE
               ELEN (I) = PN - P1 + 1
               P3 = PN
               DO 170 P = P2 + 1, P1 + LEN (I) - 1
                  J = IW (P)
                  NVJ = NV (J)
                  IF (NVJ .GT. 0) THEN
                     IF (DEGREE(J).LE.N) DEG=DEG+NVJ
                     IW (PN) = J
                     PN = PN + 1
                     HASH = HASH + J
                  ENDIF
 170           CONTINUE
               IF ((DEG .EQ. 0).AND.(DENXT(ME).EQ.NBD)) THEN
                  PE (I) = -ME
                  NVI = -NV (I)
                  DEGME = DEGME - NVI
                  NVPIV = NVPIV + NVI
                  NEL = NEL + NVI
                  NV (I) = 0
                  ELEN (I) = 0
               ELSE
                  DEGREE(I) = MIN (DEG+NBD-DENXT(ME), DEGREE(I))
                  IW (PN) = IW (P3)
                  IW (P3) = IW (P1)
                  IW (P1) = ME
                  LEN (I) = PN - P1 + 1
                  HASH = ABS(MOD (HASH, HMOD)) + 1
                  J = HEAD (HASH)
                  IF (J .LE. 0) THEN
                     DENXT (I) = -J
                     HEAD (HASH) = -I
                  ELSE
                     DENXT (I) = LAST (J)
                     LAST (J) = I
                  ENDIF
                  LAST (I) = HASH
               ENDIF
 180        CONTINUE
         ELSE
            DO 183 PME = PME1, PME2
               I = IW (PME)
               P1 = PE (I)
               P2 = P1 + ELEN (I) - 1
               PN = P1
               HASH = 0
               DEG = 0
               DO 181 P = P1, P2
                  E = IW (P)
                  DEXT = W (E) - WFLG
                  IF (DEXT .GT. 0) THEN
                     DEG = DEG + DEXT
                     IW (PN) = E
                     PN = PN + 1
                     HASH = HASH + E
                  ELSE IF (DEXT .EQ. 0) THEN
                     PE (E) = -ME
                     W (E)  = 0
                  ENDIF
 181           CONTINUE
               ELEN (I) = PN - P1 + 1
               P3 = PN
               DO 182 P = P2 + 1, P1 + LEN (I) - 1
                  J = IW (P)
                  NVJ = NV (J)
                  IF (NVJ .GT. 0) THEN
                     DEG=DEG+NVJ
                     IW (PN) = J
                     PN = PN + 1
                     HASH = HASH + J
                  ENDIF
 182           CONTINUE
               IF (DEG .EQ. 0) THEN
                  PE (I) = -ME
                  NVI = -NV (I)
                  DEGME = DEGME - NVI
                  NVPIV = NVPIV + NVI
                  NEL = NEL + NVI
                  NV (I) = 0
                  ELEN (I) = 0
               ELSE
                  DEGREE(I) = MIN (DEG,  DEGREE(I))
                  IW (PN) = IW (P3)
                  IW (P3) = IW (P1)
                  IW (P1) = ME
                  LEN (I) = PN - P1 + 1
                  HASH = ABS(MOD (HASH, HMOD)) + 1
                  J = HEAD (HASH)
                  IF (J .LE. 0) THEN
                     DENXT (I) = -J
                     HEAD (HASH) = -I
                  ELSE
                     DENXT (I) = LAST (J)
                     LAST (J) = I
                  ENDIF
                  LAST (I) = HASH
               ENDIF
 183        CONTINUE
         END IF
         DEGREE (ME) = DEGME
         DMAX = MAX (DMAX, DEGME)
         WFLG = WFLG + DMAX
         IF (WFLG .GE. IOVFLO - N) THEN
            DO 190 X = 1, N
               IF (W (X) .NE. 0) W (X) = 1
 190        CONTINUE
            WFLG = 2
         ENDIF
         DO 250 PME = PME1, PME2
            I = IW (PME)
            IF ( (NV(I).GE.0) .OR. (DEGREE(I).GT.N) ) GO TO 250
            HASH = LAST (I)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
               I = -J
               HEAD (HASH) = 0
            ELSE
               I = LAST (J)
               LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
            DO 247 JDUMMY = 1,N
               IF (DENXT (I) .EQ. 0) GO TO 250
               LN = LEN (I)
               ELN = ELEN (I)
               DO 210 P = PE (I) + 1, PE (I) + LN - 1
                  W (IW (P)) = WFLG
 210           CONTINUE
               JLAST = I
               J = DENXT (I)
               DO 245 IDUMMY=1,N
                  IF (J .EQ. 0) GO TO 246
                  IF (LEN (J) .NE. LN) GO TO 240
                  IF (ELEN (J) .NE. ELN) GO TO 240
                  DO 230 P = PE (J) + 1, PE (J) + LN - 1
                     IF (W (IW (P)) .NE. WFLG) GO TO 240
 230              CONTINUE
                  PE (J) = -I
                  NV (I) = NV (I) + NV (J)
                  NV (J) = 0
                  ELEN (J) = 0
                  J = DENXT (J)
                  DENXT (JLAST) = J
                  GO TO 245
 240              CONTINUE
                  JLAST = J
                  J = DENXT (J)
 245           CONTINUE
 246           WFLG = WFLG + 1
               I = DENXT (I)
               IF (I .EQ. 0) GO TO 250
 247         CONTINUE
 250      CONTINUE
          P = PME1
          NLEFT = N - NEL
          DO 260 PME = PME1, PME2
             I = IW (PME)
             NVI = -NV (I)
             IF (NVI .LE. 0) GO TO 260
             NV (I) = NVI
             IF (DEGREE(I).GT.N) GO TO 258
             DEG = MIN (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
             DEGREE (I) = DEG
             IDENSE = .FALSE.
             IF (THRESM.GE.0) THEN
                IF ((DEG+NVI .GE. THRESM).OR.
     &               (DEG+NVI .GE. NLEFT)) THEN
                   IF (THRESM.EQ.N) THEN
                      IF ((ELEN(I).LE.2) .AND.((DEG+NVI).EQ.NLEFT)
     &                     .AND. NBD.EQ.NFULL ) THEN
                         DEGREE(I) = N+1
                         IDENSE = .TRUE.
                      ENDIF
                   ELSE
                      IDENSE = .TRUE.
                      IF ((ELEN(I).LE.2).AND. ((DEG+NVI).EQ.NLEFT)
     &                     .AND. NBD.EQ.NFULL ) THEN
                         DEGREE(I) = N+1
                      ELSE
                         DEGREE(I) = N+1+DEGREE(I)
                      ENDIF
                   ENDIF
                ENDIF
                IF (IDENSE) THEN
                   P1 = PE(I)
                   P2 = P1 + ELEN(I) - 1
                   DO 255 PJ=P1,P2
                      E= IW(PJ)
                      DENXT (E) = DENXT(E) + NVI
 255               CONTINUE
                   NBD = NBD+NVI
                   DEG = N
                   IF (DEGREE(I).EQ.N+1) THEN
                      NFULL = NFULL +NVI
                      IF (LASTD.EQ.0) THEN
                         LASTD     = I
                         HEAD(DEG) = I
                         DENXT(I)   = 0
                         LAST(I)   = 0
                      ELSE
                            DENXT(LASTD) = I
                            LAST(I)     = LASTD
                            LASTD       = I
                            DENXT(I)     = 0
                      ENDIF
                   ELSE
                      INEXT = HEAD(DEG)
                      IF (INEXT .NE. 0) LAST (INEXT) = I
                      DENXT (I) = INEXT
                      HEAD (DEG) = I
                      LAST(I)    = 0
                      IF (LASTD.EQ.0) LASTD=I
                   ENDIF
                ENDIF
             ENDIF
             IF (.NOT.IDENSE) THEN
                INEXT = HEAD (DEG)
                IF (INEXT .NE. 0) LAST (INEXT) = I
                DENXT (I) = INEXT
                LAST (I) = 0
                HEAD (DEG) = I
             ENDIF
             MINDEG = MIN (MINDEG, DEG)
 258         CONTINUE
             IW (P) = I
             P = P + 1
 260      CONTINUE
          OPS = OPS + DEGME*NVPIV + DEGME * NVPIV*NVPIV +
     *         DEGME*DEGME*NVPIV + NVPIV*NVPIV*NVPIV/3 +
     *         NVPIV*NVPIV/2 + NVPIV/6 + NVPIV
          NRLADU = NRLADU + (NVPIV*(NVPIV+1))/2 + (DEGME*NVPIV)
          NV (ME) = NVPIV + DEGME
          LEN (ME) = P - PME1
          IF (LEN (ME) .EQ. 0) THEN
             PE (ME) = 0
             W (ME) = 0
          ENDIF
          IF (NEWMEM .NE. 0) THEN
             PFREE = P
             MEM = MEM - NEWMEM + LEN (ME)
          ENDIF
          GO TO 30
       ENDIF
       GO TO 265
 263   NELME    = -(NEL+1)
       DO 264 X=1,N
          IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = -ME
          ELSEIF (DEGREE(X).EQ.N+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = -ME
             ELEN(X) = 0
             NV(X) = 0
          ENDIF
 264   CONTINUE
       ELEN(ME) = NELME
       NV(ME)   = NBD
       NRLADU = NRLADU + (NBD*(NBD+1))/2
       OPS = OPS + NBD*NBD*NBD/3 + NBD*NBD/2 + NBD/6 + NBD
       PE(ME)   = 0
 265   CONTINUE
       DO 290 I = 1, N
          IF (ELEN (I) .EQ. 0) THEN
             J = -PE (I)
             DO 270 JDUMMY = 1,N
                IF (ELEN (J) .LT. 0) GO TO 275
                J = -PE (J)
 270         CONTINUE
 275         E = J
             K = -ELEN (E)
             J = I
             DO 280 IDUMMY = 1,N
                IF (ELEN (J) .LT. 0) GO TO 285
                JNEXT = -PE (J)
                PE (J) = -E
                IF (ELEN (J) .EQ. 0) THEN
                   ELEN (J) = K
                   K = K + 1
                ENDIF
                J = JNEXT
 280         CONTINUE
 285         ELEN (E) = -K
          ENDIF
 290   CONTINUE
       DO 300 I = 1, N
          K = ABS (ELEN (I))
          LAST (K) = I
          ELEN (I) = K
 300   CONTINUE
       RJNFO(1) = OPS
       RJNFO(2) = NRLADU
       JNFO(1) = NCMPA
       JNFO(2) = RSTRT
       PFREE = MAXMEM
       RETURN
       END
* *******************************************************************
* COPYRIGHT (c) 1988 Hyprotech UK and
* Council for the Central Laboratory of the Research Councils
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
C Original date 14 June 2001
C  June 2001: threadsafe version of MC41
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC71AD(N,KASE,X,EST,W,IW,KEEP)
      INTEGER ITMAX
      PARAMETER (ITMAX=5)
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION EST
      INTEGER KASE,N
      DOUBLE PRECISION W(*),X(*)
      INTEGER IW(*),KEEP(5)
      DOUBLE PRECISION ALTSGN,TEMP
      INTEGER I,ITER,J,JLAST,JUMP
      INTEGER IDAMAX
      EXTERNAL IDAMAX
      INTRINSIC ABS,SIGN,NINT,DBLE
      IF (N.LE.0) THEN
        KASE = -1
        RETURN
      END IF
      IF (KASE.EQ.0) THEN
        DO 10 I = 1,N
          X(I) = ONE/DBLE(N)
   10   CONTINUE
        KASE = 1
        JUMP = 1
        KEEP(1) = JUMP
        KEEP(2) = 0
        KEEP(3) = 0
        KEEP(4) = 0
        RETURN
      END IF
      JUMP  = KEEP(1)
      ITER  = KEEP(2)
      J     = KEEP(3)
      JLAST = KEEP(4)
      GO TO (100,200,300,400,500) JUMP
  100 CONTINUE
      IF (N.EQ.1) THEN
        W(1) = X(1)
        EST = ABS(W(1))
        GO TO 510
      END IF
      DO 110 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  110 CONTINUE
      KASE = 2
      JUMP = 2
      GO TO 1010
  200 CONTINUE
      J = IDAMAX(N,X,1)
      ITER = 2
  220 CONTINUE
      DO 230 I = 1,N
        X(I) = ZERO
  230 CONTINUE
      X(J) = ONE
      KASE = 1
      JUMP = 3
      GO TO 1010
  300 CONTINUE
      DO 310 I = 1,N
        W(I) = X(I)
  310 CONTINUE
      DO 320 I = 1,N
        IF (NINT(SIGN(ONE,X(I))).NE.IW(I)) GO TO 330
  320 CONTINUE
      GO TO 410
  330 CONTINUE
      DO 340 I = 1,N
        X(I) = SIGN(ONE,X(I))
        IW(I) = NINT(X(I))
  340 CONTINUE
      KASE = 2
      JUMP = 4
      GO TO 1010
  400 CONTINUE
      JLAST = J
      J = IDAMAX(N,X,1)
      IF ((ABS(X(JLAST)).NE.ABS(X(J))) .AND. (ITER.LT.ITMAX)) THEN
        ITER = ITER + 1
        GO TO 220
      END IF
  410 CONTINUE
      EST = ZERO
      DO 420 I = 1,N
        EST = EST + ABS(W(I))
  420 CONTINUE
      ALTSGN = ONE
      DO 430 I = 1,N
        X(I) = ALTSGN* (ONE+DBLE(I-1)/DBLE(N-1))
        ALTSGN = -ALTSGN
  430 CONTINUE
      KASE = 1
      JUMP = 5
      GO TO 1010
  500 CONTINUE
      TEMP = ZERO
      DO 520 I = 1,N
        TEMP = TEMP + ABS(X(I))
  520 CONTINUE
      TEMP = 2.0*TEMP/DBLE(3*N)
      IF (TEMP.GT.EST) THEN
        DO 530 I = 1,N
          W(I) = X(I)
  530   CONTINUE
        EST = TEMP
      END IF
  510 KASE = 0
 1010 CONTINUE
      KEEP(1) = JUMP
      KEEP(2) = ITER
      KEEP(3) = J
      KEEP(4) = JLAST
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1987 Hyprotech UK
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
* Original date 10 Feb 1993
C       Toolpack tool decs employed.
C 20/2/02 Cosmetic changes applied to reduce single/double differences

C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC34AD(N,IRN,JCOLST,YESA,A,IW)
      INTEGER N
      LOGICAL YESA
      DOUBLE PRECISION A(*)
      INTEGER IRN(*),IW(*),JCOLST(*)
      INTEGER CKP1,I,I1,I2,II,IPKP1,IPOS,J,JSTART,LENK,NDIAG,NEWTAU,
     +        OLDTAU
      OLDTAU = JCOLST(N+1) - 1
      DO 5 I = 1,N
        IW(I) = 0
    5 CONTINUE
      NDIAG = 0
      DO 20 J = 1,N
        I1 = JCOLST(J)
        I2 = JCOLST(J+1) - 1
        IW(J) = IW(J) + I2 - I1 + 1
        DO 10 II = I1,I2
          I = IRN(II)
          IF (I.NE.J) THEN
            IW(I) = IW(I) + 1
          ELSE
            NDIAG = NDIAG + 1
          END IF
   10   CONTINUE
   20 CONTINUE
      NEWTAU = 2*OLDTAU - NDIAG
      IPKP1 = OLDTAU + 1
      CKP1 = NEWTAU + 1
      DO 40 J = N,1,-1
        I1 = JCOLST(J)
        I2 = IPKP1
        LENK = I2 - I1
        JSTART = CKP1
        IPKP1 = I1
        I2 = I2 - 1
        DO 30 II = I2,I1,-1
          JSTART = JSTART - 1
          IF (YESA) A(JSTART) = A(II)
          IRN(JSTART) = IRN(II)
   30   CONTINUE
        JCOLST(J) = JSTART
        CKP1 = CKP1 - IW(J)
        IW(J) = LENK
   40 CONTINUE
      DO 80 J = N,1,-1
        I1 = JCOLST(J)
        I2 = JCOLST(J) + IW(J) - 1
        DO 60 II = I1,I2
          I = IRN(II)
          IF (I.EQ.J) GO TO 60
          JCOLST(I) = JCOLST(I) - 1
          IPOS = JCOLST(I)
          IF (YESA) A(IPOS) = A(II)
          IRN(IPOS) = J
   60   CONTINUE
   80 CONTINUE
      JCOLST(N+1) = NEWTAU + 1
      RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1993 Council for the Central Laboratory
*                    of the Research Councils
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

C Original date 29 Jan 2001
C 29 January 2001. Modified from MC49 to be threadsafe.

C 12th July 2004 Version 1.0.0. Version numbering added.
C 28 February 2008. Version 1.0.1. Comments flowed to column 72.
C 21 September 2009. Version 1.0.2. Minor change to documentation.

      SUBROUTINE MC59AD(ICNTL,NC,NR,NE,IRN,LJCN,JCN,LA,A,LIP,IP,
     &                  LIW,IW,INFO)
      INTEGER LA,LIP,LIW,LJCN,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER ICNTL(10),IP(LIP),INFO(10),IRN(NE),IW(LIW),JCN(LJCN)
      INTEGER I,ICNTL1,ICNTL2,ICNTL3,ICNTL6,LAA
      INTEGER IDUP,IOUT,IUP,JOUT,LP,MP,KNE,PART
      LOGICAL LCHECK
      EXTERNAL MC59BD,MC59CD,MC59DD,MC59ED,MC59FD
      INTRINSIC MAX
      DO 10 I = 1,10
         INFO(I) = 0
   10 CONTINUE
      ICNTL1 = ICNTL(1)
      ICNTL2 = ICNTL(2)
      ICNTL3 = ICNTL(3)
      ICNTL6 = ICNTL(6)
      LCHECK = (ICNTL1.EQ.0)
      LP = ICNTL(4)
      MP = ICNTL(5)
      IF (ICNTL2.GT.2 .OR. ICNTL2.LT.0) THEN
         INFO(1) = -1
         INFO(2) = ICNTL2
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9010) ICNTL2
         END IF
         GO TO 70
      END IF
      IF (ICNTL6.GT.2 .OR. ICNTL6.LT.-2) THEN
         INFO(1) = -11
         INFO(2) = ICNTL6
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9150) ICNTL6
         END IF
         GO TO 70
      END IF
      IF (NC.LT.1) THEN
        INFO(1) = -2
        INFO(2) = NC
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9020) NC
        END IF
        GO TO 70
      END IF
      IF (NR.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9030) NR
        END IF
        GO TO 70
      END IF
      IF (ICNTL6.NE.0 .AND. NR.NE.NC) THEN
        INFO(1) = -3
        INFO(2) = NR
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9035) NC,NR
        END IF
        GO TO 70
      END IF
      IF (NE.LT.1) THEN
        INFO(1) = -4
        INFO(2) = NE
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9040) NE
        END IF
        GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.1) THEN
        IF (LJCN.LT.NE) THEN
          INFO(1) = -5
          INFO(2) = NE
        END IF
      ELSE
        IF (LJCN.LT.1) THEN
          INFO(1) = -5
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-5) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9050) LJCN,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL3.EQ.0) THEN
        IF (LA.LT.NE) THEN
          INFO(1) = -6
          INFO(2) = NE
        END IF
      ELSE
        IF (LA.LT.1) THEN
          INFO(1) = -6
          INFO(2) = 1
        END IF
      END IF
      IF (INFO(1).EQ.-6) THEN
         IF (LP.GT.0) THEN
            WRITE (LP,FMT=9000) INFO(1)
            WRITE (LP,FMT=9060) LA,INFO(2)
         END IF
         GO TO 70
      END IF
      IF (ICNTL2.EQ.0 .OR. ICNTL2.EQ.2) THEN
        IF (LIP.LT.NC+1) THEN
          INFO(1) = -7
          INFO(2) = NC+1
        END IF
      ELSE IF (LIP.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -7
        INFO(2) = MAX(NR,NC)+1
      END IF
      IF (INFO(1).EQ.-7) THEN
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9065) LIP,INFO(2)
        END IF
        GO TO 70
      END IF
      IF (LIW.LT.MAX(NR,NC)+1) THEN
        INFO(1) = -8
        INFO(2) = MAX(NR,NC)+1
        IF (LP.GT.0) THEN
          WRITE (LP,FMT=9000) INFO(1)
          WRITE (LP,FMT=9070) LIW,INFO(2)
        END IF
        GO TO 70
      END IF
      LAA = NE
      IF (ICNTL3.NE.0) LAA = 1
      IOUT = 0
      JOUT = 0
      IDUP = 0
      IUP = 0
      PART = 0
      IF (ICNTL6.NE.0) PART = 1
      IF (ICNTL2.EQ.0) THEN
        CALL MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LAA,A,IP,IW,
     +              IOUT,JOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59ED(NC,NR,NE,IRN,LIP,IP,LAA,A,IW,IDUP,
     &                          KNE,ICNTL6)
      ELSE IF (ICNTL2.EQ.1) THEN
        IF (ICNTL6.NE.0) PART = -1
        CALL MC59BD(LCHECK,PART,NR,NC,NE,JCN,IRN,LAA,A,IW,IP,
     +              JOUT,IOUT,KNE)
        IF (KNE.EQ.0) GO TO 50
        IF (LCHECK) CALL MC59ED(NR,NC,NE,JCN,NR+1,IW,LAA,A,IP,
     &                          IDUP,KNE,ICNTL6)
        CALL MC59CD(NC,NR,KNE,IRN,JCN,LAA,A,IP,IW)
      ELSE IF (ICNTL2.EQ.2) THEN
        IF (LCHECK) THEN
          CALL MC59FD(NC,NR,NE,IRN,NC+1,IP,LAA,A,LIW,IW,IDUP,
     +                IOUT,IUP,KNE,ICNTL6,INFO)
          IF (INFO(1).EQ.-9) GO TO 40
          IF (KNE.EQ.0) GO TO 50
        ELSE
           KNE = NE
        END IF
        CALL MC59DD(NC,KNE,IRN,IP,LAA,A)
      END IF
      INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(6) = KNE
      INFO(7) = IUP
      IF (IDUP.GT.0) INFO(1) = INFO(1) + 1
      IF (IOUT.GT.0) INFO(1) = INFO(1) + 2
      IF (JOUT.GT.0) INFO(1) = INFO(1) + 4
      IF (INFO(1).GT.0 .AND. MP.GT.0) THEN
        WRITE (MP,FMT=9080) INFO(1)
        IF (IOUT.GT.0) WRITE (MP,FMT=9090) IOUT
        IF (JOUT.GT.0) WRITE (MP,FMT=9110) JOUT
        IF (IDUP.GT.0) WRITE (MP,FMT=9100) IDUP
        IF (IUP.GT.0)  WRITE (MP,FMT=9130) IUP
      END IF
      GO TO 70
   40 INFO(3) = IDUP
      INFO(4) = IOUT
      INFO(7) = IUP
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9140)
      END IF
      GO TO 70
   50 INFO(1) = -10
      INFO(4) = IOUT
      INFO(5) = JOUT
      INFO(2) = IOUT + JOUT
      IF (LP.GT.0) THEN
        WRITE (LP,FMT=9000) INFO(1)
        WRITE (LP,FMT=9120)
      END IF
   70 RETURN
 9000 FORMAT (/,' *** Error return from MC59AD *** INFO(1) = ',I3)
 9010 FORMAT (1X,'ICNTL(2) = ',I2,' is out of range')
 9020 FORMAT (1X,'NC = ',I6,' is out of range')
 9030 FORMAT (1X,'NR = ',I6,' is out of range')
 9035 FORMAT (1X,'Symmetric case. NC = ',I6,' but NR = ',I6)
 9040 FORMAT (1X,'NE = ',I10,' is out of range')
 9050 FORMAT (1X,'Increase LJCN from ',I10,' to at least ',I10)
 9060 FORMAT (1X,'Increase LA from ',I10,' to at least ',I10)
 9065 FORMAT (1X,'Increase LIP from ',I8,' to at least ',I10)
 9070 FORMAT (1X,'Increase LIW from ',I8,' to at least ',I10)
 9080 FORMAT (/,' *** Warning message from MC59AD *** INFO(1) = ',I3)
 9090 FORMAT (1X,I8,' entries in IRN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9100 FORMAT (1X,I8,' duplicate entries were supplied by the user')
 9110 FORMAT (1X,I8,' entries in JCN supplied by the user were ',
     +       /,'       out of range and were ignored by the routine')
 9120 FORMAT (1X,'All entries out of range')
 9130 FORMAT (1X,I8,' of these entries were in the upper triangular ',
     +       /,'       part of matrix')
 9140 FORMAT (1X,'Entries in IP are not monotonic increasing')
 9150 FORMAT (1X,'ICNTL(6) = ',I2,' is out of range')
      END
C***********************************************************************
      SUBROUTINE MC59BD(LCHECK,PART,NC,NR,NE,IRN,JCN,LA,A,IP,IW,IOUT,
     +                  JOUT,KNE)
      INTEGER LA,NC,NE,NR,IOUT,JOUT,KNE,PART
      LOGICAL LCHECK
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NC+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,JCE,JCEP,K,L,LOC
      DO 10 J = 1,NC + 1
        IW(J) = 0
   10 CONTINUE
      KNE = 0
      IOUT = 0
      JOUT = 0
      IF (LCHECK) THEN
        IF (LA.GT.1) THEN
          IF (PART.EQ.0) THEN
            DO 20 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                A(KNE) = A(K)
                IW(J) = IW(J) + 1
              END IF
   20       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 21 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   21       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 22 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
                A(KNE) = A(K)
              END IF
   22       CONTINUE
          END IF
        ELSE
          IF (PART.EQ.0) THEN
            DO 25 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IRN(KNE) = I
                JCN(KNE) = J
                IW(J) = IW(J) + 1
              END IF
   25       CONTINUE
          ELSE IF (PART.EQ.1) THEN
            DO 26 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.LT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   26       CONTINUE
          ELSE IF (PART.EQ.-1) THEN
            DO 27 K = 1,NE
              I = IRN(K)
              J = JCN(K)
              IF (I.GT.NR .OR. I.LT.1) THEN
                IOUT = IOUT + 1
                IF (J.GT.NC .OR. J.LT.1)  JOUT = JOUT + 1
              ELSE IF (J.GT.NC .OR. J.LT.1) THEN
                JOUT = JOUT + 1
              ELSE
                KNE = KNE + 1
                IF (I.GT.J) THEN
                  IRN(KNE) = J
                  JCN(KNE) = I
                  IW(I) = IW(I) + 1
                ELSE
                  IRN(KNE) = I
                  JCN(KNE) = J
                  IW(J) = IW(J) + 1
                END IF
              END IF
   27       CONTINUE
          END IF
        END IF
        IF (KNE.EQ.0) GO TO 130
      ELSE
        KNE = NE
        IF (PART.EQ.0) THEN
          DO 30 K = 1,NE
            J = JCN(K)
            IW(J) = IW(J) + 1
   30     CONTINUE
        ELSE IF (PART.EQ.1) THEN
          DO 35 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.LT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   35     CONTINUE
        ELSE IF (PART.EQ.-1) THEN
          DO 36 K = 1,NE
            I = IRN(K)
            J = JCN(K)
            IF (I.GT.J) THEN
               IRN(K) = J
               JCN(K) = I
               IW(I) = IW(I) + 1
            ELSE
              IW(J) = IW(J) + 1
            END IF
   36     CONTINUE
        END IF
      END IF
      IP(1) = 1
      DO 37 J = 2,NC + 1
        IP(J) = IW(J-1) + IP(J-1)
        IW(J-1) = IP(J-1)
   37 CONTINUE
      IF (LA.EQ.1) THEN
        DO 70 L = 1,NC
          DO 60 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            DO 40 J = 1,NE
              IF (JCE.EQ.L) GO TO 50
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
   40       CONTINUE
   50       JCN(K) = JCE
            IRN(K) = ICE
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 120 L = 1,NC
          DO 110 K = IW(L),IP(L+1) - 1
            ICE = IRN(K)
            JCE = JCN(K)
            ACE = A(K)
            DO 90 J = 1,NE
              IF (JCE.EQ.L) GO TO 100
              LOC = IW(JCE)
              JCEP = JCN(LOC)
              ICEP = IRN(LOC)
              IW(JCE) = LOC + 1
              JCN(LOC) = JCE
              IRN(LOC) = ICE
              JCE = JCEP
              ICE = ICEP
              ACEP = A(LOC)
              A(LOC) = ACE
              ACE = ACEP
   90       CONTINUE
  100       JCN(K) = JCE
            IRN(K) = ICE
            A(K) = ACE
  110     CONTINUE
  120   CONTINUE
      END IF
  130 CONTINUE
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59CD(NC,NR,NE,IRN,JCN,LA,A,IP,IW)
      INTEGER LA,NC,NE,NR
      DOUBLE PRECISION A(LA)
      INTEGER IP(NC+1),IRN(NE),IW(NR+1),JCN(NE)
      DOUBLE PRECISION ACE,ACEP
      INTEGER I,ICE,ICEP,J,J1,J2,K,L,LOC,LOCP
      DO 10 J = 1,NC
        IP(J) = 0
   10 CONTINUE
      IF (LA.GT.1) THEN
        DO 20 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
          IRN(K) = JCN(K)
   20   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 30 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
   30   CONTINUE
        DO 50 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 40 J = J1,J2
            K = IRN(J)
            L = IP(K) - 1
            JCN(J) = L
            IRN(J) = I
            IP(K) = L
   40     CONTINUE
   50   CONTINUE
        IP(NC+1) = NE + 1
        DO 70 J = 1,NE
          LOC = JCN(J)
          IF (LOC.EQ.0) GO TO 70
          ICE = IRN(J)
          ACE = A(J)
          JCN(J) = 0
          DO 60 K = 1,NE
            LOCP = JCN(LOC)
            ICEP = IRN(LOC)
            ACEP = A(LOC)
            JCN(LOC) = 0
            IRN(LOC) = ICE
            A(LOC) = ACE
            IF (LOCP.EQ.0) GO TO 70
            ICE = ICEP
            ACE = ACEP
            LOC = LOCP
   60     CONTINUE
   70   CONTINUE
      ELSE
        DO 90 K = 1,NE
          I = JCN(K)
          IP(I) = IP(I) + 1
   90   CONTINUE
        IP(NC+1) = NE + 1
        IP(1) = IP(1) + 1
        DO 100 J = 2,NC
          IP(J) = IP(J) + IP(J-1)
  100   CONTINUE
        DO 120 I = NR,1,-1
          J1 = IW(I)
          J2 = IW(I+1) - 1
          DO 110 J = J1,J2
            K = JCN(J)
            L = IP(K) - 1
            IRN(L) = I
            IP(K) = L
  110     CONTINUE
  120   CONTINUE
      END IF
      RETURN
      END
C**********************************************************
      SUBROUTINE MC59DD(NC,NE,IRN,IP,LA,A)
      INTEGER LA,NC,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(NC)
      DOUBLE PRECISION ACE
      INTEGER ICE,IK,J,JJ,K,KDUMMY,KLO,KMAX,KOR
      INTRINSIC ABS
      IF (LA.GT.1) THEN
        KMAX = NE
        DO 50 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 40
          KOR = KMAX
          DO 30 KDUMMY = KLO,KMAX
            ACE = A(KOR-1)
            ICE = IRN(KOR-1)
            DO 10 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 20
              IRN(K-1) = IK
              A(K-1) = A(K)
   10       CONTINUE
            K = KMAX + 1
   20       IRN(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
   30     CONTINUE
   40     KMAX = KLO - 2
   50   CONTINUE
      ELSE
        KMAX = NE
        DO 150 JJ = 1,NC
          J = NC + 1 - JJ
          KLO = IP(J) + 1
          IF (KLO.GT.KMAX) GO TO 140
          KOR = KMAX
          DO 130 KDUMMY = KLO,KMAX
            ICE = IRN(KOR-1)
            DO 110 K = KOR,KMAX
              IK = IRN(K)
              IF (ABS(ICE).LE.ABS(IK)) GO TO 120
              IRN(K-1) = IK
  110       CONTINUE
            K = KMAX + 1
  120       IRN(K-1) = ICE
            KOR = KOR - 1
  130     CONTINUE
  140     KMAX = KLO - 2
  150   CONTINUE
      END IF
      END
C***********************************************************************
      SUBROUTINE MC59ED(NC,NR,NE,IRN,LIP,IP,LA,A,IW,IDUP,KNE,ICNTL6)
      INTEGER ICNTL6,IDUP,KNE,LIP,LA,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(NR)
      INTEGER I,J,K,KSTART,KSTOP,NZJ
      IDUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          KSTOP = IP(J+1)
          IP(J+1) = IP(J)
          DO 40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE MC59FD(NC,NR,NE,IRN,LIP,IP,LA,A,LIW,IW,IDUP,IOUT,
     +                  IUP,KNE,ICNTL6,INFO)
      INTEGER ICNTL6,IDUP,IOUT,IUP,KNE,LA,LIP,LIW,NC,NR,NE
      DOUBLE PRECISION A(LA)
      INTEGER IRN(NE),IP(LIP),IW(LIW),INFO(2)
      INTEGER I,J,K,KSTART,KSTOP,NZJ,LOWER
      IDUP = 0
      IOUT = 0
      IUP = 0
      KNE = 0
      DO 10 I = 1,NR
        IW(I) = 0
   10 CONTINUE
      KSTART = IP(1)
      LOWER = 1
      IF (LA.GT.1) THEN
        NZJ = 0
        DO 30 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO 20 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.LT.J) IUP = IUP + 1
            ELSE IF (IW(I).LE.NZJ) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              A(KNE) = A(K)
              IP(J+1) = IP(J+1) + 1
              IW(I) = KNE
            ELSE
              IDUP = IDUP + 1
              IF (ICNTL6.GE.0) A(IW(I)) = A(IW(I)) + A(K)
            END IF
   20     CONTINUE
          KSTART = KSTOP
          NZJ = KNE
   30   CONTINUE
      ELSE
        DO 50 J = 1,NC
          IF (ICNTL6.NE.0) LOWER = J
          KSTOP = IP(J+1)
          IF (KSTART.GT.KSTOP) THEN
            INFO(1) = -9
            INFO(2) = J
            RETURN
          END IF
          IP(J+1) = IP(J)
          DO  40 K = KSTART,KSTOP - 1
            I = IRN(K)
            IF (I.GT.NR .OR. I.LT.LOWER) THEN
              IOUT = IOUT + 1
              IF (ICNTL6.NE.0 .AND. I.GT.1) IUP = IUP + 1
            ELSE IF (IW(I).LT.J) THEN
              KNE = KNE + 1
              IRN(KNE) = I
              IP(J+1) = IP(J+1) + 1
              IW(I) = J
            ELSE
              IDUP = IDUP + 1
            END IF
   40     CONTINUE
          KSTART = KSTOP
   50   CONTINUE
      END IF
      RETURN
      END
C *******************************************************************
C COPYRIGHT (c) 1999 Council for the Central Laboratory
*                    of the Research Councils
C All rights reserved.
C
C None of the comments in this Copyright notice between the lines
C of asterisks shall be removed or altered in any way.
C
C This Package is intended for compilation without modification,
C so most of the embedded comments have been removed.
C
C ALL USE IS SUBJECT TO LICENCE. For full details of the ACADEMIC
C SOFTWARE LICENCE, see http://hsl.rl.ac.uk/hsl2007/cou/academic.html
C
C Please note that for an ACADEMIC Licence:
C
C 1. The Packages may only be used for academic research or teaching
C    purposes by the Licensee, and must not be copied by the Licensee for
C    use by any other persons. Use of the Packages in any commercial
C    application shall be subject to prior written agreement between
C    Hyprotech UK Limited and the Licensee on suitable terms and
C    conditions, which will include financial conditions.
C 2. All information on the Package is provided to the Licensee on the
C    understanding that the details thereof are confidential.
C 3. All publications issued by the Licensee that include results obtained
C    with the help of one or more of the Packages shall acknowledge the
C    use of the Packages. The Licensee will notify the Numerical Analysis
C    Group at Rutherford Appleton Laboratory (STFC) of any such publication.
C 4. The Packages may be modified by or on behalf of the Licensee
C    for such use in research applications but at no time shall such
C    Packages or modifications thereof become the property of the
C    Licensee. The Licensee shall make available free of charge to the
C    copyright holder for any purpose all information relating to
C    any modification.
C 5. Neither STFC nor Hyprotech UK Limited shall be liable for any
C    direct or consequential loss or damage whatsoever arising out of
C    the use of Packages by the Licensee.
C *******************************************************************
C
C Original date July 1999
CCCCC PACKAGE MC64A/AD
CCCCC AUTHORS Iain Duff (i.duff@rl.ac.uk) and
CCCCC         Jacko Koster (jacko.koster@uninett.no)

C 12th July 2004 Version 1.0.0. Version numbering added.

C 30/07/04  Version 1.1.0. Permutation array flagged negative to
C           indicate dependent columns in singular case.  Calls to
C           MC64F/FD changed to avoid unsafe reference to array L.
C 21st February 2005 Version 1.2.0. FD05 dependence changed to FD15.
C  7th March 2005 Version 1.3.0. Scan of dense columns avoided in
C           the cheap assignment phase in MC64W/WD.
C 15th May 2007 Version 1.4.0. Minor change made in MC64W/WD to avoid
C           crash when the input array contains NaNs.
C 28th June 2007  Version 1.5.0. Permutation array flagged negative to
C           indicate dependent columns in singular case for all values
C           of the JOB parameter (previously was only for JOB 4 and 5).
C           Also some very minor changes concerning tests in DO loops
C           being moved from beginning to end of DO loop.


      SUBROUTINE MC64ID(ICNTL)
      IMPLICIT NONE
      INTEGER ICNTL(10)
      INTEGER I
      ICNTL(1) = 6
      ICNTL(2) = 6
      ICNTL(3) = -1
      DO 10 I = 4,10
        ICNTL(I) = 0
   10 CONTINUE
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64AD(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,
     &           ICNTL,INFO)
      IMPLICIT NONE
      INTEGER JOB,N,NE,NUM,LIW,LDW
      INTEGER IP(N+1),IRN(NE),CPERM(N),IW(LIW),ICNTL(10),INFO(10)
      DOUBLE PRECISION A(NE),DW(LDW)
      INTEGER I,J,K
      DOUBLE PRECISION FACT,ZERO,RINF
      PARAMETER (ZERO=0.0D+00)
      EXTERNAL FD15AD,MC21AD,MC64BD,MC64RD,MC64SD,MC64WD
      DOUBLE PRECISION FD15AD
      INTRINSIC ABS,LOG
      RINF = FD15AD('H')
      IF (JOB.LT.1 .OR. JOB.GT.5) THEN
        INFO(1) = -1
        INFO(2) = JOB
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'JOB',JOB
        GO TO 99
      ENDIF
      IF (N.LT.1) THEN
        INFO(1) = -2
        INFO(2) = N
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'N',N
        GO TO 99
      ENDIF
      IF (NE.LT.1) THEN
        INFO(1) = -3
        INFO(2) = NE
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9001) INFO(1),'NE',NE
        GO TO 99
      ENDIF
      IF (JOB.EQ.1) K = 5*N
      IF (JOB.EQ.2) K = 4*N
      IF (JOB.EQ.3) K = 10*N + NE
      IF (JOB.EQ.4) K = 5*N
      IF (JOB.EQ.5) K = 5*N
      IF (LIW.LT.K) THEN
        INFO(1) = -4
        INFO(2) = K
        IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9004) INFO(1),K
        GO TO 99
      ENDIF
      IF (JOB.GT.1) THEN
        IF (JOB.EQ.2) K = N
        IF (JOB.EQ.3) K = NE
        IF (JOB.EQ.4) K = 2*N + NE
        IF (JOB.EQ.5) K = 3*N + NE
        IF (LDW.LT.K) THEN
          INFO(1) = -5
          INFO(2) = K
          IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9005) INFO(1),K
          GO TO 99
        ENDIF
      ENDIF
      IF (ICNTL(4).EQ.0) THEN
        DO 3 I = 1,N
          IW(I) = 0
    3   CONTINUE
        DO 6 J = 1,N
          DO 4 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (I.LT.1 .OR. I.GT.N) THEN
              INFO(1) = -6
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9006) INFO(1),J,I
              GO TO 99
            ENDIF
            IF (IW(I).EQ.J) THEN
              INFO(1) = -7
              INFO(2) = J
              IF (ICNTL(1).GE.0) WRITE(ICNTL(1),9007) INFO(1),J,I
              GO TO 99
            ELSE
              IW(I) = J
            ENDIF
    4     CONTINUE
    6   CONTINUE
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9020) JOB,N,NE
        WRITE(ICNTL(3),9021) (IP(J),J=1,N+1)
        WRITE(ICNTL(3),9022) (IRN(J),J=1,NE)
        IF (JOB.GT.1) WRITE(ICNTL(3),9023) (A(J),J=1,NE)
      ENDIF
      DO 8 I=1,10
        INFO(I) = 0
    8 CONTINUE
      IF (JOB.EQ.1) THEN
        DO 10 J = 1,N
          IW(J) = IP(J+1) - IP(J)
   10   CONTINUE
        CALL MC21AD(N,IRN,NE,IP,IW(1),CPERM,NUM,IW(N+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.2) THEN
        CALL MC64BD(N,NE,IP,IRN,A,CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),DW)
        GO TO 90
      ENDIF
      IF (JOB.EQ.3) THEN
        DO 20 K = 1,NE
          IW(K) = IRN(K)
          DW(K) = ABS(A(K))
   20   CONTINUE
        CALL MC64RD(N,NE,IP,IW,DW)
        CALL MC64SD(N,NE,IP,IW(1),DW,CPERM,NUM,IW(NE+1),
     &     IW(NE+N+1),IW(NE+2*N+1),IW(NE+3*N+1),IW(NE+4*N+1),
     &     IW(NE+5*N+1),IW(NE+6*N+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.4) THEN
        DO 50 J = 1,N
          FACT = ZERO
          DO 30 K = IP(J),IP(J+1)-1
            IF (ABS(A(K)).GT.FACT) FACT = ABS(A(K))
   30     CONTINUE
          DO 40 K = IP(J),IP(J+1)-1
            DW(2*N+K) = FACT - ABS(A(K))
   40     CONTINUE
   50   CONTINUE
        CALL MC64WD(N,NE,IP,IRN,DW(2*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        GO TO 90
      ENDIF
      IF (JOB.EQ.5) THEN
        DO 75 J = 1,N
          FACT = ZERO
          DO 60 K = IP(J),IP(J+1)-1
            DW(3*N+K) = ABS(A(K))
            IF (DW(3*N+K).GT.FACT) FACT = DW(3*N+K)
   60     CONTINUE
          DW(2*N+J) = FACT
          IF (FACT.NE.ZERO) THEN
            FACT = LOG(FACT)
          ELSE
            FACT = RINF/N
          ENDIF
          DO 70 K = IP(J),IP(J+1)-1
            IF (DW(3*N+K).NE.ZERO) THEN
              DW(3*N+K) = FACT - LOG(DW(3*N+K))
            ELSE
              DW(3*N+K) = RINF/N
            ENDIF
   70     CONTINUE
   75   CONTINUE
        CALL MC64WD(N,NE,IP,IRN,DW(3*N+1),CPERM,NUM,
     &     IW(1),IW(N+1),IW(2*N+1),IW(3*N+1),IW(4*N+1),
     &     DW(1),DW(N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
            IF (DW(2*N+J).NE.ZERO) THEN
              DW(N+J) = DW(N+J) - LOG(DW(2*N+J))
            ELSE
              DW(N+J) = ZERO
            ENDIF
   80     CONTINUE
        ENDIF
        FACT = 0.5*LOG(RINF)
        DO 86 J = 1,N
          IF (DW(J).LT.FACT .AND. DW(N+J).LT.FACT) GO TO 86
          INFO(1) = 2
          GO TO 90
   86   CONTINUE
      ENDIF
   90 IF (INFO(1).EQ.0 .AND. NUM.LT.N) THEN
        INFO(1) = 1
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9011) INFO(1)
      ENDIF
      IF (INFO(1).EQ.2) THEN
        IF (ICNTL(2).GE.0) WRITE(ICNTL(2),9012) INFO(1)
      ENDIF
      IF (ICNTL(3).GE.0) THEN
        WRITE(ICNTL(3),9030) (INFO(J),J=1,2)
        WRITE(ICNTL(3),9031) NUM
        WRITE(ICNTL(3),9032) (CPERM(J),J=1,N)
        IF (JOB.EQ.5) THEN
          WRITE(ICNTL(3),9033) (DW(J),J=1,N)
          WRITE(ICNTL(3),9034) (DW(N+J),J=1,N)
        ENDIF
      ENDIF
   99 RETURN
 9001 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2,
     &        ' because ',(A),' = ',I10)
 9004 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LIW too small, must be at least ',I8)
 9005 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        LDW too small, must be at least ',I8)
 9006 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains an entry with invalid row index ',I8)
 9007 FORMAT (' ****** Error in MC64A/AD. INFO(1) = ',I2/
     &        '        Column ',I8,
     &        ' contains two or more entries with row index ',I8)
 9011 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        The matrix is structurally singular.')
 9012 FORMAT (' ****** Warning from MC64A/AD. INFO(1) = ',I2/
     &        '        Some scaling factors may be too large.')
 9020 FORMAT (' ****** Input parameters for MC64A/AD:'/
     &        ' JOB = ',I8/' N   = ',I8/' NE  = ',I8)
 9021 FORMAT (' IP(1:N+1)  = ',8I8/(14X,8I8))
 9022 FORMAT (' IRN(1:NE)  = ',8I8/(14X,8I8))
 9023 FORMAT (' A(1:NE)    = ',4(1PD14.4)/(14X,4(1PD14.4)))
 9030 FORMAT (' ****** Output parameters for MC64A/AD:'/
     &        ' INFO(1:2)  = ',2I8)
 9031 FORMAT (' NUM        = ',I8)
 9032 FORMAT (' CPERM(1:N) = ',8I8/(14X,8I8))
 9033 FORMAT (' DW(1:N)    = ',5(F11.3)/(14X,5(F11.3)))
 9034 FORMAT (' DW(N+1:2N) = ',5(F11.3)/(14X,5(F11.3)))
      END
C**********************************************************************
      SUBROUTINE MC64BD(N,NE,IP,IRN,A,IPERM,NUM,JPERM,PR,Q,L,D)
      IMPLICIT NONE
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),JPERM(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),D(N)
      INTEGER I,II,J,JJ,JORD,Q0,QLEN,IDUM,JDUM,ISP,JSP,
     &        K,KK,KK1,KK2,I0,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DNEW,DQ0,AI,A0,BV
      DOUBLE PRECISION RINF,ZERO,MINONE
      PARAMETER (ZERO=0.0D+0,MINONE=-1.0D+0)
      INTRINSIC ABS,MIN
      EXTERNAL FD15AD,MC64DD,MC64ED,MC64FD
      DOUBLE PRECISION FD15AD
      RINF = FD15AD('H')
      NUM = 0
      BV = RINF
      DO 10 K = 1,N
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        D(K) = ZERO
   10 CONTINUE
      DO 20 J = 1,N
        A0 = MINONE
        DO 30 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.GT.D(I)) D(I) = AI
          IF (JPERM(J).NE.0) GO TO 30
          IF (AI.GE.BV) THEN
            A0 = BV
            IF (IPERM(I).NE.0) GO TO 30
            JPERM(J) = I
            IPERM(I) = J
            NUM = NUM + 1
          ELSE
            IF (AI.LE.A0) GO TO 30
            A0 = AI
            I0 = I
          ENDIF
   30   CONTINUE
        IF (A0.NE.MINONE .AND. A0.LT.BV) THEN
          BV = A0
          IF (IPERM(I0).NE.0) GO TO 20
          IPERM(I0) = J
          JPERM(J) = I0
          NUM = NUM + 1
        ENDIF
   20 CONTINUE
      DO 25 I = 1,N
        BV = MIN(BV,D(I))
   25 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        DO 50 K = IP(J),IP(J+1)-1
          I = IRN(K)
          AI = ABS(A(K))
          IF (AI.LT.BV) GO TO 50
          IF (IPERM(I).EQ.0) GO TO 90
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 50
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).NE.0) GO TO 70
            IF (ABS(A(KK)).GE.BV) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   50   CONTINUE
        GO TO 95
   80   JPERM(JJ) = II
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = I
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,N
        D(I) = MINONE
        L(I) = 0
   99 CONTINUE
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        QLEN = 0
        LOW = N + 1
        UP = N + 1
        CSP = MINONE
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = ABS(A(K))
          IF (CSP.GE.DNEW) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = I
            JSP = J
            IF (CSP.GE.BV) GO TO 160
          ELSE
            D(I) = DNEW
            IF (DNEW.GE.BV) THEN
              LOW = LOW - 1
              Q(LOW) = I
            ELSE
              QLEN = QLEN + 1
              L(I) = QLEN
              CALL MC64DD(I,N,Q,D,L,1)
            ENDIF
            JJ = IPERM(I)
            PR(JJ) = J
          ENDIF
  115   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (CSP.GE.D(I)) GO TO 160
            BV = D(I)
            DO 152 IDUM = 1,N
              CALL MC64ED(QLEN,N,Q,D,L,1)
              L(I) = 0
              LOW = LOW - 1
              Q(LOW) = I
              IF (QLEN.EQ.0) GO TO 153
              I = Q(1)
              IF (D(I).NE.BV) GO TO 153
  152       CONTINUE
          ENDIF
  153     UP = UP - 1
          Q0 = Q(UP)
          DQ0 = D(Q0)
          L(Q0) = UP
          J = IPERM(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
            DNEW = MIN(DQ0,ABS(A(K)))
            IF (CSP.GE.DNEW) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = I
              JSP = J
              IF (CSP.GE.BV) GO TO 160
            ELSE
              DI = D(I)
              IF (DI.GE.BV .OR. DI.GE.DNEW) GO TO 155
              D(I) = DNEW
              IF (DNEW.GE.BV) THEN
                IF (DI.NE.MINONE) THEN
                  LPOS = L(I)
                  CALL MC64FD(LPOS,QLEN,N,Q,D,L,1)
                ENDIF
                L(I) = 0
                LOW = LOW - 1
                Q(LOW) = I
              ELSE
                IF (DI.EQ.MINONE) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,1)
              ENDIF
              JJ = IPERM(I)
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.MINONE) GO TO 190
        BV = MIN(BV,CSP)
        NUM = NUM + 1
        I = ISP
        J = JSP
        DO 170 JDUM = 1,NUM+1
          I0 = JPERM(J)
          JPERM(J) = I
          IPERM(I) = J
          J = PR(J)
          IF (J.EQ.-1) GO TO 190
          I = I0
  170   CONTINUE
  190   DO 191 KK = UP,N
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  191   CONTINUE
        DO 192 KK = LOW,UP-1
          I = Q(KK)
          D(I) = MINONE
  192   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = MINONE
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          PR(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 I = 1,N
        IF (JPERM(I).NE.0) GO TO 320
        K = K + 1
        JDUM = PR(K)
        IPERM(JDUM) = - I
  320 CONTINUE
 1000 RETURN
      END
C**********************************************************************
      SUBROUTINE MC64DD(I,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER I,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DI
      POS = L(I)
      IF (POS.LE.1) GO TO 20
      DI = D(I)
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64ED(QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = 1
      IF (IWAY.EQ.1) THEN
        DO 10 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   10   CONTINUE
      ELSE
        DO 15 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 20
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 20
          Q(POS) = Q(POSK)
          L(Q(POS)) = POS
          POS = POSK
   15   CONTINUE
      ENDIF
   20 Q(POS) = I
      L(I) = POS
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64FD(POS0,QLEN,N,Q,D,L,IWAY)
      IMPLICIT NONE
      INTEGER POS0,QLEN,N,IWAY
      INTEGER Q(N),L(N)
      DOUBLE PRECISION D(N)
      INTEGER I,IDUM,K,POS,POSK,QK
      PARAMETER (K=2)
      DOUBLE PRECISION DK,DR,DI
      IF (QLEN.EQ.POS0) THEN
        QLEN = QLEN - 1
        RETURN
      ENDIF
      I = Q(QLEN)
      DI = D(I)
      QLEN = QLEN - 1
      POS = POS0
      IF (IWAY.EQ.1) THEN
        IF (POS.LE.1) GO TO 20
        DO 10 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.LE.D(QK)) GO TO 20
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 20
   10   CONTINUE
   20   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 30 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.LT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.GE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   30   CONTINUE
      ELSE
        IF (POS.LE.1) GO TO 34
        DO 32 IDUM = 1,N
          POSK = POS/K
          QK = Q(POSK)
          IF (DI.GE.D(QK)) GO TO 34
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
          IF (POS.LE.1) GO TO 34
   32   CONTINUE
   34   Q(POS) = I
        L(I) = POS
        IF (POS.NE.POS0) RETURN
        DO 36 IDUM = 1,N
          POSK = K*POS
          IF (POSK.GT.QLEN) GO TO 40
          DK = D(Q(POSK))
          IF (POSK.LT.QLEN) THEN
            DR = D(Q(POSK+1))
            IF (DK.GT.DR) THEN
              POSK = POSK + 1
              DK = DR
            ENDIF
          ENDIF
          IF (DI.LE.DK) GO TO 40
          QK = Q(POSK)
          Q(POS) = QK
          L(QK) = POS
          POS = POSK
   36   CONTINUE
      ENDIF
   40 Q(POS) = I
      L(I) = POS
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64RD(N,NE,IP,IRN,A)
      IMPLICIT NONE
      INTEGER N,NE
      INTEGER IP(N+1),IRN(NE)
      DOUBLE PRECISION A(NE)
      INTEGER THRESH,TDLEN
      PARAMETER (THRESH=15,TDLEN=50)
      INTEGER J,IPJ,K,LEN,R,S,HI,FIRST,MID,LAST,TD
      DOUBLE PRECISION HA,KEY
      INTEGER TODO(TDLEN)
      DO 100 J = 1,N
        LEN = IP(J+1) - IP(J)
        IF (LEN.LE.1) GO TO 100
        IPJ = IP(J)
        IF (LEN.LT.THRESH) GO TO 400
        TODO(1) = IPJ
        TODO(2) = IPJ + LEN
        TD = 2
  500   CONTINUE
        FIRST = TODO(TD-1)
        LAST = TODO(TD)
        KEY = A((FIRST+LAST)/2)
        DO 475 K = FIRST,LAST-1
          HA = A(K)
          IF (HA.EQ.KEY) GO TO 475
          IF (HA.GT.KEY) GO TO 470
          KEY = HA
          GO TO 470
  475   CONTINUE
        TD = TD - 2
        GO TO 425
  470   MID = FIRST
        DO 450 K = FIRST,LAST-1
          IF (A(K).LE.KEY) GO TO 450
          HA = A(MID)
          A(MID) = A(K)
          A(K) = HA
          HI = IRN(MID)
          IRN(MID) = IRN(K)
          IRN(K) = HI
          MID = MID + 1
  450   CONTINUE
        IF (MID-FIRST.GE.LAST-MID) THEN
          TODO(TD+2) = LAST
          TODO(TD+1) = MID
          TODO(TD) = MID
        ELSE
          TODO(TD+2) = MID
          TODO(TD+1) = FIRST
          TODO(TD) = LAST
          TODO(TD-1) = MID
        ENDIF
        TD = TD + 2
  425   CONTINUE
        IF (TD.EQ.0) GO TO 400
        IF (TODO(TD)-TODO(TD-1).GE.THRESH) GO TO 500
        TD = TD - 2
        GO TO 425
  400   DO 200 R = IPJ+1,IPJ+LEN-1
          IF (A(R-1) .LT. A(R)) THEN
            HA = A(R)
            HI = IRN(R)
            A(R) = A(R-1)
            IRN(R) = IRN(R-1)
            DO 300 S = R-1,IPJ+1,-1
              IF (A(S-1) .LT. HA) THEN
                A(S) = A(S-1)
                IRN(S) = IRN(S-1)
              ELSE
                A(S) = HA
                IRN(S) = HI
                GO TO 200
              END IF
  300       CONTINUE
            A(IPJ) = HA
            IRN(IPJ) = HI
          END IF
  200   CONTINUE
  100 CONTINUE
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64SD(N,NE,IP,IRN,A,IPERM,NUMX,
     &           W,LEN,LENL,LENH,FC,IW,IW4)
      IMPLICIT NONE
      INTEGER N,NE,NUMX
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        W(N),LEN(N),LENL(N),LENH(N),FC(N),IW(N),IW4(4*N)
      DOUBLE PRECISION A(NE)
      INTEGER NUM,NVAL,WLEN,II,I,J,K,L,CNT,MOD,IDUM1,IDUM2,IDUM3
      DOUBLE PRECISION BVAL,BMIN,BMAX,RINF
      EXTERNAL FD15AD,MC64QD,MC64UD
      DOUBLE PRECISION FD15AD
      RINF = FD15AD('H')
      DO 20 J = 1,N
        FC(J) = J
        IW(J) = 0
        LEN(J) = IP(J+1) - IP(J)
   20 CONTINUE
      CNT = 1
      MOD = 1
      NUMX = 0
      CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUMX,N,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))
      NUM = NUMX
      IF (NUM.NE.N) THEN
        BMAX = RINF
      ELSE
        BMAX = RINF
        DO 30 J = 1,N
          BVAL = 0.0
          DO 25 K = IP(J),IP(J+1)-1
            IF (A(K).GT.BVAL) BVAL = A(K)
   25     CONTINUE
          IF (BVAL.LT.BMAX) BMAX = BVAL
   30   CONTINUE
        BMAX = 1.001 * BMAX
      ENDIF
      BVAL = 0.0
      BMIN = 0.0
      WLEN = 0
      DO 48 J = 1,N
        L = IP(J+1) - IP(J)
        LENH(J) = L
        LEN(J) = L
        DO 45 K = IP(J),IP(J+1)-1
          IF (A(K).LT.BMAX) GO TO 46
   45   CONTINUE
        K = IP(J+1)
   46   LENL(J) = K - IP(J)
        IF (LENL(J).EQ.L) GO TO 48
        WLEN = WLEN + 1
        W(WLEN) = J
   48 CONTINUE
      DO 90 IDUM1 = 1,NE
        IF (NUM.EQ.NUMX) THEN
          DO 50 I = 1,N
            IPERM(I) = IW(I)
   50     CONTINUE
          DO 80 IDUM2 = 1,NE
            BMIN = BVAL
            IF (BMAX .EQ. BMIN) GO TO 99
            CALL MC64QD(IP,LENL,LEN,W,WLEN,A,NVAL,BVAL)
            IF (NVAL.LE.1) GO TO 99
            K = 1
            DO 70 IDUM3 = 1,N
              IF (K.GT.WLEN) GO TO 71
              J = W(K)
              DO 55 II = IP(J)+LEN(J)-1,IP(J)+LENL(J),-1
                IF (A(II).GE.BVAL) GO TO 60
                I = IRN(II)
                IF (IW(I).NE.J) GO TO 55
                IW(I) = 0
                NUM = NUM - 1
                FC(N-NUM) = J
   55         CONTINUE
   60         LENH(J) = LEN(J)
              LEN(J) = II - IP(J) + 1
              IF (LENL(J).EQ.LENH(J)) THEN
                W(K) = W(WLEN)
                WLEN = WLEN - 1
              ELSE
                K = K + 1
              ENDIF
   70       CONTINUE
   71       IF (NUM.LT.NUMX) GO TO 81
   80     CONTINUE
   81     MOD = 1
        ELSE
          BMAX = BVAL
          CALL MC64QD(IP,LEN,LENH,W,WLEN,A,NVAL,BVAL)
          IF (NVAL.EQ.0. OR. BVAL.EQ.BMIN) GO TO 99
          K = 1
          DO 87 IDUM3 = 1,N
            IF (K.GT.WLEN) GO TO 88
            J = W(K)
            DO 85 II = IP(J)+LEN(J),IP(J)+LENH(J)-1
              IF (A(II).LT.BVAL) GO TO 86
   85       CONTINUE
   86       LENL(J) = LEN(J)
            LEN(J) = II - IP(J)
            IF (LENL(J).EQ.LENH(J)) THEN
              W(K) = W(WLEN)
              WLEN = WLEN - 1
            ELSE
              K = K + 1
            ENDIF
   87     CONTINUE
   88     MOD = 0
        ENDIF
        CNT = CNT + 1
        CALL MC64UD(CNT,MOD,N,IRN,NE,IP,LEN,FC,IW,NUM,NUMX,
     &            IW4(1),IW4(N+1),IW4(2*N+1),IW4(3*N+1))
   90 CONTINUE
   99 IF (NUMX.EQ.N) GO TO 1000
      DO 300 J = 1,N
        W(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          IW(K) = I
        ELSE
          J = IPERM(I)
          W(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (W(J).NE.0) GO TO 320
        K = K + 1
        IDUM1 = IW(K)
        IPERM(IDUM1) =  - J
  320 CONTINUE
 1000 RETURN
      END
C**********************************************************************
      SUBROUTINE MC64QD(IP,LENL,LENH,W,WLEN,A,NVAL,VAL)
      IMPLICIT NONE
      INTEGER WLEN,NVAL
      INTEGER IP(*),LENL(*),LENH(*),W(*)
      DOUBLE PRECISION A(*),VAL
      INTEGER XX,J,K,II,S,POS
      PARAMETER (XX=10)
      DOUBLE PRECISION SPLIT(XX),HA
      NVAL = 0
      DO 10 K = 1,WLEN
        J = W(K)
        DO 15 II = IP(J)+LENL(J),IP(J)+LENH(J)-1
          HA = A(II)
          IF (NVAL.EQ.0) THEN
            SPLIT(1) = HA
            NVAL = 1
          ELSE
            DO 20 S = NVAL,1,-1
              IF (SPLIT(S).EQ.HA) GO TO 15
              IF (SPLIT(S).GT.HA) THEN
                POS = S + 1
                GO TO 21
              ENDIF
  20        CONTINUE
            POS = 1
  21        DO 22 S = NVAL,POS,-1
              SPLIT(S+1) = SPLIT(S)
  22        CONTINUE
            SPLIT(POS) = HA
            NVAL = NVAL + 1
          ENDIF
          IF (NVAL.EQ.XX) GO TO 11
  15    CONTINUE
  10  CONTINUE
  11  IF (NVAL.GT.0) VAL = SPLIT((NVAL+1)/2)
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64UD(ID,MOD,N,IRN,LIRN,IP,LENC,FC,IPERM,NUM,NUMX,
     &           PR,ARP,CV,OUT)
      IMPLICIT NONE
      INTEGER ID,MOD,N,LIRN,NUM,NUMX
      INTEGER ARP(N),CV(N),IRN(LIRN),IP(N),
     &        FC(N),IPERM(N),LENC(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,J,J1,JORD,K,KK,LAST,NFC,
     &        NUM0,NUM1,NUM2,ID0,ID1
      IF (ID.EQ.1) THEN
        DO 5 I = 1,N
          CV(I) = 0
          ARP(I) = 0
    5   CONTINUE
        NUM1 = N
        NUM2 = N
      ELSE
        IF (MOD.EQ.1) THEN
          DO 8 I = 1,N
            ARP(I) = 0
    8     CONTINUE
        ENDIF
        NUM1 = NUMX
        NUM2 = N - NUMX
      ENDIF
      NUM0 = NUM
      NFC = 0
      ID0 = (ID-1)*N
      DO 100 JORD = NUM0+1,N
        ID1 = ID0 + JORD
        J = FC(JORD-NUM0)
        PR(J) = -1
        DO 70 K = 1,JORD
          IF (ARP(J).GE.LENC(J)) GO TO 30
          IN1 = IP(J) + ARP(J)
          IN2 = IP(J) + LENC(J) - 1
          DO 20 II = IN1,IN2
            I = IRN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = LENC(J)
   30     OUT(J) = LENC(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENC(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = IRN(II)
              IF (CV(I).EQ.ID1) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = ID1
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       J1 = PR(J)
            IF (J1.EQ.-1) THEN
              NFC = NFC + 1
              FC(NFC) = J
              IF (NFC.GT.NUM2) THEN
                LAST = JORD
                GO TO 101
              ENDIF
              GO TO 100
            ENDIF
            J = J1
   60     CONTINUE
   70   CONTINUE
   80   IPERM(I) = J
        ARP(J) = II - IP(J) + 1
        NUM = NUM + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 95
          II = IP(J) + LENC(J) - OUT(J) - 2
          I = IRN(II)
          IPERM(I) = J
   90   CONTINUE
   95   IF (NUM.EQ.NUM1) THEN
          LAST = JORD
          GO TO 101
        ENDIF
  100 CONTINUE
      LAST = N
  101 DO 110 JORD = LAST+1,N
        NFC = NFC + 1
        FC(NFC) = FC(JORD-NUM0)
  110 CONTINUE
      RETURN
      END
C**********************************************************************
      SUBROUTINE MC64WD(N,NE,IP,IRN,A,IPERM,NUM,
     &           JPERM,OUT,PR,Q,L,U,D)
      IMPLICIT NONE
      INTEGER N,NE,NUM
      INTEGER IP(N+1),IRN(NE),IPERM(N),
     &        JPERM(N),OUT(N),PR(N),Q(N),L(N)
      DOUBLE PRECISION A(NE),U(N),D(N)
      INTEGER I,I0,II,J,JJ,JORD,Q0,QLEN,JDUM,ISP,JSP,
     &        K,K0,K1,K2,KK,KK1,KK2,UP,LOW,LPOS
      DOUBLE PRECISION CSP,DI,DMIN,DNEW,DQ0,VJ
      DOUBLE PRECISION RINF,ZERO
      PARAMETER (ZERO=0.0D+0)
      EXTERNAL FD15AD,MC64DD,MC64ED,MC64FD
      DOUBLE PRECISION FD15AD
      RINF = FD15AD('H')
      NUM = 0
      DO 10 K = 1,N
        U(K) = RINF
        D(K) = ZERO
        IPERM(K) = 0
        JPERM(K) = 0
        PR(K) = IP(K)
        L(K) = 0
   10 CONTINUE
      DO 30 J = 1,N
        DO 20 K = IP(J),IP(J+1)-1
          I = IRN(K)
          IF (A(K).GT.U(I)) GO TO 20
          U(I) = A(K)
          IPERM(I) = J
          L(I) = K
   20   CONTINUE
   30 CONTINUE
      DO 40 I = 1,N
        J = IPERM(I)
        IF (J.EQ.0) GO TO 40
        IPERM(I) = 0
        IF (JPERM(J).NE.0) GO TO 40
        IF (IP(J+1)-IP(J) .GT. N/10 .AND. N.GT.50) GO TO 40
        NUM = NUM + 1
        IPERM(I) = J
        JPERM(J) = L(I)
   40 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 95 J = 1,N
        IF (JPERM(J).NE.0) GO TO 95
        K1 = IP(J)
        K2 = IP(J+1) - 1
        IF (K1.GT.K2) GO TO 95
        I0 = IRN(K1)
        VJ = A(K1) - U(I0)
        K0 = K1
        DO 50 K = K1+1,K2
          I = IRN(K)
          DI = A(K) - U(I)
          IF (DI.GT.VJ) GO TO 50
          IF (DI.LT.VJ .OR. DI.EQ.RINF) GO TO 55
          IF (IPERM(I).NE.0 .OR. IPERM(I0).EQ.0) GO TO 50
   55     VJ = DI
          I0 = I
          K0 = K
   50   CONTINUE
        D(J) = VJ
        K = K0
        I = I0
        IF (IPERM(I).EQ.0) GO TO 90
        DO 60 K = K0,K2
          I = IRN(K)
          IF (A(K)-U(I).GT.VJ) GO TO 60
          JJ = IPERM(I)
          KK1 = PR(JJ)
          KK2 = IP(JJ+1) - 1
          IF (KK1.GT.KK2) GO TO 60
          DO 70 KK = KK1,KK2
            II = IRN(KK)
            IF (IPERM(II).GT.0) GO TO 70
            IF (A(KK)-U(II).LE.D(JJ)) GO TO 80
   70     CONTINUE
          PR(JJ) = KK2 + 1
   60   CONTINUE
        GO TO 95
   80   JPERM(JJ) = KK
        IPERM(II) = JJ
        PR(JJ) = KK + 1
   90   NUM = NUM + 1
        JPERM(J) = K
        IPERM(I) = J
        PR(J) = K + 1
   95 CONTINUE
      IF (NUM.EQ.N) GO TO 1000
      DO 99 I = 1,N
        D(I) = RINF
        L(I) = 0
   99 CONTINUE
      DO 100 JORD = 1,N
        IF (JPERM(JORD).NE.0) GO TO 100
        DMIN = RINF
        QLEN = 0
        LOW = N + 1
        UP = N + 1
        CSP = RINF
        J = JORD
        PR(J) = -1
        DO 115 K = IP(J),IP(J+1)-1
          I = IRN(K)
          DNEW = A(K) - U(I)
          IF (DNEW.GE.CSP) GO TO 115
          IF (IPERM(I).EQ.0) THEN
            CSP = DNEW
            ISP = K
            JSP = J
          ELSE
            IF (DNEW.LT.DMIN) DMIN = DNEW
            D(I) = DNEW
            QLEN = QLEN + 1
            Q(QLEN) = K
          ENDIF
  115   CONTINUE
        Q0 = QLEN
        QLEN = 0
        DO 120 KK = 1,Q0
          K = Q(KK)
          I = IRN(K)
          IF (CSP.LE.D(I)) THEN
            D(I) = RINF
            GO TO 120
          ENDIF
          IF (D(I).LE.DMIN) THEN
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
          ELSE
            QLEN = QLEN + 1
            L(I) = QLEN
            CALL MC64DD(I,N,Q,D,L,2)
          ENDIF
          JJ = IPERM(I)
          OUT(JJ) = K
          PR(JJ) = J
  120   CONTINUE
        DO 150 JDUM = 1,NUM
          IF (LOW.EQ.UP) THEN
            IF (QLEN.EQ.0) GO TO 160
            I = Q(1)
            IF (D(I).GE.CSP) GO TO 160
            DMIN = D(I)
  152       CALL MC64ED(QLEN,N,Q,D,L,2)
            LOW = LOW - 1
            Q(LOW) = I
            L(I) = LOW
            IF (QLEN.EQ.0) GO TO 153
            I = Q(1)
            IF (D(I).GT.DMIN) GO TO 153
            GO TO 152
          ENDIF
  153     Q0 = Q(UP-1)
          DQ0 = D(Q0)
          IF (DQ0.GE.CSP) GO TO 160
          UP = UP - 1
          J = IPERM(Q0)
          VJ = DQ0 - A(JPERM(J)) + U(Q0)
          DO 155 K = IP(J),IP(J+1)-1
            I = IRN(K)
            IF (L(I).GE.UP) GO TO 155
            DNEW = VJ + A(K)-U(I)
            IF (DNEW.GE.CSP) GO TO 155
            IF (IPERM(I).EQ.0) THEN
              CSP = DNEW
              ISP = K
              JSP = J
            ELSE
              DI = D(I)
              IF (DI.LE.DNEW) GO TO 155
              IF (L(I).GE.LOW) GO TO 155
              D(I) = DNEW
              IF (DNEW.LE.DMIN) THEN
                LPOS = L(I)
                IF (LPOS.NE.0)
     *            CALL MC64FD(LPOS,QLEN,N,Q,D,L,2)
                LOW = LOW - 1
                Q(LOW) = I
                L(I) = LOW
              ELSE
                IF (L(I).EQ.0) THEN
                  QLEN = QLEN + 1
                  L(I) = QLEN
                ENDIF
                CALL MC64DD(I,N,Q,D,L,2)
              ENDIF
              JJ = IPERM(I)
              OUT(JJ) = K
              PR(JJ) = J
            ENDIF
  155     CONTINUE
  150   CONTINUE
  160   IF (CSP.EQ.RINF) GO TO 190
        NUM = NUM + 1
        I = IRN(ISP)
        IPERM(I) = JSP
        JPERM(JSP) = ISP
        J = JSP
        DO 170 JDUM = 1,NUM
          JJ = PR(J)
          IF (JJ.EQ.-1) GO TO 180
          K = OUT(J)
          I = IRN(K)
          IPERM(I) = JJ
          JPERM(JJ) = K
          J = JJ
  170   CONTINUE
  180   DO 185 KK = UP,N
          I = Q(KK)
          U(I) = U(I) + D(I) - CSP
  185   CONTINUE
  190   DO 191 KK = LOW,N
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  191   CONTINUE
        DO 193 KK = 1,QLEN
          I = Q(KK)
          D(I) = RINF
          L(I) = 0
  193   CONTINUE
  100 CONTINUE
 1000 DO 200 J = 1,N
        K = JPERM(J)
        IF (K.NE.0) THEN
          D(J) = A(K) - U(IRN(K))
        ELSE
          D(J) = ZERO
        ENDIF
        IF (IPERM(J).EQ.0) U(J) = ZERO
  200 CONTINUE
      IF (NUM.EQ.N) GO TO 1100
      DO 300 J = 1,N
        JPERM(J) = 0
  300 CONTINUE
      K = 0
      DO 310 I = 1,N
        IF (IPERM(I).EQ.0) THEN
          K = K + 1
          OUT(K) = I
        ELSE
          J = IPERM(I)
          JPERM(J) = I
        ENDIF
  310 CONTINUE
      K = 0
      DO 320 J = 1,N
        IF (JPERM(J).NE.0) GO TO 320
        K = K + 1
        JDUM = OUT(K)
        IPERM(JDUM) = - J
  320 CONTINUE
 1100 RETURN
      END
* *******************************************************************
* COPYRIGHT (c) 1977 Hyprotech UK
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
* Original date 8 Oct 1992
C######8/10/92 Toolpack tool decs employed.
C######8/10/92 D version created by name change only.
C 13/3/02 Cosmetic changes applied to reduce single/double differences
C
C 12th July 2004 Version 1.0.0. Version numbering added.

      SUBROUTINE MC21AD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
      INTEGER LICN,N,NUMNZ
      INTEGER ICN(LICN),IP(N),IPERM(N),IW(N,4),LENR(N)
      EXTERNAL MC21BD
      CALL MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),
     +            IW(1,3),IW(1,4))
      RETURN
      END
      SUBROUTINE MC21BD(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
      INTEGER LICN,N,NUMNZ
      INTEGER ARP(N),CV(N),ICN(LICN),IP(N),IPERM(N),LENR(N),OUT(N),PR(N)
      INTEGER I,II,IN1,IN2,IOUTK,J,J1,JORD,K,KK
      DO 10 I = 1,N
        ARP(I) = LENR(I) - 1
        CV(I) = 0
        IPERM(I) = 0
   10 CONTINUE
      NUMNZ = 0
      DO 100 JORD = 1,N
        J = JORD
        PR(J) = -1
        DO 70 K = 1,JORD
          IN1 = ARP(J)
          IF (IN1.LT.0) GO TO 30
          IN2 = IP(J) + LENR(J) - 1
          IN1 = IN2 - IN1
          DO 20 II = IN1,IN2
            I = ICN(II)
            IF (IPERM(I).EQ.0) GO TO 80
   20     CONTINUE
          ARP(J) = -1
   30     CONTINUE
          OUT(J) = LENR(J) - 1
          DO 60 KK = 1,JORD
            IN1 = OUT(J)
            IF (IN1.LT.0) GO TO 50
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
            DO 40 II = IN1,IN2
              I = ICN(II)
              IF (CV(I).EQ.JORD) GO TO 40
              J1 = J
              J = IPERM(I)
              CV(I) = JORD
              PR(J) = J1
              OUT(J1) = IN2 - II - 1
              GO TO 70
   40       CONTINUE
   50       CONTINUE
            J = PR(J)
            IF (J.EQ.-1) GO TO 100
   60     CONTINUE
   70   CONTINUE
   80   CONTINUE
        IPERM(I) = J
        ARP(J) = IN2 - II - 1
        NUMNZ = NUMNZ + 1
        DO 90 K = 1,JORD
          J = PR(J)
          IF (J.EQ.-1) GO TO 100
          II = IP(J) + LENR(J) - OUT(J) - 2
          I = ICN(II)
          IPERM(I) = J
   90   CONTINUE
  100 CONTINUE
      IF (NUMNZ.EQ.N) RETURN
      DO 110 I = 1,N
        ARP(I) = 0
  110 CONTINUE
      K = 0
      DO 130 I = 1,N
        IF (IPERM(I).NE.0) GO TO 120
        K = K + 1
        OUT(K) = I
        GO TO 130
  120   CONTINUE
        J = IPERM(I)
        ARP(J) = I
  130 CONTINUE
      K = 0
      DO 140 I = 1,N
        IF (ARP(I).NE.0) GO TO 140
        K = K + 1
        IOUTK = OUT(K)
        IPERM(IOUTK) = I
  140 CONTINUE
      RETURN
      END
* COPYRIGHT (c) 1988 AEA Technology
* Original date 17 Feb 2005

C 17th February 2005 Version 1.0.0. Replacement for FD05.

      DOUBLE PRECISION FUNCTION FD15AD(T)
C----------------------------------------------------------------
C  Fortran 77 implementation of the Fortran 90 intrinsic
C    functions: EPSILON, TINY, HUGE and RADIX.  Note that
C    the RADIX result is returned as DOUBLE PRECISION.
C
C  The CHARACTER argument specifies the type of result:
C       
C   'E'  smallest positive real number: 1.0 + DC(1) > 1.0, i.e.
C          EPSILON(DOUBLE PRECISION)
C   'T'  smallest full precision positive real number, i.e.
C          TINY(DOUBLE PRECISION)
C   'H'  largest finite positive real number, i.e.
C          HUGE(DOUBLE PRECISION)
C   'R'  the base of the floating point arithematic, i.e.
C          RADIX(DOUBLE PRECISION)
C
C    any other value gives a result of zero.
C----------------------------------------------------------------
      CHARACTER T

      IF ( T.EQ.'E' ) THEN
         FD15AD = EPSILON(1.0D0)
      ELSE IF ( T.EQ.'T' ) THEN
         FD15AD = TINY(1.0D0)
      ELSE IF ( T.EQ.'H' ) THEN
         FD15AD = HUGE(1.0D0)
      ELSE IF ( T.EQ.'R' ) THEN
         FD15AD = DBLE(RADIX(1.0D0))
      ELSE
         FD15AD = 0.0D0
      ENDIF
      RETURN
      END

