* *******************************************************************
* COPYRIGHT (c) 2000 Council for the Central Laboratory
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
* Original date 30 March 2000

c JAS: edited to conform to pfort
c JAS: commented changes
C 20/3/02 Cosmetic changes applied to reduce single/double differences
c Also comments standardised (second column blank) so that they
c are stripped effectively for HSL distribution.
c   LAPACK : DSBEV, DLAMCH, DLASET, DLACPY, DLARTG
c            DLARFG, DRSCL, DLARNV, DGESVD, DSYEV
c   BLAS   : DNRM2, IDAMAX, DSCAL, DCOPY, DGEMM, DAXPY, DGEMV,
c            DGER, DTBSV
c   HSL    : KB07AD, KB08AD

c 12th July 2004 Version 1.0.0. Version numbering added.
c 1st March 2005 Version 1.1.0. ZA02 dependence changed to ZA12.
c 27 February 2008 Version 1.1.1. Comments flowed to column 72.
c 9 April 2008 Version 1.2.1. EXTERNAL ZA12AD corrected in sl version.

      subroutine EA16AD(N, BLK, NWANT, NV, LIWORK, LWORK, ICNTL, INFO)
      integer N, NV, BLK, NWANT
      integer LIWORK, LWORK
      integer ICNTL(20), INFO(20)
      integer ISAVE(1)
      external EA16GD
      call EA16GD(N, NV, NWANT, BLK, LIWORK, LWORK, ISAVE, 1, ICNTL,
     &            INFO)
      return
      end
      subroutine EA16BD(N, BLK, NWANT, NV, MODE, WHICH, IDO, IPOS,
     &                  V, LDV, BV, LDBV, RANGE, SIGMA, NEINEG, IWORK,
     &                  LIWORK, WORK, LWORK, ICNTL, CNTL, INFO)
      integer IDO, N, NV, BLK, MODE, WHICH, NWANT, LDV, LDBV
      integer NEINEG, LIWORK, LWORK
      integer IPOS(10), IWORK(LIWORK), ICNTL(20), INFO(20)
      double precision SIGMA
      double precision V(LDV,NV), BV(LDBV,BLK), CNTL(15)
      double precision WORK(LWORK), RANGE(2)
      integer INFO8, NLEJA2, NDONE, IDOGS2, NCOKE2, MODE2
      integer IDUM1, IDUM2, IDUM3, IISTE2, LURES2, LUWAR2, LUERR2
      integer NCONV2, NCOOL2, NMVA2, NMVB2, VINIT2, NWACO2, HARMON
      integer I, DIFF, LU, SHIFT, NSHIFT, NCOLEF, NCORIG, LDERV2
      INTEGER IPOS7
      double precision NORMT, GAMMA2, OLDSIG, LARGE
      logical KEEP, FIRST, ISINF, CVGED, INFLAN, INFLOC, EXPRST
      double precision ZERO, ONE
      parameter       (ZERO=0.0D0, ONE=1.0D0)
      integer IERR, MEMO, NMVA, NMVB, NRESTA, NPOLE
      parameter (IERR=1, MEMO=2, NMVA=3, NMVB=4, NRESTA=5, NPOLE=6)
      integer INCONV
      parameter (INCONV=7)
      integer TOLLAN, TOL, TOLLOC, PLGROW, RATNOR
      parameter (TOLLAN=1, TOLLOC=1, TOL=2, PLGROW=4, RATNOR=5)
      integer DIAGNO, VINIT, CNWPOL, ISEED, MXNMVA
      integer IRSTRT, IPREFI, CHAMOD, PRIRES, POLTIM, POLRTR
      parameter (DIAGNO=2, VINIT=3, CNWPOL=5, ISEED=10,
     &           MXNMVA=6, IRSTRT=9, IPREFI = 15,
     &           CHAMOD=7,PRIRES=16, POLTIM=17, POLRTR=18)
      integer PEPS, PSAFM, PERRM, PEPSS, PRANGE, PEA16E, PINTER
      integer TRUST, GAMMA, RKSGRO, GAMGRO, PEA17N, PNORMT, TIMFA1
      integer TIMFA2
      parameter (TRUST=1, PEPS=5, PSAFM=6, PERRM=7, PEPSS=8, PRANGE=9,
     &           PINTER=11, PEA17N=15, PEA16E=25, GAMMA=29,
     &           RKSGRO=30, GAMGRO=31)
      parameter (PNORMT=32, TIMFA1=33, TIMFA2=34)
      integer SEED, JUMP, LUERR, LUWARN, LUREST
      parameter (SEED=1, JUMP=5, LUERR=6, LUWARN=7, LUREST=9)
      integer NWACO, ISETPL, IBGS, IF, ILDF, ILCF
      parameter (NWACO=10, ISETPL=11, IBGS=12, IF=13, ILDF=14, ILCF=15)
      integer ILCF2, ITEMP3, ILTMP3, IL, ILDL
      parameter (ILCF2=16, ITEMP3=17, ILTMP3=18, IL=19, ILDL=20)
      integer ILCL, ITAU, ILTAU, ILTMP1, ITEMP1
      parameter (ILCL=21, ITAU=22, ILTAU=23, ILTMP1=24, ITEMP1=25)
      integer IRANDM, GSFAIL, ILDT, IT, INORM
      parameter (IRANDM=26, GSFAIL=27, ILDT=28, IT=29, INORM=30)
      integer ICOEF, ILDCOE, IPIPE, ITEMP2, INPURG
      parameter (ICOEF=31, ILDCOE=32, IPIPE=33, ITEMP2=34, INPURG=35)
      integer IEIGEN, IRESID, IRITZ, IEIGVE
      parameter (IEIGEN=36, IRESID=38, IRITZ=39, IEIGVE=40)
      integer ILDEIG, ILCEIG, INCOOL, ILSTEP, LLSTEP
      parameter (ILDEIG=41, ILCEIG=42, INCOOL=43, ILSTEP=44, LLSTEP=45)
      integer ERMOVE, LDERVE, IRMOLA, ILERLA
      parameter (ERMOVE=46, LDERVE=47, IRMOLA=48, ILERLA=49)
      integer IISTEP, IMLANC, INKEEP, NBACK
      parameter (IISTEP=50, IMLANC=51, INKEEP=52, NBACK=53)
      integer INSTEP, IKKRYL, INLANC, ILLANC
      parameter (INSTEP=54, IKKRYL=55, INLANC=56, ILLANC=57)
      integer GSSTRT, GSNUMB, NCOKEE, NWPOLE
      parameter (GSSTRT=58, GSNUMB=59, NCOKEE=60, NWPOLE=61)
      integer IDOGS, IVINIT, IOFFST, INEIGV, INPREF
      parameter (IDOGS=62, IVINIT=63, IOFFST=64, INEIGV=65, INPREF=66)
      integer ILEJA, NLEJA, RLEJA
      parameter (ILEJA=67, NLEJA=68, RLEJA=69)
      integer IMV, IEA16E, MSTEP, INDONE, IEA17N
      parameter (IMV=70, IEA16E=71, MSTEP=87, INDONE=88, IEA17N=89)
      integer IMODE, INSHIF, INCOLE, INCORI
      parameter (IMODE=95, INSHIF=96, INCOLE=97, INCORI=98)
      integer ICVG, LMLEJA, IWHRES, IFILT, IHARM
      parameter (ICVG=99, LMLEJA=100, IWHRES=101, IFILT=102, IHARM=103)
      integer RANDOM, LDT, T, NORM, LTEMP1, TEMP1
      integer COEF, LDCOEF, PIPE, TEMP2, EIGEN
      integer RESID, RITZ, EIGVEC, LDEIGV, LCEIGV
      integer NCOOLD, NPURG, NEIGV
      integer F, LDF, LCF, LCF2, TEMP3, LTEMP3, L, LDL, LCL, TAU, LTAU
      integer LSTEP, MLEJA, LLEJA
      integer ISTEP, NKEEP
      integer NSTEP, KKRYL, NLANC, LLANC, MLANC, NNSTEP
      double precision SAFMIN, EPS, EPSS, ERRMAX
      integer NCONV, OFFSET
      logical BGS, SETPOL
      integer MV
      integer MVAL
      logical NORESC
      integer LIW, LW, NVECS, IERR2
      double precision NOW
      integer ERMOLA, LDERLA
      double precision TOLREC, GROW
      double precision RATNO2
      integer NEINE2, FILT, MXA
      double precision DUMMY
      double precision   DLAMCH, EA18ED
      external           DLAMCH, EA18ED
      double precision   ZA12AD
      external ZA12AD
      external DLASET, DCOPY
      external EA16ED, EA16FD, EA16GD, EA16HD, EA16JD, EA16PD, EA16QD
      external EA16RD, EA16SD, EA16TD, EA16WD, EA16XD, EA16YD, EA16ZD
      external EA17AD, EA17BD, EA17CD, EA17ED, EA17FD, EA17LD, EA17MD
      external EA17ND, EA17PD, EA17QD
      external EA18GD, EA18CD, EA18HD, EA17DD, EA18XD, EA18YD, EA18WD
      intrinsic abs, sqrt
      if (IDO.eq.0) then
         do 10 I = 1,20
            INFO(I) = 0
 10      continue
         do 12 I = 1,10
            IPOS(I) = 0
 12      continue
         call EA16GD(N, NV, NWANT, BLK, LIW, LW, IWORK, LIWORK,
     &               ICNTL, INFO)
         if (info(IERR).lt.0) then
            IDO = 100
            return
         end if
         if (WHICH.eq.3 .or. WHICH.eq.-3. .or. WHICH.eq.1) then
            if (MODE.eq.2 .or. MODE.eq.4 .or. MODE.eq.5) then
               info(IERR) =  -6
            end if
            if (info(IERR).lt.0 .and. iwork(LUERR).ge.0)
     &         write(iwork(LUERR),9018)
         end if
 9018    format(
     &     'ERROR : WHICH and MODE are not compatible',/)
         if (LIWORK .lt. LIW) then
            if (iwork(LUERR).ge.0) then
               write(iwork(LUERR),9002) LIWORK, LIW
            end if
            info(IERR) = -12
            info(MEMO) = LIW
            IDO = 100
            return
         end if
 9002    format('ERROR : LIWORK is too small.  LIWORK = ', I10,'.',/,
     &          '        Expand to LIWORK = ', I10,/)
         if (LWORK .lt. LW) then
            if (iwork(LUERR).ge.0) then
               write(iwork(LUERR),9003) LWORK, LW
            end if
            info(IERR) = -13
            info(MEMO) = LW
            IDO = 100
            return
         end if
 9003    format('ERROR : LWORK is too small.  LWORK = ', I10,'.',/,
     &          '        Expand to LWORK = ', I10,/)
         if (MODE.lt.1 .or. MODE.gt.5) then
            info(IERR) = -1
            if (iwork(LUERR).ge.0) then
               write(iwork(LUERR),9005) MODE
            end if
            IDO = 100
            return
         end if
 9005    format('ERROR : MODE is out of range  (MODE=', I5, ').',/)
         if (WHICH.ne.1 .and. WHICH.ne.-1. and. WHICH.ne.2 .and.
     &       WHICH.ne.-2 .and. WHICH.ne.3 .and. WHICH.ne.-3 .and.
     &       WHICH.ne.4 .and. WHICH.ne.-4 .and. WHICH.ne.5 .and.
     &       WHICH.ne.10) then
            info(IERR) = -7
            if (iwork(LUERR).ge.0) then
               write(iwork(LUERR),9007) WHICH
            end if
            IDO = 100
            return
         end if
 9007    format('ERROR : WHICH is out of range  (WHICH=', I5, ').',/)
         if (LDV.lt.N) then
            info(IERR) = -8
            if (iwork(LUERR).ge.0) then
               write(iwork(LUERR),9008) LDV
            end if
            IDO = 100
            return
         end if
 9008    format('ERROR : LDV is too small  (LDV=', I5, ').',/)
         I = N
         if (MODE.eq.1 .or. MODE.eq.2) I = 1
         if (LDBV.lt.I) then
            info(IERR) = -9
            if (iwork(LUERR).ge.0) then
               write(iwork(LUERR),9009) LDBV, I
            end if
         end if
 9009    format('ERROR : LDBV is too small  (LDBV = ', I5, ').',/,
     &          '        Use LDBV >= ', I5,/)
         if (WHICH.eq.5) then
            if (RANGE(1).ge.RANGE(2)) then
               info(IERR) = -10
               if (iwork(LUERR).ge.0) then
                  write(iwork(LUERR),9010) RANGE(1), RANGE(2)
               end if
            end if
            if (MODE.eq.5 .and. icntl(IPREFI).gt.0
     &         .and. RANGE(1).le.ZERO .and. RANGE(2).ge.ZERO) then
               info(IERR) = -10
               if (iwork(LUERR).ge.0) then
                  write(iwork(LUERR),9110)
               end if
            end if
         end if
 9010    format('ERROR : RANGE(1) >= RANGE(2) (RANGE=',
     &          1PE12.5,',',1PE12.5, ').',/)
 9110    format('ERROR : RANGE(1) and RANGE(2) have different signs.')
         if (info(IERR).lt.0) then
            IDO = 100
            return
         end if
         work(RKSGRO) = abs(cntl(PLGROW))
         if (work(RKSGRO).lt.ONE) then
            call EA16ZD(info(IERR), 2**3)
            work(RKSGRO) = 500.0D0
            if (iwork(LUWARN).ge.0) then
               write(iwork(LUWARN),8003) work(RKSGRO)
            end if
         end if
 8003    format('WARNING : CNTL(4) is smaller than or equal to one.',/,
     &          '          The value used is ', 1PE10.2)
         work(RKSGRO) = sqrt(work(RKSGRO))
         LUWAR2 = iwork(LUWARN)
         call EA17ED(WHICH, MODE, icntl(IRSTRT),  iwork(IWHRES), LUWAR2)
         if (icntl(IRSTRT).gt.0 .and.
     &       icntl(IRSTRT).ne.iwork(IWHRES)) then
            call EA16ZD(info(IERR), 2**5)
         end if
         iwork(NWPOLE) = abs(ICNTL(CNWPOL))
         if (iwork(NWPOLE).gt.3) iwork(NWPOLE) = 0
         if (MODE.eq.1 .or. MODE.eq.3) then
            iwork(NWPOLE) = 0
         end if
         EPS = DLAMCH('Epsilon')
         SAFMIN = DLAMCH('SAFE MINIMUM')
         iwork(SEED)   = icntl(ISEED)
         iwork(SEED+1) = icntl(ISEED+1)
         iwork(SEED+2) = icntl(ISEED+2)
         iwork(SEED+3) = icntl(ISEED+3)
         do 20 I = SEED,SEED+3
            if (IWORK(I).lt.0 .or. IWORK(I).gt.4095) then
               if (iwork(LUWARN).ge.0)
     &            write(iwork(LUWARN),8004) I-SEED+ISEED, IWORK(I)
               call EA16ZD(info(IERR), 2**2)
               if (IWORK(I).lt.0) then
                  IWORK(I) = 0
               else if (IWORK(I).gt.4095) then
                  IWORK(I) = 4095
               end if
            end if
 20      continue
 8004    format(
     &     'WARNING : seed for random generation is out of bounds',/,
     &     '          ICNTL(',I2,') = ', I10,/)
         if (mod(iwork(SEED+3),2).eq.0) then
            if (iwork(LUWARN).ge.0) write(iwork(LUWARN),8044)
            call EA16ZD(info(IERR), 2**2)
            iwork(SEED+3) = iwork(SEED+3) + 1
         end if
 8044    format('WARNING : ICNTL(13) must be odd.',/)
         LU = iwork(LUERR)
         if (icntl(DIAGNO).gt.2 .and. LU.ge.0) then
            write(LU,6100)
            write(LU,6000)
            write(LU,6001) MODE
            write(LU,6002) IDO
            write(LU,6003) LDV
            write(LU,6004) LDBV
            write(LU,6005) N
            write(LU,6006) NV
            write(LU,6007) NWANT
            write(LU,6008) BLK
            write(LU,6010) WHICH
            if (WHICH.eq.1.or.WHICH.eq.-1.or.WHICH.eq.4.or.WHICH.eq.-4)
     &         write(LU,6111) RANGE(1)
            if (WHICH.eq.5)
     &         write(LU,6011) RANGE(1), RANGE(2)
            if (MODE.ne.1 .and. MODE.ne.3 .and.
     &          (WHICH.eq.2 .or. WHICH.eq.-2 .or. WHICH.eq.10) )
     &           write(LU,6012) SIGMA
            write(LU,6013) LIWORK
            write(LU,6014) LWORK
            write(LU,6015) (ICNTL(I),I = 1,18)
            write(LU,6016) (CNTL(I),I = 1,5)
            write(LU,6100)
         end if
 6100    format('----------------------------------------------------')
 6000    format('Arguments of EA16BD on input :')
 6001    format('  MODE     : ', I10)
 6002    format('  IDO      : ', I10)
 6003    format('  LDV      : ', I10)
 6004    format('  LDBV     : ', I10)
 6005    format('  N        : ', I10)
 6006    format('  NV       : ', I10)
 6007    format('  NWANT    : ', I10)
 6008    format('  BLK      : ', I10)
 6010    format('  WHICH    : ', I10)
 6011    format('  RANGE    : ', 1PE13.6,1X,1PE13.6)
 6111    format('  RANGE(1) : ', 1PE13.6)
 6012    format('  SIGMA    : ', 1PE13.6)
 6013    format('  LIWORK   : ', I10)
 6014    format('  LWORK    : ', I10)
 6015    format('  ICNTL    :', 5(I11),/,
     &         ('            ', 5(I11)))
 6016    format('  CNTL     :', 5(1PE11.3))
         BGS = MODE.gt.2
         if (BGS) then
            iwork(IBGS) = 1
         else
            iwork(IBGS) = 0
         end if
         EPSS = EPS*sqrt(real(N))*BLK
         work(PRANGE) = ZERO
         work(PRANGE+1) = ZERO
         if (WHICH.eq.-1 .or. WHICH.eq.1 .or. WHICH.eq.5
     &       .or. WHICH.eq.4 .or. WHICH.eq.-4) then
            work(PRANGE) = RANGE(1)
            if (WHICH.eq.5) then
               work(PRANGE+1) = RANGE(2)
            end if
         end if
         iwork(IFILT) = 0
         if (MODE.ge.4 .and. icntl(IPREFI).gt.0) then
            iwork(IFILT) = icntl(IPREFI)
         end if
         iwork(IMODE) = MODE
         info(NMVA) = 0
         info(NMVB) = 0
         iwork(NBACK) = -1
         info(NRESTA) = 0
         ISTEP = 1
         NSTEP = min(iwork(MSTEP),iwork(IMV)/BLK-1)
         NCONV = 0
         iwork(NWACO) = 0
         CVGED = .false.
         iwork(ICVG) = 0
         OFFSET = 0
         ERRMAX = ZERO
         work(GAMMA) = ZERO
         call EA16XD(NSTEP, 0, work(iwork(IRMOLA)), iwork(ILERLA),
     &               ZERO)
         info(NPOLE) = 0
         iwork(JUMP) = 0
         iwork(IVINIT) = icntl(VINIT)
         iwork(NLEJA) = 0
         iwork(INDONE) = 0
         work(PINTER) = ZERO
         work(PINTER+1) = ZERO
         work(PINTER+2) = ZERO
         work(PINTER+3) = ZERO
         iwork(JUMP)  = 0
         iwork(IEIGEN) = 0
         iwork(INCOOL) = 0
         iwork(ILSTEP) = 0
         iwork(IISTEP) = ISTEP
         iwork(INSTEP) = NSTEP
         iwork(IKKRYL) = 0
         iwork(INLANC) = 0
         iwork(NCOKEE) = 0
         iwork(IOFFST) = OFFSET
         iwork(INEIGV) = 0
         iwork(INPREF) = 0
         info(INCONV)  = NCONV
         work(PEPS)    = EPS
         work(PEPSS)   = EPSS
         work(PSAFM)   = SAFMIN
         work(PERRM)   = ERRMAX
      else if (iwork(JUMP).eq.0) then
         go to 1
      else if (iwork(JUMP).eq.4) then
         go to 4
      else if (iwork(JUMP).eq.5) then
         go to 5
      else
         go to 2
      end if
 5    continue
      if (iwork(IMODE).ne.1 .and. iwork(IMODE).ne.3) then
         iwork(JUMP) = 5
         if (abs(IDO).eq.4) then
            work(TIMFA2) = ZA12AD(DUMMY)
            if (iwork(LUREST).ge.0)
     &         write(iwork(LUREST),4805) work(TIMFA2)-work(TIMFA1)
         end if
         FIRST = .true.
         ISINF = .false.
         NVECS = 0
         MODE2 = iwork(IMODE)
         LUWAR2 = iwork(LUWARN)
         LURES2 = iwork(LUREST)
         if (IDO.eq.4) then
            NEINE2 = NEINEG
         else
            NEINE2 = -1
         end if
         FILT = iwork(IFILT)
         SAFMIN = work(PSAFM)
         call EA17ND(IDO, FIRST, icntl(CNWPOL).LT.0, ISINF,
     &               MODE2, 0, FILT, WHICH, SIGMA,
     &               work(PRANGE), work(RKSGRO), GROW,
     &               CVGED, work(iwork(IRITZ)), work(iwork(IRESID)),
     &               work(iwork(IRITZ)), 0, 0, 0, NWANT,
     &               NVECS, info(IERR), NEINE2, work(TRUST),
     &               iwork(iwork(IRANDM)), iwork(IEA17N),
     &               work(PEA17N), EXPRST, SAFMIN, EPS, N,
     &               LUWAR2, LURES2)
         if (IDO.eq.4) then
            info(NPOLE) = info(NPOLE) + 1
            work(TIMFA1) = ZA12AD(DUMMY)
            INFO(11) = 0
            IPOS(5) = iwork(IRITZ)
            IPOS(6) = IPOS(5)-1
            IPOS(7) = iwork(IRESID)
            IPOS(8) = IPOS(7)-1
            return
         end if
         if (info(IERR).eq.-16) then
            IDO = 100
            if (iwork(LUERR).ge.0) write(iwork(LUERR),9016)
            return
         else if (info(IERR).eq.-15) then
            IDO = 100
            if (iwork(LUERR).ge.0) write(iwork(LUERR),9017)
            return
         end if
 9016    format('ERROR : the user returned IDO=-4.',/,
     &          '        EA16 is not able to suggest a new pole.',/)
 9017    format('ERROR : the initial pole is too close to zero',/)
      else
         go to 1
      end if
      IDO = 0
      iwork(IISTEP) = 0
 4    continue
         MODE2 = iwork(IMODE)
         LURES2 = iwork(LUREST)
         VINIT2 = iwork(IVINIT)
         IISTE2 = iwork(IISTEP)
         NMVA2 = info(NMVA)
         NMVB2 = info(NMVB)
         FILT = iwork(IFILT)
         call EA16JD(IDO, MODE2, IPOS, N, BLK, VINIT2,
     &               FILT, iwork(SEED), V,LDV, IISTE2,
     &               NMVA2, NMVB2, LURES2)
         info(NMVA) = NMVA2
         info(NMVB) = NMVB2
         iwork(IISTEP) = IISTE2
         iwork(IVINIT) = VINIT2
         iwork(JUMP) = 4
         if (IDO.ne.100) return
      continue
      IDO = 0
      iwork(IISTEP) = 1
 1    continue
            BGS = iwork(IBGS).ne.0
            TOLREC = max(abs(cntl(TOLLAN)), abs(cntl(TOLLOC)))
            IERR2 = 0
            NSTEP = iwork(INSTEP)
            OFFSET = iwork(IOFFST)
            RANDOM = iwork(IRANDM)
            NORM = iwork(INORM)
            LDT = iwork(ILDT)
            IISTE2 = iwork(IISTEP)
            VINIT2 = iwork(IVINIT)
            LDCOEF = iwork(ILDCOE)
            LTEMP1 = iwork(ILTMP1)
            NMVA2 = info(NMVA)
            NMVB2 = info(NMVB)
            LUWAR2 = iwork(LUWARN)
            LUERR2 = iwork(LUERR)
            NCONV = info(INCONV)
            LDERV2 = iwork(LDERVE)
            LDERLA = iwork(ILERLA)
            LSTEP = iwork(ILSTEP)
            SAFMIN = work(PSAFM)
            EPSS = work(PEPSS)
            ERRMAX = work(PERRM)
            RATNO2 = ZERO
            if (iwork(IFILT).gt.0) RATNO2 = abs(cntl(RATNOR))
            NNSTEP = NSTEP
            NNSTEP = MAX(1,NSTEP)
            call EA16ED(IDO, IPOS, V, LDV, NV, BV, LDBV, IISTE2, LSTEP,
     &               NSTEP, NNSTEP, N, BLK, BGS, NCONV, OFFSET,
     &               VINIT2, WORK(iwork(ICOEF)), LDCOEF,
     &               WORK(iwork(IT)), LDT, IWORK(RANDOM),
     &               WORK(NORM), IWORK(iwork(GSSTRT)),
     &               LIWORK-iwork(GSSTRT)+1,
     &               IWORK(iwork(GSNUMB)), iwork(IMV)-BLK+1,
     &               WORK(iwork(ITEMP1)), LTEMP1, iwork(SEED),
     &               SAFMIN, EPSS, TOLREC, TOLREC, RATNO2,
     &               WORK(iwork(ERMOVE)), LDERV2,
     &               WORK(iwork(IRMOLA)), LDERLA,
     &               LUERR2, IERR2, NMVA2, NMVB2, iwork(IEA16E),
     &               work(PEA16E), ERRMAX)
            iwork(IISTEP) = IISTE2
            info(NMVA) = NMVA2
            info(NMVB) = NMVB2
            work(PERRM) = ERRMAX
            iwork(INSTEP) = NSTEP
            iwork(JUMP) = 0
            if (IDO.ne.100) return
 2       continue
         EPS    = work(PEPS)
         SAFMIN = work(PSAFM)
         ERRMAX = work(PERRM)
         EPSS   = work(PEPSS)
         F      = iwork(IF)
         LDF    = iwork(ILDF)
         LCF    = iwork(ILCF)
         LCF2   = iwork(ILCF2)
         TEMP3  = iwork(ITEMP3)
         LTEMP3 = iwork(ILTMP3)
         L      = iwork(IL)
         LDL    = iwork(ILDL)
         LCL    = iwork(ILCL)
         TAU    = iwork(ITAU)
         LTAU   = iwork(ILTAU)
         RANDOM = iwork(IRANDM)
         T      = iwork(IT)
         LDT    = iwork(ILDT)
         NORM   = iwork(INORM)
         TEMP1  = iwork(ITEMP1)
         LTEMP1 = iwork(ILTMP1)
         COEF   = iwork(ICOEF)
         LDCOEF = iwork(ILDCOE)
         PIPE   = iwork(IPIPE)
         TEMP2  = iwork(ITEMP2)
         EIGEN  = iwork(IEIGEN)
         RESID  = iwork(IRESID)
         RITZ   = iwork(IRITZ)
         EIGVEC = iwork(IEIGVE)
         LDEIGV = iwork(ILDEIG)
         LCEIGV = iwork(ILCEIG)
         LSTEP  = iwork(ILSTEP)
         ERMOLA = iwork(IRMOLA)
         LDERLA = iwork(ILERLA)
         ISTEP  = iwork(IISTEP)
         MLANC  = iwork(IMLANC)
         NKEEP  = iwork(INKEEP)
         NPURG  = iwork(INPURG)
         NSTEP  = iwork(INSTEP)
         KKRYL  = iwork(IKKRYL)
         NLANC  = iwork(INLANC)
         LLANC  = iwork(ILLANC)
         OFFSET = iwork(IOFFST)
         NEIGV  = iwork(INEIGV)
         MV     = iwork(IMV)
         CVGED  = iwork(ICVG).eq.1
         LLEJA  = iwork(ILEJA)
         MLEJA  = iwork(LMLEJA)
         NCONV  = info(INCONV)
         NCOOLD  = iwork(INCOOL)
         if (iwork(JUMP).eq.4000) then
            go to 4000
         else if (iwork(JUMP).eq.4100) then
            go to 4100
         else if (iwork(JUMP).eq.4200) then
            go to 4200
         else if (iwork(JUMP).eq.4300) then
            go to 4300
         else if (iwork(JUMP).eq.4800) then
            go to 4800
         else if (iwork(JUMP).eq.8700) then
            go to 8700
         else if (iwork(JUMP).gt.8000) then
            go to 8002
         else if (iwork(JUMP).eq.0) then
         end if
         iwork(GSFAIL) = 0
         if (IERR2.eq.-4) iwork(GSFAIL) = 1
         NLANC = ISTEP*BLK + OFFSET - NCONV
         if (NLANC.lt.BLK) then
            NLANC = BLK
            NCONV = OFFSET
            info(INCONV) = NCONV
         end if
         NVECS = NLANC + NCONV
         iwork(INLANC) = NLANC
         if (IERR2.ne.0 .and. IERR2.ne.-4) then
            info(IERR) = IERR2
            IDO = 100
            go to 7995
         end if
         NCOOLD = NCONV
         iwork(INCOOL) = NCOOLD
         NEIGV = NLANC + BLK
         iwork(INEIGV) = NEIGV
         if (iwork(IFILT).gt.0) then
            SHIFT = PIPE
            NSHIFT = max(0, min(iwork(IFILT), ISTEP-2))
            call DLASET('All', NSHIFT, 1, ZERO, ZERO,
     &                  work(SHIFT), NSHIFT)
            if (NSHIFT.gt.0) then
               if (iwork(LUREST).ge.0)
     &            write(iwork(LUREST),8001) NSHIFT
 8001             format('Performing ', I2, ' implicit filter steps.',/)
               call DLASET('All', NEIGV, NEIGV, ZERO, ONE,
     &                     WORK(EIGVEC), LDEIGV)
               iwork(INPREF) = NLANC
               KEEP = .TRUE.
               call EA17LD(work(T+LDT*(NCONV-OFFSET)), LDT, NLANC,
     &                  BLK, WORK(SHIFT), NSHIFT, WORK(EIGVEC),
     &                  LDEIGV, 0, LCEIGV, WORK(L), LDL, LCL, KEEP,
     &                  WORK(F), BLK, (L-F)/BLK, WORK(TAU), LTAU,
     &                  WORK(TEMP3), LTEMP3)
            end if
            NLANC = NLANC - NSHIFT*BLK
            ISTEP = ISTEP - NSHIFT
            iwork(INLANC) = NLANC
            iwork(IISTEP) = ISTEP
            iwork(INSHIF) = NSHIFT
         end if
         NVECS = NLANC + NCOOLD
         NORMT = EA18ED(NCOOLD, NVECS-NCOOLD, BLK, work(RITZ),
     &                  work(T+(NCOOLD-OFFSET)*LDT), LDT)
         work(PNORMT) = NORMT
         call DLASET('All', NLANC+BLK, NLANC+BLK, ZERO, ONE,
     &               work(EIGVEC), LDEIGV)
         call EA16FD(work(T+LDT*(NCOOLD-OFFSET)), LDT, NLANC, BLK,
     &               WORK(RITZ+NCOOLD), WORK(EIGVEC), LDEIGV,
     &               WORK(RESID+NCOOLD), WORK(PIPE), BLK,
     &               WORK(TEMP3), LTEMP3, IERR2)
         iwork(IHARM) = 0
         if (IERR2.ne.0) then
            iwork(IVINIT) = 0
            LSTEP = 0
            IDO = 0
            go to 11
         end if
         call EA16PD(NCOOLD, NCONV, NVECS, WORK(RITZ),
     &               WORK(RESID), WORK(EIGVEC), LDEIGV, NLANC+BLK,
     &               work(PIPE), BLK, CNTL(TOL), cntl(TOLLAN),
     &               IWORK(RANDOM), WORK(TEMP3), EPSS)
         info(INCONV) = NCONV
         if (iwork(IFILT).gt.0) then
            call EA16QD(NCONV-NCOOLD, BLK, WORK(RITZ+NCOOLD),
     &                  WORK(EIGVEC), LDEIGV, NLANC+BLK, WORK(PIPE),
     &                  BLK, WORK(EIGVEC+LDEIGV*NLANC), LDEIGV, SAFMIN)
         end if
         if (iwork(IFILT).gt.0) then
            call EA18CD(WORK(F), BLK, WORK(TAU), BLK, iwork(INPREF),
     &                  iwork(INSHIF), WORK(EIGVEC), LDEIGV, LCEIGV,
     &                  NLANC+BLK, WORK(TEMP3), LTEMP3)
         end if
         HARMON = iwork(IHARM)
         NVECS = NCOOLD + NLANC
         call EA18YD(MODE, WHICH, RANGE, HARMON,
     &              work(T), LDT, MLANC+BLK, OFFSET, NLANC,
     &              NCOOLD, NCONV, BLK, SIGMA, work(RITZ), NVECS,
     &              work(PIPE), work(RESID), work(EIGVEC), LDEIGV,
     &              LCEIGV, NEIGV,
     &              work(L), LDL, LCL, work(F), LDF, LCF, work(TAU),
     &              LTAU, work(TEMP3), LTEMP3, SAFMIN, IERR2)
         if (IERR2.ne.0) then
            call EA16SD(V(1,NCOOLD+1), LDV, MV-NCOOLD, N, NEIGV,
     &                  WORK(EIGVEC), LDEIGV, LCEIGV,
     &                  NCONV-NCOOLD, WORK(F), LDF, LCF2)
            IDO = 0
            LSTEP = 0
            OFFSET = NCONV
            NSTEP = min(NV/BLK-1, (NVECS-OFFSET)/BLK)
            iwork(IOFFST) = OFFSET
            iwork(INSTEP) = NSTEP
            iwork(IVINIT) = 0
            go to 11
         end if
         iwork(IHARM) = HARMON
         IDO = 0
         EIGEN = TEMP3
         iwork(IEIGEN) = EIGEN
         iwork(JUMP) = 4100
 4100    continue
         NVECS  = NLANC+NCOOLD
         LURES2 = iwork(LUREST)
         MODE2  = iwork(IMODE)
         HARMON = iwork(IHARM)
         IDUM1 = iwork(NWACO)
         NCOLEF = 0
         NCORIG = 0
         call EA17AD(IDO, WHICH, MODE2, HARMON, SIGMA, NCOLEF,
     &               NCORIG, IDUM1, NCONV, NVECS, WORK(RITZ),
     &               work(PRANGE), IWORK(RANDOM), WORK(EIGEN),
     &               SAFMIN)
         iwork(NWACO) = IDUM1
         iwork(INCOLE) = NCOLEF
         iwork(INCORI) = NCORIG
         if (IDO.ne.100) then
            IPOS(5) = EIGEN
            IPOS(6) = EIGEN+NVECS-1
            IPOS(9) = RANDOM
            IPOS(10) = RANDOM+NVECS-1
            return
         end if
         NCOLEF = iwork(INCOLE)
         NCORIG = iwork(INCORI)
         HARMON = iwork(IHARM)
         call EA17MD(CVGED, INFO8, WHICH, iwork(IMODE), HARMON, SIGMA,
     &              work(PRANGE), NWANT, iwork(NWACO), NCOLEF, NCORIG,
     &              iwork(RANDOM), NCOOLD+NLANC, WORK(RITZ),
     &              WORK(TEMP3), SAFMIN, iwork(LUWARN), info(IERR),
     &              MVAL)
         INFO(8) = INFO8
         CVGED = CVGED .and. iwork(GSFAIL).eq.0
         if (CVGED) iwork(NWACO) = min(NWANT,iwork(NWACO))
         if (CVGED) then
            iwork(ICVG) = 1
         else
            iwork(ICVG) = 0
         end if
         iwork(NCOKEE) = min(NWANT+LLANC, NCOOLD+LLANC, MV-3*BLK)
         iwork(NCOKEE) = min(NWANT+LLANC, MV-3*BLK)
         iwork(NCOKEE) = max(iwork(NCOKEE),min(NWANT,iwork(NWACO)))
         NCOKE2 = iwork(NCOKEE)
         NVECS = NLANC+NCOOLD
         call EA17BD(NCONV, NCOOLD, NCOKE2, NVECS, iwork(RANDOM))
         info(INCONV) = NCONV
         INFO(9) = min(NCONV,MVAL)
         iwork(NCOKEE) = NCOKE2
         call EA17CD(3, NCOOLD, NVECS, work(RITZ), work(RESID+NCOOLD),
     &               work(EIGVEC), LDEIGV, NEIGV, work(PIPE), BLK,
     &               iwork(RANDOM), work(TEMP3))
         if (iwork(LUREST).ge.0) then
            write(iwork(LUREST),9900) info(NRESTA)
           write(iwork(LUREST),9901) NCONV
            write(iwork(LUREST),9902) INFO(8)
            write(iwork(LUREST),9904) INFO(9)
            if (NCONV.lt.NCOOLD+NLANC) then
                write(iwork(LUREST),9903) NCONV+1, WORK(RESID+NCONV)
            end if
         end if
 9900    format(/'Restart ', I6, ' :')
 9901    format('  Number of converged Ritz values                  : ',
     +   I8)
 9902    format('  Number of wanted converged Ritz values           : ',
     +   I8)
 9903    format('  Residual norm of the ', I4, 'th Ritz value : ',
     &          1PE10.2)
 9904    format('  Number of converged Ritz values satisfying WHICH : ',
     +   I8)
c
         GROW = work(RKSGRO)
         NORMT = work(PNORMT)
         LARGE = DLAMCH('Overflow')
         call EA18GD(work(GAMMA), .false., EPSS, LARGE, NORMT, GROW)
         if (NCONV.gt.iwork(NCOKEE)) then
            LURES2 = iwork(LUREST)
            NCONV2 = NCONV
            NCOOL2 = NCOOLD
            call EA16HD(NCONV2, NCOOL2, iwork(NCOKEE), NLANC,
     &                  iwork(RANDOM), work(RITZ), work(RESID),
     &                  V, LDV, N, NCOOLD+NEIGV, LURES2, NCONV, NCOOLD)
            iwork(INCOOL) = NCOOLD
            info(INCONV) = NCONV
            iwork(NWACO) = min(iwork(NWACO), NCONV)
         end if
         NKEEP = min(LLANC,NCOOLD+NLANC,MV-2*BLK)
         LSTEP = min(iwork(LLSTEP), max(0,(NKEEP-NCONV)/BLK))
         LSTEP = max(0,(NKEEP-NCONV)/BLK)
         if (iwork(GSFAIL).ne.0) LSTEP = 0
         NKEEP = NCONV+LSTEP*BLK
         call EA18HD(NLANC+NCOOLD-NKEEP, NPURG, WORK(RITZ+NKEEP),
     &               WORK(RESID+NKEEP),
     &               WORK(EIGVEC+LDEIGV*(NKEEP-NCOOLD)), LDEIGV, NEIGV,
     &               WORK(PIPE+BLK*(NKEEP-NCOOLD)), BLK, BLK,
     &               sqrt(EPS)*NORMT, WORK(TEMP3), IWORK(RANDOM))
         NPURG = NPURG + NKEEP
         NCONV = min(NCONV,NPURG)
         OFFSET = mod(NPURG-NCONV,BLK)
         OFFSET = NCONV - OFFSET
         if (OFFSET.lt.0) then
            NPURG = min(NLANC, NPURG - OFFSET)
            OFFSET = 0
         end if
         NKEEP = min(NKEEP, NPURG)
         LSTEP = min(iwork(LLSTEP), max(0,(NKEEP-OFFSET)/BLK))
         LSTEP = max(0,(NKEEP-NCONV)/BLK)
         if (iwork(GSFAIL).ne.0) LSTEP = 0
         if (LSTEP.eq.0) then
            OFFSET = NCONV
         end if
         NKEEP = OFFSET + LSTEP*BLK
         NCONV = min(NCONV,NKEEP)
         KKRYL = OFFSET + iwork(MSTEP)*BLK
         DIFF = max(0,KKRYL+BLK-MV)
         NSTEP = iwork(MSTEP) - (DIFF+BLK-1)/BLK
         KKRYL = OFFSET + NSTEP*BLK
         if (iwork(IWHRES).eq.1) then
            NPURG = NKEEP
            NCONV = min(NCONV, NPURG)
         end if
         iwork(ILSTEP) = LSTEP
         iwork(INKEEP) = NKEEP
         iwork(IKKRYL) = KKRYL
         iwork(INSTEP) = NSTEP
         iwork(IOFFST) = OFFSET
         iwork(INPURG) = NPURG
         info(INCONV) = NCONV
         if (LSTEP.eq.0) NPURG = OFFSET
         if (LSTEP.gt.0) then
            call DLASET('All', BLK, NCONV-OFFSET, ZERO, ZERO,
     &               work(T+1),LDT)
            call DCOPY(NCONV-OFFSET, work(RITZ+OFFSET),1, work(T),LDT)
            call EA16RD(NCOOLD+NLANC-NCONV, BLK, 0,
     &               NPURG-NCONV, work(T+LDT*(NCONV-OFFSET)),
     &               LDT, WORK(RITZ+NCONV),
     &               WORK(EIGVEC+LDEIGV*(NCONV-NCOOLD)), LDEIGV,
     &               LCEIGV, NEIGV, WORK(PIPE+BLK*(NCONV-NCOOLD)),
     &               BLK, WORK(L), LDL, LCL)
         end if
         if (iwork(NWPOLE).gt.0 .and. CVGED .and. WHICH.ne.10) then
            if (NEINEG.ge.0) then
               NVECS = NLANC+NCOOLD
               call DCOPY(NVECS, WORK(RITZ),1, WORK(F),1)
               IDO = 0
               EIGEN = TEMP3
               iwork(IEIGEN) = EIGEN
               LURES2 = iwork(LUREST)
               MODE2 = iwork(IMODE)
               HARMON = iwork(IHARM)
               call EA17AD(IDO, WHICH, MODE2, HARMON, SIGMA,
     &                     IDUM1, IDUM2, IDUM3,
     &                     0, NCONV, WORK(F), work(PRANGE),
     &                     IWORK(RANDOM), WORK(EIGEN), SAFMIN)
               call EA17CD(1, 0, NCONV, work(F), work(RESID),
     &                    work(EIGVEC), LDEIGV, 1, work(PIPE), 1,
     &                    iwork(RANDOM), work(EIGEN))
               EIGEN = TEMP3
               NVECS = NLANC+NCOOLD
               IPOS(5) = EIGEN
               IPOS(6) = EIGEN+NVECS-1
               IPOS(7) = RESID
               IPOS(8) = RESID+NVECS-1
               call DCOPY(NVECS, WORK(F),1, WORK(EIGEN),1)
            else
               go to 4801
            end if
         else
            go to 4801
         end if
         iwork(IDOGS) = 0
         IDO = 0
         iwork(JUMP) = 4800
 4800    continue
            if (abs(IDO).eq.4) then
               work(TIMFA2) = ZA12AD(DUMMY)
               if (iwork(LUREST).ge.0)
     &            write(iwork(LUREST),4805) work(TIMFA2)-work(TIMFA1)
            end if
 4805       format('  Factorization time : ', 1PE9.2, ' seconds.')
            MODE2 = iwork(IMODE)
            LURES2 = iwork(LUREST)
            FILT = iwork(IFILT)
            HARMON = iwork(IHARM)
            NVECS = NLANC+NCOOLD
            if (IDO.eq.-4) iwork(IDOGS) = -4
            IDOGS2 = iwork(IDOGS)
            GROW = work(RKSGRO)
            NORESC = ICNTL(PRIRES).eq.0
            call EA17DD(IDOGS2, MODE2, HARMON, FILT, N, WHICH,
     &                  SIGMA, NEINEG, work(EIGEN), NVECS, NCONV,
     &                  NWANT, work(TRUST), RANGE, work(PEA17N),
     &                  iwork(IEA17N), CVGED, GROW, EPS, SAFMIN,
     &                  NORESC, LURES2)
            iwork(IDOGS) = IDOGS2
            if (CVGED) then
               iwork(ICVG) = 1
            else
               iwork(ICVG) = 0
            end if
            if (iwork(IDOGS).ge.4 .and. iwork(IDOGS).le.6) then
               work(TIMFA1) = ZA12AD(DUMMY)
               info(NPOLE) = info(NPOLE) + 1
               IDO = 4
               return
            end if
 4801    continue
         if (CVGED) then
            NKEEP = NPURG
            iwork(INKEEP) = NKEEP
            go to 4600
         end if
         iwork(IDOGS) = 0
         if (iwork(IWHRES).gt.1) then
            SHIFT = F
            if (iwork(IDOGS).eq.0) then
               NSHIFT = LDF
               NDONE = iwork(INDONE)
               IDOGS2 = iwork(IDOGS)
               NLEJA2 = iwork(NLEJA)
               IDUM1 = iwork(IWHRES)
               call EA17FD(IDOGS2, info(NRESTA), NPURG-NCONV,
     &                     NKEEP-NCONV, BLK, WHICH,
     &                     work(PRANGE), work(PINTER),
     &                     WORK(SHIFT), LDF, WORK(RITZ+NCONV),
     &                     IDUM1, NSHIFT, MLEJA, NLEJA2, NDONE,
     &                     IWORK(LLEJA), WORK(iwork(RLEJA)),
     &                     SAFMIN)
               iwork(IDOGS) = IDOGS2
               iwork(NLEJA) = NLEJA2
               iwork(INDONE) = NDONE
               if (IDOGS2.ne.100) then
                  IDO = 7
                  IPOS(5) = RITZ
                  IPOS(6) = RITZ+NPURG-1
                  IPOS(7) = SHIFT
                  IPOS(8) = SHIFT+NSHIFT-1
                  info(INCONV) = NCONV
                  iwork(JUMP) = 4200
                  iwork(IDOGS) = 1
                  return
               end if
            end if
         end if
 4200    continue
         if (iwork(IWHRES).gt.1) then
            if (IDO.eq.7) then
               NSHIFT = IPOS(8)-IPOS(7)+1
               IDO = 0
               SHIFT = F
            end if
            KEEP = .FALSE.
            call EA17LD(work(T+LDT*(NCONV-OFFSET)), LDT, NPURG-NCONV,
     &                  BLK, WORK(SHIFT), NSHIFT,
     &                  WORK(EIGVEC+(NCONV-NCOOLD)*LDEIGV), LDEIGV,
     &                  NEIGV, LCEIGV-NCONV+NCOOLD, WORK(L),
     &                  LDL, LCL, KEEP, WORK(PIPE), BLK, MLANC,
     &                  WORK(TAU), LTAU, WORK(TEMP3), LTEMP3)
            iwork(IDOGS) = 0
         end if
         SETPOL = .false.
         if (iwork(NWPOLE).gt.0) then
          if (abs(ICNTL(CNWPOL)).eq.2 .or. abs(ICNTL(CNWPOL)).eq.3)
     &      SETPOL = info(NRESTA)-iwork(NBACK).ge.abs(icntl(POLRTR))
          NOW = ZA12AD(DUMMY)
          if (iwork(LUREST).ge.0) then
            write(iwork(LUREST),4010) NOW-work(TIMFA2)
 4010       format('  Time since factorization : ', 1PE9.2)
          end if
          if (abs(ICNTL(CNWPOL)).eq.1 .or. abs(ICNTL(CNWPOL)).eq.3)
     &      SETPOL = SETPOL .or.
     &       (NOW-work(TIMFA2) .gt.
     &       abs(icntl(POLTIM)) * 0.01 * (work(TIMFA2)-work(TIMFA1)))
         end if
         if (icntl(CHAMOD).ne.0
     &      .and. (iwork(IMODE).eq.1 .or. iwork(IMODE).eq.3)
     &      .and. (abs(WHICH).eq.2 .or. WHICH.eq.10) ) then
            SETPOL = .true.
            iwork(NWPOLE) = abs(ICNTL(CNWPOL))
         end if
         if (SETPOL) then
            iwork(ISETPL) = 1
         else
            iwork(ISETPL) = 0
         end if
         if (iwork(ISETPL).ne.0) then
            EIGEN = TEMP3
            NVECS = NLANC+NCOOLD
         end if
         if (iwork(IMODE).eq.1 .or.  iwork(IMODE).eq.3) then
            IDO = 0
         else
            IDO = 2
         end if
         iwork(JUMP) = 4000
 4000    continue
         if (iwork(ISETPL).ne.0) then
            if (IDO.eq.5) then
                if (MODE.eq.iwork(IMODE)) then
                   go to 4600
                else if ( (iwork(IMODE).eq.1 .or. iwork(IMODE).eq.3)
     &                    .and. MODE.eq.iwork(IMODE)+1) then
                   IDO = 3
                else
                   if (iwork(LUERR).ge.0) write(iwork(LUERR),9015)
     &                iwork(IMODE), iwork(IMODE)+1
                   IDO = 100
                   info(IERR) = -1
                   return
                end if
            end if
 9015       format('ERROR : MODE is wrong. Choose MODE=',I1,' or ',I1,
     &             '.')
            FIRST = iwork(IMODE).eq.1 .or. iwork(IMODE).eq.3
            ISINF = FIRST
            if (abs(IDO).eq.4) then
               work(TIMFA2) = ZA12AD(DUMMY)
               if (iwork(LUREST).ge.0)
     &            write(iwork(LUREST),4805) work(TIMFA2)-work(TIMFA1)
            end if
 4700       continue
            MODE2 = iwork(IMODE)
            LUWAR2 = iwork(LUWARN)
            LURES2 = iwork(LUREST)
            NWACO2 = iwork(NWACO)
            if (MODE2.eq.2.or.MODE2.eq.4.or.MODE2.eq.5) then
               NEINE2 = NEINEG
            end if
            NVECS = NLANC+NCOOLD
            FILT = iwork(IFILT)
            HARMON = iwork(IHARM)
            EXPRST = .false.
            call EA17ND(IDO, FIRST, icntl(CNWPOL).LT.0, ISINF, MODE2,
     &               HARMON, FILT, WHICH, SIGMA, work(PRANGE),
     &               work(RKSGRO), GROW, CVGED, work(RITZ),
     &               work(RESID), work(EIGEN), NCONV, NWACO2, NVECS,
     &               NWANT, NKEEP, info(IERR), NEINEG, work(TRUST),
     &               iwork(RANDOM), iwork(IEA17N), work(PEA17N),
     &               EXPRST, SAFMIN, EPS, N, LUWAR2, LURES2)
            if (IDO.eq.4) work(GAMGRO) = GROW
            if (EXPRST) then
               LSTEP = 0
               OFFSET = NCONV
               NSTEP = min(NV/BLK-1, (NVECS-OFFSET)/BLK)
               iwork(ILSTEP) = 0
               iwork(IOFFST) = OFFSET
               iwork(INSTEP) = NSTEP
               iwork(IVINIT) = 0
            end if
            IPOS(5) = EIGEN
            IPOS(6) = EIGEN+NVECS-1
            IPOS(7) = RESID
            IPOS(8) = RESID+NVECS-1
            if (IDO.eq.4) then
               info(NPOLE) = info(NPOLE) + 1
               work(TIMFA1) = ZA12AD(DUMMY)
            end if
            if (IDO.eq.3 .and. FIRST .and.
     &           (MODE.eq.1 .or. MODE.eq.3) ) IDO = 5
            if (IDO.eq.3 .and. ICNTL(CNWPOL).ge.0) go to 4700
            if (IDO.ne.100 .and. IDO.lt.6) go to 8000
            if (IDO.eq.6) then
               NVECS = NLANC + NCOOLD
               OLDSIG = work(PEA17N+4)
               INFLAN = FIRST .or. iwork(IHARM).ne.0
               INFLOC = FIRST
               call EA17QD(OLDSIG, SIGMA, INFLAN, INFLOC, iwork(IMODE),
     &                  NVECS, NKEEP, BLK, NCONV, OFFSET,
     &                  WORK(RITZ), WORK(RESID), MV-BLK,
     &                  work(T), LDT, MV-OFFSET-BLK,
     &                  WORK(L), LDL, LCL, WORK(F), LDF, LCF,
     &                  WORK(EIGVEC+(NCONV-NCOOLD)*LDEIGV), LDEIGV,
     &                  LCEIGV-NCONV+NCOOLD, NEIGV, WORK(TAU), LTAU,
     &                  WORK(EIGEN), LTEMP3, SAFMIN)
               iwork(NBACK) = info(NRESTA)
               if (iwork(IMODE).ne.MODE) then
                  LUWAR2 = iwork(LUWARN)
                  call EA17ED(WHICH, MODE, icntl(IRSTRT), iwork(IWHRES),
     &                        LUWAR2)
                  if (icntl(IRSTRT).gt.0 .and.
     &                icntl(IRSTRT).ne.iwork(IWHRES)) then
                     call EA16ZD(info(IERR), 2**5)
                  end if
                  iwork(IMODE) = MODE
               end if
               iwork(IHARM) = 0
               GROW = work(GAMGRO)
               LARGE = DLAMCH('Overflow')
               call EA18GD(work(GAMMA), .true., EPSS, LARGE,
     &                     ZERO, GROW)
            else
               if (info(IERR).eq.-16) then
                  if (iwork(LUERR).ge.0) write(iwork(LUERR),9016)
                  IDO = 100
                  go to 7995
               end if
            end if
         end if
         CVGED = .false.
         iwork(ICVG) = 0
 4600    continue
         if (iwork(IHARM).ne.0) then
            call EA18XD(MODE, work(T),
     &              LDT, NKEEP-OFFSET, OFFSET, NKEEP, NCOOLD, NCONV,
     &              BLK, SIGMA, work(RITZ),
     &              work(EIGVEC), LDEIGV, LCEIGV, NEIGV,
     &              work(L), LDL, LCL, work(F), LDF, LCF, work(TAU),
     &              LTAU, work(TEMP3), LTEMP3, SAFMIN, IERR2)
            iwork(IHARM) = 0
            if (IERR2.ne.0) then
               call EA16SD(V(1,NCOOLD+1), LDV, MV-NCOOLD, N, NEIGV,
     &                     WORK(EIGVEC), LDEIGV, LCEIGV,
     &                     NCONV-NCOOLD, WORK(F), LDF, LCF2)
               IDO = 0
               LSTEP = 0
               OFFSET = NCONV
               NSTEP = min(NV/BLK-1, (NVECS-OFFSET)/BLK)
               iwork(IOFFST) = OFFSET
               iwork(INSTEP) = NSTEP
               iwork(IVINIT) = 0
               go to 11
            end if
         end if
         call EA16SD(V(1,NCOOLD+1), LDV, MV-NCOOLD, N, NEIGV,
     &               WORK(EIGVEC), LDEIGV, LCEIGV,
     &               NKEEP+BLK-NCOOLD, WORK(F), LDF, LCF2)
         if (CVGED .or. ISTEP.eq.0) then
            IDO = 100
            NVECS = NKEEP
            goto 8500
         end if
 4300    continue
         if (info(IERR).eq.-14) then
            info(IERR) = 0
            call DCOPY(NCONV, work(F),1, work(RITZ),1)
         end if
         MXA = icntl(MXNMVA)
         if (MXA.lt.1) MXA = 2
         if (info(NMVA)+(iwork(MSTEP)-LSTEP).ge.MXA*N) then
            IDO = 100
            if (iwork(LUERR).ge.0) write(iwork(LUERR),9014)
            NVECS = NKEEP
            info(IERR) = -14
            iwork(JUMP) = 4300
            go to 8500
         end if
 9014    format('ERROR : The maximum number of operations with OP',
     &          ' is reached',/)
         call EA16TD(work(T), LDT, BLK, LSTEP, NCONV,
     &               OFFSET, WORK(ERMOLA), LDERLA,
     &               WORK(COEF), LDCOEF, WORK(NORM), WORK(TEMP1),
     &               LTEMP1, IERR2)
 11      continue
         GAMMA2 = work(GAMMA)
         call EA16WD(work(RITZ), NCONV, WORK(RESID), ERRMAX,
     &               GAMMA2, WORK(iwork(ERMOVE)), iwork(LDERVE))
         if (ERRMAX.ge.ONE) ERRMAX = ONE
         call EA16YD(LSTEP, 3, WORK(ERMOLA), LDERLA, ERRMAX)
         GAMMA2 = work(GAMMA)
         call EA16XD(iwork(MSTEP), LSTEP, work(ERMOLA),
     &               LDERLA, GAMMA2)
         ISTEP = LSTEP + 1
         NVECS = KKRYL
         iwork(IISTEP) = ISTEP
         IDO = 0
         info(NRESTA) = info(NRESTA) + 1
         iwork(IVINIT) = 1
         if (iwork(GSFAIL).ne.0) then
            iwork(GSFAIL) = 0
            iwork(IVINIT) = 0
         end if
      go to 1
 8500 continue
 7995 continue
      iwork(IDOGS) = IDO
      IDO = 0
      EIGEN = TEMP2
      iwork(IEIGEN) = EIGEN
      iwork(JUMP) = 8700
 8700 continue
         LURES2 = iwork(LUREST)
         MODE2  = iwork(IMODE)
         HARMON = iwork(IHARM)
         if (ISTEP.ne.0) then
           call EA17AD(IDO, WHICH, MODE2, HARMON, SIGMA, IDUM1, IDUM2,
     &               IDUM3, 0, NCONV, WORK(RITZ), work(PRANGE),
     &               IWORK(RANDOM), WORK(EIGEN), SAFMIN)
            IDUM1 = iwork(NWACO)
         else
            NCOLEF = 0
            NCORIG = 0
            call EA17AD(IDO, WHICH, MODE2, HARMON, SIGMA, NCOLEF,
     &                  NCORIG, IDUM1, NCONV, NCONV, WORK(RITZ),
     &                  work(PRANGE), IWORK(RANDOM), WORK(EIGEN),
     &                  SAFMIN)
            iwork(NWACO) = IDUM1
         call EA17MD(CVGED, INFO8, WHICH, iwork(IMODE), HARMON, SIGMA,
     &              work(PRANGE), NWANT, iwork(NWACO), NCOLEF, NCORIG,
     &              iwork(RANDOM), NCOOLD+NLANC, WORK(RITZ),
     &              WORK(TEMP3), SAFMIN, iwork(LUWARN), info(IERR),
     &              MVAL)
            INFO(8) = INFO8
            INFO(9) = MIN (NCONV,MVAL)
         end if
      if (IDO.ne.100) then
         IPOS(5) = EIGEN
         IPOS(6) = EIGEN+NCONV-1
         IPOS(9) = RANDOM
         IPOS(10) = RANDOM+NCONV-1
         IDO = 6
         return
      end if
      call EA17CD(2, 0, NCONV, work(RITZ), work(RESID), V, LDV, N,
     &            work(PIPE),1, iwork(RANDOM), work(EIGEN))
      if (ICNTL(PRIRES).ne.0) IDO = 0
 8002 continue
      if (ICNTL(PRIRES).ne.0) then
         LARGE = ONE/SAFMIN
         IPOS7 = IPOS(7)
         call EA18WD(IDO, MODE, IPOS, NCONV, N, BLK, V, LDV, BV, LDBV,
     &               work(RITZ), work(RESID), NV-BLK+1, SAFMIN,
     &               LARGE, iwork(JUMP), IPOS7)
         IPOS(7) = IPOS7
         if (IDO.ne.100) return
      end if
      iwork(JUMP) = 4300
      IDO = iwork(IDOGS)
      if (info(IERR).eq.-14) call DCOPY(NCONV, work(RITZ),1, work(F),1)
      call EA17PD(MODE, iwork(IHARM), NCONV, WORK(RITZ),
     &            SIGMA, RANGE, SAFMIN)
 8000 continue
      if (IDO.eq.100) then
         IPOS(5) = RITZ
         IPOS(6) = RITZ+NCONV-1
         IPOS(7) = RESID
         IPOS(8) = RESID+NCONV-1
      end if
      LU = iwork(LUERR)
      if (info(IERR).ge.0 .and. INFO(8).lt.NWANT .and. IDO.eq.100) then
         if (LU.ge.0 .and. info(IERR).ne.2**4) write(LU,9099)
         call EA16ZD(info(IERR), 2**4)
      end if
 9099 format(/'WARNING : the number of eigenvalues computed',
     &   ' is smaller than NWANT.',/)
      if (icntl(DIAGNO).gt.2 .and. LU.ge.0 .and. IDO.eq.100) then
         write(LU,7100)
         write(LU,7000)
         write(LU,7001) (IPOS(I),I = 1,10)
         write(LU,7002) (INFO(I),I = 1,9)
         if (MODE.ne.1 .and. MODE.ne.3) write(LU,7003) SIGMA
         write(LU,7100)
      end if
 7100 format('--------------------------------------------------------')
 7000 format('Arguments of EA16BD on output :')
 7001 format('IPOS   :', 5(1X,I10),/,
     &       '        ', 5(1X,I10))
 7002 format('INFO   :', 5(1X,I10),/
     &       '        ', 4(1X,I10))
 7003 format('SIGMA  :', 1PE13.6)
      if (IDO.eq.100.and.ICNTL(2).ge.4.and.ICNTL(1).ge.0) then
         NWACO2 = iwork(NWACO)
         if (ICNTL(PRIRES).ne.0) then
            write(ICNTL(1),7053)
            do 8911 I = 1, NCONV
               write(ICNTL(1), 7054 ) I, WORK(RITZ-1+I), WORK(RESID-1+I)
 8911       continue
         else
            write(ICNTL(1),7050)
            do 8910 I = 1, NCONV
               write(ICNTL(1), 7051 ) I, WORK(RITZ-1+I)
 8910       continue
         end if
         write(ICNTL(1),7052)
      end if
 7050 format(
     + '+--------------- Ritz values on output ----------------+',/)
 7051 FORMAT( ' ', I6, ':  ', 1PE13.5 )
 7052 format('+------------------------------------------------------+')
 7053 format(
     + '+------- Ritz values / residual norms on output -------+',/)
 7054 FORMAT( ' ', I6, ':  ', 1PE13.5,' /  ',  1PE13.5 )
      return
      end
      subroutine EA16ED(IDO, IPOS, V, LDV, NV, BV, LDBV, ISTEP, LSTEP,
     &                  NSTEP, NNSTEP, N, BLK, BGS, NCONV, OFFSET,
     &                  VINIT, COEF, LDCOEF, T, LDT, RANDOM, NORM,
     &                  GSSTRT, LGSSTR, GSNUMB, LGSNUM, TEMP, LTEMP,
     &                  ISEED, SAFMIN, EPSS, TOLLAN, TOLLOC, RATNOR,
     &                  ERMOVE, LDERVE, ERMOLA, LDERLA, LUERR,
     &                  IERR, NMVA, NMVB, ISAVE, RSAVE, ERRMAX)
      logical BGS
      integer IDO, ISTEP, LSTEP, NSTEP, NNSTEP, N, NV, BLK, VINIT, LDV
      integer LDBV, LDCOEF, LDT
      integer NCONV, OFFSET, LDERVE, LDERLA, LTEMP, LUERR
      integer LGSSTR, LGSNUM
      integer IERR, NMVA, NMVB
      double precision SAFMIN, EPSS, TOLLAN, TOLLOC, RATNOR, ERRMAX
      integer ISEED(4), IPOS(4), RANDOM(BLK)
      integer GSSTRT(LGSSTR), GSNUMB(LGSNUM), ISAVE(16)
      double precision V(LDV,NV), BV(LDBV,BLK)
      double precision COEF(LDCOEF,BLK), T(LDT,NCONV-OFFSET+NNSTEP*BLK)
      double precision NORM(BLK), TEMP(LTEMP)
      double precision ERMOVE(LDERVE,4), ERMOLA(LDERLA,6)
      double precision RSAVE(4)
      integer GSMAX
      parameter (GSMAX=2)
      double precision DGKS, ZERO, ONE
      parameter       (DGKS=0.5D0, ZERO=0.0D0, ONE=1.0D0)
      integer ERRLAN, MINRAT, REFRAT
      parameter (ERRLAN=1, MINRAT=3, REFRAT=4)
      integer SECOND, THIRD
      parameter (SECOND=2, THIRD=3)
      integer I, K, MTASKS
      integer DIAGT, OFFT
      integer NTASK, IDOGS
      double precision ERRLOC
      integer PIDOGS, JUMP, PVECGS, PVECI, PVECIM, PVECIP, PNTASK
      integer HASBV, IGS, JSTEP
      parameter (PIDOGS=1, JUMP=2, PVECGS=3, PVECI=4, PVECIM=5,
     &           PVECIP=6, PNTASK=7, HASBV=8, IGS=9, JSTEP=16)
      integer VECGS, VECI, VECIM1, VECIP1, NTASKS
      double precision EA18ID
      external         EA18ID
      external DLARNV, DGEMM, DSCAL
      external EA16KD, EA16LD, EA18LD, EA18KD, EA16VD, EA16UD
      external EA16MD, EA16ND, EA16OD, EA18JD
      if (IDO.gt.0) then
         VECI = isave(PVECI)
         VECIM1 = isave(PVECIM)
         VECIP1 = isave(PVECIP)
         if (isave(JUMP).eq.500) then
            go to 500
         else if (isave(JUMP).eq.1000) then
            go to 1000
         else if (isave(JUMP).eq.1500) then
            go to 1500
         else if (isave(JUMP).eq.1200) then
            go to 1200
         else if (isave(JUMP).eq.1300) then
            go to 1300
         else if (isave(JUMP).eq.7000) then
            go to 7000
         end if
      end if
      isave(PVECI) = 0
      isave(PVECIM) = 0
      isave(PVECIP) = 0
      isave(HASBV) = 0
      isave(JUMP) = 0
      if (ISTEP.eq.1) then
         if (VINIT.eq.0) then
            do 10 I = NCONV+1,NCONV+BLK
               call DLARNV(2, ISEED, N, V(1,I))
 10         continue
         end if
         K = 1
         ERMOLA(1,SECOND) = ZERO
         ERMOLA(1,THIRD) = ZERO
         rsave(ERRLAN+1) = ZERO
 501     continue
         if (.not.BGS) then
            call EA16KD(N, OFFSET, BLK, V,LDV, COEF, LDCOEF,
     &                NORM, TEMP, LTEMP, RANDOM, SAFMIN, DGKS,
     &                GSMAX*BLK, ISEED, IERR)
            rsave(REFRAT) = ONE
            if (IERR.gt.OFFSET+1) then
               call DLARNV(2, ISEED, N, V(1,IERR))
               go to 501
            else if (IERR.gt.0) then
               IERR = -4
               IDO = 100
               ISTEP = 0
               return
            end if
         else
            isave(JUMP) = 500
            isave(PIDOGS) = 0
            IDO = 2
            NMVB = NMVB + 1
            IPOS(1) = OFFSET+1
            IPOS(2) = OFFSET+BLK
            IPOS(3) = 1
            IPOS(4) = BLK
            return
         end if
      end if
 500  continue
      if (isave(JUMP).eq.500) then
         IDOGS = isave(PIDOGS)
         call EA16LD(IDOGS, IPOS, N, OFFSET, BLK, V,LDV,
     &            BV,LDBV, COEF, LDCOEF, NORM, TEMP, LTEMP,
     &            RANDOM, SAFMIN, DGKS, GSMAX*BLK, ISEED, isave(IGS),
     &            rsave(REFRAT), IERR)
         isave(PIDOGS) = IDOGS
         if (IDOGS.ne.100) then
            isave(JUMP) = 500
            IDO = 2
            NMVB = NMVB + 1
            return
         end if
         if (IERR.gt.OFFSET+1) then
            call DSCAL(N, ZERO, V(1,IERR),1)
            call DSCAL(N, ZERO, BV(1,IERR-OFFSET),1)
            isave(PIDOGS) = 0
            go to 500
         else if (IERR.gt.0) then
            IERR = -4
            ISTEP = 0
            IDO = 100
            return
         end if
         isave(HASBV) = 1
      end if
      VECI = OFFSET+(ISTEP-1)*BLK
      isave(PVECI) = VECI
      VECIM1 = max(OFFSET,VECI-BLK)
      isave(PVECIM) = VECIM1
      VECIP1 = VECI+BLK
      isave(PVECIP) = VECIP1
      if (BGS .and. isave(HASBV).eq.0) then
         IDO = 2
         IPOS(1) = VECI+1
         IPOS(2) = VECI+BLK
         IPOS(3) = 1
         IPOS(4) = BLK
         isave(JUMP) = 1000
         return
      end if
 1000 continue
         IDO = 1
         IPOS(1) = VECI+1
         IPOS(2) = VECI+BLK
         IPOS(3) = VECIP1+1
         IPOS(4) = VECIP1+BLK
         NMVA = NMVA + 1
         isave(JUMP) = 1200
         return
 1200    continue
         if (VECI.ge.OFFSET+BLK) then
            call EA18LD(T(1,VECIM1-OFFSET+1), LDT, BLK, BLK,
     &                  COEF, LDCOEF)
            call DGEMM('Notranspose', 'Transpose', N, BLK, BLK, -ONE,
     &                 V(1,VECIM1+1),LDV, COEF(BLK+1,1),LDCOEF,
     &                 ONE, V(1,VECIP1+1),LDV)
            VECGS = VECI+1
         else
            VECGS = VECIM1+1
         end if
         isave(PVECGS) = VECGS
         if (.not.BGS) then
            call EA16KD(N, VECIP1-VECGS+1, BLK, V(1,VECGS),LDV, COEF,
     &                LDCOEF, NORM, TEMP, LTEMP, RANDOM, SAFMIN,
     &                DGKS, GSMAX, ISEED, IERR)
            rsave(MINRAT) = ONE
            go to 1600
         else
            IDO = 2
            IPOS(1) = IPOS(3)
            IPOS(2) = IPOS(4)
            IPOS(3) = 1
            IPOS(4) = BLK
            NMVB = NMVB + 1
            isave(JUMP) = 1300
            return
         end if
 1300    continue
         isave(PIDOGS) = 0
 1500    continue
         VECGS = isave(PVECGS)
         IDOGS = isave(PIDOGS)
         call EA16LD(IDOGS, IPOS, N, VECIP1-VECGS+1, BLK,
     &            V(1,VECGS),LDV, BV,LDBV, COEF, LDCOEF, NORM, TEMP,
     &            LTEMP, RANDOM, SAFMIN, DGKS, GSMAX, ISEED, isave(IGS),
     &            rsave(MINRAT), IERR)
         isave(PIDOGS) = IDOGS
         if (IDOGS.ne.100) then
            IDO = 2
            IPOS(1) = IPOS(1) + VECGS-1
            IPOS(2) = IPOS(2) + VECGS-1
            NMVB = NMVB + 1
            isave(JUMP) = 1500
            return
         end if
 1600    continue
         isave(HASBV) = 1
         DIAGT = 1
         OFFT = DIAGT+BLK
         call EA18KD(T(1,VECI-OFFSET+1), LDT, BLK, BLK,
     &               COEF(DIAGT,1), LDCOEF)
         if (IERR.lt.0) then
           if (ISTEP.gt.LSTEP+1) then
               ISTEP = ISTEP-1
               NSTEP = ISTEP
               ERRMAX = max(ERRMAX,
     &                  EA18ID(2, 1, NSTEP-1, ERMOLA, LDERLA),
     &                  EA18ID(3, 1, NSTEP, ERMOLA, LDERLA),
     &                  EA18ID(2, 1, NCONV, ERMOVE, LDERVE) )
               IERR = 0
               IDO = 100
               return
            else
               IERR = -11
               if (LUERR.ge.0) write(LUERR,9011) ISTEP
               IDO = 100
               return
            end if
         end if
 9011    format('ERROR : innerproduct was negative at Lanczos',
     &          ' iteration ',I4,/)
         if (IERR.gt.0 .and. IERR.lt.VECIP1-VECGS+2) then
            NSTEP = ISTEP
            ERRMAX = max(ERRMAX,
     &               EA18ID(2, 1, NSTEP-1, ERMOLA, LDERLA),
     &               EA18ID(3, 1, NSTEP, ERMOLA, LDERLA),
     &               EA18ID(2, 1, NCONV, ERMOVE, LDERVE) )
            IERR = -4
            IDO = 100
            return
         end if
         IERR = 0
         isave(JSTEP) = ISTEP
         rsave(ERRLAN) = rsave(ERRLAN+1)
         call EA16VD(ISTEP, isave(JSTEP), BLK, 0, ERMOLA, LDERLA,
     &               COEF(DIAGT,1), LDCOEF,
     &               COEF(OFFT,1), LDCOEF, NORM, TEMP, LTEMP,
     &               SAFMIN, EPSS, rsave(ERRLAN+1), IERR)
         IERR = 0
         call EA16UD(ISTEP, BLK, ERMOLA, LDERLA, ERMOVE, LDERVE,
     &               NCONV, NORM, SAFMIN, ERRLOC)
         if (ERRLOC.le.TOLLOC .and. rsave(ERRLAN+1).le.TOLLAN)
     &   go to 8000
            isave(JSTEP) = ISTEP - 2
            VECGS = VECI-BLK
 7100       continue
               isave(JSTEP) = isave(JSTEP) + 1
               if (isave(JSTEP).gt.ISTEP) go to 8000
               VECGS = VECGS + BLK
               NTASKS = NCONV+1
               do 75 I = 1,NCONV+1
                  GSSTRT(I) = I
                  GSNUMB(I) = 0
 75            continue
               do 77 K = 1,NCONV
                  if (ERMOVE(K,2).gt.TOLLOC .and.
     &                ERMOVE(K,isave(JSTEP)-ISTEP+2).gt.EPSS) then
                     GSNUMB(K) = 1
                     ERMOVE(K,isave(JSTEP)-ISTEP+2) = EPSS
                  end if
 77            continue
               if (rsave(ERRLAN+1).gt.TOLLAN .and.
     &             rsave(ERRLAN+isave(JSTEP)-ISTEP+1).gt.EPSS) then
                  GSNUMB(NTASKS) = isave(JSTEP)*BLK+OFFSET-NCONV
                  do 76 K = 1,isave(JSTEP)
                     ERMOLA(K,isave(JSTEP)-ISTEP+3) = EPSS
 76               continue
                  do 78 K = OFFSET+1,NCONV
                     GSNUMB(K) = 1
                     ERMOVE(K,isave(JSTEP)-ISTEP+2) = EPSS
 78               continue
                  rsave(ERRLAN+isave(JSTEP)-ISTEP+1) = EPSS
               end if
               MTASKS = NTASKS
               call EA16MD(GSSTRT, GSNUMB, MTASKS, NTASKS)
               if (NTASKS.eq.0) go to 7100
               if (.not.BGS) then
                  call EA16ND(N, GSSTRT, GSNUMB, NTASKS, VECGS, BLK,
     &                        V,LDV, NORM, COEF, LDCOEF, TEMP, LTEMP,
     &                        RANDOM, SAFMIN, DGKS, GSMAX, ISEED,
     &                        IERR)
                  go to 7200
               end if
               isave(PNTASK) = NTASKS
               isave(PIDOGS) = 0
               if (isave(JSTEP).lt.ISTEP) isave(HASBV) = 0
               isave(PVECGS) = VECGS
               if (isave(HASBV).eq.0) then
                  IDO = 2
                  NMVB = NMVB + 1
                  IPOS(1) = VECGS+1
                  IPOS(2) = VECGS+BLK
                  IPOS(3) = 1
                  IPOS(4) = BLK
                  isave(JUMP) = 7000
                  return
               end if
 7000          continue
                  VECGS = isave(PVECGS)
                  NTASK = isave(PNTASK)
                  IDOGS = isave(PIDOGS)
                  call EA16OD(IDOGS, IPOS, N, GSSTRT, GSNUMB,
     &                 NTASK, VECGS, BLK, V,LDV, BV,LDBV,
     &                 COEF,LDCOEF, NORM, TEMP, LTEMP, RANDOM, SAFMIN,
     &                 DGKS, GSMAX, ISEED, isave(IGS), rsave(MINRAT),
     &                 IERR)
                  isave(PIDOGS) = IDOGS
               if (isave(PIDOGS).ne.100) then
                  NMVB = NMVB + 1
                  IDO = 2
                  isave(JUMP) = 7000
                  return
               end if
               if (IERR.lt.0) then
                  ISTEP = ISTEP-1
                  NSTEP = ISTEP
                  ERRMAX = max(ERRMAX,
     &                     EA18ID(2, 1, NSTEP-1, ERMOLA, LDERLA),
     &                     EA18ID(3, 1, NSTEP, ERMOLA, LDERLA),
     &                     EA18ID(2, 1, NCONV, ERMOVE, LDERVE) )
                  IERR = 0
                  IDO = 100
                  return
               end if
 7200          continue
               if (VECGS.ge.OFFSET+BLK)
     &            call EA18JD(T(1,VECGS-OFFSET-BLK+1), LDT, BLK,
     &                        COEF(VECGS-BLK+1,1), LDCOEF)
               if (IERR.gt.0) then
                  if (IERR.lt.VECGS+1) ISTEP = ISTEP-1
                  NSTEP = ISTEP
                  ERRMAX = max(ERRMAX,
     &                     EA18ID(2, 1, NSTEP-1, ERMOLA, LDERLA),
     &                     EA18ID(3, 1, NSTEP, ERMOLA, LDERLA),
     &                     EA18ID(2, 1, NCONV, ERMOVE, LDERVE) )
                  IERR = -4
                  IDO = 100
                  return
               end if
               IERR = 0
               isave(HASBV) = 0
               if (isave(JSTEP).eq.ISTEP) isave(HASBV) = 1
            go to 7100
 8000    continue
         ERRMAX = max(ERRMAX,
     &                EA18ID(1, 1, ISTEP-2, ERMOLA, LDERLA),
     &                EA18ID(1, 1, NCONV, ERMOVE, LDERVE) )
         if (rsave(MINRAT).lt.rsave(REFRAT)*RATNOR.and.
     &       ISTEP.gt.LSTEP+1) then
            NSTEP = ISTEP
         end if
         VECIM1 = VECI
         VECI   = VECIP1
         VECIP1 = VECIP1 + BLK
         isave(PVECI) = VECI
         isave(PVECIM) = VECIM1
         isave(PVECIP) = VECIP1
         ISTEP  = ISTEP + 1
      if (ISTEP.le.NSTEP) go to 1000
      ISTEP = NSTEP
      ERRMAX = max(ERRMAX,
     &             EA18ID(2, 1, NSTEP-1, ERMOLA, LDERLA),
     &             EA18ID(3, 1, NSTEP, ERMOLA, LDERLA),
     &             EA18ID(2, 1, NCONV, ERMOVE, LDERVE) )
      IDO = 100
      end
      subroutine EA16FD(T, LDT, NT, BLK, RITZ, EIGVEC, LDEIGV,
     &                  RESID, PIPE, LDPIPE, WORK, LWORK, IERR)
      integer LDT, NT, BLK, LDEIGV, LDPIPE, LWORK, IERR
      double precision T(LDT,NT), RITZ(NT)
      double precision EIGVEC(LDEIGV,NT), RESID(NT), PIPE(LDPIPE,NT)
      double precision WORK(LWORK)
      integer KD, I, J, K, JT, IE
      double precision DNRM2
      external DNRM2
      external DSBEV
      do 10 I = 1,BLK
         do 100 J = 1,BLK
            PIPE(J,I) = T(J+1,NT-BLK+I)
 100     continue
 10   continue
      KD = LDT-1
      KD = min(LDT-1,NT-1)
      call DSBEV('Vectors', 'Lower', NT, KD, T, LDT,
     &           RITZ, EIGVEC, LDEIGV, WORK, IERR)
      if (IERR.ne.0) return
      do 20 I = 1,BLK
         do 200 J = 1,BLK
            T(J+1,NT-BLK+I) = PIPE(J,I)
 200     continue
 20   continue
      do 30 I = 1,NT
         do 300 J = 1,BLK
            PIPE(J,I) = 0.
 300     continue
 30   continue
      do 40 I = 1,BLK
         IE = NT - BLK + I
         JT = BLK - I + 1
         do 400 J = 1,I
            JT = JT + 1
            do 4000 K = 1,NT
               PIPE(J,K) = PIPE(J,K)
     &                       + T(JT,IE) * EIGVEC(IE,K)
 4000       continue
 400     continue
 40   continue
      do 50 I = 1,NT
         RESID(I) = DNRM2(BLK, PIPE(1,I), 1)
 50   continue
      return
      end
      subroutine EA16GD(N, NV, NWANT, BLK, LIWORK, LWORK,
     &                  ISAVE, LISAVE, ICNTL, INFO)
      integer N, NV, BLK, NWANT
      integer LIWORK, LWORK, LISAVE
      integer ICNTL(20), INFO(20), ISAVE(LISAVE)
      integer MRITZ, SHIFT, LILEJA
      integer IERR, MEMO
      parameter (IERR=1, MEMO=2)
      integer DIAGNO, IMSTEP, REDUC, IRSTRT, IMLEJA, IPREFI
      parameter (DIAGNO=2, IMSTEP=4, REDUC=8, IRSTRT=9, IMLEJA=14,
     &           IPREFI = 15)
      integer ILUERR, ILUWAR, ILURES
      parameter (ILUERR=6, ILUWAR=7, ILURES=9)
      integer IF, ILDF, ILCF
      parameter (IF=13, ILDF=14, ILCF=15)
      integer ILCF2, ITEMP3, ILTMP3, IL, ILDL
      parameter (ILCF2=16, ITEMP3=17, ILTMP3=18, IL=19, ILDL=20)
      integer ILCL, ITAU, ILTAU, ILTMP1, ITEMP1
      parameter (ILCL=21, ITAU=22, ILTAU=23, ILTMP1=24, ITEMP1=25)
      integer IRANDM, ILDT, IT, INORM
      parameter (IRANDM=26, ILDT=28, IT=29, INORM=30)
      integer ICOEF, ILDCOE, IPIPE, ITEMP2, INPURG
      parameter (ICOEF=31, ILDCOE=32, IPIPE=33, ITEMP2=34, INPURG=35)
      integer IRESID, IRITZ, IEIGVE
      parameter (IRESID=38, IRITZ=39, IEIGVE=40)
      integer ILDEIG, ILCEIG, ILLSTE
      parameter (ILDEIG=41, ILCEIG=42, ILLSTE=45)
      integer IRMOVE, ILERVE, IRMOLA, ILERLA
      parameter (IRMOVE=46, ILERVE=47, IRMOLA=48, ILERLA=49)
      integer IMLANC, INKEEP
      parameter (IMLANC=51, INKEEP=52)
      integer ILLANC
      parameter (ILLANC=57)
      integer IGSSTR, IGSNUM
      parameter (IGSSTR=58, IGSNUM=59)
      integer ILEJA, IRLEJA
      parameter (ILEJA=67, IRLEJA=69)
      integer IMV, JMSTEP, LMLEJA
      parameter (IMV=70, JMSTEP=87, LMLEJA=100)
      integer LUERR, LUWARN, LUREST
      integer RANDOM, LDT, LCT, T, NORM, LTEMP1, TEMP1
      integer COEF, LDCOEF, PIPE, TEMP2, LTEMP2
      integer RESID, RITZ, EIGVEC, LDEIGV, LCEIGV
      integer F, LDF, LCF, LCF2, TEMP3, LTEMP3, L, LDL, LCL, TAU, LTAU
      integer LLSTEP, LEJA, GSNUMB, GSSTRT, RLEJA
      integer LLANC, MLANC, MV
      integer MISAVE, MRSAVE
      parameter (MISAVE=103, MRSAVE=34)
      integer ERMOVE, ERMOLA, LDERLA, LDERVE
      integer MINNV, NLEJA, JRESTR, MSTEP, MPREF
      external EA16ZD, EA17ED
      LUERR  = ICNTL(1)
      LUWARN = ICNTL(1)
      LUREST = ICNTL(1)
      if (ICNTL(DIAGNO).lt.1) LUERR  = -1
      if (ICNTL(DIAGNO).lt.2) LUWARN = -1
      if (ICNTL(DIAGNO).lt.5) LUREST = -1
      LWORK = 0
      LIWORK = 0
      info(IERR) = 0
      if (N.le.2) then
         if (LUERR.ge.0) then
            write(LUERR,9001) N
         end if
         info(IERR) = -2
         return
      end if
 9001 format('ERROR : N is smaller than 3  (N = ', I5, ').',/)
      if (BLK.lt.1 .or. BLK.gt.N/4) then
         info(IERR) = -5
         if (LUERR.ge.0) then
            write(LUERR,9002) BLK, N/4
         end if
         return
      end if
 9002 format('ERROR : BLK is out of range  (BLK = ', I4, ').',/,
     &      'Ensure BLK is in the range [1,',I5,']',/)
      MV = min(NV,N)
      if (NV.gt.MV) then
         call EA16ZD(info(IERR), 2**0)
         if (LUWARN.ge.0) write(LUWARN,8002) NV, MV
      end if
 8002 format('WARNING : NV is larger than necessary (NV = ', I10, ').'
     &    ,/,'          The user can reduce NV to ', I10,'.',/)
      if (NWANT.le.0 .or. NWANT.gt.N-3*BLK) then
         info(IERR) = -4
         if (LUERR.ge.0) write(LUERR,9007) NWANT, N-3*BLK
         return
      end if
 9007 format('ERROR : NWANT is out of range  (NWANT = ', I10,').',/,
     &       '        Ensure NWANT is in the range [1,',I10,']',/)
      MINNV =  max(MV,((NWANT+BLK-1)/BLK)*BLK+3*BLK)
      if (MV.lt.MINNV) then
         info(MEMO) = MINNV
         info(IERR) = -3
         if (LUERR.ge.0) write(LUERR,9003) NV, MINNV
         return
      end if
 9003 format('ERROR : NV is too small  (NV = ', I10,').',/,
     &       '        Ensure NV is at least ',I10,/)
      MSTEP = icntl(IMSTEP)
      if (MSTEP.lt.2 .or. MSTEP.gt.MV/BLK-1) then
         MSTEP = MV/BLK-1
      end if
      MLANC = MSTEP * BLK
      MRITZ = MV-BLK
      LLSTEP = (MSTEP*icntl(REDUC))/100
      if (icntl(REDUC).lt.1) then
         call EA16ZD(info(IERR), 2**1)
         if (LUWARN.ge.0) then
            write(LUWARN,8011) ICNTL(REDUC), 50
         end if
         LLSTEP = 0
      end if
      if (icntl(REDUC).gt.99) then
         call EA16ZD(info(IERR), 2**1)
         if (LUWARN.ge.0) then
            write(LUWARN,8011) ICNTL(REDUC), 50
         end if
         LLSTEP = MSTEP
      end if
 8011 format('WARNING : ICNTL(8) is out of bounds [1,99].',/,
     &       '          Continue with the default = ', I4, '.', /)
      if (LLSTEP.le.0) LLSTEP = 1
      if (LLSTEP.ge.MSTEP) LLSTEP = MSTEP-1
      LLANC = LLSTEP * BLK
      call EA17ED(1, 1, icntl(IRSTRT),  JRESTR, -1)
      NLEJA = max(6,MINNV/BLK,icntl(IMLEJA))
      LEJA = MISAVE + 1
      LILEJA = 0
      if (JRESTR.ne.3) NLEJA=1
      if (JRESTR.eq.3) LILEJA = 2*NLEJA+2
      LILEJA = 2*NLEJA+2
      RANDOM = LEJA + LILEJA
      GSNUMB = RANDOM + BLK
      GSSTRT = GSNUMB + MV
      LIWORK = GSSTRT + MV
      LIWORK = max(LIWORK,RANDOM-1+MV-BLK)
      RLEJA  = MRSAVE + 1
      LILEJA = 0
      if (JRESTR.eq.3) LILEJA = 3*NLEJA
      LILEJA = 3*NLEJA
      T      = RLEJA + LILEJA
      LDT    = BLK + 1
      LCT    = MLANC + BLK
      RESID  = T + LDT*LCT
      RITZ   = RESID + (MV-BLK)
      ERMOVE = RITZ + MV-BLK
      LDERVE = MRITZ
      ERMOLA = ERMOVE + LDERVE*4
      LDERLA = MSTEP
      NORM   = ERMOLA + 6*LDERLA
      COEF   = NORM + BLK
      LDCOEF = MV
      TEMP1  = COEF + LDCOEF*BLK
      LTEMP1 = max(3*BLK,5*BLK-4,MV*BLK)
      LWORK  = TEMP1 + LTEMP1 - 1
      PIPE   = ERMOVE
      EIGVEC = PIPE + BLK*MLANC
      LDEIGV = MLANC + BLK
      LCEIGV = MLANC + BLK
      TEMP2  = EIGVEC + LDEIGV*LCEIGV
      LTEMP2 = max(BLK,3*MLANC-2,MV-BLK)
      LWORK  = max(LWORK,TEMP2+LTEMP2-1)
      MPREF = MLANC*min(MSTEP-1, max(0,ICNTL(IPREFI)))
      F = EIGVEC + LDEIGV*LCEIGV
      LDF = MLANC + BLK
      LCF = MLANC + BLK
      LCF2 = MLANC + BLK
      SHIFT = F
      L = max(F + LDF*LCF, F + BLK*MPREF, SHIFT+MLANC)
      LDL = max(3*BLK+1,2*BLK+3)
      LCL = MLANC
      TAU = max( L + LDL*LCL, F + LDF*abs(ICNTL(IPREFI)) )
      LTAU = max(MLANC+BLK, MPREF)
      TEMP3 = max(F + LDF*LCF2, F + MV-BLK, TAU + LTAU)
      LWORK = max(LWORK, TEMP3+MLANC+BLK-1, TEMP3+MV-BLK-1,
     &            TEMP3+3*MLANC-2-1)
      LTEMP3 = LWORK - TEMP3 + 1
      if (LISAVE.lt.MISAVE) return
      ISAVE(ILUERR) = LUERR
      ISAVE(ILUWAR) = LUWARN
      ISAVE(ILURES) = LUREST
      ISAVE(IF)     = F
      ISAVE(ILDF)   = LDF
      ISAVE(ILCF)   = LCF
      ISAVE(ILCF2)  = LCF2
      ISAVE(ITEMP3) = TEMP3
      ISAVE(ILTMP3) = LTEMP3
      ISAVE(IL)     = L
      ISAVE(ILDL)   = LDL
      ISAVE(ILCL)   = LCL
      ISAVE(ITAU)   = TAU
      ISAVE(ILTAU)  = LTAU
      ISAVE(IRANDM) = RANDOM
      ISAVE(ILDT)   = LDT
      ISAVE(IT)     = T
      ISAVE(INORM)  = NORM
      ISAVE(ILTMP1) = LTEMP1
      ISAVE(ITEMP1) = TEMP1
      ISAVE(ICOEF)  = COEF
      ISAVE(ILDCOE) = LDCOEF
      ISAVE(IPIPE)  = PIPE
      ISAVE(ITEMP2) = TEMP2
      ISAVE(IRESID) = RESID
      ISAVE(IRITZ)  = RITZ
      ISAVE(IEIGVE) = EIGVEC
      ISAVE(ILDEIG) = LDEIGV
      ISAVE(ILCEIG) = LCEIGV
      ISAVE(ILLSTE) = LLSTEP
      ISAVE(IRMOVE) = ERMOVE
      ISAVE(ILERVE) = LDERVE
      ISAVE(IRMOLA) = ERMOLA
      ISAVE(ILERLA) = LDERLA
      ISAVE(IMLANC) = MLANC
      ISAVE(INKEEP) = 0
      ISAVE(INPURG) = 0
      ISAVE(ILLANC) = LLANC
      ISAVE(IMV)    = MV
      ISAVE(ILEJA)  = LEJA
      ISAVE(IGSNUM) = GSNUMB
      ISAVE(IGSSTR) = GSSTRT
      ISAVE(IRLEJA) = RLEJA
      ISAVE(JMSTEP) = MSTEP
      ISAVE(LMLEJA) = NLEJA
      return
      end
      subroutine EA16HD(NCONV, NCOOLD, MRITZ, NLANC, PERM, RITZ, RESID,
     &                  V, LDV, N, NVECS, LU, NCONV2, NCOOL2)
      integer NCONV, NCOOLD, MRITZ, NLANC, LDV, N, NVECS, LU
      integer NCONV2, NCOOL2
      integer PERM(NCOOLD)
      double precision RITZ(NCOOLD+NLANC), RESID(NCOOLD+NLANC)
      double precision V(LDV,NVECS)
      integer IRITZ, IPOS, TOREM, NOLD
      double precision ZERO, ONE
      parameter       (ZERO=0.0D0, ONE=1.0D0)
      external DCOPY
      TOREM = NCONV-MRITZ
      do 10 IRITZ = NCOOLD-TOREM+1,NCOOLD
         IPOS = PERM(IRITZ)
         if (IPOS.le.NCOOLD) then
            RESID(IPOS) = -ONE
         end if
 10   continue
      NOLD = 0
      do 20 IRITZ = 1,NCOOLD
         if (NOLD.gt.0) then
            RITZ(IRITZ-NOLD) = RITZ(IRITZ)
            RESID(IRITZ-NOLD) = RESID(IRITZ)
            call DCOPY(N, V(1,IRITZ),1, V(1,IRITZ-NOLD),1)
         end if
         if (RESID(IRITZ).lt.ZERO) then
            NOLD = NOLD + 1
         end if
 20   continue
      if (NOLD.gt.0) then
         do 30 IRITZ = NCOOLD+1,NCOOLD+NLANC
            RITZ(IRITZ-NOLD) = RITZ(IRITZ)
            RESID(IRITZ-NOLD) = RESID(IRITZ)
            call DCOPY(N, V(1,IRITZ),1, V(1,IRITZ-NOLD),1)
 30      continue
         do 40 IRITZ = NCOOLD+NLANC+1,NVECS
            call DCOPY(N, V(1,IRITZ),1, V(1,IRITZ-NOLD),1)
 40      continue
      end if
      if (LU.ge.0) write(LU,9000) NOLD
 9000 format('  Removing ', I4, ' unwanted locked Ritz vectors',/)
      NCOOL2 = NCOOLD - NOLD
      NCONV2 = NCONV - NOLD
      return
      end
      subroutine EA16ID(ICNTL, CNTL)
      integer ICNTL(20)
      double precision CNTL(15)
      double precision EPS
      double precision DLAMCH
      external DLAMCH
      intrinsic sqrt
      EPS = DLAMCH('EPSILON')
      CNTL(1) = sqrt(EPS)
      CNTL(2) = 0.
      CNTL(3) = sqrt(EPS)
      CNTL(4) = 500.0D0
      CNTL(5) = sqrt(sqrt(EPS))
      ICNTL(1) = 6
      ICNTL(2) = 2
      ICNTL(3) = 0
      ICNTL(4) = 0
      ICNTL(5) = 0
      ICNTL(6) = 2
      ICNTL(7) = 0
      ICNTL(8) = 50
      ICNTL(9) = 1
      ICNTL(10) = 0
      ICNTL(11) = 0
      ICNTL(12) = 0
      ICNTL(13) = 1
      ICNTL(14) = 100
      ICNTL(15) = 0
      ICNTL(16) = 0
      ICNTL(17) = 200
      ICNTL(18) = 1
      ICNTL(19) = 0
      ICNTL(20) = 0
      return
      end
      subroutine EA16JD(IDO, MODE, IPOS, N, BLK, VINIT, INIVEC, ISEED,
     &                  V,LDV, ISTEP, NMVA, NMVB, LU)
      integer IDO, MODE, N, BLK, VINIT, INIVEC, LDV, ISTEP
      integer NMVA, NMVB, LU
      integer IPOS(4), ISEED(4)
      double precision V(LDV,2*BLK)
      integer I
      external DLARNV, DLACPY
      if (VINIT.eq.0) then
         do 10 I = 1,BLK
            call DLARNV(2, ISEED, N, V(1,I))
 10      continue
         VINIT = 1
      end if
      if (MODE.lt.4) then
         IDO = 100
         return
      end if
      if (IDO.eq.0) then
         ISTEP = 0
         if (INIVEC.gt.0 .and. LU.ge.0) write(LU,9000) INIVEC
 9000    format('Pre-processing the initial vector (',I2,' steps).',/)
      end if
      if (ISTEP.ge.INIVEC) then
         IDO = 100
         return
      end if
      if (IDO.ne.2 .and. MODE.ge.4) then
         IDO = 2
         IPOS(1) = 1
         IPOS(2) = BLK
         IPOS(3) = 1
         IPOS(4) = BLK
         NMVB = NMVB + 1
         return
      end if
      call DLACPY('All', N, BLK, V,LDV, V(1,BLK+1),LDV)
      IDO = 1
      IPOS(1) = BLK+1
      IPOS(2) = 2*BLK
      IPOS(3) = 1
      IPOS(4) = BLK
      NMVA = NMVA + 1
      ISTEP = ISTEP + 1
      return
      end
      subroutine EA16KD(N, K, M, V, LDV, COEF, LDCOEF, NORM, TEMP,
     &                  LTEMP, RANDOM, UDFL, DGKS, MAXIT, ISEED, INFO)
      integer N, K, M, LDV, LDCOEF, MAXIT, LTEMP, INFO
      integer RANDOM(M), ISEED(4)
      double precision UDFL, DGKS
      double precision V(LDV,K+M), COEF(LDCOEF,M), NORM(M)
      double precision TEMP(LTEMP)
      double precision ONE, ZERO
      parameter (ZERO=0.0D0, ONE=1.0D0)
      double precision VNORM, VINVNR
      integer I, J, IJ, II
      integer KK, MM, BASE, ITER, REORTH
      double precision DNRM2
      external DNRM2
      external DLASET, DGEMM, DGER, DLARNV, DSCAL, DGEMV, DRSCL
      KK = K
      MM = M
      BASE = 0
      ITER = 0
      do 15 I = 1,M
         do 150 J = 1,K+M
            COEF(J,I) = ZERO
 150     continue
         COEF(K+I,I) = ONE
 15   continue
      do 25 I = 1,M
         RANDOM(I) = 1
 25   continue
      INFO = 0
 1000 continue
      do 35 I = 1,MM
         NORM(BASE+I) = DNRM2(N, V(1,KK+I),1)
 35   continue
 1100 continue
         if (ITER.ge.MAXIT) then
            INFO = KK+1
            call DLASET('All', MM, MM, ZERO, ZERO,
     &                  COEF(KK+1,BASE+1),LDCOEF)
            return
         end if
         ITER = ITER + 1
         if (KK.gt.0) then
            call DGEMM('T', 'N', KK, MM, N, ONE, V,LDV, V(1,KK+1),LDV,
     &                 ZERO, TEMP,KK)
            call DGEMM('N', 'N', N, MM, KK, -ONE, V,LDV, TEMP,KK, ONE,
     &                 V(1,KK+1),LDV)
            IJ = 1
            do 30 I = BASE+1,M
               VNORM = COEF(K+I,I)
               do 300 J = 1,KK
                  COEF(J,I) = COEF(J,I) + VNORM*TEMP(IJ)
                  IJ = IJ + 1
 300           continue
 30         continue
         end if
         REORTH = KK
         do 40 I = 1,MM
            VNORM = DNRM2(N, V(1,KK+I),1)
            if (VNORM.eq.ZERO) then
               if (I.lt.MM) then
                  call DGER(N, MM-I, ONE, V(1,KK+I),1,
     &                      COEF(KK+I,BASE+I+1),LDCOEF, V(1,KK+I+1),LDV)
                  call DGER(KK+I, MM-I, ONE, COEF(1,BASE+I),1,
     &                      COEF(KK+I,BASE+I+1),LDCOEF,
     &                      COEF(1,BASE+I+1),LDCOEF)
               end if
               call DLARNV(2, ISEED, N, V(1,KK+I))
               RANDOM(BASE+I) = 0
               COEF(KK+I,BASE+I) = ZERO
               ITER = 0
               KK = REORTH
               BASE = KK - K
               MM = M - BASE
               go to 1000
            end if
            VINVNR = VNORM
            IJ = BASE+I
            COEF(KK+I,IJ) = COEF(KK+I,IJ)*VINVNR
            TEMP(BASE+I) = VINVNR
            if (VINVNR.ge.UDFL) then
               VINVNR = ONE/VINVNR
               call DSCAL(N, VINVNR, V(1,KK+I),1)
            else
               call DRSCL(N, VINVNR, V(1,KK+I),1)
            end if
            if (I.lt.MM) then
               call DGEMV('T', N, MM-I, ONE, V(1,KK+I+1),LDV,
     &                    V(1,KK+I),1, ZERO, TEMP(BASE+I+1),1)
               do 430 II = BASE+I+1,M
                  if (RANDOM(II).eq.1) then
                     COEF(KK+I,II) = COEF(KK+I,II) * TEMP(BASE+I)
                     do 4300 J = BASE+I+1,II
                        COEF(KK+I,II) = COEF(KK+I,II)
     &                                + TEMP(J) * COEF(K+J,II)
 4300                continue
                  end if
 430           continue
               call DGER(N, MM-I, -ONE, V(1,KK+I),1, TEMP(BASE+I+1),1,
     &                   V(1,KK+I+1),LDV)
            end if
            if (VNORM.ge.DGKS*NORM(BASE+I) .and. REORTH.eq.KK+I-1) then
               REORTH = KK+I
            end if
            NORM(BASE+I) = ONE
 40      continue
      if (REORTH.lt.K+M) then
         KK = REORTH
         BASE = KK - K
         MM = M - BASE
         go to 1100
      end if
      return
      end
      subroutine EA16LD(IDO, IPOS, N, K, M, V, LDV, BV, LDBV, COEF,
     &                LDCOEF, NORM, TEMP, LTEMP, RANDOM, UDFL, DGKS,
     &                MAXIT, ISEED, ISAVE, MINRAT, INFO)
      integer IDO, N, K, M, LDV, LDBV, LDCOEF, MAXIT, LTEMP, INFO
      integer IPOS(4), RANDOM(M), ISEED(4), ISAVE(6)
      double precision UDFL, DGKS, MINRAT
      double precision V(LDV,K+M), BV(LDBV,M), COEF(LDCOEF,M), NORM(M)
      double precision TEMP(LTEMP)
      double precision ONE, ZERO
      parameter (ZERO=0.0D0, ONE=1.0D0)
      double precision VNORM, VINVNR, RATIO
      integer I, J, IJ, II
      integer   IJUMP, IKK, IMM, IBASE, IRESTA, IREORT
      parameter (IJUMP=1, IKK=2, IMM=3, IBASE=4, IRESTA=5, IREORT=6)
      integer   JUMP, KK, MM, BASE, ITER, REORTH
      external EA18ND, DLASET, DGEMM, DGER, DLARNV, DGEMV, EA17JD
      if (IDO.eq.0) then
         KK = K
         MM = M
         BASE = 0
         MINRAT = ONE/UDFL
         ITER = 0
         do 15 I = 1,M
            do 150 J = 1,K+M
               COEF(J,I) = ZERO
 150        continue
            COEF(K+I,I) = ONE
 15      continue
         do 25 I = 1,M
            RANDOM(I) = 1
 25      continue
         INFO = 0
         REORTH = 0
         JUMP = 0
         IDO = 1
      else
         JUMP = isave(IJUMP)
         KK = isave(IKK)
         MM = isave(IMM)
         BASE = isave(IBASE)
         ITER = isave(IRESTA)
         REORTH = isave(IREORT)
         if (JUMP.eq.1000) then
            go to 1000
         else if (JUMP.eq.2000) then
            go to 2000
         end if
      end if
 1000 continue
      do 35 I = 1,MM
         call EA18ND(N, BV(1,BASE+I), V(1,KK+I), NORM(BASE+I),
     &               RATIO, UDFL, INFO)
         if (INFO.ne.0) then
            INFO = -(BASE+I)
            IDO = 100
            go to 8000
         end if
 35   continue
 1100 continue
      if (ITER.ge.MAXIT) then
         IDO = 100
         INFO = KK+1
         call DLASET('All', MM, MM, ZERO, ZERO,
     &               COEF(KK+1,BASE+1),LDCOEF)
         go to 8000
      end if
      ITER = ITER + 1
      if (KK.gt.0) then
      call DGEMM('T', 'N', KK, MM, N, ONE, V,LDV, BV(1,BASE+1),LDBV,
     &           ZERO, TEMP,KK)
      call DGEMM('N', 'N', N, MM, KK, -ONE, V,LDV, TEMP,KK, ONE,
     &           V(1,KK+1),LDV)
      IJ = 1
      do 30 I = BASE+1,M
         VNORM = COEF(K+I,I)
         do 300 J = 1,KK
            COEF(J,I) = COEF(J,I) + VNORM*TEMP(IJ)
            IJ = IJ + 1
 300     continue
 30   continue
      IPOS(1) = KK+1
      IPOS(2) = KK+MM
      IPOS(3) = BASE+1
      IPOS(4) = BASE+MM
      IDO = 1
      JUMP = 2000
      go to 8000
      end if
 2000 continue
      REORTH = KK
      do 40 I = 1,MM
         call EA18ND(N, BV(1,BASE+I), V(1,KK+I), VNORM,
     &               RATIO, UDFL, INFO)
         if (INFO.ne.0) then
            IDO = 100
            INFO = -(BASE+I)
            go to 8000
         end if
         if (VNORM.eq.ZERO) then
            if (I.lt.MM) then
               call DGER(N, MM-I, ONE, V(1,KK+I),1,
     &                   COEF(KK+I,BASE+I+1),LDCOEF, V(1,KK+I+1),LDV)
               call DGER(KK+I, MM-I, ONE, COEF(1,BASE+I),1,
     &                   COEF(KK+I,BASE+I+1),LDCOEF,
     &                   COEF(1,BASE+I+1),LDCOEF)
            end if
            call DLARNV(2, ISEED, N, V(1,KK+I))
            RANDOM(BASE+I) = 0
            COEF(KK+I,BASE+I) = ZERO
            ITER = 0
            IDO = 1
            IPOS(1) = KK+I
            IPOS(2) = KK+I
            IPOS(3) = BASE+I
            IPOS(4) = BASE+I
            JUMP = 1000
            KK = REORTH
            BASE = KK - K
            MM = M - BASE
            go to 8000
         end if
         MINRAT = min(MINRAT, RATIO)
         VINVNR = VNORM
         IJ = BASE+I
         COEF(KK+I,IJ) = COEF(KK+I,IJ)*VINVNR
         TEMP(BASE+I) = VINVNR
         call EA17JD(N, V(1,KK+I), VINVNR, UDFL, INFO)
         if (INFO.ne.0) then
            INFO = M+1
            IDO = 100
            return
         end if
         call EA17JD(N, BV(1,BASE+I), VINVNR, UDFL, INFO)
         if (INFO.ne.0) then
            INFO = M+1
            IDO = 100
            return
         end if
         if (I.lt.MM) then
            call DGEMV('T', N, MM-I, ONE, V(1,KK+I+1),LDV,
     &                 BV(1,BASE+I),1, ZERO, TEMP(BASE+I+1),1)
            do 430 II = BASE+I+1,M
               if (RANDOM(II).eq.1) then
                  COEF(KK+I,II) = COEF(KK+I,II) * TEMP(BASE+I)
                  do 4300 J = BASE+I+1,II
                     COEF(KK+I,II) = COEF(KK+I,II)
     &                             + TEMP(J) * COEF(K+J,II)
 4300             continue
               end if
 430        continue
            call DGER(N, MM-I, -ONE, V(1,KK+I),1, TEMP(BASE+I+1),1,
     &                V(1,KK+I+1),LDV)
            call DGER(N, MM-I, -ONE, BV(1,BASE+I),1, TEMP(BASE+I+1),1,
     &                BV(1,BASE+I+1),LDBV)
         end if
         if (VNORM.ge.DGKS*NORM(BASE+I) .and.
     &       REORTH.eq.KK+I-1) then
            REORTH = KK+I
         end if
         NORM(BASE+I) = ONE
 40   continue
      if (REORTH.lt.K+M) then
         KK = REORTH
         BASE = KK - K
         MM = M - BASE
         go to 1100
      end if
      IDO = 100
 8000 continue
      isave(IJUMP)  = JUMP
      isave(IKK)    = KK
      isave(IMM)    = MM
      isave(IBASE)  = BASE
      isave(IRESTA)  = ITER
      isave(IREORT) = REORTH
      return
      end
      subroutine EA16MD(FIRST, NUMBER, NNUMB, NTASKS)
      integer NTASKS, NNUMB
      integer FIRST(NNUMB), NUMBER(NNUMB)
      logical NEW
      integer ITASK, NT
      NT = 0
      do 10 ITASK = 1,NTASKS
         if (NUMBER(ITASK).eq.0) go to 10
         if (NT.eq.0) then
            NEW = .true.
         else
            NEW = FIRST(NT)+NUMBER(NT).lt.FIRST(ITASK)
         end if
         if (NEW) then
            NT = NT + 1
            FIRST(NT) = FIRST(ITASK)
            NUMBER(NT) = NUMBER(ITASK)
         else
            NUMBER(NT) = NUMBER(NT) + NUMBER(ITASK)
         end if
 10   continue
      NTASKS = NT
      return
      end
      subroutine EA16ND(N, FIRST, NUMBER, NTASKS, K, M, V,
     &                 LDV, NORM, COEF, LDCOEF, TEMP, LTEMP,
     &                 RANDOM, UDFL, DGKS, MAXIT, ISEED, INFO)
      integer N, NTASKS, K, M, LDV, LDCOEF, MAXIT, LTEMP
      integer INFO
      integer FIRST(NTASKS+1), NUMBER(NTASKS+1), RANDOM(M)
      integer ISEED(4)
      double precision UDFL, DGKS
      double precision V(LDV,K+M), COEF(LDCOEF,M), NORM(M)
      double precision TEMP(LTEMP)
      double precision ONE, ZERO
      parameter (ZERO=0.0D0, ONE=1.0D0)
      double precision VNORM, VINVNR
      integer I, J, IJ, II, ITASK, NN, FF
      integer KK, MM, BASE, ITER, REORTH, NT
      double precision DNRM2
      external DNRM2
      external DGEMM, DGER, DLARNV, DSCAL, DRSCL, DGEMV
      KK = K
      MM = M
      BASE = 0
      NT = NTASKS
      if (FIRST(NT)+NUMBER(NT).lt.K) then
         NT = NT + 1
         FIRST(NT) = K+1
         NUMBER(NT) = 0
      end if
      ITER = 0
      do 15 I = 1,M
         do 150 J = 1,KK+MM
            COEF(J,I) = ZERO
 150     continue
         COEF(K+I,I) = ONE
 15   continue
      do 25 I = 1,MM
         RANDOM(I) = 1
 25   continue
      INFO = 0
 1000 continue
      do 35 I = 1,MM
         NORM(BASE+I) = DNRM2(N, V(1,KK+I),1)
 35   continue
 1100 continue
      if (ITER.ge.MAXIT) then
         INFO = KK+1
         return
      end if
      ITER = ITER + 1
      do 30 ITASK = 1,NT
         NN = NUMBER(ITASK)
         FF = FIRST(ITASK)
         if (NN.gt.0) then
            call DGEMM('T', 'N', NN, MM, N, ONE,
     &                 V(1,FF),LDV, V(1,KK+1),LDV,
     &                 ZERO, TEMP,NN)
            call DGEMM('N', 'N', N, MM, NN, -ONE, V(1,FF),LDV, TEMP,NN,
     &                 ONE, V(1,KK+1),LDV)
            IJ = 1
            do 300 I = BASE+1,M
               VNORM = COEF(K+I,I)
               do 3000 J = FF,FF+NN-1
                  COEF(J,I) = COEF(J,I) + VNORM*TEMP(IJ)
                  IJ = IJ + 1
 3000          continue
 300        continue
         end if
 30   continue
      REORTH = KK
      do 40 I = 1,MM
         VNORM = DNRM2(N, V(1,KK+I),1)
         if (VNORM.eq.ZERO) then
            if (I.lt.MM) then
               call DGER(N, MM-I, ONE, V(1,KK+I),1,
     &                   COEF(KK+I,BASE+I+1),LDCOEF, V(1,KK+I+1),LDV)
               call DGER(KK+I, MM-I, ONE, COEF(1,BASE+I),1,
     &                   COEF(KK+I,BASE+I+1),LDCOEF,
     &                   COEF(1,BASE+I+1),LDCOEF)
            end if
            call DLARNV(2, ISEED, N, V(1,KK+I))
            RANDOM(BASE+I) = 0
            COEF(KK+I,BASE+I) = ZERO
            ITER = 0
            go to 1000
         end if
         VINVNR = VNORM
         IJ = BASE+I
         COEF(KK+I,IJ) = COEF(KK+I,IJ)*VINVNR
         TEMP(BASE+I) = VINVNR
         if (VINVNR.ge.UDFL) then
            VINVNR = ONE/VINVNR
            call DSCAL(N, VINVNR, V(1,KK+I),1)
         else
            call DRSCL(N, VINVNR, V(1,KK+I),1)
         end if
         if (I.lt.MM) then
            call DGEMV('T', N, MM-I, ONE, V(1,KK+I+1),LDV,
     &                 V(1,KK+I),1, ZERO, TEMP(BASE+I+1),1)
            do 430 II = BASE+I+1,M
               if (RANDOM(II).eq.1) then
                  COEF(KK+I,II) = COEF(KK+I,II) * TEMP(BASE+I)
                  do 4300 J = BASE+I+1,II
                     COEF(KK+I,II) = COEF(KK+I,II)
     &                             + TEMP(J) * COEF(K+J,II)
 4300             continue
               end if
 430        continue
            call DGER(N, MM-I, -ONE, V(1,KK+I),1, TEMP(BASE+I+1),1,
     &                V(1,KK+I+1),LDV)
         end if
         if (VNORM.ge.DGKS*NORM(BASE+I) .and. REORTH.eq.KK+I-1) then
            REORTH = KK+I
            ITER = 0
         end if
         NORM(BASE+I) = ONE
 40   continue
      if (REORTH.lt.K+M) then
         NUMBER(NT) = NUMBER(NT) + REORTH-KK
         KK = REORTH
         BASE = KK - K
         MM = M - BASE
         go to 1100
      end if
      return
      end
      subroutine EA16OD(IDO, IPOS, N, FIRST, NUMBER, NTASKS, K, M, V,
     &                 LDV, BV, LDBV, COEF, LDCOEF, NORM, TEMP, LTEMP,
     &                 RANDOM, UDFL, DGKS, MAXIT, ISEED, ISAVE, MINRAT,
     &                 INFO)
      integer IDO, N, NTASKS, K, M, LDV, LDBV, LDCOEF, MAXIT, LTEMP
      integer INFO
      integer IPOS(4), FIRST(NTASKS+1), NUMBER(NTASKS+1), RANDOM(M)
      integer ISEED(4), ISAVE(7)
      double precision UDFL, DGKS, MINRAT
      double precision V(LDV,K+M), BV(LDBV,M), COEF(LDCOEF,M), NORM(M)
      double precision TEMP(LTEMP)
      double precision ONE, ZERO
      parameter (ZERO=0.0D0, ONE=1.0D0)
      double precision VNORM, VINVNR, RATIO
      integer I, J, IJ, II, ITASK, NN, FF
      integer JUMP, KK, MM, BASE, ITER, REORTH, NT
      integer IJUMP, IKK, IMM, IBASE, IRESTA, IREORT, INT
      parameter (IJUMP=1, IKK=2, IMM=3, IBASE=4, IRESTA=5, IREORT=6,
     &           INT=7)
      external EA18ND, DGEMM, DGER, DLARNV, DGEMV, EA17JD
      if (IDO.eq.0) then
         KK = K
         MM = M
         BASE = 0
         MINRAT = ONE/UDFL
         NT = NTASKS
         if (FIRST(NT)+NUMBER(NT).lt.K) then
            NT = NT + 1
            FIRST(NT) = K+1
            NUMBER(NT) = 0
         end if
         ITER = 0
         do 15 I = 1,M
            do 150 J = 1,KK+MM
               COEF(J,I) = ZERO
 150        continue
            COEF(K+I,I) = ONE
 15      continue
         do 25 I = 1,MM
            RANDOM(I) = 1
 25      continue
         INFO = 0
         JUMP = 0
         REORTH = 0
         IDO = 1
      else
         JUMP   = isave(IJUMP)
         KK     = isave(IKK)
         MM     = isave(IMM)
         BASE   = isave(IBASE)
         ITER   = isave(IRESTA)
         REORTH = isave(IREORT)
         NT     = isave(INT)
         if (JUMP.eq.1000) then
            go to 1000
         else if (JUMP.eq.2000) then
            go to 2000
         end if
      end if
 1000 continue
      do 35 I = 1,MM
         call EA18ND(N, BV(1,BASE+I), V(1,KK+I), NORM(BASE+I),
     &               RATIO, UDFL, INFO)
         if (INFO.ne.0) then
            INFO = -(BASE+I)
            IDO = 100
            go to 8000
         end if
 35   continue
 1100 continue
      if (ITER.ge.MAXIT) then
         IDO = 100
         INFO = KK+1
         go to 8000
      end if
      ITER = ITER + 1
      if (KK.gt.0) then
      do 30 ITASK = 1,NT
         NN = NUMBER(ITASK)
         FF = FIRST(ITASK)
         if (NN.gt.0) then
            call DGEMM('T', 'N', NN, MM, N, ONE,
     &                 V(1,FF),LDV, BV(1,BASE+1),LDBV,
     &                 ZERO, TEMP,NN)
            call DGEMM('N', 'N', N, MM, NN, -ONE, V(1,FF),LDV, TEMP,NN,
     &                 ONE, V(1,KK+1),LDV)
            IJ = 1
            do 300 I = BASE+1,M
               VNORM = COEF(K+I,I)
               do 3000 J = FF,FF+NN-1
                  COEF(J,I) = COEF(J,I) + VNORM*TEMP(IJ)
                  IJ = IJ + 1
 3000          continue
 300        continue
         end if
 30   continue
      IPOS(1) = KK+1
      IPOS(2) = KK+MM
      IPOS(3) = BASE+1
      IPOS(4) = BASE+MM
      IDO = 1
      JUMP = 2000
      go to 8000
      end if
 2000 continue
      REORTH = KK
      do 40 I = 1,MM
         call EA18ND(N, BV(1,BASE+I), V(1,KK+I), VNORM,
     &               RATIO, UDFL, INFO)
         if (INFO.ne.0) then
            INFO = -(BASE+I)
            IDO = 100
            go to 8000
         end if
         if (VNORM.eq.ZERO) then
            if (I.lt.MM) then
               call DGER(N, MM-I, ONE, V(1,KK+I),1,
     &                   COEF(KK+I,BASE+I+1),LDCOEF, V(1,KK+I+1),LDV)
               call DGER(KK+I, MM-I, ONE, COEF(1,BASE+I),1,
     &                   COEF(KK+I,BASE+I+1),LDCOEF,
     &                   COEF(1,BASE+I+1),LDCOEF)
            end if
            call DLARNV(2, ISEED, N, V(1,KK+I))
            RANDOM(BASE+I) = 0
            COEF(KK+I,BASE+I) = ZERO
            ITER = 0
            IDO = 1
            IPOS(1) = KK+I
            IPOS(2) = KK+I
            IPOS(3) = BASE+I
            IPOS(4) = BASE+I
            JUMP = 1000
            KK = REORTH
            BASE = KK - K
            MM = M - BASE
            NUMBER(NT) = NUMBER(NT) + REORTH-KK
            go to 8000
         end if
         MINRAT = min(MINRAT, RATIO)
         VINVNR = VNORM
         IJ = BASE+I
         COEF(KK+I,IJ) = COEF(KK+I,IJ)*VINVNR
         TEMP(BASE+I) = VINVNR
         call EA17JD(N, V(1,KK+I), VINVNR, UDFL, INFO)
         if (INFO.ne.0) then
            INFO = M+1
            IDO = 100
            return
         end if
         call EA17JD(N, BV(1,BASE+I), VINVNR, UDFL, INFO)
         if (INFO.ne.0) then
            INFO = M+1
            IDO = 100
            return
         end if
         if (I.lt.MM) then
            call DGEMV('T', N, MM-I, ONE, V(1,KK+I+1),LDV,
     &                 BV(1,BASE+I),1, ZERO, TEMP(BASE+I+1),1)
            do 430 II = BASE+I+1,M
               if (RANDOM(II).eq.1) then
                  COEF(KK+I,II) = COEF(KK+I,II) * TEMP(BASE+I)
                  do 4300 J = BASE+I+1,II
                     COEF(KK+I,II) = COEF(KK+I,II)
     &                             + TEMP(J) * COEF(K+J,II)
 4300             continue
               end if
 430        continue
            call DGER(N, MM-I, -ONE, V(1,KK+I),1, TEMP(BASE+I+1),1,
     &                V(1,KK+I+1),LDV)
            call DGER(N, MM-I, -ONE, BV(1,BASE+I),1, TEMP(BASE+I+1),1,
     &                BV(1,BASE+I+1),LDBV)
         end if
         if (VNORM.ge.DGKS*NORM(BASE+I) .and.
     &       REORTH.eq.KK+I-1) then
            REORTH = KK+I
         end if
         NORM(BASE+I) = ONE
 40   continue
      if (REORTH.lt.K+M) then
         NUMBER(NT) = NUMBER(NT) + REORTH-KK
         KK = REORTH
         BASE = KK - K
         MM = M - BASE
         go to 1100
      end if
      IDO = 100
 8000 continue
      isave(IJUMP)  = JUMP
      isave(IKK)    = KK
      isave(IMM)    = MM
      isave(IBASE)  = BASE
      isave(IRESTA)  = ITER
      isave(IREORT) = REORTH
      isave(INT)    = NT
      return
      end
      subroutine EA16PD(NCOOLD, NCONV, NKRYL, RITZ, RESID,
     &                  EIGVEC, LDEIGV, NEIGV, PIPE, LDPIPE, TOL,
     &                  TOLLAN, PERM, TEMP, EPSS)
      integer NCOOLD, NCONV, NKRYL, LDEIGV, NEIGV, LDPIPE
      integer PERM(NKRYL-NCOOLD)
      double precision TOLLAN, EPSS
      double precision RITZ(NKRYL), RESID(NKRYL)
      double precision TOL(2), EIGVEC(LDEIGV,NKRYL-NCOOLD)
      double precision PIPE(LDPIPE,NKRYL-NCOOLD), TEMP(NKRYL-NCOOLD)
      integer I, NVAL
      double precision NORMT, RTOL, MAXRES, ATOL
      double precision ZERO
      parameter (ZERO=0.0D0)
      external EA18OD
      MAXRES = ZERO
      do 30 I = 1,NCOOLD
         MAXRES = max(MAXRES, RESID(I))
 30   continue
      NORMT = ZERO
      do 35 I = NCOOLD+1,NKRYL
         NORMT = max(NORMT, abs(RITZ(I)))
 35   continue
      ATOL = abs(TOL(1)) + NORMT*EPSS
      RTOL = abs(TOL(2))
      NCONV = 0
      NVAL = NCOOLD
      do 20 I = NCOOLD+1,NKRYL
         if (RESID(I) .le. ATOL + RTOL*abs(RITZ(I))) then
            NCONV = NCONV + 1
            PERM(NCONV) = I - NCOOLD
         else
            PERM(NKRYL-NVAL) = I - NCOOLD
            NVAL = NVAL + 1
         end if
 20   continue
      NCONV = NCOOLD + NCONV
      call EA18OD(NKRYL-NCOOLD, RITZ(NCOOLD+1), 1, 1, PERM, TEMP)
      call EA18OD(NKRYL-NCOOLD, RESID(NCOOLD+1), 1, 1, PERM, TEMP)
      call EA18OD(NKRYL-NCOOLD, EIGVEC, LDEIGV, NEIGV, PERM, TEMP)
      call EA18OD(NKRYL-NCOOLD, PIPE,LDPIPE, LDPIPE, PERM, TEMP)
      MAXRES = max(MAXRES, NORMT*max(abs(TOLLAN),EPSS))
      do 40 I = NCOOLD+1,NCONV
         RESID(I) = max(RESID(I), MAXRES)
 40   continue
      return
      end
      subroutine EA16QD(NRITZ, BLK, RITZ, RITVEC, LDRITV, N, PIPE,
     &                  LDPIPE, RES, LDRES, SAFMIN)
      integer NRITZ, BLK, LDRITV, N, LDPIPE, LDRES
      double precision SAFMIN
      double precision RITZ(NRITZ), RITVEC(LDRITV,NRITZ)
      double precision PIPE(LDPIPE,NRITZ), RES(LDRES,BLK)
      integer IRITZ
      double precision VALUE
      double precision ONE
      parameter (ONE=1.0D0)
      external DGEMV
      intrinsic SIGN, ABS
      do 10 IRITZ = 1,NRITZ
         VALUE = RITZ(IRITZ)
         if (abs(VALUE).ge.SAFMIN) then
            VALUE = sign(ONE/abs(VALUE),VALUE)
            call DGEMV('No transpose', N, BLK, VALUE, RES,LDRES,
     &                 PIPE(1,IRITZ),1, ONE, RITVEC(1,IRITZ),1)
         end if
 10   continue
      return
      end
      subroutine EA16RD(NKRYL, BLK, NCONV, NPURG, T, LDT, RITZ, EIGVEC,
     &                  LDEIGV, LCEIGV, NEIGV, PIPE, LDPIPE, TT, LDTT,
     &                  LCTT)
      integer NKRYL, BLK, NCONV, NPURG, NEIGV, LDT, LDEIGV, LCEIGV, LDTT
      integer LCTT, LDPIPE
      double precision T(LDT,NKRYL), RITZ(NPURG)
      double precision EIGVEC(LDEIGV,LCEIGV), TT(LDTT,LCTT)
      double precision PIPE(LDPIPE,NPURG)
      integer I, J
      double precision ZERO
      parameter (ZERO=0.0D0)
      external DLASET, DCOPY, EA18PD
      if (NPURG.lt.NKRYL) then
         do 40 I = 1,BLK
            do 400 J = 1,NEIGV
               EIGVEC(J,NPURG+I) = EIGVEC(J,NKRYL+I)
 400        continue
 40      continue
      end if
      if (NCONV.gt.0) then
         call DLASET('All', BLK, NCONV, ZERO, ZERO, T(2,1), LDT)
         call DCOPY(NCONV, RITZ,1, T,LDT)
      end if
      if (NCONV.ge.NPURG) return
      call EA18PD(T(1,NCONV+1), LDT, BLK, NPURG-NCONV, RITZ(NCONV+1),
     &            PIPE(1,NCONV+1), LDPIPE,
     &            EIGVEC(1,NCONV+1), LDEIGV, NEIGV, TT, LDTT)
      return
      end
      subroutine EA16SD(V, LDV, LCV, N, M, Z, LDZ, LCZ, K, WORK, LDW,
     &                  LCW)
      integer LDV, LCV, N, M, LDZ, LCZ, K, LDW, LCW
      double precision V(LDV,LCV), Z(LDZ,LCZ), WORK(LDW,*)
      double precision ZERO, ONE
      parameter (ZERO=0.0D0, ONE=1.0D0)
      integer I, STEP
      external DGEMM, DLACPY
      do 10 I = 1,N,LDW
         STEP = min(LDW,N-I+1)
         call DGEMM('No transpose', 'No transpose', STEP, K, M,
     &              ONE, V(I,1),LDV, Z,LDZ, ZERO, WORK, LDW)
         call DLACPY('All', STEP, K, WORK, LDW, V(I,1),LDV)
 10   continue
      return
      end
      subroutine EA16TD(T, LDT, BLK, NSTEP, NCONV, OFFSET,
     &                  ERMOLA, LDERLA, F, LDF, SVD, WORK, LWORK,
     &                  INFO)
      integer BLK, LDT, NSTEP, NCONV, OFFSET, LDERLA, LDF, LWORK
      integer INFO
      double precision T(LDT,NSTEP*BLK+NCONV-OFFSET), ERMOLA(LDERLA,6)
      double precision F(LDF,BLK), SVD(BLK), WORK(LWORK)
      integer ALPHA, BETA
      parameter (ALPHA=4, BETA=5)
      integer ICOL, ISTEP, I
      double precision U(1,1), VT(1,1), RTEMP
      double precision ZERO
      parameter       (ZERO=0.0D0)
      double precision EA18FD
      external         EA18FD
      external EA18LD, DSYEV, DGESVD
      ICOL = -BLK + 1
      do 10 ISTEP=1,NSTEP
         ICOL = ICOL + BLK
         call EA18LD(T(1,ICOL),LDT, BLK, BLK, F, LDF)
         call DSYEV('N', 'L', BLK, F(1,1), LDF, SVD, WORK, LWORK, INFO)
         if (INFO.eq.0) then
            RTEMP = ZERO
            do 100 I = 1,BLK
               RTEMP = max(RTEMP, abs(SVD(I)))
 100        continue
         else
            call EA18LD(T(1,ICOL),LDT, BLK, BLK, F, LDF)
            RTEMP = EA18FD('L', BLK, BLK, F, LDF)
            INFO = 4
         end if
         ermola(ISTEP,ALPHA) =  RTEMP
      if (BLK.EQ.1) then
         SVD(1) = ABS(F(BLK+1,1))
      else
         call DGESVD('N', 'N', BLK, BLK, F(BLK+1,1), LDF,
     &               SVD, U, 1, VT, 1, WORK, LWORK, INFO )
      end if
         if (INFO.eq.0) then
            RTEMP = SVD(1)
         else
            call EA18LD(T(1,ICOL),LDT, BLK, BLK, F, LDF)
            RTEMP = EA18FD('B', BLK, BLK, F(BLK+1,1), LDF)
            INFO = 4
            return
         end if
         ermola(ISTEP,BETA) = RTEMP
 10   continue
      return
      end
      subroutine EA16UD(ISTEP, BLK, ERMOLA, LDERLA, ERMOVE, LDERVE,
     &                  NCONV, SVD, SAFMIN, ERRLOC)
      integer ISTEP, BLK, NCONV, LDERLA, LDERVE
      double precision SAFMIN, ERRLOC
      double precision ERMOLA(LDERLA,5), ERMOVE(LDERVE,4), SVD(BLK)
      double precision ZERO, ONE
      parameter (ZERO=0.0D0, ONE=1.0D0)
      integer ALPHA, BETA, FIRST, SECOND
      parameter (ALPHA=4,BETA=5,FIRST=1,SECOND=2)
      integer TAU, ZNEW
      parameter (TAU=3, ZNEW=4)
      integer K
      double precision INVBET, THIRD, LARLOG
      double precision THIRD1, THIRD2, THIRD3
      INVBET = SVD(BLK)
      if (INVBET.lt.SAFMIN) then
         do 50 K = 1,NCONV
            ermove(K,FIRST) = ermove(K,SECOND)
            ermove(K,SECOND) = ONE
 50      continue
      else
         INVBET = ONE/INVBET
         LARLOG = -log(SAFMIN*4.0D0)
         do 60 K = 1,NCONV
            THIRD = ZERO
            if (ermola(ISTEP,ALPHA)+ermove(K,TAU).gt.ZERO
     &          .and. ermove(K,SECOND).gt.ZERO) then
               THIRD1 = log(INVBET)
     &               + log(ermola(ISTEP,ALPHA) + ermove(K,TAU))
     &               + log(ermove(K,SECOND))
               THIRD1 = min(THIRD1, LARLOG)
               THIRD = THIRD + exp(THIRD1)
            end if
            if (ermola(ISTEP,BETA).gt.ZERO
     &          .and. ermove(K,FIRST).gt.ZERO) then
               THIRD2 = log(INVBET)
     &                + log(ermola(ISTEP,BETA))
     &                + log(ermove(K,FIRST))
               THIRD2 = min(THIRD2, LARLOG)
               THIRD = THIRD + exp(THIRD2)
            end if
            if (ermove(K,ZNEW).gt.ZERO) then
               THIRD3 = log(INVBET)
     &                + log(ermove(K,ZNEW))
               THIRD3 = min(THIRD3, LARLOG)
               THIRD = THIRD + exp(THIRD3)
            end if
            ermove(K,FIRST) = ermove(K,SECOND)
            ermove(K,SECOND) = THIRD
 60      continue
      end if
      ERRLOC = ZERO
      do 70 K = 1,NCONV
         ERRLOC = max(ERRLOC, abs(ERMOVE(K,SECOND)))
 70   continue
      return
      end
      subroutine EA16VD(ISTEP, JSTEP, BLK, NLOCK, ERMOLA, LDERLA,
     &                  TDIA, LDTDIA, TOFF, LDTOFF,
     &                  SVD, WORK, LWORK, SAFMIN, EPSS, ERRLAN,
     &                  INFO)
      integer          ISTEP, JSTEP, BLK, NLOCK, LDERLA, LDTDIA, LDTOFF
      integer          LWORK, INFO
      double precision SAFMIN, EPSS, ERRLAN
      double precision ERMOLA(LDERLA,6), TDIA(LDTDIA,*)
      double precision TOFF(LDTOFF,*), SVD(BLK), WORK(LWORK)
      double precision ZERO, ONE
      parameter       (ZERO=0.0D0, ONE=1.0D0)
      integer    ALPHA, BETA, GAMMA, FIRST, SECOND, THIRD
      parameter (ALPHA=4,BETA=5,GAMMA=6,FIRST=1,SECOND=2,THIRD=3)
      integer          K, KMAX
      double precision U(1,1), VT(1,1), INVBET, RTEMP
      double precision EA18FD
      external         EA18FD
      external DGESVD, DSYEV
      RTEMP = EA18FD('L', BLK, BLK, TDIA, LDTDIA)
      call DSYEV('N', 'L', BLK, TDIA, LDTDIA, SVD, WORK, LWORK, INFO)
      if (INFO.ne.0) then
         INFO = 4
      else
         RTEMP = ZERO
         do 10 K = 1,BLK
            RTEMP = max(RTEMP, abs(SVD(K)))
 10      continue
      end if
      ermola(ISTEP,ALPHA) =  RTEMP
      RTEMP = EA18FD('B', BLK, BLK, TOFF, LDTOFF)
      if (BLK.EQ.1) then
         SVD(1) = ABS(TOFF(1,1))
      else
         call DGESVD('N', 'N', BLK, BLK, TOFF, LDTOFF,
     &                SVD, U, 1, VT, 1, WORK, LWORK, INFO )
      end if
      if (INFO.ne.0) then
         INFO = 4
         INVBET = ZERO
         ermola(ISTEP,BETA) = RTEMP
      else
         ermola(ISTEP,BETA) =  SVD(1)
         INVBET = SVD(BLK)
      end if
      do 30 K = NLOCK+1,ISTEP-2
         ermola(K,FIRST) = ermola(K,SECOND)
         ermola(K,SECOND) = ermola(K,THIRD)
 30   continue
      if (ISTEP.gt.1) ermola(ISTEP-1,SECOND) = ermola(ISTEP-1,THIRD)
      do 40 K = JSTEP,ISTEP
         ermola(K,THIRD) = EPSS
 40   continue
      if (JSTEP.gt.1) then
         if (INVBET.lt.SAFMIN) then
            KMAX = min(JSTEP-1,ISTEP)
            do 50 K = 1,KMAX
               ermola(K,THIRD) = ONE
 50         continue
         else if (ISTEP.gt.2) then
            INVBET = ONE/INVBET
            KMAX = min(JSTEP-1,ISTEP-2)
            do 60 K = NLOCK+1,KMAX
               if (K.eq.NLOCK+1) then
                  ermola(K,THIRD) = ZERO + ermola(K,GAMMA)
               else
                  ermola(K,THIRD) = ermola(K-1,BETA)*ermola(K-1,SECOND)
     &                            + ermola(K,GAMMA)
               end if
               ermola(K,THIRD) = INVBET * (
     &             + ermola(K,THIRD)
     &             + ermola(K,BETA)*ermola(K+1,SECOND)
     &             + ermola(ISTEP-1,BETA)*ermola(K,FIRST)
     &             + (ermola(ISTEP,ALPHA) + ermola(K,ALPHA))
     &                       *ermola(K,SECOND)
     &             )
 60         continue
            if (JSTEP.gt.ISTEP-1) then
               ermola(ISTEP-1,THIRD) = INVBET * (
     &                ( 2*ermola(ISTEP-1,BETA) + ermola(ISTEP,ALPHA)
     &                  + ermola(ISTEP-1,ALPHA)
     &                ) * EPSS
     &                + ermola(ISTEP-2,BETA)*ermola(ISTEP-2,SECOND)
     &                + ermola(ISTEP-1,GAMMA)
     &            )
            end if
         else
            INVBET = ONE/INVBET
            if (JSTEP.gt.ISTEP-1) then
               ermola(ISTEP-1,THIRD) = INVBET * (
     &                (  2*ermola(ISTEP-1,BETA)
     &                   + ermola(ISTEP,ALPHA) + ermola(ISTEP-1,ALPHA)
     &                ) * EPSS
     &                + ermola(ISTEP-1,GAMMA)
     &            )
            end if
         end if
      end if
      ERRLAN = ZERO
      do 80 K = NLOCK+1,ISTEP-1
         ERRLAN = max(ERRLAN, abs(ermola(K,THIRD)))
 80   continue
      return
      end
      subroutine EA16WD(RITZ, NCONV, RESID, ERRMAX, GAMMA,
     &                  ERMOVE, LDERVE)
      integer NCONV, LDERVE
      double precision RITZ(NCONV), ERMOVE(LDERVE,4)
      double precision RESID(NCONV)
      double precision GAMMA, ERRMAX
      integer TAU, ZNEW
      parameter (TAU=3, ZNEW=4)
      integer I
      double precision VALUE
      do 10 I = 1,NCONV
         ermove(I,TAU) = abs(RITZ(I))
 10   continue
      do 20 I = 1,NCONV
         ermove(I,ZNEW) = RESID(I) + GAMMA
 20   continue
      do 30 I = 1,NCONV
         VALUE = max(ERRMAX, RESID(I))
         ermove(I,1) = VALUE
         ermove(I,2) = VALUE
 30   continue
      return
      end
      subroutine EA16XD(NSTEP, LSTEP, ERMO, LDER, ERRLEV)
      integer NSTEP, LSTEP, LDER
      double precision ERRLEV
      double precision ERMO(LDER,6)
      integer GAMMA
      parameter (GAMMA=6)
      integer I
      double precision ZERO
      parameter       (ZERO=0.0D0)
      do 10 I = 1,LSTEP
         ERMO(I,GAMMA) = ERRLEV
 10   continue
      do 20 I = LSTEP+1,NSTEP
         ERMO(I,GAMMA) = ZERO
 20   continue
      return
      end
      subroutine EA16YD(NSTEP, KBLOCK, ERMO, LDER, ERRLEV)
      integer NSTEP, KBLOCK, LDER
      double precision ERRLEV
      double precision ERMO(LDER,KBLOCK)
      integer I, J
      double precision ZERO
      parameter       (ZERO=0.0D0)
      do 10 I = 1,KBLOCK
         do 100 J = 1,NSTEP
            ERMO(J,I) = ERRLEV
 100     continue
         do 120 J = NSTEP+1,LDER
            ERMO(J,I) = ZERO
 120     continue
 10   continue
      return
      end
      subroutine EA16ZD(IERR, WARN)
      integer IERR, WARN
      if (mod(IERR/WARN,2).eq.0) IERR = IERR + WARN
      return
      end
      subroutine EA17AD(IDO, WHICH, MODE, HARMON, POLE, NCOLEF, NCORIG,
     &                  NWACO, NCONV, NKRYL, RITZ, RANGE, PERM, RITZSO,
     &                  SAFMIN)
      integer IDO, WHICH, MODE, HARMON, NCORIG, NCOLEF, NWACO, NCONV
      integer NKRYL
      double precision RANGE(2), POLE
      double precision SAFMIN
      integer PERM(NKRYL)
      double precision RITZ(NKRYL)
      double precision RITZSO(NKRYL)
      integer I, J, ITEMP
      integer NVAL(2)
      double precision TARGET, VALUE
      external EA17PD, KB08AD, EA18MD
      if (NKRYL.le.0) then
         IDO = 100
         return
      end if
      if (IDO.eq.0) then
         do 100 I = 1,NKRYL
            PERM(I) = I
            RITZSO(I) = RITZ(I)
 100     continue
         call EA17PD(MODE, HARMON, NKRYL, RITZSO, POLE, RANGE, SAFMIN)
         TARGET = RANGE(1)
         if (abs(WHICH).eq.5) then
            TARGET = ( RANGE(1) + RANGE(2) ) * 0.5D0
         end if
      end if
      if (WHICH.eq.2 .or. WHICH.eq.-3) then
         do 210 I = 1,NKRYL
            RITZSO(I) = -RITZSO(I)
 210     continue
      else if (WHICH.eq.-2 .or. WHICH.eq.3) then
         do 220 I = 1,NKRYL
            RITZSO(I) = RITZSO(I)
 220     continue
      else if (WHICH.eq.-1 .or. WHICH.eq.5) then
         do 230 I = 1,NKRYL
            RITZSO(I) = -abs(RITZSO(I)-TARGET)
 230     continue
      else if (WHICH.eq.1) then
         do 240 I = 1,NKRYL
            RITZSO(I) = abs(RITZSO(I)-TARGET)
 240     continue
      else if (WHICH.eq.-4) then
         VALUE = RANGE(1)
         do 250 I = 1,NKRYL
            VALUE = max(VALUE,RITZSO(I))
 250     continue
         do 255 I = 1,NKRYL
            if (RITZSO(I).ge.RANGE(1)) then
               RITZSO(I) = VALUE - RITZSO(I)
            else
               RITZSO(I) = -abs(RITZSO(I)-TARGET)
            end if
 255     continue
      else if (WHICH.eq.4) then
         VALUE = RANGE(1)
         do 260 I = 1,NKRYL
            VALUE = min(VALUE,RITZSO(I))
 260     continue
         do 265 I = 1,NKRYL
            if (RITZSO(I).le.RANGE(1)) then
               RITZSO(I) = RITZSO(I) - VALUE
            else
               RITZSO(I) = -abs(RITZSO(I)-TARGET)
            end if
 265     continue
      else if (WHICH.eq.10) then
         if (IDO.eq.0) then
            IDO = 6
            return
         else
            do 270 I = 1,NKRYL
               RITZSO(I) = PERM(I)
 270        continue
         end if
      end if
      if (NKRYL.gt.0) call KB08AD(RITZSO, NKRYL, PERM)
      if (WHICH.eq.5) then
         do 72 I = 1,NKRYL
            RITZSO(I) = RITZ(PERM(I))
 72      continue
         call EA17PD(MODE, HARMON, NKRYL, RITZSO, POLE, RANGE, SAFMIN)
         NVAL(1) = NKRYL
         do 62 I = 1,NKRYL
            VALUE = RITZSO(I)
            if (VALUE.lt.RANGE(1) .or.
     &          VALUE.gt.RANGE(2)) then
               NVAL(1) = I-1
               go to 63
            end if
 62      continue
 63      continue
         do 71 I = 1,NVAL(1)
            RITZSO(I) = -RITZSO(I)
 71      continue
         call EA18MD(RITZSO, NVAL(1), PERM)
      else if (abs(WHICH).eq.3) then
         do 68 I = 2,NKRYL,2
            ITEMP = PERM(NKRYL)
            do 680 J = NKRYL-1,I,-1
               PERM(J+1) = PERM(J)
 680        continue
            PERM(I) = ITEMP
 68      continue
      end if
      NWACO = 0
      do 80 I = 1,NKRYL
         if (PERM(I).le.NCONV) then
            NWACO = NWACO + 1
         else
            go to 81
         end if
 80   continue
 81   continue
      if (WHICH.eq.3 .or. WHICH.eq.-3) then
         NCOLEF = 0
         do 30 I = 1,NKRYL,2
            if (PERM(I).le.NCONV) then
               NCOLEF = NCOLEF + 1
            else
               go to 31
            end if
 30      continue
 31      continue
         NCORIG = 0
         do 35 I = 2,NKRYL,2
            if (PERM(I).le.NCONV) then
               NCORIG = NCORIG + 1
            else
               go to 36
            end if
 35      continue
 36      continue
         if (WHICH.eq.3) then
            I = NCORIG
            NCORIG = NCOLEF
            NCOLEF = I
         end if
      end if
      IDO = 100
      return
      end
      subroutine EA17BD(NCONV, NCOOLD, NCOKEE,  NKRYL, PERM)
      integer NCONV, NCOOLD, NCOKEE, NKRYL
      integer PERM(NKRYL)
      integer I, J, JCONV, KCONV, LCONV, ITEMP
      JCONV = 0
      KCONV = NCOOLD
      do 10 I = 1,NKRYL
         if (PERM(I).le.NCOOLD) then
            JCONV = JCONV + 1
         else if (PERM(I).le.NCONV) then
            JCONV = JCONV + 1
            KCONV = KCONV + 1
         end if
         if (JCONV.ge.NCOKEE) go to 11
 10   continue
 11   continue
      LCONV = KCONV
      JCONV = 0
      do 20 I = 1,NKRYL
         if (PERM(I).le.NCOOLD) then
            ITEMP = PERM(I)
            JCONV = JCONV + 1
            do 200 J = I-1,JCONV,-1
               PERM(J+1) = PERM(J)
 200        continue
            PERM(JCONV) = ITEMP
         else
            PERM(I) = PERM(I) - NCOOLD
         end if
 20   continue
      JCONV = NCOOLD
      do 30 I = NCOOLD+1,NKRYL
         if (PERM(I).le.NCONV-NCOOLD) then
            if (JCONV.ge.LCONV) go to 31
            JCONV = JCONV+1
            ITEMP = PERM(I)
            do 300 J = I-1,JCONV,-1
               PERM(J+1) = PERM(J)
 300        continue
            PERM(JCONV) = ITEMP
         end if
 30   continue
 31   continue
      NCONV = LCONV
      return
      end
      subroutine EA17CD(WHAT, NCOOLD, NKRYL, RITZ, RESID,
     &                  V, LDV, NV, PIPE, LDPIPE, PERM, WORK)
      integer WHAT, NCOOLD
      integer NKRYL, LDV, NV, LDPIPE
      integer PERM(NKRYL)
      double precision RITZ(NKRYL), RESID(NKRYL), V(LDV,NKRYL-NCOOLD)
      double precision PIPE(LDPIPE,NKRYL-NCOOLD), WORK(NKRYL)
      external EA18OD
      if (NKRYL.le.0) return
      if (WHAT.gt.0) then
         call EA18OD(NKRYL-NCOOLD, RITZ(NCOOLD+1),1,1, PERM(NCOOLD+1),
     &               WORK)
      end if
      if (WHAT.gt.1) then
         call EA18OD(NKRYL-NCOOLD, RESID,1,1, PERM(NCOOLD+1), WORK)
         call EA18OD(NKRYL-NCOOLD, V,LDV,NV, PERM(NCOOLD+1), WORK)
      end if
      if (WHAT.gt.2) then
         call EA18OD(NKRYL-NCOOLD, PIPE,LDPIPE,LDPIPE,
     &               PERM(NCOOLD+1), WORK)
      end if
      return
      end
      subroutine EA17DD(IDO, MODE, HARMON, PURIF, N, WHICH, SIGMA,
     &                  NEINEG, EIGEN, NRITZ, NCONV, NWANT, TRUST,
     &                  RANGE, RSAVE, ISAVE, OK, MAXCND, EPS, SAFMIN,
     &                  NORESC, LU)
      logical OK, NORESC
      integer IDO, MODE, HARMON, PURIF, N, WHICH, NWANT, NEINEG
      integer NCONV, LU, NRITZ
      integer ISAVE(6)
      double precision SIGMA, MAXCND, EPS, SAFMIN
      double precision EIGEN(NRITZ)
      double precision RANGE(2), TRUST(4), RSAVE(8)
      logical ALLOW2, PURIF5
      integer INERT1, INERT2, NBETST, NEITST, NEINE2, I
      double precision POLE1, POLE2, LBOUND, RBOUND
      double precision SIGOLD, SIGTST, ABSDIS, RELDIS, LARGE
      double precision RANLOC(2)
      integer LPOLE1, LPOLE2, LLBOUN, LRBOUN, LSIGOL, LSIGTS, LABSDI
      parameter (LPOLE1=1, LPOLE2=2, LLBOUN=3, LRBOUN=4, LSIGOL=5,
     &           LSIGTS=6, LABSDI=7)
      integer LINER1, LINER2, LNBETS, LNEITS, IALLOW
      parameter (LINER1=1, LINER2=2, LNBETS=4, LNEITS=5,
     &           IALLOW=6)
      double precision ONE
      parameter (ONE=1.0D0)
      integer EA18ZD
      double precision DLAMCH
      external EA18ZD, EA18RD, DLAMCH
      external EA17PD
      intrinsic MIN, MAX
      RELDIS = (ONE / MAXCND) * (ONE+sqrt(EPS))
      if (IDO.eq.0) then
         call EA17PD(MODE, HARMON, NRITZ, EIGEN, SIGMA, RANGE, SAFMIN)
      end if
      ALLOW2 = isave(IALLOW).eq.1
      PURIF5 = MODE.eq.5 .and. PURIF.gt.0
      ABSDIS = rsave(LABSDI)
      if (WHICH.eq.-1 .or. WHICH.eq.4 .or. WHICH.eq.-4 .or.
     &    WHICH.eq.5) then
         RANLOC(1) = RANGE(1)
         if (WHICH.eq.5) RANLOC(2) = RANGE(2)
      end if
      if (.not.ALLOW2) then
         if (WHICH.eq.2) then
            RANLOC(1) = ONE/SAFMIN
            do 10 I = 1,NRITZ
               RANLOC(1) = min(RANLOC(1),
     &                  EIGEN(I)-RELDIS*abs(SIGMA-EIGEN(I)))
 10         continue
         else if (WHICH.eq.-2) then
            RANLOC(1) = -(ONE/SAFMIN)
            do 20 I = 1,NRITZ
               RANLOC(1) = max(RANLOC(1),
     &                  EIGEN(I)+RELDIS*abs(SIGMA-EIGEN(I)))
 20         continue
         end if
      end if
      NEINE2 = EA18ZD(NEINEG, MODE, SIGMA, N)
      LARGE = DLAMCH('Overflow')
      INERT1 = isave(LINER1)
      INERT2 = isave(LINER2)
      NEITST = isave(LNEITS)
      NBETST = isave(LNBETS)
      POLE1  = rsave(LPOLE1)
      POLE2  = rsave(LPOLE2)
      LBOUND = rsave(LLBOUN)
      RBOUND = rsave(LRBOUN)
      SIGOLD = rsave(LSIGOL)
      SIGTST = rsave(LSIGTS)
      ABSDIS = rsave(LABSDI)
      call EA18RD(IDO, WHICH, ALLOW2, PURIF5, N, RANLOC, OK, INERT1,
     &            INERT2, POLE1, POLE2, LBOUND, RBOUND, SIGOLD, SIGMA,
     &            NEINE2, EIGEN, NRITZ, NCONV, NWANT, TRUST, ABSDIS,
     &            RELDIS, SIGTST, NEITST, NBETST, LARGE, NORESC, LU)
      isave(LINER1) = INERT1
      isave(LINER2) = INERT2
      isave(LNEITS) = NEITST
      isave(LNBETS) = NBETST
      rsave(LPOLE1) = POLE1
      rsave(LPOLE2) = POLE2
      rsave(LLBOUN) = LBOUND
      rsave(LRBOUN) = RBOUND
      rsave(LSIGOL) = SIGOLD
      rsave(LSIGTS) = SIGTST
      return
      end
      subroutine EA17ED(WHICH, MODE, IRORIG,  IRFINA, LU)
      integer WHICH, MODE, IRORIG, IRFINA, LU
      integer ABSWHI
      IRFINA = IRORIG
      if (IRFINA.eq.4) then
         ABSWHI = abs(WHICH)
         if ((MODE.eq.1 .or. MODE.eq.3) .and. (
     &                       ABSWHI.eq.3 .or. ABSWHI.eq.2
     &                       .or. WHICH.eq.1 )
     &      ) then
         else
            if (LU.ge.0) write(LU,9001)
            IRFINA = 1
         end if
      else if (IRFINA.ne.1 .and. IRFINA.ne.2 .and. IRFINA.ne.10) then
         IRFINA = 1
         if (LU.ge.0) write(LU,9003)
      end if
 9001 format('WARNING : Chebyshev shifts cannot be selected for this',/,
     &       '          combination of MODE and WHICH.',/,
     &       '          The code uses purging.',/)
 9003 format('WARNING : the value of ICNTL(9) is out of bounds',/,
     &       '          The code uses purging.',/)
      return
      end
      subroutine EA17FD(IDO, NRESTA, NKRYL, NKEEP, BLK, WHICH,
     &                  RANGE, INTERV, SHIFT, MSHIFT, RITZ,
     &                  IRSTRT, NSHIFT, MLEJA, NLEJA, NDONE, ILEJA,
     &                  RLEJA, SAFMIN)
      integer IDO, NKRYL, NKEEP, BLK, WHICH, MSHIFT, IRSTRT
      integer NRESTA, NSHIFT, MLEJA, NLEJA, NDONE
      integer ILEJA(2*MLEJA)
      double precision SAFMIN
      double precision RANGE(2), INTERV(4), SHIFT(MSHIFT), RITZ(NKRYL)
      double precision RLEJA(3*MLEJA)
      logical FIRST, ODDIT, LEFT, RIGHT
      integer I, STEP, LKEEP
      double precision LRANGE(2)
      double precision INTER1, INTER2, INTER3, INTER4
      double precision ZERO
      parameter       (ZERO=0.0D0)
      external EA17GD, EA17HD, EA17KD
      FIRST = NRESTA.eq.0
      ODDIT = mod(NRESTA,2).eq.0
      IDO = 100
      NSHIFT = (NKRYL-NKEEP)/BLK
      if (NSHIFT.le.0) return
      if (IRSTRT.eq.2) then
         STEP = 1
         do 10 I = 1,NSHIFT
            SHIFT(I) = RITZ(NKRYL+STEP-I*STEP)
 10      continue
         return
      else if (IRSTRT.eq.10) then
         IDO = 7
         return
      end if
      LRANGE(1) = RANGE(1)
      LRANGE(2) = RANGE(2)
      INTER1 = INTERV(1)
      INTER2 = INTERV(2)
      LKEEP = max(1,NKEEP+(NKRYL-NKEEP)/3)
      if (WHICH.eq.1 .or. WHICH.eq.3 .or. WHICH.eq.-3
     &    .or. WHICH.eq.-2 .or. WHICH.eq.2) then
         LEFT = .true.
         RIGHT = .true.
         if (WHICH.eq.1 .and..not.FIRST) then
            LEFT = INTERV(1).lt.RANGE(1)
            RIGHT = INTERV(1).gt.RANGE(1)
         else if (WHICH.eq.2) then
            RIGHT = .false.
         elseif (WHICH.eq.-2) then
            LEFT = .false.
         end if
         call EA17GD(FIRST, IRSTRT, LEFT, RIGHT, INTER1,
     &               INTER2, RITZ, NKRYL, LKEEP, NSHIFT, SHIFT,
     &               MLEJA, NLEJA, NDONE, ILEJA, RLEJA)
      else if (WHICH.eq.-1) then
         INTER3 = INTERV(3)
         INTER4 = INTERV(4)
         call EA17HD(FIRST, ODDIT, IRSTRT, LRANGE(1), INTER1,
     &               INTER2, INTER3, INTER4, RITZ,
     &               NKRYL, NKEEP, NSHIFT, SHIFT, MLEJA, NLEJA, NDONE,
     &               ILEJA, RLEJA, SAFMIN)
         INTERV(3) = INTER3
         INTERV(4) = INTER4
      else if (WHICH.eq.5) then
         INTER3 = INTERV(3)
         INTER4 = INTERV(4)
         call EA17KD(FIRST, ODDIT, IRSTRT, ZERO, RANGE(1),
     &               RANGE(2), INTER1, INTER2, INTER3,
     &               INTER4, RITZ, NKRYL, NKEEP, NSHIFT,
     &               SHIFT, MLEJA, NLEJA, NDONE, ILEJA, RLEJA,
     &               SAFMIN)
         INTERV(3) = INTER3
         INTERV(4) = INTER4
      else if (WHICH.eq.4) then
         INTER3 = INTERV(3)
         INTER4 = INTERV(4)
         call EA17KD(FIRST, ODDIT, IRSTRT, ZERO, RANGE(1),
     &               RANGE(1), INTER1, INTER2, INTER3,
     &               INTER4, RITZ, NKRYL, NKEEP, NSHIFT,
     &               SHIFT, MLEJA, NLEJA, NDONE, ILEJA, RLEJA,
     &               SAFMIN)
         INTERV(3) = INTER3
         INTERV(4) = INTER4
      else if (WHICH.eq.-4) then
         INTER3 = INTERV(3)
         INTER4 = INTERV(4)
         call EA17KD(FIRST, ODDIT, IRSTRT, ZERO, RANGE(1),
     &               RANGE(1), INTER1, INTER2, INTER3,
     &               INTER4, RITZ, NKRYL, NKEEP, NSHIFT,
     &               SHIFT, MLEJA, NLEJA, NDONE, ILEJA, RLEJA,
     &               SAFMIN)
         INTERV(3) = INTER3
         INTERV(4) = INTER4
      end if
      INTERV(1) = INTER1
      INTERV(2) = INTER2
      return
      end
      subroutine EA17GD(FIRST, IRST, LEFT, RIGHT, MINVAL, MAXVAL,
     &                 RITZ, NRITZ, NKEEP, NSHIFT, SHIFT, MLEJA, NLEJA,
     &                 NDONE, IW, RW)
      logical FIRST, LEFT, RIGHT
      integer IRST, NKEEP, NRITZ, NSHIFT, MLEJA, NLEJA, NDONE
      double precision MINVAL, MAXVAL
      integer IW(2*MLEJA)
      double precision RITZ(NRITZ), SHIFT(NSHIFT), RW(3*MLEJA)
      integer I
      double precision MEAN, LENGT
      double precision HALFPI, AHALF
      parameter (HALFPI = 1.5707963267949D0, AHALF=0.5D0)
      if (FIRST) then
         NLEJA = 0
         MINVAL = RITZ(NKEEP+1)
         MAXVAL = MINVAL
      end if
      do 30 I = 1,NKEEP
         if (RITZ(I).le.MAXVAL .and. RITZ(I).ge.MINVAL) then
            if (LEFT) then
               MINVAL = RITZ(NKEEP+1)
               NLEJA = 0
            end if
            if (RIGHT) then
               MAXVAL = RITZ(NKEEP+1)
               NLEJA = 0
            end if
            go to 31
         end if
 30   continue
 31   continue
      do 10 I = NKEEP+1,NRITZ
         if (RITZ(I).lt.MINVAL) then
            MINVAL = RITZ(I)
         end if
         if (RITZ(I).gt.MAXVAL) then
            MAXVAL = RITZ(I)
         end if
 10   continue
      if (IRST.eq.4) then
         MEAN = (MAXVAL+MINVAL)*AHALF
         LENGT = (MAXVAL-MINVAL)*AHALF
         do 20 I = 1,NSHIFT
            SHIFT(I) = MEAN + LENGT * COS((-1+2*I)*HALFPI/NSHIFT)
 20      continue
      end if
      return
      end
      subroutine EA17HD(FIRST, ODDIT, IRST, MIDDLE, A, B, C, D,
     &                  RITZ, NRITZ, NKEEP, NSHIFT, SHIFT, MLEJA,
     &                  NLEJA, NDONE, IW, RW, SAFMIN)
      logical FIRST
      logical ODDIT
      integer IRST, NRITZ, NSHIFT, NKEEP, MLEJA, NLEJA, NDONE
      integer IW(2*MLEJA)
      double precision MIDDLE, A, B, C, D, SAFMIN
      double precision RITZ(NRITZ), SHIFT(NSHIFT), RW(3*MLEJA)
      integer I
      double precision VALUE
      double precision LENGT
      if (FIRST) then
         A = 0.0D0
         B = 0.0D0
         C = 0.0D0
         D = 0.0D0
      end if
      if (A.eq.B) then
         A = MIDDLE
         do 10 I = NKEEP+1,NRITZ
            A = min(A, RITZ(I))
 10      continue
         B = A
      end if
      if (C.eq.D) then
         D = MIDDLE
         do 15 I = NKEEP+1,NRITZ
            D = max(D, RITZ(I))
 15      continue
         C = D
      end if
      do 20 I = 1,NKEEP
         if (RITZ(I).le.B) B = A
         if (RITZ(I).ge.C) C = D
 20   continue
      do 30 I = NKEEP+1,NRITZ
         VALUE = RITZ(I)
         if (VALUE.ge.MIDDLE) then
            C = min(C,VALUE)
         else
            B = max(B,VALUE)
         end if
         A = min(A,VALUE)
         D = max(D,VALUE)
 30   continue
      LENGT = max(C-MIDDLE, MIDDLE-B)
      C = MIDDLE + LENGT
      B = MIDDLE - LENGT
      A = min(A,B)
      D = max(C,D)
      return
      end
      subroutine EA17JD(N, V, VALUE, SAFMIN, INFO)
      integer N, INFO
      double precision VALUE, SAFMIN
      double precision V(N)
      double precision SCALAR, MAXV
      double precision ONE
      parameter (ONE=1.0D0)
      integer IDAMAX
      external IDAMAX
      external DRSCL, DSCAL
      INFO = 0
      MAXV = abs(V(IDAMAX(N, V, 1)))
      if (VALUE.lt.MAXV*SAFMIN) then
         INFO = -1
         return
      else if (VALUE.lt.SAFMIN) then
         call DRSCL(N, VALUE, V,1)
      else
         SCALAR = ONE / VALUE
         call DSCAL(N, SCALAR, V,1)
      end if
      return
      end
      subroutine EA17KD(FIRST, ODDIT, IRST, TOL, R1, R2,
     &                  A, B, C, D, RITZ, NRITZ, NKEEP, NSHIFT, SHIFT,
     &                  MLEJA, NLEJA, NDONE, IW, RW, SAFMIN)
      logical FIRST, ODDIT
      integer IRST, NRITZ, NKEEP, NSHIFT, MLEJA, NLEJA, NDONE
      integer IW(2*MLEJA)
      double precision TOL, R1, R2, A, B, C, D, SAFMIN
      double precision RITZ(NRITZ), SHIFT(NSHIFT), RW(3*MLEJA)
      integer I
      double precision VALUE, c cc, B1, BOLD, COLD, AOLD, DOLD
      logical ALLIN
      B1 = R1-TOL
      c cc = R2+TOL
      ALLIN = .true.
      if (FIRST .or. A.eq.B) then
         A = B1
         do 10 I = NKEEP+1,NRITZ
            VALUE = RITZ(I)
            A = min(A,VALUE)
 10      continue
      else
         AOLD = A
         BOLD = B
      end if
      if (.not. FIRST .and. A.eq.B) then
         AOLD = B
         BOLD = B
      end if
      if (FIRST .or. C.eq.D) then
         D = c cc
         do 15 I = NKEEP+1,NRITZ
            VALUE = RITZ(I)
            D = max(D,VALUE)
 15      continue
      else
         COLD = C
         DOLD = D
      end if
      if (.not. FIRST .and. C.eq.D) then
         COLD = C
         DOLD = C
      end if
      B = A
      C = D
      do 20 I = NKEEP+1,NRITZ
         VALUE = RITZ(I)
         if (VALUE.ge.c cc) then
            C = min(C,VALUE)
         elseif (VALUE.le.B1) then
            B = max(B,VALUE)
         else
            ALLIN = .false.
         end if
         A = min(A,VALUE)
         D = max(D,VALUE)
 20   continue
      c cc = R2 + 0.01D0*(D-R2)
      B1 = R1 - 0.01D0*(R1-A)
      if (.not.ALLIN) then
         C = c cc
         B = B1
      else
         C = max(C,c cc)
         B = min(B,B1)
      end if
      if (.not.FIRST) then
         if (COLD.lt.DOLD) C = min(COLD, C)
         if (AOLD.lt.BOLD) B = max(BOLD, B)
      end if
      return
      end
      subroutine EA17LD(T, LDT, NKRYL, BLK, SHIFT, NSHIFT, Z, LDZ, NZ,
     &                  LCZ, TT, LDTT, LCTT, KEEP, HV, LDHV, LCHV, TAU,
     &                  LTAU, WORK, LWORK)
      logical KEEP
      integer LDT, BLK, NKRYL, NSHIFT, NZ, LDZ, LCZ, LDTT
      integer LCTT, LDHV, LCHV, LTAU, LWORK
      double precision SHIFT(NSHIFT)
      double precision T(LDT,NKRYL), TT(LDTT,LCTT), HV(LDHV,LCHV)
      double precision TAU(LTAU), Z(LDZ,LCZ), WORK(LWORK)
      integer ISTEP, IKRYL, IKEEP
      external EA18DD
      IKRYL = NKRYL
      IKEEP = 1
      do 10 ISTEP=1,NSHIFT
         call EA18DD(T, LDT, IKRYL, BLK, SHIFT(ISTEP), Z, LDZ, NZ, LCZ,
     &               TT, LDTT, LCTT, KEEP, HV(1,IKEEP), LDHV,
     &               LCHV-IKEEP+1, TAU(IKEEP), LTAU-IKEEP+1, WORK,
     &               LWORK)
         if (KEEP) IKEEP = IKEEP + NKRYL
         IKRYL = IKRYL - BLK
 10   continue
      return
      end
      subroutine EA17ND(IDO, FIRST, USER, ISINF, MODE, HARMON, PURIF,
     &                  WHICH, SIGMA, RANGE, MAXCND, GROW,
     &                  CVGED, RITZ, RESID, EIGEN, NCONV, NWACO, NRITZ,
     &                  NWANT, NKEEP, IERR, NEINEG, TRUST, IWORK, ISAVE,
     &                  RSAVE, EXPRST, SAFMIN, EPS, N, LUWARN, LU)
      logical FIRST, USER, ISINF, CVGED, EXPRST
      integer IDO, MODE, HARMON, PURIF, WHICH, NCONV
      integer NWACO, NRITZ, NWANT, NKEEP, N, NEINEG, LU, LUWARN
      integer IERR
      integer ISAVE(6), IWORK(*)
      double precision SIGMA, MAXCND, EPS, SAFMIN, GROW
      double precision RANGE(2), RITZ(*), RESID(*), EIGEN(*)
      double precision TRUST(4), RSAVE(8)
      double precision ZERO,ONE
      parameter (ZERO=0.0D0, ONE=1.0D0)
      integer IPOLE1, IPOLE2, ILBOUN, IRBOUN, ISIGOL, SIGTST, IABSDI
      parameter (IPOLE1=1, IPOLE2=2, ILBOUN=3, IRBOUN=4, ISIGOL=5)
      parameter (SIGTST=6, IABSDI=7)
      integer LINER1, LINER2, LTRIAL, NEITST, NBETST, IALLOW
      parameter (LINER1=1, LINER2=2, LTRIAL=3, NEITST=4, NBETST=5,
     &           IALLOW=6)
      logical ALLOW2
      double precision ALPHA, BETA, GAMMA, DELTA
      double precision LARGE, MINVAL, MAXVAL
      integer IERR2, NWAC2, I
      integer INERT1, INERT2, TRIAL
      double precision POLE1, POLE2, LBOUND, RBOUND, SIGOLD, ABSDIS
      double precision MAXDST
      integer EA18ZD
      double precision DLAMCH
      external DLAMCH, EA18ZD
      external DCOPY, EA17WD, EA18SD, EA18TD
      external EA16ZD, EA17OD, EA17PD
      intrinsic max, min, sqrt
      CVGED = .false.
      LARGE = ONE/SAFMIN
      if (IDO.ge.0 .and. IDO.le.2) isave(LTRIAL) = 0
 1    continue
      if (IDO.lt.3 .and. IDO.ge.0) then
         if (FIRST) rsave(IABSDI) = SAFMIN
         if (NRITZ.gt.1) then
            MAXVAL = EIGEN(1)
            MINVAL = EIGEN(1)
            do 10 I = 1,NRITZ
               MINVAL = min(MINVAL, EIGEN(I))
               MAXVAL = max(MAXVAL, EIGEN(I))
 10         continue
            rsave(IABSDI) = max(rsave(IABSDI),
     &                          (MAXVAL-MINVAL) * sqrt(EPS))
         end if
         if (IDO.eq.0) then
            isave(LINER1) = -N-1
            isave(LINER2) = -N-1
            rsave(IPOLE1) = ZERO
            rsave(IPOLE2) = ZERO
            rsave(ISIGOL) = ZERO
            rsave(SIGTST) = DLAMCH('Overflow')
            isave(NEITST) = -1
            isave(NBETST) = -1
            rsave(IRBOUN) =  DLAMCH('Overflow')
            rsave(ILBOUN) = -DLAMCH('Overflow')
            TRUST(1) = -DLAMCH('Overflow')
            TRUST(4) =  DLAMCH('Overflow')
            TRUST(2) =  TRUST(4)
            TRUST(3) =  TRUST(1)
            isave(IALLOW) = 1
            if (MODE.eq.5 .or. PURIF.gt.0) isave(IALLOW) = 0
         end if
         if (NRITZ.eq.0) then
            if (WHICH.eq.-1 .or. WHICH.eq.4 .or. WHICH.eq.-4
     &         .or. WHICH.eq.5) SIGMA = RANGE(1)
            rsave(IPOLE1) = SIGMA
            rsave(IPOLE2) = SIGMA
            rsave(ISIGOL) = SIGMA
            isave(IALLOW) = 1
            if (MODE.eq.5 .or. PURIF.gt.0) isave(IALLOW) = 0
            if (MODE.eq.5 .and. abs(SIGMA).lt.SAFMIN) then
               IERR = -15
               IDO = 100
               return
            end if
            IDO = 4
            return
         end if
         call DCOPY(NRITZ, RITZ,1, EIGEN,1)
         call EA17PD(MODE, HARMON, NRITZ, EIGEN, SIGMA, RANGE, SAFMIN)
         if (MODE.eq.2 .or. MODE.eq.4 .or. MODE.eq.5)
     &      rsave(ISIGOL) = SIGMA
         POLE1  = rsave(IPOLE1)
         POLE2  = rsave(IPOLE2)
         SIGOLD = rsave(ISIGOL)
         RBOUND = rsave(IRBOUN)
         LBOUND = rsave(ILBOUN)
         ABSDIS = rsave(IABSDI)
         INERT1 = isave(LINER1)
         INERT2 = isave(LINER2)
         TRIAL  = isave(LTRIAL)
         ALLOW2 = isave(IALLOW).eq.1
         if (WHICH.eq.2 .or. WHICH.eq.-2) then
            call EA18SD(FIRST, MODE, ALLOW2, WHICH, RANGE, SIGMA,
     &                  POLE1, POLE2, ISINF, SIGOLD, INERT1, INERT2,
     &                  LBOUND, RBOUND, RITZ, NRITZ, NCONV, NWACO,
     &                  NKEEP-NCONV, NWANT, RESID, EIGEN, IWORK,
     &                  N, TRUST, TRIAL, MAXCND, ABSDIS, EPS,
     &                  SAFMIN, IERR, LU)
         else if (WHICH.eq.5 .or. WHICH.eq.-4 .or. WHICH.eq.4) then
            NWAC2 = NWACO
            call EA18TD(MODE, WHICH, RANGE, SIGMA, POLE1,
     &                  POLE2, LBOUND, RBOUND, ISINF, SIGOLD, INERT1,
     &                  INERT2, CVGED, RITZ, NRITZ, NCONV, NWAC2,
     &                  NWANT, NKEEP-NCONV, EIGEN, IWORK,
     &                  N, TRUST, MAXCND, ABSDIS, EPS,
     &                  SAFMIN, LU)
         else if (WHICH.eq.-1) then
         end if
         rsave(IPOLE1) = POLE1
         rsave(IPOLE2) = POLE2
         rsave(ISIGOL) = SIGOLD
         rsave(IRBOUN) = RBOUND
         rsave(ILBOUN) = LBOUND
         isave(LINER1) = INERT1
         isave(LINER2) = INERT2
         isave(LTRIAL) = TRIAL
         if (ALLOW2) then
            isave(IALLOW) = 1
         else
            isave(IALLOW) = 0
         end if
         if (SIGMA.eq.rsave(ISIGOL)) then
            IDO = 100
         else
            IDO = 3
         end if
         return
      else if (IDO.eq.3) then
         call EA17OD(MODE, HARMON, NRITZ, RITZ, rsave(ISIGOL),
     &               SIGMA, ISINF, GROW, MAXDST, SAFMIN, LARGE)
         EXPRST = .false.
         ABSDIS = rsave(IABSDI)
         if (GROW.gt.MAXCND) then
            if (MAXDST.lt.ABSDIS) then
               if (USER) then
                  if (LUWARN.ge.0) write(LUWARN,8004)
                  call EA16ZD(IERR, 2**6)
               end if
               SIGMA = rsave(ISIGOL)
               if (.not.FIRST) then
                  rsave(IPOLE2) = SIGMA
                  isave(LINER2) = EA18ZD(NEINEG, MODE, SIGMA, N)
                  IDO = 100
                  return
               end if
            else
               EXPRST = .true.
            end if
         end if
 8004    format(/,'WARNING : The pole is not altered.',/)
         call EA17WD(rsave(ISIGOL), SIGMA, MODE, ISINF, ALPHA, BETA,
     &               GAMMA, DELTA, SAFMIN, IERR2)
         if (IERR2.eq.0 .or. FIRST) then
            IDO = 4
         else
            if (USER) then
               if (LUWARN.ge.0) write(LUWARN,8004)
               call EA16ZD(IERR, 2**6)
            end if
            SIGMA = rsave(ISIGOL)
            IDO = 100
         end if
      else if (IDO.eq.-4) then
         if (FIRST.and..not.ISINF) then
            IDO = 100
            IERR = -16
         else if (ISINF) then
            IDO = 2
            go to 1
         else
            if (SIGMA.eq.rsave(ISIGOL)) then
               IDO = 100
               IERR = -16
               return
            end if
            SIGMA = rsave(ISIGOL)
            GROW = ONE
            IDO = 4
         end if
      else if (IDO.eq.4) then
         call EA17OD(MODE, HARMON, NRITZ, RITZ, rsave(ISIGOL),
     &               SIGMA, ISINF, GROW, MAXDST, SAFMIN, LARGE)
         if (FIRST .or. SIGMA.ne.rsave(ISIGOL)) then
            IDO = 6
         else
            IDO = 100
         end if
         if (FIRST .or. SIGMA.ne.rsave(ISIGOL)) then
            isave(LINER1) = isave(LINER2)
            isave(LINER2) = EA18ZD(NEINEG, MODE, SIGMA, N)
            rsave(IPOLE1) = rsave(IPOLE2)
            rsave(IPOLE2) = SIGMA
         else if (SIGMA.eq.rsave(ISIGOL)) then
            rsave(IPOLE2) = SIGMA
            isave(LINER2) = EA18ZD(NEINEG, MODE, SIGMA, N)
         end if
         if (LU.ge.0) write(LU,8000) SIGMA
 8000    format('  Selected pole : ', 1PE14.6)
      end if
      return
      end
      subroutine EA17MD(CHECK, NVAL, WHICH, MODE, HARMON, POLE, RANGE,
     &                  NWANT, NWACO, NCOLEF, NCORIG, COUNT, NKRYL,
     &                  RITZ, WORK, SAFMIN, LU, IERR, MVAL)
      logical CHECK
      integer NVAL, WHICH, MODE, HARMON, NWANT, NWACO, NCOLEF, NCORIG
      integer LU, IERR, NKRYL, MVAL
      integer COUNT(NWACO)
      double precision SAFMIN
      double precision RANGE(2), POLE, RITZ(NKRYL), WORK(NKRYL)
      logical LEFT, RIGHT
      integer I, ILEFT, IRIGHT
      external EA17PD, EA16ZD
      intrinsic min
      if (WHICH.eq.1 .or. WHICH.eq.-1 .or. WHICH.eq.2 .or.
     &    WHICH.eq.-2 .or. WHICH.eq.10) then
         NVAL = min(NWANT,NWACO)
         MVAL = NVAL
         CHECK = NWACO.ge.NWANT
         return
      else if (WHICH.eq.3) then
         NVAL = min(NWANT,NWACO)
         MVAL = NVAL
         CHECK = NCOLEF.ge.NWANT/2 .and. NCORIG.ge.(NWANT+1)/2
         return
      else if (WHICH.eq.-3) then
         NVAL = min(NWANT,NWACO)
         MVAL = NVAL
         CHECK = NCORIG.ge.NWANT/2 .and. NCOLEF.ge.(NWANT+1)/2
         return
      end if
      do 10 I = 1,NWACO
         WORK(I) = RITZ(COUNT(I))
 10   continue
      call EA17PD(MODE, HARMON, NWACO, WORK, POLE, RANGE, SAFMIN)
      if (WHICH.eq.5) then
         LEFT = .false.
         RIGHT = .false.
         NVAL = 0
         do 50 I = 1,NWACO
            if (WORK(I).lt.RANGE(1)) then
               LEFT = .true.
            else if (WORK(I).gt.RANGE(2)) then
               RIGHT = .true.
            else
               NVAL = NVAL + 1
            end if
 50      continue
         ILEFT = 0
         IRIGHT = 0
         do 52 I = NWACO+1,NKRYL
            if (WORK(I).lt.RANGE(1)) then
               ILEFT = ILEFT + 1
            else if (WORK(I).gt.RANGE(2)) then
               IRIGHT = IRIGHT + 1
            end if
 52      continue
         if (ILEFT.eq.0) LEFT = .true.
         if (IRIGHT.eq.0) RIGHT = .true.
         CHECK = LEFT .and. RIGHT
         MVAL = NVAL
      else if (WHICH.eq.4) then
         CHECK = .false.
         NVAL = 0
         do 60 I = 1,NWACO
            if (WORK(I).le.RANGE(1)) then
               NVAL = NVAL + 1
            else
               CHECK = .true.
            end if
 60      continue
         MVAL = NVAL
      else if (WHICH.eq.-4) then
         CHECK = .false.
         NVAL = 0
         do 65 I = 1,NWACO
            if (WORK(I).ge.RANGE(1)) then
               NVAL = NVAL + 1
            else
               CHECK = .true.
            end if
 65      continue
         MVAL = NVAL
      end if
      NVAL = min(NVAL,NWANT)
      CHECK = CHECK .or. MVAL.ge.NWANT
      CHECK = CHECK .and. NWACO.ge.NWANT
      if (CHECK .and. NVAL.lt.NWANT) then
         if (LU.ge.0) write(LU,9000)
         call EA16ZD(IERR, 2**4)
      end if
 9000 format('WARNING : the required number of eigenvalues computed.',/,
     &   'But the number satisfying WHICH is smaller than NWANT.')
      return
      end
      subroutine EA17OD(MODE, HARMON, NKRYL, RITZ, OLDPOL,
     &                  NEWPOL, ISINF, GROW, MAXDST, SAFMIN, LARGE)
      logical ISINF
      integer MODE, HARMON, NKRYL
      double precision OLDPOL, NEWPOL, GROW, MAXDST, SAFMIN, LARGE
      double precision RITZ(*)
      double precision ZERO, ONE
      parameter (ZERO=0.0D0, ONE=1.0D0)
      external EA17WD
      logical DOINF
      integer I, IERR
      double precision MAXVAL, ALPHA, BETA, GAMMA, DELTA
      double precision NOMIN, DENOM
      intrinsic abs, max
      if (NKRYL.le.0) then
         GROW = ONE
         return
      end if
      call EA17WD(OLDPOL, NEWPOL, MODE, ISINF, ALPHA, BETA,
     &            GAMMA, DELTA, SAFMIN, IERR)
      DOINF = .false.
      if (MODE.eq.2 .or. MODE.eq.4 .or. MODE.eq.5) then
         DOINF = HARMON.ne.0
      end if
      if (IERR.gt.2) then
         GROW = LARGE
         return
      end if
      MAXVAL = ZERO
      if (ISINF) then
         do 10 I = 1,NKRYL
            DENOM = abs(RITZ(I)-NEWPOL)
            NOMIN = abs(RITZ(I))
            if (DENOM.le.SAFMIN*NOMIN) then
               MAXVAL = ONE / SAFMIN
            else
               MAXVAL = max(NOMIN/DENOM, MAXVAL)
            end if
 10      continue
      else if (DOINF) then
         do 11 I = 1,NKRYL
            DENOM = abs(RITZ(I)-NEWPOL)
            NOMIN = abs(RITZ(I)-OLDPOL)
            if (DENOM.le.SAFMIN*NOMIN) then
               MAXVAL = ONE / SAFMIN
            else
               MAXVAL = max(NOMIN/DENOM, MAXVAL)
            end if
 11      continue
      else
         do 12 I = 1,NKRYL
            DENOM = max(SAFMIN, abs( ALPHA + BETA * RITZ(I) ))
            MAXVAL = max(ONE/DENOM, MAXVAL)
 12      continue
      end if
      GROW = MAXVAL
      MAXDST = ZERO
      do 13 I = 1,NKRYL
         MAXDST = max(MAXDST, abs(NEWPOL-RITZ(I)))
 13   continue
      return
      end
      subroutine EA17PD(MODE, HARMON, N, RITZ, SIGMA, RANGE, SAFMIN)
      integer MODE, HARMON, N
      double precision SIGMA, SAFMIN
      double precision RANGE(2), RITZ(N)
      logical BUCK, SHIINV
      integer I
      double precision VALUE, TARGET
      double precision ONE
      parameter (ONE=1.0D0)
      intrinsic sign, abs
      if (MODE.eq.1 .or. MODE.eq.3) then
         return
      else if (MODE.eq.2 .or. MODE.eq.4) then
         if (HARMON.ne.0) return
         SHIINV = .true.
         BUCK = .false.
         TARGET = SIGMA
      else if (MODE.eq.5) then
         if (HARMON.ne.0) return
         SHIINV = .false.
         BUCK = .true.
      end if
      if (SHIINV) then
         do 10 I = 1,N
            VALUE = RITZ(I)
            if (abs(VALUE).lt.SAFMIN) then
               RITZ(I) = sign(ONE/SAFMIN, VALUE)
            else
               RITZ(I) = TARGET + sign(ONE/abs(VALUE), VALUE)
            end if
 10      continue
      elseif (BUCK) then
         do 20 I = 1,N
            VALUE = RITZ(I) - ONE
            if (abs(VALUE).le.SAFMIN*abs(SIGMA)) then
               RITZ(I) = sign(ONE/SAFMIN, VALUE) * sign(ONE, SIGMA)
            else
               RITZ(I) = SIGMA + SIGMA/VALUE
            end if
 20      continue
      end if
      return
      end
      subroutine EA17QD(OLDPOL, NEWPOL, INFLAN, INFLOC, MODE, NKRYL,
     &                  NKEEP, BLK, NCONV, OFFSET, RITZ, RESID, LRESID,
     &                  T, LDT, LCT, L, LDL, LCL, F, LDF, LCF, Z, LDZ,
     &                  LCZ, NZ, TAU, LTAU, WORK, LWORK, SAFMIN)
      logical INFLAN, INFLOC
      integer MODE, NKRYL, NKEEP, NCONV, OFFSET, BLK, LDT, LDL, LDF
      integer LDZ, LCZ, LRESID, LCL, LCF, LTAU, LWORK, LCT, NZ
      double precision OLDPOL, NEWPOL, SAFMIN
      double precision RITZ(LRESID), RESID(LRESID), T(LDT,LCT)
      double precision L(LDL,LCL), F(LDF,LCF)
      double precision Z(LDZ,LCZ), TAU(LTAU)
      double precision WORK(LWORK)
      integer IERR2
      double precision ALPHA, BETA, GAMMA, DELTA
      external DCOPY, EA17VD, EA17UD, EA17WD
      IERR2 = 0
      if (NKEEP.gt.NCONV) then
         call EA17WD(OLDPOL, NEWPOL, MODE, INFLAN, ALPHA, BETA,
     &            GAMMA, DELTA, SAFMIN, IERR2)
         if (IERR2.eq.0) then
            call EA17VD(T(1,NCONV-OFFSET+1), LDT, LCT-NCONV+OFFSET,
     &               NKEEP-NCONV, BLK, Z, LDZ, LCZ, NZ, ALPHA, BETA,
     &               GAMMA, DELTA, L, LDL, LCL, F, LDF, LCF, .true.,
     &               TAU, LTAU, WORK, LWORK, SAFMIN, IERR2)
         end if
      end if
      if (IERR2.eq.0) then
         call EA17UD(RITZ, 1, NKRYL, OLDPOL, NEWPOL, MODE, INFLAN,
     &               SAFMIN)
         call DCOPY(NCONV-OFFSET, RITZ(OFFSET+1),1, T,LDT)
      end if
      return
      end
      subroutine EA17RD(MODE, HARMON, DIR, RANGE, RITZ, NRITZ,
     &                  RESID, RESFAC, EIGEN, POLE, MAXCND, OLDPOL,
     &                  ISINF, SAFMIN, EPS, FIRST, LAST, STEP, ABSDIS,
     &                  TRIAL, IERR)
      logical ISINF
      integer MODE, HARMON, DIR, NRITZ, FIRST, LAST, STEP, TRIAL
      integer IERR
      double precision RESFAC, POLE, MAXCND, OLDPOL, ABSDIS, SAFMIN, EPS
      double precision RANGE(2), RITZ(NRITZ), EIGEN(NRITZ), RESID(NRITZ)
      integer I, I1
      double precision BOUND, LARGE, RELDIS
      double precision GROW, MAXDST
      double precision ONE
      parameter (ONE=1.0D0)
      external EA17OD
      IERR = 0
      LARGE = ONE/SAFMIN
      RELDIS = (ONE / MAXCND) * (ONE+sqrt(EPS))
      I1 = FIRST + TRIAL
      do 10 I = I1,LAST,STEP
         if (ISINF) then
            BOUND = RESFAC*RESID(I)
            BOUND = BOUND + abs(EIGEN(I)) * RELDIS
         else if (MODE.eq.2 .or. MODE.eq.4) then
            BOUND = abs(EIGEN(I)-OLDPOL) * RELDIS
         else if (MODE.eq.5) then
            BOUND = abs(EIGEN(I)-OLDPOL) * RELDIS
         end if
         BOUND = max(BOUND, ABSDIS) * (ONE+10*EPS)
         TRIAL = TRIAL + 1
         POLE = EIGEN(I) - DIR*BOUND
         if (MODE.eq.5 .and. abs(POLE).le.SAFMIN) go to 10
         call EA17OD(MODE, HARMON, NRITZ, RITZ, OLDPOL, POLE,
     &               ISINF, GROW, MAXDST, SAFMIN, LARGE)
         if (GROW.le.MAXCND) go to 11
 10   continue
      POLE = OLDPOL
      IERR = -16
 11   continue
      return
      end
      subroutine EA17SD(MODE, HARMON, OLDPOL, POLE, ISINF, RITZ, NRITZ,
     &                  EIGEN, REF, FIRST, LAST, STEP, DIR,
     &                  LBOUND, ABSDIS, MAXCND, SAFMIN, EPS)
      logical ISINF
      integer MODE, HARMON, NRITZ, REF, FIRST, LAST, STEP, DIR
      double precision OLDPOL, POLE, ABSDIS, MAXCND, LBOUND, SAFMIN, EPS
      double precision RITZ(NRITZ), EIGEN(NRITZ)
      integer I, IMAX
      double precision MINPOL, MAXPOL, GROW, LARGE
      double precision DIST, EIG1, EIG2, DISTMX, DISFAC
      double precision ONE, ZETMIN, AHALF
      parameter (ONE=1.0D0, ZETMIN=100.0D0, AHALF=0.5D0)
      external EA17OD
      intrinsic sqrt
      LARGE = ONE/SAFMIN
      DISFAC = ONE + sqrt(EPS)
      do 10 I = FIRST,LAST-STEP,STEP
         if (STEP.gt.0) then
            EIG1 = EIGEN(I)
            EIG2 = EIGEN(I+STEP)
         else
            EIG1 = EIGEN(I+STEP)
            EIG2 = EIGEN(I)
         end if
         DIST = max(abs(OLDPOL-EIG1) / MAXCND, ABSDIS)
         DIST = DIST*(ONE+2*EPS)
         MINPOL = EIG1+DIR*DIST
         DIST = max(abs(OLDPOL-EIG2) / MAXCND, ABSDIS)
         DIST = DIST*DISFAC
         MAXPOL = EIG2-DIR*DIST
         MINPOL = max(MINPOL,(ZETMIN*EIG1-EIGEN(REF))/(ZETMIN-1))
         MAXPOL = min(MAXPOL,
     &                (ZETMIN*EIG2+EIGEN(REF))/(ZETMIN+1))
         if (DIR*MINPOL.gt.DIR*MAXPOL) go to 10
         POLE = AHALF * (MINPOL+MAXPOL)
         if (DIR*POLE.lt.DIR*LBOUND) go to 10
         call EA17OD(MODE, HARMON, NRITZ, RITZ, OLDPOL, POLE,
     &               ISINF, GROW, DISTMX, SAFMIN, LARGE)
         if (GROW .le. MAXCND) go to 11
 10   continue
      MINPOL = LBOUND
      MAXPOL = EIGEN(FIRST)
      do 30 I = 1,NRITZ
         EIG1 = min(LBOUND, EIGEN(I))
         EIG2 = max(LBOUND, EIGEN(I))
         DIST = max(abs(OLDPOL-EIG2) / MAXCND, ABSDIS)
         DIST = DIST*DISFAC
         MAXPOL = min(MAXPOL, EIG2-DIR*DIST)
 30   continue
      if (DIR*LBOUND.le.DIR*MAXPOL) then
         POLE = MAXPOL
         MINPOL = POLE
         call EA17OD(MODE, HARMON, NRITZ, RITZ, OLDPOL, POLE,
     &               ISINF, GROW, DISTMX, SAFMIN, LARGE)
         if (GROW .le. MAXCND) go to 11
      end if
      DISTMX = abs(EIGEN(FIRST)-POLE)
      IMAX = FIRST-1
      do 20 I = FIRST,LAST-STEP,STEP
         DIST = abs(EIGEN(I)-EIGEN(I+1))
         if (DIST.gt.DISTMX) then
            DISTMX = DIST
            IMAX = I
         end if
 20   continue
      if (DISTMX.gt.0) then
         if (IMAX.eq.FIRST-1) then
            POLE = AHALF * ( POLE + EIGEN(FIRST) )
         else
            POLE = AHALF * ( EIGEN(IMAX) + EIGEN(IMAX+1) )
         end if
         if (DIR*POLE.lt.DIR*LBOUND) then
            POLE = OLDPOL
         end if
         go to 11
      end if
         MINPOL = LBOUND
         MAXPOL = LBOUND
         POLE = LBOUND
 11   continue
      if (DIR*OLDPOL.lt.DIR*MAXPOL .and. DIR*OLDPOL.gt.DIR*MINPOL) then
         POLE = OLDPOL
      end if
      return
      end
      subroutine EA17TD(NKRYL, NCONV, RITZ, RESID, OLDPOL,
     &                  NEWPOL, MODE, ISINF)
      logical ISINF
      integer NKRYL, NCONV, MODE
      double precision OLDPOL, NEWPOL
      double precision RITZ(NKRYL), RESID(*)
      double precision ONE
      parameter (ONE=1.0D0)
      integer I
      double precision VALUE, DIFF
      DIFF = NEWPOL-OLDPOL
      do 20 I = 1,NKRYL
         if (ISINF) then
            VALUE = abs(RITZ(I))
         else if (MODE.eq.4 .or. MODE.eq.2) then
            VALUE = abs(ONE + DIFF * RITZ(I))
         else if (MODE.eq.5) then
            VALUE = abs( (NEWPOL/OLDPOL) * RITZ(I) - RITZ(I) + ONE )
         end if
         if (I.lt.NCONV) RESID(I) = RESID(I) * VALUE
 20   continue
      return
      end
      subroutine EA17UD(E, LDE, N, OLDPOL, NEWPOL, MODE, ISINF, SAFMIN)
      logical ISINF
      integer LDE, N, MODE
      double precision OLDPOL, NEWPOL, SAFMIN
      double precision E(LDE,N)
      integer I
      double precision VALUE, DIFF, THETA, NOMIN, DENOM
      double precision ONE
      parameter (ONE=1.0D0)
      intrinsic max, sign, abs
      if (ISINF) then
         do 20 I = 1,N
            VALUE = (E(1,I) - NEWPOL)
            THETA = sign( ONE / max(SAFMIN,abs(VALUE)) , VALUE )
            if (MODE.eq.5) THETA = THETA * E(1,I)
            E(1,I ) = THETA
 20      continue
      else
         if (MODE.eq.5) then
            NOMIN = OLDPOL
            DENOM = NEWPOL
         else
            NOMIN = ONE
            DENOM = ONE
         end if
         DIFF = OLDPOL - NEWPOL
         do 10 I = 1,N
            THETA = E(1,I)
            VALUE = DENOM + THETA * DIFF
            THETA = THETA*NOMIN
            if (abs(VALUE).lt.abs(THETA)) then
               VALUE = VALUE / THETA
               THETA = sign(ONE/max(SAFMIN,abs(VALUE)), VALUE)
            else
               THETA = THETA / VALUE
            end if
            E(1,I ) = THETA
 10      continue
      end if
      return
      end
      subroutine EA17VD(T, LDT, LCT, NKRYL, BLK, Z, LDZ, LCZ, NZ, ALPHA,
     &                  BETA, GAMMA, DELTA, L, LDL, LCL, F, LDF, LCF,
     &                  KEEPZ, TAU, LTAU, WORK, LWORK, SAFMIN, IERR)
      logical KEEPZ
      integer LDT, LCT, NKRYL, BLK, LDZ, LCZ, NZ, LDL, LCL, LDF, LCF
      integer LWORK, LTAU, IERR
      double precision ALPHA, BETA, GAMMA, DELTA, SAFMIN
      double precision T(LDT,LCT), Z(LDZ,LCZ), L(LDL,LCL)
      double precision F(LDF,LCF)
      double precision TAU(LTAU), WORK(LWORK)
      double precision ONE, ZERO
      parameter (ONE=1.0D0, ZERO=0.0D0)
      integer ICOL, IROW, J, JT, LDIAG, ZLEN, FLEN, HLEN, LLEN
      external EA17YD, EA17ZD, DLARFG, DTBSV
      LDIAG = 2*BLK+1
      do 10 ICOL = 1,NKRYL
         do 100 J = 1,BLK
            L(J,ICOL) = ZERO
 100     continue
         do 150 J = 0,BLK
            L(LDIAG+J,ICOL) = BETA * T(1+J,ICOL)
 150     continue
         L(LDIAG,ICOL) = L(LDIAG,ICOL) + ALPHA
 10   continue
      do 20 ICOL = 1,NKRYL
         do 200 J = max(1,ICOL-BLK),ICOL-1
            L(J-ICOL+LDIAG,ICOL) = L(ICOL-J+LDIAG,J)
 200     continue
 20   continue
      do 30 ICOL = 1,NKRYL
         HLEN = BLK+1
         call DLARFG(HLEN, L(LDIAG,ICOL), L(LDIAG+1,ICOL),1, TAU(ICOL))
         if (ICOL.lt.NKRYL) then
            LLEN = min(2*BLK,NKRYL-ICOL)
            call EA17YD(HLEN, LLEN, L(LDIAG+1,ICOL), 1,
     &                  TAU(ICOL), L(LDIAG-1,ICOL+1), LDL-1,
     &                  L(LDIAG,ICOL+1), LDL-1, WORK)
         end if
 30   continue
      FLEN = NKRYL+BLK
      do 40 ICOL = 1,NKRYL
         do 400 IROW=1,FLEN
            F(IROW,ICOL) = ZERO
 400     continue
         F(ICOL,ICOL) = ONE
 40   continue
      IERR = 0
      do 45 ICOL = 1,NKRYL
         if (abs(L(LDIAG,ICOL)).lt.SAFMIN) IERR=ICOL
 45   continue
      if (IERR.ne.0) return
      HLEN = BLK+1
      do 50 ICOL = 1,NKRYL
         FLEN =  min(ICOL+BLK,NKRYL)
         call EA17YD(HLEN, FLEN, L(LDIAG+1,ICOL), 1, TAU(ICOL),
     &               F(ICOL,1), LDF, F(ICOL+1,1), LDF, WORK)
 50   continue
      if (KEEPZ) then
         ZLEN =  NZ
         do 55 ICOL = 1,NKRYL
            call EA17ZD(BLK+1, ZLEN, L(LDIAG+1,ICOL),1, TAU(ICOL),
     &                  Z(1,ICOL),1, Z(1,ICOL+1),LDZ, WORK)
 55      continue
      end if
      do 60 IROW=1,NKRYL+BLK
         call DTBSV('Upper', 'Tranpose', 'Non-unit_diag', NKRYL, 2*BLK,
     &              L,LDL, F(IROW,1),LDF)
 60   continue
      do 70 IROW=NKRYL+BLK,BLK+2,-1
         HLEN = IROW-BLK
         call DLARFG(HLEN, F(IROW,HLEN), F(IROW,1),LDF, TAU(HLEN))
         FLEN = min(NKRYL,IROW-1)
         call EA17YD(HLEN, FLEN, F(IROW,1),LDF, TAU(HLEN),
     &               F(HLEN,1), LDF, F,LDF, WORK)
         FLEN = IROW-1
         call EA17ZD(HLEN, FLEN, F(IROW,1),LDF, TAU(HLEN),
     &               F(1,HLEN),1, F(1,1),LDF, WORK)
 70   continue
      do 80 ICOL = 1,NKRYL
         JT = 1
         do 800 J = ICOL,ICOL+BLK
            T(JT,ICOL) = GAMMA * F(J,ICOL)
            JT = JT + 1
 800     continue
         T(1,ICOL) = T(1,ICOL) + DELTA
 80   continue
      if (KEEPZ) then
         FLEN = NZ
         do 90 IROW=NKRYL+BLK,BLK+2,-1
            HLEN = IROW-BLK
            call EA17ZD(HLEN, FLEN, F(IROW,1),LDF, TAU(HLEN),
     &                  Z(1,HLEN),1, Z(1,1),LDZ, WORK)
 90      continue
      end if
      return
      end
      subroutine EA17WD(OLDPOL, NEWPOL, MODE, ISINF, ALPHA, BETA,
     &                  GAMMA, DELTA, SAFMIN, IERR)
      logical ISINF
      integer MODE, IERR
      double precision OLDPOL, NEWPOL, SAFMIN, ALPHA, BETA, GAMMA, DELTA
      double precision ZERO, ONE
      parameter (ZERO=0.0D0, ONE=1.0D0)
      double precision EA17XD
      external EA17XD
      IERR = 0
      if (ISINF) then
         if (MODE.eq.5) then
            ALPHA = -NEWPOL
            BETA = ONE
            GAMMA = NEWPOL
            DELTA = ONE
         else
            ALPHA = -NEWPOL
            BETA = ONE
            GAMMA = ONE
            DELTA = ZERO
         end if
      else
         if (MODE.eq.5) then
            ALPHA = EA17XD(NEWPOL, OLDPOL, SAFMIN, IERR)
            if (IERR.ne.0) then
               IERR = -1
               return
            end if
            BETA = EA17XD(OLDPOL-NEWPOL, OLDPOL, SAFMIN, IERR)
            if (IERR.ne.0) then
               IERR = -2
               return
            end if
         else
            ALPHA = ONE
            BETA = OLDPOL-NEWPOL
         end if
         GAMMA = -EA17XD(ALPHA, BETA, SAFMIN, IERR)
         if (IERR.ne.0) then
            IERR = -3
            return
         end if
         DELTA = EA17XD(ONE, BETA, SAFMIN, IERR)
         if (IERR.ne.0) then
            IERR = -4
            return
         end if
      end if
      return
      end
      double precision function EA17XD(NOMIN, DENOM, SAFMIN, IERR)
      double precision NOMIN, DENOM, SAFMIN
      integer IERR
      double precision T
      double precision ONE
      parameter (ONE=1.0D0)
      IERR = 0
      if (abs(DENOM).ge.abs(NOMIN)) then
         T = NOMIN / DENOM
      else
         T = DENOM / NOMIN
         if (abs(T).lt.SAFMIN) then
            EA17XD = ONE
            IERR = -1
            return
         end if
         T = ONE / T
      end if
      EA17XD = T
      return
      end
      subroutine EA17YD(HLEN, NCOLS, V, INCV, TAU, A1, INCA, A2,LDA,
     &                  WORK)
      integer HLEN, NCOLS, INCV,  INCA, LDA
      double precision TAU
      double precision V(INCV*(HLEN-2)+1), A1(INCA*(NCOLS-1)+1)
      double precision A2(LDA,*), WORK(NCOLS)
      double precision ONE
      parameter (ONE=1.0D0)
      integer I, IA
      external DGEMV, DGER
      IA = 1
      do 10 I = 1,NCOLS
         WORK(I) = A1(IA)
         IA = IA + INCA
 10   continue
      call DGEMV('Transpose', HLEN-1, NCOLS, ONE, A2,LDA,
     &           V,INCV, ONE, WORK,1)
      IA = 1
      do 30 I = 1,NCOLS
         A1(IA) = A1(IA) - TAU * WORK(I)
         IA = IA + INCA
 30   continue
      call DGER(HLEN-1, NCOLS, -TAU, V,INCV, WORK,1, A2,LDA)
      return
      end
      subroutine EA17ZD(HLEN, NCOLS, V, INCV, TAU, A1, INCA, A2,LDA,
     &                  WORK)
      integer HLEN, NCOLS, INCV,  INCA, LDA
      double precision TAU
      double precision V(INCV*(HLEN-2)+1), A1(INCA*(NCOLS-1)+1)
      double precision A2(LDA,HLEN-1), WORK(*)
      double precision ONE
      parameter (ONE=1.0D0)
      integer I, IA
      external DGEMV, DGER
      IA = 1
      do 10 I = 1,NCOLS
         WORK(I) = A1(IA)
         IA = IA + INCA
 10   continue
      call DGEMV('No transpose', NCOLS, HLEN-1, ONE,
     &           A2,LDA, V,INCV, ONE, WORK,1)
      IA = 1
      do 30 I = 1,NCOLS
         A1(IA) = A1(IA) - TAU * WORK(I)
         IA = IA + INCA
 30   continue
      call DGER(NCOLS, HLEN-1, -TAU, WORK,1, V,INCV, A2,LDA)
      return
      end
      subroutine EA18CD(HV, LDVH, TAU, BLK, NT, NSTEPS, Z, LDZ, LCZ, NZ,
     &                  WORK, LWORK)
      integer LDVH, BLK, NT, NSTEPS, LDZ, LCZ, NZ, LWORK
      double precision HV(LDVH,NT*NSTEPS), TAU(NT*NSTEPS), Z(LDZ,LCZ),
     &                 WORK(LWORK)
      integer ISTEP, ICOL, ITR, HLEN, ICOL1
      external EA17YD
      HLEN = BLK+1
      ICOL1 = NT - (NSTEPS-1)*BLK
      do 10 ISTEP=NSTEPS,1,-1
         ITR = NT*ISTEP + ICOL1 - NT
         do 100 ICOL = ICOL1,1,-1
            call EA17YD(HLEN, NZ, HV(1,ITR),1, TAU(ITR),
     &                  Z(ICOL,1),LDZ, Z(ICOL+1,1),LDZ, WORK)
            ITR = ITR - 1
 100     continue
         ICOL1 = ICOL1 + BLK
 10   continue
      return
      end
      subroutine EA18DD(T, LDT, N, BLK, SHIFT, Z, LDZ, NZ, LCZ, TT,
     &                  LDTT, LCTT, KEEP, HV, LDHV, LCHV, TAU, LTAU,
     &                  WORK, LWORK)
      logical KEEP
      integer N, LDT, BLK, NZ, LDZ, LCZ, LDTT, LCTT, LDHV, LCHV, LTAU
      integer LWORK
      double precision SHIFT
      double precision T(LDT,N), TT(LDTT,LCTT), HV(LDHV,LCHV)
      double precision TAU(LTAU), Z(LDZ,LCZ), WORK(LWORK)
      integer I, ICOL, IDIAG, ITR, J, HLEN, LLEN
      double precision ZERO
      parameter (ZERO=0.0D0)
      external DLASET, DCOPY, DLARFG, EA17YD, EA17ZD, DSCAL
      IDIAG = 2*BLK+1
      call DLASET('All', LDTT, N, ZERO, ZERO, TT,LDTT)
      do 30 ICOL = 1,N
         do 300 I = 1,BLK
            TT(IDIAG+I,ICOL) = T(1+I,ICOL)
 300     continue
         TT(IDIAG,ICOL) = T(1,ICOL) - SHIFT
         do 350 I = 1,min(BLK,N-ICOL)
            TT(IDIAG-I,ICOL+I) = T(1+I,ICOL)
 350     continue
 30   continue
      call DLASET('All', BLK, BLK, ZERO, ZERO, HV,LDHV)
      do 50 I = 1,BLK
         TAU(I) = ZERO
 50   continue
      if (KEEP) then
         ITR = 1
      else
         ITR = BLK+1
      end if
      do 10 ICOL = 1,N
         HLEN = BLK+1
         call DCOPY(HLEN-1, TT(IDIAG+1,ICOL),1, HV(1,ITR),1)
         call DLARFG(HLEN, TT(IDIAG,ICOL), HV(1,ITR),1, TAU(ITR))
         call DSCAL(BLK, ZERO, TT(IDIAG+1,ICOL),1)
         if (ICOL.lt.N) then
            LLEN = min(2*BLK,N-ICOL)
            call EA17YD(HLEN, LLEN, HV(1,ITR),1, TAU(ITR),
     &                  TT(IDIAG-1,ICOL+1),LDTT-1, TT(IDIAG,ICOL+1),
     &                  LDTT-1, WORK)
         end if
         call EA17ZD(HLEN, NZ, HV(1,ITR),1, TAU(ITR),
     &               Z(1,ICOL),1, Z(1,ICOL+1),LDZ, WORK)
         if (ICOL.ge.BLK+1) then
            call EA17ZD(HLEN, 2*BLK+1, HV(1,ITR-BLK),1, TAU(ITR-BLK),
     &                  TT(BLK+1,ICOL-BLK),1,
     &                  TT(BLK,ICOL-BLK+1),LDTT-1, WORK)
         end if
         if (KEEP) then
            ITR = ITR + 1
         else
            do 100 I = 2,ITR
               TAU(I-1) = TAU(I)
               do 1000 J = 1,HLEN-1
                  HV(J,I-1) = HV(J,I)
 1000          continue
 100        continue
         end if
 10   continue
      do 40 ICOL = 1,N-BLK
         T(1,ICOL) = TT(IDIAG,ICOL) + SHIFT
         do 400 I = 1,BLK
            T(1+I,ICOL) = TT(IDIAG+I,ICOL)
 400     continue
 40   continue
      return
      end
      double precision function EA18ED(M, N, NB, E, A, LDA)
      integer M, N, NB, LDA
      double precision E(M), A(LDA,N)
      integer I, J
      double precision NORM, MMAX
      double precision ZERO
      parameter       (ZERO=0.0D0)
      intrinsic abs, max, sqrt
      MMAX = ZERO
      do 15 I = 1,M
         MMAX = max(MMAX, abs(E(I)))
 15   continue
      do 10 I = M+1,N
         do 100 J = 1,NB+1
            MMAX = max(MMAX, abs(A(J,I)))
 100     continue
 10   continue
      if (MMAX.eq.ZERO) then
         EA18ED = ZERO
         return
      end if
      NORM = ZERO
      do 20 I = 1,M
         NORM = NORM + ( E(I)/MMAX ) ** 2
 20   continue
      do 25 I = M+1,N-NB
         NORM = NORM + ( A(1,I)/MMAX ) ** 2
         do 250 J = 2,NB+1
            NORM = NORM + 2 * ( A(J,I)/MMAX ) ** 2
 250     continue
 25   continue
      do 26 I = max(M+1,N-NB+1),N
         do 260 J = 1,NB+1
            NORM = NORM + ( A(J,I)/MMAX ) ** 2
 260     continue
 26   continue
      EA18ED = sqrt(NORM) * MMAX
      return
      end
      double precision function EA18FD(UPLO, M, N, A, LDA)
      integer M, N, LDA
      character UPLO
      double precision A(LDA,*)
      integer I, J, MN
      double precision NORM, MMAX
      double precision ZERO
      parameter       (ZERO=0.0D0)
      intrinsic max, sqrt
      if (UPLO.eq.'L') then
         MN = min(M,N)
         MMAX = ZERO
         do 20 I = 1,MN
            do 200 J = I,MN
               MMAX = max(MMAX, abs(A(J,I)))
 200        continue
 20      continue
         if (MMAX.eq.ZERO) then
            EA18FD = ZERO
            return
         end if
         NORM = ZERO
         do 25 I = 1,MN
            do 250 J = I,MN
               NORM = NORM + ( A(J,I)/MMAX ) **2
 250        continue
 25      continue
         NORM = sqrt(NORM) * MMAX
      else if (UPLO.eq.'B') then
         MMAX = ZERO
         do 30 I = 1,N
            do 300 J = 1,M
               MMAX = max(MMAX, abs(A(J,I)))
 300        continue
 30      continue
         if (MMAX.eq.ZERO) then
            EA18FD = ZERO
            return
         end if
         NORM = ZERO
         do 45 I = 1,N
            do 450 J = 1,M
               NORM = NORM + ( A(J,I)/MMAX ) **2
 450        continue
 45      continue
         NORM = sqrt(NORM) * MMAX
      end if
      EA18FD = NORM
      return
      end
      subroutine EA18GD(GAMMA, POLCHA, EPS, LARGE, NORMT, GROW)
      logical POLCHA
      double precision GAMMA, EPS, LARGE, NORMT, GROW
      intrinsic max, sqrt
      if (POLCHA) then
         if (GAMMA.lt.1) then
            if (GROW*sqrt(GAMMA).lt.sqrt(LARGE)) then
               GAMMA = (GAMMA * GROW) * GROW
            else
               GAMMA = LARGE
            end if
         else if (GROW.lt.sqrt(LARGE/GAMMA)) then
               GAMMA = GAMMA * GROW**2
         else
            GAMMA = LARGE
         end if
      else
         GAMMA = max(GAMMA, EPS*NORMT)
         GAMMA = EPS*NORMT
      end if
      return
      end
      subroutine EA18HD(NKRYL, NPURG, RITZ, RESID, V1, LDV1, NV1,
     &                  V2, LDV2, NV2, TOL, WORK, IWORK)
      integer NKRYL, NPURG, LDV1, NV1, LDV2, NV2
      integer IWORK(NKRYL)
      double precision TOL
      double precision RITZ(NKRYL), RESID(NKRYL), WORK(NKRYL)
      double precision V1(LDV1,NKRYL), V2(LDV2,NKRYL)
      integer I, J, IL, IS
      if (NKRYL.le.0) then
         NPURG = 0
         return
      end if
      NPURG = 0
      do 10 I = 1,NKRYL
         if (RESID(I).gt.TOL) NPURG = NPURG + 1
 10   continue
      IL = 1
      IS = NPURG+1
      do 20 I = 1,NKRYL
         if (RESID(I).gt.TOL) then
            IWORK(I) = IL
            IL = IL + 1
         else
            IWORK(I) = IS
            IS = IS + 1
         end if
 20   continue
      do 30 I = 1,NKRYL
         WORK(I) = RITZ(I)
 30   continue
      do 35 I = 1,NKRYL
         RITZ(I) = WORK(IWORK(I))
 35   continue
      do 50 J = 1,NV1
         do 500 I = 1,NKRYL
            WORK(I) = V1(J,I)
 500     continue
         do 550 I = 1,NKRYL
            V1(J,I) = WORK(IWORK(I))
 550     continue
 50   continue
      do 70 J = 1,NV2
         do 700 I = 1,NKRYL
            WORK(I) = V2(J,I)
 700     continue
         do 750 I = 1,NKRYL
            V2(J,I) = WORK(IWORK(I))
 750     continue
 70   continue
      do 60 I = NKRYL,1,-1
         if (RESID(I).gt.TOL) then
            NPURG = I
            go to 61
         end if
 60   continue
      NPURG = 0
 61   continue
      return
      end
      double precision function EA18ID(ISTEP, LSTEP, MSTEP, ERRMOD,
     &                                 LDERRM)
      integer ISTEP, LSTEP, MSTEP, LDERRM
      double precision ERRMOD(LDERRM,ISTEP)
      integer K
      EA18ID = 0.0D0
      do 70 K = LSTEP,MSTEP
         EA18ID = max(EA18ID, ERRMOD(K,ISTEP))
 70   continue
      return
      end
      subroutine EA18JD(T, LDT, BLK, F, LDF)
      integer LDT, BLK, LDF
      double precision T(LDT,*), F(LDF,*)
      integer ICOL, J, K, JROW
      do 10 ICOL = 1,BLK
         JROW = ICOL
         do 100 J = 1,BLK+1-ICOL
            do 1000 K = 1,ICOL
               T(J,ICOL) = T(J,ICOL) + F(JROW,K) * T(BLK+1+K-ICOL,ICOL)
 1000       continue
            JROW = JROW + 1
 100     continue
 10   continue
      do 20 JROW = 1,BLK
         do 200 ICOL = JROW,BLK
            T(BLK+1-ICOL+JROW,ICOL) =
     &              T(BLK+1-ICOL+JROW,ICOL) * F(BLK+JROW,JROW)
            do 2000 K = JROW+1,ICOL
               T(BLK+1-ICOL+JROW,ICOL) = T(BLK+1-ICOL+JROW,ICOL)
     &              + F(BLK+JROW,K) * T(BLK+1-ICOL+K,ICOL)
 2000       continue
 200     continue
 20   continue
      return
      end
      subroutine EA18KD(T, LDT, BLK, NCOLS, F, LDF)
      integer LDT, BLK, NCOLS, LDF
      double precision T(LDT,NCOLS), F(LDF,NCOLS)
      integer ICOL, J, JROW
      do 20 ICOL = 1,NCOLS
         JROW = ICOL
         do 200 J = 1,BLK+1
            T(J,ICOL) = F(JROW,ICOL)
            JROW = JROW + 1
 200     continue
 20   continue
      return
      end
      subroutine EA18LD(T, LDT, BLK, NCOLS, F, LDF)
      integer LDT, BLK, NCOLS, LDF
      double precision T(LDT,NCOLS), F(LDF,NCOLS)
      integer ICOL, J, JROW
      do 10 ICOL = 1,NCOLS
         do 100 JROW = 1,ICOL-1
            F(JROW,ICOL) = 0.0D0
 100     continue
         do 110 JROW = ICOL+BLK+1,NCOLS+BLK
            F(JROW,ICOL) = 0.0D0
 110     continue
 10   continue
      do 20 ICOL = 1,NCOLS
         JROW = ICOL
         do 200 J = 1,BLK+1
            F(JROW,ICOL) = T(J,ICOL)
            JROW = JROW + 1
 200     continue
 20   continue
      return
      end
      subroutine EA18MD(COUNT,N,INDEX)
      integer N
      double precision COUNT(N)
      integer INDEX(N)
      double precision AV,X
      integer I,IF,IFK,IFKA,INT,INTEST,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M,
     +        MLOOP
      integer MARK(50)
      IF (N.LE.1) GO TO 200
      M = 12
      LA = 2
      IS = 1
      IF = N
      DO 190 MLOOP = 1,N
        IFKA = IF - IS
        IF ((IFKA+1).GT.M) GO TO 70
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
      end
      subroutine EA18ND(N, V, W, NORM, RATIO, SAFMIN, INFO)
      integer N, INFO
      double precision NORM, RATIO, SAFMIN
      double precision V(N), W(N)
      integer IMAX, I
      double precision VSCAL, WSCAL
      double precision ZERO, ONE
      parameter (ZERO=0.0D0, ONE=1.0D0)
      integer IDAMAX
      external IDAMAX
      intrinsic abs, min, sqrt
      INFO = 0
      IMAX = IDAMAX(N, V, 1)
      VSCAL = abs(V(IMAX))
      if (VSCAL.eq.ZERO) then
         NORM = ZERO
         RATIO = ZERO
         return
      end if
      IMAX = IDAMAX(N, W, 1)
      WSCAL = abs(W(IMAX))
      if (WSCAL.eq.ZERO) then
         NORM = ZERO
         RATIO = ONE / SAFMIN
         return
      end if
      if (WSCAL.le.VSCAL*SAFMIN) then
         RATIO = ONE/SAFMIN
      else
         RATIO = VSCAL/WSCAL
      end if
      NORM = ZERO
      do 10 I = 1,N
         NORM = NORM + (V(I)/VSCAL) * (W(I)/WSCAL)
 10   continue
      if (NORM.lt.ZERO) then
         INFO = 1
         return
      end if
      NORM = sqrt(NORM) * sqrt(VSCAL) * sqrt(WSCAL)
      if (NORM.lt.SAFMIN*min(VSCAL,WSCAL)) NORM=ZERO
      return
      end
      subroutine EA18OD(NORD, V, LDV, NV, PERM, TEMP)
      integer NORD, LDV, NV
      integer PERM(NORD)
      double precision V(LDV,NORD)
      double precision TEMP(NORD)
      integer I, J
      do 50 J = 1,NV
         do 500 I = 1,NORD
            TEMP(I) = V(J,I)
 500     continue
         do 550 I = 1,NORD
            V(J,I) = TEMP(PERM(I))
 550     continue
 50   continue
      return
      end
      subroutine EA18PD(T, LDT, BLK, N, RITZ, PIPE, LDPIPE,
     &                  Q, LDQ, NQ, TT, LDTT)
      integer LDT, BLK, N, LDPIPE, LDQ, NQ, LDTT
      double precision T(LDT,N), RITZ(N), PIPE(LDPIPE,N)
      double precision Q(LDQ,N), TT(LDTT,N)
      integer I, ICOL, JCOL, IBLK, IDIAG, ICOL1, ICOL2, I1
      double precision SS, CC, RR, WW1, WW2
      double precision ONE, ZERO
      parameter (ONE=1.0D0, ZERO=0.0D0)
      external DLARTG
      IDIAG = BLK+2
      do 20 ICOL = 1,N
         do 200 I = 1,2*BLK+3
            TT(I,ICOL) = ZERO
 200     continue
         TT(IDIAG,ICOL) = RITZ(ICOL)
 20   continue
      do 10 ICOL = 1-BLK,N-BLK-1
         do 100 IBLK = BLK,max(1,-ICOL+1),-1
            ICOL2 = ICOL+IBLK+1
            ICOL1 = ICOL2-1
            call DLARTG(PIPE(IBLK,ICOL2),PIPE(IBLK,ICOL1),
     &                  CC, SS, RR)
            PIPE(IBLK,ICOL1) = ZERO
            PIPE(IBLK,ICOL2) = RR
            do 1000 I = 1,2*BLK+2
               WW1 = TT(I+1,ICOL1)
               WW2 = TT(I,ICOL2)
               TT(I+1,ICOL1) = WW1 * CC - WW2 * SS
               TT(I,ICOL2) = WW1 * SS + WW2 * CC
 1000       continue
            do 1010 I = 1,2*BLK+2
               I1 = ICOL1+BLK+2-I
               if (I1.le.0 .or. I1.gt.N) go to 1010
               WW1 = TT(I,I1)
               WW2 = TT(I+1,I1)
               TT(I,I1) = WW1 * CC - WW2 * SS
               TT(I+1,I1) = WW1 * SS + WW2 * CC
 1010       continue
            do 1020 I = 1,IBLK-1
               WW1 = PIPE(I,ICOL1)
               WW2 = PIPE(I,ICOL2)
               PIPE(I,ICOL1) = WW1 * CC - WW2 * SS
               PIPE(I,ICOL2)   = WW1 * SS + WW2 * CC
 1020       continue
            do 1030 I = 1,NQ
               WW1 = Q(I,ICOL1)
               WW2 = Q(I,ICOL2)
               Q(I,ICOL1) = WW1 * CC - WW2 * SS
               Q(I,ICOL2) = WW1 * SS + WW2 * CC
 1030       continue
 100     continue
         do 150 JCOL=ICOL,1,-1
            ICOL1 = JCOL
            ICOL2 = ICOL1+1
            call DLARTG(TT(2*BLK+2,ICOL2),TT(2*BLK+3,ICOL1), CC, SS, RR)
            if (JCOL.lt.ICOL-1 .and. CC.eq.ONE .and. SS.eq.ZERO) then
               go to 151
            end if
            TT(2*BLK+2,ICOL2) = RR
            TT(2*BLK+3,ICOL1) = ZERO
            do 1050 I = 1,2*BLK+1
               WW1 = TT(I+1,ICOL1)
               WW2 = TT(I,ICOL2)
               TT(I+1,ICOL1) = WW1 * CC - WW2 * SS
               TT(I,ICOL2) = WW1 * SS + WW2 * CC
 1050       continue
            do 1040 I = 1,2*BLK+2
               I1 = ICOL1+BLK+2-I
               if (I1.le.0 .or. I1.gt.N) go to 1040
               WW1 = TT(I,I1)
               WW2 = TT(I+1,I1)
               TT(I,I1) = WW1 * CC - WW2 * SS
               TT(I+1,I1) = WW1 * SS + WW2 * CC
 1040       continue
            do 1060 I = 1,NQ
               WW1 = Q(I,ICOL1)
               WW2 = Q(I,ICOL2)
               Q(I,ICOL1) = WW1 * CC - WW2 * SS
               Q(I,ICOL2) = WW1 * SS + WW2 * CC
 1060       continue
 150     continue
 151     continue
 10   continue
      IDIAG = BLK+1
      do 50 ICOL = 1,N
         do 500 I = 1,BLK+1
            T(I,ICOL) = TT(IDIAG+I,ICOL)
 500     continue
 50   continue
      do 60 ICOL = max(1,BLK-N+1),BLK
         do 600 IBLK = 1,ICOL
            JCOL = N-BLK+ICOL
            T(BLK+1-ICOL+IBLK,JCOL) = PIPE(IBLK,JCOL)
 600     continue
 60   continue
      return
      end
      subroutine EA18RD(IDO, WHICH, ALLOW2, PURIF5, N, RANGE, OK,
     &                  INERT1, INERT2, POLE1, POLE2, RBOUND, LBOUND,
     &                  SIGOLD, SIGMA, NEINEG, EIGEN, NEIGEN,
     &                  NCONV, NWANT, TRUST, ABSDIS, RELDIS, SIGTST,
     &                  NEITST, NBETST, LARGE, NORESC, LU)
      logical ALLOW2, OK, PURIF5, NORESC
      integer IDO, WHICH, INERT1, INERT2, NEINEG, NEIGEN, NCONV
      integer N, LU, NEITST, NBETST, NWANT
      double precision POLE1, POLE2, LBOUND, RBOUND, SIGOLD, SIGMA
      double precision ABSDIS, RELDIS, SIGTST, LARGE
      double precision EIGEN(NEIGEN), TRUST(4), RANGE(2)
      logical DO2, DOM2, DO4, DOM4
      integer NEIG, NBETWE, LAST, I
      double precision BRIGHT, BLEFT, EXTREM
      double precision ZERO, AHALF
      parameter (ZERO=0.0D0, AHALF=0.5D0)
      if (IDO.eq.-4) then
         if (SIGOLD.eq.SIGMA) then
            IDO = 100
         else
            SIGMA = SIGOLD
            IDO = 6
         end if
         OK = .true.
         return
      end if
      if (.not.ALLOW2) then
         DO2  = .false.
         DOM2 = .false.
         DO4  = WHICH.eq.4 .or. WHICH.eq.-2
         DOM4 = WHICH.eq.-4 .or. WHICH.eq.5 .or. WHICH.eq.2
      else
         DO2  = WHICH.eq.2
         DOM2 = WHICH.eq.-2
         DO4  = WHICH.eq.4
         DOM4 = WHICH.eq.-4 .or. WHICH.eq.5
      end if
      if (IDO.eq.0) then
      if (WHICH.eq.-1) then
         IDO = 100
         return
      end if
      if (INERT1.lt.-N .or. INERT1.gt.N) then
         POLE1 = POLE2
         INERT1 = INERT2
      end if
      BRIGHT = max(POLE1,POLE2)
      BLEFT = min(POLE1,POLE2)
      LAST = 0
      NBETWE = 0
      if (DO4) then
         EXTREM = RANGE(1)
         do 40 I = 1,NCONV
            if (EIGEN(I).le.RANGE(1) .and. EIGEN(I).gt.POLE1)
     &         LAST = I
            if (EIGEN(I).le.BRIGHT .and. EIGEN(I).ge.BLEFT)
     &         NBETWE = NBETWE + 1
            if (EIGEN(I).lt.EXTREM.and.I.le.NWANT) EXTREM = EIGEN(I)
 40      continue
         NEIG = abs(INERT2-INERT1)
      else if (DOM4) then
         EXTREM = RANGE(1)
         do 50 I = 1,NCONV
            if (EIGEN(I).ge.RANGE(1) .and. EIGEN(I).lt.POLE1)
     &         LAST = I
            if (EIGEN(I).le.BRIGHT .and. EIGEN(I).ge.BLEFT) then
               NBETWE = NBETWE + 1
            end if
            if (EIGEN(I).gt.EXTREM.and.I.le.NWANT) EXTREM = EIGEN(I)
 50      continue
         NEIG = abs(INERT2-INERT1)
      else if (DO2) then
         do 20 I = 1,NCONV
            if (EIGEN(I).le.BRIGHT) NBETWE = NBETWE + 1
 20      continue
         NEIG = max(INERT2,INERT1)
      else if (DOM2) then
         do 25 I = 1,NCONV
            if (EIGEN(I).ge.BLEFT) NBETWE = NBETWE + 1
 25      continue
         NEIG = N-min(INERT2,INERT1)
      end if
      if (LU.ge.0) write(LU,9000) NEIG, min(NEIG,NBETWE)
      OK = NEIG .eq. NBETWE
      if (PURIF5 .and. BRIGHT*BLEFT.lt.ZERO) then
         IDO = 100
         if (SIGMA.ne.SIGOLD) IDO = 6
         SIGMA = SIGOLD
         OK = .true.
         return
      end if
      if (WHICH.eq.5) then
         if (.not.OK .and. LAST+NBETWE.lt.NWANT .and.
     &       BRIGHT.le.RANGE(2)) then
            IDO = 100
            return
         end if
      else
         if (.not.OK .and. LAST+NBETWE.lt.NWANT) then
            IDO = 100
            return
         end if
      end if
      if (ALLOW2) then
         if (WHICH.eq.4) then
            if (LAST+NEIG.ge.NWANT) LBOUND = max(LBOUND,BLEFT)
         else if (WHICH.eq.-4 .or. WHICH.eq.5) then
            if (LAST+NEIG.ge.NWANT) RBOUND = min(RBOUND,BRIGHT)
         end if
      end if
      if (OK) then
         if (DO4 .or. DOM4 .or. WHICH.eq.-1) then
            TRUST(2) = min(BLEFT, TRUST(2))
            TRUST(3) = max(BRIGHT, TRUST(3))
            if (LAST+NEIG.ge.NWANT) then
               IDO = 100
               return
            end if
         else if (DO2) then
            TRUST(1) = max(BRIGHT, TRUST(1))
            TRUST(2) = LARGE
            TRUST(3) = -LARGE
            if (NEIG.ge.NWANT) then
               IDO = 100
               return
            end if
         else if (DOM2) then
            TRUST(4) = min(BLEFT, TRUST(4))
            TRUST(2) = LARGE
            TRUST(3) = -LARGE
            if (NEIG.ge.NWANT) then
               IDO = 100
               return
            end if
         end if
      end if
      SIGOLD = SIGMA
      if (DO2) then
         BLEFT = EIGEN(NWANT)
         BRIGHT = BLEFT
         do 402 I = NWANT+1,NEIGEN
            if (EIGEN(I).gt.BLEFT) then
               if (BRIGHT.gt.BLEFT) then
                  BRIGHT = min(BRIGHT, EIGEN(I))
               else
                  BRIGHT = EIGEN(I)
               end if
            end if
 402     continue
      else if (DOM2) then
         BRIGHT = EIGEN(NWANT)
         BLEFT = BRIGHT
         do 403 I = NWANT+1,NEIGEN
            if (EIGEN(I).lt.BRIGHT) then
               if (BLEFT.lt.BRIGHT) then
                  BLEFT = max(BLEFT, EIGEN(I))
               else
                  BLEFT = EIGEN(I)
               end if
            end if
 403     continue
         BLEFT = EIGEN(NWANT+1)
      else if (DO4) then
         if (LAST+INERT1.le.NWANT .and. ALLOW2) then
            SIGMA = RANGE(1)
            IDO = 4
            return
         end if
         BRIGHT = EXTREM
         BLEFT = BRIGHT
         do 401 I = NWANT+1,NEIGEN
            if (EIGEN(I).lt.BRIGHT) then
               if (BLEFT.lt.BRIGHT) then
                  BLEFT = max(BLEFT, EIGEN(I))
               else
                  BLEFT = EIGEN(I)
               end if
            end if
 401     continue
         if (BRIGHT.eq.BLEFT) then
            SIGMA = EXTREM - max(ABSDIS, RELDIS*abs(SIGOLD-EXTREM))
            IDO = 4
            return
         end if
      else if (DOM4) then
         if (LAST+(N-INERT1).le.NWANT .and. ALLOW2) then
            SIGMA = RANGE(1)
            IDO = 4
            return
         end if
         if (WHICH.eq.5 .and. EXTREM.ge.RANGE(2)) then
            SIGMA = RANGE(2)
            IDO = 4
            return
         end if
         BLEFT = EXTREM
         BRIGHT = BLEFT
         do 405 I = NWANT+1,NEIGEN
            if (EIGEN(I).gt.BLEFT) then
               if (BRIGHT.gt.BLEFT) then
                  BRIGHT = min(BRIGHT, EIGEN(I))
               else
                  BRIGHT = EIGEN(I)
               end if
            end if
 405     continue
         if (WHICH.ne.5) then
            if (BRIGHT.eq.BLEFT) then
               SIGMA = EXTREM + max(ABSDIS, RELDIS*abs(SIGOLD-EXTREM))
               IDO = 4
               return
            end if
         else
            if (BLEFT.gt.RANGE(2).or.BLEFT.lt.RANGE(1) .or.
     &          BRIGHT.gt.RANGE(2).or.BRIGHT.lt.RANGE(1)) then
               SIGMA = RANGE(2)
               IDO = 4
               return
            end if
         end if
      end if
      SIGMA = AHALF * (BLEFT+BRIGHT)
      if (SIGOLD.ge.BLEFT .and. SIGOLD.le.BRIGHT) then
         SIGMA = SIGOLD
         IDO = 100
         return
      end if
      IDO = 4
      else if (IDO.eq.4) then
         if (.not.ALLOW2) then
            DO2  = .false.
            DOM2 = .false.
            DO4  = WHICH.eq.4 .or. WHICH.eq.-2
            DOM4 = WHICH.eq.-4 .or. WHICH.eq.5 .or. WHICH.eq.2
         else
            DO2 = WHICH.eq.2
            DOM2 = WHICH.eq.-2
            if (WHICH.eq.4) then
               DO2 = SIGMA.eq.RANGE(1)
            else if (WHICH.eq.-4 .or. WHICH.eq.5) then
               DOM2 = SIGMA.eq.RANGE(1)
            end if
         end if
         LAST = 0
         NBETWE = 0
         if (DO2) then
            do 120 I = 1,NCONV
               if (EIGEN(I).le.SIGMA) NBETWE = NBETWE + 1
 120        continue
            NEIG = NEINEG
         else if (DOM2) then
            do 125 I = 1,NCONV
               if (EIGEN(I).ge.SIGMA) NBETWE = NBETWE + 1
 125        continue
            NEIG = N-NEINEG
         else if (DO4) then
            BRIGHT = POLE1
            BLEFT = SIGMA
            if (PURIF5 .and. BRIGHT*BLEFT.lt.ZERO) then
               IDO = 100
               if (SIGMA.ne.SIGOLD) IDO = 6
               SIGMA = SIGOLD
               OK = .true.
               return
            end if
            do 140 I = 1,NCONV
               if (EIGEN(I).le.RANGE(1) .and. EIGEN(I).gt.BLEFT)
     &            LAST = I
               if (EIGEN(I).le.BRIGHT .and. EIGEN(I).ge.BLEFT)
     &            NBETWE = NBETWE + 1
 140        continue
            NEIG = abs(INERT1-NEINEG)
         else if (DOM4) then
            BLEFT = POLE1
            BRIGHT = SIGMA
            if (PURIF5 .and. BRIGHT*BLEFT.lt.ZERO) then
               IDO = 100
               if (SIGMA.ne.SIGOLD) IDO = 6
               SIGMA = SIGOLD
               OK = .true.
               return
            end if
            do 150 I = 1,NCONV
               if (EIGEN(I).ge.RANGE(1) .and. EIGEN(I).lt.BRIGHT)
     &            LAST = I
               if (EIGEN(I).le.BRIGHT .and. EIGEN(I).ge.BLEFT)
     &            NBETWE = NBETWE + 1
 150        continue
            NEIG = abs(INERT1-NEINEG)
         end if
         OK = NEIG.le.NBETWE
         if (LU.ge.0) write(LU,9000) NEIG, min(NEIG,NBETWE)
         if (WHICH.eq.5) then
            if (SIGMA.eq.RANGE(2) .and. .not.OK) then
               if (SIGMA.eq.SIGOLD) then
                  IDO = 100
               else
                  SIGMA = SIGOLD
                  IDO = 6
               end if
               go to 1111
            end if
         else if (DO4 .or. DOM4) then
            if (SIGMA.eq.RANGE(1) .and. .not.OK) then
               if (SIGMA.eq.SIGOLD) then
                  IDO = 100
               else
                  SIGMA = SIGOLD
                  IDO = 6
               end if
               go to 1111
            end if
         end if
         if (OK) then
            if (DO2) then
               TRUST(1) = max(SIGMA, TRUST(1))
               TRUST(2) = LARGE
               TRUST(3) = -LARGE
            else if (DOM2) then
               TRUST(4) = min(SIGMA, TRUST(4))
               TRUST(2) = LARGE
               TRUST(3) = -LARGE
            else
               TRUST(2) = min(BLEFT, TRUST(2))
               TRUST(3) = max(BRIGHT, TRUST(3))
            end if
         end if
         if (DO2.or.DOM2.or..not.ALLOW2) then
         else if (WHICH.eq.4 .or. WHICH.eq.-4) then
            if (LAST+NEIG.lt.NWANT) then
               OK = .false.
            end if
         else if (WHICH.eq.5) then
            if (LAST+NEIG.lt.NWANT .and. SIGMA.lt.RANGE(2)) then
               OK = .false.
            end if
         else
            if (LAST+NEIG.lt.NWANT) OK = .false.
         end if
 1111    continue
         if (NEITST.ge.0) then
            if (NEITST.eq.NEIG .and. NBETST.eq.NBETWE .and.
     &         SIGTST.eq.SIGMA) OK = .true.
         end if
         NEITST = NEIG
         NBETST = NBETWE
         SIGTST = SIGMA
         if (OK) SIGMA = SIGOLD
         if (SIGMA.eq.SIGOLD .and.NORESC) then
            IDO = 100
         else
            SIGMA = SIGOLD
            IDO = 6
         end if
      else if (IDO.eq.6) then
         IDO = 100
      end if
      return
 9000 format('  ',
     &       'The number of eigenvalues in the Sturm sequence test : ',
     &       I6,/,'  ',
     &       'The number of Ritz values satisfying the test        : ',
     &       I6)
      end
      subroutine EA18SD(FIRST, MODE, ALLOW2, WHICH, RANGE, SIGMA, POLE1,
     &                  POLE2, ISINF, SIGOLD, INERT1, INERT2, LBOUND,
     &                  RBOUND, RITZ,
     &                  NRITZ, NCONV, NWACO, NLANC, NWANT, RESID,
     &                  EIGEN, IWORK, N, TRUST, TRIAL, MAXCND,
     &                  ABSDIS, EPS, SAFMIN, IERR, LU)
      logical ALLOW2, FIRST, ISINF
      integer MODE, WHICH, NRITZ, NCONV, NWACO
      integer NLANC, LU, N, TRIAL, IERR, INERT1, INERT2, NWANT
      integer IWORK(NRITZ)
      double precision SIGMA, POLE1, POLE2, SIGOLD, MAXCND, EPS
      double precision ABSDIS, SAFMIN, RBOUND, LBOUND
      double precision RANGE(2), RITZ(NRITZ), EIGEN(NRITZ), RESID(NRITZ)
      double precision TRUST(4)
      double precision RESFAC
      parameter (RESFAC=1.0D0)
      integer HARMON
      parameter (HARMON=1)
      logical SWAP
      logical DOFRS
      integer DIR, IERR2, NEIG, NBETWE, I, EFIRST, ELAST
      double precision BRIGHT, BLEFT, MINPOL, MAXPOL, LARGE
      double precision DLAMCH
      external DLAMCH
      external EA17RD, EA17SD, KB07AD, KB08AD
      LARGE = DLAMCH('O')
      if (WHICH.gt.0) then
         DIR = 1
      else
         DIR = -1
      end if
      if (DIR.gt.0) then
         NEIG = INERT2
      else
         NEIG = N-INERT2
      end if
      if (MODE.ne.5 .and. NEIG.eq.0) ALLOW2 = .true.
      if (ALLOW2 .and. INERT2.ge.0) then
         if (DIR.gt.0) then
            NEIG = INERT2
         else
            NEIG = N-INERT2
         end if
         NBETWE = 0
         do 10 I = 1,NCONV
            if (DIR*EIGEN(I).lt.DIR*POLE2) NBETWE = NBETWE + 1
 10      continue
         if (NEIG.le.NBETWE) then
            if (DIR.gt.0) then
               TRUST(1) = max(TRUST(1), POLE2)
               LBOUND = max(LBOUND, POLE2)
            else
               TRUST(4) = min(TRUST(4), POLE2)
               RBOUND = min(RBOUND, POLE2)
            end if
         end if
         DOFRS = FIRST .or. NEIG.gt.NBETWE
         if (DIR.gt.0) then
            DOFRS = DOFRS .and. LBOUND.le.-LARGE
         else
            DOFRS = DOFRS .and. RBOUND.ge.LARGE
         end if
      else
         DOFRS = NWACO.eq.0
      end if
      if (DOFRS .or. FIRST) then
         EFIRST = 1
         ELAST = min(NWANT+1,NRITZ)
      end if
      if (FIRST) then
         if (DIR.gt.0) then
            call KB07AD(EIGEN, NRITZ, IWORK)
         else
            call KB08AD(EIGEN, NRITZ, IWORK)
         end if
         call EA17RD(MODE, HARMON, DIR, RANGE, RITZ, NRITZ,
     &               RESID, RESFAC, EIGEN, SIGMA, MAXCND,
     &               SIGOLD, ISINF, SAFMIN, EPS, EFIRST,
     &               ELAST, 1, ABSDIS, TRIAL, IERR2)
         if (IERR2.ne.0) then
            IERR = -16
            return
         end if
         if (SIGMA.ne.SIGOLD) INERT2 = -1-N
      else if (DOFRS) then
         if (DIR.gt.0) then
            call KB07AD(EIGEN, NRITZ, IWORK)
         else
            call KB08AD(EIGEN, NRITZ, IWORK)
         end if
         call EA17RD(MODE, HARMON, DIR, RANGE, RITZ, NRITZ,
     &               RESID, RESFAC, EIGEN, SIGMA, MAXCND,
     &               SIGOLD, ISINF, SAFMIN, EPS, EFIRST,
     &               NRITZ, 1, ABSDIS, TRIAL, IERR2)
         if (SIGMA.lt.LBOUND) SIGMA = LBOUND
         if (SIGMA.gt.RBOUND) SIGMA = RBOUND
         if (SIGMA.ne.SIGOLD) INERT2 = -1-N
      else
         MAXPOL = DIR*LARGE
         MINPOL = SIGOLD
         SWAP = .false.
         if (INERT1.le.N .and. INERT1.ge.-N .and.
     &       INERT2.le.N .and. INERT2.ge.-N) then
            if (.not.ALLOW2) then
               BLEFT = min(POLE1,POLE2)
               BRIGHT = max(POLE1,POLE2)
               NEIG = abs(INERT2-INERT1)
               NBETWE = 0
               do 40 I = 1,NCONV
                  if (EIGEN(I).le.BRIGHT .and.
     &                EIGEN(I).ge.BLEFT) NBETWE = NBETWE + 1
 40            continue
            else if (DIR.eq.1) then
               NEIG = INERT2
               BLEFT = -LARGE
               BRIGHT = POLE2
               if (NEIG.eq.0) LBOUND = max(POLE2, LBOUND)
               if (NEIG.gt.NWANT) RBOUND = min(RBOUND,POLE2)
            else
               NEIG = N - INERT2
               BRIGHT = LARGE
               BLEFT = POLE2
               if (NEIG.eq.0) RBOUND = min(POLE2, RBOUND)
               if (NEIG.gt.NWANT) LBOUND = max(POLE2, LBOUND)
            end if
            if (.not.ALLOW2 .and. NEIG.eq.NBETWE) then
               TRUST(2) = min(TRUST(2), BLEFT)
               TRUST(3) = max(TRUST(3), BRIGHT)
            end if
            if (LU.ge.0) then
               write(LU,9001) NEIG
               write(LU,9000) NBETWE
            end if
            if (NEIG-NBETWE.gt.NLANC) then
               SWAP = .true.
               MINPOL = POLE1
               MAXPOL = POLE2
            else if (NEIG.gt.NBETWE) then
               SIGMA = SIGOLD
               return
            end if
         end if
         if (DIR.gt.0) then
            call KB07AD(EIGEN, NRITZ, IWORK)
         else
            call KB08AD(EIGEN, NRITZ, IWORK)
         end if
         if (SWAP) then
            EFIRST = 1
         else
            EFIRST = max(1,NWACO)
         end if
         ELAST = min(NWANT+1,NRITZ)
         EFIRST = min(ELAST,EFIRST)
         call EA17SD(MODE, HARMON, SIGOLD, SIGMA, ISINF, RITZ, NRITZ,
     &               EIGEN, EFIRST, EFIRST, ELAST, 1,
     &               DIR, MINPOL, ABSDIS, MAXCND, SAFMIN, EPS)
         if (SIGMA.eq.SIGOLD) then
            call EA17SD(MODE, HARMON, SIGOLD, SIGMA, ISINF, RITZ,
     &                  NRITZ, EIGEN, EFIRST, ELAST, 1, -1, DIR,
     &                  MINPOL, ABSDIS, MAXCND, SAFMIN, EPS)
         end if
         if (DIR*SIGMA.gt.DIR*MAXPOL) SIGMA = MAXPOL
         if (SIGMA.ne.SIGOLD .and. SWAP) then
            INERT2 = INERT1
            POLE2 = POLE1
         end if
      end if
 9000 format('  ',
     &    'The number of Ritz values satisfying the test        : ',
     &    I6)
 9001 format('  ',
     &       'The number of eigenvalues in the Sturm sequence test : ',
     &       I6)
      end
      subroutine EA18TD(MODE, WHICH, RANGE, SIGMA, POLE1, POLE2,
     &                  LBOUND, RBOUND, ISINF, SIGOLD, INERT1, INERT2,
     &                  CVGED, RITZ, NRITZ, NCONV, NWACO, NWANT,
     &                  LLANC, EIGEN, IWORK, N, TRUST,
     &                  MAXCND, ABSDIS, EPS, SAFMIN, LU)
      logical CVGED, ISINF
      integer MODE, WHICH, NRITZ, NCONV, INERT1, INERT2, N
      integer LLANC, LU, NWANT
      integer IWORK(NRITZ)
      double precision SIGMA, POLE1, POLE2, LBOUND, RBOUND, SIGOLD
      double precision ABSDIS, SAFMIN, MAXCND, EPS
      double precision RANGE(2), RITZ(NRITZ), EIGEN(NRITZ)
      double precision TRUST(4)
      integer ESTEP
      parameter (ESTEP=1)
      integer HARMON
      parameter (HARMON=0)
      logical SWAP
      integer EFIRST, ELAST, DIR, NEIG, NBETWE, I, LAST, NLOCK
      integer EREF, NWACO
      double precision BRIGHT, BLEFT, MINPOL, MAXPOL, LARGE
      double precision SIG, ONEEPS, GROW, MAXDST
      double precision DLAMCH
      external EA17OD, DLAMCH
      external EA17SD, KB07AD, KB08AD
      double precision ONE
      parameter (ONE=1.0D0)
      LARGE = DLAMCH('O')
      ONEEPS = ONE + 100*EPS
      DIR = 1
      if (WHICH.eq.4) DIR = -1
      MAXPOL = DIR*max(DIR*RBOUND, DIR*LBOUND)
      MINPOL = SIGOLD
      SWAP = .false.
      NEIG = -1
      if (INERT1.le.N .and. INERT2.le.N .and.
     &    INERT1.gt.-N .and. INERT2.gt.-N) then
         BRIGHT = max(POLE1,POLE2)
         BLEFT = min(POLE1,POLE2)
         NEIG = abs(INERT2-INERT1)
         NBETWE = 0
         do 40 I = 1,NCONV
            if (EIGEN(I).ge.BLEFT .and. EIGEN(I).le.BRIGHT) then
               NBETWE = NBETWE + 1
            end if
 40      continue
         if (NEIG.eq.NBETWE) then
            TRUST(2) = min(TRUST(2),min(POLE1,POLE2))
            TRUST(3) = max(TRUST(3),max(POLE1,POLE2))
         end if
         if (LU.ge.0) write(LU,9000) NEIG, NBETWE
         if (NEIG.eq.NBETWE) then
            if (WHICH.eq.4 .and. INERT2.eq.0) then
               CVGED = .true.
            else if (WHICH.eq.-4 .and. INERT2.eq.N) then
               CVGED = .true.
            else if (WHICH.eq.5) then
               if (POLE2.ge.RANGE(2)) CVGED = .true.
            end if
         end if
         if (NEIG-NBETWE.gt.LLANC .or. NEIG.ne.NBETWE) then
            SWAP = .true.
            MINPOL = POLE1
            MAXPOL = POLE2
         end if
      else
         MINPOL = RANGE(1)
         MAXPOL = DIR*LARGE
      end if
      if (DIR.eq.1) then
         call KB07AD(EIGEN, NRITZ, IWORK)
      else
         call KB08AD(EIGEN, NRITZ, IWORK)
      end if
      EREF = NRITZ
      LAST = 0
      NLOCK = 0
      NWACO = NRITZ
      do 10 I = 1,NRITZ
         if (DIR*EIGEN(I).lt.DIR*RANGE(1)) then
            LAST = I
            NLOCK = I
            go to 10
         end if
         if (DIR*EIGEN(I).le.DIR*POLE1) NLOCK = I
         if (IWORK(I).gt.NCONV) then
            EREF = I
            go to 11
         end if
         NWACO = I
 10   continue
 11   continue
      NLOCK = NLOCK - LAST
      NWACO = NWACO - LAST
      if (NEIG.ge.0) then
      if (NLOCK+NEIG.ge.NWANT) then
         if (DIR.gt.0) then
            RBOUND = min(RBOUND, POLE2)
         else
            LBOUND = max(LBOUND, POLE2)
         end if
      end if
      end if
      MAXPOL = DIR*min(DIR*MAXPOL,max(DIR*RBOUND, DIR*LBOUND))
      if (SWAP) then
         EFIRST = LAST+1
      else
         EFIRST = LAST+max(1,NWACO)
      end if
      ELAST = min(LAST+NWANT+1,NRITZ)
      EFIRST = min(EFIRST,ELAST)
      call EA17SD(MODE, HARMON, SIGOLD, SIGMA, ISINF, RITZ, NRITZ,
     &            EIGEN, EREF, EFIRST, ELAST, ESTEP,
     &            DIR, MINPOL, ABSDIS, MAXCND, SAFMIN, EPS)
      if (DIR*SIGMA.gt.DIR*MAXPOL) SIGMA = MAXPOL
      if (WHICH.eq.5) then
         if (SIGMA.gt.RANGE(2)) then
            call EA17OD(MODE, HARMON, NRITZ, RITZ, SIGOLD,
     &                  RANGE(2), ISINF, GROW, MAXDST, SAFMIN, LARGE)
            if (GROW.le.MAXCND) then
               SIGMA = RANGE(2)
            else
               SIG = RANGE(2)+(DIR/MAXCND)*(ONEEPS)
               call EA17OD(MODE, HARMON, NRITZ, RITZ, SIGOLD,
     &                     SIG, ISINF, GROW, MAXDST, SAFMIN, LARGE)
               if (GROW.le.MAXCND .and. MAXDST.le.ABSDIS) then
                  SIG = SIGMA
               end if
            end if
            if (SIGMA.gt.RANGE(2) .and. SIGOLD.ge.RANGE(2))
     &         SIGMA = SIGOLD
         end if
      end if
      if (SIGMA.ne.SIGOLD .and. SWAP) then
         INERT2 = INERT1
         POLE2 = POLE1
      end if
 9000 format('  ',
     &       'The number of eigenvalues in the Sturm sequence test : ',
     &       I6,/,'  ',
     &       'The number of Ritz values satisfying the test        : ',
     &       I6)
      end
      subroutine EA18WD(IDO, MODE, IPOS, NRITZ, N, BLK, V, LDV, BV,
     &                  LDBV, RITZ, RESID, VR, SAFMIN, LARGE, JUMP,
     &                  IVEC)
      integer MODE, IDO, NRITZ, N, BLK, LDV, LDBV, VR, JUMP, IVEC
      integer IPOS(4)
      double precision V(LDV,*), RITZ(NRITZ), RESID(NRITZ), BV(LDBV,BLK)
      double precision SAFMIN, LARGE
      double precision RATIO
      integer I, INFO
      double precision DNRM2
      external DNRM2, DAXPY, EA18ND
      if (IDO.eq.0) then
         IVEC = 1
      else if (JUMP.eq.8001) then
         go to 8001
      elseif (JUMP.eq.8002) then
         go to 8002
      elseif (JUMP.eq.8003) then
         go to 8003
      elseif (JUMP.eq.8004) then
         go to 8004
      end if
 1    continue
      if (IVEC.gt.NRITZ) then
         IDO = 100
         return
      end if
 8001 continue
      if (MODE.gt.2) then
         IPOS(1) = IVEC
         IPOS(2) = min(NRITZ,IVEC+BLK-1)
            IPOS(3) = 1
         IPOS(4) = 1
         IDO = 2
         JUMP = 8002
         return
      end if
 8002 continue
      IPOS(1) = IVEC
      IPOS(2) = min(NRITZ,IVEC+BLK-1)
      IPOS(3) = VR
      IPOS(4) = VR + IPOS(2) - IPOS(1)
      IDO = 1
      JUMP = 8003
      return
 8003 continue
      do 10 I = 0,IPOS(2)-IPOS(1)
         call DAXPY(N, -RITZ(IVEC+I), V(1,IVEC+I),1, V(1,VR+I),1)
 10   continue
      if (MODE.gt.2) then
         JUMP = 8004
         IPOS(2) = VR + IPOS(2)-IPOS(1)
         IPOS(1) = VR
         IPOS(3) = 1
         IPOS(4) = IPOS(2)-IPOS(1) + 1
         IDO = 2
         return
      end if
 8004 continue
      if (MODE.le.2) then
         do 20 I = 0,IPOS(2)-IPOS(1)
            RESID(IVEC+I) = DNRM2(N, V(1,VR+I),1)
 20      continue
      else
         do 30 I = 0,IPOS(2)-IPOS(1)
            call EA18ND(N, BV(1,I+1), V(1,VR+I), RESID(IVEC+I),
     &                  RATIO, SAFMIN, INFO)
            if (INFO.ne.0) RESID(IVEC+I) = LARGE
            if (INFO.eq.1) RESID(IVEC+I) = 0.0D0
 30      continue
      end if
      do 40 I = IVEC,IVEC-IPOS(1)+IPOS(2)
         if (abs(RITZ(I)).ge.SAFMIN*RESID(I)) then
            RESID(I) = RESID(I) / abs(RITZ(I))
         else
            RESID(I) = LARGE
         end if
 40   continue
      IVEC = IVEC + 1
      go to 1
      end
      subroutine EA18XD(MODE, T, LDT, LCT, OFFSET,
     &                 NKEEP, NCOOLD, NCONV, BLK, SIGMA, RITZ,
     &                 Z, LDZ, LCZ, NZ, L, LDL, LCL, F, LDF, LCF,
     &                 TAU, LTAU, WORK, LWORK, SAFMIN, IERR)
      integer MODE
      integer LDT, LCT, OFFSET, NKEEP, NCOOLD, NCONV, BLK, LDZ, LCZ
      integer NZ, LDL, LCL
      integer LDF, LCF, LWORK, LTAU, IERR
      double precision SIGMA, SAFMIN
      double precision T(LDT,LCT), RITZ(NKEEP)
      double precision Z(LDZ,LCZ), L(LDL,LCL)
      double precision F(LDF,LCF)
      double precision TAU(LTAU), WORK(LWORK)
      double precision ONE, ZERO
      parameter (ONE=1.0D0, ZERO=0.0D0)
      logical SHIINV
      double precision ALPHA, BETA, GAMMA, DELTA
      external EA17VD, EA17UD, DCOPY
      if (MODE.eq.2 .or. MODE.eq.4) then
         SHIINV = .true.
      else if (MODE.eq.5) then
         SHIINV = .true.
      end if
      if (SHIINV) then
         if (MODE.eq.5) then
            ALPHA = -SIGMA
            BETA = ONE
            GAMMA = SIGMA
            DELTA = ONE
         else
            ALPHA = -SIGMA
            BETA = ONE
            GAMMA = ONE
            DELTA = ZERO
         end if
      end if
      IERR = 0
      if (NCONV.lt.NKEEP) then
         call EA17VD(T(1,NCONV-OFFSET+1), LDT, LCT-NCONV+OFFSET,
     &               NKEEP-NCONV, BLK, Z(1,NCONV-NCOOLD+1), LDZ,
     &               LCZ-(NCONV-NCOOLD), NZ, ALPHA, BETA, GAMMA,
     &               DELTA, L, LDL, LCL, F, LDF, LCF, .true.,
     &               TAU, LTAU, WORK, LWORK, SAFMIN, IERR)
         if (IERR.ne.0) return
      end if
      if (SHIINV) then
         call EA17UD(RITZ, 1, NCONV, SIGMA, SIGMA, MODE, .true.,
     &               SAFMIN)
      end if
      if (NCONV.gt.OFFSET)
     &   call DCOPY(NCONV-OFFSET, RITZ(OFFSET+1),1, T,LDT)
      return
      end
      subroutine EA18YD(MODE, WHICH, RANGE, HARMON, T, LDT, LCT, OFFSET,
     &                 NLANC, NCOOLD, NCONV, BLK, SIGMA, RITZ, NRITZ,
     &                 PIPE, RESID, Z, LDZ, LCZ, NZ, L, LDL, LCL, F,
     &                 LDF, LCF, TAU, LTAU, WORK, LWORK, SAFMIN, IERR)
      integer MODE, WHICH, HARMON
      integer LDT, LCT, NLANC, NCOOLD, NCONV, BLK, LDZ, LCZ, NZ
      integer LDF, LCF, OFFSET, NRITZ, LDL, LCL, LWORK, LTAU, IERR
      double precision SIGMA, SAFMIN
      double precision RANGE(2), T(LDT,LCT), RITZ(NRITZ), RESID(NRITZ)
      double precision PIPE(BLK,NLANC), Z(LDZ,LCZ), L(LDL,LCL)
      double precision F(LDF,LCF)
      double precision TAU(LTAU), WORK(LWORK)
      double precision ONE, ZERO
      parameter (ONE=1.0D0, ZERO=0.0D0)
      logical SHIINV
      integer I, J, JT
      integer IE, K
      double precision ALPHA, BETA, GAMMA, DELTA
      external DSBEV, EA16RD, EA17VD, EA17PD, EA16SD
      double precision DNRM2
      external DNRM2
      IERR = 0
      HARMON = 0
      if (MODE.eq.2 .or. MODE.eq.4) then
         if (WHICH.eq.2 .or. WHICH.eq.-2) then
            HARMON = 1
         end if
         SHIINV = .true.
      else if (MODE.eq.5) then
         SHIINV = .true.
         if (WHICH.eq.2 .or. WHICH.eq.-2) then
            HARMON = 1
         end if
      end if
      if (HARMON.eq.0) return
      call EA16RD(NLANC, BLK, NCONV-NCOOLD, NLANC,
     &            T(1,NCOOLD-OFFSET+1), LDT, RITZ(NCOOLD+1), Z,
     &            LDZ, LCZ, NZ, PIPE, BLK, L, LDL, LCL)
      if (SHIINV) then
         if (MODE.eq.5) then
            ALPHA = ONE
            BETA = -ONE
            GAMMA = -SIGMA
            DELTA = SIGMA
         else
            ALPHA = ZERO
            BETA = ONE
            GAMMA = ONE
            DELTA = SIGMA
         end if
      end if
      if (NCONV.lt.NCOOLD+NLANC) then
          call EA17VD(T(1,NCONV-OFFSET+1), LDT, LCT-NCONV+OFFSET,
     &                NLANC+NCOOLD-NCONV, BLK, Z(1,NCONV-NCOOLD+1),
     &                LDZ, LCZ-(NCONV-NCOOLD), NZ, ALPHA, BETA,
     &                GAMMA, DELTA, L, LDL, LCL, F, LDF, LCF,
     &                .true., TAU, LTAU, WORK, LWORK, SAFMIN, IERR)
         if (IERR.ne.0) return
      end if
      do 11 I = OFFSET+1,NCONV
         WORK(I-OFFSET) = RITZ(I)
 11   continue
      if (SHIINV) then
         call EA17PD(MODE, 0, NCONV-OFFSET, WORK, SIGMA, RANGE,
     &               SAFMIN)
      end if
      do 12 I = 1,NCONV-OFFSET
         T(1,I) = WORK(I)
 12   continue
      do 20 I = 1,BLK
         do 200 J = 1,BLK
            PIPE(J,I) = T(J+1,NLANC+NCOOLD-OFFSET-BLK+I)
 200     continue
 20   continue
      if (NCONV.lt.NLANC+NCOOLD) then
         call DSBEV('V', 'L', NLANC+NCOOLD-NCONV, BLK,
     &              T(1,NCONV-OFFSET+1), LDT,
     &              RITZ(NCONV+1), F, LDF, WORK, IERR)
         if (IERR.ne.0) then
            return
         end if
      end if
      if (SHIINV) then
         call EA17PD(MODE, 0, NCONV, RITZ, SIGMA, RANGE, SAFMIN)
      end if
      call EA16SD(Z(1,NCONV-NCOOLD+1), LDZ, LCZ-(NCONV-NCOOLD),
     &            NZ, NLANC+NCOOLD-NCONV, F, LDF,
     &            LCF, NLANC+NCOOLD-NCONV, WORK, 1, NLANC+NCOOLD-NCONV)
      do 30 I = 1,BLK
         do 300 J = 1,BLK
            T(J+1,NLANC+NCOOLD-OFFSET-BLK+I) = PIPE(J,I)
 300     continue
 30   continue
      do 40 I = 1,NLANC
         do 400 J = 1,BLK
            PIPE(J,I) = 0.
 400     continue
 40   continue
      do 50 I = max(NCONV+BLK-NLANC-NCOOLD+1,1),BLK
         IE = NLANC+NCOOLD-BLK+I
         JT = BLK-I+1
         do 500 J = 1,I
            JT = JT + 1
            do 5000 K = NCONV+1,NCOOLD+NLANC
               PIPE(J,K-NCOOLD) = PIPE(J,K-NCOOLD)
     &                       + T(JT,IE-OFFSET) * F(IE-NCONV,K-NCONV)
 5000       continue
 500     continue
 50   continue
      do 60 I = NCONV-NCOOLD+1,NLANC
         RESID(I+NCOOLD) = DNRM2(BLK, PIPE(1,I),1)
 60   continue
      return
      end
      integer function EA18ZD(NEINEG, MODE, SIGMA, N)
      integer NEINEG, MODE, N
      double precision SIGMA
      double precision ZERO
      parameter (ZERO=0.0D0)
      if (NEINEG.lt.0 .or. NEINEG.gt.N) then
         EA18ZD = -N-1
      else
         if (SIGMA.lt.ZERO.and.MODE.eq.5) then
            EA18ZD = -NEINEG
         else
            EA18ZD = NEINEG
         end if
      end if
      return
      end

