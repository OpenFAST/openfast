 MODULE AeroGenSubs


   USE NWTC_LIBRARY

! SUBROUTINE AllocArrays( Arg )
! SUBROUTINE ElemOpen( ElemFile )
! SUBROUTINE ElemOut( )

 CONTAINS
 ! ****************************************************
   SUBROUTINE AllocArrays ( Arg )
 !  Allocates space to the phenomenal number of arrays
 !  we use in this program
 ! ****************************************************


USE             Airfoil
USE             Bedoes
USE             Blade
USE             DynInflow
USE             Element
USE             ElOutParams
USE             Switch


IMPLICIT        NONE


   ! Passed Variables:

CHARACTER(*) :: Arg


   ! Local Variables:

INTEGER(4)   :: Sttus

Sttus = 0.0

IF (Arg(1:7) == 'Element') THEN

   IF (.NOT. ALLOCATED(ElPrList)) ALLOCATE ( ElPrList(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for ElPrList array.' )
   ElPrList ( : ) = 0

   IF (.NOT. ALLOCATED(WndElPrList)) ALLOCATE ( WndElPrList(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for WndElPrList array.' )
   WndElPrList ( : ) = 0

   IF (.NOT. ALLOCATED(A)) ALLOCATE ( A(NELM,NB) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for A array.' )
   A ( :, : ) = 0.0

   IF (.NOT. ALLOCATED(AP)) ALLOCATE ( AP(NELM,NB) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AP array.' )

! Beddoes arrays
   IF (DSTALL) THEN

      IF (.NOT. ALLOCATED(ADOT)) ALLOCATE ( ADOT(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ADOT array.' )
      ADOT ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(ADOT1)) ALLOCATE ( ADOT1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ADOT1 array.' )

      IF (.NOT. ALLOCATED(AFE)) ALLOCATE ( AFE(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AFE array.' )
      AFE(:,:) = 0.0


      IF (.NOT. ALLOCATED(AFE1)) ALLOCATE ( AFE1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AFE1 array.' )

      IF (.NOT. ALLOCATED(ANE)) ALLOCATE ( ANE(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ANE array.' )
      ANE ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(ANE1)) ALLOCATE ( ANE1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ANE1 array.' )

      IF (.NOT. ALLOCATED(AOD)) ALLOCATE ( AOD(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AOD array.' )

      IF (.NOT. ALLOCATED(AOL)) ALLOCATE ( AOL(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AOL array.' )

      IF (.NOT. ALLOCATED(BEDSEP)) ALLOCATE ( BEDSEP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for BEDSEP array.' )
      BEDSEP(:,:) = .FALSE.

      IF (.NOT. ALLOCATED(CDO)) ALLOCATE ( CDO(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CDO array.' )

      IF (.NOT. ALLOCATED(CNA)) ALLOCATE ( CNA(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNA array.' )

      IF (.NOT. ALLOCATED(CNP)) ALLOCATE ( CNP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNP array.' )
      CNP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(CNP1)) ALLOCATE ( CNP1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNP1 array.' )

      IF (.NOT. ALLOCATED(CNPD)) ALLOCATE ( CNPD(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPD array.' )
      CNPD ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(CNPD1)) ALLOCATE ( CNPD1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPD1 array.' )

      IF (.NOT. ALLOCATED(CNPOT)) ALLOCATE ( CNPOT(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPOT array.' )
      CNPOT ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(CNPOT1)) ALLOCATE ( CNPOT1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPOT1 array.' )

      IF (.NOT. ALLOCATED(CNS)) ALLOCATE ( CNS(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNS array.' )

      IF (.NOT. ALLOCATED(CNSL)) ALLOCATE ( CNSL(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNSL array.' )

      IF (.NOT. ALLOCATED(CNV)) ALLOCATE ( CNV(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNV array.' )
      CNV ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(CVN)) ALLOCATE ( CVN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CVN array.' )
      CVN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(CVN1)) ALLOCATE ( CVN1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CVN1 array.' )

      IF (.NOT. ALLOCATED(DF)) ALLOCATE ( DF(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DF array.' )
      DF( :, : ) = 0.0

      IF (.NOT. ALLOCATED(DFAFE)) ALLOCATE ( DFAFE(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFAFE array.' )
      DFAFE ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(DFAFE1)) ALLOCATE ( DFAFE1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFAFE1 array.' )

      IF (.NOT. ALLOCATED(DFC)) ALLOCATE ( DFC(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFC array.' )
      DFC ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(DN)) ALLOCATE ( DN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DN array.' )
      DN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(DPP)) ALLOCATE ( DPP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DP array.' )
      DPP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(DQ)) ALLOCATE ( DQ(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DQ array.' )
      DQ ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(DQP)) ALLOCATE ( DQP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DQP array.' )
      DQP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(DQP1)) ALLOCATE ( DQP1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DQP1 array.' )

      IF (.NOT. ALLOCATED(FSP)) ALLOCATE ( FSP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSP array.' )
      FSP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(FSP1)) ALLOCATE ( FSP1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSP1 array.' )

      IF (.NOT. ALLOCATED(FSPC)) ALLOCATE ( FSPC(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSPC array.' )
      FSPC ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(FSPC1)) ALLOCATE ( FSPC1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSPC1 array.' )

      IF (.NOT. ALLOCATED(OLDCNV)) ALLOCATE ( OLDCNV(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDCNV array.' )

      IF (.NOT. ALLOCATED(OLDDF)) ALLOCATE ( OLDDF(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDF array.' )

      IF (.NOT. ALLOCATED(OLDDFC)) ALLOCATE ( OLDDFC(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDFC array.' )

      IF (.NOT. ALLOCATED(OLDDN)) ALLOCATE ( OLDDN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDN array.' )

      IF (.NOT. ALLOCATED(OLDDPP)) ALLOCATE ( OLDDPP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDP array.' )

      IF (.NOT. ALLOCATED(OLDDQ)) ALLOCATE ( OLDDQ(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDQ array.' )

      IF (.NOT. ALLOCATED(OLDTAU)) ALLOCATE ( OLDTAU(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDTAU array.' )

      IF (.NOT. ALLOCATED(OLDXN)) ALLOCATE ( OLDXN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDXN array.' )

      IF (.NOT. ALLOCATED(OLDYN)) ALLOCATE ( OLDYN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDYN array.' )

      IF (.NOT. ALLOCATED(OLDSEP)) ALLOCATE ( OLDSEP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDSEP array.' )

      IF (.NOT. ALLOCATED(QX)) ALLOCATE ( QX(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for QX array.' )
      QX ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(QX1)) ALLOCATE ( QX1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for QX1 array.' )

      IF (.NOT. ALLOCATED(TAU)) ALLOCATE ( TAU(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for TAU array.' )
      TAU(:,:) = 0.0

      IF (.NOT. ALLOCATED(XN)) ALLOCATE ( XN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for XN array.' )
      XN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(YN)) ALLOCATE ( YN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for YN array.' )
      YN ( :, : ) = 0.0

   ENDIF ! Beddoes arrays

   IF (.NOT. ALLOCATED(C)) ALLOCATE ( C(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for C array.' )

   IF (.NOT. ALLOCATED(DR)) ALLOCATE ( DR(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DR array.' )

   IF (.NOT. ALLOCATED(RELM)) ALLOCATE ( RELM(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for RELM array.' )

   IF (.NOT. ALLOCATED(TWIST)) ALLOCATE ( TWIST(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for TWIST array.' )

   IF (.NOT. ALLOCATED(TLCNST)) ALLOCATE ( TLCNST(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for TLCNST array.' )
   TLCNST = 99.0

   IF (.NOT. ALLOCATED(HLCNST)) ALLOCATE ( HLCNST(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for HLCNST array.' )
   HLCNST = 99.0

   IF (.NOT. ALLOCATED(NFOIL)) ALLOCATE ( NFOIL(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for NFOIL array.' )

   IF (.NOT. ALLOCATED(NLIFT)) ALLOCATE ( NLIFT(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for NLIFT array.' )

   IF (.NOT. ALLOCATED(NTables)) ALLOCATE ( NTables(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for NTables array.' )

   IF (DYNINFL .OR. DYNINIT) THEN
      IF (.NOT. ALLOCATED(RMC_SAVE)) ALLOCATE ( RMC_SAVE ( NB, NELM, MAXINFL ) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for RMC_SAVE array.' )
      RMC_SAVE = 0.0

      IF (.NOT. ALLOCATED(RMS_SAVE)) ALLOCATE ( RMS_SAVE ( NB, NELM, MAXINFL ) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for RMS_SAVE array.' )
      RMS_SAVE = 0.0
   ENDIF


ELSEIF (Arg(1:7) == 'ElPrint') THEN

   IF ( NumElOut > 0 ) THEN
      IF (.NOT. ALLOCATED(AAA)) ALLOCATE ( AAA(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AAA array.' )

      IF (.NOT. ALLOCATED(AAP)) ALLOCATE ( AAP(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AAP array.' )

      IF (.NOT. ALLOCATED(ALF)) ALLOCATE ( ALF(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ALF array.' )

      IF (.NOT. ALLOCATED(CDD)) ALLOCATE ( CDD(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CDD array.' )

      IF (.NOT. ALLOCATED(CLL)) ALLOCATE ( CLL(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CLL array.' )

      IF (.NOT. ALLOCATED(CMM)) ALLOCATE ( CMM(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CMM array.' )

      IF (.NOT. ALLOCATED(CNN)) ALLOCATE ( CNN(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNN array.' )

      IF (.NOT. ALLOCATED(CTT)) ALLOCATE ( CTT(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CTT array.' )

      IF (.NOT. ALLOCATED(DFNSAV)) ALLOCATE ( DFNSAV(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFNSAV array.' )

      IF (.NOT. ALLOCATED(DFTSAV)) ALLOCATE ( DFTSAV(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFTSAV array.' )

      IF (.NOT. ALLOCATED(DynPres)) ALLOCATE ( DynPres(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DynPres array.' )

      IF (.NOT. ALLOCATED(PMM)) ALLOCATE ( PMM(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for PMM array.' )

      IF (.NOT. ALLOCATED(PITSAV)) ALLOCATE ( PITSAV(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for PITSAV array.' )

      IF (.NOT. ALLOCATED(ReyNum)) ALLOCATE ( ReyNum(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ReyNum array.' )

      IF (.NOT. ALLOCATED(ElPrNum)) ALLOCATE ( ElPrNum(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for ElPrNum array.' )
      ElPrNum ( : ) = 0

   END IF

   IF ( NumWndElOut > 0 ) THEN

      IF (.NOT. ALLOCATED(WndElPrNum)) ALLOCATE ( WndElPrNum(NumWndElOut) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for WndElPrNum array.' )
      WndElPrNum ( : ) = 0

      IF (.NOT. ALLOCATED(SaveVX)) ALLOCATE ( SaveVX(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVX array.' )

      IF (.NOT. ALLOCATED(SaveVY)) ALLOCATE ( SaveVY(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVY array.' )

      IF (.NOT. ALLOCATED(SaveVZ)) ALLOCATE ( SaveVZ(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVZ array.' )

   END IF

ELSEIF (Arg(1:8) == 'Aerodata') THEN

   IF (.NOT. ALLOCATED(AL)) ALLOCATE ( AL(NumFoil,NumCL) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AL array.' )

   IF (.NOT. ALLOCATED(CD)) ALLOCATE ( CD(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CD array.' )

   IF (.NOT. ALLOCATED(CL)) ALLOCATE ( CL(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CL array.' )

   IF (.NOT. ALLOCATED(CM)) ALLOCATE ( CM(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CM array.' )

   IF (.NOT. ALLOCATED(MulTabMet)) ALLOCATE ( MulTabMet(NumFoil,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for MulTabMet array.' )

   IF (DSTALL) THEN

      IF (.NOT. ALLOCATED(FTB)) ALLOCATE ( FTB(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FTB array.' )

      IF (.NOT. ALLOCATED(FTBC)) ALLOCATE ( FTBC(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FTBC array.' )

   ENDIF ! Beddoes arrays

ENDIF



RETURN
END SUBROUTINE AllocArrays

 ! *****************************************************
   SUBROUTINE ElemOpen (ElemFile)
 !  This subroutine opens the element output file and writes
 !   column headings separated by tab characters
 !  ElemFile = file name
 ! *****************************************************

USE               Blade    ! We need the number of blades for the wind output

USE               ElOutParams
USE               Element
USE               Switch


IMPLICIT          NONE


   ! Passed Variables:

CHARACTER(  *) :: ElemFile


   ! Local Variables:

INTEGER(4)     :: JE

INTEGER(4)     :: JB       ! Counter for number of blades

CHARACTER(  2) :: Dst_Unit
CHARACTER(  2) :: Frc_Unit
CHARACTER(140) :: Frmt
CHARACTER(  3) :: Prs_Unit


IF (NumWndElOut > 0) THEN
   CALL OpenFOutFile (UnWndOut, TRIM(ElemFile)//'.wind')
   WRITE (UnWndOut,"( 'This file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' )")  &
        TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()
ENDIF


 ! Open the Element Print file if requested
IF (ELEMPRN) THEN
   CALL OpenFOutFile (UnElem, TRIM(ElemFile))
   WRITE (UnElem,"( 'This file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' )")  &
        TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()
ELSE
   RETURN
ENDIF

 ! Set the units labels
IF (SIUNIT) THEN
  Dst_Unit = 'm'
  Frc_Unit = 'N'
  Prs_Unit = 'Pa'
ELSE
  Dst_Unit = 'ft'
  Frc_Unit = 'lb'
  Prs_Unit = 'psf'
ENDIF


Frmt = '( A4, 3(A1,A2,I2.2),    (: A1, A, I2.2 ) )'

IF ( PMOMENT ) THEN
   WRITE(Frmt(22:24), '(I3)') 14*NumElOut
   WRITE(UnElem, Frmt) 'Time',                    &
               TAB,    'VX',       NELM,          &
               TAB,    'VY',       NELM,          &
               TAB,    'VZ',       NELM,          &
             ( TAB,    'Alpha',    ElPrNum(JE),  &
               TAB,    'DynPres',  ElPrNum(JE),  &
               TAB,    'CLift',    ElPrNum(JE),  &
               TAB,    'CDrag',    ElPrNum(JE),  &
               TAB,    'CNorm',    ElPrNum(JE),  &
               TAB,    'CTang',    ElPrNum(JE),  &
               TAB,    'CMomt',    ElPrNum(JE),  &
               TAB,    'Pitch',    ElPrNum(JE),  &
               TAB,    'AxInd',    ElPrNum(JE),  &
               TAB,    'TanInd',   ElPrNum(JE),  &
               TAB,    'ForcN',    ElPrNum(JE),  &
               TAB,    'ForcT',    ElPrNum(JE),  &
               TAB,    'Pmomt',    ElPrNum(JE),  &
               TAB,    'ReNum',    ElPrNum(JE),  &
                         JE = 1, NumElOut )

   Frmt = '( A5, 3(A1,A8),    (: A1, A ) )'
   WRITE(Frmt(17:19), '(I3)') 14*NumElOut
   WRITE(UnElem, Frmt) '(sec)',                       &
               TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
               TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
               TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
             ( TAB,    '(deg)',                       &
               TAB,    '('//TRIM(Prs_Unit)//')',      &
               TAB,    '(-)',                         &
               TAB,    '(-)',                         &
               TAB,    '(-)',                         &
               TAB,    '(-)',                         &
               TAB,    '(-)',                         &
               TAB,    '(deg)',                       &
               TAB,    '(-)',                         &
               TAB,    '(-)',                         &
               TAB,    '('//TRIM(Frc_Unit)//')',      &
               TAB,    '('//TRIM(Frc_Unit)//')',      &
               TAB,    '('//TRIM(Frc_Unit)//'-'//TRIM(Dst_Unit)//')', &
               TAB,    '(x10^6)',                     &
                         JE = 1, NumElOut )

ELSE
   WRITE(Frmt(22:24), '(I3)') 12*NumElOut
   WRITE(UnElem, Frmt) 'Time',                    &
               TAB,    'VX',       NELM,          &
               TAB,    'VY',       NELM,          &
               TAB,    'VZ',       NELM,          &
             ( TAB,    'Alpha',    ElPrNum(JE),  &
               TAB,    'DynPres',  ElPrNum(JE),  &
               TAB,    'CLift',    ElPrNum(JE),  &
               TAB,    'CDrag',    ElPrNum(JE),  &
               TAB,    'CNorm',    ElPrNum(JE),  &
               TAB,    'CTang',    ElPrNum(JE),  &
               TAB,    'Pitch',    ElPrNum(JE),  &
               TAB,    'AxInd',    ElPrNum(JE),  &
               TAB,    'TanInd',   ElPrNum(JE),  &
               TAB,    'ForcN',    ElPrNum(JE),  &
               TAB,    'ForcT',    ElPrNum(JE),  &
               TAB,    'ReNum',    ElPrNum(JE),  &
                         JE = 1, NumElOut )

   Frmt = '( A5, 3(A1,A8),    (: A1, A ) )'
   WRITE(Frmt(17:19), '(I3)') 12*NumElOut
   WRITE(UnElem, Frmt) '(sec)',                       &
               TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
               TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
               TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
             ( TAB,    '(deg)',                       &
               TAB,    '('//TRIM(Prs_Unit)//')',      &
               TAB,    '(-)',                         &
               TAB,    '(-)',                         &
               TAB,    '(-)',                         &
               TAB,    '(-)',                         &
               TAB,    '(deg)',                       &
               TAB,    '(-)',                         &
               TAB,    '(-)',                         &
               TAB,    '('//TRIM(Frc_Unit)//')',      &
               TAB,    '('//TRIM(Frc_Unit)//')',      &
               TAB,    '(x10^6)',                     &
                         JE = 1, NumElOut )
ENDIF

IF ( NumWndElOut > 0 ) THEN
   Frmt = '( A4, XXX(A1,A2,I2.2,"-B",I1.1) )'
   WRITE(Frmt(7:9), '(I3)') 3*NumWndElOut*NB

   WRITE(UnWndOut, Frmt) 'Time',          &
             ( ( TAB, 'VX',  WndElPrNum(JE),  JB,   &
                 TAB, 'VY',  WndElPrNum(JE),  JB,   &
                 TAB, 'VZ',  WndElPrNum(JE),  JB,   &
                   JE = 1, NumWndElOut ) , &
                   JB = 1, NB )

   Frmt = '( A5, XXX(A1,A8) )'
   WRITE(Frmt(7:9), '(I3)') 3*NumWndElOut*NB

   WRITE(UnWndOut, Frmt)   '(sec)',                       &
               ( TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
                 TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
                 TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
                                  JE = 1, NumWndElOut*NB )
ENDIF

RETURN
END SUBROUTINE ElemOpen


 ! *****************************************************
   SUBROUTINE ElemOut
 !  This subroutine writes the element output values
 !   for the desired elements
 ! *****************************************************

USE               Blade       ! Get the number of blades

USE               ElOutParams
USE               AeroTime
USE               Switch


IMPLICIT          NONE


   ! Local Variables:

INTEGER(4)     :: JE

INTEGER(4)     :: JB    ! Counter for number of blades

CHARACTER(30)  :: Frmt



 ! Write the element data if requested
IF (ELEMPRN) THEN

   Frmt = '( F10.3,    ( : A1, ES12.5 ) )'

   IF ( PMOMENT ) THEN
      WRITE(Frmt(10:12), '(I3)') 14*NumElOut + 3
      WRITE(UnElem,Frmt) TIME,      TAB,           &
                     VXSAV,         TAB,           &
                     VYSAV,         TAB,           &
                     VZSAV,                        &
                   ( TAB,   ALF    (JE), &
                     TAB,   DynPres(JE), &
                     TAB,   CLL    (JE), &
                     TAB,   CDD    (JE), &
                     TAB,   CNN    (JE), &
                     TAB,   CTT    (JE), &
                     TAB,   CMM    (JE), &
                     TAB,   PITSAV (JE), &
                     TAB,   AAA    (JE), &
                     TAB,   AAP    (JE), &
                     TAB,   DFNSAV (JE), &
                     TAB,   DFTSAV (JE), &
                     TAB,   PMM    (JE), &
                     TAB,   ReyNum (JE), &
                            JE= 1, NumElOut )


   ELSE
      WRITE(Frmt(10:12), '(I3)') 12*NumElOut + 3
      WRITE(UnElem,Frmt) TIME,      TAB,           &
                     VXSAV,         TAB,           &
                     VYSAV,         TAB,           &
                     VZSAV,                        &
                   ( TAB,   ALF    (JE), &
                     TAB,   DynPres(JE), &
                     TAB,   CLL    (JE), &
                     TAB,   CDD    (JE), &
                     TAB,   CNN    (JE), &
                     TAB,   CTT    (JE), &
                     TAB,   PITSAV (JE), &
                     TAB,   AAA    (JE), &
                     TAB,   AAP    (JE), &
                     TAB,   DFNSAV (JE), &
                     TAB,   DFTSAV (JE), &
                     TAB,   ReyNum (JE), &
                            JE= 1, NumElOut )
   ENDIF ! PMOMENT

IF (NumWndElOut > 0) THEN

   WRITE(Frmt(10:12), '(I3)') 3*NumWndElOut*NB
   WRITE(UnWndOut,Frmt) TIME,                   &
                    ( (  TAB, SaveVX( JE, JB ), &
                         TAB, SaveVY( JE, JB ), &
                         TAB, SaveVZ( JE, JB ), &
                              JE = 1,NumWndElOut ), JB = 1,NB )
ENDIF

ENDIF ! ELEMPRN



RETURN
END SUBROUTINE ElemOut
!=======================================================================

END MODULE AeroGenSubs
