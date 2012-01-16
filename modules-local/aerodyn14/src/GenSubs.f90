 !bjj Start of proposed change v12.70b
 MODULE AeroGenSubs


   USE NWTC_LIBRARY

!bjj rm SUBROUTINE SetProgName( )
! FUNCTION Cross_Product(Vector1, Vector2)
! SUBROUTINE AllocArrays( Arg )
! SUBROUTINE ElemOpen( ElemFile )
! SUBROUTINE ElemOut( )

 CONTAINS
!bjj: this function could be moved to NWTC Library... I am not using it here anymore, either.
! FUNCTION Cross_Product(Vector1, Vector2)
!   ! Cross_Product = Vector1 X Vector2 (resulting in a vector)
!   
!      IMPLICIT NONE
!
!
!      ! Passed variables:
!
!   REAL(ReKi), INTENT(IN )         :: Vector1       (3)
!   REAL(ReKi), INTENT(IN )         :: Vector2       (3)
!
!      ! Function definition
!   REAL(ReKi)                      :: Cross_Product (3)   ! = Vector1 X Vector2 (resulting in a vector)
!
!
!   Cross_Product(1) = Vector1(2)*Vector2(3) - Vector1(3)*Vector2(2)
!   Cross_Product(2) = Vector1(3)*Vector2(1) - Vector1(1)*Vector2(3)
!   Cross_Product(3) = Vector1(1)*Vector2(2) - Vector1(2)*Vector2(1)
!
!
!   RETURN
!END FUNCTION Cross_Product
 !Bjj end of proposed change

!bjj start of proposed change v13.00
!bjj remove this routine:
!rm ! *****************************************************
!rm   SUBROUTINE SetProgName
!rm !  This subroutine writes the dynamic and aerodynamic
!rm !   program names and versions to a string.
!rm ! *****************************************************
!rm
!rm
!rmUSE               Identify
!rm
!rm
!rmIMPLICIT          NONE
!rm
!rm
!rm   ! Local Variables:
!rm
!rm!bjj rm CHARACTER(80)  :: Frmt
!rm!bjj rm CHARACTER(80)  :: MESAGE
!rm
!rm
!rm
!rmAeroProg = 'AeroDyn '
!rm!bjj Start of proposed change v12.70-bjj
!rm!rm AeroVer  = '(12.58, 28-Jun-2005)'
!rmAeroVer  = '(12.70-bjj, 10-Jul-2009)'
!rm!bjj End of proposed change
!rmCreator  = 'Windward Engineering LC'
!rm
!rmSELECT CASE (TRIM(DynProg))
!rm   CASE ('ADAMS','YawDyn','SymDyn','FAST') ! These are expected
!rm   CASE DEFAULT
!rm!bjj Start of proposed change v12.70-bjj
!rm!bjj: is this necessary? I'd like to remove all instances of DynProg...
!rm!rm   MESAGE = 'AeroDyn does not recognize the dynamics program name: '//Trim(DynProg)
!rm!rm   Frmt   = '(A/, ''  Unexpected problems may occur during simulation.'')'
!rm!rm   CALL ErrLog( MESAGE, Frmt, 'SetProgName', 801, 'WARN' )
!rm      CALL ProgWarn( 'AeroDyn does not recognize the dynamics program name: '//Trim(DynProg)// &
!rm                       '. Unexpected problems may occur during simulation.' )
!rm!bjj End of proposed change v12.70-bjj
!rmEND SELECT
!rm
!rmWRITE(Prog,'(A)') TRIM(AeroProg)//TRIM(AeroVer)//' in '// &
!rm                     TRIM(DynProg)//TRIM(DynVer)
!rm
!rm
!rm
!rmRETURN
!rmEND SUBROUTINE SetProgName
!bjj end of proposed change
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
!bjj Start of proposed change v12.70-bjj
!rm USE             InducedVel  !bjj THIS ISN'T USED?
!bjj End of proposed change
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
!bjj start of proposed change v13.00b
      AFE(:,:) = 0.0
!bjj start of proposed change v13.00b


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
!bjj start of proposed change v13.00b
      BEDSEP(:,:) = .FALSE.
!bjj end of proposed change v13.00b

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
!bjj start of proposed change v13.00b
      TAU(:,:) = 0.0
!bjj end of proposed change

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

!BJJ START OF PROPOSED CHANGE
   IF ( NumElOut > 0 ) THEN
!BJJ END OF PROPOSED CHANGE
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

!BJJ START OF PROPOSED CHANGE
   END IF

   IF ( NumWndElOut > 0 ) THEN
!BJJ END OF PROPOSED CHANGE

      IF (.NOT. ALLOCATED(WndElPrNum)) ALLOCATE ( WndElPrNum(NumWndElOut) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for WndElPrNum array.' )
      WndElPrNum ( : ) = 0

      IF (.NOT. ALLOCATED(SaveVX)) ALLOCATE ( SaveVX(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVX array.' )

      IF (.NOT. ALLOCATED(SaveVY)) ALLOCATE ( SaveVY(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVY array.' )

      IF (.NOT. ALLOCATED(SaveVZ)) ALLOCATE ( SaveVZ(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVZ array.' )

!BJJ START OF PROPOSED CHANGE
   END IF
!BJJ END OF PROPOSED CHANGE

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
!bjj rm NWTC_Lib:USE               Constant
USE               Element
USE               Switch


IMPLICIT          NONE


   ! Passed Variables:

CHARACTER(  *) :: ElemFile


   ! Local Variables:

INTEGER(4)     :: JE
!rm not used:INTEGER(4)     :: mode

INTEGER(4)     :: JB       ! Counter for number of blades

CHARACTER(  2) :: Dst_Unit
CHARACTER(  2) :: Frc_Unit
CHARACTER(140) :: Frmt
CHARACTER(  3) :: Prs_Unit


!bjj:   ! YES, the file name should be fixed, but for now
!bjj    ! this works....
IF (NumWndElOut > 0) THEN
!bjj Start of proposed change v12.70b-bjj
!rm      CALL OpenOutFile (UnWndOut, TRIM(ElemFile)//'.wind')
   CALL OpenFOutFile (UnWndOut, TRIM(ElemFile)//'.wind')
   WRITE (UnWndOut,"( 'This file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' )")  &
        TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()
!bjj End of proposed change
ENDIF


 ! Open the Element Print file if requested
IF (ELEMPRN) THEN
!bjj Start of proposed change v12.70b-bjj
!rm      CALL OpenOutFile (UnElem, TRIM(ElemFile))
   CALL OpenFOutFile (UnElem, TRIM(ElemFile))
   WRITE (UnElem,"( 'This file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' )")  &
        TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()
!bjj End of proposed change
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
!bjj rm NWTC_Lib:USE               Constant
USE               Switch
!bjj Start of proposed change v12.70
!rm USE               DynInflow !NOT USED
!bjj End of proposed change


IMPLICIT          NONE


   ! Local Variables:

INTEGER(4)     :: JE
!rm not used: INTEGER(4)     :: mode

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


!bjj start of proposed change
!rm! ***********************************************************
!rm   SUBROUTINE ErrLog ( MESAGE, ErrFrmt, Sub, IDADAMS, ErrType )
!rm!  Logs errors from YawDyn to a file for future access.
!rm
!rm!  D. Laino  11/00 Windward Engineering
!rm! ***********************************************************
!rm
!rmUSE                          AeroTime
!rmUSE                          ErrCount
!rmUSE                          Identify
!rm
!rmUSE                          AD_IOParams, ONLY: UnErrLog
!rm
!rmIMPLICIT                     NONE
!rm
!rm
!rm   ! Passed Variables:
!rm
!rmINTEGER(4)                :: IDADAMS
!rm
!rmCHARACTER( *)             :: ErrFrmt
!rmCHARACTER( *)             :: ErrType
!rmCHARACTER( *)             :: MESAGE
!rmCHARACTER( *)             :: Sub
!rm
!rm
!rm   ! Local Variables:
!rm
!rmINTEGER(4)                :: IERR
!rm
!rmLOGICAL                   :: Exists
!rmLOGICAL,    SAVE          :: FrstPass = .TRUE.
!rmLOGICAL,    SAVE          :: LogFlag  = .TRUE.
!rm
!rmCHARACTER(60)             :: Frmt
!rmCHARACTER(80)             :: Text
!rm
!rm
!rm   ! Global Functions:
!rm
!rm!bjj rm 12.70b-bjj CHARACTER(11), EXTERNAL   :: CurDate
!rm!bjj rm 12.70b-bjj CHARACTER( 8), EXTERNAL   :: CurTime
!rm!bjj rm 12.70b-bjj CHARACTER(15), EXTERNAL   :: Flt2LStr
!rm!bjj rm 12.70b-bjj CHARACTER(11), EXTERNAL   :: Int2LStr
!rm
!rm
!rm
!rmIF (LogFlag .AND. ErrType /= 'INFO') THEN
!rm   INQUIRE ( FILE='Error.log' , EXIST=Exists )
!rm
!rm   IF ( Exists )  THEN
!rm      Frmt = '(////,A)'
!rm   ELSE
!rm      Frmt = '(''Error log started on             at       '',/A)'
!rm      WRITE(Frmt(24:35),'(A11)') CurDate()
!rm      WRITE(Frmt(39:44),'(A5 )') CurTime()
!rm   ENDIF
!rm
!rm   OPEN(UNIT = UnErrLog, FILE = 'error.log', STATUS = 'UNKNOWN',  &
!rm                   POSITION = 'APPEND', IOSTAT = IERR )
!rm
!rm   IF (IERR /= 0)  THEN
!rm      CALL WrScr ( 'ERROR OPENING FILE error.log')
!rm      Text = 'IOSTAT='//TRIM(Int2LStr(IERR))//' FILE INDEX= 99'
!rm      CALL WrScr ( Text)
!rm      CALL WrScr ('CHECK FILE')
!rm      CALL USRMES ( .TRUE. ,  &
!rm                   ' The error log could not be accessed.  '// &
!rm                   '  Errors will no longer be logged for this simulation.' &
!rm                     , IDADAMS , 'WARN' )
!rm      LogFlag = .FALSE.
!rm   ENDIF
!rmENDIF
!rm
!rmIF (LogFlag .AND. ErrType /= 'INFO') THEN
!rm   IF (FrstPass) THEN
!rm      NumWarn = 0
!rm      NumErr  = 0
!rm      WRITE(UnErrLog,Frmt)'============================================================ '
!rm      WRITE(UnErrLog,'(A/A)')'Errors logged by '//TRIM(Prog) , &
!rm                     ' on '//CurDate()//' at '//CurTime()//'.'
!rm      FrstPass = .FALSE.
!rm   ENDIF
!rm
!rm   IF (ErrType == 'WARN') THEN
!rm      NumWarn = NumWarn + 1
!rm      WRITE(UnErrLog,'(//,A,I2)')  'Simulation Warning #',NumWarn
!rm   ELSE
!rm      NumErr = NumErr + 1
!rm      WRITE(UnErrLog,'(//,A,I2)')  'Simulation Error #',NumErr
!rm   ENDIF
!rm   WRITE(UnErrLog,'(A)')' Simulation Time = '//TRIM(Flt2LStr(REAL( TIME, ReKi )))//' seconds.'
!rm   WRITE(UnErrLog,'(A)')  ' Error detected in Subroutine '//TRIM(Sub)//'.'
!rm   WRITE(UnErrLog, ErrFrmt) TRIM(MESAGE)
!rm   SELECT CASE (ErrType)
!rm    CASE ('WARN')
!rm      WRITE(UnErrLog,'(/1x,A)')'Warning Generated! Simulation continuing...'
!rm    CASE ('ERROR')
!rm      WRITE(UnErrLog,'(/1x,A)')'Non-Fatal Error. Simulation continuing...'
!rm    CASE ('FAULT')
!rm      WRITE(UnErrLog,'(/1x,A)')'FATAL ERROR. SIMULATION ABORTED.'
!rm   END SELECT !ErrType
!rm   CLOSE (UnErrLog)
!rmENDIF
!rm
!rmSELECT CASE (ErrType)
!rm CASE ('INFO')
!rm   SELECT CASE (TRIM(DynProg))
!rm      CASE ('ADAMS')
!rm         CALL USRMES (.TRUE., MESAGE, IDADAMS, ErrType)
!rm      CASE DEFAULT
!rm!bjj Start of proposed change v12.70b
!rm!RM   CALL WrScr1 (TRIM(MESAGE), TRIM(ErrFrmt))
!rm!bjj: this doesn't work for all formats      WRITE( Text, ErrFrmt ) TRIM(MESAGE)
!rm      CALL WrScr1 ( Mesage )
!rm!bjj End of proposed change
!rm   END SELECT
!rm CASE ('WARN')
!rm   Text = ' Program warning in Subroutine '//TRIM(Sub)//'.'
!rm   CALL WrScr1 (Text)
!rm!bjj Start of proposed change v12.70b
!rm!RM   CALL WrScr (TRIM(MESAGE), ErrFrmt)
!rm!this doesn't work for all formats:   WRITE( Text, ErrFrmt ) TRIM(MESAGE)
!rm   CALL WrScr ( Mesage )
!rm!bjj End of proposed change
!rm   CALL USRMES (.TRUE., 'Warning in '//TRIM(Sub), IDADAMS, ErrType)
!rm CASE ('ERROR')
!rm   Text = ' Program error in Subroutine '//TRIM(Sub)//'.'
!rm   CALL WrScr1 (Text)
!rm!bjj Start of proposed change v12.70b
!rm!RM   CALL WrScr (TRIM(MESAGE), ErrFrmt)
!rm!this doesn't work for all formats:   WRITE( Text, ErrFrmt ) TRIM(MESAGE)
!rm   CALL WrScr ( Mesage )
!rm!bjj End of proposed change
!rm   CALL USRMES (.TRUE., 'Non-fatal error in '//TRIM(Sub), IDADAMS, ErrType)
!rm CASE ('FAULT')
!rm!bjj Start of proposed change v12.70b
!rm!rm   Text = ' FATAL ERROR in Subroutine '//TRIM(Sub)//'.'
!rm!rm   CALL WrScr (Text,'(/,50(''*''),/A)')
!rm!rm   CALL WrScr (TRIM(MESAGE), ErrFrmt)
!rm!rm   CALL WrScr (' ','(A,50(''*''))')
!rm   CALL WrScr1('**************************************************')
!rm   Text = ' FATAL ERROR in Subroutine '//TRIM(Sub)//'.'
!rm   CALL WrScr (Text)
!rm   CALL WrScr ( MESAGE )
!rm   CALL WrScr( '**************************************************')
!rm!bjj End of proposed change
!rm   Text = TRIM(Int2LStr(NumWarn))// &
!rm         ' warnings recorded in error.log file for this simulation.'
!rm   CALL WrScr1 (Text)
!rm   Text = TRIM(Int2LStr(NumErr))// &
!rm         ' errors recorded in error.log file for this simulation.'
!rm   CALL WrScr (Text)
!rm   CALL USRMES (.TRUE., 'FATAL ERROR in '//TRIM(Sub), IDADAMS, ErrType)
!rmEND SELECT !ErrType
!rm
!rm
!rmRETURN
!rmEND SUBROUTINE ErrLog
!=======================================================================

END MODULE AeroGenSubs
!BJJ end of proposed change
