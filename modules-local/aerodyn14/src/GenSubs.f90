!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
 MODULE AeroGenSubs


   USE NWTC_LIBRARY
   USE AeroDyn_Types

IMPLICIT        NONE

! SUBROUTINE AllocArrays( Arg )
! SUBROUTINE ElemOpen( ElemFile )
! SUBROUTINE ElemOut( )

   INTEGER(IntKi) , PARAMETER :: MAXINFL = 6 


 CONTAINS
 ! ****************************************************
   SUBROUTINE AllocArrays ( InitInp, P, xc, xd, z, O, y, Arg  )
 !  Allocates space to the phenomenal number of arrays
 !  we use in this program
 ! ****************************************************

   TYPE(AD_InitInputType),       INTENT(INOUT)  :: InitInp
   TYPE(AD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD_ContinuousStateType), INTENT(INOUT)  :: xc          ! Initial continuous states
   TYPE(AD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(AD_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O!therState ! Initial other/optimization states
   TYPE(AD_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
   CHARACTER(*) :: Arg



   ! Passed Variables:



   ! Local Variables:

INTEGER(4)   :: Sttus

INTEGER :: NElm, NB, MaxTable, NumCl, NumFoil, NumElOut, NumWndElOut

NB          = P%NumBl
Nelm        = P%Element%Nelm
NumFoil     = P%AirFoil%NumFoil
NumCl       = P%AirFoil%NumCl
MaxTable    = P%AirFoil%MaxTable
NumElOut    = O%ElOut%NumElOut
NumWndElOut = O%ElOut%NumWndElOut


Sttus = 0.0

IF (Arg(1:7) == 'Element') THEN

   IF (.NOT. ALLOCATED(O%ElOut%ElPrList)) ALLOCATE ( O%ElOut%ElPrList(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for ElPrList array.' )
   O%ElOut%ElPrList ( : ) = 0

   IF (.NOT. ALLOCATED(O%ElOut%WndElPrList)) ALLOCATE ( O%ElOut%WndElPrList(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for WndElPrList array.' )
   O%ElOut%WndElPrList ( : ) = 0

   IF (.NOT. ALLOCATED(O%Element%A)) ALLOCATE ( O%Element%A(NELM,NB) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for A array.' )
   O%Element%A ( :, : ) = 0.0

   IF (.NOT. ALLOCATED(O%Element%AP)) ALLOCATE ( O%Element%AP(NELM,NB) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AP array.' )

! Beddoes arrays
   IF (P%Dstall) THEN

      IF (.NOT. ALLOCATED(O%Beddoes%ADOT)) ALLOCATE ( O%Beddoes%ADOT(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ADOT array.' )
      O%Beddoes%ADOT ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%ADOT1)) ALLOCATE ( O%Beddoes%ADOT1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ADOT1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%AFE)) ALLOCATE ( O%Beddoes%AFE(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AFE array.' )
      O%Beddoes%AFE(:,:) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%AFE1)) ALLOCATE ( O%Beddoes%AFE1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AFE1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%ANE)) ALLOCATE ( O%Beddoes%ANE(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ANE array.' )
      O%Beddoes%ANE ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%ANE1)) ALLOCATE ( O%Beddoes%ANE1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ANE1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%AOD)) ALLOCATE ( O%Beddoes%AOD(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AOD array.' )
      O%Beddoes%AOD = 0.0_ReKi
      
      IF (.NOT. ALLOCATED(O%Beddoes%AOL)) ALLOCATE ( O%Beddoes%AOL(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AOL array.' )
      O%Beddoes%AOL = 0.0_ReKi


      IF (.NOT. ALLOCATED(O%Beddoes%CDO)) ALLOCATE ( O%Beddoes%CDO(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CDO array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%CNA)) ALLOCATE ( O%Beddoes%CNA(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNA array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%CNP)) ALLOCATE ( O%Beddoes%CNP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNP array.' )
      O%Beddoes%CNP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%CNP1)) ALLOCATE ( O%Beddoes%CNP1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNP1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%CNPD)) ALLOCATE ( O%Beddoes%CNPD(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPD array.' )
      O%Beddoes%CNPD ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%CNPD1)) ALLOCATE ( O%Beddoes%CNPD1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPD1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%CNPOT)) ALLOCATE ( O%Beddoes%CNPOT(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPOT array.' )
      O%Beddoes%CNPOT ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%CNPOT1)) ALLOCATE ( O%Beddoes%CNPOT1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPOT1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%CNS)) ALLOCATE ( O%Beddoes%CNS(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNS array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%CNSL)) ALLOCATE ( O%Beddoes%CNSL(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNSL array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%CNV)) ALLOCATE ( O%Beddoes%CNV(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNV array.' )
      O%Beddoes%CNV ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%CVN)) ALLOCATE ( O%Beddoes%CVN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CVN array.' )
      O%Beddoes%CVN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%CVN1)) ALLOCATE ( O%Beddoes%CVN1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CVN1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%DF)) ALLOCATE ( O%Beddoes%DF(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DF array.' )
      O%Beddoes%DF( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%DFAFE)) ALLOCATE ( O%Beddoes%DFAFE(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFAFE array.' )
      O%Beddoes%DFAFE ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%DFAFE1)) ALLOCATE ( O%Beddoes%DFAFE1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFAFE1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%DFC)) ALLOCATE ( O%Beddoes%DFC(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFC array.' )
      O%Beddoes%DFC ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%DN)) ALLOCATE ( O%Beddoes%DN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DN array.' )
      O%Beddoes%DN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%DPP)) ALLOCATE ( O%Beddoes%DPP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DP array.' )
      O%Beddoes%DPP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%DQ)) ALLOCATE ( O%Beddoes%DQ(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DQ array.' )
      O%Beddoes%DQ ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%DQP)) ALLOCATE ( O%Beddoes%DQP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DQP array.' )
      O%Beddoes%DQP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%DQP1)) ALLOCATE ( O%Beddoes%DQP1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DQP1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%FSP)) ALLOCATE ( O%Beddoes%FSP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSP array.' )
      O%Beddoes%FSP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%FSP1)) ALLOCATE ( O%Beddoes%FSP1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSP1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%FSPC)) ALLOCATE ( O%Beddoes%FSPC(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSPC array.' )
      O%Beddoes%FSPC ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%FSPC1)) ALLOCATE ( O%Beddoes%FSPC1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSPC1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%OLDCNV)) ALLOCATE ( O%Beddoes%OLDCNV(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDCNV array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%OLDDF)) ALLOCATE ( O%Beddoes%OLDDF(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDF array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%OLDDFC)) ALLOCATE ( O%Beddoes%OLDDFC(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDFC array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%OLDDN)) ALLOCATE ( O%Beddoes%OLDDN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDN array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%OLDDPP)) ALLOCATE ( O%Beddoes%OLDDPP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDP array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%OLDDQ)) ALLOCATE ( O%Beddoes%OLDDQ(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDQ array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%OLDTAU)) ALLOCATE ( O%Beddoes%OLDTAU(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDTAU array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%OLDXN)) ALLOCATE ( O%Beddoes%OLDXN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDXN array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%OLDYN)) ALLOCATE ( O%Beddoes%OLDYN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDYN array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%QX)) ALLOCATE ( O%Beddoes%QX(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for QX array.' )
      O%Beddoes%QX ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%QX1)) ALLOCATE ( O%Beddoes%QX1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for QX1 array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%TAU)) ALLOCATE ( O%Beddoes%TAU(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for TAU array.' )
      O%Beddoes%TAU(:,:) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%XN)) ALLOCATE ( O%Beddoes%XN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for XN array.' )
      O%Beddoes%XN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%YN)) ALLOCATE ( O%Beddoes%YN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for YN array.' )
      O%Beddoes%YN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(O%Beddoes%OLDSEP)) ALLOCATE ( O%Beddoes%OLDSEP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDSEP array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%BEDSEP)) ALLOCATE ( O%Beddoes%BEDSEP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for BEDSEP array.' )
      O%Beddoes%BEDSEP(:,:) = .FALSE.

   ENDIF ! Beddoes arrays

   IF (.NOT. ALLOCATED(P%Blade%C)) ALLOCATE ( P%Blade%C(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for C array.' )

   IF (.NOT. ALLOCATED(P%Blade%DR)) ALLOCATE ( P%Blade%DR(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DR array.' )

   IF (.NOT. ALLOCATED(P%Element%RELM)) ALLOCATE ( P%Element%RELM(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for RELM array.' )

   IF (.NOT. ALLOCATED(P%Element%TWIST)) ALLOCATE ( P%Element%TWIST(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for TWIST array.' )

   IF (.NOT. ALLOCATED(P%Element%TLCNST)) ALLOCATE ( P%Element%TLCNST(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for TLCNST array.' )
   P%Element%TLCNST = 99.0

   IF (.NOT. ALLOCATED(P%Element%HLCNST)) ALLOCATE ( P%Element%HLCNST(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for HLCNST array.' )
   P%Element%HLCNST = 99.0

   IF (.NOT. ALLOCATED(P%AirFoil%NFOIL)) ALLOCATE ( P%AirFoil%NFOIL(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for NFOIL array.' )

   IF (.NOT. ALLOCATED(P%AirFoil%NLIFT)) ALLOCATE ( P%AirFoil%NLIFT(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for NLIFT array.' )

   IF (.NOT. ALLOCATED(P%AirFoil%NTables)) ALLOCATE ( P%AirFoil%NTables(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for NTables array.' )

   IF (P%Dyninfl .OR. O%Dyninit) THEN
      IF (.NOT. ALLOCATED(O%DynInflow%RMC_SAVE)) ALLOCATE ( O%DynInflow%RMC_SAVE ( NB, NELM, MAXINFL ) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for RMC_SAVE array.' )
      O%DynInflow%RMC_SAVE = 0.0

      IF (.NOT. ALLOCATED(O%DynInflow%RMS_SAVE)) ALLOCATE ( O%DynInflow%RMS_SAVE ( NB, NELM, MAXINFL ) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for RMS_SAVE array.' )
      O%DynInflow%RMS_SAVE = 0.0
   ENDIF

ELSEIF (Arg(1:7) == 'ElPrint') THEN

   IF ( O%ElOut%NumElOut > 0 ) THEN
      IF (.NOT. ALLOCATED(O%ElOut%AAA)) ALLOCATE ( O%ElOut%AAA(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AAA array.' )

      IF (.NOT. ALLOCATED(O%ElOut%AAP)) ALLOCATE ( O%ElOut%AAP(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AAP array.' )

      IF (.NOT. ALLOCATED(O%ElOut%ALF)) ALLOCATE ( O%ElOut%ALF(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ALF array.' )

      IF (.NOT. ALLOCATED(O%ElOut%CDD)) ALLOCATE ( O%ElOut%CDD(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CDD array.' )

      IF (.NOT. ALLOCATED(O%ElOut%CLL)) ALLOCATE ( O%ElOut%CLL(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CLL array.' )

      IF (.NOT. ALLOCATED(O%ElOut%CMM)) ALLOCATE ( O%ElOut%CMM(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CMM array.' )

      IF (.NOT. ALLOCATED(O%ElOut%CNN)) ALLOCATE ( O%ElOut%CNN(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNN array.' )

      IF (.NOT. ALLOCATED(O%ElOut%CTT)) ALLOCATE ( O%ElOut%CTT(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CTT array.' )

      IF (.NOT. ALLOCATED(O%ElOut%DFNSAV)) ALLOCATE ( O%ElOut%DFNSAV(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFNSAV array.' )

      IF (.NOT. ALLOCATED(O%ElOut%DFTSAV)) ALLOCATE ( O%ElOut%DFTSAV(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFTSAV array.' )

      IF (.NOT. ALLOCATED(O%ElOut%DynPres)) ALLOCATE ( O%ElOut%DynPres(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DynPres array.' )

      IF (.NOT. ALLOCATED(O%ElOut%PMM)) ALLOCATE ( O%ElOut%PMM(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for PMM array.' )

      IF (.NOT. ALLOCATED(O%ElOut%PITSAV)) ALLOCATE ( O%ElOut%PITSAV(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for PITSAV array.' )

      IF (.NOT. ALLOCATED(O%ElOut%ReyNum)) ALLOCATE ( O%ElOut%ReyNum(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ReyNum array.' )

      IF (.NOT. ALLOCATED(O%ElOut%ElPrNum)) ALLOCATE ( O%ElOut%ElPrNum(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for ElPrNum array.' )
      O%ElOut%ElPrNum ( : ) = 0

   END IF

   IF ( NumWndElOut > 0 ) THEN

      IF (.NOT. ALLOCATED(O%ElOut%WndElPrNum)) ALLOCATE ( O%ElOut%WndElPrNum(NumWndElOut) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for WndElPrNum array.' )
      O%ElOut%WndElPrNum ( : ) = 0

      IF (.NOT. ALLOCATED(O%ElOut%SaveVX)) ALLOCATE ( O%ElOut%SaveVX(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVX array.' )

      IF (.NOT. ALLOCATED(O%ElOut%SaveVY)) ALLOCATE ( O%ElOut%SaveVY(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVY array.' )

      IF (.NOT. ALLOCATED(O%ElOut%SaveVZ)) ALLOCATE ( O%ElOut%SaveVZ(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVZ array.' )

   END IF

ELSEIF (Arg(1:8) == 'Aerodata') THEN

   IF (.NOT. ALLOCATED(O%AirFoil%AL)) ALLOCATE ( O%AirFoil%AL(NumFoil,NumCL) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AL array.' )

   IF (.NOT. ALLOCATED(O%AirFoil%CD)) ALLOCATE ( O%AirFoil%CD(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CD array.' )

   IF (.NOT. ALLOCATED(O%AirFoil%CL)) ALLOCATE ( O%AirFoil%CL(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CL array.' )

   IF (.NOT. ALLOCATED(O%AirFoil%CM)) ALLOCATE ( O%AirFoil%CM(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CM array.' )

   IF (.NOT. ALLOCATED(p%AirFoil%MulTabMet)) ALLOCATE ( p%AirFoil%MulTabMet(NumFoil,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for MulTabMet array.' )

   IF (P%DSTALL) THEN

      IF (.NOT. ALLOCATED(O%Beddoes%FTB)) ALLOCATE ( O%Beddoes%FTB(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FTB array.' )

      IF (.NOT. ALLOCATED(O%Beddoes%FTBC)) ALLOCATE ( O%Beddoes%FTBC(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FTBC array.' )

   ENDIF ! Beddoes arrays

!jm we did not recognize the argument consider that an error
ELSE

   CALL ProgAbort( 'Unknown switch argument to AllocArrays' )

ENDIF



RETURN
END SUBROUTINE AllocArrays

 ! *****************************************************
   SUBROUTINE ElemOpen (ElemFile, P, O, ErrStat, ErrMsg, AD_Ver )
 !  This subroutine opens the element output file and writes
 !   column headings separated by tab characters
 !  ElemFile = file name
 ! *****************************************************

   ! Passed Variables:
   CHARACTER(*), INTENT(IN) :: ElemFile
   TYPE(AD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: O !therState  Initial other/optimization states
   INTEGER(IntKi),               INTENT(OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   TYPE(ProgDesc)    ,           INTENT(IN)     :: AD_Ver

   ! Local Variables:

INTEGER(4)     :: JE

INTEGER(4)     :: JB       ! Counter for number of blades

CHARACTER(  2) :: Dst_Unit
CHARACTER(  2) :: Frc_Unit
CHARACTER(140) :: Frmt
CHARACTER(  3) :: Prs_Unit


ErrStat = ErrID_None
ErrMsg = ""



IF (O%ElOut%NumWndElOut > 0) THEN
   CALL GetNewUnit(P%UnWndOut)
   CALL OpenFOutFile (P%UnWndOut, TRIM(ElemFile)//'.wind', ErrStat, ErrMsg)
   IF (ErrStat /= ErrID_None) RETURN
   WRITE (P%UnWndOut,"( 'This file was generated by ' , A , ' on ' , A , ' at ' , A , '.' )")  &
        GETNVD(AD_Ver), CurDate(), CurTime()
ENDIF


 ! Open the Element Print file if requested
IF (p%ELEMPRN) THEN
   CALL GetNewUnit(P%UnElem)
   CALL OpenFOutFile (p%UnElem, TRIM(ElemFile), ErrStat, ErrMsg)
   IF (ErrStat /= ErrID_None) RETURN
   WRITE (p%UnElem,"( 'This file was generated by ' , A , ' on ' , A , ' at ' , A , '.' )")  &
        GETNVD(AD_Ver), CurDate(), CurTime()
ELSE
   RETURN
ENDIF

 ! Set the units labels
IF (p%SIUNIT) THEN
  Dst_Unit = 'm'
  Frc_Unit = 'N'
  Prs_Unit = 'Pa'
ELSE
  Dst_Unit = 'ft'
  Frc_Unit = 'lb'
  Prs_Unit = 'psf'
ENDIF


Frmt = '( A4, 3(A1,A2,I2.2),    (: A1, A, I2.2 ) )'

IF ( p%PMOMENT ) THEN
   WRITE(Frmt(22:24), '(I3)') 14*O%ElOut%NumElOut
   WRITE(p%UnElem, Frmt) 'Time',                    &
               TAB,    'VX',       p%Element%NELM,          &
               TAB,    'VY',       p%Element%NELM,          &
               TAB,    'VZ',       p%Element%NELM,          &
             ( TAB,    'Alpha',    O%ElOut%ElPrNum(JE),  &
               TAB,    'DynPres',  O%ElOut%ElPrNum(JE),  &
               TAB,    'CLift',    O%ElOut%ElPrNum(JE),  &
               TAB,    'CDrag',    O%ElOut%ElPrNum(JE),  &
               TAB,    'CNorm',    O%ElOut%ElPrNum(JE),  &
               TAB,    'CTang',    O%ElOut%ElPrNum(JE),  &
               TAB,    'CMomt',    O%ElOut%ElPrNum(JE),  &
               TAB,    'Pitch',    O%ElOut%ElPrNum(JE),  &
               TAB,    'AxInd',    O%ElOut%ElPrNum(JE),  &
               TAB,    'TanInd',   O%ElOut%ElPrNum(JE),  &
               TAB,    'ForcN',    O%ElOut%ElPrNum(JE),  &
               TAB,    'ForcT',    O%ElOut%ElPrNum(JE),  &
               TAB,    'Pmomt',    O%ElOut%ElPrNum(JE),  &
               TAB,    'ReNum',    O%ElOut%ElPrNum(JE),  &
                         JE = 1, O%ElOut%NumElOut )

   Frmt = '( A5, 3(A1,A8),    (: A1, A ) )'
   WRITE(Frmt(17:19), '(I3)') 14*O%ElOut%NumElOut
   WRITE(p%UnElem, Frmt) '(sec)',                       &
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
                         JE = 1, O%ElOut%NumElOut )

ELSE
   WRITE(Frmt(22:24), '(I3)') 12*O%ElOut%NumElOut
   WRITE(p%UnElem, Frmt) 'Time',                    &
               TAB,    'VX',       p%Element%NELM,          &
               TAB,    'VY',       p%Element%NELM,          &
               TAB,    'VZ',       p%Element%NELM,          &
             ( TAB,    'Alpha',    O%ElOut%ElPrNum(JE),  &
               TAB,    'DynPres',  O%ElOut%ElPrNum(JE),  &
               TAB,    'CLift',    O%ElOut%ElPrNum(JE),  &
               TAB,    'CDrag',    O%ElOut%ElPrNum(JE),  &
               TAB,    'CNorm',    O%ElOut%ElPrNum(JE),  &
               TAB,    'CTang',    O%ElOut%ElPrNum(JE),  &
               TAB,    'Pitch',    O%ElOut%ElPrNum(JE),  &
               TAB,    'AxInd',    O%ElOut%ElPrNum(JE),  &
               TAB,    'TanInd',   O%ElOut%ElPrNum(JE),  &
               TAB,    'ForcN',    O%ElOut%ElPrNum(JE),  &
               TAB,    'ForcT',    O%ElOut%ElPrNum(JE),  &
               TAB,    'ReNum',    O%ElOut%ElPrNum(JE),  &
                         JE = 1, O%ElOut%NumElOut )

   Frmt = '( A5, 3(A1,A8),    (: A1, A ) )'
   WRITE(Frmt(17:19), '(I3)') 12*O%ElOut%NumElOut
   WRITE(p%UnElem, Frmt) '(sec)',                       &
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
                         JE = 1, O%ElOut%NumElOut )
ENDIF

IF ( O%ElOut%NumWndElOut > 0 ) THEN
   Frmt = '( A4, XXX(A1,A2,I2.2,"-B",I1.1) )'
   WRITE(Frmt(7:9), '(I3)') 3*O%ElOut%NumWndElOut*p%NumBl

   WRITE(p%UnWndOut, Frmt) 'Time',          &
             ( ( TAB, 'VX',  O%ElOut%WndElPrNum(JE),  JB,   &
                 TAB, 'VY',  O%ElOut%WndElPrNum(JE),  JB,   &
                 TAB, 'VZ',  O%ElOut%WndElPrNum(JE),  JB,   &
                   JE = 1, O%ElOut%NumWndElOut ) , &
                   JB = 1, p%NumBl )

   Frmt = '( A5, XXX(A1,A8) )'
   WRITE(Frmt(7:9), '(I3)') 3*O%ElOut%NumWndElOut*p%NumBl

   WRITE(p%UnWndOut, Frmt)   '(sec)',                       &
               ( TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
                 TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
                 TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
                                  JE = 1, O%ElOut%NumWndElOut*p%NumBl )
ENDIF

RETURN
END SUBROUTINE ElemOpen


 ! *****************************************************
   SUBROUTINE ElemOut( time, P, O )
 !  This subroutine writes the element output values
 !   for the desired elements
 ! *****************************************************

 
   REAL(DbKi), INTENT(IN) :: time
   TYPE(AD_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD_OtherStateType),      INTENT(IN)  :: O !therState  Initial other/optimization states


   ! Local Variables:

   INTEGER(IntKi)     :: JE
   INTEGER(IntKi)     :: JB    ! Counter for number of blades
   CHARACTER(30)  :: Frmt


 ! Write the element data if requested
IF (p%ELEMPRN) THEN

   Frmt = '( F10.3,    ( : A1, ES12.5 ) )'

   IF ( P%PMOMENT ) THEN
      WRITE(Frmt(10:12), '(I3)') 14*O%ElOut%NumElOut + 3
      WRITE(p%UnElem,Frmt) TIME,                 &
                     TAB,   O%ElOut%VXSAV,       &
                     TAB,   O%ElOut%VYSAV,       &
                     TAB,   O%ElOut%VZSAV,       &
                   ( TAB,   O%ElOut%ALF    (JE), &
                     TAB,   O%ElOut%DynPres(JE), &
                     TAB,   O%ElOut%CLL    (JE), &
                     TAB,   O%ElOut%CDD    (JE), &
                     TAB,   O%ElOut%CNN    (JE), &
                     TAB,   O%ElOut%CTT    (JE), &
                     TAB,   O%ElOut%CMM    (JE), &
                     TAB,   O%ElOut%PITSAV (JE), &
                     TAB,   O%ElOut%AAA    (JE), &
                     TAB,   O%ElOut%AAP    (JE), &
                     TAB,   O%ElOut%DFNSAV (JE), &
                     TAB,   O%ElOut%DFTSAV (JE), &
                     TAB,   O%ElOut%PMM    (JE), &
                     TAB,   O%ElOut%ReyNum (JE), &
                            JE= 1, O%ElOut%NumElOut )


   ELSE
      WRITE(Frmt(10:12), '(I3)') 12*O%ElOut%NumElOut + 3
      WRITE(p%UnElem,Frmt) TIME,                 &
                     TAB,   O%ElOut%VXSAV,       &
                     TAB,   O%ElOut%VYSAV,       &
                     TAB,   O%ElOut%VZSAV,       &
                   ( TAB,   O%ElOut%ALF    (JE), &
                     TAB,   O%ElOut%DynPres(JE), &
                     TAB,   O%ElOut%CLL    (JE), &
                     TAB,   O%ElOut%CDD    (JE), &
                     TAB,   O%ElOut%CNN    (JE), &
                     TAB,   O%ElOut%CTT    (JE), &
                     TAB,   O%ElOut%PITSAV (JE), &
                     TAB,   O%ElOut%AAA    (JE), &
                     TAB,   O%ElOut%AAP    (JE), &
                     TAB,   O%ElOut%DFNSAV (JE), &
                     TAB,   O%ElOut%DFTSAV (JE), &
                     TAB,   O%ElOut%ReyNum (JE), &
                            JE= 1, O%ElOut%NumElOut )
   ENDIF ! PMOMENT

IF (O%ElOut%NumWndElOut > 0) THEN

   WRITE(Frmt(10:12), '(I3)') 3*O%ElOut%NumWndElOut*p%NumBl
   WRITE(p%UnWndOut,Frmt) TIME,                   &
                    ( (  TAB, O%ElOut%SaveVX( JE, JB ), &
                         TAB, O%ElOut%SaveVY( JE, JB ), &
                         TAB, O%ElOut%SaveVZ( JE, JB ), &
                              JE = 1,O%ElOut%NumWndElOut ), JB = 1,p%NumBl )
ENDIF

ENDIF ! ELEMPRN



RETURN
END SUBROUTINE ElemOut
!=======================================================================

END MODULE AeroGenSubs
