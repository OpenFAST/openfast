!**********************************************************************************************************************************
 MODULE AeroGenSubs


   USE NWTC_LIBRARY
   USE AeroDyn14_Types

IMPLICIT        NONE

! SUBROUTINE AllocArrays( Arg )
! SUBROUTINE ElemOpen( ElemFile )
! SUBROUTINE ElemOut( )

   INTEGER(IntKi) , PARAMETER :: MAXINFL = 6 


 CONTAINS
 ! ****************************************************
   SUBROUTINE AllocArrays ( InitInp, P, xc, xd, z, m, y, Arg  )
 !  Allocates space to the phenomenal number of arrays
 !  we use in this program
 ! ****************************************************
   ! Passed Variables:
   TYPE(AD14_InitInputType),       INTENT(INOUT)  :: InitInp
   TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD14_ContinuousStateType), INTENT(INOUT)  :: xc          ! Initial continuous states
   TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(AD14_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Initial misc/optimization variables
   TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
   CHARACTER(*),                   INTENT(IN   )  :: Arg

   ! Local Variables:

   INTEGER(4)   :: Sttus

   INTEGER :: NElm, NB, MaxTable, NumCl, NumFoil, NumElOut, NumWndElOut


!bjj: I really don't understand why these aren't 3 separate subroutines...


NB          = P%NumBl
Nelm        = P%Element%Nelm
NumFoil     = P%AirFoil%NumFoil
NumCl       = P%AirFoil%NumCl
MaxTable    = P%AirFoil%MaxTable
NumElOut    = m%ElOut%NumElOut
NumWndElOut = m%ElOut%NumWndElOut


Sttus = 0.0

IF (Arg(1:7) == 'Element') THEN

   IF (.NOT. ALLOCATED(m%ElOut%ElPrList)) ALLOCATE ( m%ElOut%ElPrList(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for ElPrList array.' )
   m%ElOut%ElPrList ( : ) = 0

   IF (.NOT. ALLOCATED(m%ElOut%WndElPrList)) ALLOCATE ( m%ElOut%WndElPrList(NELM) , STAT=Sttus )
   IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for WndElPrList array.' )
   m%ElOut%WndElPrList ( : ) = 0

   IF (.NOT. ALLOCATED(m%Element%A)) ALLOCATE ( m%Element%A(NELM,NB) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for A array.' )
   m%Element%A ( :, : ) = 0.0

   IF (.NOT. ALLOCATED(m%Element%AP)) ALLOCATE ( m%Element%AP(NELM,NB) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AP array.' )

! Beddoes arrays
   IF (P%Dstall) THEN

      IF (.NOT. ALLOCATED(m%Beddoes%ADOT)) ALLOCATE ( m%Beddoes%ADOT(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ADOT array.' )
      m%Beddoes%ADOT ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%ADOT1)) ALLOCATE ( m%Beddoes%ADOT1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ADOT1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%AFE)) ALLOCATE ( m%Beddoes%AFE(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AFE array.' )
      m%Beddoes%AFE(:,:) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%AFE1)) ALLOCATE ( m%Beddoes%AFE1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AFE1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%ANE)) ALLOCATE ( m%Beddoes%ANE(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ANE array.' )
      m%Beddoes%ANE ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%ANE1)) ALLOCATE ( m%Beddoes%ANE1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ANE1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%AOD)) ALLOCATE ( m%Beddoes%AOD(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AOD array.' )
      m%Beddoes%AOD = 0.0_ReKi
      
      IF (.NOT. ALLOCATED(m%Beddoes%AOL)) ALLOCATE ( m%Beddoes%AOL(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AOL array.' )
      m%Beddoes%AOL = 0.0_ReKi


      IF (.NOT. ALLOCATED(m%Beddoes%CDO)) ALLOCATE ( m%Beddoes%CDO(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CDO array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%CNA)) ALLOCATE ( m%Beddoes%CNA(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNA array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%CNP)) ALLOCATE ( m%Beddoes%CNP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNP array.' )
      m%Beddoes%CNP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%CNP1)) ALLOCATE ( m%Beddoes%CNP1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNP1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%CNPD)) ALLOCATE ( m%Beddoes%CNPD(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPD array.' )
      m%Beddoes%CNPD ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%CNPD1)) ALLOCATE ( m%Beddoes%CNPD1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPD1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%CNPOT)) ALLOCATE ( m%Beddoes%CNPOT(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPOT array.' )
      m%Beddoes%CNPOT ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%CNPOT1)) ALLOCATE ( m%Beddoes%CNPOT1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNPOT1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%CNS)) ALLOCATE ( m%Beddoes%CNS(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNS array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%CNSL)) ALLOCATE ( m%Beddoes%CNSL(NELM,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNSL array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%CNV)) ALLOCATE ( m%Beddoes%CNV(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNV array.' )
      m%Beddoes%CNV ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%CVN)) ALLOCATE ( m%Beddoes%CVN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CVN array.' )
      m%Beddoes%CVN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%CVN1)) ALLOCATE ( m%Beddoes%CVN1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CVN1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%DF)) ALLOCATE ( m%Beddoes%DF(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DF array.' )
      m%Beddoes%DF( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%DFAFE)) ALLOCATE ( m%Beddoes%DFAFE(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFAFE array.' )
      m%Beddoes%DFAFE ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%DFAFE1)) ALLOCATE ( m%Beddoes%DFAFE1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFAFE1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%DFC)) ALLOCATE ( m%Beddoes%DFC(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFC array.' )
      m%Beddoes%DFC ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%DN)) ALLOCATE ( m%Beddoes%DN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DN array.' )
      m%Beddoes%DN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%DPP)) ALLOCATE ( m%Beddoes%DPP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DP array.' )
      m%Beddoes%DPP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%DQ)) ALLOCATE ( m%Beddoes%DQ(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DQ array.' )
      m%Beddoes%DQ ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%DQP)) ALLOCATE ( m%Beddoes%DQP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DQP array.' )
      m%Beddoes%DQP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%DQP1)) ALLOCATE ( m%Beddoes%DQP1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DQP1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%FSP)) ALLOCATE ( m%Beddoes%FSP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSP array.' )
      m%Beddoes%FSP ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%FSP1)) ALLOCATE ( m%Beddoes%FSP1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSP1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%FSPC)) ALLOCATE ( m%Beddoes%FSPC(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSPC array.' )
      m%Beddoes%FSPC ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%FSPC1)) ALLOCATE ( m%Beddoes%FSPC1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FSPC1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%OLDCNV)) ALLOCATE ( m%Beddoes%OLDCNV(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDCNV array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%OLDDF)) ALLOCATE ( m%Beddoes%OLDDF(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDF array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%OLDDFC)) ALLOCATE ( m%Beddoes%OLDDFC(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDFC array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%OLDDN)) ALLOCATE ( m%Beddoes%OLDDN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDN array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%OLDDPP)) ALLOCATE ( m%Beddoes%OLDDPP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDP array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%OLDDQ)) ALLOCATE ( m%Beddoes%OLDDQ(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDDQ array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%OLDTAU)) ALLOCATE ( m%Beddoes%OLDTAU(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDTAU array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%OLDXN)) ALLOCATE ( m%Beddoes%OLDXN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDXN array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%OLDYN)) ALLOCATE ( m%Beddoes%OLDYN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDYN array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%QX)) ALLOCATE ( m%Beddoes%QX(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for QX array.' )
      m%Beddoes%QX ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%QX1)) ALLOCATE ( m%Beddoes%QX1(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for QX1 array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%TAU)) ALLOCATE ( m%Beddoes%TAU(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for TAU array.' )
      m%Beddoes%TAU(:,:) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%XN)) ALLOCATE ( m%Beddoes%XN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for XN array.' )
      m%Beddoes%XN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%YN)) ALLOCATE ( m%Beddoes%YN(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for YN array.' )
      m%Beddoes%YN ( :, : ) = 0.0

      IF (.NOT. ALLOCATED(m%Beddoes%OLDSEP)) ALLOCATE ( m%Beddoes%OLDSEP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLDSEP array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%BEDSEP)) ALLOCATE ( m%Beddoes%BEDSEP(NELM,NB) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for BEDSEP array.' )
      m%Beddoes%BEDSEP(:,:) = .FALSE.

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

   IF (P%Dyninfl .OR. m%Dyninit) THEN
      IF (.NOT. ALLOCATED(m%DynInflow%RMC_SAVE)) ALLOCATE ( m%DynInflow%RMC_SAVE ( NB, NELM, MAXINFL ) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for RMC_SAVE array.' )
      m%DynInflow%RMC_SAVE = 0.0

      IF (.NOT. ALLOCATED(m%DynInflow%RMS_SAVE)) ALLOCATE ( m%DynInflow%RMS_SAVE ( NB, NELM, MAXINFL ) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for RMS_SAVE array.' )
      m%DynInflow%RMS_SAVE = 0.0
   ENDIF

ELSEIF (Arg(1:7) == 'ElPrint') THEN

   IF ( m%ElOut%NumElOut > 0 ) THEN
      IF (.NOT. ALLOCATED(m%ElOut%AAA)) ALLOCATE ( m%ElOut%AAA(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AAA array.' )

      IF (.NOT. ALLOCATED(m%ElOut%AAP)) ALLOCATE ( m%ElOut%AAP(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AAP array.' )

      IF (.NOT. ALLOCATED(m%ElOut%ALF)) ALLOCATE ( m%ElOut%ALF(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ALF array.' )

      IF (.NOT. ALLOCATED(m%ElOut%CDD)) ALLOCATE ( m%ElOut%CDD(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CDD array.' )

      IF (.NOT. ALLOCATED(m%ElOut%CLL)) ALLOCATE ( m%ElOut%CLL(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CLL array.' )

      IF (.NOT. ALLOCATED(m%ElOut%CMM)) ALLOCATE ( m%ElOut%CMM(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CMM array.' )

      IF (.NOT. ALLOCATED(m%ElOut%CNN)) ALLOCATE ( m%ElOut%CNN(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CNN array.' )

      IF (.NOT. ALLOCATED(m%ElOut%CTT)) ALLOCATE ( m%ElOut%CTT(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CTT array.' )

      IF (.NOT. ALLOCATED(m%ElOut%DFNSAV)) ALLOCATE ( m%ElOut%DFNSAV(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFNSAV array.' )

      IF (.NOT. ALLOCATED(m%ElOut%DFTSAV)) ALLOCATE ( m%ElOut%DFTSAV(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DFTSAV array.' )

      IF (.NOT. ALLOCATED(m%ElOut%DynPres)) ALLOCATE ( m%ElOut%DynPres(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for DynPres array.' )

      IF (.NOT. ALLOCATED(m%ElOut%PMM)) ALLOCATE ( m%ElOut%PMM(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for PMM array.' )

      IF (.NOT. ALLOCATED(m%ElOut%PITSAV)) ALLOCATE ( m%ElOut%PITSAV(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for PITSAV array.' )

      IF (.NOT. ALLOCATED(m%ElOut%ReyNum)) ALLOCATE ( m%ElOut%ReyNum(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for ReyNum array.' )

      IF (.NOT. ALLOCATED(m%ElOut%Gamma)) ALLOCATE ( m%ElOut%Gamma(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for Gamma array.' )

      IF (.NOT. ALLOCATED(m%ElOut%ElPrNum)) ALLOCATE ( m%ElOut%ElPrNum(NumElOut) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for ElPrNum array.' )
      m%ElOut%ElPrNum ( : ) = 0

   END IF

   IF ( NumWndElOut > 0 ) THEN

      IF (.NOT. ALLOCATED(m%ElOut%WndElPrNum)) ALLOCATE ( m%ElOut%WndElPrNum(NumWndElOut) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for WndElPrNum array.' )
      m%ElOut%WndElPrNum ( : ) = 0

      IF (.NOT. ALLOCATED(m%ElOut%SaveVX)) ALLOCATE ( m%ElOut%SaveVX(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVX array.' )

      IF (.NOT. ALLOCATED(m%ElOut%SaveVY)) ALLOCATE ( m%ElOut%SaveVY(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVY array.' )

      IF (.NOT. ALLOCATED(m%ElOut%SaveVZ)) ALLOCATE ( m%ElOut%SaveVZ(NumWndElOut,NB) , STAT=Sttus )
      IF ( Sttus /= 0 )  CALL ProgAbort ( ' Error allocating memory for SaveVZ array.' )

   END IF

ELSEIF (Arg(1:8) == 'Aerodata') THEN

   IF (.NOT. ALLOCATED(m%AirFoil%AL)) ALLOCATE ( m%AirFoil%AL(NumFoil,NumCL) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for AL array.' )

   IF (.NOT. ALLOCATED(m%AirFoil%CD)) ALLOCATE ( m%AirFoil%CD(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CD array.' )

   IF (.NOT. ALLOCATED(m%AirFoil%CL)) ALLOCATE ( m%AirFoil%CL(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CL array.' )

   IF (.NOT. ALLOCATED(m%AirFoil%CM)) ALLOCATE ( m%AirFoil%CM(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CM array.' )

   IF (.NOT. ALLOCATED(p%AirFoil%MulTabMet)) ALLOCATE ( p%AirFoil%MulTabMet(NumFoil,MAXTABLE) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for MulTabMet array.' )

   IF (P%DSTALL) THEN

      IF (.NOT. ALLOCATED(m%Beddoes%FTB)) ALLOCATE ( m%Beddoes%FTB(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FTB array.' )

      IF (.NOT. ALLOCATED(m%Beddoes%FTBC)) ALLOCATE ( m%Beddoes%FTBC(NumFoil,NumCL,MAXTABLE) , STAT=Sttus )
      IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for FTBC array.' )

   ENDIF ! Beddoes arrays

!jm we did not recognize the argument consider that an error
ELSE

   CALL ProgAbort( 'Unknown switch argument to AllocArrays' )

ENDIF



RETURN
END SUBROUTINE AllocArrays

 ! *****************************************************
   SUBROUTINE ElemOpen (ElemFile, P, m, ErrStat, ErrMsg, AD14_Ver )
 !  This subroutine opens the element output file and writes
 !   column headings separated by tab characters
 !  ElemFile = file name
 ! *****************************************************

   ! Passed Variables:
   CHARACTER(*), INTENT(IN) :: ElemFile
   TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! misc/optimization variables
   INTEGER(IntKi),                 INTENT(OUT)    :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(OUT)    :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   TYPE(ProgDesc)    ,             INTENT(IN)     :: AD14_Ver

   ! Local Variables:

INTEGER(4)     :: JE

INTEGER(4)     :: JB       ! Counter for number of blades

CHARACTER(  2) :: Dst_Unit
CHARACTER(  2) :: Frc_Unit
CHARACTER(140) :: Frmt
CHARACTER(  3) :: Prs_Unit


ErrStat = ErrID_None
ErrMsg = ""



IF (m%ElOut%NumWndElOut > 0) THEN
   CALL GetNewUnit(P%UnWndOut)
   CALL OpenFOutFile (P%UnWndOut, TRIM(ElemFile)//'.wind', ErrStat, ErrMsg)
   IF (ErrStat /= ErrID_None) RETURN
   WRITE (P%UnWndOut,"( 'This file was generated by ' , A , ' on ' , A , ' at ' , A , '.' )")  &
        TRIM(GETNVD(AD14_Ver)), CurDate(), CurTime()
ENDIF


 ! Open the Element Print file if requested
IF (p%ELEMPRN) THEN
   CALL GetNewUnit(P%UnElem)
   CALL OpenFOutFile (p%UnElem, TRIM(ElemFile), ErrStat, ErrMsg)
   IF (ErrStat >= AbortErrLev) RETURN
   WRITE (p%UnElem,"(/, 'This file was generated by ' , A , ' on ' , A , ' at ' , A , '.' )")  &
        TRIM(GETNVD(AD14_Ver)), CurDate(), CurTime()
   WRITE (p%UnElem,'(/,/,/)') ! write some blank lines so the output file so the headers are on the same line as FAST's
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
   WRITE(Frmt(22:24), '(I3)') 15*m%ElOut%NumElOut
   WRITE(p%UnElem, Frmt) 'Time',                    &
               TAB,    'VX',       p%Element%NELM,          &
               TAB,    'VY',       p%Element%NELM,          &
               TAB,    'VZ',       p%Element%NELM,          &
             ( TAB,    'Alpha',    m%ElOut%ElPrNum(JE),  &
               TAB,    'DynPres',  m%ElOut%ElPrNum(JE),  &
               TAB,    'CLift',    m%ElOut%ElPrNum(JE),  &
               TAB,    'CDrag',    m%ElOut%ElPrNum(JE),  &
               TAB,    'CNorm',    m%ElOut%ElPrNum(JE),  &
               TAB,    'CTang',    m%ElOut%ElPrNum(JE),  &
               TAB,    'CMomt',    m%ElOut%ElPrNum(JE),  &
               TAB,    'Pitch',    m%ElOut%ElPrNum(JE),  &
               TAB,    'AxInd',    m%ElOut%ElPrNum(JE),  &
               TAB,    'TanInd',   m%ElOut%ElPrNum(JE),  &
               TAB,    'ForcN',    m%ElOut%ElPrNum(JE),  &
               TAB,    'ForcT',    m%ElOut%ElPrNum(JE),  &
               TAB,    'Pmomt',    m%ElOut%ElPrNum(JE),  &
               TAB,    'ReNum',    m%ElOut%ElPrNum(JE),  &
               TAB,    'Gamma',    m%ElOut%ElPrNum(JE),  &
                         JE = 1, m%ElOut%NumElOut )

   Frmt = '( A5, 3(A1,A8),    (: A1, A ) )'
   WRITE(Frmt(17:19), '(I3)') 15*m%ElOut%NumElOut
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
               TAB,    '(m^2/sec)',                   &
                         JE = 1, m%ElOut%NumElOut )

ELSE
   WRITE(Frmt(22:24), '(I3)') 13*m%ElOut%NumElOut
   WRITE(p%UnElem, Frmt) 'Time',                    &
               TAB,    'VX',       p%Element%NELM,          &
               TAB,    'VY',       p%Element%NELM,          &
               TAB,    'VZ',       p%Element%NELM,          &
             ( TAB,    'Alpha',    m%ElOut%ElPrNum(JE),  &
               TAB,    'DynPres',  m%ElOut%ElPrNum(JE),  &
               TAB,    'CLift',    m%ElOut%ElPrNum(JE),  &
               TAB,    'CDrag',    m%ElOut%ElPrNum(JE),  &
               TAB,    'CNorm',    m%ElOut%ElPrNum(JE),  &
               TAB,    'CTang',    m%ElOut%ElPrNum(JE),  &
               TAB,    'Pitch',    m%ElOut%ElPrNum(JE),  &
               TAB,    'AxInd',    m%ElOut%ElPrNum(JE),  &
               TAB,    'TanInd',   m%ElOut%ElPrNum(JE),  &
               TAB,    'ForcN',    m%ElOut%ElPrNum(JE),  &
               TAB,    'ForcT',    m%ElOut%ElPrNum(JE),  &
               TAB,    'ReNum',    m%ElOut%ElPrNum(JE),  &
               TAB,    'Gamma',    m%ElOut%ElPrNum(JE),  &
                         JE = 1, m%ElOut%NumElOut )

   Frmt = '( A5, 3(A1,A8),    (: A1, A ) )'
   WRITE(Frmt(17:19), '(I3)') 13*m%ElOut%NumElOut
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
               TAB,    '(m^2/sec)',                   &
                         JE = 1, m%ElOut%NumElOut )
ENDIF

IF ( m%ElOut%NumWndElOut > 0 ) THEN
   Frmt = '( A4, XXX(A1,A2,I2.2,"-B",I1.1) )'
   WRITE(Frmt(7:9), '(I3)') 3*m%ElOut%NumWndElOut*p%NumBl

   WRITE(p%UnWndOut, Frmt) 'Time',          &
             ( ( TAB, 'VX',  m%ElOut%WndElPrNum(JE),  JB,   &
                 TAB, 'VY',  m%ElOut%WndElPrNum(JE),  JB,   &
                 TAB, 'VZ',  m%ElOut%WndElPrNum(JE),  JB,   &
                   JE = 1, m%ElOut%NumWndElOut ) , &
                   JB = 1, p%NumBl )

   Frmt = '( A5, XXX(A1,A8) )'
   WRITE(Frmt(7:9), '(I3)') 3*m%ElOut%NumWndElOut*p%NumBl

   WRITE(p%UnWndOut, Frmt)   '(sec)',                       &
               ( TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
                 TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
                 TAB,    '('//TRIM(Dst_Unit)//'/sec)',  &
                                  JE = 1, m%ElOut%NumWndElOut*p%NumBl )
ENDIF

RETURN
END SUBROUTINE ElemOpen


 ! *****************************************************
   SUBROUTINE ElemOut( time, P, m )
 !  This subroutine writes the element output values
 !   for the desired elements
 ! *****************************************************

 
   REAL(DbKi), INTENT(IN) :: time
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(IN)  :: m           ! misc/optimization variables


   ! Local Variables:

   INTEGER(IntKi)     :: JE
   INTEGER(IntKi)     :: JB    ! Counter for number of blades
   CHARACTER(30)  :: Frmt


 ! Write the element data if requested
IF (p%ELEMPRN) THEN

   Frmt = '( F10.3,    ( : A1, ES12.5 ) )'

   IF ( P%PMOMENT ) THEN
      WRITE(Frmt(10:12), '(I3)') 15*m%ElOut%NumElOut + 3
      WRITE(p%UnElem,Frmt) TIME,                 &
                     TAB,   m%ElOut%VXSAV,       &
                     TAB,   m%ElOut%VYSAV,       &
                     TAB,   m%ElOut%VZSAV,       &
                   ( TAB,   m%ElOut%ALF    (JE), &
                     TAB,   m%ElOut%DynPres(JE), &
                     TAB,   m%ElOut%CLL    (JE), &
                     TAB,   m%ElOut%CDD    (JE), &
                     TAB,   m%ElOut%CNN    (JE), &
                     TAB,   m%ElOut%CTT    (JE), &
                     TAB,   m%ElOut%CMM    (JE), &
                     TAB,   m%ElOut%PITSAV (JE), &
                     TAB,   m%ElOut%AAA    (JE), &
                     TAB,   m%ElOut%AAP    (JE), &
                     TAB,   m%ElOut%DFNSAV (JE), &
                     TAB,   m%ElOut%DFTSAV (JE), &
                     TAB,   m%ElOut%PMM    (JE), &
                     TAB,   m%ElOut%ReyNum (JE), &
                     TAB,   m%ElOut%Gamma  (JE), &
                            JE= 1, m%ElOut%NumElOut )


   ELSE
      WRITE(Frmt(10:12), '(I3)') 13*m%ElOut%NumElOut + 3
      WRITE(p%UnElem,Frmt) TIME,                 &
                     TAB,   m%ElOut%VXSAV,       &
                     TAB,   m%ElOut%VYSAV,       &
                     TAB,   m%ElOut%VZSAV,       &
                   ( TAB,   m%ElOut%ALF    (JE), &
                     TAB,   m%ElOut%DynPres(JE), &
                     TAB,   m%ElOut%CLL    (JE), &
                     TAB,   m%ElOut%CDD    (JE), &
                     TAB,   m%ElOut%CNN    (JE), &
                     TAB,   m%ElOut%CTT    (JE), &
                     TAB,   m%ElOut%PITSAV (JE), &
                     TAB,   m%ElOut%AAA    (JE), &
                     TAB,   m%ElOut%AAP    (JE), &
                     TAB,   m%ElOut%DFNSAV (JE), &
                     TAB,   m%ElOut%DFTSAV (JE), &
                     TAB,   m%ElOut%ReyNum (JE), &
                     TAB,   m%ElOut%Gamma  (JE), &
                            JE= 1, m%ElOut%NumElOut )
   ENDIF ! PMOMENT

IF (m%ElOut%NumWndElOut > 0) THEN

   WRITE(Frmt(10:12), '(I3)') 3*m%ElOut%NumWndElOut*p%NumBl
   WRITE(p%UnWndOut,Frmt) TIME,                   &
                    ( (  TAB, m%ElOut%SaveVX( JE, JB ), &
                         TAB, m%ElOut%SaveVY( JE, JB ), &
                         TAB, m%ElOut%SaveVZ( JE, JB ), &
                              JE = 1,m%ElOut%NumWndElOut ), JB = 1,p%NumBl )
ENDIF

ENDIF ! ELEMPRN



RETURN
END SUBROUTINE ElemOut
!=======================================================================

END MODULE AeroGenSubs
