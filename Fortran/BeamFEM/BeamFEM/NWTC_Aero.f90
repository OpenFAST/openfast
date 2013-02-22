MODULE NWTC_Aero


   ! This module contains aerodynamics routines with non-system-specific logic and references.


   ! It contains the following routines:

   !     SUBROUTINE AeroInt  ( ISeg, Alpha, Re, AF_Table, IntData, DoCl, DoCd, DoCm, DoCpmin [, ErrStat] )
   !     SUBROUTINE CompDR   ( NumSeg, RLoc, HubRad, RotorRad, DimenInp, DelRLoc [, ErrStat] )
   !     SUBROUTINE GetAF    ( AF_File, AF_Table, ISeg )
   !     FUNCTION   GetCoef  ( ISeg, Alpha, AlfaTab, CoefTab, NumRows, Ind [, ErrStat] )
   !     SUBROUTINE GetCoefs ( ISeg, Alpha, Re, AF_Table, ClInt, CdInt, CmInt, CpminInt, DoCl, DoCd, DoCm, DoCpmin [, ErrStat] )


   USE                             NWTC_IO
   USE                             NWTC_Num
   
   IMPLICIT  NONE


!=======================================================================


      ! Global aerodynamics-related variables.

   TYPE                            :: AeroData                                  ! Declare new type that holds the interpolated aero data from the big tables.
      REAL(ReKi)                   :: AlfaStal                                  ! The stall AoA for this table.
      REAL(ReKi)                   :: AOD                                       ! The AoA for minimum CD.
      REAL(ReKi)                   :: AOL                                       ! The zero-lift AoA.
      REAL(ReKi)                   :: Cd0                                       ! The minimum Cd value.
      REAL(ReKi)                   :: CnA                                       ! The Cn slope for zero-lift.
      REAL(ReKi)                   :: CnS                                       ! The Cn at stall value for positive AoA.
      REAL(ReKi)                   :: CnSL                                      ! Cn at stall value for negative AoA.
      REAL(ReKi)                   :: Cl                                        ! The lift coefficient.
      REAL(ReKi)                   :: Cd                                        ! The drag coefficient.
      REAL(ReKi)                   :: Cm                                        ! The pitching-moment coefficient.
      REAL(ReKi)                   :: Cpmin                                     ! The minimum pressure coefficient.
      REAL(ReKi)                   :: FTB                                       ! The normal-coefficient divided by the Cn slope at zero lift.
      REAL(ReKi)                   :: FTBC                                      ! The chordwise-coefficient divided by the Cn slope at zero lift.
   ENDTYPE AeroData

   TYPE                            :: AeroTable                                 ! Declare new type that is an allocatable table of data.
      REAL(ReKi)                   :: AlfaStal                                  ! The stall AoA for this table.
      REAL(ReKi)                   :: AOD                                       ! The AoA for minimum CD.
      REAL(ReKi)                   :: AOL                                       ! The zero-lift AoA.
      REAL(ReKi)                   :: Cd0                                       ! The minimum Cd value.
      REAL(ReKi)                   :: CnA                                       ! The Cn slope for zero-lift.
      REAL(ReKi)                   :: CnS                                       ! The Cn at stall value for positive AoA.
      REAL(ReKi)                   :: CnSL                                      ! Cn at stall value for negative AoA.
      REAL(ReKi)                   :: Re                                        ! The Re for this table.
      REAL(ReKi)                   :: Ctrl                                      ! The control setting for this table.
      INTEGER                      :: Ind      = 0                              ! Last-used index into table.  Zero at beginning.
      INTEGER                      :: NumAlf                                    ! Number of angles of attack in the table.
      REAL(ReKi), ALLOCATABLE      :: Alpha    (:)                              ! The angle of attack vector.
      REAL(ReKi), ALLOCATABLE      :: Cl       (:)                              ! The lift-coefficient vector.
      REAL(ReKi), ALLOCATABLE      :: Cd       (:)                              ! The drag-coefficient vector.
      REAL(ReKi), ALLOCATABLE      :: Cm       (:)                              ! The pitching-moment-coefficient vector.
      REAL(ReKi), ALLOCATABLE      :: Cpmin    (:)                              ! The minimum-pressure-coefficient vector.
      REAL(ReKi), ALLOCATABLE      :: FTB      (:)                              ! The normal-coefficient divided by the Cn slope at zero lift.
      REAL(ReKi), ALLOCATABLE      :: FTBC     (:)                              ! The chordwise-coefficient divided by the Cn slope at zero lift.
   ENDTYPE AeroTable

   TYPE                            :: AlfIndx                                   ! Declare new type that is an allocatable table of alpha indices.
      INTEGER                      :: NumBld                                    ! Number of blades in the table.
      INTEGER                      :: NumElm                                    ! Number of segments in the table.
      INTEGER, ALLOCATABLE         :: Ind      (:,:)                            ! The tables in this supertable.
   ENDTYPE AlfIndx

   TYPE                            :: ElmTable                                  ! Declare new type that is an allocatable table of data.
      INTEGER                      :: NumTabs                                   ! Number of tables in the supertable for an element.
      TYPE(AeroTable), ALLOCATABLE :: Tab      (:)                              ! The tables in this supertable.
   ENDTYPE ElmTable

   LOGICAL                         :: UseCm    = .FALSE.                        ! Flag to tell if there are Cm data in the airfoil files.
   LOGICAL                         :: UseCpmin = .FALSE.                        ! Flag to tell if there are Cp,min data in the airfoil files.


CONTAINS

!=======================================================================
   SUBROUTINE AeroInt ( ISeg, Alpha, Re, AF_Table, IntData, DoCl, DoCd, DoCm, DoCpmin, ErrStat )

      ! This routine finds the Re-bounding tables and then calls GetCoef() to get the
      ! desired coefficients for the two tables and then interpolates between them.

!NOTE: This routine needs to be modified to account for various control settings.  mlb  1-May-2010

      ! Argument declarations.

   REAL(ReKi), INTENT(IN)            :: Alpha                                   ! Angle of attack to get the coefficient for.
   REAL(ReKi), INTENT(IN)            :: Re                                      ! Reynolds number.

   INTEGER, INTENT(OUT), OPTIONAL    :: ErrStat                                 ! Error status; if present, program does not abort on error
   INTEGER, INTENT(IN)               :: ISeg                                    ! The current segment.

   LOGICAL, INTENT(IN)               :: DoCd                                    ! Get Cd.
   LOGICAL, INTENT(IN)               :: DoCl                                    ! Get Cl.
   LOGICAL, INTENT(IN)               :: DoCm                                    ! Get Cm.
   LOGICAL, INTENT(IN)               :: DoCpmin                                 ! Get Cp,min.

   TYPE (ElmTable), INTENT(INOUT)    :: AF_Table                                ! The table of airfoil data for the current segment.
   TYPE (AeroData), INTENT(OUT)      :: IntData                                 ! The interpolated airfoil data for the current segment.


      ! Local declarations.

   REAL(ReKi)                        :: CdHi                                    ! The drag coefficient for the higher Re.
   REAL(ReKi)                        :: ClHi                                    ! The lift coefficient for the higher Re.
   REAL(ReKi)                        :: CmHi                                    ! The pitching-moment coefficient for the higher Re.
   REAL(ReKi)                        :: CpminHi                                 ! The minimum pressure coefficient for the higher Re.
   REAL(ReKi)                        :: Fract                                   ! The fractional distance between tables.

   INTEGER                           :: ITab                                    ! An index for table number.
   INTEGER                           :: ITabLo                                  ! The table number that is the lower bound for Re.
   INTEGER                           :: ITabHi                                  ! The table number that is the lower bound for Re.

   LOGICAL                           :: OneTable                                ! Flag that tells if we need to read only one table (no interpolation).



      ! Find the bounding tables (if multiple) for this Re.  If there is only one table
      ! or if we are outside the range of tables, we won't need to interpolate.

   IF ( Re <= AF_Table%Tab(1)%Re )  THEN
      ITabLo   = 1
      OneTable = .TRUE.
   ELSE IF ( Re >= AF_Table%Tab(AF_Table%NumTabs)%Re )  THEN
      ITabLo   = AF_Table%NumTabs
      OneTable = .TRUE.
   ELSE IF ( AF_Table%NumTabs > 1 )  THEN
      DO ITab=1,AF_Table%NumTabs-1
         IF ( Re <= AF_Table%Tab(ITab+1)%Re )  THEN
            ITabLo = ITab
            ITabHi = ITab + 1
            EXIT
         END IF
      END DO
      OneTable = .FALSE.
   ELSE
      ITabLo   = 1
      OneTable = .TRUE.
   END IF


      ! Get the coefficients for ITabLo.

   IF ( DoCl )  THEN
      IF ( PRESENT( ErrStat ) ) THEN
         IntData%Cl = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cl, AF_Table%Tab(ITabLo)%NumAlf &
                           , AF_Table%Tab(ITabLo)%Ind, ErrStat )
         IF (ErrStat /= 0) RETURN
      ELSE
         IntData%Cl = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cl, AF_Table%Tab(ITabLo)%NumAlf &
                           , AF_Table%Tab(ITabLo)%Ind )
      END IF
   END IF

   IF ( DoCd )  THEN
      IF ( PRESENT( ErrStat ) ) THEN
         IntData%Cd  = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cd, AF_Table%Tab(ITabLo)%NumAlf &
                              , AF_Table%Tab(ITabLo)%Ind, ErrStat )
         IF (ErrStat /= 0) RETURN
      ELSE
         IntData%Cd  = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cd, AF_Table%Tab(ITabLo)%NumAlf &
                              , AF_Table%Tab(ITabLo)%Ind )
      END IF
   END IF

   IF ( DoCm )  THEN
      IF ( PRESENT( ErrStat ) ) THEN
         IntData%Cm = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cm, AF_Table%Tab(ITabLo)%NumAlf &
                             , AF_Table%Tab(ITabLo)%Ind, ErrStat )
         IF (ErrStat /= 0) RETURN
      ELSE
         IntData%Cm = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cm, AF_Table%Tab(ITabLo)%NumAlf &
                             , AF_Table%Tab(ITabLo)%Ind )
      END IF
   END IF

   IF ( DoCpmin )  THEN
      IF ( PRESENT( ErrStat ) ) THEN
         IntData%Cpmin = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cpmin, AF_Table%Tab(ITabLo)%NumAlf &
                                , AF_Table%Tab(ITabLo)%Ind, ErrStat )
         IF (ErrStat /= 0) RETURN
      ELSE
         IntData%Cpmin = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cpmin, AF_Table%Tab(ITabLo)%NumAlf &
                                , AF_Table%Tab(ITabLo)%Ind )
      END IF
   END IF

   IntData%AlfaStal = AF_Table%Tab(ITabLo)%AlfaStal
   IntData%AOD      = AF_Table%Tab(ITabLo)%AOD
   IntData%AOL      = AF_Table%Tab(ITabLo)%AOL
   IntData%Cd0      = AF_Table%Tab(ITabLo)%Cd0
   IntData%CnA      = AF_Table%Tab(ITabLo)%CnA
   IntData%CnS      = AF_Table%Tab(ITabLo)%CnS
   IntData%CnSL     = AF_Table%Tab(ITabLo)%CnSL


      ! If we don't need to interpolate, we don't need to make a second call and we are done.

   IF ( OneTable )  RETURN


      ! Get the coefficients for ITabHi.  Use step-wise interpolation for all but the first coefficient called.

   Fract = ( Re - AF_Table%Tab(ITabLo)%Re )/( AF_Table%Tab(ITabHi)%Re - AF_Table%Tab(ITabLo)%Re )


   IF ( DoCl )  THEN
      IF ( PRESENT( ErrStat ) ) THEN
         ClHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cl, AF_Table%Tab(ITabHi)%NumAlf, &
                                      AF_Table%Tab(ITabHi)%Ind, ErrStat )
         IF (ErrStat /= 0) RETURN
      ELSE
         ClHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cl, AF_Table%Tab(ITabHi)%NumAlf, &
                                      AF_Table%Tab(ITabHi)%Ind )
      END IF
      IntData%Cl  = IntData%Cl + Fract*( ClHi - IntData%Cl )
   END IF

   IF ( DoCd )  THEN
      IF ( PRESENT( ErrStat ) ) THEN
         CdHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cd, AF_Table%Tab(ITabHi)%NumAlf, &
                                      AF_Table%Tab(ITabHi)%Ind, ErrStat )
         IF (ErrStat /= 0) RETURN
      ELSE
         CdHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cd, AF_Table%Tab(ITabHi)%NumAlf, &
                                      AF_Table%Tab(ITabHi)%Ind )
      END IF
      IntData%Cd = IntData%Cd + Fract*( CdHi - IntData%Cd )
   END IF

   IF ( DoCm )  THEN
      IF ( PRESENT( ErrStat ) ) THEN
         CmHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cm, AF_Table%Tab(ITabHi)%NumAlf, &
                                      AF_Table%Tab(ITabHi)%Ind, ErrStat )
         IF (ErrStat /= 0) RETURN
      ELSE
         CmHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cm, AF_Table%Tab(ITabHi)%NumAlf, &
                                      AF_Table%Tab(ITabHi)%Ind )
      END IF
      IntData%Cm = IntData%Cm + Fract*( CmHi - IntData%Cm )
   END IF

   IF ( DoCpmin )  THEN
      IF ( PRESENT( ErrStat ) ) THEN
         CpminHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cpmin, AF_Table%Tab(ITabHi)%NumAlf, &
                                         AF_Table%Tab(ITabHi)%Ind, ErrStat )
         IF (ErrStat /= 0) RETURN
      ELSE
         CpminHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cpmin, AF_Table%Tab(ITabHi)%NumAlf, &
                                         AF_Table%Tab(ITabHi)%Ind )
      END IF
      IntData%Cpmin = IntData%Cpmin + Fract*( CpminHi - IntData%Cpmin )
   END IF

   IntData%AlfaStal = IntData%AlfaStal + Fract*( AF_Table%Tab(ITabHi)%AlfaStal - IntData%AlfaStal )
   IntData%AOD      = IntData%AOD      + Fract*( AF_Table%Tab(ITabHi)%AOD      - IntData%AOD      )
   IntData%AOL      = IntData%AOL      + Fract*( AF_Table%Tab(ITabHi)%AOL      - IntData%AOL      )
   IntData%Cd0      = IntData%Cd0      + Fract*( AF_Table%Tab(ITabHi)%Cd0      - IntData%Cd0      )
   IntData%CnA      = IntData%CnA      + Fract*( AF_Table%Tab(ITabHi)%CnA      - IntData%CnA      )
   IntData%CnS      = IntData%CnS      + Fract*( AF_Table%Tab(ITabHi)%CnS      - IntData%CnS      )
   IntData%CnSL     = IntData%CnSL     + Fract*( AF_Table%Tab(ITabHi)%CnSL     - IntData%CnSL     )


   RETURN
   END SUBROUTINE AeroInt ! ( ISeg, Alpha, Re, AF_Table, IntData, ClInt, CdInt, CmInt )
!=======================================================================
   SUBROUTINE CompDR ( NumSeg, RLoc, HubRad, RotorRad, DimenInp, DelRLoc, ErrStat )


      ! This routine computes the segment lengths from the local radii and the rotor radius.
      ! It prints and error if the list of radii is not realizable.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: NumSeg                                       ! Number of blade segments.
   INTEGER, INTENT(OUT),OPTIONAL:: ErrStat                                      ! Error status; if present, program does not abort on error

   REAL(ReKi), INTENT(OUT)      :: DelRLoc (NumSeg)                             ! The array of segment lengths.
   REAL(ReKi), INTENT(IN)       :: HubRad                                       ! The hub radius.
   REAL(ReKi), INTENT(IN)       :: RLoc    (NumSeg)                             ! The array of radii (segment centers).
   REAL(ReKi), INTENT(IN)       :: RotorRad                                     ! The rotor radius.

   LOGICAL, INTENT(IN)          :: DimenInp                                     ! Flag that tells if input is dimensional or not.


      ! Local declarations.

   REAL(ReKi)                   :: CompRad                                      ! The computed radius of the rotor.
   REAL(ReKi)                   :: ErrFact                                      ! The conversion to non-dimensional form if needed.
   REAL(ReKi)                   :: SegBeg                                       ! The beginning of the current segment.

   INTEGER                      :: ISeg                                         ! Segment index



   IF ( PRESENT(ErrStat) ) ErrStat = 0


      ! Determine the correct units for error messages.

   IF ( DimenInp )  THEN
      ErrFact = 1.0
   ELSE
      ErrFact = RotorRad
   END IF


      ! We will work our way from the root to the tip.

   SegBeg = HubRad

   DO ISeg=1,NumSeg

      IF ( RLoc(ISeg) <= SegBeg )  THEN
         ! START v1.03.00b-dcm  21-Jun-2010  D. Maniaci
         CALL WrScr1('Your analysis nodes are incorrectly defined.  Please review the following forum topic for an explaination' &
                    //' of this error:')
         CALL WrScr1('  https://wind.nrel.gov/forum/wind/viewtopic.php?f=4&t=241')
         ! END v1.03.00b-dcm  21-Jun-2010  D. Maniaci
         CALL ProgAbort ( ' The radius for blade segment #'//Trim( Int2LStr( ISeg ) )//' is too far inboard for a physically' &
                    //' realizable blade.  It must be greater than '//Trim( Num2LStr( SegBeg/ErrFact ) )//'.', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = 1
         RETURN
      END IF

      DelRLoc(ISeg) = 2.0*( RLoc(ISeg) - SegBeg )
      SegBeg        = SegBeg + DelRLoc(ISeg)

   END DO ! ISeg


      ! Ensure that the segments (almost) exactly fill the blade.

   CompRad = RLoc(NumSeg) + 0.5*DelRLoc(NumSeg)

   IF ( ABS( CompRad - RotorRad )/RotorRad > 0.005 )  THEN
         ! START v1.03.00b-dcm  21-Jun-2010  D. Maniaci
         CALL WrScr1('Your analysis nodes are incorrectly defined.  Please review the following forum topic for an explaination' &
                    //' of this error:')
         CALL WrScr1('  https://wind.nrel.gov/forum/wind/viewtopic.php?f=4&t=241')
         ! END v1.03.00b-dcm  21-Jun-2010  D. Maniaci
         CALL ProgAbort ( ' The sum of the lengths of the blade segments does not match the rotor radius.  The segments add up' &
                    //' to a rotor radius of '//Trim( Num2LStr( CompRad ) )//' instead of the specified radius of ' &
                    //Trim( Num2LStr( RotorRad ) )//'.  They must agree within 0.5%', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = 1
         RETURN
   ELSE IF ( ABS( CompRad - RotorRad )/RotorRad > 0.001 )  THEN
! Nice message, Marshall! ;-)  Thank you!  :-)  I don't even remember writing this.
         ! START v1.03.00b-dcm  21-Jun-2010  D. Maniaci
         CALL WrScr1('Your analysis nodes are incorrectly defined.  Please review the following forum topic for an explaination' &
                    //' of this error:')
         CALL WrScr1('  https://wind.nrel.gov/forum/wind/viewtopic.php?f=4&t=241')
         ! END v1.03.00b-dcm  21-Jun-2010  D. Maniaci
         CALL WrScr1 ( ' The sum of the lengths of the blade segments does not match the rotor radius.  The segments add up to a' &
                    //' rotor radius of '//Trim( Num2LStr( CompRad ) )//' instead of the specified radius of ' &
                    //Trim( Num2LStr( RotorRad ) )//'.  They really should agree within 0.1%, but I''ll let you slide.' )
!      IF ( Beep ) &
         CALL UsrAlarm
   END IF


   RETURN
   END SUBROUTINE CompDR ! ( NumSeg, RLoc, RotorRad, DimenInp, DelRLoc [, ErrStat] )
!=======================================================================
   SUBROUTINE GetAF ( AF_File, AF_Table, ISeg )
!bjj: note that this routine aborts instead of allowing an optional returned error code.

      ! Routine to get airfoil data from either a new NWTC-style or an old AeroDyn-style airfoil file.


      ! Argument declarations.

   TYPE (ElmTable), INTENT(OUT) :: AF_Table                                  ! The table of airfoil data for the current segment.

   INTEGER, INTENT(IN)          :: ISeg                                      ! The segment number.

   CHARACTER(*), INTENT(IN)     :: AF_File                                   ! Name of file containing AeroDyn-style airfoil data.


      ! Local declarations.

      ! Because of what seems to be a compiler bug, we cannot dynamically allocate the data arrays for the new-style
      ! airfoil files.  We really need to do it for the old-style files because there is no limit on the number of points.

!   TYPE                            :: DataRowO                                  ! Declare new type that is an allocatable table of data using a linked list.
!      REAL(ReKi), ALLOCATABLE      :: Data      (:)
!      TYPE(DataRowO), POINTER      :: Next            => NULL()
!   ENDTYPE DataRowO

   REAL(ReKi)                      :: AF_Data   (5)                             ! The values from one line of airfol data.
   REAL(ReKi), ALLOCATABLE         :: AF_DataO  (:)                             ! The values from one line of airfol data.
   REAL(ReKi), ALLOCATABLE         :: RnAry     (:)                             ! The temporary array for Re.
   REAL(ReKi), ALLOCATABLE         :: ASAry     (:)                             ! The temporary array for Stall AoA.
   REAL(ReKi), ALLOCATABLE         :: AOLAry    (:)                             ! The temporary array for zero-lift AoA.
   REAL(ReKi)                      :: Cc                                        ! The chordwise force coefficient.
   REAL(ReKi)                      :: Cn                                        ! The normal force coefficient.
   REAL(ReKi), ALLOCATABLE         :: CnAAry    (:)                             ! The temporary array for Cn slope for zero lift.
   REAL(ReKi), ALLOCATABLE         :: CnSAry    (:)                             ! The temporary array for Cn at stall value for positive AoA.
   REAL(ReKi), ALLOCATABLE         :: CnSLAry   (:)                             ! The temporary array for Cn at stall value for negative AoA.
   REAL(ReKi), ALLOCATABLE         :: AODAry    (:)                             ! The temporary array for AoA for minimum Cd.
   REAL(ReKi), ALLOCATABLE         :: CDOAry    (:)                             ! The temporary array for minimum Cd value.

   INTEGER                         :: IAlf                                      ! A generic array index for angle of attack.
   INTEGER                         :: Ind                                       ! A generic array index.
   INTEGER                         :: IOS                                       ! The status of an I/O operation.
   INTEGER                         :: ITab                                      ! The table index.
   INTEGER                         :: NumAlf                                    ! The number of lines in an old-style airfoil table.
   INTEGER                         :: NumAlpha                                  ! The number of non--blank lines in an old-style airfoil table.
   INTEGER                         :: NumCoef                                   ! The number of coefficiants in an airfoil table.
   INTEGER                         :: NumVals                                   ! The total number of values on one line of airfoil data.
   INTEGER                         :: Sttus                                     ! The status returned from the allocation.
   INTEGER                         :: UnAF     = 20                             ! I/O unit number for the airfoil file.

   CHARACTER( 15)                  :: Frmt = "(1000(F11.4,:))"                  ! Output format for a line of airfoil data.
   CHARACTER(999)                  :: Line                                      ! A line of text.
   CHARACTER(  3)                  :: Line3                                     ! The first three characters of a line of text.



      ! Open the airfoil data file.

   CALL OpenFInpFile ( UnAF, AF_File )


      ! Read the header block of the airfoil file.  Look to see if this is a new-format file.

   READ (UnAF,'(A)',IOSTAT=IOS)  Line

   CALL CheckIOS ( IOS, AF_File, 'FirstHead', StrType )

   IF ( Echo )  THEN
      WRITE (UnEc,"(15X,A,T30,' - ',A,/,2X,A)")  'FirstHead', 'First line in the airfoil file.', TRIM( Line )
   END IF

   CALL Conv2UC  ( Line )

   IF ( Line(:21) == 'AERODYN AIRFOIL FILE.' )  THEN


         ! This is new style of AeroDyn file.

      CALL ReadCom  ( UnAF, AF_File, 'the first title' )
      CALL ReadCom  ( UnAF, AF_File, 'the second title' )
      CALL ReadIVar ( UnAF, AF_File, AF_Table%NumTabs, 'NumTabs', 'Number of airfoil tables for segment #' &
                                                                 //TRIM( Int2LStr( ISeg ) )//'.' )

      IF ( AF_Table%NumTabs < 1 )  CALL ProgAbort ( ' Number of tables in airfoil file, "'//TRIM( AF_File ) &
                                              //'", must be > 0 for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )


         ! Are we expecting Cm data in the file?  Allocate the temporary data array.

      IF ( UseCm )  THEN
         NumVals = 4
      ELSE
         NumVals = 3
      END IF


         ! Are we expecting Cp,min data in the file?  Allocate the temporary data array.

      IF ( UseCpmin )  THEN
         NumVals = NumVals + 1
      END IF


         ! Allocate the AF_Table of pointers for this element.

      ALLOCATE ( AF_Table%Tab(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the Tab subtable of the AF_Table of pointers for segment #' &
                    //TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF


         ! Read the NumTabs airfoil tables.

      DO ITab=1,AF_Table%NumTabs


            ! Read in the Table ID (Re), control setting, and stall parameters for this table.

         CALL ReadRVar ( UnAF, AF_File, AF_Table%Tab(ITab)%Re      , 'Re('      //TRIM( Int2LStr( ITab ) )//')'   &
                                                                   , 'Reynolds number for this airfoil table.'    )
         CALL ReadRVar ( UnAF, AF_File, AF_Table%Tab(ITab)%Ctrl    , 'Ctrl('    //TRIM( Int2LStr( ITab ) )//')'   &
                                                                   , 'Control setting for this airfoil table.'    )
         CALL ReadRVar ( UnAF, AF_File, AF_Table%Tab(ITab)%AlfaStal, 'AlfaStal('//TRIM( Int2LStr( ITab ) )//')'   &
                                                                   , 'stall AoA for this airfoil table.'          )
         CALL ReadRVar ( UnAF, AF_File, AF_Table%Tab(ITab)%AOL     , 'AOL('     //TRIM( Int2LStr( ITab ) )//')'   &
                                                                   , 'zero-lift AoA.'                             )
         CALL ReadRVar ( UnAF, AF_File, AF_Table%Tab(ITab)%CnA     , 'CnA('     //TRIM( Int2LStr( ITab ) )//')'   &
                                                                   , 'Cn slope for zero-lift.'                    )
         CALL ReadRVar ( UnAF, AF_File, AF_Table%Tab(ITab)%CnS     , 'CnS('     //TRIM( Int2LStr( ITab ) )//')'   &
                                                                   , 'Cn at stall value for positive AoA.'        )
         CALL ReadRVar ( UnAF, AF_File, AF_Table%Tab(ITab)%CnSL    , 'CnSL('    //TRIM( Int2LStr( ITab ) )//')'   &
                                                                   , 'Cn at stall value for negative AoA.'        )
         CALL ReadRVar ( UnAF, AF_File, AF_Table%Tab(ITab)%AOD     , 'AOD('     //TRIM( Int2LStr( ITab ) )//')'   &
                                                                   , 'AoA for minimum CD.'                        )
         CALL ReadRVar ( UnAF, AF_File, AF_Table%Tab(ITab)%Cd0     , 'Cd0('     //TRIM( Int2LStr( ITab ) )//')'   &
                                                                   , 'minimum Cd value.'                          )


            ! Convert to proper units.

         AF_Table%Tab(ITab)%AlfaStal = AF_Table%Tab(ITab)%AlfaStal*D2R
         AF_Table%Tab(ITab)%AOD      = AF_Table%Tab(ITab)%AOD     *D2R
         AF_Table%Tab(ITab)%AOL      = AF_Table%Tab(ITab)%AOL     *D2R
         AF_Table%Tab(ITab)%Re       = AF_Table%Tab(ITab)%Re      *1.0e6


            ! Find the length of this table.

         AF_Table%Tab(ITab)%NumAlf = 0

         DO

            READ (UnAF,'(A)',IOSTAT=IOS)  Line3

            IF ( IOS < 0 )  THEN
               CALL PremEOF ( AF_File , 'the "EOT" end-of-table mark for airfoil table #'//TRIM( Int2LStr( ITab ) ) &
                                      //' and segment #'//TRIM( Int2LStr( ISeg ) ) )
            ELSE IF ( IOS > 0 )  THEN
               CALL WrScr1 ( ' Invalid character input for file "'//TRIM( AF_File )//'.' )
               CALL ProgAbort  ( ' The error occurred while trying to read line #'//TRIM( Int2LStr( AF_Table%Tab(ITab)%NumAlf+1 ) )&
                           //' of airfoil table #'//TRIM( Int2LStr( ITab ) )//' for segment #'//TRIM( Int2LStr( ISeg ) )//'.' )
            END IF

            CALL Conv2UC ( Line3 )
            IF ( Line3 == 'EOT' )  EXIT
            AF_Table%Tab(ITab)%NumAlf = AF_Table%Tab(ITab)%NumAlf + 1

         END DO


            ! Rewind the file to the beginning of this table.

         DO IAlf=1,AF_Table%Tab(ITab)%NumAlf+1
            BACKSPACE UnAF
         END DO ! IAlf


            ! Let's allocate the permanent table.

         ALLOCATE ( AF_Table%Tab(ITab)%Alpha(AF_Table%Tab(ITab)%NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the Alpha subtable for segment #'//TRIM( Int2LStr( ISeg) ) &
                       //' and table #'//TRIM( Int2LStr( ITab) )//').' )
         END IF

         ALLOCATE ( AF_Table%Tab(ITab)%Cl(AF_Table%Tab(ITab)%NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the Cl subtable for segment #'//TRIM( Int2LStr( ISeg) ) &
                       //' and table #'//TRIM( Int2LStr( ITab) )//').' )
         END IF

         ALLOCATE ( AF_Table%Tab(ITab)%Cd(AF_Table%Tab(ITab)%NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the Cd subtable for segment #'//TRIM( Int2LStr( ISeg) ) &
                       //' and table #'//TRIM( Int2LStr( ITab) )//').' )
         END IF

         IF ( UseCm )  THEN
            ALLOCATE ( AF_Table%Tab(ITab)%Cm(AF_Table%Tab(ITab)%NumAlf) , STAT=Sttus )
            IF ( Sttus /= 0 )  THEN
               CALL ProgAbort ( ' Error allocating memory for the Cm subtable for segment #'//TRIM( Int2LStr( ISeg) ) &
                       //' and table #'//TRIM( Int2LStr( ITab) )//').' )
            END IF
         END IF

         IF ( UseCpmin )  THEN
            ALLOCATE ( AF_Table%Tab(ITab)%Cpmin(AF_Table%Tab(ITab)%NumAlf) , STAT=Sttus )
            IF ( Sttus /= 0 )  THEN
               CALL ProgAbort ( ' Error allocating memory for the Cpmin subtable for segment #'//TRIM( Int2LStr( ISeg) ) &
                       //' and table #'//TRIM( Int2LStr( ITab) )//').' )
            END IF
         END IF

         ALLOCATE ( AF_Table%Tab(ITab)%FTB(AF_Table%Tab(ITab)%NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the FTB subtable for segment #'//TRIM( Int2LStr( ISeg) ) &
                       //' and table #'//TRIM( Int2LStr( ITab) )//').' )
         END IF

         ALLOCATE ( AF_Table%Tab(ITab)%FTBC(AF_Table%Tab(ITab)%NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the FTBC subtable for segment #'//TRIM( Int2LStr( ISeg) ) &
                       //' and table #'//TRIM( Int2LStr( ITab) )//').' )
         END IF


            ! Read in the airfoil data for this table.

         DO IAlf=1,AF_Table%Tab(ITab)%NumAlf

            READ (UnAF,*,IOSTAT=IOS)  ( AF_Data(Ind), Ind=1,NumVals )

            CALL CheckIOS ( IOS, AF_File, 'AF_Data', NumType )

            IF ( Echo )  WRITE (UnEc,Frmt)  ( AF_Data(Ind), Ind=1,NumVals )

            AF_Table%Tab(ITab)%Alpha(IAlf) = AF_Data(1)

            AF_Table%Tab(ITab)%Cl(IAlf) = AF_Data(2)
            AF_Table%Tab(ITab)%Cd(IAlf) = AF_Data(3)

            IF ( ABS( AF_Table%Tab(ITab)%CnA ) .GT. 1.0e-6 )  THEN

               Cn = AF_Data(2)*COS( AF_Data(1) ) + ( AF_Data(3) - AF_Table%Tab(ITab)%Cd0 )*SIN( AF_Data(1) )
               Cc = AF_Data(2)*SIN( AF_Data(1) ) - ( AF_Data(3) - AF_Table%Tab(ITab)%Cd0 )*COS( AF_Data(1) )

               AF_Table%Tab(ITab)%FTB (IAlf) = Cn/AF_Table%Tab(ITab)%CnA
               AF_Table%Tab(ITab)%FTBC(IAlf) = Cc/AF_Table%Tab(ITab)%CnA

            ELSE

               AF_Table%Tab(ITab)%FTB (IAlf) = 1.0
               AF_Table%Tab(ITab)%FTBC(IAlf) = 1.0

            END IF

            IF ( UseCm )  THEN
               AF_Table%Tab(ITab)%Cm(IAlf) = AF_Data(4)
            END IF

            IF ( UseCpmin )  THEN
               AF_Table%Tab(ITab)%Cpmin(IAlf) = AF_Data(NumVals)
            END IF

         END DO ! IAlf


            ! Check AoA range.

         IF ( ( AF_Table%Tab(ITab)%Alpha(1)                         > -180.0 ) .OR. &
              ( AF_Table%Tab(ITab)%Alpha(AF_Table%Tab(ITab)%NumAlf) <  180.0 ) )  THEN
            CALL ProgAbort ( 'Angle of attack range for airfoil table #'//TRIM( Int2LStr( ITab ) )//' of segment #' &
                       //TRIM( Int2LStr( ISeg ) )//' must be from -180 to 180.' )
         END IF


            ! Skip this EOT mark.

         READ (UnAF,'()')

      END DO ! ITab

   ELSE


         ! This is old style of AeroDyn file.

      CALL ReadCom  ( UnAF, AF_File, 'the second title' )
      CALL ReadIVar ( UnAF, AF_File, AF_Table%NumTabs, 'NumTabs', 'Number of airfoil tables for segment #' &
                                                                //TRIM( Int2LStr( ISeg ) )//'.' )

      ALLOCATE ( RnAry(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the RnAry array for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF

      ALLOCATE ( ASAry(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the ASAry array for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF

      ALLOCATE ( AOLAry(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the AOLAry array for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF

      ALLOCATE ( CnAAry(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the CnAAry array for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF

      ALLOCATE ( CnSAry(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the CnSAry array for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF

      ALLOCATE ( CnSLAry(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the CnSLAry array for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF

      ALLOCATE ( AODAry(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the AODAry array for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF

      ALLOCATE ( CDOAry(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the CDOAry array for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF

      CALL ReadRAry ( UnAF, AF_File, RnAry  , AF_Table%NumTabs, 'RnAry'  , 'Reynolds number values for the airfoil tables.' )
      CALL ReadRAry ( UnAF, AF_File, ASAry  , AF_Table%NumTabs, 'ASAry'  , 'Stall AoA for this airfoil table.'              )
      CALL ReadCom  ( UnAF, AF_File,                                       'the unused first obsolete stall parameter'      )
      CALL ReadCom  ( UnAF, AF_File,                                       'the unused second obsolete stall parameter'     )
      CALL ReadCom  ( UnAF, AF_File,                                       'the unused third obsolete stall parameter'      )
      CALL ReadRAry ( UnAF, AF_File, AOLAry , AF_Table%NumTabs, 'AOLAry' , 'zero-lift AoA'                                  )
      CALL ReadRAry ( UnAF, AF_File, CnAAry , AF_Table%NumTabs, 'CnAAry' , 'Cn slope for zero lift'                         )
      CALL ReadRAry ( UnAF, AF_File, CnSAry , AF_Table%NumTabs, 'CnSAry' , 'Cn at stall value for positive AoA'             )
      CALL ReadRAry ( UnAF, AF_File, CnSLAry, AF_Table%NumTabs, 'CnSLAry', 'Cn at stall value for negative AoA'             )
      CALL ReadRAry ( UnAF, AF_File, AODAry , AF_Table%NumTabs, 'AODAry' , 'AoA for minimum Cd'                             )
      CALL ReadRAry ( UnAF, AF_File, CDOAry , AF_Table%NumTabs, 'CDOAry' , 'minimum Cd value'                               )


         ! Are we expecting Cm data in the file?  Allocate the temporary data array.

      IF ( UseCm )  THEN
         NumCoef = 3
      ELSE
         NumCoef = 2
      END IF

      NumVals = 1 + NumCoef*AF_Table%NumTabs

      ALLOCATE ( AF_DataO(NumVals) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the AF_DataO array for segment #'//TRIM( Int2LStr( ISeg ) )//' in GetAF.' )
      END IF


         ! Allocate the AF_Table of pointers for this element.

      ALLOCATE ( AF_Table%Tab(AF_Table%NumTabs) , STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the AF_Table of pointers in GetAF.' )
      END IF


            ! Find the length of this table.

      NumAlf = 0

      DO
         READ (UnAF,'()',IOSTAT=IOS)
         IF ( IOS < 0 )  EXIT
         NumAlf = NumAlf + 1
      END DO


         ! Rewind the file to the beginning of this table.

      DO IAlf=1,NumAlf+1
         BACKSPACE UnAF
      END DO ! IAlf


         ! Let's allocate the tables.

      DO ITab=1,AF_Table%NumTabs

         ALLOCATE ( AF_Table%Tab(ITab)%Alpha(NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the Alpha vector of airfoil table #'//TRIM( Int2LStr( ITab) ) &
                       //' for element #'//TRIM( Int2LStr( ISeg) )//'.' )
         END IF

         ALLOCATE ( AF_Table%Tab(ITab)%Cl(NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the Cl vector of airfoil table #'//TRIM( Int2LStr( ITab) ) &
                       //' for element #'//TRIM( Int2LStr( ISeg) )//'.' )
         END IF

         ALLOCATE ( AF_Table%Tab(ITab)%Cd(NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the Cd vector of airfoil table #'//TRIM( Int2LStr( ITab) ) &
                       //' for element #'//TRIM( Int2LStr( ISeg) )//'.' )
         END IF

         IF ( UseCm )  THEN
            ALLOCATE ( AF_Table%Tab(ITab)%Cm(NumAlf) , STAT=Sttus )
            IF ( Sttus /= 0 )  THEN
               CALL ProgAbort ( ' Error allocating memory for the Cm vector of airfoil table #'//TRIM( Int2LStr( ITab) ) &
                          //' for element #'//TRIM( Int2LStr( ISeg) )//'.' )
            END IF
         END IF

!         ALLOCATE ( AF_Table%Tab(ITab)%FTB(AF_Table%Tab(ITab)%NumAlf) , STAT=Sttus )
         ALLOCATE ( AF_Table%Tab(ITab)%FTB(NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the FTB subtable for segment #'//TRIM( Int2LStr( ISeg) ) &
                       //' and table #'//TRIM( Int2LStr( ITab) )//').' )
         END IF

         ALLOCATE ( AF_Table%Tab(ITab)%FTBC(NumAlf) , STAT=Sttus )
         IF ( Sttus /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the FTBC subtable for segment #'//TRIM( Int2LStr( ISeg) ) &
                      //' and table #' //TRIM( Int2LStr( ITab) )//').' )
         END IF

      END DO ! ITab


         ! Let's read the data this time through.

      NumAlpha = NumAlf

      DO IAlf=1,NumAlf


            ! Let's skip blank lines.  Decrement the number of alphas when we find them.

         READ (UnAF,'(A)')  Line

         IF ( LEN_TRIM( Line ) == 0 )  THEN
            NumAlpha = NumAlpha - 1
            CYCLE
         END IF


            ! Let's get the data from the non-blank line.

         READ (Line,*,IOSTAT=IOS)  ( AF_DataO(Ind), Ind=1,NumVals )

         CALL CheckIOS ( IOS, AF_File, 'AF_DataO', NumType )

         IF ( Echo )  THEN
            WRITE (UnEc,Frmt)  ( AF_DataO(Ind), Ind=1,NumVals )
         END IF


            ! Let's move this good data into permanent storage.

         DO ITab=1,AF_Table%NumTabs

            AF_Table%Tab(ITab)%Alpha(IAlf) = AF_DataO(1)
            AF_Table%Tab(ITab)%Cl   (IAlf) = AF_DataO(NumCoef*(ITab-1)+2)
            AF_Table%Tab(ITab)%Cd   (IAlf) = AF_DataO(NumCoef*(ITab-1)+3)

            IF ( UseCm )  THEN
               AF_Table%Tab(ITab)%Cm(IAlf) = AF_DataO(NumCoef*(ITab-1)+4)
            END IF

            IF ( ABS( AF_Table%Tab(ITab)%CnA ) .GT. 1.0e-6 )  THEN

               Cn = AF_Data(2)*COS( AF_Data(1) ) + ( AF_Data(3) - AF_Table%Tab(ITab)%Cd0 )*SIN( AF_Data(1) )
               Cc = AF_Data(2)*SIN( AF_Data(1) ) - ( AF_Data(3) - AF_Table%Tab(ITab)%Cd0 )*COS( AF_Data(1) )

               AF_Table%Tab(ITab)%FTB (IAlf) = Cn/AF_Table%Tab(ITab)%CnA
               AF_Table%Tab(ITab)%FTBC(IAlf) = Cc/AF_Table%Tab(ITab)%CnA

            ELSE

               AF_Table%Tab(ITab)%FTB (IAlf) = 1.0
               AF_Table%Tab(ITab)%FTBC(IAlf) = 1.0

            END IF

         END DO ! ITab

      END DO ! IAlf


         ! Check AoA range.  AoAs are the same for all tables in a given segment.

      IF ( ( AF_Table%Tab(1)%Alpha(1)      > -180.0 ) .OR. &
           ( AF_Table%Tab(1)%Alpha(NumAlf) <  180.0 ) )  THEN
!           ( AF_Table%Tab(1)%Alpha(AF_Table%Tab(ITab)%NumAlf) <  180.0 ) )  THEN
         CALL ProgAbort ( 'Angle of attack range of airfoil tables for segment #'//TRIM( Int2LStr( ISeg ) ) &
                    //' must be from -180 to 180.' )
      END IF


           ! Store the header data in the permanent structure.

      DO ITab=1,AF_Table%NumTabs
         AF_Table%Tab(ITab)%AlfaStal = ASAry  (ITab)*D2R
         AF_Table%Tab(ITab)%Re       = RnAry  (ITab)*1.0e6
         AF_Table%Tab(ITab)%AOD      = AODAry (ITab)*D2R
         AF_Table%Tab(ITab)%AOL      = AOLAry (ITab)*D2R
         AF_Table%Tab(ITab)%Cd0      = CDOAry (ITab)
         AF_Table%Tab(ITab)%CnA      = CnAAry (ITab)
         AF_Table%Tab(ITab)%CnS      = CnSAry (ITab)
         AF_Table%Tab(ITab)%CnSL     = CnSLAry(ITab)
         AF_Table%Tab(ITab)%NumAlf   = NumAlpha
      END DO ! ITab


         ! Deallocate the temporary Re array.

      DEALLOCATE ( RnAry, STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error deallocating memory for the RnAry array in GetAF.' )
      END IF


         ! Deallocate the temporary data array.

      DEALLOCATE ( AF_DataO, STAT=Sttus )
      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( ' Error deallocating memory for the AF_DataO array in GetAF.' )
      END IF

   END IF

   CLOSE ( UnAF )


   RETURN
   END SUBROUTINE GetAF !  ( AF_File, AF_Table, ISeg )
!=======================================================================
   FUNCTION GetCoef( ISeg, Alpha, AlfaTab, CoefTab, NumRows, Ind, ErrStat )


      ! Interpolation routine for airfoil section coefficients.


      ! Function declaration.

   REAL(ReKi)                        :: GetCoef                                 ! The value returned by this function.


      ! Argument declarations.

   INTEGER, INTENT(INOUT)            :: Ind                                     ! The starting/resulting index into the tables.
   INTEGER, INTENT(IN)               :: ISeg                                    ! The current segment.
   INTEGER, INTENT(IN)               :: NumRows                                 ! The length of the arrays.
   INTEGER, INTENT(OUT), OPTIONAL    :: ErrStat                                 ! Error status; if present, program does not abort on error

   REAL(ReKi), INTENT(IN)            :: AlfaTab   (NumRows)                     ! Table of AoAs.
   REAL(ReKi), INTENT(IN)            :: Alpha                                   ! Angle of attack to get the coefficient for.
   REAL(ReKi), INTENT(IN)            :: CoefTab   (NumRows)                     ! Table of coefficients.


   IF ( PRESENT(ErrStat) ) ErrStat = 1


      ! If Alpha is to the outside the table, the user needs to make up some data.  Warn the user and stop the program.

   IF ( ( Alpha < AlfaTab(1) ) .OR. ( AlfaTab(NumRows) < Alpha ) )  THEN

      CALL ProgAbort ( ' For segment '//TRIM( Int2LStr( ISeg ) )//', the current angle of attack ('//TRIM( Num2LStr( Alpha ) ) &
                 //' degrees) is outside the domain of your data table (' //TRIM( Num2LStr( AlfaTab(1) ) )//' to ' &
                 //TRIM( Num2LStr( AlfaTab(NumRows) ) )//' degrees).  Please extend your data table.', PRESENT(ErrStat) )
      IF ( PRESENT(ErrStat) ) ErrStat = 1
      RETURN

   END IF


      ! Alpha is in range.  Interpolate.  Use binary interpolation if this is the first time to access this table.

   IF ( Ind == 0 )  THEN
      GetCoef = InterpBin( Alpha, AlfaTab, CoefTab, Ind, NumRows )
   ELSE
      GetCoef = InterpStp( Alpha, AlfaTab, CoefTab, Ind, NumRows )
   END IF


   RETURN
   END FUNCTION GetCoef ! ( ISeg, Alpha, AlfaTab, CoefTab, NumRows, Ind )
!=======================================================================
   SUBROUTINE GetCoefs ( ISeg, Alpha, Re, AF_Table, ClInt, CdInt, CmInt, CpminInt, DoCl, DoCd, DoCm, DoCpmin, ErrStat )


      ! This routine finds the Re-bounding tables and then calls GetCoef() to get the
      ! desired coefficients for the two tables and then interpolates between them.


      ! Argument declarations.

   TYPE (ElmTable), INTENT(INOUT)    :: AF_Table                                ! The table of airfoil data for the current segment.

   INTEGER, INTENT(IN)               :: ISeg                                    ! The current segment.
   INTEGER, INTENT(OUT), OPTIONAL    :: ErrStat                                 ! Error status; if present, program does not abort on error

   LOGICAL, INTENT(IN)               :: DoCd                                    ! Get Cd.
   LOGICAL, INTENT(IN)               :: DoCl                                    ! Get Cl.
   LOGICAL, INTENT(IN)               :: DoCm                                    ! Get Cm.
   LOGICAL, INTENT(IN)               :: DoCpmin                                 ! Get Cp,min.

   REAL(ReKi), INTENT(IN)            :: Alpha                                   ! Angle of attack to get the coefficient for.
   REAL(ReKi), INTENT(OUT)           :: CdInt                                   ! Interpolated drag coefficient.
   REAL(ReKi), INTENT(OUT)           :: ClInt                                   ! Interpolated lift coefficient.
   REAL(ReKi), INTENT(OUT)           :: CmInt                                   ! Interpolated pitching-moment coefficient.
   REAL(ReKi), INTENT(OUT)           :: CpminInt                                ! Interpolated minimum-pressure coefficient.
   REAL(ReKi), INTENT(IN)            :: Re                                      ! Reynolds number.


      ! Local declarations.

   REAL(ReKi)                        :: CdHi                                    ! The drag coefficient for the higher Re.
   REAL(ReKi)                        :: ClHi                                    ! The lift coefficient for the higher Re.
   REAL(ReKi)                        :: CmHi                                    ! The pitching-moment coefficient for the higher Re.
   REAL(ReKi)                        :: CpminHi                                 ! The minimum-pressure coefficient for the higher Re.
   REAL(ReKi)                        :: Fract                                   ! The fractional distance between tables.

   INTEGER                           :: ITab                                    ! An index for table number.
   INTEGER                           :: ITabLo                                  ! The table number that is the lower bound for Re.
   INTEGER                           :: ITabHi                                  ! The table number that is the lower bound for Re.

   LOGICAL                           :: OneTable                                ! Flag that tells if we need to read only one table (no interpolation).


   IF ( PRESENT(ErrStat) ) ErrStat = 0

      ! Find the bounding tables (if multiple) for this Re.  If there is only one table
      ! or if we are outside the range of tables, we won't need to interpolate.

   IF ( Re <= AF_Table%Tab(1)%Re )  THEN
      ITabLo   = 1
      OneTable = .TRUE.
   ELSE IF ( Re >= AF_Table%Tab(AF_Table%NumTabs)%Re )  THEN
      ITabLo   = AF_Table%NumTabs
      OneTable = .TRUE.
   ELSE IF ( AF_Table%NumTabs > 1 )  THEN
      DO ITab=1,AF_Table%NumTabs-1
         IF ( Re <= AF_Table%Tab(ITab+1)%Re )  THEN
            ITabLo = ITab
            ITabHi = ITab + 1
            EXIT
         END IF
      END DO
      OneTable = .FALSE.
   ELSE
      ITabLo   = 1
      OneTable = .TRUE.
   END IF


      ! Get the coefficients for ITabLo.


   IF ( PRESENT(ErrStat) ) THEN
      IF ( DoCl )  THEN
         ClInt = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cl, AF_Table%Tab(ITabLo)%NumAlf, &
                                       AF_Table%Tab(ITabLo)%Ind, ErrStat )
         IF (ErrStat > 0) RETURN
      END IF

      IF ( DoCd )  THEN
         CdInt = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cd, AF_Table%Tab(ITabLo)%NumAlf, &
                                       AF_Table%Tab(ITabLo)%Ind, ErrStat )
         IF (ErrStat > 0) RETURN
      END IF

      IF ( DoCm )  THEN
         CmInt = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cm, AF_Table%Tab(ITabLo)%NumAlf, &
                                       AF_Table%Tab(ITabLo)%Ind, ErrStat )
         IF (ErrStat > 0) RETURN
      END IF

      IF ( DoCpmin )  THEN
         CpminInt = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cpmin, AF_Table%Tab(ITabLo)%NumAlf, &
                                          AF_Table%Tab(ITabLo)%Ind, ErrStat )
         IF (ErrStat > 0) RETURN
      END IF

   ELSE  ! Abort the program when errors are found

      IF ( DoCl )  THEN
         ClInt = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cl, AF_Table%Tab(ITabLo)%NumAlf, &
                                       AF_Table%Tab(ITabLo)%Ind )
      END IF

      IF ( DoCd )  THEN
         CdInt = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cd, AF_Table%Tab(ITabLo)%NumAlf, &
                                       AF_Table%Tab(ITabLo)%Ind )
      END IF

      IF ( DoCm )  THEN
         CmInt = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cm, AF_Table%Tab(ITabLo)%NumAlf, &
                                       AF_Table%Tab(ITabLo)%Ind )
      END IF

      IF ( DoCpmin )  THEN
         CpminInt = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabLo)%Alpha, AF_Table%Tab(ITabLo)%Cpmin, AF_Table%Tab(ITabLo)%NumAlf, &
                                          AF_Table%Tab(ITabLo)%Ind )
      END IF

   END IF


      ! If we don't need to interpolate, we don't need to make a second call and we are done.

   IF ( OneTable )  RETURN                 ! We probably shouldn't do this.  We should probably do a block IF.  mlb


      ! Get the coefficients for ITabHi.  Use step-wise interpolation for all but the first coefficient called.

   Fract = ( Re - AF_Table%Tab(ITabLo)%Re )/( AF_Table%Tab(ITabHi)%Re - AF_Table%Tab(ITabLo)%Re )

   IF ( DoCl )  THEN
      IF ( PRESENT(ErrStat) ) THEN
         ClHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cl, AF_Table%Tab(ITabHi)%NumAlf, &
                                       AF_Table%Tab(ITabHi)%Ind, ErrStat )
         IF (ErrStat > 0) RETURN
      ELSE
         ClHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cl, AF_Table%Tab(ITabHi)%NumAlf, &
                                       AF_Table%Tab(ITabHi)%Ind )
      END IF
      ClInt = ClInt + Fract*( ClHi - ClInt )
   END IF

   IF ( DoCd )  THEN
      IF ( PRESENT(ErrStat) ) THEN
         CdHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cd, AF_Table%Tab(ITabHi)%NumAlf, &
                                       AF_Table%Tab(ITabHi)%Ind, ErrStat )
         IF (ErrStat > 0) RETURN
      ELSE
         CdHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cd, AF_Table%Tab(ITabHi)%NumAlf, &
                                       AF_Table%Tab(ITabHi)%Ind )
      END IF
      CdInt = CdInt + Fract*( CdHi - CdInt )
   END IF

   IF ( DoCm )  THEN
      IF ( PRESENT(ErrStat) ) THEN
         CmHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cm, AF_Table%Tab(ITabHi)%NumAlf, &
                                       AF_Table%Tab(ITabHi)%Ind, ErrStat )
         IF (ErrStat > 0) RETURN
      ELSE
         CmHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cm, AF_Table%Tab(ITabHi)%NumAlf, &
                                       AF_Table%Tab(ITabHi)%Ind )
      END IF
      CmInt = CmInt + Fract*( CmHi - CmInt )
   END IF

   IF ( DoCpmin )  THEN
      IF ( PRESENT(ErrStat) ) THEN
         CpminHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cpmin, AF_Table%Tab(ITabHi)%NumAlf, &
                                          AF_Table%Tab(ITabHi)%Ind, ErrStat )
         IF (ErrStat > 0) RETURN
      ELSE
         CpminHi = GetCoef( ISeg, Alpha, AF_Table%Tab(ITabHi)%Alpha, AF_Table%Tab(ITabHi)%Cpmin, AF_Table%Tab(ITabHi)%NumAlf, &
                                          AF_Table%Tab(ITabHi)%Ind )
      END IF
      CpminInt = CpminInt + Fract*( CpminHi - CpminInt )
   END IF


   RETURN
   END SUBROUTINE GetCoefs ! ( ISeg, Alpha, Re, AF_Table, ClInt, CdInt, CmInt, CpminInt, DoCl, DoCd, DoCm, DoCpmin, ErrStat )
!=======================================================================

END MODULE NWTC_Aero
