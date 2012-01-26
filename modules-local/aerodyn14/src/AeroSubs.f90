MODULE AeroSubs

   USE                     NWTC_Library
   USE                     InflowWind
   USE                     SharedInflowDefns


    IMPLICIT                NONE


CONTAINS
! ************************************************
!  AeroDyn Subroutines for YawDyn, ADAMS,
!   SymDyn and FAST.
! ************************************************
!  UNIVERSITY OF UTAH, MECHANICAL ENGINEERING DEPARTMENT
!====================================================================================================
SUBROUTINE AD_GetInput(UnIn, AeroInFile, WindFileName, Title, ErrStat )
!====================================================================================================

   USE                           AeroGenSubs, ONLY: AllocArrays

   USE                           AD_IOParams
   USE                           AeroTime
   USE                           ElOutParams
   USE                           Airfoil
   USE                           Blade
   USE                           Element
   USE                           InducedVel
   USE                           Rotor
   USE                           Switch
   USE                           TwrProps
   USE                           Wind


   IMPLICIT                      NONE

      ! Passed Variables:
   INTEGER, INTENT(IN)       :: UnIn
   CHARACTER(*), INTENT(IN)  :: AeroInFile
   CHARACTER(*), INTENT(OUT) :: WindFileName
   CHARACTER(*), INTENT(OUT) :: Title           ! used for ADOUT() -- bjj: get rid of this when ADOUT is gone!!!!
   INTEGER, INTENT(OUT)      :: ErrStat

      ! Local Variables:

   INTEGER                    :: ElIndex
   INTEGER                    :: IELM
   INTEGER                    :: IndPrint
   INTEGER                    :: K

   CHARACTER(1024)            :: LINE
   CHARACTER(1024)            :: FilePath             ! The path name of the AeroDyn input file (so files listed in it can be defined relative to the main input file location)

   !-------------------------------------------------------------------------------------------------
   ! Open the AeroDyn input file
   !-------------------------------------------------------------------------------------------------
   CALL OpenFInpFile(UnIn, TRIM(AeroInFile), ErrStat)
   IF (ErrStat /= 0 ) RETURN

   CALL GetPath( AeroInFile, FilePath )

   !-------------------------------------------------------------------------------------------------
   ! If the echo file is open, write the header...
   !-------------------------------------------------------------------------------------------------
   IF ( Echo ) THEN      
      WRITE( UnEc, '(// A /)' ) 'AeroDyn input data from file "'//TRIM( AeroInFile )//'":'  
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Read the AeroDyn input file
   !-------------------------------------------------------------------------------------------------

      ! Read in the title line
   CALL ReadStr( UnIn, AeroInFile, Title, 'Title', 'File title', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   CALL WrScr( ' Heading of the AeroDyn input file: '//TRIM(Title) )


      ! Read in the units specifier  - REMOVE SOON!!!!!
   CALL ReadVar( UnIn, AeroInFile, LINE, 'Units', 'Units option', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   LINE = ADJUSTL( LINE )
   CALL Conv2UC(LINE(1:3))

   IF (LINE(1:2) /= 'SI') THEN
      CALL ProgWarn( 'English units are no longer allowed in AeroDyn. Please modify your input files to use the metric system.')
      ErrStat = 1
      RETURN
   END IF

   SIunit = .TRUE.


      ! Read in the stall model
   CALL ReadVar( UnIn, AeroInFile, LINE, 'DStall', 'Stall model', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   CALL Conv2UC(LINE(1:7))

   SELECT CASE ( TRIM(Line) )
      CASE ('STEADY')
         DStall = .FALSE.
      CASE ('BEDDOES')
         DStall = .TRUE.
      CASE DEFAULT
         CALL ProgWarn( ' Error: Expecting "STEADY" or "BEDDOES" stall model option.')
         ErrStat = 1
         RETURN
   END SELECT


      ! Read in the CM option
   CALL ReadVar( UnIn, AeroInFile, LINE, 'PMoment', 'Pitching moment option', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   CALL Conv2UC(LINE(1:6))

   SELECT CASE ( TRIM(Line) )
      CASE ('USE_CM')
         PMoment = .TRUE.
      CASE ('NO_CM')
         PMoment = .FALSE.
      CASE DEFAULT
         CALL ProgWarn( ' Error: Expecting "USE_CM" or "NO_CM" pitching moment option.')
         ErrStat = 1
         RETURN
   END SELECT


      ! Read in the inflow model option
   CALL ReadVar( UnIn, AeroInFile, LINE, 'DynInfl', 'Inflow model', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   CALL Conv2UC(LINE(1:7))

   SELECT CASE ( Line(1:5) )
      CASE ('EQUIL')
         DynInfl = .FALSE.
         DynInit = .FALSE.

         IF (Line(6:7) == 'DA') THEN
            EqAIDmult = 1.0
            EquilDA   = .TRUE.
            EquilDT   = .FALSE.
         ELSEIF (LINE(6:7) == 'DB') THEN
            EqAIDmult = 1.0
            EquilDA   = .TRUE.
            EquilDT   = .TRUE.
         ELSEIF (LINE(6:7) == 'DT') THEN
            EqAIDmult = 0.0
            EquilDA   = .FALSE.
            EquilDT   = .TRUE.
         ELSE
            EqAIDmult = 0.0
            EquilDA   = .FALSE.
            EquilDT   = .FALSE.
         ENDIF

      CASE ('DYNIN')
         DynInfl = .TRUE.
         DynInit = .TRUE.
      CASE DEFAULT
         CALL ProgWarn( ' Error: Expecting "EQUIL" or "DYNIN" inflow model option.')
         ErrStat = 1
         RETURN
   END SELECT


      ! Read in the wake model
   CALL ReadVar( UnIn, AeroInFile, LINE, 'Wake', 'Wake model', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   CALL Conv2UC(LINE(1:5))

   SELECT CASE ( TRIM(Line) )
      CASE ('NONE')
         Wake = .FALSE.
         Swirl = .FALSE.

         CALL ProgWarn( ' All wake calculations are turned off! This option is recommended only '// &
                        'in high winds or for debugging.' )
      CASE ('WAKE')
         Wake = .TRUE.
         Swirl = .FALSE.
      CASE ('SWIRL')
         Wake = .TRUE.
         Swirl = .TRUE.
      CASE DEFAULT
         CALL ProgWarn( ' Error: Expecting "NONE", "WAKE", or "SWIRL" wake model option.')
         ErrStat = 1
         RETURN
   END SELECT

      ! Read in the tolerance for the wake model
   CALL ReadVar( UnIn, AeroInFile, AToler, 'AToler', 'Induction factor tolerance', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

       ! Read in the tip-loss model for EQUIL inflow
   CALL ReadVar( UnIn, AeroInFile, LINE, 'TLoss', 'Tip-loss model', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   IF ( DynInfl ) THEN          ! Initialize these variables, though they shouldn't be used
      TLoss = .FALSE.
      GTech = .FALSE.
   ELSE
      CALL Conv2UC(LINE(1:5))

      SELECT CASE ( LINE(1:5) )
         CASE ('NONE ')
            TLoss = .FALSE.
            GTech = .FALSE.
         CASE ('PRAND')
            TLoss = .TRUE.
            GTech = .FALSE.
         CASE ('GTECH')
            TLoss = .TRUE.
            GTech = .TRUE.
         CASE DEFAULT
            CALL ProgWarn( ' Error: Expecting "NONE", "PRAND", or "GTECH" tip-loss model option.')
            ErrStat = 1
            RETURN
      END SELECT
   END IF


       ! Read in the hub-loss model for EQUIL inflow
   CALL ReadVar( UnIn, AeroInFile, LINE, 'HLoss', 'Hub-loss model', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   IF ( DynInfl ) THEN          ! Initialize these variables, though they shouldn't be used
      HLoss = .FALSE.
   ELSE
      CALL Conv2UC( LINE(1:5) )

      SELECT CASE ( LINE(1:5) )
         CASE ('NONE')
            HLoss = .FALSE.
         CASE ('PRAND')
            HLoss = .TRUE.
         CASE DEFAULT
            CALL ProgWarn( ' Error: Expecting "NONE" or "PRAND" hub-loss model option.')
            ErrStat = 1
            RETURN
      END SELECT
   END IF


      ! Read in the wind file name
   CALL ReadVar( UnIn, AeroInFile, WindFileName, 'WindFileName', 'Name of the wind file', ErrStat)
   IF ( ErrStat /= 0 ) RETURN
   IF ( PathIsRelative( WindFileName ) ) WindFileName = TRIM(FilePath)//TRIM(WindFileName)

      ! Read in the wind reference (hub) height above ground
   CALL ReadVar( UnIn, AeroInFile, HH, 'RefHt', 'Wind reference height', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

!bjj: this is a hack job to allow both the new tower influence and the old tower wake models to be used
   CALL ReadStr( UnIn, AeroInFile, Line, 'NewTowerModel?', 'Check for tower influence model', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

      ! Check if this is the "special string" to indicate the new tower influence model
   CALL Conv2UC( Line )
   IF ( INDEX(Line, "NEWTOWER" ) > 0 ) THEN

      !----------------------------------------------------------------------------------------------
      ! New tower influence model, as implemented by PJM
      !----------------------------------------------------------------------------------------------
      PJM_Version = .TRUE.


         ! Read in the tower potential flow switch
      CALL ReadVar( UnIn, AeroInFile, TwrPotent, 'TwrPotent', 'Tower influence model', ErrStat)
      IF ( ErrStat /= 0 ) RETURN

         ! Read in the tower shadow switch
      CALL ReadVar( UnIn, AeroInFile, TwrShadow, 'TwrShadow', 'Tower shadow model', ErrStat)
      IF ( ErrStat /= 0 ) RETURN

         ! Read in the tower drag file name
      CALL ReadVar( UnIn, AeroInFile, TwrFile, 'TwrFile', 'Tower data file name', ErrStat)
      IF ( ErrStat /= 0 ) RETURN
      
      IF ( PathIsRelative( TwrFile ) ) TwrFile = TRIM(FilePath)//TRIM(TwrFile)

   ELSE
      !----------------------------------------------------------------------------------------------
      ! Old tower influence model, read TwrShad from Line for now
      !----------------------------------------------------------------------------------------------
      PJM_Version = .FALSE.

      TwrPotent = .FALSE.     ! We don't want to read the tower file!
      TwrShadow = .FALSE.     ! We don't want to read the tower file!

         ! Read in the tower shadow deficit
      IF ( INDEX( 'FTft', Line(:1) ) > 0 )  THEN
         CALL ProgWarn( ' Invalid numeric input. "'//TRIM( Line )//'" found when trying to read TwrShad.' )

         ErrStat = 1
         RETURN
      ELSE
         READ (Line,*,IOSTAT=ErrStat)  TwrShad
         CALL CheckIOS ( ErrStat, AeroInFile, 'TwrShad', NumType, .TRUE. )

         IF ( Echo )  THEN
            WRITE (UnEc,"( 2X, ES11.4e2, 2X, A, T30, ' - ', A )")  TwrShad, "TwrShad", 'Tower shadow deficit'
         END IF

      END IF
!----------------

      IF ( TwrShad >= 1.0 ) THEN
         CALL ProgWarn( ' Tower shadow deficit cannot be >= 1.  Setting default deficit = 0.3' )
         TwrShad = 0.3
      END IF


         ! Read in the tower shadow width
      CALL ReadVar( UnIn, AeroInFile, ShadHWid, 'ShadHWid', 'Tower shadow half width', ErrStat)
      IF ( ErrStat /= 0 ) RETURN

      IF ( ShadHWid <= 0.0 ) THEN
         CALL ProgWarn( ' Tower shadow width cannot be <= zero.  Setting default half width = 1.0' )
         ShadHWid = 1.0
      END IF


         ! Read in the tower shadow reference point (distance from yaw axis to hub)
      CALL ReadVar( UnIn, AeroInFile, T_Shad_Refpt, 'T_Shad_Refpt', 'Tower shadow reference point', ErrStat)
      IF ( ErrStat /= 0 ) RETURN

         ! Constants used in tower shadow calculations
      TShadC1 = ShadHWid / SQRT ( ABS( T_Shad_Refpt ) )
      TShadC2 = TwrShad  * SQRT ( ABS( T_Shad_Refpt ) )

   END IF


      ! Read in the air density
   CALL ReadVar( UnIn, AeroInFile, Rho, 'Rho', 'Air density', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   !bjj do we need to check the the air density is non-negative?

   IF ( Rho == 0.0 .AND. DynInfl ) THEN  ! Turn off the GDW if RHO = 0.  It will crash
      CALL ProgWarn( 'Air density is zero. Dynamic Inflow will be turned off to avoid program crash.' )
      DynInfl = .FALSE.
   ENDIF

      ! Read in the kinematic viscosity
   CALL ReadVar( UnIn, AeroInFile, KinVisc, 'KinVisc', 'Kinematic viscosity', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

      ! Aero calculation time interval
   CALL ReadVar( UnIn, AeroInFile, DtAero, 'DtAero', 'Aero calculation time step', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

      ! Read the number of airfoil files
   CALL ReadVar( UnIn, AeroInFile, NumFoil, 'NumFoil', 'Number of airfoil files', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   IF ( NumFoil < 1 ) THEN
      CALL ProgWarn( ' Error: Number of airfoil files must be a positive integer.')
      ErrStat = 1
      RETURN
   END IF

      !..............................................................................................
      ! Allocate space for the airfoil data file name(s), then read them
      !..............................................................................................
   IF (.NOT. ALLOCATED(FoilNm)) THEN
      ALLOCATE ( FoilNm( NumFoil ) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn(' Error allocating space for the FoilNm array.')
         RETURN
      END IF
   END IF

   CALL ReadAryLines( UnIn, AeroInFile, FoilNm(1:NumFoil), NumFoil, 'FoilNm', 'Airfoil file names', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   DO K=1,NumFoil
      IF ( PathIsRelative( FoilNm(K) ) ) FoilNm(K) = TRIM(FilePath)//TRIM( FoilNm(K) )
   END DO      


      ! Read in the number of blade elements
   CALL ReadVar( UnIn, AeroInFile, NElm, 'NElm', 'Number of blade elements', ErrStat)
   IF ( ErrStat /= 0 ) RETURN


      !..............................................................................................
      ! Allocate space for blade element data and read the arrays
      ! Read blade element data, check some inputs, convert twist to radians
      !..............................................................................................
   CALL AllocArrays ('Element')

   NumElOut    = 0 ! Initialize the element print array index
   NumWndElOut = 0

   CALL ReadCom( UnIn, AeroInFile, 'Element table headers', ErrStat)
   IF ( ErrStat /= 0 ) RETURN

   DO IElm = 1, NElm

      READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line to see if print/no print is enabled

      IF (ErrStat == 0) THEN
         READ(Line,*,IOSTAT=ErrStat) RElm(IElm), Twist(IElm), DR(IElm), C(IElm), NFoil(IElm)
      END IF

      IF ( ErrStat == 0 ) THEN

         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! Check if AeroDyn will print out the element and/or wind data
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

         CALL Conv2UC(LINE)

         ElPrList(IElm) = 0               ! INITIALIZE
         IndPrint = INDEX(LINE,"PRINT")
         IF (IndPrint > 0) THEN
            IF (LINE(IndPrint-2:IndPrint+4) /= "NOPRINT") THEN
               NumElOut = NumElOut + 1
               ElPrList(IElm) = NumElOut
            END IF
         ENDIF


         WndElPrList(IElm) = 0            ! INITIALIZE
         IndPrint = INDEX(LINE,"WIND")
         IF (IndPrint > 0) THEN
            IF (LINE(IndPrint-2:IndPrint-1) /= "NO") THEN
               NumWndElOut = NumWndElOut + 1
               WndElPrList(IElm) = NumWndElOut
            END IF
         ENDIF

         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! Echo data to the file NWTC_Library echo file, if requested
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

         IF ( Echo ) THEN     ! NWTC_Library echo file
            WRITE (UnEc,'(4(2X,ES11.4e2),2X,I11,1(2X,L11))')  RElm(IElm), Twist(IElm), DR(IElm), C(IElm), NFoil(IElm), &
                   ElPrList(IElm) /= 0  !,  WndElPrList(IElm) == 0
         END IF

      ELSE IF ( ErrStat < 0 ) THEN
         CALL ProgWarn(' Error reading from line '//TRIM(Int2Lstr(IElm))//' of the AeroDyn element table.' )

         CALL ProgWarn( ' Premature end of file while reading line '//TRIM(Int2Lstr(IElm))// &
                     ' of the AeroDyn element table in file "'//TRIM(AeroInFile)//'."' )
         RETURN
      ELSE
         CALL ProgWarn(' Error reading from line '//TRIM(Int2Lstr(IElm))// &
                       ' of the AeroDyn element table in file "'//TRIM(AeroInFile)//'."' )
         RETURN
      END IF


         ! Check for valid data:

      IF( C(IElm) <= 0.0 ) THEN
         CALL ProgWarn(' Error reading from line '//TRIM(Int2Lstr(IElm))//' of the AeroDyn element table.'// &
                       ' Chord length must be larger than 0.' )
         ErrStat = 1
         RETURN
      ENDIF

      IF (NFoil(IElm) < 1 .OR. NFoil(IElm) > NumFoil) THEN
         CALL ProgWarn(' Error reading from line '//TRIM(Int2Lstr(IElm))//' of the AeroDyn element table.'// &
                       ' Airfoil file ID must be a number between 1 and '//TRIM(Int2Lstr(NumFoil))//'.' )
         ErrStat = 1
         RETURN
      END IF

         ! Convert Twist to radians:

      Twist(IElm) = Twist(IElm)*D2R

   ENDDO ! IELM


      !..............................................................................................
      ! Read multiple airfoil table option
      !..............................................................................................
!bjj, this would be ideal, but it's annoying to read EOF for all these files...
!   CALL ReadVar( UnIn, AeroInFile, Line, 'MultiTab', 'Multiple airfoil table option', ErrStat)
!   IF ( ErrStat > 0 ) RETURN

      READ(UnIn,*,IOSTAT=ErrStat) Line         !read MultiTab -- it may not exist
      IF (ErrStat > 0 ) THEN
         CALL WrScr1 ( ' Invalid character input for file "'//TRIM( AeroInFile )//'".' )
         CALL ProgWarn ( ' The error occurred while trying to read "MultiTab".' )
         RETURN
      ELSE IF (ErrStat == 0) THEN
         IF ( Echo )  THEN
            WRITE (UnEc, "( 15X, A, T30, ' - ', A, /, 2X, A )" )  &
                         'MultiTab', 'Multiple airfoil table option', '"'//TRIM( Line )//'"'
         END IF
      ELSE
         MultiTab = .FALSE.
         Reynolds = .FALSE.
      !  CALL PremEOF ( TRIM( Fil ), Variable, TrapThisError )
      END IF

   !-------------------------------------------------------------------------------------------------
   ! Close AeroDyn input file
   !-------------------------------------------------------------------------------------------------
   CLOSE(UnIn)

   !-------------------------------------------------------------------------------------------------
   ! Read airfoil data and check for MultiTab values using LINE, which was read above
   !-------------------------------------------------------------------------------------------------
   CALL READFL()                                      ! Read airfoil files; bjj: make sure there isn't a conflict with ErrStat when this gets rewritten; also make it use the same UnIn

   MulTabLoc = 0.0                                    ! Initialize this value


            ! Read in the type of airfoil data table in each file
   IF ( ErrStat < 0 ) THEN             ! If we hit the end of the file without MultiTab, use only 1 airfoil table
      IF ( ANY( NTables(1:NumFoil) > 1 ) ) THEN
         CALL ProgWarn( ' Error reading multiple airfoil table option. Only one table for each file will be used.' )
      END IF
      MultiTab = .FALSE.
      Reynolds = .FALSE.
      ErrStat  = 0
   ELSE ! ErrStat = 0

      IF ( ANY( NTables(1:NumFoil) > 1 ) ) THEN
         CALL Conv2UC(LINE(1:6))

         SELECT CASE ( TRIM(Line) )
            CASE ( 'USER' )
               MultiTab = .TRUE.
               Reynolds = .FALSE.

               DO K = 1, NumFoil
                  IF ( NTables(K) > 1 ) THEN
                     IF ( ( MulTabLoc < MulTabMet(K,1) ) .OR. ( MulTabLoc > MulTabMet(K,NTables(K) ) ))THEN  !bjj: why don't we have this error elsewhere??? can't MulTabLoc change?
                        CALL ProgWarn( 'Error interpolating between airfoil tables. '// &
                                 ' Initial interpolation value = '//TRIM(Flt2LStr(MulTabLoc))// &
                                 ' is outside table range of '//TRIM(Flt2LStr(MulTabMet(K,1)))// &
                                 ' to '//TRIM(Flt2LStr(MulTabMet(K,NTables(K))))// &
                                 ' in airfoil file #'//TRIM(Int2LStr(K))//'.' )
                        ErrStat = 1
                        RETURN
                     END IF
                  END IF ! NTables(K) > 1
               ENDDO ! K

            CASE ( 'RENUM' )
               MultiTab = .TRUE.
               Reynolds = .TRUE.
            CASE ( 'SINGLE' )
               MultiTab = .FALSE.
               Reynolds = .FALSE.
            CASE DEFAULT
               CALL WrScr( ' Error: control model option must be "SINGLE" or "MULTI".' )
         END SELECT
      ELSE
         MultiTab = .FALSE.
         Reynolds = .FALSE.
      END IF

   ENDIF

   !-------------------------------------------------------------------------------------------------
   ! Read tower drag input file, if necessary
   !-------------------------------------------------------------------------------------------------
   IF (TwrPotent .OR. TwrShadow) THEN                 ! Read in the tower drag file
      CALL READTwr(UnIn, TwrFile, ErrStat)
      IF (ErrStat /= 0) RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Initialize variables for printing
   !-------------------------------------------------------------------------------------------------
   IF ( NumElOut > 0 .OR. NumWndElOut > 0 ) THEN
      ElemPrn = .TRUE.
      CALL AllocArrays ('ElPrint')     ! Allocate arrays with dimension NumElOut

      ElIndex = 0                      ! Re-Initialize the element print array index for wind
      DO IElm = 1, NElm
         IF (WndElPrList(IElm) > 0) THEN
            ElIndex = ElIndex + 1
            WndElPrNum(ElIndex) = IElm
         END IF
      END DO ! IELM

      ElIndex = 0 ! Re-Initialize the element print array index
      DO IElm = 1, NElm
         IF (ElPrList(IElm) > 0) THEN
            ElIndex = ElIndex + 1
            ElPrNum(ElIndex) = IElm
         END IF
      END DO ! IELM
   ELSE
      ElemPrn = .FALSE.
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Initialize BEDDOES dynamic stall data
   !-------------------------------------------------------------------------------------------------
   IF ( DStall ) CALL BEDDAT()


   RETURN

END SUBROUTINE AD_GetInput
!====================================================================================================
   SUBROUTINE ADOut ( TITLE, HubRad, WindFileName )
 !  used to output data pertinent
 !  to the yaw dynamics routines
 ! *****************************************************

USE                         AD_IOParams
USE                         ElOutParams
USE                         AeroTime
USE                         Airfoil
USE                         Blade
USE                         Element
USE                         InducedVel
USE                         Rotor
USE                         Switch
USE                         Wind

USE                         TwrProps



IMPLICIT                    NONE


   ! Passed Variables:

CHARACTER(*),  INTENT(IN):: TITLE
REAL(ReKi),    INTENT(IN):: HubRad
CHARACTER(*),  INTENT(IN):: WindFileName

   ! Local Variables:

INTEGER                  :: IElm
INTEGER                  :: IFoil

CHARACTER(  2)           :: Dst_Unit
CHARACTER(150)           :: Frmt
CHARACTER(  4)           :: Mass_Unit
CHARACTER( 35)           :: MESAGE
CHARACTER(  3)           :: Vel_Unit

CHARACTER(1)             :: Delim


IF (SIUNIT) THEN
   Dst_Unit  = 'm'
   Mass_Unit = 'kg'
   Vel_Unit  = 'mps'
ELSE
   Dst_Unit  = 'ft'
   Mass_Unit = 'slug'
   Vel_Unit  = 'fps'
ENDIF

 ! Reiterate the input file
WRITE(UnADopt,'(/A)') 'Inputs read in from the AeroDyn input file:'
WRITE(UnADopt,'(A/)') TRIM(TITLE)


MESAGE = 'Units for input and output'
IF ( SIUNIT ) THEN
   WRITE(UnADopt,'(A)') 'SI'//TAB//MESAGE
ELSE
   WRITE(UnADopt,'(A)') 'ENGLISH'//TAB//MESAGE
ENDIF

MESAGE = 'Dynamic stall model'
IF ( DSTALL ) THEN
   WRITE(UnADopt,'(A)') 'BEDDOES'//TAB//MESAGE//' [Beddoes]'
ELSE
   WRITE(UnADopt,'(A)') 'STEADY'//TAB//MESAGE//' [NO Dynamic stall]'
ENDIF

MESAGE = 'Aerodynamic pitching moment model'
IF ( PMOMENT ) THEN
   WRITE(UnADopt,'(A)') 'USE_CM'//TAB//MESAGE//' [Pitching Moments calculated]'
ELSE
   WRITE(UnADopt,'(A)') 'NO_CM'//TAB//MESAGE//' [NO Pitching Moments calculated]'
ENDIF

MESAGE = 'Inflow model'
IF ( DYNINFL ) THEN
   WRITE(UnADopt,'(A)') 'DYNIN'//TAB//MESAGE//' [Dynamic Inflow]'
ELSE
   IF ( EquilDA .AND. EquilDT )  THEN
      WRITE(UnADopt,'(A)') 'EQUILDB'//TAB//MESAGE//' [Equilibrium w/ axial and tangential drag]'
   ELSEIF ( EquilDA )  THEN
      WRITE(UnADopt,'(A)') 'EQUILDA'//TAB//MESAGE//' [Equilibrium w/ axial drag]'
   ELSEIF ( EquilDT )  THEN
      WRITE(UnADopt,'(A)') 'EQUILDA'//TAB//MESAGE//' [Equilibrium w/ tangential drag]'
   ELSE
      WRITE(UnADopt,'(A)') 'EQUIL'//TAB//MESAGE//' [Equilibrium]'
   ENDIF
ENDIF


MESAGE = 'Induction factor model'
IF ( WAKE ) THEN
   IF (SWIRL) THEN
      WRITE(UnADopt,'(A)') 'SWIRL'//TAB//MESAGE//' [Normal and Radial flow induction factors calculated]'
   ELSE
      WRITE(UnADopt,'(A)') 'WAKE'//TAB//MESAGE//' [Normal flow induction factors calculated]'
   ENDIF
   WRITE(UnADopt,'(A)') TRIM(Flt2LStr(ATOLER))//TAB// 'Convergence tolerance for induction factor'
ELSE
   WRITE(UnADopt,'(A)') 'NONE'//TAB//MESAGE//' [NO induction factors calculated]'
   WRITE(UnADopt,'(A)') '[Not Used]'//TAB//'Convergence tolerance for induction factor'
ENDIF

MESAGE = 'Tip-loss model'
IF (.NOT. DYNINFL) THEN
   IF ( TLOSS ) THEN
      IF (GTECH) THEN
         WRITE(UnADopt,'(A)') 'GTECH'//TAB//MESAGE//' [Georgia Tech correction to Prandtl model]'
      ELSE
         WRITE(UnADopt,'(A)') 'PRAND'//TAB//MESAGE//' [Prandtl model]'
      ENDIF
   ELSE
      WRITE(UnADopt,'(A)') 'NONE'//TAB//MESAGE//' [NO tip-loss calculated]'
   ENDIF
ELSE
      WRITE(UnADopt,'(A)') '[Not Used]'//TAB//MESAGE
ENDIF

MESAGE = 'Hub-loss model'
IF (.NOT. DYNINFL) THEN
   IF ( HLOSS ) THEN
      WRITE(UnADopt,'(A)') 'PRAND'//TAB//MESAGE//' [Prandtl model]'
   ELSE
      WRITE(UnADopt,'(A)') 'NONE'//TAB//MESAGE//' [NO hub-loss calculated]'
   ENDIF
ELSE
      WRITE(UnADopt,'(A)') '[Not Used]'//TAB//MESAGE
ENDIF


WRITE(UnADopt, '(A)') '"'//TRIM(WindFileName)//'"'//TAB//'  Wind file name'

WRITE(UnADopt,'(A)') TRIM(Flt2LStr(HH))//TAB// 'Wind reference (hub) height, '//TRIM(Dst_Unit)

IF ( PJM_Version ) THEN
   WRITE(UnADopt,'(L2, A)') TwrPotent, TAB//'Calculate tower potential flow [T or F]'
   WRITE(UnADopt,'(L2, A)') TwrShadow, TAB//'Calculate tower shadow [T or F]'
   IF ( TwrPotent .OR. TwrShadow ) THEN
      WRITE(UnADopt,'(A)') '"'//TRIM( TwrFile )//'"'//TAB//'Tower drag file name'
   ELSE
      WRITE(UnADopt,'(A)') '[none]'//TAB//'No tower drag properties file'
   ENDIF
ELSE
   WRITE(UnADopt,'(A)') TRIM(Flt2LStr(TwrShad))//TAB// 'Tower shadow centerline velocity deficit'
   WRITE(UnADopt,'(A)') TRIM(Flt2LStr(ShadHWid))//TAB// 'Tower shadow half width, '//TRIM(Dst_Unit)
   WRITE(UnADopt,'(A)') TRIM(Flt2LStr(T_Shad_Refpt))//TAB// 'Tower shadow reference point, '//TRIM(Dst_Unit)
END IF


WRITE(UnADopt,'(A)') TRIM(Flt2LStr(RHO))//TAB// 'Air density, '//TRIM(Mass_Unit)//'/'//TRIM(Dst_Unit)//'^3'
WRITE(UnADopt,'(A)') TRIM(Flt2LStr(KinVisc))//TAB// 'Kinematic air viscosity, '//TRIM(Dst_Unit)//'^2/sec'

WRITE(UnADopt,'(A)') TRIM(Flt2LStr(DTAERO))//TAB// 'Time interval for aerodynamic calculations, sec'

WRITE(UnADopt,'(A)') TRIM(Int2LStr(NUMFOIL))//TAB// 'Number of airfoil files used. Files listed below:'
DO IFoil = 1, NUMFOIL
   WRITE(UnADopt,'(A)') '"'//TRIM(FOILNM(IFoil))//'"'
END DO ! IFoil

WRITE(UnADopt,'(A)') TRIM(Int2LStr(NELM))//TAB// 'Number of blade elements per blade'

   Delim = ' ' !or Delim = TAB

   !-------------------------------------------------------------------------------------------------
   ! write out element information
   !-------------------------------------------------------------------------------------------------
   Frmt = '(3X,A10,8("'//Delim//'",A10))'

   WRITE(UnADopt,'( )')

      ! column names

   WRITE(UnADopt,Frmt) '  Element ', &
                       '   RELM   ', &
                       '   Twist  ', &
                       '    DR    ', &
                       '   Chord  ', &
                       '   NFoil  ', &
                       '  Print?  ', &
                       ' Tip-loss ', &
                       ' Hub-loss '

      ! column units

   WRITE(UnADopt,Frmt) '    (-)   ', &
                       '    (m)   ', &
                       '   (deg)  ', &
                       '    (m)   ', &
                       '    (m)   ', &
                       '    (-)   ', &
                       ' (Yes/No) ', &
                       ' constant ', &
                       ' constant '

   WRITE(UnADopt,Frmt) '----------', &
                       '----------', &
                       '----------', &
                       '----------', &
                       '----------', &
                       '----------', &
                       '----------', &
                       '----------', &
                       '----------'

      ! column data
   Frmt = '(3X, I10, 4("'//Delim//'",F10.5),"'//Delim//'",I10,"'//Delim//'",A10, 2("'//Delim//'",F10.5) )'

   DO IElm = 1, NELM

      IF (ElPrList(IElm) /= 0) THEN
         MESAGE = 'Yes'
      ELSE
         MESAGE = 'No'
      ENDIF

      WRITE(UnADopt, Frmt) IElm, RELM(IElm), TWIST(IElm)*R2D, DR(IElm),  C(IElm), &
                           NFOIL(IElm), TRIM(Mesage), TLCNST(IElm), HLCNST(IElm)
   END DO


IF ( MultiTab ) THEN
   WRITE(UnADopt,'(A)') 'MULTI    Multiple airfoil tables used'
ELSE
   WRITE(UnADopt,'( )')
ENDIF

      WRITE(UnADopt,"(/' Rotor radius     = ',F7.3,' m')") R
      WRITE(UnADopt,"( ' Hub radius       = ',F7.3,' m')") HubRad
      WRITE(UnADopt,"( ' Number of blades = ',I3       )") NB

IF ( DSTALL ) CALL BEDWRT

IF ( ELEMPRN ) WRITE(UnADopt,'(/A/)')'Blade element aerodynamic time series data written to file.' !bjj: what? isn't that what "Print? Y/N" is for?

CLOSE (UnADopt )


RETURN
END SUBROUTINE ADOut

 ! ****************************************************
   SUBROUTINE READFL
 !  Reads a data file containing airfoil angle of attack,
 !   CL and CD, and dynamic stall parameters
 ! ****************************************************

USE                             AD_IOParams
USE                             Airfoil
USE                             Bedoes
USE                             AeroGenSubs, ONLY: AllocArrays
USE                             Switch


IMPLICIT                        NONE


   ! Local Variables:

REAL(ReKi), ALLOCATABLE      :: CDNegPI(:)
REAL(ReKi), ALLOCATABLE      :: CDPosPI(:)
REAL(ReKi), ALLOCATABLE      :: CLNegPI(:)
REAL(ReKi), ALLOCATABLE      :: CLPosPI(:)
REAL(ReKi), ALLOCATABLE      :: CMNegPI(:)
REAL(ReKi), ALLOCATABLE      :: CMPosPI(:)

INTEGER                      :: IPHI
INTEGER                      :: I
INTEGER                      :: K
INTEGER                      :: Sttus
INTEGER                      :: NFOILID
INTEGER                      :: NumLines
INTEGER                      :: NUNIT
INTEGER                      :: IOS

LOGICAL                      :: ALPosPI
LOGICAL                      :: ALNegPI

CHARACTER( 40)               :: TITLE  (2)
CHARACTER(1024)              :: LINE


NUNIT    = UnAirfl
NumCL    = 0

 ! The first loop checks existence and file length to set NumCL
DO NFOILID = 1, NUMFOIL

 ! Open the file for reading # of lines
   CALL OpenFInpFile (NUNIT, TRIM(FOILNM(NFOILID)))

 ! Determine the maximum number of aerodata points in all files

   NumLines = 0
   IOS = 0
   DO WHILE (IOS == 0)
      READ ( NUNIT, '()', IOSTAT=IOS )
      NumLines = NumLines + 1
   END DO

   NumCL = MAX(NumLines - 14, NumCL)

   CLOSE (NUNIT)

END DO ! NFOILID

 ! Allocate the arrays

Call AllocArrays ('Aerodata')

 ! The second loop reads the files
DO NFOILID = 1, NUMFOIL

 ! Open the file for reading inputs
   CALL OpenFInpFile (NUNIT, TRIM(Adjustl(FOILNM(NFOILID))) )

 ! Set up the file to read the aerodata
   READ(NUNIT,'( A )') TITLE(1)
   READ(NUNIT,'( A )') TITLE(2)

 ! Read in airfoil table dimension parameters:
 !   NTables = number of airfoil data tables

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , '# of tables' )

   READ(LINE,*,ERR=205) NTables( NFOILID )

 ! Allocate local arrays with NTables dimension

   Sttus = 0
   IF (.NOT. ALLOCATED(CLPosPI)) ALLOCATE ( CLPosPI(NTables(NFOILID)) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CLPosPI array.' )

   IF (.NOT. ALLOCATED(CDPosPI)) ALLOCATE ( CDPosPI(NTables(NFOILID)) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CDPosPI array.' )

   IF (.NOT. ALLOCATED(CMPosPI)) ALLOCATE ( CMPosPI(NTables(NFOILID)) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CMPosPI array.' )

   IF (.NOT. ALLOCATED(CLNegPI)) ALLOCATE ( CLNegPI(NTables(NFOILID)) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CLNegPI array.' )

   IF (.NOT. ALLOCATED(CDNegPI)) ALLOCATE ( CDNegPI(NTables(NFOILID)) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CDNegPI array.' )

   IF (.NOT. ALLOCATED(CMNegPI)) ALLOCATE ( CMNegPI(NTables(NFOILID)) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for CMNegPI array.' )


 ! Read in airfoil data table identification array

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , 'multi-table metric' )
   READ(LINE,*,ERR=205)  (MulTabMet ( NFOILID, K ), K = 1, NTables(NFOILID))

 ! Read in four lines that are no longer used
 ! These are retained for future USE and backwards compatibility only

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , '5th line' )
   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , '6th line' )
   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , '7th line' )
   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , '8th line' )

 ! Read Beddoes stall parameters

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , 'Angle of zero lift (AOL)' )
   IF (DSTALL) READ(LINE,*,ERR=205)  (AOL( NFOILID, K ), K = 1, NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , 'CNA' )
   IF (DSTALL) READ(LINE,*,ERR=205)  (CNA   ( NFOILID, K ), K = 1, NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , 'CNS' )
   IF (DSTALL) READ(LINE,*,ERR=205)  (CNS   ( NFOILID, K ), K = 1, NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , 'CNSL' )
   IF (DSTALL) READ(LINE,*,ERR=205)  (CNSL  ( NFOILID, K ), K = 1, NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , 'AOD' )
   IF (DSTALL) READ(LINE,*,ERR=205)  (AOD   ( NFOILID, K ), K = 1, NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) CALL PremEOF ( Trim(FOILNM(NFOILID)) , 'CDO' )
   IF (DSTALL) READ(LINE,*,ERR=205)  (CDO   ( NFOILID, K ), K = 1, NTables(NFOILID))

 ! Convert angles to radians
   IF (DSTALL) THEN
      AOD   ( NFOILID, : ) = AOD( NFOILID, : )*D2R
      AOL   ( NFOILID, : ) = AOL( NFOILID, : )*D2R
   ENDIF


 ! Read airfoil data tables to end of file

   NLIFT(NFOILID) = 0
   ALPosPI = .FALSE.
   ALNegPI = .FALSE.

   DO I = 1, NumCL

      IF ( PMOMENT ) THEN

         READ( NUNIT,*,END=150 ) AL(NFOILID,I), &
             (CL(NFOILID,I,IPHI), CD(NFOILID,I,IPHI), &
              CM(NFOILID,I,IPHI), IPHI = 1, NTables(NFOILID))

      ELSE

         READ( NUNIT,*,END=150 ) AL(NFOILID,I), &
             (CL(NFOILID,I,IPHI), CD(NFOILID,I,IPHI), &
              IPHI = 1, NTables(NFOILID))

         CM(NFOILID,I,:) = 0.

      ENDIF

 ! Check to see if values look reasonable

      DO IPHI = 1, NTables(NFOILID)
        IF ( ABS( AL( NFOILID, I ) ) > 185.) THEN
           CALL ProgAbort( 'Probable error in airfoil data table number '//TRIM(Int2LStr(NFOILID))// &
                           ' Angle of attack exceeds 185 degrees.')
        ELSEIF (ABS( CL( NFOILID, I, IPHI ) ) > 3. ) THEN
           CALL ProgAbort( 'Probable error in airfoil data table number '//TRIM(Int2LStr(NFOILID))// &
                           ' Coefficient of Lift exceeds 3.0.')
        ELSEIF (ABS( CD( NFOILID, I, IPHI ) ) > 3. ) THEN
           CALL ProgAbort( 'Probable error in airfoil data table number '//TRIM(Int2LStr(NFOILID))// &
                           ' Coefficient of Drag exceeds 3.0.')
        ELSEIF (ABS( CM( NFOILID, I, IPHI ) ) > 3. ) THEN
           CALL ProgAbort( 'Probable error in airfoil data table number '//TRIM(Int2LStr(NFOILID))// &
                           ' Coefficient of Moment exceeds 3.0.')
        ENDIF
      ENDDO ! IPHI

 ! Store the values at 180 deg. and -180 deg. for check
      IF ( AL (NFOILID, I ) == 180. ) THEN
         ALPosPI = .TRUE.
         Do IPHI = 1, NTables(NFOILID)
            CLPosPI(IPHI) = CL(NFOILID,I,IPHI)
            CDPosPI(IPHI) = CD(NFOILID,I,IPHI)
            CMPosPI(IPHI) = CM(NFOILID,I,IPHI)
         END Do ! IPHI

      ELSEIF ( AL (NFOILID, I ) == -180. ) THEN
         ALNegPI = .TRUE.
         Do IPHI = 1, NTables(NFOILID)
            CLNegPI(IPHI) = CL(NFOILID,I,IPHI)
            CDNegPI(IPHI) = CD(NFOILID,I,IPHI)
            CMNegPI(IPHI) = CM(NFOILID,I,IPHI)
          END Do ! IPHI
      ENDIF

      AL ( NFOILID, I ) = AL(NFOILID,I) * D2R
      NLIFT ( NFOILID ) = NLIFT(NFOILID) + 1

   ENDDO ! I

   150 CLOSE( NUNIT )

 ! Check to see if values at 180 deg. equal those at -180 deg.
   IF (ALPosPI .AND. ALNegPI) THEN
      Do IPHI = 1, NTables(NFOILID)
         IF (CLPosPI(IPHI) /= CLNegPI(IPHI) .OR. &
             CDPosPI(IPHI) /= CDNegPI(IPHI) .OR. &
             CMPosPI(IPHI) /= CMNegPI(IPHI)) THEN
            CALL ProgAbort( ' The airfoil data at +180 deg is different from -180 deg in file :'//Trim(FOILNM(NFOILID)) )
         ENDIF
      END Do ! IPHI
   ENDIF

 ! Deallocate arrays to make them available for the next file

   IF ( ALLOCATED(CLPosPI) ) DEALLOCATE ( CLPosPI )
   IF ( ALLOCATED(CDPosPI) ) DEALLOCATE ( CDPosPI )
   IF ( ALLOCATED(CMPosPI) ) DEALLOCATE ( CMPosPI )
   IF ( ALLOCATED(CLNegPI) ) DEALLOCATE ( CLNegPI )
   IF ( ALLOCATED(CDNegPI) ) DEALLOCATE ( CDNegPI )
   IF ( ALLOCATED(CMNegPI) ) DEALLOCATE ( CMNegPI )


END DO !NUMFOIL
RETURN

!bjj FIX THIS GOTO!!!!!!!!
205 CALL ProgAbort( ' Error reading line: "'//TRIM(Line)//'" in file : "'//TRIM(FOILNM(NFOILID))//'"' )


RETURN
END SUBROUTINE READFL

!====================================================================================================
SUBROUTINE READTwr(UnIn, FilName, ErrStat)
! This subroutine reads the tower properties input file, allocating TwrProps variables to do so.
! The tower data file contains radius and Re vs CD data as well as the tower wake constant.
!====================================================================================================

   USE                             TwrProps

   IMPLICIT                        NONE

      ! Passed variables:

   INTEGER,      INTENT(IN)     :: UnIn      ! unit number for tower input file
   CHARACTER(*), INTENT(IN)     :: FilName   ! name of the tower input file
   INTEGER,      INTENT(OUT)    :: ErrStat   ! returns 0 if no errors were encountered; non-zero otherwise


      ! Local Variables:

   INTEGER                      :: I         ! loop counter for rows in the data tables
   INTEGER                      :: J         ! loop counter for columns in the data tables

   CHARACTER(99)                :: Fmt       ! format for printing to an echo file

   !-------------------------------------------------------------------------------------------------
   ! Open the file for reading
   !-------------------------------------------------------------------------------------------------
   CALL OpenFInpFile (UnIn, TRIM(FilName), ErrStat )
   IF ( ErrStat /= 0 ) RETURN


   !-------------------------------------------------------------------------------------------------
   ! Read the heading, section 1
   !-------------------------------------------------------------------------------------------------

      ! Read in 2 header/comment lines
   CALL ReadCom( UnIn, FilName, 'Title line 1', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   CALL ReadCom( UnIn, FilName, 'Title line 2', ErrStat )
   IF ( ErrStat /= 0 ) RETURN


      ! Read in number of tower height entries, NTwrHt
   CALL ReadVar( UnIn, FilName, NTwrHt, 'NTwrHt', 'Number of tower stations', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   IF (NTwrHt < 1) THEN
      CALL ProgWarn( 'Number of tower height entries, NTwrHt, must be greater than zero.' )
      ErrStat = 1
      RETURN
   ENDIF

      ! Read in number of tower Reynolds number entries, NTwrRe
   CALL ReadVar( UnIn, FilName, NTwrRe, 'NTwrRe', 'Number of tower Reynolds number rows', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   IF (NTwrRe < 1) THEN
      CALL ProgWarn( 'Number of tower Reynolds number entries, NTwrRe, must be greater than zero.' )
      ErrStat = 1
      RETURN
   ENDIF


      ! Read in number of tower CD entries, NTwrCD
   CALL ReadVar( UnIn, FilName, NTwrCD, 'NTwrCD', 'Number of tower CD columns', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   IF (NTwrCD < 1) THEN
      CALL ProgWarn( 'Number of tower CD entries, NTwrCD, must be greater than zero.' )
      ErrStat = 1
      RETURN
   ENDIF


      ! Read in constant for tower wake model = 0 full potential flow = 0.1 model of Bak et al.
   CALL ReadVar( UnIn, FilName, Tower_Wake_Constant, 'Tower_Wake_Constant', 'Constant for tower wake model', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   ! bjj: should there be a sanity check here, too?

   !-------------------------------------------------------------------------------------------------
   ! Allocate TwrProps arrays with NTwrHt, NTwrRe, and NTwrCD dimensions; these arrays are
   ! read in the next 2 sections of this file.
   !-------------------------------------------------------------------------------------------------

   IF ( .NOT. ALLOCATED( TwrHtFr ) ) THEN
      ALLOCATE ( TwrHtFr(NTwrHt) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for TwrHtFr array.' )
         ErrStat = 1
         RETURN
      END IF
   END IF

   IF ( .NOT. ALLOCATED( TwrWid ) ) THEN
      ALLOCATE ( TwrWid(NTwrHt) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for TwrWid array.' )
         ErrStat = 1
         RETURN
      END IF
   END IF

   IF ( .NOT. ALLOCATED( NTwrCDCol ) ) THEN
      ALLOCATE ( NTwrCDCol(NTwrHt) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for NTwrCDCol array.' )
         ErrStat = 1
         RETURN
      END IF
   END IF

   IF ( .NOT. ALLOCATED( TwrRe ) ) THEN
      ALLOCATE ( TwrRe(NTwrRe) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for TwrRe array.' )
         ErrStat = 1
         RETURN
      END IF
   END IF

   IF ( .NOT. ALLOCATED( TwrCD ) ) THEN
      ALLOCATE ( TwrCD(NTwrRe, NTwrCD) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for TwrCD array.' )
         ErrStat = 1
         RETURN
      END IF
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Read section 2, DISTRIBUTED TOWER PROPERTIES;
   ! section contains 2 heading lines in addition to NTwrHt rows of data with 3 columns
   !-------------------------------------------------------------------------------------------------
      ! Read in 2 header/comment lines
   CALL ReadCom( UnIn, FilName, 'Distributed Tower Properties header 1', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   CALL ReadCom( UnIn, FilName, 'Distributed Tower Properties header 2', ErrStat )
   IF ( ErrStat /= 0 ) RETURN


      ! Read tower radius data table

   DO I = 1, NTwrHt     ! 1 line per fraction of height


      READ( UnIn,*,IOSTAT=ErrStat ) TwrHtFr(I), TwrWid(I), NTwrCDCol(I)

      IF ( ErrStat == 0 ) THEN
         IF ( Echo ) THEN
            WRITE (UnEc,'(2X,ES11.4e2, 2X,ES11.4e2, 2X,I11)')  TwrHtFr(I), TwrWid(I), NTwrCDCol(I)
         END IF
      ELSE IF ( ErrStat < 0 ) THEN
         CALL ProgWarn( ' Premature end of file while reading line '//TRIM(Int2Lstr(I))// &
                     ' of the distributed tower properties in file "'//TRIM(FilName)//'."' )
         RETURN
      ELSE
         CALL ProgWarn( ' Error reading line '//TRIM(Int2Lstr(I))// &
                     ' of the distributed tower properties in file "'//TRIM(FilName)//'."' )
         RETURN
      END IF


      !..............................................................................................
      ! Check to see if values look reasonable
      !..............................................................................................

         ! Make sure tower height fractions are between 0 and 1
      IF ( TwrHtFr( I ) < 0.0 .OR. TwrHtFr( I ) > 1.0 ) THEN
         CALL ProgWarn( ' Error on line '//TRIM(Int2Lstr(I))//' of the distributed tower properties in file "' &
                           //TRIM(FilName)//'." Tower height fractions must be between 0.0 and 1.0.' )
         ErrStat = 1
         RETURN
      END IF


         ! Make sure the tower height increases for each entry
      IF (I > 1) THEN
         IF (TwrHtFr(I) <= TwrHtFr(I-1)) THEN
            CALL ProgWarn( ' Error on line '//TRIM(Int2Lstr(I))//' of the distributed tower properties in file "' &
                           //TRIM(FilName)//'." Tower height fraction entries must be in order of increasing height.' )
            ErrStat = 1
            RETURN
         ENDIF
      ENDIF

         ! Make sure tower width is positive
      IF ( TwrWid( I ) <= 0.0 ) THEN
         CALL ProgWarn( ' Error on line '//TRIM(Int2Lstr(I))//' of the distributed tower properties in file "' &
                           //TRIM(FilName)//'." Tower width must be positive.' )
         ErrStat = 1
         RETURN
      ENDIF

         ! Make sure the tower CD column is within range
      IF ( NTwrCDCol(I) < 1 .OR. NTwrCDCol(I) > NTwrCD ) THEN
         CALL ProgWarn( ' Error on line '//TRIM(Int2Lstr(I))//' of the distributed tower properties in file "' &
                           //TRIM(FilName)//'." Tower height CD column must be between 1 and '//TRIM(Int2Lstr(NTwrCD))//'.' )
         ErrStat = 1
         RETURN
      END IF

   END DO ! I

   !-------------------------------------------------------------------------------------------------
   ! Read section 3, Re vs CD PROPERTIES;
   ! this section has 2 header lines plus NTwrRe rows of data with NTwrCD+1 columns
   !-------------------------------------------------------------------------------------------------

      ! Read in 2 header/comment lines
   CALL ReadCom( UnIn, FilName, 'Re vs CD Properties header 1', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   CALL ReadCom( UnIn, FilName, 'Re vs CD Properties header 2', ErrStat )
   IF ( ErrStat /= 0 ) RETURN

   Fmt = '('//TRIM(Int2Lstr(NTwrCD+1))//'(2X,ES11.4e2))'

   DO I = 1, NTwrRe

      READ( UnIn,*,IOSTAT=ErrStat ) TwrRe(I), (TwrCD(I,J), J = 1, NTwrCD)

      IF ( ErrStat == 0 ) THEN
         IF ( Echo ) THEN
            WRITE (UnEc,Fmt)  TwrRe(I), (TwrCD(I,J), J = 1, NTwrCD)
         END IF
      ELSE IF ( ErrStat < 0 ) THEN
         CALL ProgWarn( ' Premature end of file while reading line '//TRIM(Int2Lstr(I))// &
                     ' of the tower Re vs CD properties in file "'//TRIM(FilName)//'."' )
         RETURN
      ELSE
         CALL ProgWarn( ' Error reading line '//TRIM(Int2Lstr(I))// &
                     ' of the tower Re vs CD properties in file "'//TRIM(FilName)//'."' )
         RETURN
      END IF

   END DO ! I

   !-------------------------------------------------------------------------------------------------
   ! close the file and return
   !-------------------------------------------------------------------------------------------------

   CLOSE( UnIn )

   RETURN

END SUBROUTINE READTwr
!====================================================================================================
! Dynamics Program aerodynamics force interface gateway
   SUBROUTINE ELEMFRC (PSI, RLOCAL, J, IBlade, VNROTOR2, VT, VNW, &
                       VNB, DFN, DFT, PMA, Initial)
 !  calculates the aerodynamic forces on one
 !  blade element.  Inputs include all velocities.
 !  Normal and tangential forces and 'A' are returned.
 ! ************************************************

USE                           Airfoil
USE                           ElOutParams
USE                           Blade
USE                           Element
USE                           ElemInflow
USE                           InducedVel
USE                           Rotor
USE                           Switch
USE                           Wind

IMPLICIT                      NONE


   ! Passed Variables:

REAL(ReKi),INTENT(OUT)     :: DFN
REAL(ReKi),INTENT(OUT)     :: DFT
REAL(ReKi),INTENT(OUT)     :: PMA
REAL(ReKi),INTENT(IN)      :: PSI
REAL(ReKi),INTENT(IN)      :: RLOCAL
REAL(ReKi),INTENT(IN)      :: VNB
REAL(ReKi),INTENT(IN)      :: VNROTOR2
REAL(ReKi),INTENT(IN)      :: VNW
REAL(ReKi),INTENT(INOUT)   :: VT

INTEGER, INTENT(IN)        :: J
INTEGER, INTENT(IN)        :: IBlade

LOGICAL,   INTENT(IN)      :: Initial

   ! Local Variables:

INTEGER                    :: ErrStat
REAL(ReKi)                 :: CDA
REAL(ReKi)                 :: CLA
REAL(ReKi)                 :: CMA
REAL(ReKi)                 :: CPHI
REAL(ReKi)                 :: PHI
REAL(ReKi)                 :: QA
REAL(ReKi)                 :: ReNum
REAL(ReKi)                 :: SPHI
REAL(ReKi)                 :: Vinduced
REAL(ReKi)                 :: VN

!DJL Start of proposed change (commented out for release)
!DJL LOGICAL(1)                 :: DYN ! Flag used in debugging dyn-equil switch
!DJL LOGICAL(1)                 :: FirstPass = .TRUE.  ! Flag to indicate the first pass
!DJL LOGICAL(1)                 :: DYNwait   = .FALSE. ! Flag to wait for GDW reactivation

!DJL REAL(ReKi)                 :: SwitchTime = 0.0 ! The last time the GDW was switched on or off
!DJL REAL(ReKi)                 :: Dyn_WT     = 0.0 ! The last time we passed GDW on threshold

!DJL CHARACTER( 15)             :: Flt2LStr
!DJL CHARACTER(100)             :: MESAGE

!DJL SAVE ! We need to save values for the next time


! Inititalize the DYN flag
!DJL IF (FirstPass) THEN
!DJL    DYN = DYNINIT ! USE DYNINIT because DYNINFLOW may be FALSE on initialization
!DJL    FirstPass = .FALSE.
!DJL ENDIF
!DJL END of proposed change


   ! initialize TanInd variables
A (J,IBLADE) = 0.0
AP(J,IBLADE) = 0.0


 !-mlb  Check for being at the center of rotation.
 ! If we are at the center of rotation, the induction equations
 !  are undefined, so let's just USE zeros.


IF ( RLOCAL < 0.01 )  THEN
   A (J,IBLADE) = 0.0
   AP(J,IBLADE) = 0.0
ELSEIF( DYNINFL .AND. R * REVS < 2.0 )  THEN   !ACH 3/10/03 This block deals with dyn. inflow problems at low tip speed
   A (J,IBLADE) = 0.0
   AP(J,IBLADE) = 0.0
   DYNINIT = .TRUE.    !Re-initialize if we begin using dynamic inflow again
ELSE

 ! Turn wake off when using dynamic inflow and tip speed goes low.  Wake will remain off.
! Eliminated with the addtion of the ELSEIF above - ACH 3/10/03
!   IF( WAKE .AND. DYNINFL .AND. (TIME > 20D0) .AND. (R * REVS < 5.0) ) THEN
!      WAKE = .FALSE.
!      WRITE(*,*) 'Wake turned off because tip speed < 5'
!   ENDIF

 ! Get induction factor = A using static airfoil coefficients
   IF ( WAKE .AND. .NOT. Initial) THEN

!DJL Start of proposed change (commented out in release version)
!DJL Testing of possible fix to GDW problem at low wind speed - 05/30/03
!DJL       IF ( DYNINFL .AND. DYN) THEN ! The GDW routines are active

!DJL          IF ( TotalInf <= 0.1 ) THEN ! Deactivate GDW
!DJL             MESAGE = " TotalInf has dropped below 0.1; GDW is being turned off. Time = " &
!DJL                      //Flt2LStr(REAL(TIME, ReKi))
!DJL             CALL ErrLog ( MESAGE, '(A)', 'ELEMFRC', 301, 'WARN' )
!DJL        DYN = .FALSE.
!DJL             SwitchTime = REAL(TIME, ReKi)
!DJL          ENDIF

!DJL       ELSEIF ( DYNINFL .AND. .NOT. DYN ) THEN ! GDW inactivated

!DJL          IF ( TotalInf >= 0.1 ) THEN ! Passed threshold

!DJL             IF ( TIME - SwitchTime > 5.0 .AND. .NOT. DYNwait) THEN ! 5 seconds elapsed since last trigger - prepare to reactivate
!DJL                DYNwait = .TRUE.
!DJL                DYN_WT  = REAL(TIME, ReKi)
!DJL             ENDIF

!DJL          ELSE

!DJL             IF (DYNwait) DYNwait = .FALSE. ! dropped below threshold - no activation of GDW

!DJL          ENDIF
!DJL       ENDIF

!DJL       IF (DYNwait .AND. TIME - DYN_WT > 0.1) THEN ! Activate GDW
!DJL          MESAGE = " TotalInf has returned above 0.1; GDW is being turned on. Time = " &
!DJL                  //Flt2LStr(REAL(TIME, ReKi))
!DJL          CALL ErrLog ( MESAGE, '(A)', 'ELEMFRC', 302, 'WARN' )
!DJL          DYNwait = .FALSE.
!DJL          DYN = .TRUE.
!DJL          SwitchTime = REAL(TIME, ReKi)
!DJL       ENDIF

      IF ( DYNINFL ) THEN
!DJL USE statement below in place of above IF - Not yet ready for prime-time (commented out for release)
!DJL      IF ( DYN ) THEN
!DJL End of proposed change

 !       USE dynamic inflow model to find A
         CALL VINDINF( J, IBlade, RLOCAL, VNW, VNB, VT, PSI ) !possibly changes VT, A, and AP
      ELSE
 !       USE momentum balance to find A
         CALL VIND( J, IBlade, RLOCAL, VNROTOR2, VNW, VNB, VT )  !changes VT, A, and AP
 !       Apply skewed-wake correction, if applicable
         IF( SKEW ) CALL VNMOD( J, IBlade, RLOCAL, PSI ) !changes A
      ENDIF
   ELSE
 !    Ignore the wake calculation entirely
      A (J,IBLADE) = 0.0
      AP(J,IBLADE) = 0.0
   ENDIF

ENDIF

Vinduced = VNW  * A(J,IBLADE)
VN = VNW + VNB - Vinduced

SumInfl = SumInfl + Vinduced * RLOCAL * DR(J)

 ! Get the angle of attack

PHI   = ATAN2( VN, VT )
ALPHA(J,IBlade) = PHI - PITNOW

CALL MPI2PI ( ALPHA(J,IBlade) )

W2(J,IBlade) = VN * VN + VT * VT

 ! Get the Reynold's number for the element
 !  Returns Reynold's number x 10^6    !bjj: Reynold's number x 10^-6 ?
ReNum = GetReynolds( SQRT(W2(J,IBlade)), C(J) )
IF (Reynolds) MulTabLoc = ReNum

 ! Get lift coefficient from dynamic stall routine if desired
 !  note that the induced velocity was calculated
 !  using the static CL, not the dynamic CL

IF ( DSTALL ) THEN
 ! USE BEDDOES dynamic stall model
   IF (Initial) THEN ! USE static data on first pass
      CALL BEDINIT (J, IBlade, ALPHA(J,IBlade))
      CALL CLCD( ALPHA(J,IBlade), CLA, CDA, CMA, NFOIL(J), ErrStat )
   ELSE
      CALL BEDDOES( W2(J,IBlade), J, IBlade, ALPHA(J,IBlade), CLA, CDA, CMA)
   ENDIF
ELSE
 ! Don't USE dynamic stall model
   CALL CLCD( ALPHA(J,IBlade), CLA, CDA, CMA, NFOIL(J), ErrStat )
ENDIF

QA       = 0.5 * RHO * W2(J,IBlade) * DR(J) * C(J)
CPHI     = COS( PHI )
SPHI     = SIN( PHI )
DFN      = ( CLA * CPHI + CDA * SPHI ) * QA
DFT      = ( CLA * SPHI - CDA * CPHI ) * QA

IF ( PMOMENT ) THEN
   PMA  = CMA * QA * C(J)
ELSE
   PMA  = 0.
   CMA  = 0.
ENDIF


 ! Save values at appropriate station

IF ( IBLADE == 1 ) THEN
   IF ( ElPrList(J) > 0 )  THEN
      AAA    ( ElPrList(J) )    = A (J,IBLADE)
      AAP    ( ElPrList(J) )    = AP(J,IBLADE)
      ALF    ( ElPrList(J) )    = ALPHA(J,IBlade) * R2D
      CDD    ( ElPrList(J) )    = CDA
      CLL    ( ElPrList(J) )    = CLA
      CMM    ( ElPrList(J) )    = CMA
      CNN    ( ElPrList(J) )    = CLA * COS(ALPHA(J,IBlade)) + CDA * SIN(ALPHA(J,IBlade))
      CTT    ( ElPrList(J) )    =-CDA * COS(ALPHA(J,IBlade)) + CLA * SIN(ALPHA(J,IBlade))
      DFNSAV ( ElPrList(J) )    = DFN
      DFTSAV ( ElPrList(J) )    = DFT
      DynPres( ElPrList(J) )    = 0.5 * RHO * W2(J,IBlade)
      PITSAV ( ElPrList(J) )    = PITNOW * R2D
      PMM    ( ElPrList(J) )    = PMA
      ReyNum ( ElPrList(J) )    = ReNum
   ENDIF

ENDIF

RETURN
END SUBROUTINE ELEMFRC

!======================================================
   SUBROUTINE VIND( J, IBlade, RLOCAL, VNROTOR2, VNW, VNB, VT )
 !  calculates the axial induction factor for each
 !  annular segment and time step.
 ! ***************************************************

USE                              Blade
USE                              Element
USE                              InducedVel


IMPLICIT                         NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)         :: RLOCAL
REAL(ReKi),INTENT(IN)         :: VNB
REAL(ReKi),INTENT(IN)         :: VNROTOR2
REAL(ReKi),INTENT(IN)         :: VNW
REAL(ReKi),INTENT(INOUT)      :: VT

INTEGER, INTENT(IN)           :: J
INTEGER, INTENT(IN)           :: IBlade

   ! Local Variables:

REAL(ReKi)                    :: A2
REAL(ReKi)                    :: A2P
REAL(ReKi)                    :: AI
REAL(ReKi)                    :: ALPHA
REAL(ReKi)                    :: ASTEP
REAL(ReKi)                    :: ATOLER2
REAL(ReKi)                    :: ATOLERBY10
REAL(ReKi)                    :: CDA
REAL(ReKi)                    :: CLA
REAL(ReKi)                    :: CMA
REAL(ReKi)                    :: DAI
REAL(ReKi)                    :: DAI1
REAL(ReKi)                    :: DELAI
REAL(ReKi)                    :: PHI
REAL(ReKi), ALLOCATABLE, SAVE :: OLD_A_NS  ( :, : )
REAL(ReKi), ALLOCATABLE, SAVE :: OLD_AP_NS ( :, : )
REAL(ReKi)                    :: SOLFACT
REAL(ReKi)                    :: VNA
REAL(ReKi)                    :: VT2_Inv
REAL(ReKi)                    :: VTA

INTEGER                       :: ICOUNT
INTEGER                       :: MAXICOUNT
INTEGER                       :: Sttus


 ! Allocate and initialize the local array on the first pass
IF ( .NOT. ALLOCATED (OLD_A_NS) ) THEN
   ALLOCATE ( OLD_A_NS ( NELM, NB) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLD_A_NS array.' )
   OLD_A_NS(:,:) = 0.0
ENDIF

 ! Allocate and initialize the local array on the first pass
IF ( .NOT. ALLOCATED (OLD_AP_NS) ) THEN
   ALLOCATE ( OLD_AP_NS ( NELM, NB) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLD_AP_NS array.' )
   OLD_AP_NS(:,:) = 0.0
ENDIF

 ! Set maximum iterations
MAXICOUNT = 1000     !bjj why isn't this a parameter?

 ! CH--  Alternate convergence criteria
ATOLER2    =  2.0 * ATOLER
ATOLERBY10 =  0.1 * ATOLER

 ! Bypass calculations for low wind speed, assume no induced velocity.

IF ( VNROTOR2 < 0.1 ) THEN
   A(J,IBLADE) = 0.0
   RETURN
ENDIF

 ! SOLFACT is solidity factor divided by 2*VNROTOR2
 ! VT2_Inv is 1./VT**2 to save computation time

IF ( RLOCAL == 0.0 ) THEN   ! Avoid div/0 in FAST2
   SOLFACT = 1.0/VNROTOR2
ELSE
   SOLFACT = NB * C(J) / ( TWOPI * RLOCAL * VNROTOR2)
ENDIF
VT2_Inv = 1. / ( VT * VT )

 !-mlb  Let's USE the old value of the A from before it was corrected for skew.
AI      = OLD_A_NS( J, IBLADE )
DAI1    = 0.05
A2P     = OLD_AP_NS( J, IBLADE )
ASTEP   = 0.5
ICOUNT  = 0

 ! Check for extremely high VN and bypass calculations if necessary

IF ( ABS( VNB ) > 100. ) THEN
   A( J, IBLADE ) = 0.0
   CALL VINDERR( VNW, VNB, 'VNB', J, IBLADE )
   RETURN
ELSEIF ( ABS( VT ) > 400. ) THEN
   A( J, IBLADE ) = 0.0
   CALL VINDERR( VNW, VT, 'VT', J, IBLADE )
   RETURN
ENDIF

A2 = AI

CALL AXIND ( VNW, VNB, VNA, VTA, VT, VT2_Inv, VNROTOR2, A2, A2P, &
            J, SOLFACT, ALPHA, PHI, CLA, CDA, CMA, RLOCAL )

DAI = A2  - AI

DELAI = ASTEP * DAI

 !CH--  Modification of mlb's proposed change
 ! Must pass two criteria. If we have crossed zero many times
 !  then the first criterion will be easier to meet than the second
 !  because ASTEP will be small (but the second is relaxed to ATOLER2)

DO WHILE ( ABS( DELAI ) > ATOLERBY10 .AND. ABS(DAI) > ATOLER2 )

   ICOUNT = ICOUNT + 1

   A2 = AI

   CALL AXIND ( VNW, VNB, VNA, VTA, VT, VT2_Inv, VNROTOR2, A2, A2P, &
             J, SOLFACT, ALPHA, PHI, CLA, CDA, CMA, RLOCAL )

   DAI = A2  - AI

   DELAI = ASTEP * DAI

 ! Test for convergence, program warning after 1000 iterations

   IF ( ICOUNT > MAXICOUNT ) THEN
      CALL ProgWarn( 'Induction factor calculation did not converge after'//TRIM(Int2LStr(MAXICOUNT))// &
                     ' iterations. AeroDyn will continue using induction factors from previous successful time step.' )
      A2  = OLD_A_NS (J,IBLADE)
      A2P = OLD_AP_NS(J,IBLADE)
      EXIT
   ENDIF

 ! Reduce step size after a zero crossing
 !CH--  Put floor under ASTEP to keep it reasonable after many zero crossings

   IF( NINT( SIGN(1., DAI) ) /= NINT( SIGN(1., DAI1) ) ) ASTEP = MAX( 1.0E-4, 0.5*ASTEP )

   AI   = AI + DELAI
   DAI1 = DELAI

END DO

 ! Passed test, we're done
A (J,IBLADE) = A2
AP(J,IBLADE) = A2P
VT = VT * ( 1. + A2P )  !bjj: why are we changing the total velocity?
OLD_A_NS  (J,IBLADE) = A2
OLD_AP_NS (J,IBLADE) = A2P



RETURN
END SUBROUTINE VIND



 ! ***************************************************
   SUBROUTINE VINDERR( VNW, VX, VID, J, IBLADE )
 !  used to write warning messages to the screen
 !  when VN or VT is high.
 ! ***************************************************



IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: VNW
REAL(ReKi),INTENT(IN)      :: VX

INTEGER   ,INTENT(IN)      :: IBLADE
INTEGER   ,INTENT(IN)      :: J

CHARACTER(  *),INTENT(IN)  :: VID

   ! Local Variables:

INTEGER   , SAVE           :: NERRORS = 0

LOGICAL,    SAVE           :: AFLAG   = .FALSE.


 ! Don't write messages if we've already done it 5 times

IF ( AFLAG ) RETURN

NERRORS = NERRORS + 1

   CALL ProgWarn( ' High '//TRIM(VID)//' velocity encountered during induction factor calculation.' )
   CALL WrScr( '  Blade number '//TRIM(Int2LStr(IBLADE))//', Element number '//TRIM(Int2LStr(J )) )
   CALL WrScr( '  VNW = '       //TRIM(Flt2LStr(VNW))//', '//TRIM(VID)//' = '//TRIM(Flt2LStr(VX)) )

IF ( NERRORS >= 5 ) THEN
   AFLAG = .TRUE.
   CALL ProgWarn( ' Induced velocity warning written 5 times. '//&
                 ' The message will not be repeated, though the condition may persist.' )
ENDIF



RETURN
END SUBROUTINE VINDERR



 ! ******************************************************
   SUBROUTINE AXIND ( VNW, VNB, VNA, VTA, VT, VT2_Inv, VNROTOR2, A2, &
                      A2P, J, SOLFACT, ALPHA, PHI, CLA, CDA, CMA, RLOCAL )
 !  calculates a new axial induction factor from
 !   given values of velocities and geometry.  This routine
 !   is called by vind as part of the iteration process
 ! ******************************************************

USE                           Element
USE                           InducedVel
USE                           Airfoil, ONLY: NFoil
USE                           Switch


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(INOUT)   :: A2
REAL(ReKi),INTENT(INOUT)   :: A2P
REAL(ReKi),INTENT(OUT)     :: ALPHA
REAL(ReKi),INTENT(OUT)     :: CDA
REAL(ReKi),INTENT(OUT)     :: CLA
REAL(ReKi),INTENT(OUT)     :: CMA
REAL(ReKi),INTENT(OUT)     :: PHI
REAL(ReKi),INTENT(IN)      :: RLOCAL
REAL(ReKi),INTENT(IN)      :: SOLFACT
REAL(ReKi),INTENT(OUT)     :: VNA
REAL(ReKi),INTENT(IN)      :: VNB
REAL(ReKi),INTENT(IN)      :: VNROTOR2
REAL(ReKi),INTENT(IN)      :: VNW
REAL(ReKi),INTENT(IN)      :: VT
REAL(ReKi),INTENT(IN)      :: VT2_Inv
REAL(ReKi),INTENT(OUT)     :: VTA

INTEGER   ,INTENT(IN)      :: J

   ! Local Variables:

REAL(ReKi)                 :: CH
REAL(ReKi)                 :: CPhi                                                     ! COS( PHI )
REAL(ReKi), SAVE           :: HUBLOSS = 1
REAL(ReKi), SAVE           :: LOSS = 1
REAL(ReKi)                 :: SPHI
REAL(ReKi)                 :: SWRLARG
REAL(ReKi), SAVE           :: TIPLOSS = 1
REAL(ReKi)                 :: W2

INTEGER                    :: ErrStat


VNA    = VNW * ( 1. - A2 ) + VNB
VTA    = VT  * ( 1. + A2P )

 ! Get airfoil CL and CD
PHI    = ATAN2( VNA, VTA )
ALPHA  = PHI - PITNOW

CALL MPI2PI ( ALPHA )

CALL CLCD ( ALPHA, CLA, CDA, CMA, NFoil(J), ErrStat )

W2   = VNA * VNA + VTA * VTA
SPHI = VNA/SQRT( W2 )
CPhi = COS( Phi )

 ! Calculate new value of A.  Optionally include normal force due to drag.

CH = W2*SOLFACT*( CLA*CPhi + EqAIDmult*CDA*SPhi )


 ! Get the tip loss values for the element (if they change)
IF (TLOSS) CALL GetTipLoss (J, SPHI, TIPLOSS, RLOCAL)

 ! Get the hub loss values for the element (if they change)
IF (HLOSS) CALL GetPrandtlLoss (HLCNST(J), SPHI, HUBLOSS)

 ! Get the total loss for the element
LOSS = TIPLOSS * HUBLOSS

 ! Check for diverging CH and correct if necessary

IF ( ABS( CH ) > 2. ) CH = SIGN( 2., CH )

IF ( CH < 0.96*LOSS ) THEN
   A2 = 0.5*( 1 - SQRT( 1.0 - CH/LOSS ) )
ELSE
   A2 = 0.1432 + SQRT( -0.55106 + .6427*CH/LOSS)
ENDIF

 ! Calculate induced swirl (a') if desired.
 !  From C. Ross Harmon's paper on PROPX.

IF ( SWIRL ) THEN
   IF ( EquilDT )  THEN     ! USE PROP-PC style tangential induction equation with the addition of the drag term.
         ! Because of the singularity that occurs when phi approaches zero,
         ! let's test for small phi and set a' equal to a small, negative number.
      IF ( ( ABS( SPhi ) > 0.01 ) .AND. ( ABS( CPhi ) > 0.01 ) )  THEN
         A2P = SOLFACT*( CLA*SPhi - CDA*CPhi )*( 1.0 + A2P )*VNROTOR2/( 4.0*LOSS*SPhi*CPhi )
      ELSEIF ( ABS( SPhi ) > 0.01 )  THEN   ! Tangential velocity near zero, phi near 90 degrees.
         A2P = SOLFACT*( CLA*SPhi - CDA*SIGN( 0.01, CPhi ) )*( 1.0 + A2P )*VNROTOR2/( 4.0*LOSS*SPhi*SIGN( 0.01, CPhi ) )
      ELSE   ! Normal velocity near zero, phi near 0 degrees.
         A2P = SOLFACT*( CLA*SIGN( 0.01, SPhi ) - CDA*CPhi )*( 1.0 + A2P )*VNROTOR2/( 4.0*LOSS*SIGN( 0.01, SPhi )*CPhi )
      ENDIF
   ELSE
      SWRLARG = 1.0 + 4.0*LOSS*A2*VNW*VNA*VT2_Inv
      IF ( SWRLARG < 0.0 ) THEN
         A2P = 0.0
      ELSE
         A2P = 0.5*( -1.0 + SQRT( SWRLARG ) )
      ENDIF
   ENDIF
ELSE
   A2P = 0.0
ENDIF



RETURN
END SUBROUTINE AXIND


 ! ***************************************************
   FUNCTION GetReynolds( WindSpd, ChordLen )
 !  computes the Reynolds number for the element, divided by 1.0E6
 ! ***************************************************

   USE                           Wind,       ONLY: KinVisc

IMPLICIT                      NONE


   ! Passed Variables:

REAL(ReKi),INTENT(IN)      :: WindSpd
REAL(ReKi),INTENT(IN)      :: ChordLen

   ! function definition
REAL(ReKi)     :: GetReynolds

GetReynolds = 1.0E-6 * WindSpd * ChordLen / KinVisc


RETURN
END FUNCTION GetReynolds

 ! ***************************************************
   SUBROUTINE GetTipLoss( J, SPHI, TIPLOSS, RLOCAL )
 !  computes the tip loss constant for element J
 !  TIPLOSS is returned to AXIND
 ! Uses the Prandtl tip loss model with a correction
 !  from Georgia Tech (2002 ASME Wind Energy Symposium)
 ! ***************************************************


USE                           Blade
USE                           Element
USE                           Switch


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi), INTENT(IN)     :: SPHI
REAL(ReKi), INTENT(IN)     :: RLOCAL
REAL(ReKi), INTENT(OUT)    :: TIPLOSS

INTEGER   , INTENT(IN)     :: J


   ! Local Variables:

REAL(ReKi)                 :: Dist2pt7 = 0.7 ! current element distance to r/R = 0.7
REAL(ReKi)                 :: OLDDist7       ! previous element distance to r/R = 0.7
REAL(ReKi)                 :: percentR
REAL(ReKi), SAVE           :: TLpt7 ! Tiploss factor at r/R = 0.7

INTEGER                    :: Jpt7 = 0 ! The element closest to r/R = 0.7

LOGICAL,    SAVE           :: FirstPass = .TRUE.



 ! Calculate PRANDTL tip loss model
CALL GetPrandtlLoss( TLCNST(J), SPHI, TIPLOSS )


 ! Apply Georgia Tech correction to Prandtl model if activated
IF (GTECH) THEN
   percentR = RLOCAL/R

 ! Search for the element closest to r/R = 0.7
   IF (FirstPass) THEN
    ! If the current element is closer than the previous, update values
      IF ( ABS(percentR - 0.7) < Dist2pt7 ) THEN
         OLDDist7 = Dist2pt7
         Dist2pt7 = ABS(percentR - 0.7)
         Jpt7 = J
         TLpt7 = TIPLOSS
      ENDIF
      IF (J == NELM) THEN ! We're done after one pass through the blades
         FirstPass = .FALSE.
      ELSE
         RETURN ! Don't do the correction until we calculate the correct TLpt7
      ENDIF
   ENDIF

   IF ( J == Jpt7 ) TLpt7 = TIPLOSS ! Update the value of TLpt7 at the proper element

 ! Do the actual Georgia Tech correction to the Prandtl model
   IF (percentR >= 0.7) THEN
      TIPLOSS = (TIPLOSS**0.85 + 0.5 ) / 2.0
   ELSE
      TIPLOSS = 1.0 - percentR*(1.0 - TLpt7)/0.7
   ENDIF
ENDIF



RETURN
END SUBROUTINE GetTipLoss


 ! ***************************************************

   SUBROUTINE GetPrandtlLoss( LCnst, SPHI, PrLOSS )
 !  computes the hub loss constant for element J
 !  HUBLOSS is returned to AXIND
 ! Uses the Prandtl loss model
 ! ***************************************************

IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: LCnst
REAL(ReKi),INTENT(OUT)     :: PrLOSS
REAL(ReKi),INTENT(IN)      :: SPHI


   ! Local Variables:

REAL(ReKi)                 :: F


 ! Calculate PRANDTL loss model
 !  Check values of SPHI to save runtime.
IF ( ABS( SPHI ) < 1.E-4 ) THEN
   PrLOSS = 1.0
ELSE
 ! USE ABS function to account for unusual PHI.
   F = ABS( LCnst / SPHI )
 !  Check values of F to avoid underflow of EXP function.
   IF ( F < 7. ) THEN
      PrLOSS = ACOS( EXP( -F ) ) / PIBY2
   ELSE
      PrLOSS = 1.0
   ENDIF
ENDIF


RETURN
END SUBROUTINE GetPrandtlLoss

!====================================================================================================
SUBROUTINE GetTwrInfluence (VX, VY, InputPosition)
!  Computes tower shadow or dam influence on the blade
!  Note that this routine assumes there are NO tower deflections.
!
!  Use the Riso tower dam model of Bak, Madsen and Johansen
!  This model is based on potential flow and is applicable in front of and behind the tower,
!  although in the wake we use the method of Powles
!  PJM, NREL
!
! bjj, jmj: this function should return the influence parameters, which will be based on some mean
!  wind direction (or possibly the direction at the tower) for a given height, instead of using the
!  local velocity/direction so that all points on a horizontal slice of air can use the same
!  deficit at each given time.  Will need the tower position, too.
!====================================================================================================

   USE                           TwrProps
   USE                           Rotor,      ONLY: HH


   IMPLICIT                      NONE


      ! Passed Variables:

   REAL(ReKi), INTENT(INOUT)  :: VX                   ! on input, U-velocity without tower effect; on output, U-velocity including tower effect
   REAL(ReKi), INTENT(INOUT)  :: VY                   ! on input, V-velocity without tower effect; on output, V-velocity including tower effect
   REAL(ReKi), INTENT(IN)     :: InputPosition(3)     !velocities in global coordinates with tower effect

      ! Local Variables:

   LOGICAL, SAVE              :: AFLAG   = .FALSE.    ! set to .TRUE. on possible tower strike

   REAL(ReKi)                 :: ANGLE                ! Angle determining whether blade is upwind or downwind of tower
   REAL(ReKi)                 :: CenterDist           ! Distance from blade element to wake centerline
   REAL(ReKi)                 :: phi                  ! Angle between x-axis and horizontal wind direction based upon instantaneous velocities VY and VX (at the blade element)
   REAL(ReKi)                 :: CosPhi               ! COS(phi)
   REAL(ReKi)                 :: SinPhi               ! SIN(phi)
   REAL(ReKi)                 :: Distance             ! Normalized horizontal distance from tower to blade element
   REAL(ReKi)                 :: SHADOW               ! Value of the tower shadow deficit at the blade element
   REAL(ReKi)                 :: THETA                ! Angle between x-axis and line from tower to blade element
   REAL(ReKi)                 :: TwrCD_Station        ! Drag coefficient of the tower
   REAL(ReKi)                 :: TwrRad               ! Radius of the tower at the height of interest
   REAL(ReKi)                 :: WIDTH                ! Half width of the wake after accounting for wake expansion proportional to square root of Distance

   REAL(ReKi)                 :: V_total              ! total freestream wind speed
   REAL(ReKi)                 :: VX_wind              ! tower influenced wind speeds in wind coordinates
   REAL(ReKi)                 :: VY_wind              ! tower influenced wind speeds in wind coordinates
   REAL(ReKi)                 :: WindXInf             ! Influence of the tower on X wind velocity in wind reference frame
   REAL(ReKi)                 :: WindYInf             ! Influence of the tower on Y wind velocity in wind reference frame
   REAL(ReKi)                 :: Xtemp                ! Temporary variable used in tower dam calculations
   REAL(ReKi)                 :: Xtemp2               ! Temporary variable used in tower dam calculations
   REAL(ReKi)                 :: Xwind                ! X Location of element in a wind-based coordinate system
   REAL(ReKi)                 :: Yg                   ! Variable used to smooth dam effect above the tower
   REAL(ReKi)                 :: Ytemp2               ! Temporary variable used in tower dam calculations
   REAL(ReKi)                 :: Ywind                ! Y Location of element in a wind-based coordinate system

   REAL(ReKi)                 :: ZGrnd                ! distance between position and undeflected hub


   !-------------------------------------------------------------------------------------------------
   ! This subroutine is only valid for TwrPotent and TwrShadow features
   !-------------------------------------------------------------------------------------------------
   IF (.NOT. TwrPotent .AND. .NOT. TwrShadow) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Initialize some variables
   !-------------------------------------------------------------------------------------------------
   ZGrnd   = InputPosition(3) - HH                    ! distance between position and hub      !BJJ: this should really be the tower height (position), not HH
   V_total = SQRT( VX**2 + VY**2 )                    ! total wind speed

   !-------------------------------------------------------------------------------------------------
   ! Tower influence calculations aren't necessary for zero velocity
   !-------------------------------------------------------------------------------------------------
   IF ( V_total <= 0.0 ) RETURN

   !-------------------------------------------------------------------------------------------------
   ! For the current element location, get the appropriate tower properties and calculate the
   ! element distance from from the tower.
   ! BJJ: If we're above the tower top, is the radius zero?
   !-------------------------------------------------------------------------------------------------

   CALL GetTwrSectProp (InputPosition(:), V_total, TwrRad, TwrCD_Station)     ! Get the tower properties for the current element location

   Distance = SQRT ( InputPosition(1)**2 + InputPosition(2)**2 ) / TwrRad     ! normalized distance to tower center

      ! Check for tower strike
   IF ( Distance < 1.0 ) THEN    ! potentially inside the tower            !bjj: only if we're not ABOVE the tower, though....

      IF (ZGrnd < 0.0) THEN      !bjj added this condition.... check that it's correct

         IF( .NOT. AFLAG) THEN
            CALL ProgWarn( ' Tower model temporarily disabled due to possible tower strike.'// &
                           ' This message will not be repeated though the condition may persist.' )
               !write a blank line (so FAST doesn't write over it)
            CALL WrScr( ' ' ) 
            AFLAG = .TRUE.
         ENDIF

         RETURN

      END IF

   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Store the wind direction for later
   !-------------------------------------------------------------------------------------------------

   phi    = ATAN2( VY, VX )                                 ! angle between x-axis and instantaneous horizontal wind direction
   CosPhi = COS( phi )
   SinPhi = SIN( phi )


   !-------------------------------------------------------------------------------------------------
   !  Calculate the influence due to potential flow around tower based on velocity at element
   !-------------------------------------------------------------------------------------------------

   IF ( TwrPotent ) THEN

         ! When above the tower, smooth the transition to the free-stream
      IF ( ZGrnd > 0 )  THEN
         Yg = SQRT( InputPosition(2)**2 + ZGrnd**2 )
      ELSE
         Yg = InputPosition(2)
      ENDIF

         ! Get the element location in the wind reference frame
      Xwind = InputPosition(1) * CosPhi + Yg                * SinPhi
      Ywind = Yg               * CosPhi - InputPosition(1)  * SinPhi

         ! Normalize the location coordinates
      Xwind = Xwind / TwrRad
      Ywind = Ywind / TwrRad


      Xtemp  = Xwind + Tower_Wake_Constant !PJM fixed this error 3/30/06
      Xtemp2 = Xtemp * Xtemp
      Ytemp2 = Ywind * Ywind


         ! Calculate the tower influence
      WindXInf =  1.0 - (Xtemp2 - Ytemp2)/( (Xtemp2 + Ytemp2)**2 ) &
                      + TwrCD_Station/TwoPi * Xtemp/(Xtemp2 + Ytemp2)
      WindYInf = -2.0 * (Xtemp * Ywind)/( (Xtemp2 + Ytemp2)**2 ) &   !PJM fixed sign error 3/30/06 added minus sign from Bak
                      + TwrCD_Station/TwoPi * Ywind/(Xtemp2 + Ytemp2)
   ELSE
      WindXInf = 1.0
      WindYInf = 0.0
   ENDIF

   !-------------------------------------------------------------------------------------------------
   ! Calculate the influence of tower shadow if user specifies and we are downwind of the tower
   !-------------------------------------------------------------------------------------------------

   IF ( TwrShadow ) THEN

      theta  = ATAN2( InputPosition(2), InputPosition(1) )     ! angle between x-axis and line from tower to blade element
      angle  = ABS( theta - phi )

      CALL mPi2Pi( angle )                                     ! angle difference between -pi and pi
      angle  = ABS(angle)        !BJJ: SHOULDN'T THIS BE ABS(angle)?  I'm adding it here....

      IF ( angle <= PiBy2 ) THEN ! We are downwind of the tower in shadow territory

         width = SQRT ( Distance )

            ! Calculate how far we are from the free-stream centerline
         CenterDist = Distance * SIN ( angle )                 ! bjj: non-negative because angle is non-negative (now)

            ! If above hub height apply shadow in arc above hub to maintain a continuous deficit function.  Somewhat of a nacelle deficit.
         IF ( ZGrnd > 0.0 )  THEN
            CenterDist = SQRT( CenterDist**2 + ZGrnd**2 )
         END IF

            ! See if we are in the wake. If not, then no velocity deficit

         IF ( CenterDist < width ) THEN                        ! We are in the wake

            shadow = ( COS( PiBy2 * CenterDist / width ) )**2 * TwrCD_Station / width
            WindXInf = 1.0 - shadow                               ! Overwrites the potential flow solution in the x direction only (BJJ: longitudinal? not x?)

            WindXInf = MAX( WindXInf, REAL(0.0, ReKi) )           ! Assume tower does not reverse flow direction

         END IF

      END IF ! angle <= PiBy2

   END IF !TwrShadow


   !-------------------------------------------------------------------------------------------------
   ! Apply the tower influence to the input wind speeds
   !-------------------------------------------------------------------------------------------------

   VX_wind = WindXInf*V_total
   VY_wind = WindYInf*V_total

      ! Need to transpose these back to the global reference frame
   VX = VX_wind * CosPhi - VY_wind * SinPhi
   VY = VY_wind * CosPhi + VX_wind * SinPhi

   RETURN

END SUBROUTINE GetTwrInfluence

!====================================================================================================
SUBROUTINE GetTwrSectProp (InputPosition, VelHor, TwrElRad, TwrElCD)
!  Returns the tower radius and CD for the vertical location
!   of the element currently being evaluated for tower influence.
!====================================================================================================

   USE                           Rotor,      ONLY: HH
   USE                           TwrProps

   IMPLICIT                      NONE


      ! Passed Variables:
   REAL(ReKi), INTENT(IN)     :: InputPosition(3)  ! Location where tower properties are desired
   REAL(ReKi), INTENT(IN)     :: VelHor            ! The horizontal wind speed, used to get Reynolds number, if necessary
   REAL(ReKi), INTENT(OUT)    :: TwrElRad          ! Radius of the tower element
   REAL(ReKi), INTENT(OUT)    :: TwrElCD           ! Drag coefficient of the tower element


      ! Local Variables:
   REAL(ReKi)                 :: P1        ! Interpolation weighting factor
   REAL(ReKi)                 :: P2        ! Interpolation weighting factor
   REAL(ReKi)                 :: TwrElCD1  ! Dummy variable for 2-D interpolation
   REAL(ReKi)                 :: TwrElCD2  ! Dummy variable for 2-D interpolation
   REAL(ReKi)                 :: TwrElHt   ! Non-dimensional height of the tower element
   REAL(ReKi)                 :: TwrElRe   ! Reynold's # of the tower element

   INTEGER                    :: N1        ! Index position in table for interpolation
   INTEGER                    :: N2        ! Index position in table for interpolation
   INTEGER                    :: N1P1      ! Index position + 1 in table for interpolation
   INTEGER                    :: N2P1      ! Index position + 1 in table for interpolation


   !-------------------------------------------------------------------------------------------------
   ! Get the tower radius, TwrElRad, by interpolating into the TwrWid(:) array
   !-------------------------------------------------------------------------------------------------

   TwrElHt  = InputPosition(3) / HH                                  !!!!BJJ!!!! HH????
   TwrElRad = 0.5*InterpBin( TwrElHt, TwrHtFr, TwrWid, N2, NTwrHt )

   !-------------------------------------------------------------------------------------------------
   ! Get the section CD, TwrElCD, by interpolating into the TwrCD(:,:) array
   !-------------------------------------------------------------------------------------------------

   IF ( NTwrRe == 1 ) THEN                                  ! There is only one Re row

      IF ( NTwrCD == 1 ) THEN                               ! There is only one CD column
         TwrElCD = TwrCD(1,1)
      ELSE IF ( NTwrHt == 1 ) THEN                          ! There is more than one column of CD, but only one used
         TwrElCD = TwrCD(1,NTwrCDCol(1))
      ELSE                                                  ! Interpolate;  this will be the same Indx as before...
         TwrElCD = InterpStp( TwrElHt, TwrHtFr, TwrCD(1,:), N2, NTwrHt )
      END IF

   ELSE                                                     ! There are multiple Re rows

      TwrElRe = GetReynolds( VelHor, 2.0*TwrElRad )

      IF ( NTwrCD == 1 ) THEN                               ! There is only one CD column
         TwrElCD = InterpBin( TwrElRe, TwrRe, TwrCD(:,1), N1, NTwrRe )
      ELSE IF ( NTwrHt == 1 ) THEN                          ! Interpolate over Re only
         TwrElCD = InterpBin( TwrElRe, TwrRe, TwrCD(:,NTwrCDCol(1)), N1, NTwrRe )
      ELSE                                                  ! A 2-D interpolation is needed
         CALL LocateBin( TwrElRe, TwrRe, N1, NTwrRe )

            ! Let's use nearest-neighbor extrapolation with bi-linear interpolation:

         N1   = MIN( MAX( N1, 1 ), NTwrRe-1 )
         N1P1 = N1+1

         P1   = MIN( MAX( (TwrElRe - TwrRe(N1))   / (TwrRe(N1P1)   - TwrRe(N1))  , REAL(0.0, ReKi) ), REAL(1.0, ReKi) )

         N2P1 = N2 + 1
         P2   = MIN( MAX( (TwrElHt - TwrHtFr(N2)) / (TwrHtFr(N2P1) - TwrHtFr(N2)), REAL(0.0, ReKi) ), REAL(1.0, ReKi) )


         TwrElCD1 = TwrCD(N1,N2  ) + P1 * ( TwrCD(N1P1,N2  ) - TwrCD(N1,N2  ) )
         TwrElCD2 = TwrCD(N1,N2P1) + P1 * ( TwrCD(N1P1,N2P1) - TwrCD(N1,N2P1) )


         TwrElCD = TwrElCD1 + P2 * ( TwrElCD2 - TwrElCD1 )
      END IF

   END IF


RETURN
END SUBROUTINE GetTwrSectProp
!====================================================================================================
FUNCTION AD_WindVelocityWithDisturbance( InputPosition, ErrStat  )
!  This function computes the (dimensional) wind velocity components at the location InputPosition
!  in the inertial frame of reference, including any tower shadow defecit.
!  ** Formerly SUBROUTINE VEL and later SUBROUTINE VWrel2G( VNRotor2, At_Hub ) **
!----------------------------------------------------------------------------------------------------

   USE                              AeroTime,   ONLY: Time  ! input time, used to get the wind speed
   USE                              Rotor,      ONLY: HH
   USE                              TwrProps,   ONLY: PJM_Version, TShadC1, TShadC2



   IMPLICIT                         NONE


      ! Passed Variables:

   REAL(ReKi),INTENT(IN)            :: InputPosition(3)
   INTEGER,   INTENT(OUT),OPTIONAL  :: ErrStat


      ! function definition

   REAL(ReKi)                      :: AD_WindVelocityWithDisturbance(3)


      ! Local variables

   TYPE(InflIntrpOut)               :: InflowVel

   REAL(ReKi)                       :: angle    ! absolute difference between theta and phi
   REAL(ReKi)                       :: dist     ! distance from blade element to wake centerline
   REAL(ReKi)                       :: phi      ! angle between x-axis and instantaneous horizontal wind direction
   REAL(ReKi)                       :: RADIUS   ! horizontal distance from tower to blade element ** BJJ NOTE: in actuality, it's the distance from the undeflected tower centerline, not the actual tower
   REAL(ReKi)                       :: ROOTR    ! SQRT(radius)
   REAL(ReKi)                       :: SHADOW
   REAL(ReKi)                       :: TEMP
   REAL(ReKi)                       :: theta    ! Angle between x-axis and line from tower to blade element
   REAL(ReKi)                       :: width    ! half width of the wake after accounting for wake expansion proportional to square root of RADIUS

   INTEGER                          :: Sttus


      ! Get the undisturbed velocity

   InflowVel = WindInf_GetVelocity( REAL(Time, ReKi), InputPosition, Sttus)

   IF ( PRESENT(ErrStat) ) THEN  !BJJ: This is a hack job used for A2AD
      ErrStat = Sttus
      IF (Sttus /= 0) CALL ProgWarn( ' Error getting velocity in AeroDyn/AD_WindVelocityWithDisturbance().' )
   ELSE
      IF (Sttus /=0) CALL ProgAbort( ' Error getting velocity in AeroDyn/AD_WindVelocityWithDisturbance().' )
   END IF

   AD_WindVelocityWithDisturbance(:) = InflowVel%Velocity(:)


         ! Add the tower influence to the undisturbed velocity.

   IF ( PJM_Version ) THEN

      CALL GetTwrInfluence ( AD_WindVelocityWithDisturbance(1), AD_WindVelocityWithDisturbance(2), InputPosition(:) )

   ELSE !Old tower version

         ! Apply tower shadow if the blade element is in the wake

      IF ( TShadC2 > 0.0 ) THEN ! Perform calculations only if the wake strength is positive

            ! Bypass tower shadow for zero horizontal wind, check U-component first to save time.

         IF ( AD_WindVelocityWithDisturbance(1) /= 0.0 .OR. AD_WindVelocityWithDisturbance(2) /= 0.0 ) THEN

            phi    = ATAN2( AD_WindVelocityWithDisturbance(2), AD_WindVelocityWithDisturbance(1) )    ! angle between x-axis and instantaneous horizontal wind direction
            theta  = ATAN2( InputPosition(2),                  InputPosition(1)                  )    ! angle between x-axis and line from tower to blade element
            angle  = ABS( theta - phi )

            CALL MPi2Pi( angle )
            angle = ABS( angle )

            IF ( angle <= PiBy2 ) THEN  ! Skip cases where we are upwind of the tower -- bjj: DOES THIS ACTUALLY WORK? WHAT ABOUT

               radius = SQRT( InputPosition(1)**2 + InputPosition(2)**2 )                             ! bjj: shouldn't this be relative to the hub position?

               RootR  = SQRT( radius )
               width  = TShadC1 * RootR                                                               ! half width of the wake after accounting for wake expansion proportional to square root of RADIUS

               IF ( width > 0 ) THEN   ! Skip cases with zero width or radius so we don't divide by zero

                  dist = radius * SIN( angle )                                                        ! distance from blade element to wake centerline
                  IF ( InputPosition(3) > HH )  THEN                                                  ! Apply shadow in arc above hub to maintain a continuous deficit function.  Somewhat of a nacelle deficit.
                     dist = SQRT( dist**2 + (InputPosition(3)-HH)**2 )                                !bjj: I think this should use hub position, not HH
                  END IF

                  IF ( width > dist ) THEN   ! There is velocity deficit in the wake only.
                     temp   = COS ( PiBy2 * dist/width )
                     shadow = TShadC2/RootR * temp * temp

                        ! Adjust only the horizontal components; vertical wind is not changed

                     AD_WindVelocityWithDisturbance(1:2) = AD_WindVelocityWithDisturbance(1:2) * ( 1. - shadow )

                  END IF ! width > dist
               END IF ! width > 0

            END IF ! angle <= PiBy2

         END IF !AD_WindVelocityWithDisturbance(1) /= 0.0 .OR. AD_WindVelocityWithDisturbance(2) /= 0.0
      END IF ! TShadC2 > 0.0

   END IF

RETURN

END FUNCTION AD_WindVelocityWithDisturbance
!====================================================================================================
   SUBROUTINE DiskVel
 !  calculates the mean velocities relative to the rotor disk
 !  calls routine to get wind velocity at a specified location
 !
 !  Updated on 08/12/97 xyz-direction changed
 !  Combined VELD and GETSKEW 04/24/01
 !  Updated 12/1/09 to use new inflow module; WindInf_ADhack_diskVel MUST be replaced!
 ! ********************************************

USE                           AeroTime, ONLY: time
USE                           InflowWind
USE                           Rotor
USE                           Switch
USE                           Wind


IMPLICIT                      NONE


REAL(ReKi)                 :: Vinplane
REAL(ReKi)                 :: VXY
REAL(ReKi)                 :: AvgInfVel(3)
REAL(ReKi)                 :: Position(3)
INTEGER                    :: ErrStat

Position = (/0.0, 0.0, HH /)
AvgInfVel(:) = WindInf_ADhack_diskVel( REAL(Time, ReKi), Position, ErrStat )

VXY   = AvgInfVel(1) * CYaw - AvgInfVel(2) * SYaw

 ! Mean velocities in rotor disk coord. Includes yaw rate.
 ! X = Normal to plane, parallel to axis of rotation DOWNWIND
 ! Y = Inplane, horizontal to left (looking downwind)
 ! Z = Inplane, vertical up

VROTORX = VXY * CTILT + AvgInfVel(3) * STILT

VROTORY = AvgInfVel(1) * SYaw + AvgInfVel(2) * CYaw + YAWVEL

VROTORZ = -1.* VXY * STILT + AvgInfVel(3) * CTILT

 ! Skewed wake correction not needed for GDW
IF (.NOT. DYNINFL) THEN
 !  Set SKEW flag and assign values to related parameters
 !   used in the skewed wake correction.

 ! Vinplane is the resultant in-plane velocity
   Vinplane = SQRT( VROTORY * VROTORY + VROTORZ * VROTORZ )

 ! SKEW is TRUE if there is a cross flow, FALSE otherwise.
   IF ( Vinplane >= 1.0E-3 ) THEN
      SKEW   = .TRUE.
      SDEL   =  VROTORY/Vinplane
      CDEL   = -VROTORZ/Vinplane
      ANGFLW = ATAN2( ABS( VROTORX - AVGINFL ), Vinplane )
   ELSE
      SKEW   = .FALSE.
   ENDIF

ENDIF



RETURN
END SUBROUTINE DiskVel


 ! ****************************************
   SUBROUTINE VNMOD( J, IBlade, RLOCAL, PSI )
 !  applies the skewed wake correction
 !   to the axial induction factor A.
 ! ****************************************

USE                           Blade
USE                           Element
USE                           Wind


IMPLICIT                      NONE


   ! Passed Variables:

REAL(ReKi),INTENT(IN)      :: PSI
REAL(ReKi),INTENT(IN)      :: RLOCAL

INTEGER,   INTENT(IN)      :: J
INTEGER,   INTENT(IN)      :: IBlade

   ! Local Variables:

REAL(ReKi)                 :: BB
REAL(ReKi)                 :: SANG


SANG = SIN( ANGFLW )
BB   = 0.7363 * SQRT( ( 1. - SANG )/(1. + SANG) )

A(J,IBLADE) = A(J,IBLADE) * ( 1. + 2. * RLOCAL/R * BB *  &
             ( SDEL * SIN( PSI ) + CDEL * COS( PSI ) )  )



RETURN
END SUBROUTINE VNMOD


 ! **********************************************************
   SUBROUTINE BEDINIT( J, IBlade, ALPHA )
 !  calculates initial values of Beddoes 'f' arrays
 ! **********************************************************

USE                           Airfoil
USE                           Bedoes


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(INOUT)   :: ALPHA

INTEGER,   INTENT(IN)      :: J
INTEGER,   INTENT(IN)      :: IBlade


   ! Local Variables:

REAL(ReKi)                 :: AOL1
REAL(ReKi)                 :: CNA1
REAL(ReKi)                 :: FSPA
REAL(ReKi)                 :: FSPB
REAL(ReKi)                 :: FSPCA
REAL(ReKi)                 :: FSPCB
REAL(ReKi)                 :: P
REAL(ReKi)                 :: P1
REAL(ReKi)                 :: P2
REAL(ReKi)                 :: SRFP
REAL(ReKi)                 :: TEMP

INTEGER                    :: I
INTEGER                    :: I1
INTEGER                    :: I1P1
INTEGER                    :: I2
INTEGER                    :: I2P1
INTEGER                    :: N
INTEGER                    :: NP1


ANE(J,IBLADE) = ALPHA
AFE(J,IBLADE) = ALPHA

I = NFOIL(J)

IF ( NTables(I) > 1 ) THEN
   MulTabLoc = MIN( MAX( MulTabLoc, MulTabMet(I,1) ), MulTabMet(I,NTables(I)) )
   CALL LocateBin( MulTabLoc, MulTabMet(I,1:NTables(I)), N, NTables(I) )

   IF (N == 0 ) THEN
      CNA1 = CNA(I,1)
      AOL1 = AOL(I,1)
   ELSE IF( N == NTables(I) ) THEN
      CNA1 = CNA(I,N)
      AOL1 = AOL(I,N)
   ELSE
      NP1   = N+1
      P     = (MulTabLoc-MulTabMet(I,N))/(MulTabMet(I,NP1)-MulTabMet(I,N))
      CNA1  = CNA(I,N) + P * ( CNA(I,NP1) - CNA(I,N) )
      AOL1  = AOL(I,N) + P * ( AOL(I,NP1) - AOL(I,N) )
   END IF
ELSE
   CNA1  = CNA(I,1)
   AOL1  = AOL(I,1)
ENDIF

CNPOT(J,IBLADE) = CNA1 * (ALPHA - AOL1)
CNP(J,IBLADE)   = CNPOT(J,IBLADE)


ALPHA = MIN( MAX( ALPHA, AL(I,1) ), AL(I,NLIFT(I)) )
CALL LocateBin( ALPHA, AL(I,1:NLIFT(I)), I1, NLIFT(I) )

IF ( I1 == 0 ) THEN
   I1   = 1
   I1P1 = 2
   P1   = 0.0
ELSEIF ( I1 == NLIFT(I) ) THEN
   I1P1 = I1
   I1   = I1 - 1
   P1   = 1.0
ELSE
   I1P1 = I1 + 1
   P1   = ( AL(I,I1) - ALPHA ) / ( AL(I,I1) - AL(I,I1P1) )
ENDIF


IF ( NTables(I) > 1 ) THEN

 ! Locate the multiple airfoil table position in the table


   MulTabLoc = MIN( MAX( MulTabLoc, MulTabMet(I,1) ), MulTabMet(I,NTables(I)) )
   CALL LocateBin( MulTabLoc, MulTabMet(I,1:NTables(I)), I2, NTables(I) )

   IF ( I2 == 0 ) THEN
      I2P1 = 2
      I2   = 1
      P2   = 0.0
   ELSE IF( I2 == NTables(I) ) THEN
      I2P1 = I2
      I2   = I2 - 1
      P2   = 1.0
   ELSE
      I2P1 = I2 + 1
      P2 = (MulTabLoc-MulTabMet(I,I2))/(MulTabMet(I,I2P1)-MulTabMet(I,I2))
   ENDIF

 ! Interpolate the F-table values

   FSPB  = FTB( I,I1,I2P1) - (FTB( I,I1,I2P1) - FTB( I,I1P1,I2P1))*P1
   FSPCB = FTBC(I,I1,I2P1) - (FTBC(I,I1,I2P1) - FTBC(I,I1P1,I2P1))*P1
   FSPA  = FTB( I,I1,I2  ) - (FTB( I,I1,I2  ) - FTB( I,I1P1,I2  ))*P1
   FSPCA = FTBC(I,I1,I2  ) - (FTBC(I,I1,I2  ) - FTBC(I,I1P1,I2  ))*P1

   FSP( J,IBLADE) = FSPA  + P2 * ( FSPB - FSPA )
   FSPC(J,IBLADE) = FSPCA + P2 * ( FSPCB - FSPCA )

ELSE

   FSP( J,IBLADE) = FTB( I,I1,1) - ( FTB( I,I1,1) - FTB( I,I1P1,1) )*P1
   FSPC(J,IBLADE) = FTBC(I,I1,1) - ( FTBC(I,I1,1) - FTBC(I,I1P1,1) )*P1

ENDIF

IF ( ABS( AFE(J,IBLADE) - AOL1 ) < 1.E-10 ) THEN

   FSP(J,IBLADE)  = 1.0
   FSPC(J,IBLADE) = 1.0

ELSE

   TEMP = 2.*SQRT(ABS(FSP(J,IBLADE)/(AFE(J,IBLADE)-AOL1)))-1.
   FSP(J,IBLADE) = TEMP * TEMP * SIGN ( 1., TEMP )
   IF ( FSP(J,IBLADE) >  1.0 ) FSP(J,IBLADE) =  1.0
   IF ( FSP(J,IBLADE) < -1.0 ) FSP(J,IBLADE) = -1.0

   IF ( ABS( AFE(J,IBLADE) ) < 1.E-10 ) THEN
      FSPC(J,IBLADE) = 1.0
   ELSE
      TEMP = FSPC(J,IBLADE)/((AFE(J,IBLADE)-AOL1)*AFE(J,IBLADE))
      FSPC(J,IBLADE) = TEMP * TEMP * SIGN ( 1., TEMP )
      IF ( FSPC(J,IBLADE) >  1.0 ) FSPC(J,IBLADE) =  1.0
      IF ( FSPC(J,IBLADE) < -1.0 ) FSPC(J,IBLADE) = -1.0
   ENDIF

ENDIF

SRFP = SQRT( ABS( FSP(J,IBLADE) ) ) * SIGN( 1., FSP(J,IBLADE) ) + 1.
FK   = 0.25 * SRFP * SRFP
CVN(J,IBLADE) = CNPOT(J,IBLADE) * ( 1. - FK )



RETURN
END SUBROUTINE BEDINIT


 ! *****************************************************
   SUBROUTINE BedUpdate
 !  Update old Beddoes parameters at new time step
 ! *****************************************************


USE            Bedoes


ANE1    = ANE
ADOT1   = ADOT
OLDXN   = XN
OLDYN   = YN
CNPOT1  = CNPOT
OLDDPP  = DPP
FSP1    = FSP
FSPC1   = FSPC
OLDTAU  = TAU
OLDDF   = DF
OLDDFC  = DFC
OLDDN   = DN
OLDCNV  = CNV
CVN1    = CVN
CNP1    = CNP
CNPD1   = CNPD
OLDSEP  = BEDSEP
QX1     = QX
OLDDQ   = DQ
AFE1    = AFE
DQP1    = DQP
DFAFE1  = DFAFE



RETURN
END SUBROUTINE BedUpdate


 ! *****************************************************
   SUBROUTINE BEDDAT
 !  USED TO INPUT PARAMETERS FOR THE
 !   BEDDOES DYNAMIC STALL MODEL
 ! *****************************************************


USE                           Airfoil
USE                           Bedoes
USE                           Switch


IMPLICIT                      NONE


   ! Local Variables:

REAL(ReKi)                 :: ETA
REAL(ReKi)                 :: CA
REAL(ReKi)                 :: SA

INTEGER                    :: I
INTEGER                    :: J
INTEGER                    :: K



 ! Empirical constants for the Beddoes model

 ! TVL   = Non-dimensional time of transit for the
 !          vortex moving across the airfoil surface
 ! TP    = Time constant for pressure lag
 ! TV    = Time constant for strength of shed vortex
 ! TF    = Time constant applied to location of
 !          the separation point
 ! AS    = Speed of sound for Mach number calculation

TVL   = 11.0
TP    = 1.7
TV    = 6.0
TF    = 3.0
ETA   = .99  !bjj: this doesn't seem to be used for anything....

IF ( SIUNIT ) THEN
 ! SI UNITS--m/sec
   AS = 335.
ELSE
 ! ENGLISH UNITS--ft/sec
   AS = 1100.
ENDIF

 ! Generate table of F values from airfoil data table

DO J =1,NUMFOIL
   DO K =1,NTables(J)
      DO I = 1, NLIFT(J)

         CA = COS( AL(J,I) )
         SA = SIN( AL(J,I) )
         CN = CL(J,I,K) * CA + ( CD(J,I,K) - CDO(J,K) ) * SA
         CC = CL(J,I,K) * SA - ( CD(J,I,K) - CDO(J,K) ) * CA

         IF ( ABS( CNA(J,K) ) .GT. 1.E-6 ) THEN
            FTB(J,I,K)  = CN / CNA(J,K)
            FTBC(J,I,K) = CC / CNA(J,K)
         ELSE
            FTB(J,I,K)  = 1.0
            FTBC(J,I,K) = 1.0
         ENDIF

      END DO !I
   END DO !K
END DO !J

VOR   = .FALSE.
SHIFT = .FALSE.

BEDSEP = .FALSE.
ANE1   = 0.
OLDCNV = 0.
CVN1   = 0.
CNPOT1 = 0.
CNP1   = 0.
CNPD1  = 0.
OLDDF  = 0.
OLDDFC = 0.
OLDDPP = 0.
FSP1   = 0.
FSPC1  = 0.
TAU    = 0.
OLDTAU = 0.
OLDXN  = 0.
OLDYN  = 0.



RETURN
END SUBROUTINE BEDDAT


 ! *****************************************************
   SUBROUTINE BEDWRT
 !  USED TO OUTPUT PARAMETERS FOR THE
 !   BEDDOES DYNAMIC STALL MODEL
 ! *****************************************************

USE                           AD_IOParams
USE                           Airfoil
USE                           Bedoes

IMPLICIT                      NONE


   ! Local Variables:

INTEGER                    :: I
INTEGER                    :: K

CHARACTER(70)              :: Frmt

Frmt = '(3X,A, 21(:F8.4,3X) )'
DO K = 1, NTables(1)
   WRITE(UnADopt,'(/A/)') '  BEDDOES DYNAMIC STALL PARAMETERS:'
   WRITE(UnADopt, Frmt) 'CN SLOPE         ', ( CNA(I,K),      I = 1, NUMFOIL )
   WRITE(UnADopt, Frmt) 'STALL CN (UPPER) ', ( CNS(I,K),      I = 1, NUMFOIL )
   WRITE(UnADopt, Frmt) 'STALL CN (LOWER) ', ( CNSL(I,K),     I = 1, NUMFOIL )
   WRITE(UnADopt, Frmt) 'ZERO LIFT AOA    ', ( AOL(I,K)*R2D,  I = 1, NUMFOIL )
   WRITE(UnADopt, Frmt) 'MIN DRAG AOA     ', ( AOD(I,K)*R2D,  I = 1, NUMFOIL )
   WRITE(UnADopt, Frmt) 'MIN DRAG COEFF   ', ( CDO(I,K),      I = 1, NUMFOIL )
   WRITE(UnADopt,'(/)')
ENDDO !K

WRITE(UnADopt,*) '    VORTEX TRANSIT TIME FROM LE TO TE ', TVL
WRITE(UnADopt,*) '    PRESSURE TIME CONSTANT            ', TP
WRITE(UnADopt,*) '    VORTEX TIME CONSTANT              ', TV
WRITE(UnADopt,*) '    F-PARAMETER TIME CONSTANT         ', TF



RETURN
END SUBROUTINE BEDWRT


 ! ******************************************************
   SUBROUTINE BEDDOES( W2, J, IBlade, ALPHA, CLA, CDA ,CMA )
 !  uses the Beddoes dynamic stall model
 !   the routine is entered with an angle of attack
 !   and returns CL and CD.
 !  This routine is used regardless of whether the element
 !   is actually in dynamic stall state.
 !
 !  VARIABLES:
 !   W2    = Relative velocity squared over blade element
 !   J     = Index which identifies the blade element
 !   ALPHA = Angle of attack in radians
 !   CLA   = Lift coeff. which is calculated by the routine
 !   CDA   = Drag coeff. which is calculated by the routine
 !   CMA   = Moment coeff. which is calculated by the routine
 ! ******************************************************

USE                           Airfoil
USE                           Bedoes

IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: ALPHA
REAL(ReKi),INTENT(OUT)     :: CDA
REAL(ReKi),INTENT(OUT)     :: CLA
REAL(ReKi),INTENT(OUT)     :: CMA
REAL(ReKi),INTENT(IN)      :: W2

INTEGER,   INTENT(IN)      :: J
INTEGER,   INTENT(IN)      :: IBlade

   ! Local Variables:

REAL(ReKi)                 :: AE
REAL(ReKi)                 :: AOD1
REAL(ReKi)                 :: AOL1
REAL(ReKi)                 :: CA
REAL(ReKi)                 :: CDO1
REAL(ReKi)                 :: CNA1
REAL(ReKi)                 :: CNS1
REAL(ReKi)                 :: CNSL1
REAL(ReKi)                 :: P
REAL(ReKi)                 :: SA
REAL(ReKi)                 :: VREL

INTEGER                    :: I
INTEGER                    :: N
INTEGER                    :: NP1



 ! Check to see if element has multiple airfoil tables, then interpolate values
 !  of constants based on the current location.

I = NFOIL(J)

IF (NTables(I) > 1)THEN

   MulTabLoc = MIN( MAX( MulTabLoc, MulTabMet(I,1) ), MulTabMet(I,NTables(I)) )
   CALL LocateBin( MulTabLoc, MulTabMet(I,1:NTables(I)), N, NTables(I) )

   IF ( N == 0 ) THEN
      CNA1  = CNA( I,1)
      AOL1  = AOL( I,1)
      CNS1  = CNS( I,1)
      CNSL1 = CNSL(I,1)
      AOD1  = AOD( I,1)
      CDO1  = CDO( I,1)
   ELSE IF ( N == NTables(I) ) THEN
      CNA1  = CNA( I,N)
      AOL1  = AOL( I,N)
      CNS1  = CNS( I,N)
      CNSL1 = CNSL(I,N)
      AOD1  = AOD( I,N)
      CDO1  = CDO( I,N)
   ELSE
      NP1   = N+1
      P     = (MulTabLoc-MulTabMet(I,N))/(MulTabMet(I,NP1)-MulTabMet(I,N))

      CNA1  = CNA(I,N) + P * ( CNA(I,NP1) - CNA(I,N) )
      AOL1  = AOL(I,N) + P * ( AOL(I,NP1) - AOL(I,N) )
      CNS1  = CNS(I,N) + P * ( CNS(I,NP1) - CNS(I,N) )
      CNSL1 = CNSL(I,N)+ P * ( CNSL(I,NP1)- CNSL(I,N))
      AOD1  = AOD(I,N) + P * ( AOD(I,NP1) - AOD(I,N) )
      CDO1  = CDO(I,N) + P * ( CDO(I,NP1) - CDO(I,N) )
   END IF

ELSE
   CNA1  = CNA(I,1)
   AOL1  = AOL(I,1)
   CNS1  = CNS(I,1)
   CNSL1 = CNSL(I,1)
   AOD1  = AOD(I,1)
   CDO1  = CDO(I,1)
ENDIF

 ! Jump back if lift-curve slope is zero

IF ( CNA1 == 0.0 ) THEN
   CLA = 0.0
   CDA = CDO1
   RETURN
ENDIF

AN   = ALPHA
VREL = SQRT( W2 )

CALL ATTACH( VREL, J, IBlade, CNA1, AOL1, AE )

CALL SEPAR( NLIFT(I), J, IBlade, I, CNA1, AOL1, CNS1, CNSL1 )

CALL VORTEX( J, IBlade, AE )

CA  = COS( AN )
SA  = SIN( AN )
CLA = CN * CA + CC * SA
CDA = CN * SA - CC * CA + CDO1
CMA = PMC



RETURN
END SUBROUTINE BEDDOES


 ! ******************************************************
   SUBROUTINE ATTACH( VREL, J, IBlade, CNA1, AOL1, AE )
 !  PART OF THE BEDDOES DYNAMIC STALL MODEL.
 ! ******************************************************

USE                           AeroTime
USE                           Bedoes
USE                           Blade  !C


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(OUT)     :: AE
REAL(ReKi),INTENT(IN)      :: AOL1
REAL(ReKi),INTENT(IN)      :: CNA1
REAL(ReKi),INTENT(IN)      :: VREL

INTEGER,   INTENT(IN)      :: J
INTEGER,   INTENT(IN)      :: IBlade


   ! Local Variables:

REAL(ReKi)                 :: B2
REAL(ReKi)                 :: BS
REAL(ReKi)                 :: CNI
REAL(ReKi)                 :: CNQ
REAL(ReKi)                 :: CO
REAL(ReKi)                 :: DA
REAL(ReKi)                 :: PRP
REAL(ReKi)                 :: X
REAL(ReKi)                 :: XKA
REAL(ReKi)                 :: XM

LOGICAL,    SAVE           :: SuperSonic = .FALSE.


IF ( ABS( AN ) <= PIBY2 ) THEN
   ANE(J,IBLADE) = AN
ELSEIF ( AN > PiBy2 ) THEN
   ANE(J,IBLADE) = PI - AN
ELSE
   ANE(J,IBLADE) = - PI - AN
ENDIF


XM  = VREL / AS

 ! Check to see that the element is not supersonic
IF ( .NOT. SuperSonic .AND. XM >= 1.0 ) THEN
   XM = 0.7
   SuperSonic = .TRUE.

   CALL ProgWarn( ' Blade #'//TRIM(Int2LStr(IBLADE))//' element #'//TRIM(Int2LStr(J))//' is supersonic! '//&
                  ' Other elements are likely supersonic as well. Supersonic mach nos. will be set to '//&
                  TRIM(Flt2LStr(XM))//' to attempt continuation.' )
ELSEIF (SuperSonic .AND. XM < 1.0) THEN
   SuperSonic = .FALSE.
   CALL ProgWarn( ' Supersonic condition has subsided with Blade #'// TRIM(Int2LStr(IBLADE))// &
                  ' element #'//TRIM(Int2LStr(J))//'.')
ENDIF

B2  = 1.0 - XM * XM
DS  = 2. * DT * VREL/C(J)
BS  = B2 * DS
XKA = .75/( ( 1. - XM ) + PI * B2 * XM * XM * 0.413 )
X   = DT * AS / C(J) / XKA
CO  = XKA * C(J) / AS / XM

DA  = ANE(J,IBLADE) - ANE1(J,IBLADE)
ADOT(J,IBLADE) = DA / DT

PRP = ADOT(J,IBLADE) * C(J) / VREL
PRP = SAT( PRP, 0.03, 0.1 )
ADOT(J,IBLADE) = PRP * VREL / C(J)

DN(J,IBLADE) = OLDDN(J,IBLADE) * EXP(-X) + &
               (ADOT(J,IBLADE) - ADOT1(J,IBLADE)) * EXP(-.5*X)
CNI = 4. * CO * ( ADOT(J,IBLADE) - DN(J,IBLADE) )
CMI = -.25 * CNI

QX(J,IBLADE) = (ADOT(J,IBLADE) - ADOT1(J,IBLADE)) * C(J)/VREL/DT
DQ(J,IBLADE) = OLDDQ(J,IBLADE)*EXP(-X) + &
               ( QX(J,IBLADE) - QX1(J,IBLADE) ) * EXP(-.5*X)
CNQ = -CO * (QX(J,IBLADE) - DQ(J,IBLADE))
DQP(J,IBLADE) = DQP1(J,IBLADE) * EXP(-X/XKA) + (QX(J,IBLADE) &
                - QX1(J,IBLADE)) * EXP(-.5*X/XKA)

CMQ = -.25 * CNQ - (XKA*CO/3.) * (QX(J,IBLADE) - DQP(J,IBLADE))

CNIQ = MIN( ABS( CNI+CNQ ), 1. ) * SIGN( 1., CNI+CNQ )

XN(J,IBLADE) = OLDXN(J,IBLADE)*EXP(-.14*BS) + .3*DA*EXP(-.07*BS)
YN(J,IBLADE) = OLDYN(J,IBLADE)*EXP(-.53*BS) + .7*DA*EXP(-.265*BS)

AE   = ANE(J,IBLADE) - YN(J,IBLADE) - XN(J,IBLADE)
CNCP = CNA1 * ( AE - AOL1 )
CNPOT(J,IBLADE) = CNCP + CNIQ
CC   =  CNPOT(J,IBLADE) * AE



RETURN
END SUBROUTINE ATTACH


 ! ******************************************************
   SUBROUTINE SEPAR( NFT, J, IBlade, IFOIL, CNA1, AOL1, CNS1, CNSL1 )
 !  PART OF THE BEDDOES DYNAMIC STALL MODEL
 ! ******************************************************

USE                           Airfoil
USE                           Bedoes

IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: AOL1
REAL(ReKi),INTENT(IN)      :: CNA1
REAL(ReKi),INTENT(IN)      :: CNS1
REAL(ReKi),INTENT(IN)      :: CNSL1

INTEGER   ,INTENT(IN)      :: IFOIL
INTEGER,   INTENT(IN)      :: J
INTEGER,   INTENT(IN)      :: IBlade
INTEGER   ,INTENT(IN)      :: NFT

   ! Local Variables:

REAL(ReKi)                 :: AFEP
REAL(ReKi)                 :: AFF
REAL(ReKi)                 :: CMPA
REAL(ReKi)                 :: CMPB
REAL(ReKi)                 :: FSPA
REAL(ReKi)                 :: FSPB
REAL(ReKi)                 :: FSPCA
REAL(ReKi)                 :: FSPCB
REAL(ReKi)                 :: P1
REAL(ReKi)                 :: P2
REAL(ReKi)                 :: SRFP
REAL(ReKi)                 :: SRFPC
REAL(ReKi)                 :: TEMP
REAL(ReKi)                 :: TFE

INTEGER                    :: I1
INTEGER                    :: I1P1
INTEGER                    :: I2
INTEGER                    :: I2P1


TFE  = TF
DPP(J,IBLADE)   = OLDDPP(J,IBLADE) * EXP(-DS/TP) &
                 + ( CNPOT(J,IBLADE) - CNPOT1(J,IBLADE) ) * EXP(-.5*DS/TP)
CNP(J,IBLADE)  = CNPOT(J,IBLADE) -   DPP(J,IBLADE)
CNPD(J,IBLADE) = CNP(J,IBLADE) - CNP1(J,IBLADE)

 ! USE CNPD to determine if AOA is increasing or decreasing.
 ! Vortex lift decays more rapidly for decreasing AOA.

IF ( ANE(J,IBLADE) * CNPD(J,IBLADE) < 0. ) THEN
   SHIFT = .TRUE.
ELSE
   SHIFT = .FALSE.
ENDIF

AFF = CNP(J,IBLADE)/CNA1 + AOL1

IF ( ABS( AN ) <= PIBY2 ) THEN
   AFE(J,IBLADE) =  AFF
ELSEIF ( AN > PIBY2 ) THEN
   AFE(J,IBLADE) = PI -  AFF
ELSE
   AFE(J,IBLADE) = -PI -  AFF
ENDIF

CALL MPI2PI ( AFE(J,IBLADE) )

IF ( ( AFE(J,IBLADE) < AL(IFOIL,1) ) .OR.  &
     ( AFE(J,IBLADE) > AL(IFOIL,NLIFT(IFOIL) ) ) ) THEN !bjj compare w/ MIN(MAX()) below.  Is NLIFT(IFOIL)=NFT ? yes!!!
   CALL ProgAbort( ' Angle of attack = '//Flt2LStr(REAL( AFE(J,IBLADE)*R2D, ReKi ))// &
                   ' is outside table.' )
ENDIF

AFE(J,IBLADE) = MIN( MAX( AFE(J,IBLADE), AL(IFOIL,1) ), AL(IFOIL,NFT) )
CALL LocateBin( AFE(J,IBLADE), AL(IFOIL,1:NFT), I1, NFT )

IF (I1 == 0) THEN
   I1 = 1
ELSE IF ( I1 == NFT ) THEN
   I1 = I1 - 1
END IF

I1P1 = I1 + 1

P1 = ( AL(IFOIL,I1) - AFE(J,IBLADE) ) / ( AL(IFOIL,I1) - AL(IFOIL,I1P1) )

IF ( NTables(IFOIL) > 1 ) THEN

 ! Locate the multiple airfoil position in the table

   MulTabLoc = MIN( MAX( MulTabLoc, MulTabMet(IFOIL,1) ), MulTabMet(IFOIL,NTables(IFOIL)) )
   CALL LocateBin( MulTabLoc, MulTabMet(IFOIL,1:NTables(IFOIL)),I2,NTables(IFOIL) )

   IF ( I2 == 0 ) THEN
      I2   = 1
      I2P1 = 2
      P2   = 0.0
   ELSE IF ( I2 == NTables(IFOIL) ) THEN
      I2P1 = I2
      I2   = I2 - 1

      P2   = 1.0
   ELSE
      I2P1 = I2 + 1

      P2=(MulTabLoc-MulTabMet(IFOIL,I2))/(MulTabMet(IFOIL,I2P1)-MulTabMet(IFOIL,I2))
   END IF

 ! Interpolate the F-table values


   FSPB  = FTB( IFOIL,I1,I2P1) - (FTB( IFOIL,I1,I2P1) - FTB( IFOIL,I1P1,I2P1) ) * P1
   FSPCB = FTBC(IFOIL,I1,I2P1) - (FTBC(IFOIL,I1,I2P1) - FTBC(IFOIL,I1P1,I2P1) ) * P1
   FSPA  = FTB( IFOIL,I1,I2  ) - (FTB( IFOIL,I1,I2  ) - FTB( IFOIL,I1P1,I2  ) ) * P1
   FSPCA = FTBC(IFOIL,I1,I2  ) - (FTBC(IFOIL,I1,I2  ) - FTBC(IFOIL,I1P1,I2  ) ) * P1

   FSP(J,IBLADE) = FSPA  + P2 * (FSPB-FSPA)

   FSPC(J,IBLADE)= FSPCA + P2 * (FSPCB-FSPCA)

ELSE

   FSP(J,IBLADE) = FTB(IFOIL,I1,1) - (FTB(IFOIL,I1,1) - FTB(IFOIL,I1P1,1) ) * P1

   FSPC(J,IBLADE)= FTBC(IFOIL,I1,1) - (FTBC(IFOIL,I1,1) - FTBC(IFOIL,I1P1,1) ) * P1

ENDIF

IF ( ABS( AFE(J,IBLADE) - AOL1 ) < 1.E-10 ) THEN
   FSP(J,IBLADE)  = 1.0
   FSPC(J,IBLADE) = 1.0
ELSE
   TEMP = 2.*SQRT(ABS(FSP(J,IBLADE)/(AFE(J,IBLADE)-AOL1)))-1.
   FSP(J,IBLADE) = TEMP * TEMP * SIGN ( 1., TEMP )
   IF ( FSP(J,IBLADE) >  1.0 ) FSP(J,IBLADE) =  1.0
   IF ( FSP(J,IBLADE) < -1.0 ) FSP(J,IBLADE) = -1.0

   IF ( ABS( AFE(J,IBLADE) ) < 1.E-10 ) THEN
      FSPC(J,IBLADE) = 1.0
   ELSE
      TEMP = FSPC(J,IBLADE)/((AFE(J,IBLADE)-AOL1)*AFE(J,IBLADE))
      FSPC(J,IBLADE) = TEMP * TEMP * SIGN ( 1., TEMP )
      IF ( FSPC(J,IBLADE) >  1.0 ) FSPC(J,IBLADE) =  1.0
      IF ( FSPC(J,IBLADE) < -1.0 ) FSPC(J,IBLADE) = -1.0
   ENDIF

ENDIF

IF ( CNP(J,IBLADE) > CNS1  ) BEDSEP(J,IBLADE) = .TRUE.
IF ( CNP(J,IBLADE) < CNSL1 ) BEDSEP(J,IBLADE) = .TRUE.

IF ( BEDSEP(J,IBLADE) ) TAU(J,IBLADE) = OLDTAU(J,IBLADE) + DS/TVL

IF (SHIFT) TFE = 1.5*TFE

DF(J,IBLADE) = OLDDF(J,IBLADE) * EXP(-DS/TFE)  &
              + (FSP(J,IBLADE) - FSP1(J,IBLADE)) * EXP(-.5*DS/TFE)
DFC(J,IBLADE)= OLDDFC(J,IBLADE) * EXP(-DS/TFE) &
              + (FSPC(J,IBLADE) - FSPC1(J,IBLADE)) * EXP(-.5*DS/TFE)

FP   = FSP(J,IBLADE) - DF(J,IBLADE)
FPC  = FSPC(J,IBLADE) - DFC(J,IBLADE)
SRFP = SQRT( ABS(FP) )  * SIGN( 1., FP ) + 1.
SRFPC= SQRT( ABS(FPC) ) * SIGN( 1.,FPC )

FK   = 0.25 * SRFP * SRFP
CN   = CNCP * FK + CNIQ

CC   =  CC * SRFPC

DFAFE(J,IBLADE) = DFAFE1(J,IBLADE) * EXP(-DS/(.1*TFE)) &
                 + (AFE(J,IBLADE) - AFE1(J,IBLADE)) * EXP(-.5*DS/(.1*TFE))

AFEP=AFE(J,IBLADE) - DFAFE(J,IBLADE)


AFEP = MIN( MAX( AFEP, AL(IFOIL,1) ), AL(IFOIL,NFT) )
CALL LocateBin( AFEP, AL(IFOIL,1:NFT), I1, NFT )

IF (I1 == 0) THEN
   I1   = 1
   I1P1 = 2
   P1   = 0.0
ELSEIF ( I1 == NFT ) THEN
   I1P1 = I1
   I1   = I1 - 1
   P1   = 1.0
ELSE
   I1P1 = I1 + 1
   P1 = (AL(IFOIL,I1) - AFEP)/(AL(IFOIL,I1) - AL(IFOIL,I1P1))
END IF


IF (NTables(IFOIL) > 1) THEN

 ! Interpolate the F-table values

   CMPB = CM(IFOIL,I1,I2P1) - (CM(IFOIL,I1,I2P1) - CM(IFOIL,I1P1,I2P1) ) * P1

   CMPA = CM(IFOIL,I1,I2) - (CM(IFOIL,I1,I2) - CM(IFOIL,I1P1,I2) ) * P1

   PMC = CMPA  + P2*(CMPB-CMPA)

ELSE

   PMC = CM(IFOIL,I1,1) - ( (CM(IFOIL,I1,1) - CM(IFOIL,I1P1,1) ) * P1 )

ENDIF



RETURN
END SUBROUTINE SEPAR


 ! ******************************************************
   SUBROUTINE VORTEX( J, IBlade, AE )
 !  PART OF THE BEDDOES DYNAMIC STALL MODEL
 ! ******************************************************

USE                           Airfoil
USE                           Bedoes


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: AE

INTEGER,   INTENT(IN)      :: J
INTEGER,   INTENT(IN)      :: IBlade

   ! Local Variables:

REAL(ReKi)                 :: CMV
REAL(ReKi)                 :: TSH
REAL(ReKi)                 :: TVE



TVE = TV
CVN(J,IBLADE) = CNCP * ( 1. - FK )

IF ( TAU(J,IBLADE) < 1. ) THEN
   VOR = .TRUE.
   IF (SHIFT) VOR = .FALSE.
ELSE
   VOR = .FALSE.
ENDIF

IF (VOR) THEN
   CNV(J,IBLADE) = OLDCNV(J,IBLADE) * EXP(-DS/TVE)  &
                 + (CVN(J,IBLADE) - CVN1(J,IBLADE)) * EXP(-.5*DS/TVE)
ELSE
   CNV(J,IBLADE) = OLDCNV(J,IBLADE) * EXP(-DS/(TVE*.5))
ENDIF

CN  = CN + CNV(J,IBLADE)
CC  = CC + CNV(J,IBLADE) * AE * (1.- TAU(J,IBLADE))
CMV = -0.2 * (1. - COS(PI*TAU(J,IBLADE)) ) * CNV(J,IBLADE)
PMC = PMC + CMV + CMI + CMQ

TSH = 2.*(1.- FP)/.19

IF ( TAU(J,IBLADE) .GT. 1. + TSH/TVL .AND. .NOT. SHIFT) THEN
   TAU(J,IBLADE) = 0.
   BEDSEP(J,IBLADE) = .FALSE.
ENDIF

IF ( TAU(J,IBLADE) .GT. 1. ) THEN
   IF ( ANE(J,IBLADE) .LT. 0. ) THEN

      IF (CNPD(J,IBLADE) .LE. 0. .AND. CNPD1(J,IBLADE) .GE. 0.) THEN
         BEDSEP(J,IBLADE) = .FALSE.
         TAU(J,IBLADE) = 0.
      ENDIF

      IF (ANE1(J,IBLADE) .GT. 0.) THEN
         BEDSEP(J,IBLADE) = .FALSE.
         TAU(J,IBLADE) = 0.
      ENDIF

   ELSE

      IF (CNPD(J,IBLADE) .GE. 0. .AND. CNPD1(J,IBLADE) .LE. 0.) THEN
         BEDSEP(J,IBLADE) = .FALSE.
         TAU(J,IBLADE) = 0.
      ENDIF

      IF (ANE1(J,IBLADE) .LT. 0.) THEN
         BEDSEP(J,IBLADE) = .FALSE.
         TAU(J,IBLADE) = 0.
      ENDIF

   ENDIF
ENDIF



RETURN
END SUBROUTINE VORTEX


 ! ******************************************************
   SUBROUTINE CLCD( ALPHA, CLA, CDA, CMA, I, ErrStat )
 !  returns values of lift and drag coeffs.
 !   This subroutine interpolates airfoil coefficients
 !   from a table of airfoil data.  The table must consist
 !   of ALPHA, CL and CD over the entire range of angles
 !   that will be encountered.
 !
 ! VARIABLES:
 !    CLA      = Returned value of lift coefficient
 !    CDA      = Returned value of drag coeff
 !    CMA      = Returned value of pitching moment coeff
 !    ALPHA    = Angle of attack (radians)
 !    AL       = Array containing the angle of attack
 !    CL       = Array containing the lift coeffs. at AL(I)
 !    CD       = Array containing the drag coeffs. at AL(I)
 !    CM       = Array containing the moment coeffs. at AL(I)
 !    I        = Airfoil ID for this element, equal to NFoil(J), where J is the index identifying the blade element
 !    MulTabLoc= Multiple airfoil table location for this element
 !    MulTabMet= Array containing the multiple airfoil table metric
 ! ******************************************************


USE                           Airfoil

IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(INOUT)   :: ALPHA
REAL(ReKi),INTENT(OUT)     :: CDA
REAL(ReKi),INTENT(OUT)     :: CLA
REAL(ReKi),INTENT(OUT)     :: CMA

INTEGER   ,INTENT(IN)      :: I      ! NFOIL(J)
INTEGER,   INTENT(OUT)     :: ErrStat

   ! Local Variables:

REAL(ReKi)                 :: CDA1
REAL(ReKi)                 :: CDA2
REAL(ReKi)                 :: CLA1
REAL(ReKi)                 :: CLA2
REAL(ReKi)                 :: CMA1
REAL(ReKi)                 :: CMA2
REAL(ReKi)                 :: P1
REAL(ReKi)                 :: P2

INTEGER                    :: N1
INTEGER                    :: N1P1
INTEGER                    :: N2
INTEGER                    :: N2P1
INTEGER                    :: NTAB


IF (.NOT. ALLOCATED(NFoil) ) THEN
   CDA = 0
   CLA = 0
   CMA = 0
   ErrStat = 1
   RETURN
ELSE
   ErrStat = 0
END IF

NTAB = NLIFT(I)

IF ( ( ALPHA < AL(I,1) ) .OR. ( ALPHA > AL(I,NTAB) ) )   THEN
!bjj: This error message isn't necessarially accurate:
   CALL ProgAbort( ' Angle of attack = '//TRIM(Flt2LStr(ALPHA*R2D))// &
                   ' deg is outside data table range. '// & !Blade #'//TRIM(Int2LStr(IBLADE))//&
                   ' Airfoil '//TRIM(Int2LStr(I))//'.' )
!                   ' element '//TRIM(Int2LStr(J))//'.' )

   ErrStat = 1
ENDIF

ALPHA = MIN( MAX( ALPHA, AL(I,1) ), AL(I,NTAB) )
CALL LocateBin (ALPHA, AL(I,1:NTAB), N1, NTAB )

IF (N1 == 0) THEN
   N1   = 1
   N1P1 = 2
   P1   = 0.0
ELSEIF(N1 == NTAB) THEN
   N1P1 = N1
   N1   = N1 - 1
   P1   = 1.0
ELSE
   N1P1 = N1 + 1
   P1   = ( ALPHA - AL(I, N1) )/( AL(I, N1P1) - AL(I, N1) )
END IF

 ! If the element has multiple airfoil tables, do a 2-D linear interpolation
 !  for Cl and CD

IF (NTables(I) > 1) THEN

   MulTabLoc = MIN( MAX( MulTabLoc, MulTabMet(I,1) ), MulTabMet(I,NTables(I)) )
   CALL LocateBin (MulTabLoc, MulTabMet(I,1:NTables(I)),N2,NTables(I))

   IF (N2 == 0) THEN
      N2   = 1
      N2P1 = 2
      P2   = 0.0
   ELSE IF ( N2 == NTables(I) ) THEN
      N2P1 = N2
      N2   = N2 - 1
      P2   = 1.0
   ELSE
      N2P1 = N2 + 1
      P2   = (MulTabLoc - MulTabMet(I,N2))/(MulTabMet(I,N2P1)-MulTabMet(I,N2))
   END IF

   CLA1 = CL(I,N1,N2) + P1 * ( CL(I,N1P1,N2) - CL(I,N1,N2) )
   CDA1 = CD(I,N1,N2) + P1 * ( CD(I,N1P1,N2) - CD(I,N1,N2) )
   CMA1 = CM(I,N1,N2) + P1 * ( CM(I,N1P1,N2) - CM(I,N1,N2) )

   CLA2 = CL(I,N1,N2P1) + P1 * ( CL(I,N1P1,N2P1) - CL(I,N1,N2P1) )
   CDA2 = CD(I,N1,N2P1) + P1 * ( CD(I,N1P1,N2P1) - CD(I,N1,N2P1) )
   CMA2 = CM(I,N1,N2P1) + P1 * ( CM(I,N1P1,N2P1) - CM(I,N1,N2P1) )

   CLA = CLA1 + P2 * ( CLA2 - CLA1 )
   CDA = CDA1 + P2 * ( CDA2 - CDA1 )
   CMA = CMA1 + P2 * ( CMA2 - CMA1 )

ELSE

   CLA  = CL(I,N1,1) + P1 * ( CL(I,N1P1,1) - CL(I,N1,1) )
   CDA  = CD(I,N1,1) + P1 * ( CD(I,N1P1,1) - CD(I,N1,1) )
   CMA  = CM(I,N1,1) + P1 * ( CM(I,N1P1,1) - CM(I,N1,1) )

ENDIF


RETURN
END SUBROUTINE CLCD

 ! **************************************************
   FUNCTION SAT( X, VAL, SLOPE )
 !  AOA saturation function 02/15/98
 ! **************************************************

IMPLICIT                      NONE


   ! Passed Variables:

REAL(ReKi)                 :: SAT
REAL(ReKi),INTENT(IN)       :: SLOPE
REAL(ReKi),INTENT(IN)       :: VAL
REAL(ReKi),INTENT(IN)       :: X


IF ( ABS(X) <= VAL )  THEN
    SAT = X
ELSEIF ( X > VAL)  THEN
    SAT = SLOPE * X + VAL * ( 1. - SLOPE )
ELSE
    SAT = SLOPE * X - VAL * ( 1. - SLOPE )
ENDIF


RETURN
END FUNCTION SAT


 ! **************** Dynamic Inflow Subroutines ***************
 !
 !  Generalized Dynamic Wake (GDW) Model replaced the Pitt &
 !   Peters model.
 !   A. Suzuki, 06/23/00

 ! Subroutine VG2ROTOR was added again to v11.31
 !  A. Suzuki, 01/24/00

 !
 !  Modified for FFWIND. Subroutine VG2ROTOR was added.
 !   A. Suzuki, 11/05/99.

 !  This model is based on Pitt & Peters model.
 !   AB predictor-corrector is used for integration.
 !   A. Suzuki, 07/22/98.
 !
 ! *************************************************************

 ! **************************************

   SUBROUTINE Inflow()

 ! Gateway to the dynamic inflow routines.
 ! Called by NEWTIME after a time step.
 ! **************************************


USE               AeroTime
USE               DynInflow
USE               Switch


IMPLICIT          NONE

      IF ( DYNINFL ) THEN

 ! INITIALIZE DYNAMIC INFLOW PARAMETERS
         IF ( DYNINIT .AND. ( TIME > 0.0D0 ) ) THEN
            CALL INFINIT
            !WRITE(*,*) 'Activating dynamic inflow calculation'
            old_Alph = 0.0
            old_Beta = 0.0
            DYNINIT = .FALSE.
            SKEW = .FALSE.
         ENDIF

 ! Update the dynamic inflow parameters
         CALL INFUPDT
         CALL GetPhiLq
 ! Calculate the dynamic inflow paremeters for the new time step
         IF( TIME > 1.0D0 ) CALL INFDIST

      ENDIF   ! DYNINFL

RETURN
END SUBROUTINE Inflow


 ! **************************************

   SUBROUTINE GetRM (rLocal, DFN, DFT, psi, J, IBlade)

 ! Returns RM(MODE), the [mode]-th moment of the blade normal force.
 !  Here, the force is in [N/m], while the moment arm is
 !  non-dimensional (RLOCAL/R).  Also see FUNCTION XPHI.
 ! Called as each element is processed by AeroDyn.
 ! **************************************


USE                           Blade
USE                           DynInflow
USE                           Switch


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: DFN
REAL(ReKi),INTENT(IN)      :: DFT
REAL(ReKi),INTENT(IN)      :: psi
REAL(ReKi),INTENT(IN)      :: rLocal

INTEGER,   INTENT(IN)      :: J
INTEGER,   INTENT(IN)      :: IBlade


   ! Local Variables:

REAL(ReKi)                 :: fElem
! psiBar is Suzuki's, WindPsi is Shawler's
!REAL(ReKi)                :: psiBar
REAL(ReKi)                 :: Rzero
REAL(ReKi)                 :: WindPsi

INTEGER                    :: mode



IF ( SWIRL ) THEN
   fElem = SQRT( DFN * DFN + DFT * DFT )
ELSE
   fElem = DFN
ENDIF
fElem = fElem / R

Rzero = rLocal / R
! Suzukis inflow azimuth measure
!psiBar = - psi - piBy2
! Shawler: wind based inflow azimuth measure
CALL WindAzimuthZero (psi,WindPsi)

! Save values rotor loads for USE in Newtime (to accumulate rotor loads)
DO mode = 1, MAXINFL0
   RMC_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode )
END DO ! mode

!+++++++++++++++++++++++++++++++++++++++++++++++++++++
!Suzuki's method
!DO mode = MAXINFL0+1, maxInfl
!   RMC_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode )  * COS( REAL(MRvector(mode)) * psiBar )
!   RMS_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode )  * SIN( REAL(MRvector(mode)) * psiBar )
!END DO ! mode
! Shawler's method
DO mode = MAXINFL0+1, maxInfl
   RMC_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode )  * COS( REAL(MRvector(mode)) * Windpsi )
   RMS_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode )  * SIN( REAL(MRvector(mode)) * Windpsi )
END DO ! mode
!++++++++++++++++++++++++++++++++++++++++++++++++++++++



RETURN
END SUBROUTINE GetRM


 ! **************************************

   SUBROUTINE GetPhiLq

 ! Accumulate the rotor forces for dynamic inflow calculations
 !  PhiLqC is Lq times PHI (shape function) in COS equation.
 !  PhiLqS is Lq times PHI (shape function) in SIN equation.
 !  RM?_SAVE, which were calculated for each element,
 !   are summed here after a time step.
 ! **************************************

USE            Blade
USE            DynInflow
USE            Element


IMPLICIT       NONE


INTEGER     :: iblad
INTEGER     :: ielem
INTEGER     :: mode



PhiLqC = 0.
PhiLqS = 0.

DO mode = 1, maxInfl
   DO iblad = 1, NB
      DO ielem = 1, NELM
         PhiLqC(mode) = PhiLqC(mode) + RMC_SAVE(iblad, ielem, mode)
         IF (mode >= maxInfl0+1) &
         PhiLqS(mode) = PhiLqS(mode) + RMS_SAVE(iblad, ielem, mode)
      END DO !ielem
   END DO !iblad
END DO !mode



RETURN
END SUBROUTINE GetPhiLq


 ! *************************************************************
   SUBROUTINE infinit
 !  Initializes the variables in the dynamic inflow equation.
 !  Called only once to initialize the GDW parameters.
 ! *************************************************************


USE                           AeroTime
USE                           Blade
USE                           DynInflow
USE                           Element
USE                           Rotor
USE                           Wind


IMPLICIT                      NONE


   ! Local Variables:

REAL(ReKi)                 :: tauCos ( maxInfl )
REAL(ReKi)                 :: tauSin ( maxInfl0+1 : maxInfl )
REAL(ReKi)                 :: v1
REAL(ReKi)                 :: v2
REAL(ReKi)                 :: v3
REAL(ReKi)                 :: Vplane2

INTEGER                    :: i
INTEGER                    :: iElem
INTEGER                    :: irow
INTEGER                    :: jBlade
INTEGER                    :: jcol
INTEGER                    :: k
INTEGER                    :: mode


 ! Initialize the MRvector & NJVector.
 !  MRvector is the azimuthal mode of the inflow distribution.
 !  NJvector is the radial mode of the inflow distribution.

MRVector(:) = (/ 0, 0, 1, 1, 2, 3 /)  !bjj why aren't these parameters?
NJVector(:) = (/ 1, 3, 2, 4, 3, 4 /)

 ! For your reference,
 !  Marray(irow,jcol) = MRvector(jcol)
 !  Rarray(irow,jcol) = MRvector(irow)
 !  Narray(irow,jcol) = NJvector(jcol)
 !  Jarray(irow,jcol) = NJvector(irow)

 ! Initialize the time derivatives.
dAlph_dt(:,:) = 0.
dBeta_dt(:,:) = 0.

 ! Set up the constants.
 !  xMinv is [M]^-1.  Because [M]^-1 is just a diagonal matrix,
 !  it is stored as a column vector.
DO irow = 1, maxInfl
   xMinv(irow) = PIBY2 / hfunc(MRvector(irow), NJvector(irow))   !bjj: this is really just a parameter, too.
END DO !irow

 ! Set up the GAMMA matrix which is used to calculate [L] matrix.
 !  FUNCTION FGAMMA is called.
DO irow = 1, maxInfl
   DO jcol = 1, maxInfl
      gamma( irow, jcol )  &
           = fgamma( MRvector( irow ), NJvector( irow ),  &
                     MRvector( jcol ), NJvector( jcol )  )        !bjj: this is really just a parameter, too.
   END DO !jcol
END DO !irow

 ! calculate and store the M-R matrices, which are used in Subroutine LMATRIX.
DO irow = 1, maxInfl
   DO jcol = 1, maxInfl
      MminR  (irow,jcol) = MIN( MRvector(jcol) , MRvector(irow) ) !bjj: this is really just a parameter, too.
      MplusR (irow,jcol) =      MRvector(jcol) + MRvector(irow)   !bjj: this is really just a parameter, too.
      MminusR(irow,jcol) = ABS( MRvector(jcol) - MRvector(irow) ) !bjj: this is really just a parameter, too.
   END DO !jcol
END DO !irow

 ! Calculate the tip speed of the rotor. This isn't constant in ADAMS.
 !  Thus, it will be updated at every time step in ADAMS.
TipSpeed = MAX(r * revs, 1.0e-6)

 ! Calculate the disk loading normalization factor.
 !  This is not exactly pressure but let's call it P0.
 !  The actual unit is [N/m] or [Pa*m].
 !    Pzero = PI * AirDensity * (Rotational Speed)^2 * (Radius)^3
 !    Pzero = pi * rho * revs**2 * r**3
 !    Pzero = pi * rho * revs * revs * r * r * r
Pzero = pi * rho * TipSpeed * TipSpeed * r
 ! Non-dimensional time
DT0   = DT * REVS !bjj: this isn't used in this subroutine?

 ! Calculate the initial values of inflow distribution parameters

 ! Calculate the non-dimensional wind velocity components.

v1 =   VrotorZ / TipSpeed     !inplane, upward
v2 =   VrotorY / TipSpeed     !inplane, right looking downwind
v3 = - VrotorX / TipSpeed     !out-of-plane, normal to the rotor

 ! Calculate the initial value of lambda_m by taking the average
 !  of the induction factors(A). The A's are calculated by
 !  momentum balance during the first rotation of the trim solution.
xLambda_M = 0.
DO jBlade=1,nb
   DO iElem =1,nelm
      xLambda_M = xLambda_M + a(iElem,jBlade)
   END DO !iElem
END DO !jBlade
xLambda_M = xLambda_M / ( nb * nelm )

 ! A's are normalized by the normal wind speed, while Lambda's are
 !  mormalized by the tip speed. Make the conversion.
 !  xLambda_M = xLambda_M * (-VrotorX/TipSpeed)
xLambda_M = xLambda_M * v3

 ! totalInf is the non-dimensional total wind speed normal to the rotor.
totalInf = - v3 + xLambda_M
 ! Vplane2 is the square of in-plane component of the non-dimensional
 !  wind velocity.
Vplane2 = v1 * v1 + v2 * v2
 ! VTOTAL is the total wind speed at the rotor.
Vtotal  = SQRT( totalInf * totalInf + Vplane2 )
 ! VPARAM is the velocity parameter.
Vparam  =( Vplane2 + ( totalInf + old_LmdM ) * totalInf ) / Vtotal

 ! Calculate the disk skew angle function using the effective disk angle
 !  of attack.
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! Shawler - USE the single wake skew angle, squared variable means its always positive
 ! and because blade azimuth position is measured using windpsi the wake skew
 ! is defined in the right direction relative to the oncoming wind,
 ! ie directly downwind.
 !Suzuki:
!xKaiC = TAN( .5 * ATAN( -v2 / totalInf ) )
!xKaiS = TAN( .5 * ATAN(  v1 / totalInf ) )
 !xkai = TAN( .5 * SIGN( ATAN( SQRT( vplane2 ) / totalInf ), v2 ) )
 !Shawler:
xkai = TAN( .5 * ATAN( SQRT( Vplane2 ) / totalInf ) )
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ! To calculate the initial values of xAlpha & xBeta, get [L] matrices
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! Suzuki:
!CALL LMATRIX ( xKaiC, 1 )
!CALL LMATRIX ( xKaiS, 2 )
 ! Shawler:
CALL LMATRIX ( xKai, 1 )
CALL LMATRIX ( xKai, 2 )
! CALL LMATRIX ( xkai )
 ! Here we need [L_cos] & [L_sin], not [L_***}^-1.
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ! Get TAU's (pressure coefficient). Refer to Subroutine INFDIST.
DO mode = 1, maxInfl0
   tauCos(mode) = - PhiLqC(mode) / Pzero * .5
END DO !mode

DO mode = maxInfl0+1, maxInfl
   tauCos(mode) = - PhiLqC(mode) / Pzero
   tauSin(mode) = - PhiLqS(mode) / Pzero
END DO !mode

 ! Get the steady values of alpha(1)
 !  If m=0 and n=1, USE VTOTAL.
xAlpha(1) = 0.
DO k = 1, maxInfl
   xAlpha(1) = xAlpha(1) + xLcos(1,k) * tauCos(k)
END DO !k
xAlpha(1) = .5 * xAlpha(1) / Vtotal

 ! If m=0 but NOT n=1, USE VPARAM.
DO i = 2, maxInfl0
   xAlpha(i) = 0.
   DO k = 1, maxInfl
      xAlpha(i) = xAlpha(i) + xLcos(i,k) * tauCos(k)
   END DO !k
   xAlpha(i) = .5 * xAlpha(i) / Vparam
END DO !i

 ! Get the steady values of alpha's & beta's
DO i = maxInfl0+1, maxInfl
   xAlpha(i) = 0.
 ! akihiro
!   DO k = 1, maxInfl
   DO k = maxInfl0+1, maxInfl
      xAlpha(i) = xAlpha(i) + xLcos(i,k) * tauCos(k)
   END DO !k
   xAlpha(i) = .5 * xAlpha(i) / Vparam

   xBeta (i) = 0.
   DO k = maxInfl0+1, maxInfl
      xBeta (i) = xBeta (i) + xLsin(i,k) * tauSin(k)
   END DO !k
   xBeta (i) = .5 * xBeta (i) / Vparam
END DO !i

 ! Invert [L_cos] & [L_sin] matrices in case the XKAI is constant
 !  and the same [L]'s are used later.
CALL MATINV  ( xLcos, xLsin, maxInfl, maxInfl0, 1 )
CALL MATINV  ( xLcos, xLsin, maxInfl, maxInfl0, 2 )



RETURN
END SUBROUTINE infinit


 ! **********************************************************************
   SUBROUTINE infupdt
 !  INFUPDT updates the OLD variables of inflow distribution parameters.
 !   The program must call this subroutine before calling
 !   the subroutine INFDIST.
 ! **********************************************************************


USE            DynInflow


IMPLICIT       NONE


INTEGER     :: i
!rm not used:INTEGER(4)  :: mode



!+++++++++++++++++++++++++++
!Suzuki:
!oldKaiC     = xKaiC
!oldKaiS     = xKaiS
!Shawler:
oldKai = xKai
!+++++++++++++++++++++++++++

old_LmdM    = xLambda_M

DO i = 1, maxInfl
   old_Alph(i)   = xAlpha  (i)
   dAlph_dt(i,4) = dAlph_dt(i,3)
   dAlph_dt(i,3) = dAlph_dt(i,2)
   dAlph_dt(i,2) = dAlph_dt(i,1)
END DO !i

DO i = maxInfl0+1, maxInfl
   old_Beta(i)   = xBeta   (i)
   dBeta_dt(i,4) = dBeta_dt(i,3)
   dBeta_dt(i,3) = dBeta_dt(i,2)
   dBeta_dt(i,2) = dBeta_dt(i,1)
END DO !i



RETURN
END SUBROUTINE infupdt


 ! **********************************************************************
   SUBROUTINE DynDebug (RHScos, RHSsin)
 !  Write out debugging information
 ! **********************************************************************

USE                           AeroTime
USE                           DynInflow


IMPLICIT                      NONE


   ! Passed Variables:

REAL(ReKi)                 :: RHScos ( maxInfl )
REAL(ReKi)                 :: RHSsin ( maxInfl0+1:maxInfl )


   ! Local Variables:

INTEGER                    :: i
INTEGER                    :: NumOut
INTEGER                    :: UnDyn = 80

LOGICAL                    :: OnePass = .TRUE.

CHARACTER(50)              :: Frmt

SAVE                                                   ! Save *all* local variables.  Is this necessary, or is OnePass enough.


NumOut = maxInfl+(maxInfl-maxInfl0) + 1

IF (OnePass) THEN

   CALL OpenFOutFile (UnDyn, 'DynDebug.plt')

   Frmt = '( A4,    (: A1, A, I1.1 ) )'

   WRITE(Frmt(7:9), '(I3)') NumOut
   WRITE(UnDyn, Frmt) 'Time',                    &
            ( TAB,    'dAlph_dt',       i,         &
                       i = 1, maxInfl ),         &
            ( TAB,    'dBeta_dt',       i,         &
                     i = maxInfl0+1, maxInfl ),  &
              TAB,    'TotalInf'

   OnePass = .FALSE.

ENDIF

Frmt = '( F10.3,    ( : A1, ES12.5 ) )'

IF (TIME > 0.0D0) THEN

   WRITE(Frmt(10:12), '(I3)') NumOut
   WRITE(UnDyn,Frmt) TIME,                    &
        ( TAB,       dAlph_dt(i,1),               &
                   i = 1, maxInfl ),          &
        ( TAB,       dBeta_dt(i,1),               &
                   i = maxInfl0+1, maxInfl ), &
              TAB,   totalInf

ENDIF



RETURN
END SUBROUTINE DynDebug


 ! **********************************************************************
   SUBROUTINE infdist
 !  INFDIST calculates the inflow (induced flow) distribution
 !  parameters using Generalized Dynamic Wake Theory.
 ! **********************************************************************

USE                           AeroTime
USE                           Blade
USE                           DynInflow
USE                           Rotor
USE                           Wind


IMPLICIT                      NONE


   ! Local Variables:

REAL(ReKi)                 :: RHScos ( maxInfl )                     ! The cosine terms go from 1 to 6.  The sine terms go from 3 to 6.
REAL(ReKi)                 :: RHSsin ( maxInfl0+1:maxInfl )          ! The right hand side of the governing equation.
REAL(ReKi)                 :: tauCos ( maxInfl )                     ! Forcing Functions
REAL(ReKi)                 :: tauSin ( maxInfl0+1:maxInfl )
REAL(ReKi)                 :: v1
REAL(ReKi)                 :: v2
REAL(ReKi)                 :: v3
REAL(ReKi)                 :: Vplane2

INTEGER                    :: i
INTEGER                    :: k
INTEGER                    :: mode



      TipSpeed = MAX(r  * revs, 1.0e-6)   !bjj: why is this here AND in InfInit()?

      Pzero    = pi * rho * TipSpeed * TipSpeed * r

 ! Non-dimensional time
      DT0      = DT * REVS

 ! Calculate the wind velocity components in rotor-fixed
 !  coordinates(1-2-3), which are normalized by the tipspeed.
v1 =   VrotorZ / TipSpeed     !inplane, upward
v2 =   VrotorY / TipSpeed     !inplane, right looking downwind
v3 = - VrotorX / TipSpeed     !out-of-plane (normal to the rotor)
 ! Vplane2 is the square of in-plane component of the non-dimensional
 !  wind velocity.
Vplane2  = v1 * v1 + v2 * v2

 ! Note: Direction of non-dimensional velocity (All normal to the rotor plane).
 !  totalInf:  positive downwind. This is the total inflow to the rotor.
 !        v3:        positive upwind (thus, in normal condition, v3 < 0 )
 ! xLambda_M: positive downwind (opposite to A(i,j) )
 !  old_LmdM:  positive downwind (opposite to A(i,j) )
totalInf = - v3 + old_LmdM
! if ( totalInf .le. 0. ) then
!     call usrmes( .true. , &
!        'In SUBROUTINE INFDIST. totalInf =< 0.', &
!        27, 'WARN' )
! endif

 !  VTOTAL is the speed of the total inflow to the rotor.
Vtotal = SQRT( totalInf * totalInf + Vplane2 )
IF (vtotal <= 1.0e-6) vtotal=1.0e-6

 ! VPARAM is the inflow velocity parameter.
Vparam = ( Vplane2 + ( totalInf + old_LmdM ) * totalInf ) / Vtotal
 ! Calculate the disk skew angle function
 !  using the effective disk angle of attack.
IF ( totalInf == 0. ) THEN
!   WRITE(*,*) v3, old_LmdM
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Suzuki:
!   xKaiC = 0.
!   xKaiS = 0.
! Shawler:
    xKai = 0
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ELSE
 ! Note the definition of Yaw Angle is around -Z axis.
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Suzuki:
!   xKaiC = TAN( .5 * ATAN( -v2 / totalInf ) )
!   xKaiS = TAN( .5 * ATAN(  v1 / totalInf ) )
!!   xKaiC = TAN( .5 * SIGN( ATAN( SQRT( vplane2 ) / totalInf ), v2 ) )
! Shawler:
   xkai  = TAN( .5 *       ATAN( SQRT( vplane2 ) / totalInf )       )
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ENDIF

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Suzuki:
 ! Update [L_cos] and [L_sin] matrices only if 'xkai' has changed.
 !  Then invert [L_cos] & [L_sin] matrices.
!IF ( xKaiC /= oldKaiC ) THEN
!   CALL LMATRIX ( xKaiC, 1 )
!   CALL MATINV  ( xLcos, xLsin, maxInfl, maxInfl0, 1 )
!ENDIF

!IF ( xKaiS /= oldKaiS ) THEN
!   CALL LMATRIX ( xKaiS, 2 )
!   CALL MATINV  ( xLcos, xLsin, maxInfl, maxInfl0, 2 )
!ENDIF
! Shawler:
!IF ( xKai /= oldKai ) THEN
   CALL LMATRIX ( xKai, 1 )
   CALL LMATRIX ( xKai, 2 )
   CALL MATINV  ( xLcos, xLsin, maxInfl, maxInfl0, 1 )
   CALL MATINV  ( xLcos, xLsin, maxInfl, maxInfl0, 2 )
!ENDIF
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 !  [L_***] is now [L_***]^-1.

 ! Calculate the forcing functions (lift or normal forces).
 !  In the generalized dynamic wake model, the normal forces on the blades
 !  are positve upwind (just like a helicopter rotor). On the other hand,
 !  they are positve downwind in YawDyn, which is normal in wind turbine rotor.
 !  Put a minus sign for each forcing function.

 ! The modes for r=0.
DO mode = 1, maxInfl0
   tauCos(mode) = - PhiLqC(mode) / Pzero * .5
END DO !mode

 ! The modes for r>0.
DO mode = maxInfl0+1, maxInfl
   tauCos(mode) = - PhiLqC(mode) / Pzero
   tauSin(mode) = - PhiLqS(mode) / Pzero
END DO !mode

 ! Solve for the time derivatives of xAlpha's, {d(alpha)_dt}.
 !  Calculate the right hand side of the governing equation.
 !  {rhs} = 0.5*{tau} - [V][L]^-1*{alpha}

 ! First, calculate {rhs} = V[L]^-1*{alpha}
 !  USE "VTOTAL" for the first row of r=0. Cosine Matrix only.
RHScos(1) = 0.
DO k = 1, maxInfl
   RHScos(1) = RHScos(1) + xLcos(1,k) * old_Alph(k)
END DO !k
RHScos(1) = .5 * tauCos(1) - vtotal * RHScos(1)

 ! USE "VPARAM" for the second row or higher of r=0.
DO i = 2, maxInfl0
   RHScos(i) = 0.
   DO k = 1, maxInfl
      RHScos(i) = RHScos(i) + xLcos(i,k) * old_Alph(k)
   END DO !k
   RHScos(i) = .5 * tauCos(i) - Vparam * RHScos(i)
END DO !i

 ! Rows with r=1 or greater. Both cosine and sine matrices.
DO i = maxInfl0+1, maxInfl
   RHScos(i) = 0.
   RHSsin(i) = 0.
 ! First, calculate {rhs} = V[L]^-1*{alpha}
   DO k = 1, maxInfl
!   DO 260 k = maxInfl0+1, maxInfl
      RHScos(i) = RHScos(i) + xLcos(i,k) * old_Alph(k)
   END DO !k
   DO k = maxInfl0+1, maxInfl
      RHSsin(i) = RHSsin(i) + xLsin(i,k) * old_Beta(k)
   END DO !k
 ! Second, calculate {rhs} = 0.5*{tau} - [V]{rhs}
 !                         = 0.5*{tau} - [V][L]^-1*{alpha}
 !  USE "VPARAM" for m>0
   RHScos(i) = .5 * tauCos(i) - Vparam * RHScos(i)
   RHSsin(i) = .5 * tauSin(i) - Vparam * RHSsin(i)
END DO !i

DO i = 1, maxInfl0
   dAlph_dt(i,1) = xMinv(i) * RHScos(i)
END DO !i
DO i = maxInfl0+1, maxInfl
   dAlph_dt(i,1) = xMinv(i) * RHScos(i)
   dBeta_dt(i,1) = xMinv(i) * RHSsin(i)
END DO !i

 ! Integration using Adams-Bashford predictor corrector method.
 !  Note: The time step is nondimensional. t0=Omega*t
 !  Thus, in YawDyn, DT0=DT*REVS=(DELPSI/REVS)*REVS=DELPSI
 ! USE DT*REVS for compatibility with AeroDyn (DT0 not constant).
 !  DT is in module 'Lift'
 ! Determines xAlpha and xBeta
CALL ABPRECOR ( xAlpha, old_Alph, dAlph_dt, DT0, maxInfl, 1 )
CALL ABPRECOR ( xBeta,  old_Beta, dBeta_dt, DT0, maxInfl, maxInfl0+1 )

 ! Calculate the new lambda_m.
!bjj: why are there 2 do loops? can't they be combined???? (assuming maxInfl0 < maxInfl) or is one of these supposed to be sin?
xLambda_M= 0.
DO k = 1, maxInfl0-1
   xLambda_M = xLambda_M + xLcos(1,k) * xAlpha(k)
END DO !k
DO k = maxInfl0, maxInfl
   xLambda_M = xLambda_M + xLcos(1,k) * xAlpha(k)
END DO !k
! xLambda_M = 2. / sqrt(3.) * xLambda_M
xLambda_M = 1.1547005 * xLambda_M

! Added additional output for GDW debugging - comment out for distribution
!CALL DynDebug (RHScos, RHSsin)



RETURN
END SUBROUTINE infdist

 ! *************************************************************
   SUBROUTINE vindinf( iradius, iblade, rlocal, vnw, VNB, VT, psi )
 !  vindinf calculates the axial induction factor for each
 !   element position using the calculated inflow parameters.
 !  Called by ElemFrc for each element at a new time step.
 ! *************************************************************

USE                           Blade
USE                           DynInflow
USE                           Element
USE                           Switch


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: psi
REAL(ReKi),INTENT(IN)      :: rlocal
REAL(ReKi),INTENT(IN)      :: VNB
REAL(ReKi),INTENT(IN)      :: vnw
REAL(ReKi),INTENT(INOUT)   :: VT

INTEGER,   INTENT(IN)      :: iradius
INTEGER,   INTENT(IN)      :: iblade

   ! Local Variables:

REAL(ReKi)                 :: A2P
!Suzuki uses psiBar, Shawler - WindPsi
!REAL(ReKi)                 :: psibar
REAL(ReKi)                 :: Rzero
REAL(ReKi)                 :: SWRLARG
REAL(ReKi)                 :: Windpsi

INTEGER                    :: mode


 ! Rzero is the non-dimensional radius.
Rzero  = rlocal / r

 ! PSIBAR is the azimuth angle measured counterclockwise from
 !  the horizontal to the left looking downwind.  The directions
 !  of PSI and PSIBAR are opposite and there is 90 deg difference.
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Suzuki:
! psiBar = - psi - piBy2
! Shawler:
CALL WindAzimuthZero (psi,WindPsi)
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ! Calculate the induction factor using the inflow parameters.
 !  Here it's normalized by the tipspeed.

A(iRadius,iBlade) =  0.

DO mode = 1, maxInfl0
   A(iRadius,iBlade) = A(iRadius,iBlade)  &
                     + xphi(Rzero,mode) * xAlpha(mode)
!  &  + phis(Rzero, MRvector(mode), NJvector(mode) )* xAlpha(mode)
END DO !mode

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Suzuki:
!DO mode = maxInfl0+1, maxInfl
!   A(iRadius,iBlade) = A(iRadius,iBlade) + xphi(Rzero,mode) *      &
!!  &    + phis(Rzero, MRvector(mode), NJvector(mode) ) *
!            ( xAlpha(mode) * COS( REAL(MRvector(MODE)) * psibar )  &
!            + xBeta (mode) * SIN( REAL(MRvector(MODE)) * psibar ) )
!END DO !mode
! Shawler:
DO mode = maxInfl0+1, maxInfl
   A(iRadius,iBlade) = A(iRadius,iBlade) + xphi(Rzero,mode) *      &
!  &     + phis(Rzero, MRvector(mode), NJvector(mode) ) *
            ( xAlpha(mode) * COS( REAL(MRvector(MODE)) * Windpsi )  &
            + xBeta (mode) * SIN( REAL(MRvector(MODE)) * Windpsi ) )
END DO !mode
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ! A is positive upwind, while alpha's & beta's are positve downwind.
 !  Also, they are normalized by different factors.
A(iRadius,iBlade) = - A(iRadius,iBlade) * TipSpeed / VNW

 ! Calculate induced swirl (a') if desired.

IF ( SWIRL ) THEN
 ! akihiro 10/26/99
   SWRLARG = 1.0 + 4.0 *  A(iradius,iblade) * VNW * &
!   SWRLARG = 1.0 + 4.0 * TIPLOSS * A(iradius,iblade) * VNW *
           ( (1.0 - A(iradius,iblade)) * VNW + VNB ) / VT / VT
   IF ( SWRLARG > 0.0 ) THEN
      A2P = 0.5 * ( -1.0 + SQRT( SWRLARG ) )
      VT  = VT * ( 1.0 + A2P)
! bjj: this value was not properly set before. We could also just replace the local A2P variable with AP() instead.
      AP(iRadius,iBlade) = A2P
   ENDIF

ENDIF



RETURN
END SUBROUTINE vindinf

 ! ***********************************************************************
   SUBROUTINE ABPRECOR( F, OLDF, DFDT, DT, N, N0 )
 !  Integrates Function F by Adams-Bashford Predictor and Adams-Moulton
 !   Corrector using four previous values of dF/dt.
 ! ***********************************************************************


IMPLICIT                      NONE


   ! Passed Variables:
INTEGER   ,INTENT(IN)      :: N
INTEGER   ,INTENT(IN)      :: N0

REAL(ReKi),INTENT(IN)      :: DT
REAL(ReKi),INTENT(IN)      :: DFDT ( N0:N, 4 )
REAL(ReKi),INTENT(OUT)     :: F    ( N0:N )
REAL(ReKi),INTENT(IN)      :: OLDF ( N0:N )

   ! Local Variables:

REAL(ReKi)                 :: DFDT0

INTEGER                    :: I



DO I = N0, N
 ! Adams-Bashford Predictor
   F(I) = OLDF(I) + ( 55. * DFDT(I,1) - 59. * DFDT(I,2)  &
                  + 37. * DFDT(I,3) -  9. * DFDT(I,4) ) * DT / 24.

 ! New time derivative for corrector
   DFDT0 = ( F(I) - OLDF(I) ) / DT

 ! Adams-Moulton Corrector
   F(I) = OLDF(I) + (  9. * DFDT0     + 19. * DFDT(I,1)  &
                  -  5. * DFDT(I,2) +       DFDT(I,3) ) * DT / 24.

END DO !I



RETURN
END SUBROUTINE ABPRECOR

 ! ***********************************************************************
   SUBROUTINE LMATRIX ( X, matrixMode )
 !  LMATRIX calculates the [L_***] matrix using Gamma matrix and xkai=X.
 !   matrixMode = 1 : Calculate [L_cos] matrix
 !   matrixMode = 2 : Calculate [L_sin] matrix
 ! ***********************************************************************

USE                           DynInflow


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: X

INTEGER   ,INTENT(IN)      :: matrixMode


   ! Local Variables:

REAL(ReKi)                 :: XM

INTEGER                    :: IROW
INTEGER                    :: JCOL


 ! Check the value of X
IF ( ( X < -1.) .OR. ( X > 1.) ) THEN
!   IF ( ( X < 0.) .OR. ( X > 1.) ) THEN
   CALL ProgAbort ( 'Value of X = '//TRIM(Flt2LStr(X))//' must be between -1 and 1.')
ENDIF

SELECT CASE ( matrixMode )

 ! Calculate the rows for r=0 of [L_cos] matrix,
 !  which needs a separate formula.
   CASE (1)
    DO JCOL = 1, maxInfl
       XM = X ** MRvector(JCOL)
       DO IROW = 1, maxInfl0
          xLcos( IROW, JCOL ) = GAMMA( IROW, JCOL ) * XM
       END DO !IROW
    END DO !JCOL

 ! Calculate the [L_cos] matrix,
 !  the rows for r=1 and higher.
    DO IROW = maxInfl0+1, maxInfl
       DO JCOL = 1, maxInfl
          xLcos( IROW, JCOL ) = GAMMA( IROW, JCOL )  &
                  *(    X  ** MminusR( IROW, JCOL )  &
                      + X  ** MplusR ( IROW, JCOL )  &
                  *  (-1.) ** MminR  ( IROW, JCOL )  )
       END DO !JCOL
    END DO !IROW

 ! Calculate [L_sin] matrix.
   CASE (2)
    DO IROW = maxInfl0+1, maxInfl
       DO JCOL = maxInfl0+1, maxInfl
          xLsin( IROW, JCOL ) = GAMMA( IROW, JCOL )  &
                  *(    X  ** MminusR( IROW, JCOL )  &
                      - X  ** MplusR ( IROW, JCOL )  &
                  *  (-1.) ** MminR  ( IROW, JCOL )  )
      END DO !JCOL
    END DO !IROW

   CASE DEFAULT
    CALL ProgAbort( 'Value of matrixMode = '//TRIM(Int2LStr(matrixMode))//' must be 1 or 2.')
END SELECT



RETURN
END SUBROUTINE LMATRIX

 ! ***********************************************************************
   SUBROUTINE MATINV( A0, A1, N, N0, invMode )
 !  Inverts the [L_cos] and [L_sin] matrices.
 !   Subroutine GAUSSJ (modified) from "Numerical Recipe" is needed.
 !   invMode = 1 : Invert [L_cos] matrix
 !   invMode = 2 : Invert [L_sin] matrix
 !**********************************************************************


IMPLICIT                      NONE


   ! Passed Variables:
INTEGER   ,INTENT(IN)      :: invMode
INTEGER   ,INTENT(IN)      :: N
INTEGER   ,INTENT(IN)      :: N0

REAL(ReKi),INTENT(INOUT)   :: A0   (      N   ,      N    )
REAL(ReKi),INTENT(INOUT)   :: A1   ( N0+1:N   , N0+1:N    )

   ! Local Variables:

REAL(ReKi)                 :: DUMMY(      N-N0,      N-N0 )

INTEGER                    :: I
INTEGER                    :: J



 ! Invert [L_cos] in all cases
 !  [L_cos] matrix can be inverted without a dummy.
 !  Invert the [A0]=[L_cos] matrix by Gauss-Jordan Method
SELECT CASE ( invMode )

   CASE (1)
    CALL GAUSSJ(A0,N)

 ! [L_sin] matrix needs a dummy array, because the index goes
 !  from maxInfl0(=N0) to maxInfl(=N),
 !  which is incompatible with SUBROUTINE GAUSSJ.
 !BJJ: IS THIS REALLY NECESSARY?  Aren't the indicies an abstraction?? if you pass an array (1:3) and use it as (0:2) in another subroutine, it's really okay???
   CASE (2)
    DO I=1,N-N0
       DO J=1,N-N0
          DUMMY(I,J) = A1(I+N0,J+N0)
       END DO !J
    END DO !I

 ! Invert the [A1]=[L_sin] matrix by Gauss-Jordan Method.
    CALL GAUSSJ(DUMMY,N-N0)

 ! Put the dummy back into [A1]=[L_sin].
    DO I=1,N-N0
       DO J=1,N-N0
          A1(I+N0,J+N0) = DUMMY(I,J)
       END DO !J
    END DO !I

   CASE DEFAULT
   CALL ProgAbort( 'Value of invMode = '//TRIM(Int2LStr(invMode))//' must be 1 or 2.')

END SELECT



RETURN
END SUBROUTINE MATINV

 ! ***********************************************************************
   SUBROUTINE gaussj(a,n)
 !  Invert a matrix by Gauss-Jordan Method. The original source code
 !   from "Numerical Recipe" was slightly modified.
 ! ***********************************************************************


IMPLICIT                      NONE


   ! Passed Variables:
INTEGER   ,INTENT(IN)      :: n

REAL(ReKi),INTENT(INOUT)   :: a(n,n)

   ! Local Variables:

INTEGER   , PARAMETER      :: NMAX = 6

REAL(ReKi)                 :: big
REAL(ReKi)                 :: dum
REAL(ReKi)                 :: pivinv

INTEGER                    :: i
INTEGER                    :: icol
INTEGER                    :: indxc(NMAX)
INTEGER                    :: indxr(NMAX)
INTEGER                    :: ipiv(NMAX)
INTEGER                    :: irow
INTEGER                    :: j
INTEGER                    :: k
INTEGER                    :: l
INTEGER                    :: ll



DO j=1,n
   ipiv(j)=0
END DO !j
DO i=1,n
   big=0.
   DO j=1,n
      IF (ipiv(j) /= 1) THEN
         DO k=1,n
            IF (ipiv(k) == 0) THEN
               IF (ABS(a(j,k)) >= big) THEN
                  big=ABS(a(j,k))
                  irow=j
                  icol=k
                ENDIF

            ELSE IF (ipiv(k) > 1) THEN
               CALL ProgAbort( 'Singular matrix encountered.' )
            ENDIF
         END DO !k
      ENDIF
   END DO !j
   ipiv(icol)=ipiv(icol)+1
   IF (irow /= icol) THEN
      DO l=1,n
         dum=a(irow,l)
         a(irow,l)=a(icol,l)
         a(icol,l)=dum
      END DO !l
   ENDIF
   indxr(i)=irow
   indxc(i)=icol
   IF (a(icol,icol) == 0.) THEN
      CALL ProgAbort( 'Singular matrix encountered.' )
   ENDIF
   pivinv=1./a(icol,icol)
   a(icol,icol)=1.
   DO l=1,n
      a(icol,l)=a(icol,l)*pivinv
   END DO !l
   DO ll=1,n
      if (ll /= icol) THEN
         dum=a(ll,icol)
         a(ll,icol)=0.
         DO l=1,n
            a(ll,l)=a(ll,l)-a(icol,l)*dum
         END DO !l
      ENDIF
   END DO !ll
END DO !i
DO l=n,1,-1
   IF (indxr(l) /= indxc(l)) THEN

      DO k=1,n
         dum=a(k,indxr(l))
         a(k,indxr(l))=a(k,indxc(l))
         a(k,indxc(l))=dum
      END DO !k
   ENDIF
END DO !l



RETURN
END SUBROUTINE gaussj

! ***********************************************************************
  FUNCTION FGAMMA( R, J, M, N )
!  Calculate the GAMMA matrix. It is NOT the statistical function.
! ***********************************************************************


IMPLICIT                      NONE


   ! Passed Variables:

REAL(ReKi)                  :: FGAMMA
INTEGER   ,INTENT(IN)       :: J
INTEGER   ,INTENT(IN)       :: M
INTEGER   ,INTENT(IN)       :: N
INTEGER   ,INTENT(IN)       :: R


IF ( MOD(R+M,2) == 0 ) THEN
   FGAMMA = (-1)**((N+J-2*R)*.5) * 2.        &
          * SQRT( REAL( (2*N+1) * (2*J+1) ) ) &
          / SQRT( HFUNC(M,N) * HFUNC(R,J) )   &
          / REAL( (J+N) * (J+N+2) * ((J-N)*(J-N)-1) )

ELSE IF ( ABS(J-N) == 1 ) THEN  !bjj: why don't we use the pi() variable?
   FGAMMA = 3.14159265 * SIGN(1., REAL(R-M) ) * .5 &
          / SQRT( HFUNC(M,N) * HFUNC(R,J) )        &
          / SQRT( REAL( (2*N+1) * (2*J+1) ) )

ELSE
   FGAMMA = 0.

ENDIF



RETURN
END FUNCTION FGAMMA

 ! ***********************************************************************
   FUNCTION HFUNC( M, N )
 !  Calculates the value of function H(m,n).
 !   Warning: This subroutine is not optimized for large m or n.
 !   Although H(m,n) is a well behaving function, it may
 !   cause math overflow, if m or n is large.
 !bjj: we only call with with the MRvector and NJvector (parameter) values.  This could
 !  possibly increase performance if implemented differently....
 ! ***********************************************************************


IMPLICIT                      NONE


   ! Passed Variables:

REAL(ReKi)                 :: HFUNC

INTEGER   ,INTENT(IN)      :: M
INTEGER   ,INTENT(IN)      :: N


   ! Local Variables:

INTEGER                    :: NMM   ! n minus M
INTEGER                    :: NPM   ! N plus M



IF ( N <= M ) THEN
   CALL ProgAbort( 'Value of N = '//TRIM(Int2LStr(N))//' must be geater than M = '//TRIM(Int2LStr(M))//'.' )
ENDIF

NPM = N + M
NMM = N - M

HFUNC = ( REAL( IDUBFACT(NPM-1) ) / REAL( IDUBFACT(NPM) ) ) &
      * ( REAL( IDUBFACT(NMM-1) ) / REAL( IDUBFACT(NMM) ) )



RETURN
END FUNCTION HFUNC

 ! ***********************************************************************
   INTEGER FUNCTION IDUBFACT( I )
 !  Calculates the double factorial of an integer I
 !   IDUBFACT( I ) = I!! = I*(I-2)*(I-4)*...*4*2 for I = even
 !                    or = I*(I-2)*(I-4)*...*3*1 for I = odd
 ! ***********************************************************************


IMPLICIT                   NONE


   ! Passed Variables:
INTEGER   ,INTENT(IN)   :: I

   ! Local Variables:

INTEGER                 :: K



IF ( I >= 1 ) THEN
   IDUBFACT = 1

   DO K = I, 1, -2
      IDUBFACT = IDUBFACT * K
   END DO !K

ELSE IF ( I == 0 .OR. I == -1 ) THEN
   IDUBFACT = 1
ELSE IF ( I == -3 ) THEN  ! use definition of n!! for odd negative numbers
   IDUBFACT = -1
ELSE
   CALL ProgAbort( 'Double factorial is NOT defined for '//TRIM(Int2LStr(I))//' in FUNCTION IDUBFACT.')
ENDIF



RETURN
END FUNCTION IDUBFACT

 ! ***********************************************************************
   FUNCTION xphi( Rzero, mode )
 !  Set up the PHI coefficients. They are the results from Mathematica.
 !  phi(1)  = sqrt(   3.)                    ! m=0, n=1
 !  phi(2)  = 2*sqrt(7) (1.5 - 3.75 Rzero**2 ) / 3.! m=0, n=3
 !  phi(3)  = sqrt(  15./ 2.) *Rzero            ! m=1, n=2
 !  phi(4)  = 4*(15/4 * Rzero - 105/16 *Rzero**3 )/sqrt(5)! m=1, n=4
 !  phi(5)  = sqrt( 105./ 2.) / 2.  *Rzero**2      ! m=2, n=3
 !  phi(6)  = sqrt(  35.) *3. / 4.  *Rzero**3      ! m=3, n=4
 ! ***********************************************************************


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: Rzero
REAL(ReKi)                 :: xphi

INTEGER   ,INTENT(IN)      :: mode




IF ( Rzero < 0. ) THEN
   CALL ProgAbort( 'Value of Rzero = '//TRIM(Flt2LStr(Rzero))//' must be larger than 0 in xphi().')
ELSE IF ( Rzero > 1. ) THEN
   CALL ProgAbort( 'Value of Rzero = '//TRIM(Flt2LStr(Rzero))//' must be smaller than 1 in xphi().')
ENDIF

SELECT CASE ( mode )
   CASE (1)
    xphi =   1.732051
   CASE (2)
    xphi =   2.645751 - 6.6143783 * Rzero * Rzero
   CASE (3)
    xphi =   2.738613 * Rzero
   CASE (4)
    xphi = ( 6.708204 - 11.73936 * Rzero * Rzero ) * Rzero
   CASE (5)
    xphi =   3.622844 * Rzero * Rzero
   CASE (6)
    xphi =   4.437060 * Rzero * Rzero * Rzero
   CASE DEFAULT
      CALL ProgAbort('Integer MODE = '//TRIM(Int2LStr(MODE))//' must be 1 through 6.' )
END SELECT



RETURN
END FUNCTION xphi

 !***********************************************************************
 ! akihiro 06/25/00
 ! This subroutine is not currently used, may be used in the future.
 !
   FUNCTION phis( Rzero, r, j )
 !  Calculates the PHI coefficients. This function is not used unless
 !   the # of inflow states is increased greater than 10 (Mode 7 and higher)
 ! ***********************************************************************


IMPLICIT                     NONE


   ! Passed Variables:

REAL(ReKi)                :: phis
REAL(ReKi),INTENT(IN)     :: Rzero

INTEGER   ,INTENT(IN)     :: j
INTEGER   ,INTENT(IN)     :: r

   ! Local Variables:

INTEGER                   :: q


IF ( Rzero < 0. ) THEN
   CALL ProgAbort('Value of Rzero = '//TRIM(Flt2LStr(Rzero))//' must be larger than 0 in phis().' )
ELSE IF ( Rzero > 1. ) THEN
   CALL ProgAbort('Value of Rzero = '//TRIM(Flt2LStr(Rzero))//' must be smaller than 1 in phis().' )
ENDIF

phis = 0.

DO q = r, j-1, 2
   phis = phis  &
        + Rzero ** q * (-1.) **((q-r)/2) * REAL( idubfact(j+q) ) &
        / REAL( idubfact(q-r) * idubfact(q+r) * idubfact(j-q-1) )
END DO !q

phis = phis * SQRT( REAL( 2*j+1 ) * hfunc(r,j) )


RETURN
END FUNCTION phis


!***********************************************************
SUBROUTINE WindAzimuthZero (psi,WindPsi)
! Subroutine added by JRS to define the zero azimuth datum in
! a wind based co-ordinate system, for USE in the dynamic inflow
! routines.
! Calculates the rotational measurement in radians
! of the resultant of two vectors, VrotorZ and VrotorY.
! VrotorZ is positive vertically upwards and negative vertically
! downwards, VrotorY is positive to the left and negative
! to the right, both when looking towards the rotor from upwind.
! Zero degrees azimuth is defined when veritically down
! and rises with a clockwise rotation.


USE                           Wind


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: psi
REAL(ReKi),INTENT(OUT)     :: WindPsi

WindPsi = psi - ATAN2(VrotorY,-VrotorZ)

RETURN
END SUBROUTINE WindAzimuthZero

 !     ************* END OF FILE ***************

! ********************************************************************
  SUBROUTINE AeroDyn_Terminate ( )
! ********************************************************************

   USE   AD_IOParams
   USE   Airfoil
   USE   Bedoes
   USE   Blade
   USE   DynInflow
   USE   Element
   USE   ElOutParams
   USE   ElemInflow
   USE   TwrProps


      ! Deallocate memory for arrays stored in AeroDyn modules

        ! Airfoil
   IF ( ALLOCATED( FOILNM     ) )   DEALLOCATE( FOILNM      )
   IF ( ALLOCATED( AL         ) )   DEALLOCATE( AL          )
   IF ( ALLOCATED( CD         ) )   DEALLOCATE( CD          )
   IF ( ALLOCATED( CL         ) )   DEALLOCATE( CL          )
   IF ( ALLOCATED( CM         ) )   DEALLOCATE( CM          )
   IF ( ALLOCATED( MulTabMet  ) )   DEALLOCATE( MulTabMet   )
   IF ( ALLOCATED( NFOIL      ) )   DEALLOCATE( NFOIL       )
   IF ( ALLOCATED( NLIFT      ) )   DEALLOCATE( NLIFT       )
   IF ( ALLOCATED( NTables    ) )   DEALLOCATE( NTables     )
   IF ( ALLOCATED( FOILNM     ) )   DEALLOCATE( FOILNM      )


        ! Bedoes
   IF ( ALLOCATED( ADOT       ) )   DEALLOCATE( ADOT        )
   IF ( ALLOCATED( ADOT1      ) )   DEALLOCATE( ADOT1       )
   IF ( ALLOCATED( AFE        ) )   DEALLOCATE( AFE         )
   IF ( ALLOCATED( AFE1       ) )   DEALLOCATE( AFE1        )
   IF ( ALLOCATED( ANE        ) )   DEALLOCATE( ANE         )
   IF ( ALLOCATED( ANE1       ) )   DEALLOCATE( ANE1        )
   IF ( ALLOCATED( AOD        ) )   DEALLOCATE( AOD         )
   IF ( ALLOCATED( AOL        ) )   DEALLOCATE( AOL         )
   IF ( ALLOCATED( CDO        ) )   DEALLOCATE( CDO         )
   IF ( ALLOCATED( CNA        ) )   DEALLOCATE( CNA         )
   IF ( ALLOCATED( CNP        ) )   DEALLOCATE( CNP         )
   IF ( ALLOCATED( CNP1       ) )   DEALLOCATE( CNP1        )
   IF ( ALLOCATED( CNPD       ) )   DEALLOCATE( CNPD        )
   IF ( ALLOCATED( CNPD1      ) )   DEALLOCATE( CNPD1       )
   IF ( ALLOCATED( CNPOT      ) )   DEALLOCATE( CNPOT       )
   IF ( ALLOCATED( CNPOT1     ) )   DEALLOCATE( CNPOT1      )
   IF ( ALLOCATED( CNS        ) )   DEALLOCATE( CNS         )
   IF ( ALLOCATED( CNSL       ) )   DEALLOCATE( CNSL        )
   IF ( ALLOCATED( CNV        ) )   DEALLOCATE( CNV         )
   IF ( ALLOCATED( CVN        ) )   DEALLOCATE( CVN         )
   IF ( ALLOCATED( CVN1       ) )   DEALLOCATE( CVN1        )
   IF ( ALLOCATED( DF         ) )   DEALLOCATE( DF          )
   IF ( ALLOCATED( DFAFE      ) )   DEALLOCATE( DFAFE       )
   IF ( ALLOCATED( DFAFE1     ) )   DEALLOCATE( DFAFE1      )
   IF ( ALLOCATED( DFC        ) )   DEALLOCATE( DFC         )
   IF ( ALLOCATED( DN         ) )   DEALLOCATE( DN          )
   IF ( ALLOCATED( DPP        ) )   DEALLOCATE( DPP         )
   IF ( ALLOCATED( DQ         ) )   DEALLOCATE( DQ          )
   IF ( ALLOCATED( DQP        ) )   DEALLOCATE( DQP         )
   IF ( ALLOCATED( DQP1       ) )   DEALLOCATE( DQP1        )
   IF ( ALLOCATED( FSP        ) )   DEALLOCATE( FSP         )
   IF ( ALLOCATED( FSP1       ) )   DEALLOCATE( FSP1        )
   IF ( ALLOCATED( FSPC       ) )   DEALLOCATE( FSPC        )
   IF ( ALLOCATED( FSPC1      ) )   DEALLOCATE( FSPC1       )
   IF ( ALLOCATED( FTB        ) )   DEALLOCATE( FTB         )
   IF ( ALLOCATED( FTBC       ) )   DEALLOCATE( FTBC        )
   IF ( ALLOCATED( OLDCNV     ) )   DEALLOCATE( OLDCNV      )
   IF ( ALLOCATED( OLDDF      ) )   DEALLOCATE( OLDDF       )
   IF ( ALLOCATED( OLDDFC     ) )   DEALLOCATE( OLDDFC      )
   IF ( ALLOCATED( OLDDN      ) )   DEALLOCATE( OLDDN       )
   IF ( ALLOCATED( OLDDPP     ) )   DEALLOCATE( OLDDPP      )
   IF ( ALLOCATED( OLDDQ      ) )   DEALLOCATE( OLDDQ       )
   IF ( ALLOCATED( OLDTAU     ) )   DEALLOCATE( OLDTAU      )
   IF ( ALLOCATED( OLDXN      ) )   DEALLOCATE( OLDXN       )
   IF ( ALLOCATED( OLDYN      ) )   DEALLOCATE( OLDYN       )
   IF ( ALLOCATED( QX         ) )   DEALLOCATE( QX          )
   IF ( ALLOCATED( QX1        ) )   DEALLOCATE( QX1         )
   IF ( ALLOCATED( TAU        ) )   DEALLOCATE( TAU         )
   IF ( ALLOCATED( XN         ) )   DEALLOCATE( XN          )
   IF ( ALLOCATED( YN         ) )   DEALLOCATE( YN          )
   IF ( ALLOCATED( BEDSEP     ) )   DEALLOCATE( BEDSEP      )
   IF ( ALLOCATED( OLDSEP     ) )   DEALLOCATE( OLDSEP      )


        ! Blade
   IF ( ALLOCATED( C          ) )   DEALLOCATE( C           )
   IF ( ALLOCATED( DR         ) )   DEALLOCATE( DR          )


        ! DynInflow
   IF ( ALLOCATED( RMC_SAVE   ) )   DEALLOCATE( RMC_SAVE    )
   IF ( ALLOCATED( RMS_SAVE   ) )   DEALLOCATE( RMS_SAVE    )


        ! Element
   IF ( ALLOCATED( A          ) )   DEALLOCATE( A           )
   IF ( ALLOCATED( AP         ) )   DEALLOCATE( AP          )
   IF ( ALLOCATED( HLCNST     ) )   DEALLOCATE( HLCNST      )
   IF ( ALLOCATED( RELM       ) )   DEALLOCATE( RELM        )
   IF ( ALLOCATED( TLCNST     ) )   DEALLOCATE( TLCNST      )
   IF ( ALLOCATED( TWIST      ) )   DEALLOCATE( TWIST       )


        ! ElOutParams
   IF ( ALLOCATED( AAA        ) )   DEALLOCATE( AAA         )
   IF ( ALLOCATED( AAP        ) )   DEALLOCATE( AAP         )
   IF ( ALLOCATED( ALF        ) )   DEALLOCATE( ALF         )
   IF ( ALLOCATED( CDD        ) )   DEALLOCATE( CDD         )
   IF ( ALLOCATED( CLL        ) )   DEALLOCATE( CLL         )
   IF ( ALLOCATED( CMM        ) )   DEALLOCATE( CMM         )
   IF ( ALLOCATED( CNN        ) )   DEALLOCATE( CNN         )
   IF ( ALLOCATED( CTT        ) )   DEALLOCATE( CTT         )
   IF ( ALLOCATED( DFNSAV     ) )   DEALLOCATE( DFNSAV      )
   IF ( ALLOCATED( DFTSAV     ) )   DEALLOCATE( DFTSAV      )
   IF ( ALLOCATED( DynPres    ) )   DEALLOCATE( DynPres     )
   IF ( ALLOCATED( PITSAV     ) )   DEALLOCATE( PITSAV      )
   IF ( ALLOCATED( PMM        ) )   DEALLOCATE( PMM         )
   IF ( ALLOCATED( ReyNum     ) )   DEALLOCATE( ReyNum      )
   IF ( ALLOCATED( SaveVX     ) )   DEALLOCATE( SaveVX      )
   IF ( ALLOCATED( SaveVY     ) )   DEALLOCATE( SaveVY      )
   IF ( ALLOCATED( SaveVZ     ) )   DEALLOCATE( SaveVZ      )
   IF ( ALLOCATED( WndElPrList) )   DEALLOCATE( WndElPrList )
   IF ( ALLOCATED( WndElPrNum ) )   DEALLOCATE( WndElPrNum  )
   IF ( ALLOCATED( ElPrList   ) )   DEALLOCATE( ElPrList    )
   IF ( ALLOCATED( ElPrNum    ) )   DEALLOCATE( ElPrNum     )

      ! ElemInflow

   IF ( ALLOCATED( W2         ) )   DEALLOCATE( W2          )
   IF ( ALLOCATED( Alpha      ) )   DEALLOCATE( Alpha       )

         ! TwrProps
   IF ( ALLOCATED( TwrHtFr ) )      DEALLOCATE ( TwrHtFr    )
   IF ( ALLOCATED( TwrWid  ) )      DEALLOCATE ( TwrWid     )
   IF ( ALLOCATED( TwrCD  ) )       DEALLOCATE ( TwrCD      )
   IF ( ALLOCATED( TwrRe ) )        DEALLOCATE ( TwrRe      )
   IF ( ALLOCATED( NTwrCDCol ) )    DEALLOCATE ( NTwrCDCol  )

   ! DEALLOCATE ( OLD_A_NS )     !FH these are local variables of SUBROUTINE VIND
   ! DEALLOCATE ( OLD_AP_NS )       !FH which cannot be DEALLOCATED here ... still a small memory leak ...

      ! Close files that were opened in AeroDyn

      ! AD_IOParams
   CLOSE(UnADin)              ! ipt file
   CLOSE(UnADopt)             ! opt file
   CLOSE(UnAirfl)             ! Airfoil data file
   CLOSE(UnWind)              ! HH or FF wind file


      ! NWTC_Library
   CLOSE(UnEc)                ! I/O unit number for the echo file.


!      ! Close wind inflow
!   CALL WindInf_Terminate( ErrStat )

   RETURN
END SUBROUTINE AeroDyn_Terminate
!====================================================================================================
SUBROUTINE CheckRComp( ADFile, HubRadius, TipRadius, ErrStat )
! This routine checks to see if RElm(:) and DR(:) are compatible within a millimeter;
!----------------------------------------------------------------------------------------------------

   USE                           Blade,      ONLY: DR
   USE                           Element,    ONLY: NElm, RElm

   IMPLICIT                      NONE

      ! Passed variables

   CHARACTER(*), INTENT(IN)   :: ADFile                     ! Name of the AeroDyn input file, used for printing error message
   REAL(ReKi),   INTENT(IN)   :: HubRadius                  ! Hub radius, used to verify that RElm and DR are input correctly
   REAL(ReKi),   INTENT(IN)   :: TipRadius                  ! Tip radius, used to verify that RElm and DR are input correctly
   INTEGER,      INTENT(OUT)  :: ErrStat


      ! Local variables.

   REAL(ReKi)                 :: DRNodesNew(NElm)           ! Length of variable-spaced blade elements--calculated from input RElm(:).
   REAL(ReKi)                 :: DRSum                      ! Sum of DRs--should be close to TipRadius
   REAL(ReKi), PARAMETER      :: EPS = EPSILON(HubRadius)   ! A small value used to compare two numbers

   INTEGER                    :: I                          ! Generic index.

   CHARACTER(33)              :: DRChange                   ! A string showing how to change DR to get campatibility



      ! Initialize ErrStat to 0 (no error):

   ErrStat = 0


      ! Calculate DRNodesNew(:) based on input RElm(:) and compare to input DR(:):

!   AssumedHubRadius = RElm(1) - 0.5*DR(1)
   DRNodesNew(1)    = 2.0*( RElm(1) - HubRadius )
   DRSum            = DRNodesNew(1) + HubRadius

   IF ( DRNodesNew(1) <= EPS )  THEN                              ! Check to see if RElm(1) > HubRad; if not, ProgAbort program
      CALL WrScr(' RElm(1) must be larger than the hub radius (HubRadius = RElm(1) - 0.5*DR(1)). ')
      ErrStat = 1
      RETURN
   ELSEIF ( ABS( DRNodesNew(1) - DR(1) ) > 0.001 )  THEN     ! Check to see if the calculated DRNodes(1) is close to the inputted DRNodes(1); if not, set flag--this will cause the program to ProgAbort later.
      ErrStat = 1
   ENDIF

   DO I = 2,NElm ! Loop through all but the innermost blade element

      DRNodesNew(I) = 2.0*( RElm(I) - RElm(I-1) ) - DRNodesNew(I-1)
      DRSum         = DRSum + DRNodesNew(I)

      IF ( DRNodesNew(I) <= EPS )  THEN                           ! Check to see if it is possible to have compatible DR(:) with the input RElm(:); if not, abort program
         CALL WrScr( 'RElm('//TRIM( Int2LStr(I) )//') produces ill-conditioned DR(:)' )
         ErrStat = 1
         RETURN
      ELSEIF ( ABS( DRNodesNew(I) - DR(I) ) > 0.001 )  THEN  ! Check to see if the calculated DRNodes(I) is close to the inputted DRNodes(I); if not, set flag--this will cause the program to Abort later.
         ErrStat = 1
      ENDIF

   END DO             ! I - all but the innermost blade element


      ! Abort program if necessary

   IF ( ErrStat /= 0 )  THEN

         ! Write error message since the input DR(:) are not close to the calculated DRNodesNew(:)

      CALL WrScr1(' Input values for DR(:) are not compatible with input RElm(:). ')
      CALL WrScr (' To make them compatible, please modify DR in the AeroDyn input file, '// TRIM( ADFile ) //', as follows:')
      CALL WrScr1(' DR (Old) --> DR (New) ')

      DO I = 1,NElm
         WRITE( DRChange, "(' ', F13.4, ' --> ', F13.4, ' ')" ) DR(I), DRNodesNew(I)
         CALL WrScr( DRChange )
      ENDDO !I

      RETURN

   ELSEIF ( ABS( DRSum - TipRadius ) > 0.001 )  THEN

         ! Abort program since SUM( DRNodes(:) ) /= ( TipRadius - HubRadius)

      CALL WrScr(' TipRadius must be equal to HubRadius + SUM( DR(:) ). ')
      ErrStat = 1

      RETURN

   ENDIF


   RETURN
END SUBROUTINE CheckRComp


END MODULE AeroSubs
