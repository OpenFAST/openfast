!**********************************************************************************************************************************
MODULE AeroSubs

   USE                     NWTC_Library
   USE                     AeroDyn14_Types

   USE AeroGenSubs, ONLY : MaxInfl


   IMPLICIT                NONE

    
   INTEGER(IntKi) , PARAMETER, DIMENSION(1:6)  :: MRvector = (/ 0, 0, 1, 1, 2, 3 /)  !bjj why aren't these parameters? Now they are.
   INTEGER(IntKi) , PARAMETER, DIMENSION(1:6)  :: NJVector = (/ 1, 3, 2, 4, 3, 4 /)    
                     

CONTAINS
! ************************************************
!  AeroDyn Subroutines for YawDyn, ADAMS,
!   SymDyn and FAST.
! ************************************************
!  UNIVERSITY OF UTAH, MECHANICAL ENGINEERING DEPARTMENT


! Updated version that uses FAST Interface types
!====================================================================================================
SUBROUTINE AD14_GetInput(InitInp, P, x, xd, z, m, y, ErrStat, ErrMess )
!====================================================================================================
   USE                           AeroGenSubs, ONLY: AllocArrays
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_InitInputType),       INTENT(INOUT)  :: InitInp
   TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD14_ContinuousStateType), INTENT(INOUT)  :: x           ! Initial continuous states
   TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(AD14_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Initial misc/optimization variables
   TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMess

      ! Local Variables:

   INTEGER                    :: ElIndex
   INTEGER                    :: IELM
   INTEGER                    :: IndPrint
   INTEGER                    :: K
   INTEGER                    :: UnIn
   INTEGER                    :: ErrStatLcl

   LOGICAL                    :: PremEOF_indicator
   CHARACTER(1024)            :: LINE
   CHARACTER(1024)            :: FilePath             ! The path name of the AeroDyn input file (so files listed in it can be defined relative to the main input file location)
   CHARACTER(ErrMsgLen)       :: ErrMessLcl
   character(*), parameter    :: RoutineName = 'AD14_GetInput'

   !bjj: error handling here needs to be fixed! (we overwrite any non-AbortErrLev errors)
   
   ErrStat = ErrID_None
   ErrMess = ''
   
      ! Function definition

   call GetNewUnit(UnIn, ErrStatLcl, ErrMessLcl)

   !-------------------------------------------------------------------------------------------------
   ! Open the AeroDyn input file
   !-------------------------------------------------------------------------------------------------
   CALL OpenFInpFile(UnIn, TRIM(InitInp%ADFileName), ErrStatLcl, ErrMessLcl)
   CALL SetErrStat(ErrStatLcl, ErrMessLcl,ErrStat, ErrMess,RoutineName)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(UnIn)
      RETURN
   END IF

   CALL GetPath( InitInp%ADFileName, FilePath )

   !-------------------------------------------------------------------------------------------------
   ! If the echo file is open, write the header...
   !-------------------------------------------------------------------------------------------------
   IF ( p%Echo ) THEN
      WRITE( p%UnEc, '(// A /)' ) 'AeroDyn input data from file "'//TRIM( InitInp%ADFileName )//'":'
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Read the AeroDyn input file
   !-------------------------------------------------------------------------------------------------

   CALL ReadCom( UnIn, InitInp%ADFileName, 'Header', ErrStatLcl, ErrMessLcl )
      CALL SetErrStat(ErrStatLcl, ErrMessLcl,ErrStat, ErrMess,RoutineName)   
   
      ! Read in the title line
   CALL ReadStr( UnIn, InitInp%ADFileName, InitInp%Title, VarName='Title', VarDescr='File title', ErrStat=ErrStatLcl, ErrMsg=ErrMessLcl)
   CALL SetErrStat(ErrStatLcl, ErrMessLcl,ErrStat, ErrMess,RoutineName)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(UnIn)
      RETURN
   END IF

   IF (NWTC_VerboseLevel == NWTC_Verbose) THEN
      CALL WrScr( ' Heading of the AeroDyn input file: '//NewLine//'   '//TRIM(InitInp%Title) )
   END IF
   p%TITLE = InitInp%TITLE

   p%SIunit = .TRUE.

      ! Read in the stall model
   CALL ReadVar( UnIn, InitInp%ADFileName, LINE, VarName='DStall', VarDescr='Stall model', ErrStat=ErrStat, ErrMsg=ErrMess)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(UnIn)
      RETURN
   END IF

   CALL Conv2UC(LINE(1:7))

   SELECT CASE ( TRIM(Line) )
      CASE ('STEADY')
         P%DStall = .FALSE.
      CASE ('BEDDOES')   !  added -- maybe the input format changed? jm
         P%DStall = .TRUE.
      CASE DEFAULT
         CALL ProgWarn( ' Error: Expecting "STEADY" or "BEDDOES" stall model option.')
         ErrStat = ErrID_Fatal
         CLOSE(UnIn)
         RETURN
   END SELECT


      ! Read in the CM option
   CALL ReadVar( UnIn, InitInp%ADFileName, LINE, VarName='PMoment', VarDescr='Pitching moment option', ErrStat=ErrStat, ErrMsg=ErrMess)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(UnIn)
      RETURN
   END IF

   CALL Conv2UC(LINE(1:6))

   SELECT CASE ( TRIM(Line) )
      CASE ('USE_CM')
         P%PMoment = .TRUE.
      CASE ('NO_CM')
         P%PMoment = .FALSE.
      CASE DEFAULT
         CALL ProgWarn( ' Error: Expecting "USE_CM" or "NO_CM" pitching moment option.')
         ErrStat = ErrID_Fatal
         CLOSE(UnIn)
         RETURN
   END SELECT


      ! Read in the inflow model option
   CALL ReadVar( UnIn, InitInp%ADFileName, LINE, VarName='DynInfl', VarDescr='Inflow model', ErrStat=ErrStat, ErrMsg=ErrMess)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(UnIn)
      RETURN
   END IF
   
   CALL Conv2UC(LINE(1:7))

   SELECT CASE ( Line(1:5) )
      CASE ('EQUIL')
         P%DynInfl = .FALSE.
         m%DynInit = .FALSE.

         IF (Line(6:7) == 'DA') THEN
            P%Inducedvel%EqAIDmult = 1.0
            P%Inducedvel%EquilDA   = .TRUE.
            P%Inducedvel%EquilDT   = .FALSE.
         ELSEIF (LINE(6:7) == 'DB') THEN
            P%Inducedvel%EqAIDmult = 1.0
            P%Inducedvel%EquilDA   = .TRUE.
            P%Inducedvel%EquilDT   = .TRUE.
         ELSEIF (LINE(6:7) == 'DT') THEN
            P%Inducedvel%EqAIDmult = 0.0
            P%Inducedvel%EquilDA   = .FALSE.
            P%Inducedvel%EquilDT   = .TRUE.
         ELSE
            P%Inducedvel%EqAIDmult = 0.0
            P%Inducedvel%EquilDA   = .FALSE.
            P%Inducedvel%EquilDT   = .FALSE.
         ENDIF

      CASE ('DYNIN')
         P%DynInfl = .TRUE.
         m%DynInit = .TRUE.
      CASE DEFAULT
         CALL ProgWarn( ' Error: Expecting "EQUIL" or "DYNIN" inflow model option.')
         ErrStat = ErrID_Fatal
         CLOSE(UnIn)
         RETURN
   END SELECT


      ! Read in the wake model
   CALL ReadVar( UnIn, InitInp%ADFileName, LINE, VarName='Wake', VarDescr='Wake model', ErrStat=ErrStat, ErrMsg=ErrMess)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(UnIn)
      RETURN
   END IF
   CALL Conv2UC(LINE(1:5))

   SELECT CASE ( TRIM(Line) )
      CASE ('NONE')
         P%Wake = .FALSE.
         P%Swirl = .FALSE.

         CALL ProgWarn( ' All wake calculations are turned off! This option is recommended only '// &
                        'in high winds or for debugging.' )
      CASE ('WAKE')
         P%Wake = .TRUE.
         P%Swirl = .FALSE.
      CASE ('SWIRL')
         P%Wake = .TRUE.
         P%Swirl = .TRUE.
      CASE DEFAULT
         CALL ProgWarn( ' Error: Expecting "NONE", "WAKE", or "SWIRL" wake model option.')
         ErrStat = ErrID_Fatal
         CLOSE(UnIn)
         RETURN
   END SELECT

      ! Read in the tolerance for the wake model
   CALL ReadVar( UnIn, InitInp%ADFileName, P%Inducedvel%AToler, VarName='AToler', VarDescr='Induction factor tolerance', ErrStat=ErrStat, ErrMsg=ErrMess)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(UnIn)
      RETURN
   END IF

       ! Read in the tip-loss model for EQUIL inflow
   CALL ReadVar( UnIn, InitInp%ADFileName, LINE, VarName='TLoss', VarDescr='Tip-loss model', ErrStat=ErrStat, ErrMsg=ErrMess)
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(UnIn)
      RETURN
   END IF
   
   IF ( P%DynInfl ) THEN          ! Initialize these variables, though they shouldn't be used
      P%Inducedvel%TLoss = .FALSE.
      P%Inducedvel%GTech = .FALSE.
   ELSE
      CALL Conv2UC(LINE(1:5))

      SELECT CASE ( LINE(1:5) )
         CASE ('NONE ')
            P%Inducedvel%TLoss = .FALSE.
            P%Inducedvel%GTech = .FALSE.
         CASE ('PRAND')
            P%Inducedvel%TLoss = .TRUE.
            P%Inducedvel%GTech = .FALSE.
         CASE ('GTECH')
            P%Inducedvel%TLoss = .TRUE.
            P%Inducedvel%GTech = .TRUE.
         CASE DEFAULT
            CALL ProgWarn( ' Error: Expecting "NONE", "PRAND", or "GTECH" tip-loss model option.')
            ErrStat = ErrID_Fatal
            CLOSE(UnIn)
            RETURN
      END SELECT
   END IF


       ! Read in the hub-loss model for EQUIL inflow
   CALL ReadVar( UnIn, InitInp%ADFileName, LINE, VarName='HLoss', VarDescr='Hub-loss model', ErrStat=ErrStat, ErrMsg=ErrMess)
   IF (ErrStat >= AbortErrLev) RETURN

   IF ( P%DynInfl ) THEN          ! Initialize these variables, though they shouldn't be used
      P%Inducedvel%HLoss = .FALSE.
   ELSE
      CALL Conv2UC( LINE(1:5) )

      SELECT CASE ( LINE(1:5) )
         CASE ('NONE')
            P%Inducedvel%HLoss = .FALSE.
         CASE ('PRAND')
            P%Inducedvel%HLoss = .TRUE.
         CASE DEFAULT
            CALL ProgWarn( ' Error: Expecting "NONE" or "PRAND" hub-loss model option.')
            ErrStat = ErrID_Fatal
            CLOSE(UnIn)
            RETURN
      END SELECT
   END IF

   p%Rotor%HH = InitInp%HubHt   
   
!bjj: this is a hack job to allow both the new tower influence and the old tower wake models to be used
!   CALL ReadStr( UnIn, InitInp%ADFileName, LINE, VarName='NewTowerModel?', VarDescr='Check for tower influence model', ErrStat=ErrStat )
   CALL ReadVar( UnIn, InitInp%ADFileName, LINE, VarName='NewTowerModel?', VarDescr='Check for tower influence model', ErrStat=ErrStat, ErrMsg=ErrMess )
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(UnIn)
      RETURN
   END IF

      ! Check if this is the "special string" to indicate the new tower influence model
   CALL Conv2UC( Line )
   IF ( INDEX(Line, "NEWTOWER" ) > 0 ) THEN

      !----------------------------------------------------------------------------------------------
      ! New tower influence model, as implemented by PJM
      !----------------------------------------------------------------------------------------------
      P%TwrProps%PJM_Version = .TRUE.

         ! Read in the tower potential flow switch
      CALL ReadVar( UnIn, InitInp%ADFileName, P%TwrProps%TwrPotent, VarName='TwrPotent', VarDescr='Tower influence model', ErrStat=ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

         ! Read in the tower shadow switch
      CALL ReadVar( UnIn, InitInp%ADFileName, P%TwrProps%TwrShadow, VarName='TwrShadow', VarDescr='Tower shadow model', ErrStat=ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

         ! Read in the tower drag file name
      CALL ReadVar( UnIn, InitInp%ADFileName, P%TwrProps%TwrFile, VarName='TwrFile', VarDescr='Tower data file name', ErrStat=ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

      IF ( PathIsRelative( P%TwrProps%TwrFile ) ) P%TwrProps%TwrFile = TRIM(FilePath)//TRIM(P%TwrProps%TwrFile)

         ! Read in the flag to tell AeroDyn to compute tower aerodynamics.
      CALL ReadVar( UnIn, InitInp%ADFileName, P%TwrProps%CalcTwrAero, VarName='CalcTwrAero', VarDescr='Flag to calculate tower aerodynamics', ErrStat=ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

   ELSE
      !----------------------------------------------------------------------------------------------
      ! Old tower influence model, read TwrShad from Line for now
      !----------------------------------------------------------------------------------------------
      P%TwrProps%PJM_Version = .FALSE.

      P%TwrProps%TwrPotent   = .FALSE.     ! We don't want to read the tower file!
      P%TwrProps%TwrShadow   = .FALSE.     ! We don't want to read the tower file!
      P%TwrProps%CalcTwrAero = .FALSE.     ! We don't want to read the tower file!

         ! Read in the tower shadow deficit
      IF ( INDEX( 'FTft', Line(:1) ) > 0 )  THEN
         CALL ProgWarn( ' Invalid numeric input. "'//TRIM( Line )//'" found when trying to read TwrShad.' )
         close(unin)
         ErrStat = ErrID_Fatal
         RETURN
      ELSE
         READ (Line,*,IOSTAT=ErrStat)  P%TwrProps%TwrShad
         CALL CheckIOS ( ErrStat, InitInp%ADFileName, 'TwrShad', NumType, ErrStat, ErrMess )
!bjj: is this aborting?
         IF ( p%Echo )  THEN
            WRITE (p%UnEc,"( 2X, ES11.4e2, 2X, A, T30, ' - ', A )")  P%TwrProps%TwrShad, "TwrShad", 'Tower shadow deficit'
         END IF

      END IF
!----------------

      IF ( P%TwrProps%TwrShad >= 1.0 ) THEN
         CALL ProgWarn( ' Tower shadow deficit cannot be >= 1.  Setting default deficit = 0.3' )
         P%TwrProps%TwrShad = 0.3
      END IF


         ! Read in the tower shadow width
      CALL ReadVar( UnIn, InitInp%ADFileName, P%TwrProps%ShadHWid, VarName='ShadHWid', VarDescr='Tower shadow half width', ErrStat=ErrStat, ErrMsg=ErrMess)
      IF ( ErrStat /= ErrID_None ) RETURN

      IF ( P%TwrProps%ShadHWid <= 0.0 ) THEN
         CALL ProgWarn( ' Tower shadow width cannot be <= zero.  Setting default half width = 1.0' )
         P%TwrProps%ShadHWid = 1.0
      END IF


         ! Read in the tower shadow reference point (distance from yaw axis to hub)
      CALL ReadVar( UnIn, InitInp%ADFileName, P%TwrProps%T_Shad_Refpt, VarName='T_Shad_Refpt', VarDescr='Tower shadow reference point', ErrStat=ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

         ! Constants used in tower shadow calculations
      P%TwrProps%TShadC1 = P%TwrProps%ShadHWid / SQRT ( ABS( P%TwrProps%T_Shad_Refpt ) )
      P%TwrProps%TShadC2 = P%TwrProps%TwrShad  * SQRT ( ABS( P%TwrProps%T_Shad_Refpt ) )

   END IF

      ! Read in the air density
   CALL ReadVar( UnIn, InitInp%ADFileName, P%Wind%Rho, VarName='Rho', VarDescr='Air density', ErrStat=ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

   !bjj do we need to check the the air density is non-negative?

   IF ( P%Wind%Rho == 0.0 .AND. P%DynInfl ) THEN  ! Turn off the GDW if RHO = 0.  It will crash
      CALL ProgWarn( 'Air density is zero. Dynamic Inflow will be turned off to avoid program crash.' )
      P%DynInfl = .FALSE.
   ENDIF

      ! Read in the kinematic viscosity
   CALL ReadVar( UnIn, InitInp%ADFileName, P%Wind%KinVisc, 'KinVisc', 'Kinematic viscosity', ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

      ! Aero calculation time interval
   CALL ReadVar( UnIn, InitInp%ADFileName, Line, 'DtAero', 'Aero calculation time step', ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "DEFAULT" ) /= 1 ) THEN ! If it's not "default", read this variable; otherwise use the value already stored in P%DtAero
         READ( Line, *, IOSTAT=ErrStatLcl) P%DtAero
         IF ( ErrStatLcl /= 0 ) THEN
            CALL CheckIOS ( ErrStatLcl, InitInp%ADFileName, "DT", NumType, ErrStat, ErrMess )
            RETURN
         ELSE
            ErrStat = ErrID_None
         END IF
      END IF  
      
            
      ! Read the number of airfoil files
   CALL ReadVar( UnIn, InitInp%ADFileName, P%AirFoil%NumFoil, 'NumFoil', 'Number of airfoil files', ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

   IF ( P%AirFoil%NumFoil < 1 ) THEN
      CALL ProgWarn( ' Error: Number of airfoil files must be a positive integer.')
      ErrStat = ErrID_Fatal
      close(unin)
      RETURN
   END IF

      !..............................................................................................
      ! Allocate space for the airfoil data file name(s), then read them
      !..............................................................................................
   IF (.NOT. ALLOCATED(P%AirFoil%FoilNm)) THEN
      ALLOCATE ( P%AirFoil%FoilNm( P%AirFoil%NumFoil ) , STAT=ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL ProgWarn(' Error allocating space for the FoilNm array.')
         close(unin)
         RETURN
      END IF
   END IF

   CALL ReadAryLines( UnIn, InitInp%ADFileName, P%AirFoil%FoilNm, P%AirFoil%NumFoil, AryName='FoilNm', AryDescr='Airfoil file names', ErrStat=ErrStat, ErrMsg=ErrMess )
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

   DO K=1,P%AirFoil%NumFoil
      IF ( PathIsRelative( P%AirFoil%FoilNm(K) ) ) P%AirFoil%FoilNm(K) = TRIM(FilePath)//TRIM( P%AirFoil%FoilNm(K) )
   END DO


      ! Read in the number of blade elements
   CALL ReadVar( UnIn, InitInp%ADFileName, P%Element%NElm, VarName='NElm', VarDescr='Number of blade elements', ErrStat=ErrStat, ErrMsg=ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF


      !..............................................................................................
      ! Allocate space for blade element data and read the arrays
      ! Read blade element data, check some inputs, convert twist to radians
      !..............................................................................................
   CALL AllocArrays (InitInp, P, x, xd, z, m, y, 'Element')

   m%ElOut%NumElOut    = 0 ! Initialize the element print array index
   m%ElOut%NumWndElOut = 0

   CALL ReadCom( UnIn, InitInp%ADFileName, 'Element table headers', ErrStat, ErrMess)
      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(UnIn)
         RETURN
      END IF

   DO IElm = 1, p%Element%NElm

      READ(UnIn,'(A)',IOSTAT=ErrStat) Line      !read into a line to see if print/no print is enabled

      IF (ErrStat == 0) THEN
         READ(Line,*,IOSTAT=ErrStat) P%Element%RElm(IElm), P%Element%Twist(IElm), P%Blade%DR(IElm), P%Blade%C(IElm), P%AirFoil%NFoil(IElm)
      END IF

      IF ( ErrStat == 0 ) THEN

         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! Check if AeroDyn will print out the element and/or wind data
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

         CALL Conv2UC(LINE)

         m%ElOut%ElPrList(IElm) = 0               ! INITIALIZE
         IndPrint = INDEX(LINE,"PRINT")
         IF (IndPrint > 0) THEN
            IF (LINE(IndPrint-2:IndPrint+4) /= "NOPRINT") THEN
               m%ElOut%NumElOut = m%ElOut%NumElOut + 1
               m%ElOut%ElPrList(IElm) = m%ElOut%NumElOut
            END IF
         ENDIF


         m%ElOut%WndElPrList(IElm) = 0            ! INITIALIZE
         IndPrint = INDEX(LINE,"WIND")
         IF (IndPrint > 0) THEN
            IF (LINE(IndPrint-2:IndPrint-1) /= "NO") THEN
               m%ElOut%NumWndElOut = m%ElOut%NumWndElOut + 1
               m%ElOut%WndElPrList(IElm) = m%ElOut%NumWndElOut
            END IF
         ENDIF

         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! Echo data to the file NWTC_Library echo file, if requested
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

         IF ( p%Echo ) THEN     ! NWTC_Library echo file
            WRITE (p%UnEc,'(4(2X,ES11.4e2),2X,I11,1(2X,L11))')  P%Element%RElm(IElm), P%Element%Twist(IElm), P%Blade%DR(IElm), P%Blade%C(IElm), P%AirFoil%NFoil(IElm), &
                   m%ElOut%ElPrList(IElm) /= 0  !, m%ElOut%WndElPrList(IElm) == 0
         END IF

      ELSE IF ( ErrStat < 0 ) THEN

         CALL ProgWarn( ' Premature end of file while reading line '//TRIM(Int2Lstr(IElm))// &
                     ' of the AeroDyn element table in file "'//TRIM(InitInp%ADFileName)//'."' )
         close(unin)
         ErrStat = ErrID_Fatal
         ErrMess = 'Error reading from line '//TRIM(Int2Lstr(IElm))//' of the AeroDyn element table.'
         RETURN
      ELSE
         close(unin)
         ErrStat = ErrID_Fatal
         ErrMess = 'Error reading from line '//TRIM(Int2Lstr(IElm))// &
                       ' of the AeroDyn element table in file "'//TRIM(InitInp%ADFileName)//'."'
         RETURN
      END IF


         ! Check for valid data:

      IF( P%Blade%C(IElm) <= 0.0 ) THEN
         CALL ProgWarn(' Error reading from line '//TRIM(Int2Lstr(IElm))//' of the AeroDyn element table.'// &
                       ' Chord length must be larger than 0.' )
         ErrStat = ErrID_Fatal
         close(unin)
         RETURN
      ENDIF

      IF (p%AirFoil%NFoil(IElm) < 1 .OR. p%AirFoil%NFoil(IElm) > p%AirFoil%NumFoil) THEN
         CALL ProgWarn(' Error reading from line '//TRIM(Int2Lstr(IElm))//' of the AeroDyn element table.'// &
                       ' Airfoil file ID must be a number between 1 and '//TRIM(Int2Lstr(p%AirFoil%NumFoil))//'.' )
         ErrStat = ErrID_Fatal
         close(unin)
         RETURN
      END IF

         ! Convert Twist to radians:

      P%Element%Twist(IElm) = P%Element%Twist(IElm)*D2R

   ENDDO ! IELM


      !..............................................................................................
      ! Read multiple airfoil table option
      !..............................................................................................

      PremEOF_indicator = .FALSE.
      READ(UnIn,*,IOSTAT=ErrStatLcl) Line         !read MultiTab -- it may not exist
      IF (ErrStatLcl > 0 ) THEN
         CALL WrScr1 ( ' Invalid character input for file "'//TRIM( InitInp%ADFileName )//'".' )
         CALL ProgWarn ( ' The error occurred while trying to read "MultiTab".' )
         ErrStat=ErrID_Fatal
         close(unin)      
         RETURN
      ELSE IF (ErrStatLcl == 0) THEN
         IF ( p%Echo )  THEN
            WRITE (p%UnEc, "( 15X, A, T30, ' - ', A, /, 2X, A )" )  &
                         'MultiTab', 'Multiple airfoil table option', '"'//TRIM( Line )//'"'
         END IF
      ELSE
         p%MultiTab = .FALSE.
         p%Reynolds = .FALSE.
         PremEOF_indicator = .TRUE.
      !  CALL PremEOF ( TRIM( Fil ), Variable, TrapThisError )
      END IF

   !-------------------------------------------------------------------------------------------------
   ! Close AeroDyn input file
   !-------------------------------------------------------------------------------------------------
   CLOSE(UnIn)

   !-------------------------------------------------------------------------------------------------
   ! Read airfoil data and check for MultiTab values using LINE, which was read above
   !-------------------------------------------------------------------------------------------------
   CALL READFL(InitInp, P, x, xd, z, m, y, ErrStatLcl, ErrMessLcl)     
   CALL SetErrStat(ErrStatLcl, ErrMessLcl, ErrStat, ErrMess,'AD_GetInput')   
   IF ( ErrStat >= AbortErrLev ) THEN
      close(unin)
      RETURN
   END IF
      
   m%AirFoil%MulTabLoc = 0.0                                    ! Initialize this value


            ! Read in the type of airfoil data table in each file
   IF ( PremEOF_indicator ) THEN             ! If we hit the end of the file without MultiTab, use only 1 airfoil table
      IF ( ANY( p%AirFoil%NTables(1:p%AirFoil%NumFoil) > 1 ) ) THEN
         CALL ProgWarn( ' Error reading multiple airfoil table option. Only one table for each file will be used.' )
      END IF
      p%MultiTab = .FALSE.
      p%Reynolds = .FALSE.
   ELSE ! not PremEOF_indicator

      IF ( ANY( p%AirFoil%NTables(1:p%AirFoil%NumFoil) > 1 ) ) THEN
         CALL Conv2UC(LINE(1:6))

         SELECT CASE ( TRIM(Line) )
            CASE ( 'USER' )
               p%MultiTab = .TRUE.
               p%Reynolds = .FALSE.

               DO K = 1, p%AirFoil%NumFoil
                  IF ( p%AirFoil%NTables(K) > 1 ) THEN
                     IF ( ( m%AirFoil%MulTabLoc < p%AirFoil%MulTabMet(K,1) ) .OR. &
                          ( m%AirFoil%MulTabLoc > p%AirFoil%MulTabMet(K,p%AirFoil%NTables(K) ) ))THEN
                        CALL ProgWarn( 'Error interpolating between airfoil tables. '// &
                                 ' Initial interpolation value = '//TRIM(Num2LStr(m%AirFoil%MulTabLoc))// &
                                 ' is outside table range of '//TRIM(Num2LStr(p%AirFoil%MulTabMet(K,1)))// &
                                 ' to '//TRIM(Num2LStr(p%AirFoil%MulTabMet(K,p%AirFoil%NTables(K))))// &
                                 ' in airfoil file #'//TRIM(Int2LStr(K))//'.' )
                        ErrStat = ErrID_Fatal
                        RETURN
                     END IF
                  END IF ! NTables(K) > 1
               ENDDO ! K

            CASE ( 'RENUM' )
               p%MultiTab = .TRUE.
               p%Reynolds = .TRUE.
            CASE ( 'SINGLE' )
               p%MultiTab = .FALSE.
               p%Reynolds = .FALSE.
            CASE DEFAULT
               CALL WrScr( ' Error: control model option must be "USER", "RENUM" or "SINGLE".' )
         END SELECT
      ELSE
         p%MultiTab = .FALSE.
         p%Reynolds = .FALSE.
      END IF

   ENDIF

   !-------------------------------------------------------------------------------------------------
   ! Read tower drag input file, if necessary
   !-------------------------------------------------------------------------------------------------
   IF (p%TwrProps%TwrPotent .OR. p%TwrProps%TwrShadow .OR. p%TwrProps%CalcTwrAero) THEN                 ! Read in the tower drag file
      CALL READTwr(UnIn, InitInp, P, x, xd, z, m, y, ErrStat, ErrMess )
      IF (ErrStat /= ErrID_None) RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Initialize variables for printing
   !-------------------------------------------------------------------------------------------------
   IF ( m%ElOut%NumElOut > 0 .OR. m%ElOut%NumWndElOut > 0 ) THEN
      p%ElemPrn = .TRUE.
      CALL AllocArrays (InitInp, P, x, xd, z, m, y, 'ElPrint')

      ElIndex = 0                      ! Re-Initialize the element print array index for wind
      DO IElm = 1, p%Element%NElm
         IF (m%ElOut%WndElPrList(IElm) > 0) THEN
            ElIndex = ElIndex + 1
            m%ElOut%WndElPrNum(ElIndex) = IElm
         END IF
      END DO ! IELM

      ElIndex = 0 ! Re-Initialize the element print array index
      DO IElm = 1, p%Element%NElm
         IF (m%ElOut%ElPrList(IElm) > 0) THEN
            ElIndex = ElIndex + 1
            m%ElOut%ElPrNum(ElIndex) = IElm
         END IF
      END DO ! IELM
   ELSE
      p%ElemPrn = .FALSE.
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Initialize Beddoes dynamic stall data
   !-------------------------------------------------------------------------------------------------
   IF ( p%DStall ) CALL BEDDAT( P, x, xd, z, m, y, ErrStat, ErrMess )


   RETURN
   
contains 
   SUBROUTINE CleanUp()
   
      CLOSE(UnIn)
   
   END SUBROUTINE CleanUp

END SUBROUTINE AD14_GetInput
!====================================================================================================
   SUBROUTINE ADOut(InitInp, P, m, AD14_Ver, FileName, ErrStat, ErrMess )
 !  used to output data to a summary file
 ! *****************************************************
   IMPLICIT                    NONE
      ! Passed Variables:
   TYPE(AD14_InitInputType),       INTENT(IN   )  :: InitInp
   TYPE(AD14_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(IN   )  :: m           ! Misc/optimization variables
   TYPE(ProgDesc),                 INTENT(IN   )  :: AD14_ver
   INTEGER,                        INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMess
   CHARACTER(*),                   INTENT(IN   )  :: FileName

   

   ! Local Variables:

   INTEGER                                :: IElm
   INTEGER                                :: IFoil
   INTEGER                                :: UnOut
   INTEGER                                :: I,K
   INTEGER(IntKi)                         :: ErrStatLcl
                                          
   CHARACTER(  2)                         :: Dst_Unit
   CHARACTER(150)                         :: Frmt
   CHARACTER(  4)                         :: Mass_Unit
   CHARACTER( 35)                         :: MESAGE
   CHARACTER(  3)                         :: Vel_Unit
                                          
   CHARACTER(1),PARAMETER                 :: Delim = ' '  ! bjj: made this a parameter because I don't think tabs work very well in a summary file
   CHARACTER(ErrMsgLen)                   :: ErrMessLcl


   ErrStat = ErrID_None
   ErrMess = ""

   ! Function definition

   CALL GetNewUnit( UnOut, ErrStat, ErrMess )
   CALL OpenFOutFile( UnOut, FileName, ErrStatLcl, ErrMessLcl)
      CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,'ADOut' )
      IF ( ErrStat >= AbortErrLev )  THEN
         CLOSE(UnOut)
         RETURN
      END IF


      WRITE (UnOut,"(/A)")  'This file was generated by '//TRIM(GetNVD(AD14_Ver))//&
                                ' on '//CurDate()//' at '//CurTime()//'.'


IF (p%SIUNIT) THEN
   Dst_Unit  = 'm'
   Mass_Unit = 'kg'
   Vel_Unit  = 'mps'
ELSE
   Dst_Unit  = 'ft'
   Mass_Unit = 'slug'
   Vel_Unit  = 'fps'
ENDIF

 ! Reiterate the input file
WRITE(UnOut,'(/A)') 'Inputs read in from the AeroDyn input file:'
WRITE(UnOut,'(A/)') TRIM(p%TITLE)

! Ec_Ch11Frmt is a parameter defined for echo output in the NWTC Subroutine Library
MESAGE = 'Units for input and output'
IF ( p%SIUNIT ) THEN
   WRITE(UnOut,Ec_Ch11Frmt) 'SI','SysUnits',MESAGE
ELSE
   WRITE(UnOut,Ec_Ch11Frmt) 'ENGLISH','SysUnits',MESAGE
ENDIF

MESAGE = 'Dynamic stall model'
IF ( p%DSTALL ) THEN
   WRITE(UnOut,Ec_Ch11Frmt) 'BEDDOES',"StallMod", MESAGE//' [Beddoes]'
ELSE
   WRITE(UnOut,Ec_Ch11Frmt) 'STEADY',"StallMod", MESAGE//' [NO Dynamic stall]'
ENDIF

MESAGE = 'Aerodynamic pitching moment model'
IF ( p%PMOMENT ) THEN
   WRITE(UnOut,Ec_Ch11Frmt) 'USE_CM','UseCm',MESAGE//' [Pitching Moments calculated]'
ELSE
   WRITE(UnOut,Ec_Ch11Frmt) 'NO_CM','UseCm',MESAGE//' [NO Pitching Moments calculated]'
ENDIF

MESAGE = 'Inflow model'
IF ( p%DYNINFL ) THEN
   WRITE(UnOut,Ec_Ch11Frmt) 'DYNIN','InfModel',MESAGE//' [Dynamic Inflow]'
ELSE
   IF ( p%InducedVel%EquilDA .AND. p%InducedVel%EquilDT )  THEN
      WRITE(UnOut,Ec_Ch11Frmt) 'EQUILDB','InfModel',MESAGE//' [Equilibrium w/ axial and tangential drag]'
   ELSEIF ( p%InducedVel%EquilDA )  THEN
      WRITE(UnOut,Ec_Ch11Frmt) 'EQUILDA','InfModel',MESAGE//' [Equilibrium w/ axial drag]'
   ELSEIF ( p%InducedVel%EquilDT )  THEN
      WRITE(UnOut,Ec_Ch11Frmt) 'EQUILDT','InfModel',MESAGE//' [Equilibrium w/ tangential drag]'
   ELSE
      WRITE(UnOut,Ec_Ch11Frmt) 'EQUIL','InfModel',MESAGE//' [Equilibrium]'
   ENDIF
ENDIF


MESAGE = 'Induction factor model'
IF ( p%WAKE ) THEN
   IF (p%SWIRL) THEN
      WRITE(UnOut,Ec_Ch11Frmt) 'SWIRL','IndModel',MESAGE//' [Normal and Radial flow induction factors calculated]'
   ELSE
      WRITE(UnOut,Ec_Ch11Frmt) 'WAKE','IndModel',MESAGE//' [Normal flow induction factors calculated]'
   ENDIF
   WRITE(UnOut,Ec_ReFrmt) p%InducedVel%ATOLER,'AToler','Convergence tolerance for induction factor'
ELSE
   WRITE(UnOut,Ec_Ch11Frmt) 'NONE','IndModel',MESAGE//' [NO induction factors calculated]'
   WRITE(UnOut,Ec_Ch11Frmt) '[Not Used]','AToler','Convergence tolerance for induction factor'
ENDIF

MESAGE = 'Tip-loss model'
IF (.NOT. p%DYNINFL) THEN
   IF ( p%InducedVel%TLOSS ) THEN
      IF (p%InducedVel%GTECH) THEN
         WRITE(UnOut,Ec_Ch11Frmt) 'GTECH','TLModel',MESAGE//' [Georgia Tech correction to Prandtl model]'
      ELSE
         WRITE(UnOut,Ec_Ch11Frmt) 'PRAND','TLModel',MESAGE//' [Prandtl model]'
      ENDIF
   ELSE
      WRITE(UnOut,Ec_Ch11Frmt) 'NONE','TLModel',MESAGE//' [NO tip-loss calculated]'
   ENDIF
ELSE
   WRITE(UnOut,Ec_Ch11Frmt) '[Not Used]','TLModel',MESAGE
ENDIF

MESAGE = 'Hub-loss model'
IF (.NOT. p%DYNINFL) THEN
   IF ( p%InducedVel%HLOSS ) THEN
      WRITE(UnOut,Ec_Ch11Frmt) 'PRAND','HLModel',MESAGE//' [Prandtl model]'
   ELSE
      WRITE(UnOut,Ec_Ch11Frmt) 'NONE','HLModel',MESAGE//' [NO hub-loss calculated]'
   ENDIF
ELSE
   WRITE(UnOut,Ec_Ch11Frmt) '[Not Used]','HLModel',MESAGE
ENDIF


WRITE(UnOut,Ec_ReFrmt) p%Rotor%HH,'HH','Wind reference (hub) height, '//TRIM(Dst_Unit)

IF ( p%TwrProps%PJM_Version ) THEN
   WRITE(UnOut,Ec_LgFrmt) p%TwrProps%TwrPotent,'TwrPotent','Calculate tower potential flow [T or F]'
   WRITE(UnOut,Ec_LgFrmt) p%TwrProps%TwrShadow,'TwrShadow','Calculate tower shadow [T or F]'
   IF ( p%TwrProps%TwrPotent .OR. p%TwrProps%TwrShadow ) THEN
      WRITE(UnOut,Ec_StrFrmt) 'TwrFile','Tower drag file name',TRIM(P%TwrProps%TwrFile)
   ELSE
      WRITE(UnOut,Ec_Ch11Frmt) '[none]','TwrFile','No tower drag properties file'
   ENDIF
ELSE
   WRITE(UnOut,Ec_ReFrmt) p%TwrProps%TwrShad,'TwrShad','Tower shadow centerline velocity deficit'
   WRITE(UnOut,Ec_ReFrmt) p%TwrProps%ShadHWid,'ShadHWid','Tower shadow half width, '//TRIM(Dst_Unit)
   WRITE(UnOut,Ec_ReFrmt) p%TwrProps%T_Shad_Refpt,'T_Shad_Refpt','Tower shadow reference point, '//TRIM(Dst_Unit)
END IF


WRITE(UnOut,Ec_ReFrmt)  p%Wind%RHO,'AirDens','Air density, '//TRIM(Mass_Unit)//'/'//TRIM(Dst_Unit)//'^3'
WRITE(UnOut,Ec_ReFrmt)  p%Wind%KinVisc,'KinVisc','Kinematic air viscosity, '//TRIM(Dst_Unit)//'^2/sec'
WRITE(UnOut,Ec_ReFrmt)  p%DTAERO,'DTAERO','Time interval for aerodynamic calculations, sec'
WRITE(UnOut,Ec_IntFrmt) p%AirFoil%NUMFOIL,'NumFoil','Number of airfoil files used. Files listed below:'

DO IFoil = 1, p%AirFoil%NUMFOIL
   WRITE(UnOut,'(A)') '"'//TRIM(p%AirFoil%FOILNM(IFoil))//'"'
END DO ! IFoil

WRITE(UnOut,Ec_IntFrmt) p%Element%NELM,'BldNodes','Number of blade elements per blade'

   !-------------------------------------------------------------------------------------------------
   ! write out element information
   !-------------------------------------------------------------------------------------------------
   Frmt = '(3X,A10,8("'//Delim//'",A10))'

   WRITE(UnOut,'( )')

      ! column names

   WRITE(UnOut,Frmt) '  Element ', &
                     '   RELM   ', &
                     '   Twist  ', &
                     '    DR    ', &
                     '   Chord  ', &
                     '   NFoil  ', &
                     '  Print?  ', &
                     ' Tip-loss ', &
                     ' Hub-loss '

      ! column units

   WRITE(UnOut,Frmt) '    (-)   ', &
                     '    (m)   ', &
                     '   (deg)  ', &
                     '    (m)   ', &
                     '    (m)   ', &
                     '    (-)   ', &
                     ' (Yes/No) ', &
                     ' constant ', &
                     ' constant '

   WRITE(UnOut,Frmt) '----------', &
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

   DO IElm = 1, p%Element%NELM

      IF (m%ElOut%ElPrList(IElm) /= 0) THEN
         MESAGE = 'Yes'
      ELSE
         MESAGE = 'No'
      ENDIF

      WRITE(UnOut, Frmt) IElm, p%Element%RELM(IElm), p%Element%TWIST(IElm)*R2D, p%Blade%DR(IElm),  p%Blade%C(IElm), &
                         p%AirFoil%NFOIL(IElm), TRIM(Mesage), p%Element%TLCNST(IElm), p%Element%HLCNST(IElm)
   END DO


IF ( p%MultiTab ) THEN
   WRITE(UnOut,'(A)') 'MULTI    Multiple airfoil tables used'
ELSE
   WRITE(UnOut,'( )')
ENDIF

      WRITE(UnOut,"(/' Rotor radius     = ',F7.3,' m')") p%Blade%R
      WRITE(UnOut,"( ' Hub radius       = ',F7.3,' m')") p%HubRad
      WRITE(UnOut,"( ' Number of blades = ',I3       )") p%NumBl

IF ( p%DSTALL ) THEN
   Frmt = '(3X,A, 21(:F8.4,3X) )'
   DO K = 1, P%AirFoil%NTables(1)
      WRITE(UnOut,'(/A/)') '  BEDDOES DYNAMIC STALL PARAMETERS:'
      WRITE(UnOut, Frmt) 'CN SLOPE         ', ( m%Beddoes%CNA(I,K),      I = 1, P%AirFoil%NUMFOIL )
      WRITE(UnOut, Frmt) 'STALL CN (UPPER) ', ( m%Beddoes%CNS(I,K),      I = 1, P%AirFoil%NUMFOIL )
      WRITE(UnOut, Frmt) 'STALL CN (LOWER) ', ( m%Beddoes%CNSL(I,K),     I = 1, P%AirFoil%NUMFOIL )
      WRITE(UnOut, Frmt) 'ZERO LIFT AOA    ', ( m%Beddoes%AOL(I,K)*R2D,  I = 1, P%AirFoil%NUMFOIL )
      WRITE(UnOut, Frmt) 'MIN DRAG AOA     ', ( m%Beddoes%AOD(I,K)*R2D,  I = 1, P%AirFoil%NUMFOIL )
      WRITE(UnOut, Frmt) 'MIN DRAG COEFF   ', ( m%Beddoes%CDO(I,K),      I = 1, P%AirFoil%NUMFOIL )
      WRITE(UnOut,'(/)')
   ENDDO !K

   WRITE(UnOut,*) '    VORTEX TRANSIT TIME FROM LE TO TE ', P%Beddoes%TVL
   WRITE(UnOut,*) '    PRESSURE TIME CONSTANT            ', P%Beddoes%TP
   WRITE(UnOut,*) '    VORTEX TIME CONSTANT              ', P%Beddoes%TV
   WRITE(UnOut,*) '    F-PARAMETER TIME CONSTANT         ', P%Beddoes%TF   
END IF

IF ( p%ELEMPRN ) WRITE(UnOut,'(/A/)')'Blade element aerodynamic time series data written to file.'

CLOSE (UnOut )

RETURN
END SUBROUTINE ADOut

 ! ****************************************************
   SUBROUTINE READFL(InitInp, P, x, xd, z, m, y, ErrStat, ErrMess )
 !  Reads a data file containing airfoil angle of attack,
 !   CL and CD, and dynamic stall parameters
 ! ****************************************************
!====================================================================================================
   USE                           AeroGenSubs, ONLY: AllocArrays
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_InitInputType),       INTENT(INOUT)  :: InitInp
   TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD14_ContinuousStateType), INTENT(INOUT)  :: x           ! Initial continuous states
   TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(AD14_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMess

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

   INTEGER                    :: ErrStatLcL        ! Error status returned by called routines.
   CHARACTER(ErrMsgLen)       :: ErrMessLcl          ! Error message returned by called routines.



   ErrStat = ErrID_None
   ErrMess = ""

   CALL GetNewUnit(NUNIT, ErrStat, ErrMess) 
   p%AirFoil%NumCL    = 0

 ! The first loop checks existence and file length to set NumCL
DO NFOILID = 1, p%AirFoil%NUMFOIL

 ! Open the file for reading # of lines
   CALL OpenFInpFile (NUNIT, TRIM(p%AirFoil%FOILNM(NFOILID)), ErrStatLcL, ErrMessLcl)
   CALL SetErrStat( ErrStatLcL, ErrMessLcl, ErrStat, ErrMess, 'READFL')
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(NUNIT)
      RETURN
   END IF

 ! Determine the maximum number of aerodata points in all files

   NumLines = 0
   IOS = 0
   DO WHILE (IOS == 0)
      READ ( NUNIT, '()', IOSTAT=IOS )
      NumLines = NumLines + 1
   END DO

   p%AirFoil%NumCL = MAX(NumLines - 14, p%AirFoil%NumCL)

   CLOSE (NUNIT)

END DO ! NFOILID

 ! Allocate the arrays

CALL AllocArrays (InitInp, P, x, xd, z, m, y, 'Aerodata')
   !CALL SetErrStat( ErrStatLcL, ErrMessLcl, ErrStat, ErrMess, 'READFL')
   !IF (ErrStat >= AbortErrLev) RETURN

 ! The second loop reads the files
DO NFOILID = 1, p%AirFoil%NUMFOIL

 ! Open the file for reading inputs
   CALL OpenFInpFile (NUNIT, TRIM(Adjustl(p%AirFoil%FOILNM(NFOILID))), ErrStatLcL, ErrMessLcl )
   CALL SetErrStat( ErrStatLcL, ErrMessLcl, ErrStat, ErrMess, 'READFL')
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(NUNIT)
      RETURN
   END IF
   
 ! Set up the file to read the aerodata
   READ(NUNIT,'( A )',IOSTAT=IOS) TITLE(1)
   READ(NUNIT,'( A )',IOSTAT=IOS) TITLE(2)

 ! Read in airfoil table dimension parameters:
 !   NTables = number of airfoil data tables

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) THEN
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , '# of tables', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
   END IF
   

   READ(LINE,*,ERR=205) p%AirFoil%NTables( NFOILID )

 ! Allocate local arrays with NTables dimension

   Sttus = 0
   IF (.NOT. ALLOCATED(CLPosPI)) THEN
      ALLOCATE ( CLPosPI(P%AirFoil%NTables(NFOILID)) , STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Error allocating memory for CLPosPI array.', ErrStat, ErrMess, 'READFL' )
         close(NUNIT)
         RETURN
      END IF
   END IF
   

   IF (.NOT. ALLOCATED(CDPosPI)) THEN
      ALLOCATE ( CDPosPI(P%AirFoil%NTables(NFOILID)) , STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Error allocating memory for CDPosPI array.', ErrStat, ErrMess, 'READFL' )
         close(NUNIT)
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(CMPosPI)) THEN
      ALLOCATE ( CMPosPI(P%AirFoil%NTables(NFOILID)) , STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Error allocating memory for CMPosPI array.', ErrStat, ErrMess, 'READFL' )
         close(NUNIT)
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(CLNegPI)) THEN
      ALLOCATE ( CLNegPI(P%AirFoil%NTables(NFOILID)) , STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Error allocating memory for CLNegPI array.', ErrStat, ErrMess, 'READFL' )
         close(NUNIT)
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(CDNegPI)) THEN
      ALLOCATE ( CDNegPI(P%AirFoil%NTables(NFOILID)) , STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Error allocating memory for CDNegPI array.', ErrStat, ErrMess, 'READFL' )
         close(NUNIT)
         RETURN
      END IF
   END IF

   IF (.NOT. ALLOCATED(CMNegPI)) THEN
      ALLOCATE ( CMNegPI(P%AirFoil%NTables(NFOILID)) , STAT=Sttus )
      IF ( Sttus /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Error allocating memory for CMNegPI array.', ErrStat, ErrMess, 'READFL' )
         close(NUNIT)
         RETURN
      END IF
   END IF


 ! Read in airfoil data table identification array

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) then
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , 'multi-table metric', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   READ(LINE,*,ERR=205)  (p%AirFoil%MulTabMet ( NFOILID, K ), K = 1, p%AirFoil%NTables(NFOILID))

 ! Read in four lines that are no longer used
 ! These are retained for future USE and backwards compatibility only

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) then
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , '5th line', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) then
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , '6th line', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) then
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , '7th line', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) then
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , '8th line', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF

 ! Read Beddoes stall parameters

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) THEN
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , 'Angle of zero lift (AOL)', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   IF (p%DSTALL) READ(LINE,*,ERR=205)  (m%Beddoes%AOL( NFOILID, K ), K = 1, p%AirFoil%NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) THEN
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , 'CNA', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   IF (p%DSTALL) READ(LINE,*,ERR=205)  (m%Beddoes%CNA   ( NFOILID, K ), K = 1, p%AirFoil%NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) THEN
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , 'CNS', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   IF (p%DSTALL) READ(LINE,*,ERR=205)  (m%Beddoes%CNS   ( NFOILID, K ), K = 1, p%AirFoil%NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) THEN 
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , 'CNSL', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   IF (p%DSTALL) READ(LINE,*,ERR=205)  (m%Beddoes%CNSL  ( NFOILID, K ), K = 1, p%AirFoil%NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) THEN
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , 'AOD', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   IF (p%DSTALL) READ(LINE,*,ERR=205)  (m%Beddoes%AOD   ( NFOILID, K ), K = 1, p%AirFoil%NTables(NFOILID))

   READ(NUNIT,'( A )',IOSTAT=IOS) LINE
   IF ( IOS < 0 ) THEN
      CALL PremEOF ( Trim(p%AirFoil%FOILNM(NFOILID)) , 'CDO', .TRUE., ErrMessLcl )
      close(NUNIT)
      CALL SetErrStat( ErrID_Fatal, ErrMessLcl, ErrStat, ErrMess, 'READFL')
      RETURN
   END IF
   IF (p%DSTALL) READ(LINE,*,ERR=205)  (m%Beddoes%CDO   ( NFOILID, K ), K = 1, p%AirFoil%NTables(NFOILID))

 ! Convert angles to radians
   IF (p%DSTALL) THEN
      m%Beddoes%AOD   ( NFOILID, : ) = m%Beddoes%AOD( NFOILID, : )*D2R
      m%Beddoes%AOL   ( NFOILID, : ) = m%Beddoes%AOL( NFOILID, : )*D2R
   ENDIF


 ! Read airfoil data tables to end of file

   p%AirFoil%NLIFT(NFOILID) = 0
   ALPosPI = .FALSE.
   ALNegPI = .FALSE.

   DO I = 1, p%AirFoil%NumCL

      IF ( p%PMOMENT ) THEN

         READ( NUNIT,*,END=150,ERR=150 ) m%AirFoil%AL(NFOILID,I), &
             (m%AirFoil%CL(NFOILID,I,IPHI), m%AirFoil%CD(NFOILID,I,IPHI), &
              m%AirFoil%CM(NFOILID,I,IPHI), IPHI = 1, p%AirFoil%NTables(NFOILID))

      ELSE

         READ( NUNIT,*,END=150,ERR=150 ) m%AirFoil%AL(NFOILID,I), &
             (m%AirFoil%CL(NFOILID,I,IPHI), m%AirFoil%CD(NFOILID,I,IPHI), &
              IPHI = 1, p%AirFoil%NTables(NFOILID))

         m%AirFoil%CM(NFOILID,I,:) = 0.

      ENDIF

 ! Check to see if values look reasonable

      DO IPHI = 1, p%AirFoil%NTables(NFOILID)
        IF ( ABS( m%AirFoil%AL( NFOILID, I ) ) > 185.) THEN
            CALL SetErrStat( ErrID_Fatal, 'Probable error in airfoil data table number '//TRIM(Int2LStr(NFOILID))// &
                           ' Angle of attack exceeds 185 degrees.', ErrStat, ErrMess, 'READFL')   
            CLOSE(NUNIT)
            RETURN

        ELSEIF (ABS( m%AirFoil%CL( NFOILID, I, IPHI ) ) > 3. ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Probable error in airfoil data table number '//TRIM(Int2LStr(NFOILID))// &
                           ' Coefficient of Lift exceeds 3.0.', ErrStat, ErrMess, 'READFL')   
            CLOSE(NUNIT)
            RETURN
        ELSEIF (ABS( m%AirFoil%CD( NFOILID, I, IPHI ) ) > 3. ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Probable error in airfoil data table number '//TRIM(Int2LStr(NFOILID))// &
                           ' Coefficient of Drag exceeds 3.0.', ErrStat, ErrMess, 'READFL')   
            CLOSE(NUNIT)
            RETURN
        ELSEIF (ABS( m%AirFoil%CM( NFOILID, I, IPHI ) ) > 3. ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Probable error in airfoil data table number '//TRIM(Int2LStr(NFOILID))// &
                           ' Coefficient of Moment exceeds 3.0.', ErrStat, ErrMess, 'READFL')   
            CLOSE(NUNIT)
            RETURN
           
           
        ENDIF
      ENDDO ! IPHI

 ! Store the values at 180 deg. and -180 deg. for check
      IF ( m%AirFoil%AL (NFOILID, I ) == 180. ) THEN
         ALPosPI = .TRUE.
         Do IPHI = 1, p%AirFoil%NTables(NFOILID)
            CLPosPI(IPHI) = m%AirFoil%CL(NFOILID,I,IPHI)
            CDPosPI(IPHI) = m%AirFoil%CD(NFOILID,I,IPHI)
            CMPosPI(IPHI) = m%AirFoil%CM(NFOILID,I,IPHI)
         END Do ! IPHI

      ELSEIF ( m%AirFoil%AL (NFOILID, I ) == -180. ) THEN
         ALNegPI = .TRUE.
         Do IPHI = 1, p%AirFoil%NTables(NFOILID)
            CLNegPI(IPHI) = m%AirFoil%CL(NFOILID,I,IPHI)
            CDNegPI(IPHI) = m%AirFoil%CD(NFOILID,I,IPHI)
            CMNegPI(IPHI) = m%AirFoil%CM(NFOILID,I,IPHI)
          END Do ! IPHI
      ENDIF

      m%AirFoil%AL ( NFOILID, I ) = m%AirFoil%AL(NFOILID,I) * D2R
      p%AirFoil%NLIFT ( NFOILID ) = p%AirFoil%NLIFT(NFOILID) + 1

   ENDDO ! I

   150 CLOSE( NUNIT )

 ! Check to see if values at 180 deg. equal those at -180 deg.
   IF (ALPosPI .AND. ALNegPI) THEN
      Do IPHI = 1, p%AirFoil%NTables(NFOILID)
         IF (CLPosPI(IPHI) /= CLNegPI(IPHI) .OR. &
             CDPosPI(IPHI) /= CDNegPI(IPHI) .OR. &
             CMPosPI(IPHI) /= CMNegPI(IPHI)) THEN
            CALL SetErrStat( ErrID_Fatal, ' The airfoil data at +180 deg is different from -180 deg in file :'//Trim(P%AirFoil%FOILNM(NFOILID)), ErrStat, ErrMess, 'READFL')            
            return
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

205 CALL SetErrStat( ErrID_Fatal, ' Error reading line: "'//TRIM(Line)//'" in file : "'//TRIM(P%AirFoil%FOILNM(NFOILID))//'"', ErrStat, ErrMess, 'READFL')   
   CLOSE(NUNIT)
   RETURN

RETURN
END SUBROUTINE READFL

 ! ****************************************************
   SUBROUTINE READTwr(UnIn, InitInp, P, x, xd, z, m, y, ErrStat, ErrMess )
! This subroutine reads the tower properties input file, allocating TwrProps variables to do so.
! The tower data file contains radius and Re vs CD data as well as the tower wake constant.
 ! ****************************************************
!====================================================================================================
   IMPLICIT                      NONE
      ! Passed Variables:
   INTEGER,                      INTENT(IN)     :: UnIn
   TYPE(AD14_InitInputType),       INTENT(INOUT)  :: InitInp
   TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD14_ContinuousStateType), INTENT(INOUT)  :: x           ! Initial continuous states
   TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(AD14_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
   INTEGER, INTENT(OUT)                          :: ErrStat
   CHARACTER(*), INTENT(OUT)                     :: ErrMess

      ! Local Variables:

   INTEGER                      :: I         ! loop counter for rows in the data tables
   INTEGER                      :: J         ! loop counter for columns in the data tables

   CHARACTER(99)                :: Fmt       ! format for printing to an echo file
   CHARACTER(1024)              :: FilName   ! file name

   !-------------------------------------------------------------------------------------------------
   ! Open the file for reading
   !-------------------------------------------------------------------------------------------------
   FilName = p%TwrProps%TwrFile
   CALL OpenFInpFile (UnIn, TRIM(FilName), ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN


   !-------------------------------------------------------------------------------------------------
   ! Read the heading, section 1
   !-------------------------------------------------------------------------------------------------

      ! Read in 2 header/comment lines
   CALL ReadCom( UnIn, FilName, 'Title line 1', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN

   CALL ReadCom( UnIn, FilName, 'Title line 2', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN


      ! Read in number of tower height entries, NTwrHt
   CALL ReadVar( UnIn, FilName, p%TwrProps%NTwrHt, 'NTwrHt', 'Number of tower stations', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN

   IF (p%TwrProps%NTwrHt < 1) THEN
      CALL ProgWarn( 'Number of tower height entries, NTwrHt, must be greater than zero.' )
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF

      ! Read in number of tower Reynolds number entries, NTwrRe
   CALL ReadVar( UnIn, FilName, p%TwrProps%NTwrRe, 'NTwrRe', 'Number of tower Reynolds number rows', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN

   IF (p%TwrProps%NTwrRe < 1) THEN
      CALL ProgWarn( 'Number of tower Reynolds number entries, NTwrRe, must be greater than zero.' )
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF


      ! Read in number of tower CD entries, NTwrCD
   CALL ReadVar( UnIn, FilName, p%TwrProps%NTwrCD, 'NTwrCD', 'Number of tower CD columns', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN

   IF (p%TwrProps%NTwrCD < 1) THEN
      CALL ProgWarn( 'Number of tower CD entries, NTwrCD, must be greater than zero.' )
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF


      ! Read in constant for tower wake model = 0 full potential flow = 0.1 model of Bak et al.
   CALL ReadVar( UnIn, FilName, p%TwrProps%Tower_Wake_Constant, 'Tower_Wake_Constant', 'Constant for tower wake model', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN

   ! bjj: should there be a sanity check here, too?

   !-------------------------------------------------------------------------------------------------
   ! Allocate TwrProps arrays with NTwrHt, NTwrRe, and NTwrCD dimensions; these arrays are
   ! read in the next 2 sections of this file.
   !-------------------------------------------------------------------------------------------------

   IF ( .NOT. ALLOCATED( p%TwrProps%TwrHtFr ) ) THEN
      ALLOCATE ( p%TwrProps%TwrHtFr(p%TwrProps%NTwrHt) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for TwrHtFr array.' )
         ErrStat = ErrID_Fatal
         RETURN
      END IF
   END IF

   IF ( .NOT. ALLOCATED( p%TwrProps%TwrWid ) ) THEN
      ALLOCATE ( p%TwrProps%TwrWid(p%TwrProps%NTwrHt) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for TwrWid array.' )
         ErrStat = ErrID_Fatal
         RETURN
      END IF
   END IF

   IF ( .NOT. ALLOCATED( p%TwrProps%NTwrCDCol ) ) THEN
      ALLOCATE ( p%TwrProps%NTwrCDCol(p%TwrProps%NTwrHt) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for NTwrCDCol array.' )
         ErrStat = ErrID_Fatal
         RETURN
      END IF
   END IF

   IF ( .NOT. ALLOCATED( p%TwrProps%TwrRe ) ) THEN
      ALLOCATE ( p%TwrProps%TwrRe(p%TwrProps%NTwrRe) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for TwrRe array.' )
         ErrStat = ErrID_Fatal
         RETURN
      END IF
   END IF

   IF ( .NOT. ALLOCATED( p%TwrProps%TwrCD ) ) THEN
      ALLOCATE ( p%TwrProps%TwrCD(p%TwrProps%NTwrRe, p%TwrProps%NTwrCD) , STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgWarn( ' Error allocating memory for TwrCD array.' )
         ErrStat = ErrID_Fatal
         RETURN
      END IF
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Read section 2, DISTRIBUTED TOWER PROPERTIES;
   ! section contains 2 heading lines in addition to NTwrHt rows of data with 3 columns
   !-------------------------------------------------------------------------------------------------
      ! Read in 2 header/comment lines
   CALL ReadCom( UnIn, FilName, 'Distributed Tower Properties header 1', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN

   CALL ReadCom( UnIn, FilName, 'Distributed Tower Properties header 2', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN


      ! Read tower height data table

   DO I = 1, p%TwrProps%NTwrHt     ! 1 line per fraction of height


      READ( UnIn,*,IOSTAT=ErrStat ) p%TwrProps%TwrHtFr(I), p%TwrProps%TwrWid(I), p%TwrProps%NTwrCDCol(I)

      IF ( ErrStat == 0 ) THEN
         IF ( p%Echo ) THEN
            WRITE (p%UnEc,'(2X,ES11.4e2, 2X,ES11.4e2, 2X,I11)')  p%TwrProps%TwrHtFr(I), p%TwrProps%TwrWid(I), p%TwrProps%NTwrCDCol(I)
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
      IF ( p%TwrProps%TwrHtFr( I ) < 0.0 .OR. p%TwrProps%TwrHtFr( I ) > 1.0 ) THEN
         CALL ProgWarn( ' Error on line '//TRIM(Int2Lstr(I))//' of the distributed tower properties in file "' &
                           //TRIM(FilName)//'." Tower height fractions must be between 0.0 and 1.0.' )
         ErrStat = ErrID_Fatal
         RETURN
      END IF


         ! Make sure the tower height increases for each entry
      IF (I > 1) THEN
         IF (p%TwrProps%TwrHtFr(I) <= p%TwrProps%TwrHtFr(I-1)) THEN
            CALL ProgWarn( ' Error on line '//TRIM(Int2Lstr(I))//' of the distributed tower properties in file "' &
                           //TRIM(FilName)//'." Tower height fraction entries must be in order of increasing height.' )
            ErrStat = ErrID_Fatal
            RETURN
         ENDIF
      ENDIF

         ! Make sure tower width is positive
      IF ( p%TwrProps%TwrWid( I ) <= 0.0 ) THEN
         CALL ProgWarn( ' Error on line '//TRIM(Int2Lstr(I))//' of the distributed tower properties in file "' &
                           //TRIM(FilName)//'." Tower width must be positive.' )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

         ! Make sure the tower CD column is within range
      IF ( p%TwrProps%NTwrCDCol(I) < 1 .OR. p%TwrProps%NTwrCDCol(I) > P%TwrProps%NTwrCD ) THEN
         CALL ProgWarn( ' Error on line '//TRIM(Int2Lstr(I))//' of the distributed tower properties in file "' &
                           //TRIM(FilName)//'." Tower height CD column must be between 1 and '//TRIM(Int2Lstr(P%TwrProps%NTwrCD))//'.' )
         ErrStat = ErrID_Fatal
         RETURN
      END IF

   END DO ! I

   !-------------------------------------------------------------------------------------------------
   ! Read section 3, Re vs CD PROPERTIES;
   ! this section has 2 header lines plus NTwrRe rows of data with NTwrCD+1 columns
   !-------------------------------------------------------------------------------------------------

      ! Read in 2 header/comment lines
   CALL ReadCom( UnIn, FilName, 'Re vs CD Properties header 1', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN

   CALL ReadCom( UnIn, FilName, 'Re vs CD Properties header 2', ErrStat, ErrMess )
   IF ( ErrStat /= ErrID_None ) RETURN

   Fmt = '('//TRIM(Int2Lstr(p%TwrProps%NTwrCD+1))//'(2X,ES11.4e2))'

   DO I = 1, p%TwrProps%NTwrRe

      READ( UnIn,*,IOSTAT=ErrStat ) p%TwrProps%TwrRe(I), (p%TwrProps%TwrCD(I,J), J = 1, p%TwrProps%NTwrCD)

      IF ( ErrStat == 0 ) THEN
         IF ( p%Echo ) THEN
            WRITE (p%UnEc,Fmt)  p%TwrProps%TwrRe(I), (p%TwrProps%TwrCD(I,J), J = 1, p%TwrProps%NTwrCD)
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
!> Calculates the axial and tangential induction factor for each annular segment
! and time step (i.e. sets m%Element%A and m%Element%AP)
   SUBROUTINE ELEM_INDUCTIONS( p, m, ErrStat, ErrMess, &
                      PSI, RLOCAL, J, IBlade, VNROTOR2, VT, VNW, &
                      VNB, Initial)

   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   REAL(ReKi),INTENT(IN)      :: PSI
   REAL(ReKi),INTENT(IN)      :: RLOCAL
   REAL(ReKi),                    INTENT(IN   ) :: VNB      ! Normal (relative) velocity of the element 
   REAL(ReKi),INTENT(IN)      :: VNROTOR2
   REAL(ReKi),INTENT(IN)      :: VNW
   REAL(ReKi),                    INTENT(IN   ) :: VT
   INTEGER, INTENT(IN)        :: J
   INTEGER, INTENT(IN)        :: IBlade
   LOGICAL,   INTENT(IN)      :: Initial

   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.
   CHARACTER(ErrMsgLen)                      :: ErrMessLcl          ! Error message returned by called routines.

   ErrStat = ErrID_None
   ErrMess = ""

 !-mlb  Check for being at the center of rotation.
 ! If we are at the center of rotation, the induction equations
 !  are undefined, so let's just USE zeros.

! initialize AxInd and TanInd variables
   m%Element%A (J,IBLADE) = 0.0
   m%Element%AP(J,IBLADE) = 0.0

IF ( RLOCAL < 0.01 )  THEN
    ! Already set to 0
ELSEIF( P%DYNINFL .AND. P%Blade%R * m%Rotor%REVS < 2.0 )  THEN   !ACH 3/10/03 This block deals with dyn. inflow problems at low tip speed
    ! Already set to 0
   m%DYNINIT = .TRUE.    !Re-initialize if we begin using dynamic inflow again
ELSE

 ! Turn wake off when using dynamic inflow and tip speed goes low.  Wake will remain off.
 ! Get induction factor = A using static airfoil coefficients
   IF ( P%WAKE .AND. .NOT. Initial) THEN

      IF ( P%DYNINFL ) THEN
 !       USE dynamic inflow model to find A
         CALL VINDINF( P, m, ErrStatLcl, ErrMessLcl, &
                       J, IBlade, RLOCAL, VNW, VNB, VT, PSI ) !possibly changes A, and AP
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'ELEM_INDUCTIONS' )
            IF (ErrStat >= AbortErrLev) RETURN
      ELSE
 !       USE momentum balance to find A
         CALL VIND( P, m, ErrStatLcl, ErrMessLcl, &
                    J, IBlade, RLOCAL, VNROTOR2, VNW, VNB, VT )  !changes A, and AP
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'ELEM_INDUCTIONS' )
            IF (ErrStat >= AbortErrLev) RETURN
 !       Apply skewed-wake correction, if applicable
         IF( m%SKEW ) CALL VNMOD( P, m, ErrStatLcl, ErrMessLcl,&
                                  J, IBlade, RLOCAL, PSI ) !changes A
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'ELEM_INDUCTIONS' )
            IF (ErrStat >= AbortErrLev) RETURN
      ENDIF
   ENDIF
ENDIF

END SUBROUTINE ELEM_INDUCTIONS

SUBROUTINE ELEMFRC2( p, m, ErrStat, ErrMess, J, IBlade, &
                      DFN, DFT, PMA, Initial, phi )

   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),  INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),    INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER,                   INTENT(  OUT)  :: ErrStat
   CHARACTER(*),              INTENT(  OUT)  :: ErrMess

   REAL(ReKi),                INTENT(  OUT)  :: DFN
   REAL(ReKi),                INTENT(  OUT)  :: DFT
   REAL(ReKi),                INTENT(  OUT)  :: PMA
   INTEGER,                   INTENT(IN)     :: J
   INTEGER,                   INTENT(IN)     :: IBlade
   LOGICAL,                   INTENT(IN)     :: Initial

   ! Local Variables:

   REAL(ReKi)                 :: CDA
   REAL(ReKi)                 :: CLA
   REAL(ReKi)                 :: CMA
   REAL(ReKi)                 :: CPHI
   REAL(ReKi), intent(in)     :: PHI
   REAL(ReKi)                 :: QA
   REAL(ReKi)                 :: ReNum
   REAL(ReKi)                 :: SPHI

   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.
   CHARACTER(ErrMsgLen)                      :: ErrMessLcl          ! Error message returned by called routines.

   ErrStat = ErrID_None
   ErrMess = ""

 ! Get the Reynold's number for the element
 !  Returns Reynold's number x 10^6    !bjj: Reynold's number x 10^-6 ?
ReNum = GetReynolds( SQRT(m%Element%W2(J,IBlade)), P%Blade%C(J), P%Wind%KinVisc )
IF (P%Reynolds) m%AirFoil%MulTabLoc = ReNum

 ! Get lift coefficient from dynamic stall routine if desired
 !  note that the induced velocity was calculated
 !  using the static CL, not the dynamic CL

IF ( P%DSTALL ) THEN
 ! USE Beddoes dynamic stall model
   IF (Initial) THEN ! USE static data on first pass
      CALL BEDINIT ( P, m, ErrStatLcl, ErrMessLcl, &
                     J, IBlade, m%Element%ALPHA(J,IBlade))
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'ELEMFRC' )
            IF (ErrStat >= AbortErrLev) RETURN
      
      CALL CLCD( P, m, ErrStatLcl, ErrMessLcl, &
                 m%Element%ALPHA(J,IBlade), CLA, CDA, CMA, P%AirFoil%NFOIL(J) )
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'ELEMFRC' )
            IF (ErrStat >= AbortErrLev) RETURN
   ELSE
      CALL BeddoesModel( P,  m,  ErrStatLcl, ErrMessLcl, &
                        m%Element%W2(J,IBlade), J, IBlade, m%Element%ALPHA(J,IBlade), CLA, CDA, CMA)
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'ELEMFRC' )
            IF (ErrStat >= AbortErrLev) RETURN
   ENDIF
ELSE
 ! Don't USE dynamic stall model
   CALL CLCD( P,  m,  ErrStatLcl, ErrMessLcl, &
              m%Element%ALPHA(J,IBlade), CLA, CDA, CMA, P%AirFoil%NFOIL(J) )
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'ELEMFRC' )
      IF (ErrStat >= AbortErrLev) RETURN
ENDIF

QA       = 0.5 * P%Wind%RHO * m%Element%W2(J,IBlade) * P%Blade%DR(J) * P%Blade%C(J)
CPHI     = COS( PHI )
SPHI     = SIN( PHI )
DFN      = ( CLA * CPHI + CDA * SPHI ) * QA
DFT      = ( CLA * SPHI - CDA * CPHI ) * QA

IF ( P%PMOMENT ) THEN
   PMA  = CMA * QA * P%Blade%C(J)
ELSE
   PMA  = 0.
   CMA  = 0.
ENDIF


 ! Save values at appropriate station

IF ( IBLADE == 1 ) THEN
   IF ( m%ElOut%ElPrList(J) > 0 )  THEN
      m%ElOut%AAA    ( m%ElOut%ElPrList(J) )    = m%Element%A (J,IBLADE)
      m%ElOut%AAP    ( m%ElOut%ElPrList(J) )    = m%Element%AP(J,IBLADE)
      m%ElOut%ALF    ( m%ElOut%ElPrList(J) )    = m%Element%ALPHA(J,IBlade) * R2D
      m%ElOut%CDD    ( m%ElOut%ElPrList(J) )    = CDA
      m%ElOut%CLL    ( m%ElOut%ElPrList(J) )    = CLA
      m%ElOut%CMM    ( m%ElOut%ElPrList(J) )    = CMA
      m%ElOut%CNN    ( m%ElOut%ElPrList(J) )    = CLA * COS(m%Element%ALPHA(J,IBlade)) + CDA * SIN(m%Element%ALPHA(J,IBlade))
      m%ElOut%CTT    ( m%ElOut%ElPrList(J) )    =-CDA * COS(m%Element%ALPHA(J,IBlade)) + CLA * SIN(m%Element%ALPHA(J,IBlade))
      m%ElOut%DFNSAV ( m%ElOut%ElPrList(J) )    = DFN
      m%ElOut%DFTSAV ( m%ElOut%ElPrList(J) )    = DFT
      m%ElOut%DynPres( m%ElOut%ElPrList(J) )    = 0.5 * P%Wind%RHO * m%Element%W2(J,IBlade)
      m%ElOut%PITSAV ( m%ElOut%ElPrList(J) )    = m%Element%PitNow(J,IBlade) * R2D
      m%ElOut%PMM    ( m%ElOut%ElPrList(J) )    = PMA
      m%ElOut%ReyNum ( m%ElOut%ElPrList(J) )    = ReNum
      m%ElOut%Gamma  ( m%ElOut%ElPrList(J) )    = 0.5 * P%Blade%C(J) * sqrt(m%Element%W2(J,IBlade)) * CLA ! 1/2 c Urel Cl [m^2/s]
   ENDIF

ENDIF

RETURN
END SUBROUTINE ELEMFRC2

!======================================================
   SUBROUTINE VIND( p, m, ErrStat, ErrMess, &
                      J, IBlade, RLOCAL, VNROTOR2, VNW, VNB, VT )
 ! Calculates the axial and tangential induction factor for each annular segment
 ! and time step (i.e. sets m%Element%A and m%Element%AP)
 ! ***************************************************
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),     INTENT(IN   )  :: p            ! Parameters
   TYPE(AD14_MiscVarType),       INTENT(INOUT)  :: m            ! Misc/optimization variables

   REAL(ReKi),                   INTENT(IN   )  :: RLOCAL
   REAL(ReKi),                   INTENT(IN   )  :: VNB
   REAL(ReKi),                   INTENT(IN   )  :: VNROTOR2
   REAL(ReKi),                   INTENT(IN   )  :: VNW
   REAL(ReKi),                   INTENT(IN   )  :: VT ! tangential velocity from relative blade motion and wind, no induction
                                 
   INTEGER,                      INTENT(IN   )  :: J
   INTEGER,                      INTENT(IN   )  :: IBlade
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMess


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
   REAL(ReKi)                    :: SOLFACT
   REAL(ReKi)                    :: VNA
   REAL(ReKi)                    :: VT2_Inv
   REAL(ReKi)                    :: VTA

   INTEGER                       :: ICOUNT
   INTEGER                       :: MAXICOUNT
   INTEGER                       :: Sttus
   INTEGER(IntKi)                :: ErrStatLcl
   character(ErrMsgLen)          :: ErrMessLcl

   ErrStat = ErrID_None
   ErrMess = ""

 ! Allocate and initialize the local array on the first pass
IF ( .NOT. ALLOCATED (m%Element%OLD_A_NS) ) THEN
   ALLOCATE ( m%Element%OLD_A_NS ( P%Element%NELM, P%NumBl) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLD_A_NS array.' )
   m%Element%OLD_A_NS(:,:) = 0.0
ENDIF

 ! Allocate and initialize the local array on the first pass
IF ( .NOT. ALLOCATED (m%Element%OLD_AP_NS) ) THEN
   ALLOCATE ( m%Element%OLD_AP_NS ( P%Element%NELM, P%NumBl) , STAT=Sttus )
   IF ( Sttus /= 0 ) CALL ProgAbort ( ' Error allocating memory for OLD_AP_NS array.' )
   m%Element%OLD_AP_NS(:,:) = 0.0
ENDIF

 ! Set maximum iterations
MAXICOUNT = P%MAXICOUNT

 ! CH--  Alternate convergence criteria
ATOLER2    =  2.0 * P%InducedVel%ATOLER
ATOLERBY10 =  0.1 * P%InducedVel%ATOLER

 ! Bypass calculations for low wind speed, assume no induced velocity.

IF ( VNROTOR2 < 0.1 ) THEN
   m%Element%A(J,IBLADE) = 0.0
   RETURN
ENDIF

 ! SOLFACT is solidity factor divided by 2*VNROTOR2
 ! VT2_Inv is 1./VT**2 to save computation time

IF ( RLOCAL == 0.0 ) THEN   ! Avoid div/0 in FAST2
   SOLFACT = 1.0/VNROTOR2
ELSE
   SOLFACT = P%NumBl * P%Blade%C(J) / ( TWOPI * RLOCAL * VNROTOR2)
ENDIF
VT2_Inv = 1. / ( VT * VT )

 !-mlb  Let's USE the old value of the A from before it was corrected for skew.
AI      = m%Element%OLD_A_NS( J, IBLADE )
DAI1    = 0.05
A2P     = m%Element%OLD_AP_NS( J, IBLADE )
ASTEP   = 0.5
ICOUNT  = 0

 ! Check for extremely high VN and bypass calculations if necessary

IF ( ABS( VNB ) > 100. ) THEN
   m%Element%A( J, IBLADE ) = 0.0
   CALL VINDERR( m, ErrStat, ErrMess, &
                 VNW, VNB, 'VNB', J, IBLADE )
   RETURN   
ELSEIF ( ABS( VT ) > 400. ) THEN
   m%Element%A( J, IBLADE ) = 0.0
   CALL VINDERR( m, ErrStat, ErrMess, &
                 VNW, VT, 'VT', J, IBLADE )
   RETURN
ENDIF

A2 = AI

CALL AXIND ( P, m, ErrStat, ErrMess, &
             VNW, VNB, VNA, VTA, VT, VT2_Inv, VNROTOR2, A2, A2P, &
             J, IBlade, SOLFACT, ALPHA, PHI, CLA, CDA, CMA, RLOCAL )
      IF (ErrStat >= AbortErrLev) RETURN

DAI = A2  - AI

DELAI = ASTEP * DAI

 !CH--  Modification of mlb's proposed change
 ! Must pass two criteria. If we have crossed zero many times
 !  then the first criterion will be easier to meet than the second
 !  because ASTEP will be small (but the second is relaxed to ATOLER2)

DO WHILE ( ABS( DELAI ) > ATOLERBY10 .AND. ABS(DAI) > ATOLER2 )

   ICOUNT = ICOUNT + 1

   A2 = AI

   CALL AXIND ( P, m, ErrStatLcl, ErrMessLcl,       &
                VNW, VNB, VNA, VTA, VT, VT2_Inv, VNROTOR2, A2, A2P, &
                J, IBlade, SOLFACT, ALPHA, PHI, CLA, CDA, CMA, RLOCAL )
      CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,'VIND' )
      IF (ErrStat >= AbortErrLev) RETURN

   DAI = A2  - AI

   DELAI = ASTEP * DAI

 ! Test for convergence, program warning after 1000 iterations

   IF ( ICOUNT > MAXICOUNT ) THEN
      CALL ProgWarn( 'Induction factor calculation did not converge after'//TRIM(Int2LStr(MAXICOUNT))// &
                     ' iterations. AeroDyn will continue using induction factors from previous successful time step.' )
      A2  = m%Element%OLD_A_NS (J,IBLADE)
      A2P = m%Element%OLD_AP_NS(J,IBLADE)
      !CALL SetErrStat(ErrID_Warn,ErrMessLcl,ErrStat,ErrMess,'VIND' )

      EXIT
   ENDIF

 ! Reduce step size after a zero crossing
 !CH--  Put floor under ASTEP to keep it reasonable after many zero crossings

   IF( NINT( SIGN(1.0_ReKi, DAI) ) /= NINT( SIGN(1.0_ReKi, DAI1) ) ) ASTEP = MAX( 1.0E-4_ReKi, 0.5_ReKi*ASTEP )

   AI   = AI + DELAI
   DAI1 = DELAI

END DO

 ! Passed test, we're done
m%Element%A (J,IBLADE) = A2
m%Element%AP(J,IBLADE) = A2P
m%Element%OLD_A_NS  (J,IBLADE) = A2
m%Element%OLD_AP_NS (J,IBLADE) = A2P

RETURN
END SUBROUTINE VIND


 ! ***************************************************
   SUBROUTINE VINDERR( m, ErrStat, ErrMess, &
                      VNW, VX, VID, J, IBLADE )
!   SUBROUTINE VINDERR( VNW, VX, VID, J, IBLADE )
 !  used to write warning messages to the screen
 !  when VN or VT is high.
 ! ***************************************************
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   REAL(ReKi),INTENT(IN)      :: VNW
   REAL(ReKi),INTENT(IN)      :: VX

   INTEGER   ,INTENT(IN)      :: IBLADE
   INTEGER   ,INTENT(IN)      :: J

   CHARACTER(  *),INTENT(IN)  :: VID


   ErrStat = ErrID_None
   ErrMess = ""

 ! Don't write messages if we've already done it 5 times

IF ( m%AFLAGVinderr ) RETURN

m%NERRORS = m%NERRORS + 1

   CALL ProgWarn( ' High '//TRIM(VID)//' velocity encountered during induction factor calculation.' )
   CALL WrScr( '  Blade number '//TRIM(Int2LStr(IBLADE))//', Element number '//TRIM(Int2LStr(J )) )
   CALL WrScr( '  VNW = '       //TRIM(Num2LStr(VNW))//', '//TRIM(VID)//' = '//TRIM(Num2LStr(VX)) )

IF ( m%NERRORS >= 5 ) THEN
   m%AFLAGVinderr = .TRUE.
   CALL ProgWarn( ' Induced velocity warning written 5 times. '//&
                 ' The message will not be repeated, though the condition may persist.' )
ENDIF


RETURN
END SUBROUTINE VINDERR

 ! ******************************************************
   SUBROUTINE  AXIND (P, m, ErrStat, ErrMess,      &
                      VNW, VNB, VNA, VTA, VT, VT2_Inv, VNROTOR2, A2,     &
                      A2P, J, IBlade, SOLFACT, ALPHA, PHI, CLA, CDA, CMA, RLOCAL )
!   SUBROUTINE AXIND ( VNW, VNB, VNA, VTA, VT, VT2_Inv, VNROTOR2, A2, &
!                      A2P, J, SOLFACT, ALPHA, PHI, CLA, CDA, CMA, RLOCAL )
 !  calculates a new axial induction factor from
 !   given values of velocities and geometry.  This routine
 !   is called by vind as part of the iteration process
 ! ******************************************************
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

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
   INTEGER   ,INTENT(IN)      :: IBlade

   ! Local Variables:

   REAL(ReKi)                 :: CH
   REAL(ReKi)                 :: CPhi                                                     ! COS( PHI )
   REAL(ReKi)                 :: SPHI
   REAL(ReKi)                 :: SWRLARG
   REAL(ReKi)                 :: W2
   
   
   ErrStat = ErrID_None
   ErrMess = ""


VNA    = VNW * ( 1. - A2 ) + VNB
VTA    = VT  * ( 1. + A2P )

 ! Get airfoil CL and CD
PHI    = ATAN2( VNA, VTA )
ALPHA  = PHI - m%Element%PitNow(J,IBlade)

CALL MPI2PI ( ALPHA )

CALL CLCD ( P,  m,  ErrStat, ErrMess, &
            ALPHA, CLA, CDA, CMA, P%AirFoil%NFoil(J) )
      IF (ErrStat >= AbortErrLev) RETURN

W2   = VNA * VNA + VTA * VTA
SPHI = VNA/SQRT( W2 )
CPhi = COS( Phi )

 ! Calculate new value of A.  Optionally include normal force due to drag.

CH = W2*SOLFACT*( CLA*CPhi + P%InducedVel%EqAIDmult*CDA*SPhi )


 ! Get the tip loss values for the element (if they change)
IF (p%InducedVel%TLOSS) CALL GetTipLoss ( P, m, ErrStat, ErrMess,      &
                                          J, SPHI, m%TIPLOSS, RLOCAL)

 ! Get the hub loss values for the element (if they change)
IF (p%InducedVel%HLOSS) CALL GetPrandtlLoss ( P%Element%HLCNST(J), SPHI, m%HUBLOSS)

 ! Get the total loss for the element
m%LOSS = m%TIPLOSS * m%HUBLOSS

 ! Check for diverging CH and correct if necessary

IF ( ABS( CH ) > 2. ) CH = SIGN( 2.0_ReKi, CH )

IF ( CH < 0.96*m%LOSS ) THEN
   A2 = 0.5*( 1 - SQRT( 1.0 - CH/m%LOSS ) )
ELSE
   A2 = 0.1432 + SQRT( -0.55106 + .6427*CH/m%LOSS)
ENDIF

 ! Calculate induced swirl (a') if desired.
 !  From C. Ross Harmon's paper on PROPX.

IF ( p%SWIRL ) THEN
   IF ( p%InducedVel%EquilDT )  THEN     ! USE PROP-PC style tangential induction equation with the addition of the drag term.
         ! Because of the singularity that occurs when phi approaches zero,
         ! let's test for small phi and set a' equal to a small, negative number.
      IF ( ( ABS( SPhi ) > 0.01 ) .AND. ( ABS( CPhi ) > 0.01 ) )  THEN
         A2P = SOLFACT*( CLA*SPhi - CDA*CPhi )*( 1.0 + A2P )*VNROTOR2/( 4.0*m%LOSS*SPhi*CPhi )
      ELSEIF ( ABS( SPhi ) > 0.01 )  THEN   ! Tangential velocity near zero, phi near 90 degrees.
         A2P = SOLFACT*( CLA*SPhi - CDA*SIGN( 0.01_ReKi, CPhi ) )*( 1.0 + A2P )*VNROTOR2/( 4.0*m%LOSS*SPhi*SIGN( 0.01_ReKi, CPhi ) )
      ELSE   ! Normal velocity near zero, phi near 0 degrees.
         A2P = SOLFACT*( CLA*SIGN( 0.01_ReKi, SPhi ) - CDA*CPhi )*( 1.0 + A2P )*VNROTOR2/( 4.0*m%LOSS*SIGN( 0.01_ReKi, SPhi )*CPhi )
      ENDIF
   ELSE
      SWRLARG = 1.0 + 4.0*m%LOSS*A2*VNW*VNA*VT2_Inv
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
   FUNCTION GetReynolds( WindSpd, ChordLen, KinVisc )
 !  computes the Reynolds number for the element, divided by 1.0E6
 ! ***************************************************

IMPLICIT                      NONE

   ! Passed Variables:

REAL(ReKi),INTENT(IN)      :: WindSpd
REAL(ReKi),INTENT(IN)      :: ChordLen
REAL(ReKi),INTENT(IN)      :: KinVisc

   ! function definition
REAL(ReKi)     :: GetReynolds

GetReynolds = 1.0E-6 * WindSpd * ChordLen / KinVisc


RETURN
END FUNCTION GetReynolds

 ! ***************************************************
   SUBROUTINE GetTipLoss( P, m, ErrStat, ErrMess,      &
                          J, SPHI, TIPLOSS, RLOCAL )
!   SUBROUTINE GetTipLoss( J, SPHI, TIPLOSS, RLOCAL )
 !  computes the tip loss constant for element J
 !  TIPLOSS is returned to AXIND
 ! Uses the Prandtl tip loss model with a correction
 !  from Georgia Tech (2002 ASME Wind Energy Symposium)
 ! ***************************************************
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   ! Passed Variables:
   REAL(ReKi), INTENT(IN)     :: SPHI
   REAL(ReKi), INTENT(IN)     :: RLOCAL
   REAL(ReKi), INTENT(OUT)    :: TIPLOSS

   INTEGER   , INTENT(IN)     :: J


   ! Local Variables:

   REAL(ReKi)                 :: Dist2pt7 = 0.7 ! current element distance to r/R = 0.7
   REAL(ReKi)                 :: OLDDist7       ! previous element distance to r/R = 0.7
   REAL(ReKi)                 :: percentR

   INTEGER                    :: Jpt7 = 0 ! The element closest to r/R = 0.7


   ErrStat = ErrID_None
   ErrMess = ""

 ! Calculate PRANDTL tip loss model
CALL GetPrandtlLoss( P%Element%TLCNST(J), SPHI, TIPLOSS )


 ! Apply Georgia Tech correction to Prandtl model if activated
IF (p%InducedVel%GTECH) THEN
   percentR = RLOCAL/P%Blade%R

 ! Search for the element closest to r/R = 0.7
   IF (m%FirstPassGTL) THEN
    ! If the current element is closer than the previous, update values
      IF ( ABS(percentR - 0.7) < Dist2pt7 ) THEN
         OLDDist7 = Dist2pt7
         Dist2pt7 = ABS(percentR - 0.7)
         Jpt7 = J
         m%TLpt7 = TIPLOSS
      ENDIF
      IF (J == P%Element%NELM) THEN ! We're done after one pass through the blades
         m%FirstPassGTL = .FALSE.
      ELSE
         RETURN ! Don't do the correction until we calculate the correct TLpt7
      ENDIF
   ENDIF

   IF ( J == Jpt7 ) m%TLpt7 = TIPLOSS ! Update the value of TLpt7 at the proper element

 ! Do the actual Georgia Tech correction to the Prandtl model
   IF (percentR >= 0.7) THEN
      TIPLOSS = (TIPLOSS**0.85 + 0.5 ) / 2.0
   ELSE
      TIPLOSS = 1.0 - percentR*(1.0 - m%TLpt7)/0.7
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
SUBROUTINE GetTwrInfluence ( P, m, ErrStat, ErrMess,      &
                             VX, VY, InputPosition)
!SUBROUTINE GetTwrInfluence (VX, VY, InputPosition)
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

   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

      ! Passed Variables:

   REAL(ReKi), INTENT(INOUT)  :: VX                   ! on input, U-velocity without tower effect; on output, U-velocity including tower effect
   REAL(ReKi), INTENT(INOUT)  :: VY                   ! on input, V-velocity without tower effect; on output, V-velocity including tower effect
   REAL(ReKi), INTENT(IN)     :: InputPosition(3)     !velocities in global coordinates with tower effect

      ! Local Variables:

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

   
   ErrStat = ErrID_None
   ErrMess = ""

   !-------------------------------------------------------------------------------------------------
   ! This subroutine is only valid for TwrPotent and TwrShadow features
   !-------------------------------------------------------------------------------------------------
   IF (.NOT. p%TwrProps%TwrPotent .AND. .NOT. p%TwrProps%TwrShadow) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Initialize some variables
   !-------------------------------------------------------------------------------------------------
   ZGrnd   = InputPosition(3) - p%Rotor%HH            ! distance between position and hub !BJJ: this should really be the tower height (position), not HH
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

   CALL GetTwrSectProp (P, m, ErrStat, ErrMess,      &
                        InputPosition(:), V_total, TwrRad, TwrCD_Station)     ! Get the tower properties for the current element location

   Distance = SQRT ( InputPosition(1)**2 + InputPosition(2)**2 ) / TwrRad     ! normalized distance to tower center

      ! Check for tower strike
   IF ( Distance < 1.0 ) THEN    ! potentially inside the tower            !bjj: only if we're not ABOVE the tower, though....

      IF (ZGrnd < 0.0) THEN      !bjj added this condition.... check that it's correct
!bjj: what does this exactly mean? "temporarily disabled?"
         IF( .NOT. m%AFLAGTwrInflu) THEN
            CALL ProgWarn( ' Tower model temporarily disabled due to possible tower strike.'// &
                           ' This message will not be repeated though the condition may persist.' )
               !write a blank line (so FAST doesn't write over it)
            CALL WrScr( ' ' )
            m%AFLAGTwrInflu = .TRUE.
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

   IF ( P%TwrProps%TwrPotent ) THEN

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


      Xtemp  = Xwind + P%TwrProps%Tower_Wake_Constant !PJM fixed this error 3/30/06
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

   IF ( P%TwrProps%TwrShadow ) THEN

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
SUBROUTINE GetTwrSectProp ( P, m, ErrStat, ErrMess,      &
                            InputPosition, VelHor, TwrElRad, TwrElCD )
!SUBROUTINE GetTwrSectProp (InputPosition, VelHor, TwrElRad, TwrElCD)
!  Returns the tower radius and CD for the vertical location
!   of the element currently being evaluated for tower influence.
!====================================================================================================

   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess


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


   ErrStat = ErrID_None
   ErrMess = ""

   !-------------------------------------------------------------------------------------------------
   ! Get the tower radius, TwrElRad, by interpolating into the TwrWid(:) array
   !-------------------------------------------------------------------------------------------------

   TwrElHt  = InputPosition(3) / P%Rotor%HH                                  !!!!BJJ!!!! HH????   !MLB: It's the hub height.  It really should be the tower height.
   TwrElRad = 0.5*InterpBin( TwrElHt, p%TwrProps%TwrHtFr, p%TwrProps%TwrWid, N2, p%TwrProps%NTwrHt )

   !-------------------------------------------------------------------------------------------------
   ! Get the section CD, TwrElCD, by interpolating into the TwrCD(:,:) array
   !-------------------------------------------------------------------------------------------------

   IF ( p%TwrProps%NTwrRe == 1 ) THEN                                  ! There is only one Re row

      IF ( p%TwrProps%NTwrCD == 1 ) THEN                               ! There is only one CD column
         TwrElCD = p%TwrProps%TwrCD(1,1)
      ELSE IF ( p%TwrProps%NTwrHt == 1 ) THEN                          ! There is more than one column of CD, but only one used
         TwrElCD = p%TwrProps%TwrCD(1,p%TwrProps%NTwrCDCol(1))
      ELSE                                                  ! Interpolate;  this will be the same Indx as before...
         TwrElCD = InterpStp( TwrElHt, p%TwrProps%TwrHtFr, p%TwrProps%TwrCD(1,:), N2, p%TwrProps%NTwrHt )
      END IF

   ELSE                                                     ! There are multiple Re rows

      TwrElRe = GetReynolds( VelHor, 2.0_ReKi*TwrElRad, P%Wind%KinVisc )

      IF ( p%TwrProps%NTwrCD == 1 ) THEN                               ! There is only one CD column
         TwrElCD = InterpBin( TwrElRe, p%TwrProps%TwrRe, p%TwrProps%TwrCD(:,1), N1, p%TwrProps%NTwrRe )
      ELSE IF ( p%TwrProps%NTwrHt == 1 ) THEN                          ! Interpolate over Re only
         TwrElCD = InterpBin( TwrElRe, p%TwrProps%TwrRe, p%TwrProps%TwrCD(:,p%TwrProps%NTwrCDCol(1)), N1, p%TwrProps%NTwrRe )
      ELSE                                                  ! A 2-D interpolation is needed
         CALL LocateBin( TwrElRe, p%TwrProps%TwrRe, N1, p%TwrProps%NTwrRe )

            ! Let's use nearest-neighbor extrapolation with bi-linear interpolation:

         N1   = MIN( MAX( N1, 1 ), p%TwrProps%NTwrRe-1 )
         N1P1 = N1+1

         P1   = MIN( MAX( (TwrElRe - p%TwrProps%TwrRe(N1))   / (p%TwrProps%TwrRe(N1P1)   - p%TwrProps%TwrRe(N1))  , REAL(0.0, ReKi) ), REAL(1.0, ReKi) )

         N2P1 = N2 + 1
         P2   = MIN( MAX( (TwrElHt - p%TwrProps%TwrHtFr(N2)) / (p%TwrProps%TwrHtFr(N2P1) - p%TwrProps%TwrHtFr(N2)), REAL(0.0, ReKi) ), REAL(1.0, ReKi) )


         TwrElCD1 = p%TwrProps%TwrCD(N1,N2  ) + P1 * ( p%TwrProps%TwrCD(N1P1,N2  ) - p%TwrProps%TwrCD(N1,N2  ) )
         TwrElCD2 = p%TwrProps%TwrCD(N1,N2P1) + P1 * ( p%TwrProps%TwrCD(N1P1,N2P1) - p%TwrProps%TwrCD(N1,N2P1) )


         TwrElCD = TwrElCD1 + P2 * ( TwrElCD2 - TwrElCD1 )
      END IF

   END IF


RETURN
END SUBROUTINE GetTwrSectProp


!====================================================================================================
FUNCTION AD_WindVelocityWithDisturbance(  Time, u, p, x, xd, z, m, y, ErrStat, ErrMsg,      &
                                          InputPosition, InputVelocity )
!                                          InputPosition, TShadC1, TShadC2, PJM_Version  )
!FUNCTION AD_WindVelocityWithDisturbance( InputPosition, TShadC1, TShadC2, PJM_Version, LeStat  )
!  This function computes the (dimensional) wind velocity components at the location InputPosition
!  in the inertial frame of reference, including any tower shadow defecit.
!  ** Formerly SUBROUTINE VEL and later SUBROUTINE VWrel2G( VNRotor2, At_Hub ) **
!----------------------------------------------------------------------------------------------------


   IMPLICIT                         NONE

      ! Passed Variables:

      REAL(DbKi),                     INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD14_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD14_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD14_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD14_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(AD14_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
      TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat  ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg   ! Error message if ErrStat /= ErrID_None

      REAL(ReKi),INTENT(IN)            :: InputPosition(3)        ! 
      REAL(ReKi),INTENT(IN)            :: InputVelocity(3)        ! undisturbed velocity

   !   REAL(ReKi),INTENT(IN)            :: TShadC1
!   REAL(ReKi),INTENT(IN)            :: TShadC2
!   LOGICAL,INTENT(IN)               :: PJM_Version

      ! function definition

   REAL(ReKi)                      :: AD_WindVelocityWithDisturbance(3)

      ! Local variables

   REAL(ReKi)                       :: angle    ! absolute difference between theta and phi
   REAL(ReKi)                       :: dist     ! distance from blade element to wake centerline
   REAL(ReKi)                       :: phi      ! angle between x-axis and instantaneous horizontal wind direction
   REAL(ReKi)                       :: RADIUS   ! horizontal distance from tower to blade element ** BJJ NOTE: in actuality, it's the distance from the undeflected tower centerline, not the actual tower
   REAL(ReKi)                       :: ROOTR    ! SQRT(radius)
   REAL(ReKi)                       :: SHADOW
   REAL(ReKi)                       :: TEMP
   REAL(ReKi)                       :: theta    ! Angle between x-axis and line from tower to blade element
   REAL(ReKi)                       :: width    ! half width of the wake after accounting for wake expansion proportional to square root of RADIUS

!   INTEGER                          :: Sttus

!   INTEGER                          :: TmpErrStat
!   CHARACTER(ErrMsgLen)             :: TmpErrMsg


   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Get the undisturbed velocity

   AD_WindVelocityWithDisturbance(:) = InputVelocity

         ! Add the tower influence to the undisturbed velocity.

   IF ( p%TwrProps%PJM_Version ) THEN

      CALL GetTwrInfluence ( P, m, ErrStat, ErrMsg, &
                             AD_WindVelocityWithDisturbance(1), AD_WindVelocityWithDisturbance(2), InputPosition(:) )

   ELSE !Old tower version

         ! Apply tower shadow if the blade element is in the wake

      IF ( p%TwrProps%TShadC2 > 0.0 ) THEN ! Perform calculations only if the wake strength is positive

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
               width  = p%TwrProps%TShadC1 * RootR                                                               ! half width of the wake after accounting for wake expansion proportional to square root of RADIUS

               IF ( width > 0 ) THEN   ! Skip cases with zero width or radius so we don't divide by zero

                  dist = radius * SIN( angle )                                                        ! distance from blade element to wake centerline
                  IF ( InputPosition(3) > p%Rotor%HH )  THEN                                                  ! Apply shadow in arc above hub to maintain a continuous deficit function.  Somewhat of a nacelle deficit.
                     dist = SQRT( dist**2 + (InputPosition(3)-p%Rotor%HH)**2 )                                !bjj: I think this should use hub position, not HH
                  END IF

                  IF ( width > dist ) THEN   ! There is velocity deficit in the wake only.
                     temp   = COS ( PiBy2 * dist/width )
                     shadow = p%TwrProps%TShadC2/RootR * temp * temp

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
SUBROUTINE DiskVel ( Time, P, m, AvgInfVel, ErrStat, ErrMess )
!   SUBROUTINE DiskVel
 !  calculates the mean velocities relative to the rotor disk
 !  calls routine to get wind velocity at a specified location
 !
 !  Updated on 08/12/97 xyz-direction changed
 !  Combined VELD and GETSKEW 04/24/01
 !  Updated 12/1/09 to use new inflow module; WindInf_ADhack_diskVel MUST be replaced!
 ! ********************************************
   IMPLICIT                      NONE
      ! Passed Variables:
   REAL(DbKi), INTENT(IN) :: Time
   TYPE(AD14_ParameterType),     INTENT(IN)     :: p            ! Parameters
   TYPE(AD14_MiscVarType),       INTENT(INOUT)  :: m            ! Misc/optimization variables
   REAL(ReKi),                   INTENT(IN)     :: AvgInfVel(3) !some sort of averaged wind speed (currently depends on wind file type)
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess


   REAL(ReKi)                 :: Vinplane
   REAL(ReKi)                 :: VXY

   ErrStat = ErrID_None
   ErrMess = ""
   
!Position = (/0.0_ReKi, 0.0_ReKi, P%Rotor%HH /)
!AvgInfVel(:) = WindInf_ADhack_diskVel( REAL(Time, ReKi), Position, ErrStat )


VXY   = AvgInfVel(1) * m%Rotor%CYaw - AvgInfVel(2) * m%Rotor%SYaw

 ! Mean velocities in rotor disk coord. Includes yaw rate.
 ! X = Normal to plane, parallel to axis of rotation DOWNWIND
 ! Y = Inplane, horizontal to left (looking downwind)
 ! Z = Inplane, vertical up

m%Wind%VROTORX = VXY * m%Rotor%CTILT + AvgInfVel(3) * m%Rotor%STILT

m%Wind%VROTORY = AvgInfVel(1) * m%Rotor%SYaw + AvgInfVel(2) * m%Rotor%CYaw + m%Rotor%YAWVEL

m%Wind%VROTORZ = -1.* VXY * m%Rotor%STILT + AvgInfVel(3) * m%Rotor%CTILT

 ! Skewed wake correction not needed for GDW
IF (.NOT. P%DYNINFL) THEN
 !  Set SKEW flag and assign values to related parameters
 !   used in the skewed wake correction.

 ! Vinplane is the resultant in-plane velocity
   Vinplane = SQRT( m%Wind%VROTORY * m%Wind%VROTORY + m%Wind%VROTORZ * m%Wind%VROTORZ )

 ! SKEW is TRUE if there is a cross flow, FALSE otherwise.
   IF ( Vinplane >= 1.0E-3 ) THEN
      m%SKEW   = .TRUE.
      m%Wind%SDEL   =  m%Wind%VROTORY/Vinplane
      m%Wind%CDEL   = -m%Wind%VROTORZ/Vinplane
      m%Wind%ANGFLW = ATAN2( ABS( m%Wind%VROTORX - m%Rotor%AVGINFL ), Vinplane )
   ELSE
      m%SKEW   = .FALSE.
   ENDIF

ENDIF



RETURN
END SUBROUTINE DiskVel
!=======================================================================
SUBROUTINE TwrAeroLoads ( p, Node, NodeDCMGbl, NodeVelGbl, NodeWindVelGbl, NodeFrcGbl )


      ! This routine calcualtes the aeroynamic loads of all tower nodes above the mean sea level.
      ! It doesn't worry about whether or not a node is below water.  The aero loads will be far less than the hydro loads.

   IMPLICIT                                     NONE


      ! Arguments:

   REAL(R8Ki), INTENT(IN )                   :: NodeDCMGbl    (3,3)           ! The direction-cosine matrix used to transform from the global system to the node system.
   REAL(ReKi), INTENT(OUT)                   :: NodeFrcGbl    (3)             ! The forces per unit length at the current tower element.
   REAL(ReKi), INTENT(IN )                   :: NodeVelGbl    (3)             ! The 3 components of the translational velocity at the node in the global system.
   REAL(ReKi), INTENT(IN )                   :: NodeWindVelGbl(3)             ! The 3 components of the wind at the node in the global system.

   INTEGER,    INTENT(IN )                   :: Node                          ! Tower node index.

   TYPE(AD14_ParameterType), INTENT(IN)      :: p                             ! The AeroDyn parameters.


      ! Local variables.

   REAL(ReKi)                                :: NodeFrcLcl    (3)             ! The forces per unit length on the tower node in the local system.
!   REAL(ReKi)                                :: NodeLocTwr    (3)             ! The location of the node in the tower coordinate system.  This is used to get the tower section properties.
   REAL(ReKi)                                :: NodeVelRelGbl (3)             ! The relative wind velocity in the global reference frame..
   REAL(ReKi)                                :: NodeVelRelLcl (3)             ! The relative wind velocity in the local reference frame..
   REAL(ReKi)                                :: RelNmlWndSpd                  ! The relative wind speed normal to the tower axis.  sqrt(u^2+v^2)
!   REAL(ReKi)                                :: TwrFrcLcl     (3)             ! The drag coefficient corresponding to the computed Reynolds Number.
   REAL(ReKi)                                :: TwrNodeCd                     ! The drag coefficient corresponding to the computed Reynolds Number.
   REAL(ReKi)                                :: TwrNodeRe                     ! The Reynolds Number computed using the wind speed normal to the tower axis.
!   REAL(ReKi)                                :: WndDirLcl                     ! The wind direction relative to the local node coordinate system using only the u and v components.

   INTEGER(IntKi)                            :: IndLo                         ! The index pointing to the lower of the two points bounding an interpolated value.



      ! Calculate the relative local wind velocity.

   NodeVelRelGbl(:) = NodeWindVelGbl(:) - NodeVelGbl(:)


      ! Transform the relative velocity into the node coordinate system.

   NodeVelRelLcl(:) = MATMUL( NodeDCMGbl, NodeVelRelGbl )


      ! Compute the relative normal wind speed and its square.

   RelNmlWndSpd  = SQRT( DOT_PRODUCT( NodeVelRelLcl(1:2), NodeVelRelLcl(1:2) ) )


      ! Compute the Reynolds Number.  Because interpolation is expensive, we will compute the magnitude of the drag and resolve it into components.

   TwrNodeRe = GetReynolds( RelNmlWndSpd, p%TwrProps%TwrNodeWidth(Node), p%Wind%KinVisc )


      ! Get the local value of the drag coefficient.

   IndLo = 1

!MLB: Why have two different calls?  Can't the second accommodate the first?  Is the first method more efficient than the second method if there is only one Cd column?
!MLB: This logic was stolen from AeroSubs.f90\GetTwrSectProp().

   IF ( p%TwrProps%NTwrCD == 1 ) THEN                               ! There is only one CD column
      TwrNodeCd = InterpBin( TwrNodeRe, p%TwrProps%TwrRe, p%TwrProps%TwrCD(:,1), IndLo, p%TwrProps%NTwrRe )
   ELSE IF ( p%TwrProps%NTwrHt == 1 ) THEN                          ! Interpolate over Re only
      TwrNodeCd = InterpBin( TwrNodeRe, p%TwrProps%TwrRe, p%TwrProps%TwrCD(:,p%TwrProps%NTwrCDCol(1)), IndLo, p%TwrProps%NTwrRe )  !MLB Why have two different calls?  Can't the second accommodate the first?
   END IF


      ! Compute the forces on the tower in the local system.

   NodeFrcLcl(1) = 0.5*TwrNodeCd*p%Wind%Rho*p%TwrProps%TwrNodeWidth(Node)*RelNmlWndSpd*NodeVelRelLcl(1)
   NodeFrcLcl(2) = 0.5*TwrNodeCd*p%Wind%Rho*p%TwrProps%TwrNodeWidth(Node)*RelNmlWndSpd*NodeVelRelLcl(2)
   NodeFrcLcl(3) = 0.0


      ! Convert the force to global coordinates.

   NodeFrcGbl = MATMUL( TRANSPOSE( NodeDCMGbl ), NodeFrcLcl )


! Temporarily set the returned force to zero until we figure why the forces are not being applied correctly.
!NodeFrcGbl(:) = 0.0


   RETURN

END SUBROUTINE TwrAeroLoads
!=======================================================================


 ! ****************************************
   SUBROUTINE VNMOD( P, m, ErrStat, ErrMess, &
                     J, IBlade, RLOCAL, PSI )
!   SUBROUTINE VNMOD( J, IBlade, RLOCAL, PSI )
 !  applies the skewed wake correction
 !   to the axial induction factor A.
 ! ****************************************
!USE                           Blade
!USE                           Element
!USE                           Wind
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   ! Passed Variables:

   REAL(ReKi),INTENT(IN)      :: PSI
   REAL(ReKi),INTENT(IN)      :: RLOCAL

   INTEGER,   INTENT(IN)      :: J
   INTEGER,   INTENT(IN)      :: IBlade

   ! Local Variables:

   REAL(ReKi)                 :: BB
   REAL(ReKi)                 :: SANG


   ErrStat = ErrID_None
   ErrMess = ""

SANG = SIN( m%Wind%ANGFLW )
BB   = 0.7363 * SQRT( ( 1. - SANG )/(1. + SANG) )

m%Element%A(J,IBLADE) = m%Element%A(J,IBLADE) * ( 1. + 2. * RLOCAL/P%Blade%R * BB *  &
             ( m%Wind%SDEL * SIN( PSI ) + m%Wind%CDEL * COS( PSI ) )  )



RETURN
END SUBROUTINE VNMOD


 ! **********************************************************
   SUBROUTINE BEDINIT( P, m, ErrStat, ErrMess, &
                       J, IBlade, ALPHA )
!   SUBROUTINE BEDINIT( J, IBlade, ALPHA )
 !  calculates initial values of Beddoes 'f' arrays
 ! **********************************************************
!USE                           Airfoil
!USE                           Beddoes
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

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
   REAL(ReKi)                 :: P0
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

   ErrStat = ErrID_None
   ErrMess = ""

m%Beddoes%ANE(J,IBLADE) = ALPHA
m%Beddoes%AFE(J,IBLADE) = ALPHA

I = P%AirFoil%NFOIL(J)

IF ( P%AirFoil%NTables(I) > 1 ) THEN
   m%AirFoil%MulTabLoc = MIN( MAX( m%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1) ), P%AirFoil%MulTabMet(I,P%AirFoil%NTables(I)) )
   CALL LocateBin( m%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1:P%AirFoil%NTables(I)), N, P%AirFoil%NTables(I) )

   IF (N == 0 ) THEN
      CNA1 = m%Beddoes%CNA(I,1)
      AOL1 = m%Beddoes%AOL(I,1)
   ELSE IF( N == P%AirFoil%NTables(I) ) THEN
      CNA1 = m%Beddoes%CNA(I,N)
      AOL1 = m%Beddoes%AOL(I,N)
   ELSE
      NP1   = N+1
      P0    = (m%AirFoil%MulTabLoc-P%AirFoil%MulTabMet(I,N))/(P%AirFoil%MulTabMet(I,NP1)-P%AirFoil%MulTabMet(I,N))
      CNA1  = m%Beddoes%CNA(I,N) + P0 * ( m%Beddoes%CNA(I,NP1) - m%Beddoes%CNA(I,N) )
      AOL1  = m%Beddoes%AOL(I,N) + P0 * ( m%Beddoes%AOL(I,NP1) - m%Beddoes%AOL(I,N) )
   END IF
ELSE
   CNA1  = m%Beddoes%CNA(I,1)
   AOL1  = m%Beddoes%AOL(I,1)
ENDIF

m%Beddoes%CNPOT(J,IBLADE) = CNA1 * (ALPHA - AOL1)
m%Beddoes%CNP(J,IBLADE)   = m%Beddoes%CNPOT(J,IBLADE)


ALPHA = MIN( MAX( ALPHA, m%AirFoil%AL(I,1) ), m%AirFoil%AL(I,P%AirFoil%NLIFT(I)) )
CALL LocateBin( ALPHA, m%AirFoil%AL(I,1:P%AirFoil%NLIFT(I)), I1, P%AirFoil%NLIFT(I) )

IF ( I1 == 0 ) THEN
   I1   = 1
   I1P1 = 2
   P1   = 0.0
ELSEIF ( I1 == P%AirFoil%NLIFT(I) ) THEN
   I1P1 = I1
   I1   = I1 - 1
   P1   = 1.0
ELSE
   I1P1 = I1 + 1
   !bjj: check for division by zero?
   P1   = (  m%AirFoil%AL(I,I1) - ALPHA ) / (  m%AirFoil%AL(I,I1) -  m%AirFoil%AL(I,I1P1) )
ENDIF


IF ( P%AirFoil%NTables(I) > 1 ) THEN

 ! Locate the multiple airfoil table position in the table


   m%AirFoil%MulTabLoc = MIN( MAX( m%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1) ), P%AirFoil%MulTabMet(I,P%AirFoil%NTables(I)) )
   CALL LocateBin( m%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1:P%AirFoil%NTables(I)), I2, P%AirFoil%NTables(I) )

   IF ( I2 == 0 ) THEN
      I2P1 = 2
      I2   = 1
      P2   = 0.0
   ELSE IF( I2 == P%AirFoil%NTables(I) ) THEN
      I2P1 = I2
      I2   = I2 - 1
      P2   = 1.0
   ELSE
      I2P1 = I2 + 1
      P2 = (m%AirFoil%MulTabLoc-P%AirFoil%MulTabMet(I,I2))/(P%AirFoil%MulTabMet(I,I2P1)-P%AirFoil%MulTabMet(I,I2))
   ENDIF

 ! Interpolate the F-table values

   FSPB  = m%Beddoes%FTB( I,I1,I2P1) - (m%Beddoes%FTB( I,I1,I2P1) - m%Beddoes%FTB( I,I1P1,I2P1))*P1
   FSPCB = m%Beddoes%FTBC(I,I1,I2P1) - (m%Beddoes%FTBC(I,I1,I2P1) - m%Beddoes%FTBC(I,I1P1,I2P1))*P1
   FSPA  = m%Beddoes%FTB( I,I1,I2  ) - (m%Beddoes%FTB( I,I1,I2  ) - m%Beddoes%FTB( I,I1P1,I2  ))*P1
   FSPCA = m%Beddoes%FTBC(I,I1,I2  ) - (m%Beddoes%FTBC(I,I1,I2  ) - m%Beddoes%FTBC(I,I1P1,I2  ))*P1

   m%Beddoes%FSP( J,IBLADE) = FSPA  + P2 * ( FSPB - FSPA )
   m%Beddoes%FSPC(J,IBLADE) = FSPCA + P2 * ( FSPCB - FSPCA )

ELSE

   m%Beddoes%FSP( J,IBLADE) = m%Beddoes%FTB( I,I1,1) - ( m%Beddoes%FTB( I,I1,1) - m%Beddoes%FTB( I,I1P1,1) )*P1
   m%Beddoes%FSPC(J,IBLADE) = m%Beddoes%FTBC(I,I1,1) - ( m%Beddoes%FTBC(I,I1,1) - m%Beddoes%FTBC(I,I1P1,1) )*P1

ENDIF

IF ( ABS( m%Beddoes%AFE(J,IBLADE) - AOL1 ) < 1.E-10 ) THEN

   m%Beddoes%FSP(J,IBLADE)  = 1.0
   m%Beddoes%FSPC(J,IBLADE) = 1.0

ELSE

   TEMP = 2.*SQRT(ABS(m%Beddoes%FSP(J,IBLADE)/(m%Beddoes%AFE(J,IBLADE)-AOL1)))-1.
   m%Beddoes%FSP(J,IBLADE) = TEMP * TEMP * SIGN ( 1.0_ReKi, TEMP )
   IF ( m%Beddoes%FSP(J,IBLADE) >  1.0 ) m%Beddoes%FSP(J,IBLADE) =  1.0
   IF ( m%Beddoes%FSP(J,IBLADE) < -1.0 ) m%Beddoes%FSP(J,IBLADE) = -1.0

   IF ( ABS( m%Beddoes%AFE(J,IBLADE) ) < 1.E-10 ) THEN
      m%Beddoes%FSPC(J,IBLADE) = 1.0
   ELSE
      TEMP = m%Beddoes%FSPC(J,IBLADE)/((m%Beddoes%AFE(J,IBLADE)-AOL1)*m%Beddoes%AFE(J,IBLADE))
      m%Beddoes%FSPC(J,IBLADE) = TEMP * TEMP * SIGN ( 1.0_ReKi, TEMP )
      IF ( m%Beddoes%FSPC(J,IBLADE) >  1.0 ) m%Beddoes%FSPC(J,IBLADE) =  1.0
      IF ( m%Beddoes%FSPC(J,IBLADE) < -1.0 ) m%Beddoes%FSPC(J,IBLADE) = -1.0
   ENDIF

ENDIF

SRFP = SQRT( ABS( m%Beddoes%FSP(J,IBLADE) ) ) * SIGN( 1.0_ReKi, m%Beddoes%FSP(J,IBLADE) ) + 1.
m%Beddoes%FK   = 0.25 * SRFP * SRFP
m%Beddoes%CVN(J,IBLADE) = m%Beddoes%CNPOT(J,IBLADE) * ( 1. - m%Beddoes%FK )

RETURN
END SUBROUTINE BEDINIT


 ! *****************************************************
   SUBROUTINE BedUpdate( m )
 !  Update old Beddoes parameters at new time step
 ! *****************************************************
!USE            Beddoes
   IMPLICIT                      NONE
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables

m%Beddoes%ANE1    = m%Beddoes%ANE
m%Beddoes%ADOT1   = m%Beddoes%ADOT
m%Beddoes%OLDXN   = m%Beddoes%XN
m%Beddoes%OLDYN   = m%Beddoes%YN
m%Beddoes%CNPOT1  = m%Beddoes%CNPOT
m%Beddoes%OLDDPP  = m%Beddoes%DPP
m%Beddoes%FSP1    = m%Beddoes%FSP
m%Beddoes%FSPC1   = m%Beddoes%FSPC
m%Beddoes%OLDTAU  = m%Beddoes%TAU
m%Beddoes%OLDDF   = m%Beddoes%DF
m%Beddoes%OLDDFC  = m%Beddoes%DFC
m%Beddoes%OLDDN   = m%Beddoes%DN
m%Beddoes%OLDCNV  = m%Beddoes%CNV
m%Beddoes%CVN1    = m%Beddoes%CVN
m%Beddoes%CNP1    = m%Beddoes%CNP
m%Beddoes%CNPD1   = m%Beddoes%CNPD
m%Beddoes%OLDSEP  = m%Beddoes%BEDSEP
m%Beddoes%QX1     = m%Beddoes%QX
m%Beddoes%OLDDQ   = m%Beddoes%DQ
m%Beddoes%AFE1    = m%Beddoes%AFE
m%Beddoes%DQP1    = m%Beddoes%DQP
m%Beddoes%DFAFE1  = m%Beddoes%DFAFE


RETURN
END SUBROUTINE BedUpdate


 ! *****************************************************
   SUBROUTINE BEDDAT ( P, x, xd, z, m, y, ErrStat, ErrMess )
 !  USED TO INPUT PARAMETERS FOR THE
 !   Beddoes DYNAMIC STALL MODEL
 ! *****************************************************
!USE                           Airfoil
!USE                           Beddoes
!USE                           Switch
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(AD14_ContinuousStateType), INTENT(INOUT)  :: x           ! Initial continuous states
   TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(AD14_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   ! Local Variables:

   REAL(ReKi)                 :: ETA
   REAL(ReKi)                 :: CA
   REAL(ReKi)                 :: SA

   INTEGER                    :: I
   INTEGER                    :: J
   INTEGER                    :: K


   ErrStat = ErrID_None
   ErrMess = ""

 ! Empirical constants for the Beddoes model

 ! TVL   = Non-dimensional time of transit for the
 !          vortex moving across the airfoil surface
 ! TP    = Time constant for pressure lag
 ! TV    = Time constant for strength of shed vortex
 ! TF    = Time constant applied to location of
 !          the separation point
 ! AS    = Speed of sound for Mach number calculation

P%Beddoes%TVL   = 11.0
P%Beddoes%TP    = 1.7
P%Beddoes%TV    = 6.0
P%Beddoes%TF    = 3.0
ETA   = .99  !bjj: this doesn't seem to be used for anything....

IF ( P%SIUNIT ) THEN
 ! SI UNITS--m/sec
   P%Beddoes%AS = 335.
ELSE
 ! ENGLISH UNITS--ft/sec
   P%Beddoes%AS = 1100.
ENDIF

 ! Generate table of F values from airfoil data table

DO J =1,P%AirFoil%NUMFOIL
   DO K =1,P%AirFoil%NTables(J)
      DO I = 1, P%AirFoil%NLIFT(J)

         CA = COS( m%AirFoil%AL(J,I) )
         SA = SIN( m%AirFoil%AL(J,I) )
         m%Beddoes%CN = m%AirFoil%CL(J,I,K) * CA + ( m%AirFoil%CD(J,I,K) - m%Beddoes%CDO(J,K) ) * SA
         m%Beddoes%CC = m%AirFoil%CL(J,I,K) * SA - ( m%AirFoil%CD(J,I,K) - m%Beddoes%CDO(J,K) ) * CA

         IF ( ABS( m%Beddoes%CNA(J,K) ) .GT. 1.E-6 ) THEN
            m%Beddoes%FTB(J,I,K)  = m%Beddoes%CN / m%Beddoes%CNA(J,K)
            m%Beddoes%FTBC(J,I,K) = m%Beddoes%CC / m%Beddoes%CNA(J,K)
         ELSE
            m%Beddoes%FTB(J,I,K)  = 1.0
            m%Beddoes%FTBC(J,I,K) = 1.0
         ENDIF

      END DO !I
   END DO !K
END DO !J

m%Beddoes%VOR   = .FALSE.
m%Beddoes%SHIFT = .FALSE.

m%Beddoes%BEDSEP = .FALSE.
m%Beddoes%ANE1   = 0.
m%Beddoes%OLDCNV = 0.
m%Beddoes%CVN1   = 0.
m%Beddoes%CNPOT1 = 0.
m%Beddoes%CNP1   = 0.
m%Beddoes%CNPD1  = 0.
m%Beddoes%OLDDF  = 0.
m%Beddoes%OLDDFC = 0.
m%Beddoes%OLDDPP = 0.
m%Beddoes%FSP1   = 0.
m%Beddoes%FSPC1  = 0.
m%Beddoes%TAU    = 0.
m%Beddoes%OLDTAU = 0.
m%Beddoes%OLDXN  = 0.
m%Beddoes%OLDYN  = 0.

RETURN
END SUBROUTINE BEDDAT


 ! ******************************************************
   SUBROUTINE BeddoesModel( P,  m,  ErrStat, ErrMess, &
                      W2, J, IBlade, ALPHA, CLA, CDA ,CMA )
 !  uses the Beddoes dynamic stall model
 !   the routine is entered with an angle of attack
 !   and returns CL, CD, and CM.
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

   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

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
   REAL(ReKi)                 :: P1
   REAL(ReKi)                 :: SA
   REAL(ReKi)                 :: VREL

   INTEGER                    :: I
   INTEGER                    :: N
   INTEGER                    :: NP1

   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.
   CHARACTER(ErrMsgLen)                      :: ErrMessLcl          ! Error message returned by called routines.
   
   ErrStat = ErrID_None
   ErrMess = ""
   
 ! Check to see if element has multiple airfoil tables, then interpolate values
 !  of constants based on the current location.

I = P%AirFoil%NFOIL(J)

IF (P%AirFoil%NTables(I) > 1)THEN

   m%AirFoil%MulTabLoc = MIN( MAX( m%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1) ), P%AirFoil%MulTabMet(I,P%AirFoil%NTables(I)) )
   CALL LocateBin( m%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1:P%AirFoil%NTables(I)), N, P%AirFoil%NTables(I) )

   IF ( N == 0 ) THEN
      CNA1  = m%Beddoes%CNA( I,1)
      AOL1  = m%Beddoes%AOL( I,1)
      CNS1  = m%Beddoes%CNS( I,1)
      CNSL1 = m%Beddoes%CNSL(I,1)
      AOD1  = m%Beddoes%AOD( I,1)
      CDO1  = m%Beddoes%CDO( I,1)
   ELSE IF ( N == P%AirFoil%NTables(I) ) THEN
      CNA1  = m%Beddoes%CNA( I,N)
      AOL1  = m%Beddoes%AOL( I,N)
      CNS1  = m%Beddoes%CNS( I,N)
      CNSL1 = m%Beddoes%CNSL(I,N)
      AOD1  = m%Beddoes%AOD( I,N)
      CDO1  = m%Beddoes%CDO( I,N)
   ELSE
      NP1   = N+1
   !bjj: check for division by zero?
      P1     = (m%AirFoil%MulTabLoc-P%AirFoil%MulTabMet(I,N))/(P%AirFoil%MulTabMet(I,NP1)-P%AirFoil%MulTabMet(I,N))

      CNA1  = m%Beddoes%CNA(I,N) + P1 * ( m%Beddoes%CNA(I,NP1) - m%Beddoes%CNA(I,N) )
      AOL1  = m%Beddoes%AOL(I,N) + P1 * ( m%Beddoes%AOL(I,NP1) - m%Beddoes%AOL(I,N) )
      CNS1  = m%Beddoes%CNS(I,N) + P1 * ( m%Beddoes%CNS(I,NP1) - m%Beddoes%CNS(I,N) )
      CNSL1 = m%Beddoes%CNSL(I,N)+ P1 * ( m%Beddoes%CNSL(I,NP1)- m%Beddoes%CNSL(I,N))
      AOD1  = m%Beddoes%AOD(I,N) + P1 * ( m%Beddoes%AOD(I,NP1) - m%Beddoes%AOD(I,N) )
      CDO1  = m%Beddoes%CDO(I,N) + P1 * ( m%Beddoes%CDO(I,NP1) - m%Beddoes%CDO(I,N) )
   END IF

ELSE
   CNA1  = m%Beddoes%CNA(I,1)
   AOL1  = m%Beddoes%AOL(I,1)
   CNS1  = m%Beddoes%CNS(I,1)
   CNSL1 = m%Beddoes%CNSL(I,1)
   AOD1  = m%Beddoes%AOD(I,1)
   CDO1  = m%Beddoes%CDO(I,1)
ENDIF

 ! Jump back if lift-curve slope is zero

IF ( EqualRealNos(CNA1, 0.0_ReKi) ) THEN
   CLA = 0.0
   CDA = CDO1
   CMA = 0.0
   RETURN
ENDIF

m%Beddoes%AN   = ALPHA
VREL = SQRT( W2 )

CALL ATTACH( P, m, ErrStatLcl, ErrMessLcl, VREL, J, IBlade, CNA1, AOL1, AE )
   CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'BeddoesModel' )
   IF (ErrStat >= AbortErrLev) RETURN


CALL SEPAR(  P, m, ErrStatLcl, ErrMessLcl, P%AirFoil%NLIFT(I), J, IBlade, I, CNA1, AOL1, CNS1, CNSL1 )
   CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'BeddoesModel' )
   IF (ErrStat >= AbortErrLev) RETURN

CALL VORTEX( P, m, ErrStatLcl, ErrMessLcl, J, IBlade, AE )
   CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'BeddoesModel' )
   IF (ErrStat >= AbortErrLev) RETURN

CA  = COS( m%Beddoes%AN )
SA  = SIN( m%Beddoes%AN )
CLA = m%Beddoes%CN * CA + m%Beddoes%CC * SA
CDA = m%Beddoes%CN * SA - m%Beddoes%CC * CA + CDO1
CMA = m%AirFoil%PMC

RETURN
END SUBROUTINE BeddoesModel


 ! ******************************************************
   SUBROUTINE ATTACH( P, m, ErrStat, ErrMess, &
                      VREL, J, IBlade, CNA1, AOL1, AE )
 !  PART OF THE Beddoes DYNAMIC STALL MODEL.
 ! ******************************************************

   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),   INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),     INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER(IntKi),           INTENT(  OUT)   :: ErrStat
   CHARACTER(*),             INTENT(  OUT)   :: ErrMess

   REAL(ReKi),               INTENT(  OUT)   :: AE
   REAL(ReKi),               INTENT(IN)      :: AOL1
   REAL(ReKi),               INTENT(IN)      :: CNA1
   REAL(ReKi),               INTENT(IN)      :: VREL
                             
   INTEGER,                  INTENT(IN)      :: J
   INTEGER,                  INTENT(IN)      :: IBlade

   ! Local Variables:

   REAL(ReKi)                 :: B2
   REAL(ReKi)                 :: BS
   REAL(ReKi)                 :: CNI
   REAL(ReKi)                 :: CNQ
   REAL(ReKi)                 :: CO
   REAL(ReKi)                 :: DA
   REAL(ReKi)                 :: PRP
   REAL(ReKi)                 :: X0
   REAL(ReKi)                 :: XKA
   REAL(ReKi)                 :: XM


   ErrStat = ErrID_None
   ErrMess = ""

IF ( ABS( m%Beddoes%AN ) <= PIBY2 ) THEN
   m%Beddoes%ANE(J,IBLADE) = m%Beddoes%AN
ELSEIF ( m%Beddoes%AN > PiBy2 ) THEN
   m%Beddoes%ANE(J,IBLADE) = PI - m%Beddoes%AN
ELSE
   m%Beddoes%ANE(J,IBLADE) = - PI - m%Beddoes%AN
ENDIF


XM  = VREL / p%Beddoes%AS

 ! Check to see that the element is not supersonic
IF ( .NOT. m%SuperSonic .AND. XM >= 1.0 ) THEN
   XM = 0.7
   m%SuperSonic = .TRUE.
   ErrMess = 'ATTACH: Blade #'//TRIM(Int2LStr(IBLADE))//' element #'//TRIM(Int2LStr(J))//' is supersonic! '//&
             ' Other elements are likely supersonic as well. Supersonic mach nos. will be set to '//&
             TRIM(Num2LStr(XM))//' to attempt continuation.' 
   ErrStat = ErrID_Warn
ELSEIF (m%SuperSonic .AND. XM < 1.0) THEN
   m%SuperSonic = .FALSE.
   ErrMess = 'ATTACH: Supersonic condition has subsided with Blade #'// TRIM(Int2LStr(IBLADE))// &
                  ' element #'//TRIM(Int2LStr(J))//'.'
   ErrStat = ErrID_Info
ENDIF
IF (ErrStat >= AbortErrLev) RETURN

B2  = 1.0 - XM * XM
m%Beddoes%DS  = 2. * m%DT * VREL/p%Blade%C(J)
BS  = B2 * m%Beddoes%DS
XKA = .75/( ( 1. - XM ) + PI * B2 * XM * XM * 0.413 )
X0  = m%DT * p%Beddoes%AS / p%Blade%C(J) / XKA
CO  = XKA * p%Blade%C(J) / p%Beddoes%AS / XM

DA  = m%Beddoes%ANE(J,IBLADE) - m%Beddoes%ANE1(J,IBLADE)
m%Beddoes%ADOT(J,IBLADE) = DA / m%DT

PRP = m%Beddoes%ADOT(J,IBLADE) * p%Blade%C(J) / VREL
PRP = SAT( PRP, 0.03_ReKi, 0.1_ReKi )
m%Beddoes%ADOT(J,IBLADE) = PRP * VREL / p%Blade%C(J)

m%Beddoes%DN(J,IBLADE) = m%Beddoes%OLDDN(J,IBLADE) * EXP(-X0) + &
               (m%Beddoes%ADOT(J,IBLADE) - m%Beddoes%ADOT1(J,IBLADE)) * EXP(-.5*X0)
CNI = 4._ReKi * CO * ( m%Beddoes%ADOT(J,IBLADE) - m%Beddoes%DN(J,IBLADE) )
m%Beddoes%CMI = -.25_ReKi * CNI

m%Beddoes%QX(J,IBLADE) = (m%Beddoes%ADOT(J,IBLADE) - m%Beddoes%ADOT1(J,IBLADE)) * P%Blade%C(J)/VREL/m%DT
m%Beddoes%DQ(J,IBLADE) = m%Beddoes%OLDDQ(J,IBLADE)*EXP(-X0) + &
               ( m%Beddoes%QX(J,IBLADE) - m%Beddoes%QX1(J,IBLADE) ) * EXP(-.5*X0)
CNQ = -CO * (m%Beddoes%QX(J,IBLADE) - m%Beddoes%DQ(J,IBLADE))
m%Beddoes%DQP(J,IBLADE) = m%Beddoes%DQP1(J,IBLADE) * EXP(-X0/XKA) + (m%Beddoes%QX(J,IBLADE) &
                - m%Beddoes%QX1(J,IBLADE)) * EXP(-.5*X0/XKA)

m%Beddoes%CMQ = -.25 * CNQ - (XKA*CO/3.) * (m%Beddoes%QX(J,IBLADE) - m%Beddoes%DQP(J,IBLADE))

m%Beddoes%CNIQ = MIN( ABS( CNI+CNQ ), 1.0_ReKi ) * SIGN( 1.0_ReKi, CNI+CNQ )

m%Beddoes%XN(J,IBLADE) = m%Beddoes%OLDXN(J,IBLADE)*EXP(-.14*BS) + .3*DA*EXP(-.07*BS)
m%Beddoes%YN(J,IBLADE) = m%Beddoes%OLDYN(J,IBLADE)*EXP(-.53*BS) + .7*DA*EXP(-.265*BS)

AE   = m%Beddoes%ANE(J,IBLADE) - m%Beddoes%YN(J,IBLADE) - m%Beddoes%XN(J,IBLADE)
m%Beddoes%CNCP = CNA1 * ( AE - AOL1 )
m%Beddoes%CNPOT(J,IBLADE) = m%Beddoes%CNCP + m%Beddoes%CNIQ
m%Beddoes%CC   =  m%Beddoes%CNPOT(J,IBLADE) * AE


RETURN
END SUBROUTINE ATTACH


 ! ******************************************************
   SUBROUTINE SEPAR( P, m, ErrStat, ErrMess, &
                     NFT, J, IBlade, IFOIL, CNA1, AOL1, CNS1, CNSL1 )
!   SUBROUTINE SEPAR( NFT, J, IBlade, IFOIL, CNA1, AOL1, CNS1, CNSL1 )
 !  PART OF THE Beddoes DYNAMIC STALL MODEL
 ! ******************************************************

   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

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

   ErrStat = ErrID_None
   ErrMess = ""


TFE  = p%Beddoes%TF
m%Beddoes%DPP(J,IBLADE)   = m%Beddoes%OLDDPP(J,IBLADE) * EXP(-m%Beddoes%DS/P%Beddoes%TP) &
                 + ( m%Beddoes%CNPOT(J,IBLADE) - m%Beddoes%CNPOT1(J,IBLADE) ) * EXP(-.5*m%Beddoes%DS/P%Beddoes%TP)
m%Beddoes%CNP(J,IBLADE)  =  m%Beddoes%CNPOT(J,IBLADE) -   m%Beddoes%DPP(J,IBLADE)
m%Beddoes%CNPD(J,IBLADE) =  m%Beddoes%CNP(J,IBLADE) -  m%Beddoes%CNP1(J,IBLADE)

 ! USE CNPD to determine if AOA is increasing or decreasing.
 ! Vortex lift decays more rapidly for decreasing AOA.

IF ( m%Beddoes%ANE(J,IBLADE) * m%Beddoes%CNPD(J,IBLADE) < 0. ) THEN
   m%Beddoes%SHIFT = .TRUE.
ELSE
   m%Beddoes%SHIFT = .FALSE.
ENDIF

AFF = m%Beddoes%CNP(J,IBLADE)/CNA1 + AOL1

IF ( ABS( m%Beddoes%AN ) <= PIBY2 ) THEN
   m%Beddoes%AFE(J,IBLADE) =  AFF
ELSEIF ( m%Beddoes%AN > PIBY2 ) THEN
   m%Beddoes%AFE(J,IBLADE) = PI -  AFF
ELSE
   m%Beddoes%AFE(J,IBLADE) = -PI -  AFF
ENDIF

CALL MPI2PI ( m%Beddoes%AFE(J,IBLADE) )

IF ( ( m%Beddoes%AFE(J,IBLADE) < m%AirFoil%AL(IFOIL,1) ) .OR.  &
     ( m%Beddoes%AFE(J,IBLADE) > m%AirFoil%AL(IFOIL,p%AirFoil%NLIFT(IFOIL) ) ) ) THEN
   ErrMess = 'SEPAR: Angle of attack = '//TRIM(Num2LStr( m%Beddoes%AFE(J,IBLADE)*R2D ))//' is outside table.'
   ErrStat = ErrID_Fatal
   RETURN
ENDIF

m%Beddoes%AFE(J,IBLADE) = MIN( MAX( m%Beddoes%AFE(J,IBLADE), m%AirFoil%AL(IFOIL,1) ), m%AirFoil%AL(IFOIL,NFT) )
CALL LocateBin( m%Beddoes%AFE(J,IBLADE), m%AirFoil%AL(IFOIL,1:NFT), I1, NFT )

IF (I1 == 0) THEN
   I1 = 1
ELSE IF ( I1 == NFT ) THEN
   I1 = I1 - 1
END IF

I1P1 = I1 + 1

   !bjj: check for division by zero?
P1 = ( m%AirFoil%AL(IFOIL,I1) - m%Beddoes%AFE(J,IBLADE) ) / ( m%AirFoil%AL(IFOIL,I1) - m%AirFoil%AL(IFOIL,I1P1) )

IF ( p%AirFoil%NTables(IFOIL) > 1 ) THEN

 ! Locate the multiple airfoil position in the table

   m%AirFoil%MulTabLoc = MIN( MAX( m%AirFoil%MulTabLoc, p%AirFoil%MulTabMet(IFOIL,1) ), p%AirFoil%MulTabMet(IFOIL,p%AirFoil%NTables(IFOIL)) )
   CALL LocateBin( m%AirFoil%MulTabLoc, p%AirFoil%MulTabMet(IFOIL,1:p%AirFoil%NTables(IFOIL)),I2,p%AirFoil%NTables(IFOIL) )

   IF ( I2 == 0 ) THEN
      I2   = 1
      I2P1 = 2
      P2   = 0.0
   ELSE IF ( I2 == p%AirFoil%NTables(IFOIL) ) THEN
      I2P1 = I2
      I2   = I2 - 1

      P2   = 1.0
   ELSE
      I2P1 = I2 + 1

      P2=(m%AirFoil%MulTabLoc-p%AirFoil%MulTabMet(IFOIL,I2))/(p%AirFoil%MulTabMet(IFOIL,I2P1)-p%AirFoil%MulTabMet(IFOIL,I2))
   END IF

 ! Interpolate the F-table values


   FSPB  = m%Beddoes%FTB( IFOIL,I1,I2P1) - (m%Beddoes%FTB( IFOIL,I1,I2P1) - m%Beddoes%FTB( IFOIL,I1P1,I2P1) ) * P1
   FSPCB = m%Beddoes%FTBC(IFOIL,I1,I2P1) - (m%Beddoes%FTBC(IFOIL,I1,I2P1) - m%Beddoes%FTBC(IFOIL,I1P1,I2P1) ) * P1
   FSPA  = m%Beddoes%FTB( IFOIL,I1,I2  ) - (m%Beddoes%FTB( IFOIL,I1,I2  ) - m%Beddoes%FTB( IFOIL,I1P1,I2  ) ) * P1
   FSPCA = m%Beddoes%FTBC(IFOIL,I1,I2  ) - (m%Beddoes%FTBC(IFOIL,I1,I2  ) - m%Beddoes%FTBC(IFOIL,I1P1,I2  ) ) * P1

   m%Beddoes%FSP(J,IBLADE) = FSPA  + P2 * (FSPB-FSPA)

   m%Beddoes%FSPC(J,IBLADE)= FSPCA + P2 * (FSPCB-FSPCA)

ELSE

   m%Beddoes%FSP(J,IBLADE) = m%Beddoes%FTB(IFOIL,I1,1) - (m%Beddoes%FTB(IFOIL,I1,1) - m%Beddoes%FTB(IFOIL,I1P1,1) ) * P1

   m%Beddoes%FSPC(J,IBLADE)= m%Beddoes%FTBC(IFOIL,I1,1) - (m%Beddoes%FTBC(IFOIL,I1,1) - m%Beddoes%FTBC(IFOIL,I1P1,1) ) * P1

ENDIF

IF ( ABS( m%Beddoes%AFE(J,IBLADE) - AOL1 ) < 1.E-10 ) THEN
   m%Beddoes%FSP(J,IBLADE)  = 1.0
   m%Beddoes%FSPC(J,IBLADE) = 1.0
ELSE
   TEMP = 2.*SQRT(ABS(m%Beddoes%FSP(J,IBLADE)/(m%Beddoes%AFE(J,IBLADE)-AOL1)))-1.
   m%Beddoes%FSP(J,IBLADE) = TEMP * TEMP * SIGN ( 1.0_ReKi, TEMP )
   IF ( m%Beddoes%FSP(J,IBLADE) >  1.0 ) m%Beddoes%FSP(J,IBLADE) =  1.0
   IF ( m%Beddoes%FSP(J,IBLADE) < -1.0 ) m%Beddoes%FSP(J,IBLADE) = -1.0

   IF ( ABS( m%Beddoes%AFE(J,IBLADE) ) < 1.E-10 ) THEN
      m%Beddoes%FSPC(J,IBLADE) = 1.0
   ELSE
      TEMP = m%Beddoes%FSPC(J,IBLADE)/((m%Beddoes%AFE(J,IBLADE)-AOL1)*m%Beddoes%AFE(J,IBLADE))
      m%Beddoes%FSPC(J,IBLADE) = TEMP * TEMP * SIGN ( 1.0_ReKi, TEMP )
      IF ( m%Beddoes%FSPC(J,IBLADE) >  1.0 ) m%Beddoes%FSPC(J,IBLADE) =  1.0
      IF ( m%Beddoes%FSPC(J,IBLADE) < -1.0 ) m%Beddoes%FSPC(J,IBLADE) = -1.0
   ENDIF

ENDIF

IF ( m%Beddoes%CNP(J,IBLADE) > CNS1  ) m%Beddoes%BEDSEP(J,IBLADE) = .TRUE.
IF ( m%Beddoes%CNP(J,IBLADE) < CNSL1 ) m%Beddoes%BEDSEP(J,IBLADE) = .TRUE.

IF ( m%Beddoes%BEDSEP(J,IBLADE) ) m%Beddoes%TAU(J,IBLADE) = m%Beddoes%OLDTAU(J,IBLADE) + m%Beddoes%DS/P%Beddoes%TVL

IF (m%Beddoes%SHIFT) TFE = 1.5*TFE

m%Beddoes%DF(J,IBLADE) = m%Beddoes%OLDDF(J,IBLADE) * EXP(-m%Beddoes%DS/TFE)  &
              + (m%Beddoes%FSP(J,IBLADE) - m%Beddoes%FSP1(J,IBLADE)) * EXP(-.5*m%Beddoes%DS/TFE)
m%Beddoes%DFC(J,IBLADE)= m%Beddoes%OLDDFC(J,IBLADE) * EXP(-m%Beddoes%DS/TFE) &
              + (m%Beddoes%FSPC(J,IBLADE) - m%Beddoes%FSPC1(J,IBLADE)) * EXP(-.5*m%Beddoes%DS/TFE)

m%Beddoes%FP   = m%Beddoes%FSP(J,IBLADE) - m%Beddoes%DF(J,IBLADE)
m%Beddoes%FPC  = m%Beddoes%FSPC(J,IBLADE) - m%Beddoes%DFC(J,IBLADE)
SRFP = SQRT( ABS(m%Beddoes%FP) )  * SIGN( 1.0_ReKi, m%Beddoes%FP ) + 1.
SRFPC= SQRT( ABS(m%Beddoes%FPC) ) * SIGN( 1.0_ReKi,m%Beddoes%FPC )

m%Beddoes%FK   = 0.25 * SRFP * SRFP
m%Beddoes%CN   = m%Beddoes%CNCP * m%Beddoes%FK + m%Beddoes%CNIQ

m%Beddoes%CC   =  m%Beddoes%CC * SRFPC

m%Beddoes%DFAFE(J,IBLADE) = m%Beddoes%DFAFE1(J,IBLADE) * EXP(-m%Beddoes%DS/(.1*TFE)) &
                 + (m%Beddoes%AFE(J,IBLADE) - m%Beddoes%AFE1(J,IBLADE)) * EXP(-.5*m%Beddoes%DS/(.1*TFE))

AFEP=m%Beddoes%AFE(J,IBLADE) - m%Beddoes%DFAFE(J,IBLADE)


AFEP = MIN( MAX( AFEP, m%AirFoil%AL(IFOIL,1) ), m%AirFoil%AL(IFOIL,NFT) )
CALL LocateBin( AFEP, m%AirFoil%AL(IFOIL,1:NFT), I1, NFT )

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
   P1 = (m%AirFoil%AL(IFOIL,I1) - AFEP)/(m%AirFoil%AL(IFOIL,I1) - m%AirFoil%AL(IFOIL,I1P1))
END IF


IF (p%AirFoil%NTables(IFOIL) > 1) THEN

 ! Interpolate the F-table values

   CMPB = m%AirFoil%CM(IFOIL,I1,I2P1) - (m%AirFoil%CM(IFOIL,I1,I2P1) - m%AirFoil%CM(IFOIL,I1P1,I2P1) ) * P1

   CMPA = m%AirFoil%CM(IFOIL,I1,I2) - (m%AirFoil%CM(IFOIL,I1,I2) - m%AirFoil%CM(IFOIL,I1P1,I2) ) * P1

   m%AirFoil%PMC = CMPA  + P2*(CMPB-CMPA)

ELSE

   m%AirFoil%PMC = m%AirFoil%CM(IFOIL,I1,1) - ( (m%AirFoil%CM(IFOIL,I1,1) - m%AirFoil%CM(IFOIL,I1P1,1) ) * P1 )

ENDIF



RETURN
END SUBROUTINE SEPAR


 ! ******************************************************
   SUBROUTINE VORTEX( P, m, ErrStat, ErrMess, &
                      J, IBlade, AE )
 !  PART OF THE Beddoes DYNAMIC STALL MODEL
 ! ******************************************************

   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   ! Passed Variables:
   REAL(ReKi),INTENT(IN)      :: AE

   INTEGER,   INTENT(IN)      :: J
   INTEGER,   INTENT(IN)      :: IBlade

   ! Local Variables:
   REAL(ReKi)                 :: CMV
   REAL(ReKi)                 :: TSH
   REAL(ReKi)                 :: TVE

   ErrStat = ErrID_None
   ErrMess = ""

TVE = P%Beddoes%TV
m%Beddoes%CVN(J,IBLADE) = m%Beddoes%CNCP * ( 1. - m%Beddoes%FK )

IF ( m%Beddoes%TAU(J,IBLADE) < 1. ) THEN
   m%Beddoes%VOR = .TRUE.
   IF (m%Beddoes%SHIFT) m%Beddoes%VOR = .FALSE.
ELSE
   m%Beddoes%VOR = .FALSE.
ENDIF

IF (m%Beddoes%VOR) THEN
   m%Beddoes%CNV(J,IBLADE) = m%Beddoes%OLDCNV(J,IBLADE) * EXP(-m%Beddoes%DS/TVE)  &
                 + (m%Beddoes%CVN(J,IBLADE) - m%Beddoes%CVN1(J,IBLADE)) * EXP(-.5*m%Beddoes%DS/TVE)
ELSE
   m%Beddoes%CNV(J,IBLADE) = m%Beddoes%OLDCNV(J,IBLADE) * EXP(-m%Beddoes%DS/(TVE*.5))
ENDIF

m%Beddoes%CN  = m%Beddoes%CN + m%Beddoes%CNV(J,IBLADE)
m%Beddoes%CC  = m%Beddoes%CC + m%Beddoes%CNV(J,IBLADE) * AE * (1.- m%Beddoes%TAU(J,IBLADE))
CMV = -0.2 * (1. - COS(PI*m%Beddoes%TAU(J,IBLADE)) ) * m%Beddoes%CNV(J,IBLADE)
m%AirFoil%PMC = m%AirFoil%PMC + CMV + m%Beddoes%CMI + m%Beddoes%CMQ

TSH = 2.*(1.- m%Beddoes%FP)/.19

IF ( m%Beddoes%TAU(J,IBLADE) .GT. 1. + TSH/p%Beddoes%TVL .AND. .NOT. m%Beddoes%SHIFT) THEN
   m%Beddoes%TAU(J,IBLADE) = 0.
   m%Beddoes%BEDSEP(J,IBLADE) = .FALSE.
ENDIF

IF ( m%Beddoes%TAU(J,IBLADE) .GT. 1. ) THEN
   IF ( m%Beddoes%ANE(J,IBLADE) .LT. 0. ) THEN

      IF (m%Beddoes%CNPD(J,IBLADE) .LE. 0. .AND. m%Beddoes%CNPD1(J,IBLADE) .GE. 0.) THEN
         m%Beddoes%BEDSEP(J,IBLADE) = .FALSE.
         m%Beddoes%TAU(J,IBLADE) = 0.
      ENDIF

      IF (m%Beddoes%ANE1(J,IBLADE) .GT. 0.) THEN
         m%Beddoes%BEDSEP(J,IBLADE) = .FALSE.
         m%Beddoes%TAU(J,IBLADE) = 0.
      ENDIF

   ELSE

      IF (m%Beddoes%CNPD(J,IBLADE) .GE. 0. .AND. m%Beddoes%CNPD1(J,IBLADE) .LE. 0.) THEN
         m%Beddoes%BEDSEP(J,IBLADE) = .FALSE.
         m%Beddoes%TAU(J,IBLADE) = 0.
      ENDIF

      IF (m%Beddoes%ANE1(J,IBLADE) .LT. 0.) THEN
         m%Beddoes%BEDSEP(J,IBLADE) = .FALSE.
         m%Beddoes%TAU(J,IBLADE) = 0.
      ENDIF

   ENDIF
ENDIF



RETURN
END SUBROUTINE VORTEX


 ! ******************************************************
   SUBROUTINE CLCD( P,  m,  ErrStat, ErrMess, &
                    ALPHA, CLA, CDA, CMA, I )
!   SUBROUTINE CLCD( ALPHA, CLA, CDA, CMA, I, ErrStat )
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
!USE                           Airfoil
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   ! Passed Variables:
   REAL(ReKi),INTENT(INOUT)   :: ALPHA
   REAL(ReKi),INTENT(OUT)     :: CDA
   REAL(ReKi),INTENT(OUT)     :: CLA
   REAL(ReKi),INTENT(OUT)     :: CMA

   INTEGER   ,INTENT(IN)      :: I      ! NFOIL(J)

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

ErrStat = ErrID_None
ErrMess = ""

IF (.NOT. ALLOCATED(P%AirFoil%NFoil) ) THEN
   CDA = 0.
   CLA = 0.
   CMA = 0.
   ErrStat = ErrID_Fatal
   RETURN
ELSE
   ErrStat = ErrID_None
END IF

NTAB = P%AirFoil%NLIFT(I)

IF ( ( ALPHA < m%AirFoil%AL(I,1) ) .OR. ( ALPHA > m%AirFoil%AL(I,NTAB) ) )   THEN
!bjj: This error message isn't necessarially accurate:
   CDA = 0.
   CLA = 0.
   CMA = 0.
  ErrMess = ' Angle of attack = '//TRIM(Num2LStr(ALPHA*R2D))// &
            ' deg is outside data table range. '// & !Blade #'//TRIM(Int2LStr(IBLADE))//&
            ' Airfoil '//TRIM(Int2LStr(I))//'.' 
!                   ' element '//TRIM(Int2LStr(J))//'.' )

   ErrStat = ErrID_Fatal
   RETURN
ENDIF

ALPHA = MIN( MAX( ALPHA, m%AirFoil%AL(I,1) ), m%AirFoil%AL(I,NTAB) )
CALL LocateBin (ALPHA, m%AirFoil%AL(I,1:NTAB), N1, NTAB )

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
   P1   = ( ALPHA - m%AirFoil%AL(I, N1) )/( m%AirFoil%AL(I, N1P1) - m%AirFoil%AL(I, N1) )
END IF

 ! If the element has multiple airfoil tables, do a 2-D linear interpolation
 !  for Cl and CD

IF (P%AirFoil%NTables(I) > 1) THEN

   m%AirFoil%MulTabLoc = MIN( MAX( m%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1) ), P%AirFoil%MulTabMet(I,P%AirFoil%NTables(I)) )
   CALL LocateBin (m%AirFoil%MulTabLoc, P%AirFoil%MulTabMet(I,1:P%AirFoil%NTables(I)),N2,P%AirFoil%NTables(I))

   IF (N2 == 0) THEN
      N2   = 1
      N2P1 = 2
      P2   = 0.0
   ELSE IF ( N2 == P%AirFoil%NTables(I) ) THEN
      N2P1 = N2
      N2   = N2 - 1
      P2   = 1.0
   ELSE
      N2P1 = N2 + 1
      P2   = (m%AirFoil%MulTabLoc - P%AirFoil%MulTabMet(I,N2))/(P%AirFoil%MulTabMet(I,N2P1)-P%AirFoil%MulTabMet(I,N2))
   END IF

   CLA1 = m%AirFoil%CL(I,N1,N2) + P1 * ( m%AirFoil%CL(I,N1P1,N2) - m%AirFoil%CL(I,N1,N2) )
   CDA1 = m%AirFoil%CD(I,N1,N2) + P1 * ( m%AirFoil%CD(I,N1P1,N2) - m%AirFoil%CD(I,N1,N2) )
   CMA1 = m%AirFoil%CM(I,N1,N2) + P1 * ( m%AirFoil%CM(I,N1P1,N2) - m%AirFoil%CM(I,N1,N2) )

   CLA2 = m%AirFoil%CL(I,N1,N2P1) + P1 * ( m%AirFoil%CL(I,N1P1,N2P1) - m%AirFoil%CL(I,N1,N2P1) )
   CDA2 = m%AirFoil%CD(I,N1,N2P1) + P1 * ( m%AirFoil%CD(I,N1P1,N2P1) - m%AirFoil%CD(I,N1,N2P1) )
   CMA2 = m%AirFoil%CM(I,N1,N2P1) + P1 * ( m%AirFoil%CM(I,N1P1,N2P1) - m%AirFoil%CM(I,N1,N2P1) )

   CLA = CLA1 + P2 * ( CLA2 - CLA1 )
   CDA = CDA1 + P2 * ( CDA2 - CDA1 )
   CMA = CMA1 + P2 * ( CMA2 - CMA1 )

ELSE

   CLA  = m%AirFoil%CL(I,N1,1) + P1 * ( m%AirFoil%CL(I,N1P1,1) - m%AirFoil%CL(I,N1,1) )
   CDA  = m%AirFoil%CD(I,N1,1) + P1 * ( m%AirFoil%CD(I,N1P1,1) - m%AirFoil%CD(I,N1,1) )
   CMA  = m%AirFoil%CM(I,N1,1) + P1 * ( m%AirFoil%CM(I,N1P1,1) - m%AirFoil%CM(I,N1,1) )

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

   SUBROUTINE Inflow( Time, P, m, ErrStat, ErrMess )

 ! Gateway to the dynamic inflow routines.
 ! Called by NEWTIME after a time step.
 ! **************************************

   IMPLICIT                      NONE
      ! Passed Variables:
   REAL(DbKi),                   INTENT(IN)  :: Time
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.
   CHARACTER(ErrMsgLen)                      :: ErrMessLcl          ! Error message returned by called routines.
   
   ErrStat = ErrID_None
   ErrMess = ""
   
      IF ( P%DYNINFL ) THEN

! INITIALIZE DYNAMIC INFLOW PARAMETERS
        IF ( m%DYNINIT .AND. ( TIME > 0.0D0 ) ) THEN
           CALL INFINIT( P, m, ErrStatLcl, ErrMessLcl )
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'Inflow' )
            IF (ErrStat >= AbortErrLev) RETURN
           !WRITE(*,*) 'Activating dynamic inflow calculation'
           m%DynInflow%old_Alph = 0.0
           m%DynInflow%old_Beta = 0.0
           m%DYNINIT = .FALSE.
           m%SKEW = .FALSE.
        ENDIF

 ! Update the dynamic inflow parameters
         CALL INFUPDT(P, m, ErrStatLcl, ErrMessLcl)
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'Inflow' )
            IF (ErrStat >= AbortErrLev) RETURN
            
         CALL GetPhiLq(P, m, ErrStatLcl, ErrMessLcl)
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'Inflow' )
            IF (ErrStat >= AbortErrLev) RETURN
            
 ! Calculate the dynamic inflow paremeters for the new time step
         IF( TIME > 1.0D0 ) THEN
            CALL INFDIST(P, m, ErrStatLcl, ErrMessLcl)
               CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'Inflow' )
               IF (ErrStat >= AbortErrLev) RETURN
         END IF
         
      ENDIF   ! DYNINFL

RETURN
END SUBROUTINE Inflow


 ! **************************************

   SUBROUTINE GetRM ( P,  m,  ErrStat, ErrMess, &
                      rLocal, DFN, DFT, psi, J, IBlade)

 ! Returns RM(MODE), the [mode]-th moment of the blade normal force.
 !  Here, the force is in [N/m], while the moment arm is
 !  non-dimensional (RLOCAL/R).  Also see FUNCTION XPHI.
 ! Called as each element is processed by AeroDyn.
 ! **************************************
   IMPLICIT                      NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables

   REAL(ReKi),INTENT(IN)      :: DFN
   REAL(ReKi),INTENT(IN)      :: DFT
   REAL(ReKi),INTENT(IN)      :: psi
   REAL(ReKi),INTENT(IN)      :: rLocal

   INTEGER,   INTENT(IN)      :: J
   INTEGER,   INTENT(IN)      :: IBlade
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   ! Local Variables:
   REAL(ReKi)                 :: fElem
   ! psiBar is Suzuki's, WindPsi is Shawler's
   !REAL(ReKi)                :: psiBar
   REAL(ReKi)                 :: Rzero
   REAL(ReKi)                 :: WindPsi

   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.
   CHARACTER(ErrMsgLen)                      :: ErrMessLcl          ! Error message returned by called routines.
   
   
   INTEGER                    :: mode

   ErrStat = ErrID_None
   ErrMess = ""


IF ( P%SWIRL ) THEN
   fElem = SQRT( DFN * DFN + DFT * DFT )
ELSE
   fElem = DFN
ENDIF
fElem = fElem / P%Blade%R

Rzero = rLocal / P%Blade%R
! Suzukis inflow azimuth measure
!psiBar = - psi - piBy2
! Shawler: wind based inflow azimuth measure
CALL WindAzimuthZero (psi,m%Wind%VrotorY,m%Wind%VrotorZ,WindPsi)

! Save values rotor loads for USE in Newtime (to accumulate rotor loads)
DO mode = 1, P%DynInflow%MAXINFLO
   m%DynInflow%RMC_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode, ErrStatLcl, ErrMessLcl )
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'GetRM' )
      IF (ErrStat >= AbortErrLev) RETURN
END DO ! mode

!+++++++++++++++++++++++++++++++++++++++++++++++++++++
!Suzuki's method
!DO mode = MaxInflo+1, maxInfl
!   m%DynInflow%RMC_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode )  * COS( REAL(MRvector(mode), ReKi) * psiBar )
!   m%DynInflow%RMS_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode )  * SIN( REAL(MRvector(mode), ReKi) * psiBar )
!END DO ! mode
! Shawler's method
DO mode = p%DynInflow%MaxInflo+1, maxInfl
   m%DynInflow%RMC_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode, ErrStatLcl, ErrMessLcl )  * COS( REAL(MRvector(mode), ReKi) * WindPsi )
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'GetRM' )
      IF (ErrStat >= AbortErrLev) RETURN
   
   m%DynInflow%RMS_SAVE(IBLADE, J, mode) = fElem * XPHI( Rzero, mode, ErrStatLcl, ErrMessLcl )  * SIN( REAL(MRvector(mode), ReKi) * WindPsi )
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'GetRM' )
      IF (ErrStat >= AbortErrLev) RETURN

END DO ! mode
!++++++++++++++++++++++++++++++++++++++++++++++++++++++

RETURN
END SUBROUTINE GetRM


 ! **************************************

   SUBROUTINE GetPhiLq( P, m, ErrStat, ErrMess )

 ! Accumulate the rotor forces for dynamic inflow calculations
 !  PhiLqC is Lq times PHI (shape function) in COS equation.
 !  PhiLqS is Lq times PHI (shape function) in SIN equation.
 !  RM?_SAVE, which were calculated for each element,
 !   are summed here after a time step.
 ! **************************************

!USE            Blade
!USE            DynInflow
!USE            Element

   IMPLICIT       NONE
      ! Passed Variables:
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess


   INTEGER     :: iblad
   INTEGER     :: ielem
   INTEGER     :: mode

   ErrStat = ErrID_None
   ErrMess = ""

m%DynInflow%PhiLqC = 0.
m%DynInflow%PhiLqS = 0.

DO mode = 1, maxInfl
   DO iblad = 1, P%NumBl
      DO ielem = 1, P%Element%NELM
         m%DynInflow%PhiLqC(mode) = m%DynInflow%PhiLqC(mode) + m%DynInflow%RMC_SAVE(iblad, ielem, mode)
         IF (mode >= p%DynInflow%MaxInflo+1) &
         m%DynInflow%PhiLqS(mode) = m%DynInflow%PhiLqS(mode) + m%DynInflow%RMS_SAVE(iblad, ielem, mode)
      END DO !ielem
   END DO !iblad
END DO !mode

RETURN
END SUBROUTINE GetPhiLq


 ! *************************************************************
   SUBROUTINE infinit( P, m, ErrStat, ErrMess )
 !  Initializes the variables in the dynamic inflow equation.
 !  Called only once to initialize the GDW parameters.
 ! *************************************************************
      IMPLICIT                      NONE

      TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
      INTEGER, INTENT(OUT)                   :: ErrStat
      CHARACTER(*), INTENT(OUT)              :: ErrMess

   ! Local Variables:

      REAL(ReKi)                 :: tauCos ( maxInfl )
      REAL(ReKi)                 :: tauSin ( P%DynInflow%MaxInflo+1 : maxInfl )
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

   ErrStat = ErrID_None
   ErrMess = ""

 ! Initialize the MRvector & NJVector.
 !  MRvector is the azimuthal mode of the inflow distribution.
 !  NJvector is the radial mode of the inflow distribution.

!P%DynInflow%MRVector(:) = (/ 0, 0, 1, 1, 2, 3 /)  !bjj why aren't these parameters? Now they are.
!P%DynInflow%NJVector(:) = (/ 1, 3, 2, 4, 3, 4 /)

 ! For your reference,
 !  Marray(irow,jcol) = MRvector(jcol)
 !  Rarray(irow,jcol) = MRvector(irow)
 !  Narray(irow,jcol) = NJvector(jcol)
 !  Jarray(irow,jcol) = NJvector(irow)

 ! Initialize the time derivatives.
m%DynInflow%dAlph_dt(:,:) = 0.
m%DynInflow%dBeta_dt(:,:) = 0.

! ! Set up the constants.
! !  xMinv is [M]^-1.  Because [M]^-1 is just a diagonal matrix,
! !  it is stored as a column vector.
!bjj: this is a parameter that we'll set at initialization
!DO irow = 1, maxInfl
!   p%DynInflow%xMinv(irow) = PIBY2 / hfunc(MRvector(irow), NJvector(irow))   !bjj: this is really just a parameter, too.
!END DO !irow

 ! Set up the GAMMA matrix which is used to calculate [L] matrix.
 !  FUNCTION FGAMMA is called.
DO irow = 1, maxInfl
   DO jcol = 1, maxInfl
      m%DynInflow%gamma( irow, jcol )  &
           = fgamma( MRvector( irow ), NJvector( irow ),  &
                     MRvector( jcol ), NJvector( jcol )  )  !bjj: this is really just a parameter, too.
   END DO !jcol
END DO !irow

 ! calculate and store the M-R matrices, which are used in Subroutine LMATRIX.
DO irow = 1, maxInfl
   DO jcol = 1, maxInfl
      m%DynInflow%MminR  (irow,jcol) = MIN( MRvector(jcol) , MRvector(irow) )
      m%DynInflow%MplusR (irow,jcol) =      MRvector(jcol) + MRvector(irow)
      m%DynInflow%MminusR(irow,jcol) = ABS( MRvector(jcol) - MRvector(irow) )
   END DO !jcol
END DO !irow

 ! Calculate the tip speed of the rotor. This isn't constant in ADAMS.
 !  Thus, it will be updated at every time step in ADAMS.
m%DynInflow%TipSpeed = MAX(P%Blade%r * m%Rotor%revs, 1.0e-6_ReKi)

 ! Calculate the disk loading normalization factor.
 !  This is not exactly pressure but let's call it P0.
 !  The actual unit is [N/m] or [Pa*m].
 !    Pzero = PI * AirDensity * (Rotational Speed)^2 * (Radius)^3
 !    Pzero = pi * rho * revs**2 * P%Blade%r**3
 !    Pzero = pi * rho * revs * revs * P%Blade%r * P%Blade%r * P%Blade%r
m%DynInflow%Pzero = pi * p%Wind%rho * m%DynInflow%TipSpeed * m%DynInflow%TipSpeed * P%Blade%r
 ! Non-dimensional time
m%DynInflow%DTO   = m%DT * m%Rotor%REVS !bjj: this isn't used in this subroutine?

 ! Calculate the initial values of inflow distribution parameters

 ! Calculate the non-dimensional wind velocity components.

v1 =   m%Wind%VrotorZ / m%DynInflow%TipSpeed     !inplane, upward
v2 =   m%Wind%VrotorY / m%DynInflow%TipSpeed     !inplane, right looking downwind
v3 = - m%Wind%VrotorX / m%DynInflow%TipSpeed     !out-of-plane, normal to the rotor

 ! Calculate the initial value of lambda_m by taking the average
 !  of the induction factors(A). The A's are calculated by
 !  momentum balance during the first rotation of the trim solution.
m%DynInflow%xLambda_M = 0.
DO jBlade=1,P%NumBl
   DO iElem =1,P%Element%nelm
      m%DynInflow%xLambda_M = m%DynInflow%xLambda_M + m%Element%a(iElem,jBlade)
   END DO !iElem
END DO !jBlade
m%DynInflow%xLambda_M = m%DynInflow%xLambda_M / ( P%NumBl * P%Element%nelm )

 ! A's are normalized by the normal wind speed, while Lambda's are
 !  mormalized by the tip speed. Make the conversion.
 !  xLambda_M = xLambda_M * (-VrotorX/TipSpeed)
m%DynInflow%xLambda_M = m%DynInflow%xLambda_M * v3

 ! totalInf is the non-dimensional total wind speed normal to the rotor.
m%DynInflow%totalInf = - v3 + m%DynInflow%xLambda_M
 ! Vplane2 is the square of in-plane component of the non-dimensional
 !  wind velocity.
Vplane2 = v1 * v1 + v2 * v2
 ! VTOTAL is the total wind speed at the rotor.
m%DynInflow%Vtotal  = SQRT( m%DynInflow%totalInf * m%DynInflow%totalInf + Vplane2 )
 ! VPARAM is the velocity parameter.
m%DynInflow%Vparam  =( Vplane2 + ( m%DynInflow%totalInf + m%DynInflow%old_LmdM ) * m%DynInflow%totalInf ) / m%DynInflow%Vtotal

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
m%DynInflow%xkai = TAN( .5 * ATAN( SQRT( Vplane2 ) / m%DynInflow%totalInf ) )
 !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ! To calculate the initial values of xAlpha & xBeta, get [L] matrices
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ! Suzuki:
!CALL LMATRIX ( m,p,xKaiC, 1 )
!CALL LMATRIX ( m,p,xKaiS, 2 )
 ! Shawler:
CALL LMATRIX ( m,p, 1, ErrStat, ErrMess )
   IF (ErrStat /= ErrID_None) THEN
      ErrMess='infinit:'//TRIM(ErrMess)
      RETURN
   END IF
CALL LMATRIX ( m,p, 2, ErrStat, ErrMess )
   IF (ErrStat /= ErrID_None) THEN
      ErrMess='infinit:'//TRIM(ErrMess)
      RETURN
   END IF
   
! CALL LMATRIX ( xkai )
 ! Here we need [L_cos] & [L_sin], not [L_***}^-1.
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ! Get TAU's (pressure coefficient). Refer to Subroutine INFDIST.
DO mode = 1, P%DynInflow%MaxInflo
   tauCos(mode) = - m%DynInflow%PhiLqC(mode) / m%DynInflow%Pzero * .5
END DO !mode

DO mode = P%DynInflow%MaxInflo+1, maxInfl
   tauCos(mode) = - m%DynInflow%PhiLqC(mode) / m%DynInflow%Pzero
   tauSin(mode) = - m%DynInflow%PhiLqS(mode) / m%DynInflow%Pzero
END DO !mode

 ! Get the steady values of alpha(1)
 !  If m=0 and n=1, USE VTOTAL.
m%DynInflow%xAlpha(1) = 0.
DO k = 1, maxInfl
   m%DynInflow%xAlpha(1) = m%DynInflow%xAlpha(1) + m%DynInflow%xLcos(1,k) * tauCos(k)
END DO !k
m%DynInflow%xAlpha(1) = .5 * m%DynInflow%xAlpha(1) / m%DynInflow%Vtotal

 ! If m=0 but NOT n=1, USE VPARAM.
DO i = 2, P%DynInflow%MaxInflo
   m%DynInflow%xAlpha(i) = 0.
   DO k = 1, maxInfl
      m%DynInflow%xAlpha(i) = m%DynInflow%xAlpha(i) + m%DynInflow%xLcos(i,k) * tauCos(k)
   END DO !k
   m%DynInflow%xAlpha(i) = .5 * m%DynInflow%xAlpha(i) / m%DynInflow%Vparam
END DO !i

 ! Get the steady values of alpha's & beta's
DO i = P%DynInflow%MaxInflo+1, maxInfl
   m%DynInflow%xAlpha(i) = 0.
 ! akihiro
!   DO k = 1, maxInfl
   DO k = P%DynInflow%MaxInflo+1, maxInfl
      m%DynInflow%xAlpha(i) = m%DynInflow%xAlpha(i) + m%DynInflow%xLcos(i,k) * tauCos(k)
   END DO !k
   m%DynInflow%xAlpha(i) = .5 * m%DynInflow%xAlpha(i) / m%DynInflow%Vparam

   m%DynInflow%xBeta (i) = 0.
   DO k = P%DynInflow%MaxInflo+1, maxInfl
      m%DynInflow%xBeta (i) = m%DynInflow%xBeta (i) + m%DynInflow%xLsin(i,k) * tauSin(k)
   END DO !k
   m%DynInflow%xBeta (i) = .5 * m%DynInflow%xBeta (i) / m%DynInflow%Vparam
END DO !i

 ! Invert [L_cos] & [L_sin] matrices in case the XKAI is constant
 !  and the same [L]'s are used later.
CALL MATINV  ( m%DynInflow%xLcos, m%DynInflow%xLsin, maxInfl, P%DynInflow%MaxInflo, 1, ErrStat, ErrMess )
   IF (ErrStat /= ErrID_None) THEN
      ErrMess='infinit:'//TRIM(ErrMess)
      RETURN
   END IF
CALL MATINV  ( m%DynInflow%xLcos, m%DynInflow%xLsin, maxInfl, P%DynInflow%MaxInflo, 2, ErrStat, ErrMess )
   IF (ErrStat /= ErrID_None) THEN
      ErrMess='infinit:'//TRIM(ErrMess)
      RETURN
   END IF

RETURN
END SUBROUTINE infinit


 ! **********************************************************************
   SUBROUTINE infupdt( P, m, ErrStat, ErrMess )
 !  INFUPDT updates the OLD variables of inflow distribution parameters.
 !   The program must call this subroutine before calling
 !   the subroutine INFDIST.
 ! **********************************************************************

!USE            DynInflow

   IMPLICIT       NONE

   TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess
   INTEGER     :: i


   ErrStat = ErrID_None
   ErrMess = ""

!+++++++++++++++++++++++++++
!Suzuki:
!oldKaiC     = xKaiC
!oldKaiS     = xKaiS
!Shawler:
m%DynInflow%oldKai = m%DynInflow%xKai
!+++++++++++++++++++++++++++

m%DynInflow%old_LmdM    = m%DynInflow%xLambda_M

DO i = 1, maxInfl
   m%DynInflow%old_Alph(i)   = m%DynInflow%xAlpha  (i)
   m%DynInflow%dAlph_dt(i,4) = m%DynInflow%dAlph_dt(i,3)
   m%DynInflow%dAlph_dt(i,3) = m%DynInflow%dAlph_dt(i,2)
   m%DynInflow%dAlph_dt(i,2) = m%DynInflow%dAlph_dt(i,1)
END DO !i

DO i = P%DynInflow%MaxInflo+1, maxInfl
   m%DynInflow%old_Beta(i)   = m%DynInflow%xBeta   (i)
   m%DynInflow%dBeta_dt(i,4) = m%DynInflow%dBeta_dt(i,3)
   m%DynInflow%dBeta_dt(i,3) = m%DynInflow%dBeta_dt(i,2)
   m%DynInflow%dBeta_dt(i,2) = m%DynInflow%dBeta_dt(i,1)
END DO !i

RETURN
END SUBROUTINE infupdt


 ! **********************************************************************
   SUBROUTINE DynDebug (Time, P, x, xd, z, m, y, ErrStat, ErrMess, RHScos, RHSsin)
 !  Write out debugging information
 ! **********************************************************************

   IMPLICIT                      NONE
   REAL(DbKi),                   INTENT(IN)  :: Time
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_ContinuousStateType), INTENT(INOUT)  :: x           ! Initial continuous states
   TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(AD14_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   ! Passed Variables:

   REAL(ReKi)                 :: RHScos ( maxInfl )
   REAL(ReKi)                 :: RHSsin ( P%DynInflow%MaxInflo+1:maxInfl )


   ! Local Variables:

   INTEGER                    :: i
   INTEGER                    :: NumOut
   INTEGER, PARAMETER         :: UnDyn = 80

   CHARACTER(50)              :: Frmt

   ErrStat = ErrID_None
   ErrMess = ""

!SAVE                                                   ! Save *all* local variables.  Is this necessary, or is OnePass enough.

NumOut = maxInfl+(maxInfl-P%DynInflow%MaxInflo) + 1

IF (m%OnePassDynDbg) THEN

   CALL OpenFOutFile (UnDyn, 'DynDebug.plt', ErrStat, ErrMess)
   IF (ErrStat>=AbortErrLev) RETURN

   Frmt = '( A4,    (: A1, A, I1.1 ) )'

   WRITE(Frmt(7:9), '(I3)') NumOut
   WRITE(UnDyn, Frmt) 'Time',                    &
            ( TAB,    'dAlph_dt',       i,         &
                       i = 1, maxInfl ),         &
            ( TAB,    'dBeta_dt',       i,         &
                     i = P%DynInflow%MaxInflo+1, maxInfl ),  &
              TAB,    'TotalInf'

   m%OnePassDynDbg = .FALSE.

ENDIF

Frmt = '( F10.3,    ( : A1, ES12.5 ) )'

IF (TIME > 0.0D0) THEN

   WRITE(Frmt(10:12), '(I3)') NumOut
   WRITE(UnDyn,Frmt) TIME,                    &
        ( TAB,       m%DynInflow%dAlph_dt(i,1),               &
                   i = 1, maxInfl ),          &
        ( TAB,       m%DynInflow%dBeta_dt(i,1),               &
                   i = P%DynInflow%MaxInflo+1, maxInfl ), &
              TAB,   m%DynInflow%totalInf

ENDIF

RETURN
END SUBROUTINE DynDebug


 ! **********************************************************************
   SUBROUTINE infdist( P, m, ErrStat, ErrMess )
 !  INFDIST calculates the inflow (induced flow) distribution
 !  parameters using Generalized Dynamic Wake Theory.
 ! **********************************************************************
   IMPLICIT                      NONE
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

   ! Local Variables:

   REAL(ReKi)                 :: RHScos ( maxInfl )                     ! The cosine terms go from 1 to 6.  The sine terms go from 3 to 6.
   REAL(ReKi)                 :: RHSsin ( P%DynInflow%MaxInflo+1:maxInfl )          ! The right hand side of the governing equation.
   REAL(ReKi)                 :: tauCos ( maxInfl )                     ! Forcing Functions
   REAL(ReKi)                 :: tauSin ( P%DynInflow%MaxInflo+1:maxInfl )
   REAL(ReKi)                 :: v1
   REAL(ReKi)                 :: v2
   REAL(ReKi)                 :: v3
   REAL(ReKi)                 :: Vplane2

   INTEGER                    :: i
   INTEGER                    :: k
   INTEGER                    :: mode

   ErrStat = ErrID_None
   ErrMess = ""
m%DynInflow%TipSpeed = MAX(P%Blade%r  * m%Rotor%revs, 1.0e-6_ReKi)   !bjj: why is this here AND in InfInit()?

m%DynInflow%Pzero    = pi * p%Wind%rho * m%DynInflow%TipSpeed * m%DynInflow%TipSpeed * P%Blade%r

 ! Non-dimensional time
m%DynInflow%DTO      = m%DT * m%Rotor%REVS

 ! Calculate the wind velocity components in rotor-fixed
 !  coordinates(1-2-3), which are normalized by the tipspeed.
v1 =   m%Wind%VrotorZ / m%DynInflow%TipSpeed     !inplane, upward
v2 =   m%Wind%VrotorY / m%DynInflow%TipSpeed     !inplane, right looking downwind
v3 = - m%Wind%VrotorX / m%DynInflow%TipSpeed     !out-of-plane (normal to the rotor)
 ! Vplane2 is the square of in-plane component of the non-dimensional
 !  wind velocity.
Vplane2  = v1 * v1 + v2 * v2

 ! Note: Direction of non-dimensional velocity (All normal to the rotor plane).
 !  totalInf:  positive downwind. This is the total inflow to the rotor.
 !        v3:        positive upwind (thus, in normal condition, v3 < 0 )
 ! xLambda_M: positive downwind (opposite to A(i,j) )
 !  old_LmdM:  positive downwind (opposite to A(i,j) )
m%DynInflow%totalInf = - v3 + m%DynInflow%old_LmdM
! if ( m%DynInflow%totalInf .le. 0. ) then
!     call usrmes( .true. , &
!        'In SUBROUTINE INFDIST. totalInf =< 0.', &
!        27, 'WARN' )
! endif

 !  VTOTAL is the speed of the total inflow to the rotor.
m%DynInflow%Vtotal = SQRT( m%DynInflow%totalInf * m%DynInflow%totalInf + Vplane2 )
IF (m%DynInflow%vtotal <= 1.0e-6) m%DynInflow%vtotal=1.0e-6

 ! VPARAM is the inflow velocity parameter.
m%DynInflow%Vparam = ( Vplane2 + ( m%DynInflow%totalInf + m%DynInflow%old_LmdM ) * m%DynInflow%totalInf ) / m%DynInflow%Vtotal
 ! Calculate the disk skew angle function
 !  using the effective disk angle of attack.
IF (m%DynInflow% totalInf == 0. ) THEN
!   WRITE(*,*) v3, m%DynInflow%old_LmdM
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Suzuki:
!   xKaiC = 0.
!   xKaiS = 0.
! Shawler:
    m%DynInflow%xKai = 0
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ELSE
 ! Note the definition of Yaw Angle is around -Z axis.
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Suzuki:
!   xKaiC = TAN( .5 * ATAN( -v2 / totalInf ) )
!   xKaiS = TAN( .5 * ATAN(  v1 / totalInf ) )
!!   xKaiC = TAN( .5 * SIGN( ATAN( SQRT( vplane2 ) / totalInf ), v2 ) )
! Shawler:
   m%DynInflow%xkai  = TAN( .5 *       ATAN( SQRT( vplane2 ) / m%DynInflow%totalInf )       )
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ENDIF

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Suzuki:
 ! Update [L_cos] and [L_sin] matrices only if 'xkai' has changed.
 !  Then invert [L_cos] & [L_sin] matrices.
!IF ( xKaiC /= oldKaiC ) THEN
!   CALL LMATRIX ( xKaiC, 1 )
!   CALL MATINV  ( xLcos, xLsin, maxInfl, MaxInflo, 1 )
!ENDIF

!IF ( xKaiS /= oldKaiS ) THEN
!   CALL LMATRIX ( m%DynInflow%xKaiS, 2 )
!   CALL MATINV  ( m%DynInflow%xLcos, m%DynInflow%xLsin, maxInfl, P%DynInflow%MaxInflo, 2 )
!ENDIF
! Shawler:
!IF ( xKai /= oldKai ) THEN
   CALL LMATRIX ( m,p, 1, ErrStat, ErrMess )
      IF (ErrStat /= ErrID_None) THEN
         ErrMess='infdist:'//TRIM(ErrMess)
         RETURN
      END IF   
   CALL LMATRIX ( m,p, 2, ErrStat, ErrMess )
      IF (ErrStat /= ErrID_None) THEN
         ErrMess='infdist:'//TRIM(ErrMess)
         RETURN
      END IF      
   CALL MATINV  ( m%DynInflow%xLcos, m%DynInflow%xLsin, maxInfl, P%DynInflow%MaxInflo, 1, ErrStat, ErrMess )
      IF (ErrStat /= ErrID_None) THEN
         ErrMess='infdist:'//TRIM(ErrMess)
         RETURN
      END IF   
   CALL MATINV  ( m%DynInflow%xLcos, m%DynInflow%xLsin, maxInfl, P%DynInflow%MaxInflo, 2, ErrStat, ErrMess )
      IF (ErrStat /= ErrID_None) THEN
         ErrMess='infdist:'//TRIM(ErrMess)
         RETURN
      END IF   
!ENDIF
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 !  [L_***] is now [L_***]^-1.

 ! Calculate the forcing functions (lift or normal forces).
 !  In the generalized dynamic wake model, the normal forces on the blades
 !  are positve upwind (just like a helicopter rotor). On the other hand,
 !  they are positve downwind in YawDyn, which is normal in wind turbine rotor.
 !  Put a minus sign for each forcing function.

 ! The modes for r=0.
DO mode = 1, P%DynInflow%MaxInflo
   tauCos(mode) = - m%DynInflow%PhiLqC(mode) / m%DynInflow%Pzero * .5
END DO !mode

 ! The modes for r>0.
DO mode = P%DynInflow%MaxInflo+1, maxInfl
   tauCos(mode) = - m%DynInflow%PhiLqC(mode) / m%DynInflow%Pzero
   tauSin(mode) = - m%DynInflow%PhiLqS(mode) / m%DynInflow%Pzero
END DO !mode

 ! Solve for the time derivatives of xAlpha's, {d(alpha)_dt}.
 !  Calculate the right hand side of the governing equation.
 !  {rhs} = 0.5*{tau} - [V][L]^-1*{alpha}

 ! First, calculate {rhs} = V[L]^-1*{alpha}
 !  USE "VTOTAL" for the first row of r=0. Cosine Matrix only.
RHScos(1) = 0.
DO k = 1, maxInfl
   RHScos(1) = RHScos(1) + m%DynInflow%xLcos(1,k) * m%DynInflow%old_Alph(k)
END DO !k
RHScos(1) = .5 * tauCos(1) - m%DynInflow%vtotal * RHScos(1)

 ! USE "VPARAM" for the second row or higher of r=0.
DO i = 2, P%DynInflow%MaxInflo
   RHScos(i) = 0.
   DO k = 1, maxInfl
      RHScos(i) = RHScos(i) + m%DynInflow%xLcos(i,k) * m%DynInflow%old_Alph(k)
   END DO !k
   RHScos(i) = .5 * tauCos(i) - m%DynInflow%Vparam * RHScos(i)
END DO !i

 ! Rows with r=1 or greater. Both cosine and sine matrices.
DO i = P%DynInflow%MaxInflo+1, maxInfl
   RHScos(i) = 0.
   RHSsin(i) = 0.
 ! First, calculate {rhs} = V[L]^-1*{alpha}
   RHScos(i) = sum(m%DynInflow%xLcos(i,1:maxInfl) * m%DynInflow%old_Alph(1:maxInfl))
   RHSsin(i) = sum(m%DynInflow%xLsin(i,P%DynInflow%MaxInflo+1:maxInfl) * m%DynInflow%old_Beta(P%DynInflow%MaxInflo+1:maxInfl))
 ! Second, calculate {rhs} = 0.5*{tau} - [V]{rhs}
 !                         = 0.5*{tau} - [V][L]^-1*{alpha}
 !  USE "VPARAM" for m>0
   RHScos(i) = .5 * tauCos(i) - m%DynInflow%Vparam * RHScos(i)
   RHSsin(i) = .5 * tauSin(i) - m%DynInflow%Vparam * RHSsin(i)
END DO !i

DO i = 1, P%DynInflow%MaxInflo
   m%DynInflow%dAlph_dt(i,1) = p%DynInflow%xMinv(i) * RHScos(i)
END DO !i
DO i = P%DynInflow%MaxInflo+1, maxInfl
   m%DynInflow%dAlph_dt(i,1) = p%DynInflow%xMinv(i) * RHScos(i)
   m%DynInflow%dBeta_dt(i,1) = p%DynInflow%xMinv(i) * RHSsin(i)
END DO !i

 ! Integration using Adams-Bashford predictor corrector method.
 !  Note: The time step is nondimensional. t0=Omega*t
 !  Thus, in YawDyn, DT0=DT*REVS=(DELPSI/REVS)*REVS=DELPSI
 ! USE DT*REVS for compatibility with AeroDyn (DT0 not constant).
 !  DT is in module 'Lift'
 ! Determines xAlpha and xBeta
CALL ABPRECOR ( m%DynInflow%xAlpha, m%DynInflow%old_Alph, m%DynInflow%dAlph_dt, m%DynInflow%DTO, maxInfl, 1 )
CALL ABPRECOR ( m%DynInflow%xBeta,  m%DynInflow%old_Beta, m%DynInflow%dBeta_dt, m%DynInflow%DTO, maxInfl, P%DynInflow%MaxInflo+1 )

 ! Calculate the new lambda_m.
!bjj: why are there 2 do loops? can't they be combined???? (assuming MaxInflo < maxInfl) or is one of these supposed to be sin?
m%DynInflow%xLambda_M= 0.
DO k = 1, P%DynInflow%MaxInflo-1
   m%DynInflow%xLambda_M = m%DynInflow%xLambda_M + m%DynInflow%xLcos(1,k) * m%DynInflow%xAlpha(k)
END DO !k
DO k = P%DynInflow%MaxInflo, maxInfl
   m%DynInflow%xLambda_M = m%DynInflow%xLambda_M + m%DynInflow%xLcos(1,k) * m%DynInflow%xAlpha(k)
END DO !k
! xLambda_M = 2. / sqrt(3.) * xLambda_M
m%DynInflow%xLambda_M = 1.1547005 * m%DynInflow%xLambda_M

! Added additional output for GDW debugging - comment out for distribution
!CALL DynDebug (Time, P, x, xd, z, m, y, ErrStat, ErrMess, RHScos, RHSsin)

RETURN
END SUBROUTINE infdist

 ! *************************************************************
   SUBROUTINE VINDINF( P, m, ErrStat, ErrMess, &
                       iradius, iblade, rlocal, vnw, VNB, VT, psi )
 ! Calculates the axial and tangential induction factor for each annular segment
 ! and time step (i.e. sets m%Element%A and m%Element%AP)
 ! Uses the calculated inflow parameters
 !  Called by ElemFrc for each element at a new time step.
 ! *************************************************************

   IMPLICIT                      NONE
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess


   ! Passed Variables:
   REAL(ReKi),INTENT(IN)      :: psi
   REAL(ReKi),INTENT(IN)      :: rlocal
   REAL(ReKi),INTENT(IN)      :: VNB
   REAL(ReKi),INTENT(IN)      :: vnw
   REAL(ReKi),INTENT(IN   )   :: VT     ! Tangential velocity from relative blade motion and wind, no induction

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
   
   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.
   CHARACTER(ErrMsgLen)                      :: ErrMessLcl          ! Error message returned by called routines.


ErrStat = ErrID_None
ErrMess = ""

 ! Rzero is the non-dimensional radius.
Rzero  = rlocal / P%Blade%r

 ! PSIBAR is the azimuth angle measured counterclockwise from
 !  the horizontal to the left looking downwind.  The directions
 !  of PSI and PSIBAR are opposite and there is 90 deg difference.
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Suzuki:
! psiBar = - psi - piBy2
! Shawler:
CALL WindAzimuthZero (psi,m%Wind%VrotorY,m%Wind%VrotorZ,WindPsi)
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ! Calculate the induction factor using the inflow parameters.
 !  Here it's normalized by the tipspeed.

m%Element%A(iRadius,iBlade) =  0.

DO mode = 1, p%DynInflow%MaxInflo
   m%Element%A(iRadius,iBlade) = m%Element%A(iRadius,iBlade)  &
                     + xphi(Rzero,mode,ErrStatLcl, ErrMessLcl) * m%DynInflow%xAlpha(mode)
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'VINDINF' )
      IF (ErrStat >= AbortErrLev) RETURN
   
!  &  + phis(Rzero, MRvector(mode), NJvector(mode) )* xAlpha(mode)
END DO !mode

 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!Suzuki:
!DO mode = MaxInflo+1, maxInfl
!   A(iRadius,iBlade) = A(iRadius,iBlade) + xphi(Rzero,mode) *      &
!!  &    + phis(Rzero, MRvector(mode), NJvector(mode) ) *
!            ( xAlpha(mode) * COS( REAL(MRvector(MODE), ReKi) * psibar )  &
!            + xBeta (mode) * SIN( REAL(MRvector(MODE), ReKi) * psibar ) )
!END DO !mode
! Shawler:
DO mode = p%DynInflow%MaxInflo+1, maxInfl
   m%Element%A(iRadius,iBlade) = m%Element%A(iRadius,iBlade) + xphi(Rzero,mode,ErrStatLcl, ErrMessLcl) *      &
            ( m%DynInflow%xAlpha(mode) * COS( REAL(MRvector(MODE), ReKi) * Windpsi )  &
            + m%DynInflow%xBeta (mode) * SIN( REAL(MRvector(MODE), ReKi) * Windpsi ) )
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'VINDINF' )
      IF (ErrStat >= AbortErrLev) RETURN
   
END DO !mode
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 ! A is positive upwind, while alpha's & beta's are positve downwind.
 !  Also, they are normalized by different factors.
m%Element%A(iRadius,iBlade) = - m%Element%A(iRadius,iBlade) * m%DynInflow%TipSpeed / VNW

 ! Calculate induced swirl (a') if desired.

m%Element%AP(iRadius,iBlade) = 0.0_ReKi ! Default value

IF ( P%SWIRL ) THEN
 ! akihiro 10/26/99
   SWRLARG = 1.0 + 4.0 *  m%Element%A(iradius,iblade) * VNW * &
           ( (1.0 - m%Element%A(iradius,iblade)) * VNW + VNB ) / VT / VT
   IF ( SWRLARG > 0.0 ) THEN
      A2P = 0.5 * ( -1.0 + SQRT( SWRLARG ) )
! bjj: this value was not properly set before. We could also just replace the local A2P variable with AP() instead.
      m%Element%AP(iRadius,iBlade) = A2P
   ENDIF

ENDIF

RETURN
END SUBROUTINE VINDINF

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
   SUBROUTINE LMATRIX ( m, p, matrixMode,ErrStat, ErrMess)
 !  LMATRIX calculates the [L_***] matrix using Gamma matrix and xkai=X.
 !   matrixMode = 1 : Calculate [L_cos] matrix
 !   matrixMode = 2 : Calculate [L_sin] matrix
 ! ***********************************************************************

IMPLICIT                      NONE
   ! Passed Variables:
   TYPE(AD14_MiscVarType),      INTENT(INOUT)  :: m           ! Misc/optimization variables
TYPE(AD14_ParameterType),       INTENT(IN)     :: p           ! Parameters
INTEGER   ,INTENT(IN)      :: matrixMode

INTEGER(IntKi),   INTENT(OUT)  :: ErrStat      ! Error status of the operation
CHARACTER(*),     INTENT(OUT)  :: ErrMess      ! Error message if ErrStat /= ErrID_None


REAL(ReKi)      :: X 


   ! Local Variables:

REAL(ReKi)                 :: XM

INTEGER                    :: IROW
INTEGER                    :: JCOL

   ErrStat = ErrID_None
   ErrMess = ""


X = m%DynInflow%xKai

 ! Check the value of X
IF ( ( X < -1.) .OR. ( X > 1.) ) THEN
!   IF ( ( X < 0.) .OR. ( X > 1.) ) THEN
   ErrStat = ErrID_Fatal
   ErrMess = 'LMATRIX: Value of X = '//TRIM(Num2LStr(X))//' must be between -1 and 1.'
   RETURN
ENDIF

SELECT CASE ( matrixMode )

 ! Calculate the rows for r=0 of [L_cos] matrix,
 !  which needs a separate formula.
   CASE (1)
    DO JCOL = 1, maxInfl
       XM = X ** MRvector(JCOL)
       DO IROW = 1, p%DynInflow%MaxInflo
          m%DynInflow%xLcos( IROW, JCOL ) = m%DynInflow%gamma( IROW, JCOL ) * XM
       END DO !IROW
    END DO !JCOL

 ! Calculate the [L_cos] matrix,
 !  the rows for r=1 and higher.
    DO IROW = p%DynInflow%MaxInflo+1, maxInfl
       DO JCOL = 1, maxInfl
          m%DynInflow%xLcos( IROW, JCOL ) = m%DynInflow%GAMMA( IROW, JCOL )  &
                  *(    X  ** m%DynInflow%MminusR( IROW, JCOL )  &
                      + X  ** m%DynInflow%MplusR ( IROW, JCOL )  &
                  *  (-1.) ** m%DynInflow%MminR  ( IROW, JCOL )  )
       END DO !JCOL
    END DO !IROW

 ! Calculate [L_sin] matrix.
   CASE (2)
    DO IROW = p%DynInflow%MaxInflo+1, maxInfl
       DO JCOL = p%DynInflow%MaxInflo+1, maxInfl
          m%DynInflow%xLsin( IROW, JCOL ) = m%DynInflow%GAMMA( IROW, JCOL )  &
                  *(    X  ** m%DynInflow%MminusR( IROW, JCOL )  &
                      - X  ** m%DynInflow%MplusR ( IROW, JCOL )  &
                  *  (-1.) ** m%DynInflow%MminR  ( IROW, JCOL )  )
      END DO !JCOL
    END DO !IROW

   CASE DEFAULT
    CALL ProgAbort( 'Value of matrixMode = '//TRIM(Int2LStr(matrixMode))//' must be 1 or 2.')
END SELECT

RETURN
END SUBROUTINE LMATRIX

 ! ***********************************************************************
   SUBROUTINE MATINV( A0, A1, N, N0, invMode,ErrStat,ErrMsg )
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

INTEGER(IntKi),   INTENT(OUT)  :: ErrStat      ! Error status of the operation
CHARACTER(*),     INTENT(OUT)  :: ErrMsg       ! Error message if ErrStat /= ErrID_None


   ! Local Variables:

REAL(ReKi)                 :: DUMMY(      N-N0,      N-N0 )

INTEGER                    :: I
INTEGER                    :: J

   ErrStat = ErrID_None
   ErrMsg  = ""


 ! Invert [L_cos] in all cases
 !  [L_cos] matrix can be inverted without a dummy.
 !  Invert the [A0]=[L_cos] matrix by Gauss-Jordan Method
SELECT CASE ( invMode )

   CASE (1)
    CALL GAUSSJ(A0,N,ErrStat,ErrMsg)
    IF (ErrStat /= ErrID_None) THEN
      ErrMsg = 'MATINV:'//TRIM(ErrMsg)
      RETURN
    END IF

 ! [L_sin] matrix needs a dummy array, because the index goes
 !  from MaxInflo(=N0) to maxInfl(=N),
 !  which is incompatible with SUBROUTINE GAUSSJ.
 !BJJ: IS THIS REALLY NECESSARY?  Aren't the indicies an abstraction?? if you pass an array (1:3) and use it as (0:2) in another subroutine, it's really okay???
   CASE (2)
    DO I=1,N-N0
       DO J=1,N-N0
          DUMMY(I,J) = A1(I+N0,J+N0)
       END DO !J
    END DO !I

 ! Invert the [A1]=[L_sin] matrix by Gauss-Jordan Method.
    CALL GAUSSJ(DUMMY,N-N0,ErrStat,ErrMsg)
    IF (ErrStat /= ErrID_None) THEN
      ErrMsg = 'MATINV:'//TRIM(ErrMsg)
      RETURN
    END IF

 ! Put the dummy back into [A1]=[L_sin].
    DO I=1,N-N0
       DO J=1,N-N0
          A1(I+N0,J+N0) = DUMMY(I,J)
       END DO !J
    END DO !I

   CASE DEFAULT
   CALL ProgAbort( 'Value of invMode = '//TRIM(Int2LStr(invMode))//' must be 1 or 2.')
    IF (ErrStat >= AbortErrLev) RETURN    
   

END SELECT



RETURN
END SUBROUTINE MATINV

 ! ***********************************************************************
   SUBROUTINE gaussj(a,n, ErrStat, ErrMsg)
 !  Invert a matrix by Gauss-Jordan Method. The original source code
 !   from "Numerical Recipe" was slightly modified.
 ! ***********************************************************************


IMPLICIT                      NONE


   ! Passed Variables:
INTEGER   ,    INTENT(IN)      :: n
REAL(ReKi),    INTENT(INOUT)   :: a(n,n)

INTEGER(IntKi),   INTENT(OUT)  :: ErrStat      ! Error status of the operation
CHARACTER(*),     INTENT(OUT)  :: ErrMsg       ! Error message if ErrStat /= ErrID_None


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


ErrStat=ErrID_None
ErrMsg=""

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
               ErrStat = ErrID_Fatal
               ErrMsg  = "gaussj: Singular matrix encountered." 
               RETURN
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
      ErrStat = ErrID_Fatal
      ErrMsg  = "gaussj: Singular matrix encountered." 
      RETURN
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
          * SQRT( REAL( (2*N+1) * (2*J+1), ReKi ) ) &
          / SQRT( HFUNC(M,N) * HFUNC(R,J) )   &
          / REAL( (J+N) * (J+N+2) * ((J-N)*(J-N)-1), ReKi )

ELSE IF ( ABS(J-N) == 1 ) THEN  !bjj: why don't we use the pi() variable? or PibyTwo
   FGAMMA = Pi * SIGN(1.0_ReKi, REAL(R-M, ReKi) ) * .5 &
          / SQRT( HFUNC(M,N) * HFUNC(R,J) )        &
          / SQRT( REAL( (2*N+1) * (2*J+1) , ReKi) )

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

HFUNC = ( REAL( IDUBFACT(NPM-1), ReKi ) / REAL( IDUBFACT(NPM), ReKi ) ) &
      * ( REAL( IDUBFACT(NMM-1), ReKi ) / REAL( IDUBFACT(NMM), ReKi ) )



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
   FUNCTION xphi( Rzero, mode,ErrStat, ErrMsg )
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
INTEGER(IntKi),   INTENT(OUT)  :: ErrStat      ! Error status of the operation
CHARACTER(*),     INTENT(OUT)  :: ErrMsg       ! Error message if ErrStat /= ErrID_None




IF ( Rzero < 0. ) THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Value of Rzero = '//TRIM(Num2LStr(Rzero))//' must be larger than 0 in xphi().'
   RETURN
ELSE IF ( Rzero > 1. ) THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Value of Rzero = '//TRIM(Num2LStr(Rzero))//' must be smaller than 1 in xphi().'
   RETURN
ELSE
   ErrStat = ErrID_None
   ErrMsg  = ''   
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
   CALL ProgAbort('Value of Rzero = '//TRIM(Num2LStr(Rzero))//' must be larger than 0 in phis().' )
ELSE IF ( Rzero > 1. ) THEN
   CALL ProgAbort('Value of Rzero = '//TRIM(Num2LStr(Rzero))//' must be smaller than 1 in phis().' )
ENDIF

phis = 0.

DO q = r, j-1, 2
   phis = phis  &
        + Rzero ** q * (-1.) **((q-r)/2) * REAL( idubfact(j+q), ReKi ) &
        / REAL( idubfact(q-r) * idubfact(q+r) * idubfact(j-q-1), ReKi )
END DO !q

phis = phis * SQRT( REAL( 2*j+1, ReKi ) * hfunc(r,j) )


RETURN
END FUNCTION phis


!***********************************************************
SUBROUTINE WindAzimuthZero (psi,VrotorY,VrotorZ,WindPsi)
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


IMPLICIT                      NONE


   ! Passed Variables:
REAL(ReKi),INTENT(IN)      :: psi,VrotorY,VrotorZ
REAL(ReKi),INTENT(OUT)     :: WindPsi

WindPsi = psi - ATAN2(VrotorY,-VrotorZ)

RETURN
END SUBROUTINE WindAzimuthZero

 !     ************* END OF FILE ***************

! ********************************************************************
!====================================================================================================
SUBROUTINE CheckRComp( P, x, xd, z, m, y, ErrStat, ErrMess, &
                       ADFile, HubRadius, TipRadius )
! This routine checks to see if RElm(:) and DR(:) are compatible within a millimeter;
!----------------------------------------------------------------------------------------------------
   IMPLICIT                      NONE
   TYPE(AD14_ParameterType),       INTENT(IN)  :: p           ! Parameters
   TYPE(AD14_ContinuousStateType), INTENT(INOUT)  :: x           ! Initial continuous states
   TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(AD14_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD14_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables
   TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess

      ! Passed variables

   CHARACTER(*), INTENT(IN)   :: ADFile                     ! Name of the AeroDyn input file, used for printing error message
   REAL(ReKi),   INTENT(IN)   :: HubRadius                  ! Hub radius, used to verify that RElm and DR are input correctly
   REAL(ReKi),   INTENT(IN)   :: TipRadius                  ! Tip radius, used to verify that RElm and DR are input correctly


      ! Local variables.

   REAL(ReKi)                 :: DRNodesNew(P%Element%NElm) ! Length of variable-spaced blade elements--calculated from input RElm(:).
   REAL(ReKi)                 :: DRSum                      ! Sum of DRs--should be close to TipRadius
   REAL(ReKi), PARAMETER      :: EPS = EPSILON(HubRadius)   ! A small value used to compare two numbers

   INTEGER                    :: I                          ! Generic index.

   CHARACTER(33)              :: DRChange                   ! A string showing how to change DR to get campatibility



      ! Initialize ErrStat to 0 (no error):

   ErrStat = ErrID_None
   ErrMess = ""


      ! Calculate DRNodesNew(:) based on input RElm(:) and compare to input DR(:):

!   AssumedHubRadius = P%Element%RElm(1) - 0.5*P%Blade%DR(1)
   DRNodesNew(1)    = 2.0*( P%Element%RElm(1) - HubRadius )
   DRSum            = DRNodesNew(1) + HubRadius

   IF ( DRNodesNew(1) <= EPS )  THEN                              ! Check to see if RElm(1) > HubRad; if not, ProgAbort program
      ErrMess = 'CheckRComp: RElm(1) must be larger than the hub radius (HubRadius = RElm(1) - 0.5*DR(1)). '
      ErrStat = ErrID_Fatal
      RETURN
   ELSEIF ( ABS( DRNodesNew(1) - P%Blade%DR(1) ) > 0.001 )  THEN     ! Check to see if the calculated DRNodes(1) is close to the inputted DRNodes(1); if not, set flag--this will cause the program to ProgAbort later.
      !ErrMess = ' this error message will be written at the end of the routine '
      ErrStat = ErrID_Fatal
   ENDIF

   DO I = 2,P%Element%NElm ! Loop through all but the innermost blade element

      DRNodesNew(I) = 2.0*( P%Element%RElm(I) - P%Element%RElm(I-1) ) - DRNodesNew(I-1)
      DRSum         = DRSum + DRNodesNew(I)


      IF ( DRNodesNew(I) <= EPS )  THEN                           ! Check to see if it is possible to have compatible DR(:) with the input RElm(:); if not, abort program
         ErrMess = 'CheckRComp: RElm('//TRIM( Int2LStr(I) )//') produces ill-conditioned DR(:)' 
         ErrStat = ErrID_Fatal
         RETURN
      ELSEIF ( ABS( DRNodesNew(I) - P%Blade%DR(I) ) > 0.001 )  THEN  ! Check to see if the calculated DRNodes(I) is close to the inputted DRNodes(I); if not, set flag--this will cause the program to Abort later.
         !ErrMess = ' this error message will be written at the end of the routine '
         ErrStat = ErrID_Fatal
      ENDIF

   END DO             ! I - all but the innermost blade element


      ! Abort program if necessary

   IF ( ErrStat /= ErrID_None )  THEN

         ! Write error message since the input DR(:) are not close to the calculated DRNodesNew(:)

      ErrMess = ' Input values for DR(:) are not compatible with input RElm(:).'
      CALL WrScr1(TRIM(ErrMess))
      CALL WrScr( ' To make them compatible, please modify DR in the AeroDyn input file, '// TRIM( ADFile ) //', as follows:' )
      CALL WrScr1(' DR (Old) --> DR (New) ')

      DO I = 1,P%Element%NElm
         WRITE( DRChange, "(' ', F13.4, ' --> ', F13.4, ' ')" ) P%Blade%DR(I), DRNodesNew(I)
         CALL WrScr( DRChange )
      ENDDO !I

      ErrMess = 'CheckRComp: '//TRIM(ErrMess)
      RETURN

   ELSEIF ( ABS( DRSum - TipRadius ) > 0.001 )  THEN

         ! Abort program since SUM( DRNodes(:) ) /= ( TipRadius - HubRadius)

      ErrMess = 'CheckRComp: TipRadius must be equal to HubRadius + SUM( P%Blade%DR(:) ).'
      ErrStat = ErrID_Fatal

      RETURN

   ENDIF


   RETURN
END SUBROUTINE CheckRComp

                             
END MODULE AeroSubs
