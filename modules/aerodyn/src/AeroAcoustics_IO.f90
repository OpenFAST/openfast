!**********************************************************************************************************************************
! File last committed: 2020-02-12
!**********************************************************************************************************************************
MODULE AeroAcoustics_IO

   use NWTC_Library
   use AeroAcoustics_Types

   implicit none

   type(ProgDesc), parameter  :: AA_Ver = ProgDesc( 'AeroAcoustics', 'v1.00.00', '18-Aug-2016' )
   character(*),   parameter  :: AA_Nickname = 'AA'



   INTEGER(IntKi), PARAMETER      :: Time      =    0

     ! Parameters related to output length (number of characters allowed in the output data headers):

   INTEGER(IntKi), PARAMETER      :: OutStrLenM1 = ChanLen - 1

   INTEGER(IntKi), PARAMETER        :: MaxBl    =  3                                   ! Maximum number of blades allowed in simulation

   ! model identifiers
   integer(intKi), parameter        :: ModelUnknown  = -1

! FLAG TO COMPUTE BLUNTNESS NOISE           = 0 No, =1 Yes
   integer(intKi), parameter        :: IBLUNT_None  = 0
   integer(intKi), parameter        :: IBLUNT_BPM  = 1

! FLAG TO COMPUTE Laminar Boundary Layer Noise          = 0 No, =1 Yes
   integer(intKi), parameter        :: ILAM_None      = 0  ! steady model
   integer(intKi), parameter        :: ILAM_BPM       = 1  !

! FLAG TO COMPUTE Tip  Noise          = 0 No, =1 Yes
   integer(intKi), parameter        :: ITIP_None  = 0  !
   integer(intKi), parameter        :: ITIP_On      = 1  !

   integer(intKi), parameter        :: ITRIP_None     = 0  ! not tripped boundary layer
   integer(intKi), parameter        :: ITRIP_Heavy    = 1  ! heavily tripped boundary layer
   integer(intKi), parameter        :: ITRIP_Light      = 2  ! light tripped boundary layer

! calculation method for boundary layer properties,  = 1 BPM = 2 Pretabulated BL values
   integer(intKi), parameter        :: X_BLMethod_BPM     = 1  !
   integer(intKi), parameter        :: X_BLMethod_Tables  = 2  !

   integer(intKi), parameter        :: TICalc_Interp  = 1  ! interpolate from pretabulated
   integer(intKi), parameter        :: TICalc_Every   = 2  ! calculate ti automatically

   integer(intKi), parameter        :: ITURB_None           = 0  ! TBLTE noise is not calculated
   integer(intKi), parameter        :: ITURB_BPM            = 1  ! TBLTE noise is calculated with BPM
   integer(intKi), parameter        :: ITURB_TNO            = 2  ! TBLTE noise is calculated with TNO

   integer(intKi), parameter        :: IInflow_None                 = 0  ! IInflow noise is not calculated
   integer(intKi), parameter        :: IInflow_BPM                  = 1  ! IInflow noise is calculated with BPM
   integer(intKi), parameter        :: IInflow_FullGuidati         = 2  ! IInflow noise is calculated with FullGuidati
   integer(intKi), parameter        :: IInflow_SimpleGuidati    = 3  ! IInflow noise is calculated with SimpleGuidati

contains
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadInputFiles( InputFileName, BL_Files, InputFileData, Default_DT, OutFileRoot, NumBlades, UnEcho, ErrStat, ErrMsg )
    ! This subroutine reads the input file and stores all the data in the AA_InputFile structure.
    ! It does not perform data validation.
    !..................................................................................................................................
    ! Passed variables
    REAL(DbKi),              INTENT(IN)    :: Default_DT      ! The default DT (from glue code)
    CHARACTER(*),            INTENT(IN)    :: InputFileName   ! Name of the aeroacoustics input file
    CHARACTER(*), dimension(:),            INTENT(IN)    :: BL_Files         ! Name of the BL input file
    CHARACTER(*),            INTENT(IN)    :: OutFileRoot     ! The rootname of all the output files written by this routine.
    TYPE(AA_InputFile),      INTENT(OUT)   :: InputFileData   ! Data stored in the module's input file
    INTEGER(IntKi),          INTENT(OUT)   :: UnEcho          ! Unit number for the echo file
    INTEGER(IntKi),          INTENT(IN)    :: NumBlades       ! Number of blades for this model
    INTEGER(IntKi),          INTENT(OUT)   :: ErrStat         ! The error status code
    CHARACTER(*),            INTENT(OUT)   :: ErrMsg          ! The error message, if an error occurred
    ! local variables
    INTEGER(IntKi)                         :: I
    INTEGER(IntKi)                         :: ErrStat2        ! The error status code
    CHARACTER(ErrMsgLen)                   :: ErrMsg2         ! The error message, if an error occurred
    CHARACTER(*), PARAMETER                :: RoutineName = 'ReadInputFiles'
    ! initialize values:
    ErrStat = ErrID_None
    ErrMsg  = ''
    UnEcho  = -1


    ! Reads the module input-file data
    CALL ReadPrimaryFile( InputFileName, InputFileData, Default_DT,  OutFileRoot, UnEcho, ErrStat2, ErrMsg2 )
    if(Failed()) return

    ! get the blade input-file data
    ALLOCATE( InputFileData%BladeProps( size(BL_Files) ), STAT = ErrStat2 )
    IF (ErrStat2 /= 0) THEN
        CALL SetErrStat(ErrID_Fatal,"Error allocating memory for BladeProps.", ErrStat, ErrMsg, RoutineName)
        return
    END IF

    if ((InputFileData%ITURB==2) .or. (InputFileData%X_BLMethod==2) .or. (InputFileData%IBLUNT==1)) then
        ! We need to read the BL tables
        CALL ReadBLTables( InputFileName, BL_Files, InputFileData, ErrStat2, ErrMsg2 )
        if (Failed())return
    endif

    IF(   (InputFileData%TICalcMeth==1) ) THEN
        CALL REadTICalcTables(InputFileName,InputFileData,  ErrStat2, ErrMsg2); if(Failed()) return
    ENDIF

CONTAINS
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
    end function Failed

END SUBROUTINE ReadInputFiles
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads in the primary Noise input file and places the values it reads in the InputFileData structure.
!   It opens and prints to an echo file if requested.
SUBROUTINE ReadPrimaryFile( InputFile, InputFileData, Default_DT, OutFileRoot, UnEc, ErrStat, ErrMsg )
    integer(IntKi),     intent(out)     :: UnEc                                ! I/O unit for echo file. If > 0, file is open for writing.
    integer(IntKi),     intent(out)     :: ErrStat                             ! Error status
    REAL(DbKi),         INTENT(IN)      :: Default_DT                          ! The default DT (from glue code)
    character(*),       intent(in)      :: InputFile                           ! Name of the file containing the primary input data
    character(*),       intent(out)     :: ErrMsg                              ! Error message
    character(*),       intent(in)      :: OutFileRoot                         ! The rootname of the echo file, possibly opened in this routine
    type(AA_InputFile), intent(inout)   :: InputFileData                       ! All the data in the Noise input file
    ! Local variables:
    real(ReKi)                    :: TmpAry(3)                                 ! array to help read tower properties table
    integer(IntKi)                :: I                                         ! loop counter
    integer(IntKi)                :: UnIn,UnIn2                                ! Unit number for reading file
    integer(IntKi)                :: loop1                                     ! loop counter
    character(1024)               :: ObserverFile                              ! name of the files containing obesever location
    integer(IntKi)                :: ErrStat2, IOS,cou                             ! Temporary Error status
    logical                       :: Echo                                      ! Determines if an echo file should be written
    character(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
    character(1024)               :: PriPath                                   ! Path name of the primary file
    character(200)                :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
    character(*), parameter       :: RoutineName = 'ReadPrimaryFile'
    integer(IntKi)                :: n                                         ! dummy integer
    ! Initialize some variables:
    ErrStat = ErrID_None
    ErrMsg  = ""

    UnEc = -1
    Echo = .FALSE.
    CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.

    ! Open the Primary input file.
    CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 ); call check
    CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 ); call check
    IF ( ErrStat >= AbortErrLev ) THEN
        CALL Cleanup()
        RETURN
    END IF

    ! Read the lines up/including to the "Echo" simulation control variable
    ! If echo is FALSE, don't write these lines to the echo file.
    ! If Echo is TRUE, rewind and write on the second try.
    I = 1 !set the number of times we've read the file
    DO
        !----------- HEADER -------------------------------------------------------------
        CALL ReadCom( UnIn, InputFile, 'File header: Module Version (line 1)', ErrStat2, ErrMsg2, UnEc ); call check
        CALL ReadStr( UnIn, InputFile, InputFileData%FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc ); call check
        IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
        END IF

        !----------- GENERAL OPTIONS ----------------------------------------------------
        CALL ReadCom( UnIn, InputFile, 'Section Header: General Options', ErrStat2, ErrMsg2, UnEc ); call check
        ! Echo - Echo input to "<RootName>.AD.ech".
        CALL ReadVar( UnIn, InputFile, Echo, 'Echo',   'Echo flag', ErrStat2, ErrMsg2, UnEc); call check
        IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop
        ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
        I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
        CALL OpenEcho ( UnEc, TRIM(OutFileRoot)//'.ech', ErrStat2, ErrMsg2, AA_Ver ); call check
        IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
        END IF
        IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(AA_Ver%Name)//' primary input file "'//TRIM( InputFile )//'":'
        REWIND( UnIn, IOSTAT=ErrStat2 )
        IF (ErrStat2 /= 0_IntKi ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
        END IF
    END DO

    IF (NWTC_VerboseLevel == NWTC_Verbose) THEN
        CALL WrScr( ' Heading of the '//TRIM(AA_Ver%Name)//' input file: ' )
        CALL WrScr( '   '//TRIM( InputFileData%FTitle ) )
    END IF

    ! DT_AA - Time interval for aerodynamic calculations {or default} (s):
    Line = ""
    CALL ReadVar( UnIn, InputFile, Line, "DT_AA", "Time interval for aeroacoustics calculations {or default} (s)", ErrStat2, ErrMsg2, UnEc); call check
    CALL Conv2UC( Line )

    IF ( INDEX(Line, "DEFAULT" ) /= 1 ) THEN ! If DT_AA is not "default", read it and make sure it is a multiple of DTAero from AeroDyn. Else, just use DTAero
        READ( Line, *, IOSTAT=IOS) InputFileData%DT_AA
        CALL CheckIOS ( IOS, InputFile, 'DT_AA', NumType, ErrStat2, ErrMsg2 ); call check

        IF (abs(InputFileData%DT_AA / Default_DT - NINT(InputFileData%DT_AA / Default_DT)) .gt. 1E-10) THEN
            CALL SetErrStat(ErrID_Fatal,"The Aeroacoustics input DT_AA must be a multiple of DTAero.", ErrStat, ErrMsg, RoutineName)
            return
        END IF
    ELSE
        InputFileData%DT_AA = Default_DT
    END IF

    CALL ReadVar(UnIn,InputFile,InputFileData%AAStart      ,"AAStart"      ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%AA_Bl_Prcntge,"BldPrcnt"     ,"-",ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadCom( UnIn, InputFile, 'Section Header: Aeroacoustic Models', ErrStat2, ErrMsg2, UnEc ); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%IInflow      ,"InflowMod"    ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%TICalcMeth   ,"TICalcMeth"   ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVAr(UnIn,InputFile,InputFileData%TICalcTabFile,"TICalcTabFile","" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%Lturb        ,"Lturb"        ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%ITURB        ,"TurbMod"      ,"" ,ErrStat2,ErrMsg2,UnEc); call check ! ITURB - TBLTE NOISE
    CALL ReadVar(UnIn,InputFile,InputFileData%X_BLMethod   ,"BLMod"        ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%ITRIP        ,"TripMod"      ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%ILAM         ,"LamMod"       ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%ITIP         ,"TipMod"       ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%ROUND        ,"RoundTip"     ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%ALPRAT       ,"ALPRAT"       ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar(UnIn,InputFile,InputFileData%IBLUNT       ,"BluntMod"     ,"" ,ErrStat2,ErrMsg2,UnEc); call check

    ! Return on error at end of section
    IF ( ErrStat >= AbortErrLev ) THEN
        CALL Cleanup()
        RETURN
    END IF

    !----------- OBSERVER INPUT  ------------------------------
    CALL ReadCom( UnIn, InputFile, 'Section Header: Observer Input ', ErrStat2, ErrMsg2, UnEc ); call check
    !----- read from observer file
    CALL ReadVar ( UnIn, InputFile, ObserverFile, ObserverFile, 'Name of file  observer locations', ErrStat2, ErrMsg2, UnEc ); call check
    IF ( PathIsRelative( ObserverFile ) ) ObserverFile = TRIM(PriPath)//TRIM(ObserverFile)

    CALL GetNewUnit( UnIn2, ErrStat2, ErrMsg2 ); call check

    CALL OpenFInpFile ( UnIn2, ObserverFile, ErrStat2, ErrMsg2 ); call check
    IF ( ErrStat >= AbortErrLev ) RETURN
    
    ! NrObsLoc  - Nr of Observers (-):
    CALL ReadVar( UnIn2, ObserverFile, InputFileData%NrObsLoc, "NrObsLoc", "Nr of Observers (-)", ErrStat2, ErrMsg2, UnEc); call check

    ! Observer location in tower-base coordinate  (m):
    CALL AllocAry( InputFileData%ObsX,InputFileData%NrObsLoc, 'ObsX', ErrStat2, ErrMsg2); call check
    CALL AllocAry( InputFileData%ObsY,InputFileData%NrObsLoc, 'ObsY', ErrStat2, ErrMsg2); call check
    CALL AllocAry( InputFileData%ObsZ,InputFileData%NrObsLoc, 'ObsZ', ErrStat2, ErrMsg2); call check

    CALL ReadCom( UnIn2, InputFile, ' Header', ErrStat2, ErrMsg2, UnEc ); call check

    DO cou=1,InputFileData%NrObsLoc
        READ( UnIn2, *, IOStat=IOS )  InputFileData%ObsX(cou), InputFileData%ObsY(cou), InputFileData%ObsZ(cou)
        CALL CheckIOS( IOS, ObserverFile, 'Obeserver Locations '//TRIM(Num2LStr(cou)), NumType, ErrStat2, ErrMsg2 ); call check
        ! Return on error if we couldn't read this line
        IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
        END IF
    ENDDO
    CLOSE ( UnIn2 )
    !----- end read from observer file

    !----------- OUTPUTS  -----------------------------------------------------------
    CALL ReadCom( UnIn, InputFile, 'Section Header: Outputs', ErrStat2, ErrMsg2, UnEc); call check
    CALL ReadVar( UnIn,InputFile,InputFileData%aweightflag  ,"AWeighting"   ,"" ,ErrStat2,ErrMsg2,UnEc); call check
    CALL ReadVar( UnIn, InputFile, InputFileData%NrOutFile, "NrOutFile", "Nr of Output Files (-)", ErrStat2, ErrMsg2, UnEc); call check
    CALL AllocAry( InputFileData%AAOutFile,InputFileData%NrOutFile, 'AAOutFile', ErrStat2, ErrMsg2); call check
    CALL ReadVar ( UnIn, InputFile, InputFileData%AAOutFile(1), 'AAOutFile', 'Name of output file ', ErrStat2, ErrMsg2, UnEc ); call check
    DO I=InputFileData%NrOutFile,1,-1
        ! one file name is given by the user and the XXFile1.out XXFile2.out XXFile3.out is generated
        IF ( PathIsRelative( InputFileData%AAOutFile(I) ) ) InputFileData%AAOutFile(I) = TRIM(PriPath)//TRIM(InputFileData%AAOutFile(1))//TRIM(Num2Lstr(I))//".out"
    ENDDO

    ! Return on error at end of section
    IF ( ErrStat >= AbortErrLev ) THEN
        CALL Cleanup()
        RETURN
    END IF
    !---------------------- END OF FILE -----------------------------------------
    CALL Cleanup( )

CONTAINS
   SUBROUTINE Check()
       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END SUBROUTINE Check

   !...............................................................................................................................
   SUBROUTINE Cleanup()
       IF (UnIn > 0) CLOSE ( UnIn )
   END SUBROUTINE Cleanup
   !...............................................................................................................................
END SUBROUTINE ReadPrimaryFile
!----------------------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

subroutine ReadRealMatrix(fid, FileName, Mat, VarName, nLines,nRows, iStat, Msg, iLine )
    integer, intent(in)                     :: fid
    real(DbKi), dimension(:,:), allocatable :: Mat
    character(len=*), intent(in)            :: FileName
    character(len=*), intent(in)            :: VarName
    integer, intent(in)                     :: nLines
    integer, intent(in)                     :: nRows
    integer, intent(out)                    :: iStat
    integer, intent(inout)                  :: iLine
    character(len=*), intent(inout)         :: Msg
    ! local variables
    integer :: i
    if (allocated(Mat)) deallocate(Mat)
    call allocAry( Mat, nLines, nRows, VarName,  iStat, Msg);
    if (iStat /= 0) return
    !Read Stiffness
    DO I =1,nLines
       iLine=iLine+1
       ! TODO use ReadCAryFromStr when available in the NWTCIO, it performs more checks
       CALL ReadAry( fid, FileName, Mat(I,:), nRows, trim(VarName)//' Line '//Num2LStr(iLine), VarName, iStat, Msg) ! From NWTC_Library
       if (iStat /= 0) return
    ENDDO
end subroutine



SUBROUTINE ReadBLTables( InputFile, BL_Files, InputFileData, ErrStat, ErrMsg )
    ! Passed variables
    character(*),       intent(in)      :: InputFile                           ! Name of the file containing the primary input data
    character(*), dimension(:),       intent(in)      :: BL_Files                           ! Name of the file containing the primary input data
type(AA_InputFile), intent(inout)       :: InputFileData                       ! All the data in the Noise input file
    integer(IntKi),     intent(out)     :: ErrStat                             ! Error status
    character(*),       intent(out)     :: ErrMsg                              ! Error message
    ! Local variables:
    integer(IntKi)                :: UnIn,UnIn2                                ! Unit number for reading file
    character(1024)               :: FileName                              ! name of the files containing obesever location
    integer(IntKi)                :: ErrStat2                                ! Temporary Error status
    logical                       :: Echo                                      ! Determines if an echo file should be written
    character(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
    character(1024)               :: PriPath                                   ! Path name of the primary file
    character(1024)               :: FTitle                                    ! "File Title": the 2nd line of the input file, which contains a description of its contents
    character(200)                :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
    character(*), parameter       :: RoutineName = 'readbltable'
    integer(IntKi)                :: nRe, nAoA, nAirfoils  !  Number of Reynolds number, angle of attack, and number of airfoils listed
    integer(IntKi)                :: iAF , iRe, iAoA, iDummy, iBuffer ! loop counters
    real(DbKi),dimension(:,:),ALLOCATABLE :: Buffer
    integer                 :: iLine
    ! Initialize some variables:
    ErrStat = ErrID_None
    ErrMsg  = ""

    CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
    nAirfoils = size(BL_Files)
    do iAF=1,nAirfoils

        FileName = trim(BL_Files(iAF))

        print*,'AeroAcoustics_IO: reading BL table:'//trim(Filename)

        CALL GetNewUnit(UnIn, ErrStat2, ErrMsg2); if(Failed()) return
        CALL OpenFInpFile(UnIn, FileName, ErrStat2, ErrMsg2); if(Failed()) return

        CALL ReadCom(UnIn, FileName, "! Boundary layer", ErrStat2, ErrMsg2); if(Failed()) return
        CALL ReadCom(UnIn, FileName, "! Legend: aoa", ErrStat2, ErrMsg2); if(Failed()) return

        CALL ReadVar(UnIn, FileName, nRe,  "ReListBL",   "", ErrStat2, ErrMsg2); if(Failed()) return
        CALL ReadVar(UnIn, FileName, nAoA, "aoaListBL",  "", ErrStat2, ErrMsg2); if(Failed()) return

        if (iAF==1) then
            CALL AllocAry(InputFileData%Pres_DispThick ,nAoA,nRe,nAirfoils,'InputFileData%Pres_DispThick' ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Suct_DispThick ,nAoA,nRe,nAirfoils,'InputFileData%Suct_DispThick' ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Pres_BLThick   ,nAoA,nRe,nAirfoils,'InputFileData%Pres_BLThick'   ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Suct_BLThick   ,nAoA,nRe,nAirfoils,'InputFileData%Suct_BLThick'   ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Pres_Cf        ,nAoA,nRe,nAirfoils,'InputFileData%Pres_Cf'        ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Suct_Cf        ,nAoA,nRe,nAirfoils,'InputFileData%Suct_Cf'        ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Pres_EdgeVelRat,nAoA,nRe,nAirfoils,'InputFileData%Pres_EdgeVelRat',ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Suct_EdgeVelRat,nAoA,nRe,nAirfoils,'InputFileData%Suct_EdgeVelRat',ErrStat2,ErrMsg2); if (Failed())return

            CALL AllocAry(InputFileData%ReListBL,nRe,'InputFileData%ReListBL',ErrStat2,ErrMsg2); if (Failed())return


            CALL AllocAry(Buffer,nAoA,9, 'Buffer', ErrStat2, ErrMsg2); if(Failed()) return
         endif
        iLine=8
        do iRe=1,nRe
            CALL ReadVar(UnIn, FileName, InputFileData%ReListBL(iRe), 'InputFileData%ReListBL','ReListBL', ErrStat2, ErrMsg2); if(Failed()) return
            InputFileData%ReListBL(iRe) = InputFileData%ReListBL(iRe)  * 1.e+006
            CALL ReadCom(UnIn, FileName, "aoa     Ue_Vinf_SS     Ue_Vinf_PS      Dstar_SS     Dstar_PS   Theta_SS   Theta_PS    Cf_SS   Cf_PS", ErrStat2, ErrMsg2); if(Failed()) return
            CALL ReadCom(UnIn, FileName, "(deg)   (-)            (-)             (-)          (-)        (-)        (-)         (-)     (-)",   ErrStat2, ErrMsg2); if(Failed()) return

            call ReadRealMatrix(UnIn, FileName, Buffer, 'BL Matrix', nAoA, 9, ErrStat2, ErrMsg2, iLine)

            if(Failed()) return
            do iAoA=1,nAoA
                InputFileData%Suct_EdgeVelRat(iAoA,iRe,iAF)= Buffer(iAoA, 2) ! EdgeVelRat1 Suction
                InputFileData%Pres_EdgeVelRat(iAoA,iRe,iAF)= Buffer(iAoA, 3) ! EdgeVelRat2 Pressure
                InputFileData%Suct_DispThick (iAoA,iRe,iAF)= Buffer(iAoA, 4) ! dStarAll1 Suction
                InputFileData%Pres_DispThick (iAoA,iRe,iAF)= Buffer(iAoA, 5) ! dStarAll2 Pressure
                InputFileData%Suct_BLThick   (iAoA,iRe,iAF)= Buffer(iAoA, 6) ! d99All1 Suction
                InputFileData%Pres_BLThick   (iAoA,iRe,iAF)= Buffer(iAoA, 7) ! d99All2 Pressure
                InputFileData%Suct_Cf        (iAoA,iRe,iAF)= Buffer(iAoA, 8) ! CfAll1 Suction
                InputFileData%Pres_Cf        (iAoA,iRe,iAF)= Buffer(iAoA, 9) ! CfAll2 Pressure
            enddo
        enddo

        if (iAF == 1) then
            CALL AllocAry(InputFileData%AoAListBL,nAoA, 'InputFileData%AoAListBL', ErrStat2, ErrMsg2); if(Failed()) return
                do iAoA=1,nAoA
                    InputFileData%AoAListBL(iAoA)= Buffer(iAoA, 1) ! AoA
                enddo
        endif
        
        if (InputFileData%IBLUNT==1) then
            call ReadCom(UnIn, FileName, 'Comment' , ErrStat2, ErrMsg2)
            call ReadCom(UnIn, FileName, 'Comment' , ErrStat2, ErrMsg2)
            call ReadVar(UnIn, FileName, InputFileData%BladeProps(iAF)%TEAngle, 'TEAngle', 'TE Angle',ErrStat2, ErrMsg2); if(Failed()) return
            call ReadVar(UnIn, FileName, InputFileData%BladeProps(iAF)%TEThick, 'TEThick', 'TE Thick',ErrStat2, ErrMsg2); if(Failed()) return
        else
            InputFileData%BladeProps(iAF)%TEAngle = 0._ReKi
            InputFileData%BladeProps(iAF)%TEThick = 0._ReKi
        endif
        
        if (UnIn > 0) CLOSE(UnIn)

    enddo
    CALL Cleanup( )
CONTAINS
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if(Failed) call cleanup()
    end function Failed
    SUBROUTINE Cleanup()
        IF (UnIn > 0) CLOSE ( UnIn )
    END SUBROUTINE Cleanup
END SUBROUTINE ReadBLTables
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadTICalcTables(InputFile, InputFileData, ErrStat, ErrMsg)
    ! Passed variables
    integer(IntKi),     intent(out)     :: ErrStat                             ! Error status
    character(*),       intent(out)     :: ErrMsg                              ! Error message
    type(AA_InputFile), intent(inout)   :: InputFileData                       ! All the data in the Noise input file
    character(*),       intent(in)      :: InputFile                           ! Name of the file containing the primary input data
    ! Local variables:
    integer(IntKi)                :: I                                         ! loop counter
    integer(IntKi)                :: UnIn,UnIn2                                ! Unit number for reading file
    integer(IntKi)                :: loop1                                     ! loop counter
    character(1024)               :: FileName                              ! name of the files containing obesever location
    integer(IntKi)                :: ErrStat2, IOS,cou                             ! Temporary Error status
    logical                       :: Echo                                      ! Determines if an echo file should be written
    character(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
    character(1024)               :: PriPath                                   ! Path name of the primary file
    character(1024)               :: FTitle                                    ! "File Title": the 2nd line of the input file, which contains a description of its contents
    character(200)                :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
    character(*), parameter       :: RoutineName = 'REadTICalcTables'
    integer(IntKi)                :: GridY                                     !
    integer(IntKi)                :: GridZ                                    !
    integer(IntKi)                :: cou1
    ! Initialize some variables:
    ErrStat = ErrID_None
    ErrMsg  = ""

    CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.

    FileName = TRIM(PriPath)//InputFileData%TICalcTabFile

    CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2); call check()
    CALL OpenFInpFile ( UnIn, FileName, ErrStat2, ErrMsg2 ); if(Failed()) return
    CALL ReadCom(UnIn, FileName, 'Text Line', ErrStat2, ErrMsg2); call check
    CALL ReadVar(UnIn, FileName, InputFileData%AvgV, 'AvgV',   'Echo flag', ErrStat2, ErrMsg2); call check
    CALL ReadCom(UnIn, FileName, 'Text Line', ErrStat2, ErrMsg2); call check
    CALL ReadVar(UnIn, FileName, GridY, 'GridY',   'Echo flag', ErrStat2, ErrMsg2); call check
    CALL ReadCom(UnIn, FileName, 'Text Line', ErrStat2, ErrMsg2);call check
    CALL ReadVar(UnIn, FileName, GridZ, 'GridZ',   'Echo flag', ErrStat2, ErrMsg2); call check
    CALL ReadCom(UnIn, FileName, 'Text Line', ErrStat2, ErrMsg2); call check
    CALL ReadVar(UnIn, FileName, InputFileData%dy_turb_in, 'InputFileData%dy_turb_in',   'Echo flag', ErrStat2, ErrMsg2); call check
    CALL ReadCom(UnIn, FileName, 'Text Line', ErrStat2, ErrMsg2); call check
    CALL ReadVar(UnIn, FileName, InputFileData%dz_turb_in, 'InputFileData%dz_turb_in',   'Echo flag', ErrStat2, ErrMsg2); call check
    if(Failed()) return

    CALL AllocAry( InputFileData%TI_Grid_In,GridZ,GridY,'InputFileData%TI_Grid_In', ErrStat2, ErrMsg2);
    if(Failed()) return
    DO cou1=1,size(InputFileData%TI_Grid_In,1)
        read(UnIn,*)  InputFileData%TI_Grid_In(cou1,:)
    ENDDO
    !---------------------- END OF FILE -----------------------------------------
    CALL Cleanup( )

CONTAINS
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if(Failed) call cleanup()
    end function Failed
   SUBROUTINE Check()
       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   END SUBROUTINE Check
   SUBROUTINE Cleanup()
       IF (UnIn > 0) CLOSE ( UnIn )
   END SUBROUTINE Cleanup
END SUBROUTINE REadTICalcTables
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine validates the inputs from the AeroDyn input files.
SUBROUTINE ValidateInputData( InputFileData, NumBl, ErrStat, ErrMsg )
   type(AA_InputFile),       intent(in)     :: InputFileData                       !< All the data in the AeroDyn input file
   integer(IntKi),           intent(in)     :: NumBl                               !< Number of blades
   integer(IntKi),           intent(out)    :: ErrStat                             !< Error status
   character(*),             intent(out)    :: ErrMsg                              !< Error message
   ! local variables
   integer(IntKi)                           :: k                                   ! Blade number
   integer(IntKi)                           :: j                                   ! node number
   character(*), parameter                  :: RoutineName = 'ValidateInputData'
   ErrStat = ErrID_None
   ErrMsg  = ""
   if (NumBl > MaxBl .or. NumBl < 1) call SetErrStat( ErrID_Fatal, 'Number of blades must be between 1 and '//trim(num2lstr(MaxBl))//'.', ErrSTat, ErrMsg, RoutineName )
   if (InputFileData%DT_AA <= 0.0)  call SetErrStat ( ErrID_Fatal, 'DT_AA must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%IBLUNT /= IBLUNT_None .and. InputFileData%IBLUNT /= IBLUNT_BPM) then
       call SetErrStat ( ErrID_Fatal, &
           'IBLUNT must '//trim(num2lstr(IBLUNT_None))//' (none) or '//trim(num2lstr(IBLUNT_BPM))//' (Bluntness noise calculated).', ErrStat, ErrMsg, RoutineName )
   endif
   if (InputFileData%ILAM /= ILAM_None .and. InputFileData%ilam /= ILAM_BPM) then
       call SetErrStat ( ErrID_Fatal, 'ILAM must be '//trim(num2lstr(ILAM_None))//' No calculation '//&
           trim(num2lstr(ILAM_BPM))//' (ILAM Calculated).', ErrStat, ErrMsg, RoutineName )
   end if
   if (InputFileData%ITIP /= ITIP_None .and. InputFileData%ITIP /= ITIP_ON) then
       call SetErrStat ( ErrID_Fatal, 'ITIP must be '//trim(num2lstr(ITIP_None))//' (Off) or '//&
           trim(num2lstr(ITIP_On))//' (ITIP On).', ErrStat, ErrMsg, RoutineName )
   end if
   if (InputFileData%ITRIP /= ITRIP_None .and. InputFileData%ITRIP /= ITRIP_Heavy .and. InputFileData%ITRIP /= ITRIP_Light) then
       call SetErrStat ( ErrID_Fatal,'ITRIP must be '//trim(num2lstr(ITRIP_None))//' (none) or '//trim(num2lstr(ITRIP_Heavy))//&
           ' (heavily tripped BL Calculation) or '//trim(num2lstr(ITRIP_Light))//' (lightly tripped BL)' ,ErrStat, ErrMsg, RoutineName )
   end if
   if (InputFileData%ITURB /= ITURB_None .and. InputFileData%ITURB /= ITURB_BPM .and. InputFileData%ITURB /= ITURB_TNO) then
       call SetErrStat ( ErrID_Fatal, 'ITURB must be 0 (off) or 1 (BPM) or 2 (TNO) .', ErrStat, ErrMsg, RoutineName )
   end if
   if (InputFileData%IInflow /= IInflow_None .and. InputFileData%IInflow /= IInflow_BPM &
       .and. InputFileData%IInflow /= IInflow_FullGuidati .and. InputFileData%IInflow /= IInflow_SimpleGuidati ) then
       call SetErrStat ( ErrID_Fatal, 'IInflow must be 0 (off) or 1 (only Amiet)  or 2 (Full Guidati)'//&
           'or 3 (Simple Guidati).', ErrStat, ErrMsg, RoutineName )
   end if
   if (InputFileData%TICalcMeth /= TICalc_Every .and. InputFileData%TICalcMeth /= TICalc_Interp ) then
       call SetErrStat ( ErrID_Fatal, 'TICalcMeth must be '//trim(num2lstr(TICalc_Every))//' TICalc automatic  or '//&
           trim(num2lstr(TICalc_Interp))//' (TICalcMeth interp).', ErrStat, ErrMsg, RoutineName )
   end if

   if (InputFileData%X_BLMethod /= X_BLMethod_BPM .and. InputFileData%X_BLMethod /= X_BLMethod_Tables) then
       call SetErrStat ( ErrID_Fatal, 'X_BLMethod must be '//trim(num2lstr(X_BLMethod_BPM))//' X_BLMethod_ with BPM or '//&
           trim(num2lstr(X_BLMethod_Tables))//' (X_BLMethod with BL tables).', ErrStat, ErrMsg, RoutineName )
   end if
   if (InputFileData%NrObsLoc <= 0.0) call SetErrStat ( ErrID_Fatal, 'Number of Observer Locations should be greater than zero', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%NrOutFile /= 1 .and. InputFileData%NrOutFile /= 2 .and. InputFileData%NrOutFile /= 3 &
       .and. InputFileData%NrOutFile /= 4) then
       call SetErrStat ( ErrID_Fatal, ' NrOutFile must be 1 or 2 or 3 or 4', ErrStat, ErrMsg, RoutineName )
   end if
END SUBROUTINE ValidateInputData

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AA_PrintSum( InputFileData, p, u, y, ErrStat, ErrMsg )
    ! This routine generates the summary file, which contains a summary of input file options.
    ! passed variables
    TYPE(AA_InputFile),        INTENT(IN)  :: InputFileData                        ! Input-file data
    TYPE(AA_ParameterType),    INTENT(IN)  :: p                                    ! Parameters
    TYPE(AA_InputType),        INTENT(IN)  :: u                                    ! inputs
    TYPE(AA_OutputType),       INTENT(IN)  :: y                                    ! outputs
    INTEGER(IntKi),            INTENT(OUT) :: ErrStat
    CHARACTER(*),              INTENT(OUT) :: ErrMsg
    ! Local variables.
    INTEGER(IntKi)               :: I                                               ! Index for the nodes.
    INTEGER(IntKi)               :: UnSu                                            ! I/O unit number for the summary output file
    CHARACTER(*), PARAMETER      :: FmtDat    = '(A,T35,1(:,F13.3))'                ! Format for outputting mass and modal data.
    CHARACTER(*), PARAMETER      :: FmtDatT   = '(A,T35,1(:,F13.8))'                ! Format for outputting time steps.
    CHARACTER(30)                :: OutPFmt                                         ! Format to print list of selected output channels to summary file
    CHARACTER(100)               :: Msg                                             ! temporary string for writing appropriate text to summary file
    ! Open the summary file and give it a heading.
    ErrStat = ErrID_None
    ErrMsg  = ""
    RETURN
END SUBROUTINE AA_PrintSum
!..................................................................................................................................
!> This subroutine sets the initialization output data structure, which contains data to be returned to the calling program (e.g.,
!! FAST or AeroAcoustics_Driver)
subroutine AA_SetInitOut(p, InputFileData, InitOut, errStat, errMsg)
    type(AA_InitOutputType),       intent(  out)  :: InitOut          ! output data
    type(AA_InputFile),            intent(in   )  :: InputFileData    ! input file data (for setting airfoil shape outputs)
    type(AA_ParameterType),        intent(in   )  :: p                ! Parameters
    integer(IntKi),                intent(  out)  :: errStat          ! Error status of the operation
    character(*),                  intent(  out)  :: errMsg           ! Error message if ErrStat /= ErrID_None
    ! Local variables
    integer(intKi)                               :: ErrStat2          ! temporary Error status
    character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
    character(*), parameter                      :: RoutineName = 'AA_SetInitOut'
    integer(IntKi)                               :: i, j, k,m,oi
    integer(IntKi)                               :: NumCoords
    character(500)                               :: chanPrefix
    ! Initialize variables for this routine
    errStat = ErrID_None
    errMsg  = ""
    InitOut%AirDens = p%AirDens
    ! FIRST  FILE HEADER,UNIT
    call AllocAry(InitOut%WriteOutputHdr, p%numOuts, 'WriteOutputHdr', errStat2, errMsg2); if(Failed()) return
    call AllocAry(InitOut%WriteOutputUnt, p%numOuts, 'WriteOutputUnt', errStat2, errMsg2); if(Failed()) return
    do j=1,p%NrObsLoc
        InitOut%WriteOutputHdr(j)="Obs"//trim(num2lstr(j))
        InitOut%WriteOutputUnt(j) = "OASPL"
    enddo

    ! SECOND FILE HEADER,UNIT
    call AllocAry(InitOut%WriteOutputHdrforPE, p%numOutsforPE, 'WriteOutputHdrforPE', errStat2, errMsg2); if(Failed()) return
    call AllocAry(InitOut%WriteOutputUntforPE, p%numOutsforPE, 'WriteOutputUntforPE', errStat2, errMsg2); if(Failed()) return
    i=0
    do j=1,p%NrObsLoc
        do k=1,size(p%FreqList)
         i=i+1
         InitOut%WriteOutputHdrforPE(i) = "Obs"//trim(num2lstr(j))//"_Freq"//trim(num2lstr(p%FreqList(k)))
         if(p%aweightflag .eqv. .TRUE.) then
            InitOut%WriteOutputUntforPE(i) = "SPL_A"
         else
            InitOut%WriteOutputUntforPE(i) = "SPL"
         endif
        end do
    enddo
    ! THIRD FILE HEADER,UNIT
    call AllocAry(InitOut%WriteOutputHdrSep, p%NumOutsForSep, 'WriteOutputHdrSep', errStat2, errMsg2); if(Failed()) return
    call AllocAry(InitOut%WriteOutputUntSep, p%NumOutsForSep, 'WriteOutputUntSep', errStat2, errMsg2); if(Failed()) return
    i=0
    do j=1,p%NrObsLoc
      do k=1,size(p%FreqList)
         do oi=1,7
            i=i+1
            InitOut%WriteOutputHdrSep(i) = "Obs"//trim(num2lstr(j))//"_Freq"//trim(num2lstr(p%FreqList(k)))//"_Type"//trim(num2lstr(oi))
            InitOut%WriteOutputHdrSep(i)=trim(InitOut%WriteOutputHdrSep(i))
            if(p%aweightflag .eqv. .TRUE.) then
               InitOut%WriteOutputUntSep(i) = "SPL_A"
            else
               InitOut%WriteOutputUntSep(i) = "SPL"
            endif
         enddo
      enddo
    enddo

    ! FOURTH FILE HEADER,UNIT
    call AllocAry(InitOut%WriteOutputHdrNodes,p%numBlades*p%NumBlNds*p%NrObsLoc, 'InitOut%WriteOutputHdrNodes', errStat2, errMsg2); if(Failed()) return
    call AllocAry(InitOut%WriteOutputUntNodes,p%numBlades*p%NumBlNds*p%NrObsLoc, 'InitOut%WriteOutputUntNodes', errStat2, errMsg2); if(Failed()) return
    i=0
    do oi = 1,p%numBlades
        do k = 1,p%NumBlNds
            do j = 1,p%NrObsLoc
                i=i+1
                InitOut%WriteOutputHdrNodes(i) = "Bld"//trim(num2lstr(oi))//"Node"//trim(num2lstr(k))//"Obs"//trim(num2lstr(j))
                InitOut%WriteOutputUntNodes(i) = "SPL"
            enddo
        enddo
    enddo
    InitOut%Ver = AA_Ver
    InitOut%delim = Tab 

contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
    end function Failed
end subroutine AA_SetInitOut
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AA_InitializeOutputFile(p, InputFileData,InitOut,errStat, errMsg)
   type(AA_InputFile),       intent(in   ) :: InputFileData    !< All the data in the AeroDyn input file
   type(AA_ParameterType) ,  intent(inout) :: p                !<
   type(AA_InitOutputType),  intent(in  )  :: InitOut          !< output data
   integer(IntKi)         ,  intent(inout) :: errStat          !< Status of error message
   character(*)           ,  intent(inout) :: errMsg           !< Error message if ErrStat /= ErrID_None
   ! locals
   integer(IntKi) :: i
   integer(IntKi) :: numOuts
   character(200) :: frmt                                      ! A string to hold a format specifier
   character(15)  :: tmpStr                                    ! temporary string to print the time output as text
   ! FIRST FILE
   IF (InputFileData%NrOutFile .gt.0) THEN
      call GetNewUnit( p%unOutFile, ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) then
         p%unOutFile = -1
         return
      end if

      call OpenFOutFile ( p%unOutFile, trim(InputFileData%AAOutFile(1)), ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) return

      write (p%unOutFile,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using AA '//trim(GetNVD(InitOut%ver))
      write (p%unOutFile,'(A)')  ''
      write( p%unOutFile,'(A,I5)' )      'Number of observers      :', p%NrObsLoc
      write (p%unOutFile,'(A)')  'Description from AA input file, line2: '//trim(InputFileData%FTitle)
      write (p%unOutFile,'(A)')  ''
      numOuts = size(InitOut%WriteOutputHdr)
      !......................................................
      ! Write the names of the output parameters on one line: line 7
      !......................................................
      call WrFileNR ( p%unOutFile, '     Time           ' )
      do i=1,NumOuts
         call WrFileNR ( p%unOutFile, InitOut%delim//InitOut%WriteOutputHdr(i) )
      end do ! i
      write (p%unOutFile,'()')
      !......................................................
      ! Write the units of the output parameters on one line: line 8
      !......................................................
      call WrFileNR ( p%unOutFile, '      (s)           ' )
      do i=1,NumOuts
         call WrFileNR ( p%unOutFile, InitOut%delim//InitOut%WriteOutputUnt(i) )
      end do ! i
      write (p%unOutFile,'()')
   ENDIF
   ! SECOND FILE
   IF (InputFileData%NrOutFile .gt. 1) THEN
      call GetNewUnit( p%unOutFile2, ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) then
         p%unOutFile = -1
         return
      end if
      call OpenFOutFile ( p%unOutFile2, trim(InputFileData%AAOutFile(2)), ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) return
      write (p%unOutFile2,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using AA '//trim(GetNVD(InitOut%Ver))
      write (p%unOutFile2,'(A)')  ''
      write( p%unOutFile2,'(A,I5,A,I5)' )      'Number of observers      :', p%NrObsLoc,';  Number of frequencies    :', size(p%FreqList)
      write (p%unOutFile2,'(A)')  'Description from AA input file, line2: '//trim(InputFileData%FTitle)
      write (p%unOutFile2,'(A)')  ''
      numOuts = size(InitOut%WriteOutputHdrforPE)
      !......................................................
      ! Write the names of the output parameters on one line: line 7
      !......................................................
      call WrFileNR ( p%unOutFile2, '     Time           ' )
      do i=1,NumOuts
         call WrFileNR ( p%unOutFile2, InitOut%delim//InitOut%WriteOutputHdrforPE(i) )
      end do ! i
      write (p%unOutFile2,'()')
      !......................................................
      ! Write the units of the output parameters on one line: line 8
      !......................................................
      call WrFileNR ( p%unOutFile2, '      (s)           ' )
      do i=1,NumOuts
         call WrFileNR ( p%unOutFile2, InitOut%delim//InitOut%WriteOutputUntforPE(i) )
      end do ! i
      write (p%unOutFile2,'()')
      !frmt = '"'//p%delim//'"'//trim(p%outFmt)      ! format for array elements from individual modules
      !call WrNumAryFileNR ( p%unOutFile2, p%FreqList,  frmt, errStat, errMsg )
      !if ( errStat >= AbortErrLev ) return
      !write (p%unOutFile2,'()')
   ENDIF
   ! THIRD FILE
   IF (InputFileData%NrOutFile .gt. 2) THEN
      call GetNewUnit( p%unOutFile3, ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) then
         p%unOutFile = -1
         return
      end if
      call OpenFOutFile ( p%unOutFile3, trim(InputFileData%AAOutFile(3)), ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) return
      write (p%unOutFile3,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using AA '//trim(GetNVD(InitOut%Ver))
      write (p%unOutFile3,'(A)')  ''
      write( p%unOutFile3,'(A,I5,A,I5)' )      'Number of observers      :', p%NrObsLoc,';  Number of frequencies    :', size(p%FreqList)
      write (p%unOutFile3,'(A)')  'Description from AA input file, line2: '//trim(InputFileData%FTitle)
      numOuts = size(InitOut%WriteOutputHdrSep)
      !......................................................
      ! Write the names of the output parameters on one line: line 7
      !......................................................
      call WrFileNR ( p%unOutFile3,  "1-LBL 2-TBLPres 3-TBLSuc 4-Sep  5-BLUNT 6-TIP 7-Inflow")
      write (p%unOutFile3,'()')
      call WrFileNR ( p%unOutFile3, '     Time           ' )
      do i=1,NumOuts
         call WrFileNR ( p%unOutFile3, InitOut%delim//InitOut%WriteOutputHdrSep(i) )
      end do ! i
      write (p%unOutFile3,'()')
      !......................................................
      ! Write the units of the output parameters on one line: line 8
      !......................................................
      call WrFileNR ( p%unOutFile3, '      (s)           ' )

      do i=1,NumOuts
         call WrFileNR ( p%unOutFile3, InitOut%delim//InitOut%WriteOutputUntSep(i) )
      end do ! i
      write (p%unOutFile3,'()')
   ENDIF
   ! FOURTH FILE
   IF (InputFileData%NrOutFile .gt. 3) THEN
      call GetNewUnit( p%unOutFile4, ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) then
         p%unOutFile = -1
         return
      end if
      call OpenFOutFile ( p%unOutFile4, trim(InputFileData%AAOutFile(4)), ErrStat, ErrMsg )
      if ( ErrStat >= AbortErrLev ) return
      write (p%unOutFile4,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using AA '//trim(GetNVD(InitOut%Ver))
      write (p%unOutFile4,'()')
      write( p%unOutFile4,'(A,I5)' )      'Number of observers      :', p%NrObsLoc, ';        Number of blades         :', p%numBlades,'     Number of nodes per blade:', p%NumBlNds
      write (p%unOutFile4,'(A)')  'Description from AA input file, line2: '//trim(InputFileData%FTitle)
      numOuts = size(InitOut%WriteOutputHdrNodes)
      !......................................................
      ! Write the names of the output parameters on one line: line 7
      !......................................................
      write (p%unOutFile4,'()')
      call WrFileNR ( p%unOutFile4, '     Time           ' )
      do i=1,NumOuts
         call WrFileNR ( p%unOutFile4, InitOut%delim//InitOut%WriteOutputHdrNodes(i) )
      end do ! i
      write (p%unOutFile4,'()')
      !......................................................
      ! Write the units of the output parameters on one line: line 8
      !......................................................
      call WrFileNR ( p%unOutFile4, '      (s)           ' )
      do i=1,NumOuts
         call WrFileNR ( p%unOutFile4, InitOut%delim//InitOut%WriteOutputUntNodes(i) )
      end do ! i
      write (p%unOutFile4,'()')
   ENDIF
end subroutine AA_InitializeOutputFile
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AA_WriteOutputLine(y, t, p, errStat, errMsg)
    real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
    type(AA_OutputType)    ,  intent(in   )   :: y
    type(AA_ParameterType) ,  intent(in   )   :: p
    integer(IntKi)         ,  intent(inout)   :: errStat              ! Status of error message
    character(*)           ,  intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
    ! Local variables.
    character(200)                   :: frmt                                      ! A string to hold a format specifier
    character(15)                    :: tmpStr                                    ! temporary string to print the time output as text
    integer :: numOuts
    errStat = ErrID_None
    errMsg  = ''
    ! FIRST FILE
    IF (p%NrOutFile .gt. 0) THEN
       numOuts = size(y%WriteOutput)
       frmt = '"'//p%delim//'"'//trim(p%outFmt)      ! format for array elements from individual modules
       ! time
       write( tmpStr, '(F15.4)' ) t
       call WrFileNR( p%unOutFile, tmpStr )
       call WrNumAryFileNR ( p%unOutFile, y%WriteOutput,  frmt, errStat, errMsg )
       if ( errStat >= AbortErrLev ) return
       ! write a new line (advance to the next line)
       write (p%unOutFile,'()')
    ENDIF

    !! SECOND FILE
    IF (p%NrOutFile .gt. 1) THEN
       numOuts = size(y%WriteOutputforPE)
       frmt = '"'//p%delim//'"'//trim(p%outFmt)      ! format for array elements from individual modules
       ! time
       write( tmpStr, '(F15.4)' ) t
       call WrFileNR( p%unOutFile2, tmpStr )
       call WrNumAryFileNR ( p%unOutFile2, y%WriteOutputforPE,  frmt, errStat, errMsg )
       if ( errStat >= AbortErrLev ) return
       ! write a new line (advance to the next line)
       write (p%unOutFile2,'()')
    ENDIF
    ! THIRD FILE
    IF (p%NrOutFile .gt. 2) THEN
       numOuts = size(y%WriteOutputSep)
       frmt = '"'//p%delim//'"'//trim(p%outFmt)      ! format for array elements from individual modules
       ! time
       write( tmpStr, '(F15.4)' ) t
       call WrFileNR( p%unOutFile3, tmpStr )
       call WrNumAryFileNR ( p%unOutFile3, y%WriteOutputSep,  frmt, errStat, errMsg )
       if ( errStat >= AbortErrLev ) return
       ! write a new line (advance to the next line)
       write (p%unOutFile3,'()')
    ENDIF
    ! Fourth FILE
    IF (p%NrOutFile .gt. 3) THEN
       numOuts = size(y%WriteOutputNode)
       frmt = '"'//p%delim//'"'//trim(p%outFmt)      ! format for array elements from individual modules
       ! time
       write( tmpStr, '(F15.4)' ) t
       call WrFileNR( p%unOutFile4, tmpStr )
       call WrNumAryFileNR ( p%unOutFile4, y%WriteOutputNode,  frmt, errStat, errMsg )
       if ( errStat >= AbortErrLev ) return
       ! write a new line (advance to the next line)
       write (p%unOutFile4,'()')
    ENDIF
end subroutine AA_WriteOutputLine
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Calc_WriteOutput( p, u, m, y, ErrStat, ErrMsg )
   TYPE(AA_ParameterType),    INTENT(IN   )  :: p                                 ! The module parameters
   TYPE(AA_InputType),        INTENT(IN   )  :: u                                 ! inputs
   TYPE(AA_MiscVarType),      INTENT(INOUT)  :: m                                 ! misc variables
   TYPE(AA_OutputType),       INTENT(INOUT)  :: y                                 ! outputs
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat                           ! The error status code
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg                            ! The error message, if an error occurred
   ! local variables
   CHARACTER(*), PARAMETER                   :: RoutineName = 'Calc_WriteOutput'
   INTEGER(intKi)                            :: ErrStat2
   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   INTEGER(IntKi)                            :: j,k,counter,i,oi,III
   ! start routine:
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! FOR THE FIRST OUTPUT FILE
   IF (p%NrOutFile .gt. 0) THEN
      y%WriteOutput(1:p%NrObsLoc)=y%DirectiviOutput
   endif

   ! FOR THE SECOND OUTPUT FILE
   IF (p%NrOutFile .gt. 1) THEN
      counter=0
      DO K = 1,p%NrObsLoc
         DO III = 1,size(p%FreqList)
            counter=counter+1
            y%WriteOutputforPE(counter) = y%PtotalFreq(K,III)
         END DO !
      END DO !
   ENDIF

   ! FOR THE THIRD OUTPUT FILE
   IF (p%NrOutFile .gt. 2) THEN
      counter=0
      do K = 1,p%NrObsLoc
         do III = 1,size(p%FreqList)
            do oi=1,size(y%OASPL_Mech,1)
               counter=counter+1
               y%WriteOutputSep(counter) = y%SumSpecNoiseSep(oi,K,III)
            enddo
         enddo
      enddo
   ENDIF

    ! FOR THE FOURTH OUTPUT FILE
    IF (p%NrOutFile .gt. 3) THEN
        counter=0
        DO I = 1,p%numBlades
            DO J = 1,p%NumBlNds
                DO K = 1,p%NrObsLoc
                    counter=counter+1
                    y%WriteOutputNode(counter) = y%OASPL(K,J,I)
                END DO !
            END DO !
        ENDDO
    ENDIF
END SUBROUTINE Calc_WriteOutput
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE AeroAcoustics_IO
