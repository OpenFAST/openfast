!**********************************************************************************************************************************
! File last committed: 2020-02-12
!**********************************************************************************************************************************
MODULE AeroAcoustics_IO

   use NWTC_Library
   use AeroAcoustics_Types

   implicit none

   type(ProgDesc), parameter  :: AA_Ver = ProgDesc( 'AeroAcoustics', '', '' )
   character(*),   parameter  :: AA_Nickname = 'AA'
   character(*),   parameter  :: delim = Tab

   LOGICAL,        parameter  :: AA_OutputToSeparateFile = .true.

   integer(intKi), parameter  :: nNoiseMechanism = 7  ! number of noise mechanisms

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

   integer(intKi), parameter        :: TICalc_Interp  = 1  ! interpolate from pretabulated (TICalcMethod)
   integer(intKi), parameter        :: TICalc_Every   = 2  ! calculate ti automatically (TICalcMethod)

   integer(intKi), parameter        :: ITURB_None           = 0  ! TBLTE noise is not calculated
   integer(intKi), parameter        :: ITURB_BPM            = 1  ! TBLTE noise is calculated with BPM
   integer(intKi), parameter        :: ITURB_TNO            = 2  ! TBLTE noise is calculated with TNO

   integer(intKi), parameter        :: IInflow_None             = 0  ! IInflow noise is not calculated
   integer(intKi), parameter        :: IInflow_BPM              = 1  ! IInflow noise is calculated with BPM
   integer(intKi), parameter        :: IInflow_FullGuidati      = 2  ! IInflow noise is calculated with FullGuidati
   integer(intKi), parameter        :: IInflow_SimpleGuidati    = 3  ! IInflow noise is calculated with SimpleGuidati

contains
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadInputFiles( InputFileName, AFInfo, InputFileData, Default_DT, OutFileRoot, ErrStat, ErrMsg )
    ! This subroutine reads the input file and stores all the data in the AA_InputFile structure.
    ! It does not perform data validation.
    !..................................................................................................................................
    ! Passed variables
    REAL(DbKi),              INTENT(IN)    :: Default_DT      ! The default DT (from glue code)
    CHARACTER(*),            INTENT(IN)    :: InputFileName   ! Name of the aeroacoustics input file
    TYPE(AFI_ParameterType), INTENT(IN)    :: AFInfo(:)          ! airfoil array: contains names of the BL input file
    CHARACTER(*),            INTENT(IN)    :: OutFileRoot     ! The rootname of all the output files written by this routine.
    TYPE(AA_InputFile),      INTENT(OUT)   :: InputFileData   ! Data stored in the module's input file
    INTEGER(IntKi),          INTENT(OUT)   :: ErrStat         ! The error status code
    CHARACTER(*),            INTENT(OUT)   :: ErrMsg          ! The error message, if an error occurred
    ! local variables
    INTEGER(IntKi)                         :: UnEcho          ! Unit number for the echo file
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
    ALLOCATE( InputFileData%BladeProps( size(AFInfo) ), STAT = ErrStat2 )
    IF (ErrStat2 /= 0) THEN
        CALL SetErrStat(ErrID_Fatal,"Error allocating memory for BladeProps.", ErrStat, ErrMsg, RoutineName)
        call cleanup()
        return
    END IF

    if (InputFileData%ITURB==ITURB_TNO .or. InputFileData%X_BLMethod==X_BLMethod_Tables .or. InputFileData%IBLUNT==IBLUNT_BPM) then
        ! We need to read the BL tables
        CALL ReadBLTables( InputFileName, AFInfo, InputFileData, UnEcho, ErrStat2, ErrMsg2 )
        if (Failed()) return
    endif

CONTAINS
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call cleanup()
    end function Failed
    subroutine cleanup()
       if (UnEcho > 0) close(UnEcho)
    end subroutine

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
    integer(IntKi)                :: I                                         ! loop counter
    integer(IntKi)                :: UnIn,UnIn2                                ! Unit number for reading file
    character(1024)               :: ObserverFile                              ! name of the files containing obesever location
    integer(IntKi)                :: ErrStat2, cou                             ! Temporary Error status
    logical                       :: Echo                                      ! Determines if an echo file should be written
    character(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
    character(1024)               :: PriPath                                   ! Path name of the primary file
    character(1024)               :: OutPath                                   ! Path name of the default output file
    character(200)                :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
    character(*), parameter       :: RoutineName = 'ReadPrimaryFile'
    real(ReKi)                    :: TmpArray(3)
    ! Initialize some variables:
    ErrStat = ErrID_None
    ErrMsg  = ""

    UnEc = -1
    UnIn = -1
    UnIn2 = -1
    Echo = .FALSE.
    CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.

    ! Open the Primary input file.
    CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 ); if (Failed()) return;
    CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 ); if (Failed()) return;

    ! Read the lines up/including to the "Echo" simulation control variable
    ! If echo is FALSE, don't write these lines to the echo file.
    ! If Echo is TRUE, rewind and write on the second try.
    I = 1 !set the number of times we've read the file
    DO
        !----------- HEADER -------------------------------------------------------------
        CALL ReadCom( UnIn, InputFile, 'File header: Module Version (line 1)', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return;
        CALL ReadStr( UnIn, InputFile, InputFileData%FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return;

        !----------- GENERAL OPTIONS ----------------------------------------------------
        CALL ReadCom( UnIn, InputFile, 'Section Header: General Options', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return;
        ! Echo - Echo input to "<RootName>.AD.AA.ech".
        CALL ReadVar( UnIn, InputFile, Echo, 'Echo',   'Echo flag', ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
        
        IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop
        ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
        I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
        CALL OpenEcho ( UnEc, TRIM(OutFileRoot)//'.ech', ErrStat2, ErrMsg2, AA_Ver ); if (Failed()) return;

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
    CALL ReadVarWDefault( UnIn, InputFile, InputFileData%DT_AA, "DT_AA", "Time interval for aeroacoustics calculations {or default} (s)", Default_DT, ErrStat2, ErrMsg2, UnEc )
    if (Failed()) return;
    
    IF (.NOT. EqualRealNos( InputFileData%DT_AA, NINT(InputFileData%DT_AA / Default_DT)*Default_DT ) ) THEN
       CALL SetErrStat(ErrID_Fatal,"The Aeroacoustics input DT_AA must be a multiple of DTAero.", ErrStat, ErrMsg, RoutineName)
       call Cleanup()
       return
    END IF
    

    CALL ReadVar(UnIn,InputFile,InputFileData%AAStart      ,"AAStart"      ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%AA_Bl_Prcntge,"BldPrcnt"     ,"-",ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadCom( UnIn, InputFile, 'Section Header: Aeroacoustic Models', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%IInflow      ,"InflowMod"    ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%TICalcMeth   ,"TICalcMeth"   ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVAr(UnIn,InputFile,InputFileData%TI           ,"TI"           ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVAr(UnIn,InputFile,InputFileData%avgV         ,"avgV"         ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%Lturb        ,"Lturb"        ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%ITURB        ,"TurbMod"      ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return; ! ITURB - TBLTE NOISE
    CALL ReadVar(UnIn,InputFile,InputFileData%X_BLMethod   ,"BLMod"        ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%ITRIP        ,"TripMod"      ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%ILAM         ,"LamMod"       ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%ITIP         ,"TipMod"       ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%ROUND        ,"RoundTip"     ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%ALPRAT       ,"ALPRAT"       ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar(UnIn,InputFile,InputFileData%IBLUNT       ,"BluntMod"     ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;

    !----------- OBSERVER INPUT  ------------------------------
    CALL ReadCom( UnIn, InputFile, 'Section Header: Observer Input ', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return;
    !----- read from observer file
    CALL ReadVar ( UnIn, InputFile, ObserverFile, ObserverFile, 'Name of file  observer locations', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return;
    IF ( PathIsRelative( ObserverFile ) ) ObserverFile = TRIM(PriPath)//TRIM(ObserverFile)

    CALL GetNewUnit( UnIn2, ErrStat2, ErrMsg2 ); if (Failed()) return;
    CALL OpenFInpFile ( UnIn2, ObserverFile, ErrStat2, ErrMsg2 ); if (Failed()) return;
    
    ! NrObsLoc  - Nr of Observers (-):
    CALL ReadVar( UnIn2, ObserverFile, InputFileData%NrObsLoc, "NrObsLoc", "Nr of Observers (-)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
       if (InputFileData%NrObsLoc < 1) then
          call SetErrStat(ErrID_Fatal,"NrObsLoc must be a positive number", ErrStat, ErrMsg, RoutineName)
          call Cleanup()
          return
       end if
    
    CALL ReadCom( UnIn2, ObserverFile, ' Header', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return;

    ! Observer location in tower-base coordinate  (m):
    CALL AllocAry( InputFileData%ObsXYZ,3,InputFileData%NrObsLoc, 'ObsX', ErrStat2, ErrMsg2); if (Failed()) return;
    DO cou=1,InputFileData%NrObsLoc
        CALL ReadAry( UnIn2, ObserverFile, InputFileData%ObsXYZ(:,cou), SIZE(TmpArray), 'Observer Locations Line '//trim(Num2LStr(cou)), 'Observer Locations', ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
    ENDDO
    CLOSE ( UnIn2 )
    UnIn2 = -1
    !----- end read from observer file

    !----------- OUTPUTS  -----------------------------------------------------------
    CALL ReadCom( UnIn, InputFile, 'Section Header: Outputs', ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
    CALL ReadVar( UnIn,InputFile,InputFileData%aweightflag  ,"AWeighting"   ,"" ,ErrStat2,ErrMsg2,UnEc); if (Failed()) return;
    CALL ReadVar( UnIn, InputFile, InputFileData%NrOutFile, "NrOutFile", "Nr of Output Files (-)", ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
    if (InputFileData%NrOutFile < 1 .OR. InputFileData%NrOutFile > 4) then
       call SetErrStat(ErrID_Fatal, "NrOutFile must be a value between 1 and 4.", ErrStat, ErrMsg, RoutineName)
       CALL Cleanup( )
       return
    end if
    
    CALL ReadVar ( UnIn, InputFile, InputFileData%AAOutFile(1), 'AAOutFile', 'Name of output file ', ErrStat2, ErrMsg2, UnEc ); if (Failed()) return;
    Line = InputFileData%AAOutFile(1)
    call Conv2UC(Line)
    IF ( INDEX(Line, "DEFAULT" ) /= 1 ) THEN 
       IF ( PathIsRelative( InputFileData%AAOutFile(1) ) ) then
          CALL GetPath( OutFileRoot, OutPath )   ! Output files will be relative to the path where the primary output file is located.
          InputFileData%AAOutFile(1) = TRIM(OutPath)//TRIM(InputFileData%AAOutFile(1))
       END IF
    ELSE ! use default program root
       InputFileData%AAOutFile(1) = TRIM(OutFileRoot)
    ENDIF
    
    DO I=InputFileData%NrOutFile,1,-1
       ! one file name is given by the user and the XXFile1.out XXFile2.out XXFile3.out is generated
       InputFileData%AAOutFile(I) = TRIM(InputFileData%AAOutFile(1))//TRIM(Num2Lstr(I))//".out"
    ENDDO

    !---------------------- END OF FILE -----------------------------------------
    CALL Cleanup( )

CONTAINS
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call cleanup()
    end function Failed
   !...............................................................................................................................
   SUBROUTINE Cleanup()
       IF (UnIn > 0) CLOSE ( UnIn )
       IF (UnIn2 > 0) CLOSE ( UnIn2 )
   END SUBROUTINE Cleanup
   !...............................................................................................................................
END SUBROUTINE ReadPrimaryFile
!----------------------------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadBLTables( InputFile, AFInfo, InputFileData, UnEc, ErrStat, ErrMsg )
    ! Passed variables
    character(*),       intent(in)      :: InputFile                           ! Name of the file containing the primary input data
    TYPE(AFI_ParameterType), INTENT(IN) :: AFInfo(:)                              ! airfoil array: contains names of the BL input file
    type(AA_InputFile), intent(inout)   :: InputFileData                       ! All the data in the Noise input file
    integer(IntKi),     intent(in)      :: UnEc                                ! I/O unit for echo file. If > 0, file is open for writing.
    integer(IntKi),     intent(out)     :: ErrStat                             ! Error status
    character(*),       intent(out)     :: ErrMsg                              ! Error message
    
    ! Local variables:
    integer(IntKi)                :: UnIn                                      ! Unit number for reading file
    character(1024)               :: FileName                                  ! name of the files containing obesever location
    integer(IntKi)                :: ErrStat2                                  ! Temporary Error status
    character(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
    character(1024)               :: PriPath                                   ! Path name of the primary file
    character(*), parameter       :: RoutineName = 'ReadBLTables'
    integer(IntKi)                :: nRe, nAoA, nAirfoils                      !  Number of Reynolds number, angle of attack, and number of airfoils listed
    integer(IntKi)                :: iAF , iRe, iAoA                           ! loop counters
    real(ReKi)                    :: Buffer(9)
    real(ReKi)                    :: TempRe
    
    ! Initialize some variables:
    ErrStat = ErrID_None
    ErrMsg  = ""
    UnIn = -1

    CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
    nAirfoils = size(AFInfo)
    do iAF=1,nAirfoils

        FileName = trim(AFInfo(iAF)%BL_file)

        call WrScr('AeroAcoustics_IO: reading BL table:'//trim(Filename))

        CALL GetNewUnit(UnIn, ErrStat2, ErrMsg2); if(Failed()) return
        CALL OpenFInpFile(UnIn, FileName, ErrStat2, ErrMsg2); if(Failed()) return

        CALL ReadCom(UnIn, FileName, "! Boundary layer", ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
        CALL ReadCom(UnIn, FileName, "! Legend: aoa", ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

        CALL ReadVar(UnIn, FileName, nRe,  "ReListBL",   "", ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
        CALL ReadVar(UnIn, FileName, nAoA, "aoaListBL",  "", ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

        if (iAF==1) then
           if (nAoA < 1 .OR. nRe < 1 ) call SetErrStat(ErrID_Fatal,"ReListBL and aoaListBL must be positive numbers.", ErrStat, ErrMsg, RoutineName)
           
            CALL AllocAry(InputFileData%Pres_DispThick ,nAoA,nRe,nAirfoils,'Pres_DispThick' ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Suct_DispThick ,nAoA,nRe,nAirfoils,'Suct_DispThick' ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Pres_BLThick   ,nAoA,nRe,nAirfoils,'Pres_BLThick'   ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Suct_BLThick   ,nAoA,nRe,nAirfoils,'Suct_BLThick'   ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Pres_Cf        ,nAoA,nRe,nAirfoils,'Pres_Cf'        ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Suct_Cf        ,nAoA,nRe,nAirfoils,'Suct_Cf'        ,ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Pres_EdgeVelRat,nAoA,nRe,nAirfoils,'Pres_EdgeVelRat',ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%Suct_EdgeVelRat,nAoA,nRe,nAirfoils,'Suct_EdgeVelRat',ErrStat2,ErrMsg2); if (Failed())return

            CALL AllocAry(InputFileData%AoAListBL,      nAoA,              'AoAListBL',      ErrStat2,ErrMsg2); if (Failed())return
            CALL AllocAry(InputFileData%ReListBL,       nRe,               'ReListBL',       ErrStat2,ErrMsg2); if (Failed())return
        else
            if (nAoA /= SIZE(InputFileData%Pres_DispThick,1) .OR. &
                 nRe /= SIZE(InputFileData%Pres_DispThick,2) ) then
               call SetErrStat(ErrID_Fatal,'All aeroacoustics airfoils must have the same number of angles of attack and reynolds numbers', ErrStat, ErrMsg, RoutineName)
               call cleanup()
               return
            end if
        endif
        
        do iRe=1,nRe
            CALL ReadVar(UnIn, FileName, TempRe, 'InputFileData%ReListBL','ReListBL', ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
            if (iAF == 1) then
               InputFileData%ReListBL(iRe) = TempRe  * 1.e+006
               
               if (iRe > 1) then
                  if (InputFileData%ReListBL(iRe) <= InputFileData%ReListBL(iRe-1) ) then
                     call SetErrStat(ErrID_Fatal,'All aeroacoustics BL tables must have Reynolds Numbers entered in increasing order.',ErrStat, ErrMsg, RoutineName)
                     call cleanup()
                     return
                  end if
               end if
               
            else
               if ( nRe > 1 .AND. .NOT. EqualRealNos(InputFileData%ReListBL(iRe), TempRe * 1.e+006 ) ) then
                  call SetErrStat(ErrID_Fatal,'All aeroacoustics BL tables must have the same Reynolds Numbers.',ErrStat, ErrMsg, RoutineName)
                  call cleanup()
                  return
               end if
            end if

            CALL ReadCom(UnIn, FileName, "aoa     Ue_Vinf_SS     Ue_Vinf_PS      Dstar_SS     Dstar_PS   Theta_SS   Theta_PS    Cf_SS   Cf_PS", ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
            CALL ReadCom(UnIn, FileName, "(deg)   (-)            (-)             (-)          (-)        (-)        (-)         (-)     (-)",   ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

            do iAoA=1,nAoA
                CALL ReadAry( UnIn, FileName, Buffer, SIZE(Buffer), 'BL Table Line '//Num2LStr(iAoA+8), 'BL Table for suction and pressure', ErrStat2, ErrMsg2, UnEc) ! From NWTC_Library
                if(Failed()) return
               
                Buffer(1) = Buffer(1)*D2R ! convert to radians
                call MPi2Pi( Buffer(1)  ) ! convert to radians between -pi and pi
                Buffer(1) = Buffer(1)*R2D ! convert back to degrees
               
                if (iAF == 1 .AND. iRe == 1) then
                   InputFileData%AoAListBL(iAoA) = Buffer( 1) ! AoA in degrees
                   
                  if (iAoA > 1) then

                     if (InputFileData%AoAListBL(iAoA) <= InputFileData%AoAListBL(iAoA-1) ) then
                        call SetErrStat(ErrID_Fatal,'All aeroacoustics BL tables angles of attack must be entered in increasing order.',ErrStat, ErrMsg, RoutineName)
                        call cleanup()
                        return
                     end if
                  end if

                else
                   if ( .NOT. EqualRealNos(InputFileData%AoAListBL(iAoA), Buffer( 1) ) ) then
                      call SetErrStat(ErrID_Fatal,'All aeroacoustics BL tables must have the same angles of attack.',ErrStat, ErrMsg, RoutineName)
                      call cleanup()
                      return
                   end if
                end if
                
                InputFileData%Suct_EdgeVelRat(iAoA,iRe,iAF)= Buffer(2) ! EdgeVelRat1 Suction
                InputFileData%Pres_EdgeVelRat(iAoA,iRe,iAF)= Buffer(3) ! EdgeVelRat2 Pressure
                InputFileData%Suct_DispThick (iAoA,iRe,iAF)= Buffer(4) ! dStarAll1 Suction
                InputFileData%Pres_DispThick (iAoA,iRe,iAF)= Buffer(5) ! dStarAll2 Pressure
                InputFileData%Suct_BLThick   (iAoA,iRe,iAF)= Buffer(6) ! d99All1 Suction
                InputFileData%Pres_BLThick   (iAoA,iRe,iAF)= Buffer(7) ! d99All2 Pressure
                InputFileData%Suct_Cf        (iAoA,iRe,iAF)= Buffer(8) ! CfAll1 Suction
                InputFileData%Pres_Cf        (iAoA,iRe,iAF)= Buffer(9) ! CfAll2 Pressure
            enddo
        enddo


        if (InputFileData%IBLUNT==IBLUNT_BPM) then
            call ReadCom(UnIn, FileName, 'Comment' , ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
            call ReadCom(UnIn, FileName, 'Comment' , ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
            call ReadVar(UnIn, FileName, InputFileData%BladeProps(iAF)%TEAngle, 'TEAngle', 'TE Angle',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
            call ReadVar(UnIn, FileName, InputFileData%BladeProps(iAF)%TEThick, 'TEThick', 'TE Thick',ErrStat2, ErrMsg2, UnEc ); if(Failed()) return
        else
            InputFileData%BladeProps(iAF)%TEAngle = 0._ReKi
            InputFileData%BladeProps(iAF)%TEThick = 0._ReKi
        endif
        
        call Cleanup()

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
!> This routine validates the inputs from the AeroDyn input files.
SUBROUTINE ValidateInputData( InputFileData, NumBl, ErrStat, ErrMsg )
   type(AA_InputFile),       intent(in)     :: InputFileData                       !< All the data in the AeroDyn input file
   integer(IntKi),           intent(in)     :: NumBl                               !< Number of blades
   integer(IntKi),           intent(out)    :: ErrStat                             !< Error status
   character(*),             intent(out)    :: ErrMsg                              !< Error message
   ! local variables
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
      call SetErrStat ( ErrID_Fatal, 'TICalcMeth must be '//trim(num2lstr(TICalc_Every))//' (TICalc automatic)  or '//&
           trim(num2lstr(TICalc_Interp))//' (TICalcMeth interp).', ErrStat, ErrMsg, RoutineName )
   end if

   if (InputFileData%X_BLMethod /= X_BLMethod_BPM .and. InputFileData%X_BLMethod /= X_BLMethod_Tables) then
      call SetErrStat ( ErrID_Fatal, 'X_BLMethod must be '//trim(num2lstr(X_BLMethod_BPM))//' X_BLMethod_ with BPM or '//&
           trim(num2lstr(X_BLMethod_Tables))//' (X_BLMethod with BL tables).', ErrStat, ErrMsg, RoutineName )
   end if
   
   if (InputFileData%NrObsLoc <= 0.0) call SetErrStat ( ErrID_Fatal, 'Number of Observer Locations should be greater than zero', ErrStat, ErrMsg, RoutineName )
   
   if (InputFileData%NrOutFile /= 1 .and. InputFileData%NrOutFile /= 2 .and. InputFileData%NrOutFile /= 3 .and. InputFileData%NrOutFile /= 4) then
      call SetErrStat ( ErrID_Fatal, ' NrOutFile must be 1 or 2 or 3 or 4', ErrStat, ErrMsg, RoutineName )
   end if

   if (InputFileData%AA_Bl_Prcntge > 100.0 .or. InputFileData%AA_Bl_Prcntge < 0.0) then
      call SetErrStat ( ErrID_Fatal, ' AA_Bl_Prcntge must be between 0 and 100%', ErrStat, ErrMsg, RoutineName )
   end if
   
END SUBROUTINE ValidateInputData

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets the initialization output data structure, which contains data to be returned to the calling program (e.g.,
!! FAST or AeroAcoustics_Driver)
subroutine AA_SetInitOut(p, InitOut, errStat, errMsg)
    type(AA_InitOutputType),       intent(  out)  :: InitOut          ! output data
    type(AA_ParameterType),        intent(in   )  :: p                ! Parameters
    integer(IntKi),                intent(  out)  :: errStat          ! Error status of the operation
    character(*),                  intent(  out)  :: errMsg           ! Error message if ErrStat /= ErrID_None
    ! Local variables
    integer(intKi)                               :: ErrStat2          ! temporary Error status
    character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
    character(*), parameter                      :: RoutineName = 'AA_SetInitOut'
    integer(IntKi)                               :: i, j, k,oi
    CHARACTER(16)                                :: ChanBNPrefix                     ! Name prefix (AeroB#_Z######y_)
    CHARACTER(6)                                 :: TmpChar                          ! Temporary char array to hold the node digits (2 places only!!!!)
    
    ! Initialize variables for this routine
    
    errStat = ErrID_None
    errMsg  = ""
    
    ! FIRST  FILE HEADER,UNIT
    call AllocAry(InitOut%WriteOutputHdr, p%numOutsAll(1), 'WriteOutputHdr', errStat2, errMsg2); if(Failed()) return
    call AllocAry(InitOut%WriteOutputUnt, p%numOutsAll(1), 'WriteOutputUnt', errStat2, errMsg2); if(Failed()) return
    do j=1,p%NrObsLoc
        InitOut%WriteOutputHdr(j)="Obs"//trim(num2lstr(j))
        InitOut%WriteOutputUnt(j) = "OASPL"
    enddo

    ! SECOND FILE HEADER,UNIT
    if (p%NrOutFile>1) then
       i=0
       call AllocAry(InitOut%WriteOutputHdrforPE, p%numOutsAll(2), 'WriteOutputHdrPE', errStat2, errMsg2); if(Failed()) return
       call AllocAry(InitOut%WriteOutputUntforPE, p%numOutsAll(2), 'WriteOutputUntPE', errStat2, errMsg2); if(Failed()) return
       do j=1,p%NrObsLoc
          do k=1,size(p%FreqList)
             i=i+1
             InitOut%WriteOutputHdrforPE(i) = "Obs"//trim(num2lstr(j))//"_Freq"//trim(num2lstr(p%FreqList(k)))
          end do
       enddo
       if(p%aweightflag) then  ! whole array
          InitOut%WriteOutputUntforPE = "SPL_A"
       else
          InitOut%WriteOutputUntforPE = "SPL"
       endif
       
       
       if (p%NrOutFile>2) then
          ! THIRD FILE HEADER,UNIT
          call AllocAry(InitOut%WriteOutputHdrSep, p%numOutsAll(3), 'WriteOutputHdrSep', errStat2, errMsg2); if(Failed()) return
          call AllocAry(InitOut%WriteOutputUntSep, p%numOutsAll(3), 'WriteOutputUntSep', errStat2, errMsg2); if(Failed()) return
          i=0
          do j=1,p%NrObsLoc
             do k=1,size(p%FreqList)
                 do oi=1,nNoiseMechanism
                    i=i+1
                    InitOut%WriteOutputHdrSep(i) = "Obs"//trim(num2lstr(j))//"_Freq"//trim(num2lstr(p%FreqList(k)))//"_Type"//trim(num2lstr(oi))
                    InitOut%WriteOutputHdrSep(i)=trim(InitOut%WriteOutputHdrSep(i))
                 enddo
             enddo
          enddo
          if(p%aweightflag) then ! whole array
             InitOut%WriteOutputUntSep = "SPL_A"
          else
             InitOut%WriteOutputUntSep = "SPL"
          endif

          if (p%NrOutFile>3) then
             ! FOURTH FILE HEADER,UNIT
             call AllocAry(InitOut%WriteOutputHdrNodes,p%numOutsAll(4), 'InitOut%WriteOutputHdrNodes', errStat2, errMsg2); if(Failed()) return
             call AllocAry(InitOut%WriteOutputUntNodes,p%numOutsAll(4), 'InitOut%WriteOutputUntNodes', errStat2, errMsg2); if(Failed()) return
             i=0
             do oi = 1,p%numBlades
                 do k = 1,p%NumBlNds
                     do j = 1,p%NrObsLoc
                         i=i+1
                         ChanBNPrefix = setChannelBldNdPrefix(oi,k)
                         InitOut%WriteOutputHdrNodes(i) = trim(ChanBNPrefix)//"Obs"//trim(num2lstr(j))
                     enddo
                 enddo
             enddo
             InitOut%WriteOutputUntNodes = "SPL"
             
         end if ! file 4
      end if ! file 3
    end if ! file 2
    
contains

   function setChannelBldNdPrefix(IdxBlade,IdxNode) result(ChanPrefix)
      INTEGER(IntKi), intent(in)  :: IdxBlade                         ! Counter to which blade we are on
      INTEGER(IntKi), intent(in)  :: IdxNode                          ! Counter to the blade node we ae on
      CHARACTER(16)               :: ChanPrefix                       ! Name prefix (AeroB#_Z######y_)
    
         ! Create the name prefix:
      WRITE (TmpChar,'(I3.3)')  IdxNode         ! 3 digit number
      ChanPrefix = 'AB' // TRIM(Num2LStr(IdxBlade)) // 'N' // TRIM(TmpChar)
    end function

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

   ! Local variables
   integer(IntKi) :: i
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'AA_InitializeOutputFile'
   

   p%unOutFile = -1
   ErrStat  = ErrID_None
   ErrMsg  = ""

   ErrStat2 = ErrID_None
   ErrMsg2 = ""
      
   
   ! FIRST FILE
   call WriteHeader(1,InitOut%WriteOutputHdr,InitOut%WriteOutputUnt,'Number of observers      :'//TRIM(num2lstr(p%NrObsLoc)) )
      if (Failed()) return
   
   ! SECOND FILE
   IF (InputFileData%NrOutFile > 1) THEN
      call WriteHeader(2,InitOut%WriteOutputHdrforPE,InitOut%WriteOutputUntforPE,'Number of observers      :'//TRIM(num2lstr(p%NrObsLoc))//';  Number of frequencies    :'//TRIM(num2lstr(size(p%FreqList)))  )
         if (Failed()) return

      ! THIRD FILE
      IF (InputFileData%NrOutFile > 2) THEN
         call WriteHeader(3,InitOut%WriteOutputHdrSep,InitOut%WriteOutputUntSep,'Number of observers      :'//TRIM(num2lstr(p%NrObsLoc))//';  Number of frequencies    :'//TRIM(num2lstr(size(p%FreqList)))//"; 1-LBL 2-TBLPres 3-TBLSuc 4-Sep  5-BLUNT 6-TIP 7-Inflow"  )
            if (Failed()) return
         
         ! FOURTH FILE
         IF (InputFileData%NrOutFile > 3) THEN
            call WriteHeader(4,InitOut%WriteOutputHdrNodes,InitOut%WriteOutputUntNodes,'Number of observers      :'//TRIM(num2lstr(p%NrObsLoc))//';  Number of blades         :'//TRIM(num2lstr(p%numBlades))//';     Number of nodes per blade:'//TRIM(num2lstr(p%NumBlNds))  )
               if (Failed()) return
         ENDIF
      
      ENDIF
   ENDIF

   
contains
   !-------------------------------------------------------------------------------------------------
   subroutine WriteHeader(iFile,WrOutHdr,WrOutUnt,LineTxt)
      integer(IntKi), intent(in) :: iFile
      CHARACTER(*), intent(in)     :: WrOutHdr(:)
      CHARACTER(*), intent(in)     :: WrOutUnt(:)
      character(*), intent(in)   :: LineTxt   ! text description to write to line 3 of the file
      
      call GetNewUnit( p%unOutFile(iFile), ErrStat2, ErrMsg2 )
      if (Failed()) return

      call OpenFOutFile ( p%unOutFile(iFile), trim(InputFileData%AAOutFile(iFile)), ErrStat2, ErrMsg2 )
      if (Failed()) return

      write (p%unOutFile(iFile),'(A)')  ''
      write (p%unOutFile(iFile),'(A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using AA '//trim(GetNVD(AA_ver))
      write (p%unOutFile(iFile),'(A)')  ''
      write( p%unOutFile(iFile),'(A)')  TRIM(LineTxt)
      write (p%unOutFile(iFile),'(A)')  'Description from AA input file, line2: '//trim(InputFileData%FTitle)
      write (p%unOutFile(iFile),'(A)')  ''

      !......................................................
      ! Write the names of the output parameters on one line: line 7
      !......................................................
      call WrFileNR ( p%unOutFile(iFile), '     Time           ' )
      do i=1,p%NumOutsAll(iFile)
         call WrFileNR ( p%unOutFile(iFile), delim//WrOutHdr(i) )
      end do ! i
      write (p%unOutFile(iFile),'()')
      
      !......................................................
      ! Write the units of the output parameters on one line: line 8
      !......................................................
      call WrFileNR ( p%unOutFile(iFile), '      (s)           ' )
      do i=1,p%NumOutsAll(iFile)
         call WrFileNR ( p%unOutFile(iFile), delim//WrOutUnt(i) )
      end do ! i
      write (p%unOutFile(iFile),'()')
      
   end subroutine
   !-------------------------------------------------------------------------------------------------
   subroutine cleanup()

   end subroutine
   !-------------------------------------------------------------------------------------------------
   LOGICAL function Failed()
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      Failed = ErrStat >= AbortErrLev
      if (Failed) call cleanup()
   end function
   !-------------------------------------------------------------------------------------------------
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

    errStat = ErrID_None
    errMsg  = ''
    
    frmt = '"'//delim//'"'//trim(p%outFmt)      ! format for array elements from individual modules
    write( tmpStr, '(F15.4)' ) t
    
    ! FIRST FILE
    IF (p%NrOutFile .gt. 0) THEN
       ! time
       call WrFileNR( p%unOutFile(1), tmpStr )
       call WrNumAryFileNR ( p%unOutFile(1), y%WriteOutput,  frmt, errStat, errMsg )
       if ( errStat >= AbortErrLev ) return
       ! write a new line (advance to the next line)
       write (p%unOutFile(1),'()')
    ENDIF

    !! SECOND FILE
    IF (p%NrOutFile .gt. 1) THEN
       ! time
       call WrFileNR( p%unOutFile(2), tmpStr )
       call WrNumAryFileNR ( p%unOutFile(2), y%WriteOutputforPE,  frmt, errStat, errMsg )
       if ( errStat >= AbortErrLev ) return
       ! write a new line (advance to the next line)
       write (p%unOutFile(2),'()')
    ENDIF
    
    ! THIRD FILE
    IF (p%NrOutFile .gt. 2) THEN
       ! time
       call WrFileNR( p%unOutFile(3), tmpStr )
       call WrNumAryFileNR ( p%unOutFile(3), y%WriteOutputSep,  frmt, errStat, errMsg )
       if ( errStat >= AbortErrLev ) return
       ! write a new line (advance to the next line)
       write (p%unOutFile(3),'()')
    ENDIF
    
    ! Fourth FILE
    IF (p%NrOutFile .gt. 3) THEN
       ! time
       call WrFileNR( p%unOutFile(4), tmpStr )
       call WrNumAryFileNR ( p%unOutFile(4), y%WriteOutputNodes,  frmt, errStat, errMsg )
       if ( errStat >= AbortErrLev ) return
       ! write a new line (advance to the next line)
       write (p%unOutFile(4),'()')
    ENDIF
end subroutine AA_WriteOutputLine
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Calc_WriteOutput( p, m, y, ErrStat, ErrMsg )
   TYPE(AA_ParameterType),    INTENT(IN   )  :: p                                 ! The module parameters
   TYPE(AA_MiscVarType),      INTENT(INOUT)  :: m                                 ! misc variables
   TYPE(AA_OutputType),       INTENT(INOUT)  :: y                                 ! outputs
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat                           ! The error status code
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg                            ! The error message, if an error occurred
   ! local variables
   CHARACTER(*), PARAMETER                   :: RoutineName = 'Calc_WriteOutput'
!   INTEGER(intKi)                            :: ErrStat2
!   CHARACTER(ErrMsgLen)                      :: ErrMsg2
   INTEGER(IntKi)                            :: j,k,counter,i,oi,III
   ! start routine:
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! FOR THE FIRST OUTPUT FILE
   IF (p%NrOutFile .gt. 0) THEN
      y%WriteOutput(1:p%NrObsLoc)=m%DirectiviOutput

      ! FOR THE SECOND OUTPUT FILE
      IF (p%NrOutFile .gt. 1) THEN
         counter=0
         DO K = 1,p%NrObsLoc
            DO III = 1,size(p%FreqList)
               counter=counter+1
               y%WriteOutputforPE(counter) = m%PtotalFreq(III,K)
            END DO !
         END DO !

         ! FOR THE THIRD OUTPUT FILE
         IF (p%NrOutFile .gt. 2) THEN
            counter=0
            do K = 1,p%NrObsLoc
               do III = 1,size(p%FreqList)
                  do oi=1,nNoiseMechanism
                     counter=counter+1
                     y%WriteOutputSep(counter) = m%SumSpecNoiseSep(oi,III,K)
                  enddo
               enddo
            enddo

             ! FOR THE FOURTH OUTPUT FILE
             IF (p%NrOutFile .gt. 3) THEN
                 counter=0
                 DO I = 1,p%numBlades
                     DO J = 1,p%NumBlNds
                         DO K = 1,p%NrObsLoc
                             counter=counter+1
                             y%WriteOutputNodes(counter) = m%OASPL(K,J,I)
                         END DO !
                     END DO !
                 ENDDO
             ENDIF
         ENDIF
      ENDIF
   ENDIF
    
END SUBROUTINE Calc_WriteOutput
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE AeroAcoustics_IO
