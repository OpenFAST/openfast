!**********************************************************************************************************************************
! The ExtPtfm_MCKF.f90, ExtPtfm_MCKF_IO.f90 and  ExtPtfm_MCKF_Types.f90 make up the ExtPtfm_MCKF module of the
! FAST Modularization Framework. ExtPtfm_MCKF_Types is auto-generated based on FAST_Registry.txt.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2016 National Renewable Energy Laboratory
!
!    This file is part of ExtPtfm_MCKF.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
!> This module contains definitions of compile-time PARAMETERS for the ExtPtfm_MCKF module.
!! Every variable defined here MUST have the PARAMETER attribute.
MODULE ExtPtfm_MCKF_Parameters
   USE NWTC_Library

   TYPE(ProgDesc), PARAMETER :: ExtPtfm_Ver = ProgDesc( 'ExtPtfm_MCKF', '', '' ) !< module date/version information

   CHARACTER(len=4), DIMENSION(3), PARAMETER :: StrIntMethod = (/'RK4 ','AB4 ','ABM4'/)

   ! Variables for output channels
   INTEGER(IntKi), PARAMETER :: FILEFORMAT_GUYANASCII = 0
   INTEGER(IntKi), PARAMETER :: FILEFORMAT_FLEXASCII  = 1


   ! Variables for output channels
   INTEGER(IntKi), PARAMETER :: MaxOutChs   = 9 + 3*200 ! Maximum number of output channels
                                                        ! Harcoded to outputs of 200 CB modes
   INTEGER(IntKi), PARAMETER :: OutStrLenM1 = ChanLen - 1
   INTEGER(IntKi), PARAMETER :: ID_Time     = 0
   INTEGER(IntKi), PARAMETER :: ID_PtfFx    = 1
   INTEGER(IntKi), PARAMETER :: ID_PtfFy    = 2
   INTEGER(IntKi), PARAMETER :: ID_PtfFz    = 3
   INTEGER(IntKi), PARAMETER :: ID_PtfMx    = 4
   INTEGER(IntKi), PARAMETER :: ID_PtfMy    = 5
   INTEGER(IntKi), PARAMETER :: ID_PtfMz    = 6
   INTEGER(IntKi), PARAMETER :: ID_WaveElev = 7
   INTEGER(IntKi), PARAMETER :: ID_QStart   = 8
END MODULE ExtPtfm_MCKF_Parameters

!**********************************************************************************************************************************
!> This module contains file I/O routines and data validation routines.
MODULE ExtPtfm_MCKF_IO

   USE ExtPtfm_MCKF_Parameters
   USE ExtPtfm_MCKF_Types

   IMPLICIT NONE
   private

   public :: ReadPrimaryFile
   public :: SetOutParam 
   public :: SetOutParamLin
   public :: ExtPtfm_PrintSum
   
CONTAINS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Helper functions for the module
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets the error status and error message for a routine, it's a simplified version of SetErrStat from NWTC_Library
subroutine SetErrStatSimple(ErrStat, ErrMess, RoutineName, LineNumber)
  INTEGER(IntKi), INTENT(INOUT)        :: ErrStat      ! Error status of the operation
  CHARACTER(*),   INTENT(INOUT)        :: ErrMess      ! Error message if ErrStat /= ErrID_None
  CHARACTER(*),   INTENT(IN   )        :: RoutineName  ! Name of the routine error occurred in
  INTEGER(IntKi), INTENT(IN), OPTIONAL :: LineNumber   ! Line of input file 
  if (ErrStat /= ErrID_None) then
     ErrMess = TRIM(RoutineName)//':'//TRIM(ErrMess)
     if (present(LineNumber)) then
         ErrMess = TRIM(ErrMess)//' Line: '//TRIM(Num2LStr(LineNumber))//'.'
     endif
  end if
end subroutine SetErrStatSimple

subroutine disp2r8(u,varname,a)
    integer,intent(in) ::u
    character(len=*),intent(in)::varname
    real(ReKi),intent(in),dimension(:,:) ::a
    integer :: n, m,i
    character(len=20) :: fmt
    character(len=*),parameter :: RFMT='EN13.3E2'
    n=size(a,1)
    m=size(a,2)
    if (n>0 .and. m>0) then
        write(u,"(A,A)") varname,"=["
        write(fmt,*) m
        do i=1,n-1
            write(u,"("//adjustl(fmt)//RFMT//")") a(i,:)
        enddo
        i=n
        write(u,"("//trim(fmt)//RFMT//",A)") a(i,:), "  ];"
    else
        write(u,'(A,A)') varname,'=[];'
    endif
end subroutine

!----------------------------------------------------------------------------------------------------------------------------------
!> Helper functions to read primary file
real(ReKi) function ReadFloatFromStr(s, VarName, iStat, Msg ) result(myfloat)
   character(len=*), intent(in)    :: s
   character(len=*), intent(in)    :: VarName
   character(len=*), intent(inout) :: Msg
   integer, intent(out)            :: iStat
   read(s,*, iostat=iStat ) myfloat 
   if (iStat /= 0) then
      iStat=ErrID_Fatal
      Msg = trim(Msg)//'Error extracting float while reading '//VarName
   endif
end function ReadFloatFromStr
integer function ReadIntFromStr(s, VarName, iStat, Msg ) result(myint)
   character(len=*), intent(in)    :: s
   character(len=*), intent(in)    :: VarName
   character(len=*), intent(inout) :: Msg
   integer, intent(out)            :: iStat
   read(s,*, iostat=iStat ) myint 
   if (iStat /= 0) then
      iStat=ErrID_Fatal
      Msg = trim(Msg)//'Error extracting integer while reading '//VarName
   endif
end function ReadIntFromStr
subroutine ReadRealMatrix(fid, FileName, Mat, VarName, nLines,nRows, iStat, Msg, iLine )
   integer, intent(in)                     :: fid
   real(ReKi), dimension(:,:), allocatable :: Mat
   character(len=*), intent(in)            :: FileName
   character(len=*), intent(in)            :: VarName
   integer, intent(in)                     :: nLines
   integer, intent(in)                     :: nRows
   integer, intent(out)                    :: iStat
   integer, intent(inout)                  :: iLine
   character(len=*), intent(inout)         :: Msg
   ! local variables
   integer :: i
   call allocAry( Mat, nLines, nRows, VarName,  iStat, Msg); 
   if (iStat /= 0) return
   !Read Stiffness
   DO I =1,nLines
      iLine=iLine+1
      ! TODO use ReadCAryFromStr when available in the NWTCIO, it performs more checks
      CALL ReadAry( fid, FileName, Mat(I,:), nRows, trim(VarName)//' Line '//Num2LStr(iLine), VarName, iStat, Msg)
      if (iStat /= 0) return
   ENDDO
end subroutine


SUBROUTINE SetOutParam(OutList, NumOuts_in, p, ErrStat, ErrMsg )
! This routine checks to see if any requested output channel names (stored in the OutList(:)) are invalid. It returns a 
! warning if any of the channels are not available outputs from the module.
!  It assigns the settings for OutParam(:) (i.e, the index, name, and units of the output channels, WriteOutput(:)).
!  the sign is set to 0 if the channel is invalid.
! It sets assumes the value p%NumOuts has been set before this routine has been called, and it sets the values of p%OutParam here.
!..................................................................................................................................
   CHARACTER(ChanLen),           INTENT(IN)     :: OutList(:)         !< The list out user-requested outputs
   INTEGER(IntKi),               INTENT(IN)     :: NumOuts_in         !< Effective number of output channels
   TYPE(ExtPtfm_ParameterType),  INTENT(INOUT)  :: p                  !< The module parameters
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat            !< The error status code
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg             !< The error message, if an error occurred
   ! Local variables
   INTEGER                      :: ErrStat2                                        ! temporary (local) error status
   INTEGER                      :: I                                               ! Generic loop-counting index
   INTEGER                      :: INDX                                            ! Index for valid arrays
   CHARACTER(ChanLen)           :: OutListTmp                                      ! A string to temporarily hold OutList(I)
   CHARACTER(*), PARAMETER      :: RoutineName = "SetOutParam"

   CHARACTER(OutStrLenM1), PARAMETER  :: ValidParamAry(7) =  (/ & ! This lists the names of the allowed parameters, which must be sorted alphabetically
                               "INTRFFX  ","INTRFFY  ","INTRFFZ  ","INTRFMX  ","INTRFMY  ","INTRFMZ  ","WAVELEV  "/) 
   CHARACTER(OutStrLenM1), PARAMETER :: ParamUnitsAry(7) =  (/ &                     ! This lists the units corresponding to the allowed parameters
                               "(N)      ","(N)      ","(N)      ","(Nm)     ","(Nm)     ","(Nm)     ","(m)      "/)
   INTEGER(IntKi), PARAMETER :: ParamIndxAry(7) =  (/ &                            ! This lists the index into AllOuts(:) of the allowed parameters ValidParamAry(:)
                              ID_PtfFx, ID_PtfFy, ID_PtfFz, ID_PtfMx, ID_PtfMy, ID_PtfMz, ID_WaveElev /)
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   p%NumOuts = NumOuts_in
   allocate(p%OutParam(0:p%NumOuts) , stat=ErrStat )
   if ( ErrStat /= 0_IntKi )  THEN
      CALL SetErrStat(ErrID_Fatal,"Error allocating memory for the InflowWind OutParam array.", ErrStat, ErrMsg, RoutineName)
      return
   endif

   ! Set index, name, and units for the time output channel:
   p%OutParam(0)%Indx  = 0
   p%OutParam(0)%Name  = "Time"    ! OutParam(0) is the time channel by default.
   p%OutParam(0)%Units = "(s)"
   p%OutParam(0)%SignM = 1

  ! Set index, name, and units for all of the output channels.
   do I = 1,p%NumOuts
      p%OutParam(I)%Name  = OutList(I)
      OutListTmp          = OutList(I)
      p%OutParam(I)%Indx   = 0
      p%OutParam(I)%Units  = "(NA)"
      CALL Conv2UC( OutListTmp )  ! Convert OutListTmp to upper case
      ! Reverse the sign of the channel if the prefix is "-", "_" or "M"
      if  ( index( "-_M", OutListTmp(1:1) ) > 0 ) then
         p%OutParam(I)%SignM = -1 
         OutListTmp          = OutListTmp(2:)
      else
         p%OutParam(I)%SignM = 1
      end if
      ! Find the index of the channel in the AllOut list
      Indx = IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )
      if (Indx>0) then
          p%OutParam(I)%Indx  = ParamIndxAry(Indx)
          p%OutParam(I)%Units = ParamUnitsAry(Indx)
      else if (index(OutListTmp,'CBQ_') > 0 ) then
          call setDOFChannel(5,ID_QStart+0*p%nCB-1); if(Failed()) return
      else if (index(OutListTmp,'CBQD_') > 0 ) then
          call setDOFChannel(6,ID_QStart+1*p%nCB-1); if(Failed()) return
      else if (index(OutListTmp,'CBF_') > 0 ) then
          call setDOFChannel(5,ID_QStart+2*p%nCB-1); if(Failed()) return
      else
          call setInvalidChannel() ! INVALID
      endif
      !write(*,*) p%OutParam(I)%Name, p%OutParam(I)%Indx, p%OutParam(I)%Units
   end do
   return
contains
    logical function Failed()
        CALL SetErrStatSimple(ErrStat, ErrMsg, 'ExtPtfm_SetOutParam')
        Failed =  ErrStat >= AbortErrLev
    end function Failed
    subroutine setDOFChannel(nCharBefore,nOffset)
        !> Sets channel when the channel name has the form "YYYY_XXX" where XXX is a DOF number
        integer, intent(in) :: nCharBefore !< Number of characters to ignore in OutListTmp
        integer, intent(in) :: nOffset     !< Index offset to add to iDOF
        integer             :: idof ! index of CB DOF extracted from 
        iDOF = ReadIntFromStr(OutListTmp(nCharBefore:), 'Output channel '//trim(OutList(I)), ErrStat, ErrMsg);
        if(ErrStat/=0) return
        if ((iDOF> p%nCB) .or. (iDOF<1)) then
            call setInvalidChannel() ! INVALID
        else
            p%OutParam(I)%Indx  = nOffset+iDOF
            p%OutParam(I)%Units = '(-)'
        endif
    end subroutine
    subroutine setInvalidChannel()
        ! If a selected output channel is not available by this module set ErrStat = ErrID_Warn.
        p%OutParam(I)%Units = "INVALID"
        p%OutParam(I)%Indx = 0
        call SetErrStat(ErrID_Warn, TRIM(p%OutParam(I)%Name)//" is not an available output channel.",ErrStat,ErrMsg,'ExtPtfm_SetOutParam')
        write(*,*)TRIM(p%OutParam(I)%Name)//" is not an available output channel."
    end subroutine
END SUBROUTINE SetOutParam
!----------------------------------------------------------------------------------------------------------------------------------
!..................................................................................................................................
!> This routine checks to see if any requested output channel names are to be output in linearization analysis.
!! note that we output all WriteOutput values and assume that none of them depend on inputs (so I don't need this mapping any more)
SUBROUTINE SetOutParamLin( p, ErrStat, ErrMsg )
   ! Passed variables
   TYPE(ExtPtfm_ParameterType),  INTENT(INOUT)  :: p                  !< The module parameters
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat            !< The error status code
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg             !< The error message, if an error occurred
   ! Local variables
   !INTEGER                   :: ErrStat2                                        ! temporary (local) error status
   !INTEGER                   :: I                                               ! Generic loop-counting index
   !INTEGER                   :: J                                               ! Generic loop-counting index
   !CHARACTER(ErrMsgLen)      :: ErrMsg2
   !CHARACTER(*), PARAMETER   :: RoutineName = "SetOutParamLin"
   ErrStat = ErrID_None
   ErrMsg  = ""
!    call AllocAry(p%OutParamLinIndx, 2, p%NumOuts, 'OutParamLinIndx', ErrStat2, ErrMsg2)
!    call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
!    if (ErrStat >= AbortErrLev) return
   !do i = 1,p%NumOuts
   !   if (p%OutParam(i)%SignM /= 0 ) then
   !      do j=1,size(WindVelX)
   !         if ( p%OutParam(i)%Indx == WindVelX(j) ) then
   !            p%OutParamLinIndx(1,i) = j
   !            p%OutParamLinIndx(2,i) = 1
   !            exit !exit j loop; move to next parameter
   !         elseif ( p%OutParam(i)%Indx == WindVelY(j) ) then
   !            p%OutParamLinIndx(1,i) = j
   !            p%OutParamLinIndx(2,i) = 2
   !            exit !exit j loop; move to next parameter
   !         elseif ( p%OutParam(i)%Indx == WindVelZ(j) ) then
   !            p%OutParamLinIndx(1,i) = j
   !            p%OutParamLinIndx(2,i) = 3
   !            exit !exit j loop; move to next parameter
   !         end if
   !      end do
   !      
   !   end if      
   !end do
END SUBROUTINE SetOutParamLin
!----------------------------------------------------------------------------------------------------------------------------------
!> Checks that all inputs were correctly read
subroutine CheckAllInputsRead(p,ErrStat,ErrMsg)
    TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p              !< All the parameter matrices stored in this input file
    INTEGER(IntKi),              INTENT(OUT)   :: ErrStat        !< Error status                              
    CHARACTER(*),                INTENT(OUT)   :: ErrMsg         !< Error message
    ErrStat = ErrID_None
    ErrMsg  = ""
    if (ErrStat/=0) return
    if (p%nTot<0)                   then ; ErrStat=ErrID_Fatal; ErrMsg='The total number of DOF was not set'; endif
    if (.not.allocated(p%PtfmAM))   then ; ErrStat=ErrID_Fatal; ErrMsg='The mass matrix was not allocated.' ; endif
    if (.not.allocated(p%Stff))     then ; ErrStat=ErrID_Fatal; ErrMsg='The stiffness matrix was not allocated.' ; endif
    if (.not.allocated(p%Damp))     then ; ErrStat=ErrID_Fatal; ErrMsg='The damping matrix was not allocated.' ; endif
    if (.not.allocated(p%PtfmFt))   then ; ErrStat=ErrID_Fatal; ErrMsg='The loads were not allocated.';endif
    if (.not.allocated(p%PtfmFt_t)) then ; ErrStat=ErrID_Fatal; ErrMsg='The time vector was not allocated.'; endif
end subroutine CheckAllInputsRead
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadPrimaryFile(InputFile, p, OutFileRoot, InputFileData, ErrStat, ErrMsg)
!..................................................................................................................................
   ! Passed variables
   CHARACTER(*),                INTENT(IN)    :: InputFile      !< Name of the file containing the primary input data
   TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p              !< All the parameter matrices stored in this input file
   CHARACTER(*),                INTENT(IN)    :: OutFileRoot    !< The rootname of all the output files written by this routine.
   TYPE(ExtPtfm_InputFile),     INTENT(OUT)   :: InputFileData ! Data stored in the module's input file
   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat        !< Error status                              
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg         !< Error message
   ! Local variables:
   INTEGER(IntKi)                       :: I                    ! loop counter
   INTEGER(IntKi)                       :: UnIn                 ! Unit number for reading file
   INTEGER(IntKi)                       :: UnEc                 ! Unit number for echo
   INTEGER(IntKi)                       :: iLine                ! Current position in file
   CHARACTER(200)                       :: Line                 ! Temporary storage of a line from the input file (to compare with "default")
   CHARACTER(1024)                      :: PriPath              ! Path name of the primary file
   LOGICAL                              :: Echo
   ! --- Initialization
   ErrStat = ErrID_None
   ErrMsg  = ""
   Echo = .FALSE.
   UnEc = -1                             ! Echo file not opened, yet
   CALL GetPath(InputFile, PriPath)     ! Input files will be relative to the path where the primary input file is located.
   CALL AllocAry(InputFileData%OutList, MaxOutChs, "ExtPtfm Input File's Outlist", ErrStat, ErrMsg); if(Failed()) return
   
   ! Get an available unit number for the file.
   CALL GetNewUnit(UnIn, ErrStat, ErrMsg);              if(Failed()) return
   ! Open the Primary input file.
   CALL OpenFInpFile(UnIn, InputFile, ErrStat, ErrMsg); if(Failed()) return
   
   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file.
   ! If Echo is TRUE, rewind and write on the second try.
   I    = 1 ! the number of times we've read the file (used for the Echo variable)
   DO
       iLine=1
       !-------------------------- HEADER ---------------------------------------------
       CALL ReadCom(UnIn, InputFile, 'File Header', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return; 
       CALL ReadStr(UnIn, InputFile, Line, 'Header', 'File Header: File Description', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return 
       !---------------------- SIMULATION CONTROL --------------------------------------
       CALL ReadCom(UnIn, InputFile, 'Section Header: Simulation Control', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return
       ! Echo - Echo input to "<RootName>.ech".
       CALL ReadVar(UnIn, InputFile, Echo, 'Echo','Echo switch', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return
       IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop
       ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
       I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
       CALL OpenEcho(UnEc, TRIM(OutFileRoot)//'.ech', ErrStat, ErrMsg, ExtPtfm_Ver ); if(Failed()) return;
       IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(ExtPtfm_Ver%Name)//' primary input file "'//TRIM( InputFile )//'":'
       REWIND( UnIn, IOSTAT=ErrStat )
       IF (ErrStat /= 0_IntKi ) THEN
           CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".',ErrStat,ErrMsg,'ExtPtfm_ReadPrimaryFile' );
           IF (ErrStat >= AbortErrLev) RETURN
       END IF
   END DO

   IF (NWTC_VerboseLevel == NWTC_Verbose) THEN
      CALL WrScr(' Heading of the '//TRIM(ExtPtfm_Ver%Name)//' input file: ')
      CALL WrScr('   '//TRIM( Line ))
   END IF

   ! DT - Requested integration time for ElastoDyn (seconds):
   InputFileData%DT=-1
   CALL ReadVar( UnIn, InputFile, Line, "DT", "Integration time for ExtPtfm (s)", ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   CALL Conv2UC( Line )
   IF ( INDEX(Line, "DEFAULT" ) /= 1 ) THEN ! If it's not "default", read this variable; otherwise use the value already stored in InputFileData%DT
      READ(Line, *, IOSTAT=ErrStat) InputFileData%DT
      IF ( ErrStat /= 0 ) THEN
         CALL CheckIOS(ErrStat, InputFile, "DT", NumType, ErrStat, ErrMsg); if(Failed()) return
      END IF
   END IF
   ! Method - Integration method for loose coupling
   CALL ReadVar( UnIn, InputFile, InputFileData%IntMethod, "IntMethod", "Integration method for ExtPtfm {1: RK4, 2: AB4, or 3: ABM4}", ErrStat, ErrMsg, UnEc); if(LineFailed()) return

   !---------------------- REDUCTION INPUTS ---------------------------------------------------
   CALL ReadCom(UnIn, InputFile, 'Section Header: ReductionInputs', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   ! File Format switch
   CALL ReadVar(UnIn, InputFile, InputFileData%FileFormat, "FileFormat", "File format switch", ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   ! Reduction Filename
   CALL ReadVar(UnIn, InputFile, InputFileData%RedFile   , 'Red_FileName', 'Path containing Guyan/Craig-Bampton inputs', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   IF ( PathIsRelative(InputFileData%RedFile) ) InputFileData%RedFile = TRIM(PriPath)//TRIM(InputFileData%RedFile)
   CALL ReadVar(UnIn, InputFile, InputFileData%RedFileCst, 'RedCst_FileName', 'Path containing Guyan/Craig-Bampton constant inputs', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   IF ( PathIsRelative(InputFileData%RedFileCst) ) InputFileData%RedFileCst = TRIM(PriPath)//TRIM(InputFileData%RedFileCst)
   !---------------------- OUTPUT --------------------------------------------------
   CALL ReadCom(UnIn, InputFile, 'Section Header: Output', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   ! SumPrint - Print summary data to <RootName>.sum (flag):
   CALL ReadVar(UnIn, InputFile, InputFileData%SumPrint, "SumPrint", "Print summary data to <RootName>.sum (flag)", ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   ! OutFile - Switch to determine where output will be placed: (1: in module output file only; 2: in glue code output file only; 3: both) (-):
   CALL ReadVar(UnIn, InputFile, InputFileData%OutFile , "OutFile", "Switch to determine where output will be placed: (1: in module output file only; 2: in glue code output file only; 3: both) (-)", ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   ! TabDelim - Flag to cause tab-delimited text output (delimited by space otherwise) (flag):
   CALL ReadVar(UnIn, InputFile, InputFileData%TabDelim, "TabDelim", "Flag to cause tab-delimited text output (delimited by space otherwise) (flag)", ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   ! OutFmt - Format used for module's text tabular output (except time); resulting field should be 10 characters (-):
   CALL ReadVar(UnIn, InputFile, InputFileData%OutFmt  , "OutFmt", "Format used for module's text tabular output (except time); resulting field should be 10 characters (-)", ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   ! Tstart - Time to start module's tabular output (seconds):
   CALL ReadVar(UnIn, InputFile, InputFileData%Tstart  , "Tstart", "Time to start module's tabular output (seconds)", ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   !---------------------- OUTLIST  --------------------------------------------
   CALL ReadCom(UnIn, InputFile, 'Section Header: OutList', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   ! OutList - List of user-requested output channels (-):
   CALL ReadOutputList(UnIn, InputFile, InputFileData%OutList, InputFileData%NumOuts, 'OutList', "List of user-requested output channels", ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   !---------------------- END OF FILE -----------------------------------------
   call cleanup()

   ! --- Reading Reduced file
   call ReadReducedFile(InputFileData%RedFile, p, InputFileData%FileFormat, ErrStat, ErrMsg); if(Failed()) return;
   ! Checking that everyting was correctly read and set
   call CheckAllInputsRead(p, ErrStat, ErrMsg);  if(Failed()) return

   return

CONTAINS
    logical function Failed()
        CALL SetErrStatSimple(ErrStat, ErrMsg, 'ExtPtfm_ReadPrimaryFile')
        Failed =  ErrStat >= AbortErrLev
        if(Failed) call cleanup()
    end function Failed
    logical function LineFailed()
        CALL SetErrStatSimple(ErrStat, ErrMsg, 'ExtPtfm_ReadPrimaryFile',iLine)
        LineFailed =  ErrStat >= AbortErrLev
        if(LineFailed) call cleanup()
        iLine=iLine+1  ! Increase line number
    end function LineFailed
   subroutine cleanup()
        if (UnIn>0) close(UnIn)
        if (UnEc>0) close(UnEc)
   end subroutine cleanup
END SUBROUTINE ReadPrimaryFile
!..................................................................................................................................
SUBROUTINE ReadReducedFile( InputFile, p, FileFormat, ErrStat, ErrMsg )
!..................................................................................................................................
   ! Passed variables
   CHARACTER(*),                INTENT(IN)    :: InputFile                           !< Name of the file containing the primary input data
   TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p                                   !< All the parameter matrices stored in this input file
   INTEGER(IntKi),              INTENT(IN)    :: FileFormat                          !< File format for reduction inputs
   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                             !< Error status                              
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                              !< Error message
   ! Local variables:
   REAL(ReKi), dimension(:),allocatable :: TmpAry                                 ! temporary array for reading row from file
   INTEGER(IntKi)                       :: I                                         ! loop counter
   INTEGER(IntKi)                       :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)                       :: iLine                                     ! Current position in file
   CHARACTER(200)                       :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
   ErrStat = ErrID_None
   ErrMsg  = ""
   if     (FileFormat==FILEFORMAT_GUYANASCII) then
       call ReadGuyanASCII()
   elseif (FileFormat==FILEFORMAT_FLEXASCII) then
       call ReadFlexASCII()
   else
       call SetErrStat(ErrID_Fatal, 'FileFormat not implemented: '//trim(Num2LStr(FileFormat)), ErrStat, ErrMsg, 'ExtPtfm_ReadReducedFile')
       return
   endif
   ! --- The code below can detect between FlexASCII and GuyanASCII format by looking at the two first lines
   ! Get an available unit number for the file.
   !CALL GetNewUnit( UnIn, ErrStat, ErrMsg );               if(Failed()) return
   !! Open the Primary input file.
   !CALL OpenFInpFile ( UnIn, InputFile, ErrStat, ErrMsg ); if(Failed()) return
   !iLine=1
   !!-------------------------- Read the first two lines
   !CALL ReadStr( UnIn, InputFile, Line, 'Line'//Num2LStr(iLine), 'External Platform MCKF file', ErrStat, ErrMsg)
   !if(Failed()) return
   !iLine=iLine+1
   !CALL ReadStr( UnIn, InputFile, Line2, 'Line'//Num2LStr(iLine), 'External Platform MCKF file', ErrStat, ErrMsg)
   !if(Failed()) return
   !iLine=iLine+1
   !call CONV2UC(Line)
   !call CONV2UC(Line2)
   !call cleanup()
   !!-------------------------- Detecting file format
   !if (index(Line2,'#MASS')==1) then
   !    write(*,*) 'File detected as Guyan ASCII file format: '//trim(InputFile)
   !    call ReadGuyanASCII()
   !else if (index(Line2,'FLEX 5 FORMAT')>=1) then
   !    write(*,*) 'File detected as FLEX ASCII file format: '//trim(InputFile)
   !    call ReadFlexASCII()
   !endif

CONTAINS
    !> 
    logical function Failed()
        CALL SetErrStatSimple(ErrStat, ErrMsg, 'ExtPtfm_ReadReducedFile')
        Failed =  ErrStat >= AbortErrLev
        if(Failed) call cleanup()
    end function Failed
    !> 
    subroutine cleanup()
        close( UnIn )
        if (allocated(TmpAry)) deallocate(TmpAry)
    end subroutine cleanup

   !> Reads a FLEX ASCII file for Guyan or CraigBampton reductions
   SUBROUTINE ReadFlexASCII()
       REAL(ReKi) :: dt !< time step
       REAL(ReKi) :: T  !< total simulation time

       T=-1
       dt=-1
       ! Get an available unit number for the file.
       CALL GetNewUnit( UnIn, ErrStat, ErrMsg );            if ( ErrStat /= 0 ) return
       ! Open the Primary input file.
       CALL OpenFInpFile(UnIn, InputFile, ErrStat, ErrMsg); if ( ErrStat /= 0 ) return

       ! --- Reading file line by line
       ErrStat=0
       iLine=0
       do while (ErrStat==0)
           iLine=iLine+1
           read(UnIn,'(A)', iostat=ErrStat) Line
           if (ErrStat/=0) then
               if (ErrStat < 0) then
                   ErrStat=0 ! End of file is fine
               else
                   ErrMsg='Error while reading file '//trim(InputFile)// ' line '//Num2LStr(iLine)
               endif
               exit
           endif
           ! Line content is analyzed as case incensitive 
           call Conv2UC(Line)
           if (index(Line,'!DIMENSION')==1) then
               p%nTot = ReadIntFromStr(Line(12:), '`dimension`, file '//trim(InputFile)//', line '//Num2LStr(iLine), ErrStat, ErrMsg); if (ErrStat /= 0) exit
               p%nCB=p%nTot-6

           else if (index(Line,'!TIME INCREMENT IN SIMULATION:')==1) then
               dt =  ReadFloatFromStr(Line(31:), '`time increment`, file '//trim(InputFile)//', line '//Num2LStr(iLine), ErrStat, ErrMsg); if (ErrStat /= 0) exit

           else if (index(Line,'!TOTAL SIMULATION TIME IN FILE:')==1) then
               T =  ReadFloatFromStr(Line(32:), '`total simulation time`, file '//trim(InputFile)//', line '//Num2LStr(iLine), ErrStat, ErrMsg ); if (ErrStat /= 0) exit

           else if (index(Line,'!MASS MATRIX')==1) then
               iLine=iLine+1
               CALL ReadCom( UnIn, InputFile, 'Comment - Line '//Num2LStr(iLine), ErrStat, ErrMsg); if (ErrStat /= 0) exit
               if (p%nTot<0) exit
               call ReadRealMatrix(UnIn, InputFile, p%PtfmAM, 'Mass Matrix', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)

           else if (index(Line,'!STIFFNESS MATRIX')==1) then
               iLine=iLine+1
               CALL ReadCom( UnIn, InputFile, 'Comment - Line '//Num2LStr(iLine), ErrStat, ErrMsg);  if (ErrStat /= 0) exit
               if (p%nTot<0) exit
               call ReadRealMatrix(UnIn, InputFile, p%Stff, 'Stiffness Matrix', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)

           else if (index(Line,'!DAMPING MATRIX')==1) then
               iLine=iLine+1
               CALL ReadCom( UnIn, InputFile, 'Comment - Line '//Num2LStr(iLine), ErrStat, ErrMsg); if (ErrStat /= 0) exit
               if (p%nTot<0) exit
               call ReadRealMatrix(UnIn, InputFile, p%Damp, 'Damping Matrix', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)

           else if (index(Line,'!LOADING AND WAVE ELEVATION')==1) then
               iLine=iLine+1
               CALL ReadCom( UnIn, InputFile, 'Comment - Line '//Num2LStr(iLine), ErrStat, ErrMsg)
               if (ErrStat /= 0) exit
               p%nPtfmFt = nint(T/dt)+1
               if (p%nTot<0) exit
               call allocAry( p%PtfmFt,   max(1,p%nPtfmFt), p%nTot, 'p%PtfmFt',   ErrStat, ErrMsg); if (ErrStat /= 0) exit
               call allocAry( p%PtfmFt_t, max(1,p%nPtfmFt),         'p%PtfmFt_t', ErrStat, ErrMsg); if (ErrStat /= 0) exit
               if (p%nPtfmFt == 0) then
                  p%PtfmFt   = 0.0_ReKi
                  p%PtfmFt_t = 0.0_ReKi
                  p%nPtfmFt  = 1
               else
                  allocate(TmpAry(1:p%nTot+1))
                  do i=1,p%nPtfmFt
                     iLine=iLine+1
                     call ReadAry( UnIn, InputFile, TmpAry, p%nTot+1, 'PtfmFt - Line: '//Num2LStr(iLine)//' Value: '//trim(Num2LStr(i))//'/'//Num2LStr(p%nPtfmFt), 'PtfmFt time-history', ErrStat, ErrMsg)
                     if (ErrStat /= 0) exit
                     p%PtfmFt_t(i) = TmpAry(1)
                     p%PtfmFt(i,:) = TmpAry(2:p%nTot+1)
                  end do
               end if

           elseif (index(Line,'!')==1) then
               !write(*,*) 'Ignored comment: '//trim(Line)
           else
               ! Ignore unsupported lines
               !write(*,*) 'Ignored line: '//trim(Line)
           endif
       enddo
       close( UnIn )
   END SUBROUTINE ReadFlexASCII

   !> Reads a Guyan ASCII file 
   SUBROUTINE ReadGuyanASCII()
       ! Guyan reduction has 6 DOF, 0 CB DOFs
       p%nCB  = 0
       p%nTot = 6
       ! Get an available unit number for the file.
       CALL GetNewUnit( UnIn, ErrStat, ErrMsg );               if ( ErrStat /= 0 ) return
       ! Open the Primary input file.
       CALL OpenFInpFile ( UnIn, InputFile, ErrStat, ErrMsg ); if ( ErrStat /= 0 ) return

       !-------------------------- HEADER ---------------------------------------------
       CALL ReadStr( UnIn, InputFile, Line, 'Header line', 'File Header: External Platform MCKF Matrices (line 1)', ErrStat, ErrMsg)
       if ( ErrStat /= 0 ) return
       !---------------------- MASS MATRIX --------------------------------------
       CALL ReadCom( UnIn, InputFile, 'Section Header: Mass Matrix', ErrStat, ErrMsg)
       if ( ErrStat /= 0 ) return
       CALL ReadRealMatrix(UnIn, InputFile, p%PtfmAM, 'Mass Matrix', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)
       if ( ErrStat /= 0 ) return
       !---------------------- DAMPING MATRIX --------------------------------------
       CALL ReadCom( UnIn, InputFile, 'Section Header: Damping Matrix', ErrStat, ErrMsg)
       if ( ErrStat /= 0 ) return
       CALL ReadRealMatrix(UnIn, InputFile, p%Damp, 'Damping Matrix', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)
       if ( ErrStat /= 0 ) return
       !---------------------- STIFFNESS MATRIX --------------------------------------
       CALL ReadCom( UnIn, InputFile, 'Section Header: Stiffness Matrix', ErrStat, ErrMsg)
       if ( ErrStat /= 0 ) return
       CALL ReadRealMatrix(UnIn, InputFile, p%Stff, 'Stiffness Matrix', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)
       if ( ErrStat /= 0 ) return
       !---------------------- LOAD time-history --------------------------------------
       p%nPtfmFt = 0
       CALL ReadCom( UnIn, InputFile, 'Section Header: Loads time-history', ErrStat, ErrMsg)
       CALL ReadCom( UnIn, InputFile, 'Loads time-history table channel names', ErrStat, ErrMsg)
       CALL ReadCom( UnIn, InputFile, 'Loads time-history table channel units', ErrStat, ErrMsg)
       allocate(TmpAry(1:p%nTot+1))
       if (ErrStat < AbortErrLev) then
          ! let's figure out how many rows of data are in the time-history table:
          read( UnIn, *, IOSTAT=ErrStat ) TmpAry
          do while (ErrStat==0)
             p%nPtfmFt = p%nPtfmFt + 1
             read( UnIn, *, IOSTAT=ErrStat ) TmpAry
          end do
       end if
       call allocAry( p%PtfmFt,   max(1,p%nPtfmFt), p%nTot, 'p%PtfmFt',   ErrStat, ErrMsg); if ( ErrStat /= 0 ) return
       call allocAry( p%PtfmFt_t, max(1,p%nPtfmFt),         'p%PtfmFt_t', ErrStat, ErrMsg); if ( ErrStat /= 0 ) return
       if (p%nPtfmFt == 0) then
          p%PtfmFt = 0.0_ReKi
          p%PtfmFt_t = 0.0_ReKi
          p%nPtfmFt = 1
       else
          rewind(UnIn)
          do i=1,25 ! skip the first 25 rows of the file until we get to the data for the time-history table
             read(UnIn,*,IOSTAT=ErrStat) line
          end do
          do i=1,p%nPtfmFt
             call ReadAry( UnIn, InputFile, TmpAry, p%nTot+1, 'PtfmFt', 'PtfmFt time-history', ErrStat, ErrMsg)
             if ( ErrStat /= 0 ) return
             p%PtfmFt_t(i) = TmpAry(1)
             p%PtfmFt(i,:) = TmpAry(2:p%nTot+1)
          end do
       end if
       !---------------------- END OF FILE -----------------------------------------
       close( UnIn )
   END SUBROUTINE ReadGuyanASCII
END SUBROUTINE ReadReducedFile

!> This routine generates the summary file, which contains a regurgitation of  the input data and interpolated flexible body data.
SUBROUTINE ExtPtfm_PrintSum(p, OtherState, RootName, ErrStat, ErrMsg)
   ! passed variables
   TYPE(ExtPtfm_ParameterType),    INTENT(IN   )  :: p           !< Parameters of the structural dynamics module
   TYPE(ExtPtfm_OtherStateType),   INTENT(IN   )  :: OtherState  !< Other states of the structural dynamics module 
   CHARACTER(*),                   INTENT(IN   )  :: RootName    !< Root Name to write the summary file
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   ! Local variables.
   INTEGER(IntKi)               :: I                                               ! Loop counter
   INTEGER(IntKi)               :: UnSu                                            ! I/O unit number for the summary output file
!    CHARACTER(*), PARAMETER      :: Fmt1      = "(34X,3(6X,'Blade',I2,:))"          ! Format for outputting blade headings.
!    CHARACTER(*), PARAMETER      :: Fmt2      = "(34X,3(6X,A,:))"                   ! Format for outputting blade headings.
!    CHARACTER(*), PARAMETER      :: FmtDat    = '(A,T35,3(:,F13.3))'                ! Format for outputting mass and modal data.
!    CHARACTER(*), PARAMETER      :: FmtDatT   = '(A,T35,1(:,F13.8))'                ! Format for outputting time steps.
   CHARACTER(30)                :: OutPFmtS                                        ! Format to print list of selected output channel names to summary file
   CHARACTER(30)                :: OutPFmt                                         ! Format to print list of selected output channels to summary file
   CHARACTER(ChanLen),PARAMETER :: TitleStr(2) = (/ 'Parameter', 'Units    ' /)
   CHARACTER(ChanLen),PARAMETER :: TitleStrLines(2) = (/ '---------------', '---------------' /)
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Open the summary file and give it a heading.
   CALL GetNewUnit(UnSu, ErrStat, ErrMsg);
   CALL OpenFOutFile(UnSu, TRIM( RootName )//'.sum', ErrStat, ErrMsg)
   IF ( ErrStat /= ErrID_None ) RETURN
   ! Heading:
   WRITE (UnSu,'(/,A)')  '!This summary information was generated by '//TRIM( GetNVD(ExtPtfm_Ver) )// &
                         ' on '//CurDate()//' at '//CurTime()//'.'

   write(UnSu,'(A)')      '!Module input file'
   write(UnSu,'(A,A)')    'Time integration method: ',StrIntMethod(p%IntMethod)
   write(UnSu,'(A,F13.8)')'Integration time step  : ',p%EP_DeltaT
   write(UnSu,'(A)')      '!Reduction input file'
   write(UnSu,'(A,I0)')   'Number of time steps   : ',p%nPtfmFt
   write(UnSu,'(A,F13.8)')'Start time             : ',p%PtfmFt_t(1)
   write(UnSu,'(A,F13.8)')'End time               : ',p%PtfmFt_t(p%nPtfmFt)
   write(UnSu,'(A)')      '!Derived parameters'
   write(UnSu,'(A,I0)')   'Total number of DOF    : ',p%nTot
   write(UnSu,'(A,I0)')   'Number of CB modes     : ',p%nCB
! 
   write(UnSu,'(A)')'!State matrices'
   call disp2r8(UnSu, 'A',p%AMat)
   call disp2r8(UnSu, 'B',p%BMat)
   call disp2r8(UnSu, 'C',p%CMat)
   call disp2r8(UnSu, 'D',p%DMat)
   write(UnSu,'(A)')'!Input matrices'
   call disp2r8(UnSu, 'M',p%PtfmAM)
   call disp2r8(UnSu, 'K',p%Stff)
   call disp2r8(UnSu, 'C',p%Damp)
!    call disp2r8(UnSu, 'F',p%PtfmFt)
!    call disp2r8(UnSu, 'M11',p%M11)
!    call disp2r8(UnSu, 'M12',p%M12)
!    call disp2r8(UnSu, 'M21',p%M21)
!    call disp2r8(UnSu, 'M22',p%M22)
!    call disp2r8(UnSu, 'K11',p%K11)
!    call disp2r8(UnSu, 'K22',p%K22)
!    call disp2r8(UnSu, 'C11',p%C11)
!    call disp2r8(UnSu, 'C12',p%C12)
!    call disp2r8(UnSu, 'C21',p%C21)
!    call disp2r8(UnSu, 'C22',p%C22)


   OutPFmt  = '( I4, 3X,A '//TRIM(Num2LStr(ChanLen))//',1 X, A'//TRIM(Num2LStr(ChanLen))//' )'
   OutPFmtS = '( A4, 3X,A '//TRIM(Num2LStr(ChanLen))//',1 X, A'//TRIM(Num2LStr(ChanLen))//' )'
   write(UnSu,'(//,A,//)')  '!Requested Outputs:'
   write(UnSu,OutPFmtS)  "Col", TitleStr
   write(UnSu,OutPFmtS)  "---", TitleStrLines
   DO I = 0,p%NumOuts
      write (UnSu,OutPFmt)  I, p%OutParam(I)%Name, p%OutParam(I)%Units
   END DO             

   call cleanup()

CONTAINS
    !> 
    logical function Failed()
        CALL SetErrStatSimple(ErrStat, ErrMsg, 'ExtPtfm_PrintSum')
        Failed =  ErrStat >= AbortErrLev
        if(Failed) call cleanup()
    end function Failed
    !> 
    subroutine cleanup()
        if (UnSu>0) close(UnSu)
    end subroutine cleanup
END SUBROUTINE ExtPtfm_PrintSum

END MODULE ExtPtfm_MCKF_IO
