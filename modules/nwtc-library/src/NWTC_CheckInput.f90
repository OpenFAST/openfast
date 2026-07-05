!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2026 National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
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
!**********************************************************************************************************************************
!> Reusable pieces for every OpenFAST executable's `-CheckInput` mode: an in-memory collector for
!! per-component check results (built on the framework's SetErrStat/ErrStat/ErrMsg convention), a
!! console summary printer built on WrScr, and a crash-survivable incremental YAML report writer.
!!
!! Report liveness contract: CkIn_OpenReport writes a `check_status: crashed` placeholder that is never
!! rewritten (sequential formatted files cannot edit earlier lines); CkIn_CloseReport appends the
!! authoritative `overall_status:` block. A reader must treat a file with no trailing `overall_status:`
!! block as a crashed run.
MODULE NWTC_CheckInput

   USE NWTC_Library
   USE YAML, ONLY: yaml_write_comm

   IMPLICIT NONE

   PRIVATE

   INTEGER(IntKi), PARAMETER, PUBLIC :: CkIn_St_Passed      = 1  !< all checks for this component passed
   INTEGER(IntKi), PARAMETER, PUBLIC :: CkIn_St_Failed      = 2  !< one or more Severe/Fatal messages collected
   INTEGER(IntKi), PARAMETER, PUBLIC :: CkIn_St_NotUsed     = 3  !< component not enabled by the input file
   INTEGER(IntKi), PARAMETER, PUBLIC :: CkIn_St_Crashed     = 4  !< report-level liveness placeholder only
   INTEGER(IntKi), PARAMETER, PUBLIC :: CkIn_St_Unavailable = 5  !< attempted with fallback data (upstream failed)
   INTEGER(IntKi), PARAMETER, PUBLIC :: CkIn_St_Skipped     = 6  !< stage not run (blocked by upstream failures)

   INTEGER(IntKi), PARAMETER :: CkIn_NameLen = 64
   INTEGER(IntKi), PARAMETER :: CkIn_MsgLen  = 512

   !> One collected message. Severity is the scalar ErrStat CkIn_Collect was called with: SetErrStat's
   !! concatenation does not preserve per-line severity, so every line split from one call shares it.
   TYPE, PUBLIC :: CkIn_MsgType
      CHARACTER(CkIn_NameLen) :: Component = ''
      INTEGER(IntKi)          :: Severity  = ErrID_None
      CHARACTER(CkIn_NameLen) :: Source    = ''   !< inner RoutineName parsed from "RoutineName:text"
      CHARACTER(CkIn_MsgLen)  :: Text      = ''   !< message text with the "RoutineName:" prefix stripped
   END TYPE CkIn_MsgType

   !> The collector. Arrays grow on demand; default-initialize and pass INTENT(INOUT) everywhere —
   !! INTENT(OUT) would re-apply the default initializers and silently wipe collected data.
   TYPE, PUBLIC :: CheckInputCollectorType
      TYPE(CkIn_MsgType), ALLOCATABLE      :: Msgs(:)
      INTEGER(IntKi)                       :: NumMsgs      = 0
      CHARACTER(CkIn_NameLen), ALLOCATABLE :: CompNames(:) !< distinct names, in order first seen
      INTEGER(IntKi),          ALLOCATABLE :: CompStat(:)  !< aggregate CkIn_St_*, parallel to CompNames
      INTEGER(IntKi)                       :: NumComps     = 0
      INTEGER(IntKi)                       :: NumErrors    = 0
      INTEGER(IntKi)                       :: NumWarnings  = 0
      INTEGER(IntKi)                       :: UnYaml       = -1
      CHARACTER(1024)                      :: YamlFileName = ''
   END TYPE CheckInputCollectorType

   PUBLIC :: CkIn_Collect
   PUBLIC :: CkIn_ComponentStatus
   PUBLIC :: CkIn_StatusName
   PUBLIC :: CkIn_WrSummary
   PUBLIC :: CkIn_OpenReport
   PUBLIC :: CkIn_ReportComponent
   PUBLIC :: CkIn_CloseReport
   PUBLIC :: CkIn_ExitCode

CONTAINS

   !=======================================================================
   !> Records one component's check result. ErrStat/ErrMsg are exactly what a module Init returns via
   !! the SetErrStat convention: ErrMsg may be multiple "RoutineName:text" lines joined by NewLine.
   !! Status optionally forces a state that cannot be inferred from ErrStat alone
   !! ('not_used' | 'unavailable' | 'skipped' | 'failed' | 'passed'). 'crashed' is deliberately not
   !! accepted: a crashed component never returns to call this routine.
   SUBROUTINE CkIn_Collect(collector, component, ErrStat, ErrMsg, Status)

      TYPE(CheckInputCollectorType), INTENT(INOUT) :: collector
      CHARACTER(*),                  INTENT(IN)    :: component
      INTEGER(IntKi),                INTENT(IN)    :: ErrStat
      CHARACTER(*),                  INTENT(IN)    :: ErrMsg
      CHARACTER(*), OPTIONAL,        INTENT(IN)    :: Status

      CHARACTER(:), ALLOCATABLE :: Remaining, Line, Src, Txt
      INTEGER(IntKi)            :: NLPos, ColonPos, CompStat

      IF ( PRESENT(Status) ) THEN
         SELECT CASE ( TRIM(Status) )
         CASE ('not_used');    CompStat = CkIn_St_NotUsed
         CASE ('unavailable'); CompStat = CkIn_St_Unavailable
         CASE ('skipped');     CompStat = CkIn_St_Skipped
         CASE ('failed');      CompStat = CkIn_St_Failed
         CASE ('passed');      CompStat = CkIn_St_Passed
         CASE DEFAULT
            IF ( ErrStat >= ErrID_Severe ) THEN
               CompStat = CkIn_St_Failed
            ELSE
               CompStat = CkIn_St_Passed
            END IF
         END SELECT
      ELSE IF ( ErrStat >= ErrID_Severe ) THEN
         CompStat = CkIn_St_Failed
      ELSE
         CompStat = CkIn_St_Passed
      END IF

      CALL CkIn_SetComponentStatus( collector, component, CompStat )

      IF ( ErrStat /= ErrID_None .AND. LEN_TRIM(ErrMsg) > 0 ) THEN

         Remaining = ErrMsg

         DO WHILE ( LEN_TRIM(Remaining) > 0 )

            NLPos = INDEX( Remaining, NewLine )
            IF ( NLPos > 0 ) THEN
               Line      = Remaining(1:NLPos-1)
               Remaining = Remaining( (NLPos + LEN(NewLine)): )
            ELSE
               Line      = Remaining
               Remaining = ''
            END IF

            IF ( LEN_TRIM(Line) > 0 ) THEN

               ! SetErrStat writes "RoutineName:text" -- split on the FIRST colon; a line with no
               ! colon (not produced via SetErrStat) becomes text with a blank Source.
               ColonPos = INDEX( Line, ':' )
               IF ( ColonPos > 1 ) THEN
                  Src = TRIM( ADJUSTL( Line(1:ColonPos-1) ) )
                  Txt = TRIM( ADJUSTL( Line(ColonPos+1:) ) )
               ELSE
                  Src = ''
                  Txt = TRIM( ADJUSTL( Line ) )
               END IF

               CALL CkIn_AppendMsg( collector, component, ErrStat, Src, Txt )

               IF ( ErrStat >= ErrID_Severe ) THEN
                  collector%NumErrors = collector%NumErrors + 1
               ELSE IF ( ErrStat == ErrID_Warn ) THEN
                  collector%NumWarnings = collector%NumWarnings + 1
               END IF

            END IF

         END DO

      END IF

   END SUBROUTINE CkIn_Collect

   !=======================================================================
   FUNCTION CkIn_ComponentStatus(collector, component, Found) RESULT(Status)
      TYPE(CheckInputCollectorType), INTENT(IN)  :: collector
      CHARACTER(*),                  INTENT(IN)  :: component
      LOGICAL, OPTIONAL,             INTENT(OUT) :: Found
      INTEGER(IntKi)                             :: Status
      INTEGER(IntKi) :: idx
      idx = CkIn_FindComponent( collector, component )
      IF ( idx > 0 ) THEN
         Status = collector%CompStat(idx)
         IF ( PRESENT(Found) ) Found = .TRUE.
      ELSE
         Status = CkIn_St_Unavailable
         IF ( PRESENT(Found) ) Found = .FALSE.
      END IF
   END FUNCTION CkIn_ComponentStatus

   !=======================================================================
   FUNCTION CkIn_StatusName(Status) RESULT(Name)
      INTEGER(IntKi), INTENT(IN) :: Status
      CHARACTER(11)              :: Name
      SELECT CASE (Status)
      CASE (CkIn_St_Passed);      Name = 'passed'
      CASE (CkIn_St_Failed);      Name = 'failed'
      CASE (CkIn_St_NotUsed);     Name = 'not_used'
      CASE (CkIn_St_Crashed);     Name = 'crashed'
      CASE (CkIn_St_Unavailable); Name = 'unavailable'
      CASE (CkIn_St_Skipped);     Name = 'skipped'
      CASE DEFAULT;               Name = 'unavailable'
      END SELECT
   END FUNCTION CkIn_StatusName

   !=======================================================================
   !> Prints the "INPUT CHECK SUMMARY" block: per-component status, every message with severity tag,
   !! and the final PASSED/FAILED banner. Un present => that file unit; absent => console via WrScr.
   SUBROUTINE CkIn_WrSummary(collector, Un)

      TYPE(CheckInputCollectorType), INTENT(IN)           :: collector
      INTEGER(IntKi),                INTENT(IN), OPTIONAL :: Un

      INTEGER(IntKi) :: i, j
      CHARACTER(72)  :: Bar
      CHARACTER(CkIn_NameLen + CkIn_MsgLen + 96) :: Line

      Bar = REPEAT( '=', 72 )

      CALL CkIn_WrLine( Un, Bar )
      CALL CkIn_WrLine( Un, ' INPUT CHECK SUMMARY' )
      CALL CkIn_WrLine( Un, Bar )

      DO i = 1, collector%NumComps
         Line = '  '//TRIM(collector%CompNames(i))//': '//TRIM(CkIn_StatusBanner(collector%CompStat(i)))
         CALL CkIn_WrLine( Un, Line )
         DO j = 1, collector%NumMsgs
            IF ( TRIM(collector%Msgs(j)%Component) == TRIM(collector%CompNames(i)) ) THEN
               IF ( collector%Msgs(j)%Severity >= ErrID_Severe ) THEN
                  Line = '     [error] '
               ELSE IF ( collector%Msgs(j)%Severity == ErrID_Warn ) THEN
                  Line = '     [warn]  '
               ELSE
                  Line = '     [info]  '
               END IF
               IF ( LEN_TRIM(collector%Msgs(j)%Source) > 0 ) THEN
                  Line = TRIM(Line)//TRIM(collector%Msgs(j)%Source)//': '//TRIM(collector%Msgs(j)%Text)
               ELSE
                  Line = TRIM(Line)//TRIM(collector%Msgs(j)%Text)
               END IF
               CALL CkIn_WrLine( Un, Line )
            END IF
         END DO
      END DO

      CALL CkIn_WrLine( Un, Bar )
      IF ( collector%NumErrors > 0 ) THEN
         Line = ' INPUT CHECK FAILED: '//TRIM(Num2LStr(collector%NumErrors))//' errors, '// &
                TRIM(Num2LStr(collector%NumWarnings))//' warnings'
      ELSE
         Line = ' INPUT CHECK PASSED'
         IF ( collector%NumWarnings > 0 ) THEN
            Line = TRIM(Line)//' ('//TRIM(Num2LStr(collector%NumWarnings))//' warnings)'
         END IF
      END IF
      CALL CkIn_WrLine( Un, Line )
      CALL CkIn_WrLine( Un, Bar )

   END SUBROUTINE CkIn_WrSummary

   !=======================================================================
   !> Opens <RootName>.verify.yaml and writes the header + CRASHED liveness placeholder (never
   !! rewritten -- see module header contract).
   SUBROUTINE CkIn_OpenReport(collector, RootName, ErrStat, ErrMsg)

      TYPE(CheckInputCollectorType), INTENT(INOUT) :: collector
      CHARACTER(*),                  INTENT(IN)    :: RootName
      INTEGER(IntKi),                INTENT(OUT)   :: ErrStat
      CHARACTER(*),                  INTENT(OUT)   :: ErrMsg

      CHARACTER(*), PARAMETER :: RoutineName = 'CkIn_OpenReport'
      INTEGER(IntKi)          :: ErrStat2
      CHARACTER(ErrMsgLen)    :: ErrMsg2

      ErrStat = ErrID_None
      ErrMsg  = ''

      collector%YamlFileName = TRIM(RootName)//'.verify.yaml'

      CALL GetNewUnit( collector%UnYaml, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      CALL OpenFOutFile( collector%UnYaml, TRIM(collector%YamlFileName), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      CALL yaml_write_comm( collector%UnYaml, &
         'OpenFAST -CheckInput report, auto-generated '//CurDate()//' '//CurTime()//' -- do not hand-edit', &
         ErrStat2, ErrMsg2 )

      WRITE (collector%UnYaml, '(A)') 'check_status: crashed  # liveness placeholder; see "overall_status" at EOF'
      WRITE (collector%UnYaml, '(A)') 'components:'

      FLUSH(collector%UnYaml)

   END SUBROUTINE CkIn_OpenReport

   !=======================================================================
   !> Appends one component's status + messages to the open report, then FLUSHes so a segfault in the
   !! NEXT component's init still leaves this block intact on disk. Call once per component,
   !! immediately after its CkIn_Collect call.
   SUBROUTINE CkIn_ReportComponent(collector, component, ErrStat, ErrMsg)

      TYPE(CheckInputCollectorType), INTENT(IN)  :: collector
      CHARACTER(*),                  INTENT(IN)  :: component
      INTEGER(IntKi),                INTENT(OUT) :: ErrStat
      CHARACTER(*),                  INTENT(OUT) :: ErrMsg

      CHARACTER(*), PARAMETER :: RoutineName = 'CkIn_ReportComponent'
      INTEGER(IntKi) :: Un, Stat, i, nErr, nWarn, nMsg, IOS

      ErrStat = ErrID_None
      ErrMsg  = ''
      Un      = collector%UnYaml

      IF ( Un <= 0 ) THEN
         CALL SetErrStat( ErrID_Severe, 'Report file is not open; call CkIn_OpenReport first.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      Stat  = CkIn_ComponentStatus( collector, component )
      nErr  = 0
      nWarn = 0
      nMsg  = 0
      DO i = 1, collector%NumMsgs
         IF ( TRIM(collector%Msgs(i)%Component) == TRIM(component) ) THEN
            nMsg = nMsg + 1
            IF ( collector%Msgs(i)%Severity >= ErrID_Severe ) nErr  = nErr  + 1
            IF ( collector%Msgs(i)%Severity == ErrID_Warn   ) nWarn = nWarn + 1
         END IF
      END DO

      WRITE (Un, '(2X,"- name: ",A)', IOSTAT=IOS) TRIM(component)
      WRITE (Un, '(4X,"status: ",A)')             TRIM(CkIn_StatusName(Stat))
      WRITE (Un, '(4X,"errors: ",I0)')            nErr
      WRITE (Un, '(4X,"warnings: ",I0)')          nWarn

      IF ( nMsg == 0 ) THEN
         WRITE (Un, '(4X,"messages: []")')
      ELSE
         WRITE (Un, '(4X,"messages:")')
         DO i = 1, collector%NumMsgs
            IF ( TRIM(collector%Msgs(i)%Component) == TRIM(component) ) THEN
               WRITE (Un, '(6X,"- severity: ",A)') TRIM(CkIn_SeverityName(collector%Msgs(i)%Severity))
               IF ( LEN_TRIM(collector%Msgs(i)%Source) > 0 ) THEN
                  WRITE (Un, '(8X,"source: ",A)') TRIM(collector%Msgs(i)%Source)
               END IF
               WRITE (Un, '(8X,"text: """,A,"""")') CkIn_YamlEscape( TRIM(collector%Msgs(i)%Text) )
            END IF
         END DO
      END IF

      IF ( IOS /= 0 ) THEN
         CALL SetErrStat( ErrID_Severe, 'Error writing component "'//TRIM(component)//'" to '// &
                          TRIM(collector%YamlFileName)//'.', ErrStat, ErrMsg, RoutineName )
      END IF

      FLUSH(Un)

   END SUBROUTINE CkIn_ReportComponent

   !=======================================================================
   !> Writes the authoritative overall_status/counts block and closes the report.
   SUBROUTINE CkIn_CloseReport(collector, ErrStat, ErrMsg)

      TYPE(CheckInputCollectorType), INTENT(INOUT) :: collector
      INTEGER(IntKi),                INTENT(OUT)   :: ErrStat
      CHARACTER(*),                  INTENT(OUT)   :: ErrMsg

      CHARACTER(*), PARAMETER :: RoutineName = 'CkIn_CloseReport'
      INTEGER(IntKi) :: Un, Overall, IOS

      ErrStat = ErrID_None
      ErrMsg  = ''
      Un      = collector%UnYaml

      IF ( Un <= 0 ) THEN
         CALL SetErrStat( ErrID_Severe, 'Report file is not open; call CkIn_OpenReport first.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

      IF ( collector%NumErrors > 0 ) THEN
         Overall = CkIn_St_Failed
      ELSE
         Overall = CkIn_St_Passed
      END IF

      WRITE (Un, '(A)', IOSTAT=IOS) 'overall_status: '//TRIM(CkIn_StatusName(Overall))
      WRITE (Un, '(A,I0)')          'total_errors: ',       collector%NumErrors
      WRITE (Un, '(A,I0)')          'total_warnings: ',     collector%NumWarnings
      WRITE (Un, '(A,I0)')          'components_checked: ', collector%NumComps

      FLUSH(Un)
      CLOSE(Un)

      IF ( IOS /= 0 ) THEN
         CALL SetErrStat( ErrID_Severe, 'Error finalizing '//TRIM(collector%YamlFileName)//'.', ErrStat, ErrMsg, RoutineName )
      END IF

      collector%UnYaml = -1

   END SUBROUTINE CkIn_CloseReport

   !=======================================================================
   !> 0 if no Severe/Fatal message and no failed component; 1 otherwise. Warnings never affect it.
   FUNCTION CkIn_ExitCode(collector) RESULT(Code)
      TYPE(CheckInputCollectorType), INTENT(IN) :: collector
      INTEGER(IntKi)                            :: Code
      INTEGER(IntKi)                            :: i
      Code = 0
      IF ( collector%NumErrors > 0 ) THEN
         Code = 1
         RETURN
      END IF
      DO i = 1, collector%NumComps
         IF ( collector%CompStat(i) == CkIn_St_Failed ) THEN
            Code = 1
            RETURN
         END IF
      END DO
   END FUNCTION CkIn_ExitCode

   !=======================================================================
   ! ---- private helpers ----

   SUBROUTINE CkIn_WrLine(Un, Str)
      INTEGER(IntKi), INTENT(IN), OPTIONAL :: Un
      CHARACTER(*),   INTENT(IN)           :: Str
      IF ( PRESENT(Un) ) THEN
         WRITE (Un, '(A)') TRIM(Str)
      ELSE
         CALL WrScr( TRIM(Str) )
      END IF
   END SUBROUTINE CkIn_WrLine

   FUNCTION CkIn_StatusBanner(Status) RESULT(Txt)
      INTEGER(IntKi), INTENT(IN) :: Status
      CHARACTER(11)              :: Txt
      SELECT CASE (Status)
      CASE (CkIn_St_Passed);      Txt = 'PASSED'
      CASE (CkIn_St_Failed);      Txt = 'FAILED'
      CASE (CkIn_St_NotUsed);     Txt = 'NOT USED'
      CASE (CkIn_St_Crashed);     Txt = 'CRASHED'
      CASE (CkIn_St_Unavailable); Txt = 'UNAVAILABLE'
      CASE (CkIn_St_Skipped);     Txt = 'SKIPPED'
      CASE DEFAULT;               Txt = 'UNKNOWN'
      END SELECT
   END FUNCTION CkIn_StatusBanner

   FUNCTION CkIn_SeverityName(Severity) RESULT(Name)
      INTEGER(IntKi), INTENT(IN) :: Severity
      CHARACTER(5)               :: Name
      SELECT CASE (Severity)
      CASE (ErrID_Info);   Name = 'info'
      CASE (ErrID_Warn);   Name = 'warn'
      CASE (ErrID_Severe); Name = 'error'
      CASE (ErrID_Fatal);  Name = 'fatal'
      CASE DEFAULT;        Name = 'none'
      END SELECT
   END FUNCTION CkIn_SeverityName

   !> Minimal YAML double-quoted-scalar escaping (backslash, then double-quote). Error text routinely
   !! contains ':' and sometimes '"'; YAML.f90's writers have no escaping and no list-of-mappings
   !! support, so messages are hand-written quoted scalars.
   FUNCTION CkIn_YamlEscape(Str) RESULT(Esc)
      CHARACTER(*), INTENT(IN)  :: Str
      CHARACTER(:), ALLOCATABLE :: Esc
      INTEGER(IntKi) :: i
      Esc = ''
      DO i = 1, LEN_TRIM(Str)
         SELECT CASE ( Str(i:i) )
         CASE ('\')
            Esc = Esc//'\\'
         CASE ('"')
            Esc = Esc//'\"'
         CASE DEFAULT
            Esc = Esc//Str(i:i)
         END SELECT
      END DO
   END FUNCTION CkIn_YamlEscape

   FUNCTION CkIn_FindComponent(collector, component) RESULT(Idx)
      TYPE(CheckInputCollectorType), INTENT(IN) :: collector
      CHARACTER(*),                  INTENT(IN) :: component
      INTEGER(IntKi)                            :: Idx
      INTEGER(IntKi) :: i
      Idx = 0
      DO i = 1, collector%NumComps
         IF ( TRIM(collector%CompNames(i)) == TRIM(component) ) THEN
            Idx = i
            RETURN
         END IF
      END DO
   END FUNCTION CkIn_FindComponent

   SUBROUTINE CkIn_SetComponentStatus(collector, component, Status)
      TYPE(CheckInputCollectorType), INTENT(INOUT) :: collector
      CHARACTER(*),                  INTENT(IN)    :: component
      INTEGER(IntKi),                INTENT(IN)    :: Status
      INTEGER(IntKi) :: idx
      idx = CkIn_FindComponent( collector, component )
      IF ( idx > 0 ) THEN
         ! Failed status is sticky: a later benign collect for the same component (e.g. an info-level
         ! note after a fatal) must not downgrade the aggregate back to passed.
         IF ( collector%CompStat(idx) == CkIn_St_Failed .AND. Status == CkIn_St_Passed ) RETURN
         collector%CompStat(idx) = Status
      ELSE
         CALL CkIn_GrowComps( collector )
         collector%NumComps                      = collector%NumComps + 1
         collector%CompNames(collector%NumComps) = component
         collector%CompStat(collector%NumComps)  = Status
      END IF
   END SUBROUTINE CkIn_SetComponentStatus

   SUBROUTINE CkIn_AppendMsg(collector, component, Severity, Source, Text)
      TYPE(CheckInputCollectorType), INTENT(INOUT) :: collector
      CHARACTER(*),                  INTENT(IN)    :: component
      INTEGER(IntKi),                INTENT(IN)    :: Severity
      CHARACTER(*),                  INTENT(IN)    :: Source
      CHARACTER(*),                  INTENT(IN)    :: Text
      CALL CkIn_GrowMsgs( collector )
      collector%NumMsgs = collector%NumMsgs + 1
      collector%Msgs(collector%NumMsgs)%Component = component
      collector%Msgs(collector%NumMsgs)%Severity  = Severity
      collector%Msgs(collector%NumMsgs)%Source    = Source
      collector%Msgs(collector%NumMsgs)%Text      = Text
   END SUBROUTINE CkIn_AppendMsg

   !> ALLOCATE+MOVE_ALLOC doubling growth (portable across the older compilers OpenFAST supports;
   !! whole-array realloc-on-assignment for derived-type arrays is less universally reliable).
   SUBROUTINE CkIn_GrowComps(collector)
      TYPE(CheckInputCollectorType), INTENT(INOUT) :: collector
      CHARACTER(CkIn_NameLen), ALLOCATABLE :: TmpNames(:)
      INTEGER(IntKi),          ALLOCATABLE :: TmpStat(:)
      INTEGER(IntKi) :: OldCap, NewCap
      IF ( .NOT. ALLOCATED(collector%CompNames) ) THEN
         ALLOCATE( collector%CompNames(8) )
         ALLOCATE( collector%CompStat(8) )
         RETURN
      END IF
      OldCap = SIZE(collector%CompNames)
      IF ( collector%NumComps < OldCap ) RETURN
      NewCap = OldCap * 2
      ALLOCATE( TmpNames(NewCap) ); TmpNames(1:OldCap) = collector%CompNames
      ALLOCATE( TmpStat(NewCap)  ); TmpStat(1:OldCap)  = collector%CompStat
      CALL MOVE_ALLOC( TmpNames, collector%CompNames )
      CALL MOVE_ALLOC( TmpStat,  collector%CompStat )
   END SUBROUTINE CkIn_GrowComps

   SUBROUTINE CkIn_GrowMsgs(collector)
      TYPE(CheckInputCollectorType), INTENT(INOUT) :: collector
      TYPE(CkIn_MsgType), ALLOCATABLE :: TmpMsgs(:)
      INTEGER(IntKi) :: OldCap, NewCap
      IF ( .NOT. ALLOCATED(collector%Msgs) ) THEN
         ALLOCATE( collector%Msgs(16) )
         RETURN
      END IF
      OldCap = SIZE(collector%Msgs)
      IF ( collector%NumMsgs < OldCap ) RETURN
      NewCap = OldCap * 2
      ALLOCATE( TmpMsgs(NewCap) ); TmpMsgs(1:OldCap) = collector%Msgs
      CALL MOVE_ALLOC( TmpMsgs, collector%Msgs )
   END SUBROUTINE CkIn_GrowMsgs

END MODULE NWTC_CheckInput
