!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
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
!
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE ModMesh

!======================================================================================================================
! WARNING:  Because this code uses some preprocessor directives to comment out code
!           and because some lines go beyond column 132, you must specify the following compiler options:
!
!              Intel:   /fpp
!              Gnu:     -x f95-cpp-input -ffree-line-length-none
!======================================================================================================================

   USE ModMesh_Types
   IMPLICIT NONE

   INTEGER, PARAMETER :: BUMPUP = 64  ! do not set to less than 2

   INTERFACE MeshConstructElement
      MODULE PROCEDURE MeshConstructElement_1PT ,                            &
                       MeshConstructElement_2PT , MeshConstructElement_3PT , &
                       MeshConstructElement_4PT ,                            &
                       MeshConstructElement_6PT ,                            &
                       MeshConstructElement_8PT ,                            &
                       MeshConstructElement_10PT,                            &
                                                  MeshConstructElement_15PT, &
                       MeshConstructElement_20PT
   END INTERFACE

CONTAINS

   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE MeshWrBin ( UnIn, M, ErrStat, ErrMsg, FileName)
   !...............................................................................................................................
   ! This routine writes mesh information in binary form. If UnIn is < 0, it gets a new unit number and opens the file,
   ! otherwise the file is appended. It is up to the caller of this routine to close the file when it's finished.
   !...............................................................................................................................
      INTEGER, INTENT(INOUT)                ::  UnIn     ! fortran output unit
      TYPE(MeshType),INTENT(INOUT)          ::      M    ! mesh to be reported on

      INTEGER(IntKi),  INTENT(OUT)          :: ErrStat   ! Indicates whether an error occurred (see NWTC_Library)
      CHARACTER(*),    INTENT(OUT)          :: ErrMsg    ! Error message associated with the ErrStat
      CHARACTER(*),    INTENT(IN), OPTIONAL :: FileName  ! Name of the file to write the output in

      ! local variables
      INTEGER(IntKi)                        :: ErrStat2  ! Temporary storage for local errors
      INTEGER(IntKi)                        :: I         ! loop counter

      IF (UnIn < 0) THEN
         CALL GetNewUnit( UnIn, ErrStat, ErrMsg )

         CALL OpenBOutFile ( UnIn, TRIM(FileName), ErrStat, ErrMsg )
         IF ( ErrStat >= AbortErrLev ) RETURN
      END IF


      ! Write information about mesh structure:
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(ReKi,B4Ki)
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(FIELDMASK_SIZE,B4Ki)
      WRITE (UnIn, IOSTAT=ErrStat2)   M%fieldmask           ! BJJ: do we need to verify that this is size B4Ki?
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(M%Nnodes,B4Ki)
      WRITE (UnIn, IOSTAT=ErrStat2)   INT(M%nelemlist,B4Ki)


      !...........
      ! Write nodal information:
      !...........

      WRITE (UnIn, IOSTAT=ErrStat2)   M%Position
         IF ( ErrStat2 /= 0 ) THEN
            CALL CheckError( ErrID_Fatal, 'Error writing Position to the mesh binary file.' )
            RETURN
         END IF

      WRITE (UnIn, IOSTAT=ErrStat2)   M%RefOrientation
         IF ( ErrStat2 /= 0 ) THEN
            CALL CheckError( ErrID_Fatal, 'Error writing RefOrientation to the mesh binary file.' )
            RETURN
         END IF


      ! Write fields:

      IF ( M%fieldmask(MASKID_FORCE) .AND. ALLOCATED(M%Force) ) THEN
         WRITE (UnIn, IOSTAT=ErrStat2)   M%Force
            IF ( ErrStat2 /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, 'Error writing Force to the mesh binary file.' )
               RETURN
            END IF
      END IF

      IF ( M%fieldmask(MASKID_Moment) .AND. ALLOCATED(M%Moment)) THEN
         WRITE (UnIn, IOSTAT=ErrStat2)   M%Moment
            IF ( ErrStat2 /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, 'Error writing Force to the mesh binary file.' )
               RETURN
            END IF
      END IF

      IF ( M%fieldmask(MASKID_ORIENTATION) .AND. ALLOCATED(M%Orientation)) THEN
         WRITE (UnIn, IOSTAT=ErrStat2)   M%Orientation
            IF ( ErrStat2 /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, 'Error writing Orientation to the mesh binary file.' )
               RETURN
            END IF
      END IF

      IF ( M%fieldmask(MASKID_TRANSLATIONDISP) .AND. ALLOCATED(M%TranslationDisp)) THEN
         WRITE (UnIn, IOSTAT=ErrStat2)   M%TranslationDisp
            IF ( ErrStat2 /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, 'Error writing TranslationDisp to the mesh binary file.' )
               RETURN
            END IF
      END IF

      IF ( M%fieldmask(MASKID_TRANSLATIONVEL) .AND. ALLOCATED(M%TranslationVel)) THEN
         WRITE (UnIn, IOSTAT=ErrStat2)   M%TranslationVel
            IF ( ErrStat2 /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, 'Error writing TranslationVel to the mesh binary file.' )
               RETURN
            END IF
      END IF

      IF ( M%fieldmask(MASKID_ROTATIONVEL) .AND. ALLOCATED(M%RotationVel)) THEN
         WRITE (UnIn, IOSTAT=ErrStat2)   M%RotationVel
            IF ( ErrStat2 /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, 'Error writing RotationVel to the mesh binary file.' )
               RETURN
            END IF
      END IF

      IF ( M%fieldmask(MASKID_TRANSLATIONACC) .AND. ALLOCATED(M%TranslationAcc)) THEN
         WRITE (UnIn, IOSTAT=ErrStat2)   M%TranslationAcc
            IF ( ErrStat2 /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, 'Error writing TranslationAcc to the mesh binary file.' )
               RETURN
            END IF
      END IF

      IF ( M%fieldmask(MASKID_ROTATIONACC) .AND. ALLOCATED(M%RotationAcc)) THEN
         WRITE (UnIn, IOSTAT=ErrStat2)   M%RotationAcc
            IF ( ErrStat2 /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, 'Error writing RotationAcc to the mesh binary file.' )
               RETURN
            END IF
      END IF

      IF ( M%fieldmask(MASKID_SCALAR) .AND. ALLOCATED(M%Scalars)) THEN
         WRITE (UnIn, IOSTAT=ErrStat2)   M%Scalars
            IF ( ErrStat2 /= 0 ) THEN
               CALL CheckError( ErrID_Fatal, 'Error writing Scalars to the mesh binary file.' )
               RETURN
            END IF
      END IF


      !...........
      ! Write element information:
      !...........

      DO i=1,M%nelemlist
         WRITE (UnIn, IOSTAT=ErrStat2)   INT(M%ElemList(i)%Element%Xelement       ,B4Ki) ! what kind of element
         WRITE (UnIn, IOSTAT=ErrStat2)   INT(SIZE(M%ElemList(i)%Element%ElemNodes),B4Ki) ! how many nodes
         WRITE (UnIn, IOSTAT=ErrStat2)   INT(M%ElemList(i)%Element%ElemNodes      ,B4Ki) ! which nodes
      END DO

   CONTAINS
      !............................................................................................................................
      SUBROUTINE CheckError(ErrID,Msg)
      !............................................................................................................................
      ! This subroutine sets the error message and level
      !............................................................................................................................

            ! Passed arguments
         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


         !.........................................................................................................................
         ! Set error status/message;
         !.........................................................................................................................

         IF ( ErrID /= ErrID_None ) THEN

            IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
            ErrMsg = TRIM(ErrMsg)//' '//TRIM(Msg)
            ErrStat = MAX(ErrStat, ErrID)

            !......................................................................................................................
            ! Clean up if we're going to return on error: close file, deallocate local arrays
            !......................................................................................................................
            IF ( ErrStat >= AbortErrLev ) THEN
            END IF

         END IF

      END SUBROUTINE CheckError
      !............................................................................................................................
   END SUBROUTINE MeshWrBin
   !-------------------------------------------------------------------------------------------------------------------------------


   SUBROUTINE MeshPrintInfo ( U, M, N)
     INTEGER, INTENT(IN   )                ::      U  ! fortran output unit
     TYPE(MeshType),INTENT(INOUT)          ::      M  ! mesh to be reported on
     INTEGER, OPTIONAL,INTENT(IN   )       ::      N  ! Number to print, default 5
    ! Local
     INTEGER isz,i,j,nn,CtrlCode,Ielement,Xelement
     INTEGER                    :: ErrStat
     CHARACTER(256)             :: ErrMess

     nn = 5
     IF (PRESENT(N)) nn = N

     write(U,*)'-----------  MeshPrintInfo:  -------------'

     write(U,*)  'Initialized: ', M%initialized
     write(U,*)  'Committed:   ', M%Committed
     IF ( ASSOCIATED(M%RemapFlag) ) write(U,*)  'Remap Flag: ', M%RemapFlag

     write(U,*)  'Fieldmask:   ', M%FieldMask
     IF ( M%FieldMask( MASKID_FORCE           ) )  write(U,*)  '  Defined : Force'
     IF ( M%FieldMask( MASKID_MOMENT          ) )  write(U,*)  '  Defined : Moment'
     IF ( M%FieldMask( MASKID_ORIENTATION     ) )  write(U,*)  '  Defined : Orientation'
     IF ( M%FieldMask( MASKID_TRANSLATIONDISP ) )  write(U,*)  '  Defined : TranslationDisp'
     IF ( M%FieldMask( MASKID_TRANSLATIONVEL  ) )  write(U,*)  '  Defined : TranslationVel'
     IF ( M%FieldMask( MASKID_ROTATIONVEL     ) )  write(U,*)  '  Defined : RotationVel'
     IF ( M%FieldMask( MASKID_TRANSLATIONACC  ) )  write(U,*)  '  Defined : TranslationAcc'
     IF ( M%FieldMask( MASKID_ROTATIONACC     ) )  write(U,*)  '  Defined : RotationAcc'
     IF ( M%FieldMask( MASKID_SCALAR          ) )  write(U,*)  '  Defined : Scalar'
     write(U,*)  'Ios:         ', M%Ios
     write(U,*)  'Nnodes:      ', M%Nnodes


     DO i = 1, NELEMKINDS
       IF ( M%ElemTable(i)%nelem .GT. 0 ) THEN
         WRITE(U,*)ElemNames(i),' nelem: ',M%ElemTable(i)%nelem,' max: ',M%ElemTable(i)%maxelem
         IF(M%initialized.AND.ASSOCIATED(M%ElemTable(i)%Elements))THEN
           DO j = 1,min(nn,M%ElemTable(i)%nelem)
             write(U,*)' ',j,M%ElemTable(i)%Elements(j)%ElemNodes(:)
           ENDDO
         ENDIF
       ENDIF
     END DO

! Here are some built in derived data types that can represent values at the nodes
! the last dimension of each of these has range 1:nnodes for the mesh being represented
! and they are indexed by the element arrays above
! only some of these would be allocated, depending on what's being represented
! on the mesh.
! Whether or not these are allocated is indicted in the fieldmask, which can
! be interrogated by a routine using an instance of the type. If you add a field
! here, be sure to change the table of parameters used to size and index fieldmask above.
     IF(M%initialized.AND.ASSOCIATED(M%Position))THEN
       isz=size(M%Position,2)
       write(U,*)'Position: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%Position(:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%RefOrientation))THEN
       isz=size(M%RefOrientation,3)
       write(U,*)'RefOrientation: ',isz
       DO i=1,min(nn,isz) !bjj: printing this like a matrix:
         write(U,'(1X,I3, 3(1X,F10.4))') i, M%RefOrientation(1,:,i)
         write(U,'(4X,    3(1X,F10.4))')    M%RefOrientation(2,:,i)
         write(U,'(4X,    3(1X,F10.4))')    M%RefOrientation(3,:,i)
       ENDDO
     ENDIF


     IF(ALLOCATED(M%Force))THEN
       isz=size(M%Force,2)
       write(U,*)'Force: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%Force(:,i)
       ENDDO
     ENDIF
     IF(ALLOCATED(M%Moment))THEN
       isz=size(M%Moment,2)
       write(U,*)'Moment: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%Moment(:,i)
       ENDDO
     ENDIF
     IF(ALLOCATED(M%Orientation))THEN
       isz=size(M%Orientation,3)
       write(U,*)'Orientation: ',isz,' node(s)'
       DO i=1,min(nn,isz) !bjj: printing this like a matrix:
         write(U,'(1X,I3, 3(1X,F10.4))') i, M%Orientation(1,:,i)
         write(U,'(4X,    3(1X,F10.4))')    M%Orientation(2,:,i)
         write(U,'(4X,    3(1X,F10.4))')    M%Orientation(3,:,i)
       ENDDO
     ENDIF
     IF(ALLOCATED(M%TranslationDisp))THEN
       isz=size(M%TranslationDisp,2)
       write(U,*)'TranslationDisp: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%TranslationDisp(:,i)
       ENDDO
     ENDIF
     IF(ALLOCATED(M%RotationVel))THEN
       isz=size(M%RotationVel,2)
       write(U,*)'RotationVel: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%RotationVel(:,i)
       ENDDO
     ENDIF
     IF(ALLOCATED(M%TranslationVel))THEN
       isz=size(M%TranslationVel,2)
       write(U,*)'TranslationVel: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%TranslationVel(:,i)
       ENDDO
     ENDIF
     IF(ALLOCATED(M%RotationAcc))THEN
       isz=size(M%RotationAcc,2)
       write(U,*)'RotationAcc: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%RotationAcc(:,i)
       ENDDO
     ENDIF
     IF(ALLOCATED(M%TranslationAcc))THEN
       isz=size(M%TranslationAcc,2)
       write(U,*)'TranslationAcc: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%TranslationAcc(:,i)
       ENDDO
     ENDIF
     IF(ALLOCATED(M%Scalars))THEN
       isz=size(M%Scalars,1)
       write(U,*)'Scalars: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%Scalars(:,i)
       ENDDO
     ENDIF
     write(U,*)'--------- Traverse Element List ----------'
     CtrlCode = 0
     CALL MeshNextElement( M, CtrlCode, ErrStat, ErrMess, Ielement=Ielement, Xelement=Xelement )
     IF (ErrStat >= AbortErrLev) THEN
        WRITE(U,*) ' Error in MeshNextElement(): '
        WRITE(U,*) TRIM(ErrMess)
        CtrlCode = MESH_NOMORE
      END IF

     DO WHILE ( CtrlCode .NE. MESH_NOMORE )
       WRITE(U,*)'  Ielement: ', Ielement,' ',ElemNames(Xelement), &
                 ' det_jac: ',M%ElemList(Ielement)%Element%det_jac,  &
                 ' Nodes: ',M%ElemList(Ielement)%Element%ElemNodes

       CtrlCode = MESH_NEXT
       CALL MeshNextElement( M, CtrlCode, ErrStat, ErrMess, Ielement=Ielement, Xelement=Xelement )
        IF (ErrStat >= AbortErrLev) THEN
           WRITE(U,*) ' Error in MeshNextElement(): '
           WRITE(U,*) TRIM(ErrMess)
           RETURN
         END IF

     ENDDO
     write(U,*)'---------  End of Element List  ----------'

   END SUBROUTINE MeshPrintInfo

   ! operations to create a mesh

   SUBROUTINE MeshCreate ( BlankMesh                                                       &
                          ,IOS                                                             &
                          ,Nnodes                                                          &
                          ,ErrStat                                                         &
                          ,ErrMess                                                         &
                             ! optional arguments that say whether to allocate fields      &
                             ! in the mesh. These are always dimensioned npoints           &
                          ,Force                                                           &
                          ,Moment                                                          &
                          ,Orientation                                                     &
                          ,TranslationDisp                                                 &
                          ,TranslationVel                                                  &
                          ,RotationVel                                                     &
                          ,TranslationAcc                                                  &
                          ,RotationAcc                                                     &
                          ,nScalars                                                        &
                          ,IsNewSibling                                                    &
                         )

      TYPE(MeshType), INTENT(INOUT)   :: BlankMesh ! Mesh to be created
      INTEGER,INTENT(IN)         :: IOS                  ! input (COMPONENT_INPUT), output(COMPONENT_OUTPUT), or state(COMPONENT_STATE)
      INTEGER,INTENT(IN)         :: Nnodes               ! Number of nodes in mesh
      INTEGER(IntKi),INTENT(OUT) :: ErrStat              ! error status/level
      CHARACTER(*),INTENT(OUT)   :: ErrMess              ! error message
                                   ! optional arguments from here down
                                   ! optional arguments that say whether to allocate fields
                                   ! in the mesh. These are always dimensioned npoints
      LOGICAL,OPTIONAL,INTENT(IN):: Force                ! If present and true, allocate Force field
      LOGICAL,OPTIONAL,INTENT(IN):: Moment               ! If present and true, allocate Moment field
      LOGICAL,OPTIONAL,INTENT(IN):: Orientation          ! If present and true, allocate Orientation field
      LOGICAL,OPTIONAL,INTENT(IN):: TranslationDisp      ! If present and true, allocate TranslationDisp field
      LOGICAL,OPTIONAL,INTENT(IN):: TranslationVel       ! If present and true, allocate TranslationVel field
      LOGICAL,OPTIONAL,INTENT(IN):: RotationVel          ! If present and true, allocate RotationVel field
      LOGICAL,OPTIONAL,INTENT(IN):: TranslationAcc       ! If present and true, allocate TranslationAcc field
      LOGICAL,OPTIONAL,INTENT(IN):: RotationAcc          ! If present and true, allocate RotationAcc field
!
      INTEGER,OPTIONAL,INTENT(IN):: nScalars             ! If present and > 0, allocate nScalars Scalars
      LOGICAL,OPTIONAL,INTENT(IN):: IsNewSibling         ! If present and true, this is an new sibling so don't allocate new shared fields (RemapFlag, position, RefOrientation, and ElemTable)

    ! Local
      INTEGER i
      LOGICAL                    :: IsNewSib

      LOGICAL                    :: IsMotion
      LOGICAL                    :: IsLoad


         ! Local initializations:

      ErrStat = ErrID_None
      ErrMess = ""
      IsMotion = .FALSE.
      IsLoad   = .FALSE.

      IF ( mesh_debug ) print*,'Called MeshCreate'

      CALL MeshDestroy( BlankMesh, ErrStat, ErrMess, .TRUE. )
                                                        ! make sure we're not leaving any pointers dangling
                                                        ! and nullify them for good measure
                                                        ! See comment on optional IgnoreSibling argument
                                                        ! in definition of MeshDestroy
      IF (ErrStat >= AbortErrLev) RETURN

      BlankMesh%initialized = .TRUE.
      BlankMesh%IOS         = IOS

!bjj: check that IOS is valid and Nnodes > 0?

      BlankMesh%Nnodes = Nnodes
      BlankMesh%nelemlist = 0 ; BlankMesh%maxelemlist = 0 ;


      ! These fields are shared between siblings, so we don't want to recreate space for them here.
      IsNewSib = .FALSE.
      IF ( PRESENT(IsNewSibling) ) IsNewSib = IsNewSibling

      IF ( .NOT. IsNewSib ) THEN
         CALL AllocPAry( BlankMesh%Position, 3, Nnodes, 'MeshCreate: Position' )
         CALL AllocPAry( BlankMesh%RefOrientation, 3, 3, Nnodes, 'MeshCreate: RefOrientation' )

         ALLOCATE(BlankMesh%ElemTable(NELEMKINDS))
         DO i = 1, NELEMKINDS
            BlankMesh%ElemTable(i)%nelem = 0  ; BlankMesh%ElemTable(i)%maxelem = 0
            NULLIFY(BlankMesh%ElemTable(i)%Elements )
         ENDDO

         ALLOCATE(BlankMesh%RemapFlag, Stat=ErrStat ) ! assign some space for this pointer to point to
         BlankMesh%RemapFlag = .true.

      ELSE
         NULLIFY( BlankMesh%Position )
         NULLIFY( BlankMesh%RefOrientation )
         NULLIFY( BlankMesh%ElemTable )
         NULLIFY( BlankMesh%ElemList )
         NULLIFY( BlankMesh%RemapFlag )
      END IF
      NULLIFY( BlankMesh%SiblingMesh )

   ! handle optionals
      BlankMesh%FieldMask = .FALSE.

      IF ( PRESENT(Force) ) THEN
         IF ( Force ) THEN
            CALL AllocAry( BlankMesh%Force, 3, Nnodes, 'MeshCreate: Force', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%Force = 0.
            BlankMesh%FieldMask(MASKID_FORCE) = .TRUE.
            IsLoad = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(Moment) ) THEN
         IF ( Moment ) THEN
            CALL AllocAry( BlankMesh%Moment, 3, Nnodes, 'MeshCreate: Moment', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%Moment = 0.
            BlankMesh%FieldMask(MASKID_MOMENT) = .TRUE.
            IsLoad = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(Orientation) ) THEN
         IF ( Orientation ) THEN
            CALL AllocAry( BlankMesh%Orientation, 3, 3, Nnodes, 'MeshCreate: Orientation', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            CALL Eye(BlankMesh%Orientation, ErrStat, ErrMess)  ! set this orientation to the identity matrix
            BlankMesh%FieldMask(MASKID_ORIENTATION) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(TranslationDisp) ) THEN
         IF ( TranslationDisp ) THEN
            CALL AllocAry( BlankMesh%TranslationDisp, 3, Nnodes, 'MeshCreate: TranslationDisp', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationDisp = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONDISP) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(TranslationVel) ) THEN
         IF ( TranslationVel ) THEN
            CALL AllocAry( BlankMesh%TranslationVel, 3, Nnodes, 'MeshCreate: TranslationVel', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationVel = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONVEL) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(RotationVel) ) THEN
         IF ( RotationVel ) THEN
            CALL AllocAry( BlankMesh%RotationVel, 3, Nnodes, 'MeshCreate: RotationVel', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%RotationVel = 0.
            BlankMesh%FieldMask(MASKID_ROTATIONVEL) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(TranslationAcc) ) THEN
         IF ( TranslationAcc ) THEN
            CALL AllocAry( BlankMesh%TranslationAcc, 3, Nnodes, 'MeshCreate: TranslationAcc', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationAcc = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONACC) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(RotationAcc) ) THEN
         IF ( RotationAcc ) THEN
            CALL AllocAry( BlankMesh%RotationAcc, 3, Nnodes, 'MeshCreate: RotationAcc', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%RotationAcc = 0.
            BlankMesh%FieldMask(MASKID_ROTATIONACC) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF


      BlankMesh%nScalars = 0
      IF ( PRESENT(nScalars) ) THEN
         IF ( nScalars .GT. 0 ) THEN
            CALL AllocAry( BlankMesh%Scalars, nScalars, Nnodes, 'MeshCreate: Scalars', ErrStat,ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%Scalars = 0.
            BlankMesh%FieldMask(MASKID_Scalar) = .TRUE.
            BlankMesh%nScalars = nScalars
            !IsMotion = .TRUE. !bjj: do we care about this one?
         ENDIF
      ENDIF


         !........................
         ! Let's make sure that we have all the necessary fields:
         !........................
      IF ( IsMotion .AND. IOS == COMPONENT_INPUT ) THEN

         IF ( .NOT. BlankMesh%FieldMask(MASKID_Orientation) ) THEN
            CALL AllocAry( BlankMesh%Orientation, 3, 3, Nnodes, 'MeshCreate: Orientation', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            CALL Eye(BlankMesh%Orientation, ErrStat, ErrMess)  ! set this orientation to the identity matrix
            BlankMesh%FieldMask(MASKID_ORIENTATION) = .TRUE.
         !   ErrStat = ErrID_Info
         !   ErrMess = ' MeshCreate: Meshes with motion fields must also contain the Orientation field.'
         END IF

         IF ( .NOT. BlankMesh%FieldMask(MASKID_TRANSLATIONDISP)) THEN
            CALL AllocAry( BlankMesh%TranslationDisp, 3, Nnodes, 'MeshCreate: TranslationDisp', ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationDisp = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONDISP) = .TRUE.
         !   ErrStat = ErrID_Info
         !   ErrMess = ' MeshCreate: Meshes with motion fields must also contain the TranslationDisp field.'
         ENDIF

      END IF

      RETURN

   END SUBROUTINE MeshCreate

   RECURSIVE SUBROUTINE MeshDestroy ( Mesh, ErrStat, ErrMess, IgnoreSibling )
     TYPE(MeshType),  INTENT(INOUT) :: Mesh            ! Mesh to be vaporized
     INTEGER(IntKi),  INTENT(OUT)   :: ErrStat         ! Error status/code
     CHARACTER(*),    INTENT(OUT)   :: ErrMess         ! Error message
    ! On a brand new mesh, the pointers to siblings may not be nullified and
    ! thus undefined (which may cause ASSOCIATED to report .true. erroneously)
    ! This despite use of => NULL in declaration of this fields for MeshType. Sigh.
    ! So if IgnoreSibling is present and true, don't follow the sibling pointers.
    ! Instead just unconditionally nullify these.
    ! Use this carefully, since it can leave dangling memory if used for a
    ! mesh that already exists and has existing siblings.
     LOGICAL, INTENT(IN), OPTIONAL :: IgnoreSibling

    ! Local
      LOGICAL IgSib
      INTEGER i, j, k

      ErrStat = ErrID_None

      !IF ( .NOT. Mesh%Initialized ) RETURN

         ! Deallocate/Nullify/Deinitialize values that are not shared between siblings:

      Mesh%initialized = .FALSE.
      Mesh%committed   = .FALSE.
      Mesh%fieldmask   = .FALSE.
      Mesh%ios         = 0
      Mesh%Nnodes      = 0

      IF ( ALLOCATED(Mesh%Force)          ) DEALLOCATE(Mesh%Force)
      IF ( ALLOCATED(Mesh%Moment)         ) DEALLOCATE(Mesh%Moment)
      IF ( ALLOCATED(Mesh%Orientation)    ) DEALLOCATE(Mesh%Orientation)
      IF ( ALLOCATED(Mesh%TranslationDisp)) DEALLOCATE(Mesh%TranslationDisp)
      IF ( ALLOCATED(Mesh%RotationVel)    ) DEALLOCATE(Mesh%RotationVel)
      IF ( ALLOCATED(Mesh%TranslationVel) ) DEALLOCATE(Mesh%TranslationVel)
      IF ( ALLOCATED(Mesh%RotationAcc)    ) DEALLOCATE(Mesh%RotationAcc)
      IF ( ALLOCATED(Mesh%TranslationAcc) ) DEALLOCATE(Mesh%TranslationAcc)
      IF ( ALLOCATED(Mesh%Scalars)        ) DEALLOCATE(Mesh%Scalars)
      


!#if 0
     !IF (ASSOCIATED(Mesh%element_point))     NULLIFY(Mesh%element_point)
     !IF (ASSOCIATED(Mesh%element_line2))     NULLIFY(Mesh%element_line2)
     !IF (ASSOCIATED(Mesh%element_line3))     NULLIFY(Mesh%element_line3)
     !IF (ASSOCIATED(Mesh%element_tri3))      NULLIFY(Mesh%element_tri3)
     !IF (ASSOCIATED(Mesh%element_tri6))      NULLIFY(Mesh%element_tri6)
     !IF (ASSOCIATED(Mesh%element_quad4))     NULLIFY(Mesh%element_quad4)
     !IF (ASSOCIATED(Mesh%element_quad8))     NULLIFY(Mesh%element_quad8)
     !IF (ASSOCIATED(Mesh%element_tet4))      NULLIFY(Mesh%element_tet4)
     !IF (ASSOCIATED(Mesh%element_tet10))     NULLIFY(Mesh%element_tet10)
     !IF (ASSOCIATED(Mesh%element_hex8))      NULLIFY(Mesh%element_hex8)
     !IF (ASSOCIATED(Mesh%element_hex20))     NULLIFY(Mesh%element_hex20)
     !IF (ASSOCIATED(Mesh%element_wedge6))    NULLIFY(Mesh%element_wedge6)
     !IF (ASSOCIATED(Mesh%element_wedge15))   NULLIFY(Mesh%element_wedge15)
     !IF (ASSOCIATED(Mesh%Position))          NULLIFY(Mesh%Position)
     !IF (ASSOCIATED(Mesh%Force))             NULLIFY(Mesh%Force)
     !IF (ASSOCIATED(Mesh%Moment))            NULLIFY(Mesh%Moment)
     !IF (ASSOCIATED(Mesh%Orientation))       NULLIFY(Mesh%Orientation)
     !IF (ASSOCIATED(Mesh%TranslationDisp))   NULLIFY(Mesh%TranslationDisp)
     !IF (ASSOCIATED(Mesh%RotationVel))       NULLIFY(Mesh%RotationVel)
     !IF (ASSOCIATED(Mesh%TranslationVel))    NULLIFY(Mesh%TranslationVel)
     !IF (ASSOCIATED(Mesh%RotationAcc))       NULLIFY(Mesh%RotationAcc)
     !IF (ASSOCIATED(Mesh%TranslationAcc))    NULLIFY(Mesh%TranslationAcc)
     !IF (ASSOCIATED(Mesh%Scalars))           NULLIFY(Mesh%Scalars)
!#endif

!bjj: if we keep the sibling, deleting this table is going to be a problem

      IgSib = .FALSE.
      IF ( PRESENT( IgnoreSibling ) ) THEN
         IgSib = IgnoreSibling
      ENDIF


      IF ( .NOT. ASSOCIATED( Mesh%SiblingMesh ) ) THEN ! There is no sibling mesh so we don't want to keep the data.

            ! Deallocate and Nullify all fields that can be shared between siblings

         IF ( ASSOCIATED(Mesh%RemapFlag) ) THEN
            DEALLOCATE(Mesh%RemapFlag)
            NULLIFY(Mesh%RemapFlag)
         END IF

         IF ( ASSOCIATED(Mesh%ElemTable) ) THEN
            DO i = 1, NELEMKINDS
               Mesh%ElemTable(i)%nelem = 0  ; Mesh%ElemTable(i)%maxelem = 0
               IF (ASSOCIATED(Mesh%ElemTable(i)%Elements)) THEN
                  DO j = 1, SIZE(Mesh%ElemTable(i)%Elements)
                     IF (ASSOCIATED(Mesh%ElemTable(i)%Elements(j)%ElemNodes)) THEN
                        DEALLOCATE(Mesh%ElemTable(i)%Elements(j)%ElemNodes)
                        NULLIFY(Mesh%ElemTable(i)%Elements(j)%ElemNodes)
                     ENDIF
                     IF (ASSOCIATED(Mesh%ElemTable(i)%Elements(j)%Neighbors)) THEN
                        DEALLOCATE(Mesh%ElemTable(i)%Elements(j)%Neighbors)
                        NULLIFY(Mesh%ElemTable(i)%Elements(j)%Neighbors)
                     ENDIF
                  ENDDO
                  DEALLOCATE(Mesh%ElemTable(i)%Elements)
                  NULLIFY(Mesh%ElemTable(i)%Elements)
               ENDIF
            ENDDO
            DEALLOCATE(Mesh%ElemTable)
            NULLIFY(Mesh%ElemTable)
         ENDIF

         IF (ASSOCIATED(Mesh%ElemList) ) THEN
            DEALLOCATE(Mesh%ElemList) ! These elements pointed to the ElemTable data, which we've deallocated already
            NULLIFY( Mesh%ElemList )
         END IF

         IF ( ASSOCIATED(Mesh%Position) ) THEN
            DEALLOCATE(Mesh%Position)
            NULLIFY(Mesh%Position)
         END IF

         IF ( ASSOCIATED(Mesh%RefOrientation) ) THEN
            DEALLOCATE(Mesh%RefOrientation)
            NULLIFY(Mesh%RefOrientation)
         END IF


      ELSE ! Keep the data for an existing sibling mesh (nullify but don't deallocate):

         NULLIFY( Mesh%RemapFlag )
         NULLIFY( Mesh%ElemTable )
         NULLIFY( Mesh%ElemList  )
         NULLIFY( Mesh%Position  )
         NULLIFY( Mesh%RefOrientation  )

            ! Tell the sibling that this sibling doesn't exist (avoid endless recursion):
         IF ( ASSOCIATED( Mesh%SiblingMesh%SiblingMesh ) ) THEN ! the mesh should be associated with Mesh%SiblingMesh%SiblingMesh
            NULLIFY( Mesh%SiblingMesh%SiblingMesh )
         ELSE
            ! Throw a fault here. The mesh's twin should point back to this mesh.
            ErrStat = ErrID_Fatal
            ErrMess  = 'ModMesh: MeshDestroy: estranged twin mesh.'
            CALL ProgAbort ( ErrMess ) !bjj: our handbook says we shouldn't call ProgAbort, except in the Glue code
                                       !jm:  but this is framework code, not a contributed module, and this is a fatal error
                                       !bjj: it's called from contributed modules; but, assuming we've done the rest of this correctly, we'd never get an estranged twin mesh, right?
                                       !     If we abort using the FAST-for-Matlab code, Matlab will almost certainly have to restart because files are locked (or it runs out of memory). a pain.
            RETURN
         ENDIF

         IF ( .not. IgSib ) THEN  ! don't Ignore the sibling (i.e., delete it, too)
            CALL MeshDestroy( Mesh%SiblingMesh, ErrStat, ErrMess )
            IF (ErrStat >= AbortErrLev) RETURN
         ENDIF !IgSib

      END IF

      NULLIFY( Mesh%SiblingMesh )


   END SUBROUTINE MeshDestroy

!.............
!bjj: these are unused, I assume they were replaced with parameters MASKID_FORCE and MASKID_MOMENT
!#define HDR_FIELDMASK_FORCE   1
!#define HDR_FIELDMASK_MOMENT  2
!.............

! Format of the Int buffer
!   word
!     1        Total size of Int buffer in bytes
!     2        Total size of Real buffer in bytes
!     3        Total size of Db  buffer in bytes
!     4        IOS
!     5        Number of Nodes
!     6        Number of element records
!     7        FieldMask                           FIELDMASK_SIZE
!     7+$7     Table Entries                       $5 * SIZE(ElemRecType)
!
   SUBROUTINE MeshPack ( Mesh, ReBuf, DbBuf, IntBuf , ErrStat, ErrMess, SizeOnly )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being packed
     REAL(ReKi),     ALLOCATABLE, INTENT(  OUT) :: ReBuf(:)  ! Real buffer
     REAL(DbKi),     ALLOCATABLE, INTENT(  OUT) :: DbBuf(:)  ! Double buffer
     INTEGER(IntKi), ALLOCATABLE, INTENT(  OUT) :: IntBuf(:) ! Int buffer
     INTEGER(IntKi),              INTENT(  OUT) :: ErrStat
     CHARACTER(*),                INTENT(  OUT) :: ErrMess
     LOGICAL,OPTIONAL,            INTENT(IN   ) :: SizeOnly
   ! Local
     INTEGER i,ic,nelem,n_int,n_re,n_db,l,ii,jj,CtrlCode,x
     INTEGER Ielement, Xelement
     TYPE(ElemRecType), POINTER :: ElemRec
     LOGICAL SzOnly

     SzOnly = .FALSE.
     IF ( PRESENT(SizeOnly) ) SzOnly = SizeOnly


!bjj: if we modify the code below to check if things are allocated/etc we can remove this requirement (and then we have to rethink the unpack)
     IF (.NOT. Mesh%Committed ) THEN
        ErrStat = ErrID_Fatal
        ErrMess = " MeshPack: Uncommitted meshes cannot be packed."
        RETURN
     END IF

    ! Traverse the element list and calculate size needed to store element records
     n_int = 0
     nelem = 0
     CtrlCode = 0
     CALL MeshNextElement( Mesh, CtrlCode, ErrStat, ErrMess, Ielement=Ielement, Xelement=Xelement, ElemRec=ElemRec )
     IF (ErrStat >= AbortErrLev) RETURN

     DO WHILE ( CtrlCode .NE. MESH_NOMORE )
       nelem = nelem + 1
       n_int = n_int + 1                       ! add word for element kind
       n_int = n_int + 1                       ! add word for number of nodes
       n_int = n_int + NumNodes( Xelement )    ! space for nodes in this element
#if 0
 TODO
       n_int = n_int + ElemRec%Nneighbors         ! space for neighbor list, stored as element indices
#endif
       CtrlCode = MESH_NEXT
       CALL MeshNextElement( Mesh, CtrlCode, ErrStat, ErrMess, Ielement=Ielement, Xelement=Xelement )
       IF (ErrStat >= AbortErrLev) RETURN

     ENDDO
     n_int = n_int + HDR_FIRSTELEM + FIELDMASK_SIZE  ! add space for header
     ALLOCATE( IntBuf( n_int ) )

     n_re = 0
     n_re = n_re + Mesh%Nnodes * 3 ! Position
     n_re = n_re + Mesh%Nnodes * 9 ! RefOrientation
     IF ( Mesh%FieldMask(MASKID_FORCE) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_MOMENT) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_ORIENTATION) ) n_re = n_re + Mesh%Nnodes * 9
     IF ( Mesh%FieldMask(MASKID_ROTATIONVEL) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_ROTATIONACC) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_TRANSLATIONDISP) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_TRANSLATIONVEL) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_TRANSLATIONACC) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%nScalars .GT. 0 ) n_re = n_re + Mesh%Nnodes * Mesh%nScalars

     ALLOCATE( ReBuf( n_re ) )

     n_db = 0  ! unused for now

     IF ( .NOT. SzOnly ) THEN

       ! Actually pack the Int buffer
       IntBuf(HDR_INTBUFSIZE) = BYTES_IN_INT  * n_int
       IntBuf(HDR_REALBUFSIZE) = BYTES_IN_REAL * n_re
       IntBuf(HDR_DBLBUFSIZE) = BYTES_IN_DBL  * n_db
       IntBuf(HDR_IOS) = Mesh%IOS
       IntBuf(HDR_NUMNODES) = Mesh%Nnodes
       IntBuf(HDR_NUMELEMREC) = nelem
       DO i =0,FIELDMASK_SIZE-1
         x = 0 ; IF ( Mesh%FieldMask(i+1) ) x = 1   ! convert logical to ints
         IntBuf(HDR_FIELDMASK+i) = x
       ENDDO
       !ic = FIELDMASK_SIZE+HDR_FIXEDLEN+1
       ic = HDR_FIRSTELEM
       CtrlCode = 0
       CALL MeshNextElement( Mesh, CtrlCode, ErrStat, ErrMess, Ielement=Ielement, Xelement=Xelement, ElemRec=ElemRec )
       IF (ErrStat >= AbortErrLev) RETURN

       DO WHILE ( CtrlCode .NE. MESH_NOMORE )
         IntBuf(ic) = Xelement ; ic = ic + 1
         IntBuf(ic) = NumNodes(Xelement) ; ic = ic + 1
         DO i = 1,NumNodes(Xelement)
           IntBuf(ic) = ElemRec%ElemNodes(i)
           ic = ic + 1
         ENDDO
#if 0
 TODO
         DO i = 1,ElemRec%Nneighbors
           IntBuf(ic) = 0 !   ElemRec%Neighbors(i)
           ic = ic + 1
         ENDDO
#endif
         CtrlCode = MESH_NEXT
         CALL MeshNextElement( Mesh, CtrlCode, ErrStat, ErrMess, Ielement=Ielement, Xelement=Xelement, ElemRec=ElemRec )
         IF (ErrStat >= AbortErrLev) RETURN

       ENDDO

       ! Actually pack the Real buffer
       ic = 1
       DO i = 1, Mesh%Nnodes
         ReBuf(ic) = Mesh%Position(1,i) ; ic = ic + 1
         ReBuf(ic) = Mesh%Position(2,i) ; ic = ic + 1
         ReBuf(ic) = Mesh%Position(3,i) ; ic = ic + 1
       ENDDO


       DO i = 1, Mesh%Nnodes
         DO jj = 1,3
            DO ii = 1,3
               ReBuf(ic) = Mesh%RefOrientation(ii,jj,i) ; ic = ic + 1
            ENDDO
         ENDDO
       ENDDO

       IF ( Mesh%FieldMask(MASKID_FORCE) ) THEN ! n_re = n_re + Mesh%Nnodes * 3
         DO i = 1, Mesh%Nnodes
           ReBuf(ic) = Mesh%Force(1,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%Force(2,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%Force(3,i) ; ic = ic + 1
         ENDDO
       ENDIF
       IF ( Mesh%FieldMask(MASKID_MOMENT) ) THEN ! n_re = n_re + Mesh%Nnodes * 3
         DO i = 1, Mesh%Nnodes
           ReBuf(ic) = Mesh%Moment(1,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%Moment(2,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%Moment(3,i) ; ic = ic + 1
         ENDDO
       ENDIF
       IF ( Mesh%FieldMask(MASKID_ORIENTATION) ) THEN ! n_re = n_re + Mesh%Nnodes * 9
         DO i = 1, Mesh%Nnodes
           DO jj = 1,3
             DO ii = 1,3
               ReBuf(ic) = Mesh%Orientation(ii,jj,i) ; ic = ic + 1
             ENDDO
           ENDDO
         ENDDO
       ENDIF
       IF ( Mesh%FieldMask(MASKID_TRANSLATIONDISP) ) THEN ! n_re = n_re + Mesh%Nnodes * 3
         DO i = 1, Mesh%Nnodes
           ReBuf(ic) = Mesh%TranslationDisp(1,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%TranslationDisp(2,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%TranslationDisp(3,i) ; ic = ic + 1
         ENDDO
       ENDIF
       IF ( Mesh%FieldMask(MASKID_ROTATIONVEL) ) THEN ! n_re = n_re + Mesh%Nnodes * 3
         DO i = 1, Mesh%Nnodes
           ReBuf(ic) = Mesh%RotationVel(1,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%RotationVel(2,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%RotationVel(3,i) ; ic = ic + 1
         ENDDO
       ENDIF
       IF ( Mesh%FieldMask(MASKID_TRANSLATIONVEL) ) THEN ! n_re = n_re + Mesh%Nnodes * 3
         DO i = 1, Mesh%Nnodes
           ReBuf(ic) = Mesh%TranslationVel(1,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%TranslationVel(2,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%TranslationVel(3,i) ; ic = ic + 1
         ENDDO
       ENDIF
       IF ( Mesh%FieldMask(MASKID_ROTATIONACC) ) THEN ! n_re = n_re + Mesh%Nnodes * 3
         DO i = 1, Mesh%Nnodes
           ReBuf(ic) = Mesh%RotationAcc(1,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%RotationAcc(2,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%RotationAcc(3,i) ; ic = ic + 1
         ENDDO
       ENDIF
       IF ( Mesh%FieldMask(MASKID_TRANSLATIONACC) ) THEN ! n_re = n_re + Mesh%Nnodes * 3
         DO i = 1, Mesh%Nnodes
           ReBuf(ic) = Mesh%TranslationAcc(1,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%TranslationAcc(2,i) ; ic = ic + 1
           ReBuf(ic) = Mesh%TranslationAcc(3,i) ; ic = ic + 1
         ENDDO
       ENDIF
       IF ( Mesh%nScalars .GT. 0 ) THEN ! n_re = n_re + Mesh%Nnodes * Mesh%nScalar
         DO i = 1, Mesh%Nnodes
           DO ii = 1,Mesh%nScalars
             ReBuf(ic) = Mesh%Scalars(ii,i) ; ic = ic + 1
           ENDDO
         ENDDO
       ENDIF
     ENDIF
     ALLOCATE( DbBuf( 1 ) )  ! dummy allocation

     !bjj: where are we keeping track of which ones are siblings so that we can unpack them (set pointers) properly for restart?
   END SUBROUTINE MeshPack

   SUBROUTINE MeshUnpack( Mesh, Re_Buf, Db_Buf, Int_Buf, ErrStat, ErrMess )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh
     REAL(ReKi),     ALLOCATABLE, INTENT(IN   ) :: Re_Buf(:)
     REAL(DbKi),     ALLOCATABLE, INTENT(IN   ) :: Db_Buf(:)
     INTEGER(IntKi), ALLOCATABLE, INTENT(IN   ) :: Int_Buf(:)
     INTEGER(IntKi),              INTENT(  OUT) :: ErrStat
     CHARACTER(*),                INTENT(  OUT) :: ErrMess

   ! Local
     LOGICAL Force, Moment, Orientation, TranslationDisp, TranslationVel, RotationVel, TranslationAcc, RotationAcc, AddedMass
     INTEGER nScalars
     INTEGER i,ic,nelem,n_int,n_re,n_db,l,ii,jj,CtrlCode,x,nelemnodes
     INTEGER Ielement, Xelement
     TYPE(ElemRecType), POINTER :: ElemRec

     Force = .FALSE.
     Moment = .FALSE.
     Orientation = .FALSE.
     TranslationDisp = .FALSE.
     TranslationVel = .FALSE.
     RotationVel = .FALSE.
     TranslationAcc = .FALSE.
     RotationAcc = .FALSE.
     nScalars = 0
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_FORCE-1)       .EQ. 1 ) Force = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_MOMENT-1)      .EQ. 1 ) Moment = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_ORIENTATION-1) .EQ. 1 ) Orientation = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_TRANSLATIONDISP-1) .EQ. 1 ) TranslationDisp = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_TRANSLATIONVEL-1) .EQ. 1 ) TranslationVel = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_ROTATIONVEL-1)    .EQ. 1 ) RotationVel = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_TRANSLATIONACC-1) .EQ. 1 ) TranslationAcc = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_ROTATIONACC-1)    .EQ. 1 ) RotationAcc = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_SCALAR-1)      .GT. 0 ) nScalars = Int_Buf(HDR_FIELDMASK+MASKID_SCALAR-1)

!write(0,*)'Int_Buf(HDR_INTBUFSIZE) ',Int_Buf(HDR_INTBUFSIZE)
!write(0,*)'Int_Buf(HDR_IOS) ',Int_Buf(HDR_IOS)
!write(0,*)'Int_Buf(HDR_NUMNODES) ',Int_Buf(HDR_NUMNODES)
     CALL MeshCreate( Mesh, Int_Buf(HDR_IOS), Int_Buf(HDR_NUMNODES)                    &
                     ,ErrStat=ErrStat, ErrMess=ErrMess                                 &
                     ,Force=Force, Moment=Moment, Orientation=Orientation              &
                     ,TranslationDisp=TranslationDisp                                  &
                     ,TranslationVel=TranslationVel, RotationVel=RotationVel           &
                     ,TranslationAcc=TranslationAcc, RotationAcc=RotationAcc           &
                     , nScalars = nScalars                                             &
                    )
     IF (ErrStat >= AbortErrLev) RETURN

!write(0,*)'Int_Buf(HDR_NUMELEMREC) ',Int_Buf(HDR_NUMELEMREC)
     ic = HDR_FIRSTELEM
     DO i = 1, Int_Buf(HDR_NUMELEMREC)
       Xelement = Int_Buf(ic) ; ic = ic + 1
       nelemnodes = Int_Buf(ic) ; ic = ic + 1
!write(0,*)'ic ',ic,' Xelement: ',Xelement,' nelemnodes ',nelemnodes
       SELECT CASE (nelemnodes )
         CASE (1)
           CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess                          &
            , P1 =Int_Buf(ic   )                                                                &
            )
         CASE (2)
           CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess                          &
            , P1 =Int_Buf(ic   ), P2= Int_Buf(ic+1 )                                            &
            )
         CASE (3)
           CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess                          &
            , P1 =Int_Buf(ic   ), P2 =Int_Buf(ic+1 ), P3 =Int_Buf(ic+2 )                        &
            )
         CASE (4)
           CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess                          &
            , P1 =Int_Buf(ic   ), P2 =Int_Buf(ic+1 ), P3 =Int_Buf(ic+2 ), P4 =Int_Buf(ic+3 )    &
            )
         CASE (6)
           CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess                          &
            , P1 =Int_Buf(ic   ), P2 =Int_Buf(ic+1 ), P3 =Int_Buf(ic+2 ), P4 =Int_Buf(ic+3 )    &
            , P5 =Int_Buf(ic+4 ), P6 =Int_Buf(ic+5 )                                            &
            )
         CASE (8)
           CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess                          &
            , P1 =Int_Buf(ic   ), P2 =Int_Buf(ic+1 ), P3 =Int_Buf(ic+2 ), P4 =Int_Buf(ic+3 )    &
            , P5 =Int_Buf(ic+4 ), P6 =Int_Buf(ic+5 ), P7 =Int_Buf(ic+6 ), P8 =Int_Buf(ic+7 )    &
            )
         CASE (10)
           CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess                          &
            , P1 =Int_Buf(ic   ), P2 =Int_Buf(ic+1 ), P3 =Int_Buf(ic+2 ), P4 =Int_Buf(ic+3 )    &
            , P5 =Int_Buf(ic+4 ), P6 =Int_Buf(ic+5 ), P7 =Int_Buf(ic+6 ), P8 =Int_Buf(ic+7 )    &
            , P9 =Int_Buf(ic+8 ), P10=Int_Buf(ic+9 )                                            &
            )
         CASE (15)
           CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess                          &
            , P1 =Int_Buf(ic   ), P2 =Int_Buf(ic+1 ), P3 =Int_Buf(ic+2 ), P4 =Int_Buf(ic+3 )    &
            , P5 =Int_Buf(ic+4 ), P6 =Int_Buf(ic+5 ), P7 =Int_Buf(ic+6 ), P8 =Int_Buf(ic+7 )    &
            , P9 =Int_Buf(ic+8 ), P10=Int_Buf(ic+9 ), P11=Int_Buf(ic+10), P12=Int_Buf(ic+11)    &
            , P13=Int_Buf(ic+12), P14=Int_Buf(ic+13), P15=Int_Buf(ic+14)                        &
            )
         CASE (20)
           CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess                          &
            , P1 =Int_Buf(ic   ), P2 =Int_Buf(ic+1 ), P3 =Int_Buf(ic+2 ), P4 =Int_Buf(ic+3 )    &
            , P5 =Int_Buf(ic+4 ), P6 =Int_Buf(ic+5 ), P7 =Int_Buf(ic+6 ), P8 =Int_Buf(ic+7 )    &
            , P9 =Int_Buf(ic+8 ), P10=Int_Buf(ic+9 ), P11=Int_Buf(ic+10), P12=Int_Buf(ic+11)    &
            , P13=Int_Buf(ic+12), P14=Int_Buf(ic+13), P15=Int_Buf(ic+14), P16=Int_Buf(ic+15)    &
            , P17=Int_Buf(ic+16), P18=Int_Buf(ic+17), P19=Int_Buf(ic+18), P20=Int_Buf(ic+19)    &
            )
         CASE DEFAULT
           ! Throw a fault here. There is no such element. Probably buffer mangled.
            ErrStat = ErrID_Fatal
            ErrMess = 'ModMesh: MeshUnpack: No such element.  Probably manged buffer.'
            CALL ProgAbort ( ErrMess )
            RETURN
      END SELECT
      IF (ErrStat >= AbortErrLev) RETURN

       ic = ic + nelemnodes ;

     ENDDO

     ! Unpack the Real Buffer
     ic = 1
     DO i = 1, Mesh%Nnodes
       Mesh%Position(1,i) = Re_Buf(ic) ; ic = ic + 1
       Mesh%Position(2,i) = Re_Buf(ic) ; ic = ic + 1
       Mesh%Position(3,i) = Re_Buf(ic) ; ic = ic + 1
     ENDDO

     DO i = 1, Mesh%Nnodes
        DO jj = 1,3
           DO ii = 1,3
             Mesh%RefOrientation(ii,jj,i) = Re_Buf(ic) ; ic = ic + 1
           ENDDO
        ENDDO
     ENDDO

     IF ( Mesh%FieldMask(MASKID_FORCE) ) THEN
       DO i = 1, Mesh%Nnodes
         Mesh%Force(1,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%Force(2,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%Force(3,i) = Re_Buf(ic) ; ic = ic + 1
       ENDDO
     ENDIF
     IF ( Mesh%FieldMask(MASKID_MOMENT) ) THEN
       DO i = 1, Mesh%Nnodes
         Mesh%Moment(1,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%Moment(2,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%Moment(3,i) = Re_Buf(ic) ; ic = ic + 1
       ENDDO
     ENDIF
     IF ( Mesh%FieldMask(MASKID_ORIENTATION) ) THEN
       DO i = 1, Mesh%Nnodes
         DO jj = 1,3
           DO ii = 1,3
             Mesh%Orientation(ii,jj,i) = Re_Buf(ic) ; ic = ic + 1
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     IF ( Mesh%FieldMask(MASKID_TRANSLATIONDISP) ) THEN
       DO i = 1, Mesh%Nnodes
         Mesh%TranslationDisp(1,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%TranslationDisp(2,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%TranslationDisp(3,i) = Re_Buf(ic) ; ic = ic + 1
       ENDDO
     ENDIF
     IF ( Mesh%FieldMask(MASKID_ROTATIONVEL) ) THEN
       DO i = 1, Mesh%Nnodes
         Mesh%RotationVel(1,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%RotationVel(2,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%RotationVel(3,i) = Re_Buf(ic) ; ic = ic + 1
       ENDDO
     ENDIF
     IF ( Mesh%FieldMask(MASKID_TRANSLATIONVEL) ) THEN
       DO i = 1, Mesh%Nnodes
         Mesh%TranslationVel(1,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%TranslationVel(2,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%TranslationVel(3,i) = Re_Buf(ic) ; ic = ic + 1
       ENDDO
     ENDIF
     IF ( Mesh%FieldMask(MASKID_ROTATIONACC) ) THEN
       DO i = 1, Mesh%Nnodes
         Mesh%RotationAcc(1,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%RotationAcc(2,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%RotationAcc(3,i) = Re_Buf(ic) ; ic = ic + 1
       ENDDO
     ENDIF
     IF ( Mesh%FieldMask(MASKID_TRANSLATIONACC) ) THEN
       DO i = 1, Mesh%Nnodes
         Mesh%TranslationAcc(1,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%TranslationAcc(2,i) = Re_Buf(ic) ; ic = ic + 1
         Mesh%TranslationAcc(3,i) = Re_Buf(ic) ; ic = ic + 1
       ENDDO
     ENDIF
     IF ( Mesh%nScalars .GT. 0 ) THEN
       DO i = 1, Mesh%Nnodes
         DO ii = 1,Mesh%nScalars
           Mesh%Scalars(ii,i) = Re_Buf(ic) ; ic = ic + 1
         ENDDO
       ENDDO
     ENDIF

     CALL MeshCommit(Mesh, ErrStat, ErrMess)

     RETURN



   END SUBROUTINE MeshUnpack

   SUBROUTINE MeshCopy( SrcMesh, DestMesh, CtrlCode, ErrStat , ErrMess   &
                      ,IOS, Force, Moment, Orientation, TranslationDisp, TranslationVel &
                      ,RotationVel, TranslationAcc, RotationAcc, nScalars )
     TYPE(MeshType), TARGET,      INTENT(INOUT) :: SrcMesh  ! Mesh being copied
     TYPE(MeshType), TARGET,      INTENT(INOUT) :: DestMesh ! Copy of mesh
     INTEGER(IntKi),              INTENT(IN)    :: CtrlCode ! MESH_NEWCOPY, MESH_SIBLING, or
                                                            ! MESH_UPDATECOPY
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat  ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess  ! Error message
    ! Optional arguments (used only if CtrlCode is MESH_SIBLING):
     INTEGER(IntKi), OPTIONAL,    INTENT(IN)    :: IOS               ! If present, IOS of new sibling: input (COMPONENT_INPUT), output(COMPONENT_OUTPUT), or state(COMPONENT_STATE)
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: Force             ! If present and true, allocate Force field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: Moment            ! If present and true, allocate Moment field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: Orientation       ! If present and true, allocate Orientation field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: TranslationDisp   ! If present and true, allocate TranslationDisp field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: TranslationVel    ! If present and true, allocate TranslationVel field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: RotationVel       ! If present and true, allocate RotationVel field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: TranslationAcc    ! If present and true, allocate TranslationAcc field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: RotationAcc       ! If present and true, allocate RotationAcc field
     INTEGER(IntKi), OPTIONAL,    INTENT(IN)    :: nScalars          ! If present and > 0 , alloc n Scalars
    ! Local
     INTEGER(IntKi)                             :: IOS_l               ! IOS of new sibling
     LOGICAL                                    :: Force_l           & ! If true, allocate Force field
                                                 , Moment_l          & ! If true, allocate Moment field
                                                 , Orientation_l     & ! If true, allocate Orientation field
                                                 , TranslationDisp_l & ! If true, allocate TranslationDisp field
                                                 , TranslationVel_l  & ! If true, allocate TranslationVel field
                                                 , RotationVel_l     & ! If true, allocate RotationVel field
                                                 , TranslationAcc_l  & ! If true, allocate TranslationAcc field
                                                 , RotationAcc_l       ! If true, allocate RotationAcc field
     INTEGER(IntKi)                             :: nScalars_l          ! If > 0, alloc n Scalars
     INTEGER i, j, k


      ErrStat = ErrID_None
      ErrMess = ""

      IF (.NOT. SrcMesh%Initialized) RETURN !bjj: maybe we should first CALL MeshDestroy(DestMesh,ErrStat, ErrMess)

      IF ( CtrlCode .EQ. MESH_NEWCOPY .OR. CtrlCode .EQ. MESH_SIBLING ) THEN

         IF ( CtrlCode .EQ. MESH_NEWCOPY ) THEN
            CALL MeshCreate( DestMesh, IOS=SrcMesh%IOS, Nnodes=SrcMesh%Nnodes, ErrStat=ErrStat, ErrMess=ErrMess &
                            ,Force=SrcMesh%FieldMask(MASKID_FORCE)                                              &
                            ,Moment=SrcMesh%FieldMask(MASKID_MOMENT)                                            &
                            ,Orientation=SrcMesh%FieldMask(MASKID_ORIENTATION)                                  &
                            ,TranslationDisp=SrcMesh%FieldMask(MASKID_TRANSLATIONDISP)                          &
                            ,TranslationVel=SrcMesh%FieldMask(MASKID_TRANSLATIONVEL)                            &
                            ,RotationVel=SrcMesh%FieldMask(MASKID_ROTATIONVEL)                                  &
                            ,TranslationAcc=SrcMesh%FieldMask(MASKID_TRANSLATIONACC)                            &
                            ,RotationAcc=SrcMesh%FieldMask(MASKID_ROTATIONACC)                                  &
                            ,nScalars=SrcMesh%nScalars                                                          )

            IF (ErrStat >= AbortErrLev) RETURN

            DestMesh%Position       =  SrcMesh%Position
            DestMesh%RefOrientation =  SrcMesh%RefOrientation

            DO i = 1, NELEMKINDS
               DestMesh%ElemTable(i)%nelem = SrcMesh%ElemTable(i)%nelem
               DestMesh%ElemTable(i)%maxelem = SrcMesh%ElemTable(i)%maxelem
               DestMesh%ElemTable(i)%XElement = SrcMesh%ElemTable(i)%XElement

               ALLOCATE(DestMesh%ElemTable(i)%Elements(DestMesh%ElemTable(i)%maxelem),STAT=ErrStat)
               IF (ErrStat /=0) THEN
                  ErrStat = ErrID_Fatal
                  ErrMess=' MeshCopy: Error allocating ElemTable%Elements.'
                  RETURN !Early return
               END IF

               DO j = 1,DestMesh%ElemTable(i)%nelem

                  ALLOCATE(DestMesh%ElemTable(i)%Elements(j)%ElemNodes(size(SrcMesh%ElemTable(i)%Elements(j)%ElemNodes)),STAT=ErrStat)
                  IF (ErrStat /=0) THEN
                     ErrStat = ErrID_Fatal
                     ErrMess=' MeshCopy: Error allocating ElemTable%ElemNodes.'
                     RETURN !Early return
                  END IF
                  DestMesh%ElemTable(i)%Elements(j)%ElemNodes = SrcMesh%ElemTable(i)%Elements(j)%ElemNodes
                  DestMesh%ElemTable(i)%Elements(j)%det_jac   = SrcMesh%ElemTable(i)%Elements(j)%det_jac
                  DestMesh%ElemTable(i)%Elements(j)%Xelement  = SrcMesh%ElemTable(i)%Elements(j)%Xelement
                  DestMesh%ElemTable(i)%Elements(j)%Nneighbors= SrcMesh%ElemTable(i)%Elements(j)%Nneighbors
!bjj: allocate Neighbors, too:
!                  IF (DestMesh%ElemTable(i)%Elements(j)%Nneighbors > 0) then

               ENDDO
            ENDDO


               ! Regenerate new list of elements (point to ElemTable)
   !bjj: call meshCommit?
            DestMesh%nelemlist   = SrcMesh%nelemlist
            DestMesh%maxelemlist = SrcMesh%maxelemlist
            DestMesh%nextelem    = SrcMesh%nextelem

            IF ( SrcMesh%Committed ) THEN

               ALLOCATE(DestMesh%ElemList(DestMesh%maxelemlist))
               IF (ErrStat /=0) THEN
                  ErrStat = ErrID_Fatal
                  ErrMess=' MeshCopy: Error allocating ElemList.'
                  RETURN !Early return
               END IF

               k = 0
               DO i = 1, NELEMKINDS
                  DO j = 1, DestMesh%ElemTable(i)%nelem
                     k = k + 1
                     DestMesh%elemlist(k)%Element => DestMesh%ElemTable(i)%Elements(j)
                     DestMesh%elemlist(k)%Element%Xelement = i
                  ENDDO
               END DO

            END IF ! SrcMesh%Committed

            DestMesh%RemapFlag   = SrcMesh%RemapFlag

         ELSE IF ( CtrlCode .EQ. MESH_SIBLING ) THEN
!bjj: we should make sure the mesh has been committed, otherwise the element lists haven't been created, yet (and thus not shared)
            IF ( ASSOCIATED(SrcMesh%SiblingMesh) ) THEN
               ErrStat = ErrID_Fatal
               ErrMess = ' MeshCopy: A mesh can have only one sibling.'
               RETURN !early return
            END IF

            IF (.NOT. SrcMesh%Committed ) THEN
               ErrStat = ErrID_Fatal
               ErrMess = ' MeshCopy: An uncommitted mesh cannot have a sibling.'
               RETURN !early return
            END IF

            IOS_l          = SrcMesh%IOS ; IF ( PRESENT(IOS) )                         IOS_l = IOS
            Force_l            = .FALSE. ; IF ( PRESENT(Force) )                     Force_l = Force
            Moment_l           = .FALSE. ; IF ( PRESENT(Moment) )                   Moment_l = Moment
            Orientation_l      = .FALSE. ; IF ( PRESENT(Orientation) )         Orientation_l = Orientation
            TranslationDisp_l  = .FALSE. ; IF ( PRESENT(TranslationDisp) ) TranslationDisp_l = TranslationDisp
            TranslationVel_l   = .FALSE. ; IF ( PRESENT(TranslationVel) )   TranslationVel_l = TranslationVel
            RotationVel_l      = .FALSE. ; IF ( PRESENT(RotationVel) )         RotationVel_l = RotationVel
            TranslationAcc_l   = .FALSE. ; IF ( PRESENT(TranslationAcc) )   TranslationAcc_l = TranslationAcc
            RotationAcc_l      = .FALSE. ; IF ( PRESENT(RotationAcc) )         RotationAcc_l = RotationAcc
            nScalars_l         = 0       ; IF ( PRESENT(nScalars) )               nScalars_l = nScalars
            CALL MeshCreate( DestMesh, IOS=IOS_l, Nnodes=SrcMesh%Nnodes, ErrStat=ErrStat, ErrMess=ErrMess   &
                            ,Force=Force_l                                                                  &
                            ,Moment=Moment_l                                                                &
                            ,Orientation=Orientation_l                                                      &
                            ,TranslationDisp=TranslationDisp_l                                              &
                            ,TranslationVel=TranslationVel_l                                                &
                            ,RotationVel=RotationVel_l                                                      &
                            ,TranslationAcc=TranslationAcc_l                                                &
                            ,RotationAcc=RotationAcc_l                                                      &
                            ,nScalars=nScalars_l                                                            &
                            ,IsNewSibling=.TRUE.)
            IF (ErrStat >= AbortErrLev) RETURN

            !bjj: Doesn't this logic mean we can have only one sibling?
            ! I added a check that SrcMesh%SiblingMesh isn't already associated so that we don't lose siblings
            DestMesh%SiblingMesh => SrcMesh
            SrcMesh%SiblingMesh  => DestMesh

            DestMesh%Position       => SrcMesh%Position
            DestMesh%RefOrientation => SrcMesh%RefOrientation
            DestMesh%RemapFlag      => SrcMesh%RemapFlag
            DestMesh%ElemTable      => SrcMesh%ElemTable
            DestMesh%ElemList       => SrcMesh%ElemList

            DestMesh%nelemlist   = SrcMesh%nelemlist
            DestMesh%maxelemlist = SrcMesh%maxelemlist
            DestMesh%nextelem    = SrcMesh%nextelem


         ENDIF

         DO i = 1, NELEMKINDS
            IF ( ASSOCIATED(SrcMesh%ElemTable) ) THEN
            ENDIF
            IF ( ASSOCIATED(DestMesh%ElemTable) ) THEN
            ENDIF
         ENDDO


#if 0
          IF ( CtrlCode .EQ. MESH_NEWCOPY ) THEN
             DestMesh%npoint   = SrcMesh%npoint ; DestMesh%maxpoint = SrcMesh%maxpoint
             IF ( SrcMesh%npoint .GT. 0 .AND. ASSOCIATED( SrcMesh%element_point ) ) THEN
      !         CALL AllocPAry( DestMesh%element_point, 1, SrcMesh%maxpoint, 'MeshCreate: element_point' )
               DestMesh%element_point => SrcMesh%element_point
             ENDIF
             DestMesh%nline2   = SrcMesh%nline2 ; DestMesh%maxline2 = SrcMesh%maxline2
             IF ( SrcMesh%nline2 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_line2 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_line2, 2, SrcMesh%maxline2, 'MeshCreate: element_line2' )
               DestMesh%element_line2 => SrcMesh%element_line2
             ENDIF
             DestMesh%nline3   = SrcMesh%nline3 ; DestMesh%maxline3 = SrcMesh%maxline3
             IF ( SrcMesh%nline3 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_line3 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_line3, 3, SrcMesh%maxline3, 'MeshCreate: element_line3' )
               DestMesh%element_line3 => SrcMesh%element_line3
             ENDIF
             DestMesh%ntri3   = SrcMesh%ntri3 ; DestMesh%maxtri3 = SrcMesh%maxtri3
             IF ( SrcMesh%ntri3 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_tri3 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_tri3, 3, SrcMesh%maxtri3, 'MeshCreate: element_tri3' )
               DestMesh%element_tri3 => SrcMesh%element_tri3
             ENDIF
             DestMesh%ntri6   = SrcMesh%ntri6 ; DestMesh%maxtri6 = SrcMesh%maxtri6
             IF ( SrcMesh%ntri6 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_tri6 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_tri6, 6, SrcMesh%maxtri6, 'MeshCreate: element_tri6' )
               DestMesh%element_tri6 => SrcMesh%element_tri6
             ENDIF
             DestMesh%nquad4   = SrcMesh%nquad4 ; DestMesh%maxquad4 = SrcMesh%maxquad4
             IF ( SrcMesh%nquad4 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_quad4 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_quad4, 4, SrcMesh%maxquad4, 'MeshCreate: element_quad4' )
               DestMesh%element_quad4 => SrcMesh%element_quad4
             ENDIF
             DestMesh%nquad8   = SrcMesh%nquad8 ; DestMesh%maxquad8 = SrcMesh%maxquad8
             IF ( SrcMesh%nquad8 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_quad8 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_quad8, 8, SrcMesh%maxquad8, 'MeshCreate: element_quad8' )
               DestMesh%element_quad8 => SrcMesh%element_quad8
             ENDIF
             DestMesh%ntet4   = SrcMesh%ntet4 ; DestMesh%maxtet4 = SrcMesh%maxtet4
             IF ( SrcMesh%ntet4 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_tet4 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_tet4, 4, SrcMesh%maxtet4, 'MeshCreate: element_tet4' )
               DestMesh%element_tet4 => SrcMesh%element_tet4
             ENDIF
             DestMesh%ntet10   = SrcMesh%ntet10 ; DestMesh%maxtet10 = SrcMesh%maxtet10
             IF ( SrcMesh%ntet10 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_tet10 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_tet10, 10, SrcMesh%maxtet10, 'MeshCreate: element_tet10' )
               DestMesh%element_tet10 = SrcMesh%element_tet10
             ENDIF
             DestMesh%nhex8   = SrcMesh%nhex8 ; DestMesh%maxhex8 = SrcMesh%maxhex8
             IF ( SrcMesh%nhex8 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_hex8 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_hex8, 8, SrcMesh%maxhex8, 'MeshCreate: element_hex8' )
               DestMesh%element_hex8 => SrcMesh%element_hex8
             ENDIF
             DestMesh%nhex20   = SrcMesh%nhex20 ; DestMesh%maxhex20 = SrcMesh%maxhex20
             IF ( SrcMesh%nhex20 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_hex20 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_hex20, 20, SrcMesh%maxhex20, 'MeshCreate: element_hex20' )
               DestMesh%element_hex20 => SrcMesh%element_hex20
             ENDIF
             DestMesh%nwedge6   = SrcMesh%nwedge6 ; DestMesh%maxwedge6 = SrcMesh%maxwedge6
             IF ( SrcMesh%nwedge6 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_wedge6 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_wedge6, 6, SrcMesh%maxwedge6, 'MeshCreate: element_wedge6' )
               DestMesh%element_wedge6 => SrcMesh%element_wedge6
             ENDIF
             DestMesh%nwedge15   = SrcMesh%nwedge15 ; DestMesh%maxwedge15 = SrcMesh%maxwedge15
             IF ( SrcMesh%nwedge15 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_wedge15 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_wedge15, 15, SrcMesh%maxwedge15, 'MeshCreate: element_wedge15' )
               DestMesh%element_wedge15 => SrcMesh%element_wedge15
             ENDIF
          ENDIF
#endif
      ELSE IF ( CtrlCode .EQ. MESH_UPDATECOPY ) THEN
         IF ( SrcMesh%nNodes .NE. DestMesh%nNodes ) THEN
            ErrStat = ErrID_Fatal
            ErrMess = "MeshCopy: MESH_UPDATECOPY of meshes with different numbers of nodes."
         ELSE

            ! bjj: should we update positions?
            !DestMesh%RemapFlag = SrcMesh%RemapFlag
            !DestMesh%nextelem  = SrcMesh%nextelem

         ENDIF

      ELSE
         ErrStat = ErrID_Fatal
         ErrMess  = 'MeshCopy: Invalid CtrlCode.'
         RETURN
      ENDIF

         ! These aren't shared between siblings, so they get copied, no matter what the CtrlCode:

      DestMesh%Initialized = SrcMesh%Initialized
      DestMesh%Committed   = SrcMesh%Committed
      IF ( ALLOCATED(SrcMesh%Force          ) .AND. ALLOCATED(DestMesh%Force          ) ) DestMesh%Force = SrcMesh%Force
      IF ( ALLOCATED(SrcMesh%Moment         ) .AND. ALLOCATED(DestMesh%Moment         ) ) DestMesh%Moment = SrcMesh%Moment
      IF ( ALLOCATED(SrcMesh%Orientation    ) .AND. ALLOCATED(DestMesh%Orientation    ) ) DestMesh%Orientation = SrcMesh%Orientation
      IF ( ALLOCATED(SrcMesh%TranslationDisp) .AND. ALLOCATED(DestMesh%TranslationDisp) ) DestMesh%TranslationDisp = SrcMesh%TranslationDisp
      IF ( ALLOCATED(SrcMesh%TranslationVel ) .AND. ALLOCATED(DestMesh%TranslationVel ) ) DestMesh%TranslationVel = SrcMesh%TranslationVel
      IF ( ALLOCATED(SrcMesh%RotationVel    ) .AND. ALLOCATED(DestMesh%RotationVel    ) ) DestMesh%RotationVel = SrcMesh%RotationVel
      IF ( ALLOCATED(SrcMesh%TranslationAcc ) .AND. ALLOCATED(DestMesh%TranslationAcc ) ) DestMesh%TranslationAcc = SrcMesh%TranslationAcc
      IF ( ALLOCATED(SrcMesh%RotationAcc    ) .AND. ALLOCATED(DestMesh%RotationAcc    ) ) DestMesh%RotationAcc = SrcMesh%RotationAcc
      IF ( ALLOCATED(SrcMesh%Scalars        ) .AND. ALLOCATED(DestMesh%Scalars        ) ) DestMesh%Scalars = SrcMesh%Scalars


      !DestMesh%spatial = SrcMesh%spatial !bjj: unused?

   END SUBROUTINE MeshCopy

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshPositionNode( Mesh, Inode, Pos, ErrStat, ErrMess, Orient )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh         ! Mesh being spatio-located
     INTEGER(IntKi),              INTENT(IN   ) :: Inode        ! Number of node being located
     REAL(ReKi),                  INTENT(IN   ) :: Pos(3)       ! Xi,Yi,Zi, coordinates of node
     INTEGER(IntKi),              INTENT(  OUT) :: ErrStat      ! Error code
     CHARACTER(*),                INTENT(  OUT) :: ErrMess      ! Error message
     REAL(ReKi), OPTIONAL,        INTENT(IN   ) :: Orient(3,3)  ! Orientation (direction cosine matrix) of node; identity by default

     ErrStat = ErrID_None
     ErrMess = ""
    ! Safety first
     IF ( .NOT. Mesh%Initialized ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: attempt to use uncreated mesh."
     ENDIF
     IF ( Mesh%Nnodes < 0 ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: Invalid number of nodes: "//TRIM(Num2LStr(Mesh%Nnodes))
     ENDIF
     IF (Mesh%Committed ) THEN  !bjj: perhaps this shouldn't be an error? Maybe it just sets Mesh%RemapFlag = .TRUE.
        ErrStat = ErrID_Fatal
        ErrMess = " MeshPositionNode: attempt to reposition committed mesh."
     END IF

     IF ( .NOT. ( Inode .GE. 1 .AND. Inode .LE. Mesh%Nnodes ) ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: not 1 <= Inode("//TRIM(Num2LStr(Inode))//") <= "//TRIM(Num2LStr(Mesh%Nnodes))
     ENDIF
     IF ( .NOT. ASSOCIATED(Mesh%Position) ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: Position array not associated"
     ENDIF
     IF ( .NOT. SIZE(Mesh%Position,2) .EQ. Mesh%Nnodes ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: Position array not big enough"
     ENDIF
     IF ( .NOT. ASSOCIATED(Mesh%RefOrientation) ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: RefOrientation array not associated"
     ENDIF
     IF ( .NOT. SIZE(Mesh%RefOrientation,3) .EQ. Mesh%Nnodes ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: RefOrientation array not big enough"
     ENDIF



     IF ( ErrStat .NE. ErrID_None ) RETURN   ! early return on error

    ! Business
     Mesh%Position(:,Inode) = Pos

     IF ( PRESENT(Orient) ) THEN
        Mesh%RefOrientation(:,:,Inode) = Orient
     ELSE
        Mesh%RefOrientation(:,1,Inode) = (/ 1., 0., 0. /)
        Mesh%RefOrientation(:,2,Inode) = (/ 0., 1., 0. /)
        Mesh%RefOrientation(:,3,Inode) = (/ 0., 0., 1. /)
     END IF

     RETURN

   END SUBROUTINE MeshPositionNode

   SUBROUTINE MeshCommit( Mesh, ErrStat, ErrMess )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh ! Mesh being committed
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess ! Error message
    ! Local
     REAL(ReKi)                                 :: n1_n2_vector(3)               ! vector going from node 1 to node 2 in a Line2 element
     INTEGER n0d, n1d, n2d, n3d, nn, i, j, NElem, n1,n2
     LOGICAL                                    :: NodeInElement(Mesh%NNodes)    ! determines if each node is part of an element

     !TYPE(ElemListType), POINTER :: tmp(:)

     ! Check for spatial constraints -- can't mix 1D with 2D with 3D
     n0d = Mesh%ElemTable(ELEMENT_POINT)%nelem
     n1d = Mesh%ElemTable(ELEMENT_LINE2)%nelem+Mesh%ElemTable(ELEMENT_LINE3)%nelem
     n2d = Mesh%ElemTable(ELEMENT_TRI3)%nelem+Mesh%ElemTable(ELEMENT_TRI6)%nelem +     &
           Mesh%ElemTable(ELEMENT_QUAD4)%nelem+Mesh%ElemTable(ELEMENT_QUAD8)%nelem
     n3d = Mesh%ElemTable(ELEMENT_TET4)%nelem+Mesh%ElemTable(ELEMENT_TET10)%nelem +    &
           Mesh%ElemTable(ELEMENT_HEX8)%nelem+Mesh%ElemTable(ELEMENT_HEX20)%nelem +    &
           Mesh%ElemTable(ELEMENT_WEDGE6)%nelem+Mesh%ElemTable(ELEMENT_WEDGE15)%nelem
     nn = n0d + n1d + n2d + n3d
     IF ( max(n0d,n1d,n2d,n3d) .LT. nn ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshCommit: mixing elements of different spatial dimensionality"
       RETURN  ! Early return
     ENDIF

     IF ( ALL( Mesh%ElemTable(:)%nelem /= SUM(Mesh%ElemTable(:)%nelem ) ) ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshCommit: a mesh can have only one type of element."
       RETURN  ! Early return
     ENDIF

     ! Check that every node is part of an element.
      NodeInElement = .FALSE.
      DO i = 1, NELEMKINDS
         DO j = 1, Mesh%ElemTable(i)%nelem
            DO n1 = 1,NumNodes( Mesh%ElemTable(i)%XElement )
               n2 = Mesh%ElemTable(i)%Elements(j)%ElemNodes(n1)
               NodeInElement( n2 ) = .TRUE.
            END DO
         END DO
      END DO

      IF ( .NOT. ALL(NodeInElement) ) THEN
         ErrStat = ErrID_Fatal
         ErrMess = "MeshCommit: all nodes must be part of an element."
         RETURN  ! Early return
      END IF


   !bjj: check that allocated fields make sense (i.e., Motions must include Orientation and TranslationDisp)

   !BJJ: initialize allocated orientation field


     ! Construct list of elements


      ! first determine how many elements there are
     Mesh%Nelemlist = 0
     DO i = 1, NELEMKINDS
       DO j = 1, Mesh%ElemTable(i)%nelem
         Mesh%Nelemlist = Mesh%Nelemlist + 1
       END DO
     END DO

     Mesh%maxelemlist = Mesh%Nelemlist
     ALLOCATE( Mesh%elemlist(Mesh%maxelemlist), STAT=ErrStat ) !Allocates the pointer array

     IF (ErrStat /= 0) THEN
        ErrStat = ErrID_Fatal
        ErrMess = " MeshCommit: Error allocating element list."
        RETURN
     END IF

     NElem = 0
     DO i = 1, NELEMKINDS
       DO j = 1, Mesh%ElemTable(i)%nelem
          NElem = NElem + 1
          Mesh%elemlist(NElem)%Element => Mesh%ElemTable(i)%Elements(j)
          Mesh%elemlist(NElem)%Element%Xelement = i
       END DO
     END DO


     ! calculate det_jac:

     DO j=1,Mesh%ElemTable(ELEMENT_POINT)%nelem
         Mesh%ElemTable(ELEMENT_POINT)%Elements(J)%det_jac  = 0.0_ReKi
     END DO

     DO j = 1,Mesh%ElemTable(ELEMENT_LINE2)%nelem

        n1_n2_vector = Mesh%Position(:,Mesh%ElemTable(ELEMENT_LINE2)%Elements(j)%ElemNodes(2)) &
                     - Mesh%Position(:,Mesh%ElemTable(ELEMENT_LINE2)%Elements(j)%ElemNodes(1))

        Mesh%ElemTable(ELEMENT_LINE2)%Elements(J)%det_jac  = 0.5_ReKi * SQRT( DOT_PRODUCT(n1_n2_vector,n1_n2_vector) )   ! = L / 2

     END DO

     ! we're finished:

      Mesh%Committed = .TRUE.

   END SUBROUTINE MeshCommit

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshConstructElement_1PT( Mesh, Xelement, ErrStat, ErrMess, P1 )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER,                     INTENT(IN   ) :: P1
    ! Local
     TYPE(ElemRecType), POINTER :: tmp(:)

     IF ( mesh_debug ) print*,'Called MeshConstructElement_1PT'
     ErrStat = ErrID_None
     ErrMess = ''
    ! Safety first
     IF ( .NOT. Mesh%Initialized ) THEN
       ErrStat = ErrID_Fatal
       ErrMess="MeshConstructElement_1PT: attempt to use uncreated mesh."
     ELSEIF ( P1 .LT. 1 .OR. P1 .GT. Mesh%Nnodes ) THEN !BJJ moved to ELSE
       ErrStat = ErrID_Fatal
       ErrMess="MeshConstructElement_1PT: invalid P1 ("//TRIM(Num2LStr(P1))//") for mesh with "//TRIM(Num2LStr(Mesh%Nnodes))//" nodes."
     ELSEIF (Mesh%Committed ) THEN
        ErrStat = ErrID_Fatal
        ErrMess = " MeshConstructElement_1PT: attempt to add element to committed mesh."
     ENDIF
     IF ( ErrStat .NE. ErrID_None ) THEN
        CALL WrScr(TRIM(ErrMess))
        RETURN  !  early return on error
     ENDIF
    ! Business
     IF ( Xelement .EQ. ELEMENT_POINT ) THEN
       Mesh%ElemTable(ELEMENT_POINT)%nelem = Mesh%ElemTable(ELEMENT_POINT)%nelem + 1
       Mesh%ElemTable(ELEMENT_POINT)%XElement = ELEMENT_POINT

       IF ( Mesh%ElemTable(ELEMENT_POINT)%nelem .GE. Mesh%ElemTable(ELEMENT_POINT)%maxelem ) THEN
!write(0,*)'>>>>>>>>>> bumping maxpoint',Mesh%ElemTable(ELEMENT_POINT)%maxelem
         ALLOCATE(tmp(Mesh%ElemTable(ELEMENT_POINT)%maxelem+BUMPUP))
         IF ( Mesh%ElemTable(ELEMENT_POINT)%nelem .GT. 1 ) &
           tmp(1:Mesh%ElemTable(ELEMENT_POINT)%maxelem) = Mesh%ElemTable(ELEMENT_POINT)%Elements(1:Mesh%ElemTable(ELEMENT_POINT)%maxelem)
         IF ( ASSOCIATED(Mesh%ElemTable(ELEMENT_POINT)%Elements) ) DEALLOCATE(Mesh%ElemTable(ELEMENT_POINT)%Elements)
         Mesh%ElemTable(ELEMENT_POINT)%Elements => tmp
         Mesh%ElemTable(ELEMENT_POINT)%maxelem = Mesh%ElemTable(ELEMENT_POINT)%maxelem + BUMPUP
!write(0,*)'>>>>>>>>>> bumped maxpoint',Mesh%ElemTable(ELEMENT_POINT)%maxelem
       ENDIF

       ALLOCATE(Mesh%ElemTable(ELEMENT_POINT)%Elements(Mesh%ElemTable(ELEMENT_POINT)%nelem)%ElemNodes(1))
       Mesh%ElemTable(ELEMENT_POINT)%Elements(Mesh%ElemTable(ELEMENT_POINT)%nelem)%ElemNodes(1) = P1
       Mesh%ElemTable(ELEMENT_POINT)%Elements(Mesh%ElemTable(ELEMENT_POINT)%nelem)%Xelement = ELEMENT_POINT
     ELSE
       ErrStat = ErrID_Fatal
       ErrMess = 'MeshConstructElement_1PT called for invalid element type'
     ENDIF

     RETURN

   END SUBROUTINE MeshConstructElement_1PT

   SUBROUTINE MeshConstructElement_2PT( Mesh, Xelement, ErrStat, ErrMess, P1, P2 )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER,                     INTENT(IN   ) :: P1,  P2
    ! Local
     TYPE(ElemRecType), POINTER :: tmp(:)

     IF ( mesh_debug ) print*,'Called MeshConstructElement_2PT'
     ErrStat = ErrID_None
     ErrMess = ""
    ! Safety first
     IF ( .NOT. Mesh%Initialized ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshConstructElement_1PT: attempt to use uncreated mesh."
     ELSEIF ( P1 .LT. 1 .OR. P1 .GT. Mesh%Nnodes .OR. &
          P2 .LT. 1 .OR. P2 .GT. Mesh%Nnodes ) THEN
       ErrStat = ErrID_Fatal
       ErrMess ="MeshConstructElement_2PT: invalid P1 ("//TRIM(Num2LStr(P1))//") or P2 ("//TRIM(Num2LStr(P2))//&
                     ") for mesh with "//TRIM(Num2LStr(Mesh%Nnodes))//" nodes."
     ELSEIF (Mesh%Committed ) THEN
        ErrStat = ErrID_Fatal
        ErrMess = " MeshConstructElement_2PT: attempt to add element to committed mesh."
     ENDIF
! TODO: hook the element into a list of elements stored in this mesh to
! allow traversal over elements
     IF ( ErrStat .NE. ErrID_None ) THEN
        CALL WrScr( TRIM(ErrMess) )
        RETURN  !  early return on error
     ENDIF
    ! Business
     IF ( Xelement .EQ. ELEMENT_LINE2 ) THEN
       Mesh%ElemTable(ELEMENT_LINE2)%nelem = Mesh%ElemTable(ELEMENT_LINE2)%nelem + 1
       Mesh%ElemTable(ELEMENT_LINE2)%XElement = ELEMENT_LINE2

       IF ( Mesh%ElemTable(ELEMENT_LINE2)%nelem .GE. Mesh%ElemTable(ELEMENT_LINE2)%maxelem ) THEN
!write(0,*)'>>>>>>>>>> bumping maxline2',Mesh%ElemTable(ELEMENT_LINE2)%maxelem
         ALLOCATE(tmp(Mesh%ElemTable(ELEMENT_LINE2)%maxelem+BUMPUP))
         IF ( Mesh%ElemTable(ELEMENT_LINE2)%nelem .GT. 1 ) &
           tmp(1:Mesh%ElemTable(ELEMENT_LINE2)%maxelem) = Mesh%ElemTable(ELEMENT_LINE2)%Elements(1:Mesh%ElemTable(ELEMENT_LINE2)%maxelem)
         IF ( ASSOCIATED(Mesh%ElemTable(ELEMENT_LINE2)%Elements) ) DEALLOCATE(Mesh%ElemTable(ELEMENT_LINE2)%Elements)
         Mesh%ElemTable(ELEMENT_LINE2)%Elements => tmp
         Mesh%ElemTable(ELEMENT_LINE2)%maxelem = Mesh%ElemTable(ELEMENT_LINE2)%maxelem + BUMPUP
!write(0,*)'>>>>>>>>>> bumped maxline2',Mesh%ElemTable(ELEMENT_LINE2)%maxelem
       ENDIF
       ALLOCATE(Mesh%ElemTable(ELEMENT_LINE2)%Elements(Mesh%ElemTable(ELEMENT_LINE2)%nelem)%ElemNodes(2))
       Mesh%ElemTable(ELEMENT_LINE2)%Elements(Mesh%ElemTable(ELEMENT_LINE2)%nelem)%ElemNodes(1) = P1
       Mesh%ElemTable(ELEMENT_LINE2)%Elements(Mesh%ElemTable(ELEMENT_LINE2)%nelem)%ElemNodes(2) = P2
       Mesh%ElemTable(ELEMENT_LINE2)%Elements(Mesh%ElemTable(ELEMENT_LINE2)%nelem)%Xelement = ELEMENT_LINE2

     ELSE
       ErrMess = 'MeshConstructElement_2PT called for invalid element type'
       ErrStat = ErrID_Fatal
     ENDIF
     RETURN
   END SUBROUTINE MeshConstructElement_2PT

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshConstructElement_3PT( Mesh, Xelement, ErrStat, ErrMess, P1, P2, P3 )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER,                     INTENT(IN   ) :: P1,  P2,  P3
     IF ( mesh_debug ) print*,'Called MeshConstructElement_3PT'
     ErrStat = ErrID_None
     ErrMess = ''
     ErrStat = ErrID_Fatal
     ErrMess = 'MeshConstructElement_3PT not supported'
   END SUBROUTINE MeshConstructElement_3PT

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshConstructElement_4PT( Mesh, Xelement, ErrStat, ErrMess, P1, P2, P3, P4 )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER,                     INTENT(IN   ) :: P1,  P2,  P3,  P4
     IF ( mesh_debug ) print*,'Called MeshConstructElement_4PT'
     ErrStat = ErrID_None
     ErrMess = ''
     ErrStat = ErrID_Fatal
     ErrMess = 'MeshConstructElement_4PT not supported'
   END SUBROUTINE MeshConstructElement_4PT

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshConstructElement_6PT( Mesh, Xelement, ErrStat, ErrMess, P1, P2, P3, P4, P5, P6 )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER,                     INTENT(IN   ) :: P1,  P2,  P3,  P4,  P5, P6
     IF ( mesh_debug ) print*,'Called MeshConstructElement_6PT'
     ErrStat = ErrID_None
     ErrMess = ''
     ErrStat = ErrID_Fatal
     ErrMess = 'MeshConstructElement_6PT not supported'
   END SUBROUTINE MeshConstructElement_6PT

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshConstructElement_8PT( Mesh, Xelement, ErrStat, ErrMess, P1, P2, P3, P4, P5, P6, P7, P8 )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER,                     INTENT(IN   ) :: P1,  P2,  P3,  P4,  P5,  P6,  P7,  P8
     IF ( mesh_debug ) print*,'Called MeshConstructElement_8PT'
     ErrStat = ErrID_None
     ErrMess = ''
     ErrStat = ErrID_Fatal
     ErrMess = 'MeshConstructElement_8PT not supported'
   END SUBROUTINE MeshConstructElement_8PT

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshConstructElement_10PT( Mesh, Xelement, ErrStat, ErrMess, P1, P2, P3, P4, P5, &
                                                                           P6, P7, P8, P9, P10 )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER,                     INTENT(IN   ) :: P1,  P2,  P3,  P4,  P5,  &
                                                   P6,  P7,  P8,  P9,  P10
     IF ( mesh_debug ) print*,'Called MeshConstructElement_10PT'
     ErrStat = ErrID_None
     ErrMess = ''
     ErrStat = ErrID_Fatal
     ErrMess = 'MeshConstructElement_10PT not supported'
   END SUBROUTINE MeshConstructElement_10PT

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshConstructElement_15PT( Mesh, Xelement, ErrStat, ErrMess, P1,  P2,  P3,  P4,  P5,  &
                                                                           P6,  P7,  P8,  P9,  P10, &
                                                                           P11, P12, P13, P14, P15 )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER,                     INTENT(IN   ) :: P1,  P2,  P3,  P4,  P5,  &
                                                   P6,  P7,  P8,  P9,  P10, &
                                                   P11, P12, P13, P14, P15
     IF ( mesh_debug ) print*,'Called MeshConstructElement_15PT'
     ErrStat = ErrID_None
     ErrMess = ''
     ErrStat = ErrID_Fatal
     ErrMess = 'MeshConstructElement_15PT not supported'
   END SUBROUTINE MeshConstructElement_15PT

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshConstructElement_20PT( Mesh, Xelement, ErrStat, ErrMess, P1,  P2,  P3,  P4,  P5,  &
                                                                           P6,  P7,  P8,  P9,  P10, &
                                                                           P11, P12, P13, P14, P15, &
                                                                           P16, P17, P18, P19, P20 )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER,                     INTENT(IN   ) :: P1,  P2,  P3,  P4,  P5,  &
                                                   P6,  P7,  P8,  P9,  P10, &
                                                   P11, P12, P13, P14, P15, &
                                                   P16, P17, P18, P19, P20
     IF ( mesh_debug ) print*,'Called MeshConstructElement_20PT'
     ErrStat = ErrID_None
     ErrMess = ''
     ErrStat = ErrID_Fatal
     ErrMess = 'MeshConstructElement_20PT not supported'
   END SUBROUTINE MeshConstructElement_20PT

   SUBROUTINE MeshNextElement ( Mesh, CtrlCode, ErrStat, ErrMess, Ielement, Xelement, ElemRec )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
     INTEGER(IntKi),              INTENT(INOUT) :: CtrlCode  ! CtrlCode
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
     INTEGER(IntKi),OPTIONAL,     INTENT(OUT)   :: Ielement  ! Element index
     INTEGER(IntKi),OPTIONAL,     INTENT(OUT)   :: Xelement  ! Element identifier
     TYPE(ElemRecType),POINTER,OPTIONAL,INTENT(INOUT)   :: ElemRec ! Return array of elements of kind Xelement


     ErrStat = ErrID_None
     ErrMess = ""

     IF (.not. Mesh%Committed) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshNextElement: Uncommitted mesh."
       RETURN ! Early Return
     ENDIF

     IF ( .NOT. CtrlCode .EQ. MESH_NEXT .AND. (CtrlCode .LT. 0 .OR. CtrlCode .GT. Mesh%nelemlist) ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshNextElement: Invalid CtrlCode"
       RETURN ! Early Return
     ENDIF
     IF ( CtrlCode .EQ. 0 ) THEN              ! 0 means start traversal from beginning
       Mesh%nextelem = 1
     ELSE IF ( CtrlCode .NE. MESH_NEXT ) THEN ! Use the index provided in CtrlCode
       Mesh%nextelem = CtrlCode
     ENDIF
     CtrlCode = 0
     IF ( Mesh%nextelem .GT. Mesh%nelemlist ) THEN
       CtrlCode = MESH_NOMOREELEMS
       RETURN ! Early Return
     ENDIF
     IF ( PRESENT(Ielement) ) Ielement = Mesh%nextelem
     IF ( PRESENT(ElemRec) )  ElemRec => Mesh%elemlist(Mesh%nextelem)%Element
     IF ( PRESENT(Xelement) ) Xelement = Mesh%elemlist(Mesh%nextelem)%Element%Xelement
     Mesh%nextelem = Mesh%nextelem + 1      !bjj should we put this in a modulo statement? (i.e, if we go past the end, nextelem = 1)
     RETURN

   END SUBROUTINE MeshNextElement




   !...............................................................................................................................
    SUBROUTINE MeshExtrapInterp1(u1, u2, tin, u_out, tin_out, ErrStat, ErrMsg )
   !
   ! This subroutine calculates a extrapolated (or interpolated) input u_out at time t_out, from previous/future time
   ! values of u (which has values associated with times in t).  Order of the interpolation is 1
   !
   !  expressions below based on either
   !
   !  f(t) = a
   !  f(t) = a + b * t, or
   !  f(t) = a + b * t + c * t**2
   !
   !  where a, b and c are determined as the solution to
   !  f(t1) = u1, f(t2) = u2, f(t3) = u3  (as appropriate)
   !
   !...............................................................................................................................

    TYPE(MeshType),      INTENT(IN)     :: u1            ! Inputs at t1 > t2
    TYPE(MeshType),      INTENT(IN)     :: u2            ! Inputs at t2
    REAL(DbKi),          INTENT(IN   )  :: tin(:)        ! Times associated with the inputs
    TYPE(MeshType),      INTENT(INOUT)  :: u_out         ! Inputs at tin_out
    REAL(DbKi),          INTENT(IN   )  :: tin_out       ! time to be extrap/interp'd to
    INTEGER(IntKi),      INTENT(  OUT)  :: ErrStat       ! Error status of the operation
    CHARACTER(*),        INTENT(  OUT)  :: ErrMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
    REAL(DbKi)                          :: t(SIZE(tin))  ! Times associated with the inputs
    REAL(DbKi)                          :: t_out         ! Time to which to be extrap/interpd
    INTEGER(IntKi), parameter           :: order = 1     ! order of polynomial fit (max 2)

    REAL(DbKi)                          :: scaleFactor   ! temporary for extrapolation/interpolation

       ! Initialize ErrStat
       ErrStat = ErrID_None
       ErrMsg  = ""

          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)

          ! some error checking:

       if ( size(t) .ne. order+1) then
          ErrStat = ErrID_Fatal
          ErrMsg = ' Error in MeshExtrapInterp1: size(t) must equal 2.'
          RETURN
       end if

      IF ( EqualRealNos( t(1), t(2) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in MeshExtrapInterp1: t(1) must not equal t(2) to avoid a division-by-zero error.'
         RETURN
      END IF

         ! now let's interpolate/extrapolate the fields:
      scaleFactor = t_out / t(2)

      IF ( ALLOCATED(u1%Force) ) THEN
         u_out%Force = u1%Force + (u2%Force - u1%Force) * scaleFactor
      END IF

      IF ( ALLOCATED(u1%Moment) ) THEN
         u_out%Moment = u1%Moment + (u2%Moment - u1%Moment) * scaleFactor
      END IF

      IF ( ALLOCATED(u1%TranslationDisp) ) THEN
         u_out%TranslationDisp = u1%TranslationDisp + (u2%TranslationDisp - u1%TranslationDisp) * scaleFactor
      END IF

      IF ( ALLOCATED(u1%RotationVel) ) THEN
         u_out%RotationVel = u1%RotationVel + (u2%RotationVel - u1%RotationVel) * scaleFactor
      END IF

      IF ( ALLOCATED(u1%TranslationVel) ) THEN
         u_out%TranslationVel = u1%TranslationVel + (u2%TranslationVel - u1%TranslationVel) * scaleFactor
      END IF

      IF ( ALLOCATED(u1%RotationAcc) ) THEN
         u_out%RotationAcc = u1%RotationAcc + (u2%RotationAcc - u1%RotationAcc) * scaleFactor
      END IF

      IF ( ALLOCATED(u1%TranslationAcc) ) THEN
         u_out%TranslationAcc = u1%TranslationAcc + (u2%TranslationAcc - u1%TranslationAcc) * scaleFactor
      END IF

      IF ( ALLOCATED(u1%Scalars) ) THEN
         u_out%Scalars = u1%Scalars + (u2%Scalars - u1%Scalars) * scaleFactor
      END IF

      IF ( ALLOCATED(u1%Orientation) ) THEN
         !ErrStat=ErrID_Info
         !ErrMsg='Orientations not implemented in MeshExtrapInterp1; using nearest neighbor approach instead.'

         IF (scaleFactor < 0.5) THEN
            u_out%Orientation = u1%Orientation
         ELSE
            u_out%Orientation = u2%Orientation
         END IF
      END IF

   END SUBROUTINE MeshExtrapInterp1

   !...............................................................................................................................
    SUBROUTINE MeshExtrapInterp2(u1, u2, u3, tin, u_out, tin_out, ErrStat, ErrMsg )
   !
   ! This subroutine calculates a extrapolated (or interpolated) input u_out at time t_out, from previous/future time
   ! values of u (which has values associated with times in t).  Order of the interpolation is 2.
   !
   !  expressions below based on either
   !
   !  f(t) = a
   !  f(t) = a + b * t, or
   !  f(t) = a + b * t + c * t**2
   !
   !  where a, b and c are determined as the solution to
   !  f(t1) = u1, f(t2) = u2, f(t3) = u3  (as appropriate)
   !
   !...............................................................................................................................

    TYPE(MeshType),      INTENT(IN)     :: u1            ! Inputs at t1 > t2 > t3
    TYPE(MeshType),      INTENT(IN)     :: u2            ! Inputs at t2 > t3
    TYPE(MeshType),      INTENT(IN)     :: u3            ! Inputs at t3
    REAL(DbKi),          INTENT(IN   )  :: tin(:)        ! Times associated with the inputs
    TYPE(MeshType),      INTENT(INOUT)  :: u_out         ! Inputs at tin_out
    REAL(DbKi),          INTENT(IN   )  :: tin_out       ! time to be extrap/interp'd to
    INTEGER(IntKi),      INTENT(  OUT)  :: ErrStat       ! Error status of the operation
    CHARACTER(*),        INTENT(  OUT)  :: ErrMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
    REAL(DbKi)                          :: t(SIZE(tin))  ! Times associated with the inputs
    REAL(DbKi)                          :: t_out         ! Time to which to be extrap/interpd
    INTEGER(IntKi), parameter           :: order = 2     ! order of polynomial fit (max 2)

    REAL(DbKi)                          :: scaleFactor   ! temporary for extrapolation/interpolation

         ! Initialize ErrStat
       ErrStat = ErrID_None
       ErrMsg  = ""

!bjj: TODO check that we've initialized the mesh

          ! we'll subtract a constant from the times to resolve some
          ! numerical issues when t gets large (and to simplify the equations)
       t = tin - tin(1)
       t_out = tin_out - tin(1)


         ! Some error checking:

      if ( size(t) .ne. order+1) then
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error in MeshExtrapInterp2: size(t) must equal 2.'
         RETURN
      end if

      IF ( EqualRealNos( t(1), t(2) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in MeshExtrapInterp2: t(1) must not equal t(2) to avoid a division-by-zero error.'
         RETURN
      END IF
      IF ( EqualRealNos( t(2), t(3) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in MeshExtrapInterp2: t(2) must not equal t(3) to avoid a division-by-zero error.'
         RETURN
      END IF
      IF ( EqualRealNos( t(1), t(3) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in MeshExtrapInterp2: t(1) must not equal t(3) to avoid a division-by-zero error.'
         RETURN
      END IF

         ! Now let's interpolate/extrapolate:

      scaleFactor = t_out / ( t(2) * t(3) * (t(2) - t(3)) )

      IF ( ALLOCATED(u1%Force) ) THEN

         u_out%Force =   u1%Force &
                       + ( t(3)**2 * (u1%Force - u2%Force) + t(2)**2*(-u1%Force + u3%Force) ) * scaleFactor &
                       + ( (t(2)-t(3))*u1%Force + t(3)*u2%Force - t(2)*u3%Force ) *scaleFactor * t_out

      END IF
      IF ( ALLOCATED(u1%Moment) ) THEN
         u_out%Moment =   u1%Moment &
                       + ( t(3)**2 * (u1%Moment - u2%Moment) + t(2)**2*(-u1%Moment + u3%Moment) ) * scaleFactor &
                       + ( (t(2)-t(3))*u1%Moment + t(3)*u2%Moment - t(2)*u3%Moment ) *scaleFactor * t_out
      END IF

      IF ( ALLOCATED(u1%TranslationDisp) ) THEN
         u_out%TranslationDisp =   u1%TranslationDisp &
                               + ( t(3)**2 * ( u1%TranslationDisp - u2%TranslationDisp) &
                                 + t(2)**2 * (-u1%TranslationDisp + u3%TranslationDisp) ) * scaleFactor &
                               + ( (t(2)-t(3))*u1%TranslationDisp + t(3)*u2%TranslationDisp &
                                                                  - t(2)*u3%TranslationDisp )*scaleFactor*t_out
      END IF

      IF ( ALLOCATED(u1%RotationVel) ) THEN
         u_out%RotationVel =   u1%RotationVel &
                           + ( t(3)**2 * ( u1%RotationVel - u2%RotationVel) &
                             + t(2)**2 * (-u1%RotationVel + u3%RotationVel) ) * scaleFactor &
                           + ( (t(2)-t(3))*u1%RotationVel + t(3)*u2%RotationVel - t(2)*u3%RotationVel )*scaleFactor*t_out
      END IF

      IF ( ALLOCATED(u1%TranslationVel) ) THEN
         u_out%TranslationVel =   u1%TranslationVel &
                              +( t(3)**2 * ( u1%TranslationVel - u2%TranslationVel) &
                               + t(2)**2 * (-u1%TranslationVel + u3%TranslationVel) ) * scaleFactor &
                              +( (t(2)-t(3))*u1%TranslationVel + t(3)*u2%TranslationVel - t(2)*u3%TranslationVel)*scaleFactor*t_out
      END IF

      IF ( ALLOCATED(u1%RotationAcc) ) THEN
         u_out%RotationVel =   u1%RotationVel &
                             + ( t(3)**2 * ( u1%RotationVel - u2%RotationVel) &
                               + t(2)**2 * (-u1%RotationVel + u3%RotationVel) ) * scaleFactor &
                            + ( (t(2)-t(3))*u1%RotationVel  + t(3)*u2%RotationVel - t(2)*u3%RotationVel )*scaleFactor*t_out
      END IF

      IF ( ALLOCATED(u1%TranslationAcc) ) THEN
         u_out%TranslationAcc =   u1%TranslationAcc &
                              +( t(3)**2 * ( u1%TranslationAcc - u2%TranslationAcc) &
                               + t(2)**2 * (-u1%TranslationAcc + u3%TranslationAcc) ) * scaleFactor &
                              +( (t(2)-t(3))*u1%TranslationAcc + t(3)*u2%TranslationAcc - t(2)*u3%TranslationAcc)*scaleFactor*t_out
      END IF

      IF ( ALLOCATED(u1%Scalars) ) THEN
         u_out%Scalars =   u1%Scalars &
                       + ( t(3)**2 * (u1%Scalars - u2%Scalars) + t(2)**2*(-u1%Scalars + u3%Scalars) )*scaleFactor &
                       + ( (t(2)-t(3))*u1%Scalars + t(3)*u2%Scalars - t(2)*u3%Scalars )*scaleFactor * t_out
      END IF

      IF ( ALLOCATED(u1%Orientation) ) THEN
         !ErrStat=ErrID_Info
         !ErrMsg=' Orientations are not implemented in MeshExtrapInterp2; using nearest neighbor approach instead.'

         IF ( t_out < 0.5_DbKi*(t(2)+t(1)) ) THEN
            u_out%Orientation = u1%Orientation
         ELSEIF ( t_out < 0.5_DbKi*(t(3)+t(2)) ) THEN
            u_out%Orientation = u2%Orientation
         ELSE
            u_out%Orientation = u3%Orientation
         END IF
      END IF

   END SUBROUTINE MeshExtrapInterp2

   !...............................................................................................................................
    SUBROUTINE MeshExtrapInterp(u, tin, u_out, tin_out, ErrStat, ErrMsg )
   !
   ! This subroutine calculates a extrapolated (or interpolated) input u_out at time t_out, from previous/future time
   ! values of u (which has values associated with times in t).  Order of the interpolation is given by the size of u
   !
   !  expressions below based on either
   !
   !  f(t) = a
   !  f(t) = a + b * t, or
   !  f(t) = a + b * t + c * t**2
   !
   !  where a, b and c are determined as the solution to
   !  f(t1) = u1, f(t2) = u2, f(t3) = u3  (as appropriate)
   !
   !...............................................................................................................................

    TYPE(MeshType),      INTENT(INOUT)  :: u(:)          ! Inputs at t1 > t2 > t3
    REAL(DbKi),          INTENT(IN   )  :: tin(:)        ! Times associated with the inputs
    TYPE(MeshType),      INTENT(INOUT)  :: u_out         ! Inputs at tin_out
    REAL(DbKi),          INTENT(IN   )  :: tin_out       ! time to be extrap/interp'd to
    INTEGER(IntKi),      INTENT(  OUT)  :: ErrStat       ! Error status of the operation
    CHARACTER(*),        INTENT(  OUT)  :: ErrMsg        ! Error message if ErrStat /= ErrID_None

      ! local variables
    REAL(DbKi)                          :: t(SIZE(tin))  ! Times associated with the inputs
    REAL(DbKi)                          :: t_out         ! Time to which to be extrap/interpd
    INTEGER(IntKi)                      :: order         ! order of polynomial fit (max 2)

    REAL(DbKi)                          :: scaleFactor   ! temporary for extrapolation/interpolation

    ! Initialize ErrStat
    ErrStat = ErrID_None
    ErrMsg  = ""

       ! we'll subtract a constant from the times to resolve some
       ! numerical issues when t gets large (and to simplify the equations)
    t = tin - tin(1)
    t_out = tin_out - tin(1)

    if ( size(t) .ne. size(u)) then
       ErrStat = ErrID_Fatal
       ErrMsg = ' Error in MeshExtrapInterp: size(t) must equal size(u).'
       RETURN
    end if

   order = SIZE(u) - 1

   SELECT CASE ( order )
   CASE ( 0 )
      CALL MeshCopy(u(1), u_out, MESH_UPDATECOPY, ErrStat, ErrMsg )
   CASE ( 1 )

      IF ( EqualRealNos( t(1), t(2) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in MeshExtrapInterp: t(1) must not equal t(2) to avoid a division-by-zero error.'
         RETURN
      END IF

      scaleFactor = t_out / t(2)

      IF ( ALLOCATED(u(1)%Force) ) THEN
         u_out%Force = u(1)%Force + (u(2)%Force - u(1)%Force) * scaleFactor
      END IF
      IF ( ALLOCATED(u(1)%Moment) ) THEN
         u_out%Moment = u(1)%Moment + (u(2)%Moment - u(1)%Moment) * scaleFactor
      END IF

      IF ( ALLOCATED(u(1)%TranslationDisp) ) THEN
         u_out%TranslationDisp = u(1)%TranslationDisp + (u(2)%TranslationDisp - u(1)%TranslationDisp) * scaleFactor
      END IF

      IF ( ALLOCATED(u(1)%RotationVel) ) THEN
         u_out%RotationVel = u(1)%RotationVel + (u(2)%RotationVel - u(1)%RotationVel) * scaleFactor
      END IF

      IF ( ALLOCATED(u(1)%TranslationVel) ) THEN
         u_out%TranslationVel = u(1)%TranslationVel + (u(2)%TranslationVel - u(1)%TranslationVel) * scaleFactor
      END IF

      IF ( ALLOCATED(u(1)%RotationAcc) ) THEN
         u_out%RotationAcc = u(1)%RotationAcc + (u(2)%RotationAcc - u(1)%RotationAcc) * scaleFactor
      END IF

      IF ( ALLOCATED(u(1)%TranslationAcc) ) THEN
         u_out%TranslationAcc = u(1)%TranslationAcc + (u(2)%TranslationAcc - u(1)%TranslationAcc) * scaleFactor
      END IF

      IF ( ALLOCATED(u(1)%Scalars) ) THEN
         u_out%Scalars = u(1)%Scalars + (u(2)%Scalars - u(1)%Scalars) * scaleFactor
      END IF

      IF ( ALLOCATED(u(1)%Orientation) ) THEN
         ErrStat=ErrID_Info
         ErrMsg='Orientations not implemented in MeshExtrapInterp; using nearest neighbor approach instead.'

         IF (scaleFactor < 0.5) THEN
            u_out%Orientation = u(1)%Orientation
         ELSE
            u_out%Orientation = u(2)%Orientation
         END IF
      END IF
   CASE ( 2 )

      IF ( EqualRealNos( t(1), t(2) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in MeshExtrapInterp: t(1) must not equal t(2) to avoid a division-by-zero error.'
         RETURN
      END IF
      IF ( EqualRealNos( t(2), t(3) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in MeshExtrapInterp: t(2) must not equal t(3) to avoid a division-by-zero error.'
         RETURN
      END IF
      IF ( EqualRealNos( t(1), t(3) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in MeshExtrapInterp: t(1) must not equal t(3) to avoid a division-by-zero error.'
         RETURN
      END IF

      scaleFactor = 1.0_DbKi / ( t(2) * t(3) * (t(2) - t(3)) )

      IF ( ALLOCATED(u(1)%Force) ) THEN

         u_out%Force =   u(1)%Force &
                       + ( t(3)**2 * (u(1)%Force - u(2)%Force) + t(2)**2*(-u(1)%Force + u(3)%Force) ) * scaleFactor * t_out &
                       + ( (t(2)-t(3))*u(1)%Force + t(3)*u(2)%Force - t(2)*u(3)%Force ) *scaleFactor * t_out**2

      END IF
      IF ( ALLOCATED(u(1)%Moment) ) THEN
         u_out%Moment =   u(1)%Moment &
                       + ( t(3)**2 * (u(1)%Moment - u(2)%Moment) + t(2)**2*(-u(1)%Moment + u(3)%Moment) ) * scaleFactor * t_out &
                       + ( (t(2)-t(3))*u(1)%Moment + t(3)*u(2)%Moment - t(2)*u(3)%Moment ) *scaleFactor * t_out**2
      END IF

      IF ( ALLOCATED(u(1)%TranslationDisp) ) THEN
         u_out%TranslationDisp =   u(1)%TranslationDisp &
                               + ( t(3)**2 * ( u(1)%TranslationDisp - u(2)%TranslationDisp) &
                                 + t(2)**2 * (-u(1)%TranslationDisp + u(3)%TranslationDisp) ) * scaleFactor * t_out &
                               + ( (t(2)-t(3))*u(1)%TranslationDisp + t(3)*u(2)%TranslationDisp &
                                                                    - t(2)*u(3)%TranslationDisp )*scaleFactor*t_out**2
      END IF

      IF ( ALLOCATED(u(1)%RotationVel) ) THEN
         u_out%RotationVel =   u(1)%RotationVel &
                           + ( t(3)**2 * ( u(1)%RotationVel - u(2)%RotationVel) &
                             + t(2)**2 * (-u(1)%RotationVel + u(3)%RotationVel) ) * scaleFactor * t_out &
                           + ( (t(2)-t(3))*u(1)%RotationVel + t(3)*u(2)%RotationVel &
                                                            - t(2)*u(3)%RotationVel )*scaleFactor*t_out**2
      END IF

      IF ( ALLOCATED(u(1)%TranslationVel) ) THEN
         u_out%TranslationVel =   u(1)%TranslationVel &
                              + ( t(3)**2 * ( u(1)%TranslationVel - u(2)%TranslationVel) &
                                + t(2)**2 * (-u(1)%TranslationVel + u(3)%TranslationVel) ) * scaleFactor * t_out &
                              + ( (t(2)-t(3))*u(1)%TranslationVel + t(3)*u(2)%TranslationVel &
                                                                  - t(2)*u(3)%TranslationVel )*scaleFactor*t_out**2
      END IF

      IF ( ALLOCATED(u(1)%RotationAcc) ) THEN
         u_out%RotationVel =   u(1)%RotationVel &
                             + ( t(3)**2 * ( u(1)%RotationVel  - u(2)%RotationVel) &
                                + t(2)**2 * (-u(1)%RotationVel + u(3)%RotationVel) ) * scaleFactor * t_out &
                             + ( (t(2)-t(3))*u(1)%RotationVel  + t(3)*u(2)%RotationVel &
                                                               - t(2)*u(3)%RotationVel )*scaleFactor*t_out**2
      END IF

      IF ( ALLOCATED(u(1)%TranslationAcc) ) THEN
         u_out%TranslationAcc =   u(1)%TranslationAcc &
                              + ( t(3)**2 * ( u(1)%TranslationAcc - u(2)%TranslationAcc) &
                                + t(2)**2 * (-u(1)%TranslationAcc + u(3)%TranslationAcc) ) * scaleFactor * t_out &
                              + ( (t(2)-t(3))*u(1)%TranslationAcc + t(3)*u(2)%TranslationAcc &
                                                                  - t(2)*u(3)%TranslationAcc )*scaleFactor*t_out**2
      END IF

      IF ( ALLOCATED(u(1)%Scalars) ) THEN
         u_out%Scalars =   u(1)%Scalars &
                       + ( t(3)**2 * (u(1)%Scalars - u(2)%Scalars) + t(2)**2*(-u(1)%Scalars + u(3)%Scalars) )*scaleFactor * t_out &
                       + ( (t(2)-t(3))*u(1)%Scalars + t(3)*u(2)%Scalars - t(2)*u(3)%Scalars )*scaleFactor * t_out**2
      END IF

      IF ( ALLOCATED(u(1)%Orientation) ) THEN
         ErrStat=ErrID_Info
         ErrMsg=' Orientations are not implemented in MeshExtrapInterp; using nearest neighbor approach instead.'

         IF ( t_out < 0.5_DbKi*(t(2)+t(1)) ) THEN
            u_out%Orientation = u(1)%Orientation
         ELSEIF ( t_out < 0.5_DbKi*(t(3)+t(2)) ) THEN
            u_out%Orientation = u(2)%Orientation
         ELSE
            u_out%Orientation = u(3)%Orientation
         END IF

      END IF
   CASE DEFAULT
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error in MeshExtrapInterp: size(u) must be less than 4 (order must be less than 3).'
   END SELECT


 END SUBROUTINE MeshExtrapInterp

END MODULE ModMesh


