!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of the NWTC Subroutine Library.
!
!    ServoDyn is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with ServoDyn.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
MODULE ModMesh

!======================================================================================================================
! WARNING:  Because someone used preprocessor directives to comment out code instead of temporary comment characters,
!           and because some lines go beyond column 132, you must specify the following compiler options:
!
!              Intel:   /fpp
!              Gnu:     -x f95-cpp-input -ffree-line-length-none
!======================================================================================================================

   USE ModMesh_Types
   IMPLICIT NONE
!   PRIVATE

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

   SUBROUTINE MeshPrintInfo ( U, M, N)
     INTEGER, INTENT(IN   )                ::      U  ! fortran output unit
     TYPE(MeshType),TARGET,INTENT(INOUT)   ::      M  ! mesh to be reported on
     INTEGER, OPTIONAL,INTENT(IN   )       ::      N  ! Number to print, default 5
    ! Local
     INTEGER isz,i,j,nn,CtrlCode,Ielement,Xelement
     INTEGER                    :: ErrStat
     CHARACTER(256)             :: ErrMess

     nn = 5
     IF (PRESENT(N)) nn = N

     write(U,*)  'Initialized: ', M%initialized
     write(U,*)  'Fieldmask:   ', M%FieldMask
     IF ( M%FieldMask( MASKID_FORCE           ) )  write(U,*)  '  Defined : Force'
     IF ( M%FieldMask( MASKID_MOMENT          ) )  write(U,*)  '  Defined : Moment'
     IF ( M%FieldMask( MASKID_ORIENTATION     ) )  write(U,*)  '  Defined : Orientation'
     IF ( M%FieldMask( MASKID_TRANSLATIONDISP ) )  write(U,*)  '  Defined : TranslationDisp'
     IF ( M%FieldMask( MASKID_TRANSLATIONVEL  ) )  write(U,*)  '  Defined : TranslationVel'
     IF ( M%FieldMask( MASKID_ROTATIONVEL     ) )  write(U,*)  '  Defined : RotationVel'
     IF ( M%FieldMask( MASKID_TRANSLATIONACC  ) )  write(U,*)  '  Defined : TranslationAcc'
     IF ( M%FieldMask( MASKID_ROTATIONACC     ) )  write(U,*)  '  Defined : RotationAcc'
     IF ( M%FieldMask( MASKID_ADDEDMASS       ) )  write(U,*)  '  Defined : AddedMass'
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
       write(U,*)'Position: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%Position(:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%Force))THEN
       isz=size(M%Force,2)
       write(U,*)'Force: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%Force(:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%Moment))THEN
       isz=size(M%Moment,2)
       write(U,*)'Moment: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%Moment(:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%Orientation))THEN
       isz=size(M%Orientation,3)
       write(U,*)'Orientation: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%Orientation(:,:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%TranslationDisp))THEN
       isz=size(M%TranslationDisp,2)
       write(U,*)'TranslationDisp: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%TranslationDisp(:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%RotationVel))THEN
       isz=size(M%RotationVel,2)
       write(U,*)'RotationVel: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%RotationVel(:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%TranslationVel))THEN
       isz=size(M%TranslationVel,2)
       write(U,*)'TranslationVel: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%TranslationVel(:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%RotationAcc))THEN
       isz=size(M%RotationAcc,2)
       write(U,*)'RotationAcc: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%RotationAcc(:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%TranslationAcc))THEN
       isz=size(M%TranslationAcc,2)
       write(U,*)'TranslationAcc: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%TranslationAcc(:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%AddedMass))THEN
       isz=size(M%AddedMass,3)
       write(U,*)'AddedMass: ',isz
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%AddedMass(:,:,i)
       ENDDO
     ENDIF
     IF(M%initialized.AND.ASSOCIATED(M%Scalars))THEN
       isz=size(M%Scalars,1)
       write(U,*)'Scalars: ',isz
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
        RETURN
      END IF

     DO WHILE ( CtrlCode .NE. MESH_NOMORE )
       WRITE(U,*)'  Ielement: ', Ielement,' ',ElemNames(Xelement)
       CtrlCode = MESH_NEXT
       CALL MeshNextElement( M, CtrlCode, ErrStat, ErrMess, Ielement=Ielement, Xelement=Xelement )
        IF (ErrStat >= AbortErrLev) THEN
           WRITE(U,*) ' Error in MeshNextElement(): '
           WRITE(U,*) TRIM(ErrMess)
           RETURN
         END IF

     ENDDO

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
                          ,AddedMass                                                       &
                          ,nScalars                                                        &
                         )

     TYPE(MeshType),TARGET,INTENT(INOUT)   :: BlankMesh ! Mesh to be created
     INTEGER,INTENT(IN)         :: IOS                  ! input (1), output(2), or state(3)
     INTEGER,INTENT(IN)         :: Nnodes               ! Number of nodes in mesh
     INTEGER,INTENT(INOUT)      :: ErrStat              ! error status/level
     CHARACTER(*),INTENT(INOUT) :: ErrMess              ! error message
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
     LOGICAL,OPTIONAL,INTENT(IN):: AddedMass            ! If present and true, allocate AddedMass field
!
     INTEGER,OPTIONAL,INTENT(IN):: nScalars             ! If present > 0, allocate n Scalars

    ! Local
     INTEGER i

     ErrStat = ErrID_None
     ErrMess = ""

     IF ( mesh_debug ) print*,'Called MeshCreate'


#if 0
! This is problematic called from MeshCreate for a new mesh because
! we can't assume that uninitialized pointers are uninitialized
! Just put it on the user not to recreate meshes that already exist (memory leak danger if they do)


!bjj: I found one pointer wasn't initialized with => NULL
! I fixed the default initializion and now tested this code on both IVF
! and gfortran (albeit a small case and only on Windows). It works for me.


     CALL MeshDestroy( BlankMesh, ErrStat, ErrMess, .TRUE. )
                                                        ! make sure we're not leaving any pointers dangling
                                                        ! and nullify them for good measure
                                                        ! See comment on optional IgnoreSibling argument
                                                        ! in definition of MeshDestroy
      IF (ErrStat >= AbortErrLev) RETURN
#endif

     BlankMesh%initialized = .TRUE.
     BlankMesh%IOS         = IOS

!bjj: check that IOS is valid and Nnodes > 0?

     BlankMesh%Nnodes = Nnodes
     NULLIFY( BlankMesh%ElemTable )
     NULLIFY( BlankMesh%SiblingMesh )
     NULLIFY( BlankMesh%ElemList )
     BlankMesh%nelemlist = 0 ; BlankMesh%maxelemlist = 0 ;

     CALL AllocPAry( BlankMesh%Position, 3, Nnodes, 'CreateMesh: Position' )

     ALLOCATE(BlankMesh%ElemTable(NELEMKINDS))
     DO i = 1, NELEMKINDS
       BlankMesh%ElemTable(i)%nelem = 0  ; BlankMesh%ElemTable(i)%maxelem = 0
       NULLIFY(BlankMesh%ElemTable(i)%Elements )
     ENDDO

   ! handle optionals
     BlankMesh%FieldMask = .FALSE.

     NULLIFY(BlankMesh%Force)
     IF ( PRESENT(Force) ) THEN
       IF ( Force ) THEN
         CALL AllocPAry( BlankMesh%Force, 3, Nnodes, 'CreateMesh: Force' )
         BlankMesh%FieldMask(MASKID_FORCE) = .TRUE.
       ENDIF
     ENDIF
     NULLIFY(BlankMesh%Moment)
     IF ( PRESENT(Moment) ) THEN
       IF ( Moment ) THEN
         CALL AllocPAry( BlankMesh%Moment, 3, Nnodes, 'CreateMesh: Moment' )
         BlankMesh%FieldMask(MASKID_MOMENT) = .TRUE.
       ENDIF
     ENDIF
     NULLIFY(BlankMesh%Orientation)
     IF ( PRESENT(Orientation) ) THEN
       IF ( Orientation ) THEN
         CALL AllocPAry( BlankMesh%Orientation, 3, 3, Nnodes, 'CreateMesh: Orientation' )
         BlankMesh%FieldMask(MASKID_ORIENTATION) = .TRUE.
       ENDIF
     ENDIF
     NULLIFY(BlankMesh%TranslationDisp)
     IF ( PRESENT(TranslationDisp) ) THEN
       IF ( TranslationDisp ) THEN
         CALL AllocPAry( BlankMesh%TranslationDisp, 3, Nnodes, 'CreateMesh: TranslationDisp' )
         BlankMesh%FieldMask(MASKID_TRANSLATIONDISP) = .TRUE.
       ENDIF
     ENDIF
     NULLIFY(BlankMesh%TranslationVel)
     IF ( PRESENT(TranslationVel) ) THEN
       IF ( TranslationVel ) THEN
         CALL AllocPAry( BlankMesh%TranslationVel, 3, Nnodes, 'CreateMesh: TranslationVel' )
         BlankMesh%FieldMask(MASKID_TRANSLATIONVEL) = .TRUE.
       ENDIF
     ENDIF
     NULLIFY(BlankMesh%RotationVel)
     IF ( PRESENT(RotationVel) ) THEN
       IF ( RotationVel ) THEN
         CALL AllocPAry( BlankMesh%RotationVel, 3, Nnodes, 'CreateMesh: RotationVel' )
         BlankMesh%FieldMask(MASKID_ROTATIONVEL) = .TRUE.
       ENDIF
     ENDIF
     NULLIFY(BlankMesh%TranslationAcc)
     IF ( PRESENT(TranslationAcc) ) THEN
       IF ( TranslationAcc ) THEN
         CALL AllocPAry( BlankMesh%TranslationAcc, 3, Nnodes, 'CreateMesh: TranslationAcc' )
         BlankMesh%FieldMask(MASKID_TRANSLATIONACC) = .TRUE.
       ENDIF
     ENDIF
     NULLIFY(BlankMesh%RotationAcc)
     IF ( PRESENT(RotationAcc) ) THEN
       IF ( RotationAcc ) THEN
         CALL AllocPAry( BlankMesh%RotationAcc, 3, Nnodes, 'CreateMesh: RotationAcc' )
         BlankMesh%FieldMask(MASKID_ROTATIONACC) = .TRUE.
       ENDIF
     ENDIF
     NULLIFY(BlankMesh%AddedMass)
     IF ( PRESENT(AddedMass) ) THEN
       IF ( AddedMass ) THEN
         CALL AllocPAry( BlankMesh%AddedMass, 6, 6, Nnodes, 'CreateMesh: AddedMass' )
         BlankMesh%FieldMask(MASKID_AddedMass) = .TRUE.
       ENDIF
     ENDIF
     NULLIFY(BlankMesh%Scalars)
     BlankMesh%nScalars = 0
     IF ( PRESENT(nScalars) ) THEN
       IF ( nScalars .GT. 0 ) THEN
         CALL AllocPAry( BlankMesh%Scalars, nScalars, Nnodes, 'CreateMesh: Scalars' )
         BlankMesh%FieldMask(MASKID_Scalar) = .TRUE.
         BlankMesh%nScalars = nScalars
       ENDIF
     ENDIF

     RETURN

   END SUBROUTINE MeshCreate

   RECURSIVE SUBROUTINE MeshDestroy ( Mesh, ErrStat, ErrMess, IgnoreSibling )
     TYPE(MeshType),INTENT(INOUT)   :: Mesh            ! Mesh to be vaporized
     INTEGER, INTENT(OUT)  :: ErrStat
     CHARACTER(*), INTENT(OUT)  :: ErrMess
    ! On a brand new mesh, the pointers to siblings may not be nullified and
    ! thus undefined (which may cause ASSOCIATED to report .true. erroneously)
    ! This despite use of => NULL in declaration of this fields for MeshType. Sigh.
    ! So if IgnoreSibling is present and true, don't follow the sibling pointers.
    ! Instead just unconditionally nullify these.
    ! Use this carefully, since it can leave dangling memory if used for a
    ! mesh that already exists and has existing siblings.
     LOGICAL, OPTIONAL :: IgnoreSibling

    ! Local
     LOGICAL IgSib
     INTEGER i, j, k

     ErrStat = ErrID_None

     IF ( .NOT. Mesh%Initialized ) RETURN

     Mesh%initialized = .FALSE.
     Mesh%fieldmask   = .FALSE.
     Mesh%ios         = 0
     Mesh%Nnodes      = 0

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
         ENDIF
       ENDDO
       DEALLOCATE(Mesh%ElemTable)
       NULLIFY(Mesh%ElemTable)
     ENDIF

#if 0
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
     IF (ASSOCIATED(Mesh%Position))          NULLIFY(Mesh%Position)
     IF (ASSOCIATED(Mesh%Force))             NULLIFY(Mesh%Force)
     IF (ASSOCIATED(Mesh%Moment))            NULLIFY(Mesh%Moment)
     IF (ASSOCIATED(Mesh%Orientation))       NULLIFY(Mesh%Orientation)
     IF (ASSOCIATED(Mesh%TranslationDisp))    NULLIFY(Mesh%TranslationDisp)
     IF (ASSOCIATED(Mesh%RotationVel))       NULLIFY(Mesh%RotationVel)
     IF (ASSOCIATED(Mesh%TranslationVel))    NULLIFY(Mesh%TranslationVel)
     IF (ASSOCIATED(Mesh%RotationAcc))       NULLIFY(Mesh%RotationAcc)
     IF (ASSOCIATED(Mesh%TranslationAcc))    NULLIFY(Mesh%TranslationAcc)
     IF (ASSOCIATED(Mesh%Scalars))           NULLIFY(Mesh%Scalars)
#endif

     IgSib = .FALSE.
     IF ( PRESENT( IgnoreSibling ) ) THEN
       IgSib = IgnoreSibling
     ENDIF
     IF ( .NOT. IgSib ) THEN
       IF ( ASSOCIATED( Mesh%SiblingMesh ) ) THEN
         ! Nullify the back pointer to avoid an endless recursion
         IF ( ASSOCIATED( Mesh%SiblingMesh%SiblingMesh )) THEN
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
         NULLIFY( Mesh%SiblingMesh%ElemTable )
         CALL MeshDestroy( Mesh%SiblingMesh, ErrStat, ErrMess )
         IF (ErrStat >= AbortErrLev) RETURN
       ENDIF
     ENDIF
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
     REAL(ReKi), ALLOCATABLE, INTENT(OUT)  :: ReBuf(:)           ! Real buffer
     REAL(DbKi), ALLOCATABLE, INTENT(OUT)  :: DbBuf(:)           ! Double buffer
     INTEGER(IntKi), ALLOCATABLE, INTENT(OUT) :: IntBuf(:)       ! Int buffer
     INTEGER(IntKi),   INTENT(  OUT) :: ErrStat
     CHARACTER(*),     INTENT(  OUT) :: ErrMess
     LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly
   ! Local
     INTEGER i,ic,nelem,n_int,n_re,n_db,l,ii,jj,CtrlCode,x
     INTEGER Ielement, Xelement
     TYPE(ElemRecType), POINTER :: ElemRec
     LOGICAL SzOnly

     SzOnly = .FALSE.
     IF ( PRESENT(SizeOnly) ) SzOnly = SizeOnly

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
    ! Position
     n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_FORCE) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_MOMENT) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_ORIENTATION) ) n_re = n_re + Mesh%Nnodes * 9
     IF ( Mesh%FieldMask(MASKID_ROTATIONVEL) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_ROTATIONACC) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_TRANSLATIONDISP) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_TRANSLATIONVEL) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_TRANSLATIONACC) ) n_re = n_re + Mesh%Nnodes * 3
     IF ( Mesh%FieldMask(MASKID_ADDEDMASS) ) n_re = n_re + Mesh%Nnodes * 36
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
       IF ( Mesh%FieldMask(MASKID_ADDEDMASS) ) THEN ! n_re = n_re + Mesh%Nnodes * 36
         DO i = 1, Mesh%Nnodes
           DO jj = 1,6
             DO ii = 1,6
               ReBuf(ic) = Mesh%AddedMass(ii,jj,i) ; ic = ic + 1
             ENDDO
           ENDDO
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
     REAL(ReKi), ALLOCATABLE,     INTENT(IN   ) :: Re_Buf(:)
     REAL(DbKi), ALLOCATABLE,     INTENT(IN   ) :: Db_Buf(:)
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
     AddedMass = .FALSE.
     nScalars = 0
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_FORCE-1)       .EQ. 1 ) Force = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_MOMENT-1)      .EQ. 1 ) Moment = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_ORIENTATION-1) .EQ. 1 ) Orientation = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_TRANSLATIONDISP-1) .EQ. 1 ) TranslationDisp = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_TRANSLATIONVEL-1) .EQ. 1 ) TranslationVel = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_ROTATIONVEL-1)    .EQ. 1 ) RotationVel = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_TRANSLATIONACC-1) .EQ. 1 ) TranslationAcc = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_ROTATIONACC-1)    .EQ. 1 ) RotationAcc = .TRUE.
     IF ( Int_Buf(HDR_FIELDMASK+MASKID_ADDEDMASS-1)   .EQ. 1 ) AddedMass = .TRUE.
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
                     ,AddedMass=AddedMass, nScalars = nScalars                         &
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
     IF ( Mesh%FieldMask(MASKID_ADDEDMASS) ) THEN
       DO i = 1, Mesh%Nnodes
         DO jj = 1,6
           DO ii = 1,6
             Mesh%AddedMass(ii,jj,i) = Re_Buf(ic) ; ic = ic + 1
           ENDDO
         ENDDO
       ENDDO
     ENDIF
     IF ( Mesh%nScalars .GT. 0 ) THEN
       DO i = 1, Mesh%Nnodes
         DO ii = 1,Mesh%nScalars
           Mesh%Scalars(ii,i) = Re_Buf(ic) ; ic = ic + 1
         ENDDO
       ENDDO
     ENDIF
     RETURN

   END SUBROUTINE MeshUnpack

   SUBROUTINE MeshCopy( SrcMesh, DestMesh, CtrlCode, ErrStat , ErrMess   &
                      , Force, Moment, Orientation, TranslationDisp, TranslationVel, RotationVel, TranslationAcc, RotationAcc, AddedMass, nScalars )
     TYPE(MeshType), TARGET,      INTENT(INOUT) :: SrcMesh  ! Mesh being copied
     TYPE(MeshType), TARGET,      INTENT(INOUT) :: DestMesh ! Copy of mesh
     INTEGER(IntKi),              INTENT(IN)    :: CtrlCode ! MESH_NEWCOPY, MESH_SIBLING, or
                                                            ! MESH_UPDATECOPY
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat  ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess  ! Error message
    ! Optional arguments (used only if CtrlCode is MESH_SIBLING):
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: Force             & ! If present and true, allocate Force field
                                                 , Moment            & ! If present and true, allocate Moment field
                                                 , Orientation       & ! If present and true, allocate Orientation field
                                                 , TranslationDisp   & ! If present and true, allocate TranslationDisp field
                                                 , TranslationVel    & ! If present and true, allocate TranslationVel field
                                                 , RotationVel       & ! If present and true, allocate RotationVel field
                                                 , TranslationAcc    & ! If present and true, allocate TranslationAcc field
                                                 , RotationAcc       & ! If present and true, allocate RotationAcc field
                                                 , AddedMass           ! If present and true, allocate AddedMess field
     INTEGER(IntKi), OPTIONAL,    INTENT(IN)    :: nScalars            ! If present and > 0 , alloc n Scalars
    ! Local
     LOGICAL                                    :: Force_l           & ! If present and true, allocate Force field
                                                 , Moment_l          & ! If present and true, allocate Moment field
                                                 , Orientation_l     & ! If present and true, allocate Orientation field
                                                 , TranslationDisp_l & ! If present and true, allocate TranslationDisp field
                                                 , TranslationVel_l  & ! If present and true, allocate TranslationVel field
                                                 , RotationVel_l     & ! If present and true, allocate RotationVel field
                                                 , TranslationAcc_l  & ! If present and true, allocate TranslationAcc field
                                                 , RotationAcc_l     & ! If present and true, allocate RotationAcc field
                                                 , AddedMass_l         ! If present and true, allocate AddedMess field
     INTEGER(IntKi)                             :: nScalars_l          ! If present and true, alloc n Scalars
     INTEGER i, j, k


     ErrStat = ErrID_None

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
                            ,AddedMass=SrcMesh%FieldMask(MASKID_ADDEDMASS)                                      &
                            ,nScalars=SrcMesh%nScalars                                                          )

            IF ( ASSOCIATED(SrcMesh%Force ) ) THEN
               DestMesh%Force = SrcMesh%Force
            ENDIF
            IF ( ASSOCIATED(SrcMesh%Moment ) ) THEN
              DestMesh%Moment = SrcMesh%Moment
            ENDIF
            IF ( ASSOCIATED(SrcMesh%Orientation ) ) THEN
              DestMesh%Orientation = SrcMesh%Orientation
            ENDIF
            IF ( ASSOCIATED(SrcMesh%TranslationDisp ) ) THEN
              DestMesh%TranslationDisp = SrcMesh%TranslationDisp
            ENDIF
            IF ( ASSOCIATED(SrcMesh%TranslationVel ) ) THEN
              DestMesh%TranslationVel = SrcMesh%TranslationVel
            ENDIF
            IF ( ASSOCIATED(SrcMesh%RotationVel ) ) THEN
              DestMesh%RotationVel = SrcMesh%RotationVel
            ENDIF
            IF ( ASSOCIATED(SrcMesh%TranslationAcc ) ) THEN
              DestMesh%TranslationAcc = SrcMesh%TranslationAcc
            ENDIF
            IF ( ASSOCIATED(SrcMesh%RotationAcc ) ) THEN
              DestMesh%RotationAcc = SrcMesh%RotationAcc
            ENDIF
            IF ( ASSOCIATED(SrcMesh%AddedMass ) ) THEN
              DestMesh%AddedMass = SrcMesh%AddedMass
            ENDIF
            IF ( ASSOCIATED(SrcMesh%Scalars ) ) THEN
              DestMesh%Scalars = SrcMesh%Scalars
            ENDIF

            DestMesh%Position =  SrcMesh%Position


            ALLOCATE(DestMesh%ElemTable(size(SrcMesh%ElemTable)))
            DO i = 1, NELEMKINDS
               DestMesh%ElemTable(i)%nelem = SrcMesh%ElemTable(i)%nelem
               DestMesh%ElemTable(i)%maxelem = SrcMesh%ElemTable(i)%maxelem
               ALLOCATE(DestMesh%ElemTable(i)%Elements(DestMesh%ElemTable(i)%maxelem))
               DO j = 1,DestMesh%ElemTable(i)%nelem
                  ALLOCATE(DestMesh%ElemTable(i)%Elements(j)%ElemNodes(size(SrcMesh%ElemTable(i)%Elements(j)%ElemNodes)))
                  DestMesh%ElemTable(i)%Elements(j)%ElemNodes = SrcMesh%ElemTable(i)%Elements(j)%ElemNodes
               ENDDO
            ENDDO

!bjj: aren't we also copying the other fields:

               ! Regenerate new list of elements (point to ElemTable)

            DestMesh%nelemlist   = SrcMesh%nelemlist
            DestMesh%maxelemlist = SrcMesh%maxelemlist
            ALLOCATE(DestMesh%ElemList(DestMesh%maxelemlist))

            k = 0
            DO i = 1, NELEMKINDS
               DO j = 1, DestMesh%ElemTable(i)%nelem
                  k = k + 1
                  DestMesh%elemlist(k)%Element => DestMesh%ElemTable(i)%Elements(j)
                  DestMesh%elemlist(k)%Element%Xelement = i
               ENDDO
            END DO

            DestMesh%RemapFlag   = SrcMesh%RemapFlag

         ELSE IF ( CtrlCode .EQ. MESH_SIBLING ) THEN
            Force_l            = .FALSE. ; IF ( PRESENT(Force) )                     Force_l = Force
            Moment_l           = .FALSE. ; IF ( PRESENT(Moment) )                   Moment_l = Moment
            Orientation_l      = .FALSE. ; IF ( PRESENT(Orientation) )         Orientation_l = Orientation
            TranslationDisp_l  = .FALSE. ; IF ( PRESENT(TranslationDisp) ) TranslationDisp_l = TranslationDisp
            TranslationVel_l   = .FALSE. ; IF ( PRESENT(TranslationVel) )   TranslationVel_l = TranslationVel
            RotationVel_l      = .FALSE. ; IF ( PRESENT(RotationVel) )         RotationVel_l = RotationVel
            TranslationAcc_l   = .FALSE. ; IF ( PRESENT(TranslationAcc) )   TranslationAcc_l = TranslationAcc
            RotationAcc_l      = .FALSE. ; IF ( PRESENT(RotationAcc) )         RotationAcc_l = RotationAcc
            AddedMass_l        = .FALSE. ; IF ( PRESENT(AddedMass) )             AddedMass_l = AddedMass
            nScalars_l         = 0       ; IF ( PRESENT(nScalars) )               nScalars_l = nScalars
            CALL MeshCreate( DestMesh, IOS=SrcMesh%IOS, Nnodes=SrcMesh%Nnodes, ErrStat=ErrStat, ErrMess=ErrMess &
                            ,Force=Force_l                                                                      &
                            ,Moment=Moment_l                                                                    &
                            ,Orientation=Orientation_l                                                          &
                            ,TranslationDisp=TranslationDisp_l                                                  &
                            ,TranslationVel=TranslationVel_l                                                    &
                            ,RotationVel=RotationVel_l                                                          &
                            ,TranslationAcc=TranslationAcc_l                                                    &
                            ,RotationAcc=RotationAcc_l                                                          &
                            ,AddedMass=AddedMass_l                                                              &
                            ,nScalars=nScalars_l                                                                )
            DestMesh%SiblingMesh => SrcMesh
            SrcMesh%SiblingMesh => DestMesh

            DestMesh%Position => SrcMesh%Position
            DestMesh%ElemTable => SrcMesh%ElemTable

            !bjj: what about this one:
            DestMesh%ElemList => SrcMesh%ElemList

            !bjj: and these?
            DestMesh%nelemlist   = SrcMesh%nelemlist
            DestMesh%maxelemlist = SrcMesh%maxelemlist
            DestMesh%nextelem    = SrcMesh%nextelem


            !bjj: this isn't a pointer, but possibly should be
            DestMesh%RemapFlag   = SrcMesh%RemapFlag
            !DestMesh%RemapFlag => SrcMesh%RemapFlag

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
      !         CALL AllocPAry( DestMesh%element_point, 1, SrcMesh%maxpoint, 'CreateMesh: element_point' )
               DestMesh%element_point => SrcMesh%element_point
             ENDIF
             DestMesh%nline2   = SrcMesh%nline2 ; DestMesh%maxline2 = SrcMesh%maxline2
             IF ( SrcMesh%nline2 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_line2 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_line2, 2, SrcMesh%maxline2, 'CreateMesh: element_line2' )
               DestMesh%element_line2 => SrcMesh%element_line2
             ENDIF
             DestMesh%nline3   = SrcMesh%nline3 ; DestMesh%maxline3 = SrcMesh%maxline3
             IF ( SrcMesh%nline3 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_line3 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_line3, 3, SrcMesh%maxline3, 'CreateMesh: element_line3' )
               DestMesh%element_line3 => SrcMesh%element_line3
             ENDIF
             DestMesh%ntri3   = SrcMesh%ntri3 ; DestMesh%maxtri3 = SrcMesh%maxtri3
             IF ( SrcMesh%ntri3 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_tri3 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_tri3, 3, SrcMesh%maxtri3, 'CreateMesh: element_tri3' )
               DestMesh%element_tri3 => SrcMesh%element_tri3
             ENDIF
             DestMesh%ntri6   = SrcMesh%ntri6 ; DestMesh%maxtri6 = SrcMesh%maxtri6
             IF ( SrcMesh%ntri6 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_tri6 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_tri6, 6, SrcMesh%maxtri6, 'CreateMesh: element_tri6' )
               DestMesh%element_tri6 => SrcMesh%element_tri6
             ENDIF
             DestMesh%nquad4   = SrcMesh%nquad4 ; DestMesh%maxquad4 = SrcMesh%maxquad4
             IF ( SrcMesh%nquad4 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_quad4 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_quad4, 4, SrcMesh%maxquad4, 'CreateMesh: element_quad4' )
               DestMesh%element_quad4 => SrcMesh%element_quad4
             ENDIF
             DestMesh%nquad8   = SrcMesh%nquad8 ; DestMesh%maxquad8 = SrcMesh%maxquad8
             IF ( SrcMesh%nquad8 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_quad8 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_quad8, 8, SrcMesh%maxquad8, 'CreateMesh: element_quad8' )
               DestMesh%element_quad8 => SrcMesh%element_quad8
             ENDIF
             DestMesh%ntet4   = SrcMesh%ntet4 ; DestMesh%maxtet4 = SrcMesh%maxtet4
             IF ( SrcMesh%ntet4 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_tet4 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_tet4, 4, SrcMesh%maxtet4, 'CreateMesh: element_tet4' )
               DestMesh%element_tet4 => SrcMesh%element_tet4
             ENDIF
             DestMesh%ntet10   = SrcMesh%ntet10 ; DestMesh%maxtet10 = SrcMesh%maxtet10
             IF ( SrcMesh%ntet10 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_tet10 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_tet10, 10, SrcMesh%maxtet10, 'CreateMesh: element_tet10' )
               DestMesh%element_tet10 = SrcMesh%element_tet10
             ENDIF
             DestMesh%nhex8   = SrcMesh%nhex8 ; DestMesh%maxhex8 = SrcMesh%maxhex8
             IF ( SrcMesh%nhex8 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_hex8 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_hex8, 8, SrcMesh%maxhex8, 'CreateMesh: element_hex8' )
               DestMesh%element_hex8 => SrcMesh%element_hex8
             ENDIF
             DestMesh%nhex20   = SrcMesh%nhex20 ; DestMesh%maxhex20 = SrcMesh%maxhex20
             IF ( SrcMesh%nhex20 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_hex20 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_hex20, 20, SrcMesh%maxhex20, 'CreateMesh: element_hex20' )
               DestMesh%element_hex20 => SrcMesh%element_hex20
             ENDIF
             DestMesh%nwedge6   = SrcMesh%nwedge6 ; DestMesh%maxwedge6 = SrcMesh%maxwedge6
             IF ( SrcMesh%nwedge6 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_wedge6 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_wedge6, 6, SrcMesh%maxwedge6, 'CreateMesh: element_wedge6' )
               DestMesh%element_wedge6 => SrcMesh%element_wedge6
             ENDIF
             DestMesh%nwedge15   = SrcMesh%nwedge15 ; DestMesh%maxwedge15 = SrcMesh%maxwedge15
             IF ( SrcMesh%nwedge15 .GT. 0 .AND. ASSOCIATED( SrcMesh%element_wedge15 ) ) THEN
      !         CALL AllocPAry( DestMesh%element_wedge15, 15, SrcMesh%maxwedge15, 'CreateMesh: element_wedge15' )
               DestMesh%element_wedge15 => SrcMesh%element_wedge15
             ENDIF
          ENDIF
#endif
      ELSE IF ( CtrlCode .EQ. MESH_UPDATECOPY ) THEN
         IF ( SrcMesh%nNodes .NE. DestMesh%nNodes ) THEN
            ErrStat = ErrID_Fatal
            ErrMess = "MeshCopy. MESH_UPDATECOPY of meshes with different numbers of nodes."
         ELSE
            DestMesh%Position    = SrcMesh%Position
            IF ( ASSOCIATED(SrcMesh%Force ) ) THEN
               DestMesh%Force = SrcMesh%Force
            ENDIF
            IF ( ASSOCIATED(SrcMesh%Moment ) ) THEN
               DestMesh%Moment = SrcMesh%Moment
            ENDIF
            IF ( ASSOCIATED(SrcMesh%Orientation ) ) THEN
               DestMesh%Orientation = SrcMesh%Orientation
            ENDIF
            IF ( ASSOCIATED(SrcMesh%TranslationDisp ) ) THEN
               DestMesh%TranslationDisp = SrcMesh%TranslationDisp
            ENDIF
            IF ( ASSOCIATED(SrcMesh%TranslationVel ) ) THEN
               DestMesh%TranslationVel = SrcMesh%TranslationVel
            ENDIF
            IF ( ASSOCIATED(SrcMesh%RotationVel ) ) THEN
               DestMesh%RotationVel = SrcMesh%RotationVel
            ENDIF
            IF ( ASSOCIATED(SrcMesh%TranslationAcc ) ) THEN
               DestMesh%TranslationAcc = SrcMesh%TranslationAcc
            ENDIF
            IF ( ASSOCIATED(SrcMesh%RotationAcc ) ) THEN
               DestMesh%RotationAcc = SrcMesh%RotationAcc
            ENDIF
            IF ( ASSOCIATED(SrcMesh%AddedMass ) ) THEN
               DestMesh%AddedMass = SrcMesh%AddedMass
            ENDIF
            IF ( ASSOCIATED(SrcMesh%Scalars ) ) THEN
               DestMesh%Scalars = SrcMesh%Scalars
            ENDIF


            !bjj: should we also update the Remapflag and nextelem fields?
            DestMesh%RemapFlag = SrcMesh%RemapFlag
            DestMesh%nextelem  = SrcMesh%nextelem

         ENDIF

      ENDIF

      DestMesh%Initialized = SrcMesh%Initialized
      !DestMesh%spatial = SrcMesh%spatial !bjj: unused?

   END SUBROUTINE MeshCopy

! added 20130102 as stub for AeroDyn work
   SUBROUTINE MeshPositionNode( Mesh, Inode, Pos, ErrStat, ErrMess )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being spatio-located
     INTEGER(IntKi),              INTENT(IN   ) :: Inode     ! Number of node being located
     REAL(ReKi),                  INTENT(IN   ) :: Pos(3)    ! Xi,Yi,Zi, coordinates of node
     INTEGER(IntKi),              INTENT(  OUT) :: ErrStat   ! Error code
     CHARACTER(*),                INTENT(  OUT) :: ErrMess   ! Error message

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
     IF ( ErrStat .NE. ErrID_None ) RETURN   ! early return on error
    ! Business
     Mesh%Position(:,Inode) = Pos
     RETURN

   END SUBROUTINE MeshPositionNode

   SUBROUTINE MeshCommit( Mesh, ErrStat, ErrMess )
     TYPE(MeshType),              INTENT(INOUT) :: Mesh ! Mesh being committed
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat ! Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess ! Error message
    ! Local
     INTEGER n0d, n1d, n2d, n3d, nn, i, j
     TYPE(ElemListType), POINTER :: tmp(:)

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


     ! Construct list of elements

     DO i = 1, NELEMKINDS
       DO j = 1, Mesh%ElemTable(i)%nelem
         Mesh%Nelemlist = Mesh%Nelemlist + 1
         IF ( Mesh%Nelemlist .GE. Mesh%Maxelemlist ) THEN
           ALLOCATE(tmp(Mesh%maxelemlist+BUMPUP))
           ! if Mesh%Nelemlist .EQ. 1 then this is the first time through
           IF ( Mesh%Nelemlist .GT. 1 ) tmp(1:Mesh%nelemlist) = Mesh%elemlist(1:Mesh%Nelemlist)
           IF ( ASSOCIATED(Mesh%elemlist) ) THEN
               DEALLOCATE( Mesh%elemlist )
               NULLIFY( Mesh%elemlist )
           ENDIF
           Mesh%elemlist => tmp
           Mesh%Maxelemlist = Mesh%Maxelemlist + BUMPUP
         ENDIF
         Mesh%elemlist(Mesh%Nelemlist)%Element => Mesh%ElemTable(i)%Elements(j)
         Mesh%elemlist(Mesh%Nelemlist)%Element%Xelement = i
       ENDDO
     ENDDO

!#if 0
!     DO i = 1, Mesh%npoint
!       Mesh%Spatial = 0
!       CALL AddElement(Mesh, ELEMENT_POINT, i, Mesh%element_point )
!     ENDDO
!     DO i = 1, Mesh%nline2
!       Mesh%Spatial = 1
!       CALL AddElement(Mesh, ELEMENT_LINE2, i, Mesh%element_line2 )
!     ENDDO
!     DO i = 1, Mesh%nline3
!       Mesh%Spatial = 1
!       CALL AddElement(Mesh, ELEMENT_LINE3, i, Mesh%element_line3 )
!     ENDDO
!     DO i = 1, Mesh%ntri3
!       Mesh%Spatial = 2
!       CALL AddElement(Mesh, ELEMENT_TRI3, i, Mesh%element_tri3 )
!     ENDDO
!     DO i = 1, Mesh%ntri6
!       Mesh%Spatial = 2
!       CALL AddElement(Mesh, ELEMENT_TRI6, i, Mesh%element_tri6 )
!     ENDDO
!     DO i = 1, Mesh%nquad4
!       Mesh%Spatial = 2
!       CALL AddElement(Mesh, ELEMENT_QUAD4, i, Mesh%element_quad4 )
!     ENDDO
!     DO i = 1, Mesh%nquad8
!       Mesh%Spatial = 2
!       CALL AddElement(Mesh, ELEMENT_QUAD8, i, Mesh%element_quad8 )
!     ENDDO
!     DO i = 1, Mesh%ntet4
!       Mesh%Spatial = 3
!       CALL AddElement(Mesh, ELEMENT_TET4, i, Mesh%element_tet4 )
!     ENDDO
!     DO i = 1, Mesh%ntet10
!       Mesh%Spatial = 3
!       CALL AddElement(Mesh, ELEMENT_TET10, i, Mesh%element_tet10 )
!     ENDDO
!     DO i = 1, Mesh%nhex8
!       Mesh%Spatial = 3
!       CALL AddElement(Mesh, ELEMENT_HEX8, i, Mesh%element_hex8 )
!     ENDDO
!     DO i = 1, Mesh%nhex20
!       Mesh%Spatial = 3
!       CALL AddElement(Mesh, ELEMENT_HEX20, i, Mesh%element_hex20 )
!     ENDDO
!     DO i = 1, Mesh%nwedge6
!       Mesh%Spatial = 3
!       CALL AddElement(Mesh, ELEMENT_WEDGE6, i, Mesh%element_wedge6 )
!     ENDDO
!     DO i = 1, Mesh%nwedge15
!       Mesh%Spatial = 3
!       CALL AddElement(Mesh, ELEMENT_WEDGE15, i, Mesh%element_wedge15 )
!     ENDDO
!
!     CONTAINS
!       SUBROUTINE AddElement( Mesh, Xelement, i_in_kind, element_array )
!         IMPLICIT NONE
!         TYPE(MeshType), INTENT(INOUT) :: Mesh  ! Mesh being committed
!         INTEGER, INTENT(IN) :: Xelement        ! which kind of element
!         INTEGER, INTENT(IN) ::  i_in_kind      ! index of element in Xelement kind array
!         INTEGER, POINTER, INTENT(IN) :: element_array(:,:) ! element kind array
!       ! Local
!
!         TYPE(ElemRec), POINTER :: tmp(:)
!         Mesh%Nelemlist = Mesh%Nelemlist + 1
!         IF ( Mesh%Nelemlist .GE. Mesh%Maxelemlist ) THEN
!            ALLOCATE(tmp(Mesh%maxelemlist+BUMPUP))
!            ! if Mesh%Nelemlist .EQ. 1 then this is the first time through
!            IF ( Mesh%Nelemlist .GT. 1 ) tmp(1:Mesh%nelemlist) = Mesh%elemlist(1:Mesh%Nelemlist)
!            IF ( ASSOCIATED(Mesh%elemlist) ) THEN
!               DEALLOCATE( Mesh%elemlist )
!               NULLIFY( Mesh%elemlist )
!            ENDIF
!            Mesh%elemlist => tmp
!            Mesh%Maxelemlist = Mesh%Maxelemlist + BUMPUP
!         ENDIF
!         Mesh%elemlist(Mesh%Nelemlist)%ElemAry    => element_array
!         Mesh%elemlist(Mesh%Nelemlist)%I_in_kind = i_in_kind
!
!       END SUBROUTINE AddElement
!#endif
!
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
     ELSEIF ( P1 .LT. 1 .OR. P1 .GT. Mesh%Nnodes ) THEN !BJJ moved to ELSE because
       ErrStat = ErrID_Fatal
       ErrMess="MeshConstructElement_1PT: invalid P1 ("//TRIM(Num2LStr(P1))//") for mesh with "//TRIM(Num2LStr(Mesh%Nnodes))//" nodes."
     ENDIF
     IF ( ErrStat .NE. ErrID_None ) THEN
        CALL WrScr(TRIM(ErrMess))
        RETURN  !  early return on error
     ENDIF
    ! Business
     IF ( Xelement .EQ. ELEMENT_POINT ) THEN
       Mesh%ElemTable(ELEMENT_POINT)%nelem = Mesh%ElemTable(ELEMENT_POINT)%nelem + 1

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

  ! TODO : check to make sure mesh is committed first

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
     Mesh%nextelem = Mesh%nextelem + 1
     ErrStat = ErrID_None
     ErrMess = ""
     RETURN

   END SUBROUTINE MeshNextElement



END MODULE ModMesh


