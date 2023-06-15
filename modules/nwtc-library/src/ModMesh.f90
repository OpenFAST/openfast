!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
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
!> The modules ModMesh and ModMesh_Types provide data structures and subroutines for representing and manipulating meshes
!! and meshed data in the FAST modular framework. 
!!
!! A mesh is comprised of a set of "nodes" (simple points in space) together with information specifying how they are connected 
!! to form "elements"  representing spatial boundaries between components. ModMesh and ModMesh_Types define point, line, surface, 
!! and volume elements in a standard isoparametric mapping from finite element analysis. Currently only points and straight line 
!! (line2) elements are implemented.
!!   
!! Associated with a mesh are one or more "fields" that represent the values of variables or "degrees of freedom" at each node. 
!! A mesh always has a named "Position" that specifies the location in three-dimensional space as an Xi,Yi,Zi triplet of each node 
!! and a field named "RefOrientation" that specifies the orientation (as a direction cosine matrix) of the node. 
!! The ModMesh_Types module predefines a number of other fields of triples representing velocities, forces, and moments as well as
!! a field of nine values representing a direction cosine matrix. 
!!   
!! The operations on meshes defined in the ModMesh module are creation, spatio-location of nodes, construction, committing the 
!! mesh definition, initialization of fields, accessing field data, updating field data, copying, deallocating, and destroying meshes. 
!! See https://nwtc.nrel.gov/FAST-Developers and https://nwtc.nrel.gov/system/files/ProgrammingHandbook_Mod20130717.pdf
MODULE ModMesh
   use VTK, only: WrVTK_header, WrVTK_footer

   USE ModMesh_Types
   IMPLICIT NONE
!   INTEGER :: DEBUG_UNIT = 74

   INTEGER,     PARAMETER, PRIVATE :: BUMPUP = 64                !< size element list will be increased when adding an element that does not fit in the currently allocated space; do not set to less than 2
   CHARACTER(*),PARAMETER          :: VTK_AryFmt = '(3(F30.5))'  !< text format for triplets written to VTK text files


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

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes mesh information in binary form. If UnIn is < 0, it gets a new unit number and opens the file,
!! otherwise the file is appended. It is up to the caller of this routine to close the file when it's finished.
SUBROUTINE MeshWrBin ( UnIn, M, ErrStat, ErrMsg, FileName)

      
   INTEGER, INTENT(INOUT)                ::  UnIn     !< fortran output unit
   TYPE(MeshType),  INTENT(IN)           ::  M        !< mesh to be reported on

   INTEGER(IntKi),  INTENT(OUT)          :: ErrStat   !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: ErrMsg    !< Error message associated with the ErrStat
   CHARACTER(*),    INTENT(IN), OPTIONAL :: FileName  !< Name of the file to write the output in

   ! local variables
   INTEGER(IntKi)                        :: ErrStat2  ! Temporary storage for local errors
   INTEGER(IntKi)                        :: I         ! loop counter
   CHARACTER(*), PARAMETER               :: RoutineName = 'MeshWrBin'

   ErrStat = ErrID_None
   ErrMsg  = ""
   
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
   if (M%Fieldmask(MASKID_SCALAR))  WRITE (UnIn, IOSTAT=ErrStat2)   INT(M%nScalars,B4Ki)


   !...........
   ! Write nodal information:
   !...........
      
   IF (.NOT. M%Initialized) RETURN
      
   WRITE (UnIn, IOSTAT=ErrStat2)   M%Position
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error writing Position to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF

   WRITE (UnIn, IOSTAT=ErrStat2)   M%RefOrientation
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, 'Error writing RefOrientation to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
         RETURN
      END IF


   ! Write fields:

   IF ( M%fieldmask(MASKID_FORCE) .AND. ALLOCATED(M%Force) ) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   M%Force
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error writing Force to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
   END IF

   IF ( M%fieldmask(MASKID_Moment) .AND. ALLOCATED(M%Moment)) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   M%Moment
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error writing Force to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
   END IF

   IF ( M%fieldmask(MASKID_ORIENTATION) .AND. ALLOCATED(M%Orientation)) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   M%Orientation
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error writing Orientation to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
   END IF

   IF ( M%fieldmask(MASKID_TRANSLATIONDISP) .AND. ALLOCATED(M%TranslationDisp)) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   M%TranslationDisp
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error writing TranslationDisp to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
   END IF

   IF ( M%fieldmask(MASKID_TRANSLATIONVEL) .AND. ALLOCATED(M%TranslationVel)) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   M%TranslationVel
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error writing TranslationVel to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
   END IF

   IF ( M%fieldmask(MASKID_ROTATIONVEL) .AND. ALLOCATED(M%RotationVel)) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   M%RotationVel
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error writing RotationVel to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
   END IF

   IF ( M%fieldmask(MASKID_TRANSLATIONACC) .AND. ALLOCATED(M%TranslationAcc)) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   M%TranslationAcc
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error writing TranslationAcc to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
   END IF

   IF ( M%fieldmask(MASKID_ROTATIONACC) .AND. ALLOCATED(M%RotationAcc)) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   M%RotationAcc
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error writing RotationAcc to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
            RETURN
         END IF
   END IF

   IF ( M%fieldmask(MASKID_SCALAR) .AND. ALLOCATED(M%Scalars)) THEN
      WRITE (UnIn, IOSTAT=ErrStat2)   M%Scalars
         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error writing Scalars to the mesh binary file.', ErrStat, ErrMsg, RoutineName )
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

END SUBROUTINE MeshWrBin
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes the reference position and orientations of a mesh in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
SUBROUTINE MeshWrVTKreference (RefPoint, M, FileRootName, ErrStat, ErrMsg )
   
   REAL(SiKi),      INTENT(IN)           :: RefPoint(3)   !< reference location, normally (0,0,0)
   TYPE(MeshType),  INTENT(IN)           :: M             !< mesh to be written
   CHARACTER(*),    INTENT(IN)           :: FileRootName  !< Name of the file to write the output in (excluding extension)

   INTEGER(IntKi),  INTENT(OUT)          :: ErrStat       !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: ErrMsg        !< Error message associated with the ErrStat
   
   ! local variables
   INTEGER(IntKi)                        :: Un            ! fortran unit number
   INTEGER(IntKi)                        :: I, J          ! loop counters
      
   INTEGER(IntKi)                        :: ErrStat2 
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*),PARAMETER                :: RoutineName = 'MeshWrVTKreference'
      
   CHARACTER(*),PARAMETER                :: RefOrientation(3) = (/ 'RefOrientationX','RefOrientationY','RefOrientationZ' /)
      
      
   ErrStat = ErrID_None
   ErrMsg  = ""
      
      
   ! PolyData (.vtp) - Serial vtkPolyData (unstructured)
   call WrVTK_header( TRIM(FileRootName)//'_Reference.vtp', M%Nnodes, M%ElemTable(ELEMENT_LINE2)%nelem, 0, Un, ErrStat2, ErrMsg2 )    
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
         
! points (i.e., nodes):      
      WRITE(Un,'(A)')         '      <Points>'         
      WRITE(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i)
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Points>'
   
! point data (orientation vectors):
      WRITE(Un,'(A)')         '      <PointData>'        
   DO j=1,3 
      WRITE(Un,'(A,A,A)')   '        <DataArray type="Float32" Name="', RefOrientation(j), '" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         WRITE(Un,VTK_AryFmt) RefPoint + M%RefOrientation(j,:,i)
      END DO
      WRITE(Un,'(A)')      '        </DataArray>'
   END DO
      WRITE(Un,'(A)')         '      </PointData>'
   
! lines (i.e., elements; for line2 meshes only):
   if ( M%ElemTable(ELEMENT_LINE2)%nelem > 0) then
      WRITE(Un,'(A)')         '      <Lines>'
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'
      DO i=1,M%ElemTable(ELEMENT_LINE2)%nelem
         WRITE(Un,'(2(i7))') M%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1) - 1, M%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2) - 1
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'
      DO i=1,M%ElemTable(ELEMENT_LINE2)%nelem
         WRITE(Un,'(i7)') 2*i
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Lines>'
   end if

   call WrVTK_footer( Un )
      
END SUBROUTINE MeshWrVTKreference   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes mesh information in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
SUBROUTINE MeshWrVTK ( RefPoint, M, FileRootName, VTKcount, OutputFieldData, ErrStat, ErrMsg, Twidth, Sib )
      
   REAL(SiKi),      INTENT(IN)           :: RefPoint(3)     !< reference location, normally (0,0,0)
   TYPE(MeshType),  INTENT(IN)           :: M               !< mesh to be written
   CHARACTER(*),    INTENT(IN)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   INTEGER(IntKi),  INTENT(IN)           :: VTKcount        !< Indicates number for VTK output file (when 0, the routine will also write reference information)
   LOGICAL,         INTENT(IN)           :: OutputFieldData !< flag to determine if we want to output field data or just the absolute position of this mesh
   INTEGER(IntKi),  INTENT(OUT)          :: ErrStat         !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: ErrMsg          !< Error message associated with the ErrStat
   INTEGER(IntKi),  INTENT(IN)           :: Twidth          !< Number of digits in the maximum write-out step (used to pad the VTK write-out in the filename with zeros)

   TYPE(MeshType),  INTENT(IN), OPTIONAL :: Sib             !< "functional" Sibling of M that contains translational displacement information (used to place forces at displaced positions) 

   ! local variables
   INTEGER(IntKi)                        :: Un            ! fortran unit number
   INTEGER(IntKi)                        :: i,j           ! loop counters
   CHARACTER(1024)                       :: FileName
   CHARACTER(Twidth)                     :: Tstr          ! string for current VTK write-out step (padded with zeros)
      
   INTEGER(IntKi)                        :: ErrStat2 
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*),PARAMETER                :: RoutineName = 'MeshWrVTK'
   CHARACTER(*),PARAMETER                :: Orientation(3) = (/ 'OrientationX','OrientationY','OrientationZ' /)

   
   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (.NOT. M%Initialized) RETURN
         
   !.................................................................
   !> We'll write the mesh reference fields on the first timestep only:
   !.................................................................
   if (VTKcount == 0) then
      call MeshWrVTKreference(RefPoint, M, FileRootName, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= AbortErrLev) return
   end if

   !.................................................................
   ! write the data that potentially changes each time step:
   !.................................................................

   ! construct the string for the zero-padded VTK write-out step
   write(Tstr, '(i' // trim(Num2LStr(Twidth)) //'.'// trim(Num2LStr(Twidth)) // ')') VTKcount
      
   ! PolyData (.vtp) - Serial vtkPolyData (unstructured) file
   FileName = TRIM(FileRootName)//'.'//Tstr//'.vtp'
      
   call WrVTK_header( trim(FileName), M%Nnodes, M%ElemTable(ELEMENT_LINE2)%nelem, 0, Un, ErrStat2, ErrMsg2 )    
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
            
      ! Write a VTP mesh file (Polygonal VTK file) with positions, lines, and field information
      ! (note alignment of WRITE statements to make sure spaces are lined up in XML file)
   
! points (nodes):   
      WRITE(Un,'(A)')         '      <Points>'
      WRITE(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      IF (ALLOCATED(M%TranslationDisp)) THEN
         DO i=1,M%Nnodes
            WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i) + M%TranslationDisp(:,i)
         END DO
      ELSEIF ( PRESENT(Sib) ) THEN
         if (allocated(Sib%TranslationDisp)) then
            DO i=1,M%Nnodes
               WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i) + Sib%TranslationDisp(:,i) ! @note: M%Position and Sib%Position should be the same!
            END DO
         else
            DO i=1,M%Nnodes
               WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i)
            END DO
         end if         
      ELSE         
         DO i=1,M%Nnodes
            WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i)
         END DO
      END IF
   
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Points>'
            
   if (OutputFieldData) then            ! point data for any existing mesh fields:   
      WRITE(Un,'(A)')         '      <PointData>'
      call MeshWrVTKfields ( Un, M, 1)
      
      if ( PRESENT(Sib) ) then ! write the sibling fields, too, so we don't have so many output files
         if (Sib%Nnodes == M%Nnodes .and. Sib%nelemlist == M%nelemlist ) then
            call MeshWrVTKfields ( Un, Sib, 1)
         end if         
      end if
      WRITE(Un,'(A)')         '      </PointData>'
   end if !(OutputFieldData)      

      
! lines (i.e., elements; for line2 meshes only):
   if ( M%ElemTable(ELEMENT_LINE2)%nelem > 0) then    
      WRITE(Un,'(A)')         '      <Lines>'
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'
      DO i=1,M%ElemTable(ELEMENT_LINE2)%nelem
         WRITE(Un,'(2(i7))') M%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1) - 1, M%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2) - 1
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'
      DO i=1,M%ElemTable(ELEMENT_LINE2)%nelem
         WRITE(Un,'(i7)') 2*i
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Lines>'
   end if      

      call WrVTK_footer( Un )               
      
END SUBROUTINE MeshWrVTK
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes mesh field information in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
SUBROUTINE MeshWrVTKfields ( Un, M, n )
      
   INTEGER(IntKi),  INTENT(IN)           :: Un            !< unit number of already-open vtk file in which to write the field information
   TYPE(MeshType),  INTENT(IN)           :: M             !< mesh to be written
   INTEGER(IntKi),  INTENT(IN)           :: n             !< number of times to write field value for each mesh node (> 1 when added to surface)


   ! local variables
   INTEGER(IntKi)                        :: i,j,k         ! loop counters
   CHARACTER(1024)                       :: FileName
      
   !INTEGER(IntKi)                        :: ErrStat2 
   !CHARACTER(ErrMsgLen)                  :: ErrMsg2
   !CHARACTER(*),PARAMETER                :: RoutineName = 'MeshWrVTKfields'
   CHARACTER(*),PARAMETER                :: Orientation(3) = (/ 'OrientationX','OrientationY','OrientationZ' /)

   
   
! point data for any existing mesh fields:   
      !WRITE(Un,'(A)')         '      <PointData>'
      
   IF ( M%fieldmask(MASKID_FORCE) .AND. ALLOCATED(M%Force) ) THEN
      WRITE(Un,'(A)')         '        <DataArray type="Float32" Name="Force" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         do k=1,n
            WRITE(Un,VTK_AryFmt) M%Force(:,i)
         end do !k         
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
   END IF

   IF ( M%fieldmask(MASKID_MOMENT) .AND. ALLOCATED(M%Moment) ) THEN
      WRITE(Un,'(A)')         '        <DataArray type="Float32" Name="Moment" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         do k=1,n
            WRITE(Un,VTK_AryFmt) M%Moment(:,i)
         end do !k
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
   END IF

#ifdef VTK_OUTPUT_TRANSLATIONDISP
   IF ( M%fieldmask(MASKID_TRANSLATIONDISP) .AND. ALLOCATED(M%TranslationDisp)) THEN
      WRITE(Un,'(A)')         '        <DataArray type="Float32" Name="TranslationalDisplacement" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         do k=1,n
            WRITE(Un,VTK_AryFmt) M%TranslationDisp(:,i)
         end do         
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
   END IF
#endif

   IF ( M%fieldmask(MASKID_TRANSLATIONVEL) .AND. ALLOCATED(M%TranslationVel)) THEN
      WRITE(Un,'(A)')         '        <DataArray type="Float32" Name="TranslationalVelocity" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         do k=1,n
            WRITE(Un,VTK_AryFmt) M%TranslationVel(:,i)
         end do         
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
   END IF

   IF ( M%fieldmask(MASKID_ROTATIONVEL) .AND. ALLOCATED(M%RotationVel)) THEN
      WRITE(Un,'(A)')         '        <DataArray type="Float32" Name="RotationalVelocity" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         do k=1,n
            WRITE(Un,VTK_AryFmt) M%RotationVel(:,i)
         end do !k         
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
   END IF

   IF ( M%fieldmask(MASKID_TRANSLATIONACC) .AND. ALLOCATED(M%TranslationAcc)) THEN
      WRITE(Un,'(A)')         '        <DataArray type="Float32" Name="TranslationalAcceleration" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         do k=1,n
            WRITE(Un,VTK_AryFmt) M%TranslationAcc(:,i)
         end do !k
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
   END IF

   IF ( M%fieldmask(MASKID_ROTATIONACC) .AND. ALLOCATED(M%RotationAcc)) THEN
      WRITE(Un,'(A)')         '        <DataArray type="Float32" Name="RotationalAcceleration" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         do k=1,n
            WRITE(Un,VTK_AryFmt) M%RotationAcc(:,i)
         end do         
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
   END IF

IF (M%fieldmask(MASKID_ORIENTATION) .AND. ALLOCATED(M%Orientation)) THEN
   DO j=1,3 
      WRITE(Un,'(A,A,A)')   '        <DataArray type="Float32" Name="', Orientation(j), '" NumberOfComponents="3" format="ascii">'
      DO i=1,M%Nnodes
         do k=1,n
            WRITE(Un,VTK_AryFmt) M%Orientation(j,:,i)
         end do         
      END DO
      WRITE(Un,'(A)')      '        </DataArray>'
   END DO      
END IF

   IF ( M%fieldmask(MASKID_SCALAR) .AND. ALLOCATED(M%Scalars) .AND. M%nScalars > 0) THEN
      WRITE(Un,'(A,I7,A)')         '        <DataArray type="Float32" Name="Scalars" NumberOfComponents="', M%nScalars, '" format="ascii">'
      DO i=1,M%Nnodes
         do k=1,n
            WRITE(Un,'('//trim(num2lstr(M%nScalars))//'(F20.6))') M%Scalars(:,i) ! not very efficient, but it's easy and I'm not sure anyone uses this field
         end do         
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'
   END IF

      !WRITE(Un,'(A)')         '      </PointData>'
      
               
      
END SUBROUTINE MeshWrVTKfields
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes line2 mesh surface information in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
SUBROUTINE MeshWrVTK_Ln2Surface ( RefPoint, M, FileRootName, VTKcount, OutputFieldData, ErrStat, ErrMsg, Twidth, NumSegments, Radius, verts, Sib )
      
   REAL(SiKi),      INTENT(IN)           :: RefPoint(3)     !< reference location, normally (0,0,0)
   TYPE(MeshType),  INTENT(IN)           :: M               !< mesh to be written
   CHARACTER(*),    INTENT(IN)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   INTEGER(IntKi),  INTENT(IN)           :: VTKcount        !< Indicates number for VTK output file (when 0, the routine will also write reference information)
   LOGICAL,         INTENT(IN)           :: OutputFieldData !< flag to determine if we want to output field data or just the absolute position of this mesh
   INTEGER(IntKi),  INTENT(IN)           :: Twidth          !< Number of digits in the maximum write-out step (used to pad the VTK write-out in the filename with zeros)
   INTEGER(IntKi),  INTENT(IN), OPTIONAL :: NumSegments     !< Number of segments to split the circle into
   REAL(SiKi),      INTENT(IN), OPTIONAL :: Radius(:)       !< Radius of each node
   REAL(SiKi),      INTENT(IN), OPTIONAL :: verts(:,:,:)    !< X-Y verticies (2x{NumSegs}xNNodes) of points that define a shape around each node
   TYPE(MeshType),  INTENT(IN), OPTIONAL :: Sib             !< Sibling of M that contains more field information (used only if OutputFieldData is true, to minimize number of files being written)
   
   INTEGER(IntKi),  INTENT(OUT)          :: ErrStat         !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: ErrMsg          !< Error message associated with the ErrStat


   ! local variables
   INTEGER(IntKi)                        :: Un            ! fortran unit number
   INTEGER(IntKi)                        :: i,j           ! loop counters
   INTEGER(IntKi)                        :: offset_cnt    ! counter for offsets (corresponding to number of nodes in polygons written)
   CHARACTER(1024)                       :: FileName
   REAL(SiKi)                            :: angle
   REAL(SiKi)                            :: xyz(3)
   CHARACTER(Twidth)                     :: Tstr          ! string for current write-out step (padded with zeros)

   INTEGER(IntKi)                        :: firstPntEnd, firstPntStart, secondPntStart, secondPntEnd  ! node indices for forming rectangle 
   INTEGER(IntKi)                        :: NumSegments1
   
   
   
   INTEGER(IntKi)                        :: ErrStat2 
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*),PARAMETER                :: RoutineName = 'MeshWrVTK_Ln2Surface'

   
   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (.NOT. M%Initialized) RETURN
   IF (.NOT. ALLOCATED(M%TranslationDisp) ) RETURN
   IF (.NOT. ALLOCATED(M%Orientation) ) RETURN
      
   if (present(verts)) then
      NumSegments1   = size(verts,2)            
   elseif (present(Radius) .and. present(NumSegments)) then
      NumSegments1 = NumSegments
   else
      call SetErrStat(ErrID_Fatal,'Incorrect number of arguments.',ErrStat,ErrMsg,RoutineName)
      RETURN
   end if   
   
   !.................................................................
   ! write the data that potentially changes each time step:
   !.................................................................

   ! construct the string for the zero-padded VTK write-out step
   write(Tstr, '(i' // trim(Num2LStr(Twidth)) //'.'// trim(Num2LStr(Twidth)) // ')') VTKcount
      
   ! PolyData (.vtp) - Serial vtkPolyData (unstructured) file
   FileName = TRIM(FileRootName)//'.'//Tstr//'.vtp'
       
      ! Write a VTP mesh file (Polygonal VTK file) with positions and polygons (surfaces)
      ! (note alignment of WRITE statements to make sure spaces are lined up in XML file)
   call WrVTK_header(   FileName=trim(FileName)                                        &
                      , NumberOfPoints=M%Nnodes*NumSegments1                           &
                      , NumberOfLines=M%ElemTable(ELEMENT_LINE2)%nelem                 &
                      , NumberOfPolys=M%ElemTable(ELEMENT_LINE2)%nelem*(NumSegments1+2)& 
                      , Un=Un                                                          &
                      , ErrStat=ErrStat2                                               &
                      , ErrMsg=ErrMsg2                                                 )  
   
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
                     
! points (nodes, augmented with NumSegments):   
      WRITE(Un,'(A)')         '      <Points>'
      WRITE(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      
      xyz(3) = 0.0_SiKi
      if (present(verts)) then
         
         DO i=1,M%Nnodes
            DO j=1,NumSegments1
               xyz(1:2) = verts(1:2,j,i)
               WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i) + M%TranslationDisp(:,i) + matmul(xyz,M%Orientation(:,:,i))
            END DO
         END DO         
      else               
         DO i=1,M%Nnodes
            DO j=1,NumSegments1
               angle = TwoPi*(j-1.0_ReKi)/NumSegments1
               xyz(1) = radius(i)*COS(angle)
               xyz(2) = radius(i)*SIN(angle)
               WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i) + M%TranslationDisp(:,i) + matmul(xyz,M%Orientation(:,:,i))
            END DO
         END DO
      end if
      
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Points>'
  
      
   if (OutputFieldData) then            ! point data for any existing mesh fields:   
      WRITE(Un,'(A)')         '      <PointData>'
      call MeshWrVTKfields ( Un, M, NumSegments1)
      
      if ( PRESENT(Sib) ) then ! write the sibling fields, too, so we don't have so many output files
         if (Sib%Nnodes == M%Nnodes .and. Sib%nelemlist == M%nelemlist ) then
            call MeshWrVTKfields ( Un, Sib, NumSegments1)
         end if         
      end if
      WRITE(Un,'(A)')         '      </PointData>'
   end if !(OutputFieldData)      

   if ( M%ElemTable(ELEMENT_LINE2)%nelem > 0) then   
         ! Using rectangle to render surfaces (for line2 meshes only):
      WRITE(Un,'(A)')         '      <Polys>'
      
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'      
      DO i=1,M%ElemTable(ELEMENT_LINE2)%nelem
         firstPntStart  = (M%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(1)-1)*NumSegments1
         firstPntEnd    = firstPntStart + NumSegments1 - 1
         secondPntStart = (M%ElemTable(ELEMENT_LINE2)%Elements(i)%ElemNodes(2)-1)*NumSegments1
         secondPntEnd   = secondPntStart + NumSegments1 - 1
         DO j=1,NumSegments1-1
            WRITE(Un,'(4(i7))') firstPntStart + (j-1), firstPntStart + j, &
                                secondPntStart + j,    secondPntStart + (j-1)
         END DO
         WRITE(Un,'(4(i7))')  firstPntEnd, firstPntStart, secondPntStart, secondPntEnd
         
         ! make top and bottom of this element, making sure surface normals point outward
         WRITE(Un,'('//trim(num2lstr(NumSegments1))//'(i7))') (j, j=firstPntEnd,firstPntStart,-1)               
         WRITE(Un,'('//trim(num2lstr(NumSegments1))//'(i7))') (j, j=secondPntStart,secondPntEnd)              
                  
      END DO      
      WRITE(Un,'(A)')         '        </DataArray>'

      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'
      offset_cnt = 0
      DO i=1,M%ElemTable(ELEMENT_LINE2)%nelem
         DO j=1,NumSegments1
            offset_cnt = offset_cnt + 4 ! number of nodes in polygon
            WRITE(Un,'(i7)') offset_cnt
         END DO
         DO j=1,2 ! top and bottom
            offset_cnt = offset_cnt + NumSegments1 ! number of nodes in this polygon
            WRITE(Un,'(i7)') offset_cnt
         END DO
      END DO
      WRITE(Un,'(A)')         '        </DataArray>'

      WRITE(Un,'(A)')         '      </Polys>'      
            
   end if ! do this only for line2 elements    

      call WrVTK_footer( Un )         
                     
   END SUBROUTINE MeshWrVTK_Ln2Surface
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine writes point mesh surfaces information in VTK format.
!! see VTK file information format for XML, here: http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
SUBROUTINE MeshWrVTK_PointSurface ( RefPoint, M, FileRootName, VTKcount, OutputFieldData, ErrStat, ErrMsg, Twidth, NumSegments, Radius, verts, Sib )
      
   REAL(SiKi),      INTENT(IN)           :: RefPoint(3)     !< reference location, normally (0,0,0)
   TYPE(MeshType),  INTENT(IN)           :: M               !< mesh to be written
   CHARACTER(*),    INTENT(IN)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   INTEGER(IntKi),  INTENT(IN)           :: VTKcount        !< Indicates number for VTK output file (when 0, the routine will also write reference information)
   LOGICAL,         INTENT(IN)           :: OutputFieldData !< flag to determine if we want to output field data or just the absolute position of this mesh
   INTEGER(IntKi),  INTENT(IN)           :: Twidth          !< Number of digits in the maximum write-out timestep (used to pad the VTK write-out in the filename with zeros)
   INTEGER(IntKi),  INTENT(IN), OPTIONAL :: NumSegments     !< Number of segments to split the circle into
   REAL(SiKi),      INTENT(IN), OPTIONAL :: Radius          !< Radius of each node
   REAL(SiKi),      INTENT(IN), OPTIONAL :: verts(:,:)      !< X-Y-Z verticies (3xn) of points that define a volume around each node
   !bjj: we don't need this to be limited to 8, I guess...
   TYPE(MeshType),  INTENT(IN), OPTIONAL :: Sib             !< Sibling of M that contains more field information (used only if OutputFieldData is true, to minimize number of files being written)
   
   INTEGER(IntKi),  INTENT(OUT)          :: ErrStat         !< Indicates whether an error occurred (see NWTC_Library)
   CHARACTER(*),    INTENT(OUT)          :: ErrMsg          !< Error message associated with the ErrStat


   ! local variables
   INTEGER(IntKi)                        :: Un            ! fortran unit number
   INTEGER(IntKi)                        :: i,j,k         ! loop counters
   INTEGER(IntKi)                        :: offset_cnt    ! counter for offsets (corresponding to number of nodes in polygons written)
   INTEGER(IntKi)                        :: NumberOfPoints, NumberOfPointsPerNode
   INTEGER(IntKi)                        :: NumberOfPolys
   INTEGER(IntKi)                        :: NumSegments1   
   INTEGER(IntKi)                        :: NumSegments2   
   CHARACTER(1024)                       :: FileName
   REAL(SiKi)                            :: angle, r, ratio
   REAL(SiKi)                            :: xyz(3)
   CHARACTER(Twidth)                     :: Tstr          ! string for current VTK write-out step (padded with zeros)

   INTEGER(IntKi)                        :: firstPntEnd, firstPntStart, secondPntStart, secondPntEnd  ! node indices for forming rectangle 
   
   
   
   INTEGER(IntKi)                        :: ErrStat2 
   CHARACTER(ErrMsgLen)                  :: ErrMsg2
   CHARACTER(*),PARAMETER                :: RoutineName = 'MeshWrVTK_PointSurface'

   
   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (.NOT. M%Initialized) RETURN
   IF (.NOT. ALLOCATED(M%TranslationDisp) ) RETURN
   IF (.NOT. ALLOCATED(M%Orientation) ) RETURN
      
   if (present(verts)) then
      NumberOfPointsPerNode = size(verts,2)
      if (size(verts,2)==8) then
         NumberOfPolys = 6 ! per node
      elseif (size(verts,2)==4) then
         NumberOfPolys = 1 ! per node
      else
         ! it would be nice if we could add this sometime, but ...
         call SetErrStat(ErrID_Fatal,'When verticies are specified, there must be exactly 4 or 8.',ErrStat,ErrMsg,RoutineName)
         RETURN
      end if
      
      NumberOfPolys  = M%Nnodes * NumberOfPolys
      
   elseif (present(Radius) ) then
      if (present(NumSegments)) then ! a volume
         NumSegments1 = max(1,abs(NumSegments))
         NumSegments2 = max(4,NumSegments1)         
         
         NumberOfPolys  = M%Nnodes * max(1,(NumSegments2-1))* NumSegments1
         NumberOfPointsPerNode = NumSegments2 * NumSegments1
         
      else                           ! a plane
         NumSegments1 = 20
         NumSegments2 = 1
         
         NumberOfPolys  = M%Nnodes 
         NumberOfPointsPerNode = NumSegments1
         
      end if
      
      
   else
      call SetErrStat(ErrID_Fatal,'Incorrect number of arguments.',ErrStat,ErrMsg,RoutineName)
      RETURN
   end if
   
   NumberOfPoints = M%Nnodes*NumberOfPointsPerNode
   !.................................................................
   ! write the data that potentially changes each time step:
   !.................................................................

   ! construct the string for the zero-padded VTK write-out step
   write(Tstr, '(i' // trim(Num2LStr(Twidth)) //'.'// trim(Num2LStr(Twidth)) // ')') VTKcount
      
   ! PolyData (.vtp) - Serial vtkPolyData (unstructured) file
   FileName = TRIM(FileRootName)//'.'//Tstr//'.vtp'
      
      ! Write a VTP mesh file (Polygonal VTK file) with positions and polygons (surfaces)
      ! (note alignment of WRITE statements to make sure spaces are lined up in XML file)
   call WrVTK_header( trim(FileName), M%Nnodes*NumberOfPoints, 0, NumberOfPolys, Un, ErrStat2, ErrMsg2 )   
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
                        
! points (nodes, augmented with NumSegments):   
      WRITE(Un,'(A)')         '      <Points>'
      WRITE(Un,'(A)')         '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
      
      if ( present(verts) ) then
         do i=1,M%Nnodes
            do j=1,NumberOfPointsPerNode
               WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i) + M%TranslationDisp(:,i) + MATMUL(verts(:,j),M%Orientation(:,:,i))
            end do    
         end do
         
      else                           
         DO i=1,M%Nnodes
            DO j=1,NumSegments2
               if (NumSegments2>1) then
                  ratio  = 2.0_SiKi*REAL(j-1,SiKi)/REAL(NumSegments2-1,SiKi) - 1.0_SiKi  ! where we are in [-1, 1]
               else
                  ratio = 0.0_SiKi
                  !WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i) + M%TranslationDisp(:,i)  ! write center of node
               end if               
               
               xyz(3) = radius*ratio
               
               ! now calculate the radius of the x-y circle we're going to create at this z:              
               !r = acos( radius / xyz(3) ) or r = sqrt( radius**2 - xyz(3)**2 ) = radius*sqrt(1 - ratio**2)
               r = radius*sqrt(abs(1.0_SiKi - ratio**2))  ! note the abs in case ratio**2 gets slightly larger than 1 
               
               DO k=1,NumSegments1
                  angle = TwoPi*(k-1.0_ReKi)/NumSegments1
                  xyz(1) = r*COS(angle)
                  xyz(2) = r*SIN(angle)
                  WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i) + M%TranslationDisp(:,i) + MATMUL(xyz,M%Orientation(:,:,i))
                  !WRITE(Un,VTK_AryFmt) RefPoint + M%Position(:,i) + M%TranslationDisp(:,i) + xyz
               END DO            
            END DO
         END DO         
      end if
      
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Points>'
  
      
   if (OutputFieldData) then            ! point data for any existing mesh fields:   
      WRITE(Un,'(A)')         '      <PointData>'
      call MeshWrVTKfields ( Un, M, NumberOfPointsPerNode)      
      
      if ( PRESENT(Sib) ) then ! write the sibling fields, too, so we don't have so many output files
         if (Sib%Nnodes == M%Nnodes .and. Sib%nelemlist == M%nelemlist ) then
            call MeshWrVTKfields ( Un, Sib, NumberOfPointsPerNode)         
         end if         
      end if
      WRITE(Un,'(A)')         '      </PointData>'
   end if !(OutputFieldData)            
                  
      WRITE(Un,'(A)')         '      <Polys>'
      
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="connectivity" format="ascii">'   
      
      offset_cnt = 4
      if ( present(verts) ) then
         if (size(verts,2)==8) then
            
               ! Write points for the 6 corners of the box
            
            do i=1,M%Nnodes
               firstPntStart  = (i-1)*size(verts,2)

               WRITE(Un,'(4(i7))') firstPntStart  ,firstPntStart+1,firstPntStart+2,firstPntStart+3    ! bottom
               WRITE(Un,'(4(i7))') firstPntStart+4,firstPntStart+5,firstPntStart+6,firstPntStart+7    ! top
               
               WRITE(Un,'(4(i7))') firstPntStart+7,firstPntStart+6,firstPntStart+1,firstPntStart      ! sides
               WRITE(Un,'(4(i7))') firstPntStart+6,firstPntStart+5,firstPntStart+2,firstPntStart+1
               WRITE(Un,'(4(i7))') firstPntStart+5,firstPntStart+4,firstPntStart+3,firstPntStart+2
               WRITE(Un,'(4(i7))') firstPntStart+4,firstPntStart+7,firstPntStart  ,firstPntStart+3                  
            end do
            
         elseif (size(verts,2)==4) then
            
               ! Write points for the 4 corners of the polygon
            
            do i=1,M%Nnodes
               firstPntStart  = (i-1)*size(verts,2)
               
               WRITE(Un,'(4(i7))') firstPntStart  ,firstPntStart+1,firstPntStart+2,firstPntStart+3    ! bottom
            end do
            
         end if
      else
            
         if (NumSegments2==1) then
            offset_cnt = NumSegments1 ! number of nodes in this polygon
            
            DO i=1,M%Nnodes
               WRITE(Un,'('//trim(num2lstr(NumSegments1))//'(i7))') (k, k=0,NumSegments1-1)               
            END DO             
         else
            
            DO i=1,M%Nnodes
               DO j=1,NumSegments2-1
                  firstPntStart  = (j-1)*NumSegments1 + (i-1)*NumSegments1*NumSegments2
                  firstPntEnd    = firstPntStart + NumSegments1 - 1
                  secondPntStart = j*NumSegments1 + (i-1)*NumSegments1*NumSegments2
                  secondPntEnd   = secondPntStart + NumSegments1 - 1
                  DO k=1,NumSegments1-1
                     WRITE(Un,'(4(i7))') firstPntStart + (k-1), firstPntStart + k, &
                                         secondPntStart + k,    secondPntStart + (k-1)
                  END DO
                  WRITE(Un,'(4(i7))')  firstPntEnd, firstPntStart, secondPntStart, secondPntEnd
               END DO         
            END DO      
            
         end if         
      end if
      
      WRITE(Un,'(A)')         '        </DataArray>'
      
      WRITE(Un,'(A)')         '        <DataArray type="Int32" Name="offsets" format="ascii">'      
      
      do i=1,NumberOfPolys
         WRITE(Un,'(i7)') offset_cnt*i
      end do      
      WRITE(Un,'(A)')         '        </DataArray>'
      WRITE(Un,'(A)')         '      </Polys>'      
            

      call WrVTK_footer( Un )         
                     
   END SUBROUTINE MeshWrVTK_PointSurface
   
!-------------------------------------------------------------------------------------------------------------------------------
!> This routine writes mesh information in text form. It is used for debugging.
   SUBROUTINE MeshPrintInfo ( U, M, N, MeshName)
         
     INTEGER, INTENT(IN   )                ::      U  !< fortran output unit
     TYPE(MeshType),INTENT(IN   )          ::      M  !< mesh to be reported on
     INTEGER, OPTIONAL,INTENT(IN   )       ::      N  !< Number to print, default is all nodes
     character(*), optional, intent(in   ) :: MeshName !< name of the mesh
    ! Local
     INTEGER isz,i,j,nn,Ielement,Xelement

     nn = M%Nnodes !5
     IF (PRESENT(N)) nn = min(nn,N)

     if (present(MeshName)) then
        write(U,*)'-----------  MeshPrintInfo: '//trim(MeshName)//'  -------------'
     else
        write(U,*)'-----------  MeshPrintInfo:  -------------'
     endif

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


     IF (ASSOCIATED(M%ElemTable)) THEN
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
      END IF

! Here are some built in derived data types that can represent values at the nodes
! the last dimension of each of these has range 1:nnodes for the mesh being represented
! and they are indexed by the element arrays above
! only some of these would be allocated, depending on what's being represented
! on the mesh.
! Whether or not these are allocated is indicted in the fieldmask, which can
! be interrogated by a routine using an instance of the type. If you add a field
! here, be sure to change the table of parameters used to size and index fieldmask above.
     IF(ASSOCIATED(M%Position))THEN
       isz=size(M%Position,2)
       write(U,*)'Position: ',isz,' node(s)'
       DO i=1,min(nn,isz)
         write(U,*)' ',i,M%Position(:,i)
       ENDDO
     ENDIF
     IF(ASSOCIATED(M%RefOrientation))THEN
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
     
     DO Ielement=1,M%nelemlist
        Xelement = M%ElemList(Ielement)%Element%Xelement
        
        WRITE(U,'("  Ielement: ",I10,1x,A," det_jac: ",ES15.7," Nodes:",'//&
                Num2LStr( size(M%ElemList(Ielement)%Element%ElemNodes) )//'(1x,I10))') &
                    Ielement,&
                    ElemNames(Xelement), &
                    M%ElemList(Ielement)%Element%det_jac,  &
                    M%ElemList(Ielement)%Element%ElemNodes                      
     END DO 
          
     write(U,*)'---------  End of Element List  ----------'

   END SUBROUTINE MeshPrintInfo

!----------------------------------------------------------------------------------------------------------------------------------
   ! operations to create a mesh

!> Takes a blank, uninitialized instance of Type(MeshType) and defines the number of nodes in the mesh. Optional 
!! arguments indicate the fields that will be allocated and associated with the nodes of the mesh. The fields that may 
!! be associated with the mesh nodes are Force, Moment, Orientation, Rotation, TranslationDisp, RotationVel, TranslationVel, 
!! RotationAcc, TranslationAcc, and an arbitrary number of Scalars. See the definition of ModMeshType for descriptions of these fields.  
! After the first 5 arguments, the others are optional that say whether to allocate fields in the mesh.
! These are always dimensioned npoints 
   SUBROUTINE MeshCreate ( BlankMesh                                                       &
                          ,IOS                                                             &
                          ,Nnodes                                                          &
                          ,ErrStat                                                         &
                          ,ErrMess                                                         &
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
      
      TYPE(MeshType), INTENT(INOUT)   :: BlankMesh !< Mesh to be created
      INTEGER,INTENT(IN)         :: IOS                  !< input (COMPONENT_INPUT), output(COMPONENT_OUTPUT), or state(COMPONENT_STATE)
      INTEGER,INTENT(IN)         :: Nnodes               !< Number of nodes in mesh
      INTEGER(IntKi),INTENT(OUT) :: ErrStat              !< error status/level
      CHARACTER(*),INTENT(OUT)   :: ErrMess              !< error message
                                   ! optional arguments from here down
                                   ! optional arguments that say whether to allocate fields
                                   ! in the mesh. These are always dimensioned npoints
      LOGICAL,OPTIONAL,INTENT(IN):: Force                !< If present and true, allocate Force field
      LOGICAL,OPTIONAL,INTENT(IN):: Moment               !< If present and true, allocate Moment field
      LOGICAL,OPTIONAL,INTENT(IN):: Orientation          !< If present and true, allocate Orientation field
      LOGICAL,OPTIONAL,INTENT(IN):: TranslationDisp      !< If present and true, allocate TranslationDisp field
      LOGICAL,OPTIONAL,INTENT(IN):: TranslationVel       !< If present and true, allocate TranslationVel field
      LOGICAL,OPTIONAL,INTENT(IN):: RotationVel          !< If present and true, allocate RotationVel field
      LOGICAL,OPTIONAL,INTENT(IN):: TranslationAcc       !< If present and true, allocate TranslationAcc field
      LOGICAL,OPTIONAL,INTENT(IN):: RotationAcc          !< If present and true, allocate RotationAcc field
!
      INTEGER,OPTIONAL,INTENT(IN):: nScalars             !< If present and > 0, allocate nScalars Scalars
      LOGICAL,OPTIONAL,INTENT(IN):: IsNewSibling         !< If present and true, this is an new sibling so don't allocate new shared fields (RemapFlag, position, RefOrientation, and ElemTable)

    ! Local
      INTEGER i
      LOGICAL                    :: IsNewSib

      LOGICAL                    :: IsMotion
      LOGICAL                    :: IsLoad
      INTEGER(IntKi)             :: ErrStat2 
      CHARACTER(ErrMsgLen)       :: ErrMess2
      CHARACTER(*),PARAMETER     :: RoutineName = 'MeshCreate'

         ! Local initializations:

      ErrStat = ErrID_None
      ErrMess = ""
      IsMotion = .FALSE.
      IsLoad   = .FALSE.

      IF ( mesh_debug ) print*,'Called MeshCreate'

      CALL MeshDestroy( BlankMesh, ErrStat2, ErrMess2, .TRUE. )
                                                        ! make sure we're not leaving any pointers dangling
                                                        ! and nullify them for good measure
                                                        ! See comment on optional IgnoreSibling argument
                                                        ! in definition of MeshDestroy
      CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName) 
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
         CALL AllocPAry( BlankMesh%Position, 3, Nnodes, 'MeshCreate: Position', ErrStat2, ErrMess2 )
            CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
         CALL AllocPAry( BlankMesh%RefOrientation, 3, 3, Nnodes, 'MeshCreate: RefOrientation', ErrStat2, ErrMess2 )
            CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
            IF (ErrStat >= AbortErrLev) RETURN
            ! initialize these variables:
            BlankMesh%Position = 0.0_ReKi
            CALL Eye(BlankMesh%RefOrientation, ErrStat2, ErrMess2)
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)     
            
         ALLOCATE(BlankMesh%ElemTable(NELEMKINDS),STAT=ErrStat2)
         IF (ErrStat2/=0) THEN
            CALL SetErrStat(ErrID_Fatal, "Error allocating ElemTable.", ErrStat, ErrMess,RoutineName)     
            RETURN
         END IF
         
            
         DO i = 1, NELEMKINDS
            BlankMesh%ElemTable(i)%nelem = 0  ; BlankMesh%ElemTable(i)%maxelem = 0
            NULLIFY(BlankMesh%ElemTable(i)%Elements )
         ENDDO

         ALLOCATE(BlankMesh%RemapFlag, Stat=ErrStat2 ) ! assign some space for this pointer to point to
         IF (ErrStat2/=0) THEN
            CALL SetErrStat(ErrID_Fatal, "Error allocating RemapFlag.", ErrStat, ErrMess,RoutineName)     
            RETURN
         END IF
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
            CALL AllocAry( BlankMesh%Force, 3, Nnodes, 'MeshCreate: Force', ErrStat2, ErrMess2 )
            CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%Force = 0.
            BlankMesh%FieldMask(MASKID_FORCE) = .TRUE.
            IsLoad = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(Moment) ) THEN
         IF ( Moment ) THEN
            CALL AllocAry( BlankMesh%Moment, 3, Nnodes, 'MeshCreate: Moment', ErrStat2, ErrMess2 )
            CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
            IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%Moment = 0.
            BlankMesh%FieldMask(MASKID_MOMENT) = .TRUE.
            IsLoad = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(Orientation) ) THEN
         IF ( Orientation ) THEN
            CALL AllocAry( BlankMesh%Orientation, 3, 3, Nnodes, 'MeshCreate: Orientation', ErrStat2,ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            CALL Eye(BlankMesh%Orientation, ErrStat2, ErrMess2)  ! set this orientation to the identity matrix
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
            BlankMesh%FieldMask(MASKID_ORIENTATION) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(TranslationDisp) ) THEN
         IF ( TranslationDisp ) THEN
            CALL AllocAry( BlankMesh%TranslationDisp, 3, Nnodes, 'MeshCreate: TranslationDisp', ErrStat2, ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationDisp = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONDISP) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(TranslationVel) ) THEN
         IF ( TranslationVel ) THEN
            CALL AllocAry( BlankMesh%TranslationVel, 3, Nnodes, 'MeshCreate: TranslationVel', ErrStat2, ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationVel = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONVEL) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(RotationVel) ) THEN
         IF ( RotationVel ) THEN
            CALL AllocAry( BlankMesh%RotationVel, 3, Nnodes, 'MeshCreate: RotationVel', ErrStat2, ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%RotationVel = 0.
            BlankMesh%FieldMask(MASKID_ROTATIONVEL) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(TranslationAcc) ) THEN
         IF ( TranslationAcc ) THEN
            CALL AllocAry( BlankMesh%TranslationAcc, 3, Nnodes, 'MeshCreate: TranslationAcc', ErrStat2, ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationAcc = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONACC) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF

      IF ( PRESENT(RotationAcc) ) THEN
         IF ( RotationAcc ) THEN
            CALL AllocAry( BlankMesh%RotationAcc, 3, Nnodes, 'MeshCreate: RotationAcc', ErrStat2, ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%RotationAcc = 0.
            BlankMesh%FieldMask(MASKID_ROTATIONACC) = .TRUE.
            IsMotion = .TRUE.
         ENDIF
      ENDIF


      BlankMesh%nScalars = 0
      IF ( PRESENT(nScalars) ) THEN
         IF ( nScalars .GT. 0 ) THEN
            CALL AllocAry( BlankMesh%Scalars, nScalars, Nnodes, 'MeshCreate: Scalars', ErrStat2,ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
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
      !> This routine will add any required fields that were not explicitly requested.
      !! If the mesh has motion fields and it is an input mesh, it must always have the following fields: TranslationDisp, TranslationVel, TranslationAcc
      IF ( IsMotion .AND. IOS == COMPONENT_INPUT ) THEN

         IF ( .NOT. BlankMesh%FieldMask(MASKID_TRANSLATIONDISP)) THEN
            CALL AllocAry( BlankMesh%TranslationDisp, 3, Nnodes, 'MeshCreate: TranslationDisp', ErrStat2, ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationDisp = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONDISP) = .TRUE.
            !CALL SetErrStat(ErrID_Info, 'Meshes with motion fields must also contain the TranslationDisp field.',ErrStat,ErrMsg,'MeshCreate')
         ENDIF

         IF ( .NOT. BlankMesh%FieldMask(MASKID_TRANSLATIONVEL)) THEN
            CALL AllocAry( BlankMesh%TranslationVel, 3, Nnodes, 'MeshCreate: TranslationVel', ErrStat2, ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationVel = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONVEL) = .TRUE.
            !CALL SetErrStat(ErrID_Info, 'Meshes with motion fields must also contain the TranslationVel field.',ErrStat,ErrMsg,'MeshCreate')
         ENDIF
         
         IF ( .NOT. BlankMesh%FieldMask(MASKID_TRANSLATIONACC)) THEN
            CALL AllocAry( BlankMesh%TranslationAcc, 3, Nnodes, 'MeshCreate: TranslationAcc', ErrStat2, ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%TranslationAcc = 0.
            BlankMesh%FieldMask(MASKID_TRANSLATIONACC) = .TRUE.
            !CALL SetErrStat(ErrID_Info, 'Meshes with motion fields must also contain the TranslationAcc field.',ErrStat,ErrMsg,'MeshCreate')
         ENDIF
                           
      END IF
      
      !> If the mesh has load fields and it is an input mesh, it must always have the following fields: Moment
      IF ( IsLoad .AND. IOS == COMPONENT_INPUT ) THEN
               
         IF ( .NOT. BlankMesh%FieldMask(MASKID_MOMENT)) THEN
            CALL AllocAry( BlankMesh%Moment, 3, Nnodes, 'MeshCreate: Moment', ErrStat2, ErrMess2 )
               CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess,RoutineName)          
               IF (ErrStat >= AbortErrLev) RETURN
            BlankMesh%Moment = 0.
            BlankMesh%FieldMask(MASKID_MOMENT) = .TRUE.
            !CALL SetErrStat(ErrID_Info, 'Meshes with load fields must also contain the Moment field.',ErrStat,ErrMsg,'MeshCreate')
         ENDIF               
         
      END IF
      

      RETURN

   END SUBROUTINE MeshCreate

!> Destroy the given mesh and deallocate all of its data. If the optional IgnoreSibling argument 
!! is set to TRUE, destroying a sibling in a set has no effect on the other siblings other than 
!! to remove the victim from the list of siblings. If IgnoreSibling is omitted or is set to FALSE, 
!! all of the other siblings in the set will be destroyed as well.
   RECURSIVE SUBROUTINE MeshDestroy ( Mesh, ErrStat, ErrMess, IgnoreSibling )

     TYPE(MeshType),  INTENT(INOUT) :: Mesh            !< Mesh to be vaporized
     INTEGER(IntKi),  INTENT(OUT)   :: ErrStat         !< Error status/code
     CHARACTER(*),    INTENT(OUT)   :: ErrMess         !< Error message
    ! On a brand new mesh, the pointers to siblings may not be nullified and
    ! thus undefined (which may cause ASSOCIATED to report .true. erroneously)
    ! This despite use of => NULL in declaration of this fields for MeshType. Sigh. (15-dec-2015 bjj: not sure this is true; some fields didn't use => NULL)
    ! So ...
     LOGICAL, INTENT(IN), OPTIONAL :: IgnoreSibling    !< if IgnoreSibling is present and true, don't follow the sibling pointers.
                                                       !! Instead just unconditionally nullify these.
                                                       !! Use this carefully, since it can leave dangling memory if used for a
                                                       !! mesh that already exists and has existing siblings. 

    ! Local
      LOGICAL IgSib
      INTEGER i, j

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
                     IF (ALLOCATED(Mesh%ElemTable(i)%Elements(j)%ElemNodes)) THEN
                        DEALLOCATE(Mesh%ElemTable(i)%Elements(j)%ElemNodes)
                     ENDIF
                     !IF (ASSOCIATED(Mesh%ElemTable(i)%Elements(j)%Neighbors)) THEN
                     !   DEALLOCATE(Mesh%ElemTable(i)%Elements(j)%Neighbors)
                     !   NULLIFY(Mesh%ElemTable(i)%Elements(j)%Neighbors)
                     !ENDIF
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

!----------------------------------------------------------------------------------------------------------------------------------
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
!> Given a mesh and allocatable buffers of type INTEGER(IntKi), REAL(ReKi), and REAL(DbKi), 
!! return the mesh information compacted into consecutive elements of the corresponding buffers. 
!! This would be done to allow subsequent writing of the buffers to a file for restarting later. 
!! The sense of the name is "pack the data from the mesh into buffers". IMPORTANT: MeshPack 
!! allocates the three buffers. It is incumbent upon the calling program to deallocate the 
!! buffers when they are no longer needed. For sibling meshes, MeshPack should be called 
!! separately for each sibling, because the fields allocated with the siblings are separate 
!! and unique to each sibling.
   SUBROUTINE MeshPack ( Mesh, ReKiBuf, DbKiBuf, IntKiBuf , ErrStat, ErrMess, SizeOnly )
   
     TYPE(MeshType),              INTENT(IN   ) :: Mesh        ! Mesh being packed
     REAL(ReKi),     ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)  ! Real buffer
     REAL(DbKi),     ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)  ! Double buffer
     INTEGER(IntKi), ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:) ! Int buffer
     INTEGER(IntKi),              INTENT(  OUT) :: ErrStat
     CHARACTER(*),                INTENT(  OUT) :: ErrMess
     LOGICAL,OPTIONAL,            INTENT(IN   ) :: SizeOnly
     
   ! Local
     INTEGER(IntKi)                             :: Re_BufSz      ! number of reals in the buffer
     INTEGER(IntKi)                             :: Re_Xferred    ! number of reals transferred
     INTEGER(IntKi)                             :: Db_BufSz      ! number of doubles in the buffer
     INTEGER(IntKi)                             :: Db_Xferred    ! number of doubles transferred
     INTEGER(IntKi)                             :: Int_BufSz     ! number of integers in the buffer
     INTEGER(IntKi)                             :: Int_Xferred   ! number of integers transferred
   
   
     INTEGER i,j, nelemnodes
     LOGICAL OnlySize
     INTEGER(IntKi)                             :: ErrStat2
     !CHARACTER(1024)                            :: ErrMess2
     CHARACTER(*),      PARAMETER               :: RoutineName = "MeshPack"   

     
     ErrStat = ErrID_None
     ErrMess = ""
     
     OnlySize = .FALSE.
     IF ( PRESENT(SizeOnly) ) OnlySize = SizeOnly


     ! bjj: figure out what to do about sibling meshes... (for now, I'm going to ignore them)
     
     !.........................................
     ! get number of integer values
     !.........................................
      IF (.NOT. Mesh%Initialized) THEN ! we don't need to store any data; it's a blank mesh
         Int_BufSz = 1
      ELSE ! initialized, may or may not be committed           
         Int_BufSz =  3                & ! number of logicals in MeshType (initialized, committed, RemapFlag)
                     + FIELDMASK_SIZE  & ! number of logicals in MeshType (fieldmask)
                     + 5                 ! number of non-pointer integers (ios, nnodes, nextelem, nscalars, refNode)
         
         !......
         ! we'll store the element structure (and call MeshCommit on Unpack if necessary to get the remaining fields like det_jac)
         !......
         DO i = 1, NELEMKINDS

            Int_BufSz = Int_BufSz+1 ! Mesh%ElemTable(i)%nelem
            if (Mesh%ElemTable(i)%nelem > 0) Int_BufSz = Int_BufSz+1 ! number of nodes in this kind of element            
                  
            DO j = 1, Mesh%ElemTable(i)%nelem            
               !Int_BufSz = Int_BufSz+1 ! which kind of element 
               !Int_BufSz = Int_BufSz+1 ! skip Nneighbors until that's implemented (as well as neighbor list)
               Int_BufSz = Int_BufSz + SIZE( Mesh%ElemTable(i)%Elements(j)%ElemNodes ) ! nodes in this element                     
            END DO
         
         END DO             
         
      END IF
      
     !.........................................
     ! get number of real values
     !.........................................
     Re_BufSz = 0
     IF (Mesh%Initialized) THEN
        Re_BufSz = Re_BufSz + Mesh%Nnodes * 3 ! Position
        !Re_BufSz = Re_BufSz + Mesh%Nnodes * 9 ! RefOrientation
        IF ( Mesh%FieldMask(MASKID_FORCE) ) Re_BufSz = Re_BufSz + Mesh%Nnodes * 3
        IF ( Mesh%FieldMask(MASKID_MOMENT) ) Re_BufSz = Re_BufSz + Mesh%Nnodes * 3
        !IF ( Mesh%FieldMask(MASKID_ORIENTATION) ) Re_BufSz = Re_BufSz + Mesh%Nnodes * 9
        !IF ( Mesh%FieldMask(MASKID_TRANSLATIONDISP) ) Re_BufSz = Re_BufSz + Mesh%Nnodes * 3
        IF ( Mesh%FieldMask(MASKID_ROTATIONVEL) ) Re_BufSz = Re_BufSz + Mesh%Nnodes * 3
        IF ( Mesh%FieldMask(MASKID_TRANSLATIONVEL) ) Re_BufSz = Re_BufSz + Mesh%Nnodes * 3
        IF ( Mesh%FieldMask(MASKID_ROTATIONACC) ) Re_BufSz = Re_BufSz + Mesh%Nnodes * 3
        IF ( Mesh%FieldMask(MASKID_TRANSLATIONACC) ) Re_BufSz = Re_BufSz + Mesh%Nnodes * 3
        IF ( Mesh%nScalars .GT. 0 ) Re_BufSz = Re_BufSz + Mesh%Nnodes * Mesh%nScalars
     END IF
     
     !.........................................
     ! get number of double values (none now)
     !.........................................
     Db_BufSz = 0
     IF (Mesh%Initialized) THEN
        Db_BufSz = Db_BufSz + Mesh%Nnodes * 9 ! RefOrientation
        IF ( Mesh%FieldMask(MASKID_ORIENTATION) ) Db_BufSz = Db_BufSz + Mesh%Nnodes * 9
        IF ( Mesh%FieldMask(MASKID_TRANSLATIONDISP) ) Db_BufSz = Db_BufSz + Mesh%Nnodes * 3
     END IF 
     
     !.........................................
     ! allocate buffer arrays
     !.........................................          
     IF ( Re_BufSz  .GT. 0 ) THEN 
        ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )
        IF (ErrStat2 /= 0) THEN 
          CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMess,RoutineName)
          RETURN
        END IF
     END IF
     IF ( Db_BufSz  .GT. 0 ) THEN 
        ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )
        IF (ErrStat2 /= 0) THEN 
          CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMess,RoutineName)
          RETURN
        END IF
     END IF
     IF ( Int_BufSz  .GT. 0 ) THEN 
        ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )
        IF (ErrStat2 /= 0) THEN 
          CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMess,RoutineName)
          RETURN
        END IF
     END IF
     IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)

     
     !.........................................
     ! store data in buffer arrays
     !.........................................     
     Re_Xferred  = 1
     Db_Xferred  = 1
     Int_Xferred = 1        
     
     ! ..... fill IntKiBuf .....
     
      IF (.NOT. Mesh%Initialized) THEN ! we don't need to store any data; it's a blank mesh
         IntKiBuf(Int_Xferred) = 0; ;  Int_Xferred = Int_Xferred + 1
      ELSE ! initialized, may or may not be committed                    
            ! transfer the logicals
         IntKiBuf(Int_Xferred) = 1;   Int_Xferred = Int_Xferred + 1 
         IntKiBuf(Int_Xferred) = TRANSFER( Mesh%committed, IntKiBuf(1) );  Int_Xferred = Int_Xferred + 1 
         IntKiBuf(Int_Xferred:Int_Xferred+FIELDMASK_SIZE-1) = TRANSFER( Mesh%fieldmask, IntKiBuf(Int_Xferred:Int_Xferred+FIELDMASK_SIZE-1) );  Int_Xferred = Int_Xferred + FIELDMASK_SIZE
         IntKiBuf(Int_Xferred) = TRANSFER( Mesh%RemapFlag, IntKiBuf(1) );  Int_Xferred = Int_Xferred + 1        
            ! integers
         IntKiBuf(Int_Xferred) = Mesh%ios;           Int_Xferred = Int_Xferred + 1 
         IntKiBuf(Int_Xferred) = Mesh%nnodes;        Int_Xferred = Int_Xferred + 1 
         IntKiBuf(Int_Xferred) = Mesh%refnode;       Int_Xferred = Int_Xferred + 1 
         IntKiBuf(Int_Xferred) = Mesh%nextelem;      Int_Xferred = Int_Xferred + 1 
         IntKiBuf(Int_Xferred) = Mesh%nscalars;      Int_Xferred = Int_Xferred + 1 
         
            ! element structure
         DO i = 1, NELEMKINDS
            
            IntKiBuf(Int_Xferred) = Mesh%ElemTable(i)%nelem;      Int_Xferred = Int_Xferred + 1 ! number of elements
            
            if (Mesh%ElemTable(i)%nelem > 0) then
               nelemnodes = SIZE( Mesh%ElemTable(i)%Elements(1)%ElemNodes );
               IntKiBuf(Int_Xferred) = nelemnodes; Int_Xferred = Int_Xferred + 1 ! nodes per element
                                          
                  ! nodes in this element
               DO j = 1, Mesh%ElemTable(i)%nelem            
                  IntKiBuf(Int_Xferred:Int_Xferred+nelemnodes-1) = Mesh%ElemTable(i)%Elements(j)%ElemNodes; Int_Xferred = Int_Xferred + nelemnodes 
               END DO
            end if
         
         END DO             
         
      END IF         
     
     ! ..... fill ReKiBuf and DbKiBuf .....
     IF (Mesh%Initialized) THEN
         DO i = 1, Mesh%Nnodes ! Position
            ReKiBuf(Re_Xferred:Re_Xferred+2) = Mesh%Position(:,i); Re_Xferred = Re_Xferred + 3
         END DO
         DO i = 1, Mesh%Nnodes ! RefOrientation
            DO j = 1,3
               DbKiBuf(Db_Xferred:Db_Xferred+2) = Mesh%RefOrientation(:,j,i); Db_Xferred = Db_Xferred + 3
            ENDDO
         END DO
                
         IF ( Mesh%FieldMask(MASKID_FORCE) ) THEN ! Force
            DO i = 1, Mesh%Nnodes
               ReKiBuf(Re_Xferred:Re_Xferred+2) = Mesh%Force(:,i); Re_Xferred = Re_Xferred + 3
            ENDDO
         ENDIF
         IF ( Mesh%FieldMask(MASKID_MOMENT) ) THEN ! Moment
            DO i = 1, Mesh%Nnodes
               ReKiBuf(Re_Xferred:Re_Xferred+2) = Mesh%Moment(:,i); Re_Xferred = Re_Xferred + 3
            ENDDO
         ENDIF
         IF ( Mesh%FieldMask(MASKID_ORIENTATION) ) THEN ! Orientation
            DO i = 1, Mesh%Nnodes 
               DO j = 1,3
                  DbKiBuf(Db_Xferred:Db_Xferred+2) = Mesh%Orientation(:,j,i); Db_Xferred = Db_Xferred + 3
               ENDDO
            END DO
         ENDIF
         IF ( Mesh%FieldMask(MASKID_TRANSLATIONDISP) ) THEN ! TranslationDisp
            DO i = 1, Mesh%Nnodes
               DbKiBuf(Db_Xferred:Db_Xferred+2) = Mesh%TranslationDisp(:,i); Db_Xferred = Db_Xferred + 3
            ENDDO
         ENDIF
         IF ( Mesh%FieldMask(MASKID_ROTATIONVEL) ) THEN ! RotationVel
            DO i = 1, Mesh%Nnodes
               ReKiBuf(Re_Xferred:Re_Xferred+2) = Mesh%RotationVel(:,i); Re_Xferred = Re_Xferred + 3
            ENDDO
         ENDIF
         IF ( Mesh%FieldMask(MASKID_TRANSLATIONVEL) ) THEN ! TranslationVel
            DO i = 1, Mesh%Nnodes
               ReKiBuf(Re_Xferred:Re_Xferred+2) = Mesh%TranslationVel(:,i); Re_Xferred = Re_Xferred + 3
            ENDDO
         ENDIF
         IF ( Mesh%FieldMask(MASKID_ROTATIONACC) ) THEN ! RotationAcc
            DO i = 1, Mesh%Nnodes
               ReKiBuf(Re_Xferred:Re_Xferred+2) = Mesh%RotationAcc(:,i); Re_Xferred = Re_Xferred + 3
            ENDDO
         ENDIF
         IF ( Mesh%FieldMask(MASKID_TRANSLATIONACC) ) THEN ! TranslationAcc
            DO i = 1, Mesh%Nnodes
               ReKiBuf(Re_Xferred:Re_Xferred+2) = Mesh%TranslationAcc(:,i); Re_Xferred = Re_Xferred + 3
            ENDDO
         ENDIF
         
         IF ( Mesh%nScalars .GT. 0 ) THEN ! n_re = n_re + Mesh%Nnodes * Mesh%nScalar
            DO i = 1, Mesh%Nnodes
               ReKiBuf(Re_Xferred:Re_Xferred+Mesh%nScalars-1) = Mesh%Scalars(:,i); Re_Xferred = Re_Xferred + Mesh%nScalars
            ENDDO
         ENDIF         
         
     END IF
     
     !bjj: where are we keeping track of which ones are siblings so that we can unpack them (set pointers) properly for restart?
   END SUBROUTINE MeshPack

!----------------------------------------------------------------------------------------------------------------------------------
!> Given a blank, uncreated mesh and buffers of type INTEGER(IntKi), REAL(ReKi), and 
!! REAL(DbKi), unpack the mesh information from the buffers. This would be done to 
!! recreate a mesh after reading in the buffers on a restart of the program. The sense 
!! of the name is "unpack the mesh from buffers." The resulting mesh will be returned 
!! in the exact state as when the data in the buffers was packed using MeshPack. 
   SUBROUTINE MeshUnpack( Mesh, ReKiBuf, DbKiBuf, IntKiBuf, ErrStat, ErrMess )
      
      ! bjj: not implemented yet:  
      ! If the mesh has an already recreated sibling mesh from a previous call to MeshUnpack, specify 
      ! the existing sibling as an optional argument so that the sibling relationship is also recreated.
   
      TYPE(MeshType),              INTENT(INOUT) :: Mesh
      REAL(ReKi),     ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)
      REAL(DbKi),     ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)
      INTEGER(IntKi), ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)
      INTEGER(IntKi),              INTENT(  OUT) :: ErrStat
      CHARACTER(*),                INTENT(  OUT) :: ErrMess

         ! Local
      LOGICAL committed, RemapFlag, fieldmask(FIELDMASK_SIZE)
      INTEGER nScalars, ios, nnodes, nextelem, nelemnodes, nelem, refnode
      INTEGER i,j
     
      INTEGER(IntKi)                             :: Re_Xferred    ! number of reals transferred
      INTEGER(IntKi)                             :: Db_Xferred    ! number of doubles transferred
      INTEGER(IntKi)                             :: Int_Xferred   ! number of integers transferred
      
      INTEGER(IntKi)                             :: ErrStat2
      CHARACTER(ErrMsgLen)                       :: ErrMess2
      CHARACTER(*),      PARAMETER               :: RoutineName = "MeshUnpack"        
     
      Re_Xferred  = 1
      Db_Xferred  = 1
      Int_Xferred = 1

      ErrStat = ErrID_None
      ErrMess = ""
     
      IF (IntKiBuf(Int_Xferred) == 0 ) THEN ! this is a blank mesh
         CALL MeshDestroy( Mesh, ErrStat2, ErrMess2, .TRUE. )
         CALL SetErrStat(ErrStat2,ErrMess2,ErrStat,ErrMess,RoutineName)
         RETURN
      END IF
      
      
         ! initialized, may or may not be committed           
            
      Mesh%initialized = .true.; Int_Xferred = Int_Xferred + 1
      committed        = TRANSFER( IntKiBuf(Int_Xferred), Mesh%committed );  Int_Xferred = Int_Xferred + 1 
      fieldmask        = TRANSFER( IntKiBuf(Int_Xferred:Int_Xferred+FIELDMASK_SIZE-1), fieldmask );  Int_Xferred = Int_Xferred + FIELDMASK_SIZE 
      RemapFlag        = TRANSFER( IntKiBuf(Int_Xferred), Mesh%RemapFlag ); Int_Xferred = Int_Xferred + 1 
         ! integers
      ios              = IntKiBuf(Int_Xferred) ; Int_Xferred = Int_Xferred + 1 
      nnodes           = IntKiBuf(Int_Xferred) ; Int_Xferred = Int_Xferred + 1 
      refnode          = IntKiBuf(Int_Xferred) ; Int_Xferred = Int_Xferred + 1 
      nextelem         = IntKiBuf(Int_Xferred) ; Int_Xferred = Int_Xferred + 1 
      nscalars         = IntKiBuf(Int_Xferred) ; Int_Xferred = Int_Xferred + 1 
                  
                                                            
      CALL MeshCreate( Mesh, ios, nnodes                                         &
                        ,ErrStat=ErrStat2, ErrMess=ErrMess2                      &
                        ,Force          =fieldmask(MASKID_FORCE)                 &
                        ,Moment         =fieldmask(MASKID_MOMENT)                &
                        ,Orientation    =fieldmask(MASKID_ORIENTATION)           &
                        ,TranslationDisp=fieldmask(MASKID_TRANSLATIONDISP)       &
                        ,TranslationVel =fieldmask(MASKID_TRANSLATIONVEL )       &
                        ,RotationVel    =fieldmask(MASKID_ROTATIONVEL )          &
                        ,TranslationAcc =fieldmask(MASKID_TRANSLATIONACC )       &
                        ,RotationAcc    =fieldmask(MASKID_ROTATIONACC )          &
                        ,nScalars       = nScalars                               &
                     )
         CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess, RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN

      Mesh%RefNode = refnode
      Mesh%RemapFlag = RemapFlag
      Mesh%nextelem  = nextelem
     
            ! element structure
      DO i = 1, NELEMKINDS            
         nelem = IntKiBuf(Int_Xferred);      Int_Xferred = Int_Xferred + 1    ! number of elements            
            
         if (nelem > 0) then
            nelemnodes = IntKiBuf(Int_Xferred); Int_Xferred = Int_Xferred + 1 ! nodes per element
                                                                             
               ! nodes in this element
            DO j = 1,nelem                                                                     
               
               SELECT CASE (nelemnodes)
               CASE (1)
                  CALL MeshConstructElement( Mesh, i, ErrStat, ErrMess                                      &
                  , P1 =IntKiBuf(Int_Xferred  )                                                             &
                  )
               CASE (2)
                  CALL MeshConstructElement( Mesh, i, ErrStat, ErrMess                                      &
                  , P1 =IntKiBuf(Int_Xferred   ),P2 =IntKiBuf(Int_Xferred+ 1)                               &
                  )
               CASE (3)
                  CALL MeshConstructElement( Mesh, i, ErrStat, ErrMess                                      &
                  , P1 =IntKiBuf(Int_Xferred   ),P2 =IntKiBuf(Int_Xferred+ 1),P3 =IntKiBuf(Int_Xferred+ 2)  &
                  )
               CASE (4)
                  CALL MeshConstructElement( Mesh, i, ErrStat, ErrMess                                      &
                  , P1 =IntKiBuf(Int_Xferred   ),P2 =IntKiBuf(Int_Xferred+ 1),P3 =IntKiBuf(Int_Xferred+ 2)  &
                  , P4 =IntKiBuf(Int_Xferred+ 3)                                                            &
                  )
               CASE (6)
                  CALL MeshConstructElement( Mesh, i, ErrStat, ErrMess                                      &
                  , P1 =IntKiBuf(Int_Xferred   ),P2 =IntKiBuf(Int_Xferred+ 1),P3 =IntKiBuf(Int_Xferred+ 2)  &
                  , P4 =IntKiBuf(Int_Xferred+ 3),P5 =IntKiBuf(Int_Xferred+ 4),P6 =IntKiBuf(Int_Xferred+ 5)  &
                  )
               CASE (8)
                  CALL MeshConstructElement( Mesh, i, ErrStat2, ErrMess2                                    &
                  , P1 =IntKiBuf(Int_Xferred   ),P2 =IntKiBuf(Int_Xferred+ 1),P3 =IntKiBuf(Int_Xferred+ 2)  &
                  , P4 =IntKiBuf(Int_Xferred+ 3),P5 =IntKiBuf(Int_Xferred+ 4),P6 =IntKiBuf(Int_Xferred+ 5)  &
                  , P7 =IntKiBuf(Int_Xferred+ 6),P8 =IntKiBuf(Int_Xferred+ 7)                               &
                  )
               CASE (10)
                  CALL MeshConstructElement( Mesh, i, ErrStat, ErrMess                                      &
                  , P1 =IntKiBuf(Int_Xferred   ),P2 =IntKiBuf(Int_Xferred+ 1),P3 =IntKiBuf(Int_Xferred+ 2)  &
                  , P4 =IntKiBuf(Int_Xferred+ 3),P5 =IntKiBuf(Int_Xferred+ 4),P6 =IntKiBuf(Int_Xferred+ 5)  &
                  , P7 =IntKiBuf(Int_Xferred+ 6),P8 =IntKiBuf(Int_Xferred+ 7),P9 =IntKiBuf(Int_Xferred+ 8)  &
                  , P10=IntKiBuf(Int_Xferred+ 9)                                                            &
                  )
               CASE (15)
                  CALL MeshConstructElement( Mesh, i, ErrStat, ErrMess                                      &
                  , P1 =IntKiBuf(Int_Xferred   ),P2 =IntKiBuf(Int_Xferred+ 1),P3 =IntKiBuf(Int_Xferred+ 2)  &
                  , P4 =IntKiBuf(Int_Xferred+ 3),P5 =IntKiBuf(Int_Xferred+ 4),P6 =IntKiBuf(Int_Xferred+ 5)  &
                  , P7 =IntKiBuf(Int_Xferred+ 6),P8 =IntKiBuf(Int_Xferred+ 7),P9 =IntKiBuf(Int_Xferred+ 8)  &
                  , P10=IntKiBuf(Int_Xferred+ 9),P11=IntKiBuf(Int_Xferred+10),P12=IntKiBuf(Int_Xferred+11)  &
                  , P13=IntKiBuf(Int_Xferred+12),P14=IntKiBuf(Int_Xferred+13),P15=IntKiBuf(Int_Xferred+14)  &
                  )
               CASE (20)
                  CALL MeshConstructElement( Mesh, i, ErrStat, ErrMess                                      &
                  , P1 =IntKiBuf(Int_Xferred   ),P2 =IntKiBuf(Int_Xferred+ 1),P3 =IntKiBuf(Int_Xferred+ 2)  &
                  , P4 =IntKiBuf(Int_Xferred+ 3),P5 =IntKiBuf(Int_Xferred+ 4),P6 =IntKiBuf(Int_Xferred+ 5)  &
                  , P7 =IntKiBuf(Int_Xferred+ 6),P8 =IntKiBuf(Int_Xferred+ 7),P9 =IntKiBuf(Int_Xferred+ 8)  &
                  , P10=IntKiBuf(Int_Xferred+ 9),P11=IntKiBuf(Int_Xferred+10),P12=IntKiBuf(Int_Xferred+11)  &
                  , P13=IntKiBuf(Int_Xferred+12),P14=IntKiBuf(Int_Xferred+13),P15=IntKiBuf(Int_Xferred+14)  &
                  , P16=IntKiBuf(Int_Xferred+15),P17=IntKiBuf(Int_Xferred+16),P18=IntKiBuf(Int_Xferred+17)  &
                  , P19=IntKiBuf(Int_Xferred+18),P20=IntKiBuf(Int_Xferred+19)                               &
                  )
               CASE DEFAULT
                  CALL SetErrStat(ErrID_Fatal,"No such element. Probably manged buffer.",ErrStat,ErrMess,RoutineName)
                  RETURN
               END SELECT  
               Int_Xferred = Int_Xferred + nelemnodes
            END DO   ! Elements of this kind
         end if ! if there are any elements of this kind

      END DO ! kinds of elements
      
     ! ..... fill ReKiBuf .....
      DO i = 1, Mesh%Nnodes ! Position
         Mesh%Position(:,i) = ReKiBuf(Re_Xferred:Re_Xferred+2); Re_Xferred = Re_Xferred + 3
      END DO
      DO i = 1, Mesh%Nnodes ! RefOrientation
         DO j = 1,3
            Mesh%RefOrientation(:,j,i) = DbKiBuf(Db_Xferred:Db_Xferred+2); Db_Xferred = Db_Xferred + 3
         ENDDO
      END DO
                
      IF ( FieldMask(MASKID_FORCE) ) THEN ! Force
         DO i = 1, Mesh%Nnodes
            Mesh%Force(:,i) = ReKiBuf(Re_Xferred:Re_Xferred+2); Re_Xferred = Re_Xferred + 3
         ENDDO
      ENDIF
      IF ( FieldMask(MASKID_MOMENT) ) THEN ! Moment
         DO i = 1, Mesh%Nnodes
            Mesh%Moment(:,i) = ReKiBuf(Re_Xferred:Re_Xferred+2); Re_Xferred = Re_Xferred + 3
         ENDDO
      ENDIF
      IF ( FieldMask(MASKID_ORIENTATION) ) THEN ! Orientation
         DO i = 1, Mesh%Nnodes 
            DO j = 1,3
               Mesh%Orientation(:,j,i) = DbKiBuf(Db_Xferred:Db_Xferred+2); Db_Xferred = Db_Xferred + 3
            ENDDO
         END DO
      ENDIF
      IF ( FieldMask(MASKID_TRANSLATIONDISP) ) THEN ! TranslationDisp
         DO i = 1, Mesh%Nnodes
            Mesh%TranslationDisp(:,i) = DbKiBuf(Db_Xferred:Db_Xferred+2); Db_Xferred = Db_Xferred + 3
         ENDDO
      ENDIF
      IF ( FieldMask(MASKID_ROTATIONVEL) ) THEN ! RotationVel
         DO i = 1, Mesh%Nnodes
            Mesh%RotationVel(:,i) = ReKiBuf(Re_Xferred:Re_Xferred+2); Re_Xferred = Re_Xferred + 3
         ENDDO
      ENDIF
      IF ( FieldMask(MASKID_TRANSLATIONVEL) ) THEN ! TranslationVel
         DO i = 1, Mesh%Nnodes
            Mesh%TranslationVel(:,i) = ReKiBuf(Re_Xferred:Re_Xferred+2); Re_Xferred = Re_Xferred + 3
         ENDDO
      ENDIF
      IF ( FieldMask(MASKID_ROTATIONACC) ) THEN ! RotationAcc
         DO i = 1, Mesh%Nnodes
            Mesh%RotationAcc(:,i) = ReKiBuf(Re_Xferred:Re_Xferred+2); Re_Xferred = Re_Xferred + 3
         ENDDO
      ENDIF
      IF ( FieldMask(MASKID_TRANSLATIONACC) ) THEN ! TranslationAcc
         DO i = 1, Mesh%Nnodes
            Mesh%TranslationAcc(:,i) = ReKiBuf(Re_Xferred:Re_Xferred+2); Re_Xferred = Re_Xferred + 3
         ENDDO
      ENDIF
         
      IF ( Mesh%nScalars .GT. 0 ) THEN ! n_re = n_re + Mesh%Nnodes * Mesh%nScalar
         DO i = 1, Mesh%Nnodes
            Mesh%Scalars(:,i) = ReKiBuf(Re_Xferred:Re_Xferred+Mesh%nScalars-1); Re_Xferred = Re_Xferred + Mesh%nScalars
         ENDDO
      ENDIF         
                     
      ! commit the mesh
      IF (committed) THEN
         CALL MeshCommit(Mesh, ErrStat2, ErrMess2)
         CALL SetErrStat(ErrStat2, ErrMess2, ErrStat, ErrMess, RoutineName)
      END IF
      
      RETURN

   END SUBROUTINE MeshUnpack

!----------------------------------------------------------------------------------------------------------------------------------
!> Given an existing mesh and a destination mesh, create a completely new copy, a sibling, or 
!!   update the fields of a second existing mesh from the first mesh. When CtrlCode is 
!!   MESH_NEWCOPY or MESH_SIBLING, the destination mesh must be a blank, uncreated mesh.
!! 
!! If CtrlCode is MESH_NEWCOPY, an entirely new copy of the mesh is created, including all fields, 
!!   with the same data values as the original, but as an entirely separate copy in memory. The new 
!!   copy is in the same state as the original--if the original has not been committed, neither is 
!!   the copy; in this case, an all-new copy of the mesh must be committed separately.
!!
!! If CtrlCode is MESH_SIBLING, the destination mesh is created with the same mesh and position/reference 
!!   orientation information of the source mesh, and this new sibling is added to the end of the list for 
!!   the set of siblings. Siblings may have different fields (other than Position and RefOrientation). 
!!   Therefore, for a sibling, it is necessary, as with MeshCreate, to indicate the fields the sibling 
!!   will have using optional arguments. Sibling meshes should not be created unless the original mesh 
!!   has been committed first.
!!
!! If CtrlCode is MESH_UPDATECOPY, all of the allocatable fields of the destination mesh are updated 
!!   with the values of the fields in the source. (The underlying mesh is untouched.) The mesh and field 
!!   definitions of the source and destination meshes must match and both must have been already committed. 
!!   The destination mesh may be an entirely different copy or it may be a sibling of the source mesh.
   SUBROUTINE MeshCopy( SrcMesh, DestMesh, CtrlCode, ErrStat , ErrMess   &
                      ,IOS, Force, Moment, Orientation, TranslationDisp, TranslationVel &
                      ,RotationVel, TranslationAcc, RotationAcc, nScalars )
   

   
     TYPE(MeshType), TARGET,      INTENT(INOUT) :: SrcMesh           !< Mesh being copied
     TYPE(MeshType), TARGET,      INTENT(INOUT) :: DestMesh          !< Copy of mesh
     INTEGER(IntKi),              INTENT(IN)    :: CtrlCode          !< MESH_NEWCOPY, MESH_SIBLING, or
                                                                     !! MESH_UPDATECOPY
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat           !< Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess           !< Error message
    ! Optional arguments (used only if CtrlCode is MESH_SIBLING):
     INTEGER(IntKi), OPTIONAL,    INTENT(IN)    :: IOS               !< If present, IOS of new sibling: input (COMPONENT_INPUT), output(COMPONENT_OUTPUT), or state(COMPONENT_STATE)
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: Force             !< If present and true, allocate Force field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: Moment            !< If present and true, allocate Moment field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: Orientation       !< If present and true, allocate Orientation field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: TranslationDisp   !< If present and true, allocate TranslationDisp field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: TranslationVel    !< If present and true, allocate TranslationVel field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: RotationVel       !< If present and true, allocate RotationVel field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: TranslationAcc    !< If present and true, allocate TranslationAcc field
     LOGICAL,        OPTIONAL,    INTENT(IN)    :: RotationAcc       !< If present and true, allocate RotationAcc field
     INTEGER(IntKi), OPTIONAL,    INTENT(IN)    :: nScalars          !< If present and > 0 , alloc n Scalars               
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
     INTEGER i, j, k, ErrStat2


     
     
      ErrStat = ErrID_None
      ErrMess = ""

      IF (.NOT. SrcMesh%Initialized) RETURN !bjj: maybe we should first CALL MeshDestroy(DestMesh,ErrStat, ErrMess)

      IF ( CtrlCode .EQ. MESH_NEWCOPY .OR. CtrlCode .EQ. MESH_SIBLING .OR. CtrlCode .EQ. MESH_COUSIN ) THEN
         
         IF (CtrlCode .EQ. MESH_NEWCOPY) THEN
            IOS_l              = SrcMesh%IOS
            Force_l            = SrcMesh%FieldMask(MASKID_FORCE)                     
            Moment_l           = SrcMesh%FieldMask(MASKID_MOMENT)                   
            Orientation_l      = SrcMesh%FieldMask(MASKID_ORIENTATION)         
            TranslationDisp_l  = SrcMesh%FieldMask(MASKID_TRANSLATIONDISP) 
            TranslationVel_l   = SrcMesh%FieldMask(MASKID_TRANSLATIONVEL)   
            RotationVel_l      = SrcMesh%FieldMask(MASKID_ROTATIONVEL)         
            TranslationAcc_l   = SrcMesh%FieldMask(MASKID_TRANSLATIONACC)   
            RotationAcc_l      = SrcMesh%FieldMask(MASKID_ROTATIONACC)         
            nScalars_l         = SrcMesh%nScalars          
         ELSE ! Sibling or cousin
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
         END IF
            
         IF ( CtrlCode .EQ. MESH_NEWCOPY .OR. CtrlCode .EQ. MESH_COUSIN ) THEN
                                    
            CALL MeshCreate( DestMesh, IOS=IOS_l, Nnodes=SrcMesh%Nnodes, ErrStat=ErrStat, ErrMess=ErrMess &
                            ,Force=Force_l                                                                &
                            ,Moment=Moment_l                                                              &
                            ,Orientation=Orientation_l                                                    &
                            ,TranslationDisp=TranslationDisp_l                                            &
                            ,TranslationVel=TranslationVel_l                                              &
                            ,RotationVel=RotationVel_l                                                    &
                            ,TranslationAcc=TranslationAcc_l                                              &
                            ,RotationAcc=RotationAcc_l                                                    &
                            ,nScalars=nScalars_l                                                          )

            IF (ErrStat >= AbortErrLev) RETURN

            DestMesh%Position       =  SrcMesh%Position
            DestMesh%RefOrientation =  SrcMesh%RefOrientation

            DO i = 1, NELEMKINDS
               DestMesh%ElemTable(i)%nelem = SrcMesh%ElemTable(i)%nelem
               DestMesh%ElemTable(i)%maxelem = SrcMesh%ElemTable(i)%maxelem
               DestMesh%ElemTable(i)%XElement = SrcMesh%ElemTable(i)%XElement

               ALLOCATE(DestMesh%ElemTable(i)%Elements(DestMesh%ElemTable(i)%maxelem),STAT=ErrStat2)
               IF (ErrStat2 /=0) THEN
                  ErrStat = ErrID_Fatal
                  ErrMess='MeshCopy:Error allocating ElemTable%Elements.'
                  RETURN !Early return
               END IF

               DO j = 1,DestMesh%ElemTable(i)%nelem

                  ALLOCATE(DestMesh%ElemTable(i)%Elements(j)%ElemNodes(size(SrcMesh%ElemTable(i)%Elements(j)%ElemNodes)),STAT=ErrStat2)
                  IF (ErrStat2 /=0) THEN
                     ErrStat = ErrID_Fatal
                     ErrMess='MeshCopy:Error allocating ElemTable%ElemNodes.'
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

               ALLOCATE(DestMesh%ElemList(DestMesh%maxelemlist),Stat=ErrStat2)
               IF (ErrStat2 /=0) THEN
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

      ELSE IF ( CtrlCode .EQ. MESH_UPDATECOPY ) THEN
         
         IF ( SrcMesh%nNodes .NE. DestMesh%nNodes ) THEN
            ErrStat = ErrID_Fatal
            ErrMess = "MeshCopy: MESH_UPDATECOPY of meshes with different numbers of nodes."
            RETURN
         ENDIF
                  
      ELSE IF ( CtrlCode .EQ. MESH_UPDATEREFERENCE ) THEN

         IF ( SrcMesh%nNodes .NE. DestMesh%nNodes ) THEN
            ErrStat = ErrID_Fatal
            ErrMess = "MeshCopy:MESH_UPDATEREFERENCE of meshes with different numbers of nodes."
            RETURN
         ENDIF         ! if we have a different number of nodes or different element connectivity, we'll have to redo this
         
         DestMesh%Position       = SrcMesh%Position
         DestMesh%RefOrientation = SrcMesh%RefOrientation
         DestMesh%RemapFlag      = SrcMesh%RemapFlag            
                        
      ELSE
         ErrStat = ErrID_Fatal
         ErrMess  = 'MeshCopy: Invalid CtrlCode.'
         RETURN
      ENDIF

         ! These aren't shared between siblings, so they get copied, no matter what the CtrlCode:

      DestMesh%Initialized = SrcMesh%Initialized
      DestMesh%Committed   = SrcMesh%Committed
      DestMesh%refNode = SrcMesh%refNode
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

!----------------------------------------------------------------------------------------------------------------------------------
!> For a given node in a mesh, assign the coordinates of the node in the global coordinate space. 
!! If an Orient argument is included, the node will also be assigned the specified orientation 
!! (orientation is assumed to be the identity matrix if omitted). Returns a non-zero value in  
!! ErrStat if Inode is outside the range 1..Nnodes.     
   SUBROUTINE MeshPositionNode( Mesh, Inode, Pos, ErrStat, ErrMess, Orient, Ref )
   
     TYPE(MeshType),              INTENT(INOUT) :: Mesh         !< Mesh being spatio-located
     INTEGER(IntKi),              INTENT(IN   ) :: Inode        !< Number of node being located
     REAL(ReKi),                  INTENT(IN   ) :: Pos(3)       !< Xi,Yi,Zi, coordinates of node
     INTEGER(IntKi),              INTENT(  OUT) :: ErrStat      !< Error code
     CHARACTER(*),                INTENT(  OUT) :: ErrMess      !< Error message
     REAL(R8Ki), OPTIONAL,        INTENT(IN   ) :: Orient(3,3)  !< Orientation (direction cosine matrix) of node; identity by default
     LOGICAL, OPTIONAL,           INTENT(IN   ) :: Ref
     
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
     IF ( .NOT. SIZE(Mesh%Position,2) .GE. Mesh%Nnodes ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: Position array not big enough"
     ENDIF
     IF ( .NOT. ASSOCIATED(Mesh%RefOrientation) ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: RefOrientation array not associated"
     ENDIF
     IF ( .NOT. SIZE(Mesh%RefOrientation,3) .GE. Mesh%Nnodes ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshPositionNode: RefOrientation array not big enough"
     ENDIF



     IF ( ErrStat .NE. ErrID_None ) RETURN   ! early return on error

    ! Business
     Mesh%Position(:,Inode) = Pos

     IF ( PRESENT(Orient) ) THEN
        Mesh%RefOrientation(:,:,Inode) = Orient
     ELSE
        Mesh%RefOrientation(:,1,Inode) = (/ 1._R8Ki, 0._R8Ki, 0._R8Ki /)
        Mesh%RefOrientation(:,2,Inode) = (/ 0._R8Ki, 1._R8Ki, 0._R8Ki /)
        Mesh%RefOrientation(:,3,Inode) = (/ 0._R8Ki, 0._R8Ki, 1._R8Ki /)
     END IF

     IF (PRESENT(Ref)) THEN
        Mesh%RefNode = Inode
     END IF

     RETURN

   END SUBROUTINE MeshPositionNode

!----------------------------------------------------------------------------------------------------------------------------------
!> Given a mesh that has been created, spatio-located, and constructed, 
!! commit the definition of the mesh, making it ready for initialization 
!! and use. Explicitly committing a mesh provides the opportunity to precompute 
!! traversal information, neighbor lists and other information about the mesh. 
!! Returns non-zero in value of ErrStat on error.     
   SUBROUTINE MeshCommit( Mesh, ErrStat, ErrMess )
         
     TYPE(MeshType),              INTENT(INOUT) :: Mesh              !< Mesh being committed
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat           !< Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess           !< Error message
    ! Local
     REAL(ReKi)                                 :: n1_n2_vector(3)               ! vector going from node 1 to node 2 in a Line2 element
     INTEGER n0d, n1d, n2d, n3d, nn, i, j, NElem, n1,n2
     LOGICAL                                    :: NodeInElement(Mesh%NNodes)    ! determines if each node is part of an element

     !TYPE(ElemListType), POINTER :: tmp(:)

     IF (Mesh%Committed) then
       ErrStat = ErrID_Warn
       ErrMess = "MeshCommit: mesh was already committed."
       RETURN  ! Early return
     ENDIF
     
     !> Check for spatial constraints -- can't mix 1D with 2D with 3D
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

     !bjj: Not sure Mesh%ElemTable(:)%nelem can be used on all versions of gfortran
     IF ( ALL( Mesh%ElemTable(:)%nelem /= SUM(Mesh%ElemTable(:)%nelem ) ) ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshCommit: a mesh can have only one type of element."
       RETURN  ! Early return
     ENDIF

     !> Check that every node is part of an element.
      NodeInElement = .FALSE.
      DO i = 1, NELEMKINDS
         DO j = 1, Mesh%ElemTable(i)%nelem
            DO n1 = 1,NumNodes( Mesh%ElemTable(i)%XElement )
               n2 = Mesh%ElemTable(i)%Elements(j)%ElemNodes(n1)
               NodeInElement( n2 ) = .TRUE.
            END DO
         END DO
      END DO
      
      ! maybe we should check that the RefOrientation is a DCM? (or at least non-zero?)

      IF ( .NOT. ALL(NodeInElement) ) THEN
         ErrStat = ErrID_Fatal
         ErrMess = "MeshCommit: all nodes must be part of an element."
         RETURN  ! Early return
      END IF

      IF ( .NOT. ANY(Mesh%FieldMask) ) THEN
         ErrStat = ErrID_Fatal
         ErrMess = "MeshCommit: Mesh does not contain any fields."
         RETURN
      END IF
         
      
     ! make sure the arrays are allocated properly...
      IF ( SIZE(Mesh%Position,2) < Mesh%Nnodes) THEN
         ErrStat = ErrID_Fatal
         ErrMess = "MeshCommit: Position array smaller than number of nodes."
         RETURN  ! Early return
      ELSEIF ( SIZE(Mesh%Position,2) > Mesh%Nnodes ) THEN
         
         ! bjj: need to get rid of the extra storage so that this doesn't cause errors in MeshCopy....
         
         
      END IF
        
      
      
      
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
        ErrMess = "MeshCommit: Error allocating element list."
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

        n2 = Mesh%ElemTable(ELEMENT_LINE2)%Elements(j)%ElemNodes(2)
        n1 = Mesh%ElemTable(ELEMENT_LINE2)%Elements(j)%ElemNodes(1)
        n1_n2_vector = Mesh%Position(:,n2) &
                     - Mesh%Position(:,n1)

        Mesh%ElemTable(ELEMENT_LINE2)%Elements(J)%det_jac  = 0.5_ReKi * TwoNorm( n1_n2_vector )   ! = L / 2
        
        IF ( 2.0_ReKi*Mesh%ElemTable(ELEMENT_LINE2)%Elements(J)%det_jac < MIN_LINE2_ELEMENT_LENGTH ) THEN
           ErrStat = ErrID_Fatal
           ErrMess = trim(ErrMess)//"MeshCommit: Line2 element "//TRIM(Num2Lstr(j))//" has 0 length."//NewLine// &
                     "   n2 = n("//TRIM(Num2Lstr(n2))//") = ("//TRIM(Num2Lstr(Mesh%Position(1,n2)))//','//TRIM(Num2Lstr(mesh%position(2,n2)))//','//TRIM(Num2Lstr(mesh%position(3,n2))) //')'//NewLine// &
                     "   n1 = n("//TRIM(Num2Lstr(n1))//") = ("//TRIM(Num2Lstr(Mesh%Position(1,n1)))//','//TRIM(Num2Lstr(mesh%position(2,n1)))//','//TRIM(Num2Lstr(mesh%position(3,n1))) //')'//NewLine
           RETURN
        END IF

     END DO

   
     ! we're finished:

      Mesh%Committed = .TRUE.

   END SUBROUTINE MeshCommit

!----------------------------------------------------------------------------------------------------------------------------------
!> Given a mesh and an element name, construct a point element whose vertex is the 
!! node index listed as the remaining argument of the call to MeshConstructElement.
!! Returns a non-zero ErrStat value on error.     
   SUBROUTINE MeshConstructElement_1PT( Mesh, Xelement, ErrStat, ErrMess, P1 )
          
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      !< Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  !< See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   !< Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   !< Error message
     INTEGER,                     INTENT(IN   ) :: P1        !< node index for this point element

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

       CALL BumpupElementTable( Mesh, Xelement, ErrStat, ErrMess )
       IF (ErrStat >= AbortErrLev ) RETURN

       ALLOCATE(Mesh%ElemTable(ELEMENT_POINT)%Elements(Mesh%ElemTable(ELEMENT_POINT)%nelem)%ElemNodes(1),sTAT=ErrStat)
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMess = "MeshConstructElement_1PT: Couldn't allocate space for element nodes."
            RETURN
         END IF
       Mesh%ElemTable(ELEMENT_POINT)%Elements(Mesh%ElemTable(ELEMENT_POINT)%nelem)%ElemNodes(1) = P1
       Mesh%ElemTable(ELEMENT_POINT)%Elements(Mesh%ElemTable(ELEMENT_POINT)%nelem)%Xelement = ELEMENT_POINT
     ELSE
       ErrStat = ErrID_Fatal
       ErrMess = 'MeshConstructElement_1PT called for invalid element type'
     ENDIF

     RETURN

   END SUBROUTINE MeshConstructElement_1PT

!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BumpupElementTable_New( Mesh, Xelement, ErrStat, ErrMess )
   ! bjj: I am getting weird errors with some models using gfortran (ivf is fine),
   ! so I am implementing this method which does not just set a pointer, pointing to
   ! a local variable. It actually copies the data twice, but allocates the pointer 
   ! in the dataype.
   
      TYPE(MeshType),              INTENT(INOUT) :: Mesh      ! Mesh being constructed
      INTEGER(IntKi),              INTENT(IN)    :: Xelement  ! See Element Names
      INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   ! Error code
      CHARACTER(*),                INTENT(OUT)   :: ErrMess   ! Error message
   
         ! Local
      TYPE(ElemRecType),             ALLOCATABLE :: tmp(:)
      INTEGER                                    :: i 

      ErrStat = ErrID_None
      ErrMess = ""
      
   
       IF ( Mesh%ElemTable(Xelement)%nelem .GE. Mesh%ElemTable(Xelement)%maxelem ) THEN
!write(0,*)'>>>>>>>>>> bumping maxpoint',Mesh%ElemTable(Xelement)%maxelem
         
         IF (Mesh%ElemTable(Xelement)%maxelem .GT. 0 ) THEN 
            
               ! copy data in pointer to temp copy:
            ALLOCATE ( tmp(Mesh%ElemTable(Xelement)%maxelem), STAT=ErrStat )
            IF (ErrStat /= 0) THEN
               ErrStat = ErrID_Fatal
               ErrMess = "BumpupElementTable: Couldn't allocate space for element table copy."
               RETURN
            END IF
            
            DO i=1,Mesh%ElemTable(Xelement)%maxelem
               CALL Mesh_MoveAlloc_ElemRecType( Mesh%ElemTable(Xelement)%Elements(i), tmp(i) )
            END DO
            
         END IF
         
            ! deallocate the pointer, then reallocate to a larger size:
         IF ( ASSOCIATED(Mesh%ElemTable(Xelement)%Elements) ) DEALLOCATE(Mesh%ElemTable(Xelement)%Elements)
         
         ALLOCATE( Mesh%ElemTable(Xelement)%Elements( Mesh%ElemTable(Xelement)%maxelem + BUMPUP ),STAT=ErrStat )   
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMess = "BumpupElementTable: Couldn't allocate space for element table."
            IF ( ALLOCATED(tmp) ) DEALLOCATE(tmp)
            RETURN
         END IF
         
            ! copy old data to the table pointer:         
         DO i=1,Mesh%ElemTable(Xelement)%maxelem
            CALL Mesh_MoveAlloc_ElemRecType( tmp(i), Mesh%ElemTable(Xelement)%Elements(i) )
         END DO
               
         IF ( ALLOCATED(tmp) ) DEALLOCATE(tmp)
         
            ! set the new size of the element table:
         Mesh%ElemTable(Xelement)%maxelem = Mesh%ElemTable(Xelement)%maxelem + BUMPUP  
!write(0,*)'>>>>>>>>>> bumped maxpoint',Mesh%ElemTable(Xelement)%maxelem                  
      END IF
     
   END SUBROUTINE BumpupElementTable_New
   
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine increases the allocated space for Mesh%ElemTable(Xelement)%Elements
!! if adding a new element will exceed the pre-allocated space.
   SUBROUTINE BumpupElementTable( Mesh, Xelement, ErrStat, ErrMess )

   ! bjj: this is the old method of increasing the element table size.
   ! it was duplicated in MeshConstructElement_1PT and MeshConstructElement_2PT
   ! I have made it a subroutine so we don't have to duplicate it anymore.
   
      TYPE(MeshType),              INTENT(INOUT) :: Mesh      !< Mesh being constructed
      INTEGER(IntKi),              INTENT(IN)    :: Xelement  !< type of element (See Element Names)
      INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   !< Error code
      CHARACTER(*),                INTENT(OUT)   :: ErrMess   !< Error message
   
    ! Local
     TYPE(ElemRecType),             POINTER      :: tmp(:)
     INTEGER                                     :: i

      ErrStat = ErrID_None
      ErrMess = ""
                     
       IF ( Mesh%ElemTable(Xelement)%nelem .GE. Mesh%ElemTable(Xelement)%maxelem ) THEN
!write(0,*)'>>>>>>>>>> bumping maxline2',Mesh%ElemTable(Xelement)%maxelem
         ALLOCATE(tmp(Mesh%ElemTable(Xelement)%maxelem+BUMPUP),Stat=ErrStat)
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMess = "BumpupElementTableOld: Couldn't allocate space for element table"
            RETURN
         END IF
                     
         !IF (Mesh%ElemTable(Xelement)%maxelem .GT. 0 ) &
               ! tmp(1:Mesh%ElemTable(Xelement)%maxelem) = Mesh%ElemTable(Xelement)%Elements(1:Mesh%ElemTable(Xelement)%maxelem)
         DO i=1,Mesh%ElemTable(Xelement)%maxelem
            CALL Mesh_MoveAlloc_ElemRecType( Mesh%ElemTable(Xelement)%Elements(i), tmp(i) )
         END DO
                        
         IF ( ASSOCIATED(Mesh%ElemTable(Xelement)%Elements) ) DEALLOCATE(Mesh%ElemTable(Xelement)%Elements)
         Mesh%ElemTable(Xelement)%Elements => tmp
         Mesh%ElemTable(Xelement)%maxelem = Mesh%ElemTable(Xelement)%maxelem + BUMPUP
!write(0,*)'>>>>>>>>>> bumped maxline2',Mesh%ElemTable(Xelement)%maxelem
       ENDIF       
       
       
   END SUBROUTINE BumpupElementTable
      
!----------------------------------------------------------------------------------------------------------------------------------
!> Given a mesh and an element name, construct 2-point line (line2) element whose 
!! vertices are the node indices listed as the remaining arguments of the call to 
!! MeshConstructElement. The adjacency of elements is implied when elements are 
!! created that share some of the same nodes. Returns a non-zero value on error.     
   SUBROUTINE MeshConstructElement_2PT( Mesh, Xelement, ErrStat, ErrMess, P1, P2 )
      
     TYPE(MeshType),              INTENT(INOUT) :: Mesh      !< Mesh being constructed
     INTEGER(IntKi),              INTENT(IN)    :: Xelement  !< See Element Names
     INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   !< Error code
     CHARACTER(*),                INTENT(OUT)   :: ErrMess   !< Error message
     INTEGER,                     INTENT(IN   ) :: P1        !< 1 of 2 points that make a 2-point line element
     INTEGER,                     INTENT(IN   ) :: P2        !< 1 of 2 points that make a 2-point line element
     
     IF ( mesh_debug ) print*,'Called MeshConstructElement_2PT'
     ErrStat = ErrID_None
     ErrMess = ""
    ! Safety first
     IF ( .NOT. Mesh%Initialized ) THEN
       ErrStat = ErrID_Fatal
       ErrMess = "MeshConstructElement_2PT: attempt to use uncreated mesh."
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

      CALL BumpupElementTable( Mesh, Xelement, ErrStat, ErrMess )
      IF (ErrStat >= AbortErrLev) RETURN

       ALLOCATE(Mesh%ElemTable(ELEMENT_LINE2)%Elements(Mesh%ElemTable(ELEMENT_LINE2)%nelem)%ElemNodes(2),Stat=ErrStat)
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMess = "MeshConstructElement_2PT: Couldn't allocate space for element nodes."
            RETURN
         END IF
       Mesh%ElemTable(ELEMENT_LINE2)%Elements(Mesh%ElemTable(ELEMENT_LINE2)%nelem)%ElemNodes(1) = P1
       Mesh%ElemTable(ELEMENT_LINE2)%Elements(Mesh%ElemTable(ELEMENT_LINE2)%nelem)%ElemNodes(2) = P2
       Mesh%ElemTable(ELEMENT_LINE2)%Elements(Mesh%ElemTable(ELEMENT_LINE2)%nelem)%Xelement = ELEMENT_LINE2

     ELSE
       ErrMess = 'MeshConstructElement_2PT called for invalid element type'
       ErrStat = ErrID_Fatal
     ENDIF
     RETURN
   END SUBROUTINE MeshConstructElement_2PT

!----------------------------------------------------------------------------------------------------------------------------------
!> added 20130102 as stub for AeroDyn work
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

!----------------------------------------------------------------------------------------------------------------------------------
!> added 20130102 as stub for AeroDyn work
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

!----------------------------------------------------------------------------------------------------------------------------------
!> added 20130102 as stub for AeroDyn work
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

!----------------------------------------------------------------------------------------------------------------------------------
!> added 20130102 as stub for AeroDyn work
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

!----------------------------------------------------------------------------------------------------------------------------------
!> added 20130102 as stub for AeroDyn work
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

!----------------------------------------------------------------------------------------------------------------------------------
!> added 20130102 as stub for AeroDyn work
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

!----------------------------------------------------------------------------------------------------------------------------------
!> added 20130102 as stub for AeroDyn work
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

!................................................................                                                                                                                                                      
!> This routine splits a line2 element into two separate elements, using p1 as
!! the new node connecting the two new elements formed from E1.
   SUBROUTINE MeshSplitElement_2PT( Mesh, Xelement, ErrStat, ErrMess, E1, P1  )
      
      TYPE(MeshType),              INTENT(INOUT) :: Mesh      !< Mesh being constructed
      INTEGER(IntKi),              INTENT(IN)    :: Xelement  !< See Element Names
      INTEGER(IntKi),              INTENT(OUT)   :: ErrStat   !< Error code
      CHARACTER(*),                INTENT(OUT)   :: ErrMess   !< Error message
      INTEGER,                     INTENT(IN   ) :: E1        !< number of element in Element Table
      INTEGER,                     INTENT(IN   ) :: P1        !< node number

      INTEGER                                    :: p2        ! local copy of P2, in case we end up deallocating the array pointed to in the call to construct a new element
      
      IF ( mesh_debug ) print*,'Called MeshSplitElement_2PT'
      ErrStat = ErrID_None
      ErrMess = ""
      
      ! Safety first
      IF ( Xelement .NE. ELEMENT_LINE2 ) THEN
         ErrMess = 'MeshSplitElement_2PT called for invalid element type.'
         ErrStat = ErrID_Fatal        
      ELSEIF ( .NOT. Mesh%Initialized ) THEN
         ErrStat = ErrID_Fatal
         ErrMess = "MeshSplitElement_2PT: attempt to use uncreated mesh."
      ELSEIF ( P1 .LT. 1 .OR. P1 .GT. Mesh%Nnodes ) THEN
         ErrStat = ErrID_Fatal
         ErrMess ="MeshSplitElement_2PT: invalid P1 ("//TRIM(Num2LStr(P1))//") for mesh with "//TRIM(Num2LStr(Mesh%Nnodes))//" nodes."
      ELSEIF ( E1 .LT. 1 .OR. E1 .GT. Mesh%ElemTable(ELEMENT_LINE2)%nelem ) THEN
         ErrStat = ErrID_Fatal
         ErrMess ="MeshSplitElement_2PT: invalid E1 ("//TRIM(Num2LStr(E1))//") for mesh with "//TRIM(Num2LStr(Mesh%ElemTable(ELEMENT_LINE2)%nelem))//" Line2 elements."
      ELSEIF (Mesh%Committed ) THEN
         ErrStat = ErrID_Fatal
         ErrMess = " MeshSplitElement_2PT: attempt to add element to committed mesh."
      ELSEIF ( Mesh%ElemTable(ELEMENT_LINE2)%Elements(E1)%ElemNodes(1) == P1 .OR. &
              Mesh%ElemTable(ELEMENT_LINE2)%Elements(E1)%ElemNodes(2) == P1 ) THEN
         ErrStat = ErrID_Fatal
         ErrMess ="MeshSplitElement_2PT: node P1 ("//TRIM(Num2LStr(P1))//") is already a node of element E1 ("//TRIM(Num2LStr(E1))//")."                
      ENDIF
     
      IF ( ErrStat .NE. ErrID_None ) THEN
         !CALL WrScr( TRIM(ErrMess) )
         RETURN  !  early return on fatal error
      ENDIF
               
     
    ! Business
      ! E1 currently has nodes (n1,n2):
      ! Create a new element with nodes (p1,n2):
      p2 = Mesh%ElemTable(ELEMENT_LINE2)%Elements(E1)%ElemNodes(2)
      CALL MeshConstructElement( Mesh, Xelement, ErrStat, ErrMess, p1=P1, p2=p2)
    
         ! Make element E1 now have nodes (n1,p1):
      Mesh%ElemTable(ELEMENT_LINE2)%Elements(E1)%ElemNodes(2) = P1
    
      RETURN
       
   END SUBROUTINE MeshSplitElement_2PT
!................................................................                                                                           
                                                                           
!----------------------------------------------------------------------------------------------------------------------------------
!> Given a control code and a mesh that has been committed, retrieve the next element in the mesh. 
!!   Used to traverse mesh element by element. On entry, the CtrlCode argument contains a control code: 
!!   zero indicates start from the beginning, an integer between 1 and Mesh%Nelemlist returns that element,
!!   and MESH_NEXT means return the next element in traversal. On exit, CtrlCode contains the status of the 
!!   traversal in (zero or MESH_NOMOREELEMS). The routine optionally outputs the index of the element in the
!!   mesh's element list, the name of the element (see "Element Names"), and a pointer to the element.    
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
!> This subroutine returns the names of the output rows/columns in the Jacobian matrices. It assumes both force and moment
!! fields are allocated.
   SUBROUTINE PackLoadMesh_Names(M, MeshName, Names, indx_first)
   
      TYPE(MeshType)                    , INTENT(IN   ) :: M                          !< Load mesh
      CHARACTER(*)                      , INTENT(IN   ) :: MeshName                   !< name of mesh 
      CHARACTER(LinChanLen)             , INTENT(INOUT) :: Names(:)                   !< name of row/column of jacobian 
      INTEGER(IntKi)                    , INTENT(INOUT) :: indx_first                 !< index into Names array; gives location of next array position to fill
   
         ! local variables:
      INTEGER(IntKi)                :: i, j
      character(1), parameter       :: Comp(3) = (/'X','Y','Z'/)
      character(2)                  :: UnitDesc

            
      UnitDesc = ''
      if (M%Committed) then
         if (M%ElemTable(ELEMENT_LINE2)%nelem > 0) UnitDesc = '/m'
      end if
                     
      do i=1,M%NNodes
         do j=1,3
            Names(indx_first) = trim(MeshName)//' '//Comp(j)//' force, node '//trim(num2lstr(i))//', N'//UnitDesc
            indx_first = indx_first + 1
         end do      
      end do
            ! This is needed for MAP meshes because it only contains the Force field not the Moment field
      if ( M%fieldmask(MASKID_Moment) .AND. ALLOCATED(M%Moment)) then                    
         do i=1,M%NNodes
            do j=1,3
               Names(indx_first) = trim(MeshName)//' '//Comp(j)//' moment, node '//trim(num2lstr(i))//', Nm'//UnitDesc
               indx_first = indx_first + 1
            end do
         end do
      end if

   END SUBROUTINE PackLoadMesh_Names
!...............................................................................................................................
!> This subroutine returns the operating point values of the mesh fields. It assumes both force and moment
!! fields are allocated.
   SUBROUTINE PackLoadMesh(M, Ary, indx_first)
   
      TYPE(MeshType)                    , INTENT(IN   ) :: M                          !< Load mesh
      REAL(ReKi)                        , INTENT(INOUT) :: Ary(:)                     !< array to pack this mesh into 
      INTEGER(IntKi)                    , INTENT(INOUT) :: indx_first                 !< index into Ary; gives location of next array position to fill
   
         ! local variables:
      INTEGER(IntKi)                :: i, j

                     
      do i=1,M%NNodes
         do j=1,3
            Ary(indx_first) = M%Force(j,i)
            indx_first = indx_first + 1
         end do      
      end do
      
         ! This is needed for MAP meshes because it only contains the Force field not the Moment field
      if ( M%fieldmask(MASKID_Moment) .AND. ALLOCATED(M%Moment)) then     
         do i=1,M%NNodes
            do j=1,3
               Ary(indx_first) = M%Moment(j,i)
               indx_first = indx_first + 1
            end do
         end do
      end if

   END SUBROUTINE PackLoadMesh
!...............................................................................................................................
!> This subroutine computes the differences of two meshes and packs that value into appropriate locations in the dY array.
!! Do not change this packing without making sure subroutine aerodyn::init_jacobian is consistant with this routine!
   SUBROUTINE PackLoadMesh_dY(M_p, M_m, dY, indx_first)
   
      TYPE(MeshType)                    , INTENT(IN   ) :: M_p                        !< AD outputs on given mesh at \f$ u + \Delta u \f$ (p=plus)
      TYPE(MeshType)                    , INTENT(IN   ) :: M_m                        !< AD outputs on given mesh at \f$ u - \Delta u \f$ (m=minus)   
      REAL(R8Ki)                        , INTENT(INOUT) :: dY(:)                      !< column of dYdu or dYdz \f$ \frac{\partial Y}{\partial u_i} = \frac{y_p - y_m}{2 \, \Delta u}\f$ 
      INTEGER(IntKi)                    , INTENT(INOUT) :: indx_first                 !< index into dY array; gives location of next array position to fill
   
         ! local variables:
      INTEGER(IntKi)                :: i, indx_last

   
      do i=1,M_p%NNodes
         indx_last  = indx_first + 2 
         dY(indx_first:indx_last) = M_p%Force(:,i) - M_m%Force(:,i)
         indx_first = indx_last + 1
      end do

         ! This is needed for MAP meshes because it only contains the Force field not the Moment field
      if ( M_p%fieldmask(MASKID_Moment) .AND. ALLOCATED(M_p%Moment)) then     
         do i=1,M_p%NNodes
            indx_last  = indx_first + 2 
            dY(indx_first:indx_last) = M_p%Moment(:,i) - M_m%Moment(:,i)
            indx_first = indx_last + 1
         end do
      end if

   END SUBROUTINE PackLoadMesh_dY
!...............................................................................................................................
!> This subroutine returns the names of rows/columns of motion meshes in the Jacobian matrices. It assumes all fields marked
!! by FieldMask are allocated; Some fields may be allocated by the ModMesh module and not used in
!! the linearization procedure, thus I am not using the check if they are allocated to determine if they should be included.
!...............................................................................................................................
   SUBROUTINE PackMotionMesh_Names(M, MeshName, Names, indx_first, FieldMask)
   
      TYPE(MeshType)                    , INTENT(IN   ) :: M                          !< Motion mesh
      CHARACTER(*)                      , INTENT(IN   ) :: MeshName                   !< name of mesh 
      CHARACTER(LinChanLen)             , INTENT(INOUT) :: Names(:)                   !< name of row/column of jacobian
      INTEGER(IntKi)                    , INTENT(INOUT) :: indx_first                 !< index into Names array; gives location of next array position to fill
      LOGICAL, OPTIONAL                 , INTENT(IN   ) :: FieldMask(FIELDMASK_SIZE)  !< flags to determine if this field is part of the packing
      
      
         ! local variables:
      INTEGER(IntKi)                :: i, j
      character(1), parameter       :: Comp(3) = (/'X','Y','Z'/)
      LOGICAL                       :: Mask(FIELDMASK_SIZE)               !< flags to determine if this field is part of the packing

      if (present(FieldMask)) then
         Mask = FieldMask
      else
         Mask = .true.
      end if
            
   
      if (Mask(MASKID_TRANSLATIONDISP)) then
         do i=1,M%NNodes
            do j=1,3
               Names(indx_first) = trim(MeshName)//' '//Comp(j)//' translation displacement, node '//trim(num2lstr(i))//', m'
               indx_first = indx_first + 1
            end do      
         end do
      end if
      
      if (Mask(MASKID_ORIENTATION)) then
         do i=1,M%NNodes
            do j=1,3
               Names(indx_first) = trim(MeshName)//' '//Comp(j)//' orientation angle, node '//trim(num2lstr(i))//', rad'
               indx_first = indx_first + 1
            end do      
         end do
      end if
      
      if (Mask(MASKID_TRANSLATIONVEL)) then
         do i=1,M%NNodes
            do j=1,3
               Names(indx_first) = trim(MeshName)//' '//Comp(j)//' translation velocity, node '//trim(num2lstr(i))//', m/s'
               indx_first = indx_first + 1
            end do      
         end do
      end if
      
      if (Mask(MASKID_ROTATIONVEL)) then
         do i=1,M%NNodes
            do j=1,3
               Names(indx_first) = trim(MeshName)//' '//Comp(j)//' rotation velocity, node '//trim(num2lstr(i))//', rad/s'
               indx_first = indx_first + 1
            end do      
         end do
      end if
         
      if (Mask(MASKID_TRANSLATIONACC)) then
         do i=1,M%NNodes
            do j=1,3
               Names(indx_first) = trim(MeshName)//' '//Comp(j)//' translation acceleration, node '//trim(num2lstr(i))//', m/s^2'
               indx_first = indx_first + 1
            end do      
         end do
      end if
   
      if (Mask(MASKID_ROTATIONACC)) then
         do i=1,M%NNodes
            do j=1,3
               Names(indx_first) = trim(MeshName)//' '//Comp(j)//' rotation acceleration, node '//trim(num2lstr(i))//', rad/s^2'
               indx_first = indx_first + 1
            end do      
         end do
      end if


   END SUBROUTINE PackMotionMesh_Names
!...............................................................................................................................
!> This subroutine returns the operating point values of the mesh fields. It assumes all fields marked
!! by FieldMask are allocated; Some fields may be allocated by the ModMesh module and not used in
!! the linearization procedure, thus I am not using the check if they are allocated to determine if they should be included.
   SUBROUTINE PackMotionMesh(M, Ary, indx_first, FieldMask, TrimOP)
   
      TYPE(MeshType)                    , INTENT(IN   ) :: M                          !< Motion mesh
      REAL(ReKi)                        , INTENT(INOUT) :: Ary(:)                     !< array to pack this mesh into 
      INTEGER(IntKi)                    , INTENT(INOUT) :: indx_first                 !< index into Ary; gives location of next array position to fill
      LOGICAL, OPTIONAL                 , INTENT(IN   ) :: FieldMask(FIELDMASK_SIZE)  !< flags to determine if this field is part of the packing
      LOGICAL, OPTIONAL                 , INTENT(IN   ) :: TrimOP                     !< flag to determine if the orientation should be packed as a DCM or a log map
      
      
         ! local variables:
      INTEGER(IntKi)                :: i, j, k
      LOGICAL                       :: Mask(FIELDMASK_SIZE)               !< flags to determine if this field is part of the packing
      LOGICAL                       :: PackForTrimSolution
      REAL(R8Ki)                    :: logmap(3)                          !< array to pack dcm vector representation (logmaps) into 
      INTEGER(IntKi)                :: ErrStat2
      CHARACTER(ErrMsgLen)          :: ErrMsg2
      

      if (present(FieldMask)) then
         Mask = FieldMask
      else
         Mask = .true.
      end if
            
      if (present(TrimOP)) then
         PackForTrimSolution = TrimOP
      else
         PackForTrimSolution = .false.
      end if
      
   
      if (Mask(MASKID_TRANSLATIONDISP)) then
         do i=1,M%NNodes
            do j=1,3
               Ary(indx_first) = M%TranslationDisp(j,i)
               indx_first = indx_first + 1
            end do      
         end do
      end if
      
      if (Mask(MASKID_ORIENTATION)) then
         
         if (PackForTrimSolution) then
            do i=1,M%NNodes
               call DCM_logMap(M%Orientation(:,:,i), logmap, ErrStat2, ErrMsg2) !NOTE: we cannot use GetSmllRotAngs because we CANNOT assume that all DCMs in the code are small.
               do k=1,3
                  Ary(indx_first) = logmap(k)
                  indx_first = indx_first + 1
               end do
            end do
         else
            do i=1,M%NNodes
               do j=1,3
                  do k=1,3 ! note this gives us 9 values instead of 3 for this "operating point"
                     Ary(indx_first) = M%Orientation(j,k,i)
                     indx_first = indx_first + 1
                  end do
               end do
            end do
         end if
      end if
      
      if (Mask(MASKID_TRANSLATIONVEL)) then
         do i=1,M%NNodes
            do j=1,3
               Ary(indx_first) = M%TranslationVel(j,i)
               indx_first = indx_first + 1
            end do      
         end do
      end if
      
      if (Mask(MASKID_ROTATIONVEL)) then
      
         do i=1,M%NNodes
            do j=1,3
               Ary(indx_first) = M%RotationVel(j,i)
               indx_first = indx_first + 1
            end do      
         end do
      end if
         
      if (Mask(MASKID_TRANSLATIONACC)) then
         do i=1,M%NNodes
            do j=1,3
               Ary(indx_first) = M%TranslationAcc(j,i)
               indx_first = indx_first + 1
            end do
         end do
      end if
   
      if (Mask(MASKID_ROTATIONACC)) then
         if (PackForTrimSolution) then ! these are difficult to converge in a trim solution
            do i=1,M%NNodes
               do j=1,3
                  Ary(indx_first) = 0.0_ReKi
                  indx_first = indx_first + 1
               end do
            end do
         else
            do i=1,M%NNodes
               do j=1,3
                  Ary(indx_first) = M%RotationAcc(j,i)
                  indx_first = indx_first + 1
               end do
            end do
         end if

      end if


   END SUBROUTINE PackMotionMesh
!...............................................................................................................................
!> This subroutine computes the differences of two meshes and packs that value into appropriate locations in the dY array.
   SUBROUTINE PackMotionMesh_dY(M_p, M_m, dY, indx_first, FieldMask, UseSmlAngle)
   
      TYPE(MeshType)                    , INTENT(IN   ) :: M_p                        !< ED outputs on given mesh at \f$ u + \Delta u \f$ (p=plus)
      TYPE(MeshType)                    , INTENT(IN   ) :: M_m                        !< ED outputs on given mesh at \f$ u - \Delta u \f$ (m=minus)   
      REAL(R8Ki)                        , INTENT(INOUT) :: dY(:)                      !< column of dYdu \f$ \frac{\partial Y}{\partial u_i} = \frac{y_p - y_m}{2 \, \Delta u}\f$ 
      INTEGER(IntKi)                    , INTENT(INOUT) :: indx_first                 !< index into dY array; gives location of next array position to fill
      LOGICAL, OPTIONAL                 , INTENT(IN   ) :: FieldMask(FIELDMASK_SIZE)  !< flags to determine if this field is part of the packing
      LOGICAL, OPTIONAL                 , INTENT(IN   ) :: UseSmlAngle                 !< flag to determine if the orientation should be packed as a DCM or a log map
   
         ! local variables:
      INTEGER(IntKi)                :: ErrStat2 ! we're ignoring the errors about small angles
      CHARACTER(ErrMsgLen)          :: ErrMsg2  
   
      INTEGER(IntKi)                :: i, indx_last
      LOGICAL                       :: OutputSmlAngle
      !REAL(R8Ki)                    :: lambda_m(3)
      !REAL(R8Ki)                    :: lambda_p(3)
      REAL(R8Ki)                    :: angles_m(3)
      REAL(R8Ki)                    :: angles_p(3)
      REAL(R8Ki)                    :: smallAngles(3)
      REAL(R8Ki)                    :: orientation(3,3)
      LOGICAL                       :: Mask(FIELDMASK_SIZE)               !< flags to determine if this field is part of the packing

      if (present(FieldMask)) then
         Mask = FieldMask
      else
         Mask = .true.
      end if

   
      if (Mask(MASKID_TRANSLATIONDISP)) then
         do i=1,M_p%NNodes
            indx_last  = indx_first + 2 
            dY(indx_first:indx_last) = M_p%TranslationDisp(:,i) - M_m%TranslationDisp(:,i)
            indx_first = indx_last + 1
         end do
      end if
   
      if (Mask(MASKID_ORIENTATION)) then
         if (present(UseSmlAngle)) then
            OutputSmlAngle = UseSmlAngle
         else
            OutputSmlAngle = .false.
         end if

         if (OutputSmlAngle) then
            do i=1,M_p%NNodes
               !call DCM_logMap( M_m%Orientation(:,:,i), lambda_m, ErrStat2, ErrMsg2 )
               !call DCM_logMap( M_p%Orientation(:,:,i), lambda_p, ErrStat2, ErrMsg2 )
               angles_m =  GetSmllRotAngs ( M_m%Orientation(:,:,i), ErrStat2, ErrMsg2 )
               angles_p =  GetSmllRotAngs ( M_p%Orientation(:,:,i), ErrStat2, ErrMsg2 )

               !smallAngles = lambda_p - lambda_m
               smallAngles = angles_p - angles_m

               indx_last  = indx_first + 2
               dY(indx_first:indx_last) = smallAngles
               indx_first = indx_last + 1
            end do
         else
            do i=1,M_p%NNodes
               orientation = transpose(M_m%Orientation(:,:,i))
               orientation = matmul(orientation, M_p%Orientation(:,:,i))

               smallAngles = GetSmllRotAngs( orientation, ErrStat2, ErrMsg2 )

               indx_last  = indx_first + 2
               dY(indx_first:indx_last) = smallAngles
               indx_first = indx_last + 1
            end do
         endif
      end if
      
      if (Mask(MASKID_TRANSLATIONVEL)) then
         do i=1,M_p%NNodes
            indx_last  = indx_first + 2 
            dY(indx_first:indx_last) = M_p%TranslationVel(:,i) - M_m%TranslationVel(:,i)
            indx_first = indx_last + 1
         end do
      end if
      
      if (Mask(MASKID_ROTATIONVEL)) then
         do i=1,M_p%NNodes
            indx_last  = indx_first + 2 
            dY(indx_first:indx_last) = M_p%RotationVel(:,i) - M_m%RotationVel(:,i)
            indx_first = indx_last + 1
         end do
      end if
         
      if (Mask(MASKID_TRANSLATIONACC)) then
         do i=1,M_p%NNodes
            indx_last  = indx_first + 2 
            dY(indx_first:indx_last) = M_p%TranslationAcc(:,i) - M_m%TranslationAcc(:,i)
            indx_first = indx_last + 1
         end do
      end if
   
      if (Mask(MASKID_ROTATIONACC)) then
         do i=1,M_p%NNodes
            indx_last  = indx_first + 2 
            dY(indx_first:indx_last) = M_p%RotationAcc(:,i) - M_m%RotationAcc(:,i)
            indx_first = indx_last + 1
         end do
      end if


   END SUBROUTINE PackMotionMesh_dY

!...............................................................................................................................
!> This subroutine calculates a extrapolated (or interpolated) input u_out at time t_out, from previous/future time
!! values of u (which has values associated with times in t).  Order of the interpolation is 1.
    SUBROUTINE MeshExtrapInterp1(u1, u2, tin, u_out, tin_out, ErrStat, ErrMsg )
      
    TYPE(MeshType),      INTENT(IN)     :: u1                        !< Inputs at t1 > t2
    TYPE(MeshType),      INTENT(IN)     :: u2                        !< Inputs at t2
    REAL(DbKi),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
    TYPE(MeshType),      INTENT(INOUT)  :: u_out                     !< Inputs at tin_out
    REAL(DbKi),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
    INTEGER(IntKi),      INTENT(  OUT)  :: ErrStat                   !< Error status of the operation
    CHARACTER(*),        INTENT(  OUT)  :: ErrMsg                    !< Error message if ErrStat /= ErrID_None
                                                                     
      ! local variables                                              
    INTEGER(IntKi), parameter           :: order = 1                 ! order of polynomial fit (max 2)
    REAL(DbKi)                          :: t(SIZE(tin))              ! Times associated with the inputs
    REAL(DbKi)                          :: t_out                     ! Time to which to be extrap/interpd
                                                                     
    REAL(DbKi)                          :: scaleFactor               ! temporary for extrapolation/interpolation
    REAL(DbKi)                          :: tensor(3, order+1)        ! for extrapolation of orientations 
    REAL(DbKi)                          :: tensor_interp(3)          ! for extrapolation of orientations    
    REAL(DbKi)                          :: Orient(3,3)               ! for extrapolation of orientations    
    
    INTEGER(IntKi)                      :: node                      ! node counter

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
          ErrMsg = 'MeshExtrapInterp1: size(t) must equal 2.'
          RETURN
       end if

      IF ( EqualRealNos( t(1), t(2) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'MeshExtrapInterp1: t(1) must not equal t(2) to avoid a division-by-zero error.'
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
                  
         if ( EqualRealNos(t_out, t(1)) ) then
            u_out%Orientation = u1%Orientation
         elseif ( EqualRealNos(t_out, t(2)) ) then
            u_out%Orientation = u2%Orientation
         else
                                 
            DO node=1,u_out%Nnodes
            
               Orient = u1%Orientation(:,:,node)
               CALL DCM_logmap ( Orient, tensor(:,1), ErrStat, ErrMsg )
                  IF (ErrStat >= AbortErrLev ) THEN 
                     ErrMsg = 'MeshExtrapInterp1:'//TRIM(ErrMsg)
                     RETURN
                  END IF
                  
               Orient = u2%Orientation(:,:,node)
               CALL DCM_logmap ( Orient, tensor(:,2), ErrStat, ErrMsg )
                  IF (ErrStat >= AbortErrLev ) THEN 
                     ErrMsg = 'MeshExtrapInterp1:'//TRIM(ErrMsg)
                     RETURN
                  END IF
                                    
               CALL DCM_SetLogMapForInterp( tensor )            
                      
               tensor_interp  = tensor(:,1) + (tensor(:,2) - tensor(:,1)) * scaleFactor            
                                                
               u_out%Orientation(:,:,node) = DCM_exp( tensor_interp ) 
               
            END DO
            
         end if
                  
      END IF

   END SUBROUTINE MeshExtrapInterp1

!...............................................................................................................................
!> This subroutine calculates a extrapolated (or interpolated) input u_out at time t_out, from previous/future time
!! values of u (which has values associated with times in t).  Order of the interpolation is 2.
    SUBROUTINE MeshExtrapInterp2(u1, u2, u3, tin, u_out, tin_out, ErrStat, ErrMsg )
   
    TYPE(MeshType),      INTENT(IN)     :: u1                        !< Inputs at t1 > t2 > t3
    TYPE(MeshType),      INTENT(IN)     :: u2                        !< Inputs at t2 > t3
    TYPE(MeshType),      INTENT(IN)     :: u3                        !< Inputs at t3
    REAL(DbKi),          INTENT(IN   )  :: tin(:)                    !< Times associated with the inputs
    TYPE(MeshType),      INTENT(INOUT)  :: u_out                     !< Inputs at tin_out
    REAL(DbKi),          INTENT(IN   )  :: tin_out                   !< time to be extrap/interp'd to
    INTEGER(IntKi),      INTENT(  OUT)  :: ErrStat                   !< Error status of the operation
    CHARACTER(*),        INTENT(  OUT)  :: ErrMsg                    !< Error message if ErrStat /= ErrID_None
                                                                     
      ! local variables                                              
    INTEGER(IntKi), parameter           :: order = 2                 ! order of polynomial fit (max 2)

    REAL(DbKi)                          :: t(SIZE(tin))              ! Times associated with the inputs
    REAL(DbKi)                          :: t_out                     ! Time to which to be extrap/interpd                                                                     
    REAL(DbKi)                          :: scaleFactor               ! temporary for extrapolation/interpolation    
    REAL(DbKi)                          :: tensor(3, order+1)        ! for extrapolation of orientations 
    REAL(DbKi)                          :: tensor_interp(3)          ! for extrapolation of orientations 
    REAL(DbKi)                          :: Orient(3,3)               ! for extrapolation of orientations    
    
    INTEGER(IntKi)                      :: node                      ! node counter
    
    
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
         ErrMsg = 'MeshExtrapInterp2: size(t) must equal 3.'
         RETURN
      end if

      IF ( EqualRealNos( t(1), t(2) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'MeshExtrapInterp2: t(1) must not equal t(2) to avoid a division-by-zero error.'
         RETURN
      END IF
      IF ( EqualRealNos( t(2), t(3) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'MeshExtrapInterp2: t(2) must not equal t(3) to avoid a division-by-zero error.'
         RETURN
      END IF
      IF ( EqualRealNos( t(1), t(3) ) ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = 'MeshExtrapInterp2: t(1) must not equal t(3) to avoid a division-by-zero error.'
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
         u_out%RotationAcc =   u1%RotationAcc &
                             + ( t(3)**2 * ( u1%RotationAcc - u2%RotationAcc) &
                               + t(2)**2 * (-u1%RotationAcc + u3%RotationAcc) ) * scaleFactor &
                            + ( (t(2)-t(3))*u1%RotationAcc  + t(3)*u2%RotationAcc - t(2)*u3%RotationAcc )*scaleFactor*t_out
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
                     
         if ( EqualRealNos(t_out, t(1)) ) then
            u_out%Orientation = u1%Orientation
         elseif ( EqualRealNos(t_out, t(2)) ) then
            u_out%Orientation = u2%Orientation
         elseif ( EqualRealNos(t_out, t(3)) ) then
            u_out%Orientation = u3%Orientation
         else                  
            DO node=1,u_out%Nnodes
               
               Orient = u1%Orientation(:,:,node)
               CALL DCM_logmap ( Orient, tensor(:,1), ErrStat, ErrMsg )
                  IF (ErrStat >= AbortErrLev ) THEN 
                     ErrMsg = 'MeshExtrapInterp2:'//TRIM(ErrMsg)
                     RETURN
                  END IF
                  
               Orient = u2%Orientation(:,:,node)
               CALL DCM_logmap ( Orient, tensor(:,2), ErrStat, ErrMsg )
                  IF (ErrStat >= AbortErrLev ) THEN 
                     ErrMsg = 'MeshExtrapInterp2:'//TRIM(ErrMsg)
                     RETURN
                  END IF
                  
               Orient = u3%Orientation(:,:,node)
               CALL DCM_logmap ( Orient, tensor(:,3), ErrStat, ErrMsg )
                  IF (ErrStat >= AbortErrLev ) THEN 
                     ErrMsg = 'MeshExtrapInterp2:'//TRIM(ErrMsg)
                     RETURN
                  END IF
               
               CALL DCM_SetLogMapForInterp( tensor )
                                              
               tensor_interp =   tensor(:,1) &
                                 + ( t(3)**2 * (tensor(:,1) - tensor(:,2)) + t(2)**2*(-tensor(:,1) + tensor(:,3)) )*scaleFactor &
                                 + ( (t(2)-t(3))*tensor(:,1) + t(3)*tensor(:,2) - t(2)*tensor(:,3) )*scaleFactor * t_out
               u_out%Orientation(:,:,node) = DCM_exp( tensor_interp )  

            END DO
         end if
         
                                                
      END IF

   END SUBROUTINE MeshExtrapInterp2

!...............................................................................................................................
!> High level function to easily create a point mesh with one node and one element
   SUBROUTINE CreatePointMesh(mesh, posInit, orientInit, errStat, errMsg, hasMotion, hasLoads, hasAcc)
      type(MeshType),               intent(inout) :: mesh             !< Mesh to be created
      real(ReKi),                   intent(in   ) :: PosInit(3)       !< Xi,Yi,Zi, coordinates of node
      real(R8Ki),                   intent(in   ) :: orientInit(3,3)  !< Orientation (direction cosine matrix) of node; identity by default
      logical,                      intent(in   ) :: hasMotion        !< include displacements in mesh
      logical,                      intent(in   ) :: hasLoads         !< include loads in mesh
      logical, optional,            intent(in   ) :: hasAcc           !< include acceleration (default is true)      
      integer(IntKi)              , intent(out)   :: errStat          ! Status of error message
      character(*)                , intent(out)   :: errMsg           ! Error message if ErrStat /= ErrID_None
      logical                                     :: hasAcc_loc       !< include acceleration
      
      integer(IntKi)                              :: errStat2         ! local status of error message
      character(ErrMsgLen)                        :: errMsg2          ! local error message if ErrStat /= ErrID_None
      character(*), parameter                     :: RoutineName = 'CreatePointMesh'
      
      errStat = ErrID_None
      errMsg  = ''
      
      hasAcc_loc = .true.
      if (present(hasAcc)) hasAcc_loc=hasAcc

      call MeshCreate(mesh, COMPONENT_INPUT, 1, errStat2, errMsg2,  &
         Orientation=hasMotion, TranslationDisp=hasMotion, TranslationVel=hasMotion, RotationVel=hasMotion, &
         TranslationAcc=hasAcc_loc, RotationAcc=hasAcc_loc, &
         Force = hasLoads, Moment = hasLoads)
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

      call MeshPositionNode(mesh, 1, posInit, errStat2, errMsg2, orientInit); 
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

      call MeshConstructElement(mesh, ELEMENT_POINT, errStat2, errMsg2, p1=1); 
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

      call MeshCommit(mesh, errStat2, errMsg2);
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, RoutineName)

      
      ! bjj: this initialization in done in MeshCreate already...
      ! Initialize fields
      if (hasLoads) then
         mesh%Force    = 0.0_ReKi
         mesh%Moment   = 0.0_ReKi
      endif
      if (hasMotion) then
         mesh%Orientation      = mesh%RefOrientation
         mesh%TranslationDisp  = 0.0_ReKi
         mesh%TranslationVel   = 0.0_ReKi
         mesh%RotationVel      = 0.0_ReKi
      endif
      if (hasAcc_loc) then
         mesh%TranslationAcc   = 0.0_ReKi
         mesh%RotationAcc      = 0.0_ReKi
      endif

   END SUBROUTINE CreatePointMesh
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE ModMesh


