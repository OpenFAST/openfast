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
   !
   INTEGER(IntKi), parameter :: N_INPUTS = 18
   INTEGER(IntKi), parameter :: N_OUTPUTS = 6


   CHARACTER(len=4), DIMENSION(3), PARAMETER :: StrIntMethod = (/'RK4 ','AB4 ','ABM4'/)

   ! Variables for output channels
   INTEGER(IntKi), PARAMETER :: MaxOutChs   = 9 + 3*200 ! Maximum number of output channels
                                                        ! Harcoded to outputs of 200 CB modes
   INTEGER(IntKi), PARAMETER :: ID_Time     = 0
   INTEGER(IntKi), PARAMETER :: ID_PtfFx    = 1
   INTEGER(IntKi), PARAMETER :: ID_PtfFy    = 2
   INTEGER(IntKi), PARAMETER :: ID_PtfFz    = 3
   INTEGER(IntKi), PARAMETER :: ID_PtfMx    = 4
   INTEGER(IntKi), PARAMETER :: ID_PtfMy    = 5
   INTEGER(IntKi), PARAMETER :: ID_PtfMz    = 6
   INTEGER(IntKi), PARAMETER :: ID_WaveElev = 7
   INTEGER(IntKi), PARAMETER :: ID_InpFx    = 8
   INTEGER(IntKi), PARAMETER :: ID_InpFy    = 9
   INTEGER(IntKi), PARAMETER :: ID_InpFz    = 10
   INTEGER(IntKi), PARAMETER :: ID_InpMx    = 11
   INTEGER(IntKi), PARAMETER :: ID_InpMy    = 12
   INTEGER(IntKi), PARAMETER :: ID_InpMz    = 13
   INTEGER(IntKi), PARAMETER :: ID_QStart   = 14
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
     print*,'ErrMess',TRIM(ErrMess)
     write(ErrMess,'(A)') TRIM(RoutineName)//':'//TRIM(ErrMess)
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
subroutine disp1r8(u,varname,a)
    integer,intent(in) ::u
    character(len=*),intent(in)::varname
    real(ReKi),intent(in),dimension(:) ::a
    integer :: n
    character(len=20) :: fmt
    character(len=*),parameter :: RFMT='EN13.3E2'
    n=size(a,1)
    if (n>0) then
        write(fmt,*) n
        write(u,"(A,"//adjustl(fmt)//RFMT//",A)") varname//" =[ ", a(:), " ];"
    else
        write(u,'(A,A)') varname,'=[];'
    endif
end subroutine
subroutine disp1i(u,varname,a)
    integer,intent(in) ::u
    character(len=*),intent(in)::varname
    integer(IntKi),intent(in),dimension(:) ::a
    integer :: n
    character(len=20) :: fmt
    character(len=*),parameter :: RFMT='I5'
    n=size(a,1)
    if (n>0) then
        write(fmt,*) n
        write(u,"(A,"//adjustl(fmt)//RFMT//",A)") varname//" =[ ", a(:), " ];"
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
   INTEGER                      :: I                                               ! Generic loop-counting index
   INTEGER                      :: INDX                                            ! Index for valid arrays
   CHARACTER(ChanLen)           :: OutListTmp                                      ! A string to temporarily hold OutList(I)
   CHARACTER(*), PARAMETER      :: RoutineName = "SetOutParam"

   CHARACTER(OutStrLenM1), PARAMETER  :: ValidParamAry(13) =  (/ & ! This lists the names of the allowed parameters, which must be sorted alphabetically
                               "INPF_FX  ","INPF_FY  ","INPF_FZ  ","INPF_MX  ","INPF_MY  ","INPF_MZ  ",&
                               "INTRFFX  ","INTRFFY  ","INTRFFZ  ","INTRFMX  ","INTRFMY  ","INTRFMZ  ",&
                               "WAVELEV  "/) 
   CHARACTER(OutStrLenM1), PARAMETER :: ParamUnitsAry(13) =  (/ &                     ! This lists the units corresponding to the allowed parameters
                               "(N)      ","(N)      ","(N)      ","(Nm)     ","(Nm)     ","(Nm)     ",&
                               "(N)      ","(N)      ","(N)      ","(Nm)     ","(Nm)     ","(Nm)     ","(m)      "/)
   INTEGER(IntKi), PARAMETER :: ParamIndxAry(13) =  (/ &                            ! This lists the index into AllOuts(:) of the allowed parameters ValidParamAry(:)
                              ID_InpFx, ID_InpFy, ID_InpFz, ID_InpMx, ID_InpMy, ID_InpMz,&
                              ID_PtfFx, ID_PtfFy, ID_PtfFz, ID_PtfMx, ID_PtfMy, ID_PtfMz,&
                              ID_WaveElev  /)
   character(ErrMsgLen)                         :: WarnMsg  !Warning Message
   ErrStat = ErrID_None
   ErrMsg  = ""
   WarnMsg  = ""
   
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
          call setDOFChannel(5,ID_QStart+0*p%nCBFull-1); if(Failed()) return ! NOTE: using full CB
      else if (index(OutListTmp,'CBQD_') > 0 ) then
          call setDOFChannel(6,ID_QStart+1*p%nCBFull-1); if(Failed()) return ! NOTE: using full CB
      else if (index(OutListTmp,'CBF_') > 0 ) then
          call setDOFChannel(5,ID_QStart+2*p%nCBFull-1); if(Failed()) return ! NOTE: using full CB
      else
          call setInvalidChannel() ! INVALID
      endif
      !write(*,*) p%OutParam(I)%Name, p%OutParam(I)%Indx, p%OutParam(I)%Units
   end do
   if (len(trim(WarnMsg))>0) then
       call SetErrStat(ErrID_Warn, WarnMsg,ErrStat,ErrMsg,'ExtPtfm_SetOutParam')
       write(*,'(A)')trim(WarnMsg)
   endif
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
        if ( any( p%ActiveCBDOF== iDOF ) ) then
            p%OutParam(I)%Indx  = nOffset+iDOF
            p%OutParam(I)%Units = '(-)'
        else
!         if ((iDOF> p%nCB) .or. (iDOF<1)) then
            call setInvalidChannel() ! INVALID
!         else
        endif
    end subroutine
    subroutine setInvalidChannel()
        ! If a selected output channel is not available by this module set ErrStat = ErrID_Warn.
        p%OutParam(I)%Units = "INVALID"
        p%OutParam(I)%Indx = 0
        WarnMsg=trim(WarnMsg)//TRIM(p%OutParam(I)%Name)//" is not an available output channel."//NewLine
!         call SetErrStat(ErrID_Warn, TRIM(p%OutParam(I)%Name)//" is not an available output channel.",ErrStat,ErrMsg,'ExtPtfm_SetOutParam')
!         write(*,*)TRIM(p%OutParam(I)%Name)//" is not an available output channel."
    end subroutine
END SUBROUTINE SetOutParam
!----------------------------------------------------------------------------------------------------------------------------------
!> Checks that all inputs were correctly read
subroutine CheckReducedInputs(Inp, p, ErrStat, ErrMsg)
    TYPE(ExtPtfm_InputFile),     INTENT(INOUT) :: Inp        !< Data stored in the module's input file
    TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p          !< All the parameter matrices stored in this input file
    INTEGER(IntKi),              INTENT(OUT)   :: ErrStat    !< Error status
    CHARACTER(*),                INTENT(OUT)   :: ErrMsg     !< Error message
    INTEGER(IntKi)                             :: i,j
    REAL(ReKi),                  PARAMETER     :: RelTol = 0.001_ReKi
    REAL(ReKi),                  PARAMETER     :: AbsTol = 0.0001_ReKi
    ErrStat = ErrID_None
    ErrMsg  = ""

    if (p%nTot<0)                 then ; ErrStat=ErrID_Fatal; ErrMsg='The total number of DOF was not set'; return; endif
    if (.not.allocated(p%Mass))   then ; ErrStat=ErrID_Fatal; ErrMsg='The mass matrix was not allocated.' ; return; endif
    if (.not.allocated(p%Stff))   then ; ErrStat=ErrID_Fatal; ErrMsg='The stiffness matrix was not allocated.' ; return; endif
    if (.not.allocated(p%Damp))   then ; ErrStat=ErrID_Fatal; ErrMsg='The damping matrix was not allocated.' ; return; endif
    if (.not.allocated(p%W0))     then ; ErrStat=ErrID_Fatal; ErrMsg='Constant self-weight was not allocated.';return; endif
    if (.not.allocated(p%WStff))  then ; ErrStat=ErrID_Fatal; ErrMsg='Self-weight stiffness matrix was not allocated.';return; endif
    !if (.not.allocated(p%Forces)) then ; ErrStat=ErrID_Fatal; ErrMsg='The loads were not allocated.';return; endif
    !if (.not.allocated(p%times))  then ; ErrStat=ErrID_Fatal; ErrMsg='The time vector was not allocated.'; return; endif
    if (allocated(Inp%ActiveCBDOF)) then 
        if (maxval(Inp%ActiveCBDOF)>size(p%Mass,1)-6) then
            ErrStat=ErrID_Fatal; ErrMsg='The maximum index of `ActiveCBDOF` (active CB DOF) should be less than the total number of CB DOF.'; return;
        endif
    endif

    ! Check mass matrix
    do i = 1,p%nTot
        do j = i+1,p%nTot
            if ( ABS( p%Mass(i,j) - p%Mass(j,i) ) > MAX(RelTol*ABS(p%Mass(i,j)),AbsTol) ) then
                   ErrStat=ErrID_Fatal; ErrMsg='Mass matrix is not symmetric'; return;
            end if
        end do
    end do
    if (Inp%HasRBMode) then
        ! Check stiffness matrix
        do i = 1,p%nTot
            do j = 1,6
                if ( ABS( p%Stff(i,j) ) > AbsTol ) then
                     ErrStat=ErrID_Fatal; ErrMsg='When rigid-body modes are enabled, the stiffness matrix should have only zeros in the first 6 rows and columns.'; return;
                end if
            end do
        end do
        do i = 1,6
            do j = 7,p%nTot
                if ( ABS( p%Stff(i,j) ) > AbsTol ) then
                     ErrStat=ErrID_Fatal; ErrMsg='When rigid-body modes are enabled, the stiffness matrix should have only zeros in the first 6 rows and columns.'; return;
                end if
            end do
        end do
        ! Check mass matrix
        p%RBMass   =  p%Mass(1,1)
        p%RBCoG(1) =  p%Mass(2,6)/p%RBMass
        p%RBCoG(2) = -p%Mass(1,6)/p%RBMass
        p%RBCoG(3) =  p%Mass(1,5)/p%RBMass
        p%RBInertia = p%Mass(4:6,4:6)
        do i = 1,3
            do j = i,3
                if ( i == j ) then
                    if ( ABS( p%Mass(i,j) - p%RBMass ) > MAX( RelTol*ABS(p%RBMass), AbsTol) ) then
                        ErrStat=ErrID_Fatal; ErrMsg='Mass matrix associated with the first 6 modes inconsistent with rigid-body modes.'; return;
                    end if
                else
                    if ( ABS( p%Mass(i,j) ) > AbsTol ) then
                        ErrStat=ErrID_Fatal; ErrMsg='Mass matrix associated with the first 6 modes inconsistent with rigid-body modes.'; return;
                    end if
                end if
            end do
            if ( ABS( p%Mass(i,i+3) ) > AbsTol ) then
                ErrStat=ErrID_Fatal; ErrMsg='Mass matrix associated with the first 6 modes inconsistent with rigid-body modes.'; return;
            end if
        end do
        if ( ABS( p%Mass(1,5) + p%Mass(2,4) ) > MAX( RelTol*ABS(p%Mass(1,5)), AbsTol ) ) then
            ErrStat=ErrID_Fatal; ErrMsg='Mass matrix associated with the first 6 modes inconsistent with rigid-body modes.'; return;
        end if
        if ( ABS( p%Mass(1,6) + p%Mass(3,4) ) > MAX( RelTol*ABS(p%Mass(1,6)), AbsTol ) ) then
            ErrStat=ErrID_Fatal; ErrMsg='Mass matrix associated with the first 6 modes inconsistent with rigid-body modes.'; return;
        end if
        if ( ABS( p%Mass(2,6) + p%Mass(3,5) ) > MAX( RelTol*ABS(p%Mass(2,6)), AbsTol ) ) then
            ErrStat=ErrID_Fatal; ErrMsg='Mass matrix associated with the first 6 modes inconsistent with rigid-body modes.'; return;
        end if
        ! Check W0
        if (     ABS(p%W0(1)) > AbsTol &
            .or. ABS(p%W0(2)) > AbsTol &
            .or. ABS(p%W0(6)) > AbsTol ) then
            ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, constant self-weight can only have non-zero entires in the (3) heave, (4) roll, and (5) pitch directions.'; return;
        end if
        if ( ABS( p%W0(4) - p%W0(3)*p%RBCoG(2) ) > MAX(RelTol*ABS(p%W0(4)), AbsTol) .or. &
             ABS( p%W0(5) + p%W0(3)*p%RBCoG(1) ) > MAX(RelTol*ABS(p%W0(4)), AbsTol) ) then
            ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, rigid-body center of mass inconsistent between constant self-weight and mass matrix.'; return;
        end if
        ! Check WStff
        do i = 1,p%nTot
           if (     ABS(p%WStff(i,1)) > AbsTol &
               .or. ABS(p%WStff(i,2)) > AbsTol &
               .or. ABS(p%WStff(i,3)) > AbsTol &
               .or. ABS(p%WStff(i,6)) > AbsTol ) then
               ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, self-weight stiffness matrix must contain only zeros in the first (surge), second (sway), third (heave), and sixth (yaw) columns.'; return;
            end if
        end do
        do i = 1,6
            do j = 4,5
               if (i==1 .and. j==5) then
                  if ( ABS( p%WStff(i,j) - p%W0(3)) > MAX( RelTol * ABS(p%W0(3)), AbsTol) ) then
                     ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, self-weight stiffness matrix entry (1,5) must be equal to  W0(3).'; return;
                  end if
               else if (i==2 .and. j==4) then
                  if ( ABS( p%WStff(i,j) + p%W0(3)) > MAX( RelTol * ABS(p%W0(3)), AbsTol) ) then
                     ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, self-weight stiffness matrix entry (2,4) must be equal to -W0(3).'; return;
                  end if
               else if (i==4 .and. j==4) then 
                  if ( ABS( p%WStff(i,j) - p%W0(3)*p%RBCoG(3) ) > MAX( RelTol * ABS(p%W0(3)*p%RBCoG(3)), AbsTol) ) then
                     ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, self-weight stiffness matrix entry (4,4) must be equal to  W0(3)*zCG.'; return;
                  end if
               else if (i==5 .and. j==5) then
                  if ( ABS( p%WStff(i,j) - p%W0(3)*p%RBCoG(3) ) > MAX( RelTol * ABS(p%W0(3)*p%RBCoG(3)), AbsTol) ) then
                     ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, self-weight stiffness matrix entry (5,5) must be equal to  W0(3)*zCG.'; return;
                  end if
               else if (i==6 .and. j==4) then 
                  if ( ABS( p%WStff(i,j) + p%W0(3)*p%RBCoG(1) ) > MAX( RelTol * ABS(p%W0(3)*p%RBCoG(1)), AbsTol) ) then
                     ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, self-weight stiffness matrix entry (6,4) must be equal to -W0(3)*xCG.'; return;
                  end if
               else if (i==6 .and. j==5) then
                  if ( ABS( p%WStff(i,j) + p%W0(3)*p%RBCoG(2) ) > MAX( RelTol * ABS(p%W0(3)*p%RBCoG(2)), AbsTol) ) then
                     ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, self-weight stiffness matrix entry (6,5) must be equal to -W0(3)*yCG.'; return;
                  end if
               else
                  if ( ABS( p%WStff(i,j) ) > AbsTol ) then
                     ErrStat=ErrID_Fatal; ErrMsg='With RBMod>0, self-weight stiffness matrix entry ('//TRIM(num2lstr(i))//','//TRIM(num2lstr(j))//') should be zero.'; return;
                  end if
               end if
            end do
        end do
    end if
end subroutine CheckReducedInputs
subroutine CheckConnInputs(Inp, InitInp, p, ErrStat, ErrMsg)
    TYPE(ExtPtfm_InputFile),     INTENT(INOUT) :: Inp        !< Data stored in the module's input file
    TYPE(ExtPtfm_InitInputType), INTENT(IN   ) :: InitInp     !< Input data for initialization routine
    TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p          !< All the parameter matrices stored in this input file
    INTEGER(IntKi),              INTENT(OUT)   :: ErrStat    !< Error status
    CHARACTER(*),                INTENT(OUT)   :: ErrMsg     !< Error message
    INTEGER(IntKi)                             :: i,j,iConn
    real(ReKi), dimension(3)                   :: RelPosConn0
    REAL(ReKi),                  PARAMETER     :: RelTol = 0.001_ReKi
    REAL(ReKi),                  PARAMETER     :: AbsTol = 0.0001_ReKi
    ErrStat = ErrID_None
    ErrMsg  = ""

    if ( Inp%HasRBMode ) then
        ! If rigid-body modes are present, the first six columns of PhiConn must be consistent with rigid-body kinematics
        do iConn=1,p%nConn
            RelPosConn0 = p%PosConn(iConn,:) - [InitInp%PtfmRefxt,InitInp%PtfmRefyt,InitInp%PtfmRefzt]
            do i = 1,3
               do j = 1,3
                  if ( i == j ) then
                     if ( ABS( p%PhiConn(3*(iConn-1)+i,j) - 1.0_ReKi ) > AbsTol ) then
                        ErrStat=ErrID_Fatal; ErrMsg='With rigid-body modes, entry ('//TRIM(num2lstr(3*(iConn-1)+i))//','//TRIM(num2lstr(j))//') of the connection displacement matrix should be 1.'; return;
                     end if
                  else
                     if ( ABS( p%PhiConn(3*(iConn-1)+i,j) ) > AbsTol ) then
                        ErrStat=ErrID_Fatal; ErrMsg='With rigid-body modes, entry ('//TRIM(num2lstr(3*(iConn-1)+i))//','//TRIM(num2lstr(j))//') of the connection displacement matrix should be 0.'; return;
                     end if
                  end if
               end do
               if ( ABS( p%PhiConn(3*(iConn-1)+i,i+3) ) > AbsTol ) then
                  ErrStat=ErrID_Fatal; ErrMsg='With rigid-body modes, entry ('//TRIM(num2lstr(3*(iConn-1)+i))//','//TRIM(num2lstr(i+3))//') of the connection displacement matrix should be 0.'; return;
               end if
            end do
            if ( ABS( p%PhiConn(3*(iConn-1)+1,5) - RelPosConn0(3) ) > MAX( RelTol*ABS(RelPosConn0(3)) , AbsTol) ) then
               ErrStat=ErrID_Fatal; ErrMsg='With rigid-body modes, entry ('//TRIM(num2lstr(3*(iConn-1)+1))//','//TRIM(num2lstr(5))//') of the connection displacement matrix should be '//TRIM(num2lstr( RelPosConn0(3)))//'.'; return;
            end if
            if ( ABS( p%PhiConn(3*(iConn-1)+1,6) + RelPosConn0(2) ) > MAX( RelTol*ABS(RelPosConn0(2)) , AbsTol) ) then
               ErrStat=ErrID_Fatal; ErrMsg='With rigid-body modes, entry ('//TRIM(num2lstr(3*(iConn-1)+1))//','//TRIM(num2lstr(6))//') of the connection displacement matrix should be '//TRIM(num2lstr(-RelPosConn0(2)))//'.'; return;
            end if
            if ( ABS( p%PhiConn(3*(iConn-1)+2,4) + RelPosConn0(3) ) > MAX( RelTol*ABS(RelPosConn0(3)) , AbsTol) ) then
               ErrStat=ErrID_Fatal; ErrMsg='With rigid-body modes, entry ('//TRIM(num2lstr(3*(iConn-1)+2))//','//TRIM(num2lstr(4))//') of the connection displacement matrix should be '//TRIM(num2lstr(-RelPosConn0(3)))//'.'; return;
            end if
            if ( ABS( p%PhiConn(3*(iConn-1)+2,6) - RelPosConn0(1) ) > MAX( RelTol*ABS(RelPosConn0(1)) , AbsTol) ) then
               ErrStat=ErrID_Fatal; ErrMsg='With rigid-body modes, entry ('//TRIM(num2lstr(3*(iConn-1)+2))//','//TRIM(num2lstr(6))//') of the connection displacement matrix should be '//TRIM(num2lstr( RelPosConn0(1)))//'.'; return;
            end if
            if ( ABS( p%PhiConn(3*(iConn-1)+3,4) - RelPosConn0(2) ) > MAX( RelTol*ABS(RelPosConn0(2)) , AbsTol) ) then
               ErrStat=ErrID_Fatal; ErrMsg='With rigid-body modes, entry ('//TRIM(num2lstr(3*(iConn-1)+3))//','//TRIM(num2lstr(4))//') of the connection displacement matrix should be '//TRIM(num2lstr( RelPosConn0(2)))//'.'; return;
            end if
            if ( ABS( p%PhiConn(3*(iConn-1)+3,5) + RelPosConn0(1) ) > MAX( RelTol*ABS(RelPosConn0(1)) , AbsTol) ) then
               ErrStat=ErrID_Fatal; ErrMsg='With rigid-body modes, entry ('//TRIM(num2lstr(3*(iConn-1)+3))//','//TRIM(num2lstr(5))//') of the connection displacement matrix should be '//TRIM(num2lstr(-RelPosConn0(1)))//'.'; return;
            end if
        end do
    end if

end subroutine CheckConnInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadPrimaryFile(InputFile, InitInp, p, OutFileRoot, InputFileData, ErrStat, ErrMsg)
!..................................................................................................................................
   ! Passed variables
   CHARACTER(*),                INTENT(IN)    :: InputFile      !< Name of the file containing the primary input data
   TYPE(ExtPtfm_InitInputType), INTENT(IN)    :: InitInp        !< Initialization input
   TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p              !< All the parameter matrices stored in this input file
   CHARACTER(*),                INTENT(IN)    :: OutFileRoot    !< The rootname of all the output files written by this routine.
   TYPE(ExtPtfm_InputFile),     INTENT(OUT)   :: InputFileData ! Data stored in the module's input file
   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat        !< Error status                              
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg         !< Error message
   ! Local variables:
   INTEGER(IntKi)                       :: I                    ! loop counter
   INTEGER(IntKi)                       :: N                    ! Number of list elements
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
   ! Rigid-body mode flag
   CALL ReadVar(UnIn, InputFile, InputFileData%RBMod, "RBMod", "Method for handling rigid-body motion", ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   InputFileData%hasRBMode = InputFileData%RBMod > 0_IntKi
   ! Reduction Filename
   CALL ReadVar(UnIn, InputFile, InputFileData%RedFile   , 'Red_FileName', 'Path containing Guyan/Craig-Bampton inputs', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   IF ( PathIsRelative(InputFileData%RedFile) ) InputFileData%RedFile = TRIM(PriPath)//TRIM(InputFileData%RedFile)
   CALL ReadVar(UnIn, InputFile, InputFileData%RedFileCst, 'RedCst_FileName', 'Path containing Guyan/Craig-Bampton constant inputs', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   IF ( PathIsRelative(InputFileData%RedFileCst) ) InputFileData%RedFileCst = TRIM(PriPath)//TRIM(InputFileData%RedFileCst)
   CALL ReadVar(UnIn, InputFile, N , 'NActiveCBDOF','Number of active CB mode listed in ActiveCBDOF, -1 for all modes', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return
   if (N<0) then
       CALL ReadCom(UnIn, InputFile, 'ActiveCBDOF', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   elseif (N==0) then
       ! Allocating ActiveDOF of size 0 => Guyan modes only
       CALL AllocAry(InputFileData%ActiveCBDOF, N, 'ActiveCBDOF',  ErrStat, ErrMsg ); if (Failed()) return
       CALL ReadCom(UnIn, InputFile, 'ActiveCBDOF', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   else
       CALL AllocAry(InputFileData%ActiveCBDOF, N, 'ActiveCBDOF',  ErrStat, ErrMsg ); if (Failed()) return
       CALL ReadAry(UnIn, InputFile, InputFileData%ActiveCBDOF, N, 'ActiveCBDOF', 'List of active CB modes', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   endif
   ! TODO TODO TODO CALL ReadVar(UnIn, InputFile, InputFileData%EquilStart, 'EquilStart','Find the equilibrium initial positions for the CB modes', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return
   CALL ReadVar(UnIn, InputFile, N , 'NInitPosList','Number of initial positions listed in InitPosList', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return
   if (N<=0) then
       CALL ReadCom(UnIn, InputFile, 'InitPosList', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   else
       CALL AllocAry(InputFileData%InitPosList, N, 'InitPosList',  ErrStat, ErrMsg ); if (Failed()) return
       CALL ReadAry(UnIn, InputFile, InputFileData%InitPosList, N, 'InitPosList', 'Initial positions', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   endif
   CALL ReadVar(UnIn, InputFile, N , 'NInitVelList','Number of initial velocties listed in InitVelList', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return
   if (N<=0) then
       CALL ReadCom(UnIn, InputFile, 'InitVelList', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   else
       CALL AllocAry(InputFileData%InitVelList, N, 'InitVelList',  ErrStat, ErrMsg ); if (Failed()) return
       CALL ReadAry(UnIn, InputFile, InputFileData%InitVelList, N, 'InitVelList', 'Initial velocities', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   endif
   !---------------------- CONNECTION INPUTS ---------------------------------------
   CALL ReadCom(UnIn, InputFile, 'Section Header: Connections', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   CALL ReadVar(UnIn, InputFile, InputFileData%HasConnections, 'Connections','Flag for connections', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return
   CALL ReadVar(UnIn, InputFile, InputFileData%ConnFile, 'Conn_FileName', 'Path containing connections inputs', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   !---------------------- USER FORCING INPUTS ---------------------------------------
   CALL ReadCom(UnIn, InputFile, 'Section Header: User Forcing', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
   CALL ReadVar(UnIn, InputFile, InputFileData%HasUserForcing, 'UserForcing','Flag for user prescribed forcing', ErrStat, ErrMsg, UnEc ); if(LineFailed()) return
   CALL ReadVar(UnIn, InputFile, InputFileData%ForceFile, 'Force_FileName', 'Path containing user forcing inputs', ErrStat, ErrMsg, UnEc); if(LineFailed()) return
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
   call ReadReducedFile(InputFileData%RedFile, p, ErrStat, ErrMsg); if(Failed()) return
   ! Checking that everyting was correctly read and set
   call CheckReducedInputs(InputFileData, p, ErrStat, ErrMsg);  if(Failed()) return

   ! --- Reading connection file
   if (InputFileData%hasConnections) then
      call ReadConnFile(InputFileData%ConnFile, p, ErrStat, ErrMsg); if(Failed()) return
      call CheckConnInputs(InputFileData, InitInp, p, ErrStat, ErrMsg); if(Failed()) return
   else
      p%NConn = 0_IntKi
   end if

   ! --- Reading user forcing file
   if (InputFileData%HasUserForcing) then
      call ReadForceFile(InputFileData%ForceFile, p, ErrStat, ErrMsg); if(Failed()) return
   else
      p%nTimeSteps = 1_IntKi
      call allocAry( p%Forces, 1_IntKi, p%nTot, 'p%Forces', ErrStat, ErrMsg); if(Failed()) return
      call allocAry( p%times , 1_IntKi,         'p%times' , ErrStat, ErrMsg); if(Failed()) return
      p%Forces= 0.0_ReKi
      p%times = 0.0_ReKi
   end if
 
   ! --- Reducing the number of DOF if needed
   p%nCBFull=p%nCB
   if (allocated(InputFileData%ActiveCBDOF)) then
       call allocAry(p%ActiveCBDOF, size(InputFileData%ActiveCBDOF), 'ActiveCBDOF',  ErrStat, ErrMsg); if(Failed()) return
       do I=1,size(InputFileData%ActiveCBDOF)
           p%ActiveCBDOF(I) = InputFileData%ActiveCBDOF(I);
       enddo
       call ReduceNumberOfDOF(p, ErrStat, ErrMsg);
   else
       call allocAry(p%ActiveCBDOF, p%nCBFull, 'ActiveCBDOF',  ErrStat, ErrMsg); if(Failed()) return
       do I=1,p%nCBFull
           p%ActiveCBDOF(I) = I
       enddo
   endif

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

!> Reduce the number of degrees of freedom given as input
SUBROUTINE ReduceNumberOfDOF(p, ErrStat, ErrMsg)
   TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p                                   !< All the parameter matrices stored in this input file
   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                             !< Error status                              
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                              !< Error message
   integer(IntKi) :: nActive
   integer(IntKi), dimension(:), allocatable :: FullActiveCBDOF
   integer(IntKi) :: I

   ! Preprending 1-6 to ActiveDOF
   call allocAry(FullActiveCBDOF, size(p%ActiveCBDOF)+6, 'FullActiveCBDOF',  ErrStat, ErrMsg); if(Failed()) return
   FullActiveCBDOF(1:6)=(/1,2,3,4,5,6/)
   do I=1,size(p%ActiveCBDOF);
       FullActiveCBDOF(I+6)=p%ActiveCBDOF(I)+6;
   enddo
   nActive=size(FullActiveCBDOF)

   ! Reducing matrices and load matrix
   call SquareMatRed(p%Mass)
   call SquareMatRed(p%Stff)
   call SquareMatRed(p%Damp)
   call TimeMatRed(p%Forces)
   if (allocated(p%PhiConn)) then
      call RectMatRed(p%PhiConn)
   end if

   ! Trigger
   p%nCB = size(p%ActiveCBDOF)
   p%nTot= p%nCB+6
CONTAINS
    !> Takes M and returns M(I,I) where I is a list of indexes to keep
    subroutine SquareMatRed(M)
        real(Reki), dimension(:,:), allocatable :: M
        real(Reki), dimension(:,:), allocatable :: tmp
        integer(IntKi) :: I,J
        ! Storing M to a tmp array
        call allocAry( tmp, size(M,1), size(M,2), 'Mtmp',  ErrStat, ErrMsg); if(Failed()) return
        tmp=M
        ! Reallocating M and storing only the desired DOF
        deallocate(M)
        call allocAry(M, nActive, nActive, 'M',  ErrStat, ErrMsg); if(Failed()) return
        do I=1,nActive
            do J=1,nActive
                M(I,J) = tmp(FullActiveCBDOF(I), FullActiveCBDOF(J))
            enddo
        enddo
        deallocate(tmp)
    end subroutine 
    !> Takes M and returns M(:,I) where I is a list of indexes to keep
    subroutine RectMatRed(M)
        real(Reki), dimension(:,:), allocatable :: M
        real(Reki), dimension(:,:), allocatable :: tmp
        integer(IntKi) :: J
        ! Storing M to a tmp array
        call allocAry( tmp, size(M,1), size(M,2), 'Mtmp',  ErrStat, ErrMsg); if(Failed()) return
        tmp=M
        ! Reallocating M and storing only the desired DOF
        deallocate(M)
        call allocAry(M, size(tmp,1), nActive, 'M',  ErrStat, ErrMsg); if(Failed()) return
        do J=1,nActive
            M(:,J) = tmp(:, FullActiveCBDOF(J))
        enddo
        deallocate(tmp)
    end subroutine
    !> Takes M and returns M(:,I) where I is a list of indexes to keep
    subroutine TimeMatRed(M)
        real(Reki), dimension(:,:), allocatable :: M
        real(Reki), dimension(:,:), allocatable :: tmp
        integer(IntKi) :: I,J
        ! Storing M to a tmp array
        call allocAry( tmp, size(M,1), size(M,2), 'MTimeTmp',  ErrStat, ErrMsg); if(Failed()) return
        tmp=M
        ! Reallocating M and storing only the desired DOF
        deallocate(M)
        call allocAry(M, size(tmp,1), nActive, 'MTime',  ErrStat, ErrMsg); if(Failed()) return
        do I=1,size(tmp,1)
            do J=1,nActive
                M(I,J) = tmp(I, FullActiveCBDOF(J))
            enddo
        enddo
        deallocate(tmp)
    end subroutine 
    logical function Failed()
        CALL SetErrStatSimple(ErrStat, ErrMsg, 'ExtPtfm_ReduceNumberOfDOF')
        Failed =  ErrStat >= AbortErrLev
    end function Failed
END SUBROUTINE ReduceNumberOfDOF


!..................................................................................................................................
SUBROUTINE ReadReducedFile( InputFile, p, ErrStat, ErrMsg )
!..................................................................................................................................
   ! Passed variables
   CHARACTER(*),                INTENT(IN)    :: InputFile                           !< Name of the file containing the primary input data
   TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p                                   !< All the parameter matrices stored in this input file
   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                             !< Error status
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                              !< Error message
   ! Local variables:
   REAL(ReKi), dimension(:),allocatable :: TmpAry                                 ! temporary array for reading row from file
   INTEGER(IntKi)                       :: I                                         ! loop counter
   INTEGER(IntKi)                       :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)                       :: iLine                                     ! Current position in file
   CHARACTER(4096)                      :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
   ErrStat = ErrID_None
   ErrMsg  = ""

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
        ! Line content is analyzed as case insensitive
        call Conv2UC(Line)
        if (index(Line,'!DIMENSION')==1) then
            p%nTot = ReadIntFromStr(Line(12:), '`dimension`, file '//trim(InputFile)//', line '//Num2LStr(iLine), ErrStat, ErrMsg); if (ErrStat /= 0) exit
            p%nCB=p%nTot-6

        else if (index(Line,'!MASS')==1) then
            if (p%nTot<0) exit
            call ReadRealMatrix(UnIn, InputFile, p%Mass, 'Mass Matrix', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)
        else if (index(Line,'!STIFFNESS')==1) then
            if (p%nTot<0) exit
            call ReadRealMatrix(UnIn, InputFile, p%Stff, 'Stiffness Matrix', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)
        else if (index(Line,'!DAMPING')==1) then
            if (p%nTot<0) exit
            call ReadRealMatrix(UnIn, InputFile, p%Damp, 'Damping Matrix', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)
        else if (index(Line,'!WEIGHT CONSTANT')==1) then
            if (p%nTot<0) exit
            CALL AllocAry(p%W0, p%nTot, 'W0', ErrStat, ErrMsg ); if (Failed()) return
            call ReadAry(UnIn, InputFile, p%W0, p%nTot, 'W0', 'Weight constant', ErrStat, ErrMsg)
        else if (index(Line,'!WEIGHT STIFFNESS')==1) then
            if (p%nTot<0) exit
            call ReadRealMatrix(UnIn, InputFile, p%WStff, 'Weight stiffness', p%nTot, p%nTot, ErrStat, ErrMsg, iLine)
        ! elseif (index(Line,'!')==1) then
            !write(*,*) 'Ignored comment: '//trim(Line)
        ! else
            ! Ignore unsupported lines
            !write(*,*) 'Ignored line: '//trim(Line)
        endif
    enddo
    close( UnIn )

CONTAINS
    logical function Failed()
        CALL SetErrStatSimple(ErrStat, ErrMsg, 'ExtPtfm_ReadReducedFile')
        Failed =  ErrStat >= AbortErrLev
        if(Failed) call cleanup()
    end function Failed
    subroutine cleanup()
        close( UnIn )
        if (allocated(TmpAry)) deallocate(TmpAry)
    end subroutine cleanup
END SUBROUTINE ReadReducedFile

!..................................................................................................................................
SUBROUTINE ReadConnFile( InputFile, p, ErrStat, ErrMsg )
!..................................................................................................................................
   ! Passed variables
   CHARACTER(*),                INTENT(IN)    :: InputFile                           !< Name of the file containing the primary input data
   TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p                                   !< All the parameter matrices stored in this input file
   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                             !< Error status
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                              !< Error message
   ! Local variables:
   INTEGER(IntKi)                       :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)                       :: iLine                                     ! Current position in file
   CHARACTER(4096)                      :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Get an available unit number for the file.
   CALL GetNewUnit( UnIn, ErrStat, ErrMsg );            if (Failed()) return
   ! Open the Primary input file.
   CALL OpenFInpFile(UnIn, InputFile, ErrStat, ErrMsg); if (Failed()) return

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
       ! Line content is analyzed as case insensitive
       call Conv2UC(Line)
       if (index(Line,'!NCONN:')==1) then
           p%nConn = ReadIntFromStr(Line(8:), '`Nconn`, file '//trim(InputFile)//', line '//Num2LStr(iLine), ErrStat, ErrMsg); if (Failed()) return
           if (p%nConn<=0_IntKi) return
        else if (index(Line,'!CONNECTIONS')==1) then
           call ReadRealMatrix(UnIn, InputFile, p%PosConn, 'Connections', p%nConn, 3_IntKi, ErrStat, ErrMsg, iLine)
        else if (index(Line,'!DISPLACEMENT')==1) then
           call ReadRealMatrix(UnIn, InputFile, p%PhiConn, 'Displacement', 3*p%nConn, p%nTot, ErrStat, ErrMsg, iLine)
        ! else if (index(Line,'!')==1) then
           !write(*,*) 'Ignored comment: '//trim(Line)
        ! else
           ! Ignore unsupported lines
           !write(*,*) 'Ignored line: '//trim(Line)
        end if
    end do
    close( UnIn )

CONTAINS
    logical function Failed()
        CALL SetErrStatSimple(ErrStat, ErrMsg, 'ExtPtfm_ReadConnFile')
        Failed =  ErrStat >= AbortErrLev
        if(Failed) call cleanup()
    end function Failed
    subroutine cleanup()
        close( UnIn )
    end subroutine cleanup
END SUBROUTINE ReadConnFile

!..................................................................................................................................
SUBROUTINE ReadForceFile( InputFile, p, ErrStat, ErrMsg )
!..................................................................................................................................
   ! Passed variables
   CHARACTER(*),                INTENT(IN)    :: InputFile                           !< Name of the file containing the primary input data
   TYPE(ExtPtfm_ParameterType), INTENT(INOUT) :: p                                   !< All the parameter matrices stored in this input file
   INTEGER(IntKi),              INTENT(OUT)   :: ErrStat                             !< Error status
   CHARACTER(*),                INTENT(OUT)   :: ErrMsg                              !< Error message
   ! Local variables:
   REAL(ReKi), dimension(:),allocatable :: TmpAry                                 ! temporary array for reading row from file
   INTEGER(IntKi)                       :: I                                         ! loop counter
   INTEGER(IntKi)                       :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)                       :: iLine                                     ! Current position in file
   CHARACTER(4096)                      :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
   LOGICAL                              :: foundNSteps
   LOGICAL                              :: foundForcing
   ErrStat = ErrID_None
   ErrMsg  = ""

   foundNSteps  = .false.
   foundForcing = .false.

   p%nTimeSteps = 0_IntKi

   ! Get an available unit number for the file.
   CALL GetNewUnit( UnIn, ErrStat, ErrMsg );            if (Failed()) return
   ! Open the Primary input file.
   CALL OpenFInpFile(UnIn, InputFile, ErrStat, ErrMsg); if (Failed()) return

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
        ! Line content is analyzed as case insensitive
        call Conv2UC(Line)
        if (index(Line,'!NSTEPS')==1) then
            foundNSteps = .true.
            p%nTimeSteps =  ReadIntFromStr(Line(9:), '`Nsteps`, file '//trim(InputFile)//', line '//Num2LStr(iLine), ErrStat, ErrMsg); if (Failed()) return
        else if (index(Line,'!FORCING')==1) then
            foundForcing = .true.
            if (p%nTot<0 .or. p%nTimeSteps==0) exit
            call allocAry( p%Forces, max(1,p%nTimeSteps), p%nTot, 'p%Forces', ErrStat, ErrMsg); if (Failed()) return
            call allocAry( p%times , max(1,p%nTimeSteps),         'p%times' , ErrStat, ErrMsg); if (Failed()) return
            allocate(TmpAry(1:p%nTot+1))
            do i=1,p%nTimeSteps
                iLine=iLine+1
                TmpAry(1:p%nTot+1)=-999.9E-09
                read(UnIn, fmt='(A)', iostat=ErrStat) Line
                if (ErrStat/=0) then
                ErrStat = ErrID_Fatal
                ErrMSg='Failed to read line '//trim(Num2LStr(iLine))//' (out of '//trim(Num2LStr(p%nTimeSteps))//' expected lines) in file: '//trim(InputFile)
                exit
                end if
                ! Extract fields (ReadR8AryFromStr is in NWTC_IO)
                CALL ReadAry(Line, TmpAry, p%nTot+1, 'Forces', 'Forces', ErrStat, ErrMsg)
                if (ErrStat/=0) then
                ErrStat = ErrID_Fatal
                ErrMsg='Failed to extract fields from line '//trim(Num2LStr(iLine))//'. '//trim(ErrMsg)//'. Check that the number of columns is correct in file: '//trim(InputFile)
                exit
                end if
                if (ErrStat /= 0) exit
                p%times(i)    = TmpAry(1)
                p%Forces(i,:) = TmpAry(2:p%nTot+1)
            end do
        ! elseif (index(Line,'!')==1) then
            !write(*,*) 'Ignored comment: '//trim(Line)
        ! else
            ! Ignore unsupported lines
            !write(*,*) 'Ignored line: '//trim(Line)
        end if
    enddo
    close( UnIn )

    if (.not.foundNSteps) then
        ErrStat = ErrID_Fatal
        ErrMsg  = "Did not find '!NSteps:' followed by the number of time steps in file "//trim(InputFile)//"."
        if(Failed()) return
    end if

    if (.not.foundForcing) then
        ErrStat = ErrID_Fatal
        ErrMsg  = "Did not find forcing time series after '!FORCING' in file "//trim(InputFile)//". Note that the time series should be after '!NSteps:'."
        if(Failed()) return
    end if

    if (p%nTimeSteps <= 0_IntKi) then
        p%nTimeSteps = 1_IntKi
        call allocAry( p%Forces, 1_IntKi, p%nTot, 'p%Forces', ErrStat, ErrMsg); if(Failed()) return
        call allocAry( p%times , 1_IntKi,         'p%times' , ErrStat, ErrMsg); if(Failed()) return
        p%Forces= 0.0_ReKi
        p%times = 0.0_ReKi
    end if

CONTAINS
    logical function Failed()
        CALL SetErrStatSimple(ErrStat, ErrMsg, 'ExtPtfm_ReadForceFile')
        Failed =  ErrStat >= AbortErrLev
        if(Failed) call cleanup()
    end function Failed
    subroutine cleanup()
        close( UnIn )
        if (allocated(TmpAry)) deallocate(TmpAry)
    end subroutine cleanup
END SUBROUTINE ReadForceFile

!> This routine generates the summary file, which contains a regurgitation of  the input data and interpolated flexible body data.
SUBROUTINE ExtPtfm_PrintSum(x, p, m, RootName, ErrStat, ErrMsg)
   ! passed variables
   TYPE(ExtPtfm_ContinuousStateType), INTENT(IN)  :: x           !< Initial continuous states
   TYPE(ExtPtfm_ParameterType),    INTENT(IN   )  :: p           !< Parameters of the structural dynamics module
   TYPE(ExtPtfm_MiscVarType),      INTENT(IN   )  :: m           !< Misc variables for optimization (not copied in glue code)
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
   ! TODO TODO TODO YAML FORMAT
   ! TODO TODO TODO ONLY Open Summary if no optional unit given 

   ! Open the summary file and give it a heading.
   CALL GetNewUnit(UnSu, ErrStat, ErrMsg);
   CALL OpenFOutFile(UnSu, TRIM( RootName )//'.sum', ErrStat, ErrMsg)
   IF ( ErrStat /= ErrID_None ) RETURN
   ! Heading:
   WRITE (UnSu,'(/,A)')  '!This summary information was generated by '//TRIM( GetNVD(ExtPtfm_Ver) )// &
                         ' on '//CurDate()//' at '//CurTime()//'.'

   write(UnSu,'(A)')      '!Module input file'
   write(UnSu,'(A,A)')    'Time integration method      : ',StrIntMethod(p%IntMethod)
   write(UnSu,'(A,F13.8)')'Integration time step        : ',p%EP_DeltaT
   write(UnSu,'(A)')      '!Reduction input file'
   write(UnSu,'(A,I0)')   'Number of time steps         : ',p%nTimeSteps
   write(UnSu,'(A,F13.8)')'Start time                   : ',p%times(1)
   write(UnSu,'(A,F13.8)')'End time                     : ',p%times(p%nTimeSteps)
   write(UnSu,'(A,I0)')   'Total number of DOF (input)  : ',p%nCBFull+6
   write(UnSu,'(A,I0)')   'Number of CB modes (input)   : ',p%nCBFull
   write(UnSu,'(A)')      '!Degrees of freedom'
   write(UnSu,'(A,I0)')   'Total number of DOF (active) : ',p%nTot
   write(UnSu,'(A,I0)')   'Number of CB modes (active)  : ',p%nCB
   call disp1i(UnSu, 'ActiveCBDOF',p%ActiveCBDOF)
! 
   if (m%EquilStart) then
       write(UnSu,'(A)')'!Initial conditions (before equilibrium)'
   else
       write(UnSu,'(A)')'!Initial conditions (no equilibrium will be computed)'
   endif
   call disp1r8(UnSu, 'qm'   ,x%qm)
   call disp1r8(UnSu, 'qmdot',x%qmdot)

   !write(UnSu,'(A)')'!State matrices'
   !call disp2r8(UnSu, 'A',p%AMat)
   !call disp2r8(UnSu, 'B',p%BMat)
   !call disp2r8(UnSu, 'C',p%CMat)
   !call disp2r8(UnSu, 'D',p%DMat)
   write(UnSu,'(A)')'!Input matrices'
   call disp2r8(UnSu, 'M',p%Mass)
   call disp2r8(UnSu, 'K',p%Stff)
   call disp2r8(UnSu, 'C',p%Damp)
!    call disp2r8(UnSu, 'F',p%Forces)
   write(UnSu,'(A)')'!Input sub-matrices'
   call disp2r8(UnSu, 'M11',p%M11)
   call disp2r8(UnSu, 'M12',p%M12)
   call disp2r8(UnSu, 'M21',p%M21)
   call disp2r8(UnSu, 'M22',p%M22)
   call disp2r8(UnSu, 'K11',p%K11)
   call disp2r8(UnSu, 'K22',p%K22)
   call disp2r8(UnSu, 'C11',p%C11)
   call disp2r8(UnSu, 'C12',p%C12)
   call disp2r8(UnSu, 'C21',p%C21)
   call disp2r8(UnSu, 'C22',p%C22)

   write(UnSu,'(//,A,//)')  '!Connections:'
   write(UnSu,'(A,I0)')     'Number of connections        : ',p%nConn
   call disp2r8(UnSu, 'PosConn',p%PosConn)
   call disp2r8(UnSu, 'PhiConn',p%PhiConn)

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
