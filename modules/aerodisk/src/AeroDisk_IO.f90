!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2024  National Renewable Energy Laboratory
!
!    This file is part of AeroDisk
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
MODULE AeroDisk_IO

   USE AeroDisk_Types
   USE AeroDisk_Output_Params
   USE NWTC_Library

   implicit none

   ! data storage for reading/parsing table
   type :: TableIndexType
      integer(IntKi)    :: ColTSR   !< Column number for Tip-Speed Ratio
      integer(IntKi)    :: ColRtSpd !< Column number for Rotor Speed
      integer(IntKi)    :: ColVRel  !< Column number for VRel
      integer(IntKi)    :: ColSkew  !< Column number for Skew
      integer(IntKi)    :: ColPitch !< Column number for Pitch
      integer(IntKi)    :: NumColNamesGiven  !< total number of column names given
   end type TableIndexType

contains

!---------------------------------------------------------------
!> Parse the input in the InFileInfo (FileInfo_Type data structure):
subroutine ADsk_ParsePrimaryFileData( InitInp, RootName, interval, FileInfo_In, InputFileData, UnEc, ErrStat, ErrMsg )
   type(ADsk_InitInputType),  intent(in   )  :: InitInp              !< Input data for initialization routine
   character(1024),           intent(in   )  :: RootName             !< root name for summary file
   real(DBKi),                intent(in   )  :: interval             !< timestep
   type(FileInfoType),        intent(in   )  :: FileInfo_In          !< The input file stored in a data structure
   type(ADsk_InputFile),      intent(inout)  :: InputFileData        !< The data for initialization
   integer(IntKi),            intent(  out)  :: UnEc                 !< The local unit number for this module's echo file
   integer(IntKi),            intent(  out)  :: ErrStat              !< Error status  from this subroutine
   character(*),              intent(  out)  :: ErrMsg               !< Error message from this subroutine

   ! local vars
   integer(IntKi)                            :: CurLine              !< current entry in FileInfo_In%Lines array
   integer(IntKi)                            :: i                    !< generic counter
   type(TableIndexType)                      :: TabIdx               !< indices for table columnns, for simplifying data parsing/passing
   real(SiKi)                                :: TmpRe(10)            !< temporary 10 number array for reading values in from table
   integer(IntKi)                            :: ErrStat2             !< Temporary error status  for subroutine and function calls
   character(ErrMsgLen)                      :: ErrMsg2              !< Temporary error message for subroutine and function calls
   character(*),              parameter      :: RoutineName="ADsk_ParsePrimaryFileData"

      ! Initialize ErrStat
   ErrStat  = ErrID_None
   ErrMsg   = ""
   UnEc     = -1  ! No file


   CALL AllocAry( InputFileData%OutList, MaxOutPts, "Outlist", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !======  General  ====================================================================================
   CurLine = 4    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2 )
         if (Failed()) return;

   if ( InputFileData%Echo ) then
      CALL OpenEcho ( UnEc, TRIM(RootName)//'.ech', ErrStat2, ErrMsg2 )
         if (Failed()) return;
      WRITE(UnEc, '(A)') 'Echo file for AeroDisk primary input file: '//trim(InitInp%InputFile)
      ! Write the first three lines into the echo file
      WRITE(UnEc, '(A)') FileInfo_In%Lines(1)
      WRITE(UnEc, '(A)') FileInfo_In%Lines(2)
      WRITE(UnEc, '(A)') FileInfo_In%Lines(3)

      CurLine = 4
      call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return
   endif


      ! DeltaT - Time interval for aerodynamic calculations {or default} (s):
   call ParseVarWDefault ( FileInfo_In, CurLine, "DT", InputFileData%DT, interval, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

   !======  Environmental Conditions  ===================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! AirDens - Air density {or default} (kg/m^3)
   call ParseVarWDefault( FileInfo_In, CurLine, "AirDens", InputFileData%AirDens, InitInp%defAirDens, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

   !======  Actuator Disk Properties  ===================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! RotorRad - Rotor radius (m) (or "default")
   call ParseVarWDefault( FileInfo_In, CurLine, "RotorRad", InputFileData%RotorRad, InitInp%RotorRad, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! InColNames - names for the input columns for the table.
   call Get_InColNames( FileInfo_In, CurLine, TabIdx, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! InColDims - Number of values in each column (-) (must have same number of columns as InColName) [each >=2]
   call Get_InColDims( FileInfo_In, CurLine, TabIdx, InputFileData%AeroTable, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! Column headers
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)
   CurLine = CurLine + 1

      ! Read table
   call Get_RtAeroTableData( FileInfo_In, CurLine, TabIdx, InputFileData%AeroTable, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return


   !======  Outputs  ====================================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
!      ! SumPrint - Generate a summary file listing input options and interpolated properties to "<rootname>.AD.sum"?  (flag)
!   call ParseVar( FileInfo_In, CurLine, "SumPrint", InputFileData%SumPrint, ErrStat2, ErrMsg2, UnEc )
!      if (Failed()) return

   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
   call ReadOutputListFromFileInfo( FileInfo_In, CurLine, InputFileData%OutList, &
            InputFileData%NumOuts, ErrStat2, ErrMsg2, UnEc )
         if (Failed()) return;

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed
   subroutine Cleanup()
      ! Only do this on a fault.  Leave open for calling routine in case we want to write anything else.
      if (UnEc > 0_IntKi)  close(UnEc)
   end subroutine Cleanup
   subroutine GetVarNamePos(FileName,LineNo,ChAry,NameToCheck,ColNum,ErrStat3,ErrMsg3)
      character(*),           intent(in   )  :: FileName
      integer(IntKi),         intent(in   )  :: LineNo
      character(*),           intent(in   )  :: ChAry(:)
      character(*),           intent(in   )  :: NameToCheck
      integer(IntKi),         intent(  out)  :: ColNum
      integer(IntKi),         intent(  out)  :: ErrStat3
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3
      character(len(ChAry(1)))               :: TmpNm
      character(len=len(NameToCheck))        :: TmpNmCk
      integer(IntKi)                         :: WordNum
      logical                                :: TmpFlag
      ErrStat3 = ErrID_None
      ErrMsg3  = ""
      TmpFlag  = .false.
      ColNum   = 0_IntKi      ! if not present
      TmpNmCk  = NameToCheck; call Conv2UC(TmpNmCk)
      do WordNum=1,size(ChAry)
         TmpNm = ChAry(WordNum); call Conv2UC(TmpNm)
         if (trim(TmpNm) == TmpNmCk) then
            if (TmpFlag) then
               ErrStat3 = ErrID_Fatal
               ErrMsg3  = NewLine//' >> A fatal error occurred when parsing data from '// &
                    trim(FileName)//'.'//NewLine// &
                    ' >> The column name '//trim(NameToCheck)//' occurs more than once on line '// &
                    trim(Num2LStr(LineNo))//'.'
               return
            endif
            TmpFlag = .true.
            ColNum = WordNum
         endif
      enddo
   end subroutine GetVarNamePos
   subroutine Get_InColNames(Info,LineNo,Idx,ErrStat3,ErrMsg3,UnEc)
      ! A custom routine is used here rather than using the ParseChAry as ParseChAry does not
      ! handle an unknown number of quoted strings well (i.e: "RtSpd,VRel,Skew,Pitch")
      type(FileInfoType),     intent(in   )  :: Info
      integer(IntKi),         intent(inout)  :: LineNo
      type(TableIndexType),   intent(inout)  :: Idx
      integer(IntKi),         intent(  out)  :: ErrStat3
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3
      integer(IntKi),         intent(in   )  :: UnEc
      integer(IntKi)                         :: WordPos
      character(1024)                        :: thisFile ! Simplify some error management
      integer(IntKi)                         :: thisLine ! Simplify some error management
      character(20)                          :: TmpChAry(6)
      ErrStat3 = ErrID_None
      ErrMsg3  = ""
      ! Echo if we have a file
      if ( UnEc > 0 )  write (UnEc,'(A)')  TRIM( Info%Lines(LineNo) )

      ! Parse out words from the line
      call GetWords( Info%Lines(LineNo), TmpChAry, size(TmpChAry))

      ! file info for error handling -- note that this could be in a different file than the above!!!!
      thisFile =trim(Info%FileList(Info%FileIndx(LineNo)))
      thisLine =     Info%FileLine(LineNo)

      ! Check that this is the right line in the file and InColNames exists in first 5 words of line
      call GetVarNamePos(thisFile, thisLine, TmpChAry, "InColNames", WordPos, ErrStat3, ErrMsg3)   ! ignore duplicate entry error here
         ! Error handling if wrong row
         if (WordPos <= 0_IntKi) then
            ErrStat3 = ErrID_Fatal
            ErrMsg3  = NewLine//' >> A fatal error occurred when parsing data from '// &
                       trim(thisFile)//'.'//NewLine// &
                       ' >> The variable "InColNames" was not found on line '//trim(Num2LStr(thisLine))//'.'
            return
         endif

      ! Get order of the other columns (returns 0 if column name not in file)
      call GetVarNamePos(thisFile, thisLine, TmpChAry, "TSR",   Idx%ColTSR,   ErrStat3,ErrMsg3); if (ErrStat3 >= ErrID_Fatal) return
      call GetVarNamePos(thisFile, thisLine, TmpChAry, "RtSpd", Idx%ColRtSpd, ErrStat3,ErrMsg3); if (ErrStat3 >= ErrID_Fatal) return
      call GetVarNamePos(thisFile, thisLine, TmpChAry, "VRel",  Idx%ColVRel,  ErrStat3,ErrMsg3); if (ErrStat3 >= ErrID_Fatal) return
      call GetVarNamePos(thisFile, thisLine, TmpChAry, "Skew",  Idx%ColSkew,  ErrStat3,ErrMsg3); if (ErrStat3 >= ErrID_Fatal) return
      call GetVarNamePos(thisFile, thisLine, TmpChAry, "Pitch", Idx%ColPitch, ErrStat3,ErrMsg3); if (ErrStat3 >= ErrID_Fatal) return

      ! total number of columns specified
      Idx%NumColNamesGiven = maxval( (/ Idx%ColTSR, Idx%ColRtSpd, Idx%ColVRel, Idx%ColSkew, Idx%ColPitch /) )

      ! make sure have have column names
      if ((Idx%ColTSR + Idx%ColRtSpd + Idx%ColVRel + Idx%ColSkew + Idx%ColPitch) <= 0_IntKi) then
         ErrStat3 = ErrID_Fatal
         ErrMsg3  =  NewLine//' >> A fatal error occurred when parsing data from '// &
               trim(thisFile)//'.'//NewLine// &
               ' >> At least one column of data named "TSR", "RtSpd", "VRel", "Skew", or "Pitch" must exist'// &
               ' in the table header, but none were found'//trim(Num2LStr(thisLine))//'.'
         return
      endif
      LineNo = LineNo + 1      ! Picked up column names, so increment to next line and return
      return
   end subroutine Get_InColNames
   subroutine Get_InColDims(Info,LineNo,Idx,AeroTable,ErrStat3,ErrMsg3,UnEc)
      ! A custom routine is used here rather than using the ParseChAry as ParseChAry does not
      ! handle an unknown number of quoted strings well (i.e: "RtSpd,VRel,Skew,Pitch")
      type(FileInfoType),     intent(in   )  :: Info
      integer(IntKi),         intent(inout)  :: LineNo
      type(TableIndexType),   intent(in   )  :: Idx
      type(ADsk_AeroTable),   intent(inout)  :: AeroTable
      integer(IntKi),         intent(  out)  :: ErrStat3
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3
      integer(IntKi),         intent(in   )  :: UnEc
      integer(IntKi)                         :: WordPos
      integer(IntKi)                         :: IOS
      character(1024)                        :: thisFile    ! Simplify some error management
      integer(IntKi)                         :: thisLine    ! Simplify some error management
      character(20)                          :: TmpChAry(6) ! Assume 5 columns and the name
      ErrStat3 = ErrID_None
      ErrMsg3  = ""
      ! Echo if we have a file
      if ( UnEc > 0 )  write (UnEc,'(A)')  TRIM( Info%Lines(LineNo) )

      ! Parse out words from the line
      call GetWords( Info%Lines(LineNo), TmpChAry, size(TmpChAry))

      ! file info for error handling -- note that this could be in a different file than the above!!!!
      thisFile =trim(Info%FileList(Info%FileIndx(LineNo)))
      thisLine =     Info%FileLine(LineNo)

      ! Check that this is the right line in the file and InColDims exists in first 5 words of line
      call GetVarNamePos(thisFile, thisLine, TmpChAry, "InColDims", WordPos, ErrStat3, ErrMsg3)   ! ignore duplicate entry error here
         ! Error handling if wrong row
         if (WordPos <= 0_IntKi) then
            ErrStat3 = ErrID_Fatal
            ErrMsg3  = NewLine//' >> A fatal error occurred when parsing data from '// &
                       trim(thisFile)//'.'//NewLine// &
                       ' >> The variable "InColDims" was not found on line '//trim(Num2LStr(thisLine))//'.'
            return
         endif

      ! set number of indices for each to 0
      AeroTable%N_TSR   = 0
      AeroTable%N_RtSpd = 0
      AeroTable%N_VRel  = 0
      AeroTable%N_Skew  = 0
      AeroTable%N_Pitch = 0

      ! Read numbers for number of various entries from the tmpChAry
      if (Idx%ColTSR    > 0_IntKi) then
         READ (TmpChAry(Idx%ColTSR  ),*,IOSTAT=IOS)   AeroTable%N_TSR
         if (IOS /= 0) then
            call SetErrStat(ErrID_Fatal,'Could not read "N_TSR" from column '//     &
               trim(Num2LStr(Idx%ColTSR))//' on line '//trim(Num2LStr(thisLine))//  &
               ' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
            return
         endif
      endif
      if (Idx%ColRtSpd  > 0_IntKi) then
         READ (TmpChAry(Idx%ColRtSpd),*,IOSTAT=IOS)   AeroTable%N_RtSpd
         if (IOS /= 0) then
            call SetErrStat(ErrID_Fatal,'Could not read "N_RtSpd" from column '//     &
               trim(Num2LStr(Idx%ColRtSpd))//' on line '//trim(Num2LStr(thisLine))//  &
               ' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
            return
         endif
      endif
      if (Idx%ColVRel   > 0_IntKi) then
         READ (TmpChAry(Idx%ColVRel ),*,IOSTAT=IOS)   AeroTable%N_VRel
         if (IOS /= 0) then
            call SetErrStat(ErrID_Fatal,'Could not read "N_VRel" from column '//     &
               trim(Num2LStr(Idx%ColVRel))//' on line '//trim(Num2LStr(thisLine))//  &
               ' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
            return
         endif
      endif
      if (Idx%ColSkew   > 0_IntKi) then
         READ (TmpChAry(Idx%ColSkew ),*,IOSTAT=IOS)   AeroTable%N_Skew
         if (IOS /= 0) then
            call SetErrStat(ErrID_Fatal,'Could not read "N_Skew" from column '//     &
               trim(Num2LStr(Idx%ColSkew))//' on line '//trim(Num2LStr(thisLine))//  &
               ' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
            return
         endif
      endif
      if (Idx%ColPitch  > 0_IntKi) then
         READ (TmpChAry(Idx%ColPitch),*,IOSTAT=IOS)   AeroTable%N_Pitch
         if (IOS /= 0) then
            call SetErrStat(ErrID_Fatal,'Could not read "N_Pitch" from column '//     &
               trim(Num2LStr(Idx%ColPitch))//' on line '//trim(Num2LStr(thisLine))//  &
               ' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
            return
         endif
      endif
      ! make sure all values positive
      if (AeroTable%N_TSR < 0_IntKi)   call SetErrStat(ErrID_Fatal,'Entry for "N_TSR" must be postive valued from column '//     &
               trim(Num2LStr(Idx%ColTSR  ))//' on line '//trim(Num2LStr(thisLine))//' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
      if (AeroTable%N_RtSpd < 0_IntKi)   call SetErrStat(ErrID_Fatal,'Entry for "N_RtSpd" must be postive valued from column '//     &
               trim(Num2LStr(Idx%ColRtSpd))//' on line '//trim(Num2LStr(thisLine))//' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
      if (AeroTable%N_VRel < 0_IntKi)   call SetErrStat(ErrID_Fatal,'Entry for "N_VRel" must be postive valued from column '//     &
               trim(Num2LStr(Idx%ColVRel ))//' on line '//trim(Num2LStr(thisLine))//' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
      if (AeroTable%N_Skew < 0_IntKi)   call SetErrStat(ErrID_Fatal,'Entry for "N_Skew" must be postive valued from column '//     &
               trim(Num2LStr(Idx%ColSkew ))//' on line '//trim(Num2LStr(thisLine))//' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
      if (AeroTable%N_Pitch < 0_IntKi)   call SetErrStat(ErrID_Fatal,'Entry for "N_Pitch" must be postive valued from column '//     &
               trim(Num2LStr(Idx%ColPitch))//' on line '//trim(Num2LStr(thisLine))//' in file '//trim(thisFile)//'.',ErrStat3,ErrMsg3,'')
      ! NOTE:  we are storing 0 for the dimensions that don't exist in the table. We will
      !        modify this later to make table usage simpler
      LineNo = LineNo + 1      ! Picked up column names, so increment to next line and return
      return
   end subroutine Get_InColDims

end subroutine ADsk_ParsePrimaryFileData


subroutine Get_RtAeroTableData(Info,LineNo,Idx,AeroTable,ErrStat,ErrMsg,UnEc)
   ! Table entries may not be in order.  So sort while reading in.
   ! NOTE: Sparse data files are not currently supported
   type(FileInfoType),     intent(in   )  :: Info
   integer(IntKi),         intent(inout)  :: LineNo
   type(TableIndexType),   intent(in   )  :: Idx
   type(ADsk_AeroTable),   intent(inout)  :: AeroTable
   integer(IntKi),         intent(  out)  :: ErrStat
   character(ErrMsgLen),   intent(  out)  :: ErrMsg
   integer(IntKi),         intent(in   )  :: UnEc
   integer(IntKi)                         :: WordPos
   integer(IntKi)                         :: IOS
   character(1024)                        :: thisFile    ! Simplify some error management
   integer(IntKi)                         :: thisLine    ! Simplify some error management
   integer(IntKi)                         :: i,j         ! Generic counters
   integer(IntKi)                         :: NumCols     ! number of columns expected
   integer(IntKi)                         :: NumRows     ! number of rows expected
   integer(IntKi)                         :: Sz(5)       ! Array size -- for readability of code
   real(SiKi),             allocatable    :: TmpTab(:,:)       ! temporary real array of allocatable size -- will read entire table in, then do the sorting
   logical,                allocatable    :: Mask(:,:,:,:,:)   ! to make sure we aren't missing any terms
   integer(IntKi)                         :: ErrStat2    !< Temporary error status  for subroutine and function calls
   character(ErrMsgLen)                   :: ErrMsg2     !< Temporary error message for subroutine and function calls
   character(*),           parameter      :: RoutineName='Get_RtAeroTableData'
   ErrStat  = ErrID_None
   ErrMsg   = ""


   ! Number of columns in table
   Sz = 0_IntKi
   if (AeroTable%N_TSR   > 0_IntKi)    Sz(1) = 1_IntKi
   if (AeroTable%N_RtSpd > 0_IntKi)    Sz(2) = 1_IntKi
   if (AeroTable%N_VRel  > 0_IntKi)    Sz(3) = 1_IntKi
   if (AeroTable%N_Pitch > 0_IntKi)    Sz(4) = 1_IntKi
   if (AeroTable%N_Skew  > 0_IntKi)    Sz(5) = 1_IntKi
   NumCols = Idx%NumColNamesGiven + 6_IntKi      ! Add DOF columns

   ! temporary array for sizing -- note that min dimension size is 1 so we calculate number of rows correctly
   Sz(1) = max(AeroTable%N_TSR,  1_IntKi)
   Sz(2) = max(AeroTable%N_RtSpd,1_IntKi)
   Sz(3) = max(AeroTable%N_VRel, 1_IntKi)
   Sz(4) = max(AeroTable%N_Pitch,1_IntKi)
   Sz(5) = max(AeroTable%N_Skew, 1_IntKi)

   ! Total number of rows we expect to find in the table
   NumRows= Sz(1) * Sz(2) * Sz(3) * Sz(4) * Sz(5)

   ! Allocate arrays for forces/moments
   call AllocAry(AeroTable%C_Fx, Sz(1), Sz(2), Sz(3), Sz(4), Sz(5), 'AeroTable%C_Fx', ErrStat2, ErrMsg2);   if (Failed())  return;  AeroTable%C_Fx = 0.0_SiKi
   call AllocAry(AeroTable%C_Fy, Sz(1), Sz(2), Sz(3), Sz(4), Sz(5), 'AeroTable%C_Fy', ErrStat2, ErrMsg2);   if (Failed())  return;  AeroTable%C_Fy = 0.0_SiKi
   call AllocAry(AeroTable%C_Fz, Sz(1), Sz(2), Sz(3), Sz(4), Sz(5), 'AeroTable%C_Fz', ErrStat2, ErrMsg2);   if (Failed())  return;  AeroTable%C_Fz = 0.0_SiKi
   call AllocAry(AeroTable%C_Mx, Sz(1), Sz(2), Sz(3), Sz(4), Sz(5), 'AeroTable%C_Mx', ErrStat2, ErrMsg2);   if (Failed())  return;  AeroTable%C_Mx = 0.0_SiKi
   call AllocAry(AeroTable%C_My, Sz(1), Sz(2), Sz(3), Sz(4), Sz(5), 'AeroTable%C_My', ErrStat2, ErrMsg2);   if (Failed())  return;  AeroTable%C_My = 0.0_SiKi
   call AllocAry(AeroTable%C_Mz, Sz(1), Sz(2), Sz(3), Sz(4), Sz(5), 'AeroTable%C_Mz', ErrStat2, ErrMsg2);   if (Failed())  return;  AeroTable%C_Mz = 0.0_SiKi

   ! Allocate a mask for data checks
   allocate(Mask(Sz(1),Sz(2),Sz(3),Sz(4),Sz(5)),STAT=ErrStat2)
   if (ErrStat2 /= 0) then
      call SetErrStat(ErrID_Fatal,'Could not allocate array for data mask',ErrStat,ErrMsg,RoutineName)
      call Cleanup()
      return
   endif
   Mask = .false.

   ! Slurp the entire table into real array -- order for more efficient searching (bad for reading in)
   call AllocAry(TmpTab,NumRows,NumCols,'TemporaryTable',ErrStat2,ErrMsg2); if (Failed())  return
   do i=1,NumRows
      read (Info%Lines(LineNo),*,iostat=IOS)    TmpTab(i,1:NumCols)
      if (IOS /= 0_IntKi) then
         thisFile =trim(Info%FileList(Info%FileIndx(LineNo)))
         thisLine =     Info%FileLine(LineNo)
         call SetErrStat(ErrID_Fatal,'Could not read "Actuator Disk Properties Table" row '//trim(Num2LStr(i))//' (line '//trim(Num2LStr(thisLine))//  &
                  ' in file '//trim(thisFile)//').  Expecting '//trim(Num2LStr(NumRows))//' rows and '//trim(Num2LStr(NumCols))//' columns in table.',ErrStat,ErrMsg,RoutineName)
         call Cleanup()
         return
      endif
      if (UnEc > 0_IntKi)  write(UnEc,'(A)') trim(Info%Lines(LineNo))
      LineNo = LineNo + 1
   enddo

   ! Find the unique values in the indexing columns of the table (if name given, otherwise skip)
   if (AeroTable%N_TSR   > 0_IntKi) call GetTabIndexVals( TmpTab(:,Idx%ColTSR  ),'TSR'  ,AeroTable%N_TSR  ,AeroTable%TSR  , ErrStat2, ErrMsg2); if (Failed()) return
   if (AeroTable%N_RtSpd > 0_IntKi) call GetTabIndexVals( TmpTab(:,Idx%ColRtSpd),'RtSpd',AeroTable%N_RtSpd,AeroTable%RtSpd, ErrStat2, ErrMsg2); if (Failed()) return
   if (AeroTable%N_VRel  > 0_IntKi) call GetTabIndexVals( TmpTab(:,Idx%ColVRel ),'VRel' ,AeroTable%N_VRel ,AeroTable%VRel , ErrStat2, ErrMsg2); if (Failed()) return
   if (AeroTable%N_Pitch > 0_IntKi) call GetTabIndexVals( TmpTab(:,Idx%ColPitch),'Pitch',AeroTable%N_Pitch,AeroTable%Pitch, ErrStat2, ErrMsg2); if (Failed()) return
   if (AeroTable%N_Skew  > 0_IntKi) call GetTabIndexVals( TmpTab(:,Idx%ColSkew ),'Skew' ,AeroTable%N_Skew ,AeroTable%Skew , ErrStat2, ErrMsg2); if (Failed()) return
   
   ! Now populate matrix -- read each line and put in correct table entry location
   call PopulateAeroTabs(AeroTable,Mask,Idx,TmpTab,NumRows,NumCols,ErrStat2,ErrMsg2); if (Failed()) return
   call CheckAeroTabs(AeroTable,Mask,ErrStat2,ErrMsg2);    if (Failed()) return

   ! Now convert RtSpd from rpm to rad/s, and Pitch and Skew from deg to rad
   if (AeroTable%N_RtSpd > 0_IntKi)    AeroTable%RtSpd = (AeroTable%RtSpd * Pi_S)/30.0_SiKi
   if (AeroTable%N_Pitch > 0_IntKi)    AeroTable%Pitch = (AeroTable%Pitch * Pi_S)/180.0_SiKi
   if (AeroTable%N_Skew  > 0_IntKi)    AeroTable%Skew  = (AeroTable%Skew  * Pi_S)/180.0_SiKi

   call Cleanup()

contains
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed
   subroutine Cleanup()
      if (allocated(TmpTab))  deallocate(TmpTab)
      if (allocated(Mask))    deallocate(Mask)
   end subroutine Cleanup
   subroutine GetTabIndexVals(TabCol,ColName,NumExpect,UniqueArray,ErrStat3,ErrMsg3)
      real(SiKi),             intent(in   )  :: TabCol(:)
      character(*),           intent(in   )  :: ColName
      integer(IntKi),         intent(in   )  :: NumExpect
      real(SiKi), allocatable,intent(  out)  :: UniqueArray(:)
      integer(IntKi),         intent(  out)  :: ErrStat3
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3
      integer(IntKi)                         :: NumFound
      call UniqueRealValues( TabCol, UniqueArray, NumFound, ErrStat3, ErrMsg3 )
      if (NumExpect /= NumFound) then
         call SetErrStat(ErrID_Fatal,'Expecting '//trim(Num2LStr(NumExpect))//' unique '//ColName//   &
            ' entries in "Actuator Disk Properties Table", but found '//trim(Num2LStr(NumFound))//' instead.',ErrStat3,ErrMsg3,'')
      endif
   end subroutine GetTabIndexVals
   subroutine PopulateAeroTabs(Aero,DatMask,TabIdx,Dat,nRow,nCol,ErrStat3,ErrMsg3)
      type(ADsk_AeroTable),   intent(inout)  :: Aero
      logical,    allocatable,intent(inout)  :: DatMask(:,:,:,:,:)
      type(TableIndexType),   intent(in   )  :: TabIdx
      real(SiKi), allocatable,intent(in   )  :: Dat(:,:)
      integer(IntKi),         intent(in   )  :: nCol        ! number of columns in Dat
      integer(IntKi),         intent(in   )  :: nRow        ! number of rows in Dat
      integer(IntKi),         intent(  out)  :: ErrStat3
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3
      real(SiKi)                             :: TmpR6(6)
      integer(IntKi)                         :: row
      integer(IntKi)                         :: iTSR,iRtSpd,iVRel,iPitch,iSkew
      ErrStat3 = ErrID_None
      ErrMsg3  = ''
      DatMask = .false.    ! Make sure data mask is false.  Set to true for entries we populate
      ! Step through each row of raw data, and find corresponding data location
      ! Set index for placement to 1 in the event that the variable isn't actually used
      do row=1,nRow
         TmpR6(1:6) = Dat(row,nCol-5:nCol)   ! DOF values
         ! Find TSR entry        LocateStp(XVal, XAry, Ind, AryLen)
         iTSR   = 1_IntKi
         iRtSpd = 1_IntKi
         iVRel  = 1_IntKi
         iPitch = 1_IntKi
         iSkew  = 1_IntKi
         ! if no entries, the corresponding array is unallocated, so skip and leave index at 1
         if (Aero%N_TSR   > 0_IntKi) call LocateStp( Dat(row,TabIdx%ColTSR  ),Aero%TSR  ,iTSR  ,Aero%N_TSR  )
         if (Aero%N_RtSpd > 0_IntKi) call LocateStp( Dat(row,TabIdx%ColRtSpd),Aero%RtSpd,iRtSpd,Aero%N_RtSpd)
         if (Aero%N_VRel  > 0_IntKi) call LocateStp( Dat(row,TabIdx%ColVRel ),Aero%VRel ,iVRel ,Aero%N_VRel )
         if (Aero%N_Pitch > 0_IntKi) call LocateStp( Dat(row,TabIdx%ColPitch),Aero%Pitch,iPitch,Aero%N_Pitch)
         if (Aero%N_Skew  > 0_IntKi) call LocateStp( Dat(row,TabIdx%ColSkew ),Aero%Skew ,iSkew ,Aero%N_Skew )
         Aero%C_Fx(iTSR,iRtSpd,iVRel,iPitch,iSkew) = TmpR6(1)
         Aero%C_Fy(iTSR,iRtSpd,iVRel,iPitch,iSkew) = TmpR6(2)
         Aero%C_Fz(iTSR,iRtSpd,iVRel,iPitch,iSkew) = TmpR6(3)
         Aero%C_Mx(iTSR,iRtSpd,iVRel,iPitch,iSkew) = TmpR6(4)
         Aero%C_My(iTSR,iRtSpd,iVRel,iPitch,iSkew) = TmpR6(5)
         Aero%C_Mz(iTSR,iRtSpd,iVRel,iPitch,iSkew) = TmpR6(6)
         if (DatMask(iTSR,iRtSpd,iVRel,iPitch,iSkew)) then
            call SetErrStat(ErrID_Fatal,'Duplicate data entry in "Actuator Disk Properties Table" row '//trim(Num2LStr(row))//'.',ErrStat3,ErrMsg3,'')
         else
            DatMask(iTSR,iRtSpd,iVRel,iPitch,iSkew) = .true.
         endif
      enddo
   end subroutine PopulateAeroTabs
   subroutine CheckAeroTabs(Aero,DatMask,ErrStat3,ErrMsg3)
      type(ADsk_AeroTable),   intent(in   )  :: Aero
      logical,    allocatable,intent(in   )  :: DatMask(:,:,:,:,:)
      integer(IntKi),         intent(  out)  :: ErrStat3
      character(ErrMsgLen),   intent(  out)  :: ErrMsg3
      integer(IntKi)                         :: iTSR,iRtSpd,iVRel,iPitch,iSkew
      logical                                :: DataMiss
      logical                                :: FirstMiss
      character(60)                          :: TmpChar
      real(SiKi)                             :: TmpR5(5)
      ErrStat3 = ErrID_None
      ErrMsg3  = ''
      DataMiss = .false.
      FirstMiss = .false.
      do iSkew=1,size(DatMask,DIM=5)
         do iPitch=1,size(DatMask,DIM=4)
            do iVRel=1,size(DatMask,DIM=3)
               do iRtSpd=1,size(DatMask,DIM=2)
                  do iTSR=1,size(DatMask,DIM=1)
                     if (.not. DatMask(iTSR,iRtSpd,iVRel,iPitch,iSkew)) then
                        DataMiss = .true.
                        if (.not. FirstMiss) then
                           FirstMiss = .true.
                           ErrStat3 = ErrID_Fatal
                           ErrMsg3  = NewLine//'Data missing from table!!! (note order may not be same as file)'//NewLine//   &
                                 '   TSR      RtSpd     VRel      Pitch     Skew    '//NewLine
                        endif
                        TmpR5 = NaN_S
                        if (Aero%N_TSR  >0_IntKi)  TmpR5(1) = Aero%TSR  (iTSR  )
                        if (Aero%N_RtSpd>0_IntKi)  TmpR5(2) = Aero%RtSpd(iRtSpd)
                        if (Aero%N_VRel >0_IntKi)  TmpR5(3) = Aero%VRel (iVRel )
                        if (Aero%N_Pitch>0_IntKi)  TmpR5(4) = Aero%Pitch(iPitch)
                        if (Aero%N_Skew >0_IntKi)  TmpR5(5) = Aero%Skew (iSkew )
                        write(TmpChar,'(6(f10.3))') TmpR5(1), TmpR5(2), TmpR5(3), TmpR5(4), TmpR5(5)
                        ErrMsg3 = ErrMsg3//TmpChar//NewLine
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
   end subroutine CheckAeroTabs
end subroutine Get_RtAeroTableData


!> This subroutine counts the number of unique values in an array and returns a sorted array of them.
!! NOTE: this routine is found in the WAMIT2.f90 file as well
SUBROUTINE UniqueRealValues( DataArrayIn, DataArrayOut, NumUnique, ErrStat, ErrMsg )
   IMPLICIT NONE
   REAL(SiKi),                         INTENT(IN   )  :: DataArrayIn(:)    !< Array to search
   REAL(SiKi),       ALLOCATABLE,      INTENT(  OUT)  :: DataArrayOut(:)   !< Array to return results in
   INTEGER(IntKi),                     INTENT(  OUT)  :: NumUnique         !< Number of unique values found
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           !< Error Status
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            !< Message about the error
      ! Local variables
   REAL(SiKi)                                         :: TmpReal           !< Temporary real value
   INTEGER(IntKi)                                     :: I                 !< Generic counter
   INTEGER(IntKi)                                     :: J                 !< Generic counter
   REAL(SiKi),       ALLOCATABLE                      :: TmpRealArray(:)   !< Temporary real array
   LOGICAL                                            :: Duplicate         !< If there is a duplicate value
      ! Error handling
   INTEGER(IntKi)                                     :: ErrStatTmp
   CHARACTER(2048)                                    :: ErrMsgTmp
   CHARACTER(*), PARAMETER                            :: RoutineName = 'UniqueRealValues'

      ! Initialize things
   ErrStat     =  ErrID_None
   ErrStatTmp  =  ErrID_None
   ErrMsg      =  ''
   ErrMsgTmp   =  ''

      ! Allocate the temporary array
   CALL AllocAry( TmpRealArray, SIZE(DataArrayIn,1), 'Temporary array for data sorting', ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp
      RETURN
   ENDIF

      ! Initialize the array with a large negative number.
   TmpRealArray   = -9.9e9_SiKi

      ! The first point is unique since we haven't compared it to anything yet.
   TmpRealArray(1)   = DataArrayIn(1)
   NumUnique         =  1

      ! Step through the DataArrayIn and put unique values into TmpRealArray.  Start at second point
   DO I=2,SIZE(DataArrayIn,1)
         ! Check the current value against the largest stored value (I-1).  If the current value is
         ! larger than the last stored one, then it should be stored after it.
      IF ( DataArrayIn(I) > TmpRealArray(NumUnique) ) THEN
         TmpRealArray(NumUnique + 1) = DataArrayIn(I)
         NumUnique = NumUnique + 1
      ELSE
         ! Otherwise, if the value should not be put last, then we have to find where it goes.  Before
         ! searching for the location, first make sure this isn't a duplicate value.
         Duplicate = .FALSE.
         DO J= NumUnique, 1, -1
            IF ( EqualRealNos( DataArrayIn(I), TmpRealArray(J) )) THEN
               Duplicate = .TRUE.
               EXIT     ! Stop searching
            ENDIF
         ENDDO

         ! If this is not a duplicate, the location where it goes has to be find.  To do this, we will
         ! sequentially shift each value one index further as we step backwards through the sorted
         ! array.  When we find the location between values where this goes, we put the value there.
         IF ( .NOT. Duplicate ) THEN
            DO J= NumUnique, 1, -1           ! TempRealArray only has NumUnique values.  Everything after is junk.
               IF ( DataArrayIn(I) < TmpRealArray(J) )   THEN
                  TmpRealArray(J+1) = TmpRealArray(J)       ! Shift this value further along the array
                  IF ( J == 1 )  THEN                       ! If we are at the first value, then it goes here.
                     TmpRealArray(J) = DataArrayIn(I)
                     NumUnique = NumUnique + 1
                  ELSE
                     IF ( DataArrayIn(I) > TmpRealArray(J-1) ) THEN  ! If larger than the preceeding number, it goes here
                        TmpRealArray(J) = DataArrayIn(I)
                        NumUnique = NumUnique + 1
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF
   ENDDO

      ! Now that we have the array sorted into unique values, we need to construct an array to return the values in.
   CALL AllocAry( DataArrayOut, NumUnique, 'Return array with sorted values', ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp
      RETURN
   ENDIF

      ! Copy the values over
   DataArrayOut = TmpRealArray(1:NumUnique)
   call Cleanup()
contains
   subroutine Cleanup()
      if (allocated(TmpRealArray))  deallocate(TmpRealArray)
   end subroutine Cleanup
END SUBROUTINE UniqueRealValues


!> Check inputdata
subroutine ADskInput_ValidateInput( InitInp, InputFileData, ErrStat, ErrMsg )
   type(ADsk_InitInputType),  intent(in   )  :: InitInp              !< Input data for initialization
   type(ADsk_InputFile),      intent(in   )  :: InputFileData        !< The data for initialization
   integer(IntKi),            intent(  out)  :: ErrStat              !< Error status  from this subroutine
   character(*),              intent(  out)  :: ErrMsg               !< Error message from this subroutine
   character(*),              parameter      :: RoutineName="ADskInput_ValidateInput"

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! InitInput checks
   if (InputFileData%DT       <= 0.0_DbKi)   call SetErrStat(ErrID_Fatal,'DT must not be negative.',     ErrStat,ErrMsg,RoutineName)
   if (InputFileData%AirDens  <= 0.0_ReKi)   call SetErrStat(ErrID_Fatal,'AirDens must not be negative.',ErrStat,ErrMsg,RoutineName) 
   if (InputFileData%RotorRad <= 0.0_ReKi)   call SetErrStat(ErrID_Fatal,'RotorRad must not be negative.',ErrStat,ErrMsg,RoutineName) 
   if (InitInp%Linearize)  call SetErrStat(ErrID_Fatal,'AeroDisk cannot perform linearization analysis.',ErrStat,ErrMsg,RoutineName)

   ! Some sanity checks AeroTable
   if (InputFileData%AeroTable%N_TSR > 0_IntKi) then
      if (minval(InputFileData%AeroTable%TSR) <= 0.0_SiKi) then
         call SetErrStat(ErrID_Fatal,'All TSR values in table must be postive.',ErrStat,ErrMsg,RoutineName)
      endif
   endif
   if (InputFileData%AeroTable%N_RtSpd > 0_IntKi) then
      if (minval(InputFileData%AeroTable%RtSpd) <= 0.0_SiKi) then
         call SetErrStat(ErrID_Fatal,'All RtSpd values in table must be postive.',ErrStat,ErrMsg,RoutineName)
      endif
   endif
   if (InputFileData%AeroTable%N_VRel > 0_IntKi) then
      if (minval(InputFileData%AeroTable%VRel) < 0.0_SiKi ) then
         call SetErrStat(ErrID_Fatal,'All VRel values in table must be postive.',ErrStat,ErrMsg,RoutineName)
      endif
   endif
   if (InputFileData%AeroTable%N_Pitch > 0_IntKi) then   ! input table as deg, already converted to rad
      if (minval(InputFileData%AeroTable%Pitch) <= -Pi_S .or. maxval(InputFileData%AeroTable%Pitch) >= Pi_S) then
         call SetErrStat(ErrID_Fatal,'All Pitch values in table must be between -180 and 180 degrees.',ErrStat,ErrMsg,RoutineName)
      endif
   endif
   if (InputFileData%AeroTable%N_Skew > 0_IntKi) then    ! input table as deg, already converted to rad
      if (minval(InputFileData%AeroTable%Skew) <= 0.0_SiKi .or. maxval(InputFileData%AeroTable%Skew) >= Pi_S) then
         call SetErrStat(ErrID_Fatal,'All Skew values in table must be between 0 and 180 degrees.',ErrStat,ErrMsg,RoutineName)
      endif
   endif

   ! Either TSR or RtSpd+VRel columns must be provided
   if (InputFileData%AeroTable%N_TSR > 1_IntKi) then
      if (InputFileData%AeroTable%N_RtSpd > 1_IntKi .or. InputFileData%AeroTable%N_VRel > 1_IntKi) then
         call SetErrStat(ErrID_Fatal,'TSR values present in table along with RtSpd or VRel values. '//NewLine// &
                  ' --> Either RtSpd and VRel values must be in the table, or TSR values may be present.'//NewLine// &
                  '       To skip columns, you may enter "0" for the column index and leave the table values.',ErrStat,ErrMsg,RoutineName)
         return
      endif
   else
      if (InputFileData%AeroTable%N_RtSpd < 2_IntKi .and. InputFileData%AeroTable%N_VRel < 2_IntKi) then
         call SetErrStat(ErrID_Fatal,'TSR values NOT present in table, but RtSpd and VRel values are not both present. '//NewLine// &
                  ' --> Either RtSpd and VRel values must be in the table, or TSR values may be present.'//NewLine// &
                  '       To skip columns, you may enter "0" for the column index and leave the table values.',ErrStat,ErrMsg,RoutineName)
         return
      endif
   endif
end subroutine ADskInput_ValidateInput


!> validate and process input file data (some was done during parsing of input file)
subroutine ADskInput_SetParameters( InitInp, Interval, InputFileData, p, ErrStat, ErrMsg )
   type(ADsk_InitInputType),  intent(in   )  :: InitInp              !< Input data for initialization
   real(DbKi),                intent(inout)  :: Interval             !< Coupling interval in seconds
   type(ADsk_InputFile),      intent(inout)  :: InputFileData        !< The data for initialization
   type(ADsk_ParameterType),  intent(inout)  :: p                    !<
   integer(IntKi),            intent(  out)  :: ErrStat              !< Error status  from this subroutine
   character(*),              intent(  out)  :: ErrMsg               !< Error message from this subroutine
   integer(IntKi)                            :: ErrStat2             !< Temporary error status  for subroutine and function calls
   character(ErrMsgLen)                      :: ErrMsg2              !< Temporary error message for subroutine and function calls
   character(*),              parameter      :: RoutineName="ADskInput_SetParameters"

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Set parameters
   p%DT        = InputFileData%DT
   Interval    = p%DT                        ! Tell glue code what we want for DT
   p%numOuts   = InputFileData%NumOuts
   p%RootName  = InitInp%RootName
   p%RotorRad  = InputFileData%RotorRad
   p%AirDens   = InputFileData%AirDens
   p%UseTSR    = .false.      ! Reset below if N_TSR>1

      ! Derived parameter
   p%halfRhoA  = 0.5_ReKi * p%AirDens * Pi * p%RotorRad*p%RotorRad

      ! Table of values
   p%AeroTable%N_TSR   = InputFileData%AeroTable%N_TSR
   p%AeroTable%N_RtSpd = InputFileData%AeroTable%N_RtSpd
   p%AeroTable%N_VRel  = InputFileData%AeroTable%N_VRel
   p%AeroTable%N_Pitch = InputFileData%AeroTable%N_Pitch
   p%AeroTable%N_Skew  = InputFileData%AeroTable%N_Skew
   if (allocated( InputFileData%AeroTable%TSR  ))  call move_alloc( InputFileData%AeroTable%TSR,   p%AeroTable%TSR   )
   if (allocated( InputFileData%AeroTable%RtSpd))  call move_alloc( InputFileData%AeroTable%RtSpd, p%AeroTable%RtSpd )
   if (allocated( InputFileData%AeroTable%VRel ))  call move_alloc( InputFileData%AeroTable%VRel,  p%AeroTable%VRel  )
   if (allocated( InputFileData%AeroTable%Pitch))  call move_alloc( InputFileData%AeroTable%Pitch, p%AeroTable%Pitch )
   if (allocated( InputFileData%AeroTable%Skew ))  call move_alloc( InputFileData%AeroTable%Skew,  p%AeroTable%Skew  )
   if (allocated( InputFileData%AeroTable%C_Fx ))  call move_alloc( InputFileData%AeroTable%C_Fx,  p%AeroTable%C_Fx  )
   if (allocated( InputFileData%AeroTable%C_Fy ))  call move_alloc( InputFileData%AeroTable%C_Fy,  p%AeroTable%C_Fy  )
   if (allocated( InputFileData%AeroTable%C_Fz ))  call move_alloc( InputFileData%AeroTable%C_Fz,  p%AeroTable%C_Fz  )
   if (allocated( InputFileData%AeroTable%C_Mx ))  call move_alloc( InputFileData%AeroTable%C_Mx,  p%AeroTable%C_Mx  )
   if (allocated( InputFileData%AeroTable%C_My ))  call move_alloc( InputFileData%AeroTable%C_My,  p%AeroTable%C_My  )
   if (allocated( InputFileData%AeroTable%C_Mz ))  call move_alloc( InputFileData%AeroTable%C_Mz,  p%AeroTable%C_Mz  )

      ! Use the TSR values
   if (p%AeroTable%N_TSR > 1_IntKi)    p%UseTSR = .true.

      ! Set the outputs
   call SetOutParam(InputFileData%OutList, p, ErrStat, ErrMsg )
end subroutine ADskInput_SetParameters


!> Write the table out to a whatever UnOut is (as long as > 0).
subroutine WriteAeroTab(Aero, UnOut)
   type(ADsk_AeroTable),   intent(in   )  :: Aero
   integer(IntKi),         intent(in   )  :: UnOut
   integer(IntKi)                         :: i1,i2,i3,i4,i5       !< loop counters
   character(*),           parameter      :: RoutineName="ADskInput_SetParameters"
   if (UnOut <= 0_IntKi) return
   ! Write header info
   write(UnOut,'(A)')   '======================================='
   write(UnOut,'(A)')   'AeroDisk Actuator Disk Properties table'
   write(UnOut,'(A)')   '  NOTE: the units correspond to units used internally within code, not units of input'
   if (Aero%N_TSR   > 0_IntKi) then
      write(UnOut,'(A12,I5,A7)') '     TSR    ',Aero%N_TSR,' values'
   else
      write(UnOut,'(A)')   '     TSR    ---- unused ----'
   endif
   if (Aero%N_RtSpd > 0_IntKi) then
      write(UnOut,'(A12,I5,A7)') '     RtSpd  ',Aero%N_RtSpd,' values'
   else
      write(UnOut,'(A)')   '     RtSpd  ---- unused ----'
   endif
   if (Aero%N_VRel  > 0_IntKi) then
      write(UnOut,'(A12,I5,A7)') '     VRel   ',Aero%N_VRel,' values'
   else
      write(UnOut,'(A)')   '     VRel   ---- unused ----'
   endif
   if (Aero%N_Pitch > 0_IntKi) then
      write(UnOut,'(A12,I5,A7)') '     Pitch  ',Aero%N_Pitch,' values'
   else
      write(UnOut,'(A)')   '     Pitch  ---- unused ----'
   endif
   if (Aero%N_Skew  > 0_IntKi) then
      write(UnOut,'(A12,I5,A7)') '     Skew   ',Aero%N_Skew,' values'
   else
      write(UnOut,'(A)')   '     Skew   ---- unused ----'
   endif
   ! Table header
   if (Aero%N_TSR   > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      TSR      '
   if (Aero%N_Skew  > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      RtSpd    '
   if (Aero%N_VRel  > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      VRel     '
   if (Aero%N_Pitch > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      Pitch    '
   if (Aero%N_Skew  > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      Skew     '
   write(UnOut,'(A15)',ADVANCE='NO')   '      C_fx     '
   write(UnOut,'(A15)',ADVANCE='NO')   '      C_fy     '
   write(UnOut,'(A15)',ADVANCE='NO')   '      C_fz     '
   write(UnOut,'(A15)',ADVANCE='NO')   '      C_mx     '
   write(UnOut,'(A15)',ADVANCE='NO')   '      C_my     '
   write(UnOut,'(A15)',ADVANCE='NO')   '      C_mz     '
   write(UnOut,'(A)') ''
   if (Aero%N_TSR   > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      (-)      '
   if (Aero%N_Skew  > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      (rad/s)  '
   if (Aero%N_VRel  > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      (m/s)    '
   if (Aero%N_Pitch > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      (rad)    '
   if (Aero%N_Skew  > 0_IntKi)   write(UnOut,'(A15)',ADVANCE='NO')   '      (rad)    '
   write(UnOut,'(A15)',ADVANCE='NO')   '      (-)      '
   write(UnOut,'(A15)',ADVANCE='NO')   '      (-)      '
   write(UnOut,'(A15)',ADVANCE='NO')   '      (-)      '
   write(UnOut,'(A15)',ADVANCE='NO')   '      (-)      '
   write(UnOut,'(A15)',ADVANCE='NO')   '      (-)      '
   write(UnOut,'(A15)',ADVANCE='NO')   '      (-)      '
   write(UnOut,'(A)') ''
   ! Table itself
   do i1=1,max(1,Aero%N_TSR  )
      do i2=1,max(1,Aero%N_Skew )
         do i3=1,max(1,Aero%N_VRel )
            do i4=1,max(1,Aero%N_Pitch)
               do i5=1,max(1,Aero%N_Skew )
                  if (Aero%N_TSR  > 0_IntKi) write(UnOut,'(2x,f13.6)',ADVANCE='NO')  Aero%TSR  (i1)
                  if (Aero%N_Skew > 0_IntKi) write(UnOut,'(2x,f13.6)',ADVANCE='NO')  Aero%Skew (i2)
                  if (Aero%N_VRel > 0_IntKi) write(UnOut,'(2x,f13.6)',ADVANCE='NO')  Aero%VRel (i3)
                  if (Aero%N_Pitch> 0_IntKi) write(UnOut,'(2x,f13.6)',ADVANCE='NO')  Aero%Pitch(i4)
                  if (Aero%N_Skew > 0_IntKi) write(UnOut,'(2x,f13.6)',ADVANCE='NO')  Aero%Skew (i5)
                  write(UnOut,   '(f13.6)',ADVANCE='NO')    Aero%C_Fx(i1,i2,i3,i4,i5)
                  write(UnOut,'(2x,f13.6)',ADVANCE='NO')    Aero%C_Fy(i1,i2,i3,i4,i5)
                  write(UnOut,'(2x,f13.6)',ADVANCE='NO')    Aero%C_Fz(i1,i2,i3,i4,i5)
                  write(UnOut,'(2x,f13.6)',ADVANCE='NO')    Aero%C_Mx(i1,i2,i3,i4,i5)
                  write(UnOut,'(2x,f13.6)',ADVANCE='NO')    Aero%C_My(i1,i2,i3,i4,i5)
                  write(UnOut,'(2x,f13.6)',ADVANCE='NO')    Aero%C_Mz(i1,i2,i3,i4,i5)
                  write(UnOut,'(A)') ''
               enddo
            enddo
         enddo
      enddo
   enddo
end subroutine WriteAeroTab


!----------------------------------------------------------------------------------------------------------------------------------
!> this routine fills the AllOuts array, which is used to send data to the glue code to be written to an output file.
!! NOTE: AllOuts is ReKi, but most calculations in this module are in single precision. This requires a bunch of conversions at this
!! stage.
subroutine Calc_WriteOutput( u, p, y, m, ErrStat, ErrMsg, CalcWriteOutput )
   type(ADsk_InputType),         intent(in   )  :: u                 !< The inputs at time T
   type(ADsk_ParameterType),     intent(in   )  :: p                 !< The module parameters
   type(ADsk_OutputType),        intent(in   )  :: y                 !< outputs
   type(ADsk_MiscVarType),       intent(inout)  :: m                 !< misc/optimization variables (for computing mesh transfers)
   integer(IntKi),               intent(  out)  :: ErrStat           !< The error status code
   character(*),                 intent(  out)  :: ErrMsg            !< The error message, if an error occurred
   logical,                      intent(in   )  :: CalcWriteOutput   !< flag that determines if we need to compute AllOuts (or just the reaction loads that get returned to ServoDyn)
   ! local variables
   character(*), parameter                      :: RoutineName = 'Calc_WriteOutput'
   integer(IntKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   real(ReKi)                                   :: Tmp3(3)
   real(ReKi)                                   :: Rxyz(3,3)         !< rotation matrix for x,y,z of local coordinates

   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! return if we are not providing outputs
   if (.not. CalcWriteOutput) return

   ! rotation matrix
!FIXME: make sure this is actually correct and not the transpose of what we want
   Rxyz(1:3,1) = real(m%x_hat, ReKi)
   Rxyz(1:3,2) = real(m%y_hat, ReKi)
   Rxyz(1:3,3) = real(m%z_hat, ReKi)

   ! Rotspeed etc
   m%AllOuts( ADSpeed   ) = u%RotSpeed * 30.0_ReKi / Pi
   m%AllOuts( ADPitch   ) = u%BlPitch * 180.0_ReKi / Pi

   m%AllOuts( ADTSR     ) = real(m%lambda, ReKi)      ! TSR -- tip speed ratio   (-)
   m%AllOuts( ADVRel    ) = real(m%VRel, ReKi)        ! magnitude of VRel vector (m/s)
   m%AllOuts( ADSkew    ) = real(m%Chi * 180.0_ReKi / Pi, ReKi)

   ! Wind in local frame, inertial frame
   Tmp3 = matmul(Rxyz(1:3,1:3), m%DiskAvgVel)
   m%AllOuts( ADVWindx  ) = Tmp3(1)
   m%AllOuts( ADVWindy  ) = Tmp3(2)
   m%AllOuts( ADVWindz  ) = Tmp3(3)
   m%AllOuts( ADVWindxi ) = m%DiskAvgVel(1)
   m%AllOuts( ADVWindyi ) = m%DiskAvgVel(2)
   m%AllOuts( ADVWindzi ) = m%DiskAvgVel(3)

   ! Rotor velocity in local frame, inertial frame
   Tmp3 = matmul(Rxyz(1:3,1:3), u%HubMotion%TranslationVel(1:3,1))
   m%AllOuts( ADSTVx    ) = Tmp3(1)
   m%AllOuts( ADSTVy    ) = Tmp3(2)
   m%AllOuts( ADSTVz    ) = Tmp3(3)
   m%AllOuts( ADSTVxi   ) = u%HubMotion%TranslationVel(1,1)
   m%AllOuts( ADSTVyi   ) = u%HubMotion%TranslationVel(2,1)
   m%AllOuts( ADSTVzi   ) = u%HubMotion%TranslationVel(3,1)

   ! Coefficients
   if (EqualRealNos(m%VRel_xd,0.0_SiKi)) then
      m%AllOuts( ADCp   ) = 0.0_ReKi
   else
      m%AllOuts( ADCp   ) = (real(m%Moment(1),ReKi) * u%RotSpeed) / (p%halfRhoA * real(m%Vrel_xd,ReKi)**3_IntKi )
   endif
   m%AllOuts( ADCt      ) = real(m%C_F(1),ReKi)
   m%AllOuts( ADCq      ) = real(m%C_M(1),ReKi)

   ! Power
   m%AllOuts( ADPower   ) = real(m%Moment(1),ReKi) * u%RotSpeed

   ! Resulting forces
   m%AllOuts( ADFx      ) = real(m%Force(1), ReKi)
   m%AllOuts( ADFy      ) = real(m%Force(2), ReKi)
   m%AllOuts( ADFz      ) = real(m%Force(3), ReKi)
   m%AllOuts( ADMx      ) = real(m%Moment(1),ReKi)
   m%AllOuts( ADMy      ) = real(m%Moment(2),ReKi)
   m%AllOuts( ADMz      ) = real(m%Moment(3),ReKi)
   !Tmp3 = m%Force( 1)*m%x_hat + m%Force( 2)*m%y_hat + m%Force( 3)*m%z_hat
   Tmp3 = matmul(real(m%Force(1:3),ReKi), Rxyz(1:3,1:3))
   m%AllOuts( ADFxi     ) = Tmp3(1)
   m%AllOuts( ADFyi     ) = Tmp3(2)
   m%AllOuts( ADFzi     ) = Tmp3(3)
   !Tmp3 = m%Moment(1)*m%x_hat + m%Moment(2)*m%y_hat + m%Moment(3)*m%z_hat
   Tmp3 = matmul(real(m%Force(1:3),ReKi), Rxyz(1:3,1:3))
   m%AllOuts( ADMxi     ) = Tmp3(1)
   m%AllOuts( ADMyi     ) = Tmp3(2)
   m%AllOuts( ADMzi     ) = Tmp3(3)

end subroutine Calc_WriteOutput


!**********************************************************************************************************************************
! NOTE: The following lines of code were generated by a Matlab script called "Write_ChckOutLst.m"
!      using the parameters listed in the "OutListParameters.xlsx" Excel file. Any changes to these 
!      lines should be modified in the Matlab script and/or Excel worksheet as necessary. 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine checks to see if any requested output channel names (stored in the OutList(:)) are invalid. It returns a 
!! warning if any of the channels are not available outputs from the module.
!!  It assigns the settings for OutParam(:) (i.e, the index, name, and units of the output channels, WriteOutput(:)).
!!  the sign is set to 0 if the channel is invalid.
!! It sets assumes the value p%NumOuts has been set before this routine has been called, and it sets the values of p%OutParam here.
!! 
!! This routine was generated by Write_ChckOutLst.m using the parameters listed in OutListParameters.xlsx at 24-Feb-2022 16:52:56.
SUBROUTINE SetOutParam(OutList, p, ErrStat, ErrMsg )
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables

   CHARACTER(ChanLen),        INTENT(IN)     :: OutList(:)                        !< The list out user-requested outputs
   TYPE(ADsk_ParameterType),  INTENT(INOUT)  :: p                                 !< The module parameters
   INTEGER(IntKi),            INTENT(OUT)    :: ErrStat                           !< The error status code
   CHARACTER(*),              INTENT(OUT)    :: ErrMsg                            !< The error message, if an error occurred

      ! Local variables

   INTEGER                      :: ErrStat2                                        ! temporary (local) error status
   INTEGER                      :: I                                               ! Generic loop-counting index
   INTEGER                      :: J                                               ! Generic loop-counting index
   INTEGER                      :: INDX                                            ! Index for valid arrays

   LOGICAL                      :: CheckOutListAgain                               ! Flag used to determine if output parameter starting with "M" is valid (or the negative of another parameter)
   LOGICAL                      :: InvalidOutput(0:MaxOutPts)                      ! This array determines if the output channel is valid for this configuration
   CHARACTER(ChanLen)           :: OutListTmp                                      ! A string to temporarily hold OutList(I)
   CHARACTER(*), PARAMETER      :: RoutineName = "SetOutParam"

   CHARACTER(OutStrLenM1), PARAMETER  :: ValidParamAry(34) =  (/  &   ! This lists the names of the allowed parameters, which must be sorted alphabetically
                               "ADCP     ","ADCQ     ","ADCT     ","ADFX     ","ADFXI    ","ADFY     ","ADFYI    ","ADFZ     ", &
                               "ADFZI    ","ADMX     ","ADMXI    ","ADMY     ","ADMYI    ","ADMZ     ","ADMZI    ","ADPITCH  ", &
                               "ADPOWER  ","ADSKEW   ","ADSPEED  ","ADSTVX   ","ADSTVXI  ","ADSTVY   ","ADSTVYI  ","ADSTVZ   ", &
                               "ADSTVZI  ","ADTSR    ","ADVREL   ","ADVWINDX ","ADVWINDXI","ADVWINDY ","ADVWINDYI","ADVWINDZ ", &
                               "ADVWINDZI","ADYAWERR "/)
   INTEGER(IntKi), PARAMETER :: ParamIndxAry(34) =  (/ &                            ! This lists the index into AllOuts(:) of the allowed parameters ValidParamAry(:)
                                     ADCp ,      ADCq ,      ADCt ,      ADFx ,     ADFxi ,      ADFy ,     ADFyi ,      ADFz , &
                                    ADFzi ,      ADMx ,     ADMxi ,      ADMy ,     ADMyi ,      ADMz ,     ADMzi ,   ADPitch , &
                                  ADPower ,    ADSkew ,   ADSpeed ,    ADSTVx ,   ADSTVxi ,    ADSTVy ,   ADSTVyi ,    ADSTVz , &
                                  ADSTVzi ,     ADTSR ,    ADVRel ,  ADVWindx , ADVWindxi ,  ADVWindy , ADVWindyi ,  ADVWindz , &
                                ADVWindzi ,  ADYawErr /)
   CHARACTER(ChanLen), PARAMETER :: ParamUnitsAry(34) =  (/  &  ! This lists the units corresponding to the allowed parameters
                               "(-)  ","(-)  ","(-)  ","(N)  ","(N)  ","(N)  ","(N)  ","(N)  ", &
                               "(N)  ","(N-m)","(N-m)","(N-m)","(N-m)","(N-m)","(N-m)","(deg)", &
                               "(W)  ","(deg)","(rpm)","(m/s)","(m/s)","(m/s)","(m/s)","(m/s)", &
                               "(m/s)","(-)  ","(m/s)","(m/s)","(m/s)","(m/s)","(m/s)","(m/s)", &
                               "(m/s)","(deg)"/)


      ! Initialize values
   ErrStat = ErrID_None
   ErrMsg = ""
   InvalidOutput = .FALSE.


!   ..... Developer must add checking for invalid inputs here: .....

!   ................. End of validity checking .................


   !-------------------------------------------------------------------------------------------------
   ! Allocate and set index, name, and units for the output channels
   ! If a selected output channel is not available in this module, set error flag.
   !-------------------------------------------------------------------------------------------------

   ALLOCATE ( p%OutParam(0:p%NumOuts) , STAT=ErrStat2 )
   IF ( ErrStat2 /= 0_IntKi )  THEN
      CALL SetErrStat( ErrID_Fatal,"Error allocating memory for the AeroDisk OutParam array.", ErrStat, ErrMsg, RoutineName )
      RETURN
   ENDIF

      ! Set index, name, and units for the time output channel:

   p%OutParam(0)%Indx  = Time
   p%OutParam(0)%Name  = "Time"    ! OutParam(0) is the time channel by default.
   p%OutParam(0)%Units = "(s)"
   p%OutParam(0)%SignM = 1


      ! Set index, name, and units for all of the output channels.
      ! If a selected output channel is not available by this module set ErrStat = ErrID_Warn.

   DO I = 1,p%NumOuts

      p%OutParam(I)%Name  = OutList(I)
      OutListTmp          = OutList(I)

      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !   channel name with a "-", "_", "m", or "M" character indicating "minus".


      CheckOutListAgain = .FALSE.

      IF      ( INDEX( "-_", OutListTmp(1:1) ) > 0 ) THEN
         p%OutParam(I)%SignM = -1                         ! ex, "-TipDxc1" causes the sign of TipDxc1 to be switched.
         OutListTmp          = OutListTmp(2:)
      ELSE IF ( INDEX( "mM", OutListTmp(1:1) ) > 0 ) THEN ! We'll assume this is a variable name for now, (if not, we will check later if OutListTmp(2:) is also a variable name)
         CheckOutListAgain   = .TRUE.
         p%OutParam(I)%SignM = 1
      ELSE
         p%OutParam(I)%SignM = 1
      END IF

      CALL Conv2UC( OutListTmp )    ! Convert OutListTmp to upper case


      Indx = IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )


         ! If it started with an "M" (CheckOutListAgain) we didn't find the value in our list (Indx < 1)

      IF ( CheckOutListAgain .AND. Indx < 1 ) THEN    ! Let's assume that "M" really meant "minus" and then test again
         p%OutParam(I)%SignM = -1                     ! ex, "MTipDxc1" causes the sign of TipDxc1 to be switched.
         OutListTmp          = OutListTmp(2:)

         Indx = IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )
      END IF


      IF ( Indx > 0 ) THEN ! we found the channel name
         IF ( InvalidOutput( ParamIndxAry(Indx) ) ) THEN  ! but, it isn't valid for these settings
            p%OutParam(I)%Indx  = 0                 ! pick any valid channel (I just picked "Time=0" here because it's universal)
            p%OutParam(I)%Units = "INVALID"
            p%OutParam(I)%SignM = 0
         ELSE
            p%OutParam(I)%Indx  = ParamIndxAry(Indx)
            p%OutParam(I)%Units = ParamUnitsAry(Indx) ! it's a valid output
         END IF
      ELSE ! this channel isn't valid
         p%OutParam(I)%Indx  = 0                    ! pick any valid channel (I just picked "Time=0" here because it's universal)
         p%OutParam(I)%Units = "INVALID"
         p%OutParam(I)%SignM = 0                    ! multiply all results by zero

         CALL SetErrStat(ErrID_Fatal, TRIM(p%OutParam(I)%Name)//" is not an available output channel.",ErrStat,ErrMsg,RoutineName)
      END IF

   END DO

   RETURN
END SUBROUTINE SetOutParam
!----------------------------------------------------------------------------------------------------------------------------------
!End of code generated by Matlab script
!**********************************************************************************************************************************
END MODULE AeroDisk_IO
