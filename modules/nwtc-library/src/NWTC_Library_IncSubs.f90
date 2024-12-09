!STARTOFREGISTRYGENERATEDFILE 'NWTC_Library_Subs.f90'
!
! WARNING This file is generated automatically by the FAST registry.
! Do not edit.  Your changes to this file will be lost.
!
! FAST Registry'

subroutine NWTC_Library_CopyMapType(SrcMapTypeData, DstMapTypeData, CtrlCode, ErrStat, ErrMsg)
   type(MapType), intent(in) :: SrcMapTypeData
   type(MapType), intent(inout) :: DstMapTypeData
   integer(IntKi),  intent(in   ) :: CtrlCode
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   character(*), parameter        :: RoutineName = 'NWTC_Library_CopyMapType'
   ErrStat = ErrID_None
   ErrMsg  = ''
   DstMapTypeData%OtherMesh_Element = SrcMapTypeData%OtherMesh_Element
   DstMapTypeData%distance = SrcMapTypeData%distance
   DstMapTypeData%couple_arm = SrcMapTypeData%couple_arm
   DstMapTypeData%shape_fn = SrcMapTypeData%shape_fn
end subroutine

subroutine NWTC_Library_DestroyMapType(MapTypeData, ErrStat, ErrMsg)
   type(MapType), intent(inout) :: MapTypeData
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   character(*), parameter        :: RoutineName = 'NWTC_Library_DestroyMapType'
   ErrStat = ErrID_None
   ErrMsg  = ''
end subroutine

subroutine NWTC_Library_PackMapType(RF, Indata)
   type(RegFile), intent(inout) :: RF
   type(MapType), intent(in) :: InData
   character(*), parameter         :: RoutineName = 'NWTC_Library_PackMapType'
   if (RF%ErrStat >= AbortErrLev) return
   call RegPack(RF, InData%OtherMesh_Element)
   call RegPack(RF, InData%distance)
   call RegPack(RF, InData%couple_arm)
   call RegPack(RF, InData%shape_fn)
   if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine NWTC_Library_UnPackMapType(RF, OutData)
   type(RegFile), intent(inout)    :: RF
   type(MapType), intent(inout) :: OutData
   character(*), parameter            :: RoutineName = 'NWTC_Library_UnPackMapType'
   if (RF%ErrStat /= ErrID_None) return
   call RegUnpack(RF, OutData%OtherMesh_Element); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpack(RF, OutData%distance); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpack(RF, OutData%couple_arm); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpack(RF, OutData%shape_fn); if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine NWTC_Library_CopyMeshMapLinearizationType(SrcMeshMapLinearizationTypeData, DstMeshMapLinearizationTypeData, CtrlCode, ErrStat, ErrMsg)
   type(MeshMapLinearizationType), intent(in) :: SrcMeshMapLinearizationTypeData
   type(MeshMapLinearizationType), intent(inout) :: DstMeshMapLinearizationTypeData
   integer(IntKi),  intent(in   ) :: CtrlCode
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   integer(B4Ki)                  :: LB(2), UB(2)
   integer(IntKi)                 :: ErrStat2
   character(*), parameter        :: RoutineName = 'NWTC_Library_CopyMeshMapLinearizationType'
   ErrStat = ErrID_None
   ErrMsg  = ''
   if (allocated(SrcMeshMapLinearizationTypeData%mi)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%mi)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%mi)
      if (.not. allocated(DstMeshMapLinearizationTypeData%mi)) then
         allocate(DstMeshMapLinearizationTypeData%mi(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%mi.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%mi = SrcMeshMapLinearizationTypeData%mi
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%fx_p)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%fx_p)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%fx_p)
      if (.not. allocated(DstMeshMapLinearizationTypeData%fx_p)) then
         allocate(DstMeshMapLinearizationTypeData%fx_p(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%fx_p.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%fx_p = SrcMeshMapLinearizationTypeData%fx_p
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%tv_uD)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%tv_uD)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%tv_uD)
      if (.not. allocated(DstMeshMapLinearizationTypeData%tv_uD)) then
         allocate(DstMeshMapLinearizationTypeData%tv_uD(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%tv_uD.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%tv_uD = SrcMeshMapLinearizationTypeData%tv_uD
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%tv_uS)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%tv_uS)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%tv_uS)
      if (.not. allocated(DstMeshMapLinearizationTypeData%tv_uS)) then
         allocate(DstMeshMapLinearizationTypeData%tv_uS(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%tv_uS.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%tv_uS = SrcMeshMapLinearizationTypeData%tv_uS
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%ta_uD)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%ta_uD)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%ta_uD)
      if (.not. allocated(DstMeshMapLinearizationTypeData%ta_uD)) then
         allocate(DstMeshMapLinearizationTypeData%ta_uD(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%ta_uD.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%ta_uD = SrcMeshMapLinearizationTypeData%ta_uD
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%ta_uS)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%ta_uS)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%ta_uS)
      if (.not. allocated(DstMeshMapLinearizationTypeData%ta_uS)) then
         allocate(DstMeshMapLinearizationTypeData%ta_uS(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%ta_uS.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%ta_uS = SrcMeshMapLinearizationTypeData%ta_uS
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%ta_rv)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%ta_rv)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%ta_rv)
      if (.not. allocated(DstMeshMapLinearizationTypeData%ta_rv)) then
         allocate(DstMeshMapLinearizationTypeData%ta_rv(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%ta_rv.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%ta_rv = SrcMeshMapLinearizationTypeData%ta_rv
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%li)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%li)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%li)
      if (.not. allocated(DstMeshMapLinearizationTypeData%li)) then
         allocate(DstMeshMapLinearizationTypeData%li(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%li.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%li = SrcMeshMapLinearizationTypeData%li
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%M_uS)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%M_uS)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%M_uS)
      if (.not. allocated(DstMeshMapLinearizationTypeData%M_uS)) then
         allocate(DstMeshMapLinearizationTypeData%M_uS(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%M_uS.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%M_uS = SrcMeshMapLinearizationTypeData%M_uS
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%M_uD)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%M_uD)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%M_uD)
      if (.not. allocated(DstMeshMapLinearizationTypeData%M_uD)) then
         allocate(DstMeshMapLinearizationTypeData%M_uD(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%M_uD.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%M_uD = SrcMeshMapLinearizationTypeData%M_uD
   end if
   if (allocated(SrcMeshMapLinearizationTypeData%M_f)) then
      LB(1:2) = lbound(SrcMeshMapLinearizationTypeData%M_f)
      UB(1:2) = ubound(SrcMeshMapLinearizationTypeData%M_f)
      if (.not. allocated(DstMeshMapLinearizationTypeData%M_f)) then
         allocate(DstMeshMapLinearizationTypeData%M_f(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapLinearizationTypeData%M_f.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapLinearizationTypeData%M_f = SrcMeshMapLinearizationTypeData%M_f
   end if
end subroutine

subroutine NWTC_Library_DestroyMeshMapLinearizationType(MeshMapLinearizationTypeData, ErrStat, ErrMsg)
   type(MeshMapLinearizationType), intent(inout) :: MeshMapLinearizationTypeData
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   character(*), parameter        :: RoutineName = 'NWTC_Library_DestroyMeshMapLinearizationType'
   ErrStat = ErrID_None
   ErrMsg  = ''
   if (allocated(MeshMapLinearizationTypeData%mi)) then
      deallocate(MeshMapLinearizationTypeData%mi)
   end if
   if (allocated(MeshMapLinearizationTypeData%fx_p)) then
      deallocate(MeshMapLinearizationTypeData%fx_p)
   end if
   if (allocated(MeshMapLinearizationTypeData%tv_uD)) then
      deallocate(MeshMapLinearizationTypeData%tv_uD)
   end if
   if (allocated(MeshMapLinearizationTypeData%tv_uS)) then
      deallocate(MeshMapLinearizationTypeData%tv_uS)
   end if
   if (allocated(MeshMapLinearizationTypeData%ta_uD)) then
      deallocate(MeshMapLinearizationTypeData%ta_uD)
   end if
   if (allocated(MeshMapLinearizationTypeData%ta_uS)) then
      deallocate(MeshMapLinearizationTypeData%ta_uS)
   end if
   if (allocated(MeshMapLinearizationTypeData%ta_rv)) then
      deallocate(MeshMapLinearizationTypeData%ta_rv)
   end if
   if (allocated(MeshMapLinearizationTypeData%li)) then
      deallocate(MeshMapLinearizationTypeData%li)
   end if
   if (allocated(MeshMapLinearizationTypeData%M_uS)) then
      deallocate(MeshMapLinearizationTypeData%M_uS)
   end if
   if (allocated(MeshMapLinearizationTypeData%M_uD)) then
      deallocate(MeshMapLinearizationTypeData%M_uD)
   end if
   if (allocated(MeshMapLinearizationTypeData%M_f)) then
      deallocate(MeshMapLinearizationTypeData%M_f)
   end if
end subroutine

subroutine NWTC_Library_PackMeshMapLinearizationType(RF, Indata)
   type(RegFile), intent(inout) :: RF
   type(MeshMapLinearizationType), intent(in) :: InData
   character(*), parameter         :: RoutineName = 'NWTC_Library_PackMeshMapLinearizationType'
   if (RF%ErrStat >= AbortErrLev) return
   call RegPackAlloc(RF, InData%mi)
   call RegPackAlloc(RF, InData%fx_p)
   call RegPackAlloc(RF, InData%tv_uD)
   call RegPackAlloc(RF, InData%tv_uS)
   call RegPackAlloc(RF, InData%ta_uD)
   call RegPackAlloc(RF, InData%ta_uS)
   call RegPackAlloc(RF, InData%ta_rv)
   call RegPackAlloc(RF, InData%li)
   call RegPackAlloc(RF, InData%M_uS)
   call RegPackAlloc(RF, InData%M_uD)
   call RegPackAlloc(RF, InData%M_f)
   if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine NWTC_Library_UnPackMeshMapLinearizationType(RF, OutData)
   type(RegFile), intent(inout)    :: RF
   type(MeshMapLinearizationType), intent(inout) :: OutData
   character(*), parameter            :: RoutineName = 'NWTC_Library_UnPackMeshMapLinearizationType'
   integer(B4Ki)   :: LB(2), UB(2)
   integer(IntKi)  :: stat
   logical         :: IsAllocAssoc
   if (RF%ErrStat /= ErrID_None) return
   call RegUnpackAlloc(RF, OutData%mi); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%fx_p); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%tv_uD); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%tv_uS); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%ta_uD); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%ta_uS); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%ta_rv); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%li); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%M_uS); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%M_uD); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%M_f); if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine NWTC_Library_CopyMeshMapType(SrcMeshMapTypeData, DstMeshMapTypeData, CtrlCode, ErrStat, ErrMsg)
   type(MeshMapType), intent(inout) :: SrcMeshMapTypeData
   type(MeshMapType), intent(inout) :: DstMeshMapTypeData
   integer(IntKi),  intent(in   ) :: CtrlCode
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   integer(B4Ki)   :: i1, i2, i3
   integer(B4Ki)                  :: LB(3), UB(3)
   integer(IntKi)                 :: ErrStat2
   character(ErrMsgLen)           :: ErrMsg2
   character(*), parameter        :: RoutineName = 'NWTC_Library_CopyMeshMapType'
   ErrStat = ErrID_None
   ErrMsg  = ''
   if (allocated(SrcMeshMapTypeData%MapLoads)) then
      LB(1:1) = lbound(SrcMeshMapTypeData%MapLoads)
      UB(1:1) = ubound(SrcMeshMapTypeData%MapLoads)
      if (.not. allocated(DstMeshMapTypeData%MapLoads)) then
         allocate(DstMeshMapTypeData%MapLoads(LB(1):UB(1)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%MapLoads.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      do i1 = LB(1), UB(1)
         call NWTC_Library_CopyMapType(SrcMeshMapTypeData%MapLoads(i1), DstMeshMapTypeData%MapLoads(i1), CtrlCode, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
      end do
   end if
   if (allocated(SrcMeshMapTypeData%MapMotions)) then
      LB(1:1) = lbound(SrcMeshMapTypeData%MapMotions)
      UB(1:1) = ubound(SrcMeshMapTypeData%MapMotions)
      if (.not. allocated(DstMeshMapTypeData%MapMotions)) then
         allocate(DstMeshMapTypeData%MapMotions(LB(1):UB(1)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%MapMotions.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      do i1 = LB(1), UB(1)
         call NWTC_Library_CopyMapType(SrcMeshMapTypeData%MapMotions(i1), DstMeshMapTypeData%MapMotions(i1), CtrlCode, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
      end do
   end if
   if (allocated(SrcMeshMapTypeData%MapSrcToAugmt)) then
      LB(1:1) = lbound(SrcMeshMapTypeData%MapSrcToAugmt)
      UB(1:1) = ubound(SrcMeshMapTypeData%MapSrcToAugmt)
      if (.not. allocated(DstMeshMapTypeData%MapSrcToAugmt)) then
         allocate(DstMeshMapTypeData%MapSrcToAugmt(LB(1):UB(1)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%MapSrcToAugmt.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      do i1 = LB(1), UB(1)
         call NWTC_Library_CopyMapType(SrcMeshMapTypeData%MapSrcToAugmt(i1), DstMeshMapTypeData%MapSrcToAugmt(i1), CtrlCode, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) return
      end do
   end if
   call MeshCopy(SrcMeshMapTypeData%Augmented_Ln2_Src, DstMeshMapTypeData%Augmented_Ln2_Src, CtrlCode, ErrStat2, ErrMsg2 )
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
   call MeshCopy(SrcMeshMapTypeData%Lumped_Points_Src, DstMeshMapTypeData%Lumped_Points_Src, CtrlCode, ErrStat2, ErrMsg2 )
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
   if (allocated(SrcMeshMapTypeData%LoadLn2_A_Mat_Piv)) then
      LB(1:1) = lbound(SrcMeshMapTypeData%LoadLn2_A_Mat_Piv)
      UB(1:1) = ubound(SrcMeshMapTypeData%LoadLn2_A_Mat_Piv)
      if (.not. allocated(DstMeshMapTypeData%LoadLn2_A_Mat_Piv)) then
         allocate(DstMeshMapTypeData%LoadLn2_A_Mat_Piv(LB(1):UB(1)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%LoadLn2_A_Mat_Piv.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapTypeData%LoadLn2_A_Mat_Piv = SrcMeshMapTypeData%LoadLn2_A_Mat_Piv
   end if
   if (allocated(SrcMeshMapTypeData%DisplacedPosition)) then
      LB(1:3) = lbound(SrcMeshMapTypeData%DisplacedPosition)
      UB(1:3) = ubound(SrcMeshMapTypeData%DisplacedPosition)
      if (.not. allocated(DstMeshMapTypeData%DisplacedPosition)) then
         allocate(DstMeshMapTypeData%DisplacedPosition(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%DisplacedPosition.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapTypeData%DisplacedPosition = SrcMeshMapTypeData%DisplacedPosition
   end if
   if (allocated(SrcMeshMapTypeData%LoadLn2_A_Mat)) then
      LB(1:2) = lbound(SrcMeshMapTypeData%LoadLn2_A_Mat)
      UB(1:2) = ubound(SrcMeshMapTypeData%LoadLn2_A_Mat)
      if (.not. allocated(DstMeshMapTypeData%LoadLn2_A_Mat)) then
         allocate(DstMeshMapTypeData%LoadLn2_A_Mat(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%LoadLn2_A_Mat.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapTypeData%LoadLn2_A_Mat = SrcMeshMapTypeData%LoadLn2_A_Mat
   end if
   if (allocated(SrcMeshMapTypeData%LoadLn2_F)) then
      LB(1:2) = lbound(SrcMeshMapTypeData%LoadLn2_F)
      UB(1:2) = ubound(SrcMeshMapTypeData%LoadLn2_F)
      if (.not. allocated(DstMeshMapTypeData%LoadLn2_F)) then
         allocate(DstMeshMapTypeData%LoadLn2_F(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%LoadLn2_F.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapTypeData%LoadLn2_F = SrcMeshMapTypeData%LoadLn2_F
   end if
   if (allocated(SrcMeshMapTypeData%LoadLn2_M)) then
      LB(1:2) = lbound(SrcMeshMapTypeData%LoadLn2_M)
      UB(1:2) = ubound(SrcMeshMapTypeData%LoadLn2_M)
      if (.not. allocated(DstMeshMapTypeData%LoadLn2_M)) then
         allocate(DstMeshMapTypeData%LoadLn2_M(LB(1):UB(1),LB(2):UB(2)), stat=ErrStat2)
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal, 'Error allocating DstMeshMapTypeData%LoadLn2_M.', ErrStat, ErrMsg, RoutineName)
            return
         end if
      end if
      DstMeshMapTypeData%LoadLn2_M = SrcMeshMapTypeData%LoadLn2_M
   end if
   call NWTC_Library_CopyMeshMapLinearizationType(SrcMeshMapTypeData%dM, DstMeshMapTypeData%dM, CtrlCode, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return
end subroutine

subroutine NWTC_Library_DestroyMeshMapType(MeshMapTypeData, ErrStat, ErrMsg)
   type(MeshMapType), intent(inout) :: MeshMapTypeData
   integer(IntKi),  intent(  out) :: ErrStat
   character(*),    intent(  out) :: ErrMsg
   integer(B4Ki)   :: i1, i2, i3
   integer(B4Ki)   :: LB(3), UB(3)
   integer(IntKi)                 :: ErrStat2
   character(ErrMsgLen)           :: ErrMsg2
   character(*), parameter        :: RoutineName = 'NWTC_Library_DestroyMeshMapType'
   ErrStat = ErrID_None
   ErrMsg  = ''
   if (allocated(MeshMapTypeData%MapLoads)) then
      LB(1:1) = lbound(MeshMapTypeData%MapLoads)
      UB(1:1) = ubound(MeshMapTypeData%MapLoads)
      do i1 = LB(1), UB(1)
         call NWTC_Library_DestroyMapType(MeshMapTypeData%MapLoads(i1), ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end do
      deallocate(MeshMapTypeData%MapLoads)
   end if
   if (allocated(MeshMapTypeData%MapMotions)) then
      LB(1:1) = lbound(MeshMapTypeData%MapMotions)
      UB(1:1) = ubound(MeshMapTypeData%MapMotions)
      do i1 = LB(1), UB(1)
         call NWTC_Library_DestroyMapType(MeshMapTypeData%MapMotions(i1), ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end do
      deallocate(MeshMapTypeData%MapMotions)
   end if
   if (allocated(MeshMapTypeData%MapSrcToAugmt)) then
      LB(1:1) = lbound(MeshMapTypeData%MapSrcToAugmt)
      UB(1:1) = ubound(MeshMapTypeData%MapSrcToAugmt)
      do i1 = LB(1), UB(1)
         call NWTC_Library_DestroyMapType(MeshMapTypeData%MapSrcToAugmt(i1), ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end do
      deallocate(MeshMapTypeData%MapSrcToAugmt)
   end if
   call MeshDestroy( MeshMapTypeData%Augmented_Ln2_Src, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call MeshDestroy( MeshMapTypeData%Lumped_Points_Src, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (allocated(MeshMapTypeData%LoadLn2_A_Mat_Piv)) then
      deallocate(MeshMapTypeData%LoadLn2_A_Mat_Piv)
   end if
   if (allocated(MeshMapTypeData%DisplacedPosition)) then
      deallocate(MeshMapTypeData%DisplacedPosition)
   end if
   if (allocated(MeshMapTypeData%LoadLn2_A_Mat)) then
      deallocate(MeshMapTypeData%LoadLn2_A_Mat)
   end if
   if (allocated(MeshMapTypeData%LoadLn2_F)) then
      deallocate(MeshMapTypeData%LoadLn2_F)
   end if
   if (allocated(MeshMapTypeData%LoadLn2_M)) then
      deallocate(MeshMapTypeData%LoadLn2_M)
   end if
   call NWTC_Library_DestroyMeshMapLinearizationType(MeshMapTypeData%dM, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
end subroutine

subroutine NWTC_Library_PackMeshMapType(RF, Indata)
   type(RegFile), intent(inout) :: RF
   type(MeshMapType), intent(in) :: InData
   character(*), parameter         :: RoutineName = 'NWTC_Library_PackMeshMapType'
   integer(B4Ki)   :: i1, i2, i3
   integer(B4Ki)   :: LB(3), UB(3)
   if (RF%ErrStat >= AbortErrLev) return
   call RegPack(RF, allocated(InData%MapLoads))
   if (allocated(InData%MapLoads)) then
      call RegPackBounds(RF, 1, lbound(InData%MapLoads), ubound(InData%MapLoads))
      LB(1:1) = lbound(InData%MapLoads)
      UB(1:1) = ubound(InData%MapLoads)
      do i1 = LB(1), UB(1)
         call NWTC_Library_PackMapType(RF, InData%MapLoads(i1)) 
      end do
   end if
   call RegPack(RF, allocated(InData%MapMotions))
   if (allocated(InData%MapMotions)) then
      call RegPackBounds(RF, 1, lbound(InData%MapMotions), ubound(InData%MapMotions))
      LB(1:1) = lbound(InData%MapMotions)
      UB(1:1) = ubound(InData%MapMotions)
      do i1 = LB(1), UB(1)
         call NWTC_Library_PackMapType(RF, InData%MapMotions(i1)) 
      end do
   end if
   call RegPack(RF, allocated(InData%MapSrcToAugmt))
   if (allocated(InData%MapSrcToAugmt)) then
      call RegPackBounds(RF, 1, lbound(InData%MapSrcToAugmt), ubound(InData%MapSrcToAugmt))
      LB(1:1) = lbound(InData%MapSrcToAugmt)
      UB(1:1) = ubound(InData%MapSrcToAugmt)
      do i1 = LB(1), UB(1)
         call NWTC_Library_PackMapType(RF, InData%MapSrcToAugmt(i1)) 
      end do
   end if
   call MeshPack(RF, InData%Augmented_Ln2_Src) 
   call MeshPack(RF, InData%Lumped_Points_Src) 
   call RegPackAlloc(RF, InData%LoadLn2_A_Mat_Piv)
   call RegPackAlloc(RF, InData%DisplacedPosition)
   call RegPackAlloc(RF, InData%LoadLn2_A_Mat)
   call RegPackAlloc(RF, InData%LoadLn2_F)
   call RegPackAlloc(RF, InData%LoadLn2_M)
   call NWTC_Library_PackMeshMapLinearizationType(RF, InData%dM) 
   if (RegCheckErr(RF, RoutineName)) return
end subroutine

subroutine NWTC_Library_UnPackMeshMapType(RF, OutData)
   type(RegFile), intent(inout)    :: RF
   type(MeshMapType), intent(inout) :: OutData
   character(*), parameter            :: RoutineName = 'NWTC_Library_UnPackMeshMapType'
   integer(B4Ki)   :: i1, i2, i3
   integer(B4Ki)   :: LB(3), UB(3)
   integer(IntKi)  :: stat
   logical         :: IsAllocAssoc
   if (RF%ErrStat /= ErrID_None) return
   if (allocated(OutData%MapLoads)) deallocate(OutData%MapLoads)
   call RegUnpack(RF, IsAllocAssoc); if (RegCheckErr(RF, RoutineName)) return
   if (IsAllocAssoc) then
      call RegUnpackBounds(RF, 1, LB, UB); if (RegCheckErr(RF, RoutineName)) return
      allocate(OutData%MapLoads(LB(1):UB(1)),stat=stat)
      if (stat /= 0) then 
         call SetErrStat(ErrID_Fatal, 'Error allocating OutData%MapLoads.', RF%ErrStat, RF%ErrMsg, RoutineName)
         return
      end if
      do i1 = LB(1), UB(1)
         call NWTC_Library_UnpackMapType(RF, OutData%MapLoads(i1)) ! MapLoads 
      end do
   end if
   if (allocated(OutData%MapMotions)) deallocate(OutData%MapMotions)
   call RegUnpack(RF, IsAllocAssoc); if (RegCheckErr(RF, RoutineName)) return
   if (IsAllocAssoc) then
      call RegUnpackBounds(RF, 1, LB, UB); if (RegCheckErr(RF, RoutineName)) return
      allocate(OutData%MapMotions(LB(1):UB(1)),stat=stat)
      if (stat /= 0) then 
         call SetErrStat(ErrID_Fatal, 'Error allocating OutData%MapMotions.', RF%ErrStat, RF%ErrMsg, RoutineName)
         return
      end if
      do i1 = LB(1), UB(1)
         call NWTC_Library_UnpackMapType(RF, OutData%MapMotions(i1)) ! MapMotions 
      end do
   end if
   if (allocated(OutData%MapSrcToAugmt)) deallocate(OutData%MapSrcToAugmt)
   call RegUnpack(RF, IsAllocAssoc); if (RegCheckErr(RF, RoutineName)) return
   if (IsAllocAssoc) then
      call RegUnpackBounds(RF, 1, LB, UB); if (RegCheckErr(RF, RoutineName)) return
      allocate(OutData%MapSrcToAugmt(LB(1):UB(1)),stat=stat)
      if (stat /= 0) then 
         call SetErrStat(ErrID_Fatal, 'Error allocating OutData%MapSrcToAugmt.', RF%ErrStat, RF%ErrMsg, RoutineName)
         return
      end if
      do i1 = LB(1), UB(1)
         call NWTC_Library_UnpackMapType(RF, OutData%MapSrcToAugmt(i1)) ! MapSrcToAugmt 
      end do
   end if
   call MeshUnpack(RF, OutData%Augmented_Ln2_Src) ! Augmented_Ln2_Src 
   call MeshUnpack(RF, OutData%Lumped_Points_Src) ! Lumped_Points_Src 
   call RegUnpackAlloc(RF, OutData%LoadLn2_A_Mat_Piv); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%DisplacedPosition); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%LoadLn2_A_Mat); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%LoadLn2_F); if (RegCheckErr(RF, RoutineName)) return
   call RegUnpackAlloc(RF, OutData%LoadLn2_M); if (RegCheckErr(RF, RoutineName)) return
   call NWTC_Library_UnpackMeshMapLinearizationType(RF, OutData%dM) ! dM 
end subroutine
!ENDOFREGISTRYGENERATEDFILE
