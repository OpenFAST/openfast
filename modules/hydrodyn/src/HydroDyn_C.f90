!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2021 National Renewable Energy Lab 
!
! This file is part of HydroDyn.
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
MODULE HydroDyn_C

    USE ISO_C_BINDING
    USE HydroDyn
    USE HydroDyn_Types
    USE NWTC_Library

   IMPLICIT NONE

   PUBLIC :: HydroDyn_Init_C
   PUBLIC :: HydroDyn_CalcOutput_C
   PUBLIC :: HydroDyn_UpdateStates_C
   PUBLIC :: HydroDyn_End_C

   integer(IntKi),   parameter            :: InterpOrder = 2   !< quadratic interpolation

   !  All HydroDyn data is stored within the following data structures inside this module.  No data is
   !  stored within HydroDyn itself, but is instead passed in from this module.  This data is not available
   !  to the calling code unless explicitly passed through the interface (derived types such as these are
   !  non-trivial to pass through the c-bindings).
   type(HydroDyn_InputType),  allocatable :: u(:)              !< Inputs at T, T-dt, T-2*dt  (history kept for updating states) 
   type(HydroDyn_InitInputType)           :: InitInp           !< Initialization data
   type(HydroDyn_InitOutputType)          :: InitOutData       !< Initial output data -- Names, units, and version info.
   type(HydroDyn_ParameterType)           :: p                 !< Parameters
   type(HydroDyn_ContinuousStateType)     :: x                 !< continuous states at Time t
   type(HydroDyn_DiscreteStateType)       :: xd                !< discrete states   at Time t
   type(HydroDyn_ConstraintStateType)     :: z                 !< Constraint states at Time t
   type(HydroDyn_OtherStateType)          :: OtherStates       !< Initial other/optimization states
   type(HydroDyn_OutputType)              :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
   type(HydroDyn_MiscVarType)             :: m                 !< Misc variables for optimization (not copied in glue code)

   !  If we are doing a correction step, the previous state information will be needed so we can compute with the new inputs
   type(HydroDyn_ContinuousStateType)     :: x_prev            !< continuous states at Time t of previous call
   type(HydroDyn_DiscreteStateType)       :: xd_prev           !< discrete states   at Time t of previous call
   type(HydroDyn_ConstraintStateType)     :: z_prev            !< Constraint states at Time t of previous call

   !  Mesh mapping
   type(MeshType)                         :: RefPtMesh         ! 1-node Point mesh located at (0,0,0) in global system where all PRP-related driver inputs are set
   type(MeshMapType)                      :: HD_Ref_2_WB_P     ! Mesh mapping between Reference pt mesh and WAMIT body(ies) mesh
   type(MeshMapType)                      :: HD_Ref_2_M_P      ! Mesh mapping between Reference pt mesh and Morison mesh
   real(R8Ki)                             :: theta(3)          ! mesh creation helper data
 
   !  Time tracking for when we are repeating a timestep
   integer(IntKi)                         :: T_Global          ! Time of this call 
   integer(IntKi)                         :: T_Global_prev     ! time of the previous call
   integer(IntKi)                         :: N_T_Global        ! count of which timestep we are on -- calculated internally based on the time passed into UpdateStates
   integer(IntKi)                         :: N_T_Global_prev   ! timestep of the previous call
   logical                                :: CorrectionStep    ! if we are repeating a timestep in UpdateStates,

   !  This must exactly match the value in the python-lib. If ErrMsgLen changes
   !  at some point in the nwtc-library, this should be updated, but the logic
   !  exists to correctly handle different lengths of the strings
   integer(IntKi),   parameter            :: ErrMsgLen_C=1025

CONTAINS

!> This routine sets the error status in C_CHAR for export to calling code.
!! Make absolutely certain that we do not overrun the end of ErrMsg_C.  That is hard coded to 1025,
!! but ErrMsgLen is set in the nwtc_library, and could change without updates here.  We don't want an
!! inadvertant buffer overrun -- that can lead to bad things.
subroutine SetErr(ErrStat, ErrMsg, ErrStat_C, ErrMsg_C)
   integer,                intent(in   )  :: ErrStat                 !< aggregated error message (fortran type)
   character(ErrMsgLen),   intent(in   )  :: ErrMsg                  !< aggregated error message (fortran type)
   integer(c_int),         intent(  out)  :: ErrStat_C
   character(kind=c_char), intent(  out)  :: ErrMsg_C(ErrMsgLen_C)
   ErrStat_C = ErrStat     ! We will send back the same error status that is used in OpenFAST
   if (ErrMsgLen > ErrMsgLen_C-1) then   ! If ErrMsgLen is > the space in ErrMsg_C, do not copy everything over
      ErrMsg_C = TRANSFER( trim(ErrMsg(1:ErrMsgLen_C-1))//C_NULL_CHAR, ErrMsg_C )
   else
      ErrMsg_C = TRANSFER( trim(ErrMsg)//C_NULL_CHAR, ErrMsg_C )
   endif
end subroutine SetErr


!===============================================================================================================
!--------------------------------------------- HydroDyn Init----------------------------------------------------
!===============================================================================================================
!FIXME: pass in interporder
SUBROUTINE HydroDyn_Init_c( InputFileString_C, InputFileStringLength_C,             &
               Gravity_C, defWtrDens_C, defWtrDpth_C, defMSL2SWL_C,                 &
               DT_C, TMax_C,                                                        &
               NumChannels_C, OutputChannelNames_C, OutputChannelUnits_C,           &
               ErrStat_C, ErrMsg_C) BIND (C, NAME='HydroDyn_Init_c')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_Init_c
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_Init_c
#endif

   type(c_ptr),               intent(in   )  :: InputFileString_C          !< Input file as a single string with lines deliniated by C_NULL_CHAR
   integer(c_int),            intent(in   )  :: InputFileStringLength_C    !< lenght of the input file string
   real(c_float),             intent(in   )  :: Gravity_C                  !< Gravitational constant (set by calling code)
   real(c_float),             intent(in   )  :: defWtrDens_C               !< Default value for water density (may be overridden by input file)
   real(c_float),             intent(in   )  :: defWtrDpth_C               !< Default value for water density (may be overridden by input file)
   real(c_float),             intent(in   )  :: defMSL2SWL_C               !< Default Offset between still-water level and mean sea level (m) [positive upward] (may be overridden by input file)
   real(c_double),            intent(in   )  :: DT_C                       !< Timestep used with HD for stepping forward from t to t+dt.  Must be constant.
   real(c_double),            intent(in   )  :: TMax_C                     !< Maximum time for simulation (used to set arrays for wave kinematics)
   integer(c_int),            intent(  out)  :: NumChannels_C              !< Number of output channels requested from the input file
   character(kind=c_char),    intent(  out)  :: OutputChannelNames_C(ChanLen*MaxHDOutputs+1)
   character(kind=c_char),    intent(  out)  :: OutputChannelUnits_C(ChanLen*MaxHDOutputs+1)
   integer(c_int),            intent(  out)  :: ErrStat_C                  !< Error status
   character(kind=c_char),    intent(  out)  :: ErrMsg_C(ErrMsgLen_C)      !< Error message (C_NULL_CHAR terminated)

   ! Local Variables
   character(kind=C_char, len=InputFileStringLength_C), pointer   :: InputFileString                  !< Input file as a single string with NULL chracter separating lines

   real(DbKi)                                                     :: TimeInterval                     !< timestep for HD 
   integer(IntKi)                                                 :: ErrStat                          !< aggregated error message
   character(ErrMsgLen)                                           :: ErrMsg                           !< aggregated error message
   integer(IntKi)                                                 :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)                                           :: ErrMsg2                          !< temporary error message from a call
   integer(IntKi)                                                 :: i,j,k
   character(*), parameter                                        :: RoutineName = 'HydroDyn_Init_C'  !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Get fortran pointer to C_NULL_CHAR deliniated input file as a string 
   call C_F_pointer(InputFileString_C, InputFileString)

   ! Get the data to pass to HD_Init
   call InitFileInfo(InputFileString, InitInp%PassedFileData, ErrStat2, ErrMsg2);   if (Failed())  return

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the InFileInfo data structure.
!   call Print_FileInfo_Struct( CU, InitInp%PassedFileData ) ! CU is the screen -- different number on different systems.


   ! Set other inputs for calling HydroDyn_Init
   InitInp%InputFile             = "passed_hd_file"         ! dummy
   InitInp%OutRootName           = "HDpassed"               ! used for making echo files
   InitInp%UseInputFile          = .FALSE.
   InitInp%Linearize             = .FALSE.
   InitInp%HasIce                = .FALSE.
   InitInp%Linearize             = .FALSE.                  ! we may want to pass this flag to linearize at T=0 for added mass effects

!FIXME: allocate this based on driver input file
!  InitInp%WaveElevXY
   InitInp%PtfmLocationX         = 0.0_ReKi
   InitInp%PtfmLocationY         = 0.0_ReKi

   ! Values passed in
   InitInp%Gravity               = REAL(Gravity_C, ReKi)
   InitInp%defWtrDens            = REAL(defWtrDens_C, ReKi)
   InitInp%defWtrDpth            = REAL(defWtrDpth_C, ReKi)
   InitInp%defMSL2SWL            = REAL(defMSL2SWL_C, ReKi)
   TimeInterval                  = REAL(DT_C, DbKi)
   InitInp%TMax                  = REAL(TMax_C, DbKi)


   !-----------------------
   ! Allocate input array u
   !-----------------------
   !     These inputs are used in the time stepping algorithm within HD_UpdateStates
   !     For quadratic interpolation, 3 timesteps are used.  For linear, 2 timesteps
   !     (the HD code can handle either).
   !        u(1)  inputs at t
   !        u(2)  inputs at t -   dt
   !        u(3)  inputs at t - 2*dt
   allocate(u(3), STAT=ErrStat2)
      if (ErrStat2 /= 0) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Could not allocate WriteOutputHdr array"
         if (Failed())  return
      endif

!FIXME: add handling of input mesh info
print*,'FIXME: add input platform location information'

   ! Call the main subroutine HydroDyn_Init
   !     TimeInterval and InitInp are passed into HD_Init, all the rest are set by HD_Init
   !
   !     NOTE: Pass u(1) only (this is empty and will be set inside Init).  We will copy
   !           this to u(2) and u(3) afterwards
   call HydroDyn_Init( InitInp, u(1), p, x, xd, z, OtherStates, y, m, TimeInterval, InitOutData, ErrStat2, ErrMsg2 )
      if (Failed())  return

print*,'FIXME: add mesh handling and mapping' 

!FIXME: copy u(1) to others
!FIXME: what from InitOutData should be returned
!FIXME: what error handling is necessary 



   !--------------------------------------------------
   !  Set output channel informatoion for driver code.
   !--------------------------------------------------

   ! Number of channels
   NumChannels_C = size(InitOutData%WriteOutputHdr)

   ! transfer the output channel names and units to c_char arrays for returning
   k=1
   do i=1,NumChannels_C
      do j=1,ChanLen    ! max length of channel name.  Same for units
         OutputChannelNames_C(k)=InitOutData%WriteOutputHdr(i)(j:j)
         OutputChannelUnits_C(k)=InitOutData%WriteOutputUnt(i)(j:j)
         k=k+1
      enddo
   enddo

   ! null terminate the string
   OutputChannelNames_C(k) = C_NULL_CHAR
   OutputChannelUnits_C(k) = C_NULL_CHAR


   call Cleanup()
   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed) then
         call Cleanup()
         call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
      endif
   end function Failed
   subroutine Cleanup()    ! NOTE: we are ignoring any error reporting from here
   end subroutine Cleanup
END SUBROUTINE HydroDyn_Init_c


!===============================================================================================================
!--------------------------------------------- HydroDyn CalcOutput ---------------------------------------------
!===============================================================================================================

SUBROUTINE HydroDyn_CalcOutput_c(Time_C,OutputChannelValues_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_CalcOutput_c')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_CalcOutput_c
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_CalcOutput_c
#endif
   real(c_double),          intent(in   ) :: Time_C
   real(c_float),           intent(  out) :: OutputChannelValues_C(p%NumOuts)
   integer(c_int),          intent(  out) :: ErrStat_C
   character(kind=c_char),  intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   real(DbKi)                             :: Time
   integer(IntKi)                         :: ErrStat                          !< aggregated error status 
   character(ErrMsgLen)                   :: ErrMsg                           !< aggregated error message
   integer(IntKi)                         :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)                   :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter                :: RoutineName = 'HydroDyn_CalcOutput_c' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   ! Convert the inputs from C to Fortrn
   Time = REAL(Time_C,DbKi)

!FIXME: Reshape position, velocity, acceleration and set mesh accordingly.

!FIXME:Set inputs array and extrap/interp necessary

   ! Call the main subroutine HydroDyn_CalcOutput to get the resulting forces and moments at time T 
   CALL HydroDyn_CalcOutput( Time, u, p, x, xd, z, OtherStates, y, m, ErrStat2, ErrMsg2 )
      if (Failed())  return

   ! Get the output channel info out of y
   OutputChannelValues_C = REAL(y%WriteOutput, C_FLOAT)

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

CONTAINS
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE HydroDyn_CalcOutput_c

!===============================================================================================================
!--------------------------------------------- HydroDyn UpdateStates ---------------------------------------------
!===============================================================================================================
SUBROUTINE HydroDyn_UpdateStates_c(Time_C,ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_UpdateStates_c')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_UpdateStates_c
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_UpdateStates_c
#endif
   real(c_double),          intent(in   ) :: Time_C
   integer(c_int),          intent(  out) :: ErrStat_C
   character(kind=c_char),  intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   real(DbKi)                 :: Time
   integer                    :: ErrStat                          !< aggregated error status 
   character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
   integer                    :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter    :: RoutineName = 'HydroDyn_UpdateStates_c' !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

!FIXME: set time as array
   ! Convert the inputs from C to Fortran
   Time = REAL(Time_C,DbKi)

!FIXME: do I need to have a flag for stepping forward in time? What if we are doing a correction step????
!Rotate values for inputs (do extrap interp if we are not passed new inputs here)

!FIXME: Reshape position, velocity, acceleration and set mesh inputs

!FIXME:Set inputs array and extrap/interp necessary

!   ! Call the main subroutine HydroDyn_UpdateStates to get the velocities
!   CALL HydroDyn_UpdateStates( Time, u, p, x, xd, z, OtherStates, y, m, ErrStat2, ErrMsg2 )
!      if (Failed())  return

!FIXME: what do we need to update after the call?
!     states are handled internally.  If we are doing correction steps, we may need the old states

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)
   end function Failed
END SUBROUTINE HydroDyn_UpdateStates_c

!===============================================================================================================
!--------------------------------------------------- HydroDyn End-----------------------------------------------
!===============================================================================================================
!  NOTE: the error handling in this routine is slightly different than the other routines

SUBROUTINE HydroDyn_End_C(ErrStat_C,ErrMsg_C) BIND (C, NAME='HydroDyn_End_c')
   implicit none
#ifndef IMPLICIT_DLLEXPORT
!DEC$ ATTRIBUTES DLLEXPORT :: HydroDyn_End_c
!GCC$ ATTRIBUTES DLLEXPORT :: HydroDyn_End_c
#endif
   integer(c_int),          intent(  out) :: ErrStat_C
   character(kind=c_char),  intent(  out) :: ErrMsg_C(ErrMsgLen_C)

   ! Local variables
   integer(IntKi)             :: i                                !< generic loop counter
   integer                    :: ErrStat                          !< aggregated error status 
   character(ErrMsgLen)       :: ErrMsg                           !< aggregated error message
   integer                    :: ErrStat2                         !< temporary error status  from a call
   character(ErrMsgLen)       :: ErrMsg2                          !< temporary error message from a call
   character(*), parameter    :: RoutineName = 'HydroDyn_End_c'   !< for error handling

   ! Initialize error handling
   ErrStat  =  ErrID_None
   ErrMsg   =  ""

   !  NOTE: HydroDyn_End only takes 1 instance of u, not the array.  So extra
   !        logic is required here (this isn't necessary in the fortran driver
   !        or in openfast, but may be when this code is called from C, Python,
   !        or some other code using the c-bindings.
   if (allocated(u)) then
      do i=2,size(u)    ! leave first one for passing to HD_End for destruction
         call HydroDyn_DestroyInput( u, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      enddo
   endif

   ! Call the main subroutine HydroDyn_End
   !     If u is not allocated, then we didn't get far at all in initialization,
   !     or HD_End_C got called before Init.  We don't want a segfault, so check
   !     for allocation.
   if (allocated(u)) then
      call HydroDyn_End( u(1), p, x, xd, z, OtherStates, y, m, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   endif

   ! if u is still allocated, deallocate it now
   if (allocated(u))    deallocate(u)

   call SetErr(ErrStat,ErrMsg,ErrStat_C,ErrMsg_C)

END SUBROUTINE HydroDyn_End_c

END MODULE HydroDyn_C
