!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of FAST's Controls and Electrical Drive Module, "ServoDyn".
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
MODULE REDWINinterface

!FIXME: when done, remove the ifdef NO_LibLoad checks.
   USE NWTC_Library
   USE SoilDyn_Types

   IMPLICIT NONE

   INTEGER(IntKi),   PARAMETER   :: IDtask_unkown  = 0_IntKi      ! Unknown task (placeholder for error checking)
   INTEGER(IntKi),   PARAMETER   :: IDtask_init    = 1_IntKi      ! Initialize DLL
   INTEGER(IntKi),   PARAMETER   :: IDtask_calc    = 2_IntKi      ! Calculate resultant force
   INTEGER(IntKi),   PARAMETER   :: IDtask_stiff   = 3_IntKi      ! Return stiffness 6x6

      !> Definition of the DLL Interface (from REDWIN):
   abstract interface
      subroutine REDWINdll_interface_v00(PROPSFILE, LDISPFILE, IDTask, nErrorCode, ErrorCode, Props, StVar, StVarPrint, Disp, Force, D)
         USE, INTRINSIC :: ISO_C_Binding, only : C_INT, C_CHAR, C_DOUBLE
         character(kind=c_char), intent(in   )  :: PROPSFILE(45)
         character(kind=c_char), intent(in   )  :: LDISPFILE(45)
         integer(c_int),         intent(in   )  :: IDTask
         integer(c_int),         intent(  out)  :: nErrorCode
         real(c_double),         intent(inout)  :: Props(1:100, 1:200)
         real(c_double),         intent(inout)  :: StVar(1:12, 1:100)
         integer(c_int),         intent(inout)  :: StVarPrint(1:12, 1:100)
         real(c_double),         intent(in   )  :: Disp(1:6)
         real(c_double),         intent(  out)  :: Force(1:6)
         real(c_double),         intent(  out)  :: D(1:6,1:6)
         integer(c_int),         intent(  out)  :: ErrorCode(1:100)
      end subroutine REDWINdll_interface_v00
   end interface

#ifdef STATIC_DLL_LOAD
   interface
         ! DO NOT REMOVE or MODIFY LINES starting with "!DEC$" or "!GCC$"
         ! !DEC$ specifies attributes for IVF and !GCC$ specifies attributes for gfortran
         ! NOTE: BIND(C... does not appear to be built into the DLL from REDWIN.
      subroutine INTERFACEFOUNDATION ( PROPSFILE, LDISPFILE, IDTask, nErrorCode, ErrorCode, Props, StVar, StVarPrint, Disp, Force, D )  !BIND(C, NAME='INTERFACEFOUNDATION')
         !DEC$ ATTRIBUTES DLLIMPORT :: INTERFACEFOUNDATION
         !GCC$ ATTRIBUTES DLLIMPORT :: INTERFACEFOUNDATION
         USE, INTRINSIC :: ISO_C_Binding, only : C_INT, C_CHAR, C_DOUBLE
         character(kind=c_char), intent(in   )  :: PROPSFILE(45)
         character(kind=c_char), intent(in   )  :: LDISPFILE(45)
         integer(c_int),         intent(in   )  :: IDTask
         integer(c_int),         intent(  out)  :: nErrorCode
         real(c_double),         intent(inout)  :: Props(1:100, 1:200)
         real(c_double),         intent(inout)  :: StVar(1:12, 1:100)
         integer(c_int),         intent(inout)  :: StVarPrint(1:12, 1:100)
         real(c_double),         intent(inout)  :: Disp(1:6)
         real(c_double),         intent(inout)  :: Force(1:6)
         real(c_double),         intent(inout)  :: D(1:6,1:6)
         integer(c_int),         intent(inout)  :: ErrorCode(1:100)
      end subroutine INTERFACEFOUNDATION
   end interface
#endif

   type(ProgDesc), parameter    :: REDWINinterface_Ver = ProgDesc( 'SoilDyn Interface for REDWIN soil interaction DLLs', 'using '//TRIM(OS_Desc), '99-Feb-2020' )

      ! Interface version (in case we end up with multiple different versions supported at some later date)
   INTEGER(IntKi), PARAMETER    :: RW_v00 = 0         ! Version number
   INTEGER(IntKi), PARAMETER    :: RW_ver = RW_v00    ! Current version number (read from DLL file)


CONTAINS
!==================================================================================================================================
!> This SUBROUTINE is used to call the REDWIN-style DLL.
subroutine CallREDWINdll ( u, DLL, dll_data, p, ErrStat, ErrMsg )

      ! Passed Variables:
   type(SlD_InputType),       intent(in   )  :: u              ! System inputs
   type(DLL_Type),            intent(in   )  :: DLL            ! The DLL to be called.
   type(REDWINdllType),       intent(inout)  :: dll_data       ! data type containing the avrSWAP, accINFILE, and avcOUTNAME arrays
   type(SlD_ParameterType),   intent(in   )  :: p              ! Parameters

   integer(IntKi),            intent(  out)  :: ErrStat        ! Error status of the operation
   character(*),              intent(  out)  :: ErrMsg         ! Error message if ErrStat /= ErrID_None

      ! Local Variables:
   character(len=45)                    :: PROPSFILE  ! properties input file
   character(len=45)                    :: LDISPFILE  ! displacement input file

   PROCEDURE(REDWINdll_interface_V00),POINTER:: REDWIN_Subroutine_v00                 ! The address of the procedure in the RedWin DLL

      ! Set names of DLL input files to pass
   PROPSFILE = TRIM(dll_data%PROPSfile)
   LDISPFILE = TRIM(dll_data%LDISPfile)

      ! Check existance of DLL input files.  The DLL does not check this, and will
      ! catastrophically fail if they are not found.

#ifdef STATIC_DLL_LOAD
      ! if we're statically loading the library (i.e., OpenFOAM), we can just call INTERFACEFOUNDATION();
   CALL INTERFACEFOUNDATION( PROPSFILE, LDISPFILE, &
         dll_data%IDTask, dll_data%nErrorCode, dll_data%ErrorCode, &
         dll_data%Props, dll_data%StVar, dll_data%StVarPrint, &
         dll_data%Disp, dll_data%Force, dll_data%D )
#else
      ! Call the DLL (first associate the address from the procedure in the DLL with the subroutine):
   if (RW_Ver == RW_v00) then
      CALL C_F_PROCPOINTER( transfer(DLL%ProcAddr(1),C_NULL_FUNPTR), REDWIN_Subroutine_v00)
      CALL REDWIN_Subroutine_v00 ( PROPSFILE, LDISPFILE, &
            dll_data%IDTask, dll_data%nErrorCode, dll_data%ErrorCode, &
            dll_data%Props, dll_data%StVar, dll_data%StVarPrint, &
            dll_data%Disp, dll_data%Force, dll_data%D )
   endif
#endif

      ! Call routine for error trapping the returned ErrorCodes
   call CheckREDWINerrors( dll_data, p%DLL_Model, dll_data%SuppressWarn, ErrStat, ErrMsg )
   return
end subroutine CallREDWINdll


!==================================================================================================================================
!> This routine initializes variables used in the REDWIN DLL interface.
subroutine REDWINinterface_Init(u,p,dll_data,y,InputFileData, ErrStat, ErrMsg)

   type(SlD_InputType),            intent(inout)  :: u               !< An initial guess for the input; input mesh must be defined
   type(SlD_ParameterType),        intent(inout)  :: p               !< Parameters
   type(REDWINdllType),            intent(inout)  :: dll_data
   type(SlD_OutputType),           intent(inout)  :: y               !< Initial system outputs (outputs are not calculated;
                                                                     !!   only the output mesh is initialized)
   type(SlD_InputFile),            intent(inout)  :: InputFileData   !< Data stored in the module's input file
   integer(IntKi),                 intent(  out)  :: ErrStat         !< Error status of the operation
   character(*),                   intent(  out)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)                                 :: ErrStat2       ! The error status code
   character(ErrMsgLen)                           :: ErrMsg2        ! The error message, if an error occurred
   character(*), parameter                        :: RoutineName = 'REDWINinterface_Init'


   ErrStat = ErrID_None
   ErrMsg= ''

   CALL DispNVD( REDWINinterface_Ver )  ! Display the version of this interface
#ifdef NO_Libload
   CALL SetErrStat( ErrID_Warn,'   -->  Skipping LoadDynamicLib call for '//TRIM(InputFileData%DLL_FileName),ErrStat,ErrMsg,RoutineName )
#endif


   ! Load the DLL
#ifdef STATIC_DLL_LOAD
      ! because OpenFOAM needs the MPI task to copy the library, we're not going to dynamically load it; it needs to be loaded at runtime.
   p%DLL_Trgt%FileName = ''
   p%DLL_Trgt%ProcName = ''
#else
   ! Define and load the DLL:
   p%DLL_Trgt%FileName = InputFileData%DLL_FileName
   p%DLL_Trgt%ProcName = "" ! initialize all procedures to empty so we try to load only one
   p%DLL_Trgt%ProcName(1) = InputFileData%DLL_ProcName
#ifndef NO_LibLoad
   CALL LoadDynamicLib ( p%DLL_Trgt, ErrStat2, ErrMsg2 );   if(Failed()) return;
#endif
#endif

      ! Initialize DLL
   dll_data%IDtask = IDtask_init
#ifndef NO_LibLoad
   CALL CallREDWINdll(u, p%DLL_Trgt, dll_data, p, ErrStat2, ErrMsg2);   if(Failed()) return;
#endif


!TODO: can we add a check on which type of library we actually loaded and compare to the model we set????
   ! Set status flag:
   p%UseREDWINinterface = .TRUE.

CONTAINS
   logical function Failed()
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      Failed =    ErrStat >= AbortErrLev
      if ( ErrStat >= AbortErrLev )    p%UseREDWINinterface = .FALSE.
   end function Failed
end subroutine REDWINinterface_Init


!==================================================================================================================================
!> This routine would call the DLL a final time, but there appears to be no end routine for the DLL,
!! so we don't need to make a last call.  It also frees the dynamic library (doesn't do anything on
!! static linked).
subroutine REDWINinterface_End(u, p, ErrStat, ErrMsg)

   TYPE(SlD_InputType),             INTENT(IN   )  :: u               !< System inputs
   TYPE(SlD_ParameterType),         INTENT(INOUT)  :: p               !< Parameters
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat         !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

      ! local variables:
   INTEGER(IntKi)                                 :: ErrStat2    ! The error status code
   CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! The error message, if an error occurred
   character(*), parameter                        :: RoutineName = 'REDWINinterface_End'

#ifdef NO_LibLoad
   ErrStat = ErrID_None
   ErrMsg= ''
   CALL SetErrStat( ErrID_Warn,'   -->  Skipping DynamicLib call for '//TRIM(p%DLL_Trgt%FileName),ErrStat,ErrMsg,RoutineName )
#else
      ! Free the library (note: this doesn't do anything #ifdef STATIC_DLL_LOAD  because p%DLL_Trgt is 0 (NULL))
   CALL FreeDynamicLib( p%DLL_Trgt, ErrStat, ErrMsg )
#endif
end subroutine REDWINinterface_End


!==================================================================================================================================
!> This routine sets the AVRswap array, calls the routine from the REDWIN DLL, and sets the outputs from the call to be used as
!! necessary in the main ServoDyn CalcOutput routine.
subroutine REDWINinterface_CalcOutput(t, u, p, dll_data, ErrStat, ErrMsg)

   real(DbKi),                     intent(in   )  :: t           !< Current simulation time in seconds
   type(SlD_InputType),            intent(in   )  :: u           !< Inputs at t
   type(SlD_ParameterType),        intent(in   )  :: p           !< Parameters
   type(REDWINdllType),            intent(inout)  :: dll_data
   integer(IntKi),                 intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                   intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables:
   integer(IntKi)                                 :: ErrStat2    ! The error status code
   character(ErrMsgLen)                           :: ErrMsg2     ! The error message, if an error occurred
   character(*), parameter                        :: RoutineName = 'REDWINinterface_CalcOutput'

      ! Initialize error values:
   ErrStat = ErrID_None
   ErrMsg= ''

!FIXME: coordinate transform and copy over u information here!


!FIXME: add some debugging options
#ifdef DEBUG_REDWIN_INTERFACE
!CALL WrNumAryFileNR ( 58, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat, ErrMsg )
!write(58,'()')
#endif

#ifdef NO_LibLoad
   CALL SetErrStat( ErrID_Warn,'   -->  Skipping DynamicLib call for '//TRIM(p%DLL_Trgt%FileName),ErrStat,ErrMsg,RoutineName )
#else
      ! Call the REDWIN-style DLL:
   dll_data%IDtask = IDtask_calc
   CALL CallREDWINdll(u, p%DLL_Trgt,  dll_data, p, ErrStat, ErrMsg); if(Failed()) return;
#endif


!FIXME: coordinate transform here!!!!

#ifdef DEBUG_REDWIN_INTERFACE
!CALL WrNumAryFileNR ( 59, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat, ErrMsg )
!write(59,'()')
#endif

contains
   logical function Failed()
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      Failed =    ErrStat >= AbortErrLev
   end function Failed
end subroutine REDWINinterface_CalcOutput


!==================================================================================================================================
!> This routine sets the AVRswap array, calls the routine from the REDWIN DLL, and sets the outputs from the call to be used as
!! necessary in the main ServoDyn CalcOutput routine.
subroutine REDWINinterface_GetStiffMatrix(t, u, p, dll_data, StiffMatrix, ErrStat, ErrMsg)

   real(DbKi),                   intent(in   )  :: t           !< Current simulation time in seconds
   type(SlD_InputType),          intent(in   )  :: u           !< Inputs at t
   type(SlD_ParameterType),      intent(in   )  :: p           !< Parameters
   type(REDWINdllType),          intent(inout)  :: dll_data
   real(ReKi),                   intent(  out)  :: StiffMatrix(6,6)
   integer(IntKi),               intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables:
   integer(IntKi)                               :: ErrStat2    ! The error status code
   character(ErrMsgLen)                         :: ErrMsg2     ! The error message, if an error occurred
   character(*), parameter                      :: RoutineName = 'REDWINinterface_CalcOutput'

      ! Initialize error values:
   ErrStat = ErrID_None
   ErrMsg= ''

!FIXME: coordinate transform and copy over u information here!

#ifdef NO_LibLoad
   CALL SetErrStat( ErrID_Warn,'   -->  Skipping DynamicLib call for '//TRIM(p%DLL_Trgt%FileName),ErrStat,ErrMsg,RoutineName )
   StiffMatrix = 0.0_ReKi
#else
      ! Call the REDWIN-style DLL:
   dll_data%IDtask = IDtask_stiff
   CALL CallREDWINdll(u, p%DLL_Trgt,  dll_data, p, ErrStat, ErrMsg); if(Failed()) return;
   StiffMatrix = real(dll_data%D,ReKi)    ! NOTE: converting types here
#endif

!FIXME: coordinate transform here!!!!

#ifdef DEBUG_REDWIN_INTERFACE
!CALL WrNumAryFileNR ( 59, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat, ErrMsg )
!write(59,'()')
#endif

contains
   logical function Failed()
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      Failed =    ErrStat >= AbortErrLev
   end function Failed
end subroutine REDWINinterface_GetStiffMatrix


!==================================================================================================================================
!> Check errors from REDWIN
!!    Error values taken from "20150014-11-R_Rev0_3D_Foundation Model Library.pdf"
subroutine CheckREDWINerrors( dll_data, DLL_Model, SuppressWarn, ErrStat, ErrMsg )
   type(REDWINdllType),       intent(in   )  :: dll_data       ! data type
   integer(IntKi),            intent(in   )  :: DLL_Model      ! Model type of the DLL
   logical,                   intent(inout)  :: SuppressWarn   ! from dll_data%SupressWarn
   integer(IntKi),            intent(  out)  :: ErrStat        ! Error status of the operation
   character(*),              intent(  out)  :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: i
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   ErrStat  =  ErrID_none
   ErrMsg   =  ''

   select case (DLL_Model)
      case(1)
         do i=1,dll_data%nErrorCode
            call CheckErrorsModel1(dll_data%ErrorCode(i),ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'REDWIN DLL error')
         enddo
      case(2)
         do i=1,dll_data%nErrorCode
            call CheckErrorsModel2(dll_data%ErrorCode(i),ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'REDWIN DLL error')
         enddo
      case(3)
         do i=1,dll_data%nErrorCode
            call CheckErrorsModel3(dll_data%ErrorCode(i),ErrStat2,ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'REDWIN DLL error')
         enddo
      case default
   end select

   ! Check if this is only a warning, and if we should supress further warnings (only one warning exists in each DLL model, rest are errors
   if (ErrStat == ErrID_Warn) then
      if ( SuppressWarn ) then
         ErrStat  =  ErrID_None
         ErrMsg   =  ''
      else
         SuppressWarn = .TRUE.
      endif
   endif

CONTAINS

   !> Check error codes from DLL model 1
   subroutine CheckErrorsModel1(ErrVal,ErrStat,ErrMsg)
            !  1     Warning: The number of rows in LDISDPFILE exceed the maximum number supported (200). The calibration will proceed using the first 200 values.
            !           Reduce the number of data points in the input file.
            !  2     Error in the interpolation tool used in the model calibration. The value you are trying to interpolate is outside the interpolation curve.
            !           Please inspect LDISPFILE to make sure that it covers a wide enough range and that all values are positive.
            !           Try to extend the input load-displacement curves in LDISPFILE.
      integer(IntKi),         intent(in   )  :: ErrVal
      integer(IntKi),         intent(  out)  :: ErrStat
      character(ErrMsgLen),   intent(  out)  :: ErrMSg
      integer(IntKi),         parameter      :: MaxErr=2

      if ( (ErrVal > MaxErr) .or. (ErrVal < 0) ) then
         ErrStat  =  ErrID_Fatal
         ErrMSg   =  'Unknown error from REDWIN DLL: '//trim(num2lstr(ErrVal))//'.  Only '//trim(num2lstr(MaxErr))//' values are known.' &
                     //NewLine//' --> Check that the correct REDWIN DLL model is specified and used.'
         return
      endif

      select case(ErrVal)
         case(0)
            ErrStat  =  ErrID_None
            ErrMsg   = ''
         case(1)
            ErrStat  =  ErrID_Warn
            ErrMsg   =  'The number of rows in LDISDPFILE exceed the maximum number supported (200). The calibration will proceed using the first 200 values.'  &
                        //NewLine//' --> Reduce the number of data points in the input file.'
         case(2)
            ErrStat  =  ErrID_Fatal
            ErrMsg   =  'Error in the interpolation tool used in the model calibration. The value you are trying to interpolate is outside the interpolation curve.' &
                        //NewLine//' --> Please inspect LDISPFILE to make sure that it covers a wide enough range and that all values are positive. ' &
                        //'Try to extend the input load-displacement curves in LDISPFILE.'
      end select
   end subroutine CheckErrorsModel1

   !> Check error codes from DLL model 2
   subroutine CheckErrorsModel2(ErrVal,ErrStat,ErrMsg)
            ! 1      Warning: The plastic force- displacement calibration curve has several zero-rows. The solution does not stop, but the results may be inaccurate or erroneous.
            !           Check that the provided coefficients of the elastic stiffness matrix are consistent with the load-displacement input curves.
            ! 2      Error. The iteration to find the plastic rotation increment and the plastic displacement increment did not converge.
            !           The force you are trying to apply might be outside the calibrated range. Please extend the input load-displacement curves in LDISPFILE.
            !           Alternatively, increase the number of iterations in PROPSFILE.
            ! 3      Error in the interpolation tool used in the model calibration. The value you are trying to interpolate is outside the interpolation curve.
            !           Please inspect LDISPFILE to make sure that it covers a wide enough range and that all values are positive. Try to extend the input load-displacement curves in LDISPFILE.
            ! 4      Error in the calibration tool. The contours of plastic horizontal displacement and the contours of plastic rotation are parallel.
            !           The input might be non-physical. Please check that LDISPFILE is in the correct format and that the units are consistent.
            ! 5      Error in the calibration tool. The calculation of the orientation of the yield surfaces might be wrong.
            !           The input might be non-physical. Please check that LDISPFILE is in the correct format and that the units are consistent.
            ! 6      Error in the calibration tool. The contours of plastic horizontal displacement are steeper than the contours of plastic rotation.
            !           The input might be non-physical. Please check that LDISPFILE is in the correct format and that the units are consistent.
      integer(IntKi),         intent(in   )  :: ErrVal
      integer(IntKi),         intent(  out)  :: ErrStat
      character(ErrMsgLen),   intent(  out)  :: ErrMSg
      integer(IntKi),         parameter      :: MaxErr=6

      if ( (ErrVal > MaxErr) .or. (ErrVal < 0) ) then
         ErrStat  =  ErrID_Fatal
         ErrMSg   =  'Unknown error from REDWIN DLL: '//trim(num2lstr(ErrVal))//'.  Only '//trim(num2lstr(MaxErr))//' values are known.'
         return
      endif

      select case(ErrVal)
         case(0)
            ErrStat  =  ErrID_None
            ErrMsg   = ''
         case(1)
            ErrStat  =  ErrID_Warn
            ErrMsg   =  'The plastic force- displacement calibration curve has several zero-rows. The solution does not stop, but the results may be inaccurate or erroneous.'   &
                        //NewLine//' --> Check that the provided coefficients of the elastic stiffness matrix are consistent with the load-displacement input curves.'
         case(2)
            ErrStat  =  ErrID_Fatal
            ErrMsg   =  'The iteration to find the plastic rotation increment and the plastic displacement increment did not converge.'  &
                        //NewLine//' --> The force you are trying to apply might be outside the calibrated range. Please extend the input load-displacement curves ' &
                        //'in LDISPFILE. Alternatively, increase the number of iterations in PROPSFILE.'
         case(3)
            ErrStat  =  ErrID_Fatal
            ErrMsg   =  'Error in the interpolation tool used in the model calibration. The value you are trying to interpolate is outside the interpolation curve.' &
                        //NewLine//' --> Please inspect LDISPFILE to make sure that it covers a wide enough range and that all values are positive. ' &
                        //'Try to extend the input load-displacement curves in LDISPFILE.'
         case(4)
            ErrStat  =  ErrID_Fatal
            ErrMsg   =  'Error in the calibration tool. The contours of plastic horizontal displacement and the contours of plastic rotation are parallel. '   &
                        //NewLine//' --> The input might be non-physical. Please check that LDISPFILE is in the correct format and that the units are consistent.'
         case(5)
            ErrStat  =  ErrID_Fatal
            ErrMsg   =  'Error in the calibration tool. The calculation of the orientation of the yield surfaces might be wrong.'   &
                        //NewLine//' --> The input might be non-physical. Please check that LDISPFILE is in the correct format and that the units are consistent.'
         case(6)
            ErrStat  =  ErrID_Fatal
            ErrMsg   =  'Error in the calibration tool. The contours of plastic horizontal displacement are steeper than the contours of plastic rotation.' &
                        //NewLine//' --> The input might be non-physical. Please check that LDISPFILE is in the correct format and that the units are consistent.'
      end select
   end subroutine CheckErrorsModel2

   !> Check error codes from DLL model 3
   subroutine CheckErrorsModel3(ErrVal,ErrStat,ErrMsg)
            ! 1      Warning. The solution in the current sub-step seems to be diverging. Will attempt to reduce the step size.
            !           The step size may be too large for convergence to be reached. The model will attempt to try again with a smaller step size.
            ! 2      Error. The sub-stepping algorithm in the multi-surface plasticity model did not converge.
            !           The cause of divergence is usually that the applied loads exceed the calibration range, or that there are several identical spring stiffness for low load levels.
            !           Possible solutions are: reduce the number of yield surfaces (Ns), increase the number of substeps (nsub), increase the range of the input load-displacement files.
            ! 3      Error in the calibration tool. The input file cannot be found.
            !           Check that the file name and path of the input files PROPSFILE and LDISPFILE are correctly specified.
            ! 4      Error in the calibration tool during read of PROPSFILE or LDISPFILE.
            !           Check that the format of the input files are correct.
      integer(IntKi),         intent(in   )  :: ErrVal
      integer(IntKi),         intent(  out)  :: ErrStat
      character(ErrMsgLen),   intent(  out)  :: ErrMSg
      integer(IntKi),         parameter      :: MaxErr=4

      if ( (ErrVal > MaxErr) .or. (ErrVal < 0) ) then
         ErrStat  =  ErrID_Fatal
         ErrMSg   =  'Unknown error from REDWIN DLL: '//trim(num2lstr(ErrVal))//'.  Only '//trim(num2lstr(MaxErr))//' values are known.'
         return
      endif

      select case(ErrVal)
         case(0)
            ErrStat  =  ErrID_None
            ErrMsg   = ''
         case(1)
            ErrStat  =  ErrID_Warn
            ErrMsg   =  'The solution in the current sub-step seems to be diverging. Will attempt to reduce the step size.' &
                        //NewLine//' --> The step size may be too large for convergence to be reached.' &
                        //' The model will attempt to try again with a smaller step size.'
         case(2)
            ErrMsg   =  'The sub-stepping algorithm in the multi-surface plasticity model did not converge.' &
                        //NewLine//' --> The cause of divergence is usually that the applied loads exceed the calibration range, or that there' &
                        //' are several identical spring stiffness for low load levels. Possible solutions are: reduce the number of yield surfaces' &
                        //'(Ns), increase the number of substeps (nsub), increase the range of the input load-displacement files.'
         case(3)
            ErrMsg   =  'Error in the calibration tool. The input file cannot be found.' &
                        //NewLine//' --> Check that the file name and path of the input files PROPSFILE and LDISPFILE are correctly specified.'
         case(4)
            ErrMsg   =  'Error in the calibration tool during read of PROPSFILE or LDISPFILE.' &
                        //NewLine//' --> Check that the format of the input files are correct.'
      end select
   end subroutine CheckErrorsModel3
end subroutine CheckREDWINerrors









end module REDWINinterface
