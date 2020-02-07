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

   USE NWTC_Library  
   USE SoilDyn_Types

   IMPLICIT NONE

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
   
#ifdef STATIC_DLL_LOAD
      ! because OpenFOAM needs the MPI task to copy the library, we're not going to dynamically load it; it needs to be loaded at runtime.
   p%DLL_Trgt%FileName = ''
   p%DLL_Trgt%ProcName = ''
#else
   ! Define and load the DLL:
   p%DLL_Trgt%FileName = InputFileData%DLL_FileName
   p%DLL_Trgt%ProcName = "" ! initialize all procedures to empty so we try to load only one
   p%DLL_Trgt%ProcName(1) = InputFileData%DLL_ProcName
#ifdef NO_LibLoad
   CALL SetErrStat( ErrID_Warn,'   -->  Skipping LoadDynamicLib call for '//TRIM(InputFileData%DLL_FileName),ErrStat,ErrMsg,RoutineName )
#else
   CALL LoadDynamicLib ( p%DLL_Trgt, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
#endif
#endif
      
   ! Set status flag:
   p%UseREDWINinterface = .TRUE.

!FIXME: do I need to call it once to get it to actually initialize?  And when?
#ifdef NO_LibLoad
   CALL SetErrStat( ErrID_Warn,'   -->  Skipping DynamicLib call for '//TRIM(p%DLL_Trgt%FileName),ErrStat,ErrMsg,RoutineName )
#else
      ! Initialize DLL 
   CALL CallREDWINdll(u, p%DLL_Trgt, dll_data, p, ErrStat2, ErrMsg2)
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
#endif

CONTAINS   
   ! Sets the error message and level and cleans up if the error is >= AbortErrLev
   subroutine CheckError(ErrID,Msg)

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN
         IF ( ErrStat /= ErrID_None ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'REDWINinterface_Init:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            p%UseREDWINinterface = .FALSE.
         END IF
      END IF


   end subroutine CheckError      
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
   CALL FreeDynamicLib( p%DLL_Trgt, ErrStat2, ErrMsg2 )
   IF (ErrStat2 /= ErrID_None) THEN  
      ErrStat = MAX(ErrStat, ErrStat2)      
      ErrMsg = TRIM(ErrMsg)//NewLine//TRIM(ErrMsg2)
   END IF
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

!FIXME: add some debugging options   
#ifdef DEBUG_REDWIN_INTERFACE
!CALL WrNumAryFileNR ( 58, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat, ErrMsg )
!write(58,'()')
#endif
   
#ifdef NO_LibLoad
   CALL SetErrStat( ErrID_Warn,'   -->  Skipping DynamicLib call for '//TRIM(p%DLL_Trgt%FileName),ErrStat,ErrMsg,RoutineName )
#else
      ! Call the REDWIN-style DLL:
   CALL CallREDWINdll(u, p%DLL_Trgt,  dll_data, p, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN
#endif

#ifdef DEBUG_REDWIN_INTERFACE
!CALL WrNumAryFileNR ( 59, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat, ErrMsg )
!write(59,'()')
#endif

end subroutine REDWINinterface_CalcOutput  
end module REDWINinterface
