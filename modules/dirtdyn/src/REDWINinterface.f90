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
   USE DirtDyn_Types
   USE, INTRINSIC :: ISO_C_Binding

   IMPLICIT NONE

   type(ProgDesc), parameter    :: REDWINinterface_Ver = ProgDesc( 'DirtDyn Interface for REDWIN soil interaction DLLs', 'using '//TRIM(OS_Desc), '99-Feb-2020' )
   
      !> Definition of the DLL Interface (from REDWIN):
      !! Note that aviFAIL and avcMSG should be used as INTENT(OUT), but I'm defining them INTENT(INOUT) just in case the compiler decides to reinitialize something that's INTENT(OUT)
   abstract interface
      subroutine REDWINdll_interface(PROPSFILE, LDISPFILE, IDTask, nErrorCode, ErrorCode, Props, StVar, StVarPrint, Disp, Force, D)
         USE, INTRINSIC :: ISO_C_Binding
         ! Define standard arguments for 'InterfaceFoundation'
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
      end subroutine REDWINdll_interface
   end interface

#ifdef STATIC_DLL_LOAD
   interface
      subroutine INTERFACEFOUNDATION ( PROPSFILE, LDISPFILE, IDTask, nErrorCode, ErrorCode, Props, StVar, StVarPrint, Disp, Force, D )  BIND(C, NAME='INTERFACEFOUNDATION')
         USE, INTRINSIC :: ISO_C_Binding
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


      ! Some constants for the Interface:
!   INTEGER(IntKi), PARAMETER    :: R_v36 = 85         !< Start of below-rated torque-speed look-up table (record no.) for REDWIN version 3.6
!   INTEGER(IntKi), PARAMETER    :: R_v4  = 145        !< Start of below-rated torque-speed look-up table (record no.) for REDWIN version 3.8 and later
!
!   INTEGER(IntKi), PARAMETER    :: R = R_v4           !< start of the generator speed look-up table  
            

CONTAINS
!==================================================================================================================================
!> This SUBROUTINE is used to call the REDWIN-style DLL.
subroutine CallREDWINdll ( u, DLL, dll_data, p, ErrStat, ErrMsg )

      ! Passed Variables:
   type(DirtD_InputType),     intent(in   )  :: u              ! System inputs
   type(DLL_Type),            intent(in   )  :: DLL            ! The DLL to be called.
   type(REDWINdllType),       intent(inout)  :: dll_data       ! data type containing the avrSWAP, accINFILE, and avcOUTNAME arrays 
   type(DirtD_ParameterType), intent(in   )  :: p              ! Parameters

   integer(IntKi),            intent(  out)  :: ErrStat        ! Error status of the operation
   character(*),              intent(  out)  :: ErrMsg         ! Error message if ErrStat /= ErrID_None
   
!FIXME: update these
      ! Local Variables:
   character(kind=c_char)                    :: accINFILE(LEN_TRIM(p%DLL_InFile)+1)  ! INFILE
   character(kind=c_char)                    :: avcOUTNAME(LEN_TRIM(p%RootName)+1)   ! OUTNAME (Simulation RootName)
 
   PROCEDURE(REDWINdll_interface),   POINTER :: REDWIN_Subroutine                 ! The address of the procedure in the RedWin DLL

!FIXME: not sure if this step is needed for the REDWIN DLLs
      !Convert to C-type characters: the "C_NULL_CHAR" converts the Fortran string to a C-type string (i.e., adds //CHAR(0) to the end)
   accINFILE  = TRANSFER( TRIM(p%DLL_InFile)//C_NULL_CHAR, accINFILE  )
   avcOUTNAME = TRANSFER( TRIM(p%RootName)//C_NULL_CHAR,   avcOUTNAME )
   
#ifdef STATIC_DLL_LOAD
      ! if we're statically loading the library (i.e., OpenFOAM), we can just call INTERFACEFOUNDATION(); 
   CALL INTERFACEFOUNDATION( accINFILE, avcOUTNAME, & !PROPSFILE, LDISPFILE, &
         dll_data%IDTask, dll_data%nErrorCode, dll_data%ErrorCode, &
         dll_data%Props, dll_data%StVar, dll_data%StVarPrint, &
         dll_data%Disp, dll_data%Force, dll_data%D )
#else
      ! Call the DLL (first associate the address from the procedure in the DLL with the subroutine):
   CALL C_F_PROCPOINTER( DLL%ProcAddr(1), REDWIN_Subroutine) 
   CALL REDWIN_Subroutine ( accINFILE, avcOUTNAME, & !PROPSFILE, LDISPFILE, &
         dll_data%IDTask, dll_data%nErrorCode, dll_data%ErrorCode, &
         dll_data%Props, dll_data%StVar, dll_data%StVarPrint, &
         dll_data%Disp, dll_data%Force, dll_data%D ) 
#endif


!FIXME: add error trapping routine that parses the ErrorCodes.

   
   return
end subroutine CallREDWINdll
!==================================================================================================================================
!> This routine initializes variables used in the REDWIN DLL interface.
subroutine REDWINinterface_Init(u,p,m,y,InputFileData, ErrStat, ErrMsg)
   
   type(DirtD_InputType),          intent(inout)  :: u               !< An initial guess for the input; input mesh must be defined
   type(DirtD_ParameterType),      intent(inout)  :: p               !< Parameters
   type(DirtD_MiscVarType),        intent(inout)  :: m               !< Initial misc (optimization) variables
   type(DirtD_OutputType),         intent(inout)  :: y               !< Initial system outputs (outputs are not calculated;
                                                                    !!   only the output mesh is initialized)
   type(DirtD_InputFile),          intent(inout)  :: InputFileData   !< Data stored in the module's input file
   integer(IntKi),                 intent(  out)  :: ErrStat         !< Error status of the operation
   character(*),                   intent(  out)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

      ! local variables
   integer(IntKi)                                  :: ErrStat2       ! The error status code
   character(ErrMsgLen)                            :: ErrMsg2        ! The error message, if an error occurred
      

   
   ErrStat = ErrID_None
   ErrMsg= ''
   
   CALL DispNVD( REDWINinterface_Ver )  ! Display the version of this interface
   
!   p%Ptch_Cntrl        = InputFileData%Ptch_Cntrl
   p%DLL_InFile        = InputFileData%DLL_InFile
   



!FIXME: do we need to set one?
!   p%DLL_DT            = InputFileData%DLL_DT
!   IF ( .NOT. EqualRealNos( NINT( p%DLL_DT / p%DT ) * p%DT, p%DLL_DT ) ) THEN
!      CALL CheckError( ErrID_Fatal, 'DLL_DT must be an integer multiple of DT.' )
!   END IF
!   IF ( p%DLL_DT < EPSILON( p%DLL_DT ) ) THEN 
!      CALL CheckError( ErrID_Fatal, 'DLL_DT must be larger than zero.' )
!   END IF
!   IF ( ErrStat >= AbortErrLev ) RETURN

!FIXME: initialize whatever here (mesh mappings etc)
   
#ifdef STATIC_DLL_LOAD
      ! because OpenFOAM needs the MPI task to copy the library, we're not going to dynamically load it; it needs to be loaded at runtime.
   p%DLL_Trgt%FileName = ''
   p%DLL_Trgt%ProcName = ''
#else
   ! Define and load the DLL:
   p%DLL_Trgt%FileName = InputFileData%DLL_FileName
   p%DLL_Trgt%ProcName = "" ! initialize all procedures to empty so we try to load only one
   p%DLL_Trgt%ProcName(1) = InputFileData%DLL_ProcName
   CALL LoadDynamicLib ( p%DLL_Trgt, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN
#endif
      
   ! Set status flag:
   p%UseREDWINinterface = .TRUE.

!FIXME: do I need to call it once to get it to actually initialize?  And when?
   !CALL CallREDWINdll(p%DLL_Trgt,  m%dll_data, ErrStat2, ErrMsg2)
   !   CALL CheckError(ErrStat2,ErrMsg2)
   !   IF ( ErrStat >= AbortErrLev ) RETURN
CONTAINS   
   !...............................................................................................................................
   subroutine CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

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
!> This routine calls the DLL for the final time (if it was previously called), and frees the dynamic library.
subroutine REDWINinterface_End(u, p, m, ErrStat, ErrMsg)
   
   TYPE(DirtD_InputType),           INTENT(IN   )  :: u               !< System inputs
   TYPE(DirtD_ParameterType),       INTENT(INOUT)  :: p               !< Parameters
   TYPE(DirtD_MiscVarType),         INTENT(INOUT)  :: m               !< misc (optimization) variables
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat         !< Error status of the operation
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

      ! local variables:
   INTEGER(IntKi)                                 :: ErrStat2    ! The error status code
   CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! The error message, if an error occurred
   
      ! call DLL final time, but skip if we've never called it
!FIXME: is there an end routine for the REDWIN DLLs????
!         CALL CallREDWINdll(u, p%DLL_Trgt,  m%dll_data, p, ErrStat, ErrMsg)
 
   CALL FreeDynamicLib( p%DLL_Trgt, ErrStat2, ErrMsg2 )  ! this doesn't do anything #ifdef STATIC_DLL_LOAD  because p%DLL_Trgt is 0 (NULL)
   IF (ErrStat2 /= ErrID_None) THEN  
      ErrStat = MAX(ErrStat, ErrStat2)      
      ErrMsg = TRIM(ErrMsg)//NewLine//TRIM(ErrMsg2)
   END IF
   
end subroutine REDWINinterface_End
!==================================================================================================================================
!> This routine sets the AVRswap array, calls the routine from the REDWIN DLL, and sets the outputs from the call to be used as
!! necessary in the main ServoDyn CalcOutput routine.
subroutine REDWINinterface_CalcOutput(t, u, p, m, ErrStat, ErrMsg)

   real(DbKi),                     intent(in   )  :: t           !< Current simulation time in seconds
   type(DirtD_InputType),          intent(in   )  :: u           !< Inputs at t
   type(DirtD_ParameterType),      intent(in   )  :: p           !< Parameters
   type(DirtD_MiscVarType),        intent(inout)  :: m           !< misc (optimization) variables
   integer(IntKi),                 intent(  out)  :: ErrStat     !< Error status of the operation
   character(*),                   intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      
      ! local variables:
   integer(IntKi)                                 :: ErrStat2    ! The error status code
   character(ErrMsgLen)                           :: ErrMsg2     ! The error message, if an error occurred
   
   
      ! Initialize error values:   
   ErrStat = ErrID_None
   ErrMsg= ''
   
   
#ifdef DEBUG_REDWIN_INTERFACE
!CALL WrNumAryFileNR ( 58, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat, ErrMsg )
!write(58,'()')
#endif
   
      ! Call the REDWIN-style DLL:
   CALL CallREDWINdll(u, p%DLL_Trgt,  m%dll_data, p, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN

#ifdef DEBUG_REDWIN_INTERFACE
!CALL WrNumAryFileNR ( 59, m%dll_data%avrSWAP,'1x,ES15.6E2', ErrStat, ErrMsg )
!write(59,'()')
#endif



end subroutine REDWINinterface_CalcOutput  
end module REDWINinterface
