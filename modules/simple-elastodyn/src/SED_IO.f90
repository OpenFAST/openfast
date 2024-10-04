!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2024  National Renewable Energy Laboratory
!
!    This file is part of Simplified-ElastoDyn (SED)
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
MODULE SED_IO

   USE SED_Types
   USE SED_Output_Params
   USE NWTC_Library

   implicit none

   integer(IntKi),   parameter   :: Method_RK4  = 1
   integer(IntKi),   parameter   :: Method_AB4  = 2
   integer(IntKi),   parameter   :: Method_ABM4 = 3


   real(ReKi),       parameter   :: SmallAngleLimit_Deg  =  15.0     ! Largest input angle considered "small" (used as a check on input data), degrees
   integer(IntKi),   parameter   :: MaxBl    = 3                     ! Maximum number of blades allowed in simulation

   integer(IntKi),   parameter   :: DOF_Az   = 1                     ! Rotor azimuth



contains

!---------------------------------------------------------------
!> Parse the input in the InFileInfo (FileInfo_Type data structure):
subroutine SED_ParsePrimaryFileData( InitInp, RootName, interval, FileInfo_In, InputFileData, UnEc, ErrStat, ErrMsg )
   type(SED_InitInputType),   intent(in   )  :: InitInp              !< Input data for initialization routine
   character(1024),           intent(in   )  :: RootName             !< root name for summary file
   real(DBKi),                intent(in   )  :: interval             !< timestep
   type(FileInfoType),        intent(in   )  :: FileInfo_In          !< The input file stored in a data structure
   type(SED_InputFile),       intent(inout)  :: InputFileData        !< The data for initialization
   integer(IntKi),            intent(  out)  :: UnEc                 !< The local unit number for this module's echo file
   integer(IntKi),            intent(  out)  :: ErrStat              !< Error status  from this subroutine
   character(*),              intent(  out)  :: ErrMsg               !< Error message from this subroutine

   ! local vars
   integer(IntKi)                            :: CurLine              !< current entry in FileInfo_In%Lines array
   integer(IntKi)                            :: i                    !< generic counter
   real(SiKi)                                :: TmpRe(10)            !< temporary 10 number array for reading values in from table
   integer(IntKi)                            :: ErrStat2             !< Temporary error status  for subroutine and function calls
   character(ErrMsgLen)                      :: ErrMsg2              !< Temporary error message for subroutine and function calls
   character(*),              parameter      :: RoutineName="SED_ParsePrimaryFileData"

      ! Initialize ErrStat
   ErrStat  = ErrID_None
   ErrMsg   = ""
   UnEc     = -1  ! No file


   CALL AllocAry( InputFileData%OutList, MaxOutPts, "Outlist", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !======  Simulation control  =========================================================================
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

      ! IntMethod - Integration method: {1: RK4, 2: AB4, or 3: ABM4} (-):
   call ParseVar ( FileInfo_In, CurLine, "IntMethod", InputFileData%IntMethod, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! DT - Time interval for aerodynamic calculations {or default} (s):
   call ParseVarWDefault ( FileInfo_In, CurLine, "DT", InputFileData%DT, interval, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return


   !======  Degrees of Freedom  =========================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
      ! GenDOF -  Generator DOF (flag)
   call ParseVar( FileInfo_In, CurLine, "GenDOF", InputFileData%GenDOF, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! YawDOF -  Yaw DOF (flag)
   call ParseVar( FileInfo_In, CurLine, "YawDOF", InputFileData%YawDOF, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return


   !======  Iniital Conditions  =========================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      ! Azimuth - Initial azimuth angle for blades (degrees)
   call ParseVar( FileInfo_In, CurLine, "Azimuth", InputFileData%Azimuth, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   InputFileData%Azimuth = InputFileData%Azimuth * D2R

      ! BlPitch - Blades initial pitch (degrees)
   call ParseVar( FileInfo_In, CurLine, "BlPitch", InputFileData%BlPitch, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   InputFileData%BlPitch = InputFileData%BlPitch * D2R

      ! RotSpeed - Initial or fixed rotor speed (rpm)
   call ParseVar( FileInfo_In, CurLine, "RotSpeed", InputFileData%RotSpeed, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   InputFileData%RotSpeed = InputFileData%RotSpeed * RPM2RPS

      ! NacYaw - Initial or fixed nacelle-yaw angle (degrees)
   call ParseVar( FileInfo_In, CurLine, "NacYaw", InputFileData%NacYaw, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   InputFileData%NacYaw = InputFileData%NacYaw * D2R

      ! PtfmPitch - Fixed pitch tilt rotational displacement of platform (degrees)
   call ParseVar( FileInfo_In, CurLine, "PtfmPitch", InputFileData%PtfmPitch, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   InputFileData%PtfmPitch = InputFileData%PtfmPitch * D2R


   !======  Turbine Configuration  ======================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      ! NumBl       - Number of blades (-)
   call ParseVar( FileInfo_In, CurLine, "NumBl", InputFileData%NumBl, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! TipRad      - The distance from the rotor apex to the blade tip (meters)
   call ParseVar( FileInfo_In, CurLine, "TipRad", InputFileData%TipRad, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! HubRad      - The distance from the rotor apex to the blade root (meters)
   call ParseVar( FileInfo_In, CurLine, "HubRad", InputFileData%HubRad, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! PreCone     - Blades cone angle (degrees)
   call ParseVar( FileInfo_In, CurLine, "PreCone", InputFileData%PreCone, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   InputFileData%PreCone = InputFileData%PreCone * D2R

      ! OverHang    - Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)
   call ParseVar( FileInfo_In, CurLine, "OverHang", InputFileData%OverHang, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! ShftTilt    - Rotor shaft tilt angle (degrees)
   call ParseVar( FileInfo_In, CurLine, "ShftTilt", InputFileData%ShftTilt, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return
   InputFileData%ShftTilt = InputFileData%ShftTilt * D2R

      ! Twr2Shft    - Vertical distance from the tower-top to the rotor shaft (meters)
   call ParseVar( FileInfo_In, CurLine, "Twr2Shft", InputFileData%Twr2Shft, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! TowerHt     - Height of tower above ground level [onshore] or MSL [offshore] (meters)
   call ParseVar( FileInfo_In, CurLine, "TowerHt", InputFileData%TowerHt, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return


   !======  Mass and Inertia  ===========================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      ! RotIner     - Rot inertia about rotor axis [blades + hub] (kg m^2)
   call ParseVar( FileInfo_In, CurLine, "RotIner", InputFileData%RotIner, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return

      ! GenIner     - Generator inertia about HSS (kg m^2)
   call ParseVar( FileInfo_In, CurLine, "GenIner", InputFileData%GenIner, ErrStat2, ErrMsg2, UnEc )
      if (Failed()) return


   !======  Drivetrain  =================================================================================
   if ( InputFileData%Echo )   WRITE(UnEc, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      ! GBoxRatio     - Gearbox ratio (-)
   call ParseVar( FileInfo_In, CurLine, "GBoxRatio", InputFileData%GBoxRatio, ErrStat2, ErrMsg2, UnEc )
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
end subroutine SED_ParsePrimaryFileData


!> Check inputdata
subroutine SEDInput_ValidateInput( InitInp, InputFileData, ErrStat, ErrMsg )
   type(SED_InitInputType),   intent(in   )  :: InitInp              !< Input data for initialization
   type(SED_InputFile),       intent(in   )  :: InputFileData        !< The data for initialization
   integer(IntKi),            intent(  out)  :: ErrStat              !< Error status  from this subroutine
   character(*),              intent(  out)  :: ErrMsg               !< Error message from this subroutine
   character(*),              parameter      :: RoutineName="SEDInput_ValidateInput"

      ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! InitInput checks
   if (InitInp%Linearize)  call SetErrStat(ErrID_Fatal,'AeroDisk cannot perform linearization analysis.',ErrStat,ErrMsg,RoutineName)
   if (InputFileData%DT       <= 0.0_DbKi)   call SetErrStat(ErrID_Fatal,'DT must not be negative.',     ErrStat,ErrMsg,RoutineName)
   if ((InputFileData%IntMethod /= Method_RK4) .and. (InputFileData%IntMethod /= Method_AB4) .and. (InputFileData%IntMethod /= Method_ABM4))   &
      call SetErrStat(ErrID_Fatal,'IntMethod must be '//trim(Num2LStr(Method_RK4))//': RK4, '//trim(Num2LStr(Method_AB4))//': AB4, or '//trim(Num2LStr(Method_ABM4))//': ABM4',     ErrStat,ErrMsg,RoutineName)

   ! initial settings check
   if (abs(InputFileData%Azimuth) > TwoPi)     &
      call SetErrStat(ErrID_Fatal,'Starting Azimuth must be between -360 and 360 degrees.',     ErrStat,ErrMsg,RoutineName)
   if ((InputFileData%BlPitch  <= -pi ) .or. (InputFileData%BlPitch > pi))    &
      call SetErrStat( ErrID_Fatal, 'BlPitch must be greater than -pi radians and '// &
                                   'less than or equal to pi radians (i.e., in the range (-180, 180] degrees).',ErrStat,ErrMsg,RoutineName)
   if (InputFileData%RotSpeed < 0_ReKi)   call SetErrStat(ErrID_Fatal,'RotSpeed must not be negative',   ErrStat,ErrMsg,RoutineName)
   IF ((InputFileData%NacYaw <= -pi) .or. (InputFileData%NacYaw > pi)) &
      call SetErrStat( ErrID_Fatal, 'NacYaw must be in the range (-pi, pi] radians (i.e., (-180, 180] degrees).',ErrStat,ErrMsg,RoutineName)
   if ( ABS( InputFileData%PtfmPitch ) > SmallAngleLimit_Deg*D2R )  &
      call SetErrStat( ErrID_Fatal, 'PtfmPitch must be between -'//TRIM(Num2LStr(SmallAngleLimit_Deg))//' and ' &
                                    //TRIM(Num2LStr(SmallAngleLimit_Deg))//' degrees.',ErrStat,ErrMsg,RoutineName)

   ! turbine configuration
   if ((InputFileData%NumBl < 1) .or. (InputFileData%NumBl > MaxBl))    &
      call SetErrStat( ErrID_Fatal, 'NumBl must be 1, 2, or 3.',ErrStat,ErrMsg,RoutineName)
   if (InputFileData%TipRad < 0.0_ReKi)   call SetErrStat(ErrID_Fatal,'TipRad must not be negative.',ErrStat,ErrMsg,RoutineName)
   if (InputFileData%HubRad < 0.0_ReKi)   call SetErrStat(ErrID_Fatal,'HubRad must not be negative.',ErrStat,ErrMsg,RoutineName)
   if (abs(InputFileData%PreCone) >= PiBy2)     &
      call SetErrStat( ErrID_Fatal, 'PreCone must be in the range (-pi/2, pi/2) '//&
                                   'radians (i.e., (-90, 90) degrees).',ErrStat,ErrMsg,RoutineName)
   if (abs(InputFileData%OverHang) > InputFileData%TipRad)     &
      call SetErrStat( ErrID_Fatal, 'Overhang larger than tip-radius.  Does your model make sense?',ErrStat,ErrMsg,RoutineName)
   if (abs(InputFileData%ShftTilt) > PiBy2)    &
      call SetErrStat(ErrID_Fatal,'ShftTilt must be between -pi/2 and pi/2 radians (i.e., in the range [-90, 90] degrees).',ErrStat,ErrMsg,RoutineName)
   if (InputFileData%Twr2Shft < 0.0_ReKi)    call SetErrStat(ErrID_Fatal,'Twr2Shft must not be negative.',ErrStat,ErrMsg,RoutineName)
   if ( InputFileData%TowerHt <= 0.0_ReKi)   call SetErrStat( ErrID_Fatal, 'TowerHt must be greater than zero.',ErrStat,ErrMsg,RoutineName )

   if ( InputFileData%TowerHt + InputFileData%Twr2Shft + InputFileData%OverHang*SIN(InputFileData%ShftTilt) <= InputFileData%TipRad )     &
      call SetErrStat( ErrID_Fatal, 'TowerHt + Twr2Shft + OverHang*SIN(ShftTilt) must be greater than TipRad.',ErrStat,ErrMsg,RoutineName)

   ! inertias
   if (InputFileData%RotIner < 0.0_ReKi)     call SetErrStat(ErrID_Fatal,'RotIner must not be negative.',ErrStat,ErrMsg,RoutineName)
   if (InputFileData%GenIner < 0.0_ReKi)     call SetErrStat(ErrID_Fatal,'GenIner must not be negative.',ErrStat,ErrMsg,RoutineName)

   !GBRatio  -- no sanity checks on this
end subroutine SEDInput_ValidateInput


!----------------------------------------------------------------------------------------------------------------------------------
!> this routine fills the AllOuts array, which is used to send data to the glue code to be written to an output file.
!! NOTE: AllOuts is ReKi, but most calculations in this module are in single precision. This requires a bunch of conversions at this
!! stage.
subroutine Calc_WriteOutput( u, p, x, dxdt, y, m, ErrStat, ErrMsg, CalcWriteOutput )
   type(SED_InputType),          intent(in   )  :: u                 !< The inputs at time T
   type(SED_ParameterType),      intent(in   )  :: p                 !< The module parameters
   type(SED_ContinuousStateType),intent(in   )  :: x                 !< Continuous states at t
   type(SED_ContinuousStateType),intent(in   )  :: dxdt              !< Derivative of continuous states at t
   type(SED_OutputType),         intent(in   )  :: y                 !< outputs
   type(SED_MiscVarType),        intent(inout)  :: m                 !< misc/optimization variables (for computing mesh transfers)
   integer(IntKi),               intent(  out)  :: ErrStat           !< The error status code
   character(*),                 intent(  out)  :: ErrMsg            !< The error message, if an error occurred
   logical,                      intent(in   )  :: CalcWriteOutput   !< flag that determines if we need to compute AllOuts (or just the reaction loads that get returned to ServoDyn)
   ! local variables
   character(*), parameter                      :: RoutineName = 'Calc_WriteOutput'
   integer(IntKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   real(ReKi)                                   :: Tmp3(3)
   real(ReKi)                                   :: Rxyz(3,3)         !< rotation matrix for x,y,z of local coordinates

   ! Initialize
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! return if we are not providing outputs
   if (.not. CalcWriteOutput) return

   ! Azimuth
   m%AllOuts( Azimuth  ) = x%QT( DOF_Az)
   call Zero2TwoPi(m%AllOuts( Azimuth  ))    ! modulo
   m%AllOuts( Azimuth  ) = m%AllOuts( Azimuth  ) * R2D

   ! speed
   m%AllOuts( RotSpeed ) = x%QDT(DOF_Az) * RPS2RPM
   m%AllOuts( GenSpeed ) = x%QDT(DOF_Az) * RPS2RPM    * p%GBoxRatio

   ! accel
   m%AllOuts( RotAcc   ) = dxdt%QDT(DOF_Az) * RPS2RPM
   m%AllOuts( GenAcc   ) = dxdt%QDT(DOF_Az) * RPS2RPM * p%GBoxRatio

   ! Yaw commands
   m%AllOuts( Yaw      ) = y%Yaw     * R2D
   m%AllOuts( YawRate  ) = y%YawRate * R2D

   ! BlPitch1
   m%AllOuts( BlPitch1 ) = u%BlPitchCom(1) * R2D
   if (p%NumBl > 1)  m%AllOuts( BlPitch2 ) = u%BlPitchCom(2) * R2D
   if (p%NumBl > 2)  m%AllOuts( BlPitch3 ) = u%BlPitchCom(3) * R2D

   ! LLS torqque and power
   m%AllOuts( RotTorq  ) = y%RotTrq * 0.001_ReKi ! LSShftTq  (kN-m)
   m%AllOuts( RotPwr   ) = y%RotPwr * 0.001_ReKi ! LSShftPwr (kN-m)
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
!! This routine was generated by Write_ChckOutLst.m using the parameters listed in OutListParameters.xlsx at 10-Aug-2022 08:44:32.
SUBROUTINE SetOutParam(OutList, p, ErrStat, ErrMsg )
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables

   CHARACTER(ChanLen),        INTENT(IN)     :: OutList(:)                        !< The list out user-requested outputs
   TYPE(SED_ParameterType),    INTENT(INOUT)  :: p                                 !< The module parameters
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

   CHARACTER(OutStrLenM1), PARAMETER  :: ValidParamAry(25) =  (/  &   ! This lists the names of the allowed parameters, which must be sorted alphabetically
                               "AZIMUTH  ","BLDPITCH1","BLDPITCH2","BLDPITCH3","BLPITCH1 ","BLPITCH2 ","BLPITCH3 ","GENACC   ", &
                               "GENSPEED ","HSSHFTA  ","HSSHFTV  ","LSSHFTPWR","LSSHFTTQ ","LSSTIPA  ","LSSTIPAXA","LSSTIPAXS", &
                               "LSSTIPV  ","LSSTIPVXA","LSSTIPVXS","ROTACC   ","ROTPWR   ","ROTSPEED ","ROTTORQ  ","YAW      ", &
                               "YAWRATE  "/)
   INTEGER(IntKi), PARAMETER :: ParamIndxAry(25) =  (/ &                            ! This lists the index into AllOuts(:) of the allowed parameters ValidParamAry(:)
                                  Azimuth ,  BlPitch1 ,  BlPitch2 ,  BlPitch3 ,  BlPitch1 ,  BlPitch2 ,  BlPitch3 ,    GenAcc , &
                                 GenSpeed ,    GenAcc ,  GenSpeed ,    RotPwr ,   RotTorq ,    RotAcc ,    RotAcc ,    RotAcc , &
                                 RotSpeed ,  RotSpeed ,  RotSpeed ,    RotAcc ,    RotPwr ,  RotSpeed ,   RotTorq ,       Yaw , &
                                  YawRate /)
   CHARACTER(ChanLen), PARAMETER :: ParamUnitsAry(25) =  (/  &  ! This lists the units corresponding to the allowed parameters
                               "(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg)    ","(deg/s^2)", &
                               "(rpm)    ","(deg/s^2)","(rpm)    ","(kW)     ","(kN-m)   ","(deg/s^2)","(deg/s^2)","(deg/s^2)", &
                               "(rpm)    ","(rpm)    ","(rpm)    ","(deg/s^2)","(kW)     ","(rpm)    ","(kN-m)   ","(deg)    ", &
                               "(deg/s)  "/)


      ! Initialize values
   ErrStat = ErrID_None
   ErrMsg = ""
   InvalidOutput = .FALSE.

   if (p%NumBl < 3)  InvalidOutput( BlPitch3 ) = .true.
   if (p%NumBl < 2)  InvalidOutput( BlPitch2 ) = .true.

   !-------------------------------------------------------------------------------------------------
   ! Allocate and set index, name, and units for the output channels
   ! If a selected output channel is not available in this module, set error flag.
   !-------------------------------------------------------------------------------------------------

   ALLOCATE ( p%OutParam(0:p%NumOuts) , STAT=ErrStat2 )
   IF ( ErrStat2 /= 0_IntKi )  THEN
      CALL SetErrStat( ErrID_Fatal,"Error allocating memory for the SimpleElastoDyn OutParam array.", ErrStat, ErrMsg, RoutineName )
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
END MODULE SED_IO
