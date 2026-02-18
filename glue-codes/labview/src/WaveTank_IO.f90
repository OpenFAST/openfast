!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2025  National Renewable Energy Laboratory
!
!    This file is a module specific to an experimental wave tank at NREL.
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
!
!  This code is designed to connect with LabView for a specific wave tank test case and likely will not work for other purposes.
!
!  For this test, a physical platform is deployed in a wave tank with cable acutators that are controlled through LabView.  This
!  module is called to provide some loads that are not present in the physical tank setup.  These include the following:
!     - rotor loading from a fixed RPM MHK rotor from AeroDyn.  This is calculated from either steady current provided by SeaState
!     - Mooring loads from MoorDyn
!
!
!**********************************************************************************************************************************
MODULE WaveTank_IO
   use ISO_C_BINDING
   use NWTC_Library
   use NWTC_IO
   use WaveTank_Types

   implicit none
   private

   public :: ParseInputFile
   public :: ValidateInputFile
   public :: TransferOutChanNamesUnits
   public :: InitOutputFile
   public :: WriteOutputLine

   ! These channels are output by default
   integer(IntKi), parameter :: NumDefChans = 42
   character(OutStrLenM1), parameter  :: DefChanNames(NumDefChans) =  (/ "Time          ",   &
                              "Ptfm_x        ","Ptfm_y        ","Ptfm_z        ",   &  ! position (absolute global)
                              "Ptfm_Rx       ","Ptfm_Ry       ","Ptfm_Rz       ",   &  ! Euler angles phi,theta,psi
                              "Ptfm_Vx       ","Ptfm_Vy       ","Ptfm_Vz       ",   &  ! translation vel
                              "Ptfm_RVx      ","Ptfm_RVy      ","Ptfm_RVz      ",   &  ! rotation vel
                              "Ptfm_Ax       ","Ptfm_Ay       ","Ptfm_Az       ",   &  ! translation acc
                              "Ptfm_RAx      ","Ptfm_RAy      ","Ptfm_RAz      ",   &  ! rotation acc
                              "Ptfm_Fx       ","Ptfm_Fy       ","Ptfm_Fz       ",   &  ! Forces total
                              "Ptfm_Mx       ","Ptfm_My       ","Ptfm_Mz       ",   &  ! Moments total
                              "Ptfm_MD_Fx    ","Ptfm_MD_Fy    ","Ptfm_MD_Fz    ",   &  ! Forces from MD
                              "Ptfm_MD_Mx    ","Ptfm_MD_My    ","Ptfm_MD_Mz    ",   &  ! Moments from MD
                              "Ptfm_ADI_Fx   ","Ptfm_ADI_Fy   ","Ptfm_ADI_Fz   ",   &  ! Forces from ADI
                              "Ptfm_ADI_Mx   ","Ptfm_ADI_My   ","Ptfm_ADI_Mz   ",   &  ! Moments from ADI
                              "Azimuth       ","RotSpeed      ","BlPitch       ",   &
                              "NacYaw        ","BuoyElev      "                     &
                           /)
   character(OutStrLenM1), parameter  :: DefChanUnits(NumDefChans) =  (/ "(s)           ",   &
                              "(m)           ","(m)           ","(m)           ",   &
                              "(rad)         ","(rad)         ","(rad)         ",   &
                              "(m/s)         ","(m/s)         ","(m/s)         ",   &
                              "(rad/s)       ","(rad/s)       ","(rad/s)       ",   &
                              "(m/s^2)       ","(m/s^2)       ","(m/s^2)       ",   &
                              "(rad/s^2)     ","(rad/s^2)     ","(rad/s^2)     ",   &
                              "(N)           ","(N)           ","(N)           ",   &  ! Forces total
                              "(N-m)         ","(N-m)         ","(N-m)         ",   &  ! Moments total
                              "(N)           ","(N)           ","(N)           ",   &  ! Forces from MD
                              "(N-m)         ","(N-m)         ","(N-m)         ",   &  ! Moments from MD
                              "(N)           ","(N)           ","(N)           ",   &  ! Forces from ADI
                              "(N-m)         ","(N-m)         ","(N-m)         ",   &  ! Moments from ADI
                              "(deg)         ","(RPM)         ","(deg)         ",   &
                              "(deg)         ","(m)           "                     &
                           /)

contains

subroutine ParseInputFile(FileInfo_In, SimSettings, ErrStat, ErrMsg)
   type(FileInfoType),        intent(in   )  :: FileInfo_In       !< The derived type for holding the file information.
   type(SimSettingsType),     intent(  out)  :: SimSettings
   integer(IntKi),            intent(  out)  :: ErrStat
   character(*),              intent(  out)  :: ErrMsg

   ! Local variables
   integer                                   :: CurLine
   character(1024), target                   :: TmpPath
   character(1024)                           :: FileName
   integer(IntKi)                            :: ErrStat2             ! local status of error message
   character(ErrMsgLen)                      :: ErrMsg2              ! local error message if errStat /= ErrID_None
   character(*), parameter                   :: RoutineName = 'WaveTankTesting.ParseInputFile'

   ErrStat = ErrID_None
   ErrMsg  = " "

   CurLine = 1
   ! Separator/header lines skipped automatically
   ! ----- Simulation control -------------
   call ParseVar( FileInfo_In, CurLine, 'DT',               SimSettings%Sim%DT,                    ErrStat2, ErrMsg2); if(Failed()) return;  ! timestep (unused)
   call ParseVar( FileInfo_In, CurLine, 'TMax',             SimSettings%Sim%TMax,                  ErrStat2, ErrMsg2); if(Failed()) return;  ! Max sim time (used only with SeaState wavemod 5)
   call ParseVar( FileInfo_In, CurLine, 'MHK',              SimSettings%Sim%MHK,                   ErrStat2, ErrMsg2); if(Failed()) return;  ! MHK turbine type (switch) {0=Not an MHK turbine; 1=Fixed MHK turbine; 2=Floating MHK turbine}
   call ParseVar( FileInfo_In, CurLine, 'InterpOrd',        SimSettings%Sim%InterpOrd,             ErrStat2, ErrMsg2); if(Failed()) return;  ! Interpolation order [unused]
!TODO: These are placeholders for later use.  Some of the logic is incomplete which is why this has been commented out.
!   call ParseVar( FileInfo_In, CurLine, 'ScaleFact',        SimSettings%Sim%ScaleFact,             ErrStat2, ErrMsg2); if(Failed()) return;  ! scaling factor for scaling full size model to wavetank scale results (Froude scaling: lambda = full_dimension / scale_dimension) [>1 expected] (-)
!   call ParseVar( FileInfo_In, CurLine, 'DensFact',         SimSettings%Sim%DensFact,              ErrStat2, ErrMsg2); if(Failed()) return;  ! ratio of density - Density_full/Density_model (rho_F/rho_M).  Used with Froude scaling of forces/moments" (-)
   call ParseVar( FileInfo_In, CurLine, 'DebugLevel',       SimSettings%Sim%DebugLevel,            ErrStat2, ErrMsg2); if(Failed()) return;  ! 0: none, 1: I/O summary, 2: +positions/orientations passed, 3:, 4: +all meshes
   call ParseVar( FileInfo_In, CurLine, 'OutRootName',      SimSettings%Sim%OutRootName,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Root name for any summary or other files
   ! -------- Environment ----------------
   call ParseVar( FileInfo_In, CurLine, 'Gravity',          SimSettings%Env%Gravity,               ErrStat2, ErrMsg2); if(Failed()) return;  ! Gravitational acceleration (m/s^2)
   call ParseVar( FileInfo_In, CurLine, 'WtrDens',          SimSettings%Env%WtrDens,               ErrStat2, ErrMsg2); if(Failed()) return;  ! Water density (kg/m^3)
   call ParseVar( FileInfo_In, CurLine, 'WtrVisc',          SimSettings%Env%WtrVisc,               ErrStat2, ErrMsg2); if(Failed()) return;  ! Kinematic viscosity of working fluid (m^2/s)
   call ParseVar( FileInfo_In, CurLine, 'SpdSound',         SimSettings%Env%SpdSound,              ErrStat2, ErrMsg2); if(Failed()) return;  ! Speed of sound in working fluid (m/s)
   call ParseVar( FileInfo_In, CurLine, 'Patm',             SimSettings%Env%Patm,                  ErrStat2, ErrMsg2); if(Failed()) return;  ! Atmospheric pressure (Pa) [used only for an MHK turbine cavitation check]
   call ParseVar( FileInfo_In, CurLine, 'Pvap',             SimSettings%Env%Pvap,                  ErrStat2, ErrMsg2); if(Failed()) return;  !  Vapour pressure of working fluid (Pa) [used only for an MHK turbine cavitation check]
   call ParseVar( FileInfo_In, CurLine, 'WtrDpth',          SimSettings%Env%WtrDpth,               ErrStat2, ErrMsg2); if(Failed()) return;  ! Water depth (m)
   call ParseVar( FileInfo_In, CurLine, 'MSL2SWL',          SimSettings%Env%MSL2SWL,               ErrStat2, ErrMsg2); if(Failed()) return;  ! Offset between still-water level and mean sea level (m) [positive upward]
   ! -------- SeaState -------------------
   call ParseVar( FileInfo_In, CurLine, 'SS_InputFile',     SimSettings%ModSettings%SS_InputFile,  ErrStat2, ErrMsg2); if(Failed()) return;  ! SeaState input file
   call ParseVar( FileInfo_In, CurLine, 'WaveTimeShift',    SimSettings%ModSettings%WaveTimeShift, ErrStat2, ErrMsg2); if(Failed()) return;  ! Shift the SeaState wavetime by this amount (for phase shifting waves to match tank)
   ! -------- MoorDyn --------------------
   call ParseVar( FileInfo_In, CurLine, 'MD_InputFile',     SimSettings%ModSettings%MD_InputFile,  ErrStat2, ErrMsg2); if(Failed()) return;  ! MoorDyn input file
   ! -------- AeroDyn + InflowWind -------
   call ParseVar( FileInfo_In, CurLine, 'AD_InputFile',     SimSettings%ModSettings%AD_InputFile,  ErrStat2, ErrMsg2); if(Failed()) return;  ! AeroDyn input file
   call ParseVar( FileInfo_In, CurLine, 'IfW_InputFile',    SimSettings%ModSettings%IfW_InputFile, ErrStat2, ErrMsg2); if(Failed()) return;  ! InflowWind input file
   ! -------- Turbine Configuration ------
   call ParseVar( FileInfo_In, CurLine, 'NumBl',            SimSettings%TrbCfg%NumBl,              ErrStat2, ErrMsg2); if(Failed()) return;  ! Number of blades (-)
   call ParseVar( FileInfo_In, CurLine, 'HubRad',           SimSettings%TrbCfg%HubRad,             ErrStat2, ErrMsg2); if(Failed()) return;  ! The distance from the rotor apex to the blade root (meters)
   call ParseVar( FileInfo_In, CurLine, 'PreCone',          SimSettings%TrbCfg%PreCone,            ErrStat2, ErrMsg2); if(Failed()) return;  ! Blade cone angle (degrees)
   SimSettings%TrbCfg%PreCone    = D2R * SimSettings%TrbCfg%PreCone
   call ParseVar( FileInfo_In, CurLine, 'OverHang',         SimSettings%TrbCfg%OverHang,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)
   call ParseVar( FileInfo_In, CurLine, 'ShftTilt',         SimSettings%TrbCfg%ShftTilt,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Rotor shaft tilt angle (degrees)
   SimSettings%TrbCfg%ShftTilt    = D2R * SimSettings%TrbCfg%ShftTilt
   call ParseVar( FileInfo_In, CurLine, 'Twr2Shft',         SimSettings%TrbCfg%Twr2Shft,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Vertical distance from the tower-top to the rotor shaft, center of nacelle (meters)
   call ParseVar( FileInfo_In, CurLine, 'TowerHt',          SimSettings%TrbCfg%TowerHt,            ErrStat2, ErrMsg2); if(Failed()) return;  ! Height of tower relative MSL
   call ParseAry( FileInfo_In, CurLine, 'TowerBsPt',        SimSettings%TrbCfg%TowerBsPt,     3,   ErrStat2, ErrMsg2); if(Failed()) return;  ! Height of tower base relative to PtfmRefPos in x,y, and water surface in z (meters)
   call ParseAry( FileInfo_In, CurLine, 'PtfmRefPos',       SimSettings%TrbCfg%PtfmRefPos,    3,   ErrStat2, ErrMsg2); if(Failed()) return;  ! Location of platform reference point, relative to MSL.  Motions and loads all connect to this point
   call ParseAry( FileInfo_In, CurLine, 'PtfmRefOrient',    SimSettings%TrbCfg%PtfmRefOrient, 3,   ErrStat2, ErrMsg2); if(Failed()) return;  ! Orientation of platform reference point, Euler angle set of roll,pitch,yaw" (deg)
   SimSettings%TrbCfg%PtfmRefOrient    = D2R * SimSettings%TrbCfg%PtfmRefOrient
   ! -------- Turbine Operating Point ----
   call ParseVar( FileInfo_In, CurLine, 'RotSpeed',         SimSettings%TrbInit%RotSpeed,          ErrStat2, ErrMsg2); if(Failed()) return;  ! Rotational speed of rotor in rotor coordinates (rpm)
   call ParseVar( FileInfo_In, CurLine, 'NacYaw',           SimSettings%TrbInit%NacYaw,            ErrStat2, ErrMsg2); if(Failed()) return;  ! Initial or fixed nacelle-yaw angle (deg read)
   call ParseVar( FileInfo_In, CurLine, 'BldPitch',         SimSettings%TrbInit%BldPitch,          ErrStat2, ErrMsg2); if(Failed()) return;  ! Blade 1 pitch (deg read)
   call ParseVar( FileInfo_In, CurLine, 'Azimuth',          SimSettings%TrbInit%Azimuth,           ErrStat2, ErrMsg2); if(Failed()) return;  ! Initial azimuth (deg read)
   ! angles are read as degrees, store as radians internally
   SimSettings%TrbInit%NacYaw    = D2R * SimSettings%TrbInit%NacYaw
   SimSettings%TrbInit%BldPitch  = D2R * SimSettings%TrbInit%BldPitch
   SimSettings%TrbInit%Azimuth   = D2R * SimSettings%TrbInit%Azimuth
   ! -------- Wave Buoy------------------
   call ParseAry( FileInfo_In, CurLine, 'WaveBuoyLoc',      SimSettings%WaveBuoy%XYLoc,      2,    ErrStat2, ErrMsg2); if(Failed()) return;  ! Location of the wave elevation measurement buoy. SeaState data returned at every timestep at this location (m)
   ! -------- Output ---------------------
   call ParseVar( FileInfo_In, CurLine, 'SendScreenToFile', SimSettings%Outs%SendScreenToFile,     ErrStat2, ErrMsg2); if(Failed()) return;  ! send to file <OutRootName>.screen.log if true
   call ParseVar( FileInfo_In, CurLine, 'OutFile',          SimSettings%Outs%OutFile,              ErrStat2, ErrMsg2); if(Failed()) return;  ! 0: no output file of channels, 1: output file in text format (at default DT) 
   call ParseVar( FileInfo_In, CurLine, 'OutFmt',           SimSettings%Outs%OutFmt,               ErrStat2, ErrMsg2); if(Failed()) return;  ! Format used for text tabular output, excluding the time channel. (quoted string)
   ! -------- VTK output -----------------
   call ParseVar( FileInfo_In, CurLine, 'WrVTK_Dir',        SimSettings%Viz%WrVTK_Dir,             ErrStat2, ErrMsg2); if(Failed()) return;  ! output directory for visualization
   call ParseVar( FileInfo_In, CurLine, 'WrVTK',            SimSettings%Viz%WrVTK,                 ErrStat2, ErrMsg2); if(Failed()) return;  ! VTK visualization data output: (switch) {0=none; 1=initialization data only; 2=animation; 3=mode shapes}
   call ParseVar( FileInfo_In, CurLine, 'WrVTK_type',       SimSettings%Viz%WrVTK_type,            ErrStat2, ErrMsg2); if(Failed()) return;  ! Type of VTK visualization data: (switch) {1=surfaces; 2=basic meshes (lines/points); 3=all meshes (debug)} [unused if WrVTK=0]
   call ParseVar( FileInfo_In, CurLine, 'WrVTK_DT',         SimSettings%Viz%WrVTK_DT,              ErrStat2, ErrMsg2); if(Failed()) return;  ! DT for writing VTK files
   call ParseAry( FileInfo_In, CurLine, 'VTKNacDim',        SimSettings%Viz%VTKNacDim,        6,   ErrStat2, ErrMsg2); if(Failed()) return;  !  Nacelle dimension passed in for VTK surface rendering [0,y0,z0,Lx,Ly,Lz] (m)

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine



subroutine ValidateInputFile(SimSettings, ErrStat, ErrMsg)
   type(SimSettingsType),     intent(inout)  :: SimSettings
   integer(IntKi),            intent(  out)  :: ErrStat
   character(*),              intent(  out)  :: ErrMsg
   integer(IntKi)                            :: ErrStat2
   character(ErrMsgLen)                      :: ErrMsg2
   character(*), parameter                   :: RoutineName = 'WaveTankTesting.ValidateInputFile'
   logical                                   :: file_exists

   ErrStat = ErrID_None
   ErrMsg  = ""

   !------------------------
   ! Sim Control
   !------------------------
   if (SimSettings%Sim%MHK /= 2_c_int)          call SetErrStat(ErrID_Fatal, "WaveTank module only works for floating MHK turbines at present (MHK=2).",ErrStat,ErrMsg,RoutineName)
   if (SimSettings%Sim%ScaleFact < 1.0_c_float) call SetErrStat(ErrID_Fatal, "ScaleFact should be > 1", ErrStat,ErrMsg,RoutineName)
   if (SimSettings%Sim%ScaleFact > 1.0_c_float) then
      call SetErrStat(ErrID_Warn,  "ScaleFact should be == 1 for now.  Scaling is untested and incomplete!!!!", ErrStat,ErrMsg,RoutineName)
      call WrScr("/---------------------------------------------------------\")
      call WrScr("|---  WARNING ----  WARNING  ---- WARNING ---- WARNING ---|")
      call WrScr("|---------------------------------------------------------|")
      call WrScr("|                                                         |")
      call WrScr("|  Froude scaling is not complete, nor is it tested!!!!   |")
      call WrScr("|                                                         |")
      call WrScr("|  At present, only some inputs are scaled, but equations |")
      call WrScr("|  have not been verified yet.  This is useful just for   |")
      call WrScr("|  observing motions are occuring, but will corrupt your  |")
      call WrScr("|  simulation.                                            |")
      call WrScr("|                                                         |")
      call WrScr("|  Set ScaleFact=1.0 in your input file.                  |")
      call WrScr("|                                                         |")
      call WrScr("|  TODO:                                                  |")
      call WrScr("|    - verify equations in FroudeScaling* functions       |")
      call WrScr("|    - scale resulting forces / moments                   |")
      call WrScr("|    - add scaled time, pos, vel, acc, frc, mom to output |")
      call WrScr("|      channels and add subscripting to differentiate     |")
      call WrScr("|                                                         |")
      call WrScr("\---------------------------------------------------------/")
   endif


   !------------------------
   ! Environment
   !------------------------

   !------------------------
   ! Turbine Config
   !------------------------

   !------------------------
   ! Input files
   !------------------------
   inquire(file=SimSettings%ModSettings%SS_InputFile,  exist=file_exists);  if (.not. file_exists) call SetErrStat(ErrID_Fatal,"Cannot find SeaState input file "  //trim(SimSettings%ModSettings%SS_InputFile ),ErrStat,ErrMsg,RoutineName)
   inquire(file=SimSettings%ModSettings%MD_InputFile,  exist=file_exists);  if (.not. file_exists) call SetErrStat(ErrID_Fatal,"Cannot find MoorDyn input file "   //trim(SimSettings%ModSettings%MD_InputFile ),ErrStat,ErrMsg,RoutineName) 
   inquire(file=SimSettings%ModSettings%AD_InputFile,  exist=file_exists);  if (.not. file_exists) call SetErrStat(ErrID_Fatal,"Cannot find AeroDyn input file "   //trim(SimSettings%ModSettings%AD_InputFile ),ErrStat,ErrMsg,RoutineName)
   inquire(file=SimSettings%ModSettings%IfW_InputFile, exist=file_exists);  if (.not. file_exists) call SetErrStat(ErrID_Fatal,"Cannot find InflowWind input file "//trim(SimSettings%ModSettings%IfW_InputFile),ErrStat,ErrMsg,RoutineName)
   if (ErrStat >= AbortErrLev) return

   ! Wave time shift must be positive
   if (SimSettings%ModSettings%WaveTimeShift < 0.0_DbKi) call SetErrStat(ErrID_Fatal, "WaveTimeShift must be >= 0",Errstat,ErrMsg,RoutineName)

   !------------------------
   ! Turbine Operating point
   !------------------------

   !------------------------
   ! Output
   !------------------------

   !------------------------
   ! VTK
   !------------------------

contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine


!> Transfer names or units from c character arrays into fortran character arrays for output file writing
subroutine TransferOutChanNamesUnits(NumChans, name, NamesUnits_C, NamesUnits, ErrStat, ErrMsg)
   integer(IntKi),                  intent(in   )  :: NumChans
   character(*),                    intent(in   )  :: name
   character(c_char),               intent(in   )  :: NamesUnits_C(:)
   character(ChanLen), allocatable, intent(  out)  :: NamesUnits(:)
   integer(IntKi),                  intent(  out)  :: ErrStat
   character(ErrMsgLen),            intent(  out)  :: ErrMsg
   integer(IntKi) :: i,idxStart,idxEnd
   call AllocAry(NamesUnits, NumChans, name, ErrStat, ErrMsg)
   if (ErrStat >= AbortErrLev) return
   do i=1,NumChans
      idxStart = (i-1)*ChanLen + 1
      idxEnd   = i*ChanLen
      NamesUnits(i) = transfer(NamesUnits_C(idxStart:idxEnd),NamesUnits(i))
   enddo
end subroutine


!> open the output file and populate the header
subroutine  InitOutputFile(WrOutputData, ErrStat, ErrMsg)
   type(WrOutputDataType), intent(inout)  :: WrOutputData
   integer(IntKi),         intent(  out)  :: ErrStat
   character(ErrMsgLen),   intent(  out)  :: ErrMsg
   integer(IntKi)          :: i
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'InitOutputFile'
   if (WrOutputData%OutUn > 0) then
      ErrStat = ErrID_Warn
      ErrMsg   = "Output file "//trim(WrOutputData%OutName)//" already open"
      return
   endif
   
   call GetNewUnit(WrOutputData%OutUn,ErrStat2,ErrMsg2)
   call OpenFOutFile (WrOutputData%OutUn, trim(WrOutputData%OutName), ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return

!FIXME: add version info here
   ! write header line
   write (WrOutputData%OutUn,'(/,A)')  'Simulation run on '//CurDate()//' at '//CurTime()//' using the WaveTank c-binding interface'
   write (WrOutputData%OutUn,'()')
   write (WrOutputData%OutUn,'()')
   write (WrOutputData%OutUn,'()')
   write (WrOutputData%OutUn,'()')
   
   !......................................................
   ! Write the names of the output parameters on one line:
   !......................................................
   ! default outputs
   call WrFileNR(WrOutputData%OutUn, DefChanNames(1))   ! time channel
   do i=2,NumDefChans
      call WrFileNR(WrOutputData%OutUn, tab//DefChanNames(i))
   enddo
   ! SS
   do i=1,WrOutputData%NumChans_SS
      call WrFileNR(WrOutputData%OutUn, tab//WrOutputData%WriteOutputHdr_SS(i))
   enddo
   ! MD
   do i=1,WrOutputData%NumChans_MD
      call WrFileNR(WrOutputData%OutUn, tab//WrOutputData%WriteOutputHdr_MD(i))
   enddo
   ! ADI
   do i=1,WrOutputData%NumChans_ADI
      call WrFileNR(WrOutputData%OutUn, tab//WrOutputData%WriteOutputHdr_ADI(i))
   enddo
   write (WrOutputData%OutUn,'()')

   !......................................................
   ! Write the units of the output parameters on one line:
   !......................................................
   ! default outputs
   call WrFileNR(WrOutputData%OutUn, DefChanUnits(1))   ! time channel
   do i=2,NumDefChans
      call WrFileNR(WrOutputData%OutUn, tab//DefChanUnits(i))
   enddo
   ! SS
   do i=1,WrOutputData%NumChans_SS
      call WrFileNR(WrOutputData%OutUn, tab//WrOutputData%WriteOutputUnt_SS(i))
   enddo
   ! MD
   do i=1,WrOutputData%NumChans_MD
      call WrFileNR(WrOutputData%OutUn, tab//WrOutputData%WriteOutputUnt_MD(i))
   enddo
   ! ADI
   do i=1,WrOutputData%NumChans_ADI
      call WrFileNR(WrOutputData%OutUn, tab//WrOutputData%WriteOutputUnt_ADI(i))
   enddo
   write (WrOutputData%OutUn,'()')

end subroutine


subroutine WriteOutputLine(OutFmt, CalcStepIO, StructTmp, WrOutputData, ErrStat, ErrMsg)
   character(*),              intent(in   )  :: OutFmt
   type(CalcStepIOdataType),  intent(in   )  :: CalcStepIO
   type(StructTmpType),       intent(in   )  :: StructTmp            ! operating states are in here
   type(WrOutputDataType),    intent(in   )  :: WrOutputData
   integer(IntKi),            intent(  out)  :: ErrStat
   character(ErrMsgLen),      intent(  out)  :: ErrMsg
   integer(IntKi)                            :: OutUnit
   integer(IntKi)                            :: errStat2             ! Status of error message (we're going to ignore errors in writing to the file)
   character(ErrMsgLen)                      :: errMsg2              ! Error message if ErrStat /= ErrID_None
   character(200)                            :: frmt                 ! A string to hold a format specifier
   character(15)                             :: tmpStr               ! temporary string to print the time output as text
   real(ReKi)                                :: TmpAry5(5)           ! temporary array for RotSpeed, BlPitch, NacYaw, Azimuth, BuoyWaveElev

   ErrStat = ErrID_None
   ErrMSg  = ""

   frmt = '"'//tab//'"'//trim(OutFmt)      ! format for array elements from individual modules
   OutUnit = WrOutputData%OutUn
   if (OutUnit <= 0_IntKi) then
      ErrStat = ErrID_Severe
      ErrMSg  = 'Cannot write to output file '//trim(WrOutputData%OutName)
      return
   endif
   
   ! time
   write( tmpStr, '(F15.6)' ) CalcStepIO%Time_c
   call WrFileNR( OutUnit, tmpStr )
   ! position / orientation euler angles, velocity, accel, resulting force/moment
   call WrNumAryFileNR(OutUnit, CalcStepIO%PosAng_c,          frmt, errStat2, errMsg2); if (Failed()) return
   call WrNumAryFileNR(OutUnit, CalcStepIO%Vel_c,             frmt, errStat2, errMsg2); if (Failed()) return
   call WrNumAryFileNR(OutUnit, CalcStepIO%Acc_c,             frmt, errStat2, errMsg2); if (Failed()) return
   call WrNumAryFileNR(OutUnit, CalcStepIO%FrcMom_c,          frmt, errStat2, errMsg2); if (Failed()) return ! total
   call WrNumAryFileNR(OutUnit, StructTmp%FrcMom_MD_at_Ptfm,  frmt, errStat2, errMsg2); if (Failed()) return
   call WrNumAryFileNR(OutUnit, StructTmp%FrcMom_ADI_at_Ptfm, frmt, errStat2, errMsg2); if (Failed()) return
   TmpAry5 = (/ R2D*StructTmp%Azimuth, StructTmp%RotSpeed, R2D*StructTmp%BldPitch, R2D*StructTmp%NacYaw, CalcStepIO%BuoyWaveElev /)
   call WrNumAryFileNR(OutUnit, TmpAry5,                      frmt, errStat2, errMsg2); if (Failed()) return
   ! channels from modules
   call WrNumAryFileNR(OutUnit, WrOutputData%OutData_SS,      frmt, errStat2, ErrMsg2); if (Failed()) return
   call WrNumAryFileNR(OutUnit, WrOutputData%OutData_MD,      frmt, errStat2, ErrMsg2); if (Failed()) return
   call WrNumAryFileNR(OutUnit, WrOutputData%OutData_ADI,     frmt, errStat2, ErrMsg2); if (Failed()) return
   ! write a new line (advance to the next line)
   write (OutUnit,'()')
contains
   logical function Failed()
        call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'WriteOutputLine') 
        Failed =  errStat >= AbortErrLev
   end function Failed
end subroutine
 
END MODULE WaveTank_IO
