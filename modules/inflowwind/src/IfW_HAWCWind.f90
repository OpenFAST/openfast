!>  This module uses full-field binary wind files to determine the wind inflow.
!!  This module assumes that the origin, (0,0,0), is located at the tower centerline at ground level,
!!  and that all units are specified in the metric system (using meters and seconds).
!!  Data is assumed periodic in the X direction (and thus not shifted like FFWind files are).
MODULE IfW_HAWCWind
!!
!!  Created 25-June-2010 by B. Jonkman, National Renewable Energy Laboratory
!!     using subroutines and modules from AeroDyn v12.58
!!
!!----------------------------------------------------------------------------------------------------
!! Updated 8-Aug-2015 for InflowWind v3.0 in the FAST v8 Framework
!
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
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
   USE                                          NWTC_Library
   USE                                          IfW_HAWCWind_Types
   USE                                          IfW_FFWind_Base

   IMPLICIT                                     NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_HAWCWind_Ver = ProgDesc( 'IfW_HAWCWind', '', '' )

   PUBLIC                                    :: IfW_HAWCWind_Init
   PUBLIC                                    :: IfW_HAWCWind_End
   PUBLIC                                    :: IfW_HAWCWind_CalcOutput

   INTEGER(IntKi), PARAMETER  :: nc = 3                           !< number of wind components
   
CONTAINS
!====================================================================================================
!>  This routine is used to initialize the parameters for using HAWC wind format files.
SUBROUTINE IfW_HAWCWind_Init(InitInp, p, MiscVars, Interval, InitOut, ErrStat, ErrMsg)

      ! Passed Variables
   TYPE(IfW_HAWCWind_InitInputType),         INTENT(IN   )  :: InitInp           !< Initialization data passed to the module
   TYPE(IfW_HAWCWind_ParameterType),         INTENT(  OUT)  :: p                 !< Parameters
   TYPE(IfW_HAWCWind_MiscVarType),           INTENT(  OUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)
   TYPE(IfW_HAWCWind_InitOutputType),        INTENT(  OUT)  :: InitOut           !< Initialization output

   REAL(DbKi),                               INTENT(IN   )  :: Interval          !< Time Interval to use (passed through here)


      ! Error Handling
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Message about errors


      ! Temporary variables for error handling
   INTEGER(IntKi)                                           :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg         ! temporary error message
   CHARACTER(*),                             PARAMETER      :: RoutineName = 'IfW_HAWCWind_Init'

      ! Local Variables:



   ErrStat = ErrID_None
   ErrMsg  = ""

   ! validate init input data:
   call ValidateInput(InitInp, TmpErrStat, TmpErrMsg)
      CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN
  
      
   !-------------------------------------------------------------------------------------------------
   ! Set some internal module parameters based on input file values
   !-------------------------------------------------------------------------------------------------
   p%FF%WindFileFormat  = 0
   p%FF%Periodic        = .true.
   p%FF%InterpTower     = .true.
   
   p%FF%NFFComp  = 3
   p%FF%NFFSteps = InitInp%nx
   p%FF%NYGrids  = InitInp%ny
   p%FF%NZGrids  = InitInp%nz
   p%FF%NTGrids  = 0

   p%FF%MeanFFWS        = InitInp%FF%URef
   p%FF%FFDTime         = InitInp%dx / InitInp%FF%URef
   p%FF%FFRate          = 1.0 / p%FF%FFDTime
   p%FF%InvFFYD         = 1.0 / InitInp%dy
   p%FF%InvFFZD         = 1.0 / InitInp%dz
   p%FF%InvMFFWS        = 1.0 / p%FF%MeanFFWS
   
   p%FF%TotalTime       = InitInp%nx * InitInp%dx / InitInp%FF%URef
   p%FF%FFYHWid         = 0.5 * InitInp%dy * (InitInp%ny-1)
   p%FF%FFZHWid         = 0.5 * InitInp%dz * (InitInp%nz-1)
   p%FF%GridBase        = InitInp%FF%RefHt - p%FF%FFZHWid
   p%FF%RefHt           = InitInp%FF%RefHt

   p%FF%WindProfileType = InitInp%FF%WindProfileType
   p%FF%Z0              = InitInp%FF%Z0
   p%FF%PLExp           = InitInp%FF%PLExp
   p%FF%AddMeanAfterInterp = .true.

   
   p%FF%InitXPosition   = InitInp%FF%XOffset

   IF ( p%FF%GridBase < 0.0_ReKi ) THEN
      call SetErrStat( ErrID_Severe, 'WARNING: The bottom of the grid is located at a height of '//&
                      TRIM( Num2LStr(p%FF%GridBase) )//' meters, which is below the ground.'//&
                      ' Winds below the ground will be set to 0.', ErrStat,ErrMsg, RoutineName)
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Read data files:
   !-------------------------------------------------------------------------------------------------
   call ReadTurbulenceData( p, InitInp, TmpErrStat, TmpErrMsg)
      CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN
   
   !-------------------------------------------------------------------------------------------------
   ! scale to requested TI (or use requested scale factors)
   !-------------------------------------------------------------------------------------------------
   call ScaleTurbulence(InitInp%FF, p%FF%FFData, InitOut%sf, TmpErrStat, TmpErrMsg)
      CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN      
      

   !-------------------------------------------------------------------------------------------------
   ! Add the mean wind speed to the u component.
   !-------------------------------------------------------------------------------------------------
   if (InitInp%FF%ScaleMethod /= ScaleMethod_None) call SubtractMeanVelocity(p%FF%FFData)
   if (.not. p%FF%AddMeanAfterInterp) call AddMeanVelocity(InitInp%FF, p%FF%GridBase, 1.0_ReKi/p%FF%InvFFZD, p%FF%FFData)
   
   !-------------------------------------------------------------------------------------------------
   ! write info to summary file, if necessary
   !-------------------------------------------------------------------------------------------------
      
   IF ( InitInp%SumFileUnit > 0 ) THEN
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)    'HAWC wind type.  Read by InflowWind sub-module '//TRIM(GetNVD(IfW_HAWCWind_Ver))      
      
      WRITE(InitInp%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Reference height (m):        ',InitInp%FF%RefHt
      WRITE(InitInp%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Timestep (s):                ',p%FF%FFDTime
      WRITE(InitInp%SumFileUnit,'(A34,I12)',  IOSTAT=TmpErrStat)    '     Number of timesteps:         ',p%FF%NFFSteps
      WRITE(InitInp%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Mean windspeed (m/s):        ',p%FF%MeanFFWS
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Time range (s):              [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr( p%FF%TotalTime ))//' ]'
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     X range (m):                 [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr( p%FF%TotalTime * p%FF%MeanFFWS ))//' ]'
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Y range (m):                 [ '// &
                     TRIM(Num2LStr(-p%FF%FFYHWid))//' : '//TRIM(Num2LStr(p%FF%FFYHWid))//' ]'
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Z range (m):                 [ '// &
                     TRIM(Num2LStr(p%FF%GridBase))//' : '//TRIM(Num2LStr(p%FF%GridBase + p%FF%FFZHWid*2.0))//' ]'
      
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)    'Scaling factors used:'
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)    '  u           v           w       '
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)    '----------  ----------  ----------'
      WRITE(InitInp%SumFileUnit,'(F10.3,2x,F10.3,2x,F10.3)',IOSTAT=TmpErrStat)   InitOut%sf      
   ENDIF 
   
   MiscVars%DummyMiscVar = 0
      
   RETURN

END SUBROUTINE IfW_HAWCWind_Init
!====================================================================================================
!>  This routine is used to make sure the initInp data is valid.
SUBROUTINE ValidateInput(InitInp, ErrStat, ErrMsg)

      ! Passed Variables
   TYPE(IfW_HAWCWind_InitInputType),         INTENT(IN   )  :: InitInp           !< Initialization input data passed to the module

      ! Error Handling
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Message about errors
   character(*), parameter                                  :: RoutineName = 'ValidateInput'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! validate the inputs for scaling turbulence:
   CALL FFWind_ValidateInput(InitInp%FF, nc, ErrStat, ErrMsg)   
   
   IF ( InitInp%nx < 1 ) CALL SetErrStat( ErrID_Fatal, 'Number of grid points in the X direction must be at least 1.', ErrStat, ErrMsg, RoutineName )
   IF ( InitInp%ny < 1 ) CALL SetErrStat( ErrID_Fatal, 'Number of grid points in the Y direction must be at least 1.', ErrStat, ErrMsg, RoutineName )
   IF ( InitInp%nz < 1 ) CALL SetErrStat( ErrID_Fatal, 'Number of grid points in the Z direction must be at least 1.', ErrStat, ErrMsg, RoutineName )

   IF ( InitInp%dx < 0.0_ReKi .or. EqualRealNos( InitInp%dx, 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal, 'The grid spacing in the X direction must be larger than 0.', ErrStat, ErrMsg, RoutineName )
   IF ( InitInp%dy < 0.0_ReKi .or. EqualRealNos( InitInp%dy, 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal, 'The grid spacing in the Y direction must be larger than 0.', ErrStat, ErrMsg, RoutineName )
   IF ( InitInp%dz < 0.0_ReKi .or. EqualRealNos( InitInp%dz, 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal, 'The grid spacing in the Z direction must be larger than 0.', ErrStat, ErrMsg, RoutineName )
      
   
END SUBROUTINE ValidateInput
!====================================================================================================
!>  This routine is used read the full-field turbulence data stored in HAWC format.
SUBROUTINE ReadTurbulenceData(p, InitInp, ErrStat, ErrMsg)

      ! Passed Variables
   TYPE(IfW_HAWCWind_ParameterType),         INTENT(INOUT)  :: p                 !< Parameters
   TYPE(IfW_HAWCWind_InitInputType),         INTENT(IN   )  :: InitInp           !< Initialization input data passed to the module
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Message about errors
   
      ! Local Variables:
   INTEGER                                                  :: IC                ! Loop counter for the number of wind components
   INTEGER                                                  :: IX                ! Loop counter for the number of grid points in the X direction
   INTEGER                                                  :: IY                ! Loop counter for the number of grid points in the Y direction
   INTEGER                                                  :: unWind            ! unit number for reading binary files

   INTEGER(IntKi)                                           :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg         ! temporary error message
   CHARACTER(*),                             PARAMETER      :: RoutineName = 'ReadTurbulenceData'

      
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      
   !-------------------------------------------------------------------------------------------------
   ! Allocate space for the wind array.
   !-------------------------------------------------------------------------------------------------

   CALL AllocAry( p%FF%FFData, p%FF%NZGrids,p%FF%NYGrids,p%FF%NFFComp, p%FF%NFFSteps, 'p%FF%FFData', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Read the 3 files containg the turbulent wind speeds.
   !-------------------------------------------------------------------------------------------------
!bjj: check these indices... they do not seem to be very consistant between the WAsP IEC Turbulence
!     simulator and documentation of OC3 file formats... the current implementation is from the
!     OC3/Kenneth Thompson documentation.

      
      ! this could take a while, so we'll write a message indicating what's going on:
         
   CALL WrScr( NewLine//'   Reading HAWC wind files with grids of '//&
      TRIM( Num2LStr(p%FF%NFFSteps) )//' x '//TRIM( Num2LStr(p%FF%NYGrids) )//' x '//TRIM( Num2LStr(p%FF%NZGrids) )//' points'// &
      ' ('//TRIM( Num2LStr(p%FF%FFYHWid*2) )//' m wide, '// TRIM( Num2LStr(p%FF%GridBase) )//' m to '// &
               TRIM( Num2LStr(p%FF%GridBase+p%FF%FFZHWid*2) )//&
               ' m above ground) with a characteristic wind speed of '//TRIM( Num2LStr(p%FF%MeanFFWS) )//' m/s. ' )
            
      
   CALL GetNewUnit( UnWind, TmpErrStat, TmpErrMsg )    
      CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN
                     
      ! The array must be filled so that x(i) < x(i+1), y(i) < y(i+1), and z(i) < z(i+1)
      ! Also, note that the time axis is the negative x axis.      
      
   DO IC = 1,NC

      CALL OpenBInpFile ( UnWind, InitInp%WindFileName(IC), TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
         IF (ErrStat >= AbortErrLev) RETURN

      DO IX = 1,p%FF%NFFSteps
         DO IY = p%FF%NYGrids,1,-1
            !DO IZ = 1,p%FF%NZGrids

               READ( UnWind, IOSTAT=TmpErrStat ) p%FF%FFData(:,iy,ic,ix)  ! note that FFData is SiKi (4-byte reals, not default kinds)

               IF (TmpErrStat /= 0) THEN
                  TmpErrMsg = ' Error reading binary data from "'//TRIM(InitInp%WindFileName(IC))//'". I/O error ' &
                                       //TRIM(Num2LStr(TmpErrStat))//' occurred at IY='//TRIM(Num2LStr(IY))//', IX='//TRIM(Num2LStr(IX))//'.'
                  CLOSE ( UnWind )
                  CALL SetErrStat(ErrID_Fatal, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
                  RETURN
               END IF

            !END DO
         END DO
      END DO

      CLOSE ( UnWind )

   END DO  
   
   
END SUBROUTINE ReadTurbulenceData
!====================================================================================================
!> This routine acts as a wrapper for the GetWindSpeed routine. It steps through the array of input
!! positions and calls the GetWindSpeed routine to calculate the velocities at each point.
!!
!! There are inefficiencies in how this set of routines is coded, but that is a problem for another
!! day. For now, it merely needs to be functional. It can be fixed up and made all pretty later.
!!
!!   16-Apr-2013 - A. Platt, NREL.  Converted to modular framework. Modified for NWTC_Library 2.0
SUBROUTINE IfW_HAWCWind_CalcOutput(Time, PositionXYZ, p, Velocity, DiskVel, MiscVars, ErrStat, ErrMsg)

   IMPLICIT NONE

      ! Passed Variables
   REAL(DbKi),                               INTENT(IN   )  :: Time              !< time from the start of the simulation
   REAL(ReKi),                               INTENT(IN   )  :: PositionXYZ(:,:)  !< Array of XYZ coordinates, 3xN
   TYPE(IfW_HAWCWind_ParameterType),         INTENT(IN   )  :: p                 !< Parameters
   REAL(ReKi),                               INTENT(INOUT)  :: Velocity(:,:)     !< Velocity output at Time    (Set to INOUT so that array does not get deallocated)
   REAL(ReKi),                               INTENT(  OUT)  :: DiskVel(3)        !< HACK for AD14: disk velocity output at Time
   TYPE(IfW_HAWCWind_MiscVarType),           INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)

      ! Error handling
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< The error message


      ! temporary variables
   CHARACTER(*), PARAMETER                                  :: RoutineName = 'IfW_HAWCWind_CalcOutput'


      !-------------------------------------------------------------------------------------------------
      ! Check that the module has been initialized.
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ''

      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------

   CALL IfW_FFWind_CalcOutput(Time, PositionXYZ, p%FF, Velocity, DiskVel, ErrStat, ErrMsg)
   RETURN

END SUBROUTINE IfW_HAWCWind_CalcOutput
!====================================================================================================
!>  This subroutine cleans up any data that is still allocated.  The (possibly) open files are
!!  closed in InflowWindMod.
SUBROUTINE IfW_HAWCWind_End( p, MiscVars, ErrStat, ErrMsg)

   IMPLICIT                                                 NONE

   CHARACTER(*),           PARAMETER                     :: RoutineName="IfW_HAWCWind_End"



      ! Passed Variables
   TYPE(IfW_HAWCWind_ParameterType),      INTENT(INOUT)  :: p                 !< Parameters
   TYPE(IfW_HAWCWind_MiscVarType),        INTENT(  OUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)


      ! Error Handling
   INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg            !< Message about errors


      ! Local Variables
   INTEGER(IntKi)                                        :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                  :: TmpErrMsg      ! temporary error message


      !-=- Initialize the routine -=-

   ErrMsg   = ''
   ErrStat  = ErrID_None




      ! Destroy parameter data

   CALL IfW_HAWCWind_DestroyParam(       p,     TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)


      ! Destroy the state data

   CALL IfW_HAWCWind_DestroyMisc(  MiscVars,   TmpErrStat, TmpErrMsg )
   CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)



END SUBROUTINE IfW_HAWCWind_End

!====================================================================================================
END MODULE IfW_HAWCWind
