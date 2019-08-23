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

   IMPLICIT                                     NONE
   PRIVATE

   TYPE(ProgDesc),   PARAMETER               :: IfW_HAWCWind_Ver = ProgDesc( 'IfW_HAWCWind', '', '' )

   PUBLIC                                    :: IfW_HAWCWind_Init
   PUBLIC                                    :: IfW_HAWCWind_End
   PUBLIC                                    :: IfW_HAWCWind_CalcOutput

   INTEGER(IntKi), PARAMETER  :: nc = 3                           !< number of wind components
   INTEGER(IntKi), PARAMETER  :: WindProfileType_Constant = 0     !< constant wind
   INTEGER(IntKi), PARAMETER  :: WindProfileType_Log      = 1     !< logarithmic
   INTEGER(IntKi), PARAMETER  :: WindProfileType_PL       = 2     !< power law

   INTEGER(IntKi), PARAMETER  :: ScaleMethod_None         = 0     !< no scaling
   INTEGER(IntKi), PARAMETER  :: ScaleMethod_Direct       = 1     !< direct scaling factors
   INTEGER(IntKi), PARAMETER  :: ScaleMethod_StdDev       = 2     !< requested standard deviation
   
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
   p%nx           = InitInp%nx   
   p%ny           = InitInp%ny   
   p%nz           = InitInp%nz   
   p%RefHt        = InitInp%RefHt
   p%URef         = InitInp%URef
   p%InitPosition = 0.0_ReKi  ! bjj: someday we may want to let the users give an offset time/position
   p%InitPosition(1) = InitInp%dx

   p%deltaXInv   = 1.0 / InitInp%dx
   p%deltaYInv   = 1.0 / InitInp%dy
   p%deltaZInv   = 1.0 / InitInp%dz
     
   
   p%LengthX     = InitInp%dx * p%nx !(nx-1)   !because the turbulence box is periodic in the X direction, we need to consider the length between point 1 and the next point 1 (instead of between points 1 and nx)
   p%LengthYHalf = 0.5*InitInp%dy * (p%ny-1)
   p%GridBase    = p%RefHt - 0.5*(p%nz-1)*InitInp%dz

   IF ( p%GridBase < 0.0_ReKi ) THEN
      call SetErrStat( ErrID_Fatal, 'The bottom of the grid is located at a height of '//&
                      TRIM( Num2LStr(p%GridBase) )//' meters, which is below the ground.', ErrStat,ErrMsg, RoutineName)
      RETURN
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
   call ScaleTurbulence( p, InitInp, Interval, InitOut, MiscVars, TmpErrStat, TmpErrMsg)
      CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN      
      

   !-------------------------------------------------------------------------------------------------
   ! Add the mean wind speed to the u component.
   !-------------------------------------------------------------------------------------------------
   call AddMeanVelocity(p, InitInp, TmpErrStat, TmpErrMsg)
      CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN  

   !-------------------------------------------------------------------------------------------------
   ! write info to summary file, if necessary
   !-------------------------------------------------------------------------------------------------
      
   IF ( InitInp%SumFileUnit > 0 ) THEN
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)    'HAWC wind type.  Read by InflowWind sub-module '//TRIM(GetNVD(IfW_HAWCWind_Ver))      
      
      WRITE(InitInp%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Reference height (m):        ',p%RefHt
      WRITE(InitInp%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Timestep (s):                ',p%deltaXInv / p%URef
      WRITE(InitInp%SumFileUnit,'(A34,I12)',  IOSTAT=TmpErrStat)    '     Number of timesteps:         ',p%nx
      WRITE(InitInp%SumFileUnit,'(A34,G12.4)',IOSTAT=TmpErrStat)    '     Mean windspeed (m/s):        ',p%URef
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Time range (s):              [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr( p%LengthX / p%URef ))//' ]'
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     X range (m):                 [ '// &
                     TRIM(Num2LStr(0.0_ReKi))//' : '//TRIM(Num2LStr( p%LengthX ))//' ]'
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Y range (m):                 [ '// &
                     TRIM(Num2LStr(-p%LengthYHalf))//' : '//TRIM(Num2LStr(p%LengthYHalf))//' ]'
      WRITE(InitInp%SumFileUnit,'(A)',        IOSTAT=TmpErrStat)    '     Z range (m):                 [ '// &
                     TRIM(Num2LStr(p%GridBase))//' : '//TRIM(Num2LStr(p%GridBase + p%nz / p%deltaZInv))//' ]'
      
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)    'Scaling factors used:'
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)    '  u           v           w       '
      WRITE(InitInp%SumFileUnit,'(A)', IOSTAT=TmpErrStat)    '----------  ----------  ----------'
      WRITE(InitInp%SumFileUnit,'(F10.3,2x,F10.3,2x,F10.3)',IOSTAT=TmpErrStat)   InitOut%sf      
   ENDIF 
      
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
   
   integer(intki)                                           :: ic                ! loop counter
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF ( InitInp%nx < 1 ) CALL SetErrStat( ErrID_Fatal, 'Number of grid points in the X direction must be at least 1.', ErrStat, ErrMsg, RoutineName )
   IF ( InitInp%ny < 1 ) CALL SetErrStat( ErrID_Fatal, 'Number of grid points in the Y direction must be at least 1.', ErrStat, ErrMsg, RoutineName )
   IF ( InitInp%nz < 1 ) CALL SetErrStat( ErrID_Fatal, 'Number of grid points in the Z direction must be at least 1.', ErrStat, ErrMsg, RoutineName )

   IF ( InitInp%dx < 0.0_ReKi .or. EqualRealNos( InitInp%dx, 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal, 'The grid spacing in the X direction must be larger than 0.', ErrStat, ErrMsg, RoutineName )
   IF ( InitInp%dy < 0.0_ReKi .or. EqualRealNos( InitInp%dy, 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal, 'The grid spacing in the Y direction must be larger than 0.', ErrStat, ErrMsg, RoutineName )
   IF ( InitInp%dz < 0.0_ReKi .or. EqualRealNos( InitInp%dz, 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal, 'The grid spacing in the Z direction must be larger than 0.', ErrStat, ErrMsg, RoutineName )
   
   IF ( InitInp%RefHt < 0.0_ReKi .or. EqualRealNos( InitInp%RefHt, 0.0_ReKi ) ) call SetErrStat( ErrID_Fatal, 'The grid reference height must be larger than 0.', ErrStat, ErrMsg, RoutineName )

   if ( InitInp%ScaleMethod == ScaleMethod_Direct) then
      do ic=1,nc
         if ( InitInp%sf(ic) < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'Turbulence scaling factors must not be negative.', ErrStat, ErrMsg, RoutineName ) 
      end do
   elseif ( InitInp%ScaleMethod == ScaleMethod_StdDev ) then
      do ic=1,nc
         if ( InitInp%sigmaf(ic) < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'Turbulence standard deviations must not be negative.', ErrStat, ErrMsg, RoutineName ) 
      end do
#ifdef UNUSED_INPUTFILE_LINES      
      if ( InitInp%TStart < 0.0_ReKi )  CALL SetErrStat( ErrID_Fatal, 'TStart for turbulence standard deviation calculations must not be negative.', ErrStat, ErrMsg, RoutineName ) 
      if ( InitInp%TEnd <= InitInp%TStart )  CALL SetErrStat( ErrID_Fatal, 'TEnd for turbulence standard deviation calculations must be after TStart.', ErrStat, ErrMsg, RoutineName )       
#endif      
   elseif ( InitInp%ScaleMethod /= ScaleMethod_None ) then
      CALL SetErrStat( ErrID_Fatal, 'HAWC scaling method must be 0 (none), 1 (direct scaling factors), or 2 (target standard deviation).', ErrStat, ErrMsg, RoutineName )             
   end if

   
   if (InitInp%WindProfileType == WindProfileType_Log) then
      if ( InitInp%z0 < 0.0_ReKi .or. EqualRealNos( InitInp%z0, 0.0_ReKi ) ) &
         call SetErrStat( ErrID_Fatal, 'The surface roughness length, Z0, must be greater than zero', ErrStat, ErrMsg, RoutineName )
   elseif ( InitInp%WindProfileType < WindProfileType_Constant .or. InitInp%WindProfileType > WindProfileType_PL)  then                              
       call SetErrStat( ErrID_Fatal, 'The WindProfile type must be 0 (constant), 1 (logarithmic) or 2 (power law).', ErrStat, ErrMsg, RoutineName )
   end if

   IF ( InitInp%URef < 0.0_ReKi ) call SetErrStat( ErrID_Fatal, 'The reference wind speed must not be negative.', ErrStat, ErrMsg, RoutineName )
   
   
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
   INTEGER                                                   :: IC               ! Loop counter for the number of wind components
   INTEGER                                                   :: IX               ! Loop counter for the number of grid points in the X direction
   INTEGER                                                   :: IY               ! Loop counter for the number of grid points in the Y direction
   INTEGER                                                   :: unWind           ! unit number for reading binary files

   INTEGER(IntKi)                                           :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg         ! temporary error message
   CHARACTER(*),                             PARAMETER      :: RoutineName = 'ReadTurbulenceData'

      
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      
   !-------------------------------------------------------------------------------------------------
   ! Allocate space for the wind array.
   !-------------------------------------------------------------------------------------------------

   CALL AllocAry( p%HAWCData, p%nz, p%ny, p%nx, nc, 'p%HAWCData', TmpErrStat, TmpErrMsg )
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
      TRIM( Num2LStr(p%nx) )//' x '//TRIM( Num2LStr(p%ny) )//' x '//TRIM( Num2LStr(p%nz) )//' points.' )
            
      
   CALL GetNewUnit( UnWind, TmpErrStat, TmpErrMsg )    
      CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
      IF (ErrStat >= AbortErrLev) RETURN
                     
      ! The array must be filled so that x(i) < x(i+1), y(i) < y(i+1), and z(i) < z(i+1)
      ! Also, note that the time axis is the negative x axis.      
      
   DO IC = 1,NC

      CALL OpenBInpFile ( UnWind, InitInp%WindFileName(IC), TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName) 
         IF (ErrStat >= AbortErrLev) RETURN

      DO IX = p%nx,1,-1                  ! Time is the opposite of X ....
         DO IY = p%ny,1,-1
            !DO IZ = 1,p%nz

               READ( UnWind, IOSTAT=TmpErrStat ) p%HAWCData(:,iy,ix,ic)  ! note that HAWCData is SiKi (4-byte reals, not default kinds)

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
!>  This routine is used read scale the full-field turbulence data stored in HAWC format.
SUBROUTINE ScaleTurbulence(p, InitInp, Interval, InitOut, MiscVars, ErrStat, ErrMsg)

      ! Passed Variables
   TYPE(IfW_HAWCWind_ParameterType),         INTENT(INOUT)  :: p                 !< Parameters
   TYPE(IfW_HAWCWind_InitInputType),         INTENT(IN   )  :: InitInp           !< Initialization input data passed to the module
   TYPE(IfW_HAWCWind_InitOutputType),        INTENT(INOUT)  :: InitOut           !< Initialization output data passed from the module
   REAL(DbKi),                               INTENT(IN   )  :: Interval          !< Time Interval to use (passed through here)
   TYPE(IfW_HAWCWind_MiscVarType),           INTENT(INOUT)  :: MiscVars          !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Message about errors
   
      ! Local Variables:
   REAL(DbKi)                                               :: v(3)              ! instanteanous wind speed at target position   
   REAL(DbKi)                                               :: vMean(3)          ! average wind speeds over time at target position
   REAL(DbKi)                                               :: vSum(3)           ! sum over time of wind speeds at target position
   REAL(DbKi)                                               :: vSum2(3)          ! sum of wind speeds squared
   REAL(ReKi)                                               :: ActualSigma(3)    ! computed standard deviation
   
   INTEGER                                                  :: ic                ! Loop counter for wind component
   INTEGER                                                  :: ix                ! Loop counter for x
   INTEGER                                                  :: iy                ! Loop counter for y
   INTEGER                                                  :: iz                ! Loop counter for z

   INTEGER(IntKi)                                           :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg         ! temporary error message
   CHARACTER(*),                             PARAMETER      :: RoutineName = 'CalcScaleFactors'

      
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   if ( InitInp%ScaleMethod == ScaleMethod_None ) then
      InitOut%sf = 1.0_ReKi      
   else ! ScaleMethod_Direct or ScaleMethod_StdDev
      if ( InitInp%ScaleMethod == ScaleMethod_Direct ) then
         InitOut%sf = InitInp%sf
      else !if ( InitInp%ScaleMethod == ScaleMethod_StdDev ) then
         
#ifdef UNUSED_INPUTFILE_LINES         
         ! calculate actual standard deviations, then compute scaling factor based on ratio of target to actual
         position(1) = 0.0_ReKi
         position(2) = 0.0_ReKi
         position(3) = p%RefHt
         
         vSum  = 0.0
         vSum2 = 0.0
                       
         n = nint( ( InitInp%TEnd - InitInp%TStart ) / Interval )
         DO it=1,n 
            t = InitInp%TStart + it*Interval
            v = FF_Interp(t,position,p,OtherStates,TmpErrStat,TmpErrMsg)
               CALL SetErrStat( TmpErrStat,TmpErrMsg, ErrStat, ErrMsg, RoutineName )             

            vSum  = vSum  + v
            vSum2 = vSum2 + v**2
         ENDDO ! IT
               
         vMean = vSum/n 
         ActualSigma = SQRT( ABS( (vSum2/n) - vMean**2 ) )
                     
         !InitOut%sf = InitInp%SigmaF / ActualSigma  ! factor = Target / actual
        
         do ic=1,nc
            if ( EqualRealNos( ActualSigma(ic), 0.0_ReKi ) ) then
               InitOut%sf(ic) = 0.0_ReKi
               if ( .not. EqualRealNos( InitInp%SigmaF(ic), 0.0_ReKi ) ) then
                  call SetErrStat( ErrID_Fatal,"Computed standard deviation is zero; cannot scale to achieve target standard deviation.", ErrStat, ErrMsg, RoutineName )                  
               end if         
            else
               InitOut%sf(ic) = InitInp%SigmaF(ic) / ActualSigma(ic)
            end if                           
         end do                  
#else

            ! roughly the point in the center of the grid
         iz = (p%nz + 1) / 2 ! integer division
         iy = (p%ny + 1) / 2 ! integer division
         
         DO ix=1,p%nx 
            v = p%HAWCData(iz,iy,ix,:)
            
            vSum  = vSum  + v
            vSum2 = vSum2 + v**2
         ENDDO ! IT
               
         vMean = vSum/p%nx 
         ActualSigma = SQRT( ABS( (vSum2/p%nx) - vMean**2 ) )
                     
         !InitOut%sf = InitInp%SigmaF / ActualSigma  ! factor = Target / actual        
         do ic=1,nc
            if ( EqualRealNos( ActualSigma(ic), 0.0_ReKi ) ) then
               InitOut%sf(ic) = 0.0_ReKi
               if ( .not. EqualRealNos( InitInp%SigmaF(ic), 0.0_ReKi ) ) then
                  call SetErrStat( ErrID_Fatal,"Computed standard deviation is zero; cannot scale to achieve target standard deviation.", ErrStat, ErrMsg, RoutineName )                  
               end if         
            else
               InitOut%sf(ic) = InitInp%SigmaF(ic) / ActualSigma(ic)
            end if                           
         end do
#endif        
      end if
      
         ! scale using our scaling factors:
      do ic = 1,nc
         p%HAWCData( :, :, :, ic ) = InitOut%sf(IC)*p%HAWCData( :, :, :, ic )   
      end do !IC 
                  
   end if
               
END SUBROUTINE ScaleTurbulence   
!====================================================================================================
!>  This routine is used to add a mean wind profile to the HAWC format turbulence data.
SUBROUTINE AddMeanVelocity(p, InitInp, ErrStat, ErrMsg)

      ! Passed Variables
   TYPE(IfW_HAWCWind_ParameterType),         INTENT(INOUT)  :: p                 !< Parameters
   TYPE(IfW_HAWCWind_InitInputType),         INTENT(IN   )  :: InitInp           !< Initialization input data passed to the module
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< determines if an error has been encountered
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Message about errors
   
      ! Local Variables:
   REAL(ReKi)                                               :: Z                 ! height
   REAL(ReKi)                                               :: U                 ! mean wind speed
   INTEGER(IntKi)                                           :: iz                ! loop counter
   INTEGER(IntKi)                                           :: TmpErrStat        ! temporary error status
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg         ! temporary error message
   CHARACTER(*),                             PARAMETER      :: RoutineName = 'AddMeanVelocity'

      
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   

   DO IZ = 1,p%NZ

      Z = p%GridBase  + ( IZ - 1 )*InitInp%dz

      SELECT CASE ( InitInp%WindProfileType )

      CASE ( WindProfileType_PL )
            
            U = p%URef*( Z / p%RefHt )**InitInp%PLExp      ! [IEC 61400-1 6.3.1.2 (10)]

      CASE ( WindProfileType_Log )

            IF ( .not. EqualRealNos( p%RefHt, InitInp%Z0 ) .and. z > 0.0_ReKi ) THEN
               U = p%URef*( LOG( Z / InitInp%Z0 ) )/( LOG( p%RefHt / InitInp%Z0 ) )
            ELSE
               U = 0.0_ReKi
            ENDIF

      CASE ( WindProfileType_Constant )
         
           U = p%URef            
            
      CASE DEFAULT
         
            U = 0.0_ReKi

      END SELECT

      p%HAWCData( IZ, :, :, 1 ) = p%HAWCData( IZ, :, :, 1 ) + U


   END DO ! IZ
   
               
END SUBROUTINE AddMeanVelocity   
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

      ! local variables
   INTEGER(IntKi)                                           :: NumPoints      ! Number of points specified by the PositionXYZ array

      ! local counters
   INTEGER(IntKi)                                           :: PointNum       ! a loop counter for the current point

      ! temporary variables
   INTEGER(IntKi)                                           :: TmpErrStat     ! temporary error status
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg      ! temporary error message
   CHARACTER(*), PARAMETER                                  :: RoutineName = 'IfW_HAWCWind_CalcOutput'


      !-------------------------------------------------------------------------------------------------
      ! Check that the module has been initialized.
      !-------------------------------------------------------------------------------------------------

   ErrStat     = ErrID_None
   ErrMsg      = ''

      !-------------------------------------------------------------------------------------------------
      ! Initialize some things
      !-------------------------------------------------------------------------------------------------


      ! The array is transposed so that the number of points is the second index, x/y/z is the first.
      ! This is just in case we only have a single point, the SIZE command returns the correct number of points.
   NumPoints   =  SIZE(PositionXYZ,2)

      ! Step through all the positions and get the velocities
   DO PointNum = 1, NumPoints

         ! Calculate the velocity for the position
      Velocity(:,PointNum) = FF_Interp(Time,PositionXYZ(:,PointNum),p,MiscVars,TmpErrStat,TmpErrMsg)


         ! Error handling
      IF (TmpErrStat /= ErrID_None) THEN  !  adding this so we don't have to convert numbers to strings every time
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName//" [position=("//   &
                                                      TRIM(Num2LStr(PositionXYZ(1,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(2,PointNum)))//", "// &
                                                      TRIM(Num2LStr(PositionXYZ(3,PointNum)))//") in wind-file coordinates]" )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF

   ENDDO



      !REMOVE THIS for AeroDyn 15
      ! Return the average disk velocity values needed by AeroDyn 14.  This is the WindInf_ADhack_diskVel routine.
   DiskVel(1)   =  p%URef
   DiskVel(2:3) =  0.0_ReKi


   RETURN

END SUBROUTINE IfW_HAWCWind_CalcOutput
!====================================================================================================
   !+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
   !>    This function is used to interpolate into the full-field wind array.  It receives X, Y, Z and
   !!    TIME from the calling routine.  It then computes a time shift due to a nonzero X based upon
   !!    the average windspeed.  The modified time is used to decide which pair of time slices to interpolate
   !!    within and between.  After finding the two time slices, it decides which four grid points bound the
   !!    (Y,Z) pair.  It does a 3-d linear interpolation.  This routine assumes that X is downwind, Y is to the left when
   !!    looking downwind and Z is up.  It also assumes that no extrapolation will be needed.
   FUNCTION FF_Interp(Time, Position, p, MiscVars, ErrStat, ErrMsg)

      REAL(DbKi),                         INTENT(IN   )  :: Time           !< Time at which to find wind speed
      REAL(ReKi),                         INTENT(IN   )  :: Position(3)    !< takes the place of XGrnd, YGrnd, ZGrnd
      TYPE(IfW_HAWCWind_ParameterType),   INTENT(IN   )  :: p              !< Parameters
      TYPE(IfW_HAWCWind_MiscVarType),     INTENT(INOUT)  :: MiscVars       !< Misc variables for optimization (not copied in glue code)
      REAL(ReKi)                                         :: FF_Interp(3)   !< The U, V, W velocities

      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< error status
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< error message   

         ! Local Variables:
      
      CHARACTER(*),           PARAMETER                  :: RoutineName="FF_Interp"


      REAL(ReKi)                                         :: ShiftedXPosition
      REAL(ReKi),PARAMETER                               :: Tol = 1.0E-3   ! a tolerance for determining if two reals are the same (for extrapolation)
      REAL(ReKi)                                         :: X              ! value between -1 and 1 (where we are between 2 x grid points)
      REAL(ReKi)                                         :: XGRID
      REAL(ReKi)                                         :: Y
      REAL(ReKi)                                         :: YGRID
      REAL(ReKi)                                         :: Z
      REAL(ReKi)                                         :: ZGRID
      REAL(ReKi)                                         :: N(8)           ! array for holding scaling factors for the interpolation algorithm
      REAL(ReKi)                                         :: u(8)           ! array for holding the corner values for the interpolation algorithm across a cubic volume

      INTEGER(IntKi)                                     :: IDIM
      INTEGER(IntKi)                                     :: IXHI
      INTEGER(IntKi)                                     :: IXLO
      INTEGER(IntKi)                                     :: IYHI
      INTEGER(IntKi)                                     :: IYLO
      INTEGER(IntKi)                                     :: IZHI
      INTEGER(IntKi)                                     :: IZLO


      !-------------------------------------------------------------------------------------------------
      ! Initialize variables
      !-------------------------------------------------------------------------------------------------

      FF_Interp(:)         = 0.0_ReKi                         ! the output velocities (in case p%NFFComp /= 3)

      ErrStat              = ErrID_None
      ErrMsg               = ""
      
      !-------------------------------------------------------------------------------------------------
      ! Find the bounding time slices.
      !-------------------------------------------------------------------------------------------------

   ! bjj: should we shift by MIN(YHalfWid,FFZHWid)?

            ! Assume Taylor's Frozen Turbulence Hypothesis applies: u(X,Y,Z,t) = u( X-U*t, Y, Z, 0)

      ShiftedXPosition = Position(1) - TIME*p%URef - p%InitPosition(1)      !this puts the first X grid point at the undeflected tower centerline (or p%InitXPosition)

         ! The wind file is periodic so we'll translate this position to ( 0 <= ShiftedXPosition < p%LengthX )

      ShiftedXPosition = MODULO( ShiftedXPosition, p%LengthX )
      ! If ShiftedXPosition is a very small negative number, modulo returns the incorrect value due to internal rounding errors.
      ! See bug report #471
      IF (ShiftedXPosition == p%LengthX) ShiftedXPosition = 0.0_ReKi

      XGrid = ShiftedXPosition * p%deltaXInv

      IXLO = INT( XGrid )              ! convert REAL to INTEGER
      X = 2.0_ReKi * ( XGrid - REAL(IXLO, ReKi) ) - 1.0_ReKi     ! a value between -1 and 1 that indicates a relative position between IXLO and IXHI
      IXLO = IXLO + 1  ! Add 1 because our grids start at 1, not zero
      
      IF ( IXLO == p%NX ) THEN
         IXHI = 1
      ELSE
         IXHI = IXLO + 1
      ENDIF
      
      !-------------------------------------------------------------------------------------------------
      ! Find the bounding rows for the Z position. [The lower-left corner is (1,1) when looking upwind.]
      !-------------------------------------------------------------------------------------------------

      ZGRID = ( Position(3) - p%GridBase - p%InitPosition(3) )*p%deltaZInv

         ! Index for start and end slices
      IZLO = INT( ZGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
      IZHI = IZLO + 1

         ! Set Z as a value between -1 and 1 for the relative location between IZLO and IZHI.
         ! Subtract 1_IntKi from Z since the indices are starting at 1, not 0
      Z = 2.0_ReKi * (ZGRID - REAL(IZLO - 1_IntKi, ReKi)) - 1.0_ReKi

      IF ( IZLO < 1 ) THEN
         IF ( IZLO == 0 .AND. Z >= 1.0-TOL ) THEN
            Z    = -1.0_ReKi
            IZLO = 1
         ELSE
            ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                        TRIM(Num2LStr(Position(3)))//' m is below the grid).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF
      ELSEIF ( IZLO >= p%nz ) THEN
         IF ( IZLO == p%nz .AND. Z <= TOL ) THEN
            Z    = -1.0_ReKi
            IZHI = IZLO                   ! We're right on the last point, which is still okay
         ELSE
            ErrMsg   = ' FF wind array boundaries violated. Grid too small in Z direction (Z='//&
                        TRIM(Num2LStr(Position(3)))//' m is above the grid).'
            ErrStat  = ErrID_Fatal
            RETURN
         ENDIF
      ENDIF


      !-------------------------------------------------------------------------------------------------
      ! Find the bounding columns for the Y position. [The lower-left corner is (1,1) when looking upwind.]
      !-------------------------------------------------------------------------------------------------
      YGRID = ( Position(2) + p%LengthYHalf - p%InitPosition(2) )*p%deltaYInv    ! really, it's (Position(2) - -1.0*p%LengthYHalf)

      IYLO = INT( YGRID ) + 1             ! convert REAL to INTEGER, then add one since our grids start at 1, not 0
      IYHI = IYLO + 1

         ! Set Y as a value between -1 and 1 for the relative location between IYLO and IYHI.  Used in the interpolation.
         ! Subtract 1_IntKi from IYLO since grids start at index 1, not 0
      Y = 2.0_ReKi * (YGRID - REAL(IYLO - 1_IntKi, ReKi)) - 1.0_ReKi

      IF ( IYLO >= p%ny .OR. IYLO < 1 ) THEN
         IF ( IYLO == 0 .AND. Y >= 1.0-TOL ) THEN
            Y    = -1.0_ReKi
            IYLO = 1
         ELSE IF ( IYLO == p%ny .AND. Y <= TOL ) THEN
            Y    = -1.0_ReKi
            IYHI = IYLO                   ! We're right on the last point, which is still okay
         ELSE
            ErrMsg   = ' FF wind array boundaries violated: Grid too small in Y direction. Y='// &
                        TRIM(Num2LStr(Position(2)))//'; Y boundaries = ['//TRIM(Num2LStr(-1.0*p%LengthYHalf))// &
                        ', '//TRIM(Num2LStr(p%LengthYHalf))//']'
            ErrStat = ErrID_Fatal         ! we don't return anything
            RETURN
         ENDIF
      ENDIF

      !-------------------------------------------------------------------------------------------------
      ! Interpolate on the grid
      !-------------------------------------------------------------------------------------------------

      DO IDIM=1,nc       ! all the components


         N(1)  = ( 1.0_ReKi + Z )*( 1.0_ReKi - Y )*( 1.0_ReKi - X )
         N(2)  = ( 1.0_ReKi + Z )*( 1.0_ReKi + Y )*( 1.0_ReKi - X )
         N(3)  = ( 1.0_ReKi - Z )*( 1.0_ReKi + Y )*( 1.0_ReKi - X )
         N(4)  = ( 1.0_ReKi - Z )*( 1.0_ReKi - Y )*( 1.0_ReKi - X )
         N(5)  = ( 1.0_ReKi + Z )*( 1.0_ReKi - Y )*( 1.0_ReKi + X )
         N(6)  = ( 1.0_ReKi + Z )*( 1.0_ReKi + Y )*( 1.0_ReKi + X )
         N(7)  = ( 1.0_ReKi - Z )*( 1.0_ReKi + Y )*( 1.0_ReKi + X )
         N(8)  = ( 1.0_ReKi - Z )*( 1.0_ReKi - Y )*( 1.0_ReKi + X )
         N     = N / REAL( SIZE(N), ReKi )  ! normalize


         u(1)  = p%HAWCData( IZHI, IYLO, IXLO, IDIM )
         u(2)  = p%HAWCData( IZHI, IYHI, IXLO, IDIM )
         u(3)  = p%HAWCData( IZLO, IYHI, IXLO, IDIM )
         u(4)  = p%HAWCData( IZLO, IYLO, IXLO, IDIM )
         u(5)  = p%HAWCData( IZHI, IYLO, IXHI, IDIM )
         u(6)  = p%HAWCData( IZHI, IYHI, IXHI, IDIM )
         u(7)  = p%HAWCData( IZLO, IYHI, IXHI, IDIM )
         u(8)  = p%HAWCData( IZLO, IYLO, IXHI, IDIM )
            
         FF_Interp(IDIM)  =  SUM ( N * u ) 


      END DO !IDIM

      RETURN

   END FUNCTION FF_Interp
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
