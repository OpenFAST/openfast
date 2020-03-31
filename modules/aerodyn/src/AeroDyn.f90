!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2018  Envision Energy USA, LTD
!
!    This file is part of AeroDyn.
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
!> AeroDyn is a time-domain aerodynamics module for horizontal-axis wind turbines.
module AeroDyn
    
   use NWTC_Library
   use AeroDyn_Types
   use AeroDyn_IO
   use BEMT
   use AirfoilInfo
   use NWTC_LAPACK
   use UnsteadyAero
   
   
   implicit none

   private
         

   ! ..... Public Subroutines ...................................................................................................

   public :: AD_Init                           ! Initialization routine
   public :: AD_ReInit                         ! Routine to reinitialize driver (re-initializes the states)
   public :: AD_End                            ! Ending routine (includes clean up)
   public :: AD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   public :: AD_CalcOutput                     ! Routine for computing outputs
   public :: AD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   
   
   PUBLIC :: AD_JacobianPInput                 ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the inputs(u)
   PUBLIC :: AD_JacobianPContState             ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the continuous
                                               !   states(x)
   PUBLIC :: AD_JacobianPDiscState             ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the discrete
                                               !   states(xd)
   PUBLIC :: AD_JacobianPConstrState           ! Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                               !   (Xd), and constraint - state(Z) functions all with respect to the constraint
                                               !   states(z)
   PUBLIC :: AD_GetOP                          !< Routine to pack the operating point values (for linearization) into arrays
   
  
contains    
!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine sets the initialization output data structure, which contains data to be returned to the calling program (e.g.,
!! FAST or AeroDyn_Driver)   
subroutine AD_SetInitOut(p, InputFileData, InitOut, errStat, errMsg)

   type(AD_InitOutputType),       intent(  out)  :: InitOut          ! output data
   type(AD_InputFile),            intent(in   )  :: InputFileData    ! input file data (for setting airfoil shape outputs)
   type(AD_ParameterType),        intent(in   )  :: p                ! Parameters
   integer(IntKi),                intent(  out)  :: errStat          ! Error status of the operation
   character(*),                  intent(  out)  :: errMsg           ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'AD_SetInitOut'
   
   
   
   integer(IntKi)                               :: i, j, k, f
   integer(IntKi)                               :: NumCoords
#ifdef DBG_OUTS
   integer(IntKi)                               :: m
   character(6)                                 ::chanPrefix
   character(3)                                 :: TmpChar
#endif   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
   InitOut%AirDens = p%AirDens
   
   call AllocAry( InitOut%WriteOutputHdr, p%numOuts, 'WriteOutputHdr', errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   
   call AllocAry( InitOut%WriteOutputUnt, p%numOuts, 'WriteOutputUnt', errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

   if (ErrStat >= AbortErrLev) return
      
   
#ifdef DBG_OUTS
   ! Loop over blades and nodes to populate the output channel names and units
   
   do k=1,p%numBlades
      do j=1, p%NumBlNds
         
         m = (k-1)*p%NumBlNds*23 + (j-1)*23 
         
         WRITE (TmpChar,'(I3.3)') j
         chanPrefix = "B"//trim(num2lstr(k))//"N"//TmpChar
         InitOut%WriteOutputHdr( m + 1 ) = trim(chanPrefix)//"Twst"
         InitOut%WriteOutputUnt( m + 1 ) = '  (deg)  '
         InitOut%WriteOutputHdr( m + 2 ) = trim(chanPrefix)//"Psi"
         InitOut%WriteOutputUnt( m + 2 ) = '  (deg)  '
         InitOut%WriteOutputHdr( m + 3 ) = trim(chanPrefix)//"Vx"
         InitOut%WriteOutputUnt( m + 3 ) = '  (m/s)  '
         InitOut%WriteOutputHdr( m + 4 ) = trim(chanPrefix)//"Vy"
         InitOut%WriteOutputUnt( m + 4 ) = '  (m/s)  '
         InitOut%WriteOutputHdr( m + 5 ) = ' '//trim(chanPrefix)//"AIn"
         InitOut%WriteOutputUnt( m + 5 ) = '  (-)  '
         InitOut%WriteOutputHdr( m + 6 ) = ' '//trim(chanPrefix)//"ApIn"
         InitOut%WriteOutputUnt( m + 6 ) = '  (-)  '
         InitOut%WriteOutputHdr( m + 7 ) = trim(chanPrefix)//"Vrel"
         InitOut%WriteOutputUnt( m + 7 ) = '  (m/s)  '
         InitOut%WriteOutputHdr( m + 8 ) = ' '//trim(chanPrefix)//"Phi"
         InitOut%WriteOutputUnt( m + 8 ) = '  (deg)  '
         InitOut%WriteOutputHdr( m + 9 ) = ' '//trim(chanPrefix)//"AOA"
         InitOut%WriteOutputUnt( m + 9 ) = '  (deg)  '
         InitOut%WriteOutputHdr( m + 10 ) = ' '//trim(chanPrefix)//"Cl"
         InitOut%WriteOutputUnt( m + 10 ) = '   (-)   '
         InitOut%WriteOutputHdr( m + 11 ) = ' '//trim(chanPrefix)//"Cd"
         InitOut%WriteOutputUnt( m + 11 ) = '   (-)   '
         InitOut%WriteOutputHdr( m + 12 ) = ' '//trim(chanPrefix)//"Cm"
         InitOut%WriteOutputUnt( m + 12 ) = '   (-)   '
         InitOut%WriteOutputHdr( m + 13 ) = ' '//trim(chanPrefix)//"Cx"
         InitOut%WriteOutputUnt( m + 13 ) = '   (-)   '
         InitOut%WriteOutputHdr( m + 14 ) = ' '//trim(chanPrefix)//"Cy"
         InitOut%WriteOutputUnt( m + 14 ) = '   (-)   '
         InitOut%WriteOutputHdr( m + 15 ) = ' '//trim(chanPrefix)//"Cn"
         InitOut%WriteOutputUnt( m + 15 ) = '   (-)   '
         InitOut%WriteOutputHdr( m + 16 ) = ' '//trim(chanPrefix)//"Ct"
         InitOut%WriteOutputUnt( m + 16 ) = '   (-)   '
         InitOut%WriteOutputHdr( m + 17 ) = ' '//trim(chanPrefix)//"Fl"
         InitOut%WriteOutputUnt( m + 17 ) = '  (N/m)  '
         InitOut%WriteOutputHdr( m + 18 ) = ' '//trim(chanPrefix)//"Fd"
         InitOut%WriteOutputUnt( m + 18 ) = '  (N/m)  '
         InitOut%WriteOutputHdr( m + 19 ) = ' '//trim(chanPrefix)//"M"
         InitOut%WriteOutputUnt( m + 19 ) = ' (N/m^2) '
         InitOut%WriteOutputHdr( m + 20 ) = ' '//trim(chanPrefix)//"Fx"
         InitOut%WriteOutputUnt( m + 20 ) = '  (N/m)  '
         InitOut%WriteOutputHdr( m + 21 ) = ' '//trim(chanPrefix)//"Fy"
         InitOut%WriteOutputUnt( m + 21 ) = '  (N/m)  '
         InitOut%WriteOutputHdr( m + 22 ) = ' '//trim(chanPrefix)//"Fn"
         InitOut%WriteOutputUnt( m + 22 ) = '  (N/m)  '
         InitOut%WriteOutputHdr( m + 23 ) = ' '//trim(chanPrefix)//"Ft"
         InitOut%WriteOutputUnt( m + 23 ) = '  (N/m)  '
         
      end do
   end do
#else
   do i=1,p%NumOuts
      InitOut%WriteOutputHdr(i) = p%OutParam(i)%Name
      InitOut%WriteOutputUnt(i) = p%OutParam(i)%Units
   end do
#endif
                      
   
   InitOut%Ver = AD_Ver
   
! set visualization data:
      ! this check is overly restrictive, but it would be a lot of work to ensure that only the *used* airfoil 
      ! tables have the same number of coordinates.
   if ( allocated(p%AFI) ) then  
      
      if ( p%AFI(1)%NumCoords > 0 ) then
         NumCoords = p%AFI(1)%NumCoords
         do i=2,size(p%AFI)
            if (p%AFI(i)%NumCoords /= NumCoords) then
               call SetErrStat( ErrID_Info, 'Airfoil files do not contain the same number of x-y coordinates.', ErrStat, ErrMsg, RoutineName )
               NumCoords = -1
               exit
            end if            
         end do
            
         if (NumCoords > 0) then
            if (NumCoords < 3) then
               call SetErrStat( ErrID_Info, 'Airfoil files with NumCoords > 0 must contain at least 2 coordinates.', ErrStat, ErrMsg, RoutineName )
               return
            end if     
               
            allocate( InitOut%BladeShape( p%numBlades ), STAT=ErrStat2 )
            if (ErrStat2 /= 0) then
               call SetErrStat( ErrID_Info, 'Error allocationg InitOut%AD_BladeShape', ErrStat, ErrMsg, RoutineName )
               return
            end if     
            
            do k=1,p%numBlades
               call allocAry(  InitOut%BladeShape(k)%AirfoilCoords, 2, NumCoords-1, InputFileData%BladeProps(k)%NumBlNds, 'AirfoilCoords', ErrStat2, ErrMsg2)
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  if (ErrStat >= AbortErrLev) return
                  
               do j=1,InputFileData%BladeProps(k)%NumBlNds
                  f = InputFileData%BladeProps(k)%BlAFID(j)
                  
                  do i=1,NumCoords-1                                                     
                     InitOut%BladeShape(k)%AirfoilCoords(1,i,j) = InputFileData%BladeProps(k)%BlChord(j)*( p%AFI(f)%Y_Coord(i+1) - p%AFI(f)%Y_Coord(1) )
                     InitOut%BladeShape(k)%AirfoilCoords(2,i,j) = InputFileData%BladeProps(k)%BlChord(j)*( p%AFI(f)%X_Coord(i+1) - p%AFI(f)%X_Coord(1) )
                  end do                  
               end do
                                 
            end do
            
         end if                  
      end if
      
   end if
   
   
   ! set blade properties data  ! bjj: I would probably do a move_alloc() at the end of the init routine rather than make a copy like this.... 
   ALLOCATE(InitOut%BladeProps(p%numBlades), STAT = ErrStat2)
   IF (ErrStat2 /= 0) THEN
      CALL SetErrStat(ErrID_Fatal,"Error allocating memory for BladeProps.", ErrStat, ErrMsg, RoutineName)
      RETURN
   END IF
   do k=1,p%numBlades
      ! allocate space and copy blade data:
      CALL AD_CopyBladePropsType(InputFileData%BladeProps(k), InitOut%BladeProps(k), MESH_NEWCOPY, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end do

   !Tower data
   IF ( p%NumTwrNds > 0 ) THEN
      ALLOCATE(InitOut%TwrElev(p%NumTwrNds), STAT = ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating memory for TwrElev.", ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF
      InitOut%TwrElev(:) = InputFileData%TwrElev(:)

      ALLOCATE(InitOut%TwrDiam(p%NumTwrNds), STAT = ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,"Error allocating memory for TwrDiam.", ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF   
      InitOut%TwrDiam(:) = p%TwrDiam(:)
   END IF  
   
end subroutine AD_SetInitOut

!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine AD_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   type(AD_InitInputType),       intent(in   ) :: InitInp       !< Input data for initialization routine
   type(AD_InputType),           intent(  out) :: u             !< An initial guess for the input; input mesh must be defined
   type(AD_ParameterType),       intent(  out) :: p             !< Parameters
   type(AD_ContinuousStateType), intent(  out) :: x             !< Initial continuous states
   type(AD_DiscreteStateType),   intent(  out) :: xd            !< Initial discrete states
   type(AD_ConstraintStateType), intent(  out) :: z             !< Initial guess of the constraint states
   type(AD_OtherStateType),      intent(  out) :: OtherState    !< Initial other states
   type(AD_OutputType),          intent(  out) :: y             !< Initial system outputs (outputs are not calculated;
                                                                !!   only the output mesh is initialized)
   type(AD_MiscVarType),         intent(  out) :: m             !< Initial misc/optimization variables
   real(DbKi),                   intent(inout) :: interval      !< Coupling interval in seconds: the rate that
                                                                !!   (1) AD_UpdateStates() is called in loose coupling &
                                                                !!   (2) AD_UpdateDiscState() is called in tight coupling.
                                                                !!   Input is the suggested time from the glue code;
                                                                !!   Output is the actual coupling interval that will be used
                                                                !!   by the glue code.
   type(AD_InitOutputType),      intent(  out) :: InitOut       !< Output for initialization routine
   integer(IntKi),               intent(  out) :: errStat       !< Error status of the operation
   character(*),                 intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None
   

      ! Local variables
   integer(IntKi)                              :: i             ! loop counter
   
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
      
   type(AD_InputFile)                          :: InputFileData ! Data stored in the module's input file
   integer(IntKi)                              :: UnEcho        ! Unit number for the echo file
   
   character(*), parameter                     :: RoutineName = 'AD_Init'
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   UnEcho  = -1

      ! Initialize the NWTC Subroutine Library

   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   call DispNVD( AD_Ver )
   
   
   p%NumBlades = InitInp%NumBlades ! need this before reading the AD input file so that we know how many blade files to read
   !bjj: note that we haven't validated p%NumBlades before using it below!
   p%RootName  = TRIM(InitInp%RootName)//'.AD'
   
      ! Read the primary AeroDyn input file
   call ReadInputFiles( InitInp%InputFile, InputFileData, interval, p%RootName, p%NumBlades, UnEcho, ErrStat2, ErrMsg2 )   
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
         
      
      ! Validate the inputs
   call ValidateInputData( InitInp, InputFileData, p%NumBlades, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
      
      !............................................................................................
      ! Define parameters
      !............................................................................................
      
      ! Initialize AFI module (read Airfoil tables)
   call Init_AFIparams( InputFileData, p%AFI, UnEcho, p%NumBlades, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
         
      
      ! set the rest of the parameters
   call SetParameters( InitInp, InputFileData, p, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
   
      !............................................................................................
      ! Define and initialize inputs here 
      !............................................................................................
   
   call Init_u( u, p, InputFileData, InitInp, errStat2, errMsg2 ) 
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

      ! 

      !............................................................................................
      ! Initialize the BEMT module (also sets other variables for sub module)
      !............................................................................................
      
      ! initialize BEMT after setting parameters and inputs because we are going to use the already-
      ! calculated node positions from the input meshes
      
   call Init_BEMTmodule( InputFileData, u, m%BEMT_u(1), p, x%BEMT, xd%BEMT, z%BEMT, &
                           OtherState%BEMT, m%BEMT_y, m%BEMT, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
         
   call BEMT_CopyInput( m%BEMT_u(1), m%BEMT_u(2), MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
         
      
      !............................................................................................
      ! Define outputs here
      !............................................................................................
   call Init_y(y, u, p, errStat2, errMsg2) ! do this after input meshes have been initialized
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
   
   
      !............................................................................................
      ! Initialize states and misc vars
      !............................................................................................
      
      ! many states are in the BEMT module, which were initialized in BEMT_Init()
      
   call Init_MiscVars(m, p, u, y, errStat2, errMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      
      !............................................................................................
      ! Define initialization output here
      !............................................................................................
   call AD_SetInitOut(p, InputFileData, InitOut, errStat2, errMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   
      ! after setting InitOut variables, we really don't need the airfoil coordinates taking up
      ! space in AeroDyn
   if ( allocated(p%AFI) ) then  
      do i=1,size(p%AFI)
         if (allocated(p%AFI(i)%X_Coord)) deallocate( p%AFI(i)%X_Coord) 
         if (allocated(p%AFI(i)%Y_Coord)) deallocate( p%AFI(i)%Y_Coord) 
      end do
   end if
   
      !............................................................................................
      ! Initialize Jacobian:
      !............................................................................................
   if (InitInp%Linearize) then      
      call Init_Jacobian(InputFileData, p, u, y, m, InitOut, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if
   
      !............................................................................................
      ! Print the summary file if requested:
      !............................................................................................
   if (InputFileData%SumPrint) then
      call AD_PrintSum( InputFileData, p, u, y, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   end if
      
            
   call Cleanup() 
      
contains
   subroutine Cleanup()

      CALL AD_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2 )
      IF ( UnEcho > 0 ) CLOSE( UnEcho )
      
   end subroutine Cleanup

end subroutine AD_Init
!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine reinitializes BEMT and UA, assuming that we will start the simulation over again, with only the inputs being different.
!! This allows us to bypass reading input files and allocating arrays because p is already set.
subroutine AD_ReInit(p, x, xd, z, OtherState, m, Interval, ErrStat, ErrMsg )   

   type(AD_ParameterType),       intent(in   ) :: p             !< Parameters
   type(AD_ContinuousStateType), intent(inout) :: x             !< Initial continuous states
   type(AD_DiscreteStateType),   intent(inout) :: xd            !< Initial discrete states
   type(AD_ConstraintStateType), intent(inout) :: z             !< Initial guess of the constraint states
   type(AD_OtherStateType),      intent(inout) :: OtherState    !< Initial other states
   type(AD_MiscVarType),         intent(inout) :: m             !< Initial misc/optimization variables
   real(DbKi),                   intent(in   ) :: interval      !< Coupling interval in seconds: the rate that
                                                                !!   (1) AD_UpdateStates() is called in loose coupling &
                                                                !!   (2) AD_UpdateDiscState() is called in tight coupling.
                                                                !!   Input is the suggested time from the glue code;
                                                                !!   Output is the actual coupling interval that will be used
                                                                !!   by the glue code.
   integer(IntKi),               intent(  out) :: errStat       !< Error status of the operation
   character(*),                 intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None

   character(*), parameter                     :: RoutineName = 'AD_ReInit'

   
   ErrStat = ErrID_None
   ErrMsg = ''
   
   if ( .not. EqualRealNos(p%DT, interval) ) then
      call SetErrStat( ErrID_Fatal, 'When AD is reinitialized, DT must not change.', ErrStat, ErrMsg, RoutineName )
      ! we could get around this by figuring out what needs to change when we modify the dt parameter... probably just some unused-parameters
      ! and the UA filter
   end if
      
   call BEMT_ReInit(p%BEMT,x%BEMT,xd%BEMT,z%BEMT,OtherState%BEMT,m%BEMT,p%AFI)
      
end subroutine AD_ReInit
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes (allocates) the misc variables for use during the simulation.
subroutine Init_MiscVars(m, p, u, y, errStat, errMsg)
   type(AD_MiscVarType),          intent(inout)  :: m                !< misc/optimization data (not defined in submodules)
   type(AD_ParameterType),        intent(in   )  :: p                !< Parameters
   type(AD_InputType),            intent(inout)  :: u                !< input for HubMotion mesh (create sibling mesh here)
   type(AD_OutputType),           intent(in   )  :: y                !< output (create mapping between output and otherstate mesh here)
   integer(IntKi),                intent(  out)  :: errStat          !< Error status of the operation
   character(*),                  intent(  out)  :: errMsg           !< Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: k
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_OtherStates'

      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
   call AllocAry( m%DisturbedInflow, 3_IntKi, p%NumBlNds, p%numBlades, 'OtherState%DisturbedInflow', ErrStat2, ErrMsg2 ) ! must be same size as u%InflowOnBlade
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( m%WithoutSweepPitchTwist, 3_IntKi, 3_IntKi, p%NumBlNds, p%numBlades, 'OtherState%WithoutSweepPitchTwist', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   
   call allocAry( m%SigmaCavit, p%NumBlNds, p%numBlades, 'm%SigmaCavit', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( m%SigmaCavitCrit, p%NumBlNds, p%numBlades, 'm%SigmaCavitCrit', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call allocAry( m%CavitWarnSet, p%NumBlNds, p%numBlades, 'm%CavitWarnSet', errStat2, errMsg2); call setErrStat(errStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   m%SigmaCavit     = 0.0_ReKi      !Init to zero for output files in case a cavit check isnt done but output is requested 
   m%SigmaCavitCrit = 0.0_ReKi
   m%CavitWarnSet   = .false.
         ! arrays for output
#ifdef DBG_OUTS
   allocate( m%AllOuts(0:p%NumOuts), STAT=ErrStat2 ) ! allocate starting at zero to account for invalid output channels
#else
   allocate( m%AllOuts(0:MaxOutPts), STAT=ErrStat2 ) ! allocate starting at zero to account for invalid output channels
#endif
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Fatal, "Error allocating AllOuts.", errStat, errMsg, RoutineName )
         return
      end if
   m%AllOuts = 0.0_ReKi
   
      ! save these tower calculations for output:
   call AllocAry( m%W_Twr, p%NumTwrNds, 'm%W_Twr', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( m%X_Twr, p%NumTwrNds, 'm%X_Twr', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( m%Y_Twr, p%NumTwrNds, 'm%Y_Twr', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      ! save blade calculations for output:
if (p%TwrPotent /= TwrPotent_none .or. p%TwrShadow) then
   call AllocAry( m%TwrClrnc, p%NumBlNds, p%NumBlades, 'm%TwrClrnc', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
end if            
   call AllocAry( m%Curve, p%NumBlNds, p%NumBlades, 'm%Curve', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )            
   call AllocAry( m%X, p%NumBlNds, p%NumBlades, 'm%X', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( m%Y, p%NumBlNds, p%NumBlades, 'm%Y', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( m%M, p%NumBlNds, p%NumBlades, 'm%M', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      ! mesh mapping data for integrating load over entire rotor:
   allocate( m%B_L_2_H_P(p%NumBlades), Stat = ErrStat2)
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Fatal, "Error allocating B_L_2_H_P mapping structure.", errStat, errMsg, RoutineName )
         return
      end if

   call MeshCopy (  SrcMesh  = u%HubMotion        &
                  , DestMesh = m%HubLoad          &
                  , CtrlCode = MESH_SIBLING       &
                  , IOS      = COMPONENT_OUTPUT   &
                  , force    = .TRUE.             &
                  , moment   = .TRUE.             &
                  , ErrStat  = ErrStat2           &
                  , ErrMess  = ErrMsg2            )
   
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) RETURN         
   
   do k=1,p%NumBlades
      CALL MeshMapCreate( y%BladeLoad(k), m%HubLoad, m%B_L_2_H_P(k), ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName//':B_L_2_H_P('//TRIM(Num2LStr(K))//')' )
   end do
   
   if (ErrStat >= AbortErrLev) RETURN

   ! 
   if (p%NumTwrNds > 0) then
      m%W_Twr = 0.0_ReKi
      m%X_Twr = 0.0_ReKi
      m%Y_Twr = 0.0_ReKi
   end if
   
   
   
end subroutine Init_MiscVars
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes AeroDyn meshes and output array variables for use during the simulation.
subroutine Init_y(y, u, p, errStat, errMsg)
   type(AD_OutputType),           intent(  out)  :: y               !< Module outputs
   type(AD_InputType),            intent(inout)  :: u               !< Module inputs -- intent(out) because of mesh sibling copy
   type(AD_ParameterType),        intent(in   )  :: p               !< Parameters
   integer(IntKi),                intent(  out)  :: errStat         !< Error status of the operation
   character(*),                  intent(  out)  :: errMsg          !< Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: k                 ! loop counter for blades
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_y'

      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
         
   if (p%TwrAero) then
            
      call MeshCopy ( SrcMesh  = u%TowerMotion    &
                    , DestMesh = y%TowerLoad      &
                    , CtrlCode = MESH_SIBLING     &
                    , IOS      = COMPONENT_OUTPUT &
                    , force    = .TRUE.           &
                    , moment   = .TRUE.           &
                    , ErrStat  = ErrStat2         &
                    , ErrMess  = ErrMsg2          )
   
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
         if (ErrStat >= AbortErrLev) RETURN         
         
         !y%TowerLoad%force = 0.0_ReKi  ! shouldn't have to initialize this
         !y%TowerLoad%moment= 0.0_ReKi  ! shouldn't have to initialize this
   else
      y%TowerLoad%nnodes = 0
   end if

   
   allocate( y%BladeLoad(p%numBlades), stat=ErrStat2 )
   if (errStat2 /= 0) then
      call SetErrStat( ErrID_Fatal, 'Error allocating y%BladeLoad.', ErrStat, ErrMsg, RoutineName )      
      return
   end if
   

   do k = 1, p%numBlades
   
      call MeshCopy ( SrcMesh  = u%BladeMotion(k) &
                    , DestMesh = y%BladeLoad(k)   &
                    , CtrlCode = MESH_SIBLING     &
                    , IOS      = COMPONENT_OUTPUT &
                    , force    = .TRUE.           &
                    , moment   = .TRUE.           &
                    , ErrStat  = ErrStat2         &
                    , ErrMess  = ErrMsg2          )
   
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
                           
   end do
   
   call AllocAry( y%WriteOutput, p%numOuts, 'WriteOutput', errStat2, errMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) RETURN      
   
   
   
end subroutine Init_y
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes AeroDyn meshes and input array variables for use during the simulation.
subroutine Init_u( u, p, InputFileData, InitInp, errStat, errMsg )
!..................................................................................................................................

   type(AD_InputType),           intent(  out)  :: u                 !< Input data
   type(AD_ParameterType),       intent(in   )  :: p                 !< Parameters
   type(AD_InputFile),           intent(in   )  :: InputFileData     !< Data stored in the module's input file
   type(AD_InitInputType),       intent(in   )  :: InitInp           !< Input data for AD initialization routine
   integer(IntKi),               intent(  out)  :: errStat           !< Error status of the operation
   character(*),                 intent(  out)  :: errMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   real(reKi)                                   :: position(3)       ! node reference position
   real(reKi)                                   :: positionL(3)      ! node local position
   real(R8Ki)                                   :: theta(3)          ! Euler angles
   real(R8Ki)                                   :: orientation(3,3)  ! node reference orientation
   real(R8Ki)                                   :: orientationL(3,3) ! node local orientation
   
   integer(intKi)                               :: j                 ! counter for nodes
   integer(intKi)                               :: k                 ! counter for blades
   
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_u'

      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Arrays for InflowWind inputs:
   
   call AllocAry( u%InflowOnBlade, 3_IntKi, p%NumBlNds, p%numBlades, 'u%InflowOnBlade', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( u%InflowOnTower, 3_IntKi, p%NumTwrNds, 'u%InflowOnTower', ErrStat2, ErrMsg2 ) ! could be size zero
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

   call AllocAry( u%UserProp, p%NumBlNds, p%numBlades, 'u%UserProp', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      
   if (errStat >= AbortErrLev) return      
      
   u%InflowOnBlade = 0.0_ReKi
   u%UserProp      = 0.0_ReKi
   
      ! Meshes for motion inputs (ElastoDyn and/or BeamDyn)
         !................
         ! tower
         !................
   if (p%NumTwrNds > 0) then
      
      u%InflowOnTower = 0.0_ReKi 
      
      call MeshCreate ( BlankMesh = u%TowerMotion   &
                       ,IOS       = COMPONENT_INPUT &
                       ,Nnodes    = p%NumTwrNds     &
                       ,ErrStat   = ErrStat2        &
                       ,ErrMess   = ErrMsg2         &
                       ,Orientation     = .true.    &
                       ,TranslationDisp = .true.    &
                       ,TranslationVel  = .true.    &
                      )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      if (errStat >= AbortErrLev) return
            
         ! set node initial position/orientation
      position = 0.0_ReKi
      do j=1,p%NumTwrNds         
         position(3) = InputFileData%TwrElev(j)
         
         call MeshPositionNode(u%TowerMotion, j, position, errStat2, errMsg2)  ! orientation is identity by default
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end do !j
         
         ! create line2 elements
      do j=1,p%NumTwrNds-1
         call MeshConstructElement( u%TowerMotion, ELEMENT_LINE2, errStat2, errMsg2, p1=j, p2=j+1 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end do !j
            
      call MeshCommit(u%TowerMotion, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      if (errStat >= AbortErrLev) return

      
      u%TowerMotion%Orientation     = u%TowerMotion%RefOrientation
      u%TowerMotion%TranslationDisp = 0.0_R8Ki
      u%TowerMotion%TranslationVel  = 0.0_ReKi
      
   end if ! we compute tower loads
   
         !................
         ! hub
         !................
   
      call MeshCreate ( BlankMesh = u%HubMotion     &
                       ,IOS       = COMPONENT_INPUT &
                       ,Nnodes    = 1               &
                       ,ErrStat   = ErrStat2        &
                       ,ErrMess   = ErrMsg2         &
                       ,Orientation     = .true.    &
                       ,TranslationDisp = .true.    &
                       ,RotationVel     = .true.    &
                      )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      if (errStat >= AbortErrLev) return
                     
      call MeshPositionNode(u%HubMotion, 1, InitInp%HubPosition, errStat2, errMsg2, InitInp%HubOrientation)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
      call MeshConstructElement( u%HubMotion, ELEMENT_POINT, errStat2, errMsg2, p1=1 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      call MeshCommit(u%HubMotion, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      if (errStat >= AbortErrLev) return

         
      u%HubMotion%Orientation     = u%HubMotion%RefOrientation
      u%HubMotion%TranslationDisp = 0.0_R8Ki
      u%HubMotion%RotationVel     = 0.0_ReKi   
      
   
         !................
         ! blade roots
         !................
         
      allocate( u%BladeRootMotion(p%NumBlades), STAT = ErrStat2 )
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Fatal, 'Error allocating u%BladeRootMotion array.', ErrStat, ErrMsg, RoutineName )
         return
      end if      
      
      do k=1,p%NumBlades
         call MeshCreate ( BlankMesh = u%BladeRootMotion(k)                  &
                          ,IOS       = COMPONENT_INPUT                       &
                          ,Nnodes    = 1                                     &
                          ,ErrStat   = ErrStat2                              &
                          ,ErrMess   = ErrMsg2                               &
                          ,Orientation     = .true.                          &
                         )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

         if (errStat >= AbortErrLev) return
            
         call MeshPositionNode(u%BladeRootMotion(k), 1, InitInp%BladeRootPosition(:,k), errStat2, errMsg2, InitInp%BladeRootOrientation(:,:,k))
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
                     
         call MeshConstructElement( u%BladeRootMotion(k), ELEMENT_POINT, errStat2, errMsg2, p1=1 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
         call MeshCommit(u%BladeRootMotion(k), errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
         if (errStat >= AbortErrLev) return

      
         u%BladeRootMotion(k)%Orientation     = u%BladeRootMotion(k)%RefOrientation
   
   end do !k=numBlades      
      
      
         !................
         ! blades
         !................
   
      allocate( u%BladeMotion(p%NumBlades), STAT = ErrStat2 )
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Fatal, 'Error allocating u%BladeMotion array.', ErrStat, ErrMsg, RoutineName )
         return
      end if
      
      do k=1,p%NumBlades
         call MeshCreate ( BlankMesh = u%BladeMotion(k)                     &
                          ,IOS       = COMPONENT_INPUT                      &
                          ,Nnodes    = InputFileData%BladeProps(k)%NumBlNds &
                          ,ErrStat   = ErrStat2                             &
                          ,ErrMess   = ErrMsg2                              &
                          ,Orientation     = .true.                         &
                          ,TranslationDisp = .true.                         &
                          ,TranslationVel  = .true.                         &
                         )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

         if (errStat >= AbortErrLev) return
            
                        
         do j=1,InputFileData%BladeProps(k)%NumBlNds

               ! reference position of the jth node in the kth blade, relative to the root in the local blade coordinate system:
            positionL(1) = InputFileData%BladeProps(k)%BlCrvAC(j)
            positionL(2) = InputFileData%BladeProps(k)%BlSwpAC(j)
            positionL(3) = InputFileData%BladeProps(k)%BlSpn(  j)
            
               ! reference position of the jth node in the kth blade:
            position = u%BladeRootMotion(k)%Position(:,1) + matmul(positionL,u%BladeRootMotion(k)%RefOrientation(:,:,1))  ! note that because positionL is a 1-D array, we're doing the transpose of matmul(transpose(u%BladeRootMotion(k)%RefOrientation),positionL)

            
               ! reference orientation of the jth node in the kth blade, relative to the root in the local blade coordinate system:
            theta(1)     =  0.0_R8Ki
            theta(2)     =  InputFileData%BladeProps(k)%BlCrvAng(j)
            theta(3)     = -InputFileData%BladeProps(k)%BlTwist( j)            
            orientationL = EulerConstruct( theta )
                                 
               ! reference orientation of the jth node in the kth blade
            orientation = matmul( orientationL, u%BladeRootMotion(k)%RefOrientation(:,:,1) )

            
            call MeshPositionNode(u%BladeMotion(k), j, position, errStat2, errMsg2, orientation)
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
               
         end do ! j=blade nodes
         
            ! create line2 elements
         do j=1,InputFileData%BladeProps(k)%NumBlNds-1
            call MeshConstructElement( u%BladeMotion(k), ELEMENT_LINE2, errStat2, errMsg2, p1=j, p2=j+1 )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         end do !j
            
         call MeshCommit(u%BladeMotion(k), errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
         if (errStat >= AbortErrLev) return

      
         u%BladeMotion(k)%Orientation     = u%BladeMotion(k)%RefOrientation
         u%BladeMotion(k)%TranslationDisp = 0.0_R8Ki
         u%BladeMotion(k)%TranslationVel  = 0.0_ReKi
   
   end do !k=numBlades
   
   
end subroutine Init_u
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets AeroDyn parameters for use during the simulation; these variables are not changed after AD_Init.
subroutine SetParameters( InitInp, InputFileData, p, ErrStat, ErrMsg )
   TYPE(AD_InitInputType),       intent(in   )  :: InitInp          !< Input data for initialization routine, out is needed because of copy below
   TYPE(AD_InputFile),           INTENT(INout)  :: InputFileData    !< Data stored in the module's input file -- intent(out) only for move_alloc statements
   TYPE(AD_ParameterType),       INTENT(INOUT)  :: p                !< Parameters
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat          !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg           !< Error message if ErrStat /= ErrID_None


      ! Local variables
   CHARACTER(ErrMsgLen)                          :: ErrMsg2         ! temporary Error message if ErrStat /= ErrID_None
   INTEGER(IntKi)                                :: ErrStat2        ! temporary Error status of the operation
   !INTEGER(IntKi)                                :: i, j
   character(*), parameter                       :: RoutineName = 'SetParameters'
   
      ! Initialize variables for this routine

   ErrStat  = ErrID_None
   ErrMsg   = ""

   p%DT               = InputFileData%DTAero      
   p%WakeMod          = InputFileData%WakeMod
   p%TwrPotent        = InputFileData%TwrPotent
   p%TwrShadow        = InputFileData%TwrShadow
   p%TwrAero          = InputFileData%TwrAero
   p%CavitCheck       = InputFileData%CavitCheck
   p%Gravity          = InitInp%Gravity
  

   
   if (InitInp%Linearize) then
      p%FrozenWake = InputFileData%FrozenWake 
   else
      p%FrozenWake = .FALSE.
   end if
   
   
 ! p%numBlades        = InitInp%numBlades    ! this was set earlier because it was necessary
   p%NumBlNds         = InputFileData%BladeProps(1)%NumBlNds
   if (p%TwrPotent == TwrPotent_none .and. .not. p%TwrShadow .and. .not. p%TwrAero) then
      p%NumTwrNds     = 0
   else
      p%NumTwrNds     = InputFileData%NumTwrNds
      
      call move_alloc( InputFileData%TwrDiam, p%TwrDiam )
      call move_alloc( InputFileData%TwrCd,   p%TwrCd )      
   end if
   
   p%AirDens          = InputFileData%AirDens          
   p%KinVisc          = InputFileData%KinVisc
   p%Patm             = InputFileData%Patm
   p%Pvap             = InputFileData%Pvap
   p%FluidDepth       = InputFileData%FluidDepth
   p%SpdSound         = InputFileData%SpdSound
   
  !p%AFI     ! set in call to AFI_Init() [called early because it wants to use the same echo file as AD]
  !p%BEMT    ! set in call to BEMT_Init()
      
  !p%RootName       = TRIM(InitInp%RootName)//'.AD'   ! set earlier to it could be used   
   
#ifdef DBG_OUTS
   p%NBlOuts          = 23  
   p%numOuts          = p%NumBlNds*p%NumBlades*p%NBlOuts
   p%NTwOuts          = 0
      
#else
   p%numOuts          = InputFileData%NumOuts  
   p%NBlOuts          = InputFileData%NBlOuts      
   p%BlOutNd          = InputFileData%BlOutNd
   
   if (p%NumTwrNds > 0) then
      p%NTwOuts = InputFileData%NTwOuts
      p%TwOutNd = InputFileData%TwOutNd
   else
      p%NTwOuts = 0
   end if
   
   call SetOutParam(InputFileData%OutList, p, ErrStat2, ErrMsg2 ) ! requires: p%NumOuts, p%numBlades, p%NumBlNds, p%NumTwrNds; sets: p%OutParam.
      call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return  
   
#endif  
   
end subroutine SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine AD_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(AD_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(AD_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
      TYPE(AD_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(AD_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(AD_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(AD_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(AD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:



         ! Destroy the input data:

      CALL AD_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      CALL AD_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      CALL AD_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL AD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL AD_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL AD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
      CALL AD_DestroyMisc(        m,           ErrStat, ErrMsg ) 

         ! Destroy the output data:

      CALL AD_DestroyOutput( y, ErrStat, ErrMsg )




END SUBROUTINE AD_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine AD_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, m, errStat, errMsg )
!..................................................................................................................................

   real(DbKi),                     intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),                 intent(in   ) :: n          !< Current simulation time step n = 0,1,...
   type(AD_InputType),             intent(inout) :: u(:)       !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                     intent(in   ) :: utimes(:)  !< Times associated with u(:), in seconds
   type(AD_ParameterType),         intent(in   ) :: p          !< Parameters
   type(AD_ContinuousStateType),   intent(inout) :: x          !< Input: Continuous states at t;
                                                               !!   Output: Continuous states at t + Interval
   type(AD_DiscreteStateType),     intent(inout) :: xd         !< Input: Discrete states at t;
                                                               !!   Output: Discrete states at t  + Interval
   type(AD_ConstraintStateType),   intent(inout) :: z          !< Input: Constraint states at t;
                                                               !!   Output: Constraint states at t+dt
   type(AD_OtherStateType),        intent(inout) :: OtherState !< Input: Other states at t;
                                                               !!   Output: Other states at t+dt
   type(AD_MiscVarType),           intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),                 intent(  out) :: errStat    !< Error status of the operation
   character(*),                   intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(AD_InputType)                           :: uInterp     ! Interpolated/Extrapolated input
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'AD_UpdateStates'
      
   ErrStat = ErrID_None
   ErrMsg  = ""
     

   call AD_CopyInput( u(1), uInterp, MESH_NEWCOPY, errStat2, errMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

      ! set values of m%BEMT_u(2) from inputs interpolated at t+dt:
   call AD_Input_ExtrapInterp(u,utimes,uInterp,t+p%DT, errStat2, errMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call SetInputs(p, uInterp, m, 2, errStat2, errMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      ! set values of m%BEMT_u(1) from inputs (uInterp) interpolated at t:
      ! I'm doing this second in case we want the other misc vars at t as before, but I don't think it matters      
   call AD_Input_ExtrapInterp(u,utimes,uInterp, t, errStat2, errMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   call SetInputs(p, uInterp, m, 1, errStat2, errMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
                        
      ! Call into the BEMT update states    NOTE:  This is a non-standard framework interface!!!!!  GJH
   call BEMT_UpdateStates(t, n, m%BEMT_u(1), m%BEMT_u(2),  p%BEMT, x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, p%AFI, m%BEMT, errStat2, errMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
           
   call Cleanup()
   
contains
   subroutine Cleanup()
      call AD_DestroyInput( uInterp, errStat2, errMsg2)
   end subroutine Cleanup
end subroutine AD_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine AD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
! NOTE: no matter how many channels are selected for output, all of the outputs are calculated
! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(AD_InputType),           INTENT(IN   )  :: u           !< Inputs at Time t
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(AD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(AD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(AD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(AD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   type(AD_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   integer, parameter                           :: indx = 1  ! m%BEMT_u(1) is at t; m%BEMT_u(2) is t+dt
   integer(intKi)                               :: i
   integer(intKi)                               :: j

   integer(intKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'AD_CalcOutput'
   real(ReKi)                                   :: SigmaCavitCrit, SigmaCavit

   ErrStat = ErrID_None
   ErrMsg  = ""

   
   call SetInputs(p, u, m, indx, errStat2, errMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
   ! Call the BEMT module CalcOutput.  Notice that the BEMT outputs are purposely attached to AeroDyn's MiscVar structure to
   ! avoid issues with the coupling code
   
   call BEMT_CalcOutput(t, m%BEMT_u(indx), p%BEMT, x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, p%AFI, m%BEMT_y, m%BEMT, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  
   call SetOutputsFromBEMT(p, m, y )
                          
   if ( p%TwrAero ) then
      call ADTwr_CalcOutput(p, u, m, y, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
   end if
   
   if ( p%CavitCheck ) then      ! Calculate the cavitation number for the airfoil at the node in quesiton, and compare to the critical cavitation number based on the vapour pressure and submerged depth       
      do j = 1,p%numBlades ! Loop through all blades
         do i = 1,p%NumBlNds  ! Loop through all nodes
                     
      if ( EqualRealNos( m%BEMT_y%Vrel(i,j), 0.0_ReKi ) ) call SetErrStat( ErrID_Fatal, 'Vrel cannot be zero to do a cavitation check', ErrStat, ErrMsg, RoutineName) 
         if (ErrStat >= AbortErrLev) return
      
      SigmaCavit= -1* m%BEMT_y%Cpmin(i,j) ! Local cavitation number on node j                                               
      SigmaCavitCrit= ( ( p%Patm + ( p%Gravity * (p%FluidDepth - (  u%BladeMotion(j)%Position(3,i) + u%BladeMotion(j)%TranslationDisp(3,i) - u%HubMotion%Position(3,1))) * p%airDens)  - p%Pvap ) / ( 0.5_ReKi * p%airDens * m%BEMT_y%Vrel(i,j)**2)) ! Critical value of Sigma, cavitation occurs if local cavitation number is greater than this
                                                                  
         if ( (SigmaCavitCrit < SigmaCavit) .and. (.not. (m%CavitWarnSet(i,j)) ) ) then     
              call WrScr( NewLine//'Cavitation occurred at blade '//trim(num2lstr(j))//' and node '//trim(num2lstr(i))//'.' )
              m%CavitWarnSet(i,j) = .true.
         end if 
                     
      m%SigmaCavit(i,j)= SigmaCavit                 
      m%SigmaCavitCrit(i,j)=SigmaCavitCrit  
                           
         end do   ! p%NumBlNds
      end do  ! p%numBlades
   end if   ! Cavitation check
      

   !-------------------------------------------------------   
   !     get values to output to file:  
   !-------------------------------------------------------   
   if (p%NumOuts > 0) then
#ifdef DBG_OUTS
      call Calc_WriteDbgOutput( p, u, m, y, ErrStat2, ErrMsg2 ) 
#else
      call Calc_WriteOutput( p, u, m, y, OtherState, indx, ErrStat2, ErrMsg2 )   
#endif   
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
   
      !...............................................................................................................................   
      ! Place the selected output channels into the WriteOutput(:) array with the proper sign:
      !...............................................................................................................................   

      do i = 1,p%NumOuts  ! Loop through all selected output channels
#ifdef DBG_OUTS
         y%WriteOutput(i) = m%AllOuts( i )
#else
         y%WriteOutput(i) = p%OutParam(i)%SignM * m%AllOuts( p%OutParam(i)%Indx )
#endif

      end do             ! i - All selected output channels
      
   end if
   
   
   
end subroutine AD_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
subroutine AD_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )   :: Time        !< Current simulation time in seconds
   TYPE(AD_InputType),           INTENT(IN   )   :: u           !< Inputs at Time
   TYPE(AD_ParameterType),       INTENT(IN   )   :: p           !< Parameters
   TYPE(AD_ContinuousStateType), INTENT(IN   )   :: x           !< Continuous states at Time
   TYPE(AD_DiscreteStateType),   INTENT(IN   )   :: xd          !< Discrete states at Time
   TYPE(AD_ConstraintStateType), INTENT(IN   )   :: z           !< Constraint states at Time (possibly a guess)
   TYPE(AD_OtherStateType),      INTENT(IN   )   :: OtherState  !< Other states at Time
   TYPE(AD_MiscVarType),         INTENT(INOUT)   :: m           !< Misc/optimization variables
   TYPE(AD_ConstraintStateType), INTENT(INOUT)   :: Z_residual  !< Residual of the constraint state equations using
                                                                !!     the input values described above
   INTEGER(IntKi),               INTENT(  OUT)   :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)   :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   
      ! Local variables   
   integer, parameter                            :: indx = 1  ! m%BEMT_u(1) is at t; m%BEMT_u(2) is t+dt
   integer(intKi)                                :: ErrStat2
   character(ErrMsgLen)                          :: ErrMsg2
   character(*), parameter                       :: RoutineName = 'AD_CalcConstrStateResidual'
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   if (.not. allocated(Z_residual%BEMT%phi)) then ! BEMT_CalcConstrStateResidual expects memory to be allocated, so let's make sure it is
      call AD_CopyConstrState( z, Z_residual, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   end if
   
   
   call SetInputs(p, u, m, indx, errStat2, errMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                                
      
   call BEMT_CalcConstrStateResidual( Time, m%BEMT_u(indx), p%BEMT, x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, m%BEMT, &
                                       Z_residual%BEMT, p%AFI, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
   
   
end subroutine AD_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine converts the AeroDyn inputs into values that can be used for its submodules. It calculates the disturbed inflow
!! on the blade if tower shadow or tower influence are enabled, then uses these values to set m%BEMT_u(indx).
subroutine SetInputs(p, u, m, indx, errStat, errMsg)
   type(AD_ParameterType),       intent(in   )  :: p                      !< AD parameters
   type(AD_InputType),           intent(in   )  :: u                      !< AD Inputs at Time
   type(AD_MiscVarType),         intent(inout)  :: m                      !< Misc/optimization variables
   integer,                      intent(in   )  :: indx                   !< index into m%BEMT_u(indx) array; 1=t and 2=t+dt (but not checked here)
   integer(IntKi),               intent(  out)  :: ErrStat                !< Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg                 !< Error message if ErrStat /= ErrID_None
                                 
   ! local variables             
   integer(intKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'SetInputs'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   if (p%TwrPotent /= TwrPotent_none .or. p%TwrShadow) then
      call TwrInfl( p, u, m, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   else
      m%DisturbedInflow = u%InflowOnBlade
   end if
               
      ! This needs to extract the inputs from the AD data types (mesh) and massage them for the BEMT module
   call SetInputsForBEMT(p, u, m, indx, errStat2, errMsg2)  
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   
end subroutine SetInputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets m%BEMT_u(indx).
subroutine SetInputsForBEMT(p, u, m, indx, errStat, errMsg)

   type(AD_ParameterType),  intent(in   )  :: p                               !< AD parameters
   type(AD_InputType),      intent(in   )  :: u                               !< AD Inputs at Time
   type(AD_MiscVarType),    intent(inout)  :: m                               !< Misc/optimization variables
   integer,                 intent(in   )  :: indx                            !< index into m%BEMT_u array; must be 1 or 2 (but not checked here)
   integer(IntKi),          intent(  out)  :: ErrStat                         !< Error status of the operation
   character(*),            intent(  out)  :: ErrMsg                          !< Error message if ErrStat /= ErrID_None
      
   ! local variables
   real(ReKi)                              :: x_hat(3)
   real(ReKi)                              :: y_hat(3)
   real(ReKi)                              :: z_hat(3)
   real(ReKi)                              :: x_hat_disk(3)
   real(ReKi)                              :: y_hat_disk(3)
   real(ReKi)                              :: z_hat_disk(3)
   real(ReKi)                              :: tmp(3)
   real(R8Ki)                              :: theta(3)
   real(R8Ki)                              :: orientation(3,3)
   real(R8Ki)                              :: orientation_nopitch(3,3)
   real(ReKi)                              :: tmp_sz, tmp_sz_y
   
   integer(intKi)                          :: j                      ! loop counter for nodes
   integer(intKi)                          :: k                      ! loop counter for blades
   integer(intKi)                          :: ErrStat2
   character(ErrMsgLen)                    :: ErrMsg2
   character(*), parameter                 :: RoutineName = 'SetInputsForBEMT'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! calculate disk-averaged relative wind speed, V_DiskAvg
   m%V_diskAvg = 0.0_ReKi
   do k=1,p%NumBlades
      do j=1,p%NumBlNds
         tmp = m%DisturbedInflow(:,j,k) - u%BladeMotion(k)%TranslationVel(:,j)
         m%V_diskAvg = m%V_diskAvg + tmp         
      end do
   end do
   m%V_diskAvg = m%V_diskAvg / real( p%NumBlades * p%NumBlNds, ReKi ) 
   
      ! orientation vectors:
   x_hat_disk = u%HubMotion%Orientation(1,:,1) !actually also x_hat_hub      
   
   m%V_dot_x  = dot_product( m%V_diskAvg, x_hat_disk )
   m%BEMT_u(indx)%Un_disk  = m%V_dot_x
   tmp    = m%V_dot_x * x_hat_disk - m%V_diskAvg
   tmp_sz = TwoNorm(tmp)
   if ( EqualRealNos( tmp_sz, 0.0_ReKi ) ) then
      y_hat_disk = u%HubMotion%Orientation(2,:,1)
      z_hat_disk = u%HubMotion%Orientation(3,:,1)
   else
     y_hat_disk = tmp / tmp_sz
     z_hat_disk = cross_product( m%V_diskAvg, x_hat_disk ) / tmp_sz
  end if
     
      ! "Angular velocity of rotor" rad/s
   m%BEMT_u(indx)%omega   = dot_product( u%HubMotion%RotationVel(:,1), x_hat_disk )    
   
      ! "Angle between the vector normal to the rotor plane and the wind vector (e.g., the yaw angle in the case of no tilt)" rad 
   tmp_sz = TwoNorm( m%V_diskAvg )
   if ( EqualRealNos( tmp_sz, 0.0_ReKi ) ) then
      m%BEMT_u(indx)%chi0 = 0.0_ReKi
   else
         ! make sure we don't have numerical issues that make the ratio outside +/-1
      tmp_sz_y = min(  1.0_ReKi, m%V_dot_x / tmp_sz )
      tmp_sz_y = max( -1.0_ReKi, tmp_sz_y )
      
      m%BEMT_u(indx)%chi0 = acos( tmp_sz_y )
      
   end if
   
      ! "Azimuth angle" rad
   do k=1,p%NumBlades
      z_hat = u%BladeRootMotion(k)%Orientation(3,:,1)      
      tmp_sz_y = -1.0*dot_product(z_hat,y_hat_disk)
      tmp_sz   =      dot_product(z_hat,z_hat_disk)
      if ( EqualRealNos(tmp_sz_y,0.0_ReKi) .and. EqualRealNos(tmp_sz,0.0_ReKi) ) then
         m%BEMT_u(indx)%psi(k) = 0.0_ReKi
      else
         m%BEMT_u(indx)%psi(k) = atan2( tmp_sz_y, tmp_sz )
      end if      
   end do
   
      ! theta, "Twist angle (includes all sources of twist)" rad
      ! Vx, "Local axial velocity at node" m/s
      ! Vy, "Local tangential velocity at node" m/s
   do k=1,p%NumBlades
      
         ! construct system equivalent to u%BladeRootMotion(k)%Orientation, but without the blade-pitch angle:
      
      !orientation = matmul( u%BladeRootMotion(k)%Orientation(:,:,1), transpose(u%HubMotion%Orientation(:,:,1)) )
      call LAPACK_gemm( 'n', 't', 1.0_R8Ki, u%BladeRootMotion(k)%Orientation(:,:,1), u%HubMotion%Orientation(:,:,1), 0.0_R8Ki, orientation, errStat2, errMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      theta = EulerExtract( orientation ) !hub_theta_root(k)
#ifndef DBG_OUTS
      m%AllOuts( BPitch(  k) ) = -theta(3)*R2D ! save this value of pitch for potential output
#endif
      theta(3) = 0.0_ReKi  
      m%hub_theta_x_root(k) = theta(1)   ! save this value for FAST.Farm
      
      orientation = EulerConstruct( theta )
      orientation_nopitch = matmul( orientation, u%HubMotion%Orientation(:,:,1) ) ! withoutPitch_theta_Root(k)
            
      do j=1,p%NumBlNds         
         
            ! form coordinate system equivalent to u%BladeMotion(k)%Orientation(:,:,j) but without live sweep (due to in-plane
            ! deflection), blade-pitch and twist (aerodynamic + elastic) angles:
         
         ! orientation = matmul( u%BladeMotion(k)%Orientation(:,:,j), transpose(orientation_nopitch) )
         call LAPACK_gemm( 'n', 't', 1.0_R8Ki, u%BladeMotion(k)%Orientation(:,:,j), orientation_nopitch, 0.0_R8Ki, orientation, errStat2, errMsg2)
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         theta = EulerExtract( orientation ) !root(k)WithoutPitch_theta(j)_blade(k)
         
         m%BEMT_u(indx)%theta(j,k) = -theta(3) ! local pitch + twist (aerodyanmic + elastic) angle of the jth node in the kth blade
         
         
         theta(1) = 0.0_ReKi
         theta(3) = 0.0_ReKi
         m%Curve(j,k) = theta(2)  ! save value for possible output later
         m%WithoutSweepPitchTwist(:,:,j,k) = matmul( EulerConstruct( theta ), orientation_nopitch ) ! WithoutSweepPitch+Twist_theta(j)_Blade(k)
                           
         x_hat = m%WithoutSweepPitchTwist(1,:,j,k)
         y_hat = m%WithoutSweepPitchTwist(2,:,j,k)
         tmp   = m%DisturbedInflow(:,j,k) - u%BladeMotion(k)%TranslationVel(:,j) ! rel_V(j)_Blade(k)
         
         m%BEMT_u(indx)%Vx(j,k) = dot_product( tmp, x_hat ) ! normal component (normal to the plane, not chord) of the inflow velocity of the jth node in the kth blade
         m%BEMT_u(indx)%Vy(j,k) = dot_product( tmp, y_hat ) ! tangential component (tangential to the plane, not chord) of the inflow velocity of the jth node in the kth blade
         
      end do !j=nodes
   end do !k=blades
   
   
      ! "Radial distance from center-of-rotation to node" m
   
   do k=1,p%NumBlades
      do j=1,p%NumBlNds
         
            ! displaced position of the jth node in the kth blade relative to the hub:
         tmp =  u%BladeMotion(k)%Position(:,j) + u%BladeMotion(k)%TranslationDisp(:,j) &
              - u%HubMotion%Position(:,1)      - u%HubMotion%TranslationDisp(:,1)
         
            ! local radius (normalized distance from rotor centerline)
         tmp_sz_y = dot_product( tmp, y_hat_disk )**2
         tmp_sz   = dot_product( tmp, z_hat_disk )**2
         m%BEMT_u(indx)%rLocal(j,k) = sqrt( tmp_sz + tmp_sz_y )
         
      end do !j=nodes      
   end do !k=blades  
  
   m%BEMT_u(indx)%UserProp = u%UserProp
   
end subroutine SetInputsForBEMT
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine converts outputs from BEMT (stored in m%BEMT_y) into values on the AeroDyn BladeLoad output mesh.
subroutine SetOutputsFromBEMT(p, m, y )

   type(AD_ParameterType),  intent(in   )  :: p                               !< AD parameters
   type(AD_OutputType),     intent(inout)  :: y                               !< AD outputs 
   type(AD_MiscVarType),    intent(inout)  :: m                               !< Misc/optimization variables
   !type(BEMT_OutputType),   intent(in   )  :: BEMT_y                          ! BEMT outputs
   !real(ReKi),              intent(in   )  :: WithoutSweepPitchTwist(:,:,:,:) ! modified orientation matrix

   integer(intKi)                          :: j                      ! loop counter for nodes
   integer(intKi)                          :: k                      ! loop counter for blades
   real(reki)                              :: force(3)
   real(reki)                              :: moment(3)
   real(reki)                              :: q
   
  
   
   force(3)    =  0.0_ReKi          
   moment(1:2) =  0.0_ReKi          
   do k=1,p%NumBlades
      do j=1,p%NumBlNds
                      
         q = 0.5 * p%airDens * m%BEMT_y%Vrel(j,k)**2              ! dynamic pressure of the jth node in the kth blade
         force(1) =  m%BEMT_y%cx(j,k) * q * p%BEMT%chord(j,k)     ! X = normal force per unit length (normal to the plane, not chord) of the jth node in the kth blade
         force(2) = -m%BEMT_y%cy(j,k) * q * p%BEMT%chord(j,k)     ! Y = tangential force per unit length (tangential to the plane, not chord) of the jth node in the kth blade
         moment(3)=  m%BEMT_y%cm(j,k) * q * p%BEMT%chord(j,k)**2  ! M = pitching moment per unit length of the jth node in the kth blade
         
            ! save these values for possible output later:
         m%X(j,k) = force(1)
         m%Y(j,k) = force(2)
         m%M(j,k) = moment(3)
         
            ! note: because force and moment are 1-d arrays, I'm calculating the transpose of the force and moment outputs
            !       so that I don't have to take the transpose of WithoutSweepPitchTwist(:,:,j,k)
         y%BladeLoad(k)%Force(:,j)  = matmul( force,  m%WithoutSweepPitchTwist(:,:,j,k) )  ! force per unit length of the jth node in the kth blade
         y%BladeLoad(k)%Moment(:,j) = matmul( moment, m%WithoutSweepPitchTwist(:,:,j,k) )  ! moment per unit length of the jth node in the kth blade
         
      end do !j=nodes
   end do !k=blades
   
   
end subroutine SetOutputsFromBEMT
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine validates the inputs from the AeroDyn input files.
SUBROUTINE ValidateInputData( InitInp, InputFileData, NumBl, ErrStat, ErrMsg )
!..................................................................................................................................
      
      ! Passed variables:

   type(AD_InitInputType),   intent(in   )  :: InitInp                           !< Input data for initialization routine
   type(AD_InputFile),       intent(in)     :: InputFileData                     !< All the data in the AeroDyn input file
   integer(IntKi),           intent(in)     :: NumBl                             !< Number of blades
   integer(IntKi),           intent(out)    :: ErrStat                           !< Error status
   character(*),             intent(out)    :: ErrMsg                            !< Error message

   
      ! local variables
   integer(IntKi)                           :: k                                 ! Blade number
   integer(IntKi)                           :: j                                 ! node number
   character(*), parameter                  :: RoutineName = 'ValidateInputData'
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   if (NumBl > MaxBl .or. NumBl < 1) call SetErrStat( ErrID_Fatal, 'Number of blades must be between 1 and '//trim(num2lstr(MaxBl))//'.', ErrSTat, ErrMsg, RoutineName )
   if (InputFileData%DTAero <= 0.0)  call SetErrStat ( ErrID_Fatal, 'DTAero must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%WakeMod /= WakeMod_None .and. InputFileData%WakeMod /= WakeMod_BEMT .and. InputFileData%WakeMod /= WakeMod_DBEMT) then
      call SetErrStat ( ErrID_Fatal, 'WakeMod must '//trim(num2lstr(WakeMod_None))//' (none), '//trim(num2lstr(WakeMod_BEMT))//' (BEMT),'// &
         'or '//trim(num2lstr(WakeMod_DBEMT))//' (DBEMT).', ErrStat, ErrMsg, RoutineName ) 
   end if
   
   if (InputFileData%AFAeroMod /= AFAeroMod_Steady .and. InputFileData%AFAeroMod /= AFAeroMod_BL_unsteady) then
      call SetErrStat ( ErrID_Fatal, 'AFAeroMod must be '//trim(num2lstr(AFAeroMod_Steady))//' (steady) or '//&
                        trim(num2lstr(AFAeroMod_BL_unsteady))//' (Beddoes-Leishman unsteady).', ErrStat, ErrMsg, RoutineName ) 
   end if
   if (InputFileData%TwrPotent /= TwrPotent_none .and. InputFileData%TwrPotent /= TwrPotent_baseline .and. InputFileData%TwrPotent /= TwrPotent_Bak) then
      call SetErrStat ( ErrID_Fatal, 'TwrPotent must be 0 (none), 1 (baseline potential flow), or 2 (potential flow with Bak correction).', ErrStat, ErrMsg, RoutineName ) 
   end if   
   
   if (InputFileData%AirDens <= 0.0) call SetErrStat ( ErrID_Fatal, 'The air density (AirDens) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%KinVisc <= 0.0) call SetErrStat ( ErrID_Fatal, 'The kinesmatic viscosity (KinVisc) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%SpdSound <= 0.0) call SetErrStat ( ErrID_Fatal, 'The speed of sound (SpdSound) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%Pvap <= 0.0) call SetErrStat ( ErrID_Fatal, 'The vapour pressure (Pvap) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%Patm <= 0.0) call SetErrStat ( ErrID_Fatal, 'The atmospheric pressure (Patm)  must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%FluidDepth <= 0.0) call SetErrStat ( ErrID_Fatal, 'Fluid depth (FluidDepth) must be greater than zero', ErrStat, ErrMsg, RoutineName )

                       
   
      ! BEMT/DBEMT inputs
      ! bjj: these checks should probably go into BEMT where they are used...
   if (InputFileData%WakeMod /= WakeMod_none) then
      if ( InputFileData%MaxIter < 1 ) call SetErrStat( ErrID_Fatal, 'MaxIter must be greater than 0.', ErrStat, ErrMsg, RoutineName )
      
      if ( InputFileData%IndToler < 0.0 .or. EqualRealNos(InputFileData%IndToler, 0.0_ReKi) ) &
         call SetErrStat( ErrID_Fatal, 'IndToler must be greater than 0.', ErrStat, ErrMsg, RoutineName )
   
      if ( InputFileData%SkewMod /= SkewMod_Uncoupled .and. InputFileData%SkewMod /= SkewMod_PittPeters) &  !  .and. InputFileData%SkewMod /= SkewMod_Coupled )
           call SetErrStat( ErrID_Fatal, 'SkewMod must be 1, or 2.  Option 3 will be implemented in a future version.', ErrStat, ErrMsg, RoutineName )      
      
   end if !BEMT/DBEMT checks
   
   
   if ( InputFileData%CavitCheck .and. InputFileData%AFAeroMod == AFAeroMod_BL_unsteady) then
      call SetErrStat( ErrID_Fatal, 'Cannot use unsteady aerodynamics module with a cavitation check', ErrStat, ErrMsg, RoutineName )
   end if
        
   if (InputFileData%InCol_Cpmin == 0 .and. InputFileData%CavitCheck) call SetErrStat( ErrID_Fatal, 'InCol_Cpmin must not be 0 to do a cavitation check.', ErrStat, ErrMsg, RoutineName )

         ! validate the number of airfoils
   if (InputFileData%NumAFfiles  < 1) call SetErrStat( ErrID_Fatal, 'The number of unique airfoil tables (NumAFfiles) must be greater than zero.', ErrStat, ErrMsg, RoutineName )   
   
      ! .............................
      ! check blade mesh data:
      ! .............................
   if ( InputFileData%BladeProps(1)%NumBlNds < 2 ) call SetErrStat( ErrID_Fatal, 'There must be at least two nodes per blade.',ErrStat, ErrMsg, RoutineName )
   do k=2,NumBl
      if ( InputFileData%BladeProps(k)%NumBlNds /= InputFileData%BladeProps(k-1)%NumBlNds ) then
         call SetErrStat( ErrID_Fatal, 'All blade property files must have the same number of blade nodes.', ErrStat, ErrMsg, RoutineName )
         exit  ! exit do loop
      end if
   end do
   
      ! Check the list of airfoil tables for blades to make sure they are all within limits.
   do k=1,NumBl
      do j=1,InputFileData%BladeProps(k)%NumBlNds
         if ( ( InputFileData%BladeProps(k)%BlAFID(j) < 1 ) .OR. ( InputFileData%BladeProps(k)%BlAFID(j) > InputFileData%NumAFfiles ) )  then
            call SetErrStat( ErrID_Fatal, 'Blade '//trim(Num2LStr(k))//' node '//trim(Num2LStr(j))//' must be a number between 1 and NumAFfiles (' &
               //TRIM(Num2LStr(InputFileData%NumAFfiles))//').', ErrStat, ErrMsg, RoutineName )
         end if
      end do ! j=nodes
   end do ! k=blades
            
      ! Check that the blade chord is > 0.
   do k=1,NumBl
      do j=1,InputFileData%BladeProps(k)%NumBlNds
         if ( InputFileData%BladeProps(k)%BlChord(j) <= 0.0_ReKi )  then
            call SetErrStat( ErrID_Fatal, 'The chord for blade '//trim(Num2LStr(k))//' node '//trim(Num2LStr(j)) &
                             //' must be greater than 0.', ErrStat, ErrMsg, RoutineName )
         endif
      end do ! j=nodes
   end do ! k=blades
   
   do k=1,NumBl
      if ( .not. EqualRealNos(InputFileData%BladeProps(k)%BlSpn(1), 0.0_ReKi) ) call SetErrStat( ErrID_Fatal, 'Blade '//trim(Num2LStr(k))//' span location must start at 0.0 m', ErrStat, ErrMsg, RoutineName)       
      do j=2,InputFileData%BladeProps(k)%NumBlNds
         if ( InputFileData%BladeProps(k)%BlSpn(j) <= InputFileData%BladeProps(k)%BlSpn(j-1) )  then
            call SetErrStat( ErrID_Fatal, 'Blade '//trim(Num2LStr(k))//' nodes must be entered in increasing elevation.', ErrStat, ErrMsg, RoutineName )
            exit
         end if
      end do ! j=nodes
   end do ! k=blades
   
      ! .............................
      ! check tower mesh data:
      ! .............................
   if (InputFileData%TwrPotent /= TwrPotent_none .or. InputFileData%TwrShadow .or. InputFileData%TwrAero ) then
      
      if (InputFileData%NumTwrNds < 2) call SetErrStat( ErrID_Fatal, 'There must be at least two nodes on the tower.',ErrStat, ErrMsg, RoutineName )
         
         ! Check that the tower diameter is > 0.
      do j=1,InputFileData%NumTwrNds
         if ( InputFileData%TwrDiam(j) <= 0.0_ReKi )  then
            call SetErrStat( ErrID_Fatal, 'The diameter for tower node '//trim(Num2LStr(j))//' must be greater than 0.' &
                            , ErrStat, ErrMsg, RoutineName )
         end if
      end do ! j=nodes
      
         ! check that the elevation is increasing:
      do j=2,InputFileData%NumTwrNds
         if ( InputFileData%TwrElev(j) <= InputFileData%TwrElev(j-1) )  then
            call SetErrStat( ErrID_Fatal, 'The tower nodes must be entered in increasing elevation.', ErrStat, ErrMsg, RoutineName )
            exit
         end if
      end do ! j=nodes
            
   end if
   
      ! .............................
      ! check outputs:
      ! .............................
   
   if ( ( InputFileData%NTwOuts < 0_IntKi ) .OR. ( InputFileData%NTwOuts > 9_IntKi ) )  then
      call SetErrStat( ErrID_Fatal, 'NTwOuts must be between 0 and 9 (inclusive).', ErrStat, ErrMsg, RoutineName )
   else
         ! Check to see if all TwOutNd(:) analysis points are existing analysis points:

      do j=1,InputFileData%NTwOuts
         if ( InputFileData%TwOutNd(j) < 1_IntKi .OR. InputFileData%TwOutNd(j) > InputFileData%NumTwrNds ) then
            call SetErrStat( ErrID_Fatal, ' All TwOutNd values must be between 1 and '//&
                           trim( Num2LStr( InputFileData%NumTwrNds ) )//' (inclusive).', ErrStat, ErrMsg, RoutineName )
            exit ! stop checking this loop
         end if
      end do         
   
   end if
         
         
   if ( ( InputFileData%NBlOuts < 0_IntKi ) .OR. ( InputFileData%NBlOuts > 9_IntKi ) )  then
      call SetErrStat( ErrID_Fatal, 'NBlOuts must be between 0 and 9 (inclusive).', ErrStat, ErrMsg, RoutineName )
   else 

   ! Check to see if all BlOutNd(:) analysis points are existing analysis points:

      do j=1,InputFileData%NBlOuts
         if ( InputFileData%BlOutNd(j) < 1_IntKi .OR. InputFileData%BlOutNd(j) > InputFileData%BladeProps(1)%NumBlNds ) then
            call SetErrStat( ErrID_Fatal, ' All BlOutNd values must be between 1 and '//&
                    trim( Num2LStr( InputFileData%BladeProps(1)%NumBlNds ) )//' (inclusive).', ErrStat, ErrMsg, RoutineName )
            exit ! stop checking this loop
         end if
      end do
      
   end if   
   
   !..................
   ! check for linearization
   !..................
   if (InitInp%Linearize) then
      if (InputFileData%AFAeroMod /= AFAeroMod_Steady) then
         call SetErrStat( ErrID_Fatal, 'Steady blade airfoil aerodynamics must be used for linearization. Set AFAeroMod=1.', ErrStat, ErrMsg, RoutineName )
      end if
      
      if (InputFileData%WakeMod == WakeMod_DBEMT) then
         call SetErrStat( ErrID_Fatal, 'DBEMT cannot currently be used for linearization. Set WakeMod=0 or WakeMod=1.', ErrStat, ErrMsg, RoutineName )
      end if
   end if
   
   
END SUBROUTINE ValidateInputData
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets up the data structures and initializes AirfoilInfo to get the necessary AFI parameters. It then verifies 
!! that the UA parameters are included in the AFI tables if UA is being used.
SUBROUTINE Init_AFIparams( InputFileData, p_AFI, UnEc, NumBl, ErrStat, ErrMsg )


      ! Passed variables
   type(AD_InputFile),                   intent(inout) :: InputFileData      !< All the data in the AeroDyn input file (intent(out) only because of the call to MOVE_ALLOC)
   type(AFI_ParameterType), allocatable, intent(  out) :: p_AFI(:)           !< parameters returned from the AFI (airfoil info) module
   integer(IntKi),                       intent(in   ) :: UnEc               !< I/O unit for echo file. If > 0, file is open for writing.
   integer(IntKi),                       intent(in   ) :: NumBl              !< number of blades (for performing check on valid airfoil data read in)
   integer(IntKi),                       intent(  out) :: ErrStat            !< Error status
   character(*),                         intent(  out) :: ErrMsg             !< Error message

      ! local variables
   type(AFI_InitInputType)                             :: AFI_InitInputs     ! initialization data for the AFI routines
   
   integer(IntKi)                                      :: j                  ! loop counter for nodes
   integer(IntKi)                                      :: k                  ! loop counter for blades
   integer(IntKi)                                      :: File               ! loop counter for airfoil files
   logical, allocatable                                :: fileUsed(:)
   
   integer(IntKi)                                      :: ErrStat2
   character(ErrMsgLen)                                :: ErrMsg2
   character(*), parameter                             :: RoutineName = 'Init_AFIparams'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   allocate(p_AFI( InputFileData%NumAFfiles), STAT = ErrStat2)
      if ( ErrStat2 /= 0 ) then
         call SetErrStat(ErrID_Fatal,'Error allocating p_AFI.',ErrStat,ErrMsg,RoutineName)
         return
      end if
   
   
      ! Setup Airfoil InitInput data structure:
   AFI_InitInputs%InCol_Alfa  = InputFileData%InCol_Alfa
   AFI_InitInputs%InCol_Cl    = InputFileData%InCol_Cl
   AFI_InitInputs%InCol_Cd    = InputFileData%InCol_Cd
   AFI_InitInputs%InCol_Cm    = InputFileData%InCol_Cm
   AFI_InitInputs%InCol_Cpmin = InputFileData%InCol_Cpmin
   AFI_InitInputs%AFTabMod    = InputFileData%AFTabMod !AFITable_1
   
      ! Call AFI_Init to read in and process the airfoil files.
      ! This includes creating the spline coefficients to be used for interpolation.
   
   do File = 1, InputFileData%NumAFfiles

      AFI_InitInputs%FileName = InputFileData%AFNames(File)

      call AFI_Init ( AFI_InitInputs, p_AFI(File), ErrStat2, ErrMsg2, UnEc )
         call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if (ErrStat >= AbortErrLev) exit
   end do
         
      
   call AFI_DestroyInitInput( AFI_InitInputs, ErrStat2, ErrMsg2 )
   if (ErrStat >= AbortErrLev) return
   
   
      ! check that we read the correct airfoil parameters for UA:      
   if ( InputFileData%AFAeroMod == AFAeroMod_BL_unsteady ) then
      
         ! determine which airfoil files will be used
      call AllocAry( fileUsed, InputFileData%NumAFfiles, 'fileUsed', errStat2, errMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)   
         if (errStat >= AbortErrLev) return
      fileUsed = .false.
            
      do k=1,NumBl
         do j=1,InputFileData%BladeProps(k)%NumBlNds
            fileUsed ( InputFileData%BladeProps(k)%BlAFID(j) ) = .true.
         end do ! j=nodes
      end do ! k=blades
      
         ! make sure all files in use have proper UA input parameters:
      do File = 1,InputFileData%NumAFfiles
         
         if (fileUsed(File)) then
            call UA_ValidateAFI(p_AFI(File), InputFileData%AFNames(File), ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
               if (errStat >= AbortErrLev) return
         end if
         
      end do
      
      if ( allocated(fileUsed) ) deallocate(fileUsed)
      
   end if
   
   
END SUBROUTINE Init_AFIparams
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the BEMT module from within AeroDyn.
SUBROUTINE Init_BEMTmodule( InputFileData, u_AD, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

   type(AD_InputFile),             intent(in   ) :: InputFileData  !< All the data in the AeroDyn input file
   type(AD_InputType),             intent(in   ) :: u_AD           !< AD inputs - used for input mesh node positions
   type(BEMT_InputType),           intent(  out) :: u              !< An initial guess for the input; input mesh must be defined
   type(AD_ParameterType),         intent(inout) :: p              !< Parameters ! intent out b/c we set the BEMT parameters here
   type(BEMT_ContinuousStateType), intent(  out) :: x              !< Initial continuous states
   type(BEMT_DiscreteStateType),   intent(  out) :: xd             !< Initial discrete states
   type(BEMT_ConstraintStateType), intent(  out) :: z              !< Initial guess of the constraint states
   type(BEMT_OtherStateType),      intent(  out) :: OtherState     !< Initial other states
   type(BEMT_OutputType),          intent(  out) :: y              !< Initial system outputs (outputs are not calculated;
                                                                   !!   only the output mesh is initialized)
   type(BEMT_MiscVarType),         intent(  out) :: m              !< Initial misc/optimization variables
   integer(IntKi),                 intent(  out) :: errStat        !< Error status of the operation
   character(*),                   intent(  out) :: errMsg         !< Error message if ErrStat /= ErrID_None


      ! Local variables
   real(DbKi)                                    :: Interval       ! Coupling interval in seconds: the rate that
                                                                   !   (1) BEMT_UpdateStates() is called in loose coupling &
                                                                   !   (2) BEMT_UpdateDiscState() is called in tight coupling.
                                                                   !   Input is the suggested time from the glue code;
                                                                   !   Output is the actual coupling interval that will be used
                                                                   !   by the glue code.
   type(BEMT_InitInputType)                      :: InitInp        ! Input data for initialization routine
   type(BEMT_InitOutputType)                     :: InitOut        ! Output for initialization routine
                                                 
   integer(intKi)                                :: j              ! node index
   integer(intKi)                                :: k              ! blade index
   real(ReKi)                                    :: tmp(3), tmp_sz_y, tmp_sz
   real(ReKi)                                    :: y_hat_disk(3)
   real(ReKi)                                    :: z_hat_disk(3)
   integer(IntKi)                                :: ErrStat2
   character(ErrMsgLen)                          :: ErrMsg2
   character(*), parameter                       :: RoutineName = 'Init_BEMTmodule'

   ! note here that each blade is required to have the same number of nodes
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! set initialization data here:   
   Interval                 = p%DT   
   InitInp%numBlades        = p%NumBlades
   
   InitInp%airDens          = InputFileData%AirDens 
   InitInp%kinVisc          = InputFileData%KinVisc                  
   InitInp%skewWakeMod      = InputFileData%SkewMod
   InitInp%yawCorrFactor    = InputFileData%SkewModFactor
   InitInp%aTol             = InputFileData%IndToler
   InitInp%useTipLoss       = InputFileData%TipLoss
   InitInp%useHubLoss       = InputFileData%HubLoss
   InitInp%useInduction     = InputFileData%WakeMod /= WakeMod_none
   InitInp%useTanInd        = InputFileData%TanInd
   InitInp%useAIDrag        = InputFileData%AIDrag        
   InitInp%useTIDrag        = InputFileData%TIDrag  
   InitInp%numBladeNodes    = p%NumBlNds
   InitInp%numReIterations  = 1                              ! This is currently not available in the input file and is only for testing  
   InitInp%maxIndIterations = InputFileData%MaxIter 
   
   call AllocAry(InitInp%chord, InitInp%numBladeNodes,InitInp%numBlades,'chord', ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   call AllocAry(InitInp%AFindx,InitInp%numBladeNodes,InitInp%numBlades,'AFindx',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   call AllocAry(InitInp%zHub,                        InitInp%numBlades,'zHub',  ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(InitInp%zLocal,InitInp%numBladeNodes,InitInp%numBlades,'zLocal',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   call AllocAry(InitInp%rLocal,InitInp%numBladeNodes,InitInp%numBlades,'rLocal',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   call AllocAry(InitInp%zTip,                        InitInp%numBlades,'zTip',  ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
   
   if ( ErrStat >= AbortErrLev ) then
      call Cleanup()
      return
   end if  

   
   do k=1,p%numBlades
      
      InitInp%zHub(k) = TwoNorm( u_AD%BladeRootMotion(k)%Position(:,1) - u_AD%HubMotion%Position(:,1) )  
      if (EqualRealNos(InitInp%zHub(k),0.0_ReKi) ) &
         call SetErrStat( ErrID_Fatal, "zHub for blade "//trim(num2lstr(k))//" is zero.", ErrStat, ErrMsg, RoutineName)
      
      InitInp%zLocal(1,k) = InitInp%zHub(k) + TwoNorm( u_AD%BladeMotion(k)%Position(:,1) - u_AD%BladeRootMotion(k)%Position(:,1) )
      do j=2,p%NumBlNds
         InitInp%zLocal(j,k) = InitInp%zLocal(j-1,k) + TwoNorm( u_AD%BladeMotion(k)%Position(:,j) - u_AD%BladeMotion(k)%Position(:,j-1) ) 
      end do !j=nodes
      
      InitInp%zTip(k) = InitInp%zLocal(p%NumBlNds,k)
      
      y_hat_disk = u_AD%HubMotion%Orientation(2,:,1)
      z_hat_disk = u_AD%HubMotion%Orientation(3,:,1)
      
      do j=1,p%NumBlNds
               ! displaced position of the jth node in the kth blade relative to the hub:
         tmp =  u_AD%BladeMotion(k)%Position(:,j)  - u_AD%HubMotion%Position(:,1) 
            ! local radius (normalized distance from rotor centerline)
         tmp_sz_y = dot_product( tmp, y_hat_disk )**2
         tmp_sz   = dot_product( tmp, z_hat_disk )**2
         InitInp%rLocal(j,k) = sqrt( tmp_sz + tmp_sz_y )
      end do !j=nodes   
   end do !k=blades
   
               
  do k=1,p%numBlades
     do j=1,p%NumBlNds
        InitInp%chord (j,k)  = InputFileData%BladeProps(k)%BlChord(j)
        InitInp%AFindx(j,k)  = InputFileData%BladeProps(k)%BlAFID(j)
     end do
  end do
   
   InitInp%UA_Flag    = InputFileData%AFAeroMod == AFAeroMod_BL_unsteady
   InitInp%UAMod      = InputFileData%UAMod
   InitInp%Flookup    = InputFileData%Flookup
   InitInp%a_s        = InputFileData%SpdSound
   
   if (InputFileData%WakeMod == WakeMod_DBEMT) then
      InitInp%DBEMT_Mod  = InputFileData%DBEMT_Mod
   else
      InitInp%DBEMT_Mod  = DBEMT_none
   end if
   InitInp%tau1_const = InputFileData%tau1_const
   
   if (ErrStat >= AbortErrLev) then
      call cleanup()
      return
   end if
   
   
   call BEMT_Init(InitInp, u, p%BEMT,  x, xd, z, OtherState, p%AFI, y, m, Interval, InitOut, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)   
         
   if (.not. equalRealNos(Interval, p%DT) ) &
      call SetErrStat( ErrID_Fatal, "DTAero was changed in Init_BEMTmodule(); this is not allowed.", ErrStat2, ErrMsg2, RoutineName)
   
   !m%UseFrozenWake = .FALSE. !BJJ: set this in BEMT
   
   call Cleanup()
   return
      
contains   
   subroutine Cleanup()
      call BEMT_DestroyInitInput( InitInp, ErrStat2, ErrMsg2 )   
      call BEMT_DestroyInitOutput( InitOut, ErrStat2, ErrMsg2 )   
   end subroutine Cleanup
   
END SUBROUTINE Init_BEMTmodule
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the tower loads for the AeroDyn TowerLoad output mesh.
SUBROUTINE ADTwr_CalcOutput(p, u, m, y, ErrStat, ErrMsg )

   TYPE(AD_InputType),           INTENT(IN   )  :: u           !< Inputs at Time t
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(AD_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
   TYPE(AD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   INTEGER(IntKi)                               :: j
   real(ReKi)                                   :: q
   real(ReKi)                                   :: V_rel(3)    ! relative wind speed on a tower node
   real(ReKi)                                   :: VL(2)       ! relative local x- and y-components of the wind speed on a tower node
   real(ReKi)                                   :: tmp(3)
   
   !integer(intKi)                               :: ErrStat2
   !character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'ADTwr_CalcOutput'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""

   
   do j=1,p%NumTwrNds
      
      V_rel = u%InflowOnTower(:,j) - u%TowerMotion%TranslationVel(:,j) ! relative wind speed at tower node
   
      tmp   = u%TowerMotion%Orientation(1,:,j)
      VL(1) = dot_product( V_Rel, tmp )            ! relative local x-component of wind speed of the jth node in the tower
      tmp   = u%TowerMotion%Orientation(2,:,j)
      VL(2) = dot_product( V_Rel, tmp )            ! relative local y-component of wind speed of the jth node in the tower
      
      m%W_Twr(j)  =  TwoNorm( VL )            ! relative wind speed normal to the tower at node j      
      q     = 0.5 * p%TwrCd(j) * p%AirDens * p%TwrDiam(j) * m%W_Twr(j)
      
         ! force per unit length of the jth node in the tower
      tmp(1) = q * VL(1)
      tmp(2) = q * VL(2)
      tmp(3) = 0.0_ReKi
      
      y%TowerLoad%force(:,j) = matmul( tmp, u%TowerMotion%Orientation(:,:,j) ) ! note that I'm calculating the transpose here, which is okay because we have 1-d arrays
      m%X_Twr(j) = tmp(1)
      m%Y_Twr(j) = tmp(2)
      
      
         ! moment per unit length of the jth node in the tower
      y%TowerLoad%moment(:,j) = 0.0_ReKi
      
   end do
   

END SUBROUTINE ADTwr_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine checks for invalid inputs to the tower influence models.
SUBROUTINE CheckTwrInfl(u, ErrStat, ErrMsg )

   TYPE(AD_InputType),           INTENT(IN   )  :: u           !< Inputs at Time t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
   
   ! local variables
   real(reKi)                                   :: ElemSize
   real(reKi)                                   :: tmp(3)
   integer(intKi)                               :: j
   character(*), parameter                      :: RoutineName = 'CheckTwrInfl'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !! the tower-influence models (tower potential flow and tower shadow) are valid only for small tower deflections;
   !! so, first throw an error to avoid a division-by-zero error if any line2 elements on the tower mesh are colocated.
   
   do j = 2,u%TowerMotion%Nnodes
      tmp =   u%TowerMotion%Position(:,j  ) + u%TowerMotion%TranslationDisp(:,j  ) &
            - u%TowerMotion%Position(:,j-1) - u%TowerMotion%TranslationDisp(:,j-1)
   
      ElemSize = TwoNorm(tmp)
      if ( EqualRealNos(ElemSize,0.0_ReKi) ) then
         call SetErrStat(ErrID_Fatal, "Division by zero:Elements "//trim(num2lstr(j))//' and '//trim(num2lstr(j-1))//' are colocated.', ErrStat, ErrMsg, RoutineName )
         exit
      end if
   end do
      
   
END SUBROUTINE CheckTwrInfl
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates m%DisturbedInflow, the influence of tower shadow and/or potential flow on the inflow velocities
SUBROUTINE TwrInfl( p, u, m, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(AD_InputType),           INTENT(IN   )  :: u                       !< Inputs at Time t
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p                       !< Parameters
   type(AD_MiscVarType),         intent(inout)  :: m                       !< Misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat                 !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg                  !< Error message if ErrStat /= ErrID_None

   ! local variables
   real(ReKi)                                   :: xbar                    ! local x^ component of r_TowerBlade (distance from tower to blade) normalized by tower radius
   real(ReKi)                                   :: ybar                    ! local y^ component of r_TowerBlade (distance from tower to blade) normalized by tower radius
   real(ReKi)                                   :: zbar                    ! local z^ component of r_TowerBlade (distance from tower to blade) normalized by tower radius
   real(ReKi)                                   :: theta_tower_trans(3,3)  ! transpose of local tower orientation expressed as a DCM
   real(ReKi)                                   :: TwrCd                   ! local tower drag coefficient
   real(ReKi)                                   :: W_tower                 ! local relative wind speed normal to the tower

   real(ReKi)                                   :: BladeNodePosition(3)    ! local blade node position
   
   
   real(ReKi)                                   :: u_TwrShadow             ! axial velocity deficit fraction from tower shadow
   real(ReKi)                                   :: u_TwrPotent             ! axial velocity deficit fraction from tower potential flow
   real(ReKi)                                   :: v_TwrPotent             ! transverse velocity deficit fraction from tower potential flow
   
   real(ReKi)                                   :: denom                   ! denominator
   real(ReKi)                                   :: v(3)                    ! temp vector
   
   integer(IntKi)                               :: j, k                    ! loop counters for elements, blades
   integer(intKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'TwrInfl'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""   
   
   
      ! these models are valid for only small tower deflections; check for potential division-by-zero errors:   
   call CheckTwrInfl( u, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return
      
   do k = 1, p%NumBlades
      do j = 1, u%BladeMotion(k)%NNodes
         
         ! for each line2-element node of the blade mesh, a nearest-neighbor line2 element or node of the tower 
         ! mesh is found in the deflected configuration, returning theta_tower, W_tower, xbar, ybar, zbar, and TowerCd:
         
         BladeNodePosition = u%BladeMotion(k)%Position(:,j) + u%BladeMotion(k)%TranslationDisp(:,j)
         
         call getLocalTowerProps(p, u, BladeNodePosition, theta_tower_trans, W_tower, xbar, ybar, zbar, TwrCd, m%TwrClrnc(j,k), ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            if (ErrStat >= AbortErrLev) return
         
      
         ! calculate tower influence:
         if ( abs(zbar) < 1.0_ReKi .and. p%TwrPotent /= TwrPotent_none ) then
            if ( p%TwrPotent == TwrPotent_baseline ) then
               
               denom = (xbar**2 + ybar**2)**2
               
               u_TwrPotent = ( -1.0*xbar**2 + ybar**2 ) / denom
               v_TwrPotent = ( -2.0*xbar    * ybar    ) / denom      
               
            elseif (p%TwrPotent == TwrPotent_Bak) then
               
               xbar = xbar + 0.1
               
               denom = (xbar**2 + ybar**2)**2               
               u_TwrPotent = ( -1.0*xbar**2 + ybar**2 ) / denom
               v_TwrPotent = ( -2.0*xbar    * ybar    ) / denom        
               
               denom = TwoPi*(xbar**2 + ybar**2)               
               u_TwrPotent = u_TwrPotent + TwrCd*xbar / denom
               v_TwrPotent = v_TwrPotent + TwrCd*ybar / denom                       
               
            end if
         else
            u_TwrPotent = 0.0_ReKi
            v_TwrPotent = 0.0_ReKi
         end if
         
         if ( p%TwrShadow .and. xbar > 0.0_ReKi .and. abs(zbar) < 1.0_ReKi) then
            denom = sqrt( sqrt( xbar**2 + ybar**2 ) )
            if ( abs(ybar) < denom ) then
               u_TwrShadow = -TwrCd / denom * cos( PiBy2*ybar / denom )**2
            else
               u_TwrShadow = 0.0_ReKi
            end if
         else            
            u_TwrShadow = 0.0_ReKi
         end if
                     
         v(1) = (u_TwrPotent + u_TwrShadow)*W_tower
         v(2) = v_TwrPotent*W_tower
         v(3) = 0.0_ReKi
         
         m%DisturbedInflow(:,j,k) = u%InflowOnBlade(:,j,k) + matmul( theta_tower_trans, v ) 
      
      end do !j=NumBlNds
   end do ! NumBlades
   
   
END SUBROUTINE TwrInfl 
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine returns the tower constants necessary to compute the tower influence. 
!! if u%TowerMotion does not have any nodes there will be serious problems. I assume that has been checked earlier.
SUBROUTINE getLocalTowerProps(p, u, BladeNodePosition, theta_tower_trans, W_tower, xbar, ybar, zbar, TwrCd, TwrClrnc, ErrStat, ErrMsg)
!..................................................................................................................................
   TYPE(AD_InputType),           INTENT(IN   )  :: u                       !< Inputs at Time t
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p                       !< Parameters
   REAL(ReKi)                   ,INTENT(IN   )  :: BladeNodePosition(3)    !< local blade node position
   REAL(ReKi)                   ,INTENT(  OUT)  :: theta_tower_trans(3,3)  !< transpose of local tower orientation expressed as a DCM
   REAL(ReKi)                   ,INTENT(  OUT)  :: W_tower                 !< local relative wind speed normal to the tower
   REAL(ReKi)                   ,INTENT(  OUT)  :: xbar                    !< local x^ component of r_TowerBlade normalized by tower radius
   REAL(ReKi)                   ,INTENT(  OUT)  :: ybar                    !< local y^ component of r_TowerBlade normalized by tower radius
   REAL(ReKi)                   ,INTENT(  OUT)  :: zbar                    !< local z^ component of r_TowerBlade normalized by tower radius
   REAL(ReKi)                   ,INTENT(  OUT)  :: TwrCd                   !< local tower drag coefficient
   REAL(ReKi)                   ,INTENT(  OUT)  :: TwrClrnc                !< tower clearance for potential output 
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat                 !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg                  !< Error message if ErrStat /= ErrID_None

   ! local variables
   real(ReKi)                                   :: r_TowerBlade(3)         ! distance vector from tower to blade
   real(ReKi)                                   :: TwrDiam                 ! local tower diameter  
   logical                                      :: found   
   character(*), parameter                      :: RoutineName = 'getLocalTowerProps'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""   
   
   ! ..............................................
   ! option 1: nearest line2 element
   ! ..............................................
   call TwrInfl_NearestLine2Element(p, u, BladeNodePosition, r_TowerBlade, theta_tower_trans, W_tower, xbar, ybar, zbar, TwrCd, TwrDiam, found)
   
   if ( .not. found) then 
      ! ..............................................
      ! option 2: nearest node
      ! ..............................................
      call TwrInfl_NearestPoint(p, u, BladeNodePosition, r_TowerBlade, theta_tower_trans, W_tower, xbar, ybar, zbar, TwrCd, TwrDiam)
         
   end if
   
   TwrClrnc = TwoNorm(r_TowerBlade) - 0.5_ReKi*TwrDiam
   if ( TwrClrnc <= 0.0_ReKi ) then
      call SetErrStat(ErrID_Fatal, "Tower strike.", ErrStat, ErrMsg, RoutineName)
   end if
   
   
END SUBROUTINE getLocalTowerProps
!----------------------------------------------------------------------------------------------------------------------------------
!> Option 1: Find the nearest-neighbor line2 element of the tower mesh for which the blade line2-element node projects orthogonally onto
!!   the tower line2-element domain (following an approach similar to the line2_to_line2 mapping search for motion and scalar quantities). 
!!   That is, for each node of the blade mesh, an orthogonal projection is made onto all possible Line2 elements of the tower mesh and 
!!   the line2 element of the tower mesh that is the minimum distance away is found.
!! Adapted from modmesh_mapping::createmapping_projecttoline2()
SUBROUTINE TwrInfl_NearestLine2Element(p, u, BladeNodePosition, r_TowerBlade, theta_tower_trans, W_tower, xbar, ybar, zbar, TwrCd, TwrDiam, found)
!..................................................................................................................................
   TYPE(AD_InputType),              INTENT(IN   )  :: u                             !< Inputs at Time t
   TYPE(AD_ParameterType),          INTENT(IN   )  :: p                             !< Parameters
   REAL(ReKi)                      ,INTENT(IN   )  :: BladeNodePosition(3)          !< local blade node position
   REAL(ReKi)                      ,INTENT(  OUT)  :: r_TowerBlade(3)               !< distance vector from tower to blade
   REAL(ReKi)                      ,INTENT(  OUT)  :: theta_tower_trans(3,3)        !< transpose of local tower orientation expressed as a DCM
   REAL(ReKi)                      ,INTENT(  OUT)  :: W_tower                       !< local relative wind speed normal to the tower
   REAL(ReKi)                      ,INTENT(  OUT)  :: xbar                          !< local x^ component of r_TowerBlade normalized by tower radius
   REAL(ReKi)                      ,INTENT(  OUT)  :: ybar                          !< local y^ component of r_TowerBlade normalized by tower radius
   REAL(ReKi)                      ,INTENT(  OUT)  :: zbar                          !< local z^ component of r_TowerBlade normalized by tower radius
   REAL(ReKi)                      ,INTENT(  OUT)  :: TwrCd                         !< local tower drag coefficient
   REAL(ReKi)                      ,INTENT(  OUT)  :: TwrDiam                       !< local tower diameter
   logical                         ,INTENT(  OUT)  :: found                         !< whether a mapping was found with this option 
      
      ! local variables
   REAL(ReKi)      :: denom
   REAL(ReKi)      :: dist
   REAL(ReKi)      :: min_dist
   REAL(ReKi)      :: elem_position, elem_position2
   REAL(SiKi)      :: elem_position_SiKi

   REAL(ReKi)      :: p1(3), p2(3)        ! position vectors for nodes on tower line 2 element
   
   REAL(ReKi)      :: V_rel_tower(3)
   
   REAL(ReKi)      :: n1_n2_vector(3)     ! vector going from node 1 to node 2 in Line2 element
   REAL(ReKi)      :: n1_Point_vector(3)  ! vector going from node 1 in Line 2 element to Destination Point
   REAL(ReKi)      :: tmp(3)              ! temporary vector for cross product calculation

   INTEGER(IntKi)  :: jElem               ! do-loop counter for elements on tower mesh

   INTEGER(IntKi)  :: n1, n2              ! nodes associated with an element

   LOGICAL         :: on_element
   
      
   found = .false.
   min_dist = HUGE(min_dist)

   do jElem = 1, u%TowerMotion%ElemTable(ELEMENT_LINE2)%nelem   ! number of elements on TowerMesh
         ! grab node numbers associated with the jElem_th element
      n1 = u%TowerMotion%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(1)
      n2 = u%TowerMotion%ElemTable(ELEMENT_LINE2)%Elements(jElem)%ElemNodes(2)

      p1 = u%TowerMotion%Position(:,n1) + u%TowerMotion%TranslationDisp(:,n1)
      p2 = u%TowerMotion%Position(:,n2) + u%TowerMotion%TranslationDisp(:,n2)

         ! Calculate vectors used in projection operation
      n1_n2_vector    = p2 - p1
      n1_Point_vector = BladeNodePosition - p1

      denom           = DOT_PRODUCT( n1_n2_vector, n1_n2_vector ) ! we've already checked that these aren't zero

         ! project point onto line defined by n1 and n2

      elem_position = DOT_PRODUCT(n1_n2_vector,n1_Point_vector) / denom

            ! note: i forumlated it this way because Fortran doesn't necessarially do shortcutting and I don't want to call EqualRealNos if we don't need it:
      if ( elem_position .ge. 0.0_ReKi .and. elem_position .le. 1.0_ReKi ) then !we're ON the element (between the two nodes)
         on_element = .true.
      else
         elem_position_SiKi = REAL( elem_position, SiKi )
         if (EqualRealNos( elem_position_SiKi, 1.0_SiKi )) then !we're ON the element (at a node)
            on_element = .true.
            elem_position = 1.0_ReKi
         elseif (EqualRealNos( elem_position_SiKi,  0.0_SiKi )) then !we're ON the element (at a node)
            on_element = .true.
            elem_position = 0.0_ReKi
         else !we're not on the element
            on_element = .false.
         end if
         
      end if

      if (on_element) then

         ! calculate distance between point and line (note: this is actually the distance squared);
         ! will only store information once we have determined the closest element
         elem_position2 = 1.0_ReKi - elem_position
         
         r_TowerBlade  = BladeNodePosition - elem_position2*p1 - elem_position*p2
         dist = dot_product( r_TowerBlade, r_TowerBlade )

         if (dist .lt. min_dist) then
            found = .true.
            min_dist = dist

            V_rel_tower =   ( u%InflowOnTower(:,n1) - u%TowerMotion%TranslationVel(:,n1) ) * elem_position2  &
                          + ( u%InflowOnTower(:,n2) - u%TowerMotion%TranslationVel(:,n2) ) * elem_position
            
            TwrDiam     = elem_position2*p%TwrDiam(n1) + elem_position*p%TwrDiam(n2)
            TwrCd       = elem_position2*p%TwrCd(  n1) + elem_position*p%TwrCd(  n2)
            
            
            ! z_hat
            theta_tower_trans(:,3) = n1_n2_vector / sqrt( denom ) ! = n1_n2_vector / twoNorm( n1_n2_vector )
            
            tmp = V_rel_tower - dot_product(V_rel_tower,theta_tower_trans(:,3)) * theta_tower_trans(:,3)
            denom = TwoNorm( tmp )
            if (.not. EqualRealNos( denom, 0.0_ReKi ) ) then
               ! x_hat
               theta_tower_trans(:,1) = tmp / denom
               
               ! y_hat
               tmp = cross_product( theta_tower_trans(:,3), V_rel_tower )
               theta_tower_trans(:,2) = tmp / denom  
               
               W_tower = dot_product( V_rel_tower,theta_tower_trans(:,1) )
               xbar    = 2.0/TwrDiam * dot_product( r_TowerBlade, theta_tower_trans(:,1) )
               ybar    = 2.0/TwrDiam * dot_product( r_TowerBlade, theta_tower_trans(:,2) )
               zbar    = 0.0_ReKi
                                             
            else
                  ! there is no tower influence because dot_product(V_rel_tower,x_hat) = 0
                  ! thus, we don't need to set the other values (except we don't want the sum of xbar^2 and ybar^2 to be 0)
               theta_tower_trans = 0.0_ReKi
               W_tower           = 0.0_ReKi
               xbar              = 1.0_ReKi
               ybar              = 0.0_ReKi  
               zbar              = 0.0_ReKi
            end if
   
            
         end if !the point is closest to this line2 element

      end if

   end do !jElem

END SUBROUTINE TwrInfl_NearestLine2Element
!----------------------------------------------------------------------------------------------------------------------------------
!> Option 2: used when the blade node does not orthogonally intersect a tower element.
!!  Find the nearest-neighbor node in the tower Line2-element domain (following an approach similar to the point_to_point mapping
!!  search for motion and scalar quantities). That is, for each node of the blade mesh, the node of the tower mesh that is the minimum 
!!  distance away is found.
SUBROUTINE TwrInfl_NearestPoint(p, u, BladeNodePosition, r_TowerBlade, theta_tower_trans, W_tower, xbar, ybar, zbar, TwrCd, TwrDiam)
!..................................................................................................................................
   TYPE(AD_InputType),              INTENT(IN   )  :: u                             !< Inputs at Time t
   TYPE(AD_ParameterType),          INTENT(IN   )  :: p                             !< Parameters
   REAL(ReKi)                      ,INTENT(IN   )  :: BladeNodePosition(3)          !< local blade node position
   REAL(ReKi)                      ,INTENT(  OUT)  :: r_TowerBlade(3)               !< distance vector from tower to blade
   REAL(ReKi)                      ,INTENT(  OUT)  :: theta_tower_trans(3,3)        !< transpose of local tower orientation expressed as a DCM
   REAL(ReKi)                      ,INTENT(  OUT)  :: W_tower                       !< local relative wind speed normal to the tower
   REAL(ReKi)                      ,INTENT(  OUT)  :: xbar                          !< local x^ component of r_TowerBlade normalized by tower radius
   REAL(ReKi)                      ,INTENT(  OUT)  :: ybar                          !< local y^ component of r_TowerBlade normalized by tower radius
   REAL(ReKi)                      ,INTENT(  OUT)  :: zbar                          !< local z^ component of r_TowerBlade normalized by tower radius
   REAL(ReKi)                      ,INTENT(  OUT)  :: TwrCd                         !< local tower drag coefficient
   REAL(ReKi)                      ,INTENT(  OUT)  :: TwrDiam                       !< local tower diameter
      
      ! local variables
   REAL(ReKi)      :: denom
   REAL(ReKi)      :: dist
   REAL(ReKi)      :: min_dist
   REAL(ReKi)      :: cosTaper

   REAL(ReKi)      :: p1(3)                     ! position vectors for nodes on tower   
   REAL(ReKi)      :: V_rel_tower(3)
   
   REAL(ReKi)      :: tmp(3)                    ! temporary vector for cross product calculation

   INTEGER(IntKi)  :: n1                        ! node
   INTEGER(IntKi)  :: node_with_min_distance    

   
   
      !.................
      ! find the closest node
      !.................
      
   min_dist = HUGE(min_dist)
   node_with_min_distance = 0

   do n1 = 1, u%TowerMotion%NNodes   ! number of nodes on TowerMesh
      
      p1 = u%TowerMotion%Position(:,n1) + u%TowerMotion%TranslationDisp(:,n1)
      
         ! calculate distance between points (note: this is actually the distance squared);
         ! will only store information once we have determined the closest node
      r_TowerBlade  = BladeNodePosition - p1         
      dist = dot_product( r_TowerBlade, r_TowerBlade )

      if (dist .lt. min_dist) then
         min_dist = dist
         node_with_min_distance = n1
               
      end if !the point is (so far) closest to this blade node

   end do !n1
   
      !.................
      ! calculate the values to be returned:  
      !..................
   if (node_with_min_distance == 0) then
      node_with_min_distance = 1
      if (NWTC_VerboseLevel == NWTC_Verbose) call WrScr( 'AD:TwrInfl_NearestPoint:Error finding minimum distance. Positions may be invalid.' )
   end if
   
   n1 = node_with_min_distance
   
   r_TowerBlade = BladeNodePosition - u%TowerMotion%Position(:,n1) - u%TowerMotion%TranslationDisp(:,n1)
   V_rel_tower  = u%InflowOnTower(:,n1) - u%TowerMotion%TranslationVel(:,n1)
   TwrDiam      = p%TwrDiam(n1) 
   TwrCd        = p%TwrCd(  n1) 
                           
   ! z_hat
   theta_tower_trans(:,3) = u%TowerMotion%Orientation(3,:,n1)
            
   tmp = V_rel_tower - dot_product(V_rel_tower,theta_tower_trans(:,3)) * theta_tower_trans(:,3)
   denom = TwoNorm( tmp )
   
   if (.not. EqualRealNos( denom, 0.0_ReKi ) ) then
      
      ! x_hat
      theta_tower_trans(:,1) = tmp / denom
               
      ! y_hat
      tmp = cross_product( theta_tower_trans(:,3), V_rel_tower )
      theta_tower_trans(:,2) = tmp / denom  
               
      W_tower = dot_product( V_rel_tower,theta_tower_trans(:,1) )

      if ( n1 == 1 .or. n1 == u%TowerMotion%NNodes) then         
         ! option 2b
         zbar    = 2.0/TwrDiam * dot_product( r_TowerBlade, theta_tower_trans(:,3) )
         if (abs(zbar) < 1) then   
            cosTaper = cos( PiBy2*zbar )
            xbar = 2.0/TwrDiam * dot_product( r_TowerBlade, theta_tower_trans(:,1) ) / cosTaper
            ybar = 2.0/TwrDiam * dot_product( r_TowerBlade, theta_tower_trans(:,2) ) / cosTaper
         else ! we check that zbar < 1 before using xbar and ybar later, but I'm going to set them here anyway:
            xbar = 1.0_ReKi
            ybar = 0.0_ReKi  
         end if                                    
      else
         ! option 2a
         xbar    = 2.0/TwrDiam * dot_product( r_TowerBlade, theta_tower_trans(:,1) )
         ybar    = 2.0/TwrDiam * dot_product( r_TowerBlade, theta_tower_trans(:,2) )
         zbar    = 0.0_ReKi
      end if

   else
      
         ! there is no tower influence because W_tower = dot_product(V_rel_tower,x_hat) = 0
         ! thus, we don't need to set the other values (except we don't want the sum of xbar^2 and ybar^2 to be 0)
      W_tower           = 0.0_ReKi
      theta_tower_trans = 0.0_ReKi
      xbar              = 1.0_ReKi
      ybar              = 0.0_ReKi  
      zbar              = 0.0_ReKi
      
   end if   

END SUBROUTINE TwrInfl_NearestPoint
!----------------------------------------------------------------------------------------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in AD_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and dZ/du are returned.
SUBROUTINE AD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu)
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(AD_InputType),                   INTENT(INOUT)           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(AD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(AD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(AD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(AD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(AD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(AD_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is
                                                                               !!   available here so that mesh parameter information (i.e.,
                                                                               !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(AD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions (Y) with respect
                                                                               !!   to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) with
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) with
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z) with
                                                                               !!   respect to the inputs (u) [intent in to avoid deallocation]
      ! local variables
   TYPE(AD_OutputType)                                           :: y_p
   TYPE(AD_OutputType)                                           :: y_m
   TYPE(AD_ConstraintStateType)                                  :: z_p
   TYPE(AD_ConstraintStateType)                                  :: z_m
   TYPE(AD_InputType)                                            :: u_perturb
   REAL(R8Ki)                                                    :: delta_p, delta_m  ! delta change in input
   INTEGER(IntKi)                                                :: i, j, k, n   
   logical                                                       :: ValidInput
   
   integer, parameter                                            :: indx = 1      ! m%BEMT_u(1) is at t; m%BEMT_u(2) is t+dt
   integer, parameter                                            :: op_indx = 2   ! m%BEMT_u(1) is at t; m%BEMT_u(2) is t+dt or the input at OP
   integer(intKi)                                                :: ErrStat2
   character(ErrMsgLen)                                          :: ErrMsg2
   character(*), parameter                                       :: RoutineName = 'AD_JacobianPInput'


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   
      ! get OP values here:
   !call AD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
   call SetInputsForBEMT(p, u, m, indx, errStat2, errMsg2)  
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
   call BEMT_CopyInput( m%BEMT_u(indx), m%BEMT_u(op_indx), MESH_UPDATECOPY, ErrStat2, ErrMsg2) ! copy the BEMT OP inputs to a temporary location that won't be overwritten
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later                        
   
   
   if ( p%FrozenWake ) then            
            ! compare arguments with call to BEMT_CalcOutput   
      call computeFrozenWake(m%BEMT_u(op_indx), p%BEMT, m%BEMT_y, m%BEMT )            
      m%BEMT%UseFrozenWake = .true.
   end if
         
   
      ! make a copy of the inputs to perturb
   call AD_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
   

   IF ( PRESENT( dYdu ) ) THEN
      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
      
      ! allocate dYdu
      if (.not. allocated(dYdu) ) then
         call AllocAry(dYdu,p%Jac_ny, size(p%Jac_u_indx,1),'dYdu', ErrStat2, ErrMsg2)
         call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if
   
      
         ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call AD_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AD_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
         
      do i=1,size(p%Jac_u_indx,1)
         
            ! get u_op + delta_p u
         call AD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         call Perturb_u( p, i, 1, u_perturb, delta_p )

            ! we need to see if these perturbed inputs put us in different solution regions:
         call SetInputsForBEMT(p, u_perturb, m, indx, errStat2, errMsg2)  
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         ValidInput = CheckBEMTInputPerturbations( p, m )
      
            ! if so, we do a 1-sided difference:
         if (.not. ValidInput) then
            call AD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
            delta_p = 0
         end if
               
         
            ! compute y at u_op + delta_p u
         call AD_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get u_op - delta_m u
         call AD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later
         call Perturb_u( p, i, -1, u_perturb, delta_m )
         
            ! we need to see if these perturbed inputs put us in different solution regions:
         call SetInputsForBEMT(p, u_perturb, m, indx, errStat2, errMsg2)  
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         ValidInput = CheckBEMTInputPerturbations( p, m )
      
            ! if so, we do a 1-sided difference:
         if (.not. ValidInput) then
            call AD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
            delta_m = 0
            if (EqualRealNos(delta_p, 0.0_R8Ki)) then
               call SetErrStat(ErrID_Fatal,'Both sides of central difference equation change solution region. '// &
                  'dYdu cannot be calculated for column '//trim(num2lstr(i))//'.',ErrStat,ErrMsg,RoutineName) 
               return
            end if
         end if         
         
         
            ! compute y at u_op - delta_m u
         call AD_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         
            
            ! get central difference:            
         call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdu(:,i) )
         
      end do
      

      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      call AD_DestroyOutput( y_p, ErrStat2, ErrMsg2 ) ! we don't need this any more   
      call AD_DestroyOutput( y_m, ErrStat2, ErrMsg2 ) ! we don't need this any more   
   
      
   END IF

   IF ( PRESENT( dXdu ) ) THEN
      if (allocated(dXdu)) deallocate(dXdu)
   END IF

   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu)
   END IF

   IF ( PRESENT( dZdu ) ) THEN

      call CheckLinearizationInput(p%BEMT, m%BEMT_u(op_indx), z%BEMT, m%BEMT, OtherState%BEMT, ErrStat2, ErrMsg2)      
         call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if         
      
      
      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the inputs (u) here:

      ! allocate dZdu
      if (.not. allocated(dZdu)) then
         call AllocAry(dZdu,size(z%BEMT%phi), size(p%Jac_u_indx,1),'dZdu', ErrStat2, ErrMsg2)
         call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if         
      end if

      
      do i=1,size(p%Jac_u_indx,1)
         
            ! get u_op + delta_p u
         call AD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         call Perturb_u( p, i, 1, u_perturb, delta_p )

            ! we need to see if these perturbed inputs put us in different solution regions:
         call SetInputsForBEMT(p, u_perturb, m, indx, errStat2, errMsg2)  
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         ValidInput = CheckBEMTInputPerturbations( p, m )
      
            ! if so, we do a 1-sided difference:
         if (.not. ValidInput) then
            call AD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
            delta_p = 0
         end if
         
                  
            ! compute z_p at u_op + delta_p u
         call AD_CalcConstrStateResidual( t, u_perturb, p, x, xd, z, OtherState, m, z_p, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
                                         
            ! get u_op - delta_m u
         call AD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)             
         call Perturb_u( p, i, -1, u_perturb, delta_m )
         
            ! we need to see if these perturbed inputs put us in different solution regions:
         call SetInputsForBEMT(p, u_perturb, m, indx, errStat2, errMsg2)  
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
         ValidInput = CheckBEMTInputPerturbations( p, m )
      
            ! if so, we do a 1-sided difference:
         if (.not. ValidInput) then
            call AD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
               call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
            delta_m = 0
            if (EqualRealNos(delta_p, 0.0_R8Ki)) then
               call SetErrStat(ErrID_Fatal,'Both sides of central difference equation change solution region. '// &
                  'dYdu cannot be calculated for column '//trim(num2lstr(i))//'.',ErrStat,ErrMsg,RoutineName) 
               return
            end if
         end if         
         
                  
            ! compute z_m at u_op - delta_m u
         call AD_CalcConstrStateResidual( t, u_perturb, p, x, xd, z, OtherState, m, z_m, ErrStat2, ErrMsg2 ) 
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
            
            
            ! get central difference:            
            
            ! we may have had an error allocating memory, so we'll check
         if (ErrStat>=AbortErrLev) then 
            call cleanup()
            return
         end if         
         
         
         do k=1,p%NumBlades ! size(z%BEMT%Phi,2)
            do j=1,p%NumBlNds ! size(z%BEMT%Phi,1)
               n = (k-1)*p%NumBlNds + j
               dZdu(n,i) = z_p%BEMT%Phi(j,k) - z_m%BEMT%Phi(j,k)
            end do            
         end do
         
         dZdu(:,i) = dZdu(:,i) / (delta_p + delta_m) 
         
      end do
      
      call AD_DestroyConstrState( z_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call AD_DestroyConstrState( z_m, ErrStat2, ErrMsg2 ) ! we don't need this any more      
      
   END IF
   
   call cleanup()
contains
   subroutine cleanup()
      m%BEMT%UseFrozenWake = .false.
   
      call AD_DestroyOutput(      y_p, ErrStat2, ErrMsg2 )
      call AD_DestroyOutput(      y_m, ErrStat2, ErrMsg2 )
      call AD_DestroyConstrState( z_p, ErrStat2, ErrMsg2 )
      call AD_DestroyConstrState( z_m, ErrStat2, ErrMsg2 )
      call AD_DestroyInput( u_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup

END SUBROUTINE AD_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
SUBROUTINE AD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(AD_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(AD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(AD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(AD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(AD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(AD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(AD_OutputType),                  INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is
                                                                               !!   available here so that mesh parameter information (i.e.,
                                                                               !!   connectivity) does not have to be recalculated for dYdx.
   TYPE(AD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdx(:,:)  !< Partial derivatives of output functions
                                                                               !!   (Y) with respect to the continuous
                                                                               !!   states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdx(:,:)  !< Partial derivatives of continuous state
                                                                               !!   functions (X) with respect to
                                                                               !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddx(:,:) !< Partial derivatives of discrete state
                                                                               !!   functions (Xd) with respect to
                                                                               !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdx(:,:)  !< Partial derivatives of constraint state
                                                                               !!   functions (Z) with respect to
                                                                               !!   the continuous states (x) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''



   IF ( PRESENT( dYdx ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:

      ! allocate and set dYdx

   END IF

   IF ( PRESENT( dXdx ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:

      ! allocate and set dXdx

   END IF

   IF ( PRESENT( dXddx ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the continuous states (x) here:

      ! allocate and set dXddx

   END IF

   IF ( PRESENT( dZdx ) ) THEN


      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the continuous states (x) here:

      ! allocate and set dZdx

   END IF


END SUBROUTINE AD_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and dZ/dxd are returned.
SUBROUTINE AD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(AD_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(AD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(AD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(AD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(AD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(AD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(AD_OutputType),                  INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is
                                                                               !!   available here so that mesh parameter information (i.e.,
                                                                               !!   connectivity) does not have to be recalculated for dYdxd.
   TYPE(AD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdxd(:,:) !< Partial derivatives of output functions
                                                                               !!  (Y) with respect to the discrete
                                                                               !!  states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdxd(:,:) !< Partial derivatives of continuous state
                                                                               !!   functions (X) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddxd(:,:)!< Partial derivatives of discrete state
                                                                               !!   functions (Xd) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdxd(:,:) !< Partial derivatives of constraint state
                                                                               !!   functions (Z) with respect to the
                                                                               !!   discrete states (xd) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''


   IF ( PRESENT( dYdxd ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the discrete states (xd) here:

      ! allocate and set dYdxd

   END IF

   IF ( PRESENT( dXdxd ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the discrete states (xd) here:

      ! allocate and set dXdxd

   END IF

   IF ( PRESENT( dXddxd ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the discrete states (xd) here:

      ! allocate and set dXddxd

   END IF

   IF ( PRESENT( dZdxd ) ) THEN

      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the discrete states (xd) here:

      ! allocate and set dZdxd

   END IF


END SUBROUTINE AD_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and dZ/dz are returned.
SUBROUTINE AD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
!..................................................................................................................................

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(AD_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(AD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(AD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(AD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(AD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(AD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(AD_OutputType),                  INTENT(INOUT)           :: y          !< Output (change to inout if a mesh copy is required);
                                                                               !!   Output fields are not used by this routine, but type is
                                                                               !!   available here so that mesh parameter information (i.e.,
                                                                               !!   connectivity) does not have to be recalculated for dYdz.
   TYPE(AD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dYdz(:,:)  !< Partial derivatives of output
                                                                               !!  functions (Y) with respect to the
                                                                               !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXdz(:,:)  !< Partial derivatives of continuous
                                                                               !!  state functions (X) with respect to
                                                                               !!  the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dXddz(:,:) !< Partial derivatives of discrete state
                                                                               !!  functions (Xd) with respect to the
                                                                               !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dZdz(:,:)  !< Partial derivatives of constraint
                                                                               !! state functions (Z) with respect to
                                                                               !!  the constraint states (z) [intent in to avoid deallocation]

      ! local variables
   TYPE(AD_OutputType)                                           :: y_p
   TYPE(AD_OutputType)                                           :: y_m
   TYPE(AD_ConstraintStateType)                                  :: Z_p
   TYPE(AD_ConstraintStateType)                                  :: Z_m
   TYPE(AD_ConstraintStateType)                                  :: z_perturb
   REAL(R8Ki)                                                    :: delta_p, delta_m  ! delta change in state
   INTEGER(IntKi)                                                :: i, j, k, n, k2, j2   

   integer, parameter                                            :: indx = 1      ! m%BEMT_u(1) is at t; m%BEMT_u(2) is t+dt
   integer, parameter                                            :: op_indx = 2   ! m%BEMT_u(1) is at t; m%BEMT_u(2) is t+dt or the input at OP
   integer(intKi)                                                :: ErrStat2
   character(ErrMsgLen)                                          :: ErrMsg2
   character(*), parameter                                       :: RoutineName = 'AD_JacobianPConstrState'

   
      ! local variables
      
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''
   
!bjj: how do I figure out if F is 0??? In that case, need to se dY/dz = 0 and dZ/dz = 1 {and need to ask jmj if this is the whole matrix or just a row/column where it applies}   

      ! get OP values here:   
   !call AD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )  ! (bjj: is this necessary? if not, still need to get BEMT inputs)
   call SetInputsForBEMT(p, u, m, indx, errStat2, errMsg2)  
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
   call BEMT_CopyInput( m%BEMT_u(indx), m%BEMT_u(op_indx), MESH_UPDATECOPY, ErrStat2, ErrMsg2) ! copy the BEMT OP inputs to a temporary location that won't be overwritten
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later                        
 
      
   if ( p%FrozenWake ) then            
            ! compare arguments with call to BEMT_CalcOutput   
      call computeFrozenWake(m%BEMT_u(op_indx), p%BEMT, m%BEMT_y, m%BEMT )      
      m%BEMT%UseFrozenWake = .true.
   end if
   
   
      ! make a copy of the constraint states to perturb
   call AD_CopyConstrState( z, z_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
   
   
   IF ( PRESENT( dYdz ) ) THEN

         ! Calculate the partial derivative of the output functions (Y) with respect to the constraint states (z) here:

      ! allocate and set dYdz
      if (.not. allocated(dYdz) ) then
         call AllocAry(dYdz,p%Jac_ny, size(z%BEMT%phi),'dYdz', ErrStat2, ErrMsg2)
         call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      end if

      
         ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call AD_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AD_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if
      
         
      do k=1,p%NumBlades ! size(z%BEMT%Phi,2)
         do j=1,p%NumBlNds ! size(z%BEMT%Phi,1)                  
            i = (k-1)*p%NumBlNds + j
            
               ! need a check if F = 0 for this case:
   
            if ( ( p%BEMT%UseTipLoss .and. EqualRealNos(p%BEMT%tipLossConst(j,k),0.0_ReKi) ) .or. &
                 ( p%BEMT%useHubLoss .and. EqualRealNos(p%BEMT%hubLossConst(j,k),0.0_ReKi) ) ) then
               ! F is zero, we we need to skip this perturbation
               dYdz(:,i) = 0.0_ReKi
            else                        
            
               call Get_phi_perturbations(p%BEMT, m%BEMT, z%BEMT%phi(j,k), delta_p, delta_m)
               
                  ! get z_op + delta_p z
               z_perturb%BEMT%phi(j,k) = z%BEMT%phi(j,k) + delta_p
            
                  ! compute y at z_op + delta_p z
               call AD_CalcOutput( t, u, p, x, xd, z_perturb, OtherState, y_p, m, ErrStat2, ErrMsg2 ) 
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
            
            
                  ! get z_op - delta_m z
               z_perturb%BEMT%phi(j,k) = z%BEMT%phi(j,k) - delta_m
            
                  ! compute y at z_op - delta_m z
               call AD_CalcOutput( t, u, p, x, xd, z_perturb, OtherState, y_m, m, ErrStat2, ErrMsg2 ) 
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) ! we shouldn't have any errors about allocating memory here so I'm not going to return-on-error until later            
            

                  ! get central difference:            
               call Compute_dY( p, y_p, y_m, delta_p, delta_m, dYdz(:,i) )
               
               
                  ! put z_perturb back (for next iteration):
               z_perturb%BEMT%phi(j,k) = z%BEMT%phi(j,k)
            end if
         
         end do
      end do
      
      if (ErrStat>=AbortErrLev) then
         call cleanup()
         return
      end if
      call AD_DestroyOutput( y_p, ErrStat2, ErrMsg2 ) ! we don't need this any more   
      call AD_DestroyOutput( y_m, ErrStat2, ErrMsg2 ) ! we don't need this any more   
      
      
   END IF

   IF ( PRESENT( dXdz ) ) THEN
      if (allocated(dXdz)) deallocate(dXdz)
   END IF

   IF ( PRESENT( dXddz ) ) THEN
      if (allocated(dXddz)) deallocate(dXddz)
   END IF

   IF ( PRESENT(dZdz) ) THEN

      call CheckLinearizationInput(p%BEMT, m%BEMT_u(op_indx), z%BEMT, m%BEMT, OtherState%BEMT, ErrStat2, ErrMsg2)      
         call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if         
         
         ! Calculate the partial derivative of the constraint state functions (Z) with respect to the constraint states (z) here:

      ! allocate and set dZdz
      if (.not. allocated(dZdz)) then
         call AllocAry(dZdz,size(z%BEMT%phi), size(z%BEMT%phi),'dZdz', ErrStat2, ErrMsg2)
         call setErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if (ErrStat>=AbortErrLev) then
            call cleanup()
            return
         end if         
      end if
      
      
      call AD_CopyConstrState( z, z_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
      
      do k=1,p%NumBlades ! size(z%BEMT%Phi,2)
         do j=1,p%NumBlNds ! size(z%BEMT%Phi,1)                  
            i = (k-1)*p%NumBlNds + j
               
            if ( ( p%BEMT%UseTipLoss .and. EqualRealNos(p%BEMT%tipLossConst(j,k),0.0_ReKi) ) .or. &
                 ( p%BEMT%useHubLoss .and. EqualRealNos(p%BEMT%hubLossConst(j,k),0.0_ReKi) ) ) then
               ! F is zero, we we need to skip this perturbation
               dZdz(:,i) = 0.0_ReKi
               dZdz(i,i) = 1.0_ReKi                              
            else                        
            
               call Get_phi_perturbations(p%BEMT, m%BEMT, z%BEMT%phi(j,k), delta_p, delta_m)
            
                  ! get z_op + delta_p z
               z_perturb%BEMT%phi(j,k) = z%BEMT%phi(j,k) + delta_p

                  ! compute z_p at z_op + delta_p z
               call AD_CalcConstrStateResidual( t, u, p, x, xd, z_perturb, OtherState, m, z_p, ErrStat2, ErrMsg2 ) 
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
                                         
                  ! get z_op - delta_m z
               z_perturb%BEMT%phi(j,k) = z%BEMT%phi(j,k) - delta_m
                     
                  ! compute z_m at u_op - delta_m u
               call AD_CalcConstrStateResidual( t, u, p, x, xd, z_perturb, OtherState, m, z_m, ErrStat2, ErrMsg2 ) 
                  call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
                  if (ErrStat>=AbortErrLev) then 
                     call cleanup()
                     return
                  end if         
            
                  ! get central difference:            
                     
               do k2=1,p%NumBlades ! size(z%BEMT%Phi,2)
                  do j2=1,p%NumBlNds ! size(z%BEMT%Phi,1)
                     n = (k2-1)*p%NumBlNds + j2
                     dZdz(n,i) = z_p%BEMT%Phi(j2,k2) - z_m%BEMT%Phi(j2,k2)
                  end do            
               end do
         
               dZdz(:,i) = dZdz(:,i) / (delta_p + delta_m) 
         
                  ! put z_perturb back (for next iteration):
               z_perturb%BEMT%phi(j,k) = z%BEMT%phi(j,k)
               
            end if
            
         end do         
      end do
      
      call AD_DestroyConstrState( z_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call AD_DestroyConstrState( z_m, ErrStat2, ErrMsg2 ) ! we don't need this any more      
      
   END IF
     
   call cleanup()
   
contains
   subroutine cleanup()
      m%BEMT%UseFrozenWake = .false.

      call AD_DestroyOutput(            y_p, ErrStat2, ErrMsg2 )
      call AD_DestroyOutput(            y_m, ErrStat2, ErrMsg2 )
      call AD_DestroyConstrState(       z_p, ErrStat2, ErrMsg2 )
      call AD_DestroyConstrState(       z_m, ErrStat2, ErrMsg2 )
      call AD_DestroyConstrState( z_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup   

END SUBROUTINE AD_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE AD_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(AD_InputType),                   INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(AD_ParameterType),               INTENT(IN   )           :: p          !< Parameters
   TYPE(AD_ContinuousStateType),         INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(AD_DiscreteStateType),           INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(AD_ConstraintStateType),         INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(AD_OtherStateType),              INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(AD_OutputType),                  INTENT(IN   )           :: y          !< Output at operating point
   TYPE(AD_MiscVarType),                 INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: z_op(:)    !< values of linearized constraint states

   INTEGER(IntKi)                                                :: index, i, j, k
   INTEGER(IntKi)                                                :: nu
   INTEGER(IntKi)                                                :: ErrStat2
   CHARACTER(ErrMsgLen)                                          :: ErrMsg2
   CHARACTER(*), PARAMETER                                       :: RoutineName = 'AD_GetOP'
   LOGICAL                                                       :: FieldMask(FIELDMASK_SIZE)

   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( u_op ) ) THEN
      
      nu = size(p%Jac_u_indx,1) + u%TowerMotion%NNodes * 6 & ! Jac_u_indx has 3 orientation angles, but the OP needs the full 9 elements of the DCM
                                + u%hubMotion%NNodes * 6     ! Jac_u_indx has 3 orientation angles, but the OP needs the full 9 elements of the DCM
      do i=1,p%NumBlades
         nu = nu + u%BladeMotion(i)%NNodes * 6 & ! Jac_u_indx has 3 orientation angles, but the OP needs the full 9 elements of the DCM
             + u%BladeRootMotion(i)%NNodes * 6   ! Jac_u_indx has 3 orientation angles, but the OP needs the full 9 elements of the DCM
      end do      
                  
      if (.not. allocated(u_op)) then
         call AllocAry(u_op, nu, 'u_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      

      index = 1
      FieldMask = .false.
      FieldMask(MASKID_TRANSLATIONDISP) = .true.
      FieldMask(MASKID_Orientation) = .true.
      FieldMask(MASKID_TRANSLATIONVel) = .true.
      call PackMotionMesh(u%TowerMotion, u_op, index, FieldMask=FieldMask)
   
      FieldMask(MASKID_TRANSLATIONVel) = .false.
      FieldMask(MASKID_RotationVel) = .true.
      call PackMotionMesh(u%HubMotion, u_op, index, FieldMask=FieldMask)
   
      FieldMask = .false.
      FieldMask(MASKID_Orientation) = .true.
      do k = 1,p%NumBlades
         call PackMotionMesh(u%BladeRootMotion(k), u_op, index, FieldMask=FieldMask)
      end do
   
      FieldMask(MASKID_TRANSLATIONDISP) = .true.
      FieldMask(MASKID_TRANSLATIONVel)  = .true.
      do k=1,p%NumBlades     
         call PackMotionMesh(u%BladeMotion(k), u_op, index, FieldMask=FieldMask)
      end do
   
      do k=1,p%NumBlades
         do i=1,p%NumBlNds
            do j=1,3
               u_op(index) = u%InflowOnBlade(j,i,k)
               index = index + 1
            end do            
         end do
      end do

      do i=1,p%NumTwrNds
         do j=1,3
            u_op(index) = u%InflowOnTower(j,i)
            index = index + 1
         end do            
      end do
      
   END IF

   IF ( PRESENT( y_op ) ) THEN
      
      if (.not. allocated(y_op)) then
         call AllocAry(y_op, p%Jac_ny, 'y_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      
      
      
      index = 1               
      call PackLoadMesh(y%TowerLoad, y_op, index)
      do k=1,p%NumBlades
         call PackLoadMesh(y%BladeLoad(k), y_op, index)                  
      end do
   
      index = index - 1
      do i=1,p%NumOuts
         y_op(i+index) = y%WriteOutput(i)
      end do   
         
      
   END IF

   IF ( PRESENT( x_op ) ) THEN

   END IF

   IF ( PRESENT( dx_op ) ) THEN

   END IF

   IF ( PRESENT( xd_op ) ) THEN

   END IF
   
   IF ( PRESENT( z_op ) ) THEN

      if (.not. allocated(z_op)) then
         call AllocAry(z_op, p%NumBlades*p%NumBlNds, 'z_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      
   
      index = 1      
      do k=1,p%NumBlades ! size(z%BEMT%Phi,2)
         do i=1,p%NumBlNds ! size(z%BEMT%Phi,1)
            z_op(index) = z%BEMT%phi(i,k)
            index = index + 1
         end do            
      end do
      
   END IF

END SUBROUTINE AD_GetOP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
SUBROUTINE Init_Jacobian_y( p, y, InitOut, ErrStat, ErrMsg)

   TYPE(AD_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(AD_OutputType)               , INTENT(IN   ) :: y                     !< outputs
   TYPE(AD_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Initialization output data (for Jacobian row/column names)
   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
      ! local variables:
   INTEGER(IntKi)                :: i, j, k, indx_next, indx_last
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'Init_Jacobian_y'
   logical, allocatable                              :: AllOut(:)
                        
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! determine how many outputs there are in the Jacobians     
   p%Jac_ny = y%TowerLoad%NNodes * 6  & ! 3 forces + 3 moments at each node
            + p%NumOuts                 ! WriteOutput values 
      
   do k=1,p%NumBlades
      p%Jac_ny = p%Jac_ny + y%BladeLoad(k)%NNodes * 6  ! 3 forces + 3 moments at each node
   end do   
   
   
      ! get the names of the linearized outputs:
   call AllocAry(InitOut%LinNames_y, p%Jac_ny,'LinNames_y',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(InitOut%RotFrame_y, p%Jac_ny,'RotFrame_y',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if (ErrStat >= AbortErrLev) return
   
         
   InitOut%RotFrame_y = .false. ! default all to false, then set the true ones below (note that meshes are in the global, not rotating frame)
   indx_next = 1  
   call PackLoadMesh_Names(y%TowerLoad, 'Tower', InitOut%LinNames_y, indx_next)
   
   indx_last = indx_next
   do k=1,p%NumBlades
      call PackLoadMesh_Names(y%BladeLoad(k), 'Blade '//trim(num2lstr(k)), InitOut%LinNames_y, indx_next)
   end do
   ! InitOut%RotFrame_y(indx_last:indx_next-1) = .true. ! The mesh fields are in the global frame, so are not in the rotating frame

   do i=1,p%NumOuts
      InitOut%LinNames_y(i+indx_next-1) = trim(InitOut%WriteOutputHdr(i))//', '//trim(InitOut%WriteOutputUnt(i))  !trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
   end do    
   
      ! check for all the WriteOutput values that are functions of blade number:
   allocate( AllOut(0:MaxOutPts), STAT=ErrStat2 ) ! allocate starting at zero to account for invalid output channels
   if (ErrStat2 /=0 ) then
      call SetErrStat(ErrID_Info, 'error allocating temporary space for AllOut',ErrStat,ErrMsg,RoutineName)
      return;
   end if
   
   AllOut = .false.
   do k=1,3
      AllOut( BAzimuth(k)) = .true.
      AllOut( BPitch  (k)) = .true.
      do j=1,9
         AllOut(BNVUndx(j,k)) = .true.
         AllOut(BNVUndy(j,k)) = .true.
         AllOut(BNVUndz(j,k)) = .true.
         AllOut(BNVDisx(j,k)) = .true.
         AllOut(BNVDisy(j,k)) = .true.
         AllOut(BNVDisz(j,k)) = .true.
         AllOut(BNSTVx (j,k)) = .true.
         AllOut(BNSTVy (j,k)) = .true.
         AllOut(BNSTVz (j,k)) = .true.
         AllOut(BNVRel (j,k)) = .true.
         AllOut(BNDynP (j,k)) = .true.
         AllOut(BNRe   (j,k)) = .true.
         AllOut(BNM    (j,k)) = .true.   
         AllOut(BNVIndx(j,k)) = .true.   
         AllOut(BNVIndy(j,k)) = .true. 
         AllOut(BNAxInd(j,k)) = .true.         
         AllOut(BNTnInd(j,k)) = .true.
         AllOut(BNAlpha(j,k)) = .true.
         AllOut(BNTheta(j,k)) = .true.
         AllOut(BNPhi  (j,k)) = .true.   
         AllOut(BNCurve(j,k)) = .true.
         AllOut(BNCl   (j,k)) = .true.
         AllOut(BNCd   (j,k)) = .true.
         AllOut(BNCm   (j,k)) = .true.
         AllOut(BNCx   (j,k)) = .true.
         AllOut(BNCy   (j,k)) = .true.
         AllOut(BNCn   (j,k)) = .true.
         AllOut(BNCt   (j,k)) = .true.
         AllOut(BNFl   (j,k)) = .true.
         AllOut(BNFd   (j,k)) = .true.
         AllOut(BNMm   (j,k)) = .true.
         AllOut(BNFx   (j,k)) = .true.
         AllOut(BNFy   (j,k)) = .true.
         AllOut(BNFn   (j,k)) = .true.
         AllOut(BNFt   (j,k)) = .true.
         AllOut(BNClrnc(j,k)) = .true.
      end do
   end do
   
   
   do i=1,p%NumOuts
      InitOut%RotFrame_y(i+indx_next-1) = AllOut( p%OutParam(i)%Indx )      
   end do    
   
   deallocate(AllOut)
          
END SUBROUTINE Init_Jacobian_y
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine initializes the array that maps rows/columns of the Jacobian to specific mesh fields.
!! Do not change the order of this packing without changing corresponding parts of AD linearization !
SUBROUTINE Init_Jacobian( InputFileData, p, u, y, m, InitOut, ErrStat, ErrMsg)

   type(AD_InputFile)                , intent(in   ) :: InputFileData         !< input file data (for default blade perturbation)
   TYPE(AD_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(AD_InputType)                , INTENT(IN   ) :: u                     !< inputs
   TYPE(AD_OutputType)               , INTENT(IN   ) :: y                     !< outputs
   TYPE(AD_MiscVarType)              , INTENT(IN   ) :: m                     !< miscellaneous variable
   TYPE(AD_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Initialization output data (for Jacobian row/column names)
   
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'Init_Jacobian'
   
      ! local variables:
   INTEGER(IntKi)                :: i, j, k, index, index_last, nu, i_meshField
   REAL(ReKi)                    :: perturb, perturb_t, perturb_b(MaxBl)
   LOGICAL                       :: FieldMask(FIELDMASK_SIZE)
   CHARACTER(1), PARAMETER       :: UVW(3) = (/'U','V','W'/)
   
            
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   call Init_Jacobian_y( p, y, InitOut, ErrStat, ErrMsg)
   
      ! these matrices will be needed for linearization with frozen wake feature
   if (p%FrozenWake) then
      call AllocAry(m%BEMT%AxInd_op,p%NumBlNds,p%numBlades,'m%BEMT%AxInd_op', ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call AllocAry(m%BEMT%TnInd_op,p%NumBlNds,p%numBlades,'m%BEMT%TnInd_op', ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   end if
   
   
      
      ! determine how many inputs there are in the Jacobians
   nu = u%TowerMotion%NNodes * 9            & ! 3 Translation Displacements + 3 orientations + 3 Translation velocities at each node
      + u%hubMotion%NNodes   * 9            & ! 3 Translation Displacements + 3 orientations + 3 Rotation velocities at each node
      + size( u%InflowOnBlade)              &
      + size( u%InflowOnTower)

   do i=1,p%NumBlades
      nu = nu + u%BladeMotion(i)%NNodes * 9 & ! 3 Translation Displacements + 3 orientations + 3 Translation velocities at each node
          + u%BladeRootMotion(i)%NNodes * 3   ! 3 orientations at each node
   end do      
      
   ! all other inputs ignored

      
   !............................                     
   ! fill matrix to store index to help us figure out what the ith value of the u vector really means
   ! (see aerodyn::perturb_u ... these MUST match )
   ! column 1 indicates module's mesh and field
   ! column 2 indicates the first index (x-y-z component) of the field
   ! column 3 is the node
   !............................                     
   
   call allocAry( p%Jac_u_indx, nu, 3, 'p%Jac_u_indx', ErrStat2, ErrMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   if (ErrStat >= AbortErrLev) return                     
            
   !...............
   ! AD input mappings stored in p%Jac_u_indx:   
   !...............            
   index = 1
   !Module/Mesh/Field: u%TowerMotion%TranslationDisp  = 1;
   !Module/Mesh/Field: u%TowerMotion%Orientation      = 2;
   !Module/Mesh/Field: u%TowerMotion%TranslationVel   = 3;
   do i_meshField = 1,3
      do i=1,u%TowerMotion%NNodes
         do j=1,3
            p%Jac_u_indx(index,1) =  i_meshField
            p%Jac_u_indx(index,2) =  j !component index:  j
            p%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
   end do
   
   !Module/Mesh/Field: u%HubMotion%TranslationDisp = 4;
   !Module/Mesh/Field: u%HubMotion%Orientation     = 5;
   !Module/Mesh/Field: u%HubMotion%RotationVel     = 6;
   do i_meshField = 4,6
      do i=1,u%HubMotion%NNodes
         do j=1,3
            p%Jac_u_indx(index,1) =  i_meshField
            p%Jac_u_indx(index,2) =  j !component index:  j
            p%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
   end do
   
   !bjj: if MaxBl (max blades) changes, we need to modify this
   !Module/Mesh/Field: u%BladeRootMotion(1)%Orientation = 7;
   !Module/Mesh/Field: u%BladeRootMotion(2)%Orientation = 8;
   !Module/Mesh/Field: u%BladeRootMotion(3)%Orientation = 9;   
   do k=1,p%NumBlades         
      do i_meshField = 6,6
         do i=1,u%BladeRootMotion(k)%NNodes
            do j=1,3
               p%Jac_u_indx(index,1) =  i_meshField + k
               p%Jac_u_indx(index,2) =  j !component index:  j
               p%Jac_u_indx(index,3) =  i !Node:   i
               index = index + 1
            end do !j      
         end do !i
            
      end do !i_meshField                            
   end do !k  
      
   !bjj: if MaxBl (max blades) changes, we need to modify this
   !Module/Mesh/Field: u%BladeMotion(1)%TranslationDisp = 10;
   !Module/Mesh/Field: u%BladeMotion(1)%Orientation     = 11;
   !Module/Mesh/Field: u%BladeMotion(1)%TranslationVel  = 12;
   !Module/Mesh/Field: u%BladeMotion(2)%TranslationDisp = 13;
   !Module/Mesh/Field: u%BladeMotion(2)%Orientation     = 14;
   !Module/Mesh/Field: u%BladeMotion(2)%TranslationVel  = 15;
   !Module/Mesh/Field: u%BladeMotion(3)%TranslationDisp = 16;
   !Module/Mesh/Field: u%BladeMotion(3)%Orientation     = 17;
   !Module/Mesh/Field: u%BladeMotion(3)%TranslationVel  = 18;      
   do k=1,p%NumBlades         
      do i_meshField = 1,3
         do i=1,u%BladeMotion(k)%NNodes
            do j=1,3
               p%Jac_u_indx(index,1) =  9 + i_meshField + (k-1)*3
               p%Jac_u_indx(index,2) =  j !component index:  j
               p%Jac_u_indx(index,3) =  i !Node:   i
               index = index + 1
            end do !j      
         end do !i
            
      end do !i_meshField                            
   end do !k
   
   !Module/Mesh/Field: u%InflowOnBlade(:,:,1) = 19;
   !Module/Mesh/Field: u%InflowOnBlade(:,:,2) = 20;
   !Module/Mesh/Field: u%InflowOnBlade(:,:,3) = 21;   
   do k=1,size(u%InflowOnBlade,3) ! p%NumBlades         
      do i=1,size(u%InflowOnBlade,2) ! numNodes
         do j=1,3
            p%Jac_u_indx(index,1) =  18 + k
            p%Jac_u_indx(index,2) =  j !component index:  j
            p%Jac_u_indx(index,3) =  i !Node:   i
            index = index + 1
         end do !j      
      end do !i
   end do !k
   
   !Module/Mesh/Field: u%InflowOnTower(:,:) = 22;
   do i=1,size(u%InflowOnTower,2) ! numNodes
      do j=1,3
         p%Jac_u_indx(index,1) =  22
         p%Jac_u_indx(index,2) =  j !component index:  j
         p%Jac_u_indx(index,3) =  i !Node:   i
         index = index + 1
      end do !j      
   end do !i
   
   
      !......................................
      ! default perturbations, p%du:
      !......................................
   call allocAry( p%du, 22, 'p%du', ErrStat2, ErrMsg2) ! 22 = number of unique values in p%Jac_u_indx(:,1)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   perturb = 2*D2R
   
   do k=1,p%NumBlades
      perturb_b(k) = 0.2_ReKi*D2R * InputFileData%BladeProps(k)%BlSpn( InputFileData%BladeProps(k)%NumBlNds )
   end do

   if ( u%TowerMotion%NNodes > 0) then
      perturb_t = 0.2_ReKi*D2R * u%TowerMotion%Position( 3, u%TowerMotion%NNodes )
   else
      perturb_t = 0.0_ReKi
   end if   
   
   p%du(1) = perturb_t                    ! u%TowerMotion%TranslationDisp  = 1
   p%du(2) = perturb                      ! u%TowerMotion%Orientation      = 2
   p%du(3) = perturb_t                    ! u%TowerMotion%TranslationVel   = 3
   p%du(4) = perturb_b(1)                 ! u%HubMotion%TranslationDisp    = 4
   p%du(5) = perturb                      ! u%HubMotion%Orientation        = 5
   p%du(6) = perturb                      ! u%HubMotion%RotationVel        = 6
   do i_meshField = 7,9   
      p%du(i_meshField) = perturb         ! u%BladeRootMotion(k)%Orientation = 6+k, for k in [1, 3]
   end do
   do k=1,p%NumBlades         
      p%du(10 + (k-1)*3) = perturb_b(k)   ! u%BladeMotion(k)%TranslationDisp = 10 + (k-1)*3
      p%du(11 + (k-1)*3) = perturb        ! u%BladeMotion(k)%Orientation     = 11 + (k-1)*3
      p%du(12 + (k-1)*3) = perturb_b(k)   ! u%BladeMotion(k)%TranslationVel  = 12 + (k-1)*3
   end do
   do k=1,p%NumBlades         
      p%du(18 + k) = perturb_b(k)         ! u%InflowOnBlade(:,:,k) = 18 + k
   end do      
   p%du(22) = perturb_t                   ! u%InflowOnTower(:,:) = 22
  
         
      !.....................
      ! get names of linearized inputs
      !.....................
   call AllocAry(InitOut%LinNames_u, nu, 'LinNames_u', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(InitOut%RotFrame_u, nu, 'RotFrame_u', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(InitOut%IsLoad_u, nu, 'IsLoad_u', ErrStat2, ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return

   InitOut%IsLoad_u   = .false. ! None of AeroDyn's inputs are loads
   InitOut%RotFrame_u = .false.
      
   index = 1
   FieldMask = .false.
   FieldMask(MASKID_TRANSLATIONDISP) = .true.
   FieldMask(MASKID_Orientation) = .true.
   FieldMask(MASKID_TRANSLATIONVel) = .true.
   call PackMotionMesh_Names(u%TowerMotion, 'Tower', InitOut%LinNames_u, index, FieldMask=FieldMask)
   
   FieldMask(MASKID_TRANSLATIONVel) = .false.
   FieldMask(MASKID_RotationVel) = .true.
   call PackMotionMesh_Names(u%HubMotion, 'Hub', InitOut%LinNames_u, index, FieldMask=FieldMask)

   index_last = index   
   FieldMask = .false.
   FieldMask(MASKID_Orientation) = .true.
   do k = 1,p%NumBlades
      call PackMotionMesh_Names(u%BladeRootMotion(k), 'Blade root '//trim(num2lstr(k)), InitOut%LinNames_u, index, FieldMask=FieldMask)
   end do
   
   FieldMask(MASKID_TRANSLATIONDISP) = .true.
   FieldMask(MASKID_TRANSLATIONVel)  = .true.
   do k=1,p%NumBlades     
      call PackMotionMesh_Names(u%BladeMotion(k), 'Blade '//trim(num2lstr(k)), InitOut%LinNames_u, index, FieldMask=FieldMask)
   end do
   
   do k=1,p%NumBlades
      do i=1,p%NumBlNds
         do j=1,3
            InitOut%LinNames_u(index) = UVW(j)//'-component inflow on blade '//trim(num2lstr(k))//', node '//trim(num2lstr(i))//', m/s'
            index = index + 1
         end do            
      end do
   end do
   !InitOut%RotFrame_u(index_last:index-1) = .true. ! values on the mesh (and from IfW) are in global coordinates, thus not in the rotating frame

   do i=1,p%NumTwrNds
      do j=1,3
         InitOut%LinNames_u(index) = UVW(j)//'-component inflow on tower node '//trim(num2lstr(i))//', m/s'
         index = index + 1
      end do            
   end do
               
   

      !.....................
      ! get names of linearized constraint states (though i don't think we really need them)
      !.....................
   call AllocAry(InitOut%LinNames_z, p%NumBlades*p%NumBlNds, 'LinNames_z', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(InitOut%RotFrame_z, p%NumBlades*p%NumBlNds, 'RotFrame_z', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) return
   InitOut%RotFrame_z = .true.
      
   index = 1      
   do k=1,p%NumBlades ! size(z%BEMT%Phi,2)
      do i=1,p%NumBlNds ! size(z%BEMT%Phi,1)
         InitOut%LinNames_z(index) = 'phi at blade '//trim(num2lstr(k))//', node '//trim(num2lstr(i))//', rad'
         index = index + 1
      end do            
   end do
         
END SUBROUTINE Init_Jacobian
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the u array (and mesh/field it corresponds to)
!! Do not change this without making sure subroutine aerodyn::init_jacobian is consistant with this routine!
SUBROUTINE Perturb_u( p, n, perturb_sign, u, du )

   TYPE(AD_ParameterType)              , INTENT(IN   ) :: p                      !< parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: n                      !< number of array element to use 
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign           !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(AD_InputType)                  , INTENT(INOUT) :: u                      !< perturbed ED inputs
   REAL( R8Ki )                        , INTENT(  OUT) :: du                     !< amount that specific input was perturbed
   

   ! local variables
   integer(intKi)                                      :: ErrStat2
   character(ErrMsgLen)                                :: ErrMsg2
   
   INTEGER                                             :: fieldIndx
   INTEGER                                             :: node
   REAL(R8Ki)                                          :: orientation(3,3)
   REAL(R8Ki)                                          :: angles(3)
      
   fieldIndx = p%Jac_u_indx(n,2) 
   node      = p%Jac_u_indx(n,3) 
   
   du = p%du(  p%Jac_u_indx(n,1) )
   
      ! determine which mesh we're trying to perturb and perturb the input:
   SELECT CASE( p%Jac_u_indx(n,1) )
      
   CASE ( 1) !Module/Mesh/Field: u%TowerMotion%TranslationDisp = 1;      
      u%TowerMotion%TranslationDisp( fieldIndx,node) = u%TowerMotion%TranslationDisp( fieldIndx,node) + du * perturb_sign       
   CASE ( 2) !Module/Mesh/Field: u%TowerMotion%Orientation = 2;
      angles = 0.0_R8Ki
      angles(fieldIndx) = du * perturb_sign
      call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), orientation, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )            
      u%TowerMotion%Orientation(:,:,node) = matmul(u%TowerMotion%Orientation(:,:,node), orientation)      
   CASE ( 3) !Module/Mesh/Field: u%TowerMotion%TranslationVel = 3;
      u%TowerMotion%TranslationVel( fieldIndx,node) = u%TowerMotion%TranslationVel( fieldIndx,node) + du * perturb_sign       
      
   CASE ( 4) !Module/Mesh/Field: u%HubMotion%TranslationDisp = 4;
      u%HubMotion%TranslationDisp(fieldIndx,node) = u%HubMotion%TranslationDisp(fieldIndx,node) + du * perturb_sign            
   CASE ( 5) !Module/Mesh/Field: u%HubMotion%Orientation = 5;
      angles = 0.0_R8Ki
      angles(fieldIndx) = du * perturb_sign
      call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), orientation, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )            
      u%HubMotion%Orientation(:,:,node) = matmul(u%HubMotion%Orientation(:,:,node), orientation)
   CASE ( 6) !Module/Mesh/Field: u%HubMotion%RotationVel = 6;
      u%HubMotion%RotationVel(fieldIndx,node) = u%HubMotion%RotationVel(fieldIndx,node) + du * perturb_sign
   
   CASE ( 7) !Module/Mesh/Field: u%BladeRootMotion(1)%Orientation = 7;
      angles = 0.0_R8Ki
      angles(fieldIndx) = du * perturb_sign
      call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), orientation, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )            
      u%BladeRootMotion(1)%Orientation(:,:,node) = matmul(u%BladeRootMotion(1)%Orientation(:,:,node), orientation)
   CASE ( 8) !Module/Mesh/Field: u%BladeRootMotion(2)%Orientation = 8;
      angles = 0.0_R8Ki
      angles(fieldIndx) = du * perturb_sign
      call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), orientation, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )            
      u%BladeRootMotion(2)%Orientation(:,:,node) = matmul(u%BladeRootMotion(2)%Orientation(:,:,node), orientation)
   CASE ( 9) !Module/Mesh/Field: u%BladeRootMotion(3)%Orientation = 9;
      angles = 0.0_R8Ki
      angles(fieldIndx) = du * perturb_sign
      call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), orientation, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )            
      u%BladeRootMotion(3)%Orientation(:,:,node) = matmul(u%BladeRootMotion(3)%Orientation(:,:,node), orientation)
      
   CASE (10) !Module/Mesh/Field: u%BladeMotion(1)%TranslationDisp = 10;
      u%BladeMotion(1)%TranslationDisp(fieldIndx,node) = u%BladeMotion(1)%TranslationDisp(fieldIndx,node) + du * perturb_sign   
   CASE (11) !Module/Mesh/Field: u%BladeMotion(1)%Orientation = 11;
      angles = 0.0_R8Ki
      angles(fieldIndx) = du * perturb_sign
      call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), orientation, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )            
      u%BladeMotion(1)%Orientation(:,:,node) = matmul(u%BladeMotion(1)%Orientation(:,:,node), orientation)
   CASE (12) !Module/Mesh/Field: u%BladeMotion(1)%TranslationVel  = 12;
      u%BladeMotion(1)%TranslationVel(fieldIndx,node) = u%BladeMotion(1)%TranslationVel(fieldIndx,node) + du * perturb_sign            
      
   CASE (13) !Module/Mesh/Field: u%BladeMotion(2)%TranslationDisp = 13;
      u%BladeMotion(2)%TranslationDisp( fieldIndx,node) = u%BladeMotion(2)%TranslationDisp( fieldIndx,node) + du * perturb_sign       
   CASE (14) !Module/Mesh/Field: u%BladeMotion(2)%Orientation = 14;
      angles = 0.0_R8Ki
      angles(fieldIndx) = du * perturb_sign
      call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), orientation, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )            
      u%BladeMotion(2)%Orientation(:,:,node) = matmul(u%BladeMotion(2)%Orientation(:,:,node), orientation)
   CASE (15) !Module/Mesh/Field: u%BladeMotion(2)%TranslationVel = 15;
      u%BladeMotion(2)%TranslationVel(fieldIndx,node) = u%BladeMotion(2)%TranslationVel(fieldIndx,node) + du * perturb_sign
      
   CASE (16) !Module/Mesh/Field: u%BladeMotion(3)%TranslationDisp = 16;
      u%BladeMotion(3)%TranslationDisp( fieldIndx,node) = u%BladeMotion(3)%TranslationDisp( fieldIndx,node) + du * perturb_sign       
   CASE (17) !Module/Mesh/Field: u%BladeMotion(3)%Orientation     = 17;
      angles = 0.0_R8Ki
      angles(fieldIndx) = du * perturb_sign
      call SmllRotTrans( 'linearization perturbation', angles(1), angles(2), angles(3), orientation, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )            
      u%BladeMotion(3)%Orientation(:,:,node) = matmul(u%BladeMotion(3)%Orientation(:,:,node), orientation)
   CASE (18) !Module/Mesh/Field: u%BladeMotion(3)%TranslationVel  = 18;
      u%BladeMotion(3)%TranslationVel(fieldIndx,node) = u%BladeMotion(3)%TranslationVel(fieldIndx,node) + du * perturb_sign

   CASE (19) !Module/Mesh/Field: u%InflowOnBlade(:,:,1) = 19;
      u%InflowOnBlade(fieldIndx,node,1) = u%InflowOnBlade(fieldIndx,node,1) + du * perturb_sign
   CASE (20) !Module/Mesh/Field: u%InflowOnBlade(:,:,2) = 20;
      u%InflowOnBlade(fieldIndx,node,2) = u%InflowOnBlade(fieldIndx,node,2) + du * perturb_sign
   CASE (21) !Module/Mesh/Field: u%InflowOnBlade(:,:,3) = 21;
      u%InflowOnBlade(fieldIndx,node,3) = u%InflowOnBlade(fieldIndx,node,3) + du * perturb_sign
      
   CASE (22) !Module/Mesh/Field: u%InflowOnTower(:,:)   = 22;
      u%InflowOnTower(fieldIndx,node) = u%InflowOnTower(fieldIndx,node) + du * perturb_sign
      
   END SELECT
      
END SUBROUTINE Perturb_u
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine aerodyn::init_jacobian is consistant with this routine!
SUBROUTINE Compute_dY(p, y_p, y_m, delta_p, delta_m, dY)
   
   TYPE(AD_ParameterType)            , INTENT(IN   ) :: p         !< parameters
   TYPE(AD_OutputType)               , INTENT(IN   ) :: y_p       !< AD outputs at \f$ u + \Delta_p u \f$ or \f$ z + \Delta_p z \f$ (p=plus)
   TYPE(AD_OutputType)               , INTENT(IN   ) :: y_m       !< AD outputs at \f$ u - \Delta_m u \f$ or \f$ z - \Delta_m z \f$ (m=minus)   
   REAL(R8Ki)                        , INTENT(IN   ) :: delta_p   !< difference in inputs or states \f$ delta_p = \Delta_p u \f$ or \f$ delta_p = \Delta_p z \f$
   REAL(R8Ki)                        , INTENT(IN   ) :: delta_m   !< difference in inputs or states \f$ delta_m = \Delta_m u \f$ or \f$ delta_m = \Delta_m z \f$
   REAL(R8Ki)                        , INTENT(INOUT) :: dY(:)     !< column of dYdu or dYdz: \f$ \frac{\partial Y}{\partial u_i} = \frac{y_p - y_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial Y}{\partial z_i} = \frac{y_p - y_m}{2 \, \Delta z}\f$
   
      ! local variables:
   INTEGER(IntKi)    :: k              ! loop over blades
   INTEGER(IntKi)    :: indx_first     ! index indicating next value of dY to be filled 

   
   
   indx_first = 1               
   call PackLoadMesh_dY(y_p%TowerLoad, y_m%TowerLoad, dY, indx_first)
   do k=1,p%NumBlades
      call PackLoadMesh_dY(y_p%BladeLoad(k), y_m%BladeLoad(k), dY, indx_first)                  
   end do
   
   
   do k=1,p%NumOuts
      dY(k+indx_first-1) = y_p%WriteOutput(k) - y_m%WriteOutput(k)
   end do   
   
   
   dY = dY / (delta_p + delta_m)
   
END SUBROUTINE Compute_dY
!----------------------------------------------------------------------------------------------------------------------------------
FUNCTION CheckBEMTInputPerturbations( p, m ) RESULT(ValidPerturb)

   type(AD_ParameterType),  intent(in   )  :: p               !< AD parameters
   type(AD_MiscVarType),    intent(inout)  :: m               !< Misc/optimization variables
   logical                                 :: ValidPerturb    !< if .true., the perturbation is valid; if false, invalid (and thus don't use it) 
   
   integer                                 :: j,k
   
   integer, parameter                      :: indx    = 1  ! index of perturbed input
   integer, parameter                      :: indx_op = 2  ! index of operating point
   real(ReKi)                              :: Vx_prod, Vy_prod
   
   ValidPerturb = .true.
   
   if ( p%BEMT%UseInduction ) then
      if (p%FrozenWake ) then
         
         do k=1,p%NumBlades      
            do j=1,p%NumBlNds         
            
                  ! don't allow the input perturbations to change Vx or Vy so that Vx+AxInd_op=0 and Vy+TnInd_op=0 to 
                  ! avoid ill-conditioning in CalcConstrStateResidual:
               if ( VelocityIsZero( m%BEMT_u(indx)%Vx(j,k)+m%BEMT%AxInd_op(j,k) ) .and. &
                    VelocityIsZero( m%BEMT_u(indx)%Vy(j,k)+m%BEMT%TnInd_op(j,k) ) ) then
                  ValidPerturb = .false.
                  return
               end if

                  ! don't allow the input perturbations to change Vx or Vy so that Vx=0 or Vy=0 to 
                  ! avoid division-by-zero errors in CalcOutput:
               if ( VelocityIsZero( m%BEMT_u(indx)%Vx(j,k) ) .or. VelocityIsZero( m%BEMT_u(indx)%Vy(j,k) ) ) then
                  ValidPerturb = .false.
                  return
               end if
                                        
            end do !j=nodes
         end do !k=blades         
                                    
      else ! not FrozenWake
         
         do k=1,p%NumBlades      
            do j=1,p%NumBlNds         
            
                  ! don't allow the input perturbations to change Vx or Vy far enough to switch sign (or go to zero)
                  ! so as to change solution regions. 
               Vx_prod = m%BEMT_u(indx)%Vx(j,k) * m%BEMT_u(indx_op)%Vx(j,k)
               if ( Vx_prod <= 0.0_ReKi ) then
                  ValidPerturb = .false.
                  return
               elseif ( VelocityIsZero( m%BEMT_u(indx)%Vx(j,k) ) .or. VelocityIsZero( m%BEMT_u(indx_op)%Vx(j,k) )  ) then
                  ValidPerturb = .false.
                  return
               else
                  Vy_prod = m%BEMT_u(indx)%Vy(j,k) * m%BEMT_u(indx_op)%Vy(j,k)
                  if (Vy_prod <= 0.0_ReKi ) then
                     ValidPerturb = .false.
                     return
                  elseif ( VelocityIsZero(  m%BEMT_u(indx)%Vy(j,k) ) .or. VelocityIsZero( m%BEMT_u(indx_op)%Vy(j,k) )) then
                     ValidPerturb = .false.
                     return
                  end if
               end if
                     
            end do !j=nodes
         end do !k=blades
                              
      end if
      
   else ! not UseInduction
      
         do k=1,p%NumBlades      
            do j=1,p%NumBlNds         
            
                  ! don't allow the input perturbations to change Vx or Vy so that Vx=0 or Vy=0:
               if ( EqualRealNos( m%BEMT_u(indx)%Vx(j,k), 0.0_ReKi ) .and. EqualRealNos( m%BEMT_u(indx)%Vy(j,k), 0.0_ReKi ) ) then
                  ValidPerturb = .false.
                  return
               end if
                                        
            end do !j=nodes
         end do !k=blades         
      
   end if
                        
END FUNCTION CheckBEMTInputPerturbations
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE AeroDyn
