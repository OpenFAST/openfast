!**********************************************************************************************************************************
! UnsteadyAero_DriverCode: This code tests a stand-alone version of the UnsteadyAero module
!..................................................................................................................................
! LICENSING
! Copyright (C) 2015  National Renewable Energy Laboratory
!
!    This file is part of UnsteadyAero.
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


   
   
   
program UnsteadyAero_Driver

   use NWTC_Library
   use AirfoilInfo
   use AirfoilInfo_Types
   use UnsteadyAero_Types
   use UnsteadyAero
   use UA_Dvr_Subs
   USE VersionInfo

   implicit none

   
   
   
   
    ! Variables
   integer(IntKi), parameter                     :: NumInp = 1           ! Number of inputs sent to HydroDyn_UpdateStates
   
   real(DbKi)  :: dt, t
   integer     :: i, j, k, n 
   real(DbKi)                                    :: InputTime(NumInp)    ! Variable for storing time associated with inputs, in seconds
   type(UA_InitInputType)                        :: InitInData           ! Input data for initialization
   type(UA_InitOutputType)                       :: InitOutData          ! Output data from initialization
   type(UA_DiscreteStateType)                    :: xd                   ! Discrete states
   type(UA_OtherStateType)                       :: OtherState           ! Other/optimization states
   type(UA_MiscVarType)                          :: m                    ! Misc/optimization variables
   type(UA_ParameterType)                        :: p                    ! Parameters
   type(UA_InputType)                            :: u(NumInp)            ! System inputs
   type(UA_OutputType)                           :: y                    ! System outputs
   integer(IntKi)                                :: ErrStat, errStat2    ! Status of error message
   character(ErrMsgLen)                          :: ErrMsg, errMsg2     ! Error message if ErrStat /= ErrID_None
   
   integer, parameter                            :: NumAFfiles = 1
   character(1024)                               :: afNames(NumAFfiles)
   type(AFI_ParameterType)                       :: AFI_Params
   integer, allocatable                          :: AFIndx(:,:)
   character(1024)                               :: outFileName
   integer                                       :: unOutFile
   character(200)                                :: TimeFrmt, Frmt 
   CHARACTER(1024)                               :: dvrFilename          ! Filename and path for the driver input file.  This is passed in as a command line argument when running the Driver exe.
   TYPE(UA_Dvr_InitInput)                        :: dvrInitInp           ! Initialization data for the driver program
   integer                                       :: nSteps
   real(DbKi)                                    :: simTime  
   integer                                       :: nSimSteps
   character(1024)                               :: RoutineName
   real(DbKi), allocatable                       :: timeArr(:)
   real(ReKi), allocatable                       :: AOAarr(:)
   real(ReKi), allocatable                       :: Uarr(:) !RRD
   CHARACTER(200)                                :: git_commit
   TYPE(ProgDesc), PARAMETER   :: version   = ProgDesc( 'UnsteadyAero Driver', '', '' )  ! The version number of this program.
      ! Initialize the NWTC library
   call NWTC_Init()
   
      ! Initialize error handling variables
   ErrMsg  = ''
   ErrStat = ErrID_None
   
   RoutineName = 'UnsteadyAero_Driver'
    
   
      ! Display the copyright notice
   CALL DispCopyrightLicense( version )   
      ! Obtain OpenFAST git commit hash
   git_commit = QueryGitVersion()
      ! Tell our users what they're running
   CALL WrScr( ' Running '//GetNVD( version )//' a part of OpenFAST - '//TRIM(git_Commit)//NewLine//' linked with '//TRIM( GetNVD( NWTC_Ver ))//NewLine )
   
   
   
   
      ! Parse the driver file if one was provided, if not, then set driver parameters using hardcoded values
   if ( command_argument_count() > 1 ) then
      call print_help()
      stop
   end if
  
   
      ! Establish initialization inputs which are fixed for the stand-alone driver, but would be
      ! variable for a coupled simulation
   InitInData%nNodesPerBlade  = 1 
   InitInData%numBlades       = 1
   
      ! Set up initialization data
   allocate(AFIndx(InitInData%nNodesPerBlade,InitInData%numBlades), STAT = ErrStat)
      if ( ErrStat /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error trying to allocate InitInData%AFIndx.', ErrStat, ErrMsg, 'UnsteadyAeroTest')  
         call Cleanup()
         stop       
      end if
   
   allocate(InitInData%c(InitInData%nNodesPerBlade,InitInData%numBlades), STAT = ErrStat)
      if ( ErrStat /= 0 ) then
         call SetErrStat( ErrID_Fatal, 'Error trying to allocate InitInData%c.', ErrStat, ErrMsg, 'UnsteadyAeroTest')  
         call Cleanup()
         stop       
      end if
   
      
      ! Parse the driver input file and run the simulation based on that file
      
   if ( command_argument_count() == 1 ) then
      
      call get_command_argument(1, dvrFilename)
      call ReadDriverInputFile( dvrFilename, dvrInitInp, errStat2, errMsg2 )
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            stop
         end if
      InitInData%a_s          = dvrInitInp%SpdSound
      InitInData%c(1,1)       = dvrInitInp%Chord
      InitInData%UAMod        = dvrInitInp%UAMod 
      InitInData%Flookup      = dvrInitInp%Flookup
   
   else
      
      dvrInitInp%OutRootName  = '.\TestingUA_Driver'
      InitInData%UAMod        = 1  
      InitInData%Flookup      = .FALSE.
      InitInData%a_s          = 340.29 ! m/s  
      InitInData%c(1,1)       = 1.0
      dvrInitInp%InflowVel    = 30.0 ! m/s
      dvrInitInp%Re           = 75  ! million
      dvrInitInp%UAMod        = 1
      dvrInitInp%Flookup      = .FALSE.
      dvrInitInp%AirFoil1     = './OSU075_FAST.txt'
      dvrInitInp%SimMod       = 1
      dvrInitInp%NCycles      = 3.0
      dvrInitInp%Frequency    = 1.2 ! Hz
      dvrInitInp%StepsPerCycle= 180
      dvrInitInp%Amplitude    = 10.0 ! deg
      dvrInitInp%Mean         =  2.0 ! deg
      dvrInitInp%Phase        =  0 ! steps of a cycle
      dvrInitInp%InputsFile   = ''
      
   end if
   
   if ( dvrInitInp%SimMod == 1 ) then
         ! Using the frequency and NCycles, determine how long the simulation needs to run
      simTime   = dvrInitInp%NCycles/dvrInitInp%Frequency
      nSimSteps = dvrInitInp%StepsPerCycle*dvrInitInp%NCycles
      dt        = simTime / nSimSteps
      
   else
         ! Read time-series data file with a 1 line header and then each row contains time-step data with 4, white-space-separated columns
         ! time  Angle-fo-attack  
      call ReadTimeSeriesData( dvrInitInp%InputsFile, nSimSteps, timeArr, AOAarr, Uarr, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if ( ErrStat >= AbortErrLev ) then
            call Cleanup()
            stop
         end if
      dt = (timeArr(nSimSteps) - timeArr(1)) / nSimSteps
   end if
      
     
      ! All nodes/blades are using the same 2D airfoil
   afNames(1)  = dvrInitInp%AirFoil1
   AFIndx(1,1) = 1
   
      ! Initialize the Airfoil Info Params
   call Init_AFI( NumAFfiles, afNames, InitInData%Flookup, dvrInitInp%UseCm, AFI_Params, errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
      if ( ErrStat >= AbortErrLev ) then
         call Cleanup()
         stop
      end if
   
    ! Initialize UnsteadyAero
   call UA_Init( InitInData, u(1), p, xd, OtherState, y, m, dt, InitOutData, errStat2, errMsg2 ) 
      call SetErrStat( errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         stop
      end if
   
   if (p%NumOuts > 0) then
         ! Initialize the output file
         ! Open the file for output
      outFileName = trim(dvrInitInp%OutRootName)//'.out'
      call GetNewUnit( unOutFile )
   
      call OpenFOutFile ( unOutFile, outFileName, errStat2, errMsg2 ) 
         call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            stop
         end if
      
      
         ! Write the output file header
      p%OutSFmt = 'A19'
      p%OutFmt  = 'ES19.5e2'
      p%Delim   =''

      Frmt = '('//trim(Int2LStr(p%NumOuts*p%numBlades*p%nNodesPerBlade))//'(:,A,'//trim( p%OutSFmt )//'))'
      
      write (unOutFile,'(/,A/)', IOSTAT=ErrStat)  'These predictions were generated by UnSteadyAero on '//CurDate()//' at '//CurTime()//'.'
      write (unOutFile,'(/,A/)', IOSTAT=ErrStat)  'Driver file name: '//trim(dvrFilename)
      
         ! Write the names of the output parameters:
      write(unOutFile, '(A15)', ADVANCE='no')  trim( 'Time' )
      write(unOutFile, Frmt, ADVANCE='no')   ( p%Delim, trim( InitOutData%WriteOutputHdr(I)   ), i=1,p%NumOuts*p%numBlades*p%nNodesPerBlade )   
      write (unOutFile,'()', IOSTAT=ErrStat)          ! write the line return
      
         ! Write the units of the output parameters:  
      write(unOutFile, '(A15)', ADVANCE='no')  trim( '(sec)' ) 
      write(unOutFile, Frmt, ADVANCE='no')   ( p%Delim, trim( InitOutData%WriteOutputUnt(I)   ), i=1,p%NumOuts*p%numBlades*p%nNodesPerBlade )
      write (unOutFile,'()', IOSTAT=ErrStat)          ! write the line return
   
      TimeFrmt = '(F15.4)'
      Frmt     = '('//trim(Int2LStr(p%NumOuts*p%numBlades*p%nNodesPerBlade))//'(:,A,'//trim( p%OutFmt )//'))'
   end if
   
      ! Set inputs which do not vary with node or time
   u(1)%U  = dvrInitInp%InflowVel  ! m/s
   u(1)%Re = dvrInitInp%Re  ! not used at the moment
            
      ! time marching loop
   do n = 1, nSimSteps
      if ( dvrInitInp%SimMod == 1 ) then
         t            = (n-1)*dt
         InputTime(1) = t
         u(1)%alpha =   (dvrInitInp%Amplitude * sin((n+dvrInitInp%Phase-1)*2*pi/dvrInitInp%StepsPerCycle) + dvrInitInp%Mean)*pi/180.0   ! This needs to be in radians
      
      else
         ! Load timestep data from the time-series inputs which were previous read from input file
         InputTime(1) = timeArr(n)
         u(1)%alpha   = AOAarr (n)*pi/180.0   ! This needs to be in radians
         u(1)%U = Uarr(n) !RRD 
         t            = timeArr(n)
      end if
      
         ! set the inputs for the node
      
      
      do j = 1,InitInData%numBlades
         do i = 1,InitInData%nNodesPerBlade
     
               ! Need to use MiscVar to store which element we are operating on
            m%iBladeNode = i
            m%iBlade     = j
            
               ! Use existing states to compute the outputs
            call UA_CalcOutput(u(1),  p, xd, OtherState, AFI_Params%AFInfo(AFIndx(i,j)), y, m, errStat2, errMsg2 )
               call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
               if (ErrStat >= AbortErrLev) then
                  call Cleanup()
                  stop
               end if
            
 
               ! Prepare states for next time step
            call UA_UpdateStates(i,j,u(1), p, xd, OtherState, AFI_Params%AFInfo(AFIndx(i,j)), m, errStat2, errMsg2 )
               call SetErrStat(errStat2, errMsg2, ErrStat, ErrMsg, RoutineName )
               if (ErrStat >= AbortErrLev) then
                  call Cleanup()
                  stop
               end if
               
         end do
      end do
      
         ! Generate file outputs
      if (p%NumOuts > 0) then
         write (unOutFile,TimeFrmt,ADVANCE='no')  t  
         write (unOutFile,Frmt,ADVANCE='no')   ( p%Delim,  y%WriteOutput(k)  , k=1,p%NumOuts*p%numBlades*p%nNodesPerBlade )   
         write (unOutFile,'()', IOSTAT=ErrStat)          ! write the line return
      end if      
   end do
   
  ! write (unOutFile,'(/,A/)', IOSTAT=ErrStat)  'This output file was closed on '//CurDate()//' at '//CurTime()//'.'
   
   !-------------------------------------------------------------------------------------------------
   ! Close our output file
   !-------------------------------------------------------------------------------------------------
   

   call Cleanup()
   
   
   contains
   
   !====================================================================================================
   subroutine Cleanup()
   !     The routine cleans up the module echo file and resets the NWTC_Library, reattaching it to 
   !     any existing echo information
   !----------------------------------------------------------------------------------------------------  
      
      if ( ErrStat /= ErrID_None ) print *, ErrMsg
      
      if (p%NumOuts > 0) close( unOutFile, IOSTAT = ErrStat )
      
      
   end subroutine Cleanup

   
   
   
   subroutine print_help()
    print '(a)', 'usage: '
    print '(a)', ''
    print '(a)', 'UnsteadyAero_Driver.exe [driverfilename]'
    print '(a)', ''
    print '(a)', 'Where the optional argument, driverfilename, is the name of the UnsteadyAero driver input file.'
    print '(a)', ''

   end subroutine print_help
   
end program UnsteadyAero_Driver

