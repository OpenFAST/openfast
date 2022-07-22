!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014  DNV KEMA Renewables, Inc.
!
!    This file is part of the IceFloe suite of subroutines
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
!************************************************************************
! modified 8-Jan-2016 by B. Jonkman, NREL to conform to changes in FAST Modularization framework (added MiscVars)

!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************

!  Test calling program for IceFloe

program main

!  rename some things to avoid conflicts
   use IceFloe_Types, dum1=>errID_none, dum2=>aborterrlev, dum3=>errID_warn, dum4=>errID_fatal, dum5=>newline
   use IceFloe
   use IceInputParams
   implicit none

   TYPE(IceFloe_InitInputType)    :: InitInp     ! Input data from FAST for initialization routine
   TYPE(IceFloe_InputType)        :: u(1)        ! An initial guess for the input
   TYPE(IceFloe_ParameterType)    :: p           ! Parameters
   TYPE(IceFloe_OutputType)       :: y           ! Initial system outputs (outputs are not calculated;
   TYPE(IceFloe_InitOutputType)   :: InitOut     ! Output for initialization routine
   INTEGER(IntKi)                :: ErrStat = ErrID_None    ! Error status of the operation
   CHARACTER(1024)               :: ErrMsg = ""             ! Error message if ErrStat /= ErrID_None
   REAL(DbKi)                    :: Interval    ! Coupling interval in seconds: the rate that
   integer(IntKi)                :: outUnitNum

! Below are not currently used in IceFloe - dummy variables and stubbed routines are included however
   TYPE(IceFloe_ContinuousStateType)  :: x           ! Initial continuous states
   TYPE(IceFloe_DiscreteStateType)    :: xd          ! Initial discrete states
   TYPE(IceFloe_ConstraintStateType)  :: z           ! Initial guess of the constraint states
   TYPE(IceFloe_OtherStateType)       :: OtherState  ! Initial other states
   TYPE(IceFloe_MiscVarType)          :: m           ! misc/optimization variables

   character(132) :: outFile
   integer(IntKi) :: n, i, numArg, nSteps, nL
   real(ReKi)     :: pos = 0.0
   real(ReKi)     :: lenTime

   numArg = COMMAND_ARGUMENT_COUNT()
   if (numArg > 0) then
      CALL GET_COMMAND_ARGUMENT( 1, InitInp%InputFile )
   else
      InitInp%InputFile = 'testIce.inp'
   endif
   
   write(*,*) ' Initializing IceFloe...'

! Initialize the NWTC Subroutine Library
   CALL NWTC_Init( )


   InitInp%rootname = InitInp%InputFile
   call IceFloe_Init( InitInp, u(1), p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
   if (ErrStat >= AbortErrLev) then
      call progAbort( ErrMsg )
      stop
   elseif (ErrStat == ErrID_Warn) then
      call progWarn( ErrMsg )
   endif

   call GetRoot( InitInp%InputFile, outFile )
   outFile = trim(outFile)//'.plt'
   call OpenFOutFile ( outUnitNum, outFile )

   pos = 0.0
   
   nSteps = size(p%loadSeries,1)
   if(p%iceType == 5) nSteps = floor(600.0/p%dt)
   write(*,*) ' Now time marching'
   
   nL = p%numLegs
   if (p%singleLoad) nL = 1
   do i = 0, nSteps-1
      call IceFloe_CalcOutput( dble(i)*p%dt, u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
      if (ErrStat >= AbortErrLev) then
         write(*,*) ErrMsg
         stop
      endif
      write(outUnitNum,'(1x,10(1pe12.5,A))') sngl(i)*p%dt, TAB, pos, TAB, u(1)%iceMesh%TranslationVel(1,1),   &
                                             TAB, y%writeoutput(1), TAB, y%writeoutput(2)
!                                             (TAB, y%iceMesh%Force(1,n), TAB, y%iceMesh%Force(2,n), TAB, y%iceMesh%Moment(3,n), n = 1,nL)
      if(p%iceType == 5) call sdof(sngl(p%dt), y%iceMesh%Force(1,1), pos, u(1)%iceMesh%TranslationVel(1,1))
   enddo

   write(*,*) ' Wrap up'

   call IceFloe_End( u(1), p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )

end program main

! single degree of system integrated via 4th order RK
subroutine sdof (dt, F, pos, vel)
   use precision
   implicit none
   real(ReKi), intent(in)  :: dt, F
   real(ReKi), intent(out) :: pos, vel
   real(ReKi), save        :: x   = 0.0
   real(ReKi), save        :: xd  = 0.0
   real(ReKi), save        :: xdd = 0.0
   real(ReKi)              :: k11, k12, k13, k14
   real(ReKi)              :: k21, k22, k23, k24
   real(ReKi), save        :: rampForce = 0.0


! From NREL 5MW monopile in 25m depth (Carbon Trust)
   real(ReKi), parameter   :: mass  = 10.6E6
   real(ReKi), parameter   :: stiff = 30.0E6  
! Tom's example
!   real(ReKi), parameter   :: mass  = 166000.0
!   real(ReKi), parameter   :: stiff = 1.5E9   
   real(ReKi), parameter   :: damp = 0.01
   real(ReKi), save        :: wn

   wn = sqrt(stiff/mass)
   
!  Ramp up the force
   rampForce = (dt/0.1)*F + (1.0-dt/0.1)*rampForce

!  first integrate xdd (acceleration)
   pos = x
   vel = xd

   k11 = rampForce/mass - wn*wn*x - 2.0*damp*wn*xd
   k21 = xd
   x  = pos + 0.5*dt*k21
   xd = vel + 0.5*dt*k11

   k12 = rampForce/mass - wn*wn*x - 2.0*damp*wn*xd
   k22 = xd
   x  = pos + 0.5*dt*k22
   xd = vel + 0.5*dt*k12

   k13 = rampForce/mass - wn*wn*x - 2.0*damp*wn*xd
   k23 = xd
   x  = pos + 0.5*dt*k23
   xd = vel + 0.5*dt*k13

   k14 = rampForce/mass - wn*wn*x - 2.0*damp*wn*xd
   k24 = xd

   pos = pos + dt*(k21 + 2.0*k22 + 2.0*k23 + k24)/6.0
   vel = vel + dt*(k11 + 2.0*k12 + 2.0*k13 + k14)/6.0
   x  = pos
   xd = vel

end subroutine sdof

!**********************************************************************************************************************************
