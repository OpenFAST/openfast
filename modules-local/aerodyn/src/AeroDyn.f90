module AeroDyn
    
   use NWTC_Library
   use BEMT
   use BEMT_Types
   use AeroDyn_Types
   use AirfoilInfo_Types
   use AirfoilInfo
  

   implicit none

   private
   
   type(ProgDesc), parameter  :: AD_Ver = ProgDesc( 'AeroDyn', 'v15.00.00a-gjh', '12-May-2015' )
   character(*),   parameter  :: AD_Nickname = 'AD'
   integer(IntKi), parameter  :: AD_numChanPerNode = 15  ! TODO This needs to be set dynamically ?? 9/18/14 GJH
   
   REAL(ReKi), PARAMETER     ::epsilon = 1e-6
   
   INTEGER(IntKi), PARAMETER        :: MaxBl    =  3                                   ! Maximum number of blades allowed in simulation
   INTEGER(IntKi), PARAMETER        :: MaxOutPts = 1  !bjj: fix me!!!!!
   
   ! ..... Public Subroutines ...................................................................................................

   public :: AD_Init                           ! Initialization routine
   public :: AD_End                            ! Ending routine (includes clean up)

   public :: AD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   public :: AD_CalcOutput                     ! Routine for computing outputs

   public :: AD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   
  
   contains
  
!subroutine Set_BEMT_InitInp(AD_InitInp, BEMT_InitInp, errStat, errMsg)
!      type(AD_InitInputType) , intent(in   )   :: AD_InitInp           ! Input data for initialization
!      type(BEMT_InitInputType), intent(  out)   :: BEMT_InitInp           ! Input data for initialization
!      integer(IntKi)         , intent(inout)   :: errStat              ! Status of error message
!      character(*)           , intent(inout)   :: errMsg               ! Error message if ErrStat /= ErrID_None
!
!         ! locals
!      integer(IntKi)                           :: i, j
!      integer(IntKi)                           :: errStat2             ! local status of error message
!      character(len(errMsg))                   :: errMsg2              ! local error message if ErrStat /= ErrID_None
!      real(ReKi)                               :: rHub, zTip, deltar
!      errStat2 = ErrID_None
!      errMsg2  = ''
!   
!      AD_InitInp%numBladeNodes  = WTP_Data%numSeg
!      AD_InitInp%numBlades      = WTP_Data%numBlade
!   
!      allocate ( AD_InitInp%chord(AD_InitInp%numBladeNodes, AD_InitInp%numBlades), STAT = errStat2 )
!      if ( errStat2 /= 0 ) then
!         errStat2 = ErrID_Fatal
!         errMsg2  = 'Error allocating memory for AD_InitInp%chord.'
!         call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'Set_WTP_InitInp' )
!         return
!      end if 
!   
!      allocate ( AD_InitInp%AFindx(AD_InitInp%numBladeNodes), STAT = errStat2 )
!      if ( errStat2 /= 0 ) then
!         errStat2 = ErrID_Fatal
!         errMsg2  = 'Error allocating memory for AD_InitInp%chord array.'
!         call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'Set_WTP_InitInp' )
!         return
!      end if 
!      
!      allocate ( AD_InitInp%zLocal(AD_InitInp%numBladeNodes, AD_InitInp%numBlades), STAT = errStat2 )
!      if ( errStat2 /= 0 ) then
!         errStat2 = ErrID_Fatal
!         errMsg2  = 'Error allocating memory for AD_InitInp%zLocal array.'
!         call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'Set_WTP_InitInp' )
!         return
!      end if 
!   
!      allocate ( AD_InitInp%zTip(AD_InitInp%numBlades), STAT = errStat2 )
!      if ( errStat2 /= 0 ) then
!         errStat2 = ErrID_Fatal
!         errMsg2  = 'Error allocating memory for AD_InitInp%zTip array.'
!         call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'Set_WTP_InitInp' )
!         return
!      end if 
!      
!      allocate ( AD_InitInp%zHub(AD_InitInp%numBlades), STAT = errStat2 )
!      if ( errStat2 /= 0 ) then
!         errStat2 = ErrID_Fatal
!         errMsg2  = 'Error allocating memory for AD_InitInp%rHub array.'
!         call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'Set_WTP_InitInp' )
!         return
!      end if 
!   
!      do i=1,AD_InitInp%numBladeNodes
!         
!         AD_InitInp%AFindx(i) = WTP_Data%BladeData(i)%AFfile  
!      end do
!   
!      do j=1,AD_InitInp%numBlades
!         BEMT_InitInp%zTip(j)  = AD_InitInp%zTip(j)
!         BEMT_InitInp%zHub(j)  = AD_InitInp%zHub(j)
!         
!         do i=1,AD_InitInp%numBladeNodes
!            BEMT_InitInp%chord (i,j)  = AD_InitInp%chord (i,j) 
!            BEMT_InitInp%zLocal(i,j)  = AD_InitInp%zLocal(i,j)
!         end do
!      end do
!      
!      AD_InitInp%DT                  = BEMT_InitInp%DT                    
!      AD_InitInp%airDens             = BEMT_InitInp%airDens          
!      AD_InitInp%kinVisc             = BEMT_InitInp%kinVisc          
!      AD_InitInp%skewWakeMod         = BEMT_InitInp%skewWakeMod      
!      AD_InitInp%useTipLoss          = BEMT_InitInp%useTipLoss       
!      AD_InitInp%useHubLoss          = BEMT_InitInp%useHubLoss       
!      AD_InitInp%useTanInd           = BEMT_InitInp%useTanInd        
!      AD_InitInp%useAIDrag           = BEMT_InitInp%useAIDrag        
!      AD_InitInp%useTIDrag           = BEMT_InitInp%useTIDrag        
!      AD_InitInp%numReIterations     = BEMT_InitInp%numReIterations  
!      AD_InitInp%maxIndIterations    = BEMT_InitInp%maxIndIterations 
!
!
!end subroutine Set_BEMT_InitInp
  
subroutine AD_SetInitOut(p, InitOut, errStat, errMsg)

   type(AD_InitOutputType),       intent(  out)  :: InitOut     ! output data
   type(AD_ParameterType),        intent(in   )  :: p           ! Parameters
   integer(IntKi),                intent(inout)  :: errStat     ! Error status of the operation
   character(*),                  intent(inout)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(len(errMsg))                        :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
   
   integer(IntKi)                                :: i, j
   character(len=10) :: chanPrefix
      ! Initialize variables for this routine

   errStat2 = ErrID_None
   errMsg2  = ""
   
   call AllocAry( InitOut%WriteOutputHdr, p%numOuts, 'WriteOutputHdr', errStat2, errMsg2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for InitOut%WriteOutputHdr.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_SetInitOut' )
      return
   end if 
   
   call AllocAry( InitOut%WriteOutputUnt, p%numOuts, 'WriteOutputUnt', errStat2, errMsg2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for InitOut%WriteOutputUnt.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_SetInitOut' )
      return
   end if
   
      ! Loop over blades and nodes to populate the output channel names and units
   
   do j=1,p%numBlades
      do i=1,p%numBladeNodes
         

         chanPrefix = "B"//trim(num2lstr(j))//"N"//trim(num2lstr(i))
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 1 ) = trim(chanPrefix)//"Theta"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 1 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 2 ) = trim(chanPrefix)//"Psi"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 2 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 3 ) = trim(chanPrefix)//"Vx"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 3 ) = '  (m/s)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 4 ) = trim(chanPrefix)//"Vy"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 4 ) = '  (m/s)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 5 ) = ' '//trim(chanPrefix)//"AxInd"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 5 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 6 ) = ' '//trim(chanPrefix)//"TanInd"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 6 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 7 ) = trim(chanPrefix)//"IndVel"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 7 ) = '  (m/s)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 8 ) = ' '//trim(chanPrefix)//"Phi"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 8 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 9 ) = ' '//trim(chanPrefix)//"AOA"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 9 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 10 ) = ' '//trim(chanPrefix)//"Cl"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 10 ) = '   (-)   '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 11 ) = ' '//trim(chanPrefix)//"Cd"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 11 ) = '   (-)   '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 12 ) = ' '//trim(chanPrefix)//"Cx"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 12 ) = '   (-)   '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 13 ) = ' '//trim(chanPrefix)//"Cy"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 13 ) = '   (-)   '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 14 ) = ' '//trim(chanPrefix)//"Fx"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 14 ) = '  (N/m)  '
         InitOut%WriteOutputHdr( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 15 ) = ' '//trim(chanPrefix)//"Fy"
         InitOut%WriteOutputUnt( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 15 ) = '  (N/m)  '
      end do
   end do
   
   
   
   end subroutine AD_SetInitOut
   
!----------------------------------------------------------------------------------------------------------------------------------   
subroutine AD_SetParameters( InitInp, p, errStat, errMsg )
! This routine is called from AD_Init.
! The parameters are set here and not changed during the simulation.
!..................................................................................................................................
   type(AD_InitInputType),       intent(inout)  :: InitInp     ! Input data for initialization routine, out is needed because of copy below
   type(AD_ParameterType),       intent(inout)  :: p           ! Parameters
   integer(IntKi),                intent(inout)  :: errStat     ! Error status of the operation
   character(*),                  intent(inout)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(len(errMsg))                        :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
   integer(IntKi)                                :: i, j
   INTEGER(IntKi)                                :: i1
   INTEGER(IntKi)                                :: i1_l  ! lower bounds for an array dimension
   INTEGER(IntKi)                                :: i1_u  ! upper bounds for an array dimension
   
      ! Initialize variables for this routine

   errStat2 = ErrID_None
   errMsg2  = ""

   p%numBladeNodes  = InitInp%numBladeNodes 
   p%numBlades      = InitInp%numBlades    
   p%numOuts        = p%numBladeNodes*p%numBlades*AD_numChanPerNode
   ! p%BEMT_SkewWakeMod = InitInp%BEMT_SkewWakeMod
   IF (ALLOCATED(InitInp%AFInfo)) THEN
      i1_l = LBOUND(InitInp%AFInfo,1)
      i1_u = UBOUND(InitInp%AFInfo,1)
      IF (.NOT. ALLOCATED(p%AFInfo)) THEN 
         ALLOCATE(p%AFInfo(i1_l:i1_u),STAT=ErrStat)
         IF (ErrStat /= 0) THEN 
            ErrStat = ErrID_Fatal 
            ErrMsg = 'AFI_CopyParam: Error allocating p%AFInfo.'
            RETURN
         END IF
      END IF
      DO i1 = i1_l, i1_u
         CALL AFI_Copyafinfotype( InitInp%AFInfo(i1), p%AFInfo(i1), MESH_NEWCOPY, ErrStat, ErrMsg )
      ENDDO
   ENDIF
   !call AFI_Copyafinfotype( InitInp%AFInfo, p%AFInfo, MESH_NEWCOPY, ErrStat, ErrMsg )
   
   allocate ( p%chord(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for p%chord.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_SetParameters' )
      return
   end if 
   
   allocate ( p%AFindx(p%numBladeNodes), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for p%chord array.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_SetParameters' )
      return
   end if 
   
   allocate ( p%tipLossConst(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for p%tipLossConst array.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_SetParameters' )
      return
   end if 
   
   allocate ( p%hubLossConst(p%numBladeNodes, p%numBlades), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for p%hubLossConst array.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_SetParameters' )
      return
   end if 
   
   do i=1,p%numBladeNodes
         
      p%AFindx(i) = InitInp%AFindx(i) 
   end do
   
      ! Compute the tip and hub loss constants using the distances along the blade (provided as input for now) TODO: See how this really needs to be handled. GJH 
   do j=1,p%numBlades
      do i=1,p%numBladeNodes
         p%chord(i,j)        = InitInp%chord(i,j)
         p%tipLossConst(i,j) = p%numBlades*(InitInp%zTip    (j) - InitInp%zLocal(i,j)) / (2.0*InitInp%zLocal(i,j))
         p%hubLossConst(i,j) = p%numBlades*(InitInp%zLocal(i,j) - InitInp%zHub    (j)) / (2.0*InitInp%zHub    (j))
      end do
   end do
   
   
   p%DT               = InitInp%DT                             
   p%airDens          = InitInp%airDens          
   p%kinVisc          = InitInp%kinVisc          
   p%skewWakeMod      = InitInp%skewWakeMod     
   p%useTipLoss       = InitInp%useTipLoss       
   p%useHubLoss       = InitInp%useHubLoss       
   p%useTanInd        = InitInp%useTanInd        
   p%useAIDrag        = InitInp%useAIDrag           
   p%useTIDrag        = InitInp%useTIDrag           
   p%numReIterations  = InitInp%numReIterations  
   p%maxIndIterations = InitInp%maxIndIterations 
   !p%NumOuts          = 2 + p%numBlades*p%numBladeNodes*AD_numChanPerNode  ! TODO This needs to be computed some other way 9/18/14 GJH   
   
end subroutine AD_SetParameters


subroutine AD_AllocOutput(y, p, errStat, errMsg)
   type(AD_OutputType),           intent(  out)  :: y           ! Input data
   type(AD_ParameterType),       intent(in   )  :: p           ! Parameters
   integer(IntKi),                intent(inout)  :: errStat     ! Error status of the operation
   character(*),                  intent(inout)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(len(errMsg))                        :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation

      ! Initialize variables for this routine

   errStat2 = ErrID_None
   errMsg2  = ""
   
   
   call AllocAry( y%WriteOutput, p%numOuts, 'WriteOutput', errStat2, errMsg2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for y%WriteOutput.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_AllocOutput' )
      return
   end if 
   y%WriteOutput = 0.0_ReKi
   
end subroutine AD_AllocOutput


!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_AllocInput( u, p, errStat, errMsg )
! This routine is called from AD_Init.
!  
!  
!..................................................................................................................................

   type(AD_InputType),           intent(  out)  :: u           ! Input data
   type(AD_ParameterType),       intent(in   )  :: p           ! Parameters
   integer(IntKi),                intent(inout)  :: errStat     ! Error status of the operation
   character(*),                  intent(inout)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(len(errMsg))                        :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation

      ! Initialize variables for this routine

   errStat2 = ErrID_None
   errMsg2  = ""

   allocate ( u%theta( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for u%theta.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_AllocInput' )
      return
   end if 
   u%theta = 0.0_ReKi
   
   allocate ( u%psi( p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for u%psi.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_AllocInput' )
      return
   end if 
   u%psi = 0.0_ReKi
   
   allocate ( u%Vx( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for u%Vx.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_AllocInput' )
      return
   end if 
   u%Vx = 0.0_ReKi
   
   allocate ( u%Vy( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for u%Vy.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_AllocInput' )
      return
   end if 
   u%Vy = 0.0_ReKi
  
   allocate ( u%Vinf( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for u%Vinf.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_AllocInput' )
      return
   end if 
   u%Vinf = 0.0_ReKi
      
   allocate ( u%rTip( p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for u%rTip.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_AllocInput' )
      return
   end if 
   u%rTip = 0.0_ReKi
   
   allocate ( u%rLocal( p%numBladeNodes, p%numBlades ), STAT = errStat2 )
   if ( errStat2 /= 0 ) then
      errStat2 = ErrID_Fatal
      errMsg2  = 'Error allocating memory for u%rLocal.'
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'AD_AllocInput' )
      return
   end if 
   u%rLocal = 0.0_ReKi
   
   
   u%lambda = 0.0_ReKi
   u%omega  = 0.0_ReKi
   
end subroutine AD_AllocInput

!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

   type(AD_InitInputType),       intent(inout)  :: InitInp     ! Input data for initialization routine, needs to be inout because there is a copy of some data in InitInp in AD_SetParameters()
   type(AD_InputType),           intent(  out)  :: u           ! An initial guess for the input; input mesh must be defined
   type(AD_ParameterType),       intent(  out)  :: p           ! Parameters
   type(AD_ContinuousStateType), intent(  out)  :: x           ! Initial continuous states
   type(AD_DiscreteStateType),   intent(  out)  :: xd          ! Initial discrete states
   type(AD_ConstraintStateType), intent(  out)  :: z           ! Initial guess of the constraint states
   type(AD_OtherStateType),      intent(  out)  :: OtherState  ! Initial other/optimization states
   type(AD_OutputType),          intent(  out)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                !   only the output mesh is initialized)
   real(DbKi),                    intent(inout)  :: interval    ! Coupling interval in seconds: the rate that
                                                                !   (1) AD_UpdateStates() is called in loose coupling &
                                                                !   (2) AD_UpdateDiscState() is called in tight coupling.
                                                                !   Input is the suggested time from the glue code;
                                                                !   Output is the actual coupling interval that will be used
                                                                !   by the glue code.
   type(AD_InitOutputType),      intent(  out)  :: InitOut     ! Output for initialization routine
   integer(IntKi),                intent(  out)  :: errStat     ! Error status of the operation
   character(*),                  intent(  out)  :: errMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables
   character(len(errMsg))                        :: errMsg2                 ! temporary Error message if ErrStat /= ErrID_None
   integer(IntKi)                                :: errStat2                ! temporary Error status of the operation
   type(AFI_InitInputType)                       :: AFI_InitInputs
   !type(AFI_ParameterType)                       :: AFI_Params
   integer                                      :: i, unEC, i1, i1_l, i1_u
   character(1024)                              :: rootName, echoFile
   !type(BEMT_InitInputType)                       :: BEMT_InitInData
   type(BEMT_InitOutputType)                      :: BEMT_InitOutData
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""


      ! Initialize the NWTC Subroutine Library

   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   call DispNVD( AD_Ver )
   
      ! Initialize the Airfoil Info module
      ! Setup Airfoil info
   AFI_InitInputs%NumAFfiles = InitInp%NumAF
   
   allocate ( AFI_InitInputs%FileNames( AFI_InitInputs%NumAFfiles ), STAT=ErrStat )
   if ( ErrStat /= 0 )  then
      ErrStat = ErrID_Fatal 
      ErrMsg = 'AD_Init: Error allocating AFI_InitInputs%FileNames.'
      call AD_InitCleanup()
   endif
   
   do i=1,AFI_InitInputs%NumAFfiles
      AFI_InitInputs%FileNames(i) = InitInp%AF_File(i)
   end do
   
   AFI_InitInputs%UA_Model    = 0
   AFI_InitInputs%NumCoefs    = 3
   AFI_InitInputs%InCol_Alfa  = 1
   AFI_InitInputs%InCol_Cl    = 2
   AFI_InitInputs%InCol_Cd    = 3
   AFI_InitInputs%InCol_Cm    = 4
   AFI_InitInputs%InCol_Cpmin = 0

     ! Open the echo file.

   CALL GetNewUnit ( UnEc, ErrStat, ErrMsg )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort ( NewLine//TRIM( ADJUSTL( ErrMsg ) ), .FALSE., 10.0, ErrStat )
   ENDIF ! ( ErrStatLcl /= 0 )

   CALL GetRoot ( AFI_InitInputs%FileNames(1), RootName )
   EchoFile = TRIM( RootName )//'.ech' 

   CALL OpenEcho ( UnEc, EchoFile, ErrStat, ErrMsg, AD_Ver )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort ( NewLine//TRIM( ADJUSTL( ErrMsg ) ), .FALSE., 10.0, ErrStat )
   ENDIF ! ( ErrStatLcl /= 0 )
   
      ! Call AFI_Init to read in and process the airfoil files.
      ! This includes creating the spline coefficients to be used for interpolation.

   CALL AFI_Init ( AFI_InitInputs, p%AFI_Params, ErrStat, ErrMsg, UnEc )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort ( NewLine//TRIM( ADJUSTL( ErrMsg ) ), .FALSE., 10.0, ErrStat )
   ENDIF ! ( ErrStatLcl /= 0 )
   
   !call Set_BEMT_InitInp(InitInp, BEMT_InitInData, errStat, errMsg)
   
   ! TODO: Look into MOVE_ALLOC()  GJH
   
      !augment the initialization inputs with the data from AFI_Init
   !IF (ALLOCATED(p%AFI_Params%AFInfo)) THEN
   !   i1_l = LBOUND(p%AFI_Params%AFInfo,1)
   !   i1_u = UBOUND(p%AFI_Params%AFInfo,1)
   !   IF (.NOT. ALLOCATED(InitInp%BEMT%AFInfo)) THEN 
   !      ALLOCATE(InitInp%BEMT%AFInfo(i1_l:i1_u),STAT=ErrStat)
   !      IF (ErrStat /= 0) THEN 
   !         ErrStat = ErrID_Fatal 
   !         ErrMsg = 'AD_Init: Error allocating InitInp%BEMT%AFInfo.'
   !         call AD_InitCleanup()
   !      END IF
   !   END IF
   !   DO i1 = LBOUND(p%AFI_Params%AFInfo,1), UBOUND(p%AFI_Params%AFInfo,1)
   !      CALL AFI_Copyafinfotype( p%AFI_Params%AFInfo(i1), InitInp%BEMT%AFInfo(i1), MESH_NEWCOPY, ErrStat, ErrMsg )
   !   ENDDO
   !ELSE
   !   ErrMsg = 'Invalid AirfoilInfo parameters returned from AD_Init()'
   !   ErrStat = ErrID_Fatal
   !   call AD_InitCleanup()
   !ENDIF
   !call AFI_Copyafinfotype( AFI_Params%AFInfo, AD_InitInData%AFInfo, MESH_NEWCOPY, ErrStat, ErrMsg )
   
     !............................................................................................
      ! Define parameters here
      !............................................................................................
   call AD_SetParameters( InitInp, p, errStat, errMsg )
   if (errStat >= AbortErrLev) return
   
   call AD_SetInitOut(p, InitOut, errStat, errMsg)
   
   call AD_AllocInput( u, p, errStat, errMsg ) 
   call AD_AllocOutput(y, p, errStat, errMsg)
   
      ! Initialized the BEM module
   call BEMT_Init(InitInp%BEMT, OtherState%BEMT_u, p%BEMT,  x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, OtherState%BEMT_y, interval, BEMT_InitOutData, errStat, errMsg )
   if (errStat >= AbortErrLev) then
      ! Clean up and exit
      call AD_InitCleanup()
   end if
   
      ! Open an output file and write the header information
   !call AD_InitializeOutputFile(AD_Ver, delim, outFmtS, outFileRoot, unOutFile, errStat, errMsg)
   
   
   if (errStat >= AbortErrLev) then
      ! Clean up and exit
      call AD_InitCleanup()
   end if
   
   ! ---------------------------------------------------
   ! DO NOT destroy this data because the driver will be using some of the initinp quantities each time step
   !
   !   ! Destroy initialization data
   !
   !call BEMT_DestroyInitInput(  BEMT_InitInData,  errStat, errMsg )
   !if (errStat >= AbortErrLev) then
   !   ! Clean up and exit
   !   call AD_InitCleanup()
   !end if
   ! ---------------------------------------------------
   call BEMT_DestroyInitOutput( BEMT_InitOutData, errStat, errMsg )
   if (errStat >= AbortErrLev) then
      ! Clean up and exit
      call AD_InitCleanup()
   end if

contains

subroutine AD_InitCleanup()
end subroutine AD_InitCleanup

end subroutine AD_Init


!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(AD_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(AD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(AD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(AD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(AD_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(AD_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



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


         ! Destroy the output data:

      CALL AD_DestroyOutput( y, ErrStat, ErrMsg )




END SUBROUTINE AD_End
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, errStat, errMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time t; Continuous and discrete states are updated for t + Interval
!..................................................................................................................................

      real(DbKi),                         intent(in   ) :: t          ! Current simulation time in seconds
      integer(IntKi),                     intent(in   ) :: n          ! Current simulation time step n = 0,1,...
      type(AD_InputType),                intent(inout) :: u(:)       ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
      real(DbKi),                         intent(in   ) :: utimes(:)  ! Times associated with u(:), in seconds
      type(AD_ParameterType),            intent(in   ) :: p          ! Parameters
      type(AD_ContinuousStateType),      intent(inout) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      type(AD_DiscreteStateType),        intent(inout) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      type(AD_ConstraintStateType),      intent(inout) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      type(AD_OtherStateType),           intent(inout) :: OtherState ! Other/optimization states
      integer(IntKi),                     intent(  out) :: errStat    ! Error status of the operation
      character(*),                       intent(  out) :: errMsg     ! Error message if ErrStat /= ErrID_None

      integer(IntKi)                                  :: nTime
      
             
     
      nTime = size(u)
      
         ! This needs to extract the inputs from the AD data types (mesh) and massage them for the BEMT module
      call SetInputsForBEMT(u(nTime), OtherState%BEMT_u, errStat, errMsg)  
      
      
         ! Call into the BEMT update states    NOTE:  This is a non-standard framework interface!!!!!  GJH
      call BEMT_UpdateStates(t, n, OtherState%BEMT_u,  p%BEMT, x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, p%AFI_Params%AFInfo, errStat, errMsg)
      
      
end subroutine AD_UpdateStates


!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! Routine for computing outputs, used in both loose and tight coupling.
! This SUBROUTINE is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
! NOTE: the descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
! for a complete description of each output parameter.
! NOTE: no matter how many channels are selected for output, all of the outputs are calcalated
! All of the calculated output channels are placed into the OtherState%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(AD_InputType),           INTENT(IN   )  :: u           ! Inputs at Time t
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(AD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(AD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                               !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   INTEGER(IntKi)                               :: i,j
   character(10)                                :: chanPrefix
   real(ReKi)                                   :: q
   
   call SetInputsForBEMT(u,OtherState%BEMT_u, ErrStat, ErrMsg)
   
!-------------------------------------------------------
!     Start of BEMT-level calculations of outputs
!-------------------------------------------------------
      ! Call the BEMT module CalcOutput.  Notice that the BEMT outputs are purposely attached to AeroDyn's OtherState structure to
      ! avoid issues with the coupling code
      ! Also we are forced to copy the mesh-based input data into the BEMT input data structure (which doesn't use a mesh)
      !
      ! NOTE: the x%BEMT, xd%BEMT, and OtherState%BEMT are simply dummy variables because the BEMT module does not use them or return them
      !
   
   call BEMT_CalcOutput(t, OtherState%BEMT_u, p%BEMT, x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, p%AFI_Params%AFInfo, OtherState%BEMT_y, ErrStat, ErrMsg )
   ! TODO Check error status
   
!-------------------------------------------------------   
!     End of BEMT calculations  
!-------------------------------------------------------
   
!-------------------------------------------------------
!     Start AeroDyn-level calculations of outputs
!-------------------------------------------------------
   
!-------------------------------------------------------   
!     End of AeroDyn calculations  
!-------------------------------------------------------   
 
   
         
         
            ! Loop over blades and nodes to populate the output channel names and units
   
   do j=1,p%numBlades
      do i=1,p%numBladeNodes
         
         
         chanPrefix = "B"//trim(num2lstr(j))//"N"//trim(num2lstr(i))
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 1 ) = OtherState%BEMT_u%theta(i,j)*R2D
      
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 2 ) = OtherState%BEMT_u%psi(j)*R2D
        
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 3 ) = OtherState%BEMT_u%Vx(i,j)
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 4 ) = OtherState%BEMT_u%Vy(i,j)
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 5 ) = OtherState%BEMT_y%AxInduction(i,j)
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 6 ) = OtherState%BEMT_y%TanInduction(i,j)
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 7 ) = OtherState%BEMT_y%inducedVel(i,j)
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 8 ) = OtherState%BEMT_y%phi(i,j)*R2D
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 9 ) = OtherState%BEMT_y%AOA(i,j)*R2D
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 10 ) = OtherState%BEMT_y%Cl(i,j)
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 11 ) = OtherState%BEMT_y%Cd(i,j)
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 12 ) = OtherState%BEMT_y%Cx(i,j)
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 13 ) = OtherState%BEMT_y%Cy(i,j)
         
         q = 0.5*p%AirDens*p%chord(i,j)*OtherState%BEMT_y%inducedVel(i,j)**2
         
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 14 ) =q*OtherState%BEMT_y%Cx(i,j)
         y%WriteOutput( (j-1)*p%numBladeNodes*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 15 ) = -q *OtherState%BEMT_y%Cy(i,j)
         
      end do
   end do
!   
!-------------------------------------------------------
   
   
end subroutine AD_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                    INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD_InputType),           INTENT(IN   )   :: u           ! Inputs at Time
      TYPE(AD_ParameterType),       INTENT(IN   )   :: p           ! Parameters
      TYPE(AD_ContinuousStateType), INTENT(IN   )   :: x           ! Continuous states at Time
      TYPE(AD_DiscreteStateType),   INTENT(IN   )   :: xd          ! Discrete states at Time
      TYPE(AD_ConstraintStateType), INTENT(INOUT)   :: z           ! Constraint states at Time (possibly a guess)
      TYPE(AD_OtherStateType),      INTENT(INOUT)   :: OtherState  ! Other/optimization states
      TYPE(AD_ConstraintStateType), INTENT(  OUT)   :: z_residual  ! Residual of the constraint state equations using
                                                                   !     the input values described above
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      !# set epsilon
   !REAL(ReKi), PARAMETER     ::epsilon = 1e-6
   
      ! Local variables
   INTEGER    :: i,j
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   call SetInputsForBEMT(u,OtherState%BEMT_u, ErrStat, ErrMsg)
   
   call BEMT_CalcConstrStateResidual( Time, OtherState%BEMT_u, p%BEMT, x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, z_residual%BEMT, p%AFI_Params%AFInfo, ErrStat, ErrMsg )
   
END SUBROUTINE AD_CalcConstrStateResidual

subroutine SetInputsForBEMT(u,BEMT_u, errStat, errMsg)
   type(AD_InputType),            intent(in   )  :: u           ! Inputs at Time
   type(BEMT_InputType),           intent(inout)  :: BEMT_u           ! Inputs at Time
   integer(IntKi),                intent(  out)  :: ErrStat     ! Error status of the operation
   character(*),                  intent(  out)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      
      ! In the actually version of AeroDyn, we will need to copy data from a mesh into the BEMT arrays.
   BEMT_u%theta    = u%theta  
   BEMT_u%chi0     = u%gamma   !TODO,  this needs to change to account for tilt GJH 11/17/14
   BEMT_u%psi      = u%psi    
   BEMT_u%omega    = u%omega  
   BEMT_u%Vx       = u%Vx     
   BEMT_u%Vy       = u%Vy  
   BEMT_u%Vinf     = u%Vinf
   BEMT_u%lambda   = u%lambda 
   BEMT_u%rTip     = u%rTip   
   BEMT_u%rLocal   = u%rLocal   
   
end subroutine SetInputsForBEMT

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_ReadInput( InputFileName, InputFileData, Default_DT, OutFileRoot, NumBl, ErrStat, ErrMsg )
! This subroutine reads the input file and stores all the data in the AD_InputFile structure.
! It does not perform data validation.
!..................................................................................................................................

      ! Passed variables
   REAL(DbKi),           INTENT(IN)       :: Default_DT      ! The default DT (from glue code)

   CHARACTER(*), INTENT(IN)               :: InputFileName   ! Name of the input file
   CHARACTER(*), INTENT(IN)               :: OutFileRoot     ! The rootname of all the output files written by this routine.

   TYPE(AD_InputFile),   INTENT(OUT)      :: InputFileData   ! Data stored in the module's input file

   INTEGER(IntKi),       INTENT(IN)       :: NumBl           ! Number of blades for this model
   INTEGER(IntKi),       INTENT(OUT)      :: ErrStat         ! The error status code
   CHARACTER(*),         INTENT(OUT)      :: ErrMsg          ! The error message, if an error occurred

      ! local variables

   INTEGER(IntKi)                         :: I
   INTEGER(IntKi)                         :: UnEcho          ! Unit number for the echo file
   INTEGER(IntKi)                         :: ErrStat2        ! The error status code
   CHARACTER(ErrMsgLen)                   :: ErrMsg2         ! The error message, if an error occurred

   CHARACTER(1024)                        :: ADBlFile(MaxBl) ! File that contains the blade information (specified in the primary input file)
   CHARACTER(*), PARAMETER                :: RoutineName = 'AD_ReadInput'
   
   
      ! initialize values:

   ErrStat = ErrID_None
   ErrMsg  = ''
   InputFileData%DTAero = Default_DT  ! the glue code's suggested DT for the module (may be overwritten in ReadPrimaryFile())

      ! get the primary/platform input-file data
      ! sets UnEcho, ADBlFile
   
   CALL ReadPrimaryFile( InputFileName, InputFileData, ADBlFile, OutFileRoot, UnEcho, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      

      ! get the blade input-file data
      
   ALLOCATE( InputFileData%BladeProps( NumBl ), STAT = ErrStat2 )
   IF (ErrStat2 /= 0) THEN
      CALL SetErrStat(ErrID_Fatal,"Error allocating memory for BladeProps.", ErrStat, ErrMsg, RoutineName)
      CALL Cleanup()
      RETURN
   END IF
      
   DO I=1,NumBl
      CALL ReadBladeInputs ( ADBlFile(I), InputFileData%BladeProps(I), InputFileData%NumBlNds, UnEcho, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName//TRIM(':Blade')//TRIM(Num2LStr(I)))
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   END DO
   

      ! cleanup

   CALL Cleanup ( )


CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up before exiting this subroutine
   !...............................................................................................................................

      IF ( UnEcho > 0 ) CLOSE( UnEcho )


   END SUBROUTINE Cleanup

END SUBROUTINE AD_ReadInput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadPrimaryFile( InputFile, InputFileData, ADBlFile, OutFileRoot, UnEc, ErrStat, ErrMsg )
! This routine reads in the primary AeroDyn input file and places the values it reads in the InputFileData structure.
!   It opens and prints to an echo file if requested.
!..................................................................................................................................


   IMPLICIT                        NONE

      ! Passed variables
   INTEGER(IntKi),     INTENT(OUT)     :: UnEc                                ! I/O unit for echo file. If > 0, file is open for writing.
   INTEGER(IntKi),     INTENT(OUT)     :: ErrStat                             ! Error status

   CHARACTER(*),       INTENT(OUT)     :: ADBlFile(MaxBl)                     ! name of the files containing blade inputs
   CHARACTER(*),       INTENT(IN)      :: InputFile                           ! Name of the file containing the primary input data
   CHARACTER(*),       INTENT(OUT)     :: ErrMsg                              ! Error message
   CHARACTER(*),       INTENT(IN)      :: OutFileRoot                         ! The rootname of the echo file, possibly opened in this routine

   TYPE(AD_InputFile), INTENT(INOUT)   :: InputFileData                       ! All the data in the AeroDyn input file
   
      ! Local variables:
   REAL(ReKi)                    :: TmpRAry(2)                                ! A temporary array to read a table from the input file
   INTEGER(IntKi)                :: I                                         ! loop counter
   INTEGER(IntKi)                :: UnIn                                      ! Unit number for reading file
     
   INTEGER(IntKi)                :: ErrStat2, IOS                             ! Temporary Error status
   LOGICAL                       :: Echo                                      ! Determines if an echo file should be written
   CHARACTER(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)               :: PriPath                                   ! Path name of the primary file
   CHARACTER(1024)               :: FTitle                                    ! "File Title": the 2nd line of the input file, which contains a description of its contents
   CHARACTER(200)                :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
   CHARACTER(*), PARAMETER       :: RoutineName = 'ReadPrimaryFile'
   
      ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   UnEc = -1
   Echo = .FALSE.   
   CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
   

   CALL AllocAry( InputFileData%OutList, MaxOutPts, "Outlist", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        
      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      
                  
      
   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file. 
   ! If Echo is TRUE, rewind and write on the second try.
   
   I = 1 !set the number of times we've read the file
   DO 
   !----------- HEADER -------------------------------------------------------------
   
      CALL ReadCom( UnIn, InputFile, 'File header: Module Version (line 1)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
      CALL ReadStr( UnIn, InputFile, FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   
   
   !----------- GENERAL OPTIONS ----------------------------------------------------
   
      CALL ReadCom( UnIn, InputFile, 'Section Header: General Options', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
         ! Echo - Echo input to "<RootName>.AD.ech".
   
      CALL ReadVar( UnIn, InputFile, Echo, 'Echo',   'Echo flag', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   
      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop
   
         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
      
      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
   
      CALL OpenEcho ( UnEc, TRIM(OutFileRoot)//'.ech', ErrStat2, ErrMsg2, AD_Ver )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   
      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(AD_Ver%Name)//' primary input file "'//TRIM( InputFile )//'":'
         
      REWIND( UnIn, IOSTAT=ErrStat2 )  
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         END IF         
      
   END DO    

   IF (NWTC_VerboseLevel == NWTC_Verbose) THEN
      CALL WrScr( ' Heading of the '//TRIM(aD_Ver%Name)//' input file: ' )      
      CALL WrScr( '   '//TRIM( FTitle ) )
   END IF
   
   
      ! DTAero - Time interval for aerodynamic calculations {or default} (s):
   CALL ReadVar( UnIn, InputFile, InputFileData%DTAero, "DTAero", "Time interval for aerodynamic calculations {or default} (s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "DEFAULT" ) /= 1 ) THEN ! If it's not "default", read this variable; otherwise use the value already stored in InputFileData%DTAero
         READ( Line, *, IOSTAT=IOS) InputFileData%DTAero
            CALL CheckIOS ( IOS, InputFile, 'DTAero', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF   
      
      ! WakeMod - Type of wake/induction model {0=none, 1=BEMT} (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%WakeMod, "WakeMod", "Type of wake/induction model {0=none, 1=BEMT} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! AFAeroBladeMod - Type of blade airfoil aerodynamics model {1=steady model, 2=Beddoes-Leishman unsteady model} (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%AFAeroBladeMod, "AFAeroBladeMod", "Type of blade airfoil aerodynamics model {1=steady model, 2=Beddoes-Leishman unsteady model} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! AFAeroTwrMod - Type of tower airfoil aerodynamics model {1=steady model, 2=Beddoes-Leishman unsteady model} (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%AFAeroTwrMod, "AFAeroTwrMod", "Type of tower airfoil aerodynamics model {1=steady model, 2=Beddoes-Leishman unsteady model} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TwrPotent - Calculate tower influence on wind based on potential flow around the tower? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrPotent, "TwrPotent", "Calculate tower influence on wind based on potential flow around the tower? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TwrShadow - Calculate downstream tower shadow? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrShadow, "TwrShadow", "Calculate downstream tower shadow? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TwrAero - Calculate tower aerodynamic loads? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrAero, "TwrAero", "Calculate tower aerodynamic loads? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
            
   !----------- ENVIRONMENTAL CONDITIONS -------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Environmental Conditions', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! AirDens - Air density (kg/m^3):
   CALL ReadVar( UnIn, InputFile, InputFileData%AirDens, "AirDens", "Air density (kg/m^3)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! KinVisc - Kinematic air viscosity (m^2/s):
   CALL ReadVar( UnIn, InputFile, InputFileData%KinVisc, "KinVisc", "Kinematic air viscosity (m^2/s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! SpdSound - Speed of sound (m/s):
   CALL ReadVar( UnIn, InputFile, InputFileData%SpdSound, "SpdSound", "Speed of sound (m/s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF

   !----------- BLADE-ELEMENT/MOMENTUM THEORY OPTIONS ------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Blade-Element/Momentum Theory Options', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! SkewMod - Type of skewed-wake correction model {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%SkewMod, "SkewMod", "Type of skewed-wake correction model {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TipLoss - Use the Prandtl tip-loss model? [used only when WakeMod=1] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TipLoss, "TipLoss", "Use the Prandtl tip-loss model? [used only when WakeMod=1] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! HubLoss - Use the Prandtl hub-loss model? [used only when WakeMod=1] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubLoss, "HubLoss", "Use the Prandtl hub-loss model? [used only when WakeMod=1] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TanInd - Include tangential induction in BEMT calculations? [used only when WakeMod=1] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TanInd, "TanInd", "Include tangential induction in BEMT calculations? [used only when WakeMod=1] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! AIDrag - Include the drag term in the axial-induction calculation? [used only when WakeMod=1] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%AIDrag, "AIDrag", "Include the drag term in the axial-induction calculation? [used only when WakeMod=1] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TIDrag - Include the drag term in the tangential-induction calculation? [used only when WakeMod=1 and TanInd=TRUE] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TIDrag, "TIDrag", "Include the drag term in the tangential-induction calculation? [used only when WakeMod=1 and TanInd=TRUE] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! IndToler - Convergence tolerance for BEM induction factors [used only when WakeMod=1] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%IndToler, "IndToler", "Convergence tolerance for BEM induction factors [used only when WakeMod=1] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! MaxIter - Maximum number of iteration steps [used only when WakeMod=1] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%MaxIter, "MaxIter", "Maximum number of iteration steps [used only when WakeMod=1] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
         
   !----------- BEDDOES-LEISHMAN UNSTEADY AIRFOIL AERODYNAMICS OPTIONS -------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Beddoes-Leishman Unsteady Airfoil Aerodynamics Options', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! DSMod - Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez’s variant (changes in Cn,Cc,Cm), 3=Minemma/Pierce variant (changes in Cc and Cm)} [used only when AFAreoMod=2] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%DSMod, "DSMod", "Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez’s variant (changes in Cn,Cc,Cm), 3=Minemma/Pierce variant (changes in Cc and Cm)} [used only when AFAreoMod=2] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! FLookup - Flag to indicate whether a lookup for f’ will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files [used only when AFAreoMod=2] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%FLookup, "FLookup", "Flag to indicate whether a lookup for f’ will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files [used only when AFAreoMod=2] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )     
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
                  
   !----------- AIRFOIL INFORMATION ------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Airfoil Information', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! InCol_Alfa - The column in the airfoil tables that contains the angle of attack (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Alfa, "InCol_Alfa", "The column in the airfoil tables that contains the angle of attack (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! InCol_Cl - The column in the airfoil tables that contains the lift coefficient (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Cl, "InCol_Cl", "The column in the airfoil tables that contains the lift coefficient (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! InCol_Cd - The column in the airfoil tables that contains the drag coefficient (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Cd, "InCol_Cd", "The column in the airfoil tables that contains the drag coefficient (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! InCol_Cm - The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Cm, "InCol_Cm", "The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! InCol_Cpmin - The column in the airfoil tables that contains the drag coefficient; use zero if there is no Cpmin column (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Cpmin, "InCol_Cpmin", "The column in the airfoil tables that contains the drag coefficient; use zero if there is no Cpmin column (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NumAFfiles - Number of airfoil files used (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NumAFfiles, "NumAFfiles", "Number of airfoil files used (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

         ! Allocate space to hold AFNames
      ALLOCATE( InputFileData%AFNames(InputFileData%NumAFfiles), STAT=ErrStat2)
         IF (ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, "Error allocating AFNames.", ErrStat, ErrMsg, RoutineName)
            CALL Cleanup()
            RETURN
         END IF
               
      ! AFNames - Airfoil file names (NumAFfiles lines) (quoted strings):
   DO I = 1,InputFileData%NumAFfiles            
      CALL ReadVar ( UnIn, InputFile, InputFileData%AFNames(I), 'AFNames('//TRIM(Num2Lstr(I))//')', 'Airfoil '//TRIM(Num2Lstr(I))//' file name', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( PathIsRelative( InputFileData%AFNames(I) ) ) InputFileData%AFNames(I) = TRIM(PriPath)//TRIM(InputFileData%AFNames(I))
   END DO      
             
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF

   !----------- ROTOR/BLADE PROPERTIES  --------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Rotor/Blade Properties', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! UseBlCm - Include aerodynamic pitching moment in calculations? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%UseBlCm, "UseBlCm", "Include aerodynamic pitching moment in calculations? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN
            
      ! NumBlNds - Number of blade nodes used in the analysis (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NumBlNds, "NumBlNds", "Number of blade nodes used in the analysis (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      ! ADBlFile - Names of files containing distributed aerodynamic properties for each blade (see AD_BladeInputFile type):
   DO I = 1,MaxBl            
      CALL ReadVar ( UnIn, InputFile, ADBlFile(I), 'ADBlFile('//TRIM(Num2Lstr(I))//')', 'Name of file containing distributed aerodynamic properties for blade '//TRIM(Num2Lstr(I)), ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( PathIsRelative( ADBlFile(I) ) ) ADBlFile(I) = TRIM(PriPath)//TRIM(ADBlFile(I))
   END DO      
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF

   !----------- TOWER INFLUENCE AND AERODYNAMICS  ----------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Tower Influence and Aerodynamics', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! TwrWakeCnst - Tower wake constant {0.0 - full potential flow, 0.1 - Bak model} [used only when TwrPotent=True] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrWakeCnst, "TwrWakeCnst", "Tower wake constant {0.0 - full potential flow, 0.1 - Bak model} [used only when TwrPotent=True] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TwrUseCm - Include aerodynamic pitching moment in tower aerodynamic load calculations? [used only when TwrAero=True] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrUseCm, "TwrUseCm", "Include aerodynamic pitching moment in tower aerodynamic load calculations? [used only when TwrAero=True] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NumTwrNds - Number of tower nodes used in the analysis (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NumTwrNds, "NumTwrNds", "Number of tower nodes used in the analysis (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      
   !....... tower properties ...................
   CALL ReadCom( UnIn, InputFile, 'Section Header: Tower Property Channels', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ReadCom( UnIn, InputFile, 'Section Header: Tower Property Units', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! allocate space for tower inputs:
   CALL AllocAry( InputFileData%TwrElev,  InputFileData%NumTwrNds, 'TwrElev',  ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( InputFileData%TwrTwist, InputFileData%NumTwrNds, 'TwrTwist', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( InputFileData%TwrChord, InputFileData%NumTwrNds, 'TwrChord', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( InputFileData%TwrAFID,  InputFileData%NumTwrNds, 'TwrAFID',  ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error if we didn't allocate space for the next inputs
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
            
   DO I=1,InputFileData%NumTwrNds
      READ( UnIn, *, IOStat=IOS ) InputFileData%TwrElev(I), InputFileData%TwrTwist(I), InputFileData%TwrChord(I), InputFileData%TwrAFID(I)
         CALL CheckIOS( IOS, InputFile, 'Tower properties row '//TRIM(Num2LStr(I)), NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
         IF (UnEc > 0 .AND. IOS == 0) THEN
            WRITE( UnEc, "(3(F9.4,1x),I9)", IOStat=IOS) InputFileData%TwrElev(I), InputFileData%TwrTwist(I), InputFileData%TwrChord(I), InputFileData%TwrAFID(I)
         END IF         
   END DO
   InputFileData%TwrTwist = InputFileData%TwrTwist*D2R
               
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
                  
   !----------- OUTPUTS  -----------------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Outputs', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! SumPrint - Generate a summary file listing input options and interpolated properties to <rootname>.AD.sum? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%SumPrint, "SumPrint", "Generate a summary file listing input options and interpolated properties to <rootname>.AD.sum? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! NBlOuts - Number of blade node outputs [0 - 9] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NBlOuts, "NBlOuts", "Number of blade node outputs [0 - 9] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      IF ( InputFileData%NBlOuts > SIZE(InputFileData%BlOutNd) ) THEN
         CALL SetErrStat( ErrID_Warn, ' Warning: number of blade output nodes exceeds '//&
                           TRIM(Num2LStr(SIZE(InputFileData%BlOutNd))) //'.', ErrStat, ErrMsg, RoutineName )
         InputFileData%NBlOuts = SIZE(InputFileData%BlOutNd)
      END IF
      
      ! BlOutNd - Blade nodes whose values will be output (-):
   CALL ReadAry( UnIn, InputFile, InputFileData%BlOutNd, InputFileData%NBlOuts, "BlOutNd", "Blade nodes whose values will be output (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      ! NTwOuts - Number of tower node outputs [0 - 9] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NTwOuts, "NTwOuts", "Number of tower node outputs [0 - 9] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      IF ( InputFileData%NTwOuts > SIZE(InputFileData%TwOutNd) ) THEN
         CALL SetErrStat( ErrID_Warn, ' Warning: number of tower output nodes exceeds '//&
                           TRIM(Num2LStr(SIZE(InputFileData%TwOutNd))) //'.', ErrStat, ErrMsg, RoutineName )
         InputFileData%NTwOuts = SIZE(InputFileData%TwOutNd)
      END IF
      
      ! TwOutNd - Tower nodes whose values will be output (-):
   CALL ReadAry( UnIn, InputFile, InputFileData%TwOutNd, InputFileData%NTwOuts, "TwOutNd", "Tower nodes whose values will be output (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
                  
   !----------- OUTLIST  -----------------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: OutList', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! OutList - List of user-requested output channels (-):
   CALL ReadOutputList ( UnIn, InputFile, InputFileData%OutList, InputFileData%NumOuts, 'OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc  )     ! Routine in NWTC Subroutine Library
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !---------------------- END OF FILE -----------------------------------------
      
   CALL Cleanup( )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up any local variables and closes input files
   !...............................................................................................................................

   IF (UnIn > 0) CLOSE ( UnIn )

   END SUBROUTINE Cleanup
   !...............................................................................................................................
END SUBROUTINE ReadPrimaryFile      
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadBladeInputs ( ADBlFile, BladeKInputFileData, NumBlNds, UnEc, ErrStat, ErrMsg )
! This routine reads a blade input file.
!..................................................................................................................................


      ! Passed variables:

   TYPE(AD_BladePropsType),  INTENT(INOUT)  :: BladeKInputFileData                 ! Data for Blade K stored in the module's input file
   CHARACTER(*),             INTENT(IN)     :: ADBlFile                            ! Name of the blade input file data
   INTEGER(IntKi),           INTENT(IN)     :: NumBlNds                            ! Number of blade nodes to read from this file
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                ! I/O unit for echo file. If present and > 0, write to UnEc

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message


      ! Local variables:

   REAL(ReKi)                   :: TmpRAry(17)                                     ! Temporary variable to read table from file (up to 17 columns)

   INTEGER(IntKi)               :: I                                               ! A generic DO index.
   INTEGER( IntKi )             :: UnIn                                            ! Unit number for reading file
   INTEGER(IntKi)               :: ErrStat2 , IOS                                  ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                         ! Temporary Err msg
   CHARACTER(*), PARAMETER      :: RoutineName = 'ReadBladeInputs'

   ErrStat = ErrID_None
   ErrMsg  = ""
   UnIn = -1
      
   ! Allocate space for these variables
   
   
   
   
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


      ! Open the input file for blade K.

   CALL OpenFInpFile ( UnIn, ADBlFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- HEADER -------------------------------------------------------

      ! Skip the header.

   CALL ReadCom ( UnIn, ADBlFile, 'unused blade file header line 1', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom ( UnIn, ADBlFile, 'unused blade file header line 2', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   !  -------------- Blade properties table ------------------------------------------                                    
   CALL ReadCom ( UnIn, ADBlFile, 'Section header: Blade Properties', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom ( UnIn, ADBlFile, 'Table header: names', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom ( UnIn, ADBlFile, 'Table header: units', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   IF ( ErrStat>= AbortErrLev ) THEN 
      CALL Cleanup()
      RETURN
   END IF
   
      
      ! allocate space for blade inputs:
   CALL AllocAry( BladeKInputFileData%BlSpn,   NumBlNds, 'BlSpn',   ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlCrvAC, NumBlNds, 'BlCrvAC', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlSwpAC, NumBlNds, 'BlSwpAC', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlCrvAng,NumBlNds, 'BlCrvAng',ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlTwist, NumBlNds, 'BlTwist', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlChord, NumBlNds, 'BlChord', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlAFID,  NumBlNds, 'BlAFID',  ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error if we didn't allocate space for the next inputs
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
            
   DO I=1,NumBlNds
      READ( UnIn, *, IOStat=IOS ) BladeKInputFileData%BlSpn(I), BladeKInputFileData%BlCrvAC(I), BladeKInputFileData%BlSwpAC(I), &
                                  BladeKInputFileData%BlCrvAng(I), BladeKInputFileData%BlTwist(I), BladeKInputFileData%BlChord(I), &
                                  BladeKInputFileData%BlAFID(I)  
         CALL CheckIOS( IOS, ADBlFile, 'Blade properties row '//TRIM(Num2LStr(I)), NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               ! Return on error if we couldn't read this line
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL Cleanup()
               RETURN
            END IF
         
         IF (UnEc > 0) THEN
            WRITE( UnEc, "(6(F9.4,1x),I9)", IOStat=IOS) BladeKInputFileData%BlSpn(I), BladeKInputFileData%BlCrvAC(I), BladeKInputFileData%BlSwpAC(I), &
                                  BladeKInputFileData%BlCrvAng(I), BladeKInputFileData%BlTwist(I), BladeKInputFileData%BlChord(I), &
                                  BladeKInputFileData%BlAFID(I)
         END IF         
   END DO
   BladeKInputFileData%BlCrvAng = BladeKInputFileData%BlCrvAng*D2R
   BladeKInputFileData%BlTwist  = BladeKInputFileData%BlTwist*D2R
                  
   !  -------------- END OF FILE --------------------------------------------

   CALL Cleanup()
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up local variables and closes files
   !...............................................................................................................................

      IF (UnIn > 0) RETURN

   END SUBROUTINE Cleanup

END SUBROUTINE ReadBladeInputs      
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE AeroDyn