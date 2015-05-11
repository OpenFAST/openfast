module AeroDyn
    
   use NWTC_Library
   use BEMT
   use BEMT_Types
   use AeroDyn_Types
   use AirfoilInfo_Types
   use AirfoilInfo
  

   implicit none

   private
   
   type(ProgDesc), parameter  :: AD_Ver = ProgDesc( 'AeroDyn', 'v14.00.00a-gjh', '01-Oct-2014' )
   character(*),   parameter  :: AD_Nickname = 'AD'
   integer(IntKi), parameter  :: AD_numChanPerNode = 15  ! TODO This needs to be set dynamically ?? 9/18/14 GJH
   REAL(ReKi), PARAMETER     ::epsilon = 1e-6
   
   
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

end module AeroDyn