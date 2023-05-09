!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
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
!
! References:
!  [1] Brooks, T. F.; Pope, D. S. & Marcolini, M. A., Airfoil self-noise and prediction, 
!      NASA, NASA, 1989. https://ntrs.nasa.gov/search.jsp?R=19890016302
module AeroAcoustics
    
   use NWTC_Library
   use AeroAcoustics_Types
   use AeroAcoustics_IO
   use NWTC_LAPACK
   USE NWTC_FFTPACK
   implicit none

   private
   ! ..... Public Subroutines ...................................................................................................
   public :: AA_Init                           ! Initialization routine
   public :: AA_End                            ! Ending routine (includes clean up)
   public :: AA_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   public :: AA_CalcOutput                     ! Routine for computing outputs

   contains    
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine AA_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
   type(AA_InitInputType),       intent(in   ) :: InitInp       !< Input data for initialization routine
   type(AA_InputType),           intent(  out) :: u             !< An initial guess for the input; input mesh must be defined
   type(AA_ParameterType),       intent(  out) :: p             !< Parameters
   type(AA_ContinuousStateType), intent(  out) :: x             !< Initial continuous states
   type(AA_DiscreteStateType),   intent(  out) :: xd            !< Initial discrete states
   type(AA_ConstraintStateType), intent(  out) :: z             !< Initial guess of the constraint states
   type(AA_OtherStateType),      intent(  out) :: OtherState    !< Initial other states
   type(AA_OutputType),          intent(  out) :: y             !< Initial system outputs (outputs are not calculated;
                                                                !!   only the output mesh is initialized)
   type(AA_MiscVarType),         intent(  out) :: m             !< Initial misc/optimization variables
   real(DbKi),                   intent(inout) :: interval      !< Coupling interval in seconds: the rate that
                                                                !!   (1) AA_UpdateStates() is called in loose coupling &
                                                                !!   (2) AA_UpdateDiscState() is called in tight coupling.
                                                                !!   Input is the suggested time from the glue code;
                                                                !!   Output is the actual coupling interval that will be used
                                                                !!   by the glue code.
   type(AA_InitOutputType),      intent(  out) :: InitOut       !< Output for initialization routine
   integer(IntKi),               intent(  out) :: errStat       !< Error status of the operation
   character(*),                 intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
   type(AA_InputFile)                          :: InputFileData ! Data stored in the module's input file
   integer(IntKi)                              :: UnEcho        ! Unit number for the echo file
   character(*), parameter                     :: RoutineName = 'AA_Init'
   
   ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""
   UnEcho  = -1
   ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )
   ! Display the module information
   call DispNVD( AA_Ver )

   ! To get rid of a compiler warning.
   x%DummyContState           = 0.0_SiKi
   z%DummyConstrState         = 0.0_SiKi
   OtherState%DummyOtherState = 0.0_SiKi

   !bjj: note that we haven't validated p%NumBlades before using it below!
   p%NumBlades = InitInp%NumBlades ! need this before reading the AD input file so that we know how many blade files to read
   p%RootName  = TRIM(InitInp%RootName)//'.NN'
   
   ! Read the primary AeroAcoustics input file in AeroAcoustics_IO
   call ReadInputFiles( InitInp%InputFile, InitInp%AFInfo, InputFileData, interval, p%RootName, UnEcho, ErrStat2, ErrMsg2 )
   if (Failed()) return
      
   ! Validate the inputs
   call ValidateInputData(InputFileData, p%NumBlades, ErrStat2, ErrMsg2); if (Failed()) return
    
   ! Validate Initialization Input data ( not found in the AeroAcoustics input file )
   if (InitInp%AirDens <= 0.0)  call SetErrStat ( ErrID_Fatal, 'The air density (AirDens) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InitInp%KinVisc <= 0.0)  call SetErrStat ( ErrID_Fatal, 'The kinesmatic viscosity (KinVisc) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InitInp%SpdSound <= 0.0) call SetErrStat ( ErrID_Fatal, 'The speed of sound (SpdSound) must be greater than zero.', ErrStat, ErrMsg, RoutineName )  
   if (Failed()) return

   ! Define parameters
   call SetParameters( InitInp, InputFileData, p, ErrStat2, ErrMsg2 ); if(Failed()) return
   ! Define and initialize inputs 
   call Init_u( u, p, errStat2, errMsg2 ); if(Failed()) return

   ! Define outputs here
   call Init_y(y, u, p, errStat2, errMsg2); if(Failed()) return

   ! Initialize states and misc vars
   call Init_MiscVars(m, p, u, y, errStat2, errMsg2); if(Failed()) return
   call Init_States(xd, p,  errStat2, errMsg2); if(Failed()) return

   ! Define initialization output here
   call AA_SetInitOut(p, InitOut, errStat2, errMsg2); if(Failed()) return
   call AA_InitializeOutputFile(p, InputFileData,InitOut,errStat2, errMsg2); if(Failed()) return
   call Cleanup() 
      
contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call Cleanup()
    end function Failed

   subroutine Cleanup()
      CALL AA_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2 )
      IF ( UnEcho > 0 ) CLOSE( UnEcho )
   end subroutine Cleanup
end subroutine AA_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets AeroAcoustics parameters for use during the simulation; these variables are not changed after AA_Init.
subroutine SetParameters( InitInp, InputFileData, p, ErrStat, ErrMsg )
    TYPE(AA_InitInputType),       INTENT(IN   ) :: InitInp        !< Input data for initialization routine, out is needed because of copy below
    TYPE(AA_InputFile),           INTENT(IN   ) :: InputFileData  !< Data stored in the module's input file -- intent(out) only for move_alloc statements
    TYPE(AA_ParameterType),       INTENT(INOUT) :: p              !< Parameters
    INTEGER(IntKi),               INTENT(  OUT) :: ErrStat        !< Error status of the operation
    CHARACTER(*),                 INTENT(  OUT) :: ErrMsg         ! Error message if ErrStat /= ErrID_None
    ! Local variables
    CHARACTER(ErrMsgLen)    :: ErrMsg2         ! temporary Error message if ErrStat /    = ErrID_None
    INTEGER(IntKi)          :: ErrStat2        ! temporary Error status of the operation
!    INTEGER(IntKi)          :: simcou,coun     ! simple loop  counter
    INTEGER(IntKi)          :: I,J,whichairfoil,K,i1_1,i10_1,i1_2,i10_2,iLE
    character(*), parameter :: RoutineName = 'SetParameters'
    LOGICAL                 :: tri,LE_flag
    REAL(ReKi)              :: val1,val10,f2,f4,lefttip,rightip,jumpreg, dist1, dist10
    ! Initialize variables for this routine
    ErrStat  = ErrID_None
    ErrMsg   = ""
    !!Assign input fiel data to parameters
    p%DT               = InputFileData%DT_AA         ! seconds
    p%AA_Bl_Prcntge    = InputFileData%AA_Bl_Prcntge  ! %
    p%fsample          = 1/p%DT     ! Hz
    p%total_sample     = 2**( ceiling(log(1*p%fsample)/log(2.0d0)))! 1 stands for the 1 seconds. Every 1 second Vrel spectra will be calculated for the dissipation calculation (change if more needed & recompile )
    p%total_sampleTI   = 5/p%DT  ! 10 seconds for TI sampling
    p%AAStart          = InputFileData%AAStart
    p%IBLUNT           = InputFileData%IBLUNT
    p%ILAM             = InputFileData%ILAM
    p%ITIP             = InputFileData%ITIP
    p%ITRIP            = InputFileData%ITRIP
    p%ITURB            = InputFileData%ITURB
    p%IInflow          = InputFileData%IInflow
    p%X_BLMethod       = InputFileData%X_BLMethod
    p%TICalcMeth       = InputFileData%TICalcMeth
    p%AweightFlag      = InputFileData%AweightFlag
    p%ROUND            = InputFileData%ROUND
    p%alprat           = InputFileData%ALPRAT
    p%NrOutFile        = InputFileData%NrOutFile
    p%delim            = Tab 
    p%outFmt           = "ES15.6E3" 
    p%NumBlNds         = InitInp%NumBlNds
    p%AirDens          = InitInp%AirDens          
    p%KinVisc          = InitInp%KinVisc
    p%SpdSound         = InitInp%SpdSound
    p%HubHeight        = InitInp%HubHeight
    p%Lturb            = InputFileData%Lturb
    p%dy_turb_in       = InputFileData%dy_turb_in
    p%dz_turb_in       = InputFileData%dz_turb_in
    p%NrObsLoc         = InputFileData%NrObsLoc
    p%FTitle           = InputFileData%FTitle

    IF ((InputFileData%TICalcMeth==1)) THEN
        call AllocAry(p%TI_Grid_In,size(InputFileData%TI_Grid_In,1), size(InputFileData%TI_Grid_In,2),  'p%TI_Grid_In', errStat2, errMsg2); if(Failed()) return
        p%TI_Grid_In=InputFileData%TI_Grid_In
    ENDIF

    p%AvgV=InputFileData%AvgV

    ! Copy AFInfo into AA module
    ! TODO Allocate AFInfo   and AFindx variables (DONE AND DONE) 
    ALLOCATE(p%AFInfo( size(InitInp%AFInfo) ), STAT=ErrStat2)
    IF ( ErrStat2 /= 0 )  THEN
        CALL SetErrStat(ErrID_Fatal, 'Error allocating memory for the InitInp%AFInfo array.', ErrStat2, ErrMsg2, RoutineName)
        RETURN
    ENDIF

    do i=1,size(InitInp%AFInfo)
        call AFI_CopyParam(InitInp%AFInfo(i), p%AFInfo(i), MESH_NEWCOPY, errStat2, errMsg2); if(Failed()) return
    end do

    ! Check 1
    tri=.true.
    IF( (p%ITURB.eq.2) .or. (p%IInflow.gt.1) )then
        ! if tno is on or one of the guidati models is on, check if we have airfoil coordinates
        DO k=1,size(p%AFInfo) ! if any of the airfoil coordinates are missing change calculation method
            IF( (size(p%AFInfo(k)%X_Coord) .lt. 5) .or. (size(p%AFInfo(k)%Y_Coord).lt.5) )then
                IF (tri) then ! Print the message for once only
                    CALL WrScr( 'Airfoil coordinates are missing: If Full or Simplified Guidati or Bl Calculation is on coordinates are needed ' )
                    CALL WrScr( 'Calculation methods enforced as BPM for TBLTE and only Amiet for inflow ' )
                    p%ITURB   = 1
                    p%IInflow = 1
                    tri=.false.
                ENDIF
            ENDIF
        ENDDO
    ENDIF
    
    ! Check 2
    ! if passed the first check and if tno, turn on boundary layer calculation
    IF( (p%ITURB.eq.2)) then
        p%X_BLMethod=X_BLMethod_Tables
    ENDIF
    
    ! Check 3
    ! if boundary layer is tripped then laminar b.l. vortex shedding mechanism is turned off
    IF( p%ITRIP.gt.0  )then
        p%ILAM=0
    ENDIF

    ! set 1/3 octave band frequency as parameter and A weighting.
    CALL AllocAry( p%FreqList, 34, 'FreqList', ErrStat2, ErrMsg2); if(Failed()) return
    p%FreqList = (/10.,12.5,16.,20.,25.,31.5,40.,50.,63.,80., &
        100.,125.,160.,200.,250.,315.,400.,500.,630.,800., & 
        1000.,1250.,1600.,2000.,2500.,3150.,4000.,5000.,6300.,8000., &
        10000.,12500.,16000.,20000./)  ! TODO this should be fortran parameter


    CALL AllocAry(p%Aweight, size(p%Freqlist), 'Aweight', ErrStat2, ErrMsg2); if(Failed()) return
    Do I=1,size(p%Freqlist)
        f2 = p%Freqlist(I)**2;
        f4 = p%Freqlist(I)**4;
        p%Aweight(I)=  10 * log(1.562339 * f4 / ((f2 + 107.65265**2) &
            * (f2 + 737.86223 **2))) / log(10.0_Reki) &
            + 10 * log(2.242881E+16 * f4 / ((f2 + 20.598997**2)**2 &
            * (f2 + 12194.22**2)**2)) / log(10.0_Reki)
    enddo

    ! Observer Locations
    call AllocAry(p%ObsX,  p%NrObsLoc, 'p%ObsX', ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(p%ObsY,  p%NrObsLoc, 'p%ObsY', ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(p%ObsZ,  p%NrObsLoc, 'p%ObsZ', ErrStat2, ErrMsg2); if(Failed()) return
    p%ObsX = InputFileData%ObsX
    p%ObsY = InputFileData%ObsY
    p%ObsZ = InputFileData%ObsZ
    ! 
    call AllocAry(p%BlAFID,      p%NumBlNds, p%numBlades, 'p%BlAFID' , ErrStat2, ErrMsg2); if(Failed()) return
    p%BlAFID=InitInp%BlAFID
    
    ! Blade Characteristics chord,span,trailing edge angle and thickness,airfoil ID for each segment
    call AllocAry(p%TEThick   ,p%NumBlNds,p%NumBlades,'p%TEThick'   ,ErrStat2,ErrMsg2); if(Failed()) return
    call AllocAry(p%TEAngle   ,p%NumBlNds,p%NumBlades,'p%TEAngle'   ,ErrStat2,ErrMsg2); if(Failed()) return
    call AllocAry(p%StallStart,p%NumBlNds,p%NumBlades,'p%StallStart',ErrStat2,ErrMsg2); if(Failed()) return
    p%StallStart(:,:) = 0.0_ReKi

     do i=1,p%NumBlades    
        do j=1,p%NumBlNds
            whichairfoil = p%BlAFID(j,i)
            p%TEThick(j,i) = InputFileData%BladeProps(whichairfoil)%TEThick
            p%TEAngle(j,i) = InputFileData%BladeProps(whichairfoil)%TEAngle
            
            if(p%AFInfo(whichairfoil)%NumTabs /=1 ) then
                call SetErrStat(ErrID_Fatal, 'Number of airfoil tables within airfoil file different than 1, which is not supported.', ErrStat2, ErrMsg2, RoutineName )
                if(Failed()) return
            endif
            p%StallStart(j,i)  = p%AFInfo(whichairfoil)%Table(1)%UA_BL%alpha1*180/PI ! approximate stall angle of attack [deg] (alpha1 in [rad])
        enddo
    enddo

    call AllocAry(p%BlSpn,       p%NumBlNds, p%NumBlades, 'p%BlSpn'  , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(p%BlChord,     p%NumBlNds, p%NumBlades, 'p%BlChord', ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(p%AerCent,  2, p%NumBlNds, p%NumBlades, 'p%AerCent', ErrStat2, ErrMsg2); if(Failed()) return
    p%BlSpn   = InitInp%BlSpn
    p%BlChord = InitInp%BlChord

    do j=p%NumBlNds,2,-1
        IF ( p%BlSpn(j,1) .lt. p%BlSpn(p%NumBlNds,1)*(100-p%AA_Bl_Prcntge)/100 )THEN ! assuming 
            p%startnode=j
            exit ! exit the loop
        endif
    enddo

    IF (p%startnode.lt.2) THEN
        p%startnode=2
    ENDIF

    !print*, 'AeroAcoustics Module is using the blade nodes starting from ' ,p%startnode,' Radius in meter ',p%BlSpn(p%startnode,1)
    !AerodYnamic center extraction for each segment 
    do i=1,p%numBlades
        do j=1,p%NumBlNds
            whichairfoil         = p%BlAFID(j,i)  ! just a temporary variable for clear coding
            ! airfoil coordinates read by AeroDyn. First value is the aerodynamic center
            p%AerCent(1,J,I)  = p%AFInfo(whichairfoil)%X_Coord(1)  ! assigned here corresponding airfoil.
            p%AerCent(2,J,I)  = p%AFInfo(whichairfoil)%Y_Coord(1)  ! assigned here corresponding airfoil.
        enddo
    enddo

    ! Dimensionalize Leading and trailing edge coordinates for later usage
    call AllocAry( p%AFTeCo, 3, p%NumBlNds,p%numBlades, 'p%AFTeCo', errStat2, errMsg2 ); if(Failed())return
    call AllocAry( p%AFLeCo, 3, p%NumBlNds,p%numBlades, 'p%AFLeCo', errStat2, errMsg2 ); if(Failed())return
    p%AFTeCo=0.0_Reki
    p%AFLeCo=0.0_Reki

    ! Normalized Leading edge coordinates (0,0,0)  
    ! Normalized Trailing edge coordinates  (1,0,0) -- > changed to 0,1,0 
    DO i=1,p%numBlades 
        DO j=1,p%NumBlNds
            p%AFLeCo(1,j,i) = ( 0.0_Reki -  p%AerCent(2,J,I)  ) * p%BlChord(j,i) ! (y_LE - y_AC) *Chord
            p%AFLeCo(2,j,i) = ( 0.0_Reki -  p%AerCent(1,J,I)  ) * p%BlChord(j,i) ! (x_LE - x_AC) *Chord
            p%AFLeCo(3,j,i) = ( 0.0_Reki -       0.0_Reki     ) * p%BlChord(j,i) ! this is always zero at the moment ( kept for 3d consistency )
            p%AFTeCo(1,j,i) = ( 0.0_Reki -  p%AerCent(2,J,I)  ) * p%BlChord(j,i) ! (y_TE - y_AC) *Chord
            p%AFTeCo(2,j,i) = ( 1.0_Reki -  p%AerCent(1,J,I)  ) * p%BlChord(j,i) ! (x_TE - x_AC) *Chord
            p%AFTeCo(3,j,i) = ( 0.0_Reki -       0.0_Reki     ) * p%BlChord(j,i) ! this is always zero at the moment  ( kept for 3d consistency )
        ENDDO
    ENDDO

    if (p%X_BLMethod .eq. X_BLMethod_Tables) then

        ! Copying inputdata list of AOA and Reynolds to parameters
        call AllocAry( p%AOAListBL, size(InputFileData%AOAListBL), 'p%AOAListBL', errStat2, errMsg2); if(Failed()) return
        call AllocAry( p%ReListBL,  size(InputFileData%ReListBL) , 'p%ReListBL' , errStat2, errMsg2); if(Failed()) return
        p%AOAListBL=InputFileData%AOAListBL
        p%ReListBL=InputFileData%ReListBL
        ! Allocate the suction and pressure side boundary layer parameters for output - will be used as tabulated data
        call AllocAry(p%dstarall1  ,size(p%AOAListBL), size(p%ReListBL),size(p%AFInfo),'p%dstarall1'  , errStat2, errMsg2); if(Failed()) return
        call AllocAry(p%dstarall2  ,size(p%AOAListBL), size(p%ReListBL),size(p%AFInfo),'p%dstarall2'  , errStat2, errMsg2); if(Failed()) return
        call AllocAry(p%d99all1    ,size(p%AOAListBL), size(p%ReListBL),size(p%AFInfo),'p%d99all1'    , errStat2, errMsg2); if(Failed()) return
        call AllocAry(p%d99all2    ,size(p%AOAListBL), size(p%ReListBL),size(p%AFInfo),'p%d99all2'    , errStat2, errMsg2); if(Failed()) return
        call AllocAry(p%Cfall1     ,size(p%AOAListBL), size(p%ReListBL),size(p%AFInfo),'p%Cfall1'     , errStat2, errMsg2); if(Failed()) return
        call AllocAry(p%Cfall2     ,size(p%AOAListBL), size(p%ReListBL),size(p%AFInfo),'p%Cfall2'     , errStat2, errMsg2); if(Failed()) return
        call AllocAry(p%EdgeVelRat1,size(p%AOAListBL), size(p%ReListBL),size(p%AFInfo),'p%EdgeVelRat1', errStat2, errMsg2); if(Failed()) return
        call AllocAry(p%EdgeVelRat2,size(p%AOAListBL), size(p%ReListBL),size(p%AFInfo),'p%EdgeVelRat2', errStat2, errMsg2); if(Failed()) return
        p%dstarall1   =0.0_ReKi
        p%dstarall2   =0.0_ReKi
        p%d99all1     =0.0_ReKi
        p%d99all2     =0.0_ReKi
        p%Cfall1      =0.0_ReKi
        p%Cfall2      =0.0_ReKi
        p%EdgeVelRat1 =0.0_ReKi
        p%EdgeVelRat2 =0.0_ReKi


        ! --- BL data are read from files and just copy what was read from the files
        p%dstarall1   = InputFileData%Suct_DispThick
        p%dstarall2   = InputFileData%Pres_DispThick
        p%d99all1     = InputFileData%Suct_BLThick
        p%d99all2     = InputFileData%Pres_BLThick
        p%Cfall1      = InputFileData%Suct_Cf
        p%Cfall2      = InputFileData%Pres_Cf
        p%EdgeVelRat1 = InputFileData%Suct_EdgeVelRat
        p%EdgeVelRat2 = InputFileData%Pres_EdgeVelRat
        
        if(Failed()) return
    endif

    ! If simplified guidati is on, calculate the airfoil thickness at 1% and at 10% chord from input airfoil coordinates
    IF (p%IInflow .EQ. 2) THEN
        call AllocAry(p%AFThickGuida,2,size(p%AFInfo),  'p%AFThickGuida', errStat2, errMsg2); if(Failed()) return
        p%AFThickGuida=0.0_Reki

        DO k=1,size(p%AFInfo) ! for each airfoil interpolation 

            ! IF ((MIN(p%AFInfo(k)%X_Coord) < 0.) .or. (MAX(p%AFInfo(k)%X_Coord) > 0.)) THEN
            !     call SetErrStat ( ErrID_Fatal,'The coordinates of airfoil '//trim(num2lstr(k))//' are mot defined between x=0 and x=1. Code stops.' ,ErrStat, ErrMsg, RoutineName )
            ! ENDIF
            
            ! Flip the flag when LE is found and find index
            LE_flag = .False.
            DO i=3,size(p%AFInfo(k)%X_Coord)
                IF (LE_flag .eqv. .False.) THEN
                    IF (p%AFInfo(k)%X_Coord(i) - p%AFInfo(k)%X_Coord(i-1) > 0.) THEN
                        LE_flag = .TRUE.
                        iLE = i
                    ENDIF
                ENDIF
            ENDDO

            ! From LE toward TE
            dist1  = ABS( p%AFInfo(k)%X_Coord(iLE) - 0.01)
            dist10 = ABS( p%AFInfo(k)%X_Coord(iLE) - 0.10)
            DO i=iLE+1,size(p%AFInfo(k)%X_Coord)
                IF (ABS(p%AFInfo(k)%X_Coord(i) - 0.01) < dist1) THEN
                    i1_1 = i
                    dist1 = ABS(p%AFInfo(k)%X_Coord(i) - 0.01)
                ENDIF
                IF (ABS(p%AFInfo(k)%X_Coord(i) - 0.1) < dist10) THEN
                    i10_1 = i
                    dist10 = ABS(p%AFInfo(k)%X_Coord(i) - 0.1)
                ENDIF
            ENDDO

            ! From TE to LE
            dist1  = 0.99
            dist10 = 0.90
            DO i=1,iLE-1
                IF (ABS(p%AFInfo(k)%X_Coord(i) - 0.01) < dist1) THEN
                    i1_2 = i
                    dist1 = ABS(p%AFInfo(k)%X_Coord(i) - 0.01)
                ENDIF
                IF (ABS(p%AFInfo(k)%X_Coord(i) - 0.1) < dist10) THEN
                    i10_2 = i
                    dist10 = ABS(p%AFInfo(k)%X_Coord(i) - 0.1)
                ENDIF
            ENDDO

            val1  = p%AFInfo(k)%Y_Coord(i1_1) - p%AFInfo(k)%Y_Coord(i1_2)
            val10 = p%AFInfo(k)%Y_Coord(i10_1) - p%AFInfo(k)%Y_Coord(i10_2)

            p%AFThickGuida(1,k)=val1  ! 1  % chord thickness
            p%AFThickGuida(2,k)=val10 ! 10  % chord thickness
        ENDDO
    ENDIF

    !! for turbulence intensity calculations on the fly every 5 meter the whole rotor area is divided vertically to store flow fields in each region
    jumpreg=7
    p%toptip           = CEILING(p%HubHeight+maxval(p%BlSpn(:,1)))+2 !Top Tip Height = Hub height plus radius
    p%bottip           = FLOOR(p%HubHeight-maxval(p%BlSpn(:,1)))-2 !Bottom Tip Height = Hub height minus radius
    call AllocAry(p%rotorregionlimitsVert,ceiling(((p%toptip)-(p%bottip))/jumpreg),  'p%rotorregionlimitsVert', errStat2, errMsg2); if(Failed()) return
    do i=0,size(p%rotorregionlimitsVert)-1
        p%rotorregionlimitsVert(i+1)=(p%bottip)+jumpreg*i
    enddo
    !! for turbulence intensity calculations on the fly every 5 meter the whole rotor area is divided horizontally to store flow fields in each region
    jumpreg=7
    lefttip           = 2*maxval(p%BlSpn(:,1))+5 !
    rightip           = 0 !
    call AllocAry( p%rotorregionlimitsHorz,ceiling(((lefttip)-(rightip))/jumpreg),  'p%rotorregionlimitsHorz', errStat2, errMsg2); if(Failed()) return
    do i=0,size(p%rotorregionlimitsHorz)-1
        p%rotorregionlimitsHorz(i+1)=rightip+jumpreg*i
    enddo
    jumpreg=60 ! 10 ! must be divisable to 360
    call AllocAry(p%rotorregionlimitsalph,INT((360/jumpreg)+1),  'p%rotorregionlimitsalph', errStat2, errMsg2); if(Failed()) return
    do i=0,size(p%rotorregionlimitsalph)-1
        p%rotorregionlimitsalph(i+1)=jumpreg*i
    enddo
    jumpreg=5
    call AllocAry( p%rotorregionlimitsrad, (CEILING( maxval(p%BlSpn(:,1))/jumpreg )+2),  'p%rotorregionlimitsrad', errStat2, errMsg2); if(Failed()) return
    do i=1,size(p%rotorregionlimitsrad)-1
        p%rotorregionlimitsrad(i+1)=jumpreg*i
    enddo
    p%rotorregionlimitsrad(1)=0.0_reki
    p%rotorregionlimitsrad(size(p%rotorregionlimitsrad)-1)=p%rotorregionlimitsrad(size(p%rotorregionlimitsrad)-1)+3

contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
    end function Failed
end subroutine SetParameters
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes AeroAcoustics module input array variables for use during the simulation.
subroutine Init_u( u, p, errStat, errMsg )
   type(AA_InputType),           intent(  out)  :: u                 !< Input data
   type(AA_ParameterType),       intent(in   )  :: p                 !< Parameters
   integer(IntKi),               intent(  out)  :: errStat           !< Error status of the operation
   character(*),                 intent(  out)  :: errMsg            !< Error message if ErrStat /= ErrID_None
   !local variables
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_u'

   call AllocAry(u%AoANoise  , p%NumBlNds, p%numBlades, 'u%AoANoise', errStat2      , errMsg2); if(Failed()) return
   call AllocAry(u%Vrel      , p%NumBlNds, p%numBlades, 'u%Vrel'    , errStat2      , errMsg2); if(Failed()) return
   call AllocAry(u%AeroCent_G, 3         , p%NumBlNds , p%numBlades , 'u%AeroCent_G', errStat2    , errMsg2); if(Failed()) return
   call AllocAry(u%Inflow    , 3_IntKi   , p%NumBlNds , p%numBlades , 'u%Inflow'    , ErrStat2    , ErrMsg2); if(Failed()) return
   call AllocAry(u%RotGtoL   , 3         , 3          , p%NumBlNds  , p%numBlades   , 'u%RotGtoL' , errStat2  , errMsg2); if(Failed()) return
   u%AoANoise   = 0.0_Reki
   u%Vrel       = 0.0_Reki
   u%RotGtoL    = 0.0_Reki
   u%AeroCent_G = 0.0_Reki
   u%Inflow     = 0.0_Reki
contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
    end function Failed
end subroutine Init_u
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes AeroAcoustics  output array variables for use during the simulation.
subroutine Init_y(y, u, p, errStat, errMsg)
    type(AA_OutputType),           intent(  out)  :: y               !< Module outputs
    type(AA_InputType),            intent(inout)  :: u               !< Module inputs -- intent(out) because of mesh sibling copy
    type(AA_ParameterType),        intent(inout)  :: p               !< Parameters
    integer(IntKi),                intent(  out)  :: errStat         !< Error status of the operation
    character(*),                  intent(  out)  :: errMsg          !< Error message if ErrStat /= ErrID_None
    ! Local variables
    integer(intKi)                               :: ErrStat2          ! temporary Error status
    character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
    character(*), parameter                      :: RoutineName = 'Init_y'
    integer(intKi)                               :: nNoiseMechanism   ! loop counter for blades
    ! Initialize variables for this routine
    errStat = ErrID_None
    errMsg  = ""
    nNoiseMechanism = 7! 7 noise mechanisms
    p%numOuts       = p%NrObsLoc
    p%NumOutsForSep = p%NrObsLoc*size(p%FreqList)*nNoiseMechanism
    p%NumOutsForPE  = p%NrObsLoc*size(p%Freqlist)
    p%NumOutsForNodes = p%NrObsLoc*p%NumBlNds*p%NumBlades
    call AllocAry(y%WriteOutput        , p%numOuts              , 'y%WriteOutput'           , errStat2                   , errMsg2); if(Failed()) return
    call AllocAry(y%WriteOutputSep     , p%NumOutsForSep        , 'y%WriteOutputSep'        , errStat2                   , errMsg2); if(Failed()) return
    call AllocAry(y%WriteOutputForPE   , p%numOutsForPE         , 'y%WriteOutputForPE'      , errStat2                   , errMsg2); if(Failed()) return
    call AllocAry(y%DirectiviOutput    , p%NrObsLoc             , 'y%DirectiviOutput'       , errStat2                   , errMsg2); if(Failed()) return
    call AllocAry(y%WriteOutputNode    , p%NumOutsForNodes      , 'y%WriteOutputSepFreq' , errStat2  , errMsg2); if(Failed()) return
    call AllocAry(y%OASPL              , p%NrObsLoc             , p%NumBlNds                , p%NumBlades                , 'y%OASPL'           , errStat2               , errMsg2); if(Failed()) return
    call AllocAry(y%SumSpecNoise       , size(p%FreqList)       , p%NrObsLoc                , p%NumBlades                , 'y%SumSpecNoise'    , errStat2               , errMsg2); if(Failed()) return
    call AllocAry(y%SumSpecNoiseSep    , 7                      , p%NrObsLoc                , size(p%FreqList)           , 'y%SumSpecNoiseSep' , errStat2               , errMsg2); if(Failed()) return
    call AllocAry(y%OASPL_Mech         , nNoiseMechanism        , p%NrObsLoc                , p%NumBlNds                 , p%NumBlades         , 'y%OASPL_Mech'         , errStat2  , errMsg2); if(Failed()) return
    call AllocAry(y%OutLECoords        , 3                      , size(p%FreqList)          , p%NrObsLoc                 , p%NumBlades         , 'y%OutLECoords'        , errStat2  , errMsg2); if(Failed()) return
    call AllocAry(y%PtotalFreq         , p%NrObsLoc             , size(p%FreqList)          , 'y%PtotalFreq'             , errStat2            , errMsg2); if(Failed()) return

    y%WriteOutput        = 0.0_reki
    y%WriteOutputSep     = 0.0_reki
    y%WriteOutputForPE   = 0.0_reki
    y%DirectiviOutput    = 0.0_reki
    y%WriteOutputNode    = 0.0_reki
    y%OASPL              = 0.0_reki
    y%OASPL_Mech         = 0.0_reki
    y%SumSpecNoise       = 0.0_reki
    y%SumSpecNoiseSep    = 0.0_reki
    y%OutLECoords        = 0.0_reki
    y%PtotalFreq         = 0.0_reki

contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
    end function Failed
end subroutine Init_y
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes (allocates) the misc variables for use during the simulation.
subroutine Init_MiscVars(m, p, u, y, errStat, errMsg)
    type(AA_MiscVarType),          intent(inout)  :: m                !< misc/optimization data (not defined in submodules)
    type(AA_ParameterType),        intent(in   )  :: p                !< Parameters
    type(AA_InputType),            intent(inout)  :: u                !< input for HubMotion mesh (create sibling mesh here)
    type(AA_OutputType),           intent(in   )  :: y                !< output (create mapping between output and otherstate mesh here)
    integer(IntKi),                intent(  out)  :: errStat          !< Error status of the operation
    character(*),                  intent(  out)  :: errMsg           !< Error message if ErrStat /= ErrID_None
    ! Local variables
    integer(intKi)                               :: ErrStat2          ! temporary Error status
    character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
    character(*), parameter                      :: RoutineName = 'Init_MiscVars'
    ! Initialize variables for this routine
    errStat = ErrID_None
    errMsg  = ""
    call AllocAry(m%ChordAngleLE, p%NrObsLoc, p%NumBlNds, p%numBlades, 'ChordAngleLE', ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(m%SpanAngleLE , p%NrObsLoc, p%NumBlNds, p%numBlades, 'SpanAngleLE' , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(m%ChordAngleTE, p%NrObsLoc, p%NumBlNds, p%numBlades, 'ChordAngleTE', ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(m%SpanAngleTE , p%NrObsLoc, p%NumBlNds, p%numBlades, 'SpanAngleTE' , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(m%rTEtoObserve, p%NrObsLoc, p%NumBlNds, p%numBlades, 'rTEtoObserve', ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(m%rLEtoObserve, p%NrObsLoc, p%NumBlNds, p%numBlades, 'rLEtoObserve', ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(m%SPLLBL      , size(p%FreqList), 'SPLLBL'    , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLP        , size(p%FreqList), 'SPLP'      , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLS        , size(p%FreqList), 'SPLS'      , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLALPH     , size(p%FreqList), 'SPLALPH'   , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLTBL      , size(p%FreqList), 'SPLTBL'    , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLBLUNT    , size(p%FreqList), 'SPLBLUNT'  , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLTIP      , size(p%FreqList), 'SPLTIP'    , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLTI       , size(p%FreqList), 'SPLTI'     , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLTIGui    , size(p%FreqList), 'SPLTIGui'  , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%CfVar       , 2               , 'CfVar'     , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%d99Var      , 2               , 'd99Var'    , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%dstarVar    , 2               , 'dstarVar'  , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%EdgeVelVar  , 2               , 'EdgeVelVar', errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%LE_Location,  3, p%NumBlNds, p%numBlades, 'LE_Location', ErrStat2, ErrMsg2); if(Failed()) return
    m%ChordAngleLE = 0.0_ReKi
    m%SpanAngleLE  = 0.0_ReKi
    m%ChordAngleTE = 0.0_ReKi
    m%SpanAngleTE  = 0.0_ReKi
    m%rTEtoObserve = 0.0_ReKi
    m%rLEtoObserve = 0.0_ReKi
    m%SPLLBL       = 0.0_ReKi
    m%SPLP         = 0.0_ReKi
    m%SPLS         = 0.0_ReKi
    m%SPLALPH      = 0.0_ReKi
    m%SPLTBL       = 0.0_ReKi
    m%SPLBLUNT     = 0.0_ReKi
    m%SPLTIP       = 0.0_ReKi
    m%SPLTI        = 0.0_ReKi
    m%SPLTIGui     = 0.0_ReKi
    m%CfVar        = 0.0_ReKi
    m%d99Var       = 0.0_ReKi
    m%dstarVar     = 0.0_ReKi
    m%EdgeVelVar   = 0.0_ReKi
    m%LE_Location  = 0.0_ReKi
    m%speccou      = 0
    m%filesopen    = 0
contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
    end function Failed
end subroutine Init_MiscVars
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes (allocates) the misc variables for use during the simulation.
subroutine Init_states(xd, p, errStat, errMsg)
    type(AA_DiscreteStateType),    intent(inout)  :: xd               !
    type(AA_ParameterType),        intent(in   )  :: p                !< Parameters
    integer(IntKi),                intent(  out)  :: errStat          !< Error status of the operation
    character(*),                  intent(  out)  :: errMsg           !< Error message if ErrStat /= ErrID_None
    ! Local variables
    integer(intKi)                               :: k,ji
    integer(intKi)                               :: ErrStat2          ! temporary Error status
    character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
    character(*), parameter                      :: RoutineName = 'Init_DiscrStates'
    ! Initialize variables for this routine
    errStat = ErrID_None
    errMsg  = ""

    call AllocAry(xd%MeanVrel,   p%NumBlNds, p%numBlades, 'xd%MeanVrel'  , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(xd%VrelSq,     p%NumBlNds, p%numBlades, 'xd%VrelSq'    , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(xd%TIVrel,     p%NumBlNds, p%numBlades, 'xd%TIVrel'    , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(xd%MeanVxVyVz, p%NumBlNds, p%numBlades, 'xd%MeanVxVyVz', ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(xd%TIVx,       p%NumBlNds, p%numBlades, 'xd%TIVx'      , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(xd%VxSq,       p%NumBlNds, p%numBlades, 'xd%VxSq'      , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(xd%VrelStore,  p%total_sample+1, p%NumBlNds, p%numBlades,'xd%VrelStore', ErrStat2, ErrMsg2) ! plus one just in case
    if(Failed()) return
    DO ji=1,size(xd%MeanVrel,2)
        DO k=1,size(xd%MeanVrel,1)
            xd%VrelSq (k,ji)     = 0.0_ReKi  ! Relative Velocity Squared for TI calculation (on the fly)
            xd%MeanVrel (k,ji)   = 0.0_ReKi  ! Relative Velocity Mean  calculation (on the fly)
            xd%TIVrel(k,ji)      = 0.0_ReKi  ! Turbulence Intensity (for on the fly calculation)
            xd%MeanVxVyVz (k,ji) = 0.0_ReKi  ! 
            xd%TIVx  (k,ji)      = 0.0_ReKi  ! 
            xd%VxSq (k,ji)       = 0.0_ReKi  !
            xd%VrelStore (1:size(xd%VrelStore,1),k,ji)       = 0.0_ReKi  !
        ENDDO
    ENDDO
    call AllocAry(xd%RegVxStor,p%total_sampleTI,size(p%rotorregionlimitsrad)-1,size(p%rotorregionlimitsalph)-1,'xd%Vxst',ErrStat2,ErrMsg2)
    if(Failed()) return
    call AllocAry(xd%allregcounter ,size(p%rotorregionlimitsrad)-1,size(p%rotorregionlimitsalph)-1,'xd%allregcounter',ErrStat2,ErrMsg2 )
    if(Failed()) return
    call AllocAry(xd%VxSqRegion    ,size(p%rotorregionlimitsrad)-1,size(p%rotorregionlimitsalph)-1,'xd%VxSqRegion'   , ErrStat2, ErrMsg2)
    if(Failed()) return
    call AllocAry(xd%RegionTIDelete,size(p%rotorregionlimitsrad)-1,size(p%rotorregionlimitsalph)-1,'xd%RegionTIDelete', ErrStat2, ErrMsg2)
    do ji=1,size(xd%allregcounter,2)
        do k=1,size(xd%allregcounter,1)
            xd%allregcounter(k,ji)      = 2.0_Reki  ! 
            xd%VxSqRegion(k,ji)         = 0.0_ReKi  ! 
            xd%RegionTIDelete(k,ji)     = 0.0_ReKi  ! 
            xd%RegVxStor(1:size(xd%RegVxStor,1),k,ji)=0.0_reki
        enddo
    enddo
contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
    end function Failed
end subroutine Init_states
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AA_UpdateStates( t, n, m, u, p,  xd,  errStat, errMsg )
   real(DbKi),                     intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),                 intent(in   ) :: n          !< Current simulation time step n = 0,1,...
   type(AA_InputType),             intent(in   ) :: u          !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   TYPE(AA_ParameterType),         INTENT(IN   ) :: p          !< Parameters
   type(AA_DiscreteStateType),     intent(inout) :: xd         !< Input: Discrete states at t;
   type(AA_MiscVarType),           intent(inout) :: m          !< misc/optimization data 
   integer(IntKi),                 intent(  out) :: errStat    !< Error status of the operation
   character(*),                   intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None
   ! local variables
!   integer(intKi)                               :: ErrStat2          ! temporary Error status
!   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'AA_UpdateStates'
   REAL(ReKi),DIMENSION(p%NumBlNds,p%numBlades) :: TEMPSTD  ! temporary standard deviation variable
   REAL(ReKi)                                   :: tempsingle,tempmean,angletemp,abs_le_x,ti_vx,U1,U2   ! temporary standard deviation variable
   integer(intKi)                               :: i,j,k,rco, y0_a,y1_a,z0_a,z1_a
   REAL(ReKi) :: yi_a,zi_a,yd_a,zd_a,c00_a,c10_a

   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Cumulative mean and standard deviation, states are updated as Vx Vy Vz changes at each time step
   TEMPSTD       =  sqrt( u%Inflow(1,:,:)**2+u%Inflow(2,:,:)**2+u%Inflow(3,:,:)**2 )
   xd%MeanVxVyVz = (TEMPSTD + xd%MeanVxVyVz*n) / (n+1)  
   !   xd%VxSq       = TEMPSTD**2 + xd%VxSq
   !   TEMPSTD     = sqrt(  (xd%VxSq/(n+1)) - (xd%MeanVxVyVz**2)   )
   !   xd%TIVx  = (TEMPSTD / xd%MeanVxVyVz ) ! check inflow noise input for multiplication with 100 or not

   m%speccou= m%speccou+1
   IF(   (p%TICalcMeth.eq.2) ) THEN
       do i=1,p%NumBlades
           do j=1,p%NumBlNds
               abs_le_x=m%LE_Location(3,j,i)-p%hubheight
               IF ((abs_le_x.lt.0).and.(m%LE_Location(2,j,i).lt.0)) THEN
                   angletemp=180+ATAN(  ABS( m%LE_Location(2,j,i)/abs_le_x )  ) * R2D_D 
               ELSEIF ((abs_le_x.lt.0).and.(m%LE_Location(2,j,i).gt.0)) THEN
                   angletemp=180-ATAN(  ABS( m%LE_Location(2,j,i)/abs_le_x )  ) * R2D_D
               ELSEIF ((abs_le_x.gt.0).and.(m%LE_Location(2,j,i).lt.0)) THEN
                   angletemp=360-ATAN(  ABS( m%LE_Location(2,j,i)/abs_le_x )  ) * R2D_D 
               ELSEIF ((abs_le_x.gt.0).and.(m%LE_Location(2,j,i).gt.0)) THEN
                   angletemp=ATAN(   m%LE_Location(2,j,i)/abs_le_x  ) * R2D_D
               ELSE
                   CALL WrScr( 'problem in angletemp Aeroacoustics module' )
               ENDIF
               !abs_le_x=ABS(abs_le_x)
               do k=1,size(p%rotorregionlimitsrad)
                   IF (p%BlSpn(j,i)-p%rotorregionlimitsrad(k).lt.0) THEN ! it means location is in the k-1 region
                       !print*, abs_le_x,p%rotorregionlimitsrad(k),k-1
                       GOTO 4758
                   ENDIF
               enddo
               4758  do rco=1,size(p%rotorregionlimitsalph)
                   IF (angletemp-p%rotorregionlimitsalph(rco).lt.0) THEN ! it means location is in the k-1 region
                       GOTO 9815
                   ENDIF
               enddo
               9815 xd%allregcounter(k-1,rco-1)=CEILING(xd%allregcounter(k-1,rco-1)+1.0_Reki)    ! increase the sample amount in that specific 5 meter height vertical region
               tempsingle         = sqrt( u%Inflow(1,j,i)**2+u%Inflow(2,j,i)**2+u%Inflow(3,j,i)**2 )  ! 
               ! with storage region dependent moving average and TI
               IF  (INT(xd%allregcounter(k-1,rco-1)) .lt. (size(xd%RegVxStor,1)+1)) THEN
                   xd%RegVxStor(INT(xd%allregcounter(k-1,rco-1)),k-1,rco-1)=tempsingle
                   xd%TIVx(j,i)       =  0
                   xd%RegionTIDelete(k-1,rco-1)=0
               ELSE
                   xd%RegVxStor((mod(INT(xd%allregcounter(k-1,rco-1))-size(xd%RegVxStor,1),size(xd%RegVxStor,1)))+1,k-1,rco-1)=tempsingle
                   tempmean=SUM(xd%RegVxStor(:,k-1,rco-1))
                   tempmean=tempmean/size(xd%RegVxStor,1)
                   xd%RegionTIDelete(k-1,rco-1)=SQRT((SUM((xd%RegVxStor(:,k-1,rco-1)-tempmean)**2)) /  size(xd%RegVxStor,1) )
                   xd%TIVx(j,i)       =  xd%RegionTIDelete(k-1,rco-1) ! only the fluctuation 
               ENDIF
           enddo
       enddo
       
   ELSE! interpolate from the user given ti values
       do i=1,p%NumBlades
           do j=1,p%NumBlNds
               zi_a=ABS(m%LE_Location(3,j,i)  -  (FLOOR(p%HubHeight-maxval(p%BlSpn(:,1))))  )   /p%dz_turb_in
               z0_a=floor(zi_a)
               z1_a=ceiling(zi_a)
               zd_a=zi_a-z0_a  
               yi_a=ABS(m%LE_Location(2,j,i)  + maxval(p%BlSpn(:,1)) )  /p%dy_turb_in
               y0_a=floor(yi_a)
               y1_a=ceiling(yi_a)
               yd_a=yi_a-y0_a
               c00_a=(1.0_ReKi-yd_a)*p%TI_Grid_In(z0_a+1,y0_a+1)+yd_a*p%TI_Grid_In(z0_a+1,y1_a+1)
               c10_a=(1.0_ReKi-yd_a)*p%TI_Grid_In(z1_a+1,y0_a+1)+yd_a*p%TI_Grid_In(z1_a+1,y1_a+1)
               
               ! This is the turbulence intensity of the wind at the location of the blade i at node j
               ti_vx = (1.0_ReKi-zd_a)*c00_a+zd_a*c10_a
               ! With some velocity triangles, we convert it into the incident turbulence intensity, i.e. the TI used by the Amiet model
               U1 = u%Vrel(J,I) 
               U2 = SQRT((p%AvgV*(1.+ti_vx))**2 + U1**2 - p%AvgV**2)
               ! xd%TIVx(j,i)=(U2-U1)/U1
               xd%TIVx(j,i)=p%AvgV*ti_vx/U1
               
               
               if (i.eq.p%NumBlades) then 
                   if (j.eq.p%NumBlNds) then 
                   endif
               endif
           enddo
       enddo
   endif
end subroutine AA_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine AA_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
    TYPE(AA_InputType),           INTENT(INOUT)  :: u           !< System inputs
    TYPE(AA_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
    TYPE(AA_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
    TYPE(AA_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
    TYPE(AA_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
    TYPE(AA_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
    TYPE(AA_OutputType),          INTENT(INOUT)  :: y           !< System outputs
    TYPE(AA_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
    INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
    CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
    ! Initialize ErrStat
    ErrStat = ErrID_None
    ErrMsg  = ""
    ! Destroy the input data:
    CALL AA_DestroyInput( u, ErrStat, ErrMsg )
    ! Destroy the parameter data:
    CALL AA_DestroyParam( p, ErrStat, ErrMsg )
    ! Destroy the state data:
    CALL AA_DestroyContState(   x,           ErrStat, ErrMsg )
    CALL AA_DestroyDiscState(   xd,          ErrStat, ErrMsg )
    CALL AA_DestroyConstrState( z,           ErrStat, ErrMsg )
    CALL AA_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
    CALL AA_DestroyMisc(        m,           ErrStat, ErrMsg ) 
    ! Destroy the output data:
    CALL AA_DestroyOutput( y, ErrStat, ErrMsg )

END SUBROUTINE AA_End

!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine AA_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg)
    ! NOTE: no matter how many channels are selected for output, all of the outputs are calcalated
    ! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
    ! placed in the y%WriteOutput(:) array.
    !..................................................................................................................................
    REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
    TYPE(AA_InputType),           INTENT(IN   )  :: u           !< Inputs at Time t
    TYPE(AA_ParameterType),       INTENT(IN   )  :: p           !< Parameters
    TYPE(AA_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
    TYPE(AA_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
    TYPE(AA_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
    TYPE(AA_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
    TYPE(AA_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
    type(AA_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
    INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
    CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
    ! Local variables
    integer, parameter      :: indx = 1  ! m%BEMT_u(1) is at t; m%BEMT_u(2) is t+dt
    integer(intKi)          :: ErrStat2
    character(ErrMsgLen)    :: ErrMsg2
    character(*), parameter :: RoutineName = 'AA_CalcOutput'
    ErrStat = ErrID_None
    ErrMsg  = ""
    ! assume integer divide is possible
    call CalcObserve(t,p,m,u,xd,errStat2, errMsg2)
    call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
    IF (t >= p%AAStart) THEN        
        IF (mod(t + 1E-10,p%DT) .lt. 1E-6) THEN
            call CalcAeroAcousticsOutput(u,p,m,xd,y,errStat2,errMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

            call Calc_WriteOutput( p, u, m, y,  ErrStat2, ErrMsg2 )     
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      

            call AA_WriteOutputLine(y, t, p, ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
        ENDIF
    ENDIF
end subroutine AA_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE CalcObserve(t,p,m,u,xd,errStat,errMsg)
    REAL(DbKi),                          INTENT(IN   )  :: t      !< Current simulation time in seconds
    TYPE(AA_DiscreteStateType),          INTENT(IN   ) :: xd      !< discrete state type
    TYPE(AA_ParameterType),              intent(in   ) :: p       !< Parameters
    TYPE(AA_InputType),                  intent(in   ) :: u       !< NN Inputs at Time
    TYPE(AA_MiscVarType),                intent(inout) :: m       !< misc/optimization data (not defined in submodules)
    INTEGER(IntKi),                      intent(  out) :: errStat !< Error status of the operation
    CHARACTER(*),                        intent(  out) :: errMsg  !< Error message if ErrStat /= ErrID_None
    ! Local variables.
    REAL(ReKi)     :: RLEObserve (3)              ! Position vector from leading edge to observer in trailing edge coordinate system
    REAL(ReKi)     :: RTEObserve (3)              ! Position vector from trailing edge to observer in trailing edge coordinate system
    REAL(ReKi)     :: RTEObserveG (3)             ! Position vector from trailing edge to observer in the coordinate system located at the trailing edge and rotated as the global
    REAL(ReKi)     :: RLEObserveG (3)             ! Position vector from leading edge to observer in the coordinate system located at the leading edge and rotated as the global
    REAL(ReKi)     :: RTEObservereal (3)          ! Location of trailing edge in global coordinate system
    REAL(ReKi)     :: RLEObservereal (3)          ! Location of leading edge in global coordinate system
    REAL(ReKi)     :: LocalToGlobal(3,3)          ! Transformation matrix
    REAL(ReKi)     :: timeLE                      ! Time of sound propagation from leading edge to observer
    REAL(ReKi)     :: timeTE                      ! Time of sound propagation from trailing edge to observer
    REAL(ReKi)     :: phi_e                       ! Spanwise directivity angle
    REAL(ReKi)     :: theta_e                     ! Chordwise directivity angle 
    INTEGER(intKi) :: I                           ! I A generic index for DO loops.
    INTEGER(intKi) :: J                           ! J A generic index for DO loops.
    INTEGER(intKi) :: K                           ! K A generic index for DO loops.
!    INTEGER(intKi)               :: ErrStat2
!    CHARACTER(ErrMsgLen)         :: ErrMsg2
    CHARACTER(*), parameter      :: RoutineName = 'CalcObserveDist'

    ErrStat = ErrID_None
    ErrMsg  = ""
    ! Loop through the blades
    DO I = 1,p%numBlades
        ! Loop through the nodes along blade span
        DO J = 1,p%NumBlNds
            ! Transpose the rotational vector GlobalToLocal to obtain the rotation LocalToGlobal
            LocalToGlobal  = TRANSPOSE(u%RotGtoL(:,:,J,I))
            ! Rotate the coordinates of leading and trailing edge from the local reference system to the global. Then add the coordinates of the aerodynamic center in the global coordinate system
            ! The global coordinate system is located on the ground, has x pointing downwind, y pointing laterally, and z pointing vertically upwards
            RTEObservereal = MATMUL(LocalToGlobal, p%AFTeCo(:,J,I)) + u%AeroCent_G(:,J,I)
            RLEObservereal = MATMUL(LocalToGlobal, p%AFLeCo(:,J,I)) + u%AeroCent_G(:,J,I)
            ! Compute the coordinates of the leading edge in the global coordinate system
            m%LE_Location(1,J,I) = RLEObservereal(1)
            m%LE_Location(2,J,I) = RLEObservereal(2)
            m%LE_Location(3,J,I) = RLEObservereal(3)
            ! If the time step is set to generate AA outputs
            IF (t >= p%AAStart) THEN
                IF ( mod(t + 1E-10,p%DT) .lt. 1E-6)  THEN
                    ! Loop through the observers
                    DO K = 1,p%NrObsLoc
                        ! Calculate the position of the observer K in a reference system located at the trailing edge and oriented as the global reference system
                        RTEObserveG(1)=p%Obsx(K)-RTEObservereal(1)
                        RTEObserveG(2)=p%Obsy(K)-RTEObservereal(2)
                        RTEObserveG(3)=p%Obsz(K)-RTEObservereal(3)
                        ! Calculate the position of the observer K in a reference system located at the leading edge and oriented as the global reference system
                        RLEObserveG(1)=p%Obsx(K)-RLEObservereal(1)
                        RLEObserveG(2)=p%Obsy(K)-RLEObservereal(2)
                        RLEObserveG(3)=p%Obsz(K)-RLEObservereal(3)
                        ! Rotate back the two reference systems from global to local. 
                        RTEObserve = MATMUL(u%RotGtoL(:,:,J,I), RTEObserveG)
                        RLEObserve = MATMUL(u%RotGtoL(:,:,J,I), RLEObserveG)

                        ! Calculate absolute distance between node and observer
                        m%rTEtoObserve(K,J,I) = SQRT (RTEObserve(1)**2+RTEObserve(2)**2+RTEObserve(3)**2)
                        m%rLEtoObserve(K,J,I) = SQRT (RLEObserve(1)**2+RLEObserve(2)**2+RLEObserve(3)**2)

                        ! Calculate time of noise propagation to observer
                        timeTE = m%rTEtoObserve(K,J,I) / p%SpdSound
                        timeLE = m%rLEtoObserve(K,J,I) / p%SpdSound
                        
                        ! The local system has y alinged with the chord, x pointing towards the airfoil suction side, and z aligned with blade span from root towards tip 
                        ! x ---> z_e
                        ! y ---> x_e
                        ! z ---> y_e

                        ! Compute spanwise directivity angle phi for the trailing edge
                        phi_e = ATAN2 (RTEObserve(1) , RTEObserve(3))
                        m%SpanAngleTE(K,J,I)  = phi_e * R2D

                        ! Compute chordwise directivity angle theta for the trailing edge
                        theta_e = ATAN2 ((RTEObserve(3) * COS (phi_e) + RTEObserve(1) * SIN (phi_e) ) , RTEObserve(2))
                        m%ChordAngleTE(K,J,I) = theta_e * R2D
                        
                        ! Compute spanwise directivity angle phi  for the leading edge (it's the same angle for the trailing edge)
                        phi_e = ATAN2 (RLEObserve(1) , RLEObserve(3))
                        m%SpanAngleLE(K,J,I)  = phi_e * R2D

                        ! Compute chordwise directivity angle theta for the leading edge
                        theta_e = ATAN2 ((RLEObserve(3) * COS (phi_e) + RLEObserve(1) * SIN (phi_e) ) , RLEObserve(2))
                        m%ChordAngleLE(K,J,I) = theta_e * R2D

                    ENDDO !K, observers
                ENDIF !  every Xth time step or so..
            ENDIF ! only if the time step is more than user input value run this part
        ENDDO  !J, blade nodes
    ENDDO  !I , number of blades
END SUBROUTINE CalcObserve
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE CalcAeroAcousticsOutput(u,p,m,xd,y,errStat,errMsg)
    TYPE(AA_InputType),                 INTENT(IN   )       :: u       !< Inputs at Time t
    TYPE(AA_OutputType),                INTENT(INOUT)       :: y       !< 
    TYPE(AA_ParameterType),                   INTENT(IN   ) :: p       !< Parameters
    TYPE(AA_MiscVarType),                   INTENT(INOUT)   :: m       !< misc/optimization data (not defined in submodules)
    TYPE(AA_DiscreteStateType),             INTENT(IN   )   :: xd      !< discrete state type
    integer(IntKi),                           INTENT(  OUT) :: errStat !< Error status of the operation
    character(*),                             INTENT(  OUT) :: errMsg  !< Error message if ErrStat /= ErrID_None
    ! Local variables.
    integer(intKi)                :: III                                             !III A generic index for DO loops.
    integer(intKi)                :: I                                               !I   A generic index for DO loops.
    integer(intKi)                :: J                                               !J   A generic index for DO loops.
    integer(intKi)                :: K !,liop,cou ,JTEMP                                   !K   A generic index for DO loops.
    integer(intKi)                :: oi                                              !K   A generic index for DO loops.
    REAL(ReKi)                    :: AlphaNoise                                 ! 
    REAL(ReKi)                    :: UNoise                                     ! 
    REAL(ReKi)                    :: elementspan                                ! 
!    REAL(ReKi),DIMENSION(p%NumBlNds)       ::tempdel
!    REAL(ReKi),DIMENSION(p%NrObsLoc,p%NumBlNds,p%numBlades)    ::OASPLTBLAll
    REAL(ReKi),DIMENSION(p%NrObsLoc,p%NumBlNds,p%numBlades,size(p%FreqList))    ::ForMaxLoc
    REAL(ReKi),DIMENSION(size(y%OASPL_Mech,1),size(p%FreqList),p%NrObsLoc,p%NumBlNds,p%numBlades)    :: ForMaxLoc3
!    REAL(ReKi),DIMENSION(size(p%FreqList),p%NrObsLoc,p%numBlades)               ::SPL_Out
    REAL(ReKi),DIMENSION(p%NumBlNds,p%numBlades)    ::temp_dispthick
    REAL(ReKi),DIMENSION(p%NumBlNds,p%numBlades)    ::temp_dispthickchord

    real(ReKi)                                                 ::  Ptotal
    real(ReKi)                                                 :: PtotalLBL    
    real(ReKi)                                                 :: PtotalTBLP   
    real(ReKi)                                                 :: PtotalTBLS   
    real(ReKi)                                                 :: PtotalSep    
    real(ReKi)                                                 :: PtotalTBLAll 
    real(ReKi)                                                 :: PtotalBlunt  
    real(ReKi)                                                 :: PtotalTip    
    real(ReKi)                                                 :: PtotalInflow 
    real(ReKi)                                                 :: PLBL
    real(ReKi)                                                 :: PTBLP
    real(ReKi)                                                 :: PTBLS
    real(ReKi)                                                 :: PTBLALH
    real(ReKi)                                                 :: PTip
    real(ReKi)                                                 :: PTI
    real(ReKi)                                                 :: PBLNT !,adforma
!    REAL(ReKi),DIMENSION(2)                                    :: Cf ,d99, d_star
!    TYPE(FFT_DataType)                                         :: FFT_Data             !< the instance of the FFT module we're using
!    REAL(ReKi),DIMENSION(p%total_sample)                     :: spect_signal
!    REAL(ReKi),DIMENSION(p%total_sample/2)                   :: spectra
!    real(ReKi),ALLOCATABLE     ::  fft_freq(:)  
    integer(intKi)                                             :: ErrStat2
    character(ErrMsgLen)                                       :: ErrMsg2
    character(*), parameter                                    :: RoutineName = 'CalcAeroAcousticsOutput'

    ErrStat = ErrID_None
    ErrMsg  = ""

    !------------------- Fill with zeros -------------------------!
    DO I = 1,p%numBlades;DO J = 1,p%NumBlNds;DO K = 1,p%NrObsLoc; 
        y%OASPL(k,j,i)        = 0.0_Reki
        DO oi=1,size(y%OASPL_Mech,1)
            y%OASPL_Mech(oi,k,j,i)= 0.0_Reki
        ENDDO;
    ENDDO;ENDDO;ENDDO

    DO K = 1,p%NrObsLoc;     
        y%DirectiviOutput(K)  = 0.0_Reki
        DO I=1,p%NumBlades;DO III=1,size(p%FreqList);
            y%SumSpecNoise(III,K,I) = 0.0_Reki
            ForMaxLoc(K,1:p%NumBlNds,I,III)=0.0_Reki
            DO oi=1,size(y%OASPL_Mech,1)
                y%SumSpecNoiseSep(oi,K,III) = 0.0_Reki
                ForMaxLoc3(oi,III,K,1:p%NumBlNds,I)=0.0_Reki
                m%SPLLBL(III)=0.0_Reki
                m%SPLP(III)=0.0_Reki
                m%SPLS(III)=0.0_Reki
                m%SPLALPH(III)=0.0_Reki
                m%SPLBLUNT(III)=0.0_Reki
                m%SPLTIP(III)=0.0_Reki
                m%SPLti(III)=0.0_Reki
            ENDDO
        ENDDO;ENDDO
    ENDDO

    DO K = 1,p%NrObsLoc;
       DO III = 1,size(p%FreqList);
          y%PtotalFreq(K,III) = 0.0_ReKi
       ENDDO
    ENDDO

    !------------------- initialize FFT  -------------------------!
    !!!IF (m%speccou .eq. p%total_sample)THEN
    !!!CALL InitFFT ( p%total_sample, FFT_Data, ErrStat=ErrStat2 )
    !!! CALL SetErrStat(ErrStat2, 'Error in InitFFT', ErrStat, ErrMsg, 'CalcAeroAcousticsOutput' )
    !!!CALL AllocAry( fft_freq,  size(spect_signal)/2-1, 'fft_freq', ErrStat2, ErrMsg2 ) 
    !!!         CALL SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
    !!!do liop=1,size(fft_freq)
    !!!    fft_freq(liop)=p%fsample*liop ! fRequncy x axis
    !!!    fft_freq(liop)=fft_freq(liop)/size(spect_signal)
    !!!enddo
    !!!ENDIF


    
    DO I = 1,p%numBlades
        DO J = p%startnode,p%NumBlNds  ! starts loop from startnode. 
            !------------------------------!!------------------------------!!------------------------------!!------------------------------!
            !------------------------------!!------------------------------!!------------------------------!!------------------------------!
            !------------------------------!!------------------------------!!------------------------------!!------------------------------!
            !--------Calculate Spectrum for dissipation calculation-------------------------!
            !IF (m%speccou .eq. p%total_sample)THEN
            !spect_signal=xd%VrelStore(  1:p%total_sample,J,I  )
            !        CALL ApplyFFT_f( spect_signal, FFT_Data, ErrStat2 )
            ! IF (ErrStat2 /= ErrID_None ) THEN
            ! CALL SetErrStat(ErrStat2, 'Error in ApplyFFT .', ErrStat, ErrMsg, 'CalcAeroAcousticsOutput' )
            ! ENDIF
            !cou=1
            !O liop=2,size(spect_signal)-1,2
            !cou=cou+1
            !spectra(cou) = spect_signal(liop)*spect_signal(liop) + spect_signal(1+liop)*spect_signal(1+liop)
            !ENDDO
            !spectra(1)=spect_signal(1)*spect_signal(1)
            !spectra=spectra/(size(spectra)*2)
            !          m%speccou=0
            !ENDIF

            Unoise =  u%Vrel(J,I) 
            IF (EqualRealNos(Unoise,0.0_ReKi)) then
                Unoise = 0.1 ! TODO TODO a value consistent with the test above should be used
            ENDIF
            IF (J .EQ. p%NumBlNds) THEN
                elementspan =   (p%BlSpn(J,I)-p%BlSpn(J-1,I))/2 
            ELSE
                elementspan =   (p%BlSpn(J,I)-p%BlSpn(J-1,I))/2    +    (p%BlSpn(J+1,I)-p%BlSpn(J,I))/2  
            ENDIF
            AlphaNoise= u%AoANoise(J,I) * R2D_D 


            !--------Read in Boundary Layer Data-------------------------!
            IF (p%X_BLMethod .EQ. X_BLMethod_Tables) THEN
                call BL_Param_Interp(p,m,Unoise,AlphaNoise,p%BlChord(J,I),p%BlAFID(J,I), errStat2, errMsg2)
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                temp_dispthick(J,I) = m%d99Var(1)
                m%d99Var            = m%d99Var*p%BlChord(J,I)
                m%dstarVar          = m%dstarVar*p%BlChord(J,I)
                temp_dispthickchord(J,I)=m%d99Var(1)
            ENDIF

            !------------------------------!!------------------------------!!------------------------------!!------------------------------!
            !------------------------------!!------------------------------!!------------------------------!!------------------------------!
            !------------------------------!!------------------------------!!------------------------------!!------------------------------!
            DO K = 1,p%NrObsLoc 
                !--------Laminar Boundary Layer Vortex Shedding Noise----------------------------!
                IF ( (p%ILAM .EQ. 1) .AND. (p%ITRIP .EQ. 0) )    THEN
                    CALL LBLVS(AlphaNoise,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
                        elementspan,m%rTEtoObserve(K,J,I), &
                        p,m%d99Var(2),m%dstarVar(1),m%dstarVar(2),m%SPLLBL,p%StallStart(J,I),errStat2,errMsg2)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
                ENDIF
                !--------Turbulent Boundary Layer Trailing Edge Noise----------------------------!
                IF (   (p%ITURB .EQ. 1) .or. (p%ITURB .EQ. 2) )   THEN
                    CALL TBLTE(AlphaNoise,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
                    elementspan,m%rTEtoObserve(K,J,I), p, m%d99Var(2),m%dstarVar(1),m%dstarVar(2),p%StallStart(J,I), &
                    m%SPLP,m%SPLS,m%SPLALPH,m%SPLTBL,errStat2,errMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
                    IF (p%ITURB .EQ. 2)  THEN
                        m%SPLP=0.0_ReKi;m%SPLS=0.0_ReKi;m%SPLTBL=0.0_ReKi;
                        m%EdgeVelVar(1)=1.000d0;m%EdgeVelVar(2)=m%EdgeVelVar(1);
                        CALL TBLTE_TNO(UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
                            elementspan,m%rTEtoObserve(K,J,I),m%CfVar,m%d99var,m%EdgeVelVar ,p, &
                            m%SPLP,m%SPLS,m%SPLALPH,m%SPLTBL,errStat2 ,errMsg2)
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
                    ENDIF
                ENDIF
                !--------Blunt Trailing Edge Noise----------------------------------------------!
                IF ( p%IBLUNT .EQ. 1 )   THEN                                          
                    CALL BLUNT(AlphaNoise,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
                    elementspan,m%rTEtoObserve(K,J,I),p%TEThick(J,I),p%TEAngle(J,I), &
                    p, m%d99Var(2),m%dstarVar(1),m%dstarVar(2),m%SPLBLUNT,p%StallStart(J,I),errStat2,errMsg2 )
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
                ENDIF
                !--------Tip Noise--------------------------------------------------------------!
                IF (  (p%ITIP .EQ. 1) .AND. (J .EQ. p%NumBlNds)  ) THEN 
                    CALL TIPNOIS(AlphaNoise,p%ALpRAT,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
                        m%rTEtoObserve(K,J,I), p, m%SPLTIP,errStat2,errMsg2)
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
                ENDIF
                !--------Inflow Turbulence Noise ------------------------------------------------!
                ! important checks to be done inflow tubulence inputs
                IF (p%IInflow.gt.0) then

                    ! Amiet's Inflow Noise Model is Calculated as long as InflowNoise is On
                    CALL InflowNoise(AlphaNoise,p%BlChord(J,I),Unoise,m%ChordAngleLE(K,J,I),m%SpanAngleLE(K,J,I),&
                        elementspan,m%rLEtoObserve(K,J,I),xd%TIVx(J,I),p,m%SPLti,errStat2,errMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
                    ! If Guidati model (simplified or full version) is also on then the 'SPL correction' to Amiet's model will be added 
                    IF ( p%IInflow .EQ. 2 )   THEN      
                        CALL Simple_Guidati(UNoise,p%BlChord(J,I),p%AFThickGuida(2,p%BlAFID(J,I)), &
                            p%AFThickGuida(1,p%BlAFID(J,I)),p,m%SPLTIGui,errStat2,errMsg2 )
                        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
                        m%SPLti=m%SPLti+m%SPLTIGui + 10. ! +10 is fudge factor to match NLR data
                    ELSEIF ( p%IInflow .EQ. 3 )   THEN                                     
                       CALL WrScr('Full Guidati removed')
                       STOP
                    ENDIF    
                ENDIF
                !----------------------------------------------------------------------------------------------------------------------------------!
                !      ADD IN THIS SEGMENT'S CONTRIBUTION ON A MEAN-SQUARE
                !      PRESSURE BASIS
                !----------------------------------------------------------------------------------------------------------------------------------!
                Ptotal = 0.0_ReKi         ! Total Sound Pressure - All (7) mechanisms, All Frequencies
               PtotalLBL= 0.0_ReKi        ! Total Sound Pressure - Laminar Boundary Layer, All Frequencies
               PtotalTBLP= 0.0_ReKi       ! Total Sound Pressure - Turbulent Boundary Layer, Pressure Contribution, All Frequencies
               PtotalTBLS= 0.0_ReKi       ! Total Sound Pressure - Turbulent Boundary Layer, Suction Contribution, All Frequencies
               PtotalSep= 0.0_ReKi        ! Total Sound Pressure - Separation, All Frequencies
               PtotalTBLAll = 0.0_ReKi    ! Total Sound Pressure - Turbulent Boundary Layer, All Frequencies
               PtotalBlunt= 0.0_ReKi      ! Total Sound Pressure - Blunt Trailing Edge, All Frequencies
               PtotalTip= 0.0_ReKi        ! Total Sound Pressure - Tip Noise, All Frequencies
               PtotalInflow= 0.0_ReKi     ! Total Sound Pressure - Turbulent Inflow, All Frequencies
               PLBL= 0.0_ReKi             ! Laminar Boundary Layer - Current Iteration
               PTBLP= 0.0_ReKi            ! Turbulent Boundary Layer, Pressure Contribution - Current Iteration
               PTBLS= 0.0_ReKi            ! Turbulent Boundary Layer, Suction Contribution - Current Iteration
               PTBLALH= 0.0_ReKi          ! Turbulent Boundary Layer, Angle of Attack Contribution - Current Iteration (Feeds into PTotalSep. Consider renaming.)
               PTip= 0.0_ReKi             ! Tip Noise - Current Iteration
               PTI= 0.0_ReKi              ! Turbulent Inflow - Current Iteration
               PBLNT= 0.0_ReKi            ! Blunt Trailing Edge - Current Iteration

            
               DO III=1,size(p%FreqList)   ! Loops through each 1/3rd octave center frequency 
         
                  ! If flag for LBL is ON and Boundary Layer Trip is OFF, then compute LBL
                  IF ( (p%ILAM .EQ. 1) .AND. (p%ITRIP .EQ. 0) )  THEN
                     IF (p%AweightFlag .eqv. .TRUE.) THEN
                         m%SPLLBL(III) = m%SPLLBL(III) + p%Aweight(III)                ! A-weighting
                     ENDIF
                        
                     PLBL = 10.0_ReKi**(m%SPLLBL(III)/10.0_ReKi)                       ! SPL to Sound Pressure (P) Conversion for III Frequency
                        
                     PtotalLBL = PtotalLBL + PLBL                                      ! Sum of Current LBL with LBL Running Total
                        Ptotal = Ptotal + PLBL                                         ! Sum of Current LBL with Overall Running Total
                     y%PtotalFreq(K,III) = y%PtotalFreq(K,III) + PLBL                  ! Running sum of observer and frequency dependent sound pressure
                  
                     y%SumSpecNoiseSep(1,K,III) = PLBL + y%SumSpecNoiseSep(1,K,III)    ! Assigns Current LBL to Appropriate Mechanism (1), Observer (K), and Frequency (III)
                  ENDIF

                  ! If flag for TBL is ON, compute Pressure, Suction, and AoA contributions
                  IF ( p%ITURB .GT. 0 )  THEN
                     IF (p%AweightFlag .eqv. .TRUE.) THEN
                        m%SPLP(III) = m%SPLP(III) + p%Aweight(III)                     ! A-weighting
                        m%SPLS(III) = m%SPLS(III) + p%Aweight(III)                     ! A-weighting
                        m%SPLALPH(III) = m%SPLALPH(III) + p%Aweight(III)               ! A-weighting
                     ENDIF

                     PTBLP = 10.0_ReKi**(m%SPLP(III)/10.0_ReKi)                        ! SPL to P Conversion for III Frequency
                     PTBLS = 10.0_ReKi**(m%SPLS(III)/10.0_ReKi)                        ! SPL to P Conversion for III Frequency
                     PTBLALH = 10.0_ReKi**(m%SPLALPH(III)/10.0_ReKi)                   ! SPL to P Conversion for III Frequency
                        
                     PtotalTBLP = PtotalTBLP + PTBLP                                   ! Sum of Current TBLP with TBLP Running Total
                     PtotalTBLS = PtotalTBLS + PTBLS                                   ! Sum of Current TBLS with TBLS Running Total         
                     PtotalSep  = PtotalSep  + PTBLALH                                 ! Sum of Current TBLALH with TBLALH Running Total
                  
                     Ptotal = Ptotal + PTBLP + PTBLS + PTBLALH                         ! Sum of Current TBL with Overall Running Total
                     y%PtotalFreq(K,III) = y%PtotalFreq(K,III) + PTBLP + PTBLS + PTBLALH  ! Running sum of observer and frequency dependent sound pressure
                     PtotalTBLAll = PtotalTBLAll + 10.0_ReKi**(m%SPLTBL(III)/10.0_ReKi)   ! SPLTBL from comment on line 1794 is the mean-square sum of SPLP, SPLS, and SPLALPH.
                                                                                          !   So this should be equal to PTBLP+PTBLS+TBLALH
                     y%SumSpecNoiseSep(2,K,III) = PTBLP   + y%SumSpecNoiseSep(2,K,III)    ! Assigns Current TBLP to Appropriate Mechanism (2), Observer (K), and Frequency (III)
                     y%SumSpecNoiseSep(3,K,III) = PTBLS   + y%SumSpecNoiseSep(3,K,III)    ! Assigns Current TBLS to Appropriate Mechanism (2), Observer (K), and Frequency (III)
                     y%SumSpecNoiseSep(4,K,III) = PTBLALH + y%SumSpecNoiseSep(4,K,III)    ! Assigns Current TBLALH to Appropriate Mechanism (2), Observer (K), and Frequency (III)
                  ENDIF

                  ! If flag for Blunt TE is ON, compute Blunt contribution
                  IF ( p%IBLUNT .GT. 0 )  THEN                                            ! NOTE: .EQ. 1 would be more accurate since only options are 0 and 1
                     IF (p%AweightFlag .eqv. .TRUE.) THEN
                        m%SPLBLUNT(III) = m%SPLBLUNT(III) + p%Aweight(III)                ! A-weighting
                     ENDIF
                        
                     PBLNT = 10.0_ReKi**(m%SPLBLUNT(III)/10.0_ReKi)                       ! SPL to P Conversion for III Frequency
                        
                     PtotalBlunt = PtotalBlunt + PBLNT                                    ! Sum of Current Blunt with Blunt Running Total
                     Ptotal = Ptotal + PBLNT                                              ! Sum of Current Blunt with Overall Running Total
                     y%PtotalFreq(K,III) = y%PtotalFreq(K,III) + PBLNT                    ! Running sum of observer and frequency dependent sound pressure
                        
                     y%SumSpecNoiseSep(5,K,III) = PBLNT + y%SumSpecNoiseSep(5,K,III)      ! Assigns Current Blunt to Appropriate Mechanism (5), Observer (K), and Frequency (III)
                  ENDIF

                  ! If flag for Tip is ON and the current blade node (J) is the last node (tip), compute Tip contribution
                  IF ( (p%ITIP .GT. 0) .AND. (J .EQ. p%NumBlNds) )  THEN                  ! NOTE: .EQ. 1 would again be more accurate
                     IF (p%AweightFlag .eqv. .TRUE.) THEN
                        m%SPLTIP(III) = m%SPLTIP(III) + p%Aweight(III)                    ! A-weighting
                     ENDIF
                        
                     PTip = 10.0_ReKi**(m%SPLTIP(III)/10.0_ReKi)                          ! SPL to P Conversion for III Frequency
                        
                     PtotalTip = PtotalTip + PTip                                         ! Sum of Current Tip with Tip Running Total
                     Ptotal = Ptotal + PTip                                               ! Sum of Current Tip with Overall Running Total
                     y%PtotalFreq(K,III) = y%PtotalFreq(K,III) + PTip                     ! Running sum of observer and frequency dependent sound pressure
                  
                     y%SumSpecNoiseSep(6,K,III) = PTip + y%SumSpecNoiseSep(6,K,III)       ! Assigns Current Tip to Appropriate Mechanism (6), Observer (K), and Frequency (III)
                  ENDIF

                  ! If flag for TI is ON, compute Turbulent Inflow contribution
                  IF ( (p%IInflow .GT. 0)  )  THEN
                     IF (p%AweightFlag .eqv. .TRUE.) THEN
                        m%SPLti(III) = m%SPLti(III) + p%Aweight(III)                      ! A-weighting
                     ENDIF
                        
                     PTI = 10.0_ReKi**(m%SPLti(III)/10.0_ReKi)                            ! SPL to P Conversion for III Frequency
                        
                     PtotalInflow = PtotalInflow + PTI                                    ! Sum of Current TI with TI Running Total
                     Ptotal = Ptotal + PTI                                                ! Sum of Current TI with Overall Running Total
                     y%PtotalFreq(K,III) = y%PtotalFreq(K,III) + PTI                      ! Running sum of observer and frequency dependent sound pressure
                  
                     y%SumSpecNoiseSep(7,K,III) = PTI + y%SumSpecNoiseSep(7,K,III)        ! Assigns Current TI to Appropriate Mechanism (7), Observer (K), and Frequency (III)
                  ENDIF
                    
                ENDDO ! III = 1, size(p%FreqList)
              
              y%DirectiviOutput(K)         = Ptotal + y%DirectiviOutput(K)            ! Assigns Overall Pressure to Appropriate Observer for Directivity   
              IF (y%DirectiviOutput(K)  .EQ. 0.)      y%DirectiviOutput(K) = 1        ! Since these will all be converted via LOG10, they will produce an error if .EQ. 0. 
                                                                                          !    Set .EQ. to 1 instead (LOG10(1)=0)
              y%OASPL(K,J,I) = Ptotal + y%OASPL(K,J,I)               ! Assigns Overall Pressure to Appropriate Observer/Blade/Node for Directivity
           ENDDO ! Loop on observers
       ENDDO ! Loop on blade nodes
   ENDDO ! Loop on blades

    ! If any Output file is wanted, convert DirectiviOutput from Directivity Factor to Directivity Index
    ! Ref: Fundamentals of Acoustics by Colin Hansen (1951)
    y%DirectiviOutput = 10.*LOG10(y%DirectiviOutput)        !! DirectiviOutput is used as total observer OASPL for Output File 1   
    ! Since these will all be converted via LOG10, they will produce an error if .EQ. 0., Set .EQ. to 1 instead (LOG10(1)=0)
    DO I = 1,p%numBlades
        DO J = 1,p%NumBlNds
            DO K = 1,p%NrObsLoc 
                IF (y%OASPL(K,J,I)  .EQ. 0.)      y%OASPL(K,J,I) = 1   
            ENDDO
        ENDDO
    ENDDO
    IF  (p%NrOutFile .gt. 0) y%OASPL = 10.*LOG10(y%OASPL)                            !! OASPL is used as observer/blade/node OASPL for Output File 4

    ! Procedure for Output file 2
    IF  (p%NrOutFile .gt. 1) THEN 
      DO K = 1,p%NrObsLoc
         DO III=1,size(p%FreqList)
            IF (y%PtotalFreq(K,III) .EQ. 0.)       y%PtotalFreq(K,III) = 1
                y%PtotalFreq(K,III)    = 10.*LOG10(y%PtotalFreq(K,III))               ! P to SPL conversion
         ENDDO
      ENDDO
    ENDIF

   ! If 3rd Output file is needed, these will need to be converted via LOG10. Change to equal 1 to avoid error.
   DO K = 1,p%NrObsLoc
      DO III = 1,size(p%FreqList)
         DO oi = 1,7
            IF (y%SumSpecNoiseSep(oi,K,III)  .EQ. 0.) y%SumSpecNoiseSep(oi,K,III) = 1 
         ENDDO
      ENDDO
   ENDDO
   
   ! Procedure for Output file 3
   IF  (p%NrOutFile .gt. 2) THEN 
        y%SumSpecNoiseSep = 10.*LOG10(y%SumSpecNoiseSep)      ! P to SPL Conversion
    ENDIF
   
END SUBROUTINE CalcAeroAcousticsOutput
!==================================================================================================================================!
SUBROUTINE LBLVS(ALPSTAR,C,U,THETA,PHI,L,R,p,d99Var2,dstarVar1,dstarVar2,SPLLAM,StallVal,errStat,errMsg)
    REAL(ReKi),                                 INTENT(IN   ) :: ALPSTAR        ! AOA
    REAL(ReKi),                                 INTENT(IN   ) :: C              ! Chord Length
    REAL(ReKi),                                 INTENT(IN   ) :: U              ! Unoise FREESTREAM VELOCITY                METERS/SEC
    REAL(ReKi),                                 INTENT(IN   ) :: THETA          ! DIRECTIVITY ANGLE                  DEGREES
    REAL(ReKi),                                 INTENT(IN   ) :: PHI            ! DIRECTIVITY ANGLE                  DEGREES
    REAL(ReKi),                                 INTENT(IN   ) :: L              ! SPAN                               METERS
    REAL(ReKi),                                 INTENT(IN   ) :: R              !  OBSERVER DISTANCE FROM SEGMENT     METERS
    REAL(ReKi),                                 INTENT(IN   ) :: d99Var2        !
    REAL(ReKi),                                 INTENT(IN   ) :: dstarVar1              !
    REAL(ReKi),                                 INTENT(IN   ) :: dstarVar2              !
    REAL(ReKi),                                 INTENT(IN   ) :: StallVal               !  
    TYPE(AA_ParameterType),                     INTENT(IN   ) :: p          ! Noise module Parameters
    REAL(ReKi),DIMENSION(size(p%FreqList)),     INTENT(  OUT) :: SPLLAM         !
    INTEGER(IntKi),                             INTENT(  OUT) :: errStat        ! Error status of the operation
    character(*),                               INTENT(  OUT) :: errMsg         ! Error message if ErrStat /= ErrID_None
    integer(intKi)                                             :: ErrStat2           ! temporary Error status
    character(ErrMsgLen)                                       :: ErrMsg2            ! temporary Error message
    character(*), parameter                                    :: RoutineName = 'LBLVS'
    ! Local variables
    real(ReKi)     :: STPRIM   !  STROUHAL NUMBER BASED ON PRESSURE SIDE BOUNDARY LAYER THICKNESS    ---
    real(ReKi)     :: M        ! MACH NUMBER
    real(ReKi)     :: RC       ! REYNOLDS NUMBER BASED ON  CHORD
    real(ReKi)     :: DELTAP   ! PRESSURE SIDE BOUNDARY LAYER THICKNESS METERS
    real(ReKi)     :: DSTRS    ! SUCTION SIDE BOUNDARY LAYER DISPLACEMENT THICKNESS           METERS
    real(ReKi)     :: DSTRP    ! PRESSURE SIDE BOUNDARY LAYER DISPLACEMENT THICKNESS           METERS
    real(ReKi)     :: DBARH    ! HIGH FREQUENCY DIRECTIVITY             ---
    real(ReKi)     :: ST1PRIM  ! REFERENCE STROUHAL NUMBER          ---
    real(ReKi)     :: STPKPRM  ! PEAK STROUHAL NUMBER               ---
    real(ReKi)     :: RC0      ! REFERENCE REYNOLDS NUMBER          ---
    real(ReKi)     :: D        ! REYNOLDS NUMBER RATIO              ---
    real(ReKi)     :: G1       ! SOUND PRESSURE LEVEL FUNCTION      DB
    real(ReKi)     :: G2       ! OVERALL SOUND PRESSURE LEVEL FUNCTION    DB
    real(ReKi)     :: G3       ! OVERALL SOUND PRESSURE LEVEL FUNCTION    DB
    real(ReKi)     :: E        ! STROUHAL NUMBER RATIO              ---
    real(ReKi)     :: SCALE    ! GEOMETRIC SCALING TERM
    integer(intKi) :: I        ! I A generic index for DO loops.
    ErrStat = ErrID_None
    ErrMsg  = ""
    !compute reynolds number and mach number
    M          = U  / p%SpdSound        ! MACH NUMBER
    RC         = U  * C/p%KinVisc       ! REYNOLDS NUMBER BASED ON  CHORD
    ! compute boundary layer thicknesses
    IF (p%X_BLMethod .eq. X_BLMethod_Tables) THEN
        DELTAP = d99Var2
        DSTRS  = dstarVar1
        DSTRP  = dstarVar2
    ELSE
        CALL THICK(C,RC,ALPSTAR,p,DELTAP,DSTRS,DSTRP,StallVal,errStat2,errMsg2)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
    ENDIF
    ! compute directivity function
    CALL DIRECTH_TE(M,THETA,PHI,DBARH,errStat2,errMsg2)
    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 

    IF (DBARH <= 0) THEN
        SPLLAM = 0.
        RETURN
    ENDIF
    ! compute reference strouhal number                                 ! Eq 55 from BPM Airfoil Self-noise and Prediction paper
    IF (RC .LE. 1.3E+05) ST1PRIM = .18
    IF((RC .GT. 1.3E+05).AND.(RC.LE.4.0E+05))ST1PRIM=.001756*RC**.3931
    IF (RC .GT. 4.0E+05) ST1PRIM = .28
    STPKPRM  = 10.**(-.04*ALPSTAR) * ST1PRIM                            ! Eq 56 from BPM Airfoil Self-noise and Prediction paper

    ! compute reference reynolds number                                 ! Eq 59 from BPM Airfoil Self-noise and Prediction paper
    IF (ALPSTAR .LE. 3.0) RC0=10.**(.215*ALPSTAR+4.978)
    IF (ALPSTAR .GT. 3.0) RC0=10.**(.120*ALPSTAR+5.263)
    ! compute peak scaled spectrum level
    D   = RC / RC0                                                      ! Used in Eq 58 from BPM Airfoil Self-noise and Prediction paper
    IF (D .LE. .3237)                        G2 =77.852*LOG10(D)+15.328       ! Begin Eq 58 from BPM Airfoil Self-noise and Prediction paper
    IF ((D .GT. .3237).AND.(D .LE. .5689))   G2 = 65.188*LOG10(D) + 9.125
    IF ((D .GT. .5689).AND.(D .LE. 1.7579))  G2 = -114.052 * LOG10(D)**2
    IF ((D .GT. 1.7579).AND.(D .LE. 3.0889)) G2 = -65.188*LOG10(D)+9.125
    IF (D .GT. 3.0889)                       G2 =-77.852*LOG10(D)+15.328      ! end
    ! compute angle-dependent level for shape curve
   G3      = 171.04 - 3.03 * ALPSTAR                                    ! Eq 60 from BPM Airfoil Self-noise and Prediction paper
    SCALE   = 10. * LOG10(DELTAP*M**5*DBARH*L/R**2)                     ! From Eq 53 from BPM Airfoil Self-noise and Prediction paper
    ! Compute scaled sound pressure levels for each strouhal number
    DO I=1,SIZE(p%FreqList)
        STPRIM  = p%FreqList(I) * DELTAP / U                            ! Eq 54 from BPM Airfoil Self-noise and Prediction paper
        E          = STPRIM / STPKPRM                                   ! Used in Eq 57 from BPM Airfoil Self-noise and Prediction paper
        IF (E .LE. .5974)                      G1 = 39.8*LOG10(E)-11.12                   ! Begin Eq 57 from BPM Airfoil Self-noise and Prediction paper   
        IF ((E .GT. .5974).AND.(E .LE. .8545)) G1 = 98.409 * LOG10(E) + 2.0
        IF ((E .GT. .8545).AND.(E .LE. 1.17))  G1 = -5.076+SQRT(2.484-506.25*(LOG10(E))**2)
        IF ((E .GT. 1.17).AND.(E .LE. 1.674))  G1 = -98.409 * LOG10(E) + 2.0
        IF (E .GT. 1.674)                      G1 = -39.80*LOG10(E)-11.12                 ! end
        SPLLAM(I) = G1 + G2 + G3 + SCALE                                      ! Eq 53 from BPM Airfoil Self-noise and Prediction paper
    ENDDO
END SUBROUTINE LBLVS
!==================================================================================================================================!
SUBROUTINE TBLTE(ALPSTAR,C,U,THETA,PHI,L,R,p,d99Var2,dstarVar1,dstarVar2,StallVal,SPLP,SPLS,SPLALPH,SPLTBL,errStat,errMsg)
    REAL(ReKi),                             INTENT(IN   )  :: ALPSTAR        ! AOA(deg)
    REAL(ReKi),                             INTENT(IN   )  :: C              ! Chord Length           (m)
!    REAL(ReKi),                             INTENT(IN   )  :: U              ! Unoise(m/s)
!    REAL(ReKi),                             INTENT(IN   )  :: THETA          ! DIRECTIVITY ANGLE      (deg)
!    REAL(ReKi),                             INTENT(IN   )  :: PHI            ! DIRECTIVITY ANGLE      (deg) 
    REAL(ReKi),                             INTENT(IN   )  :: L              ! SPAN(m)
    REAL(ReKi),                             INTENT(IN   )  :: R              ! SOURCE TO OBSERVER DISTANCE (m)

!    REAL(ReKi)                             :: ALPSTAR        ! AOA(deg)
!    REAL(ReKi)                               :: C              ! Chord Length           (m)
    REAL(ReKi)                               :: U              ! Unoise(m/s)
    REAL(ReKi)                               :: THETA          ! DIRECTIVITY ANGLE      (deg)
    REAL(ReKi)                               :: PHI            ! DIRECTIVITY ANGLE      (deg) 
!    REAL(ReKi)                               :: L              ! SPAN(m)
!    REAL(ReKi)                               :: R              ! SOURCE TO OBSERVER DISTANCE (m)

    REAL(ReKi),                             INTENT(IN   )  :: d99Var2        !  
    REAL(ReKi),                             INTENT(IN   )  :: dstarVar1              !  
    REAL(ReKi),                             INTENT(IN   )  :: dstarVar2              !  
    REAL(ReKi),                             INTENT(IN   )  :: StallVal              !  

    TYPE(AA_ParameterType),                         INTENT(IN   )  :: p              ! Noise Module Parameters
    REAL(ReKi),DIMENSION(size(p%FreqList)),  INTENT(  OUT)  :: SPLP           ! SOUND PRESSURE LEVEL DUE TO PRESSURE SIDE OF AIRFOIL (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),  INTENT(  OUT)  :: SPLS           ! SOUND PRESSURE LEVEL DUE TO SUCTION SIDE OF AIRFOIL  (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),  INTENT(  OUT)  :: SPLTBL         ! TOTAL SOUND PRESSURE LEVEL DUE TO TBLTE MECHANISM    (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),  INTENT(  OUT)  :: SPLALPH        ! SOUND PRESSURE LEVEL DUE TO ANGLE OF ATTACK CONTRIBUTION (db)
    INTEGER(IntKi),                     INTENT(  OUT)  :: errStat        ! Error status of the operation
    character(*),                                   INTENT(  OUT)  :: errMsg         ! Error message if ErrStat /= ErrID_None
    integer(intKi)                                                 :: ErrStat2       ! temporary Error status
    character(ErrMsgLen)                                           :: ErrMsg2        ! temporary Error message
    character(*), parameter                                        :: RoutineName = 'TBLTE'
    ! Local variables
    real(ReKi)   :: STP        ! PRESSURE SIDE STROUHAL NUMBER          --- 
    real(ReKi)   :: STS        ! SUCTION SIDE STROUHAL NUMBER           ---
    real(ReKi)   :: DSTRS      ! SUCTION SIDE DISPLACEMENT THICKNESS   METERS
    real(ReKi)   :: DSTRP      ! PRESSURE SIDE DISPLACEMENT THICKNESS  METERS
    real(ReKi)   :: RDSTRS     ! REYNOLDS NUMBER BASED ON SUCTION  SIDE DISPLACEMENT THICKNESS 
    real(ReKi)   :: RDSTRP     ! REYNOLDS NUMBER BASED ON PRESSURE SIDE DISPLACEMENT THICKNESS
    real(ReKi)   :: ST1        ! PEAK STROUHAL NUMBER                   ---
    real(ReKi)   :: ST2        ! PEAK STROUHAL NUMBER                   ---
    real(ReKi)   :: ST1PRIM    ! PEAK STROUHAL NUMBER                   ---
    real(ReKi)   :: A0         ! FUNCTION USED IN 'A' CALCULATION
    real(ReKi)   :: A02        ! FUNCTION USED IN 'A' CALCULATION
    real(ReKi)   :: ARA0       ! INTERPOLATION FACTOR
    real(ReKi)   :: ARA02      ! INTERPOLATION FACTOR
    real(ReKi)   :: B0         ! FUNCTION USED IN 'B' CALCULATION
    real(ReKi)   :: BMINB0     ! MINIMUM 'B' EVALUATED AT B0            DB
    real(ReKi)   :: BMINB      ! MINIMUM 'B' EVALUATED AT B             DB
    real(ReKi)   :: BMAXB0     ! MAXIMUM 'B' EVALUATED AT B0            DB
    real(ReKi)   :: BMAXB      ! MAXIMUM 'B' EVALUATED AT B             DB
    real(ReKi)   :: BRB0       ! INTERPOLATION FACTOR                   DB
    real(ReKi)   :: STPEAK     ! PEAK STROUHAL NUMBER                   ---
    real(ReKi)   :: AMINA      ! MINIMUM 'A' CURVE EVALUATED AT STROUHAL NUMBER RATIO                DB
    real(ReKi)   :: AMINB      ! MINIMUM 'A' CURVE EVALUATED AT B       DB
    real(ReKi)   :: AMAXA      ! MAXIMUM 'A' CURVE EVALUATED AT STROUHAL NUMBER RATIO    (DB)
    real(ReKi)   :: AMAXB      ! MAXIMUM 'A' CURVE EVALUATED AT B       DB
    real(ReKi)   :: AMINA0     ! MAXIMUM 'B' EVALUATED AT B0            DB
    real(ReKi)   :: AMINA02    ! MINIMUM 'A' CURVE EVALUATED AT A02     DB
    real(ReKi)   :: AMAXA0     ! MAXIMUM 'A' CURVE EVALUATED AT A0      DB
    real(ReKi)   :: AMAXA02    ! MAXIMUM 'A' CURVE EVALUATED AT A02      DB
    real(ReKi)   :: A          ! STROUHAL NUMBER RATIO                 ---
    real(ReKi)   :: B          ! STROUHAL NUMBER RATIO                 ---
    real(ReKi)   :: AA         ! 'A' SPECTRUM SHAPE EVALUATED AT STROUHAL NUMBER RATIO         DB
    real(ReKi)   :: BB         ! 'B' SPECTRUM SHAPE EVALUATED AT STROUHAL NUMBER RATIO                DB
    real(ReKi)   :: DELK1      ! CORRECTION TO AMPLITUDE FUNCTION       DB
    real(ReKi)   :: GAMMA      ! USED IN 'B' COMPUTATION                ---
    real(ReKi)   :: BETA       ! USED IN 'B' COMPUTATION               ---
    real(ReKi)   :: GAMMA0     ! USED IN 'B' COMPUTATION                ---
    real(ReKi)   :: BETA0      ! USED IN 'B' COMPUTATION               ---
    real(ReKi)   :: K1         ! AMPLITUDE FUNCTION (DB)
    real(ReKi)   :: K2         ! AMPLITUDE FUNCTION (DB)
    real(ReKi)   :: P1         ! PRESSURE SIDE PRESSURE (NT/M2)
    real(ReKi)   :: P2         ! SUCTION SIDE PRESSURE                      (NT/M2)
    real(ReKi)   :: P4         ! PRESSURE FROM ANGLE OF ATTACK CONTRIBUTION (NT/M2)
    real(ReKi)   :: M          ! MACH NUMBER
    real(ReKi)   :: RC         ! REYNOLDS NUMBER BASED ON  CHORD
    real(ReKi)   :: DELTAP     ! PRESSURE SIDE BOUNDARY LAYER THICKNESS METERS
    real(ReKi)   :: XCHECK     ! USED TO CHECK FOR ANGLE OF ATTACK CONTRIBUTION     
    real(ReKi)   :: DBARH      ! HIGH FREQUENCY DIRECTIVITY             ---
    real(ReKi)   :: DBARL      ! LOW FREQUENCY DIRECTIVITY              ---

    integer(intKi)      :: I          ! I A generic index for DO loops.

    LOGICAL     :: SWITCH  !!LOGICAL FOR COMPUTATION OF ANGLE OF ATTACK CONTRIBUTION  



    ErrStat = ErrID_None
    ErrMsg  = ""
    ! Compute reynolds number and mach number
    M          = U  / p%SpdSound
    RC         = U  * C/p%KinVisc
    ! Compute boundary layer thicknesses
    IF (p%X_BLMethod .eq. X_BLMethod_Tables) THEN
        DELTAP = d99Var2
        DSTRS  = dstarVar1
        DSTRP  = dstarVar2
    ELSE
        CALL THICK(C,RC,ALPSTAR,p,DELTAP,DSTRS,DSTRP,StallVal,errStat2,errMsg2)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
    ENDIF
    ! Compute directivity function
    CALL DIRECTL(M,THETA,PHI,DBARL,errStat2,errMsg2)
    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
    CALL DIRECTH_TE(M,THETA,PHI,DBARH,errStat2,errMsg2)
    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
    !      IF (DBARH <= 0) THEN
    !          SPLP = 0.
    !          SPLS = 0.
    !          SPLALPH = 0.
    !          RETURN
    !      ENDIF
    ! Calculate the reynolds numbers based on pressure and suction displacement thickness
    RDSTRS = DSTRS * U  / p%KinVisc
    RDSTRP = DSTRP * U  / p%KinVisc
    ! Determine peak strouhal numbers to be used for 'a' and 'b' curve calculations
    ST1    = .02 * M ** (-.6)                                                          ! Eq 32 from BPM Airfoil Self-noise and Prediction paper
     ! Eq 34 from BPM Airfoil Self-noise and Prediction paper
    IF  (ALPSTAR .LE. 1.333)                          ST2 = ST1
    IF ((ALPSTAR .GT. 1.333).AND.(ALPSTAR .LE. StallVal)) ST2 = ST1*10.**(.0054*(ALPSTAR-1.333)**2)
    IF (ALPSTAR .GT. StallVal)                           ST2 = 4.72 * ST1
    ST1PRIM = (ST1+ST2)/2.                                                             ! Eq 33 from BPM Airfoil Self-noise and Prediction paper
    CALL A0COMP(RC,A0)      ! compute -20 dB dropout   (returns A0)
    CALL A0COMP(3.0_ReKi*RC,A02)   ! compute -20 dB dropout for AoA > AoA_0   (returns A02)
    ! Evaluate minimum and maximum 'a' curves at a0
    CALL AMIN(A0,AMINA0)
    CALL AMAX(A0,AMAXA0)
    CALL AMIN(A02,AMINA02)
    CALL AMAX(A02,AMAXA02)
    ! Compute 'a' max/min ratio                                                        ! Eq 39 from BPM Airfoil Self-noise and Prediction paper   
    ARA0  = (20. + AMINA0) / (AMINA0 - AMAXA0)
    ARA02 = (20. + AMINA02)/ (AMINA02- AMAXA02)
    ! Compute b0 to be used in 'b' curve calculations                                  ! Eq 44 from BPM Airfoil Self-noise and Prediction paper
    IF (RC .LT. 9.52E+04) B0 = .30
    IF ((RC .GE. 9.52E+04).AND.(RC .LT. 8.57E+05)) &
        B0 = (-4.48E-13)*(RC-8.57E+05)**2 + .56
    IF (RC .GE. 8.57E+05) B0 = .56
    ! Evaluate minimum and maximum 'b' curves at b0
    CALL BMIN(B0,BMINB0)
    CALL BMAX(B0,BMAXB0)
    ! Compute 'b' max/min ratio
    BRB0  = (20. + BMINB0) / (BMINB0 - BMAXB0)

    ! For each center frequency, compute an 'a' prediction for the pressure side
    STPEAK = ST1
    IF (RC .LT. 2.47E+05)                        K1 = -4.31 * LOG10(RC) + 156.3        ! Begin Eq 47 from BPM Airfoil Self-noise and Prediction paper         
    IF((RC .GE. 2.47E+05).AND.(RC .LE. 8.0E+05)) K1 = -9.0 * LOG10(RC) + 181.6
    IF (RC .GT. 8.0E+05)                         K1 = 128.5                            ! end
   IF (RDSTRP .LE. 5000.) DELK1 = -ALPSTAR*(5.29-1.43*LOG10(RDSTRP))                   ! Begin Eq 48 from BPM Airfoil Self-noise and Prediction paper
    IF (RDSTRP .GT. 5000.) DELK1 = 0.0                                                 ! end      

    GAMMA   = 27.094 * M +  3.31                                                       ! Begin Eq 49 from BPM Airfoil Self-noise and Prediction paper
    BETA    = 72.650 * M + 10.74
    GAMMA0  = 23.430 * M +  4.651
    BETA0   =-34.190 * M - 13.820                                                      ! end

    IF (ALPSTAR .LE. (GAMMA0-GAMMA)) K2 = -1000.0                                      ! Begin Eq 49 from BPM Airfoil Self-noise and Prediction paper
    IF ((ALPSTAR.GT.(GAMMA0-GAMMA)).AND.(ALPSTAR.LE.(GAMMA0+GAMMA))) &
        K2=SQRT(BETA**2-(BETA/GAMMA)**2*(ALPSTAR-GAMMA0)**2)+BETA0
    IF (ALPSTAR .GT. (GAMMA0+GAMMA)) K2 = -12.0
    K2 = K2 + K1                                                                       ! end
    ! Check for 'a' computation for suction side
    XCHECK = GAMMA0
    SWITCH = .FALSE.
    !older version:
    !      IF ((ALPSTAR .GE. XCHECK).OR.(ALPSTAR .GT. 12.5))SWITCH=.TRUE. 
    ! newer version
    IF ((ALPSTAR .GE. XCHECK).OR.(ALPSTAR .GT. StallVal))SWITCH=.TRUE. 
    DO  I=1,size(p%FreqList)
        STP= p%FreqList(I) * DSTRP / U                                     ! Eq 31 from BPM Airfoil Self-noise and Prediction paper   
        A      = LOG10( STP / STPEAK )                                     ! Eq 37 from BPM Airfoil Self-noise and Prediction paper
        CALL AMIN(A,AMINA)
        CALL AMAX(A,AMAXA)
        AA     = AMINA + ARA0 * (AMAXA - AMINA)                            ! Eq 40 from BPM Airfoil Self-noise and Prediction paper

        SPLP(I)=AA+K1-3.+10.*LOG10(DSTRP*M**5*DBARH*L/R**2)+DELK1        ! Eq 25 from BPM Airfoil Self-noise and Prediction paper
        STS = p%FreqList(I) * DSTRS / U                                    ! Eq 31 from BPM Airfoil Self-noise and Prediction paper

        IF (.NOT. SWITCH) THEN
            A      = LOG10( STS / ST1PRIM )
            CALL AMIN(A,AMINA)
            CALL AMAX(A,AMAXA)
            AA = AMINA + ARA0 * (AMAXA - AMINA)
            SPLS(I) = AA+K1-3.+10.*LOG10(DSTRS*M**5*DBARH* L/R**2)       ! Eq 26 from BPM Airfoil Self-noise and Prediction paper
            !  'B' CURVE COMPUTATION
            !        B = ABS(LOG10(STS / ST2))
            B = LOG10(STS / ST2) ! abs not needed absolute taken in the AMAX,AMIN   ! Eq 43 from BPM Airfoil Self-noise and Prediction paper
            CALL BMIN(B,BMINB)
            CALL BMAX(B,BMAXB)
            BB = BMINB + BRB0 * (BMAXB-BMINB)                              ! Eq 46 from BPM Airfoil Self-noise and Prediction paper
            SPLALPH(I)=BB+K2+10.*LOG10(DSTRS*M**5*DBARH*L/R**2)          ! Eq 27 from BPM Airfoil Self-noise and Prediction paper
        ELSE
            ! The 'a' computation is dropped if 'switch' is true
            SPLS(I) = 10.*LOG10(DSTRS*M**5*DBARL*L/R**2)
            !    SPLP(I) = 0.0 + 10.*LOG10(DSTRS*M**5*DBARL*L/R**2) ! changed the line below because the SPLP should be calculatd with DSTRP not with DSTRS
            SPLP(I) = 10.*LOG10(DSTRP*M**5*DBARL*L/R**2) ! this is correct
            !        B = ABS(LOG10(STS / ST2))
            B = LOG10(STS / ST2) ! abs not needed absolute taken in the AMAX,AMIN
            CALL AMIN(B,AMINB)
            CALL AMAX(B,AMAXB)
            BB = AMINB + ARA02 * (AMAXB-AMINB)
            SPLALPH(I)=BB+K2+10.*LOG10(DSTRS*M**5*DBARL*L/R**2)            
        ENDIF
        ! Sum all contributions from 'a' and 'b' on both pressure and suction side on a mean-square pressure basis
        IF (SPLP(I)    .LT. -100.) SPLP(I)    = -100.                      ! Similar to Eq 28 of BPM Airfoil Self-noise and Prediction paper
        IF (SPLS(I)    .LT. -100.) SPLS(I)    = -100.                      ! Similar to Eq 29 of BPM Airfoil Self-noise and Prediction paper      
        IF (SPLALPH(I) .LT. -100.) SPLALPH(I) = -100.                      ! Eq 30 of BPM Airfoil Self-noise and Prediction paper recommends SPLALPH = 10log(stuff) + A' + K2, where A' is calculated same as A but with x3 Rc   

        P1  = 10.**(SPLP(I) / 10.)            ! SPL_Pressure
        P2  = 10.**(SPLS(I) / 10.)            ! SPL_Suction
        P4  = 10.**(SPLALPH(I) / 10.)         ! SPL_AoA   
        SPLTBL(I) = 10. * LOG10(P1 + P2 + P4)                              ! Eq 24 from BPM Airfoil Self-noise and Prediction paper



   ENDDO

END SUBROUTINE TBLTE
!==================================================================================================================================!
SUBROUTINE TIPNOIS(ALPHTIP,ALPRAT2,C,U ,THETA,PHI, R,p,SPLTIP, errStat, errMsg)
    REAL(ReKi),                               INTENT(IN   )  :: ALPHTIP        !< AOA
    REAL(ReKi),                               INTENT(IN   )  :: ALPRAT2        !< TIP LIFT CURVE SLOPE                 ---
    REAL(ReKi),                               INTENT(IN   )  :: C              !< Chord Length
    REAL(ReKi),                               INTENT(IN   )  :: U              !< FREESTREAM VELOCITY               METERS/SEC
    REAL(ReKi),                               INTENT(IN   )  :: THETA          !< DIRECTIVITY ANGLE                  DEGREES 
    REAL(ReKi),                               INTENT(IN   )  :: PHI            !< DIRECTIVITY ANGLE                  DEGREES
    REAL(ReKi),                               INTENT(IN   )  :: R              !< SOURCE TO OBSERVER DISTANCE        METERS
    TYPE(AA_ParameterType) ,                  INTENT(IN   )  :: p              !< Parameters
    REAL(ReKi),DIMENSION(size(p%FreqList)),   INTENT(  OUT)  :: SPLTIP         !<
    INTEGER(IntKi),                           INTENT(  OUT)  :: errStat        !< Error status of the operation
    character(*),                             INTENT(  OUT)  :: errMsg         !< Error message if ErrStat /= ErrID_None
    ! local variables
    integer(intKi)                                             :: ErrStat2       ! temporary Error status
    character(ErrMsgLen)                                       :: ErrMsg2        ! temporary Error message
    character(*), parameter                                    :: RoutineName = 'tipnoise'
    REAL(ReKi)        :: M        ! MACH NUMBER                         ---
    REAL(ReKi)        :: MM       ! MAXIMUM MACH NUMBER                 ---
    REAL(ReKi)        :: ALPTIPP  ! CORRECTED TIP ANGLE OF ATTACK      DEGREES
    REAL(ReKi)        :: DBARH    ! DIRECTIVITY                         ---
    REAL(ReKi)        :: SCALE    ! SCALING TERM                        ---
    REAL(ReKi)        :: STPP     ! STROUHAL NUMBER                     ---
    REAL(ReKi)        :: UM       ! MAXIMUM VELOCITY                  METERS/SEC
    REAL(ReKi)        :: L        ! CHARACTERISTIC LENGTH FOR TIP      METERS
    REAL(ReKi)        :: TERM     ! SCALING TERM                        ---
    integer(intKi)    :: I        !I A generic index for DO loops.
    ErrStat = ErrID_None
    ErrMsg  = ""
    IF (alphtip.eq.0.) THEN
        SPLTIP= 0
        RETURN
    ELSEIF (alphtip.lt.0.) THEN
        !         alphtip = ABS (alphtip) !  (EB_DTU) NOT possible to change inten(in) variable, INSTEAD 
        !  ALPTIPP is equal to abs(alphtip) - see next equation 
    ENDIF
    !! used to be  ALPTIPP = ALPHTIP * ALPRAT2
    ALPTIPP = ABS(ALPHTIP) * ALPRAT2
    M          = U  / p%SpdSound ! MACH NUMBER
    ! Compute directivity function
    CALL DIRECTH_TE(M,THETA,PHI,DBARH,errStat2,errMsg2)
    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
    IF (p%ROUND) THEN
        L = .008 * ALPTIPP * C                                    ! Eq 63 from BPM Airfoil Self-noise and Prediction paper   
    ELSE
        IF (ABS(ALPTIPP) .LE. 2.) THEN                            ! not sure where this comes from   
            L = (.023 + .0169*ALPTIPP) * C
        ELSE
            L = (.0378 + .0095*ALPTIPP) * C
        ENDIF
    ENDIF
    MM    = (1. + .036*ALPTIPP) * M                               ! Eq 64 from BPM Airfoil Self-noise and Prediction paper
    UM    = MM * p%SpdSound                                       ! Eq 65 from BPM Airfoil Self-noise and Prediction paper   
    TERM  = M*M*MM**3*L**2*DBARH/R**2                          ! TERM = M^2 * M_max^5 *l^2 *D / r^2 according to Semi-Empirical Aeroacoustic Noise Prediction Code for Wind Turbines paper
                                                                  ! Term is correct according to Eq 61 from BPM Airfoil self-noise and Prediction paper
    IF (TERM .NE. 0.0) THEN                                       
        SCALE = 10.*LOG10(TERM)
    ELSE
        SCALE = 0.0
    ENDIF
    DO I=1,size(p%FreqList)
        STPP      = p%FreqList(I) * L / UM                        ! Eq 62 from BPM Airfoil Self-noise and Prediction paper   
        SPLTIP(I) = 126.-30.5*(LOG10(STPP)+.3)**2 + SCALE        ! Eq 61 from BPM Airfoil Self-noise and Prediction paper
    ENDDO
END SUBROUTINE TipNois
!==================================================================================================================================!
SUBROUTINE InflowNoise(AlphaNoise,Chord,U,THETA,PHI,d,RObs,TINoise,p,SPLti,errStat,errMsg)
  REAL(ReKi),                                 INTENT(IN   ) :: AlphaNoise     ! AOA
  REAL(ReKi),                                 INTENT(IN   ) :: Chord          ! Chord Length
  REAL(ReKi),                                 INTENT(IN   ) :: U              !
  REAL(ReKi),                                 INTENT(IN   ) :: THETA          !
  REAL(ReKi),                                 INTENT(IN   ) :: PHI            ! Spanwise directivity angle
  REAL(ReKi),                                 INTENT(IN   ) :: d              ! element span
  REAL(ReKi),                                 INTENT(IN   ) :: RObs           ! distance to observer
!  REAL(ReKi),                                 INTENT(IN   ) :: MeanVNoise     !
  REAL(ReKi),                                 INTENT(IN   ) :: TINoise        !
!  REAL(ReKi),                                 INTENT(IN   ) :: LE_Location    !

!  REAL(ReKi),                                 INTENT(IN   ) :: dissip         !
  TYPE(AA_ParameterType),                     INTENT(IN   ) :: p              ! Parameters
  REAL(ReKi),DIMENSION(size(p%FreqList)),     INTENT(  OUT) :: SPLti          !
  INTEGER(IntKi),                             INTENT(  OUT) :: errStat        ! Error status of the operation
  character(*),                               INTENT(  OUT) :: errMsg         ! Error message if ErrStat /= ErrID_None
  integer(intKi)                                            :: ErrStat2       ! temporary Error status
  character(ErrMsgLen)                                      :: ErrMsg2        ! temporary Error message
  character(*), parameter                                   :: RoutineName = 'InflowNoise'
! local variables
  REAL(ReKi)                   :: Beta2                                           ! Prandtl-Glauert correction factor
  REAL(ReKi)                   :: DBARH                                           ! High-frequency directivity correction factor
  REAL(ReKi)                   :: DBARL                                           ! Low-frequency directivity correction factor
  REAL(ReKi)                   :: Directivity                                     ! Directivity correction factor
  REAL(ReKi)                   :: Frequency_cutoff                                ! Cutoff frequency between
  REAL(ReKi)                   :: LFC                                             ! low-frequency correction factor
  REAL(ReKi)                   :: Mach                                            ! local mach number
  REAL(ReKi)                   :: Sears                                           ! Sears function
  REAL(ReKi)                   :: SPLhigh                                         ! predicted high frequency sound pressure level
!  REAL(ReKi)                   :: Ums                                             ! mean square turbulence level
  REAL(ReKi)                   :: WaveNumber                                      ! wave number - non-dimensional frequency
  REAL(ReKi)                   :: Kbar                                      ! nafnoise 
  REAL(ReKi)                   :: khat                                      ! nafnoise 
!  REAL(ReKi)                   :: Kh                                        ! nafnoise 
  REAL(ReKi)                   :: ke                                        ! nafnoise 
  REAL(ReKi)                   :: alpstar                                   ! nafnoise 
!  REAL(ReKi)                   :: mu                                        ! nafnoise 
  REAL(ReKi)                   :: tinooisess                                ! nafnoise 
  ! REAL(ReKi)                   :: L_Gammas                                  ! nafnoise 

  INTEGER(intKi)           :: I        !I A generic index for DO loops.
   ErrStat = ErrID_None
   ErrMsg  = ""

   !!!--- NAF NOISE IDENTICAL
   Mach = U/p%SpdSound
   
   ! This part is recently added for height and surface roughness dependent estimation of turbulence intensity and turbulence scales
   !%Lturb=300*(Z/300)^(0.46+0.074*log(p%z0_aa));              !% Gives larger  length scale
   ! Lturb=25.d0*LE_Location**(0.35)*p%z0_aa**(-0.063)               !% Gives smaller length scale        ! Wei Jun Zhu, Modeling of Aerodynamically generated Noise From Wind Turbines
   ! L_Gammas=0.24+0.096*log10(p%z0_aa)+0.016*(log10(p%z0_aa))**2;   !% Can be computed or just give it a value.    ! Wei Jun Zhu, Modeling of Aerodynamically generated Noise From Wind Turbines
   !tinooisess=L_Gammas*log(30.d0/p%z0_aa)/log(LE_Location/p%z0_aa) !% F.E. 16% is 0.16 which is the correct input for SPLhIgh, no need to divide 100   ! ! Wei Jun Zhu, Modeling of Aerodynamically generated Noise From Wind Turbines
   tinooisess=TINoise

   !tinooisess=0.1
   !Ums = (tinooisess*U)**2
   !Ums = (tinooisess*8)**2
   CALL DIRECTL(Mach,THETA,PHI,DBARL,errStat2,errMsg2) ! assume that noise is low-freq in nature because turbulence length scale is large
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
   CALL DIRECTH_LE(Mach,THETA,PHI,DBARH,errStat2,errMsg2) ! Directivity for the leading edge at high frequencies
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   IF (DBARH <= 0) THEN
      SPLti = 0.
      RETURN
   ENDIF
    
   ! In the following lines, bibliography will be referenced as:  a) Moriarty, Guidati, Migliore, Recent Improvement of a Semi-Empirical Aeroacoustic 
   ! Prediction Code for Wind Turbines
   ! ref b) Lowson, Assessment and Prediction of Wind Turbine Noise
   
   !*********************************************** Model 1:
   !!! Nafnoise source code version see below 
   Frequency_cutoff = 10*U/PI/Chord
   Ke = 3.0/(4.0*p%Lturb) 
   Beta2 = 1-Mach*Mach
   ALPSTAR = AlphaNoise*PI/180.

   DO I=1,size(p%FreqList)
      IF (p%FreqList(I) <= Frequency_cutoff) THEN
         Directivity = DBARL
      ELSE
         Directivity = DBARH 
      ENDIF

      WaveNumber = 2.0*PI*p%FreqList(I)/U
      Kbar = WaveNumber*Chord/2.0
      Khat = WaveNumber/Ke
      ! mu = Mach*WaveNumber*Chord/2.0/Beta2

      SPLhigh = 10.*LOG10(p%AirDens*p%AirDens*p%SpdSound**4*p%Lturb*(d/2.)/ &
               (RObs*RObs)*(Mach**5)*tinooisess*tinooisess*(Khat**3)* &
               (1+Khat**2)**(-7./3.)*Directivity) + 78.4   ! ref a)
   !!!   SPLhigh = 10.*LOG10(p%Lturb*(d/2.)/ &
   !!!                  (RObs*RObs)*(Mach**5)*tinooisess*tinooisess*(WaveNumber**3) &
   !!!                  *(1+WaveNumber**2)**(-7./3.)*Directivity) + 181.3  
   
      SPLhigh = SPLhigh + 10.*LOG10(1+ 9.0*ALPSTAR*ALPSTAR)  ! Component due to angles of attack, ref a)   

      Sears = 1/(2.*PI*Kbar/Beta2+1/(1+2.4*Kbar/Beta2))      ! ref a)

   !!!   Sears = 1/(2.*PI*WaveNumber/Beta2+1/(1+2.4*WaveNumber/Beta2))  ! ref b) 
   
      LFC = 10*Sears*Mach*Kbar*Kbar/Beta2  ! ref a)
   !!!   LFC = 10*Sears*Mach*WaveNumber*WaveNumber/Beta2  ! ref b)
   
   !!!   IF (mu<(PI/4.0)) THEN                     ! ref b)
   !!!      SPLti(I) = SPLhigh + 10.*ALOG10(LFC)   ! ref b)
   !!!   ELSE                                      ! ref b)
   !!!      SPLti(I) = SPLhigh                     ! ref b)
   !!!ENDIF
      SPLti(I) = SPLhigh + 10.*LOG10(LFC/(1+LFC))
   
   ENDDO
   !!!*********************************************** end of Model 1

!  ! ********************************* Model 2:     
!  !Wei Jun Zhu et al - !Modeling of Aerodynamically Generated Noise From Wind Turbines 2005 paper
!        Beta2 = 1.d0-Mach**2; !  corresponding line:  Bsq = 1.d0 - Ma**2;
!         DO I=1,size(p%FreqList)
!       WaveNumber = PI*p%FreqList(I)*p%SpdSound/U !corresponding line: K = pi*Freq(i)*c/Vrel;  ! CarloS: This is a Mistake, c in this case is the Local Chord
!       Sears = (2.d0*PI*WaveNumber/Beta2 + (1.d0+2.4d0*WaveNumber/Beta2)**(-1))**(-1);
!     ! corresponding line: Ssq = (2.d0*pi*K/Bsq + (1.d0+2.4d0*K/Bsq)**(-1))**(-1);
!        LFC = 10.d0 * Sears*Mach*WaveNumber**2*Beta2**(-1);
!      ! corresponding line:  LFC = 10.d0 * Ssq*Ma*K**2*Bsq**(-1);
!     SPLti(I)=(p%AirDens*p%AirDens*p%SpdSound*p%SpdSound*p%Lturb*d)/(2*RObs*RObs)
!  !   SPLti(I)=SPLti(I)*(Mach**3)*(MeanVnoise**2)*(tinooisess**2)  
!     SPLti(I)=SPLti(I)*(Mach**3)*(tinooisess**2)                 
!  !   SPLti(I)=SPLti(I)*(Mach**3)*ufluct**2
!     SPLti(I)=(SPLti(I)*(WaveNumber**3)) / ((1+WaveNumber**2)**(7/3))
!     SPLti(I)=SPLti(I)*DBARH
!     SPLti(I)=10*log10(SPLti(I))+58.4
!      SPLti(I) = SPLti(I) + 10.*LOG10(LFC/(1+LFC))
!  ! SPLti(I)=10.d0*log10(DBARH*p%AirDens**2*p%SpdSound**2*p%Lturb*d/2.0*Mach**3*tinooisess**2* &
!  !WaveNumber**3*(1.d0+WaveNumber**2)**(-7.d0/3.d0)/RObs**2)+58.4d0 + 10.d0*log10(LFC/(1+LFC))
!  ! corresponding line:    SPLti(i)=10.d0*log10(Di_hi_fr*Density**2*co**2*Tbscale*L/2.0*Ma
!  !     & **3*Tbinten**2*K**3*(1.d0+K**2)**(-7.d0/3.d0)/Distance**2)+58.4d0
!  !     &    + 10.d0*log10(LFC/(1+LFC));  
!  !            !% ver2.!
!  !    Kh = 8.d0*pi*p%FreqList(i)*p%Lturb/(3.d0*U);
!  !          SPLti(i) = 10*log10(DBARH*p%Lturb*0.5*d*Mach**5*tinooisess**2*Kh**3*(1+Kh**2)**(-7/3)/RObs**2) +&
!  !              10*log10(10**18.13) + 10*log10(DBARH*LFC/(1+LFC));   
!  
!  ENDDO
!  ! ********************************* End of Model 2/ CarloSucameli: I think this model is wrong 



!!!!  ! ********************************* Model 3:     
!!!!  ! ref b) Lowson, Assessment and Prediction of Wind Turbine Noise
!!!!     Beta2 = 1.d0-Mach**2; !  corresponding line:  Bsq = 1.d0 - Ma**2;
!!!!      DO I=1,size(p%FreqList)
!!!!       WaveNumber = PI*p%FreqList(I)*Chord/U !corresponding line: K = pi*Freq(i)*c/Vrel;  
!!!!       Sears = (2.d0*PI*WaveNumber/Beta2 + (1.d0+2.4d0*WaveNumber/Beta2)**(-1))**(-1);
!!!!     ! corresponding line: Ssq = (2.d0*pi*K/Bsq + (1.d0+2.4d0*K/Bsq)**(-1))**(-1);
!!!!        LFC = 10.d0 * Sears*Mach*WaveNumber**2*Beta2**(-1);
!!!!      ! corresponding line:  LFC = 10.d0 * Ssq*Ma*K**2*Bsq**(-1);
!!!!     SPLti(I)=(p%AirDens*p%AirDens*p%SpdSound*p%SpdSound*p%Lturb*d)/(2*RObs*RObs)
!!!!     SPLti(I)=SPLti(I)*(Mach**3)*(MeanVnoise**2)*(tinooisess**2)  
!!!!     SPLti(I)=(SPLti(I)*(WaveNumber**3)) / ((1+WaveNumber**2)**(7./3.))
!!!!     SPLti(I)=SPLti(I)*DBARH
!!!!     SPLti(I)=10*log10(SPLti(I))+58.4
!!!!     SPLti(I) = SPLti(I) + 10.*LOG10(LFC/(1+LFC))
!!!!  
!!!!
!!!!     ENDDO
!!!!  ! ********************************* End of Model 3 
  
!!Buck&Oerlamans&Palo - !Experimental validation of a wind turbine turbulent inflow noise prediction code 2016 paper
!DO I=1,size(p%FreqList)
        ! IF (p%FreqList(I) <= Frequency_cutoff) THEN
        !    Directivity = DBARL
        ! ELSE
        !    Directivity = DBARH 
        ! ENDIF
 ! WaveNumber = 2.0*PI*p%FreqList(I)/U   ! (K)
        ! Kbar    = WaveNumber*Chord/2.0
        ! Khat    = WaveNumber/Ke
        ! SPLhigh = (  (p%AirDens**2) * (p%SpdSound**2) *d ) / (2*RObs*RObs)
! SPLhigh = SPLhigh * (Mach**3) * (dissip**(2/3)) * (WaveNumber**(-5/3)) * Directivity
        ! SPLhigh = 10.*LOG10(SPLhigh) + 77.6
        ! Sears   = 1/(2.*PI*Kbar/Beta2+1/(1+2.4*Kbar/Beta2))
        ! LFC = 10*Sears*(1+9.0*ALPSTAR*ALPSTAR)*Mach*Kbar*Kbar/Beta2
        ! SPLti(I) = SPLhigh + 10.*LOG10(LFC/(1+LFC))
        !ENDDO

! double commented lines are from FAST v4.0 aeroacoustics module. But Nafnoise version is used see above
!!   Mach = U/p%SpdSound
!!
!!IF (TINoise > 0) THEN
!!    Ums = (TINoise*MeanVNoise/100.)**2 ! mean square turbulence level
!!ELSE
!!    SPLti = 0.
!!    RETURN
!!ENDIF
!!
!! LTurb=60
!! LTurb=0.06
!!!  temporarily commented 
!!!  IF (FASTHH < 30.0) THEN
!!!      LTurb = 3.5*0.7*FASTHH ! Prediction sensitive to this parameter!
!!!  ELSE
!!!      LTurb = 3.5*21.
!!! ENDIF
!!
!!!LTurb = LTurb/100
!!
!!! Calculate directivity...?
!!!!!     ----------------------------
!!      CALL DIRECTL(Mach,THETA,PHI,DBARL,errStat2,errMsg2) !yes, assume that noise is low-freq in nature because turbulence length scale is large
!!    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )   
!!      CALL DIRECTH_LE(Mach,THETA,PHI,DBARH,errStat2,errMsg2)
!!    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
!!      IF (DBARH <= 0) THEN
!!    SPLti = 0.
!!          RETURN
!!      ENDIF
!!
!! Frequency_cutoff = 10*U/PI/Chord
!!
!!    IF (DBARL <= 0.) THEN
!!        SPLti = 0.
!!        RETURN
!!    ENDIF
!!
!!DO I=1,size(p%FreqList)
!!   IF (p%FreqList(I) <= Frequency_cutoff) THEN
!!       Directivity = DBARL
!!   ELSE
!!       Directivity = DBARH
!!   ENDIF
!!   WaveNumber = PI*p%FreqList(I)*Chord/U
!!  Beta2 = 1-Mach*Mach
!!   SPLhigh = 10.*LOG10(p%AirDens*p%AirDens*p%SpdSound*p%SpdSound*p%Lturb*(d/2.)/(RObs*RObs)*(Mach**3)*Ums* &
!!             (WaveNumber**3)*(1+WaveNumber**2)**(-7./3.)*Directivity) + 58.4
!!   Sears = 1/(2*PI*WaveNumber/Beta2+1/(1+2.4*WaveNumber/Beta2))
!!   LFC = 10*Sears*Mach*WaveNumber*WaveNumber/Beta2
!!   SPLti(I) = SPLhigh + 10.*LOG10(LFC/(1+LFC))
!!
!!ENDDO

END SUBROUTINE InflowNoise
!====================================================================================================
SUBROUTINE BLUNT(ALPSTAR,C,U ,THETA,PHI,L,R,H,PSI,p,d99Var2,dstarVar1,dstarVar2,SPLBLUNT,StallVal,errStat,errMsg)
  REAL(ReKi),                             INTENT(IN   )  :: ALPSTAR        ! AOA
  REAL(ReKi),                             INTENT(IN   )  :: C              ! Chord Length
  REAL(ReKi),                             INTENT(IN   )  :: U              ! Unoise
  REAL(ReKi),                             INTENT(IN   )  :: THETA          ! DIRECTIVITY ANGLE                     ---
  REAL(ReKi),                             INTENT(IN   )  :: PHI            ! DIRECTIVITY ANGLE                     ---
  REAL(ReKi),                             INTENT(IN   )  :: L              ! SPAN                                  METERS
  REAL(ReKi),                             INTENT(IN   )  :: R              ! SOURCE TO OBSERVER DISTANCE           METERS 
  REAL(ReKi),                             INTENT(IN   )  :: H              ! TRAILING EDGE BLUNTNESS              METERS
  REAL(ReKi),                             INTENT(IN   )  :: PSI            ! TRAILING EDGE ANGLE                  DEGREES 
  REAL(ReKi),                             INTENT(IN   )  :: d99Var2        !  
  REAL(ReKi),                             INTENT(IN   )  :: dstarVar1              !  
  REAL(ReKi),                             INTENT(IN   )  :: dstarVar2              !
  REAL(ReKi),                             INTENT(IN   )  :: StallVal       !< Stall angle at station i
  TYPE(AA_ParameterType),                 INTENT(IN   )  :: p              ! Parameters
  REAL(ReKi),DIMENSION(size(p%FreqList)), INTENT(  OUT)  :: SPLBLUNT       !
  INTEGER(IntKi),                         INTENT(  OUT)  :: errStat        ! Error status of the operation
  character(*),                           INTENT(  OUT)  :: errMsg         ! Error message if ErrStat /= ErrID_None
   ! Local variables
  integer(intKi)                         :: ErrStat2           ! temporary Error status
  character(ErrMsgLen)                   :: ErrMsg2            ! temporary Error message
  character(*), parameter                :: RoutineName = 'BLUNT'
  real(ReKi)                             :: STPPP    ! STROUHAL NUMBER                       ---
  real(ReKi)                             :: M        ! MACH NUMBER                           ---
  real(ReKi)                             :: RC       ! REYNOLDS NUMBER BASED ON CHORD        ---
  integer(intKi)                         :: I        ! I A generic index for DO loops.
  real(ReKi)                             :: DELTAP   ! PRESSURE SIDE BOUNDARY LAYER THICKNESS METERS
  real(ReKi)                             :: DSTRS    ! SUCTION SIDE DISPLACEMENT THICKNESS  METERS
  real(ReKi)                             :: DSTRP    ! PRESSURE SIDE DISPLACEMENT THICKNESS METERS
  real(ReKi)                             :: DBARH    ! HIGH FREQUENCY DIRECTIVITY           ---
  real(ReKi)                             :: DSTRAVG  ! AVERAGE DISPLACEMENT THICKNESS       METERS
  real(ReKi)                             :: HDSTAR   ! BLUNTNESS OVER AVERAGE DISPLACEMENT THICKNESS   ---
  real(ReKi)                             :: DSTARH   ! AVERAGE DISPLACEMENT THICKNESS OVER TRAILING EDGE BLUNTNESS       ---
  real(ReKi)                             :: ATERM    ! USED TO COMPUTE PEAK STROUHAL NO.    ---
  real(ReKi)                             :: STPEAK   ! PEAK STROUHAL NUMBER                  ---
  real(ReKi)                             :: ETA      ! RATIO OF STROUHAL NUMBERS             ---
  real(ReKi)                             :: HDSTARL  ! MINIMUM ALLOWED VALUE OF HDSTAR       ---
  real(ReKi)                             :: G514     ! G5 EVALUATED AT PSI=14.0              DB
  real(ReKi)                             :: HDSTARP  ! MODIFIED VALUE OF HDSTAR              ---
  real(ReKi)                             :: G50      ! G5 EVALUATED AT PSI=0.0               DB
  real(ReKi)                             :: G4       ! SCALED SPECTRUM LEVEL                 DB
  !   real(ReKi)                         :: G5       ! SPECTRUM SHAPE FUNCTION               DB
  REAL(ReKi),DIMENSION(size(p%FreqList)) :: G5       ! SPECTRUM SHAPE FUNCTION               DB ! corrected (EB_DTU)
  real(ReKi)                             :: G5Sum       ! SPECTRUM SHAPE FUNCTION               DB
  real(ReKi)                             :: SCALE    ! SCALING FACTOR                        ---

   ErrStat = ErrID_None
   ErrMsg  = ""

    ! Reynolds number and mach number
        M          = U  / p%SpdSound
        RC         = U  * C/p%KinVisc
    ! Compute boundary layer thicknesses
    IF (p%X_BLMethod .eq. X_BLMethod_Tables) THEN
        DELTAP = d99Var2
        DSTRS  = dstarVar1
        DSTRP  = dstarVar2
    ELSE
        CALL THICK(C,RC,ALPSTAR,p,DELTAP,DSTRS,DSTRP,StallVal,errStat2,errMsg2)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
    ENDIF
    ! Compute average displacement thickness
    DSTRAVG = (DSTRS + DSTRP) / 2.
    HDSTAR  = H / DSTRAVG
    DSTARH = 1. /HDSTAR
    ! Compute directivity function
    CALL DIRECTH_TE(M,THETA,PHI,DBARH,errStat2,errMsg2)
    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
    IF (DBARH <= 0) THEN
        SPLBLUNT = 0.
        RETURN
    ENDIF
    ! Compute peak strouhal number                                               eq 72 in BPM Airfoil Self-noise and Prediction paper
    ATERM  = .212 - .0045 * PSI
    IF (HDSTAR .GE. .2) &
        STPEAK    = ATERM / (1.+.235*DSTARH-.0132*DSTARH**2)                    ! this is what it used to be in nafnoise and fast noise module
    !!  STPEAK    = ATERM / (1+0.235*(DSTARH)**(-1)-0.0132*DSTARH**(-2)) ! check if this one is correct (EB_DTU) 
    IF (HDSTAR .LT. .2) &
        STPEAK    = .1 * HDSTAR + .095 - .00243 * PSI
    ! Compute scaled spectrum level                                              eq 74 of BPM Airfoil Self-noise and Prediction paper
    IF (HDSTAR .LE. 5.) G4=17.5*LOG10(HDSTAR)+157.5-1.114*PSI
    IF (HDSTAR .GT. 5.) G4=169.7 - 1.114 * PSI
    ! For each frequency, compute spectrum shape referenced to 0 db
    SCALE = 10. * LOG10(M**5.5*H*DBARH*L/R**2)
    G5Sum=0.0_Reki
    DO I=1,SIZE(p%FreqList)
        STPPP    = p%FreqList(I) * H / U
        ETA      = LOG10(STPPP/STPEAK)
        HDSTARL = HDSTAR
        CALL G5COMP(HDSTARL,ETA,G514,errStat2,errMsg2 )                          ! compute G5 for Phi=14deg
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
        HDSTARP = 6.724 * HDSTAR **2-4.019*HDSTAR+1.107                         ! eq 82 from BPM Airfoil Self-noise and Prediction paper
        CALL G5COMP(HDSTARP,ETA,G50,errStat2,errMsg2 )                           ! recompute G5 for Phi=0deg
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
        G5(I) = G50 + .0714 * PSI * (G514-G50)                                   ! interpolate G5 from G50 and G514
        IF (G5(I) .GT. 0.) G5(I) = 0.
        G5Sum = 10**(G5(I)/10)+G5Sum     ! to be subtracted
        SPLBLUNT(I) = G4 + G5(I) + SCALE - 10*log10(1/G5Sum) ! equation mentioned there is plus but it is stated subtract, thus ''- 10*log10(1/G5Sum)'' 
    end do
END SUBROUTINE Blunt
!====================================================================================================
SUBROUTINE G5COMP(HDSTAR,ETA,G5,errStat,errMsg)
    REAL(ReKi),          INTENT(IN   )  :: HDSTAR        !<
    REAL(ReKi),          INTENT(IN   )  :: ETA           !< 
    REAL(ReKi),          INTENT(  OUT)  :: G5            !< 
    INTEGER(IntKi),      INTENT(  OUT)  :: errStat       !< Error status of the operation
    CHARACTER(*),        INTENT(  OUT)  :: errMsg        !< Error message if ErrStat /= ErrID_None
    ! Local variables
!    INTEGER(intKi)                                                 :: ErrStat2           ! temporary Error status
!    CHARACTER(ErrMsgLen)                                           :: ErrMsg2            ! temporary Error message
    CHARACTER(*), parameter                                        :: RoutineName = 'BLUNT'
    real(ReKi)                                    :: K 
    real(ReKi)                                    :: M
    real(ReKi)                                    :: MU
    real(ReKi)                                    :: ETALIMIT
    real(ReKi)                                    :: ETA0
    ErrStat = ErrID_None
    ErrMsg  = ""
    IF ( HDSTAR .LT. .25)                          MU = .1211                    ! begin eq 78 from BPM Airfoil Self-noise and Prediction paper
    IF ((HDSTAR .GT. .25).AND.(HDSTAR .LE. .62))   MU =-.2175*HDSTAR + .1755
    IF ((HDSTAR .GT. .62).AND.(HDSTAR .LT. 1.15))  MU =-.0308*HDSTAR + .0596
    IF ( HDSTAR .GE. 1.15)                         MU = .0242                    ! end
    IF ( HDSTAR .LE. .02 )                         M = 0.0                       ! begin eq 79 from BPM Airfoil Self-noise and Prediction paper
    IF ((HDSTAR .GE. .02 ).AND.(HDSTAR .LT. .5))   M = 68.724*HDSTAR - 1.35   
    IF ((HDSTAR .GT. .5  ).AND.(HDSTAR .LE. .62))  M = 308.475*HDSTAR - 121.23
    IF ((HDSTAR .GT. .62 ).AND.(HDSTAR .LE. 1.15)) M = 224.811*HDSTAR - 69.354
    IF ((HDSTAR .GT. 1.15).AND.(HDSTAR .LT. 1.2))  M = 1583.28*HDSTAR - 1631.592
    IF ( HDSTAR .GT. 1.2 )                         M = 268.344
    IF ( M      .LT. 0.0 )                         M = 0.0                       ! end
    ETA0 = -SQRT((M*M*MU**4)/(6.25+M*M*MU*MU))                                   ! eq 80 from BPM Airfoil Self-noise and Prediction paper
    K    = 2.5*SQRT(1.-(ETA0/MU)**2)-2.5-M*ETA0                                  ! eq 81 from BPM Airfoil Self-noise and Prediction paper
    ETALIMIT = 0.03615995                                                        ! one of the bounds given in eq 76 of BPM Airfoil Self-noise and Prediction paper
    IF (ETA .LE. ETA0)                      G5 = M * ETA + K                     ! begin eq 76 from BPM Airfoil Self-noise and Prediction paper
    IF((ETA.GT.ETA0).AND.(ETA .LE. 0.))     G5 = 2.5*SQRT(1.-(ETA/MU)**2)-2.5
    IF((ETA.GT.0.  ).AND.(ETA.LE.ETALIMIT)) G5 = SQRT(1.5625-1194.99*ETA**2)-1.25
    IF (ETA.GT.ETALIMIT)                    G5 = -155.543 * ETA + 4.375          ! end
END SUBROUTINE G5Comp
!====================================================================================================
!> This subroutine defines the curve fit corresponding to the a-curve for the minimum allowed reynolds number.
SUBROUTINE AMIN(A,AMINA)
    REAL(ReKi),                             INTENT(IN   )  :: A
    REAL(ReKi),                             INTENT(OUT  )  :: AMINA
    REAL(ReKi) :: X1
    X1 = ABS(A)
    IF (X1 .LE. .204) AMINA=SQRT(67.552-886.788*X1**2)-8.219
    IF((X1 .GT. .204).AND.(X1 .LE. .244))AMINA=-32.665*X1+3.981
    IF (X1 .GT. .244)AMINA=-142.795*X1**3+103.656*X1**2-57.757*X1+6.006
END SUBROUTINE AMIN
!====================================================================================================
!> This subroutine defines the curve fit corresponding to the a-curve for the maximum allowed reynolds number.
SUBROUTINE AMAX(A,AMAXA)
    REAL(ReKi),                             INTENT(IN   )  :: A
    REAL(ReKi),                             INTENT(OUT  )  :: AMAXA
    REAL(ReKi) :: X1
    X1 = ABS(A)
    IF (X1 .LE. .13)AMAXA=SQRT(67.552-886.788*X1**2)-8.219
    IF((X1 .GT. .13).AND.(X1 .LE. .321))AMAXA=-15.901*X1+1.098
    IF (X1 .GT. .321)AMAXA=-4.669*X1**3+3.491*X1**2-16.699*X1+1.149
END SUBROUTINE AMAX
!====================================================================================================
!> This subroutine defines the curve fit corresponding to the b-curve for the minimum allowed reynolds number.
SUBROUTINE BMIN(B,BMINB)
    REAL(ReKi),                             INTENT(IN   )  :: B
    REAL(ReKi),                             INTENT(OUT  )  :: BMINB
    REAL(ReKi) :: X1
    X1 = ABS(B)
    IF (X1 .LE. .13)BMINB=SQRT(16.888-886.788*X1**2)-4.109
    IF((X1 .GT. .13).AND.(X1 .LE. .145))BMINB=-83.607*X1+8.138
    IF (X1.GT..145)BMINB=-817.81*X1**3+355.21*X1**2-135.024*X1+10.619
END SUBROUTINE BMin
!====================================================================================================
!> Define the curve fit corresponding to the b-curve for the maximum allowed reynolds number.
SUBROUTINE BMAX(B,BMAXB)
    REAL(ReKi),   INTENT(IN   )  :: B
    REAL(ReKi),   INTENT(OUT  )  :: BMAXB
    REAL(ReKi) :: X1
    X1 = ABS(B)
    IF (X1 .LE. .1) BMAXB=SQRT(16.888-886.788*X1**2)-4.109
    IF((X1 .GT. .1).AND.(X1 .LE. .187))BMAXB=-31.313*X1+1.854
    IF (X1.GT..187)BMAXB=-80.541*X1**3+44.174*X1**2-39.381*X1+2.344
END SUBROUTINE BMax
!====================================================================================================
!> Determine where the a-curve takes on a value of -20 db.
SUBROUTINE A0COMP(RC,A0)
    REAL(ReKi),   INTENT(IN   )  :: RC
    REAL(ReKi),   INTENT(OUT  )  :: A0
    IF (RC .LT. 9.52E+04) A0 = .57
    IF ((RC .GE. 9.52E+04).AND.(RC .LT. 8.57E+05)) &
        A0 = (-9.57E-13)*(RC-8.57E+05)**2 + 1.13
    IF (RC .GE. 8.57E+05) A0 = 1.13
END SUBROUTINE A0COMP
!====================================================================================================
!> Compute zero angle of attack boundary layer thickness (meters) and reynolds number
SUBROUTINE THICK(C,RC,ALPSTAR,p,DELTAP,DSTRS,DSTRP,StallVal,errStat,errMsg)
!!       VARIABLE NAME               DEFINITION                  UNITS
!!       -------------               ----------                  -----
!!       ALPSTAR            ANGLE OF ATTACK                    DEGREES
!!       C                  CHORD LENGTH                        METERS
!!       C0                 SPEED OF SOUND                    METERS/SEC
!!       DELTA0             BOUNDARY LAYER THICKNESS AT
!!                            ZERO ANGLE OF ATTACK              METERS
!!       DELTAP             PRESSURE SIDE BOUNDARY LAYER
!!                            THICKNESS                         METERS
!!       DSTR0              DISPLACEMENT THICKNESS AT ZERO
!!                            ANGLE OF ATTACK                   METERS
!!       DSTRP              PRESSURE SIDE DISPLACEMENT
!!                            THICKNESS                         METERS
!!       DSTRS              SUCTION SIDE DISPLACEMENT
!!                            THICKNESS                         METERS
!!       ITRIP              TRIGGER FOR BOUNDARY LAYER TRIPPING  ---
!!       RC                 REYNOLDS NUMBER BASED ON CHORD       ---
!!       U                  FREESTREAM VELOCITY                METERS/SEC
!!       KinViscosity       KINEMATIC VISCOSITY                M2/SEC
    REAL(ReKi),                INTENT(IN   )  :: ALPSTAR        !< AOA
    REAL(ReKi),                INTENT(IN   )  :: C              !< Chord Length
    REAL(ReKi),                INTENT(IN   )  :: RC             !< RC= U*C/KinViscosity
    TYPE(AA_ParameterType),    INTENT(IN   )  :: p              !< Parameters
    REAL(ReKi),                INTENT(  OUT)  :: DELTAP         !<
    REAL(ReKi),                INTENT(  OUT)  :: DSTRS          !<
    REAL(ReKi),                INTENT(  OUT)  :: DSTRP          !<
    REAL(ReKi),                INTENT(IN   )  :: StallVal       !< Stall angle at station i
    INTEGER(IntKi),            INTENT(  OUT)  :: errStat        !< Error status of the operation
    character(*),              INTENT(  OUT)  :: errMsg         !< Error message if ErrStat /= ErrID_None
    ! Local variables
!    integer(intKi)          :: ErrStat2           ! temporary Error status
!    character(ErrMsgLen)    :: ErrMsg2            ! temporary Error message
    character(*), parameter :: RoutineName = 'Thick'
    real(ReKi)              :: DELTA0              ! BOUNDARY LAYER THICKNESS AT ZERO ANGLE OF ATTACK METERS
    real(ReKi)              :: DSTR0      ! DISPLACEMENT THICKNESS AT ZERO   ANGLE OF ATTACK METERS
    ErrStat = ErrID_None
    ErrMsg  = ""
    ! Boundary layer thickness
    DELTA0                     = 10.**(1.6569-0.9045*LOG10(RC)+0.0596*LOG10(RC)**2)*C ! (untripped)         Eq. (5) of [1]
    IF (p%ITRIP .GT. 0) DELTA0 = 10.**(1.892 -0.9045*LOG10(RC)+0.0596*LOG10(RC)**2)*C ! (heavily tripped)   Eq. (2) of [1]
    IF (p%ITRIP .EQ. 2) DELTA0=.6*DELTA0
    ! Pressure side boundary layer thickness, Eq (8) of [1]
    DELTAP   = 10.**(-.04175*ALPSTAR+.00106*ALPSTAR**2)*DELTA0
    ! Compute zero angle of attack displacement thickness
    IF ((p%ITRIP .EQ. 1) .OR. (p%ITRIP .EQ. 2)) THEN
        ! Heavily tripped, Eq. (3) of [1]
        IF (RC .LE. .3E+06) DSTR0 = .0601 * RC **(-.114)*C
        IF (RC .GT. .3E+06) &
            DSTR0=10.**(3.411-1.5397*LOG10(RC)+.1059*LOG10(RC)**2)*C
        ! Lightly tripped
        IF (p%ITRIP .EQ. 2) DSTR0 = DSTR0 * .6
    ELSE
        ! Untripped, Eq. (6) of [1]
        DSTR0=10.**(3.0187-1.5397*LOG10(RC)+.1059*LOG10(RC)**2)*C
    ENDIF
    ! Pressure side displacement thickness, Eq. (9) of [1]
    DSTRP   = 10.**(-.0432*ALPSTAR+.00113*ALPSTAR**2)*DSTR0
    !      IF (p%ITRIP .EQ. 3) DSTRP = DSTRP * 1.48 ! commented since itrip is never 3 check if meant 2.(EB_DTU)
    ! Suction side displacement thickness
    IF (p%ITRIP .EQ. 1) THEN
        ! Heavily tripped, Eq. (12) of [1]
        IF (ALPSTAR .LE. 5.) DSTRS=10.**(.0679*ALPSTAR)*DSTR0
        IF((ALPSTAR .GT. 5.).AND.(ALPSTAR .LE. StallVal)) &
            DSTRS = .381*10.**(.1516*ALPSTAR)*DSTR0
        IF (ALPSTAR .GT. StallVal)DSTRS=14.296*10.**(.0258*ALPSTAR)*DSTR0
    ELSE
        ! Untripped or lightly tripped, Eq. (15) of [1]
        IF (ALPSTAR .LE. 7.5)DSTRS =10.**(.0679*ALPSTAR)*DSTR0
        IF((ALPSTAR .GT. 7.5).AND.(ALPSTAR .LE. StallVal)) &
            DSTRS = .0162*10.**(.3066*ALPSTAR)*DSTR0
        IF (ALPSTAR .GT. StallVal) DSTRS = 52.42*10.**(.0258*ALPSTAR)*DSTR0
    ENDIF
END SUBROUTINE Thick
!====================================================================================================
!> This subroutine computes the high frequency directivity function for the trailing edge
SUBROUTINE DIRECTH_TE(M,THETA,PHI,DBAR, errStat, errMsg)
    REAL(ReKi),        INTENT(IN   ) :: THETA      !
    REAL(ReKi),        INTENT(IN   ) :: PHI        !
    REAL(ReKi),        INTENT(IN   ) :: M          !
    REAL(ReKi),        INTENT(  OUT) :: DBAR       !
    INTEGER(IntKi),    INTENT(  OUT) :: errStat    ! Error status of the operation
    character(*),      INTENT(  OUT) :: errMsg     ! Error message if ErrStat /= ErrID_None
    ! Local variables
    character(*), parameter :: RoutineName = 'Directh_te'
    real(ReKi)              :: MC
    real(ReKi)              :: DEGRAD
    real(ReKi)              :: PHIR
    real(ReKi)              :: THETAR
    ErrStat = ErrID_None
    ErrMsg  = ""
    DEGRAD = .017453
    MC     = .8 * M
    THETAR = THETA * DEGRAD
    PHIR   = PHI * DEGRAD
    DBAR   = 2.*SIN(THETAR/2.)**2*SIN(PHIR)**2/((1.+M*COS(THETAR))* (1.+(M-MC)*COS(THETAR))**2)    ! eq B1 in BPM Airfoil Self-noise and Prediction paper
END SUBROUTINE DIRECTH_TE

!====================================================================================================
!> This subroutine computes the high frequency directivity function for the leading edge
SUBROUTINE DIRECTH_LE(M,THETA,PHI,DBAR, errStat, errMsg)
    REAL(ReKi),        INTENT(IN   ) :: THETA      !
    REAL(ReKi),        INTENT(IN   ) :: PHI        !
    REAL(ReKi),        INTENT(IN   ) :: M          !
    REAL(ReKi),        INTENT(  OUT) :: DBAR       !
    INTEGER(IntKi),    INTENT(  OUT) :: errStat    ! Error status of the operation
    character(*),      INTENT(  OUT) :: errMsg     ! Error message if ErrStat /= ErrID_None
    ! Local variables
    character(*), parameter :: RoutineName = 'Directh_le'
    real(ReKi)              :: DEGRAD
    real(ReKi)              :: PHIR
    real(ReKi)              :: THETAR
    ErrStat = ErrID_None
    ErrMsg  = ""
    DEGRAD = .017453
    THETAR = THETA * DEGRAD
    PHIR   = PHI * DEGRAD
    DBAR   = 2.*COS(THETAR/2.)**2*SIN(PHIR)**2/(1.+M*COS(THETAR))**3 
END SUBROUTINE DIRECTH_LE

!====================================================================================================
!> This subroutine computes the high frequency directivity function for the input observer location
! Paper: 
SUBROUTINE DIRECTL(M,THETA,PHI,DBAR, errStat, errMsg)
    REAL(ReKi),           INTENT(IN   ) :: THETA      !<
    REAL(ReKi),           INTENT(IN   ) :: PHI        !<
    REAL(ReKi),           INTENT(IN   ) :: M          !<
    REAL(ReKi),           INTENT(  OUT) :: DBAR       !<
    INTEGER(IntKi),       INTENT(  OUT) :: errStat    !< Error status of the operation
    character(*),         INTENT(  OUT) :: errMsg     !< Error message if ErrStat /= ErrID_None
    ! Local variables
    character(*), parameter :: RoutineName = 'DirectL'
    real(ReKi)              :: MC
    real(ReKi)              :: DEGRAD
    real(ReKi)              :: PHIR
    real(ReKi)              :: THETAR
    ErrStat = ErrID_None
    ErrMsg  = "" 
    !   This subroutine computes the low frequency directivity function for the input observer location
    DEGRAD  = .017453
    MC     = .8 * M
    THETAR = THETA * DEGRAD
    PHIR   = PHI * DEGRAD
    DBAR = (SIN(THETAR)*SIN(PHIR))**2/(1.+M*COS(THETAR))**4                   ! eq B2 in BPM Airfoil Self-noise and Prediction paper
END SUBROUTINE DIRECTL
!==================================================================================================================================!
!===============================  Simplified Guidati Inflow Turbulence Noise Addition =============================================!
!==================================================================================================================================!
! Uses simple correction for turbulent inflow noise from Moriarty et. al 2005
! Paper: Prediction of Turbulent Inflow and Trailing-Edge Noise for Wind Turbines, by Moriarty, Guidati, and Migliore
SUBROUTINE Simple_Guidati(U,Chord,thick_10p,thick_1p,p,SPLti,errStat,errMsg)
    REAL(ReKi),                             INTENT(IN   )  :: U              ! Vrel
    REAL(ReKi),                             INTENT(IN   )  :: Chord          ! Chord Length
    REAL(ReKi),                             INTENT(IN   )  :: thick_10p      ! 
    REAL(ReKi),                             INTENT(IN   )  :: thick_1p       ! 
    TYPE(AA_ParameterType),                 INTENT(IN   )  :: p              ! Parameters
    REAL(ReKi),DIMENSION(size(p%FreqList)), INTENT(  OUT)  :: SPLti          !
    INTEGER(IntKi),                         INTENT(  OUT)  :: errStat        ! Error status of the operation
    character(*),                           INTENT(  OUT)  :: errMsg         ! Error message if ErrStat /= ErrID_None
    ! local variables
!    integer(intKi)          :: ErrStat2       ! temporary Error status
!    character(ErrMsgLen)    :: ErrMsg2        ! temporary Error message
    character(*), parameter :: RoutineName = 'Simple_Guidati'
    INTEGER(intKi)          :: loop1       ! temporary
    REAL(ReKi)              :: TI_Param    ! Temporary variable thickness ratio dependent
    REAL(ReKi)              :: slope       ! Temporary variable thickness ratio dependent
    
    ErrStat = ErrID_None
    ErrMsg  = "" 

    TI_Param = thick_1p + thick_10p                                     ! Eq 2 
    slope = 1.123*TI_Param + 5.317*TI_Param*TI_Param                    ! Eq 3 
    do loop1 =1,size(p%FreqList)
        SPLti(loop1) = -slope*(2*PI*p%FreqList(loop1)*chord/U + 5.0d0)  ! Eq 4 
    enddo   ! Outputs Delta_SPL, the difference in SPL between the airfoil and a flat plate.
END SUBROUTINE Simple_Guidati
!==================================================================================================================================!
!================================ Turbulent Boundary Layer Trailing Edge Noise ====================================================!
!=================================================== TNO START ====================================================================!
SUBROUTINE TBLTE_TNO(U,THETA,PHI,D,R,Cfall,d99all,EdgeVelAll,p,SPLP,SPLS,SPLALPH,SPLTBL,errStat,errMsgn)
   USE TNO, only: SPL_integrate
    REAL(ReKi),                               INTENT(IN   ) :: U          !< Unoise                 (m/s)
    REAL(ReKi),                               INTENT(IN   ) :: THETA      !< DIRECTIVITY ANGLE      (deg)
    REAL(ReKi),                               INTENT(IN   ) :: PHI        !< DIRECTIVITY ANGLE      (deg)
    REAL(ReKi),                               INTENT(IN   ) :: D          !< SPAN                   (m)
    REAL(ReKi),                               INTENT(IN   ) :: R          !< SOURCE TO OBSERVER DISTANCE (m)
    REAL(ReKi),DIMENSION(2),                  INTENT(IN   ) :: Cfall      !< Skin friction coefficient   (-)
    REAL(ReKi),DIMENSION(2),                  INTENT(IN   ) :: d99all     !< 
    REAL(ReKi),DIMENSION(2),                  INTENT(IN   ) :: EdgeVelAll !< 
    TYPE(AA_ParameterType),                   INTENT(IN   ) :: p          !< Noise Module Parameters
    REAL(ReKi),DIMENSION(size(p%FreqList)),   INTENT(IN   ) :: SPLALPH    !< SOUND PRESSURE LEVEL DUE TO ANGLE OF ATTACK CONTRIBUTION (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),   INTENT(  OUT) :: SPLP       !< SOUND PRESSURE LEVEL DUE TO PRESSURE SIDE OF AIRFOIL (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),   INTENT(  OUT) :: SPLS       !< SOUND PRESSURE LEVEL DUE TO SUCTION SIDE OF AIRFOIL  (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),   INTENT(  OUT) :: SPLTBL     !< TOTAL SOUND PRESSURE LEVEL DUE TO TBLTE MECHANISM    (db)
    INTEGER(IntKi),                           INTENT(  OUT) :: errStat    !< Error status of the operation
    character(*),                             INTENT(  OUT) :: errMsgn    !< Error message if ErrStat /= ErrID_None
    ! Local variables
    integer(intKi)          :: ErrStat2                                                                                  ! temporary Error status
    character(ErrMsgLen)    :: ErrMsg2                                                                                   ! temporary Error message
    character(*), parameter :: RoutineName                              = 'TBLTE_TNO'
    REAL(ReKi) :: answer
    REAL(ReKi) :: Spectrum
    REAL(ReKi) :: freq(size(p%FreqList))
    REAL(ReKi) :: SPL_press,SPL_suction
    REAL(ReKi) :: band_width,band_ratio
    REAL(ReKi) :: DBARH
    REAL(ReKi) :: P1,P2,P4
    INTEGER (4)  :: n_freq
    INTEGER (4)  :: i_omega

      ! Variables passed to integration routine
   real(ReKi)  :: int_limits(2)  !< Lower and upper integration limits
   real(ReKi)  :: Mach        !< Mach number
   real(ReKi)  :: omega

    ! Init
    n_freq  = size(p%FreqList)
    freq    = p%FreqList
    ErrStat = ErrID_None
    ErrMsgn = ""
    ! Body of TNO 
    band_ratio = 2.**(1./3.)

    ! Mach number
    Mach = U  / p%SpdSound

    ! Directivity function
    CALL DIRECTH_TE(REAL(Mach,ReKi),THETA,PHI,DBARH,errStat2,errMsg2)
    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsgn, RoutineName )
 
    do i_omega = 1,n_freq
        omega = 2.*pi*freq(i_omega)
        !integration limits
        int_limits(1) = 0.0e0
        int_limits(2) = 10*omega/(Mach*p%SpdSound)
        ! Convert to third octave
        band_width = freq(i_omega)*(sqrt(band_ratio)-1./sqrt(band_ratio)) * 4. * pi
        IF (Cfall(1) .GT. 0.) THEN
            answer = SPL_integrate(omega=omega,limits=int_limits,ISSUCTION=.true.,        &
                     Mach=Mach,SpdSound=p%SpdSound,AirDens=p%AirDens,KinVisc=p%KinVisc,   &
                     Cfall=Cfall,d99all=d99all,EdgeVelAll=EdgeVelAll)
            Spectrum = D/(4.*pi*R**2)*answer
            SPL_suction = 10.*log10(Spectrum*DBARH/2.e-5/2.e-5)
            SPLS(i_omega) = SPL_suction + 10.*log10(band_width)
        ENDIF

        IF (Cfall(2) .GT. 0.) THEN
            answer = SPL_integrate(omega=omega,limits=int_limits,ISSUCTION=.FALSE.,       &
                     Mach=Mach,SpdSound=p%SpdSound,AirDens=p%AirDens,KinVisc=p%KinVisc,   &
                     Cfall=Cfall,d99all=d99all,EdgeVelAll=EdgeVelAll)
            Spectrum = D/(4.*pi*R**2)*answer
            SPL_press = 10.*log10(Spectrum*DBARH/2.e-5/2.e-5)
            SPLP(i_omega) = SPL_press + 10.*log10(band_width)
        ENDIF

        ! Sum the noise sources SPLALPH is BPM value
        IF (SPLP(i_omega)    .LT. -100.) SPLP(i_omega)    = -100.
        IF (SPLS(i_omega)    .LT. -100.) SPLS(i_omega)    = -100.

        P1  = 10.**(SPLP(i_omega) / 10.)
        P2  = 10.**(SPLS(i_omega) / 10.)
        P4  = 10.**(SPLALPH(i_omega) / 10.)
        SPLTBL(i_omega) = 10. * LOG10(P1 + P2 + P4)
    enddo
END SUBROUTINE TBLTE_TNO


!====================================================================================================
SUBROUTINE BL_Param_Interp(p,m,U,AlphaNoise,C,whichairfoil, errStat, errMsg)
  TYPE(AA_ParameterType),                INTENT(IN   ) :: p              !< Parameters
  TYPE(AA_MiscVarType),              INTENT(INOUT)     :: m              !< misc/optimization data (not defined in submodules)
  REAL(ReKi),                        INTENT(IN   )     :: U              !< METERS/SEC
  REAL(ReKi),                        INTENT(IN   )     :: AlphaNoise     !< Angle of Attack                           DEG
  REAL(ReKi),                        INTENT(IN   )     :: C              !< Chord                                     METERS
  integer(intKi),                        INTENT(IN   ) :: whichairfoil   !< whichairfoil
  integer(IntKi),                intent(  out)         :: ErrStat        !< Error status of the operation
  character(*),                  intent(  out)         :: ErrMsg         !< Error message if ErrStat /= ErrID_None
  character(*), parameter :: RoutineName = 'BL_Param_Interp'
  REAL(ReKi)              :: redif1,redif2,aoadif1,aoadif2,xx1,xx2,RC
  INTEGER(intKi)          :: loop1,loop2
  logical                 :: re_flag
  ErrStat = ErrID_None
  ErrMsg  = ""

  !!!! this if is not used but if necessary two sets of tables can be populated for tripped and untripped cases
  RC = U  * C/p%KinVisc       ! REYNOLDS NUMBER BASED ON  CHORD

  re_flag = .FALSE.
  DO loop1=1,size(p%ReListBL)-1
      IF (   (RC.le.p%ReListBL(loop1+1)) .and. (RC.gt.p%ReListBL(loop1))  ) then
          re_flag = .TRUE.
          redif1=abs(RC-p%ReListBL(loop1+1))
          redif2=abs(RC-p%ReListBL(loop1))
          DO loop2=1,size(p%AOAListBL)-1

              if (  (AlphaNoise.le.p%AOAListBL(loop2+1)) .and. (AlphaNoise.gt.p%AOAListBL(loop2))  ) then
                  aoadif1=abs(AlphaNoise-p%AOAListBL(loop2+1))
                  aoadif2=abs(AlphaNoise-p%AOAListBL(loop2))

                  xx1=( p%dstarall1(loop2,loop1+1,whichairfoil)*redif2+p%dstarall1(loop2,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  xx2=( p%dstarall1(loop2+1,loop1+1,whichairfoil)*redif2+p%dstarall1(loop2+1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  m%dstarVar(1)=(xx1*aoadif1+xx2*aoadif2) / (aoadif1+aoadif2)

                  xx1=( p%dstarall2(loop2,loop1+1,whichairfoil)*redif2+p%dstarall2(loop2,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  xx2=( p%dstarall2(loop2+1,loop1+1,whichairfoil)*redif2+p%dstarall2(loop2+1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  m%dstarVar(2)=(xx1*aoadif1+xx2*aoadif2) / (aoadif1+aoadif2)

                  xx1=( p%d99all1(loop2,loop1+1,whichairfoil)*redif2+p%d99all1(loop2,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  xx2=( p%d99all1(loop2+1,loop1+1,whichairfoil)*redif2+p%d99all1(loop2+1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  m%d99Var(1)=(xx1*aoadif1+xx2*aoadif2) / (aoadif1+aoadif2)

                  xx1=( p%d99all2(loop2,loop1+1,whichairfoil)*redif2+p%d99all2(loop2,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  xx2=( p%d99all2(loop2+1,loop1+1,whichairfoil)*redif2+p%d99all2(loop2+1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  m%d99Var(2)=(xx1*aoadif1+xx2*aoadif2) / (aoadif1+aoadif2)

                  xx1=( p%Cfall1(loop2,loop1+1,whichairfoil)*redif2+p%Cfall1(loop2,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  xx2=( p%Cfall1(loop2+1,loop1+1,whichairfoil)*redif2+p%Cfall1(loop2+1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  m%CfVar(1)=(xx1*aoadif1+xx2*aoadif2) / (aoadif1+aoadif2)

                  xx1=( p%Cfall2(loop2,loop1+1,whichairfoil)*redif2+p%Cfall2(loop2,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  xx2=( p%Cfall2(loop2+1,loop1+1,whichairfoil)*redif2+p%Cfall2(loop2+1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  m%CfVar(2)=(xx1*aoadif1+xx2*aoadif2) / (aoadif1+aoadif2)

                  xx1=( p%EdgeVelRat1(loop2,loop1+1,whichairfoil)*redif2+p%EdgeVelRat1(loop2,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  xx2=( p%EdgeVelRat1(loop2+1,loop1+1,whichairfoil)*redif2+p%EdgeVelRat1(loop2+1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  m%EdgeVelVar(1)=(xx1*aoadif1+xx2*aoadif2) / (aoadif1+aoadif2)

                  xx1=( p%EdgeVelRat2(loop2,loop1+1,whichairfoil)*redif2+p%EdgeVelRat2(loop2,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  xx2=( p%EdgeVelRat2(loop2+1,loop1+1,whichairfoil)*redif2+p%EdgeVelRat2(loop2+1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  m%EdgeVelVar(2)=(xx1*aoadif1+xx2*aoadif2) / (aoadif1+aoadif2)

                  return ! We exit the routine !
              endif
              if (loop2 .eq. (size(p%AOAListBL)-1) ) then

                  if (AlphaNoise .gt. p%AOAListBL(size(p%AOAListBL))) then
                      CALL WrScr( 'Warning AeroAcoustics Module - Angle of attack (AoA) range is not in the range provided by the user')
                      CALL WrScr( 'Station '// trim(num2lstr(whichairfoil)) )
                      CALL WrScr( 'Airfoil AoA '//trim(num2lstr(AlphaNoise))//'; Using the closest AoA '//trim(num2lstr(p%AOAListBL(loop2+1))))
                      m%dStarVar  (1) = ( p%dstarall1  (loop2+1,loop1+1,whichairfoil)*redif2 + p%dstarall1  (loop2+1,loop1,whichairfoil)*redif1 )/(redif1+redif2)
                      m%dStarVar  (2) = ( p%dstarall2  (loop2+1,loop1+1,whichairfoil)*redif2 + p%dstarall2  (loop2+1,loop1,whichairfoil)*redif1 )/(redif1+redif2)
                      m%d99Var    (1) = ( p%d99all1    (loop2+1,loop1+1,whichairfoil)*redif2 + p%d99all1    (loop2+1,loop1,whichairfoil)*redif1 )/(redif1+redif2)
                      m%d99Var    (2) = ( p%d99all2    (loop2+1,loop1+1,whichairfoil)*redif2 + p%d99all2    (loop2+1,loop1,whichairfoil)*redif1 )/(redif1+redif2)
                      m%CfVar     (1) = ( p%Cfall1     (loop2+1,loop1+1,whichairfoil)*redif2 + p%Cfall1     (loop2+1,loop1,whichairfoil)*redif1 )/(redif1+redif2)
                      m%CfVar     (2) = ( p%Cfall2     (loop2+1,loop1+1,whichairfoil)*redif2 + p%Cfall2     (loop2+1,loop1,whichairfoil)*redif1 )/(redif1+redif2)
                      m%EdgeVelVar(1) = ( p%EdgeVelRat1(loop2+1,loop1+1,whichairfoil)*redif2 + p%EdgeVelRat1(loop2+1,loop1,whichairfoil)*redif1 )/(redif1+redif2)
                      m%EdgeVelVar(2) = ( p%EdgeVelRat2(loop2+1,loop1+1,whichairfoil)*redif2 + p%EdgeVelRat2(loop2+1,loop1,whichairfoil)*redif1 )/(redif1+redif2)
                  elseif (AlphaNoise .lt. p%AOAListBL(1)) then
                      CALL WrScr( 'Warning AeroAcoustics Module - Angle of attack (AoA) range is not in the range provided by the user')
                      CALL WrScr( 'Station '// trim(num2lstr(whichairfoil)) )
                      CALL WrScr( 'Airfoil AoA '//trim(num2lstr(AlphaNoise))//'; Using the closest AoA '//trim(num2lstr(p%AOAListBL(1))) )
                      m%dStarVar(1)   = ( p%dstarall1  (1,loop1+1,whichairfoil)*redif2 + p%dstarall1  (1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                      m%dStarVar(2)   = ( p%dstarall2  (1,loop1+1,whichairfoil)*redif2 + p%dstarall2  (1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                      m%d99Var(1)     = ( p%d99all1    (1,loop1+1,whichairfoil)*redif2 + p%d99all1    (1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                      m%d99Var(2)     = ( p%d99all2    (1,loop1+1,whichairfoil)*redif2 + p%d99all2    (1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                      m%CfVar(1)      = ( p%Cfall1     (1,loop1+1,whichairfoil)*redif2 + p%Cfall1     (1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                      m%CfVar(2)      = ( p%Cfall2     (1,loop1+1,whichairfoil)*redif2 + p%Cfall2     (1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                      m%EdgeVelVar(1) = ( p%EdgeVelRat1(1,loop1+1,whichairfoil)*redif2 + p%EdgeVelRat1(1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                      m%EdgeVelVar(2) = ( p%EdgeVelRat2(1,loop1+1,whichairfoil)*redif2 + p%EdgeVelRat2(1,loop1,whichairfoil)*redif1 ) / (redif1+redif2)
                  endif
              endif
          enddo
      endif    
  enddo 
  if (.not. re_flag) then
    call SetErrStat( ErrID_Fatal, 'Warning AeroAcoustics Module - the Reynolds number is not in the range provided by the user. Code stopping.', ErrStat, ErrMsg, RoutineName )
  stop
  endif 
END SUBROUTINE BL_Param_Interp


SUBROUTINE Aero_Tests()
    !--------Laminar Boundary Layer Vortex Shedding Noise----------------------------!
    !CALL LBLVS(AlphaNoise,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
    !    elementspan,m%rTEtoObserve(K,J,I), &
    !    p,m%d99Var(2),m%dstarVar(1),m%dstarVar(2),m%SPLLBL,ErrStat2,errMsg2)
    !--------Turbulent Boundary Layer Trailing Edge Noise----------------------------!
    !CALL TBLTE(3.0d0,0.22860d0,63.920d0,90.0d0,90.0d0,0.5090d0,1.220d0, &
    !    p, m%d99Var(2),m%dstarVar(1),m%dstarVar(2),p%StallStart(J,I),m%SPLP,m%SPLS,m%SPLALPH,m%SPLTBL,ErrStat2,errMsg2 )
    !m%SPLP=0.0_ReKi;m%SPLS=0.0_ReKi;m%SPLTBL=0.0_ReKi;
    !m%EdgeVelVar(1)=1.000d0;m%EdgeVelVar(2)=m%EdgeVelVar(1);
    !m%CfVar(1) = 0.0003785760d0;m%CfVar(2) = 0.001984380d0;m%d99var(1)= 0.01105860d0; m%d99var(2)= 0.007465830d0;m%EdgeVelVar(1)=1.000d0;m%EdgeVelVar(2)=m%EdgeVelVar(1);
    !CALL TBLTE_TNO(0.22860_Reki,63.9200_Reki,90.00_Reki,90.0_Reki,0.5090_Reki,1.220_Reki, &
    !    m%CfVar,m%d99var,m%EdgeVelVar, p, m%SPLP,m%SPLS,m%SPLALPH,m%SPLTBL,ErrStat2 ,errMsg2)
    !--------Blunt Trailing Edge Noise----------------------------------------------!
    !CALL BLUNT(3.0d0,0.22860d0,63.920d0,90.0d0,90.0d0,0.5090d0,1.220d0,&
    !    p%TEThick(J,I),p%TEAngle(J,I),p, m%d99Var(2),m%dstarVar(1),m%dstarVar(2),m%SPLBLUNT,ErrStat2,errMsg2 )
    !--------Tip Noise--------------------------------------------------------------!
    !CALL TIPNOIS(AlphaNoise,p%ALpRAT,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
    !    m%rTEtoObserve(K,J,I), p, m%SPLTIP,ErrStat2,errMsg2)
    !--------Inflow Turbulence Noise ------------------------------------------------!
    !CALL InflowNoise(3.0d0,0.22860d0,63.920d0,90.0d0,90.0d0,0.5090d0,1.220d0, xd%TIVx(J,I),0.050d0,p,m%SPLti,ErrStat2,errMsg2 )
    !CALL FullGuidati(3.0d0,63.920d0,0.22860d0,0.5090d0,1.220d0,90.0d0,90.0d0,xd%MeanVrel(J,I),xd%TIVrel(J,I), &
    !    p,p%BlAFID(J,I),m%SPLTIGui,ErrStat2 )
    !CALL Simple_Guidati(UNoise,0.22860d0,0.120d0,0.020d0,p,m%SPLTIGui,ErrStat2,errMsg2 )
END SUBROUTINE 
END MODULE AeroAcoustics

