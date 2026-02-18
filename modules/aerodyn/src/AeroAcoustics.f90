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
! NOTE: This paper is also known as "BPM Airfoil Self-noise and Prediction paper" in the code documentation.
! NOTE: curve fit equations in the Brooks, Pope, and Marcolini paper use AoA in **degrees** (not radians).
   
!  [2] Moriarty, Guidati, Migliore, Recent Improvement of a Semi-Empirical Aeroacoustic 
!      Prediction Code for Wind Turbines, 2003, NREL/TP-500-34478 (https://docs.nrel.gov/docs/fy04osti/34478.pdf)
!  [3] Lowson, M.V.; Assessment and Prediction of Wind Turbine Noise, Volumes 13-284 of ETSU W. 1993. https://books.google.com/books?id=IgVKGwAACAAJ

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
   
   REAL(ReKi), parameter :: AA_u_min = 0.1_ReKi
   REAL(ReKi), parameter :: AA_EPSILON = 1.E-16 ! EPSILON(AA_EPSILON)
   
   REAL(ReKi), parameter :: RotorRegionAlph_delta = 60.0_ReKi ! degrees : size of bin, must be a number that evenly divides 360 degrees
   REAL(ReKi), parameter :: RotorRegionRad_delta  =  5.0_ReKi ! meters : size of bin along blade span (rotor radius)
   REAL(ReKi), parameter :: RotorRegionTimeSampling  =  5.0_ReKi ! seconds (for Num_total_sampleTI)

   contains    
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine AA_Init( InitInp, u, p, xd, OtherState, y, m, Interval, AFInfo, InitOut, ErrStat, ErrMsg )
   type(AA_InitInputType),       intent(inout) :: InitInp       !< Input data for initialization routine; out because we move allocated array
   type(AA_InputType),           intent(  out) :: u             !< An initial guess for the input; input mesh must be defined
   type(AA_ParameterType),       intent(  out) :: p             !< Parameters
   !type(AA_ContinuousStateType), intent(  out) :: x             !< Initial continuous states
   type(AA_DiscreteStateType),   intent(  out) :: xd            !< Initial discrete states
   !type(AA_ConstraintStateType), intent(  out) :: z             !< Initial guess of the constraint states
   type(AA_OtherStateType),      intent(  out) :: OtherState    !< Initial other states
   type(AA_OutputType),          intent(  out) :: y             !< Initial system outputs (outputs are not calculated;
                                                                !!   only the output mesh is initialized)
   type(AA_MiscVarType),         intent(  out) :: m             !< Initial misc/optimization variables
   real(DbKi),                   intent(in   ) :: interval      !< Coupling interval in seconds: the rate that
                                                                !!   (1) AA_UpdateStates() is called in loose coupling &
                                                                !!   (2) AA_UpdateDiscState() is called in tight coupling.
                                                                !!   Input is the suggested time from the glue code;
                                                                !!   Output is the actual coupling interval that will be used
                                                                !!   by the glue code.
   type(AA_InitOutputType),      intent(  out) :: InitOut       !< Output for initialization routine
   integer(IntKi),               intent(  out) :: errStat       !< Error status of the operation
   character(*),                 intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None
   type(AFI_ParameterType),      intent(in   ) :: AFInfo(:)   !< The airfoil parameter data
!   integer(IntKi),               intent(in   ) :: AFIndx(:,:)
   
   ! Local variables
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
   type(AA_InputFile)                          :: InputFileData ! Data stored in the module's input file
   character(*), parameter                     :: RoutineName = 'AA_Init'
   
   ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""
   ! Initialize the NWTC Subroutine Library
   call NWTC_Init( EchoLibVer=.FALSE. )
   ! Display the module information
   call DispNVD( AA_Ver )

   ! To get rid of a compiler warning.
   !x%DummyContState           = 0.0_SiKi
   !z%DummyConstrState         = 0.0_SiKi

   !bjj: note that we haven't validated p%NumBlades before using it below!
   p%NumBlades = InitInp%NumBlades ! need this before reading the AD input file so that we know how many blade files to read
   p%RootName  = TRIM(InitInp%RootName)//'.'//trim(AA_Nickname)
   
   ! Read the primary AeroAcoustics input file in AeroAcoustics_IO
   call ReadInputFiles( InitInp%InputFile, AFInfo, InputFileData, interval, p%RootName, ErrStat2, ErrMsg2 )
   if (Failed()) return
      
   ! Validate the inputs
   call ValidateInputData(InputFileData, p%NumBlades, ErrStat2, ErrMsg2); if (Failed()) return
    
   ! Validate Initialization Input data ( not found in the AeroAcoustics input file )
   if (InitInp%AirDens <= 0.0)  call SetErrStat ( ErrID_Fatal, 'The air density (AirDens) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InitInp%KinVisc <= 0.0)  call SetErrStat ( ErrID_Fatal, 'The kinesmatic viscosity (KinVisc) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InitInp%SpdSound <= 0.0) call SetErrStat ( ErrID_Fatal, 'The speed of sound (SpdSound) must be greater than zero.', ErrStat, ErrMsg, RoutineName )  
   if (InitInp%NumBlNds < 1) call SetErrStat ( ErrID_Fatal, 'AeroAcoustics requires at least 1 node.', ErrStat, ErrMsg, RoutineName )
   if (Failed()) return

   ! Define parameters
   call SetParameters( InitInp, InputFileData, p, AFInfo, ErrStat2, ErrMsg2 ); if(Failed()) return
   ! Define and initialize inputs 
   call Init_u( u, p, errStat2, errMsg2 ); if(Failed()) return

   ! Initialize states and misc vars
   call Init_MiscVars(m, p, errStat2, errMsg2); if(Failed()) return
   call Init_States(xd, OtherState, p,  errStat2, errMsg2); if(Failed()) return

   ! Define write outputs here (must initialize AFTER Init_MiscVars)
   call Init_y(y, m, p, errStat2, errMsg2); if(Failed()) return

   ! Define initialization output here
   call AA_SetInitOut(p, InitOut, errStat2, errMsg2); if(Failed()) return
   if (AA_OutputToSeparateFile) then
      call AA_InitializeOutputFile(p, InputFileData,InitOut,errStat2, errMsg2); if(Failed()) return
   end if
   call Cleanup() 
      
contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call Cleanup()
    end function Failed

   subroutine Cleanup()
      CALL AA_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2 )
   end subroutine Cleanup
end subroutine AA_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine sets AeroAcoustics parameters for use during the simulation; these variables are not changed after AA_Init.
subroutine SetParameters( InitInp, InputFileData, p, AFInfo, ErrStat, ErrMsg )
    TYPE(AA_InitInputType),       INTENT(INOUT) :: InitInp        !< Input data for initialization routine, out is needed because of copy below
    TYPE(AA_InputFile),           INTENT(INOUT) :: InputFileData  !< Data stored in the module's input file -- intent(out) only for move_alloc statements
    TYPE(AA_ParameterType),       INTENT(INOUT) :: p              !< Parameters
    type(AFI_ParameterType),      intent(in   ) :: AFInfo(:)      !< The airfoil parameter data
    INTEGER(IntKi),               INTENT(  OUT) :: ErrStat        !< Error status of the operation
    CHARACTER(*),                 INTENT(  OUT) :: ErrMsg         ! Error message if ErrStat /= ErrID_None
    ! Local variables
    CHARACTER(ErrMsgLen)    :: ErrMsg2         ! temporary Error message if ErrStat /    = ErrID_None
    INTEGER(IntKi)          :: ErrStat2        ! temporary Error status of the operation
!    INTEGER(IntKi)          :: simcou,coun     ! simple loop  counter
    INTEGER(IntKi)          :: I,J,whichairfoil,K,i1_1,i10_1,i1_2,i10_2,iLE
    character(*), parameter :: RoutineName = 'SetParameters'
    REAL(ReKi)              :: val1,val10,f2,f4, dist1, dist10
    REAL(ReKi)              :: BladeSpanUsedForNoise
    
    ! Initialize variables for this routine
    ErrStat  = ErrID_None
    ErrMsg   = ""
    !!Assign input file data to parameters
    p%DT               = InputFileData%DT_AA         ! seconds
    p%Num_total_sampleTI   = max( NINT(RotorRegionTimeSampling / InputFileData%DT_AA), 1 )
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
    p%outFmt           = "ES15.6E3" 
    p%NumBlNds         = InitInp%NumBlNds
    p%AirDens          = InitInp%AirDens
    p%KinVisc          = InitInp%KinVisc
    p%SpdSound         = InitInp%SpdSound
    p%HubHeight        = InitInp%HubHeight
    p%Lturb            = InputFileData%Lturb
    p%NrObsLoc         = InputFileData%NrObsLoc
    p%FTitle           = InputFileData%FTitle
    p%TI               = InputFileData%TI
    p%avgV             = InputFileData%avgV


    ! Check 1
    IF( (p%ITURB.eq.ITURB_TNO) .or. p%IInflow == IInflow_FullGuidati .OR. p%IInflow == IInflow_SimpleGuidati )then
        ! if tno is on or one of the guidati models is on, check if we have airfoil coordinates
        DO k=1,size(AFInfo) ! if any of the airfoil coordinates are missing change calculation method
            IF( AFInfo(k)%NumCoords .lt. 5 )then
               CALL WrScr( 'Airfoil coordinates are missing: If Full or Simplified Guidati or Bl Calculation is on coordinates are needed ' )
               CALL WrScr( 'Calculation methods enforced as BPM for TBLTE and only Amiet for inflow ' )
               p%ITURB   = ITURB_BPM
               p%IInflow = IInflow_BPM
               exit ! stop checking do loop
            ENDIF
        ENDDO
    ENDIF
    
    ! Check 2
    ! if passed the first check and if tno, turn on boundary layer calculation
    IF( (p%ITURB.eq.ITURB_TNO)) p%X_BLMethod=X_BLMethod_Tables
    
    ! Check 3
    ! if boundary layer is tripped then laminar b.l. vortex shedding mechanism is turned off
    IF( p%ITRIP /= ITRIP_None  ) p%ILAM=ILAM_None

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
    call MOVE_ALLOC(InputFileData%ObsXYZ,p%ObsXYZ)

    ! 
    call MOVE_ALLOC(InitInp%BlAFID,p%BlAFID)
    
    ! Blade Characteristics chord,span,trailing edge angle and thickness,airfoil ID for each segment
    call AllocAry(p%TEThick   ,p%NumBlNds,p%NumBlades,'p%TEThick'   ,ErrStat2,ErrMsg2); if(Failed()) return
    call AllocAry(p%TEAngle   ,p%NumBlNds,p%NumBlades,'p%TEAngle'   ,ErrStat2,ErrMsg2); if(Failed()) return
    call AllocAry(p%StallStart,p%NumBlNds,p%NumBlades,'p%StallStart',ErrStat2,ErrMsg2); if(Failed()) return
    p%StallStart = 0.0_ReKi

     do i=1,p%NumBlades
        do j=1,p%NumBlNds
            whichairfoil = p%BlAFID(j,i)
            p%TEThick(j,i) = InputFileData%BladeProps(whichairfoil)%TEThick
            p%TEAngle(j,i) = InputFileData%BladeProps(whichairfoil)%TEAngle
            
            if(AFInfo(whichairfoil)%NumTabs /=1 ) then
                call SetErrStat(ErrID_Fatal, 'Number of airfoil tables within airfoil file different than 1, which is not supported.', ErrStat2, ErrMsg2, RoutineName )
                if(Failed()) return
            endif
            p%StallStart(j,i)  = AFInfo(whichairfoil)%Table(1)%UA_BL%alpha1*180/PI ! approximate stall angle of attack [deg] (alpha1 in [rad])
        enddo
    enddo

    call AllocAry(p%BlSpn,       p%NumBlNds, p%NumBlades, 'p%BlSpn'    , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(p%BlElemSpn,   p%NumBlNds, p%NumBlades, 'p%BlElemSpn', ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(p%BlChord,     p%NumBlNds, p%NumBlades, 'p%BlChord'  , ErrStat2, ErrMsg2); if(Failed()) return
    call AllocAry(p%AerCent,  2, p%NumBlNds, p%NumBlades, 'p%AerCent'  , ErrStat2, ErrMsg2); if(Failed()) return
    p%BlSpn   = InitInp%BlSpn
    p%BlChord = InitInp%BlChord

    p%startnode = max(1, p%NumBlNds - 1)
    BladeSpanUsedForNoise = p%BlSpn(p%NumBlNds,1)*(1.0 - InputFileData%AA_Bl_Prcntge/100.0)
    do j=p%NumBlNds-1,2,-1
        IF ( p%BlSpn(j,1) .lt. BladeSpanUsedForNoise )THEN
            p%startnode=j
            exit ! exit the loop
        endif
    enddo
    p%startnode = max(min(p%NumBlNds,2),p%startnode)
    
   p%BlElemSpn = 0;
   DO I = 1,p%numBlades
      DO J = p%startnode,p%NumBlNds  ! starts loop from startnode. 
         IF (J < 2) THEN
            p%BlElemSpn(J,I) = p%BlSpn(J,I) !assume this is the innermost node
         ELSEIF (J .EQ. p%NumBlNds) THEN
            p%BlElemSpn(J,I) =   p%BlSpn(J,I)-p%BlSpn(J-1,I)
         ELSE
            p%BlElemSpn(J,I) =   (p%BlSpn(J,I)-p%BlSpn(J-1,I))/2 + (p%BlSpn(J+1,I)-p%BlSpn(J,I))/2 ! this is the average element size around this node, equivalent to (p%BlSpn(J+1,I) - p%BlSpn(J-1,I))/2
         ENDIF
      end do
   end do

    !print*, 'AeroAcoustics Module is using the blade nodes starting from ' ,p%startnode,' Radius in meter ',p%BlSpn(p%startnode,1)
    !AerodYnamic center extraction for each segment 
    do i=1,p%numBlades
        do j=1,p%NumBlNds
            whichairfoil         = p%BlAFID(j,i)  ! just a temporary variable for clear coding
            ! airfoil coordinates read by AeroDyn. First value is the aerodynamic center
            if (AFInfo(whichairfoil)%NumCoords > 0) then
               p%AerCent(1,J,I)  = AFInfo(whichairfoil)%X_Coord(1)  ! assigned here corresponding airfoil.
               p%AerCent(2,J,I)  = AFInfo(whichairfoil)%Y_Coord(1)  ! assigned here corresponding airfoil.
            else
               p%AerCent(1,J,I)  = 0.0_ReKi
               p%AerCent(2,J,I)  = 0.0_ReKi
            end if
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
        call MOVE_ALLOC(InputFileData%AOAListBL,p%AOAListBL)
        call MOVE_ALLOC(InputFileData%ReListBL,p%ReListBL)
        
        ! --- BL data are read from files and just copy what was read from the files
        call MOVE_ALLOC(InputFileData%Suct_DispThick  , p%dstarall1   )
        call MOVE_ALLOC(InputFileData%Pres_DispThick  , p%dstarall2   )
        call MOVE_ALLOC(InputFileData%Suct_BLThick    , p%d99all1     )
        call MOVE_ALLOC(InputFileData%Pres_BLThick    , p%d99all2     )
        call MOVE_ALLOC(InputFileData%Suct_Cf         , p%Cfall1      )
        call MOVE_ALLOC(InputFileData%Pres_Cf         , p%Cfall2      )
        call MOVE_ALLOC(InputFileData%Suct_EdgeVelRat , p%EdgeVelRat1 )
        call MOVE_ALLOC(InputFileData%Pres_EdgeVelRat , p%EdgeVelRat2 )
    endif

    ! If guidati is on, calculate the airfoil thickness at 1% and at 10% chord from input airfoil coordinates
    IF (p%IInflow .EQ. IInflow_FullGuidati) THEN
        call AllocAry(p%AFThickGuida,2,size(AFInfo),  'p%AFThickGuida', errStat2, errMsg2); if(Failed()) return
        p%AFThickGuida=0.0_Reki

        DO k=1,size(AFInfo) ! for each airfoil interpolation 

            ! find index where LE is found
            DO i=3,size(AFInfo(k)%X_Coord)
               IF (AFInfo(k)%X_Coord(i) - AFInfo(k)%X_Coord(i-1) > 0.) THEN
                  iLE = i
                  exit ! end the innermost do loop (i)
               ENDIF
            ENDDO

            ! From LE toward TE
            dist1  = ABS( AFInfo(k)%X_Coord(iLE) - 0.01)
            dist10 = ABS( AFInfo(k)%X_Coord(iLE) - 0.10)
            DO i=iLE+1,size(AFInfo(k)%X_Coord)
                IF (ABS(AFInfo(k)%X_Coord(i) - 0.01) < dist1) THEN
                    i1_1 = i
                    dist1 = ABS(AFInfo(k)%X_Coord(i) - 0.01)
                ENDIF
                IF (ABS(AFInfo(k)%X_Coord(i) - 0.1) < dist10) THEN
                    i10_1 = i
                    dist10 = ABS(AFInfo(k)%X_Coord(i) - 0.1)
                ENDIF
            ENDDO

            ! From TE to LE
            dist1  = 0.99
            dist10 = 0.90
            DO i=1,iLE-1
                IF (ABS(AFInfo(k)%X_Coord(i) - 0.01) < dist1) THEN
                    i1_2 = i
                    dist1 = ABS(AFInfo(k)%X_Coord(i) - 0.01)
                ENDIF
                IF (ABS(AFInfo(k)%X_Coord(i) - 0.1) < dist10) THEN
                    i10_2 = i
                    dist10 = ABS(AFInfo(k)%X_Coord(i) - 0.1)
                ENDIF
            ENDDO

            val1  = AFInfo(k)%Y_Coord(i1_1 ) - AFInfo(k)%Y_Coord(i1_2)
            val10 = AFInfo(k)%Y_Coord(i10_1) - AFInfo(k)%Y_Coord(i10_2)

            p%AFThickGuida(1,k)=val1  ! 1  % chord thickness
            p%AFThickGuida(2,k)=val10 ! 10  % chord thickness
        ENDDO
    ENDIF

   p%NumRotorRegionLimitsAlph = NINT(360./RotorRegionAlph_delta) + 1
   p%NumRotorRegionLimitsRad = CEILING( maxval(p%BlSpn)/RotorRegionRad_delta )+2
   
   call AllocAry( p%RotorRegion_k_minus1, p%NumBlNds, p%NumBlades,  'p%RotorRegion_k_minus1', errStat2, errMsg2); if(Failed()) return
   p%RotorRegion_k_minus1 = 0
   do i=1,p%NumBlades
      do j=1,p%NumBlNds
         p%RotorRegion_k_minus1(j,i) = CEILING( p%BlSpn(j,i) / RotorRegionRad_delta )
         p%RotorRegion_k_minus1(j,i) = MIN( p%NumRotorRegionLimitsRad - 1, MAX( 1, p%RotorRegion_k_minus1(j,i) ) ) !safety
      end do
   enddo

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
subroutine Init_y(y, m, p, errStat, errMsg)
    type(AA_OutputType),           intent(  out)  :: y               !< Module outputs
    type(AA_MiscVarType),          intent(in   )  :: m               !< misc/optimization data
    type(AA_ParameterType),        intent(inout)  :: p               !< Parameters
    integer(IntKi),                intent(  out)  :: errStat         !< Error status of the operation
    character(*),                  intent(  out)  :: errMsg          !< Error message if ErrStat /= ErrID_None

    ! Local variables
    integer(intKi)                               :: ErrStat2          ! temporary Error status
    character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
    character(*), parameter                      :: RoutineName = 'Init_y'

    ! Initialize variables for this routine
    errStat = ErrID_None
    errMsg  = ""

    p%numOutsAll = 0
    
    p%numOutsAll(1) = SIZE(m%DirectiviOutput)
    if (p%NrOutFile > 1) p%numOutsAll(2) = SIZE(m%PtotalFreq) ! SIZE returns total size, including all dimensions of the multi-dimensional array
    if (p%NrOutFile > 2) p%numOutsAll(3) = SIZE(m%SumSpecNoiseSep)
    if (p%NrOutFile > 3) p%numOutsAll(4) = SIZE(m%OASPL)

    if (AA_OutputToSeparateFile) then
        p%numOuts = 0
    else
       p%numOuts = SUM(p%numOutsAll)
    end if

    call AllocAry(y%WriteOutput      , p%numOutsAll(1), 'y%WriteOutput'        , errStat2 , errMsg2); if(Failed()) return
    call AllocAry(y%WriteOutputSep   , p%numOutsAll(3), 'y%WriteOutputSep'     , errStat2 , errMsg2); if(Failed()) return
    call AllocAry(y%WriteOutputForPE , p%numOutsAll(2), 'y%WriteOutputForPE'   , errStat2 , errMsg2); if(Failed()) return
    call AllocAry(y%WriteOutputNodes , p%numOutsAll(4), 'y%WriteOutputSepFreq' , errStat2 , errMsg2); if(Failed()) return

    y%WriteOutput      = 0.0_reki
    y%WriteOutputSep   = 0.0_reki
    y%WriteOutputForPE = 0.0_reki
    y%WriteOutputNodes = 0.0_reki

contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
    end function Failed
end subroutine Init_y
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes (allocates) the misc variables for use during the simulation.
subroutine Init_MiscVars(m, p, errStat, errMsg)
    type(AA_MiscVarType),          intent(inout)  :: m                !< misc/optimization data (not defined in submodules)
    type(AA_ParameterType),        intent(in   )  :: p                !< Parameters
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
    call AllocAry(m%SPLBLUNT    , size(p%FreqList), 'SPLBLUNT'  , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLTIP      , size(p%FreqList), 'SPLTIP'    , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLTI       , size(p%FreqList), 'SPLTI'     , errStat2, errMsg2); if(Failed()) return
    call AllocAry(m%SPLTIGui    , size(p%FreqList), 'SPLTIGui'  , errStat2, errMsg2); if(Failed()) return

    call AllocAry(m%LE_Location,  3, p%NumBlNds, p%numBlades, 'LE_Location', ErrStat2, ErrMsg2); if(Failed()) return
    
   ! arrays for computing WriteOutput values
    call AllocAry(m%DirectiviOutput    , p%NrObsLoc                                                                                            , 'm%DirectiviOutput'    , errStat2 , errMsg2); if(Failed()) return
    call AllocAry(m%SumSpecNoiseSep    , nNoiseMechanism        , size(p%FreqList)          , p%NrObsLoc                                       , 'm%SumSpecNoiseSep'    , errStat2 , errMsg2); if(Failed()) return
    call AllocAry(m%PtotalFreq         , size(p%FreqList)       , p%NrObsLoc                                                                   , 'm%PtotalFreq'         , errStat2 , errMsg2); if(Failed()) return
    call AllocAry(m%OASPL              , p%NrObsLoc             , p%NumBlNds                , p%NumBlades                                      , 'm%OASPL'              , errStat2 , errMsg2); if(Failed()) return

    m%ChordAngleTE = 0.0_ReKi
    m%SpanAngleTE  = 0.0_ReKi
    m%rTEtoObserve = 0.0_ReKi
    m%rLEtoObserve = 0.0_ReKi

    m%SPLTIGui     = 0.0_ReKi
    m%CfVar        = 0.0_ReKi
    m%d99Var       = 0.0_ReKi
    m%dstarVar     = 0.0_ReKi
    m%EdgeVelVar   = 0.0_ReKi
    m%LE_Location  = 0.0_ReKi
contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
    end function Failed
end subroutine Init_MiscVars
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine initializes (allocates) the misc variables for use during the simulation.
subroutine Init_states(xd, OtherState, p, errStat, errMsg)
   type(AA_DiscreteStateType),   intent(inout)  :: xd               !
   type(AA_OtherStateType),      intent(inout)  :: OtherState       !< Initial other states
   type(AA_ParameterType),       intent(in   )  :: p                !< Parameters
   integer(IntKi),               intent(  out)  :: errStat          !< Error status of the operation
   character(*),                 intent(  out)  :: errMsg           !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_states'
    
   ! Initialize variables for this routine
   errStat = ErrID_None
   errMsg  = ""

   call AllocAry(xd%TIVx,       p%NumBlNds, p%numBlades, 'xd%TIVx'      , ErrStat2, ErrMsg2); if(Failed()) return
   xd%TIVx       = 0.0_ReKi

   if (p%TICalcMeth == TICalc_Every) then
      call AllocAry(xd%RegVxStor,         p%Num_total_sampleTI, p%NumRotorRegionLimitsRad-1,p%NumRotorRegionLimitsAlph-1,'xd%Vxst',                  ErrStat2,ErrMsg2); if(Failed()) return
      call AllocAry(OtherState%allregcounter ,                  p%NumRotorRegionLimitsRad-1,p%NumRotorRegionLimitsAlph-1,'OtherState%allregcounter', ErrStat2,ErrMsg2); if(Failed()) return
    
      xd%RegVxStor  = 0.0_reki
      OtherState%allregcounter  = 0
   endif

contains
    logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
    end function Failed
end subroutine Init_states
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AA_UpdateStates( t, n, m, u, p,  xd, OtherState, errStat, errMsg )
   real(DbKi),                     intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),                 intent(in   ) :: n          !< Current simulation time step n = 0,1,...
   type(AA_InputType),             intent(in   ) :: u          !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   TYPE(AA_ParameterType),         INTENT(IN   ) :: p          !< Parameters
   type(AA_DiscreteStateType),     intent(inout) :: xd         !< Input: Discrete states at t;
    type(AA_OtherStateType),       intent(inout) :: OtherState !< Other states (integers)
   type(AA_MiscVarType),           intent(inout) :: m          !< misc/optimization data 
   integer(IntKi),                 intent(  out) :: errStat    !< Error status of the operation
   character(*),                   intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None
   ! local variables
!   integer(intKi)                               :: ErrStat2          ! temporary Error status
!   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'AA_UpdateStates'
!   REAL(ReKi),DIMENSION(p%NumBlNds,p%numBlades) :: TEMPSTD  ! temporary standard deviation variable
   REAL(ReKi)                                   :: InflowNorm,meanInflow,angletemp,abs_le_x   ! temporary standard deviation variable
   integer(intKi)                               :: i,j
   integer(intKi)                               :: k_minus1,rco_minus1

   ErrStat = ErrID_None
   ErrMsg  = ""

   !! Cumulative mean and standard deviation, states are updated as Vx Vy Vz changes at each time step
   !TEMPSTD       =  sqrt( u%Inflow(1,:,:)**2+u%Inflow(2,:,:)**2+u%Inflow(3,:,:)**2 )
   !xd%MeanVxVyVz = (TEMPSTD + xd%MeanVxVyVz*n) / (n+1)  
   !!   xd%VxSq       = TEMPSTD**2 + xd%VxSq
   !!   TEMPSTD     = sqrt(  (xd%VxSq/(n+1)) - (xd%MeanVxVyVz**2)   )
   !!   xd%TIVx  = (TEMPSTD / xd%MeanVxVyVz ) ! check inflow noise input for multiplication with 100 or not

               
   IF( p%TICalcMeth == TICalc_Every ) THEN
       call Calc_LE_Location_Array(p,m,u) ! sets m%LE_Location(:,:,:)
   
       do i=1,p%NumBlades
           do j=1,p%NumBlNds
               abs_le_x=m%LE_Location(3,j,i)-p%hubheight
               
               if (EqualRealNos(abs_le_x, 0.0_ReKi)) then
                  rco_minus1 = 1
               else
                  angletemp = ATAN2(m%LE_Location(2,j,i), abs_le_x) * R2D ! returns angles in the range [-180, 180] degrees
                  if (angletemp<0.) angletemp = angletemp + 360. ! in calculation for rco_minus1 below, we compare angles in the range [0, 360] degrees
                  rco_minus1 = ceiling(angletemp / RotorRegionAlph_delta)
                  rco_minus1 = MIN( p%NumRotorRegionLimitsAlph-1, MAX(1, rco_minus1) ) ! safety
               end if

               k_minus1 = p%RotorRegion_k_minus1(j,i)
               
               OtherState%allregcounter(k_minus1,rco_minus1) = OtherState%allregcounter(k_minus1,rco_minus1) + 1    ! increase the sample amount in that specific bin
               
               InflowNorm = TwoNorm( u%Inflow(:,j,i) )
               !note: p%Num_total_sampleTI = size(xd%RegVxStor,1)
               ! with storage region dependent moving average and TI
               IF  ( OtherState%allregcounter(k_minus1,rco_minus1) <= p%Num_total_sampleTI ) THEN
                   xd%RegVxStor(OtherState%allregcounter(k_minus1,rco_minus1),k_minus1,rco_minus1) = InflowNorm
                   xd%TIVx(j,i) = 0
               ELSE
                   xd%RegVxStor( mod( OtherState%allregcounter(k_minus1,rco_minus1), p%Num_total_sampleTI )+1, k_minus1, rco_minus1)=InflowNorm
                   meanInflow = SUM( xd%RegVxStor(:,k_minus1,rco_minus1) ) /p%Num_total_sampleTI

                   if ( EqualRealNos(meanInflow,0.0_ReKi)) then
                      xd%TIVx(j,i) = 0.0_ReKi
                   else
                      xd%TIVx(j,i) = SQRT( SUM((xd%RegVxStor(:,k_minus1,rco_minus1)-meanInflow)**2) /  p%Num_total_sampleTI ) ! only the fluctuation  (this is the population standard deviation, not TI)
                      xd%TIVx(j,i) = xd%TIVx(j,i) / meanInflow ! this is TI as a fraction (std(U)/mean(U))
                   end if
               ENDIF
           enddo
       enddo
       
   ELSE! interpolate from the user given ti values
       do i=1,p%NumBlades
           do j=1,p%NumBlNds
               ! We scale the incident turbulence intensity by the ratio of average to incident wind speed
                ! The scaled TI is used by the Amiet model
                xd%TIVx(j,i)=p%TI * p%avgV/u%Vrel(J,I) 
           enddo
       enddo
   endif
   
end subroutine AA_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine AA_End( u, p, xd, OtherState, y, m, ErrStat, ErrMsg )
    TYPE(AA_InputType),           INTENT(INOUT)  :: u           !< System inputs
    TYPE(AA_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
    !TYPE(AA_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
    TYPE(AA_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
    !TYPE(AA_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
    TYPE(AA_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
    TYPE(AA_OutputType),          INTENT(INOUT)  :: y           !< System outputs
    TYPE(AA_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables
    INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
    CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      
    integer(IntKi)                               :: j

      ! Initialize ErrStat
    ErrStat = ErrID_None
    ErrMsg  = ""
    
      
   do j=1,SIZE(p%unOutFile)
      if (p%unOutFile(j) > 0) then
         close(p%unOutFile(j))
         p%unOutFile(j) = -1
      end if
   end do
       
       
    !! Destroy the input data:
    !CALL AA_DestroyInput( u, ErrStat, ErrMsg )
    !
    !! Destroy the parameter data:
    !CALL AA_DestroyParam( p, ErrStat, ErrMsg )
    !
    !! Destroy the state data:
    !CALL AA_DestroyContState(   x,           ErrStat, ErrMsg )
    !CALL AA_DestroyDiscState(   xd,          ErrStat, ErrMsg )
    !CALL AA_DestroyConstrState( z,           ErrStat, ErrMsg )
    !CALL AA_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
    !CALL AA_DestroyMisc(        m,           ErrStat, ErrMsg ) 
    !! Destroy the output data:
    !CALL AA_DestroyOutput( y, ErrStat, ErrMsg )

END SUBROUTINE AA_End

!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine AA_CalcOutput( t, u, p, xd, OtherState, y, m, ErrStat, ErrMsg)
    ! NOTE: no matter how many channels are selected for output, all of the outputs are calcalated
    ! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
    ! placed in the y%WriteOutput(:) array.
    !..................................................................................................................................
    REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
    TYPE(AA_InputType),           INTENT(IN   )  :: u           !< Inputs at Time t
    TYPE(AA_ParameterType),       INTENT(IN   )  :: p           !< Parameters
    !TYPE(AA_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
    TYPE(AA_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
    !TYPE(AA_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
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
    
    IF (t >= p%AAStart) THEN
       
        IF (.NOT. AA_OutputToSeparateFile .or. mod(t + 1E-10,p%DT) .lt. 1E-6) THEN !bjj: should check NINT(t/p%DT)?
            call CalcObserve(p,m,u,errStat2, errMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (ErrStat >= AbortErrLev) return
           
            call CalcAeroAcousticsOutput(u,p,m,xd,errStat2,errMsg2)
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (ErrStat >= AbortErrLev) return

            call Calc_WriteOutput( p, m, y,  ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (ErrStat >= AbortErrLev) return

            if (AA_OutputToSeparateFile) then
               call AA_WriteOutputLine(y, t, p, ErrStat2, ErrMsg2)
               call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName); if (ErrStat >= AbortErrLev) return
            end if
        ENDIF

    ENDIF
    
end subroutine AA_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
REAL(ReKi) FUNCTION Log10AA(X) RESULT(F)
   REAL(ReKi),INTENT(IN) :: X
   
    F = LOG10( MAX(AA_EPSILON, X) )
    
END FUNCTION Log10AA
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE Calc_LE_Location_Array(p,m,u)
    TYPE(AA_ParameterType),              intent(in   ) :: p       !< Parameters
    TYPE(AA_InputType),                  intent(in   ) :: u       !< NN Inputs at Time
    TYPE(AA_MiscVarType),                intent(inout) :: m       !< misc/optimization data (not defined in submodules)
    ! Local variables.
    INTEGER(intKi) :: I                           ! I A generic index for DO loops.
    INTEGER(intKi) :: J                           ! J A generic index for DO loops.


    ! Loop through the blades
    DO I = 1,p%numBlades
        ! Loop through the nodes along blade span
        DO J = 1,p%NumBlNds
            ! Transpose the rotational vector GlobalToLocal to obtain the rotation LocalToGlobal
            ! LocalToGlobal  = TRANSPOSE(u%RotGtoL(:,:,J,I))

            ! Rotate the coordinates of leading and trailing edge from the local reference system to the global. Then add the coordinates of the aerodynamic center in the global coordinate system
            ! The global coordinate system is located on the ground, has x pointing downwind, y pointing laterally, and z pointing vertically upwards

            !m%LE_Location(:,J,I) = RLEObservereal = MATMUL(LocalToGlobal, p%AFLeCo(:,J,I)) + u%AeroCent_G(:,J,I)
            m%LE_Location(:,J,I) = MATMUL(p%AFLeCo(:,J,I), u%RotGtoL(:,:,J,I) ) + u%AeroCent_G(:,J,I) ! = because this is a matrix times a vector, we can do the transpose of the actual equation: MATMUL(TRANSPOSE(u%RotGtoL(:,:,J,I)), p%AFLeCo(:,J,I)) + u%AeroCent_G(:,J,I)
            
        ENDDO  !J, blade nodes
    ENDDO  !I , number of blades
    
END SUBROUTINE Calc_LE_Location_Array
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE CalcObserve(p,m,u,errStat,errMsg)
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
     
    call Calc_LE_Location_Array(p,m,u) ! sets m%LE_Location(:,:,:)
    
   ! Loop through the blades
   DO I = 1,p%numBlades
      ! Loop through the nodes along blade span
      DO J = 1,p%NumBlNds
         ! Rotate the coordinates of leading and trailing edge from the local reference system to the global. Then add the coordinates of the aerodynamic center in the global coordinate system
         ! The global coordinate system is located on the ground, has x pointing downwind, y pointing laterally, and z pointing vertically upwards
         RTEObservereal = MATMUL(p%AFTeCo(:,J,I), u%RotGtoL(:,:,J,I)) + u%AeroCent_G(:,J,I) ! Note that with the vector math, this is equivalent to MATMUL(TRANSPOSE(p%RotGtoL(:,:,J,I)), p%AFTeCo(:,J,I)) + u%AeroCent_G(:,J,I)
            
         ! Loop through the observers
         DO K = 1,p%NrObsLoc
            
            RTEObserveG=p%ObsXYZ(:,K)-RTEObservereal            ! Calculate the position of the observer K in a reference system located at the trailing edge and oriented as the global reference system
            RLEObserveG=p%ObsXYZ(:,K)-m%LE_Location(:,J,I)      ! Calculate the position of the observer K in a reference system located at the leading edge and oriented as the global reference system
            ! Rotate back the two reference systems from global to local. 
            RTEObserve = MATMUL(u%RotGtoL(:,:,J,I), RTEObserveG)
            RLEObserve = MATMUL(u%RotGtoL(:,:,J,I), RLEObserveG)

            ! Calculate absolute distance between node and observer
            m%rTEtoObserve(K,J,I) = max(AA_Epsilon, TwoNorm(RTEObserve) )
            m%rLEtoObserve(K,J,I) = max(AA_Epsilon, TwoNorm(RLEObserve) )

            ! Calculate time of noise propagation to observer
            !timeTE = m%rTEtoObserve(K,J,I) / p%SpdSound
            !timeLE = m%rLEtoObserve(K,J,I) / p%SpdSound
                        
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
      ENDDO  !J, blade nodes
   ENDDO  !I , number of blades

END SUBROUTINE CalcObserve
!----------------------------------------------------------------------------------------------------------------------------------!
SUBROUTINE CalcAeroAcousticsOutput(u,p,m,xd,errStat,errMsg)
    TYPE(AA_InputType),                     INTENT(IN   )   :: u       !< Inputs at Time t
    TYPE(AA_ParameterType),                 INTENT(IN   )   :: p       !< Parameters
    TYPE(AA_MiscVarType),                   INTENT(INOUT)   :: m       !< misc/optimization data (not defined in submodules)
    TYPE(AA_DiscreteStateType),             INTENT(IN   )   :: xd      !< discrete state type
    integer(IntKi),                         INTENT(  OUT)   :: errStat !< Error status of the operation
    character(*),                           INTENT(  OUT)   :: errMsg  !< Error message if ErrStat /= ErrID_None
    ! Local variables.
    integer(intKi)                :: III                              ! III A generic index for DO loops (frequency)
    integer(intKi)                :: I                                ! I   A generic index for DO loops (blade)
    integer(intKi)                :: J                                ! J   A generic index for DO loops (blade node)
    integer(intKi)                :: K                                ! K   A generic index for DO loops (NrObsLoc)
    integer(intKi)                :: oi                               ! oi   A generic index for DO loops (NoiseMechanism)
    REAL(ReKi)                    :: AlphaNoise                       
    REAL(ReKi)                    :: AlphaNoise_Deg                   ! 
    REAL(ReKi)                    :: UNoise                           ! 

    real(ReKi)                    ::  Ptotal
    character(*), parameter       :: RoutineName = 'CalcAeroAcousticsOutput'
    
    ErrStat = ErrID_None
    ErrMsg  = ""

    !------------------- Initialize arrays with zeros -------------------------!
   ! values for WriteOutput
   m%OASPL = 0.0_Reki
   m%DirectiviOutput = 0.0_Reki
   m%SumSpecNoiseSep = 0.0_Reki
   !----------------
   m%SPLLBL=0.0_Reki
   m%SPLP=0.0_Reki
   m%SPLS=0.0_Reki
   m%SPLALPH=0.0_Reki
   m%SPLBLUNT=0.0_Reki
   m%SPLTIP=0.0_Reki
   m%SPLti=0.0_Reki


   DO I = 1,p%numBlades
      DO J = p%startnode,p%NumBlNds  ! starts loop from startnode. 
         !------------------------------!!------------------------------!!------------------------------!!------------------------------!

         Unoise =  u%Vrel(J,I) 
         IF (abs(Unoise) < AA_u_min) then
            Unoise = SIGN(AA_u_min, Unoise)
         ENDIF
            
         AlphaNoise= u%AoANoise(J,I)
         call MPi2Pi(AlphaNoise) ! make sure this is in an appropriate range [-pi,pi]
         AlphaNoise_Deg = AlphaNoise * R2D_D ! convert to degrees since that is how this code is set up.

         !--------Read in Boundary Layer Data-------------------------!
         IF (p%X_BLMethod .EQ. X_BLMethod_Tables) THEN
            call BL_Param_Interp(p, m, Unoise, AlphaNoise_Deg, p%BlChord(J,I), p%BlAFID(J,I))

            m%d99Var     = m%d99Var*p%BlChord(J,I)
            m%dstarVar   = m%dstarVar*p%BlChord(J,I)
         ENDIF

            !------------------------------!!------------------------------!!------------------------------!!------------------------------!
            !------------------------------!!------------------------------!!------------------------------!!------------------------------!
            !------------------------------!!------------------------------!!------------------------------!!------------------------------!
         DO K = 1,p%NrObsLoc
            Ptotal = 0.0_ReKi         ! Total Sound Pressure - All (7) mechanisms, All Frequencies

            !--------Laminar Boundary Layer Vortex Shedding Noise----------------------------!
            IF ( (p%ILAM .EQ. ILAM_BPM) .AND. (p%ITRIP .EQ. ITRIP_None) )    THEN
               CALL LBLVS(AlphaNoise_Deg,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
                  p%BlElemSpn(J,I),m%rTEtoObserve(K,J,I), p,m%d99Var(2),m%dstarVar(1),m%dstarVar(2),m%SPLLBL,p%StallStart(J,I))
               
               call TotalContributionFromType(m%SPLLBL,Ptotal,NoiseMech=1)
            ENDIF
            
            !--------Turbulent Boundary Layer Trailing Edge Noise----------------------------!
            IF ( p%ITURB /= ITURB_None )   THEN
               !returns  m%SPLP, m%SPLS, m%SPLALPH
               CALL TBLTE(AlphaNoise_Deg,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
                      p%BlElemSpn(J,I),m%rTEtoObserve(K,J,I), p, m%d99Var(2),m%dstarVar(1),m%dstarVar(2),p%StallStart(J,I), &
                      m%SPLP,m%SPLS,m%SPLALPH )
                    
               IF (p%ITURB .EQ. ITURB_TNO)  THEN
                  m%EdgeVelVar=1.0_ReKi
                  !returns m%SPLP, m%SPLS from TBLTE
                  CALL TBLTE_TNO(UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
                        p%BlElemSpn(J,I),m%rTEtoObserve(K,J,I),m%CfVar,m%d99var,m%EdgeVelVar ,p, &
                        m%SPLP,m%SPLS)
               ENDIF
               
               ! If flag for TBL is ON, compute Pressure, Suction, and AoA contributions
               call TotalContributionFromType(m%SPLP,Ptotal,NoiseMech=2)
               call TotalContributionFromType(m%SPLS,Ptotal,NoiseMech=3)
               call TotalContributionFromType(m%SPLALPH,Ptotal,NoiseMech=4)
            ENDIF
            
                
            !--------Blunt Trailing Edge Noise----------------------------------------------!
            IF ( p%IBLUNT == IBLUNT_BPM )   THEN   ! calculate m%SPLBLUNT(1:nFreq)
               CALL BLUNT(AlphaNoise_Deg,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
               p%BlElemSpn(J,I),m%rTEtoObserve(K,J,I),p%TEThick(J,I),p%TEAngle(J,I), &
               p, m%d99Var(2),m%dstarVar(1),m%dstarVar(2),m%SPLBLUNT,p%StallStart(J,I) )
               
               call TotalContributionFromType(m%SPLBLUNT,Ptotal,NoiseMech=5)
            ENDIF
            
            
            !--------Tip Noise--------------------------------------------------------------!
            IF ( (p%ITIP == ITIP_ON) .AND. (J .EQ. p%NumBlNds) ) THEN ! calculate m%SPLTIP(1:nFreq)
               CALL TIPNOIS(AlphaNoise_Deg,p%ALpRAT,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
                  m%rTEtoObserve(K,J,I), p, m%SPLTIP)
               
               ! If flag for Tip is ON and the current blade node (J) is the last node (tip), compute Tip contribution
               call TotalContributionFromType(m%SPLTIP,Ptotal,NoiseMech=6)
            ENDIF

                
            !--------Inflow Turbulence Noise ------------------------------------------------!
            ! important checks to be done inflow tubulence inputs
            IF (p%IInflow /= IInflow_None) then

               ! Amiet's Inflow Noise Model is Calculated as long as InflowNoise is On
               CALL InflowNoise(AlphaNoise,p%BlChord(J,I),Unoise,m%ChordAngleLE(K,J,I),m%SpanAngleLE(K,J,I),&
                  p%BlElemSpn(J,I),m%rLEtoObserve(K,J,I),xd%TIVx(J,I),p,m%SPLti )
                    
               ! If Guidati model (simplified or full version) is also on then the 'SPL correction' to Amiet's model will be added 
               IF ( p%IInflow .EQ. IInflow_FullGuidati )   THEN      
                  CALL Simple_Guidati(UNoise,p%BlChord(J,I),p%AFThickGuida(2,p%BlAFID(J,I)), p%AFThickGuida(1,p%BlAFID(J,I)),p,m%SPLTIGui )
                  m%SPLti = m%SPLti+m%SPLTIGui + 10. ! +10 is fudge factor to match NLR data
               ELSEIF ( p%IInflow .EQ. IInflow_SimpleGuidati )   THEN
                  call setErrStat(ErrID_Fatal,'Full Guidati removed',ErrStat, ErrMsg,RoutineName)
                  return
               ENDIF
                    
               call TotalContributionFromType(m%SPLti,Ptotal,NoiseMech=7) ! compute Turbulent Inflow contribution
            ENDIF
            !m%DirectiviOutput(K) = Ptotal + m%DirectiviOutput(K)   ! Assigns Overall Pressure to Appropriate Observer for Directivity
            
            m%OASPL(K,J,I) = Ptotal + m%OASPL(K,J,I)               ! Assigns Overall Pressure to Appropriate Observer/Blade/Node for Directivity
         ENDDO ! Loop on observers (K)
            
      ENDDO ! Loop on blade nodes (J)
   ENDDO ! Loop on blades (I)

    ! If any Output file is wanted, convert DirectiviOutput from Directivity Factor to Directivity Index
    ! Ref: Fundamentals of Acoustics by Colin Hansen (1951)
    
    ! Since these will all be converted via LOG10, they will produce an error if .EQ. 0., Set .EQ. to 1 instead (LOG10(1)=0)
   DO K = 1,p%NrObsLoc 
      m%DirectiviOutput(K) = SUM(m%SumSpecNoiseSep(:,:,K))
      
      IF (m%DirectiviOutput(K)  .NE. 0.)      m%DirectiviOutput(K) = 10.*LOG10(m%DirectiviOutput(K))        !! DirectiviOutput is used as total observer OASPL for Output File 1
   ENDDO ! Loop on observers
   
   IF  (p%NrOutFile .gt. 1) THEN

      ! Procedure for Output file 2
      DO K = 1,p%NrObsLoc
         DO III=1,size(p%FreqList)
            m%PtotalFreq(III,K) = SUM( m%SumSpecNoiseSep(:,III,K) )

            IF (m%PtotalFreq(III,K) .NE. 0.)  m%PtotalFreq(III,K)    = 10.*LOG10(m%PtotalFreq(III,K))               ! P to SPL conversion
         ENDDO
      ENDDO

      ! Procedure for Output file 3; If 3rd Output file is needed, convert P to SPL (skip values = 0).
      IF  (p%NrOutFile .gt. 2) THEN 
         DO K = 1,p%NrObsLoc
            DO III = 1,size(p%FreqList)
               DO oi = 1,nNoiseMechanism
                  IF (m%SumSpecNoiseSep(oi,III,K)  .NE. 0.) m%SumSpecNoiseSep(oi,III,K) = 10.*LOG10(m%SumSpecNoiseSep(oi,III,K))      ! P to SPL Conversion
               ENDDO
            ENDDO
         ENDDO
            
      ! Procedure for Output file 3; If 4th Output file is needed, convert P to SPL (skip values = 0).
        IF  (p%NrOutFile .gt. 3) THEN 
            DO I = 1,p%numBlades
               DO J = 1,p%NumBlNds
                  DO K = 1,p%NrObsLoc 
                     IF (m%OASPL(K,J,I)  .NE. 0.) m%OASPL(K,J,I) = 10.*LOG10(m%OASPL(K,J,I))
                  ENDDO
               ENDDO
            ENDDO
         END IF ! file 4

      ENDIF ! file 3

   END IF ! file 2

   
contains

   subroutine TotalContributionFromType(SPL,Ptotal,NoiseMech)
      REAL(ReKi),     intent(inout) :: SPL(:)
      INTEGER(IntKi), intent(in   ) :: NoiseMech ! number of noise mechanism (index into SumSpecNoiseSep)
      REAL(ReKi),     intent(inout) :: Ptotal
      REAL(ReKi)                    :: Pt
      REAL(ReKi)                    :: P_SumAllFreq

      IF (p%AweightFlag) THEN
         SPL = SPL + p%Aweight                     ! A-weighting for all frequencies
      ENDIF

      P_SumAllFreq = 0.0_ReKi
                  
      do III=1,size(p%FreqList)   ! Loops through each 1/3rd octave center frequency 

         Pt = 10.0_ReKi**(SPL(III)/10.0_ReKi)                                           ! SPL to P Conversion for III Frequency
                        
         P_SumAllFreq = P_SumAllFreq + Pt                                               ! Sum for Running Total
         m%SumSpecNoiseSep(NoiseMech,III,K) = m%SumSpecNoiseSep(NoiseMech,III,K) + Pt   ! Running sum of observer and frequency dependent sound pressure

      end do
      Ptotal = Ptotal + P_SumAllFreq
   end subroutine
   
END SUBROUTINE CalcAeroAcousticsOutput
!==================================================================================================================================!
SUBROUTINE LBLVS(ALPSTAR,C,U,THETA,PHI,L,R,p,d99Var2,dstarVar1,dstarVar2,SPLLAM,StallVal)
    REAL(ReKi),                                 INTENT(IN   ) :: ALPSTAR        ! AOA, deg
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
    
    !compute reynolds number and mach number
    M          = U  / p%SpdSound        ! MACH NUMBER
    RC         = U  * C/p%KinVisc       ! REYNOLDS NUMBER BASED ON  CHORD
    
    ! compute boundary layer thicknesses
    IF (p%X_BLMethod .eq. X_BLMethod_Tables) THEN
        DELTAP = d99Var2
        DSTRS  = dstarVar1
        DSTRP  = dstarVar2
    ELSE
        CALL THICK(C,RC,ALPSTAR,p,DELTAP,DSTRS,DSTRP,StallVal)
    ENDIF
    
    ! compute directivity function
    DBARH = DIRECTH_TE(M,THETA,PHI)

    IF (DBARH <= 0) THEN
        SPLLAM = 0.
        RETURN
    ENDIF
    
    ! compute reference strouhal number                                 ! Eq 55 from BPM Airfoil Self-noise and Prediction paper
    if (RC .LE. 1.3E+05) then
       ST1PRIM = .18
    elseif (RC.LE.4.0E+05) then
       ST1PRIM=.001756*RC**.3931
    else
       ST1PRIM = .28
    end if
    STPKPRM  = 10.**(-.04*ALPSTAR) * ST1PRIM                            ! Eq 56 from BPM Airfoil Self-noise and Prediction paper

    ! compute reference reynolds number                                 ! Eq 59 from BPM Airfoil Self-noise and Prediction paper
    IF (ALPSTAR .LE. 3.0) then
       RC0=10.**(.215*ALPSTAR+4.978)
    else
       RC0=10.**(.120*ALPSTAR+5.263)
    end if
    
    ! compute peak scaled spectrum level
    D   = RC / RC0                                                      ! Used in Eq 58 from BPM Airfoil Self-noise and Prediction paper
    if (D .LE. .3237) then
       G2 =77.852*LOG10AA(D)+15.328       ! Begin Eq 58 from BPM Airfoil Self-noise and Prediction paper
    elseif (D .LE. .5689) then
       G2 = 65.188*LOG10(D) + 9.125
    elseif (D .LE. 1.7579) then
       G2 = -114.052 * LOG10(D)**2
    elseif (D .LE. 3.0889) then
       G2 = -65.188*LOG10(D)+9.125
    else 
       G2 =-77.852*LOG10(D)+15.328
    end if
    
    ! compute angle-dependent level for shape curve
    G3      = 171.04 - 3.03 * ALPSTAR                                    ! Eq 60 from BPM Airfoil Self-noise and Prediction paper
    SCALE   = 10. * Log10AA(DELTAP*M**5*DBARH*L/R**2)                    ! From Eq 53 from BPM Airfoil Self-noise and Prediction paper
    
    ! Compute scaled sound pressure levels for each strouhal number
    DO I=1,SIZE(p%FreqList)
        STPRIM  = p%FreqList(I) * DELTAP / U                            ! Eq 54 from BPM Airfoil Self-noise and Prediction paper
        E       = STPRIM / STPKPRM                                      ! Used in Eq 57 from BPM Airfoil Self-noise and Prediction paper
        IF (E .LE. .5974) then
           G1 = 39.8*LOG10AA(E)-11.12                   ! Begin Eq 57 from BPM Airfoil Self-noise and Prediction paper   
        ELSEIF(E .LE. .8545) then
           G1 = 98.409 * LOG10(E) + 2.0
        ELSEIF (E .LE. 1.17) then
           G1 = -5.076+SQRT(2.484-506.25*(LOG10(E))**2)
        ELSEIF (E .LE. 1.674) then
           G1 = -98.409 * LOG10(E) + 2.0
        ELSE
           G1 = -39.80*LOG10(E)-11.12
        END IF
        SPLLAM(I) = G1 + G2 + G3 + SCALE                                      ! Eq 53 from BPM Airfoil Self-noise and Prediction paper
    ENDDO
END SUBROUTINE LBLVS
!==================================================================================================================================!
SUBROUTINE TBLTE(ALPSTAR,C,U,THETA,PHI,L,R,p,d99Var2,dstarVar1,dstarVar2,StallVal,SPLP,SPLS,SPLALPH)
    REAL(ReKi),                             INTENT(IN   )  :: ALPSTAR        ! AOA(deg)
    REAL(ReKi),                             INTENT(IN   )  :: C              ! Chord Length           (m)
    REAL(ReKi),                             INTENT(IN   )  :: U              ! Unoise(m/s)
    REAL(ReKi),                             INTENT(IN   )  :: THETA          ! DIRECTIVITY ANGLE      (deg)
    REAL(ReKi),                             INTENT(IN   )  :: PHI            ! DIRECTIVITY ANGLE      (deg) 
    REAL(ReKi),                             INTENT(IN   )  :: L              ! SPAN(m)
    REAL(ReKi),                             INTENT(IN   )  :: R              ! SOURCE TO OBSERVER DISTANCE (m)
    TYPE(AA_ParameterType),                 INTENT(IN   )  :: p              ! Noise Module Parameters

!    REAL(ReKi)                               :: L              ! SPAN(m)
!    REAL(ReKi)                               :: R              ! SOURCE TO OBSERVER DISTANCE (m)

    REAL(ReKi),                             INTENT(IN   )  :: d99Var2        !  
    REAL(ReKi),                             INTENT(IN   )  :: dstarVar1              !  
    REAL(ReKi),                             INTENT(IN   )  :: dstarVar2              !  
    REAL(ReKi),                             INTENT(IN   )  :: StallVal              !  

    REAL(ReKi),DIMENSION(size(p%FreqList)),  INTENT(  OUT)  :: SPLP           ! SOUND PRESSURE LEVEL DUE TO PRESSURE SIDE OF AIRFOIL (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),  INTENT(  OUT)  :: SPLS           ! SOUND PRESSURE LEVEL DUE TO SUCTION SIDE OF AIRFOIL  (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),  INTENT(  OUT)  :: SPLALPH        ! SOUND PRESSURE LEVEL DUE TO ANGLE OF ATTACK CONTRIBUTION (db)
    
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
    !real(ReKi)   :: P1         ! PRESSURE SIDE PRESSURE (NT/M2)
    !real(ReKi)   :: P2         ! SUCTION SIDE PRESSURE                      (NT/M2)
    !real(ReKi)   :: P4         ! PRESSURE FROM ANGLE OF ATTACK CONTRIBUTION (NT/M2)
    real(ReKi)   :: M          ! MACH NUMBER
    real(ReKi)   :: RC         ! REYNOLDS NUMBER BASED ON  CHORD
    real(ReKi)   :: DELTAP     ! PRESSURE SIDE BOUNDARY LAYER THICKNESS METERS
    real(ReKi)   :: XCHECK     ! USED TO CHECK FOR ANGLE OF ATTACK CONTRIBUTION     
    real(ReKi)   :: DBARH      ! HIGH FREQUENCY DIRECTIVITY             ---
    real(ReKi)   :: DBARL      ! LOW FREQUENCY DIRECTIVITY              ---
    
    integer(intKi)      :: I          ! I A generic index for DO loops.

    LOGICAL     :: SWITCH  !!LOGICAL FOR COMPUTATION OF ANGLE OF ATTACK CONTRIBUTION  

    ! Compute reynolds number and mach number
    M          = U  / p%SpdSound
    RC         = U  * C/p%KinVisc
    
    ! Compute boundary layer thicknesses
    IF (p%X_BLMethod .eq. X_BLMethod_Tables) THEN
        DELTAP = d99Var2
        DSTRS  = dstarVar1
        DSTRP  = dstarVar2
    ELSE
        CALL THICK(C,RC,ALPSTAR,p,DELTAP,DSTRS,DSTRP,StallVal)
    ENDIF
    
    ! Compute directivity function
    DBARL = DIRECTL(M,THETA,PHI)
    DBARH = DIRECTH_TE(M,THETA,PHI)
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
    IF  (ALPSTAR .LE. 1.333) then
       ST2 = ST1
    elseif (ALPSTAR .LE. StallVal) then
       ST2 = ST1*10.**(.0054*(ALPSTAR-1.333)**2)
    else
       ST2 = 4.72 * ST1
    end if
    
    ST1PRIM = (ST1+ST2)/2.                                                             ! Eq 33 from BPM Airfoil Self-noise and Prediction paper
    A0 = A0COMP(RC)      ! compute -20 dB dropout   (returns A0)
    A02 = A0COMP(3.0_ReKi*RC)   ! compute -20 dB dropout for AoA > AoA_0   (returns A02)
    ! Evaluate minimum and maximum 'a' curves at a0
    AMINA0 = AMIN(A0)
    AMAXA0 = AMAX(A0)
    AMINA02 = AMIN(A02)
    AMAXA02 = AMAX(A02)
    ! Compute 'a' max/min ratio                                                        ! Eq 39 from BPM Airfoil Self-noise and Prediction paper   
    ARA0  = (20. + AMINA0) / (AMINA0 - AMAXA0)
    ARA02 = (20. + AMINA02)/ (AMINA02- AMAXA02)
    
    ! Compute b0 to be used in 'b' curve calculations                                  ! Eq 44 from BPM Airfoil Self-noise and Prediction paper
    IF (RC .LT. 9.52E+04) then
       B0 = .30
    elseif (RC .LT. 8.57E+05) then
       B0 = (-4.48E-13)*(RC-8.57E+05)**2 + .56
    else
       B0 = .56
    end if
    
    ! Evaluate minimum and maximum 'b' curves at b0
    BMINB0 = BMIN(B0)
    BMAXB0 = BMAX(B0)
    ! Compute 'b' max/min ratio
    BRB0  = (20. + BMINB0) / (BMINB0 - BMAXB0)

    ! For each center frequency, compute an 'a' prediction for the pressure side
    STPEAK = ST1
    IF (RC .LT. 2.47E+05) then
       K1 = -4.31 * LOG10AA(RC) + 156.3        ! Begin Eq 47 from BPM Airfoil Self-noise and Prediction paper         
    elseif (RC .LE. 8.0E+05) then
       K1 = -9.0 * LOG10(RC) + 181.6
    else
       K1 = 128.5
    end if
    
    IF (RDSTRP .LE. 5000.) then
       DELK1 = -ALPSTAR*(5.29-1.43*LOG10AA(RDSTRP))                   ! Begin Eq 48 from BPM Airfoil Self-noise and Prediction paper
    else
       DELK1 = 0.0
    end if
    
    GAMMA   = 27.094 * M +  3.31                                                       ! Begin Eq 49 from BPM Airfoil Self-noise and Prediction paper
    BETA    = 72.650 * M + 10.74
    GAMMA0  = 23.430 * M +  4.651
    BETA0   =-34.190 * M - 13.820                                                      ! end

    if (ALPSTAR .LE. (GAMMA0-GAMMA)) then
       K2 = -1000.0                                      ! Begin Eq 49 from BPM Airfoil Self-noise and Prediction paper
    else if (ALPSTAR.LE.(GAMMA0+GAMMA)) then
       K2=SQRT(BETA**2-(BETA/GAMMA)**2*(ALPSTAR-GAMMA0)**2)+BETA0
    else
       K2 = -12.0
    end if
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
        A      = LOG10AA( STP / STPEAK )                                           ! Eq 37 from BPM Airfoil Self-noise and Prediction paper
        AMINA = AMIN(A)
        AMAXA = AMAX(A)
        AA     = AMINA + ARA0 * (AMAXA - AMINA)                            ! Eq 40 from BPM Airfoil Self-noise and Prediction paper

        SPLP(I)=AA+K1-3.+10.*LOG10AA(DSTRP*M**5*DBARH*L/R**2)+DELK1        ! Eq 25 from BPM Airfoil Self-noise and Prediction paper
        STS = p%FreqList(I) * DSTRS / U                                    ! Eq 31 from BPM Airfoil Self-noise and Prediction paper

        IF (.NOT. SWITCH) THEN
            A      = LOG10AA( STS / ST1PRIM )
            AMINA = AMIN(A)
            AMAXA = AMAX(A)
            AA = AMINA + ARA0 * (AMAXA - AMINA)
            SPLS(I) = AA+K1-3.+10.*LOG10AA(DSTRS*M**5*DBARH* L/R**2)       ! Eq 26 from BPM Airfoil Self-noise and Prediction paper
            !  'B' CURVE COMPUTATION
            !        B = ABS(LOG10(STS / ST2))
            B = LOG10AA(STS / ST2) ! abs not needed absolute taken in the BMAX,BMIN   ! Eq 43 from BPM Airfoil Self-noise and Prediction paper
            BMINB = BMIN(B)
            BMAXB = BMAX(B)
            BB = BMINB + BRB0 * (BMAXB-BMINB)                              ! Eq 46 from BPM Airfoil Self-noise and Prediction paper
            SPLALPH(I)=BB+K2+10.*LOG10AA(DSTRS*M**5*DBARH*L/R**2)          ! Eq 27 from BPM Airfoil Self-noise and Prediction paper
        ELSE
            ! The 'a' computation is dropped if 'switch' is true
            SPLS(I) = 10.*LOG10AA(DSTRS*M**5*DBARL*L/R**2)
            
            !    SPLP(I) = 0.0 + 10.*LOG10(DSTRS*M**5*DBARL*L/R**2) ! changed the line below because the SPLP should be calculatd with DSTRP not with DSTRS
            SPLP(I) = 10.*LOG10AA(DSTRP*M**5*DBARL*L/R**2) ! this is correct
            
            !        B = ABS(LOG10(STS / ST2))
            B = LOG10AA(STS / ST2) ! abs not needed absolute taken in the AMAX,AMIN
            AMINB = AMIN(B)
            AMAXB = AMAX(B)
            BB = AMINB + ARA02 * (AMAXB-AMINB)
            SPLALPH(I)=BB+K2+10.*LOG10AA(DSTRS*M**5*DBARL*L/R**2)
        ENDIF
        ! Sum all contributions from 'a' and 'b' on both pressure and suction side on a mean-square pressure basis
        IF (SPLP(I)    .LT. -100.) SPLP(I)    = -100.                      ! Similar to Eq 28 of BPM Airfoil Self-noise and Prediction paper
        IF (SPLS(I)    .LT. -100.) SPLS(I)    = -100.                      ! Similar to Eq 29 of BPM Airfoil Self-noise and Prediction paper      
        IF (SPLALPH(I) .LT. -100.) SPLALPH(I) = -100.                      ! Eq 30 of BPM Airfoil Self-noise and Prediction paper recommends SPLALPH = 10log(stuff) + A' + K2, where A' is calculated same as A but with x3 Rc   

        !P1  = 10.**(SPLP(I) / 10.)            ! SPL_Pressure
        !P2  = 10.**(SPLS(I) / 10.)            ! SPL_Suction
        !P4  = 10.**(SPLALPH(I) / 10.)         ! SPL_AoA   
        !SPLTBL(I) = 10. * LOG10AA(P1 + P2 + P4)                                     ! Eq 24 from BPM Airfoil Self-noise and Prediction paper



   ENDDO

END SUBROUTINE TBLTE
!==================================================================================================================================!
SUBROUTINE TIPNOIS(ALPHTIP,ALPRAT2,C,U ,THETA,PHI, R,p,SPLTIP)
    REAL(ReKi),                               INTENT(IN   )  :: ALPHTIP        !< AOA, deg
    REAL(ReKi),                               INTENT(IN   )  :: ALPRAT2        !< TIP LIFT CURVE SLOPE                 ---
    REAL(ReKi),                               INTENT(IN   )  :: C              !< Chord Length
    REAL(ReKi),                               INTENT(IN   )  :: U              !< FREESTREAM VELOCITY               METERS/SEC
    REAL(ReKi),                               INTENT(IN   )  :: THETA          !< DIRECTIVITY ANGLE                  DEGREES 
    REAL(ReKi),                               INTENT(IN   )  :: PHI            !< DIRECTIVITY ANGLE                  DEGREES
    REAL(ReKi),                               INTENT(IN   )  :: R              !< SOURCE TO OBSERVER DISTANCE        METERS
    TYPE(AA_ParameterType) ,                  INTENT(IN   )  :: p              !< Parameters
    REAL(ReKi),DIMENSION(size(p%FreqList)),   INTENT(  OUT)  :: SPLTIP         !<
    
    ! local variables
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
    DBARH = DIRECTH_TE(M,THETA,PHI)
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
    TERM  = M*M*MM**3*L**2*DBARH/R**2                             ! TERM = M^2 * M_max^5 *l^2 *D / r^2 according to Semi-Empirical Aeroacoustic Noise Prediction Code for Wind Turbines paper
                                                                  ! Term is correct according to Eq 61 from BPM Airfoil self-noise and Prediction paper
    IF (TERM .NE. 0.0) THEN                                       
        SCALE = 10.*LOG10(TERM)
    ELSE
        SCALE = 0.0
    ENDIF
    DO I=1,size(p%FreqList)
        STPP      = p%FreqList(I) * L / UM                       ! Eq 62 from BPM Airfoil Self-noise and Prediction paper   
        SPLTIP(I) = 126.-30.5*(LOG10AA(STPP)+.3)**2 + SCALE        ! Eq 61 from BPM Airfoil Self-noise and Prediction paper
    ENDDO
END SUBROUTINE TipNois
!==================================================================================================================================!
SUBROUTINE InflowNoise(AlphaNoise,Chord,U,THETA,PHI,d,RObs,TINoise,p,SPLti)
  REAL(ReKi),                                 INTENT(IN   ) :: AlphaNoise     ! AOA, radians
  REAL(ReKi),                                 INTENT(IN   ) :: Chord          ! Chord Length
  REAL(ReKi),                                 INTENT(IN   ) :: U              !
  REAL(ReKi),                                 INTENT(IN   ) :: THETA          !
  REAL(ReKi),                                 INTENT(IN   ) :: PHI            ! Spanwise directivity angle
  REAL(ReKi),                                 INTENT(IN   ) :: d              ! element span
  REAL(ReKi),                                 INTENT(IN   ) :: RObs           ! distance to observer
!  REAL(ReKi),                                 INTENT(IN   ) :: MeanVNoise     !
  REAL(ReKi),                                 INTENT(IN   ) :: TINoise        ! turbulence intensity (NOT in percent)
!  REAL(ReKi),                                 INTENT(IN   ) :: LE_Location    !

!  REAL(ReKi),                                 INTENT(IN   ) :: dissip         !
  TYPE(AA_ParameterType),                     INTENT(IN   ) :: p              ! Parameters
  REAL(ReKi),DIMENSION(size(p%FreqList)),     INTENT(  OUT) :: SPLti          !

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
  REAL(ReKi)                   :: WaveNumber                                      ! wave number - non-dimensional frequency
  REAL(ReKi)                   :: Kbar                                      ! nafnoise 
  REAL(ReKi)                   :: khat                                      ! nafnoise 
  REAL(ReKi)                   :: ke                                        ! nafnoise 
  REAL(ReKi)                   :: tinooisess                                ! nafnoise 

  INTEGER(intKi)           :: I        !I A generic index for DO loops.

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
   DBARL = DIRECTL(Mach,THETA,PHI) ! assume that noise is low-freq in nature because turbulence length scale is large
   DBARH = DIRECTH_LE(Mach,THETA,PHI) ! Directivity for the leading edge at high frequencies
   
   IF (DBARH <= 0) THEN
      SPLti = 0.
      RETURN
   ENDIF
   
   ! In the following lines, bibliography will be referenced as:  a) Moriarty, Guidati, Migliore, Recent Improvement of a Semi-Empirical Aeroacoustic 
   ! Prediction Code for Wind Turbines (https://docs.nrel.gov/docs/fy04osti/34478.pdf)
   ! ref b) Lowson, Assessment and Prediction of Wind Turbine Noise ()
   
   !*********************************************** Model 1:
   !!! Nafnoise source code version see below 
   Frequency_cutoff = 10*U/PI/Chord
   Ke = 3.0/(4.0*p%Lturb) 
   Beta2 = 1-Mach*Mach

   DO I=1,size(p%FreqList)
      IF (p%FreqList(I) <= Frequency_cutoff) THEN
         Directivity = DBARL
      ELSE
         Directivity = DBARH 
      ENDIF

      WaveNumber = TwoPi*p%FreqList(I)/U
      Kbar = WaveNumber*Chord/2.0
      Khat = WaveNumber/Ke
      ! mu = Mach*WaveNumber*Chord/2.0/Beta2
      
      !Note: when we set RObs in CalcObserve(), we make sure it is >= AA_EPSILON ! avoid divide-by-zero
      ! tinooisess could be 0, especially on the first step, so we need to check (use LOG10AA instead of LOG10)
      SPLhigh = 10.*LOG10AA(p%AirDens**2 * p%SpdSound**4 * p%Lturb * (d/2.) / (RObs**2) *(Mach**5) * &
                            tinooisess**2 *(Khat**3)* (1+Khat**2)**(-7./3.) * Directivity) + 78.4   ! ref a; [2] )
      !bjj 01-13-2026: comparing with Eq 8 in ref [2], 
      ! (1) The paper uses "Kbar" instead of Khat (which the code uses).
      ! (2) In the paper, "I" is in percent and it adds the constant 58.4. In the code, we have "I" as a fraction and I is squared, so 
      ! 10*log10(x*100^2)+58.4 = 10*(log10(x)+log10(100^2)) + 58.4 = 10*log10(x) + 10*log10(100^2) + 58.4 = 10*log10(x) + 40 + 58.4
      ! Seems like we should be adding 98.4 instead of 78.4 in this code. However, I also haven't found documentation for the "component due to angles of attack" below,
      ! so maybe this isn't wrong.
      
   !!!   SPLhigh = 10.*LOG10(p%Lturb*(d/2.)/ &
   !!!                  (RObs*RObs)*(Mach**5)*tinooisess*tinooisess*(WaveNumber**3) &
   !!!                  *(1+WaveNumber**2)**(-7./3.)*Directivity) + 181.3  
   
      SPLhigh = SPLhigh + 10.*LOG10(1+ 9.0*AlphaNoise**2)  ! Component due to angles of attack, ref a [2])   

      Sears = 1./(TwoPi*Kbar/Beta2 + 1./(1.+2.4*Kbar/Beta2))      ! ref a [2])

   !!!   Sears = 1/(2.*PI*WaveNumber/Beta2+1/(1+2.4*WaveNumber/Beta2))  ! ref b [3]) 
   
      LFC = MAX(AA_Epsilon, 10*Sears*Mach*Kbar**2/Beta2)  ! ref a)
   !!!   LFC = 10*Sears*Mach*WaveNumber**2/Beta2  ! ref b [3])
   
   !!!   IF (mu<(PI/4.0)) THEN                     ! ref b [3])
   !!!      SPLti(I) = SPLhigh + 10.*ALOG10(LFC)   ! ref b [3])
   !!!   ELSE                                      ! ref b [3])
   !!!      SPLti(I) = SPLhigh                     ! ref b [3])
   !!!ENDIF
      SPLti(I) = SPLhigh + 10.*LOG10AA(LFC/(1+LFC))
   
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
!!      DBARL = DIRECTL(Mach,THETA,PHI) !yes, assume that noise is low-freq in nature because turbulence length scale is large
!!      DBARH = DIRECTH_LE(Mach,THETA,PHI)
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
SUBROUTINE BLUNT(ALPSTAR,C,U ,THETA,PHI,L,R,H,PSI,p,d99Var2,dstarVar1,dstarVar2,SPLBLUNT,StallVal)
  REAL(ReKi),                             INTENT(IN   )  :: ALPSTAR        ! AOA, deg
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
   ! Local variables
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
  real(ReKi)                             :: G514     ! G5 EVALUATED AT PSI=14.0              DB
  real(ReKi)                             :: HDSTARP  ! MODIFIED VALUE OF HDSTAR              ---
  real(ReKi)                             :: G50      ! G5 EVALUATED AT PSI=0.0               DB
  real(ReKi)                             :: G4       ! SCALED SPECTRUM LEVEL                 DB
  !   real(ReKi)                         :: G5       ! SPECTRUM SHAPE FUNCTION               DB
  REAL(ReKi),DIMENSION(size(p%FreqList)) :: G5       ! SPECTRUM SHAPE FUNCTION               DB ! corrected (EB_DTU)
  real(ReKi)                             :: G5Sum       ! SPECTRUM SHAPE FUNCTION               DB
  real(ReKi)                             :: SCALE    ! SCALING FACTOR                        ---
  real(ReKi)                             :: LogVal   ! temp variable to help us not take log10(0)    ---

    ! Reynolds number and mach number
        M          = U  / p%SpdSound
        RC         = U  * C/p%KinVisc
    ! Compute boundary layer thicknesses
    IF (p%X_BLMethod .eq. X_BLMethod_Tables) THEN
        DELTAP = d99Var2
        DSTRS  = dstarVar1
        DSTRP  = dstarVar2
    ELSE
        CALL THICK(C,RC,ALPSTAR,p,DELTAP,DSTRS,DSTRP,StallVal)
    ENDIF
    
    ! Compute average displacement thickness
    DSTRAVG = (DSTRS + DSTRP) / 2.
    HDSTAR  = H / DSTRAVG
    DSTARH = 1. /HDSTAR
    ! Compute directivity function
    DBARH = DIRECTH_TE(M,THETA,PHI)
    IF (DBARH <= 0) THEN
        SPLBLUNT = 0.
        RETURN
    ENDIF
    
    ! Compute peak strouhal number                                               eq 72 in BPM Airfoil Self-noise and Prediction paper
    ATERM  = .212 - .0045 * PSI
    IF (HDSTAR .GE. .2) then
       STPEAK    = ATERM / (1.+.235*DSTARH-.0132*DSTARH**2)                    ! this is what it used to be in nafnoise and fast noise module
    !!  STPEAK    = ATERM / (1+0.235*(DSTARH)**(-1)-0.0132*DSTARH**(-2)) ! check if this one is correct (EB_DTU) 
    else
       STPEAK    = .1 * HDSTAR + .095 - .00243 * PSI
    end if
    
    ! Compute scaled spectrum level                                              eq 74 of BPM Airfoil Self-noise and Prediction paper
    if (HDSTAR .LE. 5.) then
       G4=17.5*LOG10AA(HDSTAR)+157.5-1.114*PSI
    else
       G4=169.7 - 1.114 * PSI
    end if
    
    ! For each frequency, compute spectrum shape referenced to 0 db
    SCALE = 10. * LOG10AA(M**5.5 * H * DBARH * L / R**2)
    G5Sum=0.0_Reki
    DO I=1,SIZE(p%FreqList)
        STPPP    = p%FreqList(I) * H / U
        ETA      = LOG10AA(STPPP/STPEAK)
        G514 = G5COMP(HDSTAR,ETA)                          ! compute G5 for Phi=14deg
        
        HDSTARP = 6.724 * HDSTAR **2-4.019*HDSTAR+1.107                         ! eq 82 from BPM Airfoil Self-noise and Prediction paper
        G50 = G5COMP(HDSTARP,ETA)                           ! recompute G5 for Phi=0deg
        
        G5(I) = G50 + .0714 * PSI * (G514-G50)                                   ! interpolate G5 from G50 and G514
        IF (G5(I) .GT. 0.) G5(I) = 0.
        G5Sum = 10**(G5(I)/10)+G5Sum     ! to be subtracted
        if ( G5Sum .ne. 0) then
            LogVal = MAX(AA_EPSILON,1/G5Sum)
        else
           LogVal = 1
        end if
        SPLBLUNT(I) = G4 + G5(I) + SCALE - 10*log10(LogVal)  ! equation mentioned there is plus but it is stated subtract, thus ''- 10*log10(1/G5Sum)'' 
    end do
END SUBROUTINE Blunt
!====================================================================================================
REAL(ReKi) FUNCTION G5COMP(HDSTAR,ETA) result(G5)
    REAL(ReKi),          INTENT(IN   )  :: HDSTAR        !<
    REAL(ReKi),          INTENT(IN   )  :: ETA           !< 

    ! Local variables
    real(ReKi)                                    :: K 
    real(ReKi)                                    :: M
    real(ReKi)                                    :: MU
    real(ReKi)                                    :: ETA0
    
    IF ( HDSTAR .LT. .25) then
       MU = .1211                    ! begin eq 78 from BPM Airfoil Self-noise and Prediction paper
    elseif (HDSTAR .LE. .62) then
       MU =-.2175*HDSTAR + .1755
    elseif (HDSTAR .LT. 1.15) then
       MU =-.0308*HDSTAR + .0596
    else
       MU = .0242
    end if
    
    IF ( HDSTAR .LE. .02 ) then
       M = 0.0                       ! begin eq 79 from BPM Airfoil Self-noise and Prediction paper
    elseif (HDSTAR .LT. 0.5) then
       M = 68.724*HDSTAR - 1.35
    elseif (HDSTAR .LE. .62) then
       M = 308.475*HDSTAR - 121.23
    elseif (HDSTAR .LE. 1.15) then
       M = 224.811*HDSTAR - 69.354
    elseif (HDSTAR .LT. 1.2) then
       M = 1583.28*HDSTAR - 1631.592
    else
       M = 268.344
    end if
    M = MAX(M, 0.0_ReKi) !bjj: not sure this is necessary... previous iterations of this statement missed some of the cases so may have had uninitialized values; otherwise, it's not possible to get M<0
    
    ETA0 = -SQRT((M*M*MU**4)/(6.25+M*M*MU*MU))                                   ! eq 80 from BPM Airfoil Self-noise and Prediction paper
    
    IF (ETA .LE. ETA0) then
       K  = 2.5*SQRT(1.-(ETA0/MU)**2)-2.5-M*ETA0                                 ! eq 81 from BPM Airfoil Self-noise and Prediction paper
       G5 = M * ETA + K                     ! begin eq 76 from BPM Airfoil Self-noise and Prediction paper
    elseif (ETA .LE. 0.) then
       G5 = 2.5*SQRT(1.-(ETA/MU)**2)-2.5
    elseif (ETA .LE. 0.03615995) then
       G5 = SQRT(1.5625-1194.99*ETA**2)-1.25
    else
       G5 = -155.543 * ETA + 4.375
    end if
    
END FUNCTION G5Comp
!====================================================================================================
!> This subroutine defines the curve fit corresponding to the a-curve for the minimum allowed reynolds number.
REAL(ReKi) FUNCTION AMIN(A) result(AMINA)
    REAL(ReKi),                             INTENT(IN   )  :: A
    REAL(ReKi) :: X1
    
    X1 = ABS(A)
    IF (X1 .LE. .204) then
       AMINA=SQRT(67.552-886.788*X1**2)-8.219
    elseif (X1 .LE. .244) then
       AMINA=-32.665*X1+3.981
    else
       AMINA=-142.795*X1**3+103.656*X1**2-57.757*X1+6.006
    end if
    
END FUNCTION AMIN
!====================================================================================================
!> This subroutine defines the curve fit corresponding to the a-curve for the maximum allowed reynolds number.
REAL(ReKi) FUNCTION AMAX(A) result(AMAXA)
    REAL(ReKi),                             INTENT(IN   )  :: A
    REAL(ReKi) :: X1

    X1 = ABS(A)
    IF (X1 .LE. .13) then
       AMAXA=SQRT(67.552-886.788*X1**2)-8.219
    elseif (X1 .LE. .321) then
       AMAXA=-15.901*X1+1.098
    else
       AMAXA=-4.669*X1**3+3.491*X1**2-16.699*X1+1.149
    end if
    
END FUNCTION AMAX
!====================================================================================================
!> This subroutine defines the curve fit corresponding to the b-curve for the minimum allowed reynolds number.
REAL(ReKi) FUNCTION BMIN(B) result(BMINB)
    REAL(ReKi),                             INTENT(IN   )  :: B
    REAL(ReKi) :: X1

    X1 = ABS(B)
    IF (X1 .LE. .13) then
       BMINB=SQRT(16.888-886.788*X1**2)-4.109
    elseif (X1 .LE. .145) then
       BMINB=-83.607*X1+8.138
    else
       BMINB=-817.81*X1**3+355.21*X1**2-135.024*X1+10.619
    end if
    
END FUNCTION BMin
!====================================================================================================
!> Define the curve fit corresponding to the b-curve for the maximum allowed reynolds number.
REAL(ReKi) FUNCTION BMAX(B) result(BMAXB)
    REAL(ReKi),   INTENT(IN   )  :: B
    REAL(ReKi) :: X1
    X1 = ABS(B)
    IF (X1 .LE. .1) then
       BMAXB=SQRT(16.888-886.788*X1**2)-4.109
    else if (X1 .LE. .187) then
       BMAXB=-31.313*X1+1.854
    else
       BMAXB=-80.541*X1**3+44.174*X1**2-39.381*X1+2.344
    end if
END FUNCTION BMax
!====================================================================================================
!> Determine where the a-curve takes on a value of -20 db.
REAL(ReKi) FUNCTION A0COMP(RC) result(A0)
    REAL(ReKi),   INTENT(IN   )  :: RC
    IF (RC .LT. 9.52E+04) then
       A0 = .57
    elseif (RC .LT. 8.57E+05) then
       A0 = (-9.57E-13)*(RC-8.57E+05)**2 + 1.13
    else
       A0 = 1.13
    end if
END FUNCTION A0COMP
!====================================================================================================
!> Compute zero angle of attack boundary layer thickness (meters) and reynolds number
SUBROUTINE THICK(C,RC,ALPSTAR,p,DELTAP,DSTRS,DSTRP,StallVal)
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
    REAL(ReKi),                INTENT(IN   )  :: ALPSTAR        !< AOA, deg
    REAL(ReKi),                INTENT(IN   )  :: C              !< Chord Length
    REAL(ReKi),                INTENT(IN   )  :: RC             !< RC= U*C/KinViscosity
    TYPE(AA_ParameterType),    INTENT(IN   )  :: p              !< Parameters
    REAL(ReKi),                INTENT(  OUT)  :: DELTAP         !< Pressure side boundary layer thickness
    REAL(ReKi),                INTENT(  OUT)  :: DSTRS          !< Suction side displacement thickness
    REAL(ReKi),                INTENT(  OUT)  :: DSTRP          !< Pressure side displacement thickness
    REAL(ReKi),                INTENT(IN   )  :: StallVal       !< Stall angle at station i

    ! Local variables
!    integer(intKi)          :: ErrStat2           ! temporary Error status
!    character(ErrMsgLen)    :: ErrMsg2            ! temporary Error message
    character(*), parameter :: RoutineName = 'Thick'
    real(ReKi)              :: DELTA0              ! BOUNDARY LAYER THICKNESS AT ZERO ANGLE OF ATTACK METERS
    real(ReKi)              :: DSTR0      ! DISPLACEMENT THICKNESS AT ZERO   ANGLE OF ATTACK METERS
    real(ReKi)              :: LogRC      ! LOG10(RC)
    
    LogRC = LOG10AA( RC )
    
    ! Boundary layer thickness
    DELTA0                     = 10.**(1.6569-0.9045*LogRC+0.0596*LogRC**2)*C ! (untripped)         Eq. (5) of [1]
    IF (p%ITRIP /= ITRIP_None) DELTA0 = 10.**(1.892 -0.9045*LogRC+0.0596*LogRC**2)*C ! (heavily tripped)   Eq. (2) of [1]
    IF (p%ITRIP .EQ. ITRIP_Light) DELTA0=.6*DELTA0
    
    ! Pressure side boundary layer thickness, Eq (8) of [1]
    DELTAP   = 10.**(-.04175*ALPSTAR+.00106*ALPSTAR**2)*DELTA0
    
    ! Compute zero angle of attack displacement thickness
    IF (p%ITRIP /= ITRIP_None) THEN
        ! Heavily tripped, Eq. (3) of [1]
        IF (RC .LE. .3E+06) THEN
           DSTR0 = .0601 * RC **(-.114)*C
        ELSE
           DSTR0=10.**(3.411-1.5397*LogRC+.1059*LogRC**2)*C
        END IF
        ! Lightly tripped
        IF (p%ITRIP .EQ. ITRIP_Light) DSTR0 = DSTR0 * .6
    ELSE
        ! Untripped, Eq. (6) of [1]
        DSTR0=10.**(3.0187-1.5397*LogRC+.1059*LogRC**2)*C
    ENDIF
    
    ! Pressure side displacement thickness, Eq. (9) of [1]
   DSTRP   = 10.**(-.0432*ALPSTAR+.00113*ALPSTAR**2)*DSTR0
    !      IF (p%ITRIP .EQ. 3) DSTRP = DSTRP * 1.48 ! commented since itrip is never 3 check if meant 2.(EB_DTU)

    ! Suction side displacement thickness
   IF (p%ITRIP .EQ. ITRIP_Heavy) THEN
      ! Heavily tripped, Eq. (12) of [1]
      IF (ALPSTAR .LE. 5.) THEN
         DSTRS=10.**(.0679*ALPSTAR)*DSTR0
      ELSEIF (ALPSTAR .LE. StallVal) THEN
         DSTRS = 0.381 * 10.**(.1516*ALPSTAR)*DSTR0
      ELSE
         DSTRS = 14.296 * 10.**(.0258*ALPSTAR)*DSTR0
      ENDIF
   ELSE
        ! Untripped or lightly tripped, Eq. (15) of [1]
      IF (ALPSTAR .LE. 7.5) THEN
         DSTRS =10.**(.0679*ALPSTAR)*DSTR0
      ELSEIF(ALPSTAR .LE. StallVal) THEN
         DSTRS = .0162*10.**(.3066*ALPSTAR)*DSTR0
      ELSE
         DSTRS = 52.42*10.**(.0258*ALPSTAR)*DSTR0
      ENDIF
   ENDIF
   
END SUBROUTINE Thick
!====================================================================================================
!> This subroutine computes the high frequency directivity function for the trailing edge
REAL(ReKi) FUNCTION DIRECTH_TE(M,THETA,PHI) result(DBAR)
    REAL(ReKi),        INTENT(IN   ) :: THETA      !
    REAL(ReKi),        INTENT(IN   ) :: PHI        !
    REAL(ReKi),        INTENT(IN   ) :: M          !
    ! Local variables
    real(ReKi)              :: MC
    real(ReKi), parameter   :: DEGRAD = .017453
    real(ReKi)              :: PHIR
    real(ReKi)              :: THETAR

    MC     = .8 * M
    THETAR = THETA * DEGRAD
    PHIR   = PHI * DEGRAD
    DBAR   = 2.*SIN(THETAR/2.)**2 * SIN(PHIR)**2 / ((1.+M*COS(THETAR))* (1.+(M-MC)*COS(THETAR))**2)    ! eq B1 in BPM Airfoil Self-noise and Prediction paper
END FUNCTION DIRECTH_TE

!====================================================================================================
!> This subroutine computes the high frequency directivity function for the leading edge
REAL(ReKi) FUNCTION DIRECTH_LE(M,THETA,PHI) result(DBAR)
    REAL(ReKi),        INTENT(IN   ) :: THETA      !
    REAL(ReKi),        INTENT(IN   ) :: PHI        !
    REAL(ReKi),        INTENT(IN   ) :: M          !

    ! Local variables
    real(ReKi), parameter   :: DEGRAD = .017453
    real(ReKi)              :: PHIR
    real(ReKi)              :: THETAR

    THETAR = THETA * DEGRAD
    PHIR   = PHI * DEGRAD
    DBAR   = 2.*COS(THETAR/2.)**2*SIN(PHIR)**2/(1.+M*COS(THETAR))**3 
END FUNCTION DIRECTH_LE

!====================================================================================================
!> This subroutine computes the high frequency directivity function for the input observer location
! Paper: 
REAL(ReKi) FUNCTION DIRECTL(M,THETA,PHI) result(DBAR)
    REAL(ReKi),           INTENT(IN   ) :: THETA      !<
    REAL(ReKi),           INTENT(IN   ) :: PHI        !<
    REAL(ReKi),           INTENT(IN   ) :: M          !<

    ! Local variables
    real(ReKi)              :: MC
    real(ReKi), parameter   :: DEGRAD  = .017453
    real(ReKi)              :: PHIR
    real(ReKi)              :: THETAR

    !   This subroutine computes the low frequency directivity function for the input observer location
    
    MC     = .8 * M
    THETAR = THETA * DEGRAD
    PHIR   = PHI * DEGRAD
    DBAR = (SIN(THETAR)*SIN(PHIR))**2/(1.+M*COS(THETAR))**4                   ! eq B2 in BPM Airfoil Self-noise and Prediction paper
END FUNCTION DIRECTL
!==================================================================================================================================!
!===============================  Simplified Guidati Inflow Turbulence Noise Addition =============================================!
!==================================================================================================================================!
! Uses simple correction for turbulent inflow noise from Moriarty et. al 2005
! Paper: Prediction of Turbulent Inflow and Trailing-Edge Noise for Wind Turbines, by Moriarty, Guidati, and Migliore
SUBROUTINE Simple_Guidati(U,Chord,thick_10p,thick_1p,p,SPLti)
    REAL(ReKi),                             INTENT(IN   )  :: U              ! Vrel
    REAL(ReKi),                             INTENT(IN   )  :: Chord          ! Chord Length
    REAL(ReKi),                             INTENT(IN   )  :: thick_10p      ! 
    REAL(ReKi),                             INTENT(IN   )  :: thick_1p       ! 
    TYPE(AA_ParameterType),                 INTENT(IN   )  :: p              ! Parameters
    REAL(ReKi),DIMENSION(size(p%FreqList)), INTENT(  OUT)  :: SPLti          !
    ! local variables
    character(*), parameter :: RoutineName = 'Simple_Guidati'
    INTEGER(intKi)          :: loop1       ! temporary
    REAL(ReKi)              :: TI_Param    ! Temporary variable thickness ratio dependent
    REAL(ReKi)              :: slope       ! Temporary variable thickness ratio dependent
    REAL(ReKi)              :: const1      ! Temporary variable
    REAL(ReKi)              :: const2      ! Temporary variable
    
    TI_Param = thick_1p + thick_10p                                     ! Eq 2 
    slope = 1.123*TI_Param + 5.317*TI_Param*TI_Param                    ! Eq 3 
    const1 = -slope*TwoPi*chord/U
    const2 = -slope*5.0d0
    
    do loop1 =1,size(p%FreqList)
!       SPLti(loop1) = -slope*(TwoPi * chord/U * p%FreqList(loop1) + 5.0d0)  ! Eq 4 
        SPLti(loop1) = const1 * p%FreqList(loop1) + const2  ! Eq 4 
    enddo   ! Outputs Delta_SPL, the difference in SPL between the airfoil and a flat plate.
    
END SUBROUTINE Simple_Guidati
!==================================================================================================================================!
!================================ Turbulent Boundary Layer Trailing Edge Noise ====================================================!
!=================================================== TNO START ====================================================================!
SUBROUTINE TBLTE_TNO(U,THETA,PHI,D,R,Cfall,d99all,EdgeVelAll,p,SPLP,SPLS)
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
!   REAL(ReKi),DIMENSION(size(p%FreqList)),   INTENT(IN   ) :: SPLALPH    !< SOUND PRESSURE LEVEL DUE TO ANGLE OF ATTACK CONTRIBUTION (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),   INTENT(  OUT) :: SPLP       !< SOUND PRESSURE LEVEL DUE TO PRESSURE SIDE OF AIRFOIL (db)
    REAL(ReKi),DIMENSION(size(p%FreqList)),   INTENT(  OUT) :: SPLS       !< SOUND PRESSURE LEVEL DUE TO SUCTION SIDE OF AIRFOIL  (db)

    ! Local variables
    REAL(ReKi) :: answer
    REAL(ReKi) :: Spectrum
    REAL(ReKi) :: freq(size(p%FreqList))
    REAL(ReKi) :: SPL_press,SPL_suction
    REAL(ReKi) :: band_width,band_ratio
    REAL(ReKi) :: DBARH
    !REAL(ReKi) :: P1,P2,P4
    INTEGER (4)  :: n_freq
    INTEGER (4)  :: i_omega

      ! Variables passed to integration routine
   real(ReKi)  :: int_limits(2)  !< Lower and upper integration limits
   real(ReKi)  :: Mach        !< Mach number
   real(ReKi)  :: omega

    ! Init
    n_freq  = size(p%FreqList)
    freq    = p%FreqList
    
    SPLS = 0.0_ReKi ! initialize in case Cfall(1) <= 0
    SPLP = 0.0_ReKi ! initialize in case Cfall(2) <= 0
    
    ! Body of TNO 
    band_ratio = 2.**(1./3.)

    ! Mach number
    Mach = U  / p%SpdSound

    ! Directivity function
    DBARH = DIRECTH_TE(REAL(Mach,ReKi),THETA,PHI)
 
    do i_omega = 1,n_freq
        omega = TwoPi*p%FreqList(i_omega)
        !integration limits
        int_limits(1) = 0.0e0
        int_limits(2) = 10*omega/(Mach*p%SpdSound)
        ! Convert to third octave
        band_width = 2. * omega * (sqrt(band_ratio)-1./sqrt(band_ratio))
        
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

        !P1  = 10.**(SPLP(i_omega) / 10.)
        !P2  = 10.**(SPLS(i_omega) / 10.)
        !P4  = 10.**(SPLALPH(i_omega) / 10.)
        !
        !SPLTBL(i_omega) = 10. * LOG10(P1 + P2 + P4)
    enddo
END SUBROUTINE TBLTE_TNO


!====================================================================================================
SUBROUTINE BL_Param_Interp(p,m,U,AlphaNoise_Deg,C,whichAirfoil)
   TYPE(AA_ParameterType),                INTENT(IN   ) :: p              !< Parameters
   TYPE(AA_MiscVarType),                  INTENT(INOUT) :: m              !< misc/optimization data (not defined in submodules)
   REAL(ReKi),                            INTENT(IN   ) :: U              !< METERS/SEC
   REAL(ReKi),                            INTENT(IN   ) :: AlphaNoise_Deg     !< Angle of Attack                           DEG
   REAL(ReKi),                            INTENT(IN   ) :: C              !< Chord                                     METERS
   integer(intKi),                        INTENT(IN   ) :: whichAirfoil   !< whichairfoil
   
   character(*), parameter :: RoutineName = 'BL_Param_Interp'
   REAL(ReKi)              :: RC
   INTEGER(intKi)          :: i
  
   INTEGER, PARAMETER :: NumDimensions = 2
   INTEGER(IntKi)                                :: MaxIndx(NumDimensions)                       ! max sizes associated with each dimension of array
   INTEGER(IntKi)                                :: Indx_Lo(NumDimensions)                       ! index associated with lower bound of dimension 1,2 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   INTEGER(IntKi)                                :: Indx_Hi(NumDimensions)                       ! index associated with upper bound of dimension 1,2 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   REAL(ReKi)                                    :: Pos_Lo(NumDimensions)                        ! coordinate value with lower bound of dimension 1,2
   REAL(ReKi)                                    :: Pos_Hi(NumDimensions)                        ! coordinate value with upper bound of dimension 1,2

   REAL(ReKi)                                    :: isopc(NumDimensions)                         ! isoparametric coordinates
   REAL(ReKi)                                    :: N(2**NumDimensions)                          ! size 2^n
   REAL(ReKi)                                    :: InCoord(NumDimensions)                       !< Arranged as (x, y)


  !!!! this if is not used but if necessary two sets of tables can be populated for tripped and untripped cases
   RC = U  * C/p%KinVisc       ! REYNOLDS NUMBER BASED ON  CHORD


      ! find the indices into the arrays representing coordinates of each dimension:
      !  (by using LocateStp, we do not require equally spaced arrays)
   InCoord = (/ AlphaNoise_Deg, RC /)
   
   MaxIndx(1) = SIZE(p%AOAListBL)
   MaxIndx(2) = SIZE(p%ReListBL)

   CALL LocateStp( InCoord(1), p%AOAListBL, m%LastIndex(1), MaxIndx(1) )
   CALL LocateStp( InCoord(2), p%ReListBL,  m%LastIndex(2), MaxIndx(2) )

   Indx_Lo = m%LastIndex  ! at this point, 0 <= Indx_Lo(i) <= n(i) for all i

   ! RE (indx 2)
   do i = 1,2
      IF (Indx_Lo(i) == 0) THEN
         Indx_Lo(i) = 1
      ELSEIF (Indx_Lo(i) == MaxIndx(i) ) THEN
         Indx_Lo(i) = max( MaxIndx(i) - 1, 1 )                ! make sure it's a valid index
      END IF
      Indx_Hi(i) = min( Indx_Lo(i) + 1 , MaxIndx(i) )         ! make sure it's a valid index
   end do
  
   ! calculate the bounding box; the positions of all dimensions:

   pos_Lo(1) = p%AOAListBL( Indx_Lo(1) )
   pos_Hi(1) = p%AOAListBL( Indx_Hi(1) )

   pos_Lo(2) = p%ReListBL( Indx_Lo(2) )
   pos_Hi(2) = p%ReListBL( Indx_Hi(2) )


      ! 2-D linear interpolation:

   CALL IsoparametricCoords( InCoord, pos_Lo, pos_Hi, isopc )      ! Calculate iospc

   N(1)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )
   N(2)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )
   N(3)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )
   N(4)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )
   N     = N / REAL( SIZE(N), ReKi )  ! normalize
          
   m%dStarVar  (1) = InterpData( p%dstarall1(  :,:,whichAirfoil) )
   m%dStarVar  (2) = InterpData( p%dstarall2(  :,:,whichAirfoil) )
   m%d99Var    (1) = InterpData( p%d99all1(    :,:,whichAirfoil) )
   m%d99Var    (2) = InterpData( p%d99all2(    :,:,whichAirfoil) )
   m%CfVar     (1) = InterpData( p%Cfall1(     :,:,whichAirfoil) )
   m%CfVar     (2) = InterpData( p%Cfall2(     :,:,whichAirfoil) )
   m%EdgeVelVar(1) = InterpData( p%EdgeVelRat1(:,:,whichAirfoil) )
   m%EdgeVelVar(2) = InterpData( p%EdgeVelRat2(:,:,whichAirfoil) )
  
contains
   real(ReKi) function InterpData(Dataset)
      REAL(ReKi),                  INTENT(IN   ) :: Dataset(:,:)                                 !< Arranged as (x, y)
      REAL(ReKi)                                 :: u(2**NumDimensions)                          ! size 2^n

      u(1)  = Dataset( Indx_Hi(1), Indx_Lo(2) )
      u(2)  = Dataset( Indx_Hi(1), Indx_Hi(2) )
      u(3)  = Dataset( Indx_Lo(1), Indx_Hi(2) )
      u(4)  = Dataset( Indx_Lo(1), Indx_Lo(2) )

      InterpData = SUM ( N * u )
   
   end function
END SUBROUTINE BL_Param_Interp



SUBROUTINE Aero_Tests()
    !--------Laminar Boundary Layer Vortex Shedding Noise----------------------------!
    !CALL LBLVS(AlphaNoise_Deg,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
    !    elementspan,m%rTEtoObserve(K,J,I), &
    !    p,m%d99Var(2),m%dstarVar(1),m%dstarVar(2),m%SPLLBL)
    !--------Turbulent Boundary Layer Trailing Edge Noise----------------------------!
    !CALL TBLTE(3.0d0,0.22860d0,63.920d0,90.0d0,90.0d0,0.5090d0,1.220d0, &
    !    p, m%d99Var(2),m%dstarVar(1),m%dstarVar(2),p%StallStart(J,I),m%SPLP,m%SPLS,m%SPLALPH,m%SPLTBL )
    !m%SPLP=0.0_ReKi;m%SPLS=0.0_ReKi;m%SPLTBL=0.0_ReKi;
    !m%EdgeVelVar(1)=1.000d0;m%EdgeVelVar(2)=m%EdgeVelVar(1);
    !m%CfVar(1) = 0.0003785760d0;m%CfVar(2) = 0.001984380d0;m%d99var(1)= 0.01105860d0; m%d99var(2)= 0.007465830d0;m%EdgeVelVar(1)=1.000d0;m%EdgeVelVar(2)=m%EdgeVelVar(1);
    !CALL TBLTE_TNO(0.22860_Reki,63.9200_Reki,90.00_Reki,90.0_Reki,0.5090_Reki,1.220_Reki, &
    !    m%CfVar,m%d99var,m%EdgeVelVar, p, m%SPLP,m%SPLS,m%SPLALPH,m%SPLTBL,ErrStat2 ,errMsg2)
    !--------Blunt Trailing Edge Noise----------------------------------------------!
    !CALL BLUNT(3.0d0,0.22860d0,63.920d0,90.0d0,90.0d0,0.5090d0,1.220d0,&
    !    p%TEThick(J,I),p%TEAngle(J,I),p, m%d99Var(2),m%dstarVar(1),m%dstarVar(2),m%SPLBLUNT )
    !--------Tip Noise--------------------------------------------------------------!
    !CALL TIPNOIS(AlphaNoise_Deg,p%ALpRAT,p%BlChord(J,I),UNoise,m%ChordAngleTE(K,J,I),m%SpanAngleTE(K,J,I), &
    !    m%rTEtoObserve(K,J,I), p, m%SPLTIP,ErrStat2,errMsg2)
    !--------Inflow Turbulence Noise ------------------------------------------------!
    !CALL InflowNoise(3.0d0,0.22860d0,63.920d0,90.0d0,90.0d0,0.5090d0,1.220d0, xd%TIVx(J,I),0.050d0,p,m%SPLti )
    !CALL FullGuidati(3.0d0,63.920d0,0.22860d0,0.5090d0,1.220d0,90.0d0,90.0d0,xd%MeanVrel(J,I),xd%TIVrel(J,I), &
    !    p,p%BlAFID(J,I),m%SPLTIGui,ErrStat2 )
    !CALL Simple_Guidati(UNoise,0.22860d0,0.120d0,0.020d0,p,m%SPLTIGui)
END SUBROUTINE 
END MODULE AeroAcoustics

