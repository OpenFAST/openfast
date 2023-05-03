!**********************************************************************************************************************************
! The WAMIT and WAMIT_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. WAMITs_Types will be auto-generated based on a description of the variables for the module.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
!
!    This file is part of WAMIT.
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
!    
!**********************************************************************************************************************************
MODULE WAMIT

   USE Waves
   USE WAMIT_Types 
   USE WAMIT_Interp
   USE NWTC_Library
  ! USE Waves_Types
   USE Conv_Radiation
   USE SS_Radiation
   USE SS_Excitation
   USE NWTC_FFTPACK
   use SeaState_Interp
   
   IMPLICIT NONE
   
   PRIVATE
   
   REAL(DbKi), PARAMETER, PRIVATE       :: OnePlusEps  = 1.0 + EPSILON(OnePlusEps)   ! The number slighty greater than unity in the precision of DbKi.

   TYPE(ProgDesc), PARAMETER            :: WAMIT_ProgDesc = ProgDesc( 'WAMIT', '', '' )

   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: WAMIT_Init                           ! Initialization routine
   PUBLIC :: WAMIT_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: WAMIT_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: WAMIT_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: WAMIT_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: WAMIT_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: WAMIT_UpdateDiscState                ! Tight coupling routine for updating discrete states
        
   
CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine transforms  WAMIT input file data from a local (heading-angle, based) coordinate system to the global system. 
   subroutine TransformWAMITMatrices( NBody, RotZ, M )
!..................................................................................................................................
      integer(IntKi), intent( in    ) :: NBody   ! Number of WAMIT bodies in this WAMIT object ( = 1 if NBodyMod > 1)
      real(R8Ki),     intent( in    ) :: RotZ(:) ! NBody heading angles (radians)
      real(SiKi),     intent( inout ) :: M(:,:)  ! Matrix data to be transformed, if NBodyMOD = 1 and NBody > 1 then we will be transforming the individual sub 6x6 matrices
      
      integer(IntKi)   :: i,j,ii,jj,iSub,jSub
      real(R8Ki)       :: Rj(3,3)
      real(R8Ki)       :: Ri(3,3)
      
      do j = 1, NBody
         Rj(1,:) = (/ cos(RotZ(j)), sin(RotZ(j)), 0.0_R8Ki/)
         Rj(2,:) = (/-sin(RotZ(j)), cos(RotZ(j)), 0.0_R8Ki/)
         Rj(3,:) = (/ 0.0_R8Ki    , 0.0_R8Ki    , 1.0_R8Ki/)
         do i = 1, NBody
            if ( (.not. EqualRealNos(RotZ(i), 0.0_R8Ki)) .or. (.not. EqualRealNos(RotZ(j), 0.0_R8Ki)) ) then
               Ri(1,:) = (/ cos(RotZ(i)), sin(RotZ(i)), 0.0_R8Ki/)
               Ri(2,:) = (/-sin(RotZ(i)), cos(RotZ(i)), 0.0_R8Ki/)
               Ri(3,:) = (/ 0.0_R8Ki    , 0.0_R8Ki    , 1.0_R8Ki/)
               do jj = 1,2
                  jSub = (j-1)*6 + (jj-1)*3 + 1  
                  do ii = 1,2
                     iSub = (i-1)*6 + (ii-1)*3 + 1
                     M(iSub:iSub+2,jSub:jSub+2) = matmul( transpose(Ri), matmul( M(iSub:iSub+2,jSub:jSub+2), Rj ) )
                  end do
               end do 
            end if
         end do

      end do
   end subroutine TransformWAMITMatrices
   
   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE WAMIT_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(WAMIT_InitInputType),       INTENT(INOUT)  :: InitInp       !< Input data for initialization routine.  NOTE: we need INOUT because we may be moving the allocation of SS_Excitation data
      TYPE(WAMIT_InputType),           INTENT(  OUT)  :: u             !< An initial guess for the input; input mesh must be defined
      TYPE(WAMIT_ParameterType),       INTENT(  OUT)  :: p             !< Parameters      
      TYPE(WAMIT_ContinuousStateType), INTENT(  OUT)  :: x             !< Initial continuous states
      TYPE(WAMIT_DiscreteStateType),   INTENT(  OUT)  :: xd            !< Initial discrete states
      TYPE(WAMIT_ConstraintStateType), INTENT(  OUT)  :: z             !< Initial guess of the constraint states
      TYPE(WAMIT_OtherStateType),      INTENT(  OUT)  :: OtherState    !< Initial other states            
      TYPE(WAMIT_OutputType),          INTENT(  OUT)  :: y             !< Initial system outputs (outputs are not calculated; 
                                                                       !!   only the output mesh is initialized)
      TYPE(WAMIT_MiscVarType),         INTENT(  OUT)  :: m             !< Initial misc/optimization variables            
      REAL(DbKi),                      INTENT(IN   )  :: Interval      !< Coupling interval in seconds: the rate that 
                                                                       !!   (1) WAMIT_UpdateStates() is called in loose coupling &
                                                                       !!   (2) WAMIT_UpdateDiscState() is called in tight coupling.
                                                                       !!   Input is the suggested time from the glue code; 
                                                                       !!   Output is the actual coupling interval that will be used 
                                                                       !!   by the glue code.
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat       !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg        !< Error message if ErrStat /= ErrID_None


     
         
      
         ! These are dummy variables to satisfy the framework, but are not used 
      TYPE(Conv_Rdtn_InitInputType)            :: Conv_Rdtn_InitInp                     ! Local version of the intialization data for the radiation module
      TYPE(Conv_Rdtn_InitOutputType)           :: Conv_Rdtn_InitOut                     ! Initialization Outputs from the Conv_Rdtn module initialization
      !TYPE(Conv_Rdtn_InitOutputType)          :: Conv_RdtnInitOutData                     
      TYPE(SS_Rad_InitInputType)               :: SS_Rdtn_InitInp                       ! Local version of the intialization data for the radiation module
      TYPE(SS_Rad_InitOutputType)              :: SS_Rdtn_InitOut                       ! Initialization Outputs from the SS_Rdtn module initialization
      TYPE(SS_Exc_InitInputType)               :: SS_Exctn_InitInp                       ! Local version of the intialization data for the SS wave excitation module
      TYPE(SS_Exc_InitOutputType)              :: SS_Exctn_InitOut                       ! Initialization Outputs from the SS wave excitation module initialization
     
       
         ! Local Variables
      REAL(DbKi)                             :: Interval_Sub                         ! Local timestep for the SS_Rad and SS_Exc modules, based on RdtnDT   
      COMPLEX(SiKi), ALLOCATABLE             :: HdroExctn (:,:,:)                    ! Frequency- and direction-dependent complex hydrodynamic wave excitation force per unit wave amplitude vector (kg/s^2, kg-m/s^2)
      COMPLEX(SiKi), ALLOCATABLE             :: WaveExctnC(:,:)                      ! Discrete Fourier transform of the instantaneous value of the total excitation force on the support platfrom from incident waves (N, N-m)
      COMPLEX(SiKi), ALLOCATABLE             :: WaveExctnCGrid(:,:,:)                ! Discrete Fourier transform of the instantaneous value of the total excitation force on the grid points from incident waves (N, N-m)
      REAL(ReKi)                             :: DffrctDim (6)                        ! Matrix used to redimensionalize WAMIT hydrodynamic wave excitation force  output (kg/s^2, kg-m/s^2            )
      REAL(SiKi), ALLOCATABLE                :: HdroAddMs (:,:,:)                    ! The frequency-dependent hydrodynamic added mass matrix from the radiation problem (kg  , kg-m  , kg-m^2  )
      REAL(SiKi), ALLOCATABLE                :: HdroDmpng (:,:,:)                    ! The frequency-dependent hydrodynamic damping    matrix from the radiation problem (kg/s, kg-m/s, kg-m^2/s)
      REAL(SiKi), ALLOCATABLE                :: HdroFreq  (:)                        ! Frequency components inherent in the hydrodynamic added mass matrix, hydrodynamic daming matrix, and complex wave excitation force per unit wave amplitude vector (rad/s)
      REAL(SiKi), ALLOCATABLE                :: HdroWvDir (:)                        ! Incident wave propagation heading direction components inherent in the complex wave excitation force per unit wave amplitude vector (degrees)
      REAL(ReKi)                             :: HighFreq                             ! The highest frequency component in the WAMIT file, not counting infinity.
      REAL(SiKi)                             :: Omega                                ! Wave frequency (rad/s)
      REAL(ReKi)                             :: PrvDir                               ! The value of TmpDir from the previous line (degrees)
      REAL(ReKi)                             :: PrvPer                               ! The value of TmpPer from the previous line (sec    )
      REAL(ReKi)                             :: SttcDim   (6,6)                      ! Matrix used to redimensionalize WAMIT hydrostatic  restoring              output (kg/s^2, kg-m/s^2, kg-m^2/s^2)
      REAL(ReKi)                             :: RdtnDim   (6,6)                      ! Matrix used to redimensionalize WAMIT hydrodynamic added mass and damping output (kg    , kg-m    , kg-m^2    )
      REAL(ReKi)                             :: TmpData1                             ! A temporary           value  read in from a WAMIT file (-      )
      REAL(ReKi)                             :: TmpData2                             ! A temporary           value  read in from a WAMIT file (-      )
      REAL(ReKi)                             :: TmpDir                               ! A temporary direction        read in from a WAMIT file (degrees)
      REAL(ReKi)                             :: TmpIm                                ! A temporary imaginary value  read in from a WAMIT file (-      ) - stored as a REAL value
      REAL(ReKi)                             :: TmpPer                               ! A temporary period           read in from a WAMIT file (sec    )
      REAL(ReKi)                             :: TmpRe                                ! A temporary real      value  read in from a WAMIT file (-      )
      REAL(SiKi)                             :: TmpCoord(2)                          ! A temporary real array to hold the (Omega,WaveDir) pair for interpolation
      COMPLEX(SiKi),ALLOCATABLE              :: tmpComplexArr(:)                     ! A temporary array (0:NStepWave2-1) for FFT use. 
      REAL(ReKi), ALLOCATABLE                :: WAMITFreq (:)                        ! Frequency      components as ordered in the WAMIT output files (rad/s  )
      REAL(ReKi), ALLOCATABLE                :: WAMITPer  (:)                        ! Period         components as ordered in the WAMIT output files (sec    )
      REAL(ReKi), ALLOCATABLE                :: WAMITWvDir(:)                        ! Wave direction components as ordered in the WAMIT output files (degrees)

      INTEGER                                :: I,iGrid,iX,iY                        ! Generic index
      INTEGER                                :: InsertInd                            ! The lowest sorted index whose associated frequency component is higher than the current frequency component -- this is to sort the frequency components from lowest to highest
      INTEGER                                :: J                                    ! Generic index
      INTEGER                                :: K                                    ! Generic index
      INTEGER                                :: LastInd                              ! Index into the arrays saved from the last call as a starting point for this call
      INTEGER                                :: LastInd2(2)                          ! Index into the arrays saved from the last call as a starting point for this call. 2D
      INTEGER                                :: NInpFreq                             ! Number of input frequency components inherent in the hydrodynamic added mass matrix, hydrodynamic daming matrix, and complex wave excitation force per unit wave amplitude vector (-)
      INTEGER                                :: NInpWvDir                            ! Number of input incident wave propagation heading direction components inherent in the complex wave excitation force per unit wave amplitude vector (-)
      INTEGER,    ALLOCATABLE                :: SortFreqInd (:)                      ! The array of indices such that WAMITFreq (SortFreqInd (:)) is sorted from lowest to highest frequency (-)
      INTEGER,    ALLOCATABLE                :: SortWvDirInd(:)                      ! The array of indices such that WAMITWvDir(SortWvDirInd(:)) is sorted from lowest to highest agnle     (-)
      INTEGER                                :: Sttus                                ! Status returned by an attempted allocation or READ.
      INTEGER                                :: UnW1                                 ! I/O unit number for the WAMIT output file with the .1   extension; this file contains the linear, nondimensionalized, frequency-dependent solution to the radiation   problem.
      INTEGER                                :: UnW3                                 ! I/O unit number for the WAMIT output file with the .3   extension; this file contains the linear, nondimensionalized, frequency-dependent solution to the diffraction problem.
      INTEGER                                :: UnWh                                 ! I/O unit number for the WAMIT output file with the .hst extension; this file contains the linear, nondimensionalized hydrostatic restoring matrix.

      LOGICAL                                :: FirstFreq                            ! When .TRUE., indicates we're still looping through the first frequency component.
      LOGICAL                                :: FirstPass                            ! When .TRUE., indicates we're on the first pass through a loop.
      LOGICAL                                :: InfFreq                              ! When .TRUE., indicates that the infinite-frequency limit of added mass is contained within the WAMIT output files.
      LOGICAL                                :: NewPer                               ! When .TRUE., indicates that the period has just changed.
      LOGICAL                                :: ZeroFreq                             ! When .TRUE., indicates that the zero    -frequency limit of added mass is contained within the WAMIT output files.
      
      CHARACTER(MaxFileInfoLineLen)          :: Line                                 ! String to temporarily hold the value of a line within a WAMIT output file.

      TYPE(FFT_DataType)                     :: FFT_Data                             ! the instance of the FFT module we're using
      integer(IntKi)                         :: iSub, jSub                           ! indices into the 6x6 sub-matrices used to redimensionalize the WAMIT data (Needed because NBodyMod=1 could have WAMIT matrices which are 6N x 6N)
      integer(IntKi)                         :: iBody                                ! WAMIT body index
      real(R8Ki)                             :: orientation(3,3)                     ! Initial orientation of the WAMIT body 
      real(R8Ki)                             :: theta(3)                             ! Euler angle rotations of the WAMIT body
      real(ReKi)                             :: WaveNmbr                             ! Frequency-dependent wave number
      COMPLEX(SiKi)                          :: Fxy                                  ! Phase correction term for Wave excitation forces
      real(ReKi)                             :: tmpAngle                             ! Frequency and heading and platform offset dependent phase shift angle for Euler's Equation e^(-j*tmpAngle)
      COMPLEX(SiKi), ALLOCATABLE             :: HdroExctn_Local (:,:,:)              ! Temporary Frequency- and direction-dependent complex hydrodynamic wave excitation force per unit wave amplitude vector (kg/s^2, kg-m/s^2)
         ! Error handling
      CHARACTER(ErrMsgLen)                   :: ErrMsg2                              ! Temporary error message for calls
      INTEGER(IntKi)                         :: ErrStat2                             ! Temporary error status for calls
      COMPLEX(SiKi)                          :: Ctmp1, Ctmp2, Ctmp4, Ctmp5           ! Temporary COMPLEX transformation terms
      character(*), parameter                :: RoutineName = 'WAMIT_Init'


         ! Initialize data
         
      HighFreq    = 0.0
      UnW1        = 31
      UnW3        = 32
      UnWh        = 33
      LastInd     = 1
      LastInd2    = 0
      InfFreq     = .FALSE.    
      ZeroFreq    = .FALSE.
      
      
         ! Initialize ErrStat
         
      ErrStat  = ErrID_None         
      ErrMsg   = ""               
      
         ! Copy Output Init data from Waves Module Init call
         
      p%NStepWave    = InitInp%NStepWave
      p%ExctnMod     = InitInp%ExctnMod
      p%ExctnDisp    = InitInp%ExctnDisp
      p%ExctnCutOff  = InitInp%ExctnCutOff
      p%NBodyMod     = InitInp%NBodyMod
      p%NBody        = InitInp%NBody            ! In the context of this WAMIT object NBody is 1 if NBodyMod > 1 [there are NBody different WAMIT objects in this case]
      p%seast_interp_p = InitInp%seast_interp_p
         ! This module's implementation requires that if NBodyMod = 2 or 3, then there is one instance of a WAMIT module for each body, therefore, HydroDyn may have NBody > 1, but this WAMIT module will have NBody = 1
      if ( (p%NBodyMod > 1) .and. (p%NBody > 1) ) then
         CALL SetErrStat( ErrID_Fatal, "DEVELOPER ERROR: If NBodyMod = 2 or 3, then NBody for the a WAMIT object must be equal to 1", ErrStat, ErrMsg, RoutineName) 
         return
      end if     
      
         ! Allocate misc var and parameter vectors/matrices
      call AllocAry( p%F_HS_Moment_Offset,  6, p%NBody, 'p%F_HS_Moment_Offset', ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      call AllocAry( m%F_HS              ,  6*p%NBody, 'm%F_HS'              , ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      call AllocAry( m%F_Waves1          ,  6*p%NBody, 'm%F_Waves1'          , ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      call AllocAry( m%F_Rdtn            ,  6*p%NBody, 'm%F_Rdtn'            , ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      call AllocAry( m%F_PtfmAM          ,  6*p%NBody, 'm%F_PtfmAM'          , ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      call AllocAry( p%HdroAdMsI, 6*p%NBody,6*p%NBody, 'p%HdroAdMsI'         , ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      call AllocAry( p%HdroSttc , 6*p%NBody,6*p%NBody, 'p%HdroSttc'          , ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      
      do iBody = 1, p%NBody     
         p%F_HS_Moment_Offset(1,iBody) = 0.0_ReKi
         p%F_HS_Moment_Offset(2,iBody) = 0.0_ReKi
         p%F_HS_Moment_Offset(3,iBody) =  InitInp%RhoXg*InitInp%PtfmVol0(iBody)                                             ! except for the hydrostatic buoyancy force from Archimede's Principle when the support platform is in its undisplaced position
         p%F_HS_Moment_Offset(4,iBody) =  InitInp%RhoXg*InitInp%PtfmVol0(iBody)*( InitInp%PtfmCOByt(iBody) - InitInp%PtfmRefyt(iBody)  )  ! and the moment about X due to the COB being offset from the local WAMIT reference point
         p%F_HS_Moment_Offset(5,iBody) = -InitInp%RhoXg*InitInp%PtfmVol0(iBody)*( InitInp%PtfmCOBxt(iBody) - InitInp%PtfmRefxt(iBody)  )  ! and the moment about Y due to the COB being offset from the localWAMIT reference point
         p%F_HS_Moment_Offset(6,iBody) = 0.0_ReKi
      end do 

         
         ! Tell our nice users what is about to happen that may take a while:

      CALL WrScr ( ' Reading in WAMIT output with root name "'//TRIM(InitInp%WAMITFile)//'".' )

         ! Let's set up the matrices used to redimensionalize the hydrodynamic data
         !   from WAMIT; all these matrices are symmetric and need to be used with
         !   element-by-element multiplication, instead of matrix-by-matrix
         !   multiplication:

      SttcDim(1,1) = InitInp%RhoXg  *InitInp%WAMITULEN**2  ! Force-translation
      SttcDim(1,4) = InitInp%RhoXg  *InitInp%WAMITULEN**3  ! Force-rotation/Moment-translation - Hydrostatic restoring
      SttcDim(4,4) = InitInp%RhoXg  *InitInp%WAMITULEN**4  ! Moment-rotation

      RdtnDim(1,1) = InitInp%WtrDens*InitInp%WAMITULEN**3  ! Force-translation
      RdtnDim(1,4) = InitInp%WtrDens*InitInp%WAMITULEN**4  ! Force-rotation/Moment-translation - Hydrodynamic added mass and damping
      RdtnDim(4,4) = InitInp%WtrDens*InitInp%WAMITULEN**5  ! Moment-rotation

      DffrctDim(1) = InitInp%RhoXg  *InitInp%WAMITULEN**2  ! Force-translation - Hydrodynamic wave excitation force
      DffrctDim(4) = InitInp%RhoXg  *InitInp%WAMITULEN**3  ! Moment-rotation

      DO I = 1,3     ! Loop through all force-translation elements (rows)

         DO J = 1,3  ! Loop through all force-translation elements (columns)

            SttcDim(I,J) = SttcDim(1,1)

            RdtnDim(I,J) = RdtnDim(1,1)

         END DO       ! J - All force-translation elements (columns)

         DffrctDim (I  ) = DffrctDim(1)

      END DO          ! I - All force-translation elements (rows)

      DO I = 1,3     ! Loop through all force-rotation/moment-translation elements (rows/columns)

         DO J = 4,6  ! Loop through all force-rotation/moment-translation elements (columns/rows)

            SttcDim(I,J) = SttcDim(1,4)
            SttcDim(J,I) = SttcDim(1,4)

            RdtnDim(I,J) = RdtnDim(1,4)
            RdtnDim(J,I) = RdtnDim(1,4)

         END DO       ! J - All force-rotation/moment-translation elements (rows/columns)

      END DO          ! I - All force-rotation/moment-translation elements (columns/rows)

      DO I = 4,6     ! Loop through all moment-rotation elements (rows)

         DO J = 4,6  ! Loop through all moment-rotation elements (columns)

            SttcDim(I,J) = SttcDim(4,4)

            RdtnDim(I,J) = RdtnDim(4,4)

         END DO       ! J - All moment-rotation elements (columns)

         DffrctDim (I  ) = DffrctDim(4)

      END DO          ! I - All moment-rotation elements (rows)




         ! Let's read in and redimensionalize the hydrodynamic data from the WAMIT
         !   output files:



         ! Linear restoring from the hydrostatics problem:

      CALL OpenFInpFile ( UnWh, TRIM(InitInp%WAMITFile)//'.hst', ErrStat2, ErrMsg2 )  ! Open file.
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
         IF ( ErrStat >= AbortErrLev )  THEN
            CALL Cleanup()
            RETURN
         END IF
      p%HdroSttc (:,:) = 0.0 ! Initialize to zero

      DO    ! Loop through all rows in the file


         READ (UnWh,*,IOSTAT=Sttus)  I, J, TmpData1   ! Read in the row index, column index, and nondimensional data from the WAMIT file

         IF ( Sttus == 0 )  THEN                ! .TRUE. when data is read in successfully

         ! In case NBodyMod = 1, we now have WAMIT matrices which are potentially larger than 6x6, so we need to determine how the SttcDim multiplier matrix (a 6x6)
         !   should be applied to the larger WAMIT matrix.  
            
            iSub = mod(I-1,6)+1                         ! Finds the 6x6 sub-matrix indexing for the SttcDim multiplier matrix
            jSub = mod(J-1,6)+1  

            p%HdroSttc (I,J) = TmpData1*SttcDim(iSub,jSub)    ! Redimensionalize the data and place it at the appropriate location within the array

         ELSE                                           ! We must have reached the end of the file, so stop reading in data

            EXIT

         END IF


      END DO ! End loop through all rows in the file

      CLOSE ( UnWh ) ! Close file.

         ! need to transform p%HdroSttc when PtfmRefztRot is nonzero per plan
      call TransformWAMITMatrices( p%NBody, InitInp%PtfmRefztRot, p%HdroSttc )

         ! Linear, frequency-dependent hydrodynamic added mass and damping from the
         !   radiation problem:

      CALL OpenFInpFile ( UnW1, TRIM(InitInp%WAMITFile)//'.1', ErrStat2, ErrMsg2   )  ! Open file.
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
         IF ( ErrStat >= AbortErrLev )  THEN
            CALL Cleanup()
            RETURN
         END IF


         ! First find the number of input frequency components inherent in the
         !   hydrodynamic added mass matrix, hydrodynamic daming matrix, and complex
         !   wave excitation force per unit wave amplitude vector:

      NInpFreq  = 0        ! Initialize to zero
      PrvPer    = 0.0      ! Initialize to a don't care
      FirstPass = .TRUE.   ! Initialize to .TRUE. for the first pass

      DO    ! Loop through all rows in the file


         READ (UnW1,*,IOSTAT=Sttus)  TmpPer  ! Read in only the period from the WAMIT file

         IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully

            IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period currently read in is different than the previous period read in; thus we found a new frequency in the WAMIT file!
               NInpFreq  = NInpFreq + 1      ! Since we found a new frequency, count it in the total
               PrvPer    = TmpPer            ! Store the current period as the previous period for the next pass
               FirstPass = .FALSE.           ! Sorry, you can only have one first pass
            END IF

         ELSE                    ! We must have reached the end of the file, so stop reading in data  !bjj -- this isn't necessarially true....

            EXIT

         END IF


      END DO ! End loop through all rows in the file


      REWIND (UNIT=UnW1)   ! REWIND the file so we can read it in a second time.


      ! Now that we know how many frequencies there are, we can ALLOCATE the arrays
      !   to store the frequencies and frequency-dependent hydrodynamic added mass
      !   and damping matrices:

      CALL AllocAry( WAMITFreq,    NInpFreq,    'WAMITFreq',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL AllocAry( WAMITPer,     NInpFreq,    'WAMITPer',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL AllocAry( SortFreqInd,  NInpFreq,    'SortFreqInd',  ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL AllocAry( HdroFreq,     NInpFreq,    'HdroFreq',     ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL AllocAry( HdroAddMs,    NInpFreq, 6*p%NBody, 6*p%NBody, 'HdroAddMs',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL AllocAry( HdroDmpng,    NInpFreq, 6*p%NBody, 6*p%NBody, 'HdroDmpng',    ErrStat2, ErrMsg2 ); CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
         IF ( ErrStat >= AbortErrLev )  THEN
            CALL Cleanup()
            RETURN
         END IF



         ! Now find out how the frequencies are ordered in the file.  When we read in
         !   the added mass and damping matrices, we need to have them sorted by
         !   increasing frequency.  Thus, find the array of indices, SortFreqInd(),
         !   such that WAMITFreq(SortFreqInd(:)) is sorted from lowest to highest
         !   frequency:

      K         = 0        ! Initialize to zero
      PrvPer    = 0.0      ! Initialize to a don't care
      FirstPass = .TRUE.   ! Initialize to .TRUE. for the first pass

      DO    ! Loop through all rows in the file


         READ (UnW1,*,IOSTAT=Sttus)  TmpPer  ! Read in only the period from the WAMIT file

         IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully

            IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period currently read in is different than the previous period read in; thus we found a new frequency in the WAMIT file!

               K               = K + 1       ! This is current count of which frequency component we are on
               PrvPer          = TmpPer      ! Store the current period as the previous period for the next pass
               FirstPass       = .FALSE.     ! Sorry, you can only have one first pass

               WAMITPer    (K) = TmpPer         ! Store the periods                         in the order they appear in the WAMIT file
               IF (     TmpPer <  0.0 )  THEN   ! Periods less than zero in WAMIT represent infinite period = zero frequency
                  WAMITFreq(K) = 0.0
                  ZeroFreq     = .TRUE.
               ELSEIF ( TmpPer == 0.0 )  THEN   ! Periods equal to  zero in WAMIT represent infinite frequency
                  WAMITFreq(K) = HUGE(TmpPer)   ! Use HUGE() to approximate infinity in the precision of ReKi
                  InfFreq      = .TRUE.
               ELSE                             ! We must have positive, non-infinite frequency
                  WAMITFreq(K) = TwoPi/TmpPer   ! Store the periods as frequencies in rad/s in the order they appear in the WAMIT file
                  HighFreq     = MAX( HighFreq, WAMITFreq(K) ) ! Find the highest frequency (HighFreq) in the WAMIT output file, not counting infinity (even if the infinite frequency limit is in the file).
               END IF

               InsertInd       = K           ! Initialize as the K'th component
               DO I = 1,K-1   ! Loop throuh all previous frequencies
                  IF ( ( WAMITFreq(I) > WAMITFreq(K) ) )  THEN ! .TRUE. if a previous frequency component is higher than the current frequency component
                     InsertInd      = MIN( InsertInd, SortFreqInd(I) )  ! Store the lowest sorted index whose associated frequency component is higher than the current frequency component
                     SortFreqInd(I) = SortFreqInd(I) + 1                ! Shift all of the sorted indices up by 1 whose associated frequency component is higher than the current frequency component
                  END IF
               END DO          ! I - All previous frequencies
               SortFreqInd(K)  = InsertInd   ! Store the index such that WAMITFreq(SortFreqInd(:)) is sorted from lowest to highest frequency

            END IF

         ELSE                    ! We must have reached the end of the file, so stop reading in data

            EXIT

         END IF


      END DO ! End loop through all rows in the file


      REWIND (UNIT=UnW1)   ! REWIND the file so we can read it in a third time.  (This is getting ridiculous!)


         ! Now we can finally read in the frequency-dependent added mass and damping
         !   matrices; only store the upper-triangular portions (diagonal and above)
         !   of these matrices:

      K              = 0      ! Initialize to zero
      PrvPer         = 0.0    ! Initialize to a don't care
      FirstPass      = .TRUE. ! Initialize to .TRUE. for the first pass

      HdroAddMs = 0.0    ! Initialize to zero
      HdroDmpng = 0.0    ! Initialize to zero

      DO    ! Loop through all rows in the file


         READ (UnW1,'(A)',IOSTAT=Sttus)  Line   ! Read in the entire line

         IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully


            READ (Line,*)  TmpPer               ! Read in only the period from the WAMIT file


            IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period currently read in is different than the previous period read in; thus we found a new frequency in the WAMIT file!

               K              = K + 1           ! This is current count of which frequency component we are on
               PrvPer         = TmpPer          ! Store the current period as the previous period for the next pass
               FirstPass      = .FALSE.         ! Sorry, you can only have one first pass

               IF (     TmpPer <  0.0 )  THEN   ! Periods less than zero in WAMIT represent infinite period = zero frequency
                  HdroFreq (SortFreqInd(K)) = 0.0
               ELSEIF ( TmpPer == 0.0 )  THEN   ! Periods equal to  zero in WAMIT represent infinite frequency; a value slightly larger than HighFreq is returned to approximate infinity while still maintaining an effective interpolation later on.
                  HdroFreq (SortFreqInd(K)) = HighFreq*OnePlusEps ! Set the infinite frequency to a value slightly larger than HighFreq
               ELSE                             ! We must have positive, non-infinite frequency
                  HdroFreq (SortFreqInd(K)) = TwoPi/TmpPer  ! Convert the period in seconds to a frequency in rad/s and store them sorted from lowest to highest
               END IF

            END IF


            IF ( TmpPer <= 0.0 )  THEN          ! .TRUE. if the current period is less than or equal to zero, which in WAMIT represents the zero and infinite frequency limits, respectively; in these cases, only the added mass matrix is computed and output by WAMIT (and based on hydrodynamic theory, the damping matrix is zero as initialized above)

               READ (Line,*,IOSTAT=Sttus)  TmpPer, I, J, TmpData1           ! Read in the period, row index, column index, and nondimensional data from the WAMIT file

               IF ( Sttus /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, "Error reading line from WAMIT file", ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN
               END IF              
!bjj: verify that I and J are valid indices for RdtnDim                  
                  
                  
              !  IF ( J >= I )  THEN  ! .TRUE. if we are on or above the diagonal
              !    Indx = 6*( I - 1 ) + J - ( I*( I - 1 ) )/2                                       ! Convert from row/column indices to an index in the format used to save only the upper-triangular portion of the matrix.  NOTE: ( I*( I - 1 ) )/2 = SUM(I,START=1,END=I-1).
                  iSub = mod(I-1,6)+1                                                              ! Finds the 6x6 sub-matrix indexing for the SttcDim multiplier matrix
                  jSub = mod(J-1,6)+1  

                  HdroAddMs(SortFreqInd(K),I,J) = TmpData1*RdtnDim(iSub,jSub)                     ! Redimensionalize the data and place it at the appropriate location within the array
              ! END IF

            ELSE                                ! We must have a positive, non-infinite frequency.

               READ (Line,*,IOSTAT=Sttus)  TmpPer, I, J, TmpData1, TmpData2                        ! Read in the period, row index, column index, and nondimensional data from the WAMIT file
               IF ( Sttus /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, "Error reading line from WAMIT file", ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN
               END IF              
                  

               !IF ( J >= I )  THEN  ! .TRUE. if we are on or above the diagonal
                !  Indx = 6*( I - 1 ) + J - ( I*( I - 1 ) )/2                                             ! Convert from row/column indices to an index in the format used to save only the upper-triangular portion of the matrix.  NOTE: ( I*( I - 1 ) )/2 = SUM(I,START=1,END=I-1).
                  iSub = mod(I-1,6)+1                                                                    ! Finds the 6x6 sub-matrix indexing for the SttcDim multiplier matrix
                  jSub = mod(J-1,6)+1  
                  HdroAddMs(SortFreqInd(K),I,J) = TmpData1*RdtnDim(iSub,jSub)                           ! Redimensionalize the data and place it at the appropriate location within the array
                  HdroDmpng(SortFreqInd(K),I,J) = TmpData2*RdtnDim(iSub,jSub)*HdroFreq(SortFreqInd(K))  ! Redimensionalize the data and place it at the appropriate location within the array
              ! END IF

            END IF


         ELSE                    ! We must have reached the end of the file, so stop reading in data


            EXIT


         END IF


      END DO ! End loop through all rows in the file


      CLOSE ( UnW1 ) ! Close file.

         ! need to transform HdroAddMs and HdroDmpng when PtfmRefztRot is nonzero per plan
      do I = 1, NInpFreq
         call TransformWAMITMatrices( p%NBody, InitInp%PtfmRefztRot, HdroAddMs(I,:,:) )
         call TransformWAMITMatrices( p%NBody, InitInp%PtfmRefztRot, HdroDmpng(I,:,:) )
      end do
      
      
         ! Linear, frequency- and direction-dependent complex hydrodynamic wave
         !   excitation force per unit wave amplitude vector from the diffraction
         !   problem:

      CALL OpenFInpFile ( UnW3, TRIM(InitInp%WAMITFile)//'.3', ErrStat2, ErrMsg2   )  ! Open file.
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
         IF ( ErrStat >= AbortErrLev )  THEN
            CALL Cleanup()
            RETURN
         END IF
            

         ! First find the number of input incident wave propagation heading direction
         !   components inherent in the complex wave excitation force per unit wave
         !   amplitude vector:

      NInpWvDir = 0        ! Initialize to zero
      PrvDir    = 0.0      ! Initialize to a don't care
      FirstPass = .TRUE.   ! Initialize to .TRUE. for the first pass

      DO    ! Loop through all rows in the file


         READ (UnW3,'(A)',IOSTAT=Sttus)  Line   ! Read in the entire line

         IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully


            READ (Line,*)  TmpPer, TmpDir ! Read in only the period and direction from the WAMIT file !bjj why don't we check IOSTAT here, too????


            IF ( FirstPass                           )  THEN   ! .TRUE. if we are on the first pass
               PrvPer = TmpPer            ! Store the current period    as the previous period    for the next pass
            END IF


            IF (                  TmpPer /= PrvPer   )  THEN   ! .TRUE.                                if the period    currently read in is different than the previous period    read in; thus we found a new period    in the WAMIT file, so stop reading in data
               EXIT
            END IF


            IF ( FirstPass .OR. ( TmpDir /= PrvDir ) )  THEN   ! .TRUE. if we are on the first pass or if the direction currently read in is different than the previous direction read in; thus we found a new direction in the WAMIT file!
               NInpWvDir = NInpWvDir + 1  ! Since we found a new direction, count it in the total
               PrvDir    = TmpDir         ! Store the current direction as the previous direction for the next pass
               FirstPass = .FALSE.        ! Sorry, you can only have one first pass
            END IF


         ELSE                    ! We must have reached the end of the file, so stop reading in data


            EXIT


         END IF


      END DO ! End loop through all rows in the file


      REWIND (UNIT=UnW3)   ! REWIND the file so we can read it in a second time.


      ! Now that we know how many directions there are, we can ALLOCATE the arrays to
      !   to store the directions and frequency- and direction-dependent complex wave
      !   excitation force per unit wave amplitude vector:

      CALL AllocAry(  WAMITWvDir,    NInpWvDir, 'WAMITWvDir',   ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
      CALL AllocAry(  SortWvDirInd,  NInpWvDir, 'SortWvDirInd', ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
      CALL AllocAry(  HdroWvDir,     NInpWvDir, 'HdroWvDir',    ErrStat2, ErrMsg2 );  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
         IF ( ErrStat >= AbortErrLev )  THEN
            CALL Cleanup()
            RETURN
         END IF
            
      ALLOCATE ( HdroExctn   (NInpFreq,NInpWvDir,6*p%NBody) , STAT=ErrStat2 ) ! complex so we don't have a built in subroutine
      IF ( ErrStat2 /= 0 )  THEN
         CALL SetErrStat( ErrID_Fatal, 'Error allocating space for HdroExctn array', ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
         RETURN
      END IF


         ! Now find out how the directions are ordered in the file.  When we read in
         !   the wave excitation force vector, we need to have them sorted by
         !   increasing angle.  Thus, find the array of indices, SortWvDirInd(),
         !   such that WAMITWvDir(SortWvDirInd(:)) is sorted from lowest to highest
         !   angle.  At the same time, make sure that the frequencies in the .3 file are
         !   ordered in the same way they are in the .1 file and make sure that the
         !   directions are the same for each frequency component:

      K         = 0        ! Initialize to zero
      PrvPer    = 0.0      ! Initialize to a don't care
      PrvDir    = 0.0      ! Initialize to a don't care
      FirstPass = .TRUE.   ! Initialize to .TRUE. for the first pass

      DO    ! Loop through all rows in the file


         READ (UnW3,'(A)',IOSTAT=Sttus)  Line   ! Read in the entire line

         IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully


            READ (Line,*,IOSTAT=Sttus)  TmpPer, TmpDir ! Read in only the period and direction from the WAMIT file
               IF ( Sttus /= 0 )  THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error reading period and direction from WAMIT file.', ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN
               END IF


            IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period    currently read in is different than the previous period    read in; thus we found a new period    in the WAMIT file!

               J         = 0           ! Reset the count of directions to zero
               K         = K + 1       ! This is current count of which frequency component we are on
               PrvPer    = TmpPer      ! Store the current period    as the previous period    for the next pass
               FirstFreq = FirstPass   ! Sorry, you can only loop through the first frequency once
               NewPer    = .TRUE.      ! Reset the new period flag

               DO WHILE ( WAMITPer(K) <= 0.0 )  ! Periods less than or equal to zero in WAMIT represent infinite period = zero frequency and infinite frequency, respectively.  However, only the added mass is output by WAMIT at these limits.  The damping and wave excitation are left blank, so skip them!
                  K = K + 1
               END DO

               IF ( TmpPer /= WAMITPer(K) )  THEN  ! Abort if the .3 and .1 files do not contain the same frequency components (not counting zero and infinity)
                  ErrMsg2  = ' Other than zero and infinite frequencies, "'   //TRIM(InitInp%WAMITFile)//'.3",' // &
                               ' contains different frequency components than "'//TRIM(InitInp%WAMITFile)//'.1". '// &
                               ' Both WAMIT output files must be generated from the same run.'
                  CALL SetErrStat( ErrID_Fatal, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN
               END IF

            END IF


            IF ( FirstPass .OR. ( TmpDir /= PrvDir ) .OR. NewPer )  THEN   ! .TRUE. if we are on the first pass, or if this is new period, or if the direction currently read in is different than the previous direction read in; thus we found a new direction in the WAMIT file!

               J         = J + 1       ! This is current count of which direction component we are on
               PrvDir    = TmpDir      ! Store the current direction as the previous direction for the next pass
               FirstPass = .FALSE.     ! Sorry, you can only have one first pass
               NewPer    = .FALSE.     ! Disable the new period flag

               IF ( FirstFreq )  THEN                    ! .TRUE. while we are still looping through all directions for the first frequency component
                  WAMITWvDir(J)   = TmpDir      ! Store the directions in the order they appear in the WAMIT file

                  InsertInd       = J           ! Initialize as the J'th component
                  DO I = 1,J-1   ! Loop throuh all previous directions
                     IF ( ( WAMITWvDir(I) > WAMITWvDir(J) ) )  THEN  ! .TRUE. if a previous direction component is higher than the current direction component
                        InsertInd       = MIN( InsertInd, SortWvDirInd(I) )   ! Store the lowest sorted index whose associated direction component is higher than the current direction component
                        SortWvDirInd(I) = SortWvDirInd(I) + 1                 ! Shift all of the sorted indices up by 1 whose associated direction component is higher than the current direction component
                     END IF
                  END DO          ! I - All previous directions
                  SortWvDirInd(J) = InsertInd   ! Store the index such that WAMITWvDir(SortWvDirInd(:)) is sorted from lowest to highest direction
               ELSEIF ( TmpDir /= WAMITWvDir(J) )  THEN  ! We must have looped through all directions at least once; so check to make sure all subsequent directions are consistent with the directions from the first frequency component, otherwise Abort
                  ErrMsg2  = ' Not every frequency component in "'//TRIM(InitInp%WAMITFile)//'.3"'// &
                               ' contains the same listing of direction angles.  Check for' // &
                               ' errors in the WAMIT output file.'
                  CALL SetErrStat( ErrID_Fatal, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN
               END IF

            END IF


         ELSE                    ! We must have reached the end of the file, so stop reading in data


            EXIT


         END IF


      END DO ! End loop through all rows in the file

if (p%ExctnMod == 1 ) then
      REWIND (UNIT=UnW3)   ! REWIND the file so we can read it in a third time.  (This is getting ridiculous!)


         ! Now we can finally read in the frequency- and direction-dependent complex
         !   wave excitation force per unit wave amplitude vector:

      K                = 0       ! Initialize to zero
      PrvPer           = 0.0     ! Initialize to a don't care
      PrvDir           = 0.0     ! Initialize to a don't care
      FirstPass        = .TRUE.  ! Initialize to .TRUE. for the first pass

      HdroExctn(:,:,:) = 0.0     ! Initialize to zero

      DO    ! Loop through all rows in the file


         READ (UnW3,'(A)',IOSTAT=Sttus)  Line   ! Read in the entire line

         IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully


            READ (Line,*,IOSTAT=Sttus)  TmpPer, TmpDir, I, TmpData1, TmpData2, TmpRe, TmpIm   ! Read in the period, direction, row index, and nondimensional data from the WAMIT file
               IF ( Sttus /= 0 )  THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error reading period and direction, row index, and nondimensional data from the WAMIT file.', ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN
               END IF


            IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period    currently read in is different than the previous period    read in; thus we found a new period    in the WAMIT file!

               J            = 0           ! Reset the count of directions to zero
               K            = K + 1       ! This is current count of which frequency component we are on
               PrvPer       = TmpPer      ! Store the current period    as the previous period    for the next pass
               FirstFreq    = FirstPass   ! Sorry, you can only loop through the first frequency once
               NewPer       = .TRUE.      ! Reset the new period flag

               DO WHILE ( WAMITPer(K) <= 0.0 )  ! Periods less than or equal to zero in WAMIT represent infinite period = zero frequency and infinite frequency, respectively.  However, only the added mass is output by WAMIT at these limits.  The damping and wave excitation are left blank, so skip them!
                  K = K + 1
               END DO

            END IF


            IF ( FirstPass .OR. ( TmpDir /= PrvDir ) .OR. NewPer )  THEN   ! .TRUE. if we are on the first pass, or if this is new period, or if the direction currently read in is different than the previous direction read in; thus we found a new direction in the WAMIT file!

               J            = J + 1       ! This is current count of which direction component we are on
               PrvDir       = TmpDir      ! Store the current direction as the previous direction for the next pass
               FirstPass    = .FALSE.     ! Sorry, you can only have one first pass
               NewPer       = .FALSE.     ! Disable the new period flag

               IF ( FirstFreq )  THEN  ! .TRUE. while we are still looping through all directions for the first frequency component
                  HdroWvDir(SortWvDirInd(J)) = TmpDir ! Store the directions sorted from lowest to highest
               END IF

            END IF

            iSub = mod(I-1,6)+1                                                                    ! Finds the 6x6 sub-matrix indexing for the SttcDim multiplier matrix 
            HdroExctn(SortFreqInd(K),SortWvDirInd(J),I) = CMPLX( TmpRe, TmpIm )*DffrctDim(iSub) ! Redimensionalize the data and place it at the appropriate location within the array


         ELSE                    ! We must have reached the end of the file, so stop reading in data


            EXIT


         END IF


      END DO ! End loop through all rows in the file
end if
      CLOSE ( UnW3 ) ! Close file.

      ! For some reason, WAMIT computes the zero- and infinite- frequency limits for
      !   only the added mass.  Based on hydrodynamic theory, the damping is zero at
      !   these limits (as initialized).  Hydrodynamic theory also says that the
      !   infinite-frequency limit of the diffraction force is zero (as initialized);
      !   however, the zero-frequency limit need not be zero.  Thus, if necessary
      !   (i.e., if we have read in a WAMIT output file that contains the
      !   zero-frequency limit of the added mass), compute the zero-frequency limit
      !   of the diffraction problem using the known values at the lowest
      !   nonzero-valued frequency available:

      DO I = 1,NInpFreq       ! Loop through all input frequency components

         IF ( HdroFreq(I) > 0.0 )  THEN ! .TRUE. at the lowest nonzero-valued frequency component

            DO J = I-1,1,-1   ! Loop through all zero-valued frequency components
               HdroExctn(J,:,:) = HdroExctn(I,:,:) ! Set the zero-frequency limits to equal the known values at the lowest nonzero-valued frequency available
            END DO             ! J - All zero-valued frequency components

            EXIT  ! Since HdroFreq(:) is sorted from lowest to highest frequency, there is no reason to continue on once we have found the lowest nonzero-valued frequency component

         END IF

      END DO                   ! I - All input frequency components


      
            ! Tell our nice users what is about to happen that may take a while:

      CALL WrScr ( ' Computing radiation impulse response functions and wave diffraction forces.' )



         ! Abort if the WAMIT files do not contain both the zero- and and infinite-
         !   frequency limits of added mass.  

      IF ( .NOT. ( ZeroFreq .AND. InfFreq ) )  THEN   ! .TRUE. if both the zero- and infinite-frequency limits of added mass are contained within the WAMIT file
         ErrMsg2  = ' "'//TRIM(InitInp%WAMITFile)// &
                          '.1" must contain both the zero- and infinite-frequency limits of added mass.'
         CALL SetErrStat( ErrID_Fatal, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
         RETURN         
      END IF



         ! Set the infinite-frequency limit of the frequency-dependent hydrodynamic
         !   added mass matrix, HdroAdMsI, based on the highest frequency available:
!TODO: Is this index order correct for computational speed? GJH 9/5/19
      !Indx = 0
      DO J = 1,6*p%NBody        ! Loop through all rows    of HdroAdMsI
         DO K = 1,6*p%NBody     ! Loop through all columns of HdroAdMsI above and including the diagonal
           ! Indx = Indx + 1
            p%HdroAdMsI(J,K) = HdroAddMs(NInpFreq,J,K)
         END DO          ! K - All columns of HdroAdMsI above and including the diagonal
      !   DO K = J+1,6   ! Loop through all rows    of HdroAdMsI below the diagonal
      !      p%HdroAdMsI(K,J) = p%HdroAdMsI(J,K)
      !   END DO          ! K - All rows    of HdroAdMsI below the diagonal
      END DO             ! J - All rows    of HdroAdMsI


      if ( ( p%ExctnMod == 0 ) ) then
         
         ! no need to allocate the p%WaveExctn array because it won't be used
         
      else
           ! Initialize the variables associated with the incident wave:

         SELECT CASE ( InitInp%WaveMod ) ! Which incident wave kinematics model are we using?
         CASE ( 0 )  ! No waves, NOTE: for this case we are forcing ExctnDisp = 0, so only p%WaveExctn needs to be allocated, not p%WaveExctnGrid
            if ( p%ExctnMod == 1 ) then

                  ! Initialize everything to zero:

               ALLOCATE ( p%WaveExctn (0:InitInp%NStepWave,6*p%NBody) , STAT=ErrStat2 )
               IF ( ErrStat2 /= 0 )  THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the WaveExctn array.', ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN
               END IF

               p%WaveExctn = 0.0   
         
            else if ( p%ExctnMod == 2 ) then
               Interval_Sub                  = InitInp%Conv_Rdtn%RdtnDT
               SS_Exctn_InitInp%InputFile    = InitInp%WAMITFile    
               SS_Exctn_InitInp%WaveDir      = InitInp%WaveDir
               SS_Exctn_InitInp%NStepWave    = p%NStepWave
               SS_Exctn_InitInp%NBody        = InitInp%NBody
               SS_Exctn_InitInp%PtfmRefztRot = InitInp%PtfmRefztRot 
               SS_Exctn_InitInp%ExctnDisp    = InitInp%ExctnDisp
               
                  ! No other modules need this WaveElev0 array so we will simply move the allocation over to the SS_Exctn module
               IF (ALLOCATED(InitInp%WaveElev0)) CALL MOVE_ALLOC(InitInp%WaveElev0, SS_Exctn_InitInp%WaveElev0) 
               !SS_Exctn_InitInp%WaveElev0 => InitInp%WaveElev0
               SS_Exctn_InitInp%WaveElev1 => InitInp%WaveElev1
!TODO: Verify what happens within SS_Exctn when we have no waves. 
               
                  ! We need the WaveTime array to stay intact for use in other modules, so we will make a copy instead of moving the allocation
               !ALLOCATE ( SS_Exctn_InitInp%WaveTime (0:InitInp%NStepWave) , STAT=ErrStat2 )
               !IF ( ErrStat2 /= 0 )  THEN
               !   CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the SS_Exctn_InitInp%WaveTime array.', ErrStat, ErrMsg, RoutineName)
               !   CALL Cleanup()
               !   RETURN            
               !END IF
               SS_Exctn_InitInp%WaveTime => InitInp%WaveTime 
               
               call SS_Exc_Init(SS_Exctn_InitInp, m%SS_Exctn_u, p%SS_Exctn, x%SS_Exctn, xd%SS_Exctn, z%SS_Exctn, OtherState%SS_Exctn, &
                                      m%SS_Exctn_y, m%SS_Exctn, Interval_Sub, SS_Exctn_InitOut, ErrStat2, ErrMsg2)
            
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  if ( ErrStat >= AbortErrLev ) then
                     call Cleanup()
                     return
                  end if   
            end if
         CASE ( 1, 2, 3, 4, 5, 7, 10 )    ! Plane progressive (regular) wave, JONSWAP/Pierson-Moskowitz spectrum (irregular) wave, white-noise wave,  or user-defined spectrum (irregular) wave.


            if ( p%ExctnMod == 1 ) then
               
               ! Abort if we have chosen a wave heading direction that is outside the range
               !   of directions where the complex wave excitation force per unit wave
               !   amplitude vector has been defined, else interpolate to find the complex
               !   wave excitation force per unit wave amplitude vector at the chosen wave
               !   heading direction:
               ! NOTE: we may end up inadvertantly aborting if the wave direction crosses
               !   the -Pi / Pi boundary (-180/180 degrees).

               IF ( ( InitInp%WaveDirMin < HdroWvDir(1) ) .OR. ( InitInp%WaveDirMax > HdroWvDir(NInpWvDir) ) )  THEN
                  ErrMsg2  = 'All Wave directions must be within the wave heading angle range available in "' &
                                 //TRIM(InitInp%WAMITFile)//'.3" (inclusive).'
                  CALL SetErrStat( ErrID_Fatal, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN
               END IF

                  ! Calculate the WaveExctn data from WAMIT data if ExctnMod = 1
               
                  ! ALLOCATE the arrays:

               ALLOCATE ( WaveExctnC(0:InitInp%NStepWave2 ,6*p%NBody) , STAT=ErrStat2 )
               IF ( ErrStat2 /= 0 )  THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the WaveExctnC array.', ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN            
               END IF

               if (p%ExctnDisp > 0 ) then
                  ALLOCATE ( WaveExctnCGrid(0:InitInp%NStepWave2 ,p%SeaSt_Interp_p%n(2)*p%SeaSt_Interp_p%n(3),6*p%NBody) , STAT=ErrStat2 )
                  IF ( ErrStat2 /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the WaveExctnC array.', ErrStat, ErrMsg, RoutineName)
                     CALL Cleanup()
                     RETURN            
                  END IF
                  ALLOCATE ( p%WaveExctnGrid (0:InitInp%NStepWave,p%SeaSt_Interp_p%n(2),p%SeaSt_Interp_p%n(3), 6*p%NBody) , STAT=ErrStat2 )
                  IF ( ErrStat2 /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the WaveExctn array.', ErrStat, ErrMsg, RoutineName)
                     CALL Cleanup()
                     RETURN            
                  END IF
               else
                  ALLOCATE ( p%WaveExctn (0:InitInp%NStepWave,6*p%NBody) , STAT=ErrStat2 )
                  IF ( ErrStat2 /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the WaveExctn array.', ErrStat, ErrMsg, RoutineName)
                     CALL Cleanup()
                     RETURN            
                  END IF
               end if
              

               !====================================
               ! Transform the wave excitation coefs
               !====================================
               
               if ( p%NBodyMod == 2 ) then
                  
                  ! Since NBodyMod = 2, then NBody = 1 for this WAMIT object (this requirement is encoded at the HydroDyn module level)
                  
                  allocate (  HdroExctn_Local(NInpFreq, NInpWvDir, 6),   STAT=ErrStat2 )
                  IF ( ErrStat2 /= 0 )  THEN
                     CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the HdroExctn_Local array.', ErrStat, ErrMsg, RoutineName)
                     CALL Cleanup()
                     RETURN            
                  END IF

                  do K = 1,6           ! Loop through all wave excitation forces and moments
                     do J = 1, NInpWvDir                       
                        TmpCoord(2) = HdroWvDir(J) - InitInp%PtfmRefztRot(1)*R2D  ! apply locale Z rotation to heading angle (degrees)
                        do I = 1, NInpFreq
                           TmpCoord(1) = HdroFreq(I)
                           ! Iterpolate to find new coef
                           call WAMIT_Interp2D_Cplx( TmpCoord, HdroExctn(:,:,K), HdroFreq, HdroWvDir, LastInd2, HdroExctn_Local(I,J,K), ErrStat2, ErrMsg2 )
                        end do
                     end do
                  end do

                  ! Now apply rotation and phase shift 

                  do J = 1, NInpWvDir  
                     do I = 1, NInpFreq
                           ! Fxy = exp(-j * k(w) * ( X*cos(Beta(w)) + Y*sin(Beta(w)) )
                        WaveNmbr   = WaveNumber ( HdroFreq(I), InitInp%Gravity, InitInp%WtrDpth )
                        tmpAngle   = WaveNmbr * ( InitInp%PtfmRefxt(1)*cos(HdroWvDir(J)*D2R) + InitInp%PtfmRefyt(1)*sin(HdroWvDir(J)*D2R) )
                        TmpRe =  cos(tmpAngle)
                        TmpIm = -sin(tmpAngle)
                        Fxy   = CMPLX( TmpRe, TmpIm )

                        HdroExctn(I,J,1) = Fxy*( HdroExctn_Local(I,J,1)*cos(InitInp%PtfmRefztRot(1)) -  HdroExctn_Local(I,J,2)*sin(InitInp%PtfmRefztRot(1)) )
                        HdroExctn(I,J,2) = Fxy*( HdroExctn_Local(I,J,1)*sin(InitInp%PtfmRefztRot(1)) +  HdroExctn_Local(I,J,2)*cos(InitInp%PtfmRefztRot(1)) )
                        HdroExctn(I,J,3) = Fxy*( HdroExctn_Local(I,J,3) )
                        HdroExctn(I,J,4) = Fxy*( HdroExctn_Local(I,J,4)*cos(InitInp%PtfmRefztRot(1)) -  HdroExctn_Local(I,J,5)*sin(InitInp%PtfmRefztRot(1)) )
                        HdroExctn(I,J,5) = Fxy*( HdroExctn_Local(I,J,4)*sin(InitInp%PtfmRefztRot(1)) +  HdroExctn_Local(I,J,5)*cos(InitInp%PtfmRefztRot(1)) )
                        HdroExctn(I,J,6) = Fxy*( HdroExctn_Local(I,J,6) )

                     end do
                  end do  
                  deallocate(HdroExctn_Local)
               else
                  
                     ! Apply rotation only for NBodyMod = 1,3
                  do J = 1, NInpWvDir  
                     do I = 1, NInpFreq

                        Ctmp1 = ( HdroExctn(I,J,1)*cos(InitInp%PtfmRefztRot(1)) ) - ( HdroExctn(I,J,2)*sin(InitInp%PtfmRefztRot(1)) )
                        Ctmp2 = ( HdroExctn(I,J,1)*sin(InitInp%PtfmRefztRot(1)) ) + ( HdroExctn(I,J,2)*cos(InitInp%PtfmRefztRot(1)) )  
                        Ctmp4 = ( HdroExctn(I,J,4)*cos(InitInp%PtfmRefztRot(1)) ) - ( HdroExctn(I,J,5)*sin(InitInp%PtfmRefztRot(1)) )
                        Ctmp5 = ( HdroExctn(I,J,4)*sin(InitInp%PtfmRefztRot(1)) ) + ( HdroExctn(I,J,5)*cos(InitInp%PtfmRefztRot(1)) )

                        HdroExctn(I,J,1) = Ctmp1
                        HdroExctn(I,J,2) = Ctmp2
                        HdroExctn(I,J,4) = Ctmp4
                        HdroExctn(I,J,5) = Ctmp5
                     end do
                  end do  
                  
               end if
             if (p%ExctnDisp == 0 ) then  
               ! Compute the positive-frequency components (including zero) of the discrete
               !   Fourier transform of the wave excitation force:

               DO I = 0,InitInp%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transform

                     ! Compute the frequency of this component:

                  Omega = I*InitInp%WaveDOmega
           

                     ! Compute the discrete Fourier transform of the instantaneous value of the
                     !   total excitation force on the support platfrom from incident waves:

                  DO J = 1,6*p%NBody           ! Loop through all wave excitation forces and moments
                     TmpCoord(1) = Omega
                     TmpCoord(2) = InitInp%WaveDirArr(I)
                     CALL WAMIT_Interp2D_Cplx( TmpCoord, HdroExctn(:,:,J), HdroFreq, HdroWvDir, LastInd2, WaveExctnC(I,J), ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                     IF ( ErrStat >= AbortErrLev ) THEN
                        CALL Cleanup()
                        RETURN
                     END IF
                     WaveExctnC(I,J) = WaveExctnC(I,J) * CMPLX(InitInp%WaveElevC0(1,I), InitInp%WaveElevC0(2,I))
                  
                  END DO                ! J - All wave excitation forces and moments


               END DO                ! I - The positive frequency components (including zero) of the discrete Fourier transform
            ! Compute the inverse discrete Fourier transform to find the time-domain
                  !   representation of the wave excitation force:

               CALL InitFFT ( InitInp%NStepWave, FFT_Data, .TRUE., ErrStat2 )
                  CALL SetErrStat( ErrStat2, 'Error in call to InitFFT.', ErrStat, ErrMsg, RoutineName)
                  IF ( ErrStat >= AbortErrLev) THEN
                     CALL Cleanup()
                     RETURN
                  END IF
         
               DO J = 1,6*p%NBody           ! Loop through all wave excitation forces and moments
                  CALL ApplyFFT_cx ( p%WaveExctn(0:InitInp%NStepWave-1,J), WaveExctnC(:,J), FFT_Data, ErrStat2 )
                  CALL SetErrStat( ErrStat2, ' An error occured while applying an FFT to WaveExctnC.', ErrStat, ErrMsg, RoutineName)
                  IF ( ErrStat >= AbortErrLev) THEN
                     CALL Cleanup()
                     RETURN
                  END IF
            
                     ! Append first datpoint as the last as aid for repeated wave data
                  p%WaveExctn(InitInp%NStepWave,J) = p%WaveExctn(0,J)
               END DO                ! J - All wave excitation forces and moments

               CALL ExitFFT(FFT_Data, ErrStat2)
                  CALL SetErrStat( ErrStat2, 'Error in call to ExitFFT.', ErrStat, ErrMsg, RoutineName)
                  IF ( ErrStat >= AbortErrLev) THEN
                     CALL Cleanup()
                     RETURN
                  END IF

            else
               DO I = 0,InitInp%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transform

                     ! Compute the frequency of this component:

                  Omega = I*InitInp%WaveDOmega
               
                     ! Compute the discrete Fourier transform of the instantaneous value of the
                     !   total excitation force on the support platfrom from incident waves:

                  DO J = 1,6*p%NBody           ! Loop through all wave excitation forces and moments
                     TmpCoord(1) = Omega
                     TmpCoord(2) = InitInp%WaveDirArr(I)
                     CALL WAMIT_Interp2D_Cplx( TmpCoord, HdroExctn(:,:,J), HdroFreq, HdroWvDir, LastInd2, WaveExctnC(I,J), ErrStat2, ErrMsg2 )
                     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                     IF ( ErrStat >= AbortErrLev ) THEN
                        CALL Cleanup()
                        RETURN
                     END IF
                     do iGrid = 1, p%SeaSt_Interp_p%n(2)*p%SeaSt_Interp_p%n(3)
                        WaveExctnCGrid(I,iGrid,J) = WaveExctnC(I,J) * CMPLX(InitInp%WaveElevC(1,I,iGrid), InitInp%WaveElevC(2,I,iGrid))                        
                     end do
                  END DO                ! J - All wave excitation forces and moments
               END DO                ! I - The positive frequency components (including zero) of the discrete Fourier transform
               
               ! Compute the inverse discrete Fourier transform to find the time-domain
                  !   representation of the wave excitation force:

               CALL InitFFT ( InitInp%NStepWave, FFT_Data, .TRUE., ErrStat2 )
                  CALL SetErrStat( ErrStat2, 'Error in call to InitFFT.', ErrStat, ErrMsg, RoutineName)
                  IF ( ErrStat >= AbortErrLev) THEN
                     CALL Cleanup()
                     RETURN
                  END IF
         
               DO J = 1,6*p%NBody           ! Loop through all wave excitation forces and moments
                  do iGrid = 1, p%SeaSt_Interp_p%n(2)*p%SeaSt_Interp_p%n(3)
                        iX = mod(iGrid-1, p%SeaSt_Interp_p%n(2)) + 1  ! 1st n index is time
                        iY = (iGrid-1) / p%SeaSt_Interp_p%n(2) + 1
                        CALL ApplyFFT_cx ( p%WaveExctnGrid(0:InitInp%NStepWave-1,iX,iY,J), WaveExctnCGrid(:,iGrid,J), FFT_Data, ErrStat2 )
                        CALL SetErrStat( ErrStat2, ' An error occured while applying an FFT to WaveExctnC.', ErrStat, ErrMsg, RoutineName)
                        IF ( ErrStat >= AbortErrLev) THEN
                           CALL Cleanup()
                           RETURN
                        END IF
                           ! Append first datpoint as the last as aid for repeated wave data
                        p%WaveExctnGrid(InitInp%NStepWave,iX,iY,J) = p%WaveExctnGrid(0,iX,iY,J)
                  end do
                  
               END DO                ! J - All wave excitation forces and moments

               CALL ExitFFT(FFT_Data, ErrStat2)
                  CALL SetErrStat( ErrStat2, 'Error in call to ExitFFT.', ErrStat, ErrMsg, RoutineName)
                  IF ( ErrStat >= AbortErrLev) THEN
                     CALL Cleanup()
                     RETURN
                  END IF
                  
            end if
            

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Dump the HdroFreq variable to a file for debugging
      ! Open and write header info to the HydroDyn Output File
      !CALL OpenFOutFile ( 66, 'C:\Dev\NREL_SVN\HydroDyn\branches\HydroDyn_Modularization\Samples\NRELOffshrBsline5MW_OC3Hywind\HdroFreq_HD.txt', ErrStat   )  ! Open motion file.
      !DO K = 1, NInpFreq
      !   WRITE ( 66, '(2(e20.9))', IOSTAT = ErrStat) REAL(K), HdroFreq(K)
      !END DO
      !CLOSE ( 66 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Dump the WaveElevCO variable to a file for debugging
      ! Open and write header info to the HydroDyn Output File
      !CALL OpenFOutFile ( 66, 'C:\Dev\NREL_SVN\HydroDyn\branches\HydroDyn_Modularization\Samples\NRELOffshrBsline5MW_OC3Hywind\WaveElevC0_HD.txt', ErrStat   )  ! Open motion file.
      !DO K = 0, InitInp%NStepWave2
      !   WRITE ( 66, '(2(e20.9))', IOSTAT = ErrStat) REAL(K), REAL(InitInp%WaveElevC0(K))
      !END DO
      !CLOSE ( 66 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Dump the WaveExctnC variable to a file for debugging
      ! Open and write header info to the HydroDyn Output File
      !CALL OpenFOutFile ( 66, 'C:\Dev\NREL_SVN\HydroDyn\branches\HydroDyn_Modularization\Samples\NRELOffshrBsline5MW_OC3Hywind\WaveExctnC_HD.txt', ErrStat   )  ! Open motion file.
      !DO K = 0, InitInp%NStepWave2 
      !   WRITE ( 66, '(7(e20.9))', IOSTAT = ErrStat) REAL(K), REAL(WaveExctnC(K,:))
      !END DO
      !CLOSE ( 66 )
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              
            else if ( p%ExctnMod == 2 ) then
               Interval_Sub                  = InitInp%Conv_Rdtn%RdtnDT
               SS_Exctn_InitInp%InputFile    = InitInp%WAMITFile    
               SS_Exctn_InitInp%WaveDir      = InitInp%WaveDir
               SS_Exctn_InitInp%NStepWave    = p%NStepWave
               SS_Exctn_InitInp%NBody        = InitInp%NBody
               SS_Exctn_InitInp%PtfmRefztRot = InitInp%PtfmRefztRot
               SS_Exctn_InitInp%SeaSt_Interp_p = InitInp%SeaSt_Interp_p
               SS_Exctn_InitInp%ExctnDisp    = InitInp%ExctnDisp
               ! We have been passed a pointer to WaveElev0 for use by the State Space excitation module.
               ! If the special case shown below is not used, then the state space model simply uses WaveElev0, as is.
               ! however, if we are using the special case, then WaveElev0 will be modified.  This is okay, because no one else
               ! is using WaveElev0 data
               if (p%ExctnDisp == 0 ) then
                  if (allocated(InitInp%WaveElev0)) then
                  
                     ! No other modules need this WaveElev0 array so we will simply move the allocation over to the SS_Exctn module
                     call MOVE_ALLOC(InitInp%WaveElev0, SS_Exctn_InitInp%WaveElev0) 
                    ! SS_Exctn_InitInp%WaveElev0 => InitInp%WaveElev0
               
                     ! Handle special case when NBodyMod=2 and (PtfmRefxt /= 0 or PtfmRefyt /= 0)  : Need to phase shift the wave elevation data for the offset body
                     if ( p%NBodyMod==2 .and. (InitInp%PtfmRefxt(1) /= 0 .or. InitInp%PtfmRefyt(1) /= 0) ) then
                  
                        ! Need to start with the DFT of the Wave Elevation data at the Platform reference point: InitInp%WaveElevC0
               
                        ! Now apply the phase shift in the frequency space

                        do J = 1, NInpWvDir  
                           do I = 0,InitInp%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transform

                        ! Compute the frequency of this component:

                              Omega = I*InitInp%WaveDOmega
                                 ! Fxy = exp(-j * k(w) * ( X*cos(Beta(w)) + Y*sin(Beta(w)) )
                              WaveNmbr   = WaveNumber ( Omega, InitInp%Gravity, InitInp%WtrDpth )
                              tmpAngle   = WaveNmbr * ( InitInp%PtfmRefxt(1)*cos(HdroWvDir(J)*D2R) + InitInp%PtfmRefyt(1)*sin(HdroWvDir(J)*D2R) )
                              TmpRe =  cos(tmpAngle)
                              TmpIm = -sin(tmpAngle)
                              Fxy   = CMPLX( TmpRe, TmpIm )

                              tmpComplexArr(I) = Fxy*CMPLX(InitInp%WaveElevC0(1,I), InitInp%WaveElevC0(2,I))
                          

                           end do
                        end do  
                     
                        ! Compute the inverse discrete Fourier transforms to find the time-domain
                        !   representations of the wave kinematics without stretching:

                        CALL InitFFT ( InitInp%NStepWave, FFT_Data, .TRUE., ErrStat2 )
                        CALL SetErrStat(ErrStat2,'Error occured while initializing the FFT.',ErrStat,ErrMsg,RoutineName)
                        IF ( ErrStat >= AbortErrLev ) THEN
                           CALL CleanUp()
                           RETURN
                        END IF
      
                           ! We'll need the following for wave stretching once we implement it.
                        CALL ApplyFFT_cx (  SS_Exctn_InitInp%WaveElev0(0:InitInp%NStepWave-1),  tmpComplexArr(:  ), FFT_Data, ErrStat2 )
                        CALL SetErrStat(ErrStat2,'Error occured while applying the FFT to WaveElev0.',ErrStat,ErrMsg,RoutineName)
                        IF ( ErrStat >= AbortErrLev ) THEN
                           CALL CleanUp()
                           RETURN
                        END IF
                  
                        CALL ExitFFT(FFT_Data, ErrStat2)
                           CALL SetErrStat( ErrStat2, 'Error in call to ExitFFT.', ErrStat, ErrMsg, RoutineName)
                           IF ( ErrStat >= AbortErrLev) THEN
                              CALL Cleanup()
                              RETURN
                           END IF
                     
                     end if 
                  else
                     !TODO: Error message because we need WaveElev0 for ExctnDisp=0
                  end if
               else
                     SS_Exctn_InitInp%WaveElev1 => InitInp%WaveElev1
               end if
                  
                     ! We need the WaveTime array to stay intact for use in other modules, so we will make a copy instead of moving the allocation
                  !ALLOCATE ( SS_Exctn_InitInp%WaveTime (0:InitInp%NStepWave) , STAT=ErrStat2 )
                  !IF ( ErrStat2 /= 0 )  THEN
                  !   CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the SS_Exctn_InitInp%WaveTime array.', ErrStat, ErrMsg, RoutineName)
                  !   CALL Cleanup()
                  !   RETURN            
                  !END IF
               SS_Exctn_InitInp%WaveTime => InitInp%WaveTime 
                  
      
               
               call SS_Exc_Init(SS_Exctn_InitInp, m%SS_Exctn_u, p%SS_Exctn, x%SS_Exctn, xd%SS_Exctn, z%SS_Exctn, OtherState%SS_Exctn, &
                                      m%SS_Exctn_y, m%SS_Exctn, Interval_Sub, SS_Exctn_InitOut, ErrStat2, ErrMsg2)
            
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  if ( ErrStat >= AbortErrLev ) then
                     call Cleanup()
                     return
                  end if   
            end if
            
            IF ( (p%ExctnMod>0) .AND. (p%ExctnDisp==2) ) THEN ! Allocate array for filtered potential-flow body positions
               p%ExctnFiltConst = exp(-2.0*Pi*p%ExctnCutOff * Interval)
               ALLOCATE ( xd%BdyPosFilt(1:2, 1:p%NBody, 1:3) , STAT=ErrStat2 )
               IF ( ErrStat2 /= 0 )  THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the BdyPosFilt array.', ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  RETURN            
               END IF
               xd%BdyPosFilt = 0.0_ReKi
            END IF

         CASE ( 6 )              ! User wave data.

            CALL SetErrStat( ErrID_Fatal, 'User input wave data not applicable for floating platforms.', ErrStat, ErrMsg, RoutineName)
            CALL Cleanup()
            RETURN

         ENDSELECT   
      end if
      
      IF ( InitInp%RdtnTMax == 0.0 )  THEN   ! .TRUE. when we don't want to model wave radiation damping; set RdtnTMax to some minimum value greater than zero to avoid an error in the calculations below.
      
         p%RdtnMod   = 0
         
      ELSE                                    ! We will be modeling wave radiation damping.
                 
         p%RdtnMod   = InitInp%RdtnMod
         
         if ( InitInp%RdtnMod == 1 ) THEN
            
            ! this check used to occur with equivalent variables after calling Conv_Rdtn_Init
            if (.not. EqualRealNos( InitInp%Conv_Rdtn%RdtnDT, Interval) ) then
               call SetErrStat(ErrID_Fatal,'RdtnDT must be the same as the HD time step', ErrStat, ErrMsg, RoutineName)
               call Cleanup()
               return
            end if
               
            ! Set Initialization data for the Conv_Rdtn submodule
            ! Would be nice if there were a copy InitInput function in the *_Types file
            ! BJJ 6/25/2014: There is a copy InitInput function.... ???
            
            CALL MOVE_ALLOC( HdroFreq,  Conv_Rdtn_InitInp%HdroFreq  )
            CALL MOVE_ALLOC( HdroAddMs, Conv_Rdtn_InitInp%HdroAddMs )
            CALL MOVE_ALLOC( HdroDmpng, Conv_Rdtn_InitInp%HdroDmpng )
            Conv_Rdtn_InitInp%NBody               = InitInp%NBody
            Conv_Rdtn_InitInp%RdtnTMax            = InitInp%RdtnTMax
            Conv_Rdtn_InitInp%RdtnDT              = InitInp%Conv_Rdtn%RdtnDT
            Conv_Rdtn_InitInp%HighFreq            = HighFreq
            Conv_Rdtn_InitInp%WAMITFile           = InitInp%WAMITFile
            Conv_Rdtn_InitInp%NInpFreq            = NInpFreq
    
         
            CALL Conv_Rdtn_Init(Conv_Rdtn_InitInp, m%Conv_Rdtn_u, p%Conv_Rdtn, x%Conv_Rdtn, xd%Conv_Rdtn, z%Conv_Rdtn, OtherState%Conv_Rdtn, &
                                   m%Conv_Rdtn_y, m%Conv_Rdtn, Conv_Rdtn_InitOut, ErrStat2, ErrMsg2)
            
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL Cleanup()
                  RETURN
               END IF
               

            
         ELSE IF ( InitInp%RdtnMod == 2 ) THEN
            
            SS_Rdtn_InitInp%InputFile    = InitInp%WAMITFile    

            call AllocAry(SS_Rdtn_InitInp%enabledDOFs, 6*p%NBody, 'SS_Rdtn_InitInp%enabledDOFs', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL Cleanup()
                  RETURN
               END IF
            SS_Rdtn_InitInp%enabledDOFs  = 1                        !  Set to 1 (True) for all DOFs, meaning each DOF is to be used in the analysis.   
            Interval_Sub                 = InitInp%Conv_Rdtn%RdtnDT
            SS_Rdtn_InitInp%NBody        = InitInp%NBody
            call AllocAry(SS_Rdtn_InitInp%PtfmRefztRot, p%NBody, 'SS_Rdtn_InitInp%PtfmRefztRot', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL Cleanup()
                  RETURN
               END IF
            SS_Rdtn_InitInp%PtfmRefztRot = InitInp%PtfmRefztRot
            CALL SS_Rad_Init(SS_Rdtn_InitInp, m%SS_Rdtn_u, p%SS_Rdtn, x%SS_Rdtn, xd%SS_Rdtn, z%SS_Rdtn, OtherState%SS_Rdtn, &
                                   m%SS_Rdtn_y, m%SS_Rdtn, Interval_Sub, SS_Rdtn_InitOut, ErrStat2, ErrMsg2)
            
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
               IF ( ErrStat >= AbortErrLev ) THEN
                  CALL Cleanup()
                  RETURN
               END IF            
            
         END IF
         
      END IF
      
         ! create the input and output meshes
         ! CALL MeshCreate(u%MeshData, COMPONENT_INPUT, 1, ErrStat2, ErrMsg2, .TRUE.)
         ! deallocate arrays

      IF ( ALLOCATED( HdroExctn    ) ) DEALLOCATE( HdroExctn    )
      IF ( ALLOCATED( WaveExctnC   ) ) DEALLOCATE( WaveExctnC   )
      IF ( ALLOCATED( HdroAddMs    ) ) DEALLOCATE( HdroAddMs    )
      IF ( ALLOCATED( HdroDmpng    ) ) DEALLOCATE( HdroDmpng    )
      IF ( ALLOCATED( HdroFreq     ) ) DEALLOCATE( HdroFreq     )
      IF ( ALLOCATED( HdroWvDir    ) ) DEALLOCATE( HdroWvDir    )
      
      IF ( ALLOCATED( WAMITFreq    ) ) DEALLOCATE( WAMITFreq    )
      IF ( ALLOCATED( WAMITPer     ) ) DEALLOCATE( WAMITPer     )
      IF ( ALLOCATED( WAMITWvDir   ) ) DEALLOCATE( WAMITWvDir   )
      IF ( ALLOCATED( SortFreqInd  ) ) DEALLOCATE( SortFreqInd  )
      IF ( ALLOCATED( SortWvDirInd ) ) DEALLOCATE( SortWvDirInd )

      
      
      ! Define system output initializations (set up meshes) here:
      
     
          ! Create the input and output meshes associated with lumped loads
     
      
     
         call MeshCreate( BlankMesh        = u%Mesh            &
                        ,IOS               = COMPONENT_INPUT   &
                        ,Nnodes            = p%NBody           &
                        ,ErrStat           = ErrStat2          &
                        ,ErrMess           = ErrMsg2           &
                        ,TranslationDisp   = .TRUE.            &
                        ,Orientation       = .TRUE.            &
                        ,TranslationVel    = .TRUE.            &
                        ,RotationVel       = .TRUE.            &
                        ,TranslationAcc    = .TRUE.            &
                        ,RotationAcc       = .TRUE.)
         
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL Cleanup()
               RETURN
            END IF
         
      do iBody = 1, p%NBody

         theta = (/ 0.0_R8Ki, 0.0_R8Ki, 0.0_R8Ki /)
         orientation = EulerConstruct(theta)
         
         
            ! Create the node on the mesh
  
         CALL MeshPositionNode (u%Mesh                                &
                                 , iBody                              &
                                 , (/InitInp%PtfmRefxt(iBody), InitInp%PtfmRefyt(iBody), InitInp%PtfmRefzt(iBody)/)   &  
                                 , ErrStat2                           &
                                 , ErrMsg2                            &
                                 , orientation )
      
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      
            ! Create the mesh element
         CALL MeshConstructElement (  u%Mesh              &
                                     , ELEMENT_POINT      &                         
                                     , ErrStat2           &
                                     , ErrMsg2            &
                                     , iBody              &
                                                 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
      end do

      CALL MeshCommit ( u%Mesh              &
                        , ErrStat2            &
                        , ErrMsg2             )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF      

        call MeshCopy ( SrcMesh   = u%Mesh           &
                       ,DestMesh  = y%Mesh           &
                       ,CtrlCode  = MESH_SIBLING     &
                       ,IOS       = COMPONENT_OUTPUT &
                       ,ErrStat   = ErrStat2         &
                       ,ErrMess   = ErrMsg2          &
                       ,Force     = .TRUE.           &
                       ,Moment    = .TRUE.           )
        
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF    
      u%Mesh%RemapFlag  = .TRUE.
      y%Mesh%RemapFlag  = .TRUE.

      
         ! Define initialization-routine output here:
         
      
     
                                               

      ! initialize misc vars:      
   m%LastIndWave = 1
       
       CALL Cleanup()
       
CONTAINS


   SUBROUTINE Cleanup()
   
      ! destroy local variables that are types in the framework:
      
      CALL Conv_Rdtn_DestroyInitInput(  Conv_Rdtn_InitInp,  ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )
      CALL Conv_Rdtn_DestroyInitOutput( Conv_Rdtn_InitOut,  ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )

      CALL SS_Rad_DestroyInitInput(     SS_Rdtn_InitInp,    ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )
      CALL SS_Rad_DestroyInitOutput(    SS_Rdtn_InitOut,    ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )
      CALL SS_Exc_DestroyInitInput(     SS_Exctn_InitInp,   ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )
      CALL SS_Exc_DestroyInitOutput(    SS_Exctn_InitOut,   ErrStat2, ErrMsg2, DEALLOCATEpointers=.false. )
      
      
      ! destroy local variables that are allocatable arrays:
      
      IF ( ALLOCATED( HdroExctn   ) ) DEALLOCATE(HdroExctn   )
      IF ( ALLOCATED( WaveExctnC  ) ) DEALLOCATE(WaveExctnC  )
      
      IF ( ALLOCATED( HdroAddMs   ) ) DEALLOCATE(HdroAddMs   )
      IF ( ALLOCATED( HdroDmpng   ) ) DEALLOCATE(HdroDmpng   )
      IF ( ALLOCATED( HdroFreq    ) ) DEALLOCATE(HdroFreq    )
      IF ( ALLOCATED( HdroWvDir   ) ) DEALLOCATE(HdroWvDir   )
      
      IF ( ALLOCATED( WAMITFreq   ) ) DEALLOCATE(WAMITFreq   )
      IF ( ALLOCATED( WAMITPer    ) ) DEALLOCATE(WAMITPer    )
      IF ( ALLOCATED( WAMITWvDir  ) ) DEALLOCATE(WAMITWvDir  )
      IF ( ALLOCATED( SortFreqInd ) ) DEALLOCATE(SortFreqInd )
      IF ( ALLOCATED( SortWvDirInd) ) DEALLOCATE(SortWvDirInd)
   
   
   END SUBROUTINE Cleanup
            
END SUBROUTINE WAMIT_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE WAMIT_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(WAMIT_InputType),             INTENT(INOUT)  :: u           !< System inputs
      TYPE(WAMIT_ParameterType),         INTENT(INOUT)  :: p           !< Parameters     
      TYPE(WAMIT_ContinuousStateType),   INTENT(INOUT)  :: x           !< Continuous states
      TYPE(WAMIT_DiscreteStateType),     INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(WAMIT_ConstraintStateType),   INTENT(INOUT)  :: z           !< Constraint states
      TYPE(WAMIT_OtherStateType),        INTENT(INOUT)  :: OtherState  !< Other states            
      TYPE(WAMIT_OutputType),            INTENT(INOUT)  :: y           !< System outputs
      TYPE(WAMIT_MiscVarType),           INTENT(INOUT)  :: m           !< Initial misc/optimization variables            
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:


         ! Close files here:     
                  
                  

         ! Destroy the input data:
         
      CALL WAMIT_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:
      CALL WAMIT_DestroyParam( p, ErrStat, ErrMsg, DEALLOCATEpointers=.false. )


         ! Destroy the state data:
         
      CALL WAMIT_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL WAMIT_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL WAMIT_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL WAMIT_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         
         ! Destroy misc vars:
      CALL WAMIT_DestroyMisc(  m,  ErrStat, ErrMsg )
      
         ! Destroy the output data:
         
      CALL WAMIT_DestroyOutput( y, ErrStat, ErrMsg )


      

END SUBROUTINE WAMIT_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE WAMIT_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
      TYPE(WAMIT_InputType),              INTENT(IN   ) :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(WAMIT_ParameterType),          INTENT(IN   ) :: p               !< Parameters
      TYPE(WAMIT_ContinuousStateType),    INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                           !!   Output: Continuous states at t + Interval
      TYPE(WAMIT_DiscreteStateType),      INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                           !!   Output: Discrete states at t + Interval
      TYPE(WAMIT_ConstraintStateType),    INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                           !!   Output: Constraint states at t + Interval
      TYPE(WAMIT_OtherStateType),         INTENT(INOUT) :: OtherState      !< Input: Other states at t;
                                                                           !!   Output: Other states at t + Interval
      TYPE(WAMIT_MiscVarType),            INTENT(INOUT) :: m               !< Misc/optimization variables            
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None

         ! Local variables

      INTEGER                                           :: I               ! Generic loop counter
      INTEGER                                           :: nTime           ! Number of inputs
      integer(IntKi)                                    :: iBody           ! WAMIT body index

      integer(IntKi)                                    :: indxStart, indxEnd                   ! Starting and ending indices for the iBody_th sub vector in an NBody long vector

!      INTEGER(IntKi)                                    :: ErrStat2        ! Error status of the operation (secondary error)
!      CHARACTER(ErrMsgLen)                              :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
      REAL(ReKi)                                        :: bodyPosition(2)
      
          ! Create dummy variables required by framework but which are not used by the module
      
      TYPE(Conv_Rdtn_InputType), ALLOCATABLE :: Conv_Rdtn_u(:)         ! Inputs
      
      TYPE(SS_Rad_InputType), ALLOCATABLE    :: SS_Rdtn_u(:)           ! Inputs
      TYPE(SS_Exc_InputType), ALLOCATABLE    :: SS_Exctn_u(:)          ! Inputs
      TYPE(WAMIT_InputType),  ALLOCATABLE    :: WAMIT_u(:)             ! Inputs
      TYPE(WAMIT_InputType)                  :: WAMIT_u_t
      
                        
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      nTime = size(Inputs)   
      
      
      IF ( p%RdtnMod == 1 ) THEN       ! Update the convolution radiation memory effect sub-module's state  
         
            ! Allocate array of Conv_Rdtn inputs
        
         ALLOCATE( Conv_Rdtn_u(nTime), STAT = ErrStat )
         IF (ErrStat /=0) THEN
            ErrMsg = ' Failed to allocate array Conv_Rdtn_u.'
            RETURN
         END IF
         
         DO I=1,nTime
            ALLOCATE( Conv_Rdtn_u(I)%Velocity(size(m%Conv_Rdtn_u%Velocity)), STAT = ErrStat  )
            IF (ErrStat /=0) THEN
               ErrMsg = ' Failed to allocate array Conv_Rdtn_u(I)%Velocity.'
               RETURN
            END IF
            do iBody=1,p%NBody
               indxStart = (iBody-1)*6+1
               indxEnd   = indxStart+5        
               Conv_Rdtn_u(I)%Velocity(indxStart:indxEnd) = (/Inputs(I)%Mesh%TranslationVel(:,iBody), Inputs(I)%Mesh%RotationVel(:,iBody)/) 
            end do
         END DO
                 
         CALL Conv_Rdtn_UpdateStates( t, n, Conv_Rdtn_u, InputTimes, p%Conv_Rdtn, x%Conv_Rdtn, xd%Conv_Rdtn, &
                                      z%Conv_Rdtn, OtherState%Conv_Rdtn, m%Conv_Rdtn, ErrStat, ErrMsg )
         
         DEALLOCATE(Conv_Rdtn_u)
         
      ELSE IF ( p%RdtnMod == 2 )  THEN       ! Update the state-space radiation memory effect sub-module's state      
          
           ! Allocate array of SS_Rdtn inputs
        
         ALLOCATE( SS_Rdtn_u(nTime), STAT = ErrStat )
         IF (ErrStat /=0) THEN
            ErrMsg = ' Failed to allocate array SS_Rdtn_u.'
            RETURN
         END IF
         
         DO I=1,nTime
            ALLOCATE( SS_Rdtn_u(I)%dq(size(m%SS_Rdtn_u%dq)), STAT = ErrStat  )
            IF (ErrStat /=0) THEN
               ErrMsg = ' Failed to allocate array SS_Rdtn_u(I)%dq.'
               RETURN
            END IF
            do iBody=1,p%NBody
               indxStart = (iBody-1)*6+1 
               SS_Rdtn_u(I)%dq(indxStart:indxStart+2)   = Inputs(I)%Mesh%TranslationVel(:,iBody)
               SS_Rdtn_u(I)%dq(indxStart+3:indxStart+5) = Inputs(I)%Mesh%RotationVel(:,iBody)
               !SS_Rdtn_u(I)%dq = reshape((/Inputs(I)%Mesh%TranslationVel(:,1), Inputs(I)%Mesh%RotationVel(:,1)/), (/6,1/)) !reshape(u%Velocity, (/6,1/)) ! dq is a 6x1 matrix
            end do
         END DO
         
         CALL SS_Rad_UpdateStates( t, n, SS_Rdtn_u, InputTimes, p%SS_Rdtn, x%SS_Rdtn, xd%SS_Rdtn, z%SS_Rdtn, OtherState%SS_Rdtn, m%SS_Rdtn, ErrStat, ErrMsg )
         
         DEALLOCATE(SS_Rdtn_u)
         
      END IF
      
      IF ( (p%ExctnMod>0).AND.(p%ExctnDisp==2) ) THEN
         ALLOCATE( WAMIT_u(nTime), STAT = ErrStat )
         IF (ErrStat /=0) THEN
            ErrMsg = ' Failed to allocate array WAMIT_u.'
            RETURN
         END IF
         DO I=1,nTime
            ALLOCATE( WAMIT_u(I)%Mesh%TranslationDisp(3,p%NBody), STAT = ErrStat  )
            IF (ErrStat /=0) THEN
               ErrMsg = ' Failed to allocate array WAMIT_u(I)%Mesh%TranslationDisp.'
               RETURN
            END IF
            DO iBody=1,p%NBody
               WAMIT_u(I)%Mesh%TranslationDisp(:,iBody) = Inputs(I)%Mesh%TranslationDisp(:,iBody)
            END DO
         END DO
         ! Interpolate WAMIT input at time t+dt
         CALL WAMIT_Input_ExtrapInterp(WAMIT_u, InputTimes, WAMIT_u_t, t+p%dt, ErrStat, ErrMsg)
         DO iBody = 1,p%NBody
            ! Current unfiltered body position at time t+dt
            bodyPosition(1) = WAMIT_u_t%Mesh%TranslationDisp(1,iBody)
            bodyPosition(2) = WAMIT_u_t%Mesh%TranslationDisp(2,iBody)
            ! Filtered body position
            xd%BdyPosFilt(:,iBody,3) = xd%BdyPosFilt(:,iBody,2)
            xd%BdyPosFilt(:,iBody,2) = xd%BdyPosFilt(:,iBody,1)
            xd%BdyPosFilt(1,iBody,1) = p%ExctnFiltConst * xd%BdyPosFilt(1,iBody,1) + (1.0_ReKi - p%ExctnFiltConst) * bodyPosition(1)
            xd%BdyPosFilt(2,iBody,1) = p%ExctnFiltConst * xd%BdyPosFilt(2,iBody,1) + (1.0_ReKi - p%ExctnFiltConst) * bodyPosition(2)  
         END DO
         CALL WAMIT_DestroyInput( WAMIT_u_t, ErrStat, ErrMsg)
         DO I=1,nTime
            CALL WAMIT_DestroyInput( WAMIT_u(I), ErrStat, ErrMsg)
         END DO
         DEALLOCATE(WAMIT_u)
      END IF
      IF ( p%ExctnMod == 2 )  THEN       ! Update the state-space wave excitation sub-module's states      
          
           ! Allocate array of dummy SS_Excitation inputs for the framework
        
         allocate( SS_Exctn_u(nTime), STAT = ErrStat )
         if (ErrStat /=0) then
            ErrMsg = ' Failed to allocate array SS_Exctn_u.'
            return
         end if
         if (p%ExctnDisp == 1) then ! Use unfiltered position
            DO I=1,nTime
               ALLOCATE( SS_Exctn_u(I)%PtfmPos(3,p%NBody), STAT = ErrStat  )
               IF (ErrStat /=0) THEN
                  ErrMsg = ' Failed to allocate array SS_Exctn_u(I)%PtfmPos.'
                  RETURN
               END IF
               if (p%NBodyMod == 2) then
                  do iBody=1,p%NBody                  
                     SS_Exctn_u(I)%PtfmPos(:,iBody)   = Inputs(I)%Mesh%TranslationDisp(:,iBody) + Inputs(I)%Mesh%Position(:,iBody)
                  end do
               else
                  do iBody=1,p%NBody                  
                     SS_Exctn_u(I)%PtfmPos(:,iBody)   = Inputs(I)%Mesh%TranslationDisp(:,iBody)
                  end do
               end if
               
            END DO
         else if (p%ExctnDisp == 2) then ! Use filtered position (only need x and y coordinates)       
            DO I=1,nTime
               ALLOCATE( SS_Exctn_u(I)%PtfmPos(3,p%NBody), STAT = ErrStat  )
               IF (ErrStat /=0) THEN
                  ErrMsg = ' Failed to allocate array SS_Exctn_u(I)%PtfmPos.'
                  RETURN
               END IF
               if (p%NBodyMod == 2) then
                  do iBody=1,p%NBody                  
                     SS_Exctn_u(I)%PtfmPos(1:2,iBody)   = xd%BdyPosFilt(:,iBody,I) + Inputs(I)%Mesh%Position(1:2,iBody)
                  end do
               else
                  do iBody=1,p%NBody                  
                     SS_Exctn_u(I)%PtfmPos(1:2,iBody)   = xd%BdyPosFilt(:,iBody,I)
                  end do
               end if
               
            END DO

         end if
         
         call SS_Exc_UpdateStates( t, n, SS_Exctn_u, InputTimes, p%SS_Exctn, x%SS_Exctn, xd%SS_Exctn, z%SS_Exctn, OtherState%SS_Exctn, m%SS_Exctn, ErrStat, ErrMsg )
         
         deallocate(SS_Exctn_u)
         
      END IF
      
END SUBROUTINE WAMIT_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE WAMIT_CalcOutput( Time, WaveTime, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                      INTENT(IN   )  :: Time        !< Current simulation time in seconds
      real(SiKi),                      intent(in   )  :: WaveTime(:) !< Array of wave kinematic time samples, (sec)
      TYPE(WAMIT_InputType),           INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(WAMIT_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(WAMIT_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(WAMIT_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(WAMIT_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(WAMIT_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(WAMIT_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                     !!   nectivity information does not have to be recalculated)
      TYPE(WAMIT_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables            
      INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      
            
         ! Local Variables:
      !REAL(ReKi)                           :: F_HS     (6)                            ! Total load contribution from hydrostatics, including the effects of waterplane area and the center of buoyancy (N, N-m)
      !REAL(ReKi)                           :: F_Waves  (6)                            ! Total load contribution from incident waves (i.e., the diffraction problem) (N, N-m)   
      !REAL(ReKi)                           :: F_Rdtn   (6)                            ! Total load contribution from wave radiation damping (i.e., the diffraction problem) (N, N-m)
      INTEGER(IntKi)                       :: I,iStart                                ! Generic index
      REAL(ReKi)                           :: q(6*p%NBody), qdot(6*p%NBody), qdotdot(6*p%NBody)  ! kinematics for all WAMIT bodies
      REAL(ReKi)                           :: rotdisp(3)                              ! small angle rotational displacements
      integer(IntKi)                       :: iBody                                   ! Counter for WAMIT bodies.  If NBodyMod > 1 then NBody = 1, and hence iBody = 1
      integer(IntKi)                       :: indxStart, indxEnd                      ! Starting and ending indices for the iBody_th sub vector in an NBody long vector
      real(ReKi)                           :: bodyPosition(2)                         ! x-y displaced location of a WAMIT body (relative to 
               ! Error handling
      CHARACTER(1024)                        :: ErrMsg2                              ! Temporary error message for calls
      INTEGER(IntKi)                         :: ErrStat2                             ! Temporary error status for calls

         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute outputs here:
         
      
      
         ! Compute the load contribution from incident waves (i.e., the diffraction problem):
      if ( p%ExctnMod == 0 ) then
         
         m%F_Waves1 = 0.0_ReKi
         
      else if ( p%ExctnMod == 1 ) then

         if ( p%ExctnDisp == 0 ) then         
               ! Abort if the wave excitation loads have not been computed yet:
            IF ( .NOT. ALLOCATED ( p%WaveExctn ) )  THEN
               ErrMsg  = ' Routine WAMIT_Init() must be called before routine WAMIT_CalcOutput().'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
         
            DO I = 1,6*p%NBody     ! Loop through all wave excitation forces and moments
               m%F_Waves1(I) = InterpWrappedStpReal ( REAL(Time, SiKi), WaveTime(:), p%WaveExctn(:,I), &
                                                        m%LastIndWave, p%NStepWave + 1       )
            END DO          ! I - All wave excitation forces and moments
         else ! p%ExctnDisp > 0
            IF ( .NOT. allocated ( p%WaveExctnGrid ) )  THEN
               ErrMsg  = ' Routine WAMIT_Init() must be called before routine WAMIT_CalcOutput().'
               ErrStat = ErrID_Fatal
               RETURN
            END IF
            ! We are using the displaced x,y location of the WAMIT bodies to determine the Wave Exication force

            DO iBody  = 1,p%NBody
               IF ( p%ExctnDisp == 1 ) THEN
                  ! Current unfiltered body position
                  bodyPosition(1) = u%Mesh%TranslationDisp(1,iBody)
                  bodyPosition(2) = u%Mesh%TranslationDisp(2,iBody)
               ELSE IF ( p%ExctnDisp == 2 ) THEN
                  ! Use filtered body position
                  bodyPosition(1) = xd%BdyPosFilt(1,iBody,1)
                  bodyPosition(2) = xd%BdyPosFilt(2,iBody,1)
               END IF
               iStart = (iBody-1)*6+1
               ! WaveExctnGrid dimensions are: 1st: wavetime, 2nd: X, 3rd: Y, 4th: Force component for each WAMIT Body
               m%F_Waves1(iStart:iStart+5) = SeaSt_Interp_3D_Vec6( Time, bodyPosition, p%WaveExctnGrid(:,:,:,iStart:iStart+5), p%SeaSt_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                  call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SeaState_CalcOutput' )
            END DO
         end if
         
      else if ( p%ExctnMod == 2 ) then
         
         call SS_Exc_CalcOutput( Time, m%SS_Exctn_u, p%SS_Exctn, x%SS_Exctn, xd%SS_Exctn,  &
                                z%SS_Exctn, OtherState%SS_Exctn, m%SS_Exctn_y, m%SS_Exctn, ErrStat, ErrMsg )
         m%F_Waves1  (:) = m%SS_Exctn_y%y
         
      end if
      
      
      do iBody = 1, p%NBody
         
            ! Determine the rotational angles from the direction-cosine matrix
         rotdisp   = GetSmllRotAngs ( u%Mesh%Orientation(:,:,iBody), ErrStat, ErrMsg )
         indxStart = (iBody-1)*6+1
         indxEnd   = indxStart+5
         q      (indxStart:indxEnd)   = reshape((/real(u%Mesh%TranslationDisp(:,iBody),ReKi),rotdisp(:)/),(/6/))
         qdot   (indxStart:indxEnd)   = reshape((/u%Mesh%TranslationVel(:,iBody),u%Mesh%RotationVel(:,iBody)/),(/6/))
         qdotdot(indxStart:indxEnd)   = reshape((/u%Mesh%TranslationAcc(:,iBody),u%Mesh%RotationAcc(:,iBody)/),(/6/))
         
      end do
      
         ! Compute the load contribution from hydrostatics:

      m%F_HS = -matmul(p%HdroSttc,q)
      
      do iBody = 1, p%NBody
         indxStart = (iBody-1)*6+1
         indxEnd   = indxStart+5
         m%F_HS(indxStart:indxEnd) =  m%F_HS(indxStart:indxEnd) + p%F_HS_Moment_Offset(:,iBody)  ! except for the hydrostatic buoyancy force from Archimede's Principle when the support platform is in its undisplaced position
      end do   
      
     
         ! If necessary, compute the load contribution from wave radiation damping
         !   (i.e., the radiation problem):

      IF ( p%RdtnMod == 1 )  THEN ! .TRUE. when we will be modeling wave radiation damping.
         m%Conv_Rdtn_u%Velocity = qdot
         CALL Conv_Rdtn_CalcOutput( Time, m%Conv_Rdtn_u, p%Conv_Rdtn, x%Conv_Rdtn, xd%Conv_Rdtn,  &
                                z%Conv_Rdtn, OtherState%Conv_Rdtn, m%Conv_Rdtn_y, m%Conv_Rdtn, ErrStat, ErrMsg )
         m%F_Rdtn  (:) = m%Conv_Rdtn_y%F_Rdtn       

      ELSE IF ( p%RdtnMod == 2 )  THEN 
         m%SS_Rdtn_u%dq = qdot
         CALL SS_Rad_CalcOutput( Time, m%SS_Rdtn_u, p%SS_Rdtn, x%SS_Rdtn, xd%SS_Rdtn,  &
                                z%SS_Rdtn, OtherState%SS_Rdtn, m%SS_Rdtn_y, m%SS_Rdtn, ErrStat, ErrMsg )
         m%F_Rdtn  (:) = m%SS_Rdtn_y%y
      ELSE ! We must not be modeling wave radiation damping.


      ! Set the total load contribution from radiation damping to zero:

         m%F_Rdtn        (:) = 0.0


      END IF       
      
      
      ! Compute Added Mass Forces
      
         ! Set the platform added mass matrix, PtfmAM, to be the infinite-frequency
         !   limit of the frequency-dependent hydrodynamic added mass matrix,
         !   HdroAdMsI:
         
         !added mass:

      m%F_PtfmAM     =   -matmul(p%HdroAdMsI, qdotdot)

      
      
         ! Compute outputs here:
      do iBody = 1, p%NBody
         indxStart = (iBody-1)*6
         
         DO I=1,3
            y%Mesh%Force(I,iBody)    = m%F_PtfmAM(indxStart+I)   + m%F_Rdtn(indxStart+I)   + m%F_Waves1(indxStart+I)   + m%F_HS(indxStart+I) 
         END DO
         DO I=1,3
            y%Mesh%Moment(I,iBody)   = m%F_PtfmAM(indxStart+I+3) + m%F_Rdtn(indxStart+I+3) + m%F_Waves1(indxStart+I+3) + m%F_HS(indxStart+I+3)
         END DO
      end do
      
      
      ! Output channels will be dealt with by the HydroDyn module
             

END SUBROUTINE WAMIT_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE WAMIT_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )  
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(WAMIT_InputType),             INTENT(IN   )  :: u           !< Inputs at Time                    
      TYPE(WAMIT_ParameterType),         INTENT(IN   )  :: p           !< Parameters                             
      TYPE(WAMIT_ContinuousStateType),   INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(WAMIT_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(WAMIT_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(WAMIT_OtherStateType),        INTENT(IN   )  :: OtherState  !< Other states                    
      TYPE(WAMIT_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables            
      TYPE(WAMIT_ContinuousStateType),   INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      integer(IntKi)                                    :: iBody                ! WAMIT body index
      integer(IntKi)                                    :: indxStart            ! Starting and ending indices for the iBody_th sub vector in an NBody long vector
      real(SiKi)   :: waveElev0(p%NBody)   
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute the first time derivatives of the continuous states here:
      if (p%RdtnMod == 2) then
         do iBody = 1, p%NBody
            indxStart = (iBody-1)*6+1 
            m%SS_Rdtn_u%dq(indxStart:indxStart+2)   = u%Mesh%TranslationVel(:,iBody) 
            m%SS_Rdtn_u%dq(indxStart+3:indxStart+5) = u%Mesh%RotationVel(:,iBody) 
         end do
      
         CALL SS_Rad_CalcContStateDeriv( Time, m%SS_Rdtn_u, p%SS_Rdtn, x%SS_Rdtn, xd%SS_Rdtn, z%SS_Rdtn, OtherState%SS_Rdtn, m%SS_Rdtn, dxdt%SS_Rdtn, ErrStat, ErrMsg )      
      end if
         
         ! NOTE: The input below (0.0) will only work as part of a linearization Get_OP call! If this routine (WAMIT_CalcContStateDeriv) is called in another context, then the following
         ! input needs to be implemented generically. As of Aug 10, 2020, this is only called for Get_OP related work. GJH
      if (p%ExctnMod == 2) then
         waveElev0 = 0.0_SiKi
         CALL SS_Exc_CalcContStateDeriv( Time, waveElev0, p%SS_Exctn, x%SS_Exctn, xd%SS_Exctn, z%SS_Exctn, OtherState%SS_Exctn, m%SS_Exctn, dxdt%SS_Exctn, ErrStat, ErrMsg )      
      end if

END SUBROUTINE WAMIT_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
SUBROUTINE WAMIT_UpdateDiscState( Time, n, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                         INTENT(IN   )  :: Time        !< Current simulation time in seconds 
      INTEGER(IntKi),                     INTENT(IN   )  :: n           !< Current step of the simulation: t = n*Interval
      TYPE(WAMIT_InputType),              INTENT(IN   )  :: u           !< Inputs at Time                       
      TYPE(WAMIT_ParameterType),          INTENT(IN   )  :: p           !< Parameters                                 
      TYPE(WAMIT_ContinuousStateType),    INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(WAMIT_DiscreteStateType),      INTENT(INOUT)  :: xd          !< Input: Discrete states at Time; 
                                                                        !<   Output: Discrete states at Time + Interval
      TYPE(WAMIT_ConstraintStateType),    INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(WAMIT_OtherStateType),         INTENT(INOUT)  :: OtherState  !< Other states at Time (THIS [intent out] VIOLATES THE FRAMEWORK)          
      TYPE(WAMIT_MiscVarType),            INTENT(INOUT)  :: m           !< Misc/optimization variables            
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None
      
      integer(IntKi)                                     :: iBody               ! WAMIT body index
      integer(IntKi)                                     :: indxStart, indxEnd  ! Starting and ending indices for the iBody_th sub vector in an NBody long vector

         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Update discrete states here:
      IF ( p%RdtnMod == 1 )  THEN ! .TRUE. when we will be modeling wave radiation damping.   
         do iBody=1,p%NBody
               indxStart = (iBody-1)*6+1
               indxEnd   = indxStart+5    
               m%Conv_Rdtn_u%Velocity(indxStart:indxEnd) = (/u%Mesh%TranslationVel(:,iBody), u%Mesh%RotationVel(:,iBody)/)
         end do
         CALL Conv_Rdtn_UpdateDiscState( Time, n, m%Conv_Rdtn_u, p%Conv_Rdtn, x%Conv_Rdtn, xd%Conv_Rdtn, z%Conv_Rdtn, &
                                         OtherState%Conv_Rdtn, m%Conv_Rdtn, ErrStat, ErrMsg )
         
      END IF
         
       
END SUBROUTINE WAMIT_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE WAMIT_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds   
      TYPE(WAMIT_InputType),             INTENT(IN   )  :: u           !< Inputs at Time                       
      TYPE(WAMIT_ParameterType),         INTENT(IN   )  :: p           !< Parameters                           
      TYPE(WAMIT_ContinuousStateType),   INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(WAMIT_DiscreteStateType),     INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(WAMIT_ConstraintStateType),   INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
      TYPE(WAMIT_OtherStateType),        INTENT(IN   )  :: OtherState  !< Other states                    
      TYPE(WAMIT_MiscVarType),           INTENT(INOUT)  :: m           !< Misc/optimization variables            
      TYPE(WAMIT_ConstraintStateType),   INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using  
                                                                       !!     the input values described above      
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Solve for the constraint states here:
      call SS_Rad_CalcConstrStateResidual( Time, m%SS_Rdtn_u, p%SS_Rdtn, x%SS_Rdtn, xd%SS_Rdtn, z%SS_Rdtn, OtherState%SS_Rdtn, m%SS_Rdtn, z_residual%SS_Rdtn, ErrStat, ErrMsg ) 

END SUBROUTINE WAMIT_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE WAMIT
!**********************************************************************************************************************************
