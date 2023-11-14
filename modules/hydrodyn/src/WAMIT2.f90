!**********************************************************************************************************************************
!  NOTE: documentation in this file is written for use with Doxygen 1.8.6 and higher.
!
!> WAMIT2 module
!!
!!  This module calculates the second order wave forces on a structure. This module is used with HydroDyn in FAST.
!!
!!  This software is written in the FAST modular framework.
!!
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
!
!    This file is part of the WAMIT2 sub-module of HydroDyn.
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
MODULE WAMIT2



         !> Known issues:
         !!
         !! 1. _Problem:_
         !!       In the Read_DataFile3D and Read_DataFile4D routines, if the wave direction information crosses the +/-pi boundary,
         !!       it will end up sorted incorrectly.  This issue has not been fixed or addressed in any way yet.
         !!    _Effects:_
         !!       During checking of the wave directions to see if the data in the WAMIT file spans the full range of the required
         !!       wave directions when multidirectional waves are used, we will have an issue.  This will occur during the error
         !!       checking within the MnDrift_InitCalc, NewmanApp_InitCalc, DiffQTF_InitCalc, and SumQTF_InitCalc routines.
         !!    _Solution:_
         !!       -- Fix the ordering of the wave direction array in the Read_DataFile3D and Read_DataFile4D arrays.  Allow these arrays
         !!       span from -2pi to 2pi so that the directions are contiguous.  Add checks here to make sure the data read in is
         !!       shifted appropriately (i.e., data crosses the +/-pi boundary with values at 175 and -175 degrees, shift everything
         !!       so that it lies between 0 and 360 degrees so that it is contiguous).
         !!       -- In the _InitCalc routines, change the handling of the WaveDirMin and WaveDirMax and shift them when testing
         !!       across the +/-pi boundary.
         !!    _Reason not implimented:_
         !!       It takes so long for WAMIT to perform the calculations for the QTF that it is unlikely that any data where this is
         !!       a problem will arise before time is found to fix it.  Right now, time is not available.
         !!


   USE Waves, ONLY: WaveNumber
   USE WAMIT2_Types
   USE WAMIT_Interp
   USE NWTC_Library
   USE NWTC_FFTPACK

   IMPLICIT NONE

   PRIVATE

!   INTEGER(IntKi), PARAMETER                             :: DataFormatID = 1  !< Update this value if the data types change (used in WAMIT_Pack)
   TYPE(ProgDesc), PARAMETER                             :: WAMIT2_ProgDesc = ProgDesc( 'WAMIT2', '', '' )
                                                                              !< This holds the name of the program, version info, and date.

   REAL(SiKi), PARAMETER, PRIVATE                        :: OnePlusEps  = 1.0 + EPSILON(OnePlusEps)   ! The number slighty greater than unity in the precision of SiKi.


      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: WAMIT2_Init                           !< Initialization routine

   PUBLIC :: WAMIT2_CalcOutput                     !< Routine for computing outputs



   ! Derived types for data storage.
   !     There are a couple of good methods that would work for storing the data from the files.  We could:
   !     1. Store the data according to which file it came from.  Seems redundant
   !     2. Store the data in generic 3D or 4D structures (types).  This makes passing data messy.
   !     3. Store the data according to difference or sum method.  The difference method data becomes larger as it has to include
   !           the 3D and 4D options.  But it makes the passing of data much simpler to subroutines.
   !     4. Store the data as QTF and nonQTF method data.  The QTF method data would be only 4D, the nonQTF could be either 3D or
   !           or 4D.  Both the NewmanApp and MnDrift methods could use either the 4D (tossing most of it) or 3D.  This option
   !           becomes unatractive since we end up mixing the data for the full QTF diff method (DiffQTF) with the sum QTF (SumQTF)
   !
   !  In light of the options, it was arbitrarily chosen to use the 3rd method.  This is favors keeping the difference methods
   !  together, and the sum method as a separate entity.  This will streamline calls from the DiffQTF method to the MnDrift method
   !  (the first term of the DiffQTF equation is the MnDrift equation).


      !> This type is only used locally for holding data during the Initialization.  It will not be used outside of the WAMIT2_Init
      !! routine.  The 3D data is of the form F_k( WvFreq1, WvDir1, WvDir2, k ) where k is coordinate index to the load components
      !! (surge, sway, heave, roll, pitch, yaw).  DataSet is of size (NumFreq1,NumWvDir1,NumWvDir2,6*NumBodies). The LoadComponents array
      !! contains flags that will indicate which of the LoadComponents were read in from the data file.
      !!
      !! The DataMask array is of the same size as the DataSet matrix and corresponds to it.  Each point in DataMask that is true
      !! indicates that the datapoint in DataSet is an actual value.  Any missing points are set to 0 in DataSet and marked as
      !! FALSE in DataMask.  This is used for determining if the matrix is sparse or not. If any points are missing, then DataIsSparse
      !! is set to TRUE.  This will determine which algorithm we use for interpolation.
   TYPE, PRIVATE  :: W2_InitData3D_Type
      INTEGER(IntKi)                   :: NumWvFreq1           !< Number of frequencies in first  frequency direction set
      INTEGER(IntKi)                   :: NumWvDir1            !< Number of wave directions in first   wave direction set
      INTEGER(IntKi)                   :: NumWvDir2            !< Number of wave directions in second  wave direction set
      INTEGER(IntKi)                   :: NumBodies            !< Number of bodies in the file (based on highest LoadComponent listed)
      LOGICAL,          ALLOCATABLE    :: DataIsSparse(:)      !< Flag to indicate if the data is sparse or complete.  One per component direction
      LOGICAL,          ALLOCATABLE    :: LoadComponents(:)    !< Which load components actually exist in the input file
      COMPLEX(SiKi),    ALLOCATABLE    :: DataSet(:,:,:,:)     !< Storage for 3D data from the 2nd order WAMIT file
      LOGICAL,          ALLOCATABLE    :: DataMask(:,:,:,:)    !< Mask for knowing which data points are complete, which are missing
      REAL(SiKi),       ALLOCATABLE    :: WvFreq1(:)           !< (1:NumFreq1)  elements -- values correspond to index 1 of DataSet
      REAL(SiKi),       ALLOCATABLE    :: WvDir1(:)            !< (1:NumWvDir1) elements -- values correspond to index 2 of DataSet
      REAL(SiKi),       ALLOCATABLE    :: WvDir2(:)            !< (1:NumWvDir2) elements -- values correspond to index 3 of DataSet
   END TYPE W2_InitData3D_Type


      !> This type is only used locally for holding data during the Initialization.  It will not be used outside of the WAMIT2_Init
      !! routine.  The 4D data is of the form F_k( WvFreq1, WvFreq2, WvDir1, WvDir2, k ) where k is coordinate index to the load components
      !! (surge, sway, heave, roll, pitch, yaw).  DataSet is of size (NumFreq1,NumFreq2,NumWvDir1,NumWvDir2,6*NumBodies).  The LoadComponents
      !! array contains flags that will indicate which of the LoadComponents were read in from the data file.
      !!
      !! The DataMask array is of the same size as the DataSet matrix and corresponds to it.  Each point in DataMask that is true
      !! indicates that the datapoint in DataSet is an actual value.  Any missing points are set to 0 in DataSet and marked as
      !! FALSE in DataMask.  This is used for determining if the matrix is sparse or not.
   TYPE, PRIVATE  :: W2_InitData4D_Type
      INTEGER(IntKi)                   :: NumWvFreq1           !< Number of frequencies in first  frequency direction set
      INTEGER(IntKi)                   :: NumWvFreq2           !< Number of frequencies in second frequency direction set
      INTEGER(IntKi)                   :: NumWvDir1            !< Number of wave directions in first   wave direction set
      INTEGER(IntKi)                   :: NumWvDir2            !< Number of wave directions in second  wave direction set
      INTEGER(IntKi)                   :: NumBodies            !< Number of bodies in the file (based on highest LoadComponent listed)
      LOGICAL,          ALLOCATABLE    :: DataIsSparse(:)      !< Flag to indicate if the data is sparse or complete.  One per component direction
      LOGICAL                          :: WvFreqDiagComplete   !< Flag to indicate if the diagonal element is complete or not (Omega_m, Omega_m).
      LOGICAL                          :: IsSumForce           !< Flag to indicate that this is a sum type array.  Used in setting the F_{nm} term from F_{mn} term.
      LOGICAL,          ALLOCATABLE    :: LoadComponents(:)    !< Which load components actually exist in the input file
      COMPLEX(SiKi),    ALLOCATABLE    :: DataSet(:,:,:,:,:)   !< Storage for 4D data from the 2nd order WAMIT file
      LOGICAL,          ALLOCATABLE    :: DataMask(:,:,:,:,:)  !< Mask for knowing which data points are complete, which are missing
      REAL(SiKi),       ALLOCATABLE    :: WvFreq1(:)           !< (1:NumFreq1)  elements -- values correspond to index 1 of DataSet
      REAL(SiKi),       ALLOCATABLE    :: WvFreq2(:)           !< (1:NumFreq2)  elements -- values correspond to index 2 of DataSet
      REAL(SiKi),       ALLOCATABLE    :: WvDir1(:)            !< (1:NumWvDir1) elements -- values correspond to index 3 of DataSet
      REAL(SiKi),       ALLOCATABLE    :: WvDir2(:)            !< (1:NumWvDir2) elements -- values correspond to index 4 of DataSet
   END TYPE W2_InitData4D_Type


      !> This type is used to hold the data for difference frequency methods.  Only one difference method can be used at a time,
      !! but the data could be 3D or 4D depending on the method and input file used.  So in order to simplify keeping track of
      !! the data, this data structure was created to contain all difference data.  Note that since only one method can be used
      !! at a time, only one input file can be used at a time, and only one of Data3D or Data4D can be populated.  The flags
      !! DataIs3D and DataIs4D can be true at a time.
      !!
      !!  |  Algorithm    |  3D data   |  4D data  |
      !!  | :-----------: | :--------: | :-------: |
      !!  | MnDrift       |    yes     |    yes    |
      !!  | NewmanApp     |    yes     |    yes    |
      !!  | DiffQTF       |    NO      |    yes    |
      !!
   TYPE, PRIVATE  :: W2_DiffData_Type
      CHARACTER(1024)                  :: FileName             !< The filename the data came from
      LOGICAL                          :: DataIs3D             !< Only one of DataIs3D and DataIs4D can be true
      LOGICAL                          :: DataIs4D             !< Only one of DataIs3D and DataIs4D can be true
      TYPE(W2_InitData3D_Type)         :: Data3D               !< The 3D type from above
      TYPE(W2_InitData4D_Type)         :: Data4D               !< The 4D type from above
   END TYPE W2_DiffData_Type



      !> This type is used to hold the data for sum frequency method.  Only 4D data can be used for this method, so we
      !! will only allow 4D data within the type.  The only reason for actually making this type is to maintain consistency
      !! with the difference methods (helps with the programming implimentation).
      !!
      !! The flag DataIs4D is a carryover from the W2_DiffData_Type, and is preserved only for constency in the programming.
      !!
      !!  |  Algorithm    |  3D data   |  4D data  |
      !!  | :-----------: | :--------: | :-------: |
      !!  | SumQTF        |    NO      |    yes    |
      !!
   TYPE, PRIVATE  :: W2_SumData_Type
      CHARACTER(1024)                  :: FileName             !< The filename the data came from
      LOGICAL                          :: DataIs4D             !< Only true if Data4D is populated
      TYPE(W2_InitData4D_Type)         :: Data4D               !< The 4D type from above
   END TYPE W2_SumData_Type


CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> @brief
!!    This routine is called at the start of the simulation to perform initialization steps.
!!    The parameters that are set here are not changed during the simulation.
!!    The initial states and initial guess for the input are defined.
SUBROUTINE WAMIT2_Init( InitInp, p, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(WAMIT2_InitInputType),         INTENT(IN   )  :: InitInp              !< Input data for initialization routine
      TYPE(WAMIT2_ParameterType),         INTENT(  OUT)  :: p                    !< Parameters
      TYPE(WAMIT2_OutputType),            INTENT(  OUT)  :: y                    !< Initial system outputs (outputs are not calculated; only the output mesh is initialized)
      TYPE(WAMIT2_MiscVarType),           INTENT(  OUT)  :: m                    !< Initial misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat              !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg               !< Error message if ErrStat /= ErrID_None


         ! Local Variables
      INTEGER(IntKi)                                     :: IBody                !< Counter for current body
      INTEGER(IntKi)                                     :: ThisDim              !< Counter to currrent dimension
      INTEGER(IntKi)                                     :: Idx                  !< Generic counter
      REAL(R8Ki)                                         :: theta(3)             !< rotation about z for the current body (0 about x,y)
      REAL(R8Ki)                                         :: orientation(3,3)     !< Orientation matrix for orientation of the current body
      REAL(ReKi)                                         :: XYZloc(3)            !< Starting position of this WAMIT2 body


         ! QTF storage -- from the data files that are read in
      TYPE(W2_DiffData_Type)                             :: MnDriftData          !< Data storage for the Mean Drift method
      TYPE(W2_DiffData_Type)                             :: NewmanAppData        !< Data storage for the Newman's Approximation method
      TYPE(W2_DiffData_Type)                             :: DiffQTFData          !< Data storage for the full difference QTF method
      TYPE(W2_SumData_Type)                              :: SumQTFData           !< Data storage for the full sum QTF method

         ! Force arrays
      REAL(SiKi),    ALLOCATABLE                         :: MnDriftForce(:)      !< MnDrift force array.   Constant for all time.  First index is force component
      REAL(SiKi),    ALLOCATABLE                         :: NewmanAppForce(:,:)  !< NewmanApp force array.  Index 1: Time,    Index 2: force component
      REAL(SiKi),    ALLOCATABLE                         :: DiffQTFForce(:,:)    !< DiffQTF force array.    Index 1: Time,    Index 2: force component
      REAL(SiKi),    ALLOCATABLE                         :: SumQTFForce(:,:)     !< SumQTF force array.     Index 1: Time,    Index 2: force component

         ! Temporary error trapping variables
      INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary variable for holding the error status  returned from a CALL statement
      CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary variable for holding the error message returned from a CALL statement
      CHARACTER(*), PARAMETER                            :: RoutineName = 'WAMIT2_Init'


         !> ### Subroutine contents

         !> Initialize Error handling variables

      ErrStat     = ErrID_None
      ErrMsg      = ""


         !> Initialize the data storage
      MnDriftData%DataIs3D    = .FALSE.
      MnDriftData%DataIs4D    = .FALSE.
      MnDriftData%Filename    = ''

      NewmanAppData%DataIs3D  = .FALSE.
      NewmanAppData%DataIs4D  = .FALSE.
      NewmanAppData%Filename  = ''

      DiffQTFData%DataIs3D    = .FALSE.
      DiffQTFData%DataIs4D    = .FALSE.
      DiffQTFData%Filename    = ''

      SumQTFData%DataIs4D     = .FALSE.
      SumQTFData%Filename     = ''


      !-----------------------------------------------------------------------------
      !> Before attempting to do any real calculations, we first check what was
      !! passed in through _InitInp_ to make sure it makes sense.  That routine will
      !! then copy over the relevant information that should be kept in parameters
      !! (_p_).
      !!
      !! _InitInp_ will also check the flags, existence of files, and set flags
      !! accordingly.
      !-----------------------------------------------------------------------------

         !-------------------------------------------------------------------------------------------------------------
         !> 1. Check the data file related values (_MnDrift_, _MnDriftF_ etc). Also copy over important things from _InitInp_ to _p_ and _InitOut_.

      CALL CheckInitInput( InitInp, p, MnDriftData, NewmanAppData, DiffQTFData, SumQTFData, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         !-------------------------------------------------------------------------------------------------------------
         !> 2. Now read in the data files and store the data and information about it using _Read_DataFile3D_ or
         !!    _Read_DataFile4D_.  Flags should have been set in the _CheckInitInput_ subroutine to indicate
         !!    what type of data is going to be read in.  The consistency of the WAMIT 2nd order data files is checked
         !!    during the reading in subroutines.

         !> If the MnDrift method will be used, read in the data for it.
      IF ( p%MnDriftF ) THEN
         IF ( MnDriftData%DataIs3D ) THEN
            CALL Read_DataFile3D( InitInp, TRIM(MnDriftData%Filename), MnDriftData%Data3D, ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         ELSEIF ( MnDriftData%DataIs4D ) THEN
            MnDriftData%Data4D%IsSumForce = .FALSE.
            CALL Read_DataFile4D( InitInp, TRIM(MnDriftData%Filename), MnDriftData%Data4D, ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         ELSE
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  MnDrift method flags incorrectly set by '// &
                  'CheckInitInput subroutine.', ErrStat, ErrMsg, RoutineName )
         ENDIF
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !> If the NewmanApp method will be used, read in the data for it.
      IF ( p%NewmanAppF ) THEN
         IF ( NewmanAppData%DataIs3D ) THEN
            CALL Read_DataFile3D( InitInp, TRIM(NewmanAppData%Filename), NewmanAppData%Data3D, ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         ELSEIF ( NewmanAppData%DataIs4D ) THEN
            NewmanAppData%Data4D%IsSumForce = .FALSE.
            CALL Read_DataFile4D( InitInp, TRIM(NewmanAppData%Filename), NewmanAppData%Data4D, ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         ELSE
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  NewmanApp method flags incorrectly set by '// &
                  'CheckInitInput subroutine.', ErrStat, ErrMsg, RoutineName )
         ENDIF
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !> If the DiffQTF method will be used, read in the data for it.  This can only be 4D data even though DiffQTF can hold
         !! 3D data.  This was an intentional choice so that DiffQTFData could be passed to the MnDrift method (the first term
         !! of the DiffQTF method equations is the MnDrift term).
      IF ( p%DiffQTFF ) THEN
         IF ( DiffQTFData%DataIs3D ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  DiffQTF method flags incorrectly set by '// &
                  'CheckInitInput subroutine. 3D data cannot be used in the DiffQTF method.', ErrStat, ErrMsg, RoutineName )
         ELSEIF ( DiffQTFData%DataIs4D ) THEN
            DiffQTFData%Data4D%IsSumForce = .FALSE.
            CALL Read_DataFile4D( InitInp, TRIM(DiffQTFData%Filename), DiffQTFData%Data4D, ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         ELSE
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  DiffQTF method flags incorrectly set by '// &
                                          'CheckInitInput subroutine.', ErrStat, ErrMsg, RoutineName )
         ENDIF
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !> If the SumQTF method will be used, read in the data for it.  This can only be 4D data (SumQTFData only holds 4D).
      IF ( p%SumQTFF ) THEN
         IF ( SumQTFData%DataIs4D ) THEN
            SumQTFData%Data4D%IsSumForce = .TRUE.
            CALL Read_DataFile4D( InitInp, TRIM(SumQTFData%Filename), SumQTFData%Data4D, ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         ELSE
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  SumQTF method flags incorrectly set by '// &
                  'CheckInitInput subroutine.', ErrStat, ErrMsg, RoutineName )
         ENDIF
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF



         !-------------------------------------------------------------------------------------------------------------
         !> 3. Now check the data to ensure that all the dimensions that were requested acually exist in the file.  At
         !!    this point, the MnDriftData and NewmanApp data will be either 3D or 4D, but not both.
         !!
         !> Check the MnDrift data: check both in case we have already copied 4D data into the 3D.  Will decide later which is used.
      IF ( p%MnDriftF ) THEN
         IF ( MnDriftData%DataIs3D ) THEN
               ! Check the dimensions used.  The LoadComponents(Idx) flag will be set to .TRUE. if data was found in the file
            DO IBody=1,MnDriftData%Data3D%NumBodies
               DO ThisDim=1,6
                  Idx =  (IBody-1)*6+ThisDim
                  IF ( p%MnDriftDims(ThisDim) .AND. ( .NOT. MnDriftData%Data3D%LoadComponents(Idx) ) ) &
                     CALL SetErrStat( ErrID_Warn, ' '//TRIM(MnDriftData%Filename)//' does not contain information for the '// &
                           TRIM(Num2LStr(ThisDim))//' force component for the MnDrift method.  Setting this component to zero.',    &
                           ErrStat,ErrMsg,RoutineName)
               ENDDO
            ENDDO
         ELSE IF ( MnDriftData%DataIs4D ) THEN
               ! Check the dimensions used.  The LoadComponents(Idx) flag will be set to .TRUE. if data was found in the file
            DO IBody=1,MnDriftData%Data4D%NumBodies
               DO ThisDim=1,6
                  Idx =  (IBody-1)*6+ThisDim
                  IF ( p%MnDriftDims(ThisDim) .AND. ( .NOT. MnDriftData%Data4D%LoadComponents(Idx) )  )&
                     CALL SetErrStat( ErrID_Warn, ' '//TRIM(MnDriftData%Filename)//' does not contain information for the '// &
                           TRIM(Num2LStr(ThisDim))//' force component for the MnDrift method.  Setting this component to zero.',    &
                           ErrStat,ErrMsg,RoutineName)
               ENDDO
            ENDDO
         ELSE  ! We didn't find any data to use...
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  MnDrift flag is set, but no data has been read in.',ErrStat,ErrMsg, RoutineName)
         ENDIF
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !> Check the NewmanApp data: check both in case we have already copied 4D data into the 3D.  Will decide later which is used.
      IF ( p%NewmanAppF ) THEN
         IF ( NewmanAppData%DataIs3D ) THEN
               ! Check the dimensions used.  The LoadComponents(Idx) flag will be set to .TRUE. if data was found in the file
            DO IBody=1,NewmanAppData%Data3D%NumBodies
               DO ThisDim=1,6
                  Idx =  (IBody-1)*6+ThisDim
                  IF ( p%NewmanAppDims(ThisDim) .AND. ( .NOT. NewmanAppData%Data3D%LoadComponents(Idx) ) ) &
                     CALL SetErrStat( ErrID_Warn, ' '//TRIM(NewmanAppData%Filename)//' does not contain information for the '// &
                           TRIM(Num2LStr(ThisDim))//' force component for the NewmanApp method.  Setting this component to zero.',    &
                           ErrStat,ErrMsg,RoutineName)
               ENDDO
            ENDDO
         ELSE IF ( NewmanAppData%DataIs4D ) THEN
               ! Check the dimensions used.  The LoadComponents(Idx) flag will be set to .TRUE. if data was found in the file
            DO IBody=1,NewmanAppData%Data4D%NumBodies
               DO ThisDim=1,6
                  Idx =  (IBody-1)*6+ThisDim
                  IF ( p%NewmanAppDims(ThisDim) .AND. ( .NOT. NewmanAppData%Data4D%LoadComponents(Idx) ) ) &
                     CALL SetErrStat( ErrID_Warn, ' '//TRIM(NewmanAppData%Filename)//' does not contain information for the '// &
                           TRIM(Num2LStr(ThisDim))//' force component for the NewmanApp method.  Setting this component to zero.',    &
                           ErrStat,ErrMsg,RoutineName)
               ENDDO
            ENDDO
         ELSE
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  NewmanApp flag is set, but no data has been read in.',ErrStat,ErrMsg, RoutineName)
         ENDIF
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !> Check the DiffQTF data: Don't check the 3D data.  We may have copied 4D into it for the MnDrift term.
      IF ( p%DiffQTFF ) THEN
         IF ( DiffQTFData%DataIs4D ) THEN
               ! Check the dimensions used.  The LoadComponents(Idx) flag will be set to .TRUE. if data was found in the file
            DO IBody=1,DiffQTFData%Data4D%NumBodies
               DO ThisDim=1,6
                  Idx =  (IBody-1)*6+ThisDim
                  IF ( p%DiffQTFDims(ThisDim) .AND. ( .NOT. DiffQTFData%Data4D%LoadComponents(Idx) ) ) &
                     CALL SetErrStat( ErrID_Warn, ' '//TRIM(DiffQTFData%Filename)//' does not contain information for the '// &
                           TRIM(Num2LStr(ThisDim))//' force component for the DiffQTF method.  Setting this component to zero.',    &
                           ErrStat,ErrMsg,RoutineName)
               ENDDO
            ENDDO
         ELSE
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  DiffQTF flag is set, but no data has been read in.',ErrStat,ErrMsg, RoutineName)
         ENDIF
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !> Check the SumQTF data
      IF ( p%SumQTFF ) THEN
         IF ( SumQTFData%DataIs4D ) THEN
               ! Check the dimensions used.  The LoadComponents(Idx) flag will be set to .TRUE. if data was found in the file
            DO IBody=1,SumQTFData%Data4D%NumBodies
               DO ThisDim=1,6
                  Idx =  (IBody-1)*6+ThisDim
                  IF ( p%SumQTFDims(ThisDim) .AND. ( .NOT. SumQTFData%Data4D%LoadComponents(Idx) ) ) &
                     CALL SetErrStat( ErrID_Warn, ' '//TRIM(SumQTFData%Filename)//' does not contain information for the '// &
                           TRIM(Num2LStr(ThisDim))//' force component for the SumQTF method.  Setting this component to zero.',    &
                           ErrStat,ErrMsg,RoutineName)
               ENDDO
            ENDDO
         ELSE
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  SumQTF flag is set, but no data has been read in.',ErrStat,ErrMsg, RoutineName)
         ENDIF
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !-------------------------------------------------------------------------------------------------------------
         !> 4. At this point, all the data from the WAMIT second order output files should be loaded.  We can now proceed
         !!    to perform the calculations for each of the methods.  For the call to the method, we pass the dataset
         !!    corresponding to each method.  The calculation methods for MnDrift and NewmanApp can handle both 3D and 4D
         !!    data, so the entire MnDriftData and NewmanAppData sets are handed to the calculation subroutines.
         !!
         !! @note The frequency range and wave direction range present in the QTF information that was read in has not
         !!    been checked.  The frequency range needed by the Newman approximation is different than that of the other
         !!    methods.  Therefore the decision has been made to check and apply the frequency check and wave directionality
         !!    checks at the start of the calculation subroutine for each method.
         !!
         !!    For the DiffQTF and SumQTF methods, only 4D data can be used in the calculations.  Again, the entire dataset
         !!    for each is passed into the respective calculation subroutine.  In the case of the DiffQTFData, only the 4D
         !!    data portion should be populated above.  If 3D data is present in it, the subroutine will exit with a
         !!    warning about a programming error.  For the SumQTF method, the data can only be 4D.
         !!
         !! @note The wamit2:mndrift_initcalc routine will be called from wamit2::diffqtf_initcalc.  The first term of the
         !!       DiffQTF method is the MnDrift method.
         !!
         !! For the MnDrift and NewmanApp methods, it is possible to convert 4D data where the diagonal is complete but
         !! off diagonal values are not.  The data in this case will have to be ported over to a 3D storage.  This will
         !! be handled in the MnDrift_InitCalc and NewmanApp_InitCalc functions.
         !!
         !! @note The MnDrift and NewmanApp can use the 4D data if the diagonal is not complete, but the full QTF is.
         !!       This would only occur when the stepsize in \f$ \omega_1 \f$ and \f$ \omega_2 \f$ are not the same
         !!       resulting in no cases where \f$ \omega_1 == \omega_2 \f$.



         !> If the MnDrift method will be used, call the subroutine to calculate the force time series
      IF ( p%MnDriftF ) THEN

            ! Tell our nice users what is about to happen that may take a while:
         CALL WrScr ( ' Calculating second order mean drift force.' )

         CALL MnDrift_InitCalc( InitInp, p, MnDriftData, MnDriftForce, ErrMsgTmp, ErrStatTmp )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !> If the NewmanApp method will be used, call the subroutine to calculate the force time series
      IF ( p%NewmanAppF ) THEN
            ! Tell our nice users what is about to happen that may take a while:
         CALL WrScr ( " Calculating second order difference-frequency force using the Newman's approximation." )

         CALL NewmanApp_InitCalc( InitInp, p, NewmanAppData, NewmanAppForce, ErrMsgTmp, ErrStatTmp )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF



         !> If the DiffQTF method will be used, call the subroutine to calculate the force time series
         !! Note that the MnDrift calculation is included.
      IF ( p%DiffQTFF ) THEN
            ! Tell our nice users what is about to happen that may take a while:
         CALL WrScr ( ' Calculating second order difference-frequency force using the full quadratic transfer function.' )

         CALL DiffQTF_InitCalc( InitInp, p, DiffQTFData, DiffQTFForce, ErrMsgTmp, ErrStatTmp )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !> If the SumQTF method will be used, call the subroutine to calculate the force time series
      IF ( p%SumQTFF ) THEN
            ! Tell our nice users what is about to happen that may take a while:
         CALL WrScr ( ' Calculating second order sum-frequency force using the full quadratic transfer function.' )

         CALL SumQTF_InitCalc( InitInp, p, SumQTFData, SumQTFForce, ErrMsgTmp, ErrStatTmp )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
      ENDIF


         !----------------------------------------------------------------------
         !> Copy output forces over to parameters as needed.
         !----------------------------------------------------------------------

         ! Initialize the second order force to zero.
      p%WaveExctn2 = 0.0_SiKi


         ! Difference method data.  Only one difference method can be calculated at a time.
      IF ( p%MnDriftF ) THEN

         DO IBody=1,p%NBody        ! Loop through load components. Set ones that are calculated.
            DO ThisDim=1,6
               Idx =  (IBody-1)*6+ThisDim
               IF ( p%MnDriftDims(ThisDim) ) THEN
                  p%WaveExctn2(:,Idx) = MnDriftForce(Idx)
               ENDIF
            ENDDO
         ENDDO

      ELSE IF ( p%NewmanAppF ) THEN

         DO IBody=1,p%NBody        ! Loop through load components. Set ones that are calculated.
            DO ThisDim=1,6
               Idx =  (IBody-1)*6+ThisDim
               IF ( p%NewmanAppDims(ThisDim) ) THEN
                  p%WaveExctn2(:,Idx) = NewmanAppForce(:,Idx)
               ENDIF
            ENDDO
         ENDDO

      ELSE IF ( p%DiffQTFF ) THEN

         DO IBody=1,p%NBody        ! Loop through load components. Set ones that are calculated.
            DO ThisDim=1,6
               Idx =  (IBody-1)*6+ThisDim
               IF ( p%DiffQTFDims(ThisDim) ) THEN
                  p%WaveExctn2(:,Idx) = DiffQTFForce(:,Idx)
               ENDIF
            ENDDO
         ENDDO

      ENDIF


         ! Summation method
      IF ( p%SumQTFF ) THEN

         DO IBody=1,p%NBody        ! Loop through load components. Set ones that are calculated.
            DO ThisDim=1,6
               Idx =  (IBody-1)*6+ThisDim
               IF ( p%SumQTFDims(ThisDim) ) THEN
                  p%WaveExctn2(:,Idx) = p%WaveExctn2(:,Idx) + SumQTFForce(:,Idx)
               ENDIF
            ENDDO
         ENDDO

      ENDIF


         ! Deallocate the force arrays since we are done with them.  Note that the MnDrift force array is
         ! not deallocated since it is not time dependent.
      IF (ALLOCATED(NewmanAppForce))         DEALLOCATE(NewmanAppForce)
      IF (ALLOCATED(DiffQTFForce))           DEALLOCATE(DiffQTFForce)
      IF (ALLOCATED(SumQTFForce))            DEALLOCATE(SumQTFForce)



         !----------------------------------------------------------------------
         !> 5. Create the mesh used for exporting the 2nd order forces from the
         !!    WAMIT2_CalcOuput routine at each timestep.  Also set the outputs
         !!    for each timestep.
         !----------------------------------------------------------------------

         ! Create the input and output meshes associated with lumped loads
      CALL MeshCreate(  BlankMesh         = y%Mesh           , &
                        IOS               = COMPONENT_OUTPUT , &
                        Nnodes            = p%NBody          , &
                        ErrStat           = ErrStatTmp       , &
                        ErrMess           = ErrMsgTmp        , &
                        Force             = .TRUE.           , &
                        Moment            = .TRUE.)

      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF

      DO IBody = 1,p%NBody

            ! Set orientation and position for each body in mesh
         theta = (/ 0.0_R8Ki, 0.0_R8Ki, 0.0_R8Ki /)
         orientation = EulerConstruct(theta)
         XYZloc      = (/InitInp%PtfmRefxt(IBody), InitInp%PtfmRefyt(IBody), InitInp%PtfmRefzt(IBody)/)

            ! Create the node on the mesh
         CALL MeshPositionNode (y%Mesh, IBody, XYZloc, ErrStatTmp, ErrMsgTmp, orientation )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)

            ! Create the mesh element
         CALL MeshConstructElement (  y%Mesh, ELEMENT_POINT, ErrStatTmp, ErrMsgTmp, IBody )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      ENDDO

      CALL MeshCommit ( y%Mesh, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      END IF


      y%Mesh%RemapFlag  = .TRUE.


         !----------------------------------------------------------------------
         !> 6. Set zero values for unused outputs.  This is mostly so that the
         !!    compiler does not complain.  Also set misc vars
         !----------------------------------------------------------------------
      CALL AllocAry( m%LastIndWave, p%NBody, 'm%LastIndWave', ErrStatTmp, ErrMsgTmp)
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      m%LastIndWave              = 1_IntKi
      call AllocAry(m%F_Waves2, 6*p%NBody, 'm%F_Waves2', ErrStatTmp, ErrMsgTmp)
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)

         ! Cleanup remaining stuff
      CALL CleanUp()

      RETURN

   CONTAINS


   !-------------------------------------------------------------------------------------------------------------------------------
   !> This subroutine cleans up anything that is still allocated.
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CleanUp()

      IF (ALLOCATED(MnDriftForce)) DEALLOCATE(MnDriftForce)
      IF (ALLOCATED(NewmanAppForce)) DEALLOCATE(NewmanAppForce)
      IF (ALLOCATED(DiffQTFForce)) DEALLOCATE(DiffQTFForce)
      IF (ALLOCATED(SumQTFForce)) DEALLOCATE(SumQTFForce)
   
      CALL Destroy_InitData3D(   MnDriftData%Data3D )
      CALL Destroy_InitData4D(   MnDriftData%Data4D )
      CALL Destroy_InitData3D( NewmanAppData%Data3D )
      CALL Destroy_InitData4D( NewmanAppData%Data4D )
      CALL Destroy_InitData3D(   DiffQTFData%Data3D )
      CALL Destroy_InitData4D(   DiffQTFData%Data4D )
      CALL Destroy_InitData4D(    SumQTFData%Data4D )

   END SUBROUTINE CleanUp

END SUBROUTINE WAMIT2_Init


   !-------------------------------------------------------------------------------------------------------------------------------
   !> This subroutine calculates the force time series using the MnDrift calculation.
   !! The data is stored in either 3D or 4D arrays depending on the file type used.
   !! At each step in the summation of the mth term, a call is made to the 3D or 4D interpolation algorithm to find the value of
   !! \f$ F^-_k(\omega_m, \omega_n) \f$ corresponding to the \f$ Z(\omega_m) \f$ term in the complex wave amplitude, _WaveElevC_.
   !! The limits of \f$ \omega_{lo-d} :math:`\le` \omega_m :math:`\le` \omega_{hi-d} \f$ are imposed during the summation with values
   !! outside this range set to zero.
   !!
   !! For multi-directional waves where the equal energy discretization is used, each frequency has a single wave direction
   !! associated with it. Since the mean drift force calculation only involves summing over terms involving only a single frequency
   !! at a time, only a single wave direction is involved at each step. So, if all the diagonal elements of the 4D matrix where
   !! \f$ \omega_1 = \omega_2 \f$ and \f$ \beta_1 = \beta_2 \f$ were present, it would be possible to simplify the 4D interpolation
   !! required to two dimensional interpolation for this particular case.  However, since the same interpolation routine is used for
   !! all the 4D data handling, it is programatically simpler to use only one interpolation algorithm for the 4D data and incur a
   !! slight penalty at the initialization step where this routine is called.  Overall it is a small price to pay for maintainable
   !! code.
   !!
   !! The single summation equation used here is given by
   !! \f$  {F_{{ex}~k}^{{-}(2)}} = \Re \left( \sum\limits_{m=0}^{N/2}
   !!                                        {a_m} {a_m^*} F_k^{-}(\omega_m, \beta_m)\right) \cdot\Delta\omega  \f$
   !!      for \f$\quad k=1,2,\ldots,6,      \f$
   !!
   !! where \f$ k \f$ indicates the index to the load component,  \f$ {F_{{ex}~k}^{{-}(2)}} \f$ is the resulting time
   !! independent mean drift force, and \f$ a_m \f$ and \f$ a_m^* \f$ are the complex wave amplitude and its complex conjugate for
   !! the \f$ m^{th} \f$ frequency.  _Note the lack of time dependence in this equation:_ the mean drift is the average drift
   !! force over the entire simulation. Note also that \f$ F_k^{-} (\omega_m, \beta_m) \f$ is the dimensionalized real valued
   !! diagonal of the QTF read from the WAMIT file and interpolated for the \f$ m^{th} \f$ wave frequency.  Note that this is
   !! is a numerical integral, so the \f$ \Delta\omega \f$ term is necessary.
   !!
   !! @note The mean drift is a static value representing the mean drift force for the entire simulation.  Therefore, any offset
   !!       in the body location can be ignored since it cancels out (it is in the complex part that is ignored in the summation).
   !!       The orientation of the platform is therefore handled after the summation.
   !!
   !! Since the frequency range of the QTF has not yet been checked, we will do that now.
   !!
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE MnDrift_InitCalc( InitInp, p, MnDriftData, MnDriftForce, ErrMsg, ErrStat )

      IMPLICIT NONE

      TYPE(WAMIT2_InitInputType),         INTENT(IN   )  :: InitInp              !< Input data for initialization routine
      TYPE(WAMIT2_ParameterType),         INTENT(IN   )  :: p                    !< Parameters
      TYPE(W2_DiffData_Type),             INTENT(INOUT)  :: MnDriftData          !< Data storage for the MnDrift method.  Set to INOUT in case we need to convert 4D to 3D
      REAL(SiKi),  ALLOCATABLE,           INTENT(  OUT)  :: MnDriftForce(:)      !< Force data.  Index 1 is the force component.  Constant for all time.
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat

         ! Local Variables
      CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary error message for calls
      INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
!      REAL(SiKi)                                         :: TmpReal1             !< Temporary real
!      REAL(SiKi)                                         :: TmpReal2             !< Temporary real
      LOGICAL                                            :: TmpFlag              !< Temporary logical flag
      INTEGER(IntKi)                                     :: ThisDim              !< Generic counter for dimension
      INTEGER(IntKi)                                     :: IBody                !< Index to which body we are on
      INTEGER(IntKi)                                     :: Idx                  !< Index to the full set of 6*NBody
      INTEGER(IntKi)                                     :: J                    !< Generic counter
!      INTEGER(IntKi)                                     :: K                    !< Generic counter
      CHARACTER(*), PARAMETER                            :: RoutineName = 'MnDrift_InitCalc'


         ! Wave information and QTF temporary
      COMPLEX(SiKi)                                      :: QTF_Value            !< Temporary complex number for QTF
      COMPLEX(SiKi)                                      :: aWaveElevC           !< Wave elevation of current frequency component.  NStepWave2 normalization is removed.
      REAL(ReKi)                                         :: Omega1               !< Wave frequency of this component

         ! Interpolation routine indices and value to search for, and smaller array to pass
      INTEGER(IntKi)                                     :: LastIndex3(3)        !< Last used index for searching in the interpolation algorithms
      INTEGER(IntKi)                                     :: LastIndex4(4)        !< Last used index for searching in the interpolation algorithms
      REAL(SiKi)                                         :: RotateZdegOffset     !< Offset to wave heading (NBodyMod==2 only)
      REAL(SiKi)                                         :: RotateZMatrixT(2,2)  !< The transpose of rotation in matrix form for rotation about z (from global to local)
      REAL(SiKi)                                         :: Coord3(3)            !< The (omega1,beta1,beta2) coordinate we want in the 3D dataset
      REAL(SiKi)                                         :: Coord4(4)            !< The (omega1,omega2,beta1,beta2) coordinate we want in the 4D dataset
      COMPLEX(SiKi),ALLOCATABLE                          :: TmpData3D(:,:,:)     !< Temporary 3D array we put the 3D data into (minus the load component indice)
      COMPLEX(SiKi),ALLOCATABLE                          :: TmpData4D(:,:,:,:)   !< Temporary 4D array we put the 4D data into (minus the load component indice)


         ! Initialize a few things
      ErrMsg      = ''
      ErrStat     = ErrID_None

         ! Initialize resulting forces
      ALLOCATE( MnDriftForce(6*p%NBody), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) THEN
         CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for the resulting mean drift force '// &
                                             'of the 2nd order force.',ErrStat, ErrMsg, RoutineName)
         RETURN
      ENDIF
      MnDriftForce = 0.0_SiKi ! initialize this subroutine return value


         !> 1. Check the data to see if low cutoff on the difference frequency is 0.  If it is above zero, that implies no mean drift
         !!    term since \f$ \omega_1=\omega_2 \f$

      IF ( InitInp%WvLowCOffD > 0.0_SiKi )  THEN
         CALL SetErrStat( ErrID_Warn, ' WvLowCOffD > 0.0, so no mean drift term is calculated (the mean drift uses only the equal '//&
                        'frequency terms of the QTF).  Setting the mean drift force to 0.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF



         !> 2. Check the data to see if the wave frequencies are present in the QTF data.  Since the mean drift term only uses
         !!    frequencies where \f$ \omega_1=\omega_2 \f$, the data read in from the files must contain the full range of frequencies
         !!    present in the waves.

      IF ( MnDriftData%DataIs3D ) THEN

            ! Check the low frequency cutoff
         IF ( MINVAL( MnDriftData%Data3D%WvFreq1 ) > InitInp%WvLowCOffD ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(MnDriftData%Data3D%WvFreq1)))// &
                           ' rad/s for first wave period) data in '//TRIM(MnDriftData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOffD.',ErrStat,ErrMsg,RoutineName)
         ENDIF

            ! Check the high frequency cutoff -- using the Difference high frequency cutoff.  The first order high frequency
            ! cutoff is typically too high for this in most cases.
         IF ( (MAXVAL(MnDriftData%Data3D%WvFreq1 ) < InitInp%WvHiCOffD) ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(MnDriftData%Data3D%WvFreq1)))// &
                           ' rad/s for first wave period) data in '//TRIM(MnDriftData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOffD.',ErrStat,ErrMsg,RoutineName)
         ENDIF

      ELSE IF ( MnDriftData%DataIs4D ) THEN   ! only check if not 3D data. If there is 3D data, we default to using it for calculations

             ! Check the low frequency cutoff
         IF ( MINVAL( MnDriftData%Data4D%WvFreq1 ) > InitInp%WvLowCOffD ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(MnDriftData%Data4D%WvFreq1)))// &
                           ' rad/s first wave period) data in '//TRIM(MnDriftData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOff.',ErrStat,ErrMsg,RoutineName)
         ENDIF
         IF ( MINVAL( MnDriftData%Data4D%WvFreq2 ) > InitInp%WvLowCOff ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(MnDriftData%Data4D%WvFreq2)))// &
                           ' rad/s for second wave period) data in '//TRIM(MnDriftData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOffD.',ErrStat,ErrMsg,RoutineName)
         ENDIF

            ! Check the high frequency cutoff -- using the Difference high frequency cutoff.  The first order high frequency
            ! cutoff is typically too high for this in most cases.
         IF ( (MAXVAL(MnDriftData%Data4D%WvFreq1) < InitInp%WvHiCOffD) ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(MnDriftData%Data4D%WvFreq1)))// &
                           ' rad/s for first wave period) data in '//TRIM(MnDriftData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOffD.',ErrStat,ErrMsg,RoutineName)
         ENDIF
         IF ( (MAXVAL(MnDriftData%Data4D%WvFreq2) < InitInp%WvHiCOffD) ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(MnDriftData%Data4D%WvFreq1)))// &
                           ' rad/s second wave period) data in '//TRIM(MnDriftData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOffD.',ErrStat,ErrMsg,RoutineName)
         ENDIF

      ELSE
            ! This is a catastrophic issue.  We should not have called this routine without data that is usable for the MnDrift calculation
         CALL SetErrStat( ErrID_Fatal, ' Mean drift calculation called without data.',ErrStat,ErrMsg,RoutineName)
      ENDIF

      IF ( ErrStat >= AbortErrLev )    RETURN



         !> 3. Check the data to see if the wave directions are present.  May need to adjust for the boundary at +/- PI
      IF ( MnDriftData%DataIs3D ) THEN

            ! If we are using multidirectional waves, then we should have more than 1 wave direction in the WAMIT file.
         IF ( InitInp%WaveMultiDir .AND. (MnDriftData%Data3D%NumWvDir1 == 1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(MnDriftData%Filename)//' only contains one wave '// &
                        'direction at '//TRIM(Num2LStr(MnDriftData%Data3D%WvDir1(1)))//' degrees (first wave direction). '// &
                        'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                        ErrStat,ErrMsg,RoutineName)
         ELSE IF ( InitInp%WaveMultiDir .AND. (MnDriftData%Data3D%NumWvDir2 == 1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(MnDriftData%Filename)//' only contains one wave '// &
                        'direction at '//TRIM(Num2LStr(MnDriftData%Data3D%WvDir2(1)))//' degrees (second wave direction). '// &
                        'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                        ErrStat,ErrMsg,RoutineName)
         ELSE

               ! See Known Issues #1 at the top of this file.  There may be problems if the data spans the +/- Pi boundary.  For
               ! now (since time is limited) we will issue a warning if any of the wave directions for multidirectional waves
               ! or data from the WAMIT file for the wavedirections is close to the +/-pi boundary (>150 degrees, <-150 degrees),
               ! we will issue a warning.
            IF ( (InitInp%WaveDirMin > 150.0_SiKi) .OR. (InitInp%WaveDirMax < -150.0_SiKi) .OR. &
                 (minval(MnDriftData%data3d%WvDir1) > 150.0_SiKi) .OR.  (maxval(MnDriftData%data3d%WvDir1) < -150.0_SiKi) .OR. &
                 (minval(MnDriftData%data3d%WvDir2) > 150.0_SiKi) .OR.  (maxval(MnDriftData%data3d%WvDir2) < -150.0_SiKi) ) THEN
               CALL SetErrStat( ErrID_Warn,' There may be issues with how the wave direction data is handled when the wave '// &
                                          'direction of interest is near the +/- 180 direction.  This is a known issue with '// &
                                          'the WAMIT2 module that has not yet been addressed.',ErrStat,ErrMsg,RoutineName)
            ENDIF

               ! Now check the limits for the first wave direction
               !   --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
            IF ( InitInp%WaveDirMin < MINVAL(MnDriftData%Data3D%WvDir1) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                     'found in the WAMIT data file '//TRIM(MnDriftData%Filename)//' for the first wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF
            IF ( InitInp%WaveDirMax > MAXVAL(MnDriftData%Data3D%WvDir1) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                     'found in the WAMIT data file '//TRIM(MnDriftData%Filename)//' for the first wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF


               ! Now check the limits for the second wave direction
               !   --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
            IF ( InitInp%WaveDirMin < MINVAL(MnDriftData%Data3D%WvDir2) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                     'found in the WAMIT data file '//TRIM(MnDriftData%Filename)//' for the second wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF
            IF ( InitInp%WaveDirMax > MAXVAL(MnDriftData%Data3D%WvDir2) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                     'found in the WAMIT data file '//TRIM(MnDriftData%Filename)//' for the second wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF

         ENDIF

      ELSEIF ( MnDriftData%DataIs4D ) THEN

            ! If we are using multidirectional waves, then we should have more than 1 wave direction in the WAMIT file.
         IF ( InitInp%WaveMultiDir .AND. (MnDriftData%Data4D%NumWvDir1 == 1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(MnDriftData%Filename)//' only contains one wave '// &
                        'direction at '//TRIM(Num2LStr(MnDriftData%Data4D%WvDir1(1)))//' degrees (first wave direction). '// &
                        'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                        ErrStat,ErrMsg,RoutineName)
         ELSE IF ( InitInp%WaveMultiDir .AND. (MnDriftData%Data4D%NumWvDir2 == 1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(MnDriftData%Filename)//' only contains one wave '// &
                        'direction at '//TRIM(Num2LStr(MnDriftData%Data4D%WvDir2(1)))//' degrees (second wave direction). '// &
                        'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                        ErrStat,ErrMsg,RoutineName)
         ELSE

               ! See Known Issues #1 at the top of this file.  There may be problems if the data spans the +/- Pi boundary.  For
               ! now (since time is limited) we will issue a warning if any of the wave directions for multidirectional waves
               ! or data from the WAMIT file for the wavedirections is close to the +/-pi boundary (>150 degrees, <-150 degrees),
               ! we will issue a warning.
            IF ( (InitInp%WaveDirMin > 150.0_SiKi) .OR. (InitInp%WaveDirMax < -150.0_SiKi) .OR. &
                 (MINVAL(MnDriftData%Data4D%WvDir1) > 150.0_SiKi) .OR.  (MAXVAL(MnDriftData%Data4D%WvDir1) < -150.0_SiKi) .OR. &
                 (MINVAL(MnDriftData%Data4D%WvDir2) > 150.0_SiKi) .OR.  (MAXVAL(MnDriftData%Data4D%WvDir2) < -150.0_SiKi) ) THEN
               CALL SetErrStat( ErrID_Warn,' There may be issues with how the wave direction data is handled when the wave '// &
                                          'direction of interest is near the +/- 180 direction.  This is a known issue with '// &
                                          'the WAMIT2 module that has not yet been addressed.',ErrStat,ErrMsg,RoutineName)
            ENDIF

               ! Now check the limits for the first wave direction
               !  --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
               !  --> FIXME: modify to allow shifting values by TwoPi before comparing
            IF ( InitInp%WaveDirMin < MINVAL(MnDriftData%Data4D%WvDir1) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                     'found in the WAMIT data file '//TRIM(MnDriftData%Filename)//' for the first wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF
            IF ( InitInp%WaveDirMax > MAXVAL(MnDriftData%Data4D%WvDir1) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                     'found in the WAMIT data file '//TRIM(MnDriftData%Filename)//' for the first wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF


               ! Now check the limits for the second wave direction
               !   --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
            IF ( InitInp%WaveDirMin < MINVAL(MnDriftData%Data4D%WvDir2) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                     'found in the WAMIT data file '//TRIM(MnDriftData%Filename)//' for the second wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF
            IF ( InitInp%WaveDirMax > MAXVAL(MnDriftData%Data4D%WvDir2) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                     'found in the WAMIT data file '//TRIM(MnDriftData%Filename)//' for the second wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF

         ENDIF

      ELSE
            ! No data. This is a catastrophic issue.  We should not have called this routine without data that is usable for the MnDrift calculation
         CALL SetErrStat( ErrID_Fatal, ' Mean drift calculation called without data.',ErrStat,ErrMsg,RoutineName)
      ENDIF

      IF ( ErrStat >= AbortErrLev )    RETURN



         !> 4. Check the data to see if we need to convert to 3D arrays before continuing (4D is sparse in any dimension we want and
         !!    frequency diagonal is complete).  Only check if we don't have 3D data.

      IF ( .NOT. MnDriftData%DataIs3D .AND. MnDriftData%Data4D%WvFreqDiagComplete ) THEN
         TmpFlag = .FALSE.    ! if this goes true, then we need to convert to 3D data
         DO IBody=1,MnDriftData%Data4D%NumBodies
            IF (TmpFlag) EXIT
            DO ThisDim=1,6
               Idx = (IBody-1)*6+ThisDim
               IF ( p%MnDriftDims(IBody) ) THEN        ! Flag indicating which dimension we are calculating for
                  IF ( MnDriftData%Data4D%DataIsSparse(Idx) .AND. MnDriftData%Data4D%LoadComponents(Idx) ) THEN
                     TmpFlag = .TRUE.
                     EXIT ! inner DO
                  END IF
               ENDIF
            ENDDO
         ENDDO

            ! If we need to create the 3D data set, then
         IF (TmpFlag) THEN
            CALL Copy_InitData4Dto3D( MnDriftData%Data4D, MnDriftData%Data3D, ErrStatTmp, ErrMsgTmp )
            CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
            IF (ErrStat >= AbortErrLev) RETURN

            MnDriftData%DataIs3D = .TRUE.       ! Set flag to indicate we now have the 3D data.
         END IF ! TmpFlag

      ENDIF



         !> 5. Now check to make sure we have data that will work.  For either 3D or 4D data, it must not be sparse.
         !!    We simplified the 4D sparse case to 3D above if we could.  To check this, we have to check the load
         !!    components that we will use.  So, we will loop through them and set the TmpFlag to true if there is
         !!    a sparse matrix for one of them.
         !! FIXME: remove this check and warning once the sparse matrix interpolation routines are implimented.
      TmpFlag = .FALSE.
      IF ( MnDriftData%DataIs3D ) THEN
         DO IBody=1,MnDriftData%Data3D%NumBodies
            IF (TmpFlag) EXIT
            DO ThisDim=1,6
               Idx = (IBody-1)*6+ThisDim
               IF ( MnDriftData%Data3D%DataIsSparse(Idx) .AND. MnDriftData%Data3D%LoadComponents(Idx) .AND. p%MnDriftDims(ThisDim) ) THEN
                  TmpFlag = .TRUE.
                  EXIT
               END IF
            ENDDO
         ENDDO
      ELSE     ! must be 4D -- we checked that we had something at the start of this routine.
         DO IBody=1,MnDriftData%Data4D%NumBodies
            IF (TmpFlag) EXIT
            DO ThisDim=1,6
               Idx = (IBody-1)*6+ThisDim
               IF ( MnDriftData%Data4D%DataIsSparse(Idx) .AND. MnDriftData%Data4D%LoadComponents(Idx) .AND. p%MnDriftDims(ThisDim) ) THEN
                  TmpFlag = .TRUE.
                  EXIT
               END IF
            ENDDO
         ENDDO
      ENDIF
      IF (TmpFlag) THEN
         CALL SetErrStat(ErrID_Fatal,' The second order WAMIT data in '//TRIM(MnDriftData%Filename)//' is too sparse '// &
               'for the interpolation routine used in the mean drift calculation.  At some later point, we will allow for '// &
               'sparse data when a different interpolation routine is implimented.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF


         !> 5.  Calculate the mean drift force
         !! The single summation equation used here is given by
         !! \f$  {F_{{ex}~k}^{{-}(2)}} = \Re \left( \sum\limits_{m=0}^{N/2}
         !!                                        {A_m} {A_m^*} F_k^{-}(\omega_m, \beta_m)\right) \cdot\Delta\omega \f$
         !!      for \f$\quad k=1,2,\ldots,6,      \f$


         ! Setup temporary arrays that we will be passing the data to the interpolation routines in.
      IF (MnDriftData%DataIs3D) THEN
         ALLOCATE( TmpData3D(MnDriftData%Data3D%NumWvFreq1, MnDriftData%Data3D%NumWvDir1,MnDriftData%Data3D%NumWvDir2),STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate temporary array for interpolation of 3D QTF data.', &
                                                ErrStat, ErrMsg, RoutineName)
      ELSE
         ALLOCATE( TmpData4D(MnDriftData%Data4D%NumWvFreq1,MnDriftData%Data4D%NumWvFreq2,MnDriftData%Data4D%NumWvDir1,MnDriftData%Data4D%NumWvDir2),STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate temporary array for interpolation of 4D QTF data.', &
                                                ErrStat, ErrMsg, RoutineName)
      ENDIF

         ! If something went wrong during allocation of the temporary arrays...
      IF ( ErrStat >= AbortErrLev ) THEN
         call cleanup()
         RETURN
      ENDIF


         ! Now loop through all the dimensions and perform the calculation
      DO IBody=1,p%NBody

            ! Heading correction, only applies to NBodyMod == 2
         if (p%NBodyMod==2) then
            RotateZdegOffset = InitInp%PtfmRefztRot(IBody)*R2D
         else
            RotateZdegOffset = 0.0_SiKi
         endif

            ! NOTE: RotateZMatrixT is the rotation from local to global.
         RotateZMatrixT(:,1) = (/  cos(InitInp%PtfmRefztRot(IBody)), -sin(InitInp%PtfmRefztRot(IBody)) /)
         RotateZMatrixT(:,2) = (/  sin(InitInp%PtfmRefztRot(IBody)),  cos(InitInp%PtfmRefztRot(IBody)) /)


         DO ThisDim=1,6

            Idx = (IBody-1)*6 + ThisDim

               ! Set the MnDrift force to 0.0 (Even ones we don't calculate)
            MnDriftForce(Idx)   = 0.0_SiKi

            IF (MnDriftData%DataIs3D) THEN
               TmpFlag = MnDriftData%Data3D%LoadComponents(Idx)
            ELSE
               TmpFlag = MnDriftData%Data4D%LoadComponents(Idx)
            END IF

               ! Only on the dimensions we requested, and if it is present in the data
            IF ( p%MnDriftDims(ThisDim) .AND. TmpFlag ) THEN

                  ! Set an initial search index for the 3D and 4D array interpolation
               LastIndex3 = (/0,0,0/)
               LastIndex4 = (/0,0,0,0/)

                  ! To make things run slightly quicker, copy the data we will be interpolating over into the temporary arrays
               IF (MnDriftData%DataIs3D) THEN
                  TmpData3D = MnDriftData%Data3D%DataSet(:,:,:,Idx)
               ELSE
                  TmpData4D = MnDriftData%Data4D%DataSet(:,:,:,:,Idx)
               END IF


               DO J=1,InitInp%NStepWave2

                     ! NOTE: since the Mean Drift only returns a static time independent average value for the drift force, we do not
                     !        need to account for any offset in the location of the WAMIT body (this term vanishes).
                     ! First get the wave amplitude -- must be reconstructed from the WaveElevC0 array.  First index is the real (1) or
                     ! imaginary (2) part.  Divide by NStepWave2 to remove the built in normalization in WaveElevC0.
                  aWaveElevC = CMPLX( InitInp%WaveElevC0(1,J), InitInp%WaveElevC0(2,J), SiKi) / InitInp%NStepWave2

                     ! Calculate the frequency
                  Omega1 = J * InitInp%WaveDOmega


                     ! Only get a QTF value if within the range of frequencies we have wave amplitudes for (first order cutoffs).  This
                     ! is done only for efficiency. 
                  
                  !BJJ: If WaveMod==1, this could result in zeroing out the wrong values... 
                  !InitInp%WvLowCOff and InitInp%WvHiCOff are not used in SeaState when WaveMod = 0,1, or 6
                  ! Probably could just remove this IF statement????
                  IF ( (Omega1 >= InitInp%WvLowCOff) .AND. (Omega1 <= InitInp%WvHiCOff) ) THEN

                        ! Now get the QTF value that corresponds to this frequency and wavedirection pair.
                     IF ( MnDriftData%DataIs3D ) THEN

                           ! Set the (omega1,beta1,beta2) point we are looking for. (angles in degrees here)
                        Coord3 = (/ REAL(Omega1,SiKi), InitInp%WaveDirArr(J), InitInp%WaveDirArr(J) /)

                           ! Apply local Z rotation to heading angle (degrees) to put wave direction into the local (rotated) body frame
                        Coord3(2) = Coord3(2) - RotateZdegOffset
                        Coord3(3) = Coord3(3) - RotateZdegOffset

                           ! get the interpolated value for F(omega1,beta1,beta2)
                        CALL WAMIT_Interp3D_Cplx( Coord3, TmpData3D, MnDriftData%Data3D%WvFreq1, &
                                             MnDriftData%Data3D%WvDir1, MnDriftData%Data3D%WvDir2, LastIndex3, QTF_Value, ErrStatTmp, ErrMsgTmp )
                           CALL SetErrStat(ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)

                     ELSE

                           ! Set the (omega1,omega2,beta1,beta2) point we are looking for. (angles in degrees here)
                        Coord4 = (/ REAL(Omega1,SiKi), REAL(Omega1,SiKi), InitInp%WaveDirArr(J), InitInp%WaveDirArr(J) /)

                           ! Apply local Z rotation to heading angle (degrees) to put wave direction into the local (rotated) body frame
                        Coord4(3) = Coord4(3) - RotateZdegOffset
                        Coord4(4) = Coord4(4) - RotateZdegOffset

                           ! get the interpolated value for F(omega1,omega2,beta1,beta2)
                        CALL WAMIT_Interp4D_Cplx( Coord4, TmpData4D, MnDriftData%Data4D%WvFreq1, MnDriftData%Data4D%WvFreq2, &
                                             MnDriftData%Data4D%WvDir1, MnDriftData%Data4D%WvDir2, LastIndex4, QTF_Value, ErrStatTmp, ErrMsgTmp )
                           CALL SetErrStat(ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)


                     ENDIF !QTF value find


                  ELSE     ! outside the frequency range

                     QTF_Value = CMPLX(0.0,0.0,SiKi)

                  ENDIF    ! frequency check


                     ! Check and make sure nothing bombed in the interpolation that we need to be aware of
                  IF ( ErrStat >= AbortErrLev ) THEN
                     call cleanup()
                     RETURN
                  ENDIF


                     ! Now we have the value of the QTF.  These values should only be real for the omega1=omega2 case of the mean drift.
                     ! However if the value came from the 4D interpolation routine, it might have some residual complex part to it.  So
                     ! we throw the complex part out.
                  QTF_Value = CMPLX(REAL(QTF_Value,SiKi),0.0,SiKi)

                     ! NOTE:  any offset in platform location vanishes when the only the REAL part is kept (the offset resides in the
                     !        phase shift, which is in the imaginary part)
                     ! Now put it all together... note the frequency stepsize is multiplied after the summation
                  MnDriftForce(Idx) = MnDriftForce(Idx) + REAL(QTF_Value * aWaveElevC * CONJG(aWaveElevC)) !bjj: put QTF_Value first so that if it's zero, the rest gets set to zero (to hopefully avoid overflow issues)

               ENDDO ! NStepWave2

            ENDIF    ! Load component to calculate


         ENDDO ! ThisDim   -- Load Component on body


            ! Now rotate the force components with platform orientation
         MnDriftForce(1:2) = MATMUL( RotateZMatrixT, MnDriftForce(1:2) )       ! Fx and Fy, rotation about z
         MnDriftForce(4:5) = MATMUL( RotateZMatrixT, MnDriftForce(4:5) )       ! Mx and My, rotation about z

      ENDDO    ! IBody



         ! Cleanup
      call cleanup()
      
   CONTAINS
      subroutine cleanup()
         IF (ALLOCATED(TmpData3D))        DEALLOCATE(TmpData3D,STAT=ErrStatTmp)
         IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
      end subroutine cleanup

   END SUBROUTINE MnDrift_InitCalc







   !-------------------------------------------------------------------------------------------------------------------------------
   !> This subroutine calculates the force time series using the NewmanApp calculation.
   !! The data is stored in either 3D or 4D arrays depending on the file type used.
   !! At each step in the summation of the mth term, a call is made to the 3D or 4D interpolation algorithm to find the value of
   !! \f$ F^-_k(\omega_m, \omega_n) \f$ corresponding to the \f$ Z(\omega_m) \f$ term in the complex wave amplitude, _aWaveElevC_.
   !! A combination of the limits of \f$ \max(\omega_{lo-d},\omega_{lo}) \le \omega_m \le  \min(\omega_{hi-d},\omega_{hi}) \f$ are
   !! imposed during the summation before the FFT and values outside this range are set to zero.
   !!
   !! For multi-directional waves where the equal energy discretization is used, each frequency has a single wave direction
   !! associated with it. Since the Newman's approximation calculation only involves summing over terms involving only a single
   !! frequency at a time, only a single wave direction is involved at each step. So, if all the diagonal elements of the 4D matrix
   !! where \f$ \omega_1 = \omega_2 \f$ and \f$ \beta_1 = \beta_2 \f$ were present, it would be possible to simplify the 4D
   !! interpolation required to two dimensional interpolation for this particular case.  However, since the same interpolation
   !! routine is used for all the 4D data handling, it is programatically simpler to use only one interpolation algorithm for the
   !! 4D data and incur a slight penalty at the initialization step where this routine is called.  Overall it is a small price to
   !! pay for maintainable code.
   !!
   !! The single summation equation used here is Standing's version of Newman's equation (written in a slightly different form than
   !! appears in Standing's paper, but mathematically equivalent).  Standing's version differs sligthly from Newman's original
   !! equation in that it allows for multidirectional waves, and does not require high frequency filtering afterwards.  Our
   !! mathematically equivalent equation to Standing's version of Newman's equation is given as,
   !!    \f$      {F_{{ex}~k}^{- (2)}} \approx
   !!                      \left.   \left[ abs \left( \sum\limits_{m=0}^{N/2}
   !!                                     {a_m} \sqrt{   F_k^{-}(\omega_m,\omega_m)}
   !!                                 \cdot \mathrm{e}^{i \omega_m t} \right) \right]^2 \right|_{F_{k}^{-}(\omega_m,\omega_m) > 0}
   !!                   -  \left.   \left[ abs \left( \sum\limits_{m=0}^{N/2}
   !!                                  {a_m} \sqrt{  - F_k^{-}(\omega_m,\omega_m)}
   !!                                  \cdot \mathrm{e}^{i \omega_m t} \right) \right]^2 \right|_{F_{k}^{-}(\omega_m,\omega_m) < 0} \quad \f$
   !!                      for     \f$ \quad k=1,2,\ldots,6,         \f$
   !!
   !! where \f$ k \f$ indicates the index to the load component,  \f$ {F_{{ex}~k}^{{-}(2)}} \f$ is the resulting Newman's
   !! approximation of the second order forces, and \f$ a_m \f$ is the complex wave amplitude for the \f$ m^{th} \f$ frequency.  Note
   !! that the two IFFTs cover different values of \f$ F_k^{-} \f$.   Note also that \f$ F_k^{-} (\omega_m, \beta_m) \f$ is the
   !! dimensionalized real valued diagonal of the QTF read from the WAMIT file and interpolated for the \f$ m^{th} \f$ wave frequency.
   !! Note that this is an IFFT and not a numerical integral.
   !!
   !! For the Newman approximation, the limits used on the calculations are a little different than with the other methods.  For
   !! this method, the limits are imposed as the maximum of the low frequency limits, and the minimum of the high frequency limits.
   !!
   !! @note The inverse Fourier transform in each term is over data that is not conjugate symmetric (where the negative frequencies
   !!       are the complex conjugate of the positive frequencies).  This means we cannot use the IFFTs in FFT_Mod.  To get around
   !!       this, the SlowFT_cx2real routine was written.  It evaluates the sum at each timestep.
   !!
   !! Since the frequency range of the QTF has not yet been checked, we will do that now.
   !!
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE NewmanApp_InitCalc( InitInp, p, NewmanAppData, NewmanAppForce, ErrMsg, ErrStat )

      IMPLICIT NONE

      TYPE(WAMIT2_InitInputType),         INTENT(IN   )  :: InitInp              !< Input data for initialization routine
      TYPE(WAMIT2_ParameterType),         INTENT(IN   )  :: p                    !< Parameters
      TYPE(W2_DiffData_Type),             INTENT(INOUT)  :: NewmanAppData        !< Data storage for the NewmanApp method.  Set to INOUT in case we need to convert 4D to 3D
      REAL(SiKi),  ALLOCATABLE,           INTENT(  OUT)  :: NewmanAppForce(:,:)  !< Force data.  Index 1 is the timestep, index 2 is the load component.
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat

         ! Local Variables
      CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary error message for calls
      INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
      REAL(SiKi)                                         :: TmpReal1             !< Temporary real
!     REAL(SiKi)                                         :: TmpReal2             !< Temporary real
      LOGICAL                                            :: TmpFlag              !< Temporary logical flag
      INTEGER(IntKi)                                     :: ThisDim              !< Generic counter for dimension
      INTEGER(IntKi)                                     :: IBody                !< Index to which body we are on
      INTEGER(IntKi)                                     :: Idx                  !< Index to the full set of 6*NBody
      INTEGER(IntKi)                                     :: J                    !< Generic counter
!      INTEGER(IntKi)                                     :: K                    !< Generic counter
      TYPE(FFT_DataType)                                 :: FFT_Data             !< Temporary array for the FFT module we're using
      CHARACTER(*), PARAMETER                            :: RoutineName = 'NewmanApp_InitCalc'


         ! Wave information and QTF temporary
      COMPLEX(SiKi)                                      :: QTF_Value            !< Temporary complex number for QTF
      COMPLEX(SiKi), ALLOCATABLE                         :: NewmanTerm1C(:,:)    !< First  term in the newman calculation, complex frequency space.  All dimensions, this body.
      COMPLEX(SiKi), ALLOCATABLE                         :: NewmanTerm2C(:,:)    !< Second term in the newman calculation, complex frequency space.  All dimensions, this body.
      COMPLEX(SiKi), ALLOCATABLE                         :: NewmanTerm1t(:)      !< First  term in the newman calculation, time domain.  Current load dimension.
      COMPLEX(SiKi), ALLOCATABLE                         :: NewmanTerm2t(:)      !< Second term in the newman calculation, time domain.  Current load dimension.
      COMPLEX(SiKi)                                      :: aWaveElevC           !< Wave elevation of current frequency component.  NStepWave2 factor removed.
      REAL(ReKi)                                         :: Omega1               !< Wave frequency of this component

         ! Interpolation routine indices and value to search for, and smaller array to pass
      INTEGER(IntKi)                                     :: LastIndex3(3)        !< Last used index for searching in the interpolation algorithms
      INTEGER(IntKi)                                     :: LastIndex4(4)        !< Last used index for searching in the interpolation algorithms
      REAL(SiKi)                                         :: Coord3(3)            !< The (omega1,beta1,beta2) coordinate we want in the 3D dataset
      REAL(SiKi)                                         :: Coord4(4)            !< The (omega1,omega2,beta1,beta2) coordinate we want in the 4D dataset
      REAL(SiKi)                                         :: RotateZdegOffset     !< Offset to wave heading (NBodyMod==2 only)
      REAL(SiKi)                                         :: RotateZMatrixT(2,2)  !< The transpose of rotation in matrix form for rotation about z (from global to local)
      COMPLEX(SiKi)                                      :: PhaseShiftXY         !< The phase shift offset to apply to the body
      REAL(SiKi)                                         :: WaveNmbr1            !< Wavenumber for this frequency
      COMPLEX(SiKi), ALLOCATABLE                         :: TmpData3D(:,:,:)     !< Temporary 3D array we put the 3D data into (minus the load component indice)
      COMPLEX(SiKi), ALLOCATABLE                         :: TmpData4D(:,:,:,:)   !< Temporary 4D array we put the 4D data into (minus the load component indice)


         ! Initialize a few things
      ErrMsg      = ''
      ErrMsgTmp   = ''
      ErrStat     = ErrID_None
      ErrStatTmp  = ErrID_None



         !> 1. Check the data to see if the wave frequencies are present in the QTF data.  Since Newman's approximation only uses
         !!    frequencies where \f$ \omega_1=\omega_2 \f$, the data read in from the files must contain the full range of frequencies
         !!    present in the waves.
!bjj: InitInp%WvLowCOff and InitInp%WvHiCOff aren't supposed to be used when WaveMod=0, 1, or 6, but they are used here regardless of those conditions.
!     Can we get rid of these checks????
      IF ( NewmanAppData%DataIs3D ) THEN

            ! Check the low frequency cutoff
         IF ( MINVAL( NewmanAppData%Data3D%WvFreq1 ) > InitInp%WvLowCOff ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(NewmanAppData%Data3D%WvFreq1)))// &
                           ' rad/s for first wave period) data in '//TRIM(NewmanAppData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOff.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF

            ! Check the high frequency cutoff -- using the Difference high frequency cutoff.  The first order high frequency
            ! cutoff is typically too high for this in most cases.
         IF ( MAXVAL(NewmanAppData%Data3D%WvFreq1 ) < InitInp%WvHiCOff ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(NewmanAppData%Data3D%WvFreq1)))// &
                           ' rad/s for first wave period) data in '//TRIM(NewmanAppData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOff.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF

      ELSE IF ( NewmanAppData%DataIs4D ) THEN   ! only check if not 3D data. If there is 3D data, we default to using it for calculations

             ! Check the low frequency cutoff
         IF ( MINVAL( NewmanAppData%Data4D%WvFreq1 ) > InitInp%WvLowCOff ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(NewmanAppData%Data4D%WvFreq1)))// &
                           ' rad/s first wave period) data in '//TRIM(NewmanAppData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOff.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF
         IF ( MINVAL( NewmanAppData%Data4D%WvFreq2 ) > InitInp%WvLowCOff ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(NewmanAppData%Data4D%WvFreq2)))// &
                           ' rad/s for second wave period) data in '//TRIM(NewmanAppData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOff.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF

            ! Check the high frequency cutoff -- using the Difference high frequency cutoff.  The first order high frequency
            ! cutoff is typically too high for this in most cases.
         IF ( MAXVAL(NewmanAppData%Data4D%WvFreq1) < InitInp%WvHiCOff ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(NewmanAppData%Data4D%WvFreq1)))// &
                           ' rad/s for first wave period) data in '//TRIM(NewmanAppData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOff.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF
         IF ( MAXVAL(NewmanAppData%Data4D%WvFreq2) < InitInp%WvHiCOff ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(NewmanAppData%Data4D%WvFreq1)))// &
                           ' rad/s second wave period) data in '//TRIM(NewmanAppData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOff.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF

      ELSE
            ! This is a catastrophic issue.  We should not have called this routine without data that is usable for the NewmanApp calculation
         CALL SetErrStat( ErrID_Fatal, ' Newman approximation calculation called without data.',ErrStat,ErrMsg,RoutineName)
      ENDIF

      IF ( ErrStat >= AbortErrLev )    RETURN



         !> 2. Check the data to see if the wave directions are present.  May need to adjust for the boundary at +/- PI
      IF ( NewmanAppData%DataIs3D ) THEN

            ! If we are using multidirectional waves, then we should have more than 1 wave direction in the WAMIT file.
         IF ( InitInp%WaveMultiDir .AND. (NewmanAppData%Data3D%NumWvDir1 == 1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(NewmanAppData%Filename)//' only contains one wave '// &
                        'direction at '//TRIM(Num2LStr(NewmanAppData%Data3D%WvDir1(1)))//' degrees (first wave direction). '// &
                        'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                        ErrStat,ErrMsg,RoutineName)
         ELSE IF ( InitInp%WaveMultiDir .AND. (NewmanAppData%Data3D%NumWvDir2 == 1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(NewmanAppData%Filename)//' only contains one wave '// &
                        'direction at '//TRIM(Num2LStr(NewmanAppData%Data3D%WvDir2(1)))//' degrees (second wave direction). '// &
                        'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                        ErrStat,ErrMsg,RoutineName)
         ELSE

               ! See Known Issues #1 at the top of this file.  There may be problems if the data spans the +/- Pi boundary.  For
               ! now (since time is limited) we will issue a warning if any of the wave directions for multidirectional waves
               ! or data from the WAMIT file for the wavedirections is close to the +/-pi boundary (>150 degrees, <-150 degrees),
               ! we will issue a warning.
            IF ( (InitInp%WaveDirMin > 150.0_SiKi) .OR. (InitInp%WaveDirMax < -150.0_SiKi) .OR. &
                 (minval(NewmanAppData%data3d%WvDir1) > 150.0_SiKi) .OR.  (maxval(NewmanAppData%data3d%WvDir1) < -150.0_SiKi) .OR. &
                 (minval(NewmanAppData%data3d%WvDir2) > 150.0_SiKi) .OR.  (maxval(NewmanAppData%data3d%WvDir2) < -150.0_SiKi) ) THEN
               CALL SetErrStat( ErrID_Warn,' There may be issues with how the wave direction data is handled when the wave '// &
                                          'direction of interest is near the +/- 180 direction.  This is a known issue with '// &
                                          'the WAMIT2 module that has not yet been addressed.',ErrStat,ErrMsg,RoutineName)
            ENDIF

               ! Now check the limits for the first wave direction
               !   --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
            IF ( InitInp%WaveDirMin < MINVAL(NewmanAppData%Data3D%WvDir1) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                     'found in the WAMIT data file '//TRIM(NewmanAppData%Filename)//' for the first wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF
            IF ( InitInp%WaveDirMax > MAXVAL(NewmanAppData%Data3D%WvDir1) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                     'found in the WAMIT data file '//TRIM(NewmanAppData%Filename)//' for the first wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
           ENDIF


               ! Now check the limits for the second wave direction
               !   --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
            IF ( InitInp%WaveDirMin < MINVAL(NewmanAppData%Data3D%WvDir2) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                     'found in the WAMIT data file '//TRIM(NewmanAppData%Filename)//' for the second wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF
            IF ( InitInp%WaveDirMax > MAXVAL(NewmanAppData%Data3D%WvDir2) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                     'found in the WAMIT data file '//TRIM(NewmanAppData%Filename)//' for the second wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF

         ENDIF

      ELSEIF ( NewmanAppData%DataIs4D ) THEN

            ! If we are using multidirectional waves, then we should have more than 1 wave direction in the WAMIT file.
         IF ( InitInp%WaveMultiDir .AND. (NewmanAppData%Data4D%NumWvDir1 == 1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(NewmanAppData%Filename)//' only contains one wave '// &
                        'direction at '//TRIM(Num2LStr(NewmanAppData%Data4D%WvDir1(1)))//' degrees (first wave direction). '// &
                        'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                        ErrStat,ErrMsg,RoutineName)
         ELSE IF ( InitInp%WaveMultiDir .AND. (NewmanAppData%Data4D%NumWvDir2 == 1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(NewmanAppData%Filename)//' only contains one wave '// &
                        'direction at '//TRIM(Num2LStr(NewmanAppData%Data4D%WvDir2(1)))//' degrees (second wave direction). '// &
                        'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                        ErrStat,ErrMsg,RoutineName)
         ELSE

               ! See Known Issues #1 at the top of this file.  There may be problems if the data spans the +/- Pi boundary.  For
               ! now (since time is limited) we will issue a warning if any of the wave directions for multidirectional waves
               ! or data from the WAMIT file for the wavedirections is close to the +/-pi boundary (>150 degrees, <-150 degrees),
               ! we will issue a warning.
            IF ( (InitInp%WaveDirMin > 150.0_SiKi) .OR. (InitInp%WaveDirMax < -150.0_SiKi) .OR. &
                 (MINVAL(NewmanAppData%Data4D%WvDir1) > 150.0_SiKi) .OR.  (MAXVAL(NewmanAppData%Data4D%WvDir1) < -150.0_SiKi) .OR. &
                 (MINVAL(NewmanAppData%Data4D%WvDir2) > 150.0_SiKi) .OR.  (MAXVAL(NewmanAppData%Data4D%WvDir2) < -150.0_SiKi) ) THEN
               CALL SetErrStat( ErrID_Warn,' There may be issues with how the wave direction data is handled when the wave '// &
                                          'direction of interest is near the +/- 180 direction.  This is a known issue with '// &
                                          'the WAMIT2 module that has not yet been addressed.',ErrStat,ErrMsg,RoutineName)
            ENDIF

               ! Now check the limits for the first wave direction
               !  --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
               !  --> FIXME: modify to allow shifting values by TwoPi before comparing
            IF ( InitInp%WaveDirMin < MINVAL(NewmanAppData%Data4D%WvDir1) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                     'found in the WAMIT data file '//TRIM(NewmanAppData%Filename)//' for the first wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF
            IF ( InitInp%WaveDirMax > MAXVAL(NewmanAppData%Data4D%WvDir1) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                     'found in the WAMIT data file '//TRIM(NewmanAppData%Filename)//' for the first wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF


               ! Now check the limits for the second wave direction
               !   --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
            IF ( InitInp%WaveDirMin < MINVAL(NewmanAppData%Data4D%WvDir2) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                     'found in the WAMIT data file '//TRIM(NewmanAppData%Filename)//' for the second wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF
            IF ( InitInp%WaveDirMax > MAXVAL(NewmanAppData%Data4D%WvDir2) ) THEN
               CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                     'found in the WAMIT data file '//TRIM(NewmanAppData%Filename)//' for the second wave direction.', &
                     ErrStat, ErrMsg, RoutineName)
            ENDIF

         ENDIF

      ELSE
            ! No data. This is a catastrophic issue.  We should not have called this routine without data that is usable for the NewmanApp calculation
         CALL SetErrStat( ErrID_Fatal, ' Newman approximation calculation called without data.',ErrStat,ErrMsg,RoutineName)
      ENDIF

      IF ( ErrStat >= AbortErrLev )    RETURN



         !> 3. Check the data to see if we need to convert to 3D arrays before continuing (4D is sparse in any dimension we want and
         !!    frequency diagonal is complete).  Only check if we don't have 3D data.

      IF ( .NOT. NewmanAppData%DataIs3D .AND. NewmanAppData%Data4D%WvFreqDiagComplete ) THEN
         TmpFlag = .FALSE.    ! if this goes true, then we need to convert to 3D data
         DO IBody=1,NewmanAppData%Data4D%NumBodies
            DO ThisDim=1,6
               Idx = (IBody-1)*6+ThisDim
               IF ( p%NewmanAppDims(ThisDim) ) THEN        ! Flag indicating which dimension we are calculating for
                  IF ( NewmanAppData%Data4D%DataIsSparse(Idx) .AND. NewmanAppData%Data4D%LoadComponents(Idx) )      TmpFlag = .TRUE.
               ENDIF
            ENDDO
         ENDDO

            ! If we need to create the 3D data set, then
         CALL Copy_InitData4Dto3D( NewmanAppData%Data4D, NewmanAppData%Data3D, ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) RETURN

         NewmanAppData%DataIs3D = .TRUE.       ! Set flag to indicate we now have the 3D data.

      ENDIF



         !> 4. Now check to make sure we have data that will work.  For either 3D or 4D data, it must not be sparse.
         !!    We simplified the 4D sparse case to 3D above if we could.  To check this, we have to check the load
         !!    components that we will use.  So, we will loop through them and set the TmpFlag to true if there is
         !!    a sparse matrix for one of them.
         !! FIXME: remove this check and warning once the sparse matrix interpolation routines are implimented.
      TmpFlag = .FALSE.
      IF ( NewmanAppData%DataIs3D ) THEN
         DO IBody=1,NewmanAppData%Data3D%NumBodies
            DO ThisDim=1,6
               Idx = (IBody-1)*6+ThisDim
               IF ( NewmanAppData%Data3D%DataIsSparse(Idx) .AND. NewmanAppData%Data3D%LoadComponents(Idx) .AND. p%NewmanAppDims(ThisDim) )    TmpFlag = .TRUE.
            ENDDO
         ENDDO
      ELSE     ! must be 4D -- we checked that we had something at the start of this routine.
         DO IBody=1,NewmanAppData%Data4D%NumBodies
            DO ThisDim=1,6
               Idx = (IBody-1)*6+ThisDim
               IF ( NewmanAppData%Data4D%DataIsSparse(Idx) .AND. NewmanAppData%Data4D%LoadComponents(Idx) .AND. p%NewmanAppDims(ThisDim) )    TmpFlag = .TRUE.
            ENDDO
         ENDDO
      ENDIF
      IF (TmpFlag) THEN
         CALL SetErrStat(ErrID_Fatal,' The second order WAMIT data in '//TRIM(NewmanAppData%Filename)//' is too sparse '// &
               'for the interpolation routine used in the Newman approximation calculation.  At some later point, we will allow for '// &
               'sparse data when a different interpolation routine is implimented.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF


         !> 5.  Calculate the Newman's approximation for the second order forces
         !! The single summation equation used here is given by
         !!    \f$      {F_{{ex}~k}^{- (2)}} \approx
         !!                      \left.   \left| \sum\limits_{m=0}^{N/2}
         !!                                     {a_m} \sqrt{   F_k^{-}(\omega_m,\omega_m)}
         !!                                 \cdot \mathrm{e}^{i \omega_m t} \right|^2 \right|_{F_{k}^{-}(\omega_m,\omega_m) > 0}
         !!                   -  \left.   \left| \sum\limits_{m=0}^{N/2}
         !!                                  {a_m} \sqrt{  - F_k^{-}(\omega_m,\omega_m)}
         !!                                  \cdot \mathrm{e}^{i \omega_m t} \right|^2 \right|_{F_{k}^{-}(\omega_m,\omega_m) < 0} \quad \f$
         !!                      for     \f$ \quad k=1,2,\ldots,6,         \f$


         ! Setup temporary arrays that we will be passing the data to the interpolation routines in.
      IF (NewmanAppData%DataIs3D) THEN
         ALLOCATE( TmpData3D(NewmanAppData%Data3D%NumWvFreq1, NewmanAppData%Data3D%NumWvDir1,NewmanAppData%Data3D%NumWvDir2),STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate temporary array for interpolation of 3D QTF data.', &
                                                ErrStat, ErrMsg, RoutineName)
      ELSE
         ALLOCATE( TmpData4D(NewmanAppData%Data4D%NumWvFreq1,NewmanAppData%Data4D%NumWvFreq1,NewmanAppData%Data4D%NumWvDir1,NewmanAppData%Data4D%NumWvDir2),STAT=ErrStatTmp )
         IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate temporary array for interpolation of 4D QTF data.', &
                                                ErrStat, ErrMsg, RoutineName)
      ENDIF


         ! Setup the arrays holding the Newman terms, both the complex frequency domain and real time domain pieces
      ALLOCATE( NewmanTerm1t( 0:InitInp%NStepWave  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for calculating the first term of the Newmans '// &
                                             'approximation in the time domain.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( NewmanTerm2t( 0:InitInp%NStepWave  ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for calculating the second term of the Newmans '// &
                                             'approximation in the time domain.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( NewmanTerm1C( 0:InitInp%NStepWave2, 6 ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for calculating the first term of the Newmans '// &
                                             'approximation in the frequency domain.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( NewmanTerm2C( 0:InitInp%NStepWave2, 6 ), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for calculating the second term of the Newmans '// &
                                             'approximation in the frequency domain.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( NewmanAppForce( 0:InitInp%NStepWave, 6*p%NBody), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for the resulting Newmans '// &
                                             'approximation of the 2nd order force.',ErrStat, ErrMsg, RoutineName)



         ! If something went wrong during allocation of the temporary arrays...
      IF ( ErrStat >= AbortErrLev ) THEN
         IF (ALLOCATED(TmpData3D))        DEALLOCATE(TmpData3D,STAT=ErrStatTmp)
         IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm1t))     DEALLOCATE(NewmanTerm1t,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm2t))     DEALLOCATE(NewmanTerm2t,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm1C))     DEALLOCATE(NewmanTerm1C,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm2C))     DEALLOCATE(NewmanTerm2C,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanAppForce))   DEALLOCATE(NewmanAppForce,STAT=ErrStatTmp)
         RETURN
      ENDIF


         ! Set the resulting force to zero.  Will calculate for the dimensions requested.
      NewmanAppForce = 0.0_SiKi


         ! Initialize the FFT library
      CALL InitCFFT ( InitInp%NStepWave, FFT_Data, .FALSE., ErrStatTmp )      ! Complex result FFT initialize
      CALL SetErrStat(ErrStatTmp,'Error occured while initializing the FFT.',ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         IF (ALLOCATED(TmpData3D))        DEALLOCATE(TmpData3D,STAT=ErrStatTmp)
         IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm1t))     DEALLOCATE(NewmanTerm1t,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm2t))     DEALLOCATE(NewmanTerm2t,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm1C))     DEALLOCATE(NewmanTerm1C,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm2C))     DEALLOCATE(NewmanTerm2C,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanAppForce))   DEALLOCATE(NewmanAppForce,STAT=ErrStatTmp)
         RETURN
      END IF


         ! Loop through all bodies
      DO IBody=1,p%NBody

            ! set all frequency terms to zero to start
         NewmanTerm1C(:,:) = CMPLX(0.0, 0.0, SiKi)
         NewmanTerm2C(:,:) = CMPLX(0.0, 0.0, SiKi)


            ! Heading correction, only applies to NBodyMod == 2
         if (p%NBodyMod==2) then
            RotateZdegOffset = InitInp%PtfmRefztRot(IBody)*R2D
         else
            RotateZdegOffset = 0.0_SiKi
         endif

         !----------------------------------------------------
         ! Populate the frequency terms for this body
         !----------------------------------------------------

         DO ThisDim=1,6

            Idx= (IBody-1)*6+ThisDim

            IF (NewmanAppData%DataIs3D) THEN
               TmpFlag = NewmanAppData%Data3D%LoadComponents(Idx)
            ELSE
               TmpFlag = NewmanAppData%Data4D%LoadComponents(Idx)
            END IF

               ! Only on the dimensions we requested, and if it is present in the data
            IF ( p%NewmanAppDims(ThisDim) .AND. TmpFlag ) THEN

                  ! Set an initial search index for the 3D and 4D array interpolation
               LastIndex3 = (/0,0,0/)
               LastIndex4 = (/0,0,0,0/)

                  ! To make things run slightly quicker, copy the data we will be interpolating over into the temporary arrays
               IF (NewmanAppData%DataIs3D) THEN
                  TmpData3D = NewmanAppData%Data3D%DataSet(:,:,:,Idx)
               ELSE
                  TmpData4D = NewmanAppData%Data4D%DataSet(:,:,:,:,Idx)
               END IF


               DO J=1,InitInp%NStepWave2

                     ! First get the wave amplitude -- must be reconstructed from the WaveElevC array.  First index is the real (1) or
                     ! imaginary (2) part.  Divide by NStepWave2 so that the wave amplitude is of the same form as the paper.
                  aWaveElevC = CMPLX( InitInp%WaveElevC0(1,J), InitInp%WaveElevC0(2,J), SiKi) / InitInp%NStepWave2

                     ! Calculate the frequency
                  Omega1 = J * InitInp%WaveDOmega


                     ! Only get a QTF value if within the range of frequencies between the cutoffs for the difference frequency
                  IF ( (Omega1 >= InitInp%WvLowCOff) .AND. (Omega1 <= InitInp%WvHiCOff) ) THEN

                        ! Now get the QTF value that corresponds to this frequency and wavedirection pair.
                     IF ( NewmanAppData%DataIs3D ) THEN

                           ! Set the (omega1,beta1,beta2) point we are looking for.
                        Coord3 = (/ REAL(Omega1,SiKi), InitInp%WaveDirArr(J), InitInp%WaveDirArr(J) /)

                           ! Apply local Z rotation to heading angle (degrees) to put wave direction into the local (rotated) body frame
                        Coord3(2) = Coord3(2) - RotateZdegOffset
                        Coord3(3) = Coord3(3) - RotateZdegOffset

                           ! get the interpolated value for F(omega1,beta1,beta2)
                        CALL WAMIT_Interp3D_Cplx( Coord3, TmpData3D, NewmanAppData%Data3D%WvFreq1, &
                                             NewmanAppData%Data3D%WvDir1, NewmanAppData%Data3D%WvDir2, LastIndex3, QTF_Value, ErrStatTmp, ErrMsgTmp )

                     ELSE

                           ! Set the (omega1,omega2,beta1,beta2) point we are looking for.
                        Coord4 = (/ REAL(Omega1,SiKi), REAL(Omega1,SiKi), InitInp%WaveDirArr(J), InitInp%WaveDirArr(J) /)

                           ! Apply local Z rotation to heading angle (degrees) to put wave direction into the local (rotated) body frame
                        Coord4(3) = Coord4(3) - RotateZdegOffset
                        Coord4(4) = Coord4(4) - RotateZdegOffset

                           ! get the interpolated value for F(omega1,omega2,beta1,beta2)
                        CALL WAMIT_Interp4D_Cplx( Coord4, TmpData4D, NewmanAppData%Data4D%WvFreq1, NewmanAppData%Data4D%WvFreq2, &
                                             NewmanAppData%Data4D%WvDir1, NewmanAppData%Data4D%WvDir2, LastIndex4, QTF_Value, ErrStatTmp, ErrMsgTmp )


                     ENDIF !QTF value find

                        ! Now we have the value of the QTF.  These values should only be real for the omega1=omega2 case of the approximation.
                        ! However if the value came from the 4D interpolation routine, it might have some residual complex part to it.  So
                        ! we throw the complex part out.  NOTE: the phase shift due to location will be added before the FFT.
                     QTF_Value = CMPLX(REAL(QTF_Value,SiKi),0.0,SiKi)


                  ELSE     ! outside the frequency range

                     QTF_Value = CMPLX(0.0,0.0,SiKi)

                  ENDIF    ! frequency check

                     ! Check and make sure nothing bombed in the interpolation that we need to be aware of
                  CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
                  IF ( ErrStat >= AbortErrLev ) THEN
                     IF (ALLOCATED(TmpData3D))        DEALLOCATE(TmpData3D,STAT=ErrStatTmp)
                     IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
                     IF (ALLOCATED(NewmanTerm1t))     DEALLOCATE(NewmanTerm1t,STAT=ErrStatTmp)
                     IF (ALLOCATED(NewmanTerm2t))     DEALLOCATE(NewmanTerm2t,STAT=ErrStatTmp)
                     IF (ALLOCATED(NewmanTerm1C))     DEALLOCATE(NewmanTerm1C,STAT=ErrStatTmp)
                     IF (ALLOCATED(NewmanTerm2C))     DEALLOCATE(NewmanTerm2C,STAT=ErrStatTmp)
                     IF (ALLOCATED(NewmanAppForce))   DEALLOCATE(NewmanAppForce,STAT=ErrStatTmp)
                     RETURN
                  ENDIF

                     ! Now calculate the Newman terms
                  IF (REAL(QTF_Value) > 0.0_SiKi) THEN

                     NewmanTerm1C(J,ThisDim) = aWaveElevC * (QTF_Value)**0.5_SiKi
                     NewmanTerm2C(J,ThisDim) = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)

                  ELSE IF (REAL(QTF_Value) < 0.0_SiKi) THEN

                     NewmanTerm1C(J,ThisDim) = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
                     NewmanTerm2C(J,ThisDim) = aWaveElevC * (-QTF_Value)**0.5_SiKi

                  ELSE ! at 0

                     NewmanTerm1C(J,ThisDim) = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)
                     NewmanTerm2C(J,ThisDim) = CMPLX(0.0_SiKi, 0.0_SiKi, SiKi)

                  ENDIF


               ENDDO    ! J=1,InitInp%NStepWave2

            ENDIF    ! Load component to calculate

         ENDDO ! ThisDim -- index to current dimension


         !----------------------------------------------------
         ! Rotate back to global frame and phase shift and set the terms for the summation
         !----------------------------------------------------

            ! Set rotation
            ! NOTE: RotateZMatrixT is the rotation from local to global.
         RotateZMatrixT(:,1) = (/  cos(InitInp%PtfmRefztRot(IBody)), -sin(InitInp%PtfmRefztRot(IBody)) /)
         RotateZMatrixT(:,2) = (/  sin(InitInp%PtfmRefztRot(IBody)),  cos(InitInp%PtfmRefztRot(IBody)) /)

            ! Loop through all the frequencies
         DO J=1,InitInp%NStepWave2

               ! Frequency
            Omega1 = J * InitInp%WaveDOmega

            !> Phase shift due to offset in location, only for NBodyMod==2
            if (p%NBodyMod == 2) then

               !> The phase shift due to an (x,y) offset is of the form
               !! \f$  exp[-\imath k(\omega) ( X cos(\beta(w)) + Y sin(\beta(w)) )] \f$
               !  NOTE: the phase shift applies to the aWaveElevC of the incoming wave.  Including it here instead
               !        of above is mathematically equivalent, but only because each frequency has only one wave
               !        direction associated with it through the equal energy approach used in multidirectional waves.

               WaveNmbr1   = WaveNumber ( REAL(Omega1,SiKi), InitInp%Gravity, InitInp%WtrDpth )    ! SiKi returned
               TmpReal1    = WaveNmbr1 * ( InitInp%PtfmRefxt(1)*cos(InitInp%WaveDirArr(J)*D2R) + InitInp%PtfmRefyt(1)*sin(InitInp%WaveDirArr(J)*D2R) )
               PhaseShiftXY = CMPLX( cos(TmpReal1), -sin(TmpReal1) )

               ! Apply the phase shift
               DO ThisDim=1,6
                  NewmanTerm1C(J,ThisDim) = NewmanTerm1C(J,ThisDim)*PhaseShiftXY       ! Newman term 1
                  NewmanTerm2C(J,ThisDim) = NewmanTerm2C(J,ThisDim)*PhaseShiftXY       ! Newman term 2
               ENDDO
            endif


               ! Apply the rotation to get back to global frame  -- Term 1
            NewmanTerm1C(J,1:2) = MATMUL(RotateZMatrixT, NewmanTerm1C(J,1:2))
            NewmanTerm1C(J,4:5) = MATMUL(RotateZMatrixT, NewmanTerm1C(J,4:5))

               ! Apply the rotation to get back to global frame  -- Term 2
            NewmanTerm2C(J,1:2) = MATMUL(RotateZMatrixT, NewmanTerm2C(J,1:2))
            NewmanTerm2C(J,4:5) = MATMUL(RotateZMatrixT, NewmanTerm2C(J,4:5))

         ENDDO    ! J=1,InitInp%NStepWave2



         !----------------------------------------------------
         ! Apply the FFT to get time domain results
         !----------------------------------------------------

         DO ThisDim=1,6    ! Loop through all dimensions

            Idx= (IBody-1)*6+ThisDim

               ! Now we apply the FFT to the first piece.
            CALL ApplyCFFT(  NewmanTerm1t(:), NewmanTerm1C(:,ThisDim), FFT_Data, ErrStatTmp )
            CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               IF (ALLOCATED(TmpData3D))        DEALLOCATE(TmpData3D,STAT=ErrStatTmp)
               IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanTerm1t))     DEALLOCATE(NewmanTerm1t,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanTerm2t))     DEALLOCATE(NewmanTerm2t,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanTerm1C))     DEALLOCATE(NewmanTerm1C,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanTerm2C))     DEALLOCATE(NewmanTerm2C,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanAppForce))   DEALLOCATE(NewmanAppForce,STAT=ErrStatTmp)
               RETURN
            END IF

               ! Now we apply the FFT to the second piece.
            CALL ApplyCFFT( NewmanTerm2t(:), NewmanTerm2C(:,ThisDim), FFT_Data, ErrStatTmp )
            CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               IF (ALLOCATED(TmpData3D))        DEALLOCATE(TmpData3D,STAT=ErrStatTmp)
               IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanTerm1t))     DEALLOCATE(NewmanTerm1t,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanTerm2t))     DEALLOCATE(NewmanTerm2t,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanTerm1C))     DEALLOCATE(NewmanTerm1C,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanTerm2C))     DEALLOCATE(NewmanTerm2C,STAT=ErrStatTmp)
               IF (ALLOCATED(NewmanAppForce))   DEALLOCATE(NewmanAppForce,STAT=ErrStatTmp)
               RETURN
            ENDIF


               ! Now square the real part of the resulting time domain pieces and add them together to get the final force time series.
            DO J=0,InitInp%NStepWave-1
               NewmanAppForce(J,Idx) = (abs(NewmanTerm1t(J)))**2 - (abs(NewmanTerm2t(J)))**2
            ENDDO

               ! Copy the last first term to the last so that it is cyclic
            NewmanAppForce(InitInp%NStepWave,Idx) = NewmanAppForce(0,Idx)

         ENDDO ! ThisDim -- index to current dimension

      ENDDO    ! IBody -- current body


         ! Done with the FFT library routines, so end them.
      CALL  ExitCFFT(FFT_Data, ErrStatTmp)
      CALL  SetErrStat(ErrStatTmp,'Error occured while cleaning up after the FFTs.', ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         IF (ALLOCATED(TmpData3D))        DEALLOCATE(TmpData3D,STAT=ErrStatTmp)
         IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm1t))     DEALLOCATE(NewmanTerm1t,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm2t))     DEALLOCATE(NewmanTerm2t,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm1C))     DEALLOCATE(NewmanTerm1C,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanTerm2C))     DEALLOCATE(NewmanTerm2C,STAT=ErrStatTmp)
         IF (ALLOCATED(NewmanAppForce))   DEALLOCATE(NewmanAppForce,STAT=ErrStatTmp)
         RETURN
      END IF


         ! Cleanup
      IF (ALLOCATED(TmpData3D))        DEALLOCATE(TmpData3D,STAT=ErrStatTmp)
      IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
      IF (ALLOCATED(NewmanTerm1t))     DEALLOCATE(NewmanTerm1t,STAT=ErrStatTmp)
      IF (ALLOCATED(NewmanTerm2t))     DEALLOCATE(NewmanTerm2t,STAT=ErrStatTmp)
      IF (ALLOCATED(NewmanTerm1C))     DEALLOCATE(NewmanTerm1C,STAT=ErrStatTmp)
      IF (ALLOCATED(NewmanTerm2C))     DEALLOCATE(NewmanTerm2C,STAT=ErrStatTmp)



   END SUBROUTINE NewmanApp_InitCalc





   !-------------------------------------------------------------------------------------------------------------------------------
   !> This subroutine calculates the force time series using the full DiffQTF calculation.
   !! The data is stored in a 4D array.  A 3D array also exists in the dataset and may be used by the MnDrift calculation when it
   !! is called from here.
   !! At each step in the summation of the mth term, a call is made to the 4D interpolation algorithm to find the value of
   !! \f$ F^-_k(\omega_m, \omega_n, \beta_1, \beta_2) \f$ corresponding to the \f$ Z(\omega_m) \f$ term in the complex wave
   !! amplitude, _aWaveElevC_. A combination of the limits of \f$ \omega_{lo-d} \le \omega_\mu \le  \omega_{hi-d} \f$ are
   !! imposed during the summation before the FFT and values outside this range are set to zero.  \f$ \omega_{\mu^-} \f$ is the
   !! difference frequency.
   !!
   !! For multi-directional waves where the equal energy discretization is used, each frequency has a single wave direction
   !! associated with it. In the full QTF calculation, the combination of the wave frequencies and their associated wave
   !! directions are used.
   !!
   !! The single summation equation used here is given by
   !!    \f$   {F_{{ex}~k}^{(-)(2)}} = \Re \left(
   !!                 \sum\limits_{m=0}^{N/2}   a_m a_m^* F_k^{-}(\omega_m,\omega_m) +
   !!              2 \cdot \sum\limits_{\mu^{-}=1}^{N/2-1} H_{\mu^{-}} \mathrm{e}^{i(\omega_{\mu^{-}})t}  \right)\qquad  \f$,
   !! for \f$  \qquad k=1,2,\ldots,6 \qquad\f$
   !! where  \f$ \qquad
   !!             H_{\mu^{-}} = \frac{1}{2} \sum\limits_{n=0}^{N/2-\mu^{-}} A_{\mu^{-}+n} A_n^* F_k^{-}(\omega_{\mu^{-}+n},\omega_n), \qquad \f$
   !!   for \f$ \qquad\quad 1\le\mu^{-}\le N-1 \f$.
   !!
   !! Note that \f$ k \f$ indicates the index to the load component,  \f$ {F_{{ex}~k}^{{-}(2)}} \f$ is the resulting second order
   !! force, and \f$ a_m \f$ and \f$ a_m^* \f$ is the complex wave amplitude and its conjugate for the \f$ m^{th} \f$ frequency.
   !! Note also that \f$ F_k^{-} (\omega_m, \omega_n, \beta_m, \beta_n) \f$ is the dimensionalized the QTF value read from the
   !! WAMIT file and interpolated for the \f$ m^{th} \f$ wave frequency. The value of \f$ \mu^{-}=m-n \f$, the difference index.
   !!
   !! The first term in this equation is the MnDrift force, so we will call that routine to get it. The second term is an IFFT
   !! of the summation of \f$ H_{\mu^{-}} \f$.  The limits that are imposed on this equation are the WvLowCOffD and WvHiCOffD
   !! values read in from the input file.  These limits are placed on \f$ \mu \f$ in the second term of the equation.
   !!
   !! Since the frequency range of the QTF has not yet been checked, we will do that now as well.
   !!
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE DiffQTF_InitCalc( InitInp, p, DiffQTFData, DiffQTFForce, ErrMsg, ErrStat )

      IMPLICIT NONE

      TYPE(WAMIT2_InitInputType),         INTENT(IN   )  :: InitInp              !< Input data for initialization routine
      TYPE(WAMIT2_ParameterType),         INTENT(IN   )  :: p                    !< Parameters
      TYPE(W2_DiffData_Type),             INTENT(INOUT)  :: DiffQTFData          !< Data storage for the DiffQTF method.  Set to INOUT in case we need to convert 4D to 3D
      REAL(SiKi),  ALLOCATABLE,           INTENT(  OUT)  :: DiffQTFForce(:,:)    !< Force data.  Index 1 is the timestep, index 2 is the load component.
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat

         ! Local Variables
      CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary error message for calls
      INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
      REAL(SiKi)                                         :: TmpReal1             !< Temporary real
      REAL(SiKi)                                         :: TmpReal2             !< Temporary real
      LOGICAL                                            :: TmpFlag              !< Temporary logical flag
      INTEGER(IntKi)                                     :: ThisDim              !< Generic counter for dimension
      INTEGER(IntKi)                                     :: IBody                !< Index to which body we are on
      INTEGER(IntKi)                                     :: Idx                  !< Index to the full set of 6*NBody
      INTEGER(IntKi)                                     :: J                    !< Generic counter
      INTEGER(IntKi)                                     :: K                    !< Generic counter
      TYPE(FFT_DataType)                                 :: FFT_Data             !< Temporary array for the FFT module we're using
      CHARACTER(*), PARAMETER                            :: RoutineName = 'DiffQTF_InitCalc'


         ! Wave information and QTF temporary
      COMPLEX(SiKi)                                      :: QTF_Value            !< Temporary complex number for QTF
      COMPLEX(SiKi), ALLOCATABLE                         :: TmpComplexArr(:,:)   !< Temporary complex array for frequency domain of one complete load component
      COMPLEX(SiKi)                                      :: TmpHMinusC           !< Temporary variable for holding the current value of \f$ H^- \f$
      COMPLEX(SiKi)                                      :: aWaveElevC1          !< Wave elevation of the first  frequency component.  NStepWave2 factor removed.
      COMPLEX(SiKi)                                      :: aWaveElevC2          !< Wave elevation of the second frequency component.  NStepWave2 factor removed.
      REAL(ReKi)                                         :: OmegaDiff            !< Wave difference frequency
      REAL(SiKi),    ALLOCATABLE                         :: TmpDiffQTFForce(:)   !< The resulting diffQTF force for this load component
      REAL(ReKi)                                         :: Omega1               !< First  wave frequency
      REAL(ReKi)                                         :: Omega2               !< Second wave frequency
      REAL(SiKi),    ALLOCATABLE                         :: MnDriftForce(:)      !< Mean drift force (first term).  MnDrift_InitCalc routine will return this.
      REAL(SiKi)                                         :: RotateZdegOffset     !< Offset to wave heading (NBodyMod==2 only)
      REAL(SiKi)                                         :: RotateZMatrixT(2,2)  !< The transpose of rotation in matrix form for rotation about z (from global to local)
      COMPLEX(SiKi)                                      :: PhaseShiftXY         !< The phase shift offset to apply to the body
      REAL(SiKi)                                         :: WaveNmbr1            !< Wavenumber for this frequency
      REAL(SiKi)                                         :: WaveNmbr2            !< Wavenumber for this frequency

         ! Interpolation routine indices and value to search for, and smaller array to pass
      INTEGER(IntKi)                                     :: LastIndex4(4)        !< Last used index for searching in the interpolation algorithms.  First  wave freq
      REAL(SiKi)                                         :: Coord4(4)            !< The (omega1,omega2,beta1,beta2) coordinate we want in the 4D dataset. First  wave freq.
      COMPLEX(SiKi), ALLOCATABLE                         :: TmpData4D(:,:,:,:)   !< Temporary 4D array we put the 4D data into (minus the load component indice)


         ! Initialize a few things
      ErrMsg      = ''
      ErrMsgTmp   = ''
      ErrStat     = ErrID_None
      ErrStatTmp  = ErrID_None


         !> 1. Check the data to see if the wave frequencies are present in the QTF data.

      IF ( DiffQTFData%DataIs4D ) THEN   ! We must have a 4D data set

             ! Check the low frequency cutoff
         IF ( MINVAL( DiffQTFData%Data4D%WvFreq1 ) > InitInp%WvLowCOffD ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(DiffQTFData%Data4D%WvFreq1)))// &
                           ' rad/s first wave period) data in '//TRIM(DiffQTFData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOffD.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF
         IF ( MINVAL( DiffQTFData%Data4D%WvFreq2 ) > InitInp%WvLowCOffD ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(DiffQTFData%Data4D%WvFreq2)))// &
                           ' rad/s for second wave period) data in '//TRIM(DiffQTFData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOffD.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF

            ! Check the high frequency cutoff -- using the Difference high frequency cutoff.  The first order high frequency
            ! cutoff is typically too high for this in most cases.
         IF ( MAXVAL(DiffQTFData%Data4D%WvFreq1) < InitInp%WvHiCOffD ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(DiffQTFData%Data4D%WvFreq1)))// &
                           ' rad/s for first wave period) data in '//TRIM(DiffQTFData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOffD.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF
         IF ( MAXVAL(DiffQTFData%Data4D%WvFreq2) < InitInp%WvHiCOffD ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(DiffQTFData%Data4D%WvFreq1)))// &
                           ' rad/s second wave period) data in '//TRIM(DiffQTFData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOffD.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF

      ELSE
            ! This is a catastrophic issue.  We should not have called this routine without data that is usable for the DiffQTF calculation
         CALL SetErrStat( ErrID_Fatal, ' The full Difference QTF method requires 4D data, and was not passed any.',ErrStat,ErrMsg,RoutineName)
      ENDIF

      IF ( ErrStat >= AbortErrLev )    RETURN




         ! If we are using multidirectional waves, then we should have more than 1 wave direction in the WAMIT file.
      IF ( InitInp%WaveMultiDir .AND. (DiffQTFData%Data4D%NumWvDir1 == 1) ) THEN
         CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(DiffQTFData%Filename)//' only contains one wave '// &
                     'direction at '//TRIM(Num2LStr(DiffQTFData%Data4D%WvDir1(1)))//' degrees (first wave direction). '// &
                     'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                     ErrStat,ErrMsg,RoutineName)
      ELSE IF ( InitInp%WaveMultiDir .AND. (DiffQTFData%Data4D%NumWvDir2 == 1) ) THEN
         CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(DiffQTFData%Filename)//' only contains one wave '// &
                     'direction at '//TRIM(Num2LStr(DiffQTFData%Data4D%WvDir2(1)))//' degrees (second wave direction). '// &
                     'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                     ErrStat,ErrMsg,RoutineName)
      ELSE

            ! See Known Issues #1 at the top of this file.  There may be problems if the data spans the +/- Pi boundary.  For
            ! now (since time is limited) we will issue a warning if any of the wave directions for multidirectional waves
            ! or data from the WAMIT file for the wavedirections is close to the +/-pi boundary (>150 degrees, <-150 degrees),
            ! we will issue a warning.
         IF ( (InitInp%WaveDirMin > 150.0_SiKi) .OR. (InitInp%WaveDirMax < -150.0_SiKi) .OR. &
              (MINVAL(DiffQTFData%Data4D%WvDir1) > 150.0_SiKi) .OR.  (MAXVAL(DiffQTFData%Data4D%WvDir1) < -150.0_SiKi) .OR. &
              (MINVAL(DiffQTFData%Data4D%WvDir2) > 150.0_SiKi) .OR.  (MAXVAL(DiffQTFData%Data4D%WvDir2) < -150.0_SiKi) ) THEN
            CALL SetErrStat( ErrID_Warn,' There may be issues with how the wave direction data is handled when the wave '// &
                                       'direction of interest is near the +/- 180 direction.  This is a known issue with '// &
                                       'the WAMIT2 module that has not yet been addressed.',ErrStat,ErrMsg,RoutineName)
         ENDIF

            ! Now check the limits for the first wave direction
            !  --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
            !  --> FIXME: modify to allow shifting values by TwoPi before comparing
         IF ( InitInp%WaveDirMin < MINVAL(DiffQTFData%Data4D%WvDir1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                  'found in the WAMIT data file '//TRIM(DiffQTFData%Filename)//' for the first wave direction.', &
                  ErrStat, ErrMsg, RoutineName)
         ENDIF
         IF ( InitInp%WaveDirMax > MAXVAL(DiffQTFData%Data4D%WvDir1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                  'found in the WAMIT data file '//TRIM(DiffQTFData%Filename)//' for the first wave direction.', &
                  ErrStat, ErrMsg, RoutineName)
         ENDIF


            ! Now check the limits for the second wave direction
            !   --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
         IF ( InitInp%WaveDirMin < MINVAL(DiffQTFData%Data4D%WvDir2) ) THEN
            CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                  'found in the WAMIT data file '//TRIM(DiffQTFData%Filename)//' for the second wave direction.', &
                  ErrStat, ErrMsg, RoutineName)
         ENDIF
         IF ( InitInp%WaveDirMax > MAXVAL(DiffQTFData%Data4D%WvDir2) ) THEN
            CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                  'found in the WAMIT data file '//TRIM(DiffQTFData%Filename)//' for the second wave direction.', &
                    ErrStat, ErrMsg, RoutineName)
         ENDIF

      ENDIF

      IF ( ErrStat >= AbortErrLev )    RETURN




         !> 4. Now check to make sure we have data that will work.  For the 4D data, it must not be sparse.
         !!    To check this, we have to check the load components that we will use.  So, we will loop through them
         !!    and set the TmpFlag to true if there is a sparse matrix for one of them.
         !! FIXME: remove this check and warning once the sparse matrix interpolation routines are implemented.
      TmpFlag = .FALSE.
      DO IBody=1,DiffQTFData%Data4D%NumBodies
         DO ThisDim=1,6
            Idx = (IBody-1)*6+ThisDim
            IF ( DiffQTFData%Data4D%DataIsSparse(Idx) .AND. DiffQTFData%Data4D%LoadComponents(Idx) .AND. p%DiffQTFDims(ThisDim) )    TmpFlag = .TRUE.
         ENDDO
      ENDDO
      IF (TmpFlag) THEN
         CALL SetErrStat(ErrID_Fatal,' The second order WAMIT data in '//TRIM(DiffQTFData%Filename)//' is too sparse '// &
               'for the interpolation routine used in the full Difference QTF calculation.  At some later point, we will allow for '// &
               'sparse data when a different interpolation routine is implemented.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF


         !> 5.  Calculate the second order forces
         !! The single summation equation used here is given by


         ! Setup temporary arrays that we will be passing the data to the interpolation routines in.
      ALLOCATE( TmpData4D(DiffQTFData%Data4D%NumWvFreq1,DiffQTFData%Data4D%NumWvFreq1,DiffQTFData%Data4D%NumWvDir1,DiffQTFData%Data4D%NumWvDir2),STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate temporary array for interpolation of 4D QTF data.', &
                                             ErrStat, ErrMsg, RoutineName)


         ! Setup the arrays holding the DiffQTF terms, both the complex frequency domain and real time domain pieces
      ALLOCATE( TmpDiffQTFForce( 0:InitInp%NStepWave), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for one load component of the full difference '// &
                                             'QTF 2nd order force time series.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( TmpComplexArr( 0:InitInp%NStepWave2, 6), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for one load component of the full difference '// &
                                             'QTF 2nd order force in the frequency domain.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( DiffQTFForce( 0:InitInp%NStepWave, 6*p%NBody), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for the full difference '// &
                                             'QTF 2nd order force time series.',ErrStat, ErrMsg, RoutineName)

         ! If something went wrong during allocation of the temporary arrays...
      IF ( ErrStat >= AbortErrLev ) THEN
         call cleanup()
         RETURN
      ENDIF


         ! Set the resulting force to zero.  Will calculate for the dimensions requested.
      DiffQTFForce = 0.0_SiKi


         ! Initialize the FFT library.  Do not apply normalization.
      CALL InitFFT ( InitInp%NStepWave, FFT_Data, .FALSE., ErrStatTmp )
      CALL SetErrStat(ErrStatTmp,'Error occured while initializing the FFT.',ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         call cleanup()
         RETURN
      END IF




         ! Before we continue, we will get the MnDriftForce results.
         !  --> Note that we can pass the DiffQTFData directly since we are using the same type for both
      CALL MnDrift_InitCalc( InitInp, p, DiffQTFData, MnDriftForce, ErrMsgTmp, ErrStatTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         call cleanup()
         RETURN
      ENDIF

         ! Make sure we have a value for the mean drift for each of the dimensions we were requested to calculate. To
         ! find out we will just check the status of the flags for each dimension in the MnDriftDims against the ones
         ! for the DiffQTFDims
      DO ThisDim=1,6
         IF ( p%DiffQTFDims(ThisDim) .AND. (.NOT. p%MnDriftDims(ThisDim)) ) &
            CALL SetErrStat( ErrID_Fatal,' The DiffQTF method requires the use of the MnDrift method for the first term. '// &
                     'Something went wrong and the MnDrift method returned a different number of load components.', &
                     ErrStat,ErrMsg,RoutineName)
      ENDDO
      IF ( ErrStat >= AbortErrLev ) THEN
         call cleanup()
         RETURN
      ENDIF



         ! Now loop through all the dimensions and perform the calculation
      DO IBody=1,p%NBody

            ! Initialize the temporary array to zero.
         TmpComplexArr = CMPLX(0.0_SiKi,0.0_SiKi,SiKi)

            ! Heading correction, only applies to NBodyMod == 2
         if (p%NBodyMod==2) then
            RotateZdegOffset = InitInp%PtfmRefztRot(IBody)*R2D
         else
            RotateZdegOffset = 0.0_SiKi
         endif

         !----------------------------------------------------
         ! Populate the frequency terms for this body
         !     -- with phase shift for NBodyMod == 2
         !----------------------------------------------------

         DO ThisDim=1,6
            Idx = (IBody-1)*6+ThisDim

               ! Only on the dimensions we requested, and it exists in the dataset
            IF ( p%DiffQTFDims(ThisDim) .AND. DiffQTFData%Data4D%LoadComponents(Idx) ) THEN


                  ! Set an initial search index for the 4D array interpolation
               LastIndex4 = (/0,0,0,0/)

                  ! To make things run slightly quicker, copy the data we will be interpolating over into the temporary arrays
               TmpData4D = DiffQTFData%Data4D%DataSet(:,:,:,:,Idx)

                  ! Outer loop to create the TmpComplexArr
               DO J=1,InitInp%NStepWave2-1

                     ! Calculate the frequency  -- This is the difference frequency.
                  OmegaDiff = J * InitInp%WaveDOmega


                     ! Only perform calculations if the difference frequency is in the right range
                  IF ( (OmegaDiff >= InitInp%WvLowCOffD) .AND. (OmegaDiff <= InitInp%WvHiCOffD) ) THEN

                        ! Set the \f$ H^- \f$ term to zero before we start
                     TmpHMinusC = CMPLX(0.0_SiKi,0.0_SiKi,SiKi)


                       ! Do the sum over H^-
                     DO K=1,InitInp%NStepWave2-J        ! note the funny upper limit.  This is because we are doing a summation on a triangular area.

                           ! set the two frequencies that the difference frequency comes from
                        Omega1 = (J + K) * InitInp%WaveDOmega        ! the mth frequency -- \mu^- + n = m
                        Omega2 = K * InitInp%WaveDOmega              ! the nth frequency

                           ! Find the Wave amplitudes 1 and 2
                        aWaveElevC1 = CMPLX( InitInp%WaveElevC0(1,J+K), InitInp%WaveElevC0(2,J+K), SiKi)  / InitInp%NStepWave2
                        aWaveElevC2 = CMPLX( InitInp%WaveElevC0(1,K),   InitInp%WaveElevC0(2,K),   SiKi)  / InitInp%NStepWave2

                           ! Set the (omega1,omega2,beta1,beta2) point we are looking for.
                        Coord4 = (/ REAL(Omega1,SiKi), REAL(Omega2,SiKi), InitInp%WaveDirArr(J+K), InitInp%WaveDirArr(K) /)

                           ! Apply local Z rotation to heading angle (degrees) to put wave direction into the local (rotated) body frame
                        Coord4(3) = Coord4(3) - RotateZdegOffset
                        Coord4(4) = Coord4(4) - RotateZdegOffset

                           ! get the interpolated value for F(omega1,omega2,beta1,beta2)  --> QTF_Value
                        CALL WAMIT_Interp4D_Cplx( Coord4, TmpData4D, DiffQTFData%Data4D%WvFreq1, DiffQTFData%Data4D%WvFreq2, &
                                             DiffQTFData%Data4D%WvDir1, DiffQTFData%Data4D%WvDir2, LastIndex4, QTF_Value, ErrStatTmp, ErrMsgTmp )
                        CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
                        IF (ErrStat >= AbortErrLev ) THEN
                           call cleanup()
                           RETURN
                        ENDIF

                        !--------------------------
                        ! Phase shift due to offset
                        !--------------------------
                        if (p%NBodyMod == 2) then
                           !> The phase shift due to an (x,y) offset for second order difference frequencies is of the form
                           !! \f$  exp[-\imath ( k(\omega_1) ( X cos(\beta(w_1)) + Y sin(\beta(w_1)) )
                           !!                  - k(\omega_2) ( X cos(\beta(w_2)) + Y sin(\beta(w_2)) ) ) ]\f$
                           !  NOTE: the phase shift applies to the aWaveElevC of the incoming wave.  Including it here instead
                           !        of above is mathematically equivalent, but only because each frequency has only one wave
                           !        direction associated with it through the equal energy approach used in multidirectional waves.

                           WaveNmbr1   = WaveNumber ( REAL(Omega1,SiKi), InitInp%Gravity, InitInp%WtrDpth )    ! SiKi returned
                           WaveNmbr2   = WaveNumber ( REAL(Omega2,SiKi), InitInp%Gravity, InitInp%WtrDpth )    ! SiKi returned
                           TmpReal1    = WaveNmbr1 * ( InitInp%PtfmRefxt(1)*cos(InitInp%WaveDirArr(J+K)*D2R) + InitInp%PtfmRefyt(1)*sin(InitInp%WaveDirArr(J+K)*D2R) )
                           TmpReal2    = WaveNmbr2 * ( InitInp%PtfmRefxt(1)*cos(InitInp%WaveDirArr(K)*D2R)   + InitInp%PtfmRefyt(1)*sin(InitInp%WaveDirArr(K)*D2R)   )

                           ! Set the phase shift for the set of difference frequencies
                           PhaseShiftXY = CMPLX( cos(TmpReal1 - TmpReal2), -sin(TmpReal1 - TmpReal2) )

                           ! For similicity, apply to the QTF_Value (mathematically equivalent to applying to the wave elevations)
                           QTF_Value = QTF_Value*PhaseShiftXY

                        endif ! Phaseshift for NBodyMod==2

                           ! Calculate this value and add it to what we have so far.
                        TmpHMinusC = TmpHMinusC + aWaveElevC1 * CONJG(aWaveElevC2) * QTF_Value

                     ENDDO

                        ! Copy this value difference frequency information over to the array we will take the IFFT.  Divide
                        ! by two for the single sided FFT given in the documentation.
                     TmpComplexArr(J,ThisDim) = TmpHMinusC / 2.0_SiKi

                  ELSE     ! outside the frequency range, so

                     TmpComplexArr(J,ThisDim) = CMPLX(0.0_SiKi,0.0_SiKi,SiKi)

                  ENDIF    ! frequency check


               ENDDO
            ENDIF    ! Load component to calculate
         ENDDO ! ThisDim -- The current dimension


         !----------------------------------------------------
         ! Rotate back to global frame
         !----------------------------------------------------

            ! Set rotation
            ! NOTE: RotateZMatrixT is the rotation from local to global.
         RotateZMatrixT(:,1) = (/  cos(InitInp%PtfmRefztRot(IBody)), -sin(InitInp%PtfmRefztRot(IBody)) /)
         RotateZMatrixT(:,2) = (/  sin(InitInp%PtfmRefztRot(IBody)),  cos(InitInp%PtfmRefztRot(IBody)) /)

            ! Loop through all the frequencies
         DO J=1,InitInp%NStepWave2

               ! Apply the rotation to get back to global frame
            TmpComplexArr(J,1:2) = MATMUL(RotateZMatrixT, TmpComplexArr(J,1:2))
            TmpComplexArr(J,4:5) = MATMUL(RotateZMatrixT, TmpComplexArr(J,4:5))

         ENDDO    ! J=1,InitInp%NStepWave2



         !----------------------------------------------------
         ! Apply the FFT to get time domain results
         !----------------------------------------------------

         DO ThisDim=1,6
            Idx = (IBody-1)*6+ThisDim

               ! Now we apply the FFT to the result of the sum
            CALL ApplyFFT_cx(  TmpDiffQTFForce(:),  TmpComplexArr(:,ThisDim), FFT_Data, ErrStatTmp )
            CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to the second term of the difference QTF.', &
                           ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               call cleanup()
               RETURN
            END IF


               ! Now we multiply the result by 2 and save it to the DiffQTFForce array and add the MnDrift term
               ! NOTE: phase shift and orientations on the MnDriftForce term have already been applied
               ! NOTE: the "-1" since TmpDiffQTFForce(InitInp%NStepWave) is not set and DiffQTFForce(InitInp%NStepWave,Idx) gets overwritten
            DO K=0,InitInp%NStepWave-1
               DiffQTFForce(K,Idx) = 2.0_SiKi * TmpDiffQTFForce(K) + MnDriftForce(Idx)
            ENDDO

               ! Copy the last first term to the first so that it is cyclic
            DiffQTFForce(InitInp%NStepWave,Idx) = DiffQTFForce(0,Idx)

         ENDDO ! ThisDim -- The current dimension
      ENDDO    ! IBody -- This WAMIT body



         ! Done with the FFT library routines, so end them.
      CALL  ExitFFT(FFT_Data, ErrStatTmp)
      CALL  SetErrStat(ErrStatTmp,'Error occured while cleaning up after the FFTs.', ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         call cleanup()
         RETURN
      END IF

      call cleanup()
      
   contains
!--------------------------------------------------- 
      subroutine cleanup()

            ! Cleanup
         IF (ALLOCATED(TmpData4D))           DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
         IF (ALLOCATED(TmpDiffQTFForce))     DEALLOCATE(TmpDiffQTFForce,STAT=ErrStatTmp)
         IF (ALLOCATED(TmpComplexArr))       DEALLOCATE(TmpComplexArr,STAT=ErrStatTmp)
      end subroutine cleanup
!--------------------------------------------------- 
   END SUBROUTINE DiffQTF_InitCalc








   !-------------------------------------------------------------------------------------------------------------------------------
   !> This subroutine calculates the force time series using the full SumQTF calculation. The data is stored in a 4D array.
   !! At each step in the summation of the mth term, a call is made to the 4D interpolation algorithm to find the value of
   !! \f$ F^-_k(\omega_m, \omega_n, \beta_1, \beta_2) \f$ corresponding to the \f$ Z(\omega_m) \f$ term in the complex wave
   !! amplitude, _WaveElevC_. A combination of the limits of \f$ \omega_{lo-s} \le \omega_\mu \le  \omega_{hi-s} \f$ are
   !! imposed during the summation before the FFT and values outside this range are set to zero.
   !!
   !! For multi-directional waves where the equal energy discretization is used, each frequency has a single wave direction
   !! associated with it. In the full QTF calculation, the combination of the wave frequencies and their associated wave
   !! directions are used.
   !!
   !! The single summation equation used here is given by
   !!    \f$   {F_{{ex}~k}^{(+)(2)}} = \Re  \left(
   !!          \sum\limits_{m=1}^{\lfloor N/4 \rfloor}   A_m A_m F_k^{+}(\omega_m,\omega_m) \mathrm{e}^{2\cdot i\omega_m t} +
   !!                2 \sum\limits_{\mu^{+}=2}^{N/2} H_{\mu^{+}}\mathrm{e}^{i (\omega_{\mu^{+}-1})t}\right)  ,\quad         \f$
   !!      for   \f$\quad k=1,2,\ldots,6 \qquad\f$
   !!
   !!   where  \f$
   !!             H_{\mu^{+}} = \sum\limits_{n=1}^{\left\lfloor\frac{\mu^{+}-n}{2}\right\rfloor}
   !!                                     A_n A_{\mu^{+}-n} F_k^{+}(\omega_n,\omega_{\mu^{+}-1}),\f$
   !!                                                                for \f$ 2\leq \mu^{+}\leq N/2+1 \f$
   !!
   !!    and   \f$
   !!             H_{\mu^{+}} = \displaystyle\sum\limits_{n=\mu^{+}-N}^{\left\lfloor\frac{\mu^{+}-n}{2}\right\rfloor}
   !!                                     A_n A_{\mu^{+}-n} F_k^{+}(\omega_n,\omega_{\mu^{+}-1}),\f$
   !!                                                                for \f$ N/2+2\leq\mu^{+}\leq N \f$
   !!
   !!    where \f$ \mu^{+} = m + n \f$, \f$ \lfloor x \rfloor \f$ represents the floor function given by
   !!          \f$  \lfloor x \rfloor \equiv \max \left\{ m \in \mathbb{Z} \middle| m\leq x \right\}.\f$
   !!
   !! Note that \f$ k \f$ indicates the index to the load component,  \f$ {F_{{ex}~k}^{{+}(2)}} \f$ is the resulting second order
   !! force, and \f$ A_m \f$ is the complex wave amplitude for the \f$ m^{th} \f$ frequency.
   !! Note also that \f$ F_k^{+} (\omega_m, \omega_n, \beta_m, \beta_n) \f$ is the dimensionalized QTF value read from the
   !! WAMIT file and interpolated for the \f$ m^{th} \f$ wave frequency.
   !!
   !! The first term in this equation is an IFFT of the diagonal elements where \f$ \omega_m = \omega_n \f$. The second term is
   !! an IFFT of the summation of \f$ H_{\mu^{+}} \f$ term.  The limits that are imposed on this equation are the WvLowCOffS and
   !! WvHiCOffS values read in from the input file.  These limits are placed on \f$ \mu^+ \f$ in the both terms of the equation.
   !!
   !! The second term is an IFFT over \f$ N \f$ terms instead of the usual \f$ N/2 \f$ terms.  This results in \f$ 2N \f$ timesteps
   !! with a timestep of \f$ \Delta t / 2 \f$.  However, since we only report at \f$ \Delta t \f$ timestep intervals, the highest
   !! reportable frequency is the Nyquist frequency of \f$ \pi / \Delta t \f$, so anything above that will be lost.  Considering
   !! this, there is no reason to calculate any summation frequencies above the Nyquist frequency.  So we will change the equation
   !! to the following form:
   !!    \f$   {F_{{ex}~k}^{(+)(2)}} = \Re  \left(
   !!          \sum\limits_{m=1}^{\lfloor N/4 \rfloor}   A_m A_m F_k^{+}(\omega_m,\omega_m) \mathrm{e}^{2\cdot i\omega_m t} +
   !!                2 \sum\limits_{\mu^{+}=2}^{N/2} H_{\mu^{+}}\mathrm{e}^{i (\omega_{\mu^{+}-1})t}\right),        \f$
   !!      for   \f$ k=1,2,\ldots,6 \f$
   !!
   !!   where  \f$
   !!             H_{\mu^{+}} = \sum\limits_{n=1}^{\left\lfloor\frac{\mu^{+}-n}{2}\right\rfloor}
   !!                                     A_n A_{\mu^{+}-n} F_k^{+}(\omega_n,\omega_{\mu^{+}-1}), \f$
   !!                                                                for \f$ 2\le \mu^{+}\le N/2 \f$
   !!
   !!
   !!
   !!
   !! Before doing all the calculations, we will first check the frequency range of the QTF and check if we will be truncating data
   !! if the cutoff for the highest sum frequency, _WvHiCOffS_ is above the Nyquist frequency.
   !!
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE SumQTF_InitCalc( InitInp, p, SumQTFData, SumQTFForce, ErrMsg, ErrStat )

      IMPLICIT NONE

      TYPE(WAMIT2_InitInputType),         INTENT(IN   )  :: InitInp              !< Input data for initialization routine
      TYPE(WAMIT2_ParameterType),         INTENT(IN   )  :: p                    !< Parameters
      TYPE(W2_SumData_Type),              INTENT(INOUT)  :: SumQTFData           !< Data storage for the SumQTF method.  Set to INOUT in case we need to convert 4D to 3D
      REAL(SiKi),  ALLOCATABLE,           INTENT(  OUT)  :: SumQTFForce(:,:)     !< Force data.  Index 1 is the timestep, index 2 is the load component.
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat

         ! Local Variables
      CHARACTER(2048)                                    :: ErrMsgTmp            !< Temporary error message for calls
      INTEGER(IntKi)                                     :: ErrStatTmp           !< Temporary error status for calls
      LOGICAL                                            :: TmpFlag              !< Temporary logical flag
      INTEGER(IntKi)                                     :: ThisDim              !< Generic counter for dimension
      INTEGER(IntKi)                                     :: IBody                !< Index to which body we are on
      INTEGER(IntKi)                                     :: Idx                  !< Index to the full set of 6*NBody
      INTEGER(IntKi)                                     :: J                    !< Generic counter
      INTEGER(IntKi)                                     :: K                    !< Generic counter
      TYPE(FFT_DataType)                                 :: FFT_Data             !< Temporary array for the FFT module we're using. For the first  term in the equation.
      CHARACTER(*), PARAMETER                            :: RoutineName = 'SumQTF_InitCalc'


         ! Wave information and QTF temporary
      COMPLEX(SiKi)                                      :: QTF_Value           !< Temporary complex number for QTF
      COMPLEX(SiKi), ALLOCATABLE                         :: Term1ArrayC(:,:)     !< Temporary complex array for frequency domain of one load component. For first  term.
      COMPLEX(SiKi), ALLOCATABLE                         :: Term2ArrayC(:,:)     !< Temporary complex array for frequency domain of one load component. For second term.
      REAL(SiKi),    ALLOCATABLE                         :: Term1Array(:)        !< Temporary complex array for time      domain of one load component. For first  term.
      REAL(SiKi),    ALLOCATABLE                         :: Term2Array(:)        !< Temporary complex array for time      domain of one load component. For second term.
      COMPLEX(SiKi)                                      :: TmpHPlusC            !< Temporary variable for holding the current value of \f$ H^+ \f$
      COMPLEX(SiKi)                                      :: aWaveElevC1          !< Wave elevation of the first  frequency component.  NStepWave2 factor removed.
      COMPLEX(SiKi)                                      :: aWaveElevC2          !< Wave elevation of the second frequency component.  NStepWave2 factor removed.
      REAL(ReKi)                                         :: OmegaSum             !< Wave difference frequency
      REAL(ReKi)                                         :: Omega1               !< First  wave frequency
      REAL(ReKi)                                         :: Omega2               !< Second wave frequency
      REAL(SiKi)                                         :: RotateZdegOffset     !< Offset to wave heading (NBodyMod==2 only)
      REAL(SiKi)                                         :: RotateZMatrixT(2,2)  !< The transpose of rotation in matrix form for rotation about z (from global to local)
      COMPLEX(SiKi)                                      :: PhaseShiftXY         !< The phase shift offset to apply to the body
      REAL(SiKi)                                         :: WaveNmbr1            !< Wavenumber for this frequency
      REAL(SiKi)                                         :: WaveNmbr2            !< Wavenumber for this frequency
      REAL(SiKi)                                         :: TmpReal1             !< Temporary real
      REAL(SiKi)                                         :: TmpReal2             !< Temporary real

         ! Interpolation routine indices and value to search for, and smaller array to pass
      INTEGER(IntKi)                                     :: LastIndex4(4)        !< Last used index for searching in the interpolation algorithms.  First  wave freq
      REAL(SiKi)                                         :: Coord4(4)            !< The (omega1,omega2,beta1,beta2) coordinate we want in the 4D dataset. First  wave freq.
      COMPLEX(SiKi), ALLOCATABLE                         :: TmpData4D(:,:,:,:)   !< Temporary 4D array we put the 4D data into (minus the load component indice)


         ! Initialize a few things
      ErrMsg      = ''
      ErrMsgTmp   = ''
      ErrStat     = ErrID_None
      ErrStatTmp  = ErrID_None


         !> 1. Check the data to see if the wave frequencies are present in the QTF data.

      IF ( SumQTFData%DataIs4D ) THEN   ! We must have a 4D data set

             ! Check the low frequency cutoff
         IF ( MINVAL( SumQTFData%Data4D%WvFreq1 ) > InitInp%WvLowCOffS ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(SumQTFData%Data4D%WvFreq1)))// &
                           ' rad/s first wave period) data in '//TRIM(SumQTFData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOffS.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF
         IF ( MINVAL( SumQTFData%Data4D%WvFreq2 ) > InitInp%WvLowCOffS ) THEN
            CALL SetErrStat( ErrID_Fatal,' The lowest frequency ( '//TRIM(Num2LStr(MINVAL(SumQTFData%Data4D%WvFreq2)))// &
                           ' rad/s for second wave period) data in '//TRIM(SumQTFData%Filename)// &
                           ' is above the low frequency cutoff set by WvLowCOffS.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF

            ! Check the high frequency cutoff -- using the Difference high frequency cutoff.  The first order high frequency
            ! cutoff is typically too high for this in most cases.
         IF ( MAXVAL(SumQTFData%Data4D%WvFreq1) < InitInp%WvHiCOffS ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(SumQTFData%Data4D%WvFreq1)))// &
                           ' rad/s for first wave period) data in '//TRIM(SumQTFData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOffS.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF
         IF ( MAXVAL(SumQTFData%Data4D%WvFreq2) < InitInp%WvHiCOffS ) THEN
            CALL SetErrStat( ErrID_Fatal,' The highest frequency ( '//TRIM(Num2LStr(MAXVAL(SumQTFData%Data4D%WvFreq1)))// &
                           ' rad/s second wave period) data in '//TRIM(SumQTFData%Filename)// &
                           ' is below the high frequency cutoff set by WvHiCOffS.', &
                           ErrStat,ErrMsg,RoutineName)
         ENDIF

      ELSE
            ! This is a catastrophic issue.  We should not have called this routine without data that is usable for the SumQTF calculation
         CALL SetErrStat( ErrID_Fatal, ' The full Sum QTF method requires 4D data, and was not passed any.',ErrStat,ErrMsg,RoutineName)
      ENDIF

      IF ( ErrStat >= AbortErrLev )    RETURN




         ! If we are using multidirectional waves, then we should have more than 1 wave direction in the WAMIT file.
      IF ( InitInp%WaveMultiDir .AND. (SumQTFData%Data4D%NumWvDir1 == 1) ) THEN
         CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(SumQTFData%Filename)//' only contains one wave '// &
                     'direction at '//TRIM(Num2LStr(SumQTFData%Data4D%WvDir1(1)))//' degrees (first wave direction). '// &
                     'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                     ErrStat,ErrMsg,RoutineName)
      ELSE IF ( InitInp%WaveMultiDir .AND. (SumQTFData%Data4D%NumWvDir2 == 1) ) THEN
         CALL SetErrStat( ErrID_Fatal,' WAMIT output file '//TRIM(SumQTFData%Filename)//' only contains one wave '// &
                     'direction at '//TRIM(Num2LStr(SumQTFData%Data4D%WvDir2(1)))//' degrees (second wave direction). '// &
                     'It cannot be used with multidirectional waves.  Set WaveDirMod to 0 to use this file.', &
                     ErrStat,ErrMsg,RoutineName)
      ELSE

            ! See Known Issues #1 at the top of this file.  There may be problems if the data spans the +/- Pi boundary.  For
            ! now (since time is limited) we will issue a warning if any of the wave directions for multidirectional waves
            ! or data from the WAMIT file for the wavedirections is close to the +/-pi boundary (>150 degrees, <-150 degrees),
            ! we will issue a warning.
         IF ( (InitInp%WaveDirMin > 150.0_SiKi) .OR. (InitInp%WaveDirMax < -150.0_SiKi) .OR. &
              (MINVAL(SumQTFData%Data4D%WvDir1) > 150.0_SiKi) .OR.  (MAXVAL(SumQTFData%Data4D%WvDir1) < -150.0_SiKi) .OR. &
              (MINVAL(SumQTFData%Data4D%WvDir2) > 150.0_SiKi) .OR.  (MAXVAL(SumQTFData%Data4D%WvDir2) < -150.0_SiKi) ) THEN
            CALL SetErrStat( ErrID_Warn,' There may be issues with how the wave direction data is handled when the wave '// &
                                       'direction of interest is near the +/- 180 direction.  This is a known issue with '// &
                                       'the WAMIT2 module that has not yet been addressed.',ErrStat,ErrMsg,RoutineName)
         ENDIF

            ! Now check the limits for the first wave direction
            !  --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
            !  --> FIXME: modify to allow shifting values by TwoPi before comparing
         IF ( InitInp%WaveDirMin < MINVAL(SumQTFData%Data4D%WvDir1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                  'found in the WAMIT data file '//TRIM(SumQTFData%Filename)//' for the first wave direction.', &
                  ErrStat, ErrMsg, RoutineName)
         ENDIF
         IF ( InitInp%WaveDirMax > MAXVAL(SumQTFData%Data4D%WvDir1) ) THEN
            CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                  'found in the WAMIT data file '//TRIM(SumQTFData%Filename)//' for the first wave direction.', &
                  ErrStat, ErrMsg, RoutineName)
         ENDIF


            ! Now check the limits for the second wave direction
            !   --> FIXME: sometime fix this to handle the above case.  See Known Issue #1 at top of file.
         IF ( InitInp%WaveDirMin < MINVAL(SumQTFData%Data4D%WvDir2) ) THEN
            CALL SetErrStat( ErrID_Fatal,' Minimum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMin))//' is not'//&
                  'found in the WAMIT data file '//TRIM(SumQTFData%Filename)//' for the second wave direction.', &
                  ErrStat, ErrMsg, RoutineName)
         ENDIF
         IF ( InitInp%WaveDirMax > MAXVAL(SumQTFData%Data4D%WvDir2) ) THEN
            CALL SetErrStat( ErrID_Fatal,' Maximum wave direction required of '//TRIM(Num2LStr(InitInp%WaveDirMax))//' is not'//&
                  'found in the WAMIT data file '//TRIM(SumQTFData%Filename)//' for the second wave direction.', &
                    ErrStat, ErrMsg, RoutineName)
         ENDIF

      ENDIF

      IF ( ErrStat >= AbortErrLev )    RETURN




         !> 4. Now check to make sure we have data that will work.  For the 4D data, it must not be sparse.
         !!    To check this, we have to check the load components that we will use.  So, we will loop through them
         !!    and set the TmpFlag to true if there is a sparse matrix for one of them.
         !! FIXME: remove this check and warning once the sparse matrix interpolation routines are implimented.
      TmpFlag = .FALSE.
      DO IBody=1,SumQTFData%Data4D%NumBodies
         DO ThisDim=1,6
            Idx = (IBody-1)*6+ThisDim
            IF ( SumQTFData%Data4D%DataIsSparse(Idx) .AND. SumQTFData%Data4D%LoadComponents(Idx) .AND. p%SumQTFDims(ThisDim) )    TmpFlag = .TRUE.
         ENDDO
      ENDDO
      IF (TmpFlag) THEN
         CALL SetErrStat(ErrID_Fatal,' The second order WAMIT data in '//TRIM(SumQTFData%Filename)//' is too sparse '// &
               'for the interpolation routine used in the full sum QTF calculation.  At some later point, we will allow for '// &
               'sparse data when a different interpolation routine is implimented.',ErrStat,ErrMsg,RoutineName)
         RETURN
      ENDIF


         !> 5.  Calculate the second order forces
         !! The single summation equation used here is given by


         ! Setup temporary arrays that we will be passing the data to the interpolation routines in.
      ALLOCATE( TmpData4D(SumQTFData%Data4D%NumWvFreq1,SumQTFData%Data4D%NumWvFreq1,SumQTFData%Data4D%NumWvDir1,SumQTFData%Data4D%NumWvDir2),STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate temporary array for interpolation of 4D QTF data.', &
                                             ErrStat, ErrMsg, RoutineName)


         ! Setup the arrays holding the SumQTF terms, both the complex frequency domain and real time domain pieces
      ALLOCATE( Term1ArrayC( 0:InitInp%NStepWave2, 6), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for the first term of one load component of the full sum '// &
                                             'QTF 2nd order force in the frequency domain.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( Term2ArrayC( 0:InitInp%NStepWave2, 6), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for the second term of one load component of the full sum '// &
                                             'QTF 2nd order force in the frequency domain.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( Term1Array( 0:InitInp%NStepWave), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for the first term of one load component of the full sum '// &
                                             'QTF 2nd order force in the time domain.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( Term2Array( 0:InitInp%NStepWave), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for the second term of one load component of the full sum '// &
                                             'QTF 2nd order force in the time domain.',ErrStat, ErrMsg, RoutineName)
      ALLOCATE( SumQTFForce( 0:InitInp%NStepWave, 6*p%NBody), STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,' Cannot allocate array for the full difference '// &
                                             'QTF 2nd order force time series.',ErrStat, ErrMsg, RoutineName)

         ! If something went wrong during allocation of the temporary arrays...
      IF ( ErrStat >= AbortErrLev ) THEN
         IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
         IF (ALLOCATED(Term1ArrayC))      DEALLOCATE(Term1ArrayC,STAT=ErrStatTmp)
         IF (ALLOCATED(Term2ArrayC))      DEALLOCATE(Term2ArrayC,STAT=ErrStatTmp)
         IF (ALLOCATED(Term1Array))       DEALLOCATE(Term1Array,STAT=ErrStatTmp)
         IF (ALLOCATED(Term2Array))       DEALLOCATE(Term2Array,STAT=ErrStatTmp)
         RETURN
      ENDIF


         ! Set the resulting force to zero.  Will calculate for the dimensions requested.
      SumQTFForce = 0.0_SiKi


         ! Initialize the FFT library.  Normalization not required in this formulation.
      CALL InitFFT ( InitInp%NStepWave, FFT_Data, .FALSE., ErrStatTmp )        ! FIXME:
      CALL SetErrStat(ErrStatTmp,'Error occured while initializing the FFT.',ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
         IF (ALLOCATED(SumQTFForce))      DEALLOCATE(SumQTFForce,STAT=ErrStatTmp)
         RETURN
      END IF



         ! Now loop through all the dimensions and perform the calculation
      DO IBody=1,p%NBody

            ! Initialize the temporary array to zero.
         Term1ArrayC = CMPLX(0.0_SiKi,0.0_SiKi,SiKi)
         Term2ArrayC = CMPLX(0.0_SiKi,0.0_SiKi,SiKi)

            ! Heading correction, only applies to NBodyMod == 2
         if (p%NBodyMod==2) then
            RotateZdegOffset = InitInp%PtfmRefztRot(IBody)*R2D
         else
            RotateZdegOffset = 0.0_SiKi
         endif

         !----------------------------------------------------
         ! Populate the frequency terms for this body
         !     -- with phase shift for NBodyMod == 2
         !----------------------------------------------------

         DO ThisDim=1,6
            Idx = (IBody-1)*6+ThisDim

               ! Only on the dimensions we requested, and if it is present in the data
            IF ( p%SumQTFDims(ThisDim) .AND. SumQTFData%Data4D%LoadComponents(Idx) ) THEN


                  ! To make things run slightly quicker, copy the data we will be interpolating over into the temporary arrays
               TmpData4D = SumQTFData%Data4D%DataSet(:,:,:,:,Idx)

               !---------------------------------------------------------------------------------
               ! Calculate the first term
               ! This term is only the FFT over the diagonal elements where omega_1 = omega_2
               ! note that the sum frequency is 2*omega.  The index for the sum frequency is
               ! therefore 2*J.  Since we are placing the calculated value for the A_m * A_m *
               ! F_k^+ term in the 2*omega location, we will only run through the first half of
               ! the frequencies (the sum frequency will exceed the bounds of the frequencies
               ! used in the FFT otherwise).
               ! The IFFT will be calculated later.

                  ! Set an initial search index for the 4D array interpolation
               LastIndex4 = (/0,0,0,0/)

                  ! The limits look a little funny.  But remember we are placing the value in the 2*J location,
                  ! so we cannot overun the end of the array, and the highest frequency must be zero.  The
                  ! floor function is just in case (NStepWave2 - 1) is an odd number
               DO J=1,FLOOR(REAL(InitInp%NStepWave2-1)/2.0_SiKi)

                     ! The frequency
                  Omega1   = REAL(J,ReKi) * InitInp%WaveDOmega
                  OmegaSum = 2.0_SiKi * Omega1           ! the sum frequency

                     ! Only perform calculations if the difference frequency is in the right range
                  IF ( (OmegaSum >= InitInp%WvLowCOffS) .AND. (OmegaSum <= InitInp%WvHiCOffS) ) THEN


                        ! Find the wave amplitude at frequency omega
                     aWaveElevC1 = CMPLX( InitInp%WaveElevC0(1,J), InitInp%WaveElevC0(2,J), SiKi )  / InitInp%NStepWave2

                        ! Set the (omega1,omega2,beta1,beta2) point we are looking for.
                     Coord4 = (/ REAL(Omega1,SiKi), REAL(Omega1,SiKi), InitInp%WaveDirArr(J), InitInp%WaveDirArr(J) /)

                        ! Apply local Z rotation to heading angle (degrees) to put wave direction into the local (rotated) body frame
                     Coord4(3) = Coord4(3) - RotateZdegOffset
                     Coord4(4) = Coord4(4) - RotateZdegOffset

                        ! get the interpolated value for F(omega1,omega2,beta1,beta2)  --> QTF_Value
                     CALL WAMIT_Interp4D_Cplx( Coord4, TmpData4D, SumQTFData%Data4D%WvFreq1, SumQTFData%Data4D%WvFreq2, &
                                          SumQTFData%Data4D%WvDir1, SumQTFData%Data4D%WvDir2, LastIndex4, QTF_Value, ErrStatTmp, ErrMsgTmp )
                     CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
                     IF (ErrStat >= AbortErrLev ) THEN
                        IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
                        RETURN
                     ENDIF

                     !--------------------------
                     ! Phase shift due to offset
                     !--------------------------
                     if (p%NBodyMod == 2) then
                        !> The phase shift due to an (x,y) offset for second order difference frequencies is of the form
                        !! \f$  exp[-\imath ( k(\omega_1) ( X cos(\beta(w_1)) + Y sin(\beta(w_1)) )
                        !!                  1 k(\omega_2) ( X cos(\beta(w_2)) + Y sin(\beta(w_2)) ) ) ]\f$.
                        !! For the first term, \f$ \omega_1 = \omega_2 \$f.
                        !  NOTE: the phase shift applies to the aWaveElevC of the incoming wave.  Including it here instead
                        !        of above is mathematically equivalent, but only because each frequency has only one wave
                        !        direction associated with it through the equal energy approach used in multidirectional waves.

                        WaveNmbr1   = WaveNumber ( REAL(Omega1,SiKi), InitInp%Gravity, InitInp%WtrDpth )    ! SiKi returned
                        TmpReal1    = WaveNmbr1 * ( InitInp%PtfmRefxt(1)*cos(InitInp%WaveDirArr(J)*D2R) + InitInp%PtfmRefyt(1)*sin(InitInp%WaveDirArr(J)*D2R) )

                        ! Set the phase shift for the set of sum frequencies
                        PhaseShiftXY = CMPLX( cos(TmpReal1 + TmpReal1), -sin(TmpReal1 + TmpReal1) )

                        ! For similicity, apply to the QTF_Value (mathematically equivalent to applying to the wave elevations)
                        QTF_Value = QTF_Value*PhaseShiftXY

                     endif ! Phaseshift for NBodyMod==2

                        ! Set the value of the first term in the frequency domain
                     Term1ArrayC(2*J,ThisDim) = aWaveElevC1 * aWaveElevC1 * QTF_Value


                  ENDIF    ! Check on the limits
               ENDDO       ! First term calculation


               !---------------------------------------------------------------------------------
               ! Calculate the second term.
               ! In this term, we are are now stepping through the sum frequencies.  The inner
               ! sum essentially covers all the off diagonal terms (omega_m /= omega_n).  The limits
               ! on the outer integral that is the FFT run through the full frequency range that
               ! we are using


                  ! Set an initial search index for the 4D array interpolation
               LastIndex4 = (/0,0,0,0/)

                  ! Check the limits for the high frequency cutoff. If WvHiCOffS is less than the
                  ! maximum frequency possible with the value of WaveDT (omega_max = pi/WaveDT = NStepWave2*WaveDOmega),
                  ! then we are good.  If the WvHiCOff > 1/2 omega_max, then we will be potentially
                  ! throwing away information.  However, remember the following:
                  !     WaveDT   Omega_max    wavelength
                  !      (s)      (rad/s)        (m)
                  !       .25       4 Pi           < 1
                  !      0.5        2 Pi           1
                  !      1.0          Pi          10
                  ! so, we don't need a really small WaveDT

         !This section has been removed since it is kind of annoying.
         !      IF ( InitInp%WvHiCOffS > InitInp%NStepWave2*InitInp%WaveDOmega ) THEN
         !         CALL SetErrStat( ErrID_Warn,' The high frequency cutoff for second order wave forces, WvHiCOffS, '// &
         !                  'is larger than the Nyquist frequency for the given time step of WaveDT. The Nyquist frequency '// &
         !                  '(highest frequency) that can be computed is OmegaMax = PI/WaveDT = '// &
         !                  TRIM(Num2LStr(InitInp%NStepWave2*InitInp%WaveDOmega))// &
         !                  ' radians/second.  If you need those frequencies, decrease WaveDT.  For reference, 2*PI '// &
         !                  'radians/second corresponds to a wavelength of ~1 meter.',&
         !                  ErrStat,ErrMsg,RoutineName)
         !      ENDIF



                  ! Outer loop to create the Term2ArrayC. This is stepwise through the sum frequencies.
               DO J=1,InitInp%NStepWave2

                     ! Calculate the frequency  -- This is the sum frequency.
                  OmegaSum = J * InitInp%WaveDOmega



                     ! Set the \f$ H^+ \f$ term to zero before we start
                  TmpHPlusC = CMPLX(0.0_SiKi,0.0_SiKi,SiKi)


                     ! Only perform calculations if the difference frequency is in the right range
                  IF ( (OmegaSum >= InitInp%WvLowCOffS) .AND. (OmegaSum <= InitInp%WvHiCOffS) ) THEN

                     !> Now do the inner sum.  We are going to perform a sum up to the maximum frequency that we
                     !! can support (Nyquist frequency) for the given WaveDOmega and NStepWave2 (WaveOmegaMax =
                     !! NStepWave2 * WaveDOmega = Pi / WaveDT rad/second). Note that this means the largest diagonal
                     !! frequency we use is \f$ \omega = \Delta\omega * \f$ _NStepwave_/4.
                     !! So, we essentially end up running into the sampling limit.  If we want higher frequency
                     !! terms, we need to use a smaller stepsize.

                     DO K=0,FLOOR(Real(J-1)/2.0_SiKi)

                           ! Calculate the frequency pair
                        Omega1 =    K  * InitInp%WaveDOmega
                        Omega2 = (J-K) * InitInp%WaveDOmega

                           ! Find the wave amplitude at frequency omega.  Remove the NStepWave2 normalization built into WaveElevC0 from Waves module
                        aWaveElevC1 = CMPLX( InitInp%WaveElevC0(1,  K), InitInp%WaveElevC0(2,  K), SiKi )    / InitInp%NStepWave2
                        aWaveElevC2 = CMPLX( InitInp%WaveElevC0(1,J-K), InitInp%WaveElevC0(2,J-K), SiKi )    / InitInp%NStepWave2

                           ! Set the (omega1,omega2,beta1,beta2) point we are looking for.
                        Coord4 = (/ REAL(Omega1,SiKi), REAL(Omega2,SiKi), InitInp%WaveDirArr(K), InitInp%WaveDirArr(J-K) /)

                           ! Apply local Z rotation to heading angle (degrees) to put wave direction into the local (rotated) body frame
                        Coord4(3) = Coord4(3) - RotateZdegOffset
                        Coord4(4) = Coord4(4) - RotateZdegOffset

                           ! get the interpolated value for F(omega1,omega2,beta1,beta2)  --> QTF_Value
                        CALL WAMIT_Interp4D_Cplx( Coord4, TmpData4D, SumQTFData%Data4D%WvFreq1, SumQTFData%Data4D%WvFreq2, &
                                             SumQTFData%Data4D%WvDir1, SumQTFData%Data4D%WvDir2, LastIndex4, QTF_Value, ErrStatTmp, ErrMsgTmp )
                        CALL SetErrStat(ErrStatTmp,ErrMsgTmp,ErrStat,ErrMsg,RoutineName)
                        IF (ErrStat >= AbortErrLev ) THEN
                           IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
                           RETURN
                        ENDIF

                        !--------------------------
                        ! Phase shift due to offset
                        !--------------------------
                        if (p%NBodyMod == 2) then
                           !> The phase shift due to an (x,y) offset for second order difference frequencies is of the form
                           !! \f$  exp[-\imath ( k(\omega_1) ( X cos(\beta(w_1)) + Y sin(\beta(w_1)) )
                           !!                  - k(\omega_2) ( X cos(\beta(w_2)) + Y sin(\beta(w_2)) ) ) ]\f$
                           !  NOTE: the phase shift applies to the aWaveElevC of the incoming wave.  Including it here instead
                           !        of above is mathematically equivalent, but only because each frequency has only one wave
                           !        direction associated with it through the equal energy approach used in multidirectional waves.

                           WaveNmbr1   = WaveNumber ( REAL(Omega1,SiKi), InitInp%Gravity, InitInp%WtrDpth )    ! SiKi returned
                           WaveNmbr2   = WaveNumber ( REAL(Omega2,SiKi), InitInp%Gravity, InitInp%WtrDpth )    ! SiKi returned
                           TmpReal1    = WaveNmbr1 * ( InitInp%PtfmRefxt(1)*cos(InitInp%WaveDirArr(K)*D2R)   + InitInp%PtfmRefyt(1)*sin(InitInp%WaveDirArr(K)*D2R)   )
                           TmpReal2    = WaveNmbr2 * ( InitInp%PtfmRefxt(1)*cos(InitInp%WaveDirArr(J-K)*D2R) + InitInp%PtfmRefyt(1)*sin(InitInp%WaveDirArr(J-K)*D2R) )

                           ! Set the phase shift for the set of sum frequencies
                           PhaseShiftXY = CMPLX( cos(TmpReal1 + TmpReal2), -sin(TmpReal1 + TmpReal2) )

                           QTF_Value = QTF_Value*PhaseShiftXY

                        endif ! Phaseshift for NBodyMod==2

                           ! Set the value of the first term in the frequency domain.
                        TmpHPlusC = TmpHPlusC + aWaveElevC1 * aWaveElevC2 * QTF_Value

                     ENDDO

                        ! Save the value from the summation.
                     Term2ArrayC(J,ThisDim) = TmpHPlusC

                  ENDIF    ! Check on the limits

               ENDDO       ! Second term calculation -- frequency step on the sum frequency
            ENDIF    ! Load component to calculate
         ENDDO ! ThisDim -- current dimension


         !----------------------------------------------------
         ! Rotate back to global frame
         !----------------------------------------------------

            ! Set rotation
            ! NOTE: RotateZMatrixT is the rotation from local to global.
         RotateZMatrixT(:,1) = (/  cos(InitInp%PtfmRefztRot(IBody)), -sin(InitInp%PtfmRefztRot(IBody)) /)
         RotateZMatrixT(:,2) = (/  sin(InitInp%PtfmRefztRot(IBody)),  cos(InitInp%PtfmRefztRot(IBody)) /)

            ! Loop through all the frequencies
         DO J=1,InitInp%NStepWave2

               ! Apply the rotation to get back to global frame -- term 1
            Term1ArrayC(J,1:2) = MATMUL(RotateZMatrixT, Term1ArrayC(J,1:2))
            Term1ArrayC(J,4:5) = MATMUL(RotateZMatrixT, Term1ArrayC(J,4:5))

               ! Apply the rotation to get back to global frame -- term 2
            Term2ArrayC(J,1:2) = MATMUL(RotateZMatrixT, Term2ArrayC(J,1:2))
            Term2ArrayC(J,4:5) = MATMUL(RotateZMatrixT, Term2ArrayC(J,4:5))

         ENDDO    ! J=1,InitInp%NStepWave2



         !----------------------------------------------------
         ! Apply the FFT to get time domain results
         !----------------------------------------------------

         DO ThisDim=1,6
            Idx = (IBody-1)*6+ThisDim

               ! Now we apply the FFT to the result of the sum.
            CALL ApplyFFT_cx(  Term1Array(:),  Term1ArrayC(:,ThisDim), FFT_Data, ErrStatTmp )
            CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to the first term of the Sum QTF.', &
                           ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
               RETURN
            END IF

               ! Now we apply the FFT to the result of the sum.
            CALL ApplyFFT_cx(  Term2Array(:), Term2ArrayC(:,ThisDim), FFT_Data, ErrStatTmp )
            CALL SetErrStat(ErrStatTmp,'Error occured while applying the FFT to the second term of the Sum QTF.', &
                           ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) THEN
               IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
               RETURN
            ENDIF

               ! Now we add the two terms together.  The 0.5 multiplier on is because the double sided FFT was used.
            DO J=0,InitInp%NStepWave-1  !bjj: Term1Array and Term2Array don't set the last element, so we can get over-flow errors here. SumQTFForce(InitInp%NStepWave,Idx) gets overwritten later, so Idx'm setting the array bounds to be -1.
               SumQTFForce(J,Idx) = 0.5_SiKi*(REAL(Term1Array(J) + 2*Term2Array(J), SiKi))
            ENDDO

               ! Copy the last first term to the first so that it is cyclic
            SumQTFForce(InitInp%NStepWave,Idx) = SumQTFForce(0,Idx)


         ENDDO ! ThisDim -- current dimension
      ENDDO    ! IBody -- current WAMIT body



         ! Done with the FFT library routines, so end them.
      CALL  ExitFFT(FFT_Data, ErrStatTmp)
      CALL  SetErrStat(ErrStatTmp,'Error occured while cleaning up after the FFTs.', ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         IF (ALLOCATED(TmpData4D))        DEALLOCATE(TmpData4D,STAT=ErrStatTmp)
         RETURN
      END IF


         ! Cleanup
      IF (ALLOCATED(TmpData4D))           DEALLOCATE(TmpData4D,STAT=ErrStatTmp)

   END SUBROUTINE SumQTF_InitCalc










   !-------------------------------------------------------------------------------------------------------------------------------
   !> This subroutine checks the input values present in InitInput. This should in theory have been done in the HydroDyn_Input
   !! routine, so it should not be necessary to do this again here.  However, we don't necessarily trust that it has been done,
   !! particularly if this module gets used without HydroDyn.
   !!
   !! This subroutine also populates the InitOut and creates the filenames for each of the calculation types.
   !!
   SUBROUTINE CheckInitInput( InitInp, p, MnDriftData, NewmanAppData, DiffQTFData, SumQTFData, ErrStat, ErrMsg )

      IMPLICIT NONE

         ! Passed variables.
      TYPE(WAMIT2_InitInputType),         INTENT(IN   )  :: InitInp        !< Input data for initialization routine
      TYPE(WAMIT2_ParameterType),         INTENT(  OUT)  :: p              !< The parameters

         ! QTF storage -- from the data files that are read in
      TYPE(W2_DiffData_Type)                             :: MnDriftData    !< Data storage for the Mean Drift method
      TYPE(W2_DiffData_Type)                             :: NewmanAppData  !< Data storage for the Newman's Approximation method
      TYPE(W2_DiffData_Type)                             :: DiffQTFData    !< Data storage for the full difference QTF method
      TYPE(W2_SumData_Type)                              :: SumQTFData     !< Data storage for the full sum QTF method

      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< The error value
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< A message about the error.

         ! Temporary variables
      LOGICAL                                            :: TmpFileExist   !< Temporary flag to indicate existence of file

         ! Temporary Error Variables
      INTEGER(IntKi)                                     :: ErrStatTmp     !< Temporary variable for the local error status
!      CHARACTER(2048)                                    :: ErrMsgTmp      !< Temporary error message variable
      CHARACTER(*), PARAMETER                            :: RoutineName = 'CheckInitInput'

      !> ## Subroutine contents

      !--------------------------------------------------------------------------------
      !> ### Initialize variables
      !--------------------------------------------------------------------------------

         !> 1. Initialize error variables
      ErrStat     = ErrID_None
      ErrMsg      = ''

         !> 2. Initialize filenames
      MnDriftData%Filename   = ''
      NewmanAppData%Filename = ''
      DiffQTFData%Filename   = ''
      SumQTFData%Filename    = ''


      !--------------------------------------------------------------------------------
      !> ### Check the algorithm flags and file values:
      !!
      !!  In the following tests, we check to make sure that a good set of values is
      !!  given for the tests to run. The following chart shows what is allowed.
      !!
      !!  |  Algorithm    |   Flag True    |   Flag False   |   Invalid combination |
      !!  | :-----------: | :------------: | :------------: | :-------------------: |
      !!  | MnDrift       |    7 - 12      |      0         |     any other         |
      !!  | NewmanApp     |    7 - 12      |      0         |     any other         |
      !!  | DiffQTF       |   10 - 12      |      0         |     any other         |
      !!  | SumQTF        |   10 - 12      |      0         |     any other         |
      !!
      !!  Also, we check that only one of MnDrift, NewmanApp, or DiffQTF are nonzero.
      !!  Rather than test on the flag value (MnDriftF etc.), we test on the input value
      !!  then verify the flag matches. It could be done based on the flag, but it isn't
      !!  for arbitrarily chosen reasons.
      !!
      !!  @note The .8 input files only support dimensions 1,2, and 6 (surge, sway,
      !!        yaw). This gets checked in this routine, and a set of parameters are set
      !!        for each method listing which dimensions to use.
      !!
      !--------------------------------------------------------------------------------


         !> 1. Check that we only specified one of MnDrift, NewmanApp, or DiffQTF
         !!        (compared pairwise -- if any two are both true, we have a problem)

      IF ( ( InitInp%MnDrift /= 0 .AND. InitInp%NewmanApp /= 0 ) .OR. &
           ( InitInp%DiffQTF /= 0 .AND. InitInp%NewmanApp /= 0 ) .OR. &
           ( InitInp%MnDrift /= 0 .AND. InitInp%DiffQTF   /= 0 ) ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Only one of MnDrift, NewmanApp, or DiffQTF can be non-zero.', ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
      END IF


         !> 2. Check that we have a valid values for MnDrift, check flag status

      IF ( InitInp%MnDrift == 0 ) THEN             ! We are not doing anything
         IF ( InitInp%MnDriftF ) THEN    ! Should be false in this case, so if true there is a problem
            CALL SetErrStat( ErrID_Fatal, ' Programming error in call to WAMIT2_Init: '//NewLine// &
                  '           MnDriftF flag should be set to false by calling program for MnDrift = 0.'//NewLine// &
                  '              --> This should have been checked by the calling program.',ErrStat, ErrMsg, RoutineName)
         END IF
      ELSE IF( InitInp%MnDrift >= 7 .AND. InitInp%MnDrift <= 12 ) THEN              ! Valid values
         IF ( .not. InitInp%MnDriftF ) THEN    ! Should be true in this case, so if false there is a problem
            CALL SetErrStat( ErrID_Fatal, ' Programming error in call to WAMIT2_Init: '//NewLine// &
                  '           MnDriftF flag should be set to true by calling program for MnDrift /= 0.'//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
         END IF
         IF ( InitInp%NBody > 1 .AND. InitInp%MnDrift == 8 ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Mean drift calculation cannot be used with input file type 8 when more than 1 WAMIT body is present.', &
                                 ErrStat, ErrMsg, RoutineName)
         END IF
      ELSE
         CALL SetErrStat( ErrID_Fatal, ' Programming Error in call to WAMIT2_Init: '//NewLine// &
                  '           MnDrift can only have values of 0, 7, 8, 9, 10, 11, or 12. '//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
      END IF
      IF ( ErrStat >= AbortErrLev ) RETURN


         !> 3. Check that we have a valid values for NewmanApp, check flag status

      IF ( InitInp%NewmanApp == 0 ) THEN             ! We are not doing anything
         IF ( InitInp%NewmanAppF ) THEN    ! Should be false in this case, so if true there is a problem
            CALL SetErrStat( ErrID_Fatal, ' Programming error in call to WAMIT2_Init: '//NewLine// &
                  '           NewmanAppF flag should be set to false by calling program for NewmanApp = 0.'//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
         END IF
      ELSE IF( InitInp%NewmanApp >= 7 .AND. InitInp%NewmanApp <= 12 ) THEN              ! Valid values
         IF ( InitInp%NewmanAppF .eqv. .FALSE. ) THEN    ! Should be true in this case, so if false there is a problem
            CALL SetErrStat( ErrID_Fatal, ' Programming error in call to WAMIT2_Init: '//NewLine// &
                  '           NewmanAppF flag should be set to true by calling program for NewmanApp /= 0.'//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
         END IF
         IF ( InitInp%NBody > 1 .AND. InitInp%NewmanApp == 8 ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Newman''s approximation cannot be used with input file type 8 when more than 1 WAMIT body is present.', &
                                 ErrStat, ErrMsg, RoutineName)
         END IF
      ELSE
         CALL SetErrStat( ErrID_Fatal, ' Programming Error in call to WAMIT2_Init: '//NewLine// &
                  '           NewmanApp can only have values of 0, 7, 8, 9, 10, 11, or 12. '//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
      END IF
      IF ( ErrStat >= AbortErrLev ) RETURN


         !> 4. Check that we have a valid values for DiffQTF, check flag status

      IF ( InitInp%DiffQTF == 0 ) THEN             ! We are not doing anything
         IF ( InitInp%DiffQTFF .eqv. .TRUE. ) THEN    ! Should be false in this case, so if true there is a problem
            CALL SetErrStat( ErrID_Fatal, ' Programming error in call to WAMIT2_Init: '//NewLine// &
                  '           DiffQTFF flag should be set to false by calling program for DiffQTF = 0.'//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
         END IF
      ELSE IF( InitInp%DiffQTF >= 10 .AND. InitInp%DiffQTF <= 12 ) THEN              ! Valid values
         IF ( InitInp%DiffQTFF .eqv. .FALSE. ) THEN    ! Should be true in this case, so if false there is a problem
            CALL SetErrStat( ErrID_Fatal, ' Programming error in call to WAMIT2_Init: '//NewLine// &
                  '           DiffQTFF flag should be set to true by calling program for DiffQTF /= 0.'//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
         END IF
      ELSE
         CALL SetErrStat( ErrID_Fatal, ' Programming Error in call to WAMIT2_Init: '//NewLine// &
                  '           DiffQTF can only have values of 0, 10, 11, or 12. '//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
      END IF
      IF ( ErrStat >= AbortErrLev ) RETURN


         !> 5. Check that we have a valid values for SumQTF, check flag status

      IF ( InitInp%SumQTF == 0 ) THEN             ! We are not doing anything
         IF ( InitInp%SumQTFF .eqv. .TRUE. ) THEN    ! Should be false in this case, so if true there is a problem
            CALL SetErrStat( ErrID_Fatal, ' Programming error in call to WAMIT2_Init: '//NewLine// &
                  '           SumQTFF flag should be set to false by calling program for SumQTF = 0.'//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
         END IF
      ELSE IF( InitInp%SumQTF >= 10 .AND. InitInp%SumQTF <= 12 ) THEN              ! Valid values
         IF ( InitInp%SumQTFF .eqv. .FALSE. ) THEN    ! Should be true in this case, so if false there is a problem
            CALL SetErrStat( ErrID_Fatal, ' Programming error in call to WAMIT2_Init: '//NewLine// &
                  '           SumQTFF flag should be set to true by calling program for SumQTF /= 0.'//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
         END IF
      ELSE
         CALL SetErrStat( ErrID_Fatal, ' Programming Error in call to WAMIT2_Init: '//NewLine// &
                  '           SumQTF can only have values of 0, 10, 11, or 12. '//NewLine// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
      END IF
      IF ( ErrStat >= AbortErrLev ) RETURN


         !--------------------------------------------------------------------------------
         !> ### Check the Min and Max frequencies for the full QTF cases
         !!
         !!  -- these checks are performed based on the DiffQTFF and SumQTFF flags
         !--------------------------------------------------------------------------------


         !> 1. Check that the min / max diff frequencies make sense if using DiffQTF

      IF ( InitInp%DiffQTFF .eqv. .TRUE. ) THEN
         IF ( ( InitInp%WvHiCOffD < InitInp%WvLowCOffD ) .OR. ( InitInp%WvLowCOffD < 0.0 ) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming Error in call to WAMIT2_Init: '//NewLine// &
                  '           WvHiCOffD must be larger than WvLowCOffD. Both must be positive.'// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
            RETURN
         END IF
      END IF


         !> 2. Check that the min / max diff frequencies make sense if using SumQTF

      IF ( InitInp%SumQTFF .eqv. .TRUE. ) THEN
         IF ( ( InitInp%WvHiCOffS < InitInp%WvLowCOffS ) .OR. ( InitInp%WvLowCOffS < 0.0 ) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming Error in call to WAMIT2_Init: '//NewLine// &
                  '           WvHiCOffS must be larger than WvLowCOffS. Both must be positive.'// &
                  '              --> This should have been checked by the calling program.', ErrStat, ErrMsg, RoutineName)
            RETURN
         END IF
      END IF



         !--------------------------------------------------------------------------------
         !> ### Assemble the names of the WAMIT data files we are using and verify existence
         !--------------------------------------------------------------------------------

         !> 1. Check MnDrift file if we are doing the mean drift calculations.  Set flag for data type.

      MnDriftData%DataIs3D = .FALSE.
      MnDriftData%DataIs4D = .FALSE.

      IF ( InitInp%MnDrift /= 0) THEN
         IF ( InitInp%MnDrift <= 9 ) THEN                                                          ! For file types 7, 8, 9
            MnDriftData%Filename  = TRIM(InitInp%WAMITFile)//'.'//TRIM(Num2LStr(InitInp%MnDrift))
            INQUIRE( file=TRIM(MnDriftData%Filename), exist=TmpFileExist )
            MnDriftData%DataIs3D = .TRUE.
         ELSE                                                                                      ! For full QTF file types 10, 11, 12
            MnDriftData%Filename  = TRIM(InitInp%WAMITFile)//'.'//TRIM(Num2LStr(InitInp%MnDrift))//'d'
            INQUIRE( file=TRIM(MnDriftData%Filename), exist=TmpFileExist )
            MnDriftData%DataIs4D = .TRUE.
         ENDIF
         IF ( .not. TmpFileExist ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Cannot find the WAMIT file '//TRIM(MnDriftData%Filename)// &
                        ' required by the MnDrift option.', ErrStat, ErrMsg, RoutineName)
            RETURN
         END IF
      END IF


         !> 2. Check NewmanApp file if we are doing the Newman's approximation calculations. Set flag for data type.

      NewmanAppData%DataIs3D  = .FALSE.
      NewmanAppData%DataIs4D  = .FALSE.

      IF ( InitInp%NewmanApp /= 0) THEN
         IF ( InitInp%NewmanApp <= 9 ) THEN                                                          ! For file types 7, 8, 9
            NewmanAppData%Filename  = TRIM(InitInp%WAMITFile)//'.'//TRIM(Num2LStr(InitInp%NewmanApp))
            INQUIRE( file=TRIM(NewmanAppData%Filename), exist=TmpFileExist )
            NewmanAppData%DataIs3D  = .TRUE.
         ELSE                                                                                      ! For full QTF file types 10, 11, 12
            NewmanAppData%Filename  = TRIM(InitInp%WAMITFile)//'.'//TRIM(Num2LStr(InitInp%NewmanApp))//'d'
            INQUIRE( file=TRIM(NewmanAppData%Filename), exist=TmpFileExist )
            NewmanAppData%DataIs4D  = .TRUE.
         ENDIF
         IF ( .not. TmpFileExist ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Cannot find the WAMIT file '//TRIM(NewmanAppData%Filename)// &
                        ' required by the NewmanApp option.', ErrStat, ErrMsg, RoutineName)
            RETURN
         END IF
      END IF


         !> 3. Check DiffQTF file if we are doing the Difference QTF calculations

      DiffQTFData%DataIs3D  = .FALSE.
      DiffQTFData%DataIs4D  = .FALSE.

      IF ( InitInp%DiffQTF /= 0) THEN
         DiffQTFData%Filename  = TRIM(InitInp%WAMITFile)//'.'//TRIM(Num2LStr(InitInp%DiffQTF))//'d'
         INQUIRE( file=TRIM(DiffQTFData%Filename), exist=TmpFileExist )
         IF ( .not. TmpFileExist ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Cannot find the WAMIT file '//TRIM(DiffQTFData%Filename)// &
                        ' required by the DiffQTF option.', ErrStat, ErrMsg, RoutineName)
            RETURN
         END IF
         DiffQTFData%DataIs4D  = .TRUE.
      END IF


         !> 4. Check SumQTF file if we are doing the Sum QTF calculations

      SumQTFData%DataIs4D  = .FALSE.

      IF ( InitInp%SumQTF /= 0) THEN
         SumQTFData%Filename  = TRIM(InitInp%WAMITFile)//'.'//TRIM(Num2LStr(InitInp%SumQTF))//'s'
         INQUIRE( file=TRIM(SumQTFData%Filename), exist=TmpFileExist )
         IF ( .not. TmpFileExist ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Cannot find the WAMIT file '//TRIM(SumQTFData%Filename)// &
                        ' required by the SumQTF option.', ErrStat, ErrMsg, RoutineName)
            RETURN
         END IF
         SumQTFData%DataIs4D  = .TRUE.
      END IF


         !--------------------------------------------------------------------------------
         !> ### Check the size of arrays that were passed in containing the wave info
         !--------------------------------------------------------------------------------


         !> 1. Check that WaveElevC0 is a 2x(NStepWave2+1) sized array (0 index start)

      IF ( SIZE( InitInp%WaveElevC0, 2 ) /= (InitInp%NStepWave2 + 1) ) THEN    ! Expect a 2x(0:NStepWave2) array
         CALL SetErrStat( ErrID_Fatal, ' Programming error in call to WAMIT2_Init:'//NewLine// &
                     '        --> Expected array for WaveElevC0 to be of size 2x'//TRIM(Num2LStr(InitInp%NStepWave2 + 1))// &
                     ' (2x(NStepWave2+1)), but instead received array of size '// &
                     TRIM(Num2LStr(SIZE(InitInp%WaveElevC0,1)))//'x'//TRIM(Num2LStr(SIZE(InitInp%WaveElevC0,2)))//'.', ErrStat, ErrMsg, RoutineName)
         RETURN
      END IF




      !--------------------------------------------------------------------------------
      !> ### Now copy over things to parameters...
      !--------------------------------------------------------------------------------
      !> 1. Wave information we need to keep
      !--------------------------------------------------------------------------------
      p%NStepWave = InitInp%NStepWave


      !--------------------------------------------------------------------------------
      !> 3. WAMIT body related information
      !--------------------------------------------------------------------------------

      p%NBody                 =  InitInp%NBody              ! Number of bodies WAMIT2 sees
      p%NBodyMod              =  InitInp%NBodyMod           ! How multiple bodys are treated

         ! This module's implementation requires that if NBodyMod = 2 or 3, then there is one instance of a WAMIT module for each body, therefore, HydroDyn may have NBody > 1, but this WAMIT module will have NBody = 1
      if ( (p%NBodyMod > 1) .and. (p%NBody > 1) ) then
         CALL SetErrStat( ErrID_Fatal, "DEVELOPER ERROR: If NBodyMod = 2 or 3, then NBody for the a WAMIT2 object must be equal to 1", ErrStat, ErrMsg, RoutineName)
         return
      end if


      !--------------------------------------------------------------------------------
      !> ### Flags indicating forces to calculate for with each method
      !!     This information is stored on a method basis since we might use a .8 file
      !!     in either the NewmanApp or the MnDrift method (it doesn't contain all the
      !!     dimensions), and use a different file with the SumQTF method.
      !!
      !! The dimension numbers map as follows:
      !!  | Index | Name     | Axis       |
      !!  | :---: | :------: | :--------: |
      !!  |   1   |  Surge   | (X)        |
      !!  |   2   |  Sway    | (Y)        |
      !!  |   3   |  Heave   | (Z)        |
      !!  |   4   |  Roll    | (about X)  |
      !!  |   5   |  Pitch   | (about Y)  |
      !!  |   6   |  Yaw     | (about Z)  |
      !!
      !! @note
      !!  If we had set flags to do calculations of Heave, Roll, or Pitch (3, 4, 5) for
      !!  either the Mean Drift or Newman's Approximation methods while trying to use
      !!  the WAMIT .8 output files, we will issue a warning.
      !--------------------------------------------------------------------------------

         !> 1. For the Mean Drift method,

      IF (InitInp%MnDriftF)  THEN                        ! if the flag is true, we are doing this calculation
         IF (InitInp%MnDrift == 8) THEN                  ! the .8 files are not complete
            p%MnDriftDims(1)     =  .TRUE.
            p%MnDriftDims(2)     =  .TRUE.
            p%MnDriftDims(3)     =  .FALSE.              ! the .8 files don't contain this dimension
            p%MnDriftDims(4)     =  .FALSE.              ! the .8 files don't contain this dimension
            p%MnDriftDims(5)     =  .FALSE.              ! the .8 files don't contain this dimension
            p%MnDriftDims(6)     =  .TRUE.
         ELSE
            p%MnDriftDims        = .TRUE.
         ENDIF
      ELSE
         p%MnDriftDims(:)        = .FALSE.               ! Set all dimensions to false unless we are actually calculating something
      ENDIF



         !> 2. For the Newman Approximation method

      IF (InitInp%NewmanAppF)  THEN                      ! if the flag is true, we are doing this calculation
         IF (InitInp%NewmanApp == 8) THEN                ! the .8 files are not complete
            p%NewmanAppDims(1)   =  .TRUE.
            p%NewmanAppDims(2)   =  .TRUE.
            p%NewmanAppDims(3)   =  .FALSE.              ! the .8 files don't contain this dimension
            p%NewmanAppDims(4)   =  .FALSE.              ! the .8 files don't contain this dimension
            p%NewmanAppDims(5)   =  .FALSE.              ! the .8 files don't contain this dimension
            p%NewmanAppDims(6)   =  .TRUE.
         ELSE
            p%NewmanAppDims      = .TRUE.
         ENDIF
      ELSE
         p%NewmanAppDims(:)      = .FALSE.               ! Set all dimensions to false unless we are actually calculating something
      ENDIF



         !> 3. For the Difference QTF method,

      IF (InitInp%DiffQTFF)  THEN               ! if the flag is true, we are doing this calculation
         p%DiffQTFDims     = .TRUE.
            ! Also set the MnDrift flags.  We will be passing data from the DiffQTF method to the MnDrift method for the first term
         p%MnDriftDims     = .TRUE.
      ELSE
         p%DiffQTFDims(:)  = .FALSE.            ! Set all dimensions to false unless we are actually calculating something
      ENDIF


         !> 4. For the Summation QTF method,

      IF (InitInp%SumQTFF)  THEN                 ! if the flag is true, we are doing this calculation
         p%SumQTFDims      = .TRUE.
      ELSE
         p%SumQTFDims(:)   = .FALSE.            ! Set all dimensions to false unless we are actually calculating something
      ENDIF


      !--------------------------------------------------------------------------------
      !> Flags to perform method calculations.  Info from input file.
      !! Flag indicates method use in calculations (ends with F).
      !--------------------------------------------------------------------------------

         ! Mean drift method
      p%MnDriftF              =  InitInp%MnDriftF           ! Flag for calculation

         ! Newman approximation method
      p%NewmanAppF            =  InitInp%NewmanAppF         ! Flag for calculation

         ! Difference QTF
      p%DiffQTFF              =  InitInp%DiffQTFF           ! Flag for calculation

         ! Summation QTF
      p%SumQTFF               =  InitInp%SumQTFF            ! Flag for calculation



      !--------------------------------------------------------------------------------
      !> Make array for holding the resulting timeseries for the second order force.
      !--------------------------------------------------------------------------------

         ! Allocate array for the WaveExtcn2.
      ALLOCATE( p%WaveExctn2(0:InitInp%NStepWave,6*p%NBody), STAT=ErrStatTmp)
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array p%WaveExctn2 to store '// &
                              'the 2nd order force data.',  ErrStat,ErrMsg,'CheckInitInp')
      IF (ErrStat >= AbortErrLev ) RETURN

!
!      !--------------------------------------------------------------------------------
!      !> FAST channel output
!      !--------------------------------------------------------------------------------
!
!      p%NumOuts               =  InitInp%NumOuts
!      p%NumOutAll             =  InitInp%NumOutAll


   END SUBROUTINE CheckInitInput




   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine reads in the datafile containing 3D data (Omega, WaveDir1, WaveDir2), stores it to Data3D, and sets flags
   !! indicating how complete the data is.
   !!
   !! The datafile is first scanned to make sure it is of the correct form, then read in in its entirety into a temporary matrix.
   !! Then we find the number of unique frequencies and unique wave direction pairs that are in it, and create the necessary
   !! data structures to hold all of it.  Once this is done, we check for the infinite and zero frequency limits and add space
   !! for them as needed.
   !!
   !! This matrix holding the data is read one line at a time and is placed into the correct location within the dataset.  The mask
   !! array is then marked to indicate valid data is at that coordinate index.
   !!
   !! Since the data may not contain the zero frequency and infinite frequency, we will create those values.  For the zero frequency,
   !! the QTF must be zero.  For the infinite frequency, we use the value from the highest wave frequency (essentially flat response).
   !!
   !! At the end of all this, we check the data for completeness and set the flags accordingly.
   !!
   SUBROUTINE Read_DataFile3D( InitInp, Filename3D, Data3D, ErrStat, Errmsg )

      IMPLICIT NONE

         ! Passed variables.
      TYPE(WAMIT2_InitInputType),         INTENT(IN   )  :: InitInp              !< Input data for initialization routine
      CHARACTER(*),                       INTENT(IN   )  :: Filename3D        !< Name of the file containing the 3D data
      TYPE(W2_InitData3D_Type),           INTENT(INOUT)  :: Data3D            !< 3D QTF data
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           !< The error value
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            !< A message about the error.

         ! Local variables
      INTEGER(IntKi)                                     :: UnitDataFile      !< The unit number for the currently open file

      INTEGER(IntKi)                                     :: NumDataColumns    !< Number of data columns in the file
      INTEGER(IntKi)                                     :: NumDataLines      !< Total number of lines in the file (excluding first text line if there is one)
      INTEGER(IntKi)                                     :: NumDataLinesKeep  !< Total number of lines in the file that are to be kept (positive force component index)
      INTEGER(IntKi)                                     :: NumHeaderLines    !< Flag to indicate if the first line of the file is text

         ! Raw file data storage
      REAL(SiKi),       ALLOCATABLE                      :: RawData3D(:,:)    !< The raw data from the entirety of the input file -- after ignoring negative force components
      REAL(SiKi),       ALLOCATABLE                      :: RawData3DTmp(:,:) !< The raw data from the entirety of the input file -- as read in.
      REAL(SiKi)                                         :: HighFreq1         !< The highest frequency found.  Needed for setting the infinite frequency data.
      INTEGER(IntKi)                                     :: WvFreq1HiIdx      !< Index to the highest wave 1 frequency
      INTEGER(IntKi)                                     :: WvFreq1LoIdx      !< Index to the  lowest wave 1 frequency
      LOGICAL                                            :: HaveZeroFreq1     !< Indicates we have a zer frequency value


         ! File reading variables
      CHARACTER(MaxFileInfoLineLen)                      :: TextLine          !< One line of text read from the file
      INTEGER(IntKi)                                     :: LineLen           !< The length of the line read in
      REAL(SiKi),       ALLOCATABLE                      :: TmpRealArr(:)     !< Temporary real array
      REAL(SiKi),       ALLOCATABLE                      :: TmpDataRow(:)     !< Single row of data
      REAL(SiKi),       ALLOCATABLE                      :: TmpWvFreq1(:)     !< Temporary array to hold the wave frequencies read in.

         ! Temporary sparseness checking flag
      LOGICAL                                            :: TmpSparseFlag     !< Temporary flag for checking sparseness

         ! Generic counters
      INTEGER(IntKi)                                     :: I                 !< Generic counter
      INTEGER(IntKi)                                     :: J                 !< Generic counter
      INTEGER(IntKi)                                     :: K                 !< Generic counter
      INTEGER(IntKi)                                     :: L                 !< Generic counter
      INTEGER(IntKi)                                     :: TmpCoord(4)       !< Temporary index coords to the Data3D matrix (3D + force dimension)

            ! Error handling temporary variables
      INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary variable for the local error status
      CHARACTER(2048)                                    :: ErrMsgTmp         !< Temporary error message variable
      INTEGER(IntKi)                                     :: LclErrStat        !< Temporary error status.  Used locally to indicate when we have reached the end of the file.
      CHARACTER(*), PARAMETER                            :: RoutineName = 'Read_DataFile3D'


      !> ## Subroutine Contents

      !--------------------------------------------------------------------------------
      !> Initialize variables
      !--------------------------------------------------------------------------------

         ! Initialize error variables
      ErrStat     = ErrID_None
      ErrMsg      = ''
      HaveZeroFreq1 = .FALSE.             ! If we find a zero frequency, we will set this to true


      !--------------------------------------------------------------------------------
      !> ### Check data file for consistency
      !--------------------------------------------------------------------------------

         !------------------------------------------------------------------------------
         !> 3D data files are only used for cases of { MnDrift || NewmanApp } = { 7 || 8 || 9}\n
         !!     .7 files are available from WAMIT version 7 only --> Not checked at present\n
         !!     .8 files only contain information for dimensions 1, 2, and 6 (x, y, yaw) -- set in CheckInitInput\n
         !!     .9 files contain all dimensions
         !------------------------------------------------------------------------------

         ! Find a unit number to use
      CALL GetNewUnit(UnitDataFile,ErrStatTmp,ErrMsgTmp)
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         UnitDataFile = -1
         CALL CleanUp()
         RETURN
      ENDIF

         ! Open the file
      CALL OpenFInpFile(  UnitDataFile, TRIM(Filename3D), ErrStat, ErrMsg )  ! Open file containing mean drift information
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         UnitDataFile = -1
         CALL CleanUp()
         RETURN
      ENDIF


         ! Do an initial readthrough and find out the length of the file, if there is a header line, and the number of columns in the file.
      CALL GetFileLength( UnitDataFile, TRIM(Filename3D), NumDataColumns, NumDataLines, NumHeaderLines, ErrStatTmp, ErrMsgTmp)
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF


         ! Make sure we have 8 columns of data in the data file.
      IF ( NumDataColumns /= 8 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' The 2nd order WAMIT data file '//TRIM(Filename3D)//' has '//TRIM(Num2LStr(NumDataColumns))// &
                        ' columns instead of the 8 columns expected.', ErrStat, ErrMsg, RoutineName)
         CALL CleanUp()
         RETURN
      ENDIF


      !----------------------------------------------------------------------------------
      !> ### Read and store the raw data from the file.  Convert all wave periods (s) into
      !!       frequency (rad/s)
      !----------------------------------------------------------------------------------

         ! Allocate the temporary array for reading in one line
      CALL AllocAry( TmpDataRow, NumDataColumns, ' Array for holding one line of 4D data for 2nd order WAMIT files', ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)

         ! Allocate an array to hold the entirety of the raw contents of the file
      CALL AllocAry( RawData3DTmp, NumDataLines, NumDataColumns, ' Array for holding raw 3D data for 2nd order WAMIT files', ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Read and discard the header lines
      DO I=1,NumHeaderLines
         CALL ReadLine( UnitDataFile, '', TextLine, LineLen, LclErrStat )
      ENDDO

         ! Read in the data one line at a time and put it in RawData3D
      DO I=1,NumDataLines
         CALL ReadAry( UnitDataFile, TRIM(Filename3D), TmpDataRow, NumDataColumns, 'RawData3D('//TRIM(Num2LStr(I))//',:)', &
                     ' Line '//TRIM(Num2LStr(NumHeaderLines+I))//' of '//TRIM(Filename3D), &
                     ErrStatTmp, ErrMsgTmp )      ! Note, not echoing this to anything.
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF
         RawData3DTmp(I,:) = TmpDataRow
      ENDDO

      CLOSE( UnitDataFile )
      UnitDataFile = -1


         !> Before continuing, we need to figure out how many actual lines of data we will
         !! be keeping (lines where the force component index is positive).  The force component
         !! is in column 4 for 3D data read in.  First find out if there are negative force
         !! components, then count the number of data rows we are keeping.  Then copy data over.
      IF ( MINVAL(RawData3DTmp(:,4)) < 0_IntKi ) THEN             ! check the 4th element (force component)
         CALL SetErrStat( ErrID_Warn,' Negative load components found (moving reference frame). Ignoring', &
               ErrStat,ErrMsg,RoutineName)

            ! Count how many lines we are keeping.
         NumDataLinesKeep = 0_IntKi
         DO I=1,NumDataLines
            IF ( RawData3DTmp(I,4) > 0 ) THEN
               NumDataLinesKeep = NumDataLinesKeep + 1
            ENDIF
         ENDDO

            ! Allocate an array to hold the data from the file that we are keeping
         CALL AllocAry( RawData3D, NumDataLinesKeep, NumDataColumns, ' Array for holding raw 3D data for 2nd order WAMIT files', ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         ENDIF

            ! Now copy over the data that we are keeping. Ignore the rest
         NumDataLinesKeep = 0_IntKi
         DO I=1,NumDataLines
            IF ( RawData3DTmp(I,4) > 0 ) THEN
               NumDataLinesKeep = NumDataLinesKeep + 1
               RawData3D( NumDataLinesKeep, : ) = RawData3DTmp(I,:)
            ENDIF
         ENDDO

            ! We no longer need the raw data from the file.
         IF (ALLOCATED(RawData3DTmp))     DEALLOCATE(RawData3DTmp,STAT=ErrStatTmp)

         ! no negative force components, so just move the array.
      ELSE
         CALL MOVE_ALLOC( RawData3DTmp, RawData3D )
         NumDataLinesKeep = NumDataLines
      ENDIF


         !> Before proceeding further, we are going to change the wave period from measured
         !! in seconds to a frequency measurement in _rad/s_.  We will step through the data
         !! matrix and convert the first column.
         !!
         !! Note that wave periods that are negative correspond to the \f$ \omega=0 \f$, and
         !! wave periods with \f$ \tau=0 \f$ mean \f$ \omega=\infty \f$.  Normally these
         !! would not be calculated for sum and difference frequencies as there is little
         !! information of interest there.
         !!
         !! In the implimentation here, the infinity frequency covers everything from just above
         !! highest frequency given in the file and infinity.  This means that everything the
         !! the highest frequency is handled exactly the same.

      DO I=1,NumDataLinesKeep
         IF ( EqualRealNos(RawData3D(I,1), 0.0_SiKi) ) THEN
            ! Leave it alone.  We will have to fix it afterwards.
         ELSE IF ( RawData3D(I,1) < 0 ) THEN
            ! Leave it alone.  We will have to fix it afterwards.
         ELSE
            RawData3D(I,1) =  TwoPi_S/RawData3D(I,1)       ! First column is Tau1
         ENDIF
      ENDDO

         ! Now we can fix the value for the frequency that was negative and set it to sligtly
         ! larger than the largest value. First find the highest value.
      HighFreq1   = MAXVAL(RawData3D(:,1))

         ! First change the values given as 0.0 to the infinite frequency and set a flag to indicate we
         ! found an infinite frequency in this dimension.  We will later copy this value to the infinite
         ! frequency.  This results in having a flat response from the highest frequency found to the
         ! infinite frequency.
      DO I=1,NumDataLinesKeep
         IF ( EqualRealNos(RawData3D(I,1), 0.0_SiKi) ) THEN
            RawData3D(I,1) = HighFreq1 * OnePlusEps
         ENDIF
      ENDDO

         ! Now change the negative values to be 0.0 (zero frequency)
      DO I=1,NumDataLinesKeep
         IF ( RawData3D(I,1) < 0.0_SiKi ) THEN
            RawData3D(I,1) = 0.0_SiKi
            HaveZeroFreq1 = .TRUE.
         ENDIF
      ENDDO


      !----------------------------------------------------------------------------------
      !> ### Read through the file and find the number of WvFreq1, WaveDir1, WaveDir2.
      !!
      !! The 2nd order 3D files are arranged as follows:
      !!       Tau1, Beta1, Beta2, k, |F|, Phase, Real, Imaginary
      !!    where:
      !!    | Column |  Variable       |  Units          |  Description                      |
      !!    | :----: | :-------------: | :-------------: | :-------------------------------- |
      !!    |  1     |\f$ \omega_1 \f$ |  rad/s          |  Wave period 1 (converted above)  |
      !!    |  2     |\f$ \beta_1  \f$ |  degrees        |  Wave direction 1                 |
      !!    |  3     |\f$ \beta_2  \f$ |  degrees        |  Wave direction 2                 |
      !!    |  4     |   k             |  --             |  Force component direction (1-6)  |
      !!    |  5     |\f$ |F|      \f$ |  Nondimensional |  Magnitude      of the force QTF  |
      !!    |  6     |   Phase         |  degrees        |  Phase          of the force QTF  |
      !!    |  7     |\f$ \Re(F)   \f$ |  Nondimensional |  Real part      of the force QTF  |
      !!    |  8     |\f$ \Im(F)   \f$ |  Nondimensional |  Imaginary part of the force QTF  |
      !!
      !! Only columns 1, 2, 3, 4, 7, 8 are used.  Columns 5 & 6 are redundant information.
      !!
      !! Note that we will be checking for the zero frequency, and adding it if it was not
      !! found in the file.
      !----------------------------------------------------------------------------------


      !----------------------------------------------------------------------------------
      !> Read through the 3D data and figure out how many unique values in each dependent
      !! variable (\f$ \omega_1, \beta_1, \beta_2,\f$ Component direction) exist in the file.
      !!
      !! Allocate the necessary storage arrays when complete.
      !----------------------------------------------------------------------------------

         ! Get the number of first frequencies
      CALL UniqueRealValues( RawData3D(:,1), TmpWvFreq1, Data3D%NumWvFreq1, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Get the number of first wave directions
      CALL UniqueRealValues( RawData3D(:,2), Data3D%WvDir1, Data3D%NumWvDir1, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Get the number of second wave directions
      CALL UniqueRealValues( RawData3D(:,3), Data3D%WvDir2, Data3D%NumWvDir2, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF


         ! Find out which load components are actually in use
      CALL UniqueRealValues( RawData3D(:,4), TmpRealArr, K, ErrStatTmp,ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )

         ! Now figure out how many bodies there are. Each body will have up to 6 components.  The component
         ! column can be up to 6*N where N is the number of bodies in the file.  We will assume that we don't
         ! skip groups of bodies.
      Data3D%NumBodies = ceiling((maxval(TmpRealArr)-0.1_ReKi) / 6.0_ReKi)    ! Account for any uncertainty in the number
      IF ( Data3D%NumBodies < 1 ) CALL SetErrStat( ErrID_Fatal, ' No WAMIT bodies found (no positive load component numbers in column 4) in '// &
                     TRIM(Filename3D)//'.', ErrStat,ErrMsg,RoutineName )
      IF ( Data3D%NumBodies > 1 ) CALL SetErrStat( ErrID_Info, ' Found data for '//TRIM(Num2LStr(Data3D%NumBodies))//' WAMIT bodies in '// &
                     TRIM(Filename3D)//'.', ErrStat,ErrMsg,RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF


         ! Now that we know how many bodies are in the file, allocate the size of the arrays
      CALL AllocAry( Data3D%DataIsSparse,   6*Data3D%NumBodies, ' Array for tracking which dimension indices are sparsely populated', ErrStatTmp, ErrMsgTmp )
      CALL AllocAry( Data3D%LoadComponents, 6*Data3D%NumBodies, ' Array for tracking which dimension indices contain information',    ErrStatTmp, ErrMsgTmp )
      Data3D%DataIsSparse = .TRUE.        ! Assume the data is sparse, then change this after checking on the dimensions of interest.


         ! Now check the values we got back and set the LoadComponents flags for those with data. The
         ! load component direction must be between 1 and 6 (translational: 1,2,3; rotational: 4,5,6).
      Data3D%LoadComponents = .FALSE.
      DO I=1,K
         IF ( NINT(TmpRealArr(I)) < 1 .OR. NINT(TmpRealArr(K)) > 6*Data3D%NumBodies ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Load components listed in column 4 of '//TRIM(Filename3D)// &
                  ' must be between 1 and '//TRIM(Num2LStr(6*Data3D%NumBodies))//' for '//TRIM(Num2LStr(Data3D%NumBodies))// &
                  ' WAMIT bodies.', ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         ENDIF
         Data3D%LoadComponents(NINT(TmpRealArr(I))) = .TRUE.
      ENDDO


      IF (ALLOCATED(TmpRealArr))       DEALLOCATE(TmpRealArr,STAT=ErrStatTmp)     ! Done with this, so throw it away.
      IF (ALLOCATED(TmpDataRow))       DEALLOCATE(TmpDataRow,STAT=ErrStatTmp)




         ! Now we need to figure out if the zero frequency was given in the file.  If so, we change NumWvFreq1 to
         ! NumWvFreq1+2.  If not, change to NumWvFreq1+4.  We will add on the inifinite frequency value and
         ! zero out all values not in the input frequency range. The inifinite frequency value will be set to HUGE
         ! and we'll add/subtract epsilon to the first non-zero frequency entered so that we can achieve a step
         ! change for zeroing the values outside the input frequency range.
      IF (HaveZeroFreq1) THEN
         Data3D%NumWvFreq1 = Data3D%NumWvFreq1+2
         WvFreq1LoIdx   = 1
      ELSE
         Data3D%NumWvFreq1 = Data3D%NumWvFreq1+4
         WvFreq1LoIdx   = 3
      ENDIF

         ! Set the index for the highest frequency stored before the cutoff
      WvFreq1HiIdx   =  Data3D%NumWvFreq1-2

         ! Now allocate the array for holding the WvFreq1
      ALLOCATE( Data3D%WvFreq1( Data3D%NumWvFreq1), STAT=ErrStatTmp)
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array Data3D%WvFreq1 to store '// &
                              'the sorted 3D 2nd order WAMIT frequency data.',  ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Populate the wave frequencies with what we read in
      Data3D%WvFreq1( WvFreq1LoIdx:WvFreq1HiIdx ) = TmpWvFreq1

         ! If no zero frequency was supplied, add the two points for step-change before first entered frequency
      IF ( .NOT. HaveZeroFreq1) THEN
         Data3D%WvFreq1( 1 )                 = 0.0_SiKi
         Data3D%WvFreq1( 2 )                 = MAX( TmpWvFreq1(1) - 10.0_SiKi*EPSILON(0.0_SiKi), 0.0_SiKi )  ! make sure the Frequencies are still monotonically increasing
      ENDIF

         ! add the two points for step-change after last entered frequency
      Data3D%WvFreq1( Data3D%NumWvFreq1-1 )     = Data3D%WvFreq1(Data3D%NumWvFreq1-2) + 10.0_SiKi*EPSILON(0.0_SiKi)
      Data3D%WvFreq1( Data3D%NumWvFreq1   )     = HUGE(1.0_SiKi)/5 ! floating overflow occurs later with arithmetic so I divided by a small constant



         ! Now that we know how many frequencies and wave directions there are, we can allocate the array
         ! for saving the sorted data.
      ALLOCATE( Data3D%DataSet( Data3D%NumWvFreq1, Data3D%NumWvDir1, Data3D%NumWvDir2, 6*Data3D%NumBodies ),  STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array Data3D%DataSet to store '// &
                              'the sorted 3D 2nd order WAMIT data.',  ErrStat,ErrMsg,RoutineName)

         ! Allocate the logical array for storing the mask for which points are valid. Set to .FALSE.
      ALLOCATE( Data3D%DataMask( Data3D%NumWvFreq1, Data3D%NumWvDir1, Data3D%NumWvDir2, 6*Data3D%NumBodies ),  STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array Data3D%DataMask to store '// &
                              'the sorted 3D 2nd order WAMIT data.',  ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

      Data3D%DataMask = .FALSE.     ! Populate this later



      !----------------------------------------------------------------------------------
      !> Now step through the raw data and place values in the correct final locations,
      !! set the flag mask array, and calculate the conjugate pairs (4D only).
      !!
      !! The way this works is that the Data3D%WvFreq1, Data3D%WvDir1, and Data3D%WaveDir2
      !! contain the order frequencies and direction from the input file.  Therefore, the
      !! index numbers of these three arrays corresponds to the location in Data3D%DataSet
      !! that holds the complex force QTF value for that frequency and wave direction pair.
      !!
      !! So, to populate the Data3D%DataSet matrix, we simply have to read one line at a
      !! time from the file (stored in RawData3D matrix), and place it in the corresponding
      !! coordinate location in Data3D%DataSet.  This involves simply searching through
      !! the WvFreq1, WvDir1, and WvDir2 arrays for the correct indices.
      !!
      !! The wave force component direction is stored slightly differently, so it only
      !! needs to be read from the raw data and coverted into an integer to use as an
      !! index.
      !----------------------------------------------------------------------------------


      TmpCoord = 1         ! Initialize the search locations.

      DO I=1,NumDataLinesKeep

            ! Error checking: The LocateStp routine will return 0 if the requested value is less than the value
            !                 of the first element in the array.  It will return the index of the last element
            !                 of the array if the value is larger than the last element.  In creating the arrays
            !                 that are being searched here, the values that we are now trying to find the index
            !                 to were used.  Therefore, if the requested value is not found in the array, then
            !                 something must have gone horribly wrong while creating it, or between then and now,
            !                 which is most likely a programming error.

            ! Find the location in the WvFreq1 array that this point corresponds to.  We will check only between the
            ! cutoffs that were added to the frequency range.  This is contained within TmpWvFreq1 from reading in.
         CALL LocateStp( RawData3D(I,1), TmpWvFreq1, TmpCoord(1),   WvFreq1HiIdx - (WvFreq1LoIdx - 1) )  ! inclusive limits
         IF ( TmpCoord(1) == 0 .OR. ( RawData3D(I,1) > Data3D%WvFreq1(Data3D%NumWvFreq1)) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  Array data point not found in Data3D%WvFreq1 array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         ENDIF
         TmpCoord(1) =  TmpCoord(1) + ( WvFreq1LoIdx - 1 )     ! shift to the point in the Data3D%WvFreq1 array by adding the zero frequency step function

            ! Find the location in the WvDir1 array that this point corresponds to.
         CALL LocateStp( RawData3D(I,2), Data3D%WvDir1,  TmpCoord(2),   Data3D%NumWvDir1 )
         IF ( TmpCoord(2) == 0 .OR. ( RawData3D(I,2) > Data3D%WvDir1(Data3D%NumWvDir1)) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  Array data point not found in Data3D%WvDir1 array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         ENDIF

            ! Find the location in the WvDir2 array that this point corresponds to.
         CALL LocateStp( RawData3D(I,3), Data3D%WvDir2,  TmpCoord(3),   Data3D%NumWvDir2 )
         IF ( TmpCoord(3) == 0 .OR. ( RawData3D(I,3) > Data3D%WvDir2(Data3D%NumWvDir2)) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  Array data point not found in Data3D%WvDir2 array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         ENDIF

            ! Find which force component this belongs to
         TmpCoord(4) = NINT(RawData3D(I,4))
            ! Check that it is a valid force component
         if (TmpCoord(4) < 1 .or. TmpCoord(4) > 6*Data3D%NumBodies) then
            CALL SetErrStat( ErrID_Fatal, ' Line '//TRIM(Num2Lstr(NumHeaderLines+I))//' of '//TRIM(Filename3D)// &
                           ' contains force component '//TRIM(Num2LStr(TmpCoord(4)))//' which is outside the expected force '// &
                           ' range of 1 to '//TRIM(Num2Lstr(6*Data3D%NumBodies))//' for a '//TRIM(Num2LStr(Data3D%NumBodies))// &
                           ' body system.', ErrStat, ErrMsg, RoutineName)
            IF (ALLOCATED(RawData3D))        DEALLOCATE(RawData3D,STAT=ErrStatTmp)
            IF (ALLOCATED(RawData3DTmp))     DEALLOCATE(RawData3DTmp,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpRealArr))       DEALLOCATE(TmpRealArr,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpDataRow))       DEALLOCATE(TmpDataRow,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpWvFreq1))       DEALLOCATE(TmpWvFreq1,STAT=ErrStatTmp)
            CALL CleanUp
            RETURN
         endif



            !> The data from the WAMIT file is non-dimensional, so we need to dimensionalize it here.  This
            !! is a partial dimensionalization since the wave amplitudes are not included (this is done later
            !! in each of the calculation methods).  To dimensionalize the data, the equation is for the
            !! partially dimensionalized force (\f$ F_k \f$) is:
            !!
            !!       \f$ F_k = \rho g \cdot L^\alpha \cdot  \bar{F}_k \f$
            !!
            !! where \f$ \bar{F}_k \f$ is the force QTF value in the file for the \f$ k \f$ component
            !! direction, \f$ L \f$ is the WAMIT unit length _WAMITULEN_, \f$ \rho \f$ is the density of
            !! water, and \f$ g \f$ is the gravitational constant (\f$ \rho g \f$ is combined as _RhoXg_).
            !! The value of \f$ \alpha \f$ is 1 for \f$ k = 1,2,3 \f$ and 2 for \f$ k = 4,5,6 \f$.

            ! Here K is used for \f$ alpha \f$ in the dimensionalization equation.
         IF ( TmpCoord(4) <=3 ) THEN
            K = 1
         ELSE
            K = 2
         ENDIF


            ! Check that the current value has not been read in already.  If it has, the corresponding
            ! flag in the mask array will be set to true.  This will indicate that there was redundant
            ! data in the file.  If the value agrees, we will politely ignore it without warning.  If
            ! it does not agree, then we will have to stop since it is ambiguous which is correct.

         IF ( Data3D%DataMask( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4) ) ) THEN
            IF ( .NOT. EqualRealNos(REAL(Data3D%DataSet( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4) ),SiKi), &
                                    REAL(InitInp%RhoXg * InitInp%WAMITULEN**K * RawData3D(I,7)               ,SiKi)) .AND. &
                 .NOT. EqualRealNos( AIMAG(Data3D%DataSet( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4) ) ), &
                                    REAL(InitInp%RhoXg * InitInp%WAMITULEN**K * RawData3D(I,8)               ,SiKi)) ) THEN
               CALL SetErrStat( ErrID_Fatal, ' Line '//TRIM(Num2Lstr(NumHeaderLines+I))//' of '//TRIM(Filename3D)// &
                        ' contains different values for the real and imaginary part (columns 7 and 8) than was '// &
                        'given earlier in the file for the same values of wave frequency and wave direction '// &
                        '(force dimension = '//TRIM(Num2LStr(TmpCoord(4)))//').', &
                        ErrStat, ErrMsg, RoutineName )
               CALL CleanUp()
               RETURN
            ENDIF
         ELSE

               ! Store the data after dimensionalizing
            Data3D%DataSet( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4) ) = &
                      REAL(InitInp%RhoXg * InitInp%WAMITULEN**K,SiKi) * CMPLX(RawData3D(I,7),RawData3D(I,8),SiKi)

               ! Set flag indicating that this value has been inserted.
            Data3D%DataMask( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4) ) = .TRUE.

         ENDIF

      ENDDO


         !> Note that in the 3D WAMIT output files (.7, .8, .9 files), only the diagonal terms where
         !! \f$ \omega_m = \omega_m \f$ are given. If multidirectional waves are used, both WvDir1
         !! and WvDir2 will need to be populated.  According to the WAMIT theory guide (page 4-7),
         !! it there is a relationship between the wave direction pairs for this case:
         !!       \f$ \bar{F}(\omega_n,\omega_m) = \bar{F}^*(\omega_m,\omega_n) \f$.
         !! So, we can use this to populate missing values of WvDir2,WvDir1 from the WvDir1,WvDir2
         !! pair.
         !!
         !! Since only one wave frequency array is present, WvFreq1, this is possible to do.  For
         !! two wave frequencies, such as in a full QTF, the calculation performed here would need
         !! to be performed only when \f$ \omega_1 = \omega_2 \f$.

         ! Loop over the wave components, but only perform calculations on the ones that have values
      DO L=1,6*Data3D%NumBodies
         IF (Data3D%LoadComponents(L)) THEN        ! Only do this for the load components that exist
            DO I=1,Data3D%NumWvFreq1               ! Frequencies

                  ! Now look through the wave direction matrix (all combinations of (beta_1,beta_2) )
                  ! and if it is missing, set it equal to the complex conjugate of the
                  ! (beta_2,beta_1) element if it exists.  We don't need to check the beta_1 = beta_2
                  ! elements.
               DO J=1,Data3D%NumWvDir1
                  DO K=J,Data3D%NumWvDir2

                        ! Only non-diagonal elements
                     IF ( .NOT. EqualRealNos( Data3D%WvDir1(J), Data3D%WvDir2(K) ) ) THEN

                           ! See if WvDir1(J) == WvDir2(J) and WvDir1(K) == WvDir2(K), because
                           ! if they are not, then the discretization along the two wave directions
                           ! is different and this won't work.
                        IF ( EqualRealNos( Data3D%WvDir1(J), Data3D%WvDir2(J) ) .AND. &
                             EqualRealNos( Data3D%WvDir1(K), Data3D%WvDir2(K) ) ) THEN

                           ! Check if filled
                           IF ( Data3D%DataMask( I, J, K, L )  ) THEN

                                 ! See if the diagonal mirror one (WvDir2,WvDir1) value is not filled,
                                 ! and fill it if empty
                              IF ( .NOT. Data3D%DataMask( I, K, J, L )  ) THEN
                                 Data3D%DataSet ( I, K, J, L ) = Data3D%DataSet( I, J, K, L )
                                 Data3D%DataMask( I, K, J, L ) = .TRUE.
                              ENDIF
                           ENDIF
                        ENDIF ! Check that wave directions will pair.
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO       ! Checking the wave directions for completeness.


      !----------------------------------------------------------------------------------
      !> We added two frequencies for the \f$ omega = 0 \f$ term if it did not exist,
      !! and added two frequencies for the infinite frequency term on the end of the array,
      !! to create step changes outside the entered frequency ranges. We need to populate
      !! the these new terms (set to zero).
      !----------------------------------------------------------------------------------

      IF (.NOT. HaveZeroFreq1) THEN
         Data3D%DataSet( 1:2,:,:,:)  = CMPLX(0.0_SiKi,0.0_SiKi)                                     ! Set the values to zero for everything before entered frequency range
         Data3D%DataMask(1:2,:,:,:)  = .TRUE.                                                       ! Set the mask for these first two frequencies
      ENDIF
      Data3D%DataSet( Data3D%NumWvFreq1-1:Data3D%NumWvFreq1,:,:,:) = CMPLX(0.0_SiKi,0.0_SiKi)       ! Set the values for the last two frequencies to zero (everything higher than the last non-infinite frequency)
      Data3D%DataMask(Data3D%NumWvFreq1-1:Data3D%NumWvFreq1,:,:,:) = .TRUE.                         ! Set the mask for the last two frequencies


      !----------------------------------------------------------------------------------
      !> Find out if the data is sparse or full.  Verification that the requested component
      !! directions were found will occur in the calling routine, not here (that information
      !! was not passed in).  We will check this by sweeping through all three dimensions
      !! for only the values of k that have data in them.
      !----------------------------------------------------------------------------------

      DO L=1,6*Data3D%NumBodies       ! Loop over force component directions
         TmpSparseFlag  = .FALSE.                  ! Change this to true if any empty values are found
         IF (Data3D%LoadComponents(L)) THEN        ! Only if we found data for that component
            DO I=1,Data3D%NumWvFreq1
               DO J=1,Data3D%NumWvDir1
                  DO K=1,Data3D%NumWvDir2
                     IF (.NOT. Data3D%DataMask( I, J, K, L ) ) THEN
                       TmpSparseFlag = .TRUE.
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

            ! If any values were missing for this force component direction, TmpSparseFlag will be true.
         Data3D%DataIsSparse(L)  = TmpSparseFlag

      ENDDO    ! Sparseness check

      !> In the MnDrift and NewmanApp calculations used in this module, the wave directions are
      !! assigned as one direction per frequency.  As a result, we will never have a situation
      !! with these methods where we need two wave directions since they involve only a single
      !! frequency.  So, if multidirectional waves are used, both WvDir1 and WvDir2 will need to
      !! be populated.  According to the WAMIT theory guide (page 4-7), it is possible to calculate
      !! the value of \f$ \bar{F}(\omega_n,\omega_m) = \bar{F}^*(\omega_m,\omega_n) \f$.



         ! Clean up
      call cleanup()
      
   CONTAINS
      subroutine cleanup()
      
         if (UnitDataFile > 0) CLOSE( UnitDataFile )
         
         IF (ALLOCATED(RawData3D))        DEALLOCATE(RawData3D,STAT=ErrStatTmp)
         IF (ALLOCATED(RawData3DTmp))     DEALLOCATE(RawData3DTmp,STAT=ErrStatTmp)
         IF (ALLOCATED(TmpRealArr))       DEALLOCATE(TmpRealArr,STAT=ErrStatTmp)
         IF (ALLOCATED(TmpDataRow))       DEALLOCATE(TmpDataRow,STAT=ErrStatTmp)
         IF (ALLOCATED(TmpWvFreq1))       DEALLOCATE(TmpWvFreq1,STAT=ErrStatTmp)
      end subroutine cleanup
   END SUBROUTINE Read_DataFile3D



   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine reads in the datafile containing 4D data (Omega1, Omega2, WaveDir1, WaveDir2), stores it to Data4D, and sets flags
   !! indicating how complete the data is.
   !!
   !! The datafile is first scanned to make sure it is of the correct form, then read in in its entirety into a temporary matrix.
   !! Then we find the number of unique frequencies and unique wave direction pairs that are in it, and create the necessary
   !! data structures to hold all of it.  Once this is done, we check for the infinite and zero frequency limits and add space
   !! for them as needed.
   !!
   !! This matrix holding the data is read one line at a time and is placed into the correct location within the dataset.  The mask
   !! array is then marked to indicate valid data is at that coordinate index.
   !!
   !! Since the data may not contain the zero frequency and infinite frequency, we will create those values.  For the zero frequency,
   !! the QTF must be zero.  For the infinite frequency, we use the value from the highest wave frequency (essentially flat response).
   !!
   !! At the end of all this, we check the data for completeness and set the flags accordingly.
   !!
   SUBROUTINE Read_DataFile4D( InitInp, Filename4D, Data4D, ErrStat, Errmsg )

      IMPLICIT NONE

         ! Passed variables.
      TYPE(WAMIT2_InitInputType),         INTENT(IN   )  :: InitInp              !< Input data for initialization routine
      CHARACTER(*),                       INTENT(IN   )  :: Filename4D        !< Name of the file containing the 4D data
      TYPE(W2_InitData4D_Type),           INTENT(INOUT)  :: Data4D            !< 4D QTF data
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           !< The error value
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            !< A message about the error.

         ! Local variables
      INTEGER(IntKi)                                     :: UnitDataFile      !< The unit number for the currently open file

      INTEGER(IntKi)                                     :: NumDataColumns    !< Number of data columns in the file
      INTEGER(IntKi)                                     :: NumDataLines      !< Total number of lines in the file (excluding first text line if there is one)
      INTEGER(IntKi)                                     :: NumDataLinesKeep  !< Total number of lines in the file that are to be kept (positive force component index)
      INTEGER(IntKi)                                     :: NumHeaderLines    !< Flag to indicate if the first line of the file is text

         ! Raw file data storage
      REAL(SiKi),       ALLOCATABLE                      :: RawData4D(:,:)    !< The raw data from the entirety of the input file -- after removing negative force components
      REAL(SiKi),       ALLOCATABLE                      :: RawData4DTmp(:,:) !< The raw data from the entirety of the input file -- as read in
      REAL(SiKi)                                         :: HighFreq1         !< The highest frequency found.  Needed for setting the infinite frequency data.
      REAL(SiKi)                                         :: HighFreq2         !< The highest frequency found.  Needed for setting the infinite frequency data.
      INTEGER(IntKi)                                     :: WvFreq1HiIdx      !< Index to the highest wave 1 frequency
      INTEGER(IntKi)                                     :: WvFreq1LoIdx      !< Index to the  lowest wave 1 frequency
      INTEGER(IntKi)                                     :: WvFreq2HiIdx      !< Index to the highest wave 2 frequency
      INTEGER(IntKi)                                     :: WvFreq2LoIdx      !< Index to the  lowest wave 2 frequency
      LOGICAL                                            :: HaveZeroFreq1     !< Indicates we have a zer frequency value
      LOGICAL                                            :: HaveZeroFreq2     !< Indicates we have a zer frequency value


         ! File reading variables
      CHARACTER(MaxFileInfoLineLen)                      :: TextLine          !< One line of text read from the file
      INTEGER(IntKi)                                     :: LineLen           !< The length of the line read in
      REAL(SiKi),       ALLOCATABLE                      :: TmpRealArr(:)     !< Temporary real array
      REAL(SiKi),       ALLOCATABLE                      :: TmpDataRow(:)     !< Single row of data
      REAL(SiKi),       ALLOCATABLE                      :: TmpWvFreq1(:)     !< Temporary array to hold the wave frequencies read in.
      REAL(SiKi),       ALLOCATABLE                      :: TmpWvFreq2(:)     !< Temporary array to hold the wave frequencies read in.

         ! Temporary sparseness checking flag
      LOGICAL                                            :: TmpSparseFlag     !< Temporary flag for checking sparseness
      LOGICAL                                            :: TmpDiagComplete   !< Temporary flag for checking sparseness

         ! Generic counters
      INTEGER(IntKi)                                     :: I                 !< Generic counter
      INTEGER(IntKi)                                     :: J                 !< Generic counter
      INTEGER(IntKi)                                     :: K                 !< Generic counter
      INTEGER(IntKi)                                     :: L                 !< Generic counter
      INTEGER(IntKi)                                     :: M                 !< Generic counter
      INTEGER(IntKi)                                     :: TmpCoord(5)       !< Temporary index coords to the Data4D matrix (4D + force dimension)

            ! Error handling temporary variables
      INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary variable for the local error status
      CHARACTER(2048)                                    :: ErrMsgTmp         !< Temporary error message variable
      INTEGER(IntKi)                                     :: LclErrStat        !< Temporary error status.  Used locally to indicate when we have reached the end of the file.
      CHARACTER(*), PARAMETER                            :: RoutineName = 'Read_DataFile4D'


      !> ## Subroutine Contents

      !--------------------------------------------------------------------------------
      !> Initialize variables
      !--------------------------------------------------------------------------------

         ! Initialize error variables
      ErrStat     = ErrID_None
      ErrMsg      = ''
      HaveZeroFreq1 = .FALSE.             ! If we find a zero frequency, we will set this to true
      HaveZeroFreq2 = .FALSE.             ! If we find a zero frequency, we will set this to true
      UnitDataFile = -1

      !--------------------------------------------------------------------------------
      !> ### Check data file for consistency
      !--------------------------------------------------------------------------------

         !------------------------------------------------------------------------------
         !> 4D data files are only used for cases of { SumQTF || DiffQTF } = { 10 || 11 || 12}\n
         !!     .10 files can contain all dimensions \n
         !!     .11 files can contain all dimensions \n
         !!     .12 files can contain all dimensions
         !------------------------------------------------------------------------------

         ! Find a unit number to use
      CALL GetNewUnit(UnitDataFile,ErrStatTmp,ErrMsgTmp)
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Open the file
      CALL OpenFInpFile(  UnitDataFile, TRIM(Filename4D), ErrStatTmp, ErrMsgTmp )  ! Open file containing mean drift information
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE( UnitDataFile )
         CALL CleanUp()
         RETURN
      ENDIF

         ! Do an initial readthrough and find out the length of the file, if there is a header line, and the number of columns in the file.
      CALL GetFileLength( UnitDataFile, TRIM(Filename4D), NumDataColumns, NumDataLines, NumHeaderLines, ErrStatTmp, ErrMsgTmp)
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE( UnitDataFile )
         CALL CleanUp()
         RETURN
      ENDIF

      REWIND( UnitDataFile )



         ! Make sure we have 9 columns of data in the data file.
      IF ( NumDataColumns /= 9 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' The 2nd order WAMIT data file '//TRIM(Filename4D)//' has '//TRIM(Num2LStr(NumDataColumns))// &
                        ' columns instead of the 9 columns expected.', ErrStat, ErrMsg, RoutineName)
         CLOSE( UnitDataFile )
         CALL CleanUp()
         RETURN
      ENDIF



      !----------------------------------------------------------------------------------
      !> ### Read and store the raw data from the file.  Convert all wave periods (s) into
      !!       frequency (rad/s)
      !----------------------------------------------------------------------------------

         ! Allocate the temporary array for reading in one line
      CALL AllocAry( TmpDataRow, NumDataColumns, ' Array for holding one line of 4D data for 2nd order WAMIT files', ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)

         ! Allocate an array to hold the entirety of the raw contents of the file
      CALL AllocAry( RawData4DTmp, NumDataLines, NumDataColumns, ' Array for holding raw 4D data for 2nd order WAMIT files', ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CLOSE( UnitDataFile )
         CALL CleanUp()
         RETURN
      ENDIF

         ! Read and discard the header lines
      DO I=1,NumHeaderLines
         CALL ReadLine( UnitDataFile, '', TextLine, LineLen, LclErrStat )
      ENDDO

         ! Read in the data one line at a time and put it in RawData4D
      DO I=1,NumDataLines
         CALL ReadAry( UnitDataFile, TRIM(Filename4D), TmpDataRow, NumDataColumns, 'RawData4DTmp('//TRIM(Num2LStr(I))//',:)', &
                     ' Line '//TRIM(Num2LStr(NumHeaderLines+I))//' of '//TRIM(Filename4D), &
                     ErrStatTmp, ErrMsgTmp )      ! Note, not echoing this to anything.
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnitDataFile )
            CALL CleanUp()
            RETURN
         ENDIF
         RawData4DTmp(I,:) = TmpDataRow
      ENDDO

      CLOSE( UnitDataFile )


         !> Before continuing, we need to figure out how many actual lines of data we will
         !! be keeping (lines where the force component index is positive).  The force component
         !! is in column 5 for 4D data read in.  First find out if there are negative force
         !! components, then count the number of data rows we are keeping.  Then copy data over.
      IF ( MINVAL(RawData4DTmp(:,5)) < 0_IntKi ) THEN          ! check the 5th element (force component)
         CALL SetErrStat( ErrID_Warn,' Negative load components found (moving reference frame). Ignoring', &
               ErrStat,ErrMsg,RoutineName)

            ! Count how many lines we are keeping.
         NumDataLinesKeep = 0_IntKi
         DO I=1,NumDataLines
            IF ( RawData4DTmp(I,5) > 0 ) THEN
               NumDataLinesKeep = NumDataLinesKeep + 1
            ENDIF
         ENDDO

            ! Allocate an array to hold the data from the file that we are keeping
         CALL AllocAry( RawData4D, NumDataLinesKeep, NumDataColumns, ' Array for holding raw 4D data for 2nd order WAMIT files', ErrStatTmp, ErrMsgTmp )
         CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnitDataFile )
            CALL CleanUp()
            RETURN
         ENDIF

            ! Now copy over the data that we are keeping. Ignore the rest
         NumDataLinesKeep = 0_IntKi
         DO I=1,NumDataLines
            IF ( RawData4DTmp(I,5) > 0 ) THEN
               NumDataLinesKeep = NumDataLinesKeep + 1
               RawData4D( NumDataLinesKeep, : ) = RawData4DTmp(I,:)
            ENDIF
         ENDDO

            ! We no longer need the raw data from the file.
         IF (ALLOCATED(RawData4DTmp))     DEALLOCATE(RawData4DTmp,STAT=ErrStatTmp)

         ! no negative force components, so just move the array.
      ELSE
         CALL MOVE_ALLOC( RawData4DTmp, RawData4D )
         NumDataLinesKeep = NumDataLines
      ENDIF


         !> Before proceeding further, we are going to change the wave period from measured
         !! in seconds to a frequency measurement in _rad/s_.  We will step through the data
         !! matrix and convert the first two columns.
         !!
         !! Note that wave periods that are negative correspond to the \f$ \omega=0 \f$, and
         !! wave periods with \f$ \tau=0 \f$ mean \f$ \omega=\infty \f$.  Normally these
         !! would not be calculated for sum and difference frequencies as there is little
         !! information of interest there.
      DO I=1,NumDataLinesKeep
         IF ( EqualRealNos(RawData4D(I,1), 0.0_SiKi) ) THEN
            ! Leave it alone.  We will have to fix it afterwards.
         ELSE IF ( RawData4D(I,1) < 0 ) THEN
            ! Leave it alone.  We will have to fix it afterwards.
         ELSE
            RawData4D(I,1) =  TwoPi_S/RawData4D(I,1)       ! First column is Tau1
         ENDIF
         IF ( EqualRealNos(RawData4D(I,2), 0.0_SiKi) ) THEN
            ! Leave it alone.  We will have to fix it afterwards.
         ELSE IF ( RawData4D(I,2) < 0 ) THEN
            ! Leave it alone.  We will have to fix it afterwards.
         ELSE
            RawData4D(I,2) =  TwoPi_S/RawData4D(I,2)       ! First column is Tau2
         ENDIF
      ENDDO

         ! Now we can fix the value for the frequency that was negative and set it to sligtly
         ! larger than the largest value. First find the highest value.
      HighFreq1   = MAXVAL(RawData4D(:,1))
      HighFreq2   = MAXVAL(RawData4D(:,2))

         ! First change the values given as 0.0 to the infinite frequency and set a flag to indicate we
         ! found an infinite frequency in this dimension
      DO I=1,NumDataLinesKeep
         IF ( EqualRealNos(RawData4D(I,1), 0.0_SiKi) ) THEN
            RawData4D(I,1) = HighFreq1 * OnePlusEps
         ENDIF
         IF ( EqualRealNos(RawData4D(I,1), 0.0_SiKi) ) THEN
            RawData4D(I,1) = HighFreq2 * OnePlusEps
         ENDIF
      ENDDO

         ! Now change the negative values to be 0.0
      DO I=1,NumDataLinesKeep
         IF ( RawData4D(I,1) < 0.0_SiKi ) THEN
            RawData4D(I,1) = 0.0_SiKi
            HaveZeroFreq1 = .TRUE.
         ENDIF
         IF ( RawData4D(I,2) < 0.0_SiKi ) THEN
            RawData4D(I,2) = 0.0_SiKi
            HaveZeroFreq2 = .TRUE.
         ENDIF
      ENDDO



      !-------------------------------------------------------------------------------------------
      !> ### Read through the file and find the number of WvFreq1, WaveDir1, WaveDir2.
      !!
      !! The 2nd order 4D files are arranged as follows:
      !!       Tau1, Beta1, Beta2, k, |F|, Phase, Real, Imaginary
      !!    where:
      !!    | Column |  Variable       |  Units          |  Description                      |
      !!    | :----: | :-------------: | :-------------: | :-------------------------------- |
      !!    |  1     |\f$ \omega_1 \f$ |  rad/s          |  Wave period 1 (converted above)  |
      !!    |  2     |\f$ \omega_2 \f$ |  rad/s          |  Wave period 2 (converted above)  |
      !!    |  3     |\f$ \beta_1  \f$ |  degrees        |  Wave direction 1                 |
      !!    |  4     |\f$ \beta_2  \f$ |  degrees        |  Wave direction 2                 |
      !!    |  5     |   k             |  --             |  Force component direction (1-6)  |
      !!    |  6     |\f$ |F|      \f$ |  Nondimensional |  Magnitude      of the force QTF  |
      !!    |  7     |   Phase         |  degrees        |  Phase          of the force QTF  |
      !!    |  8     |\f$ \Re(F)   \f$ |  Nondimensional |  Real part      of the force QTF  |
      !!    |  9     |\f$ \Im(F)   \f$ |  Nondimensional |  Imaginary part of the force QTF  |
      !!
      !! Only columns 1, 2, 3, 4, 7, 8 are used.  Columns 5 & 6 are redundant information.
      !!
      !----------------------------------------------------------------------------------


      !----------------------------------------------------------------------------------
      !> Read through the 4D data and figure out how many unique values in each dependent
      !! variable (\f$ \omega_1, \omega_2, \beta_1, \beta_2,\f$ Component direction) exist
      !! in the file.
      !!
      !! Allocate the necessary storage arrays when complete.
      !----------------------------------------------------------------------------------

         ! Get the number of first frequencies
      CALL UniqueRealValues( RawData4D(:,1), TmpWvFreq1, Data4D%NumWvFreq1, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Get the number of second frequencies
      CALL UniqueRealValues( RawData4D(:,2), TmpWvFreq2, Data4D%NumWvFreq2, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Get the number of first wave directions
      CALL UniqueRealValues( RawData4D(:,3), Data4D%WvDir1, Data4D%NumWvDir1, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Get the number of second wave directions
      CALL UniqueRealValues( RawData4D(:,4), Data4D%WvDir2, Data4D%NumWvDir2, ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF


         ! Find out which load components are actually in use
      CALL UniqueRealValues( RawData4D(:,5), TmpRealArr, K, ErrStatTmp,ErrMsgTmp )

         ! Now figure out how many bodies there are. Each body will have up to 6 components.  The component
         ! column can be up to 6*N where N is the number of bodies in the file.  We will assume that we don't
         ! skip groups of bodies.
      Data4D%NumBodies = ceiling((maxval(TmpRealArr)-0.1_ReKi) / 6.0_ReKi)    ! Account for any uncertainty in the number
      IF ( Data4D%NumBodies < 1 ) CALL SetErrStat( ErrID_Fatal, ' No WAMIT bodies found (no positive load component numbers in column 4) in '// &
                     TRIM(Filename4D)//'.', ErrStat,ErrMsg,RoutineName )
      IF ( Data4D%NumBodies > 1 ) CALL SetErrStat( ErrID_Info, ' Found data for '//TRIM(Num2LStr(Data4D%NumBodies))//' WAMIT bodies in '// &
                     TRIM(Filename4D)//'.', ErrStat,ErrMsg,RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF


         ! Now that we know how many bodies are in the file, allocate the size of the arrays
      CALL AllocAry( Data4D%DataIsSparse,   6*Data4D%NumBodies, ' Array for tracking which dimension indices are sparsely populated', ErrStatTmp, ErrMsgTmp )
      CALL AllocAry( Data4D%LoadComponents, 6*Data4D%NumBodies, ' Array for tracking which dimension indices contain information',    ErrStatTmp, ErrMsgTmp )
      Data4D%DataIsSparse = .TRUE.        ! Assume the data is sparse, then change this after checking on the dimensions of interest.

     CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )


         ! Now check the values we got back and set the LoadComponents flags for those with data. The
         ! load component direction must be between 1 and 6 (translational: 1,2,3; rotational: 4,5,6).
      Data4D%LoadComponents = .FALSE.
      DO I=1,K
         IF ( NINT(TmpRealArr(I)) < 1 .OR. NINT(TmpRealArr(K)) > 6*Data4D%NumBodies ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Load components listed in column 4 of '//TRIM(Filename4D)// &
                  ' must be between 1 and '//TRIM(Num2LStr(6*Data4D%NumBodies))//' for '//TRIM(Num2LStr(Data4D%NumBodies))// &
                  ' WAMIT bodies.', ErrStat,ErrMsg,RoutineName)
            CALL CleanUp()
            RETURN
         ENDIF
         Data4D%LoadComponents(NINT(TmpRealArr(I))) = .TRUE.
      ENDDO


      IF (ALLOCATED(TmpRealArr))       DEALLOCATE(TmpRealArr,STAT=ErrStatTmp)     ! Done with this, so throw it away.
      IF (ALLOCATED(TmpDataRow))       DEALLOCATE(TmpDataRow,STAT=ErrStatTmp)


         ! Now we need to figure out if the zero frequency was given in the file.  If so, we change NumWvFreq1 to
         ! NumWvFreq1+2.  If not, change to NumWvFreq1+4.  We will add on the inifinite frequency value and
         ! zero out all values not in the input frequency range. The inifinite frequency value will be set to HUGE
         ! and we'll add/subtract epsilon to the first non-zero frequency entered so that we can achieve a step
         ! change for zeroing the values outside the input frequency range.
      IF (HaveZeroFreq1) THEN
         Data4D%NumWvFreq1 = Data4D%NumWvFreq1+2
         WvFreq1LoIdx   = 1
      ELSE
         Data4D%NumWvFreq1 = Data4D%NumWvFreq1+4
         WvFreq1LoIdx   = 3
      ENDIF

         ! Set the index for the highest frequency stored before the cutoff
      WvFreq1HiIdx   =  Data4D%NumWvFreq1-2


         ! Do the same for NumWvFreq2 as we did for NumWvFreq2
      IF (HaveZeroFreq2) THEN
         Data4D%NumWvFreq2 = Data4D%NumWvFreq2+2
         WvFreq2LoIdx   = 1
      ELSE
         Data4D%NumWvFreq2 = Data4D%NumWvFreq2+4
         WvFreq2LoIdx   = 3
      ENDIF

         ! Set the index for the highest frequency stored before the cutoff
      WvFreq2HiIdx   =  Data4D%NumWvFreq2-2


         ! Now allocate the array for holding the WvFreq1
      ALLOCATE( Data4D%WvFreq1( Data4D%NumWvFreq1), STAT=ErrStatTmp)
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array Data4D%WvFreq1 to store '// &
                              'the sorted 4D 2nd order WAMIT frequency data.',  ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF


         ! Now allocate the array for holding the WvFreq2
      ALLOCATE( Data4D%WvFreq2( Data4D%NumWvFreq2), STAT=ErrStatTmp)
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array Data4D%WvFreq2 to store '// &
                              'the sorted 4D 2nd order WAMIT frequency data.',  ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Populate the wave frequencies with what we read in
      Data4D%WvFreq1( WvFreq1LoIdx:WvFreq1HiIdx ) = TmpWvFreq1

         ! If no zero frequency was supplied, add the two points for step-change before first entered frequency
      IF ( .NOT. HaveZeroFreq1) THEN
         Data4D%WvFreq1( 1 )                 = 0.0_SiKi
         Data4D%WvFreq1( 2 )                 = MAX( TmpWvFreq1(1) - 10.0_SiKi*EPSILON(0.0_SiKi), 0.0_SiKi )  ! make sure the Frequencies are still monotonically increasing
      ENDIF

         ! add the two points for step-change after last entered frequency
      Data4D%WvFreq1( Data4D%NumWvFreq1-1 )     = Data4D%WvFreq1(Data4D%NumWvFreq1-2) + 10.0_SiKi*EPSILON(0.0_SiKi)
      Data4D%WvFreq1( Data4D%NumWvFreq1   )     = HUGE(1.0_SiKi)/5 ! floating overflow occurs later with arithmetic so I divided by a small constant



         ! Populate the wave frequencies with what we read in
      Data4D%WvFreq2( WvFreq2LoIdx:WvFreq2HiIdx ) = TmpWvFreq2

         ! If no zero frequency was supplied, add the two points for step-change before first entered frequency
      IF ( .NOT. HaveZeroFreq2) THEN
         Data4D%WvFreq2( 1 )                 = 0.0_SiKi
         Data4D%WvFreq2( 2 )                 = MAX( TmpWvFreq2(1) - 10.0_SiKi*EPSILON(0.0_SiKi), 0.0_SiKi )  ! make sure the Frequencies are still monotonically increasing
      ENDIF

         ! add the two points for step-change after last entered frequency
      Data4D%WvFreq2( Data4D%NumWvFreq2-1 )     = Data4D%WvFreq2(Data4D%NumWvFreq2-2) + 10.0_SiKi*EPSILON(0.0_SiKi)
      Data4D%WvFreq2( Data4D%NumWvFreq2   )     = HUGE(1.0_SiKi)/5 ! floating overflow occurs later with arithmetic so I divided by a small constant


         ! Now that we know how many frequencies and wave directions there are, we can allocate the array
         ! for saving the sorted data.
      ALLOCATE( Data4D%DataSet( Data4D%NumWvFreq1, Data4D%NumWvFreq2, Data4D%NumWvDir1, Data4D%NumWvDir2, 6*Data4D%NumBodies ),  STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array Data4D%DataSet to store '// &
                              'the sorted 4D 2nd order WAMIT data.',  ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Allocate the logical array for storing the mask for which points are valid. Set to .FALSE.
      ALLOCATE( Data4D%DataMask( Data4D%NumWvFreq1, Data4D%NumWvFreq2, Data4D%NumWvDir1, Data4D%NumWvDir2, 6*Data4D%NumBodies ),  STAT=ErrStatTmp )
      IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array Data4D%DataMask to store '// &
                              'the sorted 4D 2nd order WAMIT data.',  ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

      Data4D%DataMask = .FALSE.     ! Populate this later



      !----------------------------------------------------------------------------------
      !> Now step through the raw data and place values in the correct final locations,
      !! set the flag mask array, and calculate the conjugate pairs (4D only).
      !!
      !! The way this works is that the Data4D%WvFreq1,  Data4D%WvFreq2, Data4D%WvDir1, and
      !! Data4D%WaveDir2 contain the order frequencies and direction from the input file.
      !! Therefore, the index numbers of these three arrays corresponds to the location in
      !! Data4D%DataSet that holds the complex force QTF value for that frequency and wave
      !! direction pair.
      !!
      !! So, to populate the Data4D%DataSet matrix, we simply have to read one line at a
      !! time from the file (stored in RawData4D matrix), and place it in the corresponding
      !! coordinate location in Data4D%DataSet.  This involves simply searching through
      !! the WvFreq1, WvFreq2, WvDir1, and WvDir2 arrays for the correct indices.
      !!
      !! The wave force component direction is stored slightly differently, so it only
      !! needs to be read from the raw data and coverted into an integer to use as an
      !! index.
      !----------------------------------------------------------------------------------


      TmpCoord = 1         ! Initialize the search locations.

      DO I=1,NumDataLinesKeep

            ! Error checking: The LocateStp routine will return 0 if the requested value is less than the value
            !                 of the first element in the array.  It will return the index of the last element
            !                 of the array if the value is larger than the last element.  In creating the arrays
            !                 that are being searched here, the values that we are now trying to find the index
            !                 to were used.  Therefore, if the requested value is not found in the array, then
            !                 something must have gone horribly wrong while creating it, or between then and now,
            !                 which is most likely a programming error.

            ! Find the location in the WvFreq1 array that this point corresponds to.  We will check only between the
            ! cutoffs that were added to the frequency range.  This is contained within TmpWvFreq1 from reading in.
         CALL LocateStp( RawData4D(I,1), TmpWvFreq1, TmpCoord(1),   WvFreq1HiIdx - (WvFreq1LoIdx - 1) )  ! inclusive limits
         IF ( TmpCoord(1) == 0 .OR. ( RawData4D(I,1) > Data4D%WvFreq1(Data4D%NumWvFreq1)) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  Array data point not found in Data4D%WvFreq1 array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         ENDIF
         TmpCoord(1) =  TmpCoord(1) + ( WvFreq1LoIdx - 1 )     ! shift to the point in the Data4D%WvFreq1 array by adding the zero frequency step function

            ! Find the location in the WvFreq2 array that this point corresponds to.  We will check only between the
            ! cutoffs that were added to the frequency range.  This is contained within TmpWvFreq2 from reading in.
         CALL LocateStp( RawData4D(I,2), TmpWvFreq2, TmpCoord(2),   WvFreq2HiIdx - (WvFreq2LoIdx - 1) )  ! inclusive limits
         IF ( TmpCoord(2) == 0 .OR. ( RawData4D(I,2) > Data4D%WvFreq2(Data4D%NumWvFreq2)) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  Array data point not found in Data4D%WvFreq2 array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         ENDIF
         TmpCoord(2) =  TmpCoord(2) + ( WvFreq2LoIdx - 1 )     ! shift to the point in the Data4D%WvFreq2 array by adding the zero frequency step function

            ! Find the location in the WvDir1 array that this point corresponds to.
         CALL LocateStp( RawData4D(I,3), Data4D%WvDir1,  TmpCoord(3),   Data4D%NumWvDir1 )
         IF ( TmpCoord(3) == 0 .OR. ( RawData4D(I,3) > Data4D%WvDir1(Data4D%NumWvDir1)) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  Array data point not found in Data4D%WvDir1 array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         ENDIF

            ! Find the location in the WvDir2 array that this point corresponds to.
         CALL LocateStp( RawData4D(I,4), Data4D%WvDir2,  TmpCoord(4),   Data4D%NumWvDir2 )
         IF ( TmpCoord(4) == 0 .OR. ( RawData4D(I,4) > Data4D%WvDir2(Data4D%NumWvDir2)) ) THEN
            CALL SetErrStat( ErrID_Fatal, ' Programming error.  Array data point not found in Data4D%WvDir2 array.', ErrStat, ErrMsg, RoutineName)
            CALL CleanUp()
            RETURN
         ENDIF

            ! Find which force component this belongs to
         TmpCoord(5) = NINT(RawData4D(I,5))
            ! Check that it is a valid force component
         if (TmpCoord(5) < 1 .or. TmpCoord(5) > 6*Data4D%NumBodies) then
            CALL SetErrStat( ErrID_Fatal, ' Line '//TRIM(Num2Lstr(NumHeaderLines+I))//' of '//TRIM(Filename4D)// &
                           ' contains force component '//TRIM(Num2LStr(TmpCoord(5)))//' which is outside the expected force '// &
                           ' range of 1 to '//TRIM(Num2Lstr(6*Data4D%NumBodies))//' for a '//TRIM(Num2LStr(Data4D%NumBodies))// &
                           ' body system.', ErrStat, ErrMsg, RoutineName)
            IF (ALLOCATED(RawData4D))        DEALLOCATE(RawData4D,STAT=ErrStatTmp)
            IF (ALLOCATED(RawData4DTmp))     DEALLOCATE(RawData4DTmp,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpRealArr))       DEALLOCATE(TmpRealArr,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpDataRow))       DEALLOCATE(TmpDataRow,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpWvFreq1))       DEALLOCATE(TmpWvFreq1,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpWvFreq2))       DEALLOCATE(TmpWvFreq2,STAT=ErrStatTmp)
            CALL CleanUp
            RETURN
         endif

            ! Check that it is a valid force component
         if (TmpCoord(5) < 1 .or. TmpCoord(5) > 6*Data4D%NumBodies) then
            CALL SetErrStat( ErrID_Fatal, ' Line '//TRIM(Num2Lstr(NumHeaderLines+I))//' of '//TRIM(Filename4D)// &
                           ' contains force component '//TRIM(Num2LStr(TmpCoord(5)))//' which is outside the expected force '// &
                           ' range of 1 to '//TRIM(Num2Lstr(6*Data4D%NumBodies))//' for a '//TRIM(Num2LStr(Data4D%NumBodies))// &
                           ' body system.', ErrStat, ErrMsg, RoutineName)
            IF (ALLOCATED(RawData4D))        DEALLOCATE(RawData4D,STAT=ErrStatTmp)
            IF (ALLOCATED(RawData4DTmp))     DEALLOCATE(RawData4DTmp,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpRealArr))       DEALLOCATE(TmpRealArr,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpDataRow))       DEALLOCATE(TmpDataRow,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpWvFreq1))       DEALLOCATE(TmpWvFreq1,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpWvFreq2))       DEALLOCATE(TmpWvFreq2,STAT=ErrStatTmp)
            CALL CleanUp
            RETURN
         endif

            !> The data from the WAMIT file is non-dimensional, so we need to dimensionalize it here.  This
            !! is a partial dimensionalization since the wave amplitudes are not included (this is done later
            !! in each of the calculation methods).  To dimensionalize the data, the equation is for the
            !! partially dimensionalized force (\f$ F_k \f$) is:
            !!
            !!       \f$ F_k = \rho g \cdot L^\alpha \cdot  \bar{F}_k \f$
            !!
            !! where \f$ \bar{F}_k \f$ is the force QTF value in the file for the \f$ k \f$ component
            !! direction, \f$ L \f$ is the WAMIT unit length _WAMITULEN_, \f$ \rho \f$ is the density of
            !! water, and \f$ g \f$ is the gravitational constant (\f$ \rho g \f$ is combined as _RhoXg_).
            !! The value of \f$ \alpha \f$ is 1 for \f$ k = 1,2,3 \f$ and 2 for \f$ k = 4,5,6 \f$.

            ! Here K is used for \f$ alpha \f$ in the dimensionalization equation.
         IF ( TmpCoord(5) <=3 ) THEN
            K = 1
         ELSE
            K = 2
         ENDIF


            ! Check that the current value has not been read in already.  If it has, the corresponding
            ! flag in the mask array will be set to true.  This will indicate that there was redundant
            ! data in the file.  If the value agrees, we will politely ignore it without warning.  If
            ! it does not agree, then we will have to stop since it is ambiguous which is correct.

         IF ( Data4D%DataMask( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4), TmpCoord(5) ) ) THEN
            IF ( .NOT. EqualRealNos( REAL(Data4D%DataSet( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4), TmpCoord(5) ),SiKi), &
                                     REAL(InitInp%RhoXg * InitInp%WAMITULEN**K * RawData4D(I,8)                            ,SiKi)) .AND. &
                 .NOT. EqualRealNos(AIMAG(Data4D%DataSet( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4), TmpCoord(5) )), &
                                     REAL(InitInp%RhoXg * InitInp%WAMITULEN**K * RawData4D(I,9)                            ,SiKi))) THEN
               CALL SetErrStat( ErrID_Fatal, ' Line '//TRIM(Num2Lstr(NumHeaderLines+I))//' of '//TRIM(Filename4D)// &
                        ' contains different values for the real and imaginary part (columns 8 and 9) than was '// &
                        'given earlier in the file for the same values of wave frequency and wave direction '// &
                        '(force dimension = '//TRIM(Num2LStr(TmpCoord(5)))//').', &
                        ErrStat, ErrMsg, RoutineName )
               CALL CleanUp()
               RETURN
            ENDIF
         ELSE

               ! Store the data after dimensionalizing
            Data4D%DataSet( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4), TmpCoord(5) ) = &
                      REAL(InitInp%RhoXg * InitInp%WAMITULEN**K,SiKi) * CMPLX(RawData4D(I,8),RawData4D(I,9),SiKi)

               ! Set flag indicating that this value has been inserted.
            Data4D%DataMask( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4), TmpCoord(5) ) = .TRUE.

         ENDIF

            !> There are relationships between F(Omega_m,Omega_n) and F(Omega_n,Omega_m) where F
            !! resultant excitation force calculated above and stored in Data4D%DataSet.  So, given
            !! F(Omega_m,Omega_n), we can calculate F(Omega_n,Omega_m).  WAMIT only calculates one
            !! of these, so we must generate the other one (we could calculate them later when the
            !! value is used, but that adds complexity that can be avoided by calculating it here).
            !!
            !! The flag Data4D%IsSumForce is set prior to the call to this routine.  This dictates
            !! which method we use.
            !!
            !! For sum forces:
            !!             \f$ \bar{F}^{+}(\omega_m,\omega_n)  = \bar{F}^{+}(\omega_n,\omega_m)       \f$
            !!
            !! For difference forces:
            !!             \f$ \bar{F}^{-}(\omega_m,\omega_n)  = \bar{F}^{- *}(\omega_n,\omega_m)     \f$
            !!
            !!    where \f$ * \f$ indicates the complex conjugate
            !!
            !! @note    For the special case where \f$ \omega_1 = \omega_2 \f$, there is a relationship
            !!          between the wave directions.  Information on this is given in the WAMIT manual,
            !!          page 4-7:
            !!
            !!             \f$ \bar{F}^{-}(\beta_m,\beta_n)  = \bar{F}^{- *}(\beta_n,\beta_m)     \f$
            !!
            !! Be sure to note that this does not necessarily imply that:
            !!             \f$ \bar{F}^{-}(\omega_m,\omega_n,\beta_m,\beta_n)  = \bar{F}^{-}(\omega_n,\omega_m,\beta_n,\beta_m)     \f$
            !!          when \f$ \omega_m \neq \omega_n \f$,


            ! Check if the corresponding value has been filled in already or not.  If not, check if
            ! the value of the frequencies are compatible with each other (will not be true if the
            ! stepsize in the frequencies is different).
            IF ( .NOT. Data4D%DataMask( TmpCoord(2), TmpCoord(1), TmpCoord(3), TmpCoord(4), TmpCoord(5) ) ) THEN

                  ! Equal stepsize check
               IF ( EqualRealNos( Data4D%WvFreq1(TmpCoord(1)), Data4D%WvFreq2(TmpCoord(1)) ) .AND. &
                    EqualRealNos( Data4D%WvFreq1(TmpCoord(1)), Data4D%WvFreq2(TmpCoord(1)) )) THEN

                     ! Value not filled in and the frequencies correspond, so we will now fill it in based
                     ! on what type of data this is (sum or difference).
                  IF ( Data4D%IsSumForce ) THEN
                     Data4D%DataSet( TmpCoord(2), TmpCoord(1), TmpCoord(3), TmpCoord(4), TmpCoord(5) ) = &
                                    Data4D%DataSet( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4), TmpCoord(5) )
                  ELSE  ! Must be difference force, fill with complex conjugate
                     Data4D%DataSet( TmpCoord(2), TmpCoord(1), TmpCoord(3), TmpCoord(4), TmpCoord(5) ) = &
                            CONJG( Data4D%DataSet( TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(4), TmpCoord(5) ) )
                  ENDIF

                     ! We filled in the value, so mark it as filled.
                  Data4D%DataMask( TmpCoord(2), TmpCoord(1), TmpCoord(3), TmpCoord(4), TmpCoord(5) ) = .TRUE.
               ENDIF
            ENDIF

            ! If we have the special case of being on the diagonal of the frequencies where omega_m=omega_n
            ! then there is another rule we can apply to fill in missing values of the wave direction matrix
            ! for this (omega_m,omega_m) pair.
            IF ( EqualRealNos(Data4D%WvFreq1(TmpCoord(1)), Data4D%WvFreq2(TmpCoord(2))) ) THEN
                  ! Only wave direction elements not on the wave direction diagonal
               IF ( .NOT. EqualRealNos( Data4D%WvDir1(TmpCoord(3)), Data4D%WvDir2(TmpCoord(4)) ) ) THEN

                     ! See if WvDir1(J) == WvDir2(J) and WvDir1(K) == WvDir2(K), because
                     ! if they are not, then the discretization along the two wave directions
                     ! is different and this won't work.
                  IF ( EqualRealNos( Data4D%WvDir1(TmpCoord(3)), Data4D%WvDir2(TmpCoord(3)) ) .AND. &
                       EqualRealNos( Data4D%WvDir1(TmpCoord(4)), Data4D%WvDir2(TmpCoord(4)) ) ) THEN

                        ! See if the diagonal mirror one (WvDir2,WvDir1) value is not filled, set it and its flag
                     IF ( .NOT. Data4D%DataMask(TmpCoord(1), TmpCoord(2), TmpCoord(4), TmpCoord(3), TmpCoord(5))  ) THEN
                        Data4D%DataSet(TmpCoord(1), TmpCoord(2), TmpCoord(4), TmpCoord(3), TmpCoord(5) ) = &
                              Data4D%DataSet(TmpCoord(1), TmpCoord(2), TmpCoord(3), TmpCoord(3), TmpCoord(5) )
                        Data4D%DataMask(TmpCoord(1), TmpCoord(2), TmpCoord(4), TmpCoord(3), TmpCoord(5) ) = .TRUE.
                     ENDIF

                  ENDIF ! Check that wave directions will pair.

               ENDIF

            ENDIF    ! WvFreq1(TmpCoord(1)) == WvFreq2(TmpCoord(2))



      ENDDO


      !----------------------------------------------------------------------------------
      !> We added two frequencies for the \f$ omega = 0 \f$ term if it did not exist,
      !! and added two frequencies for the infinite frequency term on the end of the array,
      !! to create step changes outside the entered frequency ranges. We need to populate
      !! the these new terms (set to zero).
      !----------------------------------------------------------------------------------

      IF (.NOT. HaveZeroFreq1) THEN
         Data4D%DataSet( 1:2,:,:,:,:)  = CMPLX(0.0,0.0,SiKi)                                          ! Set the values to zero for everything before entered frequency range
         Data4D%DataMask(1:2,:,:,:,:)  = .TRUE.                                                       ! Set the mask for these first two frequencies
      ENDIF
      Data4D%DataSet( Data4D%NumWvFreq1-1:Data4D%NumWvFreq1,:,:,:,:) = CMPLX(0.0,0.0,SiKi)            ! Set the values for the last two frequencies to zero (everything higher than the last non-infinite frequency)
      Data4D%DataMask(Data4D%NumWvFreq1-1:Data4D%NumWvFreq1,:,:,:,:) = .TRUE.                         ! Set the mask for the last two frequencies

      IF (.NOT. HaveZeroFreq2) THEN
         Data4D%DataSet( :,1:2,:,:,:)  = CMPLX(0.0,0.0,SiKi)                                          ! Set the values to zero for everything before entered frequency range
         Data4D%DataMask(:,1:2,:,:,:)  = .TRUE.                                                       ! Set the mask for these first two frequencies
      ENDIF
      Data4D%DataSet( :,Data4D%NumWvFreq2-1:Data4D%NumWvFreq2,:,:,:) = CMPLX(0.0,0.0,SiKi)            ! Set the values for the last two frequencies to zero (everything higher than the last non-infinite frequency)
      Data4D%DataMask(:,Data4D%NumWvFreq2-1:Data4D%NumWvFreq2,:,:,:) = .TRUE.                         ! Set the mask for the last two frequencies


      !----------------------------------------------------------------------------------
      !> Find out if the data is sparse or full.  Verification that the requested component
      !! directions were found will occur in the calling routine, not here (that information
      !! was not passed in).  We will check this by sweeping through all three dimensions
      !! for only the values of k that have data in them.
      !----------------------------------------------------------------------------------

      DO M=1,6*Data4D%NumBodies       ! Loop over force component directions
         TmpSparseFlag  = .FALSE.                  ! Change this to true if any empty values are found
         IF (Data4D%LoadComponents(M)) THEN        ! Only if we found data for that component
            DO I=1,Data4D%NumWvFreq1
               DO J=1,Data4D%NumWvFreq2
                  DO K=1,Data4D%NumWvDir1
                     DO L=1,Data4D%NumWvDir2
                        IF (.NOT. Data4D%DataMask( I, J, K, L, M ) ) THEN
                           TmpSparseFlag = .TRUE.
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO
         ENDIF

            ! If any values were missing for this force component direction, TmpSparseFlag will be true.
         Data4D%DataIsSparse(M)  = TmpSparseFlag

      ENDDO



      !----------------------------------------------------------------------------------
      !> If the data is sparse, it may still have a complete diagonal.  This is worth
      !! checking in case this data is used for either the MnDrift or NewmanApp methods.
      !!
      !! Diagonal terms: Elements where the first wave frequency and second wave frequency match.
      !!                 For the diagonal to be complete, we require all the wave direction
      !!                 pairs for \f$ \omega_m = \omega_n \f$ to be complete.
      !!
      !! NOTE: In the 3D WAMIT files, only the diagonal terms Omega_m = Omega_m are given.  In
      !! the 4D files, the off diagonal terms are also listed.  For calculations that involve the
      !! full QTF (such as DiffQTF and SumQTF), the data must be complete to use a normal
      !! interpolation routine (may be sparse if a sparse interpolation routine is developed).
      !! For the MnDrift and NewmanApp calculations, incomplete 4D QTFs may be used if the
      !! diagonal is complete.  So, we will do some additional testing on the 4D data to make
      !! sure that the diagonal is complete.  This will allow the data to be used easily with
      !! the MnDrift and NewmanApp.
      !!
      !! We will only be considering this for the load components requested.  Only if all the
      !! requested component diagonals are complete will we consider it complete.
      !!
      !! One consideration worth noting is that the diagonal may not be present in a complete
      !! dataset.  This would occur if the stepsize in \f$ \omega_1 \f$ and \f$ \omega_2 \f$ are
      !! different.
      !----------------------------------------------------------------------------------

      TmpDiagComplete = .TRUE.                     ! Change this to false if any empty values are found

         ! If the number of frequencies differ, than this won't work.
      IF (Data4D%NumWvFreq1 /= Data4D%NumWvFreq2) THEN
         TmpDiagComplete = .FALSE.
      ELSE  ! Same number of frequencies, so we can proceed.
         DO M=1,6*Data4D%NumBodies       ! Loop over force component directions
               ! If we have data for this component, and it is sparse, proceed to check it.
            IF (Data4D%LoadComponents(M) .AND. Data4D%DataIsSparse(M)) THEN

                  ! Step through all the components
               DO I=1,Data4D%NumWvFreq1
                     ! Make sure the frequencies actually match
                  IF ( EqualRealNos(Data4D%WvFreq1(I),Data4D%WvFreq2(I)) ) THEN
                     DO K=1,Data4D%NumWvDir1
                        DO L=1,Data4D%NumWvDir2
                           IF (.NOT. Data4D%DataMask( I, I, K, L, M ) ) THEN
                              TmpDiagComplete = .FALSE.         ! If any of the elements are missing
                           ENDIF
                        ENDDO
                     ENDDO
                  ELSE  ! Frequencies don't match, so this won't work.
                     TmpDiagComplete = .FALSE.
                  ENDIF
               ENDDO ! NumWvFreq1

            ENDIF    ! Sparse load component of interest
         ENDDO ! Loop over load components
      ENDIF

      Data4D%WvFreqDiagComplete = TmpDiagComplete

      call cleanup()
      
      contains
         subroutine cleanup()

               ! Clean up
            IF (ALLOCATED(RawData4D))        DEALLOCATE(RawData4D,STAT=ErrStatTmp)
            IF (ALLOCATED(RawData4DTmp))     DEALLOCATE(RawData4DTmp,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpRealArr))       DEALLOCATE(TmpRealArr,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpDataRow))       DEALLOCATE(TmpDataRow,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpWvFreq1))       DEALLOCATE(TmpWvFreq1,STAT=ErrStatTmp)
            IF (ALLOCATED(TmpWvFreq2))       DEALLOCATE(TmpWvFreq2,STAT=ErrStatTmp)
            
         end subroutine cleanup
         
   END SUBROUTINE Read_DataFile4D






   !> This subroutine counts the number of unique values in an array and returns a sorted array of them.
   !! This is called by the _Read_DataFile3D_ and _Read_DataFile4D_ routines.
   SUBROUTINE UniqueRealValues( DataArrayIn, DataArrayOut, NumUnique, ErrStat, ErrMsg )
      IMPLICIT NONE

      REAL(SiKi),                         INTENT(IN   )  :: DataArrayIn(:)    !< Array to search
      REAL(SiKi),       ALLOCATABLE,      INTENT(  OUT)  :: DataArrayOut(:)   !< Array to return results in
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumUnique         !< Number of unique values found
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           !< Error Status
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            !< Message about the error

         ! Local variables
!      REAL(SiKi)                                         :: TmpReal           !< Temporary real value
      INTEGER(IntKi)                                     :: I                 !< Generic counter
      INTEGER(IntKi)                                     :: J                 !< Generic counter
      REAL(SiKi),       ALLOCATABLE                      :: TmpRealArray(:)   !< Temporary real array
      LOGICAL                                            :: Duplicate         !< If there is a duplicate value

         ! Error handling
      INTEGER(IntKi)                                     :: ErrStatTmp
      CHARACTER(2048)                                    :: ErrMsgTmp
      CHARACTER(*), PARAMETER                            :: RoutineName = 'UniqueRealValues'



         ! Initialize things
      ErrStat     =  ErrID_None
      ErrMsg      =  ""


         ! Allocate the temporary array
      CALL AllocAry( TmpRealArray, SIZE(DataArrayIn,1), 'Temporary array for data sorting', ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Initialize the array with a large negative number.  Only positive frequencies and wave directions from
         ! +/-180 degrees should ever be sent here.  The way this algorithm works, only empty values will have
         ! this number.  We won't be returning it to the calling routine
      TmpRealArray   = -9.9e9_SiKi


         ! The first point is unique since we haven't compared it to anything yet.
      TmpRealArray(1)   = DataArrayIn(1)
      NumUnique         =  1

         ! Step through the DataArrayIn and put unique values into TmpRealArray.  Start at second point
      DO I=2,SIZE(DataArrayIn,1)
            ! Check the current value against the largest stored value (I-1).  If the current value is
            ! larger than the last stored one, then it should be stored after it.
         IF ( DataArrayIn(I) > TmpRealArray(NumUnique) ) THEN
            TmpRealArray(NumUnique + 1) = DataArrayIn(I)
            NumUnique = NumUnique + 1
         ELSE
            ! Otherwise, if the value should not be put last, then we have to find where it goes.  Before
            ! searching for the location, first make sure this isn't a duplicate value.
            Duplicate = .FALSE.
            DO J= NumUnique, 1, -1
               IF ( EqualRealNos( DataArrayIn(I), TmpRealArray(J) )) THEN
                  Duplicate = .TRUE.
                  EXIT     ! Stop searching
               ENDIF
            ENDDO


            ! If this is not a duplicate, the location where it goes has to be find.  To do this, we will
            ! sequentially shift each value one index further as we step backwards through the sorted
            ! array.  When we find the location between values where this goes, we put the value there.
            IF ( .NOT. Duplicate ) THEN
               DO J= NumUnique, 1, -1           ! TempRealArray only has NumUnique values.  Everything after is junk.
                  IF ( DataArrayIn(I) < TmpRealArray(J) )   THEN
                     TmpRealArray(J+1) = TmpRealArray(J)       ! Shift this value further along the array
                     IF ( J == 1 )  THEN                       ! If we are at the first value, then it goes here.
                        TmpRealArray(J) = DataArrayIn(I)
                        NumUnique = NumUnique + 1
                     ELSE
                        IF ( DataArrayIn(I) > TmpRealArray(J-1) ) THEN  ! If larger than the preceeding number, it goes here
                           TmpRealArray(J) = DataArrayIn(I)
                           NumUnique = NumUnique + 1
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO


         ! Now that we have the array sorted into unique values, we need to construct an array to return the values in.
      CALL AllocAry( DataArrayOut, NumUnique, 'Return array with sorted values', ErrStatTmp, ErrMsgTmp )
      CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL CleanUp()
         RETURN
      ENDIF

         ! Copy the values over
      DataArrayOut = TmpRealArray(1:NumUnique)

      call cleanup()
      
      contains
         subroutine cleanup()
            if (allocated(TmpRealArray)) deallocate(TmpRealArray)
         end subroutine cleanup
      
   END SUBROUTINE UniqueRealValues






   !-------------------------------------------------------------------------------------------------------------------------------
   !>    This subroutine looks at a file that has been opened and finds out how many header lines there are, how many periods
   !!    (frequencies) there are (first only if there are paired periods for second order), and how many lines of data there are in
   !!    the file.
   !!
   !!    A few things are assumed about the file:
   !!       1. Any header lines are the first thing in the file.
   !!       2. No text appears anyplace other than in the file
   !!       3. The datalines only contain numbers that can be read in as reals.
   !!
   !!    Limitations:
   !!       1. only handles up to 20 words (columns) on a line
   !!       2. empty lines are considered text lines
   !!       3. All data rows must contain the same number of columns
   !!
   !!
   SUBROUTINE GetFileLength(UnitDataFile, Filename, NumDataColumns, NumDataLines, NumHeaderLines, ErrStat, ErrMsg)

      IMPLICIT NONE

         ! Passed variables
      INTEGER(IntKi),                     INTENT(IN   )  :: UnitDataFile          !< Unit number of the file we are looking at.
      CHARACTER(*),                       INTENT(IN   )  :: Filename          !< The name of the file we are looking at.
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumDataColumns    !< The number of columns in the data file.
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumDataLines      !< Number of lines containing data
      INTEGER(IntKi),                     INTENT(  OUT)  :: NumHeaderLines    !< Number of header lines at the start of the file
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            !< Error Message to return (empty if all good)
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           !< Status flag if there were any problems (ErrID_None if all good)

         ! Local Variables
      CHARACTER(2048)                                    :: ErrMsgTmp         !< Temporary message variable.  Used in calls.
      INTEGER(IntKi)                                     :: ErrStatTmp        !< Temporary error status.  Used in calls.
      INTEGER(IntKi)                                     :: LclErrStat        !< Temporary error status.  Used locally to indicate when we have reached the end of the file.
      INTEGER(IntKi)                                     :: TmpIOErrStat      !< Temporary error status for the internal read of the first word to a real number
      LOGICAL                                            :: IsRealNum         !< Flag indicating if the first word on the line was a real number

      CHARACTER(MaxFileInfoLineLen)                      :: TextLine          !< One line of text read from the file
      INTEGER(IntKi)                                     :: LineLen           !< The length of the line read in
      CHARACTER(1024)                                    :: StrRead           !< String containing the first word read in
      REAL(SiKi)                                         :: RealRead          !< Returns value of the number (if there was one), or NaN (as set by NWTC_Num) if there wasn't
!      CHARACTER(1024)                                    :: VarName           !< Name of the variable we are trying to read from the file
      CHARACTER(NWTC_SizeOfNumWord)                      :: Words(20)         !< Array of words we extract from a line.  We shouldn't have more than 20.
      INTEGER(IntKi)                                     :: i !,j,k             !< simple integer counters
      INTEGER(IntKi)                                     :: LineNumber        !< the line I am on
      LOGICAL                                            :: LineHasText       !< Flag indicating if the line I just read has text.  If so, it is a header line.
      LOGICAL                                            :: HaveReadData      !< Flag indicating if I have started reading data.
      INTEGER(IntKi)                                     :: NumWords          !< Number of words on a line
      INTEGER(IntKi)                                     :: FirstDataLineNum  !< Line number of the first row of data in the file
      CHARACTER(*), PARAMETER                            :: RoutineName = 'GetFileLength'



         ! Initialize the error handling
      ErrStat     = ErrID_None
      ErrStatTmp  = ErrID_None
      LclErrStat  = ErrID_None
      ErrMsg      = ''
      ErrMsgTmp   = ''


         ! Set some of the flags and counters
      HaveReadData   = .FALSE.
      NumDataColumns = 0
      NumHeaderLines = 0
      NumDataLines   = 0
      LineNumber     = 0


         ! Just in case we were handed a file that we are part way through reading (should never be true), rewind to the start

      REWIND( UnitDataFile )


         !------------------------------------
         !> The variable LclErrStat is used to indicate when we have reached the end of the file or had an error from
         !! ReadLine.  Until that occurs, we read each line, and decide if it contained any non-numeric data.  The
         !! first group of lines containing non-numeric data is considered the header.  The first line of all numeric
         !! data is considered the start of the data section.  Any non-numeric containing found within the data section
         !! will be considered as an invalid file format at which point we will return a fatal error from this routine.

      DO WHILE ( LclErrStat == ErrID_None )

            !> Reset the indicator flag for the non-numeric content
         LineHasText = .FALSE.

            !> Read in a single line from the file
         CALL ReadLine( UnitDataFile, '', TextLine, LineLen, LclErrStat )

            !> If there was an error in reading the file, then exit.
            !!    Possible causes: reading beyond end of file in which case we are done so don't process it.
         IF ( LclErrStat /= ErrID_None ) EXIT

            !> Increment the line counter.
         LineNumber  = LineNumber + 1

            !> Read all the words on the line into the array called 'Words'.  Only the first words will be encountered
            !! will be stored.  The others are empty (i.e. only three words on the line, so the remaining 17 are empty).
         CALL GetWords( TextLine, Words, size(Words), NumWords )


            !> Now cycle through the first 'NumWords' of non-empty values stored in 'Words'.  Words should contain
            !! everything that is one the line.  The subroutine ReadRealNumberFromString will set a flag 'IsRealNum'
            !! when the value in Words(i) can be read as a real(SiKi).  'StrRead' will contain the string equivalent.
         DO i=1,NumWords
            CALL ReadRealNumberFromString( Words(i), RealRead, StrRead, IsRealNum, ErrStatTmp, ErrMsgTmp, TmpIOErrStat )
            IF ( .NOT. IsRealNum) THEN
               LineHasText = .TRUE.
            ENDIF
         ENDDO

            !> If all the words on that line had no text in them, then it must have been a line of data.
            !! If not, then we have either a header line, which is ok, or a line containing text in the middle of the
            !! the data section, which is not good (the flag HaveReadData tells us which case this is).
         IF ( LineHasText ) THEN
            IF ( HaveReadData ) THEN      ! Uh oh, we have already read a line of data before now, so there is a problem
               CALL SetErrStat( ErrID_Fatal, ' Found text on line '//TRIM(Num2LStr(LineNumber))//' of '//TRIM(FileName)// &
                           ' when real numbers were expected.  There may be a problem with the 2nd order WAMIT file: '// &
                           TRIM(Filename)//'.', ErrStat, ErrMsg, RoutineName)
               RETURN
            ELSE
               NumHeaderLines = NumHeaderLines + 1
            ENDIF
         ELSE     ! No text, must be data line
            NumDataLines = NumDataLines + 1
               ! If this is the first row of data, then store the number of words that were on the line
            IF ( .NOT. HaveReadData )  THEN
                  ! If this is the first line of data, keep some relevant info about it and the number of columns in it
               HaveReadData      = .TRUE.
               FirstDataLineNum  = LineNumber         ! Keep the line number of the first row of data (for error reporting)
               NumDataColumns    = NumWords
            ELSE
                  ! Make sure that the number columns on the row matches the number of columnns on the first row of data.
               IF ( NumWords /= NumDataColumns ) THEN
                  CALL SetErrStat( ErrID_Fatal, ' Error in 2nd order WAMIT file: '//TRIM(Filename)//'.'// &
                           ' The number of data columns on line '//TRIM(Num2LStr(LineNumber))// &
                           '('//TRIM(Num2LStr(NumWords))//' columns) is different than the number of columns on first row of data '// &
                           ' (line: '//TRIM(Num2LStr(FirstDataLineNum))//', '//TRIM(Num2LStr(NumDataColumns))//' columns).', &
                           ErrStat, ErrMsg, RoutineName)
                  RETURN
               ENDIF
            ENDIF
         ENDIF

      ENDDO

      IF ( NumDataLines < 2 ) THEN
         CALL SetErrStat( ErrID_Fatal, ' 2nd order WAMIT file '//TRIM(Filename)//' contains only '//TRIM(Num2LStr(NumDataLines))// &
                        ' lines of data. This does not appear to be a valid WAMIT file.', ErrStat, ErrMsg, RoutineName)
         RETURN
      ENDIF

      REWIND( UnitDataFile )

   END SUBROUTINE GetFileLength

   !-------------------------------------------------------------------------------
   !> This subroutine takes a line of text that is passed in and reads the first
   !! word to see if it is a number.  An internal read is used to do this.  If
   !! it is a number, it is started in ValueRead and returned. The flag IsRealNum
   !! is set to true.  Otherwise, ValueRead is set to NaN (value from the NWTC_Num)
   !! and the flag is set to false.
   !!
   !! The IsRealNum flag is set to indicate if we actually have a real number or
   !! not.  After calling this routine, a simple if statement can be used:
   !!
   !!       @code
   !!    IF (IsRealNum) THEN
   !!       ! do something
   !!    ELSE
   !!       ! do something else
   !!    ENDIF
   !!       @endcode
   !!
   !-------------------------------------------------------------------------------
   SUBROUTINE ReadRealNumberFromString(StringToParse, ValueRead, StrRead, IsRealNum, ErrStat, ErrMsg, IOErrStat)

      CHARACTER(*),        INTENT(IN   )           :: StringToParse  !< The string we were handed.
      REAL(SiKi),          INTENT(  OUT)           :: ValueRead      !< The variable being read.  Returns as NaN (library defined) if not a Real.
      CHARACTER(*),        INTENT(  OUT)           :: StrRead        !< A string containing what was read from the ReadNum routine.
      LOGICAL,             INTENT(  OUT)           :: IsRealNum      !< Flag indicating if we successfully read a Real
      INTEGER(IntKi),      INTENT(  OUT)           :: ErrStat        !< ErrID level returned from ReadNum
      CHARACTER(*),        INTENT(  OUT)           :: ErrMsg         !< Error message including message from ReadNum
      INTEGER(IntKi),      INTENT(  OUT)           :: IOErrStat      !< Error status from the internal read. Useful for diagnostics.

      CHARACTER(2048)                              :: ErrMsgTmp            !< Temporary variable for holding the error message returned from a CALL statement


         ! Initialize some things
      ErrStat     = ErrID_None
      ErrMsg      = ''


         ! ReadNum returns a string contained in StrRead.  So, we now try to do an internal read to VarRead and then trap errors.
      read(StringToParse,*,IOSTAT=IOErrStat)   StrRead
      read(StringToParse,*,IOSTAT=IOErrStat)   ValueRead


         ! If IOErrStat==0, then we have a real number, anything else is a problem.
      if (IOErrStat==0) then
         IsRealNum   = .TRUE.
      else
         IsRealNum   = .FALSE.
         ValueRead   = NaN_S           ! This is NaN as defined in the NWTC_Num.
         ErrMsg      = 'Not a real number. '//TRIM(ErrMsgTmp)//NewLine
         ErrSTat     = ErrID_Severe
      endif



      RETURN
   END SUBROUTINE ReadRealNumberFromString


   !-------------------------------------------------------------------------------------------------------------------------------
   !-------------------------------------------------------------------------------
   !> This subroutine works with the ReadNum routine from the library.  ReadNum is
   !! called to read a word from the input file.  An internal read is then done to
   !! convert the string to a number that is stored in VarRead and returned.
   !!
   !! The IsRealNum flag is set to indicate if we actually have a real number or
   !! not.  After calling this routine, a simple if statement can be used:
   !!
   !!       @code
   !!    IF (ISRealNum) THEN
   !!       ! do something
   !!    ELSE
   !!       ! do something else
   !!    ENDIF
   !!       @endcode
   !!
   !-------------------------------------------------------------------------------
   SUBROUTINE ReadRealNumber(UnitNum, FileName, VarName, VarRead, StrRead, IsRealNum, ErrStat, ErrMsg, IOErrStat)

      INTEGER(IntKi),      INTENT(IN   )           :: UnitNum        !< The unit number of the file being read
      CHARACTER(*),        INTENT(IN   )           :: FileName       !< The name of the file being read.  Used in the ErrMsg from ReadNum (Library routine).
      CHARACTER(*),        INTENT(IN   )           :: VarName        !< The variable we are reading.  Used in the ErrMsg from ReadNum (Library routine)'.
      REAL(SiKi),          INTENT(  OUT)           :: VarRead        !< The variable being read.  Returns as NaN (library defined) if not a Real.
      CHARACTER(*),        INTENT(  OUT)           :: StrRead        !< A string containing what was read from the ReadNum routine.
      LOGICAL,             INTENT(  OUT)           :: IsRealNum      !< Flag indicating if we successfully read a Real
      INTEGER(IntKi),      INTENT(  OUT)           :: ErrStat        !< ErrID level returned from ReadNum
      CHARACTER(*),        INTENT(  OUT)           :: ErrMsg         !< Error message including message from ReadNum
      INTEGER(IntKi),      INTENT(  OUT)           :: IOErrStat      !< Error status from the internal read. Useful for diagnostics.

         ! Local vars
      INTEGER(IntKi)                      :: ErrStatTmp
      CHARACTER(2048)                     :: ErrMsgTmp



         ! Initialize some things
      ErrStat     = ErrID_None
      ErrMsg      = ''


         ! Now call the ReadNum routine to get the number
         ! If it is a word that does not start with T or F, then ReadNum won't give any errors.
      CALL ReadNum( UnitNum, FileName, StrRead, VarName, ErrStatTmp, ErrMsgTmp)


         ! ReadNum returns a string contained in StrRead.  So, we now try to do an internal read to VarRead and then trap errors.
      read(StrRead,*,IOSTAT=IOErrStat)   VarRead


         ! If IOErrStat==0, then we have a real number, anything else is a problem.
      if (IOErrStat==0) then
         IsRealNum   = .TRUE.
      else
         IsRealNum   = .FALSE.
         VarRead     = NaN_S        ! This is NaN as defined in the NWTC_Num.
         ErrMsg      = 'Not a real number. '//TRIM(ErrMsgTmp)//NewLine
         ErrStat     = ErrStatTmp         ! The ErrStatTmp returned by the ReadNum routine is an ErrID level.
      endif



      RETURN
   END SUBROUTINE ReadRealNumber






!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE WAMIT2_CalcOutput( Time, WaveTime, p, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: Time           !< Current simulation time in seconds
      real(SiKi),                         intent(in   )  :: WaveTime(:)    !< Array of wave kinematic time samples, (sec)
      TYPE(WAMIT2_ParameterType),         INTENT(IN   )  :: p              !< Parameters
      TYPE(WAMIT2_OutputType),            INTENT(INOUT)  :: y              !< Outputs computed at Time (Input only so that mesh
                                                                           !!   connectivity information does not have to be recalculated)
      TYPE(WAMIT2_MiscVarType),           INTENT(INOUT)  :: m              !< Misc/optimization variables
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat        !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None



         ! Local Variables:
      INTEGER(IntKi)                                     :: I              ! Generic index
      INTEGER(IntKi)                                     :: IBody          ! Index to body number
      INTEGER(IntKi)                                     :: indxStart      ! Starting index


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""




         ! Abort if the wave excitation loads have not been computed yet:

      IF ( .NOT. ALLOCATED ( p%WaveExctn2 ) )  THEN
         CALL SetErrStat(ErrID_Fatal,' Routine WAMIT2_Init() must be called before routine WAMIT2_CalcOutput().',ErrStat,ErrMsg,'WAMIT2_CalcOutput')
         RETURN
      END IF


         ! Compute the 2nd order load contribution from incident waves:

      do iBody = 1, p%NBody
         indxStart = (iBody-1)*6

         DO I = 1,6     ! Loop through all wave excitation forces and moments
            m%F_Waves2(indxStart+I) = InterpWrappedStpReal ( REAL(Time, SiKi), WaveTime(:), p%WaveExctn2(:,indxStart+I), &
                                                  m%LastIndWave(IBody), p%NStepWave + 1       )
         END DO          ! I - All wave excitation forces and moments



         ! Copy results to the output point mesh
         DO I=1,3
            y%Mesh%Force(I,IBody)    = m%F_Waves2(indxStart+I)
         END DO
         DO I=1,3
            y%Mesh%Moment(I,IBody)   = m%F_Waves2(indxStart+I+3)
         END DO

      enddo


END SUBROUTINE WAMIT2_CalcOutput




!-------------------------------------------------------------------------------
!> This subroutine copies data stored from a W2_InitData4D_Type to a W2_InitData3D_Type
SUBROUTINE Copy_InitData4Dto3D( Data4D, Data3D, ErrStat, ErrMsg )
   IMPLICIT NONE

   TYPE(W2_InitData4D_Type),  INTENT(IN   )  :: Data4D         !< 4D data source
   TYPE(W2_InitData3D_Type),  INTENT(INOUT)  :: Data3D         !< 3D data where it will be copied to
   INTEGER(IntKi),            INTENT(  OUT)  :: ErrStat        !< Error status
   CHARACTER(*),              INTENT(  OUT)  :: ErrMsg         !< Error message

      ! Local variables
   LOGICAL                                   :: TmpSparseFlag  !< Flag to indicate if the data is sparse
   INTEGER(IntKi)                            :: I              !< Generic counter
   INTEGER(IntKi)                            :: J              !< Generic counter
   INTEGER(IntKi)                            :: K              !< Generic counter
   INTEGER(IntKi)                            :: L              !< Generic counter
   INTEGER(IntKi)                            :: ErrStatTmp     !< Temporary error status for calls
   CHARACTER(2048)                           :: ErrMsgTmp      !< Temporary error message for calls
   CHARACTER(*), PARAMETER                   :: RoutineName = 'Copy_InitData4Dto3D'


      ! Initialize the error handling
   ErrStat     = ErrID_none
   ErrStatTmp  = ErrID_none
   ErrMsg      = ''
   ErrMsgTmp   = ''


      ! Make sure we aren't trying to copy 4D sum frequency data in to a 3D type that only holds information for difference frequencies
   IF ( Data4D%IsSumForce ) THEN
      CALL SetErrStat( ErrID_Fatal, ' Attempted to copy 4D sum-frequency data into a 3D difference frequency type.', ErrStat, ErrMsg,RoutineName)
      RETURN
   ENDIF

      ! Make sure the dimensions work
   IF ( Data4D%NumWvFreq1 /= Data4D%NumWvFreq2 ) THEN
      CALL SetErrStat( ErrID_Fatal, ' Attempted to copy 4D data to 3D data when NumFreq1 /= NumFreq2.', ErrStat, ErrMsg,RoutineName)
      RETURN
   ENDIF

      ! Make sure that the frequencies actually match each other
   DO I=1,Data4D%NumWvFreq1
      IF ( .NOT. EqualRealNos(Data4D%WvFreq1(I), Data4D%WvFreq2(I)) ) THEN
         CALL SetErrStat( ErrID_Fatal, ' Attempted to copy 4D data to 3D data when wave frequencies are not the same.', ErrStat, ErrMsg,RoutineName)
         RETURN
      ENDIF
   ENDDO


      ! Make sure that nothing has been allocated already
   ErrStatTmp = ErrID_None
   IF (ALLOCATED(Data3D%DataSet))   ErrStatTmp = ErrID_Fatal
   IF (ALLOCATED(Data3D%DataMask))  ErrStatTmp = ErrID_Fatal
   IF (ALLOCATED(Data3D%WvFreq1))   ErrStatTmp = ErrID_Fatal
   IF (ALLOCATED(Data3D%WvDir1))    ErrStatTmp = ErrID_Fatal
   IF (ALLOCATED(Data3D%WvDir2))    ErrStatTmp = ErrID_Fatal
   IF ( ErrStatTmp >= ErrID_Fatal ) THEN
      CALL SetErrStat( ErrID_Fatal, ' Attempted to copy 4D data into a populated 3D dataset.', ErrStat, ErrMsg,RoutineName)
   ENDIF



      ! Ok.  It checks out, so start copying info over
   Data3D%NumWvFreq1       =  Data4D%NumWvFreq1
   Data3D%NumWvDir1        =  Data4D%NumWvDir1
   Data3D%NumWvDir2        =  Data4D%NumWvDir2
   Data3D%NumBodies        =  Data4D%NumBodies


      ! Now allocate the storage arrays
   ALLOCATE( Data3D%DataSet( Data3D%NumWvFreq1, Data3D%NumWvDir1, Data3D%NumWvDir2, 6*Data3D%NumBodies ),  STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array Data3D%DataSet to store '// &
                           'the 3D 2nd order WAMIT data.',  ErrStat,ErrMsg,RoutineName)
   ALLOCATE( Data3D%DataMask( Data3D%NumWvFreq1, Data3D%NumWvDir1, Data3D%NumWvDir2, 6*Data3D%NumBodies ),  STAT=ErrStatTmp )
   IF (ErrStatTmp /= 0) CALL SetErrStat(ErrID_Fatal,'Cannot allocate array Data3D%DataMask to store '// &
                           'the information on the 3D 2nd order WAMIT data.',  ErrStat,ErrMsg,RoutineName)
   CALL AllocAry( Data3D%WvFreq1, Data3D%NumWvFreq1, 'Data3D WvFreq array', ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( Data3D%WvDir1, Data3D%NumWvDir1, 'Data3D WvDir1 array', ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( Data3D%WvDir2, Data3D%NumWvDir2, 'Data3D WvDir2 array', ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( Data3D%LoadComponents, 6*Data3D%NumBodies, 'Data3D LoadComponents array', ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( Data3D%DataIsSparse, 6*Data3D%NumBodies, 'Data3D DataIsSparse array', ErrStatTmp, ErrMsgTmp )
   CALL SetErrStat( ErrStatTmp, ErrMsgTmp, ErrStat, ErrMsg, RoutineName )

   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Destroy_InitData3D( Data3D )
      RETURN
   ENDIF


      ! Copy the values over
   Data3D%WvFreq1    = Data4D%WvFreq1
   Data3D%WvDir1     = Data4D%WvDir1
   Data3D%WvDir2     = Data4D%WvDir2
   Data3D%LoadComponents   =  Data4D%LoadComponents

   DO I=1,Data3D%NumWvFreq1
      Data3D%DataSet(I,:,:,:)    = Data4D%DataSet(I,I,:,:,:)
      Data3D%DataMask(I,:,:,:)   = Data4D%DataMask(I,I,:,:,:)
   ENDDO


      !----------------------------------------------------------------------------------
      !> This routine should not have been called unless the data could be copied over,
      !! which means that the 4D data is either not sparse in the dimensions of interest,
      !! or the diagonal is complete.  We won't simply trust that this is true though, so
      !! we will run the array and find out if the data we now have is sparse for the
      !! load components we have or not.
      !!
      !! Find out if the data is sparse or full.  Verification that the requested component
      !! directions were found will occur in the calling routine, not here (that information
      !! was not passed in).  We will check this by sweeping through all three dimensions
      !! for only the values of k that have data in them.
      !----------------------------------------------------------------------------------

   DO L=1,6*Data3D%NumBodies       ! Loop over force component directions
      TmpSparseFlag  = .FALSE.                  ! Change this to true if any empty values are found
      IF (Data3D%LoadComponents(L)) THEN        ! Only if we found data for that component
         DO I=1,Data3D%NumWvFreq1
            DO J=1,Data3D%NumWvDir1
               DO K=1,Data3D%NumWvDir2
                  IF (.NOT. Data3D%DataMask( I, J, K, L ) ) THEN
                     TmpSparseFlag = .TRUE.
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF

         ! If any values were missing for this force component direction, TmpSparseFlag will be true.
      Data3D%DataIsSparse(L)  = TmpSparseFlag

   ENDDO    ! Sparseness check


      ! Cleanup.  Well, nothing here to cleanup.

END SUBROUTINE Copy_InitData4Dto3D



!> This subroutine destroys data stored in a W2_InitData3D_Type
SUBROUTINE Destroy_InitData3D(Data3D)
   IMPLICIT NONE
   TYPE(W2_InitData3D_Type),  INTENT(INOUT)  :: Data3D
   INTEGER(IntKi)                            :: ErrStatTmp
   
   IF (ALLOCATED(Data3D%DataIsSparse))   DEALLOCATE(Data3D%DataIsSparse,STAT=ErrStatTmp)
   IF (ALLOCATED(Data3D%LoadComponents)) DEALLOCATE(Data3D%LoadComponents,STAT=ErrStatTmp)
   
   IF (ALLOCATED(Data3D%DataSet))      DEALLOCATE(Data3D%DataSet,STAT=ErrStatTmp)
   IF (ALLOCATED(Data3D%DataMask))     DEALLOCATE(Data3D%DataMask,STAT=ErrStatTmp)
   IF (ALLOCATED(Data3D%WvFreq1))      DEALLOCATE(Data3D%WvFreq1,STAT=ErrStatTmp)
   IF (ALLOCATED(Data3D%WvDir1))       DEALLOCATE(Data3D%WvDir1,STAT=ErrStatTmp)
   IF (ALLOCATED(Data3D%WvDir2))       DEALLOCATE(Data3D%WvDir2,STAT=ErrStatTmp)
END SUBROUTINE Destroy_InitData3D


!> This subroutine destroys data stored in a W2_InitData4D_Type
SUBROUTINE Destroy_InitData4D(Data4D)
   IMPLICIT NONE
   TYPE(W2_InitData4D_Type),  INTENT(INOUT)  :: Data4D
   INTEGER(IntKi)                            :: ErrStatTmp
   
   IF (ALLOCATED(Data4D%DataIsSparse))   DEALLOCATE(Data4D%DataIsSparse,STAT=ErrStatTmp)
   IF (ALLOCATED(Data4D%LoadComponents)) DEALLOCATE(Data4D%LoadComponents,STAT=ErrStatTmp)

   IF (ALLOCATED(Data4D%DataSet))        DEALLOCATE(Data4D%DataSet,STAT=ErrStatTmp)
   IF (ALLOCATED(Data4D%DataMask))       DEALLOCATE(Data4D%DataMask,STAT=ErrStatTmp)
   IF (ALLOCATED(Data4D%WvFreq1))        DEALLOCATE(Data4D%WvFreq1,STAT=ErrStatTmp)
   IF (ALLOCATED(Data4D%WvFreq2))        DEALLOCATE(Data4D%WvFreq2,STAT=ErrStatTmp)
   IF (ALLOCATED(Data4D%WvDir1))         DEALLOCATE(Data4D%WvDir1,STAT=ErrStatTmp)
   IF (ALLOCATED(Data4D%WvDir2))         DEALLOCATE(Data4D%WvDir2,STAT=ErrStatTmp)
END SUBROUTINE Destroy_InitData4D


!----------------------------------------------------------------------------------------------------------------------------------

END MODULE WAMIT2
!**********************************************************************************************************************************
