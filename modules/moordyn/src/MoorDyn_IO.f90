!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020-2021 Alliance for Sustainable Energy, LLC
! Copyright (C) 2015-2019 Matthew Hall
!
!    This file is part of MoorDyn.
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
MODULE MoorDyn_IO

  ! This MODULE stores variables used for input and output and provides i/o subs

  USE                              NWTC_Library
  USE                              MoorDyn_Types
  IMPLICIT                         NONE


  PRIVATE

   INTEGER(IntKi), PARAMETER            :: wordy = 0   ! verbosity level. >1 = more console output

  INTEGER, PARAMETER :: nCoef = 30  ! maximum number of entries to allow in nonlinear coefficient lookup tables
  ! it would be nice if the above worked for everything, but I think it needs to also be matched in the Registry

  ! --------------------------- Output definitions -----------------------------------------

  ! The following are some definitions for use with the output options in MoorDyn.
  ! These are for the global output quantities specified by OutList, not line-specific outputs.
  ! Output definitions follow the structure described by the MD_OutParmType .
  ! Each output channel is described by the following fields:
  !  Name   - (string) what appears at the top of the output column
  !  Units  - (string) selected from UnitList (see below) based on index QType
  !  OType  - (int) the type of object the output is from. 1=line, 2=connect (0=invalid)
  !  ObjID  - (int) the ID number of the line or connect
  !  QType  - (int) the type of quantity to output.  0=tension, 1=x pos, etc.  see the parameters below
  !  NodeID - (int) the ID number of the node of the output quantity

  ! These are the "OTypes": 1=Line, 2=Connect, 3=Rod, 4=Body

  ! Indices for computing output channels:  - customized for the MD_OutParmType approach
  ! these are the "QTypes"
  INTEGER, PARAMETER             :: Time      =    0
  INTEGER, PARAMETER             :: PosX      =    1
  INTEGER, PARAMETER             :: PosY      =    2
  INTEGER, PARAMETER             :: PosZ      =    3
  INTEGER, PARAMETER             :: VelX      =    4
  INTEGER, PARAMETER             :: VelY      =    5
  INTEGER, PARAMETER             :: VelZ      =    6
  INTEGER, PARAMETER             :: AccX      =    7
  INTEGER, PARAMETER             :: AccY      =    8
  INTEGER, PARAMETER             :: AccZ      =    9
  INTEGER, PARAMETER             :: Ten       =   10
  INTEGER, PARAMETER             :: FX        =   11
  INTEGER, PARAMETER             :: FY        =   12
  INTEGER, PARAMETER             :: FZ        =   13
  INTEGER, PARAMETER             :: MX        =   14
  INTEGER, PARAMETER             :: MY        =   15
  INTEGER, PARAMETER             :: MZ        =   16
  INTEGER, PARAMETER             :: Pitch     =   17
  INTEGER, PARAMETER             :: Roll      =   18
  INTEGER, PARAMETER             :: Yaw       =   19
  INTEGER, PARAMETER             :: Sub       =   20

  ! List of units corresponding to the quantities parameters for QTypes
  CHARACTER(ChanLen), PARAMETER :: UnitList(0:20) =  (/ &
                               "(s)       ","(m)       ","(m)       ","(m)       ", &
                               "(m/s)     ","(m/s)     ","(m/s)     ", &
                               "(m/s2)    ","(m/s2)    ","(m/s2)    ", &
                               "(N)       ","(N)       ","(N)       ","(N)       ", &
                               "(Nm)      ","(Nm)      ","(Nm)      ", &
                               "(deg)     ","(deg)     ","(deg)     ","(frac)    "/)

  CHARACTER(28), PARAMETER  :: OutPFmt = "( I4, 3X,A 10,1 X, A10 )"   ! Output format parameter output list.
  CHARACTER(28), PARAMETER  :: OutSFmt = "ES10.3E2"


  ! output naming scheme is as
  ! examples:
  !  FairTen1, AnchTen1
  !  Con1pX
  !  Con3vY (connection 3, y velocity)
  !  L2N4pX (line 2, node 4, x position)

  ! ---------------------------------------------------------------------------------------------------------




  ! PUBLIC :: MDIO_ReadInput
   PUBLIC :: setupBathymetry
   PUBLIC :: getCoefficientOrCurve
   PUBLIC :: SplitByBars
   PUBLIC :: DecomposeString
   PUBLIC :: MDIO_OpenOutput
   PUBLIC :: MDIO_CloseOutput
   PUBLIC :: MDIO_ProcessOutList
   PUBLIC :: MDIO_WriteOutputs


CONTAINS


   SUBROUTINE setupBathymetry(inputString, defaultDepth, BathGrid, BathGrid_Xs, BathGrid_Ys, ErrStat3, ErrMsg3)
   ! SUBROUTINE getBathymetry(inputString, BathGrid, BathGrid_Xs, BathGrid_Ys, BathGrid_npoints, ErrStat3, ErrMsg3)

      CHARACTER(40),           INTENT(IN   )  :: inputString         ! string describing water depth or bathymetry filename
      REAL(ReKi),              INTENT(IN   )  :: defaultDepth        ! depth to use if inputString is empty
      REAL(DbKi), ALLOCATABLE, INTENT(INOUT)  :: BathGrid (:,:)
      REAL(DbKi), ALLOCATABLE, INTENT(INOUT)  :: BathGrid_Xs (:)
      REAL(DbKi), ALLOCATABLE, INTENT(INOUT)  :: BathGrid_Ys (:)
      INTEGER(IntKi),          INTENT(  OUT)  :: ErrStat3            ! Error status of the operation
      CHARACTER(*),            INTENT(  OUT)  :: ErrMsg3             ! Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)                   :: I
      INTEGER(IntKi)                   :: UnCoef   ! unit number for coefficient input file
      
      INTEGER(IntKi)                   :: ErrStat4
      CHARACTER(120)                   :: ErrMsg4         
      CHARACTER(120)                   :: Line2

      CHARACTER(20)                    :: nGridX_string  ! string to temporarily hold the nGridX string from Line2
      CHARACTER(20)                    :: nGridY_string  ! string to temporarily hold the nGridY string from Line3
      INTEGER(IntKi)                   :: nGridX         ! integer of the size of BathGrid_Xs
      INTEGER(IntKi)                   :: nGridY         ! integer of the size of BathGrid_Ys


      IF (LEN_TRIM(inputString) == 0) THEN
         ! If the input is empty (not provided), make the 1x1 bathymetry grid using the default depth
         ALLOCATE(BathGrid(1,1), STAT=ErrStat4)
         BathGrid(1,1) = REAL(defaultDepth,R8Ki)
         
         ALLOCATE(BathGrid_Xs(1), STAT=ErrStat4)
         BathGrid_Xs(1) = 0.0_DbKi
         
         ALLOCATE(BathGrid_Ys(1), STAT=ErrStat4)
         BathGrid_Ys(1) = 0.0_DbKi     
         
      ELSE IF (SCAN(inputString, "abcdfghijklmnopqrstuvwxyzABCDFGHIJKLMNOPQRSTUVWXYZ") == 0) THEN
         ! If the input does not have any of these string values, let's treat it as a number but store in a matrix
         ALLOCATE(BathGrid(1,1), STAT=ErrStat4)
         READ(inputString, *, IOSTAT=ErrStat4) BathGrid(1,1)
         
         ALLOCATE(BathGrid_Xs(1), STAT=ErrStat4)
         BathGrid_Xs(1) = 0.0_DbKi
         
         ALLOCATE(BathGrid_Ys(1), STAT=ErrStat4)
         BathGrid_Ys(1) = 0.0_DbKi

      ELSE ! otherwise interpret the input as a file name to load the bathymetry lookup data from
         CALL WrScr("   The depth input contains letters so will load a bathymetry file.")
         
         ! load lookup table data from file
         CALL GetNewUnit( UnCoef ) ! unit number for coefficient input file
         CALL OpenFInpFile( UnCoef, TRIM(inputString), ErrStat4, ErrMsg4 )
         cALL SetErrStat(ErrStat4, ErrMsg4, ErrStat3, ErrMsg3, 'MDIO_getBathymetry')

         READ(UnCoef,'(A)',IOSTAT=ErrStat4) Line2   ! skip the first title line
         READ(UnCoef,*,IOSTAT=ErrStat4) nGridX_string, nGridX  ! read in the second line as the number of x values in the BathGrid
         READ(UnCoef,*,IOSTAT=ErrStat4) nGridY_string, nGridY  ! read in the third line as the number of y values in the BathGrid

         ! Allocate the bathymetry matrix and associated grid x and y values
         ALLOCATE(BathGrid(nGridX, nGridY), STAT=ErrStat4)
         ALLOCATE(BathGrid_Xs(nGridX), STAT=ErrStat4)
         ALLOCATE(BathGrid_Ys(nGridY), STAT=ErrStat4)

         DO I = 1, nGridY+1  ! loop through each line in the rest of the bathymetry file

            READ(UnCoef,'(A)',IOSTAT=ErrStat4) Line2   ! read into a line and call it Line2
            IF (ErrStat4 > 0) EXIT

            IF (I==1) THEN    ! if it's the first line in the Bathymetry Grid, then it's a list of all the x values
               READ(Line2, *,IOSTAT=ErrStat4) BathGrid_Xs
            ELSE              ! if it's not the first line, then the first value is a y value and the rest are the depth values
               READ(Line2, *,IOSTAT=ErrStat4) BathGrid_Ys(I-1), BathGrid(I-1,:)
            ENDIF
         
         END DO

         IF (I < 2) THEN
            ErrStat3 = ErrID_Fatal
            ErrMsg3 = "Less than the minimum of 2 data lines found in file "//TRIM(inputString)
            CLOSE (UnCoef)
            RETURN
         ELSE 
            ! BathGrid_npoints = nGridX*nGridY       ! save the number of points in the grid
            CLOSE (UnCoef)
         END IF
      
      END IF

   END SUBROUTINE setupBathymetry
   

   ! read in stiffness/damping coefficient or load nonlinear data file if applicable
   SUBROUTINE getCoefficientOrCurve(inputString, LineProp_c, LineProp_npoints, LineProp_Xs, LineProp_Ys, ErrStat3, ErrMsg3)
   
      CHARACTER(40),    INTENT(IN   )  :: inputString
      REAL(DbKi),       INTENT(INOUT)  :: LineProp_c
      INTEGER(IntKi),   INTENT(  OUT)  :: LineProp_nPoints
      REAL(DbKi),       INTENT(  OUT)  :: LineProp_Xs (nCoef)
      REAL(DbKi),       INTENT(  OUT)  :: LineProp_Ys (nCoef)
      
      INTEGER(IntKi),   INTENT( OUT)   :: ErrStat3 ! Error status of the operation
      CHARACTER(*),     INTENT( OUT)   :: ErrMsg3  ! Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)                   :: I
      INTEGER(IntKi)                   :: UnCoef   ! unit number for coefficient input file
           
           
      INTEGER(IntKi)                   :: ErrStat4
      CHARACTER(120)                   :: ErrMsg4         
      CHARACTER(120)                   :: Line2   
      
           
      if (SCAN(inputString, "abcdfghijklmnopqrstuvwxyzABCDFGHIJKLMNOPQRSTUVWXYZ") == 0) then ! "eE" are exluded as they're used for scientific notation!
      
         ! "found NO letter in the line coefficient value so treating it as a number."
         READ(inputString, *, IOSTAT=ErrStat4) LineProp_c  ! convert the entry string into a real number
         LineProp_npoints = 0;
      
      else ! otherwise interpet the input as a file name to load stress-strain lookup data from
      
         CALL WrScr("found A letter in the line coefficient value so will try to load the filename.")
         
         LineProp_c = 0.0
         
         ! load lookup table data from file
        
         CALL GetNewUnit( UnCoef )
         CALL OpenFInpFile( UnCoef, TRIM(inputString), ErrStat4, ErrMsg4 )   ! add error handling?
         
         READ(UnCoef,'(A)',IOSTAT=ErrStat4) Line2   ! skip the first two lines (title, names, and units) then parse
         READ(UnCoef,'(A)',IOSTAT=ErrStat4) Line2
         READ(UnCoef,'(A)',IOSTAT=ErrStat4) Line2
            
         DO I = 1, nCoef
            
            READ(UnCoef,'(A)',IOSTAT=ErrStat4) Line2      !read into a line

            IF (ErrStat4 > 0) then
               CALL WrScr("Error while reading lookup table file")
               EXIT
            ELSE IF (ErrStat4 < 0) then
               CALL WrScr("Read "//trim(Int2LStr(I-1))//" data lines from lookup table file")
               EXIT
            ELSE
               READ(Line2,*,IOSTAT=ErrStat4) LineProp_Xs(I), LineProp_Ys(I)
            END IF 
         END DO
         
         if (I < 2) then
            ErrStat3 = ErrID_Fatal
            ErrMsg3  = "Less than the minimum of 2 data lines found in file "//TRIM(inputString)//" (first 3 lines are headers)."
            LineProp_npoints = 0
            Close (UnCoef)
            RETURN
         else
            LineProp_npoints = I-1
            Close (UnCoef)
         end if
      
      END IF
   
   END SUBROUTINE getCoefficientOrCurve
   
   
   ! Split a string into separate strings by the bar (|) symbol
   SUBROUTINE SplitByBars(instring, n, outstrings)
   
      CHARACTER(*),          INTENT(INOUT)  :: instring
      INTEGER(IntKi),        INTENT(  OUT)  :: n
      CHARACTER(40),         INTENT(INOUT)  :: outstrings(6)  ! array of output strings. Up to 6 strings can be read
      
      INTEGER :: pos1, pos2
 
      n = 0
      pos1=1
 
      DO
         pos2 = INDEX(instring(pos1:), "|")  ! find index of next comma
         IF (pos2 == 0) THEN                 ! if there isn't another comma, read the last entry and call it done (this could be the only entry if no commas)
            n = n + 1
            outstrings(n) = instring(pos1:)
            EXIT
         END IF
         n = n + 1
         if (n > 6) then
            CALL WrScr("ERROR - SplitByBars cannot do more than 6 entries")
         end if
         outstrings(n) = instring(pos1:pos1+pos2-2)
         pos1 = pos2+pos1
      END DO
      
   END SUBROUTINE SplitByBars


   ! Split a string into separate letter strings and integers. Letters are converted to uppercase.
   SUBROUTINE DecomposeString(outWord, let1, num1, let2, num2, let3)
   
      CHARACTER(*),          INTENT(INOUT)  :: outWord
      CHARACTER(25),         INTENT(  OUT)  :: let1
 !     INTEGER(IntKi),        INTENT(  OUT)  :: num1
      CHARACTER(25),         INTENT(  OUT)  :: num1
      CHARACTER(25),         INTENT(  OUT)  :: let2
      CHARACTER(25),         INTENT(  OUT)  :: num2
!      INTEGER(IntKi),        INTENT(  OUT)  :: num2
      CHARACTER(25),         INTENT(  OUT)  :: let3
   
!      INTEGER(IntKi)               :: I                                        ! Generic loop-counting index
      
!      CHARACTER(ChanLen)           :: OutListTmp                               ! A string to temporarily hold OutList(I), the name of each output channel
!      CHARACTER(ChanLen)           :: qVal                                     ! quantity type string to match to list of valid options
      
!      INTEGER                      :: oID                                      ! ID number of connect or line object
!      INTEGER                      :: nID                                      ! ID number of node object
      INTEGER                      :: i1 = 0                                   ! indices of start of numbers or letters in OutListTmp string, for parsing
      INTEGER                      :: i2 = 0
      INTEGER                      :: i3 = 0
      INTEGER                      :: i4 = 0

   
      CALL Conv2UC(outWord)       ! convert to all uppercase for string matching purposes

      ! start these strings as empty, and fill in only if used
      let1 = ''
      num1 = ''
      let2 = ''
      num2 = ''
      let3 = ''

      ! find indicies of changes in number-vs-letter in characters of outWord and split into segments accordingly
      
      i1 = scan( outWord , '1234567890' )              ! find index of first number in the string
      if (i1 > 0) then                                 ! if there is a number
         let1 = TRIM(outWord( 1:i1-1))
         i2 = i1+verify( outWord(i1+1:) , '1234567890' )  ! find starting index of second set of letters (if first character is a letter, i.e. i1>1), otherwise index of first letter
         if (i2 > i1) then                                ! if there is a second letter/word
            num1 = TRIM(outWord(i1:i2-1))
            i3 = i2+scan( outWord(i2+1:) , '1234567890' )    ! find starting index of second set of numbers <<<<
            if (i3 > i2) then                                ! if there is a second number
               let2 = TRIM(outWord(i2:i3-1))
               i4 = i3+verify( outWord(i3+1:) , '1234567890' )  ! third letter start
               if (i4 > i3) then                                ! if there is a third letter/word
                  num2 = TRIM(outWord(i3:i4-1))
                  let3 = TRIM(outWord(i4:   ))
               else
                  num2 = TRIM(outWord(i3:))
               end if
            else
               let2 = TRIM(outWord(i2:))
            end if
         else
            num1 = TRIM(outWord(i1:))
         end if
      else
         let1 = TRIM(outWord)
      end if
      
      
      !READ(outWord(i1:i2-1)) num1
      !READ(outWord(i3:i4-1)) num2
      
      ! print *, "Decomposed string ", outWord, " into:"
      ! print *, let1
      ! print *, num1
      ! print *, let2
      ! print *, num2
      ! print *, let3
      ! print *, "based on indices (i1-i4):"
      ! print *, i1
      ! print *, i2
      ! print *, i3
      ! print *, i4
   
   END SUBROUTINE DecomposeString
   


  ! ====================================================================================================
  SUBROUTINE MDIO_ProcessOutList(OutList, p, m, y, InitOut, ErrStat, ErrMsg )

  ! This routine processes the output channels requested by OutList, checking for validity and setting
  ! the p%OutParam structures (of type MD_OutParmType) for each valid output.
  ! It assumes the value p%NumOuts has been set beforehand, and sets the values of p%OutParam.


    IMPLICIT                        NONE

    ! Passed variables
    CHARACTER(ChanLen),        INTENT(IN)     :: OutList(:)                  ! The list of user-requested outputs
    TYPE(MD_ParameterType),    INTENT(INOUT)  :: p                           ! The module parameters
    TYPE(MD_MiscVarType),      INTENT(INOUT)  :: m
    TYPE(MD_OutputType),       INTENT(INOUT)  :: y                           ! Initial system outputs (outputs are not calculated; only the output mesh is initialized)
    TYPE(MD_InitOutputType),   INTENT(INOUT)  :: InitOut                     ! Output for initialization routine
    INTEGER(IntKi),            INTENT(OUT)    :: ErrStat                     ! The error status code
    CHARACTER(*),              INTENT(OUT)    :: ErrMsg                      ! The error message, if an error occurred

    ! Local variables
    INTEGER                      :: I                                        ! Generic loop-counting index
!    INTEGER                      :: J                                        ! Generic loop-counting index
!    INTEGER                      :: INDX                                     ! Index for valid arrays

    CHARACTER(ChanLen)           :: OutListTmp                               ! A string to temporarily hold OutList(I), the name of each output channel
    CHARACTER(ChanLen)           :: qVal                                     ! quantity type string to match to list of valid options

    INTEGER                      :: oID                                      ! ID number of connect or line object
    INTEGER                      :: nID                                      ! ID number of node object
!    INTEGER                      :: i1,i2,i3,i4                              ! indices of start of numbers or letters in OutListTmp string, for parsing
    
      CHARACTER(25)                 :: let1                ! strings used for splitting and parsing identifiers
      CHARACTER(25)                 :: num1
      CHARACTER(25)                 :: let2
      CHARACTER(25)                 :: num2
      CHARACTER(25)                 :: let3
      
    INTEGER(IntKi)                            :: LineNumOuts                 ! number of entries in LineWrOutput for each line
    INTEGER(IntKi)                            :: RodNumOuts                  !   same for Rods
      

    ! see the top of the module for info on the output labelling types

    ! Initialize values
    ErrStat = ErrID_None
    ErrMsg = ""


    ALLOCATE ( p%OutParam(1:p%NumOuts) , STAT=ErrStat )   ! note: I'm skipping the time output entry at index 0 for simplicity
    IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = "Error allocating memory for the MoorDyn OutParam array."
      RETURN
    ELSE
      ErrStat = ErrID_None
    ENDIF


    ! Set index, name, and units for the time output channel: ! note: I'm skipping the time output entry at index 0
    !p%OutParam(0)%Indx  = Time
    !p%OutParam(0)%Name  = "Time"    ! OutParam(0) is the time channel by default.
    !p%OutParam(0)%Units = "(s)"
    !p%OutParam(0)%SignM = 1


    ! Set index, name, and units for all of the output channels.
    ! If a selected output channel is not valid set ErrStat = ErrID_Warn.


    ! go through list of requested output names and process (this is a bit of a mess)

    DO I = 1,p%NumOuts

      OutListTmp          = OutList(I)  ! current requested output name
      
      call DecomposeString(OutListTmp, let1, num1, let2, num2, let3)
      
      
      
      !p%OutParam(I)%Name  = OutListTmp
      CALL Conv2UC(OutListTmp)       ! convert to all uppercase for string matching purposes

   !   ! find indicies of changes in number-vs-letter in characters of OutListTmp
   !   i1 = scan( OutListTmp , '1234567890' )              ! first number in the string
   !   i2 = i1+verify( OutListTmp(i1+1:) , '1234567890' )  ! second letter start (assuming first character is a letter, i.e. i1>1)
   !   i3 = i2+scan( OutListTmp(i2+1:) , '1234567890' )    ! second number start
   !   i4 = i3+verify( OutListTmp(i3+1:) , '1234567890' )  ! third letter start
   
      ! error check
   !   IF (i1 <= 1) THEN
   !      CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
   !      CALL WrScr('Warning: invalid output specifier '//trim(OutListTmp)//'.  Starting character must be C or L.')
   !      CYCLE    ! <<<<<<<<<<< check correct usage
   !   END IF

        p%OutParam(I)%Name = OutListTmp  ! label channel with whatever name was inputted, for now


      ! figure out what type of output it is and process accordingly

      ! fairlead tension case (updated) <<<<<<<<<<<<<<<<<<<<<<<<<<< these are not currently working - need new way to find ObjID
      IF (let1 == 'FAIRTEN') THEN
        p%OutParam(I)%OType = 1                                     ! line object type
        p%OutParam(I)%QType = Ten                                   ! tension quantity type
        p%OutParam(I)%Units = UnitList(Ten)                         ! set units according to QType
        READ (num1,*) oID                                           ! this is the line number
        p%OutParam(I)%ObjID  = oID                                  ! record the ID of the line
        p%OutParam(I)%NodeID = m%LineList(oID)%N                    ! specify node N (end B, fairlead)
        ! >>> should check validity of ObjID and NodeID <<<
        
      ! achor tension case
      ELSE IF (let1 == 'ANCHTEN') THEN
        p%OutParam(I)%OType = 1                                     ! line object type
        p%OutParam(I)%QType = Ten                                   ! tension quantity type
        p%OutParam(I)%Units = UnitList(Ten)                         ! set units according to QType
        READ (num1,*) oID                                           ! this is the line number
        p%OutParam(I)%ObjID  = oID                                  ! record the ID of the line
        p%OutParam(I)%NodeID = 0                                    ! specify node 0 (end A, anchor)

      ! more general case
      ELSE

        ! what object type?
        
        ! Line case                               
        IF (let1(1:1) == 'L') THEN      ! Look for L?N?xxxx
          p%OutParam(I)%OType = 1                ! Line object type
          ! for now we'll just assume the next character(s) are "n" to represent node number or "s" to represent segment number
          IF (num2/=" ") THEN
              READ (num2,*) nID                      ! node or segment ID
              p%OutParam(I)%NodeID = nID
          ELSE
              CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
              CALL WrScr('Warning: invalid output specifier '//trim(OutListTmp)//'. Line ID or Node ID missing.')
              CYCLE
          END IF
          qVal = let3                            ! quantity type string
        
        ! Connect case                            
        ELSE IF (let1(1:1) == 'C') THEN    ! Look for C?xxx or Con?xxx
          p%OutParam(I)%OType = 2                ! Connect object type
          qVal = let2                            ! quantity type string
          
        ! Rod case                            
        ELSE IF (let1(1:1) == 'R') THEN    ! Look for R?xxx or Rod?xxx
          p%OutParam(I)%OType = 3                ! Rod object type
          IF (LEN_TRIM(let3)== 0) THEN           ! No third character cluster indicates this is a whole-rod channel
            p%OutParam(I)%NodeID = 0
            qVal = let2                          ! quantity type string
          ELSE IF (num2/=" ") THEN
            READ (num2,*) nID                    ! rod node ID
            p%OutParam(I)%NodeID = nID
            qVal = let3                          ! quantity type string
          ELSE
            CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
            CALL WrScr('Warning: invalid output specifier '//trim(OutListTmp)//'.  Rod ID or Node ID missing.')
            CYCLE
          END IF
          
        ! Body case                            
        ELSE IF (Let1(1:1) == 'B') THEN    ! Look for B?xxx or Body?xxx
          p%OutParam(I)%OType = 4                ! Body object type
          qVal = let2                            ! quantity type string

        ! should do fairlead option also!

        ! error
        ELSE
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
          CALL WrScr('Warning: invalid output specifier '//trim(OutListTmp)//'.  Must start with L, C, R, or B')
          CYCLE
        END IF

        ! object number
        IF (num1/=" ") THEN
            READ (num1,*) oID
            p%OutParam(I)%ObjID = oID                ! line or connect ID number
        ELSE
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
          CALL WrScr('Warning: invalid output specifier '//trim(OutListTmp)//'.  Object ID missing.')
          CYCLE
        END IF

        ! which kind of quantity?
        IF (qVal == 'PX') THEN
          p%OutParam(I)%QType = PosX
          p%OutParam(I)%Units = UnitList(PosX)
        ELSE IF (qVal == 'PY') THEN
          p%OutParam(I)%QType = PosY
          p%OutParam(I)%Units = UnitList(PosY)
        ELSE IF (qVal == 'PZ') THEN
          p%OutParam(I)%QType = PosZ
          p%OutParam(I)%Units = UnitList(PosZ)
        ELSE IF (qVal == 'VX') THEN
          p%OutParam(I)%QType = VelX
          p%OutParam(I)%Units = UnitList(VelX)
        ELSE IF (qVal == 'VY') THEN
          p%OutParam(I)%QType = VelY
          p%OutParam(I)%Units = UnitList(VelY)
        ELSE IF (qVal == 'VZ') THEN
          p%OutParam(I)%QType = VelZ
          p%OutParam(I)%Units = UnitList(VelZ)
        ELSE IF (qVal == 'AX') THEN
          p%OutParam(I)%QType = AccX
          p%OutParam(I)%Units = UnitList(AccX)
        ELSE IF (qVal == 'AY') THEN   ! fixed typo Nov 24
          p%OutParam(I)%QType = AccY
          p%OutParam(I)%Units = UnitList(AccY)
        ELSE IF (qVal == 'AZ') THEN
          p%OutParam(I)%QType = AccZ
          p%OutParam(I)%Units = UnitList(AccZ)
        ELSE IF ((qVal == 'T') .or. (qVal == 'TEN')) THEN
          p%OutParam(I)%QType = Ten
          p%OutParam(I)%Units = UnitList(Ten)
        ELSE IF (qVal == 'FX') THEN
          p%OutParam(I)%QType = FX
          p%OutParam(I)%Units = UnitList(FX)
        ELSE IF (qVal == 'FY') THEN
          p%OutParam(I)%QType = FY
          p%OutParam(I)%Units = UnitList(FY)
        ELSE IF (qVal == 'FZ') THEN
          p%OutParam(I)%QType = FZ
          p%OutParam(I)%Units = UnitList(FZ)   ! <<<< should add moments as well <<<<
        ELSE IF (qVal == 'ROLL') THEN
          p%OutParam(I)%QType = Roll
          p%OutParam(I)%Units = UnitList(Roll)
        ELSE IF (qVal == 'PITCH') THEN
          p%OutParam(I)%QType = Pitch
          p%OutParam(I)%Units = UnitList(Pitch)
        ELSE IF (qVal == 'YAW') THEN
          p%OutParam(I)%QType = Yaw
          p%OutParam(I)%Units = UnitList(Yaw)
        ELSE IF (qVal == 'SUB') THEN
          p%OutParam(I)%QType = Sub
          p%OutParam(I)%Units = UnitList(Sub)
        ELSE
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
          CALL WrScr('Warning: invalid output specifier '//trim(OutListTmp)//'.  Quantity type not recognized.')
          CONTINUE
        END IF

      END IF

      ! also check whether each object index and node index (if applicable) is in range
      
      IF (p%OutParam(I)%OType==1) THEN              ! Line
        IF (p%OutParam(I)%ObjID > p%NLines) THEN
          CALL WrScr('Warning: output Line index excedes number of Lines in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF
        IF (p%OutParam(I)%NodeID > m%LineList(p%OutParam(I)%ObjID)%N) THEN
          CALL WrScr('Warning: output node index excedes number of nodes in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        ELSE IF (p%OutParam(I)%NodeID < 0) THEN
          CALL WrScr('Warning: output node index is less than zero in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF
        
      ELSE IF (p%OutParam(I)%OType==2) THEN         ! Connect
        IF (p%OutParam(I)%ObjID > p%NConnects) THEN
          CALL WrScr('Warning: output Connect index excedes number of Connects in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF
        
      ELSE IF (p%OutParam(I)%OType==3) THEN         ! Rod
        IF (p%OutParam(I)%ObjID > p%NRods) THEN
          CALL WrScr('Warning: output Rod index excedes number of Rods in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF
        IF (p%OutParam(I)%NodeID > m%RodList(p%OutParam(I)%ObjID)%N) THEN
          CALL WrScr('Warning: output node index excedes number of nodes in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        ELSE IF (p%OutParam(I)%NodeID < 0) THEN
          CALL WrScr('Warning: output node index is less than zero in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF
        
      ELSE IF (p%OutParam(I)%OType==4) THEN         ! Body
        IF (p%OutParam(I)%ObjID > p%NBodies) THEN
          CALL WrScr('Warning: output Body index excedes number of Bodies in requested output '//trim(OutListTmp)//'.')
          CALL DenoteInvalidOutput(p%OutParam(I)) ! flag as invalid
        END IF

      END IF

      ! is the reverse sign functionality necessary?
      !      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !      !   channel name with a "-", "_", "m", or "M" character indicating "minus".

     END DO  ! I ... looping through OutList


!!   ! Allocate MDWrOutput which is used to store a time step's worth of output channels, prior to writing to a file.
!    ALLOCATE( MDWrOutput( p%NumOuts),  STAT = ErrStat )
!    IF ( ErrStat /= ErrID_None ) THEN
!      ErrMsg  = ' Error allocating space for MDWrOutput array.'
!      ErrStat = ErrID_Fatal
!      RETURN
!    END IF


      ! Allocate MDWrOuput2 which is used to store a time step's worth of output data for each line, just making it really big for now <<<<<<<<<<<<<<
      ! <<<<<<<<<<< should do this for each line instead.
   !   ALLOCATE( LineWriteOutputs( 200),  STAT = ErrStat )
   !   IF ( ErrStat /= ErrID_None ) THEN
   !      ErrMsg  = ' Error allocating space for LineWriteOutputs array.'
   !      ErrStat = ErrID_Fatal
   !      RETURN
   !   END IF

      !Allocate WriteOuput
      ALLOCATE(        y%WriteOutput(  p%NumOuts), &
              InitOut%WriteOutputHdr(p%NumOuts), &
              InitOut%WriteOutputUnt(p%NumOuts),  STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for y%WriteOutput array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF

      ! allocate output array in each Line
      DO I=1,p%NLines
      
      
         ! calculate number of output entries (excluding time) to write for this line
         LineNumOuts = 3*(m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(2:6)) &
                       + (m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(7:9)) &
                             + m%LineList(I)%N*SUM(m%LineList(I)%OutFlagList(10:18))
   
         ALLOCATE(m%LineList(I)%LineWrOutput( 1 + LineNumOuts), STAT = ErrStat)  
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Error allocating space for a LineWrOutput array'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
      END DO  ! I
      
      ! allocate output array in each Rod
      DO I=1,p%NRods
      
         ! calculate number of output entries (excluding time) to write for this Rod
         RodNumOuts = 3*(m%RodList(I)%N + 1)*SUM(m%RodList(I)%OutFlagList(2:9)) &
                       + (m%RodList(I)%N + 1)*SUM(m%RodList(I)%OutFlagList(10:11)) &
                             + m%RodList(I)%N*SUM(m%RodList(I)%OutFlagList(12:18))
      
         ALLOCATE(m%RodList(I)%RodWrOutput( 1 + RodNumOuts), STAT = ErrStat)  
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Error allocating space for a RodWrOutput array'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
      END DO  ! I

      !print *, "y%WriteOutput allocated to size ", size(y%WriteOutput)

      ! These variables are to help follow the framework template, but the data in them is simply a copy of data
      ! already available in the OutParam data structure
      !  ALLOCATE ( InitOut%WriteOutputHdr(p%NumOuts+p%OutAllint*p%OutAllDims), STAT = ErrStat )
      !  ALLOCATE ( InitOut%WriteOutputUnt(p%NumOuts+p%OutAllint*p%OutAllDims), STAT = ErrStat )

      DO I = 1,p%NumOuts
         InitOut%WriteOutputHdr(I) = p%OutParam(I)%Name
         InitOut%WriteOutputUnt(I) = p%OutParam(I)%Units
      END DO


   CONTAINS

      SUBROUTINE DenoteInvalidOutput( OutParm )
         TYPE(MD_OutParmType), INTENT (INOUT)  :: OutParm

         OutParm%OType = 0  ! flag as invalid
         OutParm%Name = 'Invalid'
         OutParm%Units = ' - '

      END SUBROUTINE DenoteInvalidOutput

   END SUBROUTINE MDIO_ProcessOutList
   !----------------------------------------------------------------------------------------============





   !----------------------------------------------------------------------------------------============
   SUBROUTINE MDIO_OpenOutput( MD_ProgDesc, p, m, InitOut, ErrStat, ErrMsg )
   !----------------------------------------------------------------------------------------------------
      TYPE(ProgDesc),                INTENT( IN    ) :: MD_ProgDesc
      TYPE(MD_ParameterType),        INTENT( INOUT ) :: p
      TYPE(MD_MiscVarType),          INTENT( INOUT ) :: m
      TYPE(MD_InitOutPutType ),      INTENT( IN    ) :: InitOut              !
      INTEGER,                       INTENT(   OUT ) :: ErrStat              ! a non-zero value indicates an error occurred
      CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      INTEGER                                        :: I                    ! Generic loop counter
      INTEGER                                        :: J                    ! Generic loop counter
      CHARACTER(1024)                                :: OutFileName          ! The name of the output file  including the full path.
!      INTEGER                                        :: L                    ! counter for index in LineWrOutput
      INTEGER                                        :: LineNumOuts          ! number of entries in LineWrOutput for each line
      INTEGER                                        :: RodNumOuts           ! for Rods ... redundant <<<
      CHARACTER(200)                                 :: Frmt                 ! a string to hold a format statement
      INTEGER                                        :: ErrStat2


      ErrStat = ErrID_None
      ErrMsg  = ""

      p%Delim = ' '  ! for now

      !-------------------------------------------------------------------------------------------------
      ! Open the output file, if necessary, and write the header
      !-------------------------------------------------------------------------------------------------

      IF ( ALLOCATED( p%OutParam ) .AND. p%NumOuts > 0 ) THEN           ! Output has been requested so let's open an output file

         ! Open the file for output
         OutFileName = TRIM(p%RootName)//'.out'
         CALL GetNewUnit( p%MDUnOut )

         CALL OpenFOutFile ( p%MDUnOut, OutFileName, ErrStat, ErrMsg )
         IF ( ErrStat > ErrID_None ) THEN
            ErrMsg = ' Error opening MoorDyn-level output file: '//TRIM(ErrMsg)
            ErrStat = ErrID_Fatal
            RETURN
         END IF


         ! Write empty header lines (to match standard OF out file formats used by all other modules)
         write (p%MDUnOut,'(/,A/)', IOSTAT=ErrStat)  'These predictions were generated by '//&
                        TRIM(MD_ProgDesc%Name)//' ('//TRIM(MD_ProgDesc%Ver)//TRIM(MD_ProgDesc%Date)//&
                        ') on '//CurDate()//' at '//CurTime()//'.'
         write(p%MDUnOut,*) ''
         write(p%MDUnOut,*) ''
         write(p%MDUnOut,*) ''

         !Write the names of the output parameters:

         Frmt = '(A10,'//TRIM(Int2LStr(p%NumOuts))//'(A1,A15))'

         WRITE(p%MDUnOut,Frmt, IOSTAT=ErrStat2)  TRIM( 'Time' ), ( p%Delim, TRIM( p%OutParam(I)%Name), I=1,p%NumOuts )

         WRITE(p%MDUnOut,Frmt)  TRIM( '(s)' ), ( p%Delim, TRIM( p%OutParam(I)%Units ), I=1,p%NumOuts )

 !     ELSE  ! if no outputs requested

 !        call wrscr('note, MDIO_OpenOutput thinks that no outputs have been requested.')

      END IF

      !--------------------------------------------------------------------------
      !                    now do the same for line output files
      !--------------------------------------------------------------------------

      !! allocate UnLineOuts
      !ALLOCATE(UnLineOuts(p%NLines))  ! should add error checking

      DO I = 1,p%NLines

         
         IF (m%LineList(I)%OutFlagList(1) == 1) THEN   ! only proceed if the line is flagged to output a file
           
            ! Open the file for output
            OutFileName = TRIM(p%RootName)//'.Line'//TRIM(Int2LStr(I))//'.out'
            CALL GetNewUnit( m%LineList(I)%LineUnOut )

            CALL OpenFOutFile ( m%LineList(I)%LineUnOut, OutFileName, ErrStat, ErrMsg )
            IF ( ErrStat > ErrID_None ) THEN
               ErrMsg = ' Error opening Line output file '//TRIM(ErrMsg)
               ErrStat = ErrID_Fatal
               RETURN
            END IF

                        
            ! calculate number of output entries (excluding time) to write for this line
            LineNumOuts = 3*(m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(2:6)) &
                          + (m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(7:9)) &
                                + m%LineList(I)%N*SUM(m%LineList(I)%OutFlagList(10:18))
                                  
            if (wordy > 2) PRINT *, LineNumOuts, " output channels"

            Frmt = '(A10,'//TRIM(Int2LStr(1 + LineNumOuts))//'(A1,A15))'   ! should evenutally use user specified format?
            !Frmt = '(A10,'//TRIM(Int2LStr(3+3*m%LineList(I)%N))//'(A1,A15))'
            
            ! Write the names of the output parameters:  (these use "implied DO" loops)

            WRITE(m%LineList(I)%LineUnOut,'(A10)', advance='no', IOSTAT=ErrStat2)  TRIM( 'Time' )
            IF (m%LineList(I)%OutFlagList(2) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'px', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'py', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'pz', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(3) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'vx', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'vy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'vz', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(4) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Ux', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Uy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Uz', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(5) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Dx', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Dy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Dz', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(6) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'bx', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'by', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'bz', J=0,(m%LineList(I)%N) )
            END IF
            
            IF (m%LineList(I)%OutFlagList(7) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Wz', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(8) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Kurv', J=0,(m%LineList(I)%N) )
            END IF
            
            IF (m%LineList(I)%OutFlagList(10) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Ten', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(11) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Dmp', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(12) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Str', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(13) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'SRt', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(14)== 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Lst', J=1,(m%LineList(I)%N) )
            END IF
            
            WRITE(m%LineList(I)%LineUnOut,'(A1)', IOSTAT=ErrStat2) ' '  ! make line break at the end
            
            
            ! Now write the units line

            WRITE(m%LineList(I)%LineUnOut,'(A10)', advance='no', IOSTAT=ErrStat2)  TRIM( '(s)' )
            IF (m%LineList(I)%OutFlagList(2) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m)', p%Delim, '(m)', p%Delim, '(m)', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(3) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m/s)', p%Delim, '(m/s)', p%Delim, '(m/s)', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(4) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m/s)', p%Delim, '(m/s)', p%Delim, '(m/s)', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(5) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', p%Delim, '(N)', p%Delim, '(N)', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(6) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((3+3*m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', p%Delim, '(N)', p%Delim, '(N)', J=0,(m%LineList(I)%N) )
            END IF
            
            IF (m%LineList(I)%OutFlagList(7) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(Nup)', J=0,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(8) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(1/m)', J=0,(m%LineList(I)%N) )
            END IF
            
            IF (m%LineList(I)%OutFlagList(10) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(11) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(12) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(-)', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(13) == 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(1/s)', J=1,(m%LineList(I)%N) )
            END IF
            IF (m%LineList(I)%OutFlagList(14)== 1) THEN
               WRITE(m%LineList(I)%LineUnOut,'('//TRIM(Int2LStr((m%LineList(I)%N)))//'(A1,A10))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m)', J=1,(m%LineList(I)%N) )
            END IF
            
            WRITE(m%LineList(I)%LineUnOut,'(A1)', IOSTAT=ErrStat2) ' '  ! make line break at the end
            
         END IF  ! if line is flagged for output file
         
      END DO ! I - line number




      !--------------------------------------------------------------------------
      !                    now do the same for rod output files 
      !--------------------------------------------------------------------------

      !! allocate UnLineOuts
      !ALLOCATE(UnLineOuts(p%NLines))  ! should add error checking

      DO I = 1,p%NRods

         
         IF (m%RodList(I)%OutFlagList(1) == 1) THEN   ! only proceed if the Rod is flagged to output a file
           
            ! Open the file for output
            OutFileName = TRIM(p%RootName)//'.Rod'//TRIM(Int2LStr(I))//'.out'
            CALL GetNewUnit( m%RodList(I)%RodUnOut )

            CALL OpenFOutFile ( m%RodList(I)%RodUnOut, OutFileName, ErrStat, ErrMsg )
            IF ( ErrStat > ErrID_None ) THEN
               ErrMsg = ' Error opening Rod output file '//TRIM(ErrMsg)
               ErrStat = ErrID_Fatal
               RETURN
            END IF

                        
            ! calculate number of output entries (excluding time) to write for this Rod
            RodNumOuts = 3*(m%RodList(I)%N + 1)*SUM(m%RodList(I)%OutFlagList(2:9)) &
                          + (m%RodList(I)%N + 1)*SUM(m%RodList(I)%OutFlagList(10:11)) &
                                + m%RodList(I)%N*SUM(m%RodList(I)%OutFlagList(12:18))
                                  
            if (wordy > 2) PRINT *, RodNumOuts, " output channels"

            Frmt = '(A10,'//TRIM(Int2LStr(1 + RodNumOuts))//'(A1,A15))'   ! should evenutally use user specified format?
            !Frmt = '(A10,'//TRIM(Int2LStr(3+3*m%RodList(I)%N))//'(A1,A15))'
            
            ! >>> should functionalize the below <<<
            
            
            ! Write the names of the output parameters:  (these use "implied DO" loops)

            WRITE(m%RodList(I)%RodUnOut,'(A10)', advance='no', IOSTAT=ErrStat2)  TRIM( 'Time' )
            IF (m%RodList(I)%OutFlagList(2) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'px', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'py', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'pz', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(3) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'vx', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'vy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'vz', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(4) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Ux', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Uy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Uz', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(5) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Box', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Boy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Boz', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(6) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Dx', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Dy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Dz', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(7) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Fix', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Fiy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Fiz', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(8) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Pdx', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Pdy', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Pdz', J=0,(m%RodList(I)%N) )
            END IF            
            IF (m%RodList(I)%OutFlagList(9) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'bx', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'by', p%Delim, 'Node'//TRIM(Int2Lstr(J))//'bz', J=0,(m%RodList(I)%N) )
            END IF
            
            IF (m%RodList(I)%OutFlagList(10) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Wz', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(11) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Node'//TRIM(Int2Lstr(J))//'Kurv', J=0,(m%RodList(I)%N) )
            END IF
            
            IF (m%RodList(I)%OutFlagList(12) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Ten', J=1,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(13) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Dmp', J=1,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(14) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'Str', J=1,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(15) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, 'Seg'//TRIM(Int2Lstr(J))//'SRt', J=1,(m%RodList(I)%N) )
            END IF
            
            WRITE(m%RodList(I)%RodUnOut,'(A1)', IOSTAT=ErrStat2) ' '  ! make line break at the end
            
            
            ! Now write the units line

            WRITE(m%RodList(I)%RodUnOut,'(A10)', advance='no', IOSTAT=ErrStat2)  TRIM( '(s)' )
            IF (m%RodList(I)%OutFlagList(2) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m)', p%Delim, '(m)', p%Delim, '(m)', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(3) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m/s)', p%Delim, '(m/s)', p%Delim, '(m/s)', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(4) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(m/s)', p%Delim, '(m/s)', p%Delim, '(m/s)', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(5) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', p%Delim, '(N)', p%Delim, '(N)', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(6) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', p%Delim, '(N)', p%Delim, '(N)', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(7) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', p%Delim, '(N)', p%Delim, '(N)', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(8) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', p%Delim, '(N)', p%Delim, '(N)', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(9) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((3+3*m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', p%Delim, '(N)', p%Delim, '(N)', J=0,(m%RodList(I)%N) )
            END IF
            
            IF (m%RodList(I)%OutFlagList(10) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(Nup)', J=0,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(11) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(1/m)', J=0,(m%RodList(I)%N) )
            END IF
            
            IF (m%RodList(I)%OutFlagList(12) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', J=1,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(13) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(N)', J=1,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(14) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(-)', J=1,(m%RodList(I)%N) )
            END IF
            IF (m%RodList(I)%OutFlagList(15) == 1) THEN
               WRITE(m%RodList(I)%RodUnOut,'('//TRIM(Int2LStr((m%RodList(I)%N)))//'(A1,A15))', advance='no', IOSTAT=ErrStat2) &
                  ( p%Delim, '(1/s)', J=1,(m%RodList(I)%N) )
            END IF
            
            WRITE(m%RodList(I)%RodUnOut,'(A1)', IOSTAT=ErrStat2) ' '  ! make Rod break at the end
            
         END IF  ! if rod is flagged for output file
         
      END DO ! I - rod number

      ! need to fix error handling in this sub

   END SUBROUTINE MDIO_OpenOutput
   !----------------------------------------------------------------------------------------============


   !----------------------------------------------------------------------------------------============
   SUBROUTINE MDIO_CloseOutput ( p, m, ErrStat, ErrMsg )
      ! This function cleans up after running the MoorDyn output module.
      ! It closes the output files and releases memory.

      TYPE(MD_ParameterType),       INTENT( INOUT )  :: p                    ! data for this instance of the floating platform module
      TYPE(MD_MiscVarType),         INTENT( INOUT )  :: m                    ! data for this instance of the floating platform module
      INTEGER,                      INTENT(   OUT )  :: ErrStat              ! a non-zero value indicates an error occurred
      CHARACTER(*),                 INTENT(   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)       :: I  ! generic counter


      ErrStat = 0
      ErrMsg  = ""


!FIXME: make sure thes are actually open before trying to close them. Segfault will occur otherwise!!!!
!  This bug can be triggered by an early failure of the parsing routines, before these files were ever opened
!  which returns MD to OpenFAST as ErrID_Fatal, then OpenFAST calls MD_End, which calls this.

      ! close main MoorDyn output file
      if (p%MDUnOut > 0) then
         CLOSE( p%MDUnOut, IOSTAT = ErrStat )
         IF ( ErrStat /= 0 ) THEN
            ErrMsg = 'Error closing output file'
         END IF
      end if 
      
      ! close individual rod output files
      DO I=1,p%NRods
         if (allocated(m%RodList)) then
            if (m%RodList(I)%RodUnOut > 0) then
               CLOSE( m%RodList(I)%RodUnOut, IOSTAT = ErrStat )
               IF ( ErrStat /= 0 ) THEN
                  ErrMsg = 'Error closing rod output file'
               END IF
            end if 
         end if 
      END DO
      
      ! close individual line output files
      DO I=1,p%NLines
         if (allocated(m%LineList)) then
            if (m%LineList(I)%LineUnOut > 0) then
               CLOSE( m%LineList(I)%LineUnOut, IOSTAT = ErrStat )
               IF ( ErrStat /= 0 ) THEN
                  ErrMsg = 'Error closing line output file'
               END IF
            end if 
         end if
      END DO

      ! deallocate output arrays
      IF (ALLOCATED(m%MDWrOutput)) THEN
         DEALLOCATE(m%MDWrOutput)
      ENDIF
      DO I=1,p%NLines
         IF (ALLOCATED(m%LineList(I)%LineWrOutput)) THEN
            DEALLOCATE(m%LineList(I)%LineWrOutput)       ! this may be unnecessary and handled by Line destructor
         ENDIF
      END DO

   END SUBROUTINE MDIO_CloseOutput
   !----------------------------------------------------------------------------------------============


   !----------------------------------------------------------------------------------------============
   SUBROUTINE MDIO_WriteOutputs( Time, p, m, y, ErrStat, ErrMsg )
      ! This subroutine gathers the output data defined by the OutParams list and
      ! writes it to the output file opened in MDIO_OutInit()

      REAL(DbKi),                   INTENT( IN    ) :: Time                 ! Time for this output
      TYPE(MD_ParameterType),       INTENT( IN    ) :: p                    ! MoorDyn module's parameter data
      TYPE(MD_OutputType),          INTENT( INOUT ) :: y                    ! INTENT( OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)
      TYPE(MD_MiscVarType),         INTENT( INOUT ) :: m                    ! MoorDyn module's m data
      INTEGER,                      INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs
      CHARACTER(*),                 INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      INTEGER                                :: I                           ! Generic loop counter
      INTEGER                                :: J                           ! Generic loop counter
      INTEGER                                :: K                           ! Generic loop counter
      INTEGER                                :: L                           ! counter for index in LineWrOutput
      INTEGER                                :: LineNumOuts                 ! number of entries in LineWrOutput for each line
      INTEGER                                :: RodNumOuts                  !   same for Rods
      CHARACTER(200)                         :: Frmt                        ! a string to hold a format statement


      IF ( .NOT. ALLOCATED( p%OutParam ) .OR. p%MDUnOut < 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' To write outputs for MoorDyn there must be a valid file ID and OutParam must be allocated.'
         RETURN
      ELSE
         ErrStat = ErrID_None
         ErrMsg  = ''
      END IF
      
      ! -------------------------------- main output file --------------------------------
      
      if ( p%NumOuts > 0_IntKi ) then  

         ! gather the required output quantities (INCOMPLETE!)
         DO I = 1,p%NumOuts


            IF (p%OutParam(I)%OType == 1) THEN  ! if dealing with a Line output

               SELECT CASE (p%OutParam(I)%QType)
                  CASE (PosX)
                    y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%r(1,p%OutParam(I)%NodeID)  ! x position
                  CASE (PosY)
                    y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%r(2,p%OutParam(I)%NodeID) ! y position
                  CASE (PosZ)
                    y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%r(3,p%OutParam(I)%NodeID) ! z position
                  CASE (VelX)
                    y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%rd(1,p%OutParam(I)%NodeID) ! x velocity
                  CASE (VelY)
                    y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%rd(2,p%OutParam(I)%NodeID) ! y velocity
                  CASE (VelZ)
                    y%WriteOutput(I) = m%LineList(p%OutParam(I)%ObjID)%rd(3,p%OutParam(I)%NodeID) ! z velocity
                  CASE (Ten)
                    y%WriteOutput(I) = Line_GetNodeTen(m%LineList(p%OutParam(I)%ObjID), p%OutParam(I)%NodeID, p)  ! this is actually the segment tension ( 1 < NodeID < N )  Should deal with properly!
                    
                  CASE DEFAULT
                    y%WriteOutput(I) = 0.0_ReKi
                    ErrStat = ErrID_Warn
                    ErrMsg = ' Unsupported output quantity '//TRIM(p%OutParam(I)%Name)//' requested from Line '//TRIM(Num2Lstr(p%OutParam(I)%ObjID))//'.'
               END SELECT

            ELSE IF (p%OutParam(I)%OType == 2) THEN  ! if dealing with a Connect output
               SELECT CASE (p%OutParam(I)%QType)
                  CASE (PosX)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%r(1)  ! x position
                  CASE (PosY)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%r(2) ! y position
                  CASE (PosZ)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%r(3) ! z position
                  CASE (VelX)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%rd(1) ! x velocity
                  CASE (VelY)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%rd(2) ! y velocity
                  CASE (VelZ)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%rd(3) ! z velocity
                  CASE (AccX)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%a(1) ! x acceleration
                  CASE (AccY)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%a(2) ! y acceleration
                  CASE (AccZ)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%a(3) ! z acceleration
                  CASE (Ten)
                     y%WriteOutput(I) = TwoNorm(m%ConnectList(p%OutParam(I)%ObjID)%Fnet)  ! total force magnitude on a connect (used eg. for fairlead and anchor tensions)
                  CASE (FX)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%Fnet(1)  ! total force in x - added Nov 24
                  CASE (FY)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%Fnet(2)  ! total force in y
                  CASE (FZ)
                     y%WriteOutput(I) = m%ConnectList(p%OutParam(I)%ObjID)%Fnet(3)  ! total force in z
                  CASE DEFAULT
                     y%WriteOutput(I) = 0.0_ReKi
                     ErrStat = ErrID_Warn
                     ErrMsg = ' Unsupported output quantity '//TRIM(p%OutParam(I)%Name)//' requested from Connection '//TRIM(Num2Lstr(p%OutParam(I)%ObjID))//'.'
               END SELECT

            ELSE IF (p%OutParam(I)%OType == 3) THEN  ! if dealing with a Rod output

               SELECT CASE (p%OutParam(I)%QType)
                  CASE (PosX)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%r(1,p%OutParam(I)%NodeID)  ! x position
                  CASE (PosY)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%r(2,p%OutParam(I)%NodeID) ! y position
                  CASE (PosZ)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%r(3,p%OutParam(I)%NodeID) ! z position
                  CASE (VelX)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%rd(1,p%OutParam(I)%NodeID) ! x velocity
                  CASE (VelY)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%rd(2,p%OutParam(I)%NodeID) ! y velocity
                  CASE (VelZ)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%rd(3,p%OutParam(I)%NodeID) ! z velocity                     
                  CASE (AccX)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%a6(1) ! x acceleration <<< should this become distributed for each node?
                  CASE (AccY)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%a6(2) ! y acceleration
                  CASE (AccZ)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%a6(3) ! z acceleration
                  CASE (FX)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%F6net(1)  ! total force in x - added Nov 24
                  CASE (FY)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%F6net(2)  ! total force in y
                  CASE (FZ)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%F6net(3)  ! total force in z
                  CASE (Roll)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%roll*180.0/pi  ! rod roll
                  CASE (Pitch)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%pitch*180.0/pi ! rod pitch
                  CASE (Sub)
                     y%WriteOutput(I) = m%RodList(p%OutParam(I)%ObjID)%h0 / m%RodList(p%OutParam(I)%ObjID)%UnstrLen ! rod submergence
                  CASE DEFAULT
                     y%WriteOutput(I) = 0.0_ReKi
                     ErrStat = ErrID_Warn
                     ErrMsg = ' Unsupported output quantity '//TRIM(p%OutParam(I)%Name)//' requested from Rod '//TRIM(Num2Lstr(p%OutParam(I)%ObjID))//'.'
               END SELECT

            ELSE IF (p%OutParam(I)%OType == 4) THEN  ! if dealing with a Body output
               SELECT CASE (p%OutParam(I)%QType)
                  CASE (PosX)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%r6(1)  ! x position
                  CASE (PosY)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%r6(2) ! y position
                  CASE (PosZ)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%r6(3) ! z position
                  CASE (VelX)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%v6(1) ! x velocity
                  CASE (VelY)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%v6(2) ! y velocity
                  CASE (VelZ)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%v6(3) ! z velocity
                  CASE (FX)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%F6net(1)  ! total force in x - added Nov 24
                  CASE (FY)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%F6net(2)  ! total force in y
                  CASE (FZ)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%F6net(3)  ! total force in z
                  CASE (Roll)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%r6(4)*180.0/pi  ! roll
                  CASE (Pitch)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%r6(5)*180.0/pi  ! pitch
                  CASE (Yaw)
                     y%WriteOutput(I) = m%BodyList(p%OutParam(I)%ObjID)%r6(6)*180.0/pi  ! yaw
                  CASE DEFAULT
                     y%WriteOutput(I) = 0.0_ReKi
                     ErrStat = ErrID_Warn
                     ErrMsg = ' Unsupported output quantity '//TRIM(p%OutParam(I)%Name)//' requested from Body '//TRIM(Num2Lstr(p%OutParam(I)%ObjID))//'.'
               END SELECT


            ELSE  ! it must be an invalid output, so write zero
               y%WriteOutput(I) = 0.0_ReKi

            END IF

         END DO ! I, loop through OutParam

      END IF

         ! check if this is a repeated time step, in which case exit instead of writing a duplicate line to the output files
         if (Time <= m%LastOutTime) then
            return
         else
            m%LastOutTime = Time
         end if

         ! if using a certain output time step, check whether we should output, and exit the subroutine if not
         if (p%dtOut > 0)  then
            !if (Time < (floor((Time-p%dtCoupling)/p%dtOut) + 1.0)*p%dtOut)  then
            if ( abs(MOD( Time - 0.5*p%dtOut, p%dtOut) - 0.5*p%dtOut) >= 0.5*p%dtCoupling)  then
                return
            end if
         end if
         ! What the above does is say if ((dtOut==0) || (t >= (floor((t-dtC)/dtOut) + 1.0)*dtOut)), continue to writing files

      if ( p%NumOuts > 0_IntKi ) then  
      
         ! Write the output parameters to the file
         Frmt = '(F10.4,'//TRIM(Int2LStr(p%NumOuts))//'(A1,ES15.7E2))'   ! should evenutally use user specified format?
         
         WRITE(p%MDUnOut,Frmt)  Time, ( p%Delim, y%WriteOutput(I), I=1,p%NumOuts )
      END IF



      !------------------------------------------------------------------------
      ! now do the outputs for each line!  
      
      DO I=1,p%NLines
        
        IF (m%LineList(I)%OutFlagList(1) == 1) THEN    ! only proceed if the line is flagged to output a file
           
           ! calculate number of output entries to write for this line
           !LineNumOuts = 3*(m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(2:5)) + m%LineList(I)%N*SUM(m%LineList(I)%OutFlagList(6:9))
           
           LineNumOuts = 3*(m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(2:6)) &
                         + (m%LineList(I)%N + 1)*SUM(m%LineList(I)%OutFlagList(7:9)) &
                               + m%LineList(I)%N*SUM(m%LineList(I)%OutFlagList(10:18))
           
           if (m%LineList(I)%OutFlagList(2) == 1) THEN   ! if node positions are included, make them using a float format for higher precision
            Frmt = '(F10.4,'//TRIM(Int2LStr(3*(m%LineList(I)%N + 1)))//'(A1,ES15.7E2),'//TRIM(Int2LStr(LineNumOuts - 3*(m%LineList(I)%N - 1)))//'(A1,ES15.7E2))'  
           else
            Frmt = '(F10.4,'//TRIM(Int2LStr(LineNumOuts))//'(A1,ES15.7E2))'   ! should evenutally use user specified format?
           end if
           
           L = 1 ! start of index of line output file at first entry   12345.7890
           
           ! Time
      !     m%LineList(I)%LineWrOutput(L) = Time
      !     L = L+1
           
           ! Node positions
           IF (m%LineList(I)%OutFlagList(2) == 1) THEN
              DO J = 0,m%LineList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%LineList(I)%LineWrOutput(L) = m%LineList(I)%r(K,J)
                  L = L+1
                END DO
              END DO
           END IF         
           
           ! Node velocities
           IF (m%LineList(I)%OutFlagList(3) == 1) THEN
              DO J = 0,m%LineList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%LineList(I)%LineWrOutput(L) = m%LineList(I)%rd(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           
           ! Node wave velocities (not implemented yet)
           IF (m%LineList(I)%OutFlagList(4) == 1) THEN
              DO J = 0,m%LineList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%LineList(I)%LineWrOutput(L) = m%LineList(I)%U(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           
           ! Node total hydrodynamic forces (except added mass - just drag for now)
           IF (m%LineList(I)%OutFlagList(5) == 1) THEN
              DO J = 0,m%LineList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%LineList(I)%LineWrOutput(L) = m%LineList(I)%Dp(K,J) + m%LineList(I)%Dq(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           
           ! Node seabed contact force
           IF (m%LineList(I)%OutFlagList(6) == 1) THEN
              DO J = 0,m%LineList(I)%N  
                DO K = 1,3
                  m%LineList(I)%LineWrOutput(L) = m%LineList(I)%B(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           
           ! Node weights
           IF (m%LineList(I)%OutFlagList(7) == 1) THEN
              DO J = 0,m%LineList(I)%N
                  m%LineList(I)%LineWrOutput(L) = m%LineList(I)%W(3,J)
                  L = L+1
              END DO
           END IF
           
        !   ! Node curvatures
        !   IF (m%LineList(I)%OutFlagList(8) == 1) THEN
        !      DO J = 0,m%LineList(I)%N
        !          m%LineList(I)%LineWrOutput(L) = m%LineList(I)%W(3,J)
        !          L = L+1
        !      END DO
        !   END IF
           
           
           ! Segment tension force (excludes damping term, just EA)
           IF (m%LineList(I)%OutFlagList(10) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                m%LineList(I)%LineWrOutput(L) = TwoNorm(m%LineList(I)%T(:,J) )
                L = L+1
              END DO
           END IF
           
           ! Segment internal damping force
           IF (m%LineList(I)%OutFlagList(11) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                 IF (( m%LineList(I)%Td(3,J)*m%LineList(I)%T(3,J) ) > 0)  THEN  ! if statement for handling sign (positive = tension)
                    m%LineList(I)%LineWrOutput(L) = TwoNorm(m%LineList(I)%Td(:,J) )
                 ELSE
                    m%LineList(I)%LineWrOutput(L) = -TwoNorm(m%LineList(I)%Td(:,J) )
                 END IF
                 L = L+1
              END DO
           END IF
           
           ! Segment strain
           IF (m%LineList(I)%OutFlagList(12) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                m%LineList(I)%LineWrOutput(L) = m%LineList(I)%lstr(J)/m%LineList(I)%l(J) - 1.0 
                L = L+1
              END DO
           END IF
           
           ! Segment strain rate
           IF (m%LineList(I)%OutFlagList(13) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                m%LineList(I)%LineWrOutput(L) = m%LineList(I)%lstrd(J)/m%LineList(I)%l(J)
                L = L+1
              END DO
           END IF
           
           ! Segment length
           IF (m%LineList(I)%OutFlagList(14) == 1) THEN
              DO J = 1,m%LineList(I)%N  
                m%LineList(I)%LineWrOutput(L) = m%LineList(I)%lstr(J)
                L = L+1
              END DO
           END IF
                    
                    
           
           WRITE(m%LineList(I)%LineUnOut,Frmt) Time, ( p%Delim, m%LineList(I)%LineWrOutput(J), J=1,(LineNumOuts) )
           !WRITE(m%LineList(I)%LineUnOut,Frmt)  Time, ( p%Delim, m%LineList(I)%LineWrOutput(J), J=1,(3+3*m%LineList(I)%N) )

         END IF  ! if line output file flag is on
           
      END DO ! I
      
      
      
      !------------------------------------------------------------------------
      ! now do the outputs for each Rod!  
      
      DO I=1,p%NRods
        
        IF (m%RodList(I)%OutFlagList(1) == 1) THEN    ! only proceed if the line is flagged to output a file
           
           ! calculate number of output entries to write for this Rod
           RodNumOuts = 3*(m%RodList(I)%N + 1)*SUM(m%RodList(I)%OutFlagList(2:9)) &
                         + (m%RodList(I)%N + 1)*SUM(m%RodList(I)%OutFlagList(10:11)) &
                               + m%RodList(I)%N*SUM(m%RodList(I)%OutFlagList(12:18))
           
           
           Frmt = '(F10.4,'//TRIM(Int2LStr(RodNumOuts))//'(A1,ES15.7E2))'   ! should evenutally use user specified format?

           L = 1 ! start of index of line output file at first entry
           
           ! Time
      !     m%RodList(I)%RodWrOutput(L) = Time
      !     L = L+1
           
           ! Node positions
           IF (m%RodList(I)%OutFlagList(2) == 1) THEN
              DO J = 0,m%RodList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%RodList(I)%RodWrOutput(L) = m%RodList(I)%r(K,J)
                  L = L+1
                END DO
              END DO
           END IF         
           
           ! Node velocities
           IF (m%RodList(I)%OutFlagList(3) == 1) THEN
              DO J = 0,m%RodList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%RodList(I)%RodWrOutput(L) = m%RodList(I)%rd(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           
           ! Node wave velocities (not implemented yet)
           IF (m%RodList(I)%OutFlagList(4) == 1) THEN
              DO J = 0,m%RodList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%RodList(I)%RodWrOutput(L) = m%RodList(I)%U(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           ! Node buoyancy forces
           IF (m%RodList(I)%OutFlagList(5) == 1) THEN
              DO J = 0,m%RodList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%RodList(I)%RodWrOutput(L) = m%RodList(I)%Bo(K,J)
                  L = L+1
                END DO
              END DO
           END IF  
           
           ! Node drag forces
           IF (m%RodList(I)%OutFlagList(6) == 1) THEN
              DO J = 0,m%RodList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%RodList(I)%RodWrOutput(L) = m%RodList(I)%Dp(K,J) + m%RodList(I)%Dq(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           ! Node inertia forces
           IF (m%RodList(I)%OutFlagList(7) == 1) THEN
              DO J = 0,m%RodList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%RodList(I)%RodWrOutput(L) = m%RodList(I)%Ap(K,J) + m%RodList(I)%Aq(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           ! Node dynamic pressure forces
           IF (m%RodList(I)%OutFlagList(8) == 1) THEN
              DO J = 0,m%RodList(I)%N  ! note index starts at zero because these are nodes
                DO K = 1,3
                  m%RodList(I)%RodWrOutput(L) = m%RodList(I)%Pd(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           ! Node seabed contact force
           IF (m%RodList(I)%OutFlagList(9) == 1) THEN
              DO J = 0,m%RodList(I)%N  
                DO K = 1,3
                  m%RodList(I)%RodWrOutput(L) = m%RodList(I)%B(K,J)
                  L = L+1
                END DO
              END DO
           END IF
           
           
           ! Node weights
           IF (m%RodList(I)%OutFlagList(10) == 1) THEN
              DO J = 0,m%RodList(I)%N
                  m%RodList(I)%RodWrOutput(L) = m%RodList(I)%W(3,J)
                  L = L+1
              END DO
           END IF
           
        !   ! Node curvatures
        !   IF (m%RodList(I)%OutFlagList(8) == 1) THEN
        !      DO J = 0,m%RodList(I)%N
        !          m%RodList(I)%RodWrOutput(L) = m%RodList(I)%W(3,J)
        !          L = L+1
        !      END DO
        !   END IF
           
           
           ! Segment tension force (excludes damping term, just EA)
           ! N/A
           
           ! Segment internal damping force 
           ! N/A
           
           ! Segment strain
           ! N/A
           
           ! Segment strain rate
           ! N/A
                    
           
           WRITE(m%RodList(I)%RodUnOut,Frmt) Time, ( p%Delim, m%RodList(I)%RodWrOutput(J), J=1,(RodNumOuts) )
           
         END IF  ! if line output file flag is on
           
      END DO ! I

   END SUBROUTINE MDIO_WriteOutputs
   !----------------------------------------------------------------------------------------============


   ! get tension at any node including fairlead or anchor (accounting for weight in these latter cases)
   !--------------------------------------------------------------
   FUNCTION Line_GetNodeTen(Line, i, p) result(NodeTen)

      TYPE(MD_Line),          INTENT(IN   )  :: Line           ! label for the current line, for convenience
      INTEGER(IntKi),         INTENT(IN   )  :: i              ! node index to get tension at
      TYPE(MD_ParameterType), INTENT(IN   )  :: p              ! Parameters
      REAL(DbKi)                             :: NodeTen        ! returned calculation of tension at node   
      
      INTEGER(IntKi)                   :: J      
      REAL(DbKi)                       :: Tmag_squared   
   
      if (i==0) then
         NodeTen = sqrt( Line%Fnet(1,i)**2 + Line%Fnet(2,i)**2 + Line%Fnet(3,i)**2 )  ! if an end node, use Fnet which already includes weight
      else if (i==Line%N) then                          
         NodeTen = sqrt( Line%Fnet(1,i)**2 + Line%Fnet(2,i)**2 + Line%Fnet(3,i)**2 )
      else 
         Tmag_squared = 0.0_DbKi 
         DO J=1,3
            Tmag_squared = Tmag_squared + 0.25*(Line%T(J,i) + Line%Td(J,i) + Line%T(J,i+1) + Line%Td(J,i+1))**2   ! take average of tension in adjacent segments 
         END DO
         NodeTen = sqrt(Tmag_squared) 
      end if

   END FUNCTION Line_GetNodeTen
   !--------------------------------------------------------------


END MODULE MoorDyn_IO
