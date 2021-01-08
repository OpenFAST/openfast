!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  Matthew Hall
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
MODULE MoorDyn

   USE MoorDyn_Types
   USE MoorDyn_IO
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: MD_ProgDesc = ProgDesc( 'MoorDyn-F', 'v2.a1', '5 Jan. 2020' )


   PUBLIC :: MD_Init
   PUBLIC :: MD_UpdateStates
   PUBLIC :: MD_CalcOutput
   PUBLIC :: MD_CalcContStateDeriv
   PUBLIC :: MD_End

CONTAINS

   !=========================================   MD_Init   ===================================
   SUBROUTINE MD_Init(InitInp, u, p, x, xd, z, other, y, m, DTcoupling, InitOut, ErrStat, ErrMsg)

      IMPLICIT NONE

      TYPE(MD_InitInputType),       INTENT(INOUT)  :: InitInp     ! INTENT(INOUT) : Input data for initialization routine
      TYPE(MD_InputType),           INTENT(  OUT)  :: u           ! INTENT( OUT) : An initial guess for the input; input mesh must be defined
      TYPE(MD_ParameterType),       INTENT(  OUT)  :: p           ! INTENT( OUT) : Parameters
      TYPE(MD_ContinuousStateType), INTENT(  OUT)  :: x           ! INTENT( OUT) : Initial continuous states
      TYPE(MD_DiscreteStateType),   INTENT(  OUT)  :: xd          ! INTENT( OUT) : Initial discrete states
      TYPE(MD_ConstraintStateType), INTENT(  OUT)  :: z           ! INTENT( OUT) : Initial guess of the constraint states
      TYPE(MD_OtherStateType),      INTENT(  OUT)  :: other       ! INTENT( OUT) : Initial other states
      TYPE(MD_OutputType),          INTENT(  OUT)  :: y           ! INTENT( OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)
      TYPE(MD_MiscVarType),         INTENT(  OUT)  :: m           ! INTENT( OUT) : Initial misc/optimization variables
      REAL(DbKi),                   INTENT(INOUT)  :: DTcoupling  ! Coupling interval in seconds: the rate that Output is the actual coupling interval
      TYPE(MD_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      REAL(DbKi)                                   :: t              ! instantaneous time, to be used during IC generation
      INTEGER(IntKi)                               :: l              ! index
      INTEGER(IntKi)                               :: I              ! index
      INTEGER(IntKi)                               :: J              ! index
      INTEGER(IntKi)                               :: K              ! index
      INTEGER(IntKi)                               :: Itemp          ! index
      INTEGER(IntKi)                               :: Converged      ! flag indicating whether the dynamic relaxation has converged
      INTEGER(IntKi)                               :: N              ! convenience integer for readability: number of segments in the line
      REAL(ReKi)                                   :: Pos(3)         ! array for setting absolute fairlead positions in mesh
      REAL(ReKi)                                   :: OrMat(3,3)     ! rotation matrix for setting fairlead positions correctly if there is initial platform rotation
      REAL(DbKi), ALLOCATABLE                      :: FairTensIC(:,:)! array of size nCpldCons, 3 to store three latest fairlead tensions of each line
      CHARACTER(20)                                :: TempString     ! temporary string for incidental use
      INTEGER(IntKi)                               :: ErrStat2       ! Error status of the operation
      CHARACTER(ErrMsgLen)                         :: ErrMsg2        ! Error message if ErrStat2 /= ErrID_None
      
      REAL(DbKi)                                      :: dtM         ! actual mooring dynamics time step
      INTEGER(IntKi)                                  :: NdtM        ! number of time steps to integrate through with RK2
      
      TYPE(MD_InputType)    :: u_array(1)    ! a size-one array for u to make call to TimeStep happy
      REAL(DbKi)            :: t_array(1)    ! a size-one array saying time is 0 to make call to TimeStep happy  
      TYPE(MD_InputType)                                  :: u_interp   ! interpolated instantaneous input values to be calculated for each mooring time step



      CHARACTER(MaxWrScrLen)                       :: Message
      
      ! Local variables for reading file input (Previously in MDIO_ReadInput)
      INTEGER(IntKi)               :: UnIn                 ! Unit number for the input file
      INTEGER(IntKi)               :: UnEc                 ! The local unit number for this module's echo file
    INTEGER(IntKi)   :: UnOut    ! for outputing wave kinematics data
    CHARACTER(200)   :: Frmt     ! a string to hold a format statement

      CHARACTER(1024)              :: EchoFile             ! Name of MoorDyn echo file
      CHARACTER(1024)              :: Line                 ! String to temporarially hold value of read line
      CHARACTER(20)                :: LineOutString        ! String to temporarially hold characters specifying line output options
      CHARACTER(20)                :: OptString            ! String to temporarially hold name of option variable
      CHARACTER(20)                :: OptValue             ! String to temporarially hold value of options variable input
      INTEGER(IntKi)               :: nOpts = 0            ! number of options lines in input file
      CHARACTER(40)                :: TempString1          !
      CHARACTER(40)                :: TempString2          !
      CHARACTER(40)                :: TempString3          !
      CHARACTER(40)                :: TempString4          !
      CHARACTER(1024)              :: FileName             !
      
      
      CHARACTER(25)                 :: let1                ! strings used for splitting and parsing identifiers
      CHARACTER(25)                 :: num1
      CHARACTER(25)                 :: let2
      CHARACTER(25)                 :: num2
      CHARACTER(25)                 :: let3
      
      REAL(DbKi)                    :: tempArray(6)
      REAL(ReKi)                    :: rRef(6)
      REAL(DbKi)                    :: rRefDub(3)
      
      ! for reading output channels
      CHARACTER(ChanLen),ALLOCATABLE :: OutList(:)          ! array of output channel request (moved here from InitInput)
      INTEGER                       :: MaxAryLen = 1000    ! Maximum length of the array being read
      INTEGER                       :: NumWords            ! Number of words contained on a line
      INTEGER                       :: Nx
      CHARACTER(*), PARAMETER       :: RoutineName = 'MD_Init'

      

      ErrStat = ErrID_None
      ErrMsg  = ""
      m%zeros6 = 0.0_DbKi

      ! Initialize the NWTC Subroutine Library
      CALL NWTC_Init( )

      ! Display the module information
      CALL DispNVD( MD_ProgDesc )
      InitOut%Ver = MD_ProgDesc


      !---------------------------------------------------------------------------------------------
      !                   Get all the inputs taken care of
      !---------------------------------------------------------------------------------------------


      ! set environmental parameters from input data and error check
      ! (should remove these values as options from MoorDyn input file for consistency?)

      p%g        = InitInp%g
      p%WtrDpth  = InitInp%WtrDepth
      p%rhoW     = InitInp%rhoW

      p%RootName = TRIM(InitInp%RootName)//'.MD'  ! all files written from this module will have this root name

      !---------------------------------------------------------------------------------------------
      !            read input file and create cross-referenced mooring system objects
      !---------------------------------------------------------------------------------------------
      
      
      UnEc = -1

      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      !-------------------------------------------------------------------------------------------------
      ! Open the file
      !-------------------------------------------------------------------------------------------------
      FileName = TRIM(InitInp%FileName)

      CALL GetNewUnit( UnIn )
      CALL OpenFInpFile( UnIn, FileName, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL CleanUp()
            RETURN
         END IF


      !CALL WrScr( '  MD_Init: Opening MoorDyn input file:  '//FileName )


		
		! ----------------- go through file contents a first time, counting each entry -----------------------
		
      i = 0
      read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i       !read a line
      
      do while ( ErrStat2 == 0 ) 
      
		
         if (INDEX(Line, "---") > 0) then ! look for a header line

            if ( ( INDEX(Line, "LINE DICTIONARY") > 0) .or. ( INDEX(Line, "LINE TYPES") > 0) ) then ! if line dictionary header
print *, "line dictionary"
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! find how many elements of this type there are
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nLineTypes = p%nLineTypes + 1
                  read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i," linetype on prev line"
               END DO

            else if ( (INDEX(Line, "ROD DICTIONARY") > 0) .or. ( INDEX(Line, "ROD TYPES") > 0) ) then ! if rod dictionary header
print *, "rod dictionary"
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! find how many elements of this type there are
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nRodTypes = p%nRodTypes + 1
                  read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i," rod type on prev line"
               END DO	

            else if ((INDEX(Line, "BODY LIST") > 0 ) .or. (INDEX(Line, "BODY PROPERTIES") > 0 )) then
print *, "body list"
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! find how many elements of this type there are
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nBodies = p%nBodies + 1
                  read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i," body on prev line"
               END DO		

            else if ((INDEX(Line, "ROD LIST") > 0) .or. (INDEX(Line, "ROD PROPERTIES") > 0)) then ! if rod properties header
print *, "rod list"
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! find how many elements of this type there are
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nRods = p%nRods + 1
                  read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i," rod on prev line"
               END DO	

            else if ( (INDEX(Line, "CONNECTION PROPERTIES") > 0) .or. (INDEX(Line, "NODE PROPERTIES") > 0) ) then ! if node properties header
print *, "connections"
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! find how many elements of this type there are
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nConnects = p%nConnects + 1
                  read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i," con on prev line"
               END DO

            else if (INDEX(Line, "LINE PROPERTIES") > 0) then ! if line properties header
print *, "lines"
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! find how many elements of this type there are
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nLines = p%nLines + 1
                  read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i," line on prev line"
               END DO

            else if (INDEX(Line, "FAILURE") > 0) then ! if failure conditions header

               print *, "   Reading failure conditions: ";
               
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! find how many elements of this type there are
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nFails = p%nFails + 1
                  read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               END DO
               
               
            else if (INDEX(Line, "OPTIONS") > 0) then ! if options header
print *, "options"
               ! don't skip any lines (no column headers for the options section)

               ! find how many options have been specified
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  nOpts = nOpts + 1
                  read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i," option on prev line"
               END DO
               

            else if (INDEX(Line, "OUTPUT") > 0) then ! if output header
print *, "output"
               ! we don't need to count this section...

               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i


            else  ! otherwise ignore this line that isn't a recognized header line and read the next line
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i," --- unrecognized header"
            end if
			
         else ! otherwise ignore this line, which doesn't have the "---" or header line and read the next line
            read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i, " .."
         end if
     
      end do
		
      p%nConnectsExtra = p%nConnects + 2*p%nLines    ! set maximum number of connections, accounting for possible detachment of each line end and a connection for that
		
      print *, "Identified ", p%nLineTypes  , "LineTypes in input file."
      print *, "Identified ", p%nRodTypes   , "RodTypes in input file."
      print *, "Identified ", p%nBodies     , "Bodies in input file."
      print *, "Identified ", p%nRods       , "Rods in input file."
      print *, "Identified ", p%nConnects   , "Connections in input file."
      print *, "Identified ", p%nLines      , "Lines in input file."



      ! ----------------------------- misc checks to be sorted -----------------------------


    ! make sure nLineTypes isn't zero
    IF ( p%nLineTypes < 1 ) THEN
       CALL SetErrStat( ErrID_Fatal, 'nLineTypes parameter must be greater than zero.', ErrStat, ErrMsg, RoutineName )
       CALL CleanUp()
       RETURN
    END IF
    
    ! make sure NLines is at least one
    IF ( p%NLines < 1 ) THEN
       CALL SetErrStat( ErrID_Fatal, 'NLines parameter must be at least 1.', ErrStat, ErrMsg, RoutineName )
       CALL CleanUp()
       RETURN
    END IF


    
    

      ! ----------------------------- allocate necessary arrays ----------------------------


      ! Allocate object arrays

      ALLOCATE(m%LineTypeList(p%nLineTypes), STAT = ErrStat2 ); if(AllocateFailed("LineTypeList")) return
      ALLOCATE(m%RodTypeList( p%nRodTypes ), STAT = ErrStat2 ); if(AllocateFailed("LineTypeList")) return
               
      ALLOCATE(m%BodyList(    p%nBodies   ), STAT = ErrStat2 ); if(AllocateFailed("BodyList"    )) return
      ALLOCATE(m%RodList(     p%nRods     ), STAT = ErrStat2 ); if(AllocateFailed("RodList"     )) return
      ALLOCATE(m%ConnectList( p%nConnects ), STAT = ErrStat2 ); if(AllocateFailed("ConnectList" )) return
      ALLOCATE(m%LineList(    p%nLines    ), STAT = ErrStat2 ); if(AllocateFailed("LineList"    )) return
      
      ALLOCATE(m%FailList(    p%nFails    ), STAT = ErrStat2 ); if(AllocateFailed("FailList"    )) return
    
      
      ! Allocate associated index arrays (note: some are allocated larger than will be used, for simplicity)
      ALLOCATE(m%BodyStateIs1(p%nBodies ), m%BodyStateIsN(p%nBodies ), STAT=ErrStat2); if(AllocateFailed("BodyStateIs1/N")) return
      ALLOCATE(m%RodStateIs1(p%nRods    ), m%RodStateIsN(p%nRods    ), STAT=ErrStat2); if(AllocateFailed("RodStateIs1/N" )) return
      ALLOCATE(m%ConStateIs1(p%nConnects), m%ConStateIsN(p%nConnects), STAT=ErrStat2); if(AllocateFailed("ConStateIs1/N" )) return
      ALLOCATE(m%LineStateIs1(p%nLines)  , m%LineStateIsN(p%nLines)  , STAT=ErrStat2); if(AllocateFailed("LineStateIs1/N")) return

      ALLOCATE(m%FreeBodyIs(   p%nBodies ), m%CpldBodyIs(p%nBodies ), STAT=ErrStat2); if(AllocateFailed("BodyIs")) return
      ALLOCATE(m%FreeRodIs(    p%nRods   ), m%CpldRodIs( p%nRods   ), STAT=ErrStat2); if(AllocateFailed("RodIs")) return
      ALLOCATE(m%FreeConIs(   p%nConnects), m%CpldConIs(p%nConnects),STAT=ErrStat2); if(AllocateFailed("ConnectIs")) return


      ! ---------------------- now go through again and process file contents --------------------

      REWIND(UnIn)      ! rewind to start of input file
		
      ! note: no longer worrying about "Echo" option
      
      Nx = 0  ! set state counter to zero
      i  = 0  ! set line number counter to zero
      
      read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i       !read a line
      
      do while ( ErrStat2 == 0 ) 
      
         if (INDEX(Line, "---") > 0) then ! look for a header line

            !-------------------------------------------------------------------------------------------
            if ( ( INDEX(Line, "LINE DICTIONARY") > 0) .or. ( INDEX(Line, "LINE TYPES") > 0) ) then ! if line dictionary header
               
               print *, "Reading line types"
               
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
                ! process each line
                DO l = 1,p%nLineTypes
                   
                   !read into a line
                   READ(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i

                   ! parse out entries: Name  Diam MassDenInAir EA cIntDamp >>EI(new)<<  Can  Cat Cdn  Cdt 
                   READ(Line,*,IOSTAT=ErrStat2) m%LineTypeList(l)%name, m%LineTypeList(l)%d,  &
                      m%LineTypeList(l)%w, tempString1, tempString2, tempString3, &
                      m%LineTypeList(l)%Can, m%LineTypeList(l)%Cat, m%LineTypeList(l)%Cdn, m%LineTypeList(l)%Cdt
                   
                    IF ( ErrStat2 /= ErrID_None ) THEN
                      CALL SetErrStat( ErrID_Fatal, 'Failed to process line type inputs of entry '//trim(Num2LStr(l))//'. Check formatting and correct number of columns.', ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                      RETURN
                   END IF
                   
                   ! process stiffness, damping, and bending coefficients (which might use lookup tables)
                   CALL getCoefficientOrCurve(tempString1, m%LineTypeList(l)%EA,        &
                                                           m%LineTypeList(l)%nEApoints, &
                                                           m%LineTypeList(l)%stiffXs,   &
                                                           m%LineTypeList(l)%stiffYs,  ErrStat2, ErrMsg2)
                   CALL getCoefficientOrCurve(tempString2, m%LineTypeList(l)%BA,        &
                                                           m%LineTypeList(l)%nBpoints,  &
                                                           m%LineTypeList(l)%dampXs,    &
                                                           m%LineTypeList(l)%dampYs,   ErrStat2, ErrMsg2)
                   CALL getCoefficientOrCurve(tempString3, m%LineTypeList(l)%EI,        &
                                                           m%LineTypeList(l)%nEIpoints, &
                                                           m%LineTypeList(l)%bstiffXs,  &
                                                           m%LineTypeList(l)%bstiffYs, ErrStat2, ErrMsg2)
								
                   ! specify IdNum of line type for error checking
                   m%LineTypeList(l)%IdNum = l  


                   IF ( ErrStat2 /= ErrID_None ) THEN
                      CALL SetErrStat( ErrID_Fatal, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                      RETURN
                   END IF

                END DO


            !-------------------------------------------------------------------------------------------
            else if ( (INDEX(Line, "ROD DICTIONARY") > 0) .or. ( INDEX(Line, "ROD TYPES") > 0) ) then ! if rod dictionary header
               
               print *, "Reading rod types"
               
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
                ! process each line
                DO l = 1,p%nRodTypes
                   
                   !read into a line
                   READ(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i

                   ! parse out entries: Name  Diam MassDenInAir Can  Cat Cdn  Cdt 
                   IF (ErrStat2 == 0) THEN
                      READ(Line,*,IOSTAT=ErrStat2) m%RodTypeList(l)%name, m%RodTypeList(l)%d, m%RodTypeList(l)%w, &
                         m%RodTypeList(l)%Can, m%RodTypeList(l)%Cat, m%RodTypeList(l)%Cdn, m%RodTypeList(l)%Cdt, &
                         m%RodTypeList(l)%CaEnd, m%RodTypeList(l)%CdEnd   
                   END IF

                   ! specify IdNum of rod type for error checking
                   m%RodTypeList(l)%IdNum = l  


                   IF ( ErrStat2 /= ErrID_None ) THEN
                      CALL SetErrStat( ErrID_Fatal, 'Failed to process rod type properties for rod '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                      RETURN
                   END IF

                END DO


            !-------------------------------------------------------------------------------------------
            else if ((INDEX(Line, "BODY LIST") > 0 ) .or. (INDEX(Line, "BODY PROPERTIES") > 0 )) then

               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! process each body
               DO l = 1,p%nBodies
                  
                  !read into a line
                  READ(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i

                  ! parse out entries: Node Type X Y Z M V FX FY FZ CdA Ca 
                  IF (ErrStat2 == 0) THEN
                     READ(Line,*,IOSTAT=ErrStat2) tempString1, &
                        tempArray(1), tempArray(2), tempArray(3), tempArray(4), tempArray(5), tempArray(6), &
                        m%BodyList(l)%rCG(1), m%BodyList(l)%rCG(2), m%BodyList(l)%rCG(3), &
                        m%BodyList(l)%bodyM, m%BodyList(l)%bodyV, m%BodyList(l)%bodyI(1), m%BodyList(l)%bodyI(2), m%BodyList(l)%bodyI(3), &
                        m%BodyList(l)%bodyCdA(1), m%BodyList(l)%bodyCa(1)
                  END IF

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL WrScr('   Unable to parse Body '//trim(Num2LStr(l))//' on row '//trim(Num2LStr(i))//' in input file.')  ! Specific screen output because errors likely
                     CALL WrScr('   Ensure row has all 17 columns.')  
                        CALL SetErrStat( ErrID_Fatal, 'Failed to read bodies.' , ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  


                  !----------- process body type (not considering input fixed bodies for now, only the GroundBody) -----------------
							
                  call DecomposeString(tempString1, let1, num1, let2, num2, let3)
                  	
                  READ(num1, *) m%BodyList(l)%IdNum   ! convert to int, representing parent body index
                                          
							if ((let2 == "COUPLED") .or. (let2 == "VESSEL")) then    ! if a coupled body
                     
                     m%BodyList(l)%typeNum = -1
                     p%nCpldBodies=p%nCpldBodies+1  ! add this rod to coupled list                          
                     m%CpldBodyIs(p%nCpldBodies) = l

                     ! body initial position due to coupling will be adjusted later
                     
                  else if ((let2 == "FREE") .or. (LEN_TRIM(let2)== 0)) then    ! if a free body
                     m%BodyList(l)%typeNum = 0
                     
                     p%nFreeBodies=p%nFreeBodies+1             ! add this pinned rod to the free list because it is half free
                     
                     m%BodyStateIs1(p%nFreeBodies) = Nx+1
                     m%BodyStateIsN(p%nFreeBodies) = Nx+12                     				 
                     Nx = Nx + 12                           ! add 12 state variables for free Body
                     
                     m%FreeBodyIs(p%nFreeBodies) = l
                     
                     m%BodyList(l)%r6 = tempArray     ! set initial body position and orientation
                     
                  else 
                     CALL SetErrStat( ErrID_Fatal,  "Unidentified Body type string for Body "//trim(Num2LStr(l))//": "//trim(tempString2), ErrStat, ErrMsg, RoutineName )   
                     return
                  end if
                  

                  ! check for sequential IdNums
                  IF ( m%BodyList(l)%IdNum .NE. l ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Body numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  ! set initial velocity to zero
                  m%BodyList(l)%v6 = 0.0_DbKi
                  
                  !also set number of attached rods and points to zero initially
                  m%BodyList(l)%nAttachedC = 0
                  m%BodyList(l)%nAttachedR = 0

                  ! if there was a body setup function, it would get called here, but I don't think it's needed.                  

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read data for body '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  print *, "Set up body ", l, " of type ",  m%BodyList(l)%typeNum

               END DO
               
               
            !-------------------------------------------------------------------------------------------
            else if ((INDEX(Line, "ROD LIST") > 0) .or. (INDEX(Line, "ROD PROPERTIES") > 0)) then ! if rod properties header

               print *, "Reading rods"
               
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! process each rod
               DO l = 1,p%nRods
                  
                  !read into a line
                  READ(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i

                  ! parse out entries: RodID  Type/BodyID  RodType  Xa   Ya   Za   Xb   Yb   Zb  NumSegs  Flags/Outputs
                  IF (ErrStat2 == 0) THEN
                     READ(Line,*,IOSTAT=ErrStat2) m%RodList(l)%IdNum, tempString1, tempString2, &
                         tempArray(1), tempArray(2), tempArray(3), tempArray(4), tempArray(5), tempArray(6), &
                         m%RodList(l)%N, LineOutString
                  END IF


                  !----------- process rod type -----------------
							
                  call DecomposeString(tempString1, let1, num1, let2, num2, let3)
                  	
							if ((let1 == "ANCHOR") .or. (let1 == "FIXED") .or. (let1 == "FIX")) then
								 m%RodList(l)%typeNum = 2
                     CALL Body_AddRod(m%GroundBody, l, tempArray)   ! add rod l to Ground body
                           
								
							else if ((let1 == "PINNED") .or. (let1 == "PIN")) then
						      m%RodList(l)%typeNum = 1
                     CALL Body_AddRod(m%GroundBody, l, tempArray)   ! add rod l to Ground body
                     
                     p%nFreeRods=p%nFreeRods+1  ! add this pinned rod to the free list because it is half free
                     
                     m%RodStateIs1(p%nFreeRods) = Nx+1
                     m%RodStateIsN(p%nFreeRods) = Nx+6                     				 
                     Nx = Nx + 6                                               ! add 6 state variables for each pinned rod	
                     
                     m%FreeRodIs(p%nFreeRods) = l
                     
                  else if (let1 == "BODY") then ! attached to a body (either rididly or pinned)
                  
                     if (len_trim(num1) > 0) then
                     
                        READ(num1,*) J   ! convert to int, representing parent body index
                        
                        if ((J <= p%nBodies) .and. (J > 0)) then 
                        
                           CALL Body_AddRod(m%BodyList(J), l, tempArray)   ! add rod l to the body
                           
                           if ((INDEX(let2, "PINNED") > 0) .or. (INDEX(let2, "PIN") > 0)) then
                              m%RodList(l)%typeNum = 1
                              
                              p%nFreeRods=p%nFreeRods+1  ! add this pinned rod to the free list because it is half free
                              
                              m%RodStateIs1(p%nFreeRods) = Nx+1
                              m%RodStateIsN(p%nFreeRods) = Nx+6                     				 
                              Nx = Nx + 6                                               ! add 6 state variables for each pinned rod	
                              
                              m%FreeRodIs(p%nFreeRods) = l
                     
                           else
                              m%RodList(l)%typeNum = 2
                           end if
                        
                        else
                           CALL SetErrStat( ErrID_Severe,  "Body ID out of bounds for Rod "//trim(Num2LStr(l))//".", ErrStat, ErrMsg, RoutineName )  
                           return
                        end if
                     
                  else
                     CALL SetErrStat( ErrID_Severe,  "No number provided for Rod "//trim(Num2LStr(l))//" Body attachment.", ErrStat, ErrMsg, RoutineName )   
                         return
                  end if
                  
                  else if ((let1 == "VESSEL") .or. (let1 == "VES")) then    ! if a rigid fairlead, add to list and add 
                     m%RodList(l)%typeNum = -2                     
                     m%CpldRodIs(p%nCpldRods) = l;  p%nCpldRods=p%nCpldRods+1  ! add this rod to coupled list                 				 	
                 
                  else if ((let1 == "VESSELPINNED") .or. (let1 == "VESPIN")) then  ! if a pinned fairlead, add to list and add 
                     m%RodList(l)%typeNum = -1
                     
                     p%nCpldRods=p%nCpldRods+1  ! add
                     p%nFreeRods=p%nFreeRods+1  ! add this pinned rod to the free list because it is half free
                     
                     m%RodStateIs1(p%nFreeRods) = Nx+1
                     m%RodStateIsN(p%nFreeRods) = Nx+6                     				 
                     Nx = Nx + 6                                               ! add 6 state variables for each pinned rod	
                     
                     m%CpldRodIs(p%nCpldRods) = l
                     m%FreeRodIs(p%nFreeRods) = l
                     
                  else if ((let1 == "CONNECT") .or. (let1 == "CON") .or. (let1 == "FREE")) then
                     m%RodList(l)%typeNum = 0
                     
                     p%nFreeRods=p%nFreeRods+1  ! add this pinned rod to the free list because it is half free
                     
                     m%RodStateIs1(p%nFreeRods) = Nx+1
                     m%RodStateIsN(p%nFreeRods) = Nx+12                     				 
                     Nx = Nx + 12                                              ! add 12 state variables for free Rod
                     
                     m%FreeRodIs(p%nFreeRods) = l
                     
                  else 
                  
                     CALL SetErrStat( ErrID_Severe,  "Unidentified Type/BodyID for Rod "//trim(Num2LStr(l))//": "//trim(tempString1), ErrStat, ErrMsg, RoutineName )   
                     return
                  end if


                  ! find Rod properties index
                  DO J = 1,p%nRodTypes
                     IF (trim(tempString2) == trim(m%RodTypeList(J)%name)) THEN
                       m%RodList(l)%PropsIdNum = J
                       EXIT
                       IF (J == p%nRodTypes) THEN   ! call an error if there is no match
                           CALL SetErrStat( ErrID_Severe, 'Unable to find matching rod type name for Rod '//trim(Num2LStr(l))//": "//trim(tempString2), ErrStat, ErrMsg, RoutineName )
                       END IF
                     END IF
                  END DO
                  
                
                  ! process output flag characters (LineOutString) and set line output flag array (OutFlagList)
                  m%RodList(l)%OutFlagList = 0  ! first set array all to zero
                  ! per node, 3 component
                  IF ( scan( LineOutString, 'p') > 0 )  m%RodList(l)%OutFlagList(2 ) = 1   ! node position
                  IF ( scan( LineOutString, 'v') > 0 )  m%RodList(l)%OutFlagList(3 ) = 1   ! node velocity
                  IF ( scan( LineOutString, 'U') > 0 )  m%RodList(l)%OutFlagList(4 ) = 1   ! water velocity
                  IF ( scan( LineOutString, 'B') > 0 )  m%RodList(l)%OutFlagList(5 ) = 1   ! node buoyancy force
                  IF ( scan( LineOutString, 'D') > 0 )  m%RodList(l)%OutFlagList(6 ) = 1   ! drag force
                  IF ( scan( LineOutString, 'I') > 0 )  m%RodList(l)%OutFlagList(7 ) = 1   ! inertia force
                  IF ( scan( LineOutString, 'P') > 0 )  m%RodList(l)%OutFlagList(8 ) = 1   ! dynamic pressure force
                  IF ( scan( LineOutString, 'b') > 0 )  m%RodList(l)%OutFlagList(9 ) = 1   ! seabed contact forces
                  ! per node, 1 component
                  IF ( scan( LineOutString, 'W') > 0 )  m%RodList(l)%OutFlagList(10) = 1   ! node weight/buoyancy (positive up)
                  IF ( scan( LineOutString, 'K') > 0 )  m%RodList(l)%OutFlagList(11) = 1   ! curvature at node
                  ! per element, 1 component >>> these don't apply to a rod!! <<<
                  IF ( scan( LineOutString, 't') > 0 )  m%RodList(l)%OutFlagList(12) = 1  ! segment tension force (just EA)
                  IF ( scan( LineOutString, 'c') > 0 )  m%RodList(l)%OutFlagList(13) = 1  ! segment internal damping force
                  IF ( scan( LineOutString, 's') > 0 )  m%RodList(l)%OutFlagList(14) = 1  ! Segment strain
                  IF ( scan( LineOutString, 'd') > 0 )  m%RodList(l)%OutFlagList(15) = 1  ! Segment strain rate

                  IF (SUM(m%RodList(l)%OutFlagList) > 0)   m%RodList(l)%OutFlagList(1) = 1  ! this first entry signals whether to create any output file at all
                  ! the above letter-index combinations define which OutFlagList entry corresponds to which output type


                  ! specify IdNum of line for error checking
                  m%RodList(l)%IdNum = l  

                  ! check for sequential IdNums
                  IF ( m%RodList(l)%IdNum .NE. l ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Line numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  ! set up rod
                  CALL Rod_Setup( m%RodList(l), m%RodTypeList(m%RodList(l)%PropsIdNum), tempArray, p%rhoW, ErrStat2, ErrMsg2)
                     CALL CheckError( ErrStat2, ErrMsg2 )
                     IF (ErrStat >= AbortErrLev) RETURN
                     
                  ! note: Rod was already added to its respective parent body if type > 0
							
                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read rod data for Rod '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF

               END DO   ! l = 1,p%nRods



            !-------------------------------------------------------------------------------------------
            else if ( (INDEX(Line, "CONNECTION PROPERTIES") > 0) .or. (INDEX(Line, "NODE PROPERTIES") > 0) ) then ! if node properties header
               
               print *, "Reading connections"
               
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! process each point
               DO l = 1,p%nConnects
                  
                  !read into a line
                  READ(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i

                  ! parse out entries: Node Type X Y Z M V FX FY FZ CdA Ca 
                  IF (ErrStat2 == 0) THEN
                     READ(Line,*,IOSTAT=ErrStat2) m%ConnectList(l)%IdNum, tempString1, tempArray(1), &
                        tempArray(2), tempArray(3), m%ConnectList(l)%conM, &
                        m%ConnectList(l)%conV, m%ConnectList(l)%conFX, m%ConnectList(l)%conFY, &
                        m%ConnectList(l)%conFZ, m%ConnectList(l)%conCdA, m%ConnectList(l)%conCa
                  END IF
                  

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL WrScr('   Unable to parse Connection '//trim(Num2LStr(I))//' row in input file.')  ! Specific screen output because errors likely
                     CALL WrScr('   Ensure row has all 12 columns, including CdA and Ca.')           ! to be caused by non-updated input file formats.
                        CALL SetErrStat( ErrID_Fatal, 'Failed to read connects.' , ErrStat, ErrMsg, RoutineName ) ! would be nice to specify which line <<<<<<<<<
                     CALL CleanUp()
                     RETURN
                  END IF
                  


                  !----------- process connection type -----------------
							
                  call DecomposeString(tempString1, let1, num1, let2, num2, let3)
                  	
							if ((let1 == "ANCHOR") .or. (let1 == "FIXED") .or. (let1 == "FIX")) then
								 m%ConnectList(l)%typeNum = 1
                     
                     m%ConnectList(l)%r = tempArray(1:3)   ! set initial node position
                     
                     CALL Body_AddConnect(m%GroundBody, l, tempArray(1:3))   ! add connection l to Ground body
                     
								
						   else if (let1 == "BODY") then ! attached to a body
                     if (len_trim(num1) > 0) then                     
                        READ(num1, *) J   ! convert to int, representing parent body index
                        
                        if ((J <= p%nBodies) .and. (J > 0)) then
                           m%ConnectList(l)%typeNum = 1    

                           CALL Body_AddConnect(m%BodyList(J), l, tempArray(1:3))   ! add connection l to Ground body
                           
                        else
                           CALL SetErrStat( ErrID_Severe,  "Body ID out of bounds for Connection "//trim(Num2LStr(l))//".", ErrStat, ErrMsg, RoutineName )  
                           return
                        end if                     
                     else
                        CALL SetErrStat( ErrID_Severe,  "No number provided for Connection "//trim(Num2LStr(l))//" Body attachment.", ErrStat, ErrMsg, RoutineName )   
                            return
                     end if
                  
                  else if ((let1 == "VESSEL") .or. (let1 == "VES")) then    ! if a fairlead, add to list and add 
                     m%ConnectList(l)%typeNum = -1
                     p%nCpldCons=p%nCpldCons+1  ! add this rod to coupled list                          
                     m%CpldConIs(p%nCpldCons) = l

                     ! this is temporary for backwards compatibility >>>>> will need to update for more versatile coupling >>>>
                     CALL SmllRotTrans('PtfmInit', InitInp%PtfmInit(4),InitInp%PtfmInit(5),InitInp%PtfmInit(6), OrMat, '', ErrStat2, ErrMsg2)

                     ! set initial node position, including adjustments due to initial platform rotations and translations  <<< could convert to array math
                     m%ConnectList(l)%r(1) = InitInp%PtfmInit(1) + OrMat(1,1)*tempArray(1) + OrMat(2,1)*tempArray(2) + OrMat(3,1)*tempArray(3)
                     m%ConnectList(l)%r(2) = InitInp%PtfmInit(2) + OrMat(1,2)*tempArray(1) + OrMat(2,2)*tempArray(2) + OrMat(3,2)*tempArray(3)
                     m%ConnectList(l)%r(3) = InitInp%PtfmInit(3) + OrMat(1,3)*tempArray(1) + OrMat(2,3)*tempArray(2) + OrMat(3,3)*tempArray(3)
                 
                  else if ((let1 == "CONNECT") .or. (let1 == "CON") .or. (let1 == "FREE")) then
                     m%ConnectList(l)%typeNum = 0
                     
                     p%nFreeCons=p%nFreeCons+1             ! add this pinned rod to the free list because it is half free
                     
                     m%ConStateIs1(p%nFreeCons) = Nx+1
                     m%ConStateIsN(p%nFreeCons) = Nx+6                     				 
                     Nx = Nx + 6                           ! add 12 state variables for free Connection
                     
                     m%FreeConIs(p%nFreeCons) = l
                     
                     m%ConnectList(l)%r = tempArray(1:3)   ! set initial node position
                     
                     

                  else 
                     CALL SetErrStat( ErrID_Severe,  "Unidentified Type/BodyID for Connection "//trim(Num2LStr(l))//": "//trim(tempString2), ErrStat, ErrMsg, RoutineName )   
                     return
                  end if
                  
                  ! set initial velocity to zero
                  m%ConnectList(l)%rd(1) = 0.0_DbKi
                  m%ConnectList(l)%rd(2) = 0.0_DbKi
                  m%ConnectList(l)%rd(3) = 0.0_DbKi
                  
                  ! possibly redundant <<< should revisit                  
                  m%ConnectList(l)%conX = tempArray(1)
                  m%ConnectList(l)%conY = tempArray(2)
                  m%ConnectList(l)%conZ = tempArray(3)

         
                  !also set number of attached lines to zero initially
                  m%ConnectList(l)%nAttached = 0


                  ! check for sequential IdNums
                  IF ( m%ConnectList(l)%IdNum .NE. l ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Connection numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read rod data for Connection '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  print *, "Set up connection ", l, " of type ",  m%ConnectList(l)%typeNum

               END DO   ! l = 1,p%nRods

            !-------------------------------------------------------------------------------------------
            else if (INDEX(Line, "LINE PROPERTIES") > 0) then ! if line properties header

               print *, "Reading lines"
               
               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! process each line
               DO l = 1,p%nLines
                  
                  !read into a line
                  READ(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i

                   ! parse out entries: LineID  LineType  UnstrLen  NumSegs  NodeA  NodeB  Flags/Outputs
                  IF (ErrStat2 == 0) THEN
                     READ(Line,*,IOSTAT=ErrStat2) m%LineList(l)%IdNum, tempString1, m%LineList(l)%UnstrLen, &
                        m%LineList(l)%N, tempString2, tempString3, LineOutString
                  END IF

                  ! identify index of line type
                  DO J = 1,p%nLineTypes
                     IF (trim(tempString1) == trim(m%LineTypeList(J)%name)) THEN
                       m%LineList(l)%PropsIdNum = J
                       EXIT
                       IF (J == p%nLineTypes) THEN   ! call an error if there is no match
                           CALL SetErrStat( ErrID_Severe, 'Unable to find matching line type name for Line '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                       END IF
                     END IF
                  END DO
                  
                  ! account for states of line
                  m%LineStateIs1(l) = Nx + 1
                  m%LineStateIsN(l) = Nx + 6*m%LineList(l)%N - 6
                  Nx = Nx + 6*(m%LineList(l)%N - 1)       
                
                  ! Process attachment identfiers and attach line ends 
                  
                  ! First for the anchor (or end A)...
                  
                  call DecomposeString(tempString2, let1, num1, let2, num2, let3)
                  
                  if (len_trim(num1)<1) then
                     CALL SetErrStat( ErrID_Severe,  "Error: no number provided for line "//trim(Num2LStr(l))//" end A attachment.", ErrStat, ErrMsg, RoutineName )  
                     return
                  end if 

                  READ(num1, *) J   ! convert to int
                     
                  ! if id starts with an "R" or "Rod"
                  if ((let1 == "R") .or. (let1 == "ROD")) then   
                  
                     if ((J <= p%nRods) .and. (J > 0)) then                  
                        if (let2 == "A") then
                           CALL Rod_AddLine(m%RodList(J), l, 0, 0)   ! add line l (end A, denoted by 0) to rod J (end A, denoted by 0)
                        else if (let2 == "B") then 
                           CALL Rod_AddLine(m%RodList(J), l, 0, 1)   ! add line l (end A, denoted by 0) to rod J (end B, denoted by 1)
                        else
                           CALL SetErrStat( ErrID_Severe,  "Error: rod end (A or B) must be specified for line "//trim(Num2LStr(l))//" end A attachment. Instead seeing "//let2, ErrStat, ErrMsg, RoutineName )  
                            return
                        end if
                     else
                        CALL SetErrStat( ErrID_Severe,  "Error: rod connection ID out of bounds for line "//trim(Num2LStr(l))//" end A attachment.", ErrStat, ErrMsg, RoutineName )  
                        return							
                     end if
                  
                     ! if J starts with a "C" or "Con" or goes straight ot the number then it's attached to a Connection
                  else if ((len_trim(let1)==0) .or. (let1 == "C") .or. (let1 == "CON")) then 

                     if ((J <= p%nConnects) .and. (J > 0)) then                  
                        CALL Connect_AddLine(m%ConnectList(J), l, 0)   ! add line l (end A, denoted by 0) to connection J
                     else
                        CALL SetErrStat( ErrID_Severe,  "Error: connection out of bounds for line "//trim(Num2LStr(l))//" end A attachment.", ErrStat, ErrMsg, RoutineName )  
                        return							
                     end if
                        
                  end if


                  ! Then again for the fairlead (or end B)...

                  call DecomposeString(tempString3, let1, num1, let2, num2, let3)

                  if (len_trim(num1)<1) then
                     CALL SetErrStat( ErrID_Severe,  "Error: no number provided for line "//trim(Num2LStr(l))//" end B attachment.", ErrStat, ErrMsg, RoutineName )  
                     return
                  end if 

                  READ(num1, *) J   ! convert to int
                     
                  ! if id starts with an "R" or "Rod"
                  if ((let1 == "R") .or. (let1 == "ROD")) then   

                     if ((J <= p%nRods) .and. (J > 0)) then                  
                        if (let2 == "A") then
                           CALL Rod_AddLine(m%RodList(J), l, 1, 0)   ! add line l (end B, denoted by 1) to rod J (end A, denoted by 0)
                        else if (let2 == "B") then 
                           CALL Rod_AddLine(m%RodList(J), l, 1, 1)   ! add line l (end B, denoted by 1) to rod J (end B, denoted by 1)
                        else
                           CALL SetErrStat( ErrID_Severe,  "Error: rod end (A or B) must be specified for line "//trim(Num2LStr(l))//" end B attachment. Instead seeing "//let2, ErrStat, ErrMsg, RoutineName )  
                            return
                        end if
                     else
                        CALL SetErrStat( ErrID_Severe,  "Error: rod connection ID out of bounds for line "//trim(Num2LStr(l))//" end B attachment.", ErrStat, ErrMsg, RoutineName )  
                        return							
                     end if

                  ! if J starts with a "C" or "Con" or goes straight ot the number then it's attached to a Connection
                  else if ((len_trim(let1)==0) .or. (let1 == "C") .or. (let1 == "CON")) then 

                     if ((J <= p%nConnects) .and. (J > 0)) then                  
                        CALL Connect_AddLine(m%ConnectList(J), l, 1)   ! add line l (end B, denoted by 1) to connection J
                     else
                        CALL SetErrStat( ErrID_Severe,  "Error: connection out of bounds for line "//trim(Num2LStr(l))//" end B attachment.", ErrStat, ErrMsg, RoutineName )  
                        return							
                     end if
                        
                  end if
                
                
                  ! process output flag characters (LineOutString) and set line output flag array (OutFlagList)
                  m%LineList(l)%OutFlagList = 0  ! first set array all to zero
                  ! per node 3 component
                  IF ( scan( LineOutString, 'p') > 0 )  m%LineList(l)%OutFlagList(2) = 1 
                  IF ( scan( LineOutString, 'v') > 0 )  m%LineList(l)%OutFlagList(3) = 1
                  IF ( scan( LineOutString, 'U') > 0 )  m%LineList(l)%OutFlagList(4) = 1
                  IF ( scan( LineOutString, 'D') > 0 )  m%LineList(l)%OutFlagList(5) = 1
                  IF ( scan( LineOutString, 'b') > 0 )  m%LineList(l)%OutFlagList(6) = 1   ! seabed contact forces
                  ! per node 1 component
                  IF ( scan( LineOutString, 'W') > 0 )  m%LineList(l)%OutFlagList(7) = 1  ! node weight/buoyancy (positive up)
                  IF ( scan( LineOutString, 'K') > 0 )  m%LineList(l)%OutFlagList(8) = 1  ! curvature at node
                  ! per element 1 component
                  IF ( scan( LineOutString, 't') > 0 )  m%LineList(l)%OutFlagList(10) = 1  ! segment tension force (just EA)
                  IF ( scan( LineOutString, 'c') > 0 )  m%LineList(l)%OutFlagList(11) = 1  ! segment internal damping force
                  IF ( scan( LineOutString, 's') > 0 )  m%LineList(l)%OutFlagList(12) = 1  ! Segment strain
                  IF ( scan( LineOutString, 'd') > 0 )  m%LineList(l)%OutFlagList(13) = 1  ! Segment strain rate

                  IF (SUM(m%LineList(l)%OutFlagList) > 0)   m%LineList(l)%OutFlagList(1) = 1  ! this first entry signals whether to create any output file at all
                  ! the above letter-index combinations define which OutFlagList entry corresponds to which output type


                  ! specify IdNum of line for error checking
                  m%LineList(l)%IdNum = l  


                  ! check for sequential IdNums
                  IF ( m%LineList(l)%IdNum .NE. l ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Line numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  
                  ! setup line
                  CALL SetupLine( m%LineList(l), m%LineTypeList(m%LineList(l)%PropsIdNum), p%rhoW ,  ErrStat2, ErrMsg2)
                     CALL CheckError( ErrStat2, ErrMsg2 )
                     IF (ErrStat >= AbortErrLev) RETURN
                     

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read line data for Line '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF

               END DO   ! l = 1,p%nLines



            !-------------------------------------------------------------------------------------------
            else if (INDEX(Line, "FAILURE") > 0) then ! if failure conditions header

               print *, "   Reading failure conditions: (not implemented yet) ";
               
               ! TODO: add stuff <<<<<<<<

               ! skip following two lines (label line and unit line)
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
               
               ! process each line
               DO l = 1,p%nFails
                  
                  !read into a line
                  READ(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i
                  
               END DO

            !-------------------------------------------------------------------------------------------
            else if (INDEX(Line, "OPTIONS") > 0) then ! if options header

               print *, "Reading options"
               
               ! (don't skip any lines)
               
               ! process each line
               DO l = 1,nOpts
                  
                  !read into a line
                  READ(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i      

                  ! parse out entries:  value, option keyword
                  READ(Line,*,IOSTAT=ErrStat2) OptValue, OptString  ! look at first two entries, ignore remaining words in line, which should be comments
                  

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read options.', ErrStat, ErrMsg, RoutineName ) ! would be nice to specify which line had the error
                     CALL CleanUp()
                     RETURN
                  END IF
                              
                  CALL Conv2UC(OptString)

                  ! check all possible options types and see if OptString is one of them, in which case set the variable.
                  if ( OptString == 'DTM') THEN
                     read (OptValue,*) p%dtM0   ! InitInp%DTmooring
                  else if ( OptString == 'G') then
                     read (OptValue,*) p%g
                  else if ( OptString == 'RHOW') then
                     read (OptValue,*) p%rhoW
                  else if ( OptString == 'WTRDPTH') then
                     read (OptValue,*) p%WtrDpth
                  else if ( OptString == 'KBOT')  then
                     read (OptValue,*) p%kBot
                  else if ( OptString == 'CBOT')  then
                     read (OptValue,*) p%cBot
                  else if ( OptString == 'DTIC')  then
                     read (OptValue,*) InitInp%dtIC
                  else if ( OptString == 'TMAXIC')  then
                     read (OptValue,*) InitInp%TMaxIC
                  else if ( OptString == 'CDSCALEIC')  then
                     read (OptValue,*) InitInp%CdScaleIC
                  else if ( OptString == 'THRESHIC')  then
                     read (OptValue,*) InitInp%threshIC
                  else if ( OptString == 'WATERKIN')  then
                     read (OptValue,*) p%WaterKin
                  else
                     CALL SetErrStat( ErrID_Warn, 'unable to interpret input '//trim(OptString), ErrStat, ErrMsg, RoutineName ) 
                  end if

               END DO


            !-------------------------------------------------------------------------------------------
            else if (INDEX(Line, "OUTPUT") > 0) then ! if output header

               print *, "Reading outputs"
               
               ! (don't skip any lines)
                        
               ! allocate InitInp%Outliest (to a really big number for now...)
               CALL AllocAry( OutList, MaxAryLen, "MoorDyn Input File's Outlist", ErrStat2, ErrMsg2 ); if(Failed()) return

               ! OutList - List of user-requested output channels (-):
               !CALL ReadOutputList ( UnIn, FileName, InitInp%OutList, p%NumOuts, 'OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc ); if(Failed()) return

               ! customm implementation to avoid need for "END" keyword line
 
               ! Initialize some values
               p%NumOuts = 0    ! start counter at zero
               OutList = ''


               ! Read in all of the lines containing output parameters and store them in OutList(:)
               DO
                  ! read a line
                  READ(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i

                  CALL Conv2UC(Line)   ! convert to uppercase for easy string matching

                  if ((INDEX(Line, "---") > 0) .or. (INDEX(Line, "END") > 0)) EXIT ! stop if we hit a header line or the keyword "END"
                     
                  NumWords = CountWords( Line )    ! The number of words in Line.

                  p%NumOuts = p%NumOuts + NumWords  ! The total number of output channels read in so far.


                  IF ( p%NumOuts > MaxAryLen )  THEN  ! Check to see if the maximum # allowable in the array has been reached.

                     ErrStat = ErrID_Fatal
                     ErrMsg = 'Error while reading output channels: The maximum number of output channels allowed is '//TRIM( Int2LStr(MaxAryLen) )//'.'
                     EXIT

                  ELSE
                     CALL GetWords ( Line, OutList((p%NumOuts - NumWords + 1):p%NumOuts), NumWords )

                  END IF

               END DO

               ! process the OutList array and set up the index arrays for the requested output quantities
               CALL MDIO_ProcessOutList(OutList, p, m, y, InitOut, ErrStat2, ErrMsg2 )
                  CALL CheckError( ErrStat2, ErrMsg2 )
                  IF (ErrStat >= AbortErrLev) RETURN


            !-------------------------------------------------------------------------------------------
            else  ! otherwise ignore this line that isn't a recognized header line and read the next line
               read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i," --- unrecognized header"
            end if
			
            !-------------------------------------------------------------------------------------------
         
         else ! otherwise ignore this line, which doesn't have the "---" or header line and read the next line
            read(UnIn,'(A)',IOSTAT=ErrStat2) Line; i=i+1; print*,i, " .."
         end if
     
      end do


      ! this is the end of reading the input file, so close it now
      CALL CleanUp()
      
      
      
      
      
      
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN



      !-------------------------------------------------------------------------------------------------
      !          Connect mooring system together and make necessary allocations
      !-------------------------------------------------------------------------------------------------

      CALL WrNR( '   Created mooring system.  ' )

!     p%NAnchs = 0   ! this is the number of "fixed" type Connections. <<<<<<<<<<<<<<


!      CALL WrScr(trim(Num2LStr(p%nCpldCons))//' fairleads, '//trim(Num2LStr(p%NAnchs))//' anchors, '//trim(Num2LStr(p%nFreeCons))//' connects.')




 !     ! now go back through and record the fairlead Id numbers (this >>>WAS<<< all the "connecting" that's required) <<<<
 !     J = 1  ! counter for fairlead number
 !     K = 1  ! counter for connect number
 !     DO I = 1,p%NConnects
 !        IF (m%ConnectList(I)%typeNum == 1) THEN
 !          m%CpldConIs(J) = I             ! if a vessel connection, add ID to list
 !          J = J + 1
 !        ELSE IF (m%ConnectList(I)%typeNum == 2) THEN
 !          m%FreeConIs(K) = I             ! if a connect connection, add ID to list
 !          K = K + 1
 !        END IF
 !     END DO

   print *, "nLineTypes     = ",p%nLineTypes    
   print *, "nRodTypes      = ",p%nRodTypes     
   print *, "nConnects      = ",p%nConnects     
   print *, "nConnectsExtra = ",p%nConnectsExtra
   print *, "nBodies        = ",p%nBodies       
   print *, "nRods          = ",p%nRods         
   print *, "nLines         = ",p%nLines        
   print *, "nFails         = ",p%nFails        
   print *, "nFreeBodies    = ",p%nFreeBodies   
   print *, "nFreeRods      = ",p%nFreeRods     
   print *, "nFreeCons      = ",p%nFreeCons     
   print *, "nCpldBodies    = ",p%nCpldBodies   
   print *, "nCpldRods      = ",p%nCpldRods     
   print *, "nCpldCons      = ",p%nCpldCons     
   print *, "NConns         = ",p%NConns        
   print *, "NAnchs         = ",p%NAnchs              
      
   print *, "FreeConIs are ", m%FreeConIs
   print *, "CpldConIs are ", m%CpldConIs


      !------------------------------------------------------------------------------------
      !                          fill in state vector index record holders
      !------------------------------------------------------------------------------------

      ! allocate state vector index record holders...

      

 !    ! allocate list of starting and ending state vector indices for each free connection
 !    ALLOCATE ( m%ConStateIs1(p%nFreeCons), m%ConStateIsN(p%nFreeCons), STAT = ErrStat )
 !    IF ( ErrStat /= ErrID_None ) THEN
 !      CALL CheckError(ErrID_Fatal, ' Error allocating ConStateIs array.')
 !      RETURN
 !    END IF
 !    
 !    ! allocate list of starting and ending state vector indices for each line  - does this belong elsewhere?
 !    ALLOCATE ( m%LineStateIs1(p%nLines), m%LineStateIsN(p%nLines), STAT = ErrStat )
 !    IF ( ErrStat /= ErrID_None ) THEN
 !      CALL CheckError(ErrID_Fatal, ' Error allocating LineStateIs arrays.')
 !      RETURN
 !    END IF
 !          
 !    
 !    ! fill in values for state vector index record holders...
 !    
 !    J=0 ! start off index counter at zero
 !    
 !    ! Free Bodies...
 !    ! Free Rods...
 !    
 !    ! Free Connections...
 !    DO l = 1, p%nFreeCons
 !       J = J + 1            ! assign start index
 !       m%ConStateIs1(l) = J
 !       
 !       J = J + 5            ! assign end index (5 entries further, since nodes have 2*3 states)
 !       m%ConStateIsN(l) = J
 !    END DO
 !    
 !    ! Lines
 !    DO l = 1, p%nLines
 !       J = J + 1            ! assign start index
 !       m%LineStateIs1(l) = J
 !       
 !       J = J + 6*(m%LineList(l)%N - 1) - 1  ! !add 6 state variables for each internal node
 !       m%LineStateIsN(l) = J
 !    END DO
 !
 !
 !    ! record number of states
 !    m%Nx = J
 
 
      !------------------------------------------------------------------------------------
      !                               prepare state vector etc.
      !------------------------------------------------------------------------------------

      ! the number of states is Nx
      m%Nx = Nx
      
      print *, "allocating state vectors to size ", Nx

      ! allocate state vector and temporary state vectors based on size just calculated
      ALLOCATE ( x%states(m%Nx), m%xTemp%states(m%Nx), m%xdTemp%states(m%Nx), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
        ErrMsg  = ' Error allocating state vectors.'
        !CALL CleanUp()
        RETURN
      END IF




      ! ================================ initialize system ================================



      ! call ground body to update all the fixed things...
      CALL Body_SetDependentKin(m%GroundBody, 0.0_DbKi, m)

    ! m%GroundBody%OrMat = EulerConstruct( m%GroundBody%r6(4:6) ) ! make sure it's OrMat is set up  <<< need to check this approach
      
   !   ! first set/update the kinematics of all the fixed things (>>>> eventually do this by using a ground body <<<<)
   !   ! only doing connections so far
   !   DO J = 1,p%nConnects 
   !      if (m%ConnectList(J)%typeNum == 1) then
   !         ! set the attached line endpoint positions:
   !         CALL Connect_SetKinematics(m%ConnectList(J), m%ConnectList(J)%r, (/0.0_DbKi,0.0_DbKi,0.0_DbKi/), 0.0_DbKi, m%LineList) 
   !      end if
   !   END DO 
      
      
      
      
      ! Initialize coupled objects based on passed kinematics
      ! (set up initial condition of each coupled object based on values specified by glue code) 
      ! Also create i/o meshes 
      
      ! <<<<<<<< look here when changing to shared mooring compatibility

      ! create input mesh for all coupled objects (as mesh node points)
      CALL MeshCreate(BlankMesh=u%CoupledKinematics, IOS= COMPONENT_INPUT, &
                    Nnodes = p%nCpldBodies + p%nCpldRods + p%nCpldCons, &
                    TranslationDisp=.TRUE., TranslationVel=.TRUE., &
                    Orientation=.TRUE., RotationVel=.TRUE., &
                    TranslationAcc=.TRUE., RotationAcc= .TRUE., &
                     ErrStat=ErrStat2, ErrMess=ErrMsg2)

      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
            
      ! note: in MoorDyn-F v2, the points in the mesh correspond in order to all the coupled bodies, then rods, then connections
      ! NOTE: InitInp%PtfmInit should be replaced by specific values for each coupled body/rod/connect at some point <<<<<
      
      J = 0 ! this is the counter through the mesh points
      
      DO l = 1,p%nCpldBodies
      
         J = J + 1
      
         rRef = m%BodyList(m%CpldBodyIs(l))%r6              ! for now set reference position as per input file <<< 
         
         CALL MeshPositionNode(u%CoupledKinematics, J, rRef(1:3), ErrStat2, ErrMsg2) ! defaults to identity orientation matrix
         u%CoupledKinematics%TranslationDisp(:,J) = 0.0_ReKi   ! no displacement from reference position

         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN

      !   ! Apply offsets due to initial platform rotations and translations (fixed Jun 19, 2015)
      !   CALL SmllRotTrans('body rotation matrix', InitInp%PtfmInit(4),InitInp%PtfmInit(5),InitInp%PtfmInit(6), OrMat, '', ErrStat2, ErrMsg2)
      !   u%CoupledKinematics%TranslationDisp(1, i) = InitInp%PtfmInit(1) + OrMat(1,1)*rRef(1) + OrMat(2,1)*rRef(2) + OrMat(3,1)*rRef(3) - rRef(1)
      !   u%CoupledKinematics%TranslationDisp(2, i) = InitInp%PtfmInit(2) + OrMat(1,2)*rRef(1) + OrMat(2,2)*rRef(2) + OrMat(3,2)*rRef(3) - rRef(2)
      !   u%CoupledKinematics%TranslationDisp(3, i) = InitInp%PtfmInit(3) + OrMat(1,3)*rRef(1) + OrMat(2,3)*rRef(2) + OrMat(3,3)*rRef(3) - rRef(3)
         
         CALL MeshConstructElement(u%CoupledKinematics, ELEMENT_POINT, ErrStat2, ErrMsg2, J)      ! set node as point element
         CALL Body_InitializeUnfree( m%BodyList(m%CpldBodyIs(l)), m )

      END DO 
      
      DO l = 1,p%nCpldRods   ! keeping this one simple for now, positioning at whatever is specified in input file <<<<< should change to glue code!

         J = J + 1
         
         rRef = m%RodList(m%CpldRodIs(l))%r6    ! for now set reference position as per input file <<< 
         
         CALL MeshPositionNode(u%CoupledKinematics, J, rRef, ErrStat2, ErrMsg2)  ! defaults to identity orientation matrix
         u%CoupledKinematics%TranslationDisp(:,J) = 0.0_ReKi   ! no displacement from reference position
         CALL MeshConstructElement(u%CoupledKinematics, ELEMENT_POINT, ErrStat2, ErrMsg2, J)
         CALL Rod_SetKinematics(m%RodList(m%CpldRodIs(l)), tempArray, m%zeros6, m%zeros6, 0.0_DbKi, m%LineList) 
      END DO 

      DO l = 1,p%nCpldCons   ! keeping this one simple for now, positioning at whatever is specified by glue code <<<
      
         J = J + 1
         
         rRef(1:3) = m%ConnectList(m%CpldConIs(l))%r         ! for now set reference position as per input file <<<
         CALL MeshPositionNode(u%CoupledKinematics, J, rRef, ErrStat2, ErrMsg2)   ! defaults to identity orientation matrix
         u%CoupledKinematics%TranslationDisp(:,J) = 0.0_ReKi   ! no displacement from reference position
         CALL MeshConstructElement(u%CoupledKinematics, ELEMENT_POINT, ErrStat2, ErrMsg2, J)

         ! lastly, do this to set the attached line endpoint positions:
         rRefDub = rRef(1:3)
         CALL Connect_SetKinematics(m%ConnectList(m%CpldConIs(l)), rRefDub, m%zeros6(1:3), m%zeros6(1:3), 0.0_DbKi, m%LineList)
      END DO 

         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
         
      ! set velocities/accelerations of all mesh nodes to zero
      u%CoupledKinematics%TranslationVel = 0.0_ReKi
      u%CoupledKinematics%TranslationAcc = 0.0_ReKi
      u%CoupledKinematics%RotationVel    = 0.0_ReKi
      u%CoupledKinematics%RotationAcc    = 0.0_ReKi

      CALL MeshCommit ( u%CoupledKinematics, ErrStat, ErrMsg )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN

      ! copy the input fairlead kinematics mesh to make the output mesh for fairlead loads, PtFairleadLoad
      CALL MeshCopy ( SrcMesh  = u%CoupledKinematics, DestMesh = y%CoupledLoads, &
                      CtrlCode = MESH_SIBLING,             IOS      = COMPONENT_OUTPUT, &
                      Force=.TRUE.,  Moment=.TRUE.,  ErrStat  = ErrStat2, ErrMess=ErrMsg2 )

      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      
      ! ----------------------------- Arrays for active tensioning ---------------------------
      
      ! size active tensioning inputs arrays based on highest channel number read from input file for now <<<<<<<
      
      ! find the highest channel number
      J = 0
      DO I = 1, p%NLines
         IF ( m%LineList(I)%CtrlChan > J ) then
            J = m%LineList(I)%CtrlChan       
         END IF
      END DO   

      ! allocate the input arrays
      ALLOCATE ( u%DeltaL(J), u%DeltaLdot(J), STAT = ErrStat2 )
      
      
      
      
      ! ----------------------------- Arrays for wave kinematics -----------------------------
  !   :::::::::::::: BELOW WILL BE USED EVENTUALLY WHEN WAVE INFO IS AN INPUT ::::::::::::::::::
  !   ! The rAll array contains all nodes or reference points in the system 
  !   ! (x,y,z global coordinates for each) in the order of bodies, rods, points, internal line nodes.      
  !   
  !   ! count the number of nodes to use for passing wave kinematics
  !   J=0 
  !   ! Body reference point coordinates
  !   J = J + p%nBodies
  !   ! Rod node coordinates (including ends)
  !   DO l = 1, p%nRods
  !      J = J + (m%RodList(l)%N + 1)
  !   END DO
  !   ! Point reference point coordinates
  !   J = J + p%nConnects
  !   ! Line internal node coordinates
  !   DO l = 1, p%nLines
  !      J = J + (m%LineList(l)%N - 1)
  !   END DO
  !
  !   ! allocate all relevant arrays
  !   ! allocate state vector and temporary state vectors based on size just calculated
  !   ALLOCATE ( y%rAll(3,J), u%U(3,J), u%Ud(3,J), u%zeta(J), u%PDyn(J), STAT = ErrStat )
  !   IF ( ErrStat /= ErrID_None ) THEN
  !     ErrMsg  = ' Error allocating wave kinematics vectors.'
  !     RETURN
  !   END IF
  !
  !
  !   ! go through the nodes and fill in the data (this should maybe be turned into a global function)
  !   J=0 
  !   ! Body reference point coordinates
  !   DO I = 1, p%nBodies
  !      J = J + 1                     
  !      y%rAll(:,J) = m%BodyList(I)%r6(1:3)         
  !   END DO
  !   ! Rod node coordinates
  !   DO I = 1, p%nRods
  !      DO K = 0,m%RodList(I)%N  
  !         J = J + 1             
  !         y%rAll(:,J) = m%RodList(I)%r(:,K)
  !      END DO
  !   END DO
  !   ! Point reference point coordinates
  !   DO I = 1, p%nConnects
  !      J = J + 1
  !      y%rAll(:,J) = m%ConnectList(I)%r
  !   END DO      
  !   ! Line internal node coordinates
  !   DO I = 1, p%nLines
  !      DO K = 1,m%LineList(I)%N-1
  !         J = J + 1               
  !         y%rAll(:,J) = m%LineList(I)%r(:,K)
  !      END DO
  !   END DO  
      
      ! :::::::::::::::: the above will be used eventually. For now, let's store wave info grids within this module :::::::::::::::::
      ! allocate arrays
      I = SIZE(InitInp%WaveTime)
      ALLOCATE ( p%ux  (I,8,5,5), STAT = ErrStat2 )
      ALLOCATE ( p%uy  (I,8,5,5), STAT = ErrStat2 )
      ALLOCATE ( p%uz  (I,8,5,5), STAT = ErrStat2 )
      ALLOCATE ( p%ax  (I,8,5,5), STAT = ErrStat2 )
      ALLOCATE ( p%ay  (I,8,5,5), STAT = ErrStat2 )
      ALLOCATE ( p%az  (I,8,5,5), STAT = ErrStat2 )
      ALLOCATE ( p%PDyn(I,8,5,5), STAT = ErrStat2 )
      ALLOCATE ( p%zeta(I,5,5), STAT = ErrStat2 )    ! 2D grid over x and y only
      ALLOCATE ( p%px(5), STAT = ErrStat2 )
      ALLOCATE ( p%py(5), STAT = ErrStat2 )
      ALLOCATE ( p%pz(8), STAT = ErrStat2 )
      
      ! get grid and time info (currenltly this is hard-coded to match what's in HydroDyn_Input
      DO I=1,8
         p%pz(I) =  1.0 - 2.0**(8-I)       !  -127,  -63,  -31,  -15,   -7,   -3,   -1,    0
      END DO
      DO J = 1,5
         p%py(J) = -60.0 + 20.0*J
      END DO
      DO K = 1,5   
         p%px(K) = -60.0 + 20.0*K
      END DO
      p%dtWave = InitInp%WaveTime(2) - InitInp%WaveTime(1)
      
      ! fill in the grid data (the for loops match those in HydroDyn_Input)
      DO I=1,8
         DO J = 1,5
            DO K = 1,5   
               Itemp = (I-1)*25.0 + (J-1)*5.0 + K    ! index of actual node on 3D grid
               
               p%ux  (:,I,J,K) = InitInp%WaveVel( :,Itemp,1)  ! note: indices are t, z, y, x
               p%uy  (:,I,J,K) = InitInp%WaveVel( :,Itemp,2)
               p%uz  (:,I,J,K) = InitInp%WaveVel( :,Itemp,3)
               p%ax  (:,I,J,K) = InitInp%WaveAcc( :,Itemp,1)
               p%ay  (:,I,J,K) = InitInp%WaveAcc( :,Itemp,2)
               p%az  (:,I,J,K) = InitInp%WaveAcc( :,Itemp,3)
               p%PDyn(:,I,J,K) = InitInp%WavePDyn(:,Itemp)
            END DO
         END DO
      END DO
      
      ! fill in the grid data (the for loops match those in HydroDyn_Input)      
      DO J = 1,5
         DO K = 1,5   
            Itemp = (J-1)*5.0 + K    ! index of actual node on surface 2D grid   
            p%zeta(:,J,K) = InitInp%WaveElev(:,Itemp)
         END DO
      END DO
      
      
      ! write the date to an output file for testing purposes
      
      CALL GetNewUnit( UnOut)

      CALL OpenFOutFile ( UnOut, "waves.txt", ErrStat, ErrMsg )
      IF ( ErrStat > ErrID_None ) THEN
         ErrMsg = ' Error opening wave grid file: '//TRIM(ErrMsg)
         ErrStat = ErrID_Fatal
         RETURN
      END IF

      WRITE(UnOut, *, IOSTAT=ErrStat2)  TRIM( 'MoorDyn v2 wave/current kinematics grid file' )
      WRITE(UnOut, *, IOSTAT=ErrStat2)  TRIM( '---------------------------------------------' )
      WRITE(UnOut, *, IOSTAT=ErrStat2)  TRIM( 'The following 6 lines (4-9) specify the input type then the inputs for x, then, y, then z coordinates.' )
      
      WRITE(UnOut,*, IOSTAT=ErrStat2)  TRIM( '1  - X input type (0: not used; 1: list values in ascending order; 2: uniform specified by -xlim, xlim, num)' )
      Frmt = '('//TRIM(Int2LStr(5))//'(A1,e10.4))'      
      WRITE(UnOut,*, IOSTAT=ErrStat2)  ( " ", TRIM(Num2LStr(p%px(I))), I=1,5 )
      
      WRITE(UnOut,*, IOSTAT=ErrStat2)  TRIM( '1  - Y input type (0: not used; 1: list values in ascending order; 2: uniform specified by -xlim, xlim, num)' )
      Frmt = '('//TRIM(Int2LStr(5))//'(A1,e10.4))'      
      WRITE(UnOut,*, IOSTAT=ErrStat2)  ( " ", TRIM(Num2LStr(p%py(I))), I=1,5 )
      
      WRITE(UnOut,*, IOSTAT=ErrStat2)  TRIM( '1  - Z input type (0: not used; 1: list values in ascending order; 2: uniform specified by -xlim, xlim, num)' )
      Frmt = '('//TRIM(Int2LStr(8))//'(A1,e10.4))'      
      WRITE(UnOut,*, IOSTAT=ErrStat2)  ( " ", TRIM(Num2LStr(p%pz(I))), I=1,8 )
      
      CLOSE(UnOut, IOSTAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrMsg = 'Error closing wave grid file'
      END IF
      
      
      CALL GetNewUnit( UnOut)

      CALL OpenFOutFile ( UnOut, "wave data.txt", ErrStat, ErrMsg )
      IF ( ErrStat > ErrID_None ) THEN
         ErrMsg = ' Error opening wave grid file: '//TRIM(ErrMsg)
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
      ! write channel labels
      
      
      ! time
      WRITE(UnOut,"(A10)", IOSTAT=ErrStat2, advance="no") "Time"
   
      DO J = 1,5     !y
         DO K = 1,5  !x 
            WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " ze0", Num2Lstr(J+10*K)
         END DO
      END DO
      DO I=1,8          !z
         DO J = 1,5     !y
            DO K = 1,5  !x 
               WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " ux", Num2Lstr(I+10*J+100*K)
               WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " uy", Num2Lstr(I+10*J+100*K)
               WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " uz", Num2Lstr(I+10*J+100*K)
            END DO
         END DO
      END DO
      DO I=1,8          !z
         DO J = 1,5     !y
            DO K = 1,5  !x 
               WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " ax", Num2Lstr(I+10*J+100*K)
               WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " ay", Num2Lstr(I+10*J+100*K)
               WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " az", Num2Lstr(I+10*J+100*K)
            END DO
         END DO
      END DO
      DO I=1,8          !z
         DO J = 1,5     !y
            DO K = 1,5  !x 
               WRITE(UnOut,"(A3,A8)", IOSTAT=ErrStat2, advance="no") " pd", Num2Lstr(I+10*J+100*K)
            END DO
         END DO
      END DO
   
      ! end the line
      WRITE(UnOut, "(A1)", IOSTAT=ErrStat2) " "
      
      
      
      DO l=1, SIZE(InitInp%WaveTime)  ! loop through all time steps
      
         ! time
         WRITE(UnOut,"(F10.4)", IOSTAT=ErrStat2, advance="no") (l-1)*p%dtWave
         !WRITE(UnOut,"(F10.4)", IOSTAT=ErrStat2, advance="no") InitInp%WaveTime(l)
      
         ! wave elevation (all slices for now, to check)
         DO J = 1,5     !y
            DO K = 1,5  !x 
               Itemp = (J-1)*5.0 + K    ! index of actual node
               
               WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%zeta(l,J,K)
            END DO
         END DO
      
         ! wave velocities
         DO I=1,8          !z
            DO J = 1,5     !y
               DO K = 1,5  !x 
                  Itemp = (I-1)*25.0 + (J-1)*5.0 + K    ! index of actual node
                  
                  WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%ux(l,I,J,K)
                  WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%uy(l,I,J,K)
                  WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%uz(l,I,J,K)
               END DO
            END DO
         END DO
         
         ! wave accelerations
         DO I=1,8          !z
            DO J = 1,5     !y
               DO K = 1,5  !x 
                  Itemp = (I-1)*25.0 + (J-1)*5.0 + K    ! index of actual node
                  
                  WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%ax(l,I,J,K)
                  WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%ay(l,I,J,K)
                  WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%az(l,I,J,K)
               END DO
            END DO
         END DO
      
         ! dynamic pressure
         DO I=1,8          !z
            DO J = 1,5     !y
               DO K = 1,5  !x 
                  Itemp = (I-1)*25.0 + (J-1)*5.0 + K    ! index of actual node
                  
                  WRITE(UnOut,"(A1,e10.3)", IOSTAT=ErrStat2, advance="no") " ", p%PDyn(l,I,J,K)
               END DO
            END DO
         END DO
      
         ! end the line
         WRITE(UnOut, "(A1)", IOSTAT=ErrStat2) " "
      
      
      END DO
      
      
      CLOSE(UnOut, IOSTAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         ErrMsg = 'Error closing wave grid file'
      END IF
      
      
      
   !   Frmt = '(A10,'//TRIM(Int2LStr(p%NumOuts))//'(A1,A12))'
   !
   !   WRITE(p%MDUnOut,Frmt, IOSTAT=ErrStat2)  TRIM( 'Time' ), ( p%Delim, TRIM( p%OutParam(I)%Name), I=1,p%NumOuts )
   !
   !   WRITE(p%MDUnOut,Frmt)  TRIM( '(s)' ), ( p%Delim, TRIM( p%OutParam(I)%Units ), I=1,p%NumOuts )
   !
   !
   !
   !   ! Write the output parameters to the file
   !
   !   Frmt = '(F10.4,'//TRIM(Int2LStr(p%NumOuts))//'(A1,e10.4))' 
   !
   !   WRITE(p%MDUnOut,Frmt)  Time, ( p%Delim, y%WriteOutput(I), I=1,p%NumOuts )
   
      
      
      ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      


      ! if any of the coupled objects need initialization steps, that should have been taken care of already <<<<
      
      
      ! initialize objects with states, writing their initial states to the master state vector (x%states)
      
      
      !TODO: apply any initial adjustment of line length from active tensioning <<<<<<<<<<<<
      ! >>> maybe this should be skipped <<<<
		
      
		 ! Go through Bodys and write the coordinates to the state vector
      DO l = 1,p%nFreeBodies
         CALL Body_Initialize(m%BodyList(m%FreeBodyIs(l)), x%states(m%BodyStateIs1(l) : m%BodyStateIsN(l)), m)
      END DO
      
      ! Go through independent (including pinned) Rods and write the coordinates to the state vector
      DO l = 1,p%nFreeRods
         CALL Rod_Initialize(m%RodList(m%FreeRodIs(l)), x%states(m%RodStateIs1(l):m%RodStateIsN(l)), m%LineList)
      END DO

      ! Go through independent connections (Connects) and write the coordinates to the state vector and set positions of attached line ends
      DO l = 1, p%nFreeCons
         CALL Connect_Initialize(m%ConnectList(m%FreeConIs(l)), x%states(m%ConStateIs1(l) : m%conStateIsN(l)), m%LineList)
      END DO


      ! Lastly, go through lines and initialize internal node positions using quasi-static model
      DO l = 1, p%NLines

         N = m%LineList(l)%N ! for convenience

   !      ! set end node positions and velocities from connect objects
   !      m%LineList(l)%r(:,N) = m%ConnectList(m%LineList(l)%FairConnect)%r
   !      m%LineList(l)%r(:,0) = m%ConnectList(m%LineList(l)%AnchConnect)%r
   !      m%LineList(l)%rd(:,N) = (/ 0.0, 0.0, 0.0 /)  ! set anchor end velocities to zero
   !      m%LineList(l)%rd(:,0) = (/ 0.0, 0.0, 0.0 /)  ! set fairlead end velocities to zero

         ! set initial line internal node positions using quasi-static model or straight-line interpolation from anchor to fairlead
         CALL Line_Initialize( m%LineList(l), m%LineTypeList(m%LineList(l)%PropsIdNum), p%rhoW ,  ErrStat2, ErrMsg2)
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
            IF (ErrStat >= ErrId_Warn) CALL WrScr("   Catenary solver failed for one or more lines.  Using linear node spacing.")  ! make this statement more accurate

!         print *, "Line ", l, " with NumSegs =", N
!         print *, "its states range from index ", m%LineStateIs1(l), " to ", m%LineStateIsN(l)

         ! assign the resulting internal node positions to the integrator initial state vector! (velocities leave at 0)
         DO I = 1, N-1
!            print *, "I=", I
            DO J = 1, 3
!               print*, J, " ... writing position state to index ", 1*(m%LineStateIs1(l) + 3*N-3 + 3*I-3 + J-1)
               x%states(m%LineStateIs1(l) + 3*N-3 + 3*I-3 + J-1 ) = m%LineList(l)%r(J,I) ! assign position
               x%states(m%LineStateIs1(l)         + 3*I-3 + J-1 ) = 0.0_DbKi ! assign velocities (of zero)
            END DO
!            print *, m%LineList(l)%r(:,I)
         END DO

      END DO    !l = 1, p%NLines



      ! --------------------------------------------------------------------
      !          open output file(s) and write header lines
      CALL MDIO_OpenOutput( InitInp%FileName, p, m, InitOut, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
      ! --------------------------------------------------------------------


      
!      print *,"Done setup of the system (before any dynamic relaxation. State vector is as follows:"
      
!      DO I = 1, m%Nx
!         print *, x%states(I)


!      ! try writing output for troubleshooting purposes (TEMPORARY)
!      CALL MDIO_WriteOutputs(-1.0_DbKi, p, m, y, ErrStat, ErrMsg)
!      IF ( ErrStat >= AbortErrLev ) THEN
!         ErrMsg = ' Error in MDIO_WriteOutputs: '//TRIM(ErrMsg)
!         RETURN
!      END IF
!      END DO

      ! --------------------------------------------------------------------
      !           do dynamic relaxation to get ICs
      ! --------------------------------------------------------------------

      ! only do this if TMaxIC > 0
      if (InitInp%TMaxIC > 0.0_DbKi) then

         CALL WrScr("   Finalizing ICs using dynamic relaxation."//NewLine)  ! newline because next line writes over itself

         ! boost drag coefficient of each line type  <<<<<<<< does this actually do anything or do lines hold these coefficients???
         DO I = 1, p%nLineTypes
            m%LineTypeList(I)%Cdn = m%LineTypeList(I)%Cdn * InitInp%CdScaleIC
            m%LineTypeList(I)%Cdt = m%LineTypeList(I)%Cdt * InitInp%CdScaleIC   ! <<<<< need to update this to apply to all objects' drag
         END DO

         ! allocate array holding 10 latest fairlead tensions
         ALLOCATE ( FairTensIC(p%nLines, 10), STAT = ErrStat )
         IF ( ErrStat /= ErrID_None ) THEN
            CALL CheckError( ErrID_Fatal, ErrMsg2 )
            RETURN
         END IF

         ! initialize fairlead tension memory at changing values so things start unconverged
         DO J = 1,p%nLines
            DO I = 1, 10
               FairTensIC(J,I) = I
            END DO
         END DO


         ! round dt to integer number of time steps  
         NdtM = ceiling(InitInp%DTIC/p%dtM0)            ! get number of mooring time steps to do based on desired time step size
         dtM = InitInp%DTIC/real(NdtM, DbKi)            ! adjust desired time step to satisfy dt with an integer number of time steps

         t = 0.0_DbKi     ! start time at zero

         ! because TimeStep wants an array...
         call MD_CopyInput( u, u_array(1), MESH_NEWCOPY, ErrStat2, ErrMsg2 )  ! make a size=1 array of inputs (since MD_RK2 expects an array to InterpExtrap)
         call MD_CopyInput( u,  u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  ! also make an inputs object to interpExtrap to
         t_array(1) = t                                                       ! fill in the times "array" for u_array

         DO I = 1, ceiling(InitInp%TMaxIC/InitInp%DTIC)   ! loop through IC gen time steps, up to maximum


            !loop through line integration time steps
            DO J = 1, NdtM                                 ! for (double ts=t; ts<=t+ICdt-dts; ts+=dts)

               CALL MD_RK2(t, dtM, u_interp, u_array, t_array, p, x, xd, z, other, m, ErrStat2, ErrMsg2)
                              
               ! check for NaNs - is this a good place/way to do it?
               DO K = 1, m%Nx
                  IF (Is_NaN(REAL(x%states(K),DbKi))) THEN
                     ErrStat = ErrID_Fatal
                     ErrMsg = ' NaN state detected.'
                     EXIT
                  END IF
               END DO
               
               IF (ErrStat == ErrID_Fatal) THEN
                  print *, "NaN detected at time ", t, " during IC gen. Here is the state vector: "
                  print *, x%states
                  EXIT
               END IF

            END DO  ! J  time steps

         !   ! integrate the EOMs one DTIC s time step
         !   CALL TimeStep ( t, InitInp%DTIC, u_array, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg )
         !      CALL CheckError( ErrStat2, ErrMsg2 )
         !      IF (ErrStat >= AbortErrLev) RETURN

            ! store new fairlead tension (and previous fairlead tensions for comparison)
            DO l = 1, p%nLines
            
               DO K=0,8   ! we want to count down from 10 to 2 .
                  FairTensIC(l, 10-K) = FairTensIC(l, 9-K)   ! this pushes stored values up in the array
               END DO
                  
               ! now store latest value of each line's fairlead (end B) tension   
               FairTensIC(l,1) = TwoNorm(m%LineList(l)%Fnet(:, m%LineList(l)%N))
            END DO


            ! provide status message
            ! bjj: putting this in a string so we get blanks to cover up previous values (if current string is shorter than previous one)
            Message = '   t='//trim(Num2LStr(t))//'  FairTen 1: '//trim(Num2LStr(FairTensIC(1,1)))// &
                           ', '//trim(Num2LStr(FairTensIC(1,2)))//', '//trim(Num2LStr(FairTensIC(1,3))) 
            CALL WrOver( Message )

            ! check for convergence (compare current tension at each fairlead with previous 9 values)
            IF (I > 9) THEN
               
               Converged = 1
               
               ! check for non-convergence
               
               DO l = 1, p%nLines   
                  DO K = 1,9
                     IF ( abs( FairTensIC(l,K)/FairTensIC(l,K+1) - 1.0 ) > InitInp%threshIC ) THEN
                        Converged = 0
                        EXIT
                     END IF
                  END DO
                  
                  IF (Converged == 0) EXIT   ! make sure we exit this loop too
               END DO

               IF (Converged == 1)  THEN  ! if we made it with all cases satisfying the threshold
                  CALL WrScr('   Fairlead tensions converged to '//trim(Num2LStr(100.0*InitInp%threshIC))//'% after '//trim(Num2LStr(t))//' seconds.')
                  EXIT  ! break out of the time stepping loop
               END IF
            END IF

            IF (I == ceiling(InitInp%TMaxIC/InitInp%DTIC) ) THEN
               CALL WrScr('   Fairlead tensions did not converge within TMaxIC='//trim(Num2LStr(InitInp%TMaxIC))//' seconds.')
               !ErrStat = ErrID_Warn
               !ErrMsg = '  MD_Init: ran dynamic convergence to TMaxIC without convergence'
            END IF

         END DO ! I ... looping through time steps



         CALL MD_DestroyInput( u_array(1), ErrStat2, ErrMsg2 )

         ! UNboost drag coefficient of each line type   <<<
         DO I = 1, p%nLineTypes
            m%LineTypeList(I)%Cdn = m%LineTypeList(I)%Cdn / InitInp%CdScaleIC
            m%LineTypeList(I)%Cdt = m%LineTypeList(I)%Cdt / InitInp%CdScaleIC
         END DO

      end if ! InitInp%TMaxIC > 0
      

      p%dtCoupling = DTcoupling  ! store coupling time step for use in updatestates

      other%dummy = 0
      xd%dummy    = 0
      z%dummy     = 0      
      
      
      ! TODO: add feature for automatic water depth increase based on max anchor depth!

   CONTAINS


      LOGICAL FUNCTION AllocateFailed(arrayName)

         CHARACTER(*), INTENT(IN   )      :: arrayName     ! The array name
         
         call SetErrStat(ErrStat2, "Error allocating space for "//trim(arrayName)//" array.", ErrStat, ErrMsg, 'MD_Init') 
         AllocateFailed = ErrStat2 >= AbortErrLev
         if (AllocateFailed) call CleanUp() !<<<<<<<<<< need to fix this up
      END FUNCTION AllocateFailed
      
      
      LOGICAL FUNCTION Failed()

         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MDIO_ReadInput') 
         Failed = ErrStat >= AbortErrLev
         if (Failed) call CleanUp()
      END FUNCTION Failed


      SUBROUTINE CheckError(ErrID,Msg)
         ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev

            ! Passed arguments
         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
         CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)

         ! Set error status/message;
         IF ( ErrID /= ErrID_None ) THEN

            IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine   ! if there's a pre-existing warning/error, retain the message and start a new line

            ErrMsg = TRIM(ErrMsg)//' MD_Init:'//TRIM(Msg)
            ErrStat = MAX(ErrStat, ErrID)

            ! Clean up if we're going to return on error: close files, deallocate local arrays


            IF ( ErrStat >= AbortErrLev ) THEN                
               IF (ALLOCATED(m%CpldConIs        ))  DEALLOCATE(m%CpldConIs       )
               IF (ALLOCATED(m%FreeConIs       ))  DEALLOCATE(m%FreeConIs       )
               IF (ALLOCATED(m%LineStateIs1     ))  DEALLOCATE(m%LineStateIs1     )
               IF (ALLOCATED(m%LineStateIsN     ))  DEALLOCATE(m%LineStateIsN     )
               IF (ALLOCATED(m%ConStateIs1      ))  DEALLOCATE(m%ConStateIs1     )
               IF (ALLOCATED(m%ConStateIsN      ))  DEALLOCATE(m%ConStateIsN     )
               IF (ALLOCATED(x%states           ))  DEALLOCATE(x%states           )
               IF (ALLOCATED(FairTensIC         ))  DEALLOCATE(FairTensIC         ) 
            END IF
         END IF

      END SUBROUTINE CheckError

     SUBROUTINE CleanUp()
        ! ErrStat = ErrID_Fatal  
        CLOSE( UnIn )
 !       IF (InitInp%Echo) CLOSE( UnEc )
     END SUBROUTINE

   END SUBROUTINE MD_Init
   !----------------------------------------------------------------------------------------======




   !----------------------------------------------------------------------------------------======
   SUBROUTINE MD_UpdateStates( t, n, u, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg)

      REAL(DbKi)                      , INTENT(IN   ) :: t
      INTEGER(IntKi)                  , INTENT(IN   ) :: n
      TYPE(MD_InputType)              , INTENT(INOUT) :: u(:)       ! INTENT(INOUT) ! had to change this to INOUT
      REAL(DbKi)                      , INTENT(IN   ) :: t_array(:)
      TYPE(MD_ParameterType)          , INTENT(IN   ) :: p          ! INTENT(IN   )
      TYPE(MD_ContinuousStateType)    , INTENT(INOUT) :: x          ! INTENT(INOUT)
      TYPE(MD_DiscreteStateType)      , INTENT(INOUT) :: xd         ! INTENT(INOUT)
      TYPE(MD_ConstraintStateType)    , INTENT(INOUT) :: z          ! INTENT(INOUT)
      TYPE(MD_OtherStateType)         , INTENT(INOUT) :: other      ! INTENT(INOUT)
      TYPE(MD_MiscVarType)            , INTENT(INOUT) :: m          ! INTENT(INOUT)
      INTEGER(IntKi)                  , INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*)                    , INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)                                  :: ErrStat2   ! Error status of the operation
      CHARACTER(ErrMsgLen)                            :: ErrMsg2    ! Error message if ErrStat2 /= ErrID_None

! moved to TimeStep      TYPE(MD_InputType)                              :: u_interp   !
      INTEGER(IntKi)                                  :: nTime
  
      TYPE(MD_InputType)                                  :: u_interp   ! interpolated instantaneous input values to be calculated for each mooring time step

      REAL(DbKi)                                      :: t2          ! copy of time variable that will get advanced by the integrator (not sure this is necessary<<<)
      REAL(DbKi)                                      :: dtM         ! actual mooring dynamics time step
      INTEGER(IntKi)                                  :: NdtM        ! number of time steps to integrate through with RK2
      INTEGER(IntKi)                                  :: I
      INTEGER(IntKi)                                  :: J

      nTime = size(u) ! the number of times of input data provided? <<<<<<< not used

      t2 = t

! >>> removing this section and putting it inside loop of TimeStep (to be done at every time step) <<<
!      ! create space for arrays/meshes in u_interp
!      CALL MD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
!         CALL CheckError( ErrStat2, ErrMsg2 )
!         IF (ErrStat >= AbortErrLev) RETURN
!
!      ! interpolate input mesh to correct time
!      CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t, ErrStat2, ErrMsg2)
!         CALL CheckError( ErrStat2, ErrMsg2 )
!         IF (ErrStat >= AbortErrLev) RETURN
!
!
!      ! go through fairleads and apply motions from driver
!      DO I = 1, p%nCpldCons
!         DO J = 1,3
!            m%ConnectList(m%CpldConIs(I))%r(J)  = u_interp%PtFairleadDisplacement%Position(J,I) + u_interp%PtFairleadDisplacement%TranslationDisp(J,I)
!            m%ConnectList(m%CpldConIs(I))%rd(J) = u_interp%PtFairleadDisplacement%TranslationVel(J,I)  ! is this right? <<<
!         END DO
!      END DO
!



!      ! call function that loops through mooring model time steps
!      CALL TimeStep ( t2, p%dtCoupling, u, t_array, p, x, xd, z, other, m, ErrStat2, ErrMsg2 )
!         CALL CheckError( ErrStat2, ErrMsg2 )
!         IF (ErrStat >= AbortErrLev) RETURN

         
      ! create space for arrays/meshes in u_interp   ... is it efficient to do this every time step???
      CALL MD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg)


      ! round dt to integer number of time steps   <<<< should this be calculated only once, up front?
      NdtM = ceiling(p%dtCoupling/p%dtM0)            ! get number of mooring time steps to do based on desired time step size
      dtM = p%dtCoupling/float(NdtM)                       ! adjust desired time step to satisfy dt with an integer number of time steps


      !loop through line integration time steps
      DO I = 1, NdtM                                 ! for (double ts=t; ts<=t+ICdt-dts; ts+=dts)

         CALL MD_RK2(t2, dtM, u_interp, u, t_array, p, x, xd, z, other, m, ErrStat2, ErrMsg2)
         
         
         ! check for NaNs - is this a good place/way to do it?
         DO J = 1, m%Nx
            IF (Is_NaN(REAL(x%states(J),DbKi))) THEN
               ErrStat = ErrID_Fatal
               ErrMsg = ' NaN state detected.'
               EXIT
            END IF
         END DO
         
         IF (ErrStat == ErrID_Fatal) THEN
            print *, "NaN detected at time ", t2, ". Here is the state vector: "
            print *, x%states
            EXIT
         END IF
      
      END DO  ! I  time steps


      ! destroy dxdt and x2, and u_interp
      !CALL MD_DestroyContState( dxdt, ErrStat, ErrMsg)
      !CALL MD_DestroyContState( x2, ErrStat, ErrMsg)
      CALL MD_DestroyInput(u_interp, ErrStat, ErrMsg)
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error destroying dxdt or x2.'
      END IF
      !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_UpdateStates')


      ! check for NaNs - is this a good place/way to do it?
      DO J = 1, m%Nx
         IF (Is_NaN(REAL(x%states(J),DbKi))) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' NaN state detected.'
            EXIT
         END IF
      END DO
      
      IF (ErrStat == ErrID_Fatal) THEN
         print *, "NaN detected at time ", t2, ". Here is the state vector: "
         print *, x%states
      END IF

   CONTAINS

      SUBROUTINE CheckError(ErrId, Msg)
        ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev

         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


         IF ( ErrID /= ErrID_None ) THEN

            IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine  ! keep existing error message if there is one
            ErrMsg = TRIM(ErrMsg)//' MD_UpdateStates:'//TRIM(Msg)      ! add current error message
            ErrStat = MAX(ErrStat, ErrID)

            CALL WrScr( ErrMsg )  ! do this always or only if warning level?

            IF( ErrStat > ErrID_Warn ) THEN
       !         CALL MD_DestroyInput( u_interp, ErrStat, ErrMsg )
                RETURN
            END IF
         END IF

      END SUBROUTINE CheckError

   END SUBROUTINE MD_UpdateStates
   !----------------------------------------------------------------------------------------



   !----------------------------------------------------------------------------------------
   SUBROUTINE MD_CalcOutput( t, u, p, x, xd, z, other, y, m, ErrStat, ErrMsg )

      REAL(DbKi)                     , INTENT(IN   ) :: t
      TYPE( MD_InputType )           , INTENT(IN   ) :: u       ! INTENT(IN   )
      TYPE( MD_ParameterType )       , INTENT(IN   ) :: p       ! INTENT(IN   )
      TYPE( MD_ContinuousStateType ) , INTENT(IN   ) :: x       ! INTENT(IN   )
      TYPE( MD_DiscreteStateType )   , INTENT(IN   ) :: xd      ! INTENT(IN   )
      TYPE( MD_ConstraintStateType ) , INTENT(IN   ) :: z       ! INTENT(IN   )
      TYPE( MD_OtherStateType )      , INTENT(IN   ) :: other   ! INTENT(IN   )
      TYPE( MD_OutputType )          , INTENT(INOUT) :: y       ! INTENT(INOUT)
      TYPE(MD_MiscVarType)           , INTENT(INOUT) :: m       ! INTENT(INOUT)
      INTEGER(IntKi)                 , INTENT(INOUT) :: ErrStat
      CHARACTER(*)                   , INTENT(INOUT) :: ErrMsg

   !   TYPE(MD_ContinuousStateType)                   :: dxdt    ! time derivatives of continuous states (initialized in CalcContStateDeriv)
      INTEGER(IntKi)                                 :: I         ! counter
      INTEGER(IntKi)                                 :: J         ! counter
      INTEGER(IntKi)                                 :: K         ! counter
      INTEGER(IntKi)                                 :: l         ! index used for objects
      
      Real(DbKi)                                     :: F6net(6)  ! net force and moment calculated on coupled objects
    
      INTEGER(IntKi)                                 :: ErrStat2  ! Error status of the operation
      CHARACTER(ErrMsgLen)                           :: ErrMsg2   ! Error message if ErrStat2 /= ErrID_None


      ! below updated to make sure outputs are current (based on provided x and u)  - similar to what's in UpdateStates

    !  ! go through fairleads and apply motions from driver
    !  DO I = 1, p%nCpldCons
    !     DO J = 1,3
    !        m%ConnectList(m%CpldConIs(I))%r(J)  = u%CoupledKinematics%Position(J,I) + u%CoupledKinematics%TranslationDisp(J,I)
    !        m%ConnectList(m%CpldConIs(I))%rd(J) = u%CoupledKinematics%TranslationVel(J,I)  ! is this right? <<<
    !     END DO
    !  END DO
      
      
    ! ! go through nodes and apply wave kinematics from driver
    ! IF (p%WaterKin > 0) THEN
    ! 
    !    J=0
    !    ! Body reference point coordinates
    !    DO I = 1, p%nBodies
    !       J = J + 1                     
    !       m%BodyList(I)%U    = u%U(:,J)
    !       m%BodyList(I)%Ud   = u%Ud(:,J)
    !       m%BodyList(I)%zeta = u%zeta(J)
    !    END DO
    !    ! Rod node coordinates
    !    DO I = 1, p%nRods
    !       DO K = 0,m%RodList(I)%N  
    !          J = J + 1             
    !          m%RodList(I)%U (:,K) = u%U(:,J)
    !          m%RodList(I)%Ud(:,K) = u%Ud(:,J)
    !          m%RodList(I)%zeta(K) = u%zeta(J)
    !          m%RodList(I)%PDyn(K) = u%PDyn(J)
    !       END DO
    !    END DO
    !    ! Point reference point coordinates
    !    DO I = 1, p%nConnects
    !       J = J + 1
    !       m%ConnectList(I)%U    = u%U(:,J)
    !       m%ConnectList(I)%Ud   = u%Ud(:,J)
    !       m%ConnectList(I)%zeta = u%zeta(J)
    !    END DO      
    !    ! Line internal node coordinates
    !    DO I = 1, p%nLines
    !       DO K = 1, m%LineList(I)%N-1
    !          J = J + 1               
    !          m%LineList(I)%U (:,K) = u%U(:,J)
    !          m%LineList(I)%Ud(:,K) = u%Ud(:,J)
    !          m%LineList(I)%zeta(K) = u%zeta(J)
    !       END DO
    !    END DO   
    ! 
    ! END IF
      
      

      ! call CalcContStateDeriv in order to run model and calculate dynamics with provided x and u
      CALL MD_CalcContStateDeriv( t, u, p, x, xd, z, other, m, m%xdTemp, ErrStat, ErrMsg )

    !  ! assign net force on fairlead Connects to the fairlead force output mesh
    !  DO i = 1, p%nCpldCons
    !     DO J=1,3
    !        y%PtFairleadLoad%Force(J,I) = m%ConnectList(m%CpldConIs(I))%Fnet(J)
    !     END DO
    !  END DO
      
      ! now that forces have been updated, write them to the output mesh
      J = 0
      DO l = 1,p%nCpldBodies
         J = J + 1
         CALL Body_GetCoupledForce(m%BodyList(m%CpldBodyIs(l)), F6net, m, p)
         y%CoupledLoads%Force( :,J) = F6net(1:3)
         y%CoupledLoads%Moment(:,J) = F6net(4:6)
      END DO
            
      DO l = 1,p%nCpldRods
         J = J + 1
         CALL Rod_GetCoupledForce(m%RodList(m%CpldRodIs(l)), F6net, m%LineList, p)
         y%CoupledLoads%Force( :,J) = F6net(1:3)
         y%CoupledLoads%Moment(:,J) = F6net(4:6)
      END DO
      
      DO l = 1,p%nCpldCons
         J = J + 1
         CALL Connect_GetCoupledForce(m%ConnectList(m%CpldConIs(l)), F6net(1:3), m%LineList, p)
         y%CoupledLoads%Force(:,J) = F6net(1:3)
      END DO
      
      
      
      

   !  ! write all node positions to the node positons output array
   !  ! go through the nodes and fill in the data (this should maybe be turned into a global function)
   !  J=0
   !  ! Body reference point coordinates
   !  DO I = 1, p%nBodies
   !     J = J + 1                     
   !     y%rAll(:,J) = m%BodyList(I)%r6(1:3)         
   !  END DO
   !  ! Rod node coordinates
   !  DO I = 1, p%nRods
   !     DO K = 0,m%RodList(I)%N  
   !        J = J + 1             
   !        y%rAll(:,J) = m%RodList(I)%r(:,K)
   !     END DO
   !  END DO
   !  ! Point reference point coordinates
   !  DO I = 1, p%nConnects
   !     J = J + 1
   !     y%rAll(:,J) = m%ConnectList(I)%r
   !  END DO      
   !  ! Line internal node coordinates
   !  DO I = 1, p%nLines
   !     DO K = 1, m%LineList(I)%N-1
   !        J = J + 1               
   !        y%rAll(:,J) = m%LineList(I)%r(:,K)
   !     END DO
   !  END DO   
   

      ! calculate outputs (y%WriteOutput) for glue code and write any m outputs to MoorDyn output files
      CALL MDIO_WriteOutputs(REAL(t,DbKi) , p, m, y, ErrStat2, ErrMsg2)
      CALL CheckError(ErrStat2, 'In MDIO_WriteOutputs: '//trim(ErrMsg2))
      IF ( ErrStat >= AbortErrLev ) RETURN


  !    ! destroy dxdt
  !    CALL MD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
  !    CALL CheckError(ErrStat2, 'When destroying dxdt: '//trim(ErrMsg2))
  !    IF ( ErrStat >= AbortErrLev ) RETURN



   CONTAINS

      SUBROUTINE CheckError(ErrId, Msg)
        ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev

         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


         IF ( ErrID /= ErrID_None ) THEN

            IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine  ! keep existing error message if there is one
            ErrMsg = TRIM(ErrMsg)//' MD_CalcOutput:'//TRIM(Msg)      ! add current error message
            ErrStat = MAX(ErrStat, ErrID)

            CALL WrScr( ErrMsg )  ! do this always or only if warning level? <<<<<<<<<<<<<<<<<<<<<< probably should remove all instances

      !      IF( ErrStat > ErrID_Warn ) THEN
      !          CALL MD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
      !      END IF
         END IF

      END SUBROUTINE CheckError

   END SUBROUTINE MD_CalcOutput
   !----------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------
   SUBROUTINE MD_CalcContStateDeriv( t, u, p, x, xd, z, other, m, dxdt, ErrStat, ErrMsg )
   ! Tight coupling routine for computing derivatives of continuous states
   ! this is modelled off what used to be subroutine DoRHSmaster

      REAL(DbKi),                         INTENT(IN )    :: t       ! Current simulation time in seconds
      TYPE(MD_InputType),                 INTENT(IN )    :: u       ! Inputs at t
      TYPE(MD_ParameterType),             INTENT(IN )    :: p       ! Parameters
      TYPE(MD_ContinuousStateType),       INTENT(IN )    :: x       ! Continuous states at t
      TYPE(MD_DiscreteStateType),         INTENT(IN )    :: xd      ! Discrete states at t
      TYPE(MD_ConstraintStateType),       INTENT(IN )    :: z       ! Constraint states at t
      TYPE(MD_OtherStateType),            INTENT(IN )    :: other   ! Other states at t
      TYPE(MD_MiscVarType),               INTENT(INOUT)  :: m       ! misc/optimization variables
      TYPE(MD_ContinuousStateType),       INTENT(INOUT)  :: dxdt    ! Continuous state derivatives at t
      INTEGER(IntKi),                     INTENT( OUT)   :: ErrStat ! Error status of the operation
      CHARACTER(*),                       INTENT( OUT)   :: ErrMsg  ! Error message if ErrStat /= ErrID_None


      INTEGER(IntKi)                                     :: L       ! index
      INTEGER(IntKi)                                     :: I       ! index
      INTEGER(IntKi)                                     :: J       ! index
      INTEGER(IntKi)                                     :: K       ! index
      INTEGER(IntKi)                                     :: Istart  ! start index of line/connect in state vector
      INTEGER(IntKi)                                     :: Iend    ! end index of line/connect in state vector

      REAL(DbKi)                                         :: r6_in(6) ! temporary for passing kinematics
      REAL(DbKi)                                         :: v6_in(6) ! temporary for passing kinematics
      REAL(DbKi)                                         :: a6_in(6) ! temporary for passing kinematics
      REAL(DbKi)                                         :: r_in(3)  ! temporary for passing kinematics
      REAL(DbKi)                                         :: rd_in(3) ! temporary for passing kinematics
      REAL(DbKi)                                         :: a_in(3)  ! temporary for passing kinematics

      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg = ""

!      ! allocations of dxdt (as in SubDyn.  "INTENT(OUT) automatically deallocates the arrays on entry, we have to allocate them here"  is this right/efficient?)
!      ALLOCATE ( dxdt%states(size(x%states)), STAT = ErrStat )
!      IF ( ErrStat /= ErrID_None ) THEN
!          ErrMsg  = ' Error allocating dxdt%states array.'
!          RETURN
!      END IF

      ! clear connection force and mass values <M<<<<<<<<<<<<<<<<<<<<<<<
      DO L = 1, p%NConnects
        DO J = 1,3
          m%ConnectList(L)%Fnet(J) = 0.0_DbKi
          m%ConnectList(L)%Fnet(J) = 0.0_DbKi
          DO K = 1,3
            m%ConnectList(L)%M   (K,J) = 0.0_DbKi
            m%ConnectList(L)%M   (K,J) = 0.0_DbKi
          END DO
        END DO
      END DO


      ! call ground body to update all the fixed things...
      !GroundBody->updateFairlead( t );                  <<<< manually set anchored connection stuff for now here
      r6_in = 0.0_DbKi
      v6_in = 0.0_DbKi
      CALL Body_SetKinematics(m%GroundBody, r6_in, v6_in, m%zeros6, t, m)
      
      ! ---------------------------------- coupled things ---------------------------------
      ! Apply displacement and velocity terms here. Accelerations will be considered to calculate inertial loads at the end.      
      
      J = 0  ! J is the index of the coupling points in the input mesh CoupledKinematics
      ! any coupled bodies (type -1)
      DO l = 1,p%nCpldBodies
         J = J + 1
         r6_in(1:3) = u%CoupledKinematics%Position(:,J) + u%CoupledKinematics%TranslationDisp(:,J)
         r6_in(4:6) = EulerExtract( TRANSPOSE( u%CoupledKinematics%Orientation(:,:,J) ) )
         v6_in(1:3) = u%CoupledKinematics%TranslationVel(:,J)
         v6_in(4:6) = u%CoupledKinematics%RotationVel(:,J)
         a6_in(1:3) = u%CoupledKinematics%TranslationAcc(:,J)
         a6_in(4:6) = u%CoupledKinematics%RotationAcc(:,J)
      
      
         if ((t >= 1.0) .and. (t < 1.001)) then
            print *, "orientation matrix from mesh:"
            print *, u%CoupledKinematics%Orientation(:,:,J)
            print *, "Euler extract:"
            print *, EulerExtract( u%CoupledKinematics%Orientation(:,:,J) )
            print *, "Euler extract of transpose"
            print *, EulerExtract( transpose(u%CoupledKinematics%Orientation(:,:,J) ))
            
            print *, "done"
         end if
      
         CALL Body_SetKinematics(m%BodyList(m%CpldBodyIs(l)), r6_in, v6_in, a6_in, t, m)
      END DO
      
      ! any coupled rods (type -1 or -2)    note, rotations ignored if it's a pinned rod
      DO l = 1,p%nCpldRods
         J = J + 1

         r6_in(1:3) = u%CoupledKinematics%Position(:,J) + u%CoupledKinematics%TranslationDisp(:,J)
         r6_in(4:6) = EulerExtract( TRANSPOSE( u%CoupledKinematics%Orientation(:,:,J) ) )  ! <<<< should make sure this works <<<
         v6_in(1:3) = u%CoupledKinematics%TranslationVel(:,J)
         v6_in(4:6) = u%CoupledKinematics%RotationVel(:,J)
         a6_in(1:3) = u%CoupledKinematics%TranslationAcc(:,J)
         a6_in(4:6) = u%CoupledKinematics%RotationAcc(:,J)
         
         if ((t >= 1.0) .and. (t < 1.001)) then
            print *, "orientation matrix from mesh:"
            print *, u%CoupledKinematics%Orientation(:,:,J)
            print *, "Euler extract:"
            print *, EulerExtract( u%CoupledKinematics%Orientation(:,:,J) )
            print *, "Euler extract of transpose"
            print *, EulerExtract( transpose(u%CoupledKinematics%Orientation(:,:,J) ))
            
            print *, "done"
         end if
      
         CALL Rod_SetKinematics(m%RodList(m%CpldRodIs(l)), r6_in, v6_in, a6_in, t, m%LineList)
 		
      END DO
      
      ! any coupled points (type -1)
      DO l = 1, p%nCpldCons
         J = J + 1
         
         r_in  = u%CoupledKinematics%Position(:,J) + u%CoupledKinematics%TranslationDisp(:,J)
         rd_in = u%CoupledKinematics%TranslationVel(:,J)
         a_in(1:3) = u%CoupledKinematics%TranslationAcc(:,J)
         CALL Connect_SetKinematics(m%ConnectList(m%CpldConIs(l)), r_in, rd_in, a_in, t, m%LineList)
         
         !print *, u%PtFairleadDisplacement%Position(:,l) + u%PtFairleadDisplacement%TranslationDisp(:,l)
         !print *, u%PtFairleadDisplacement%TranslationVel(:,l)
         
      END DO
      
      
      ! apply line length changes from active tensioning if applicable
      DO L = 1, p%NLines
         IF (m%LineList(L)%CtrlChan > 0) then

            ! do a bounds check to prohibit excessive segment length changes (until a method to add/remove segments is created)
            IF ( u%DeltaL(m%LineList(L)%CtrlChan) > m%LineList(L)%UnstrLen / m%LineList(L)%N ) then
                ErrStat = ErrID_Fatal
                ErrMsg  = ' Active tension command will make a segment longer than the limit of twice its original length.'
                print *, u%DeltaL(m%LineList(L)%CtrlChan), " is an increase of more than ", (m%LineList(L)%UnstrLen / m%LineList(L)%N)
                print *, u%DeltaL
                print*, m%LineList(L)%CtrlChan
                RETURN
            END IF
            IF ( u%DeltaL(m%LineList(L)%CtrlChan) < -0.5 * m%LineList(L)%UnstrLen / m%LineList(L)%N ) then
             ErrStat = ErrID_Fatal
                ErrMsg  = ' Active tension command will make a segment shorter than the limit of half its original length.'
                print *, u%DeltaL(m%LineList(L)%CtrlChan), " is a reduction of more than half of ", (m%LineList(L)%UnstrLen / m%LineList(L)%N)
                print *, u%DeltaL
                print*, m%LineList(L)%CtrlChan
                RETURN
            END IF                

            ! for now this approach only acts on the fairlead end segment, and assumes all segment lengths are otherwise equal size
            m%LineList(L)%l( m%LineList(L)%N) = m%LineList(L)%UnstrLen/m%LineList(L)%N + u%DeltaL(m%LineList(L)%CtrlChan)       
            m%LineList(L)%ld(m%LineList(L)%N) =                                       u%DeltaLdot(m%LineList(L)%CtrlChan)       
         END IF
      END DO      
      
      
   !  ! go through nodes and apply wave kinematics from driver
   !  IF (p%WaterKin > 0) THEN
   !  
   !     J=0
   !     ! Body reference point coordinates
   !     DO I = 1, p%nBodies
   !        J = J + 1                     
   !        m%BodyList(I)%U    = u%U(:,J)
   !        m%BodyList(I)%Ud   = u%Ud(:,J)
   !        m%BodyList(I)%zeta = u%zeta(J)
   !     END DO
   !     ! Rod node coordinates
   !     DO I = 1, p%nRods
   !        DO K = 0,m%RodList(I)%N  
   !           J = J + 1             
   !           m%RodList(I)%U (:,K) = u%U(:,J)
   !           m%RodList(I)%Ud(:,K) = u%Ud(:,J)
   !           m%RodList(I)%zeta(K) = u%zeta(J)
   !           m%RodList(I)%PDyn(K) = u%PDyn(J)
   !        END DO
   !     END DO
   !     ! Point reference point coordinates
   !     DO I = 1, p%nConnects
   !        J = J + 1
   !        m%ConnectList(I)%U    = u%U(:,J)
   !        m%ConnectList(I)%Ud   = u%Ud(:,J)
   !        m%ConnectList(I)%zeta = u%zeta(J)
   !     END DO      
   !     ! Line internal node coordinates
   !     DO I = 1, p%nLines
   !        DO K = 1, m%LineList(I)%N-1
   !           J = J + 1               
   !           m%LineList(I)%U (:,K) = u%U(:,J)
   !           m%LineList(I)%Ud(:,K) = u%Ud(:,J)
   !           m%LineList(I)%zeta(K) = u%zeta(J)
   !        END DO
   !     END DO   
   !  
   !  END IF
      
      
      ! independent or semi-independent things with their own states...
      
      ! give Bodies latest state variables (kinematics will also be assigned to dependent connections and rods, and thus line ends)
      DO l = 1,p%nFreeBodies
         CALL Body_SetState(m%BodyList(m%FreeBodyIs(l)), x%states(m%BodyStateIs1(l):m%BodyStateIsN(l)), t, m)
      END DO
      
      ! give independent or pinned rods' latest state variables (kinematics will also be assigned to attached line ends)
      DO l = 1,p%nFreeRods
         CALL Rod_SetState(m%RodList(m%FreeRodIs(l)), x%states(m%RodStateIs1(l):m%RodStateIsN(l)), t, m%LineList)
      END DO
      
      ! give Connects (independent connections) latest state variable values (kinematics will also be assigned to attached line ends)
      DO l = 1,p%nFreeCons
  !       Print *, "calling SetState for free connection, con#", m%FreeConIs(l), " with state range: ", m%ConStateIs1(l), "-", m%ConStateIsN(l)
         !K=K+1
         CALL Connect_SetState(m%ConnectList(m%FreeConIs(l)), x%states(m%ConStateIs1(l):m%ConStateIsN(l)), t, m%LineList)
      END DO
      
      ! give Lines latest state variable values for internal nodes
      DO l = 1,p%nLines
         CALL Line_SetState(m%LineList(l), x%states(m%LineStateIs1(l):m%LineStateIsN(l)), t)
      END DO

      ! calculate dynamics of free objects (will also calculate forces (doRHS()) from any child/dependent objects)...
         
      ! calculate line dynamics (and calculate line forces and masses attributed to connections)
      DO l = 1,p%nLines
         CALL Line_GetStateDeriv(m%LineList(l), dxdt%states(m%LineStateIs1(l):m%LineStateIsN(l)), p)  !dt might also be passed for fancy friction models
      END DO
      
      ! calculate connect dynamics (including contributions from attached lines
      ! as well as hydrodynamic forces etc. on connect object itself if applicable)
      DO l = 1,p%nFreeCons
         CALL Connect_GetStateDeriv(m%ConnectList(m%FreeConIs(l)), dxdt%states(m%ConStateIs1(l):m%ConStateIsN(l)), m%LineList, p)
      END DO
      
      ! calculate dynamics of independent Rods 
      DO l = 1,p%nFreeRods
         CALL Rod_GetStateDeriv(m%RodList(m%FreeRodIs(l)), dxdt%states(m%RodStateIs1(l):m%RodStateIsN(l)), m%LineList, p)
      END DO
      
      ! calculate dynamics of Bodies
      DO l = 1,p%nFreeBodies
         CALL Body_GetStateDeriv(m%BodyList(m%FreeBodyIs(l)), dxdt%states(m%BodyStateIs1(l):m%BodyStateIsN(l)), m, p)
      END DO
      
      
      
      ! get dynamics/forces (doRHS()) of coupled objects, which weren't addressed in above calls (this includes inertial loads)
      ! note: can do this in any order since there are no dependencies among coupled objects
      
      
      DO l = 1,p%nCpldCons
      
 !        >>>>>>>> here we should pass along accelerations and include inertial loads in the calculation!!! <<<??
 !               in other words are the below good enough or do I need to call _getCoupledFOrce??
      
         CALL Connect_DoRHS(m%ConnectList(m%CpldConIs(l)), m%LineList, p)
      END DO
      
      DO l = 1,p%nCpldRods
         CALL Rod_DoRHS(m%RodList(m%CpldRodIs(l)), m%LineList, p)
         ! NOTE: this won't compute net loads on Rod. Need Rod_GetNetForceAndMass for that. Change? <<<<
      END DO
      
      DO l = 1,p%nCpldBodies
         CALL Body_DoRHS(m%BodyList(m%CpldBodyIs(l)), m, p)
      END DO

      ! call ground body to update all the fixed things
      CALL Body_DoRHS(m%GroundBody, m, p)
      

   
   END SUBROUTINE MD_CalcContStateDeriv
   !----------------------------------------------------------------------------------------=====


   !----------------------------------------------------------------------------------------=======
   SUBROUTINE MD_End(u, p, x, xd, z, other, y, m, ErrStat , ErrMsg)
   
      TYPE(MD_InputType) ,            INTENT(INOUT) :: u
      TYPE(MD_ParameterType) ,        INTENT(INOUT) :: p
      TYPE(MD_ContinuousStateType) ,  INTENT(INOUT) :: x
      TYPE(MD_DiscreteStateType) ,    INTENT(INOUT) :: xd
      TYPE(MD_ConstraintStateType) ,  INTENT(INOUT) :: z
      TYPE(MD_OtherStateType) ,       INTENT(INOUT) :: other
      TYPE(MD_OutputType) ,           INTENT(INOUT) :: y
      TYPE(MD_MiscVarType),           INTENT(INOUT) :: m      
      INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat
      CHARACTER(*),                   INTENT(  OUT) :: ErrMsg

!      INTEGER(IntKi)                                :: i=0

      INTEGER(IntKi)                               :: ErrStat2      ! Error status of the operation
      CHARACTER(ErrMsgLen)                         :: ErrMsg2       ! Error message if ErrStat2 /= ErrID_None

      ErrStat = ErrID_None
      ErrMsg  = ""



      ! deallocate data associated with file output
      CALL MDIO_CloseOutput ( p, m, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         !IF (ErrStat >= AbortErrLev) RETURN


      ! deallocate FAST data structures
      CALL MD_DestroyInput(u, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyParam(p, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyContState(x, ErrStat2, ErrMsg2)  ! <--- getting access violation
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyDiscState(xd, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyConstrState(z, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyOtherState(other,ErrStat2,ErrMsg2) ! <--- getting access violation
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyOutput(y, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyMisc(m, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
         
         
         !TODO: any need to specifically deallocate things like m%xTemp%states in the above? <<<<

 !     IF ( ErrStat==ErrID_None) THEN
 !        CALL WrScr('MoorDyn closed without errors')
 !     ELSE
 !        CALL WrScr('MoorDyn closed with errors')
 !     END IF


   CONTAINS

      SUBROUTINE CheckError(ErrId, Msg)
        ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev

         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


         IF ( ErrID /= ErrID_None ) THEN

            IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine  ! keep existing error message if there is one
            ErrMsg = TRIM(ErrMsg)//' MD_End:'//TRIM(Msg)      ! add current error message
            ErrStat = MAX(ErrStat, ErrID)

            CALL WrScr( ErrMsg )  ! do this always or only if warning level?

         END IF

      END SUBROUTINE CheckError


   END SUBROUTINE MD_End                                                                         !   -------+
   !----------------------------------------------------------------------------------------==================


!!==========   MD_CheckError   =======     <---------------------------------------------------------------+
! SUBROUTINE MD_CheckError(InMsg,OutMsg)
!   ! Passed arguments
!!   CHARACTER(*), INTENT(IN   ) :: InMsg       ! The input string
!   CHARACTER(*), INTENT(INOUT) :: OutMsg      ! The error message (ErrMsg)!
!
 !  OutMsg = InMsg
 !  RETURN
 !END SUBROUTINE MD_CheckError                                                                  !   -------+
 !----------------------------------------------------------------------------------------==================


   ! RK2 integrater (part of what was in TimeStep)
   !--------------------------------------------------------------
   SUBROUTINE MD_RK2 ( t, dtM, u_interp, u, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg )
   
      REAL(DbKi)                     , INTENT(INOUT)      :: t          ! intial time (s) for this integration step
      REAL(DbKi)                     , INTENT(IN   )      :: dtM        ! single time step  size (s) for this integration step
      TYPE( MD_InputType )           , INTENT(INOUT)      :: u_interp   ! interpolated instantaneous input values to be calculated for each mooring time step
      TYPE( MD_InputType )           , INTENT(INOUT)      :: u(:)       ! INTENT(IN   )
      REAL(DbKi)                     , INTENT(IN   )      :: t_array(:)  ! times corresponding to elements of u(:)?
      TYPE( MD_ParameterType )       , INTENT(IN   )      :: p          ! INTENT(IN   )
      TYPE( MD_ContinuousStateType ) , INTENT(INOUT)      :: x
      TYPE( MD_DiscreteStateType )   , INTENT(IN   )      :: xd         ! INTENT(IN   )
      TYPE( MD_ConstraintStateType ) , INTENT(IN   )      :: z          ! INTENT(IN   )
      TYPE( MD_OtherStateType )      , INTENT(IN   )      :: other      ! INTENT(INOUT)
      TYPE(MD_MiscVarType)           , INTENT(INOUT)      :: m          ! INTENT(INOUT)
      INTEGER(IntKi)                 , INTENT(  OUT)      :: ErrStat
      CHARACTER(*)                   , INTENT(  OUT)      :: ErrMsg


      INTEGER(IntKi)                                      :: I          ! counter
      INTEGER(IntKi)                                      :: J          ! counter
         
   
      ! -------------------------------------------------------------------------------
      !       RK2 integrator written here, now calling CalcContStateDeriv
      !--------------------------------------------------------------------------------

      ! step 1

      CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t          , ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t)
   
      CALL MD_CalcContStateDeriv( t, u_interp, p, x, xd, z, other, m, m%xdTemp, ErrStat, ErrMsg )
      DO J = 1, m%Nx
         m%xTemp%states(J) = x%states(J) + 0.5*dtM*m%xdTemp%states(J)                                           !x1 = x0 + dt*f0/2.0;
      END DO

      ! step 2

      CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t + 0.5_DbKi*dtM, ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t+0.5*dtM)
         
      CALL MD_CalcContStateDeriv( (t + 0.5_DbKi*dtM), u_interp, p, m%xTemp, xd, z, other, m, m%xdTemp, ErrStat, ErrMsg )       !called with updated states x2 and time = t + dt/2.0
      DO J = 1, m%Nx
         x%states(J) = x%states(J) + dtM*m%xdTemp%states(J)
      END DO

      t = t + dtM  ! update time
      
      !TODO error check? <<<<

   END SUBROUTINE MD_RK2
   !--------------------------------------------------------------


   !----------------------------------------------------------------------------------------================
   ! this would do a full (coupling) time step and is no longer used
   SUBROUTINE TimeStep ( t, dtStep, u, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg )
   
      REAL(DbKi)                     , INTENT(INOUT)      :: t
      REAL(DbKi)                     , INTENT(IN   )      :: dtStep     ! how long to advance the time for
      TYPE( MD_InputType )           , INTENT(INOUT)      :: u(:)       ! INTENT(IN   )
      REAL(DbKi)                     , INTENT(IN   )      :: t_array(:)  ! times corresponding to elements of u(:)?
      TYPE( MD_ParameterType )       , INTENT(IN   )      :: p          ! INTENT(IN   )
      TYPE( MD_ContinuousStateType ) , INTENT(INOUT)      :: x
      TYPE( MD_DiscreteStateType )   , INTENT(IN   )      :: xd         ! INTENT(IN   )
      TYPE( MD_ConstraintStateType ) , INTENT(IN   )      :: z          ! INTENT(IN   )
      TYPE( MD_OtherStateType )      , INTENT(IN   )      :: other      ! INTENT(INOUT)
      TYPE(MD_MiscVarType)           , INTENT(INOUT)      :: m          ! INTENT(INOUT)
      INTEGER(IntKi)                 , INTENT(  OUT)      :: ErrStat
      CHARACTER(*)                   , INTENT(  OUT)      :: ErrMsg


      TYPE(MD_ContinuousStateType)                        :: dxdt       ! time derivatives of continuous states (initialized in CalcContStateDeriv)
      TYPE(MD_ContinuousStateType)                        :: x2         ! temporary copy of continuous states used in RK2 calculations
      INTEGER(IntKi)                                      :: NdtM       ! the number of time steps to make with the mooring model
      Real(DbKi)                                          :: dtM        ! the actual time step size to use
      INTEGER(IntKi)                                      :: Nx         ! size of states vector
      INTEGER(IntKi)                                      :: I          ! counter
      INTEGER(IntKi)                                      :: J          ! counter
      TYPE(MD_InputType)                                  :: u_interp   ! interpolated instantaneous input values to be calculated for each mooring time step

  !    Real(DbKi)                                          :: tDbKi   ! double version because that's what MD_Input_ExtrapInterp needs.
      
      
      ! allocate space for x2
      CALL MD_CopyContState( x, x2, 0, ErrStat, ErrMsg)
         
      ! create space for arrays/meshes in u_interp   ... is it efficient to do this every time step???
      CALL MD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg)
         

      Nx = size(x%states)   ! <<<< should this be the m%Nx parameter instead?


      ! round dt to integer number of time steps
      NdtM = ceiling(dtStep/p%dtM0)                  ! get number of mooring time steps to do based on desired time step size
      dtM = dtStep/float(NdtM)                       ! adjust desired time step to satisfy dt with an integer number of time steps


      !loop through line integration time steps
      DO I = 1, NdtM                                 ! for (double ts=t; ts<=t+ICdt-dts; ts+=dts)
      
      
   !      tDbKi = t        ! get DbKi version of current time (why does ExtrapInterp except different time type than UpdateStates?)
         
      
         ! -------------------------------------------------------------------------------
         !       RK2 integrator written here, now calling CalcContStateDeriv
         !--------------------------------------------------------------------------------

         ! step 1

         CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t          , ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t)
      
         CALL MD_CalcContStateDeriv( t, u_interp, p, x, xd, z, other, m, dxdt, ErrStat, ErrMsg )
         DO J = 1, Nx
            x2%states(J) = x%states(J) + 0.5*dtM*dxdt%states(J)                                           !x1 = x0 + dt*f0/2.0;
         END DO

         ! step 2
   
         CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t + 0.5_DbKi*dtM, ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t+0.5*dtM)
            
         CALL MD_CalcContStateDeriv( (t + 0.5_DbKi*dtM), u_interp, p, x2, xd, z, other, m, dxdt, ErrStat, ErrMsg )       !called with updated states x2 and time = t + dt/2.0
         DO J = 1, Nx
            x%states(J) = x%states(J) + dtM*dxdt%states(J)
         END DO

         t = t + dtM  ! update time

         !----------------------------------------------------------------------------------

   ! >>> below should no longer be necessary thanks to using ExtrapInterp of u(:) within the mooring time stepping loop.. <<<
   !      ! update Fairlead positions by integrating velocity and last position (do this AFTER the processing of the time step rather than before)
   !      DO J = 1, p%nCpldCons
   !         DO K = 1, 3
   !          m%ConnectList(m%CpldConIs(J))%r(K) = m%ConnectList(m%CpldConIs(J))%r(K) + m%ConnectList(m%CpldConIs(J))%rd(K)*dtM
   !         END DO
   !      END DO
      
   
      END DO  ! I  time steps


      ! destroy dxdt and x2, and u_interp
      CALL MD_DestroyContState( dxdt, ErrStat, ErrMsg)
      CALL MD_DestroyContState( x2, ErrStat, ErrMsg)
      CALL MD_DestroyInput(u_interp, ErrStat, ErrMsg)
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error destroying dxdt or x2.'
      END IF
      !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_UpdateStates')

      
      ! check for NaNs - is this a good place/way to do it?
      DO J = 1, Nx
         IF (Is_NaN(REAL(x%states(J),DbKi))) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' NaN state detected.'
         END IF
      END DO
 

   END SUBROUTINE TimeStep
   !--------------------------------------------------------------



   !-----------------------------------------------------------------------
   ! >>>>>>>>>>>>>> rename/reorganize this subroutine >>>>>>>>>>>>>
   SUBROUTINE SetupLine (Line, LineProp, rhoW, ErrStat, ErrMsg)
      ! allocate arrays in line object

      TYPE(MD_Line), INTENT(INOUT)       :: Line          ! the single line object of interest
      TYPE(MD_LineProp), INTENT(INOUT)   :: LineProp      ! the single line property set for the line of interest
      REAL(DbKi),    INTENT(IN)          :: rhoW
      INTEGER,       INTENT(   INOUT )   :: ErrStat       ! returns a non-zero value when an error occurs
      CHARACTER(*),  INTENT(   INOUT )   :: ErrMsg        ! Error message if ErrStat /= ErrID_None

      INTEGER(4)                         :: J             ! Generic index
      INTEGER(4)                         :: K             ! Generic index
      INTEGER(IntKi)                     :: N

      N = Line%N  ! number of segments in this line (for code readability)

      ! -------------- save some section properties to the line object itself -----------------

      Line%d   = LineProp%d
      Line%rho = LineProp%w/(Pi/4.0 * Line%d * Line%d)
      
      Line%EA   = LineProp%EA
      
      Line%Can   = LineProp%Can
      Line%Cat   = LineProp%Cat
      Line%Cdn   = LineProp%Cdn
      Line%Cdt   = LineProp%Cdt      
      
      ! Specify specific internal damping coefficient (BA) for this line.
      ! Will be equal to inputted BA of LineType if input value is positive.
      ! If input value is negative, it is considered to be desired damping ratio (zeta)
      ! from which the line's BA can be calculated based on the segment natural frequency.
      IF (LineProp%BA < 0) THEN
         ! - we assume desired damping coefficient is zeta = -LineProp%BA
         ! - highest axial vibration mode of a segment is wn = sqrt(k/m) = 2N/UnstrLen*sqrt(EA/w)
         Line%BA = -LineProp%BA * Line%UnstrLen / Line%N * SQRT(LineProp%EA * LineProp%w)
      !  print *, 'Based on zeta, BA set to ', Line%BA
         
      !  print *, 'Negative BA input detected, treating as -zeta.  For zeta = ', -LineProp%BA, ', setting BA to ', Line%BA
         
      ELSE
         Line%BA = LineProp%BA
      !  temp = Line%N * Line%BA / Line%UnstrLen * SQRT(1.0/(LineProp%EA * LineProp%w))
      !  print *, 'BA set as input to ', Line%BA, '. Corresponding zeta is ', temp
      END IF
      
      !temp = 2*Line%N / Line%UnstrLen * sqrt( LineProp%EA / LineProp%w) / TwoPi
      !print *, 'Segment natural frequency is ', temp, ' Hz'
      
      

      ! allocate node positions and velocities (NOTE: these arrays start at ZERO)
      ALLOCATE ( Line%r(3, 0:N), Line%rd(3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating r and rd arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! allocate node tangent vectors
      ALLOCATE ( Line%q(3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating q array.'
         !CALL CleanUp()
         RETURN
      END IF

      ! allocate segment scalar quantities
      ALLOCATE ( Line%l(N), Line%ld(N), Line%lstr(N), Line%lstrd(N), Line%V(N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating segment scalar quantity arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! assign values for l, ld, and V
      DO J=1,N
         Line%l(J) = Line%UnstrLen/REAL(N, DbKi)
         Line%ld(J)= 0.0_DbKi
         Line%V(J) = Line%l(J)*0.25*Pi*LineProp%d*LineProp%d
      END DO
      
      ! allocate water related vectors
      ALLOCATE ( Line%U(3, 0:N), Line%Ud(3, 0:N), Line%zeta(0:N), Line%PDyn(0:N), STAT = ErrStat )
      ! set to zero initially (important of wave kinematics are not being used)
      Line%U    = 0.0_DbKi
      Line%Ud   = 0.0_DbKi
      Line%zeta = 0.0_DbKi
      Line%PDyn = 0.0_DbKi

      ! allocate segment tension and internal damping force vectors
      ALLOCATE ( Line%T(3, N), Line%Td(3, N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating T and Td arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! allocate node force vectors
      ALLOCATE ( Line%W(3, 0:N), Line%Dp(3, 0:N), Line%Dq(3, 0:N), Line%Ap(3, 0:N), &
         Line%Aq(3, 0:N), Line%B(3, 0:N), Line%Fnet(3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating node force arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! set gravity and bottom contact forces to zero initially (because the horizontal components should remain at zero)
      DO J = 0,N
         DO K = 1,3
            Line%W(K,J) = 0.0_DbKi
            Line%B(K,J) = 0.0_DbKi
         END DO
      END DO

      ! allocate mass and inverse mass matrices for each node (including ends)
      ALLOCATE ( Line%S(3, 3, 0:N), Line%M(3, 3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating T and Td arrays.'
         !CALL CleanUp()
         RETURN
      END IF
      
    
      
      ! need to add cleanup sub <<<


   END SUBROUTINE SetupLine
   !--------------------------------------------------------------





   !----------------------------------------------------------------------------------------=======
   SUBROUTINE Line_Initialize (Line, LineProp, rhoW, ErrStat, ErrMsg)
      ! calculate initial profile of the line using quasi-static model

      TYPE(MD_Line),     INTENT(INOUT)       :: Line        ! the single line object of interest
      TYPE(MD_LineProp), INTENT(INOUT)       :: LineProp    ! the single line property set for the line of interest
      REAL(DbKi),        INTENT(IN)          :: rhoW
      INTEGER,           INTENT(   INOUT )   :: ErrStat     ! returns a non-zero value when an error occurs
      CHARACTER(*),      INTENT(   INOUT )   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      REAL(DbKi)                             :: COSPhi      ! Cosine of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
      REAL(DbKi)                             :: SINPhi      ! Sine   of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
      REAL(DbKi)                             :: XF          ! Horizontal distance between anchor and fairlead of the current mooring line (meters)
      REAL(DbKi)                             :: ZF          ! Vertical   distance between anchor and fairlead of the current mooring line (meters)
      INTEGER(4)                             :: I           ! Generic index
      INTEGER(4)                             :: J           ! Generic index


      INTEGER(IntKi)                         :: ErrStat2      ! Error status of the operation
      CHARACTER(ErrMsgLen)                   :: ErrMsg2       ! Error message if ErrStat2 /= ErrID_None
      REAL(DbKi)                             :: WetWeight
      REAL(DbKi)                             :: SeabedCD = 0.0_DbKi
      REAL(DbKi)                             :: TenTol = 0.0001_DbKi
      REAL(DbKi), ALLOCATABLE                :: LSNodes(:)
      REAL(DbKi), ALLOCATABLE                :: LNodesX(:)
      REAL(DbKi), ALLOCATABLE                :: LNodesZ(:)
      INTEGER(IntKi)                         :: N


      N = Line%N ! for convenience

       ! try to calculate initial line profile using catenary routine (from FAST v.7)
       ! note: much of this function is adapted from the FAST source code

          ! Transform the fairlead location from the inertial frame coordinate system
          !   to the local coordinate system of the current line (this coordinate
          !   system lies at the current anchor, Z being vertical, and X directed from
          !   current anchor to the current fairlead).  Also, compute the orientation
          !   of this local coordinate system:

             XF         = SQRT( ( Line%r(1,N) - Line%r(1,0) )**2.0 + ( Line%r(2,N) - Line%r(2,0) )**2.0 )
             ZF         =         Line%r(3,N) - Line%r(3,0)

             IF ( XF == 0.0 )  THEN  ! .TRUE. if the current mooring line is exactly vertical; thus, the solution below is ill-conditioned because the orientation is undefined; so set it such that the tensions and nodal positions are only vertical
                COSPhi  = 0.0_DbKi
                SINPhi  = 0.0_DbKi
             ELSE                    ! The current mooring line must not be vertical; use simple trigonometry
                COSPhi  =       ( Line%r(1,N) - Line%r(1,0) )/XF
                SINPhi  =       ( Line%r(2,N) - Line%r(2,0) )/XF
             ENDIF

        WetWeight = LineProp%w - 0.25*Pi*LineProp%d*LineProp%d*rhoW

        !LineNodes = Line%N + 1  ! number of nodes in line for catenary model to worry about

        ! allocate temporary arrays for catenary routine
        ALLOCATE ( LSNodes(N+1), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Error allocating LSNodes array.'
          CALL CleanUp()
          RETURN
        END IF

        ALLOCATE ( LNodesX(N+1), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Error allocating LNodesX array.'
          CALL CleanUp()
          RETURN
        END IF

        ALLOCATE ( LNodesZ(N+1), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Error allocating LNodesZ array.'
          CALL CleanUp()
          RETURN
        END IF

        ! Assign node arc length locations
        LSNodes(1) = 0.0_DbKi
        DO I=2,N
          LSNodes(I) = LSNodes(I-1) + Line%l(I-1)  ! note: l index is because line segment indices start at 1
        END DO
        LSNodes(N+1) = Line%UnstrLen  ! ensure the last node length isn't longer than the line due to numerical error

          ! Solve the analytical, static equilibrium equations for a catenary (or
          !   taut) mooring line with seabed interaction in order to find the
          !   horizontal and vertical tensions at the fairlead in the local coordinate
          !   system of the current line:
          ! NOTE: The values for the horizontal and vertical tensions at the fairlead
          !       from the previous time step are used as the initial guess values at
          !       at this time step (because the LAnchHTe(:) and LAnchVTe(:) arrays
          !       are stored in a module and thus their values are saved from CALL to
          !       CALL).


             CALL Catenary ( XF           , ZF          , Line%UnstrLen, LineProp%EA  , &
                             WetWeight    , SeabedCD,    TenTol,     (N+1)     , &
                             LSNodes, LNodesX, LNodesZ , ErrStat2, ErrMsg2)

      IF (ErrStat2 == ErrID_None) THEN ! if it worked, use it
          ! Transform the positions of each node on the current line from the local
          !   coordinate system of the current line to the inertial frame coordinate
          !   system:

          DO J = 0,N ! Loop through all nodes per line where the line position and tension can be output
             Line%r(1,J) = Line%r(1,0) + LNodesX(J+1)*COSPhi
             Line%r(2,J) = Line%r(2,0) + LNodesX(J+1)*SINPhi
             Line%r(3,J) = Line%r(3,0) + LNodesZ(J+1)
          ENDDO              ! J - All nodes per line where the line position and tension can be output


      ELSE ! if there is a problem with the catenary approach, just stretch the nodes linearly between fairlead and anchor

          CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Line_Initialize')

!          print *, "Node positions: "

          DO J = 0,N ! Loop through all nodes per line where the line position and tension can be output
             Line%r(1,J) = Line%r(1,0) + (Line%r(1,N) - Line%r(1,0))*REAL(J, DbKi)/REAL(N, DbKi)
             Line%r(2,J) = Line%r(2,0) + (Line%r(2,N) - Line%r(2,0))*REAL(J, DbKi)/REAL(N, DbKi)
             Line%r(3,J) = Line%r(3,0) + (Line%r(3,N) - Line%r(3,0))*REAL(J, DbKi)/REAL(N, DbKi)
             
!             print*, Line%r(:,J)
          ENDDO
          
!          print*,"FYI line end A and B node coords are"
!          print*, Line%r(:,0)
!          print*, Line%r(:,N)
      ENDIF



      CALL CleanUp()  ! deallocate temporary arrays



   CONTAINS


      !-----------------------------------------------------------------------
      SUBROUTINE CleanUp()
           ! deallocate temporary arrays

           IF (ALLOCATED(LSNodes))  DEALLOCATE(LSNodes)
           IF (ALLOCATED(LNodesX))  DEALLOCATE(LNodesX)
           IF (ALLOCATED(LNodesZ))  DEALLOCATE(LNodesZ)

        END SUBROUTINE CleanUp
      !-----------------------------------------------------------------------


      !-----------------------------------------------------------------------
      SUBROUTINE Catenary ( XF_In, ZF_In, L_In  , EA_In, &
                            W_In , CB_In, Tol_In, N    , &
                            s_In , X_In , Z_In , ErrStat, ErrMsg    )

         ! This subroutine is copied from FAST v7 with minor modifications

         ! This routine solves the analytical, static equilibrium equations
         ! for a catenary (or taut) mooring line with seabed interaction.
         ! Stretching of the line is accounted for, but bending stiffness
         ! is not.  Given the mooring line properties and the fairlead
         ! position relative to the anchor, this routine finds the line
         ! configuration and tensions.  Since the analytical solution
         ! involves two nonlinear equations (XF and  ZF) in two unknowns
         ! (HF and VF), a Newton-Raphson iteration scheme is implemented in
         ! order to solve for the solution.  The values of HF and VF that
         ! are passed into this routine are used as the initial guess in
         ! the iteration.  The Newton-Raphson iteration is only accurate in
         ! double precision, so all of the input/output arguments are
         ! converteds to/from double precision from/to default precision.

         ! >>>> TO DO: streamline this function, if it's still to be used at all <<<<

         !     USE                             Precision


         IMPLICIT                        NONE


            ! Passed Variables:

         INTEGER(4), INTENT(IN   )    :: N                                               ! Number of nodes where the line position and tension can be output (-)

         REAL(DbKi), INTENT(IN   )    :: CB_In                                           ! Coefficient of seabed static friction drag (a negative value indicates no seabed) (-)
         REAL(DbKi), INTENT(IN   )    :: EA_In                                           ! Extensional stiffness of line (N)
     !    REAL(DbKi), INTENT(  OUT)    :: HA_In                                           ! Effective horizontal tension in line at the anchor   (N)
     !    REAL(DbKi), INTENT(INOUT)    :: HF_In                                           ! Effective horizontal tension in line at the fairlead (N)
         REAL(DbKi), INTENT(IN   )    :: L_In                                            ! Unstretched length of line (meters)
         REAL(DbKi), INTENT(IN   )    :: s_In     (N)                                    ! Unstretched arc distance along line from anchor to each node where the line position and tension can be output (meters)
     !    REAL(DbKi), INTENT(  OUT)    :: Te_In    (N)                                    ! Effective line tensions at each node (N)
         REAL(DbKi), INTENT(IN   )    :: Tol_In                                          ! Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
     !    REAL(DbKi), INTENT(  OUT)    :: VA_In                                           ! Effective vertical   tension in line at the anchor   (N)
     !    REAL(DbKi), INTENT(INOUT)    :: VF_In                                           ! Effective vertical   tension in line at the fairlead (N)
         REAL(DbKi), INTENT(IN   )    :: W_In                                            ! Weight of line in fluid per unit length (N/m)
         REAL(DbKi), INTENT(  OUT)    :: X_In     (N)                                    ! Horizontal locations of each line node relative to the anchor (meters)
         REAL(DbKi), INTENT(IN   )    :: XF_In                                           ! Horizontal distance between anchor and fairlead (meters)
         REAL(DbKi), INTENT(  OUT)    :: Z_In     (N)                                    ! Vertical   locations of each line node relative to the anchor (meters)
         REAL(DbKi), INTENT(IN   )    :: ZF_In                                           ! Vertical   distance between anchor and fairlead (meters)
             INTEGER,                      INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs
             CHARACTER(*),                 INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None


            ! Local Variables:

         REAL(DbKi)                   :: CB                                              ! Coefficient of seabed static friction (a negative value indicates no seabed) (-)
         REAL(DbKi)                   :: CBOvrEA                                         ! = CB/EA
         REAL(DbKi)                   :: DET                                             ! Determinant of the Jacobian matrix (m^2/N^2)
         REAL(DbKi)                   :: dHF                                             ! Increment in HF predicted by Newton-Raphson (N)
         REAL(DbKi)                   :: dVF                                             ! Increment in VF predicted by Newton-Raphson (N)
         REAL(DbKi)                   :: dXFdHF                                          ! Partial derivative of the calculated horizontal distance with respect to the horizontal fairlead tension (m/N): dXF(HF,VF)/dHF
         REAL(DbKi)                   :: dXFdVF                                          ! Partial derivative of the calculated horizontal distance with respect to the vertical   fairlead tension (m/N): dXF(HF,VF)/dVF
         REAL(DbKi)                   :: dZFdHF                                          ! Partial derivative of the calculated vertical   distance with respect to the horizontal fairlead tension (m/N): dZF(HF,VF)/dHF
         REAL(DbKi)                   :: dZFdVF                                          ! Partial derivative of the calculated vertical   distance with respect to the vertical   fairlead tension (m/N): dZF(HF,VF)/dVF
         REAL(DbKi)                   :: EA                                              ! Extensional stiffness of line (N)
         REAL(DbKi)                   :: EXF                                             ! Error function between calculated and known horizontal distance (meters): XF(HF,VF) - XF
         REAL(DbKi)                   :: EZF                                             ! Error function between calculated and known vertical   distance (meters): ZF(HF,VF) - ZF
         REAL(DbKi)                   :: HA                                              ! Effective horizontal tension in line at the anchor   (N)
         REAL(DbKi)                   :: HF                                              ! Effective horizontal tension in line at the fairlead (N)
         REAL(DbKi)                   :: HFOvrW                                          ! = HF/W
         REAL(DbKi)                   :: HFOvrWEA                                        ! = HF/WEA
         REAL(DbKi)                   :: L                                               ! Unstretched length of line (meters)
         REAL(DbKi)                   :: Lamda0                                          ! Catenary parameter used to generate the initial guesses of the horizontal and vertical tensions at the fairlead for the Newton-Raphson iteration (-)
         REAL(DbKi)                   :: LMax                                            ! Maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead) (meters)
         REAL(DbKi)                   :: LMinVFOvrW                                      ! = L - VF/W
         REAL(DbKi)                   :: LOvrEA                                          ! = L/EA
         REAL(DbKi)                   :: s        (N)                                    ! Unstretched arc distance along line from anchor to each node where the line position and tension can be output (meters)
         REAL(DbKi)                   :: sOvrEA                                          ! = s(I)/EA
         REAL(DbKi)                   :: SQRT1VFOvrHF2                                   ! = SQRT( 1.0_DbKi + VFOvrHF2      )
         REAL(DbKi)                   :: SQRT1VFMinWLOvrHF2                              ! = SQRT( 1.0_DbKi + VFMinWLOvrHF2 )
         REAL(DbKi)                   :: SQRT1VFMinWLsOvrHF2                             ! = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF )
         REAL(DbKi)                   :: Te       (N)                                    ! Effective line tensions at each node (N)
         REAL(DbKi)                   :: Tol                                             ! Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
         REAL(DbKi)                   :: VA                                              ! Effective vertical   tension in line at the anchor   (N)
         REAL(DbKi)                   :: VF                                              ! Effective vertical   tension in line at the fairlead (N)
         REAL(DbKi)                   :: VFMinWL                                         ! = VF - WL
         REAL(DbKi)                   :: VFMinWLOvrHF                                    ! = VFMinWL/HF
         REAL(DbKi)                   :: VFMinWLOvrHF2                                   ! = VFMinWLOvrHF*VFMinWLOvrHF
         REAL(DbKi)                   :: VFMinWLs                                        ! = VFMinWL + Ws
         REAL(DbKi)                   :: VFMinWLsOvrHF                                   ! = VFMinWLs/HF
         REAL(DbKi)                   :: VFOvrHF                                         ! = VF/HF
         REAL(DbKi)                   :: VFOvrHF2                                        ! = VFOvrHF*VFOvrHF
         REAL(DbKi)                   :: VFOvrWEA                                        ! = VF/WEA
         REAL(DbKi)                   :: W                                               ! Weight of line in fluid per unit length (N/m)
         REAL(DbKi)                   :: WEA                                             ! = W*EA
         REAL(DbKi)                   :: WL                                              ! Total weight of line in fluid (N): W*L
         REAL(DbKi)                   :: Ws                                              ! = W*s(I)
         REAL(DbKi)                   :: X        (N)                                    ! Horizontal locations of each line node relative to the anchor (meters)
         REAL(DbKi)                   :: XF                                              ! Horizontal distance between anchor and fairlead (meters)
         REAL(DbKi)                   :: XF2                                             ! = XF*XF
         REAL(DbKi)                   :: Z        (N)                                    ! Vertical   locations of each line node relative to the anchor (meters)
         REAL(DbKi)                   :: ZF                                              ! Vertical   distance between anchor and fairlead (meters)
         REAL(DbKi)                   :: ZF2                                             ! = ZF*ZF

         INTEGER(4)                   :: I                                               ! Index for counting iterations or looping through line nodes (-)
         INTEGER(4)                   :: MaxIter                                         ! Maximum number of Newton-Raphson iterations possible before giving up (-)

         LOGICAL                      :: FirstIter                                       ! Flag to determine whether or not this is the first time through the Newton-Raphson interation (flag)


         ErrStat = ERrId_None


            ! The Newton-Raphson iteration is only accurate in double precision, so
            !   convert the input arguments into double precision:

         CB     = REAL( CB_In    , DbKi )
         EA     = REAL( EA_In    , DbKi )
         HF = 0.0_DbKi  !    = REAL( HF_In    , DbKi )
         L      = REAL( L_In     , DbKi )
         s  (:) = REAL( s_In  (:), DbKi )
         Tol    = REAL( Tol_In   , DbKi )
         VF = 0.0_DbKi   ! keeping this for some error catching functionality? (at first glance)  ! VF     = REAL( VF_In    , DbKi )
         W      = REAL( W_In     , DbKi )
         XF     = REAL( XF_In    , DbKi )
         ZF     = REAL( ZF_In    , DbKi )


         
      !  HF and VF cannot be initialized to zero when a  portion of the line rests on the seabed and the anchor tension is nonzero
         
      ! Generate the initial guess values for the horizontal and vertical tensions
      !   at the fairlead in the Newton-Raphson iteration for the catenary mooring
      !   line solution.  Use starting values documented in: Peyrot, Alain H. and
      !   Goulois, A. M., "Analysis Of Cable Structures," Computers & Structures,
      !   Vol. 10, 1979, pp. 805-813:
         XF2     = XF*XF
         ZF2     = ZF*ZF

         IF     ( XF           == 0.0_DbKi    )  THEN ! .TRUE. if the current mooring line is exactly vertical
            Lamda0 = 1.0D+06
         ELSEIF ( L <= SQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
            Lamda0 = 0.2_DbKi
         ELSE                                    ! The current mooring line must be slack and not vertical
            Lamda0 = SQRT( 3.0_DbKi*( ( L**2 - ZF2 )/XF2 - 1.0_DbKi ) )
         ENDIF

         HF = ABS( 0.5_DbKi*W*  XF/     Lamda0      )
         VF =      0.5_DbKi*W*( ZF/TANH(Lamda0) + L )         
                                    

            ! Abort when there is no solution or when the only possible solution is
            !   illogical:

         IF (    Tol <= EPSILON(TOL) )  THEN   ! .TRUE. when the convergence tolerance is specified incorrectly
           ErrStat = ErrID_Warn
           ErrMsg = ' Convergence tolerance must be greater than zero in routine Catenary().'
           return
         ELSEIF ( XF <  0.0_DbKi )  THEN   ! .TRUE. only when the local coordinate system is not computed correctly
           ErrStat = ErrID_Warn
           ErrMsg =  ' The horizontal distance between an anchor and its'// &
                         ' fairlead must not be less than zero in routine Catenary().'
           return

         ELSEIF ( ZF <  0.0_DbKi )  THEN   ! .TRUE. if the fairlead has passed below its anchor
           ErrStat = ErrID_Warn
           ErrMsg =  ' A fairlead has passed below its anchor.'
           return

         ELSEIF ( L  <= 0.0_DbKi )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly
           ErrStat = ErrID_Warn
           ErrMsg =  ' Unstretched length of line must be greater than zero in routine Catenary().'
           return

         ELSEIF ( EA <= 0.0_DbKi )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly
           ErrStat = ErrID_Warn
           ErrMsg =  ' Extensional stiffness of line must be greater than zero in routine Catenary().'
           return

         ELSEIF ( W  == 0.0_DbKi )  THEN   ! .TRUE. when the weight of the line in fluid is zero so that catenary solution is ill-conditioned
           ErrStat = ErrID_Warn
           ErrMsg = ' The weight of the line in fluid must not be zero. '// &
                         ' Routine Catenary() cannot solve quasi-static mooring line solution.'
           return


         ELSEIF ( W  >  0.0_DbKi )  THEN   ! .TRUE. when the line will sink in fluid

            LMax      = XF - EA/W + SQRT( (EA/W)*(EA/W) + 2.0_DbKi*ZF*EA/W )  ! Compute the maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)

            IF ( ( L  >=  LMax   ) .AND. ( CB >= 0.0_DbKi ) )  then  ! .TRUE. if the line is as long or longer than its maximum possible value with seabed interaction
               ErrStat = ErrID_Warn
               ErrMsg =  ' Unstretched mooring line length too large. '// &
                            ' Routine Catenary() cannot solve quasi-static mooring line solution.'
               return
            END IF

         ENDIF


            ! Initialize some commonly used terms that don't depend on the iteration:

         WL      =          W  *L
         WEA     =          W  *EA
         LOvrEA  =          L  /EA
         CBOvrEA =          CB /EA
         MaxIter = INT(1.0_DbKi/Tol)   ! Smaller tolerances may take more iterations, so choose a maximum inversely proportional to the tolerance



            ! To avoid an ill-conditioned situation, ensure that the initial guess for
            !   HF is not less than or equal to zero.  Similarly, avoid the problems
            !   associated with having exactly vertical (so that HF is zero) or exactly
            !   horizontal (so that VF is zero) lines by setting the minimum values
            !   equal to the tolerance.  This prevents us from needing to implement
            !   the known limiting solutions for vertical or horizontal lines (and thus
            !   complicating this routine):

         HF = MAX( HF, Tol )
         XF = MAX( XF, Tol )
         ZF = MAX( ZF, TOl )



            ! Solve the analytical, static equilibrium equations for a catenary (or
            !   taut) mooring line with seabed interaction:

            ! Begin Newton-Raphson iteration:

         I         = 1        ! Initialize iteration counter
         FirstIter = .TRUE.   ! Initialize iteration flag

         DO


            ! Initialize some commonly used terms that depend on HF and VF:

            VFMinWL            = VF - WL
            LMinVFOvrW         = L  - VF/W
            HFOvrW             =      HF/W
            HFOvrWEA           =      HF/WEA
            VFOvrWEA           =      VF/WEA
            VFOvrHF            =      VF/HF
            VFMinWLOvrHF       = VFMinWL/HF
            VFOvrHF2           = VFOvrHF     *VFOvrHF
            VFMinWLOvrHF2      = VFMinWLOvrHF*VFMinWLOvrHF
            SQRT1VFOvrHF2      = SQRT( 1.0_DbKi + VFOvrHF2      )
            SQRT1VFMinWLOvrHF2 = SQRT( 1.0_DbKi + VFMinWLOvrHF2 )


            ! Compute the error functions (to be zeroed) and the Jacobian matrix
            !   (these depend on the anticipated configuration of the mooring line):

            IF ( ( CB <  0.0_DbKi ) .OR. ( W  <  0.0_DbKi ) .OR. ( VFMinWL >  0.0_DbKi ) )  THEN   ! .TRUE. when no portion of the line      rests on the seabed

               EXF    = (   LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                       &
                          - LOG( VFMinWLOvrHF +               SQRT1VFMinWLOvrHF2 )                                         )*HFOvrW &
                      + LOvrEA*  HF                         - XF
               EZF    = (                                     SQRT1VFOvrHF2                                              &
                          -                                   SQRT1VFMinWLOvrHF2                                           )*HFOvrW &
                      + LOvrEA*( VF - 0.5_DbKi*WL )         - ZF

               dXFdHF = (   LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                       &
                          - LOG( VFMinWLOvrHF +               SQRT1VFMinWLOvrHF2 )                                         )/     W &
                      - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      ) &
                          -    ( VFMinWLOvrHF + VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2 )/( VFMinWLOvrHF + SQRT1VFMinWLOvrHF2 )   )/     W &
                      + LOvrEA
               dXFdVF = (      ( 1.0_DbKi     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      ) &
                          -    ( 1.0_DbKi     + VFMinWLOvrHF /SQRT1VFMinWLOvrHF2 )/( VFMinWLOvrHF + SQRT1VFMinWLOvrHF2 )   )/     W
               dZFdHF = (                                     SQRT1VFOvrHF2                                              &
                          -                                   SQRT1VFMinWLOvrHF2                                           )/     W &
                      - (                       VFOvrHF2     /SQRT1VFOvrHF2                                              &
                          -                     VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2                                           )/     W
               dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                              &
                          -                     VFMinWLOvrHF /SQRT1VFMinWLOvrHF2                                           )/     W &
                      + LOvrEA


            ELSEIF (                                           -CB*VFMinWL <  HF         )  THEN   ! .TRUE. when a  portion of the line      rests on the seabed and the anchor tension is nonzero

               EXF    =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          *HFOvrW &
                      - 0.5_DbKi*CBOvrEA*W*  LMinVFOvrW*LMinVFOvrW                                                                  &
                      + LOvrEA*  HF           + LMinVFOvrW  - XF
               EZF    = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi   )*HFOvrW &
                      + 0.5_DbKi*VF*VFOvrWEA                - ZF

               dXFdHF =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          /     W &
                      - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2        ) )/     W &
                      + LOvrEA
               dXFdVF = (      ( 1.0_DbKi     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2        ) )/     W &
                      + CBOvrEA*LMinVFOvrW - 1.0_DbKi/W
               dZFdHF = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi &
                          -                     VFOvrHF2     /SQRT1VFOvrHF2                                                )/     W
               dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                                )/     W &
                      + VFOvrWEA


            ELSE                                                ! 0.0_DbKi <  HF  <= -CB*VFMinWL   !             A  portion of the line must rest  on the seabed and the anchor tension is    zero

               EXF    =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          *HFOvrW &
                      - 0.5_DbKi*CBOvrEA*W*( LMinVFOvrW*LMinVFOvrW - ( LMinVFOvrW - HFOvrW/CB )*( LMinVFOvrW - HFOvrW/CB ) )        &
                      + LOvrEA*  HF           + LMinVFOvrW  - XF
               EZF    = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi   )*HFOvrW &
                      + 0.5_DbKi*VF*VFOvrWEA                - ZF

               dXFdHF =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          /     W &
                      - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      )   )/     W &
                      + LOvrEA - ( LMinVFOvrW - HFOvrW/CB )/EA
               dXFdVF = (      ( 1.0_DbKi     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      )   )/     W &
                      + HFOvrWEA           - 1.0_DbKi/W
               dZFdHF = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi &
                          -                     VFOvrHF2     /SQRT1VFOvrHF2                                                )/     W
               dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                                )/     W &
                      + VFOvrWEA


            ENDIF


            ! Compute the determinant of the Jacobian matrix and the incremental
            !   tensions predicted by Newton-Raphson:
            
            
            DET = dXFdHF*dZFdVF - dXFdVF*dZFdHF
            
            if ( EqualRealNos( DET, 0.0_DbKi ) ) then               
!bjj: there is a serious problem with the debugger here when DET = 0
                ErrStat = ErrID_Warn
                ErrMsg =  ' Iteration not convergent (DET is 0). '// &
                          ' Routine Catenary() cannot solve quasi-static mooring line solution.'
                return
            endif

               
            dHF = ( -dZFdVF*EXF + dXFdVF*EZF )/DET    ! This is the incremental change in horizontal tension at the fairlead as predicted by Newton-Raphson
            dVF = (  dZFdHF*EXF - dXFdHF*EZF )/DET    ! This is the incremental change in vertical   tension at the fairlead as predicted by Newton-Raphson

            dHF = dHF*( 1.0_DbKi - Tol*I )            ! Reduce dHF by factor (between 1 at I = 1 and 0 at I = MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)
            dVF = dVF*( 1.0_DbKi - Tol*I )            ! Reduce dHF by factor (between 1 at I = 1 and 0 at I = MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)

            dHF = MAX( dHF, ( Tol - 1.0_DbKi )*HF )   ! To avoid an ill-conditioned situation, make sure HF does not go less than or equal to zero by having a lower limit of Tol*HF [NOTE: the value of dHF = ( Tol - 1.0_DbKi )*HF comes from: HF = HF + dHF = Tol*HF when dHF = ( Tol - 1.0_DbKi )*HF]

            ! Check if we have converged on a solution, or restart the iteration, or
            !   Abort if we cannot find a solution:

            IF ( ( ABS(dHF) <= ABS(Tol*HF) ) .AND. ( ABS(dVF) <= ABS(Tol*VF) ) )  THEN ! .TRUE. if we have converged; stop iterating! [The converge tolerance, Tol, is a fraction of tension]

               EXIT


            ELSEIF ( ( I == MaxIter )        .AND. (       FirstIter         ) )  THEN ! .TRUE. if we've iterated MaxIter-times for the first time;

            ! Perhaps we failed to converge because our initial guess was too far off.
            !   (This could happen, for example, while linearizing a model via large
            !   pertubations in the DOFs.)  Instead, use starting values documented in:
            !   Peyrot, Alain H. and Goulois, A. M., "Analysis Of Cable Structures,"
            !   Computers & Structures, Vol. 10, 1979, pp. 805-813:
            ! NOTE: We don't need to check if the current mooring line is exactly
            !       vertical (i.e., we don't need to check if XF == 0.0), because XF is
            !       limited by the tolerance above.

               XF2 = XF*XF
               ZF2 = ZF*ZF

               IF ( L <= SQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
                  Lamda0 = 0.2_DbKi
               ELSE                                ! The current mooring line must be slack and not vertical
                  Lamda0 = SQRT( 3.0_DbKi*( ( L*L - ZF2 )/XF2 - 1.0_DbKi ) )
               ENDIF

               HF  = MAX( ABS( 0.5_DbKi*W*  XF/     Lamda0      ), Tol )   ! As above, set the lower limit of the guess value of HF to the tolerance
               VF  =           0.5_DbKi*W*( ZF/TANH(Lamda0) + L )


            ! Restart Newton-Raphson iteration:

               I         = 0
               FirstIter = .FALSE.
               dHF       = 0.0_DbKi
               dVF       = 0.0_DbKi


           ELSEIF ( ( I == MaxIter )        .AND. ( .NOT. FirstIter         ) )  THEN ! .TRUE. if we've iterated as much as we can take without finding a solution; Abort
             ErrStat = ErrID_Warn
             ErrMsg =  ' Iteration not convergent. '// &
                       ' Routine Catenary() cannot solve quasi-static mooring line solution.'
             RETURN


           ENDIF


            ! Increment fairlead tensions and iteration counter so we can try again:

            HF = HF + dHF
            VF = VF + dVF

            I  = I  + 1


         ENDDO



            ! We have found a solution for the tensions at the fairlead!

            ! Now compute the tensions at the anchor and the line position and tension
            !   at each node (again, these depend on the configuration of the mooring
            !   line):

         IF ( ( CB <  0.0_DbKi ) .OR. ( W  <  0.0_DbKi ) .OR. ( VFMinWL >  0.0_DbKi ) )  THEN   ! .TRUE. when no portion of the line      rests on the seabed

            ! Anchor tensions:

            HA = HF
            VA = VFMinWL


            ! Line position and tension at each node:

            DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

               IF ( ( s(I) <  0.0_DbKi ) .OR. ( s(I) >  L ) )  THEN
                 ErrStat = ErrID_Warn
                 ErrMsg = ' All line nodes must be located between the anchor ' &
                                 //'and fairlead (inclusive) in routine Catenary().'
                 RETURN
               END IF

               Ws                  = W       *s(I)                                  ! Initialize
               VFMinWLs            = VFMinWL + Ws                                   ! some commonly
               VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
               sOvrEA              = s(I)    /EA                                    ! that depend
               SQRT1VFMinWLsOvrHF2 = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

               X (I)    = (   LOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 ) &
                            - LOG( VFMinWLOvrHF  + SQRT1VFMinWLOvrHF2  )   )*HFOvrW                     &
                        + sOvrEA*  HF
               Z (I)    = (                        SQRT1VFMinWLsOvrHF2   &
                            -                      SQRT1VFMinWLOvrHF2      )*HFOvrW                     &
                        + sOvrEA*(         VFMinWL + 0.5_DbKi*Ws    )
               Te(I)    = SQRT( HF*HF +    VFMinWLs*VFMinWLs )

            ENDDO       ! I - All nodes where the line position and tension are to be computed


         ELSEIF (                                           -CB*VFMinWL <  HF         )  THEN   ! .TRUE. when a  portion of the line      rests on the seabed and the anchor tension is nonzero

            ! Anchor tensions:

            HA = HF + CB*VFMinWL
            VA = 0.0_DbKi


            ! Line position and tension at each node:

            DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

               IF ( ( s(I) <  0.0_DbKi ) .OR. ( s(I) >  L ) )  THEN
                 ErrStat = ErrID_Warn
                 ErrMsg =  ' All line nodes must be located between the anchor ' &
                                 //'and fairlead (inclusive) in routine Catenary().'
                 RETURN
               END IF

               Ws                  = W       *s(I)                                  ! Initialize
               VFMinWLs            = VFMinWL + Ws                                   ! some commonly
               VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
               sOvrEA              = s(I)    /EA                                    ! that depend
               SQRT1VFMinWLsOvrHF2 = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

               IF (     s(I) <= LMinVFOvrW             )  THEN ! .TRUE. if this node rests on the seabed and the tension is nonzero

                  X (I) = s(I)                                                                          &
                        + sOvrEA*( HF + CB*VFMinWL + 0.5_DbKi*Ws*CB )
                  Z (I) = 0.0_DbKi
                  Te(I) =       HF    + CB*VFMinWLs

               ELSE           ! LMinVFOvrW < s <= L            !           This node must be above the seabed

                  X (I) =     LOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 )    *HFOvrW                     &
                        + sOvrEA*  HF + LMinVFOvrW                    - 0.5_DbKi*CB*VFMinWL*VFMinWL/WEA
                  Z (I) = ( - 1.0_DbKi           + SQRT1VFMinWLsOvrHF2     )*HFOvrW                     &
                        + sOvrEA*(         VFMinWL + 0.5_DbKi*Ws    ) + 0.5_DbKi*   VFMinWL*VFMinWL/WEA
                  Te(I) = SQRT( HF*HF +    VFMinWLs*VFMinWLs )

               ENDIF

            ENDDO       ! I - All nodes where the line position and tension are to be computed


         ELSE                                                ! 0.0_DbKi <  HF  <= -CB*VFMinWL   !             A  portion of the line must rest  on the seabed and the anchor tension is    zero

            ! Anchor tensions:

            HA = 0.0_DbKi
            VA = 0.0_DbKi


            ! Line position and tension at each node:

            DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

               IF ( ( s(I) <  0.0_DbKi ) .OR. ( s(I) >  L ) )  THEN
                  ErrStat = ErrID_Warn
                   ErrMsg =  ' All line nodes must be located between the anchor ' &
                                 //'and fairlead (inclusive) in routine Catenary().'
                   RETURN
               END IF

               Ws                  = W       *s(I)                                  ! Initialize
               VFMinWLs            = VFMinWL + Ws                                   ! some commonly
               VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
               sOvrEA              = s(I)    /EA                                    ! that depend
               SQRT1VFMinWLsOvrHF2 = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

               IF (     s(I) <= LMinVFOvrW - HFOvrW/CB )  THEN ! .TRUE. if this node rests on the seabed and the tension is    zero

                  X (I) = s(I)
                  Z (I) = 0.0_DbKi
                  Te(I) = 0.0_DbKi

               ELSEIF ( s(I) <= LMinVFOvrW             )  THEN ! .TRUE. if this node rests on the seabed and the tension is nonzero

                  X (I) = s(I)                     - ( LMinVFOvrW - 0.5_DbKi*HFOvrW/CB )*HF/EA          &
                        + sOvrEA*( HF + CB*VFMinWL + 0.5_DbKi*Ws*CB ) + 0.5_DbKi*CB*VFMinWL*VFMinWL/WEA
                  Z (I) = 0.0_DbKi
                  Te(I) =       HF    + CB*VFMinWLs

               ELSE           ! LMinVFOvrW < s <= L            !           This node must be above the seabed

                  X (I) =     LOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 )    *HFOvrW                     &
                        + sOvrEA*  HF + LMinVFOvrW - ( LMinVFOvrW - 0.5_DbKi*HFOvrW/CB )*HF/EA
                  Z (I) = ( - 1.0_DbKi           + SQRT1VFMinWLsOvrHF2     )*HFOvrW                     &
                        + sOvrEA*(         VFMinWL + 0.5_DbKi*Ws    ) + 0.5_DbKi*   VFMinWL*VFMinWL/WEA
                  Te(I) = SQRT( HF*HF +    VFMinWLs*VFMinWLs )

               ENDIF

            ENDDO       ! I - All nodes where the line position and tension are to be computed


         ENDIF



            ! The Newton-Raphson iteration is only accurate in double precision, so
            !   convert the output arguments back into the default precision for real
            !   numbers:

         !HA_In    = REAL( HA   , DbKi )  !mth: for this I only care about returning node positions
         !HF_In    = REAL( HF   , DbKi )
         !Te_In(:) = REAL( Te(:), DbKi )
         !VA_In    = REAL( VA   , DbKi )
         !VF_In    = REAL( VF   , DbKi )
         X_In (:) = REAL( X (:), DbKi )
         Z_In (:) = REAL( Z (:), DbKi )

      END SUBROUTINE Catenary
      !-----------------------------------------------------------------------


   END SUBROUTINE Line_Initialize
   !--------------------------------------------------------------

   
   !--------------------------------------------------------------
   SUBROUTINE Line_SetState(Line, X, t)

      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! the current Line object
      Real(DbKi),       INTENT(IN   )  :: X(:)           ! state vector section for this line
      Real(DbKi),       INTENT(IN   )  :: t              ! instantaneous time

      INTEGER(IntKi)                   :: i              ! index of segments or nodes along line
      INTEGER(IntKi)                   :: J              ! index
   

      ! store current time
      Line%time = t
      
      ! set interior node positions and velocities based on state vector
      DO I=1,Line%N-1
         DO J=1,3
         
            Line%r( J,I) = X( 3*Line%N-3 + 3*I-3 + J)  ! get positions
            Line%rd(J,I) = X(              3*I-3 + J)  ! get velocities
            
         END DO
      END DO
         
   END SUBROUTINE Line_SetState
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Line_GetStateDeriv(Line, Xd, p)  !, FairFtot, FairMtot, AnchFtot, AnchMtot)

      TYPE(MD_Line), INTENT(INOUT)     :: Line          ! the current Line object
      Real(DbKi),    INTENT(INOUT)     :: Xd(:)         ! state derivative vector section for this line
      TYPE(MD_ParameterType), INTENT(IN )   :: p       ! Parameters
      
      
   !   Real(DbKi), INTENT( IN )      :: X(:)           ! state vector, provided
   !   Real(DbKi), INTENT( INOUT )   :: Xd(:)          ! derivative of state vector, returned ! cahnged to INOUT
   !   Real(DbKi), INTENT (IN)       :: t              ! instantaneous time
   !   TYPE(MD_Line), INTENT (INOUT) :: Line           ! label for the current line, for convenience
   !   TYPE(MD_LineProp), INTENT(IN) :: LineProp       ! the single line property set for the line of interest
  !    Real(DbKi), INTENT(INOUT)     :: FairFtot(:)    ! total force on Connect top of line is attached to
  !    Real(DbKi), INTENT(INOUT)     :: FairMtot(:,:)  ! total mass of Connect top of line is attached to
  !    Real(DbKi), INTENT(INOUT)     :: AnchFtot(:)    ! total force on Connect bottom of line is attached to
  !    Real(DbKi), INTENT(INOUT)     :: AnchMtot(:,:)  ! total mass of Connect bottom of line is attached to


      INTEGER(IntKi)                   :: i              ! index of segments or nodes along line
      INTEGER(IntKi)                   :: J              ! index
      INTEGER(IntKi)                   :: K              ! index
      INTEGER(IntKi)                   :: N              ! number of segments in line
      Real(DbKi)                       :: d              ! line diameter
      Real(DbKi)                       :: rho            ! line material density [kg/m^3]
      Real(DbKi)                       :: Sum1           ! for summing squares
      Real(DbKi)                       :: dummyLength    ! 
      Real(DbKi)                       :: m_i            ! node mass
      Real(DbKi)                       :: v_i            ! node submerged volume
      Real(DbKi)                       :: Vi(3)          ! relative water velocity at a given node
      Real(DbKi)                       :: Vp(3)          ! transverse relative water velocity component at a given node
      Real(DbKi)                       :: Vq(3)          ! tangential relative water velocity component at a given node
      Real(DbKi)                       :: SumSqVp        !
      Real(DbKi)                       :: SumSqVq        !
      Real(DbKi)                       :: MagVp          !
      Real(DbKi)                       :: MagVq          !


      N = Line%N                      ! for convenience
      d = Line%d    
      rho = Line%rho

      ! note that end node kinematics should have already been set by attached objects

!      ! set end node positions and velocities from connect objects' states
!      DO J = 1, 3
!         Line%r( J,N) = m%ConnectList(Line%FairConnect)%r(J)
!         Line%r( J,0) = m%ConnectList(Line%AnchConnect)%r(J)
!         Line%rd(J,N) = m%ConnectList(Line%FairConnect)%rd(J)
!         Line%rd(J,0) = m%ConnectList(Line%AnchConnect)%rd(J)
!      END DO



      ! calculate instantaneous (stretched) segment lengths and rates << should add catch here for if lstr is ever zero
      DO I = 1, N
         Sum1 = 0.0_DbKi
         DO J = 1, 3
            Sum1 = Sum1 + (Line%r(J,I) - Line%r(J,I-1)) * (Line%r(J,I) - Line%r(J,I-1))
         END DO
         Line%lstr(I) = sqrt(Sum1)                                  ! stretched segment length

         Sum1 = 0.0_DbKi
         DO J = 1, 3
            Sum1 = Sum1 + (Line%r(J,I) - Line%r(J,I-1))*(Line%rd(J,I) - Line%rd(J,I-1))
         END DO
         Line%lstrd(I) = Sum1/Line%lstr(I)                          ! segment stretched length rate of change

 !       Line%V(I) = Pi/4.0 * d*d*Line%l(I)                        !volume attributed to segment
      END DO

      !calculate unit tangent vectors (q) for each node (including ends) note: I think these are pointing toward 0 rather than N!
      CALL UnitVector(Line%r(:,0), Line%r(:,1), Line%q(:,0), dummyLength) ! compute unit vector q
      DO I = 1, N-1
        CALL UnitVector(Line%r(:,I-1), Line%r(:,I+1), Line%q(:,I), dummyLength) ! compute unit vector q ... using adjacent two nodes!
      END DO
      CALL UnitVector(Line%r(:,N-1), Line%r(:,N), Line%q(:,N), dummyLength)     ! compute unit vector q


      ! --------------------------------- apply wave kinematics ------------------------------------

      IF (p%WaterKin > 0)  THEN ! wave kinematics interpolated from global grid in Waves object
         DO i=0,N
            CALL getWaveKin(p, Line%r(1,i), Line%r(2,i), Line%r(3,i), Line%time, Line%U(:,i), Line%Ud(:,i), Line%zeta(i), Line%PDyn(i))
         END DO
      END IF


      ! --------------- calculate mass (including added mass) matrix for each node -----------------
      DO I = 0, N
         IF (I==0) THEN
            m_i = Pi/8.0 *d*d*Line%l(1)*rho
            v_i = 0.5 *Line%V(1)
         ELSE IF (I==N) THEN
            m_i = pi/8.0 *d*d*Line%l(N)*rho;
            v_i = 0.5*Line%V(N)
         ELSE
            m_i = pi/8.0 * d*d*rho*(Line%l(I) + Line%l(I+1))
            v_i = 0.5 *(Line%V(I) + Line%V(I+1))
         END IF

         DO J=1,3
            DO K=1,3
               IF (J==K) THEN
                  Line%M(K,J,I) = m_i + p%rhoW*v_i*( Line%Can*(1 - Line%q(J,I)*Line%q(K,I)) + Line%Cat*Line%q(J,I)*Line%q(K,I) )
               ELSE
                  Line%M(K,J,I) = p%rhoW*v_i*( Line%Can*(-Line%q(J,I)*Line%q(K,I)) + Line%Cat*Line%q(J,I)*Line%q(K,I) )
               END IF
            END DO
         END DO

         CALL Inverse3by3(Line%S(:,:,I), Line%M(:,:,I))             ! invert mass matrix
      END DO


      ! ------------------  CALCULATE FORCES ON EACH NODE ----------------------------

      ! loop through the segments
      DO I = 1, N

         ! line tension, inherently including possibility of dynamic length changes in l term
         IF (Line%lstr(I)/Line%l(I) > 1.0) THEN
            DO J = 1, 3
               Line%T(J,I) = Line%EA *( 1.0/Line%l(I) - 1.0/Line%lstr(I) ) * (Line%r(J,I)-Line%r(J,I-1))
            END DO
         ELSE
            DO J = 1, 3
               Line%T(J,I) = 0.0_DbKi                              ! cable can't "push"
            END DO
         END if

         ! line internal damping force based on line-specific BA value, including possibility of dynamic length changes in l and ld terms
         DO J = 1, 3
            !Line%Td(J,I) = Line%BA* ( Line%lstrd(I) / Line%l(I) ) * (Line%r(J,I)-Line%r(J,I-1)) / Line%lstr(I)  ! note new form of damping coefficient, BA rather than Cint
            Line%Td(J,I) = Line%BA* ( Line%lstrd(I) -  Line%lstr(I)*Line%ld(I)/Line%l(I) )/Line%l(I)  * (Line%r(J,I)-Line%r(J,I-1)) / Line%lstr(I)
         END DO
      END DO



      ! loop through the nodes
      DO I = 0, N

         !submerged weight (including buoyancy)
         IF (I==0) THEN
            Line%W(3,I) = Pi/8.0*d*d* Line%l(1)*(rho - p%rhoW) *(-p%g)   ! assuming g is positive
         ELSE IF (i==N)  THEN
            Line%W(3,I) = pi/8.0*d*d* Line%l(N)*(rho - p%rhoW) *(-p%g)
         ELSE
            Line%W(3,I) = pi/8.0*d*d* (Line%l(I)*(rho - p%rhoW) + Line%l(I+1)*(rho - p%rhoW) )*(-p%g)  ! left in this form for future free surface handling
         END IF

         !relative flow velocities
         DO J = 1, 3
            Vi(J) = 0.0 - Line%rd(J,I)                               ! relative flow velocity over node -- this is where wave velicites would be added
         END DO

         ! decomponse relative flow into components
         SumSqVp = 0.0_DbKi                                         ! start sums of squares at zero
         SumSqVq = 0.0_DbKi
         DO J = 1, 3
            Vq(J) = DOT_PRODUCT( Vi , Line%q(:,I) ) * Line%q(J,I);   ! tangential relative flow component
            Vp(J) = Vi(J) - Vq(J)                                    ! transverse relative flow component
            SumSqVq = SumSqVq + Vq(J)*Vq(J)
            SumSqVp = SumSqVp + Vp(J)*Vp(J)
         END DO
         MagVp = sqrt(SumSqVp)                                      ! get magnitudes of flow components
         MagVq = sqrt(SumSqVq)

         ! transverse and tangenential drag
         IF (I==0) THEN
            Line%Dp(:,I) = 0.25*p%rhoW*Line%Cdn*    d*Line%l(1) * MagVp * Vp
            Line%Dq(:,I) = 0.25*p%rhoW*Line%Cdt* Pi*d*Line%l(1) * MagVq * Vq
         ELSE IF (I==N)  THEN
            Line%Dp(:,I) = 0.25*p%rhoW*Line%Cdn*    d*Line%l(N) * MagVp * Vp
            Line%Dq(:,I) = 0.25*p%rhoW*Line%Cdt* Pi*d*Line%l(N) * MagVq * Vq
         ELSE
            Line%Dp(:,I) = 0.25*p%rhoW*Line%Cdn*    d*(Line%l(I) + Line%l(I+1)) * MagVp * vp
            Line%Dq(:,I) = 0.25*p%rhoW*Line%Cdt* Pi*d*(Line%l(I) + Line%l(I+1)) * MagVq * vq
         END IF

         ! F-K force from fluid acceleration not implemented yet

         ! bottom contact (stiffness and damping, vertical-only for now)  - updated Nov 24 for general case where anchor and fairlead ends may deal with bottom contact forces

         IF (Line%r(3,I) < -p%WtrDpth) THEN
            IF (I==0) THEN
               Line%B(3,I) = ( (-p%WtrDpth - Line%r(3,I))*p%kBot - Line%rd(3,I)*p%cBot) * 0.5*d*(            Line%l(I+1) ) 
            ELSE IF (I==N) THEN
               Line%B(3,I) = ( (-p%WtrDpth - Line%r(3,I))*p%kBot - Line%rd(3,I)*p%cBot) * 0.5*d*(Line%l(I)               ) 
            ELSE
               Line%B(3,I) = ( (-p%WtrDpth - Line%r(3,I))*p%kBot - Line%rd(3,I)*p%cBot) * 0.5*d*(Line%l(I) + Line%l(I+1) ) 



            END IF
         ELSE
            Line%B(3,I) = 0.0_DbKi
         END IF

         ! total forces
         IF (I==0)  THEN
            Line%Fnet(:,I) = Line%T(:,1)                 + Line%Td(:,1)                  + Line%W(:,I) + Line%Dp(:,I) + Line%Dq(:,I) + Line%B(:,I)
         ELSE IF (I==N)  THEN
            Line%Fnet(:,I) =                -Line%T(:,N)                  - Line%Td(:,N) + Line%W(:,I) + Line%Dp(:,I) + Line%Dq(:,I) + Line%B(:,I)
         ELSE
            Line%Fnet(:,I) = Line%T(:,I+1) - Line%T(:,I) + Line%Td(:,I+1) - Line%Td(:,I) + Line%W(:,I) + Line%Dp(:,I) + Line%Dq(:,I) + Line%B(:,I)
         END IF

      END DO  ! I  - done looping through nodes

 !     print *, "line ", Line%IdNum, " has N=", N
 !     print *, " and rd shape", shape(Line%rd)
 !     print *, " and Xd shape", shape(Xd)
      
      ! loop through internal nodes and update their states  <<< should/could convert to matrix operations instead of all these loops
      DO I=1, N-1
         DO J=1,3

            ! calculate RHS constant (premultiplying force vector by inverse of mass matrix  ... i.e. rhs = S*Forces)
            Sum1 = 0.0_DbKi                               ! reset temporary accumulator
            DO K = 1, 3
              Sum1 = Sum1 + Line%S(K,J,I) * Line%Fnet(K,I)   ! matrix-vector multiplication [S i]{Forces i}  << double check indices
            END DO ! K
            
 !           print *, "writing Xd indices ", 3*N-3 + 3*I-3 + J, " and " , 3*I-3 + J
 !
 !           print*, Line%rd(J,I)
  
            ! update states
            Xd(3*N-3 + 3*I-3 + J) = Line%rd(J,I);       ! dxdt = V  (velocities)
            Xd(        3*I-3 + J) = Sum1                ! dVdt = RHS * A  (accelerations)

         END DO ! J
      END DO  ! I


      ! check for NaNs
      DO J = 1, 6*(N-1)
         IF (Is_NaN(REAL(Xd(J),DbKi))) THEN
            print *, "NaN detected at time ", Line%time, " in Line ", Line%IdNum, " state derivatives:"
            print *, Xd
            
            
            
            print *, "m_i  p%rhoW   v_i Line%Can  Line%Cat"
            print *, m_i 
            print *, p%rhoW
            print *, v_i
            print *, Line%Can
            print *, Line%Cat
            
            print *, "Line%q"
            print *, Line%q
            
            print *, "Line%r"
            print *, Line%r
            
            
            print *, "Here is the mass matrix set"
            print *, Line%M
            
            print *, "Here is the inverted mass matrix set"
            print *, Line%S
            
            print *, "Here is the net force set"
            print *, Line%Fnet
            
            
            
            EXIT
         END IF
      END DO


   !   ! add force and mass of end nodes to the Connects they correspond to <<<<<<<<<<<< do this from Connection instead now!
   !   DO J = 1,3
   !      FairFtot(J) = FairFtot(J) + Line%F(J,N)
   !      AnchFtot(J) = AnchFtot(J) + Line%F(J,0)
   !      DO K = 1,3
   !         FairMtot(K,J) = FairMtot(K,J) + Line%M(K,J,N)
   !         AnchMtot(K,J) = AnchMtot(K,J) + Line%M(K,J,0)
   !      END DO
   !   END DO

   END SUBROUTINE Line_GetStateDeriv
   !=====================================================================


   !--------------------------------------------------------------
   SUBROUTINE Line_SetEndKinematics(Line, r_in, rd_in, t, topOfLine)

      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! the current Line object
      Real(DbKi),       INTENT(IN   )  :: r_in( 3)       ! state vector section for this line
      Real(DbKi),       INTENT(IN   )  :: rd_in(3)       ! state vector section for this line
      Real(DbKi),       INTENT(IN   )  :: t              ! instantaneous time
      INTEGER(IntKi),   INTENT(IN   )  :: topOfLine      ! 0 for end A (Node 0), 1 for end B (node N)

      Integer(IntKi)                   :: I,J      
      INTEGER(IntKi)                   :: inode
      
      IF (topOfLine==1) THEN
         inode = Line%N      
         Line%endTypeB = 0   ! set as ball rather than rigid connection (unless changed later by SetEndOrientation)
      ELSE
         inode = 0
         Line%endTypeA = 0   ! set as ball rather than rigid connection (unless changed later by SetEndOrientation)
      END IF 
      
      !Line%r( :,inode) = r_in
      !Line%rd(:,inode) = rd_in
      
      DO J = 1,3
         Line%r( :,inode) = r_in
         Line%rd(:,inode) = rd_in
      END DO
      
   !   print *, "SetEndKinematics of line ", Line%idNum, " top?:", topOfLine
   !   print *, r_in
   !   print *, Line%r( :,inode), "  - confirming, node ", inode
   !   print *, rd_in
  
      Line%time = t
   
   END SUBROUTINE Line_SetEndKinematics
   !--------------------------------------------------------------
   

   ! get force, moment, and mass of line at line end node
   !--------------------------------------------------------------
   SUBROUTINE Line_GetEndStuff(Line, Fnet_out, Moment_out, M_out, topOfLine)
   
      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! label for the current line, for convenience
      REAL(DbKi),       INTENT(  OUT)  :: Fnet_out(3)    ! net force on end node
      REAL(DbKi),       INTENT(  OUT)  :: Moment_out(3)  ! moment on end node (future capability)
      REAL(DbKi),       INTENT(  OUT)  :: M_out(3,3)     ! mass matrix of end node
      INTEGER(IntKi),   INTENT(IN   )  :: topOfLine      ! 0 for end A (Node 0), 1 for end B (node N)
      
      Integer(IntKi)                   :: I,J
      INTEGER(IntKi)                   :: inode
      
      IF (topOfLine==1) THEN           ! end B of line
         Fnet_out   = Line%Fnet(:, Line%N)
         Moment_out = Line%endMomentB
         M_out      = Line%M(:,:, Line%N)
      ELSE                             ! end A of line
         Fnet_out   = Line%Fnet(:, 0)
         Moment_out = Line%endMomentA
         M_out      = Line%M(:,:, 0)
      END IF

   END SUBROUTINE Line_GetEndStuff
   !--------------------------------------------------------------


   ! set end kinematics of a line that's attached to a Rod, and this includes rotational information
   !--------------------------------------------------------------
   SUBROUTINE Line_GetEndSegmentInfo(Line, qEnd, EIout, dlOut, topOfLine)

      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! label for the current line, for convenience
      REAL(DbKi),       INTENT(  OUT)  :: qEnd(3)
      REAL(DbKi),       INTENT(  OUT)  :: EIout
      REAL(DbKi),       INTENT(  OUT)  :: dlOut
      INTEGER(IntKi),   INTENT(IN   )  :: topOfLine      ! 0 for end A (Node 0), 1 for end B (node N)
      
      EIout = Line%EI;
      
      if (topOfLine==1) then
         CALL UnitVector(Line%r(:,Line%N-1), Line%r(:,Line%N), qEnd, dlOut)  ! unit vector of last line segment
      else
         CALL UnitVector(Line%r(:,0       ), Line%r(:,1     ), qEnd, dlOut)  ! unit vector of first line segment
      END IF
      
   END SUBROUTINE Line_GetEndSegmentInfo
   !--------------------------------------------------------------


   ! set end node unit vector of a line (this is called by an attached to a Rod, only applicable for bending stiffness)
   !--------------------------------------------------------------
   SUBROUTINE Line_SetEndOrientation(Line, qin, topOfLine, rodEndB)

      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! label for the current line, for convenience
      REAL(DbKi),       INTENT(IN   )  :: qin(3)         ! the rod's axis unit vector
      INTEGER(IntKi),   INTENT(IN   )  :: topOfLine      ! 0 for end A (Node 0), 1 for end B (node N)
      INTEGER(IntKi),   INTENT(IN   )  :: rodEndB        ! =0 means the line is attached to Rod end A, =1 means attached to Rod end B (implication for unit vector sign)
	
      if (topOfLine==1) then
      
         Line%endTypeB = 1                  ! indicate attached to Rod (at every time step, just in case line get detached)
         
         if (rodEndB==1) then
            Line%q(:,Line%N) = -qin   ! -----line----->[B<==ROD==A]
         else
            Line%q(:,Line%N) =  qin   ! -----line----->[A==ROD==>B]
         end if
      else
         
         Line%endTypeA = 1                  ! indicate attached to Rod (at every time step, just in case line get detached)                 ! indicate attached to Rod
         
         if (rodEndB==1) then
            Line%q(:,0     ) =  qin   ! [A==ROD==>B]-----line----->
         else
            Line%q(:,0     ) = -qin   ! [B<==ROD==A]-----line----->
         end if
      end if

   END SUBROUTINE Line_SetEndOrientation
   !--------------------------------------------------------------


!--------------------------------------------------------------
!            Connection-Specific Subroutines
!--------------------------------------------------------------





   !--------------------------------------------------------------
   SUBROUTINE Connect_Initialize(Connect, states, LineList)

      Type(MD_Connect), INTENT(INOUT)  :: Connect        ! the Connection object
      Real(DbKi),       INTENT(INOUT)  :: states(6)      ! state vector section for this Connection
      Type(MD_Line),    INTENT(INOUT)  :: LineList(:)    ! array of all the Line objects

      INTEGER(IntKi)                   :: l


      if (Connect%typeNum == 0) then  ! error check
      	
         ! pass kinematics to any attached lines so they have initial positions at this initialization stage
         DO l=1,Connect%nAttached
            print *, "Connect ",  Connect%IdNum, " setting end kinematics of line ", Connect%attached(l), " to ", Connect%r
            CALL Line_SetEndKinematics(LineList(Connect%attached(l)), Connect%r, Connect%rd, 0.0_DbKi, Connect%Top(l))
         END DO


         ! assign initial node kinematics to state vector
         states(4:6) = Connect%r
         states(1:3) = Connect%rd
         
         
         print *, "Initialized Connection ", Connect%IdNum
      
      else 
         print *,"   Error: wrong connection type given to Connect_Initialize for number ", Connect%idNum
      end if
      
   END SUBROUTINE Connect_Initialize
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   SUBROUTINE Connect_SetKinematics(Connect, r_in, rd_in, a_in, t, LineList)

      Type(MD_Connect), INTENT(INOUT)  :: Connect        ! the Connection object
      Real(DbKi),       INTENT(IN   )  :: r_in( 3)       ! position
      Real(DbKi),       INTENT(IN   )  :: rd_in(3)       ! velocity
      Real(DbKi),       INTENT(IN   )  :: a_in(3)        ! acceleration (only used for coupled connects)
      Real(DbKi),       INTENT(IN   )  :: t              ! instantaneous time
      Type(MD_Line),    INTENT(INOUT)  :: LineList(:)    ! array of all the Line objects

      INTEGER(IntKi)                   :: l

      ! store current time
      Connect%time = t

      
    !  if (Connect%typeNum==0) THEN ! anchor ( <<< to be changed/expanded) ... in MoorDyn F also used for coupled connections
                        
         ! set position and velocity
         Connect%r  = r_in
         Connect%rd = rd_in
         Connect%a = a_in
                 
         ! pass latest kinematics to any attached lines
         DO l=1,Connect%nAttached
            CALL Line_SetEndKinematics(LineList(Connect%attached(l)), Connect%r, Connect%rd, t, Connect%Top(l))
         END DO
      
     ! else
     !    
     !    PRINT*,"Error: setKinematics called for wrong Connection type. Connection ", Connect%IdNum, " type ", Connect%typeNum
         
   !  END IF
      
         
   END SUBROUTINE Connect_SetKinematics
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Connect_SetState(Connect, X, t, LineList)

      Type(MD_Connect),      INTENT(INOUT)  :: Connect        ! the Connection object
      Real(DbKi),            INTENT(IN   )  :: X(:)           ! state vector section for this line
      Real(DbKi),            INTENT(IN   )  :: t              ! instantaneous time
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:)    ! array of all the Line objects

      INTEGER(IntKi)                        :: l              ! index of segments or nodes along line
      INTEGER(IntKi)                        :: J              ! index
   

      ! store current time
      Connect%time = t
      
      ! from state values, get r and rdot values
      DO J=1,3
         Connect%r( J) = X(3 + J)  ! get positions
         Connect%rd(J) = X(    J)  ! get velocities
      END DO
           
     ! pass latest kinematics to any attached lines
      DO l=1,Connect%nAttached
         CALL Line_SetEndKinematics(LineList(Connect%attached(l)), Connect%r, Connect%rd, t, Connect%Top(l))
      END DO
      
   END SUBROUTINE Connect_SetState
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Connect_GetStateDeriv(Connect, Xd, LineList, p)

      Type(MD_Connect),      INTENT(INOUT)  :: Connect          ! the Connection object
      Real(DbKi),            INTENT(INOUT)  :: Xd(:)            ! state derivative vector section for this line
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:)      ! array of all the Line objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p                ! Parameters
      
      !TYPE(MD_MiscVarType), INTENT(INOUT)  :: m       ! misc/optimization variables

      !INTEGER(IntKi)             :: l         ! index of attached lines
      INTEGER(IntKi)                        :: J                ! index
      INTEGER(IntKi)                        :: K                ! index      
      Real(DbKi)                            :: Sum1             ! for adding things
      
      Real(DbKi)                            :: S(3,3)           ! inverse mass matrix


      CALL Connect_DoRHS(Connect, LineList, p)
                  
!      // solve for accelerations in [M]{a}={f} using LU decomposition
!      double M_tot[9];                     // serialize total mass matrix for easy processing
!      for (int I=0; I<3; I++) for (int J=0; J<3; J++) M_tot[3*I+J]=M[I][J];
!      double LU[9];                        // serialized matrix that will hold LU matrices combined
!      Crout(3, M_tot, LU);                  // perform LU decomposition on mass matrix
!      double acc[3];                        // acceleration vector to solve for
!      solveCrout(3, LU, Fnet, acc);     // solve for acceleration vector

      ! solve for accelerations in [M]{a}={f} using LU decomposition
!      CALL LUsolve(6, M_out, LU_temp, Fnet_out, y_temp, acc)
   
                  
      ! invert node mass matrix
      CALL Inverse3by3(S, Connect%M)

      ! accelerations 
      Connect%a = MATMUL(S, Connect%Fnet)

      ! fill in state derivatives
      Xd(4:6) = Connect%rd           ! dxdt = V    (velocities)
      Xd(1:3) = Connect%a            ! dVdt = RHS * A  (accelerations)
      

      ! check for NaNs
      DO J = 1, 6
         IF (Is_NaN(REAL(Xd(J),DbKi))) THEN
            print *, "NaN detected at time ", Connect%time, " in Connection ",Connect%IdNum, " state derivatives:"
            print *, Xd
            EXIT
         END IF
      END DO

   END SUBROUTINE Connect_GetStateDeriv
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Connect_DoRHS(Connect, LineList, p)

      Type(MD_Connect),      INTENT(INOUT)  :: Connect        ! the Connection object
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:)    ! array of all the Line objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p       ! Parameters
      
      !TYPE(MD_MiscVarType), INTENT(INOUT)  :: m       ! misc/optimization variables

      INTEGER(IntKi)             :: l         ! index of attached lines
      INTEGER(IntKi)             :: I         ! index
      INTEGER(IntKi)             :: J         ! index
      INTEGER(IntKi)             :: K         ! index

      Real(DbKi)                 :: Fnet_i(3) ! force from an attached line
      Real(DbKi)                 :: Moment_dummy(3) ! dummy vector to hold unused line end moments
      Real(DbKi)                 :: M_i(3,3)  ! mass from an attached line


      ! start with the Connection's own forces including buoyancy and weight, and its own mass
      Connect%Fnet(1) = Connect%conFX
      Connect%Fnet(2) = Connect%conFY
      Connect%Fnet(3) = Connect%conFZ + Connect%conV*p%rhoW*p%g - Connect%conM*p%g

      Connect%M    = 0.0_DbKi  ! clear (zero) the connect mass matrix
      
      DO J = 1,3
        Connect%M   (J,J) = Connect%conM  ! set the diagonals to the self-mass (to start with)
      END DO


   !   print *, "connection number", Connect%IdNum
   !   print *, "attached lines: ", Connect%attached
   !   print *, "size of line list" , size(m%LineList)

      ! loop through attached lines, adding force and mass contributions
      DO l=1,Connect%nAttached
         
      !   print *, "  l", l
      !   print *, Connect%attached(l)
      !   print *, m%LineList(Connect%attached(l))%Fnet
      !   
      !   
      !   print *, "  attached line ID", m%LineList(Connect%attached(l))%IdNum
         
         CALL Line_GetEndStuff(LineList(Connect%attached(l)), Fnet_i, Moment_dummy, M_i, Connect%Top(l))
         
         ! sum quantitites
         Connect%Fnet = Connect%Fnet + Fnet_i
         Connect%M    = Connect%M    + M_i
         
      END DO


      ! XXXWhen this sub is called, any self weight, buoyancy, or external forcing should have already been
      ! added by the calling subroutine.  The only thing left is any added mass or drag forces from the connection (e.g. float)
      ! itself, which will be added below.XXX


   !   IF (EqualRealNos(t, 0.0_DbKi)) THEN  ! this is old: with current IC gen approach, we skip the first call to the line objects, because they're set AFTER the call to the connects
   !
   !      DO J = 1,3
   !         Xd(3+J) = X(J)        ! velocities - these are unused in integration
   !         Xd(J) = 0.0_DbKi           ! accelerations - these are unused in integration
   !      END DO
   !   ELSE
   !      ! from state values, get r and rdot values
   !      DO J = 1,3
   !         Connect%r(J)  = X(3 + J)   ! get positions
   !         Connect%rd(J) = X(J)       ! get velocities
   !      END DO
   !   END IF
      

      ! add any added mass and drag forces from the Connect body itself
      DO J = 1,3
         Connect%Fnet(J)   = Connect%Fnet(J) - 0.5 * p%rhoW * Connect%rd(J) * abs(Connect%rd(J)) * Connect%conCdA;  ! add drag forces - corrected Nov 24
         Connect%M   (J,J) = Connect%M   (J,J) + Connect%conV*p%rhoW*Connect%conCa;                               ! add added mass

      END DO

   END SUBROUTINE Connect_DoRHS
   !=====================================================================


   ! calculate the force including inertial loads on connect that is coupled
   !--------------------------------------------------------------
   SUBROUTINE Connect_GetCoupledForce(Connect,  Fnet_out, LineList, p)
   
      Type(MD_Connect),      INTENT(INOUT)  :: Connect     ! the Connect object
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(3) ! force and moment vector about rRef
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:) ! array of all the Line objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters

      Real(DbKi)                            :: F_iner(3)   ! inertial force

      IF (Connect%typeNum == -1) then
         ! calculate forces and masses of connect
         CALL Connect_DoRHS(Connect, LineList, p)
	
         ! add inertial loads as appropriate
         F_iner = -MATMUL(Connect%M, Connect%a)    ! inertial loads
         Fnet_out = Connect%Fnet + F_iner          ! add inertial loads

      ELSE
         print *, "Connect_GetCoupledForce called for wrong (uncoupled) connect type!"
      END IF
      
   END SUBROUTINE Connect_GetCoupledForce


   ! calculate the force and mass contributions of the connect on the parent body (only for type 3 connects?)
   !--------------------------------------------------------------
   SUBROUTINE Connect_GetNetForceAndMass(Connect, rRef, Fnet_out, M_out, LineList, p)
   
      Type(MD_Connect),      INTENT(INOUT)  :: Connect     ! the Connect object
      Real(DbKi),            INTENT(IN   )  :: rRef(3)     ! global coordinates of reference point (i.e. the parent body)
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(6) ! force and moment vector about rRef
      Real(DbKi),            INTENT(  OUT)  :: M_out(6,6)  ! mass and inertia matrix about rRef
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:) ! array of all the Line objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters

      Real(DbKi)                            :: rRel(  3)   ! position of connection relative to the body reference point (global orientation frame)


      CALL Connect_DoRHS(Connect, LineList, p)
	
      rRel = Connect%r - rRef    ! vector from body reference point to node
			
      ! convert net force into 6dof force about body ref point
      CALL translateForce3to6DOF(rRel, Connect%Fnet, Fnet_out)
      
      ! convert mass matrix to 6by6 mass matrix about body ref point
      CALL translateMass3to6DOF(rRel, Connect%M, M_out)
	
   END SUBROUTINE Connect_GetNetForceAndMass
   
   
 
 
   ! this function handles assigning a line to a connection node
   !--------------------------------------------------------------
   SUBROUTINE Connect_AddLine(Connect, lineID, TopOfLine)

      Type(MD_Connect), INTENT (INOUT)   :: Connect        ! the Connection object
      Integer(IntKi),   INTENT( IN )     :: lineID
      Integer(IntKi),   INTENT( IN )     :: TopOfLine

      Print*, "L", lineID, "->C", Connect%IdNum
      
      IF (Connect%nAttached <10) THEN ! this is currently just a maximum imposed by a fixed array size.  could be improved.
         Connect%nAttached = Connect%nAttached + 1  ! add the line to the number connected
         Connect%Attached(Connect%nAttached) = lineID
         Connect%Top(Connect%nAttached) = TopOfLine  ! attached to line ... 1 = top/fairlead(end B), 0 = bottom/anchor(end A)
      ELSE
         Print*, "too many lines connected!"
      END IF

   END SUBROUTINE Connect_AddLine


   ! this function handles removing a line from a connection node
   !--------------------------------------------------------------
   SUBROUTINE Connect_RemoveLine(Connect, lineID, TopOfLine, rEnd, rdEnd)

      Type(MD_Connect), INTENT (INOUT)  :: Connect        ! the Connection object
      Integer(IntKi),   INTENT( IN )     :: lineID
      Integer(IntKi),   INTENT(  OUT)    :: TopOfLine
      REAL(DbKi),       INTENT(INOUT)    :: rEnd(3)
      REAL(DbKi),       INTENT(INOUT)    :: rdEnd(3)
      
      Integer(IntKi)    :: l,m,J
      
      DO l = 1,Connect%nAttached    ! look through attached lines
      
         IF (Connect%Attached(l) == lineID) THEN   ! if this is the line's entry in the attachment list
         
            TopOfLine = Connect%Top(l);                ! record which end of the line was attached
            
            DO m = l,Connect%nAttached-1 
            
               Connect%Attached(m) = Connect%Attached(m+1)  ! move subsequent line links forward one spot in the list to eliminate this line link
               Connect%Top(     m) =      Connect%Top(m+1) 
            
               Connect%nAttached = Connect%nAttached - 1                      ! reduce attached line counter by 1
            
               ! also pass back the kinematics at the end
               DO J = 1,3
                  rEnd( J) = Connect%r( J)
                  rdEnd(J) = Connect%rd(J)
               END DO
               
               print*, "Detached line ", lineID, " from Connection ", Connect%IdNum
               
               EXIT
            END DO
            
            IF (l == Connect%nAttached) THEN   ! detect if line not found
               print *, "Error: failed to find line to remove during removeLineFromConnect call to connection ", Connect%IdNum, ". Line ", lineID
            END IF
         
         END IF
         
      END DO
      
   END SUBROUTINE Connect_RemoveLine








!--------------------------------------------------------------
!            Rod-Specific Subroutines
!--------------------------------------------------------------



   !-----------------------------------------------------------------------
   SUBROUTINE Rod_Setup(Rod, RodProp, endCoords, rhoW, ErrStat, ErrMsg)
      ! calculate initial profile of the line using quasi-static model

      TYPE(MD_Rod),       INTENT(INOUT)  :: Rod          ! the single rod object of interest
      TYPE(MD_RodProp),   INTENT(INOUT)  :: RodProp      ! the single rod property set for the line of interest
      REAL(DbKi),    INTENT(IN)          :: endCoords(6)
      REAL(DbKi),    INTENT(IN)          :: rhoW
      INTEGER,       INTENT(   INOUT )   :: ErrStat       ! returns a non-zero value when an error occurs
      CHARACTER(*),  INTENT(   INOUT )   :: ErrMsg        ! Error message if ErrStat /= ErrID_None

      INTEGER(4)                         :: J             ! Generic index
      INTEGER(4)                         :: K             ! Generic index
      INTEGER(IntKi)                     :: N

      N = Rod%N  ! number of segments in this line (for code readability)

      ! -------------- save some section properties to the line object itself -----------------

      Rod%d   = RodProp%d
      Rod%rho = RodProp%w/(Pi/4.0 * Rod%d * Rod%d)
      
      Rod%Can   = RodProp%Can
      Rod%Cat   = RodProp%Cat
      Rod%Cdn   = RodProp%Cdn
      Rod%Cdt   = RodProp%Cdt      
      Rod%CaEnd = RodProp%CaEnd      
      Rod%CdEnd = RodProp%CdEnd      
      

      ! allocate node positions and velocities (NOTE: these arrays start at ZERO)
      ALLOCATE ( Rod%r(3, 0:N), Rod%rd(3, 0:N), STAT = ErrStat )   ! <<<<<< add error checks here
      IF ( ErrStat /= ErrID_None ) print *, "Alloc error 1!!!!!!!!!!!!!" 
     
      ! allocate segment scalar quantities
      ALLOCATE ( Rod%l(N), Rod%V(N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) print *, "Alloc error 2!!!!!!!!!!!!!" 

      ! allocate water related vectors
      ALLOCATE ( Rod%U(3, 0:N), Rod%Ud(3, 0:N), Rod%zeta(0:N), Rod%PDyn(0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) print *, "Alloc error 3!!!!!!!!!!!!!" 
      ! set to zero initially (important of wave kinematics are not being used)
      Rod%U    = 0.0_DbKi
      Rod%Ud   = 0.0_DbKi
      Rod%zeta = 0.0_DbKi
      Rod%PDyn = 0.0_DbKi

      ! allocate node force vectors
      ALLOCATE ( Rod%W(3, 0:N), Rod%Bo(3, 0:N), Rod%Dp(3, 0:N), Rod%Dq(3, 0:N), Rod%Ap(3, 0:N), &
         Rod%Aq(3, 0:N), Rod%Pd(3, 0:N), Rod%B(3, 0:N), Rod%Fnet(3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) print *, "Alloc error 4!!!!!!!!!!!!!" 
      
      ! allocate mass and inverse mass matrices for each node (including ends)
      ALLOCATE ( Rod%M(3, 3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) print *, "Alloc error 5!!!!!!!!!!!!!" 



      ! ------------------------- set some geometric properties and the starting kinematics -------------------------
		
		CALL UnitVector(endCoords(1:3), endCoords(4:6), Rod%q, Rod%UnstrLen)  ! get Rod axis direction vector and Rod length
		
		! set Rod positions if applicable
		if (Rod%typeNum==0) then               ! for an independent rod, set the position right off the bat
				
         Rod%r6(1:3) = endCoords(1:3)      ! (end A coordinates) 
         Rod%v6(1:3) = 0.0_DbKi            ! (end A velocity, unrotated axes) 
   
         Rod%r6(4:6) = Rod%q               ! (Rod direction unit vector)
         Rod%v6(4:6) = 0.0_DbKi            ! (rotational velocities about unrotated axes) 
			
		
		else if (abs(Rod%typeNum)==1) then    ! for a pinned rod, just set the orientation (position will be set later by parent object)
			
         Rod%r6(4:6) = Rod%q               ! (Rod direction unit vector)
         Rod%v6(4:6) = 0.0_DbKi            ! (rotational velocities about unrotated axes) 

		end if
		! otherwise (for a fixed rod) the positions will be set by the parent body or via coupling
		


      ! save mass for future calculations >>>> should calculate I_l and I_r here in future <<<<
      Rod%mass  = Rod%UnstrLen*RodProp%w


      ! assign values for l and V
      DO J=1,N
         Rod%l(J) = Rod%UnstrLen/REAL(N, DbKi)
         Rod%V(J) = Rod%l(J)*0.25*Pi*RodProp%d*RodProp%d
      END DO
      
      

      ! set gravity and bottom contact forces to zero initially (because the horizontal components should remain at zero)
      DO J = 0,N
         DO K = 1,3
            Rod%W(K,J) = 0.0_DbKi
            Rod%B(K,J) = 0.0_DbKi
         END DO
      END DO
      
      ! >>> why are the above assignments making l V W and B appear as "undefined pointer/array"s??? <<<
      
      print *, "Set up Rod ",Rod%IdNum, ", type ", Rod%typeNum

      ! need to add cleanup sub <<<

   END SUBROUTINE Rod_Setup
   !--------------------------------------------------------------




   ! Make output file for Rod and set end kinematics of any attached lines.
   ! For free Rods, fill in the initial states into the state vector.
   ! Notes: r6 and v6 must already be set.  
   !        ground- or body-pinned rods have already had setKinematics called to set first 3 elements of r6, v6.
   !--------------------------------------------------------------
   SUBROUTINE Rod_Initialize(Rod, states, LineList)

      TYPE(MD_Rod),          INTENT(INOUT)  :: Rod          ! the rod object 
      Real(DbKi),            INTENT(INOUT)  :: states(:)    ! state vector section for this line
      TYPE(MD_Line),         INTENT(INOUT)  :: LineList(:)  ! passing along all mooring objects

      INTEGER(IntKi)                        :: l           ! index of segments or nodes along line
      REAL(DbKi)                            :: rRef(3)     ! reference position of mesh node
      REAL(DbKi)                            :: OrMat(3,3)  ! DCM for body orientation based on r6_in
   
      print *, "initializing Rod ", Rod%idNum

      ! the r6 and v6 vectors should have already been set
      ! r and rd of ends have already been set by setup function or by parent object   <<<<< right? <<<<<


      ! Pass kinematics to any attached lines (this is just like what a Connection does, except for both ends)
      ! so that they have the correct initial positions at this initialization stage.
      
      if (Rod%typeNum >- 2)  CALL Rod_SetDependentKin(Rod, 0.0_DbKi, LineList)  ! don't call this for type -2 coupled Rods as it's already been called


      ! assign the resulting kinematics to its part of the state vector (only matters if it's an independent Rod)

      if (Rod%typeNum == 0) then               ! free Rod type
      
         states(1:6)   = 0.0_DbKi     ! zero velocities for initialization
         states(7:9)   = Rod%r(:,0)   ! end A position
         states(10:12) = Rod%q        ! rod direction unit vector
      
      else if (abs(Rod%typeNum) ==1 ) then           ! pinned rod type (coupled or attached to something previously via setPinKin)
      
         states(1:3)   = 0.0_DbKi     ! zero velocities for initialization
         states(4:6)   = Rod%q        ! rod direction unit vector
         
      end if
      
      ! note: this may also be called by a coupled rod (type = -1) in which case states will be empty
      
      print *, "    states: ", states
      print *, "    r0: ", Rod%r(:,0)
      print *, "    q : ", Rod%q     

      
   END SUBROUTINE Rod_Initialize
   !--------------------------------------------------------------




   ! set kinematics for Rods ONLY if they are attached to a body (including a coupled body) or coupled (otherwise shouldn't be called)
   !--------------------------------------------------------------
   SUBROUTINE Rod_SetKinematics(Rod, r6_in, v6_in, a6_in, t, LineList)

      Type(MD_Rod),     INTENT(INOUT)  :: Rod            ! the Rod object
      Real(DbKi),       INTENT(IN   )  :: r6_in(6)       ! 6-DOF position
      Real(DbKi),       INTENT(IN   )  :: v6_in(6)       ! 6-DOF velocity
      Real(DbKi),       INTENT(IN   )  :: a6_in(6)       ! 6-DOF acceleration (only used for coupled rods)
      Real(DbKi),       INTENT(IN   )  :: t              ! instantaneous time
      Type(MD_Line),    INTENT(INOUT)  :: LineList(:)    ! array of all the Line objects

      INTEGER(IntKi)                   :: l

      Rod%time = t    ! store current time

      
      if (abs(Rod%typeNum) == 2) then ! rod rigidly coupled to a body, or ground, or coupling point
         Rod%r6 = r6_in
         Rod%v6 = v6_in
         Rod%a6 = a6_in
         
         call ScaleVector(Rod%r6(4:6), 1.0_DbKi, Rod%r6(4:6)); ! enforce direction vector to be a unit vector
         
         ! since this rod has no states and all DOFs have been set, pass its kinematics to dependent Lines
         CALL Rod_SetDependentKin(Rod, t, LineList)
      
      else if (abs(Rod%typeNum) == 1) then ! rod end A pinned to a body, or ground, or coupling point
      
         ! set Rod *end A only* kinematics based on BCs (linear model for now) 
         Rod%r6(1:3) = r6_in(1:3)
         Rod%v6(1:3) = v6_in(1:3)
         Rod%a6(1:3) = a6_in(1:3)

         
         ! Rod is pinned so only end A is specified, rotations are left alone and will be 
         ! handled, along with passing kinematics to dependent lines, by separate call to setState
      
      else
         print *, "Error: Rod_SetKinematics called for a free Rod."  ! <<<
      end if
	
   
      ! update Rod direction unit vector (simply equal to last three entries of r6, presumably these were set elsewhere for pinned Rods)
		 Rod%q = Rod%r6(4:6)
      
	   
   END SUBROUTINE Rod_SetKinematics
   !--------------------------------------------------------------

   ! pass the latest states to the rod if it has any DOFs/states (then update rod end kinematics including attached lines)
   !--------------------------------------------------------------
   SUBROUTINE Rod_SetState(Rod, X, t, LineList)

      Type(MD_Rod),          INTENT(INOUT)  :: Rod        ! the Rod object
      Real(DbKi),            INTENT(IN   )  :: X(:)           ! state vector section for this line
      Real(DbKi),            INTENT(IN   )  :: t              ! instantaneous time
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:)    ! array of all the Line objects

      INTEGER(IntKi)                        :: J              ! index
   

      ! for a free Rod, there are 12 states:
      ! [ x, y, z velocity of end A, then rate of change of u/v/w coordinates of unit vector pointing toward end B,
      ! then x, y, z coordinate of end A, u/v/w coordinates of unit vector pointing toward end B]

      ! for a pinned Rod, there are 6 states (rotational only):
      ! [ rate of change of u/v/w coordinates of unit vector pointing toward end B,
      ! then u/v/w coordinates of unit vector pointing toward end B]
      
      
      ! store current time
      Rod%time = t


      ! copy over state values for potential use during derivative calculations
      if (Rod%typeNum == 0) then                         ! free Rod type
      
         ! CALL ScaleVector(X(10:12), 1.0, X(10:12))  ! enforce direction vector to be a unit vector <<<< can't do this with FAST frameowrk, could be a problem!!
         
         ! TODO: add "controller" adjusting state derivatives of X(10:12) to artificially force X(10:12) to remain a unit vector <<<<<<<<<<<

         
         Rod%r6(1:3) = X(7:9)                         ! (end A coordinates)
         Rod%v6(1:3) = X(1:3)                         ! (end A velocity, unrotated axes) 
         CALL ScaleVector(X(10:12), 1.0_DbKi, Rod%r6(4:6)) !Rod%r6(4:6) = X(10:12)                    ! (Rod direction unit vector)
         Rod%v6(4:6) = X(4:6)                         ! (rotational velocities about unrotated axes) 
         
         
         CALL Rod_SetDependentKin(Rod, t, LineList)
      
      else if (abs(Rod%typeNum) == 1) then                       ! pinned rod type (coupled or attached to something)t previously via setPinKin)
      
         !CALL ScaleVector(X(4:6), 1.0, X(4:6))      ! enforce direction vector to be a unit vector
         
         
         CALL ScaleVector(X(4:6), 1.0_DbKi, Rod%r6(4:6)) !Rod%r6(3+J) = X(3+J) ! (Rod direction unit vector)
         Rod%v6(4:6) = X(1:3)                    ! (rotational velocities about unrotated axes) 
         
         
         CALL Rod_SetDependentKin(Rod, t, LineList)
      
      else
         print *, "Error: Rod::setState called for a non-free rod type"   ! <<<
      end if

      ! update Rod direction unit vector (simply equal to last three entries of r6)
      Rod%q = Rod%r6(4:6)
      
   END SUBROUTINE Rod_SetState
   !--------------------------------------------------------------


   ! Set the Rod end kinematics then set the kinematics of dependent objects (any attached lines).
   ! This also determines the orientation of zero-length rods.
   !--------------------------------------------------------------
   SUBROUTINE Rod_SetDependentKin(Rod, t, LineList)

      Type(MD_Rod),          INTENT(INOUT)  :: Rod            ! the Rod object
      Real(DbKi),            INTENT(IN   )  :: t              ! instantaneous time
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:)    ! array of all the Line objects

      INTEGER(IntKi)                        :: l              ! index of segments or nodes along line
      INTEGER(IntKi)                        :: J              ! index
      INTEGER(IntKi)                        :: N              ! number of segments
   
      REAL(DbKi)                            :: qEnd(3)        ! unit vector of attached line end segment, following same direction convention as Rod's q vector
      REAL(DbKi)                            :: EIend          ! bending stiffness of attached line end segment
      REAL(DbKi)                            :: dlEnd          ! stretched length of attached line end segment
      REAL(DbKi)                            :: qMomentSum(3)  ! summation of qEnd*EI/dl_stretched (with correct sign) for each attached line
         
         
      ! in future pass accelerations here too? <<<<
   
      N = Rod%N

      ! from state values, set positions of end nodes 
      ! end A
      Rod%r(:,0)  = Rod%r6(1:3)  ! positions
      Rod%rd(:,0) = Rod%v6(1:3)  ! velocities
      
      !print *, Rod%r6(1:3)
      !print *, Rod%r(:,0)
      
      if (Rod%N > 0) then  ! set end B nodes only if the rod isn't zero length
         call transformKinematicsAtoB(Rod%r6(1:3), Rod%r6(4:6), Rod%UnstrLen, Rod%v6, Rod%r(:,N), Rod%rd(:,N))   ! end B    
      end if

      ! pass end node kinematics to any attached lines (this is just like what a Connection does, except for both ends)
      DO l=1,Rod%nAttachedA
         CALL Line_SetEndKinematics(LineList(Rod%attachedA(l)), Rod%r(:,0), Rod%rd(:,0), t, Rod%TopA(l))
      END DO
      DO l=1,Rod%nAttachedB
         CALL Line_SetEndKinematics(LineList(Rod%attachedB(l)), Rod%r(:,N), Rod%rd(:,N), t, Rod%TopB(l))
      END DO


      ! if this is a zero-length Rod, get bending moment-related information from attached lines and compute Rod's equilibrium orientation
      if (N==0) then
      
         DO l=1,Rod%nAttachedA
         
            CALL Line_GetEndSegmentInfo(LineList(Rod%attachedA(l)), qEnd, EIend, dlEnd, Rod%TopA(l))
            
            qMomentSum = qMomentSum + qEnd*EIend/dlEnd  ! add each component to the summation vector
            
         END DO

         DO l=1,Rod%nAttachedB
         
            CALL Line_GetEndSegmentInfo(LineList(Rod%attachedB(l)), qEnd, EIend, dlEnd, Rod%TopB(l))
            
            qMomentSum = qMomentSum + qEnd*EIend/dlEnd  ! add each component to the summation vector
            
         END DO
         
         ! solve for line unit vector that balances all moments (unit vector of summation of qEnd*EI/dl_stretched over each line)
         CALL ScaleVector(qMomentSum, 1.0_DbKi, Rod%q)
      END IF

      ! pass Rod orientation to any attached lines (this is just like what a Connection does, except for both ends)
      DO l=1,Rod%nAttachedA
         CALL Line_SetEndOrientation(LineList(Rod%attachedA(l)), Rod%q, Rod%TopA(l), 0)
      END DO
      DO l=1,Rod%nAttachedB
         CALL Line_SetEndOrientation(LineList(Rod%attachedB(l)), Rod%q, Rod%TopB(l), 1)
      END DO
      
   END SUBROUTINE Rod_SetDependentKin
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Rod_GetStateDeriv(Rod, Xd, LineList, p)

      Type(MD_Rod),          INTENT(INOUT)  :: Rod              ! the Rod object
      Real(DbKi),            INTENT(INOUT)  :: Xd(:)            ! state derivative vector section for this line
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:)      ! array of all the Line objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p                ! Parameters
      
      !TYPE(MD_MiscVarType), INTENT(INOUT)  :: m       ! misc/optimization variables

      INTEGER(IntKi)                        :: J                ! index
      
      Real(DbKi)                            :: Fnet     (6)     ! net force and moment about reference point
      Real(DbKi)                            :: M_out    (6,6)   ! mass matrix about reference point
      
      Real(DbKi)                            :: acc(6)           ! 6DOF acceleration vector about reference point
      
      Real(DbKi)                            :: Mcpl(3)          ! moment in response to end A acceleration due to inertial coupling
      
      Real(DbKi)                            :: y_temp (6)       ! temporary vector for LU decomposition
      Real(DbKi)                            :: LU_temp(6,6)     ! temporary matrix for LU decomposition
      

      CALL Rod_GetNetForceAndMass(Rod, Rod%r(:,0), Fnet, M_out, LineList, p)
                  
                  

   ! TODO: add "controller" adjusting state derivatives of X(10:12) to artificially force X(10:12) to remain a unit vector <<<<<<<<<<<

      ! fill in state derivatives	
      IF (Rod%typeNum == 0) THEN                         ! free Rod type, 12 states  
         
         ! solve for accelerations in [M]{a}={f} using LU decomposition
         CALL LUsolve(6, M_out, LU_temp, Fnet, y_temp, acc)
         
         Xd(7:9) = Rod%v6(1:3)  !Xd[6 + I] = v6[  I];       ! dxdt = V   (velocities)
         Xd(1:6) = acc          !Xd[    I] = acc[  I];      ! dVdt = a   (accelerations) 
                                !Xd[3 + I] = acc[3+I];        ! rotational accelerations	
      
         ! rate of change of unit vector components!!  CHECK!   <<<<<
         Xd(10) =                - Rod%v6(6)*Rod%r6(5) + Rod%v6(5)*Rod%r6(6) ! i.e.  u_dot_x = -omega_z*u_y + omega_y*u_z
         Xd(11) =  Rod%v6(6)*Rod%r6(4)                 - Rod%v6(4)*Rod%r6(6) ! i.e.  u_dot_y =  omega_z*u_x - omega_x*u_z
         Xd(12) = -Rod%v6(5)*Rod%r6(4) + Rod%v6(4)*Rod%r6(5)                 ! i.e.  u_dot_z = -omega_y*u_x - omega_x*u_y

         ! store accelerations in case they're useful as output
         Rod%a6 = acc

      ELSE                            ! pinned rod, 6 states (rotational only)
      
         ! account for moment in response to end A acceleration due to inertial coupling (off-diagonal sub-matrix terms)
         !Fnet(4:6) = Fnet(4:6) - MATMUL(M_out(4:6,1:3), Rod%a6(1:3))  ! <<<check that it's the right submatrix <<<
         Fnet(4:6) = Fnet(4:6) - MATMUL(M_out(1:3,4:6), Rod%a6(1:3))  ! <<< THIS order is stable. Weird. <<<
         ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ the above line seems to be causing the stability problems for USFLOWT! <<<<
         
         ! solve for accelerations in [M]{a}={f} using LU decomposition
         CALL LUsolve(3, M_out(4:6,4:6), LU_temp(4:6,4:6), Fnet(4:6), y_temp(4:6), acc(4:6))
         ! Note: solving for rotational DOFs only - excluding translational and off-diagonal 3x3 terms -
         
         Xd(1:3) = acc(4:6)   !   Xd[    I] = acc[3+I];          ! rotational accelerations
         
         ! rate of change of unit vector components!!  CHECK!   <<<<<
         Xd(4) =                - Rod%v6(6)*Rod%r6(5) + Rod%v6(5)*Rod%r6(6) ! i.e.  u_dot_x = -omega_z*u_y + omega_y*u_z
         Xd(5) =  Rod%v6(6)*Rod%r6(4)                 - Rod%v6(4)*Rod%r6(6) ! i.e.  u_dot_y =  omega_z*u_x - omega_x*u_z
         Xd(6) = -Rod%v6(5)*Rod%r6(4) + Rod%v6(4)*Rod%r6(5)                 ! i.e.  u_dot_z = -omega_y*u_x - omega_x*u_y
      
         ! store angular accelerations in case they're useful as output
         Rod%a6(4:6) = acc(4:6)
      
      END IF	
      
      ! Note: accelerations that are dependent on parent objects) will not be known to this object 
      !       (only those of free DOFs are coupled DOFs are known in this approach).
   
      ! check for NaNs (should check all state derivatives, not just first 6)
      DO J = 1, 6
         IF (Is_NaN(REAL(Xd(J),DbKi))) THEN
            print *, "NaN detected at time ", Rod%time, " in Rod ",Rod%IdNum, " state derivatives:"
            print *, Xd
            
            print *, "r0"
            print *, Rod%r(:,0)
            print *, "F"
            print *, Fnet
            print *, "M"
            print *, M_out
            print *, "acc"
            print *, acc            
            
            EXIT
         END IF
      END DO

   END SUBROUTINE Rod_GetStateDeriv
   !--------------------------------------------------------------


   ! calculate the aggregate 3/6DOF rigid-body loads of a coupled rod including inertial loads
   !--------------------------------------------------------------
   SUBROUTINE Rod_GetCoupledForce(Rod, Fnet_out, LineList, p)

      Type(MD_Rod),          INTENT(INOUT)  :: Rod         ! the Rod object
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(6) ! force and moment vector
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:) ! array of all the Line objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters
      
      Real(DbKi)                            :: F6_iner(6)   ! inertial reaction force
      
      ! do calculations of forces and masses on each rod node
      CALL Rod_DoRHS(Rod, LineList, p)

      ! add inertial loads as appropriate (written out in a redundant way just for clarity, and to support load separation in future)
      ! fixed coupled rod
      if (Rod%typeNum == -2) then                          
      
         F6_iner  = -MATMUL(Rod%M6net, Rod%a6)    ! inertial loads      
         Fnet_out = Rod%F6net + F6_iner           ! add inertial loads
      
      ! pinned coupled rod      
      else if (Rod%typeNum == -1) then                     
         ! inertial loads ... from input translational ... and solved rotational ... acceleration
         F6_iner(4:6)  = -MATMUL(Rod%M6net(1:3,1:3), Rod%a6(1:3)) - MATMUL(Rod%M6net(1:3,4:6), Rod%a6(4:6))
         Fnet_out(1:3) = Rod%F6net(1:3) + F6_iner(4:6)     ! add translational inertial loads
         Fnet_out(4:6) = 0.0_DbKi
      else
         print *, "ERROR, Rod_GetCoupledForce called for wrong (non-coupled) rod type!"
      end if
   
   END SUBROUTINE Rod_GetCoupledForce
   !--------------------------------------------------------------
   


   ! calculate the aggregate 6DOF rigid-body force and mass data of the rod 
   !--------------------------------------------------------------
   SUBROUTINE Rod_GetNetForceAndMass(Rod, rRef, Fnet_out, M_out, LineList, p)

      Type(MD_Rod),          INTENT(INOUT)  :: Rod         ! the Rod object
      Real(DbKi),            INTENT(IN   )  :: rRef(3)     ! global coordinates of reference point (end A for free Rods)
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(6) ! force and moment vector about rRef
      Real(DbKi),            INTENT(  OUT)  :: M_out(6,6)  ! mass and inertia matrix about rRef
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:) ! array of all the Line objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters
      
      Real(DbKi)                 :: rRel(  3)              ! relative position of each node i from rRef      
      
      ! do calculations of forces and masses on each rod node
      CALL Rod_DoRHS(Rod, LineList, p)

      ! note: Some difference from MoorDyn C here. If this function is called by the Rod itself, the reference point must be end A

      ! shift everything from end A reference to rRef reference point
      
      rRel = Rod%r(:,0) - rRef   ! vector from reference point to end A            
         
      CALL translateForce3to6DOF(rRel, Rod%F6net(1:3), Fnet_out)	   ! shift net forces
      Fnet_out(4:6) = Fnet_out(4:6) + Rod%F6net(4:6)               ! add in the existing moments
         
      CALL translateMass6to6DOF(rRel, Rod%M6net, M_out)          ! shift mass matrix to be about ref point
         
      ! >>> do we need to ensure zero moment is passed if it's pinned? <<<
      !if (abs(Rod%typeNum)==1) then
      !   Fnet_out(4:6) = 0.0_DbKi
      !end if
	
   
   END SUBROUTINE Rod_GetNetForceAndMass
   !--------------------------------------------------------------
   

   ! calculate the forces on the rod, including from attached lines
   !--------------------------------------------------------------
   SUBROUTINE Rod_DoRHS(Rod, LineList, p)

      Type(MD_Rod),          INTENT(INOUT)  :: Rod            ! the Rodion object
      Type(MD_Line),         INTENT(INOUT)  :: LineList(:)    ! array of all the Line objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p              ! Parameters
      
      !TYPE(MD_MiscVarType), INTENT(INOUT)  :: m       ! misc/optimization variables

      INTEGER(IntKi)             :: l            ! index of attached lines
      INTEGER(IntKi)             :: I,J,K        ! index
      
      
      INTEGER(IntKi)             :: N            ! number of rod elements for convenience

      Real(DbKi)                 :: phi, beta, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta   ! various orientation things
      Real(DbKi)                 :: k_hat(3)     ! unit vector (redundant, not used) <<<<
      Real(DbKi)                 :: Ftemp        ! temporary force component
      Real(DbKi)                 :: Mtemp        ! temporary moment component

      Real(DbKi)                 :: m_i, v_i     ! 
      Real(DbKi)                 :: zeta         ! wave elevation above/below a given node
      Real(DbKi)                 :: h0           ! distance along rod centerline from end A to the waterplane     
      Real(DbKi)                 :: deltaL       ! submerged length of a given segment
      Real(DbKi)                 :: Lsum         ! cumulative length along rod axis from bottom
      Real(DbKi)                 :: dL           ! length attributed to node
      Real(DbKi)                 :: VOF          ! fraction of volume associated with node that is submerged
      
      Real(DbKi)                 :: Vi(3)        ! relative flow velocity over a node
      Real(DbKi)                 :: SumSqVp, SumSqVq, MagVp, MagVq
      Real(DbKi)                 :: Vp(3), Vq(3) ! transverse and axial components of water velocity at a given node     
      Real(DbKi)                 :: ap(3), aq(3) ! transverse and axial components of water acceleration at a given node
      Real(DbKi)                 :: Fnet_i(3)    ! force from an attached line
      Real(DbKi)                 :: Mnet_i(3)    ! moment from an attached line
      Real(DbKi)                 :: Mass_i(3,3)  ! mass from an attached line

      ! used in lumped 6DOF calculations:
      Real(DbKi)                 :: rRel(  3)              ! relative position of each node i from rRef      
      Real(DbKi)                 :: OrMat(3,3)             ! rotation matrix to rotate global z to rod's axis
      Real(DbKi)                 :: F6_i(6)                ! a node's contribution to the total force vector
      Real(DbKi)                 :: M6_i(6,6)              ! a node's contribution to the total mass matrix
      Real(DbKi)                 :: I_l                    ! axial inertia of rod
      Real(DbKi)                 :: I_r                    ! radial inertia of rod about CG
      Real(DbKi)                 :: Imat_l(3,3)            ! inertia about CG aligned with Rod axis
      Real(DbKi)                 :: Imat(3,3)              ! inertia about CG in global frame     
      Real(DbKi)                 :: h_c                    ! location of CG along axis
      Real(DbKi)                 :: r_c(3)                 ! 3d location of CG relative to node A      
      Real(DbKi)                 :: Fcentripetal(3)        ! centripetal force
      Real(DbKi)                 :: Mcentripetal(3)        ! centripetal moment


      N = Rod%N

      ! ------------------------------ zero some things --------------------------
      
      Rod%Mext = 0.0_DbKi  ! zero the external moment sum

      Lsum = 0.0_DbKi

      
      ! ---------------------------- initial rod and node calculations ------------------------

      ! calculate some orientation information for the Rod as a whole
      call GetOrientationAngles(Rod%r( :,0), Rod%r( :,N), phi, sinPhi, cosPhi, tanPhi, beta, sinBeta, cosBeta, k_hat)
 
      ! save to internal roll and pitch variables for use in output <<< should check these, make Euler angles isntead of independent <<<
      Rod%roll  = -180.0/Pi * phi*sinBeta
      Rod%pitch =  180.0/Pi * phi*cosBeta

      ! set interior node positions and velocities (stretch the nodes between the endpoints linearly) (skipped for zero-length Rods)
      DO i=1,N-1
         Rod%r( :,i) =  Rod%r( :,0) + (Rod%r( :,N) - Rod%r( :,0)) * (REAL(i)/REAL(N))
         Rod%rd(:,i) =  Rod%rd(:,0) + (Rod%rd(:,N) - Rod%rd(:,0)) * (REAL(i)/REAL(N))
         
      
         Rod%V(i) = 0.25*pi * Rod%d*Rod%d * Rod%l(i) ! volume attributed to segment
      END DO


   ! --------------------------------- apply wave kinematics ------------------------------------

      IF (p%WaterKin == 1)  THEN ! wave kinematics interpolated from global grid in Waves object
         DO i=0,N
            CALL getWaveKin(p, Rod%r(1,i), Rod%r(2,i), Rod%r(3,i), Rod%time, Rod%U(:,i), Rod%Ud(:,i), Rod%zeta(i), Rod%PDyn(i))
            !F(i) = 1.0 ! set VOF value to one for now (everything submerged - eventually this should be element-based!!!) <<<<
            ! <<<< currently F is not being used and instead a VOF variable is used within the node loop
         END DO
      END IF


    !  ! wave kinematics not implemented yet <<<
    !  ap = 0.0_DbKi
    !  aq = 0.0_DbKi
    !  ! set U and Ud herem as well as pDyn and zeta...
    !  Rod%U    = 0.0_DbKi
    !  Rod%Ud   = 0.0_DbKi
    !  pDyn = 0.0_DbKi
    !  zeta = 0.0_DbKi
      
      ! >>> remember to check for violated conditions, if there are any... <<<
           
      zeta = Rod%zeta(N)! just use the wave elevation computed at the location of the top node for now
      
      if ((Rod%r(3,0) < zeta) .and. (Rod%r(3,N) > zeta)) then    ! check if it's crossing the water plane (should also add some limits to avoid near-horizontals at some point)
         h0 = (zeta - Rod%r(3,0))/Rod%q(3)                       ! distance along rod centerline from end A to the waterplane
      else if (Rod%r(3,0) < zeta) then
         h0 = 2.0*Rod%UnstrLen                                   ! fully submerged case
      else
         h0 = 0.0_DbKi                                           ! fully unsubmerged case (ever applicable?)
      end if

   
      ! -------------------------- loop through all the nodes -----------------------------------
      DO I = 0, N
      
      
         ! ------------------ calculate added mass matrix for each node -------------------------
      
         ! get mass and volume considering adjacent segment lengths
         IF (I==0) THEN
            dL  = 0.5*Rod%l(1)
            m_i = 0.25*Pi * Rod%d*Rod%d * dL *Rod%rho     ! (will be zero for zero-length Rods)
            v_i = 0.5 *Rod%V(1)
         ELSE IF (I==N) THEN
            dL  = 0.5*Rod%l(N)
            m_i = 0.25*pi * Rod%d*Rod%d * dL *Rod%rho
            v_i = 0.5*Rod%V(N)
         ELSE
            dL  = 0.5*(Rod%l(I) + Rod%l(I+1))
            m_i = 0.25*pi * Rod%d*Rod%d * dL *Rod%rho
            v_i = 0.5 *(Rod%V(I) + Rod%V(I+1))
         END IF

         ! get scalar for submerged portion                  
         IF (Lsum + dL < h0) THEN    ! if fully submerged 
            VOF = 1.0_DbKi
         ELSE IF (Lsum < h0) THEN    ! if partially below waterline 
            VOF = (h0 - Lsum)/dL
         ELSE                        ! must be out of water
            VOF = 0.0_DbKi
         END IF
         
         Lsum = Lsum + dL            ! add length attributed to this node to the total

         ! build mass and added mass matrix
         DO J=1,3
            DO K=1,3
               IF (J==K) THEN
                  Rod%M(K,J,I) = m_i + VOF*p%rhoW*v_i*( Rod%Can*(1 - Rod%q(J)*Rod%q(K)) + Rod%Cat*Rod%q(J)*Rod%q(K) )
               ELSE
                  Rod%M(K,J,I) = VOF*p%rhoW*v_i*( Rod%Can*(-Rod%q(J)*Rod%q(K)) + Rod%Cat*Rod%q(J)*Rod%q(K) )
               END IF
            END DO
         END DO
         
         ! <<<< what about accounting for offset of half segment from node location for end nodes? <<<<
         
         
!         CALL Inverse3by3(Rod%S(:,:,I), Rod%M(:,:,I))             ! invert mass matrix


         ! ------------------  CALCULATE FORCES ON EACH NODE ----------------------------

         if (N > 0) then ! the following force calculations are only nonzero for finite-length rods (skipping for zero-length Rods)
         
            ! >>> no nodal axial elasticity loads calculated since it's assumed rigid, but should I calculate tension/compression due to other loads? <<<

            ! weight (now only the dry weight)
            Rod%W(:,I) = (/ 0.0_DbKi, 0.0_DbKi, -m_i * p%g /)   ! assuming g is positive
            
            ! buoyance (now calculated based on outside pressure, for submerged portion only)
            ! radial buoyancy force from sides
            Ftemp = -VOF * 0.25*Pi*dL*Rod%d*Rod%d * p%rhoW*p%g * sinPhi
            Rod%Bo(:,I) = (/ Ftemp*cosBeta*cosPhi, Ftemp*sinBeta*cosPhi, -Ftemp*sinPhi /)            

            !relative flow velocities
            DO J = 1, 3
               Vi(J) = Rod%U(J,I) - Rod%rd(J,I)                               ! relative flow velocity over node -- this is where wave velicites would be added
            END DO

            ! decomponse relative flow into components
            SumSqVp = 0.0_DbKi                                         ! start sums of squares at zero
            SumSqVq = 0.0_DbKi
            DO J = 1, 3
               Vq(J) = DOT_PRODUCT( Vi , Rod%q ) * Rod%q(J);            ! tangential relative flow component
               Vp(J) = Vi(J) - Vq(J)                                    ! transverse relative flow component
               SumSqVq = SumSqVq + Vq(J)*Vq(J)
               SumSqVp = SumSqVp + Vp(J)*Vp(J)
            END DO
            MagVp = sqrt(SumSqVp)                                       ! get magnitudes of flow components
            MagVq = sqrt(SumSqVq)

            ! transverse and tangenential drag
            Rod%Dp(:,I) = VOF * 0.5*p%rhoW*Rod%Cdn*    Rod%d* dL * MagVp * Vp
            Rod%Dq(:,I) = 0.0_DbKi ! 0.25*p%rhoW*Rod%Cdt* Pi*Rod%d* dL * MagVq * Vq <<< should these axial side loads be included?

            ! fluid acceleration components for current node
            aq = DOT_PRODUCT(Rod%Ud(:,I), Rod%q) * Rod%q  ! tangential component of fluid acceleration
            ap = Rod%Ud(:,I) - aq                         ! normal component of fluid acceleration
            ! transverse Froude-Krylov force
            Rod%Ap(:,I) = VOF * p%rhoW*(1.0+Rod%Can)*0.5* v_i * ap  ! <<< are these equations right??
            ! axial Froude-Krylov force	
            Rod%Aq(:,I) = 0.0_DbKi  ! p%rhoW*(1.0+Rod%Cat)*0.5* v_i * aq  ! <<< are these equations right??

            ! dynamic pressure
            Rod%Pd(:,I) = 0.0_DbKi  ! assuming zero for sides
            
            ! bottom contact (stiffness and damping, vertical-only for now)  - updated Nov 24 for general case where anchor and fairlead ends may deal with bottom contact forces
            IF (Rod%r(3,I) < -p%WtrDpth) THEN
               IF (I==0) THEN
                  Rod%B(3,I) = ( (-p%WtrDpth - Rod%r(3,I))*p%kBot - Rod%rd(3,I)*p%cBot) * 0.5*Rod%d*(            Rod%l(I+1) ) 
               ELSE IF (I==N) THEN
                  Rod%B(3,I) = ( (-p%WtrDpth - Rod%r(3,I))*p%kBot - Rod%rd(3,I)*p%cBot) * 0.5*Rod%d*(Rod%l(I)               ) 
               ELSE
                  Rod%B(3,I) = ( (-p%WtrDpth - Rod%r(3,I))*p%kBot - Rod%rd(3,I)*p%cBot) * 0.5*Rod%d*(Rod%l(I) + Rod%l(I+1) ) 
               END IF
            ELSE
               Rod%B(3,I) = 0.0_DbKi
            END IF
            
         ELSE    ! zero-length (N=0) Rod case
         
            ! >>>>>>>>>>>>>> still need to check handling of zero length rods <<<<<<<<<<<<<<<<<<<
         
            ! for zero-length rods, make sure various forces are zero
            Rod%W  = 0.0_DbKi
            Rod%Bo = 0.0_DbKi
            Rod%Dp = 0.0_DbKi
            Rod%Dq= 0.0_DbKi
            Rod%B = 0.0_DbKi
            Rod%Pd = 0.0_DbKi
            
         END IF
         
         
         ! ------ now add forces, moments, and added mass from Rod end effects (these can exist even if N==0) -------
         
         ! end A
         IF ((I==0) .and. (h0 > 0.0_ReKi)) THEN    ! if this is end A and it is submerged 
         
         ! >>> eventually should consider a VOF approach for the ends    hTilt = 0.5*Rod%d/cosPhi <<<
         
            ! buoyancy force
            Ftemp = -VOF * 0.25*Pi*Rod%d*Rod%d * p%rhoW*p%g*Rod%r(3,I)
            Rod%Bo(:,I) = Rod%Bo(:,I) + (/ Ftemp*cosBeta*sinPhi, Ftemp*sinBeta*sinPhi, Ftemp*cosPhi /) 
         
            ! buoyancy moment
            Mtemp = -VOF * 1.0/64.0*Pi*Rod%d**4 * p%rhoW*p%g * sinPhi 
            Rod%Mext = Rod%Mext + (/ Mtemp*sinBeta, -Mtemp*cosBeta, 0.0_DbKi /) 
         
            ! axial drag
            Rod%Dq(:,I) = Rod%Dq(:,I) + VOF * 0.25* Pi*Rod%d*Rod%d * p%rhoW*Rod%CdEnd * MagVq * Vq
            
if ((Rod%time >= 1.0) .and. (Rod%time < 1.001)) then
   print *, "at Dq end 0 of rod:"
   print *, "CdEnd is  on position vector:"
   print *, Rod%CdEnd
end if
            
            ! >>> what about rotational drag?? <<<   eqn will be  Pi* Rod%d**4/16.0 omega_rel?^2...  *0.5 * Cd...

            ! Froud-Krylov force
            Rod%Aq(:,I) = Rod%Aq(:,I) + VOF * p%rhoW*(1.0+Rod%CaEnd)*0.5* (2.0/3.0*Pi*Rod%d**3 /8) * aq
            
            ! dynamic pressure force
            Rod%Pd(:,I) = Rod%Pd(:,I) + VOF * 0.25* Pi*Rod%d*Rod%d * Rod%PDyn(I) * Rod%q
            
            ! added mass
            DO J=1,3
               DO K=1,3
                  IF (J==K) THEN
                     Rod%M(K,J,I) = Rod%M(K,J,I) + VOF*p%rhoW* (Pi*Rod%d**3/6.0) * Rod%CaEnd*Rod%q(J)*Rod%q(K) 
                  ELSE
                     Rod%M(K,J,I) = Rod%M(K,J,I) + VOF*p%rhoW* (Pi*Rod%d**3/6.0) * Rod%CaEnd*Rod%q(J)*Rod%q(K) 
                  END IF
               END DO
            END DO
         
         END IF
            
         IF ((I==N) .and. (h0 > Rod%UnstrLen)) THEN    ! if this end B and it is submerged (note, if N=0, both this and previous if statement are true)
         
            ! buoyancy force
            Ftemp = VOF * 0.25*Pi*Rod%d*Rod%d * p%rhoW*p%g*Rod%r(3,I)
            Rod%Bo(:,I) = Rod%Bo(:,I) + (/ Ftemp*cosBeta*sinPhi, Ftemp*sinBeta*sinPhi, Ftemp*cosPhi /) 
         
            ! buoyancy moment
            Mtemp = VOF * 1.0/64.0*Pi*Rod%d**4 * p%rhoW*p%g * sinPhi 
            Rod%Mext = Rod%Mext + (/ Mtemp*sinBeta, -Mtemp*cosBeta, 0.0_DbKi /) 
            
            ! axial drag
            Rod%Dq(:,I) = Rod%Dq(:,I) + VOF * 0.25* Pi*Rod%d*Rod%d * p%rhoW*Rod%CdEnd * MagVq * Vq
            
            ! Froud-Krylov force
            Rod%Aq(:,I) = Rod%Aq(:,I) + VOF * p%rhoW*(1.0+Rod%CaEnd)*0.5* (2.0/3.0*Pi*Rod%d**3 /8) * aq
            
            ! dynamic pressure force
            Rod%Pd(:,I) = Rod%Pd(:,I) - VOF * 0.25* Pi*Rod%d*Rod%d * Rod%PDyn(I) * Rod%q
            
            ! added mass
            DO J=1,3
               DO K=1,3
                  IF (J==K) THEN
                     Rod%M(K,J,I) = Rod%M(K,J,I) + VOF*p%rhoW* (Pi*Rod%d**3/6.0) * Rod%CaEnd*Rod%q(J)*Rod%q(K) 
                  ELSE
                     Rod%M(K,J,I) = Rod%M(K,J,I) + VOF*p%rhoW* (Pi*Rod%d**3/6.0) * Rod%CaEnd*Rod%q(J)*Rod%q(K) 
                  END IF
               END DO
            END DO
            
         END IF
         
         
         
         ! ---------------------------- total forces for this node -----------------------------
         
         Rod%Fnet(:,I) = Rod%W(:,I) + Rod%Bo(:,I) + Rod%Dp(:,I) + Rod%Dq(:,I) &
                         + Rod%Ap(:,I) + Rod%Aq(:,I) + Rod%Pd(:,I) + Rod%B(:,I)
         
	
      END DO  ! I  - done looping through nodes

	
      ! ----- add waterplane moment of inertia moment if applicable -----
      IF ((Rod%r(3,0) < zeta) .and. (Rod%r(3,N) > zeta)) then    ! check if it's crossing the water plane
         Mtemp = 1.0/16.0 *Pi*Rod%d**4 * p%rhoW*p%g * sinPhi * (1.0 + 0.5* tanPhi**2)
         Rod%Mext = Rod%Mext + (/ Mtemp*sinBeta, -Mtemp*cosBeta, 0.0_DbKi /)
      END IF
   
      ! ---------------- now add in forces on end nodes from attached lines ------------------
         
      ! loop through lines attached to end A
      DO l=1,Rod%nAttachedA
         
         CALL Line_GetEndStuff(LineList(Rod%attachedA(l)), Fnet_i, Mnet_i, Mass_i, Rod%TopA(l))
         
         ! sum quantitites
         Rod%Fnet(:,0)= Rod%Fnet(:,0) + Fnet_i    ! total force
         Rod%Mext     = Rod%Mext      + Mnet_i    ! externally applied moment
         Rod%M(:,:,0) = Rod%M(:,:,0)  + Mass_i    ! mass at end node
         
      END DO
   
      ! loop through lines attached to end B
      DO l=1,Rod%nAttachedB
         
         CALL Line_GetEndStuff(LineList(Rod%attachedB(l)), Fnet_i, Mnet_i, Mass_i, Rod%TopB(l))
         
         ! sum quantitites
         Rod%Fnet(:,N)= Rod%Fnet(:,N) + Fnet_i    ! total force
         Rod%Mext     = Rod%Mext      + Mnet_i    ! externally applied moment
         Rod%M(:,:,N) = Rod%M(:,:,N)  + Mass_i    ! mass at end node
         
      END DO
      
      ! ---------------- now lump everything in 6DOF about end A -----------------------------
	
      ! question: do I really want to neglect the rotational inertia/drag/etc across the length of each segment?
   
      ! make sure 6DOF quantiaties are zeroed before adding them up
      Rod%F6net = 0.0_DbKi
      Rod%M6net = 0.0_DbKi

      ! now go through each node's contributions, put them about end A, and sum them
      DO i = 0,Rod%N
      
         rRel = Rod%r(:,i) - Rod%r(:,0)   ! vector from reference point to node            
         
         ! convert segment net force into 6dof force about body ref point (if the Rod itself, end A)
         CALL translateForce3to6DOF(rRel, Rod%Fnet(:,i), F6_i)			
         
         ! convert segment mass matrix to 6by6 mass matrix about body ref point  (if the Rod itself, end A)
         CALL translateMass3to6DOF(rRel, Rod%M(:,:,i), M6_i)
                  
         ! sum contributions
         Rod%F6net = Rod%F6net + F6_i
         Rod%M6net = Rod%M6net + M6_i
         
      END DO
      
      ! ------------- Calculate some items for the Rod as a whole here -----------------
      
      ! >>> could some of these be precalculated just once? <<<
            
      ! add inertia terms for the Rod assuming it is uniform density (radial terms add to existing matrix which contains parallel-axis-theorem components only)
      I_l = 0.125*Rod%mass * Rod%d*Rod%d     ! axial moment of inertia
      I_r = Rod%mass/12 * (0.75*Rod%d*Rod%d + (Rod%UnstrLen/Rod%N)**2 ) * Rod%N     ! summed radial moment of inertia for each segment individually
      
      !h_c = [value from registry]

      Imat_l(1,1) = I_r   ! inertia about CG in local orientations (as if Rod is vertical)
      Imat_l(2,2) = I_r
      Imat_l(3,3) = I_l
      
      OrMat = CalcOrientation(phi, beta, 0.0_DbKi)        ! get rotation matrix to put things in global rather than rod-axis orientations
      
      Imat = RotateM3(Imat_l, OrMat)  ! rotate to give inertia matrix about CG in global frame
      
      ! these supplementary inertias can then be added the matrix (these are the terms ASIDE from the parallel axis terms)
      Rod%M6net(4:6,4:6) = Rod%M6net(4:6,4:6) + Imat
      

      ! now add centripetal and gyroscopic forces/moments, and that should be everything
      h_c = 0.5*Rod%UnstrLen          ! distance to center of mass
      r_c = h_c*Rod%q                 ! vector to center of mass
      
      ! note that Rod%v6(4:6) is the rotational velocity vector, omega   
      Fcentripetal = 0.0_DbKi !<<<TEMP<<< -cross_product(Rod%v6(4:6), cross_product(Rod%v6(4:6), r_c ))*Rod%mass <<<
      Mcentripetal = 0.0_DbKi !<<<TEMP<<< cross_product(r_c, Fcentripetal) - cross_product(Rod%v6(4:6), MATMUL(Imat,Rod%v6(4:6)))
      
      ! add centripetal force/moment, gyroscopic moment, and any moments applied from lines at either end (might be zero)
      Rod%F6net(1:3) = Rod%F6net(1:3) + Fcentripetal 
      Rod%F6net(4:6) = Rod%F6net(4:6) + Mcentripetal + Rod%Mext
            
      ! Note: F6net saves the Rod's net forces and moments (excluding inertial ones) for use in later output
      !       (this is what the rod will apply to whatever it's attached to, so should be zero moments if pinned).
      !       M6net saves the rod's mass matrix.
      
      
   

   END SUBROUTINE Rod_DoRHS
   !=====================================================================




   ! this function handles assigning a line to a connection node
   SUBROUTINE Rod_AddLine(Rod, lineID, TopOfLine, endB)

      Type(MD_Rod), INTENT (INOUT)   :: Rod        ! the Connection object

      Integer(IntKi),   INTENT( IN )     :: lineID
      Integer(IntKi),   INTENT( IN )     :: TopOfLine
      Integer(IntKi),   INTENT( IN )     :: endB   ! add line to end B if 1, end A if 0

      if (endB==1) then   ! attaching to end B

         Print*, "L", lineID, "->R", Rod%IdNum , "b"
         
         IF (Rod%nAttachedB <10) THEN ! this is currently just a maximum imposed by a fixed array size.  could be improved.
            Rod%nAttachedB = Rod%nAttachedB + 1  ! add the line to the number connected
            Rod%AttachedB(Rod%nAttachedB) = lineID
            Rod%TopB(Rod%nAttachedB) = TopOfLine  ! attached to line ... 1 = top/fairlead(end B), 0 = bottom/anchor(end A)
         ELSE
            Print*, "too many lines connected!"
         END IF

      else              ! attaching to end A
      
         Print*, "L", lineID, "->R", Rod%IdNum , "a"
         
         IF (Rod%nAttachedA <10) THEN ! this is currently just a maximum imposed by a fixed array size.  could be improved.
            Rod%nAttachedA = Rod%nAttachedA + 1  ! add the line to the number connected
            Rod%AttachedA(Rod%nAttachedA) = lineID
            Rod%TopA(Rod%nAttachedA) = TopOfLine  ! attached to line ... 1 = top/fairlead(end B), 0 = bottom/anchor(end A)
         ELSE
            Print*, "too many lines connected!"
         END IF
         
      end if

   END SUBROUTINE Rod_AddLine


   ! this function handles removing a line from a connection node
   SUBROUTINE Rod_RemoveLine(Rod, lineID, TopOfLine, endB,  rEnd, rdEnd)

      Type(MD_Rod), INTENT (INOUT)  :: Rod        ! the Connection object

      Integer(IntKi),   INTENT( IN )     :: lineID
      Integer(IntKi),   INTENT(  OUT)    :: TopOfLine
      Integer(IntKi),   INTENT( IN )     :: endB   ! end B if 1, end A if 0
      REAL(DbKi),       INTENT(INOUT)    :: rEnd(3)
      REAL(DbKi),       INTENT(INOUT)    :: rdEnd(3)
      
      Integer(IntKi)    :: l,m,J
      
      if (endB==1) then   ! attaching to end B
         
         DO l = 1,Rod%nAttachedB    ! look through attached lines
         
            IF (Rod%AttachedB(l) == lineID) THEN   ! if this is the line's entry in the attachment list
            
               TopOfLine = Rod%TopB(l);                ! record which end of the line was attached
               
               DO m = l,Rod%nAttachedB-1 
               
                  Rod%AttachedB(m) = Rod%AttachedB(m+1)  ! move subsequent line links forward one spot in the list to eliminate this line link
                  Rod%TopB(     m) =      Rod%TopB(m+1) 
               
                  Rod%nAttachedB = Rod%nAttachedB - 1                      ! reduce attached line counter by 1
               
                  ! also pass back the kinematics at the end
                  DO J = 1,3
                     rEnd( J) = Rod%r( J,Rod%N)
                     rdEnd(J) = Rod%rd(J,Rod%N)
                  END DO
                  
                  print*, "Detached line ", lineID, " from Rod ", Rod%IdNum, " end B"
                  
                  EXIT
               END DO
               
               IF (l == Rod%nAttachedB) THEN   ! detect if line not found
                  print *, "Error: failed to find line to remove during RemoveLine call to Rod ", Rod%IdNum, ". Line ", lineID
               END IF
            END IF
         END DO
         
      else              ! attaching to end A
              
        DO l = 1,Rod%nAttachedA    ! look through attached lines
         
            IF (Rod%AttachedA(l) == lineID) THEN   ! if this is the line's entry in the attachment list
            
               TopOfLine = Rod%TopA(l);                ! record which end of the line was attached
               
               DO m = l,Rod%nAttachedA-1 
               
                  Rod%AttachedA(m) = Rod%AttachedA(m+1)  ! move subsequent line links forward one spot in the list to eliminate this line link
                  Rod%TopA(     m) =      Rod%TopA(m+1) 
               
                  Rod%nAttachedA = Rod%nAttachedA - 1                      ! reduce attached line counter by 1
               
                  ! also pass back the kinematics at the end
                  DO J = 1,3
                     rEnd( J) = Rod%r( J,0)
                     rdEnd(J) = Rod%rd(J,0)
                  END DO
                  
                  print*, "Detached line ", lineID, " from Rod ", Rod%IdNum, " end A"
                  
                  EXIT
               END DO
               
               IF (l == Rod%nAttachedA) THEN   ! detect if line not found
                  print *, "Error: failed to find line to remove during RemoveLine call to Rod ", Rod%IdNum, ". Line ", lineID
               END IF
            END IF
         END DO
      
      end if
      
   END SUBROUTINE Rod_RemoveLine








!--------------------------------------------------------------
!            Body-Specific Subroutines
!--------------------------------------------------------------


!   ! used to initialize bodies that aren't free i.e. don't have states
!   !--------------------------------------------------------------
!   SUBROUTINE Body_InitializeUnfree(Body, r6_in, mesh, mesh_index, m)
!
!      Type(MD_Body),         INTENT(INOUT)  :: Body        ! the Body object
!      Real(DbKi),            INTENT(IN   )  :: r6_in(6)    ! state vector section for this line
!      TYPE(MeshType),        INTENT(INOUT)  :: mesh        !
!      Integer(IntKi),        INTENT(IN   )  :: mesh_index  ! index of the node in the mesh for the current object being initialized
!      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
!
!      INTEGER(IntKi)                        :: l           ! index of segments or nodes along line
!      REAL(DbKi)                            :: rRef(3)     ! reference position of mesh node
!      REAL(DbKi)                            :: OrMat(3,3)  ! DCM for body orientation based on r6_in
!      REAL(DbKi)                            :: dummyStates(12) 
!   
!   
!      rRef = 0.0_DbKi   ! <<< maybe this should be the offsets of the local platform origins from the global origins in future? And that's what's specificed by the Body input coordinates?
!      
!      CALL MeshPositionNode(mesh, mesh+index, rRef,ErrStat2,ErrMsg2)! "assign the coordinates (u%PtFairleadDisplacement%Position) of each node in the global coordinate space"
!
!      CALL CheckError( ErrStat2, ErrMsg2 )
!      IF (ErrStat >= AbortErrLev) RETURN
!
!      ! Apply offsets due to initial platform rotations and translations (fixed Jun 19, 2015)
!      CALL SmllRotTrans('body rotation matrix', r6_in(4),r6_in(5),r6_in(6), OrMat, '', ErrStat2, ErrMsg2)
!      mesh%TranslationDisp(1, mesh_index) = r6_in(1) + OrMat(1,1)*rRef(1) + OrMat(2,1)*rRef(2) + OrMat(3,1)*rRef(3) - rRef(1)
!      mesh%TranslationDisp(2, mesh_index) = r6_in(2) + OrMat(1,2)*rRef(1) + OrMat(2,2)*rRef(2) + OrMat(3,2)*rRef(3) - rRef(2)
!      mesh%TranslationDisp(3, mesh_index) = r6_in(3) + OrMat(1,3)*rRef(1) + OrMat(2,3)*rRef(2) + OrMat(3,3)*rRef(3) - rRef(3)
!
!      ! what about node point orientation ???
!
!      ! If any Rod is fixed to the body (not pinned), initialize it now because otherwise it won't be initialized
!      DO l=1, Body%nAttachedR
!         if (m%RodList(Body%attachedR(l))%typeNum == 2)  CALL Rod_Initialize(m%RodList(Body%attachedR(l)), dummyStates, m%LineList)
!      END DO
!      
!      ! Note: Connections don't need any initialization
!      
!   END SUBROUTINE Body_InitializeUnfree
!   !--------------------------------------------------------------


   ! used to initialize bodies that are free
   !--------------------------------------------------------------
   SUBROUTINE Body_Initialize(Body, states, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body            ! the Body object
      Real(DbKi),            INTENT(INOUT)  :: states(6)       ! state vector section for this Body
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m               ! passing along all mooring objects

      INTEGER(IntKi)                        :: l               ! index of segments or nodes along line
      REAL(DbKi)                            :: dummyStates(12) ! dummy vector to mimic states when initializing a rigidly attached rod
   
   
      ! assign initial body kinematics to state vector
      states(7:12) = Body%r6
      states(1:6 ) = Body%v6
      	

      ! set positions of any dependent connections and rods now (before they are initialized)
      CALL Body_SetDependentKin(Body, 0.0_DbKi, m)
            
      ! If any Rod is fixed to the body (not pinned), initialize it now because otherwise it won't be initialized
      DO l=1, Body%nAttachedR
         if (m%RodList(Body%attachedR(l))%typeNum == 2)  CALL Rod_Initialize(m%RodList(Body%attachedR(l)), dummyStates,  m%LineList)
      END DO
      
      ! Note: Connections don't need any initialization
      
   END SUBROUTINE Body_Initialize
   !--------------------------------------------------------------
   
   ! used to initialize bodies that are coupled or fixed
   !--------------------------------------------------------------
   SUBROUTINE Body_InitializeUnfree(Body, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body            ! the Body object
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m               ! passing along all mooring objects

      INTEGER(IntKi)                        :: l               ! index of segments or nodes along line
      REAL(DbKi)                            :: dummyStates(12) ! dummy vector to mimic states when initializing a rigidly attached rod
   
   
      ! set positions of any dependent connections and rods now (before they are initialized)
      CALL Body_SetDependentKin(Body, 0.0_DbKi, m)
            
      ! If any Rod is fixed to the body (not pinned), initialize it now because otherwise it won't be initialized
      DO l=1, Body%nAttachedR
         if (m%RodList(Body%attachedR(l))%typeNum == 2)  CALL Rod_Initialize(m%RodList(Body%attachedR(l)), dummyStates,  m%LineList)
      END DO
      
      ! Note: Connections don't need any initialization
      
   END SUBROUTINE Body_InitializeUnfree
   !--------------------------------------------------------------



   !--------------------------------------------------------------
   SUBROUTINE Body_SetState(Body, X, t, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body           ! the Body object
      Real(DbKi),            INTENT(IN   )  :: X(:)           ! state vector section for this line
      Real(DbKi),            INTENT(IN   )  :: t              ! instantaneous time
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m              ! passing along all mooring objects

      INTEGER(IntKi)                        :: l              ! index of segments or nodes along line
      INTEGER(IntKi)                        :: J              ! index
   
      ! store current time
      Body%time = t
      
      
      
      Body%r6 = X(7:12)   ! get positions      
      Body%v6 = X(1:6)    ! get velocities
      

      ! set positions of any dependent connections and rods
      CALL Body_SetDependentKin(Body, t, m)
      
   END SUBROUTINE Body_SetState
   !--------------------------------------------------------------


   ! set kinematics for Bodies if they are coupled (or ground)
   !--------------------------------------------------------------
   SUBROUTINE Body_SetKinematics(Body, r_in, v_in, a_in, t, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body       ! the Body object
      Real(DbKi),            INTENT(IN   )  :: r_in(6)   ! 6-DOF position
      Real(DbKi),            INTENT(IN   )  :: v_in(6)   ! 6-DOF velocity
      Real(DbKi),             INTENT(IN   ) :: a_in(6)       ! 6-DOF acceleration (only used for coupled rods)
      Real(DbKi),            INTENT(IN   )  :: t         ! instantaneous time
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m         ! passing along all mooring objects (for simplicity, since Bodies deal with Rods and Connections)


      INTEGER(IntKi)                   :: l

      ! store current time
      Body%time = t

   !   if (abs(Body%typeNum) == 2) then ! body coupled in 6 DOF, or ground
         Body%r6 = r_in
         Body%v6 = v_in
         Body%a6 = a_in
                  
         ! since this body has no states and all DOFs have been set, pass its kinematics to dependent attachments
         CALL Body_SetDependentKin(Body, t, m)
      
   !   else if (abs(Body%typeNum) == 1) then ! body pinned at reference point
   !   
   !      ! set Body *end A only* kinematics based on BCs (linear model for now) 
   !      Body%r6(1:3) = r_in(1:3)
   !      Body%v6(1:3) = v_in(1:3)
   !      
   !      ! Body is pinned so only ref point posiiton is specified, rotations are left alone and will be 
   !      ! handled, along with passing kinematics to attached objects, by separate call to setState
   !   
   !   else
   !      print *, "Error: Body_SetKinematics called for a free Body."  ! <<<
   !   end if
	 
   END SUBROUTINE Body_SetKinematics
   !--------------------------------------------------------------

	
   ! set the states (positions and velocities) of any connects or rods that are part of this body
   ! also computes the orientation matrix (never skip this sub!)
   !--------------------------------------------------------------
   SUBROUTINE Body_SetDependentKin(Body, t, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body        ! the Bodyion object
      REAL(DbKi),            INTENT(IN   )  :: t
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects (for simplicity, since Bodies deal with Rods and Connections)

      INTEGER(IntKi)                        :: l              ! index of attached objects
   
      Real(DbKi)                            :: rConnect(3)
      Real(DbKi)                            :: rdConnect(3)
      Real(DbKi)                            :: rRod(6)
      Real(DbKi)                            :: vRod(6)
      Real(DbKi)                            :: aRod(6)

      

      ! calculate orientation matrix based on latest angles
      !CALL SmllRotTrans('', Body%r6(4), Body%r6(5), Body%r6(6), Body%TransMat, '', ErrStat2, ErrMsg2)
      Body%OrMat = EulerConstruct( Body%r6(4:6) )  ! full Euler angle approach <<<< need to check order 
  
if ((t >= 1.0) .and. (t < 1.001)) then
   print *, "orientation matrix OrMat of Body:"
   print *, Body%OrMat
   print *, "based on position vector:"
   print *, Body%r6
end if

      ! set kinematics of any dependent connections
      do l = 1,Body%nAttachedC
      
         CALL transformKinematics(Body%rConnectRel(:,l), Body%r6, Body%OrMat, Body%v6, rConnect, rdConnect) !<<< should double check this function
                  
         ! >>> need to add acceleration terms here too? <<<
                  
         ! pass above to the connection and get it to calculate the forces
         CALL Connect_SetKinematics( m%ConnectList(Body%attachedC(l)), rConnect, rdConnect, m%zeros6(1:3), t, m%LineList)
      end do
      
      ! set kinematics of any dependent Rods
      do l=1,Body%nAttachedR
      
         ! calculate displaced coordinates/orientation and velocities of each rod <<<<<<<<<<<<<
         ! do 3d details of Rod ref point
         CALL TransformKinematicsA( Body%r6RodRel(1:3,l), Body%r6(1:3), Body%OrMat, Body%v6, Body%a6, rRod(1:3), vRod(1:3), aRod(1:3))  ! set first three entires (end A translation) of rRod and rdRod
         ! does the above function need to take in all 6 elements of r6RodRel??
         
         ! do rotational stuff	
         rRod(4:6) = MATMUL(Body%OrMat, Body%r6RodRel(4:6,l))    !<<<<<< correct? <<<<< rotateVector3(r6RodRel[i]+3, OrMat, rRod+3);   ! rotate rod relative unit vector by OrMat to get unit vec in reference coords
         vRod(4:6) = Body%v6(4:6)  ! transformed rotational velocity.  <<< is this okay as is? <<<<
         aRod(4:6) = Body%a6(4:6) 
         
         ! pass above to the rod and get it to calculate the forces
         CALL Rod_SetKinematics(m%RodList(Body%attachedR(l)), rRod, vRod, aRod, t, m%LineList)
      end do

   END SUBROUTINE Body_SetDependentKin
   !--------------------------------------------------------------
   
      ! calculate the aggregate 3/6DOF rigid-body loads of a coupled rod including inertial loads
   !--------------------------------------------------------------
   SUBROUTINE Body_GetCoupledForce(Body, Fnet_out, m, p)

      Type(MD_Body),         INTENT(INOUT)  :: Body        ! the Body object
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(6) ! force and moment vector
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters
      
      Real(DbKi)                            :: F6_iner(6)  ! inertial reaction force
      
      ! do calculations of forces and masses on the body
      CALL Body_DoRHS(Body, m, p)

      ! add inertial loads as appropriate 
      if (Body%typeNum == -1) then                          
      
         F6_iner = 0.0_DbKi !-MATMUL(Body%M, Body%a6)     <<<<<<<< why does including F6_iner cause instability???
         Fnet_out = Body%F6net + F6_iner        ! add inertial loads
         
      else
         print *, "ERROR, Body_GetCoupledForce called for wrong (non-coupled) body type!"
      end if
   
   END SUBROUTINE Body_GetCoupledForce
   !--------------------------------------------------------------
   

   !--------------------------------------------------------------
   SUBROUTINE Body_GetStateDeriv(Body, Xd, m, p)

      Type(MD_Body),         INTENT(INOUT)  :: Body          ! the Bodyion object
      Real(DbKi),            INTENT(INOUT)  :: Xd(:)            ! state derivative vector section for this line
      
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p                ! Parameters
      
      INTEGER(IntKi)                        :: J                ! index
      
      Real(DbKi)                            :: acc(6)           ! 6DOF acceleration vector
      
      Real(DbKi)                            :: y_temp (6)       ! temporary vector for LU decomposition
      Real(DbKi)                            :: LU_temp(6,6)     ! temporary matrix for LU decomposition
      

      CALL Body_DoRHS(Body, m, p)

      ! solve for accelerations in [M]{a}={f} using LU decomposition
      CALL LUsolve(6, Body%M, LU_temp, Body%F6net, y_temp, acc)

      ! fill in state derivatives	      
      Xd(7:12) = Body%v6       ! dxdt = V   (velocities)
      Xd(1:6)  = acc           ! dVdt = a   (accelerations) 

      ! store accelerations in case they're useful as output
      Body%a6 = acc
   
      ! check for NaNs (should check all state derivatives, not just first 6)
      DO J = 1, 6
         IF (Is_NaN(REAL(Xd(J),DbKi))) THEN
            print *, "NaN detected at time ", Body%time, " in Body ",Body%IdNum, " state derivatives:"
            print *, Xd
            EXIT
         END IF
      END DO


   END SUBROUTINE Body_GetStateDeriv
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Body_DoRHS(Body, m, p)

      Type(MD_Body),         INTENT(INOUT)  :: Body        ! the Bodyion object
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters
      
      !TYPE(MD_MiscVarType), INTENT(INOUT)  :: m       ! misc/optimization variables

      INTEGER(IntKi)             :: l         ! index of attached lines
      INTEGER(IntKi)             :: I         ! index
      INTEGER(IntKi)             :: J         ! index
      INTEGER(IntKi)             :: K         ! index

      Real(DbKi)                 :: Fgrav(3)           ! body weight force
      Real(DbKi)                 :: body_rCGrotated(3) ! instantaneous vector from body ref point to CG
      Real(DbKi)                 :: U(3)               ! water velocity - zero for now
      Real(DbKi)                 :: Ud(3)              ! water acceleration- zero for now
      Real(DbKi)                 :: vi(6)              ! relative water velocity (last 3 terms are rotatonal and will be set to zero
      Real(DbKi)                 :: F6_i(6)            ! net force and moments from an attached object
      Real(DbKi)                 :: M6_i(6,6)          ! mass and inertia from an attached object

      ! First, the body's own mass matrix must be adjusted based on its orientation so that 
      ! we have a mass matrix in the global orientation frame	
      Body%M = RotateM6(Body%M0, Body%OrMat)
	
      !gravity on core body
	
      Fgrav(3) = Body%bodyV * p%rhow * p%g - Body%bodyM * p%g ! weight+buoyancy vector
	
      body_rCGrotated = MATMUL(Body%OrMat, Body%rCG) ! rotateVector3(body_rCG, OrMat, body_rCGrotated); ! relative vector to body CG in inertial orientation
      CALL translateForce3to6DOF(body_rCGrotated, Fgrav, Body%F6net)  ! gravity forces and moments about body ref point given CG location
	
	
      ! --------------------------------- apply wave kinematics ------------------------------------
      !env->waves->getU(r6, t, U); ! call generic function to get water velocities <<<<<<<<< all needs updating
	
      !	for (int J=0; J<3; J++)		
      !		Ud[J] = 0.0;                 ! set water accelerations as zero for now
      ! ------------------------------------------------------------------------------------------

      ! viscous drag calculation (on core body)
      vi(1:3) = U - Body%v6(1:3)  ! relative flow velocity over body ref point
      vi(4:6) =   - Body%v6(4:6)  ! for rotation, this is just the negative of the body's rotation for now (not allowing flow rotation)
	
      Body%F6net = Body%F6net + 0.5*p%rhoW * vi * abs(vi) * Body%bodyCdA
      ! <<< NOTE, for body this should be fixed to account for orientation!! <<< what about drag in rotational DOFs??? <<<<<<<<<<<<<<
	
   
   
      ! Get contributions from any dependent connections
      do l = 1,Body%nAttachedC
      
         ! get net force and mass from Connection on body ref point (global orientation)
         CALL Connect_GetNetForceAndMass( m%ConnectList(Body%attachedC(l)), Body%r6(1:3), F6_i, M6_i, m%LineList, p)
         
         ! sum quantitites
         Body%F6net = Body%F6net + F6_i
         Body%M     = Body%M     + M6_i
		 end do
      
      ! Get contributions from any dependent Rods
      do l=1,Body%nAttachedR
      
         ! get net force and mass from Rod on body ref point (global orientation)
         CALL Rod_GetNetForceAndMass(m%RodList(Body%attachedR(l)), Body%r6(1:3), F6_i, M6_i, m%LineList, p)
         
         ! sum quantitites
         Body%F6net = Body%F6net + F6_i
         Body%M     = Body%M     + M6_i
      end do
	

   END SUBROUTINE Body_DoRHS
   !=====================================================================




   ! this function handles assigning a connection to a body
   !--------------------------------------------------------------
   SUBROUTINE Body_AddConnect(Body, connectID, coords)

      Type(MD_Body),      INTENT(INOUT)  :: Body        ! the Connection object
      Integer(IntKi),     INTENT(IN   )  :: connectID
      REAL(DbKi),         INTENT(IN   )  :: coords(3)


      Print*, "C", connectID, "->B", Body%IdNum
      
      IF(Body%nAttachedC < 30) THEN                ! this is currently just a maximum imposed by a fixed array size.  could be improved.
         Body%nAttachedC = Body%nAttachedC + 1     ! increment the number connected
         Body%AttachedC(Body%nAttachedC) = connectID
         Body%rConnectRel(:,Body%nAttachedC) = coords  ! store relative position of connect on body
      ELSE
         Print*, "too many connections attached!"
      END IF

   END SUBROUTINE Body_AddConnect


   ! this function handles assigning a rod to a body
   !--------------------------------------------------------------
   SUBROUTINE Body_AddRod(Body, rodID, coords)

      Type(MD_Body),      INTENT(INOUT)  :: Body        ! the Connection object
      Integer(IntKi),     INTENT(IN   )  :: rodID
      REAL(DbKi),         INTENT(IN   )  :: coords(6)  ! positions of rod ends A and B relative to body
      
      REAL(DbKi)                         :: tempUnitVec(3)
      REAL(DbKi)                         :: dummyLength

      Print*, "R", rodID, "->B", Body%IdNum
      
      IF(Body%nAttachedR < 30) THEN                ! this is currently just a maximum imposed by a fixed array size.  could be improved.
         Body%nAttachedR = Body%nAttachedR + 1     ! increment the number connected
         
         ! store rod ID
         Body%AttachedR(Body%nAttachedR) = rodID   
         
         ! store Rod end A relative position and unit vector from end A to B
         CALL UnitVector(coords(1:3), coords(4:6), tempUnitVec, dummyLength)
         Body%r6RodRel(1:3, Body%nAttachedR) = coords(1:3)
         Body%r6RodRel(4:6, Body%nAttachedR) = tempUnitVec
         
      ELSE
         Print*, "too many rods attached!"
      END IF

   END SUBROUTINE Body_AddRod


   ! :::::::::::::::::::::::::: below are some wave related functions :::::::::::::::::::::::::::::::
   
   
   ! master function to get wave/water kinematics at a given point -- called by each object fro grid-based data
   SUBROUTINE getWaveKin(p, x, y, z, t, U, Ud, zeta, PDyn)
   
      ! note, this whole approach assuems that px, py, and pz are in increasing order
   
      TYPE(MD_ParameterType),INTENT (IN   )       :: p       ! MoorDyn parameters (contains the wave info for now)
      Real(DbKi),            INTENT (IN   )       :: x
      Real(DbKi),            INTENT (IN   )       :: y
      Real(DbKi),            INTENT (IN   )       :: z
      Real(DbKi),            INTENT (IN   )       :: t
      Real(DbKi),            INTENT (INOUT)       :: U(3)
      Real(DbKi),            INTENT (INOUT)       :: Ud(3)
      Real(DbKi),            INTENT (INOUT)       :: zeta
      Real(DbKi),            INTENT (INOUT)       :: PDyn


      INTEGER(IntKi)             :: ix, iy, iz, it        ! indeces for interpolation      
      INTEGER(IntKi)             :: N                     ! number of rod elements for convenience
      Real(DbKi)                 :: fx, fy, fz, ft        ! interpolation fractions
      Real(DbKi)                 :: qt                    ! used in time step interpolation
   
      
      CALL getInterpNums(p%px, x, ix, fx)
      CALL getInterpNums(p%py, y, iy, fy)
      CALL getInterpNums(p%pz, z, iz, fz)
      
      qt = t / real(p%dtWave, DbKi)
      it = floor(qt) + 1   ! adjust by 1 for fortran's indexing starting at 1
      ft = qt - it + 1.0; !(t-(it*dtWave))/dtWave  ! //remainder(t,dtWave)/dtWave;	
      
      CALL calculate3Dinterpolation(p%zeta, ix, iy, it, fx, fy, ft, zeta)
      
      CALL calculate4Dinterpolation(p%PDyn, ix, iy, iz, it, fx, fy, fz, ft, PDyn)
      
      CALL calculate4Dinterpolation(p%ux, ix, iy, iz, it, fx, fy, fz, ft, U(1)  )
      CALL calculate4Dinterpolation(p%uy, ix, iy, iz, it, fx, fy, fz, ft, U(2)  )
      CALL calculate4Dinterpolation(p%uz, ix, iy, iz, it, fx, fy, fz, ft, U(3)  )      
      CALL calculate4Dinterpolation(p%ax, ix, iy, iz, it, fx, fy, fz, ft, Ud(1) )
      CALL calculate4Dinterpolation(p%ay, ix, iy, iz, it, fx, fy, fz, ft, Ud(2) )
      CALL calculate4Dinterpolation(p%az, ix, iy, iz, it, fx, fy, fz, ft, Ud(3) )
      
   END SUBROUTINE
   
   
   SUBROUTINE getInterpNums(xlist, xin, i, fout)
      
      Real(DbKi),    INTENT (IN   )            :: xlist(:)
      Real(DbKi),    INTENT (IN   )            :: xin
      Integer(IntKi),INTENT (  OUT)            :: i
      Real(DbKi),    INTENT (  OUT)            :: fout
      
      Integer(IntKi)                           :: nx
      ! Parameters: list of x values, number of x values, x value to be interpolated, fraction to return
      ! Returns the lower index to interpolate from.  such that  y* = y[i] + fout*(y[i+1]-y[i])
      
      nx = SIZE(xlist)
      
      if (xin <= xlist(1)) THEN                !  below lowest data point
         i = 1_IntKi
         fout = 0.0_DbKi
      
      else if (xin >= xlist(nx)) THEN        ! above highest data point
         i = nx
         fout = 0.0_DbKi
      
      else                                     ! within the data range
         DO i = 1, nx-1
         	  IF (xlist(i+1) > xin) THEN
               fout = (xin - xlist(i) )/( xlist(i+1) - xlist(i) )
               exit
            END IF
         END DO		
      END IF
      
	END SUBROUTINE
   
  
   SUBROUTINE calculate4Dinterpolation(f, ix0, iy0, iz0, it0, fx, fy, fz, ft, c)

      Real(DbKi),     INTENT (IN   )        :: f(:,:,:,:)                ! data array
      INTEGER(IntKi), INTENT (IN   )        :: ix0, iy0, iz0, it0        ! indeces for interpolation
      Real(DbKi),     INTENT (IN   )        :: fx, fy, fz, ft            ! interpolation fractions
      Real(DbKi),     INTENT (  OUT)        :: c                         ! the output value
                                         
      INTEGER(IntKi)                        :: ix1, iy1, iz1, it1        ! second indices
      REAL(DbKi)                            :: c000, c001, c010, c011, c100, c101, c110, c111
      REAL(DbKi)                            :: c00, c01, c10, c11, c0, c1  
      
      ! handle end case conditions
      if (fx == 0) then 
         ix1 = ix0
      else  
         ix1 = ix0+1
      end if
      
      if (fy == 0) then
         iy1 = iy0
      else
         iy1 = iy0+1
      end if
      
      if (fz == 0) then
         iz1 = iz0
      else         
         iz1 = iz0+1
      end if
      
      if (ft == 0) then
         it1 = it0
      else  
         it1 = it0+1
      end if
      
      c000 = f(it0,iz0,iy0,ix0)*(1-ft) + f(it1,iz0,iy0,ix0)*ft
      c001 = f(it0,iz1,iy0,ix0)*(1-ft) + f(it1,iz1,iy0,ix0)*ft
      c010 = f(it0,iz0,iy1,ix0)*(1-ft) + f(it1,iz0,iy1,ix0)*ft
      c011 = f(it0,iz1,iy1,ix0)*(1-ft) + f(it1,iz1,iy1,ix0)*ft
      c100 = f(it0,iz0,iy0,ix1)*(1-ft) + f(it1,iz0,iy0,ix1)*ft
      c101 = f(it0,iz1,iy0,ix1)*(1-ft) + f(it1,iz1,iy0,ix1)*ft
      c110 = f(it0,iz0,iy1,ix1)*(1-ft) + f(it1,iz0,iy1,ix1)*ft
      c111 = f(it0,iz1,iy1,ix1)*(1-ft) + f(it1,iz1,iy1,ix1)*ft
      
      c00 = c000*(1-fx) + c100*fx
      c01 = c001*(1-fx) + c101*fx
      c10 = c010*(1-fx) + c110*fx
      c11 = c011*(1-fx) + c111*fx

      c0  = c00 *(1-fy) + c10 *fy
      c1  = c01 *(1-fy) + c11 *fy

      c   = c0  *(1-fz) + c1  *fz
            
   END SUBROUTINE


   SUBROUTINE calculate3Dinterpolation(f, ix0, iy0, iz0, fx, fy, fz, c)

         Real(DbKi),     INTENT (IN   )        :: f(:,:,:)                  ! data array
         INTEGER(IntKi), INTENT (IN   )        :: ix0, iy0, iz0             ! indeces for interpolation
         Real(DbKi),     INTENT (IN   )        :: fx, fy, fz                ! interpolation fractions
         Real(DbKi),     INTENT (  OUT)        :: c                         ! the output value
         
         INTEGER(IntKi)                        :: ix1, iy1, iz1             ! second indices
         REAL(DbKi)                            :: c000, c001, c010, c011, c100, c101, c110, c111
         REAL(DbKi)                            :: c00, c01, c10, c11, c0, c1  
         
      ! note that "z" could also be "t" - dimension names are arbitrary

      ! handle end case conditions
      if (fx == 0) then 
         ix1 = ix0
      else  
         ix1 = ix0+1
      end if
      
      if (fy == 0) then
         iy1 = iy0
      else
         iy1 = iy0+1
      end if
      
      if (fz == 0) then
         iz1 = iz0
      else         
         iz1 = iz0+1
      end if
      
      c000 = f(iz0,iy0,ix0)
      c001 = f(iz1,iy0,ix0)
      c010 = f(iz0,iy1,ix0)
      c011 = f(iz1,iy1,ix0)
      c100 = f(iz0,iy0,ix1)
      c101 = f(iz1,iy0,ix1)
      c110 = f(iz0,iy1,ix1)
      c111 = f(iz1,iy1,ix1)
      
      c00 = c000*(1-fx) + c100*fx
      c01 = c001*(1-fx) + c101*fx
      c10 = c010*(1-fx) + c110*fx
      c11 = c011*(1-fx) + c111*fx

      c0  = c00 *(1-fy) + c10 *fy
      c1  = c01 *(1-fy) + c11 *fy

      c   = c0  *(1-fz) + c1  *fz

   END SUBROUTINE



   ! ============ below are some math convenience functions ===============
   ! should add error checking if I keep these, but hopefully there are existing NWTCLib functions to replace them


   ! return unit vector (u) and in direction from r1 to r2 and distance between points
   !-----------------------------------------------------------------------
   SUBROUTINE UnitVector( r1, r2, u, Length )          ! note: order of parameters chagned in this function
      
      REAL(DbKi),       INTENT(IN   )  :: r1(:)
      REAL(DbKi),       INTENT(IN   )  :: r2(:)
      REAL(DbKi),       INTENT(  OUT)  :: u(:)
      REAL(DbKi),       INTENT(  OUT)  :: length

      u = r2 - r1
      length = TwoNorm(u)

      if ( .NOT. EqualRealNos(length, 0.0_DbKi ) ) THEN
        u = u / Length
      END IF

   END SUBROUTINE UnitVector
   !-----------------------------------------------------------------------

   ! scale vector to desired length 
   !-----------------------------------------------------------------------
   SUBROUTINE ScaleVector( u_in, newlength, u_out )
      REAL(DbKi),       INTENT(IN   )  :: u_in(3)     ! input vector
      REAL(DbKi),       INTENT(IN   )  :: newlength   ! desired length of output vector
      REAL(DbKi),       INTENT(INOUT)  :: u_out(3)    ! output vector (hopefully can be the same as u_in without issue)

      REAL(DbKi)                       :: length_squared
      REAL(DbKi)                       :: scaler
      INTEGER(IntKi)                   :: J
      
      length_squared = 0.0;	
      DO J=1,3
         length_squared = length_squared + u_in(J)*u_in(J)				
      END DO
      
      if (length_squared > 0) then
         scaler = newlength/sqrt(length_squared)	
      else                   ! if original vector is zero, return zero
         scaler = 0_DbKi
      end if
      
      DO J=1,3
         u_out(J) = u_in(J) * scaler 
      END DO

   END SUBROUTINE ScaleVector
   !-----------------------------------------------------------------------


   ! calculate orientation angles of a cylindrical object
   !-----------------------------------------------------------------------
   subroutine GetOrientationAngles(p1, p2, phi, sinPhi, cosPhi, tanPhi, beta, sinBeta, cosBeta, k_hat)
      real(DbKi),   intent(in   ) :: p1(3),p2(3)
      real(DbKi),   intent(  out) :: phi, sinPhi, cosPhi, tanPhi, beta, sinBeta, cosBeta, k_hat(3)
            
      real(DbKi)                  :: vec(3), vecLen, vecLen2D

      ! calculate isntantaneous incline angle and heading, and related trig values
      ! the first and last NodeIndx values point to the corresponding Joint nodes idices which are at the start of the Mesh
      vec      = p2 - p1   
      vecLen   = SQRT(Dot_Product(vec,vec))
      vecLen2D = SQRT(vec(1)**2+vec(2)**2)
      if ( vecLen < 0.000001 ) then
         print *, "ERROR in GetOrientationAngles" !call SeterrStat(ErrID_Fatal, 'An element of the Morison structure has co-located endpoints!  This should never occur.  Please review your model.', errStat, errMsg, 'Morison_CalcOutput' )
         print *, p1
         print *, p2
         k_hat = 1.0/0.0
      else
         k_hat = vec / vecLen 
         phi   = atan2(vecLen2D, vec(3))  ! incline angle   
      end if
      if ( phi < 0.000001) then
         beta = 0.0_ReKi
      else
         beta = atan2(vec(2), vec(1))                    ! heading of incline     
      endif
      sinPhi  = sin(phi)
      cosPhi  = cos(phi)  
      tanPhi  = tan(phi)     
      sinBeta = sin(beta)
      cosBeta = cos(beta)
            
   end subroutine GetOrientationAngles
   !-----------------------------------------------------------------------
   

   ! calculate position and velocity of point based on its position relative to moving 6DOF body
   !-----------------------------------------------------------------------
   SUBROUTINE TransformKinematics(rRelBody, r_in, TransMat, rd_in, r_out, rd_out)
      REAL(DbKi),       INTENT(IN   )  :: rRelBody(:)  ! coordinate of end A
      REAL(DbKi),       INTENT(IN   )  :: r_in(3)      ! Rod unit vector
      REAL(DbKi),       INTENT(IN   )  :: TransMat(3,3)! 
      REAL(DbKi),       INTENT(IN   )  :: rd_in(6)     ! 6DOF velecity vector about Rod end A, in global orientation frame
      REAL(DbKi),       INTENT(  OUT)  :: r_out(3)     ! coordinates of end B
      REAL(DbKi),       INTENT(  OUT)  :: rd_out(3)    ! velocity of end B

      REAL(DbKi)                       :: rRel(3)

      ! rd_in should be in global orientation frame
      ! note: it's okay if r_out and rd_out are 6-size. Only the first 3 will be written, and 4-6 will
      !       already be correct or can be assigned seperately from r_in and rd_in (assuming orientation frames are identical)


      ! locations (unrotated reference frame) about platform reference point  (2021-01-05: just transposed TransMat, it was incorrect before)
      rRel(1) = TransMat(1,1)*rRelBody(1) + TransMat(1,2)*rRelBody(2) + TransMat(1,3)*rRelBody(3) ! x
      rRel(2) = TransMat(2,1)*rRelBody(1) + TransMat(2,2)*rRelBody(2) + TransMat(2,3)*rRelBody(3) ! y
      rRel(3) = TransMat(3,1)*rRelBody(1) + TransMat(3,2)*rRelBody(2) + TransMat(3,3)*rRelBody(3) ! z

      ! absolute locations
      r_out = rRel + r_in

      ! absolute velocities
      rd_out(1) =                   - rd_in(6)*rRel(2) + rd_in(5)*rRel(3) + rd_in(1) ! x   
      rd_out(2) =  rd_in(6)*rRel(1)                    - rd_in(4)*rRel(3) + rd_in(2) ! y
      rd_out(3) = -rd_in(5)*rRel(1) + rd_in(4)*rRel(2)                    + rd_in(3) ! z
      
      ! absolute accelerations
      rd_out(1) =                   - rd_in(6)*rRel(2) + rd_in(5)*rRel(3) + rd_in(1) ! x   
      rd_out(2) =  rd_in(6)*rRel(1)                    - rd_in(4)*rRel(3) + rd_in(2) ! y
      rd_out(3) = -rd_in(5)*rRel(1) + rd_in(4)*rRel(2)                    + rd_in(3) ! z



      !rRel = MATMUL(TransMat, rRelBody)        
      !H = getH(rRel)      
      !! absolute locations
      !r_out = rRel + r_in 
      !! absolute velocities
      !rd_out = MATMUL( H, rd_in(4:6)) + rd_in(1:3)  
      

   END SUBROUTINE TransformKinematics
   !-----------------------------------------------------------------------
   
   

   ! calculate position, velocity, and acceleration of point based on its position relative to moving 6DOF body
   !-----------------------------------------------------------------------
   SUBROUTINE TransformKinematicsA(rRelBody, r_in, TransMat, v_in, a_in, r_out, v_out, a_out)
      REAL(DbKi),       INTENT(IN   )  :: rRelBody(:)  ! relative location of point about reference point, in local/reference coordinate system
      REAL(DbKi),       INTENT(IN   )  :: r_in(3)      ! translation applied to reference point
      REAL(DbKi),       INTENT(IN   )  :: TransMat(3,3)! rotation matrix describing rotation about reference point
      REAL(DbKi),       INTENT(IN   )  :: v_in(6)      ! 6DOF velecity vector about ref point in global orientation frame
      REAL(DbKi),       INTENT(IN   )  :: a_in(6)      ! 6DOF acceleration vector
      REAL(DbKi),       INTENT(  OUT)  :: r_out(3)     ! coordinates of point of interest
      REAL(DbKi),       INTENT(  OUT)  :: v_out(3)     ! velocity of point
      REAL(DbKi),       INTENT(  OUT)  :: a_out(3)     ! acceleration of point

      REAL(DbKi)                       :: rRel(3)
      REAL(DbKi)                       :: rRel2(3)

      REAL(DbKi)                       :: r_out2(3)
      REAL(DbKi)                       :: rd_out2(3)
      REAL(DbKi)                       :: H(3,3)

      ! rd_in should be in global orientation frame
      ! note: it's okay if r_out and rd_out are 6-size. Only the first 3 will be written, and 4-6 will
      !       already be correct or can be assigned seperately from r_in and rd_in (assuming orientation frames are identical)


      ! locations about ref point in *unrotated* reference frame
      !rRel2(1) = TransMat(1,1)*rRelBody(1) + TransMat(2,1)*rRelBody(2) + TransMat(3,1)*rRelBody(3)  ! x
      !rRel2(2) = TransMat(1,2)*rRelBody(1) + TransMat(2,2)*rRelBody(2) + TransMat(3,2)*rRelBody(3)  ! y
      !rRel2(3) = TransMat(1,3)*rRelBody(1) + TransMat(2,3)*rRelBody(2) + TransMat(3,3)*rRelBody(3)  ! z
      
      rRel = MATMUL(TransMat, rRelBody)  
      
      H = getH(rRel)
      
      ! absolute locations
      r_out = rRel + r_in 

      ! absolute velocities
      !rd_out2(1) =                  - v_in(6)*rRel(2) + v_in(5)*rRel(3) + v_in(1) ! x   
      !rd_out2(2) =  v_in(6)*rRel(1)                   - v_in(4)*rRel(3) + v_in(2) ! y
      !rd_out2(3) = -v_in(5)*rRel(1) + v_in(4)*rRel(2)                   + v_in(3) ! z
      
      v_out = MATMUL( H, v_in(4:6)) + v_in(1:3)  
      
      ! absolute accelerations
      a_out = MATMUL( H, a_in(4:6)) + a_in(1:3)   ! << should add second order terms!
        

   END SUBROUTINE TransformKinematicsA
   !-----------------------------------------------------------------------
   
   ! calculate position and velocity of point along rod (distance L along direction u)
   !-----------------------------------------------------------------------
   SUBROUTINE TransformKinematicsAtoB(rA, u, L, rd_in, r_out, rd_out)
      REAL(DbKi),       INTENT(IN   )  :: rA(3)        ! coordinate of end A
      REAL(DbKi),       INTENT(IN   )  :: u(3)         ! Rod unit vector
      REAL(DbKi),       INTENT(IN   )  :: L            ! Rod length from end A to B
      REAL(DbKi),       INTENT(IN   )  :: rd_in(6)     ! 6DOF velecity vector about Rod end A, in global orientation frame
      REAL(DbKi),       INTENT(  OUT)  :: r_out(3)     ! coordinates of end B
      REAL(DbKi),       INTENT(  OUT)  :: rd_out(3)    ! velocity of end B

      REAL(DbKi)                       :: rRel(3)
      
      
      ! locations (unrotated reference frame)
      rRel = L*u              ! relative location of point B from point A
      r_out = rRel + rA	        ! absolute location of point B
      
      ! absolute velocities
      rd_out(1) =                   - rd_in(6)*rRel(2) + rd_in(5)*rRel(3) + rd_in(1)  ! x   
      rd_out(2) =  rd_in(6)*rRel(1)                    - rd_in(4)*rRel(3) + rd_in(2)  ! y
      rd_out(3) = -rd_in(5)*rRel(1) + rd_in(4)*rRel(2)                    + rd_in(3)  ! z
      
		
   END SUBROUTINE TransformKinematicsAtoB
   !-----------------------------------------------------------------------
   
   !
   !-----------------------------------------------------------------------
   SUBROUTINE TranslateForce3to6DOF(dx, F, Fout)
      REAL(DbKi),       INTENT(IN   )  :: dx(3)       ! displacement vector from ref point to point of force (F) application
      REAL(DbKi),       INTENT(IN   )  :: F(3)        ! applied force
      REAL(DbKi),       INTENT(  OUT)  :: Fout(6)     ! resultant applied force and moment about ref point
        
      Fout(1:3) = F
      
      Fout(4:6) = CROSS_PRODUCT(dx, F)
		
   END SUBROUTINE TranslateForce3to6DOF
   !-----------------------------------------------------------------------
   
   
   !
   !-----------------------------------------------------------------------
   SUBROUTINE TranslateMass3to6DOF(dx, Min, Mout)
      REAL(DbKi),       INTENT(IN   )  :: dx(3)       ! displacement vector from ref point to point of mass matrix (Min)
      REAL(DbKi),       INTENT(IN   )  :: Min( 3,3)   ! original mass matrix (assumed at center of mass, or a point mass)
      REAL(DbKi),       INTENT(  OUT)  :: Mout(6,6)   ! resultant mass and inertia matrix about ref point
  
      REAL(DbKi)                       :: H(     3,3) ! "anti-symmetric tensor components" from Sadeghi and Incecik
      REAL(DbKi)                       :: tempM( 3,3)
      REAL(DbKi)                       :: tempM2(3,3)
      REAL(DbKi)                       :: Htrans(3,3)
      Integer(IntKi)                   :: I,J
   
      ! sub-matrix definitions are accordint to  | m    J |
      !                                          | J^T  I |
   
      H = getH(dx);
		
      ! mass matrix  [m'] = [m]
      Mout(1:3,1:3) = Min
		
      ! product of inertia matrix  [J'] = [m][H] + [J]
      Mout(1:3,4:6) = MATMUL(Min, H)
      Mout(4:6,1:3) = TRANSPOSE(Mout(1:3,4:6)) 
	
      !moment of inertia matrix  [I'] = [H][m][H]^T + [J]^T [H] + [H]^T [J] + [I]
      Mout(4:6,4:6) = MATMUL(MATMUL(H, Min), TRANSPOSE(H))
		
   END SUBROUTINE TranslateMass3to6DOF
   !-----------------------------------------------------------------------
   
   !
   !-----------------------------------------------------------------------
   SUBROUTINE TranslateMass6to6DOF(dx, Min, Mout)
      REAL(DbKi),       INTENT(IN   )  :: dx(3)       ! displacement vector from ref point to point of mass matrix (Min)
      REAL(DbKi),       INTENT(IN   )  :: Min( 6,6)   ! original mass matrix 
      REAL(DbKi),       INTENT(  OUT)  :: Mout(6,6)   ! resultant mass and inertia matrix about ref point
  
      REAL(DbKi)                       :: H(     3,3) ! "anti-symmetric tensor components" from Sadeghi and Incecik
         
      H = getH(dx);
		
      ! mass matrix  [m'] = [m]
      Mout(1:3,1:3) = Min(1:3,1:3)
		
      ! product of inertia matrix  [J'] = [m][H] + [J]
      Mout(1:3,4:6) = MATMUL(Min(1:3,1:3), H) + Min(1:3,4:6)
      Mout(4:6,1:3) = TRANSPOSE(Mout(1:3,4:6))      
	
      !moment of inertia matrix  [I'] = [H][m][H]^T + [J]^T [H] + [H]^T [J] + [I]
      Mout(4:6,4:6) = MATMUL(MATMUL(H, Min(1:3,1:3)), TRANSPOSE(H)) + MATMUL(Min(4:6,1:3),H) + MATMUL(TRANSPOSE(H),Min(1:3,4:6)) + Min(4:6,4:6)
		
   END SUBROUTINE TranslateMass6to6DOF
   !-----------------------------------------------------------------------
   
   ! produce alternator matrix
   !-----------------------------------------------------------------------
   FUNCTION GetH(r)
      Real(DbKi), INTENT(IN)      :: r(3)     ! inputted vector
      Real(DbKi)                  :: GetH(3,3) ! outputted matrix
      
      GetH(2,1) = -r(3)
      GetH(1,2) =  r(3)
      GetH(3,1) =  r(2)
      GetH(1,3) = -r(2)
      GetH(3,2) = -r(1)
      GetH(2,3) =  r(1)
      
      GetH(1,1) = 0.0_DbKi
      GetH(2,2) = 0.0_DbKi
      GetH(3,3) = 0.0_DbKi
	
   END FUNCTION GetH
   !-----------------------------------------------------------------------
   
   
   
   ! apply a rotation to a 6-by-6 mass/inertia tensor (see Sadeghi and Incecik 2005 for theory)
   !-----------------------------------------------------------------------
   FUNCTION RotateM6(Min, rotMat) result(outMat)
      
      Real(DbKi), INTENT(IN)      :: Min(6,6)     ! inputted matrix to be rotated
      Real(DbKi), INTENT(IN)      :: rotMat(3,3)  ! rotation matrix (DCM)
      Real(DbKi)                  :: outMat(6,6)  ! rotated matrix
      
      Real(DbKi)                  :: tempM(3,3)  
      Real(DbKi)                  :: tempMrot(3,3)  

      ! the process for each of the following is to 
      ! 1. copy out the relevant 3x3 matrix section,
      ! 2. rotate it, and
      ! 3. paste it into the output 6x6 matrix
		
      
      ! mass matrix
      outMat(1:3,1:3) = rotateM3(Min(1:3,1:3), rotMat)

      ! product of inertia matrix	
      outMat(1:3,4:6) = rotateM3(Min(1:3,4:6), rotMat)      
      outMat(4:6,1:3) = TRANSPOSE(outMat(1:3,4:6))

      ! moment of inertia matrix
      outMat(4:6,4:6) = rotateM3(Min(4:6,4:6), rotMat)
   
   END FUNCTION RotateM6


   ! apply a rotation to a 3-by-3 mass matrix or any other second order tensor
   !-----------------------------------------------------------------------
   FUNCTION RotateM3(Min, rotMat) result(outMat)
      
      Real(DbKi), INTENT(IN)      :: Min(3,3)     ! inputted matrix to be rotated
      Real(DbKi), INTENT(IN)      :: rotMat(3,3)  ! rotation matrix (DCM)
      Real(DbKi)                  :: outMat(3,3)  ! rotated matrix
   
      ! overall operation is [m'] = [a]*[m]*[a]^T
	
      outMat = MATMUL( MATMUL(rotMat, Min), TRANSPOSE(rotMat) )
	
   END FUNCTION RotateM3
   
   
   
   
   
   ! calculates rotation matrix R to rotate from global axes to a member's local axes
   !-----------------------------------------------------------------------
   FUNCTION CalcOrientation(phi, beta, gamma) result(R)
      
      REAL(DbKi),      INTENT ( IN    )  :: phi     ! member incline angle
      REAL(DbKi),      INTENT ( IN    )  :: beta    ! member incline heading
      REAL(DbKi),      INTENT ( IN    )  :: gamma   ! member twist angle
      REAL(DbKi)                         :: R(3,3)  ! rotation matrix 
     
      INTEGER(IntKi)                     :: errStat  
      CHARACTER(100)                     :: errMsg   

      REAL(DbKi)                         :: s1, c1, s2, c2, s3, c3


      ! trig terms for Euler angles rotation based on beta, phi, and gamma
      s1 = sin(beta) 
      c1 = cos(beta)
      s2 = sin(phi) 
      c2 = cos(phi)
      s3 = sin(gamma) 
      c3 = cos(gamma)
      
      ! calculate rotation matrix based on Z1Y2Z3 Euler rotation sequence from https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix

      R(1,1) = c1*c2*c3-s1*s3
      R(1,2) = -c3*s1-c1*c2*s3
      R(1,3) = c1*s2
      R(2,1) = c1*s3+c2*c3*s1
      R(2,2) = c1*c3-c2*s1*s3
      R(2,3) = s1*s2
      R(3,1) = -c3*s2
      R(3,2) = s2*s3
      R(3,3) = c2

      ! could also calculate unit normals p1 and p2 for rectangular cross sections   
      !p1 = matmul( R, [1,0,0] )               ! unit vector that is perpendicular to the 'beta' plane if gamma is zero
      !p2 = cross( q, p1 )                     ! unit vector orthogonal to both p1 and q
      
   END FUNCTION CalcOrientation
   

   !compute the inverse of a 3-by-3 matrix m
   !-----------------------------------------------------------------------
   SUBROUTINE Inverse3by3( Minv, M )
      Real(DbKi), INTENT(OUT)   :: Minv(3,3)  ! returned inverse matrix
      Real(DbKi), INTENT(IN)    :: M(3,3)     ! inputted matrix

      Real(DbKi)                :: det        ! the determinant
      Real(DbKi)                :: invdet     ! inverse of the determinant

      det = M(1, 1) * (M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3)) - &
            M(1, 2) * (M(2, 1) * M(3, 3) - M(2, 3) * M(3, 1)) + &
            M(1, 3) * (M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1));

      invdet = 1.0 / det   ! because multiplying is faster than dividing

      Minv(1, 1) = (M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3)) * invdet
      Minv(1, 2) = (M(1, 3) * M(3, 2) - M(1, 2) * M(3, 3)) * invdet
      Minv(1, 3) = (M(1, 2) * M(2, 3) - M(1, 3) * M(2, 2)) * invdet
      Minv(2, 1) = (M(2, 3) * M(3, 1) - M(2, 1) * M(3, 3)) * invdet
      Minv(2, 2) = (M(1, 1) * M(3, 3) - M(1, 3) * M(3, 1)) * invdet
      Minv(2, 3) = (M(2, 1) * M(1, 3) - M(1, 1) * M(2, 3)) * invdet
      Minv(3, 1) = (M(2, 1) * M(3, 2) - M(3, 1) * M(2, 2)) * invdet
      Minv(3, 2) = (M(3, 1) * M(1, 2) - M(1, 1) * M(3, 2)) * invdet
      Minv(3, 3) = (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2)) * invdet

   END SUBROUTINE Inverse3by3
   !-----------------------------------------------------------------------


   ! One-function implementation of Crout LU Decomposition. Solves Ax=b for x
   SUBROUTINE LUsolve(n, A, LU, b, y, x)

      INTEGER(intKi),   INTENT(IN   )  :: n           ! size of matrices and vectors
      Real(DbKi),       INTENT(IN   )  :: A( n,n)     ! LHS matrix (e.g. mass matrix)
      Real(DbKi),       INTENT(INOUT)  :: LU(n,n)     ! stores LU matrix data
      Real(DbKi),       INTENT(IN   )  :: b(n)        ! RHS vector
      Real(DbKi),       INTENT(INOUT)  :: y(n)        ! temporary vector
      Real(DbKi),       INTENT(  OUT)  :: x(n)        ! LHS vector to solve for

      INTEGER(intKi)                   :: i,j,k,p
      Real(DbKi)                       :: sum
	
      DO k = 1,n
         DO i = k,n
         
            sum = 0.0_DbKi
            
            DO p=1,k-1   !for(int p=0; p<k; ++p)
               sum = sum + LU(i,p)*LU(p,k)
            END DO
            
            LU(i,k) = A(i,k) - sum
         END DO !I
         
         DO j=k+1,n  !for(int j=k+1;j<n;++j)
            
            sum = 0.0_DbKi

            DO p=1,k-1   !for(int p=0;p<k;++p)
               sum = sum + LU(k,p)*LU(p,j)
            END DO
            
            LU(k,j) = (A(k,j)-sum)/LU(k,k)
         END DO !j
         
      END DO !K
      
      DO i=1,n
      
         sum = 0.0_DbKi
         
         DO k=1,i-1  !for(int k=0; k<i; ++k)
            sum = sum + LU(i,k)*y(k);
         END DO
         
         y(i) = (b(i)-sum)/LU(i,i)
         
      END DO
      
      DO j=1,n       ! this is actually for looping through i in reverse
         i = n+1-j   ! for(int i=n-1; i>=0; --i)
      
         sum = 0.0_DbKi
         
         DO k=i+1, n 
            sum = sum + LU(i,k)*x(k)
         END DO
         
         x(i) = (y(i)-sum) 
         
      END DO !j (actually decrementing i)
      
   END SUBROUTINE LUsolve




END MODULE MoorDyn
