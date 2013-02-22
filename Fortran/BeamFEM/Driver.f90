PROGRAM main
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init
      
      ! LOCAL VARIABLES
   INTEGER       :: I, J
   
   CALL SubDyn_Input(Init)

         !************** debug *********************************************
         write(*,*) 'RootName : ', trim(Init%RootName)
         write(*,*) 'NDiv = ', Init%NDiv
         write(*,*) 'CBMod = ', Init%CBMod
         write(*,*) ('Nmodes = ', Init%Nmodes)
         write(*,*) ('Damping = ', Init%JDampings(1:Init%Nmodes))
         write(*,*) ('NJoints = ', Init%NJoints)
         write(*, *) (Init%Joints(i, 1:Init%JointsCol), i=1,Init%NJoints)
         write(*, *) ('NReact = ', Init%NReact)
         write(*, *) (Init%Reacts(i, 1:Init%ReactCol), i = 1,Init%NReact)
         write(*, *) ('NInterf = ', Init%NInterf)
         write(*, *) (Init%Interf(i, 1:Init%InterfCol), i=1,Init%NInterf)
         write(*, *) ('NMembers = ', Init%NMembers)
         write(*, *) (Init%Members(i, 1:Init%MembersCol), i=1,Init%NMembers)
         write(*, *) ('NPropSets = ', Init%NPropSets)
         write(*, *) (Init%PropSets(i, 1:Init%PropSetsCol), i=1,Init%NPropSets)
         write(*, *) ('NXPropSets = ', Init%NXPropSets)
         write(*, *) (Init%XPropSets(i, 1:Init%XPropSetsCol), i=1,Init%NXPropSets)
         write(*, *) ('NCOSMs = ', Init%NCOSMs)
         write(*, *) (Init%COSMs(i, 1:Init%COSMsCol), i=1,Init%NCOSMs)
         write(*, *) ('NCMass = ', Init%NCMass)
         write(*, *) (Init%CMass(i, 1:Init%CMassCol), i=1,Init%NCMass)

!-------------  Discretize the structure according to the division size -----------------
  CALL SubDyn_Discrt(Init)
   

      ! Initialize index array for system K and M
  CALL InitIAJA(Init)

      ! Assemble system stiffness and mass matrices
  CALL AssembleKM(Init)

      ! Assemble system force vector
!  CALL AssembleF()      

      ! Apply constraints to stiffness and mass matrices
  CALL ApplyConstr(Init)

      ! Solve static problem
!  CALL StaticSolve()

      ! Solve dynamic problem
  CALL EigenSolve(Init)



   

END PROGRAM main
!----------------------------------------------------------------------------
SUBROUTINE SubDyn_Input(Init)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init

! loacal variable for input and output

CHARACTER( 1024)             :: Comment                                         ! String to temporarily hold the comment line.
CHARACTER(   3)              :: EndOfFile                                       ! String read in at the end of the input file.
CHARACTER(  35)              :: Frmt      = "( 2X, L11, 2X, A, T30, ' - ', A )" ! Output format for logical parameters. (matches NWTC Subroutine Library format)
CHARACTER(1000)              :: OutLine                                         ! String to temporarily hold the output parameter list.
CHARACTER(1024)              :: PriPath                                         ! The path to the primary input file
CHARACTER(1024)              :: FTitle                                          ! The title line from the primary input file.
INTEGER(4)                   :: IOS 
INTEGER(4)                   :: Sttus
LOGICAL                      :: WrEcho
INTEGER(IntKi)               :: UnIn
INTEGER(IntKi)               :: UnOut 
INTEGER(IntKi)               :: ErrStat



INTEGER(IntKi)               :: I

Init%SDInputFile = 'testframe.txt'


CALL GetNewUnit( UnIn )   
CALL OpenFInpfile(UnIn, TRIM(Init%SDInputFile), ErrStat)

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

CALL GetPath( Init%SDInputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.
CALL GetRoot( Init%SDInputFile, Init%RootName )


!-------------------------- HEADER ---------------------------------------------

   CALL ReadCom( UnIn, Init%SDInputFile, 'SubDyn input file header line 1', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


   CALL ReadCom( UnIn, Init%SDInputFile, 'SubDyn input file header line 2', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   CALL ReadCom( UnIn, Init%SDInputFile, 'SubDyn input file header line 3', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

!-------------------------- SIMULATION CONTROL PARAMETERS ----------------------

      ! Skip the comment line.

   CALL ReadCom( UnIn, Init%SDInputFile, ' SIMULATION CONTROL PARAMETERS ', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   ! Echo - Echo input to "echo.out".

READ (UnIn,*,IOSTAT=IOS)  WrEcho
CALL CheckIOS ( IOS, Init%SDInputFile, 'Echo', FlgType )
Echo = WrEcho

IF ( Echo )  THEN
   CALL OpenEcho ( UnEc, TRIM(Init%RootName)//'.ech' )
   WRITE (UnEc,'(/,A)'   )  'This file of echoed input was generated by '//TRIM(ProgName)//' '//TRIM(ProgVer)// &
                            ' on '//CurDate()//' at '//CurTime()//'.'
   WRITE (UnEc,'(/,A,/)' )  'Substructure data from file "'//TRIM( Init%SDInputFile )//'":'
   WRITE (UnEc,Frmt      )  Echo, 'Echo', 'Echo input to "echo.out"'
ENDIF

!-------------------- FEA and CRAIG-BAMPTON PARAMETERS---------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' FEA and CRAIG-BAMPTON PARAMETERS ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! FEMMod - FEM switch: element model in the FEM: 0= Euler-Bernoulli(E-B) ; 1=Tapered E-B; 2= 2-node Timoshenko;  3= 2-node tapered Timoshenko

CALL ReadIVar ( UnIn, Init%SDInputFile, Init%FEMMod, 'FEMMod', 'FEM analysis mode' )

IF ( ( Init%FEMMod < 0 ) .OR. ( Init%FEMMod > 4 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': FEMMod must be 0, 1, 2, or 3.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF


   ! NDiv - Number sub-elements per member

CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NDiv, 'NDiv', 'Number of divisions per member' )

IF ( ( Init%NDiv < 1 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': NDiv must be a positive integer', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF



   ! CBMod - Perform C-B flag.
!READ( UnIn, *, IOSTAT=ErrStat ) Init%CBMod
CALL ReadLVar ( UnIn, Init%SDInputFile, Init%CBMod, 'CBMod', 'C-B mod flag' )

IF ( ErrStat /= 0 ) THEN
   CALL CheckIOS ( ErrStat, Init%SDInputFile, 'CBMod', NumType, .TRUE. )
   CLOSE( UnIn )
   RETURN
END IF



IF (Init%CBMod) THEN

      ! Nmodes - Number of interal modes to retain.
   CALL ReadIVar ( UnIn, Init%SDInputFile, Init%Nmodes, 'Nmodes', 'Number of internal modes' )

   IF ( ( Init%Nmodes < 1 ) )  THEN
      CALL ProgAbort ( ' FEMMod must be a positive integer.' )
   ENDIF
   
      ! Damping ratios for retained modes
   ALLOCATE(Init%JDampings(Init%Nmodes), STAT=Sttus)
   
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the damping ratio array.' )
   ENDIF

   CALL ReadRAry( UnIn, Init%SDInputFile, Init%JDampings, Init%Nmodes, 'JDamping', 'Damping ratio of the internal modes', ErrStat )
      
   DO I = 1, Init%Nmodes
      IF ( ( Init%JDampings(I) .LE. 0 ) .OR.( Init%JDampings(I) .GE. Init%Nmodes ) ) THEN
         WRITE(Comment, *) Init%NModes
         CALL ProgAbort ( ' Number of damping ratio should be larger than 0 and less than '//trim(Comment)//'.' )
      ENDIF
      
      Init%JDampings(I) = Init%JDampings(I)/100.0
      
   ENDDO
   
ELSE
      ! skip 2 lines
   READ (UnIn,'(A)',IOSTAT=IOS)  Comment
   CALL CheckIOS( IOS, Init%SDInputFile, 'Nmodes ', StrType )
   READ (UnIn,'(A)',IOSTAT=IOS)  Comment
   CALL CheckIOS( IOS, Init%SDInputFile, 'Damping ratio', StrType )
      
   Init%Nmodes = 0
ENDIF

!---- STRUCTURE JOINTS: joints connect structure members -----------------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Joints connect structure members ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of joints
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NJoints, 'NJoints', 'Number of joints' )
IF ( ( Init%NJoints < 2 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': NJoints must be greater than 1', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip two lines
CALL ReadCom( UnIn, Init%SDInputFile, ' Joints description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

CALL ReadCom( UnIn, Init%SDInputFile, ' units ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
   
   ! Joints coordinates
ALLOCATE(Init%Joints(Init%NJoints, Init%JointsCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating Joints arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

DO I = 1, Init%NJoints

   CALL ReadRAry( UnIn, Init%SDInputFile, Init%Joints(I,1:Init%JointsCol), Init%JointsCol, 'Joints', 'Joint number and coordinates', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   
ENDDO

!------------------- BASE REACTION JOINTS: T/F for Locked/Free DOF @ each Reaction Node ---------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' BASE REACTION JOINTS ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! Number of reaction joints (The joints should be all fixes for now)
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NReact, 'NReact', 'Number of joints with reaction forces' )
IF ( ( Init%NReact < 1 ) .OR. (Init%NReact > Init%NJoints) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': NReact must be greater than 0 and less than number of joints', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF
   
   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' BASE REACTION JOINTS, description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF


   ! Joints with reaction forces, joint number and locked/free dof
ALLOCATE(Init%Reacts(Init%NReact, Init%ReactCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating Reacts arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF
   
   
DO I = 1, Init%NReact

   CALL ReadAry( UnIn, Init%SDInputFile, Init%Reacts(I,1:Init%ReactCol), Init%ReactCol, 'Reacts', 'Joint number and dof', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   
ENDDO


!------- INTERFACE JOINTS: T/F for Locked (to the TP)/Free DOF @each Interface Joint (only Locked-to-TP implemented thus far (=rigid TP)) ---------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' INTERFACE JOINTS ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! Number of interface joints (The joints should be all fixes for now)
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NInterf, 'NInterf', 'Number of joints fixed to TP' )
IF ( ( Init%NInterf < 0 ).OR. (Init%NInterf > Init%NJoints)  )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': NInterf must be non-negative and less than number of joints.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' INTERFACE JOINTS, description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF


   ! Joints with reaction forces, joint number and locked/free dof
ALLOCATE(Init%Interf(Init%NInterf, Init%InterfCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating Interf arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
  
ENDIF
   
   
DO I = 1, Init%NInterf

   CALL ReadAry( UnIn, Init%SDInputFile, Init%Interf(I,1:Init%InterfCol), Init%InterfCol, 'Interf', 'Interface joint number and dof', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   
ENDDO
   
!----------------------------------- MEMBERS --------------------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Members ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of members
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NMembers, 'NMembers', 'Number of members' )
IF ( ( Init%NMembers < 1 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': NMembers must be greater than 0', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip one line
CALL ReadCom( UnIn, Init%SDInputFile, ' Members description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
   
   ! Member connection
ALLOCATE(Init%Members(Init%NMembers, Init%MembersCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating Members arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF


DO I = 1, Init%NMembers

   CALL ReadAry( UnIn, Init%SDInputFile, Init%Members(I,1:Init%MembersCol), Init%MembersCol, 'Members', 'Member number and connectivity ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   
ENDDO   

!------------------ MEMBER X-SECTION PROPERTY data 1/2 [isotropic material for now: use this table if circular-tubular elements ------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Member X-Section property 1/2 ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of property sets
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NPropSets, 'NPropSets', 'Number of property sets' )
IF ( ( Init%NPropSets < 1 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': NPropSets must be greater than 0', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip two lines
CALL ReadCom( UnIn, Init%SDInputFile, ' Property sets 1/2 description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
CALL ReadCom( UnIn, Init%SDInputFile, ' Property sets 1/2 units ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
   
   ! Property sets value
ALLOCATE(Init%PropSets(Init%NPropSets, Init%PropSetsCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating PropSets arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF


DO I = 1, Init%NPropSets

   CALL ReadRAry( UnIn, Init%SDInputFile, Init%PropSets(I,1:Init%PropSetsCol), Init%PropSetsCol, 'PropSets', 'PropSets number and values ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   
ENDDO   

!------------------ MEMBER X-SECTION PROPERTY data 2/2 [isotropic material for now: use this table if any section other than circular, however provide COSM(i,j) below) ------------------------


   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Member X-Section property 2/2 ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of property sets
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NXPropSets, 'NPropSets', 'Number of property sets' )
IF ( ( Init%NPropSets < 0 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': NPropSets must be non-negative.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip two lines
CALL ReadCom( UnIn, Init%SDInputFile, ' Property sets 2/2 description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
CALL ReadCom( UnIn, Init%SDInputFile, ' Property sets 2/2 units ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
   
   ! Property sets value
ALLOCATE(Init%XPropSets(Init%NXPropSets, Init%XPropSetsCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating PropSets arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF


DO I = 1, Init%NXPropSets

   CALL ReadRAry( UnIn, Init%SDInputFile, Init%XPropSets(I,1:Init%XPropSetsCol), Init%XPropSetsCol, 'PropSets', 'PropSets number and values ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   
ENDDO   

!---------------------- MEMBER COSINE MATRICES COSM(i,j) ------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Member direction cosine matrices ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of direction cosine matrices
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NCOSMs, 'NCOSMs', 'Number of unique direction cosine matrices' )
IF ( ( Init%NPropSets < 0 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': NCOSMs must be non-negative.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip one line
CALL ReadCom( UnIn, Init%SDInputFile, ' Property sets 2/2 description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! Direction cosine matrices value
ALLOCATE(Init%COSMs(Init%NCOSMs, Init%COSMsCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating COSMs arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF


DO I = 1, Init%NCOSMs

   CALL ReadRAry( UnIn, Init%SDInputFile, Init%COSMs(I,1:Init%COSMsCol), Init%COSMsCol, 'PropSets', 'PropSets number and values ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   
ENDDO   

!------------------------ JOINT ADDITIONAL CONCENTRATED MASSES--------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, Init%SDInputFile, ' Additional concentrated masses at joints ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of joints that have concentrated masses
CALL ReadIVar ( UnIn, Init%SDInputFile, Init%NCMass, 'NCMass', 'Number of joints that have concentrated masses' )
IF ( ( Init%NPropSets < 0 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': NCMass must be non-negative.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip two lines
CALL ReadCom( UnIn, Init%SDInputFile, ' Concentrated mass description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

CALL ReadCom( UnIn, Init%SDInputFile, ' Concentrated mass units ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! Concentrated mass value
ALLOCATE(Init%CMass(Init%NCMass, Init%CMassCol), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(Init%SDInputFile)//': Error allocating CMass arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF


DO I = 1, Init%NCMass

   CALL ReadRAry( UnIn, Init%SDInputFile, Init%CMass(I,1:Init%CMassCol), Init%CMassCol, 'CMass', 'Joint number and mass values ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   
ENDDO   

!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------
!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------
!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------


!         to be finished.......


!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------
!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------
!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------



CLOSE( UnIn )


END SUBROUTINE SubDyn_Input
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE SubDyn_Discrt(Init)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init

      ! local variable
   INTEGER                       :: I, J, Node1, Node2, Prop1, Prop2, flg, flg1, flg2
   INTEGER                       :: OldJointIndex(Init%NJoints)
   INTEGER                       :: NPE      ! node per element
   INTEGER(4)                    :: Sttus
   INTEGER                       :: TempNProp
   REAL(ReKi), ALLOCATABLE       :: TempProps(:, :)
   INTEGER, ALLOCATABLE          :: TempMembers(:, :)        
   INTEGER                       :: UnOut     
   INTEGER(IntKi)                :: ErrStat
   CHARACTER(1024)               :: OutFile
   CHARACTER(  50)               :: tempStr ! string number of nodes in member
   CHARACTER(1024)               :: OutFmt
   INTEGER                       :: knode, kelem, kprop, nprop
   REAL(ReKi)                    :: x1, y1, z1, x2, y2, z2, dx, dy, dz, dd, dt, d1, d2, t1, t2
   
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   ! number of nodes per element
   IF( ( Init%FEMMod .GE. 0 ) .and. (Init%FEMMod .LE. 3) ) THEN
      NPE = 2 
   ENDIF
   
   ! Calculate total number of nodes according to divisions
   Init%NNode = Init%NJoints + ( Init%NDiv - 1 )*Init%NMembers
   ! Total number of element
   Init%NElem = Init%NMembers*Init%NDiv
   ! Total number of property sets (temp)
   TempNProp = Init%NElem*NPE
   
   ! Calculate total number of nodes and elements according to element types
   ! for 3-node or 4-node beam elements
   Init%NNode = Init%NNode + (NPE - 2)*Init%NElem
   Init%MembersCol = Init%MembersCol + (NPE - 2) 
   
   ! check the number of interior modes
   IF ( Init%NModes .GT. 6*(Init%NNode - Init%NInterf - Init%NReact) ) THEN
      WRITE(tempStr, *) 6*(Init%NNode - Init%NInterf - Init%NReact)
      CALL ProgAbort ( ' The NModes must be less or equal to'//TRIM(tempStr), .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   
   
   
   ! Allocate for nodes and elements and membernodes
   ALLOCATE(Init%Nodes(Init%NNode, Init%JointsCol), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating Nodes arrays', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%Nodes = 0
   
   ALLOCATE(Init%Elems(Init%NElem, Init%MembersCol), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating Elems arrays', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%Elems = 0

   ! for two-node element only, otherwise the number of nodes in one element is different
   ALLOCATE(Init%MemberNodes(Init%NMembers, Init%NDiv+1), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating MemberNodes arrays', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%MemberNodes = 0

      ! Allocate Temp members and property sets
   ALLOCATE(TempMembers(Init%NMembers, Init%MembersCol), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating TempProps arrays', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   TempMembers = 0

   ALLOCATE(TempProps(TempNProp, Init%PropSetsCol), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating TempProps arrays', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   TempProps = 0

   ! Initialize Nodes
   DO I = 1,Init%NJoints
      OldJointIndex(I) = Init%Joints(I, 1)
      Init%Nodes(I, 1) = I
      Init%Nodes(I, 2) = Init%Joints(I, 2)
      Init%Nodes(I, 3) = Init%Joints(I, 3)
      Init%Nodes(I, 4) = Init%Joints(I, 4)
   ENDDO
   
   ! Initialize TempMembers and Elems
   DO I = 1, Init%NMembers
      TempMembers(I, 1) = I
      TempMembers(I, 4) = Init%Members(I, 4)
      TempMembers(I, 5) = Init%Members(I, 5)
      
      Init%Elems(I,     1) = I
      Init%Elems(I, NPE+2) = Init%Members(I, 4)
      Init%Elems(I, NPE+3) = Init%Members(I, 5)
      
      Node1 = Init%Members(I, 2)
      Node2 = Init%Members(I, 3)
      
      flg1 = 0
      flg2 = 0
      DO J = 1, Init%NJoints
         IF ( Node1 == Init%Joints(J, 1) ) THEN
            TempMembers(I, 2) = J
            Init%Elems(I, 2) = J
            flg1 = 1
         ENDIF
         IF ( Node2 == Init%Joints(J, 1) ) THEN
            TempMembers(I, 3) = J
            Init%Elems(I, NPE+1) = J
            flg2 = 1
         ENDIF
         
      ENDDO
      
      IF ( (flg1 == 0) .OR. (flg2 == 0) ) THEN
         CALL ProgAbort ( ' Member has node not in the node list !', .TRUE. )
         ErrStat = 1
         RETURN
      ENDIF

   ENDDO
   
   ! Initialize Temp property set
   TempProps(1:Init%NPropSets, :) = Init%PropSets(1:Init%NPropSets, :)   
   
   ! Initialize boundary constraint vector
   ! Change the node number
   ALLOCATE(Init%BCs(6*Init%NReact, 2), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating TempProps arrays', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%BCs = 0
      
   DO I = 1, Init%NReact
      Node1 = Init%Reacts(I, 1);
      flg = 0
      DO J = 1, Init%NJoints
         IF ( Node1 == Init%Joints(J, 1) ) THEN
            Node2 = J
            flg = 1
         ENDIF
      ENDDO
      
      IF (flg == 0) THEN
         CALL ProgAbort ( ' Interf has node not in the node list !', .TRUE. )
         ErrStat = 1
         RETURN
      ENDIF

      
      DO J = 1, 6
         Init%BCs( (I-1)*6+J, 1) = (Node2-1)*6+J;
         Init%BCs( (I-1)*6+J, 2) = Init%Reacts(I, J+1);
      ENDDO
      
   ENDDO
   
      
   ! Initialize interface constraint vector
   ! Change the node number
   ALLOCATE(Init%IntFc(6*Init%NInterf, 2), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating TempProps arrays', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%IntFc = 0
      
   DO I = 1, Init%NInterf
      Node1 = Init%Interf(I, 1);
      flg = 0
      DO J = 1, Init%NJoints
         IF ( Node1 == Init%Joints(J, 1) ) THEN
            Node2 = J
            flg = 1
         ENDIF
      ENDDO
      
      IF (flg == 0) THEN
         CALL ProgAbort ( ' Interf has node not in the node list !', .TRUE. )
         ErrStat = 1
         RETURN
      ENDIF
      
      DO J = 1, 6
         Init%IntFc( (I-1)*6+J, 1) = (Node2-1)*6+J;
         Init%IntFc( (I-1)*6+J, 2) = Init%Interf(I, J+1);
      ENDDO
   ENDDO
  
   ! Change numbering in concentrated mass matrix
   DO I = 1, Init%NCMass
      Node1 = Init%CMass(I, 1)
      DO J = 1, Init%NJoints
         IF ( Node1 == Init%Joints(J, 1) ) THEN
            Init%CMass(I, 1) = J
         ENDIF
      ENDDO
   ENDDO

! discretize structure according to NDiv 

knode = Init%NJoints
kelem = 0
kprop = Init%NPropSets

IF (Init%NDiv .GT. 1) THEN
   DO I = 1, Init%NMembers
      ! create new node
      Node1 = TempMembers(I, 2)
      Node2 = TempMembers(I, 3)
      
      IF ( Node1==Node2 ) THEN
         CALL ProgAbort ( ' Same starting and ending node in the member.', .TRUE. )
         ErrStat = 4
         RETURN
      ENDIF
    
      
      
      Prop1 = TempMembers(I, 4)
      Prop2 = TempMembers(I, 5)
      
      Init%MemberNodes(I,           1) = Node1
      Init%MemberNodes(I, Init%NDiv+1) = Node2
      
      IF  ( ( .not. EqualRealNos(TempProps(Prop1, 2),TempProps(Prop2, 2) ) ) &
       .OR. ( .not. EqualRealNos(TempProps(Prop1, 3),TempProps(Prop2, 3) ) ) &
       .OR. ( .not. EqualRealNos(TempProps(Prop1, 4),TempProps(Prop2, 4) ) ) )  THEN
      
         CALL ProgAbort ( ' Material E,G and rho in a member must be the same', .TRUE. )
         ErrStat = 3
         RETURN
      ENDIF

      x1 = Init%Nodes(Node1, 2)
      y1 = Init%Nodes(Node1, 3)
      z1 = Init%Nodes(Node1, 4)

      x2 = Init%Nodes(Node2, 2)
      y2 = Init%Nodes(Node2, 3)
      z2 = Init%Nodes(Node2, 4)
      
      dx = ( x2 - x1 )/Init%NDiv
      dy = ( y2 - y1 )/Init%NDiv
      dz = ( z2 - z1 )/Init%NDiv
      
      d1 = TempProps(Prop1, 5)
      t1 = TempProps(Prop1, 6)

      d2 = TempProps(Prop2, 5)
      t2 = TempProps(Prop2, 6)
      
      dd = ( d2 - d1 )/Init%NDiv
      dt = ( t2 - t1 )/Init%NDiv
      
      
      ! node connect to Node1
      knode = knode + 1
      Init%MemberNodes(I, 2) = knode
      CALL GetNewNode(knode, x1+dx, y1+dy, z1+dz, Init)
      
      IF ( ( .NOT.(EqualRealNos( dd , 0.0 ) ) ) .OR. &
           ( .NOT.( EqualRealNos( dt , 0.0 ) ) ) ) THEN   
           ! create a new property set 
           ! k, E, G, rho, d, t, Init
           
           kprop = kprop + 1
           CALL GetNewProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3),&
                           TempProps(Prop1, 4), d1+dd, t1+dt, TempProps, TempNProp, Init%PropSetsCol)           
           kelem = kelem + 1
           CALL GetNewElem(kelem, Node1, knode, Prop1, kprop, Init)  
           nprop = kprop              
      ELSE
           kelem = kelem + 1
           CALL GetNewElem(kelem, Node1, knode, Prop1, Prop1, Init)                
           nprop = Prop1 
      ENDIF
      
      ! interior nodes
      
      DO J = 2, (Init%NDiv-1)
         knode = knode + 1
         Init%MemberNodes(I, J+1) = knode

         CALL GetNewNode(knode, x1 + J*dx, y1 + J*dy, z1 + J*dz, Init)
         
         IF ( ( .NOT.(EqualRealNos( dd , 0.0 ) ) ) .OR. &
              ( .NOT.( EqualRealNos( dt , 0.0 ) ) ) ) THEN   
              ! create a new property set 
              ! k, E, G, rho, d, t, Init
              
              kprop = kprop + 1
              CALL GetNewProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3),&
                              Init%PropSets(Prop1, 4), d1 + J*dd, t1 + J*dt, &
                              TempProps, TempNProp, Init%PropSetsCol)           
              kelem = kelem + 1
              CALL GetNewElem(kelem, knode-1, knode, nprop, kprop, Init)
              nprop = kprop                
         ELSE
              kelem = kelem + 1
              CALL GetNewElem(kelem, knode-1, knode, nprop, nprop, Init)                
               
         ENDIF
      ENDDO
      
      ! the element connect to Node2
      kelem = kelem + 1
      CALL GetNewElem(kelem, knode, Node2, nprop, Prop2, Init)                

   ENDDO ! loop over all members

ELSE ! NDiv = 1

   Init%MemberNodes(1:Init%NMembers, 1:2) = Init%Elems(1:Init%NElem, 2:3)   

ENDIF ! if NDiv is greater than 1

! set the props in Init
ALLOCATE(Init%Props(kprop, Init%PropSetsCol), STAT=Sttus)
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error allocating TempProps arrays', .TRUE. )
   ErrStat = 1
   RETURN
ENDIF
Init%NProp = kprop
Init%Props(1:kprop, 1:Init%PropSetsCol) = TempProps

! deallocate temp matrices
IF (ALLOCATED(TempProps)) DEALLOCATE(TempProps)
IF (ALLOCATED(TempMembers)) DEALLOCATE(TempMembers)



!--------------------------------------
! write discretized data to a txt file
CALL GetNewUnit( UnOut ) 

OutFile = (trim(Init%RootName)//'_Descritize.txt' )
CALL OpenFOutFile ( UnOut, OutFile , ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnOut )
   RETURN
END IF

WRITE(UnOut, '(4(1x, I5))' ) Init%NNode, Init%NElem, Init%NProp, Init%NCMass
WRITE(UnOut, '(F6.0, E15.6, E15.6, E15.6)') ((Init%Nodes(i, j), j = 1, Init%JointsCol), i = 1, Init%NNode)
WRITE(UnOut, '(5(1x, I6))') ((Init%Elems(i, j), j = 1, Init%MembersCol), i = 1, Init%NElem)
WRITE(UnOut, '(F6.0, E15.6, E15.6, E15.6, E15.6, E15.6 ) ') ((Init%Props(i, j), j = 1, 6), i = 1, Init%NProp)
WRITE(UnOut, '(I6)')Init%NReact*6
WRITE(UnOut, '(I6, I6)') ((Init%BCs(i, j), j = 1, 2), i = 1, Init%Nreact*6)
WRITE(UnOut, '(I6)')Init%NInterf*6
WRITE(UnOut, '(I6, I6)') ((Init%IntFc(i, j), j = 1, 2), i = 1, Init%NInterf*6)
WRITE(UnOut, '(I6)')Init%NCMass
WRITE(UnOut, '(F6.0, E15.6, E15.6, E15.6, E15.6)') ((Init%Cmass(i, j), j = 1, 5), i = 1, Init%NCMass)
WRITE(UnOut, '(I6)')Init%NMembers

WRITE(tempStr, '(I10)') Init%NDiv + 1 
OutFmt = ('('//trim(tempStr)//'(I6))')
WRITE(UnOut, trim(OutFmt)) ((Init%MemberNodes(i, j), j = 1, Init%NDiv+1), i = 1, Init%NMembers)

CLOSE(UnOut)

END SUBROUTINE SubDyn_Discrt
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewNode(k, x, y, z, Init)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init
   
   ! local variables
   INTEGER          :: k
   REAL(ReKi)       :: x, y, z
   
   Init%Nodes(k, 1) = k
   Init%Nodes(k, 2) = x
   Init%Nodes(k, 3) = y
   Init%Nodes(k, 4) = z


END SUBROUTINE GetNewNode
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewElem(k, n1, n2, p1, p2, Init)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init
   
   ! local variables
   INTEGER          :: k, n1, n2, p1, p2
   
   Init%Elems(k, 1) = k
   Init%Elems(k, 2) = n1
   Init%Elems(k, 3) = n2
   Init%Elems(k, 4) = p1
   Init%Elems(k, 5) = p2

END SUBROUTINE GetNewElem
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewProp(k, E, G, rho, d, t, TempProps, NTempProps, PropCol)
   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init
   
   ! local variables
   INTEGER                  :: k, NTempProps, PropCol
   REAL(ReKi)               :: E, G, rho, d, t
   REAL(ReKi)               :: TempProps(NTempProps, PropCol)
   
   TempProps(k, 1) = k
   TempProps(k, 2) = E
   TempProps(k, 3) = G
   TempProps(k, 4) = rho
   TempProps(k, 5) = d
   TempProps(k, 6) = t

END SUBROUTINE GetNewProp
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE InitIAJA(Init)
   ! for 2-node element only
   USE NWTC_Library
   USE SubDyn_Types
   USE qsort_c_module 
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init

   ! local variables
   INTEGER                      :: I, J, k, UnDbg, Sttus, ERRSTAT, r, s
   INTEGER                      :: NNZ, TDOF
   INTEGER, ALLOCATABLE         :: IA(:), JA(:)
   CHARACTER(1024)              :: OutFile, tempStr, OutFmt
   INTEGER, ALLOCATABLE         :: Col_Arr(:)
   INTEGER                      :: ND
   INTEGER                      :: SortA(Init%MaxMemjnt,2)

   ! allocate for NodesConnE and NodesConnN
   ALLOCATE(Init%NodesConnE(Init%NNode, Init%MaxMemJnt+2), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating NodesConnE matrix', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%NodesConnE = 0

   ALLOCATE(Init%NodesConnN(Init%NNode, Init%MaxMemJnt+2), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating NodesConnE matrix', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%NodesConnN = 0

   ! find the node connectivity, nodes/elements that connect to a common node
      ! initialize the temp array for sorting
   SortA(Init%MaxMemjnt,2) = 0

   DO I = 1, Init%NNode
      Init%NodesConnE(I, 1) = Init%Nodes(I, 1)
      Init%NodesConnN(I, 1) = Init%Nodes(I, 1)
      
      k = 0
      DO J = 1, Init%NElem
         IF ( Init%Nodes(I, 1)==Init%Elems(J, 2) ) THEN
            k = k + 1
            Init%NodesConnE(I, k + 2) = Init%Elems(J, 1)
            Init%NodesConnN(I, k + 2) = Init%Elems(J, 3)
         ENDIF
         IF ( Init%Nodes(I, 1)==Init%Elems(J, 3) ) THEN
            k = k + 1
            Init%NodesConnE(I, k + 2) = Init%Elems(J, 1)
            Init%NodesConnN(I, k + 2) = Init%Elems(J, 2)
         ENDIF
      ENDDO
      
      IF( k>1 )THEN ! sort the nodes ascendingly
         SortA(1:k, 1) = Init%NodesConnN(I, 3:(k+2))
         CALL QsortC( SortA(1:k, 1:2) )
         Init%NodesConnN(I, 3:(k+2)) = SortA(1:k, 1)
      ENDIF
      
      Init%NodesConnE(I, 2) = k
      Init%NodesConnN(I, 2) = k
   ENDDO
   


! allocate the column array - column numbers that have nonzero component
   ALLOCATE(Col_Arr( 6*(Init%MaxMemJnt+1)), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating Col_Arr', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Col_Arr = 0

! allocate the row array IA
   TDOF = 6*Init%NNode
   Init%TDOF = TDOF

   ALLOCATE(IA( TDOF + 1 ), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating IA', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   IA = 0
   
! allocate the row array JA
   ALLOCATE(JA( (TDOF*TDOF+TDOF)/2 ), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating IA', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   JA = 0
   

   ! for each node, find the nodes that connect to this node
   ! find the target columns in the global K and M for these nodes
   NNZ = 0
   DO I = 1, Init%NNode
     
      DO J = 1, 6 ! the common node: first 6 columns of row I
         Col_Arr(J) = (I-1)*6 + J
      ENDDO
      
      k = 1
      DO J = 1, Init%NodesConnN(I, 2) 
         nd = Init%NodesConnN(I, J+2) ! nodes that connect to the common node
         IF ( nd > I ) THEN ! only count the node number that is greater than the common node
            k = k + 1
            DO r = 1, 6
               Col_Arr( (k-1)*6 + r ) = (nd-1)*6 + r
            ENDDO ! r
         ENDIF
         
      ENDDO ! J
      
!     write(*, *) ' Col_Arr '
!     write(*, *) 'I = ', I
!     WRITE(tempStr, '(I10)') k*6 
!     OutFmt = ('('//trim(tempStr)//'(I5))')
!     WRITE(*, trim(OutFmt)) (Col_Arr(1:k*6))
      
      
     ! total number of columns is k*6
     ! only count the diagonal and the upper triangle 
     DO r = 1, 6
        IA( (I-1)*6 + r ) = NNZ + 1
        DO s = r, k*6
            NNZ = NNZ + 1
            JA(NNZ) = Col_Arr(s)
        ENDDO ! s
                 
     ENDDO ! r
   
   ENDDO ! I
   IA( TDOF + 1 ) = NNZ+1
   
   
! allocate the row array IA
   ALLOCATE(Init%IA( TDOF +1 ), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating Init%IA', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%IA = IA
   
! allocate the row array JA
   ALLOCATE(Init%JA(NNZ), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating Init%JA', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%JA(1:NNZ) = JA(1:NNZ)
   Init%NNZ = NNZ
   
! deallocate temp matrices
IF (ALLOCATED(IA)) DEALLOCATE(IA)
IF (ALLOCATED(JA)) DEALLOCATE(JA)
IF (ALLOCATED(Col_Arr)) DEALLOCATE(Col_Arr)

! test the qsortc subroutine
!TestA(:,1) = (/1,3,4,2,6/)
!TestA(:,2) = (/1,3,4,2,6/)

!CALL QsortC(TestA(1:5, 1:1))

!write(*, '(2(I5))') ((TestA(i, j), j = 1, 2), i = 1, 5)

!--------------------------------------
! write node connectivity data and IA, JA to a txt file
CALL GetNewUnit( UnDbg ) 

OutFile = (trim(Init%RootName)//'_Connectivity.txt' )
CALL OpenFOutFile ( UnDbg, OutFile , ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnDbg )
   RETURN
END IF


WRITE(tempStr, '(I10)') Init%MaxMemJnt + 2 
OutFmt = ('('//trim(tempStr)//'(I6))')
WRITE(UnDbg, '(A)') 'Elements connect to a common node'
WRITE(UnDbg, trim(OutFmt)) ((Init%NodesConnE(i, j), j = 1, Init%MaxMemJnt+2), i = 1, Init%NNode)
WRITE(UnDbg, '(A)') 'Nodes connect to a common node'
WRITE(UnDbg, trim(OutFmt)) ((Init%NodesConnN(i, j), j = 1, Init%MaxMemJnt+2), i = 1, Init%NNode)

WRITE(UnDbg, *) 'TDOF = ', TDOF
WRITE(UnDbg, '(I6, I6)') ( (i, Init%IA(i)), i = 1, TDOF+1)
Write(UnDbg, *) 'NNZ = ', NNZ
write(Undbg, '(I6, I6)') ( (i, Init%JA(i)), i = 1, NNZ)

CLOSE(UnDbg)

END SUBROUTINE InitIAJA
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE AssembleKM(Init)

   USE NWTC_Library
   USE SubDyn_Types
   IMPLICIT NONE

   TYPE(SD_InitInputType)   ::Init
   
   INTEGER                  :: I, J, K, Ja
   INTEGER                  :: NNE        ! number of nodes in one element
   INTEGER                  :: N1, N2     ! starting node and ending node in the element
   INTEGER                  :: P1, P2     ! property set numbers for starting and ending nodes

   REAL(ReKi)               :: D1, D2, t1, t2, E, G, rho ! properties of a section
   REAL(ReKi)               :: x1, y1, z1, x2, y2, z2    ! coordinates of the nodes
   REAL(ReKi)               :: DirCos(3, 3)              ! direction cosine matrices
   REAL(ReKi)               :: L                         ! length of the element
   REAL(ReKi)               :: r1, r2, t, Iyy, Jzz, Ixx, A, kappa
   LOGICAL                  :: shear
   REAL(ReKi), ALLOCATABLE  :: Ke(:,:), Me(:, :)         ! element stiffness and mass matrices
   INTEGER, ALLOCATABLE     :: nn(:)                     ! node number in element 
   INTEGER                  :: tgt_row_sys(6), row_in_elem(6), tgt_col_sys(6), col_in_elem(6)
   
   INTEGER                  :: ei, ej, ti, tj, r, s, beg_jA, end_jA, ii, jj
   CHARACTER(1024)          :: tempstr1, tempstr2, outfile
   INTEGER                  :: UnDbg, Sttus, ErrStat
   
   
   
      ! for current application
   IF ( (Init%FEMMod .LE. 3) .and. (Init%FEMMod .GE. 0)) THEN
      NNE = 2   
   ENDIF                              
   
   ! allocate element stiffness matrix
   ALLOCATE( Ke(NNE*6, NNE*6), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating element stiffness matrix Ke', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF

   ! allocate element mass matrix
   ALLOCATE( Me(NNE*6, NNE*6), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating element stiffness matrix Ke', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   
   
   ! allocate system stiffness matrix
   ALLOCATE( Init%K(Init%NNZ), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating system stiffness matrix K', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%K = 0

   ! allocate system mass matrix
   ALLOCATE( Init%M(Init%NNZ), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating element stiffness matrix Ke', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   Init%M = 0

   
   ! allocate node number in element array
   ALLOCATE( nn(NNE), STAT=Sttus)
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating nn', .TRUE. )
      ErrStat = 1
      RETURN
   ENDIF
   
      ! loop over all elements
   DO I = 1, Init%NElem
   
      DO J = 1, NNE
         NN(J) = Init%Elems(I, J + 1)
      ENDDO
   
      N1 = Init%Elems(I,       2)
      N2 = Init%Elems(I, NNE + 1)
      
      P1 = Init%Elems(I, NNE + 2)
      P2 = Init%Elems(I, NNE + 3)
      
      
      E   = Init%Props(P1, 2)
      G   = Init%Props(P1, 3)
      rho = Init%Props(P1, 4)
      D1  = Init%Props(P1, 5)
      t1  = Init%Props(P1, 6)
      D2  = Init%Props(P2, 5)
      t2  = Init%Props(P2, 6)
      
      x1  = Init%Nodes(N1, 2)
      y1  = Init%Nodes(N1, 3)
      z1  = Init%Nodes(N1, 4)
      
      x2  = Init%Nodes(N2, 2)
      y2  = Init%Nodes(N2, 3)
      z2  = Init%Nodes(N2, 4)

      CALL GetDirCos(X1, Y1, Z1, X2, Y2, Z2, DirCos, L)
      CALL SetConstants()
      
      IF ( (Init%FEMMod == 0).OR.(Init%FEMMod == 2)) THEN ! uniform element 
         r1 = 0.25*(D1 + D2)
         t  = 0.5*(t1+t2)
         
         IF ( EqualRealNos(t, 0.0) ) THEN
            r2 = 0
         ELSE
            r2 = r1 - t
         ENDIF
         
         A = Pi_D*(r1*r1-r2*r2)
         Ixx = 0.25*Pi_D*(r1**4-r2**4)
         Iyy = Ixx
         Jzz = 2.0*Ixx
         
         IF( Init%FEMMod == 0 ) THEN ! uniform Euler-Bernoulli
            Shear = .false.
            kappa = 0
         ELSEIF( Init%FEMMod == 2 ) THEN ! uniform Timoshenko
            Shear = .true.
            kappa = 0.5
         ENDIF
         
         CALL ElemK(A, L, Ixx, Iyy, Jzz, Shear, kappa, E, G, DirCos, Ke)
         CALL ElemM(A, L, Ixx, Iyy, Jzz, rho, DirCos, Me)
         
         
      ELSEIF  (Init%FEMMod == 1) THEN ! tapered Euler-Bernoulli
      
      
      ELSEIF  (Init%FEMMod == 3) THEN ! tapered Timoshenko
      ELSE
      ENDIF  
      
      !~~~~~~~~~~ assemble to system K and M in compressed row format 
      !           for two node element
      DO J = 1, NNE ! NNE = 2, 2 nodes in one element
      
         DO ja = 1, 6
            tgt_row_sys(ja) = ( nn(J) - 1 )*6 + ja
            row_in_elem(ja) = ( J - 1 )*6 + ja 
         ENDDO !ja 
         
         DO K = 1, NNE
         
            DO ja = 1, 6
               tgt_col_sys(ja) = ( nn(K) - 1 )*6 + ja
               col_in_elem(ja) = ( K - 1 )*6 + ja 
            ENDDO !ja 
         
            DO ii = 1, 6
               ei = row_in_elem(ii)
               ti = tgt_row_sys(ii)
            
               DO jj = 1, 6
                  ej = col_in_elem(jj)
                  tj = tgt_col_sys(jj)
                  
                  IF(ti .LE. tj) THEN ! store the upper triangle and the diagonal
                     beg_jA = Init%IA(ti)
                     end_jA = Init%IA(ti+1) - 1
                     
                     s = 0
                     DO r = beg_jA, end_jA
                        IF ( Init%JA(r) == tj ) THEN
                           s = r
                        ENDIF
                     ENDDO ! r
                     
                     IF ( s == 0) THEN
                        write(tempstr1, *)  ti
                        write(tempstr2, *)  tj 
                        CALL ProgAbort ( ('A( '//trim(tempstr1)//','//trim(tempstr2)//') not found in AssembleKM'), .TRUE. )
                        ErrStat = 1
                        RETURN
                     ENDIF
                     
                     Init%K(s) = Init%K(s) + Ke(ei, ej)
                     Init%M(s) = Init%M(s) + Me(ei, ej)
                     
                     
                  ENDIF
               
               
               ENDDO ! jj
            ENDDO !ii
         
         ENDDO ! K
         
      ENDDO ! J
      
      
                 
   
   ENDDO ! I end loop over elements
   
   
   ! add concentrated mass
   DO I = 1, Init%NCMass
      DO J = 1, 3
          r = ( Init%CMass(I, 1) - 1 )*6 + J
          Init%M(Init%IA(r)) = Init%M(Init%IA(r)) + Init%CMass(I, 2)
      ENDDO
      DO J = 4, 6
          r = ( Init%CMass(I, 1) - 1 )*6 + J
          Init%M(Init%IA(r)) = Init%M(Init%IA(r)) + Init%CMass(I, J-1)
      ENDDO

   ENDDO ! I concentrated mass
 

! deallocate temp matrices
IF (ALLOCATED(Ke)) DEALLOCATE(Ke)
IF (ALLOCATED(Me)) DEALLOCATE(Me)
IF (ALLOCATED(nn)) DEALLOCATE(nn)
   
   
!--------------------------------------
! write assembed K M to a txt file
CALL GetNewUnit( UnDbg ) 

OutFile = (trim(Init%RootName)//'_Compressed_K_M.txt' )
CALL OpenFOutFile ( UnDbg, OutFile , ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnDbg )
   RETURN
END IF

WRITE(UnDbg, '(A, I6)') ('Number of Elements in K', Init%IA(Init%TDOF+1)-1 )
WRITE(UnDbg, '(I6, e15.6)') ( (i, Init%K(i)), i = 1, Init%IA(Init%TDOF+1)-1)

WRITE(UnDbg, '(A, I6)') ('Number of Elements in M', Init%IA(Init%TDOF+1)-1 )
WRITE(UnDbg, '(I6, e15.6)') ( (i, Init%M(i)), i = 1, Init%IA(Init%TDOF+1)-1)


CLOSE(UnDbg)


END SUBROUTINE AssembleKM
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE GetDirCos(X1, Y1, Z1, X2, Y2, Z2, DirCos, xyz)
   
   USE NWTC_Library
   IMPLICIT NONE

   REAL(ReKi)         :: x1, y1, z1, x2, y2, z2
   REAL(ReKi)         :: DirCos(3, 3)
   
   REAL(ReKi)         :: xz, xyz
   
   integer            :: ErrStat
   
   xz = sqrt( (x1-x2)**2 + (z1-z2)**2 )
   xyz = sqrt( (x1-x2)**2 + (z1-z2)**2 + (y1-y2)**2 )
   
   IF ( EqualRealNos(xyz, 0.0) ) THEN
      CALL ProgAbort ( ' Same starting and ending location in the element.', .TRUE. )
      ErrStat = 4
      RETURN
   ENDIF
   
   DirCos = 0
   
   IF ( EqualRealNos(xz, 0.0) ) THEN
      IF( y2 < y1) THEN
         DirCos(1, 1) = 1.0
         DirCos(2, 3) = -1.0
         DirCos(3, 2) = 1.0
      ELSE
         DirCos(1, 1) = 1.0
         DirCos(2, 3) = 1.0
         DirCos(3, 2) = -1.0

      ENDIF
   ELSE
      DirCos(1, 1) = -(z1-z2)/xz
      DirCos(1, 2) = -(x1-x2)*(y1-y2)/(xz*xyz)
      DirCos(1, 3) = (x2-x1)/xyz

      DirCos(2, 2) = xz/xyz
      DirCos(2, 3) = (y2-y1)/xyz

      DirCos(3, 1) = -(x2-x1)/xz
      DirCos(3, 2) = -(y1-y2)*(z1-z2)/(xz*xyz)
      DirCos(3, 3) = (z2-z1)/xyz
   ENDIF

END SUBROUTINE GetDirCos
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE ElemK(A, L, Ixx, Iyy, Jzz, Shear, kappa, E, G, DirCos, K)
   ! element stiffness matrix for classical beam elements
   ! shear is true  -- non-tapered Timoshenko beam 
   ! shear is false -- non-tapered Euler-Bernoulli beam 
   USE NWTC_Library
   IMPLICIT NONE

   REAL(ReKi), INTENT( IN)               :: A, L, Ixx, Iyy, Jzz, E, G, kappa
   REAL(ReKi), INTENT( IN)               :: DirCos(3,3)
   LOGICAL, INTENT( IN)                  :: Shear
   
   REAL(ReKi)             :: K(12, 12)
         
   REAL(ReKi)                            :: Ax, Ay, Kx, Ky
   REAL(ReKi)                            :: DC(12, 12)
   
   Ax = kappa*A
   Ay = kappa*A
   
   K = 0
   
   IF (Shear) THEN
      Kx = 12.0*E*Iyy / (G*Ax*L*L)
      Ky = 12.0*E*Ixx / (G*Ay*L*L)
   ELSE
      Kx = 0.0
      Ky = 0.0
   ENDIF
      
   K( 9,  9) = E*A/L
   K( 7,  7) = 12.0*E*Iyy/( L*L*L*(1.0 + Kx) )
   K( 8,  8) = 12.0*E*Ixx/( L*L*L*(1.0 + Ky) )
   K(12, 12) = G*Jzz/L
   K(10, 10) = (4.0 + Ky)*E*Ixx / ( L*(1.0+Ky) )  
   K(11, 11) = (4.0 + Kx)*E*Iyy / ( L*(1.0+Kx) )
   K( 2,  4) = -6.*E*Ixx / ( L*L*(1.0+Ky) )
   K( 1,  5) =  6.*E*Iyy / ( L*L*(1.0+Kx) )
   K( 4, 10) = (2.0-Ky)*E*Ixx / ( L*(1.0+Ky) )
   K( 5, 11) = (2.0-Kx)*E*Iyy / ( L*(1.0+Kx) )
   
   K( 3,  3)  = K(9,9)
   K( 1,  1)  = K(7,7)
   K( 2,  2)  = K(8,8)
   K( 6,  6)  = K(12,12)
   K( 4,  4)  = K(10,10)
   K(5,5)  = K(11,11)
   K(4,2)  = K(2,4)
   K(5,1)  = K(1,5)
   K(10,4) = K(4,10)
   K(11,5) = K(5,11)
   K(12,6)= -K(6,6)
   K(10,2)=  K(4,2)
   K(11,1)=  K(5,1)
   K(9,3) = -K(3,3)
   K(7,1) = -K(1,1)
   K(8,2) = -K(2,2)
   K(6, 12) = -K(6,6)
   K(2, 10) =  K(4,2)
   K(1, 11) =  K(5,1)
   K(3, 9)  = -K(3,3)
   K(1, 7)  = -K(1,1)
   K(2, 8)  = -K(2,2)
   K(11,7) = -K(5,1)
   K(10,8) = -K(4,2)
   K(7,11) = -K(5,1)
   K(8,10) = -K(4,2)
   K(7,5) = -K(5,1)
   K(5,7) = -K(5,1)
   K(8,4) = -K(4,2)
   K(4,8) = -K(4,2)
   
   DC = 0
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   K = MATMUL( MATMUL(DC, K), TRANSPOSE(DC) )
   
   !write(*, *) K - TRANSPOSE(K)

END SUBROUTINE ElemK
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE ElemM(A, L, Ixx, Iyy, Jzz, rho, DirCos, M)
   ! element mass matrix for classical beam elements

   USE NWTC_Library
   IMPLICIT NONE

   REAL(ReKi), INTENT( IN)               :: A, L, Ixx, Iyy, Jzz, rho
   REAL(ReKi), INTENT( IN)               :: DirCos(3,3)
   
   REAL(ReKi)             :: M(12, 12)
         
   REAL(ReKi)                            :: t, rx, ry, po
   REAL(ReKi)                            :: DC(12, 12)
   
   t = rho*A*L;
   rx = rho*Ixx;
   ry = rho*Iyy;
   po = rho*Jzz*L;   

   M = 0
   
      
   M( 9,  9) = t/3.0
   M( 7,  7) = 13.0*t/35.0 + 6.0*ry/(5.0*L)
   M( 8,  8) = 13.0*t/35.0 + 6.0*rx/(5.0*L)
   M(12, 12) = po/3.0
   M(10, 10) = t*L*L/105.0 + 2.0*L*rx/15.0
   M(11, 11) = t*L*L/105.0 + 2.0*L*ry/15.0
   M( 2,  4) = -11.0*t*L/210.0 - rx/10.0
   M( 1,  5) =  11.0*t*L/210.0 + ry/10.0
   M( 3,  9) = t/6.0
   M( 5,  7) =  13.*t*L/420. - ry/10.
   M( 4,  8) = -13.*t*L/420. + rx/10. 
   M( 6, 12) = po/6.
   M( 2, 10) =  13.*t*L/420. - rx/10. 
   M( 1, 11) = -13.*t*L/420. + ry/10.
   M( 8, 10) =  11.*t*L/210. + rx/10.
   M( 7, 11) = -11.*t*L/210. - ry/10. 
   M( 1,  7) =  9.*t/70. - 6.*ry/(5.*L)
   M( 2,  8) =  9.*t/70. - 6.*rx/(5.*L)
   M( 4, 10) = -L*L*t/140. - rx*L/30. 
   M( 5, 11) = -L*L*t/140. - ry*L/30.
   
   M( 3,  3) = M( 9,  9)
   M( 1,  1) = M( 7,  7)
   M( 2,  2) = M( 8,  8)
   M( 6,  6) = M(12, 12)
   M( 4,  4) = M(10, 10)
   M( 5,  5) = M(11, 11)
   M( 4,  2) = M( 2,  4)
   M( 5,  1) = M( 1,  5)
   M( 9,  3) = M( 3,  9)
   M( 7,  5) = M( 5,  7)
   M( 8,  4) = M( 4,  8)
   M(12,  6) = M( 6, 12)
   M(10,  2) = M( 2, 10)
   M(11,  1) = M( 1, 11)
   M(10,  8) = M( 8, 10)
   M(11,  7) = M( 7, 11)
   M( 7,  1) = M( 1,  7)
   M( 8,  2) = M( 2,  8)
   M(10,  4) = M( 4, 10)
   M(11,  5) = M( 5, 11)
   
   DC = 0
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   M = MATMUL( MATMUL(DC, M), TRANSPOSE(DC) )

END SUBROUTINE ElemM


!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE EigenSolve(Init)


   USE NWTC_Library
   USE SubDyn_Types
   USE HSL_ZD11_double
   USE EA16_INTERFACE

   IMPLICIT NONE

   TYPE(ZD11_TYPE) MATRIXK
   TYPE(ZD11_TYPE) MATRIXM

   TYPE(SD_InitInputType)   :: Init
   

   INTEGER                  :: NOmega
   REAL(8),ALLOCATABLE               :: Omega(:), Phi(:, :)
   
   INTEGER                  :: UnDbg, ERRSTAT, I, J, ii, jj
   CHARACTER(1024)          :: TempStr, outfile
   
   
   REAL(8)               :: KT(24, 24), MT(24, 24)
   INTEGER                  :: IIA(Init%NNZ), i1, i2
   
   NOmega = 10
   ALLOCATE(Omega(NOmega))
   ALLOCATE(Phi(Init%TDOF, NOmega))
   
   
   !===============================================================================
	!=====             Construct Matrix Structure
	!===============================================================================


    DO i = 1, Init%TDOF
        i1 = Init%IA(i)
        i2 = Init%IA(i+1) - 1
        DO j = i1, i2
            IIA(j) = i
        END DO

    END DO


    ! ------MATRIXK-------
    MATRIXK%N = Init%TDOF
    MATRIXK%NE = Init%NNZ

    ALLOCATE(MATRIXK%COL(Init%NNZ),MATRIXK%ROW(Init%NNZ),MATRIXK%val(Init%NNZ))
   
    MATRIXK%ROW = IIA
    MATRIXK%COL = Init%JA
    MATRIXK%VAL = Init%K

    ! ------MATRIXM-------
    MATRIXM%N = Init%TDOF
    MATRIXM%NE = Init%NNZ

    ALLOCATE(MATRIXM%COL(Init%NNZ),MATRIXM%ROW(Init%NNZ),MATRIXM%val(Init%NNZ))

    MATRIXM%ROW = IIA
    MATRIXM%COL = Init%JA
    MATRIXM%VAL = Init%M

! deallocate temp matrices
IF (ALLOCATED(Init%K)) DEALLOCATE(Init%K)
IF (ALLOCATED(Init%M)) DEALLOCATE(Init%M)
IF (ALLOCATED(Init%JA)) DEALLOCATE(Init%JA)


    Call EA16_double(Init%TDOF, Init%NNZ, MatrixK, Init%IA,         &
                               Init%NNZ, MatrixM, Init%IA,          &
                               2, 1, 24,                            &
                               NOmega, Omega, phi)


IF (ALLOCATED(Init%IA)) DEALLOCATE(Init%IA)

   !===============================================================================
	!=====             Finish EigenSolve
	!===============================================================================



!--------------------------------------
! write assembed K M to a txt file
CALL GetNewUnit( UnDbg ) 

OutFile = (trim(Init%RootName)//'_eigen_results.txt' )
CALL OpenFOutFile ( UnDbg, OutFile , ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnDbg )
   RETURN
END IF

!write(UnDbg, '(24(1x, e15.6))') ((KT(i, j), j= 1, 24), i = 1, 24)
!write(UnDbg, '(24(1x, e15.6))') ((MT(i, j), j= 1, 24), i = 1, 24)
DO i = 1, NOmega
   Omega(i) = sqrt(Omega(i))/2.0/Pi
ENDDO 

WRITE(UnDbg, '(A, I6)') ('Number of eigen values ', NOmega )
WRITE(UnDbg, '(I6, e15.6)') ( (i, Omega(i) ), i = 1, NOmega )


CLOSE(UnDbg)


END SUBROUTINE EigenSolve

!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE ApplyConstr(Init)

   USE NWTC_Library
  
   USE SubDyn_Types
   
   IMPLICIT NONE

   TYPE(SD_InitInputType)   :: Init
   
   INTEGER                  :: I, J, k
   INTEGER                  :: bgn_j, end_j, row_n
   
   
   DO I = 1, Init%NReact*6
      row_n = Init%BCs(I, 1)
      
      bgn_j =  Init%IA(row_n)
      end_j =  Init%IA(row_n+1) - 1
      
      Init%K(bgn_j) = 1
      Init%M(bgn_j) = 10**(-5)
      
      DO J = bgn_j + 1, end_j
      
         Init%K(J) = 0
         Init%M(J) = 0
         
      ENDDO ! J
      
      
   ENDDO ! I

   

END SUBROUTINE ApplyConstr
