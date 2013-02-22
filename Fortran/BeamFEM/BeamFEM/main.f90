PROGRAM main
   USE NWTC_Library
   IMPLICIT NONE

   CHARACTER(1024) InputFile
   
   InputFile = 'testframe.txt'
   
   CALL BEAMFEM(InputFile)

END PROGRAM main
!----------------------------------------------------------------------------
SUBROUTINE BEAMFEM(InputFile)

USE NWTC_Library
IMPLICIT NONE


CHARACTER(1024)              ::   InputFile

INTEGER(IntKi)               ::    NJoints
INTEGER(IntKi), PARAMETER    ::    JointsDm = 4
INTEGER(IntKi)               ::    NMembers
INTEGER(IntKi), PARAMETER    ::    MemberDm = 5
INTEGER(IntKi)               ::    NPropSets
INTEGER(IntKi), PARAMETER    ::    PropSetsDm = 6
INTEGER(IntKi)               ::    NXPropSets
INTEGER(IntKi), PARAMETER    ::    XPropSetDm = 10
INTEGER(IntKi)               ::    NReact
INTEGER(IntKi), PARAMETER    ::    ReactDm = 7
INTEGER(IntKi)               ::    NInterf
INTEGER(IntKi), PARAMETER    ::    InterfDm = 7
INTEGER(IntKi)               ::    NForce  
INTEGER(IntKi)               ::    NCMass ! number of concentrated mass  
INTEGER(IntKi), PARAMETER    ::    CMassDm = 5
INTEGER(IntKi)               ::    NCOSMs
INTEGER(IntKi), PARAMETER    ::    COSMsDm = 10

INTEGER(IntKi)               ::    TDOF

REAL(ReKi), ALLOCATABLE      ::     Joints(:, :)
REAL(ReKi), ALLOCATABLE      ::     PropSets(:, :)
REAL(ReKi), ALLOCATABLE      ::     XPropSets(:, :)
REAL(ReKi), ALLOCATABLE      ::     Forces(:, :)
REAL(ReKi), ALLOCATABLE      :: COSMs(:, :)
REAL(ReKi), ALLOCATABLE      :: CMass(:, :)
INTEGER(ReKi), ALLOCATABLE      ::     Members(:, :)
INTEGER(IntKi), ALLOCATABLE      ::     Reacts(:, :)
INTEGER(IntKi), ALLOCATABLE      ::     Interf(:, :)

INTEGER(IntKi)               :: FEMMod                                          ! FEM switch: element model in the FEM: 0= Euler-Bernoulli(E-B) ; 1=Tapered E-B; 2= 2-node Timoshenko;  3= 2-node tapered Timoshenko 
INTEGER(IntKi)               :: NDiv                                            ! Number of divisions per member
LOGICAL                      :: CBMod                                           ! Perform C-B flag
INTEGER(IntKi)               :: Nmodes                                          ! Number of modes to retain
REAL(ReKi), ALLOCATABLE      :: JDampings(:)                                       ! Damping coefficient associate with internal modes 

! ---------   ---------------------------------------------------------------------------------

INTEGER(IntKi)               :: NNode
INTEGER(IntKi)               :: NElem
INTEGER(IntKi)               :: NProp

REAL(ReKi), ALLOCATABLE      :: Nodes(:, :)
REAL(ReKi), ALLOCATABLE      :: Props(:, :)

INTEGER, ALLOCATABLE         :: Elems(:, :)


REAL(ReKi), ALLOCATABLE      ::     K(:), M(:), F(:)
REAL(ReKi), ALLOCATABLE      ::     IA(:), JA(:)

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
CHARACTER(1024)              :: RootName                                        ! The root name of the input and output files.
INTEGER(IntKi)               :: UnIn
INTEGER(IntKi)               :: UnOut 
INTEGER(IntKi)               :: ErrStat



INTEGER(IntKi)               :: I



CALL GetNewUnit( UnIn )   
CALL OpenFInpfile(UnIn, TRIM(InputFile), ErrStat)

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

CALL GetPath( InputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.
CALL GetRoot( InputFile, RootName )

write(*,*) ('root:', trim(rootname))

!-------------------------- HEADER ---------------------------------------------

   CALL ReadCom( UnIn, InputFile, 'SubDyn input file header line 1', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


   CALL ReadCom( UnIn, InputFile, 'SubDyn input file header line 2', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   CALL ReadCom( UnIn, InputFile, 'SubDyn input file header line 3', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

!-------------------------- SIMULATION CONTROL PARAMETERS ----------------------

      ! Skip the comment line.

   CALL ReadCom( UnIn, InputFile, ' SIMULATION CONTROL PARAMETERS ', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   ! Echo - Echo input to "echo.out".

READ (UnIn,*,IOSTAT=IOS)  WrEcho
CALL CheckIOS ( IOS, InputFile, 'Echo', FlgType )
Echo = WrEcho
write(*,*) Echo

IF ( Echo )  THEN
   CALL OpenEcho ( UnEc, TRIM(RootName)//'.ech' )
   WRITE (UnEc,'(/,A)'   )  'This file of echoed input was generated by '//TRIM(ProgName)//' '//TRIM(ProgVer)// &
                            ' on '//CurDate()//' at '//CurTime()//'.'
   WRITE (UnEc,'(/,A,/)' )  'Substructure data from file "'//TRIM( InputFile )//'":'
   WRITE (UnEc,Frmt      )  Echo, 'Echo', 'Echo input to "echo.out"'
ENDIF

!-------------------- FEA and CRAIG-BAMPTON PARAMETERS---------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' FEA and CRAIG-BAMPTON PARAMETERS ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! FEMMod - FEM switch: element model in the FEM: 0= Euler-Bernoulli(E-B) ; 1=Tapered E-B; 2= 2-node Timoshenko;  3= 2-node tapered Timoshenko

CALL ReadIVar ( UnIn, InputFile, FEMMod, 'FEMMod', 'FEM analysis mode' )

IF ( ( FEMMod < 0 ) .OR. ( FEMMod > 4 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': FEMMod must be 0, 1, 2, or 3.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF


   ! NDiv - Number sub-elements per member

CALL ReadIVar ( UnIn, InputFile, NDiv, 'NDiv', 'Number of divisions per member' )

IF ( ( NDiv < 0 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': NDiv must be a non-negative integer', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

write(*,*) 'NDiv = ', NDiv


   ! CBMod - Perform C-B flag.
READ( UnIn, *, IOSTAT=ErrStat ) CBMod

IF ( ErrStat /= 0 ) THEN
   CALL CheckIOS ( ErrStat, InputFile, 'CBMod', NumType, .TRUE. )
   CLOSE( UnIn )
   RETURN
END IF

write(*,*) 'CBMod = ', CBMod


IF (CBMod) THEN

      ! Nmodes - Number of interal modes to retain.
   CALL ReadIVar ( UnIn, InputFile, Nmodes, 'Nmodes', 'Number of internal modes' )

   IF ( ( Nmodes < 1 ) )  THEN
      CALL ProgAbort ( ' FEMMod must be a positive integer.' )
   ENDIF
   
      ! Damping ratios for retained modes
   ALLOCATE(JDampings(Nmodes), STAT=Sttus)
   
   IF ( Sttus /= 0 )  THEN
      CALL ProgAbort ( ' Error allocating memory for the damping ratio array.' )
   ENDIF

   CALL ReadRAry( UnIn, InputFile, JDampings, Nmodes, 'JDamping', 'Damping ratio of the internal modes', ErrStat )
      
   DO I = 1, Nmodes
      IF ( ( JDampings(I) .LE. 0 ) .OR.( JDampings(I) .GE. 100 ) ) THEN
         CALL ProgAbort ( ' Damping ratio should be larger than 0 and less than 100.' )
      ENDIF
      
      JDampings(I) = JDampings(I)/100.0
      
   ENDDO
   
ELSE
      ! skip 2 lines
   READ (UnIn,'(A)',IOSTAT=IOS)  Comment
   CALL CheckIOS( IOS, InputFile, 'Nmodes ', StrType )
   READ (UnIn,'(A)',IOSTAT=IOS)  Comment
   CALL CheckIOS( IOS, InputFile, 'Damping ratio', StrType )
      
   Nmodes = 0
ENDIF

write(*,*) ('Nmodes = ', Nmodes)
write(*,*) ('Damping = ', JDampings)

!---- STRUCTURE JOINTS: joints connect structure members -----------------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' Joints connect structure members ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of joints
CALL ReadIVar ( UnIn, InputFile, NJoints, 'NJoints', 'Number of joints' )
IF ( ( NJoints < 2 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': NJoints must be greater than 1', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip two lines
CALL ReadCom( UnIn, InputFile, ' Joints description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

CALL ReadCom( UnIn, InputFile, ' units ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
   
   ! Joints coordinates
ALLOCATE(Joints(NJoints, JointsDm), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': Error allocating Joints arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

DO I = 1, NJoints

   CALL ReadRAry( UnIn, InputFile, Joints(I,1:JointsDm), JointsDm, 'Joints', 'Joint number and coordinates', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   write(*, *) Joints(I, 1:JointsDm)
   
ENDDO

!------------------- BASE REACTION JOINTS: T/F for Locked/Free DOF @ each Reaction Node ---------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' BASE REACTION JOINTS ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! Number of reaction joints (The joints should be all fixes for now)
CALL ReadIVar ( UnIn, InputFile, NReact, 'NReact', 'Number of joints with reaction forces' )
IF ( ( NReact < 1 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': NReact must be greater than 0', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF
   
   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' BASE REACTION JOINTS, description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

write(*, *) ('NReact = ', NReact)

   ! Joints with reaction forces, joint number and locked/free dof
ALLOCATE(Reacts(NReact, ReactDm), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': Error allocating Reacts arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF
   
   
DO I = 1, NReact

   CALL ReadAry( UnIn, InputFile, Reacts(I,1:ReactDm), ReactDm, 'Reacts', 'Joint number and dof', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   write(*, *) Reacts(I, 1:ReactDm)
   
ENDDO


!------- INTERFACE JOINTS: T/F for Locked (to the TP)/Free DOF @each Interface Joint (only Locked-to-TP implemented thus far (=rigid TP)) ---------

   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' INTERFACE JOINTS ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! Number of interface joints (The joints should be all fixes for now)
CALL ReadIVar ( UnIn, InputFile, NInterf, 'NInterf', 'Number of joints fixed to TP' )
IF ( ( NInterf < 0 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': NInterf must be non-negative.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' INTERFACE JOINTS, description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

write(*, *) ('NInterf = ', NInterf)

   ! Joints with reaction forces, joint number and locked/free dof
ALLOCATE(Interf(NInterf, InterfDm), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': Error allocating Interf arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF
   
   
DO I = 1, NInterf

   CALL ReadAry( UnIn, InputFile, Interf(I,1:InterfDm), InterfDm, 'Interf', 'Interface joint number and dof', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   write(*, *) Reacts(I, 1:InterfDm)
   
ENDDO
   
!----------------------------------- MEMBERS --------------------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' Members ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of members
CALL ReadIVar ( UnIn, InputFile, NMembers, 'NMembers', 'Number of members' )
IF ( ( NMembers < 1 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': NMembers must be greater than 0', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip one line
CALL ReadCom( UnIn, InputFile, ' Members description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
   
   ! Member connection
ALLOCATE(Members(Nmembers, MemberDm), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': Error allocating Members arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

write(*, *) ('NMembers = ', NMembers)

DO I = 1, NMembers

   CALL ReadAry( UnIn, InputFile, Members(I,1:MemberDm), MemberDm, 'Members', 'Member number and connectivity ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   write(*, *) Members(I, 1:MemberDm)
   
ENDDO   

!------------------ MEMBER X-SECTION PROPERTY data 1/2 [isotropic material for now: use this table if circular-tubular elements ------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' Member X-Section property 1/2 ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of property sets
CALL ReadIVar ( UnIn, InputFile, NPropSets, 'NPropSets', 'Number of property sets' )
IF ( ( NPropSets < 1 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': NPropSets must be greater than 0', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip two lines
CALL ReadCom( UnIn, InputFile, ' Property sets 1/2 description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
CALL ReadCom( UnIn, InputFile, ' Property sets 1/2 units ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
   
   ! Property sets value
ALLOCATE(PropSets(NPropSets, PropSetsDm), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': Error allocating PropSets arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

write(*, *) ('NPropSets = ', NPropSets)

DO I = 1, NPropSets

   CALL ReadRAry( UnIn, InputFile, PropSets(I,1:PropSetsDm), PropSetsDm, 'PropSets', 'PropSets number and values ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   write(*, *) PropSets(I, 1:PropSetsDm)
   
ENDDO   

!------------------ MEMBER X-SECTION PROPERTY data 2/2 [isotropic material for now: use this table if any section other than circular, however provide COSM(i,j) below) ------------------------


   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' Member X-Section property 2/2 ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of property sets
CALL ReadIVar ( UnIn, InputFile, NXPropSets, 'NPropSets', 'Number of property sets' )
IF ( ( NPropSets < 0 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': NPropSets must be non-negative.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip two lines
CALL ReadCom( UnIn, InputFile, ' Property sets 2/2 description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
CALL ReadCom( UnIn, InputFile, ' Property sets 2/2 units ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF
   
   ! Property sets value
ALLOCATE(XPropSets(NXPropSets, XPropSetDm), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': Error allocating PropSets arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

write(*, *) ('NXPropSets = ', NXPropSets)

DO I = 1, NXPropSets

   CALL ReadRAry( UnIn, InputFile, XPropSets(I,1:XPropSetDm), XPropSetDm, 'PropSets', 'PropSets number and values ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   write(*, *) XPropSets(I, 1:XPropSetDm)
   
ENDDO   

!---------------------- MEMBER COSINE MATRICES COSM(i,j) ------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' Member direction cosine matrices ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of direction cosine matrices
CALL ReadIVar ( UnIn, InputFile, NCOSMs, 'NCOSMs', 'Number of unique direction cosine matrices' )
IF ( ( NPropSets < 0 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': NCOSMs must be non-negative.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip one line
CALL ReadCom( UnIn, InputFile, ' Property sets 2/2 description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! Direction cosine matrices value
ALLOCATE(COSMs(NCOSMs, COSMsDm), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': Error allocating COSMs arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

write(*, *) ('NCOSMs = ', NCOSMs)

DO I = 1, NCOSMs

   CALL ReadRAry( UnIn, InputFile, COSMs(I,1:COSMsDm), COSMsDm, 'PropSets', 'PropSets number and values ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   write(*, *) COSMs(I, 1:COSMsDm)
   
ENDDO   

!------------------------ JOINT ADDITIONAL CONCENTRATED MASSES--------------------------

   ! Skip the comment line.

CALL ReadCom( UnIn, InputFile, ' Additional concentrated masses at joints ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! number of joints that have concentrated masses
CALL ReadIVar ( UnIn, InputFile, NCMass, 'NCMass', 'Number of joints that have concentrated masses' )
IF ( ( NPropSets < 0 ) )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': NCMass must be non-negative.', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

   ! Skip two lines
CALL ReadCom( UnIn, InputFile, ' Concentrated mass description ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

CALL ReadCom( UnIn, InputFile, ' Concentrated mass units ', ErrStat )

IF ( ErrStat /= 0 ) THEN
   CLOSE( UnIn )
   RETURN
END IF

   ! Concentrated mass value
ALLOCATE(CMass(NCMass, CMassDm), STAT=Sttus)
   
IF ( Sttus /= 0 )  THEN
   CALL ProgAbort ( ' Error in file "'//TRIM(InputFile)//': Error allocating CMass arrays', .TRUE. )
   ErrStat = 1
   CLOSE( UnIn )
   RETURN
ENDIF

write(*, *) ('NCMass = ', NCMass)

DO I = 1, NCMass

   CALL ReadRAry( UnIn, InputFile, CMass(I,1:CMassDm), CMassDm, 'CMass', 'Joint number and mass values ', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   write(*, *) CMass(I, 1:CMassDm)
   
ENDDO   

!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------
!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------
!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------


!         to be finished.......


!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------
!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------
!---------------------------- OUTPUT: SUMMARY & ECHO FILE ------------------------------





!-------------  Discretize the structure according to the division size -----------------
!  CALL SubDynDiscrt(NNode, Nodes, NElem, Elems, NProp, Props, &
!                    NReact, ReactDm, React, NInterf, InterfDm, Interf, &
!                    NCmass, CmassDm, Cmass)   

      ! Initialize index array for system K and M
!  CALL InitIAJA()

      ! Assemble system stiffness and mass matrices
!  CALL AssembleKM()

      ! Assemble system force vector
!  CALL AssembleF()      

      ! Apply constraints to stiffness and mass matrices
!  CALL ApplyConstr()

      ! Solve static problem
!  CALL SolveStatic()

      ! Solve dynamic problem
!  CALL SolveDynamic()



END SUBROUTINE BEAMFEM
!-----------------------------------------------------------------------------
SUBROUTINE SubDynDiscrt()
   USE NWTC_Library
   IMPLICIT NONE


END SUBROUTINE SubDynDiscrt
