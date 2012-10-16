MODULE IO

USE GlobalDataFun
IMPLICIT NONE

PRIVATE

PUBLIC Input, WriteError, Output

! variables needed by the main program
!----------------------------------------------
PUBLIC error,nkp,coord,member,material,curvature, &
	 & EIN,analysis_flag,&
	 & nmemb,nmate,ncurv,nframe,frame

! files needed/generated 
!---------------------------------------------
INTEGER,PARAMETER,PRIVATE:: CHAR_LEN=256
INTEGER,PARAMETER:: IN  =10     ! input file: inp_name
CHARACTER(CHAR_LEN)    :: inp_name 

INTEGER,PARAMETER:: EIN =20     ! file for echoing the inputs:
inp_name.ech
CHARACTER(CHAR_LEN+3)    :: ech_name 

INTEGER,PARAMETER:: OUT =40     ! file for output: inp_name.out
CHARACTER(CHAR_LEN+3)    :: out_name 
!--------------------------------------------------------------------------------
INTEGER,PARAMETER:: INIT =50     ! file for initial conditions:
inp_name.ini
CHARACTER(CHAR_LEN+3)    :: init_name 
!------------------------------------------------------------------------

!Private variables
!----------------------------------------
INTEGER:: analysis_flag
INTEGER:: nkp
INTEGER:: nmemb
INTEGER:: norder
INTEGER:: ninteg
INTEGER:: ncurv
INTEGER:: nmate
INTEGER:: nframe

REAL(DBL),ALLOCATABLE:: coord(:,:)
REAL(DBL),ALLOCATABLE:: member(:,:)
REAL(DBL),ALLOCATABLE:: material(:,:,:)
REAL(DBL),ALLOCATABLE:: curvature(:,:)
REAL(DBL),ALLOCATABLE:: frame(:,:,:)


!Public character variables
!=================================
CHARACTER(300)::error         ! a character variable holding  error message
!===================================

CONTAINS
!=================================

SUBROUTINE Input

INTEGER:: i,j
INTEGER:tmp_no

!Get the input file name from the command line
!--------------------------------------
CALL GETARG(1,inp_name)
IF(TRIM(inp_name)=='') THEN
 error='Please provide an input file name, executing as GEBT
input_file_name'
 RETURN
ENDIF
IF(FileOpen(IN, inp_name, 'OLD', 'READ',error))	 RETURN

! Create a file name for echoing the input data
!--------------------------------------------------
ech_name=TRIM(inp_name) // ".ech" 
IF(FileOpen(EIN,  ech_name,'REPLACE','WRITE',error)) RETURN

! Input and echo analysis control parameters.
!---------------------------------------------------------
READ(IN,*,IOSTAT=in_stat) analysis_flag
IF(IOError('read analysis conontrol parameters',error)) RETURN


CALL TitlePrint(EIN, 'Analysis Control Parameters')
WRITE(EIN,*) "Analysis Type                = ", analysis_flag

! Input and echo mesh control parameters
!----------------------------------------------------------------
READ(IN,*,IOSTAT=in_stat)  nkp,norder,ninteg,nmate,nframe,ncurv
IF(IOError('read mesh conontrol parameters',error)) RETURN
nmemb=nkp-1

CALL TitlePrint(EIN, 'Mesh Control Parameters')
WRITE(EIN,*) "Number of Key Points = ", nkp
WRITE(EIN,*) "Number of Orders = ", norder
WRITE(EIN,*) "Number of Integration Points    = ", ninteg
WRITE(EIN,*) "Number of Cross-sections         = ", nmate
WRITE(EIN,*) "Number of Frames                 = ", nframe
WRITE(EIN,*) "Number of Curvatures/Twist       = ", ncurv
WRITE(EIN,*) "Number of Members                = ", nmemb


! A small check for the control data
!---------------------------------
IF(nkp<=0 .OR. norder<=0 .OR. ninteg<=0.OR.nmate<=0) THEN
   error='nkp, norder, ninte, nmate must be greater than 0'
   RETURN
ENDIF

! Input and echo coordinates for key points
!=========================================
ALLOCATE(coord(nkp,NDIM),STAT=allo_stat)
IF(MemoryError('coord',error)) GOTO 9999
coord=0.0D0

CALL TitlePrint(EIN, 'Key Point Coordinates')

DO i=1,nkp
   READ(IN,*,IOSTAT=in_stat)tmp_no,coord(tmp_no,:)
   IF(IOError('read key point coordinates',error)) GOTO 9999

   WRITE(EIN,*)"Point NO: ",tmp_no
   WRITE(EIN,*) '--------------------------------'	
   CALL WriteVec(EIN,coord(tmp_no,:))
ENDDO

! Input and echo member properties: we need 6 numbers to describe a
!member including starting and ending points, structural/inertial
!property set #, # of norder, # of ninte, and geometry property set #
!=================================================
ALLOCATE(member(nmemb,MEMB_CONST),STAT=allo_stat)
IF(MemoryError('member',error)) GOTO 9999
member=0

CALL TitlePrint(EIN, 'Member Definition')
WRITE(EIN,*) '      KP1        KP2     Sec#  #norder  #ninteg  Geom#'
WRITE(EIN,*)
'----------------------------------------------------------------'	

DO i=1,nmemb  
   READ(IN,*,IOSTAT=in_stat)tmp_no,member(tmp_no,:)
   IF(IOError('read member properties',error)) GOTO 9999
   WRITE(EIN,*)"Member NO: ",tmp_no
   WRITE(EIN,*) '--------------------------------'	
   CALL WriteVec(EIN,member(tmp_no,:))
ENDDO

! Input and echo cross-sectional properties including flexibility
!matrix and mass matrix
!====================================================================================
IF(analysis_flag==0) THEN
	ALLOCATE(material(nmate,NDOF_ND,NDOF_ND),STAT=allo_stat)
	IF(MemoryError('material',error)) GOTO 9999
ELSE
	ALLOCATE(material(nmate,NDOF_ND+NDOF_ND,NDOF_ND),STAT=allo_stat)
	IF(MemoryError('material',error)) GOTO 9999
ENDIF
material=0.0D0

CALL TitlePrint(EIN, 'Sectional Properties')

DO i=1,nmate
   
   READ(IN,*,IOSTAT=in_stat) tmp_no 
   IF(IOError('read number of sectional properties',error)) GOTO 9999

   WRITE(EIN,*) 'Section No.         =', tmp_no
   WRITE(EIN,*) '--------------------------------'
   CALL TitlePrint(EIN, 'Sectional Flexibility Matrix')
   
   DO j=1,NDOF_ND
      READ(IN,*,IOSTAT=in_stat)material(tmp_no,j, :) 
	  IF(IOError('read sectional flexibility matrix',error)) GOTO 9999
	  CALL WriteVec(EIN,material(tmp_no,j,:))
   ENDDO

   IF(analysis_flag/=0) THEN
   		CALL TitlePrint(EIN, 'Sectional Mass Matrix')

		DO j=NDOF_ND+1,NDOF_ND+NDOF_ND
			READ(IN,*,IOSTAT=in_stat)material(tmp_no,j, :)  
			IF(IOError('read sectional mass matrix',error)) GOTO 9999
			CALL WriteVec(EIN,material(tmp_no,j,:))
		ENDDO
    ENDIF
ENDDO  

! Input and echo frame 
!=============================================================================================
IF(nframe>0) THEN
	ALLOCATE(frame(nframe,NDIM,NDIM),STAT=allo_stat)
	IF(MemoryError('frame',error)) GOTO 9999
	frame=0.0D0

	CALL TitlePrint(EIN, 'Member Frames')
	
	DO i=1,nframe
   
		READ(IN,*,IOSTAT=in_stat) tmp_no 
		IF(IOError('read number of member frames',error)) GOTO 9999
   
		WRITE(EIN,*) 'Frame No.         =', tmp_no
		WRITE(EIN,*) '-----------------------------------'	
		
		DO j=1,NDIM
			READ(IN,*,IOSTAT=in_stat)frame(tmp_no,j, :) 
			IF(IOError('read member frames',error)) GOTO 9999
			CALL WriteVec(EIN,frame(tmp_no,j,:)) ! echo
		ENDDO

	ENDDO       
ENDIF

! Input and echo initial curvatures/twist
!=========================================================
IF (ncurv>0) THEN
	ALLOCATE(curvature(ncurv,NDIM),STAT=allo_stat)
	IF(MemoryError('curvature',error)) GOTO 9999
	curvature=0.0D0

    CALL TitlePrint(EIN, 'Initial Curvatures/Twist')
	
	DO i=1,ncurv
		
		READ(IN,*,IOSTAT=in_stat) tmp_no
		IF(IOError('read # for initial curvatures/twist',error)) GOTO 9999
		 
	    READ(IN,*,IOSTAT=in_stat)curvature(tmp_no, :) 
	    IF(IOError('read initial curvatures/twist',error)) GOTO 9999
	
		WRITE(EIN,*) 'Case No.=', tmp_no
		WRITE(EIN,*) '---------------------------------------'	
		CALL WriteVec(EIN,curvature(tmp_no,:))
   ENDDO
ENDIF       

9999 IF(error/='')THEN
!	    IF(ALLOCATED(sol_mb))		 DEALLOCATE(sol_mb)
!    	IF(ALLOCATED(sol_pt))	     DEALLOCATE(sol_pt)
		IF(ALLOCATED(curvature))	 DEALLOCATE(curvature)
    	IF(ALLOCATED(material))		 DEALLOCATE(material)
		IF(ALLOCATED(member))		 DEALLOCATE(member)
		IF(ALLOCATED(coord)) 		 DEALLOCATE(coord)
		IF(ALLOCATED(frame))		 DEALLOCATE(frame)
     ENDIF
     
END SUBROUTINE
