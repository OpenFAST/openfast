!---------------------------------------------------------------
!
! This module contains general-purpose global constants,
! I/O functions/subroutines and math functions/subroutines
!
! Gobal constants/variables
!----------------------------
! DBL,PI,DEG_2_RAD,RAD_2_DEG,TOLERANCE, I3,in_stat,allo_stat
!
! I/O functions/subroutines
!------------------------------
! FUNCTION FileOpen (file_unit,file_name,sta_type,rw_type,error)
! FUNCTION IOError(message,error)
! FUNCTION ItoChar(n)
! FUNCTION LowerCase(string) 
! FUNCTION  MemoryError(allo_stat,vari_name,error)
! SUBROUTINE TitlePrint(file_unit, title)
! FUNCTION UpperCase(string) 
! SUBROUTINE WriteError(EIN,error)
! SUBROUTINE WriteVec(file_unit,vec)
!
! math functions/subroutines
!------------------------------
! SUBROUTINE Invert(matrix_in,matrix,vari_name,error)
! FUNCTION Norm(vector)
! FUNCTION OuterProduct(vec1, vec2)
! FUNCTION Tilde(vect)
!---------------------------------------------------------------

MODULE GlobalDataFun

IMPLICIT NONE

PRIVATE ! everything is private except declared by PUBLIC

PUBLIC DBL,PI,DEG_2_RAD,RAD_2_DEG,TOLERANCE,I3,FMT_INT,FMT_REAL,NDIM,NDOF_ND,NSTRN,e1
PUBLIC NORDER,NINTEG
PUBLIC in_stat, allo_stat,MEMB_CONST

PUBLIC FileOpen,IOError,MemoryError,TitlePrint,WriteError,WriteVec
PUBLIC Invert,MATMUL3,Norm,OuterProduct,Tilde,CrossProduct
INTERFACE WriteVec
 	  MODULE PROCEDURE WriteIntVector,WriteRealVector ! Write a vector
END INTERFACE


!=============================================================================
!
! Global constants/variables
!
!=============================================================================
INTEGER,PARAMETER:: NDIM=3         ! All the beams could behavior in the 3D space 
INTEGER,PARAMETER:: NDOF_ND=6      ! degrees of freedom per node is 6.
INTEGER,PARAMETER:: MEMB_CONST=3   ! Number of labels needed for member properties

INTEGER,PARAMETER:: NORDER = 2
INTEGER,PARAMETER:: NINTEG = 2

INTEGER, PARAMETER:: DBL=SELECTED_REAL_KIND(15,307)

REAL(DBL),PARAMETER:: PI        =    3.1415926535897932D0    
REAL(DBL),PARAMETER:: DEG_2_RAD =    1.7453292519943296D-2 ! the ratio between radians and degrees
REAL(DBL),PARAMETER:: RAD_2_DEG =    5.7295779513082321D1  ! convert radian to degree
REAL(DBL),PARAMETER:: TOLERANCE=EPSILON(1.0_DBL)           ! a smart number of the double precision real number

REAL(DBL),PARAMETER::I3(3,3) = RESHAPE((/1.D0, 0.D0, 0.D0,& 
                                       & 0.D0, 1.D0, 0.D0,&
                                       & 0.D0, 0.D0, 1.D0/),(/3,3/))  ! The 3x3 identity matrix
REAL(DBL),PARAMETER::e1(3)=(/1._DBL,0._DBL,0._DBL/)  ! The e1 unit vector

INTEGER:: in_stat      ! flag to indicate if the I/O process is successful: if positive, an error occured;
                       ! if negative, an end-of-file or end-of-record condition occurred; 
	               ! zero, no error, end-of-file, or end-of-record condition occurred.  
INTEGER:: allo_stat    ! flag to indicate status of allocating memory 


CHARACTER(*),PARAMETER :: FMT_REAL='ES15.7' ! format for output real numbers
CHARACTER(*),PARAMETER :: FMT_INT='I8'       ! format for output integer numbers
!=========================================================================================


CONTAINS
!=================================================================
!
! The following are general purpose I/O functions/subroutines
!
!=================================================================

!************************************************************
!*                                                          *
!*    To open an old or new file for reading or writing     *
!*															*
!************************************************************
FUNCTION 	FileOpen (file_unit,file_name,sta_type,rw_type,error)

LOGICAL                   ::FileOpen
INTEGER,INTENT(IN)        ::file_unit 
CHARACTER(*),INTENT(IN)   ::file_name
CHARACTER(*),INTENT(IN)   ::sta_type
CHARACTER(*),INTENT(IN)   ::rw_type
CHARACTER(*),INTENT(OUT)  ::error

error=''
FileOpen=.FALSE.

OPEN (UNIT=file_unit, file=file_name,STATUS=sta_type,ACTION = rw_type,IOSTAT=in_stat)

IF (in_stat/=0) THEN
  IF(rw_type=='READ') error='Cannot open the file '//TRIM(file_name)//' for reading!'
  IF(rw_type=='WRITE')error='Cannot open the file '//TRIM(file_name)//' for writing!'
ENDIF

IF(error/='')FileOpen=.TRUE.

END FUNCTION FileOpen
!***********************************************************



!************************************************************
!*                                                          *
!*        Check the error of I/O processing                 *
!*															*
!************************************************************
FUNCTION  IOError(message,error)

LOGICAL                 ::IOError
CHARACTER(*),INTENT(IN) ::message        ! a character variable to hold error message
CHARACTER(*),INTENT(OUT)::error

error=''
IOError=.FALSE.

IF(in_stat/=0) THEN 
    error='I/O error: '//TRIM(message)
    IOError=.TRUE.
ENDIF

END FUNCTION IOError
!***********************************************************



!************************************************************
!*                                                          *
!*         Convert an integer to character                  *
!*                                                          *
!************************************************************
FUNCTION ItoChar(n) RESULT(char)

	INTEGER,INTENT(IN):: n
	CHARACTER(20):: char

	 WRITE(char, *) n

END FUNCTION ItoChar
!***********************************************************



!************************************************************
!*                                                          *
!*      Convert a string or character to lower case         *
!*                                                          *
!************************************************************
FUNCTION LowerCase(string) RESULT(lc_string)

CHARACTER(*),INTENT(IN):: string
CHARACTER(LEN=LEN(string)):: lc_string

CHARACTER(LEN=26),PARAMETER:: UPPER='ABCDEFGHIJKLMNOPQRSTUVWXYZ',&
	                          lower='abcdefghijklmnopqrstuvwxyz'

INTEGER:: k  ! loop counter
INTEGER::loc ! position in alphabet

lc_string=string
    
DO k=1,len(string)
	loc=INDEX(UPPER,string(k:k))
	IF(loc/=0)lc_string(k:k)=lower(loc:loc)
ENDDO

END FUNCTION LowerCase
!************************************************************



!************************************************************
!*                                                          *
!*        Check the error of memory allocation              *
!*															*
!************************************************************
FUNCTION  MemoryError(vari_name,error)

LOGICAL                 ::MemoryError
CHARACTER(*),INTENT(IN) ::vari_name         ! a character variable to hold variable name
CHARACTER(*),INTENT(OUT)  ::error

error=''
MemoryError=.FALSE.

IF(allo_stat/=0) THEN
	error='Memory error: allocate '//TRIM(vari_name)
    MemoryError=.TRUE.
ENDIF


END FUNCTION MemoryError
!************************************************************



!************************************************************
!*                                                          *
!*        To print a title for a block of data              *
!*															*
!************************************************************
SUBROUTINE TitlePrint(file_unit, title)

INTEGER,INTENT(IN)        ::file_unit 
CHARACTER(*),INTENT(IN)   ::title

WRITE(file_unit,*) 
WRITE(file_unit,*) 
WRITE(file_unit,'(1x,100A)') title
WRITE(file_unit,*)'========================================================'

END SUBROUTINE TitlePrint
!***********************************************************



!************************************************************
!*                                                          *
!*    Convert a string or character to upper case           *
!*                                                          *
!************************************************************
FUNCTION UpperCase(string) RESULT(uc_string)

CHARACTER(*),INTENT(IN):: string
CHARACTER(LEN=LEN(string)):: uc_string

CHARACTER(LEN=26),PARAMETER:: UPPER='ABCDEFGHIJKLMNOPQRSTUVWXYZ',&
		                      lower='abcdefghijklmnopqrstuvwxyz'

INTEGER:: k  ! loop counter
INTEGER::loc ! position in alphabet

uc_string=string
    
DO k=1,len(string)
	loc=INDEX(lower,string(k:k))
	IF(loc/=0)uc_string(k:k)=UPPER(loc:loc)
ENDDO

END FUNCTION UpperCase
!************************************************************



!************************************************************
!*                                                          *
!*    Write error to the echo file                          *
!*															*
!************************************************************
SUBROUTINE WriteError(EIN,error)

INTEGER,INTENT(IN)::EIN ! file unit to write the error message
CHARACTER(*),INTENT(IN)  ::error

LOGICAL file_opened

INQUIRE (EIN,  OPENED = file_opened) ! Check whether the file is already opened, if yes, then dump the error message to this file

IF(file_opened)THEN

	WRITE(EIN,*) 
	IF(error/='')THEN

		CALL TitlePrint(EIN, 'Error Message')
		WRITE(EIN,'(1x, 300A)') error 

	ELSE

		WRITE(EIN,*) 'Congratulations! No errors!'
	ENDIF

	CLOSE(EIN)
ENDIF

END SUBROUTINE WriteError
!********************************************************



!************************************************************
!*                                                          *
!*  Write an integer vector to the file_unit                *
!*															*
!************************************************************
SUBROUTINE WriteIntVector(file_unit,vec)

INTEGER,INTENT(IN)::file_unit
INTEGER,INTENT(IN)::vec(:)

WRITE(file_unit,'(1x,'//TRIM(ItoChar(SIZE(vec)))//FMT_INT//')')vec

END SUBROUTINE WriteIntVector
!********************************************************



!************************************************************
!*                                                          *
!*    Write a real vector to the file_unit                  *
!*															*
!************************************************************
SUBROUTINE WriteRealVector(file_unit,vec)

INTEGER,INTENT(IN)::file_unit
REAL(DBL),INTENT(IN)::vec(:)
!REAL(DBL)::tmp(SIZE(vec)),vec_norm

!tmp=vec
!vec_norm=Norm(vec)

!IF(vec_norm>TOLERANCE) THEN
!	WHERE(ABS(tmp/vec_norm)<TOLERANCE)tmp=0.0D0  ! not output components which are negligible comparing to the biggest term in the vector
!ELSE
!	tmp=0.0D0
!ENDIF

WRITE(file_unit,'(1x,'//TRIM(ItoChar(SIZE(vec)))//FMT_REAL//')')vec

END SUBROUTINE WriteRealVector
!********************************************************


!***********************************************
!*       Invert a small square matrix          *
!*                                             * 
!***********************************************
SUBROUTINE Invert(matrix_in,matrix,vari_name,error)
 
REAL(DBL),INTENT(IN)::matrix_in(:,:) ! the matrix to be inverted
CHARACTER(*),INTENT(IN)::vari_name

REAL(DBL),INTENT(OUT)::matrix(:,:)    ! the inverse of the matrix

CHARACTER(*),INTENT(OUT)::error

INTEGER::i,k,n 
REAL(DBL)::con,diag_sum,zero

matrix=matrix_in
  
n= SIZE(matrix,1)

diag_sum=0.0d0

DO i=1,n
	diag_sum=diag_sum+matrix(i,i)
ENDDO

zero=TOLERANCE*diag_sum/n
DO k=1,n
   con=matrix(k,k);
   
   IF(ABS(CON)<zero) THEN
      error='The matrix is not invertable '//TRIM(vari_name)
      GOTO 9999
   ENDIF

   matrix(k,k)=1.d0
   matrix(k,:)=matrix(k,:)/con
   
   DO i=1,n
      IF(i/=k) THEN
         con=matrix(i,k); matrix(i,k)=0.d0
         matrix(i,:)=matrix(i,:) - matrix(k,:)*con
      END IF
   ENDDO
   
ENDDO

9999 RETURN

END SUBROUTINE
!***********************************************************




!*************************************************************
!*                                                           *   
!*  Calculate the L2 norm of a real vector                   *
!*                                                           *
!*************************************************************
FUNCTION Norm(vector)

REAL(DBL),INTENT(IN):: vector(:)
REAL(DBL)::Norm 

Norm=SQRT(DOT_PRODUCT(vector,vector))

END FUNCTION Norm
!**********************************************************



!*************************************************************
!*                                                           *   
!*  Calculate the outer product of two vectors               *
!*                                                           *
!*************************************************************
FUNCTION OuterProduct(vec1, vec2)

REAL(DBL),INTENT(IN):: vec1(:), vec2(:)
REAL(DBL)::OuterProduct(SIZE(vec1),SIZE(vec2)) 

INTEGER:: i,j,n1,n2

n1=SIZE(vec1); n2=SIZE(vec2)

DO i=1, n1
	DO j=1, n2
		OuterProduct(i,j)=vec1(i)*vec2(j)
	ENDDO
ENDDO

END FUNCTION OuterProduct
!**********************************************************



!*************************************************************
!*                                                           *   
!*  Carry out the tilde operation for a real vector          *
!*                                                           *
!*===========================================================*
!* Input:                                                    *
!*  vect     --	  a vector with three components             *
!* Output	                                                 *
!*  tilde    --  the 3X3 antisymmetric matrix                *
!*************************************************************
FUNCTION Tilde(vect)

REAL(DBL),INTENT(IN):: vect(3)
REAL(DBL):: Tilde(3,3) 

Tilde=0.0d0

Tilde(1,2)= -vect(3)
Tilde(1,3)=  vect(2)
Tilde(2,1)=  vect(3)
Tilde(2,3)= -vect(1)
Tilde(3,1)= -vect(2)
Tilde(3,2)=  vect(1)

END FUNCTION
!**********************************************************


!*************************************************************
!*                                                           *   
!*  Carry out cross product of two real vectors              *
!*                                                           *
!*===========================================================*
!* Input:                                                    *
!*  a,b     --	  two vectors with three components          *
!* Output	                                                 *
!*  CrossProduct  --  the result vector                      *
!*************************************************************
FUNCTION CrossProduct(a,b)

REAL(DBL),INTENT(IN):: a(3), b(3)
REAL(DBL):: CrossProduct(3) 


CrossProduct(1)=a(2)*b(3)-a(3)*b(2)
CrossProduct(2)=a(3)*b(1)-a(1)*b(3)
CrossProduct(3)=a(1)*b(2)-a(2)*b(1)

END FUNCTION
!**********************************************************


END MODULE GlobalDataFun
!==============================================================





