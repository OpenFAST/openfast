MODULE InternalData

USE GlobalDataFun


IMPLICIT NONE

LOGICAL,PARAMETER:: DEBUG=.FALSE. ! flag for debugging purpose

!File needed internally for debugging purpose
!===================================================
INTEGER,PARAMETER:: IOUT=30     ! file for debugging: inp_name.deb
CHARACTER(64)    :: deb_name 


!Global variables
!============================================================================
TYPE MemberInf
!	 REAL(DBL) :: vertex(2,3)
!	 REAL(DBL) :: curv(3)
	 REAL(DBL) :: mL                         ! length of member 
	 REAL(DBL),ALLOCATABLE::mate(:,:)        ! mate(:,6): flexibity and mass properties for each division
	 REAL(DBL),ALLOCATABLE::coord_node(:,:)  !coord_node(node,3) coordinate for the nodes
END TYPE MemberInf
	
END MODULE InternalData
!============================================================================

