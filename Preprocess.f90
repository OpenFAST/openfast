!***************************************************************
!*                                                             *
!* This module preprocess the finite element model including   *
!* connectivity and member information. This information are   *
!* time step independent                                        *
!***************************************************************
MODULE PreproModule

USE InternalData

IMPLICIT NONE

PRIVATE       ! So everything is private, except declared by PUBLIC
PUBLIC Preprocess 

!=============================================
CONTAINS


SUBROUTINE Preprocess(nkp,member,material,frame,coord,curvature,&
                    & memb_info,error)

IMPLICIT NONE

INTEGER,INTENT(IN)   ::nkp
INTEGER,INTENT(IN)   ::member(:,:)
REAL(DBL),INTENT(INOUT)::material(:,:,:)
REAL(DBL),INTENT(IN)   :: frame(:,:,:),coord(:,:),curvature(:,:)

TYPE (MemberInf),INTENT(OUT)::memb_info(:)

CHARACTER(*),INTENT(OUT)::error

INTEGER:: i           

INTEGER:: nmemb


nmemb=SIZE(member,1)

DO i=1,nmemb
   
   CALL MemberProperties(i,member,material,coord,curvature,&
                       & memb_info(i),error)

ENDDO


END SUBROUTINE Preprocess
!***********************************************************



!*************************************************************
!*                                                           *   
!* Extract member properties for each division               *
!*                                                           *
!*===========================================================*
!* Inputs:                                                   *
!* Output	                                                 *
!*************************************************************
SUBROUTINE MemberProperties(memb_no,member,material,frame,coord,curvature,&
                          & memb_info_i,error)

IMPLICIT NONE

INTEGER,INTENT(IN)  ::memb_no,member(:,:)
REAL(DBL),INTENT(IN)::material(:,:,:),frame(:,:,:),coord(:,:),curvature(:,:)

TYPE (MemberInf),INTENT(OUT)::memb_info_i
CHARACTER(*),INTENT(OUT)::error

INTEGER:: j,tmpN,mat_no,frame_no,curve_no
REAL(DBL)::mCab(3,3),mCoord(2,3),mL,dL
REAL(DBL)::mCurv(3),kn
REAL(DBL)::k12,kn2,kn4,kkTkn2(3,3),tmp33_1(3,3),tmp33_2(3,3),knx1


!-------------------------------
! Member sectional properties
!-------------------------------	
tmpN=NDOF_ND
IF(anlysis_flag/=0) tmpN=2*NDOF_ND

ALLOCATE(memb_info_i%mate(tmpN,NDOF_ND),STAT=allo_stat)
IF(MemoryError('memb_info_i%mate',error)) GOTO 9999

mat_no=member(memb_no,3)
memb_info_i%mate=material(mat_no,:,:)
!---------------------------------------------------

CALL MemberLength(memb_no,mCoord,member,curvature,mL,error)
memb_info_i%Length=mL




!------------------------------------------------
! Length of the member, divisions, and ending
! arc length of the division
!------------------------------------------------
mCoord(1,:)=coord(member(memb_no,1),:) ! starting point of this member
mCoord(2,:)=coord(member(memb_no,2),:) ! ending point of this member
xyz_pt1=coord(member(1,1),:)       ! starting point of the first member





dL=mL/ndiv
memb_info_i%dL=dL

ALLOCATE(memb_info_i%Le(ndiv),STAT=allo_stat)
IF(MemoryError('memb_info_i%Le',error)) GOTO 9999

DO j=1,ndiv
	memb_info_i%Le(j)=j*dL
ENDDO
!-----------------------------------------------------


!------------------------------------------------
! Frame b at the middle point of each division
! and coordinate at the middle point of each division
!------------------------------------------------
frame_no=member(memb_no,5)
IF(frame_no/=0) THEN 
	mCab =frame(frame_no,:,:)
ELSE 
	mCab=I3 ! default global frame
ENDIF

ALLOCATE(memb_info_i%triad(ndiv,3,3),STAT=allo_stat)
IF(MemoryError('memb_info_i%triad',error)) GOTO 9999

ALLOCATE(memb_info_i%coordinate(ndiv,3),STAT=allo_stat)
IF(MemoryError('memb_info_i%coordinate',error)) GOTO 9999

IF(curve_no/=0) THEN
    kkTkn2=OuterProduct(mCurv,mCurv)/(kn*kn)
	tmp33_1=I3-kkTkn2
	tmp33_2=Tilde(mCurv)/kn
	DO j=1,ndiv
		knx1=kn*(memb_info_i%Le(j)-0.5d0*dL)
		memb_info_i%triad(j,:,:)=MATMUL(mCab,tmp33_1*Cos(knx1)+tmp33_2*Sin(knx1)+kkTkn2)
		memb_info_i%coordinate(j,:)=mCoord(1,:)+MATMUL(mCab, MATMUL(tmp33_1*Sin(knx1) & 
		                        &  +tmp33_2*(1.0D0-Cos(knx1))+kkTkn2*knx1,e1/kn))
	ENDDO
ELSE
	DO j=1,ndiv
		memb_info_i%triad(j,:,:)=mCab
		memb_info_i%coordinate(j,:)=mCoord(1,:)+(j-0.5D0)/ndiv*(mCoord(2,:)-mCoord(1,:) )
	ENDDO
ENDIF
!------------------------------------------------------

9999 IF(error/='') THEN
		  IF(ALLOCATED(memb_info_i%coordinate)) DEALLOCATE(memb_info_i%coordinate)
		  IF(ALLOCATED(memb_info_i%triad)) DEALLOCATE(memb_info_i%triad)
          IF(ALLOCATED(memb_info_i%Le)) DEALLOCATE(memb_info_i%Le)
		  IF(ALLOCATED(memb_info_i%mate)) DEALLOCATE(memb_info_i%mate)
	 ENDIF

END SUBROUTINE MemberProperties
!*************************************************************

SUBROUTINE MemberLength(memb_no,mCoord,member,curvature,mL,error)

IMPLICIT NONE

REAL(DBL),INTENT(IN)  :: mCoord(:,:),curvature(:,:)
INTEGER, INTENT(IN)   :: member(:,:),memb_no
REAL(DBL),INTENT(OUT) :: mL
CHARACTER(*),INTENT(OUT)::error

INTEGER :: curve_no
REAL(DBL) :: mCurve(3), kn,k12,kn2,kn4

mL=0.0D0
mL=Norm(mCoord(2,:)-mCoord(1,:))  !length of the member for prismatic beams, the distance between the end points: |r_a-r_0|

curve_no=member(memb_no,6)
IF(curve_no/=0) THEN
	mCurv=curvature(curve_no,:)
	kn=Norm(mCurv)
	IF(kn/=0.0d0) THEN
	   IF(mCurv(1)==0.0D0) THEN
			mL=2.0D0*ASIN(0.5D0*kn*mL)/kn
	   ELSE IF(mCurv(2)/=0.0D0.OR.mCurv(3)/=0.0D0) THEN
		   k12=mCurv(1)*mCurv(1)
		   kn2=kn*kn
		   kn4=kn2*kn2
		   mL=Rtbis(CurveBeamFun,kn,mL,kn2,k12,kn4,mL,mL*kn/mCurv(1),TOLERANCE*100,100,error)
	   ENDIF
     ENDIF
ENDIF

END SUBROUTINE MemberLength


!************************************************************
! Function for evaluating arc length of initially curved and twisted beams
!****************************************************************
FUNCTION  CurveBeamFun(kn,mL,kn2,k12,kn4,xvar)
IMPLICIT NONE
REAL(DBL),INTENT(IN)::kn,mL,kn2,k12,kn4
REAL(DBL),INTENT(IN)::xvar
REAL(DBL)			::CurveBeamFun
REAL(DBL)           ::knx
	
knx=kn*xvar
CurveBeamFun=mL*mL-( 2.0D0*(kn2-k12)*(1.0D0-COS(knx))+k12*knx*knx )/kn4

END FUNCTION  CurveBeamFun
!**************************************************************************





!*************************************************************
!*                                                           *   
!*  Use biosection to find root of a function                *
!* from the book of Numerical Recipes                        *
!*                                                           *
!*===========================================================*
!* Inputs:                                                   *
!*  func     --	  the given nonlinear equation func==0.0d0   *
!*  x1, x2   --   the root is known to lie between x1, x2    *
!*  xacc     --   a small number indicating the accuracy     *
!*  maxit    --   max number of iterations                   *
!* Output	                                                 *
!*  root     --  the root                                    *
!*************************************************************
FUNCTION Rtbis(func,kn,mL,kn2,k12,kn4,x1,x2,xacc,maxit,error) RESULT(root)

IMPLICIT NONE
REAL(DBL),INTENT(IN)::kn,mL,kn2,k12,kn4
REAL(DBL),INTENT(IN)	::x1,x2,xacc
INTEGER,INTENT(IN)		::maxit
CHARACTER(*),INTENT(OUT)::error
REAL(DBL)			    ::root

INTERFACE
    FUNCTION func(kn,mL,kn2,k12,kn4,x)
	    USE GlobalDataFun
		IMPLICIT NONE
		REAL(DBL),INTENT(IN)::kn,mL,kn2,k12,kn4
		REAL(DBL),INTENT(IN)::x
		REAL(DBL)           ::func
	END FUNCTION func
END INTERFACE

INTEGER  :: j
REAL(DBL):: dx,f,fmid,xmid

fmid=func(kn,mL,kn2,k12,kn4,x2)
f=func(kn,mL,kn2,k12,kn4,x1)
IF(f*fmid>=0.0D0) error='rtbis: root is not bracketed.'
IF(f<0.0D0)THEN
	root=x1
	dx=x2-x1
ELSE
	root=x2
	dx=x1-x2
ENDIF
DO j=1,maxit
	dx=dx*0.5D0
	xmid=root+dx
	fmid=func(kn,mL,kn2,k12,kn4,xmid)
	IF(fmid<=0.0D0)root=xmid
	IF(ABS(dx)<xacc.OR.fmid==0.0D0)RETURN
ENDDO

error='Rtbis: the number of bisection exceeds the allowed.'

END FUNCTION Rtbis
!**********************************************************


END MODULE PreproModule
!*********************************************
