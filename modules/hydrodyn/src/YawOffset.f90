MODULE YawOffset

USE NWTC_LIBRARY

IMPLICIT NONE

INTEGER(IntKi), PARAMETER        :: i2h = 1_IntKi
INTEGER(IntKi), PARAMETER        :: h2i = 2_IntKi

INTERFACE hiFrameTransform
   MODULE PROCEDURE hiFrameTransformVec3
   MODULE PROCEDURE hiFrameTransformMat
END INTERFACE hiFrameTransform

INTERFACE GetRotAngs
   MODULE PROCEDURE GetRotAngsR
   MODULE PROCEDURE GetRotAngsD
END INTERFACE GetRotAngs

CONTAINS

SUBROUTINE GetPtfmRefYOrient(PtfmRefY, Orient, ErrStat, ErrMsg)

   REAL(ReKi),     INTENT(IN   ) :: PtfmRefY
   REAL(ReKi),     INTENT(  OUT) :: Orient(3,3)
   INTEGER(IntKi), INTENT(  OUT) :: ErrStat
   CHARACTER(*),   INTENT(  OUT) :: ErrMsg

   REAL(ReKi)                    :: cosRefY
   REAL(ReKi)                    :: sinRefY

   ErrStat = ErrID_None
   ErrMsg  = ''

   call Eye(Orient, ErrStat, ErrMsg)
   cosRefY = cos(PtfmRefY)
   sinRefY = sin(PtfmRefY)
   Orient(1,1) =  cosRefY
   Orient(1,2) =  sinRefY
   Orient(2,1) = -sinRefY
   Orient(2,2) =  cosRefY

END SUBROUTINE GetPtfmRefYOrient

SUBROUTINE RotTrans(RotationType,PtfmRefY,Rotation,Orientation,ErrTxt,ErrStat,ErrMsg)
   ! Compute the orientation matrix with potentially large reference yaw offset
   ! This subroutine essentially extends SmllRotTrans to accommodate a large yaw offset
   CHARACTER(*),   INTENT(IN   ) :: RotationType
   REAL(ReKi),     INTENT(IN   ) :: PtfmRefY
   REAL(R8Ki),     INTENT(IN   ) :: Rotation(3)
   REAL(R8Ki),     INTENT(  OUT) :: Orientation(3,3)
   CHARACTER(*),   INTENT(IN   ) :: ErrTxt
   INTEGER(IntKi), INTENT(  OUT) :: ErrStat
   CHARACTER(*),   INTENT(  OUT) :: ErrMsg

   REAL(ReKi)                    :: PtfmRefYOrient(3,3)
   REAL(R8Ki)                    :: SmllOMat(3,3)
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   CHARACTER(*), PARAMETER       :: RoutineName = 'RotTrans'

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Orientation matrix associated with large reference yaw offset
   call GetPtfmRefYOrient(PtfmRefY, PtfmRefYOrient, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
        
   ! Orientation matrix for the remaining small rotation from the large reference yaw offset
   call SmllRotTrans( RotationType, Rotation(1), Rotation(2), Rotation(3)-REAL(PtfmRefY,R8Ki), SmllOMat, ErrTxt, ErrStat2, ErrMsg2 )
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   ! Combine the contributions
   Orientation = matmul(SmllOMat,PtfmRefYOrient)

END SUBROUTINE RotTrans


FUNCTION GetRotAngsR(PtfmRefY, DCMat, ErrStat, ErrMsg)
   ! Compute the intrinsic Tait-Bryan angles (yaw first, pitch second, roll last) based on large yaw offset
   ! The subroutine essentially extends GetSmllRotAngs to accommodate a large yaw offset
   REAL(ReKi),     INTENT(IN   ) :: PtfmRefY
   REAL(SiKi),     INTENT(IN   ) :: DCMat(3,3)
   INTEGER(IntKi), INTENT(  OUT) :: ErrStat
   CHARACTER(*),   INTENT(  OUT) :: ErrMsg

   REAL(SiKi)                    :: GetRotAngsR( 3 )

   REAL(ReKi)                    :: PtfmRefYOrient(3,3)
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   CHARACTER(*), PARAMETER       :: RoutineName = 'GetRotAngs'

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Orientation matrix associated with large reference yaw offset
   call GetPtfmRefYOrient(PtfmRefY, PtfmRefYOrient, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   GetRotAngsR = GetSmllRotAngsR ( matmul(DCMat,transpose(REAL(PtfmRefYOrient,SiKi))), ErrStat2, ErrMsg2 )
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   GetRotAngsR(3) = GetRotAngsR(3) + PtfmRefY   

END FUNCTION GetRotAngsR

FUNCTION GetRotAngsD(PtfmRefY, DCMat, ErrStat, ErrMsg)
   ! Compute the intrinsic Tait-Bryan angles (yaw first, pitch second, roll last) based on large yaw offset
   ! The subroutine essentially extends GetSmllRotAngs to accommodate a large yaw offset
   REAL(ReKi),     INTENT(IN   ) :: PtfmRefY
   REAL(DbKi),     INTENT(IN   ) :: DCMat(3,3)
   INTEGER(IntKi), INTENT(  OUT) :: ErrStat
   CHARACTER(*),   INTENT(  OUT) :: ErrMsg

   REAL(DbKi)                    :: GetRotAngsD( 3 )

   REAL(ReKi)                    :: PtfmRefYOrient(3,3)
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2
   CHARACTER(*), PARAMETER       :: RoutineName = 'GetRotAngs'

   ErrStat = ErrID_None
   ErrMsg  = ''

   ! Orientation matrix associated with large reference yaw offset
   call GetPtfmRefYOrient(PtfmRefY, PtfmRefYOrient, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   GetRotAngsD = GetSmllRotAngsD ( matmul(DCMat,transpose(REAL(PtfmRefYOrient,DbKi))), ErrStat2, ErrMsg2 )
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   GetRotAngsD(3) = GetRotAngsD(3) + PtfmRefY   

END FUNCTION GetRotAngsD


SUBROUTINE hiFrameTransformVec3(Mode,PtfmRefY,VecIn,VecOut,ErrStat,ErrMsg)
   INTEGER(IntKi), INTENT(IN   ) :: Mode
   REAL(ReKi),     INTENT(IN   ) :: PtfmRefY
   REAL(ReKi),     INTENT(IN   ) :: VecIn(3)
   REAL(ReKi),     INTENT(  OUT) :: VecOut(3)
   INTEGER(IntKi), INTENT(  OUT) :: ErrStat
   CHARACTER(*),   INTENT(  OUT) :: ErrMsg

   REAL(ReKi)                    :: PtfmRefYOrient(3,3)
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2

   CHARACTER(*),   PARAMETER     :: RoutineName = 'hiFrameTransformVec3'

   ErrStat = ErrID_None
   ErrMsg  = ''

   call GetPtfmRefYOrient(PtfmRefY, PtfmRefYOrient, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   if (Mode .EQ. i2h) then   ! i-frame to h-frame
      VecOut = matmul(PtfmRefYOrient,VecIn)
   else if (Mode .EQ. h2i) then  ! h-frame to i-frame
      VecOut = matmul(transpose(PtfmRefYOrient),VecIn)
   else
      call SetErrStat(ErrID_Fatal, "Mode must be 1 or 2", ErrStat, ErrMsg, RoutineName) 
   end if

END SUBROUTINE hiFrameTransformVec3

SUBROUTINE hiFrameTransformMat(Mode,PtfmRefY,MatIn,MatOut,ErrStat,ErrMsg)
   INTEGER(IntKi), INTENT(IN   ) :: Mode
   REAL(ReKi),     INTENT(IN   ) :: PtfmRefY
   REAL(ReKi),     INTENT(IN   ) :: MatIn(3,3)
   REAL(ReKi),     INTENT(  OUT) :: MatOut(3,3)
   INTEGER(IntKi), INTENT(  OUT) :: ErrStat
   CHARACTER(*),   INTENT(  OUT) :: ErrMsg

   REAL(ReKi)                    :: PtfmRefYOrient(3,3)
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(ErrMsgLen)          :: ErrMsg2

   CHARACTER(*), PARAMETER       :: RoutineName = 'hiFrameTransformMat'

   ErrStat = ErrID_None
   ErrMsg  = ''

   call GetPtfmRefYOrient(PtfmRefY, PtfmRefYOrient, ErrStat2, ErrMsg2)
   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   if (Mode .EQ. i2h) then   ! i-frame to h-frame
      MatOut = matmul(matmul(PtfmRefYOrient,MatIn),transpose(PtfmRefYOrient))
   else if (Mode .EQ. h2i) then  ! h-frame to i-frame
      MatOut = matmul(matmul(transpose(PtfmRefYOrient),MatIn),PtfmRefYOrient)
   else
      call SetErrStat(ErrID_Fatal, "Mode must be 1 or 2", ErrStat, ErrMsg, RoutineName) 
   end if   

END SUBROUTINE hiFrameTransformMat

END MODULE YawOffset