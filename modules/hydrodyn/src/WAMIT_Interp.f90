!**********************************************************************************************************************************
!  NOTE: documentation in this file is written for use with Doxygen 1.8.6 and higher.
!
!> WAMIT_Interp module
!!
!!  This module is used for interpolating the QTFs used in the WAMIT and WAMIT2 modules of HydroDyn.
!!
!..................................................................................................................................
! LICENSING
! Copyright (C) 2014-2015  National Renewable Energy Laboratory
!
!    This file is part of HydroDyn.
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

MODULE WAMIT_Interp


   USE NWTC_Library
   IMPLICIT NONE

   PRIVATE

      ! Public subroutines
   PUBLIC   :: WAMIT_Interp2D_Cplx
   PUBLIC   :: WAMIT_Interp3D_Cplx
   PUBLIC   :: WAMIT_Interp4D_Cplx


CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine takes the complex valued QTF dataset and interpolates for the desired coordinate in (Omega,WaveDir) space.
!!
!! A few important notes concerning this subroutine:
!!    1. It is complex valued.  The values represent the second order wave force as calculated by WAMIT.
!!    2. The dimenions of DataSet2D are Frequency1 (positive valued) and Wave Direction (degrees).
!!    3. The wave direction requested might be between end points of wave direction dimension (ie. at 179 degrees when 
!!          WvDir1(1)=175, WvDir(Dims(3))=-175)
!!    4. The arrays WvFreq1 and WvDir1 will give the values for each dimension that correspond to each index of DataSet2D.
!!    5. The data is not necessarily equally spaced in any direction: ie. WvFreq1 may not have uniform spacing between points.
!!    6. If a point is requested, it can be assumed that it lies within DataSet2D (this is checked before calling this subroutine)
!!    7. LastIndex contains the index numbers for the lowest bound on the indexes used on the last interpolation.  If they are 0,
!!          then assume this is the first call to this subroutine.
!!    8. DataSet2D is complete and not sparse.  There are no missing values.
!!
SUBROUTINE WAMIT_Interp2D_Cplx( InCoord, DataSet2D, WvFreq1, WvDir1, LastIndex, OutForce, ErrStat, ErrMsg )

      ! I/O variables

   REAL(SiKi),       INTENT(IN   )     :: InCoord(2)                                   !< Arranged as (Omega1, WaveDir1)
   COMPLEX(SiKi),    INTENT(IN   )     :: DataSet2D(:,:)                               !< Arranged as Index 1= Omega1, Index 2= WaveDir1.
   REAL(SiKi),       INTENT(IN   )     :: WvFreq1(:)                                   !< Frequencies associated with Index 1 of DataSet2D
   REAL(SiKi),       INTENT(IN   )     :: WvDir1(:)                                    !< Directions associated with Index 2 of DataSet2D
   INTEGER(IntKi),   INTENT(INOUT)     :: LastIndex(2)                                 !< Index for the last (Omega1, WaveDir1) used
   COMPLEX(SiKi),    INTENT(  OUT)     :: OutForce                                     !< The interpolated resulting force from DataSet2D
   INTEGER(IntKi),   INTENT(  OUT)     :: ErrStat                                      !< Error status
   CHARACTER(*),     INTENT(  OUT)     :: ErrMsg                                       !< Error message if ErrStat /= ErrID_None


      ! Local variables

   REAL(SiKi)                          :: Coords(2)                                    !< coordinates with wave directions converted to range [-180, 180)
   INTEGER(IntKi)                      :: n(2)                                         !< number of points in WvFreq1 and WvDir1, and WvDir2
   
   INTEGER(IntKi)                      :: Indx_Lo(2)                                   !< index associated with lower bound of dimension 1,2 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   INTEGER(IntKi)                      :: Indx_Hi(2)                                   !< index associated with upper bound of dimension 1,2 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   REAL(SiKi)                          :: Pos_Lo(2)                                    !< coordinate value with lower bound of dimension 1,2 
   REAL(SiKi)                          :: Pos_Hi(2)                                    !< coordinate value with upper bound of dimension 1,2 
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""   
      
   n(1) = SIZE(WvFreq1)
   n(2) = SIZE(WvDir1)

         ! Check sizes of input:
   IF ( SIZE( DataSet2D, 1 ) /= n(1) )  CALL SetErrStat( ErrID_Fatal, 'Invalid dimensions: size(DataSet2D,1) must be the size(WvFreq1)', ErrStat, ErrMsg, 'WAMIT_Interp2D_Cplx' )
   IF ( SIZE( DataSet2D, 2 ) /= n(2) )  CALL SetErrStat( ErrID_Fatal, 'Invalid dimensions: size(DataSet2D,2) must be the size(WvDir1)',  ErrStat, ErrMsg, 'WAMIT_Interp2D_Cplx' )
   IF (ErrStat == ErrID_Fatal) RETURN
   

      ! find the indices into the arrays representing coordinates of each dimension:
      
   Coords = InCoord
      
      ! make sure these requested degrees fall in the range -180 <= Coords(2) < 180
   Coords(2) = MODULO( Coords(2), 360.0_SiKi )
   IF ( Coords(2) >= 180.0_SiKi ) Coords(2) = Coords(2) - 360.0_SiKi
   
   CALL LocateStp( Coords(1), WvFreq1, LastIndex(1), n(1) )
   CALL LocateStp( Coords(2), WvDir1,  LastIndex(2), n(2) )
   
   Indx_Lo = LastIndex  ! at this point, 0 <= Indx_Lo(i) <= n(i) for all i
   
   
   ! WvFreq1 (indx 1)
   IF (Indx_Lo(1) == 0) THEN
      Indx_Lo(1) = 1
   ELSEIF (Indx_Lo(1) == n(1) ) THEN
      Indx_Lo(1) = max( n(1) - 1, 1 )                    ! make sure it's a valid index
   END IF     
   Indx_Hi(1) = min( Indx_Lo(1) + 1 , n(1) )             ! make sure it's a valid index

   
   ! WvDir1 (indx 2)   [use modular arithmetic]
   IF (Indx_Lo(2) == 0) THEN
      Indx_Hi(2) = 1                           
      Indx_Lo(2) = n(2)
   ELSEIF (Indx_Lo(2) == n(2) ) THEN
      Indx_Hi(2) = 1      
   ELSE
      Indx_Hi(2) = min( Indx_Lo(2) + 1, n(2) )        ! make sure it's a valid index
   END IF      
      
      ! calculate the positions of all dimensions:
      
   pos_Lo(1) = WvFreq1(Indx_Lo(1))
   pos_Hi(1) = WvFreq1(Indx_Hi(1))
      
   pos_Lo(2) = WvDir1(Indx_Lo(2))
   pos_Hi(2) = WvDir1(Indx_Hi(2))
   
   
      ! angles have to be adjusted so that pos_Lo(2) <= Coords(2) <= pos_Hi(2)
   IF ( Indx_Hi(2) == 1 .AND. n(2) > 1 )  THEN ! we're looping around the array [periodic]
      IF ( pos_Lo(2) < Coords(2) ) THEN
         pos_Hi(2) = pos_Hi(2) + 360.0_SiKi
      ELSEIF ( pos_Lo(2) /= Coords(2) ) THEN !bjj: I think it's okay if we don't use equalRealNos here
         pos_Lo(2) = pos_Lo(2) - 360.0_SiKi 
      END IF
   END IF   
   

   CALL Interp2D_withIndx_Cplx( Coords, DataSet2D, Indx_Lo, Indx_Hi, pos_Lo, pos_Hi, OutForce )
      
   
END SUBROUTINE WAMIT_Interp2D_Cplx



!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine takes the complex valued QTF dataset and interpolates for the desired coordinate in (Omega,WaveDir,WaveDir) space.
!!
!! A few important notes concerning this subroutine:
!!    1. It is complex valued.  The values represent the second order wave force as calculated by WAMIT.
!!    2. The dimenions of DataSet3D are Frequency1 (positive valued), Wave Direction1 (degrees), and Wave Direction2 (degrees).
!!    3. The wave direction requested might be between end points of wave direction dimension (ie. at 179 degrees when 
!!          WvDir1(1)=175, WvDir(Dims(3))=-175)
!!    4. The arrays WvFreq1, WvDir1, and WvDir2, will give the values for each dimension that correspond to each index of DataSet3D.
!!    5. The data is not necessarily equally spaced in any direction: ie. WvFreq1 may not have uniform spacing between points.
!!    6. If a point is requested, it can be assumed that it lies within DataSet3D (this is checked before calling this subroutine)
!!    7. LastIndex contains the index numbers for the lowest bound on the indexes used on the last interpolation.  If they are 0,
!!          then assume this is the first call to this subroutine.
!!    8. DataSet4D is complete and not sparse.  There are no missing values.
!!
SUBROUTINE WAMIT_Interp3D_Cplx( InCoord, DataSet3D, WvFreq1, WvDir1, WvDir2, LastIndex, OutForce, ErrStat, ErrMsg )

      ! I/O variables

   REAL(SiKi),       INTENT(IN   )     :: InCoord(3)                                   !< Arranged as (Omega1, WaveDir1, WaveDir2)
   COMPLEX(SiKi),    INTENT(IN   )     :: DataSet3D(:,:,:)                             !< Arranged as Index 1= Omega1, Index 2= WaveDir1, Index 3= WaveDir2.
   REAL(SiKi),       INTENT(IN   )     :: WvFreq1(:)                                   !< Frequencies associated with Index 1 of DataSet3D
   REAL(SiKi),       INTENT(IN   )     :: WvDir1(:)                                    !< Directions associated with Index 2 of DataSet3D
   REAL(SiKi),       INTENT(IN   )     :: WvDir2(:)                                    !< Directions associated with Index 3 of DataSet3D
   INTEGER(IntKi),   INTENT(INOUT)     :: LastIndex(3)                                 !< Index for the last (Omega1, WaveDir1, WaveDir2) used
   COMPLEX(SiKi),    INTENT(  OUT)     :: OutForce                                     !< The interpolated resulting force from DataSet3D
   INTEGER(IntKi),   INTENT(  OUT)     :: ErrStat                                      !< Error status
   CHARACTER(*),     INTENT(  OUT)     :: ErrMsg                                       !< Error message if ErrStat /= ErrID_None


      ! Local variables

   REAL(SiKi)                          :: Coords(3)                                    !< coordinates with wave directions converted to range [-180, 180)
   INTEGER(IntKi)                      :: i                                            !< generic counter
   INTEGER(IntKi)                      :: n(3)                                         !< number of points in WvFreq1, WvDir1, and WvDir2
   
   INTEGER(IntKi)                      :: Indx_Lo(3)                                   !< index associated with lower bound of dimension 1,2,3 where val(Indx_lo(1,2,3)) <= InCoord(1,2,3) <= val(Indx_hi(1,2,3))
   INTEGER(IntKi)                      :: Indx_Hi(3)                                   !< index associated with upper bound of dimension 1,2,3 where val(Indx_lo(1,2,3)) <= InCoord(1,2,3) <= val(Indx_hi(1,2,3))
   REAL(SiKi)                          :: Pos_Lo(3)                                    !< coordinate value with lower bound of dimension 1,2,3 
   REAL(SiKi)                          :: Pos_Hi(3)                                    !< coordinate value with upper bound of dimension 1,2,3 
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""   
      
   n(1) = SIZE(WvFreq1)
   n(2) = SIZE(WvDir1)
   n(3) = SIZE(WvDir2)

         ! Check sizes of input:
   IF ( SIZE( DataSet3D, 1 ) /= n(1) )  CALL SetErrStat( ErrID_Fatal, 'Invalid dimensions: size(DataSet3D,1) must be the size(WvFreq1)', ErrStat, ErrMsg, 'WAMIT_Interp3D_Cplx' )
   IF ( SIZE( DataSet3D, 2 ) /= n(2) )  CALL SetErrStat( ErrID_Fatal, 'Invalid dimensions: size(DataSet3D,2) must be the size(WvDir1)',  ErrStat, ErrMsg, 'WAMIT_Interp3D_Cplx' )
   IF ( SIZE( DataSet3D, 3 ) /= n(3) )  CALL SetErrStat( ErrID_Fatal, 'Invalid dimensions: size(DataSet3D,3) must be the size(WvDir2)',  ErrStat, ErrMsg, 'WAMIT_Interp3D_Cplx' )      
   IF (ErrStat == ErrID_Fatal) RETURN
   

      ! find the indices into the arrays representing coordinates of each dimension:

      Coords = InCoord
      
   DO i=2,3  ! make sure these requested degrees fall in the range -180 <= Coord(2:3) < 180
      Coords(i) = MODULO( Coords(i), 360.0_SiKi )
      IF ( Coords(i) >= 180.0_SiKi ) Coords(i) = Coords(i) - 360.0_SiKi
   END DO
   
   CALL LocateStp( Coords(1), WvFreq1, LastIndex(1), n(1) )
   CALL LocateStp( Coords(2), WvDir1,  LastIndex(2), n(2) )
   CALL LocateStp( Coords(3), WvDir2,  LastIndex(3), n(3) )
   
   Indx_Lo = LastIndex  ! at this point, 0 <= Indx_Lo(i) <= n(i) for all i
   
   
   ! WvFreq1 (indx 1)
   IF (Indx_Lo(1) == 0) THEN
      Indx_Lo(1) = 1
   ELSEIF (Indx_Lo(1) == n(1) ) THEN
      Indx_Lo(1) = max( n(1) - 1, 1 )                    ! make sure it's a valid index
   END IF     
   Indx_Hi(1) = min( Indx_Lo(1) + 1 , n(1) )             ! make sure it's a valid index

   
   ! WvDir1, WvDir2 (indx 2,3)   [use modular arithmetic]
   DO i=2,3
      IF (Indx_Lo(i) == 0) THEN
         Indx_Hi(i) = 1                           
         Indx_Lo(i) = n(i)
      ELSEIF (Indx_Lo(i) == n(i) ) THEN
         Indx_Hi(i) = 1      
      ELSE
         Indx_Hi(i) = min( Indx_Lo(i) + 1, n(i) )        ! make sure it's a valid index
      END IF      
   END DO
      
      ! calculate the positions of all dimensions:
      
   pos_Lo(1) = WvFreq1(Indx_Lo(1))
   pos_Hi(1) = WvFreq1(Indx_Hi(1))
      
   pos_Lo(2) = WvDir1(Indx_Lo(2))
   pos_Hi(2) = WvDir1(Indx_Hi(2))
   
   pos_Lo(3) = WvDir2(Indx_Lo(3))
   pos_Hi(3) = WvDir2(Indx_Hi(3))
   
      ! angles have to be adjusted so that pos_Lo(i) <= Coords(i) <= pos_Hi(i)
   DO i=2,3      
      IF ( Indx_Hi(i) == 1 .AND. n(i) > 1 )  THEN ! we're looping around the array [periodic]
         IF ( pos_Lo(i) < Coords(i) ) THEN
            pos_Hi(i) = pos_Hi(i) + 360.0_SiKi
         ELSEIF ( pos_Lo(i) /= Coords(i) ) THEN !bjj: I think it's okay if we don't use equalRealNos here
            pos_Lo(i) = pos_Lo(i) - 360.0_SiKi 
         END IF
      END IF   
   END DO
   


   CALL Interp3D_withIndx_Cplx( Coords, DataSet3D, Indx_Lo, Indx_Hi, pos_Lo, pos_Hi, OutForce )
      
   
END SUBROUTINE WAMIT_Interp3D_Cplx

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine takes the complex valued QTF dataset and interpolates for the desired coordinate in (Omega,Omega,WaveDir,
!!    WaveDir) space.
!!
!! A few important notes concerning this subroutine:
!!    1. It is complex valued.  The values represent the second order wave force as calculated by WAMIT.
!!    2. The dimenions of DataSet4D are Frequency1 (positive valued), Frequency2 (positive valued), Wave Direction 1 (degrees),
!!          and Wave Direction 2 (degrees).
!!    3. The wave direction requested might be between end points of wave direction dimension (ie. at 179 degrees when 
!!          WvDir1(1)=175, WvDir(Dims(3))=-175)
!!    4. The arrays WvFreq1, WvFreq2, WvDir1, and WvDir2 will give the values for each dimension that correspond to
!!          each index of DataSet4D.
!!    5. The data is not necessarily equally spaced in any direction: ie. WvFreq1 may not have uniform spacing between points.
!!    6. If a point is requested, it can be assumed that it lies within DataSet4D (this is checked before calling this subroutine)
!!    7. LastIndex contains the index numbers for the lowest bound on the indexes used on the last interpolation.  If they are 0,
!!          then assume this is the first call to this subroutine.
!!    8. DataSet4D is complete and not sparse.  There are no missing values.
!!
!!
SUBROUTINE WAMIT_Interp4D_Cplx( InCoord, DataSet4D, WvFreq1, WvFreq2, WvDir1, WvDir2, LastIndex, OutForce, ErrStat, ErrMsg )

      ! I/O variables

   REAL(SiKi),       INTENT(IN   )     :: InCoord(4)                                   !< Arranged as (Omega1, Omega2, WaveDir1)
   COMPLEX(SiKi),    INTENT(IN   )     :: DataSet4D(:,:,:,:)                           !< Arranged as Index 1= Omega1, Index 2= Omega2, Index 3= WaveDir1, Index 4= Wavedir2.
   REAL(SiKi),       INTENT(IN   )     :: WvFreq1(:)                                   !< Frequencies associated with Index 1 of DataSet4D
   REAL(SiKi),       INTENT(IN   )     :: WvFreq2(:)                                   !< Frequencies associated with Index 2 of DataSet4D
   REAL(SiKi),       INTENT(IN   )     :: WvDir1(:)                                    !< Frequencies associated with Index 3 of DataSet4D
   REAL(SiKi),       INTENT(IN   )     :: WvDir2(:)                                    !< Frequencies associated with Index 3 of DataSet4D
   INTEGER(IntKi),   INTENT(INOUT)     :: LastIndex(4)                                 !< Index for the last (Omega1, Omega2, WaveDir1) used
   COMPLEX(SiKi),    INTENT(  OUT)     :: OutForce                                     !< The interpolated resulting force from DataSet4D
   INTEGER(IntKi),   INTENT(  OUT)     :: ErrStat                                      !< Error status
   CHARACTER(*),     INTENT(  OUT)     :: ErrMsg                                       !< Error message if ErrStat /= ErrID_None


      ! Local variables

   REAL(SiKi)                          :: Coords(4)                                    !< coordinates with wave directions converted to range [-180, 180)
   INTEGER(IntKi)                      :: i                                            !< generic counter
   INTEGER(IntKi)                      :: n(4)                                         !< number of points in WvFreq1, WvFreq2, WvDir1, and WvDir2
   
   INTEGER(IntKi)                      :: Indx_Lo(4)                                   !< index associated with lower bound of dimension 1-4 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   INTEGER(IntKi)                      :: Indx_Hi(4)                                   !< index associated with upper bound of dimension 1-4 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   REAL(SiKi)                          :: Pos_Lo(4)                                    !< coordinate value with lower bound of dimension 1-4 
   REAL(SiKi)                          :: Pos_Hi(4)                                    !< coordinate value with upper bound of dimension 1-4 


   
   
   ErrStat = ErrID_None
   ErrMsg  = ""   
      
   n(1) = SIZE(WvFreq1)
   n(2) = SIZE(WvFreq2)
   n(3) = SIZE(WvDir1)
   n(4) = SIZE(WvDir2)

         ! Check sizes of input:
   IF ( SIZE( DataSet4D, 1 ) /= n(1) )  CALL SetErrStat( ErrID_Fatal, 'Invalid dimensions: size(DataSet4D,1) must be the size(WvFreq1)', ErrStat, ErrMsg, 'WAMIT_Interp4D_Cplx' )
   IF ( SIZE( DataSet4D, 2 ) /= n(2) )  CALL SetErrStat( ErrID_Fatal, 'Invalid dimensions: size(DataSet4D,2) must be the size(WvFreq2)', ErrStat, ErrMsg, 'WAMIT_Interp4D_Cplx' )
   IF ( SIZE( DataSet4D, 3 ) /= n(3) )  CALL SetErrStat( ErrID_Fatal, 'Invalid dimensions: size(DataSet4D,3) must be the size(WvDir1)',  ErrStat, ErrMsg, 'WAMIT_Interp4D_Cplx' )
   IF ( SIZE( DataSet4D, 4 ) /= n(4) )  CALL SetErrStat( ErrID_Fatal, 'Invalid dimensions: size(DataSet4D,4) must be the size(WvDir2)',  ErrStat, ErrMsg, 'WAMIT_Interp4D_Cplx' )      
   IF (ErrStat == ErrID_Fatal) RETURN
   

      ! find the indices into the arrays representing coordinates of each dimension:
      
   Coords = InCoord
      
   DO i=3,4  ! make sure these requested degrees fall in the range -180 <= Coord(3:4) < 180
      Coords(i) = MODULO( Coords(i), 360.0_SiKi )
      IF ( Coords(i) >= 180.0_SiKi ) Coords(i) = Coords(i) - 360.0_SiKi
   END DO
   
   CALL LocateStp( Coords(1), WvFreq1, LastIndex(1), n(1) )
   CALL LocateStp( Coords(2), WvFreq2, LastIndex(2), n(2) )
   CALL LocateStp( Coords(3), WvDir1,  LastIndex(3), n(3) )
   CALL LocateStp( Coords(4), WvDir2,  LastIndex(4), n(4) )
   
   Indx_Lo = LastIndex
   
   
   ! WvFreq1, WvFreq2 (indx 1, 2)
   DO i=1,2   
      IF (Indx_Lo(i) == 0) THEN
         Indx_Lo(i) = 1
      ELSEIF (Indx_Lo(i) == n(i) ) THEN
         Indx_Lo(i) = max( n(i) - 1, 1 )           ! make sure it's a valid index
      END IF      
      Indx_Hi(i) = min( Indx_Lo(i) + 1, n(i) )     ! make sure it's a valid index
   END DO
         
   
   ! WvDir1, WvDir2 (indx 3,4)   [use modular arithmetic]
   DO i=3,4
      IF (Indx_Lo(i) == 0) THEN
         Indx_Hi(i) = 1                           
         Indx_Lo(i) = n(i)
      ELSEIF (Indx_Lo(i) == n(i) ) THEN
         Indx_Hi(i) = 1      
      ELSE
         Indx_Hi(i) = min( Indx_Lo(i) + 1, n(i) )     ! make sure it's a valid index
      END IF      
   END DO
      
   
      ! calculate the positions of all dimensions:
      
   pos_Lo(1) = WvFreq1(Indx_Lo(1))
   pos_Hi(1) = WvFreq1(Indx_Hi(1))
      
   pos_Lo(2) = WvFreq2(Indx_Lo(2))
   pos_Hi(2) = WvFreq2(Indx_Hi(2))
   
   pos_Lo(3) = WvDir1(Indx_Lo(3))
   pos_Hi(3) = WvDir1(Indx_Hi(3))
   
   pos_Lo(4) = WvDir2(Indx_Lo(4))
   pos_Hi(4) = WvDir2(Indx_Hi(4))
   
      ! angles have to be adjusted so that pos_Lo(i) <= Coords(i) <= pos_Hi(i)
   DO i=3,4
      IF ( Indx_Hi(i) == 1 .AND. n(i) > 1 )  THEN ! we're looping around the array [periodic]
         IF ( pos_Lo(i) < Coords(i) ) THEN
            pos_Hi(i) = pos_Hi(i) + 360.0_SiKi
         ELSEIF ( pos_Lo(i) /= Coords(i) ) THEN !bjj: I think it's okay if we don't use equalRealNos here
            pos_Lo(i) = pos_Lo(i) - 360.0_SiKi 
         END IF
      END IF   
   END DO
   


   CALL Interp4D_withIndx_Cplx( Coords, DataSet4D, Indx_Lo, Indx_Hi, pos_Lo, pos_Hi, OutForce )   
   
   
   

END SUBROUTINE WAMIT_Interp4D_Cplx

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine takes a complex valued 2D dataset (InData2D) with the positions and indices of the bounding box and
!! performs 3-D linear interpolation to estimate the value of the dataset at InCoord.
!! It does not check that indices are valid or that the bounding box is non-degenerate (i.e., that posHi /= posLo) 
!!
!! This method is described here: http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
!!
SUBROUTINE Interp2D_withIndx_Cplx( InCoord, InData2D, Indx_Lo, Indx_Hi, posLo, posHi, InterpVal )

   REAL(SiKi),     INTENT(IN   )          :: InCoord(2)                             !< Arranged as (x, y)
   COMPLEX(SiKi),  INTENT(IN   )          :: InData2D(:,:)                          !< Arranged as Index 1= x, Index 2= y.
   COMPLEX(SiKi),  INTENT(  OUT)          :: InterpVal                              !< The interpolated value of DataSet2D at InCoord

   REAL(SiKi),     INTENT(IN   )          :: posLo(2)                               !< coordinate values associated with Indx_Lo 
   REAL(SiKi),     INTENT(IN   )          :: posHi(2)                               !< coordinate values associated with Indx_Hi

   INTEGER(IntKi), INTENT(IN   )          :: Indx_Lo(2)                             !< index associated with lower bound of dimension 1 where x(xIndx_lo) <= InCoord(1) <= x(xIndx_hi)
   INTEGER(IntKi), INTENT(IN   )          :: Indx_Hi(2)                             !< index associated with upper bound of dimension 1 where x(xIndx_lo) <= InCoord(1) <= x(xIndx_hi)
       
   ! local variables
   REAL(SiKi)                             :: isopc(2)                               ! isoparametric coordinates 
   
   REAL(SiKi)                             :: N(4)                                   ! size 2^n
   COMPLEX(SiKi)                          :: u(4)                                   ! size 2^n
   
   
   CALL CalcIsoparCoords( InCoord, posLo, posHi, isopc )      ! Calculate iospc
   

   N(1)  = ( 1.0_SiKi + isopc(1) )*( 1.0_SiKi - isopc(2) )
   N(2)  = ( 1.0_SiKi + isopc(1) )*( 1.0_SiKi + isopc(2) )
   N(3)  = ( 1.0_SiKi - isopc(1) )*( 1.0_SiKi + isopc(2) )
   N(4)  = ( 1.0_SiKi - isopc(1) )*( 1.0_SiKi - isopc(2) )
   N     = N / REAL( SIZE(N), SiKi )  ! normalize
      
   u(1)  = InData2D( Indx_Hi(1), Indx_Lo(2) )
   u(2)  = InData2D( Indx_Hi(1), Indx_Hi(2) )
   u(3)  = InData2D( Indx_Lo(1), Indx_Hi(2) )
   u(4)  = InData2D( Indx_Lo(1), Indx_Lo(2) )
   
   
   InterpVal = SUM ( N * u ) 

   
END SUBROUTINE Interp2D_withIndx_Cplx

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine takes a complex valued 3D dataset (InData3D) with the positions and indices of the bounding box and
!! performs 3-D linear interpolation to estimate the value of the dataset at InCoord.
!! It does not check that indices are valid or that the bounding box is non-degenerate (i.e., that posHi /= posLo) 
!!
!! This method is described here: http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/AFEM.Ch11.d/AFEM.Ch11.pdf
!!
SUBROUTINE Interp3D_withIndx_Cplx( InCoord, InData3D, Indx_Lo, Indx_Hi, posLo, posHi, InterpVal )

   REAL(SiKi),     INTENT(IN   )          :: InCoord(3)                             !< Arranged as (x, y, z)
   COMPLEX(SiKi),  INTENT(IN   )          :: InData3D(:,:,:)                        !< Arranged as Index 1= x, Index 2= y, Index 3= z.
   COMPLEX(SiKi),  INTENT(  OUT)          :: InterpVal                              !< The interpolated value of DataSet3D at InCoord

   REAL(SiKi),     INTENT(IN   )          :: posLo(3)                               !< xyz (coordinate) values associated with Indx_Lo 
   REAL(SiKi),     INTENT(IN   )          :: posHi(3)                               !< xyz (coordinate) values associated with Indx_Hi

   INTEGER(IntKi), INTENT(IN   )          :: Indx_Lo(3)                             !< index associated with lower bound of dimension 1 where x(xIndx_lo) <= InCoord(1) <= x(xIndx_hi)
   INTEGER(IntKi), INTENT(IN   )          :: Indx_Hi(3)                             !< index associated with upper bound of dimension 1 where x(xIndx_lo) <= InCoord(1) <= x(xIndx_hi)
       
   ! local variables
   REAL(SiKi)                             :: isopc(3)                               ! isoparametric hexahedral coordinates (natural coordinates) [ xi, eta, mu ]
   
   REAL(SiKi)                             :: N(8)                                   ! size 2^n
   COMPLEX(SiKi)                          :: u(8)                                   ! size 2^n

   
   CALL CalcIsoparCoords( InCoord, posLo, posHi, isopc )      ! Calculate iospc
      
   !eta  = ( 2.0_SiKi*InCoord(1) - posLo(1) - posHi(1) ) / ( posHi(1) - posLo(1) ) !< (2*x - x1 - x2 ) / (x2-x1) !note that this is actually negative eta from the referenced paper, but it follows the other 2 dimensions better if we use the negative here and then flip the signs in N(:) (I think there is a bug in that paper)
   !xi   = ( 2.0_SiKi*InCoord(2) - posLo(2) - posHi(2) ) / ( posHi(2) - posLo(2) ) !< (2*y - y1 - y2 ) / (y2-y1)
   !mu   = ( 2.0_SiKi*InCoord(3) - posLo(3) - posHi(3) ) / ( posHi(3) - posLo(3) ) !< (2*z - z1 - z2 ) / (z2-z1)

   N(1)  = ( 1.0_SiKi + isopc(1) )*( 1.0_SiKi - isopc(2) )*( 1.0_SiKi - isopc(3) )
   N(2)  = ( 1.0_SiKi + isopc(1) )*( 1.0_SiKi + isopc(2) )*( 1.0_SiKi - isopc(3) )
   N(3)  = ( 1.0_SiKi - isopc(1) )*( 1.0_SiKi + isopc(2) )*( 1.0_SiKi - isopc(3) )
   N(4)  = ( 1.0_SiKi - isopc(1) )*( 1.0_SiKi - isopc(2) )*( 1.0_SiKi - isopc(3) )
   N(5)  = ( 1.0_SiKi + isopc(1) )*( 1.0_SiKi - isopc(2) )*( 1.0_SiKi + isopc(3) )
   N(6)  = ( 1.0_SiKi + isopc(1) )*( 1.0_SiKi + isopc(2) )*( 1.0_SiKi + isopc(3) )
   N(7)  = ( 1.0_SiKi - isopc(1) )*( 1.0_SiKi + isopc(2) )*( 1.0_SiKi + isopc(3) )
   N(8)  = ( 1.0_SiKi - isopc(1) )*( 1.0_SiKi - isopc(2) )*( 1.0_SiKi + isopc(3) )
   N     = N / REAL( SIZE(N), SiKi )  ! normalize
      
   u(1)  = InData3D( Indx_Hi(1), Indx_Lo(2), Indx_Lo(3) )
   u(2)  = InData3D( Indx_Hi(1), Indx_Hi(2), Indx_Lo(3) )
   u(3)  = InData3D( Indx_Lo(1), Indx_Hi(2), Indx_Lo(3) )
   u(4)  = InData3D( Indx_Lo(1), Indx_Lo(2), Indx_Lo(3) )
   u(5)  = InData3D( Indx_Hi(1), Indx_Lo(2), Indx_Hi(3) )
   u(6)  = InData3D( Indx_Hi(1), Indx_Hi(2), Indx_Hi(3) )
   u(7)  = InData3D( Indx_Lo(1), Indx_Hi(2), Indx_Hi(3) )
   u(8)  = InData3D( Indx_Lo(1), Indx_Lo(2), Indx_Hi(3) )
   
   
   InterpVal = SUM ( N * u ) 

   
END SUBROUTINE Interp3D_withIndx_Cplx

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine takes a complex valued 4D dataset (InData4D) with the positions and indices of the bounding box and
!! performs 4-D linear interpolation to estimate the value of the dataset at InCoord.
!! It does not check that indices are valid or that the bounding box is non-degenerate (i.e., that posHi /= posLo) 
!!
!! This method is described here: http://rjwagner49.com/Mathematics/Interpolation.pdf
!!
SUBROUTINE Interp4D_withIndx_Cplx( InCoord, InData4D, Indx_Lo, Indx_Hi, posLo, posHi, InterpVal )

   REAL(SiKi),     INTENT(IN   )          :: InCoord(4)                             !< Arranged as (x1, x2, x3, x4 )
   COMPLEX(SiKi),  INTENT(IN   )          :: InData4D(:,:,:,:)                        !< Arranged as Index 1= x1, Index 2= x2, Index 3= x3, Index 4= x4.
   COMPLEX(SiKi),  INTENT(  OUT)          :: InterpVal                              !< The interpolated value of InData4D at InCoord

   REAL(SiKi),     INTENT(IN   )          :: posLo(4)                               !< coordinate values associated with Indx_Lo 
   REAL(SiKi),     INTENT(IN   )          :: posHi(4)                               !< coordinate values associated with Indx_Hi

   INTEGER(IntKi), INTENT(IN   )          :: Indx_Lo(4)                             !< index associated wtih lower bound of dimension 1 where x(xIndx_lo) <= InCoord(1) <= x(xIndx_hi)
   INTEGER(IntKi), INTENT(IN   )          :: Indx_Hi(4)                             !< index associated wtih upper bound of dimension 1 where x(xIndx_lo) <= InCoord(1) <= x(xIndx_hi)
       
   ! local variables
   REAL(SiKi)                             :: isopc(4)                               ! isoparametric coordinates 
   REAL(SiKi)                             :: N(16)                                  ! size 2^n
   COMPLEX(SiKi)                          :: u(16)                                  ! size 2^n
   
   
   CALL CalcIsoparCoords( InCoord, posLo, posHi, isopc )      ! Calculate iospc
         
   
   N( 1) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N( 2) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N( 3) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N( 4) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )    
   N( 5) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N( 6) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N( 7) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N( 8) = ( 1.0_SiKi - isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )    
   N( 9) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N(10) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N(11) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N(12) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi - isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N(13) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N(14) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi - isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N(15) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi - isopc(4) )
   N(16) = ( 1.0_SiKi + isopc(1) ) * ( 1.0_SiKi + isopc(2) ) * ( 1.0_SiKi + isopc(3) ) * ( 1.0_SiKi + isopc(4) )
   N     = N / REAL( SIZE(N), SiKi )  ! normalize
   
   
   u( 1) = InData4D( Indx_Lo(1), Indx_Lo(2), Indx_Lo(3), Indx_Lo(4) )
   u( 2) = InData4D( Indx_Lo(1), Indx_Lo(2), Indx_Lo(3), Indx_Hi(4) )
   u( 3) = InData4D( Indx_Lo(1), Indx_Lo(2), Indx_Hi(3), Indx_Lo(4) )
   u( 4) = InData4D( Indx_Lo(1), Indx_Lo(2), Indx_Hi(3), Indx_Hi(4) )
   u( 5) = InData4D( Indx_Lo(1), Indx_Hi(2), Indx_Lo(3), Indx_Lo(4) )
   u( 6) = InData4D( Indx_Lo(1), Indx_Hi(2), Indx_Lo(3), Indx_Hi(4) )
   u( 7) = InData4D( Indx_Lo(1), Indx_Hi(2), Indx_Hi(3), Indx_Lo(4) )
   u( 8) = InData4D( Indx_Lo(1), Indx_Hi(2), Indx_Hi(3), Indx_Hi(4) )
   u( 9) = InData4D( Indx_Hi(1), Indx_Lo(2), Indx_Lo(3), Indx_Lo(4) )
   u(10) = InData4D( Indx_Hi(1), Indx_Lo(2), Indx_Lo(3), Indx_Hi(4) )
   u(11) = InData4D( Indx_Hi(1), Indx_Lo(2), Indx_Hi(3), Indx_Lo(4) )
   u(12) = InData4D( Indx_Hi(1), Indx_Lo(2), Indx_Hi(3), Indx_Hi(4) )
   u(13) = InData4D( Indx_Hi(1), Indx_Hi(2), Indx_Lo(3), Indx_Lo(4) )
   u(14) = InData4D( Indx_Hi(1), Indx_Hi(2), Indx_Lo(3), Indx_Hi(4) )
   u(15) = InData4D( Indx_Hi(1), Indx_Hi(2), Indx_Hi(3), Indx_Lo(4) )
   u(16) = InData4D( Indx_Hi(1), Indx_Hi(2), Indx_Hi(3), Indx_Hi(4) )
      
   
   InterpVal = SUM ( N * u ) 

   
END SUBROUTINE Interp4D_withIndx_Cplx


!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine calculates the iosparametric coordinates, isopc, which is a value between -1 and 1 (for each dimension of a dataset)
!! indicating where InCoord falls between posLo and posHi.
!!
SUBROUTINE CalcIsoparCoords( InCoord, posLo, posHi, isopc )


   REAL(SiKi),     INTENT(IN   )          :: InCoord(:)                             !< Arranged as (x, y)
   REAL(SiKi),     INTENT(IN   )          :: posLo(:)                               !< coordinate values associated with Indx_Lo 
   REAL(SiKi),     INTENT(IN   )          :: posHi(:)                               !< coordinate values associated with Indx_Hi
   REAL(SiKi),     INTENT(  OUT)          :: isopc(:)                               ! isoparametric coordinates 

   ! local variables
   REAL(SiKi)                             :: dx                                     ! difference between high and low coordinates in the bounding "box"
   INTEGER(IntKi)                         :: i                                      !  loop counter
   
   
   do i=1,size(isopc)
      
      dx = posHi(i) - posLo(i) 
      if (EqualRealNos(dx, 0.0_SiKi)) then
         isopc(i) = 1.0_SiKi
      else
         isopc(i) = ( 2.0_SiKi*InCoord(i) - posLo(i) - posHi(i) ) / dx
            ! to verify that we don't extrapolate, make sure this is bound between -1 and 1 (effectively nearest neighbor)
         isopc(i) = min( 1.0_SiKi, isopc(i) )
         isopc(i) = max(-1.0_SiKi, isopc(i) )
      end if
      
   end do
            
END SUBROUTINE CalcIsoparCoords

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE WAMIT_Interp
