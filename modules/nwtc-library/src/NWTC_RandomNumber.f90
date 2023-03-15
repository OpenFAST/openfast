!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020 National Renewable Energy Laboratory
!
!    This file is part of NWTC Library.
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

MODULE NWTC_RandomNumber

USE Ran_Lux_Mod
USE NWTC_Library_Types
USE NWTC_IO

IMPLICIT NONE

INTEGER(IntKi), PARAMETER :: pRNG_RANLUX    = 1
INTEGER(IntKi), PARAMETER :: pRNG_INTRINSIC = 2

INTEGER, PARAMETER        :: LuxLevel       = 3       ! Luxury Level for RanLux RNG

!> \copydoc nwtc_randomnumber::uniformrandomnumbersr4
INTERFACE UniformRandomNumbers
MODULE PROCEDURE UniformRandomNumbersR4   ! 4-byte reals
MODULE PROCEDURE UniformRandomNumbersR8   ! 8-byte reals
END INTERFACE

CONTAINS
   
SUBROUTINE RandNum_Init(p, ErrStat, ErrMsg )
   
   ! Initialize the Random Number Generators
   
   IMPLICIT NONE
   
   TYPE(NWTC_RandomNumber_ParameterType),  INTENT(IN   ) :: p           ! PARAMETERs for random number generation
   INTEGER(IntKi)  ,             INTENT(OUT)   :: ErrStat     ! allocation status
   CHARACTER(*) ,                INTENT(OUT)   :: ErrMsg      ! error message

   INTEGER                          :: I           ! loop counter
   INTEGER(IntKi), ALLOCATABLE      :: NextSeed(:) ! The array that holds the next random seed for each component
   INTEGER                          :: NumSeeds    ! number of seeds in the intrinsic random number generator

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF (p%pRNG == pRNG_RANLUX) THEN
   
      CALL RLuxGo ( LuxLevel, ABS( p%RandSeed(1) ), 0, 0 )
      
      IF (.NOT. ALLOCATED( NextSeed ) ) THEN
         CALL AllocAry( NextSeed, 2, 'nextSeed', ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF
   
   ELSEIF (p%pRNG == pRNG_INTRINSIC) THEN
   
         ! determine the number of seeds necessary (gfortran needs 8 or 12 seeds, not just 2)
   
      CALL RANDOM_SEED ( SIZE = NumSeeds )
   
      IF ( NumSeeds /= 2 ) THEN
         CALL ProgWarn( ' The random number generator in use differs from the original code provided by NREL. This pRNG uses ' &
                           //trim(Int2LStr(NumSeeds))//' seeds instead of the 2 in the input file.')
      END IF
   
      IF ( .NOT. ALLOCATED( NextSeed ) ) THEN
         CALL AllocAry( NextSeed, NumSeeds, 'nextSeed', ErrSTat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN
      END IF
   
         ! We'll just populate this with odd seeds = Seed(1) and even seeds = Seed(2)
      DO I = 1, NumSeeds,2
         NextSeed(I) = p%RandSeed(1)
      END DO
      DO I = 2, NumSeeds,2
         NextSeed(I) = p%RandSeed(2)
      END DO
   
      CALL RANDOM_SEED( PUT=NextSeed )
   
   ELSE

      ! Invalid pRNG requested
      ErrMsg = "Invalid pRNG requested."
      ErrStat = ErrID_Fatal
      RETURN
   
   END IF

END SUBROUTINE RandNum_Init

!=======================================================================
!> This subroutine produces uniformly distributed random numbers, based on
!!   the pRNG requested. This routine assumes that the random number
!!   generator has been initialized earlier in the main program.
SUBROUTINE UniformRandomNumbersR4( pRNG_Type, RandomNumbers )

   IMPLICIT NONE

   INTEGER,    INTENT(IN   )               :: pRNG_Type
   REAL(SiKi), INTENT(  OUT), DIMENSION(:) :: RandomNumbers

   REAL(ReKi), ALLOCATABLE :: RN(:)

   IF ( pRNG_Type == pRNG_INTRINSIC ) THEN

      ! The Fortran intrinsic has an interface for various floating
      ! point types, so pass the variable directly
      CALL RANDOM_NUMBER( RandomNumbers )

   ELSEIF ( pRNG_Type == pRNG_RANLUX ) THEN
      
      ! RanLux, as implemented, uses ReKi, so cast the return value as needed
      ALLOCATE( RN( SIZE(RandomNumbers) ) )
      CALL RanLux ( RN )
      RandomNumbers = REAL(RN, KIND=SiKi)

   END IF

END SUBROUTINE UniformRandomNumbersR4
!=======================================================================
!> \copydoc nwtc_randomnumber::uniformrandomnumbersr4
SUBROUTINE UniformRandomNumbersR8( pRNG_Type, RandomNumbers )

   IMPLICIT NONE

   INTEGER,    INTENT(IN   )               :: pRNG_Type
   REAL(R8Ki), INTENT(  OUT), DIMENSION(:) :: RandomNumbers

   REAL(ReKi), ALLOCATABLE :: RN(:)

   IF ( pRNG_Type == pRNG_INTRINSIC ) THEN

      ! The Fortran intrinsic has an interface for various floating
      ! point types, so pass the variable directly
      CALL RANDOM_NUMBER( RandomNumbers )

   ELSEIF ( pRNG_Type == pRNG_RANLUX ) THEN

      ! RanLux, as implemented, uses ReKi, so cast the return value as needed
      ALLOCATE( RN( SIZE(RandomNumbers) ) )
      CALL RanLux ( RN )
      RandomNumbers = REAL(RN, KIND=R8Ki)

   END IF

END SUBROUTINE UniformRandomNumbersR8

END MODULE
