!=======================================================================
SUBROUTINE Mod_vKrm ( Ht, Ucmp, Spec )

   ! This subroutine defines the "Improved" von Karman PSD model.
   ! The use of this subroutine requires that all variables have the units of meters and seconds.
   USE                        NWTC_Library

   IMPLICIT                   NONE
      
      !Passed variables

   REAL(ReKi), INTENT(IN)  :: Ht                   ! height
   REAL(ReKi), INTENT(IN)  :: UCmp                 ! wind speed
   REAL(ReKi)              :: Spec   (:,:)
 

RETURN
END SUBROUTINE Mod_vKrm
!=======================================================================
SUBROUTINE ScaleMODVKM(Ht,UCmp, LambdaU, LambdaV, LambdaW)

! THIS SUBROUTINE DEFINES HUB SCALE PARMS FOR Modified von KARMAN PSD MODEL 
   USE                        NWTC_Library

   IMPLICIT NONE

   REAL(ReKi), INTENT(IN)  :: Ht     ! height
   REAL(ReKi), INTENT(IN)  :: UCmp   ! wind speed

   REAL(ReKi)              :: LambdaU
   REAL(ReKi)              :: LambdaV
   REAL(ReKi)              :: LambdaW
 
RETURN
END SUBROUTINE ScaleMODVKM
!=======================================================================
FUNCTION FindZ0(z, sigma, U, f)

    ! This function is used in the Modified von Karman model to 
    ! determine the necessary surface roughness length for a given sigma.
   USE                     NWTC_Library

   IMPLICIT                NONE

   REAL(ReKi)            :: FindZ0        ! Estimated surface roughness length
   REAL(ReKi),INTENT(IN) :: z             ! Hub height
   REAL(ReKi),INTENT(IN) :: sigma         ! Target std deviation
   REAL(ReKi),INTENT(IN) :: U             ! Hub height wind speed
   REAL(ReKi),INTENT(IN) :: f             ! Coriolis parameter
   
   FindZ0 = 1.0   ! a default value

RETURN
END FUNCTION FindZ0
!=======================================================================
FUNCTION CalcDiff(z0Guess, z, sigma, U, f)

      ! This function calculates the difference between the specified
      ! sigma and the calculated one.
   USE                     NWTC_Library

   IMPLICIT                NONE
   
   REAL(ReKi)             :: CalcDiff  ! Output - will be nearly zero if surface roughness is correct
   REAL(ReKi), INTENT(IN) :: z0Guess   ! estimated surface roughness
   REAL(ReKi), INTENT(IN) :: z         ! Hub height (m)
   REAL(ReKi), INTENT(IN) :: sigma     ! Target standard deviation (m/s)
   REAL(ReKi), INTENT(IN) :: U         ! Mean hub-height wind speed (m/s)
   REAL(ReKi), INTENT(IN) :: f         ! Coriolis parameter 

   CalcDiff = 0.0
   
RETURN
END FUNCTION CalcDiff
!=======================================================================