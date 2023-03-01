!**********************************************************************************************************************************
! Copyright (C) 2013-2014  National Renewable Energy Laboratory
!
! This code provides a wrapper for the ScaLAPACK routines currently used at the NWTC (mainly codes in the FAST framework).
!
!**********************************************************************************************************************************
MODULE NWTC_ScaLAPACK

   USE NWTC_Base        ! we only need the precision and error level constants
!   USE, INTRINSIC               :: ISO_C_Binding, only: C_FLOAT, C_DOUBLE          ! this is included in NWTC_Library

      ! Notes:

         ! Your project must include the following files:
         !     NWTC_Base.f90         [from NWTC Library]
         !     NWTC_ScaLAPACK.f90
         !     dlasrt2.f
         !     slasrt2.f

   IMPLICIT  NONE

   INTERFACE ScaLAPACK_LASRT  ! sort
      MODULE PROCEDURE ScaLAPACK_DLASRT2
      MODULE PROCEDURE ScaLAPACK_SLASRT2
   END INTERFACE   
   
   

CONTAINS

!=======================================================================
!INTERFACE ScaLAPACK_LASRT:
!  Sort the numbers in D in increasing order (if ID = 'I') or in decreasing order (if ID = 'D' ).
!  See documentation in  DLASRT2/SLASRT2 source code.
   SUBROUTINE ScaLAPACK_DLASRT2( ID, N, D, KEY, ErrStat, ErrMsg )
          
      ! passed parameters
 
      CHARACTER,       intent(in   ) :: ID                ! = 'I': sort D in increasing order; = 'D': sort D in decreasing order.
      INTEGER,         intent(in   ) :: N                 ! The length of the array D.
      INTEGER(IntKi),  intent(  out) :: ErrStat           ! Error level 
      CHARACTER(*),    intent(  out) :: ErrMsg            ! Message describing error

      !     .. Array Arguments ..
      INTEGER,         intent(inout) :: KEY( : )
      REAL(R8Ki)      ,intent(inout) :: D( : )            ! On entry, the array to be sorted. On exit, D has been sorted into increasing/decreasing order, depending on ID

         ! Local variable
      INTEGER                        :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value 
      
      
      ErrStat = ErrID_None
      ErrMsg  = ""
      
      CALL DLASRT2( ID, N, D, KEY, INFO )
                
      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "ScaLAPACK_DLSRT2: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'ScaLAPACK_DLSRT2: Unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF      
 
      
   RETURN
   END SUBROUTINE ScaLAPACK_DLASRT2
!=======================================================================
   SUBROUTINE ScaLAPACK_SLASRT2( ID, N, D, KEY, ErrStat, ErrMsg )
          
      ! passed parameters
 
      CHARACTER,      intent(in   ) :: ID                ! = 'I': sort D in increasing order; = 'D': sort D in decreasing order.
      INTEGER,        intent(in   ) :: N                 ! The length of the array D.
      INTEGER(IntKi), intent(  out) :: ErrStat           ! Error level 
      CHARACTER(*),   intent(  out) :: ErrMsg            ! Message describing error

      !     .. Array Arguments ..
      INTEGER,        intent(inout) :: KEY( : )
      REAL(R4Ki),     intent(inout) :: D( : )            ! On entry, the array to be sorted. On exit, D has been sorted into increasing/decreasing order, depending on ID

         ! Local variable
      INTEGER                       :: INFO              ! = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value 
      
      
      ErrStat = ErrID_None
      ErrMsg  = ""
      
      CALL SLASRT2( ID, N, D, KEY, INFO )
                
      IF (INFO /= 0) THEN
         ErrStat = ErrID_FATAL
         WRITE( ErrMsg, * ) INFO
         IF (INFO < 0) THEN
            ErrMsg  = "ScaLAPACK_SLSRT2: illegal value in argument "//TRIM(ErrMsg)//"."
         ELSE
            ErrMsg = 'ScaLAPACK_SLSRT2: Unknown error '//TRIM(ErrMsg)//'.'
         END IF
      END IF      
 
      
   RETURN
   END SUBROUTINE ScaLAPACK_SLASRT2
!=======================================================================
END MODULE NWTC_ScaLAPACK
