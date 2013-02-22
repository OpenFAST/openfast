!**********************************************************************************************************************************
! (c) 2012 National Renewable Energy Laboratory
!
! This file is used to compile the template module for the FAST modularization framework. We are currently making changes to the 
! definitions of the pack/unpack routines, so expect this file to change or be replaced/removed. 
!
! This file contains two subroutines to (1) place data from an array into a (presumably larger) array of bytes (buffer) and 
! (2) retrieve data in an array from an array of bytes (buffer). The subroutines are not placed in a module (or defined by an 
! interface) so that calling routines can declare the InData and OutData arrays of arbitrary types (reals, doubles, characters, 
! integers, etc) and their data is then placed into an array of bytes. These routines are designed for save/restart capability 
! in the FAST modular framework. These routines are a way to avoid using EQUIVALENCE statements in the calling routines.
! To retrieve the same data that was placed in the buffer array, GetFromBuffer must be called in the same sequence as PlaceInBuffer
! was called (using arrays of the same type and size).
!**********************************************************************************************************************************
SUBROUTINE PlaceInBuffer( InData, BufferData, Sz, Indx )
! This routine takes an array InData, which contains Sz bytes of data, and places it into the BufferData array starting at
! the Indx+1 byte location. Indx is updated in this routine with each byte copied, and thus represents the number of bytes placed
! in the BufferData array. The calling routine should initialize Indx to 0 before any calls to PlaceInBuffer; we do not recommend
! changing its value outside of the PlaceInBuffer routine.
!----------------------------------------------------------------------------------------------------------------------------------

   USE Precision
   IMPLICIT NONE


   INTEGER(IntKi), INTENT(IN)     :: Sz                        ! Size--in bytes--of the data to be transferred to the BufferData array
   INTEGER(IntKi), INTENT(INOUT)  :: Indx                      ! Index into the buffer/bytes transferred
   INTEGER(B1Ki),  INTENT(IN)     :: InData(Sz)                ! Data to put into a buffer, treated as 1-byte integers (though it could be another type in the calling routine)
   INTEGER(B1Ki),  INTENT(INOUT)  :: BufferData(*)             ! The buffer that stores the data, treated as 1-byte integers    
   
   INTEGER(IntKi)                 :: I                         ! Loop counter
   
   
      ! Transfer the data, byte-by-byte
         
   DO I = 1,Sz
      Indx = Indx + 1                     
      BufferData(Indx) = InData(I)      
   END DO

   RETURN

END SUBROUTINE PlaceInBuffer
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE GetFromBuffer( BufferData, OutData, Sz, Indx )
! This routine takes the BufferData array and transfers Sz bytes into the OutData array, starting at the Indx + 1 byte location.
! Indx is updated in this routine with each byte copied, and thus represents the number of bytes retreived from the BufferData
! array. The calling routine should initialize Indx to 0 before any calls to GetFromBuffer; we do not recommend
! changing its value outside of the GetFromBuffer routine.
!----------------------------------------------------------------------------------------------------------------------------------

   USE Precision
   IMPLICIT NONE

   INTEGER(B1Ki),  INTENT(IN)     :: BufferData(*)             ! The buffer that stores the data, treated as 1-byte integers (though it could be another type)      
   INTEGER(IntKi), INTENT(IN)     :: Sz                        ! Size--in bytes--of the data to be transferred from the BufferData array
   INTEGER(B1Ki),  INTENT(  OUT)  :: OutData(Sz)               ! Data returned from the buffer, treated as 1-byte integers (though it could be any other type)
   INTEGER(IntKi), INTENT(INOUT)  :: Indx                      ! Current index into the buffer (it is recommended that the calling program not modify this value)
   
   INTEGER(IntKi)                 :: I                         ! Loop counter
   
   
      ! Transfer the data, byte-by-byte
   
   DO I = 1,Sz
      Indx = Indx + 1
      OutData(I) = BufferData(Indx)
   END DO

   RETURN

END SUBROUTINE GetFromBuffer
!----------------------------------------------------------------------------------------------------------------------------------
