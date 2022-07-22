!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014  DNV KEMA Renewables, Inc.
!
!    This file is part of the IceFloe suite of subroutines
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
!************************************************************************

!******************************************************************************************
! routines for logging progress and errors in the ice
! loading code

module iceLog

   use NWTC_Base
   use NWTC_IO, only : GetNewUnit, OpenFOutFile, conv2uc, progwarn, progabort, NewLine

   implicit none

   public

   integer(IntKi), parameter  :: msgLen  = 1024

!  Data structure for message logging and error handling
   type iceFloe_LoggingType
   !  These are the status codes and error message that will be sent back to the glue code.
   !  They are not static and should be reset for every call to IceFloe_init, _update, _calcOutputs, etc.
      integer(IntKi)    :: UnitNum
      INTEGER(IntKi)    :: errID
      CHARACTER(msgLen) :: errMsg
      logical           :: warnFlag  !  Set if there have been warnings so user can fix even if other fatal error occurs
   end type

contains
!===============================================================================================
!  Error handling routines
!-----------------------------------------------------------------------------------------------
!  This routine updates the error ID, appends the message, and writes to screen and log file if requested
   subroutine iceErrorHndlr (iceLog, Err, Msg, Scrn)
      type (iceFloe_loggingType), intent(inout) :: iceLog
      integer(IntKi), intent(in)                :: Err    ! Error ID
      integer(IntKi), intent(in)                :: Scrn   ! write to screen?
      character(*), intent(in)                  :: Msg

      if ( Err /= ErrID_None ) then
!     If we have an existing message, append a new line then the new message
         if ( len_trim(iceLog%errMsg) > 0 ) iceLog%errMsg = trim(iceLog%errMsg)//NewLine
         iceLog%errMsg = trim(iceLog%errMsg)//trim(Msg)
         iceLog%errID = max(iceLog%ErrID, Err)

         if (Err == ErrID_Warn) iceLog%warnFlag = .true.
      endif

!   Post the message to the IceFloe log file (this will write to screen as well if requested)
      select case (Err)
         case (ErrID_Warn)
            call logWarn(iceLog, Msg, Scrn)
         case (ErrID_Fatal)
            call logFatal(iceLog, Msg, Scrn)
      end select

   end subroutine iceErrorHndlr

!-----------------------------------------------------------------------------------------------
!  This routine simply appends a message with a new line
   subroutine addMessage (iceLog, Msg)
      type (iceFloe_loggingType), intent(inout) :: iceLog
      character(*), intent(in)                  :: Msg
      iceLog%errMsg = trim(iceLog%errMsg)//newLine//trim(Msg)      
   end subroutine addMessage

!===============================================================================================
!  Message logging routines
!-----------------------------------------------------------------------------------------------
   subroutine openIceLog (iceLog, logFile)
      type (iceFloe_loggingType), intent(inout) :: iceLog
      character(*), optional, intent(in)        :: logFile
      INTEGER(IntKi)    :: Err        
      CHARACTER(msgLen) :: Msg         

      call GetNewUnit ( iceLog%UnitNum, Err, Msg )
      if (Err >= ErrID_Severe) then
         iceLog%ErrMsg = 'OpenIceLog: '//newLine//trim(Msg)
         return
      endif
      if (present(logFile)) then
         call OpenFOutFile ( iceLog%unitNum, logFile, Err, Msg )
      else
         call OpenFOutFile ( iceLog%unitNum, 'iceLog.txt', Err, Msg )
      endif
      if (Err >= ErrID_Severe) then
         iceLog%ErrMsg = 'OpenIceLog: '//newLine//trim(Msg)
         return
      endif
   end subroutine openIceLog

!--------------------------------------------------------
   subroutine closeIceLog (iceLog)
      type (iceFloe_loggingType), intent(in) :: iceLog
      close( iceLog%unitNum )
   end subroutine closeIceLog

!--------------------------------------------------------
   subroutine logMessage (iceLog, msg)
      type (iceFloe_loggingType), intent(in) :: iceLog
      character(*), intent(in)               :: msg
      write(iceLog%unitNum,'(A)') trim(msg)
   end subroutine logMessage

!--------------------------------------------------------
   subroutine logWarn (iceLog, msg, Scrn)
      type (iceFloe_loggingType), intent(in) :: iceLog
      integer(IntKi), intent(in)             :: scrn
      character(*), intent(in)               :: msg

!  Write to log file
      write(iceLog%unitNum,*) 
      write(iceLog%unitNum,'(A)') '************ WARNING ************************************'
      write(iceLog%unitNum,'(A)') trim(msg)
      write(iceLog%unitNum,'(A)') '*********************************************************'
      write(iceLog%unitNum,*) 
!  Write to screen (may be redundant w/ calling program)
      if (scrn==1) call ProgWarn(msg)

   end subroutine logWarn

!--------------------------------------------------------
   subroutine logFatal (iceLog, msg, scrn)
      type (iceFloe_loggingType), intent(in)  :: iceLog
      integer(IntKi), intent(in)             :: scrn
      character(*), intent(in)               :: msg

!  Write to log file
      write(iceLog%unitNum,*) 
      write(iceLog%unitNum,'(A)') '******** FATAL ERROR, PROGRAM WILL TERMINATE for the following reason ********'
      write(iceLog%unitNum,'(A)') trim(msg)
      write(iceLog%unitNum,'(A)') '******************************************************************************'
      write(iceLog%unitNum,*) 
!  Write to screen (may be redundant w/ calling program)
      if (scrn==1) call ProgAbort('Fatal error, calling program will abort: '//msg, .true.)

   end subroutine logFatal

end module iceLog