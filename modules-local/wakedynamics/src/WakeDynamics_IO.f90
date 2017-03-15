!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE WakeDynamics_IO
 
   use NWTC_Library
   use WakeDynamics_Types

   
   implicit none
   
   type(ProgDesc), parameter  :: WD_Ver = ProgDesc( 'WakeDynamics', 'v00.01.00', '21-Sep-2016' )
   character(*),   parameter  :: WD_Nickname = 'WD'
      
contains
   
   


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WD_PrintSum(  p, u, y, ErrStat, ErrMsg )
! This routine generates the summary file, which contains a summary of input file options.

      ! passed variables
   !TYPE(WD_InitInput),        INTENT(IN)  :: InputFileData                        ! Input-file data
   TYPE(WD_ParameterType),    INTENT(IN)  :: p                                    ! Parameters
   TYPE(WD_InputType),        INTENT(IN)  :: u                                    ! inputs 
   TYPE(WD_OutputType),       INTENT(IN)  :: y                                    ! outputs
   INTEGER(IntKi),            INTENT(OUT) :: ErrStat
   CHARACTER(*),              INTENT(OUT) :: ErrMsg


      ! Local variables.

   INTEGER(IntKi)               :: I                                               ! Index for the nodes.
   INTEGER(IntKi)               :: UnSu                                            ! I/O unit number for the summary output file

   CHARACTER(*), PARAMETER      :: FmtDat    = '(A,T35,1(:,F13.3))'                ! Format for outputting mass and modal data.
   CHARACTER(*), PARAMETER      :: FmtDatT   = '(A,T35,1(:,F13.8))'                ! Format for outputting time steps.

   CHARACTER(30)                :: OutPFmt                                         ! Format to print list of selected output channels to summary file
   CHARACTER(100)               :: Msg                                             ! temporary string for writing appropriate text to summary file

   ! Open the summary file and give it a heading.
      
   CALL GetNewUnit( UnSu, ErrStat, ErrMsg )
   CALL OpenFOutFile ( UnSu, TRIM( p%OutFileRoot )//'.sum', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN

 

   CLOSE(UnSu)

RETURN
END SUBROUTINE WD_PrintSum
!----------------------------------------------------------------------------------------------------------------------------------



END MODULE WakeDynamics_IO
