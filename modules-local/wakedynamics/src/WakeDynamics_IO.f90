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
   


END MODULE WakeDynamics_IO
