!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
!
!    This file is part of the NWTC Subroutine Library.
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
MODULE VersionInfo

   implicit none

contains

FUNCTION QueryGitVersion()

   CHARACTER(200) :: QueryGitVersion

! The Visual Studio project sets the path for where to find the header file with version info
#ifdef GIT_INCLUDE_FILE
#include GIT_INCLUDE_FILE
#endif

#ifdef GIT_VERSION_INFO
   QueryGitVersion = GIT_VERSION_INFO
#else
   QueryGitVersion = 'unversioned'
#endif

   RETURN
END FUNCTION QueryGitVersion

END MODULE
