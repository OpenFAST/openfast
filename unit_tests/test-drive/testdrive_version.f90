! This file is part of test-drive.
! SPDX-Identifier: Apache-2.0 OR MIT
!
! Licensed under either of Apache License, Version 2.0 or MIT license
! at your option; you may not use this file except in compliance with
! the License.
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module testdrive_version
   implicit none
   private

   public :: testdrive_version_string, testdrive_version_compact
   public :: get_testdrive_version


   !> String representation of the test-drive version
   character(len=*), parameter :: testdrive_version_string = "0.4.0"

   !> Numeric representation of the test-drive version
   integer, parameter :: testdrive_version_compact(3) = [0, 4, 0]


contains


!> Getter function to retrieve test-drive version
subroutine get_testdrive_version(major, minor, patch, string)

   !> Major version number of the test-drive version
   integer, intent(out), optional :: major

   !> Minor version number of the test-drive version
   integer, intent(out), optional :: minor

   !> Patch version number of the test-drive version
   integer, intent(out), optional :: patch

   !> String representation of the test-drive version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = testdrive_version_compact(1)
   end if
   if (present(minor)) then
      minor = testdrive_version_compact(2)
   end if
   if (present(patch)) then
      patch = testdrive_version_compact(3)
   end if
   if (present(string)) then
      string = testdrive_version_string
   end if

end subroutine get_testdrive_version


end module testdrive_version
