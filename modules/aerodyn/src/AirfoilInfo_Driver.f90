!**********************************************************************************************************************************
! AFI_Driver: This code tests a stand-alone version of the AFI module
!..................................................................................................................................
! LICENSING
! Copyright (C) 2018  Envision Energy
!
!    This file is part of AirfoilInfo.
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

   
   
program AFI_Driver

   use NWTC_Library
   use VersionInfo
   use AirfoilInfo
   use AirfoilInfo_Types
   
   implicit none

   
   TYPE(ProgDesc), PARAMETER :: AFI_Ver = ProgDesc( 'AFI_driver', '', '' )
   
   
    ! Variables
   type(AFI_InitInputType)                       :: AFI_InitInputs       ! Input data for initialization
   integer, parameter                            :: NumAFI = 1;
   type(AFI_ParameterType)                       :: AFI_p(NumAFI)        ! Parameters
   type(AFI_OutputType)                          :: AFI_interp           ! interpolated AFI output values
   integer(IntKi)                                :: ErrStat              ! Status of error message
   character(1024)                               :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   character(1024)                               :: afName
   character(1024)                               :: outFileName
   integer                                       :: unOutFile(3) = -1
   character(*), parameter                       :: NumFmt = 'ES16.9E2'
   character(*), parameter                       :: Frmt = '(1x,'//NumFmt//')'
   character(*), parameter                       :: Ext(3) = (/ '.cl','.cd','.cm' /)
   character(*), parameter                       :: delim = ' '

   integer     :: i, j, k, iFile
   
   real(ReKi), allocatable                       :: Re(:)
   real(ReKi)                                    :: alpha
   real(ReKi), parameter                         :: UserProp = 0.0_ReKi


      ! Initialize the NWTC library
   call NWTC_Init(AFI_Ver%Name)
   
      ! Initialize error handling variables
   ErrMsg  = ''
   ErrStat = ErrID_None
   

   CALL DispNVD(AFI_Ver)
   
   
      ! Check for command line arguments.
   afName = ''  ! default name for input file
   CALL CheckArgs( afName )

   CALL GetRoot( afName, outFileName )
   outFileName = trim(outFileName)//'.interp.out'


      ! Setup Airfoil InitInput data structure (should come from an input file):
   AFI_InitInputs%AFTabMod    = AFITable_1 ! AFITable_2Re !
   AFI_InitInputs%InCol_Alfa  = 1
   AFI_InitInputs%InCol_Cl    = 2
   AFI_InitInputs%InCol_Cd    = 3
   AFI_InitInputs%InCol_Cm    = 4
   AFI_InitInputs%InCol_Cpmin = 0
   AFI_InitInputs%FileName    = afName
   AFI_InitInputs%UAMod       = UA_Gonzalez

   
      ! Write UA parameters to file:
   call AFI_WrHeader(delim, trim(outFileName)//'.AFI.sum', unOutFile(1), ErrStat, ErrMsg)
   
   do iFile=1,NumAFI
      if (NumAFI > 1) then
         if (iFile < 10) then
            AFI_InitInputs%FileName    = 'af00'//trim(num2lstr(iFile))//'.dat'
         else
            AFI_InitInputs%FileName    = 'af0'//trim(num2lstr(iFile))//'.dat'
         end if
      end if
   
         ! Initialize the Airfoil Info Params
      call AFI_Init ( AFI_InitInputs, AFI_p(iFile), ErrStat, ErrMsg )
         call checkError()
   
      if (ErrStat < AbortErrLev) then
         call AFI_WrData(iFile, unOutFile(1), delim, AFI_p(iFile))
         
         call AFI_WrTables(AFI_p(iFile), AFI_InitInputs%UAMod, trim(AFI_InitInputs%FileName) )
         
      end if
   end do
   close(unOutFile(1))
         
   iFile = 1
      ! allocate Re array, based on Re in the tables
   call allocAry( Re, AFI_p(iFile)%NumTabs*2 + 1, 'Re', ErrStat, ErrMsg )
      call checkError()
      
   Re(1) = AFI_p(iFile)%Table(1)%Re / 2.0_ReKi
   do i = 1, AFI_p(iFile)%NumTabs
      Re(2*i) = AFI_p(iFile)%Table(i)%Re
   end do
   do i = 1, AFI_p(iFile)%NumTabs-1
      Re(2*i+1) = (AFI_p(iFile)%Table(i)%Re + AFI_p(iFile)%Table(i+1)%Re)/2.0_ReKi
   end do
   Re(size(Re)) = AFI_p(iFile)%Table(AFI_p(iFile)%NumTabs)%Re * 2
   

   ! ------------
   do k=1,size(Ext)
      call GetNewUnit( unOutFile(k) )
      call OpenFOutFile ( unOutFile(k), trim(outFileName)//Ext(k), errStat, errMsg )
         call checkError()
      
      write( unOutFile(k), '('//trim(num2lstr(size(Re)+1))//Frmt//')' ) NaN, Re
   end do

            
      ! time marching loop
   do i = -180, 180 ! alpha
      alpha = Real(i, ReKi) ! degrees

      do k=1,size(Ext)
         write(unOutFile(k), Frmt, ADVANCE='no')  alpha
      end do
      
      alpha = alpha*pi/180.0_ReKi ! radians
      
      do j=1, size(Re)
      
         call AFI_ComputeAirfoilCoefs( alpha, Re(j), UserProp, AFI_p(iFile), AFI_interp, ErrStat, ErrMsg)
            call checkError()
      
         write(unOutFile(1), Frmt, ADVANCE='no')  AFI_interp%cl
         write(unOutFile(2), Frmt, ADVANCE='no')  AFI_interp%cd
         write(unOutFile(3), Frmt, ADVANCE='no')  AFI_interp%cm

      end do
      
      do k=1,size(Ext)
         write (unOutFile(k),'()', IOSTAT=ErrStat)          ! write the line return
      end do
      
   end do
      
   
   !-------------------------------------------------------------------------------------------------
   ! Close our output files
   !-------------------------------------------------------------------------------------------------
   

   call Cleanup()
   call NormStop()
   
   contains
   
   !====================================================================================================
   subroutine Cleanup()
   !     The routine closes any open files.
   !----------------------------------------------------------------------------------------------------  
      integer :: ie
   
      do ie=1,size(Ext)
         if (unOutFile(ie) > 0) close( unOutFile(ie), IOSTAT = ErrStat )
      end do
      
   end subroutine Cleanup

   !----------------------------------------------------------------------------------------------------  
   subroutine checkError()
      
      if (ErrStat >= AbortErrLev) then
         
         call Cleanup()
         call ProgAbort(ErrMsg)
            
      elseif ( ErrStat /= ErrID_None ) then
         
         call WrScr( trim(ErrMsg) )
            
      end if
      
   end subroutine checkError
   !----------------------------------------------------------------------------------------------------  
   
   
end program AFI_Driver

