!**********************************************************************************************************************************
! The ElastoDyn.f90 and  ElastoDyn_Types.f90 make up the ElastoDyn module of the
! FAST Modularization Framework. ElastoDyn_Types is auto-generated based on FAST_Registry.txt.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2014  National Renewable Energy Laboratory
!
!    This file is part of ElastoDyn.
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
! File last committed: $Date: 2015-11-12 13:43:42 -0700 (Thu, 12 Nov 2015) $
! (File) Revision #: $Rev: 1172 $
! URL: $HeadURL: https://windsvn.nrel.gov/FAST/branches/BJonkman/Source/ElastoDyn.f90 $
!**********************************************************************************************************************************

MODULE ElastoDyn_Parameters

      ! This module contains definitions of compile-time PARAMETERS for the StrucyDyn module.
      ! Every variable defined here MUST have the PARAMETER attribute.


   USE NWTC_Library

   TYPE(ProgDesc), PARAMETER  :: ED_Ver = ProgDesc( 'ElastoDyn', 'v1.03.01a-bjj', '5-Nov-2015' )
   CHARACTER(*),   PARAMETER  :: ED_Nickname = 'ED'
   
   REAL(ReKi), PARAMETER      :: SmallAngleLimit_Deg  =  15.0                     ! Largest input angle considered "small" (used as a check on input data), degrees


      ! Parameters related to degrees of freedom (formerly MODULE DOFs)

   INTEGER(IntKi), PARAMETER        :: MaxBl    =  3                                   ! Maximum number of blades allowed in simulation
   INTEGER(IntKi), PARAMETER        :: NumBE    =  1                                   ! Number of blade-edge modes
   INTEGER(IntKi), PARAMETER        :: NumBF    =  2                                   ! Number of blade-flap modes

   INTEGER(IntKi), PARAMETER        :: DOF_Sg   =  1                                   ! DOF index for platform surge
   INTEGER(IntKi), PARAMETER        :: DOF_Sw   =  2                                   ! DOF index for platform sway
   INTEGER(IntKi), PARAMETER        :: DOF_Hv   =  3                                   ! DOF index for platform heave
   INTEGER(IntKi), PARAMETER        :: DOF_R    =  4                                   ! DOF index for platform roll
   INTEGER(IntKi), PARAMETER        :: DOF_P    =  5                                   ! DOF index for platform pitch
   INTEGER(IntKi), PARAMETER        :: DOF_Y    =  6                                   ! DOF index for platform yaw
   INTEGER(IntKi), PARAMETER        :: DOF_TFA1 =  7                                   ! DOF index for 1st tower fore-aft mode
   INTEGER(IntKi), PARAMETER        :: DOF_TSS1 =  8                                   ! DOF index for 1st tower side-to-side mode
   INTEGER(IntKi), PARAMETER        :: DOF_TFA2 =  9                                   ! DOF index for 2nd tower fore-aft mode
   INTEGER(IntKi), PARAMETER        :: DOF_TSS2 = 10                                   ! DOF index for 2nd tower side-to-side mode
   INTEGER(IntKi), PARAMETER        :: DOF_Yaw  = 11                                   ! DOF index for nacelle-yaw
   INTEGER(IntKi), PARAMETER        :: DOF_RFrl = 12                                   ! DOF index for rotor-furl
   INTEGER(IntKi), PARAMETER        :: DOF_GeAz = 13                                   ! DOF index for the generator azimuth
   INTEGER(IntKi), PARAMETER        :: DOF_DrTr = 14                                   ! DOF index for drivetrain rotational-flexibility
   INTEGER(IntKi), PARAMETER        :: DOF_TFrl = 15                                   ! DOF index for tail-furl

   INTEGER(IntKi), PARAMETER        :: DOF_BE (MaxBl,NumBE) = RESHAPE(  &              ! DOF indices for blade edge:
                                               (/ 17, 20, 23 /),   (/MaxBl,NumBE/) )   !    1st blade edge mode for blades 1,2, and 3, respectively 17 + 3*(K-1)
   INTEGER(IntKi), PARAMETER        :: DOF_BF (MaxBl,NumBF) = RESHAPE(  &              ! DOF indices for blade flap:
                                               (/ 16, 19, 22,           &              !    1st blade flap mode for blades 1,2, and 3, respectively 16 + 3*(K-1)
                                                  18, 21, 24 /),   (/MaxBl,NumBF/) )   !    2nd blade flap mode for blades 1,2, and 3, respectively 18 + 3*(K-1)


   INTEGER(IntKi), PARAMETER        :: DOF_Teet = 22 !DOF_TFrl + 2*(NumBE+NumBF)+ 1    ! DOF index for rotor-teeter



   INTEGER(IntKi), PARAMETER        :: NPA      =  9                                   ! Number of DOFs that contribute to the angular velocity of the tail (body A) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: NPB      =  7                                   ! Number of DOFs that contribute to the angular velocity of the tower top / baseplate (body B) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: NPF      =  7                                   ! Number of DOFs that contribute to the angular velocity of the tower elements (body F) in the inertia frame                                           (body F) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: NPG      = 10                                   ! Number of DOFs that contribute to the angular velocity of the generator (body G) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: NPL      = 11                                   ! Number of DOFs that contribute to the angular velocity of the low-speed shaft (body L) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: NPN      =  8                                   ! Number of DOFs that contribute to the angular velocity of the nacelle (body N) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: NPR      =  9                                   ! Number of DOFs that contribute to the angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: NPX      =  3                                   ! Number of DOFs that contribute to the angular velocity of the platform (body X) in the inertia frame.

   INTEGER(IntKi), PARAMETER        :: PX(NPX)  = (/ DOF_R, DOF_P, DOF_Y /)                                                                                          ! Array of DOF indices (pointers) that contribute to the angular velocity of the platform                                                  (body X) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: PF(NPF)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2 /)                                                  ! Array of DOF indices (pointers) that contribute to the angular velocity of the tower elements                                            (body F) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: PB(NPB)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2 /)                                                  ! Array of DOF indices (pointers) that contribute to the angular velocity of the tower top / baseplate                                     (body B) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: PN(NPN)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw /)                                         ! Array of DOF indices (pointers) that contribute to the angular velocity of the nacelle                                                   (body N) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: PR(NPR)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl /)                               ! Array of DOF indices (pointers) that contribute to the angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: PL(NPL)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz, DOF_DrTr /)           ! Array of DOF indices (pointers) that contribute to the angular velocity of the low-speed shaft                                           (body L) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: PG(NPG)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz /)                     ! Array of DOF indices (pointers) that contribute to the angular velocity of the generator                                                 (body G) in the inertia frame.
   INTEGER(IntKi), PARAMETER        :: PA(NPA)  = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_TFrl /)                               ! Array of DOF indices (pointers) that contribute to the angular velocity of the tail                                                      (body A) in the inertia frame.


      ! Parameters related to coupling scheme

   INTEGER(IntKi), PARAMETER        :: Method_RK4  = 1                                 
   INTEGER(IntKi), PARAMETER        :: Method_AB4  = 2                                 
   INTEGER(IntKi), PARAMETER        :: Method_ABM4 = 3


   INTEGER(IntKi), PARAMETER        :: PolyOrd  =  6                                    ! Order of the polynomial describing the mode shape



! ===================================================================================================
! NOTE: The following lines of code were generated by a Matlab script called "Write_ChckOutLst.m"
!      using the parameters listed in the "OutListParameters.xlsx" Excel file. Any changes to these 
!      lines should be modified in the Matlab script and/or Excel worksheet as necessary. 
! ===================================================================================================
! This code was generated by Write_ChckOutLst.m at 02-Mar-2015 10:37:31.


     ! Parameters related to output length (number of characters allowed in the output data headers):

   INTEGER(IntKi), PARAMETER      :: OutStrLenM1 = ChanLen - 1


     ! Indices for computing output channels:
     ! NOTES: 
     !    (1) These parameters are in the order stored in "OutListParameters.xlsx"
     !    (2) Array y%AllOuts() must be dimensioned to the value of the largest output parameter

     !  Time: 

   INTEGER(IntKi), PARAMETER      :: Time      =   0


     ! Blade 1 Tip Motions:

   INTEGER(IntKi), PARAMETER      :: TipDxc1   =   1
   INTEGER(IntKi), PARAMETER      :: TipDyc1   =   2
   INTEGER(IntKi), PARAMETER      :: TipDzc1   =   3
   INTEGER(IntKi), PARAMETER      :: TipDxb1   =   4
   INTEGER(IntKi), PARAMETER      :: TipDyb1   =   5
   INTEGER(IntKi), PARAMETER      :: TipALxb1  =   6
   INTEGER(IntKi), PARAMETER      :: TipALyb1  =   7
   INTEGER(IntKi), PARAMETER      :: TipALzb1  =   8
   INTEGER(IntKi), PARAMETER      :: TipRDxb1  =   9
   INTEGER(IntKi), PARAMETER      :: TipRDyb1  =  10
   INTEGER(IntKi), PARAMETER      :: TipRDzc1  =  11
   INTEGER(IntKi), PARAMETER      :: TipClrnc1 =  12


     ! Blade 2 Tip Motions:

   INTEGER(IntKi), PARAMETER      :: TipDxc2   =  13
   INTEGER(IntKi), PARAMETER      :: TipDyc2   =  14
   INTEGER(IntKi), PARAMETER      :: TipDzc2   =  15
   INTEGER(IntKi), PARAMETER      :: TipDxb2   =  16
   INTEGER(IntKi), PARAMETER      :: TipDyb2   =  17
   INTEGER(IntKi), PARAMETER      :: TipALxb2  =  18
   INTEGER(IntKi), PARAMETER      :: TipALyb2  =  19
   INTEGER(IntKi), PARAMETER      :: TipALzb2  =  20
   INTEGER(IntKi), PARAMETER      :: TipRDxb2  =  21
   INTEGER(IntKi), PARAMETER      :: TipRDyb2  =  22
   INTEGER(IntKi), PARAMETER      :: TipRDzc2  =  23
   INTEGER(IntKi), PARAMETER      :: TipClrnc2 =  24


     ! Blade 3 Tip Motions:

   INTEGER(IntKi), PARAMETER      :: TipDxc3   =  25
   INTEGER(IntKi), PARAMETER      :: TipDyc3   =  26
   INTEGER(IntKi), PARAMETER      :: TipDzc3   =  27
   INTEGER(IntKi), PARAMETER      :: TipDxb3   =  28
   INTEGER(IntKi), PARAMETER      :: TipDyb3   =  29
   INTEGER(IntKi), PARAMETER      :: TipALxb3  =  30
   INTEGER(IntKi), PARAMETER      :: TipALyb3  =  31
   INTEGER(IntKi), PARAMETER      :: TipALzb3  =  32
   INTEGER(IntKi), PARAMETER      :: TipRDxb3  =  33
   INTEGER(IntKi), PARAMETER      :: TipRDyb3  =  34
   INTEGER(IntKi), PARAMETER      :: TipRDzc3  =  35
   INTEGER(IntKi), PARAMETER      :: TipClrnc3 =  36


     ! Blade 1 Local Span Motions:

   INTEGER(IntKi), PARAMETER      :: Spn1ALxb1 =  37
   INTEGER(IntKi), PARAMETER      :: Spn1ALyb1 =  38
   INTEGER(IntKi), PARAMETER      :: Spn1ALzb1 =  39
   INTEGER(IntKi), PARAMETER      :: Spn2ALxb1 =  40
   INTEGER(IntKi), PARAMETER      :: Spn2ALyb1 =  41
   INTEGER(IntKi), PARAMETER      :: Spn2ALzb1 =  42
   INTEGER(IntKi), PARAMETER      :: Spn3ALxb1 =  43
   INTEGER(IntKi), PARAMETER      :: Spn3ALyb1 =  44
   INTEGER(IntKi), PARAMETER      :: Spn3ALzb1 =  45
   INTEGER(IntKi), PARAMETER      :: Spn4ALxb1 =  46
   INTEGER(IntKi), PARAMETER      :: Spn4ALyb1 =  47
   INTEGER(IntKi), PARAMETER      :: Spn4ALzb1 =  48
   INTEGER(IntKi), PARAMETER      :: Spn5ALxb1 =  49
   INTEGER(IntKi), PARAMETER      :: Spn5ALyb1 =  50
   INTEGER(IntKi), PARAMETER      :: Spn5ALzb1 =  51
   INTEGER(IntKi), PARAMETER      :: Spn6ALxb1 =  52
   INTEGER(IntKi), PARAMETER      :: Spn6ALyb1 =  53
   INTEGER(IntKi), PARAMETER      :: Spn6ALzb1 =  54
   INTEGER(IntKi), PARAMETER      :: Spn7ALxb1 =  55
   INTEGER(IntKi), PARAMETER      :: Spn7ALyb1 =  56
   INTEGER(IntKi), PARAMETER      :: Spn7ALzb1 =  57
   INTEGER(IntKi), PARAMETER      :: Spn8ALxb1 =  58
   INTEGER(IntKi), PARAMETER      :: Spn8ALyb1 =  59
   INTEGER(IntKi), PARAMETER      :: Spn8ALzb1 =  60
   INTEGER(IntKi), PARAMETER      :: Spn9ALxb1 =  61
   INTEGER(IntKi), PARAMETER      :: Spn9ALyb1 =  62
   INTEGER(IntKi), PARAMETER      :: Spn9ALzb1 =  63
   INTEGER(IntKi), PARAMETER      :: Spn1TDxb1 =  64
   INTEGER(IntKi), PARAMETER      :: Spn1TDyb1 =  65
   INTEGER(IntKi), PARAMETER      :: Spn1TDzb1 =  66
   INTEGER(IntKi), PARAMETER      :: Spn2TDxb1 =  67
   INTEGER(IntKi), PARAMETER      :: Spn2TDyb1 =  68
   INTEGER(IntKi), PARAMETER      :: Spn2TDzb1 =  69
   INTEGER(IntKi), PARAMETER      :: Spn3TDxb1 =  70
   INTEGER(IntKi), PARAMETER      :: Spn3TDyb1 =  71
   INTEGER(IntKi), PARAMETER      :: Spn3TDzb1 =  72
   INTEGER(IntKi), PARAMETER      :: Spn4TDxb1 =  73
   INTEGER(IntKi), PARAMETER      :: Spn4TDyb1 =  74
   INTEGER(IntKi), PARAMETER      :: Spn4TDzb1 =  75
   INTEGER(IntKi), PARAMETER      :: Spn5TDxb1 =  76
   INTEGER(IntKi), PARAMETER      :: Spn5TDyb1 =  77
   INTEGER(IntKi), PARAMETER      :: Spn5TDzb1 =  78
   INTEGER(IntKi), PARAMETER      :: Spn6TDxb1 =  79
   INTEGER(IntKi), PARAMETER      :: Spn6TDyb1 =  80
   INTEGER(IntKi), PARAMETER      :: Spn6TDzb1 =  81
   INTEGER(IntKi), PARAMETER      :: Spn7TDxb1 =  82
   INTEGER(IntKi), PARAMETER      :: Spn7TDyb1 =  83
   INTEGER(IntKi), PARAMETER      :: Spn7TDzb1 =  84
   INTEGER(IntKi), PARAMETER      :: Spn8TDxb1 =  85
   INTEGER(IntKi), PARAMETER      :: Spn8TDyb1 =  86
   INTEGER(IntKi), PARAMETER      :: Spn8TDzb1 =  87
   INTEGER(IntKi), PARAMETER      :: Spn9TDxb1 =  88
   INTEGER(IntKi), PARAMETER      :: Spn9TDyb1 =  89
   INTEGER(IntKi), PARAMETER      :: Spn9TDzb1 =  90
   INTEGER(IntKi), PARAMETER      :: Spn1RDxb1 =  91
   INTEGER(IntKi), PARAMETER      :: Spn1RDyb1 =  92
   INTEGER(IntKi), PARAMETER      :: Spn1RDzb1 =  93
   INTEGER(IntKi), PARAMETER      :: Spn2RDxb1 =  94
   INTEGER(IntKi), PARAMETER      :: Spn2RDyb1 =  95
   INTEGER(IntKi), PARAMETER      :: Spn2RDzb1 =  96
   INTEGER(IntKi), PARAMETER      :: Spn3RDxb1 =  97
   INTEGER(IntKi), PARAMETER      :: Spn3RDyb1 =  98
   INTEGER(IntKi), PARAMETER      :: Spn3RDzb1 =  99
   INTEGER(IntKi), PARAMETER      :: Spn4RDxb1 = 100
   INTEGER(IntKi), PARAMETER      :: Spn4RDyb1 = 101
   INTEGER(IntKi), PARAMETER      :: Spn4RDzb1 = 102
   INTEGER(IntKi), PARAMETER      :: Spn5RDxb1 = 103
   INTEGER(IntKi), PARAMETER      :: Spn5RDyb1 = 104
   INTEGER(IntKi), PARAMETER      :: Spn5RDzb1 = 105
   INTEGER(IntKi), PARAMETER      :: Spn6RDxb1 = 106
   INTEGER(IntKi), PARAMETER      :: Spn6RDyb1 = 107
   INTEGER(IntKi), PARAMETER      :: Spn6RDzb1 = 108
   INTEGER(IntKi), PARAMETER      :: Spn7RDxb1 = 109
   INTEGER(IntKi), PARAMETER      :: Spn7RDyb1 = 110
   INTEGER(IntKi), PARAMETER      :: Spn7RDzb1 = 111
   INTEGER(IntKi), PARAMETER      :: Spn8RDxb1 = 112
   INTEGER(IntKi), PARAMETER      :: Spn8RDyb1 = 113
   INTEGER(IntKi), PARAMETER      :: Spn8RDzb1 = 114
   INTEGER(IntKi), PARAMETER      :: Spn9RDxb1 = 115
   INTEGER(IntKi), PARAMETER      :: Spn9RDyb1 = 116
   INTEGER(IntKi), PARAMETER      :: Spn9RDzb1 = 117


     ! Blade 2 Local Span Motions:

   INTEGER(IntKi), PARAMETER      :: Spn1ALxb2 = 118
   INTEGER(IntKi), PARAMETER      :: Spn1ALyb2 = 119
   INTEGER(IntKi), PARAMETER      :: Spn1ALzb2 = 120
   INTEGER(IntKi), PARAMETER      :: Spn2ALxb2 = 121
   INTEGER(IntKi), PARAMETER      :: Spn2ALyb2 = 122
   INTEGER(IntKi), PARAMETER      :: Spn2ALzb2 = 123
   INTEGER(IntKi), PARAMETER      :: Spn3ALxb2 = 124
   INTEGER(IntKi), PARAMETER      :: Spn3ALyb2 = 125
   INTEGER(IntKi), PARAMETER      :: Spn3ALzb2 = 126
   INTEGER(IntKi), PARAMETER      :: Spn4ALxb2 = 127
   INTEGER(IntKi), PARAMETER      :: Spn4ALyb2 = 128
   INTEGER(IntKi), PARAMETER      :: Spn4ALzb2 = 129
   INTEGER(IntKi), PARAMETER      :: Spn5ALxb2 = 130
   INTEGER(IntKi), PARAMETER      :: Spn5ALyb2 = 131
   INTEGER(IntKi), PARAMETER      :: Spn5ALzb2 = 132
   INTEGER(IntKi), PARAMETER      :: Spn6ALxb2 = 133
   INTEGER(IntKi), PARAMETER      :: Spn6ALyb2 = 134
   INTEGER(IntKi), PARAMETER      :: Spn6ALzb2 = 135
   INTEGER(IntKi), PARAMETER      :: Spn7ALxb2 = 136
   INTEGER(IntKi), PARAMETER      :: Spn7ALyb2 = 137
   INTEGER(IntKi), PARAMETER      :: Spn7ALzb2 = 138
   INTEGER(IntKi), PARAMETER      :: Spn8ALxb2 = 139
   INTEGER(IntKi), PARAMETER      :: Spn8ALyb2 = 140
   INTEGER(IntKi), PARAMETER      :: Spn8ALzb2 = 141
   INTEGER(IntKi), PARAMETER      :: Spn9ALxb2 = 142
   INTEGER(IntKi), PARAMETER      :: Spn9ALyb2 = 143
   INTEGER(IntKi), PARAMETER      :: Spn9ALzb2 = 144
   INTEGER(IntKi), PARAMETER      :: Spn1TDxb2 = 145
   INTEGER(IntKi), PARAMETER      :: Spn1TDyb2 = 146
   INTEGER(IntKi), PARAMETER      :: Spn1TDzb2 = 147
   INTEGER(IntKi), PARAMETER      :: Spn2TDxb2 = 148
   INTEGER(IntKi), PARAMETER      :: Spn2TDyb2 = 149
   INTEGER(IntKi), PARAMETER      :: Spn2TDzb2 = 150
   INTEGER(IntKi), PARAMETER      :: Spn3TDxb2 = 151
   INTEGER(IntKi), PARAMETER      :: Spn3TDyb2 = 152
   INTEGER(IntKi), PARAMETER      :: Spn3TDzb2 = 153
   INTEGER(IntKi), PARAMETER      :: Spn4TDxb2 = 154
   INTEGER(IntKi), PARAMETER      :: Spn4TDyb2 = 155
   INTEGER(IntKi), PARAMETER      :: Spn4TDzb2 = 156
   INTEGER(IntKi), PARAMETER      :: Spn5TDxb2 = 157
   INTEGER(IntKi), PARAMETER      :: Spn5TDyb2 = 158
   INTEGER(IntKi), PARAMETER      :: Spn5TDzb2 = 159
   INTEGER(IntKi), PARAMETER      :: Spn6TDxb2 = 160
   INTEGER(IntKi), PARAMETER      :: Spn6TDyb2 = 161
   INTEGER(IntKi), PARAMETER      :: Spn6TDzb2 = 162
   INTEGER(IntKi), PARAMETER      :: Spn7TDxb2 = 163
   INTEGER(IntKi), PARAMETER      :: Spn7TDyb2 = 164
   INTEGER(IntKi), PARAMETER      :: Spn7TDzb2 = 165
   INTEGER(IntKi), PARAMETER      :: Spn8TDxb2 = 166
   INTEGER(IntKi), PARAMETER      :: Spn8TDyb2 = 167
   INTEGER(IntKi), PARAMETER      :: Spn8TDzb2 = 168
   INTEGER(IntKi), PARAMETER      :: Spn9TDxb2 = 169
   INTEGER(IntKi), PARAMETER      :: Spn9TDyb2 = 170
   INTEGER(IntKi), PARAMETER      :: Spn9TDzb2 = 171
   INTEGER(IntKi), PARAMETER      :: Spn1RDxb2 = 172
   INTEGER(IntKi), PARAMETER      :: Spn1RDyb2 = 173
   INTEGER(IntKi), PARAMETER      :: Spn1RDzb2 = 174
   INTEGER(IntKi), PARAMETER      :: Spn2RDxb2 = 175
   INTEGER(IntKi), PARAMETER      :: Spn2RDyb2 = 176
   INTEGER(IntKi), PARAMETER      :: Spn2RDzb2 = 177
   INTEGER(IntKi), PARAMETER      :: Spn3RDxb2 = 178
   INTEGER(IntKi), PARAMETER      :: Spn3RDyb2 = 179
   INTEGER(IntKi), PARAMETER      :: Spn3RDzb2 = 180
   INTEGER(IntKi), PARAMETER      :: Spn4RDxb2 = 181
   INTEGER(IntKi), PARAMETER      :: Spn4RDyb2 = 182
   INTEGER(IntKi), PARAMETER      :: Spn4RDzb2 = 183
   INTEGER(IntKi), PARAMETER      :: Spn5RDxb2 = 184
   INTEGER(IntKi), PARAMETER      :: Spn5RDyb2 = 185
   INTEGER(IntKi), PARAMETER      :: Spn5RDzb2 = 186
   INTEGER(IntKi), PARAMETER      :: Spn6RDxb2 = 187
   INTEGER(IntKi), PARAMETER      :: Spn6RDyb2 = 188
   INTEGER(IntKi), PARAMETER      :: Spn6RDzb2 = 189
   INTEGER(IntKi), PARAMETER      :: Spn7RDxb2 = 190
   INTEGER(IntKi), PARAMETER      :: Spn7RDyb2 = 191
   INTEGER(IntKi), PARAMETER      :: Spn7RDzb2 = 192
   INTEGER(IntKi), PARAMETER      :: Spn8RDxb2 = 193
   INTEGER(IntKi), PARAMETER      :: Spn8RDyb2 = 194
   INTEGER(IntKi), PARAMETER      :: Spn8RDzb2 = 195
   INTEGER(IntKi), PARAMETER      :: Spn9RDxb2 = 196
   INTEGER(IntKi), PARAMETER      :: Spn9RDyb2 = 197
   INTEGER(IntKi), PARAMETER      :: Spn9RDzb2 = 198


     ! Blade 3 Local Span Motions:

   INTEGER(IntKi), PARAMETER      :: Spn1ALxb3 = 199
   INTEGER(IntKi), PARAMETER      :: Spn1ALyb3 = 200
   INTEGER(IntKi), PARAMETER      :: Spn1ALzb3 = 201
   INTEGER(IntKi), PARAMETER      :: Spn2ALxb3 = 202
   INTEGER(IntKi), PARAMETER      :: Spn2ALyb3 = 203
   INTEGER(IntKi), PARAMETER      :: Spn2ALzb3 = 204
   INTEGER(IntKi), PARAMETER      :: Spn3ALxb3 = 205
   INTEGER(IntKi), PARAMETER      :: Spn3ALyb3 = 206
   INTEGER(IntKi), PARAMETER      :: Spn3ALzb3 = 207
   INTEGER(IntKi), PARAMETER      :: Spn4ALxb3 = 208
   INTEGER(IntKi), PARAMETER      :: Spn4ALyb3 = 209
   INTEGER(IntKi), PARAMETER      :: Spn4ALzb3 = 210
   INTEGER(IntKi), PARAMETER      :: Spn5ALxb3 = 211
   INTEGER(IntKi), PARAMETER      :: Spn5ALyb3 = 212
   INTEGER(IntKi), PARAMETER      :: Spn5ALzb3 = 213
   INTEGER(IntKi), PARAMETER      :: Spn6ALxb3 = 214
   INTEGER(IntKi), PARAMETER      :: Spn6ALyb3 = 215
   INTEGER(IntKi), PARAMETER      :: Spn6ALzb3 = 216
   INTEGER(IntKi), PARAMETER      :: Spn7ALxb3 = 217
   INTEGER(IntKi), PARAMETER      :: Spn7ALyb3 = 218
   INTEGER(IntKi), PARAMETER      :: Spn7ALzb3 = 219
   INTEGER(IntKi), PARAMETER      :: Spn8ALxb3 = 220
   INTEGER(IntKi), PARAMETER      :: Spn8ALyb3 = 221
   INTEGER(IntKi), PARAMETER      :: Spn8ALzb3 = 222
   INTEGER(IntKi), PARAMETER      :: Spn9ALxb3 = 223
   INTEGER(IntKi), PARAMETER      :: Spn9ALyb3 = 224
   INTEGER(IntKi), PARAMETER      :: Spn9ALzb3 = 225
   INTEGER(IntKi), PARAMETER      :: Spn1TDxb3 = 226
   INTEGER(IntKi), PARAMETER      :: Spn1TDyb3 = 227
   INTEGER(IntKi), PARAMETER      :: Spn1TDzb3 = 228
   INTEGER(IntKi), PARAMETER      :: Spn2TDxb3 = 229
   INTEGER(IntKi), PARAMETER      :: Spn2TDyb3 = 230
   INTEGER(IntKi), PARAMETER      :: Spn2TDzb3 = 231
   INTEGER(IntKi), PARAMETER      :: Spn3TDxb3 = 232
   INTEGER(IntKi), PARAMETER      :: Spn3TDyb3 = 233
   INTEGER(IntKi), PARAMETER      :: Spn3TDzb3 = 234
   INTEGER(IntKi), PARAMETER      :: Spn4TDxb3 = 235
   INTEGER(IntKi), PARAMETER      :: Spn4TDyb3 = 236
   INTEGER(IntKi), PARAMETER      :: Spn4TDzb3 = 237
   INTEGER(IntKi), PARAMETER      :: Spn5TDxb3 = 238
   INTEGER(IntKi), PARAMETER      :: Spn5TDyb3 = 239
   INTEGER(IntKi), PARAMETER      :: Spn5TDzb3 = 240
   INTEGER(IntKi), PARAMETER      :: Spn6TDxb3 = 241
   INTEGER(IntKi), PARAMETER      :: Spn6TDyb3 = 242
   INTEGER(IntKi), PARAMETER      :: Spn6TDzb3 = 243
   INTEGER(IntKi), PARAMETER      :: Spn7TDxb3 = 244
   INTEGER(IntKi), PARAMETER      :: Spn7TDyb3 = 245
   INTEGER(IntKi), PARAMETER      :: Spn7TDzb3 = 246
   INTEGER(IntKi), PARAMETER      :: Spn8TDxb3 = 247
   INTEGER(IntKi), PARAMETER      :: Spn8TDyb3 = 248
   INTEGER(IntKi), PARAMETER      :: Spn8TDzb3 = 249
   INTEGER(IntKi), PARAMETER      :: Spn9TDxb3 = 250
   INTEGER(IntKi), PARAMETER      :: Spn9TDyb3 = 251
   INTEGER(IntKi), PARAMETER      :: Spn9TDzb3 = 252
   INTEGER(IntKi), PARAMETER      :: Spn1RDxb3 = 253
   INTEGER(IntKi), PARAMETER      :: Spn1RDyb3 = 254
   INTEGER(IntKi), PARAMETER      :: Spn1RDzb3 = 255
   INTEGER(IntKi), PARAMETER      :: Spn2RDxb3 = 256
   INTEGER(IntKi), PARAMETER      :: Spn2RDyb3 = 257
   INTEGER(IntKi), PARAMETER      :: Spn2RDzb3 = 258
   INTEGER(IntKi), PARAMETER      :: Spn3RDxb3 = 259
   INTEGER(IntKi), PARAMETER      :: Spn3RDyb3 = 260
   INTEGER(IntKi), PARAMETER      :: Spn3RDzb3 = 261
   INTEGER(IntKi), PARAMETER      :: Spn4RDxb3 = 262
   INTEGER(IntKi), PARAMETER      :: Spn4RDyb3 = 263
   INTEGER(IntKi), PARAMETER      :: Spn4RDzb3 = 264
   INTEGER(IntKi), PARAMETER      :: Spn5RDxb3 = 265
   INTEGER(IntKi), PARAMETER      :: Spn5RDyb3 = 266
   INTEGER(IntKi), PARAMETER      :: Spn5RDzb3 = 267
   INTEGER(IntKi), PARAMETER      :: Spn6RDxb3 = 268
   INTEGER(IntKi), PARAMETER      :: Spn6RDyb3 = 269
   INTEGER(IntKi), PARAMETER      :: Spn6RDzb3 = 270
   INTEGER(IntKi), PARAMETER      :: Spn7RDxb3 = 271
   INTEGER(IntKi), PARAMETER      :: Spn7RDyb3 = 272
   INTEGER(IntKi), PARAMETER      :: Spn7RDzb3 = 273
   INTEGER(IntKi), PARAMETER      :: Spn8RDxb3 = 274
   INTEGER(IntKi), PARAMETER      :: Spn8RDyb3 = 275
   INTEGER(IntKi), PARAMETER      :: Spn8RDzb3 = 276
   INTEGER(IntKi), PARAMETER      :: Spn9RDxb3 = 277
   INTEGER(IntKi), PARAMETER      :: Spn9RDyb3 = 278
   INTEGER(IntKi), PARAMETER      :: Spn9RDzb3 = 279


     ! Blade Pitch Motions:

   INTEGER(IntKi), PARAMETER      :: PtchPMzc1 = 280
   INTEGER(IntKi), PARAMETER      :: PtchPMzc2 = 281
   INTEGER(IntKi), PARAMETER      :: PtchPMzc3 = 282


     ! Teeter Motions:

   INTEGER(IntKi), PARAMETER      :: TeetPya   = 283
   INTEGER(IntKi), PARAMETER      :: TeetVya   = 284
   INTEGER(IntKi), PARAMETER      :: TeetAya   = 285


     ! Shaft Motions:

   INTEGER(IntKi), PARAMETER      :: LSSTipPxa = 286
   INTEGER(IntKi), PARAMETER      :: LSSTipVxa = 287
   INTEGER(IntKi), PARAMETER      :: LSSTipAxa = 288
   INTEGER(IntKi), PARAMETER      :: LSSGagPxa = 289
   INTEGER(IntKi), PARAMETER      :: LSSGagVxa = 290
   INTEGER(IntKi), PARAMETER      :: LSSGagAxa = 291
   INTEGER(IntKi), PARAMETER      :: HSShftV   = 292
   INTEGER(IntKi), PARAMETER      :: HSShftA   = 293


     ! Nacelle IMU Motions:

   INTEGER(IntKi), PARAMETER      :: NcIMUTVxs = 294
   INTEGER(IntKi), PARAMETER      :: NcIMUTVys = 295
   INTEGER(IntKi), PARAMETER      :: NcIMUTVzs = 296
   INTEGER(IntKi), PARAMETER      :: NcIMUTAxs = 297
   INTEGER(IntKi), PARAMETER      :: NcIMUTAys = 298
   INTEGER(IntKi), PARAMETER      :: NcIMUTAzs = 299
   INTEGER(IntKi), PARAMETER      :: NcIMURVxs = 300
   INTEGER(IntKi), PARAMETER      :: NcIMURVys = 301
   INTEGER(IntKi), PARAMETER      :: NcIMURVzs = 302
   INTEGER(IntKi), PARAMETER      :: NcIMURAxs = 303
   INTEGER(IntKi), PARAMETER      :: NcIMURAys = 304
   INTEGER(IntKi), PARAMETER      :: NcIMURAzs = 305


     ! Rotor-Furl Motions:

   INTEGER(IntKi), PARAMETER      :: RotFurlP  = 306
   INTEGER(IntKi), PARAMETER      :: RotFurlV  = 307
   INTEGER(IntKi), PARAMETER      :: RotFurlA  = 308


     ! Tail-Furl Motions:

   INTEGER(IntKi), PARAMETER      :: TailFurlP = 309
   INTEGER(IntKi), PARAMETER      :: TailFurlV = 310
   INTEGER(IntKi), PARAMETER      :: TailFurlA = 311


     ! Nacelle Yaw Motions:

   INTEGER(IntKi), PARAMETER      :: YawPzn    = 312
   INTEGER(IntKi), PARAMETER      :: YawVzn    = 313
   INTEGER(IntKi), PARAMETER      :: YawAzn    = 314


     ! Tower-Top / Yaw Bearing Motions:

   INTEGER(IntKi), PARAMETER      :: YawBrTDxp = 315
   INTEGER(IntKi), PARAMETER      :: YawBrTDyp = 316
   INTEGER(IntKi), PARAMETER      :: YawBrTDzp = 317
   INTEGER(IntKi), PARAMETER      :: YawBrTDxt = 318
   INTEGER(IntKi), PARAMETER      :: YawBrTDyt = 319
   INTEGER(IntKi), PARAMETER      :: YawBrTDzt = 320
   INTEGER(IntKi), PARAMETER      :: YawBrTAxp = 321
   INTEGER(IntKi), PARAMETER      :: YawBrTAyp = 322
   INTEGER(IntKi), PARAMETER      :: YawBrTAzp = 323
   INTEGER(IntKi), PARAMETER      :: YawBrRDxt = 324
   INTEGER(IntKi), PARAMETER      :: YawBrRDyt = 325
   INTEGER(IntKi), PARAMETER      :: YawBrRDzt = 326
   INTEGER(IntKi), PARAMETER      :: YawBrRVxp = 327
   INTEGER(IntKi), PARAMETER      :: YawBrRVyp = 328
   INTEGER(IntKi), PARAMETER      :: YawBrRVzp = 329
   INTEGER(IntKi), PARAMETER      :: YawBrRAxp = 330
   INTEGER(IntKi), PARAMETER      :: YawBrRAyp = 331
   INTEGER(IntKi), PARAMETER      :: YawBrRAzp = 332


     ! Local Tower Motions:

   INTEGER(IntKi), PARAMETER      :: TwHt1ALxt = 333
   INTEGER(IntKi), PARAMETER      :: TwHt1ALyt = 334
   INTEGER(IntKi), PARAMETER      :: TwHt1ALzt = 335
   INTEGER(IntKi), PARAMETER      :: TwHt2ALxt = 336
   INTEGER(IntKi), PARAMETER      :: TwHt2ALyt = 337
   INTEGER(IntKi), PARAMETER      :: TwHt2ALzt = 338
   INTEGER(IntKi), PARAMETER      :: TwHt3ALxt = 339
   INTEGER(IntKi), PARAMETER      :: TwHt3ALyt = 340
   INTEGER(IntKi), PARAMETER      :: TwHt3ALzt = 341
   INTEGER(IntKi), PARAMETER      :: TwHt4ALxt = 342
   INTEGER(IntKi), PARAMETER      :: TwHt4ALyt = 343
   INTEGER(IntKi), PARAMETER      :: TwHt4ALzt = 344
   INTEGER(IntKi), PARAMETER      :: TwHt5ALxt = 345
   INTEGER(IntKi), PARAMETER      :: TwHt5ALyt = 346
   INTEGER(IntKi), PARAMETER      :: TwHt5ALzt = 347
   INTEGER(IntKi), PARAMETER      :: TwHt6ALxt = 348
   INTEGER(IntKi), PARAMETER      :: TwHt6ALyt = 349
   INTEGER(IntKi), PARAMETER      :: TwHt6ALzt = 350
   INTEGER(IntKi), PARAMETER      :: TwHt7ALxt = 351
   INTEGER(IntKi), PARAMETER      :: TwHt7ALyt = 352
   INTEGER(IntKi), PARAMETER      :: TwHt7ALzt = 353
   INTEGER(IntKi), PARAMETER      :: TwHt8ALxt = 354
   INTEGER(IntKi), PARAMETER      :: TwHt8ALyt = 355
   INTEGER(IntKi), PARAMETER      :: TwHt8ALzt = 356
   INTEGER(IntKi), PARAMETER      :: TwHt9ALxt = 357
   INTEGER(IntKi), PARAMETER      :: TwHt9ALyt = 358
   INTEGER(IntKi), PARAMETER      :: TwHt9ALzt = 359
   INTEGER(IntKi), PARAMETER      :: TwHt1TDxt = 360
   INTEGER(IntKi), PARAMETER      :: TwHt1TDyt = 361
   INTEGER(IntKi), PARAMETER      :: TwHt1TDzt = 362
   INTEGER(IntKi), PARAMETER      :: TwHt2TDxt = 363
   INTEGER(IntKi), PARAMETER      :: TwHt2TDyt = 364
   INTEGER(IntKi), PARAMETER      :: TwHt2TDzt = 365
   INTEGER(IntKi), PARAMETER      :: TwHt3TDxt = 366
   INTEGER(IntKi), PARAMETER      :: TwHt3TDyt = 367
   INTEGER(IntKi), PARAMETER      :: TwHt3TDzt = 368
   INTEGER(IntKi), PARAMETER      :: TwHt4TDxt = 369
   INTEGER(IntKi), PARAMETER      :: TwHt4TDyt = 370
   INTEGER(IntKi), PARAMETER      :: TwHt4TDzt = 371
   INTEGER(IntKi), PARAMETER      :: TwHt5TDxt = 372
   INTEGER(IntKi), PARAMETER      :: TwHt5TDyt = 373
   INTEGER(IntKi), PARAMETER      :: TwHt5TDzt = 374
   INTEGER(IntKi), PARAMETER      :: TwHt6TDxt = 375
   INTEGER(IntKi), PARAMETER      :: TwHt6TDyt = 376
   INTEGER(IntKi), PARAMETER      :: TwHt6TDzt = 377
   INTEGER(IntKi), PARAMETER      :: TwHt7TDxt = 378
   INTEGER(IntKi), PARAMETER      :: TwHt7TDyt = 379
   INTEGER(IntKi), PARAMETER      :: TwHt7TDzt = 380
   INTEGER(IntKi), PARAMETER      :: TwHt8TDxt = 381
   INTEGER(IntKi), PARAMETER      :: TwHt8TDyt = 382
   INTEGER(IntKi), PARAMETER      :: TwHt8TDzt = 383
   INTEGER(IntKi), PARAMETER      :: TwHt9TDxt = 384
   INTEGER(IntKi), PARAMETER      :: TwHt9TDyt = 385
   INTEGER(IntKi), PARAMETER      :: TwHt9TDzt = 386
   INTEGER(IntKi), PARAMETER      :: TwHt1RDxt = 387
   INTEGER(IntKi), PARAMETER      :: TwHt1RDyt = 388
   INTEGER(IntKi), PARAMETER      :: TwHt1RDzt = 389
   INTEGER(IntKi), PARAMETER      :: TwHt2RDxt = 390
   INTEGER(IntKi), PARAMETER      :: TwHt2RDyt = 391
   INTEGER(IntKi), PARAMETER      :: TwHt2RDzt = 392
   INTEGER(IntKi), PARAMETER      :: TwHt3RDxt = 393
   INTEGER(IntKi), PARAMETER      :: TwHt3RDyt = 394
   INTEGER(IntKi), PARAMETER      :: TwHt3RDzt = 395
   INTEGER(IntKi), PARAMETER      :: TwHt4RDxt = 396
   INTEGER(IntKi), PARAMETER      :: TwHt4RDyt = 397
   INTEGER(IntKi), PARAMETER      :: TwHt4RDzt = 398
   INTEGER(IntKi), PARAMETER      :: TwHt5RDxt = 399
   INTEGER(IntKi), PARAMETER      :: TwHt5RDyt = 400
   INTEGER(IntKi), PARAMETER      :: TwHt5RDzt = 401
   INTEGER(IntKi), PARAMETER      :: TwHt6RDxt = 402
   INTEGER(IntKi), PARAMETER      :: TwHt6RDyt = 403
   INTEGER(IntKi), PARAMETER      :: TwHt6RDzt = 404
   INTEGER(IntKi), PARAMETER      :: TwHt7RDxt = 405
   INTEGER(IntKi), PARAMETER      :: TwHt7RDyt = 406
   INTEGER(IntKi), PARAMETER      :: TwHt7RDzt = 407
   INTEGER(IntKi), PARAMETER      :: TwHt8RDxt = 408
   INTEGER(IntKi), PARAMETER      :: TwHt8RDyt = 409
   INTEGER(IntKi), PARAMETER      :: TwHt8RDzt = 410
   INTEGER(IntKi), PARAMETER      :: TwHt9RDxt = 411
   INTEGER(IntKi), PARAMETER      :: TwHt9RDyt = 412
   INTEGER(IntKi), PARAMETER      :: TwHt9RDzt = 413
   INTEGER(IntKi), PARAMETER      :: TwHt1TPxi = 414
   INTEGER(IntKi), PARAMETER      :: TwHt1TPyi = 415
   INTEGER(IntKi), PARAMETER      :: TwHt1TPzi = 416
   INTEGER(IntKi), PARAMETER      :: TwHt2TPxi = 417
   INTEGER(IntKi), PARAMETER      :: TwHt2TPyi = 418
   INTEGER(IntKi), PARAMETER      :: TwHt2TPzi = 419
   INTEGER(IntKi), PARAMETER      :: TwHt3TPxi = 420
   INTEGER(IntKi), PARAMETER      :: TwHt3TPyi = 421
   INTEGER(IntKi), PARAMETER      :: TwHt3TPzi = 422
   INTEGER(IntKi), PARAMETER      :: TwHt4TPxi = 423
   INTEGER(IntKi), PARAMETER      :: TwHt4TPyi = 424
   INTEGER(IntKi), PARAMETER      :: TwHt4TPzi = 425
   INTEGER(IntKi), PARAMETER      :: TwHt5TPxi = 426
   INTEGER(IntKi), PARAMETER      :: TwHt5TPyi = 427
   INTEGER(IntKi), PARAMETER      :: TwHt5TPzi = 428
   INTEGER(IntKi), PARAMETER      :: TwHt6TPxi = 429
   INTEGER(IntKi), PARAMETER      :: TwHt6TPyi = 430
   INTEGER(IntKi), PARAMETER      :: TwHt6TPzi = 431
   INTEGER(IntKi), PARAMETER      :: TwHt7TPxi = 432
   INTEGER(IntKi), PARAMETER      :: TwHt7TPyi = 433
   INTEGER(IntKi), PARAMETER      :: TwHt7TPzi = 434
   INTEGER(IntKi), PARAMETER      :: TwHt8TPxi = 435
   INTEGER(IntKi), PARAMETER      :: TwHt8TPyi = 436
   INTEGER(IntKi), PARAMETER      :: TwHt8TPzi = 437
   INTEGER(IntKi), PARAMETER      :: TwHt9TPxi = 438
   INTEGER(IntKi), PARAMETER      :: TwHt9TPyi = 439
   INTEGER(IntKi), PARAMETER      :: TwHt9TPzi = 440
   INTEGER(IntKi), PARAMETER      :: TwHt1RPxi = 441
   INTEGER(IntKi), PARAMETER      :: TwHt1RPyi = 442
   INTEGER(IntKi), PARAMETER      :: TwHt1RPzi = 443
   INTEGER(IntKi), PARAMETER      :: TwHt2RPxi = 444
   INTEGER(IntKi), PARAMETER      :: TwHt2RPyi = 445
   INTEGER(IntKi), PARAMETER      :: TwHt2RPzi = 446
   INTEGER(IntKi), PARAMETER      :: TwHt3RPxi = 447
   INTEGER(IntKi), PARAMETER      :: TwHt3RPyi = 448
   INTEGER(IntKi), PARAMETER      :: TwHt3RPzi = 449
   INTEGER(IntKi), PARAMETER      :: TwHt4RPxi = 450
   INTEGER(IntKi), PARAMETER      :: TwHt4RPyi = 451
   INTEGER(IntKi), PARAMETER      :: TwHt4RPzi = 452
   INTEGER(IntKi), PARAMETER      :: TwHt5RPxi = 453
   INTEGER(IntKi), PARAMETER      :: TwHt5RPyi = 454
   INTEGER(IntKi), PARAMETER      :: TwHt5RPzi = 455
   INTEGER(IntKi), PARAMETER      :: TwHt6RPxi = 456
   INTEGER(IntKi), PARAMETER      :: TwHt6RPyi = 457
   INTEGER(IntKi), PARAMETER      :: TwHt6RPzi = 458
   INTEGER(IntKi), PARAMETER      :: TwHt7RPxi = 459
   INTEGER(IntKi), PARAMETER      :: TwHt7RPyi = 460
   INTEGER(IntKi), PARAMETER      :: TwHt7RPzi = 461
   INTEGER(IntKi), PARAMETER      :: TwHt8RPxi = 462
   INTEGER(IntKi), PARAMETER      :: TwHt8RPyi = 463
   INTEGER(IntKi), PARAMETER      :: TwHt8RPzi = 464
   INTEGER(IntKi), PARAMETER      :: TwHt9RPxi = 465
   INTEGER(IntKi), PARAMETER      :: TwHt9RPyi = 466
   INTEGER(IntKi), PARAMETER      :: TwHt9RPzi = 467


     ! Platform Motions:

   INTEGER(IntKi), PARAMETER      :: PtfmTDxt  = 468
   INTEGER(IntKi), PARAMETER      :: PtfmTDyt  = 469
   INTEGER(IntKi), PARAMETER      :: PtfmTDzt  = 470
   INTEGER(IntKi), PARAMETER      :: PtfmTDxi  = 471
   INTEGER(IntKi), PARAMETER      :: PtfmTDyi  = 472
   INTEGER(IntKi), PARAMETER      :: PtfmTDzi  = 473
   INTEGER(IntKi), PARAMETER      :: PtfmTVxt  = 474
   INTEGER(IntKi), PARAMETER      :: PtfmTVyt  = 475
   INTEGER(IntKi), PARAMETER      :: PtfmTVzt  = 476
   INTEGER(IntKi), PARAMETER      :: PtfmTVxi  = 477
   INTEGER(IntKi), PARAMETER      :: PtfmTVyi  = 478
   INTEGER(IntKi), PARAMETER      :: PtfmTVzi  = 479
   INTEGER(IntKi), PARAMETER      :: PtfmTAxt  = 480
   INTEGER(IntKi), PARAMETER      :: PtfmTAyt  = 481
   INTEGER(IntKi), PARAMETER      :: PtfmTAzt  = 482
   INTEGER(IntKi), PARAMETER      :: PtfmTAxi  = 483
   INTEGER(IntKi), PARAMETER      :: PtfmTAyi  = 484
   INTEGER(IntKi), PARAMETER      :: PtfmTAzi  = 485
   INTEGER(IntKi), PARAMETER      :: PtfmRDxi  = 486
   INTEGER(IntKi), PARAMETER      :: PtfmRDyi  = 487
   INTEGER(IntKi), PARAMETER      :: PtfmRDzi  = 488
   INTEGER(IntKi), PARAMETER      :: PtfmRVxt  = 489
   INTEGER(IntKi), PARAMETER      :: PtfmRVyt  = 490
   INTEGER(IntKi), PARAMETER      :: PtfmRVzt  = 491
   INTEGER(IntKi), PARAMETER      :: PtfmRVxi  = 492
   INTEGER(IntKi), PARAMETER      :: PtfmRVyi  = 493
   INTEGER(IntKi), PARAMETER      :: PtfmRVzi  = 494
   INTEGER(IntKi), PARAMETER      :: PtfmRAxt  = 495
   INTEGER(IntKi), PARAMETER      :: PtfmRAyt  = 496
   INTEGER(IntKi), PARAMETER      :: PtfmRAzt  = 497
   INTEGER(IntKi), PARAMETER      :: PtfmRAxi  = 498
   INTEGER(IntKi), PARAMETER      :: PtfmRAyi  = 499
   INTEGER(IntKi), PARAMETER      :: PtfmRAzi  = 500


     ! Blade 1 Root Loads:

   INTEGER(IntKi), PARAMETER      :: RootFxc1  = 501
   INTEGER(IntKi), PARAMETER      :: RootFyc1  = 502
   INTEGER(IntKi), PARAMETER      :: RootFzc1  = 503
   INTEGER(IntKi), PARAMETER      :: RootFxb1  = 504
   INTEGER(IntKi), PARAMETER      :: RootFyb1  = 505
   INTEGER(IntKi), PARAMETER      :: RootMxc1  = 506
   INTEGER(IntKi), PARAMETER      :: RootMyc1  = 507
   INTEGER(IntKi), PARAMETER      :: RootMzc1  = 508
   INTEGER(IntKi), PARAMETER      :: RootMxb1  = 509
   INTEGER(IntKi), PARAMETER      :: RootMyb1  = 510


     ! Blade 2 Root Loads:

   INTEGER(IntKi), PARAMETER      :: RootFxc2  = 511
   INTEGER(IntKi), PARAMETER      :: RootFyc2  = 512
   INTEGER(IntKi), PARAMETER      :: RootFzc2  = 513
   INTEGER(IntKi), PARAMETER      :: RootFxb2  = 514
   INTEGER(IntKi), PARAMETER      :: RootFyb2  = 515
   INTEGER(IntKi), PARAMETER      :: RootMxc2  = 516
   INTEGER(IntKi), PARAMETER      :: RootMyc2  = 517
   INTEGER(IntKi), PARAMETER      :: RootMzc2  = 518
   INTEGER(IntKi), PARAMETER      :: RootMxb2  = 519
   INTEGER(IntKi), PARAMETER      :: RootMyb2  = 520


     ! Blade 3 Root Loads:

   INTEGER(IntKi), PARAMETER      :: RootFxc3  = 521
   INTEGER(IntKi), PARAMETER      :: RootFyc3  = 522
   INTEGER(IntKi), PARAMETER      :: RootFzc3  = 523
   INTEGER(IntKi), PARAMETER      :: RootFxb3  = 524
   INTEGER(IntKi), PARAMETER      :: RootFyb3  = 525
   INTEGER(IntKi), PARAMETER      :: RootMxc3  = 526
   INTEGER(IntKi), PARAMETER      :: RootMyc3  = 527
   INTEGER(IntKi), PARAMETER      :: RootMzc3  = 528
   INTEGER(IntKi), PARAMETER      :: RootMxb3  = 529
   INTEGER(IntKi), PARAMETER      :: RootMyb3  = 530


     ! Blade 1 Local Span Loads:

   INTEGER(IntKi), PARAMETER      :: Spn1MLxb1 = 531
   INTEGER(IntKi), PARAMETER      :: Spn1MLyb1 = 532
   INTEGER(IntKi), PARAMETER      :: Spn1MLzb1 = 533
   INTEGER(IntKi), PARAMETER      :: Spn2MLxb1 = 534
   INTEGER(IntKi), PARAMETER      :: Spn2MLyb1 = 535
   INTEGER(IntKi), PARAMETER      :: Spn2MLzb1 = 536
   INTEGER(IntKi), PARAMETER      :: Spn3MLxb1 = 537
   INTEGER(IntKi), PARAMETER      :: Spn3MLyb1 = 538
   INTEGER(IntKi), PARAMETER      :: Spn3MLzb1 = 539
   INTEGER(IntKi), PARAMETER      :: Spn4MLxb1 = 540
   INTEGER(IntKi), PARAMETER      :: Spn4MLyb1 = 541
   INTEGER(IntKi), PARAMETER      :: Spn4MLzb1 = 542
   INTEGER(IntKi), PARAMETER      :: Spn5MLxb1 = 543
   INTEGER(IntKi), PARAMETER      :: Spn5MLyb1 = 544
   INTEGER(IntKi), PARAMETER      :: Spn5MLzb1 = 545
   INTEGER(IntKi), PARAMETER      :: Spn6MLxb1 = 546
   INTEGER(IntKi), PARAMETER      :: Spn6MLyb1 = 547
   INTEGER(IntKi), PARAMETER      :: Spn6MLzb1 = 548
   INTEGER(IntKi), PARAMETER      :: Spn7MLxb1 = 549
   INTEGER(IntKi), PARAMETER      :: Spn7MLyb1 = 550
   INTEGER(IntKi), PARAMETER      :: Spn7MLzb1 = 551
   INTEGER(IntKi), PARAMETER      :: Spn8MLxb1 = 552
   INTEGER(IntKi), PARAMETER      :: Spn8MLyb1 = 553
   INTEGER(IntKi), PARAMETER      :: Spn8MLzb1 = 554
   INTEGER(IntKi), PARAMETER      :: Spn9MLxb1 = 555
   INTEGER(IntKi), PARAMETER      :: Spn9MLyb1 = 556
   INTEGER(IntKi), PARAMETER      :: Spn9MLzb1 = 557
   INTEGER(IntKi), PARAMETER      :: Spn1FLxb1 = 558
   INTEGER(IntKi), PARAMETER      :: Spn1FLyb1 = 559
   INTEGER(IntKi), PARAMETER      :: Spn1FLzb1 = 560
   INTEGER(IntKi), PARAMETER      :: Spn2FLxb1 = 561
   INTEGER(IntKi), PARAMETER      :: Spn2FLyb1 = 562
   INTEGER(IntKi), PARAMETER      :: Spn2FLzb1 = 563
   INTEGER(IntKi), PARAMETER      :: Spn3FLxb1 = 564
   INTEGER(IntKi), PARAMETER      :: Spn3FLyb1 = 565
   INTEGER(IntKi), PARAMETER      :: Spn3FLzb1 = 566
   INTEGER(IntKi), PARAMETER      :: Spn4FLxb1 = 567
   INTEGER(IntKi), PARAMETER      :: Spn4FLyb1 = 568
   INTEGER(IntKi), PARAMETER      :: Spn4FLzb1 = 569
   INTEGER(IntKi), PARAMETER      :: Spn5FLxb1 = 570
   INTEGER(IntKi), PARAMETER      :: Spn5FLyb1 = 571
   INTEGER(IntKi), PARAMETER      :: Spn5FLzb1 = 572
   INTEGER(IntKi), PARAMETER      :: Spn6FLxb1 = 573
   INTEGER(IntKi), PARAMETER      :: Spn6FLyb1 = 574
   INTEGER(IntKi), PARAMETER      :: Spn6FLzb1 = 575
   INTEGER(IntKi), PARAMETER      :: Spn7FLxb1 = 576
   INTEGER(IntKi), PARAMETER      :: Spn7FLyb1 = 577
   INTEGER(IntKi), PARAMETER      :: Spn7FLzb1 = 578
   INTEGER(IntKi), PARAMETER      :: Spn8FLxb1 = 579
   INTEGER(IntKi), PARAMETER      :: Spn8FLyb1 = 580
   INTEGER(IntKi), PARAMETER      :: Spn8FLzb1 = 581
   INTEGER(IntKi), PARAMETER      :: Spn9FLxb1 = 582
   INTEGER(IntKi), PARAMETER      :: Spn9FLyb1 = 583
   INTEGER(IntKi), PARAMETER      :: Spn9FLzb1 = 584


     ! Blade 2 Local Span Loads:

   INTEGER(IntKi), PARAMETER      :: Spn1MLxb2 = 585
   INTEGER(IntKi), PARAMETER      :: Spn1MLyb2 = 586
   INTEGER(IntKi), PARAMETER      :: Spn1MLzb2 = 587
   INTEGER(IntKi), PARAMETER      :: Spn2MLxb2 = 588
   INTEGER(IntKi), PARAMETER      :: Spn2MLyb2 = 589
   INTEGER(IntKi), PARAMETER      :: Spn2MLzb2 = 590
   INTEGER(IntKi), PARAMETER      :: Spn3MLxb2 = 591
   INTEGER(IntKi), PARAMETER      :: Spn3MLyb2 = 592
   INTEGER(IntKi), PARAMETER      :: Spn3MLzb2 = 593
   INTEGER(IntKi), PARAMETER      :: Spn4MLxb2 = 594
   INTEGER(IntKi), PARAMETER      :: Spn4MLyb2 = 595
   INTEGER(IntKi), PARAMETER      :: Spn4MLzb2 = 596
   INTEGER(IntKi), PARAMETER      :: Spn5MLxb2 = 597
   INTEGER(IntKi), PARAMETER      :: Spn5MLyb2 = 598
   INTEGER(IntKi), PARAMETER      :: Spn5MLzb2 = 599
   INTEGER(IntKi), PARAMETER      :: Spn6MLxb2 = 600
   INTEGER(IntKi), PARAMETER      :: Spn6MLyb2 = 601
   INTEGER(IntKi), PARAMETER      :: Spn6MLzb2 = 602
   INTEGER(IntKi), PARAMETER      :: Spn7MLxb2 = 603
   INTEGER(IntKi), PARAMETER      :: Spn7MLyb2 = 604
   INTEGER(IntKi), PARAMETER      :: Spn7MLzb2 = 605
   INTEGER(IntKi), PARAMETER      :: Spn8MLxb2 = 606
   INTEGER(IntKi), PARAMETER      :: Spn8MLyb2 = 607
   INTEGER(IntKi), PARAMETER      :: Spn8MLzb2 = 608
   INTEGER(IntKi), PARAMETER      :: Spn9MLxb2 = 609
   INTEGER(IntKi), PARAMETER      :: Spn9MLyb2 = 610
   INTEGER(IntKi), PARAMETER      :: Spn9MLzb2 = 611
   INTEGER(IntKi), PARAMETER      :: Spn1FLxb2 = 612
   INTEGER(IntKi), PARAMETER      :: Spn1FLyb2 = 613
   INTEGER(IntKi), PARAMETER      :: Spn1FLzb2 = 614
   INTEGER(IntKi), PARAMETER      :: Spn2FLxb2 = 615
   INTEGER(IntKi), PARAMETER      :: Spn2FLyb2 = 616
   INTEGER(IntKi), PARAMETER      :: Spn2FLzb2 = 617
   INTEGER(IntKi), PARAMETER      :: Spn3FLxb2 = 618
   INTEGER(IntKi), PARAMETER      :: Spn3FLyb2 = 619
   INTEGER(IntKi), PARAMETER      :: Spn3FLzb2 = 620
   INTEGER(IntKi), PARAMETER      :: Spn4FLxb2 = 621
   INTEGER(IntKi), PARAMETER      :: Spn4FLyb2 = 622
   INTEGER(IntKi), PARAMETER      :: Spn4FLzb2 = 623
   INTEGER(IntKi), PARAMETER      :: Spn5FLxb2 = 624
   INTEGER(IntKi), PARAMETER      :: Spn5FLyb2 = 625
   INTEGER(IntKi), PARAMETER      :: Spn5FLzb2 = 626
   INTEGER(IntKi), PARAMETER      :: Spn6FLxb2 = 627
   INTEGER(IntKi), PARAMETER      :: Spn6FLyb2 = 628
   INTEGER(IntKi), PARAMETER      :: Spn6FLzb2 = 629
   INTEGER(IntKi), PARAMETER      :: Spn7FLxb2 = 630
   INTEGER(IntKi), PARAMETER      :: Spn7FLyb2 = 631
   INTEGER(IntKi), PARAMETER      :: Spn7FLzb2 = 632
   INTEGER(IntKi), PARAMETER      :: Spn8FLxb2 = 633
   INTEGER(IntKi), PARAMETER      :: Spn8FLyb2 = 634
   INTEGER(IntKi), PARAMETER      :: Spn8FLzb2 = 635
   INTEGER(IntKi), PARAMETER      :: Spn9FLxb2 = 636
   INTEGER(IntKi), PARAMETER      :: Spn9FLyb2 = 637
   INTEGER(IntKi), PARAMETER      :: Spn9FLzb2 = 638


     ! Blade 3 Local Span Loads:

   INTEGER(IntKi), PARAMETER      :: Spn1MLxb3 = 639
   INTEGER(IntKi), PARAMETER      :: Spn1MLyb3 = 640
   INTEGER(IntKi), PARAMETER      :: Spn1MLzb3 = 641
   INTEGER(IntKi), PARAMETER      :: Spn2MLxb3 = 642
   INTEGER(IntKi), PARAMETER      :: Spn2MLyb3 = 643
   INTEGER(IntKi), PARAMETER      :: Spn2MLzb3 = 644
   INTEGER(IntKi), PARAMETER      :: Spn3MLxb3 = 645
   INTEGER(IntKi), PARAMETER      :: Spn3MLyb3 = 646
   INTEGER(IntKi), PARAMETER      :: Spn3MLzb3 = 647
   INTEGER(IntKi), PARAMETER      :: Spn4MLxb3 = 648
   INTEGER(IntKi), PARAMETER      :: Spn4MLyb3 = 649
   INTEGER(IntKi), PARAMETER      :: Spn4MLzb3 = 650
   INTEGER(IntKi), PARAMETER      :: Spn5MLxb3 = 651
   INTEGER(IntKi), PARAMETER      :: Spn5MLyb3 = 652
   INTEGER(IntKi), PARAMETER      :: Spn5MLzb3 = 653
   INTEGER(IntKi), PARAMETER      :: Spn6MLxb3 = 654
   INTEGER(IntKi), PARAMETER      :: Spn6MLyb3 = 655
   INTEGER(IntKi), PARAMETER      :: Spn6MLzb3 = 656
   INTEGER(IntKi), PARAMETER      :: Spn7MLxb3 = 657
   INTEGER(IntKi), PARAMETER      :: Spn7MLyb3 = 658
   INTEGER(IntKi), PARAMETER      :: Spn7MLzb3 = 659
   INTEGER(IntKi), PARAMETER      :: Spn8MLxb3 = 660
   INTEGER(IntKi), PARAMETER      :: Spn8MLyb3 = 661
   INTEGER(IntKi), PARAMETER      :: Spn8MLzb3 = 662
   INTEGER(IntKi), PARAMETER      :: Spn9MLxb3 = 663
   INTEGER(IntKi), PARAMETER      :: Spn9MLyb3 = 664
   INTEGER(IntKi), PARAMETER      :: Spn9MLzb3 = 665
   INTEGER(IntKi), PARAMETER      :: Spn1FLxb3 = 666
   INTEGER(IntKi), PARAMETER      :: Spn1FLyb3 = 667
   INTEGER(IntKi), PARAMETER      :: Spn1FLzb3 = 668
   INTEGER(IntKi), PARAMETER      :: Spn2FLxb3 = 669
   INTEGER(IntKi), PARAMETER      :: Spn2FLyb3 = 670
   INTEGER(IntKi), PARAMETER      :: Spn2FLzb3 = 671
   INTEGER(IntKi), PARAMETER      :: Spn3FLxb3 = 672
   INTEGER(IntKi), PARAMETER      :: Spn3FLyb3 = 673
   INTEGER(IntKi), PARAMETER      :: Spn3FLzb3 = 674
   INTEGER(IntKi), PARAMETER      :: Spn4FLxb3 = 675
   INTEGER(IntKi), PARAMETER      :: Spn4FLyb3 = 676
   INTEGER(IntKi), PARAMETER      :: Spn4FLzb3 = 677
   INTEGER(IntKi), PARAMETER      :: Spn5FLxb3 = 678
   INTEGER(IntKi), PARAMETER      :: Spn5FLyb3 = 679
   INTEGER(IntKi), PARAMETER      :: Spn5FLzb3 = 680
   INTEGER(IntKi), PARAMETER      :: Spn6FLxb3 = 681
   INTEGER(IntKi), PARAMETER      :: Spn6FLyb3 = 682
   INTEGER(IntKi), PARAMETER      :: Spn6FLzb3 = 683
   INTEGER(IntKi), PARAMETER      :: Spn7FLxb3 = 684
   INTEGER(IntKi), PARAMETER      :: Spn7FLyb3 = 685
   INTEGER(IntKi), PARAMETER      :: Spn7FLzb3 = 686
   INTEGER(IntKi), PARAMETER      :: Spn8FLxb3 = 687
   INTEGER(IntKi), PARAMETER      :: Spn8FLyb3 = 688
   INTEGER(IntKi), PARAMETER      :: Spn8FLzb3 = 689
   INTEGER(IntKi), PARAMETER      :: Spn9FLxb3 = 690
   INTEGER(IntKi), PARAMETER      :: Spn9FLyb3 = 691
   INTEGER(IntKi), PARAMETER      :: Spn9FLzb3 = 692


     ! Hub and Rotor Loads:

   INTEGER(IntKi), PARAMETER      :: LSShftFxa = 693
   INTEGER(IntKi), PARAMETER      :: LSShftFya = 694
   INTEGER(IntKi), PARAMETER      :: LSShftFza = 695
   INTEGER(IntKi), PARAMETER      :: LSShftFys = 696
   INTEGER(IntKi), PARAMETER      :: LSShftFzs = 697
   INTEGER(IntKi), PARAMETER      :: LSShftMxa = 698
   INTEGER(IntKi), PARAMETER      :: LSSTipMya = 699
   INTEGER(IntKi), PARAMETER      :: LSSTipMza = 700
   INTEGER(IntKi), PARAMETER      :: LSSTipMys = 701
   INTEGER(IntKi), PARAMETER      :: LSSTipMzs = 702
   INTEGER(IntKi), PARAMETER      :: RotPwr    = 703


     ! Shaft Strain Gage Loads:

   INTEGER(IntKi), PARAMETER      :: LSSGagMya = 704
   INTEGER(IntKi), PARAMETER      :: LSSGagMza = 705
   INTEGER(IntKi), PARAMETER      :: LSSGagMys = 706
   INTEGER(IntKi), PARAMETER      :: LSSGagMzs = 707


     ! High-Speed Shaft Loads:

   INTEGER(IntKi), PARAMETER      :: HSShftTq  = 708
   INTEGER(IntKi), PARAMETER      :: HSSBrTq   = 709
   INTEGER(IntKi), PARAMETER      :: HSShftPwr = 710


     ! Rotor-Furl Bearing Loads:

   INTEGER(IntKi), PARAMETER      :: RFrlBrM   = 711


     ! Tail-Furl Bearing Loads:

   INTEGER(IntKi), PARAMETER      :: TFrlBrM   = 712


     ! Tower-Top / Yaw Bearing Loads:

   INTEGER(IntKi), PARAMETER      :: YawBrFxn  = 713
   INTEGER(IntKi), PARAMETER      :: YawBrFyn  = 714
   INTEGER(IntKi), PARAMETER      :: YawBrFzn  = 715
   INTEGER(IntKi), PARAMETER      :: YawBrFxp  = 716
   INTEGER(IntKi), PARAMETER      :: YawBrFyp  = 717
   INTEGER(IntKi), PARAMETER      :: YawBrMxn  = 718
   INTEGER(IntKi), PARAMETER      :: YawBrMyn  = 719
   INTEGER(IntKi), PARAMETER      :: YawBrMzn  = 720
   INTEGER(IntKi), PARAMETER      :: YawBrMxp  = 721
   INTEGER(IntKi), PARAMETER      :: YawBrMyp  = 722


     ! Tower Base Loads:

   INTEGER(IntKi), PARAMETER      :: TwrBsFxt  = 723
   INTEGER(IntKi), PARAMETER      :: TwrBsFyt  = 724
   INTEGER(IntKi), PARAMETER      :: TwrBsFzt  = 725
   INTEGER(IntKi), PARAMETER      :: TwrBsMxt  = 726
   INTEGER(IntKi), PARAMETER      :: TwrBsMyt  = 727
   INTEGER(IntKi), PARAMETER      :: TwrBsMzt  = 728


     ! Local Tower Loads:

   INTEGER(IntKi), PARAMETER      :: TwHt1MLxt = 729
   INTEGER(IntKi), PARAMETER      :: TwHt1MLyt = 730
   INTEGER(IntKi), PARAMETER      :: TwHt1MLzt = 731
   INTEGER(IntKi), PARAMETER      :: TwHt2MLxt = 732
   INTEGER(IntKi), PARAMETER      :: TwHt2MLyt = 733
   INTEGER(IntKi), PARAMETER      :: TwHt2MLzt = 734
   INTEGER(IntKi), PARAMETER      :: TwHt3MLxt = 735
   INTEGER(IntKi), PARAMETER      :: TwHt3MLyt = 736
   INTEGER(IntKi), PARAMETER      :: TwHt3MLzt = 737
   INTEGER(IntKi), PARAMETER      :: TwHt4MLxt = 738
   INTEGER(IntKi), PARAMETER      :: TwHt4MLyt = 739
   INTEGER(IntKi), PARAMETER      :: TwHt4MLzt = 740
   INTEGER(IntKi), PARAMETER      :: TwHt5MLxt = 741
   INTEGER(IntKi), PARAMETER      :: TwHt5MLyt = 742
   INTEGER(IntKi), PARAMETER      :: TwHt5MLzt = 743
   INTEGER(IntKi), PARAMETER      :: TwHt6MLxt = 744
   INTEGER(IntKi), PARAMETER      :: TwHt6MLyt = 745
   INTEGER(IntKi), PARAMETER      :: TwHt6MLzt = 746
   INTEGER(IntKi), PARAMETER      :: TwHt7MLxt = 747
   INTEGER(IntKi), PARAMETER      :: TwHt7MLyt = 748
   INTEGER(IntKi), PARAMETER      :: TwHt7MLzt = 749
   INTEGER(IntKi), PARAMETER      :: TwHt8MLxt = 750
   INTEGER(IntKi), PARAMETER      :: TwHt8MLyt = 751
   INTEGER(IntKi), PARAMETER      :: TwHt8MLzt = 752
   INTEGER(IntKi), PARAMETER      :: TwHt9MLxt = 753
   INTEGER(IntKi), PARAMETER      :: TwHt9MLyt = 754
   INTEGER(IntKi), PARAMETER      :: TwHt9MLzt = 755
   INTEGER(IntKi), PARAMETER      :: TwHt1FLxt = 756
   INTEGER(IntKi), PARAMETER      :: TwHt1FLyt = 757
   INTEGER(IntKi), PARAMETER      :: TwHt1FLzt = 758
   INTEGER(IntKi), PARAMETER      :: TwHt2FLxt = 759
   INTEGER(IntKi), PARAMETER      :: TwHt2FLyt = 760
   INTEGER(IntKi), PARAMETER      :: TwHt2FLzt = 761
   INTEGER(IntKi), PARAMETER      :: TwHt3FLxt = 762
   INTEGER(IntKi), PARAMETER      :: TwHt3FLyt = 763
   INTEGER(IntKi), PARAMETER      :: TwHt3FLzt = 764
   INTEGER(IntKi), PARAMETER      :: TwHt4FLxt = 765
   INTEGER(IntKi), PARAMETER      :: TwHt4FLyt = 766
   INTEGER(IntKi), PARAMETER      :: TwHt4FLzt = 767
   INTEGER(IntKi), PARAMETER      :: TwHt5FLxt = 768
   INTEGER(IntKi), PARAMETER      :: TwHt5FLyt = 769
   INTEGER(IntKi), PARAMETER      :: TwHt5FLzt = 770
   INTEGER(IntKi), PARAMETER      :: TwHt6FLxt = 771
   INTEGER(IntKi), PARAMETER      :: TwHt6FLyt = 772
   INTEGER(IntKi), PARAMETER      :: TwHt6FLzt = 773
   INTEGER(IntKi), PARAMETER      :: TwHt7FLxt = 774
   INTEGER(IntKi), PARAMETER      :: TwHt7FLyt = 775
   INTEGER(IntKi), PARAMETER      :: TwHt7FLzt = 776
   INTEGER(IntKi), PARAMETER      :: TwHt8FLxt = 777
   INTEGER(IntKi), PARAMETER      :: TwHt8FLyt = 778
   INTEGER(IntKi), PARAMETER      :: TwHt8FLzt = 779
   INTEGER(IntKi), PARAMETER      :: TwHt9FLxt = 780
   INTEGER(IntKi), PARAMETER      :: TwHt9FLyt = 781
   INTEGER(IntKi), PARAMETER      :: TwHt9FLzt = 782


     ! Internal Degrees of Freedom:

   INTEGER(IntKi), PARAMETER      :: Q_B1E1    = 783
   INTEGER(IntKi), PARAMETER      :: Q_B2E1    = 784
   INTEGER(IntKi), PARAMETER      :: Q_B3E1    = 785
   INTEGER(IntKi), PARAMETER      :: Q_B1F1    = 786
   INTEGER(IntKi), PARAMETER      :: Q_B2F1    = 787
   INTEGER(IntKi), PARAMETER      :: Q_B3F1    = 788
   INTEGER(IntKi), PARAMETER      :: Q_B1F2    = 789
   INTEGER(IntKi), PARAMETER      :: Q_B2F2    = 790
   INTEGER(IntKi), PARAMETER      :: Q_B3F2    = 791
   INTEGER(IntKi), PARAMETER      :: Q_Teet    = 792
   INTEGER(IntKi), PARAMETER      :: Q_DrTr    = 793
   INTEGER(IntKi), PARAMETER      :: Q_GeAz    = 794
   INTEGER(IntKi), PARAMETER      :: Q_RFrl    = 795
   INTEGER(IntKi), PARAMETER      :: Q_TFrl    = 796
   INTEGER(IntKi), PARAMETER      :: Q_Yaw     = 797
   INTEGER(IntKi), PARAMETER      :: Q_TFA1    = 798
   INTEGER(IntKi), PARAMETER      :: Q_TSS1    = 799
   INTEGER(IntKi), PARAMETER      :: Q_TFA2    = 800
   INTEGER(IntKi), PARAMETER      :: Q_TSS2    = 801
   INTEGER(IntKi), PARAMETER      :: Q_Sg      = 802
   INTEGER(IntKi), PARAMETER      :: Q_Sw      = 803
   INTEGER(IntKi), PARAMETER      :: Q_Hv      = 804
   INTEGER(IntKi), PARAMETER      :: Q_R       = 805
   INTEGER(IntKi), PARAMETER      :: Q_P       = 806
   INTEGER(IntKi), PARAMETER      :: Q_Y       = 807
   INTEGER(IntKi), PARAMETER      :: QD_B1E1   = 808
   INTEGER(IntKi), PARAMETER      :: QD_B2E1   = 809
   INTEGER(IntKi), PARAMETER      :: QD_B3E1   = 810
   INTEGER(IntKi), PARAMETER      :: QD_B1F1   = 811
   INTEGER(IntKi), PARAMETER      :: QD_B2F1   = 812
   INTEGER(IntKi), PARAMETER      :: QD_B3F1   = 813
   INTEGER(IntKi), PARAMETER      :: QD_B1F2   = 814
   INTEGER(IntKi), PARAMETER      :: QD_B2F2   = 815
   INTEGER(IntKi), PARAMETER      :: QD_B3F2   = 816
   INTEGER(IntKi), PARAMETER      :: QD_Teet   = 817
   INTEGER(IntKi), PARAMETER      :: QD_DrTr   = 818
   INTEGER(IntKi), PARAMETER      :: QD_GeAz   = 819
   INTEGER(IntKi), PARAMETER      :: QD_RFrl   = 820
   INTEGER(IntKi), PARAMETER      :: QD_TFrl   = 821
   INTEGER(IntKi), PARAMETER      :: QD_Yaw    = 822
   INTEGER(IntKi), PARAMETER      :: QD_TFA1   = 823
   INTEGER(IntKi), PARAMETER      :: QD_TSS1   = 824
   INTEGER(IntKi), PARAMETER      :: QD_TFA2   = 825
   INTEGER(IntKi), PARAMETER      :: QD_TSS2   = 826
   INTEGER(IntKi), PARAMETER      :: QD_Sg     = 827
   INTEGER(IntKi), PARAMETER      :: QD_Sw     = 828
   INTEGER(IntKi), PARAMETER      :: QD_Hv     = 829
   INTEGER(IntKi), PARAMETER      :: QD_R      = 830
   INTEGER(IntKi), PARAMETER      :: QD_P      = 831
   INTEGER(IntKi), PARAMETER      :: QD_Y      = 832
   INTEGER(IntKi), PARAMETER      :: QD2_B1E1  = 833
   INTEGER(IntKi), PARAMETER      :: QD2_B2E1  = 834
   INTEGER(IntKi), PARAMETER      :: QD2_B3E1  = 835
   INTEGER(IntKi), PARAMETER      :: QD2_B1F1  = 836
   INTEGER(IntKi), PARAMETER      :: QD2_B2F1  = 837
   INTEGER(IntKi), PARAMETER      :: QD2_B3F1  = 838
   INTEGER(IntKi), PARAMETER      :: QD2_B1F2  = 839
   INTEGER(IntKi), PARAMETER      :: QD2_B2F2  = 840
   INTEGER(IntKi), PARAMETER      :: QD2_B3F2  = 841
   INTEGER(IntKi), PARAMETER      :: QD2_Teet  = 842
   INTEGER(IntKi), PARAMETER      :: QD2_DrTr  = 843
   INTEGER(IntKi), PARAMETER      :: QD2_GeAz  = 844
   INTEGER(IntKi), PARAMETER      :: QD2_RFrl  = 845
   INTEGER(IntKi), PARAMETER      :: QD2_TFrl  = 846
   INTEGER(IntKi), PARAMETER      :: QD2_Yaw   = 847
   INTEGER(IntKi), PARAMETER      :: QD2_TFA1  = 848
   INTEGER(IntKi), PARAMETER      :: QD2_TSS1  = 849
   INTEGER(IntKi), PARAMETER      :: QD2_TFA2  = 850
   INTEGER(IntKi), PARAMETER      :: QD2_TSS2  = 851
   INTEGER(IntKi), PARAMETER      :: QD2_Sg    = 852
   INTEGER(IntKi), PARAMETER      :: QD2_Sw    = 853
   INTEGER(IntKi), PARAMETER      :: QD2_Hv    = 854
   INTEGER(IntKi), PARAMETER      :: QD2_R     = 855
   INTEGER(IntKi), PARAMETER      :: QD2_P     = 856
   INTEGER(IntKi), PARAMETER      :: QD2_Y     = 857


     ! The maximum number of output channels which can be output by the code.
   INTEGER(IntKi), PARAMETER      :: MaxOutPts = 857

!End of code generated by Matlab script
! ===================================================================================================


INTEGER,  PARAMETER          :: TipDxc( 3)  = (/TipDxc1,  TipDxc2,  TipDxc3/)
INTEGER,  PARAMETER          :: TipDyc( 3)  = (/TipDyc1,  TipDyc2,  TipDyc3/)
INTEGER,  PARAMETER          :: TipDzc( 3)  = (/TipDzc1,  TipDzc2,  TipDzc3/)
INTEGER,  PARAMETER          :: TipDxb( 3)  = (/TipDxb1,  TipDxb2,  TipDxb3/)
INTEGER,  PARAMETER          :: TipDyb( 3)  = (/TipDyb1,  TipDyb2,  TipDyb3/)
INTEGER,  PARAMETER          :: TipALxb(3)  = (/TipALxb1, TipALxb2, TipALxb3/)
INTEGER,  PARAMETER          :: TipALyb(3)  = (/TipALyb1, TipALyb2, TipALyb3/)
INTEGER,  PARAMETER          :: TipALzb(3)  = (/TipALzb1, TipALzb2, TipALzb3/)
INTEGER,  PARAMETER          :: TipRDxb(3)  = (/TipRDxb1, TipRDxb2, TipRDxb3/)
INTEGER,  PARAMETER          :: TipRDyb(3)  = (/TipRDyb1, TipRDyb2, TipRDyb3/)
INTEGER,  PARAMETER          :: TipRDzc(3)  = (/TipRDzc1, TipRDzc2, TipRDzc3/)
INTEGER,  PARAMETER          :: TipClrnc(3) = (/TipClrnc1,TipClrnc2,TipClrnc3/)
INTEGER,  PARAMETER          :: PtchPMzc(3) = (/PtchPMzc1,PtchPMzc2,PtchPMzc3/)

INTEGER,  PARAMETER          :: RootFxc(3) = (/ RootFxc1,RootFxc2,RootFxc3 /)
INTEGER,  PARAMETER          :: RootFyc(3) = (/ RootFyc1,RootFyc2,RootFyc3 /)
INTEGER,  PARAMETER          :: RootFzc(3) = (/ RootFzc1,RootFzc2,RootFzc3 /)
INTEGER,  PARAMETER          :: RootFxb(3) = (/ RootFxb1,RootFxb2,RootFxb3 /)
INTEGER,  PARAMETER          :: RootFyb(3) = (/ RootFyb1,RootFyb2,RootFyb3 /)
INTEGER,  PARAMETER          :: RootMxc(3) = (/ RootMxc1,RootMxc2,RootMxc3 /)
INTEGER,  PARAMETER          :: RootMyc(3) = (/ RootMyc1,RootMyc2,RootMyc3 /)
INTEGER,  PARAMETER          :: RootMzc(3) = (/ RootMzc1,RootMzc2,RootMzc3 /)
INTEGER,  PARAMETER          :: RootMxb(3) = (/ RootMxb1,RootMxb2,RootMxb3 /)
INTEGER,  PARAMETER          :: RootMyb(3) = (/ RootMyb1,RootMyb2,RootMyb3 /)

INTEGER,  PARAMETER          :: SpnALxb(9, 3) = RESHAPE( (/ &
                                    Spn1ALxb1,Spn2ALxb1,Spn3ALxb1,Spn4ALxb1,Spn5ALxb1,Spn6ALxb1,Spn7ALxb1,Spn8ALxb1,Spn9ALxb1, &
                                    Spn1ALxb2,Spn2ALxb2,Spn3ALxb2,Spn4ALxb2,Spn5ALxb2,Spn6ALxb2,Spn7ALxb2,Spn8ALxb2,Spn9ALxb2, &
                                    Spn1ALxb3,Spn2ALxb3,Spn3ALxb3,Spn4ALxb3,Spn5ALxb3,Spn6ALxb3,Spn7ALxb3,Spn8ALxb3,Spn9ALxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnALyb(9, 3) = RESHAPE( (/ &
                                    Spn1ALyb1,Spn2ALyb1,Spn3ALyb1,Spn4ALyb1,Spn5ALyb1,Spn6ALyb1,Spn7ALyb1,Spn8ALyb1,Spn9ALyb1, &
                                    Spn1ALyb2,Spn2ALyb2,Spn3ALyb2,Spn4ALyb2,Spn5ALyb2,Spn6ALyb2,Spn7ALyb2,Spn8ALyb2,Spn9ALyb2, &
                                    Spn1ALyb3,Spn2ALyb3,Spn3ALyb3,Spn4ALyb3,Spn5ALyb3,Spn6ALyb3,Spn7ALyb3,Spn8ALyb3,Spn9ALyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnALzb(9, 3) = RESHAPE( (/ &
                                    Spn1ALzb1,Spn2ALzb1,Spn3ALzb1,Spn4ALzb1,Spn5ALzb1,Spn6ALzb1,Spn7ALzb1,Spn8ALzb1,Spn9ALzb1, &
                                    Spn1ALzb2,Spn2ALzb2,Spn3ALzb2,Spn4ALzb2,Spn5ALzb2,Spn6ALzb2,Spn7ALzb2,Spn8ALzb2,Spn9ALzb2, &
                                    Spn1ALzb3,Spn2ALzb3,Spn3ALzb3,Spn4ALzb3,Spn5ALzb3,Spn6ALzb3,Spn7ALzb3,Spn8ALzb3,Spn9ALzb3  &
                                /), (/9, 3/) )

INTEGER,  PARAMETER          :: SpnFLxb(9,3) = RESHAPE( (/ &
                                    Spn1FLxb1,Spn2FLxb1,Spn3FLxb1,Spn4FLxb1,Spn5FLxb1,Spn6FLxb1,Spn7FLxb1,Spn8FLxb1,Spn9FLxb1, &
                                    Spn1FLxb2,Spn2FLxb2,Spn3FLxb2,Spn4FLxb2,Spn5FLxb2,Spn6FLxb2,Spn7FLxb2,Spn8FLxb2,Spn9FLxb2, &
                                    Spn1FLxb3,Spn2FLxb3,Spn3FLxb3,Spn4FLxb3,Spn5FLxb3,Spn6FLxb3,Spn7FLxb3,Spn8FLxb3,Spn9FLxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnFLyb(9,3) = RESHAPE( (/ &
                                    Spn1FLyb1,Spn2FLyb1,Spn3FLyb1,Spn4FLyb1,Spn5FLyb1,Spn6FLyb1,Spn7FLyb1,Spn8FLyb1,Spn9FLyb1, &
                                    Spn1FLyb2,Spn2FLyb2,Spn3FLyb2,Spn4FLyb2,Spn5FLyb2,Spn6FLyb2,Spn7FLyb2,Spn8FLyb2,Spn9FLyb2, &
                                    Spn1FLyb3,Spn2FLyb3,Spn3FLyb3,Spn4FLyb3,Spn5FLyb3,Spn6FLyb3,Spn7FLyb3,Spn8FLyb3,Spn9FLyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnFLzb(9,3) = RESHAPE( (/ &
                                    Spn1FLzb1,Spn2FLzb1,Spn3FLzb1,Spn4FLzb1,Spn5FLzb1,Spn6FLzb1,Spn7FLzb1,Spn8FLzb1,Spn9FLzb1, &
                                    Spn1FLzb2,Spn2FLzb2,Spn3FLzb2,Spn4FLzb2,Spn5FLzb2,Spn6FLzb2,Spn7FLzb2,Spn8FLzb2,Spn9FLzb2, &
                                    Spn1FLzb3,Spn2FLzb3,Spn3FLzb3,Spn4FLzb3,Spn5FLzb3,Spn6FLzb3,Spn7FLzb3,Spn8FLzb3,Spn9FLzb3  &
                                /), (/9, 3/) )

INTEGER,  PARAMETER          :: SpnMLxb(9,3) = RESHAPE( (/ &
                                    Spn1MLxb1,Spn2MLxb1,Spn3MLxb1,Spn4MLxb1,Spn5MLxb1,Spn6MLxb1,Spn7MLxb1,Spn8MLxb1,Spn9MLxb1, &
                                    Spn1MLxb2,Spn2MLxb2,Spn3MLxb2,Spn4MLxb2,Spn5MLxb2,Spn6MLxb2,Spn7MLxb2,Spn8MLxb2,Spn9MLxb2, &
                                    Spn1MLxb3,Spn2MLxb3,Spn3MLxb3,Spn4MLxb3,Spn5MLxb3,Spn6MLxb3,Spn7MLxb3,Spn8MLxb3,Spn9MLxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnMLyb(9,3) = RESHAPE( (/ &
                                    Spn1MLyb1,Spn2MLyb1,Spn3MLyb1,Spn4MLyb1,Spn5MLyb1,Spn6MLyb1,Spn7MLyb1,Spn8MLyb1,Spn9MLyb1, &
                                    Spn1MLyb2,Spn2MLyb2,Spn3MLyb2,Spn4MLyb2,Spn5MLyb2,Spn6MLyb2,Spn7MLyb2,Spn8MLyb2,Spn9MLyb2, &
                                    Spn1MLyb3,Spn2MLyb3,Spn3MLyb3,Spn4MLyb3,Spn5MLyb3,Spn6MLyb3,Spn7MLyb3,Spn8MLyb3,Spn9MLyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnMLzb(9,3) = RESHAPE( (/ &
                                    Spn1MLzb1,Spn2MLzb1,Spn3MLzb1,Spn4MLzb1,Spn5MLzb1,Spn6MLzb1,Spn7MLzb1,Spn8MLzb1,Spn9MLzb1, &
                                    Spn1MLzb2,Spn2MLzb2,Spn3MLzb2,Spn4MLzb2,Spn5MLzb2,Spn6MLzb2,Spn7MLzb2,Spn8MLzb2,Spn9MLzb2, &
                                    Spn1MLzb3,Spn2MLzb3,Spn3MLzb3,Spn4MLzb3,Spn5MLzb3,Spn6MLzb3,Spn7MLzb3,Spn8MLzb3,Spn9MLzb3  &
                                /), (/9, 3/) )

INTEGER,  PARAMETER          :: SpnTDxb(9,3) = RESHAPE( (/ &
                                    Spn1TDxb1,Spn2TDxb1,Spn3TDxb1,Spn4TDxb1,Spn5TDxb1,Spn6TDxb1,Spn7TDxb1,Spn8TDxb1,Spn9TDxb1, &
                                    Spn1TDxb2,Spn2TDxb2,Spn3TDxb2,Spn4TDxb2,Spn5TDxb2,Spn6TDxb2,Spn7TDxb2,Spn8TDxb2,Spn9TDxb2, &
                                    Spn1TDxb3,Spn2TDxb3,Spn3TDxb3,Spn4TDxb3,Spn5TDxb3,Spn6TDxb3,Spn7TDxb3,Spn8TDxb3,Spn9TDxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnTDyb(9,3) = RESHAPE( (/ &
                                    Spn1TDyb1,Spn2TDyb1,Spn3TDyb1,Spn4TDyb1,Spn5TDyb1,Spn6TDyb1,Spn7TDyb1,Spn8TDyb1,Spn9TDyb1, &
                                    Spn1TDyb2,Spn2TDyb2,Spn3TDyb2,Spn4TDyb2,Spn5TDyb2,Spn6TDyb2,Spn7TDyb2,Spn8TDyb2,Spn9TDyb2, &
                                    Spn1TDyb3,Spn2TDyb3,Spn3TDyb3,Spn4TDyb3,Spn5TDyb3,Spn6TDyb3,Spn7TDyb3,Spn8TDyb3,Spn9TDyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnTDzb(9,3) = RESHAPE( (/ &
                                    Spn1TDzb1,Spn2TDzb1,Spn3TDzb1,Spn4TDzb1,Spn5TDzb1,Spn6TDzb1,Spn7TDzb1,Spn8TDzb1,Spn9TDzb1, &
                                    Spn1TDzb2,Spn2TDzb2,Spn3TDzb2,Spn4TDzb2,Spn5TDzb2,Spn6TDzb2,Spn7TDzb2,Spn8TDzb2,Spn9TDzb2, &
                                    Spn1TDzb3,Spn2TDzb3,Spn3TDzb3,Spn4TDzb3,Spn5TDzb3,Spn6TDzb3,Spn7TDzb3,Spn8TDzb3,Spn9TDzb3  &
                                /), (/9, 3/) )

INTEGER,  PARAMETER          :: SpnRDxb(9,3) = RESHAPE( (/ &
                                    Spn1RDxb1,Spn2RDxb1,Spn3RDxb1,Spn4RDxb1,Spn5RDxb1,Spn6RDxb1,Spn7RDxb1,Spn8RDxb1,Spn9RDxb1, &
                                    Spn1RDxb2,Spn2RDxb2,Spn3RDxb2,Spn4RDxb2,Spn5RDxb2,Spn6RDxb2,Spn7RDxb2,Spn8RDxb2,Spn9RDxb2, &
                                    Spn1RDxb3,Spn2RDxb3,Spn3RDxb3,Spn4RDxb3,Spn5RDxb3,Spn6RDxb3,Spn7RDxb3,Spn8RDxb3,Spn9RDxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnRDyb(9,3) = RESHAPE( (/ &
                                    Spn1RDyb1,Spn2RDyb1,Spn3RDyb1,Spn4RDyb1,Spn5RDyb1,Spn6RDyb1,Spn7RDyb1,Spn8RDyb1,Spn9RDyb1, &
                                    Spn1RDyb2,Spn2RDyb2,Spn3RDyb2,Spn4RDyb2,Spn5RDyb2,Spn6RDyb2,Spn7RDyb2,Spn8RDyb2,Spn9RDyb2, &
                                    Spn1RDyb3,Spn2RDyb3,Spn3RDyb3,Spn4RDyb3,Spn5RDyb3,Spn6RDyb3,Spn7RDyb3,Spn8RDyb3,Spn9RDyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnRDzb(9,3) = RESHAPE( (/ &
                                    Spn1RDzb1,Spn2RDzb1,Spn3RDzb1,Spn4RDzb1,Spn5RDzb1,Spn6RDzb1,Spn7RDzb1,Spn8RDzb1,Spn9RDzb1, &
                                    Spn1RDzb2,Spn2RDzb2,Spn3RDzb2,Spn4RDzb2,Spn5RDzb2,Spn6RDzb2,Spn7RDzb2,Spn8RDzb2,Spn9RDzb2, &
                                    Spn1RDzb3,Spn2RDzb3,Spn3RDzb3,Spn4RDzb3,Spn5RDzb3,Spn6RDzb3,Spn7RDzb3,Spn8RDzb3,Spn9RDzb3  &
                                /), (/9, 3/) )


INTEGER,  PARAMETER          :: TwHtALxt(9) = (/ &
                                    TwHt1ALxt,TwHt2ALxt,TwHt3ALxt,TwHt4ALxt,TwHt5ALxt,TwHt6ALxt,TwHt7ALxt,TwHt8ALxt,TwHt9ALxt /)
INTEGER,  PARAMETER          :: TwHtALyt(9) = (/ &
                                    TwHt1ALyt,TwHt2ALyt,TwHt3ALyt,TwHt4ALyt,TwHt5ALyt,TwHt6ALyt,TwHt7ALyt,TwHt8ALyt,TwHt9ALyt /)
INTEGER,  PARAMETER          :: TwHtALzt(9) = (/ &
                                    TwHt1ALzt,TwHt2ALzt,TwHt3ALzt,TwHt4ALzt,TwHt5ALzt,TwHt6ALzt,TwHt7ALzt,TwHt8ALzt,TwHt9ALzt /)

INTEGER,  PARAMETER          :: TwHtMLxt(9) = (/ &
                                    TwHt1MLxt,TwHt2MLxt,TwHt3MLxt,TwHt4MLxt,TwHt5MLxt,TwHt6MLxt,TwHt7MLxt,TwHt8MLxt,TwHt9MLxt /)
INTEGER,  PARAMETER          :: TwHtMLyt(9) = (/ &
                                    TwHt1MLyt,TwHt2MLyt,TwHt3MLyt,TwHt4MLyt,TwHt5MLyt,TwHt6MLyt,TwHt7MLyt,TwHt8MLyt,TwHt9MLyt /)
INTEGER,  PARAMETER          :: TwHtMLzt(9) = (/ &
                                    TwHt1MLzt,TwHt2MLzt,TwHt3MLzt,TwHt4MLzt,TwHt5MLzt,TwHt6MLzt,TwHt7MLzt,TwHt8MLzt,TwHt9MLzt /)

INTEGER,  PARAMETER          :: TwHtFLxt(9) = (/ &
                                    TwHt1FLxt,TwHt2FLxt,TwHt3FLxt,TwHt4FLxt,TwHt5FLxt,TwHt6FLxt,TwHt7FLxt,TwHt8FLxt,TwHt9FLxt /)
INTEGER,  PARAMETER          :: TwHtFLyt(9) = (/ &
                                    TwHt1FLyt,TwHt2FLyt,TwHt3FLyt,TwHt4FLyt,TwHt5FLyt,TwHt6FLyt,TwHt7FLyt,TwHt8FLyt,TwHt9FLyt /)
INTEGER,  PARAMETER          :: TwHtFLzt(9) = (/ &
                                    TwHt1FLzt,TwHt2FLzt,TwHt3FLzt,TwHt4FLzt,TwHt5FLzt,TwHt6FLzt,TwHt7FLzt,TwHt8FLzt,TwHt9FLzt /)

INTEGER,  PARAMETER          :: TwHtTDxt(9) = (/ &
                                    TwHt1TDxt,TwHt2TDxt,TwHt3TDxt,TwHt4TDxt,TwHt5TDxt,TwHt6TDxt,TwHt7TDxt,TwHt8TDxt,TwHt9TDxt /)
INTEGER,  PARAMETER          :: TwHtTDyt(9) = (/ &
                                    TwHt1TDyt,TwHt2TDyt,TwHt3TDyt,TwHt4TDyt,TwHt5TDyt,TwHt6TDyt,TwHt7TDyt,TwHt8TDyt,TwHt9TDyt /)
INTEGER,  PARAMETER          :: TwHtTDzt(9) = (/ &
                                    TwHt1TDzt,TwHt2TDzt,TwHt3TDzt,TwHt4TDzt,TwHt5TDzt,TwHt6TDzt,TwHt7TDzt,TwHt8TDzt,TwHt9TDzt /)

INTEGER,  PARAMETER          :: TwHtRDxt(9) = (/ &
                                    TwHt1RDxt,TwHt2RDxt,TwHt3RDxt,TwHt4RDxt,TwHt5RDxt,TwHt6RDxt,TwHt7RDxt,TwHt8RDxt,TwHt9RDxt /)
INTEGER,  PARAMETER          :: TwHtRDyt(9) = (/ &
                                    TwHt1RDyt,TwHt2RDyt,TwHt3RDyt,TwHt4RDyt,TwHt5RDyt,TwHt6RDyt,TwHt7RDyt,TwHt8RDyt,TwHt9RDyt /)
INTEGER,  PARAMETER          :: TwHtRDzt(9) = (/ &
                                    TwHt1RDzt,TwHt2RDzt,TwHt3RDzt,TwHt4RDzt,TwHt5RDzt,TwHt6RDzt,TwHt7RDzt,TwHt8RDzt,TwHt9RDzt /)

INTEGER,  PARAMETER          :: TwHtTPxi(9) = (/ &
                                    TwHt1TPxi,TwHt2TPxi,TwHt3TPxi,TwHt4TPxi,TwHt5TPxi,TwHt6TPxi,TwHt7TPxi,TwHt8TPxi,TwHt9TPxi /)
INTEGER,  PARAMETER          :: TwHtTPyi(9) = (/ &
                                    TwHt1TPyi,TwHt2TPyi,TwHt3TPyi,TwHt4TPyi,TwHt5TPyi,TwHt6TPyi,TwHt7TPyi,TwHt8TPyi,TwHt9TPyi /)
INTEGER,  PARAMETER          :: TwHtTPzi(9) = (/ &
                                    TwHt1TPzi,TwHt2TPzi,TwHt3TPzi,TwHt4TPzi,TwHt5TPzi,TwHt6TPzi,TwHt7TPzi,TwHt8TPzi,TwHt9TPzi /)

INTEGER,  PARAMETER          :: TwHtRPxi(9) = (/ &
                                    TwHt1RPxi,TwHt2RPxi,TwHt3RPxi,TwHt4RPxi,TwHt5RPxi,TwHt6RPxi,TwHt7RPxi,TwHt8RPxi,TwHt9RPxi /)
INTEGER,  PARAMETER          :: TwHtRPyi(9) = (/ &
                                    TwHt1RPyi,TwHt2RPyi,TwHt3RPyi,TwHt4RPyi,TwHt5RPyi,TwHt6RPyi,TwHt7RPyi,TwHt8RPyi,TwHt9RPyi /)
INTEGER,  PARAMETER          :: TwHtRPzi(9) = (/ &
                                    TwHt1RPzi,TwHt2RPzi,TwHt3RPzi,TwHt4RPzi,TwHt5RPzi,TwHt6RPzi,TwHt7RPzi,TwHt8RPzi,TwHt9RPzi /)

END MODULE ElastoDyn_Parameters
!**********************************************************************************************************************************
MODULE ElastoDyn

   USE NWTC_Library
   USE NWTC_LAPACK


   USE ElastoDyn_Parameters
   USE ElastoDyn_Types


   IMPLICIT NONE

   PRIVATE
   
!   TYPE(ProgDesc), PARAMETER  :: ED_Ver = ProgDesc( 'ElastoDyn', 'v1.01.02b-bjj', '03-Oct-2013' )



      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ED_Init                           ! Initialization routine
   PUBLIC :: ED_End                            ! Ending routine (includes clean up)

   PUBLIC :: ED_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   PUBLIC :: ED_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: ED_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: ED_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: ED_UpdateDiscState                ! Tight coupling routine for updating discrete states

   !PUBLIC :: ED_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                            !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   !PUBLIC :: ED_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                            !   (Xd), and constraint-state (Z) equations all with respect to the continuous
   !                                            !   states (x)
   !PUBLIC :: ED_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                            !   (Xd), and constraint-state (Z) equations all with respect to the discrete
   !                                            !   states (xd)
   !PUBLIC :: ED_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                            !   (Xd), and constraint-state (Z) equations all with respect to the constraint
   !                                            !   states (z)


CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

   TYPE(ED_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(ED_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(ED_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
   TYPE(ED_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(ED_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(ED_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(ED_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
   TYPE(ED_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                               !   only the output mesh is initialized)
   REAL(DbKi),                   INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                               !   (1) ED_UpdateStates() is called in loose coupling &
                                                               !   (2) ED_UpdateDiscState() is called in tight coupling.
                                                               !   Input is the suggested time from the glue code;
                                                               !   Output is the actual coupling interval that will be used
                                                               !   by the glue code.
   TYPE(ED_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables

   TYPE(ED_InputFile)                           :: InputFileData           ! Data stored in the module's input file
   INTEGER(IntKi)                               :: ErrStat2                ! temporary Error status of the operation
   INTEGER(IntKi)                               :: i, K                    ! loop counters
   LOGICAL, PARAMETER                           :: GetAdamsVals = .FALSE.  ! Determines if we should read Adams values and create (update) an Adams model
   CHARACTER(ErrMsgLen)                         :: ErrMsg2                 ! temporary Error message if ErrStat /= ErrID_None


      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   CALL DispNVD( ED_Ver )

      !............................................................................................
      ! Read the input file and validate the data
      !............................................................................................
   p%BD4Blades = .NOT. InitInp%CompElast           ! if we're not using ElastoDyn for the blades, use BeamDyn
   p%UseAD14   = LEN_TRIM(InitInp%ADInputFile) > 0 ! if we're using AD14, we need to use the AD14 input files

   p%RootName = TRIM(InitInp%RootName)//'.'//ED_Nickname ! all of the output file names from this module will contain '.ED' before the extension

   CALL ED_ReadInput( InitInp%InputFile, InitInp%ADInputFile, InputFileData, GetAdamsVals, p%BD4Blades, Interval, p%RootName, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   IF ( p%BD4Blades ) THEN
   
         ! Set DOFs to FALSE for whatever values you don't want on for BeamDyn
      InputFileData%FlapDOF1 = .FALSE.
      InputFileData%FlapDOF2 = .FALSE.
      InputFileData%EdgeDOF  = .FALSE.
      
         ! Set other values not used for BeamDyn      
      InputFileData%OoPDefl  = 0.0_ReKi
      InputFileData%IPDefl   = 0.0_ReKi
      InputFileData%TipMass  = 0.0_ReKi            
      InputFileData%TipRad   = 0.0_ReKi
      InputFileData%NBlGages = 0
      InputFileData%BldGagNd = 0
       
   END IF


   CALL ED_ValidateInput( InputFileData, p%BD4Blades, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      !............................................................................................
      ! Define parameters here:
      !............................................................................................
   CALL ED_SetParameters( InputFileData, p, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      !............................................................................................
      ! Define initial system states here:
      !............................................................................................
   xd%DummyDiscState          = 0.                                             ! we don't have discrete states
   z%DummyConstrState         = 0.                                             ! we don't have constraint states

      ! initialize the continuous states:
   CALL Init_ContStates( x, p, InputFileData, OtherState, ErrStat2, ErrMsg2 )     ! initialize the continuous states
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   
   
      ! Initialize other states:
   CALL Init_OtherStates( OtherState, p, x, InputFileData, ErrStat2, ErrMsg2 )    ! initialize the other states (must do this after ED_SetParameters)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      

      !............................................................................................
      ! Define initial guess for the system inputs here:
      !............................................................................................

         ! allocate all the arrays that store data in the input type and initialize the values:
   CALL Init_u( u, p, x, InputFileData, OtherState, ErrStat2, ErrMsg2 )      
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN                  
   

      !............................................................................................
      ! Define system output initializations (set up meshes) here:
      !............................................................................................

   CALL ED_AllocOutput(p, OtherState, u, y, ErrStat2, ErrMsg2) ! u is sent so we can create sibling meshes
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
   IF (p%BD4Blades) THEN
         ! we need the initial outputs to send BeamDyn for its states. I should probably put them in
         ! the InitOutput type, but I'm going to use the ED_Outputs instead
      ! 
      !   ! set the coordinate system variables:
      !CALL SetCoordSy( t, OtherState%CoordSys, OtherState%RtHS, u%BlPitchCom, p, x, ErrStat, ErrMsg )
      !   IF (ErrStat >= AbortErrLev) RETURN
      !
      !CALL CalculatePositions(        p, x, OtherState%CoordSys,    OtherState%RtHS ) ! calculate positions
      !CALL CalculateAngularPosVelPAcc(p, x, OtherState%CoordSys,    OtherState%RtHS ) ! calculate angular positions, velocities, and partial accelerations, including partial angular quantities
      !CALL CalculateLinearVelPAcc(    p, x, OtherState%CoordSys,    OtherState%RtHS ) ! calculate linear velocities and partial accelerations

      
      
      CALL ED_CalcOutput( 0.0_DbKi, u, p, x, xd, z, OtherState, y, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
   END IF
      
      
      !............................................................................................
      ! Define initialization-routine output here:
      !............................................................................................
   CALL AllocAry( InitOut%WriteOutputHdr, p%NumOuts, 'WriteOutputHdr', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   CALL AllocAry( InitOut%WriteOutputUnt, p%NumOuts, 'WriteOutputUnt', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   do i=1,p%NumOuts
      InitOut%WriteOutputHdr(i) = p%OutParam(i)%Name
      InitOut%WriteOutputUnt(i) = p%OutParam(i)%Units
   end do
      
   InitOut%Ver         = ED_Ver
   InitOut%NumBl       = p%NumBl
   InitOut%Gravity     = p%Gravity
   InitOut%BladeLength = p%TipRad - p%HubRad
   InitOut%PlatformPos = x%QT(1:6)
   InitOut%HubHt       = p%HubHt
   InitOut%TwrBasePos  = y%TowerLn2Mesh%Position(:,p%TwrNodes + 2)

   CALL AllocAry(InitOut%BlPitch, p%NumBl, 'BlPitch', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   InitOut%BlPitch = InputFileData%BlPitch(1:p%NumBl)

      !............................................................................................
      ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
      !   this module must be called here:
      !............................................................................................

   Interval = p%DT


       ! Print the summary file if requested:
   IF (InputFileData%SumPrint) THEN
      CALL ED_PrintSum( p, OtherState, GetAdamsVals, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
   END IF
       
       ! Destroy the InputFileData structure (deallocate arrays)

   CALL ED_DestroyInputFile(InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
         

CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN         
         
         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_Init:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL ED_DestroyInputFile(InputFileData, ErrStat3, ErrMsg3 )
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE ED_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(ED_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(ED_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(ED_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(ED_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(ED_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(ED_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ED_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:



         ! Destroy the input data:

      CALL ED_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      CALL ED_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      CALL ED_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL ED_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL ED_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL ED_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )


         ! Destroy the output data:

      CALL ED_DestroyOutput( y, ErrStat, ErrMsg )




END SUBROUTINE ED_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time t; Continuous and discrete states are updated for t + Interval
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t          ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
      TYPE(ED_InputType),                 INTENT(INOUT) :: u(:)       ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
      TYPE(ED_ParameterType),             INTENT(IN   ) :: p          ! Parameters
      TYPE(ED_ContinuousStateType),       INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      TYPE(ED_DiscreteStateType),         INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      TYPE(ED_ConstraintStateType),       INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      TYPE(ED_OtherStateType),            INTENT(INOUT) :: OtherState ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      
         ! local variables

      !TYPE(ED_InputType)            :: u_interp  ! input interpolated from given u at utimes
      !TYPE(ED_ContinuousStateType)  :: xdot      ! continuous state time derivative


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""            
      

      SELECT CASE ( p%method )
         
      CASE (Method_RK4)
      
         CALL ED_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
         
      CASE (Method_AB4)
      
         CALL ED_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
      
      CASE (Method_ABM4)
      
         CALL ED_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
         
      CASE DEFAULT  !bjj: we already checked this at initialization, but for completeness:
         
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error in ElastoDyn_UpdateStates: p%method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
         RETURN
         
      END SELECT
      
         ! Make sure the rotor azimuth is not greater or equal to 360 degrees: 
         
      ! bjj: per jmj, the subtraction of TwoPi here is so that we don't run into numerical issues with large GeAz (in large simulations)
      !   this subtraction is okay because we use x%QT(DOF_GeAz) only in equations with SIN() and/or COS() so it doesn't matter
      !   if there is a discontinunity in the channel.
      ! bjj: why don't we just do a modulo on x%QT(DOF_GeAz) instead of using x%QT(DOF_DrTr) with it?   
      
      IF ( ( x%QT(DOF_GeAz) + x%QT(DOF_DrTr) ) >= TwoPi_D )  x%QT(DOF_GeAz) = x%QT(DOF_GeAz) - TwoPi_D
            
      
END SUBROUTINE ED_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! Routine for computing outputs, used in both loose and tight coupling.
! This SUBROUTINE is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
! NOTE: the descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
! for a complete description of each output parameter.
! NOTE: no matter how many channels are selected for output, all of the outputs are calcalated
! All of the calculated output channels are placed into the OtherState%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(ED_InputType),           INTENT(IN   )  :: u           ! Inputs at Time t
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(ED_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(ED_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(ED_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                               !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables:

   REAL(ReKi)                   :: AngAccEB  (3)                                   ! Angular acceleration of the base plate                                                (body B) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccEN  (3)                                   ! Angular acceleration of the nacelle                                                   (body N) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccEH  (3)                                   ! Angular acceleration of the hub                                                   (body N) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccER  (3)                                   ! Angular acceleration of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccEX  (3)                                   ! Angular acceleration of the platform                                                  (body X) in the inertia frame (body E for earth).
!   REAL(ReKi)                   :: ComDenom                                        ! Common denominator used in several expressions.
!   REAL(ReKi)                   :: CThrstys                                        ! Estimate of the ys-location of the center of thrust.
!   REAL(ReKi)                   :: CThrstzs                                        ! Estimate of the zs-location of the center of thrust.
   REAL(R8Ki)                   :: FrcMGagB  (3)                                   ! Total force at the blade element   (body M) / blade strain gage location            (point S) due to the blade above the strain gage.
   REAL(ReKi)                   :: FrcFGagT  (3)                                   ! Total force at the tower element   (body F) / tower strain gage location            (point T) due to the nacelle and rotor and tower above the strain gage.
   REAL(ReKi)                   :: FrcONcRt  (3)                                   ! Total force at the yaw bearing (point O  ) due to the nacelle, generator, and rotor
   REAL(ReKi)                   :: FrcPRot   (3)                                   ! Total force at the teeter pin  (point P  ) due to the rotor
   REAL(ReKi)                   :: FrcT0Trb  (3)                                   ! Total force at the base of flexible portion of the tower (point T(0)) due to the entire wind turbine
   REAL(ReKi)                   :: FZHydro   (3)                                   ! Total platform hydrodynamic force at the platform reference (point Z)
!   REAL(ReKi)                   :: HHWndVec  (3)                                   ! Hub-height wind vector in the AeroDyn coordinate system
   REAL(ReKi)                   :: LinAccEIMU(3)                                   ! Total linear acceleration of the nacelle IMU (point IMU) in the inertia frame (body E for earth)
   REAL(ReKi)                   :: LinAccEO  (3)                                   ! Total linear acceleration of the base plate (point O) in the inertia frame (body E for earth)
   REAL(ReKi)                   :: LinAccEZ  (3)                                   ! Total linear acceleration of the platform refernce (point Z) in the inertia frame (body E for earth)
   REAL(ReKi)                   :: MomBNcRt  (3)                                   ! Total moment at the base plate      (body B) / yaw bearing                           (point O) due to the nacelle, generator, and rotor.
   REAL(ReKi)                   :: MomFGagT  (3)                                   ! Total moment at the tower element   (body F) / tower strain gage location            (point T) due to the nacelle and rotor and tower above the strain gage.
   REAL(ReKi)                   :: MomLPRot  (3)                                   ! Total moment at the low-speed shaft (body L) / teeter pin                            (point P) due to the rotor.
   REAL(ReKi)                   :: MomMGagB  (3)                                   ! Total moment at the blade element   (body M) / blade strain gage location            (point S) due to the blade above the strain gage.
   REAL(ReKi)                   :: MomNGnRt  (3)                                   ! Total moment at the nacelle         (body N) / specified point on rotor-furl axis    (point V) due to the structure that furls with the rotor, generator, and rotor.
   REAL(ReKi)                   :: MomNTail  (3)                                   ! Total moment at the nacelle         (body N) / specified point on  tail-furl axis    (point W) due to the tail.
   REAL(ReKi)                   :: MomX0Trb  (3)                                   ! Total moment at the tower base      (body X) / base of flexible portion of the tower (point T(0)) due to the entire wind turbine.
   REAL(ReKi)                   :: MXHydro   (3)                                   ! Total platform hydrodynamic moment acting at the platform (body X) / platform reference (point Z).
   REAL(R8Ki)                   :: rOPO      (3)                                   ! Position vector from the undeflected tower top (point O prime) to the deflected tower top (point O).
   REAL(R8Ki)                   :: rOSTip    (3)                                   ! Position vector from the deflected tower top (point O) to the deflected blade tip (point S tip).
   REAL(R8Ki)                   :: rOSTipxn                                        ! Component of rOSTip directed along the xn-axis.
   REAL(R8Ki)                   :: rOSTipyn                                        ! Component of rOSTip directed along the yn-axis.
   REAL(R8Ki)                   :: rOSTipzn                                        ! Component of rOSTip directed along the zn-axis.
   REAL(R8Ki)                   :: rTPT      (3)                                   ! Position vector from the undeflected tower node (point T prime) to the deflected node (point T)
   REAL(R8Ki)                   :: rSPS      (3)                                   ! Position vector from the undeflected blade node (point S prime) to the deflected node (point S)
   REAL(R8Ki)                   :: rSTipPSTip(3)                                   ! Position vector from the undeflected blade tip (point S tip prime) to the deflected blade tip (point S tip).
   REAL(R8Ki)                   :: TmpVec    (3)                                   ! A temporary vector used in various computations.
   REAL(R8Ki)                   :: TmpVec2   (3)                                   ! A temporary vector.

   INTEGER, PARAMETER           :: NDims = 3
   REAL(ReKi)                   :: LinAccES (NDims,0:p%TipNode,p%NumBl)            ! Total linear acceleration of a point on a   blade (point S) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: LinAccET (NDims,0:p%TwrNodes)                   ! Total linear acceleration of a point on the tower (point T) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccEF (NDims,0:p%TwrNodes)                   ! Total angular acceleration of tower element J (body F) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: FrcS0B   (NDims,p%NumBl)                        ! Total force at the blade root (point S(0)) due to the blade.
   REAL(ReKi)                   :: FTTower  (NDims,p%TwrNodes)                     ! Total hydrodynamic + aerodynamic force per unit length acting on the tower at point T.
   REAL(ReKi)                   :: MFHydro  (NDims,p%TwrNodes)                     ! Total hydrodynamic + aerodynamic moment per unit length acting on a tower element (body F) at point T.
   REAL(ReKi)                   :: MomH0B   (NDims,p%NumBl)                        ! Total moment at the hub (body H) / blade root (point S(0)) due to the blade.


   INTEGER(IntKi)               :: I                                               ! Generic index
   INTEGER(IntKi)               :: J, J2                                           ! Loops through nodes / elements
   INTEGER(IntKi)               :: K                                               ! Loops through blades
   INTEGER(IntKi)               :: NodeNum                                         ! Mesh node number for given blade/node
   
   INTEGER(IntKi)               :: ErrStat2                                        ! Temporary Error code
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                         ! Temporary error message
   

   LOGICAL, PARAMETER           :: UpdateValues  = .TRUE.                          ! determines if the OtherState values need to be updated
   TYPE(ED_ContinuousStateType) :: dxdt                                            ! Continuous state derivs at t

         ! Initialize some output values
      ErrStat = ErrID_None
      ErrMsg  = ""

      
      ! SEE IF THESE NEED TO BE CALLED (i.e., if UpdateStates was called, these values are already calculated)
   IF ( UpdateValues ) THEN    
         ! Update the OtherState data by calculating the derivative...
      !OtherState%HSSBrTrqC = SIGN( u%HSSBrTrqC, x%QDT(DOF_GeAz) )
      CALL ED_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg ) ! sets OtherState%QD2T = dxdt%QDT
      CALL ED_DestroyContState( dxdt, ErrStat2, ErrMsg2 )  
      IF (ErrStat >= AbortErrLev) RETURN
   END IF      


   
      !..............................
      ! Outputs for HydroDyn, UsrPtfm and UsrTwr
      !..............................
            
      !u_UsrPtfm%X  = x_ED%QT(1:6)
      !u_UsrPtfm%XD = x_ED%QDT(1:6)

      !u_UsrTwr%X(:,J) = (/ OtherSt_ED%RtHS%rT(1,J),       -OtherSt_ED%RtHS%rT(3,J),       OtherSt_ED%RtHS%rT( 2,J)+ p%PtfmRefzt,&
      !                     OtherSt_ED%RtHS%AngPosEF(1,J), -OtherSt_ED%RtHS%AngPosEF(3,J), OtherSt_ED%RtHS%AngPosEF(2,J)         /) 
      !u_UsrTwr%XD(:,J) = (/ OtherSt_ED%RtHS%LinVelET(1,J), -OtherSt_ED%RtHS%LinVelET(3,J), OtherSt_ED%RtHS%LinVelET(2,J),&
      !                      OtherSt_ED%RtHS%AngVelEF(1,J), -OtherSt_ED%RtHS%AngVelEF(3,J), OtherSt_ED%RtHS%AngVelEF(2,J) /)                     
   
   
   
   
   ! Array OtherState%AllOuts() is initialized to 0.0 in initialization, so we are not going to reinitialize it here.

   !...............................................................................................................................
   ! Calculate all of the total forces and moments using all of the partial forces and moments calculated in RtHS().  Also,
   !   calculate all of the total angular and linear accelerations using all of the partial accelerations calculated in RtHS().
   !   To do this, first initialize the variables using the portions not associated with the accelerations.  Then add the portions
   !   associated with the accelerations one by one:
   !...............................................................................................................................

   AngAccEB   = OtherState%RtHS%AngAccEBt
   AngAccEH   = OtherState%RtHS%AngAccEHt
   AngAccEN   = OtherState%RtHS%AngAccENt
   AngAccER   = OtherState%RtHS%AngAccERt
   AngAccEX   = OtherState%RtHS%AngAccEXt
   LinAccEIMU = OtherState%RtHS%LinAccEIMUt
   LinAccEO   = OtherState%RtHS%LinAccEOt
   LinAccEZ   = OtherState%RtHS%LinAccEZt
   FrcONcRt   = OtherState%RtHS%FrcONcRtt
   FrcPRot    = OtherState%RtHS%FrcPRott
   FrcT0Trb   = OtherState%RtHS%FrcT0Trbt

   ! was FZHydro    = OtherState%RtHS%FZHydrot
   FZHydro    = u%PlatformPtMesh%Force(DOF_Sg,1)*OtherState%RtHS%PLinVelEZ(DOF_Sg,0,:) &
              + u%PlatformPtMesh%Force(DOF_Sw,1)*OtherState%RtHS%PLinVelEZ(DOF_Sw,0,:) &
              + u%PlatformPtMesh%Force(DOF_Hv,1)*OtherState%RtHS%PLinVelEZ(DOF_Hv,0,:)

   MomBNcRt   = OtherState%RtHS%MomBNcRtt
   MomLPRot   = OtherState%RtHS%MomLPRott
   MomNGnRt   = OtherState%RtHS%MomNGnRtt
   MomNTail   = OtherState%RtHS%MomNTailt
   MomX0Trb   = OtherState%RtHS%MomX0Trbt

   ! was MXHydro = OtherState%RtHS%MXHydrot
   MXHydro    = u%PlatformPtMesh%Moment(DOF_R-3,1)*OtherState%RtHS%PAngVelEX(DOF_R ,0,:) &
              + u%PlatformPtMesh%Moment(DOF_P-3,1)*OtherState%RtHS%PAngVelEX(DOF_P ,0,:) &
              + u%PlatformPtMesh%Moment(DOF_Y-3,1)*OtherState%RtHS%PAngVelEX(DOF_Y ,0,:)

   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
      AngAccEB   = AngAccEB   + OtherState%RtHS%PAngVelEB  (p%DOFs%SrtPS(I),0,:)*OtherState%QD2T(p%DOFs%SrtPS(I))
      AngAccEH   = AngAccEH   + OtherState%RtHS%PAngVelEH  (p%DOFs%SrtPS(I),0,:)*OtherState%QD2T(p%DOFs%SrtPS(I))      
      AngAccEN   = AngAccEN   + OtherState%RtHS%PAngVelEN  (p%DOFs%SrtPS(I),0,:)*OtherState%QD2T(p%DOFs%SrtPS(I))      
      AngAccER   = AngAccER   + OtherState%RtHS%PAngVelER  (p%DOFs%SrtPS(I),0,:)*OtherState%QD2T(p%DOFs%SrtPS(I))
      LinAccEIMU = LinAccEIMU + OtherState%RtHS%PLinVelEIMU(p%DOFs%SrtPS(I),0,:)*OtherState%QD2T(p%DOFs%SrtPS(I))
      LinAccEO   = LinAccEO   + OtherState%RtHS%PLinVelEO  (p%DOFs%SrtPS(I),0,:)*OtherState%QD2T(p%DOFs%SrtPS(I))
      FrcONcRt   = FrcONcRt   + OtherState%RtHS%PFrcONcRt  (:,p%DOFs%SrtPS(I)  )*OtherState%QD2T(p%DOFs%SrtPS(I))
      FrcPRot    = FrcPRot    + OtherState%RtHS%PFrcPRot   (:,p%DOFs%SrtPS(I)  )*OtherState%QD2T(p%DOFs%SrtPS(I))
      FrcT0Trb   = FrcT0Trb   + OtherState%RtHS%PFrcT0Trb  (:,p%DOFs%SrtPS(I)  )*OtherState%QD2T(p%DOFs%SrtPS(I))
      MomBNcRt   = MomBNcRt   + OtherState%RtHS%PMomBNcRt  (:,p%DOFs%SrtPS(I)  )*OtherState%QD2T(p%DOFs%SrtPS(I))
      MomLPRot   = MomLPRot   + OtherState%RtHS%PMomLPRot  (:,p%DOFs%SrtPS(I)  )*OtherState%QD2T(p%DOFs%SrtPS(I))
      MomNGnRt   = MomNGnRt   + OtherState%RtHS%PMomNGnRt  (:,p%DOFs%SrtPS(I)  )*OtherState%QD2T(p%DOFs%SrtPS(I))
      MomNTail   = MomNTail   + OtherState%RtHS%PMomNTail  (:,p%DOFs%SrtPS(I)  )*OtherState%QD2T(p%DOFs%SrtPS(I))
      MomX0Trb   = MomX0Trb   + OtherState%RtHS%PMomX0Trb  (:,p%DOFs%SrtPS(I)  )*OtherState%QD2T(p%DOFs%SrtPS(I))
   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPYE     ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)
      AngAccEX   = AngAccEX   + OtherState%RtHS%PAngVelEX  (p%DOFs%PYE  (I),0,:)*OtherState%QD2T(p%DOFs%PYE  (I))
      LinAccEZ   = LinAccEZ   + OtherState%RtHS%PLinVelEZ  (p%DOFs%PYE  (I),0,:)*OtherState%QD2T(p%DOFs%PYE  (I))


      FZHydro    = FZHydro    + (- u%PtfmAddedMass(DOF_Sg,p%DOFs%PYE(I))*OtherState%RtHS%PLinVelEZ(DOF_Sg,0,:)   &  ! was  FZHydro = FZHydro + OtherState%RtHS%PFZHydro(p%DOFs%PYE(I),:)*OtherState%QD2T(p%DOFs%PYE  (I))
                                 - u%PtfmAddedMass(DOF_Sw,p%DOFs%PYE(I))*OtherState%RtHS%PLinVelEZ(DOF_Sw,0,:)   &
                                 - u%PtfmAddedMass(DOF_Hv,p%DOFs%PYE(I))*OtherState%RtHS%PLinVelEZ(DOF_Hv,0,:) ) &
                                *OtherState%QD2T(p%DOFs%PYE  (I))
      ! was MXHydro = MXHydro    +  OtherState%RtHS%PMXHydro   (p%DOFs%PYE  (I),  :)*OtherState%QD2T(p%DOFs%PYE  (I))
      MXHydro    = MXHydro    +  (- u%PtfmAddedMass(DOF_R ,p%DOFs%PYE(I))*OtherState%RtHS%PAngVelEX(DOF_R ,0,:)   &
                                  - u%PtfmAddedMass(DOF_P ,p%DOFs%PYE(I))*OtherState%RtHS%PAngVelEX(DOF_P ,0,:)   &
                                  - u%PtfmAddedMass(DOF_Y ,p%DOFs%PYE(I))*OtherState%RtHS%PAngVelEX(DOF_Y ,0,:) ) &
                                *OtherState%QD2T(p%DOFs%PYE  (I))

   ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)



   DO K = 1,p%NumBl ! Loop through all blades

      FrcS0B  (:          ,K) = OtherState%RtHS%FrcS0Bt  (:,K          )
      MomH0B  (:          ,K) = OtherState%RtHS%MomH0Bt  (:,K          )

      DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K
         FrcS0B  (:          ,K) = FrcS0B  (:          ,K) + OtherState%RtHS%PFrcS0B  (:,K,          p%DOFs%PSE(K,I)  )*OtherState%QD2T(p%DOFs%PSE(K,I))
         MomH0B  (:          ,K) = MomH0B  (:          ,K) + OtherState%RtHS%PMomH0B  (:,K,          p%DOFs%PSE(K,I)  )*OtherState%QD2T(p%DOFs%PSE(K,I))
      ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

      DO J = 0,p%TipNode ! Loop through the blade nodes / elements

         LinAccES(:,J,K) = OtherState%RtHS%LinAccESt(:,K,J)

         DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K
            LinAccES(:,J,K) = LinAccES(:,J,K) + OtherState%RtHS%PLinVelES(K,J,p%DOFs%PSE(K,I),0,:)*OtherState%QD2T(p%DOFs%PSE(K,I))
         ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

      ENDDO             ! J - Blade nodes / elements

   ENDDO          ! K - All blades

   DO J = 0,p%TwrNodes  ! Loop through the tower nodes / elements, starting at the tower base (0)

      LinAccET(:,J) = OtherState%RtHS%LinAccETt(:,J)
      AngAccEF(:,J) = OtherState%RtHS%AngAccEFt(:,J)

      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)
         LinAccET(:,J) = LinAccET(:,J) + OtherState%RtHS%PLinVelET(J,p%DOFs%PTE(I),0,:)*OtherState%QD2T(p%DOFs%PTE(I))
         AngAccEF(:,J) = AngAccEF(:,J) + OtherState%RtHS%PAngVelEF(J,p%DOFs%PTE(I),0,:)*OtherState%QD2T(p%DOFs%PTE(I))
      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

   ENDDO ! J - Tower nodes / elements


   DO J = 1,p%TwrNodes  ! Loop through the tower nodes / elements

      FTTower (:,J) = OtherState%RtHS%FTHydrot (:,J)
      MFHydro (:,J) = OtherState%RtHS%MFHydrot (:,J)

      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)
         FTTower (:,J) = FTTower (:,J) + OtherState%RtHS%PFTHydro (:,J,p%DOFs%PTE(I)  )*OtherState%QD2T(p%DOFs%PTE(I))
         MFHydro (:,J) = MFHydro (:,J) + OtherState%RtHS%PMFHydro (:,J,p%DOFs%PTE(I)  )*OtherState%QD2T(p%DOFs%PTE(I))
      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

   ENDDO ! J - Tower nodes / elements


      ! Convert the units of the forces and moments from N and N-m
      !    to kN and kN-m:

   FrcONcRt = 0.001*FrcONcRt
   FrcPRot  = 0.001*FrcPRot
   FrcT0Trb = 0.001*FrcT0Trb
   FZHydro  = 0.001*FZHydro
   MomBNcRt = 0.001*MomBNcRt
   MomLPRot = 0.001*MomLPRot
   MomNGnRt = 0.001*MomNGnRt
   MomNTail = 0.001*MomNTail
   MomX0Trb = 0.001*MomX0Trb
   MXHydro  = 0.001*MXHydro
   FrcS0B   = 0.001*FrcS0B
   MomH0B   = 0.001*MomH0B


   !...............................................................................................................................
   ! set the values in the AllOuts array:
   !...............................................................................................................................
      ! Define the output channel specifying the current simulation time:

   OtherState%AllOuts(  Time) = REAL( t, ReKi )


      ! Blade (1-3) Tip Motions:

   DO K = 1,p%NumBl
      rSTipPSTip = OtherState%RtHS%rS0S(:,K,p%TipNode) - p%BldFlexL*OtherState%CoordSys%j3(K,:)  ! Position vector from the undeflected blade tip (point S tip prime) to the deflected blade tip (point S tip) of blade 1.
      rOSTip     = OtherState%RtHS%rS  (:,K,p%TipNode) - OtherState%RtHS%rO                ! Position vector from the deflected tower top (point O) to the deflected blade tip (point S tip) of blade 1.
      rOSTipxn   =      DOT_PRODUCT( rOSTip, OtherState%CoordSys%d1 )                ! Component of rOSTip directed along the xn-axis.
      rOSTipyn   = -1.0*DOT_PRODUCT( rOSTip, OtherState%CoordSys%d3 )                ! Component of rOSTip directed along the yn-axis.
      rOSTipzn   =      DOT_PRODUCT( rOSTip, OtherState%CoordSys%d2 )                ! Component of rOSTip directed along the zn-axis.

      IF (.NOT. p%BD4Blades) THEN
         OtherState%AllOuts(  TipDxc(K) ) = DOT_PRODUCT(            rSTipPSTip, OtherState%CoordSys%i1(K,         :) )
         OtherState%AllOuts(  TipDyc(K) ) = DOT_PRODUCT(            rSTipPSTip, OtherState%CoordSys%i2(K,         :) )
         OtherState%AllOuts(  TipDzc(K) ) = DOT_PRODUCT(            rSTipPSTip, OtherState%CoordSys%i3(K,         :) )
         OtherState%AllOuts(  TipDxb(K) ) = DOT_PRODUCT(            rSTipPSTip, OtherState%CoordSys%j1(K,         :) )
         OtherState%AllOuts(  TipDyb(K) ) = DOT_PRODUCT(            rSTipPSTip, OtherState%CoordSys%j2(K,         :) )
      !JASON: USE TipNode HERE INSTEAD OF BldNodes IF YOU ALLOCATE AND DEFINE n1, n2, n3, m1, m2, AND m3 TO USE TipNode.  THIS WILL REQUIRE THAT THE AERODYNAMIC AND STRUCTURAL TWISTS, AeroTwst() AND ThetaS(), BE KNOWN AT THE TIP!!!
         OtherState%AllOuts( TipALxb(K) ) = DOT_PRODUCT( LinAccES(:,p%TipNode,K), OtherState%CoordSys%n1(K,p%BldNodes,:) )
         OtherState%AllOuts( TipALyb(K) ) = DOT_PRODUCT( LinAccES(:,p%TipNode,K), OtherState%CoordSys%n2(K,p%BldNodes,:) )
         OtherState%AllOuts( TipALzb(K) ) = DOT_PRODUCT( LinAccES(:,p%TipNode,K), OtherState%CoordSys%n3(K,p%BldNodes,:) )
         OtherState%AllOuts( TipRDxb(K) ) = DOT_PRODUCT( OtherState%RtHS%AngPosHM(:,K,p%TipNode), OtherState%CoordSys%j1(K,         :) )*R2D
         OtherState%AllOuts( TipRDyb(K) ) = DOT_PRODUCT( OtherState%RtHS%AngPosHM(:,K,p%TipNode), OtherState%CoordSys%j2(K,         :) )*R2D
         ! There is no sense computing AllOuts( TipRDzc(K) ) here since it is always zero for FAST simulation results.
         IF ( rOSTipzn > 0.0 )  THEN   ! Tip of blade K is above the yaw bearing.
            OtherState%AllOuts(TipClrnc(K) ) = SQRT( rOSTipxn*rOSTipxn + rOSTipyn*rOSTipyn + rOSTipzn*rOSTipzn ) ! Absolute distance from the tower top / yaw bearing to the tip of blade 1.
         ELSE                          ! Tip of blade K is below the yaw bearing.
            OtherState%AllOuts(TipClrnc(K) ) = SQRT( rOSTipxn*rOSTipxn + rOSTipyn*rOSTipyn                     ) ! Perpendicular distance from the yaw axis / tower centerline to the tip of blade 1.
         ENDIF
      END IF      

   END DO !K

      ! Blade (1-3) Local Span Motions:

   DO K = 1,p%NumBl
      DO I = 1, p%NBlGages

         OtherState%AllOuts( SpnALxb(I,K) ) = DOT_PRODUCT( LinAccES(:,p%BldGagNd(I),K), OtherState%CoordSys%n1(K,p%BldGagNd(I),:) )
         OtherState%AllOuts( SpnALyb(I,K) ) = DOT_PRODUCT( LinAccES(:,p%BldGagNd(I),K), OtherState%CoordSys%n2(K,p%BldGagNd(I),:) )
         OtherState%AllOuts( SpnALzb(I,K) ) = DOT_PRODUCT( LinAccES(:,p%BldGagNd(I),K), OtherState%CoordSys%n3(K,p%BldGagNd(I),:) )

         rSPS                      = OtherState%RtHS%rS0S(:,K,p%BldGagNd(I)) - p%RNodes(p%BldGagNd(I))*OtherState%CoordSys%j3(K,:)

         OtherState%AllOuts( SpnTDxb(I,K) ) = DOT_PRODUCT( rSPS, OtherState%CoordSys%j1(K,:) )
         OtherState%AllOuts( SpnTDyb(I,K) ) = DOT_PRODUCT( rSPS, OtherState%CoordSys%j2(K,:) )
         OtherState%AllOuts( SpnTDzb(I,K) ) = DOT_PRODUCT( rSPS, OtherState%CoordSys%j3(K,:) )

         OtherState%AllOuts( SpnRDxb(I,K) ) = DOT_PRODUCT( OtherState%RtHS%AngPosHM(:,K,p%BldGagNd(I)), OtherState%CoordSys%j1(K,:) )*R2D
         OtherState%AllOuts( SpnRDyb(I,K) ) = DOT_PRODUCT( OtherState%RtHS%AngPosHM(:,K,p%BldGagNd(I)), OtherState%CoordSys%j2(K,:) )*R2D
        !OtherState%AllOuts( SpnRDzb(I,K) ) = DOT_PRODUCT( OtherState%RtHS%AngPosHM(:,K,p%BldGagNd(I)), OtherState%CoordSys%j3(K,:) )*R2D           ! this is always zero for FAST

      END DO !I
   END DO !K



      ! Blade Pitch Motions:

   OtherState%AllOuts(PtchPMzc1) = u%BlPitchCom(1)*R2D
   OtherState%AllOuts(PtchPMzc2) = u%BlPitchCom(2)*R2D
   IF ( p%NumBl == 3 )  THEN ! 3-blader

      OtherState%AllOuts(PtchPMzc3) = u%BlPitchCom(3)*R2D

   ELSE  ! 2-blader


      ! Teeter Motions:

      OtherState%AllOuts(  TeetPya) =x%QT  (DOF_Teet)*R2D
      OtherState%AllOuts(  TeetVya) =x%QDT (DOF_Teet)*R2D
      OtherState%AllOuts(  TeetAya) =  OtherState%QD2T(DOF_Teet)*R2D

   ENDIF


      ! Shaft Motions:

   y%LSSTipPxa = x%QT (DOF_GeAz) + x%QT  (DOF_DrTr) + p%AzimB1Up + PiBy2 !bjj: this used IgnoreMod for linearization
   CALL Zero2TwoPi(y%LSSTipPxa)  ! Return value between 0 and 2pi (LSSTipPxa is used only in calculations of SIN and COS, so it's okay to take MOD/MODULO here; this wouldn't be oaky for linearization)
   OtherState%AllOuts(LSSTipPxa) = y%LSSTipPxa*R2D
   
   OtherState%AllOuts(LSSGagPxa) = MODULO( (      x%QT (DOF_GeAz)                            + p%AzimB1Up)*R2D  + 90.0_R8Ki, 360.0_R8Ki ) !bjj: this used IgnoreMod for linearization (Zero2TwoPi)
   OtherState%AllOuts(   LSSTipVxa) =      (     x%QDT (DOF_GeAz) +          x%QDT (DOF_DrTr) )*RPS2RPM
   OtherState%AllOuts(   LSSTipAxa) = ( OtherState%QD2T(DOF_GeAz) + OtherState%QD2T(DOF_DrTr) )*R2D
   OtherState%AllOuts(   LSSGagVxa) =            x%QDT (DOF_GeAz)                              *RPS2RPM
   OtherState%AllOuts(   LSSGagAxa) =   OtherState%QD2T(DOF_GeAz)                              *R2D
   OtherState%AllOuts(     HSShftV) = ABS(p%GBRatio)*OtherState%AllOuts(LSSGagVxa)
   OtherState%AllOuts(     HSShftA) = ABS(p%GBRatio)*OtherState%AllOuts(LSSGagAxa)

   !IF ( .NOT. EqualRealNos( OtherState%AllOuts(WindVxi), 0.0_ReKi ) )  THEN  ! .TRUE. if the denominator in the following equation is not zero.
   !   OtherState%AllOuts(TipSpdRat) =      ( x%QDT (DOF_GeAz) + x%QDT (DOF_DrTr) )*p%AvgNrmTpRd / OtherState%AllOuts(  WindVxi)
   !ELSE
   !   OtherState%AllOuts(TipSpdRat) = 0.0
   !ENDIF


      ! Nacelle IMU Motions:

   OtherState%AllOuts(NcIMUTVxs) =      DOT_PRODUCT( OtherState%RtHS%LinVelEIMU, OtherState%CoordSys%c1 )
   OtherState%AllOuts(NcIMUTVys) = -1.0*DOT_PRODUCT( OtherState%RtHS%LinVelEIMU, OtherState%CoordSys%c3 )
   OtherState%AllOuts(NcIMUTVzs) =      DOT_PRODUCT( OtherState%RtHS%LinVelEIMU, OtherState%CoordSys%c2 )
   OtherState%AllOuts(NcIMUTAxs) =      DOT_PRODUCT(                 LinAccEIMU, OtherState%CoordSys%c1 )
   OtherState%AllOuts(NcIMUTAys) = -1.0*DOT_PRODUCT(                 LinAccEIMU, OtherState%CoordSys%c3 )
   OtherState%AllOuts(NcIMUTAzs) =      DOT_PRODUCT(                 LinAccEIMU, OtherState%CoordSys%c2 )
   OtherState%AllOuts(NcIMURVxs) =      DOT_PRODUCT( OtherState%RtHS%AngVelER  , OtherState%CoordSys%c1 )*R2D
   OtherState%AllOuts(NcIMURVys) = -1.0*DOT_PRODUCT( OtherState%RtHS%AngVelER  , OtherState%CoordSys%c3 )*R2D
   OtherState%AllOuts(NcIMURVzs) =      DOT_PRODUCT( OtherState%RtHS%AngVelER  , OtherState%CoordSys%c2 )*R2D
   OtherState%AllOuts(NcIMURAxs) =      DOT_PRODUCT(                 AngAccER  , OtherState%CoordSys%c1 )*R2D
   OtherState%AllOuts(NcIMURAys) = -1.0*DOT_PRODUCT(                 AngAccER  , OtherState%CoordSys%c3 )*R2D
   OtherState%AllOuts(NcIMURAzs) =      DOT_PRODUCT(                 AngAccER  , OtherState%CoordSys%c2 )*R2D


      ! Rotor-Furl Motions:

   OtherState%AllOuts( RotFurlP) = x%QT  (DOF_RFrl)*R2D
   OtherState%AllOuts( RotFurlV) = x%QDT (DOF_RFrl)*R2D
   OtherState%AllOuts( RotFurlA) = OtherState%QD2T(DOF_RFrl)*R2D


      ! Tail-Furl Motions:

   OtherState%AllOuts(TailFurlP) = x%QT  (DOF_TFrl)*R2D
   OtherState%AllOuts(TailFurlV) = x%QDT (DOF_TFrl)*R2D
   OtherState%AllOuts(TailFurlA) = OtherState%QD2T(DOF_TFrl)*R2D


      ! Yaw Motions:

   OtherState%AllOuts(   YawPzn) = x%QT  (DOF_Yaw )*R2D
   OtherState%AllOuts(   YawVzn) = x%QDT (DOF_Yaw )*R2D
   OtherState%AllOuts(   YawAzn) = OtherState%QD2T(DOF_Yaw )*R2D


   ! Tower-Top / Yaw Bearing Motions:

   rOPO     = OtherState%RtHS%rT0O - p%TwrFlexL*OtherState%CoordSys%a2 ! Position vector from the undeflected tower top (point O prime) to the deflected tower top (point O).

   OtherState%AllOuts(YawBrTDxp) =  DOT_PRODUCT(     rOPO, OtherState%CoordSys%b1 )
   OtherState%AllOuts(YawBrTDyp) = -DOT_PRODUCT(     rOPO, OtherState%CoordSys%b3 )
   OtherState%AllOuts(YawBrTDzp) =  DOT_PRODUCT(     rOPO, OtherState%CoordSys%b2 )
   OtherState%AllOuts(YawBrTDxt) =  DOT_PRODUCT(     rOPO, OtherState%CoordSys%a1 )
   OtherState%AllOuts(YawBrTDyt) = -DOT_PRODUCT(     rOPO, OtherState%CoordSys%a3 )
   OtherState%AllOuts(YawBrTDzt) =  DOT_PRODUCT(     rOPO, OtherState%CoordSys%a2 )
   OtherState%AllOuts(YawBrTAxp) =  DOT_PRODUCT( LinAccEO, OtherState%CoordSys%b1 )
   OtherState%AllOuts(YawBrTAyp) = -DOT_PRODUCT( LinAccEO, OtherState%CoordSys%b3 )
   OtherState%AllOuts(YawBrTAzp) =  DOT_PRODUCT( LinAccEO, OtherState%CoordSys%b2 )
   OtherState%AllOuts(YawBrRDxt) =  DOT_PRODUCT( OtherState%RtHS%AngPosXB, OtherState%CoordSys%a1 )*R2D
   OtherState%AllOuts(YawBrRDyt) = -DOT_PRODUCT( OtherState%RtHS%AngPosXB, OtherState%CoordSys%a3 )*R2D
   ! There is no sense computing OtherState%AllOuts(YawBrRDzt) here since it is always zero for FAST simulation results.
   OtherState%AllOuts(YawBrRVxp) =  DOT_PRODUCT( OtherState%RtHS%AngVelEB, OtherState%CoordSys%b1 )*R2D
   OtherState%AllOuts(YawBrRVyp) = -DOT_PRODUCT( OtherState%RtHS%AngVelEB, OtherState%CoordSys%b3 )*R2D
   OtherState%AllOuts(YawBrRVzp) =  DOT_PRODUCT( OtherState%RtHS%AngVelEB, OtherState%CoordSys%b2 )*R2D
   OtherState%AllOuts(YawBrRAxp) =  DOT_PRODUCT( AngAccEB, OtherState%CoordSys%b1 )*R2D
   OtherState%AllOuts(YawBrRAyp) = -DOT_PRODUCT( AngAccEB, OtherState%CoordSys%b3 )*R2D
   OtherState%AllOuts(YawBrRAzp) =  DOT_PRODUCT( AngAccEB, OtherState%CoordSys%b2 )*R2D


      ! Local Tower Motions:

   DO I = 1, p%NTwGages

      OtherState%AllOuts( TwHtALxt(I) ) =      DOT_PRODUCT( LinAccET(:,p%TwrGagNd(I)), OtherState%CoordSys%t1(p%TwrGagNd(I),:) )
      OtherState%AllOuts( TwHtALyt(I) ) = -1.0*DOT_PRODUCT( LinAccET(:,p%TwrGagNd(I)), OtherState%CoordSys%t3(p%TwrGagNd(I),:) )
      OtherState%AllOuts( TwHtALzt(I) ) =      DOT_PRODUCT( LinAccET(:,p%TwrGagNd(I)), OtherState%CoordSys%t2(p%TwrGagNd(I),:) )

      rTPT                   = OtherState%RtHS%rT0T(:,p%TwrGagNd(I)) - p%HNodes(p%TwrGagNd(I))*OtherState%CoordSys%a2(:)

      OtherState%AllOuts( TwHtTDxt(I) ) =      DOT_PRODUCT( rTPT,     OtherState%CoordSys%a1 )
      OtherState%AllOuts( TwHtTDyt(I) ) = -1.0*DOT_PRODUCT( rTPT,     OtherState%CoordSys%a3 )
      OtherState%AllOuts( TwHtTDzt(I) ) =      DOT_PRODUCT( rTPT,     OtherState%CoordSys%a2 )

      OtherState%AllOuts( TwHtRDxt(I) ) =      DOT_PRODUCT( OtherState%RtHS%AngPosXF(:,p%TwrGagNd(I)), OtherState%CoordSys%a1 )*R2D  !why is this zero???
      OtherState%AllOuts( TwHtRDyt(I) ) = -1.0*DOT_PRODUCT( OtherState%RtHS%AngPosXF(:,p%TwrGagNd(I)), OtherState%CoordSys%a3 )*R2D
   !   OtherState%AllOuts( TwHtRDzt(I) ) =     DOT_PRODUCT( OtherState%RtHS%AngPosXF(:,p%TwrGagNd(I)), OtherState%CoordSys%a2 )*R2D  !this will always be 0 in FAST, so no need to calculate


      OtherState%AllOuts( TwHtTPxi(I) ) =      OtherState%RtHS%rT(1,p%TwrGagNd(I))
      OtherState%AllOuts( TwHtTPyi(I) ) = -1.0*OtherState%RtHS%rT(3,p%TwrGagNd(I))
      OtherState%AllOuts( TwHtTPzi(I) ) =      OtherState%RtHS%rT(2,p%TwrGagNd(I)) + p%PtfmRefzt

      OtherState%AllOuts( TwHtRPxi(I) ) =  OtherState%RtHS%AngPosEF(1,p%TwrGagNd(I))*R2D
      OtherState%AllOuts( TwHtRPyi(I) ) = -OtherState%RtHS%AngPosEF(3,p%TwrGagNd(I))*R2D
      OtherState%AllOuts( TwHtRPzi(I) ) =  OtherState%RtHS%AngPosEF(2,p%TwrGagNd(I))*R2D

   END DO !I

      ! Platform Motions:

   OtherState%AllOuts( PtfmTDxt) =  DOT_PRODUCT(       OtherState%RtHS%rZ, OtherState%CoordSys%a1 )
   OtherState%AllOuts( PtfmTDyt) = -DOT_PRODUCT(       OtherState%RtHS%rZ, OtherState%CoordSys%a3 )
   OtherState%AllOuts( PtfmTDzt) =  DOT_PRODUCT(       OtherState%RtHS%rZ, OtherState%CoordSys%a2 )
   OtherState%AllOuts( PtfmTDxi) = x%QT  (DOF_Sg )
   OtherState%AllOuts( PtfmTDyi) = x%QT  (DOF_Sw )
   OtherState%AllOuts( PtfmTDzi) = x%QT  (DOF_Hv )
   OtherState%AllOuts( PtfmTVxt) =  DOT_PRODUCT( OtherState%RtHS%LinVelEZ, OtherState%CoordSys%a1 )
   OtherState%AllOuts( PtfmTVyt) = -DOT_PRODUCT( OtherState%RtHS%LinVelEZ, OtherState%CoordSys%a3 )
   OtherState%AllOuts( PtfmTVzt) =  DOT_PRODUCT( OtherState%RtHS%LinVelEZ, OtherState%CoordSys%a2 )
   OtherState%AllOuts( PtfmTVxi) = x%QDT (DOF_Sg )
   OtherState%AllOuts( PtfmTVyi) = x%QDT (DOF_Sw )
   OtherState%AllOuts( PtfmTVzi) = x%QDT (DOF_Hv )
   OtherState%AllOuts( PtfmTAxt) =  DOT_PRODUCT(                 LinAccEZ, OtherState%CoordSys%a1 )
   OtherState%AllOuts( PtfmTAyt) = -DOT_PRODUCT(                 LinAccEZ, OtherState%CoordSys%a3 )
   OtherState%AllOuts( PtfmTAzt) =  DOT_PRODUCT(                 LinAccEZ, OtherState%CoordSys%a2 )
   OtherState%AllOuts( PtfmTAxi) = OtherState%QD2T(DOF_Sg  )
   OtherState%AllOuts( PtfmTAyi) = OtherState%QD2T(DOF_Sw  )
   OtherState%AllOuts( PtfmTAzi) = OtherState%QD2T(DOF_Hv  )
   OtherState%AllOuts( PtfmRDxi) = x%QT  (DOF_R )*R2D
   OtherState%AllOuts( PtfmRDyi) = x%QT  (DOF_P )*R2D
   OtherState%AllOuts( PtfmRDzi) = x%QT  (DOF_Y )*R2D
   OtherState%AllOuts( PtfmRVxt) =  DOT_PRODUCT( OtherState%RtHS%AngVelEX, OtherState%CoordSys%a1 )*R2D
   OtherState%AllOuts( PtfmRVyt) = -DOT_PRODUCT( OtherState%RtHS%AngVelEX, OtherState%CoordSys%a3 )*R2D
   OtherState%AllOuts( PtfmRVzt) =  DOT_PRODUCT( OtherState%RtHS%AngVelEX, OtherState%CoordSys%a2 )*R2D
   OtherState%AllOuts( PtfmRVxi) = x%QDT (DOF_R )*R2D
   OtherState%AllOuts( PtfmRVyi) = x%QDT (DOF_P )*R2D
   OtherState%AllOuts( PtfmRVzi) = x%QDT (DOF_Y )*R2D
   OtherState%AllOuts( PtfmRAxt) =  DOT_PRODUCT(                 AngAccEX, OtherState%CoordSys%a1 )*R2D
   OtherState%AllOuts( PtfmRAyt) = -DOT_PRODUCT(                 AngAccEX, OtherState%CoordSys%a3 )*R2D
   OtherState%AllOuts( PtfmRAzt) =  DOT_PRODUCT(                 AngAccEX, OtherState%CoordSys%a2 )*R2D
   OtherState%AllOuts( PtfmRAxi) = OtherState%QD2T(DOF_R )*R2D
   OtherState%AllOuts( PtfmRAyi) = OtherState%QD2T(DOF_P )*R2D
   OtherState%AllOuts( PtfmRAzi) = OtherState%QD2T(DOF_Y )*R2D



      ! Blade Root Loads:

   DO K=1,p%NumBl
      OtherState%AllOuts( RootFxc(K) ) = DOT_PRODUCT( FrcS0B(:,K), OtherState%CoordSys%i1(K,:) )
      OtherState%AllOuts( RootFyc(K) ) = DOT_PRODUCT( FrcS0B(:,K), OtherState%CoordSys%i2(K,:) )
      OtherState%AllOuts( RootFzc(K) ) = DOT_PRODUCT( FrcS0B(:,K), OtherState%CoordSys%i3(K,:) )
      OtherState%AllOuts( RootFxb(K) ) = DOT_PRODUCT( FrcS0B(:,K), OtherState%CoordSys%j1(K,:) )
      OtherState%AllOuts( RootFyb(K) ) = DOT_PRODUCT( FrcS0B(:,K), OtherState%CoordSys%j2(K,:) )
      OtherState%AllOuts( RootMxc(K) ) = DOT_PRODUCT( MomH0B(:,K), OtherState%CoordSys%i1(K,:) )
      OtherState%AllOuts( RootMyc(K) ) = DOT_PRODUCT( MomH0B(:,K), OtherState%CoordSys%i2(K,:) )
      OtherState%AllOuts( RootMzc(K) ) = DOT_PRODUCT( MomH0B(:,K), OtherState%CoordSys%i3(K,:) )
      OtherState%AllOuts( RootMxb(K) ) = DOT_PRODUCT( MomH0B(:,K), OtherState%CoordSys%j1(K,:) )
      OtherState%AllOuts( RootMyb(K) ) = DOT_PRODUCT( MomH0B(:,K), OtherState%CoordSys%j2(K,:) )
   END DO !K


      ! Blade Local Span Loads:

   DO K = 1,p%NumBl
      DO I = 1,p%NBlGages

      ! Initialize FrcMGagB and MomMGagB using the tip brake effects:

         FrcMGagB = OtherState%RtHS%FSTipDrag(:,K) - p%TipMass(K)*( p%Gravity*OtherState%CoordSys%z2 + LinAccES(:,p%TipNode,K) )
         MomMGagB = CROSS_PRODUCT( OtherState%RtHS%rS0S(:,K,p%TipNode) - OtherState%RtHS%rS0S(:,K,p%BldGagNd(I)), FrcMGagB )

      ! Integrate to find FrcMGagB and MomMGagB using all of the nodes / elements above the current strain gage location:
         DO J = ( p%BldGagNd(I) + 1 ),p%BldNodes ! Loop through blade nodes / elements above strain gage node

            TmpVec2  = OtherState%RtHS%FSAero(:,K,J) - p%MassB(K,J)*( p%Gravity*OtherState%CoordSys%z2 + LinAccES(:,J,K) )  ! Portion of FrcMGagB associated with element J
            FrcMGagB = FrcMGagB + TmpVec2*p%DRNodes(J)

            TmpVec = CROSS_PRODUCT( OtherState%RtHS%rS0S(:,K,J) - OtherState%RtHS%rS0S(:,K,p%BldGagNd(I)), TmpVec2 )           ! Portion of MomMGagB associated with element J
            MomMGagB = MomMGagB + ( TmpVec + OtherState%RtHS%MMAero(:,K,J) )*p%DRNodes(J)

         ENDDO ! J - Blade nodes / elements above strain gage node

      ! Add the effects of 1/2 the strain gage element:
      ! NOTE: for the radius in this calculation, assume that there is no
      !   shortening effect (due to blade bending) within the element.  Thus,
      !   the moment arm for the force is 1/4 of p%DRNodes() and the element
      !   length is 1/2 of p%DRNodes().

         TmpVec2  = OtherState%RtHS%FSAero(:,K,p%BldGagNd(I)) - p%MassB(K,p%BldGagNd(I))* ( p%Gravity*OtherState%CoordSys%z2 + LinAccES(:,p%BldGagNd(I),K) ) ! Portion of FrcMGagB associated with 1/2 of the strain gage element
         FrcMGagB = FrcMGagB + TmpVec2 * 0.5 * p%DRNodes(p%BldGagNd(I))                                                    ! Portion of FrcMGagB associated with 1/2 of the strain gage element
         FrcMGagB = 0.001*FrcMGagB           ! Convert the local force to kN


         TmpVec = CROSS_PRODUCT( ( 0.25_R8Ki*p%DRNodes(p%BldGagNd(I)) )*OtherState%CoordSys%j3(K,:), TmpVec2 )                              ! Portion of MomMGagB associated with 1/2 of the strain gage element

         MomMGagB = MomMGagB + ( TmpVec + OtherState%RtHS%MMAero(:,K,p%BldGagNd(I)) )* ( 0.5 *p%DRNodes(p%BldGagNd(I)) )
         MomMGagB = 0.001*MomMGagB           ! Convert the local moment to kN-m


         OtherState%AllOuts(SpnFLxb(I,K)) = DOT_PRODUCT( FrcMGagB, OtherState%CoordSys%n1(K,p%BldGagNd(I),:) )
         OtherState%AllOuts(SpnFLyb(I,K)) = DOT_PRODUCT( FrcMGagB, OtherState%CoordSys%n2(K,p%BldGagNd(I),:) )
         OtherState%AllOuts(SpnFLzb(I,K)) = DOT_PRODUCT( FrcMGagB, OtherState%CoordSys%n3(K,p%BldGagNd(I),:) )

         OtherState%AllOuts(SpnMLxb(I,K)) = DOT_PRODUCT( MomMGagB, OtherState%CoordSys%n1(K,p%BldGagNd(I),:) )
         OtherState%AllOuts(SpnMLyb(I,K)) = DOT_PRODUCT( MomMGagB, OtherState%CoordSys%n2(K,p%BldGagNd(I),:) )
         OtherState%AllOuts(SpnMLzb(I,K)) = DOT_PRODUCT( MomMGagB, OtherState%CoordSys%n3(K,p%BldGagNd(I),:) )
      END DO ! I
   END DO ! K



      ! Hub and Rotor Loads:

   !ComDenom = 0.5*p%AirDens*p%ProjArea*OtherState%AllOuts(  WindVxi)*OtherState%AllOuts(  WindVxi)   ! Common denominator used in several expressions

   OtherState%AllOuts(LSShftFxa) =  DOT_PRODUCT(  FrcPRot, OtherState%CoordSys%e1 )
   OtherState%AllOuts(LSShftFya) =  DOT_PRODUCT(  FrcPRot, OtherState%CoordSys%e2 )
   OtherState%AllOuts(LSShftFza) =  DOT_PRODUCT(  FrcPRot, OtherState%CoordSys%e3 )
   OtherState%AllOuts(LSShftFys) = -DOT_PRODUCT(  FrcPRot, OtherState%CoordSys%c3 )
   OtherState%AllOuts(LSShftFzs) =  DOT_PRODUCT(  FrcPRot, OtherState%CoordSys%c2 )
   OtherState%AllOuts(LSShftMxa) =  DOT_PRODUCT( MomLPRot, OtherState%CoordSys%e1 )
   OtherState%AllOuts(LSSTipMya) =  DOT_PRODUCT( MomLPRot, OtherState%CoordSys%e2 )
   OtherState%AllOuts(LSSTipMza) =  DOT_PRODUCT( MomLPRot, OtherState%CoordSys%e3 )
   OtherState%AllOuts(LSSTipMys) = -DOT_PRODUCT( MomLPRot, OtherState%CoordSys%c3 )
   OtherState%AllOuts(LSSTipMzs) =  DOT_PRODUCT( MomLPRot, OtherState%CoordSys%c2 )

!   IF ( .NOT. EqualRealNos( OtherState%AllOuts(LSShftFxa), 0.0_ReKi ) )   THEN ! .TRUE. if the denominator in the following equations is not zero.
!
!      CThrstys = -OtherState%AllOuts(LSSTipMzs)/OtherState%AllOuts(LSShftFxa)  ! Estimate of the ys-location of the center of thrust
!      CThrstzs =  OtherState%AllOuts(LSSTipMys)/OtherState%AllOuts(LSShftFxa)  ! Estimate of the zs-location of the center of thrust
!
!!      OtherState%AllOuts(CThrstAzm) = MOD( ( ATAN2( -CThrstzs, -CThrstys ) + p%AzimB1Up )*R2D + 360.0 + 90.0, 360.0 )  !bjj: IgnoreMod was used for linearization... perhaps these outputs should not use the MOD function; only WriteOutputs should have that...
!!      OtherState%AllOuts(CThrstRad) = MIN( 1.0, SQRT( CThrstys*CThrstys + CThrstzs*CThrstzs )/p%AvgNrmTpRd )
!
!   ELSE
!
!      !OtherState%AllOuts(CThrstAzm) = 0.0
!      !OtherState%AllOuts(CThrstRad) = 0.0
!
!   ENDIF

   OtherState%AllOuts(   RotPwr) = ( x%QDT(DOF_GeAz) + x%QDT(DOF_DrTr) )*OtherState%AllOuts(LSShftMxa)

   !IF ( .NOT. EqualRealNos( ComDenom, 0.0_ReKi ) )  THEN   ! .TRUE. if the denominator in the following equations is not zero.
   !
   !   OtherState%AllOuts( RotCq) = 1000.0*OtherState%AllOuts(LSShftMxa) / ( ComDenom*p%TipRad )
   !   OtherState%AllOuts( RotCp) = 1000.0*OtherState%AllOuts(   RotPwr) / ( ComDenom*OtherState%AllOuts(  WindVxi) )
   !   OtherState%AllOuts( RotCt) = 1000.0*OtherState%AllOuts(LSShftFxa) /   ComDenom
   !
   !ELSE
   !
   !   OtherState%AllOuts( RotCq) = 0.0
   !   OtherState%AllOuts( RotCp) = 0.0
   !   OtherState%AllOuts( RotCt) = 0.0
   !
   !ENDIF


      ! Shaft Strain Gage Loads:

   OtherState%AllOuts(LSSGagMya) = OtherState%AllOuts(LSSTipMya) + p%ShftGagL*OtherState%AllOuts(LSShftFza)
   OtherState%AllOuts(LSSGagMza) = OtherState%AllOuts(LSSTipMza) - p%ShftGagL*OtherState%AllOuts(LSShftFya)
   OtherState%AllOuts(LSSGagMys) = OtherState%AllOuts(LSSTipMys) + p%ShftGagL*OtherState%AllOuts(LSShftFzs)
   OtherState%AllOuts(LSSGagMzs) = OtherState%AllOuts(LSSTipMzs) - p%ShftGagL*OtherState%AllOuts(LSShftFys)


      ! Generator and High-Speed Shaft Loads:

   OtherState%AllOuts( HSShftTq)  = OtherState%AllOuts(LSShftMxa)*OtherState%RtHS%GBoxEffFac/ABS(p%GBRatio)
   OtherState%AllOuts(HSShftPwr)  = OtherState%AllOuts( HSShftTq)*ABS(p%GBRatio)*x%QDT(DOF_GeAz)
   OtherState%AllOuts(HSSBrTq)    = OtherState%HSSBrTrq*0.001_ReKi


   !IF ( .NOT. EqualRealNos( ComDenom, 0.0_ReKi ) )  THEN  ! .TRUE. if the denominator in the following equations is not zero (ComDenom is the same as it is calculated above).
   !
   !   OtherState%AllOuts( HSShftCq) = 1000.0*OtherState%AllOuts( HSShftTq) / ( ComDenom*p%TipRad )
   !   OtherState%AllOuts( HSShftCp) = 1000.0*OtherState%AllOuts(HSShftPwr) / ( ComDenom*OtherState%AllOuts(  WindVxi) )
   !   OtherState%AllOuts(    GenCq) = 1000.0*OtherState%AllOuts(    GenTq) / ( ComDenom*p%TipRad )
   !   OtherState%AllOuts(    GenCp) = 1000.0*OtherState%AllOuts(   GenPwr) / ( ComDenom*OtherState%AllOuts(  WindVxi) )
   !
   !ELSE
   !
   !   OtherState%AllOuts( HSShftCq) = 0.0
   !   OtherState%AllOuts( HSShftCp) = 0.0
   !   OtherState%AllOuts(    GenCq) = 0.0
   !   OtherState%AllOuts(    GenCp) = 0.0
   !
   !ENDIF


      ! Rotor-Furl Axis Loads:

   OtherState%AllOuts(RFrlBrM  ) =  DOT_PRODUCT( MomNGnRt, OtherState%CoordSys%rfa )


      ! Tail-Furl Axis Loads:

   OtherState%AllOuts(TFrlBrM  ) =  DOT_PRODUCT( MomNTail, OtherState%CoordSys%tfa )


      ! Tower-Top / Yaw Bearing Loads:

   OtherState%AllOuts( YawBrFxn) =  DOT_PRODUCT( FrcONcRt, OtherState%CoordSys%d1 )
   OtherState%AllOuts( YawBrFyn) = -DOT_PRODUCT( FrcONcRt, OtherState%CoordSys%d3 )
   OtherState%AllOuts( YawBrFzn) =  DOT_PRODUCT( FrcONcRt, OtherState%CoordSys%d2 )
   OtherState%AllOuts( YawBrFxp) =  DOT_PRODUCT( FrcONcRt, OtherState%CoordSys%b1 )
   OtherState%AllOuts( YawBrFyp) = -DOT_PRODUCT( FrcONcRt, OtherState%CoordSys%b3 )
   OtherState%AllOuts( YawBrMxn) =  DOT_PRODUCT( MomBNcRt, OtherState%CoordSys%d1 )
   OtherState%AllOuts( YawBrMyn) = -DOT_PRODUCT( MomBNcRt, OtherState%CoordSys%d3 )
   OtherState%AllOuts( YawBrMzn) =  DOT_PRODUCT( MomBNcRt, OtherState%CoordSys%d2 )
   OtherState%AllOuts( YawBrMxp) =  DOT_PRODUCT( MomBNcRt, OtherState%CoordSys%b1 )
   OtherState%AllOuts( YawBrMyp) = -DOT_PRODUCT( MomBNcRt, OtherState%CoordSys%b3 )


      ! Tower Base Loads:

   OtherState%AllOuts( TwrBsFxt) =  DOT_PRODUCT( FrcT0Trb, OtherState%CoordSys%a1 )
   OtherState%AllOuts( TwrBsFyt) = -DOT_PRODUCT( FrcT0Trb, OtherState%CoordSys%a3 )
   OtherState%AllOuts( TwrBsFzt) =  DOT_PRODUCT( FrcT0Trb, OtherState%CoordSys%a2 )
   OtherState%AllOuts( TwrBsMxt) =  DOT_PRODUCT( MomX0Trb, OtherState%CoordSys%a1 )
   OtherState%AllOuts( TwrBsMyt) = -DOT_PRODUCT( MomX0Trb, OtherState%CoordSys%a3 )
   OtherState%AllOuts( TwrBsMzt) =  DOT_PRODUCT( MomX0Trb, OtherState%CoordSys%a2 )


      ! Local Tower Loads:

   FrcONcRt = 1000.0*FrcONcRt ! Convert the units of these forces and moments
   MomBNcRt = 1000.0*MomBNcRt ! from kN and kN-m back to N and N-m, respectively.

   DO I=1,p%NTwGages

      ! Initialize FrcFGagT and MomFGagT using the tower-top and yaw bearing mass effects:
      FrcFGagT = FrcONcRt - p%YawBrMass*( p%Gravity*OtherState%CoordSys%z2 + LinAccEO )
      MomFGagT = CROSS_PRODUCT( OtherState%RtHS%rZO - OtherState%RtHS%rZT(:,p%TwrGagNd(I)), FrcFGagT )
      MomFGagT = MomFGagT + MomBNcRt

      ! Integrate to find FrcFGagT and MomFGagT using all of the nodes / elements above the current strain gage location:
      DO J = ( p%TwrGagNd(I) + 1 ),p%TwrNodes ! Loop through tower nodes / elements above strain gage node
         TmpVec2  = FTTower(:,J) - p%MassT(J)*( p%Gravity*OtherState%CoordSys%z2 + LinAccET(:,J) )           ! Portion of FrcFGagT associated with element J
         FrcFGagT = FrcFGagT + TmpVec2*p%DHNodes(J)

         TmpVec = CROSS_PRODUCT( OtherState%RtHS%rZT(:,J) - OtherState%RtHS%rZT(:,p%TwrGagNd(I)), TmpVec2 )                          ! Portion of MomFGagT associated with element J
         MomFGagT = MomFGagT + ( TmpVec + MFHydro(:,J) )*p%DHNodes(J)
      ENDDO ! J -Tower nodes / elements above strain gage node

      ! Add the effects of 1/2 the strain gage element:
      ! NOTE: for the radius in this calculation, assume that there is no shortening
      !   effect (due to tower bending) within the element.  Thus, the moment arm
      !   for the force is 1/4 of DHNodes() and the element length is 1/2 of DHNodes().

      TmpVec2  = FTTower(:,p%TwrGagNd(I)) - p%MassT(p%TwrGagNd(I))*( p%Gravity*OtherState%CoordSys%z2 + LinAccET(:,p%TwrGagNd(I)))

      FrcFGagT = FrcFGagT + TmpVec2 * 0.5 * p%DHNodes(p%TwrGagNd(I))
      FrcFGagT = 0.001*FrcFGagT  ! Convert the local force to kN

      TmpVec = CROSS_PRODUCT( ( 0.25_R8Ki*p%DHNodes( p%TwrGagNd(I)) )*OtherState%CoordSys%a2, TmpVec2 )              ! Portion of MomFGagT associated with 1/2 of the strain gage element
      TmpVec   = TmpVec   + MFHydro(:,p%TwrGagNd(I))
      MomFGagT = MomFGagT + TmpVec * 0.5 * p%DHNodes(p%TwrGagNd(I))
      MomFGagT = 0.001*MomFGagT  ! Convert the local moment to kN-m

      OtherState%AllOuts( TwHtFLxt(I) ) =     DOT_PRODUCT( FrcFGagT, OtherState%CoordSys%t1(p%TwrGagNd(I),:) )
      OtherState%AllOuts( TwHtFLyt(I) ) = -1.*DOT_PRODUCT( FrcFGagT, OtherState%CoordSys%t3(p%TwrGagNd(I),:) )
      OtherState%AllOuts( TwHtFLzt(I) ) =     DOT_PRODUCT( FrcFGagT, OtherState%CoordSys%t2(p%TwrGagNd(I),:) )

      OtherState%AllOuts( TwHtMLxt(I) ) =     DOT_PRODUCT( MomFGagT, OtherState%CoordSys%t1(p%TwrGagNd(I),:) )
      OtherState%AllOuts( TwHtMLyt(I) ) = -1.*DOT_PRODUCT( MomFGagT, OtherState%CoordSys%t3(p%TwrGagNd(I),:) )
      OtherState%AllOuts( TwHtMLzt(I) ) =     DOT_PRODUCT( MomFGagT, OtherState%CoordSys%t2(p%TwrGagNd(I),:) )

   END DO


   !   ! Platform Loads:
   !
   !OtherState%AllOuts(  PtfmFxt) =  DOT_PRODUCT( FZHydro, OtherState%CoordSys%a1 )
   !OtherState%AllOuts(  PtfmFyt) = -DOT_PRODUCT( FZHydro, OtherState%CoordSys%a3 )
   !OtherState%AllOuts(  PtfmFzt) =  DOT_PRODUCT( FZHydro, OtherState%CoordSys%a2 )
   !OtherState%AllOuts(  PtfmFxi) =  DOT_PRODUCT( FZHydro, OtherState%CoordSys%z1 )
   !OtherState%AllOuts(  PtfmFyi) = -DOT_PRODUCT( FZHydro, OtherState%CoordSys%z3 )
   !OtherState%AllOuts(  PtfmFzi) =  DOT_PRODUCT( FZHydro, OtherState%CoordSys%z2 )
   !OtherState%AllOuts(  PtfmMxt) =  DOT_PRODUCT( MXHydro, OtherState%CoordSys%a1 )
   !OtherState%AllOuts(  PtfmMyt) = -DOT_PRODUCT( MXHydro, OtherState%CoordSys%a3 )
   !OtherState%AllOuts(  PtfmMzt) =  DOT_PRODUCT( MXHydro, OtherState%CoordSys%a2 )
   !OtherState%AllOuts(  PtfmMxi) =  DOT_PRODUCT( MXHydro, OtherState%CoordSys%z1 )
   !OtherState%AllOuts(  PtfmMyi) = -DOT_PRODUCT( MXHydro, OtherState%CoordSys%z3 )
   !OtherState%AllOuts(  PtfmMzi) =  DOT_PRODUCT( MXHydro, OtherState%CoordSys%z2 )
   !

      ! Internal p%DOFs outputs:

   OtherState%AllOuts( Q_B1E1   ) = x%QT(   DOF_BE(1,1) )
   OtherState%AllOuts( Q_B2E1   ) = x%QT(   DOF_BE(2,1) )
   OtherState%AllOuts( Q_B1F1   ) = x%QT(   DOF_BF(1,1) )
   OtherState%AllOuts( Q_B2F1   ) = x%QT(   DOF_BF(2,1) )
   OtherState%AllOuts( Q_B1F2   ) = x%QT(   DOF_BF(1,2) )
   OtherState%AllOuts( Q_B2F2   ) = x%QT(   DOF_BF(2,2) )
   OtherState%AllOuts( Q_DrTr   ) = x%QT(   DOF_DrTr    )
   OtherState%AllOuts( Q_GeAz   ) = x%QT(   DOF_GeAz    )
   OtherState%AllOuts( Q_RFrl   ) = x%QT(   DOF_RFrl    )
   OtherState%AllOuts( Q_TFrl   ) = x%QT(   DOF_TFrl    )
   OtherState%AllOuts( Q_Yaw    ) = x%QT(   DOF_Yaw     )
   OtherState%AllOuts( Q_TFA1   ) = x%QT(   DOF_TFA1    )
   OtherState%AllOuts( Q_TSS1   ) = x%QT(   DOF_TSS1    )
   OtherState%AllOuts( Q_TFA2   ) = x%QT(   DOF_TFA2    )
   OtherState%AllOuts( Q_TSS2   ) = x%QT(   DOF_TSS2    )
   OtherState%AllOuts( Q_Sg     ) = x%QT(   DOF_Sg      )
   OtherState%AllOuts( Q_Sw     ) = x%QT(   DOF_Sw      )
   OtherState%AllOuts( Q_Hv     ) = x%QT(   DOF_Hv      )
   OtherState%AllOuts( Q_R      ) = x%QT(   DOF_R       )
   OtherState%AllOuts( Q_P      ) = x%QT(   DOF_P       )
   OtherState%AllOuts( Q_Y      ) = x%QT(   DOF_Y       )

   OtherState%AllOuts( QD_B1E1  ) = x%QDT(  DOF_BE(1,1) )
   OtherState%AllOuts( QD_B2E1  ) = x%QDT(  DOF_BE(2,1) )
   OtherState%AllOuts( QD_B1F1  ) = x%QDT(  DOF_BF(1,1) )
   OtherState%AllOuts( QD_B2F1  ) = x%QDT(  DOF_BF(2,1) )
   OtherState%AllOuts( QD_B1F2  ) = x%QDT(  DOF_BF(1,2) )
   OtherState%AllOuts( QD_B2F2  ) = x%QDT(  DOF_BF(2,2) )
   OtherState%AllOuts( QD_DrTr  ) = x%QDT(  DOF_DrTr    )
   OtherState%AllOuts( QD_GeAz  ) = x%QDT(  DOF_GeAz    )
   OtherState%AllOuts( QD_RFrl  ) = x%QDT(  DOF_RFrl    )
   OtherState%AllOuts( QD_TFrl  ) = x%QDT(  DOF_TFrl    )
   OtherState%AllOuts( QD_Yaw   ) = x%QDT(  DOF_Yaw     )
   OtherState%AllOuts( QD_TFA1  ) = x%QDT(  DOF_TFA1    )
   OtherState%AllOuts( QD_TSS1  ) = x%QDT(  DOF_TSS1    )
   OtherState%AllOuts( QD_TFA2  ) = x%QDT(  DOF_TFA2    )
   OtherState%AllOuts( QD_TSS2  ) = x%QDT(  DOF_TSS2    )
   OtherState%AllOuts( QD_Sg    ) = x%QDT(  DOF_Sg      )
   OtherState%AllOuts( QD_Sw    ) = x%QDT(  DOF_Sw      )
   OtherState%AllOuts( QD_Hv    ) = x%QDT(  DOF_Hv      )
   OtherState%AllOuts( QD_R     ) = x%QDT(  DOF_R       )
   OtherState%AllOuts( QD_P     ) = x%QDT(  DOF_P       )
   OtherState%AllOuts( QD_Y     ) = x%QDT(  DOF_Y       )

   OtherState%AllOuts( QD2_B1E1 ) = OtherState%QD2T( DOF_BE(1,1) )
   OtherState%AllOuts( QD2_B2E1 ) = OtherState%QD2T( DOF_BE(2,1) )
   OtherState%AllOuts( QD2_B1F1 ) = OtherState%QD2T( DOF_BF(1,1) )
   OtherState%AllOuts( QD2_B2F1 ) = OtherState%QD2T( DOF_BF(2,1) )
   OtherState%AllOuts( QD2_B1F2 ) = OtherState%QD2T( DOF_BF(1,2) )
   OtherState%AllOuts( QD2_B2F2 ) = OtherState%QD2T( DOF_BF(2,2) )
   OtherState%AllOuts( QD2_DrTr ) = OtherState%QD2T( DOF_DrTr    )
   OtherState%AllOuts( QD2_GeAz ) = OtherState%QD2T( DOF_GeAz    )
   OtherState%AllOuts( QD2_RFrl ) = OtherState%QD2T( DOF_RFrl    )
   OtherState%AllOuts( QD2_TFrl ) = OtherState%QD2T( DOF_TFrl    )
   OtherState%AllOuts( QD2_Yaw  ) = OtherState%QD2T( DOF_Yaw     )
   OtherState%AllOuts( QD2_TFA1 ) = OtherState%QD2T( DOF_TFA1    )
   OtherState%AllOuts( QD2_TSS1 ) = OtherState%QD2T( DOF_TSS1    )
   OtherState%AllOuts( QD2_TFA2 ) = OtherState%QD2T( DOF_TFA2    )
   OtherState%AllOuts( QD2_TSS2 ) = OtherState%QD2T( DOF_TSS2    )
   OtherState%AllOuts( QD2_Sg   ) = OtherState%QD2T( DOF_Sg      )
   OtherState%AllOuts( QD2_Sw   ) = OtherState%QD2T( DOF_Sw      )
   OtherState%AllOuts( QD2_Hv   ) = OtherState%QD2T( DOF_Hv      )
   OtherState%AllOuts( QD2_R    ) = OtherState%QD2T( DOF_R       )
   OtherState%AllOuts( QD2_P    ) = OtherState%QD2T( DOF_P       )
   OtherState%AllOuts( QD2_Y    ) = OtherState%QD2T( DOF_Y       )


   IF ( p%NumBl > 2 ) THEN
      OtherState%AllOuts( Q_B3E1   ) = x%QT(   DOF_BE(3,1) )
      OtherState%AllOuts( Q_B3F1   ) = x%QT(   DOF_BF(3,1) )
      OtherState%AllOuts( Q_B3F2   ) = x%QT(   DOF_BF(3,2) )

      OtherState%AllOuts( QD_B3E1  ) = x%QDT(  DOF_BE(3,1) )
      OtherState%AllOuts( QD_B3F1  ) = x%QDT(  DOF_BF(3,1) )
      OtherState%AllOuts( QD_B3F2  ) = x%QDT(  DOF_BF(3,2) )

      OtherState%AllOuts( QD2_B3E1 ) = OtherState%QD2T( DOF_BE(3,1) )
      OtherState%AllOuts( QD2_B3F1 ) = OtherState%QD2T( DOF_BF(3,1) )
      OtherState%AllOuts( QD2_B3F2 ) = OtherState%QD2T( DOF_BF(3,2) )
   ELSE
      OtherState%AllOuts( Q_Teet   ) = x%QT(            DOF_Teet    )
      OtherState%AllOuts( QD_Teet  ) = x%QDT(           DOF_Teet    )
      OtherState%AllOuts( QD2_Teet ) = OtherState%QD2T( DOF_Teet    )
   END IF

   !...............................................................................................................................
   ! Place the selected output channels into the WriteOutput(:) array with the proper sign:
   !...............................................................................................................................

   DO I = 1,p%NumOuts  ! Loop through all selected output channels

      y%WriteOutput(I) = p%OutParam(I)%SignM * OtherState%AllOuts( p%OutParam(I)%Indx )

   ENDDO             ! I - All selected output channels

   
   !...............................................................................................................................
   ! Outputs required for AeroDyn
   !...............................................................................................................................
   
   !JASON: WE SHOULD REALLY BE PASSING TO AERODYN THE LINEAR VELOCITIES OF THE AERODYNAMIC CENTER IN THE INERTIA FRAME, NOT SIMPLY THE LINEAR VELOCITIES OF POINT S.  
   !       IS THERE ANY WAY OF GETTING THIS VELOCITY?<--DO THIS, WHEN YOU ADD THE COUPLED MODE SHAPES!!!!
   
   !...........
   ! Blade elements:
   !...........
   IF ( ALLOCATED(y%BladeLn2Mesh) ) THEN
      DO K = 1,p%NumBl ! Loop through all blades
         DO J = 0,p%TipNode ! Loop through the blade nodes / elements
            
            J2 = J            
            if (j==0) then
                  ! blade root
               NodeNum = p%BldNodes + 2
               if (p%UseAD14) j2 = 1                  
            elseif (j==p%TipNode) then
               ! blade tip
               NodeNum = p%BldNodes + 1
               if (p%UseAD14) j2 = p%BldNodes
            else
               NodeNum = J
            end if
                                                                                                        
            if (p%UseAD14) then                  
                  ! Translational Displacement (first calculate absolute position)
               y%BladeLn2Mesh(K)%TranslationDisp(1,NodeNum) =     OtherState%RtHS%rS (1,K,J2) + OtherState%RtHS%rSAerCen(1,J2,K)               ! = the distance from the undeflected tower centerline                                     to the current blade aerodynamic center in the xi ( z1) direction
               y%BladeLn2Mesh(K)%TranslationDisp(2,NodeNum) = -1.*OtherState%RtHS%rS (3,K,J2) - OtherState%RtHS%rSAerCen(3,J2,K)               ! = the distance from the undeflected tower centerline                                     to the current blade aerodynamic center in the yi (-z3) direction
               y%BladeLn2Mesh(K)%TranslationDisp(3,NodeNum) =     OtherState%RtHS%rS (2,K,J2) + OtherState%RtHS%rSAerCen(2,J2,K) + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade aerodynamic center in the zi ( z2) direction
               
                  ! Orientation 
               y%BladeLn2Mesh(K)%Orientation(1,1,NodeNum) =     OtherState%CoordSys%te1(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(2,1,NodeNum) =     OtherState%CoordSys%te2(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(3,1,NodeNum) =     OtherState%CoordSys%te3(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(1,2,NodeNum) = -1.*OtherState%CoordSys%te1(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(2,2,NodeNum) = -1.*OtherState%CoordSys%te2(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(3,2,NodeNum) = -1.*OtherState%CoordSys%te3(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(1,3,NodeNum) =     OtherState%CoordSys%te1(K,J2,2)
               y%BladeLn2Mesh(K)%Orientation(2,3,NodeNum) =     OtherState%CoordSys%te2(K,J2,2)
               y%BladeLn2Mesh(K)%Orientation(3,3,NodeNum) =     OtherState%CoordSys%te3(K,J2,2)
               
                  ! Translational Acceleration (for water-power request for added mass calculations)
               y%BladeLn2Mesh(K)%TranslationAcc(1,NodeNum) =     LinAccES(1,J2,K)
               y%BladeLn2Mesh(K)%TranslationAcc(2,NodeNum) = -1.*LinAccES(3,J2,K)
               y%BladeLn2Mesh(K)%TranslationAcc(3,NodeNum) =     LinAccES(2,J2,K)  
               
            else         
                  ! Translational Displacement (first calculate absolute position)
               y%BladeLn2Mesh(K)%TranslationDisp(1,NodeNum) =     OtherState%RtHS%rS (1,K,J2)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
               y%BladeLn2Mesh(K)%TranslationDisp(2,NodeNum) = -1.*OtherState%RtHS%rS (3,K,J2)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
               y%BladeLn2Mesh(K)%TranslationDisp(3,NodeNum) =     OtherState%RtHS%rS (2,K,J2)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
               
                  ! Orientation
               y%BladeLn2Mesh(K)%Orientation(1,1,NodeNum) =     OtherState%CoordSys%n1(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(2,1,NodeNum) =     OtherState%CoordSys%n2(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(3,1,NodeNum) =     OtherState%CoordSys%n3(K,J2,1)
               y%BladeLn2Mesh(K)%Orientation(1,2,NodeNum) = -1.*OtherState%CoordSys%n1(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(2,2,NodeNum) = -1.*OtherState%CoordSys%n2(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(3,2,NodeNum) = -1.*OtherState%CoordSys%n3(K,J2,3)
               y%BladeLn2Mesh(K)%Orientation(1,3,NodeNum) =     OtherState%CoordSys%n1(K,J2,2)
               y%BladeLn2Mesh(K)%Orientation(2,3,NodeNum) =     OtherState%CoordSys%n2(K,J2,2)
               y%BladeLn2Mesh(K)%Orientation(3,3,NodeNum) =     OtherState%CoordSys%n3(K,J2,2)
            end if
            
               ! Translational Displacement (get displacement, not absolute position):
            y%BladeLn2Mesh(K)%TranslationDisp(:,NodeNum) = y%BladeLn2Mesh(K)%TranslationDisp(:,NodeNum) - y%BladeLn2Mesh(K)%Position(:,NodeNum)
            
           
               ! Translational Velocity
            y%BladeLn2Mesh(K)%TranslationVel(1,NodeNum) =     OtherState%RtHS%LinVelES(1,J2,K)
            y%BladeLn2Mesh(K)%TranslationVel(2,NodeNum) = -1.*OtherState%RtHS%LinVelES(3,J2,K)
            y%BladeLn2Mesh(K)%TranslationVel(3,NodeNum) =     OtherState%RtHS%LinVelES(2,J2,K)  
                                                
         END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements
                  
      END DO !K = 1,p%NumBl
   END IF
   
      
   !...........
   ! Hub (for Lidar and AeroDyn15):
   !...........   
   
         ! Translation (absolute position - starting position):
   y%HubPtMotion%TranslationDisp(1,1)  =     OtherState%RtHS%rQ(1)
   y%HubPtMotion%TranslationDisp(2,1)  = -1.*OtherState%RtHS%rQ(3)
   y%HubPtMotion%TranslationDisp(3,1)  =     OtherState%RtHS%rQ(2) + p%PtfmRefzt
   y%HubPtMotion%TranslationDisp       = y%HubPtMotion%TranslationDisp - y%HubPtMotion%Position   ! relative position
   
      ! Orientation:        
   y%HubPtMotion%Orientation(1,1,1)    =     OtherState%CoordSys%g1(1) 
   y%HubPtMotion%Orientation(2,1,1)    =     OtherState%CoordSys%g2(1)
   y%HubPtMotion%Orientation(3,1,1)    =     OtherState%CoordSys%g3(1)   
   y%HubPtMotion%Orientation(1,2,1)    = -1.*OtherState%CoordSys%g1(3)
   y%HubPtMotion%Orientation(2,2,1)    = -1.*OtherState%CoordSys%g2(3) 
   y%HubPtMotion%Orientation(3,2,1)    = -1.*OtherState%CoordSys%g3(3) 
   y%HubPtMotion%Orientation(1,3,1)    =     OtherState%CoordSys%g1(2)
   y%HubPtMotion%Orientation(2,3,1)    =     OtherState%CoordSys%g2(2)
   y%HubPtMotion%Orientation(3,3,1)    =     OtherState%CoordSys%g3(2)
   
      ! Rotational velocity:
   y%HubPtMotion%RotationVel(1,1)      =     OtherState%RtHS%AngVelEH(1)
   y%HubPtMotion%RotationVel(2,1)      = -1.*OtherState%RtHS%AngVelEH(3)
   y%HubPtMotion%RotationVel(3,1)      =     OtherState%RtHS%AngVelEH(2)   
   
   !...........
   ! Blade roots (BeamDyn/AeroDyn v15):
   !...........   
         
   DO K=1,p%NumBl
         
      ! Translation displacement  ! rS at the root      
      y%BladeRootMotion(K)%TranslationDisp(1,1) =            OtherState%RtHS%rS (1,K,0)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
      y%BladeRootMotion(K)%TranslationDisp(2,1) =        -1.*OtherState%RtHS%rS (3,K,0)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
      y%BladeRootMotion(K)%TranslationDisp(3,1) =            OtherState%RtHS%rS (2,K,0)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
      y%BladeRootMotion(K)%TranslationDisp      = y%BladeRootMotion(K)%TranslationDisp - y%BladeRootMotion(K)%Position ! make it relative
      
      
      ! Orientation 
      y%BladeRootMotion(K)%Orientation(1,1,1)   =     OtherState%CoordSys%j1(K,1)
      y%BladeRootMotion(K)%Orientation(2,1,1)   =     OtherState%CoordSys%j2(K,1)
      y%BladeRootMotion(K)%Orientation(3,1,1)   =     OtherState%CoordSys%j3(K,1)
      y%BladeRootMotion(K)%Orientation(1,2,1)   = -1.*OtherState%CoordSys%j1(K,3)
      y%BladeRootMotion(K)%Orientation(2,2,1)   = -1.*OtherState%CoordSys%j2(K,3)
      y%BladeRootMotion(K)%Orientation(3,2,1)   = -1.*OtherState%CoordSys%j3(K,3)
      y%BladeRootMotion(K)%Orientation(1,3,1)   =     OtherState%CoordSys%j1(K,2)
      y%BladeRootMotion(K)%Orientation(2,3,1)   =     OtherState%CoordSys%j2(K,2)
      y%BladeRootMotion(K)%Orientation(3,3,1)   =     OtherState%CoordSys%j3(K,2)

      ! Translation velocity 
      y%BladeRootMotion(K)%TranslationVel(1,1)  =     OtherState%RtHS%LinVelES(1,0,K)
      y%BladeRootMotion(K)%TranslationVel(2,1)  = -1.*OtherState%RtHS%LinVelES(3,0,K)
      y%BladeRootMotion(K)%TranslationVel(3,1)  =     OtherState%RtHS%LinVelES(2,0,K)

      ! Rotation velocity  
      y%BladeRootMotion(K)%RotationVel(1,1)     =      OtherState%RtHS%AngVelEH(1)
      y%BladeRootMotion(K)%RotationVel(2,1)     =  -1.*OtherState%RtHS%AngVelEH(3)
      y%BladeRootMotion(K)%RotationVel(3,1)     =      OtherState%RtHS%AngVelEH(2)
      
      ! Translation acceleration
      y%BladeRootMotion(K)%TranslationAcc(1,1)  =      LinAccES(1,0,K)
      y%BladeRootMotion(K)%TranslationAcc(2,1)  =  -1.*LinAccES(3,0,K)
      y%BladeRootMotion(K)%TranslationAcc(3,1)  =      LinAccES(2,0,K)
      
      ! Rotation acceleration  
      y%BladeRootMotion(K)%RotationAcc(1,1)     =      AngAccEH(1) 
      y%BladeRootMotion(K)%RotationAcc(2,1)     =  -1.*AngAccEH(3) 
      y%BladeRootMotion(K)%RotationAcc(3,1)     =      AngAccEH(2)
      
   END DO   
   
   !...........
   ! Hub (for AeroDyn v14):
   !...........   

      ! the hub position should use rQ instead of rP, but AeroDyn 14 treats
      ! teeter deflections like blade deflections:
   
   y%HubPtMotion14%TranslationDisp(1,1)  =     OtherState%RtHS%rP(1)
   y%HubPtMotion14%TranslationDisp(2,1)  = -1.*OtherState%RtHS%rP(3)
   y%HubPtMotion14%TranslationDisp(3,1)  =     OtherState%RtHS%rP(2) + p%PtfmRefzt
   
   y%HubPtMotion14%TranslationDisp  = y%HubPtMotion14%TranslationDisp - y%HubPtMotion14%Position   
   
        ! Hub orientation should use the g instead of e system, but the current version
        ! of AeroDyn calculates forces normal and tangential to the cone of rotation
         
   y%HubPtMotion14%Orientation(1,1,1)  =     OtherState%CoordSys%e1(1) 
   y%HubPtMotion14%Orientation(2,1,1)  =     OtherState%CoordSys%e2(1)
   y%HubPtMotion14%Orientation(3,1,1)  =     OtherState%CoordSys%e3(1)   
   y%HubPtMotion14%Orientation(1,2,1)  = -1.*OtherState%CoordSys%e1(3)
   y%HubPtMotion14%Orientation(2,2,1)  = -1.*OtherState%CoordSys%e2(3) 
   y%HubPtMotion14%Orientation(3,2,1)  = -1.*OtherState%CoordSys%e3(3) 
   y%HubPtMotion14%Orientation(1,3,1)  =     OtherState%CoordSys%e1(2)
   y%HubPtMotion14%Orientation(2,3,1)  =     OtherState%CoordSys%e2(2)
   y%HubPtMotion14%Orientation(3,3,1)  =     OtherState%CoordSys%e3(2)
   
      ! Note the hub rotational velocity should be AngVelEH instead AngVelEL, but AeroDyn (13.00.00)
      ! treats teeter deflections like blade deflections:
   
   y%HubPtMotion14%RotationVel(1,1) =     OtherState%RtHS%AngVelEL(1)
   y%HubPtMotion14%RotationVel(2,1) = -1.*OtherState%RtHS%AngVelEL(3)
   y%HubPtMotion14%RotationVel(3,1) =     OtherState%RtHS%AngVelEL(2)
      
   !...........
   ! Blade roots (AeroDyn v14):
   !...........   
   
        ! Blade root orientations should use the j instead of i system, but the current version
        ! of AeroDyn calculates forces normal and tangential to the cone of rotation
      
   DO K=1,p%NumBl
         
      y%BladeRootMotion14%Orientation(1,1,K) =     OtherState%CoordSys%i1(K,1)
      y%BladeRootMotion14%Orientation(2,1,K) =     OtherState%CoordSys%i2(K,1)
      y%BladeRootMotion14%Orientation(3,1,K) =     OtherState%CoordSys%i3(K,1)
      y%BladeRootMotion14%Orientation(1,2,K) = -1.*OtherState%CoordSys%i1(K,3)
      y%BladeRootMotion14%Orientation(2,2,K) = -1.*OtherState%CoordSys%i2(K,3)
      y%BladeRootMotion14%Orientation(3,2,K) = -1.*OtherState%CoordSys%i3(K,3)
      y%BladeRootMotion14%Orientation(1,3,K) =     OtherState%CoordSys%i1(K,2)
      y%BladeRootMotion14%Orientation(2,3,K) =     OtherState%CoordSys%i2(K,2)
      y%BladeRootMotion14%Orientation(3,3,K) =     OtherState%CoordSys%i3(K,2)
            
   END DO
   
    
   !...........
   ! Rotor furl:
   !...........   
   
      ! Rotor furl position should be rP instead of rV, but AeroDyn needs this for the HubVDue2Yaw calculation:
   
   y%RotorFurlMotion14%TranslationDisp(1,1) =     OtherState%RtHS%rV(1)
   y%RotorFurlMotion14%TranslationDisp(2,1) = -1.*OtherState%RtHS%rV(3)
   y%RotorFurlMotion14%TranslationDisp(3,1) =     OtherState%RtHS%rV(2) + p%PtfmRefzt
   
   y%RotorFurlMotion14%TranslationDisp      = y%RotorFurlMotion14%TranslationDisp - y%RotorFurlMotion14%Position   
         
        ! Rotor furl orientation (note the different order than hub and blade root!)
   
   y%RotorFurlMotion14%Orientation(1,1,1) =     OtherState%CoordSys%c1(1)
   y%RotorFurlMotion14%Orientation(2,1,1) = -1.*OtherState%CoordSys%c3(1)
   y%RotorFurlMotion14%Orientation(3,1,1) =     OtherState%CoordSys%c2(1)
   y%RotorFurlMotion14%Orientation(1,2,1) = -1.*OtherState%CoordSys%c1(3)
   y%RotorFurlMotion14%Orientation(2,2,1) =     OtherState%CoordSys%c3(3)
   y%RotorFurlMotion14%Orientation(3,2,1) = -1.*OtherState%CoordSys%c2(3)
   y%RotorFurlMotion14%Orientation(1,3,1) =     OtherState%CoordSys%c1(2)
   y%RotorFurlMotion14%Orientation(2,3,1) = -1.*OtherState%CoordSys%c3(2)
   y%RotorFurlMotion14%Orientation(3,3,1) =     OtherState%CoordSys%c2(2) 
   
      ! rotaional velocity:
   y%RotorFurlMotion14%RotationVel(1,1) =     OtherState%RtHS%AngVelER(1)
   y%RotorFurlMotion14%RotationVel(2,1) = -1.*OtherState%RtHS%AngVelER(3)
   y%RotorFurlMotion14%RotationVel(3,1) =     OtherState%RtHS%AngVelER(2)
      
   !...........
   ! Nacelle :
   !...........   
      
   y%NacelleMotion%TranslationDisp(1,1) =     OtherState%RtHS%rO(1)
   y%NacelleMotion%TranslationDisp(2,1) = -1.*OtherState%RtHS%rO(3)
   y%NacelleMotion%TranslationDisp(3,1) =     OtherState%RtHS%rO(2) + p%PtfmRefzt
               
   y%NacelleMotion%TranslationDisp      = y%NacelleMotion%TranslationDisp - y%NacelleMotion%Position   
   
      ! Nacelle orientation (note the different order than hub and blade root!)
   
   y%NacelleMotion%Orientation(1,1,1) =     OtherState%CoordSys%d1(1)
   y%NacelleMotion%Orientation(2,1,1) = -1.*OtherState%CoordSys%d3(1)
   y%NacelleMotion%Orientation(3,1,1) =     OtherState%CoordSys%d2(1)
   y%NacelleMotion%Orientation(1,2,1) = -1.*OtherState%CoordSys%d1(3)
   y%NacelleMotion%Orientation(2,2,1) =     OtherState%CoordSys%d3(3)
   y%NacelleMotion%Orientation(3,2,1) = -1.*OtherState%CoordSys%d2(3)
   y%NacelleMotion%Orientation(1,3,1) =     OtherState%CoordSys%d1(2)
   y%NacelleMotion%Orientation(2,3,1) = -1.*OtherState%CoordSys%d3(2)
   y%NacelleMotion%Orientation(3,3,1) =     OtherState%CoordSys%d2(2) 
   
   y%NacelleMotion%RotationVel(1,1)   =     OtherState%RtHS%AngVelEN(1)
   y%NacelleMotion%RotationVel(2,1)   = -1.*OtherState%RtHS%AngVelEN(3)
   y%NacelleMotion%RotationVel(3,1)   =     OtherState%RtHS%AngVelEN(2) 
      
   y%NacelleMotion%TranslationVel(1,1)  =     OtherState%RtHS%LinVelEO(1)
   y%NacelleMotion%TranslationVel(2,1)  = -1.*OtherState%RtHS%LinVelEO(3)
   y%NacelleMotion%TranslationVel(3,1)  =     OtherState%RtHS%LinVelEO(2)
      
   y%NacelleMotion%RotationAcc(   1,1)  =      AngAccEN(1) 
   y%NacelleMotion%RotationAcc(   2,1)  =  -1.*AngAccEN(3) 
   y%NacelleMotion%RotationAcc(   3,1)  =      AngAccEN(2)
   
   y%NacelleMotion%TranslationAcc(1,1)  =      LinAccEO(1)
   y%NacelleMotion%TranslationAcc(2,1)  =  -1.*LinAccEO(3)
   y%NacelleMotion%TranslationAcc(3,1)  =      LinAccEO(2)
   
   
   !...........
   ! Tower :
   !...........         
   
      ! Tower base position should be rT(0) instead of rZ, but AeroDyn needs this for  the HubVDue2Yaw calculation:
   y%TowerBaseMotion14%TranslationDisp(1,1) =     OtherState%RtHS%rZ(1)
   y%TowerBaseMotion14%TranslationDisp(2,1) = -1.*OtherState%RtHS%rZ(3)
   y%TowerBaseMotion14%TranslationDisp(3,1) =     OtherState%RtHS%rZ(2) + p%PtfmRefzt
   
   y%TowerBaseMotion14%TranslationDisp      = y%TowerBaseMotion14%TranslationDisp - y%TowerBaseMotion14%Position   
      
   y%TowerBaseMotion14%RotationVel(1,1)     =     OtherState%RtHS%AngVelEX(1)
   y%TowerBaseMotion14%RotationVel(2,1)     = -1.*OtherState%RtHS%AngVelEX(3)
   y%TowerBaseMotion14%RotationVel(3,1)     =     OtherState%RtHS%AngVelEX(2) 
   
   !...............................................................................................................................
   ! Outputs required for HydroDyn
   !...............................................................................................................................
   
   y%PlatformPtMesh%TranslationDisp(1,1) = x%QT(DOF_Sg)
   y%PlatformPtMesh%TranslationDisp(2,1) = x%QT(DOF_Sw)
   y%PlatformPtMesh%TranslationDisp(3,1) = x%QT(DOF_Hv)
   
   y%PlatformPtMesh%RotationVel(1,1)    = x%QDT(DOF_R )
   y%PlatformPtMesh%RotationVel(2,1)    = x%QDT(DOF_P )
   y%PlatformPtMesh%RotationVel(3,1)    = x%QDT(DOF_Y )
   
   y%PlatformPtMesh%TranslationVel(1,1) = x%QDT(DOF_Sg)
   y%PlatformPtMesh%TranslationVel(2,1) = x%QDT(DOF_Sw)
   y%PlatformPtMesh%TranslationVel(3,1) = x%QDT(DOF_Hv) 
   

   CALL SmllRotTrans( 'platform displacement (ED_CalcOutput)', x%QT(DOF_R ),x%QT(DOF_P ),x%QT(DOF_Y ), &
          y%PlatformPtMesh%Orientation(:,:,1), errstat=ErrStat, errmsg=ErrMsg )
      IF (ErrStat /= ErrID_None)    ErrMsg = TRIM(ErrMsg)//' (occurred at '//TRIM(Num2LStr(t))//' s)'
     !IF (ErrStat >= AbortErrLev) RETURN

   y%PlatformPtMesh%RotationAcc(1,1) = OtherState%QD2T(DOF_R )     
   y%PlatformPtMesh%RotationAcc(2,1) = OtherState%QD2T(DOF_P )     
   y%PlatformPtMesh%RotationAcc(3,1) = OtherState%QD2T(DOF_Y )     

   y%PlatformPtMesh%TranslationAcc(1,1) = OtherState%QD2T(DOF_Sg)     
   y%PlatformPtMesh%TranslationAcc(2,1) = OtherState%QD2T(DOF_Sw)    
   y%PlatformPtMesh%TranslationAcc(3,1) = OtherState%QD2T(DOF_Hv)    
      
   !...............................................................................................................................
   ! Outputs required for external tower loads
   !...............................................................................................................................
      
   DO J=1,p%TwrNodes
      y%TowerLn2Mesh%TranslationDisp(1,J) =     OtherState%RtHS%rT( 1,J) - y%TowerLn2Mesh%Position(1,J)
      y%TowerLn2Mesh%TranslationDisp(2,J) = -1.*OtherState%RtHS%rT( 3,J) - y%TowerLn2Mesh%Position(2,J)
      y%TowerLn2Mesh%TranslationDisp(3,J) =     OtherState%RtHS%rT( 2,J) - y%TowerLn2Mesh%Position(3,J) + p%PtfmRefzt
            
      y%TowerLn2Mesh%Orientation(1,1,J)   =     OtherState%CoordSys%t1(J,1)
      y%TowerLn2Mesh%Orientation(3,1,J)   =     OtherState%CoordSys%t2(J,1)
      y%TowerLn2Mesh%Orientation(2,1,J)   = -1.*OtherState%CoordSys%t3(J,1)
      y%TowerLn2Mesh%Orientation(1,2,J)   = -1.*OtherState%CoordSys%t1(J,3)
      y%TowerLn2Mesh%Orientation(3,2,J)   = -1.*OtherState%CoordSys%t2(J,3)
      y%TowerLn2Mesh%Orientation(2,2,J)   =     OtherState%CoordSys%t3(J,3)
      y%TowerLn2Mesh%Orientation(1,3,J)   =     OtherState%CoordSys%t1(J,2)
      y%TowerLn2Mesh%Orientation(3,3,J)   =     OtherState%CoordSys%t2(J,2)
      y%TowerLn2Mesh%Orientation(2,3,J)   = -1.*OtherState%CoordSys%t3(J,2)     
      
      y%TowerLn2Mesh%TranslationVel(1,J)  =     OtherState%RtHS%LinVelET(1,J)
      y%TowerLn2Mesh%TranslationVel(2,J)  = -1.*OtherState%RtHS%LinVelET(3,J)
      y%TowerLn2Mesh%TranslationVel(3,J)  =     OtherState%RtHS%LinVelET(2,J)
            
      y%TowerLn2Mesh%RotationVel(1,J)     =     OtherState%RtHS%AngVelEF(1,J)
      y%TowerLn2Mesh%RotationVel(2,J)     = -1.*OtherState%RtHS%AngVelEF(3,J)
      y%TowerLn2Mesh%RotationVel(3,J)     =     OtherState%RtHS%AngVelEF(2,J) 
            
      y%TowerLn2Mesh%TranslationAcc(1,J)  =     LinAccET(1,J)
      y%TowerLn2Mesh%TranslationAcc(2,J)  = -1.*LinAccET(3,J)
      y%TowerLn2Mesh%TranslationAcc(3,J)  =     LinAccET(2,J)
            
      y%TowerLn2Mesh%RotationAcc(1,J)     =     AngAccEF(1,J)
      y%TowerLn2Mesh%RotationAcc(2,J)     = -1.*AngAccEF(3,J)
      y%TowerLn2Mesh%RotationAcc(3,J)     =     AngAccEF(2,J) 
      
   END DO
               
   
   ! p%TwrNodes+1 is the tower top:
   J = p%TwrNodes+1
   
   y%TowerLn2Mesh%TranslationDisp(1,J) =     OtherState%RtHS%rO(1) - y%TowerLn2Mesh%Position(1,J)
   y%TowerLn2Mesh%TranslationDisp(2,J) = -1.*OtherState%RtHS%rO(3) - y%TowerLn2Mesh%Position(2,J)
   y%TowerLn2Mesh%TranslationDisp(3,J) =     OtherState%RtHS%rO(2) - y%TowerLn2Mesh%Position(3,J) + p%PtfmRefzt
   
   y%TowerLn2Mesh%Orientation(1,1,J)   =     OtherState%CoordSys%b1(1)
   y%TowerLn2Mesh%Orientation(3,1,J)   =     OtherState%CoordSys%b2(1)
   y%TowerLn2Mesh%Orientation(2,1,J)   = -1.*OtherState%CoordSys%b3(1)
   y%TowerLn2Mesh%Orientation(1,2,J)   = -1.*OtherState%CoordSys%b1(3)
   y%TowerLn2Mesh%Orientation(3,2,J)   = -1.*OtherState%CoordSys%b2(3)
   y%TowerLn2Mesh%Orientation(2,2,J)   =     OtherState%CoordSys%b3(3)
   y%TowerLn2Mesh%Orientation(1,3,J)   =     OtherState%CoordSys%b1(2)
   y%TowerLn2Mesh%Orientation(3,3,J)   =     OtherState%CoordSys%b2(2)
   y%TowerLn2Mesh%Orientation(2,3,J)   = -1.*OtherState%CoordSys%b3(2)         
          
   y%TowerLn2Mesh%TranslationVel(1,J)  =     OtherState%RtHS%LinVelEO(1)         
   y%TowerLn2Mesh%TranslationVel(2,J)  = -1.*OtherState%RtHS%LinVelEO(3)   
   y%TowerLn2Mesh%TranslationVel(3,J)  =     OtherState%RtHS%LinVelEO(2)        

   y%TowerLn2Mesh%RotationVel(1,J)     =     OtherState%RtHS%AngVelEB(1)
   y%TowerLn2Mesh%RotationVel(2,J)     = -1.*OtherState%RtHS%AngVelEB(3)
   y%TowerLn2Mesh%RotationVel(3,J)     =     OtherState%RtHS%AngVelEB(2) 

   y%TowerLn2Mesh%TranslationAcc(1,J)  =     LinAccEO(1)
   y%TowerLn2Mesh%TranslationAcc(2,J)  = -1.*LinAccEO(3)
   y%TowerLn2Mesh%TranslationAcc(3,J)  =     LinAccEO(2)
   
   y%TowerLn2Mesh%RotationAcc(1,J)     =     AngAccEB(1)
   y%TowerLn2Mesh%RotationAcc(2,J)     = -1.*AngAccEB(3)
   y%TowerLn2Mesh%RotationAcc(3,J)     =     AngAccEB(2) 

   
   ! p%TwrNodes+2 is the tower base:
   J = p%TwrNodes+2

   y%TowerLn2Mesh%TranslationDisp(1,J) =     OtherState%RtHS%rZ(1) + OtherState%RtHS%rZT0(1) - y%TowerLn2Mesh%Position(1,J)
   y%TowerLn2Mesh%TranslationDisp(2,J) = -1.*OtherState%RtHS%rZ(3) - OtherState%RtHS%rZT0(3) - y%TowerLn2Mesh%Position(2,J)
   y%TowerLn2Mesh%TranslationDisp(3,J) =     OtherState%RtHS%rZ(2) + OtherState%RtHS%rZT0(2) - y%TowerLn2Mesh%Position(3,J) + p%PtfmRefzt
      
   y%TowerLn2Mesh%Orientation(1,1,J)   =     OtherState%CoordSys%a1(1)
   y%TowerLn2Mesh%Orientation(3,1,J)   =     OtherState%CoordSys%a2(1)
   y%TowerLn2Mesh%Orientation(2,1,J)   = -1.*OtherState%CoordSys%a3(1)
   y%TowerLn2Mesh%Orientation(1,2,J)   = -1.*OtherState%CoordSys%a1(3)
   y%TowerLn2Mesh%Orientation(3,2,J)   = -1.*OtherState%CoordSys%a2(3)
   y%TowerLn2Mesh%Orientation(2,2,J)   =     OtherState%CoordSys%a3(3)
   y%TowerLn2Mesh%Orientation(1,3,J)   =     OtherState%CoordSys%a1(2)
   y%TowerLn2Mesh%Orientation(3,3,J)   =     OtherState%CoordSys%a2(2)
   y%TowerLn2Mesh%Orientation(2,3,J)   = -1.*OtherState%CoordSys%a3(2)

   y%TowerLn2Mesh%TranslationVel(1,J)  =     OtherState%RtHS%LinVelET(1,0)       
   y%TowerLn2Mesh%TranslationVel(2,J)  = -1.*OtherState%RtHS%LinVelET(3,0) 
   y%TowerLn2Mesh%TranslationVel(3,J)  =     OtherState%RtHS%LinVelET(2,0)   
   
   y%TowerLn2Mesh%RotationVel(1,J)     =     OtherState%RtHS%AngVelEF(1,0)
   y%TowerLn2Mesh%RotationVel(2,J)     = -1.*OtherState%RtHS%AngVelEF(3,0)
   y%TowerLn2Mesh%RotationVel(3,J)     =     OtherState%RtHS%AngVelEF(2,0) 
   
   y%TowerLn2Mesh%TranslationAcc(1,J)  =     LinAccET(1,0)
   y%TowerLn2Mesh%TranslationAcc(2,J)  = -1.*LinAccET(3,0)
   y%TowerLn2Mesh%TranslationAcc(3,J)  =     LinAccET(2,0)
   
   y%TowerLn2Mesh%RotationAcc(1,J)     =     AngAccEF(1,0)
   y%TowerLn2Mesh%RotationAcc(2,J)     = -1.*AngAccEF(3,0)
   y%TowerLn2Mesh%RotationAcc(3,J)     =     AngAccEF(2,0) 
   
   !...............................................................................................................................
   ! Outputs required for ServoDyn
   !...............................................................................................................................
   
   y%Yaw      = x%QT( DOF_Yaw)
   y%YawRate  = x%QDT(DOF_Yaw)
   y%YawAngle = x%QT( DOF_Yaw) + x%QT(DOF_Y)  !crude approximation for yaw error... (without subtracting it from the wind direction)   
   y%BlPitch  = u%BlPitchCom !OtherState%BlPitch
   y%LSS_Spd  = x%QDT(DOF_GeAz)
   y%HSS_Spd  = ABS(p%GBRatio)*x%QDT(DOF_GeAz)
   y%RotSpeed = x%QDT(DOF_GeAz) + x%QDT(DOF_DrTr)
   
   IF ( t > 0.0_DbKi  )  THEN

      ! Calculate tower-top acceleration (fore-aft mode only) in the tower-top system:

      LinAccEO = OtherState%RtHS%LinAccEOt
      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)
         LinAccEO = LinAccEO + OtherState%RtHS%PLinVelEO(p%DOFs%PTE(I),0,:)*OtherState%QD2T(p%DOFs%PTE(I))
      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

      y%TwrAccel = DOT_PRODUCT( LinAccEO, OtherState%CoordSys%b1 )
   ELSE
      y%TwrAccel = 0
   END IF      
   
   !Control outputs for Bladed DLL:
   y%RotPwr    = OtherState%AllOuts(  RotPwr)*1000.
   DO K=1,p%NumBl
      y%RootMxc(K)   = OtherState%AllOuts( RootMxc(K) )*1000.
      y%RootMyc(K)   = OtherState%AllOuts( RootMyc(K) )*1000.
   END DO
   y%YawBrTAxp = OtherState%AllOuts( YawBrTAxp)
   y%YawBrTAyp = OtherState%AllOuts( YawBrTAyp)
   !y%LSSTipPxa = OtherState%AllOuts( LSSTipPxa)*D2R ! bjj: did this above already

   y%LSSTipMxa = OtherState%AllOuts(LSShftMxa)*1000.
   y%LSSTipMya = OtherState%AllOuts(LSSTipMya)*1000.                ! Rotating hub My (GL co-ords) (Nm)
   y%LSSTipMza = OtherState%AllOuts(LSSTipMza)*1000.                ! Rotating hub Mz (GL co-ords) (Nm)
   y%LSSTipMys = OtherState%AllOuts(LSSTipMys)*1000.                ! Fixed hub My (GL co-ords) (Nm)
   y%LSSTipMzs = OtherState%AllOuts(LSSTipMzs)*1000.                ! Fixed hub Mz (GL co-ords) (Nm)
   y%YawBrMyn  = OtherState%AllOuts( YawBrMyn)*1000.                ! Yaw bearing My (GL co-ords) (Nm) !tower accel
   y%YawBrMzn  = OtherState%AllOuts( YawBrMzn)*1000.                ! Yaw bearing Mz (GL co-ords) (Nm)
   y%NcIMURAxs = OtherState%AllOuts(NcIMURAxs)*D2R                  ! Nacelle roll    acceleration (rad/s^2) -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system
   y%NcIMURAys = OtherState%AllOuts(NcIMURAys)*D2R                  ! Nacelle nodding acceleration (rad/s^2)
   y%NcIMURAzs = OtherState%AllOuts(NcIMURAzs)*D2R                  ! Nacelle yaw     acceleration (rad/s^2) -- this is in the shaft (tilted) coordinate system, instead of the nacelle (nontilted) coordinate system
   
               
   RETURN
   

END SUBROUTINE ED_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(ED_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(ED_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(ED_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(ED_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at t
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! LOCAL variables
   LOGICAL, PARAMETER           :: UpdateValues  = .TRUE.      ! determines if the OtherState values need to be updated
      
   INTEGER(IntKi)                         :: I                 ! Loops through some or all of the DOFs.
   INTEGER(IntKi)                         :: ErrStat2          ! The error status code
   CHARACTER(ErrMsgLen)                   :: ErrMsg2           ! The error message, if an error occurred
   CHARACTER(*), PARAMETER                :: RoutineName = 'ED_CalcContStateDeriv'
   
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !OtherState%HSSBrTrqC = SIGN( u%HSSBrTrqC, x%QDT(DOF_GeAz) ) !need correct value of x%QDT(DOF_GeAz) here

         ! Compute the first time derivatives of the continuous states here:

   ! See if the values stored in OtherState%RtHS and OtherState%CoordSys need to be updated:
   IF ( UpdateValues ) THEN       
      
       !OtherState%BlPitch = u%BlPitchCom
       
         ! set the coordinate system variables:
      CALL SetCoordSy( t, OtherState%CoordSys, OtherState%RtHS, u%BlPitchCom, p, x, ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN
   
      CALL CalculatePositions(        p, x, OtherState%CoordSys,    OtherState%RtHS ) ! calculate positions
      CALL CalculateAngularPosVelPAcc(p, x, OtherState%CoordSys,    OtherState%RtHS ) ! calculate angular positions, velocities, and partial accelerations, including partial angular quantities
      CALL CalculateLinearVelPAcc(    p, x, OtherState%CoordSys,    OtherState%RtHS ) ! calculate linear velocities and partial accelerations
      CALL CalculateForcesMoments(    p, x, OtherState%CoordSys, u, OtherState%RtHS ) ! calculate the forces and moments (requires AeroBladeForces and AeroBladeMoments)            
      
   END IF
      
   !.....................................
   !  TeetMom,  RFrlMom, TFrlMom (possibly from user routines)
   ! bjj: we will want to revisit these routines:
   !.....................................
   
      ! Compute the moments from teeter springs and dampers, rotor-furl springs and dampers, tail-furl springs and dampers

   CALL Teeter  ( t, p, OtherState%RtHS%TeetAng, OtherState%RtHS%TeetAngVel, OtherState%RtHS%TeetMom ) ! Compute moment from teeter     springs and dampers, TeetMom; NOTE: TeetMom will be zero for a 3-blader since TeetAng = TeetAngVel = 0
   CALL RFurling( t, p, x%QT(DOF_RFrl),          x%QDT(DOF_RFrl),            OtherState%RtHS%RFrlMom ) ! Compute moment from rotor-furl springs and dampers, RFrlMom
   CALL TFurling( t, p, x%QT(DOF_TFrl),          x%QDT(DOF_TFrl),            OtherState%RtHS%TFrlMom ) ! Compute moment from tail-furl  springs and dampers, TFrlMom
   
   !bjj: note OtherState%RtHS%GBoxEffFac needed in OtherState only to fix HSSBrTrq (and used in FillAugMat)
   OtherState%RtHS%GBoxEffFac  = p%GBoxEff**OtherState%SgnPrvLSTQ      ! = GBoxEff if SgnPrvLSTQ = 1 OR 1/GBoxEff if SgnPrvLSTQ = -1
   
   CALL FillAugMat( p, x, OtherState%CoordSys, u, OtherState%HSSBrTrq, OtherState%RtHS, OtherState%AugMat )
   

   ! Invert the matrix to solve for the accelerations. The accelerations are returned by Gauss() in the first NActvDOF elements
   !   of the solution vector, SolnVec(). These are transfered to the proper index locations of the acceleration vector QD2T()
   !   using the vector subscript array SrtPS(), after Gauss() has been called:
   ! NOTE: QD2T( SrtPS(1:NActvDOF) ) cannot be sent directly because arrays sections with vector subscripts must not be used 
   !   in INTENT(OUT) arguments.

   IF ( p%DOFs%NActvDOF > 0 ) THEN
      OtherState%AugMat_factor = OtherState%AugMat( p%DOFs%SrtPS( 1:p%DOFs%NActvDOF ), p%DOFs%SrtPSNAUG(1:p%DOFs%NActvDOF) )
      OtherState%SolnVec       = OtherState%AugMat( p%DOFs%SrtPS( 1:p%DOFs%NActvDOF ), p%DOFs%SrtPSNAUG(1+p%DOFs%NActvDOF) )
   

      CALL LAPACK_getrf( M=p%DOFs%NActvDOF, N=p%DOFs%NActvDOF, A=OtherState%AugMat_factor, IPIV=OtherState%AugMat_pivot, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
         IF ( ErrStat >= AbortErrLev ) RETURN
      
      CALL LAPACK_getrs( TRANS='N',N=p%DOFs%NActvDOF, A=OtherState%AugMat_factor,IPIV=OtherState%AugMat_pivot, B=OtherState%SolnVec, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
   
      !CALL GaussElim( OtherState%AugMat( p%DOFs%SrtPS(    1: p%DOFs%NActvDOF   ),     &
      !                                   p%DOFs%SrtPSNAUG(1:(p%DOFs%NActvDOF+1)) ),   &
      !                                   p%DOFs%NActvDOF,  SolnVec, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
   END IF
   

   !bjj: because the deriv is INTENT(OUT), this is reallocated each time:
IF (.NOT. ALLOCATED(dxdt%qt) ) THEN
   CALL AllocAry( dxdt%qt,  SIZE(x%qt),  'dxdt%qt',  ErrStat2, ErrMsg2 )
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   IF ( ErrStat >= AbortErrLev ) RETURN
END IF

IF (.NOT. ALLOCATED(dxdt%qdt) ) THEN
   CALL AllocAry( dxdt%qdt, SIZE(x%qdt), 'dxdt%qdt', ErrStat2, ErrMsg2 )
   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   IF ( ErrStat >= AbortErrLev ) RETURN
END IF

   dxdt%QT = x%QDT
   
   dxdt%QDT = 0.0
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
      dxdt%QDT(p%DOFs%SrtPS(I)) = OtherState%SolnVec(I)    ! dxdt%QDT = OtherState%QD2T
   ENDDO             ! I - All active (enabled) DOFs

   OtherState%QD2T = dxdt%QDT
      
   
      ! Let's calculate the sign (+/-1) of the low-speed shaft torque for this time step and store it in SgnPrvLSTQ.
      !  This will be used during the next call to RtHS (bjj: currently violates framework, but DOE wants a hack for HSS brake).
      ! need OtherState%QD2T set before calling this
      
   !OtherState%SgnPrvLSTQ = SignLSSTrq(p, OtherState)
   
   
END SUBROUTINE ED_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Tight coupling routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),               INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
      TYPE(ED_InputType),           INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ED_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                  !   Output: Discrete states at t + Interval
      TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ED_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE ED_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(ED_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(ED_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(ED_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ED_ConstraintStateType), INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using
                                                                  !     the input values described above
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Solve for the constraint states here:

      z_residual%DummyConstrState = 0.

END SUBROUTINE ED_CalcConstrStateResidual
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_ReadInput( InputFileName, MeshFile, InputFileData, ReadAdmVals, BD4Blades, Default_DT, OutFileRoot, ErrStat, ErrMsg )
! This subroutine reads the input file and stores all the data in the ED_InputFile structure.
! It does not perform data validation.
!..................................................................................................................................

      ! Passed variables
   REAL(DbKi),           INTENT(IN)       :: Default_DT     ! The default DT (from glue code)

   CHARACTER(*), INTENT(IN)               :: InputFileName  ! Name of the input file
   CHARACTER(*), INTENT(IN)               :: MeshFile       ! File that contains the blade mesh information (AeroDyn input file for now) -- later this info will be defined in one of the ED input files.
   CHARACTER(*), INTENT(IN)               :: OutFileRoot    ! The rootname of all the output files written by this routine.

   TYPE(ED_InputFile),   INTENT(OUT)      :: InputFileData  ! Data stored in the module's input file

   INTEGER(IntKi),       INTENT(OUT)      :: ErrStat        ! The error status code
   LOGICAL,              INTENT(IN)       :: ReadAdmVals    ! Determines if we should read the Adams-only values
   LOGICAL,              INTENT(IN)       :: BD4Blades      ! Determines if we should read the blade values (true=don't read this file; use BeamDyn for blades instead)
   CHARACTER(*),         INTENT(OUT)      :: ErrMsg         ! The error message, if an error occurred

      ! local variables

!   INTEGER(IntKi)                         :: UnIn           ! Unit number for the input file
   INTEGER(IntKi)                         :: UnEcho         ! Unit number for the echo file
   INTEGER(IntKi)                         :: ErrStat2       ! The error status code
   CHARACTER(ErrMsgLen)                   :: ErrMsg2        ! The error message, if an error occurred

   CHARACTER(1024)                        :: BldFile(MaxBl) ! File that contains the blade information (specified in the primary input file)
   CHARACTER(1024)                        :: FurlFile       ! File that contains the furl information (specified in the primary input file)
   CHARACTER(1024)                        :: TwrFile        ! File that contains the tower information (specified in the primary input file)

      ! initialize values:

   ErrStat = ErrID_None
   ErrMsg  = ''

   InputFileData%DT = Default_DT  ! the glue code's suggested DT for the module (may be overwritten in ReadPrimaryFile())

      ! get the primary/platform input-file data
      ! sets UnEcho, BldFile, FurlFile, TwrFile
   
   CALL ReadPrimaryFile( InputFileName, InputFileData, BldFile, FurlFile, TwrFile, OutFileRoot, UnEcho, ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! get the furling input-file data
   InputFileData%Furling = .FALSE.              ! Furling is not supported in this version of ElastoDyn

   IF ( InputFileData%Furling )  THEN
      CALL ReadFurlFile( FurlFile, InputFileData, UnEcho, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
   ELSE   ! initialize all of the data that would be read by ReadFurlFile()
      InputFileData%RFrlDOF   = .FALSE.
      InputFileData%TFrlDOF   = .FALSE.
      InputFileData%RotFurl   = 0.0_ReKi  ! Radians
      InputFileData%TailFurl  = 0.0_ReKi
      InputFileData%Yaw2Shft  = 0.0
      InputFileData%ShftSkew  = 0.0
      InputFileData%RFrlCMxn  = 0.0
      InputFileData%RFrlCMyn  = 0.0
      InputFileData%RFrlCMzn  = 0.0
      InputFileData%BoomCMxn  = 0.0
      InputFileData%BoomCMyn  = 0.0
      InputFileData%BoomCMzn  = 0.0
      InputFileData%TFinCMxn  = 0.0
      InputFileData%TFinCMyn  = 0.0
      InputFileData%TFinCMzn  = 0.0
      InputFileData%TFinCPxn  = 0.0
      InputFileData%TFinCPyn  = 0.0
      InputFileData%TFinCPzn  = 0.0
      InputFileData%TFinSkew  = 0.0
      InputFileData%TFinTilt  = 0.0
      InputFileData%TFinBank  = 0.0
      InputFileData%RFrlPntxn = 0.0
      InputFileData%RFrlPntyn = 0.0
      InputFileData%RFrlPntzn = 0.0
      InputFileData%RFrlSkew  = 0.0
      InputFileData%RFrlTilt  = 0.0
      InputFileData%TFrlPntxn = 0.0
      InputFileData%TFrlPntyn = 0.0
      InputFileData%TFrlPntzn = 0.0
      InputFileData%TFrlSkew  = 0.0
      InputFileData%TFrlTilt  = 0.0
      InputFileData%RFrlMass  = 0.0
      InputFileData%BoomMass  = 0.0
      InputFileData%TFinMass  = 0.0
      InputFileData%RFrlIner  = 0.0
      InputFileData%TFrlIner  = 0.0
      InputFileData%RFrlMod   = 0
      InputFileData%RFrlSpr   = 0.0
      InputFileData%RFrlDmp   = 0.0
      InputFileData%RFrlCDmp  = 0.0
      InputFileData%RFrlUSSP  = 0.0
      InputFileData%RFrlDSSP  = 0.0
      InputFileData%RFrlUSSpr = 0.0
      InputFileData%RFrlDSSpr = 0.0
      InputFileData%RFrlUSDP  = 0.0
      InputFileData%RFrlDSDP  = 0.0
      InputFileData%RFrlUSDmp = 0.0
      InputFileData%RFrlDSDmp = 0.0
      InputFileData%TFrlMod   = 0
      InputFileData%TFrlSpr   = 0.0
      InputFileData%TFrlDmp   = 0.0
      InputFileData%TFrlCDmp  = 0.0
      InputFileData%TFrlUSSP  = 0.0
      InputFileData%TFrlDSSP  = 0.0
      InputFileData%TFrlUSSpr = 0.0
      InputFileData%TFrlDSSpr = 0.0
      InputFileData%TFrlUSDP  = 0.0
      InputFileData%TFrlDSDP  = 0.0
      InputFileData%TFrlUSDmp = 0.0
      InputFileData%TFrlDSDmp = 0.0
   END IF


      ! get the blade input-file data (from blade and mesh files)
   IF (.NOT. BD4Blades) THEN
      CALL ReadBladeInputs ( BldFile, MeshFile, ReadAdmVals, InputFileData, UnEcho, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

      ! get the tower input-file data

   CALL ReadTowerFile( TwrFile, InputFileData, ReadAdmVals, UnEcho,  ErrStat2, ErrMsg2 )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! close the echo file (if opened)

   IF ( UnEcho > 0 ) CLOSE( UnEcho )


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_ReadInput:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            IF ( UnEcho > 0 ) CLOSE( UnEcho )
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE ED_ReadInput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_ValidateInput( InputFileData, BD4Blades, ErrStat, ErrMsg )
! This subroutine validates the input file data
!..................................................................................................................................

   TYPE(ED_InputFile),       INTENT(IN)       :: InputFileData       ! Data stored in the module's input file
   LOGICAL,                  INTENT(IN)       :: BD4Blades           ! Determines if we should validate the blade values (true=don't validate; use BeamDyn for blades instead)
   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat             ! The error status code
   CHARACTER(*),             INTENT(OUT)      :: ErrMsg              ! The error message, if an error occurred

      ! Local variables:
!   INTEGER(IntKi)                             :: I                   ! Loop counter
   INTEGER(IntKi)                             :: K                   ! Blade number
   INTEGER(IntKi)                             :: ErrStat2            ! Temporary error ID
!   LOGICAL                                    :: ReadAdmVals         ! determines if an Adams model will be created (do we read/check all the inputs?)
!   LOGICAL                                    :: ReadFile            ! determines if an input file for a blade is the same as the file for the previous blade
   CHARACTER(ErrMsgLen)                       :: ErrMsg2             ! Temporary message describing error


      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ''



      ! validate the primary input data
   CALL ValidatePrimaryData( InputFileData, BD4Blades, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )


      ! validate the furling input data
   CALL ValidateFurlData ( InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )


      ! validate the blade input data
   IF (.NOT. BD4Blades) THEN
      DO K = 1,InputFileData%NumBl
         CALL ValidateBladeData ( InputFileData%InpBl(K), ErrStat2, ErrMsg2 )
            CALL CheckError( ErrStat2, ' Errors in blade '//TRIM(Num2LStr(K))//' input data: '//NewLine//TRIM(ErrMsg2) )
      END DO
      
      !bjj: validate blade discretization, too:
   END IF
   

      ! validate the tower input data
   CALL ValidateTowerData ( InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      
      
      ! validate the Output parameters:
  ! CALL ChckOutLst( InputFileData%OutList, p, ErrStat, ErrMsg )



CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_ValidateInput:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

      END IF


   END SUBROUTINE CheckError
END SUBROUTINE ED_ValidateInput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_SetParameters( InputFileData, p, ErrStat, ErrMsg )
! This subroutine sets the parameters, based on the data stored in InputFileData
!..................................................................................................................................

   TYPE(ED_InputFile),       INTENT(IN)       :: InputFileData  ! Data stored in the module's input file
   TYPE(ED_ParameterType),   INTENT(INOUT)    :: p              ! The module's parameter data
   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat        ! The error status code
   CHARACTER(*),             INTENT(OUT)      :: ErrMsg         ! The error message, if an error occurred

      ! Local variables
!   INTEGER(IntKi)                             :: K              ! Loop counter (for blades)
   INTEGER(IntKi)                             :: ErrStat2       ! Temporary error ID
   CHARACTER(ErrMsgLen)                       :: ErrMsg2        ! Temporary message describing error

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ''



      ! Set parameters from primary input file
   CALL SetPrimaryParameters( p, InputFileData, ErrStat2, ErrMsg2  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   p%DT24 = p%DT/24.0_DbKi    ! Time-step parameter needed for Solver().


      ! Set furling parameters
   CALL SetFurlParameters( p, InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Set blade parameters
   CALL SetBladeParameters( p, InputFileData%InpBl, InputFileData%InpBlMesh, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Set tower parameters
   CALL SetTowerParameters( p, InputFileData, ErrStat2, ErrMsg2 ) ! It requires p%TwrFlexL, and p%TwrNodes to be set first.
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Set the remaining (calculuated) parameters (basically the former Coeff routine)
   CALL SetOtherParameters( p, InputFileData, ErrStat2, ErrMsg2 ) ! requires MANY things to be set first
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   
CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_SetParameters:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE ED_SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_DOFparameters( InputFileData, p, ErrStat, ErrMsg )
! This subroutine initializes the ActiveDOF data type as well as the variables related to DOFs, including p%NAug and p%NDOF.
! It assumes that p%NumBl is set.
!..................................................................................................................................

!   TYPE(ED_ActiveDOFs),      INTENT(INOUT)    :: DOFs           ! ActiveDOF data
   TYPE(ED_InputFile),       INTENT(IN)       :: InputFileData  ! Data stored in the module's input file
   TYPE(ED_ParameterType),   INTENT(INOUT)    :: p              ! The module's parameter data
   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat        ! The error status code
   CHARACTER(*),             INTENT(OUT)      :: ErrMsg         ! The error message, if an error occurred

      ! Local variables
   INTEGER(IntKi)                             :: K              ! Loop counter (for blades)

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ''


   IF ( p%NumBl == 2 )  THEN
      p%NDOF = 22
   ELSE
      p%NDOF = 24
   ENDIF

   p%NAug = p%NDOF + 1

   ! ...........................................................................................................................
   ! allocate and set DOF_Flag and DOF_Desc
   ! ...........................................................................................................................
   CALL AllocAry( p%DOF_Flag, p%NDOF,   'DOF_Flag',  ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%DOF_Desc, p%NDOF,   'DOF_Desc',  ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


   IF ( p%NumBl == 2 )  THEN ! the 3rd blade overwrites the DOF_Teet position of the array, so don't use an "ELSE" for this statement
      p%DOF_Flag(DOF_Teet) = InputFileData%TeetDOF
      p%DOF_Desc(DOF_Teet) = 'Hub teetering DOF (internal DOF index = DOF_Teet)'
   END IF !


   DO K = 1,p%NumBl
      p%DOF_Flag( DOF_BF(K,1) ) = InputFileData%FlapDOF1
      p%DOF_Desc( DOF_BF(K,1) ) = '1st flapwise bending-mode DOF of blade '//TRIM(Num2LStr( K ))// &
                                  ' (internal DOF index = DOF_BF('         //TRIM(Num2LStr( K ))//',1))'

      p%DOF_Flag( DOF_BE(K,1) ) = InputFileData%EdgeDOF
      p%DOF_Desc( DOF_BE(K,1) ) = '1st edgewise bending-mode DOF of blade '//TRIM(Num2LStr( K ))// &
                                  ' (internal DOF index = DOF_BE('         //TRIM(Num2LStr( K ))//',1))'

      p%DOF_Flag( DOF_BF(K,2) ) = InputFileData%FlapDOF2
      p%DOF_Desc( DOF_BF(K,2) ) = '2nd flapwise bending-mode DOF of blade '//TRIM(Num2LStr( K ))// &
                                  ' (internal DOF index = DOF_BF('         //TRIM(Num2LStr( K ))//',2))'
   ENDDO          ! K - All blades

   p%DOF_Flag(DOF_DrTr) = InputFileData%DrTrDOF
   p%DOF_Desc(DOF_DrTr) = 'Drivetrain rotational-flexibility DOF (internal DOF index = DOF_DrTr)'
   p%DOF_Flag(DOF_GeAz) = InputFileData%GenDOF
   p%DOF_Desc(DOF_GeAz) = 'Variable speed generator DOF (internal DOF index = DOF_GeAz)'
   p%DOF_Flag(DOF_RFrl) = InputFileData%RFrlDOF
   p%DOF_Desc(DOF_RFrl) = 'Rotor-furl DOF (internal DOF index = DOF_RFrl)'
   p%DOF_Flag(DOF_TFrl) = InputFileData%TFrlDOF
   p%DOF_Desc(DOF_TFrl) = 'Tail-furl DOF (internal DOF index = DOF_TFrl)'
   p%DOF_Flag(DOF_Yaw ) = InputFileData%YawDOF
   p%DOF_Desc(DOF_Yaw ) = 'Nacelle yaw DOF (internal DOF index = DOF_Yaw)'
   p%DOF_Flag(DOF_TFA1) = InputFileData%TwFADOF1
   p%DOF_Desc(DOF_TFA1) = '1st tower fore-aft bending mode DOF (internal DOF index = DOF_TFA1)'
   p%DOF_Flag(DOF_TSS1) = InputFileData%TwSSDOF1
   p%DOF_Desc(DOF_TSS1) = '1st tower side-to-side bending mode DOF (internal DOF index = DOF_TSS1)'
   p%DOF_Flag(DOF_TFA2) = InputFileData%TwFADOF2
   p%DOF_Desc(DOF_TFA2) = '2nd tower fore-aft bending mode DOF (internal DOF index = DOF_TFA2)'
   p%DOF_Flag(DOF_TSS2) = InputFileData%TwSSDOF2
   p%DOF_Desc(DOF_TSS2) = '2nd tower side-to-side bending mode DOF (internal DOF index = DOF_TSS2)'
   p%DOF_Flag(DOF_Sg  ) = InputFileData%PtfmSgDOF
   p%DOF_Desc(DOF_Sg  ) = 'Platform horizontal surge translation DOF (internal DOF index = DOF_Sg)'
   p%DOF_Flag(DOF_Sw  ) = InputFileData%PtfmSwDOF
   p%DOF_Desc(DOF_Sw  ) = 'Platform horizontal sway translation DOF (internal DOF index = DOF_Sw)'
   p%DOF_Flag(DOF_Hv  ) = InputFileData%PtfmHvDOF
   p%DOF_Desc(DOF_Hv  ) = 'Platform vertical heave translation DOF (internal DOF index = DOF_Hv)'
   p%DOF_Flag(DOF_R   ) = InputFileData%PtfmRDOF
   p%DOF_Desc(DOF_R   ) = 'Platform roll tilt rotation DOF (internal DOF index = DOF_R)'
   p%DOF_Flag(DOF_P   ) = InputFileData%PtfmPDOF
   p%DOF_Desc(DOF_P   ) = 'Platform pitch tilt rotation DOF (internal DOF index = DOF_P)'
   p%DOF_Flag(DOF_Y   ) = InputFileData%PtfmYDOF
   p%DOF_Desc(DOF_Y   ) = 'Platform yaw rotation DOF (internal DOF index = DOF_Y)'

   ! ...........................................................................................................................
   ! allocate the arrays stored in the p%DOFs structure:
   ! ...........................................................................................................................

      ! BJJ: note that this method will cause an error if allocating data that has already been allocated...

   ALLOCATE ( p%DOFs%NPSBE(p%NumBl), p%DOFs%NPSE(p%NumBl),  STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs NPSBE and NPSE arrays.' )
      RETURN
   ENDIF


   ALLOCATE ( p%DOFs%PCE(p%NDOF), p%DOFs%PDE(p%NDOF), p%DOFs%PIE(p%NDOF), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs PCE, PDE, and PIE arrays.' )
      RETURN
   ENDIF


   ALLOCATE (  p%DOFs%PTTE(p%NDOF), p%DOFs%PTE(p%NDOF), p%DOFs%PS(p%NDOF), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs PTTE, PTE, and PS arrays.' )
      RETURN
   ENDIF


   ALLOCATE ( p%DOFs%PUE(p%NDOF), p%DOFs%PYE(p%NDOF),  STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs PUE and PYE arrays.' )
      RETURN
   ENDIF


!bjj was   ALLOCATE ( p%DOFs%PSBE(p%NumBl,3), p%DOFs%PSE(p%NumBl,p%NDOF),  STAT=ErrStat )
   ALLOCATE ( p%DOFs%PSBE(p%NumBl,(NumBE+NumBF)), p%DOFs%PSE(p%NumBl,p%NDOF),  STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs PSBE and PSE arrays.' )
      RETURN
   ENDIF


   ALLOCATE ( p%DOFs%SrtPS(p%NDOF), p%DOFs%SrtPSNAUG(p%NAug),  p%DOFs%Diag(p%NDOF), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveAOFs SrtPS, SrtPSNAUG, and Diag arrays.' )
      RETURN
   ENDIF


   !...............................................................................................................................
   ! Allocate and Initialize arrays for DOFS that contribute to the angular velocity of the hub and blade elements
   !...............................................................................................................................
   ! Define arrays of DOF indices (pointers) that contribute to the angular
   !   velocities of each rigid body of the wind turbine in the inertia frame:
   ! NOTE: We must include ALL of the appropriate DOF indices in these arrays,
   !       not just the indices of the enabled DOFs, since disabling a DOF only
   !       implies that each DOF acceleration is zero--it does not imply
   !       that each DOF velocity is zero (for example, consider disabled
   !       generator DOF, which still spins at constant speed).



   IF ( p%NumBl == 2 )  THEN ! 2-blader
      p%NPH = 12                         ! Number of DOFs that contribute to the angular velocity of the hub            (body H) in the inertia frame.
      p%NPM = 15                         ! Number of DOFs that contribute to the angular velocity of the blade elements (body M) in the inertia frame.
   ELSE                    ! 3-blader
      p%NPH = 11                         ! Number of DOFs that contribute to the angular velocity of the hub            (body H) in the inertia frame.
      p%NPM = 14                         ! Number of DOFs that contribute to the angular velocity of the blade elements (body M) in the inertia frame.
   ENDIF


   ALLOCATE ( p%PH(p%NPH),  p%PM(p%NumBl,p%NPM), STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine( ErrID_Fatal, ' Could not allocate memory for the ActiveDOFs PH and PM arrays.' )
      RETURN
   ENDIF

      ! Array of DOF indices (pointers) that contribute to the angular velocity of the hub (body H) in the inertia frame:
   p%PH(1:11) = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz, DOF_DrTr /)

   IF ( p%NumBl == 2 )  THEN ! 2-blader (add DOF_Teet to the arrays)

      p%PH(12) = DOF_Teet

         ! Array of DOF indices (pointers) that contribute to the angular velocity of the blade elements (body M) in the inertia frame:
      DO K = 1,p%NumBl ! Loop through all blades
         p%PM(K,:) = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz, DOF_DrTr, &
                        DOF_Teet,  DOF_BF(K,1) , DOF_BE(K,1)    , DOF_BF(K,2)          /)
      ENDDO          ! K - All blades

   ELSE  ! 3-blader

         ! Array of DOF indices (pointers) that contribute to the angular velocity of the blade elements (body M) in the inertia frame:
      DO K = 1,p%NumBl ! Loop through all blades
         p%PM(K,:) = (/ DOF_R, DOF_P, DOF_Y, DOF_TFA1, DOF_TSS1, DOF_TFA2, DOF_TSS2, DOF_Yaw, DOF_RFrl, DOF_GeAz, DOF_DrTr, &
                                   DOF_BF(K,1) , DOF_BE(K,1)    , DOF_BF(K,2)         /)
      ENDDO          ! K - All blades

   ENDIF


   !...............................................................................................................................
   ! Calculate the number of active (enabled) DOFs in the model, p%DOFs%NActvDOF:
   !...............................................................................................................................
   CALL SetEnabledDOFIndexArrays( p )

   RETURN


CONTAINS
   !............................................................................................................................
   SUBROUTINE ExitThisRoutine(ErrID,Msg)
   ! This subroutine cleans up all the allocatable arrays, closes the file, and sets the error status/message
   !............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error ID (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! Set error status/message

      ErrStat = ErrID
      ErrMsg  = Msg
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg = 'Error in Init_DOFparameters: '//TRIM(ErrMsg)
      END IF


   END SUBROUTINE ExitThisRoutine

END SUBROUTINE Init_DOFparameters
!----------------------------------------------------------------------------------------------------------------------------------
FUNCTION SHP(Fract, FlexL, ModShpAry, Deriv, ErrStat, ErrMsg)
! SHP calculates the Derive-derivative of the shape function ModShpAry at Fract.
! NOTE: This function only works for Deriv = 0, 1, or 2.
!..................................................................................................................................

      ! Passed variables:

   REAL(ReKi),     INTENT(IN )    :: FlexL                     ! Length of flexible beam, (m)
   REAL(ReKi),     INTENT(IN )    :: Fract                     ! Fractional distance along flexible beam, 0<=Frac<=1
   REAL(ReKi),     INTENT(IN )    :: ModShpAry(:)              ! Array holding mode shape coefficients (2:PolyOrd)
   REAL(ReKi)                     :: SHP                       ! The shape function returned by this function.

   INTEGER(IntKi), INTENT(IN )    :: Deriv                     ! Which derivative to compute Deriv = 0 (regular function SHP), 1 (D(SHP)/DZ), 2 (D2(SHP)/DZ2)
   INTEGER(IntKi), INTENT(OUT)    :: ErrStat                   ! A error level that indicates if/what error occurred
   CHARACTER(*),   INTENT(OUT)    :: ErrMsg                    ! A message indicating the error if one occurred


      ! Local variables:

   INTEGER(IntKi)                 :: CoefTmp                   ! Temporary coefficient
   INTEGER(IntKi)                 :: I                         ! Counts through polynomial array.
   INTEGER(IntKi)                 :: J                         ! I+1
   INTEGER(IntKi)                 :: Swtch(0:2)                ! Corresponds to which derivative to compute.  Sets all portions of the coefficient = 0 except those that are relevant.


   IF ( Deriv < 0 .OR. Deriv > 2 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Function SHP input Deriv='//TRIM(Num2LStr(Deriv))//' is invalid. Deriv must be 0, 1, or 2.'
      RETURN
   ELSEIF ( Fract < 0.0_ReKi .OR. Fract > 1.0_ReKi ) THEN
      ErrStat = ErrID_Warn
      ErrMsg  = 'Function SHP input Fract='//TRIM(Num2LStr(Fract))//' does not meet the condition 0<=Fract<=1.'
   ELSE
      ErrStat = ErrID_None
   END IF

   Swtch        = 0 ! Initialize Swtch(:) to 0
   Swtch(Deriv) = 1
   SHP          = 0.0

   DO I = 1,SIZE(ModShpAry,DIM=1,KIND=IntKi) ! =2,PolyOrd
      J = I + 1
      CoefTmp = Swtch(0) + Swtch(1)*J + Swtch(2)*I*J

      IF ( (J == 2) .AND. (Deriv == 2) ) THEN !bjj this could be removed as Fract**0 = 1 (0**0 = 1 in Fortran)
         SHP =       ModShpAry(I)*CoefTmp                         /( FlexL**Deriv )
      ELSE
         SHP = SHP + ModShpAry(I)*CoefTmp*( Fract**( J - Deriv ) )/( FlexL**Deriv )
      ENDIF
   ENDDO !I

   RETURN

END FUNCTION SHP
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Alloc_CoordSys( CoordSys, p, ErrStat, ErrMsg )
! This subroutine allocates the coordinate systems in the ED_CoordSys type.
!..................................................................................................................................

IMPLICIT NONE

   ! passed arguments

TYPE(ED_CoordSys),        INTENT(OUT) :: CoordSys       ! The coordinate systems, with arrays to be allocated
TYPE(ED_ParameterType),   INTENT(IN)  :: p              ! Parameters of the structural dynamics module

INTEGER(IntKi),           INTENT(OUT) :: ErrStat        ! Error status
CHARACTER(*),             INTENT(OUT) :: ErrMsg         ! Err msg


   ! local variables

CHARACTER(200), PARAMETER        :: ErrTxt = 'coordinate system arrays in SUBROUTINE Alloc_CoordSys.'


   ! Initialize ErrStat and ErrMsg

ErrStat = ErrID_None
ErrMsg  = ""


  ! Allocate coordinate system arrays:

ALLOCATE ( CoordSys%i1(p%NumBl,3), CoordSys%i2(p%NumBl,3), CoordSys%i3(p%NumBl,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the i1, i2, and i3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%j1(p%NumBl,3), CoordSys%j2(p%NumBl,3), CoordSys%j3(p%NumBl,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the j1, j2, and j3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%m1(p%NumBl,p%BldNodes,3), CoordSys%m2(p%NumBl,p%BldNodes,3), &
           CoordSys%m3(p%NumBl,p%BldNodes,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the m1, m2, and m3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%n1(p%NumBl,0:p%TipNode,3), CoordSys%n2(p%NumBl,0:p%TipNode,3), &
           CoordSys%n3(p%NumBl,0:p%TipNode,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the n1, n2, and n3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%t1(p%TwrNodes,3), CoordSys%t2(p%TwrNodes,3), CoordSys%t3(p%TwrNodes,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the t1, t2, and t3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


ALLOCATE ( CoordSys%te1(p%NumBl,p%BldNodes,3), CoordSys%te2(p%NumBl,p%BldNodes,3), &
           CoordSys%te3(p%NumBl,p%BldNodes,3), STAT=ErrStat ) !this argument doesn't work in IVF 10.1: , ERRMSG=ErrMsg
IF ( ErrStat /= 0 )  THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Error allocating the te1, te2, and te3 '//TRIM(ErrTxt)//' '//TRIM(ErrMsg)
   RETURN
END IF


RETURN
END SUBROUTINE Alloc_CoordSys
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Alloc_BladeMeshInputProperties( BladeKInputFileMesh, ErrStat, ErrMsg )
! This routine allocates arrays for the blade mesh properties from the input file
!..................................................................................................................................

   TYPE(ED_BladeMeshInputData),   INTENT(INOUT)  :: BladeKInputFileMesh      ! Data for Blade K stored in the module's input file
   INTEGER(IntKi),                INTENT(OUT)    :: ErrStat                  ! Error status
   CHARACTER(*),                  INTENT(OUT)    :: ErrMsg                   ! Err msg


   IF ( BladeKInputFileMesh%BldNodes < 1 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating arrays for blade mesh input properties: BldNodes must be at least 1.'
      RETURN
   END IF

   CALL AllocAry  ( BladeKInputFileMesh%RNodes,   BladeKInputFileMesh%BldNodes, 'RNodes'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileMesh%AeroTwst, BladeKInputFileMesh%BldNodes, 'AeroTwst', ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileMesh%Chord,    BladeKInputFileMesh%BldNodes, 'Chord'   , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


END SUBROUTINE Alloc_BladeMeshInputProperties
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Alloc_BladeInputProperties( BladeKInputFileData, AllocAdams, ErrStat, ErrMsg )
! This routine allocates arrays for the blade properties from the input file
!..................................................................................................................................

   TYPE(BladeInputData),     INTENT(INOUT)  :: BladeKInputFileData      ! Data for Blade K stored in the module's input file
   LOGICAL,                  INTENT(IN)     :: AllocAdams               ! Logical to determine if we should allocate the arrays only used for Adams
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                  ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                   ! Err message


   IF ( BladeKInputFileData%NBlInpSt < 1 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating arrays for blade input properties: NBlInpSt must be at least 1.'
      RETURN
   END IF


      ! Allocate the arrays.

   CALL AllocAry  ( BladeKInputFileData%BlFract,  BladeKInputFileData%NBlInpSt, 'BlFract'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileData%PitchAx,  BladeKInputFileData%NBlInpSt, 'PitchAx'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileData%StrcTwst, BladeKInputFileData%NBlInpSt, 'StrcTwst' , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileData%BMassDen, BladeKInputFileData%NBlInpSt, 'BMassDen' , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileData%FlpStff,  BladeKInputFileData%NBlInpSt, 'FlpStff'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileData%EdgStff,  BladeKInputFileData%NBlInpSt, 'EdgStff'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


   IF ( AllocAdams ) THEN
      CALL AllocAry  ( BladeKInputFileData%GJStff,   BladeKInputFileData%NBlInpSt, 'GJStff'   , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%EAStff,   BladeKInputFileData%NBlInpSt, 'EAStff'   , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%Alpha,    BladeKInputFileData%NBlInpSt, 'Alpha'    , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%FlpIner,  BladeKInputFileData%NBlInpSt, 'FlpIner'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%EdgIner,  BladeKInputFileData%NBlInpSt, 'EdgIner'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%PrecrvRef,BladeKInputFileData%NBlInpSt, 'PrecrvRef', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%PreswpRef,BladeKInputFileData%NBlInpSt, 'PreswpRef', ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%FlpcgOf,  BladeKInputFileData%NBlInpSt, 'FlpcgOf'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%EdgcgOf,  BladeKInputFileData%NBlInpSt, 'EdgcgOf'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%FlpEAOf,  BladeKInputFileData%NBlInpSt, 'FlpEAOf'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( BladeKInputFileData%EdgEAOf,  BladeKInputFileData%NBlInpSt, 'EdgEAOf'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   END IF


      ! BJJ: note that these used to be allocated 2:PolyOrd  :

   CALL AllocAry  ( BladeKInputFileData%BldFl1Sh,  PolyOrd-1, 'BldFl1Sh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileData%BldFl2Sh,  PolyOrd-1, 'BldFl2Sh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileData%BldEdgSh,  PolyOrd-1, 'BldEdgSh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


END SUBROUTINE Alloc_BladeInputProperties
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ValidateBladeData ( BladeKInputFileData, ErrStat, ErrMsg )
! This routine checks the blade file input data for errors
!..................................................................................................................................
   TYPE(BladeInputData),     INTENT(IN)     :: BladeKInputFileData                 ! Data for Blade K stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message

      ! local variables
   INTEGER                                  :: I                                   ! Loop counter
   INTEGER(IntKi)                           :: ErrStat2                            ! Error status
   CHARACTER(ErrMsgLen)                     :: ErrMsg2                             ! Temporary error message


   ErrStat = ErrID_None
   ErrMsg= ''


      ! Check that BlFract goes from 0.0 to 1.0 in increasing order:

   IF ( .NOT. EqualRealNos( BladeKInputFileData%BlFract(1), 0.0_ReKi ) ) THEN
      CALL SetErrors( ErrID_Fatal,'BlFract(1) must be 0.0.')
   END IF

   IF ( BladeKInputFileData%NBlInpSt /= 1_IntKi .AND. &
      .NOT. EqualRealNos( BladeKInputFileData%BlFract(BladeKInputFileData%NBlInpSt), 1.0_ReKi )  ) THEN
      CALL SetErrors( ErrID_Fatal,'BlFract('//TRIM( Num2LStr( BladeKInputFileData%NBlInpSt ) )//') must be 1.0.')
   END IF

   DO I = 2,BladeKInputFileData%NBlInpSt
      IF ( BladeKInputFileData%BlFract(I) <= BladeKInputFileData%BlFract(I-1) )  THEN
         CALL SetErrors( ErrID_Fatal,'BlFract('//TRIM( Num2LStr( I ) )//') must be greater than BlFract('&
                                                      //TRIM( Num2LStr(I-1) )//').')
      ENDIF
   END DO


   DO I = 1,BladeKInputFileData%NBlInpSt

         ! Check that PitchAx is contained in [0.0, 1.0]:
      IF ( ( BladeKInputFileData%PitchAx(I) ) < 0.0_ReKi .OR. ( BladeKInputFileData%PitchAx(I) > 1.0_ReKi ) )  THEN
         CALL SetErrors( ErrID_Fatal,'PitchAx('//TRIM( Num2LStr( I ) )//') must be between 0 and 1 (inclusive).')
      END IF

         ! Check that StrcTwst is contained in (-pi,pi] radians ( i.e., (-180.0, 180.0] degrees):
      IF ( ( BladeKInputFileData%StrcTwst(I) <= -pi ) .OR. ( BladeKInputFileData%StrcTwst(I) > pi ) )  THEN
         CALL SetErrors( ErrID_Fatal,'StrcTwst('//TRIM( Num2LStr( I ) ) // &
                     ') must be greater than -180 and less than or equal to 180.')
      END IF

         ! Check that BMassDen is contained in (0.0, inf):
      IF ( BladeKInputFileData%BMassDen(I) <= 0.0_ReKi )  THEN
         CALL SetErrors( ErrID_Fatal,'BMassDen('//TRIM( Num2LStr( I ) )//') must be greater than zero.')
      END IF

         ! Check that FlpStff is contained in (0.0, inf):
      IF ( BladeKInputFileData%FlpStff (I) <= 0.0_ReKi )  THEN
         CALL SetErrors( ErrID_Fatal,'FlpStff('//TRIM( Num2LStr( I ) )//') must be greater than zero.')
      END IF

         ! Check that EdgStff is contained in (0.0, inf):
      IF ( BladeKInputFileData%EdgStff (I) <= 0.0_ReKi )  THEN
         CALL SetErrors( ErrID_Fatal,'EdgStff('//TRIM( Num2LStr( I ) )//') must be greater than zero.')
      END IF

   END DO


      ! Check values for Adams input

   IF ( ALLOCATED(BladeKInputFileData%GJStff) ) THEN  ! We assume that if GJStff is allocated, we are using ADAMS inputs

         ! The reference axis must be coincident with the pitch axis at the blade root (I == 1):
      IF ( .NOT. EqualRealNos( BladeKInputFileData%PrecrvRef(1), 0.0_ReKi ) .OR. &
            .NOT. EqualRealNos( BladeKInputFileData%PreswpRef(1), 0.0_ReKi )      )  THEN
         CALL SetErrors( ErrID_Fatal,'Both PrecrvRef(1) and PreswpRef(1) must be zero '//&
                            '(the reference axis must be coincident with the pitch axis at the blade root).')
      END IF


      DO I = 1,BladeKInputFileData%NBlInpSt

            ! Check that GJStff is contained in (0.0, inf):
         IF ( BladeKInputFileData%GJStff(I) <= 0.0_ReKi )  THEN
            CALL SetErrors( ErrID_Fatal,'GJStff('//TRIM( Num2LStr( I ) )//') must be greater than zero.')
         END IF

            ! Check that EAStff is contained in (0.0, inf):
         IF ( BladeKInputFileData%EAStff(I) <= 0.0_ReKi )  THEN
            CALL SetErrors( ErrID_Fatal,'EAStff('//TRIM( Num2LStr( I ) )//') must be greater than zero.')
         END IF

            ! Check that Alpha is contained in (-1.0, 1):
         IF ( ( BladeKInputFileData%Alpha(I) <= -1.0_ReKi ) .OR. ( BladeKInputFileData%Alpha(I) >= 1.0_ReKi ) )  THEN
            CALL SetErrors( ErrID_Fatal,'Alpha('//TRIM( Num2LStr( I ) )//') (the blade flap/twist'// &
                         ' coupling coefficient) must be between -1 and 1 (exclusive).')
         END IF

            ! Check that FlpIner is contained in [0.0, inf):
         IF ( BladeKInputFileData%FlpIner(I) <  0.0_ReKi )  THEN
            CALL SetErrors( ErrID_Fatal,'FlpIner('//TRIM( Num2LStr( I ) )//') must not be less than zero.')
         END IF

            ! Check that EdgIner is contained in [0.0, inf):
         IF ( BladeKInputFileData%EdgIner(I) <  0.0_ReKi )  THEN
            CALL SetErrors( ErrID_Fatal,'EdgIner('//TRIM( Num2LStr( I ) )//') must not be less than zero.')
         END IF

            ! Check that PrecrvRef is 0.0 for Adams models:
         IF ( .NOT. EqualRealNos( BladeKInputFileData%PrecrvRef(I), 0.0_ReKi) )  THEN
            CALL SetErrors( ErrID_Fatal,'PrecrvRef('//TRIM( Num2LStr( I ) )//') must be zero for Adams models.')
         END IF

            ! Check that GJStff is contained in (0.0, inf):
         IF ( .NOT. EqualRealNos( BladeKInputFileData%PreswpRef(I), 0.0_ReKi) )  THEN
            CALL SetErrors( ErrID_Fatal,'PreswpRef('//TRIM( Num2LStr( I ) )//') must be zero for Adams models.')
         END IF

      END DO

   END IF  ! check for Adams models


      ! Check that the blade damping is not negative:

   IF ( ANY( BladeKInputFileData%BldFlDmp < 0.0_ReKi ) ) CALL SetErrors( ErrID_Fatal,'BldFlDmp must not be negative.')
   IF ( ANY( BladeKInputFileData%BldEdDmp < 0.0_ReKi ) ) CALL SetErrors( ErrID_Fatal,'BldEdDmp must not be negative.')


      ! Check that the stiffness tuner isn't negative:

   IF ( ANY( BladeKInputFileData%FlStTunr <= 0.0_ReKi ) ) CALL SetErrors( ErrID_Fatal,'FlStTunr must be greater than zero.')


      ! Check that the mode shape coefficients are valid:

   CALL ValidateModeShapeCoeffs(BladeKInputFileData%BldFl1Sh, 'blade flap mode 1', ErrStat2, ErrMsg2 )
   CALL SetErrors( ErrStat2, ErrMsg2)

   CALL ValidateModeShapeCoeffs(BladeKInputFileData%BldFl2Sh, 'blade flap mode 2', ErrStat2, ErrMsg2 )
   CALL SetErrors( ErrStat2, ErrMsg2)

   CALL ValidateModeShapeCoeffs(BladeKInputFileData%BldEdgSh, 'blade edge', ErrStat2, ErrMsg2 )
   CALL SetErrors( ErrStat2, ErrMsg2)

CONTAINS
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE SetErrors( ErrStat3, ErrMsg3 )
   ! This routine sets the error message and flag when an error has occurred
   !...............................................................................................................................
   INTEGER(IntKi), INTENT(IN) :: ErrStat3     ! Error status for this error
   CHARACTER(*),   INTENT(IN) :: ErrMsg3      ! Error message for this error

      ErrStat = MAX( ErrStat, ErrStat3 )
      IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg  = TRIM(ErrMsg)//TRIM(ErrMsg3)

   END SUBROUTINE SetErrors
   !-------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ValidateBladeData
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ValidateModeShapeCoeffs( Coeffs, ShpDesc, ErrStat, ErrMsg )
! This routine checks that the mode shape coefficients add to 1.0, within numerical tolerance.
!..................................................................................................................................
   REAL(ReKi),               INTENT(IN )    :: Coeffs(:)                           ! Mode shape coefficients
   CHARACTER(*),             INTENT(IN)     :: ShpDesc                             ! Description of the mode shape for the error message
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message

      ! local variables
   REAL(ReKi)                               :: Displ                               ! Blade tip/tower top displacement for a mode shape


      ! Check that the mode shape coefficients add to 1.0:

   Displ = SUM( Coeffs )
! bjj this new check seems to be a bit too restrictive for the input data:
!   IF ( .NOT. EqualRealNos( Displ, 1.0_ReKi ) ) THEN
   IF ( ABS( Displ - 1.0_ReKi ) > 0.0015_ReKi ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = '  Mode shape coefficients for '//TRIM(ShpDesc)//' must add to 1.0.'
   ELSE
      ErrStat = ErrID_None
      ErrMsg  = ''
   END IF


END SUBROUTINE ValidateModeShapeCoeffs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetBladeParameters( p, BladeInData, BladeMeshData, ErrStat, ErrMsg )
! This takes the blade input file data and sets the corresponding blade parameters, performing linear interpolation of the
! input data to the specified blade mesh.
! This routine assumes p%HubRad and p%BldFlexL are already set.
!..................................................................................................................................

   TYPE(ED_ParameterType),        INTENT(INOUT)  :: p                                   ! The parameters of the structural dynamics module
   TYPE(BladeInputData),          INTENT(IN)     :: BladeInData(:)                      ! Program input data for all blades
   TYPE(ED_BladeMeshInputData),   INTENT(IN)     :: BladeMeshData(:)                    ! Program input mesh data for all blades
   INTEGER(IntKi),                INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),                  INTENT(OUT)    :: ErrMsg                              ! Error message

      ! Local variables:
   REAL(ReKi)                                    :: x                                   ! Fractional location between two points in linear interpolation
   INTEGER(IntKi )                               :: K                                   ! Blade number
   INTEGER(IntKi )                               :: J                                   ! Index for the node arrays
   INTEGER(IntKi)                                :: InterpInd                           ! Index for the interpolation routine
   LOGICAL                                       :: SetAdmVals                          ! Logical to determine if Adams inputs should be set

      ! initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''

   IF (p%BD4Blades) THEN
      SetAdmVals = .FALSE.
   ELSE
      SetAdmVals = ALLOCATED( BladeInData(1)%GJStff )
   END IF
   

   ! ..............................................................................................................................
   ! Set the blade discretization information here:
   ! ..............................................................................................................................

   DO K=1,1 ! we're going to assume the discretization is the same for all blades

      IF (p%BD4Blades) THEN
         p%BldNodes = 0
      ELSE         
         p%BldNodes = BladeMeshData(K)%BldNodes
      END IF

      p%TipNode  = p%BldNodes + 1    ! The index for the blade tip and tower top nodes

   END DO

      ! .......... Allocate arrays for the blade parameters being set in this routine ..........:

   CALL Alloc_BladeParameters( p, SetAdmVals, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   
   IF ( .not. p%BD4Blades) then
      
      DO K=1,1 ! we're going to assume the discretization is the same for all blades

         IF ( allocated( BladeMeshData(K)%Chord ) ) THEN      
      
            p%RNodes   = BladeMeshData(K)%RNodes - p%HubRad   ! Radius to blade analysis nodes relative to root ( 0 < RNodes(:) < p%BldFlexL ) (Convert RNodes to be relative to the hub)

            p%DRNodes(1) = 2.0*p%RNodes(1)
            DO J = 2,p%BldNodes
               p%DRNodes(J) = 2.0*( p%RNodes(J) - p%RNodes(J-1) ) - p%DRNodes(J-1)
            END DO

            p%Chord     = BladeMeshData(K)%Chord
            p%AeroTwst  = BladeMeshData(K)%AeroTwst
            p%CAeroTwst = COS(p%AeroTwst)
            p%SAeroTwst = SIN(p%AeroTwst)
         
         ELSE
         
         
               ! DRNodes (Let's use constant-spaced nodes for now, but the rest of the code is written to handle variable-spaced nodes--
               !          this will be a future input!):
            p%DRNodes = p%BldFlexL/p%BldNodes !array

               ! RNodes:
            p%RNodes(1) = 0.5*p%DRNodes(1)
            DO J=2,p%BldNodes
               p%RNodes(J) = p%RNodes( J - 1 ) + 0.5*( p%DRNodes(J) + p%DRNodes( J - 1 ) )
            END DO
         
               ! these values aren't used (at least they shouldn't be):
            p%Chord     = 0.0_ReKi
            p%AeroTwst  = 0.0_ReKi
            p%CAeroTwst = 1.0_ReKi
            p%SAeroTwst = 0.0_ReKi
                        
         END IF
         

      END DO


      ! ..............................................................................................................................
      ! Interpolate the blade properties to this discretization:
      ! ..............................................................................................................................

      ! Array definitions:

      !    Input      Interp    Description
      !    -----      ------    -----------
      !    BlFract    RNodesNorm Fractional radius (0 at root, 1 at tip)
      !    PitchAx    PitchAxis  Pitch axis (0 at LE, 1 at TE)
      !    StrcTwst   ThetaS     Structural twist
      !    BMassDen   MassB      Lineal mass density
      !    FlpStff    StiffBF    Flapwise stiffness
      !    EdgStff    StiffBE    Edgewise stiffness
      !    GJStff     StiffBGJ   Blade torsional stiffness
      !    EAStff     StiffBEA   Blade extensional stiffness
      !    Alpha      BAlpha     Blade flap/twist coupling coefficient
      !    FlpIner    InerBFlp   Blade flap (about local structural yb-axis) mass inertia per unit length
      !    EdgIner    InerBEdg   Blade edge (about local structural xb-axis) mass inertia per unit length
      !    PrecrvRef  RefAxisxb  Blade offset for defining the reference axis from the pitch axis for precurved blades (along xb-axis)
      !    PreswpRef  RefAxisyb  Blade offset for defining the reference axis from the pitch axis for preswept  blades (along yb-axis)
      !    FlpcgOf    cgOffBFlp  Blade flap mass cg offset
      !    EdgcgOf    cgOffBEdg  Blade edge mass cg offset
      !    FlpEAOf    EAOffBFlp  Blade flap elastic axis offset
      !    EdgEAOf    EAOffBEdg  Blade edge elastic axis offset


         ! Define RNodesNorm() which is common to all the blades:

      p%RNodesNorm = p%RNodes/p%BldFlexL  ! Normalized radius to analysis nodes relative to hub ( 0 < RNodesNorm(:) < 1 )
      
      

         ! Perform a linear interpolation of the input data to map to the meshed data for simulation:

      DO K=1,p%NumBl
         InterpInd = 1

         p%ThetaS  (K,0)         = BladeInData(K)%StrcTwst(1)
         p%ThetaS  (K,p%TipNode) = BladeInData(K)%StrcTwst(BladeInData(K)%NBlInpSt)
      
      
         DO J=1,p%BldNodes

               ! Get the index into BlFract for all of the arrays, using the NWTC Subroutine Library
            !p%ThetaS  (K,J) = InterpStp( p%RNodesNorm(J), BladeInData(K)%BlFract, BladeInData(K)%StrcTwst, &
            !                             InterpInd, BladeInData(K)%NBlInpSt )
            p%PitchAxis(K,J) = InterpStp( p%RNodesNorm(J), BladeInData(K)%BlFract, BladeInData(K)%PitchAx, &
                                         InterpInd, BladeInData(K)%NBlInpSt )


               ! The remaining arrays will have the same x value for the linear interpolation,
               ! so we'll do it manually (with a local subroutine) instead of calling the InterpStp routine again
            IF ( BladeInData(K)%NBlInpSt < 2_IntKi ) THEN
               x         = 1.0
               InterpInd = 0
            ELSE
               x = ( p%RNodesNorm(J)                     - BladeInData(K)%BlFract(InterpInd) ) / &
                   ( BladeInData(K)%BlFract(InterpInd+1) - BladeInData(K)%BlFract(InterpInd) )
            END IF

            p%ThetaS  (K,J) = InterpAry( x, BladeInData(K)%StrcTwst, InterpInd )
            p%MassB   (K,J) = InterpAry( x, BladeInData(K)%BMassDen, InterpInd )
            p%StiffBF (K,J) = InterpAry( x, BladeInData(K)%FlpStff , InterpInd )
            p%StiffBE (K,J) = InterpAry( x, BladeInData(K)%EdgStff , InterpInd )

            IF ( SetAdmVals ) THEN
               p%StiffBGJ (K,J) = InterpAry( x, BladeInData(K)%GJStff   , InterpInd )
               p%StiffBEA (K,J) = InterpAry( x, BladeInData(K)%EAStff   , InterpInd )
               p%BAlpha   (K,J) = InterpAry( x, BladeInData(K)%Alpha    , InterpInd )
               p%InerBFlp (K,J) = InterpAry( x, BladeInData(K)%FlpIner  , InterpInd )
               p%InerBEdg (K,J) = InterpAry( x, BladeInData(K)%EdgIner  , InterpInd )
               p%RefAxisxb(K,J) = InterpAry( x, BladeInData(K)%PrecrvRef, InterpInd )
               p%RefAxisyb(K,J) = InterpAry( x, BladeInData(K)%PreswpRef, InterpInd )
               p%cgOffBFlp(K,J) = InterpAry( x, BladeInData(K)%FlpcgOf  , InterpInd )
               p%cgOffBEdg(K,J) = InterpAry( x, BladeInData(K)%EdgcgOf  , InterpInd )
               p%EAOffBFlp(K,J) = InterpAry( x, BladeInData(K)%FlpEAOf  , InterpInd )
               p%EAOffBEdg(K,J) = InterpAry( x, BladeInData(K)%EdgEAOf  , InterpInd )
            END IF



         END DO ! J (Blade nodes)

         IF ( SetAdmVals ) THEN
               ! Set the valus for the tip node
            p%RefAxisxb(K,p%TipNode) = BladeInData(K)%PrecrvRef( BladeInData(K)%NBlInpSt )
            p%RefAxisyb(K,p%TipNode) = BladeInData(K)%PreswpRef( BladeInData(K)%NBlInpSt )
         END IF


            ! Set the blade damping and stiffness tuner
         p%BldFDamp(K,:) = BladeInData(K)%BldFlDmp
         p%BldEDamp(K,:) = BladeInData(K)%BldEdDmp
         p%FStTunr (K,:) = BladeInData(K)%FlStTunr



            ! Set the mode shape arrays
         p%BldEdgSh(:,K) = BladeInData(K)%BldEdgSh
         p%BldFl1Sh(:,K) = BladeInData(K)%BldFl1Sh
         p%BldFl2Sh(:,K) = BladeInData(K)%BldFl2Sh


      END DO ! ( Blades )

            
   else
      
      p%ThetaS  = 0.0_ReKi
      
         ! Set the blade damping and stiffness tuner
      p%BldFDamp = 0.0_ReKi
      p%BldEDamp = 0.0_ReKi
      p%FStTunr  = 0.0_ReKi

         ! Set the mode shape arrays
      p%BldEdgSh = 0.0_ReKi
      p%BldFl1Sh = 0.0_ReKi
      p%BldFl2Sh = 0.0_ReKi      
      
   end if
   
   p%CThetaS = COS(p%ThetaS)
   p%SThetaS = SIN(p%ThetaS)
   

RETURN


CONTAINS
!..................................................................................................................................
   FUNCTION InterpAry( x, YAry, Ind )
      ! This subroutine is used to interpolate the arrays more efficiently (all arrays have the same X value)
      ! See InterpStpReal() for comparison. This assumes we already know Ind and that
      ! x = ( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) )


      REAL(ReKi),      INTENT(IN) :: x                ! the relative distance between Ind and Ind+ 1
      REAL(ReKi),      INTENT(IN) :: YAry (:)         ! Array of Y values to be interpolated.
      INTEGER(IntKi) , INTENT(IN) :: Ind              ! the index into the array

      REAL(ReKi)                  :: InterpAry        ! the value calculated in this function

      IF ( Ind >= SIZE(YAry) ) THEN
         InterpAry = YAry( SIZE(YAry) )
      ELSE
         InterpAry = ( YAry(Ind+1) - YAry(Ind) ) * x  + YAry(Ind)
      END IF

   END FUNCTION InterpAry
!..................................................................................................................................
END SUBROUTINE SetBladeParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Alloc_BladeParameters( p, AllocAdams, ErrStat, ErrMsg )
! This routine allocates arrays for the blade parameters.
!..................................................................................................................................

   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                                   ! The parameters of the structural dynamics module
   LOGICAL,                  INTENT(IN)     :: AllocAdams                          ! Logical to determine if Adams inputs should be allocated
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Err msg


      ! Allocate arrays to hold the blade analysis nodes.
   CALL AllocAry  ( p%RNodes,             p%BldNodes, 'RNodes'   , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%DRNodes,            p%BldNodes, 'DRNodes'  , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%Chord,              p%BldNodes, 'Chord'    , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%AeroTwst,           p%BldNodes, 'AeroTwst' , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%CAeroTwst,          p%BldNodes, 'CAeroTwst', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%SAeroTwst,          p%BldNodes, 'SAeroTwst', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN


      ! Allocate arrays to hold blade data at the analysis nodes.
   CALL AllocAry  ( p%RNodesNorm,              p%BldNodes, 'RNodesNorm' , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%PitchAxis,   p%NumBl,    p%BldNodes, 'PitchAxis'  , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   
   ALLOCATE( p%ThetaS( p%NumBl,0:P%TipNode) &
           , p%CThetaS(p%NumBl,0:P%TipNode) &
           , p%SThetaS(p%NumBl,0:P%TipNode), STAT=ErrStat ) 
   IF (ErrStat /= 0) then
      ErrStat = ErrID_Fatal
      ErrMsg = 'Error allocating ThetaS, CThetaS, and SThetaS'
   END IF
   
      
   CALL AllocAry  ( p%MassB,       p%NumBl,    p%BldNodes, 'MassB'      , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%StiffBF,     p%NumBl,    p%BldNodes, 'StiffBF'    , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%StiffBE,     p%NumBl,    p%BldNodes, 'StiffBE'    , ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN


   IF ( AllocAdams ) THEN
      CALL AllocAry  ( p%StiffBGJ,    p%NumBl,    p%BldNodes, 'StiffBGJ'   , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%StiffBEA,    p%NumBl,    p%BldNodes, 'StiffBEA'   , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%BAlpha,      p%NumBl,    p%BldNodes, 'BAlpha'     , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%InerBFlp,    p%NumBl,    p%BldNodes, 'InerBFlp'   , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%InerBEdg,    p%NumBl,    p%BldNodes, 'InerBEdg'   , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%RefAxisxb,   p%NumBl,    p%TipNode,  'RefAxisxb'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%RefAxisyb,   p%NumBl,    p%TipNode,  'RefAxisyb'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%cgOffBFlp,   p%NumBl,    p%BldNodes, 'cgOffBFlp'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%cgOffBEdg,   p%NumBl,    p%BldNodes, 'cgOffBEdg'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%EAOffBFlp,   p%NumBl,    p%BldNodes, 'EAOffBFlp'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%EAOffBEdg,   p%NumBl,    p%BldNodes, 'EAOffBEdg'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   END IF

   CALL AllocAry  ( p%BldEDamp,    p%NumBl,    NumBE,      'BldEDamp'   , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%BldFDamp,    p%NumBl,    NumBF,      'BldFDamp'   , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%FStTunr,     p%NumBl,    NumBF,      'FStTunr'    , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


         ! Allocate space for the mode shape arrays:

   ALLOCATE( p%BldEdgSh(2:PolyOrd,p%NumBl), p%BldFl1Sh(2:PolyOrd,p%NumBl), p%BldFl2Sh(2:PolyOrd,p%NumBl), STAT = ErrStat )
   IF ( ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Error allocating BldEdgSh, BldFl1Sh, and BldFl2Sh arrays.'
      RETURN
   END IF


END SUBROUTINE Alloc_BladeParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ValidateTowerData ( InputFileData, ErrStat, ErrMsg )
! This routine checks the tower file input data for errors
!..................................................................................................................................
   TYPE(ED_InputFile),     INTENT(IN   )    :: InputFileData                       ! Data stored in the module's input file
   INTEGER(IntKi),         INTENT(  OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),           INTENT(  OUT)    :: ErrMsg                              ! Error message

      ! local variables
   INTEGER                                  :: I                                   ! Loop counter
   INTEGER(IntKi)                           :: ErrStat2                            ! Error status
   CHARACTER(ErrMsgLen)                     :: ErrMsg2                             ! Temporary error message


   ErrStat = ErrID_None
   ErrMsg= ''



      ! Check that HtFract goes from 0.0 to 1.0 in increasing order:

   IF ( .NOT. EqualRealNos( InputFileData%HtFract(1), 0.0_ReKi ) ) CALL SetErrors( ErrID_Fatal, 'HtFract(1) must be 0.0.')

   IF ( InputFileData%NTwInpSt /= 1 .AND. &
      .NOT. EqualRealNos( InputFileData%HtFract(InputFileData%NTwInpSt), 1.0_ReKi )  ) THEN
      CALL SetErrors( ErrID_Fatal, 'HtFract('//TRIM( Num2LStr( InputFileData%NTwInpSt ) )//') must be 1.0.')
   END IF

   DO I = 2,InputFileData%NTwInpSt
      IF ( InputFileData%HtFract(I) <= InputFileData%HtFract(I-1) )  THEN
         CALL SetErrors( ErrID_Fatal, 'HtFract('//TRIM( Num2LStr( I ) )//') must be greater than HtFract('&
                                                      //TRIM( Num2LStr(I-1) )//').')

      ENDIF
   END DO


      ! Check the input arrays:

   DO I = 1,InputFileData%NTwInpSt
      IF ( InputFileData%TMassDen(I) <= 0.0_ReKi ) THEN
         CALL SetErrors( ErrID_Fatal, 'TMassDen('//TRIM(Num2LStr( I ))//') must be greater than zero.')
      END IF

      IF ( InputFileData%TwFAStif(I) <= 0.0_ReKi ) THEN
         CALL SetErrors( ErrID_Fatal, 'TwFAStif('//TRIM(Num2LStr( I ))//') must be greater than zero.')
      END IF

      IF ( InputFileData%TwSSStif(I) <= 0.0_ReKi ) THEN
         CALL SetErrors( ErrID_Fatal, 'TwSSStif('//TRIM(Num2LStr( I ))//') must be greater than zero.')
      END IF
   END DO

      ! Check Adams inputs

   IF ( ALLOCATED( InputFileData%TwGJStif ) ) THEN ! Assume that all of the Adams tower data is allocated

      DO I = 1,InputFileData%NTwInpSt
         IF ( InputFileData%TwGJStif(I) <= 0.0_ReKi ) THEN
            CALL SetErrors( ErrID_Fatal, 'TwGJStif('//TRIM(Num2LStr( I ))//') must be greater than zero.')
         END IF

         IF ( InputFileData%TwEAStif(I) <= 0.0_ReKi ) THEN
            CALL SetErrors( ErrID_Fatal, 'TwEAStif('//TRIM(Num2LStr( I ))//') must be greater than zero.')
         END IF

         IF ( InputFileData%TwFAIner(I) <= 0.0_ReKi ) THEN
            CALL SetErrors( ErrID_Fatal, 'TwFAIner('//TRIM(Num2LStr( I ))//') must be greater than zero.')
         END IF

         IF ( InputFileData%TwSSIner(I) <= 0.0_ReKi ) THEN
            CALL SetErrors( ErrID_Fatal, 'TwSSIner('//TRIM(Num2LStr( I ))//') must be greater than zero.')
         END IF
      END DO

   END IF ! Check items for Adams



      ! Check that the tower damping (TwrFADmp) is contained in the range [0, 100]:

   IF ( ANY( InputFileData%TwrFADmp < 0.0_ReKi ) .OR. ANY( InputFileData%TwrFADmp > 100.0_ReKi ) ) THEN
      CALL SetErrors( ErrID_Fatal, 'TwrFADmp must be between 0 and 100 (inclusive).')
   END IF

   IF ( ANY( InputFileData%TwrSSDmp < 0.0_ReKi ) .OR. ANY( InputFileData%TwrSSDmp > 100.0_ReKi ) ) THEN
      CALL SetErrors( ErrID_Fatal, 'TwrSSDmp must be between 0 and 100 (inclusive).')
   END IF


      ! Check that the tower tuners are positive numbers:

   IF ( ANY( InputFileData%FAStTunr <= 0.0_ReKi )  ) CALL SetErrors( ErrID_Fatal, 'FAStTunr must be greater than zero.' )
   IF ( ANY( InputFileData%SSStTunr <= 0.0_ReKi )  ) CALL SetErrors( ErrID_Fatal, 'SSStTunr must be greater than zero.' )



      ! Validate the mode shape coefficients:

   CALL ValidateModeShapeCoeffs( InputFileData%TwFAM1Sh, 'tower fore-aft mode 1', ErrStat2, ErrMsg2 )
      CALL SetErrors( ErrStat2, ErrMsg2 )


   CALL ValidateModeShapeCoeffs( InputFileData%TwFAM2Sh, 'tower fore-aft mode 2', ErrStat2, ErrMsg2 )
      CALL SetErrors( ErrStat2, ErrMsg2 )


   CALL ValidateModeShapeCoeffs( InputFileData%TwSSM1Sh, 'tower side-to-side mode 1', ErrStat2, ErrMsg2 )
      CALL SetErrors( ErrStat2, ErrMsg2 )


   CALL ValidateModeShapeCoeffs( InputFileData%TwSSM2Sh, 'tower side-to-side mode 2', ErrStat2, ErrMsg2 )
      CALL SetErrors( ErrStat2, ErrMsg2 )

CONTAINS
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE SetErrors( ErrStat3, ErrMsg3 )
   ! This routine sets the error message and flag when an error has occurred
   !...............................................................................................................................
   INTEGER(IntKi), INTENT(IN) :: ErrStat3     ! Error status for this error
   CHARACTER(*),   INTENT(IN) :: ErrMsg3      ! Error message for this error

      ErrStat = MAX( ErrStat, ErrStat3 )
      IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg  = TRIM(ErrMsg)//TRIM(ErrMsg3)

   END SUBROUTINE SetErrors
   !-------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ValidateTowerData
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Alloc_TowerInputProperties( InputFileData, AllocAdams, ErrStat, ErrMsg )
! This routine allocates arrays for the tower properties from the input file
!..................................................................................................................................

   TYPE(ED_InputFile),       INTENT(INOUT)  :: InputFileData      ! All the data in the ElastoDyn input file
   LOGICAL,                  INTENT(IN)     :: AllocAdams         ! Determines if the columns for Adams data will be read
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat            ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg             ! Error message


   IF ( InputFileData%NTwInpSt < 1 )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating arrays for tower input properties: NTwInpSt must be at least 1.'
      RETURN
   END IF


      ! Allocate the arrays.

   CALL AllocAry  ( InputFileData%HtFract,   InputFileData%NTwInpSt, 'HtFract'   , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( InputFileData%TMassDen,  InputFileData%NTwInpSt, 'TMassDen'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( InputFileData%TwFAStif,  InputFileData%NTwInpSt, 'TwFAStif'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( InputFileData%TwSSStif,  InputFileData%NTwInpSt, 'TwSSStif'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   IF ( AllocAdams ) THEN
      CALL AllocAry  ( InputFileData%TwGJStif,  InputFileData%NTwInpSt, 'TwGJStif'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( InputFileData%TwEAStif,  InputFileData%NTwInpSt, 'TwEAStif'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( InputFileData%TwFAIner,  InputFileData%NTwInpSt, 'TwFAIner'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( InputFileData%TwSSIner,  InputFileData%NTwInpSt, 'TwSSIner'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( InputFileData%TwFAcgOf,  InputFileData%NTwInpSt, 'TwFAcgOf'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( InputFileData%TwSScgOf,  InputFileData%NTwInpSt, 'TwSScgOf'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   END IF


      ! BJJ: note that these used to be allocated 2:PolyOrd  :
   CALL AllocAry  ( InputFileData%TwFAM1Sh,  PolyOrd-1, 'TwFAM1Sh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( InputFileData%TwFAM2Sh,  PolyOrd-1, 'TwFAM2Sh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( InputFileData%TwSSM1Sh,  PolyOrd-1, 'TwSSM1Sh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( InputFileData%TwSSM2Sh,  PolyOrd-1, 'TwSSM2Sh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


END SUBROUTINE Alloc_TowerInputProperties
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Alloc_TowerParameters( p, AllocAdams, ErrStat, ErrMsg )
! This routine allocates arrays for the tower parameters.
!..................................................................................................................................

   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                                   ! The parameters of the structural dynamics module
   LOGICAL,                  INTENT(IN)     :: AllocAdams                          ! Logical to determine if Adams inputs should be allocated
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Err msg




      ! Allocate arrays to hold tower data at the analysis nodes.
   CALL AllocAry  ( p%HNodesNorm,    p%TwrNodes, 'HNodesNorm', ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%HNodes,        p%TwrNodes, 'HNodes'    , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%DHNodes,       p%TwrNodes, 'DHNodes'   , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%MassT,         p%TwrNodes, 'MassT'     , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%StiffTFA,      p%TwrNodes, 'StiffTFA'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( p%StiffTSS,      p%TwrNodes, 'StiffTSS'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   IF ( AllocAdams ) THEN
      CALL AllocAry  ( p%StiffTGJ,      p%TwrNodes, 'StiffTGJ'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%StiffTEA,      p%TwrNodes, 'StiffTEA'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%InerTFA,       p%TwrNodes, 'InerTFA'   , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%InerTSS,       p%TwrNodes, 'InerTSS'   , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%cgOffTFA,      p%TwrNodes, 'cgOffTFA'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
      CALL AllocAry  ( p%cgOffTSS,      p%TwrNodes, 'cgOffTSS'  , ErrStat, ErrMsg )
      IF ( ErrStat /= ErrID_None ) RETURN
   END IF

   !   ! these are for HydroDyn?
   !CALL AllocAry  ( p%DiamT,         p%TwrNodes, 'DiamT'     , ErrStat, ErrMsg )
   !IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry  ( p%CAT,           p%TwrNodes, 'CAT'       , ErrStat, ErrMsg )
   !IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry  ( p%CDT,           p%TwrNodes, 'CDT'       , ErrStat, ErrMsg )
   !IF ( ErrStat /= ErrID_None ) RETURN



   !      ! Allocate space for the mode shape arrays:
   !
   !ALLOCATE( p%BldEdgSh(2:PolyOrd,p%NumBl), p%BldFl1Sh(2:PolyOrd,p%NumBl), p%BldFl2Sh(2:PolyOrd,p%NumBl), STAT = ErrStat )
   !IF ( ErrStat /= 0 ) THEN
   !   ErrStat = ErrID_Fatal
   !   ErrMsg  = ' Error allocating BldEdgSh, BldFl1Sh, and BldFl2Sh arrays.'
   !   RETURN
   !END IF



END SUBROUTINE Alloc_TowerParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetOtherParameters( p, InputFileData, ErrStat, ErrMsg )
! This routine sets the remaining parameters (replacing the former FAST Initialize routine), first allocating necessary arrays.
! It requires p%NDOF, p%NumBl, p%TTopNode, p%TipNode to be set before calling this routine.
!..................................................................................................................................

   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                ! Data stored in the module's input file
   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                            ! The parameters of the structural dynamics module
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       ! Err msg




      ! Allocate the arrays needed in the Coeff routine:

   !CALL AllocAry( p%AxRedTFA, 2,       2_IntKi, p%TTopNode,         'AxRedTFA',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( p%AxRedTSS, 2,       2_IntKi, p%TTopNode,         'AxRedTSS',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( p%AxRedBld, p%NumBl, 3_IntKi, 3_IntKi, p%TipNode, 'AxRedBld',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%BldCG,    p%NumBl,                              'BldCG',     ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%KBF,      p%NumBl, 2_IntKi, 2_IntKi,            'KBF',       ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%KBE,      p%NumBl, 1_IntKi, 1_IntKi,            'KBE',       ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%CBF,      p%NumBl, 2_IntKi, 2_IntKi,            'CBF',       ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%CBE,      p%NumBl, 1_IntKi, 1_IntKi,            'CBE',       ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%SecondMom,p%NumBl,                              'SecondMom', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%FirstMom, p%NumBl,                              'FirstMom',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%FreqBE,   p%NumBl, NumBE, 3_IntKi,              'FreqBE',    ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%FreqBF,   p%NumBl, NumBF, 3_IntKi,              'FreqBF',    ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%BldMass,  p%NumBl,                              'BldMass',   ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%rSAerCenn1,p%NumBl,p%BldNodes,  'rSAerCenn1',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%rSAerCenn2,p%NumBl,p%BldNodes,  'rSAerCenn2',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry(p%BElmntMass, p%BldNodes, p%NumBl, 'BElmntMass', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry(p%TElmntMass, p%TwrNodes,          'TElmntMass', ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN

   !CALL AllocAry( p%AxRedBld, p%NumBl, 3_IntKi, 3_IntKi, p%TipNode, 'AxRedBld',  ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) RETURN
   ALLOCATE ( p%AxRedBld(p%NumBl, 3_IntKi, 3_IntKi, 0:p%TipNode) , STAT=ErrStat )
   IF ( ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating AxRedBld array.'
      RETURN
   END IF

   
   ALLOCATE ( p%TwrFASF(2,0:p%TTopNode,0:2) , &
              p%TwrSSSF(2,0:p%TTopNode,0:2) , & 
              p%AxRedTFA(2,2,0:p%TTopNode)  , &
              p%AxRedTSS(2,2,0:p%TTopNode)  , STAT=ErrStat )
   IF ( ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating TwrFASF, TwrSSSF, AxRedTFA, and p%AxRedTSS arrays.'
      RETURN
   END IF

   ALLOCATE ( p%TwistedSF(p%NumBl,2,3,0:p%TipNode,0:2) , STAT=ErrStat )
   IF ( ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating TwistedSF array.'
      RETURN
   END IF


   CALL Coeff(p, InputFileData, ErrStat, ErrMsg)
   IF ( ErrStat /= ErrID_None ) RETURN
   
   
END SUBROUTINE SetOtherParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Alloc_RtHS( RtHS, p, ErrStat, ErrMsg  )
! This routine allocates arrays in the RtHndSide data structure.
! It requires p%TwrNodes, p%NumBl, p%TipNode, p%NDOF, p%BldNodes to be set before calling this routine.
!..................................................................................................................................

   TYPE(ED_RtHndSide),       INTENT(INOUT)  :: RtHS                         ! RtHndSide data type
   TYPE(ED_ParameterType),   INTENT(IN)     :: p                            ! Parameters of the structural dynamics module
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       ! Error message

   ! local variables:
   INTEGER(IntKi),   PARAMETER              :: Dims = 3                     ! The position arrays all must be allocated with a dimension for X,Y,and Z
   CHARACTER(*),     PARAMETER              :: RoutineName = 'Alloc_RtHS'

      ! positions:
  !CALL AllocAry( RtHS%rZT,       Dims, p%TwrNodes,        'rZT',       ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rT,        Dims, p%TwrNodes,        'rT',        ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rT0T,      Dims, p%TwrNodes,        'rT0T',      ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
  !CALL AllocAry( RtHS%rQS,       Dims, p%NumBl,p%TipNode, 'rQS',       ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
  !CALL AllocAry( RtHS%rS,        Dims, p%NumBl,p%TipNode, 'rS',        ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rS0S,      Dims, p%NumBl,p%TipNode, 'rS0S',      ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rPS0,      Dims, p%NumBl,           'rPS0',      ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%rSAerCen,  Dims, p%TipNode, p%NumBl,'rSAerCen',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
  
      ! tower
   allocate(RtHS%rZT(      Dims,  0:p%TwrNodes), &
            RtHS%AngPosEF( Dims,  0:p%TwrNodes), &
            RtHS%AngPosXF( Dims,  0:p%TwrNodes), &
            RtHS%AngVelEF( Dims,  0:p%TwrNodes), &
            RtHS%LinVelET( Dims,  0:p%TwrNodes), &
            RtHS%AngAccEFt(Dims,  0:p%TwrNodes), &
            RtHS%LinAccETt(Dims,  0:p%TwrNodes), &
      STAT=ErrStat)
      if (ErrStat /= 0) then
         ErrStat = ErrID_Fatal
         ErrMsg  = "Error allocating rZT, AngPosEF, AngPosXF, LinVelET, AngVelEF, LinAccETt, and AngAccEFt arrays."
         RETURN
      end if
               
      
      ! blades
   allocate(RtHS%rS( Dims, p%NumBl,0:p%TipNode), &
            RtHS%rQS(Dims, p%NumBl,0:p%TipNode), STAT=ErrStat)
      if (ErrStat /= 0) then
         ErrStat = ErrID_Fatal
         ErrMsg  = "Error allocating rS and rQS."
         RETURN
      end if
   

      ! angular velocities (including partial angular velocities):
   !CALL AllocAry( RtHS%AngVelEF,  Dims, p%TwrNodes,        'AngVelEF',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( RtHS%AngPosEF,  Dims, p%TwrNodes,        'AngPosEF',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( RtHS%AngPosXF,  Dims, p%TwrNodes,        'AngPosXF',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN


         ! These angular velocities are allocated to start numbering a dimension with 0 instead of 1:
   ALLOCATE ( RtHS%PAngVelEB(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEB array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PAngVelER(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelER array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PAngVelEX(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEX array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PAngVelEA(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEA array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PAngVelEF(0:p%TwrNodes, p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEF array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEG(                  p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEG array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEH(                  p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEH array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEL(                  p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEL array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEM(p%NumBl,p%TipNode,p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEM array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PAngVelEN(                  p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PAngVelEN array.'
      RETURN
   ENDIF

      ! angular accelerations:
   !CALL AllocAry( RtHS%AngAccEFt, Dims, p%TwrNodes,         'AngAccEFt', ErrStat, ErrMsg );  IF ( ErrStat /= ErrID_None ) RETURN

      ! linear velocities (including partial linear velocities):
   !CALL AllocAry( RtHS%LinVelET,  Dims, p%TwrNodes,         'LinVelET',  ErrStat, ErrMsg );  IF ( ErrStat /= ErrID_None ) RETURN         

   !CALL AllocAry( RtHS%LinVelESm2,                 p%NumBl, 'LinVelESm2',ErrStat, ErrMsg );  IF ( ErrStat /= ErrID_None ) RETURN ! The m2-component (closest to tip) of LinVelES
   ALLOCATE( RtHS%LinVelES( Dims, 0:p%TipNode, p%NumBl ), STAT=ErrStat )
   IF (ErrStat /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = RoutineName//":Error allocating LinVelES."
      RETURN
   END IF
   
            ! These linear velocities are allocated to start numbering a dimension with 0 instead of 1:

   ALLOCATE ( RtHS%PLinVelEIMU(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEIMU array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelEO(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEO array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelES(p%NumBl,0:p%TipNode,p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelES array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelET(0:p%TwrNodes,p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelET array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelEZ(p%NDOF,0:1,Dims) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEZ array.'
      RETURN
   ENDIF

   ALLOCATE ( RtHS%PLinVelEC(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEC array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelED(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelED array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelEI(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEI array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelEJ(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEJ array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelEK(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEK array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelEP(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEP array.'
      RETURN
   ENDIF
   ALLOCATE ( RtHS%PLinVelEQ(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEQ array.'
      RETURN
   ENDIF
   
   ALLOCATE ( RtHS%PLinVelEU(p%NDOF,0:1,3) , &
              RtHS%PLinVelEV(p%NDOF,0:1,3) , &
              RtHS%PLinVelEW(p%NDOF,0:1,3) , &
              RtHS%PLinVelEY(p%NDOF,0:1,3) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for the PLinVelEU, PLinVelEV, PLinVelEW and PLinVelEY arrays.'
      RETURN
   ENDIF

   
   ALLOCATE( RtHS%LinAccESt( Dims, p%NumBl, 0:p%TipNode ), STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' Error allocating memory for LinAccESt.'
      RETURN
   ENDIF


   !CALL AllocAry( RtHS%LinAccESt, Dims, p%NumBl, p%TipNode,'LinAccESt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   !CALL AllocAry( RtHS%LinAccETt, Dims, p%TwrNodes,        'LinAccETt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcS0B,   Dims, p%NumBl,p%NDOF,    'PFrcS0B',   ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%FrcS0Bt,   Dims, p%NumBl,           'FrcS0Bt',   ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomH0B,   Dims, p%NumBl, p%NDOF,   'PMomH0B',   ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%MomH0Bt,   Dims, p%NumBl,           'MomH0Bt',   ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcPRot,  Dims, p%NDOF,            'PFrcPRot',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomLPRot, Dims, p%NDOF,            'PMomLPRot', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomNGnRt, Dims, p%NDOF,            'PMomNGnRt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomNTail, Dims, p%NDOF,            'PMomNTail', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcONcRt, Dims, p%NDOF,            'PFrcONcRt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomBNcRt, Dims, p%NDOF,            'PMomBNcRt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcT0Trb, Dims, p%NDOF,            'PFrcT0Trb', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomX0Trb, Dims, p%NDOF,            'PMomX0Trb', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%FSAero,    Dims, p%NumBl,p%BldNodes,'FSAero',    ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%MMAero,    Dims, p%NumBl,p%BldNodes,'MMAero',    ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%FSTipDrag, Dims, p%NumBl,           'FSTipDrag', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%AngPosHM,  Dims, p%NumBl,p%TipNode, 'AngPosHM',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFTHydro,  Dims, p%TwrNodes, p%NDOF,'PFTHydro',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMFHydro,  Dims, p%TwrNodes, p%NDOF,'PMFHydro',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%FTHydrot,  Dims, p%TwrNodes,        'FTHydrot',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%MFHydrot,  Dims, p%TwrNodes,        'MFHydrot',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN

   CALL AllocAry( RtHS%PFrcVGnRt, Dims, p%NDOF,            'PFrcVGnRt', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcWTail, Dims, p%NDOF,            'PFrcWTail', ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PFrcZAll,  Dims, p%NDOF,            'PFrcZAll',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( RtHS%PMomXAll,  Dims, p%NDOF,            'PMomXAll',  ErrStat, ErrMsg );   IF ( ErrStat /= ErrID_None ) RETURN

END SUBROUTINE Alloc_RtHS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetTowerParameters( p, InputFileData, ErrStat, ErrMsg  )
! This takes the tower input file data and sets the corresponding tower parameters, performing linear interpolation of the
! input data to the specified tower mesh.
! It requires p%TwrFlexL, and p%TwrNodes to be set first.
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed variables

   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                            ! Parameters of the structural dynamics module
   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                ! Data stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       ! Error message

      ! Local variables:

   REAL(ReKi)                               :: x                            ! Fractional location between two points in linear interpolation
   INTEGER(IntKi )                          :: J                            ! Index for the node arrays
   INTEGER(IntKi)                           :: InterpInd                    ! Index for the interpolation routine
   LOGICAL                                  :: SetAdmVals                   ! Logical to determine if Adams inputs should be set


      ! Initialize data
   ErrStat   = ErrID_None
   ErrMsg    = ''
   SetAdmVals = ALLOCATED( InputFileData%TwGJStif )

   CALL Alloc_TowerParameters( p, SetAdmVals, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


   !...............................................................................................................................
   ! Define the tower discretization arrays:
   !...............................................................................................................................

      ! DHNodes (Let's use constant-spaced nodes for now, but the rest of the code is written to handle variable-spaced nodes--
      !          this will be a future input!):
   p%DHNodes = p%TwrFlexL/p%TwrNodes

      ! HNodes:
   p%HNodes(1) = 0.5*p%DHNodes(1)
   DO J=2,p%TwrNodes
      p%HNodes(J) = p%HNodes( J - 1 ) + 0.5*( p%DHNodes(J) + p%DHNodes( J - 1 ) )
   END DO

      ! HNodesNorm:
   p%HNodesNorm = p%HNodes/p%TwrFlexL


   !...............................................................................................................................
   ! Interpolate the input data to the tower discretization
   !...............................................................................................................................
   ! Array definitions:

   !    Input      Interp    Description
   !    -----      ------    -----------
   !    HtFract    HNodesNorm Fractional height (0 at top of rigid section, 1 at tower top)
   !    TMassDen   MassT      Lineal mass density
   !    TwFAStif   StiffTFA   Tower fore-aft stiffness
   !    TwSSStif   StiffTSS   Tower side-to-side stiffness
   !    TwGJStif   StiffTGJ   Tower torsional stiffness
   !    TwEAStif   StiffTEA   Tower extensional stiffness
   !    TwFAIner   InerTFA    Tower fore-aft (about yt-axis) mass inertia per unit length
   !    TwSSIner   InerTSS    Tower side-to-side (about xt-axis) mass inertia per unit length
   !    TwFAcgOf   cgOffTFA   Tower fore-aft mass cg offset
   !    TwSScgOf   cgOffTSS   Tower side-to-side mass cg offset

   InterpInd = 1


   DO J=1,p%TwrNodes

         ! Get the index into HtFract for all of the arrays, using the NWTC Subroutine Library
      p%MassT     (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TMassDen, InterpInd, InputFileData%NTwInpSt )
      p%StiffTFA  (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwFAStif, InterpInd, InputFileData%NTwInpSt )
      p%StiffTSS  (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwSSStif, InterpInd, InputFileData%NTwInpSt )
   END DO ! J


   IF ( SetAdmVals )  THEN          ! An ADAMS model will be created; thus, read in all the cols.
      DO J=1,p%TwrNodes
         p%StiffTGJ  (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwGJStif, InterpInd, InputFileData%NTwInpSt )
         p%StiffTEA  (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwEAStif, InterpInd, InputFileData%NTwInpSt )
         p%InerTFA   (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwFAIner, InterpInd, InputFileData%NTwInpSt )
         p%InerTSS   (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwSSIner, InterpInd, InputFileData%NTwInpSt )
         p%cgOffTFA  (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwFAcgOf, InterpInd, InputFileData%NTwInpSt )
         p%cgOffTSS  (J) = InterpStp( p%HNodesNorm(J), InputFileData%HtFract, InputFileData%TwSScgOf, InterpInd, InputFileData%NTwInpSt )
      END DO ! J
   END IF


   !...............................................................................................................................
   ! Set other tower parameters:
   !...............................................................................................................................

   p%TTopNode = p%TwrNodes + 1

   !   ! these are for HydroDyn ?
   !p%DiamT(:) = InputFileData%TwrDiam
   !p%CAT(:)   = InputFileData%TwrCA
   !p%CDT(:)   = InputFileData%TwrCD
   !

RETURN


CONTAINS
!..................................................................................................................................
   FUNCTION InterpAry( x, YAry, Ind )
      ! This subroutine is used to interpolate the arrays more efficiently (all arrays have the same X value)
      ! See InterpStpReal() for comparison. This assumes we already know Ind and that
      ! x = ( XVal - XAry(Ind) )/( XAry(Ind+1) - XAry(Ind) )


      REAL(ReKi),      INTENT(IN) :: x                ! the relative distance between Ind and Ind+ 1
      REAL(ReKi),      INTENT(IN) :: YAry (:)         ! Array of Y values to be interpolated.
      INTEGER(IntKi) , INTENT(IN) :: Ind              ! the index into the array

      REAL(ReKi)                  :: InterpAry        ! the value calculated in this function

      InterpAry = ( YAry(Ind+1) - YAry(Ind) ) * x  + YAry(Ind)

   END FUNCTION InterpAry

END SUBROUTINE SetTowerParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ValidateFurlData( InputFileData, ErrStat, ErrMsg )
! This routine validates the furling inputs.
!..................................................................................................................................

      ! Passed variables:

   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                       ! All the data in the ElastoDyn input file

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message

      ! Local variables:
!   CHARACTER(1024)                          :: TmpMsg                              ! a temporary message (so I don't have to keep typing the same error message)
   REAL(ReKi)                               :: SmallAngleLimit_Rad                 ! Largest input angle considered "small" (check in input file), radians


      ! Initialize error status and angle limit defined locally (in correct units)

   ErrStat = ErrID_None
   ErrMsg  = ''
   SmallAngleLimit_Rad = SmallAngleLimit_Deg*D2R


      ! note that all angles are assumed to be in radians here:

      ! Check that angles are in the range (-pi, pi] radians (i.e., (-180, 180] degrees ):
      ! NOTE: these are local subroutines, with ErrStat and ErrMsg INTENT(INOUT)

   CALL CheckAngle180Range( InputFileData%RotFurl,  'RotFurl',  ErrStat, ErrMsg )
   CALL CheckAngle180Range( InputFileData%TailFurl, 'TailFurl', ErrStat, ErrMsg )
   CALL CheckAngle180Range( InputFileData%TFinSkew, 'TFinSkew', ErrStat, ErrMsg )
   CALL CheckAngle180Range( InputFileData%TFinBank, 'TFinBank', ErrStat, ErrMsg )
   CALL CheckAngle180Range( InputFileData%RFrlSkew, 'RFrlSkew', ErrStat, ErrMsg )
   CALL CheckAngle180Range( InputFileData%TFrlSkew, 'TFrlSkew', ErrStat, ErrMsg )

   CALL CheckAngle180Range( InputFileData%RFrlUSSP, 'RFrlUSSP', ErrStat, ErrMsg )
   CALL CheckAngle180Range( InputFileData%TFrlUSSP, 'TFrlUSSP', ErrStat, ErrMsg )
   CALL CheckAngle180Range( InputFileData%RFrlUSDP, 'RFrlUSDP', ErrStat, ErrMsg )
   CALL CheckAngle180Range( InputFileData%TFrlUSDP, 'TFrlUSDP', ErrStat, ErrMsg )

   CALL CheckAngle180Range( InputFileData%RFrlDSSP, 'RFrlDSSP', ErrStat, ErrMsg )
   IF ( InputFileData%RFrlDSSP > InputFileData%RFrlUSSP ) THEN
      CALL SetErrors( ErrID_Fatal,'RFrlDSSP must not be larger than RFrlUSSP.')
   END IF

   CALL CheckAngle180Range( InputFileData%RFrlDSDP, 'RFrlDSDP', ErrStat, ErrMsg )
   IF ( InputFileData%RFrlDSDP > InputFileData%RFrlUSDP ) THEN
      CALL SetErrors( ErrID_Fatal,'RFrlDSDP must not be larger than RFrlUSDP.')
   END IF

   CALL CheckAngle180Range( InputFileData%TFrlDSSP, 'TFrlDSSP', ErrStat, ErrMsg )
   IF ( InputFileData%TFrlDSSP > InputFileData%TFrlUSSP ) THEN
      CALL SetErrors( ErrID_Fatal,'TFrlDSSP must not be larger than TFrlUSSP.')
   END IF

   CALL CheckAngle180Range( InputFileData%TFrlDSDP, 'TFrlDSDP', ErrStat, ErrMsg )
   IF ( InputFileData%TFrlDSDP > InputFileData%TFrlUSDP ) THEN
      CALL SetErrors( ErrID_Fatal,'TFrlDSDP must not be larger than TFrlUSDP.')
   END IF


      ! Check that tilt angles are in the range [-pi/2, pi/2] radians (i.e., [-90, 90] degrees ):

   CALL CheckAngle90Range( InputFileData%TFinTilt, 'TFinTilt', ErrStat, ErrMsg )
   CALL CheckAngle90Range( InputFileData%RFrlTilt, 'RFrlTilt', ErrStat, ErrMsg )
   CALL CheckAngle90Range( InputFileData%TFrlTilt, 'TFrlTilt', ErrStat, ErrMsg )


      ! Check for violations of the small-angle assumption (15-degree limit, using radians):

   IF ( ABS( InputFileData%ShftSkew ) > SmallAngleLimit_Rad )  THEN
      CALL SetErrors( ErrID_Fatal,'ShftSkew should only be used to skew the shaft a few degrees away from the zero-yaw ' &
                //'position and must not be used as a replacement for the yaw angle. '&
                //'ShftSkew must be between -'//TRIM(Num2LStr(SmallAngleLimit_Rad))//' and ' &
                                              //TRIM(Num2LStr(SmallAngleLimit_Rad))//' radians.' )
   END IF


      ! Warn if tail is defined upwind of the tower:

   IF ( InputFileData%BoomCMxn < 0.0_ReKi )  THEN   ! Print out warning when tail boom CM defined upwind of the tower.
      CALL SetErrors( ErrID_Warn,'WARNING: Tail boom CM is defined upwind of the tower (BoomCMxn < 0).')
   ENDIF

   IF ( InputFileData%TFinCMxn < 0.0_ReKi )  THEN   ! Print out warning when tail fin CM defined upwind of the tower.
      CALL SetErrors( ErrID_Warn,'WARNING: Tail fin CM is defined upwind of the tower (TFinCMxn < 0).')
   ENDIF

   IF ( InputFileData%TFinCPxn < 0.0_ReKi )  THEN   ! Print out warning when tail fin CP defined upwind of the tower.
      CALL SetErrors( ErrID_Warn,'WARNING: Tail fin CP is defined upwind of the tower (TFinCPxn < 0).')
   ENDIF


      ! Check that mass, inertias, damping, etc. values aren't negative:

   CALL ErrIfNegative( InputFileData%RFrlMass, 'RFrlMass', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%BoomMass, 'BoomMass', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%TFinMass, 'TFinMass', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%RFrlIner, 'RFrlIner', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%TFrlIner, 'TFrlIner', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%RFrlSpr,  'RFrlSpr',  ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%RFrlDmp,  'RFrlDmp',  ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%RFrlCDmp, 'RFrlCDmp', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%RFrlUSSpr,'RFrlUSSpr',ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%RFrlDSSpr,'RFrlDSSpr',ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%RFrlUSDmp,'RFrlUSDmp',ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%RFrlDSDmp,'RFrlDSDmp',ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%TFrlSpr,  'TFrlSpr',  ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%TFrlDmp,  'TFrlDmp',  ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%TFrlCDmp, 'TFrlCDmp', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%TFrlUSSpr,'TFrlUSSpr',ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%TFrlDSSpr,'TFrlDSSpr',ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%TFrlUSDmp,'TFrlUSDmp',ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%TFrlDSDmp,'TFrlDSDmp',ErrStat, ErrMsg )


      ! Check that furling models are valid:

   IF ( InputFileData%TFrlMod < 0 .OR. InputFileData%TFrlMod > 2 )  THEN
      CALL SetErrors( ErrID_Fatal, 'TFrlMod must be 0, 1, or 2.')
   END IF

   IF ( InputFileData%RFrlMod < 0 .OR. InputFileData%RFrlMod > 2 )  THEN
      CALL SetErrors( ErrID_Fatal, 'RFrlMod must be 0, 1, or 2.' )
   END IF


   !   ! bjj: THESE ARE checks for tail fin aerodynamics, which should be in aerodyn, in my opinion
   !CALL ErrIfNegative( TFinArea,               'TFinArea',ErrStat, ErrMsg )
   !
   !IF ( TFinMod < 0 .OR. TFinMod > 2 )  THEN
   !   CALL SetErrors( ErrID_Fatal,'TFinMod must be 0, 1, or 2.')
   !END IF


   RETURN

CONTAINS
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE SetErrors( ErrStat3, ErrMsg3 )
   ! This routine sets the error message and flag when an error has occurred
   !...............................................................................................................................
   INTEGER(IntKi), INTENT(IN) :: ErrStat3     ! Error status for this error
   CHARACTER(*),   INTENT(IN) :: ErrMsg3      ! Error message for this error

      ErrStat = MAX( ErrStat, ErrStat3 )
      IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg  = TRIM(ErrMsg)//TRIM(ErrMsg3)

   END SUBROUTINE SetErrors
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE CheckAngle180Range( Var, VarDesc, ErrStat, ErrMsg )
   ! This routine checks that an angle is in the range (-pi, pi] radians. If not, ErrStat = ErrID_Fatal
   ! Note that all values are assumed to be in radians, even if read in degrees (-180 deg, 180 deg]
   !...............................................................................................................................
   REAL(ReKi),     INTENT(IN)    :: Var         ! Variable to check
   CHARACTER(*),   INTENT(IN)    :: VarDesc     ! Description of variable (used in error message)
   INTEGER(IntKi), INTENT(INOUT) :: ErrStat     ! Error status to update if Var is not in specified range
   CHARACTER(*),   INTENT(INOUT) :: ErrMsg      ! Error message to update if Var is not in specified range


      IF ( ( Var <= -pi ) .OR. ( Var > pi ) )  THEN
         ErrStat = ErrID_Fatal
         IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg  = TRIM(ErrMsg)// &
                   TRIM(VarDesc)//' must be greater than -pi radians and less than or equal to pi radians '// &
                                        '(i.e., in the range (-180, 180] degrees).'
      END IF

   END SUBROUTINE CheckAngle180Range
   !-------------------------------------------------------------------------------------------------------------------------------
 END SUBROUTINE ValidateFurlData
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetFurlParameters( p, InputFileData, ErrStat, ErrMsg  )
! This takes the furling input file data and sets the corresponding furling parameters.
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed variables

   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                            ! Parameters of the structural dynamics module
   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                ! Data stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       ! Error message

      ! Local variables:

   REAL(ReKi)                               :: x                            ! Fractional location between two points in linear interpolation
!   INTEGER(IntKi )                          :: J                            ! Index for the node arrays
!   INTEGER(IntKi)                           :: InterpInd                    ! Index for the interpolation routine


      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''


      ! Direct copy of InputFileData to parameters

   p%RFrlMass = InputFileData%RFrlMass
   p%BoomMass = InputFileData%BoomMass
   p%TFinMass = InputFileData%TFinMass
   p%TFrlIner = InputFileData%TFrlIner

   p%RFrlMod  = InputFileData%RFrlMod
   p%TFrlMod  = InputFileData%TFrlMod

   p%RFrlSpr  = InputFileData%RFrlSpr
   p%RFrlDmp  = InputFileData%RFrlDmp
   p%RFrlCDmp = InputFileData%RFrlCDmp
   p%RFrlUSSP = InputFileData%RFrlUSSP
   p%RFrlDSSP = InputFileData%RFrlDSSP
   p%RFrlDSSpr= InputFileData%RFrlDSSpr
   p%RFrlUSSpr= InputFileData%RFrlUSSpr
   p%RFrlUSDP = InputFileData%RFrlUSDP
   p%RFrlDSDP = InputFileData%RFrlDSDP
   p%RFrlUSDmp= InputFileData%RFrlUSDmp
   p%RFrlDSDmp= InputFileData%RFrlDSDmp

   p%TFrlSpr  = InputFileData%TFrlSpr
   p%TFrlDmp  = InputFileData%TFrlDmp
   p%TFrlCDmp = InputFileData%TFrlCDmp
   p%TFrlUSSP = InputFileData%TFrlUSSP
   p%TFrlDSSP = InputFileData%TFrlDSSP
   p%TFrlUSSpr= InputFileData%TFrlUSSpr
   p%TFrlDSSpr= InputFileData%TFrlDSSpr
   p%TFrlUSDP = InputFileData%TFrlUSDP
   p%TFrlDSDP = InputFileData%TFrlDSDP
   p%TFrlUSDmp= InputFileData%TFrlUSDmp
   p%TFrlDSDmp= InputFileData%TFrlDSDmp

   p%RFrlPntxn = InputFileData%RFrlPntxn
   p%RFrlPntyn = InputFileData%RFrlPntyn
   p%RFrlPntzn = InputFileData%RFrlPntzn

   p%TFrlPntxn = InputFileData%TFrlPntxn
   p%TFrlPntyn = InputFileData%TFrlPntyn
   p%TFrlPntzn = InputFileData%TFrlPntzn


      ! Store sine/cosine values instead of some input angles:

   p%CShftSkew = COS( InputFileData%ShftSkew )
   p%SShftSkew = SIN( InputFileData%ShftSkew )

   p%CTFinSkew = COS( InputFileData%TFinSkew )
   p%STFinSkew = SIN( InputFileData%TFinSkew )
   p%CTFinTilt = COS( InputFileData%TFinTilt )
   p%STFinTilt = SIN( InputFileData%TFinTilt )
   p%CTFinBank = COS( InputFileData%TFinBank )
   p%STFinBank = SIN( InputFileData%TFinBank )

   p%CRFrlSkew = COS( InputFileData%RFrlSkew )
   p%SRFrlSkew = SIN( InputFileData%RFrlSkew )
   p%CRFrlTilt = COS( InputFileData%RFrlTilt )
   p%SRFrlTilt = SIN( InputFileData%RFrlTilt )

   p%CTFrlSkew = COS( InputFileData%TFrlSkew )
   p%STFrlSkew = SIN( InputFileData%TFrlSkew )
   p%CTFrlTilt = COS( InputFileData%TFrlTilt )
   p%STFrlTilt = SIN( InputFileData%TFrlTilt )


      ! Common multiplications of sines and cosines:

   p%CRFrlSkw2 = p%CRFrlSkew**2
   p%SRFrlSkw2 = p%SRFrlSkew**2
   p%CSRFrlSkw = p%CRFrlSkew*p%SRFrlSkew
   p%CRFrlTlt2 = p%CRFrlTilt**2
   p%SRFrlTlt2 = p%SRFrlTilt**2
   p%CSRFrlTlt = p%CRFrlTilt*p%SRFrlTilt

   p%CTFrlSkw2 = p%CTFrlSkew**2
   p%STFrlSkw2 = p%STFrlSkew**2
   p%CSTFrlSkw = p%CTFrlSkew*p%STfrlSkew
   p%CTFrlTlt2 = p%CTFrlTilt**2
   p%STFrlTlt2 = p%STFrlTilt**2
   p%CSTFrlTlt = p%CTFrlTilt*p%STFrlTilt


      ! Calculate some positions:

   p%rWIxn     = InputFileData%BoomCMxn - p%TFrlPntxn
   p%rWIyn     = InputFileData%BoomCMyn - p%TFrlPntyn
   p%rWIzn     = InputFileData%BoomCMzn - p%TFrlPntzn

   p%rWJxn     = InputFileData%TFinCMxn - p%TFrlPntxn
   p%rWJyn     = InputFileData%TFinCMyn - p%TFrlPntyn
   p%rWJzn     = InputFileData%TFinCMzn - p%TFrlPntzn

   p%rWKxn     = InputFileData%TFinCPxn - p%TFrlPntxn
   p%rWKyn     = InputFileData%TFinCPyn - p%TFrlPntyn
   p%rWKzn     = InputFileData%TFinCPzn - p%TFrlPntzn

   p%rVDxn     = InputFileData%RFrlCMxn - p%RFrlPntxn
   p%rVDyn     = InputFileData%RFrlCMyn - p%RFrlPntyn
   p%rVDzn     = InputFileData%RFrlCMzn - p%RFrlPntzn

   p%rVPxn     =        0.0_ReKi        - p%RFrlPntxn
   p%rVPyn     = InputFileData%Yaw2Shft - p%RFrlPntyn


      ! Note: These positions are also used for non-furling machines:

   p%rVPzn     = InputFileData%Twr2Shft - p%RFrlPntzn
   p%rVIMUxn   = InputFileData%NcIMUxn  - p%RFrlPntxn
   p%rVIMUyn   = InputFileData%NcIMUyn  - p%RFrlPntyn
   p%rVIMUzn   = InputFileData%NcIMUzn  - p%RFrlPntzn

END SUBROUTINE SetFurlParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetPrimaryParameters( p, InputFileData, ErrStat, ErrMsg  )
! This takes the primary input file data and sets the corresponding parameters.
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed variables

   TYPE(ED_ParameterType),   INTENT(INOUT)  :: p                            ! Parameters of the structural dynamics module
   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                ! Data stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       ! Error message

!bjj: ERROR CHECKING!!!

      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''

   !p%Twr2Shft  = InputFileData%Twr2Shft
   !p%HubIner   = InputFileData%HubIner
   !p%NacYIner  = InputFileData%NacYIner


   !...............................................................................................................................
   ! Direct copy of variables:
   !...............................................................................................................................
   p%NumBl     = InputFileData%NumBl
   p%TipRad    = InputFileData%TipRad
   p%HubRad    = InputFileData%HubRad
   p%method    = InputFileData%method
   p%TwrNodes  = InputFileData%TwrNodes

   p%PtfmCMxt = InputFileData%PtfmCMxt
   p%PtfmCMyt = InputFileData%PtfmCMyt   
   
   p%DT        = InputFileData%DT
   p%Gravity   = InputFileData%Gravity
   p%OverHang  = InputFileData%OverHang
   p%ShftGagL  = InputFileData%ShftGagL
   p%TowerHt   = InputFileData%TowerHt
   p%TowerBsHt = InputFileData%TowerBsHt
   p%PtfmRefzt = InputFileData%PtfmRefzt
   
   p%HubMass   = InputFileData%HubMass
   p%GenIner   = InputFileData%GenIner
   p%NacMass   = InputFileData%NacMass
   p%YawBrMass = InputFileData%YawBrMass
   p%PtfmMass  = InputFileData%PtfmMass
   p%PtfmRIner = InputFileData%PtfmRIner
   p%PtfmPIner = InputFileData%PtfmPIner
   p%PtfmYIner = InputFileData%PtfmYIner
   p%GBoxEff   = InputFileData%GBoxEff
   p%GBRatio   = InputFileData%GBRatio
   p%DTTorSpr  = InputFileData%DTTorSpr
   p%DTTorDmp  = InputFileData%DTTorDmp


   p%NTwGages  = InputFileData%NTwGages
   p%TwrGagNd  = InputFileData%TwrGagNd
   p%NBlGages  = InputFileData%NBlGages
   p%BldGagNd  = InputFileData%BldGagNd
   !p%OutFile   = InputFileData%OutFile
   !p%OutFileFmt= InputFileData%OutFileFmt !wrbinoutput, wrtxtoutput???
   p%OutFmt    = InputFileData%OutFmt
   p%Tstart    = InputFileData%Tstart
   !p%DecFact   = InputFileData%DecFact
   p%NumOuts   = InputFileData%NumOuts

   IF ( p%NumBl == 2 ) THEN
      p%UndSling = InputFileData%UndSling
      p%TeetMod  = InputFileData%TeetMod
      p%TeetDmpP = InputFileData%TeetDmpP
      p%TeetDmp  = InputFileData%TeetDmp
      p%TeetCDmp = InputFileData%TeetCDmp
      p%TeetSStP = InputFileData%TeetSStP
      p%TeetHStP = InputFileData%TeetHStP
      p%TeetSSSp = InputFileData%TeetSSSp
      p%TeetHSSp = InputFileData%TeetHSSp
   ELSE ! Three-bladed turbines don't use these parameters, so set them to zero.
      p%UndSling = 0.0
      p%TeetMod  = 0
      p%TeetDmpP = 0.0
      p%TeetDmp  = 0.0
      p%TeetCDmp = 0.0
      p%TeetSStP = 0.0
      p%TeetHStP = 0.0
      p%TeetSSSp = 0.0
      p%TeetHSSp = 0.0
   END IF

   
   CALL AllocAry( p%TipMass, p%NumBl, 'TipMass', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   p%TipMass   = InputFileData%TipMass

      ! initialize all of the DOF parameters:
   CALL Init_DOFparameters( InputFileData, p, ErrStat, ErrMsg ) !sets p%NDOF and p%NAug
      IF (ErrStat >= AbortErrLev) RETURN

      ! Set parameters for output channels:
   CALL SetOutParam(InputFileData%OutList, p, ErrStat, ErrMsg ) ! requires: p%NumOuts, p%NumBl, p%NBlGages, p%NTwGages; sets: p%OutParam.
      IF (ErrStat >= AbortErrLev) RETURN

   IF ( InputFileData%TabDelim ) THEN
      p%Delim = TAB
   ELSE
      p%Delim = ' '
   !ELSE
   !   p%Delim = ','
   END IF

   !...............................................................................................................................
   ! Calculate some indirect inputs:
   !...............................................................................................................................
   p%TwoPiNB   = TwoPi_D/p%NumBl                                                   ! 2*Pi/NumBl is used in RtHS().

   p%rZT0zt    = p%TowerBsHt - p%PtfmRefzt                                         ! zt-component of position vector rZT0.
   p%RefTwrHt  = p%TowerHt   - p%PtfmRefzt                                         ! Vertical distance between ElastoDyn's undisplaced tower height (variable TowerHt) and ElastoDyn's inertia frame reference point (variable PtfmRef).
   p%TwrFlexL  = p%TowerHt   - p%TowerBsHt                                         ! Height / length of the flexible portion of the tower.
   p%BldFlexL  = p%TipRad    - p%HubRad                                            ! Length of the flexible portion of the blade.
   if (p%BD4Blades) p%BldFlexL = 0.0_ReKi
   
   p%rZYzt     = InputFileData%PtfmCMzt - p%PtfmRefzt

   !...............................................................................................................................
   ! set cosine and sine of Precone and Delta3 angles:
   !...............................................................................................................................
   CALL AllocAry( p%CosPreC,  p%NumBl,                              'CosPreC',   ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry( p%SinPreC,  p%NumBl,                              'SinPreC',   ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   p%CosPreC  = COS( InputFileData%Precone(1:p%NumBl) )
   p%SinPreC  = SIN( InputFileData%Precone(1:p%NumBl) )
   p%CosDel3  = COS( InputFileData%Delta3 )
   p%SinDel3  = SIN( InputFileData%Delta3 )

   !...............................................................................................................................

      ! Calculate the average tip radius normal to the shaft (AvgNrmTpRd)
      !   and the swept area of the rotor (ProjArea):

   p%AvgNrmTpRd = p%TipRad*SUM(p%CosPreC)/p%NumBl     ! Average tip radius normal to the saft.
   p%ProjArea   = pi*( p%AvgNrmTpRd**2 )              ! Swept area of the rotor projected onto the rotor plane (the plane normal to the low-speed shaft).

   p%RotSpeed  = InputFileData%RotSpeed               ! Rotor speed in rad/sec.
   p%CShftTilt = COS( InputFileData%ShftTilt )
   p%SShftTilt = SIN( InputFileData%ShftTilt )

   p%HubHt     = p%TowerHt + InputFileData%Twr2Shft + p%OverHang*p%SShftTilt


      ! Direct copy of InputFileData to parameters

   !p%FlapDOF1  = InputFileData%FlapDOF1
   !p%FlapDOF2  = InputFileData%FlapDOF2
   !p%EdgeDOF   = InputFileData%EdgeDOF
   !p%TeetDOF   = InputFileData%TeetDOF
   !p%DrTrDOF   = InputFileData%DrTrDOF
   !p%GenDOF    = InputFileData%GenDOF
   !p%YawDOF    = InputFileData%YawDOF
   !p%TwFADOF1  = InputFileData%TwFADOF1
   !p%TwFADOF2  = InputFileData%TwFADOF2
   !p%TwSSDOF1  = InputFileData%TwSSDOF1
   !p%TwSSDOF2  = InputFileData%TwSSDOF2
   !p%PtfmSgDOF = InputFileData%PtfmSgDOF
   !p%PtfmSwDOF = InputFileData%PtfmSwDOF
   !p%PtfmHvDOF = InputFileData%PtfmHvDOF
   !p%PtfmRDOF  = InputFileData%PtfmRDOF
   !p%PtfmPDOF  = InputFileData%PtfmPDOF
   !p%PtfmYDOF  = InputFileData%PtfmYDOF
   !p%Azimuth   = InputFileData%Azimuth
   p%RotSpeed  = InputFileData%RotSpeed
   !p%TTDspFA   = InputFileData%TTDspFA
   !p%TTDspSS   = InputFileData%TTDspSS
   !p%PtfmSurge = InputFileData%PtfmSurge
   !p%PtfmSway  = InputFileData%PtfmSway
   !p%PtfmHeave = InputFileData%PtfmHeave
   !p%PtfmRoll  = InputFileData%PtfmRoll
   !p%PtfmPitch = InputFileData%PtfmPitch
   !p%PtfmYaw   = InputFileData%PtfmYaw
   p%HubCM     = InputFileData%HubCM
   p%AzimB1Up  = InputFileData%AzimB1Up

   p%NacCMxn   = InputFileData%NacCMxn
   p%NacCMyn   = InputFileData%NacCMyn
   p%NacCMzn   = InputFileData%NacCMzn
   !p%NcIMUxn   = InputFileData%NcIMUxn
   !p%NcIMUyn   = InputFileData%NcIMUyn
   !p%NcIMUzn   = InputFileData%NcIMUzn


   ! plus everything else from FAST_Initialize



END SUBROUTINE SetPrimaryParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadBladeInputs ( BldFile, MeshFile, ReadAdmVals, InputFileData, UnEc, ErrStat, ErrMsg )
! This routine reads the data from the blade and mesh inputs files.
! This routines assumes that InputFileData%NumBl has already been set.
!..................................................................................................................................


   IMPLICIT                        NONE


      ! Passed variables:

!   TYPE(ED_ParameterType), INTENT(INOUT)  :: p                                   ! Parameters of the structural dynamics module
   TYPE(ED_InputFile),     INTENT(INOUT)  :: InputFileData                       ! Input file data Data for Blade K stored in the module's input file
   CHARACTER(*),           INTENT(IN)     :: BldFile(:)                          ! The array of file names containing blade information
   CHARACTER(*),           INTENT(IN)     :: MeshFile                            ! The file names containing blade mesh information (for now, the aerodyn primary file)
   INTEGER(IntKi),         INTENT(IN)     :: UnEc                                ! I/O unit for echo file. If present and > 0, write to UnEc

   INTEGER(IntKi),         INTENT(OUT)    :: ErrStat                             ! The error ID
   CHARACTER(*),           INTENT(OUT)    :: ErrMsg                              ! Message describing error
   LOGICAL,                INTENT(IN)     :: ReadAdmVals                         ! Logical to determine if Adams inputs should be read from file


      ! Local variables:
   INTEGER(IntKi)                         :: K                                   ! Blade number
   INTEGER(IntKi)                         :: ErrStat2                            ! Temporary error ID
   LOGICAL                                :: ReadFile                            ! determines if an input file for a blade is the same as the file for the previous blade
   CHARACTER(ErrMsgLen)                   :: ErrMsg2                             ! Temporary message describing error


      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''


      ! Allocate space for the input file data
   ALLOCATE( InputFileData%InpBlMesh( 1_IntKi ), STAT=ErrStat2 )              ! for now, we're assuming the discretization is the same on all blades
   IF ( ErrStat2 /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating InpBl array'
      RETURN
   END IF

   ALLOCATE( InputFileData%InpBl( InputFileData%NumBl ), STAT=ErrStat2 )
   IF ( ErrStat2 /= 0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = 'Error allocating InpBl array'
      RETURN
   END IF



      ! Get the blade discretization here:
   IF ( len_trim(MeshFile) == 0 ) THEN
      InputFileData%InpBlMesh(1)%BldNodes = InputFileData%BldNodes
   ELSE
         ! we will get the discretization from AeroDyn's input file
      CALL ReadBladeMeshFileAD( InputFileData%InpBlMesh(1), MeshFile, UnEc, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
   END IF


      ! Read the input file(s) for all of the blades:
   ReadFile = .TRUE.
   DO K = 1,InputFileData%NumBl

      IF ( ReadFile ) THEN

            ! Add a separator to the echo file if appropriate.

         IF ( UnEc > 0 )  THEN
            WRITE (UnEc,'(//,A,/)')  'Blade '//TRIM( Num2LStr( K ) )//' input data from file "'//TRIM( BldFile(K) )//'":'
         END IF

         CALL ReadBladeFile( BldFile(K), InputFileData%InpBl(K), ReadAdmVals, UnEc, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,'Errors reading blade '//TRIM(Num2LStr(K))//' input file: '//TRIM(ErrMsg2))
            IF ( ErrStat >= AbortErrLev ) RETURN

      ELSE
         CALL ED_CopyBladeInputData( InputFileData%InpBl(K-1), InputFileData%InpBl(K), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,'Errors copying blade '//TRIM(Num2LStr(K-1))//' input file data: '//TRIM(ErrMsg2))
            IF ( ErrStat >= AbortErrLev ) RETURN
               ! bjj: we could just read the file again...

      END IF

         ! If the next file is the same as this one, don't read it again:

      IF ( K /= InputFileData%NumBl ) ReadFile = BldFile(K) /= BldFile( K + 1 )

   END DO


   RETURN

CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      CALL SetErrStat(ErrID,Msg,ErrStat,ErrMsg,'ReadBladeInputs')

         !.........................................................................................................................
         !! Clean up if we're going to return on error: close file, deallocate local arrays
         !!.........................................................................................................................
         !IF ( ErrStat >= AbortErrLev ) THEN
         !
         !END IF


   END SUBROUTINE CheckError

END SUBROUTINE ReadBladeInputs
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadBladeFile ( BldFile, BladeKInputFileData, ReadAdmVals, UnEc, ErrStat, ErrMsg )
! This routine reads a blade input file.
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed variables:

   TYPE(BladeInputData),     INTENT(INOUT)  :: BladeKInputFileData                 ! Data for Blade K stored in the module's input file
   CHARACTER(*),             INTENT(IN)     :: BldFile                             ! Name of the blade input file data
   LOGICAL,                  INTENT(IN)     :: ReadAdmVals                         ! Logical to determine if Adams inputs should be read from file
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                ! I/O unit for echo file. If present and > 0, write to UnEc

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message


      ! Local variables:

   REAL(ReKi)                   :: AdjBlMs                                         ! Factor to adjust blade mass density.
   REAL(ReKi)                   :: AdjEdSt                                         ! Factor to adjust edge stiffness.
   REAL(ReKi)                   :: AdjFlSt                                         ! Factor to adjust flap stiffness.

   REAL(ReKi)                   :: TmpRAry(17)                                     ! Temporary variable to read table from file (up to 17 columns)

   INTEGER(IntKi)               :: I                                               ! A generic DO index.
   INTEGER( IntKi )             :: UnIn                                            ! Unit number for reading file
   INTEGER( IntKi )             :: NInputCols                                      ! Number of columns to be read from the file
   INTEGER(IntKi)               :: ErrStat2                                        ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                         ! Temporary Err msg


   UnIn = -1
   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the input file for blade K.

   CALL OpenFInpFile ( UnIn, BldFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- HEADER -------------------------------------------------------

      ! Skip the header.

   CALL ReadCom ( UnIn, BldFile, 'unused blade file header line 1', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn, BldFile, 'unused blade file header line 2', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- BLADE PARAMETERS ---------------------------------------------

      ! Skip the comment line.

   CALL ReadCom ( UnIn, BldFile, 'blade parameters', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NBlInpSt - Number of blade input stations.

   CALL ReadVar ( UnIn, BldFile, BladeKInputFileData%NBlInpSt, 'NBlInpSt', 'Number of blade input stations', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! .......... Allocate the arrays based on this NBlInpSt input ..........
   CALL Alloc_BladeInputProperties( BladeKInputFileData, ReadAdmVals, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! BldFlDmp - Blade structural damping ratios in flapwise direction.

   CALL ReadAryLines( UnIn, BldFile, BladeKInputFileData%BldFlDmp, SIZE(BladeKInputFileData%BldFlDmp), 'BldFlDmp', &
                                       'Blade structural damping ratios in flapwise direction', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! BldEdDmp - Blade structural damping ratios in edgewise direction.

   CALL ReadAryLines( UnIn, BldFile, BladeKInputFileData%BldEdDmp, SIZE(BladeKInputFileData%BldEdDmp), 'BldEdDmp', &
                                       'Blade structural damping ratios in edgewise direction', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !  -------------- BLADE ADJUSTMENT FACTORS -------------------------------------


      ! Skip the comment line.

   CALL ReadCom ( UnIn, BldFile, 'blade adjustment factors', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! FlStTunr(1) - Blade flapwise modal stiffness tuners.

   CALL ReadAryLines ( UnIn, BldFile, BladeKInputFileData%FlStTunr, SIZE(BladeKInputFileData%FlStTunr), 'FlStTunr', &
                                                  'Blade flapwise modal stiffness tuners', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! AdjBlMs - Factor to adjust blade mass density.

   CALL ReadVar ( UnIn, BldFile, AdjBlMs, 'AdjBlMs', 'Factor to adjust blade mass density', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! AdjFlSt - Factor to adjust blade flap stiffness.

   CALL ReadVar ( UnIn, BldFile, AdjFlSt, 'AdjFlSt', 'Factor to adjust blade flap stiffness', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! AdjEdSt - Factor to adjust blade edge stiffness.

   CALL ReadVar ( UnIn, BldFile, AdjEdSt, 'AdjEdSt', 'Factor to adjust blade edge stiffness', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



         ! Check the locally-defined adjustment factors: AdjBlMs, AdjFlSt, AdjEdSt

      IF ( AdjBlMs <= 0.0_ReKi ) THEN
         CALL CheckError( ErrID_Warn, ' AdjBlMs must be greater than zero.' )
         IF ( ErrStat >= AbortErrLev ) RETURN
      END IF

      IF ( AdjFlSt <= 0.0_ReKi ) THEN
         CALL CheckError( ErrID_Warn, ' AdjFlSt must be greater than zero.' )
         IF ( ErrStat >= AbortErrLev ) RETURN
      END IF

      IF ( AdjEdSt <= 0.0_ReKi ) THEN
         CALL CheckError( ErrID_Warn, ' AdjEdSt must be greater than zero.' )
         IF ( ErrStat >= AbortErrLev ) RETURN
      END IF


   !  -------------- DISTRIBUTED BLADE PROPERTIES ---------------------------------


      ! Skip the comment lines.

   CALL ReadCom ( UnIn, BldFile, 'distributed blade parameters'     , ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn, BldFile, 'distributed-blade-parameter names', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn, BldFile, 'distributed-blade-parameter units', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! Read the table.

   IF ( ReadAdmVals ) THEN
      NInputCols = 17
   ELSE
      NInputCols = 6
   END IF


   DO I=1,BladeKInputFileData%NBlInpSt

      CALL ReadAry( UnIn, BldFile, TmpRAry, NInputCols, 'Line'//TRIM(Num2LStr(I)), 'Blade input station table', &
                    ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

      BladeKInputFileData%BlFract( I) = TmpRAry(1)
      BladeKInputFileData%PitchAx( I) = TmpRAry(2)
      BladeKInputFileData%StrcTwst(I) = TmpRAry(3)*D2R      ! Input in degrees; converted to radians here
      BladeKInputFileData%BMassDen(I) = TmpRAry(4)*AdjBlMs  ! Apply the correction factors to the elemental data.
      BladeKInputFileData%FlpStff( I) = TmpRAry(5)*AdjFlSt  ! Apply the correction factors to the elemental data.
      BladeKInputFileData%EdgStff( I) = TmpRAry(6)*AdjEdSt  ! Apply the correction factors to the elemental data.

      IF ( NInputCols > 6 ) THEN
         BladeKInputFileData%GJStff(   I) = TmpRAry( 7)
         BladeKInputFileData%EAStff(   I) = TmpRAry( 8)
         BladeKInputFileData%Alpha(    I) = TmpRAry( 9)
         BladeKInputFileData%FlpIner(  I) = TmpRAry(10)
         BladeKInputFileData%EdgIner(  I) = TmpRAry(11)
         BladeKInputFileData%PrecrvRef(I) = TmpRAry(12)
         BladeKInputFileData%PreswpRef(I) = TmpRAry(13)
         BladeKInputFileData%FlpcgOf(  I) = TmpRAry(14)
         BladeKInputFileData%EdgcgOf(  I) = TmpRAry(15)
         BladeKInputFileData%FlpEAOf(  I) = TmpRAry(16)
         BladeKInputFileData%EdgEAOf(  I) = TmpRAry(17)
      END IF
   ENDDO ! I



   !  -------------- BLADE MODE SHAPES --------------------------------------------


      ! Skip the comment line.

   CALL ReadCom ( UnIn, BldFile, 'blade mode shapes', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! BldFl1Sh - Blade-flap mode-1 shape coefficients.
   CALL ReadAryLines ( UnIn, BldFile, BladeKInputFileData%BldFl1Sh, SIZE(BladeKInputFileData%BldFl1Sh), 'BldFl1Sh', &
                           'Blade-flap mode-1 shape coefficients', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! BldFl2Sh - Blade-flap mode-2 shape coefficients.

   CALL ReadAryLines ( UnIn, BldFile, BladeKInputFileData%BldFl2Sh, SIZE(BladeKInputFileData%BldFl2Sh), 'BldFl2Sh', &
                    'Blade-flap mode-2 shape coefficients', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! BldEdgSh - Blade-edge mode shape coefficients.

   CALL ReadAryLines ( UnIn, BldFile, BladeKInputFileData%BldEdgSh, SIZE(BladeKInputFileData%BldEdgSh), 'BldEdgSh', &
                     'Blade-edge mode shape coefficients', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



   !  -------------- END OF FILE --------------------------------------------

      ! Close the blade file.

   CLOSE ( UnIn )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadBladeFile:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            IF (UnIn > 0) CLOSE( UnIn )
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE ReadBladeFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadBladeMeshFileAD( BladeKInputFileMesh, MeshFile, UnEc, ErrStat, ErrMsg )
! This routine reads in the AeroDyn v14.00.00 input file to get the
!   blade discretization used in the structural dynamics module.
!..................................................................................................................................


   IMPLICIT                        NONE

      ! Passed variables

   TYPE(ED_BladeMeshInputData),   INTENT(INOUT)  :: BladeKInputFileMesh                 ! All the data in the ElastoDyn input file
   CHARACTER(*),                  INTENT(IN)     :: MeshFile                            ! Name of the AeroDyn input file data (for mesh)

   INTEGER(IntKi),                INTENT(IN)     :: UnEc                                ! I/O unit for echo file. If present and > 0, write to UnEc
   INTEGER(IntKi),                INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),                  INTENT(OUT)    :: ErrMsg                              ! Error message

      ! Local variables:
   INTEGER(IntKi), PARAMETER    :: NInputCols = 4                                       ! Number of input columns to be read from the file
   REAL(ReKi)                   :: TmpRAry(NInputCols)                                  ! Temporary variable to read table from file
   INTEGER(IntKi)               :: I                                                    ! loop counter
   INTEGER(IntKi)               :: NumLin2Skp                                           ! number of lines to read
   INTEGER(IntKi)               :: NumFoil                                              ! number of airfoil lines to skip in the AD input file.
   INTEGER(IntKi)               :: UnIn                                                 ! Unit number for reading file

   INTEGER(IntKi)               :: ErrStat2                                             ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                              ! Temporary Err msg
   CHARACTER(1024)              :: Line                                                 ! Temporary string.
!   CHARACTER(1024)              :: TmpStr(1)                                            ! Temporary string.



      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the AeroDyn input file.

   CALL OpenFInpFile ( UnIn, MeshFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Add a separator to the echo file if appropriate.

   IF ( UnEc > 0 )  WRITE (UnEc,'(//,A,/)')  'Mesh input data from (AeroDyn input) file "'//TRIM( MeshFile )//'":'


   !  -------------- HEADER -------------------------------------------------------
   ! BJJ: This file is AeroDyn's input file. Until we decide on a format for the
   ! structural dynamics input, we will get this information from AeroDyn like we
   ! used to.

   DO I = 1,9
      CALL ReadCom ( UnIn, MeshFile, 'AeroDyn input (for structural dynamics mesh)', ErrStat2, ErrMsg2  )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN
   END DO

      ! See if the next line is "NEWTOWER".  If it is, read 7 more lines.  If not, read 5 more lines.

   CALL ReadVar( UnIn, MeshFile, Line, VarName='NewTowerModel?', VarDescr='Check for tower influence model', ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Check if this is the "special string" to indicate the new tower influence model

   CALL Conv2UC( Line )
   IF ( INDEX(Line, "NEWTOWER" ) > 0 ) THEN
      NumLin2Skp = 7
   ELSE
      NumLin2Skp = 5
   END IF

   DO I = 1,NumLin2Skp
      CALL ReadCom ( UnIn, MeshFile, 'AeroDyn input (for structural dynamics mesh)', ErrStat2, ErrMsg2  )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN
   END DO

   CALL ReadVar ( UnIn, MeshFile, NumFoil, 'NumFoil', &
                  'Number of airfoil lines to skip in AeroDyn input (for structural dynamics mesh)', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   DO I = 1,NumFoil
      CALL ReadCom ( UnIn, MeshFile, 'AeroDyn input (for structural dynamics mesh)', ErrStat2, ErrMsg2  )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN
   END DO


  !  -------------- Blade Mesh Data --------------------------------------------------

      ! Read in the number of blade elements
   CALL ReadVar( UnIn, MeshFile, BladeKInputFileMesh%BldNodes, 'BldNodes', 'Number of blade elements', ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Allocate the arrays to store input
   CALL Alloc_BladeMeshInputProperties( BladeKInputFileMesh, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Read comment line for the element table
   CALL ReadCom( UnIn, MeshFile, 'Blade element table headers', ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   DO I = 1, BladeKInputFileMesh%BldNodes

      CALL ReadAry( UnIn, MeshFile, TmpRAry, NInputCols, 'Blade element line'//TRIM(Num2LStr(I)), 'Blade element input table', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

         BladeKInputFileMesh%RNodes(  I) = TmpRAry(1)
         BladeKInputFileMesh%AeroTwst(I) = TmpRAry(2)*D2R  !Convert input file data (degrees) to radians
         BladeKInputFileMesh%Chord(   I) = TmpRAry(4)

   END DO

      !bjj: move this to a validation routine:
   IF ( ANY( BladeKInputFileMesh%Chord < 0.0_ReKi ) ) THEN
      CALL CheckError( ErrID_Fatal, 'Chord length must be larger than 0 meters.' )
      RETURN
   END IF


      ! Close the input file:

   CLOSE ( UnIn )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadBladeMeshFileAD:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
         END IF

      END IF


   END SUBROUTINE CheckError
   !...............................................................................................................................

END SUBROUTINE ReadBladeMeshFileAD
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadFurlFile( FurlFile, InputFileData, UnEc, ErrStat, ErrMsg  )
! This routine reads the furling file input and converts units as appropriate.
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables:

   TYPE(ED_InputFile),       INTENT(INOUT)  :: InputFileData                       ! All the data in the ElastoDyn input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                ! I/O unit for echo file. If present and > 0, write to UnEc
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message
   CHARACTER(*),             INTENT(IN)     :: FurlFile                            ! Name of the furling input file data

      ! Local variables:

   INTEGER(IntKi)               :: UnIn                                            ! Unit number for reading file
   INTEGER(IntKi)               :: ErrStat2                                        ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                         ! Temporary Err msg


      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the furling input file.

   CALL OpenFInpFile ( UnIn, FurlFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Add a separator to the echo file if appropriate.

   IF ( UnEc > 0 )  WRITE (UnEc,'(//,A,/)')  'Furling input data from file "'//TRIM( FurlFile )//'":'


   !  -------------- FILE HEADER ---------------------------------------------------

   CALL ReadCom ( UnIn, FurlFile, 'unused tower furling header line 1', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn, FurlFile, 'unused tower furling header line 2', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn, FurlFile, 'unused tower furling header line 3', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- Furling FEATURE SWITCHES  --------------------------------------


      ! Skip the comment line.

   CALL ReadCom ( UnIn, FurlFile, 'degree of freedom switches (cont)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! RFrlDOF - Rotor-furl DOF.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlDOF, 'RFrlDOF', 'Rotor-furl DOF (flag)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlDOF - Tail-furl DOF.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlDOF, 'TFrlDOF', 'Tail-furl DOF (flag)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- Furling INITIAL CONDITIONS ------------------------------------


      ! Skip the comment line.

   CALL ReadCom ( UnIn, FurlFile, 'initial conditions (cont)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RotFurl - Initial or fixed rotor-furl angle (read in degrees, converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RotFurl, 'RotFurl', 'Initial or fixed rotor-furl angle (deg)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%RotFurl   = InputFileData%RotFurl*D2R


      ! TailFurl - Initial or fixed tail-furl angle (read in degrees, converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TailFurl, 'TailFurl', 'Initial or fixed tail-furl angle (deg)',  &
                  ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%TailFurl  = InputFileData%TailFurl*D2R


   !  -------------- TURBINE CONFIGURATION (CONT) ---------------------------------


      ! Skip the comment line.

   CALL ReadCom ( UnIn, FurlFile, 'turbine configuration (cont)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Yaw2Shft - Lateral distance from yaw axis to rotor shaft.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%Yaw2Shft, 'Yaw2Shft',  &
                  'Lateral distance from yaw axis to rotor shaft (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! ShftSkew - Rotor shaft skew angle (read in degrees and converted to radians here).

   CALL ReadVar ( UnIn, FurlFile, InputFileData%ShftSkew, 'ShftSkew', 'Rotor shaft skew angle (deg)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%ShftSkew  = InputFileData%ShftSkew *D2R


      ! RFrlCMxn - Downwind distance from tower-top to CM of structure that furls with the rotor (not including rotor).

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlCMxn, 'RFrlCMxn',  &
                  'Downwind distance from tower-top to rotor-furl CM (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlCMyn - Lateral  distance from tower-top to CM of structure that furls with the rotor (not including rotor).

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlCMyn, 'RFrlCMyn',  &
                  'Lateral  distance from tower-top to rotor-furl CM (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlCMzn - Vertical distance from tower-top to CM of structure that furls with the rotor (not including rotor).

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlCMzn, 'RFrlCMzn',  &
                  'Vertical distance from tower-top to rotor-furl CM (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! BoomCMxn - Downwind distance from tower-top to tail boom CM.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%BoomCMxn, 'BoomCMxn',  &
                  'Downwind distance from tower-top to tail boom CM (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! BoomCMyn - Lateral  distance from tower-top to tail boom CM.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%BoomCMyn, 'BoomCMyn',  &
                  'Lateral  distance from tower-top to tail boom CM (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! BoomCMzn - Vertical distance from tower-top to tail boom CM.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%BoomCMzn, 'BoomCMzn', &
                   'Vertical distance from tower-top to tail boom CM (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFinCMxn - Downwind distance from tower-top to tail fin CM.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinCMxn, 'TFinCMxn', &
                   'Downwind distance from tower-top to tail fin CM (m)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFinCMyn - Lateral  distance from tower-top to tail fin CM.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinCMyn, 'TFinCMyn', &
                   'Lateral  distance from tower-top to tail fin CM (m)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFinCMzn - Vertical distance from tower-top to tail fin CM.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinCMzn, 'TFinCMzn', &
                   'Vertical distance from tower-top to tail fin CM (m)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFinCPxn - Downwind distance from tower-top to tail fin CP.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinCPxn, 'TFinCPxn', &
                  'Downwind distance from tower-top to tail fin CP (m)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFinCPyn - Lateral  distance from tower-top to tail fin CP.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinCPyn, 'TFinCPyn', &
                  'Lateral  distance from tower-top to tail fin CP (m)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFinCPzn - Vertical distance from tower-top to tail fin CP.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinCPzn, 'TFinCPzn', &
                  'Vertical distance from tower-top to tail fin CP (m)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFinSkew - Tail fin chordline skew angle (read in degrees, converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinSkew, 'TFinSkew', 'Tail fin chordline skew angle (deg)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%TFinSkew  = InputFileData%TFinSkew*D2R


      ! TFinTilt - Tail fin chordline tilt angle (read in degrees, converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinTilt, 'TFinTilt', 'Tail fin chordline tilt angle (deg)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%TFinTilt  = InputFileData%TFinTilt *D2R


      ! TFinBank - Tail fin planform  bank angle (read in degrees, converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinBank, 'TFinBank', 'Tail fin planform  bank angle (deg)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%TFinBank  = InputFileData%TFinBank *D2R


      ! RFrlPntxn - Downwind distance from tower-top to arbitrary point on rotor-furl axis.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlPntxn, 'RFrlPntxn', &
                  'Downwind distance from tower-top to arbitrary point on rotor-furl axis (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlPntyn - Lateral  distance from tower-top to arbitrary point on rotor-furl axis.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlPntyn, 'RFrlPntyn', &
                  'Lateral  distance from tower-top to arbitrary point on rotor-furl axis (m)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlPntzn - Vertical distance from tower-top to arbitrary point on rotor-furl axis.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlPntzn, 'RFrlPntzn', &
                  'Vertical distance from tower-top to arbitrary point on rotor-furl axis (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlSkew - Rotor-furl axis skew angle (read in degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlSkew, 'RFrlSkew', 'Rotor-furl axis skew angle (deg)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%RFrlSkew  = InputFileData%RFrlSkew*D2R


      ! RFrlTilt - Rotor-furl axis tilt angle (read in degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlTilt, 'RFrlTilt', 'Rotor-furl axis tilt angle (deg)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%RFrlTilt  = InputFileData%RFrlTilt*D2R


      ! TFrlPntxn - Downwind distance from tower-top to arbitrary point on tail-furl axis.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlPntxn, 'TFrlPntxn', &
                  'Downwind distance from tower-top to arbitrary point on tail-furl axis (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlPntyn - Lateral  distance from tower-top to arbitrary point on tail-furl axis.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlPntyn, 'TFrlPntyn', &
                  'Lateral  distance from tower-top to arbitrary point on tail-furl axis (m)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlPntzn - Vertical distance from tower-top to arbitrary point on tail-furl axis.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlPntzn, 'TFrlPntzn', &
                  'Vertical distance from tower-top to arbitrary point on tail-furl axis (m)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlSkew - Tail-furl axis skew angle (read in degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlSkew, 'TFrlSkew', 'Tail-furl axis skew angle (deg)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%TFrlSkew  = InputFileData%TFrlSkew *D2R


      ! TFrlTilt - Tail-furl axis tilt angle (read in degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlTilt, 'TFrlTilt', 'Tail-furl axis tilt angle (deg)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%TFrlTilt  = InputFileData%TFrlTilt *D2R


   !  -------------- MASS AND INERTIA (CONT) --------------------------------------

      ! Skip the comment line.

   CALL ReadCom ( UnIn, FurlFile, 'mass and inertia (cont)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! RFrlMass - Mass of structure that furls with the rotor (not including rotor).

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlMass, 'RFrlMass', 'Rotor-furl mass (kg)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! BoomMass - Tail boom mass.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%BoomMass, 'BoomMass', 'Tail boom mass (kg)',ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFinMass - Tail fin mass.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFinMass, 'TFinMass', 'Tail fin mass (kg)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlIner - Inertia of structure that furls with the rotor about the rotor-furl axis (not including rotor).

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlIner, 'RFrlIner', 'Rotor-furl inertia about rotor-furl axis (kg m^2)', &
                                            ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlIner - Tail boom inertia about tail-furl axis.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlIner, 'TFrlIner', 'Tail boom inertia about tail-furl axis (kg m^2)', &
                  ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !  -------------- ROTOR-FURL ---------------------------------------------------

      ! Skip the comment line.

   CALL ReadCom ( UnIn, FurlFile, 'Section heading: Rotor-Furl', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlMod - Rotor-furl spring/damper model switch.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlMod, 'RFrlMod', 'Rotor-furl spring/damper model switch', &
                  ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlSpr - Rotor-furl spring constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlSpr, 'RFrlSpr', 'Rotor-furl spring constant (N-m/rad)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlDmp - Rotor-furl damping constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlDmp, 'RFrlDmp', 'Rotor-furl damping constant (N-m/(rad/s))', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlCDmp - Rotor-furl rate-independent Coulomb-damping moment.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlCDmp, 'RFrlCDmp', 'Rotor-furl rate-independent Coulomb-damping moment (N-m)', &
                                                         ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlUSSP - Rotor-furl up-stop spring position (read in degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlUSSP, 'RFrlUSSP', 'Rotor-furl up-stop spring position (deg)', &
                  ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%RFrlUSSP  = InputFileData%RFrlUSSP*D2R


      ! RFrlDSSP - Rotor-furl down-stop spring position (read in degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlDSSP, 'RFrlDSSP', 'Rotor-furl down-stop spring position (deg)', &
                  ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%RFrlDSSP  = InputFileData%RFrlDSSP*D2R


      ! RFrlUSSpr - Rotor-furl up-stop spring constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlUSSpr, 'RFrlUSSpr', 'Rotor-furl up-stop spring constant (N-m/rad)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlDSSpr - Rotor-furl down-stop spring constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlDSSpr, 'RFrlDSSpr', 'Rotor-furl down-stop spring constant (N-m/rad)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlUSDP - Rotor-furl up-stop damper position (read in degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlUSDP, 'RFrlUSDP', 'Rotor-furl up-stop damper position (deg)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%RFrlUSDP  = InputFileData%RFrlUSDP*D2R


      ! RFrlDSDP - Rotor-furl down-stop damper position (read in degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlDSDP, 'RFrlDSDP', 'Rotor-furl down-stop damper position (deg)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%RFrlDSDP  = InputFileData%RFrlDSDP*D2R


      ! RFrlUSDmp - Rotor-furl up-stop damping constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlUSDmp, 'RFrlUSDmp', 'Rotor-furl up-stop damping constant (N-m/(rad/s))', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! RFrlDSDmp - Rotor-furl down-stop damping constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%RFrlDSDmp, 'RFrlDSDmp', 'Rotor-furl down-stop damping constant (N-m/(rad/s))', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- TAIL-FURL ----------------------------------------------------

      ! Skip the comment line.

   CALL ReadCom ( UnIn, FurlFile, 'tail-furl', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlMod - Tail-furl spring/damper model switch.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlMod, 'TFrlMod', 'Tail-furl spring/damper model switch', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlSpr - Tail-furl spring constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlSpr, 'TFrlSpr', 'Tail-furl spring constant (N-m/rad)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlDmp - Tail-furl damping constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlDmp, 'TFrlDmp', 'Tail-furl damping constant (N-m/(rad/s))', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlCDmp - Tail-furl rate-independent Coulomb-damping moment.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlCDmp, 'TFrlCDmp', 'Tail-furl rate-independent Coulomb-damping moment (N-m)', &
                                                                                ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlUSSP - Tail-furl up-stop spring position (read as degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlUSSP, 'TFrlUSSP', 'Tail-furl up-stop spring position (deg)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%TFrlUSSP  = InputFileData%TFrlUSSP*D2R


      ! TFrlDSSP - Tail-furl down-stop spring position (read as degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlDSSP, 'TFrlDSSP', 'Tail-furl down-stop spring position (deg)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%TFrlDSSP  = InputFileData%TFrlDSSP*D2R


      ! TFrlUSSpr - Tail-furl up-stop spring constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlUSSpr, 'TFrlUSSpr', 'Tail-furl up-stop spring constant (N-m/rad)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlDSSpr - Tail-furl down-stop spring constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlDSSpr, 'TFrlDSSpr', 'Tail-furl down-stop spring constant (N-m/rad)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlUSDP - Tail-furl up-stop damper position.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlUSDP, 'TFrlUSDP', 'Tail-furl up-stop damper position (deg)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%TFrlUSDP  = InputFileData%TFrlUSDP*D2R


      ! TFrlDSDP - Tail-furl down-stop damper position (read as degrees and converted to radians here)

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlDSDP, 'TFrlDSDP', 'Tail-furl down-stop damper position (deg)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%TFrlDSDP  = InputFileData%TFrlDSDP*D2R


      ! TFrlUSDmp - Tail-furl up-stop damping constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlUSDmp, 'TFrlUSDmp', 'Tail-furl up-stop damping constant (N-m/(rad/s))', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TFrlDSDmp - Tail-furl down-stop damping constant.

   CALL ReadVar ( UnIn, FurlFile, InputFileData%TFrlDSDmp, 'TFrlDSDmp', 'Tail-furl down-stop damping constant (N-m/(rad/s))', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Close the ElastoDyn furling file:

   CLOSE ( UnIn )

   RETURN
CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadFurlFile:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
         END IF

      END IF


   END SUBROUTINE CheckError
   !...............................................................................................................................
END SUBROUTINE ReadFurlFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadTowerFile( TwrFile, InputFileData, ReadAdmVals, UnEc, ErrStat, ErrMsg )
! This routine reads the tower file  input.
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables:

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                ! I/O unit for echo file. If present and > 0, write to UnEc
   LOGICAL,                  INTENT(IN)     :: ReadAdmVals                         ! Logical to determine if Adams inputs should be read from file
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message
   CHARACTER(*),             INTENT(IN)     :: TwrFile                             ! Name of the tower input file data
   TYPE(ED_InputFile),       INTENT(INOUT)  :: InputFileData                       ! All the data in the ElastoDyn input file


      ! Local variables:

   REAL(ReKi)                   :: AdjFASt                                         ! Factor to adjust tower fore-aft stiffness
   REAL(ReKi)                   :: AdjSSSt                                         ! Factor to adjust tower side-to-side stiffness
   REAL(ReKi)                   :: AdjTwMa                                         ! Factor to adjust tower mass density

   REAL(ReKi)                   :: TmpRAry(10)                                     ! Temporary variable to read table from file (up to 10 columns)

   INTEGER(IntKi)               :: I                                               ! A generic DO index.
   INTEGER(IntKi)               :: UnIn                                            ! Unit number for reading file
   INTEGER(IntKi)               :: NInputCols                                      ! Number of columns to be read from the file
   INTEGER(IntKi)               :: ErrStat2                                        ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                         ! Temporary Err msg



   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Open the tower input file.

   CALL OpenFInpFile ( UnIn, TwrFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Add a separator to the echo file if appropriate.
   IF ( UnEc > 0 )  WRITE (UnEc,'(//,A,/)')  'Tower input data from file "'//TRIM( TwrFile )//'":'


   !  -------------- FILE HEADER ---------------------------------------------------

   CALL ReadCom ( UnIn, TwrFile, 'unused tower file header line 1', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn, TwrFile, 'unused tower file header line 2', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- TOWER PARAMETERS ---------------------------------------------

   CALL ReadCom ( UnIn, TwrFile, 'heading for tower parameters', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! NTwInpSt - Number of tower input stations.

   CALL ReadVar ( UnIn, TwrFile, InputFileData%NTwInpSt, 'NTwInpSt', 'Number of tower input stations', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Allocate the input arrays based on this NTwInpSt input
   CALL Alloc_TowerInputProperties( InputFileData, ReadAdmVals, ErrStat, ErrMsg )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TwrFADmp - Tower fore-aft structural damping ratios.

   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwrFADmp, SIZE(InputFileData%TwrFADmp), 'TwrFADmp', &
                                     'Tower fore-aft structural damping ratios (%)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TwrSSDmp - Tower side-to-side structural damping ratios.

   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwrSSDmp, SIZE(InputFileData%TwrSSDmp), 'TwrSSDmp', &
                                     'Tower side-to-side structural damping ratios (%)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- TOWER ADJUSTMENT FACTORS -------------------------------------


      ! Skip the comment line.
   CALL ReadCom ( UnIn, TwrFile, 'heading for tower adjustment factors', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! FAStTunr - Tower fore-aft modal stiffness tuners.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%FAStTunr, SIZE(InputFileData%FAStTunr), 'FAStTunr', &
                                     'Tower fore-aft modal stiffness tuners (-)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! SSStTunr - Tower side-to-side modal stiffness tuners.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%SSStTunr, SIZE(InputFileData%SSStTunr), 'SSStTunr', &
                                     'Tower side-to-side modal stiffness tuners (-)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! AdjTwMa - Factor to adjust tower mass density.

   CALL ReadVar ( UnIn, TwrFile, AdjTwMa, 'AdjTwMa', 'Factor to adjust tower mass density (-)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! AdjFASt - Factor to adjust tower fore-aft stiffness.

   CALL ReadVar ( UnIn, TwrFile, AdjFASt, 'AdjFASt', 'Factor to adjust tower fore-aft stiffness (-)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! AdjSSSt - Factor to adjust tower side-to-side stiffness.

   CALL ReadVar ( UnIn, TwrFile, AdjSSSt, 'AdjSSSt', 'Factor to adjust tower side-to-side stiffness (-)', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


         ! Check the locally-defined adjustment factors: AdjTwMa, AdjFASt, AdjSSSt

   IF ( AdjTwMa <= 0.0_ReKi ) THEN
      CALL CheckError( ErrID_Warn, ' AdjTwMa must be greater than zero.' )
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF ( AdjFASt <= 0.0_ReKi ) THEN
      CALL CheckError( ErrID_Warn, ' AdjFASt must be greater than zero.' )
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF

   IF ( AdjSSSt <= 0.0_ReKi ) THEN
      CALL CheckError( ErrID_Warn, ' AdjSSSt must be greater than zero.' )
      IF ( ErrStat >= AbortErrLev ) RETURN
   END IF


   !  -------------- DISTRIBUTED TOWER PROPERTIES ---------------------------------

      ! Skip the comment lines.
   CALL ReadCom ( UnIn, TwrFile, 'heading for distributed tower parameters', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn, TwrFile, 'distributed-tower-parameter names', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn, TwrFile, 'distributed-tower-parameter units', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! Read the table.

   IF ( ReadAdmVals ) THEN
      NInputCols = 10
   ELSE
      NInputCols = 4
   END IF


   DO I=1,InputFileData%NTwInpSt

      CALL ReadAry( UnIn, TwrFile, TmpRAry, NInputCols, 'Line'//TRIM(Num2LStr(I)), 'Tower input station table', &
                    ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

      InputFileData%HtFract( I) = TmpRAry(1)
      InputFileData%TMassDen(I) = TmpRAry(2)*AdjTwMa   ! Apply the correction factors to the elemental data.
      InputFileData%TwFAStif(I) = TmpRAry(3)*AdjFASt   ! Apply the correction factors to the elemental data.
      InputFileData%TwSSStif(I) = TmpRAry(4)*AdjSSSt   ! Apply the correction factors to the elemental data.

      IF ( NInputCols > 4 ) THEN
         InputFileData%TwGJStif(I) = TmpRAry( 5)
         InputFileData%TwEAStif(I) = TmpRAry( 6)
         InputFileData%TwFAIner(I) = TmpRAry( 7)
         InputFileData%TwSSIner(I) = TmpRAry( 8)
         InputFileData%TwFAcgOf(I) = TmpRAry( 9)
         InputFileData%TwSScgOf(I) = TmpRAry(10)
      END IF

   END DO ! I


   !  -------------- TOWER FORE-AFT MODE SHAPES -----------------------------------


      ! Skip the comment line.
   CALL ReadCom ( UnIn, TwrFile, 'heading for tower fore-aft mode shapes', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TwFAM1Sh - Tower fore-aft mode-1 shape coefficients.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwFAM1Sh, SIZE(InputFileData%TwFAM1Sh), 'TwFAM1Sh', &
                           'Tower fore-aft mode-1 shape coefficients (-)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! TwFAM2Sh - Tower fore-aft mode-2 shape coefficients.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwFAM2Sh, SIZE(InputFileData%TwFAM2Sh), 'TwFAM2Sh', &
                           'Tower fore-aft mode-2 shape coefficients  (-)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- TOWER SIDE-TO-SIDE MODE SHAPES -------------------------------


      ! Skip the comment line.
   CALL ReadCom ( UnIn, TwrFile, 'heading for tower side-to-side mode shapes', ErrStat2, ErrMsg2, UnEc  )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! TwSSM1Sh - Tower side-to-side mode-1 shape coefficients.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwSSM1Sh, SIZE(InputFileData%TwSSM1Sh), 'TwSSM1Sh', &
                           'Tower side-to-side mode-1 shape coefficients (-)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN



      ! TwSSM2Sh - Tower side-to-side mode-2 shape coefficients.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwSSM2Sh, SIZE(InputFileData%TwSSM2Sh), 'TwSSM2Sh', &
                           'Tower side-to-side mode-2 shape coefficients (-)', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Close the tower file.
   CLOSE ( UnIn )


   RETURN
CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadTowerFile:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
         END IF

      END IF


   END SUBROUTINE CheckError
   !...............................................................................................................................
END SUBROUTINE ReadTowerFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadPrimaryFile( InputFile, InputFileData, BldFile, FurlFile, TwrFile, OutFileRoot, UnEc, ErrStat, ErrMsg )
! This routine reads in the primary ElastoDyn input file and places the values it reads in the InputFileData structure.
!   It opens an echo file if requested and returns the (still-open) echo file to the calling routine.
!   It also returns the names of the BldFile, FurlFile, and TrwFile for further reading of inputs.
!..................................................................................................................................


   IMPLICIT                        NONE

      ! Passed variables
   INTEGER(IntKi),     INTENT(OUT)    :: UnEc                                ! I/O unit for echo file. If > 0, file is open for writing.
   INTEGER(IntKi),     INTENT(OUT)    :: ErrStat                             ! Error status

   CHARACTER(*),       INTENT(IN)     :: InputFile                           ! Name of the file containing the primary input data
   CHARACTER(*),       INTENT(OUT)    :: ErrMsg                              ! Error message
   CHARACTER(*),       INTENT(OUT)    :: TwrFile                             ! name of the file containing tower inputs
   CHARACTER(*),       INTENT(OUT)    :: FurlFile                            ! name of the file containing furling inputs
   CHARACTER(*),       INTENT(OUT)    :: BldFile(MaxBl)                      ! name of the files containing blade inputs
   CHARACTER(*),       INTENT(IN)     :: OutFileRoot                         ! The rootname of the echo file, possibly opened in this routine

   TYPE(ED_InputFile), INTENT(INOUT)  :: InputFileData                       ! All the data in the ElastoDyn input file

      ! Local variables:
   INTEGER(IntKi)               :: I                                         ! loop counter
!   INTEGER(IntKi)               :: NumOuts                                  ! Number of output channel names read from the file
   INTEGER(IntKi)               :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)               :: IOS
   INTEGER(IntKi)               :: ErrStat2                                  ! Temporary Error status
   LOGICAL                      :: Echo                                      ! Determines if an echo file should be written
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(1024)              :: PriPath                                   ! Path name of the primary file
   CHARACTER(1024)              :: FTitle                                    ! "File Title": the 2nd line of the input file, which contains a description of its contents
   CHARACTER(200)               :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
   
      ! Initialize some variables:
   Echo = .FALSE.
   UnEc = -1                             ! Echo file not opened, yet
   CALL GetPath( InputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.


      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! Allocate arrays for input, based on maximum allowed number of blades and outputs
   CALL AllocAry( InputFileData%BlPitch, MaxBl, 'BlPitch input array', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL AllocAry( InputFileData%PreCone, MaxBl, 'Precone input array', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL AllocAry( InputFileData%TipMass, MaxBl, 'TipMass input array', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL AllocAry( InputFileData%OutList, MaxOutPts, "ElastoDyn Input File's Outlist", ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN


   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file.
   ! If Echo is TRUE, rewind and write on the second try.

   I    = 1 ! the number of times we've read the file (used for the Echo variable)
   DO
   !-------------------------- HEADER ---------------------------------------------
      CALL ReadCom( UnIn, InputFile, 'File Header: Module Version (line 1)', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL ReadStr( UnIn, InputFile, FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- SIMULATION CONTROL --------------------------------------
      CALL ReadCom( UnIn, InputFile, 'Section Header: Simulation Control', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

         ! Echo - Echo input to "<RootName>.ech".

      CALL ReadVar( UnIn, InputFile, Echo, 'Echo',   'Echo switch', ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN


      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop

         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read

      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)

      CALL OpenEcho ( UnEc, TRIM(OutFileRoot)//'.ech', ErrStat2, ErrMsg2, ED_Ver )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN

      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(ED_Ver%Name)//' primary input file "'//TRIM( InputFile )//'":'

      REWIND( UnIn, IOSTAT=ErrStat2 )
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL CheckError( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".' )
            IF ( ErrStat >= AbortErrLev ) RETURN
         END IF

   END DO

   IF (NWTC_VerboseLevel == NWTC_Verbose) THEN
      CALL WrScr( ' Heading of the '//TRIM(ED_Ver%Name)//' input file: ' )
      CALL WrScr( '   '//TRIM( FTitle ) )
   END IF

      ! Method - Integration method for loose coupling
   CALL ReadVar( UnIn, InputFile, InputFileData%method, "Method", "Requested integration method for ElastoDyn {1: RK4, 2: AB4, or 3: ABM4}", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      ! DT - Requested integration time for ElastoDyn (seconds):
   CALL ReadVar( UnIn, InputFile, Line, "DT", "Requested integration time for ElastoDyn (seconds)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "DEFAULT" ) /= 1 ) THEN ! If it's not "default", read this variable; otherwise use the value already stored in InputFileData%DT
         READ( Line, *, IOSTAT=IOS) InputFileData%DT
         IF ( IOS /= 0 ) THEN
            CALL CheckIOS ( IOS, InputFile, "DT", NumType, ErrStat2, ErrMsg2 )
            CALL CheckError( ErrStat2, ErrMsg2 )
            RETURN
         END IF
      END IF

   !---------------------- ENVIRONMENTAL CONDITION ---------------------------------
      CALL ReadCom( UnIn, InputFile, 'Section Header: Environmental Condition', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Gravity - Gravitational acceleration (m/s^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%Gravity, "Gravity", "Gravitational acceleration (m/s^2)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- DEGREES OF FREEDOM --------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Feature Flags', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! FlapDOF1 - First flapwise blade mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%FlapDOF1, "FlapDOF1", "First flapwise blade mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! FlapDOF2 - Second flapwise blade mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%FlapDOF2, "FlapDOF2", "Second flapwise blade mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! EdgeDOF - Edgewise blade mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%EdgeDOF, "EdgeDOF", "Edgewise blade mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TeetDOF - Rotor-teeter DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetDOF, "TeetDOF", "Rotor-teeter DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! DrTrDOF - Drivetrain rotational-flexibility DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%DrTrDOF, "DrTrDOF", "Drivetrain rotational-flexibility DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! GenDOF - Generator DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%GenDOF, "GenDOF", "Generator DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! YawDOF - Nacelle-yaw DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%YawDOF, "YawDOF", "Nacelle-yaw DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TwFADOF1 - First tower fore-aft bending-mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwFADOF1, "TwFADOF1", "First tower fore-aft bending-mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TwFADOF2 - Second tower fore-aft bending-mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwFADOF2, "TwFADOF2", "Second tower fore-aft bending-mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TwSSDOF1 - First tower side-to-side bending-mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwSSDOF1, "TwSSDOF1", "First tower side-to-side bending-mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TwSSDOF2 - Second tower side-to-side bending-mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwSSDOF2, "TwSSDOF2", "Second tower side-to-side bending-mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmSgDOF - Platform horizontal surge translation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmSgDOF, "PtfmSgDOF", "Platform horizontal surge translation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmSwDOF - Platform horizontal sway translation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmSwDOF, "PtfmSwDOF", "Platform horizontal sway translation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmHvDOF - Platform vertical heave translation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmHvDOF, "PtfmHvDOF", "Platform vertical heave translation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmRDOF - Platform roll tilt rotation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmRDOF, "PtfmRDOF", "Platform roll tilt rotation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmPDOF - Platform pitch tilt rotation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmPDOF, "PtfmPDOF", "Platform pitch tilt rotation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmYDOF - Platform yaw rotation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmYDOF, "PtfmYDOF", "Platform yaw rotation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- INITIAL CONDITIONS --------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Initial Conditions', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! OoPDefl - Initial out-of-plane blade-tip displacement (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%OoPDefl, "OoPDefl", "Initial out-of-plane blade-tip displacement (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! IPDefl - Initial in-plane blade-tip deflection (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%IPDefl, "IPDefl", "Initial in-plane blade-tip deflection (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! BlPitch - Initial blade pitch angles (deg) (read from file in degrees and converted to radians here):
   CALL ReadAryLines( UnIn, InputFile, InputFileData%BlPitch, MaxBl, "BlPitch", "Initial blade pitch angles (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%BlPitch = InputFileData%BlPitch*D2R

      ! TeetDefl - Initial teeter angle (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetDefl, "TeetDefl", "Initial teeter angle (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%TeetDefl = InputFileData%TeetDefl*D2R

      ! Azimuth - Initial azimuth angle for blade 1 (degrees) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%Azimuth, "Azimuth", "Initial azimuth angle for blade 1 (degrees)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%Azimuth = InputFileData%Azimuth*D2R

      ! RotSpeed - Initial rotor speed (RPM) (read in RPM and converted to rad/sec here):
   CALL ReadVar( UnIn, InputFile, InputFileData%RotSpeed, "RotSpeed", "Initial rotor speed (RPM)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%RotSpeed = InputFileData%RotSpeed*RPM2RPS

      ! NacYaw - Initial nacelle-yaw angle (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacYaw, "RotSpeed", "Initial nacelle-yaw angle (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%NacYaw = InputFileData%NacYaw*D2R

      ! TTDspFA - Initial fore-aft tower-top displacement (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TTDspFA, "TTDspFA", "Initial fore-aft tower-top displacement (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TTDspSS - Initial side-to-side tower-top displacement (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TTDspSS, "TTDspSS", "Initial side-to-side tower-top displacement (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmSurge - Initial horizontal surge translational displacement of platform (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmSurge, "PtfmSurge", "Initial horizontal surge translational displacement of platform (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmSway - Initial horizontal sway translational displacement of platform (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmSway, "PtfmSway", "Initial horizontal sway translational displacement of platform (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmHeave - Initial vertical heave translational displacement of platform (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmHeave, "PtfmHeave", "Initial vertical heave translational displacement of platform (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmRoll - Initial roll tilt rotational displacement of platform (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmRoll, "PtfmRoll", "Initial roll tilt rotational displacement of platform (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%PtfmRoll = InputFileData%PtfmRoll*D2R

      ! PtfmPitch - Initial pitch tilt rotational displacement of platform (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmPitch, "PtfmPitch", "Initial pitch tilt rotational displacement of platform (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%PtfmPitch = InputFileData%PtfmPitch*D2R

      ! PtfmYaw - Initial yaw rotational displacement of platform (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmYaw, "PtfmYaw", "Initial yaw rotational displacement of platform (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%PtfmYaw = InputFileData%PtfmYaw*D2R

   !---------------------- TURBINE CONFIGURATION -----------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Turbine Configuration', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NumBl - Number of blades (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NumBl, "NumBl", "Number of blades (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TipRad - Preconed blade-tip radius (distance from the rotor apex to the blade tip) (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TipRad, "TipRad", "Preconed blade-tip radius (distance from the rotor apex to the blade tip) (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! HubRad - Preconed hub radius (distance from the rotor apex to the blade root) (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubRad, "HubRad", "Preconed hub radius (distance from the rotor apex to the blade root) (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PreCone - Rotor precone angles (deg) (read from file in degrees and converted to radians here):
   CALL ReadAryLines( UnIn, InputFile, InputFileData%PreCone, MaxBl, "PreCone", "Rotor precone angles (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%PreCone = InputFileData%PreCone*D2R

      ! HubCM - Distance from rotor apex to hub mass (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubCM, "HubCM", "Distance from rotor apex to hub mass (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! UndSling - Undersling length (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%UndSling, "UndSling", "Undersling length (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Delta3 - Delta-3 angle for teetering rotors (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%Delta3, "Delta3", "Delta-3 angle for teetering rotors (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%Delta3 = InputFileData%Delta3*D2R

      ! AzimB1Up - Azimuth value to use for I/O when blade 1 points up (degrees) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%AzimB1Up, "AzimB1Up", "Azimuth value to use for I/O when blade 1 points up (degrees)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%AzimB1Up = InputFileData%AzimB1Up*D2R

      ! OverHang - Distance from yaw axis to rotor apex or teeter pin (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%OverHang, "OverHang", "Distance from yaw axis to rotor apex or teeter pin (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! ShftGagL - Distance from hub or teeter pin to shaft strain gages (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%ShftGagL, "ShftGagL", "Distance from hub or teeter pin to shaft strain gages (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! ShftTilt - Rotor shaft tilt angle (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%ShftTilt, "ShftTilt", "Rotor shaft tilt angle (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%ShftTilt = InputFileData%ShftTilt*D2R

      ! NacCMxn - Downwind distance from tower-top to nacelle CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacCMxn, "NacCMxn", "Downwind distance from tower-top to nacelle CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NacCMyn - Lateral distance from tower-top to nacelle CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacCMyn, "NacCMyn", "Lateral distance from tower-top to nacelle CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NacCMzn - Vertical distance from tower-top to nacelle CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacCMzn, "NacCMzn", "Vertical distance from tower-top to nacelle CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NcIMUxn - Downwind distance from the tower-top to the nacelle IMU (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NcIMUxn, "NcIMUxn", "Downwind distance from the tower-top to the nacelle IMU (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NcIMUyn - Lateral distance from the tower-top to the nacelle IMU (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NcIMUyn, "NcIMUyn", "Lateral distance from the tower-top to the nacelle IMU (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NcIMUzn - Vertical distance from the tower-top to the nacelle IMU (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NcIMUzn, "NcIMUzn", "Vertical distance from the tower-top to the nacelle IMU (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Twr2Shft - Vertical distance from the tower-top to the rotor shaft (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%Twr2Shft, "Twr2Shft", "Vertical distance from the tower-top to the rotor shaft (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TowerHt - Height of tower above ground level [onshore] or MSL [offshore] (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TowerHt, "TowerHt", "Height of tower above ground level [onshore] or MSL [offshore] (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TowerBsHt - Height of tower base above ground level [onshore] or MSL [offshore] (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TowerBsHt, "TowerBsHt", "Height of tower base above ground level [onshore] or MSL [offshore] (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      ! PtfmCMxt - Downwind distance from the ground [onshore] or MSL [offshore] to the platform CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmCMxt, "PtfmCMxt", "Downwind distance from the ground [onshore] or MSL [offshore] to the platform CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      ! PtfmCMyt - Lateral distance from the ground [onshore] or MSL [offshore] to the platform CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmCMyt, "PtfmCMzt", "Lateral distance from the ground [onshore] or MSL [offshore] to the platform CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN      
      
      ! PtfmCMzt - Vertical distance from the ground [onshore] or MSL [offshore] to the platform CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmCMzt, "PtfmCMzt", "Vertical distance from the ground [onshore] or MSL [offshore] to the platform CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmRefzt - Vertical distance from the ground [onshore] or MSL [offshore] to the platform reference point (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmRefzt, "PtfmRefzt", "Vertical distance from the ground [onshore] or MSL [offshore] to the platform reference point (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- MASS AND INERTIA ----------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Mass and Inertia', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TipMass - Tip-brake masses (kg):
   CALL ReadAryLines( UnIn, InputFile, InputFileData%TipMass, MaxBl, "TipMass", "Tip-brake masses (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! HubMass - Hub mass (kg):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubMass, "HubMass", "Hub mass (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! HubIner - Hub inertia about teeter axis (2-blader) or rotor axis (3-blader) (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubIner, "HubIner", "Hub inertia about teeter axis (2-blader) or rotor axis (3-blader) (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! GenIner - Generator inertia about HSS (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%GenIner, "GenIner", "Generator inertia about HSS (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NacMass - Nacelle mass (kg):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacMass, "NacMass", "Nacelle mass (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NacYIner - Nacelle yaw inertia (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacYIner, "NacYIner", "Nacelle yaw inertia (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! YawBrMass - Yaw bearing mass (kg):
   CALL ReadVar( UnIn, InputFile, InputFileData%YawBrMass, "YawBrMass", "Yaw bearing mass (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmMass - Platform mass (kg):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmMass, "PtfmMass", "Platform mass (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmRIner - Platform inertia for roll tilt rotation about the platform CM (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmRIner, "PtfmRIner", "Platform inertia for roll tilt rotation about the platform CM (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmPIner - Platform inertia for pitch tilt rotation about the platform CM (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmPIner, "PtfmPIner", "Platform inertia for pitch tilt rotation about the platform CM (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! PtfmYIner - Platform inertia for yaw rotation about the platform CM (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmYIner, "PtfmYIner", "Platform inertia for yaw rotation about the platform CM (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- BLADE ---------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Blade', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! BldNodes - Number of blade nodes (per blade) used for analysis (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%BldNodes, "BldNodes", "Number of blade nodes (per blade) used for analysis (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      
      ! InpBl - Blade file Input data for individual blades (see BladeInputData type):
   DO I = 1,MaxBl
      CALL ReadVar ( UnIn, InputFile, BldFile(I), 'BldFile('//TRIM(Num2Lstr(I))//')', 'Name of the file containing properties for blade '//TRIM(Num2Lstr(I)), ErrStat2, ErrMsg2, UnEc )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF ( ErrStat >= AbortErrLev ) RETURN
      IF ( PathIsRelative( BldFile(I) ) ) BldFile(I) = TRIM(PriPath)//TRIM(BldFile(I))
   END DO

   !---------------------- ROTOR-TEETER --------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Rotor-Teeter', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TeetMod - Rotor-teeter spring/damper model switch (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetMod, "TeetMod", "Rotor-teeter spring/damper model switch (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TeetDmpP - Rotor-teeter damper position (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetDmpP, "TeetDmpP", "Rotor-teeter damper position (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%TeetDmpP = InputFileData%TeetDmpP*D2R

      ! TeetDmp - Rotor-teeter damping constant (N-m/(rad/s)):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetDmp, "TeetDmp", "Rotor-teeter damping constant (N-m/(rad/s))", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TeetCDmp - Rotor-teeter rate-independent Coulomb-damping (N-m):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetCDmp, "TeetCDmp", "Rotor-teeter rate-independent Coulomb-damping (N-m)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TeetSStP - Rotor-teeter soft-stop position (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetSStP, "TeetSStP", "Rotor-teeter soft-stop position (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%TeetSStP = InputFileData%TeetSStP*D2R

      ! TeetHStP - Rotor-teeter hard-stop position (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetHStP, "TeetHStP", "Rotor-teeter hard-stop position (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%TeetHStP = InputFileData%TeetHStP*D2R

      ! TeetSSSp - Rotor-teeter soft-stop linear-spring constant (N-m/rad):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetSSSp, "TeetSSSp", "Rotor-teeter soft-stop linear-spring constant (N-m/rad)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TeetHSSp - Rotor-teeter hard-stop linear-spring constant (N-m/rad):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetHSSp, "TeetHSSp", "Rotor-teeter hard-stop linear-spring constant (N-m/rad)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- DRIVETRAIN ----------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Drivetrain', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! GBoxEff - Gearbox efficiency (%) (read from file in % and converted to fraction here):
   CALL ReadVar( UnIn, InputFile, InputFileData%GBoxEff, "GBoxEff", "Gearbox efficiency (%)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   InputFileData%GBoxEff = InputFileData%GBoxEff*0.01_ReKi

      ! GBRatio - Gearbox ratio (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%GBRatio, "GBRatio", "Gearbox ratio (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! DTTorSpr - Drivetrain torsional spring (N-m/rad):
   CALL ReadVar( UnIn, InputFile, InputFileData%DTTorSpr, "DTTorSpr", "Drivetrain torsional spring (N-m/rad)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! DTTorDmp - Drivetrain torsional damper (N-m/(rad/s)):
   CALL ReadVar( UnIn, InputFile, InputFileData%DTTorDmp, "DTTorDmp", "Drivetrain torsional damper (N-m/(rad/s))", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- FURLING -------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Furling', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Furling - Use Additional Furling parameters? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%Furling, "Furling", "Use Additional Furling parameters? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! FurlFile - Name of the file containing furling parameters:
   CALL ReadVar ( UnIn, InputFile, FurlFile, 'FurlFile', 'Name of the file containing furling parameters', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( FurlFile ) ) FurlFile = TRIM(PriPath)//TRIM(FurlFile)

   !---------------------- TOWER ---------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Tower', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TwrNodes - Number of tower nodes used in the analysis (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrNodes, "TwrNodes", "Number of tower nodes used in the analysis (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! TwrFile - Name of the file containing tower properties:
   CALL ReadVar ( UnIn, InputFile, TwrFile, 'TwrFile', 'Name of the file containing tower properties', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN
   IF ( PathIsRelative( TwrFile ) ) TwrFile = TRIM(PriPath)//TRIM(TwrFile)

   !---------------------- OUTPUT --------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Output', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! SumPrint - Print summary data to <RootName>.sum (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%SumPrint, "SumPrint", "Print summary data to <RootName>.sum (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! OutFile - Switch to determine where output will be placed: (1: in module output file only; 2: in glue code output file only; 3: both) (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%OutFile, "OutFile", "Switch to determine where output will be placed: (1: in module output file only; 2: in glue code output file only; 3: both) (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !    OutFileFmt - Format for module tabular (time-marching) output: (1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both):
   !CALL ReadVar( UnIn, InputFile, InputFileData%OutFileFmt, "OutFileFmt", "Format for module tabular (time-marching) output: (1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both)", ErrStat2, ErrMsg2, UnEc)
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN

      ! TabDelim - Flag to cause tab-delimited text output (delimited by space otherwise) (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TabDelim, "TabDelim", "Flag to cause tab-delimited text output (delimited by space otherwise) (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! OutFmt - Format used for module's text tabular output (except time); resulting field should be 10 characters (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%OutFmt, "OutFmt", "Format used for module's text tabular output (except time); resulting field should be 10 characters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! Tstart - Time to start module's tabular output (seconds):
   CALL ReadVar( UnIn, InputFile, InputFileData%Tstart, "Tstart", "Time to start module's tabular output (seconds)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! DecFact - Decimation factor for module's tabular output (1=output every step) (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%DecFact, "DecFact", "Decimation factor for module's tabular output (1=output every step) (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NTwGages - Number of tower strain gages (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NTwGages, "NTwGages", "Number of tower strain gages (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      IF ( InputFileData%NTwGages > SIZE(InputFileData%TwrGagNd) ) THEN
         CALL CheckError( ErrID_Warn, ' Warning: number of tower strain gages exceeds '// &
                                      TRIM(Num2LStr(SIZE(InputFileData%TwrGagNd)))//'.')
         InputFileData%NTwGages = SIZE(InputFileData%TwrGagNd)
      END IF

      ! TwrGagNd - Nodes closest to the tower strain gages (-):
   CALL ReadAry( UnIn, InputFile, InputFileData%TwrGagNd, InputFileData%NTwGages, "TwrGagNd", "Nodes closest to the tower strain gages (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NBlGages - Number of blade strain gages (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NBlGages, "NBlGages", "Number of blade strain gages (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      IF ( InputFileData%NBlGages > SIZE(InputFileData%BldGagNd) ) THEN
         CALL CheckError( ErrID_Warn, ' Warning: number of blade strain gages exceeds '//&
                                        TRIM(Num2LStr(SIZE(InputFileData%BldGagNd))) //'.')
         InputFileData%NBlGages = SIZE(InputFileData%BldGagNd)
      END IF

      ! BldGagNd - Nodes closest to the blade strain gages (-):
   CALL ReadAry( UnIn, InputFile, InputFileData%BldGagNd, InputFileData%NBlGages, "BldGagNd", "Nodes closest to the blade strain gages (-)", ErrStat2, ErrMsg2, UnEc)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- OUTLIST  --------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: OutList', ErrStat2, ErrMsg2, UnEc )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! OutList - List of user-requested output channels (-):
   CALL ReadOutputList ( UnIn, InputFile, InputFileData%OutList, InputFileData%NumOuts, 'OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc  )     ! Routine in NWTC Subroutine Library
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   !---------------------- END OF FILE -----------------------------------------

   CLOSE ( UnIn )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ReadPrimaryFile:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close file, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CLOSE( UnIn )
         END IF

      END IF


   END SUBROUTINE CheckError
   !...............................................................................................................................
END SUBROUTINE ReadPrimaryFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ValidatePrimaryData( InputFileData, BD4Blades, ErrStat, ErrMsg )
! This routine validates the inputs from the primary input file.
! note that all angles are assumed to be in radians in this routine:
!..................................................................................................................................

      ! Passed variables:

   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                       ! All the data in the ElastoDyn input file
   LOGICAL,                  INTENT(IN)     :: BD4Blades                           ! Use BeamDyn for blades, thus ignore ElastoDyn blade info 
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message

      ! Local variables:
   REAL(ReKi)                               :: SmallAngleLimit_Rad                 ! Largest input angle considered "small" (check in input file), radians
   INTEGER(IntKi)                           :: I                                   ! loop counter
   INTEGER(IntKi)                           :: K                                   ! blade number
   INTEGER(IntKi)                           :: FmtWidth                            ! width of the field returned by the specified OutFmt
   INTEGER(IntKi)                           :: ErrStat2                            ! Temporary error status
   CHARACTER(ErrMsgLen)                     :: ErrMsg2                             ! Temporary rror message


      ! Initialize error status and angle limit defined locally (in correct units)

   ErrStat = ErrID_None
   ErrMsg  = ''
   SmallAngleLimit_Rad = SmallAngleLimit_Deg*D2R

      ! Make sure the number of blades is valid:
   IF ( ( InputFileData%NumBl < 2 ) .OR. ( InputFileData%NumBl > 3 ) ) THEN
      CALL SetErrors( ErrID_Fatal, 'NumBl must be either 2 or 3.')
   END IF

      ! Make sure the specified integration method makes sense:
   IF ( InputFileData%method .ne. Method_RK4) THEN
      IF ( InputFileData%method .ne. Method_AB4) THEN
         IF ( InputFileData%method .ne. Method_ABM4) THEN
            CALL SetErrors( ErrID_Fatal, 'Integration method must be 1 (RK4), 2 (AB4), or 3 (ABM4)' )
         END IF
      END IF
   END IF
   
   
      ! make sure GBoxEff is 100% for now
   IF ( .NOT. EqualRealNos( InputFileData%GBoxEff, 1.0_ReKi ) .and. InputFileData%method == method_rk4 ) CALL SetErrors( ErrID_Fatal, 'GBoxEff must be 1 (i.e., 100%) when using RK4.')
      
   
      ! Don't allow these parameters to be negative (i.e., they must be in the range (0,inf)):
   CALL ErrIfNegative( InputFileData%Gravity,   'Gravity',   ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%RotSpeed,  'RotSpeed',  ErrStat, ErrMsg )
   IF (.NOT. BD4Blades) CALL ErrIfNegative( InputFileData%TipRad,    'TipRad',    ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%HubRad,    'HubRad',    ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%DTTorSpr,  'DTTorSpr',  ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%DTTorDmp,  'DTTorDmp',  ErrStat, ErrMsg )

   CALL ErrIfNegative( InputFileData%PtfmMass,  'PtfmMass',  ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%PtfmRIner, 'PtfmRIner', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%PtfmPIner, 'PtfmPIner', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%PtfmYIner, 'PtfmYIner', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%YawBrMass, 'YawBrMass', ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%NacMass,   'NacMass',   ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%HubMass,   'HubMass',   ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%Twr2Shft,  'Twr2Shft',  ErrStat, ErrMsg )

   DO K=1,InputFileData%NumBl
      CALL ErrIfNegative( InputFileData%TipMass(K), 'TipMass('//TRIM( Num2LStr( K ) )//')',   ErrStat, ErrMsg )
   ENDDO ! K

   CALL ErrIfNegative( InputFileData%NacYIner,  'NacYIner',  ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%GenIner,   'GenIner',   ErrStat, ErrMsg )
   CALL ErrIfNegative( InputFileData%HubIner,   'HubIner',   ErrStat, ErrMsg )

      ! Check that TowerHt is in the range [0,inf):
   IF ( InputFileData%TowerHt <= 0.0_ReKi )     CALL SetErrors( ErrID_Fatal, 'TowerHt must be greater than zero.' )

      ! Check that these integers are in appropriate ranges:
   IF ( InputFileData%TwrNodes < 1_IntKi ) CALL SetErrors( ErrID_Fatal, 'TwrNodes must not be less than 1.' )

   IF ( InputFileData%BldNodes < 1_IntKi ) CALL SetErrors( ErrID_Fatal, 'BldNodes must not be less than 1.' )
   
      ! Check that the gearbox efficiency is valid:
   IF ( ( InputFileData%GBoxEff <= 0.0_ReKi ) .OR. ( InputFileData%GBoxEff > 100.0 ) ) THEN
      CALL SetErrors( ErrID_Fatal, 'GBoxEff must be in the range (0,1] (i.e., (0,100] percent).' )
   END IF

      ! warn if 2nd modes are enabled without their corresponding 1st modes

   IF ( InputFileData%FlapDOF2 .AND. ( .NOT. InputFileData%FlapDOF1 ) )  THEN  ! Print out warning when flap mode 1 is not enabled and flap mode 2 is enabled
      CALL SetErrors( ErrID_Warn, '2nd blade flap mode is enabled without the 1st. '//&
                    ' This designation is recommended only for debugging purposes.')
   ENDIF

   IF ( InputFileData%TwFADOF2 .AND. ( .NOT. InputFileData%TwFADOF1 ) )  THEN  ! Print out warning when tower fore-aft mode 1 is not enabled and fore-aft mode 2 is enabled
      CALL SetErrors( ErrID_Warn, '2nd tower fore-aft mode is enabled without the 1st. '//&
                    ' This designation is recommended only for debugging purposes.')
   ENDIF

   IF ( InputFileData%TwSSDOF2 .AND. ( .NOT. InputFileData%TwSSDOF1 ) )  THEN  ! Print out warning when tower side-to-side mode 1 is not enabled and side-to-side mode 2 is enabled
      CALL SetErrors( ErrID_Warn, '2nd tower side-to-side mode is enabled without the 1st. '//&
                    ' This designation is recommended only for debugging purposes.')
   ENDIF


      ! Check that turbine configuration makes sense:

   IF ( InputFileData%Furling ) THEN
      IF ( InputFileData%OverHang > 0.0_ReKi )  THEN   ! Print out warning when downwind turbine is modeled with furling.
         CALL SetErrors( ErrID_Warn, 'Furling designation (Furling = True) specified for downwind rotor '// &
                    'configuration (OverHang > 0). Check for possible errors in the input file(s).')
      END IF
   ENDIF

   IF ( InputFileData%TowerBsHt >= InputFileData%TowerHt ) CALL SetErrors( ErrID_Fatal, 'TowerBsHt must be less than TowerHt.')

   IF ( InputFileData%PtfmCMzt  > InputFileData%TowerBsHt ) &
      CALL SetErrors( ErrID_Fatal, 'PtfmCMzt must not be greater than TowerBsHt.')
   
   IF ( InputFileData%PtfmRefzt  > InputFileData%TowerBsHt ) &
      CALL SetErrors( ErrID_Fatal, 'PtfmRefzt must not be greater than TowerBsHt.')
     
   IF (.NOT. BD4Blades ) THEN
      IF (InputFileData%HubRad >= InputFileData%TipRad ) &
      CALL SetErrors( ErrID_Fatal, 'HubRad must be less than TipRad.' )

      IF ( InputFileData%TowerHt + InputFileData%Twr2Shft + InputFileData%OverHang*SIN(InputFileData%ShftTilt) &
                                 <= InputFileData%TipRad )  THEN
         CALL SetErrors( ErrID_Fatal, 'TowerHt + Twr2Shft + OverHang*SIN(ShftTilt) must be greater than TipRad.' )
      END IF
   END IF


   IF ( InputFileData%NumBl == 2_IntKi )  THEN
      IF ( ( InputFileData%TeetDefl <= -pi ) .OR. ( InputFileData%TeetDefl > pi ) )  &
         CALL SetErrors( ErrID_Fatal, 'TeetDefl must be in the range (-pi, pi] radians (i.e., [-180,180] degrees).' )

      IF ( ABS( InputFileData%Delta3 ) >= PiBy2 )  &
         CALL SetErrors( ErrID_Fatal, 'Delta3 must be in the range (pi/2, pi/2) radians (i.e., (-90, 90) degrees).' )

      IF ( ( InputFileData%TeetSStP < 0.0_ReKi ) .OR. ( InputFileData%TeetSStP > pi ) )  &
         CALL SetErrors( ErrID_Fatal, 'TeetSStP must be in the range [0, pi] radians (i.e., [0,180] degrees).' )

      IF ( ( InputFileData%TeetDmpP < 0.0_ReKi ) .OR. ( InputFileData%TeetDmpP > pi ) )  &
         CALL SetErrors( ErrID_Fatal, 'TeetDmpP must be in the range [0, pi] radians (i.e., [0,180] degrees).' )

      IF ( ( InputFileData%TeetHStP < InputFileData%TeetSStP ) .OR. ( InputFileData%TeetHStP > pi ) )  &
         CALL SetErrors( ErrID_Fatal, 'TeetHStP must be in the range [TeetSStP, pi] radians (i.e., [TeetSStP, 180] degrees).' )

      IF ( ( InputFileData%TeetMod /= 0_IntKi ) .AND. ( InputFileData%TeetMod /= 1_IntKi ) .AND. &
           ( InputFileData%TeetMod /= 2_IntKi ) )  &
         CALL SetErrors( ErrID_Fatal, 'TeetMod must be 0, 1, or 2.' )

      CALL ErrIfNegative( InputFileData%TeetDmp,   'TeetDmp',   ErrStat, ErrMsg )
      CALL ErrIfNegative( InputFileData%TeetCDmp,  'TeetCDmp',  ErrStat, ErrMsg )
      CALL ErrIfNegative( InputFileData%TeetSSSp,  'TeetSSSp',  ErrStat, ErrMsg )
      CALL ErrIfNegative( InputFileData%TeetHSSp,  'TeetHSSp',  ErrStat, ErrMsg )
      
      ! TeetCDmp isn't allowed to be non-zero in this verison.
      IF (.NOT. EqualRealNos( InputFileData%TeetCDmp, 0.0_ReKi )) &
         CALL SetErrors( ErrID_Fatal, 'TeetCDmp must be 0 in this version of ElastoDyn.'  )
      
   ENDIF

      ! check these angles for appropriate ranges:
   IF ( ( InputFileData%NacYaw <= -pi ) .OR. ( InputFileData%NacYaw > pi ) ) &
      CALL SetErrors( ErrID_Fatal, 'NacYaw must be in the range (-pi, pi] radians (i.e., (-180, 180] degrees).' )
   IF ( ( InputFileData%Azimuth  < 0.0_ReKi ) .OR. ( InputFileData%Azimuth  >= TwoPi ) ) &
      CALL SetErrors( ErrID_Fatal, 'Azimuth must be in the range [0, 2pi) radians (i.e., [0, 360) degrees).' )
   IF ( ( InputFileData%AzimB1Up < 0.0_ReKi ) .OR. ( InputFileData%AzimB1Up >= TwoPi ) )  &
      CALL SetErrors( ErrID_Fatal, 'AzimB1Up must be in the range [0, 2pi) radians (i.e., [0, 360) degrees).' )

   DO K=1,InputFileData%NumBl
      IF ( ( InputFileData%BlPitch(K) <= -pi ) .OR. ( InputFileData%BlPitch(K) > pi ) )  THEN
         CALL SetErrors( ErrID_Fatal, 'BlPitch('//TRIM(Num2LStr(K))//')'//' must be greater than -pi radians and '// &
                                      'less than or equal to pi radians (i.e., in the range (-180, 180] degrees).' )
      END IF
      IF ( ABS( InputFileData%PreCone(K) ) >= PiBy2 )  &
         CALL SetErrors( ErrID_Fatal, 'PreCone('//TRIM( Num2LStr( K ) )//') must be in the range (-pi/2, pi/2) '//&
                                      'radians (i.e., (-90, 90) degrees).' )
   END DO

      ! Check that these angles are in the range [-pi/2, pi/2] radians (i.e., [-90, 90] degrees ):
   CALL CheckAngle90Range( InputFileData%ShftTilt, 'ShftTilt', ErrStat, ErrMsg )


      ! Check for violations of the small-angle assumption (15-degree limit, using radians):
   IF ( ABS( InputFileData%PtfmRoll ) > SmallAngleLimit_Rad ) THEN
      CALL SetErrors( ErrID_Fatal, 'PtfmRoll must be between -'//TRIM(Num2LStr(SmallAngleLimit_Rad))//' and ' &
                                                               //TRIM(Num2LStr(SmallAngleLimit_Rad))//' radians.' )
   END IF

   IF ( ABS( InputFileData%PtfmPitch ) > SmallAngleLimit_Rad ) THEN
      CALL SetErrors( ErrID_Fatal, 'PtfmPitch must be between -'//TRIM(Num2LStr(SmallAngleLimit_Rad))//' and ' &
                                                                //TRIM(Num2LStr(SmallAngleLimit_Rad))//' radians.' )
   END IF

   IF ( ABS( InputFileData%PtfmYaw ) > SmallAngleLimit_Rad ) THEN
      CALL SetErrors( ErrID_Fatal, 'PtfmYaw must be between -'//TRIM(Num2LStr(SmallAngleLimit_Rad))//' and ' &
                                                              //TRIM(Num2LStr(SmallAngleLimit_Rad))//' radians.' )
   END IF

      ! Check the output parameters:
   IF ( InputFileData%DecFact < 1_IntKi )  CALL SetErrors( ErrID_Fatal, 'DecFact must be greater than 0.' )

   IF ( ( InputFileData%NTwGages < 0_IntKi ) .OR. ( InputFileData%NTwGages > 9_IntKi ) )  THEN
      CALL SetErrors( ErrID_Fatal, 'NTwGages must be between 0 and 9 (inclusive).' )
   ELSE
         ! Check to see if all TwrGagNd(:) analysis points are existing analysis points:

      DO I=1,InputFileData%NTwGages
         IF ( InputFileData%TwrGagNd(I) < 1_IntKi .OR. InputFileData%TwrGagNd(I) > InputFileData%TwrNodes ) THEN
            CALL SetErrors( ErrID_Fatal, ' All TwrGagNd values must be between 1 and '//&
                           TRIM( Num2LStr( InputFileData%TwrNodes ) )//' (inclusive).' )
            EXIT ! stop checking this loop
         END IF
      END DO         
   
   END IF
         
         
   IF ( ( InputFileData%NBlGages < 0_IntKi ) .OR. ( InputFileData%NBlGages > 9_IntKi ) )  THEN
      CALL SetErrors( ErrID_Fatal, 'NBlGages must be between 0 and 9 (inclusive).' )
   ELSE 

   ! Check to see if all BldGagNd(:) analysis points are existing analysis points:

      DO I=1,InputFileData%NBlGages
         IF ( InputFileData%BldGagNd(I) < 1_IntKi .OR. InputFileData%BldGagNd(I) > InputFileData%InpBlMesh(1)%BldNodes ) THEN
            CALL SetErrors( ErrID_Fatal, ' All BldGagNd values must be between 1 and '//&
                              TRIM( Num2LStr( InputFileData%InpBlMesh(1)%BldNodes ) )//' (inclusive).' )
            EXIT ! stop checking this loop
         END IF
      END DO
      
   END IF

      ! Check that InputFileData%OutFmt is a valid format specifier and will fit over the column headings
   CALL ChkRealFmtStr( InputFileData%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
   IF ( ErrStat2 /= ErrID_None ) CALL SetErrors(ErrStat2, ErrMsg2 )
   IF ( FmtWidth /= ChanLen ) CALL SetErrors(ErrID_Warn, 'OutFmt produces a column width of '//TRIM(Num2LStr(FmtWidth))//&
                                                           ' instead of '//TRIM(Num2LStr(ChanLen))//' characters.' )


   RETURN

CONTAINS
   !-------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE SetErrors( ErrStat3, ErrMsg3 )
   ! This routine sets the error message and flag when an error has occurred
   !...............................................................................................................................
   INTEGER(IntKi), INTENT(IN) :: ErrStat3     ! Error status for this error
   CHARACTER(*),   INTENT(IN) :: ErrMsg3      ! Error message for this error

      ErrStat = MAX( ErrStat, ErrStat3 )
      IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg  = TRIM(ErrMsg)//TRIM(ErrMsg3)

   END SUBROUTINE SetErrors
   !-------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE ValidatePrimaryData
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ErrIfNegative( Var, VarDesc, ErrStat, ErrMsg )
! This routine checks that an value is in the range [0, inf). If not, ErrStat = ErrID_Fatal.
! Note that ErrStat and ErrMsg are INTENT(INOUT).
!..................................................................................................................................
REAL(ReKi),     INTENT(IN)    :: Var         ! Variable to check
CHARACTER(*),   INTENT(IN)    :: VarDesc     ! Description of variable (used in error message)
INTEGER(IntKi), INTENT(INOUT) :: ErrStat     ! Error status to update if Var is not in specified range
CHARACTER(*),   INTENT(INOUT) :: ErrMsg      ! Error message to update if Var is not in specified range

   IF (  Var < 0.0_ReKi )  THEN
      ErrStat = ErrID_Fatal
      IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg  = TRIM(ErrMsg)//TRIM(VarDesc)//' must not be negative.'
   END IF

END SUBROUTINE ErrIfNegative
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CheckAngle90Range( Var, VarDesc, ErrStat, ErrMsg )
! This routine checks that an angle is in the range [-pi/2, pi/2] radians. If not, ErrStat = ErrID_Fatal
! Note that all values are assumed to be in radians, even if read in degrees ( [-90 deg, 90 deg] )
! Note that ErrStat and ErrMsg are INTENT(INOUT).
!...............................................................................................................................
REAL(ReKi),     INTENT(IN)    :: Var         ! Variable to check
CHARACTER(*),   INTENT(IN)    :: VarDesc     ! Description of variable (used in error message)
INTEGER(IntKi), INTENT(INOUT) :: ErrStat     ! Error status to update if Var is not in specified range
CHARACTER(*),   INTENT(INOUT) :: ErrMsg      ! Error message to update if Var is not in specified range

   IF ( ABS( Var ) > PiBy2 )  THEN
      ErrStat = ErrID_Fatal
      IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
      ErrMsg  = TRIM(ErrMsg)// &
                  TRIM(VarDesc)//' must be between -pi/2 and pi/2 radians (i.e., in the range [-90, 90] degrees).'
   END IF

END SUBROUTINE CheckAngle90Range
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_ContStates( x, p, InputFileData, OtherState, ErrStat, ErrMsg  )
! This routine initializes the continuous states of the module.
! It assumes the parameters are set and that InputFileData contains initial conditions for the continuous states.
!..................................................................................................................................
   TYPE(ED_ContinuousStateType), INTENT(OUT)    :: x                 ! Initial continuous states
   TYPE(ED_ParameterType),       INTENT(IN)     :: p                 ! Parameters of the structural dynamics module
   TYPE(ED_InputFile),           INTENT(IN)     :: InputFileData     ! Data stored in the module's input file
   TYPE(ED_OtherStateType),      INTENT(IN)     :: OtherState        ! Initial other states
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat           ! Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg            ! Error message

      ! local variables
   REAL(ReKi)                                   :: InitQE1(p%NumBl)  ! Initial value of the 1st blade edge DOF
   REAL(ReKi)                                   :: InitQF1(p%NumBl)  ! Initial value of the 1st blade flap DOF
   REAL(ReKi)                                   :: InitQF2(p%NumBl)  ! Initial value of the 2nd blade flap DOF
!   INTEGER(IntKi)                               :: I                 ! loop counter

      
      ! First allocate the arrays stored here:

   CALL AllocAry( x%QT, p%NDOF,   'QT',   ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   CALL AllocAry( x%QDT, p%NDOF,  'QDT',  ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   
   
      ! Calculate/apply the initial blade DOF values to the corresponding DOFs.
   IF (.NOT. p%BD4Blades) THEN  !Skipping subroutine if BeamDyn = TRUE
      CALL InitBlDefl ( p, InputFileData, InitQF1, InitQF2, InitQE1, ErrStat, ErrMsg  )
   ELSE
      InitQF1 = 0.0_ReKi
      InitQF2 = 0.0_ReKi
      InitQE1 = 0.0_ReKi
   END IF
   
      
   x%QT ( DOF_BF(1:p%NumBl,1) ) = InitQF1   ! These come from InitBlDefl().
   x%QT ( DOF_BF(1:p%NumBl,2) ) = InitQF2   ! These come from InitBlDefl().
   x%QT ( DOF_BE(1:p%NumBl,1) ) = InitQE1   ! These come from InitBlDefl().
   x%QDT( DOF_BF(1:p%NumBl,1) ) = 0.0
   x%QDT( DOF_BF(1:p%NumBl,2) ) = 0.0
   x%QDT( DOF_BE(1:p%NumBl,1) ) = 0.0

      ! Teeter Motion

   IF ( p%NumBl == 2 )  THEN !note, DOF_Teet doesn't exist for 3-bladed turbine, so don't include an ELSE here

      ! Set initial teeter angle to TeetDefl and initial teeter angular velocity to 0.

      x%QT (DOF_Teet) = InputFileData%TeetDefl
      x%QDT(DOF_Teet) = 0.0
   ENDIF

      ! Generator azimuth

      ! Set initial generator azimuth angle.  Turn rotor on, whether it is
      !   fixed or variable speed.  If it is fixed speed, set up the
      !   fixed rpm.

   !JASON: CHANGE THESE MOD() FUNCTIONS INTO MODULO() FUNCTIONS SO THAT YOU CAN ELIMINATE ADDING 360:
!   x%QT (DOF_GeAz) = MOD( (InputFileData%Azimuth - p%AzimB1Up)*R2D + 270.0 + 360.0, 360.0 )*D2R   ! Internal position of blade 1
   
   x%QT (DOF_GeAz) = InputFileData%Azimuth - p%AzimB1Up - Piby2
   CALL Zero2TwoPi( x%QT (DOF_GeAz) )
   x%QDT(DOF_GeAz) = p%RotSpeed                                               ! Rotor speed in rad/sec.


      ! Shaft compliance

   ! The initial shaft compliance displacements and velocities are all zero.
   !   They will remain zero if the drivetrain DOF is disabled:

   x%QT (DOF_DrTr) = 0.0
   x%QDT(DOF_DrTr) = 0.0




      ! Rotor-furl motion

      ! Set initial rotor-furl angle to RotFurl.  If rotor-furl is off, this
      !   becomes a fixed rotor-furl angle.

   x%QT (DOF_RFrl) = InputFileData%RotFurl
   x%QDT(DOF_RFrl) = 0.0



      ! Tail-furl motion

      ! Set initial tail-furl angle to TailFurl.  If tail-furl is off, this becomes a fixed tail-furl angle.

   x%QT (DOF_TFrl) = InputFileData%TailFurl
   x%QDT(DOF_TFrl) = 0.0



      ! Yaw Motion

      ! Set initial yaw angle to NacYaw.  If yaw is off, this becomes a fixed yaw angle.

   x%QT (DOF_Yaw) = InputFileData%NacYaw
   x%QDT(DOF_Yaw) = 0.0



      ! Tower motion

      ! Assign all the displacements to mode 1 unless it is disabled.  If mode 1
      !   is disabled and mode 2 is enabled, assign all displacements to mode 2.
      ! If both modes are disabled, set the displacements to zero.

   x%QT   (DOF_TFA1) =  0.0
   x%QT   (DOF_TSS1) =  0.0
   x%QT   (DOF_TFA2) =  0.0
   x%QT   (DOF_TSS2) =  0.0

   IF (    InputFileData%TwFADOF1 )  THEN   ! First fore-aft tower mode is enabled.
      x%QT(DOF_TFA1) =  InputFileData%TTDspFA
   ELSEIF( InputFileData%TwFADOF2 )  THEN   ! Second fore-aft tower mode is enabled, but first is not.
      x%QT(DOF_TFA2) =  InputFileData%TTDspFA
   ENDIF

   IF (    InputFileData%TwSSDOF1 )  THEN   ! First side-to-side tower mode is enabled.
      x%QT(DOF_TSS1) = -InputFileData%TTDspSS
   ELSEIF( InputFileData%TwSSDOF2 )  THEN   ! Second side-to-side tower mode is enabled, but first is not.
      x%QT(DOF_TSS2) = -InputFileData%TTDspSS
   ENDIF

   x%QDT  (DOF_TFA1) =  0.0
   x%QDT  (DOF_TSS1) =  0.0
   x%QDT  (DOF_TFA2) =  0.0
   x%QDT  (DOF_TSS2) =  0.0



      ! Platform Motion

      ! Set initial platform displacements.  If platform DOFs are off, these
      !   become fixed platform displacements.

   x%QT (DOF_Sg) = InputFileData%PtfmSurge
   x%QT (DOF_Sw) = InputFileData%PtfmSway
   x%QT (DOF_Hv) = InputFileData%PtfmHeave
   x%QT (DOF_R ) = InputFileData%PtfmRoll
   x%QT (DOF_P ) = InputFileData%PtfmPitch
   x%QT (DOF_Y ) = InputFileData%PtfmYaw
   x%QDT(DOF_Sg) = 0.0
   x%QDT(DOF_Sw) = 0.0
   x%QDT(DOF_Hv) = 0.0
   x%QDT(DOF_R ) = 0.0
   x%QDT(DOF_P ) = 0.0
   x%QDT(DOF_Y ) = 0.0

   


END SUBROUTINE Init_ContStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_OtherStates( OtherState, p, x, InputFileData, ErrStat, ErrMsg  )
! This routine initializes the other states of the module.
! It assumes the parameters are set and that InputFileData contains initial conditions for the continuous states. and p%NumBl is set
!..................................................................................................................................
   TYPE(ED_OtherStateType),      INTENT(OUT)    :: OtherState        ! Initial other states
   TYPE(ED_ParameterType),       INTENT(IN)     :: p                 ! Parameters of the structural dynamics module
   TYPE(ED_ContinuousStateType), INTENT(IN)     :: x                 ! Initial continuous states
   TYPE(ED_InputFile),           INTENT(IN)     :: InputFileData     ! Data stored in the module's input file
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat           ! Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg            ! Error message

   INTEGER(IntKi)               :: I                                 ! Generic loop counter.

      ! First allocate the arrays stored here:

   CALL Alloc_RtHS( OtherState%RtHS, p, ErrStat, ErrMsg  )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL Alloc_CoordSys( OtherState%CoordSys, p, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   
   CALL AllocAry( OtherState%QD2T, p%NDOF,   'OtherState%QD2T',  ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

   
   ALLOCATE ( OtherState%AllOuts(0:MaxOutPts) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         ErrStat = ErrID_Fatal
         ErrMsg  = ' Error allocating memory for the AllOuts array.'
         RETURN
      ENDIF   
   OtherState%AllOuts = 0.0_ReKi
   
      ! for loose coupling:
   CALL AllocAry( OtherState%IC,  ED_NMX,   'IC',   ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   
   
   
   !CALL AllocAry(OtherState%BlPitch, p%NumBl, 'BlPitch', ErrStat, ErrMsg )
   !      IF ( ErrStat >= AbortErrLev ) RETURN
   !   OtherState%BlPitch = InputFileData%BlPitch(1:p%NumBl)

   CALL AllocAry( OtherState%AugMat,       p%NDOF,          p%NAug,          'AugMat',       ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL AllocAry( OtherState%OgnlGeAzRo,                    p%NAug,          'OgnlGeAzRo',   ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL AllocAry( OtherState%SolnVec,      p%DOFs%NActvDOF,                  'SolnVec',      ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL AllocAry( OtherState%AugMat_pivot, p%DOFs%NActvDOF,                  'AugMat_pivot', ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN
   CALL AllocAry( OtherState%AugMat_factor,p%DOFs%NActvDOF, p%DOFs%NActvDOF, 'AugMat_factor',ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) RETURN

   
      ! Now initialize the IC array = (/NMX, NMX-1, ... , 1 /)
      ! this keeps track of the position in the array of continuous states (stored in other states)

   OtherState%IC(1) = ED_NMX
   DO I = 2,ED_NMX
      OtherState%IC(I) = OtherState%IC(I-1) - 1
   ENDDO


   !   ! Initialize the accelerations to zero.
   !
   !OtherState%QD2 = 0.0_ReKi


   OtherState%n   = -1  ! we haven't updated OtherState%xdot, yet
   
   DO i = LBOUND(OtherState%xdot,1), UBOUND(OtherState%xdot,1)
      CALL ED_CopyContState( x, OtherState%xdot(i), MESH_NEWCOPY, ErrStat, ErrMsg)
         IF ( ErrStat >= AbortErrLev ) RETURN 
   ENDDO
   
      ! hacks for HSS brake function:
   
   OtherState%HSSBrTrq   = 0.0_ReKi
   OtherState%HSSBrTrqC  = 0.0_ReKi
   OtherState%SgnPrvLSTQ = 1
   OtherState%SgnLSTQ    = 1
   
   
END SUBROUTINE Init_OtherStates

!----------------------------------------------------------------------------------------------------------------------------------

!**********************************************************************************************************************************
! NOTE: The following lines of code were generated by a Matlab script called "Write_ChckOutLst.m"
!      using the parameters listed in the "OutListParameters.xlsx" Excel file. Any changes to these 
!      lines should be modified in the Matlab script and/or Excel worksheet as necessary. 
! This code was generated by Write_ChckOutLst.m at 02-Mar-2015 10:37:31.
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetOutParam(OutList, p, ErrStat, ErrMsg )
! This routine checks to see if any requested output channel names (stored in the OutList(:)) are invalid. It returns a 
! warning if any of the channels are not available outputs from the module.
!  It assigns the settings for OutParam(:) (i.e, the index, name, and units of the output channels, WriteOutput(:)).
!  the sign is set to 0 if the channel is invalid.
! It sets assumes the value p%NumOuts has been set before this routine has been called, and it sets the values of p%OutParam here.
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables

   CHARACTER(ChanLen),        INTENT(IN)     :: OutList(:)                        ! The list out user-requested outputs
   TYPE(ED_ParameterType),    INTENT(INOUT)  :: p                                 ! The module parameters
   INTEGER(IntKi),            INTENT(OUT)    :: ErrStat                           ! The error status code
   CHARACTER(*),              INTENT(OUT)    :: ErrMsg                            ! The error message, if an error occurred

      ! Local variables

   INTEGER                      :: I                                               ! Generic loop-counting index
   INTEGER                      :: J                                               ! Generic loop-counting index
   INTEGER                      :: INDX                                            ! Index for valid arrays
   INTEGER                      :: startIndx                                       ! Index for BeamDyn

   LOGICAL                      :: CheckOutListAgain                               ! Flag used to determine if output parameter starting with "M" is valid (or the negative of another parameter)
   LOGICAL                      :: InvalidOutput(0:MaxOutPts)                      ! This array determines if the output channel is valid for this configuration
   CHARACTER(ChanLen)           :: OutListTmp                                      ! A string to temporarily hold OutList(I)

   CHARACTER(OutStrLenM1), PARAMETER  :: ValidParamAry(972) =  (/ &                  ! This lists the names of the allowed parameters, which must be sorted alphabetically
                               "AZIMUTH  ","BLDPITCH1","BLDPITCH2","BLDPITCH3","BLPITCH1 ","BLPITCH2 ","BLPITCH3 ", &
                               "GENACCEL ","GENSPEED ","HSSBRTQ  ","HSSHFTA  ","HSSHFTPWR","HSSHFTTQ ","HSSHFTV  ", &
                               "IPDEFL1  ","IPDEFL2  ","IPDEFL3  ","LSSGAGA  ","LSSGAGAXA","LSSGAGAXS","LSSGAGFXA", &
                               "LSSGAGFXS","LSSGAGFYA","LSSGAGFYS","LSSGAGFZA","LSSGAGFZS","LSSGAGMXA","LSSGAGMXS", &
                               "LSSGAGMYA","LSSGAGMYS","LSSGAGMZA","LSSGAGMZS","LSSGAGP  ","LSSGAGPXA","LSSGAGPXS", &
                               "LSSGAGV  ","LSSGAGVXA","LSSGAGVXS","LSSHFTFXA","LSSHFTFXS","LSSHFTFYA","LSSHFTFYS", &
                               "LSSHFTFZA","LSSHFTFZS","LSSHFTMXA","LSSHFTMXS","LSSHFTPWR","LSSHFTTQ ","LSSTIPA  ", &
                               "LSSTIPAXA","LSSTIPAXS","LSSTIPMYA","LSSTIPMYS","LSSTIPMZA","LSSTIPMZS","LSSTIPP  ", &
                               "LSSTIPPXA","LSSTIPPXS","LSSTIPV  ","LSSTIPVXA","LSSTIPVXS","NACYAW   ","NACYAWA  ", &
                               "NACYAWP  ","NACYAWV  ","NCIMURAXS","NCIMURAYS","NCIMURAZS","NCIMURVXS","NCIMURVYS", &
                               "NCIMURVZS","NCIMUTAXS","NCIMUTAYS","NCIMUTAZS","NCIMUTVXS","NCIMUTVYS","NCIMUTVZS", &
                               "OOPDEFL1 ","OOPDEFL2 ","OOPDEFL3 ","PTCHDEFL1","PTCHDEFL2","PTCHDEFL3","PTCHPMZB1", &
                               "PTCHPMZB2","PTCHPMZB3","PTCHPMZC1","PTCHPMZC2","PTCHPMZC3","PTFMHEAVE","PTFMPITCH", &
                               "PTFMRAXI ","PTFMRAXT ","PTFMRAYI ","PTFMRAYT ","PTFMRAZI ","PTFMRAZT ","PTFMRDXI ", &
                               "PTFMRDYI ","PTFMRDZI ","PTFMROLL ","PTFMRVXI ","PTFMRVXT ","PTFMRVYI ","PTFMRVYT ", &
                               "PTFMRVZI ","PTFMRVZT ","PTFMSURGE","PTFMSWAY ","PTFMTAXI ","PTFMTAXT ","PTFMTAYI ", &
                               "PTFMTAYT ","PTFMTAZI ","PTFMTAZT ","PTFMTDXI ","PTFMTDXT ","PTFMTDYI ","PTFMTDYT ", &
                               "PTFMTDZI ","PTFMTDZT ","PTFMTVXI ","PTFMTVXT ","PTFMTVYI ","PTFMTVYT ","PTFMTVZI ", &
                               "PTFMTVZT ","PTFMYAW  ","QD2_B1E1 ","QD2_B1F1 ","QD2_B1F2 ","QD2_B2E1 ","QD2_B2F1 ", &
                               "QD2_B2F2 ","QD2_B3E1 ","QD2_B3F1 ","QD2_B3F2 ","QD2_DRTR ","QD2_GEAZ ","QD2_HV   ", &
                               "QD2_P    ","QD2_R    ","QD2_RFRL ","QD2_SG   ","QD2_SW   ","QD2_TEET ","QD2_TFA1 ", &
                               "QD2_TFA2 ","QD2_TFRL ","QD2_TSS1 ","QD2_TSS2 ","QD2_Y    ","QD2_YAW  ","QD_B1E1  ", &
                               "QD_B1F1  ","QD_B1F2  ","QD_B2E1  ","QD_B2F1  ","QD_B2F2  ","QD_B3E1  ","QD_B3F1  ", &
                               "QD_B3F2  ","QD_DRTR  ","QD_GEAZ  ","QD_HV    ","QD_P     ","QD_R     ","QD_RFRL  ", &
                               "QD_SG    ","QD_SW    ","QD_TEET  ","QD_TFA1  ","QD_TFA2  ","QD_TFRL  ","QD_TSS1  ", &
                               "QD_TSS2  ","QD_Y     ","QD_YAW   ","Q_B1E1   ","Q_B1F1   ","Q_B1F2   ","Q_B2E1   ", &
                               "Q_B2F1   ","Q_B2F2   ","Q_B3E1   ","Q_B3F1   ","Q_B3F2   ","Q_DRTR   ","Q_GEAZ   ", &
                               "Q_HV     ","Q_P      ","Q_R      ","Q_RFRL   ","Q_SG     ","Q_SW     ","Q_TEET   ", &
                               "Q_TFA1   ","Q_TFA2   ","Q_TFRL   ","Q_TSS1   ","Q_TSS2   ","Q_Y      ","Q_YAW    ", &
                               "RFRLBRM  ","ROLLDEFL1","ROLLDEFL2","ROLLDEFL3","ROOTFXB1 ","ROOTFXB2 ","ROOTFXB3 ", &
                               "ROOTFXC1 ","ROOTFXC2 ","ROOTFXC3 ","ROOTFYB1 ","ROOTFYB2 ","ROOTFYB3 ","ROOTFYC1 ", &
                               "ROOTFYC2 ","ROOTFYC3 ","ROOTFZB1 ","ROOTFZB2 ","ROOTFZB3 ","ROOTFZC1 ","ROOTFZC2 ", &
                               "ROOTFZC3 ","ROOTMEDG1","ROOTMEDG2","ROOTMEDG3","ROOTMFLP1","ROOTMFLP2","ROOTMFLP3", &
                               "ROOTMIP1 ","ROOTMIP2 ","ROOTMIP3 ","ROOTMOOP1","ROOTMOOP2","ROOTMOOP3","ROOTMXB1 ", &
                               "ROOTMXB2 ","ROOTMXB3 ","ROOTMXC1 ","ROOTMXC2 ","ROOTMXC3 ","ROOTMYB1 ","ROOTMYB2 ", &
                               "ROOTMYB3 ","ROOTMYC1 ","ROOTMYC2 ","ROOTMYC3 ","ROOTMZB1 ","ROOTMZB2 ","ROOTMZB3 ", &
                               "ROOTMZC1 ","ROOTMZC2 ","ROOTMZC3 ","ROTACCEL ","ROTFURL  ","ROTFURLA ","ROTFURLP ", &
                               "ROTFURLV ","ROTPWR   ","ROTSPEED ","ROTTEETA ","ROTTEETP ","ROTTEETV ","ROTTHRUST", &
                               "ROTTORQ  ","SPN1ALXB1","SPN1ALXB2","SPN1ALXB3","SPN1ALYB1","SPN1ALYB2","SPN1ALYB3", &
                               "SPN1ALZB1","SPN1ALZB2","SPN1ALZB3","SPN1FLXB1","SPN1FLXB2","SPN1FLXB3","SPN1FLYB1", &
                               "SPN1FLYB2","SPN1FLYB3","SPN1FLZB1","SPN1FLZB2","SPN1FLZB3","SPN1MLXB1","SPN1MLXB2", &
                               "SPN1MLXB3","SPN1MLYB1","SPN1MLYB2","SPN1MLYB3","SPN1MLZB1","SPN1MLZB2","SPN1MLZB3", &
                               "SPN1RDXB1","SPN1RDXB2","SPN1RDXB3","SPN1RDYB1","SPN1RDYB2","SPN1RDYB3","SPN1RDZB1", &
                               "SPN1RDZB2","SPN1RDZB3","SPN1TDXB1","SPN1TDXB2","SPN1TDXB3","SPN1TDYB1","SPN1TDYB2", &
                               "SPN1TDYB3","SPN1TDZB1","SPN1TDZB2","SPN1TDZB3","SPN2ALXB1","SPN2ALXB2","SPN2ALXB3", &
                               "SPN2ALYB1","SPN2ALYB2","SPN2ALYB3","SPN2ALZB1","SPN2ALZB2","SPN2ALZB3","SPN2FLXB1", &
                               "SPN2FLXB2","SPN2FLXB3","SPN2FLYB1","SPN2FLYB2","SPN2FLYB3","SPN2FLZB1","SPN2FLZB2", &
                               "SPN2FLZB3","SPN2MLXB1","SPN2MLXB2","SPN2MLXB3","SPN2MLYB1","SPN2MLYB2","SPN2MLYB3", &
                               "SPN2MLZB1","SPN2MLZB2","SPN2MLZB3","SPN2RDXB1","SPN2RDXB2","SPN2RDXB3","SPN2RDYB1", &
                               "SPN2RDYB2","SPN2RDYB3","SPN2RDZB1","SPN2RDZB2","SPN2RDZB3","SPN2TDXB1","SPN2TDXB2", &
                               "SPN2TDXB3","SPN2TDYB1","SPN2TDYB2","SPN2TDYB3","SPN2TDZB1","SPN2TDZB2","SPN2TDZB3", &
                               "SPN3ALXB1","SPN3ALXB2","SPN3ALXB3","SPN3ALYB1","SPN3ALYB2","SPN3ALYB3","SPN3ALZB1", &
                               "SPN3ALZB2","SPN3ALZB3","SPN3FLXB1","SPN3FLXB2","SPN3FLXB3","SPN3FLYB1","SPN3FLYB2", &
                               "SPN3FLYB3","SPN3FLZB1","SPN3FLZB2","SPN3FLZB3","SPN3MLXB1","SPN3MLXB2","SPN3MLXB3", &
                               "SPN3MLYB1","SPN3MLYB2","SPN3MLYB3","SPN3MLZB1","SPN3MLZB2","SPN3MLZB3","SPN3RDXB1", &
                               "SPN3RDXB2","SPN3RDXB3","SPN3RDYB1","SPN3RDYB2","SPN3RDYB3","SPN3RDZB1","SPN3RDZB2", &
                               "SPN3RDZB3","SPN3TDXB1","SPN3TDXB2","SPN3TDXB3","SPN3TDYB1","SPN3TDYB2","SPN3TDYB3", &
                               "SPN3TDZB1","SPN3TDZB2","SPN3TDZB3","SPN4ALXB1","SPN4ALXB2","SPN4ALXB3","SPN4ALYB1", &
                               "SPN4ALYB2","SPN4ALYB3","SPN4ALZB1","SPN4ALZB2","SPN4ALZB3","SPN4FLXB1","SPN4FLXB2", &
                               "SPN4FLXB3","SPN4FLYB1","SPN4FLYB2","SPN4FLYB3","SPN4FLZB1","SPN4FLZB2","SPN4FLZB3", &
                               "SPN4MLXB1","SPN4MLXB2","SPN4MLXB3","SPN4MLYB1","SPN4MLYB2","SPN4MLYB3","SPN4MLZB1", &
                               "SPN4MLZB2","SPN4MLZB3","SPN4RDXB1","SPN4RDXB2","SPN4RDXB3","SPN4RDYB1","SPN4RDYB2", &
                               "SPN4RDYB3","SPN4RDZB1","SPN4RDZB2","SPN4RDZB3","SPN4TDXB1","SPN4TDXB2","SPN4TDXB3", &
                               "SPN4TDYB1","SPN4TDYB2","SPN4TDYB3","SPN4TDZB1","SPN4TDZB2","SPN4TDZB3","SPN5ALXB1", &
                               "SPN5ALXB2","SPN5ALXB3","SPN5ALYB1","SPN5ALYB2","SPN5ALYB3","SPN5ALZB1","SPN5ALZB2", &
                               "SPN5ALZB3","SPN5FLXB1","SPN5FLXB2","SPN5FLXB3","SPN5FLYB1","SPN5FLYB2","SPN5FLYB3", &
                               "SPN5FLZB1","SPN5FLZB2","SPN5FLZB3","SPN5MLXB1","SPN5MLXB2","SPN5MLXB3","SPN5MLYB1", &
                               "SPN5MLYB2","SPN5MLYB3","SPN5MLZB1","SPN5MLZB2","SPN5MLZB3","SPN5RDXB1","SPN5RDXB2", &
                               "SPN5RDXB3","SPN5RDYB1","SPN5RDYB2","SPN5RDYB3","SPN5RDZB1","SPN5RDZB2","SPN5RDZB3", &
                               "SPN5TDXB1","SPN5TDXB2","SPN5TDXB3","SPN5TDYB1","SPN5TDYB2","SPN5TDYB3","SPN5TDZB1", &
                               "SPN5TDZB2","SPN5TDZB3","SPN6ALXB1","SPN6ALXB2","SPN6ALXB3","SPN6ALYB1","SPN6ALYB2", &
                               "SPN6ALYB3","SPN6ALZB1","SPN6ALZB2","SPN6ALZB3","SPN6FLXB1","SPN6FLXB2","SPN6FLXB3", &
                               "SPN6FLYB1","SPN6FLYB2","SPN6FLYB3","SPN6FLZB1","SPN6FLZB2","SPN6FLZB3","SPN6MLXB1", &
                               "SPN6MLXB2","SPN6MLXB3","SPN6MLYB1","SPN6MLYB2","SPN6MLYB3","SPN6MLZB1","SPN6MLZB2", &
                               "SPN6MLZB3","SPN6RDXB1","SPN6RDXB2","SPN6RDXB3","SPN6RDYB1","SPN6RDYB2","SPN6RDYB3", &
                               "SPN6RDZB1","SPN6RDZB2","SPN6RDZB3","SPN6TDXB1","SPN6TDXB2","SPN6TDXB3","SPN6TDYB1", &
                               "SPN6TDYB2","SPN6TDYB3","SPN6TDZB1","SPN6TDZB2","SPN6TDZB3","SPN7ALXB1","SPN7ALXB2", &
                               "SPN7ALXB3","SPN7ALYB1","SPN7ALYB2","SPN7ALYB3","SPN7ALZB1","SPN7ALZB2","SPN7ALZB3", &
                               "SPN7FLXB1","SPN7FLXB2","SPN7FLXB3","SPN7FLYB1","SPN7FLYB2","SPN7FLYB3","SPN7FLZB1", &
                               "SPN7FLZB2","SPN7FLZB3","SPN7MLXB1","SPN7MLXB2","SPN7MLXB3","SPN7MLYB1","SPN7MLYB2", &
                               "SPN7MLYB3","SPN7MLZB1","SPN7MLZB2","SPN7MLZB3","SPN7RDXB1","SPN7RDXB2","SPN7RDXB3", &
                               "SPN7RDYB1","SPN7RDYB2","SPN7RDYB3","SPN7RDZB1","SPN7RDZB2","SPN7RDZB3","SPN7TDXB1", &
                               "SPN7TDXB2","SPN7TDXB3","SPN7TDYB1","SPN7TDYB2","SPN7TDYB3","SPN7TDZB1","SPN7TDZB2", &
                               "SPN7TDZB3","SPN8ALXB1","SPN8ALXB2","SPN8ALXB3","SPN8ALYB1","SPN8ALYB2","SPN8ALYB3", &
                               "SPN8ALZB1","SPN8ALZB2","SPN8ALZB3","SPN8FLXB1","SPN8FLXB2","SPN8FLXB3","SPN8FLYB1", &
                               "SPN8FLYB2","SPN8FLYB3","SPN8FLZB1","SPN8FLZB2","SPN8FLZB3","SPN8MLXB1","SPN8MLXB2", &
                               "SPN8MLXB3","SPN8MLYB1","SPN8MLYB2","SPN8MLYB3","SPN8MLZB1","SPN8MLZB2","SPN8MLZB3", &
                               "SPN8RDXB1","SPN8RDXB2","SPN8RDXB3","SPN8RDYB1","SPN8RDYB2","SPN8RDYB3","SPN8RDZB1", &
                               "SPN8RDZB2","SPN8RDZB3","SPN8TDXB1","SPN8TDXB2","SPN8TDXB3","SPN8TDYB1","SPN8TDYB2", &
                               "SPN8TDYB3","SPN8TDZB1","SPN8TDZB2","SPN8TDZB3","SPN9ALXB1","SPN9ALXB2","SPN9ALXB3", &
                               "SPN9ALYB1","SPN9ALYB2","SPN9ALYB3","SPN9ALZB1","SPN9ALZB2","SPN9ALZB3","SPN9FLXB1", &
                               "SPN9FLXB2","SPN9FLXB3","SPN9FLYB1","SPN9FLYB2","SPN9FLYB3","SPN9FLZB1","SPN9FLZB2", &
                               "SPN9FLZB3","SPN9MLXB1","SPN9MLXB2","SPN9MLXB3","SPN9MLYB1","SPN9MLYB2","SPN9MLYB3", &
                               "SPN9MLZB1","SPN9MLZB2","SPN9MLZB3","SPN9RDXB1","SPN9RDXB2","SPN9RDXB3","SPN9RDYB1", &
                               "SPN9RDYB2","SPN9RDYB3","SPN9RDZB1","SPN9RDZB2","SPN9RDZB3","SPN9TDXB1","SPN9TDXB2", &
                               "SPN9TDXB3","SPN9TDYB1","SPN9TDYB2","SPN9TDYB3","SPN9TDZB1","SPN9TDZB2","SPN9TDZB3", &
                               "TAILFURL ","TAILFURLA","TAILFURLP","TAILFURLV","TEETAYA  ","TEETDEFL ","TEETPYA  ", &
                               "TEETVYA  ","TFRLBRM  ","TIP2TWR1 ","TIP2TWR2 ","TIP2TWR3 ","TIPALXB1 ","TIPALXB2 ", &
                               "TIPALXB3 ","TIPALYB1 ","TIPALYB2 ","TIPALYB3 ","TIPALZB1 ","TIPALZB2 ","TIPALZB3 ", &
                               "TIPCLRNC1","TIPCLRNC2","TIPCLRNC3","TIPDXB1  ","TIPDXB2  ","TIPDXB3  ","TIPDXC1  ", &
                               "TIPDXC2  ","TIPDXC3  ","TIPDYB1  ","TIPDYB2  ","TIPDYB3  ","TIPDYC1  ","TIPDYC2  ", &
                               "TIPDYC3  ","TIPDZB1  ","TIPDZB2  ","TIPDZB3  ","TIPDZC1  ","TIPDZC2  ","TIPDZC3  ", &
                               "TIPRDXB1 ","TIPRDXB2 ","TIPRDXB3 ","TIPRDYB1 ","TIPRDYB2 ","TIPRDYB3 ","TIPRDZB1 ", &
                               "TIPRDZB2 ","TIPRDZB3 ","TIPRDZC1 ","TIPRDZC2 ","TIPRDZC3 ","TTDSPAX  ","TTDSPFA  ", &
                               "TTDSPPTCH","TTDSPROLL","TTDSPSS  ","TTDSPTWST","TWHT1ALXT","TWHT1ALYT","TWHT1ALZT", &
                               "TWHT1FLXT","TWHT1FLYT","TWHT1FLZT","TWHT1MLXT","TWHT1MLYT","TWHT1MLZT","TWHT1RDXT", &
                               "TWHT1RDYT","TWHT1RDZT","TWHT1RPXI","TWHT1RPYI","TWHT1RPZI","TWHT1TDXT","TWHT1TDYT", &
                               "TWHT1TDZT","TWHT1TPXI","TWHT1TPYI","TWHT1TPZI","TWHT2ALXT","TWHT2ALYT","TWHT2ALZT", &
                               "TWHT2FLXT","TWHT2FLYT","TWHT2FLZT","TWHT2MLXT","TWHT2MLYT","TWHT2MLZT","TWHT2RDXT", &
                               "TWHT2RDYT","TWHT2RDZT","TWHT2RPXI","TWHT2RPYI","TWHT2RPZI","TWHT2TDXT","TWHT2TDYT", &
                               "TWHT2TDZT","TWHT2TPXI","TWHT2TPYI","TWHT2TPZI","TWHT3ALXT","TWHT3ALYT","TWHT3ALZT", &
                               "TWHT3FLXT","TWHT3FLYT","TWHT3FLZT","TWHT3MLXT","TWHT3MLYT","TWHT3MLZT","TWHT3RDXT", &
                               "TWHT3RDYT","TWHT3RDZT","TWHT3RPXI","TWHT3RPYI","TWHT3RPZI","TWHT3TDXT","TWHT3TDYT", &
                               "TWHT3TDZT","TWHT3TPXI","TWHT3TPYI","TWHT3TPZI","TWHT4ALXT","TWHT4ALYT","TWHT4ALZT", &
                               "TWHT4FLXT","TWHT4FLYT","TWHT4FLZT","TWHT4MLXT","TWHT4MLYT","TWHT4MLZT","TWHT4RDXT", &
                               "TWHT4RDYT","TWHT4RDZT","TWHT4RPXI","TWHT4RPYI","TWHT4RPZI","TWHT4TDXT","TWHT4TDYT", &
                               "TWHT4TDZT","TWHT4TPXI","TWHT4TPYI","TWHT4TPZI","TWHT5ALXT","TWHT5ALYT","TWHT5ALZT", &
                               "TWHT5FLXT","TWHT5FLYT","TWHT5FLZT","TWHT5MLXT","TWHT5MLYT","TWHT5MLZT","TWHT5RDXT", &
                               "TWHT5RDYT","TWHT5RDZT","TWHT5RPXI","TWHT5RPYI","TWHT5RPZI","TWHT5TDXT","TWHT5TDYT", &
                               "TWHT5TDZT","TWHT5TPXI","TWHT5TPYI","TWHT5TPZI","TWHT6ALXT","TWHT6ALYT","TWHT6ALZT", &
                               "TWHT6FLXT","TWHT6FLYT","TWHT6FLZT","TWHT6MLXT","TWHT6MLYT","TWHT6MLZT","TWHT6RDXT", &
                               "TWHT6RDYT","TWHT6RDZT","TWHT6RPXI","TWHT6RPYI","TWHT6RPZI","TWHT6TDXT","TWHT6TDYT", &
                               "TWHT6TDZT","TWHT6TPXI","TWHT6TPYI","TWHT6TPZI","TWHT7ALXT","TWHT7ALYT","TWHT7ALZT", &
                               "TWHT7FLXT","TWHT7FLYT","TWHT7FLZT","TWHT7MLXT","TWHT7MLYT","TWHT7MLZT","TWHT7RDXT", &
                               "TWHT7RDYT","TWHT7RDZT","TWHT7RPXI","TWHT7RPYI","TWHT7RPZI","TWHT7TDXT","TWHT7TDYT", &
                               "TWHT7TDZT","TWHT7TPXI","TWHT7TPYI","TWHT7TPZI","TWHT8ALXT","TWHT8ALYT","TWHT8ALZT", &
                               "TWHT8FLXT","TWHT8FLYT","TWHT8FLZT","TWHT8MLXT","TWHT8MLYT","TWHT8MLZT","TWHT8RDXT", &
                               "TWHT8RDYT","TWHT8RDZT","TWHT8RPXI","TWHT8RPYI","TWHT8RPZI","TWHT8TDXT","TWHT8TDYT", &
                               "TWHT8TDZT","TWHT8TPXI","TWHT8TPYI","TWHT8TPZI","TWHT9ALXT","TWHT9ALYT","TWHT9ALZT", &
                               "TWHT9FLXT","TWHT9FLYT","TWHT9FLZT","TWHT9MLXT","TWHT9MLYT","TWHT9MLZT","TWHT9RDXT", &
                               "TWHT9RDYT","TWHT9RDZT","TWHT9RPXI","TWHT9RPYI","TWHT9RPZI","TWHT9TDXT","TWHT9TDYT", &
                               "TWHT9TDZT","TWHT9TPXI","TWHT9TPYI","TWHT9TPZI","TWRBSFXT ","TWRBSFYT ","TWRBSFZT ", &
                               "TWRBSMXT ","TWRBSMYT ","TWRBSMZT ","TWRCLRNC1","TWRCLRNC2","TWRCLRNC3","TWSTDEFL1", &
                               "TWSTDEFL2","TWSTDEFL3","YAWACCEL ","YAWAZN   ","YAWAZP   ","YAWBRFXN ","YAWBRFXP ", &
                               "YAWBRFYN ","YAWBRFYP ","YAWBRFZN ","YAWBRFZP ","YAWBRMXN ","YAWBRMXP ","YAWBRMYN ", &
                               "YAWBRMYP ","YAWBRMZN ","YAWBRMZP ","YAWBRRAXP","YAWBRRAYP","YAWBRRAZP","YAWBRRDXT", &
                               "YAWBRRDYT","YAWBRRDZT","YAWBRRVXP","YAWBRRVYP","YAWBRRVZP","YAWBRTAXP","YAWBRTAYP", &
                               "YAWBRTAZP","YAWBRTDXP","YAWBRTDXT","YAWBRTDYP","YAWBRTDYT","YAWBRTDZP","YAWBRTDZT", &
                               "YAWPOS   ","YAWPZN   ","YAWPZP   ","YAWRATE  ","YAWVZN   ","YAWVZP   "/)
   INTEGER(IntKi), PARAMETER :: ParamIndxAry(972) =  (/ &                            ! This lists the index into AllOuts(:) of the allowed parameters ValidParamAry(:)
                                LSSTipPxa , PtchPMzc1 , PtchPMzc2 , PtchPMzc3 , PtchPMzc1 , PtchPMzc2 , PtchPMzc3 , &
                                  HSShftA ,   HSShftV ,   HSSBrTq ,   HSShftA , HSShftPwr ,  HSShftTq ,   HSShftV , &
                                  TipDyc1 ,   TipDyc2 ,   TipDyc3 , LSSGagAxa , LSSGagAxa , LSSGagAxa , LSShftFxa , &
                                LSShftFxa , LSShftFya , LSShftFys , LSShftFza , LSShftFzs , LSShftMxa , LSShftMxa , &
                                LSSGagMya , LSSGagMys , LSSGagMza , LSSGagMzs , LSSGagPxa , LSSGagPxa , LSSGagPxa , &
                                LSSGagVxa , LSSGagVxa , LSSGagVxa , LSShftFxa , LSShftFxa , LSShftFya , LSShftFys , &
                                LSShftFza , LSShftFzs , LSShftMxa , LSShftMxa ,    RotPwr , LSShftMxa , LSSTipAxa , &
                                LSSTipAxa , LSSTipAxa , LSSTipMya , LSSTipMys , LSSTipMza , LSSTipMzs , LSSTipPxa , &
                                LSSTipPxa , LSSTipPxa , LSSTipVxa , LSSTipVxa , LSSTipVxa ,    YawPzn ,    YawAzn , &
                                   YawPzn ,    YawVzn , NcIMURAxs , NcIMURAys , NcIMURAzs , NcIMURVxs , NcIMURVys , &
                                NcIMURVzs , NcIMUTAxs , NcIMUTAys , NcIMUTAzs , NcIMUTVxs , NcIMUTVys , NcIMUTVzs , &
                                  TipDxc1 ,   TipDxc2 ,   TipDxc3 ,  TipRDyb1 ,  TipRDyb2 ,  TipRDyb3 , PtchPMzc1 , &
                                PtchPMzc2 , PtchPMzc3 , PtchPMzc1 , PtchPMzc2 , PtchPMzc3 ,  PtfmTDzi ,  PtfmRDyi , &
                                 PtfmRAxi ,  PtfmRAxt ,  PtfmRAyi ,  PtfmRAyt ,  PtfmRAzi ,  PtfmRAzt ,  PtfmRDxi , &
                                 PtfmRDyi ,  PtfmRDzi ,  PtfmRDxi ,  PtfmRVxi ,  PtfmRVxt ,  PtfmRVyi ,  PtfmRVyt , &
                                 PtfmRVzi ,  PtfmRVzt ,  PtfmTDxi ,  PtfmTDyi ,  PtfmTAxi ,  PtfmTAxt ,  PtfmTAyi , &
                                 PtfmTAyt ,  PtfmTAzi ,  PtfmTAzt ,  PtfmTDxi ,  PtfmTDxt ,  PtfmTDyi ,  PtfmTDyt , &
                                 PtfmTDzi ,  PtfmTDzt ,  PtfmTVxi ,  PtfmTVxt ,  PtfmTVyi ,  PtfmTVyt ,  PtfmTVzi , &
                                 PtfmTVzt ,  PtfmRDzi ,  QD2_B1E1 ,  QD2_B1F1 ,  QD2_B1F2 ,  QD2_B2E1 ,  QD2_B2F1 , &
                                 QD2_B2F2 ,  QD2_B3E1 ,  QD2_B3F1 ,  QD2_B3F2 ,  QD2_DrTr ,  QD2_GeAz ,    QD2_Hv , &
                                    QD2_P ,     QD2_R ,  QD2_RFrl ,    QD2_Sg ,    QD2_Sw ,  QD2_Teet ,  QD2_TFA1 , &
                                 QD2_TFA2 ,  QD2_TFrl ,  QD2_TSS1 ,  QD2_TSS2 ,     QD2_Y ,   QD2_Yaw ,   QD_B1E1 , &
                                  QD_B1F1 ,   QD_B1F2 ,   QD_B2E1 ,   QD_B2F1 ,   QD_B2F2 ,   QD_B3E1 ,   QD_B3F1 , &
                                  QD_B3F2 ,   QD_DrTr ,   QD_GeAz ,     QD_Hv ,      QD_P ,      QD_R ,   QD_RFrl , &
                                    QD_Sg ,     QD_Sw ,   QD_Teet ,   QD_TFA1 ,   QD_TFA2 ,   QD_TFrl ,   QD_TSS1 , &
                                  QD_TSS2 ,      QD_Y ,    QD_Yaw ,    Q_B1E1 ,    Q_B1F1 ,    Q_B1F2 ,    Q_B2E1 , &
                                   Q_B2F1 ,    Q_B2F2 ,    Q_B3E1 ,    Q_B3F1 ,    Q_B3F2 ,    Q_DrTr ,    Q_GeAz , &
                                     Q_Hv ,       Q_P ,       Q_R ,    Q_RFrl ,      Q_Sg ,      Q_Sw ,    Q_Teet , &
                                   Q_TFA1 ,    Q_TFA2 ,    Q_TFrl ,    Q_TSS1 ,    Q_TSS2 ,       Q_Y ,     Q_Yaw , &
                                  RFrlBrM ,  TipRDxb1 ,  TipRDxb2 ,  TipRDxb3 ,  RootFxb1 ,  RootFxb2 ,  RootFxb3 , &
                                 RootFxc1 ,  RootFxc2 ,  RootFxc3 ,  RootFyb1 ,  RootFyb2 ,  RootFyb3 ,  RootFyc1 , &
                                 RootFyc2 ,  RootFyc3 ,  RootFzc1 ,  RootFzc2 ,  RootFzc3 ,  RootFzc1 ,  RootFzc2 , &
                                 RootFzc3 ,  RootMxb1 ,  RootMxb2 ,  RootMxb3 ,  RootMyb1 ,  RootMyb2 ,  RootMyb3 , &
                                 RootMxc1 ,  RootMxc2 ,  RootMxc3 ,  RootMyc1 ,  RootMyc2 ,  RootMyc3 ,  RootMxb1 , &
                                 RootMxb2 ,  RootMxb3 ,  RootMxc1 ,  RootMxc2 ,  RootMxc3 ,  RootMyb1 ,  RootMyb2 , &
                                 RootMyb3 ,  RootMyc1 ,  RootMyc2 ,  RootMyc3 ,  RootMzc1 ,  RootMzc2 ,  RootMzc3 , &
                                 RootMzc1 ,  RootMzc2 ,  RootMzc3 , LSSTipAxa ,  RotFurlP ,  RotFurlA ,  RotFurlP , &
                                 RotFurlV ,    RotPwr , LSSTipVxa ,   TeetAya ,   TeetPya ,   TeetVya , LSShftFxa , &
                                LSShftMxa , Spn1ALxb1 , Spn1ALxb2 , Spn1ALxb3 , Spn1ALyb1 , Spn1ALyb2 , Spn1ALyb3 , &
                                Spn1ALzb1 , Spn1ALzb2 , Spn1ALzb3 , Spn1FLxb1 , Spn1FLxb2 , Spn1FLxb3 , Spn1FLyb1 , &
                                Spn1FLyb2 , Spn1FLyb3 , Spn1FLzb1 , Spn1FLzb2 , Spn1FLzb3 , Spn1MLxb1 , Spn1MLxb2 , &
                                Spn1MLxb3 , Spn1MLyb1 , Spn1MLyb2 , Spn1MLyb3 , Spn1MLzb1 , Spn1MLzb2 , Spn1MLzb3 , &
                                Spn1RDxb1 , Spn1RDxb2 , Spn1RDxb3 , Spn1RDyb1 , Spn1RDyb2 , Spn1RDyb3 , Spn1RDzb1 , &
                                Spn1RDzb2 , Spn1RDzb3 , Spn1TDxb1 , Spn1TDxb2 , Spn1TDxb3 , Spn1TDyb1 , Spn1TDyb2 , &
                                Spn1TDyb3 , Spn1TDzb1 , Spn1TDzb2 , Spn1TDzb3 , Spn2ALxb1 , Spn2ALxb2 , Spn2ALxb3 , &
                                Spn2ALyb1 , Spn2ALyb2 , Spn2ALyb3 , Spn2ALzb1 , Spn2ALzb2 , Spn2ALzb3 , Spn2FLxb1 , &
                                Spn2FLxb2 , Spn2FLxb3 , Spn2FLyb1 , Spn2FLyb2 , Spn2FLyb3 , Spn2FLzb1 , Spn2FLzb2 , &
                                Spn2FLzb3 , Spn2MLxb1 , Spn2MLxb2 , Spn2MLxb3 , Spn2MLyb1 , Spn2MLyb2 , Spn2MLyb3 , &
                                Spn2MLzb1 , Spn2MLzb2 , Spn2MLzb3 , Spn2RDxb1 , Spn2RDxb2 , Spn2RDxb3 , Spn2RDyb1 , &
                                Spn2RDyb2 , Spn2RDyb3 , Spn2RDzb1 , Spn2RDzb2 , Spn2RDzb3 , Spn2TDxb1 , Spn2TDxb2 , &
                                Spn2TDxb3 , Spn2TDyb1 , Spn2TDyb2 , Spn2TDyb3 , Spn2TDzb1 , Spn2TDzb2 , Spn2TDzb3 , &
                                Spn3ALxb1 , Spn3ALxb2 , Spn3ALxb3 , Spn3ALyb1 , Spn3ALyb2 , Spn3ALyb3 , Spn3ALzb1 , &
                                Spn3ALzb2 , Spn3ALzb3 , Spn3FLxb1 , Spn3FLxb2 , Spn3FLxb3 , Spn3FLyb1 , Spn3FLyb2 , &
                                Spn3FLyb3 , Spn3FLzb1 , Spn3FLzb2 , Spn3FLzb3 , Spn3MLxb1 , Spn3MLxb2 , Spn3MLxb3 , &
                                Spn3MLyb1 , Spn3MLyb2 , Spn3MLyb3 , Spn3MLzb1 , Spn3MLzb2 , Spn3MLzb3 , Spn3RDxb1 , &
                                Spn3RDxb2 , Spn3RDxb3 , Spn3RDyb1 , Spn3RDyb2 , Spn3RDyb3 , Spn3RDzb1 , Spn3RDzb2 , &
                                Spn3RDzb3 , Spn3TDxb1 , Spn3TDxb2 , Spn3TDxb3 , Spn3TDyb1 , Spn3TDyb2 , Spn3TDyb3 , &
                                Spn3TDzb1 , Spn3TDzb2 , Spn3TDzb3 , Spn4ALxb1 , Spn4ALxb2 , Spn4ALxb3 , Spn4ALyb1 , &
                                Spn4ALyb2 , Spn4ALyb3 , Spn4ALzb1 , Spn4ALzb2 , Spn4ALzb3 , Spn4FLxb1 , Spn4FLxb2 , &
                                Spn4FLxb3 , Spn4FLyb1 , Spn4FLyb2 , Spn4FLyb3 , Spn4FLzb1 , Spn4FLzb2 , Spn4FLzb3 , &
                                Spn4MLxb1 , Spn4MLxb2 , Spn4MLxb3 , Spn4MLyb1 , Spn4MLyb2 , Spn4MLyb3 , Spn4MLzb1 , &
                                Spn4MLzb2 , Spn4MLzb3 , Spn4RDxb1 , Spn4RDxb2 , Spn4RDxb3 , Spn4RDyb1 , Spn4RDyb2 , &
                                Spn4RDyb3 , Spn4RDzb1 , Spn4RDzb2 , Spn4RDzb3 , Spn4TDxb1 , Spn4TDxb2 , Spn4TDxb3 , &
                                Spn4TDyb1 , Spn4TDyb2 , Spn4TDyb3 , Spn4TDzb1 , Spn4TDzb2 , Spn4TDzb3 , Spn5ALxb1 , &
                                Spn5ALxb2 , Spn5ALxb3 , Spn5ALyb1 , Spn5ALyb2 , Spn5ALyb3 , Spn5ALzb1 , Spn5ALzb2 , &
                                Spn5ALzb3 , Spn5FLxb1 , Spn5FLxb2 , Spn5FLxb3 , Spn5FLyb1 , Spn5FLyb2 , Spn5FLyb3 , &
                                Spn5FLzb1 , Spn5FLzb2 , Spn5FLzb3 , Spn5MLxb1 , Spn5MLxb2 , Spn5MLxb3 , Spn5MLyb1 , &
                                Spn5MLyb2 , Spn5MLyb3 , Spn5MLzb1 , Spn5MLzb2 , Spn5MLzb3 , Spn5RDxb1 , Spn5RDxb2 , &
                                Spn5RDxb3 , Spn5RDyb1 , Spn5RDyb2 , Spn5RDyb3 , Spn5RDzb1 , Spn5RDzb2 , Spn5RDzb3 , &
                                Spn5TDxb1 , Spn5TDxb2 , Spn5TDxb3 , Spn5TDyb1 , Spn5TDyb2 , Spn5TDyb3 , Spn5TDzb1 , &
                                Spn5TDzb2 , Spn5TDzb3 , Spn6ALxb1 , Spn6ALxb2 , Spn6ALxb3 , Spn6ALyb1 , Spn6ALyb2 , &
                                Spn6ALyb3 , Spn6ALzb1 , Spn6ALzb2 , Spn6ALzb3 , Spn6FLxb1 , Spn6FLxb2 , Spn6FLxb3 , &
                                Spn6FLyb1 , Spn6FLyb2 , Spn6FLyb3 , Spn6FLzb1 , Spn6FLzb2 , Spn6FLzb3 , Spn6MLxb1 , &
                                Spn6MLxb2 , Spn6MLxb3 , Spn6MLyb1 , Spn6MLyb2 , Spn6MLyb3 , Spn6MLzb1 , Spn6MLzb2 , &
                                Spn6MLzb3 , Spn6RDxb1 , Spn6RDxb2 , Spn6RDxb3 , Spn6RDyb1 , Spn6RDyb2 , Spn6RDyb3 , &
                                Spn6RDzb1 , Spn6RDzb2 , Spn6RDzb3 , Spn6TDxb1 , Spn6TDxb2 , Spn6TDxb3 , Spn6TDyb1 , &
                                Spn6TDyb2 , Spn6TDyb3 , Spn6TDzb1 , Spn6TDzb2 , Spn6TDzb3 , Spn7ALxb1 , Spn7ALxb2 , &
                                Spn7ALxb3 , Spn7ALyb1 , Spn7ALyb2 , Spn7ALyb3 , Spn7ALzb1 , Spn7ALzb2 , Spn7ALzb3 , &
                                Spn7FLxb1 , Spn7FLxb2 , Spn7FLxb3 , Spn7FLyb1 , Spn7FLyb2 , Spn7FLyb3 , Spn7FLzb1 , &
                                Spn7FLzb2 , Spn7FLzb3 , Spn7MLxb1 , Spn7MLxb2 , Spn7MLxb3 , Spn7MLyb1 , Spn7MLyb2 , &
                                Spn7MLyb3 , Spn7MLzb1 , Spn7MLzb2 , Spn7MLzb3 , Spn7RDxb1 , Spn7RDxb2 , Spn7RDxb3 , &
                                Spn7RDyb1 , Spn7RDyb2 , Spn7RDyb3 , Spn7RDzb1 , Spn7RDzb2 , Spn7RDzb3 , Spn7TDxb1 , &
                                Spn7TDxb2 , Spn7TDxb3 , Spn7TDyb1 , Spn7TDyb2 , Spn7TDyb3 , Spn7TDzb1 , Spn7TDzb2 , &
                                Spn7TDzb3 , Spn8ALxb1 , Spn8ALxb2 , Spn8ALxb3 , Spn8ALyb1 , Spn8ALyb2 , Spn8ALyb3 , &
                                Spn8ALzb1 , Spn8ALzb2 , Spn8ALzb3 , Spn8FLxb1 , Spn8FLxb2 , Spn8FLxb3 , Spn8FLyb1 , &
                                Spn8FLyb2 , Spn8FLyb3 , Spn8FLzb1 , Spn8FLzb2 , Spn8FLzb3 , Spn8MLxb1 , Spn8MLxb2 , &
                                Spn8MLxb3 , Spn8MLyb1 , Spn8MLyb2 , Spn8MLyb3 , Spn8MLzb1 , Spn8MLzb2 , Spn8MLzb3 , &
                                Spn8RDxb1 , Spn8RDxb2 , Spn8RDxb3 , Spn8RDyb1 , Spn8RDyb2 , Spn8RDyb3 , Spn8RDzb1 , &
                                Spn8RDzb2 , Spn8RDzb3 , Spn8TDxb1 , Spn8TDxb2 , Spn8TDxb3 , Spn8TDyb1 , Spn8TDyb2 , &
                                Spn8TDyb3 , Spn8TDzb1 , Spn8TDzb2 , Spn8TDzb3 , Spn9ALxb1 , Spn9ALxb2 , Spn9ALxb3 , &
                                Spn9ALyb1 , Spn9ALyb2 , Spn9ALyb3 , Spn9ALzb1 , Spn9ALzb2 , Spn9ALzb3 , Spn9FLxb1 , &
                                Spn9FLxb2 , Spn9FLxb3 , Spn9FLyb1 , Spn9FLyb2 , Spn9FLyb3 , Spn9FLzb1 , Spn9FLzb2 , &
                                Spn9FLzb3 , Spn9MLxb1 , Spn9MLxb2 , Spn9MLxb3 , Spn9MLyb1 , Spn9MLyb2 , Spn9MLyb3 , &
                                Spn9MLzb1 , Spn9MLzb2 , Spn9MLzb3 , Spn9RDxb1 , Spn9RDxb2 , Spn9RDxb3 , Spn9RDyb1 , &
                                Spn9RDyb2 , Spn9RDyb3 , Spn9RDzb1 , Spn9RDzb2 , Spn9RDzb3 , Spn9TDxb1 , Spn9TDxb2 , &
                                Spn9TDxb3 , Spn9TDyb1 , Spn9TDyb2 , Spn9TDyb3 , Spn9TDzb1 , Spn9TDzb2 , Spn9TDzb3 , &
                                TailFurlP , TailFurlA , TailFurlP , TailFurlV ,   TeetAya ,   TeetPya ,   TeetPya , &
                                  TeetVya ,   TFrlBrM , TipClrnc1 , TipClrnc2 , TipClrnc3 ,  TipALxb1 ,  TipALxb2 , &
                                 TipALxb3 ,  TipALyb1 ,  TipALyb2 ,  TipALyb3 ,  TipALzb1 ,  TipALzb2 ,  TipALzb3 , &
                                TipClrnc1 , TipClrnc2 , TipClrnc3 ,   TipDxb1 ,   TipDxb2 ,   TipDxb3 ,   TipDxc1 , &
                                  TipDxc2 ,   TipDxc3 ,   TipDyb1 ,   TipDyb2 ,   TipDyb3 ,   TipDyc1 ,   TipDyc2 , &
                                  TipDyc3 ,   TipDzc1 ,   TipDzc2 ,   TipDzc3 ,   TipDzc1 ,   TipDzc2 ,   TipDzc3 , &
                                 TipRDxb1 ,  TipRDxb2 ,  TipRDxb3 ,  TipRDyb1 ,  TipRDyb2 ,  TipRDyb3 ,  TipRDzc1 , &
                                 TipRDzc2 ,  TipRDzc3 ,  TipRDzc1 ,  TipRDzc2 ,  TipRDzc3 , YawBrTDzt , YawBrTDxt , &
                                YawBrRDyt , YawBrRDxt , YawBrTDyt , YawBrRDzt , TwHt1ALxt , TwHt1ALyt , TwHt1ALzt , &
                                TwHt1FLxt , TwHt1FLyt , TwHt1FLzt , TwHt1MLxt , TwHt1MLyt , TwHt1MLzt , TwHt1RDxt , &
                                TwHt1RDyt , TwHt1RDzt , TwHt1RPxi , TwHt1RPyi , TwHt1RPzi , TwHt1TDxt , TwHt1TDyt , &
                                TwHt1TDzt , TwHt1TPxi , TwHt1TPyi , TwHt1TPzi , TwHt2ALxt , TwHt2ALyt , TwHt2ALzt , &
                                TwHt2FLxt , TwHt2FLyt , TwHt2FLzt , TwHt2MLxt , TwHt2MLyt , TwHt2MLzt , TwHt2RDxt , &
                                TwHt2RDyt , TwHt2RDzt , TwHt2RPxi , TwHt2RPyi , TwHt2RPzi , TwHt2TDxt , TwHt2TDyt , &
                                TwHt2TDzt , TwHt2TPxi , TwHt2TPyi , TwHt2TPzi , TwHt3ALxt , TwHt3ALyt , TwHt3ALzt , &
                                TwHt3FLxt , TwHt3FLyt , TwHt3FLzt , TwHt3MLxt , TwHt3MLyt , TwHt3MLzt , TwHt3RDxt , &
                                TwHt3RDyt , TwHt3RDzt , TwHt3RPxi , TwHt3RPyi , TwHt3RPzi , TwHt3TDxt , TwHt3TDyt , &
                                TwHt3TDzt , TwHt3TPxi , TwHt3TPyi , TwHt3TPzi , TwHt4ALxt , TwHt4ALyt , TwHt4ALzt , &
                                TwHt4FLxt , TwHt4FLyt , TwHt4FLzt , TwHt4MLxt , TwHt4MLyt , TwHt4MLzt , TwHt4RDxt , &
                                TwHt4RDyt , TwHt4RDzt , TwHt4RPxi , TwHt4RPyi , TwHt4RPzi , TwHt4TDxt , TwHt4TDyt , &
                                TwHt4TDzt , TwHt4TPxi , TwHt4TPyi , TwHt4TPzi , TwHt5ALxt , TwHt5ALyt , TwHt5ALzt , &
                                TwHt5FLxt , TwHt5FLyt , TwHt5FLzt , TwHt5MLxt , TwHt5MLyt , TwHt5MLzt , TwHt5RDxt , &
                                TwHt5RDyt , TwHt5RDzt , TwHt5RPxi , TwHt5RPyi , TwHt5RPzi , TwHt5TDxt , TwHt5TDyt , &
                                TwHt5TDzt , TwHt5TPxi , TwHt5TPyi , TwHt5TPzi , TwHt6ALxt , TwHt6ALyt , TwHt6ALzt , &
                                TwHt6FLxt , TwHt6FLyt , TwHt6FLzt , TwHt6MLxt , TwHt6MLyt , TwHt6MLzt , TwHt6RDxt , &
                                TwHt6RDyt , TwHt6RDzt , TwHt6RPxi , TwHt6RPyi , TwHt6RPzi , TwHt6TDxt , TwHt6TDyt , &
                                TwHt6TDzt , TwHt6TPxi , TwHt6TPyi , TwHt6TPzi , TwHt7ALxt , TwHt7ALyt , TwHt7ALzt , &
                                TwHt7FLxt , TwHt7FLyt , TwHt7FLzt , TwHt7MLxt , TwHt7MLyt , TwHt7MLzt , TwHt7RDxt , &
                                TwHt7RDyt , TwHt7RDzt , TwHt7RPxi , TwHt7RPyi , TwHt7RPzi , TwHt7TDxt , TwHt7TDyt , &
                                TwHt7TDzt , TwHt7TPxi , TwHt7TPyi , TwHt7TPzi , TwHt8ALxt , TwHt8ALyt , TwHt8ALzt , &
                                TwHt8FLxt , TwHt8FLyt , TwHt8FLzt , TwHt8MLxt , TwHt8MLyt , TwHt8MLzt , TwHt8RDxt , &
                                TwHt8RDyt , TwHt8RDzt , TwHt8RPxi , TwHt8RPyi , TwHt8RPzi , TwHt8TDxt , TwHt8TDyt , &
                                TwHt8TDzt , TwHt8TPxi , TwHt8TPyi , TwHt8TPzi , TwHt9ALxt , TwHt9ALyt , TwHt9ALzt , &
                                TwHt9FLxt , TwHt9FLyt , TwHt9FLzt , TwHt9MLxt , TwHt9MLyt , TwHt9MLzt , TwHt9RDxt , &
                                TwHt9RDyt , TwHt9RDzt , TwHt9RPxi , TwHt9RPyi , TwHt9RPzi , TwHt9TDxt , TwHt9TDyt , &
                                TwHt9TDzt , TwHt9TPxi , TwHt9TPyi , TwHt9TPzi ,  TwrBsFxt ,  TwrBsFyt ,  TwrBsFzt , &
                                 TwrBsMxt ,  TwrBsMyt ,  TwrBsMzt , TipClrnc1 , TipClrnc2 , TipClrnc3 ,  TipRDzc1 , &
                                 TipRDzc2 ,  TipRDzc3 ,    YawAzn ,    YawAzn ,    YawAzn ,  YawBrFxn ,  YawBrFxp , &
                                 YawBrFyn ,  YawBrFyp ,  YawBrFzn ,  YawBrFzn ,  YawBrMxn ,  YawBrMxp ,  YawBrMyn , &
                                 YawBrMyp ,  YawBrMzn ,  YawBrMzn , YawBrRAxp , YawBrRAyp , YawBrRAzp , YawBrRDxt , &
                                YawBrRDyt , YawBrRDzt , YawBrRVxp , YawBrRVyp , YawBrRVzp , YawBrTAxp , YawBrTAyp , &
                                YawBrTAzp , YawBrTDxp , YawBrTDxt , YawBrTDyp , YawBrTDyt , YawBrTDzp , YawBrTDzt , &
                                   YawPzn ,    YawPzn ,    YawPzn ,    YawVzn ,    YawVzn ,    YawVzn /)
   CHARACTER(ChanLen), PARAMETER :: ParamUnitsAry(972) =  (/ &                     ! This lists the units corresponding to the allowed parameters
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg/s^2) ","(rpm)     ","(kN·m)    ","(deg/s^2) ","(kW)      ","(kN·m)    ","(rpm)     ", &
                               "(m)       ","(m)       ","(m)       ","(deg/s^2) ","(deg/s^2) ","(deg/s^2) ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ","(deg)     ","(deg)     ", &
                               "(rpm)     ","(rpm)     ","(rpm)     ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kW)      ","(kN·m)    ","(deg/s^2) ", &
                               "(deg/s^2) ","(deg/s^2) ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(rpm)     ","(rpm)     ","(rpm)     ","(deg)     ","(deg/s^2) ", &
                               "(deg)     ","(deg/s)   ","(deg/s^2) ","(deg/s^2) ","(deg/s^2) ","(deg/s)   ","(deg/s)   ", &
                               "(deg/s)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s)     ","(m/s)     ","(m/s)     ", &
                               "(m)       ","(m)       ","(m)       ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(deg)     ", &
                               "(deg/s^2) ","(deg/s^2) ","(deg/s^2) ","(deg/s^2) ","(deg/s^2) ","(deg/s^2) ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg/s)   ","(deg/s)   ","(deg/s)   ","(deg/s)   ", &
                               "(deg/s)   ","(deg/s)   ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m/s)     ","(m/s)     ","(m/s)     ","(m/s)     ","(m/s)     ", &
                               "(m/s)     ","(deg)     ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(rad/s^2) ","(rad/s^2) ","(m/s^2)   ", &
                               "(rad/s^2) ","(rad/s^2) ","(rad/s^2) ","(m/s^2)   ","(m/s^2)   ","(rad/s^2) ","(m/s^2)   ", &
                               "(m/s^2)   ","(rad/s^2) ","(m/s^2)   ","(m/s^2)   ","(rad/s^2) ","(rad/s^2) ","(m/s)     ", &
                               "(m/s)     ","(m/s)     ","(m/s)     ","(m/s)     ","(m/s)     ","(m/s)     ","(m/s)     ", &
                               "(m/s)     ","(rad/s)   ","(rad/s)   ","(m/s)     ","(rad/s)   ","(rad/s)   ","(rad/s)   ", &
                               "(m/s)     ","(m/s)     ","(rad/s)   ","(m/s)     ","(m/s)     ","(rad/s)   ","(m/s)     ", &
                               "(m/s)     ","(rad/s)   ","(rad/s)   ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(rad)     ","(rad)     ", &
                               "(m)       ","(rad)     ","(rad)     ","(rad)     ","(m)       ","(m)       ","(rad)     ", &
                               "(m)       ","(m)       ","(rad)     ","(m)       ","(m)       ","(rad)     ","(rad)     ", &
                               "(kN·m)    ","(deg)     ","(deg)     ","(deg)     ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg/s^2) ","(deg)     ","(deg/s^2) ","(deg)     ", &
                               "(deg/s)   ","(kW)      ","(rpm)     ","(deg/s^2) ","(deg)     ","(deg/s)   ","(kN)      ", &
                               "(kN·m)    ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(deg)     ","(deg/s^2) ","(deg)     ","(deg/s)   ","(deg/s^2) ","(deg)     ","(deg)     ", &
                               "(deg/s)   ","(kN·m)    ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(deg)     ","(deg)     ","(m)       ","(deg)     ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(m/s^2)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg)     ","(deg)     ","(m)       ","(m)       ", &
                               "(m)       ","(m)       ","(m)       ","(m)       ","(kN)      ","(kN)      ","(kN)      ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(m)       ","(m)       ","(m)       ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg/s^2) ","(deg/s^2) ","(deg/s^2) ","(kN)      ","(kN)      ", &
                               "(kN)      ","(kN)      ","(kN)      ","(kN)      ","(kN·m)    ","(kN·m)    ","(kN·m)    ", &
                               "(kN·m)    ","(kN·m)    ","(kN·m)    ","(deg/s^2) ","(deg/s^2) ","(deg/s^2) ","(deg)     ", &
                               "(deg)     ","(deg)     ","(deg/s)   ","(deg/s)   ","(deg/s)   ","(m/s^2)   ","(m/s^2)   ", &
                               "(m/s^2)   ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ","(m)       ", &
                               "(deg)     ","(deg)     ","(deg)     ","(deg/s)   ","(deg/s)   ","(deg/s)   "/)


      ! Initialize values
   ErrStat = ErrID_None
   ErrMsg = ""
   InvalidOutput = .FALSE.


!   ..... Developer must add checking for invalid inputs here: .....
if (p%BD4Blades) then
   startIndx = 1
else
   startIndx = p%NumBl+1
end if

   DO I = startIndx,3  ! Invalid blades

         ! motions

      InvalidOutput(   TipDxc(  I) ) = .TRUE.
      InvalidOutput(   TipDyc(  I) ) = .TRUE.
      InvalidOutput(   TipDzc(  I) ) = .TRUE.
      InvalidOutput(   TipDxb(  I) ) = .TRUE.
      InvalidOutput(   TipDyb(  I) ) = .TRUE.
      InvalidOutput(  TipALxb(  I) ) = .TRUE.
      InvalidOutput(  TipALyb(  I) ) = .TRUE.
      InvalidOutput(  TipALzb(  I) ) = .TRUE.
      InvalidOutput(  TipRDxb(  I) ) = .TRUE.
      InvalidOutput(  TipRDyb(  I) ) = .TRUE.
      InvalidOutput(  TipRDzc(  I) ) = .TRUE.
      InvalidOutput( TipClrnc(  I) ) = .TRUE.

         ! loads

      InvalidOutput(  RootFxc(  I) ) = .TRUE.
      InvalidOutput(  RootFyc(  I) ) = .TRUE.
      InvalidOutput(  RootFzc(  I) ) = .TRUE.
      InvalidOutput(  RootFxb(  I) ) = .TRUE.
      InvalidOutput(  RootFyb(  I) ) = .TRUE.
      InvalidOutput(  RootMxc(  I) ) = .TRUE.
      InvalidOutput(  RootMyc(  I) ) = .TRUE.
      InvalidOutput(  RootMzc(  I) ) = .TRUE.
      InvalidOutput(  RootMxb(  I) ) = .TRUE.
      InvalidOutput(  RootMyb(  I) ) = .TRUE.

         ! Blade node motions

      InvalidOutput(  SpnALxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnALyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnALzb(:,I) ) = .TRUE.

      InvalidOutput(  SpnTDxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnTDyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnTDzb(:,I) ) = .TRUE.

      InvalidOutput(  SpnRDxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnRDyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnRDzb(:,I) ) = .TRUE.

         ! Blade node loads

      InvalidOutput(  SpnMLxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnMLyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnMLzb(:,I) ) = .TRUE.

      InvalidOutput(  SpnFLxb(:,I) ) = .TRUE.
      InvalidOutput(  SpnFLyb(:,I) ) = .TRUE.
      InvalidOutput(  SpnFLzb(:,I) ) = .TRUE.

   END DO


   DO I = 1,p%NumBl

      DO J = p%NBlGages+1,9 ! Invalid blade gages

         InvalidOutput(  SpnALxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnALyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnALzb(J,I) ) = .TRUE.

         InvalidOutput(  SpnTDxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnTDyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnTDzb(J,I) ) = .TRUE.

         InvalidOutput(  SpnRDxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnRDyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnRDzb(J,I) ) = .TRUE.

            ! Loads

         InvalidOutput(  SpnMLxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnMLyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnMLzb(J,I) ) = .TRUE.

         InvalidOutput(  SpnFLxb(J,I) ) = .TRUE.
         InvalidOutput(  SpnFLyb(J,I) ) = .TRUE.
         InvalidOutput(  SpnFLzb(J,I) ) = .TRUE.


      END DO !J

   END DO !I

   DO J = p%NTwGages+1,9 !Invalid tower gages

         ! Motions

      InvalidOutput( TwHtALxt(J) ) = .TRUE.
      InvalidOutput( TwHtALyt(J) ) = .TRUE.
      InvalidOutput( TwHtALzt(J) ) = .TRUE.

      InvalidOutput( TwHtTDxt(J) ) = .TRUE.
      InvalidOutput( TwHtTDyt(J) ) = .TRUE.
      InvalidOutput( TwHtTDzt(J) ) = .TRUE.

      InvalidOutput( TwHtRDxt(J) ) = .TRUE.
      InvalidOutput( TwHtRDyt(J) ) = .TRUE.
      InvalidOutput( TwHtRDzt(J) ) = .TRUE.

      InvalidOutput( TwHtTPxi(J) ) = .TRUE.
      InvalidOutput( TwHtTPyi(J) ) = .TRUE.
      InvalidOutput( TwHtTPzi(J) ) = .TRUE.

      InvalidOutput( TwHtRPxi(J) ) = .TRUE.
      InvalidOutput( TwHtRPyi(J) ) = .TRUE.
      InvalidOutput( TwHtRPzi(J) ) = .TRUE.

         ! Loads

      InvalidOutput( TwHtMLxt(J) ) = .TRUE.
      InvalidOutput( TwHtMLyt(J) ) = .TRUE.
      InvalidOutput( TwHtMLzt(J) ) = .TRUE.

      InvalidOutput( TwHtFLxt(J) ) = .TRUE.
      InvalidOutput( TwHtFLyt(J) ) = .TRUE.
      InvalidOutput( TwHtFLzt(J) ) = .TRUE.

   END DO

   IF ( p%NumBl < 3_IntKi ) THEN
      InvalidOutput(PtchPMzc3) = .TRUE.

      InvalidOutput(   Q_B3E1) = .TRUE.
      InvalidOutput(   Q_B3F1) = .TRUE.
      InvalidOutput(   Q_B3F2) = .TRUE.

      InvalidOutput(  QD_B3E1) = .TRUE.
      InvalidOutput(  QD_B3F1) = .TRUE.
      InvalidOutput(  QD_B3F2) = .TRUE.

      InvalidOutput( QD2_B3E1) = .TRUE.
      InvalidOutput( QD2_B3F1) = .TRUE.
      InvalidOutput( QD2_B3F2) = .TRUE.
   ELSE IF ( p%NumBl > 2_IntKi ) THEN
      InvalidOutput(  TeetPya) = .TRUE.
      InvalidOutput(  TeetVya) = .TRUE.
      InvalidOutput(  TeetAya) = .TRUE.

      InvalidOutput(   Q_Teet) = .TRUE.
      InvalidOutput(  QD_Teet) = .TRUE.
      InvalidOutput( QD2_Teet) = .TRUE.
   END IF
   
   InvalidOutput(HSSBrTq) = p%method == Method_RK4

    IF ( p%BD4Blades ) THEN
      InvalidOutput(   Q_B1E1) = .TRUE.
      InvalidOutput(   Q_B1F1) = .TRUE.
      InvalidOutput(   Q_B1F2) = .TRUE.

      InvalidOutput(  QD_B1E1) = .TRUE.
      InvalidOutput(  QD_B1F1) = .TRUE.
      InvalidOutput(  QD_B1F2) = .TRUE.

      InvalidOutput( QD2_B1E1) = .TRUE.
      InvalidOutput( QD2_B1F1) = .TRUE.
      InvalidOutput( QD2_B1F2) = .TRUE.      
      
      InvalidOutput(   Q_B2E1) = .TRUE.
      InvalidOutput(   Q_B2F1) = .TRUE.
      InvalidOutput(   Q_B2F2) = .TRUE.

      InvalidOutput(  QD_B2E1) = .TRUE.
      InvalidOutput(  QD_B2F1) = .TRUE.
      InvalidOutput(  QD_B2F2) = .TRUE.

      InvalidOutput( QD2_B2E1) = .TRUE.
      InvalidOutput( QD2_B2F1) = .TRUE.
      InvalidOutput( QD2_B2F2) = .TRUE.
      
      InvalidOutput(   Q_B3E1) = .TRUE.
      InvalidOutput(   Q_B3F1) = .TRUE.
      InvalidOutput(   Q_B3F2) = .TRUE.

      InvalidOutput(  QD_B3E1) = .TRUE.
      InvalidOutput(  QD_B3F1) = .TRUE.
      InvalidOutput(  QD_B3F2) = .TRUE.

      InvalidOutput( QD2_B3E1) = .TRUE.
      InvalidOutput( QD2_B3F1) = .TRUE.
      InvalidOutput( QD2_B3F2) = .TRUE.      
   END IF
!   ................. End of validity checking .................


   !-------------------------------------------------------------------------------------------------
   ! Allocate and set index, name, and units for the output channels
   ! If a selected output channel is not available in this module, set error flag.
   !-------------------------------------------------------------------------------------------------

   ALLOCATE ( p%OutParam(0:p%NumOuts) , STAT=ErrStat )
   IF ( ErrStat /= 0_IntKi )  THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = "SetOutParam:Error allocating memory for the ElastoDyn OutParam array."
      RETURN
   ELSE
      ErrStat = ErrID_None
   ENDIF

      ! Set index, name, and units for the time output channel:

   p%OutParam(0)%Indx  = Time
   p%OutParam(0)%Name  = "Time"    ! OutParam(0) is the time channel by default.
   p%OutParam(0)%Units = "(s)"
   p%OutParam(0)%SignM = 1


      ! Set index, name, and units for all of the output channels.
      ! If a selected output channel is not available by this module set ErrStat = ErrID_Warn.

   DO I = 1,p%NumOuts

      p%OutParam(I)%Name  = OutList(I)
      OutListTmp          = OutList(I)

      ! Reverse the sign (+/-) of the output channel if the user prefixed the
      !   channel name with a "-", "_", "m", or "M" character indicating "minus".


      CheckOutListAgain = .FALSE.

      IF      ( INDEX( "-_", OutListTmp(1:1) ) > 0 ) THEN
         p%OutParam(I)%SignM = -1                         ! ex, "-TipDxc1" causes the sign of TipDxc1 to be switched.
         OutListTmp          = OutListTmp(2:)
      ELSE IF ( INDEX( "mM", OutListTmp(1:1) ) > 0 ) THEN ! We'll assume this is a variable name for now, (if not, we will check later if OutListTmp(2:) is also a variable name)
         CheckOutListAgain   = .TRUE.
         p%OutParam(I)%SignM = 1
      ELSE
         p%OutParam(I)%SignM = 1
      END IF

      CALL Conv2UC( OutListTmp )    ! Convert OutListTmp to upper case


      Indx = IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )


         ! If it started with an "M" (CheckOutListAgain) we didn't find the value in our list (Indx < 1)

      IF ( CheckOutListAgain .AND. Indx < 1 ) THEN    ! Let's assume that "M" really meant "minus" and then test again
         p%OutParam(I)%SignM = -1                     ! ex, "MTipDxc1" causes the sign of TipDxc1 to be switched.
         OutListTmp          = OutListTmp(2:)

         Indx = IndexCharAry( OutListTmp(1:OutStrLenM1), ValidParamAry )
      END IF


      IF ( Indx > 0 ) THEN ! we found the channel name
         p%OutParam(I)%Indx     = ParamIndxAry(Indx)
         IF ( InvalidOutput( ParamIndxAry(Indx) ) ) THEN  ! but, it isn't valid for these settings
            p%OutParam(I)%Units = "INVALID"
            p%OutParam(I)%SignM = 0
         ELSE
            p%OutParam(I)%Units = ParamUnitsAry(Indx) ! it's a valid output
         END IF
      ELSE ! this channel isn't valid
         p%OutParam(I)%Indx  = Time                 ! pick any valid channel (I just picked "Time" here because it's universal)
         p%OutParam(I)%Units = "INVALID"
         p%OutParam(I)%SignM = 0                    ! multiply all results by zero

         ErrStat = ErrID_Fatal
         ErrMsg  = "SetOutParam:"//trim(p%OutParam(I)%Name)//" is not an available output channel. "//TRIM(ErrMsg)
      END IF

   END DO

   RETURN
END SUBROUTINE SetOutParam
!----------------------------------------------------------------------------------------------------------------------------------
!End of code generated by Matlab script
!**********************************************************************************************************************************
SUBROUTINE Coeff(p,InputFileData, ErrStat, ErrMsg)
! This routine is used to compute rotor (blade and hub) properties:
!   KBF(), KBE(), CBF(), CBE(), FreqBF(), FreqBE(), AxRedBld(),
!   TwistedSF(), BldMass(), FirstMom(), SecondMom(), BldCG(),
!   RotMass, RotIner, Hubg1Iner, Hubg2Iner, rSAerCenn1(), and
!   rSAerCenn2(), BElmtMass()
! tower properties:
!   KTFA(), KTSS(), CTFA(), CTSS(), FreqTFA(), FreqTSS(),
!   AxRedTFA(), AxRedTSS(), TwrFASF(), TwrSSSF(), TwrMass, and
!   TwrTpMass, TElmtMass()
! structure that furls with the rotor (not including rotor) properties:
!   RrfaIner
! tail boom properties:
!   AtfaIner
! and nacelle properties:
!   Nacd2Iner
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed variables

   TYPE(ED_ParameterType),        INTENT(INOUT)    :: p                             ! Parameters of the structural dynamics module
   TYPE(ED_InputFile),            INTENT(IN)       :: InputFileData                 ! all the data in the ElastoDyn input file
   INTEGER(IntKi),                INTENT(OUT)      :: ErrStat                       ! Error status
   CHARACTER(*),                  INTENT(OUT)      :: ErrMsg                        ! Error message when ErrStat =/ ErrID_None


      ! Local variables.

   REAL(ReKi)                   :: AxRdBld   (3,3)                                 ! Temporary result holding the current addition to the p%AxRedBld() array.
   REAL(ReKi)                   :: AxRdBldOld(3,3)                                 ! Previous AxRdBld (i.e., AxRdBld from the previous node)
   REAL(ReKi)                   :: AxRdTFA   (2,2)                                 ! Temporary result holding the current addition to the AxRedTFA() array.
   REAL(ReKi)                   :: AxRdTFAOld(2,2)                                 ! Previous AxRdTFA (i.e., AxRdTFA from the previous node)
   REAL(ReKi)                   :: AxRdTSS   (2,2)                                 ! Temporary result holding the current addition to the AxRedTSS() array.
   REAL(ReKi)                   :: AxRdTSSOld(2,2)                                 ! Previous AxRdTSS (i.e., AxRdTSS from the previous node)
   REAL(ReKi)                   :: TmpDist                                         ! Temporary distance used in the calculation of the aero center locations.
   REAL(ReKi)                   :: TmpDistj1                                       ! Temporary distance used in the calculation of the aero center locations.
   REAL(ReKi)                   :: TmpDistj2                                       ! Temporary distance used in the calculation of the aero center locations.
   REAL(ReKi)                   :: ElmntStff                                       ! (Temporary) stiffness of an element.
   REAL(ReKi)                   :: ElStffFA                                        ! (Temporary) tower fore-aft stiffness of an element
   REAL(ReKi)                   :: ElStffSS                                        ! (Temporary) tower side-to-side  stiffness of an element
   REAL(ReKi)                   :: FMomAbvNd (p%NumBl,p%BldNodes)                  ! FMomAbvNd(K,J) = portion of the first moment of blade K about the rotor centerline (not root, like FirstMom(K)) associated with everything above node J (including tip brake masses).
   REAL(ReKi)                   :: KBECent   (p%NumBl,1,1)                         ! Centrifugal-term of generalized edgewise stiffness of the blades.
   REAL(ReKi)                   :: KBFCent   (p%NumBl,2,2)                         ! Centrifugal-term of generalized flapwise stiffness of the blades.
   REAL(ReKi)                   :: KTFAGrav  (2,2)                                 ! Gravitational-term of generalized fore-aft stiffness of the tower.
   REAL(ReKi)                   :: KTSSGrav  (2,2)                                 ! Gravitational-term of generalized side-to-side stiffness of the tower.
   REAL(ReKi)                   :: MBE       (p%NumBl,1,1)                         ! Generalized edgewise mass of the blades.
   REAL(ReKi)                   :: MBF       (p%NumBl,2,2)                         ! Generalized flapwise mass of the blades.
   REAL(ReKi)                   :: MTFA      (2,2)                                 ! Generalized fore-aft mass of the tower.
   REAL(ReKi)                   :: MTSS      (2,2)                                 ! Generalized side-to-side mass of the tower.
   REAL(ReKi)                   :: Shape                                           ! Temporary result holding a value from the SHP function
   REAL(ReKi)                   :: Shape1                                          ! Temporary result holding a value from the SHP function
   REAL(ReKi)                   :: Shape2                                          ! Temporary result holding a value from the SHP function
   REAL(ReKi)                   :: TMssAbvNd (p%TwrNodes)                          ! Portion of the tower mass associated with everything above node J (including tower-top effects)
   REAL(ReKi)                   :: TwstdSF   (2,3,0:1)                             ! Temperory result holding the current addition to the TwistedSF() array.
   REAL(ReKi)                   :: TwstdSFOld(2,3,0:1)                             ! Previous TwstdSF (i.e., TwstdSF from the previous node)

   INTEGER(IntKi)               :: I                                               ! Generic index.
   INTEGER(IntKi)               :: J                                               ! Loops through nodes / elements.
   INTEGER(IntKi)               :: K                                               ! Loops through blades.
   INTEGER(IntKi)               :: L                                               ! Generic index


   ErrStat = ErrID_None
   ErrMsg  = ''

   !...............................................................................................................................
   ! Calculate the distances from point S on a blade to the aerodynamic center in the j1 and j2 directions:
   !...............................................................................................................................

   DO K = 1,p%NumBl          ! Loop through the blades

      DO J = 1,p%BldNodes    ! Loop through the blade nodes / elements

         TmpDist           = ( 0.25 - p%PitchAxis(K,J) )*p%Chord(J)   ! Distance along the chordline from point S (25% chord) to the aerodynamic center of the blade element J--positive towards the trailing edge.
         TmpDistj1         = TmpDist*p%SAeroTwst(J)                   ! Distance along the j1-axis   from point S (25% chord) to the aerodynamic center of the blade element J
         TmpDistj2         = TmpDist*p%CAeroTwst(J)                   ! Distance along the j2-axis   from point S (25% chord) to the aerodynamic center of the blade element J
         p%rSAerCenn1(K,J) = TmpDistj1*p%CThetaS(K,J) - TmpDistj2*p%SThetaS(K,J)
         p%rSAerCenn2(K,J) = TmpDistj1*p%SThetaS(K,J) + TmpDistj2*p%CThetaS(K,J)

      ENDDO ! J - Blade nodes / elements

   ENDDO    ! K - Blades


   !...............................................................................................................................
   ! Calculate the structure that furls with the rotor inertia term:
   !...............................................................................................................................

   p%RrfaIner  = InputFileData%RFrlIner - p%RFrlMass*(      (p%rVDxn**2    )*( 1.0 - p%CRFrlSkw2*p%CRFrlTlt2 ) &
                                     +    (p%rVDzn**2    )*                    p%CRFrlTlt2   &
                                     +    (p%rVDyn**2    )*( 1.0 - p%SRFrlSkw2*p%CRFrlTlt2 ) &
                                     - 2.0*p%rVDxn*p%rVDzn*        p%CRFrlSkew*p%CSRFrlTlt   &
                                     - 2.0*p%rVDxn*p%rVDyn*        p%CSRFrlSkw*p%CRFrlTlt2   &
                                     - 2.0*p%rVDzn*p%rVDyn*        p%SRFrlSkew*p%CSRFrlTlt     )
   IF ( p%RrfaIner < 0.0 )   THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' RFrlIner must not be less than RFrlMass*( perpendicular distance between rotor-furl'// &
               ' axis and CM of the structure that furls with the rotor [not including rotor] )^2.'
      RETURN
   END IF

   !...............................................................................................................................
   ! Calculate the tail boom inertia term:
   !...............................................................................................................................

   p%AtfaIner  = p%TFrlIner - p%BoomMass*(   p%rWIxn*p%rWIxn*( 1.0 - p%CTFrlSkw2*p%CTFrlTlt2 ) &
                                       +     p%rWIzn*p%rWIzn*                    p%CTFrlTlt2   &
                                       +     p%rWIyn*p%rWIyn*( 1.0 - p%STFrlSkw2*p%CTFrlTlt2 ) &
                                       - 2.0*p%rWIxn*p%rWIzn*        p%CTFrlSkew*p%CSTFrlTlt   &
                                       - 2.0*p%rWIxn*p%rWIyn*        p%CSTFrlSkw*p%CTFrlTlt2   &
                                       - 2.0*p%rWIzn*p%rWIyn*        p%STFrlSkew*p%CSTFrlTlt     )
   IF ( p%AtfaIner < 0.0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' TFrlIner must not be less than BoomMass*( perpendicular distance between tail-furl'// &
                                        ' axis and tail boom CM )^2.'
      RETURN
   ENDIF

   !...............................................................................................................................
   ! Calculate the nacelle inertia terms:
   !...............................................................................................................................

   p%Nacd2Iner = InputFileData%NacYIner - p%NacMass*( p%NacCMxn**2 + p%NacCMyn**2 ) ! Nacelle inertia about the d2-axis
   IF ( p%Nacd2Iner < 0.0 ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg = ' NacYIner must not be less than NacMass*( NacCMxn^2 + NacCMyn^2 ).'
      RETURN
   END IF

      ! Calculate hub inertia about its centerline passing through its c.g..
      !   This calculation assumes that the hub for a 2-blader is essentially
      !   a uniform cylinder whose centerline is transverse through the cylinder
      !   passing through its c.g..  That is, for a 2-blader, Hubg1Iner =
      !   Hubg2Iner is the inertia of the hub about both the g1- and g2- axes.  For
      !   3-bladers, Hubg1Iner is simply equal to HubIner and Hubg2Iner is zero.
      ! Also, Initialize RotMass and RotIner to associated hub properties:

   IF ( p%NumBl == 2 )  THEN ! 2-blader
      p%Hubg1Iner = ( InputFileData%HubIner - p%HubMass*( ( p%UndSling - p%HubCM )**2 ) )/( p%CosDel3**2 )
      p%Hubg2Iner = p%Hubg1Iner
      IF ( p%Hubg1Iner < 0.0 ) THEN
         ErrStat = ErrID_Fatal
         ErrMsg = ' HubIner must not be less than HubMass*( UndSling - HubCM )^2 for 2-blader.'
         RETURN
      END IF
   ELSE                    ! 3-blader
      p%Hubg1Iner = InputFileData%HubIner
      p%Hubg2Iner = 0.0
   ENDIF

   p%RotMass   = p%HubMass
   p%RotIner   = p%Hubg1Iner


   !...............................................................................................................................

      ! Initialize several variables to 0.0:

   p%KBF     = 0.0
   p%KBE     = 0.0
   KBFCent   = 0.0
   KBECent   = 0.0

   p%TwrMass = 0.0
   p%KTFA    = 0.0
   p%KTSS    = 0.0
   KTFAGrav  = 0.0
   KTSSGrav  = 0.0



   DO K = 1,p%NumBl          ! Loop through the blades


      ! Initialize BldMass(), FirstMom(), and SecondMom() using TipMass() effects:

      p%BldMass  (K) = p%TipMass(K)
      p%FirstMom (K) = p%TipMass(K)*p%BldFlexL
      p%SecondMom(K) = p%TipMass(K)*p%BldFlexL*p%BldFlexL


      DO J = p%BldNodes,1,-1 ! Loop through the blade nodes / elements in reverse


      ! Calculate the mass of the current element

         p%BElmntMass(J,K) = p%MassB(K,J)*p%DRNodes(J)                        ! Mass of blade element J


      ! Integrate to find some blade properties which will be output in .fsm

         p%BldMass  (K) = p%BldMass  (K) + p%BElmntMass(J,K)
         p%FirstMom (K) = p%FirstMom (K) + p%BElmntMass(J,K)*p%RNodes(J)
         p%SecondMom(K) = p%SecondMom(K) + p%BElmntMass(J,K)*p%RNodes(J)*p%RNodes(J)


      ! Integrate to find FMomAbvNd:

         FMomAbvNd   (K,J) = ( 0.5*p%BElmntMass(J,K) )*( p%HubRad + p%RNodes(J  ) + 0.5*p%DRNodes(J  ) )

         IF ( J == p%BldNodes )  THEN ! Outermost blade element
      ! Add the TipMass() effects:

            FMomAbvNd(K,J) = FmomAbvNd(K,J) + p%TipMass(K)*p%TipRad
         ELSE                       ! All other blade elements
      ! Add to FMomAbvNd(K,J) the effects from the (not yet used) portion of element J+1

            FMomAbvNd(K,J) = FMomAbvNd(K,J) + FMomAbvNd(K,J+1) &
                           + ( 0.5*p%BElmntMass(J+1,K) )*( p%HubRad + p%RNodes(J+1) - 0.5*p%DRNodes(J+1) )
         ENDIF


      ENDDO ! J - Blade nodes / elements in reverse

      IF (.NOT. p%BD4Blades) THEN
         ! Calculate BldCG() using FirstMom() and BldMass(); and calculate
         !   RotMass and RotIner:

         p%BldCG    (K) = p%FirstMom (K) / p%BldMass    (K)
         p%RotMass      = p%RotMass      + p%BldMass    (K)
         p%RotIner      = p%RotIner      + ( p%SecondMom(K) + p%BldMass  (K)*p%HubRad*( 2.0*p%BldCG(K) + p%HubRad ) )*( p%CosPreC(K)**2 )
      END IF

   ENDDO ! K - Blades



   DO K = 1,p%NumBl          ! Loop through the blades


      ! Initialize the generalized blade masses using tip mass effects:

      MBF(K,1,1) = p%TipMass(K)
      MBF(K,2,2) = p%TipMass(K)
      MBE(K,1,1) = p%TipMass(K)


      DO J = 1,p%BldNodes    ! Loop through the blade nodes / elements


      ! Integrate to find the generalized mass of the blade (including tip mass effects).
      !   Ignore the cross-correlation terms of MBF (i.e. MBF(i,j) where i /= j) since
      !   these terms will never be used.

         Shape1 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl1Sh(:,K), 0, ErrStat, ErrMsg )
         Shape2 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl2Sh(:,K), 0, ErrStat, ErrMsg )
         MBF    (K,1,1) = MBF    (K,1,1) + p%BElmntMass(J,K)*Shape1*Shape1
         MBF    (K,2,2) = MBF    (K,2,2) + p%BElmntMass(J,K)*Shape2*Shape2

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldEdgSh(:,K), 0, ErrStat, ErrMsg )
         MBE    (K,1,1) = MBE    (K,1,1) + p%BElmntMass(J,K)*Shape *Shape


      ! Integrate to find the generalized stiffness of the blade (not including centrifugal
      !    effects).

         ElmntStff      = p%StiffBF(K,J)*p%DRNodes(J)                       ! Flapwise stiffness of blade element J
         Shape1 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl1Sh(:,K), 2, ErrStat, ErrMsg )
         Shape2 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl2Sh(:,K), 2, ErrStat, ErrMsg )
         p%KBF    (K,1,1) = p%KBF    (K,1,1) + ElmntStff*Shape1*Shape1
         p%KBF    (K,1,2) = p%KBF    (K,1,2) + ElmntStff*Shape1*Shape2
         p%KBF    (K,2,1) = p%KBF    (K,2,1) + ElmntStff*Shape2*Shape1
         p%KBF    (K,2,2) = p%KBF    (K,2,2) + ElmntStff*Shape2*Shape2

         ElmntStff      = p%StiffBE(K,J)*p%DRNodes(J)                       ! Edgewise stiffness of blade element J
         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldEdgSh(:,K), 2, ErrStat, ErrMsg )
         p%KBE    (K,1,1) = p%KBE    (K,1,1) + ElmntStff*Shape *Shape


      ! Integrate to find the centrifugal-term of the generalized flapwise and edgewise
      !   stiffness of the blades.  Ignore the cross-correlation terms of KBFCent (i.e.
      !   KBFCent(i,j) where i /= j) since these terms will never be used.

         ElmntStff      = FMomAbvNd(K,J)*p%DRNodes(J)*p%RotSpeed**2   ! Centrifugal stiffness of blade element J

         Shape1 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl1Sh(:,K), 1, ErrStat, ErrMsg )
         Shape2 = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl2Sh(:,K), 1, ErrStat, ErrMsg )
         KBFCent(K,1,1) = KBFCent(K,1,1) + ElmntStff*Shape1*Shape1
         KBFCent(K,2,2) = KBFCent(K,2,2) + ElmntStff*Shape2*Shape2

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldEdgSh(:,K), 1, ErrStat, ErrMsg )
         KBECent(K,1,1) = KBECent(K,1,1) + ElmntStff*Shape *Shape


      ! Calculate the 2nd derivatives of the twisted shape functions:

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl1Sh(:,K), 2, ErrStat, ErrMsg )
         p%TwistedSF(K,1,1,J,2) =  Shape*p%CThetaS(K,J)                  ! 2nd deriv. of Phi1(J) for blade K
         p%TwistedSF(K,2,1,J,2) = -Shape*p%SThetaS(K,J)                  ! 2nd deriv. of Psi1(J) for blade K

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldFl2Sh(:,K), 2, ErrStat, ErrMsg )
         p%TwistedSF(K,1,2,J,2) =  Shape*p%CThetaS(K,J)                  ! 2nd deriv. of Phi2(J) for blade K
         p%TwistedSF(K,2,2,J,2) = -Shape*p%SThetaS(K,J)                  ! 2nd deriv. of Psi2(J) for blade K

         Shape  = SHP( p%RNodesNorm(J), p%BldFlexL, p%BldEdgSh(:,K), 2, ErrStat, ErrMsg )
         p%TwistedSF(K,1,3,J,2) =  Shape*p%SThetaS(K,J)                  ! 2nd deriv. of Phi3(J) for blade K
         p%TwistedSF(K,2,3,J,2) =  Shape*p%CThetaS(K,J)                  ! 2nd deriv. of Psi3(J) for blade K


      ! Integrate to find the 1st derivatives of the twisted shape functions:

         DO I = 1,2     ! Loop through Phi and Psi
            DO L = 1,3  ! Loop through all blade DOFs
               TwstdSF     (  I,L,  1) = p%TwistedSF(K,I,L,J,2)*0.5*p%DRNodes(J)
               p%TwistedSF   (K,I,L,J,1) = TwstdSF   ( I,L,  1)
            ENDDO       ! L - All blade DOFs
         ENDDO          ! I - Phi and Psi

         IF ( J /= 1 )  THEN  ! All but the innermost blade element
      ! Add the effects from the (not yet used) portion of element J-1

            DO I = 1,2     ! Loop through Phi and Psi
               DO L = 1,3  ! Loop through all blade DOFs
                  p%TwistedSF(K,I,L,J,1) = p%TwistedSF(K,I,L,J,1) + p%TwistedSF(K,I,L,J-1,1) &
                                       + TwstdSFOld( I,L,  1)
               ENDDO       ! L - All blade DOFs
            ENDDO          ! I - Phi and Psi
         ENDIF


      ! Integrate to find the twisted shape functions themselves (i.e., their zeroeth derivative):

         DO I = 1,2     ! Loop through Phi and Psi
            DO L = 1,3  ! Loop through all blade DOFs
               TwstdSF     (  I,L,  0) = p%TwistedSF(K,I,L,J,1)*0.5*p%DRNodes(J)
               p%TwistedSF   (K,I,L,J,0) = TwstdSF   ( I,L,  0)
            ENDDO       ! L - All blade DOFs
         ENDDO          ! I - Phi and Psi

         IF ( J /= 1 )  THEN  ! All but the innermost blade element
      ! Add the effects from the (not yet used) portion of element J-1

            DO I = 1,2     ! Loop through Phi and Psi
               DO L = 1,3  ! Loop through all blade DOFs
                  p%TwistedSF(K,I,L,J,0) = p%TwistedSF(K,I,L,J,0) + p%TwistedSF(K,I,L,J-1,0) &
                                       + TwstdSFOld( I,L,  0)
               ENDDO       ! L - All blade DOFs
            ENDDO          ! I - Phi and Psi
         ENDIF


      ! Integrate to find the blade axial reduction shape functions:

         DO I = 1,3     ! Loop through all blade DOFs
            DO L = 1,3  ! Loop through all blade DOFs
               AxRdBld    (  I,L  ) = 0.5*p%DRNodes(J)*(                          &
                                      p%TwistedSF(K,1,I,J,1)*p%TwistedSF(K,1,L,J,1) &
                                    + p%TwistedSF(K,2,I,J,1)*p%TwistedSF(K,2,L,J,1) )
               p%AxRedBld   (K,I,L,J) = AxRdBld(I,L)
            ENDDO       ! L - All blade DOFs
         ENDDO          ! I - All blade DOFs

         IF ( J /= 1 )  THEN  ! All but the innermost blade element
      ! Add the effects from the (not yet used) portion of element J-1

            DO I = 1,3     ! Loop through all blade DOFs
               DO L = 1,3  ! Loop through all blade DOFs
                  p%AxRedBld(K,I,L,J) = p%AxRedBld(K,I,L,J) + p%AxRedBld(K,I,L,J-1)   &
                                    + AxRdBldOld(I,L)
               ENDDO       ! L - All blade DOFs
            ENDDO          ! I - All blade DOFs
         ENDIF


      ! Store the TwstdSF and AxRdBld terms of the current element (these will be used for the next element)

         TwstdSFOld = TwstdSF
         AxRdBldOld = AxRdBld


      ENDDO ! J - Blade nodes / elements




   IF (p%BD4Blades) THEN

      !p%KBF     ( K,:,:    ) = 0.0_ReKi
      
         ! the 1st and zeroeth derivatives of the twisted shape functions at the blade root:
      p%TwistedSF(K,:,:,:,1) = 0.0_ReKi
      p%TwistedSF(K,:,:,:,0) = 0.0_ReKi 
      p%AxRedBld( K,:,:,:  ) = 0.0_ReKi
   ELSE
            
      ! Apply the flapwise modal stiffness tuners of the blades to KBF():

      DO I = 1,2     ! Loop through flap DOFs
         DO L = 1,2  ! Loop through flap DOFs
            p%KBF(K,I,L) = SQRT( p%FStTunr(K,I)*p%FStTunr(K,L) )*p%KBF(K,I,L)
         ENDDO       ! L - Flap DOFs
      ENDDO          ! I - Flap DOFs
      
      ! Calculate the blade natural frequencies:
      
      DO I = 1,NumBF     ! Loop through flap DOFs
         p%FreqBF(K,I,1) = Inv2Pi*SQRT(   p%KBF(K,I,I)                   /( MBF(K,I,I) - p%TipMass(K) ) )   ! Natural blade I-flap frequency w/o centrifugal stiffening nor     tip mass effects
         p%FreqBF(K,I,2) = Inv2Pi*SQRT(   p%KBF(K,I,I)                   /  MBF(K,I,I)                )     ! Natural blade I-flap frequency w/o centrifugal stiffening, but w/ tip mass effects
         p%FreqBF(K,I,3) = Inv2Pi*SQRT( ( p%KBF(K,I,I) + KBFCent(K,I,I) )/  MBF(K,I,I)                )     ! Natural blade I-flap frequency w/  centrifugal stiffening and     tip mass effects
      ENDDO          ! I - Flap DOFs

      p%FreqBE   (K,1,1) = Inv2Pi*SQRT(   p%KBE(K,1,1)                   /( MBE(K,1,1) - p%TipMass(K) ) )   ! Natural blade 1-edge frequency w/o centrifugal stiffening nor      tip mass effects
      p%FreqBE   (K,1,2) = Inv2Pi*SQRT(   p%KBE(K,1,1)                   /  MBE(K,1,1)                )     ! Natural Blade 1-edge frequency w/o  centrifugal stiffening, but w/ tip mass effects
      p%FreqBE   (K,1,3) = Inv2Pi*SQRT( ( p%KBE(K,1,1) + KBECent(K,1,1) )/  MBE(K,1,1)                )     ! Natural Blade 1-edge frequency w/  centrifugal stiffening and      tip mass effects


      ! Calculate the generalized damping of the blades:

      DO I = 1,NumBF     ! Loop through flap DOFs
         DO L = 1,NumBF  ! Loop through flap DOFs
            p%CBF(K,I,L) = ( 0.01*p%BldFDamp(K,L) )*p%KBF(K,I,L)/( Pi*p%FreqBF(K,L,1) )
         ENDDO       ! L - Flap DOFs
      ENDDO          ! I - Flap DOFs

      p%CBE      (K,1,1) = ( 0.01*p%BldEDamp(K,1) )*p%KBE(K,1,1)/( Pi*p%FreqBE(K,1,1) )


      ! Calculate the 2nd derivatives of the twisted shape functions at the blade root:

      Shape  = SHP( 0.0_ReKi, p%BldFlexL, p%BldFl1Sh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,1,0,2) =  Shape*p%CThetaS(K,0)        ! 2nd deriv. of Phi1(0) for blade K
      p%TwistedSF(K,2,1,0,2) = -Shape*p%SThetaS(K,0)        ! 2nd deriv. of Psi1(0) for blade K

      Shape  = SHP( 0.0_ReKi, p%BldFlexL, p%BldFl2Sh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,2,0,2) =  Shape*p%CThetaS(K,0)        ! 2nd deriv. of Phi2(0) for blade K
      p%TwistedSF(K,2,2,0,2) = -Shape*p%SThetaS(K,0)        ! 2nd deriv. of Psi2(0) for blade K

      Shape  = SHP( 0.0_ReKi, p%BldFlexL, p%BldEdgSh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,3,0,2) =  Shape*p%SThetaS(K,0)        ! 2nd deriv. of Phi3(0) for blade K
      p%TwistedSF(K,2,3,0,2) =  Shape*p%CThetaS(K,0)        ! 2nd deriv. of Psi3(0) for blade K
      
      
      ! Calculate the 2nd derivatives of the twisted shape functions at the tip:

      Shape  = SHP( 1.0_ReKi, p%BldFlexL, p%BldFl1Sh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,1,p%TipNode,2) =  Shape*p%CThetaS(K,p%TipNode)        ! 2nd deriv. of Phi1(p%TipNode) for blade K
      p%TwistedSF(K,2,1,p%TipNode,2) = -Shape*p%SThetaS(K,p%TipNode)        ! 2nd deriv. of Psi1(p%TipNode) for blade K

      Shape  = SHP( 1.0_ReKi, p%BldFlexL, p%BldFl2Sh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,2,p%TipNode,2) =  Shape*p%CThetaS(K,p%TipNode)        ! 2nd deriv. of Phi2(p%TipNode) for blade K
      p%TwistedSF(K,2,2,p%TipNode,2) = -Shape*p%SThetaS(K,p%TipNode)        ! 2nd deriv. of Psi2(p%TipNode) for blade K

      Shape  = SHP( 1.0_ReKi, p%BldFlexL, p%BldEdgSh(:,K), 2, ErrStat, ErrMsg )
      p%TwistedSF(K,1,3,p%TipNode,2) =  Shape*p%SThetaS(K,p%TipNode)        ! 2nd deriv. of Phi3(p%TipNode) for blade K
      p%TwistedSF(K,2,3,p%TipNode,2) =  Shape*p%CThetaS(K,p%TipNode)        ! 2nd deriv. of Psi3(p%TipNode) for blade K


      ! Integrate to find the 1st and zeroeth derivatives of the twisted shape functions
      !   at the tip:

      DO I = 1,2     ! Loop through Phi and Psi
         DO L = 1,3  ! Loop through all blade DOFs
            p%TwistedSF(K,I,L,p%TipNode,1) = p%TwistedSF(K,I,L,p%BldNodes,1) + TwstdSFOld(I,L,1)
            p%TwistedSF(K,I,L,p%TipNode,0) = p%TwistedSF(K,I,L,p%BldNodes,0) + TwstdSFOld(I,L,0)
         ENDDO       ! L - All blade DOFs
      ENDDO          ! I - Phi and Psi

         ! the 1st and zeroeth derivatives of the twisted shape functions at the blade root:
      p%TwistedSF(K,:,:,0,1) = 0.0_ReKi
      p%TwistedSF(K,:,:,0,0) = 0.0_ReKi 
      p%AxRedBld( K,:,:,0  ) = 0.0_ReKi

      ! Integrate to find the blade axial reduction shape functions at the tip:

      DO I = 1,3     ! Loop through all blade DOFs
         DO L = 1,3  ! Loop through all blade DOFs
            p%AxRedBld(K,I,L,p%TipNode) = p%AxRedBld(K,I,L,p%BldNodes) + AxRdBldOld(I,L)
         ENDDO       ! L - All blade DOFs
      ENDDO          ! I - All blade DOFs
   END IF ! p%BD4Blades


   ENDDO ! K - Blades



      ! Calculate the tower-top mass:

   p%TwrTpMass = p%RotMass + p%RFrlMass + p%BoomMass + p%TFinMass + p%NacMass + p%YawBrMass


   DO J = p%TwrNodes,1,-1 ! Loop through the tower nodes / elements in reverse


      ! Calculate the mass of the current element

      p%TElmntMass(J)    = p%MassT(J)*p%DHNodes(J)     ! Mass of tower element J


      ! Integrate to find the tower mass which will be output in .fsm

      p%TwrMass      = p%TwrMass + p%TElmntMass(J)


      ! Integrate to find TMssAbvNd:

      TMssAbvNd   (J) = 0.5*p%TElmntMass(J)

      IF ( J == p%TwrNodes )  THEN ! Uppermost tower element
      ! Add the TwrTpMass effects:

         TMssAbvNd(J) = TMssAbvNd(J) + p%TwrTpMass
      ELSE                       ! All other tower elements
      ! Add to TMssAbvNd(J) the effects from the (not yet used) portion of element J+1

         TMssAbvNd(J) = 0.5*p%TElmntMass(J+1) + TMssAbvNd(J) + TMssAbvNd(J+1)
      ENDIF


   ENDDO ! J - Tower nodes / elements in reverse



      ! Initialize the generalized tower masses using tower-top mass effects:

   DO I = 1,2  ! Loop through all tower modes in a single direction
      MTFA(I,I) = p%TwrTpMass
      MTSS(I,I) = p%TwrTpMass
   ENDDO       ! I - All tower modes in a single direction

      ! set values for tower base (note that we haven't corrctly defined the values for (:,0,2) in the arrays below):
   p%TwrFASF(   :,0,0:1) = 0.0_ReKi
   p%TwrSSSF(   :,0,0:1) = 0.0_ReKi
   p%AxRedTFA(:,:,0)     = 0.0_ReKi
   p%AxRedTSS(:,:,0)     = 0.0_ReKi
   
   DO J = 1,p%TwrNodes    ! Loop through the tower nodes / elements


      ! Calculate the tower shape functions (all derivatives):

      p%TwrFASF(1,J,2) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM1Sh(:), 2, ErrStat, ErrMsg )
      p%TwrFASF(2,J,2) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM2Sh(:), 2, ErrStat, ErrMsg )
      p%TwrFASF(1,J,1) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM1Sh(:), 1, ErrStat, ErrMsg )
      p%TwrFASF(2,J,1) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM2Sh(:), 1, ErrStat, ErrMsg )
      p%TwrFASF(1,J,0) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM1Sh(:), 0, ErrStat, ErrMsg )
      p%TwrFASF(2,J,0) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwFAM2Sh(:), 0, ErrStat, ErrMsg )

      p%TwrSSSF(1,J,2) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM1Sh(:), 2, ErrStat, ErrMsg )
      p%TwrSSSF(2,J,2) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM2Sh(:), 2, ErrStat, ErrMsg )
      p%TwrSSSF(1,J,1) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM1Sh(:), 1, ErrStat, ErrMsg )
      p%TwrSSSF(2,J,1) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM2Sh(:), 1, ErrStat, ErrMsg )
      p%TwrSSSF(1,J,0) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM1Sh(:), 0, ErrStat, ErrMsg )
      p%TwrSSSF(2,J,0) = SHP( p%HNodesNorm(J), p%TwrFlexL, InputFileData%TwSSM2Sh(:), 0, ErrStat, ErrMsg )


      ! Integrate to find the generalized mass of the tower (including tower-top mass effects).
      !   Ignore the cross-correlation terms of MTFA (i.e. MTFA(i,j) where i /= j) and MTSS
      !   since these terms will never be used.


      DO I = 1,2     ! Loop through all tower DOFs in one direction
         MTFA  (I,I) = MTFA  (I,I) + p%TElmntMass(J)*p%TwrFASF(I,J,0)**2
         MTSS  (I,I) = MTSS  (I,I) + p%TElmntMass(J)*p%TwrSSSF(I,J,0)**2
      ENDDO          ! I - through all tower DOFs in one direction


      ! Integrate to find the generalized stiffness of the tower (not including gravitational
      !    effects).

      ElStffFA       = p%StiffTFA(J)*p%DHNodes(J)                        ! Fore-aft stiffness of tower element J
      ElStffSS       = p%StiffTSS(J)*p%DHNodes(J)                        ! Side-to-side stiffness of tower element J

      DO I = 1,2     ! Loop through all tower DOFs in one direction
         DO L = 1,2  ! Loop through all tower DOFs in one direction
            p%KTFA (I,L) = p%KTFA    (I,L) + ElStffFA *p%TwrFASF(I,J,2)*p%TwrFASF(L,J,2)
            p%KTSS (I,L) = p%KTSS    (I,L) + ElStffSS *p%TwrSSSF(I,J,2)*p%TwrSSSF(L,J,2)
         ENDDO       ! L - All tower DOFs in one direction
      ENDDO          ! I - through all tower DOFs in one direction


      ! Integrate to find the gravitational-term of the generalized stiffness of the tower.
      !   Ignore the cross-correlation terms of KTFAGrav (i.e. KTFAGrav(i,j) where i /= j)
      !   and KTSSGrav since these terms will never be used.

      ElmntStff      = -TMssAbvNd(J)*p%DHNodes(J)*p%Gravity              ! Gravitational stiffness of tower element J

      DO I = 1,2     ! Loop through all tower DOFs in one direction
         KTFAGrav(I,I) = KTFAGrav(I,I) + ElmntStff*p%TwrFASF(I,J,1)**2
         KTSSGrav(I,I) = KTSSGrav(I,I) + ElmntStff*p%TwrSSSF(I,J,1)**2
      ENDDO


      ! Integrate to find the tower axial reduction shape functions:

      DO I = 1,2     ! Loop through all tower DOFs in one direction
         DO L = 1,2  ! Loop through all tower DOFs in one direction
            AxRdTFA (I,L) = 0.5*p%DHNodes(J)*p%TwrFASF(I,J,1)*p%TwrFASF(L,J,1)
            AxRdTSS (I,L) = 0.5*p%DHNodes(J)*p%TwrSSSF(I,J,1)*p%TwrSSSF(L,J,1)

            p%AxRedTFA(I,L,J) = AxRdTFA(I,L)
            p%AxRedTSS(I,L,J) = AxRdTSS(I,L)
         ENDDO       ! L - All tower DOFs in one direction
      ENDDO

      IF ( J /= 1 )  THEN  ! All but the lowermost tower element
      ! Add the effects from the (not yet used) portion of element J-1

         DO I = 1,2     ! Loop through all tower DOFs in one direction
            DO L = 1,2  ! Loop through all tower DOFs in one direction
               p%AxRedTFA(I,L,J) = p%AxRedTFA(I,L,J) + p%AxRedTFA(I,L,J-1)+ AxRdTFAOld(I,L)
               p%AxRedTSS(I,L,J) = p%AxRedTSS(I,L,J) + p%AxRedTSS(I,L,J-1)+ AxRdTSSOld(I,L)
            ENDDO       ! L - All tower DOFs in one direction
         ENDDO
      ENDIF


      ! Store the AxRdTFA and AxRdTSS terms of the current element (these will be used for the next element)

      AxRdTFAOld = AxRdTFA
      AxRdTSSOld = AxRdTSS


   ENDDO ! J - Tower nodes / elements


   ! Apply the modal stiffness tuners of the tower to KTFA() and KTSS():

   DO I = 1,2     ! Loop through all tower DOFs in one direction
      DO L = 1,2  ! Loop through all tower DOFs in one direction
         p%KTFA(I,L) = SQRT( InputFileData%FAStTunr(I)*InputFileData%FAStTunr(L) )*p%KTFA(I,L)

         p%KTSS(I,L) = SQRT( InputFileData%SSStTunr(I)*InputFileData%SSStTunr(L) )*p%KTSS(I,L)
      ENDDO       ! L - All tower DOFs in one direction
   ENDDO          ! I - through all tower DOFs in one direction


      ! Calculate the tower natural frequencies:

   DO I = 1,2     ! Loop through all tower DOFs in one direction
      p%FreqTFA(I,1) = Inv2Pi*SQRT(   p%KTFA(I,I)                  /( MTFA(I,I) - p%TwrTpMass ) )  ! Natural tower I-fore-aft frequency w/o gravitational destiffening nor tower-top mass effects
      p%FreqTFA(I,2) = Inv2Pi*SQRT( ( p%KTFA(I,I) + KTFAGrav(I,I) )/  MTFA(I,I)               )  ! Natural tower I-fore-aft frequency w/  gravitational destiffening and tower-top mass effects
      p%FreqTSS(I,1) = Inv2Pi*SQRT(   p%KTSS(I,I)                  /( MTSS(I,I) - p%TwrTpMass ) )  ! Natural tower I-side-to-side frequency w/o gravitational destiffening nor tower-top mass effects
      p%FreqTSS(I,2) = Inv2Pi*SQRT( ( p%KTSS(I,I) + KTSSGrav(I,I) )/  MTSS(I,I)               )  ! Natural tower I-side-to-side frequency w/  gravitational destiffening and tower-top mass effects
   ENDDO          ! I - All tower DOFs in one direction


      ! Calculate the generalized damping of the tower:

   DO I = 1,2     ! Loop through all tower DOFs in one direction
      DO L = 1,2  ! Loop through all tower DOFs in one direction
         p%CTFA(I,L) = ( 0.01*InputFileData%TwrFADmp(L) )*p%KTFA(I,L)/( Pi*p%FreqTFA(L,1) )

         p%CTSS(I,L) = ( 0.01*InputFileData%TwrSSDmp(L) )*p%KTSS(I,L)/( Pi*p%FreqTSS(L,1) )
      ENDDO       ! L - All tower DOFs in one direction
   ENDDO          ! I - All tower DOFs in one direction


      ! Calculate the tower shape functions (all derivatives) at the tower-top:

   p%TwrFASF(1,p%TTopNode,2) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM1Sh(:), 2, ErrStat, ErrMsg )
   p%TwrFASF(2,p%TTopNode,2) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM2Sh(:), 2, ErrStat, ErrMsg )
   p%TwrFASF(1,p%TTopNode,1) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM1Sh(:), 1, ErrStat, ErrMsg )
   p%TwrFASF(2,p%TTopNode,1) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM2Sh(:), 1, ErrStat, ErrMsg )
   p%TwrFASF(1,p%TTopNode,0) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM1Sh(:), 0, ErrStat, ErrMsg )
   p%TwrFASF(2,p%TTopNode,0) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwFAM2Sh(:), 0, ErrStat, ErrMsg )

   p%TwrSSSF(1,p%TTopNode,2) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM1Sh(:), 2, ErrStat, ErrMsg )
   p%TwrSSSF(2,p%TTopNode,2) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM2Sh(:), 2, ErrStat, ErrMsg )
   p%TwrSSSF(1,p%TTopNode,1) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM1Sh(:), 1, ErrStat, ErrMsg )
   p%TwrSSSF(2,p%TTopNode,1) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM2Sh(:), 1, ErrStat, ErrMsg )
   p%TwrSSSF(1,p%TTopNode,0) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM1Sh(:), 0, ErrStat, ErrMsg )
   p%TwrSSSF(2,p%TTopNode,0) = SHP( 1.0_ReKi, p%TwrFlexL, InputFileData%TwSSM2Sh(:), 0, ErrStat, ErrMsg )


      ! Integrate to find the tower axial reduction shape functions at the tower-top:

   DO I = 1,2     ! Loop through all tower DOFs in one direction
      DO L = 1,2  ! Loop through all tower DOFs in one direction
         p%AxRedTFA(I,L,p%TTopNode) = p%AxRedTFA(I,L,p%TwrNodes)+ AxRdTFAOld(I,L)
         p%AxRedTSS(I,L,p%TTopNode) = p%AxRedTSS(I,L,p%TwrNodes)+ AxRdTSSOld(I,L)
      ENDDO       ! L - All tower DOFs in one direction
   ENDDO


      ! Calculate the turbine mass:

   p%TurbMass  = p%TwrTpMass + p%TwrMass


   RETURN
END SUBROUTINE Coeff
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE InitBlDefl ( p, InputFileData, InitQF1, InitQF2, InitQE1, ErrStat, ErrMsg )
! This routine calculates the initial blade deflections.
! Base the intial values of the blade DOFs, INITQF1, INITQF2, and
!   INITQE1, on OoPDefl and IPDefl.
! Write messages to the screen if the specified initial tip displacements
!  are incompatible with the enabled DOFs.
!..................................................................................................................................


      IMPLICIT                        NONE

      ! Passed variables:
   TYPE(ED_ParameterType),  INTENT(IN)  :: p                                       ! parameters of the structural dynamics module
   TYPE(ED_InputFile),      INTENT(IN)  :: InputFileData                           ! all the data in the ElastoDyn input file

   REAL(ReKi),              INTENT(OUT) :: InitQE1(p%NumBl)                        ! Initial edge deflection (output).
   REAL(ReKi),              INTENT(OUT) :: InitQF1(p%NumBl)                        ! Initial flap deflection for mode 1 (output).
   REAL(ReKi),              INTENT(OUT) :: InitQF2(p%NumBl)                        ! Initial flap deflection for mode 2 (output).

   INTEGER(IntKi),          INTENT(OUT) :: ErrStat                                 ! Error status
   CHARACTER(1024),         INTENT(OUT) :: ErrMsg                                  ! Error message when ErrStat =/ ErrID_None


      ! Local variables:
   REAL(ReKi)                   :: A(2,3)                                          ! Augmented matrix for solution of initial deflections.
   REAL(ReKi)                   :: CosPitch                                        ! Cosine of the pitch for this blade.
   REAL(ReKi)                   :: Det                                             ! Determinate of right-hand side of A.
   REAL(ReKi)                   :: SinPitch                                        ! Sine of the pitch for this blade.
   REAL(ReKi)                   :: TotResid                                        ! Generator torque.

   INTEGER(IntKi)               :: K                                               ! Blade number

      ! some warning messages
   CHARACTER(*), PARAMETER      :: Approx   = ' An approximate characterization of the specified blade deflection will be made.'
   CHARACTER(*), PARAMETER      :: BadIP    = ' Initial blade in-plane tip displacement will be ignored.'
   CHARACTER(*), PARAMETER      :: BadOoP   = ' Initial blade out-of-plane tip displacement will be ignored.'
   CHARACTER(*), PARAMETER      :: Ignore   = ' All initial blade tip displacements will be ignored.'


      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''

   InitQE1 = 0.0
   InitQF1 = 0.0
   InitQF2 = 0.0
   !bjj: replace InitQF1 and InitQF2 with an array to avoid so much duplication of logic here...

   DO K=1,p%NumBl

         ! Calculate the array of deflections(???).

      CosPitch = COS( InputFileData%BlPitch(K) )
      SinPitch = SIN( InputFileData%BlPitch(K) )

      A(1,2) =  p%TwistedSF(K,1,3,p%TipNode,0)*CosPitch + p%TwistedSF(K,2,3,p%TipNode,0)*SinPitch
      A(2,2) = -p%TwistedSF(K,1,3,p%TipNode,0)*SinPitch + p%TwistedSF(K,2,3,p%TipNode,0)*CosPitch
      A(1,3) =  InputFileData%OoPDefl
      A(2,3) =  InputFileData%IPDefl

      IF ( InputFileData%FlapDOF1 )  THEN                                ! Blade flap mode 1 is enabled

         A(1,1) =  p%TwistedSF(K,1,1,p%TipNode,0)*CosPitch + p%TwistedSF(K,2,1,p%TipNode,0)*SinPitch
         A(2,1) = -p%TwistedSF(K,1,1,p%TipNode,0)*SinPitch + p%TwistedSF(K,2,1,p%TipNode,0)*CosPitch

         DET = ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )

         IF ( .NOT. EqualRealNos( DET, 0.0_ReKi ) ) THEN                  ! Apply all flap deflection to mode 1

            InitQF1(K) = ( A(1,3)*A(2,2) - A(1,2)*A(2,3) )/DET
            InitQE1(K) = ( A(1,1)*A(2,3) - A(1,3)*A(2,1) )/DET

         ELSEIF ( .NOT. InputFileData%EdgeDOF )  THEN                     ! Blade edge mode 1 is not enabled which caused DET = 0.

            InitQE1(K) = 0.0

            IF ( .NOT. EqualRealNos( A(1,1), 0.0_ReKi ) )  THEN
               IF ( .NOT. EqualRealNos( A(2,1), 0.0_ReKi ) )  THEN        ! Find a solution of the 2 equations in 1 variable that
                                                                          !  minimizes the sum of the squares of the equation's residuals.

                  InitQF1(K) = ( A(1,1)*A(1,3) + A(2,1)*A(2,3) )/( A(1,1)**2 + A(2,1)**2 )

                  TotResid = SQRT( ( A(1,1)*InitQF1(K) - A(1,3) )**2 + ( A(2,1)*InitQF1(K) - A(2,3) )**2 )

                  IF ( .NOT. EqualRealNos( TotResid, 0.0_ReKi ) ) THEN
                     CALL CheckError( ErrID_Warn, Approx )
                  ENDIF

               ELSE !A(1,1) /= 0; A(2,1) == 0

                  InitQF1(K) = A(1,3)/A(1,1)

                  IF ( .NOT. EqualRealNos( InputFileData%IPDefl,  0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, BadIP )
                  ENDIF

               ENDIF

            ELSE ! A(1,1) == 0

               IF ( .NOT. EqualRealNos( InputFileData%OoPDefl, 0.0_ReKi ) ) THEN
                  CALL CheckError( ErrID_Warn, BadOoP )
               END IF

               IF ( .NOT. EqualRealNos( A(2,1), 0.0_ReKi ) )   THEN
                  InitQF1(K) = A(2,3)/A(2,1)
               ELSE
                  InitQF1(K) = 0.0

                  IF ( .NOT. EqualRealNos( InputFileData%IPDefl,  0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, BadIP )
                  ENDIF

               ENDIF
            ENDIF

         ELSE                                     ! It is impossible to find any "good" solution, so ignore the initial tip displacements

            InitQF1(K) = 0.0
            InitQE1(K) = 0.0

            IF ( ( InputFileData%OoPDefl /= 0.0 ) .OR. ( InputFileData%IPDefl /= 0.0 ) )  THEN
               CALL CheckError( ErrID_Warn, Ignore )
            ENDIF

         ENDIF

      ELSE                                        ! Blade flap mode 1 is not enabled.

         InitQF1(K) = 0.0

         IF ( InputFileData%FlapDOF2 )  THEN                    ! Blade flap mode 2 is enabled.

            A(1,1) =  p%TwistedSF(K,1,2,p%TipNode,0)*CosPitch + p%TwistedSF(K,2,2,p%TipNode,0)*SinPitch
            A(2,1) = -p%TwistedSF(K,1,2,p%TipNode,0)*SinPitch + p%TwistedSF(K,2,2,p%TipNode,0)*CosPitch

            DET = ( A(1,1)*A(2,2) - A(1,2)*A(2,1) )

            IF ( .NOT. EqualRealNos( DET, 0.0_ReKi ) ) THEN      ! Apply all flap deflection to mode 2
               InitQF2 = ( A(1,3)*A(2,2) - A(1,2)*A(2,3) )/DET
               InitQE1 = ( A(1,1)*A(2,3) - A(1,3)*A(2,1) )/DET

            ELSEIF ( .NOT. InputFileData%EdgeDOF )  THEN          ! Blade edge mode 1 is not enabled which caused DET = 0.

               InitQE1(K) = 0.0

               IF ( .NOT. EqualRealNos( A(1,1), 0.0_ReKi ) )  THEN
                  IF ( .NOT. EqualRealNos( A(2,1), 0.0_ReKi ) )   THEN      ! Find a solution of the 2 equations in 1 variable that
                                                                            !  minimizes the sum of the squares of the equation's residuals
                     InitQF2(K) = ( A(1,1)*A(1,3) + A(2,1)*A(2,3) )/( A(1,1)**2 + A(2,1)**2 )

                     TotResid = SQRT( ( A(1,1)*InitQF2(K) - A(1,3))**2 + ( A(2,1)*InitQF2(K) - A(2,3) )**2 )

                     IF ( .NOT. EqualRealNos( TotResid, 0.0_ReKi ) )  THEN
                        CALL CheckError( ErrID_Warn, Approx )
                     ENDIF
                  ELSE
                     InitQF2(K) = A(1,3)/A(1,1)

                     IF ( .NOT. EqualRealNos( InputFileData%IPDefl, 0.0_ReKi ) ) THEN
                        CALL CheckError( ErrID_Warn, BadIP )
                     ENDIF
                  ENDIF
               ELSE
                  IF ( .NOT. EqualRealNos( InputFileData%OoPDefl, 0.0_ReKi ) ) THEN
                     CALL CheckError( ErrID_Warn, BadOoP )
                  END IF

                  IF ( .NOT. EqualRealNos( A(2,1), 0.0_ReKi ) )  THEN
                     InitQF2(K) = A(2,3)/A(2,1)
                  ELSE
                     InitQF2(K) = 0.0

                     IF ( .NOT. EqualRealNos( InputFileData%IPDefl, 0.0_ReKi ) )  THEN
                        CALL CheckError( ErrID_Warn, BadIP )
                     ENDIF

                  ENDIF
               ENDIF

            ELSE                                  ! It is impossible to find any "good" solution, so ignore
                                                  ! the initial tip displacements.
               InitQF2(K) = 0.0
               InitQE1(K) = 0.0

               IF ( .NOT. EqualRealNos( InputFileData%OoPDefl,  0.0_ReKi ) .OR. &
                    .NOT. EqualRealNos( InputFileData%IPDefl,   0.0_ReKi ) )  THEN
                  CALL CheckError( ErrID_Warn, Ignore )
                ENDIF
            ENDIF

         ELSE                                     ! Blade flap mode 2 is not enabled.

            InitQF2(K) = 0.0

            IF ( .NOT. EqualRealNos( A(1,2), 0.0_ReKi ) )  THEN

               IF ( .NOT. EqualRealNos( A(2,2), 0.0_ReKi ) )  THEN         ! Find a solution of the 2 equations in 1 variable that minimizes
                                                                           !  the sum of the squares of the equation's residuals.
                  InitQE1(K) = ( A(1,2)*A(1,3) + A(2,2)*A(2,3) )/( A(1,2)**2 + A(2,2)**2 )

                  TotResid = SQRT( ( A(1,2)*InitQE1(K) - A(1,3) )**2 + ( A(2,2)*InitQE1(K) - A(2,3) )**2)

                  IF ( .NOT. EqualRealNos( TotResid, 0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, Approx )
                  ENDIF

               ELSE

                  InitQE1(K) = A(1,3)/A(1,2)

                  IF ( .NOT. EqualRealNos( InputFileData%IPDefl, 0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, BadIP )
                  ENDIF

               ENDIF

            ELSE

               IF ( .NOT. EqualRealNos( InputFileData%OoPDefl, 0.0_ReKi ) ) THEN
                  CALL CheckError( ErrID_Warn, BadOoP )
               END IF

               IF ( .NOT. EqualRealNos( A(2,2), 0.0_ReKi ) )  THEN
                  InitQE1(K) = A(2,3)/A(2,2)
               ELSE
                  InitQE1(K) = 0.0

                  IF ( .NOT. EqualRealNos( InputFileData%IPDefl,  0.0_ReKi ) )  THEN
                     CALL CheckError( ErrID_Warn, BadIP )
                  ENDIF

               ENDIF

            ENDIF

         ENDIF

      ENDIF

   END DO !K

   RETURN
CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'InitBlDefl:Blade '//TRIM(Num2LStr(K))// &
                     ' initial blade tip displacements are Incompat with enabled DOFs: '//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE InitBlDefl
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetEnabledDOFIndexArrays( p )
! This routine is used create arrays of DOF indices (pointers / (vector susbscript arrays) that contribute to the QD2T-related
!   linear accelerations of various points within the system in the inertia frame, based on which DOFs are presently enabled.
! NOTE: The order in which the DOFs are tested within this routine and hence the order in which the DOF indices appear in the
!       vector subscript arrays, determines the order in which the states will appear in the linearized model created by FAST
!       when AnalMode == 2.  This order is not necessarily sorted from smallest to largest DOF index.
! bjj: note that this routine is now called only in the initialization routine. It is not available during time simulation.
!----------------------------------------------------------------------------------------------------------------------------------

   IMPLICIT                        NONE


      ! passed variables
   TYPE(ED_ParameterType), INTENT(INOUT)   :: p                                  ! Parameters of the structural dynamics module

      ! Local Variables:
   INTEGER(IntKi)                   :: I                                          ! Loops through all DOFs.
   INTEGER(IntKi)                   :: K                                          ! Loops through blades.



      ! Initialize total counts to zero.

   p%DOFs%NActvDOF = 0
   p%DOFs%NPCE     = 0
   p%DOFs%NPDE     = 0
   p%DOFs%NPIE     = 0
   p%DOFs%NPTTE    = 0
   p%DOFs%NPTE     = 0
   p%DOFs%NPSBE(:) = 0
   p%DOFs%NPSE (:) = 0
   p%DOFs%NPUE     = 0
   p%DOFs%NPYE     = 0


      ! Test each DOF and include the appropriate indices in the subscript arrays
      !  and total counts:

   IF ( p%DOF_Flag(DOF_Sg  ) )  THEN  ! Platform surge.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Sg
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Sg
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_Sg
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_Sg
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_Sg
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Sg
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_Sg
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_Sg

   ENDIF


   IF ( p%DOF_Flag(DOF_Sw  ) )  THEN  ! Platform sway.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

       p%DOFs%PS     (  p%DOFs%NActvDOF) = DOF_Sw
       p%DOFs%PCE    (  p%DOFs%NPCE    ) = DOF_Sw
       p%DOFs%PDE    (  p%DOFs%NPDE    ) = DOF_Sw
       p%DOFs%PIE    (  p%DOFs%NPIE    ) = DOF_Sw
       p%DOFs%PTE    (  p%DOFs%NPTE    ) = DOF_Sw
       p%DOFs%PSE    (:,p%DOFs%NPSE (:)) = DOF_Sw
       p%DOFs%PUE    (  p%DOFs%NPUE    ) = DOF_Sw
       p%DOFs%PYE    (  p%DOFs%NPYE    ) = DOF_Sw

   ENDIF


   IF ( p%DOF_Flag(DOF_Hv  ) )  THEN  ! Platform heave.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Hv
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Hv
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_Hv
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_Hv
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_Hv
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Hv
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_Hv
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_Hv

   ENDIF


   IF ( p%DOF_Flag(DOF_R   ) )  THEN  ! Platform roll.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_R
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_R
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_R
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_R
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_R
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_R
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_R
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_R

   ENDIF


   IF ( p%DOF_Flag(DOF_P   ) )  THEN  ! Platform pitch.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_P
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_P
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_P
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_P
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_P
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_P
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_P
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_P

   ENDIF


   IF ( p%DOF_Flag(DOF_Y   ) )  THEN  ! Platform yaw.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1
      p%DOFs%NPYE     = p%DOFs%NPYE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Y
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Y
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_Y
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_Y
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_Y
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Y
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_Y
      p%DOFs%PYE     (  p%DOFs%NPYE    ) = DOF_Y

   ENDIF


   IF ( p%DOF_Flag(DOF_TFA1) )  THEN  ! 1st tower fore-aft.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTTE    = p%DOFs%NPTTE    + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TFA1
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_TFA1
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_TFA1
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TFA1
      p%DOFs%PTTE    (  p%DOFs%NPTTE   ) = DOF_TFA1
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_TFA1
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_TFA1
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_TFA1

   ENDIF


   IF ( p%DOF_Flag(DOF_TSS1) )  THEN  ! 1st tower side-to-side.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTTE    = p%DOFs%NPTTE    + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TSS1
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_TSS1
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_TSS1
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TSS1
      p%DOFs%PTTE    (  p%DOFs%NPTTE   ) = DOF_TSS1
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_TSS1
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_TSS1
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_TSS1

   ENDIF


   IF ( p%DOF_Flag(DOF_TFA2) )  THEN  ! 2nd tower fore-aft.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTTE    = p%DOFs%NPTTE    + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TFA2
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_TFA2
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_TFA2
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TFA2
      p%DOFs%PTTE    (  p%DOFs%NPTTE   ) = DOF_TFA2
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_TFA2
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_TFA2
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_TFA2

   ENDIF


   IF ( p%DOF_Flag(DOF_TSS2) )  THEN  ! 2nd tower side-to-side.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPTTE    = p%DOFs%NPTTE    + 1
      p%DOFs%NPTE     = p%DOFs%NPTE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TSS2
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_TSS2
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_TSS2
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TSS2
      p%DOFs%PTTE    (  p%DOFs%NPTTE   ) = DOF_TSS2
      p%DOFs%PTE     (  p%DOFs%NPTE    ) = DOF_TSS2
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_TSS2
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_TSS2

   ENDIF


   IF ( p%DOF_Flag(DOF_Yaw ) )  THEN  ! Nacelle yaw.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1
      p%DOFs%NPUE     = p%DOFs%NPUE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Yaw
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Yaw
      p%DOFs%PDE     (  p%DOFs%NPDE    ) = DOF_Yaw
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_Yaw
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Yaw
      p%DOFs%PUE     (  p%DOFs%NPUE    ) = DOF_Yaw

   ENDIF


   IF ( p%DOF_Flag(DOF_TFrl) )  THEN  ! Tail-furl.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPIE     = p%DOFs%NPIE     + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_TFrl
      p%DOFs%PIE     (  p%DOFs%NPIE    ) = DOF_TFrl

   ENDIF


   IF ( p%DOF_Flag(DOF_RFrl) )  THEN  ! Rotor-furl.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPDE     = p%DOFs%NPDE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1

      p%DOFs%PS     (  p%DOFs%NActvDOF) = DOF_RFrl
      p%DOFs%PCE    (  p%DOFs%NPCE    ) = DOF_RFrl
      p%DOFs%PDE    (  p%DOFs%NPDE    ) = DOF_RFrl
      p%DOFs%PSE    (:,p%DOFs%NPSE (:)) = DOF_RFrl

   ENDIF


   IF ( p%DOF_Flag(DOF_GeAz) )  THEN  ! Generator azimuth.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_GeAz
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_GeAz
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_GeAz

   ENDIF


   IF ( p%DOF_Flag(DOF_DrTr) )  THEN  ! Drivetrain torsion.

      p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
      p%DOFs%NPCE     = p%DOFs%NPCE     + 1
      p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1

      p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_DrTr
      p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_DrTr
      p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_DrTr

   ENDIF


   IF ( p%NumBl == 2 )  THEN
      IF ( p%DOF_Flag(DOF_Teet   ) )  THEN  ! Rotor-teeter.

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
         p%DOFs%NPCE     = p%DOFs%NPCE     + 1
         p%DOFs%NPSE (:) = p%DOFs%NPSE (:) + 1

         p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_Teet
         p%DOFs%PCE     (  p%DOFs%NPCE    ) = DOF_Teet
         p%DOFs%PSE     (:,p%DOFs%NPSE (:)) = DOF_Teet

      ENDIF
   ENDIF


   DO K = 1,p%NumBl ! Loop through all blades
      IF ( p%DOF_Flag(DOF_BF(K,1)) )  THEN  ! 1st blade flap.

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
         p%DOFs%NPSBE(K) = p%DOFs%NPSBE(K) + 1
         p%DOFs%NPSE (K) = p%DOFs%NPSE (K) + 1

         p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_BF(K,1)
         p%DOFs%PSBE    (K,p%DOFs%NPSBE(K)) = DOF_BF(K,1)
         p%DOFs%PSE     (K,p%DOFs%NPSE (K)) = DOF_BF(K,1)

      ENDIF
   ENDDO          ! K - Blades


   DO K = 1,p%NumBl ! Loop through all blades
      IF ( p%DOF_Flag(DOF_BE(K,1)) )  THEN  ! 1st blade edge.

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
         p%DOFs%NPSBE(K) = p%DOFs%NPSBE(K) + 1
         p%DOFs%NPSE (K) = p%DOFs%NPSE (K) + 1

         p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_BE(K,1)
         p%DOFs%PSBE    (K,p%DOFs%NPSBE(K)) = DOF_BE(K,1)
         p%DOFs%PSE     (K,p%DOFs%NPSE (K)) = DOF_BE(K,1)

      ENDIF
   ENDDO          ! K - Blades


   DO K = 1,p%NumBl ! Loop through all blades
      IF ( p%DOF_Flag(DOF_BF(K,2)) )  THEN  ! 2nd blade flap.

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1
         p%DOFs%NPSBE(K) = p%DOFs%NPSBE(K) + 1
         p%DOFs%NPSE (K) = p%DOFs%NPSE (K) + 1

         p%DOFs%PS      (  p%DOFs%NActvDOF) = DOF_BF(K,2)
         p%DOFs%PSBE    (K,p%DOFs%NPSBE(K)) = DOF_BF(K,2)
         p%DOFs%PSE     (K,p%DOFs%NPSE (K)) = DOF_BF(K,2)

      ENDIF
   ENDDO          ! K - Blades



      ! Compute the sorted (from smallest to largest p%DOFs index) version of PS(),
      !   SrtPS(), and SrtPSNAUG().  At the same time compute Diag(), which is an
      !   array containing the indices of SrtPS() associated with each enabled
      !   DOF; that is, SrtPS(Diag(I)) = I:
      ! NOTE: This calculation is recomputing NActvDOF as computed above.  This is
      !       of no concern however, since the resulting value will be the same.

   p%DOFs%NActvDOF = 0
   DO I = 1,p%NDOF  ! Loop through all DOFs
      IF ( p%DOF_Flag(I) )  THEN   ! .TRUE. if the corresponding DOF is enabled

         p%DOFs%NActvDOF = p%DOFs%NActvDOF + 1

         p%DOFs%SrtPS     (p%DOFs%NActvDOF) = I
         p%DOFs%SrtPSNAUG (p%DOFs%NActvDOF) = I
         p%DOFs%Diag      (I           ) = p%DOFs%NActvDOF

      ENDIF
   ENDDO          ! I - All DOFs

   p%DOFs%SrtPSNAUG ( p%DOFs%NActvDOF + 1 ) = p%NAug


   RETURN
END SUBROUTINE SetEnabledDOFIndexArrays
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetCoordSy( t, CoordSys, RtHSdat, BlPitch, p, x, ErrStat, ErrMsg )


   ! This routine is used to define the internal coordinate systems for this particular time step.
   ! It also sets the TeeterAng and TeetAngVel for this time step.


   IMPLICIT                        NONE

      ! Subroutine arguments (passed variables)

   REAL(DbKi),                   INTENT(IN)    :: t                             ! Current simulation time, in seconds (used only for SmllRotTrans error messages)
   REAL(ReKi),                   INTENT(IN)    :: BlPitch (:)                   ! The current blade pitch
   TYPE(ED_CoordSys),            INTENT(INOUT) :: CoordSys                      ! The coordinate systems to be set
   TYPE(ED_RtHndSide),           INTENT(INOUT) :: RtHSdat                       ! data from the RtHndSid module
   TYPE(ED_ParameterType),       INTENT(IN)    :: p                             ! The module's parameters
   TYPE(ED_ContinuousStateType), INTENT(IN)    :: x                             ! The module's continuous states

   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat                      ! Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg                       ! Error message

      ! Local variables

   REAL(R8Ki)                   :: CAzimuth                                        ! COS( rotor azimuth angle ).
   REAL(R8Ki)                   :: CgRotAng                                        ! COS( gRotAng ).
   REAL(R8Ki)                   :: CNacYaw                                         ! COS( nacelle yaw angle ).
   REAL(R8Ki)                   :: CosPitch                                        ! COS( the current pitch angle ).
   REAL(R8Ki)                   :: CPitPTwstA                                      ! COS( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of cosine.
   REAL(R8Ki)                   :: CPitPTwstS                                      ! COS( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of cosine.
   REAL(R8Ki)                   :: CRotFurl                                        ! COS( rotor-furl angle ).
   REAL(R8Ki)                   :: CTailFurl                                       ! COS( tail-furl angle ).
   REAL(R8Ki)                   :: CTeetAng                                        ! COS( TeetAng ).
   REAL(R8Ki)                   :: g1Prime   (3)                                   ! = g1.
   REAL(R8Ki)                   :: g2Prime   (3)                                   ! completes the right-handed gPrime-vector triad
   REAL(R8Ki)                   :: g3Prime   (3)                                   ! = g3 rotated about g1 so that parallel to the pitching axis of blade K (i.e., the current blade in the blade loop).
   REAL(R8Ki)                   :: gRotAng                                         ! Angle of rotation about g1 to get from the g to the gPrime system.
   REAL(R8Ki)                   :: Lj1       (3)                                   ! vector / direction Lj1 at node J for blade K.
   REAL(R8Ki)                   :: Lj2       (3)                                   ! vector / direction Lj2 at node J for blade K.
   REAL(R8Ki)                   :: Lj3       (3)                                   ! vector / direction Lj3 at node J for blade K.
   REAL(R8Ki)                   :: SAzimuth                                        ! SIN( rotor azimuth angle ).
   REAL(R8Ki)                   :: SgRotAng                                        ! SIN( gRotAng ).
   REAL(R8Ki)                   :: SinPitch                                        ! SIN( the current pitch angle ).
   REAL(R8Ki)                   :: SNacYaw                                         ! SIN( nacelle yaw angle ).
   REAL(R8Ki)                   :: SPitPTwstA                                      ! SIN( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of sine.
   REAL(R8Ki)                   :: SPitPTwstS                                      ! SIN( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of sine.
   REAL(R8Ki)                   :: SRotFurl                                        ! SIN( rotor-furl angle ).
   REAL(R8Ki)                   :: STailFurl                                       ! SIN( tail-furl angle ).
   REAL(R8Ki)                   :: STeetAng                                        ! SIN( TeetAng ).
   REAL(R8Ki)                   :: ThetaFA                                         ! Tower fore-aft tilt deflection angle.
   REAL(R8Ki)                   :: ThetaIP                                         ! Blade in-plane deflection angle at node J for blade K.
   REAL(R8Ki)                   :: ThetaLxb                                        ! Blade deflection angle about the Lxb (n1) -axis at node J for blade K.
   REAL(R8Ki)                   :: ThetaLyb                                        ! Blade deflection angle about the Lyb (n2) -axis at node J for blade K.
   REAL(R8Ki)                   :: ThetaOoP                                        ! Blade out-of-plane deflection angle at node J for blade K.
   REAL(R8Ki)                   :: ThetaSS                                         ! Tower side-to-side tilt deflection angle.
   REAL(R8Ki)                   :: TransMat  (3,3)                                 ! The resulting transformation matrix due to three orthogonal rotations, (-).

   INTEGER(IntKi)               :: J                                               ! Loops through nodes / elements.
   INTEGER(IntKi)               :: K                                               ! Loops through blades.


   INTEGER(IntKi)               :: ErrStat2                      ! Temporary error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                       ! Temporary error message


   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Inertial frame coordinate system:

   CoordSys%z1 = (/ 1.0, 0.0, 0.0 /)   ! Vector / direction z1 (=  xi from the IEC coord. system).
   CoordSys%z2 = (/ 0.0, 1.0, 0.0 /)   ! Vector / direction z2 (=  zi from the IEC coord. system).
   CoordSys%z3 = (/ 0.0, 0.0, 1.0 /)   ! Vector / direction z3 (= -yi from the IEC coord. system).


      ! Tower base / platform coordinate system:

   CALL SmllRotTrans( 'platform displacement (ElastoDyn SetCoordSy)', x%QT(DOF_R), x%QT(DOF_Y), -x%QT(DOF_P), TransMat, TRIM(Num2LStr(t))//' s', ErrStat2, ErrMsg2 )  ! Get the transformation matrix, TransMat, from inertial frame to tower base / platform coordinate systems.
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   CoordSys%a1 = TransMat(1,1)*CoordSys%z1 + TransMat(1,2)*CoordSys%z2 + TransMat(1,3)*CoordSys%z3 ! Vector / direction a1 (=  xt from the IEC coord. system).
   CoordSys%a2 = TransMat(2,1)*CoordSys%z1 + TransMat(2,2)*CoordSys%z2 + TransMat(2,3)*CoordSys%z3 ! Vector / direction a2 (=  zt from the IEC coord. system).
   CoordSys%a3 = TransMat(3,1)*CoordSys%z1 + TransMat(3,2)*CoordSys%z2 + TransMat(3,3)*CoordSys%z3 ! Vector / direction a3 (= -yt from the IEC coord. system).


   DO J = 1,p%TwrNodes ! Loop through the tower nodes / elements


      ! Tower element-fixed coordinate system:

      ThetaFA = -p%TwrFASF(1,J       ,1)*x%QT(DOF_TFA1) - p%TwrFASF(2,J       ,1)*x%QT(DOF_TFA2)
      ThetaSS =  p%TwrSSSF(1,J       ,1)*x%QT(DOF_TSS1) + p%TwrSSSF(2,J       ,1)*x%QT(DOF_TSS2)

      CALL SmllRotTrans( 'tower deflection (ElastoDyn SetCoordSy)', ThetaSS, 0.0_R8Ki, ThetaFA, TransMat, TRIM(Num2LStr(t))//' s', ErrStat2, ErrMsg2 )   ! Get the transformation matrix, TransMat, from tower-base to tower element-fixed coordinate systems.
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN

      CoordSys%t1(J,:) = TransMat(1,1)*CoordSys%a1 + TransMat(1,2)*CoordSys%a2 + TransMat(1,3)*CoordSys%a3  ! Vector / direction t1 for tower node J (=  Lxt from the IEC coord. system).
      CoordSys%t2(J,:) = TransMat(2,1)*CoordSys%a1 + TransMat(2,2)*CoordSys%a2 + TransMat(2,3)*CoordSys%a3  ! Vector / direction t2 for tower node J (=  Lzt from the IEC coord. system).
      CoordSys%t3(J,:) = TransMat(3,1)*CoordSys%a1 + TransMat(3,2)*CoordSys%a2 + TransMat(3,3)*CoordSys%a3  ! Vector / direction t3 for tower node J (= -Lyt from the IEC coord. system).


   ENDDO ! J - Tower nodes / elements


      ! Tower-top / base plate coordinate system:

   ThetaFA    = -p%TwrFASF(1,p%TTopNode,1)*x%QT(DOF_TFA1) - p%TwrFASF(2,p%TTopNode,1)*x%QT(DOF_TFA2)
   ThetaSS    =  p%TwrSSSF(1,p%TTopNode,1)*x%QT(DOF_TSS1) + p%TwrSSSF(2,p%TTopNode,1)*x%QT(DOF_TSS2)

   CALL SmllRotTrans( 'tower deflection (ElastoDyn SetCoordSy)', ThetaSS, 0.0_R8Ki, ThetaFA, TransMat, TRIM(Num2LStr(t))//' s', ErrStat2, ErrMsg2 )   ! Get the transformation matrix, TransMat, from tower-base to tower-top/base-plate coordinate systems.
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   CoordSys%b1 = TransMat(1,1)*CoordSys%a1 + TransMat(1,2)*CoordSys%a2 + TransMat(1,3)*CoordSys%a3 ! Vector / direction b1 (=  xp from the IEC coord. system).
   CoordSys%b2 = TransMat(2,1)*CoordSys%a1 + TransMat(2,2)*CoordSys%a2 + TransMat(2,3)*CoordSys%a3 ! Vector / direction b2 (=  zp from the IEC coord. system).
   CoordSys%b3 = TransMat(3,1)*CoordSys%a1 + TransMat(3,2)*CoordSys%a2 + TransMat(3,3)*CoordSys%a3 ! Vector / direction b3 (= -yp from the IEC coord. system).


      ! Nacelle / yaw coordinate system:

   CNacYaw  = COS( x%QT(DOF_Yaw ) )
   SNacYaw  = SIN( x%QT(DOF_Yaw ) )

   CoordSys%d1 = CNacYaw*CoordSys%b1 - SNacYaw*CoordSys%b3     ! Vector / direction d1 (=  xn from the IEC coord. system).
   CoordSys%d2 = CoordSys%b2                                   ! Vector / direction d2 (=  zn from the IEC coord. system).
   CoordSys%d3 = SNacYaw*CoordSys%b1 + CNacYaw*CoordSys%b3     ! Vector / direction d3 (= -yn from the IEC coord. system).


      ! Rotor-furl coordinate system:

   CRotFurl = COS( x%QT(DOF_RFrl) )
   SRotFurl = SIN( x%QT(DOF_RFrl) )

   CoordSys%rf1 = ( (   1.0 - p%CRFrlSkw2*p%CRFrlTlt2 )*CRotFurl   + p%CRFrlSkw2*p%CRFrlTlt2          )*CoordSys%d1 &
                + ( p%CRFrlSkew*p%CSRFrlTlt*( 1.0 -     CRotFurl ) - p%SRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d2 &
                + ( p%CSRFrlSkw*p%CRFrlTlt2*( CRotFurl - 1.0     ) -             p%SRFrlTilt*SRotFurl )*CoordSys%d3
   CoordSys%rf2 = ( p%CRFrlSkew*p%CSRFrlTlt*( 1.0 -     CRotFurl ) + p%SRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d1 &
                + (             p%CRFrlTlt2*            CRotFurl   +             p%SRFrlTlt2          )*CoordSys%d2 &
                + ( p%SRFrlSkew*p%CSRFrlTlt*( CRotFurl - 1.0     ) + p%CRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d3
   CoordSys%rf3 = ( p%CSRFrlSkw*p%CRFrlTlt2*( CRotFurl - 1.0     ) +             p%SRFrlTilt*SRotFurl )*CoordSys%d1 &
                + ( p%SRFrlSkew*p%CSRFrlTlt*( CRotFurl - 1.0     ) - p%CRFrlSkew*p%CRFrlTilt*SRotFurl )*CoordSys%d2 &
                + ( (   1.0 - p%SRFrlSkw2*p%CRFrlTlt2 )*CRotFurl   + p%SRFrlSkw2*p%CRFrlTlt2          )*CoordSys%d3
   CoordSys%rfa = p%CRFrlSkew*p%CRFrlTilt*CoordSys%d1 + p%SRFrlTilt*CoordSys%d2 - p%SRFrlSkew*p%CRFrlTilt*CoordSys%d3


      ! Shaft coordinate system:

   CoordSys%c1 =  p%CShftSkew*p%CShftTilt*CoordSys%rf1 + p%SShftTilt*CoordSys%rf2 - p%SShftSkew*p%CShftTilt*CoordSys%rf3  ! Vector / direction c1 (=  xs from the IEC coord. system).
   CoordSys%c2 = -p%CShftSkew*p%SShftTilt*CoordSys%rf1 + p%CShftTilt*CoordSys%rf2 + p%SShftSkew*p%SShftTilt*CoordSys%rf3  ! Vector / direction c2 (=  zs from the IEC coord. system).
   CoordSys%c3 =  p%SShftSkew*            CoordSys%rf1                            + p%CShftSkew*            CoordSys%rf3  ! Vector / direction c3 (= -ys from the IEC coord. system).


      ! Azimuth coordinate system:

   CAzimuth = COS( x%QT(DOF_DrTr) + x%QT(DOF_GeAz) )
   SAzimuth = SIN( x%QT(DOF_DrTr) + x%QT(DOF_GeAz) )

   CoordSys%e1 =  CoordSys%c1                                  ! Vector / direction e1 (=  xa from the IEC coord. system).
   CoordSys%e2 =  CAzimuth*CoordSys%c2 + SAzimuth*CoordSys%c3  ! Vector / direction e2 (=  ya from the IEC coord. system).
   CoordSys%e3 = -SAzimuth*CoordSys%c2 + CAzimuth*CoordSys%c3  ! Vector / direction e3 (=  za from the IEC coord. system).


      ! Teeter coordinate system:

      ! Lets define TeetAng, which is the current teeter angle (= QT(DOF_Teet) for
      !   2-blader or 0 for 3-blader) and is used in place of QT(DOF_Teet)
      !   throughout SUBROUTINE RtHS().  Doing it this way, we can run the same
      !   equations of motion for both the 2 and 3-blader configurations even
      !   though a 3-blader does not have a teetering DOF.

   IF ( p%NumBl == 2 )  THEN ! 2-blader
      RtHSdat%TeetAng    = x%QT (DOF_Teet)
      RtHSdat%TeetAngVel = x%QDT(DOF_Teet)
   ELSE                    ! 3-blader
      RtHSdat%TeetAng    = 0.0  ! Teeter is not an available DOF for a 3-blader
      RtHSdat%TeetAngVel = 0.0  ! Teeter is not an available DOF for a 3-blader
   ENDIF
   CTeetAng = COS( RtHSdat%TeetAng )
   STeetAng = SIN( RtHSdat%TeetAng )

   CoordSys%f1 = CTeetAng*CoordSys%e1 - STeetAng*CoordSys%e3       ! Vector / direction f1.
   CoordSys%f2 = CoordSys%e2                                       ! Vector / direction f2.
   CoordSys%f3 = STeetAng*CoordSys%e1 + CTeetAng*CoordSys%e3       ! Vector / direction f3.


      ! Hub / delta-3 coordinate system:

   CoordSys%g1 =  CoordSys%f1                                      ! Vector / direction g1 (=  xh from the IEC coord. system).
   CoordSys%g2 =  p%CosDel3*CoordSys%f2 + p%SinDel3*CoordSys%f3    ! Vector / direction g2 (=  yh from the IEC coord. system).
   CoordSys%g3 = -p%SinDel3*CoordSys%f2 + p%CosDel3*CoordSys%f3    ! Vector / direction g3 (=  zh from the IEC coord. system).


   DO K = 1,p%NumBl ! Loop through all blades


      ! Hub (Prime) coordinate system rotated to match blade K.

       gRotAng = p%TwoPiNB*(K-1)
      CgRotAng = COS( gRotAng )
      SgRotAng = SIN( gRotAng )

      g1Prime =  CoordSys%g1
      g2Prime =  CgRotAng*CoordSys%g2 + SgRotAng*CoordSys%g3
      g3Prime = -SgRotAng*CoordSys%g2 + CgRotAng*CoordSys%g3


      ! Coned coordinate system:

      CoordSys%i1(K,:) = p%CosPreC(K)*g1Prime - p%SinPreC(K)*g3Prime  ! i1(K,:) = vector / direction i1 for blade K (=  xcK from the IEC coord. system).
      CoordSys%i2(K,:) = g2Prime                                      ! i2(K,:) = vector / direction i2 for blade K (=  ycK from the IEC coord. system).
      CoordSys%i3(K,:) = p%SinPreC(K)*g1Prime + p%CosPreC(K)*g3Prime  ! i3(K,:) = vector / direction i3 for blade K (=  zcK from the IEC coord. system).


      ! Blade / pitched coordinate system:

      CosPitch = COS( BlPitch(K) )
      SinPitch = SIN( BlPitch(K) )

      CoordSys%j1(K,:) = CosPitch*CoordSys%i1(K,:) - SinPitch*CoordSys%i2(K,:)      ! j1(K,:) = vector / direction j1 for blade K (=  xbK from the IEC coord. system).
      CoordSys%j2(K,:) = SinPitch*CoordSys%i1(K,:) + CosPitch*CoordSys%i2(K,:)      ! j2(K,:) = vector / direction j2 for blade K (=  ybK from the IEC coord. system).
      CoordSys%j3(K,:) = CoordSys%i3(K,:)                                           ! j3(K,:) = vector / direction j3 for blade K (=  zbK from the IEC coord. system).


      DO J = 0,p%TipNode ! Loop through the blade nodes / elements


      ! Blade coordinate system aligned with local structural axes (not element fixed):

         Lj1 = p%CThetaS(K,J)*CoordSys%j1(K,:) - p%SThetaS(K,J)*CoordSys%j2(K,:)  ! vector / direction Lj1 at node J for blade K
         Lj2 = p%SThetaS(K,J)*CoordSys%j1(K,:) + p%CThetaS(K,J)*CoordSys%j2(K,:)  ! vector / direction Lj2 at node J for blade K
         Lj3 = CoordSys%j3(K,:)                                               ! vector / direction Lj3 at node J for blade K


      ! Blade element-fixed coordinate system aligned with local structural axes:

         ThetaOoP =   p%TwistedSF(K,1,1,J,1)*x%QT( DOF_BF(K,1) ) &
                    + p%TwistedSF(K,1,2,J,1)*x%QT( DOF_BF(K,2) ) &
                    + p%TwistedSF(K,1,3,J,1)*x%QT( DOF_BE(K,1) )
         ThetaIP  = - p%TwistedSF(K,2,1,J,1)*x%QT( DOF_BF(K,1) ) &
                    - p%TwistedSF(K,2,2,J,1)*x%QT( DOF_BF(K,2) ) &
                    - p%TwistedSF(K,2,3,J,1)*x%QT( DOF_BE(K,1) )

         ThetaLxb = p%CThetaS(K,J)*ThetaIP - p%SThetaS(K,J)*ThetaOoP
         ThetaLyb = p%SThetaS(K,J)*ThetaIP + p%CThetaS(K,J)*ThetaOoP

         CALL SmllRotTrans( 'blade deflection (ElastoDyn SetCoordSy)', ThetaLxb, ThetaLyb, 0.0_R8Ki, TransMat, TRIM(Num2LStr(t))//' s', ErrStat2, ErrMsg2 ) ! Get the transformation matrix, TransMat, from blade coordinate system aligned with local structural axes (not element fixed) to blade element-fixed coordinate system aligned with local structural axes.
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN

         CoordSys%n1(K,J,:) = TransMat(1,1)*Lj1 + TransMat(1,2)*Lj2 + TransMat(1,3)*Lj3   ! Vector / direction n1 for node J of blade K (= LxbK from the IEC coord. system).
         CoordSys%n2(K,J,:) = TransMat(2,1)*Lj1 + TransMat(2,2)*Lj2 + TransMat(2,3)*Lj3   ! Vector / direction n2 for node J of blade K (= LybK from the IEC coord. system).
         CoordSys%n3(K,J,:) = TransMat(3,1)*Lj1 + TransMat(3,2)*Lj2 + TransMat(3,3)*Lj3   ! Vector / direction n3 for node J of blade K (= LzbK from the IEC coord. system).

      ! skip these next CoordSys variables at the root and the tip; they are required only for AD14:
         
         if (j == 0 .or. j==p%TipNode) cycle  
      
         
      ! Blade element-fixed coordinate system used for calculating and returning
      !    aerodynamics loads:
      ! This coordinate system is rotated about positive n3 by the angle
      !    BlPitch(K) + ThetaS(K,J) and is coincident with the i-vector triad
      !    when the blade is undeflected.

         CPitPTwstS = CosPitch*p%CThetaS(K,J) - SinPitch*p%SThetaS(K,J)  ! = COS( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of cosine.
         SPitPTwstS = CosPitch*p%SThetaS(K,J) + SinPitch*p%CThetaS(K,J)  ! = SIN( BlPitch(K) + ThetaS(K,J) ) found using the sum of angles formulae of   sine.

         CoordSys%m1(K,J,:)  =  CPitPTwstS*CoordSys%n1(K,J,:) + SPitPTwstS*CoordSys%n2(K,J,:)   ! m1(K,J,:) = vector / direction m1 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
         CoordSys%m2(K,J,:)  = -SPitPTwstS*CoordSys%n1(K,J,:) + CPitPTwstS*CoordSys%n2(K,J,:)   ! m2(K,J,:) = vector / direction m2 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
         CoordSys%m3(K,J,:)  =  CoordSys%n3(K,J,:)                                              ! m3(K,J,:) = vector / direction m3 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).


      ! Calculate the trailing edge coordinate system used in noise calculations.
      ! This coordinate system is blade element-fixed and oriented with the local
      !   aerodynamic axes (te2 points toward trailing edge, te1 points toward
      !   suction surface):

         CPitPTwstA = CosPitch*p%CAeroTwst(J) - SinPitch*p%SAeroTwst(J)  ! = COS( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of cosine.
         SPitPTwstA = CosPitch*p%SAeroTwst(J) + SinPitch*p%CAeroTwst(J)  ! = SIN( BlPitch(K) + AeroTwst(J) ) found using the sum of angles formulae of   sine.

         CoordSys%te1(K,J,:) =  CPitPTwstA*CoordSys%m1(K,J,:) - SPitPTwstA*CoordSys%m2(K,J,:)   ! te1(K,J,:) = vector / direction te1 for node J of blade K (used to calc. noise and to calc. and return aerodynamic loads from AeroDyn).
         CoordSys%te2(K,J,:) =  SPitPTwstA*CoordSys%m1(K,J,:) + CPitPTwstA*CoordSys%m2(K,J,:)   ! te2(K,J,:) = vector / direction te2 for node J of blade K (used to calc. noise and to calc. and return aerodynamic loads from AeroDyn).
         CoordSys%te3(K,J,:) =  CoordSys%m3(K,J,:)                                              ! te3(K,J,:) = vector / direction te3 for node J of blade K (used to calc. noise and to calc. and return aerodynamic loads from AeroDyn).


      ENDDO ! J - Blade nodes / elements


   ENDDO ! K - Blades


      ! Tail-furl coordinate system:

   CTailFurl = COS( x%QT(DOF_TFrl) )
   STailFurl = SIN( x%QT(DOF_TFrl) )

   CoordSys%tf1 = ( ( 1.0 - p%CTFrlSkw2*p%CTFrlTlt2 )*CTailFurl  + p%CTFrlSkw2*p%CTFrlTlt2           )*CoordSys%d1 &
                + ( p%CTFrlSkew*p%CSTFrlTlt*(  1.0 - CTailFurl ) - p%STFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d2 &
                + ( p%CSTFrlSkw*p%CTFrlTlt2*( CTailFurl - 1.0  ) -             p%STFrlTilt*STailFurl )*CoordSys%d3
   CoordSys%tf2 = ( p%CTFrlSkew*p%CSTFrlTlt*(  1.0 - CTailFurl ) + p%STFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d1 &
                + (             p%CTFrlTlt2*         CTailFurl +               p%STFrlTlt2           )*CoordSys%d2 &
                + ( p%STFrlSkew*p%CSTFrlTlt*( CTailFurl - 1.0  ) + p%CTFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d3
   CoordSys%tf3 = ( p%CSTFrlSkw*p%CTFrlTlt2*( CTailFurl - 1.0  ) +             p%STFrlTilt*STailFurl )*CoordSys%d1 &
                + ( p%STFrlSkew*p%CSTFrlTlt*( CTailFurl - 1.0  ) - p%CTFrlSkew*p%CTFrlTilt*STailFurl )*CoordSys%d2 &
                + ( ( 1.0 - p%STFrlSkw2*p%CTFrlTlt2 )*CTailFurl  + p%STFrlSkw2*p%CTFrlTlt2           )*CoordSys%d3
   CoordSys%tfa = p%CTFrlSkew*p%CTFrlTilt*CoordSys%d1 + p%STFrlTilt*CoordSys%d2 - p%STFrlSkew*p%CTFrlTilt*CoordSys%d3


      ! Tail fin coordinate system:

   CoordSys%p1 = (                           p%CTFinSkew*p%CTFinTilt             )*CoordSys%tf1 &   ! Vector / direction p1 (= tail fin  x).
               + (                                       p%STFinTilt             )*CoordSys%tf2 &
               + (                         - p%STFinSkew*p%CTFinTilt             )*CoordSys%tf3
   CoordSys%p2 = ( p%STFinSkew*p%STFinBank - p%CTFinSkew*p%STFinTilt*p%CTFinBank )*CoordSys%tf1 &   ! Vector / direction p2 (= tail fin  z).
               + (                                       p%CTFinTilt*p%CTFinBank )*CoordSys%tf2 &
               + ( p%CTFinSkew*p%STFinBank + p%STFinSkew*p%STFinTilt*p%CTFinBank )*CoordSys%tf3
   CoordSys%p3 = ( p%STFinSkew*p%CTFinBank + p%CTFinSkew*p%STFinTilt*p%STFinBank )*CoordSys%tf1 &   ! Vector / direction p3 (= tail fin -y).
               + (                         -             p%CTFinTilt*p%STFinBank )*CoordSys%tf2 &
               + ( p%CTFinSkew*p%CTFinBank - p%STFinSkew*p%STFinTilt*p%STFinBank )*CoordSys%tf3

   RETURN
CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'SetCoordSy:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
         END IF

      END IF


   END SUBROUTINE CheckError
!----------------------------------------------------------------------------------------------------------------------------------
END SUBROUTINE SetCoordSy
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE RFurling( t, p, RFrlDef, RFrlRate, RFrlMom )
! This routine computes the rotor-furl moment due to rotor-furl deflection and rate.
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed Variables:
   REAL(DbKi), INTENT(IN)              :: t                                   ! simulation time
   TYPE(ED_ParameterType), INTENT(IN)  :: p                                   ! parameters from the structural dynamics module

   REAL(R8Ki), INTENT(IN )             :: RFrlDef                             ! The rotor-furl deflection, x%QT(DOF_RFrl)
   REAL(ReKi), INTENT(OUT)             :: RFrlMom                             ! The total moment supplied by the springs, and dampers
   REAL(R8Ki), INTENT(IN )             :: RFrlRate                            ! The rotor-furl rate, x%QDT(DOF_RFrl)


      ! Local variables:
   REAL(ReKi)                   :: RFrlDMom                                   ! The moment supplied by the rotor-furl dampers
   REAL(ReKi)                   :: RFrlSMom                                   ! The moment supplied by the rotor-furl springs


   SELECT CASE ( p%RFrlMod ) ! Which rotor-furl model are we using?

      CASE ( 0_IntKi )       ! None!


         RFrlMom = 0.0


      CASE ( 1_IntKi )        ! Standard (using inputs from the FAST furling input file).


         ! Linear spring:

         RFrlSMom = -p%RFrlSpr*RFrlDef


         ! Add spring-stops:

         IF ( RFrlDef > p%RFrlUSSP )  THEN       ! Up-stop
            RFrlSMom = RFrlSMom - p%RFrlUSSpr*( RFrlDef - p%RFrlUSSP )
         ELSEIF ( RFrlDef < p%RFrlDSSP )  THEN   ! Down-stop
            RFrlSMom = RFrlSMom - p%RFrlDSSpr*( RFrlDef - p%RFrlDSSP )
         ENDIF


         ! Linear damper:

         RFrlDMom = -p%RFrlDmp*RFrlRate


         ! Add coulomb friction:

         IF ( RFrlRate /= 0.0 )  THEN
            RFrlDMom = RFrlDMom - SIGN( p%RFrlCDmp, real(RFrlRate,ReKi) )
         ENDIF


         ! Add damper-stops:

         IF ( RFrlDef > p%RFrlUSDP )  THEN       ! Up-stop
            RFrlDMom = RFrlDMom - p%RFrlUSDmp*RFrlRate
         ELSEIF ( RFrlDef < p%RFrlDSDP )  THEN   ! Down-stop
            RFrlDMom = RFrlDMom - p%RFrlDSDmp*RFrlRate
         ENDIF


         ! Total up all the moments.

         RFrlMom = RFrlSMom + RFrlDMom


      CASE ( 2_IntKi )              ! User-defined rotor-furl spring/damper model.


         CALL UserRFrl ( RFrlDef, RFrlRate, t, p%RootName, RFrlMom )


   END   SELECT

   RETURN
END SUBROUTINE RFurling
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Teeter( t, p, TeetDef, TeetRate, TeetMom )
! This routine computes the teeter moment due to teeter deflection and rate.
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed Variables:
   REAL(DbKi), INTENT(IN)             :: t                                       ! simulation time
   TYPE(ED_ParameterType), INTENT(IN) :: p                                       ! parameters from the structural dynamics module
   REAL(R8Ki), INTENT(IN )            :: TeetDef                                 ! The teeter deflection, x%QT(DOF_Teet).
   REAL(ReKi), INTENT(OUT)            :: TeetMom                                 ! The total moment supplied by the stop, spring, and damper.
   REAL(R8Ki), INTENT(IN )            :: TeetRate                                ! The teeter rate, x%QDT(DOF_Teet).


      ! Local variables:
   REAL(ReKi)                         :: AbsDef                                   ! Absolute value of the teeter deflection.
   REAL(ReKi)                         :: SprgDef                                  ! Deflection past the spring.
   REAL(ReKi)                         :: StopDef                                  ! Deflection past the stop.
   REAL(ReKi)                         :: TeetDMom                                 ! The moment supplied by the damper.
   REAL(ReKi)                         :: TeetFMom                                 ! The moment supplied by Coulomb-friction damping.
   REAL(ReKi)                         :: TeetKMom                                 ! The moment supplied by the spring.
   REAL(ReKi)                         :: TeetSMom                                 ! The moment supplied by the stop.



   SELECT CASE ( p%TeetMod ) ! Which teeter model are we using?

   CASE ( 0_IntKi )              ! None!


      TeetMom = 0.0_ReKi


   CASE ( 1_IntKi )              ! Standard (using inputs from the primary FAST input file).


      ! Compute the absulute value of the deflection.

      AbsDef  = ABS( TeetDef )


      ! Linear teeter spring.

      SprgDef = AbsDef - p%TeetSStP

      IF ( SprgDef > 0.0_ReKi )  THEN
         TeetKMom = -SIGN( SprgDef*p%TeetSSSp, real(TeetDef,ReKi) )
      ELSE
         TeetKMom = 0.0_ReKi
      ENDIF


      ! Compute teeter-stop moment if hard stop has been contacted.

      StopDef = AbsDef - p%TeetHStP

      IF ( StopDef > 0.0_ReKi )  THEN
         TeetSMom = -p%TeetHSSp*SIGN( StopDef, real(TeetDef,reKi) )
      ELSE
         TeetSMom = 0.0_ReKi
      ENDIF


      ! Compute linear teeter-damper moment.

      IF ( ABS(TeetDef) > p%TeetDmpP ) THEN
         TeetDMom = -p%TeetDmp*TeetRate
      ELSE
         TeetDMom = 0.0_ReKi
      END IF
      
      

      ! Add coulomb friction to the teeter hinge.

      IF ( .NOT. EqualRealNos( TeetRate, 0.0_R8Ki ) )  THEN
         TeetFMom = 0.0_ReKi
      ELSE
         TeetFMom = -SIGN( p%TeetCDmp, real(TeetRate,reKi) )
      ENDIF


      ! Total up all the moments.

      TeetMom = TeetSMom + TeetDMom + TeetKMom + TeetFMom


   CASE ( 2_IntKi )              ! User-defined teeter spring/damper model.


      CALL UserTeet ( TeetDef, TeetRate, t, p%RootName, TeetMom )


   END SELECT


   RETURN
END SUBROUTINE Teeter
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE TFurling( t, p, TFrlDef, TFrlRate, TFrlMom )
! This routine computes the tail-furl moment due to tail-furl deflection and rate.
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed Variables:
   REAL(DbKi), INTENT(IN)             :: t ! simulation time
   TYPE(ED_ParameterType), INTENT(IN) :: p                                       ! parameters from the structural dynamics module

   REAL(R8Ki), INTENT(IN )            :: TFrlDef                                 ! The tail-furl deflection, QT(DOF_TFrl).
   REAL(ReKi), INTENT(OUT)            :: TFrlMom                                 ! The total moment supplied by the springs, and dampers.
   REAL(R8Ki), INTENT(IN )            :: TFrlRate                                ! The tail-furl rate, QDT(DOF_TFrl).


      ! Local variables:

   REAL(ReKi)                         :: TFrlDMom                                ! The moment supplied by the tail-furl dampers.
   REAL(ReKi)                         :: TFrlSMom                                ! The moment supplied by the tail-furl springs.



   SELECT CASE ( p%TFrlMod ) ! Which tail-furl model are we using?

      CASE ( 0_IntKi )              ! None!


         TFrlMom = 0.0


      CASE ( 1_IntKi )              ! Standard (using inputs from the FAST furling input file).


         ! Linear spring:

         TFrlSMom = -p%TFrlSpr*TFrlDef


         ! Add spring-stops:

         IF ( TFrlDef > p%TFrlUSSP )  THEN      ! Up-stop
            TFrlSMom = TFrlSMom - p%TFrlUSSpr*( TFrlDef - p%TFrlUSSP )
         ELSEIF ( TFrlDef < p%TFrlDSSP )  THEN  ! Down-stop
            TFrlSMom = TFrlSMom - p%TFrlDSSpr*( TFrlDef - p%TFrlDSSP )
         ENDIF


         ! Linear damper:

         TFrlDMom = -p%TFrlDmp*TFrlRate


         ! Add coulomb friction:

         IF ( .NOT. EqualRealNos( TFrlRate, 0.0_R8Ki) )  THEN
            TFrlDMom = TFrlDMom - SIGN( p%TFrlCDmp, real(TFrlRate,reKi) )
         ENDIF


         ! Add damper-stops:

         IF ( TFrlDef > p%TFrlUSDP )  THEN      ! Up-stop
            TFrlDMom = TFrlDMom - p%TFrlUSDmp*TFrlRate
         ELSEIF ( TFrlDef < p%TFrlDSDP )  THEN  ! Down-stop
            TFrlDMom = TFrlDMom - p%TFrlDSDmp*TFrlRate
         ENDIF


         ! Total up all the moments.

         TFrlMom = TFrlSMom + TFrlDMom


      CASE ( 2 )              ! User-defined tail-furl spring/damper model.


         CALL UserTFrl ( TFrlDef, TFrlRate, t, p%RootName, TFrlMom )


   END SELECT


   RETURN
END SUBROUTINE TFurling
!----------------------------------------------------------------------------------------------------------------------------------
FUNCTION SignLSSTrq( p, OtherState )
! This function calculates the sign (+/-1) of the low-speed shaft torque for
   !   this time step.  MomLPRot is the moment on the
   !   low-speed shaft at the teeter pin caused by the rotor.

      ! Passed variables

   TYPE(ED_ParameterType),  INTENT(IN)  :: p                 ! Parameters
   TYPE(ED_OtherStateType), INTENT(IN)  :: OtherState        ! Initial other/optimization states

   INTEGER(IntKi)                       :: SignLSSTrq        ! The sign of the LSS_Trq, output from this function

      ! Local variables

   REAL(ReKi)                           :: MomLPRot  (3)     ! The total moment on the low-speed shaft at point P caused by the rotor.
   INTEGER(IntKi)                       :: I                 ! loop counter


   MomLPRot = OtherState%RtHS%MomLPRott ! Initialize MomLPRot using MomLPRott
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      MomLPRot = MomLPRot + OtherState%RtHS%PMomLPRot(:,p%DOFs%SrtPS(I))*OtherState%QD2T(p%DOFs%SrtPS(I))  ! Add the moments associated with the accelerations of the DOFs

   ENDDO             ! I - All active (enabled) DOFs

      ! MomLProt has now been found.  Now dot this with e1 to get the
      !   low-speed shaft torque and take the SIGN of the result:

   SignLSSTrq = NINT( SIGN( 1.0_R8Ki, DOT_PRODUCT( MomLPRot, OtherState%CoordSys%e1 ) ) )

END FUNCTION SignLSSTrq
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalculatePositions( p, x, CoordSys, RtHSdat )
! This routine is used to calculate the positions stored in other states that are used in both the
! CalcOutput and CalcContStateDeriv routines.
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    ! The coordinate systems that have been set for these states/time
   TYPE(ED_RtHndSide),           INTENT(INOUT)  :: RtHSdat     ! data from the RtHndSid module (contains positions to be set)

      !Local variables
   REAL(R8Ki)                   :: rK        (3)                                   ! Position vector from inertial frame origin to tail fin center of pressure (point K).
   !REAL(R8Ki)                   :: rQ        (3)                                   ! Position vector from inertial frame origin to apex of rotation (point Q).

   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades

      !-------------------------------------------------------------------------------------------------
      ! Positions
      !-------------------------------------------------------------------------------------------------

      ! Define the position vectors between the various points on the wind turbine
      !   that are not dependent on the distributed tower or blade parameters:

   RtHSdat%rZ    = x%QT(DOF_Sg)* CoordSys%z1 + x%QT(DOF_Hv)* CoordSys%z2 - x%QT(DOF_Sw)* CoordSys%z3                          ! Position vector from inertia frame origin to platform reference (point Z).
   RtHSdat%rZY   = p%rZYzt*  CoordSys%a2 + p%PtfmCMxt*CoordSys%a1 - p%PtfmCMyt*CoordSys%a3                                    ! Position vector from platform reference (point Z) to platform mass center (point Y).      
   RtHSdat%rZT0  = p%rZT0zt* CoordSys%a2                                                                                      ! Position vector from platform reference (point Z) to tower base (point T(0))
   RtHSdat%rZO   = ( x%QT(DOF_TFA1) + x%QT(DOF_TFA2)                                                        )*CoordSys%a1 &   ! Position vector from platform reference (point Z) to tower-top / base plate (point O).
                    + ( p%RefTwrHt - 0.5*(      p%AxRedTFA(1,1,p%TTopNode)*x%QT(DOF_TFA1)*x%QT(DOF_TFA1) &
                                          +     p%AxRedTFA(2,2,p%TTopNode)*x%QT(DOF_TFA2)*x%QT(DOF_TFA2) &
                                          + 2.0*p%AxRedTFA(1,2,p%TTopNode)*x%QT(DOF_TFA1)*x%QT(DOF_TFA2) &
                                          +     p%AxRedTSS(1,1,p%TTopNode)*x%QT(DOF_TSS1)*x%QT(DOF_TSS1) &
                                          +     p%AxRedTSS(2,2,p%TTopNode)*x%QT(DOF_TSS2)*x%QT(DOF_TSS2) &
                                          + 2.0*p%AxRedTSS(1,2,p%TTopNode)*x%QT(DOF_TSS1)*x%QT(DOF_TSS2)   ) )*CoordSys%a2 &
                    + ( x%QT(DOF_TSS1) + x%QT(DOF_TSS2)                                                      )*CoordSys%a3
   RtHSdat%rOU   =   p%NacCMxn*CoordSys%d1  +  p%NacCMzn  *CoordSys%d2  -  p%NacCMyn  *CoordSys%d3                            ! Position vector from tower-top / base plate (point O) to nacelle center of mass (point U).
   RtHSdat%rOV   = p%RFrlPntxn*CoordSys%d1  +  p%RFrlPntzn*CoordSys%d2  -  p%RFrlPntyn*CoordSys%d3                            ! Position vector from tower-top / base plate (point O) to specified point on rotor-furl axis (point V).
   RtHSdat%rVIMU =   p%rVIMUxn*CoordSys%rf1 +  p%rVIMUzn  *CoordSys%rf2 -   p%rVIMUyn *CoordSys%rf3                           ! Position vector from specified point on rotor-furl axis (point V) to nacelle IMU (point IMU).
   RtHSdat%rVD   =     p%rVDxn*CoordSys%rf1 +    p%rVDzn  *CoordSys%rf2 -     p%rVDyn *CoordSys%rf3                           ! Position vector from specified point on rotor-furl axis (point V) to center of mass of structure that furls with the rotor (not including rotor) (point D).
   RtHSdat%rVP   =     p%rVPxn*CoordSys%rf1 +    p%rVPzn  *CoordSys%rf2 -     p%rVPyn *CoordSys%rf3 + p%OverHang*CoordSys%c1  ! Position vector from specified point on rotor-furl axis (point V) to teeter pin (point P).
   RtHSdat%rPQ   = -p%UndSling*CoordSys%g1                                                                                    ! Position vector from teeter pin (point P) to apex of rotation (point Q).
   RtHSdat%rQC   =     p%HubCM*CoordSys%g1                                                                                    ! Position vector from apex of rotation (point Q) to hub center of mass (point C).
   RtHSdat%rOW   = p%TFrlPntxn*CoordSys%d1  + p%TFrlPntzn *CoordSys%d2 -  p%TFrlPntyn*CoordSys%d3                             ! Position vector from tower-top / base plate (point O) to specified point on  tail-furl axis (point W).
   RtHSdat%rWI   =     p%rWIxn*CoordSys%tf1 +      p%rWIzn*CoordSys%tf2 -     p%rWIyn*CoordSys%tf3                            ! Position vector from specified point on  tail-furl axis (point W) to tail boom center of mass     (point I).
   RtHSdat%rWJ   =     p%rWJxn*CoordSys%tf1 +      p%rWJzn*CoordSys%tf2 -     p%rWJyn*CoordSys%tf3                            ! Position vector from specified point on  tail-furl axis (point W) to tail fin  center of mass     (point J).
   RtHSdat%rWK   =     p%rWKxn*CoordSys%tf1 +      p%rWKzn*CoordSys%tf2 -     p%rWKyn*CoordSys%tf3                            ! Position vector from specified point on  tail-furl axis (point W) to tail fin  center of pressure (point K).
   RtHSdat%rPC   = RtHSdat%rPQ + RtHSdat%rQC                                                                                  ! Position vector from teeter pin (point P) to hub center of mass (point C).
   RtHSdat%rT0O  = RtHSdat%rZO - RtHSdat%rZT0                                                                                 ! Position vector from the tower base (point T(0)) to tower-top / base plate (point O).
   RtHSdat%rO    = RtHSdat%rZ  + RtHSdat%rZO                                                                                  ! Position vector from inertial frame origin to tower-top / base plate (point O).
   RtHSdat%rV    = RtHSdat%rO  + RtHSdat%rOV                                                                                  ! Position vector from inertial frame origin to specified point on rotor-furl axis (point V)
   !RtHSdat%rP    = RtHSdat%rO  + RtHSdat%rOV + RtHSdat%rVP                                                                   ! Position vector from inertial frame origin to teeter pin (point P).
   RtHSdat%rP    = RtHSdat%rV  + RtHSdat%rVP                                                                                  ! Position vector from inertial frame origin to teeter pin (point P).
   RtHSdat%rQ    = RtHSdat%rP  + RtHSdat%rPQ                                                                                  ! Position vector from inertial frame origin to apex of rotation (point Q).
           rK    = RtHSdat%rO  + RtHSdat%rOW + RtHSdat%rWK                                                                    ! Position vector from inertial frame origin to tail fin center of pressure (point K).


   DO K = 1,p%NumBl ! Loop through all blades

      ! Calculate the position vector of the tip:
      RtHSdat%rS0S(:,K,p%TipNode) = ( p%TwistedSF(K,1,1,p%TipNode,0)*x%QT( DOF_BF(K,1) ) &                                       ! Position vector from the blade root (point S(0)) to the blade tip (point S(p%BldFlexL)).
                                    + p%TwistedSF(K,1,2,p%TipNode,0)*x%QT( DOF_BF(K,2) ) &
                                    + p%TwistedSF(K,1,3,p%TipNode,0)*x%QT( DOF_BE(K,1) )                     )*CoordSys%j1(K,:) &
                                  + ( p%TwistedSF(K,2,1,p%TipNode,0)*x%QT( DOF_BF(K,1) ) &
                                    + p%TwistedSF(K,2,2,p%TipNode,0)*x%QT( DOF_BF(K,2) ) &
                                    + p%TwistedSF(K,2,3,p%TipNode,0)*x%QT( DOF_BE(K,1) )                     )*CoordSys%j2(K,:) &
                                  + ( p%BldFlexL - 0.5* &
                                  (      p%AxRedBld(K,1,1,p%TipNode)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BF(K,1) ) &
                                    +    p%AxRedBld(K,2,2,p%TipNode)*x%QT( DOF_BF(K,2) )*x%QT( DOF_BF(K,2) ) &
                                    +    p%AxRedBld(K,3,3,p%TipNode)*x%QT( DOF_BE(K,1) )*x%QT( DOF_BE(K,1) ) &
                                    + 2.*p%AxRedBld(K,1,2,p%TipNode)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BF(K,2) ) &
                                    + 2.*p%AxRedBld(K,2,3,p%TipNode)*x%QT( DOF_BF(K,2) )*x%QT( DOF_BE(K,1) ) &
                                    + 2.*p%AxRedBld(K,1,3,p%TipNode)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BE(K,1) ) ) )*CoordSys%j3(K,:)
      RtHSdat%rQS (:,K,p%TipNode) = RtHSdat%rS0S(:,K,p%TipNode) + p%HubRad*CoordSys%j3(K,:)                                      ! Position vector from apex of rotation (point Q) to the blade tip (point S(p%BldFlexL)).
      RtHSdat%rS  (:,K,p%TipNode) = RtHSdat%rQS (:,K,p%TipNode) + RtHSdat%rQ                                                     ! Position vector from inertial frame origin      to the blade tip (point S(p%BldFlexL)).
      
      ! position vectors for blade root node:
      RtHSdat%rQS (:,K,0) = p%HubRad*CoordSys%j3(K,:)    
      RtHSdat%rS  (:,K,0) = p%HubRad*CoordSys%j3(K,:) + RtHSdat%rQ
      
      
         ! Calculate the position vector from the teeter pin to the blade root:
   
      RtHSdat%rPS0(:,K) = RtHSdat%rPQ + p%HubRad*CoordSys%j3(K,:)   ! Position vector from teeter pin (point P) to blade root (point S(0)).

      
      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements


      ! Calculate the position vector of the current node:

         RtHSdat%rS0S(:,K,J) = (  p%TwistedSF(K,1,1,J,0)*x%QT( DOF_BF(K,1) ) &                                                   ! Position vector from the blade root (point S(0)) to the current node (point S(RNodes(J)).
                                + p%TwistedSF(K,1,2,J,0)*x%QT( DOF_BF(K,2) ) &
                                + p%TwistedSF(K,1,3,J,0)*x%QT( DOF_BE(K,1) )                          )*CoordSys%j1(K,:) &
                            + (   p%TwistedSF(K,2,1,J,0)*x%QT( DOF_BF(K,1) ) &
                                + p%TwistedSF(K,2,2,J,0)*x%QT( DOF_BF(K,2) ) &
                                + p%TwistedSF(K,2,3,J,0)*x%QT( DOF_BE(K,1) )                          )*CoordSys%j2(K,:) &
                            + (  p%RNodes(J) - 0.5* &
                              (      p%AxRedBld(K,1,1,J)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BF(K,1) ) &
                               +     p%AxRedBld(K,2,2,J)*x%QT( DOF_BF(K,2) )*x%QT( DOF_BF(K,2) ) &
                               +     p%AxRedBld(K,3,3,J)*x%QT( DOF_BE(K,1) )*x%QT( DOF_BE(K,1) ) &
                               + 2.0*p%AxRedBld(K,1,2,J)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BF(K,2) ) &
                               + 2.0*p%AxRedBld(K,2,3,J)*x%QT( DOF_BF(K,2) )*x%QT( DOF_BE(K,1) ) &
                               + 2.0*p%AxRedBld(K,1,3,J)*x%QT( DOF_BF(K,1) )*x%QT( DOF_BE(K,1) )    ) )*CoordSys%j3(K,:)
         RtHSdat%rQS (:,K,J) = RtHSdat%rS0S(:,K,J) + p%HubRad*CoordSys%j3(K,:)                                                ! Position vector from apex of rotation (point Q) to the current node (point S(RNodes(J)).
         RtHSdat%rS  (:,K,J) = RtHSdat%rQS (:,K,J) + RtHSdat%rQ                                                               ! Position vector from inertial frame origin      to the current node (point S(RNodes(J)).


      END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements



   END DO !K = 1,p%NumBl
 

   

   !----------------------------------------------------------------------------------------------------
   ! Get the tower element positions
   !----------------------------------------------------------------------------------------------------
   RtHSdat%rZT (:,0) = RtHSdat%rZT0
   DO J = 1,p%TwrNodes  ! Loop through the tower nodes / elements


      ! Calculate the position vector of the current node:

      RtHSdat%rT0T(:,J) = ( p%TwrFASF(1,J,0)*x%QT(DOF_TFA1) + p%TwrFASF(2,J,0)*x%QT(DOF_TFA2)           )*CoordSys%a1 &       ! Position vector from base of flexible portion of tower (point T(0)) to current node (point T(J)).
                        + ( p%HNodes(J) - 0.5*(     p%AxRedTFA(1,1,J)*x%QT(DOF_TFA1)*x%QT(DOF_TFA1) &
                                              +     p%AxRedTFA(2,2,J)*x%QT(DOF_TFA2)*x%QT(DOF_TFA2) &
                                              + 2.0*p%AxRedTFA(1,2,J)*x%QT(DOF_TFA1)*x%QT(DOF_TFA2) &
                                              +     p%AxRedTSS(1,1,J)*x%QT(DOF_TSS1)*x%QT(DOF_TSS1) &
                                              +     p%AxRedTSS(2,2,J)*x%QT(DOF_TSS2)*x%QT(DOF_TSS2) &
                                              + 2.0*p%AxRedTSS(1,2,J)*x%QT(DOF_TSS1)*x%QT(DOF_TSS2)   ) )*CoordSys%a2 &
                        + ( p%TwrSSSF(1,J,0)*x%QT(DOF_TSS1) + p%TwrSSSF(2,J,0)*x%QT(DOF_TSS2)           )*CoordSys%a3
      RtHSdat%rZT (:,J) = RtHSdat%rZT0 + RtHSdat%rT0T(:,J)                                                                    ! Position vector from platform reference (point Z) to the current node (point T(HNodes(J)).


      RtHSdat%rT(:,J)   = RtHSdat%rZ   + RtHSdat%rZT (:,J)                                                                    ! Position vector from inertial frame origin        to the current node (point T(HNodes(J)).

   END DO


END SUBROUTINE CalculatePositions
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalculateAngularPosVelPAcc( p, x, CoordSys, RtHSdat )
! This routine is used to calculate the angular positions, velocities, and partial accelerations stored in other states that are used in
! both the CalcOutput and CalcContStateDeriv routines.
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    ! The coordinate systems that have been set for these states/time
   TYPE(ED_RtHndSide),           INTENT(INOUT)  :: RtHSdat     ! data from the RtHndSid module (contains positions to be set)

      !Local variables
!   REAL(ReKi)                   :: AngVelEN  (3)                                   ! Angular velocity of the nacelle (body N) in the inertia frame (body E for earth).
   REAL(ReKi)                   :: AngAccELt (3)                                   ! Portion of the angular acceleration of the low-speed shaft (body L) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades

   !-------------------------------------------------------------------------------------------------
   ! Angular and partial angular velocities
   !-------------------------------------------------------------------------------------------------

   ! Define the angular and partial angular velocities of all of the rigid bodies in the inertia frame:
   ! NOTE: PAngVelEN(I,D,:) = the Dth-derivative of the partial angular velocity of DOF I for body N in body E.

   RtHSdat%PAngVelEX(       :,0,:) = 0.0
   RtHSdat%PAngVelEX(DOF_R   ,0,:) =  CoordSys%z1
   RtHSdat%PAngVelEX(DOF_P   ,0,:) = -CoordSys%z3
   RtHSdat%PAngVelEX(DOF_Y   ,0,:) =  CoordSys%z2
   RtHSdat%AngVelEX                =                     x%QDT(DOF_R   )*RtHSdat%PAngVelEX(DOF_R   ,0,:) &
                                                       + x%QDT(DOF_P   )*RtHSdat%PAngVelEX(DOF_P   ,0,:) &
                                                       + x%QDT(DOF_Y   )*RtHSdat%PAngVelEX(DOF_Y   ,0,:)
   RtHSdat%AngPosEX                =                     x%QT (DOF_R   )*RtHSdat%PAngVelEX(DOF_R   ,0,:) &
                                                       + x%QT (DOF_P   )*RtHSdat%PAngVelEX(DOF_P   ,0,:) &
                                                       + x%QT (DOF_Y   )*RtHSdat%PAngVelEX(DOF_Y   ,0,:)

   RtHSdat%PAngVelEB(       :,0,:) =  RtHSdat%PAngVelEX(:,0,:)
   RtHSdat%PAngVelEB(DOF_TFA1,0,:) = -p%TwrFASF(1,p%TTopNode,1)*CoordSys%a3
   RtHSdat%PAngVelEB(DOF_TSS1,0,:) =  p%TwrSSSF(1,p%TTopNode,1)*CoordSys%a1
   RtHSdat%PAngVelEB(DOF_TFA2,0,:) = -p%TwrFASF(2,p%TTopNode,1)*CoordSys%a3
   RtHSdat%PAngVelEB(DOF_TSS2,0,:) =  p%TwrSSSF(2,p%TTopNode,1)*CoordSys%a1
   RtHSdat%AngVelEB                =  RtHSdat%AngVelEX + x%QDT(DOF_TFA1)*RtHSdat%PAngVelEB(DOF_TFA1,0,:) &
                                                       + x%QDT(DOF_TSS1)*RtHSdat%PAngVelEB(DOF_TSS1,0,:) &
                                                       + x%QDT(DOF_TFA2)*RtHSdat%PAngVelEB(DOF_TFA2,0,:) &
                                                       + x%QDT(DOF_TSS2)*RtHSdat%PAngVelEB(DOF_TSS2,0,:)
   RtHSdat%AngPosXB                =                     x%QT (DOF_TFA1)*RtHSdat%PAngVelEB(DOF_TFA1,0,:) &
                                                       + x%QT (DOF_TSS1)*RtHSdat%PAngVelEB(DOF_TSS1,0,:) &
                                                       + x%QT (DOF_TFA2)*RtHSdat%PAngVelEB(DOF_TFA2,0,:) &
                                                       + x%QT (DOF_TSS2)*RtHSdat%PAngVelEB(DOF_TSS2,0,:)

   RtHSdat%PAngVelEN(       :,0,:)= RtHSdat%PAngVelEB(:,0,:)
   RtHSdat%PAngVelEN(DOF_Yaw ,0,:)= CoordSys%d2
   RtHSdat%AngVelEN               = RtHSdat%AngVelEB + x%QDT(DOF_Yaw )*RtHSdat%PAngVelEN(DOF_Yaw ,0,:)

   RtHSdat%PAngVelER(       :,0,:)= RtHSdat%PAngVelEN(:,0,:)
   RtHSdat%PAngVelER(DOF_RFrl,0,:)= CoordSys%rfa
   RtHSdat%AngVelER               = RtHSdat%AngVelEN + x%QDT(DOF_RFrl)*RtHSdat%PAngVelER(DOF_RFrl,0,:)

   RtHSdat%PAngVelEL(       :,0,:)= RtHSdat%PAngVelER(:,0,:)
   RtHSdat%PAngVelEL(DOF_GeAz,0,:)= CoordSys%c1
   RtHSdat%PAngVelEL(DOF_DrTr,0,:)= CoordSys%c1
   RtHSdat%AngVelEL               = RtHSdat%AngVelER + x%QDT(DOF_GeAz)*RtHSdat%PAngVelEL(DOF_GeAz,0,:) &
                                                           + x%QDT(DOF_DrTr)*RtHSdat%PAngVelEL(DOF_DrTr,0,:)

   RtHSdat%PAngVelEH(       :,0,:)= RtHSdat%PAngVelEL(:,0,:)
   RtHSdat%AngVelEH               = RtHSdat%AngVelEL
IF ( p%NumBl == 2 )  THEN ! 2-blader
   RtHSdat%PAngVelEH(DOF_Teet,0,:)= CoordSys%f2
   RtHSdat%AngVelEH               = RtHSdat%AngVelEH + x%QDT(DOF_Teet)*RtHSdat%PAngVelEH(DOF_Teet,0,:)
ENDIF

   RtHSdat%PAngVelEG(       :,0,:) = RtHSdat%PAngVelER(:,0,:)
   RtHSdat%PAngVelEG(DOF_GeAz,0,:) = p%GBRatio*CoordSys%c1
   RtHSdat%AngVelEG                = RtHSdat%AngVelER + x%QDT(DOF_GeAz)*RtHSdat%PAngVelEG(DOF_GeAz,0,:)

   RtHSdat%PAngVelEA(       :,0,:) = RtHSdat%PAngVelEN(:,0,:)
   RtHSdat%PAngVelEA(DOF_TFrl,0,:) = CoordSys%tfa
   RtHSdat%AngVelEA                = RtHSdat%AngVelEN + x%QDT(DOF_TFrl)*RtHSdat%PAngVelEA(DOF_TFrl,0,:)



   ! Define the 1st derivatives of the partial angular velocities of all
   !   of the rigid bodies in the inertia frame and the portion of the angular
   !   acceleration of the rigid bodies in the inertia frame associated with
   !   everything but the QD2T()'s:

   RtHSdat%PAngVelEX(       :,1,:) = 0.0
   RtHSdat%AngAccEXt               = 0.0

   RtHSdat%PAngVelEB(       :,1,:) =                  RtHSdat%PAngVelEX(:,1,:)
   RtHSdat%PAngVelEB(DOF_TFA1,1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEX,                   RtHSdat%PAngVelEB(DOF_TFA1,0,:) )
   RtHSdat%PAngVelEB(DOF_TSS1,1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEX,                   RtHSdat%PAngVelEB(DOF_TSS1,0,:) )
   RtHSdat%PAngVelEB(DOF_TFA2,1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEX,                   RtHSdat%PAngVelEB(DOF_TFA2,0,:) )
   RtHSdat%PAngVelEB(DOF_TSS2,1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEX,                   RtHSdat%PAngVelEB(DOF_TSS2,0,:) )
   RtHSdat%AngAccEBt               =                  RtHSdat%AngAccEXt + x%QDT(DOF_TFA1)*RtHSdat%PAngVelEB(DOF_TFA1,1,:) &
                                                                        + x%QDT(DOF_TSS1)*RtHSdat%PAngVelEB(DOF_TSS1,1,:) &
                                                                        + x%QDT(DOF_TFA2)*RtHSdat%PAngVelEB(DOF_TFA2,1,:) &
                                                                        + x%QDT(DOF_TSS2)*RtHSdat%PAngVelEB(DOF_TSS2,1,:)

   RtHSdat%PAngVelEN(       :,1,:) =                 RtHSdat%PAngVelEB(:,1,:)
   RtHSdat%PAngVelEN(DOF_Yaw ,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEB,                    RtHSdat%PAngVelEN(DOF_Yaw ,0,:) )
   RtHSdat%AngAccENt               =                 RtHSdat%AngAccEBt  + x%QDT(DOF_Yaw )*RtHSdat%PAngVelEN(DOF_Yaw ,1,:)

   RtHSdat%PAngVelER(       :,1,:) =                 RtHSdat%PAngVelEN(:,1,:)
   RtHSdat%PAngVelER(DOF_RFrl,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEN,                    RtHSdat%PAngVelER(DOF_RFrl,0,:) )
   RtHSdat%AngAccERt               =                 RtHSdat%AngAccENt  + x%QDT(DOF_RFrl)*RtHSdat%PAngVelER(DOF_RFrl,1,:)

   RtHSdat%PAngVelEL(       :,1,:) =                 RtHSdat%PAngVelER(:,1,:)
   RtHSdat%PAngVelEL(DOF_GeAz,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelER,                    RtHSdat%PAngVelEL(DOF_GeAz,0,:) )
   RtHSdat%PAngVelEL(DOF_DrTr,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelER,                    RtHSdat%PAngVelEL(DOF_DrTr,0,:) )
           AngAccELt               =                 RtHSdat%AngAccERt  + x%QDT(DOF_GeAz)*RtHSdat%PAngVelEL(DOF_GeAz,1,:) &
                                                                        + x%QDT(DOF_DrTr)*RtHSdat%PAngVelEL(DOF_DrTr,1,:)

   RtHSdat%PAngVelEH(       :,1,:) = RtHSdat%PAngVelEL(:,1,:)
   RtHSdat%AngAccEHt               =                  AngAccELt
IF ( p%NumBl == 2 )  THEN ! 2-blader
   RtHSdat%PAngVelEH(DOF_Teet,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEH,                    RtHSdat%PAngVelEH(DOF_Teet,0,:) )
   RtHSdat%AngAccEHt               =                 RtHSdat%AngAccEHt   + x%QDT(DOF_Teet)*RtHSdat%PAngVelEH(DOF_Teet,1,:)
ENDIF

   RtHSdat%PAngVelEG(       :,1,:) = RtHSdat%PAngVelER(:,1,:)
   RtHSdat%PAngVelEG(DOF_GeAz,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelER,                    RtHSdat%PAngVelEG(DOF_GeAz,0,:) )
   RtHSdat%AngAccEGt               =                 RtHSdat%AngAccERt  + x%QDT(DOF_GeAz)*RtHSdat%PAngVelEG(DOF_GeAz,1,:)

   RtHSdat%PAngVelEA(       :,1,:) = RtHSdat%PAngVelEN(:,1,:)
   RtHSdat%PAngVelEA(DOF_TFrl,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEN,                    RtHSdat%PAngVelEA(DOF_TFrl,0,:) )
   RtHSdat%AngAccEAt               =                 RtHSdat%AngAccENt  + x%QDT(DOF_TFrl)*RtHSdat%PAngVelEA(DOF_TFrl,1,:)



   DO K = 1,p%NumBl ! Loop through all blades

      ! Define the partial angular velocities of the tip (body M(p%BldFlexL)) in the  inertia frame:
      ! NOTE: PAngVelEM(K,J,I,D,:) = the Dth-derivative of the partial angular velocity of DOF I for body M of blade K, element J in body E.

      RtHSdat%PAngVelEM(K,p%TipNode,          :,0,:) = RtHSdat%PAngVelEH(:,0,:)
      RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,1),0,:) = - p%TwistedSF(K,2,1,p%TipNode,1)*CoordSys%j1(K,:) &
                                                       + p%TwistedSF(K,1,1,p%TipNode,1)*CoordSys%j2(K,:)
      RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,2),0,:) = - p%TwistedSF(K,2,2,p%TipNode,1)*CoordSys%j1(K,:) &
                                                       + p%TwistedSF(K,1,2,p%TipNode,1)*CoordSys%j2(K,:)
      RtHSdat%PAngVelEM(K,p%TipNode,DOF_BE(K,1),0,:) = - p%TwistedSF(K,2,3,p%TipNode,1)*CoordSys%j1(K,:) &
                                                       + p%TwistedSF(K,1,3,p%TipNode,1)*CoordSys%j2(K,:)
   !           AngVelHM(K,p%TipNode              ,:) =  RtHSdat%AngVelEH + x%QDT(DOF_BF(K,1))*RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,1),0,:) & ! Currently
   !                                                                     + x%QDT(DOF_BF(K,2))*RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,2),0,:) & ! unused
   !                                                                     + x%QDT(DOF_BE(K,1))*RtHSdat%PAngVelEM(K,p%TipNode,DOF_BE(K,1),0,:)   ! calculations
       RtHSdat%AngPosHM(:,K,p%TipNode) =        x%QT (DOF_BF(K,1))*RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,1),0,:) &
                                              + x%QT (DOF_BF(K,2))*RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,2),0,:) &
                                              + x%QT (DOF_BE(K,1))*RtHSdat%PAngVelEM(K,p%TipNode,DOF_BE(K,1),0,:)


      ! Define the 1st derivatives of the partial angular velocities of the tip
      !   (body M(p%BldFlexL)) in the inertia frame:

   ! NOTE: These are currently unused by the code, therefore, they need not
   !       be calculated.  Thus, they are currently commented out.  If it
   !       turns out that they are ever needed (i.e., if inertias of the
   !       blade elements are ever added, etc...) simply uncomment out these
   !       computations:
   !   RtHSdat%PAngVelEM(K,p%TipNode,          :,1,:) = RtHSdat%PAngVelEH(:,1,:)
   !   RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,1),1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEH, RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,1),0,:)    )
   !   RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,2),1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEH, RtHSdat%PAngVelEM(K,p%TipNode,DOF_BF(K,2),0,:)    )
   !   RtHSdat%PAngVelEM(K,p%TipNode,DOF_BE(K,1),1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEH, RtHSdat%PAngVelEM(K,p%TipNode,DOF_BE(K,1),0,:)    )


      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements
      ! Define the partial angular velocities of the current node (body M(RNodes(J))) in the inertia frame:
      ! NOTE: PAngVelEM(K,J,I,D,:) = the Dth-derivative of the partial angular velocity
      !   of DOF I for body M of blade K, element J in body E.

         RtHSdat%PAngVelEM(K,J,          :,0,:) = RtHSdat%PAngVelEH(:,0,:)
         RtHSdat%PAngVelEM(K,J,DOF_BF(K,1),0,:) = - p%TwistedSF(K,2,1,J,1)*CoordSys%j1(K,:) &
                                                + p%TwistedSF(K,1,1,J,1)*CoordSys%j2(K,:)
         RtHSdat%PAngVelEM(K,J,DOF_BF(K,2),0,:) = - p%TwistedSF(K,2,2,J,1)*CoordSys%j1(K,:) &
                                                + p%TwistedSF(K,1,2,J,1)*CoordSys%j2(K,:)
         RtHSdat%PAngVelEM(K,J,DOF_BE(K,1),0,:) = - p%TwistedSF(K,2,3,J,1)*CoordSys%j1(K,:) &
                                                + p%TwistedSF(K,1,3,J,1)*CoordSys%j2(K,:)
!                 AngVelHM(K,J              ,:) =  RtHSdat%AngVelEH + x%QDT(DOF_BF(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,1),0,:) & ! Currently
!                                                                   + x%QDT(DOF_BF(K,2))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,2),0,:) & ! unused
!                                                                   + x%QDT(DOF_BE(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BE(K,1),0,:)   ! calculations
          RtHSdat%AngPosHM(:,K,J              ) =             x%QT (DOF_BF(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,1),0,:) &
                                                    + x%QT (DOF_BF(K,2))*RtHSdat%PAngVelEM(K,J,DOF_BF(K,2),0,:) &
                                                    + x%QT (DOF_BE(K,1))*RtHSdat%PAngVelEM(K,J,DOF_BE(K,1),0,:)


      ! Define the 1st derivatives of the partial angular velocities of the current node (body M(RNodes(J))) in the inertia frame:

   ! NOTE: These are currently unused by the code, therefore, they need not
   !       be calculated.  Thus, they are currently commented out.  If it
   !       turns out that they are ever needed (i.e., if inertias of the
   !       blade elements are ever added, etc...) simply uncomment out these computations:
   !      RtHSdat%PAngVelEM(K,J,          :,1,:) = RtHSdat%PAngVelEH(:,1,:)
   !      RtHSdat%PAngVelEM(K,J,DOF_BF(K,1),1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEH, PAngVelEM(K,J,DOF_BF(K,1),0,:) )
   !      RtHSdat%PAngVelEM(K,J,DOF_BF(K,2),1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEH, PAngVelEM(K,J,DOF_BF(K,2),0,:) )
   !      RtHSdat%PAngVelEM(K,J,DOF_BE(K,1),1,:) = CROSS_PRODUCT(   RtHSdat%AngVelEH, PAngVelEM(K,J,DOF_BE(K,1),0,:) )


      END DO !J = 1,p%BldNodes ! Loop through the blade nodes / elements

   END DO !K = 1,p%NumBl


   !...............
   ! tower values:
   !...............

   DO J = 0,p%TwrNodes  ! Loop through the tower nodes / elements

      ! Define the partial angular velocities (and their 1st derivatives) of the
      !   current node (body F(HNodes(J))  in the inertia frame.
      ! Also define the overall angular velocity of the current node in the inertia frame.
      !   Also, define the portion of the angular acceleration of the current node
      !   in the inertia frame associated with everything but the QD2T()'s:

      ! NOTE: PAngVelEF(J,I,D,:) = the Dth-derivative of the partial angular velocity
      !   of DOF I for body F of element J in body E.

      RtHSdat%PAngVelEF (J,       :,0,:) = RtHSdat%PAngVelEX(:,0,:)
      RtHSdat%PAngVelEF (J,DOF_TFA1,0,:) = -p%TwrFASF(1,J,1)*CoordSys%a3
      RtHSdat%PAngVelEF (J,DOF_TSS1,0,:) =  p%TwrSSSF(1,J,1)*CoordSys%a1
      RtHSdat%PAngVelEF (J,DOF_TFA2,0,:) = -p%TwrFASF(2,J,1)*CoordSys%a3
      RtHSdat%PAngVelEF (J,DOF_TSS2,0,:) =  p%TwrSSSF(2,J,1)*CoordSys%a1

      RtHSdat%PAngVelEF (J,       :,1,:) = RtHSdat%PAngVelEX(:,1,:)
      RtHSdat%PAngVelEF (J,DOF_TFA1,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEX  ,  RtHSdat%PAngVelEF(J,DOF_TFA1,0,:) )
      RtHSdat%PAngVelEF (J,DOF_TSS1,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEX  ,  RtHSdat%PAngVelEF(J,DOF_TSS1,0,:) )
      RtHSdat%PAngVelEF (J,DOF_TFA2,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEX  ,  RtHSdat%PAngVelEF(J,DOF_TFA2,0,:) )
      RtHSdat%PAngVelEF (J,DOF_TSS2,1,:) = CROSS_PRODUCT(  RtHSdat%AngVelEX  ,  RtHSdat%PAngVelEF(J,DOF_TSS2,0,:) )


      RtHSdat%AngVelEF (:,J)            =  RtHSdat%AngVelEX  + x%QDT(DOF_TFA1)*RtHSdat%PAngVelEF(J,DOF_TFA1,0,:) &
                                                             + x%QDT(DOF_TSS1)*RtHSdat%PAngVelEF(J,DOF_TSS1,0,:) &
                                                             + x%QDT(DOF_TFA2)*RtHSdat%PAngVelEF(J,DOF_TFA2,0,:) &
                                                             + x%QDT(DOF_TSS2)*RtHSdat%PAngVelEF(J,DOF_TSS2,0,:)

      RtHSdat%AngPosXF (:,J)            =                      x%QT (DOF_TFA1)*RtHSdat%PAngVelEF(J,DOF_TFA1,0,:) &
                                                             + x%QT (DOF_TSS1)*RtHSdat%PAngVelEF(J,DOF_TSS1,0,:) &
                                                             + x%QT (DOF_TFA2)*RtHSdat%PAngVelEF(J,DOF_TFA2,0,:) &
                                                             + x%QT (DOF_TSS2)*RtHSdat%PAngVelEF(J,DOF_TSS2,0,:)
      RtHSdat%AngPosEF (:,J)            =  RtHSdat%AngPosEX  + RtHSdat%AngPosXF(:,J)
      RtHSdat%AngAccEFt(:,J)            =  RtHSdat%AngAccEXt + x%QDT(DOF_TFA1)*RtHSdat%PAngVelEF(J,DOF_TFA1,1,:) &
                                                             + x%QDT(DOF_TSS1)*RtHSdat%PAngVelEF(J,DOF_TSS1,1,:) &
                                                             + x%QDT(DOF_TFA2)*RtHSdat%PAngVelEF(J,DOF_TFA2,1,:) &
                                                             + x%QDT(DOF_TSS2)*RtHSdat%PAngVelEF(J,DOF_TSS2,1,:)

   END DO ! J


END SUBROUTINE CalculateAngularPosVelPAcc
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalculateLinearVelPAcc( p, x, CoordSys, RtHSdat )
! This routine is used to calculate the linear velocities and accelerations stored in other states that are used in
! both the CalcOutput and CalcContStateDeriv routines.
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    ! The coordinate systems that have been set for these states/time
   TYPE(ED_RtHndSide),           INTENT(INOUT)  :: RtHSdat     ! data from the RtHndSid module (contains positions to be set)

      ! Local variables
   REAL(ReKi)                   :: LinAccEKt (3)                                   ! "Portion of the linear acceleration of the tail fin  center of pressure (point K) in the inertia frame (body E for earth) associated with everything but the QD2T()'s"
   REAL(ReKi)                   :: LinAccEPt (3)                                   ! "Portion of the linear acceleration of the teeter pin (point P) in the inertia frame (body E for earth) associated with everything but the QD2T()'s"
   REAL(ReKi)                   :: LinAccEQt (3)                                   ! "Portion of the linear acceleration of the apex of rotation (point Q) in the inertia frame (body E for earth) associated with everything but the QD2T()'s"
   REAL(ReKi)                   :: LinAccEVt (3)                                   ! "Portion of the linear acceleration of the selected point on the rotor-furl axis (point V) in the inertia frame (body E for earth) associated with everything but the QD2T()'s"
   REAL(ReKi)                   :: LinAccEWt (3)                                   ! "Portion of the linear acceleration of the selected point on the  tail-furl axis (point W) in the inertia frame (body E for earth) associated with everything but the QD2T()'s"
   REAL(ReKi)                   :: LinVelEK  (3)                                   ! "Linear velocity of tail fin center-of-pressure (point K) in the inertia frame"
   REAL(ReKi)                   :: LinVelHS  (3)                                   ! "Relative linear velocity of the current point on the current blade (point S) in the hub frame (body H)"
   REAL(ReKi)                   :: LinVelXO  (3)                                   ! "Relative linear velocity of the tower-top / base plate (point O) in the platform (body X)"
   REAL(ReKi)                   :: LinVelXT  (3)                                   ! "Relative linear velocity of the current point on the tower (point T) in the platform (body X)"

   REAL(ReKi)                   :: EwAXrWI   (3)                                   ! = AngVelEA X rWI
   REAL(ReKi)                   :: EwAXrWJ   (3)                                   ! = AngVelEA X rWJ
   REAL(ReKi)                   :: EwAXrWK   (3)                                   ! = AngVelEA X rWK
   REAL(ReKi)                   :: EwHXrPQ   (3)                                   ! = AngVelEH X rPQ
   REAL(ReKi)                   :: EwHXrQC   (3)                                   ! = AngVelEH X rQC
   REAL(ReKi)                   :: EwHXrQS   (3)                                   ! = AngVelEH X rQS of the current blade point S.
   REAL(ReKi)                   :: EwNXrOU   (3)                                   ! = AngVelEN X rOU
   REAL(ReKi)                   :: EwNXrOV   (3)                                   ! = AngVelEN X rOV
   REAL(ReKi)                   :: EwNXrOW   (3)                                   ! = AngVelEN X rOW
   REAL(ReKi)                   :: EwRXrVD   (3)                                   ! = AngVelER X rVD
   REAL(ReKi)                   :: EwRXrVIMU (3)                                   ! = AngVelER X rVIMU
   REAL(ReKi)                   :: EwRXrVP   (3)                                   ! = AngVelER X rVP
   REAL(ReKi)                   :: EwXXrZO   (3)                                   ! = AngVelEX X rZO
   REAL(ReKi)                   :: EwXXrZT   (3)                                   ! = AngVelEX X rZT
   REAL(ReKi)                   :: EwXXrZY   (3)                                   ! = AngVelEX X rZY

   REAL(ReKi)                   :: TmpVec0   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec1   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec2   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec3   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec4   (3)                                   ! A temporary vector used in various computations.

   INTEGER(IntKi)               :: I                                               ! Loops through some or all of the DOFs
   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades


      ! Initializations:

   RtHSdat%LinAccECt   = 0.0
   RtHSdat%LinAccEDt   = 0.0
   RtHSdat%LinAccEIMUt = 0.0
   RtHSdat%LinAccEIt   = 0.0
   RtHSdat%LinAccEJt   = 0.0
           LinAccEKt   = 0.0
   RtHSdat%LinAccEOt   = 0.0
           LinAccEPt   = 0.0
           LinAccEQt   = 0.0
   RtHSdat%LinAccESt   = 0.0
   RtHSdat%LinAccETt   = 0.0
   RtHSdat%LinAccEUt   = 0.0
           LinAccEVt   = 0.0
           LinAccEWt   = 0.0
   RtHSdat%LinAccEYt   = 0.0
   RtHSdat%LinAccEZt   = 0.0


   !-------------------------------------------------------------------------------------------------
   ! Partial linear velocities and accelerations
   !-------------------------------------------------------------------------------------------------

      ! Define the partial linear velocities (and their 1st derivatives) of all of
      !   the points on the wind turbine in the inertia frame that are not
      !   dependent on the distributed tower or blade parameters.  Also, define
      !   the portion of the linear acceleration of the points in the inertia
      !   frame associated with everything but the QD2T()'s:
      ! NOTE: PLinVelEX(I,D,:) = the Dth-derivative of the partial linear velocity
      !   of DOF I for point X in body E.

   EwXXrZY   = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%rZY   ) !
   EwXXrZO   = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%rZO   ) !
   EwNXrOU   = CROSS_PRODUCT( RtHSdat%AngVelEN, RtHSdat%rOU   ) !
   EwNXrOV   = CROSS_PRODUCT( RtHSdat%AngVelEN, RtHSdat%rOV   ) !
   EwRXrVD   = CROSS_PRODUCT( RtHSdat%AngVelER, RtHSdat%rVD   ) ! Cross products
   EwRXrVIMU = CROSS_PRODUCT( RtHSdat%AngVelER, RtHSdat%rVIMU ) ! that are used
   EwRXrVP   = CROSS_PRODUCT( RtHSdat%AngVelER, RtHSdat%rVP   ) ! in the following
   EwHXrPQ   = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%rPQ   ) ! DO...LOOPs
   EwHXrQC   = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%rQC   ) !
   EwNXrOW   = CROSS_PRODUCT( RtHSdat%AngVelEN, RtHSdat%rOW   ) !
   EwAXrWI   = CROSS_PRODUCT( RtHSdat%AngVelEA, RtHSdat%rWI   ) !
   EwAXrWJ   = CROSS_PRODUCT( RtHSdat%AngVelEA, RtHSdat%rWJ   ) !
   EwAXrWK   = CROSS_PRODUCT( RtHSdat%AngVelEA, RtHSdat%rWK   ) !


   RtHSdat%PLinVelEZ(       :,:,:) = 0.0
   RtHSdat%PLinVelEZ(DOF_Sg  ,0,:) =  CoordSys%z1
   RtHSdat%PLinVelEZ(DOF_Sw  ,0,:) = -CoordSys%z3
   RtHSdat%PLinVelEZ(DOF_Hv  ,0,:) =  CoordSys%z2

   RtHSdat%LinVelEZ                =   x%QDT(DOF_Sg  )*RtHSdat%PLinVelEZ(DOF_Sg  ,0,:) &
                                     + x%QDT(DOF_Sw  )*RtHSdat%PLinVelEZ(DOF_Sw  ,0,:) &
                                     + x%QDT(DOF_Hv  )*RtHSdat%PLinVelEZ(DOF_Hv  ,0,:)


   RtHSdat%PLinVelEY(       :,:,:) = RtHSdat%PLinVelEZ(:,:,:)
   DO I = 1,NPX   ! Loop through all DOFs associated with the angular motion of the platform (body X)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I)   ,0,:), RtHSdat%rZY  )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I)   ,0,:),     EwXXrZY  )

      RtHSdat%PLinVelEY(PX(I),0,:) = TmpVec0   +                       RtHSdat%PLinVelEY(PX(I)   ,0,:)
      RtHSdat%PLinVelEY(PX(I),1,:) = TmpVec1   +                       RtHSdat%PLinVelEY(PX(I)   ,1,:)

       RtHSdat%LinAccEYt           = RtHSdat%LinAccEYt + x%QDT(PX(I) )*RtHSdat%PLinVelEY(PX(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the platform (body X)


   RtHSdat%PLinVelEO(       :,:,:) = RtHSdat%PLinVelEZ(:,:,:)
   RtHSdat%PLinVelEO(DOF_TFA1,0,:) = CoordSys%a1 - (   p%AxRedTFA(1,1,p%TTopNode)* x%QT(DOF_TFA1) &
                                                     + p%AxRedTFA(1,2,p%TTopNode)* x%QT(DOF_TFA2)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TSS1,0,:) = CoordSys%a3 - (   p%AxRedTSS(1,1,p%TTopNode)* x%QT(DOF_TSS1) &
                                                     + p%AxRedTSS(1,2,p%TTopNode)* x%QT(DOF_TSS2)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TFA2,0,:) = CoordSys%a1 - (   p%AxRedTFA(2,2,p%TTopNode)* x%QT(DOF_TFA2) &
                                                     + p%AxRedTFA(1,2,p%TTopNode)* x%QT(DOF_TFA1)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TSS2,0,:) = CoordSys%a3 - (   p%AxRedTSS(2,2,p%TTopNode)* x%QT(DOF_TSS2) &
                                                     + p%AxRedTSS(1,2,p%TTopNode)* x%QT(DOF_TSS1)   )*CoordSys%a2

   TmpVec1 = CROSS_PRODUCT(   RtHSdat%AngVelEX   , RtHSdat%PLinVelEO(DOF_TFA1,0,:) )
   TmpVec2 = CROSS_PRODUCT(   RtHSdat%AngVelEX   , RtHSdat%PLinVelEO(DOF_TSS1,0,:) )
   TmpVec3 = CROSS_PRODUCT(   RtHSdat%AngVelEX   , RtHSdat%PLinVelEO(DOF_TFA2,0,:) )
   TmpVec4 = CROSS_PRODUCT(   RtHSdat%AngVelEX   , RtHSdat%PLinVelEO(DOF_TSS2,0,:) )

   RtHSdat%PLinVelEO(DOF_TFA1,1,:) = TmpVec1 - (   p%AxRedTFA(1,1,p%TTopNode)*x%QDT(DOF_TFA1) &
                                                 + p%AxRedTFA(1,2,p%TTopNode)*x%QDT(DOF_TFA2)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TSS1,1,:) = TmpVec2 - (   p%AxRedTSS(1,1,p%TTopNode)*x%QDT(DOF_TSS1) &
                                                 + p%AxRedTSS(1,2,p%TTopNode)*x%QDT(DOF_TSS2)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TFA2,1,:) = TmpVec3 - (   p%AxRedTFA(2,2,p%TTopNode)*x%QDT(DOF_TFA2) &
                                                 + p%AxRedTFA(1,2,p%TTopNode)*x%QDT(DOF_TFA1)   )*CoordSys%a2
   RtHSdat%PLinVelEO(DOF_TSS2,1,:) = TmpVec4 - (   p%AxRedTSS(2,2,p%TTopNode)*x%QDT(DOF_TSS2) &
                                                 + p%AxRedTSS(1,2,p%TTopNode)*x%QDT(DOF_TSS1)   )*CoordSys%a2

    LinVelXO               =              x%QDT(DOF_TFA1)*RtHSdat%PLinVelEO(DOF_TFA1,0,:) &
                                        + x%QDT(DOF_TSS1)*RtHSdat%PLinVelEO(DOF_TSS1,0,:) &
                                        + x%QDT(DOF_TFA2)*RtHSdat%PLinVelEO(DOF_TFA2,0,:) &
                                        + x%QDT(DOF_TSS2)*RtHSdat%PLinVelEO(DOF_TSS2,0,:)
    RtHSdat%LinAccEOt              =      x%QDT(DOF_TFA1)*RtHSdat%PLinVelEO(DOF_TFA1,1,:) &
                                        + x%QDT(DOF_TSS1)*RtHSdat%PLinVelEO(DOF_TSS1,1,:) &
                                        + x%QDT(DOF_TFA2)*RtHSdat%PLinVelEO(DOF_TFA2,1,:) &
                                        + x%QDT(DOF_TSS2)*RtHSdat%PLinVelEO(DOF_TSS2,1,:)
    
   RtHSdat%LinVelEO = LinVelXO + RtHSdat%LinVelEZ
   DO I = 1,NPX   ! Loop through all DOFs associated with the angular motion of the platform (body X)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I)   ,0,:), RtHSdat%rZO                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I)   ,0,:),     EwXXrZO + LinVelXO      )

      RtHSdat%PLinVelEO(PX(I),0,:) = TmpVec0    +                       RtHSdat%PLinVelEO(PX(I)   ,0,:)
      RtHSdat%PLinVelEO(PX(I),1,:) = TmpVec1    +                       RtHSdat%PLinVelEO(PX(I)   ,1,:)

      RtHSdat%LinVelEO             =  RtHSdat%LinVelEO  + x%QDT(PX(I) )*RtHSdat%PLinVelEO(PX(I)   ,0,:)
      RtHSdat%LinAccEOt            =  RtHSdat%LinAccEOt + x%QDT(PX(I) )*RtHSdat%PLinVelEO(PX(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the platform (body X)
                     

   RtHSdat%PLinVelEU(       :,:,:) = RtHSdat%PLinVelEO(:,:,:)
   DO I = 1,NPN   ! Loop through all DOFs associated with the angular motion of the nacelle (body N)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:), RtHSdat%rOU                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:),     EwNXrOU                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,1,:), RtHSdat%rOU                 )

      RtHSdat%PLinVelEU(PN(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEU(PN(I)   ,0,:)
      RtHSdat%PLinVelEU(PN(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEU(PN(I)   ,1,:)

       RtHSdat%LinAccEUt           =  RtHSdat%LinAccEUt + x%QDT(PN(I) )*RtHSdat%PLinVelEU(PN(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the nacelle (body N)


   RtHSdat%PLinVelEV(       :,:,:) = RtHSdat%PLinVelEO(:,:,:)
   DO I = 1,NPN   ! Loop through all DOFs associated with the angular motion of the nacelle (body N)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:), RtHSdat%rOV                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:),     EwNXrOV                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,1,:), RtHSdat%rOV                 )

      RtHSdat%PLinVelEV(PN(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEV(PN(I)   ,0,:)
      RtHSdat%PLinVelEV(PN(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEV(PN(I)   ,1,:)

       LinAccEVt                   =  LinAccEVt + x%QDT(PN(I) )*RtHSdat%PLinVelEV(PN(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the nacelle (body N)


   RtHSdat%PLinVelED(       :,:,:) = RtHSdat%PLinVelEV(:,:,:)
   DO I = 1,NPR   ! Loop through all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,0,:), RtHSdat%rVD                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,0,:),     EwRXrVD                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,1,:), RtHSdat%rVD                 )

      RtHSdat%PLinVelED(PR(I),0,:) = TmpVec0    +                       RtHSdat%PLinVelED(PR(I)   ,0,:)
      RtHSdat%PLinVelED(PR(I),1,:) = TmpVec1    + TmpVec2 +             RtHSdat%PLinVelED(PR(I)   ,1,:)

      RtHSdat%LinAccEDt            =  RtHSdat%LinAccEDt + x%QDT(PR(I) )*RtHSdat%PLinVelED(PR(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)


   RtHSdat%PLinVelEIMU(     :,:,:) = RtHSdat%PLinVelEV(:,:,:)
    RtHSdat%LinVelEIMU             =  RtHSdat%LinVelEZ
   DO I = 1,NPR   ! Loop through all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,0,:), RtHSdat%rVIMU               )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,0,:),     EwRXrVIMU               )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelER(PR(I)   ,1,:), RtHSdat%rVIMU               )

      RtHSdat%PLinVelEIMU(PR(I),0,:) = TmpVec0    +                         RtHSdat%PLinVelEIMU(PR(I) ,0,:)
      RtHSdat%PLinVelEIMU(PR(I),1,:) = TmpVec1    + TmpVec2 +               RtHSdat%PLinVelEIMU(PR(I) ,1,:)

      RtHSdat%LinVelEIMU             =  RtHSdat%LinVelEIMU  + x%QDT(PR(I) )*RtHSdat%PLinVelEIMU(PR(I) ,0,:)
      RtHSdat%LinAccEIMUt            =  RtHSdat%LinAccEIMUt + x%QDT(PR(I) )*RtHSdat%PLinVelEIMU(PR(I) ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)


   RtHSdat%PLinVelEP(       :,:,:) = RtHSdat%PLinVelEV(:,:,:)
   DO I = 1,NPR   ! Loop through all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)

      TmpVec0 = CROSS_PRODUCT(             RtHSdat%PAngVelER(PR(I)   ,0,:),     RtHSdat%rVP                 )
      TmpVec1 = CROSS_PRODUCT(             RtHSdat%PAngVelER(PR(I)   ,0,:), EwRXrVP                 )
      TmpVec2 = CROSS_PRODUCT(             RtHSdat%PAngVelER(PR(I)   ,1,:),     RtHSdat%rVP                 )

      RtHSdat%PLinVelEP(PR(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEP(PR(I)   ,0,:)
      RtHSdat%PLinVelEP(PR(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEP(PR(I)   ,1,:)

       LinAccEPt           =  LinAccEPt + x%QDT(PR(I) )*RtHSdat%PLinVelEP(PR(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the structure that furls with the rotor (not including rotor) (body R)


   RtHSdat%PLinVelEQ(       :,:,:) = RtHSdat%PLinVelEP(:,:,:)
    RtHSdat%LinVelEQ               =  RtHSdat%LinVelEZ
   DO I = 1,p%NPH   ! Loop through all DOFs associated with the angular motion of the hub (body H)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,0,:),   RtHSdat%rPQ  )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,0,:),       EwHXrPQ  )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,1,:),   RtHSdat%rPQ  )

      RtHSdat%PLinVelEQ(p%PH(I),0,:) = TmpVec0    +                 RtHSdat%PLinVelEQ(p%PH(I)   ,0,:)
      RtHSdat%PLinVelEQ(p%PH(I),1,:) = TmpVec1    + TmpVec2 +       RtHSdat%PLinVelEQ(p%PH(I)   ,1,:)

      RtHSdat%LinVelEQ               =  RtHSdat%LinVelEQ  + x%QDT(p%PH(I) )*RtHSdat%PLinVelEQ(p%PH(I)   ,0,:)
      LinAccEQt                      =          LinAccEQt + x%QDT(p%PH(I) )*RtHSdat%PLinVelEQ(p%PH(I)   ,1,:)
      
   ENDDO          ! I - all DOFs associated with the angular motion of the hub (body H)


   RtHSdat%PLinVelEC(       :,:,:) = RtHSdat%PLinVelEQ(:,:,:)
   DO I = 1,p%NPH   ! Loop through all DOFs associated with the angular motion of the hub (body H)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,0,:), RtHSdat%rQC )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,0,:),     EwHXrQC )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEH(p%PH(I)   ,1,:), RtHSdat%rQC )

      RtHSdat%PLinVelEC(p%PH(I),0,:) = TmpVec0    +                         RtHSdat%PLinVelEC(p%PH(I)   ,0,:)
      RtHSdat%PLinVelEC(p%PH(I),1,:) = TmpVec1    + TmpVec2 +               RtHSdat%PLinVelEC(p%PH(I)   ,1,:)

      RtHSdat%LinAccECt              =  RtHSdat%LinAccECt + x%QDT(p%PH(I) )*RtHSdat%PLinVelEC(p%PH(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the hub (body H)




   DO K = 1,p%NumBl ! Loop through all blades

      DO J = 0,p%TipNode ! Loop through the blade nodes / elements

      ! Define the partial linear velocities (and their 1st derivatives) of the
      !   current node (point S(RNodes(J))) in the inertia frame.  Also define
      !   the overall linear velocity of the current node in the inertia frame.
      !   Also, define the portion of the linear acceleration of the current node
      !   in the inertia frame associated with everything but the QD2T()'s:

         EwHXrQS = CROSS_PRODUCT(  RtHSdat%AngVelEH, RtHSdat%rQS(:,K,J) )

         RtHSdat%PLinVelES(K,J,          :,:,:) = RtHSdat%PLinVelEQ(:,:,:)
         RtHSdat%PLinVelES(K,J,DOF_BF(K,1),0,:) = p%TwistedSF(K,1,1,J,0)                          *CoordSys%j1(K,:) &  !bjj: this line can be optimized
                                                + p%TwistedSF(K,2,1,J,0)                          *CoordSys%j2(K,:) &
                                                - (   p%AxRedBld(K,1,1,J)*x%QT ( DOF_BF(K,1) ) &
                                                    + p%AxRedBld(K,1,2,J)*x%QT ( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,3,J)*x%QT ( DOF_BE(K,1) )   )*CoordSys%j3(K,:)
         RtHSdat%PLinVelES(K,J,DOF_BE(K,1),0,:) = p%TwistedSF(K,1,3,J,0)                          *CoordSys%j1(K,:) &
                                                + p%TwistedSF(K,2,3,J,0)                          *CoordSys%j2(K,:) &
                                                - (   p%AxRedBld(K,3,3,J)*x%QT ( DOF_BE(K,1) ) &
                                                    + p%AxRedBld(K,2,3,J)*x%QT ( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,3,J)*x%QT ( DOF_BF(K,1) )   )*CoordSys%j3(K,:)
         RtHSdat%PLinVelES(K,J,DOF_BF(K,2),0,:) = p%TwistedSF(K,1,2,J,0)                          *CoordSys%j1(K,:) &
                                                + p%TwistedSF(K,2,2,J,0)                          *CoordSys%j2(K,:) &
                                                - (   p%AxRedBld(K,2,2,J)*x%QT ( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,2,J)*x%QT ( DOF_BF(K,1) ) &
                                                    + p%AxRedBld(K,2,3,J)*x%QT ( DOF_BE(K,1) )   )*CoordSys%j3(K,:)

         TmpVec1 = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%PLinVelES(K,J,DOF_BF(K,1),0,:) )
         TmpVec2 = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%PLinVelES(K,J,DOF_BE(K,1),0,:) )
         TmpVec3 = CROSS_PRODUCT( RtHSdat%AngVelEH, RtHSdat%PLinVelES(K,J,DOF_BF(K,2),0,:) )

         RtHSdat%PLinVelES(K,J,DOF_BF(K,1),1,:) = TmpVec1 &
                                                - (   p%AxRedBld(K,1,1,J)*x%QDT( DOF_BF(K,1) ) &
                                                    + p%AxRedBld(K,1,2,J)*x%QDT( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,3,J)*x%QDT( DOF_BE(K,1) )   )*CoordSys%j3(K,:)
         RtHSdat%PLinVelES(K,J,DOF_BE(K,1),1,:) = TmpVec2 &
                                                - (   p%AxRedBld(K,3,3,J)*x%QDT( DOF_BE(K,1) ) &
                                                    + p%AxRedBld(K,2,3,J)*x%QDT( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,3,J)*x%QDT( DOF_BF(K,1) )   )*CoordSys%j3(K,:)
         RtHSdat%PLinVelES(K,J,DOF_BF(K,2),1,:) = TmpVec3 &
                                                - (   p%AxRedBld(K,2,2,J)*x%QDT( DOF_BF(K,2) ) &
                                                    + p%AxRedBld(K,1,2,J)*x%QDT( DOF_BF(K,1) ) &
                                                    + p%AxRedBld(K,2,3,J)*x%QDT( DOF_BE(K,1) )   )*CoordSys%j3(K,:)

         LinVelHS                 = x%QDT( DOF_BF(K,1) )*RtHSdat%PLinVelES(K,J,DOF_BF(K,1),0,:) &
                                  + x%QDT( DOF_BE(K,1) )*RtHSdat%PLinVelES(K,J,DOF_BE(K,1),0,:) &
                                  + x%QDT( DOF_BF(K,2) )*RtHSdat%PLinVelES(K,J,DOF_BF(K,2),0,:)
         RtHSdat%LinAccESt(:,K,J) = x%QDT( DOF_BF(K,1) )*RtHSdat%PLinVelES(K,J,DOF_BF(K,1),1,:) &
                                  + x%QDT( DOF_BE(K,1) )*RtHSdat%PLinVelES(K,J,DOF_BE(K,1),1,:) &
                                  + x%QDT( DOF_BF(K,2) )*RtHSdat%PLinVelES(K,J,DOF_BF(K,2),1,:)

         RtHSdat%LinVelES(:,J,K)  = LinVelHS + RtHSdat%LinVelEZ
         DO I = 1,p%NPH   ! Loop through all DOFs associated with the angular motion of the hub (body H)

            TmpVec0 = CROSS_PRODUCT(   RtHSdat%PAngVelEH(p%PH(I),0,:), RtHSdat%rQS(:,K,J)            )  !bjj: this line can be optimized
            TmpVec1 = CROSS_PRODUCT(   RtHSdat%PAngVelEH(p%PH(I),0,:),     EwHXrQS        + LinVelHS )  !bjj: this line can be optimized
            TmpVec2 = CROSS_PRODUCT(   RtHSdat%PAngVelEH(p%PH(I),1,:), RtHSdat%rQS(:,K,J)            )  !bjj: this line can be optimized

            RtHSdat%PLinVelES(K,J,p%PH(I),0,:) = RtHSdat%PLinVelES(K,J,p%PH(I),0,:) + TmpVec0            !bjj: this line can be optimized
            RtHSdat%PLinVelES(K,J,p%PH(I),1,:) = RtHSdat%PLinVelES(K,J,p%PH(I),1,:) + TmpVec1 + TmpVec2  !bjj: this line can be optimized

            RtHSdat%LinVelES(:,J,K)          = RtHSdat%LinVelES(:,J,K)   + x%QDT(p%PH(I))*RtHSdat%PLinVelES(K,J,p%PH(I),0,:)  !bjj: this line can be optimized
            RtHSdat%LinAccESt(:,K,J)         = RtHSdat%LinAccESt(:,K,J)  + x%QDT(p%PH(I))*RtHSdat%PLinVelES(K,J,p%PH(I),1,:)  !bjj: this line can be optimized

         END DO ! I - all DOFs associated with the angular motion of the hub (body H)

      END DO !J = 0,p%TipNodes ! Loop through the blade nodes / elements
      
      
   !JASON: USE TipNode HERE INSTEAD OF BldNodes IF YOU ALLOCATE AND DEFINE n1, n2, n3, m1, m2, AND m3 TO USE TipNode.  THIS WILL REQUIRE THAT THE AERODYNAMIC AND STRUCTURAL TWISTS, AeroTwst() AND ThetaS(), BE KNOWN AT THE TIP!!!
      !IF (.NOT. p%BD4Blades) THEN
      !   RtHSdat%LinVelESm2(K) = DOT_PRODUCT( RtHSdat%LinVelES(:,p%TipNode,K), CoordSys%m2(K,p%BldNodes,:) )
      !END IF
            
   END DO !K = 1,p%NumBl


   RtHSdat%PLinVelEW(       :,:,:) = RtHSdat%PLinVelEO(:,:,:)
   DO I = 1,NPN   ! Loop through all DOFs associated with the angular motion of the nacelle (body N)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:), RtHSdat%rOW                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,0,:),     EwNXrOW                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEN(PN(I)   ,1,:), RtHSdat%rOW                 )

      RtHSdat%PLinVelEW(PN(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEW(PN(I)   ,0,:)
      RtHSdat%PLinVelEW(PN(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEW(PN(I)   ,1,:)

       LinAccEWt                   =  LinAccEWt + x%QDT(PN(I) )*RtHSdat%PLinVelEW(PN(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the nacelle (body N)


   RtHSdat%PLinVelEI(       :,:,:) = RtHSdat%PLinVelEW(:,:,:)
   DO I = 1,NPA   ! Loop through all DOFs associated with the angular motion of the tail (body A)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:), RtHSdat%rWI                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:),     EwAXrWI                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,1,:), RtHSdat%rWI                 )

      RtHSdat%PLinVelEI(PA(I),0,:) = TmpVec0    +                       RtHSdat%PLinVelEI(PA(I)   ,0,:)
      RtHSdat%PLinVelEI(PA(I),1,:) = TmpVec1    + TmpVec2 +             RtHSdat%PLinVelEI(PA(I)   ,1,:)

      RtHSdat%LinAccEIt            =  RtHSdat%LinAccEIt + x%QDT(PA(I) )*RtHSdat%PLinVelEI(PA(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the tail (body A)


   RtHSdat%PLinVelEJ(       :,:,:) = RtHSdat%PLinVelEW(:,:,:)
   DO I = 1,NPA   ! Loop through all DOFs associated with the angular motion of the tail (body A)

      TmpVec0 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:), RtHSdat%rWJ                 )
      TmpVec1 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:),     EwAXrWJ                 )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,1,:), RtHSdat%rWJ                 )

      RtHSdat%PLinVelEJ(PA(I),0,:) = TmpVec0    +               RtHSdat%PLinVelEJ(PA(I)   ,0,:)
      RtHSdat%PLinVelEJ(PA(I),1,:) = TmpVec1    + TmpVec2 +     RtHSdat%PLinVelEJ(PA(I)   ,1,:)

       RtHSdat%LinAccEJt           =  RtHSdat%LinAccEJt + x%QDT(PA(I) )*RtHSdat%PLinVelEJ(PA(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the tail (body A)

   RtHSdat%PLinVelEK(       :,:,:) = RtHSdat%PLinVelEW(:,:,:)
    LinVelEK               =  RtHSdat%LinVelEZ
   DO I = 1,NPA   ! Loop through all DOFs associated with the angular motion of the tail (body A)

      TmpVec0  = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:), RtHSdat%rWK                 )
      TmpVec1  = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,0,:),         EwAXrWK             )
      TmpVec2  = CROSS_PRODUCT( RtHSdat%PAngVelEA(PA(I)   ,1,:), RtHSdat%rWK                 )

      RtHSdat%PLinVelEK(PA(I),0,:) = TmpVec0    +                RtHSdat%PLinVelEK(PA(I)   ,0,:)
      RtHSdat%PLinVelEK(PA(I),1,:) = TmpVec1    + TmpVec2 +      RtHSdat%PLinVelEK(PA(I)   ,1,:)

       LinVelEK                    =   LinVelEK  + x%QDT(PA(I) )*RtHSdat%PLinVelEK(PA(I)   ,0,:)
       LinAccEKt                   =   LinAccEKt + x%QDT(PA(I) )*RtHSdat%PLinVelEK(PA(I)   ,1,:)

   ENDDO          ! I - all DOFs associated with the angular motion of the tail (body A)



   DO J = 0,p%TwrNodes  ! Loop through the tower nodes / elements


      ! Define the partial linear velocities (and their 1st derivatives) of the current node (point T(HNodes(J))) in the inertia frame.
      !  Also define the overall linear velocity of the current node in the inertia frame.
      !  Also, define the portion of the linear acceleration of the current node in the inertia frame associated with
      !    everything but the QD2T()'s:

      EwXXrZT                   = CROSS_PRODUCT(  RtHSdat%AngVelEX, RtHSdat%rZT(:,J) )

      RtHSdat%PLinVelET(J,       :,:,:) = RtHSdat%PLinVelEZ(:,:,:)  !bjj: can this line be optimized
      RtHSdat%PLinVelET(J,DOF_TFA1,0,:) = p%TwrFASF(1,J,0)*CoordSys%a1 - (   p%AxRedTFA(1,1,J)* x%QT(DOF_TFA1) &
                                                                           + p%AxRedTFA(1,2,J)* x%QT(DOF_TFA2)   )*CoordSys%a2  
      RtHSdat%PLinVelET(J,DOF_TSS1,0,:) = p%TwrSSSF(1,J,0)*CoordSys%a3 - (   p%AxRedTSS(1,1,J)* x%QT(DOF_TSS1) &
                                                                           + p%AxRedTSS(1,2,J)* x%QT(DOF_TSS2)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TFA2,0,:) = p%TwrFASF(2,J,0)*CoordSys%a1 - (   p%AxRedTFA(2,2,J)* x%QT(DOF_TFA2) &
                                                                           + p%AxRedTFA(1,2,J)* x%QT(DOF_TFA1)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TSS2,0,:) = p%TwrSSSF(2,J,0)*CoordSys%a3 - (   p%AxRedTSS(2,2,J)* x%QT(DOF_TSS2) &
                                                                           + p%AxRedTSS(1,2,J)* x%QT(DOF_TSS1)   )*CoordSys%a2

      TmpVec1 = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%PLinVelET(J,DOF_TFA1,0,:) )
      TmpVec2 = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%PLinVelET(J,DOF_TSS1,0,:) )
      TmpVec3 = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%PLinVelET(J,DOF_TFA2,0,:) )
      TmpVec4 = CROSS_PRODUCT( RtHSdat%AngVelEX, RtHSdat%PLinVelET(J,DOF_TSS2,0,:) )

      RtHSdat%PLinVelET(J,DOF_TFA1,1,:) = TmpVec1 - (   p%AxRedTFA(1,1,J)*x%QDT(DOF_TFA1) &
                                                      + p%AxRedTFA(1,2,J)*x%QDT(DOF_TFA2)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TSS1,1,:) = TmpVec2 - (   p%AxRedTSS(1,1,J)*x%QDT(DOF_TSS1) &
                                                      + p%AxRedTSS(1,2,J)*x%QDT(DOF_TSS2)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TFA2,1,:) = TmpVec3 - (   p%AxRedTFA(2,2,J)*x%QDT(DOF_TFA2) &
                                                      + p%AxRedTFA(1,2,J)*x%QDT(DOF_TFA1)   )*CoordSys%a2
      RtHSdat%PLinVelET(J,DOF_TSS2,1,:) = TmpVec4 - (   p%AxRedTSS(2,2,J)*x%QDT(DOF_TSS2) &
                                                      + p%AxRedTSS(1,2,J)*x%QDT(DOF_TSS1)   )*CoordSys%a2

              LinVelXT       = x%QDT(DOF_TFA1)*RtHSdat%PLinVelET(J,DOF_TFA1,0,:) &
                             + x%QDT(DOF_TSS1)*RtHSdat%PLinVelET(J,DOF_TSS1,0,:) &
                             + x%QDT(DOF_TFA2)*RtHSdat%PLinVelET(J,DOF_TFA2,0,:) &
                             + x%QDT(DOF_TSS2)*RtHSdat%PLinVelET(J,DOF_TSS2,0,:)
      RtHSdat%LinAccETt(:,J) = x%QDT(DOF_TFA1)*RtHSdat%PLinVelET(J,DOF_TFA1,1,:) &
                             + x%QDT(DOF_TSS1)*RtHSdat%PLinVelET(J,DOF_TSS1,1,:) &
                             + x%QDT(DOF_TFA2)*RtHSdat%PLinVelET(J,DOF_TFA2,1,:) &
                             + x%QDT(DOF_TSS2)*RtHSdat%PLinVelET(J,DOF_TSS2,1,:)

      RtHSdat%LinVelET(:,J)  = LinVelXT + RtHSdat%LinVelEZ
      DO I = 1,NPX   ! Loop through all DOFs associated with the angular motion of the platform (body X)

         TmpVec0   = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I),0,:), RtHSdat%rZT(:,J)            )
         TmpVec1   = CROSS_PRODUCT( RtHSdat%PAngVelEX(PX(I),0,:), EwXXrZT      + LinVelXT )

         RtHSdat%PLinVelET(J,PX(I),0,:) = RtHSdat%PLinVelET(J,PX(I),0,:) + TmpVec0
         RtHSdat%PLinVelET(J,PX(I),1,:) = RtHSdat%PLinVelET(J,PX(I),1,:) + TmpVec1

         RtHSdat%LinVelET( :,        J) = RtHSdat%LinVelET( :,        J) + x%QDT(PX(I))*RtHSdat%PLinVelET(J,PX(I),0,:)
         RtHSdat%LinAccETt(:,        J) = RtHSdat%LinAccETt(:,        J) + x%QDT(PX(I))*RtHSdat%PLinVelET(J,PX(I),1,:)

      ENDDO          ! I - all DOFs associated with the angular motion of the platform (body X)


   END DO ! J


END SUBROUTINE CalculateLinearVelPAcc
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CalculateForcesMoments( p, x, CoordSys, u, RtHSdat )
! This routine is used to calculate the forces and moments stored in other states that are used in
! both the CalcOutput and CalcContStateDeriv routines.
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    ! The coordinate systems that have been set for these states/time
   TYPE(ED_InputType),           INTENT(IN   )  :: u           ! The aero (blade) & nacelle forces/moments
   TYPE(ED_RtHndSide),           INTENT(INOUT)  :: RtHSdat     ! data from the RtHndSid module (contains positions to be set)

      ! Local variables
   REAL(ReKi)                   :: TmpVec    (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec1   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec2   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec3   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec4   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec5   (3)                                   ! A temporary vector used in various computations.
      
!REAL(ReKi)                   :: rSAerCen  (3)                                   ! Position vector from a blade analysis node (point S) on the current blade to the aerodynamic center associated with the element.
   REAL(ReKi), PARAMETER        :: FKAero   (3) = 0.0                              ! The tail fin aerodynamic force acting at point K, the center-of-pressure of the tail fin. (bjj: should be an input)
   REAL(ReKi), PARAMETER        :: MAAero   (3) = 0.0                              ! The tail fin aerodynamic moment acting at point K, the center-of-pressure of the tail fin. (bjj: should be an input)   
   
   
   INTEGER(IntKi)               :: I                                               ! Loops through some or all of the DOFs
   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades
   INTEGER(IntKi)               :: NodeNum                                         ! Node number for blade element (on a single mesh)
      
!.....................................
! Compute forces and moments from properties related to Aero inputs   
!FSTipDrag
!FSAero and MMAero
!.....................................
   DO K = 1,p%NumBl ! Loop through all blades
      
         ! Calculate the tip drag forces if necessary:
      !bjj: add this back when we've figured out how to handle the tip brakes:
      !RtHSdat%FSTipDrag(:,K) = OtherState%CoordSys%m2(K,p%BldNodes,:)*SIGN( 0.5*p%AirDens*(RtHSdat%LinVelESm2(K)**2)*u%TBDrCon(K), -1.*RtHSdat%LinVelESm2(K) )
      RtHSdat%FSTipDrag = 0.0_ReKi         ! Calculate the tip drag forces if necessary

      
   ! Calculate the normal and tangential aerodynamic forces and the aerodynamic
   !   pitching moment at the current element per unit span by calling AeroDyn,
   !   if necessary:
      
      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements


   ! Calculate the aerodynamic pitching moment arm (i.e., the position vector
   !   from point S on the blade to the aerodynamic center of the element):

         RtHSdat%rSAerCen(:,J,K) = p%rSAerCenn1(K,J)*CoordSys%n1(K,J,:) + p%rSAerCenn2(K,J)*CoordSys%n2(K,J,:)   

!        rPAerCen     = OtherState%RtHS%rPQ + OtherState%RtHS%rQS(:,K,J) + OtherState%RtHS%rSAerCen(:,J,K)     ! Position vector from teeter pin (point P)  to blade analysis node aerodynamic center.
!        rAerCen      =                       OtherState%RtHS%rS (:,K,J) + OtherState%RtHS%rSAerCen(:,J,K)     ! Position vector from inertial frame origin to blade analysis node aerodynamic center.
         

   ! fill FSAero() and MMAero() with the forces resulting from inputs u%BladeLn2Mesh(K)%Force(1:2,:) and u%BladeLn2Mesh(K)%Moment(3,:):
   ! [except, we're ignoring the additional nodes we added on the mesh end points]
   
         NodeNum = J ! we're ignoring the root and tip
         
         if (p%UseAD14) then
            RtHSdat%FSAero(:,K,J) = ( u%BladePtLoads(K)%Force(1,NodeNum) * CoordSys%te1(K,J,:) &
                                    + u%BladePtLoads(K)%Force(2,NodeNum) * CoordSys%te2(K,J,:) ) / p%DRNodes(J)

            RtHSdat%MMAero(:,K,J) = CROSS_PRODUCT( RtHSdat%rSAerCen(:,J,K), RtHSdat%FSAero(:,K,J) )&
                                  + u%BladePtLoads(K)%Moment(3,NodeNum)/p%DRNodes(J) * CoordSys%te3(K,J,:)        
         else
            RtHSdat%FSAero(1,K,J) =  u%BladePtLoads(K)%Force(1,NodeNum) / p%DRNodes(J)
            RtHSdat%FSAero(2,K,J) =  u%BladePtLoads(K)%Force(3,NodeNum) / p%DRNodes(J) 
            RtHSdat%FSAero(3,K,J) = -u%BladePtLoads(K)%Force(2,NodeNum) / p%DRNodes(J)

            RtHSdat%MMAero(1,K,J) =  u%BladePtLoads(K)%Moment(1,NodeNum) / p%DRNodes(J)
            RtHSdat%MMAero(2,K,J) =  u%BladePtLoads(K)%Moment(3,NodeNum) / p%DRNodes(J)
            RtHSdat%MMAero(3,K,J) = -u%BladePtLoads(K)%Moment(2,NodeNum) / p%DRNodes(J)
         end if
                     
         
      END DO !J
   END DO  ! K 
   
   
!.....................................
! PFrcS0B and PMomH0B  
!.....................................
DO K = 1,p%NumBl ! Loop through all blades

      ! Initialize the partial forces and moments (including those associated
      !   with the QD2T()'s and those that are not) at the blade root (point S(0))
      !   using the tip brake effects:

      RtHSdat%PFrcS0B(:,K,:) = 0.0 ! Initialize these partial
      RtHSdat%PMomH0B(:,K,:) = 0.0 ! forces and moments to zero
      DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

         TmpVec1 = -p%TipMass(K)*RtHSdat%PLinVelES(K,p%TipNode,p%DOFs%PSE(K,I),0,:)                            ! The portion of PFrcS0B associated with the tip brake

         RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I)) = TmpVec1
         RtHSdat%PMomH0B(:,K,p%DOFs%PSE(K,I)) = CROSS_PRODUCT( RtHSdat%rS0S(:,K,p%TipNode), TmpVec1 )          ! The portion of PMomH0B associated with the tip brake

      ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K  
   
   
      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements

      ! Integrate to find the partial forces and moments (including those associated
      !   with the QD2T()'s and those that are not) at the blade root (point S(0)):

         DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

            TmpVec1 = -p%BElmntMass(J,K)*RtHSdat%PLinVelES(K,J,p%DOFs%PSE(K,I),0,:)   ! The portion of PFrcS0B associated with blade element J

            RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I)) = RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I)) + TmpVec1
            RtHSdat%PMomH0B(:,K,p%DOFs%PSE(K,I)) = RtHSdat%PMomH0B(:,K,p%DOFs%PSE(K,I)) + &
                                                    CROSS_PRODUCT( RtHSdat%rS0S(:,K,J), TmpVec1 )                   ! The portion of PMomH0B associated with blade element J

         ENDDO             ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K
      END DO
      
      
   END DO     
   
 
!.....................................
! FrcS0Bt and MomH0Bt
!.....................................
   DO K = 1,p%NumBl ! Loop through all blades
   
      TmpVec1 = RtHSdat%FSTipDrag(:,K) - p%TipMass(K)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccESt(:,K,p%TipNode) ) ! The portion of FrcS0Bt associated with the tip brake
      RtHSdat%FrcS0Bt(:,K) = TmpVec1
      RtHSdat%MomH0Bt(:,K) = CROSS_PRODUCT(  RtHSdat%rS0S(:,K,p%TipNode), TmpVec1 )                                 ! The portion of MomH0Bt associated with the tip brake

      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements      
      
         TmpVec1 = RtHSdat%FSAero(:,K,J)*p%DRNodes(J) - p%BElmntMass(J,K)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccESt(:,K,J) ) ! The portion of FrcS0Bt associated with blade element J
         TmpVec2 = CROSS_PRODUCT( RtHSdat%rS0S(:,K,J), TmpVec1 )                                    ! The portion of MomH0Bt associated with blade element J
         TmpVec3 = RtHSdat%MMAero(:,K,J)*p%DRNodes(J)                                               ! The total external moment applied to blade element J

         RtHSdat%FrcS0Bt(:,K) = RtHSdat%FrcS0Bt(:,K) + TmpVec1
         RtHSdat%MomH0Bt(:,K) = RtHSdat%MomH0Bt(:,K) + TmpVec2 + TmpVec3
      
      END DO !J
      
   END DO !K   
         
      
!.....................................
! PFrcPRot AND PMomLPRot:  
!   ( requires PFrcS0B and PMomH0B)
!.....................................
      ! Initialize the partial forces and moments (including those associated
      !   with the QD2T()'s and those that are not) at the teeter pin (point P) using the hub mass effects:    

   RtHSdat%PFrcPRot  = 0.0   ! Initialize these partial
   RtHSdat%PMomLPRot = 0.0   ! forces and moments to zero
   DO I = 1,p%DOFs%NPCE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the hub center of mass (point C)

      TmpVec1 = -p%HubMass*RtHSdat%PLinVelEC(p%DOFs%PCE(I),0,:)     ! The portion of PFrcPRot  associated with the HubMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rPC, TmpVec1 )      ! The portion of PMomLPRot associated with the HubMass

      RtHSdat%PFrcPRot (:,p%DOFs%PCE(I)) = TmpVec1
      RtHSdat%PMomLPRot(:,p%DOFs%PCE(I)) = TmpVec2 - p%Hubg1Iner*CoordSys%g1*DOT_PRODUCT( CoordSys%g1, RtHSdat%PAngVelEH(p%DOFs%PCE(I),0,:) ) &
                                                   - p%Hubg2Iner*CoordSys%g2*DOT_PRODUCT( CoordSys%g2, RtHSdat%PAngVelEH(p%DOFs%PCE(I),0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the hub center of mass (point C)

   
   DO K = 1,p%NumBl ! Loop through all blades
   
         ! Calculate the position vector from the teeter pin to the blade root:
   
      !rPS0(:,K) = RtHSdat%rPQ + p%HubRad*CoordSys%j3(K,:)   ! Position vector from teeter pin (point P) to blade root (point S(0)).
            
      ! Add the blade effects to the partial forces and moments (including those associated with the QD2T()'s and those that are 
      !   not) at the teeter pin (point P):

      DO I = 1,p%DOFs%NPSE(K)  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K

         TmpVec = CROSS_PRODUCT( RtHSdat%rPS0(:,K), RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I)) ) ! The portion of PMomLPRot associated with PFrcS0B.

         RtHSdat%PFrcPRot (:,p%DOFs%PSE(K,I)) = RtHSdat%PFrcPRot (:,p%DOFs%PSE(K,I)) + RtHSdat%PFrcS0B(:,K,p%DOFs%PSE(K,I))
         RtHSdat%PMomLPRot(:,p%DOFs%PSE(K,I)) = RtHSdat%PMomLPRot(:,p%DOFs%PSE(K,I)) + RtHSdat%PMomH0B(:,K,p%DOFs%PSE(K,I))+TmpVec

      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of blade K


   END DO   ! K   

!.....................................
! FrcPRott and MomLPRott:
!   (requires FrcS0Bt and MomH0Bt)
!.....................................

   TmpVec1 = -p%HubMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccECt )                     ! The portion of FrcPRott  associated with the HubMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rPC, TmpVec1 )                                        ! The portion of MomLPRott associated with the HubMass
   TmpVec  = p%Hubg1Iner*CoordSys%g1*DOT_PRODUCT( CoordSys%g1, RtHSdat%AngVelEH ) &       ! = ( Hub inertia dyadic ) dot ( angular velocity of hub in the inertia frame )
           + p%Hubg2Iner*CoordSys%g2*DOT_PRODUCT( CoordSys%g2, RtHSdat%AngVelEH )
   TmpVec3 = CROSS_PRODUCT( -RtHSdat%AngVelEH, TmpVec )                                   ! = ( -angular velocity of hub in the inertia frame ) cross ( TmpVec )

   RtHSdat%FrcPRott(1)  = TmpVec1(1) + u%HubPtLoad%Force(1,1)
   RtHSdat%FrcPRott(2)  = TmpVec1(2) + u%HubPtLoad%Force(3,1)
   RtHSdat%FrcPRott(3)  = TmpVec1(3) - u%HubPtLoad%Force(2,1)
   
   RtHSdat%MomLPRott    = TmpVec2 + TmpVec3 - p%Hubg1Iner*CoordSys%g1*DOT_PRODUCT( CoordSys%g1, RtHSdat%AngAccEHt ) &
                                            - p%Hubg2Iner*CoordSys%g2*DOT_PRODUCT( CoordSys%g2, RtHSdat%AngAccEHt )                                          
      
   RtHSdat%MomLPRott(1) = RtHSdat%MomLPRott(1) + u%HubPtLoad%Moment(1,1)
   RtHSdat%MomLPRott(2) = RtHSdat%MomLPRott(2) + u%HubPtLoad%Moment(3,1)
   RtHSdat%MomLPRott(3) = RtHSdat%MomLPRott(3) - u%HubPtLoad%Moment(2,1)
   
   DO K = 1,p%NumBl ! Loop through all blades
   
         ! Calculate the position vector from the teeter pin to the blade root:
   
      !rPS0 = RtHSdat%rPQ + p%HubRad*OtherState%CoordSys%j3(K,:)   ! Position vector from teeter pin (point P) to blade root (point S(0)).
            
      TmpVec = CROSS_PRODUCT( RtHSdat%rPS0(:,K), RtHSdat%FrcS0Bt(:,K) )       ! The portion of MomLPRott associated with FrcS0Bt.

      RtHSdat%FrcPRott  = RtHSdat%FrcPRott  + RtHSdat%FrcS0Bt(:,K)
      RtHSdat%MomLPRott = RtHSdat%MomLPRott + RtHSdat%MomH0Bt(:,K) + TmpVec

   END DO   ! K

   
   ! Define the partial forces and moments (including those associated with
   !   the QD2T()'s and those that are not) at the specified point on the
   !   rotor-furl axis (point V) / nacelle (body N) using the structure that
   !   furls with the rotor, generator, and rotor effects.
   
!.....................................
! PMomNGnRt and PFrcVGnRt
!  (requires PMomLPRot and PFrcPRot)
!..................................... 
   RtHSdat%PFrcVGnRt = RtHSdat%PFrcPRot    ! Initialize these partial forces and
   RtHSdat%PMomNGnRt = RtHSdat%PMomLPRot   ! moments using the rotor effects
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      TmpVec = CROSS_PRODUCT( RtHSdat%rVP, RtHSdat%PFrcPRot(:,p%DOFs%SrtPS(I)) )  ! The portion of PMomNGnRt associated with the PFrcPRot

      RtHSdat%PMomNGnRt(:,p%DOFs%SrtPS(I)) = RtHSdat%PMomNGnRt(:,p%DOFs%SrtPS(I)) + TmpVec

   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPDE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the center of mass of the structure that furls with the rotor (not including rotor) (point D)

      TmpVec1 = -p%RFrlMass*RtHSdat%PLinVelED(p%DOFs%PDE(I)  ,0,:)           ! The portion of PFrcVGnRt associated with the RFrlMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rVD,              TmpVec1 )  ! The portion of PMomNGnRt associated with the RFrlMass

      RtHSdat%PFrcVGnRt(:,p%DOFs%PDE(I)  ) = RtHSdat%PFrcVGnRt(:,p%DOFs%PDE(I)  ) + TmpVec1

      RtHSdat%PMomNGnRt(:,p%DOFs%PDE(I)  ) = RtHSdat%PMomNGnRt(:,p%DOFs%PDE(I)  ) + TmpVec2                                   &
                                           - p%RrfaIner*CoordSys%rfa*DOT_PRODUCT( CoordSys%rfa, RtHSdat%PAngVelER(p%DOFs%PDE(I) ,0,:) ) &
                                           - p%GenIner*CoordSys%c1 *DOT_PRODUCT(  CoordSys%c1 , RtHSdat%PAngVelEG(p%DOFs%PDE(I) ,0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the center of mass of the structure that furls with the rotor (not including rotor) (point D)
   IF ( p%DOF_Flag(DOF_GeAz) )  THEN

      RtHSdat%PMomNGnRt(:,DOF_GeAz) = RtHSdat%PMomNGnRt(:,DOF_GeAz)                                             &     ! The previous loop (DO I = 1,NPDE) misses the DOF_GeAz-contribution to: ( Generator inertia dyadic ) dot ( partial angular velocity of the generator in the inertia frame )
                            -  p%GenIner*CoordSys%c1 *DOT_PRODUCT( CoordSys%c1, RtHSdat%PAngVelEG(DOF_GeAz,0,:) )     ! Thus, add this contribution if necessary.

   ENDIF   
   
!.....................................
! FrcVGnRtt and MomNGnRtt
!  (requires FrcPRott and MomLPRott)
!.....................................
   TmpVec1 = -p%RFrlMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEDt )                ! The portion of FrcVGnRtt associated with the RFrlMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rVD      ,  TmpVec1 )                             ! The portion of MomNGnRtt associated with the RFrlMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rVP      , RtHSdat%FrcPRott )                     ! The portion of MomNGnRtt associated with the FrcPRott
   TmpVec  = p%RrfaIner*CoordSys%rfa*DOT_PRODUCT( CoordSys%rfa, RtHSdat%AngVelER )    ! = ( R inertia dyadic ) dot ( angular velocity of structure that furls with the rotor in the inertia frame )
   TmpVec4 = CROSS_PRODUCT( -RtHSdat%AngVelER, TmpVec )                               ! = ( -angular velocity of structure that furls with the rotor in the inertia frame ) cross ( TmpVec )
   TmpVec  =  p%GenIner*CoordSys%c1* DOT_PRODUCT( CoordSys%c1 , RtHSdat%AngVelEG )    ! = ( Generator inertia dyadic ) dot ( angular velocity of generator in the inertia frame )
   TmpVec5 = CROSS_PRODUCT( -RtHSdat%AngVelEG, TmpVec )                               ! = ( -angular velocity of generator in the inertia frame ) cross ( TmpVec )

   RtHSdat%FrcVGnRtt = RtHSdat%FrcPRott  + TmpVec1
   RtHSdat%MomNGnRtt = RtHSdat%MomLPRott + TmpVec2 + TmpVec3 + TmpVec4 + TmpVec5            &
                     - p%RrfaIner*CoordSys%rfa*DOT_PRODUCT( CoordSys%rfa, RtHSdat%AngAccERt ) &
                     -  p%GenIner*CoordSys%c1 *DOT_PRODUCT( CoordSys%c1 , RtHSdat%AngAccEGt )


!.....................................
! PFrcWTail and PMomNTail
!.....................................
      ! Define the partial forces and moments (including those associated with the QD2T()'s and 
      !   those that are not) at the specified point on the tail-furl axis (point W) / nacelle (body N) using the tail effects.

   RtHSdat%PFrcWTail = 0.0   ! Initialize these partial
   RtHSdat%PMomNTail = 0.0   ! forces and moments to zero
   DO I = 1,p%DOFs%NPIE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tail boom center of mass (point I)

      TmpVec1 = -p%BoomMass*RtHSdat%PLinVelEI(p%DOFs%PIE(I),0,:)    ! The portion of PFrcWTail associated with the BoomMass
      TmpVec2 = -p%TFinMass*RtHSdat%PLinVelEJ(p%DOFs%PIE(I),0,:)    ! The portion of PFrcWTail associated with the TFinMass
      TmpVec3 = CROSS_PRODUCT( RtHSdat%rWI, TmpVec1 )                      ! The portion of PMomNTail associated with the BoomMass
      TmpVec4 = CROSS_PRODUCT( RtHSdat%rWJ, TmpVec2 )                      ! The portion of PMomNTail associated with the TFinMass

      RtHSdat%PFrcWTail(:,p%DOFs%PIE(I)) = TmpVec1 + TmpVec2
      RtHSdat%PMomNTail(:,p%DOFs%PIE(I)) = TmpVec3 + TmpVec4 - p%AtfaIner*CoordSys%tfa* &
                                                             DOT_PRODUCT( CoordSys%tfa, RtHSdat%PAngVelEA(p%DOFs%PIE(I),0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tail boom center of mass (point I)

!.....................................
! FrcWTailt and MomNTailt
!  (requires FKAero and MAAero)
!.....................................

   TmpVec1 = -p%BoomMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEIt )                 ! The portion of FrcWTailt associated with the BoomMass
   TmpVec2 = -p%TFinMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEJt )                 ! The portion of FrcWTailt associated with the TFinMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rWI      , TmpVec1 )                           ! The portion of MomNTailt associated with the BoomMass
   TmpVec4 = CROSS_PRODUCT( RtHSdat%rWJ      , TmpVec2 )                           ! The portion of MomNTailt associated with the TFinMass
   TmpVec  = p%AtfaIner*CoordSys%tfa*DOT_PRODUCT( CoordSys%tfa, RtHSdat%AngVelEA )   ! = ( A inertia dyadic ) dot ( angular velocity of the tail in the inertia frame )
   TmpVec5 = CROSS_PRODUCT( -RtHSdat%AngVelEA, TmpVec  )                           ! = ( -angular velocity of the tail in the inertia frame ) cross ( TmpVec )

   RtHSdat%FrcWTailt = FKAero + TmpVec1 + TmpVec2
   RtHSdat%MomNTailt = MAAero + TmpVec3 + TmpVec4 + TmpVec5         &
                     + CROSS_PRODUCT( RtHSdat%rWK      , FKAero  )  &                         ! The portion of MomNTailt associated with FKAero
                     - p%AtfaIner*CoordSys%tfa*DOT_PRODUCT( CoordSys%tfa, RtHSdat%AngAccEAt )   
   
!.....................................
! PFrcONcRt and PMomBNcRt
!  (requires PFrcVGnRt, PMomNGnRt, PFrcWTail, PMomNTail, )
!.....................................

   ! Define the partial forces and moments (including those associated with
   !   the QD2T()'s and those that are not) at the yaw bearing (point O) /
   !   base plate (body B) using the nacelle, generator, rotor, and tail effects.

   RtHSdat%PFrcONcRt = RtHSdat%PFrcVGnRt + RtHSdat%PFrcWTail   ! Initialize these partial forces and moments using
   RtHSdat%PMomBNcRt = RtHSdat%PMomNGnRt + RtHSdat%PMomNTail   ! the rotor, rotor-furl, generator, and tail effects
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      TmpVec = CROSS_PRODUCT( RtHSdat%rOV, RtHSdat%PFrcVGnRt(:,p%DOFs%SrtPS(I)) ) ! The portion of PMomBNcRt associated with the PFrcVGnRt

      RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) = RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) + TmpVec

   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPIE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tail boom center of mass (point I)

      TmpVec = CROSS_PRODUCT( RtHSdat%rOW, RtHSdat%PFrcWTail(:,p%DOFs%PIE(I)  ) ) ! The portion of PMomBNcRt associated with the PFrcWTail

      RtHSdat%PMomBNcRt(:,p%DOFs%PIE(I) ) = RtHSdat%PMomBNcRt(:,p%DOFs%PIE(I) ) + TmpVec

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tail boom center of mass (point I)
   DO I = 1,p%DOFs%NPUE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the nacelle center of mass (point U)

      TmpVec1 = -p%NacMass*RtHSdat%PLinVelEU(p%DOFs%PUE(I),0,:)              ! The portion of PFrcONcRt associated with the NacMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rOU,               TmpVec1 ) ! The portion of PMomBNcRt associated with the NacMass

      RtHSdat%PFrcONcRt(:,p%DOFs%PUE(I)  ) = RtHSdat%PFrcONcRt(:,p%DOFs%PUE(I) ) + TmpVec1
      RtHSdat%PMomBNcRt(:,p%DOFs%PUE(I)  ) = RtHSdat%PMomBNcRt(:,p%DOFs%PUE(I) ) + TmpVec2 - p%Nacd2Iner*CoordSys%d2* &
                                             DOT_PRODUCT( CoordSys%d2, RtHSdat%PAngVelEN(p%DOFs%PUE(I),0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the nacelle center of mass (point U)


!.....................................
! FrcONcRtt and MomBNcRtt
!  (requires FrcVGnRtt, MomNGnRtt, FrcWTailt, MomNTailt)
!.....................................

   TmpVec1 = -p%NacMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEUt )                ! The portion of FrcONcRtt associated with the NacMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rOU,           TmpVec1 )                         ! The portion of MomBNcRtt associated with the NacMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rOV, RtHSdat%FrcVGnRtt )                         ! The portion of MomBNcRtt associated with the FrcVGnRtt
   TmpVec4 = CROSS_PRODUCT( RtHSdat%rOW, RtHSdat%FrcWTailt )                         ! The portion of MomBNcRtt associated with the FrcWTailt
   TmpVec  = p%Nacd2Iner*CoordSys%d2*DOT_PRODUCT( CoordSys%d2, RtHSdat%AngVelEN )    ! = ( Nacelle inertia dyadic ) dot ( angular velocity of nacelle in the inertia frame )
    
   RtHSdat%FrcONcRtt = RtHSdat%FrcVGnRtt + RtHSdat%FrcWTailt + TmpVec1 + (/ u%NacelleLoads%Force(1,1), u%NacelleLoads%Force(3,1), -u%NacelleLoads%Force(2,1) /)
   
   RtHSdat%MomBNcRtt = RtHSdat%MomNGnRtt + RtHSdat%MomNTailt + TmpVec2 + TmpVec3 + TmpVec4  &
                        + CROSS_PRODUCT( -RtHSdat%AngVelEN, TmpVec    )                     &    ! = ( -angular velocity of nacelle in the inertia frame ) cross ( TmpVec ) &
                        - p%Nacd2Iner*CoordSys%d2*DOT_PRODUCT( CoordSys%d2, RtHSdat%AngAccENt ) &
                        + (/ u%NacelleLoads%Moment(1,1), u%NacelleLoads%Moment(3,1), -u%NacelleLoads%Moment(2,1) /)

!.....................................
! PFTHydro and PMFHydro   
!  (requires TwrAddedMass)   
!.....................................

   ! Compute the partial hydrodynamic forces and moments per unit length
   !   (including those associated with the QD2T()'s and those that are not) at the current tower element (point T) / (body F):

   ! NOTE: These forces are named PFTHydro, PMFHydro, FTHydrot, and MFHydrot. However, the names should not imply that the 
   !       forces are a result of hydrodynamic contributions only.  These tower forces contain contributions from any external 
   !       load acting on the tower other than loads transmitted from aerodynamics.  For example, these tower forces contain 
   !       contributions from foundation stiffness and damping [not floating] or mooring line restoring and damping,
   !       as well as hydrostatic and hydrodynamic contributions [offshore].

   DO J=1,p%TwrNodes
   
      RtHSdat%PFTHydro(:,J,:) = 0.0
      RtHSdat%PMFHydro(:,J,:) = 0.0
      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tower

         RtHSdat%PFTHydro(:,J,p%DOFs%PTE(I)) = &
                                CoordSys%z1*( - u%TwrAddedMass(DOF_Sg,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Sg,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Sg,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_Sg,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Sg,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Sg,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   ) &
                              - CoordSys%z3*( - u%TwrAddedMass(DOF_Sw,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Sw,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Sw,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_Sw,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Sw,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Sw,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   ) &
                              + CoordSys%z2*( - u%TwrAddedMass(DOF_Hv,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Hv,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Hv,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_Hv,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Hv,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Hv,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   )
         RtHSdat%PMFHydro(:,J,p%DOFs%PTE(I)) = &
                                CoordSys%z1*( - u%TwrAddedMass(DOF_R ,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_R ,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_R ,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_R ,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_R ,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_R ,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   ) &
                              - CoordSys%z3*( - u%TwrAddedMass(DOF_P ,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_P ,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_P ,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_P ,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_P ,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_P ,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   ) &
                              + CoordSys%z2*( - u%TwrAddedMass(DOF_Y ,DOF_Sg,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Y ,DOF_Sw,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Y ,DOF_Hv,J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,2) &
                                              - u%TwrAddedMass(DOF_Y ,DOF_R ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,1) &
                                              + u%TwrAddedMass(DOF_Y ,DOF_P ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,3) &
                                              - u%TwrAddedMass(DOF_Y ,DOF_Y ,J)*RtHSdat%PAngVelEF(J,p%DOFs%PTE(I),0,2)   )

      END DO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tower
   END DO !J

!.....................................
! FTHydrot and MFHydrot   
!  (requires TwrAddedMass)   
!.....................................
   
   DO J=1,p%TwrNodes
      RtHSdat%FTHydrot(:,J) = CoordSys%z1*( u%TowerPtLoads%Force(DOF_Sg,J)/p%DHNodes(J) &
                                                  - u%TwrAddedMass(DOF_Sg,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_Sg,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_Sg,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_Sg,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_Sg,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_Sg,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   ) &
                            - CoordSys%z3*( u%TowerPtLoads%Force(DOF_Sw,J)/p%DHNodes(J) &
                                                  - u%TwrAddedMass(DOF_Sw,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_Sw,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_Sw,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_Sw,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_Sw,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_Sw,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   ) &
                             + CoordSys%z2*( u%TowerPtLoads%Force(DOF_Hv,J)/p%DHNodes(J) &
                                                  - u%TwrAddedMass(DOF_Hv,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_Hv,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_Hv,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_Hv,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_Hv,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_Hv,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   )
      RtHSdat%MFHydrot(:,J) = CoordSys%z1*( u%TowerPtLoads%Moment(DOF_R-3,J)/p%DHNodes(J) &
                                                  - u%TwrAddedMass(DOF_R ,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_R ,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_R ,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_R ,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_R ,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_R ,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   ) &
                            - CoordSys%z3*( u%TowerPtLoads%Moment(DOF_P-3 ,J)/p%DHNodes(J) &
                                                  - u%TwrAddedMass(DOF_P ,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_P ,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_P ,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_P ,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_P ,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_P ,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   ) &
                            + CoordSys%z2*( u%TowerPtLoads%Moment(DOF_Y-3 ,J)/p%DHNodes(J) &
                                                  - u%TwrAddedMass(DOF_Y ,DOF_Sg,J)*RtHSdat%LinAccETt(1,J) &
                                                  + u%TwrAddedMass(DOF_Y ,DOF_Sw,J)*RtHSdat%LinAccETt(3,J) &
                                                  - u%TwrAddedMass(DOF_Y ,DOF_Hv,J)*RtHSdat%LinAccETt(2,J) &
                                                  - u%TwrAddedMass(DOF_Y ,DOF_R ,J)*RtHSdat%AngAccEFt(1,J) &
                                                  + u%TwrAddedMass(DOF_Y ,DOF_P ,J)*RtHSdat%AngAccEFt(3,J) &
                                                  - u%TwrAddedMass(DOF_Y ,DOF_Y ,J)*RtHSdat%AngAccEFt(2,J)   )

   END DO !J
   
!.....................................
! PFrcT0Trb and PMomX0Trb
!  (requires PFrcONcRt, PMomBNcRt, PFrcT0Trb, PMomX0Trb, PFTHydro, PMFHydro)
!.....................................

      ! Initialize the partial forces and moments (including those associated
      !   with the QD2T()'s and those that are not) at the tower base (point T(0))  using everything but the tower:

   RtHSdat%PFrcT0Trb = RtHSdat%PFrcONcRt   ! Initialize these partial forces and moments
   RtHSdat%PMomX0Trb = RtHSdat%PMomBNcRt   ! using all of the effects above the yaw bearing
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      TmpVec  = CROSS_PRODUCT(  RtHSdat%rT0O, RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) )   ! The portion of PMomX0Trb associated with the PFrcONcRt

      RtHSdat%PMomX0Trb(:,p%DOFs%SrtPS(I)) = RtHSdat%PMomX0Trb(:,p%DOFs%SrtPS(I)) + TmpVec

   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

      TmpVec1 = -p%YawBrMass*RtHSdat%PLinVelEO(p%DOFs%PTE(I),0,:)               ! The portion of PFrcT0Trb associated with the YawBrMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rT0O,               TmpVec1 )   ! The portion of PMomX0Trb associated with the YawBrMass

      RtHSdat%PFrcT0Trb(:,p%DOFs%PTE(I)  ) = RtHSdat%PFrcT0Trb(:,p%DOFs%PTE(I)  ) + TmpVec1
      RtHSdat%PMomX0Trb(:,p%DOFs%PTE(I)  ) = RtHSdat%PMomX0Trb(:,p%DOFs%PTE(I)  ) + TmpVec2

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing center of mass (point O)

   TmpVec1 = -p%YawBrMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEOt ) ! The portion of FrcT0Trbt associated with the YawBrMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rT0O,   TmpVec1 )               ! The portion of MomX0Trbt associated with the YawBrMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rT0O, RtHSdat%FrcONcRtt )               ! The portion of MomX0Trbt associated with the FrcONcRtt

   RtHSdat%FrcT0Trbt = RtHSdat%FrcONcRtt + TmpVec1
   RtHSdat%MomX0Trbt = RtHSdat%MomBNcRtt + TmpVec2 + TmpVec3   
   
   ! Integrate to find the total partial forces and moments (including those
   !   associated with the QD2T()'s and those that are not) at the tower base (point T(0)):
   DO J=1,p%TwrNodes

      DO I = 1,p%DOFs%NPTE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tower

         TmpVec1 = RtHSdat%PFTHydro(:,J,p%DOFs%PTE(I))*p%DHNodes(J) - p%TElmntMass(J)*RtHSdat%PLinVelET(J,p%DOFs%PTE(I),0,:)           ! The portion of PFrcT0Trb associated with tower element J
         TmpVec2 = CROSS_PRODUCT( RtHSdat%rT0T(:,J), TmpVec1 )                 ! The portion of PMomX0Trb associated with tower element J
         TmpVec3 = RtHSdat%PMFHydro(:,J,p%DOFs%PTE(I))*p%DHNodes(J)             ! The added moment applied at tower element J

         RtHSdat%PFrcT0Trb(:,p%DOFs%PTE(I)) = RtHSdat%PFrcT0Trb(:,p%DOFs%PTE(I)) + TmpVec1
         RtHSdat%PMomX0Trb(:,p%DOFs%PTE(I)) = RtHSdat%PMomX0Trb(:,p%DOFs%PTE(I)) + TmpVec2 + TmpVec3

      ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the tower

      TmpVec1 = ( RtHSdat%FTHydrot(:,J) )*p%DHNodes(J) &
              - p%TElmntMass(J)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccETt(:,J) )          ! The portion of FrcT0Trbt associated with tower element J
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rT0T(:,J), TmpVec1 )                                 ! The portion of MomX0Trbt associated with tower element J
      TmpVec3 = ( RtHSdat%MFHydrot(:,J) )*p%DHNodes(J)                                      ! The external moment applied to tower element J

      RtHSdat%FrcT0Trbt = RtHSdat%FrcT0Trbt + TmpVec1

      RtHSdat%MomX0Trbt = RtHSdat%MomX0Trbt + TmpVec2 + TmpVec3

   END DO !J

!.....................................
! PFZHydro and  PMXHydro  
!  ( requires PtfmAddedMass )
!.....................................   
   
   !..................................................................................................................................
   ! Compute the partial platform forces and moments (including those associated with the QD2T()'s and those that are not) at the
   ! platform reference (point Z) / (body X).
   !
   ! NOTE: These forces are named PFZHydro, PMXHydro, FZHydrot, and MXHydrot. However, the names should not imply that the forces
   !   are a result of hydrodynamic contributions only. These platform forces contain contributions from any external load acting
   !   on the platform other than loads transmitted from the wind turbine. For example, these platform forces contain contributions
   !   from foundation stiffness and damping [not floating] or mooring line restoring and damping [floating], as well as hydrostatic
   !   and hydrodynamic contributions [offshore].
   !bjj: OtherState%RtHS%PFZHydro, %PMXHydro, %FZHydrot, and %MXHydrot are not used in the output routine anymore
   !      (because of their dependence on inputs, u)

   RtHSdat%PFZHydro = 0.0
   RtHSdat%PMXHydro = 0.0
   DO I = 1,p%DOFs%NPYE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)

      RtHSdat%PFZHydro(p%DOFs%PYE(I),:) = - u%PtfmAddedMass(DOF_Sg,p%DOFs%PYE(I))*RtHSdat%PLinVelEZ(DOF_Sg,0,:) &
                                          - u%PtfmAddedMass(DOF_Sw,p%DOFs%PYE(I))*RtHSdat%PLinVelEZ(DOF_Sw,0,:) &
                                          - u%PtfmAddedMass(DOF_Hv,p%DOFs%PYE(I))*RtHSdat%PLinVelEZ(DOF_Hv,0,:)
      RtHSdat%PMXHydro(p%DOFs%PYE(I),:) = - u%PtfmAddedMass(DOF_R ,p%DOFs%PYE(I))*RtHSdat%PAngVelEX(DOF_R ,0,:) &
                                          - u%PtfmAddedMass(DOF_P ,p%DOFs%PYE(I))*RtHSdat%PAngVelEX(DOF_P ,0,:) &
                                          - u%PtfmAddedMass(DOF_Y ,p%DOFs%PYE(I))*RtHSdat%PAngVelEX(DOF_Y ,0,:)

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)

   RtHSdat%FZHydrot = u%PlatformPtMesh%Force(DOF_Sg,1)*RtHSdat%PLinVelEZ(DOF_Sg,0,:) &
                    + u%PlatformPtMesh%Force(DOF_Sw,1)*RtHSdat%PLinVelEZ(DOF_Sw,0,:) &
                    + u%PlatformPtMesh%Force(DOF_Hv,1)*RtHSdat%PLinVelEZ(DOF_Hv,0,:)
   RtHSdat%MXHydrot = u%PlatformPtMesh%Moment(DOF_R-3,1)*RtHSdat%PAngVelEX(DOF_R ,0,:) &
                    + u%PlatformPtMesh%Moment(DOF_P-3,1)*RtHSdat%PAngVelEX(DOF_P ,0,:) &
                    + u%PlatformPtMesh%Moment(DOF_Y-3,1)*RtHSdat%PAngVelEX(DOF_Y ,0,:)
   
!.....................................
! PFrcZAll and PMomXAll  
!  (requires PFrcT0Trb, PMomX0Trb, PFZHydro, PMXHydro )
!.....................................   

   ! Define the partial forces and moments (including those associated with the QD2T()'s and those that are not) at the
   !   platform reference (point Z) / (body X) using the turbine and platform effects:

   RtHSdat%PFrcZAll = RtHSdat%PFrcT0Trb ! Initialize these partial forces and moments
   RtHSdat%PMomXAll = RtHSdat%PMomX0Trb ! using the effects from the wind turbine
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      TmpVec = CROSS_PRODUCT( RtHSdat%rZT0, RtHSdat%PFrcT0Trb(:,p%DOFs%SrtPS(I)) )   ! The portion of PMomXAll associated with the PFrcT0Trb

      RtHSdat%PMomXAll(:,p%DOFs%SrtPS(I)) = RtHSdat%PMomXAll(:,p%DOFs%SrtPS(I)) + TmpVec

   ENDDO             ! I - All active (enabled) DOFs
   DO I = 1,p%DOFs%NPYE  ! Loop through all active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)

      TmpVec1 = -p%PtfmMass*RtHSdat%PLinVelEY(p%DOFs%PYE(I),0,:)                ! The portion of PFrcZAll associated with the PtfmMass
      TmpVec2 = CROSS_PRODUCT( RtHSdat%rZY ,               TmpVec1 )   ! The portion of PMomXAll associated with the PtfmMass

      RtHSdat%PFrcZAll(:,p%DOFs%PYE(I)) = RtHSdat%PFrcZAll(:,p%DOFs%PYE(I)  )        + RtHSdat%PFZHydro(p%DOFs%PYE(I),:) + TmpVec1
      RtHSdat%PMomXAll(:,p%DOFs%PYE(I)) = RtHSdat%PMomXAll(:,p%DOFs%PYE(I)  )        + RtHSdat%PMXHydro(p%DOFs%PYE(I),:) + TmpVec2 &
                                    - p%PtfmRIner*CoordSys%a1*DOT_PRODUCT( CoordSys%a1, RtHSdat%PAngVelEX(p%DOFs%PYE(I),0,:) )   &
                                    - p%PtfmYIner*CoordSys%a2*DOT_PRODUCT( CoordSys%a2, RtHSdat%PAngVelEX(p%DOFs%PYE(I),0,:) )   &
                                    - p%PtfmPIner*CoordSys%a3*DOT_PRODUCT( CoordSys%a3, RtHSdat%PAngVelEX(p%DOFs%PYE(I),0,:) )

   ENDDO          ! I - All active (enabled) DOFs that contribute to the QD2T-related linear accelerations of the platform center of mass (point Y)

!.....................................
! FrcZAllt and MomXAllt
!  (requires FrcT0Trbt, MomX0Trbt)
!.....................................

   TmpVec1 = -p%PtfmMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEYt  )                                              ! The portion of FrcZAllt associated with the PtfmMass
   TmpVec2 = CROSS_PRODUCT( RtHSdat%rZY      ,   TmpVec1 )                                                                      ! The portion of MomXAllt associated with the PtfmMass
   TmpVec3 = CROSS_PRODUCT( RtHSdat%rZT0     , RtHSdat%FrcT0Trbt )                                                      ! The portion of MomXAllt associated with the FrcT0Trbt
   TmpVec  = p%PtfmRIner*CoordSys%a1*DOT_PRODUCT( CoordSys%a1, RtHSdat%AngVelEX  ) &      ! = ( Platform inertia dyadic ) dot ( angular velocity of platform in the inertia frame )
           + p%PtfmYIner*CoordSys%a2*DOT_PRODUCT( CoordSys%a2, RtHSdat%AngVelEX  ) &
           + p%PtfmPIner*CoordSys%a3*DOT_PRODUCT( CoordSys%a3, RtHSdat%AngVelEX  )
   TmpVec4 = CROSS_PRODUCT( -RtHSdat%AngVelEX,   TmpVec  )                                                      ! = ( -angular velocity of platform in the inertia frame ) cross ( TmpVec )

   RtHSdat%FrcZAllt = RtHSdat%FrcT0Trbt + RtHSdat%FZHydrot + TmpVec1
   RtHSdat%MomXAllt = RtHSdat%MomX0Trbt + RtHSdat%MXHydrot + TmpVec2 + TmpVec3 + TmpVec4   
   
   
END SUBROUTINE CalculateForcesMoments
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FillAugMat( p, x, CoordSys, u, HSSBrTrq, RtHSdat, AugMat )
! This routine is used to populate the AugMat matrix for RtHS (CalcContStateDeriv)
!..................................................................................................................................

      ! Passed variables
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
   TYPE(ED_CoordSys),            INTENT(IN   )  :: CoordSys    ! The coordinate systems that have been set for these states/time
   TYPE(ED_InputType),           INTENT(IN   )  :: u           ! The aero blade forces/moments
   TYPE(ED_RtHndSide),           INTENT(INOUT)  :: RtHSdat     ! data from the RtHndSid module (contains positions to be set)
   REAL(ReKi),                   INTENT(IN )    :: HSSBrTrq    !  SIGN( u%HSSBrTrqC, x%QDT(DOF_GeAz) ) or corrected value from FixHSS
   REAL(R8Ki),                   INTENT(OUT)    :: AugMat(:,:) ! 
   
      ! Local variables
   REAL(ReKi)                   :: TmpVec    (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec1   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: TmpVec3   (3)                                   ! A temporary vector used in various computations.
   REAL(ReKi)                   :: GBoxTrq                                         ! Gearbox torque on the LSS side in N-m (calculated from inputs and parameters).
   REAL(ReKi)                   :: GBoxEffFac2                                     ! A second gearbox efficiency factor = ( 1 / GBoxEff^SgnPrvLSTQ - 1 )

   INTEGER(IntKi)               :: I                                               ! Loops through some or all of the DOFs
   INTEGER(IntKi)               :: J                                               ! Counter for elements
   INTEGER(IntKi)               :: K                                               ! Counter for blades
   INTEGER(IntKi)               :: L                                               ! Generic index

   
      ! Initialize the matrix:
      
   AugMat      = 0.0
   GBoxTrq    = ( u%GenTrq + HSSBrTrq )*ABS(p%GBRatio) ! bjj: do we use HSSBrTrqC or HSSBrTrq?
   
   DO K = 1,p%NumBl ! Loop through all blades
   

      ! Initialize the portions of the mass matrix on and below the diagonal associated with purely blade DOFs (these portions can't
      !   be calculated using partial loads) using the tip mass effects. 
      ! Also, initialize the portions of the forcing vector associated with purely blade DOFs (these portions can't be calculated 
      !   using partial loads) using the tip mass effects:
      ! NOTE: The vector subscript array, PSBE(), used in the following loops must be sorted from smallest to largest DOF index in 
      !       order for the loops to work to enter values only on and below the diagonal of AugMat():

   
      DO L = 1,p%DOFs%NPSBE(K)    ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the tip of blade K (point S(p%BldFlexL))
         DO I = L,p%DOFs%NPSBE(K) ! Loop through all active (enabled) blade DOFs greater than or equal to L
            AugMat(p%DOFs%PSBE(K,I),p%DOFs%PSBE(K,L)) = p%TipMass(K)*&
                                        DOT_PRODUCT( RtHSdat%PLinVelES(K, p%TipNode, p%DOFs%PSBE(K,I),0,:), &   ! [C(q,t)]B
                                                     RtHSdat%PLinVelES(K, p%TipNode, p%DOFs%PSBE(K,L),0,:)    )
         ENDDO             ! I - All active (enabled) blade DOFs greater than or equal to L
      ENDDO                ! L - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the tip of blade K (point S(p%BldFlexL))

      TmpVec1 = RtHSdat%FSTipDrag(:,K) - p%TipMass(K)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccESt(:,K,p%TipNode) ) ! The portion of FrcS0Bt associated with the tip brake
      DO I = 1,p%DOFs%NPSBE(K)    ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the tip of blade K (point S(p%BldFlexL))
            AugMat(p%DOFs%PSBE(K,I), p%NAug) = DOT_PRODUCT( RtHSdat%PLinVelES(K,p%TipNode,p%DOFs%PSBE(K,I),0,:), &   ! {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                              TmpVec1                               ) ! NOTE: TmpVec1 is still the portion of FrcS0Bt associated with the tip brake
      ENDDO                ! I - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the tip of blade K (point S(p%BldFlexL))
   
      

      DO J = 1,p%BldNodes ! Loop through the blade nodes / elements


      ! Integrate to find the portions of the mass matrix on and below the diagonal associated with purely blade DOFs (these portions
      !   can't be calculated using partial loads).  Also, integrate to find the portions of the forcing vector associated with
      !    purely blade DOFs (these portions can't be calculated using partial loads):
      ! NOTE: The vector subscript array, PSBE(), used in the following loops must
      !       be sorted from smallest to largest DOF index in order for the loops
      !       to work to enter values only on and below the diagonal of AugMat():
  
         DO L = 1,p%DOFs%NPSBE(K)    ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
            DO I = L,p%DOFs%NPSBE(K) ! Loop through all active (enabled) blade DOFs greater than or equal to L
               AugMat(p%DOFs%PSBE(K,I),p%DOFs%PSBE(K,L)) = AugMat(p%DOFs%PSBE(K,I),p%DOFs%PSBE(K,L)) + p%BElmntMass(J,K)*&
                                             DOT_PRODUCT( RtHSdat%PLinVelES(K,J,p%DOFs%PSBE(K,I),0,:), &           ! [C(q,t)]B
                                                          RtHSdat%PLinVelES(K,J,p%DOFs%PSBE(K,L),0,:)   )
            ENDDO             ! I - All active (enabled) blade DOFs greater than or equal to L
         ENDDO                ! L - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
      
         TmpVec1 = RtHSdat%FSAero(:,K,J)*p%DRNodes(J) - p%BElmntMass(J,K)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccESt(:,K,J) ) ! The portion of FrcS0Bt associated with blade element J
         TmpVec3 = RtHSdat%MMAero(:,K,J)*p%DRNodes(J)                                               ! The total external moment applied to blade element J
         DO I = 1,p%DOFs%NPSBE(K)    ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
               AugMat(p%DOFs%PSBE(K,I), p%NAug) = AugMat(p%DOFs%PSBE(K,I),     p%NAug)                      & ! {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                           + DOT_PRODUCT( RtHSdat%PLinVelES(K,J,p%DOFs%PSBE(K,I),0,:), TmpVec1 ) & ! NOTE: TmpVec1 is still the portion of FrcS0Bt associated with blade element J
                                           + DOT_PRODUCT( RtHSdat%PAngVelEM(K,J,p%DOFs%PSBE(K,I),0,:), TmpVec3 )   !       and TmpVec3 is still the total external moment applied to blade element J
         ENDDO                ! I - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade


      ENDDO ! J - Blade nodes / elements      
      
      
      

      ! Initialize the portions of the mass matrix below the diagonal associated
      !   with the teeter and pure blade DOFs using the partial loads at the teeter pin; only do this if necessary:

      IF ( ( p%NumBl == 2 ) .AND. ( p%DOF_Flag(DOF_Teet) ) )  THEN
         DO L = 1,p%DOFs%NPSBE(K) ! Loop through all active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
            AugMat(DOF_Teet,p%DOFs%PSBE(K,L)) = -DOT_PRODUCT( RtHSdat%PAngVelEH(DOF_Teet,0,:), &
                                                              RtHSdat%PMomLPRot(:,p%DOFs%PSBE(K,L)) )  ! [C(q,t)]B
         ENDDO             ! L - All active (enabled) blade DOFs that contribute to the QD2T-related linear accelerations of the blade
      ENDIF



      ! If the associated DOFs are enabled, add the blade elasticity and damping
      !   forces to the forcing vector (these portions can't be calculated using
      !   partial loads):

      IF ( p%DOF_Flag(DOF_BF(K,1)) )  THEN
         AugMat(    DOF_BF(K,1),p%NAug) = AugMat(DOF_BF(K,1),p%NAug)      & !
                                        - p%KBF(K,1,1)*x%QT( DOF_BF(K,1)) &
                                        - p%KBF(K,1,2)*x%QT( DOF_BF(K,2)) &
                                        - p%CBF(K,1,1)*x%QDT(DOF_BF(K,1)) &
                                        - p%CBF(K,1,2)*x%QDT(DOF_BF(K,2))
      ENDIF
      IF ( p%DOF_Flag(DOF_BF(K,2)) )  THEN
         AugMat(    DOF_BF(K,2),p%NAug) = AugMat(DOF_BF(K,2),p%NAug)      & ! {-f(qd,q,t)}ElasticB + {-f(qd,q,t)}DampB
                                        - p%KBF(K,2,1)*x%QT( DOF_BF(K,1)) &
                                        - p%KBF(K,2,2)*x%QT( DOF_BF(K,2)) &
                                        - p%CBF(K,2,1)*x%QDT(DOF_BF(K,1)) &
                                        - p%CBF(K,2,2)*x%QDT(DOF_BF(K,2))
      ENDIF
      IF ( p%DOF_Flag(DOF_BE(K,1)) )  THEN
         AugMat(    DOF_BE(K,1),p%NAug) = AugMat(DOF_BE(K,1),p%NAug)      & !
                                        - p%KBE(K,1,1)*x%QT( DOF_BE(K,1)) &
                                        - p%CBE(K,1,1)*x%QDT(DOF_BE(K,1))
      ENDIF
                  
      
   END DO !k
         
   
      ! Initialize the portions of the mass matrix on and below the diagonal
      !   associated with purely tower DOFs (these portions can't be calculated
      !   using partial loads) using the yaw bearing mass effects.
      !   Also, initialize the portions of the forcing vector associated with
      !   purely blade DOFs (these portions can't be calculated using partial
      !   loads) using the yaw bearing mass effects:
      ! NOTE: The vector subscript array, PTTE(), used in the following loops must
      !       be sorted from smallest to largest DOF index in order for the loops
      !       to work to enter values only on and below the diagonal of AugMat():

   DO L = 1,p%DOFs%NPTTE    ! Loop through all active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing (point O)
      DO I = L,p%DOFs%NPTTE ! Loop through all active (enabled) tower DOFs greater than or equal to L
         AugMat(p%DOFs%PTTE(I),p%DOFs%PTTE(L)) = p%YawBrMass*DOT_PRODUCT( RtHSdat%PLinVelEO(p%DOFs%PTTE(I),0,:), &     ! [C(q,t)]T of YawBrMass
                                                      RtHSdat%PLinVelEO(p%DOFs%PTTE(L),0,:)    )
      ENDDO          ! I - All active (enabled) tower DOFs greater than or equal to L
   ENDDO             ! L - All active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing (point O)

   TmpVec1 = -p%YawBrMass*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccEOt ) ! The portion of FrcT0Trbt associated with the YawBrMass
   DO I = 1,p%DOFs%NPTTE    ! Loop through all active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing (point O)
         AugMat(p%DOFs%PTTE(I),   p%NAug) =           DOT_PRODUCT( RtHSdat%PLinVelEO(p%DOFs%PTTE(I),0,:), &     ! {-f(qd,q,t)}T + {-f(qd,q,t)}GravT of YawBrMass
                                                      TmpVec1                   )   ! NOTE: TmpVec1 is still the portion of FrcT0Trbt associated with YawBrMass
   ENDDO             ! I - All active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the yaw bearing (point O)
   
   
   
   

   DO J = 1,p%TwrNodes

   !..................................................................................................................................
   ! Integrate to find the portions of the mass matrix on and below the diagonal associated with purely tower DOFs (these portions
   !   can't be calculated using partial loads).  Also, integrate to find the portions of the forcing vector associated with purely
   !   tower DOFs (these portions can't be calculated using partial loads).
   ! NOTE: The vector subscript array, PTTE(), used in the following loops must be sorted from smallest to largest DOF index in order
   !   for the loops to work to enter values only on and below the diagonal of AugMat():
   !..................................................................................................................................

      DO L = 1,p%DOFs%NPTTE    ! Loop through all active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the tower
         DO I = L,p%DOFs%NPTTE ! Loop through all active (enabled) tower DOFs greater than or equal to L
            AugMat(p%DOFs%PTTE(I),p%DOFs%PTTE(L)) = AugMat(p%DOFs%PTTE(I),p%DOFs%PTTE(L))  &
                                                  + p%TElmntMass(J)   *DOT_PRODUCT( RtHSdat%PLinVelET(J,p%DOFs%PTTE(I),0,:),  &
                                                                              RtHSdat%PLinVelET(J,p%DOFs%PTTE(L),0,:) ) &   ! [C(q,t)]T + [C(q,t)]HydroT
                                                  - p%DHNodes(J)*DOT_PRODUCT( RtHSdat%PLinVelET(J,p%DOFs%PTTE(I),0,:),  &
                                                                              RtHSdat%PFTHydro (:,J,p%DOFs%PTTE(L)  ) ) &
                                                  - p%DHNodes(J)*DOT_PRODUCT( RtHSdat%PAngVelEF(J,p%DOFs%PTTE(I),0,:),  &
                                                                              RtHSdat%PMFHydro (:,J,p%DOFs%PTTE(L)  ) )
         ENDDO                 ! I - All active (enabled) tower DOFs greater than or equal to L
      ENDDO                    ! L - All active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the tower

      TmpVec1 = ( RtHSdat%FTHydrot(:,J) )*p%DHNodes(J) &
              - p%TElmntMass(J)*( p%Gravity*CoordSys%z2 + RtHSdat%LinAccETt(:,J) )          ! The portion of FrcT0Trbt associated with tower element J
      TmpVec3 = ( RtHSdat%MFHydrot(:,J) )*p%DHNodes(J)             ! The external moment applied to tower element J
      DO I = 1,p%DOFs%NPTTE    ! Loop through all active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the tower
            AugMat(p%DOFs%PTTE(I),        p%NAug) = AugMat(p%DOFs%PTTE(I),   p%NAug)                         &                 ! {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT
                                                  +  DOT_PRODUCT( RtHSdat%PLinVelET(J,p%DOFs%PTTE(I),0,:), TmpVec1        ) &  ! NOTE: TmpVec1 is still the portion of FrcT0Trbt associated with tower element J
                                                  +  DOT_PRODUCT( RtHSdat%PAngVelEF(J,p%DOFs%PTTE(I),0,:), TmpVec3        )    !       and TmpVec3 is still the total external moment to tower element J
      ENDDO                    ! I - All active (enabled) tower DOFs that contribute to the QD2T-related linear accelerations of the tower

   ENDDO ! J - Tower nodes / elements   
   
   !..................................................................................................................................
   ! If the associated DOFs are enabled, add the tower elasticity and damping forces to the forcing vector (these portions can't be
   !   calculated using partial loads):
   !..................................................................................................................................

   IF ( p%DOF_Flag(DOF_TFA1) )  THEN
      AugMat(    DOF_TFA1,p%NAug) = AugMat(DOF_TFA1,p%NAug)                                   &
                                  - p%KTFA(1,1)*x%QT( DOF_TFA1) - p%KTFA(1,2)*x%QT( DOF_TFA2) &                                     !
                                  - p%CTFA(1,1)*x%QDT(DOF_TFA1) - p%CTFA(1,2)*x%QDT(DOF_TFA2)
   ENDIF
   IF ( p%DOF_Flag(DOF_TSS1) )  THEN
      AugMat(    DOF_TSS1,p%NAug) = AugMat(DOF_TSS1,p%NAug)                                   &
                                  - p%KTSS(1,1)*x%QT( DOF_TSS1) - p%KTSS(1,2)*x%QT( DOF_TSS2) &                                     ! {-f(qd,q,t)}ElasticT + {-f(qd,q,t)}DampT
                                  - p%CTSS(1,1)*x%QDT(DOF_TSS1) - p%CTSS(1,2)*x%QDT(DOF_TSS2)
   ENDIF
   IF ( p%DOF_Flag(DOF_TFA2) )  THEN
      AugMat(    DOF_TFA2,p%NAug) = AugMat(DOF_TFA2,p%NAug)                                   &
                                  - p%KTFA(2,1)*x%QT( DOF_TFA1) - p%KTFA(2,2)*x%QT( DOF_TFA2) &                                     !
                                  - p%CTFA(2,1)*x%QDT(DOF_TFA1) - p%CTFA(2,2)*x%QDT(DOF_TFA2)
   ENDIF
   IF ( p%DOF_Flag(DOF_TSS2) )  THEN
      AugMat(    DOF_TSS2,p%NAug) = AugMat(DOF_TSS2,p%NAug)                                   &
                                  - p%KTSS(2,1)*x%QT( DOF_TSS1) - p%KTSS(2,2)*x%QT( DOF_TSS2) &                                     !
                                  - p%CTSS(2,1)*x%QDT(DOF_TSS1) - p%CTSS(2,2)*x%QDT(DOF_TSS2)
   ENDIF
   
   
   
!..................................................................................................................................
! Now that all of the partial loads have been found, let's fill in the portions of the mass matrix on and below the diagonal that
! may be calculated with the help of the partial loads.
! Also, let's fill in the portions of the forcing vector that may be calculated with the help of the partial loads.
! Also let's add in additional terms to the forcing function that can't be added with the help of the partial loads.
!
! NOTE: The vector subscript array, SrtPS(), used in the following loops must be sorted from smallest to largest DOF index in order
!   for the loops to work to enter values only on and below the diagonal of AugMat():
!..................................................................................................................................

   IF ( p%DOF_Flag (DOF_Sg  ) )  THEN
      DO I = p%DOFs%Diag(DOF_Sg  ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Sg  ) = -1.*DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Sg ,0,:), RtHSdat%PFrcZAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Sg  ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Sg ,0,:), RtHSdat%FrcZAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}R + {-f(qd,q,t)}H + {-f(qd,q,t)}B + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_Sw  ) )  THEN
      DO I = p%DOFs%Diag(DOF_Sw  ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Sw  ) = -1.*DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Sw ,0,:), RtHSdat%PFrcZAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Sw  ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Sw ,0,:), RtHSdat%FrcZAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}R + {-f(qd,q,t)}H + {-f(qd,q,t)}B + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_Hv  ) )  THEN
      DO I = p%DOFs%Diag(DOF_Hv  ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Hv  ) = -1.*DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Hv ,0,:), RtHSdat%PFrcZAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Hv  ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PLinVelEZ(DOF_Hv ,0,:), RtHSdat%FrcZAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}GravX + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_R   ) )  THEN
      DO I = p%DOFs%Diag(DOF_R   ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal

         AugMat(p%DOFs%SrtPS(I),DOF_R   ) = -1.*DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_R  ,0,:), RtHSdat%PMomXAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_R   ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_R  ,0,:), RtHSdat%MomXAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}GravX + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_P   ) )  THEN
      DO I = p%DOFs%Diag(DOF_P   ),p%DOFs%NActvDOF    ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_P   ) = -1.*DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_P  ,0,:), RtHSdat%PMomXAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                                                             ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_P            ,p%NAug) =     DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_P  ,0,:), RtHSdat%MomXAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}GravX + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
   END IF

   IF ( p%DOF_Flag (DOF_Y   ) )  THEN
      DO I = p%DOFs%Diag(DOF_Y   ),p%DOFs%NActvDOF    ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Y   ) = -1.*DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_Y  ,0,:), RtHSdat%PMomXAll (:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]X + [C(q,t)]HydroX + [C(q,t)]T + [C(q,t)]HydroT + [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                                                             ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Y   ,         p%NAug) =     DOT_PRODUCT( RtHSdat%PAngVelEX(DOF_Y  ,0,:), RtHSdat%MomXAllt              )        ! {-f(qd,q,t)}X + {-f(qd,q,t)}GravX + {-f(qd,q,t)}HydroX + {-f(qd,q,t)}T + {-f(qd,q,t)}GravT + {-f(qd,q,t)}AeroT + {-f(qd,q,t)}HydroT + {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
   ENDIF

   IF ( p%DOF_Flag (DOF_TFA1) )  THEN
      DO I = p%DOFs%Diag(DOF_TFA1),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_TFA1) = AugMat(p%DOFs%SrtPS(I),DOF_TFA1)                             &
                                          -  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TFA1,0,:),       &
                                                          RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) ) &                          ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
                                          -  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TFA1,0,:),       &
                                                          RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_TFA1,         p%NAug) = AugMat(DOF_TFA1,    p%NAug)                                  &
                                          +  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TFA1,0,:), RtHSdat%FrcONcRtt  ) &   ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                          +  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TFA1,0,:), RtHSdat%MomBNcRtt  )
   ENDIF

   IF ( p%DOF_Flag (DOF_TSS1) )  THEN
      DO I = p%DOFs%Diag(DOF_TSS1),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_TSS1) = AugMat(p%DOFs%SrtPS(I),DOF_TSS1)                             &
                                          -  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TSS1,0,:),       &
                                                          RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) ) &                          ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
                                          -  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TSS1,0,:),       &
                                                          RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_TSS1,         p%NAug) = AugMat(DOF_TSS1,    p%NAug)                                  &
                                          +  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TSS1,0,:), RtHSdat%FrcONcRtt  ) &   ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                          +  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TSS1,0,:), RtHSdat%MomBNcRtt  )
   ENDIF

   IF ( p%DOF_Flag (DOF_TFA2) )  THEN
      DO I = p%DOFs%Diag(DOF_TFA2),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_TFA2) = AugMat(p%DOFs%SrtPS(I),DOF_TFA2)                             &
                                          -  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TFA2,0,:),       &
                                                          RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) ) &                          ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
                                          -  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TFA2,0,:),       &
                                                          RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_TFA2,         p%NAug) = AugMat(DOF_TFA2,    p%NAug)                                  &
                                          +  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TFA2,0,:), RtHSdat%FrcONcRtt  ) &   ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                          +  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TFA2,0,:), RtHSdat%MomBNcRtt  )
   ENDIF

   IF ( p%DOF_Flag (DOF_TSS2) )  THEN
      DO I = p%DOFs%Diag(DOF_TSS2),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_TSS2) = AugMat(p%DOFs%SrtPS(I),DOF_TSS2)                             &
                                          -  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TSS2,0,:),       &
                                                          RtHSdat%PFrcONcRt(:,p%DOFs%SrtPS(I)) ) &                          ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
                                          -  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TSS2,0,:),       &
                                                          RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_TSS2,         p%NAug) = AugMat(DOF_TSS2,    p%NAug)                                  &
                                          +  DOT_PRODUCT( RtHSdat%PLinVelEO(DOF_TSS2,0,:), RtHSdat%FrcONcRtt  ) &   ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                          +  DOT_PRODUCT( RtHSdat%PAngVelEB(DOF_TSS2,0,:), RtHSdat%MomBNcRtt  )
   ENDIF
   
   IF ( p%DOF_Flag (DOF_Yaw ) )  THEN
      DO I = p%DOFs%Diag(DOF_Yaw ),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_Yaw ) = -DOT_PRODUCT( RtHSdat%PAngVelEN(DOF_Yaw ,0,:), RtHSdat%PMomBNcRt(:,p%DOFs%SrtPS(I)) )   ! [C(q,t)]N + [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B + [C(q,t)]A
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_Yaw ,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEN(DOF_Yaw ,0,:), RtHSdat%MomBNcRtt             ) &        ! {-f(qd,q,t)}N + {-f(qd,q,t)}GravN + {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB + {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                                              + u%YawMom                                                            ! + {-f(qd,q,t)}SpringYaw  + {-f(qd,q,t)}DampYaw; NOTE: The neutral yaw rate, YawRateNeut, defaults to zero.  It is only used for yaw control.
   ENDIF
   
   
   IF ( p%DOF_Flag (DOF_RFrl) )  THEN
      DO I = p%DOFs%Diag(DOF_RFrl),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_RFrl) = -DOT_PRODUCT( RtHSdat%PAngVelER(DOF_RFrl,0,:),       &
                                                          RtHSdat%PMomNGnRt(:,p%DOFs%SrtPS(I)) )                            ! [C(q,t)]R + [C(q,t)]G + [C(q,t)]H + [C(q,t)]B
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_RFrl,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelER(DOF_RFrl,0,:), RtHSdat%MomNGnRtt  ) &   ! {-f(qd,q,t)}R + {-f(qd,q,t)}GravR + {-f(qd,q,t)}G + {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                                              +  RtHSdat%RFrlMom                            ! + {-f(qd,q,t)}SpringRF + {-f(qd,q,t)}DampRF
   ENDIF

   TmpVec = p%GenIner*CoordSys%c1*DOT_PRODUCT( CoordSys%c1, RtHSdat%PAngVelEG(DOF_GeAz,0,:) )  ! = ( generator inertia dyadic ) Dot ( partial angular velocity of G in E for DOF_GeAz )

   IF ( p%DOF_Flag (DOF_GeAz) )  THEN
      DO I = p%DOFs%Diag(DOF_GeAz),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_GeAz) = -DOT_PRODUCT( RtHSdat%PAngVelEL(DOF_GeAz,0,:), RtHSdat%PMomLPRot(:,p%DOFs%SrtPS(I)) )! [C(q,t)]H + [C(q,t)]B
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_GeAz,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEL(DOF_GeAz,0,:), RtHSdat%MomLPRott             ) &     ! {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                              -  GBoxTrq                                                         ! + {-f(qd,q,t)}Gen + {-f(qd,q,t)}Brake


      ! The previous loop (DO I = p%DOFs%Diag(DOF_GeAz),p%DOFs%NActvDOF) misses the
      !   generator inertia-contribution to the mass matrix and forcing function.
      !   Thus, add these in as well:


         AugMat(DOF_GeAz,       DOF_GeAz) = AugMat(DOF_GeAz,DOF_GeAz)                                    &
                                            +  DOT_PRODUCT( RtHSdat%PAngVelEG(DOF_GeAz,0,:), TmpVec                )             ! [C(q,t)]G
         AugMat(DOF_GeAz,         p%NAug) = AugMat(DOF_GeAz,  p%NAug)                                    &
                                            -  DOT_PRODUCT( RtHSdat%AngAccEGt              , TmpVec                )             ! {-f(qd,q,t)}G


   ENDIF

   IF ( p%DOF_Flag (DOF_DrTr) )  THEN
      DO I = p%DOFs%Diag(DOF_DrTr),p%DOFs%NActvDOF   ! Loop through all active (enabled) DOFs on or below the diagonal
         AugMat(p%DOFs%SrtPS(I),DOF_DrTr) = -DOT_PRODUCT( RtHSdat%PAngVelEL(DOF_DrTr,0,:), RtHSdat%PMomLPRot(:,p%DOFs%SrtPS(I)) ) ! [C(q,t)]H + [C(q,t)]B
      ENDDO                            ! I - All active (enabled) DOFs on or below the diagonal
         AugMat(DOF_DrTr,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEL(DOF_DrTr,0,:), RtHSdat%MomLPRott             ) &      ! {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                          -  p%DTTorSpr*x%QT (DOF_DrTr)                                    &      ! + {-f(qd,q,t)}ElasticDrive
                                                          -  p%DTTorDmp*x%QDT(DOF_DrTr)                                           ! + {-f(qd,q,t)}DampDrive
   ENDIF

   IF ( p%DOF_Flag (DOF_TFrl) )  THEN
      ! The tail-furl DOF does not affect any DOF index larger than DOF_TFrl.  Therefore, there is no need to perform the loop: DO I = Diag(DOF_TFrl),NActvDOF
         AugMat(DOF_TFrl,       DOF_TFrl) = -DOT_PRODUCT( RtHSdat%PAngVelEA(DOF_TFrl,0,:), RtHSdat%PMomNTail(:,DOF_TFrl) )        ! [C(q,t)]A
         AugMat(DOF_TFrl,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEA(DOF_TFrl,0,:), RtHSdat%MomNTailt             ) &      ! {-f(qd,q,t)}A + {-f(qd,q,t)}GravA + {-f(qd,q,t)}AeroA
                                                              +  RtHSdat%TFrlMom                                                  ! + {-f(qd,q,t)}SpringTF + {-f(qd,q,t)}DampTF
   ENDIF

   IF ( ( p%NumBl == 2 ) .AND. ( p%DOF_Flag(DOF_Teet) ) )  THEN
      ! The teeter DOF does not affect any DOF index larger than DOF_Teet.  Therefore, there is no need to perform the loop: DO I = Diag(DOF_Teet),NActvDOF
         AugMat(DOF_Teet,       DOF_Teet) = -DOT_PRODUCT( RtHSdat%PAngVelEH(DOF_Teet,0,:), RtHSdat%PMomLPRot(:,DOF_Teet) )        ! [C(q,t)]H + [C(q,t)]B
         AugMat(DOF_Teet,         p%NAug) =  DOT_PRODUCT( RtHSdat%PAngVelEH(DOF_Teet,0,:), RtHSdat%MomLPRott             ) &      ! {-f(qd,q,t)}H + {-f(qd,q,t)}GravH + {-f(qd,q,t)}B + {-f(qd,q,t)}GravB + {-f(qd,q,t)}AeroB
                                                              +  RtHSdat%TeetMom                                                  ! + {-f(qd,q,t)}SpringTeet + {-f(qd,q,t)}DampTeet
   ENDIF
   !..................................................................................................................................
   ! So far, we have only filled in the portions of the mass matrix on and below the diagonal.  Because the mass matrix is symmetric
   !   up to this point, let's fill in the portion above the diagonal by mirroring the values from below:
   ! NOTE: The vector subscript array, SrtPS(), used in the following loops must be sorted from smallest to largest DOF index in order
   !   for the loops to work to enter values only on and below the diagonal of AugMat():
   !..................................................................................................................................

      DO L = 2,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs above the diagonal (columns)
         DO I = 1,L-1   ! Loop through all active (enabled) DOFs above the diagonal (rows)
            AugMat(p%DOFs%SrtPS(I),p%DOFs%SrtPS(L)) = AugMat(p%DOFs%SrtPS(L),p%DOFs%SrtPS(I))
         ENDDO          ! I - All active (enabled) DOFs above the diagonal (rows)
      ENDDO             ! L - All active (enabled) DOFs above the diagonal (columns)

   ! Let's add the gearbox friction terms to the mass matrix and forcing
   !   function.  These only effect the equation for the generator azimuth DOF.
   ! NOTE: the MASS MATRIX WILL NO LONGER BE SYMMETRIC after adding these
   !       terms, unless the gearbox efficiency, GBoxEff, was set to 100%:
   
   
   GBoxEffFac2 = ( 1.0/RtHSdat%GBoxEffFac - 1.0 ) ! = ( 1 / GBoxEff^SgnPrvLSTQ - 1 )
   !TmpVec = p%GenIner*CoordSys%c1*DOT_PRODUCT( CoordSys%c1, RtHSdat%PAngVelEG(DOF_GeAz,0,:) )  ! = ( generator inertia dyadic ) Dot ( partial angular velocity of G in E for DOF_GeAz )

   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      AugMat(DOF_GeAz,p%DOFs%SrtPS(I)) = AugMat(DOF_GeAz,p%DOFs%SrtPS(I)) &                                          ! NOTE: TmpVec is still = ( generator inertia dyadic ) Dot ( partial angular velocity of G in E for DOF_GeAz ) in the following equation
                                + GBoxEffFac2*  DOT_PRODUCT( RtHSdat%PAngVelEG(p%DOFs%SrtPS(I),0,:), TmpVec )        ! [C(q,t)]GBFric

   ENDDO             ! I - All active (enabled) DOFs

   AugMat(   DOF_GeAz,    p%NAug) = AugMat(DOF_GeAz,    p%NAug) &                                                    ! NOTE: TmpVec is still = ( generator inertia dyadic ) Dot ( partial angular velocity of G in E for DOF_GeAz ) in the following equation
                                - GBoxEffFac2*( DOT_PRODUCT( RtHSdat%AngAccEGt              , TmpVec ) + GBoxTrq )   ! {-f(qd,q,t)}GBFric

   
   
   
END SUBROUTINE FillAugMat
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_AllocOutput( p, OtherState, u, y, ErrStat, ErrMsg )
! This routine allocates the arrays and meshes stored in the ED_OutputType data structure (y), based on the parameters (p). 
! Inputs (u) are included only so that output meshes can be siblings of the inputs.
! The routine assumes that the arrays/meshes are not currently allocated (It will produce a fatal error otherwise.)
!..................................................................................................................................

   TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(ED_OtherStateType),      INTENT(IN   )  :: OtherState  ! Other states (initial positions, set in Init_Inputs())
   TYPE(ED_InputType),           INTENT(INOUT)  :: u           ! Input meshes (sibling)
   TYPE(ED_OutputType),          INTENT(INOUT)  :: y           ! Outputs to be allocated
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      
   
   ! local variables
   REAL(R8Ki)                                   :: Orientation(3,3) 
   REAL(ReKi)                                   :: Position(3) 
   INTEGER(IntKi)                               :: NodeNum     ! node number
   INTEGER(IntKi)                               :: J, K        ! loop counters
   INTEGER(IntKi)                               :: ErrStat2    ! The error identifier (ErrStat)
   CHARACTER(1024)                              :: ErrMsg2     ! The error message (ErrMsg)
   
   
      ! initialize variables:
      
   ErrStat = ErrID_None
   ErrMsg  = ""
      

   CALL AllocAry( y%WriteOutput, p%NumOuts, 'WriteOutput', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   CALL AllocAry( y%BlPitch, p%NumBl, 'BlPitch', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
   !.......................................................
   ! Create Line2 Mesh for motion outputs on blades:
   !.......................................................
   IF ( .NOT. p%BD4Blades) THEN
      ALLOCATE( y%BladeLn2Mesh(p%NumBl), Stat=ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL CheckError( ErrID_Fatal, 'ED: Could not allocate space for y%BladeLn2Mesh{p%NumBl}' )
         RETURN
      END IF
   
      DO K = 1,p%NumBl
         
         CALL MeshCreate( BlankMesh          = y%BladeLn2Mesh(K)      &
                           , NNodes          = p%BldNodes+2           &
                           , IOS             = COMPONENT_OUTPUT       &
                           , TranslationDisp = .TRUE.                 &
                           , Orientation     = .TRUE.                 &
                           , RotationVel     = .TRUE.                 &
                           , TranslationVel  = .TRUE.                 &
                           , RotationAcc     = .TRUE.                 &
                           , TranslationAcc  = .TRUE.                 &
                           , ErrStat         = ErrStat2               &
                           , ErrMess         = ErrMsg2                )
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
      
         DO J = 1,p%BldNodes               
            CALL MeshPositionNode ( y%BladeLn2Mesh(K), J, u%BladePtLoads(K)%Position(:,J), ErrStat2, ErrMsg2, Orient=u%BladePtLoads(K)%RefOrientation(:,:,J) )
               CALL CheckError( ErrStat2, ErrMsg2 )
               IF (ErrStat >= AbortErrLev) RETURN
         END DO
            
            ! now add position/orientation of nodes for AD14 or AD15
         if (p%UseAD14) then     ! position/orientation of nodes for AeroDyn v14 or v15    
         
               ! Use orientation at p%BldNodes for the extra node at the blade tip
            CALL MeshPositionNode ( y%BladeLn2Mesh(K), p%BldNodes + 1, (/0.0_ReKi, 0.0_ReKi, p%BldFlexL /), ErrStat2, ErrMsg2, Orient=u%BladePtLoads(K)%RefOrientation(:,:,p%BldNodes) )
               CALL CheckError( ErrStat2, ErrMsg2 )
               IF (ErrStat >= AbortErrLev) RETURN
            
               ! Use orientation at node 1 for the blade root            
            CALL MeshPositionNode ( y%BladeLn2Mesh(K), p%BldNodes + 2, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi /), ErrStat2, ErrMsg2, Orient=u%BladePtLoads(K)%RefOrientation(:,:,1) )
               CALL CheckError( ErrStat2, ErrMsg2 )
               IF (ErrStat >= AbortErrLev) RETURN
               
         else
         
            ! position the nodes on the blade root and blade tip:
            DO J = 0,p%TipNode,p%TipNode
               if (j==0) then ! blade root
                  NodeNum = p%BldNodes + 2
               elseif (j==p%TipNode) then ! blade tip
                  NodeNum = p%BldNodes + 1
               end if
         
               Orientation(1,1) =     OtherState%CoordSys%n1(K,J,1)
               Orientation(2,1) =     OtherState%CoordSys%n2(K,J,1)
               Orientation(3,1) =     OtherState%CoordSys%n3(K,J,1)
               Orientation(1,2) = -1.*OtherState%CoordSys%n1(K,J,3)
               Orientation(2,2) = -1.*OtherState%CoordSys%n2(K,J,3)
               Orientation(3,2) = -1.*OtherState%CoordSys%n3(K,J,3)
               Orientation(1,3) =     OtherState%CoordSys%n1(K,J,2)
               Orientation(2,3) =     OtherState%CoordSys%n2(K,J,2)
               Orientation(3,3) =     OtherState%CoordSys%n3(K,J,2) 
               
                  ! Translational Displacement 
               position(1) =     OtherState%RtHS%rS (1,K,J)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
               position(2) = -1.*OtherState%RtHS%rS (3,K,J)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
               position(3) =     OtherState%RtHS%rS (2,K,J)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
               
               
               CALL MeshPositionNode ( y%BladeLn2Mesh(K), NodeNum, position, ErrStat2, ErrMsg2, Orient=Orientation )
                  CALL CheckError( ErrStat2, ErrMsg2 )
                  IF (ErrStat >= AbortErrLev) RETURN
                                    
            END DO ! nodes 
            
         end if ! position/orientation of nodes for AeroDyn v14 or v15
         
         ! create elements:      
         DO J = 2,p%TipNode !p%BldNodes + 1
            
            CALL MeshConstructElement ( Mesh      = y%BladeLn2Mesh(K)  &
                                       , Xelement = ELEMENT_LINE2      &
                                       , P1       = J-1                &   ! node1 number
                                       , P2       = J                  &   ! node2 number
                                       , ErrStat  = ErrStat2           &
                                       , ErrMess  = ErrMsg2            )
               CALL CheckError( ErrStat2, ErrMsg2 )
               IF (ErrStat >= AbortErrLev) RETURN
      
         END DO ! J (blade nodes)

            ! add the other extra element, connecting the first node on the blade:
         CALL MeshConstructElement ( Mesh      = y%BladeLn2Mesh(K)  &
                                    , Xelement = ELEMENT_LINE2      &
                                    , P1       = p%BldNodes + 2     &   ! node1 number (extra node at root)
                                    , P2       = 1                  &   ! node2 number (first node on blade)
                                    , ErrStat  = ErrStat2           &
                                    , ErrMess  = ErrMsg2            )         
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
      
         
            ! that's our entire mesh:
         CALL MeshCommit ( y%BladeLn2Mesh(K), ErrStat2, ErrMsg2 )   
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
   
      END DO
      
   END IF
   
   !.......................................................
   ! Create Point Mesh for Motions Output at Platform Reference Point:
   !.......................................................
      
   CALL MeshCopy ( SrcMesh  = u%PlatformPtMesh &
                 , DestMesh = y%PlatformPtMesh &
                 , CtrlCode = MESH_SIBLING     &
                 , IOS      = COMPONENT_OUTPUT &
                 , TranslationDisp = .TRUE.    &
                 , Orientation     = .TRUE.    &
                 , RotationVel     = .TRUE.    &
                 , TranslationVel  = .TRUE.    &
                 , RotationAcc     = .TRUE.    &
                 , TranslationAcc  = .TRUE.    &
                 , ErrStat  = ErrStat2         &
                 , ErrMess  = ErrMsg2          )  ! automatically sets    y%PlatformPtMesh%RemapFlag = .TRUE.
   
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   
   
   !.......................................................
   ! Create Line2 Mesh for Motions Output on Tower Line2 Mesh:
   !  first p%TwrNodes nodes are the same as the input TowerPtLoads mesh
   !.......................................................
      
   CALL MeshCreate( BlankMesh = y%TowerLn2Mesh           &
                    , IOS             = COMPONENT_INPUT  &
                    , NNodes          = p%TwrNodes + 2   &
                    , TranslationDisp = .TRUE.           &
                    , Orientation     = .TRUE.           &
                    , RotationVel     = .TRUE.           &
                    , TranslationVel  = .TRUE.           &  
                    , RotationAcc     = .TRUE.           &  
                    , TranslationAcc  = .TRUE.           &
                    , ErrStat         = ErrStat2         &
                    , ErrMess         = ErrMsg2          )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
   
      ! position the nodes on the tower:
      DO J = 1,p%TwrNodes      
         CALL MeshPositionNode ( y%TowerLn2Mesh, J, u%TowerPtLoads%Position(:,J), ErrStat2, ErrMsg2, &
                                 orient = u%TowerPtLoads%RefOrientation(:,:,J) )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF (ErrStat >= AbortErrLev) RETURN
      END DO

   ! for now, we're going to add two nodes, one at the beginning and the other at the end
   ! they're numbered this way so that I don't have to redo all the loops in the computational part of the code
   ! I am not going to use them in my input.
      CALL MeshPositionNode ( y%TowerLn2Mesh, p%TwrNodes + 1, (/0.0_ReKi, 0.0_ReKi, p%TowerHt /), ErrStat2, ErrMsg2 ) 
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
         
      CALL MeshPositionNode ( y%TowerLn2Mesh, p%TwrNodes + 2, (/0.0_ReKi, 0.0_ReKi, p%TowerBsHt  /), ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN

            
      ! create elements:
      
   IF ( p%TwrNodes < 2_IntKi ) THEN  ! if there are less than 2 nodes, we'll throw an error:
      CALL CheckError(ErrID_Fatal,"Tower Line2 Mesh cannot be created with less than 2 elements.")
      RETURN
   ELSE ! create line2 elements from the tower nodes:
      DO J = 2,p%TwrNodes+1  !the plus 1 includes one of the end nodes
         CALL MeshConstructElement ( Mesh      = y%TowerLn2Mesh     &
                                    , Xelement = ELEMENT_LINE2      &
                                    , P1       = J-1                &   ! node1 number
                                    , P2       = J                  &   ! node2 number
                                    , ErrStat  = ErrStat2           &
                                    , ErrMess  = ErrMsg2            )
         
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
      END DO
      
   ! add the other extra element, connecting the first node:
      
         CALL MeshConstructElement ( Mesh      = y%TowerLn2Mesh     &
                                    , Xelement = ELEMENT_LINE2      &
                                    , P1       = p%TwrNodes + 2     &   ! node1 number
                                    , P2       = 1                  &   ! node2 number
                                    , ErrStat  = ErrStat2           &
                                    , ErrMess  = ErrMsg2            )
         
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
                                          
   END IF   
   
      ! that's our entire mesh:
   CALL MeshCommit ( y%TowerLn2Mesh, ErrStat2, ErrMsg2 )   
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
                                                          
   !.......................................................
   ! Create Point Meshes for motions AeroDyn/BeamDyn needs:
   !.......................................................
   
   ! -------------- Hub -----------------------------------
      !BJJ: sibling of u%HubPtLoads
   CALL MeshCopy (     SrcMesh  = u%HubPtLoad             &
                     , DestMesh = y%HubPtMotion           &
                     , CtrlCode = MESH_SIBLING            &
                     , IOS      = COMPONENT_OUTPUT        &      
                     , TranslationDisp = .TRUE.           &
                     , Orientation     = .TRUE.           &
                     , RotationVel     = .TRUE.           &
                     ,ErrStat          = ErrStat2         &
                     ,ErrMess          = ErrMsg2          )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
      
   ! -------------- pseudo-Hub (for AD v14)  -----------------------------------
   CALL MeshCreate( BlankMesh          = y%HubPtMotion14  &
                     ,IOS              = COMPONENT_OUTPUT &
                     ,NNodes           = 1                &
                     , TranslationDisp = .TRUE.           &
                     , Orientation     = .TRUE.           &
                     , RotationVel     = .TRUE.           &
                     , ErrStat         = ErrStat2         &
                     , ErrMess         = ErrMsg2          )      
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
      
      ! pseudo-Hub position and orientation (relative here as before, but should not be)
      
   CALL MeshPositionNode ( y%HubPtMotion14, 1, (/0.0_ReKi, 0.0_ReKi, p%HubHt /), ErrStat, ErrMsg ) !orientation is identity by default
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
            
   CALL CommitPointMesh( y%HubPtMotion14 )
      IF (ErrStat >= AbortErrLev) RETURN
 
      
   ! -------------- Blade Roots -----------------------------------
   ALLOCATE( y%BladeRootMotion(p%NumBl), Stat=ErrStat2 )
   IF ( ErrStat2 /= 0 ) THEN
      CALL CheckError( ErrID_Fatal, 'ED: Could not allocate space for y%BladeRootMotions{p%NumBl}' )
      RETURN
   END IF
                  
      
   DO k=1,p%NumBl      
      CALL MeshCreate( BlankMesh       = y%BladeRootMotion(k)   &
                     ,IOS              = COMPONENT_OUTPUT       &
                     ,NNodes           = 1                      &
                     ,TranslationDisp  = .TRUE.                 & 
                     ,Orientation      = .TRUE.                 & 
                     ,TranslationVel   = .TRUE.                 & 
                     ,TranslationAcc   = .TRUE.                 & 
                     ,RotationVel      = .TRUE.                 & 
                     ,RotationAcc      = .TRUE.                 & 
                     ,ErrStat          = ErrStat2               &
                     ,ErrMess          = ErrMsg2                )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN            
   END DO
   
      
   CALL MeshCreate( BlankMesh          = y%BladeRootMotion14    &
                     ,IOS              = COMPONENT_OUTPUT       &
                     ,NNodes           = p%NumBl                &
                     , Orientation     = .TRUE.                 &
                     ,ErrStat          = ErrStat2               &
                     ,ErrMess          = ErrMsg2                )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN

   DO K=1,p%NumBl      
      
      Orientation(1,1) =               p%CosPreC(K)
      Orientation(2,1) =  0.0_R8Ki
      Orientation(3,1) =  1.0_R8Ki *   p%SinPreC(K)

      Orientation(1,2) =  0.0_R8Ki
      Orientation(2,2) =  1.0_R8Ki
      Orientation(3,2) =  0.0_R8Ki

      Orientation(1,3) = -1.0_R8Ki *    p%SinPreC(K)
      Orientation(2,3) =  0.0_R8Ki
      Orientation(3,3) =                p%CosPreC(K)
                  
      Position(1) = p%HubRad*p%SinPreC(K)
      Position(2) = 0.0_ReKi
      Position(3) = p%HubRad*p%CosPreC(K)      
      
      CALL MeshPositionNode ( y%BladeRootMotion14, K, Position, &
                            ErrStat, ErrMsg, Orient=Orientation ) 
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
                  
               
      position(1) =     OtherState%RtHS%rS (1,K,0)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
      position(2) = -1.*OtherState%RtHS%rS (3,K,0)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
      position(3) =     OtherState%RtHS%rS (2,K,0)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
      
      
      Orientation(1,1) =     OtherState%CoordSys%j1(K,1)
      Orientation(2,1) =     OtherState%CoordSys%j2(K,1)
      Orientation(3,1) =     OtherState%CoordSys%j3(K,1)
      Orientation(1,2) = -1.*OtherState%CoordSys%j1(K,3)
      Orientation(2,2) = -1.*OtherState%CoordSys%j2(K,3)
      Orientation(3,2) = -1.*OtherState%CoordSys%j3(K,3)
      Orientation(1,3) =     OtherState%CoordSys%j1(K,2)
      Orientation(2,3) =     OtherState%CoordSys%j2(K,2)
      Orientation(3,3) =     OtherState%CoordSys%j3(K,2)
      
      CALL MeshPositionNode ( y%BladeRootMotion(K), 1, Position, &
                            ErrStat, ErrMsg, Orient=Orientation ) 
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN
         
   END DO
                     
   CALL CommitPointMesh( y%BladeRootMotion14 )
      IF (ErrStat >= AbortErrLev) RETURN
   
   DO k=1,p%NumBl      
      CALL CommitPointMesh( y%BladeRootMotion(K) )
         IF (ErrStat >= AbortErrLev) RETURN
   END DO
   
      
   ! -------------- Rotor Furl -----------------------------------
   CALL MeshCreate( BlankMesh          = y%RotorFurlMotion14    &
                     ,IOS              = COMPONENT_OUTPUT       &
                     ,NNodes           = 1                      &
                     , TranslationDisp = .TRUE.                 &
                     , Orientation     = .TRUE.                 &
                     , RotationVel     = .TRUE.                 &
                     ,ErrStat          = ErrStat2               &
                     ,ErrMess          = ErrMsg2                )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN

!bjj: FIX THIS>>>>     
!call wrscr(newline//'fix RotorFurlMotion initialization')
   CALL MeshPositionNode ( y%RotorFurlMotion14, 1, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi /), ErrStat, ErrMsg ) !orientation is identity by default
!<<<<<FIX THIS
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
      
   CALL CommitPointMesh( y%RotorFurlMotion14 )
      IF (ErrStat >= AbortErrLev) RETURN      
      
   ! -------------- Nacelle -----------------------------------      
   CALL MeshCopy ( SrcMesh  = u%NacelleLoads   &
                 , DestMesh = y%NacelleMotion  &
                 , CtrlCode = MESH_SIBLING     &
                 , IOS      = COMPONENT_OUTPUT &
                 , TranslationDisp = .TRUE.    &
                 , Orientation     = .TRUE.    &
                 , TranslationVel  = .TRUE.    &
                 , RotationVel     = .TRUE.    &
                 , TranslationAcc  = .TRUE.    &
                 , RotationAcc     = .TRUE.    &   
                 , ErrStat         = ErrStat2  &
                 , ErrMess         = ErrMsg2   )      ! automatically sets    y%NacelleMotion%RemapFlag   = .TRUE.
   
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
      
     
   ! -------------- Tower Base-----------------------------------
   CALL MeshCreate( BlankMesh          = y%TowerBaseMotion14    &
                     ,IOS              = COMPONENT_OUTPUT       &
                     ,NNodes           = 1                      &
                     , TranslationDisp = .TRUE.                 &
                     , RotationVel     = .TRUE.                 &
                     ,ErrStat          = ErrStat2               &
                     ,ErrMess          = ErrMsg2                )
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN

!bjj: FIX THIS>>>>     
!call wrscr(newline//'fix TowerBaseMotion14 initialization')
   CALL MeshPositionNode ( y%TowerBaseMotion14, 1, (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi /), ErrStat, ErrMsg ) !orientation is identity by default
!<<<<<FIX THIS
      CALL CheckError(ErrStat2,ErrMsg2)
      IF (ErrStat >= AbortErrLev) RETURN
      
   CALL CommitPointMesh( y%TowerBaseMotion14 )
      IF (ErrStat >= AbortErrLev) RETURN      
      
      
CONTAINS
   !...............................................................................................................................
   SUBROUTINE CommitPointMesh(NewMesh)
      ! This routine makes every node a point element and then commits the mesh.
      
      TYPE(MeshType), INTENT(INOUT)  :: NewMesh  
      
      INTEGER(IntKi) :: Node
      
      DO Node = 1,NewMesh%Nnodes
         
            ! create an element from this point      
         CALL MeshConstructElement ( Mesh = NewMesh                 &
                                    , Xelement = ELEMENT_POINT      &
                                    , P1       = Node               &   ! node number
                                    , ErrStat  = ErrStat            &
                                    , ErrMess  = ErrMsg             )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF (ErrStat >= AbortErrLev) RETURN

      END DO
      
         ! that's our entire mesh:
      CALL MeshCommit ( NewMesh, ErrStat2, ErrMsg2 )   
         CALL CheckError(ErrStat2,ErrMsg2)
         IF (ErrStat >= AbortErrLev) RETURN                

   END SUBROUTINE CommitPointMesh
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................
         ! Passed arguments
         
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)      
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_AllocOutput:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................

      END IF

   END SUBROUTINE CheckError 
   !...............................................................................................................................
END SUBROUTINE ED_AllocOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_u( u, p, x, InputFileData, OtherState, ErrStat, ErrMsg )
! This routine allocates the arrays stored in the ED_InputType data structure (u), based on the parameters (p). 
! The routine assumes that the arrays are not currently allocated (It will produce a fatal error otherwise.) It also initializes the inputs
!..................................................................................................................................

   TYPE(ED_InputType),           INTENT(INOUT)  :: u                 ! Inputs to be allocated
   TYPE(ED_ParameterType),       INTENT(IN   )  :: p                 ! Parameters
   TYPE(ED_ContinuousStateType), INTENT(IN   )  :: x                 ! Continuous states
   TYPE(ED_InputFile),           INTENT(IN   )  :: InputFileData     ! Data stored in the module's input file
   TYPE(ED_OtherStateType),      INTENT(INOUT)  :: OtherState        ! Other states
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat           ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg            ! Error message if ErrStat /= ErrID_None
      
   
   ! local variables
   REAL(R8Ki)                                   :: Orientation(3,3)  ! reference orientation matrix
   REAL(ReKi)                                   :: Position(3)       ! position vector
   TYPE(ED_ContinuousStateType)                 :: x_tmp             ! continuous states (set to 0)
   INTEGER(IntKi)                               :: J, K              ! loop counters
   INTEGER(IntKi)                               :: NodeNum           ! number of current blade node
   INTEGER(IntKi)                               :: ErrStat2          ! The error identifier (ErrStat)
   CHARACTER(ErrMsgLen)                         :: ErrMsg2           ! The error message (ErrMsg)
   CHARACTER(*), PARAMETER                      :: RoutineName = 'Init_u'
   
      ! initialize variables:
      
   ErrStat = ErrID_None
   ErrMsg  = ""

   !.......................................................
   ! allocate the u%BlPitchCom array    
   !.......................................................

   CALL AllocAry( u%BlPitchCom, p%NumBl, 'BlPitchCom', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   ! will initialize u%BlPitchCom later, after getting undisplaced positions    
   
   !.......................................................
   ! we're going to calculate the non-displaced positions of
   ! several variables so we can set up meshes properly later.
   ! want inputs and states initialized to 0 first.
   !.......................................................
   CALL ED_CopyContState( x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      x_tmp%qt  = 0.0_ReKi
      x_tmp%qdt = 0.0_ReKi
      x_tmp%QT (DOF_GeAz) = - p%AzimB1Up - Piby2
         CALL Zero2TwoPi( x_tmp%QT (DOF_GeAz) )

      u%BlPitchCom = 0.0_ReKi
      
      ! set the coordinate system variables:
   CALL SetCoordSy( -p%DT, OtherState%CoordSys, OtherState%RtHS, u%BlPitchCom, p, x_tmp, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   CALL CalculatePositions( p, x_tmp, OtherState%CoordSys, OtherState%RtHS ) ! calculate positions
   
   CALL ED_DestroyContState(x_tmp, ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !.......................................................
   ! initialize the u%BlPitchCom array    
   !.......................................................
      ! was allocated above to call SetCoordSy and CalculatePositions with undisplaced values
   u%BlPitchCom = InputFileData%BlPitch(1:p%NumBl)

   
   !.......................................................
   ! Create Line2 Meshes for loads input on blades:
   !.......................................................
      
   IF (.not. p%BD4Blades) THEN
      ALLOCATE( u%BladePtLoads(p%NumBl), STAT=ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         CALL SetErrStat( ErrID_Fatal, "Could not allocate u%BladePtLoads", ErrStat, ErrMsg, RoutineName )
         CALL Cleanup()
         RETURN
      END IF
   
   
      DO K=1,p%NumBl
      
         CALL MeshCreate( BlankMesh         = u%BladePtLoads(K)      &
                           ,IOS             = COMPONENT_INPUT        &
                           ,NNodes          = p%BldNodes             &
                           ,Force           = .TRUE.                 &
                           ,Moment          = .TRUE.                 &
                           ,ErrStat         = ErrStat2               &
                           ,ErrMess         = ErrMsg2                )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
      
         if (p%UseAD14) then
            ! position the nodes on the blades:
            DO J = 1,p%BldNodes
         
               NodeNum = J
         
               Orientation(1,1) =  p%CAeroTwst(J)
               Orientation(2,1) =  p%SAeroTwst(J)
               Orientation(3,1) =  0.0_ReKi

               Orientation(1,2) = -p%SAeroTwst(J)
               Orientation(2,2) =  p%CAeroTwst(J)
               Orientation(3,2) =  0.0_ReKi

               Orientation(1,3) =  0.0_ReKi
               Orientation(2,3) =  0.0_ReKi
               Orientation(3,3) =  1.0_ReKi
                           
               CALL MeshPositionNode ( u%BladePtLoads(K), NodeNum, (/0.0_ReKi, 0.0_ReKi, p%RNodes(J) /), ErrStat2, ErrMsg2, Orient=Orientation )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  IF (ErrStat >= AbortErrLev) THEN
                     CALL Cleanup()
                     RETURN
                  END IF
                                                               
            END DO ! nodes  
         else
            ! position the nodes on the blades:
            DO J = 1,p%BldNodes
               NodeNum = J
         
               Orientation(1,1) =     OtherState%CoordSys%n1(K,J,1)
               Orientation(2,1) =     OtherState%CoordSys%n2(K,J,1)
               Orientation(3,1) =     OtherState%CoordSys%n3(K,J,1)
               Orientation(1,2) = -1.*OtherState%CoordSys%n1(K,J,3)
               Orientation(2,2) = -1.*OtherState%CoordSys%n2(K,J,3)
               Orientation(3,2) = -1.*OtherState%CoordSys%n3(K,J,3)
               Orientation(1,3) =     OtherState%CoordSys%n1(K,J,2)
               Orientation(2,3) =     OtherState%CoordSys%n2(K,J,2)
               Orientation(3,3) =     OtherState%CoordSys%n3(K,J,2) 
               
                  ! Translational Displacement 
               position(1) =     OtherState%RtHS%rS (1,K,J)                ! = the distance from the undeflected tower centerline to the current blade node in the xi ( z1) direction
               position(2) = -1.*OtherState%RtHS%rS (3,K,J)                ! = the distance from the undeflected tower centerline to the current blade node in the yi (-z3) direction
               position(3) =     OtherState%RtHS%rS (2,K,J)  + p%PtfmRefzt ! = the distance from the nominal tower base position (i.e., the undeflected position of the tower base) to the current blade node in the zi ( z2) direction
               
               
               CALL MeshPositionNode ( u%BladePtLoads(K), NodeNum, position, ErrStat2, ErrMsg2, Orient=Orientation )
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  IF (ErrStat >= AbortErrLev) THEN
                     CALL Cleanup()
                     RETURN
                  END IF                           
                                    
            END DO ! nodes              
         end if ! position/orientation of nodes for AeroDyn v14 or v15
         
         ! create elements:      
         DO J = 1,p%BldNodes !p%BldNodes + 1
            
            CALL MeshConstructElement ( Mesh      = u%BladePtLoads(K)  &
                                       , Xelement = ELEMENT_POINT      &
                                       , P1       = J                  &   ! node1 number
                                       , ErrStat  = ErrStat2           &
                                       , ErrMess  = ErrMsg2            )
         
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               END IF
      
         END DO ! J (blade nodes)

            ! that's our entire mesh:
         CALL MeshCommit ( u%BladePtLoads(K), ErrStat2, ErrMsg2 )   
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF

   
            ! initialize it
         u%BladePtLoads(K)%Moment   = 0.0_ReKi
         u%BladePtLoads(K)%Force    = 0.0_ReKi         
         
                     
      END DO ! blades
   END IF ! p%BD4Blades
   
                     
      !.......................................................
      ! Create Point Mesh for loads input at hub point (from BeamDyn):
      !....................................................... 
    
   CALL MeshCreate( BlankMesh      = u%HubPtLoad            &
                  ,IOS             = COMPONENT_INPUT        &
                  ,NNodes          = 1                      &
                  ,Force           = .TRUE.                 &
                  ,Moment          = .TRUE.                 &
                  ,ErrStat         = ErrStat2               &
                  ,ErrMess         = ErrMsg2                )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

   ! place single node at hub; position affects mapping/coupling with other modules      
   Position(1)  =     OtherState%RtHS%rQ(1)
   Position(2)  = -1.*OtherState%RtHS%rQ(3)
   Position(3)  =     OtherState%RtHS%rQ(2) + p%PtfmRefzt
   
   Orientation(1,1) =     OtherState%CoordSys%g1(1)
   Orientation(2,1) =     OtherState%CoordSys%g2(1)
   Orientation(3,1) =     OtherState%CoordSys%g3(1)
   Orientation(1,2) = -1.*OtherState%CoordSys%g1(3)
   Orientation(2,2) = -1.*OtherState%CoordSys%g2(3)
   Orientation(3,2) = -1.*OtherState%CoordSys%g3(3)
   Orientation(1,3) =     OtherState%CoordSys%g1(2)
   Orientation(2,3) =     OtherState%CoordSys%g2(2)
   Orientation(3,3) =     OtherState%CoordSys%g3(2) 
      
   CALL MeshPositionNode ( u%HubPtLoad, 1, Position, ErrStat2, ErrMsg2, orient=Orientation )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! create an element from this point      
      CALL MeshConstructElement ( Mesh = u%HubPtLoad       &
                           , Xelement = ELEMENT_POINT      &
                           , P1       = 1                  &   ! node number
                           , ErrStat  = ErrStat2           &
                           , ErrMess  = ErrMsg2            )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF

         ! that's our entire mesh:
      CALL MeshCommit ( u%HubPtLoad, ErrStat2, ErrMsg2 )   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF

         ! initailize it
      u%HubPtLoad%Moment      = 0.0_ReKi
      u%HubPtLoad%Force       = 0.0_ReKi
         
                     
   !.......................................................
   ! Create Point Mesh for loads input at Platform Reference Point:
   !.......................................................
      
   CALL MeshCreate( BlankMesh         = u%PlatformPtMesh       &
                     ,IOS             = COMPONENT_INPUT        &
                     ,NNodes          = 1                      &
                     ,Force           = .TRUE.                 &
                     ,Moment          = .TRUE.                 &
                     ,ErrStat         = ErrStat2               &
                     ,ErrMess         = ErrMsg2                )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! place single node at platform reference point; position affects mapping/coupling with other modules
   CALL MeshPositionNode ( u%PlatformPtMesh, 1, (/0.0_ReKi, 0.0_ReKi, p%PtfmRefzt /), ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! create an element from this point      
   CALL MeshConstructElement ( Mesh = u%PlatformPtMesh        &
                              , Xelement = ELEMENT_POINT      &
                              , P1       = 1                  &   ! node number
                              , ErrStat  = ErrStat2           &
                              , ErrMess  = ErrMsg2            )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! that's our entire mesh:
   CALL MeshCommit ( u%PlatformPtMesh, ErrStat2, ErrMsg2 )   
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   
      ! initialize fields
   u%PlatformPtMesh%Moment = 0.0_ReKi
   u%PlatformPtMesh%Force  = 0.0_ReKi
      
   !.......................................................
   ! Create Point Mesh for loads input at nacelle:
   !.......................................................
         
   CALL MeshCreate( BlankMesh          = u%NacelleLoads      &
                     ,IOS              = COMPONENT_OUTPUT    &
                     ,NNodes           = 1                   &
                     ,Force            = .TRUE.              &
                     ,Moment           = .TRUE.              &   
                     ,ErrStat          = ErrStat2            &
                     ,ErrMess          = ErrMsg2             )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL MeshPositionNode ( u%NacelleLoads,  1, (/0.0_ReKi, 0.0_ReKi, p%TowerHt /), ErrStat2, ErrMsg2 ) ! orientation is identity by default
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! create an element from this point      
   CALL MeshConstructElement ( Mesh = u%NacelleLoads          &
                              , Xelement = ELEMENT_POINT      &
                              , P1       = 1                  &   ! node number
                              , ErrStat  = ErrStat2           &
                              , ErrMess  = ErrMsg2            )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      
   CALL MeshCommit ( u%NacelleLoads, ErrStat2, ErrMsg2 )   
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
            
      ! initialize fields
   u%NacelleLoads%Force    = 0.0_ReKi
   u%NacelleLoads%Moment   = 0.0_ReKi
      
   !.......................................................
   ! Create u%TwrAddedMass for loads input on tower:
   ! SHOULD REMOVE EVENTUALLY
   !.......................................................
         
   CALL AllocAry( u%TwrAddedMass,  6_IntKi, 6_IntKi, p%TwrNodes,   'TwrAddedMass',    ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! initialize it
   u%TwrAddedMass          = 0.0_ReKi  
      
   !.......................................................
   ! Create point Mesh for lumped load input on tower:
   ! note that this does not contain the end points
   !.......................................................
      
   CALL MeshCreate( BlankMesh      = u%TowerPtLoads         &
                     ,IOS          = COMPONENT_INPUT        &
                     ,NNodes       = p%TwrNodes             &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
   
      ! position the nodes on the tower:
   DO J = 1,p%TwrNodes      
      CALL MeshPositionNode ( u%TowerPtLoads, J, (/0.0_ReKi, 0.0_ReKi, p%HNodes(J) + p%TowerBsHt /), ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
   END DO
   
      ! create elements:      
   DO J = 1,p%TwrNodes
      CALL MeshConstructElement ( Mesh      = u%TowerPtLoads     &
                                 , Xelement = ELEMENT_POINT      &
                                 , P1       = J                  &   ! node1 number
                                 , ErrStat  = ErrStat2           &
                                 , ErrMess  = ErrMsg2            )
         
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
   END DO
      
   
      ! that's our entire mesh:
   CALL MeshCommit ( u%TowerPtLoads, ErrStat2, ErrMsg2 )   
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! initialize fields
   u%TowerPtLoads%Moment   = 0.0_ReKi
   u%TowerPtLoads%Force    = 0.0_ReKi   
         
   !.......................................................
   ! initialize all remaining inputs (non-allocatable):
   !.......................................................
         
   u%PtfmAddedMass         = 0.0_ReKi      
   u%YawMom                = 0.0_ReKi
   u%GenTrq                = 0.0_ReKi
   u%HSSBrTrqC             = 0.0_ReKi      
         
   CALL Cleanup()
   
CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up if the error for returning to calling routine
   !...............................................................................................................................      
         !.........................................................................................................................
         ! close files, deallocate local arrays
         !.........................................................................................................................
         CALL ED_DestroyContState( x_tmp, ErrStat2, ErrMsg2 )
         
   END SUBROUTINE Cleanup   
            
END SUBROUTINE Init_u
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!   Define constants k1, k2, k3, and k4 as 
!        k1 = dt * f(t        , x_t        )
!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!        k4 = dt * f(t + dt   , x_t + k3   ).
!   Then the continuous states at t = t + dt are
!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!
! For details, see:
! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!   Runge-Kutta." §16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!   Cambridge University Press, pp. 704-716, 1992.
!
!..................................................................................................................................

      REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),               INTENT(IN   )  :: n           ! time step number
      TYPE(ED_InputType),           INTENT(INOUT)  :: u(:)        ! Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                   INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(ED_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(ED_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(ED_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ED_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ED_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
         
      TYPE(ED_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(ED_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(ED_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(ED_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(ED_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(ED_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(ED_InputType)                           :: u_interp    ! interpolated value of inputs 

      INTEGER(IntKi)                               :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                         :: ErrMsg2     ! local error message (ErrMsg)
      
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL ED_CopyContState( x, k1, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL ED_CopyContState( x, k2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL ED_CopyContState( x, k3, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL ED_CopyContState( x, k4,    MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL ED_CopyContState( x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN


      CALL ED_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
                     
      ! interpolate u to find u_interp = u(t)
      CALL ED_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
!      HSSBrTrq_at_t = u_interp%HSSBrTrqC
!      OtherState%HSSBrTrqC = SIGN( u_interp%HSSBrTrqC, x%QDT(DOF_GeAz) )         
!      OtherState%HSSBrTrq  = OtherState%HSSBrTrqC         

      ! find xdot at t
      CALL ED_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k1%qt  = p%dt * xdot%qt
      k1%qdt = p%dt * xdot%qdt
  
      x_tmp%qt  = x%qt  + 0.5 * k1%qt
      x_tmp%qdt = x%qdt + 0.5 * k1%qdt

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL ED_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
!      u_interp%HSSBrTrqC = max(0.0_ReKi, min(u_interp%HSSBrTrqC, HSSBrTrq_at_t )) ! hack for extrapolation of limits       
!      OtherState%HSSBrTrqC = SIGN( u_interp%HSSBrTrqC, x_tmp%QDT(DOF_GeAz) )         
!      OtherState%HSSBrTrq  = OtherState%HSSBrTrqC         

      ! find xdot at t + dt/2
      CALL ED_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k2%qt  = p%dt * xdot%qt
      k2%qdt = p%dt * xdot%qdt

      x_tmp%qt  = x%qt  + 0.5 * k2%qt
      x_tmp%qdt = x%qdt + 0.5 * k2%qdt

      ! find xdot at t + dt/2
!      u_interp%HSSBrTrqC = max(0.0_ReKi, min(u_interp%HSSBrTrqC, HSSBrTrq_at_t )) ! hack for extrapolation of limits       
!      OtherState%HSSBrTrqC = SIGN( u_interp%HSSBrTrqC, x_tmp%QDT(DOF_GeAz) )         
!      OtherState%HSSBrTrq  = OtherState%HSSBrTrqC         
      CALL ED_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k3%qt  = p%dt * xdot%qt
      k3%qdt = p%dt * xdot%qdt

      x_tmp%qt  = x%qt  + k3%qt
      x_tmp%qdt = x%qdt + k3%qdt

      ! interpolate u to find u_interp = u(t + dt)
      CALL ED_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN
!      u_interp%HSSBrTrqC = max(0.0_ReKi, min(u_interp%HSSBrTrqC, HSSBrTrq_at_t )) ! hack for extrapolation of limits       
!      OtherState%HSSBrTrqC = SIGN( u_interp%HSSBrTrqC, x_tmp%QDT(DOF_GeAz) )         
!      OtherState%HSSBrTrq  = OtherState%HSSBrTrqC         

      ! find xdot at t + dt
      CALL ED_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      k4%qt  = p%dt * xdot%qt
      k4%qdt = p%dt * xdot%qdt

      x%qt  = x%qt  +  ( k1%qt  + 2. * k2%qt  + 2. * k3%qt  + k4%qt  ) / 6.      
      x%qdt = x%qdt +  ( k1%qdt + 2. * k2%qdt + 2. * k3%qdt + k4%qdt ) / 6.      

         ! clean up local variables:
      CALL ExitThisRoutine(  )
         
CONTAINS      
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)
   
   
      CALL ED_DestroyContState( xdot,     ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )

      CALL ED_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
         
   END SUBROUTINE ExitThisRoutine      
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_RK4:'//TRIM(Msg)         
         ErrStat = MAX(ErrStat,ErrID)
         
         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
                  
         
      END IF

   END SUBROUTINE CheckError                    
      
END SUBROUTINE ED_RK4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
! equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(ED_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(ED_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ED_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(ED_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ED_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ED_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! local variables
      TYPE(ED_InputType)                             :: u_interp
      TYPE(ED_ContinuousStateType)                   :: xdot
         
      INTEGER(IntKi)                                 :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! local error message (ErrMsg)


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      
      if (OtherState%n .lt. n) then

         OtherState%n = n
            
         ! Update IC() index so IC(1) is the location of xdot values at n.
         ! (this allows us to shift the indices into the array, not copy all of the values)
         OtherState%IC = CSHIFT( OtherState%IC, -1 ) ! circular shift of all values to the right
            
      elseif (OtherState%n .gt. n) then
 
         CALL CheckError( ErrID_Fatal, ' Backing up in time is not supported with a multistep method.')
         RETURN

      endif        
      
      
      ! Allocate the input arrays
      CALL ED_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      
      ! need xdot at t
      CALL ED_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN                  
      IF (EqualRealNos( x%qdt(DOF_GeAz) ,0.0_R8Ki ) ) THEN
         OtherState%HSSBrTrqC = u_interp%HSSBrTrqC
      ELSE
         OtherState%HSSBrTrqC  = SIGN( u_interp%HSSBrTrqC, real(x%qdt(DOF_GeAz),ReKi) ) ! hack for HSS brake (need correct sign)
      END IF
      OtherState%HSSBrTrq   = OtherState%HSSBrTrqC
      OtherState%SgnPrvLSTQ = OtherState%SgnLSTQ(OtherState%IC(2))
      
      CALL ED_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         
         CALL ED_CopyContState(xdot, OtherState%xdot ( OtherState%IC(1) ), MESH_NEWCOPY, ErrStat2, ErrMsg2)
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

                                                    
      if (n .le. 2) then
                                               
         CALL ED_RK4(t, n, u, utimes, p, x, xd, z, OtherState, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

      else
         
         x%qt  = x%qt  + p%DT24 * ( 55.*OtherState%xdot(OtherState%IC(1))%qt  - 59.*OtherState%xdot(OtherState%IC(2))%qt   &
                                  + 37.*OtherState%xdot(OtherState%IC(3))%qt   - 9.*OtherState%xdot(OtherState%IC(4))%qt )

         x%qdt = x%qdt + p%DT24 * ( 55.*OtherState%xdot(OtherState%IC(1))%qdt - 59.*OtherState%xdot(OtherState%IC(2))%qdt  &
                                  + 37.*OtherState%xdot(OtherState%IC(3))%qdt  - 9.*OtherState%xdot(OtherState%IC(4))%qdt )
         
         
            ! Make sure the HSS brake will not reverse the direction of the HSS
            !   for the next time step.  Do this by computing the predicted value
            !   of x%qt(); QD(DOF_GeAz,IC(NMX)) as will be done during the next time step.
            ! Only do this after the first few time steps since it doesn't work
            !   for the Runga-Kutta integration scheme.
   
         
         CALL FixHSSBrTq ( 'P', p, x, OtherState, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
            
      endif
            
      OtherState%SgnPrvLSTQ = SignLSSTrq(p, OtherState)   
      OtherState%SgnLSTQ(OtherState%IC(1)) = OtherState%SgnPrvLSTQ 
      
      
         ! clean up local variables:
      CALL ExitThisRoutine()
      
CONTAINS      
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)
   
   
      CALL ED_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( xdot,     ErrStat2, ErrMsg3 )
      
   END SUBROUTINE ExitThisRoutine    
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_AB4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
         
      END IF

   END SUBROUTINE CheckError            
         
END SUBROUTINE ED_AB4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ED_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
! differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   Adams-Bashforth Predictor:
!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!   Adams-Moulton Corrector:
!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(ED_InputType),             INTENT(INOUT)  :: u(:)        ! Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(ED_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ED_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(ED_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ED_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ED_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(ED_InputType)                             :: u_interp    ! Inputs at t
      TYPE(ED_ContinuousStateType)                   :: x_pred      ! Continuous states at t
      TYPE(ED_ContinuousStateType)                   :: xdot_pred   ! Derivative of continuous states at t

      INTEGER(IntKi)                                 :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                           :: ErrMsg2     ! local error message (ErrMsg)
      
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
         ! predict:

      CALL ED_CopyContState(x, x_pred, MESH_NEWCOPY, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL ED_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

         
      if (n .gt. 2_IntKi) then
         
            ! correct:
         
            ! allocate the arrays in u_interp
         CALL ED_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         
         CALL ED_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
            
         u_interp%HSSBrTrqC = max(0.0_ReKi, min(u_interp%HSSBrTrqC, ABS( OtherState%HSSBrTrqC) )) ! hack for extrapolation of limits  (OtherState%HSSBrTrqC is HSSBrTrqC at t)     
         IF (EqualRealNos( x_pred%qdt(DOF_GeAz) ,0.0_R8Ki ) ) THEN
            OtherState%HSSBrTrqC = u_interp%HSSBrTrqC
         ELSE
            OtherState%HSSBrTrqC  = SIGN( u_interp%HSSBrTrqC, real(x_pred%qdt(DOF_GeAz),ReKi) ) ! hack for HSS brake (need correct sign)
         END IF
         OtherState%HSSBrTrq  = OtherState%HSSBrTrqC

         CALL ED_CalcContStateDeriv(t + p%dt, u_interp, p, x_pred, xd, z, OtherState, xdot_pred, ErrStat2, ErrMsg2 )
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN

         
         x%qt  = x%qt  + p%DT24 * ( 9. * xdot_pred%qt +  19. * OtherState%xdot(OtherState%IC(1))%qt &
                                                        - 5. * OtherState%xdot(OtherState%IC(2))%qt &
                                                        + 1. * OtherState%xdot(OtherState%IC(3))%qt )

         x%qdt = x%qdt + p%DT24 * ( 9. * xdot_pred%qdt + 19. * OtherState%xdot(OtherState%IC(1))%qdt &
                                                       -  5. * OtherState%xdot(OtherState%IC(2))%qdt &
                                                       +  1. * OtherState%xdot(OtherState%IC(3))%qdt )
         
                  
          ! Make sure the HSS brake has not reversed the direction of the HSS:
         
         CALL FixHSSBrTq ( 'C', p, x, OtherState, ErrStat2, ErrMsg2 )      
            CALL CheckError(ErrStat2,ErrMsg2)
            IF ( ErrStat >= AbortErrLev ) RETURN
         OtherState%SgnPrvLSTQ = SignLSSTrq(p, OtherState)
         OtherState%SgnLSTQ(OtherState%IC(1)) = OtherState%SgnPrvLSTQ 
                                                                    
      else

         x%qt  = x_pred%qt
         x%qdt = x_pred%qdt

      endif
      
      
         ! clean up local variables:
      CALL ExitThisRoutine()
      
CONTAINS      
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)
   
   
      CALL ED_DestroyContState( xdot_pred,  ErrStat3, ErrMsg3 )
      CALL ED_DestroyContState( x_pred,     ErrStat3, ErrMsg3 )
      CALL ED_DestroyInput(     u_interp,   ErrStat3, ErrMsg3 )               
      
   END SUBROUTINE ExitThisRoutine    
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'ED_ABM4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )                  
         
      END IF

   END SUBROUTINE CheckError                 

END SUBROUTINE ED_ABM4
!----------------------------------------------------------------------------------------------------------------------------------

SUBROUTINE ED_PrintSum( p, OtherState, GenerateAdamsModel, ErrStat, ErrMsg )
! This routine generates the summary file, which contains a regurgitation of  the input data and interpolated flexible body data.

      ! passed variables
   TYPE(ED_ParameterType),    INTENT(IN)  :: p                                    ! Parameters of the structural dynamics module
   TYPE(ED_OtherStateType),   INTENT(IN)  :: OtherState                           ! Other/optimization states of the structural dynamics module 
   LOGICAL,                   INTENT(IN)  :: GenerateAdamsModel                   ! Logical to determine if Adams inputs were read/should be summarized
   INTEGER(IntKi),            INTENT(OUT) :: ErrStat
   CHARACTER(*),              INTENT(OUT) :: ErrMsg


      ! Local variables.

   INTEGER(IntKi)               :: I                                               ! Index for the nodes.
   INTEGER(IntKi)               :: K                                               ! Generic index (also for the blade number).
   INTEGER(IntKi)               :: UnSu                                            ! I/O unit number for the summary output file

   CHARACTER(*), PARAMETER      :: Fmt1      = "(34X,3(6X,'Blade',I2,:))"          ! Format for outputting blade headings.
   CHARACTER(*), PARAMETER      :: Fmt2      = "(34X,3(6X,A,:))"                   ! Format for outputting blade headings.
   CHARACTER(*), PARAMETER      :: FmtDat    = '(A,T35,3(:,F13.3))'                ! Format for outputting mass and modal data.
   CHARACTER(*), PARAMETER      :: FmtDatT   = '(A,T35,1(:,F13.8))'                ! Format for outputting time steps.
   CHARACTER(100)               :: RotorType                                       ! Text description of rotor.

   CHARACTER(30)                :: OutPFmt                                         ! Format to print list of selected output channels to summary file
   CHARACTER(10)                :: DOFEnabled                                      ! String to say if a DOF is enabled or disabled

   ! Open the summary file and give it a heading.
   
   CALL GetNewUnit( UnSu, ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   CALL OpenFOutFile ( UnSu, TRIM( p%RootName )//'.sum', ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN

   
      ! Heading:
   WRITE (UnSu,'(/,A)')  'This summary information was generated by '//TRIM( GetNVD(ED_Ver) )// &
                         ' on '//CurDate()//' at '//CurTime()//'.'
   !WRITE (UnSu,'(//,1X,A)')  TRIM( p%FTitle )


   !..................................
   ! Turbine features.
   !..................................

   WRITE (UnSu,'(//,A,/)')  'Turbine features:'

   IF ( p%OverHang > 0.0 )  THEN
      RotorType = 'Downwind,'
   ELSE
      RotorType = 'Upwind,'
   ENDIF
   IF ( p%NumBl == 2 )  THEN
      RotorType = TRIM(RotorType)//' two-bladed rotor'
      IF ( p%DOF_Flag(DOF_Teet) )  THEN
         RotorType = TRIM(RotorType)//' with teetering hub.'
      ELSE
         RotorType = TRIM(RotorType)//' with rigid hub.'
      ENDIF
   ELSE
      RotorType = TRIM(RotorType)//' three-bladed rotor with rigid hub.'
   ENDIF

   WRITE    (UnSu,'(A)')  '            '//TRIM(RotorType)

   WRITE    (UnSu,'(A)')  '            The model has '//TRIM(Num2LStr( p%DOFs%NActvDOF ))//' of '// &
                                        TRIM(Num2LStr( p%NDOF ))//' DOFs active (enabled) at start-up.'

   DO I=1,p%NDOF
   
      IF ( p%DOF_Flag( I ) ) THEN
         DOFEnabled = 'Enabled'
      ELSE
         DOFEnabled = 'Disabled'
      END IF
      
      K = INDEX( p%DOF_Desc(I), ' (internal DOF index' )
      IF (K == 0) K = LEN_TRIM(p%DOF_Desc(I))
      
      WRITE ( UnSu, '(1x,A10,1x,A)' ) DOFEnabled, p%DOF_Desc(I)(:K)
         
      
   END DO
   

      ! Time steps.

   WRITE (UnSu,'(//,A,/)')  'Time steps:'

   WRITE (UnSu,FmtDatT) '    Structural            (s)     ', p%DT

      ! Some calculated parameters.

   WRITE (UnSu,'(//,A,/)')  'Some calculated parameters:'

   WRITE (UnSu,FmtDat ) '    Hub-Height            (m)     ', p%HubHt
   WRITE (UnSu,FmtDat ) '    Flexible Tower Length (m)     ', p%TwrFlexL
IF (.NOT. p%BD4Blades) THEN
   WRITE (UnSu,FmtDat ) '    Flexible Blade Length (m)     ', p%BldFlexL
ELSE
   WRITE (UnSu,'(A)' )  '    Flexible Blade Length (m)               N/A'
END IF

      ! Rotor properties:

   WRITE (UnSu,'(//,A,/)')  'Rotor mass properties:'

   WRITE (UnSu,FmtDat ) '    Rotor Mass            (kg)    ', p%RotMass
   WRITE (UnSu,FmTDat ) '    Rotor Inertia         (kg-m^2)', p%RotINer

   WRITE (UnSu,Fmt1   ) ( K,         K=1,p%NumBl )
   WRITE (UnSu,Fmt2   ) ( '-------', K=1,p%NumBl )

IF (.NOT. p%BD4Blades) THEN

   WRITE (UnSu,FmtDat ) '    Mass                  (kg)    ', ( p%BldMass  (K), K=1,p%NumBl )
   WRITE (UnSu,FmtDat ) '    Second Mass Moment    (kg-m^2)', ( p%SecondMom(K), K=1,p%NumBl )
   WRITE (UnSu,FmtDat ) '    First Mass Moment     (kg-m)  ', ( p%FirstMom (K), K=1,p%NumBl )
   WRITE (UnSu,FmtDat ) '    Center of Mass        (m)     ', ( p%BldCG    (K), K=1,p%NumBl )

ELSE

   WRITE (UnSu,'(A)' ) '    Mass                  (kg)              N/A'
   WRITE (UnSu,'(A)' ) '    Second Mass Moment    (kg-m^2)          N/A'
   WRITE (UnSu,'(A)' ) '    First Mass Moment     (kg-m)            N/A'
   WRITE (UnSu,'(A)' ) '    Center of Mass        (m)               N/A'
   
END IF 

      ! Output additional masses:

   WRITE (UnSu,'(//,A,/)')  'Additional mass properties:'

   WRITE (UnSu,FmtDat ) '    Tower-top Mass        (kg)    ', p%TwrTpMass
   WRITE (UnSu,FmtDat ) '    Tower Mass            (kg)    ', p%TwrMass
   !WRITE (UnSu,FmtDat ) '    Turbine Mass          (kg)    ', p%TurbMass
   WRITE (UnSu,FmtDat ) '    Platform Mass         (kg)    ', p%PtfmMass
   WRITE (UnSu,FmtDat ) '    Mass Incl. Platform   (kg)    ', p%TurbMass + p%PtfmMass !TotalMass !bjj TotalMass not used anywhere else so removed it


      ! Interpolated tower properties.

   WRITE (UnSu,"(//,'Interpolated tower properties:',/)")
   IF ( GenerateAdamsModel )  THEN  ! An ADAMS model will be created; thus, print out all the cols.

      WRITE (UnSu,'(A)')  'Node  TwFract   HNodes  DHNodes  TMassDen    FAStiff    SSStiff'// &
                           '    GJStiff    EAStiff    FAIner    SSIner  FAcgOff  SScgOff'
      WRITE (UnSu,'(A)')  ' (-)      (-)      (m)      (m)    (kg/m)     (Nm^2)     (Nm^2)'// &
                           '     (Nm^2)        (N)    (kg m)    (kg m)      (m)      (m)'

      DO I=1,p%TwrNodes
         WRITE(UnSu,'(I4,3F9.3,F10.3,4ES11.3,2F10.3,2F9.3)')  I, p%HNodesNorm(I), p%HNodes(I), p%DHNodes(I), p%MassT(I), &
                                                                 p%StiffTFA(I), p%StiffTSS(I), p%StiffTGJ(I), p%StiffTEA(I),       &
                                                                 p%InerTFA(I), p%InerTSS(I), p%cgOffTFA(I), p%cgOffTSS(I)
      ENDDO ! I

   ELSE                                                     ! Only FAST will be run; thus, only print out the necessary cols.

      WRITE (UnSu,'(A)')  'Node  TwFract   HNodes  DHNodes  TMassDen    FAStiff    SSStiff'
      WRITE (UnSu,'(A)')  ' (-)      (-)      (m)      (m)    (kg/m)     (Nm^2)     (Nm^2)'

      DO I=1,p%TwrNodes
         WRITE(UnSu,'(I4,3F9.3,F10.3,2ES11.3)')  I, p%HNodesNorm(I), p%HNodes(I), p%DHNodes(I), p%MassT(I), &
                                                    p%StiffTFA(I),  p%StiffTSS(I)
      ENDDO ! I

   ENDIF


      ! Interpolated blade properties.

   DO K=1,p%NumBl

      WRITE (UnSu,'(//,A,I1,A,/)')  'Interpolated blade ', K, ' properties:'
      IF ( GenerateAdamsModel )  THEN  ! An ADAMS model will be created; thus, print out all the cols.

         WRITE (UnSu,'(A)')  'Node  BlFract   RNodes  DRNodes  PitchAxis  StrcTwst  BMassDen    FlpStff    EdgStff'//       &
                              '     GJStff     EAStff    Alpha   FlpIner   EdgIner PrecrvRef PreswpRef  FlpcgOf  EdgcgOf'// &
                              '  FlpEAOf  EdgEAOf'
         WRITE (UnSu,'(A)')  ' (-)      (-)      (m)      (m)     (-)       (deg)    (kg/m)     (Nm^2)     (Nm^2)'//       &
                              '     (Nm^2)     (Nm^2)      (-)    (kg m)    (kg m)       (m)       (m)      (m)      (m)'// &
                              '      (m)      (m)'

         DO I=1,p%BldNodes

            WRITE(UnSu,'(I4,3F9.3,3F10.3,4ES11.3,F9.3,4F10.3,4F9.3)')  I, p%RNodesNorm(I),  p%RNodes(I) + p%HubRad, p%DRNodes(I), &
                                                                          p%PitchAxis(K,I), p%ThetaS(K,I)*R2D,      p%MassB(K,I), &
                                                                          p%StiffBF(K,I),   p%StiffBE(K,I),                       &
                                                                          p%StiffBGJ(K,I),  p%StiffBEA(K,I),                      &
                                                                          p%BAlpha(K,I),    p%InerBFlp(K,I), p%InerBEdg(K,I),     &
                                                                          p%RefAxisxb(K,I), p%RefAxisyb(K,I),                     &
                                                                          p%cgOffBFlp(K,I), p%cgOffBEdg(K,I),                     &
                                                                          p%EAOffBFlp(K,I), p%EAOffBEdg(K,I)
         ENDDO ! I

      ELSE                                                     ! Only FAST will be run; thus, only print out the necessary cols.

         WRITE (UnSu,'(A)')  'Node  BlFract   RNodes  DRNodes PitchAxis  StrcTwst  BMassDen    FlpStff    EdgStff'
         WRITE (UnSu,'(A)')  ' (-)      (-)      (m)      (m)       (-)     (deg)    (kg/m)     (Nm^2)     (Nm^2)'

         DO I=1,p%BldNodes
            WRITE(UnSu,'(I4,3F9.3,3F10.3,2ES11.3)')  I, p%RNodesNorm(I), p%RNodes(I) + p%HubRad, p%DRNodes(I), &
                                                        p%PitchAxis(K,I),p%ThetaS(K,I)*R2D, p%MassB(K,I), &
                                                        p%StiffBF(K,I), p%StiffBE(K,I)
         ENDDO ! I

      ENDIF

   ENDDO ! K


   OutPFmt = '( I4, 3X,A '//TRIM(Num2LStr(ChanLen))//',1 X, A'//TRIM(Num2LStr(ChanLen))//' )'
   WRITE (UnSu,'(//,A,/)')  'Requested Outputs:'
   WRITE (UnSu,"(/, '  Col  Parameter  Units', /, '  ---  ---------  -----')")
   DO I = 0,p%NumOuts
      WRITE (UnSu,OutPFmt)  I, p%OutParam(I)%Name, p%OutParam(I)%Units
   END DO             

   CLOSE(UnSu)

RETURN
END SUBROUTINE ED_PrintSum
!=======================================================================

!=======================================================================
SUBROUTINE FixHSSBrTq ( Integrator, p, x, OtherState, ErrStat, ErrMsg )


   ! This routine is used to adjust the HSSBrTrq value if the absolute
   !   magnitudue of the HSS brake torque was strong enough to reverse
   !   the direction of the HSS, which is a physically impossible
   !   situation.  The problem arises since we are integrating in
   !   discrete time, not continuous time.

   IMPLICIT                        NONE


   ! Passed variables:

   TYPE(ED_ParameterType),      INTENT(IN   ):: p                                 ! Parameters of the structural dynamics module
   TYPE(ED_OtherStateType),     INTENT(INOUT):: OtherState                        ! Other/optimization states of the structural dynamics module 
   TYPE(ED_ContinuousStateType),INTENT(INOUT):: x                                 ! Continuous states of the structural dynamics module at n+1
   CHARACTER(1),                INTENT(IN   ):: Integrator                        ! A string holding the current integrator being used.
   INTEGER(IntKi),              INTENT(  OUT):: ErrStat
   CHARACTER(*),                INTENT(  OUT):: ErrMsg


   ! Local variables:

   REAL(ReKi)                             :: RqdFrcGeAz                           ! The force term required to produce RqdQD2GeAz.
   REAL(ReKi)                             :: RqdQD2GeAz                           ! The required QD2T(DOF_GeAz) to cause the HSS to stop rotating.

   INTEGER                                :: I                                    ! Loops through all DOFs.
   INTEGER(IntKi)                         :: ErrStat2
   CHARACTER(ErrMsgLen)                   :: ErrMsg2
   CHARACTER(*), PARAMETER                :: RoutineName = 'FixHSSBrTq'

   

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF ( .NOT. p%DOF_Flag(DOF_GeAz) .OR. EqualRealNos(OtherState%HSSBrTrqC, 0.0_ReKi ) )  RETURN


      ! The absolute magnitude of the HSS brake must have been too great
      !   that the HSS direction was reversed.  What should have happened
      !   is that the HSS should have stopped rotating.  In other words,
      !   QD(DOF_GeAz,IC(NMX)) should equal zero!  Determining what
      !   QD2T(DOF_GeAz) will make QD(DOF_GeAz,IC(NMX)) = 0, depends on
      !   which integrator we are using.

   
   SELECT CASE (Integrator)

   CASE ('C')   ! Corrector

      ! Find the required QD2T(DOF_GeAz) to cause the HSS to stop rotating (RqdQD2GeAz).
      ! This is found by solving the corrector formula for QD2(DOF_GeAz,IC(NMX))
      !   when QD(DOF_GeAz,IC(NMX)) equals zero.

      RqdQD2GeAz = ( -      OtherState%xdot(OtherState%IC(1))%qt (DOF_GeAz)/ p%DT24 &
                     - 19.0*OtherState%xdot(OtherState%IC(1))%qdt(DOF_GeAz)         &
                     +  5.0*OtherState%xdot(OtherState%IC(2))%qdt(DOF_GeAz)         &
                     -      OtherState%xdot(OtherState%IC(3))%qdt(DOF_GeAz)         ) / 9.0
      
   CASE ('P')   ! Predictor

      ! Find the required QD2T(DOF_GeAz) to cause the HSS to stop rotating (RqdQD2GeAz).
      ! This is found by solving the predictor formula for QD2(DOF_GeAz,IC(1))
      !   when QD(DOF_GeAz,IC(NMX)) equals zero.

      RqdQD2GeAz = ( -      OtherState%xdot(OtherState%IC(1))%qt( DOF_GeAz)  / p%DT24 &
                     + 59.0*OtherState%xdot(OtherState%IC(2))%qdt(DOF_GeAz) &
                     - 37.0*OtherState%xdot(OtherState%IC(3))%qdt(DOF_GeAz) &
                     +  9.0*OtherState%xdot(OtherState%IC(4))%qdt(DOF_GeAz)   )/55.0
            
   END SELECT


   ! Rearrange the augmented matrix of equations of motion to account
   !   for the known acceleration of the generator azimuth DOF.  To
   !   do this, make the known inertia like an applied force to the
   !   system.  Then set force QD2T(DOF_GeAz) to equal the known
   !   acceleration in the augmented matrix of equations of motion:
   ! Here is how the new equations are derived.  First partition the
   !   augmented matrix as follows, where Qa are the unknown
   !   accelerations, Qb are the known accelerations, Fa are the
   !   known forces, and Fb are the unknown forces:
   !      [Caa Cab]{Qa}={Fa}
   !      [Cba Cbb]{Qb}={Fb}
   !   By rearranging, the equations for the unknown and known
   !   accelerations are as follows:
   !      [Caa]{Qa}={Fa}-[Cab]{Qb} and [I]{Qb}={Qb}
   !   Combining these two sets of equations into one set yields:
   !      [Caa 0]{Qa}={{Fa}-[Cab]{Qb}}
   !      [  0 I]{Qb}={          {Qb}}
   !   Once this equation is solved, the unknown force can be found from:
   !      {Fb}=[Cba]{Qa}+[Cbb]{Qb}

   OtherState%OgnlGeAzRo    = OtherState%AugMat(DOF_GeAz,:)  ! used for HSS Brake hack; copy this row before modifying the old matrix
   
  
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs

      OtherState%AugMat(p%DOFs%SrtPS(I),    p%NAUG) = OtherState%AugMat(p%DOFs%SrtPS(I),p%NAUG) &
                                                    - OtherState%AugMat(p%DOFs%SrtPS(I),DOF_GeAz)*RqdQD2GeAz  ! {{Fa}-[Cab]{Qb}}
      OtherState%AugMat(p%DOFs%SrtPS(I),DOF_GeAz)   = 0.0                                                     ! [0]
      OtherState%AugMat(DOF_GeAz, p%DOFs%SrtPS(I))  = 0.0                                                     ! [0]

   ENDDO             ! I - All active (enabled) DOFs

   OtherState%AugMat(DOF_GeAz,DOF_GeAz) = 1.0                                                           ! [I]{Qb}={Qb}
   OtherState%AugMat(DOF_GeAz,  p%NAUG) = RqdQD2GeAz                                                    !


   ! Invert the matrix to solve for the new (updated) accelerations.  Like in
   !   CalcContStateDeriv(), the accelerations are returned by Gauss() in the first NActvDOF
   !   elements of the solution vector, SolnVec().  These are transfered to the
   !   proper index locations of the acceleration vector QD2T() using the
   !   vector subscript array SrtPS(), after Gauss() has been called:

   ! Invert the matrix to solve for the accelerations. The accelerations are returned by Gauss() in the first NActvDOF elements
   !   of the solution vector, SolnVec(). These are transfered to the proper index locations of the acceleration vector QD2T()
   !   using the vector subscript array SrtPS(), after Gauss() has been called:

      OtherState%AugMat_factor = OtherState%AugMat( p%DOFs%SrtPS( 1:p%DOFs%NActvDOF ), p%DOFs%SrtPSNAUG(1:p%DOFs%NActvDOF) )
      OtherState%SolnVec       = OtherState%AugMat( p%DOFs%SrtPS( 1:p%DOFs%NActvDOF ), p%DOFs%SrtPSNAUG(1+p%DOFs%NActvDOF) )
   
      CALL LAPACK_getrf( M=p%DOFs%NActvDOF, N=p%DOFs%NActvDOF, A=OtherState%AugMat_factor, IPIV=OtherState%AugMat_pivot, ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)   
         IF ( ErrStat >= AbortErrLev ) RETURN
      
      CALL LAPACK_getrs( TRANS='N',N=p%DOFs%NActvDOF, A=OtherState%AugMat_factor,IPIV=OtherState%AugMat_pivot, B=OtherState%SolnVec, ErrStat=ErrStat2, ErrMsg=ErrMsg2)
   
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) RETURN
   
   
      ! Find the force required to produce RqdQD2GeAz from the equations of
      !   motion using the new accelerations:

   RqdFrcGeAz = 0.0
   DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
      ! bjj: use OtherState%SolnVec(I) instead of OtherState%QD2T(p%DOFs%SrtPS(I)) here; then update OtherState%QD2T(p%DOFs%SrtPS(I))
      !      later if necessary
      !RqdFrcGeAz = RqdFrcGeAz + OgnlGeAzRo(SrtPS(I))*OtherState%QD2T(p%DOFs%SrtPS(I))  ! {Fb}=[Cba]{Qa}+[Cbb]{Qb}
      RqdFrcGeAz = RqdFrcGeAz + OtherState%OgnlGeAzRo(p%DOFs%SrtPS(I))*OtherState%SolnVec(I)  ! {Fb}=[Cba]{Qa}+[Cbb]{Qb}
   ENDDO             ! I - All active (enabled) DOFs


      ! Find the HSSBrTrq necessary to bring about this force:

   OtherState%HSSBrTrq = OtherState%HSSBrTrqC & 
                       + ( ( OtherState%OgnlGeAzRo(p%NAUG) - RqdFrcGeAz )*OtherState%RtHS%GBoxEffFac/ABS(p%GBRatio) )


      ! Make sure this new HSSBrTrq isn't larger in absolute magnitude than
      !   the original HSSBrTrq.  Indeed, the new HSSBrTrq can't be larger than
      !   the old HSSBrTrq, since the old HSSBrTrq was found solely as a
      !   function of time--and is thus the maximum possible at the current
      !   time.  If the new HSSBrTrq is larger, then the reversal in direction
      !   was caused by factors other than the HSS brake--thus the original HSS
      !   brake torque values were OK to begin with.  Thus, restore the
      !   variables changed by this subroutine, back to their original values:

   IF ( ABS( OtherState%HSSBrTrq ) > ABS( OtherState%HSSBrTrqC ) )  THEN

      OtherState%HSSBrTrq = OtherState%HSSBrTrqC !OtherState%HSSBrTrqC = SIGN( u%HSSBrTrqC, x%QDT(DOF_GeAz) )
      !OtherState%QD2T     = QD2TC

   ELSE

      ! overwrite QD2T with the new values
      OtherState%QD2T = 0.0
      DO I = 1,p%DOFs%NActvDOF ! Loop through all active (enabled) DOFs
         OtherState%QD2T(p%DOFs%SrtPS(I)) = OtherState%SolnVec(I)
      ENDDO             ! I - All active (enabled) DOFs
      
            
      ! Use the new accelerations to update the DOF values.  Again, this
      !   depends on the integrator type:

      SELECT CASE (Integrator)

      CASE ('C')  ! Corrector

      ! Update QD and QD2 with the new accelerations using the corrector.
      ! This will make QD(DOF_GeAz,IC(NMX)) equal to zero and adjust all
      !    of the other QDs as necessary.
      ! The Q's are unnaffected by this change.     
      
         x%qdt =                   OtherState%xdot(OtherState%IC(1))%qt &  ! qd at n
                 + p%DT24 * ( 9. * OtherState%QD2T &                 ! the value we just changed
                           + 19. * OtherState%xdot(OtherState%IC(1))%qdt &
                           -  5. * OtherState%xdot(OtherState%IC(2))%qdt &
                           +  1. * OtherState%xdot(OtherState%IC(3))%qdt )
            

         
      CASE ('P')  ! Predictor

      ! Update QD and QD2 with the new accelerations using predictor.  
         
         x%qdt =                OtherState%xdot(OtherState%IC(1))%qt + &  ! qd at n
                 p%DT24 * ( 55.*OtherState%QD2T &                    ! the value we just changed
                          - 59.*OtherState%xdot(OtherState%IC(2))%qdt  &
                          + 37.*OtherState%xdot(OtherState%IC(3))%qdt  &
                           - 9.*OtherState%xdot(OtherState%IC(4))%qdt )
         
         OtherState%xdot ( OtherState%IC(1) )%qdt = OtherState%QD2T        ! fix the history

         
      END SELECT
      
   ENDIF

   RETURN
END SUBROUTINE FixHSSBrTq
!=======================================================================


END MODULE ElastoDyn
!**********************************************************************************************************************************
