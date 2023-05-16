!**********************************************************************************************************************************
! The ElastoDyn.f90, ElastoDyn_IO.f90 and  ElastoDyn_Types.f90 make up the ElastoDyn module of the
! FAST Modularization Framework. ElastoDyn_Types is auto-generated based on FAST_Registry.txt.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2016 National Renewable Energy Laboratory
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
!> This module contains definitions of compile-time PARAMETERS for the ElastoDyn module.
!! Every variable defined here MUST have the PARAMETER attribute.
MODULE ElastoDyn_Parameters

   USE NWTC_Library

   USE ElastoDyn_AllBldNdOuts_IO

   TYPE(ProgDesc), PARAMETER  :: ED_Ver = ProgDesc( 'ElastoDyn', '', '' )
   
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
   INTEGER(IntKi), PARAMETER        :: ED_MaxDOFs  = 24


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
! This code was generated by Write_ChckOutLst.m at 25-Jan-2021 13:23:51.


     ! Indices for computing output channels:
     ! NOTES:
     !    (1) These parameters are in the order stored in "OutListParameters.xlsx"
     !    (2) Array AllOuts() must be dimensioned to the value of the largest output parameter

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

   INTEGER(IntKi), PARAMETER      :: TwrTpTDxi = 315
   INTEGER(IntKi), PARAMETER      :: TwrTpTDyi = 316
   INTEGER(IntKi), PARAMETER      :: TwrTpTDzi = 317
   INTEGER(IntKi), PARAMETER      :: YawBrTDxp = 318
   INTEGER(IntKi), PARAMETER      :: YawBrTDyp = 319
   INTEGER(IntKi), PARAMETER      :: YawBrTDzp = 320
   INTEGER(IntKi), PARAMETER      :: YawBrTDxt = 321
   INTEGER(IntKi), PARAMETER      :: YawBrTDyt = 322
   INTEGER(IntKi), PARAMETER      :: YawBrTDzt = 323
   INTEGER(IntKi), PARAMETER      :: YawBrTVxp = 324
   INTEGER(IntKi), PARAMETER      :: YawBrTVyp = 325
   INTEGER(IntKi), PARAMETER      :: YawBrTVzp = 326
   INTEGER(IntKi), PARAMETER      :: YawBrTAxp = 327
   INTEGER(IntKi), PARAMETER      :: YawBrTAyp = 328
   INTEGER(IntKi), PARAMETER      :: YawBrTAzp = 329
   INTEGER(IntKi), PARAMETER      :: YawBrRDxt = 330
   INTEGER(IntKi), PARAMETER      :: YawBrRDyt = 331
   INTEGER(IntKi), PARAMETER      :: YawBrRDzt = 332
   INTEGER(IntKi), PARAMETER      :: YawBrRVxp = 333
   INTEGER(IntKi), PARAMETER      :: YawBrRVyp = 334
   INTEGER(IntKi), PARAMETER      :: YawBrRVzp = 335
   INTEGER(IntKi), PARAMETER      :: YawBrRAxp = 336
   INTEGER(IntKi), PARAMETER      :: YawBrRAyp = 337
   INTEGER(IntKi), PARAMETER      :: YawBrRAzp = 338


     ! Local Tower Motions:

   INTEGER(IntKi), PARAMETER      :: TwHt1ALxt = 339
   INTEGER(IntKi), PARAMETER      :: TwHt1ALyt = 340
   INTEGER(IntKi), PARAMETER      :: TwHt1ALzt = 341
   INTEGER(IntKi), PARAMETER      :: TwHt2ALxt = 342
   INTEGER(IntKi), PARAMETER      :: TwHt2ALyt = 343
   INTEGER(IntKi), PARAMETER      :: TwHt2ALzt = 344
   INTEGER(IntKi), PARAMETER      :: TwHt3ALxt = 345
   INTEGER(IntKi), PARAMETER      :: TwHt3ALyt = 346
   INTEGER(IntKi), PARAMETER      :: TwHt3ALzt = 347
   INTEGER(IntKi), PARAMETER      :: TwHt4ALxt = 348
   INTEGER(IntKi), PARAMETER      :: TwHt4ALyt = 349
   INTEGER(IntKi), PARAMETER      :: TwHt4ALzt = 350
   INTEGER(IntKi), PARAMETER      :: TwHt5ALxt = 351
   INTEGER(IntKi), PARAMETER      :: TwHt5ALyt = 352
   INTEGER(IntKi), PARAMETER      :: TwHt5ALzt = 353
   INTEGER(IntKi), PARAMETER      :: TwHt6ALxt = 354
   INTEGER(IntKi), PARAMETER      :: TwHt6ALyt = 355
   INTEGER(IntKi), PARAMETER      :: TwHt6ALzt = 356
   INTEGER(IntKi), PARAMETER      :: TwHt7ALxt = 357
   INTEGER(IntKi), PARAMETER      :: TwHt7ALyt = 358
   INTEGER(IntKi), PARAMETER      :: TwHt7ALzt = 359
   INTEGER(IntKi), PARAMETER      :: TwHt8ALxt = 360
   INTEGER(IntKi), PARAMETER      :: TwHt8ALyt = 361
   INTEGER(IntKi), PARAMETER      :: TwHt8ALzt = 362
   INTEGER(IntKi), PARAMETER      :: TwHt9ALxt = 363
   INTEGER(IntKi), PARAMETER      :: TwHt9ALyt = 364
   INTEGER(IntKi), PARAMETER      :: TwHt9ALzt = 365
   INTEGER(IntKi), PARAMETER      :: TwHt1TDxt = 366
   INTEGER(IntKi), PARAMETER      :: TwHt1TDyt = 367
   INTEGER(IntKi), PARAMETER      :: TwHt1TDzt = 368
   INTEGER(IntKi), PARAMETER      :: TwHt2TDxt = 369
   INTEGER(IntKi), PARAMETER      :: TwHt2TDyt = 370
   INTEGER(IntKi), PARAMETER      :: TwHt2TDzt = 371
   INTEGER(IntKi), PARAMETER      :: TwHt3TDxt = 372
   INTEGER(IntKi), PARAMETER      :: TwHt3TDyt = 373
   INTEGER(IntKi), PARAMETER      :: TwHt3TDzt = 374
   INTEGER(IntKi), PARAMETER      :: TwHt4TDxt = 375
   INTEGER(IntKi), PARAMETER      :: TwHt4TDyt = 376
   INTEGER(IntKi), PARAMETER      :: TwHt4TDzt = 377
   INTEGER(IntKi), PARAMETER      :: TwHt5TDxt = 378
   INTEGER(IntKi), PARAMETER      :: TwHt5TDyt = 379
   INTEGER(IntKi), PARAMETER      :: TwHt5TDzt = 380
   INTEGER(IntKi), PARAMETER      :: TwHt6TDxt = 381
   INTEGER(IntKi), PARAMETER      :: TwHt6TDyt = 382
   INTEGER(IntKi), PARAMETER      :: TwHt6TDzt = 383
   INTEGER(IntKi), PARAMETER      :: TwHt7TDxt = 384
   INTEGER(IntKi), PARAMETER      :: TwHt7TDyt = 385
   INTEGER(IntKi), PARAMETER      :: TwHt7TDzt = 386
   INTEGER(IntKi), PARAMETER      :: TwHt8TDxt = 387
   INTEGER(IntKi), PARAMETER      :: TwHt8TDyt = 388
   INTEGER(IntKi), PARAMETER      :: TwHt8TDzt = 389
   INTEGER(IntKi), PARAMETER      :: TwHt9TDxt = 390
   INTEGER(IntKi), PARAMETER      :: TwHt9TDyt = 391
   INTEGER(IntKi), PARAMETER      :: TwHt9TDzt = 392
   INTEGER(IntKi), PARAMETER      :: TwHt1RDxt = 393
   INTEGER(IntKi), PARAMETER      :: TwHt1RDyt = 394
   INTEGER(IntKi), PARAMETER      :: TwHt1RDzt = 395
   INTEGER(IntKi), PARAMETER      :: TwHt2RDxt = 396
   INTEGER(IntKi), PARAMETER      :: TwHt2RDyt = 397
   INTEGER(IntKi), PARAMETER      :: TwHt2RDzt = 398
   INTEGER(IntKi), PARAMETER      :: TwHt3RDxt = 399
   INTEGER(IntKi), PARAMETER      :: TwHt3RDyt = 400
   INTEGER(IntKi), PARAMETER      :: TwHt3RDzt = 401
   INTEGER(IntKi), PARAMETER      :: TwHt4RDxt = 402
   INTEGER(IntKi), PARAMETER      :: TwHt4RDyt = 403
   INTEGER(IntKi), PARAMETER      :: TwHt4RDzt = 404
   INTEGER(IntKi), PARAMETER      :: TwHt5RDxt = 405
   INTEGER(IntKi), PARAMETER      :: TwHt5RDyt = 406
   INTEGER(IntKi), PARAMETER      :: TwHt5RDzt = 407
   INTEGER(IntKi), PARAMETER      :: TwHt6RDxt = 408
   INTEGER(IntKi), PARAMETER      :: TwHt6RDyt = 409
   INTEGER(IntKi), PARAMETER      :: TwHt6RDzt = 410
   INTEGER(IntKi), PARAMETER      :: TwHt7RDxt = 411
   INTEGER(IntKi), PARAMETER      :: TwHt7RDyt = 412
   INTEGER(IntKi), PARAMETER      :: TwHt7RDzt = 413
   INTEGER(IntKi), PARAMETER      :: TwHt8RDxt = 414
   INTEGER(IntKi), PARAMETER      :: TwHt8RDyt = 415
   INTEGER(IntKi), PARAMETER      :: TwHt8RDzt = 416
   INTEGER(IntKi), PARAMETER      :: TwHt9RDxt = 417
   INTEGER(IntKi), PARAMETER      :: TwHt9RDyt = 418
   INTEGER(IntKi), PARAMETER      :: TwHt9RDzt = 419
   INTEGER(IntKi), PARAMETER      :: TwHt1TPxi = 420
   INTEGER(IntKi), PARAMETER      :: TwHt1TPyi = 421
   INTEGER(IntKi), PARAMETER      :: TwHt1TPzi = 422
   INTEGER(IntKi), PARAMETER      :: TwHt2TPxi = 423
   INTEGER(IntKi), PARAMETER      :: TwHt2TPyi = 424
   INTEGER(IntKi), PARAMETER      :: TwHt2TPzi = 425
   INTEGER(IntKi), PARAMETER      :: TwHt3TPxi = 426
   INTEGER(IntKi), PARAMETER      :: TwHt3TPyi = 427
   INTEGER(IntKi), PARAMETER      :: TwHt3TPzi = 428
   INTEGER(IntKi), PARAMETER      :: TwHt4TPxi = 429
   INTEGER(IntKi), PARAMETER      :: TwHt4TPyi = 430
   INTEGER(IntKi), PARAMETER      :: TwHt4TPzi = 431
   INTEGER(IntKi), PARAMETER      :: TwHt5TPxi = 432
   INTEGER(IntKi), PARAMETER      :: TwHt5TPyi = 433
   INTEGER(IntKi), PARAMETER      :: TwHt5TPzi = 434
   INTEGER(IntKi), PARAMETER      :: TwHt6TPxi = 435
   INTEGER(IntKi), PARAMETER      :: TwHt6TPyi = 436
   INTEGER(IntKi), PARAMETER      :: TwHt6TPzi = 437
   INTEGER(IntKi), PARAMETER      :: TwHt7TPxi = 438
   INTEGER(IntKi), PARAMETER      :: TwHt7TPyi = 439
   INTEGER(IntKi), PARAMETER      :: TwHt7TPzi = 440
   INTEGER(IntKi), PARAMETER      :: TwHt8TPxi = 441
   INTEGER(IntKi), PARAMETER      :: TwHt8TPyi = 442
   INTEGER(IntKi), PARAMETER      :: TwHt8TPzi = 443
   INTEGER(IntKi), PARAMETER      :: TwHt9TPxi = 444
   INTEGER(IntKi), PARAMETER      :: TwHt9TPyi = 445
   INTEGER(IntKi), PARAMETER      :: TwHt9TPzi = 446
   INTEGER(IntKi), PARAMETER      :: TwHt1RPxi = 447
   INTEGER(IntKi), PARAMETER      :: TwHt1RPyi = 448
   INTEGER(IntKi), PARAMETER      :: TwHt1RPzi = 449
   INTEGER(IntKi), PARAMETER      :: TwHt2RPxi = 450
   INTEGER(IntKi), PARAMETER      :: TwHt2RPyi = 451
   INTEGER(IntKi), PARAMETER      :: TwHt2RPzi = 452
   INTEGER(IntKi), PARAMETER      :: TwHt3RPxi = 453
   INTEGER(IntKi), PARAMETER      :: TwHt3RPyi = 454
   INTEGER(IntKi), PARAMETER      :: TwHt3RPzi = 455
   INTEGER(IntKi), PARAMETER      :: TwHt4RPxi = 456
   INTEGER(IntKi), PARAMETER      :: TwHt4RPyi = 457
   INTEGER(IntKi), PARAMETER      :: TwHt4RPzi = 458
   INTEGER(IntKi), PARAMETER      :: TwHt5RPxi = 459
   INTEGER(IntKi), PARAMETER      :: TwHt5RPyi = 460
   INTEGER(IntKi), PARAMETER      :: TwHt5RPzi = 461
   INTEGER(IntKi), PARAMETER      :: TwHt6RPxi = 462
   INTEGER(IntKi), PARAMETER      :: TwHt6RPyi = 463
   INTEGER(IntKi), PARAMETER      :: TwHt6RPzi = 464
   INTEGER(IntKi), PARAMETER      :: TwHt7RPxi = 465
   INTEGER(IntKi), PARAMETER      :: TwHt7RPyi = 466
   INTEGER(IntKi), PARAMETER      :: TwHt7RPzi = 467
   INTEGER(IntKi), PARAMETER      :: TwHt8RPxi = 468
   INTEGER(IntKi), PARAMETER      :: TwHt8RPyi = 469
   INTEGER(IntKi), PARAMETER      :: TwHt8RPzi = 470
   INTEGER(IntKi), PARAMETER      :: TwHt9RPxi = 471
   INTEGER(IntKi), PARAMETER      :: TwHt9RPyi = 472
   INTEGER(IntKi), PARAMETER      :: TwHt9RPzi = 473


     ! Platform Motions:

   INTEGER(IntKi), PARAMETER      :: PtfmTDxt  = 474
   INTEGER(IntKi), PARAMETER      :: PtfmTDyt  = 475
   INTEGER(IntKi), PARAMETER      :: PtfmTDzt  = 476
   INTEGER(IntKi), PARAMETER      :: PtfmTDxi  = 477
   INTEGER(IntKi), PARAMETER      :: PtfmTDyi  = 478
   INTEGER(IntKi), PARAMETER      :: PtfmTDzi  = 479
   INTEGER(IntKi), PARAMETER      :: PtfmTVxt  = 480
   INTEGER(IntKi), PARAMETER      :: PtfmTVyt  = 481
   INTEGER(IntKi), PARAMETER      :: PtfmTVzt  = 482
   INTEGER(IntKi), PARAMETER      :: PtfmTVxi  = 483
   INTEGER(IntKi), PARAMETER      :: PtfmTVyi  = 484
   INTEGER(IntKi), PARAMETER      :: PtfmTVzi  = 485
   INTEGER(IntKi), PARAMETER      :: PtfmTAxt  = 486
   INTEGER(IntKi), PARAMETER      :: PtfmTAyt  = 487
   INTEGER(IntKi), PARAMETER      :: PtfmTAzt  = 488
   INTEGER(IntKi), PARAMETER      :: PtfmTAxi  = 489
   INTEGER(IntKi), PARAMETER      :: PtfmTAyi  = 490
   INTEGER(IntKi), PARAMETER      :: PtfmTAzi  = 491
   INTEGER(IntKi), PARAMETER      :: PtfmRDxi  = 492
   INTEGER(IntKi), PARAMETER      :: PtfmRDyi  = 493
   INTEGER(IntKi), PARAMETER      :: PtfmRDzi  = 494
   INTEGER(IntKi), PARAMETER      :: PtfmRVxt  = 495
   INTEGER(IntKi), PARAMETER      :: PtfmRVyt  = 496
   INTEGER(IntKi), PARAMETER      :: PtfmRVzt  = 497
   INTEGER(IntKi), PARAMETER      :: PtfmRVxi  = 498
   INTEGER(IntKi), PARAMETER      :: PtfmRVyi  = 499
   INTEGER(IntKi), PARAMETER      :: PtfmRVzi  = 500
   INTEGER(IntKi), PARAMETER      :: PtfmRAxt  = 501
   INTEGER(IntKi), PARAMETER      :: PtfmRAyt  = 502
   INTEGER(IntKi), PARAMETER      :: PtfmRAzt  = 503
   INTEGER(IntKi), PARAMETER      :: PtfmRAxi  = 504
   INTEGER(IntKi), PARAMETER      :: PtfmRAyi  = 505
   INTEGER(IntKi), PARAMETER      :: PtfmRAzi  = 506


     ! Blade 1 Root Loads:

   INTEGER(IntKi), PARAMETER      :: RootFxc1  = 507
   INTEGER(IntKi), PARAMETER      :: RootFyc1  = 508
   INTEGER(IntKi), PARAMETER      :: RootFzc1  = 509
   INTEGER(IntKi), PARAMETER      :: RootFxb1  = 510
   INTEGER(IntKi), PARAMETER      :: RootFyb1  = 511
   INTEGER(IntKi), PARAMETER      :: RootMxc1  = 512
   INTEGER(IntKi), PARAMETER      :: RootMyc1  = 513
   INTEGER(IntKi), PARAMETER      :: RootMzc1  = 514
   INTEGER(IntKi), PARAMETER      :: RootMxb1  = 515
   INTEGER(IntKi), PARAMETER      :: RootMyb1  = 516


     ! Blade 2 Root Loads:

   INTEGER(IntKi), PARAMETER      :: RootFxc2  = 517
   INTEGER(IntKi), PARAMETER      :: RootFyc2  = 518
   INTEGER(IntKi), PARAMETER      :: RootFzc2  = 519
   INTEGER(IntKi), PARAMETER      :: RootFxb2  = 520
   INTEGER(IntKi), PARAMETER      :: RootFyb2  = 521
   INTEGER(IntKi), PARAMETER      :: RootMxc2  = 522
   INTEGER(IntKi), PARAMETER      :: RootMyc2  = 523
   INTEGER(IntKi), PARAMETER      :: RootMzc2  = 524
   INTEGER(IntKi), PARAMETER      :: RootMxb2  = 525
   INTEGER(IntKi), PARAMETER      :: RootMyb2  = 526


     ! Blade 3 Root Loads:

   INTEGER(IntKi), PARAMETER      :: RootFxc3  = 527
   INTEGER(IntKi), PARAMETER      :: RootFyc3  = 528
   INTEGER(IntKi), PARAMETER      :: RootFzc3  = 529
   INTEGER(IntKi), PARAMETER      :: RootFxb3  = 530
   INTEGER(IntKi), PARAMETER      :: RootFyb3  = 531
   INTEGER(IntKi), PARAMETER      :: RootMxc3  = 532
   INTEGER(IntKi), PARAMETER      :: RootMyc3  = 533
   INTEGER(IntKi), PARAMETER      :: RootMzc3  = 534
   INTEGER(IntKi), PARAMETER      :: RootMxb3  = 535
   INTEGER(IntKi), PARAMETER      :: RootMyb3  = 536


     ! Blade 1 Local Span Loads:

   INTEGER(IntKi), PARAMETER      :: Spn1MLxb1 = 537
   INTEGER(IntKi), PARAMETER      :: Spn1MLyb1 = 538
   INTEGER(IntKi), PARAMETER      :: Spn1MLzb1 = 539
   INTEGER(IntKi), PARAMETER      :: Spn2MLxb1 = 540
   INTEGER(IntKi), PARAMETER      :: Spn2MLyb1 = 541
   INTEGER(IntKi), PARAMETER      :: Spn2MLzb1 = 542
   INTEGER(IntKi), PARAMETER      :: Spn3MLxb1 = 543
   INTEGER(IntKi), PARAMETER      :: Spn3MLyb1 = 544
   INTEGER(IntKi), PARAMETER      :: Spn3MLzb1 = 545
   INTEGER(IntKi), PARAMETER      :: Spn4MLxb1 = 546
   INTEGER(IntKi), PARAMETER      :: Spn4MLyb1 = 547
   INTEGER(IntKi), PARAMETER      :: Spn4MLzb1 = 548
   INTEGER(IntKi), PARAMETER      :: Spn5MLxb1 = 549
   INTEGER(IntKi), PARAMETER      :: Spn5MLyb1 = 550
   INTEGER(IntKi), PARAMETER      :: Spn5MLzb1 = 551
   INTEGER(IntKi), PARAMETER      :: Spn6MLxb1 = 552
   INTEGER(IntKi), PARAMETER      :: Spn6MLyb1 = 553
   INTEGER(IntKi), PARAMETER      :: Spn6MLzb1 = 554
   INTEGER(IntKi), PARAMETER      :: Spn7MLxb1 = 555
   INTEGER(IntKi), PARAMETER      :: Spn7MLyb1 = 556
   INTEGER(IntKi), PARAMETER      :: Spn7MLzb1 = 557
   INTEGER(IntKi), PARAMETER      :: Spn8MLxb1 = 558
   INTEGER(IntKi), PARAMETER      :: Spn8MLyb1 = 559
   INTEGER(IntKi), PARAMETER      :: Spn8MLzb1 = 560
   INTEGER(IntKi), PARAMETER      :: Spn9MLxb1 = 561
   INTEGER(IntKi), PARAMETER      :: Spn9MLyb1 = 562
   INTEGER(IntKi), PARAMETER      :: Spn9MLzb1 = 563
   INTEGER(IntKi), PARAMETER      :: Spn1FLxb1 = 564
   INTEGER(IntKi), PARAMETER      :: Spn1FLyb1 = 565
   INTEGER(IntKi), PARAMETER      :: Spn1FLzb1 = 566
   INTEGER(IntKi), PARAMETER      :: Spn2FLxb1 = 567
   INTEGER(IntKi), PARAMETER      :: Spn2FLyb1 = 568
   INTEGER(IntKi), PARAMETER      :: Spn2FLzb1 = 569
   INTEGER(IntKi), PARAMETER      :: Spn3FLxb1 = 570
   INTEGER(IntKi), PARAMETER      :: Spn3FLyb1 = 571
   INTEGER(IntKi), PARAMETER      :: Spn3FLzb1 = 572
   INTEGER(IntKi), PARAMETER      :: Spn4FLxb1 = 573
   INTEGER(IntKi), PARAMETER      :: Spn4FLyb1 = 574
   INTEGER(IntKi), PARAMETER      :: Spn4FLzb1 = 575
   INTEGER(IntKi), PARAMETER      :: Spn5FLxb1 = 576
   INTEGER(IntKi), PARAMETER      :: Spn5FLyb1 = 577
   INTEGER(IntKi), PARAMETER      :: Spn5FLzb1 = 578
   INTEGER(IntKi), PARAMETER      :: Spn6FLxb1 = 579
   INTEGER(IntKi), PARAMETER      :: Spn6FLyb1 = 580
   INTEGER(IntKi), PARAMETER      :: Spn6FLzb1 = 581
   INTEGER(IntKi), PARAMETER      :: Spn7FLxb1 = 582
   INTEGER(IntKi), PARAMETER      :: Spn7FLyb1 = 583
   INTEGER(IntKi), PARAMETER      :: Spn7FLzb1 = 584
   INTEGER(IntKi), PARAMETER      :: Spn8FLxb1 = 585
   INTEGER(IntKi), PARAMETER      :: Spn8FLyb1 = 586
   INTEGER(IntKi), PARAMETER      :: Spn8FLzb1 = 587
   INTEGER(IntKi), PARAMETER      :: Spn9FLxb1 = 588
   INTEGER(IntKi), PARAMETER      :: Spn9FLyb1 = 589
   INTEGER(IntKi), PARAMETER      :: Spn9FLzb1 = 590


     ! Blade 2 Local Span Loads:

   INTEGER(IntKi), PARAMETER      :: Spn1MLxb2 = 591
   INTEGER(IntKi), PARAMETER      :: Spn1MLyb2 = 592
   INTEGER(IntKi), PARAMETER      :: Spn1MLzb2 = 593
   INTEGER(IntKi), PARAMETER      :: Spn2MLxb2 = 594
   INTEGER(IntKi), PARAMETER      :: Spn2MLyb2 = 595
   INTEGER(IntKi), PARAMETER      :: Spn2MLzb2 = 596
   INTEGER(IntKi), PARAMETER      :: Spn3MLxb2 = 597
   INTEGER(IntKi), PARAMETER      :: Spn3MLyb2 = 598
   INTEGER(IntKi), PARAMETER      :: Spn3MLzb2 = 599
   INTEGER(IntKi), PARAMETER      :: Spn4MLxb2 = 600
   INTEGER(IntKi), PARAMETER      :: Spn4MLyb2 = 601
   INTEGER(IntKi), PARAMETER      :: Spn4MLzb2 = 602
   INTEGER(IntKi), PARAMETER      :: Spn5MLxb2 = 603
   INTEGER(IntKi), PARAMETER      :: Spn5MLyb2 = 604
   INTEGER(IntKi), PARAMETER      :: Spn5MLzb2 = 605
   INTEGER(IntKi), PARAMETER      :: Spn6MLxb2 = 606
   INTEGER(IntKi), PARAMETER      :: Spn6MLyb2 = 607
   INTEGER(IntKi), PARAMETER      :: Spn6MLzb2 = 608
   INTEGER(IntKi), PARAMETER      :: Spn7MLxb2 = 609
   INTEGER(IntKi), PARAMETER      :: Spn7MLyb2 = 610
   INTEGER(IntKi), PARAMETER      :: Spn7MLzb2 = 611
   INTEGER(IntKi), PARAMETER      :: Spn8MLxb2 = 612
   INTEGER(IntKi), PARAMETER      :: Spn8MLyb2 = 613
   INTEGER(IntKi), PARAMETER      :: Spn8MLzb2 = 614
   INTEGER(IntKi), PARAMETER      :: Spn9MLxb2 = 615
   INTEGER(IntKi), PARAMETER      :: Spn9MLyb2 = 616
   INTEGER(IntKi), PARAMETER      :: Spn9MLzb2 = 617
   INTEGER(IntKi), PARAMETER      :: Spn1FLxb2 = 618
   INTEGER(IntKi), PARAMETER      :: Spn1FLyb2 = 619
   INTEGER(IntKi), PARAMETER      :: Spn1FLzb2 = 620
   INTEGER(IntKi), PARAMETER      :: Spn2FLxb2 = 621
   INTEGER(IntKi), PARAMETER      :: Spn2FLyb2 = 622
   INTEGER(IntKi), PARAMETER      :: Spn2FLzb2 = 623
   INTEGER(IntKi), PARAMETER      :: Spn3FLxb2 = 624
   INTEGER(IntKi), PARAMETER      :: Spn3FLyb2 = 625
   INTEGER(IntKi), PARAMETER      :: Spn3FLzb2 = 626
   INTEGER(IntKi), PARAMETER      :: Spn4FLxb2 = 627
   INTEGER(IntKi), PARAMETER      :: Spn4FLyb2 = 628
   INTEGER(IntKi), PARAMETER      :: Spn4FLzb2 = 629
   INTEGER(IntKi), PARAMETER      :: Spn5FLxb2 = 630
   INTEGER(IntKi), PARAMETER      :: Spn5FLyb2 = 631
   INTEGER(IntKi), PARAMETER      :: Spn5FLzb2 = 632
   INTEGER(IntKi), PARAMETER      :: Spn6FLxb2 = 633
   INTEGER(IntKi), PARAMETER      :: Spn6FLyb2 = 634
   INTEGER(IntKi), PARAMETER      :: Spn6FLzb2 = 635
   INTEGER(IntKi), PARAMETER      :: Spn7FLxb2 = 636
   INTEGER(IntKi), PARAMETER      :: Spn7FLyb2 = 637
   INTEGER(IntKi), PARAMETER      :: Spn7FLzb2 = 638
   INTEGER(IntKi), PARAMETER      :: Spn8FLxb2 = 639
   INTEGER(IntKi), PARAMETER      :: Spn8FLyb2 = 640
   INTEGER(IntKi), PARAMETER      :: Spn8FLzb2 = 641
   INTEGER(IntKi), PARAMETER      :: Spn9FLxb2 = 642
   INTEGER(IntKi), PARAMETER      :: Spn9FLyb2 = 643
   INTEGER(IntKi), PARAMETER      :: Spn9FLzb2 = 644


     ! Blade 3 Local Span Loads:

   INTEGER(IntKi), PARAMETER      :: Spn1MLxb3 = 645
   INTEGER(IntKi), PARAMETER      :: Spn1MLyb3 = 646
   INTEGER(IntKi), PARAMETER      :: Spn1MLzb3 = 647
   INTEGER(IntKi), PARAMETER      :: Spn2MLxb3 = 648
   INTEGER(IntKi), PARAMETER      :: Spn2MLyb3 = 649
   INTEGER(IntKi), PARAMETER      :: Spn2MLzb3 = 650
   INTEGER(IntKi), PARAMETER      :: Spn3MLxb3 = 651
   INTEGER(IntKi), PARAMETER      :: Spn3MLyb3 = 652
   INTEGER(IntKi), PARAMETER      :: Spn3MLzb3 = 653
   INTEGER(IntKi), PARAMETER      :: Spn4MLxb3 = 654
   INTEGER(IntKi), PARAMETER      :: Spn4MLyb3 = 655
   INTEGER(IntKi), PARAMETER      :: Spn4MLzb3 = 656
   INTEGER(IntKi), PARAMETER      :: Spn5MLxb3 = 657
   INTEGER(IntKi), PARAMETER      :: Spn5MLyb3 = 658
   INTEGER(IntKi), PARAMETER      :: Spn5MLzb3 = 659
   INTEGER(IntKi), PARAMETER      :: Spn6MLxb3 = 660
   INTEGER(IntKi), PARAMETER      :: Spn6MLyb3 = 661
   INTEGER(IntKi), PARAMETER      :: Spn6MLzb3 = 662
   INTEGER(IntKi), PARAMETER      :: Spn7MLxb3 = 663
   INTEGER(IntKi), PARAMETER      :: Spn7MLyb3 = 664
   INTEGER(IntKi), PARAMETER      :: Spn7MLzb3 = 665
   INTEGER(IntKi), PARAMETER      :: Spn8MLxb3 = 666
   INTEGER(IntKi), PARAMETER      :: Spn8MLyb3 = 667
   INTEGER(IntKi), PARAMETER      :: Spn8MLzb3 = 668
   INTEGER(IntKi), PARAMETER      :: Spn9MLxb3 = 669
   INTEGER(IntKi), PARAMETER      :: Spn9MLyb3 = 670
   INTEGER(IntKi), PARAMETER      :: Spn9MLzb3 = 671
   INTEGER(IntKi), PARAMETER      :: Spn1FLxb3 = 672
   INTEGER(IntKi), PARAMETER      :: Spn1FLyb3 = 673
   INTEGER(IntKi), PARAMETER      :: Spn1FLzb3 = 674
   INTEGER(IntKi), PARAMETER      :: Spn2FLxb3 = 675
   INTEGER(IntKi), PARAMETER      :: Spn2FLyb3 = 676
   INTEGER(IntKi), PARAMETER      :: Spn2FLzb3 = 677
   INTEGER(IntKi), PARAMETER      :: Spn3FLxb3 = 678
   INTEGER(IntKi), PARAMETER      :: Spn3FLyb3 = 679
   INTEGER(IntKi), PARAMETER      :: Spn3FLzb3 = 680
   INTEGER(IntKi), PARAMETER      :: Spn4FLxb3 = 681
   INTEGER(IntKi), PARAMETER      :: Spn4FLyb3 = 682
   INTEGER(IntKi), PARAMETER      :: Spn4FLzb3 = 683
   INTEGER(IntKi), PARAMETER      :: Spn5FLxb3 = 684
   INTEGER(IntKi), PARAMETER      :: Spn5FLyb3 = 685
   INTEGER(IntKi), PARAMETER      :: Spn5FLzb3 = 686
   INTEGER(IntKi), PARAMETER      :: Spn6FLxb3 = 687
   INTEGER(IntKi), PARAMETER      :: Spn6FLyb3 = 688
   INTEGER(IntKi), PARAMETER      :: Spn6FLzb3 = 689
   INTEGER(IntKi), PARAMETER      :: Spn7FLxb3 = 690
   INTEGER(IntKi), PARAMETER      :: Spn7FLyb3 = 691
   INTEGER(IntKi), PARAMETER      :: Spn7FLzb3 = 692
   INTEGER(IntKi), PARAMETER      :: Spn8FLxb3 = 693
   INTEGER(IntKi), PARAMETER      :: Spn8FLyb3 = 694
   INTEGER(IntKi), PARAMETER      :: Spn8FLzb3 = 695
   INTEGER(IntKi), PARAMETER      :: Spn9FLxb3 = 696
   INTEGER(IntKi), PARAMETER      :: Spn9FLyb3 = 697
   INTEGER(IntKi), PARAMETER      :: Spn9FLzb3 = 698


     ! Hub and Rotor Loads:

   INTEGER(IntKi), PARAMETER      :: LSShftFxa = 699
   INTEGER(IntKi), PARAMETER      :: LSShftFya = 700
   INTEGER(IntKi), PARAMETER      :: LSShftFza = 701
   INTEGER(IntKi), PARAMETER      :: LSShftFys = 702
   INTEGER(IntKi), PARAMETER      :: LSShftFzs = 703
   INTEGER(IntKi), PARAMETER      :: LSShftMxa = 704
   INTEGER(IntKi), PARAMETER      :: LSSTipMya = 705
   INTEGER(IntKi), PARAMETER      :: LSSTipMza = 706
   INTEGER(IntKi), PARAMETER      :: LSSTipMys = 707
   INTEGER(IntKi), PARAMETER      :: LSSTipMzs = 708
   INTEGER(IntKi), PARAMETER      :: RotPwr    = 709


     ! Shaft Strain Gage Loads:

   INTEGER(IntKi), PARAMETER      :: LSSGagMya = 710
   INTEGER(IntKi), PARAMETER      :: LSSGagMza = 711
   INTEGER(IntKi), PARAMETER      :: LSSGagMys = 712
   INTEGER(IntKi), PARAMETER      :: LSSGagMzs = 713


     ! High-Speed Shaft Loads:

   INTEGER(IntKi), PARAMETER      :: HSShftTq  = 714
   INTEGER(IntKi), PARAMETER      :: HSSBrTq   = 715
   INTEGER(IntKi), PARAMETER      :: HSShftPwr = 716


     ! Rotor-Furl Bearing Loads:

   INTEGER(IntKi), PARAMETER      :: RFrlBrM   = 717


     ! Tail-Furl Bearing Loads:

   INTEGER(IntKi), PARAMETER      :: TFrlBrM   = 718


     ! Tower-Top / Yaw Bearing Loads:

   INTEGER(IntKi), PARAMETER      :: YawBrFxn  = 719
   INTEGER(IntKi), PARAMETER      :: YawBrFyn  = 720
   INTEGER(IntKi), PARAMETER      :: YawBrFzn  = 721
   INTEGER(IntKi), PARAMETER      :: YawBrFxp  = 722
   INTEGER(IntKi), PARAMETER      :: YawBrFyp  = 723
   INTEGER(IntKi), PARAMETER      :: YawBrMxn  = 724
   INTEGER(IntKi), PARAMETER      :: YawBrMyn  = 725
   INTEGER(IntKi), PARAMETER      :: YawBrMzn  = 726
   INTEGER(IntKi), PARAMETER      :: YawBrMxp  = 727
   INTEGER(IntKi), PARAMETER      :: YawBrMyp  = 728


     ! Tower Base Loads:

   INTEGER(IntKi), PARAMETER      :: TwrBsFxt  = 729
   INTEGER(IntKi), PARAMETER      :: TwrBsFyt  = 730
   INTEGER(IntKi), PARAMETER      :: TwrBsFzt  = 731
   INTEGER(IntKi), PARAMETER      :: TwrBsMxt  = 732
   INTEGER(IntKi), PARAMETER      :: TwrBsMyt  = 733
   INTEGER(IntKi), PARAMETER      :: TwrBsMzt  = 734


     ! Local Tower Loads:

   INTEGER(IntKi), PARAMETER      :: TwHt1MLxt = 735
   INTEGER(IntKi), PARAMETER      :: TwHt1MLyt = 736
   INTEGER(IntKi), PARAMETER      :: TwHt1MLzt = 737
   INTEGER(IntKi), PARAMETER      :: TwHt2MLxt = 738
   INTEGER(IntKi), PARAMETER      :: TwHt2MLyt = 739
   INTEGER(IntKi), PARAMETER      :: TwHt2MLzt = 740
   INTEGER(IntKi), PARAMETER      :: TwHt3MLxt = 741
   INTEGER(IntKi), PARAMETER      :: TwHt3MLyt = 742
   INTEGER(IntKi), PARAMETER      :: TwHt3MLzt = 743
   INTEGER(IntKi), PARAMETER      :: TwHt4MLxt = 744
   INTEGER(IntKi), PARAMETER      :: TwHt4MLyt = 745
   INTEGER(IntKi), PARAMETER      :: TwHt4MLzt = 746
   INTEGER(IntKi), PARAMETER      :: TwHt5MLxt = 747
   INTEGER(IntKi), PARAMETER      :: TwHt5MLyt = 748
   INTEGER(IntKi), PARAMETER      :: TwHt5MLzt = 749
   INTEGER(IntKi), PARAMETER      :: TwHt6MLxt = 750
   INTEGER(IntKi), PARAMETER      :: TwHt6MLyt = 751
   INTEGER(IntKi), PARAMETER      :: TwHt6MLzt = 752
   INTEGER(IntKi), PARAMETER      :: TwHt7MLxt = 753
   INTEGER(IntKi), PARAMETER      :: TwHt7MLyt = 754
   INTEGER(IntKi), PARAMETER      :: TwHt7MLzt = 755
   INTEGER(IntKi), PARAMETER      :: TwHt8MLxt = 756
   INTEGER(IntKi), PARAMETER      :: TwHt8MLyt = 757
   INTEGER(IntKi), PARAMETER      :: TwHt8MLzt = 758
   INTEGER(IntKi), PARAMETER      :: TwHt9MLxt = 759
   INTEGER(IntKi), PARAMETER      :: TwHt9MLyt = 760
   INTEGER(IntKi), PARAMETER      :: TwHt9MLzt = 761
   INTEGER(IntKi), PARAMETER      :: TwHt1FLxt = 762
   INTEGER(IntKi), PARAMETER      :: TwHt1FLyt = 763
   INTEGER(IntKi), PARAMETER      :: TwHt1FLzt = 764
   INTEGER(IntKi), PARAMETER      :: TwHt2FLxt = 765
   INTEGER(IntKi), PARAMETER      :: TwHt2FLyt = 766
   INTEGER(IntKi), PARAMETER      :: TwHt2FLzt = 767
   INTEGER(IntKi), PARAMETER      :: TwHt3FLxt = 768
   INTEGER(IntKi), PARAMETER      :: TwHt3FLyt = 769
   INTEGER(IntKi), PARAMETER      :: TwHt3FLzt = 770
   INTEGER(IntKi), PARAMETER      :: TwHt4FLxt = 771
   INTEGER(IntKi), PARAMETER      :: TwHt4FLyt = 772
   INTEGER(IntKi), PARAMETER      :: TwHt4FLzt = 773
   INTEGER(IntKi), PARAMETER      :: TwHt5FLxt = 774
   INTEGER(IntKi), PARAMETER      :: TwHt5FLyt = 775
   INTEGER(IntKi), PARAMETER      :: TwHt5FLzt = 776
   INTEGER(IntKi), PARAMETER      :: TwHt6FLxt = 777
   INTEGER(IntKi), PARAMETER      :: TwHt6FLyt = 778
   INTEGER(IntKi), PARAMETER      :: TwHt6FLzt = 779
   INTEGER(IntKi), PARAMETER      :: TwHt7FLxt = 780
   INTEGER(IntKi), PARAMETER      :: TwHt7FLyt = 781
   INTEGER(IntKi), PARAMETER      :: TwHt7FLzt = 782
   INTEGER(IntKi), PARAMETER      :: TwHt8FLxt = 783
   INTEGER(IntKi), PARAMETER      :: TwHt8FLyt = 784
   INTEGER(IntKi), PARAMETER      :: TwHt8FLzt = 785
   INTEGER(IntKi), PARAMETER      :: TwHt9FLxt = 786
   INTEGER(IntKi), PARAMETER      :: TwHt9FLyt = 787
   INTEGER(IntKi), PARAMETER      :: TwHt9FLzt = 788


     ! Internal Degrees of Freedom:

   INTEGER(IntKi), PARAMETER      :: Q_B1E1    = 789
   INTEGER(IntKi), PARAMETER      :: Q_B2E1    = 790
   INTEGER(IntKi), PARAMETER      :: Q_B3E1    = 791
   INTEGER(IntKi), PARAMETER      :: Q_B1F1    = 792
   INTEGER(IntKi), PARAMETER      :: Q_B2F1    = 793
   INTEGER(IntKi), PARAMETER      :: Q_B3F1    = 794
   INTEGER(IntKi), PARAMETER      :: Q_B1F2    = 795
   INTEGER(IntKi), PARAMETER      :: Q_B2F2    = 796
   INTEGER(IntKi), PARAMETER      :: Q_B3F2    = 797
   INTEGER(IntKi), PARAMETER      :: Q_Teet    = 798
   INTEGER(IntKi), PARAMETER      :: Q_DrTr    = 799
   INTEGER(IntKi), PARAMETER      :: Q_GeAz    = 800
   INTEGER(IntKi), PARAMETER      :: Q_RFrl    = 801
   INTEGER(IntKi), PARAMETER      :: Q_TFrl    = 802
   INTEGER(IntKi), PARAMETER      :: Q_Yaw     = 803
   INTEGER(IntKi), PARAMETER      :: Q_TFA1    = 804
   INTEGER(IntKi), PARAMETER      :: Q_TSS1    = 805
   INTEGER(IntKi), PARAMETER      :: Q_TFA2    = 806
   INTEGER(IntKi), PARAMETER      :: Q_TSS2    = 807
   INTEGER(IntKi), PARAMETER      :: Q_Sg      = 808
   INTEGER(IntKi), PARAMETER      :: Q_Sw      = 809
   INTEGER(IntKi), PARAMETER      :: Q_Hv      = 810
   INTEGER(IntKi), PARAMETER      :: Q_R       = 811
   INTEGER(IntKi), PARAMETER      :: Q_P       = 812
   INTEGER(IntKi), PARAMETER      :: Q_Y       = 813
   INTEGER(IntKi), PARAMETER      :: QD_B1E1   = 814
   INTEGER(IntKi), PARAMETER      :: QD_B2E1   = 815
   INTEGER(IntKi), PARAMETER      :: QD_B3E1   = 816
   INTEGER(IntKi), PARAMETER      :: QD_B1F1   = 817
   INTEGER(IntKi), PARAMETER      :: QD_B2F1   = 818
   INTEGER(IntKi), PARAMETER      :: QD_B3F1   = 819
   INTEGER(IntKi), PARAMETER      :: QD_B1F2   = 820
   INTEGER(IntKi), PARAMETER      :: QD_B2F2   = 821
   INTEGER(IntKi), PARAMETER      :: QD_B3F2   = 822
   INTEGER(IntKi), PARAMETER      :: QD_Teet   = 823
   INTEGER(IntKi), PARAMETER      :: QD_DrTr   = 824
   INTEGER(IntKi), PARAMETER      :: QD_GeAz   = 825
   INTEGER(IntKi), PARAMETER      :: QD_RFrl   = 826
   INTEGER(IntKi), PARAMETER      :: QD_TFrl   = 827
   INTEGER(IntKi), PARAMETER      :: QD_Yaw    = 828
   INTEGER(IntKi), PARAMETER      :: QD_TFA1   = 829
   INTEGER(IntKi), PARAMETER      :: QD_TSS1   = 830
   INTEGER(IntKi), PARAMETER      :: QD_TFA2   = 831
   INTEGER(IntKi), PARAMETER      :: QD_TSS2   = 832
   INTEGER(IntKi), PARAMETER      :: QD_Sg     = 833
   INTEGER(IntKi), PARAMETER      :: QD_Sw     = 834
   INTEGER(IntKi), PARAMETER      :: QD_Hv     = 835
   INTEGER(IntKi), PARAMETER      :: QD_R      = 836
   INTEGER(IntKi), PARAMETER      :: QD_P      = 837
   INTEGER(IntKi), PARAMETER      :: QD_Y      = 838
   INTEGER(IntKi), PARAMETER      :: QD2_B1E1  = 839
   INTEGER(IntKi), PARAMETER      :: QD2_B2E1  = 840
   INTEGER(IntKi), PARAMETER      :: QD2_B3E1  = 841
   INTEGER(IntKi), PARAMETER      :: QD2_B1F1  = 842
   INTEGER(IntKi), PARAMETER      :: QD2_B2F1  = 843
   INTEGER(IntKi), PARAMETER      :: QD2_B3F1  = 844
   INTEGER(IntKi), PARAMETER      :: QD2_B1F2  = 845
   INTEGER(IntKi), PARAMETER      :: QD2_B2F2  = 846
   INTEGER(IntKi), PARAMETER      :: QD2_B3F2  = 847
   INTEGER(IntKi), PARAMETER      :: QD2_Teet  = 848
   INTEGER(IntKi), PARAMETER      :: QD2_DrTr  = 849
   INTEGER(IntKi), PARAMETER      :: QD2_GeAz  = 850
   INTEGER(IntKi), PARAMETER      :: QD2_RFrl  = 851
   INTEGER(IntKi), PARAMETER      :: QD2_TFrl  = 852
   INTEGER(IntKi), PARAMETER      :: QD2_Yaw   = 853
   INTEGER(IntKi), PARAMETER      :: QD2_TFA1  = 854
   INTEGER(IntKi), PARAMETER      :: QD2_TSS1  = 855
   INTEGER(IntKi), PARAMETER      :: QD2_TFA2  = 856
   INTEGER(IntKi), PARAMETER      :: QD2_TSS2  = 857
   INTEGER(IntKi), PARAMETER      :: QD2_Sg    = 858
   INTEGER(IntKi), PARAMETER      :: QD2_Sw    = 859
   INTEGER(IntKi), PARAMETER      :: QD2_Hv    = 860
   INTEGER(IntKi), PARAMETER      :: QD2_R     = 861
   INTEGER(IntKi), PARAMETER      :: QD2_P     = 862
   INTEGER(IntKi), PARAMETER      :: QD2_Y     = 863


     ! The maximum number of output channels which can be output by the code.
   INTEGER(IntKi), PARAMETER      :: MaxOutPts = 863

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
!> This module contains file I/O routines and data validation routines.
MODULE ElastoDyn_IO

   USE ElastoDyn_Parameters
   USE ElastoDyn_Types

   IMPLICIT NONE
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine reads the input file and stores all the data in the ED_InputFile structure.
!! It does not perform data validation.
SUBROUTINE ED_ReadInput( InputFileName, MeshFile, InputFileData, BD4Blades, Default_DT, OutFileRoot, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   REAL(DbKi),           INTENT(IN)       :: Default_DT     !< The default DT (from glue code)

   CHARACTER(*), INTENT(IN)               :: InputFileName  !< Name of the input file
   CHARACTER(*), INTENT(IN)               :: MeshFile       !< File that contains the blade mesh information (AeroDyn input file for now) -- later this info will be defined in one of the ED input files.
   CHARACTER(*), INTENT(IN)               :: OutFileRoot    !< The rootname of all the output files written by this routine.

   TYPE(ED_InputFile),   INTENT(OUT)      :: InputFileData  !< Data stored in the module's input file

   INTEGER(IntKi),       INTENT(OUT)      :: ErrStat        !< The error status code
   LOGICAL,              INTENT(IN)       :: BD4Blades      !< Determines if we should read the blade values (true=don't read this file; use BeamDyn for blades instead)
   CHARACTER(*),         INTENT(OUT)      :: ErrMsg         !< The error message, if an error occurred

      ! local variables

   INTEGER(IntKi)                         :: UnEcho         !  Unit number for the echo file
   INTEGER(IntKi)                         :: ErrStat2       !  The error status code
   CHARACTER(ErrMsgLen)                   :: ErrMsg2        !  The error message, if an error occurred
   CHARACTER(*), PARAMETER                :: RoutineName = 'ED_ReadInput'                
   
   CHARACTER(1024)                        :: BldFile(MaxBl) !  File that contains the blade information (specified in the primary input file)
   CHARACTER(1024)                        :: FurlFile       !  File that contains the furl information (specified in the primary input file)
   CHARACTER(1024)                        :: TwrFile        !  File that contains the tower information (specified in the primary input file)

      ! initialize values:

   ErrStat = ErrID_None
   ErrMsg  = ''

   InputFileData%DT = Default_DT  ! the glue code's suggested DT for the module (may be overwritten in ReadPrimaryFile())

      ! get the primary/platform input-file data
      ! sets UnEcho, BldFile, FurlFile, TwrFile
   
   CALL ReadPrimaryFile( InputFileName, InputFileData, BldFile, FurlFile, TwrFile, OutFileRoot, UnEcho, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call Cleanup()
         return
      end if

      ! get the furling input-file data
   IF ( InputFileData%Furling )  THEN
      CALL ReadFurlFile( FurlFile, InputFileData, UnEcho, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call Cleanup()
            return
         end if         
   ELSE   ! initialize all of the data that would be read by ReadFurlFile()
      InputFileData%RFrlDOF   = .FALSE.
      InputFileData%TFrlDOF   = .FALSE.
      InputFileData%RotFurl   = 0.0_ReKi  ! Radians
      InputFileData%TailFurl  = 0.0_ReKi
      InputFileData%Yaw2Shft  = 0.0
      InputFileData%ShftSkew  = 0.0
      InputFileData%RFrlCM_n  = 0.0
      InputFileData%BoomCM_n  = 0.0
      InputFileData%TFinCM_n  = 0.0
      InputFileData%RFrlPnt_n = 0.0
      InputFileData%RFrlSkew  = 0.0
      InputFileData%RFrlTilt  = 0.0
      InputFileData%TFrlPnt_n = 0.0
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
      CALL ReadBladeInputs ( BldFile, MeshFile, InputFileData, UnEcho, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            call Cleanup()
            return
         end if
   END IF

      ! get the tower input-file data

   CALL ReadTowerFile( TwrFile, InputFileData, UnEcho,  ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      if ( ErrStat >= AbortErrLev ) then
         call Cleanup()
         return
      end if
      


      ! close the echo file (if opened)

   CALL Cleanup()


CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()   
      IF ( UnEcho > 0 ) CLOSE( UnEcho )
   END SUBROUTINE Cleanup

END SUBROUTINE ED_ReadInput
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine validates the input file data
SUBROUTINE ED_ValidateInput( InputFileData, BD4Blades, Linearize, MHK, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_InputFile),       INTENT(IN)       :: InputFileData       !< Data stored in the module's input file
   LOGICAL,                  INTENT(IN)       :: BD4Blades           !< Determines if we should validate the blade values (true=don't validate; use BeamDyn for blades instead)
   LOGICAL,                  INTENT(IN)       :: Linearize           !< Flag indicating glue code wants to linearize this module
   INTEGER(IntKi),           INTENT(IN)       :: MHK                 !< MHK turbine type switch
   INTEGER(IntKi),           INTENT(OUT)      :: ErrStat             !< The error status code
   CHARACTER(*),             INTENT(OUT)      :: ErrMsg              !< The error message, if an error occurred

      ! Local variables:
   INTEGER(IntKi)                             :: K                   ! Blade number
   INTEGER(IntKi)                             :: ErrStat2            ! Temporary error ID
   CHARACTER(ErrMsgLen)                       :: ErrMsg2             ! Temporary message describing error
   CHARACTER(*), PARAMETER                    :: RoutineName = 'ED_ValidateInput'

      ! Initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ''

      ! validate the primary input data
   CALL ValidatePrimaryData( InputFileData, BD4Blades, Linearize, MHK, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


      ! validate the furling input data
   CALL ValidateFurlData ( InputFileData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)


      ! validate the blade input data
   IF (.NOT. BD4Blades) THEN
      DO K = 1,InputFileData%NumBl
         CALL ValidateBladeData ( InputFileData%InpBl(K), ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName//':Blade'//TRIM(Num2LStr(K)))
      END DO
      
      !bjj: validate blade discretization, too:
   END IF
   

      ! validate the tower input data
   CALL ValidateTowerData ( InputFileData, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      
      ! validate the Output parameters:
  ! CALL ChckOutLst( InputFileData%OutList, p, ErrStat, ErrMsg )

END SUBROUTINE ED_ValidateInput
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads the data from the blade and mesh inputs files.
!! This routines assumes that InputFileData%NumBl has already been set.
SUBROUTINE ReadBladeInputs ( BldFile, MeshFile, InputFileData, UnEc, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables:

   TYPE(ED_InputFile),     INTENT(INOUT)  :: InputFileData                       !< Input file data Data for Blade K stored in the module's input file
   CHARACTER(*),           INTENT(IN)     :: BldFile(:)                          !< The array of file names containing blade information
   CHARACTER(*),           INTENT(IN)     :: MeshFile                            !< The file names containing blade mesh information (for now, the aerodyn primary file)
   INTEGER(IntKi),         INTENT(IN)     :: UnEc                                !< I/O unit for echo file. If present and > 0, write to UnEc

   INTEGER(IntKi),         INTENT(OUT)    :: ErrStat                             !< The error ID
   CHARACTER(*),           INTENT(OUT)    :: ErrMsg                              !< Message describing error


      ! Local variables:
   INTEGER(IntKi)                         :: K                                   ! Blade number
   INTEGER(IntKi)                         :: ErrStat2                            ! Temporary error ID
   LOGICAL                                :: ReadFile                            ! determines if an input file for a blade is the same as the file for the previous blade
   CHARACTER(ErrMsgLen)                   :: ErrMsg2                             ! Temporary message describing error
   CHARACTER(*), PARAMETER                :: RoutineName = 'ReadBladeInputs'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''


      ! Allocate space for the input file data
   ALLOCATE( InputFileData%InpBlMesh( 1_IntKi ), STAT=ErrStat2 )              ! for now, we're assuming the discretization is the same on all blades
   IF ( ErrStat2 /= 0 ) THEN
      call SetErrStat(ErrID_Fatal,'Error allocating InpBlMesh array',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF

   ALLOCATE( InputFileData%InpBl( InputFileData%NumBl ), STAT=ErrStat2 )
   IF ( ErrStat2 /= 0 ) THEN
      call SetErrStat(ErrID_Fatal,'Error allocating InpBl array',ErrStat,ErrMsg,RoutineName)
      RETURN
   END IF


      ! Get the blade discretization here:
   IF ( len_trim(MeshFile) == 0 ) THEN
      InputFileData%InpBlMesh(1)%BldNodes = InputFileData%BldNodes
   ELSE
         ! we will get the discretization from AeroDyn's input file
      CALL ReadBladeMeshFileAD( InputFileData%InpBlMesh(1), MeshFile, UnEc, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
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

         CALL ReadBladeFile( BldFile(K), InputFileData%InpBl(K), UnEc, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

      ELSE
         CALL ED_CopyBladeInputData( InputFileData%InpBl(K-1), InputFileData%InpBl(K), MESH_UPDATECOPY, ErrStat2, ErrMsg2 )
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN
               ! bjj: we could just read the file again...

      END IF

         ! If the next file is the same as this one, don't read it again:

      IF ( K /= InputFileData%NumBl ) ReadFile = BldFile(K) /= BldFile( K + 1 )

   END DO


   RETURN

END SUBROUTINE ReadBladeInputs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads a blade input file.
SUBROUTINE ReadBladeFile ( BldFile, BladeKInputFileData, UnEc, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables:

   TYPE(BladeInputData),     INTENT(INOUT)  :: BladeKInputFileData                 !< Data for Blade K stored in the module's input file
   CHARACTER(*),             INTENT(IN)     :: BldFile                             !< Name of the blade input file data
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                !< I/O unit for echo file. If present and > 0, write to UnEc

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              !< Error message


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
   CHARACTER(*), PARAMETER      :: RoutineName='ReadBladeFile'

   ErrStat = ErrID_None
   ErrMsg = ""
   
   UnIn = -1
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )


      ! Open the input file for blade K.

   CALL OpenFInpFile ( UnIn, BldFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      

   !  -------------- HEADER -------------------------------------------------------

      ! Skip the header.

   CALL ReadCom ( UnIn, BldFile, 'unused blade file header line 1', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL ReadCom ( UnIn, BldFile, 'unused blade file header line 2', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


   !  -------------- BLADE PARAMETERS ---------------------------------------------

      ! Skip the comment line.

   CALL ReadCom ( UnIn, BldFile, 'blade parameters', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NBlInpSt - Number of blade input stations.

   CALL ReadVar ( UnIn, BldFile, BladeKInputFileData%NBlInpSt, 'NBlInpSt', 'Number of blade input stations', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! .......... Allocate the arrays based on this NBlInpSt input ..........
   CALL Alloc_BladeInputProperties( BladeKInputFileData, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF



      ! BldFlDmp - Blade structural damping ratios in flapwise direction.

   CALL ReadAryLines( UnIn, BldFile, BladeKInputFileData%BldFlDmp, SIZE(BladeKInputFileData%BldFlDmp), 'BldFlDmp', &
                                       'Blade structural damping ratios in flapwise direction', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF



      ! BldEdDmp - Blade structural damping ratios in edgewise direction.

   CALL ReadAryLines( UnIn, BldFile, BladeKInputFileData%BldEdDmp, SIZE(BladeKInputFileData%BldEdDmp), 'BldEdDmp', &
                                       'Blade structural damping ratios in edgewise direction', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   !  -------------- BLADE ADJUSTMENT FACTORS -------------------------------------


      ! Skip the comment line.

   CALL ReadCom ( UnIn, BldFile, 'blade adjustment factors', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! FlStTunr(1) - Blade flapwise modal stiffness tuners.

   CALL ReadAryLines ( UnIn, BldFile, BladeKInputFileData%FlStTunr, SIZE(BladeKInputFileData%FlStTunr), 'FlStTunr', &
                                                  'Blade flapwise modal stiffness tuners', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF



      ! AdjBlMs - Factor to adjust blade mass density.

   CALL ReadVar ( UnIn, BldFile, AdjBlMs, 'AdjBlMs', 'Factor to adjust blade mass density', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! AdjFlSt - Factor to adjust blade flap stiffness.

   CALL ReadVar ( UnIn, BldFile, AdjFlSt, 'AdjFlSt', 'Factor to adjust blade flap stiffness', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! AdjEdSt - Factor to adjust blade edge stiffness.

   CALL ReadVar ( UnIn, BldFile, AdjEdSt, 'AdjEdSt', 'Factor to adjust blade edge stiffness', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


         ! Check the locally-defined adjustment factors: AdjBlMs, AdjFlSt, AdjEdSt

      IF ( AdjBlMs <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Warn, 'AdjBlMs must be greater than zero.', ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
      END IF

      IF ( AdjFlSt <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Warn, 'AdjFlSt must be greater than zero.', ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
      END IF

      IF ( AdjEdSt <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Warn, 'AdjEdSt must be greater than zero.', ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
      END IF


   !  -------------- DISTRIBUTED BLADE PROPERTIES ---------------------------------


      ! Skip the comment lines.

   CALL ReadCom ( UnIn, BldFile, 'distributed blade parameters'     , ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadCom ( UnIn, BldFile, 'distributed-blade-parameter names', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL ReadCom ( UnIn, BldFile, 'distributed-blade-parameter units', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! Read the table.

   NInputCols = 6


   DO I=1,BladeKInputFileData%NBlInpSt

      CALL ReadAry( UnIn, BldFile, TmpRAry, NInputCols, 'Line'//TRIM(Num2LStr(I)), 'Blade input station table', &
                    ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF

      BladeKInputFileData%BlFract( I) = TmpRAry(1)
      BladeKInputFileData%PitchAx( I) = TmpRAry(2)
      BladeKInputFileData%StrcTwst(I) = TmpRAry(3)*D2R      ! Input in degrees; converted to radians here
      BladeKInputFileData%BMassDen(I) = TmpRAry(4)*AdjBlMs  ! Apply the correction factors to the elemental data.
      BladeKInputFileData%FlpStff( I) = TmpRAry(5)*AdjFlSt  ! Apply the correction factors to the elemental data.
      BladeKInputFileData%EdgStff( I) = TmpRAry(6)*AdjEdSt  ! Apply the correction factors to the elemental data.

   ENDDO ! I



   !  -------------- BLADE MODE SHAPES --------------------------------------------


      ! Skip the comment line.

   CALL ReadCom ( UnIn, BldFile, 'blade mode shapes', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! BldFl1Sh - Blade-flap mode-1 shape coefficients.
   CALL ReadAryLines ( UnIn, BldFile, BladeKInputFileData%BldFl1Sh, SIZE(BladeKInputFileData%BldFl1Sh), 'BldFl1Sh', &
                           'Blade-flap mode-1 shape coefficients', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! BldFl2Sh - Blade-flap mode-2 shape coefficients.

   CALL ReadAryLines ( UnIn, BldFile, BladeKInputFileData%BldFl2Sh, SIZE(BladeKInputFileData%BldFl2Sh), 'BldFl2Sh', &
                    'Blade-flap mode-2 shape coefficients', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! BldEdgSh - Blade-edge mode shape coefficients.

   CALL ReadAryLines ( UnIn, BldFile, BladeKInputFileData%BldEdgSh, SIZE(BladeKInputFileData%BldEdgSh), 'BldEdgSh', &
                     'Blade-edge mode shape coefficients', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF



   !  -------------- END OF FILE --------------------------------------------

      ! Close the blade file.

   call Cleanup()
   RETURN


CONTAINS

   SUBROUTINE Cleanup()
      IF (UnIn > 0) CLOSE( UnIn )
   END SUBROUTINE Cleanup
   
END SUBROUTINE ReadBladeFile
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads in the AeroDyn v14.00.00 input file to get the
!! blade discretization used in the structural dynamics module.
SUBROUTINE ReadBladeMeshFileAD( BladeKInputFileMesh, MeshFile, UnEc, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables

   TYPE(ED_BladeMeshInputData),   INTENT(INOUT)  :: BladeKInputFileMesh                 !< All the data in the ElastoDyn input file
   CHARACTER(*),                  INTENT(IN)     :: MeshFile                            !< Name of the AeroDyn input file data (for mesh)

   INTEGER(IntKi),                INTENT(IN)     :: UnEc                                !< I/O unit for echo file. If present and > 0, write to UnEc
   INTEGER(IntKi),                INTENT(OUT)    :: ErrStat                             !< Error status
   CHARACTER(*),                  INTENT(OUT)    :: ErrMsg                              !< Error message

      ! Local variables:
   INTEGER(IntKi), PARAMETER    :: NInputCols = 4                                       ! Number of input columns to be read from the file
   REAL(ReKi)                   :: TmpRAry(NInputCols)                                  ! Temporary variable to read table from file
   INTEGER(IntKi)               :: I                                                    ! loop counter
   INTEGER(IntKi)               :: NumLin2Skp                                           ! number of lines to read
   INTEGER(IntKi)               :: NumFoil                                              ! number of airfoil lines to skip in the AD input file.
   INTEGER(IntKi)               :: UnIn                                                 ! Unit number for reading file

   INTEGER(IntKi)               :: ErrStat2                                             ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                              ! Temporary Err msg
   CHARACTER(*), PARAMETER      :: RoutineName = 'ReadBladeMeshFileAD'
   CHARACTER(1024)              :: Line                                                 ! Temporary string.
!   CHARACTER(1024)              :: TmpStr(1)                                            ! Temporary string.



      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the AeroDyn input file.

   CALL OpenFInpFile ( UnIn, MeshFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! Add a separator to the echo file if appropriate.

   IF ( UnEc > 0 )  WRITE (UnEc,'(//,A,/)')  'Mesh input data from (AeroDyn input) file "'//TRIM( MeshFile )//'":'


   !  -------------- HEADER -------------------------------------------------------
   ! BJJ: This file is AeroDyn's input file. Until we decide on a format for the
   ! structural dynamics input, we will get this information from AeroDyn like we
   ! used to.

   DO I = 1,9
      CALL ReadCom ( UnIn, MeshFile, 'AeroDyn input (for structural dynamics mesh)', ErrStat2, ErrMsg2  )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   END DO
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF

      ! See if the next line is "NEWTOWER".  If it is, read 7 more lines.  If not, read 5 more lines.

   CALL ReadVar( UnIn, MeshFile, Line, VarName='NewTowerModel?', VarDescr='Check for tower influence model', ErrStat=ErrStat2, ErrMsg=ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Check if this is the "special string" to indicate the new tower influence model

   CALL Conv2UC( Line )
   IF ( INDEX(Line, "NEWTOWER" ) > 0 ) THEN
      NumLin2Skp = 7
   ELSE
      NumLin2Skp = 5
   END IF

   DO I = 1,NumLin2Skp
      CALL ReadCom ( UnIn, MeshFile, 'AeroDyn input (for structural dynamics mesh)', ErrStat2, ErrMsg2  )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   END DO

   CALL ReadVar ( UnIn, MeshFile, NumFoil, 'NumFoil', &
                  'Number of airfoil lines to skip in AeroDyn input (for structural dynamics mesh)', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   DO I = 1,NumFoil
      CALL ReadCom ( UnIn, MeshFile, 'AeroDyn input (for structural dynamics mesh)', ErrStat2, ErrMsg2  )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   END DO


  !  -------------- Blade Mesh Data --------------------------------------------------

      ! Read in the number of blade elements
   CALL ReadVar( UnIn, MeshFile, BladeKInputFileMesh%BldNodes, 'BldNodes', 'Number of blade elements', ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Allocate the arrays to store input
   CALL Alloc_BladeMeshInputProperties( BladeKInputFileMesh, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Read comment line for the element table
   CALL ReadCom( UnIn, MeshFile, 'Blade element table headers', ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   DO I = 1, BladeKInputFileMesh%BldNodes

      CALL ReadAry( UnIn, MeshFile, TmpRAry, NInputCols, 'Blade element line'//TRIM(Num2LStr(I)), 'Blade element input table', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF

         BladeKInputFileMesh%RNodes(  I) = TmpRAry(1)
         BladeKInputFileMesh%AeroTwst(I) = TmpRAry(2)*D2R  !Convert input file data (degrees) to radians
         BladeKInputFileMesh%Chord(   I) = TmpRAry(4)

   END DO

      !bjj: move this to a validation routine if we plan to keep AD14 stuff in ElastoDyn:
   IF ( ANY( BladeKInputFileMesh%Chord < 0.0_ReKi ) ) THEN
      call SetErrStat( ErrID_Fatal, 'Chord length must be larger than 0 meters.', ErrStat, ErrMsg, RoutineName )
      RETURN
   END IF


      ! Close the input file:

   CALL cleanup() 
   RETURN


CONTAINS
   SUBROUTINE Cleanup()
      CLOSE( UnIn )
   END SUBROUTINE Cleanup
END SUBROUTINE ReadBladeMeshFileAD
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads the furling file input and converts units as appropriate.
SUBROUTINE ReadFurlFile( FurlFile, InputFileData, UnEc, ErrStat, ErrMsg  )
   TYPE(ED_InputFile),       INTENT(INOUT)  :: InputFileData                       !< All the data in the ElastoDyn input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             !< Error status
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                !< I/O unit for echo file. If present and > 0, write to UnEc
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              !< Error message
   CHARACTER(*),             INTENT(IN)     :: FurlFile                            !< Name of the furling input file data
   ! Local variables:
   type(FileInfoType)   :: FileInfo_In ! < The derived type for holding the file information.
   integer(IntKi)       :: iLine     !< current entry in FileInfo_In%Lines array
   character(len=1024 ) :: DummyLine
   integer(IntKi)       :: ErrStat2    !< Temporary Error status
   character(ErrMsgLen) :: ErrMsg2     !< Temporary Error message

   ! --- Read TailFurl input file into array of strings
   call ProcessComFile( FurlFile, FileInfo_In, ErrStat2, ErrMsg2)

   ! --- Parse the array of strings

   ! Skip the first two lines as they are known to be header lines and separators
   do iLine = 1,2 
      if ( UnEc>0 )   WRITE(UnEc, '(A)') FileInfo_In%Lines(iLine)    ! Write header to echo
   enddo
   iLine = 3  
   !---------------------- FEATURE FLAGS -------------------------------------------
   call ParseCom(FileInfo_in, iLine, DummyLine                           , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlDOF'  , InputFileData%RFrlDOF  , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlDOF'  , InputFileData%TFrlDOF  , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;

   ! ---------------------- INITIAL CONDITIONS --------------------------------------
   call ParseCom(FileInfo_in, iLine, DummyLine                           , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RotFurl'  , InputFileData%RotFurl  , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TailFurl' , InputFileData%TailFurl , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;

   ! ---------------------- TURBINE CONFIGURATION -----------------------------------
   call ParseCom(FileInfo_in, iLine, DummyLine                           ,  ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'Yaw2Shft' , InputFileData%Yaw2Shft ,  ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'ShftSkew' , InputFileData%ShftSkew ,  ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_in, iLine, 'RFrlCM_n' , InputFileData%RFrlCM_n ,3,ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_in, iLine, 'BoomCM_n' , InputFileData%BoomCM_n ,3,ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_in, iLine, 'TFinCM_n' , InputFileData%TFinCM_n ,3,ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_in, iLine, 'RFrlPnt_n', InputFileData%RFrlPnt_n,3,ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlSkew' , InputFileData%RFrlSkew ,  ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlTilt' , InputFileData%RFrlTilt ,  ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseAry(FileInfo_in, iLine, 'TFrlPnt_n', InputFileData%TFrlPnt_n,3,ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlSkew' , InputFileData%TFrlSkew ,  ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlTilt' , InputFileData%TFrlTilt ,  ErrStat2, ErrMsg2, UnEc); if (Failed()) return;

   ! ---------------------- MASS AND INERTIA ----------------------------------------
   call ParseCom(FileInfo_in, iLine, DummyLine                           , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlMass' , InputFileData%RFrlMass , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'BoomMass' , InputFileData%BoomMass , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFinMass' , InputFileData%TFinMass , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlIner' , InputFileData%RFrlIner , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlIner' , InputFileData%TFrlIner , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;

   ! ---------------------- ROTOR-FURL ----------------------------------------------
   call ParseCom(FileInfo_in, iLine, DummyLine                           , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlMod'  , InputFileData%RFrlMod  , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlSpr'  , InputFileData%RFrlSpr  , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlDmp'  , InputFileData%RFrlDmp  , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlUSSP' , InputFileData%RFrlUSSP , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlDSSP' , InputFileData%RFrlDSSP , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlUSSpr', InputFileData%RFrlUSSpr, ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlDSSpr', InputFileData%RFrlDSSpr, ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlUSDP' , InputFileData%RFrlUSDP , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlDSDP' , InputFileData%RFrlDSDP , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlUSDmp', InputFileData%RFrlUSDmp, ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'RFrlDSDmp', InputFileData%RFrlDSDmp, ErrStat2, ErrMsg2, UnEc); if (Failed()) return;

   ! ---------------------- TAIL-FURL -----------------------------------------------
   call ParseCom(FileInfo_in, iLine, DummyLine                           , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlMod'  , InputFileData%TFrlMod  , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlSpr'  , InputFileData%TFrlSpr  , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlDmp'  , InputFileData%TFrlDmp  , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlUSSP' , InputFileData%TFrlUSSP , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlDSSP' , InputFileData%TFrlDSSP , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlUSSpr', InputFileData%TFrlUSSpr, ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlDSSpr', InputFileData%TFrlDSSpr, ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlUSDP' , InputFileData%TFrlUSDP , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlDSDP' , InputFileData%TFrlDSDP , ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlUSDmp', InputFileData%TFrlUSDmp, ErrStat2, ErrMsg2, UnEc); if (Failed()) return;
   call ParseVar(FileInfo_in, iLine, 'TFrlDSDmp', InputFileData%TFrlDSDmp, ErrStat2, ErrMsg2, UnEc); if (Failed()) return;


   ! --- Triggers
   InputFileData%RotFurl   = InputFileData%RotFurl  * D2R
   InputFileData%TailFurl  = InputFileData%TailFurl * D2R
   InputFileData%ShftSkew  = InputFileData%ShftSkew * D2R
   InputFileData%RFrlSkew  = InputFileData%RFrlSkew * D2R
   InputFileData%RFrlTilt  = InputFileData%RFrlTilt * D2R
   InputFileData%TFrlSkew  = InputFileData%TFrlSkew * D2R
   InputFileData%TFrlTilt  = InputFileData%TFrlTilt * D2R
   InputFileData%RFrlUSSP  = InputFileData%RFrlUSSP * D2R
   InputFileData%RFrlDSSP  = InputFileData%RFrlDSSP * D2R
   InputFileData%RFrlUSDP  = InputFileData%RFrlUSDP * D2R
   InputFileData%RFrlDSDP  = InputFileData%RFrlDSDP * D2R
   InputFileData%TFrlUSSP  = InputFileData%TFrlUSSP * D2R
   InputFileData%TFrlDSSP  = InputFileData%TFrlDSSP * D2R
   InputFileData%TFrlUSDP  = InputFileData%TFrlUSDP * D2R
   InputFileData%TFrlDSDP  = InputFileData%TFrlDSDP * D2R




contains
   logical function Failed()
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadFurlFile' )
      Failed = ErrStat >= AbortErrLev
   end function Failed

   subroutine Fatal(ErrMsg_in)
      character(len=*), intent(in) :: ErrMsg_in
      call SetErrStat(ErrID_Fatal, 'File:'//trim(FurlFile)//':'//trim(ErrMsg_in), ErrStat, ErrMsg, 'ReadFurlFile')
   end subroutine Fatal

END SUBROUTINE ReadFurlFile
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads the tower file  input.
SUBROUTINE ReadTowerFile( TwrFile, InputFileData, UnEc, ErrStat, ErrMsg )
!..................................................................................................................................

   IMPLICIT                        NONE

      ! Passed variables:

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             !< Error status
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                !< I/O unit for echo file. If present and > 0, write to UnEc
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              !< Error message
   CHARACTER(*),             INTENT(IN)     :: TwrFile                             !< Name of the tower input file data
   TYPE(ED_InputFile),       INTENT(INOUT)  :: InputFileData                       !< All the data in the ElastoDyn input file


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
   CHARACTER(*), PARAMETER      :: RoutineName = 'ReadTowerFile'


   ErrStat  =  ErrID_None
   ErrMsg   =  ""


   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN

      ! Open the tower input file.

   CALL OpenFInpFile ( UnIn, TwrFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! Add a separator to the echo file if appropriate.
   IF ( UnEc > 0 )  WRITE (UnEc,'(//,A,/)')  'Tower input data from file "'//TRIM( TwrFile )//'":'


   !  -------------- FILE HEADER ---------------------------------------------------

   CALL ReadCom ( UnIn, TwrFile, 'unused tower file header line 1', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL ReadCom ( UnIn, TwrFile, 'unused tower file header line 2', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


   !  -------------- TOWER PARAMETERS ---------------------------------------------

   CALL ReadCom ( UnIn, TwrFile, 'heading for tower parameters', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! NTwInpSt - Number of tower input stations.

   CALL ReadVar ( UnIn, TwrFile, InputFileData%NTwInpSt, 'NTwInpSt', 'Number of tower input stations', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! Allocate the input arrays based on this NTwInpSt input
   CALL Alloc_TowerInputProperties( InputFileData, ErrStat, ErrMsg )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! TwrFADmp - Tower fore-aft structural damping ratios.

   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwrFADmp, SIZE(InputFileData%TwrFADmp), 'TwrFADmp', &
                                     'Tower fore-aft structural damping ratios (%)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! TwrSSDmp - Tower side-to-side structural damping ratios.

   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwrSSDmp, SIZE(InputFileData%TwrSSDmp), 'TwrSSDmp', &
                                     'Tower side-to-side structural damping ratios (%)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


   !  -------------- TOWER ADJUSTMENT FACTORS -------------------------------------


      ! Skip the comment line.
   CALL ReadCom ( UnIn, TwrFile, 'heading for tower adjustment factors', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! FAStTunr - Tower fore-aft modal stiffness tuners.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%FAStTunr, SIZE(InputFileData%FAStTunr), 'FAStTunr', &
                                     'Tower fore-aft modal stiffness tuners (-)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! SSStTunr - Tower side-to-side modal stiffness tuners.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%SSStTunr, SIZE(InputFileData%SSStTunr), 'SSStTunr', &
                                     'Tower side-to-side modal stiffness tuners (-)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! AdjTwMa - Factor to adjust tower mass density.

   CALL ReadVar ( UnIn, TwrFile, AdjTwMa, 'AdjTwMa', 'Factor to adjust tower mass density (-)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF



      ! AdjFASt - Factor to adjust tower fore-aft stiffness.

   CALL ReadVar ( UnIn, TwrFile, AdjFASt, 'AdjFASt', 'Factor to adjust tower fore-aft stiffness (-)', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF



      ! AdjSSSt - Factor to adjust tower side-to-side stiffness.

   CALL ReadVar ( UnIn, TwrFile, AdjSSSt, 'AdjSSSt', 'Factor to adjust tower side-to-side stiffness (-)', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


         ! Check the locally-defined adjustment factors: AdjTwMa, AdjFASt, AdjSSSt

   IF ( AdjTwMa <= 0.0_ReKi ) THEN
      CALL SetErrStat( ErrID_Warn, 'AdjTwMa must be greater than zero.', ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   END IF

   IF ( AdjFASt <= 0.0_ReKi ) THEN
      CALL SetErrStat( ErrID_Warn, 'AdjFASt must be greater than zero.', ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   END IF

   IF ( AdjSSSt <= 0.0_ReKi ) THEN
      CALL SetErrStat( ErrID_Warn, 'AdjSSSt must be greater than zero.', ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   END IF


   !  -------------- DISTRIBUTED TOWER PROPERTIES ---------------------------------

      ! Skip the comment lines.
   CALL ReadCom ( UnIn, TwrFile, 'heading for distributed tower parameters', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL ReadCom ( UnIn, TwrFile, 'distributed-tower-parameter names', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL ReadCom ( UnIn, TwrFile, 'distributed-tower-parameter units', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF



      ! Read the table.
      
   NInputCols = 4


   DO I=1,InputFileData%NTwInpSt

      CALL ReadAry( UnIn, TwrFile, TmpRAry, NInputCols, 'Line'//TRIM(Num2LStr(I)), 'Tower input station table', &
                    ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF

      InputFileData%HtFract( I) = TmpRAry(1)
      InputFileData%TMassDen(I) = TmpRAry(2)*AdjTwMa   ! Apply the correction factors to the elemental data.
      InputFileData%TwFAStif(I) = TmpRAry(3)*AdjFASt   ! Apply the correction factors to the elemental data.
      InputFileData%TwSSStif(I) = TmpRAry(4)*AdjSSSt   ! Apply the correction factors to the elemental data.

   END DO ! I


   !  -------------- TOWER FORE-AFT MODE SHAPES -----------------------------------


      ! Skip the comment line.
   CALL ReadCom ( UnIn, TwrFile, 'heading for tower fore-aft mode shapes', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! TwFAM1Sh - Tower fore-aft mode-1 shape coefficients.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwFAM1Sh, SIZE(InputFileData%TwFAM1Sh), 'TwFAM1Sh', &
                           'Tower fore-aft mode-1 shape coefficients (-)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! TwFAM2Sh - Tower fore-aft mode-2 shape coefficients.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwFAM2Sh, SIZE(InputFileData%TwFAM2Sh), 'TwFAM2Sh', &
                           'Tower fore-aft mode-2 shape coefficients  (-)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


   !  -------------- TOWER SIDE-TO-SIDE MODE SHAPES -------------------------------


      ! Skip the comment line.
   CALL ReadCom ( UnIn, TwrFile, 'heading for tower side-to-side mode shapes', ErrStat2, ErrMsg2, UnEc  )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! TwSSM1Sh - Tower side-to-side mode-1 shape coefficients.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwSSM1Sh, SIZE(InputFileData%TwSSM1Sh), 'TwSSM1Sh', &
                           'Tower side-to-side mode-1 shape coefficients (-)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! TwSSM2Sh - Tower side-to-side mode-2 shape coefficients.
   CALL ReadAryLines ( UnIn, TwrFile, InputFileData%TwSSM2Sh, SIZE(InputFileData%TwSSM2Sh), 'TwSSM2Sh', &
                           'Tower side-to-side mode-2 shape coefficients (-)', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! Close the tower file.
   CALL Cleanup()


   RETURN
CONTAINS
   SUBROUTINE Cleanup()
      CLOSE( UnIn )
   END SUBROUTINE Cleanup
END SUBROUTINE ReadTowerFile
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine reads in the primary ElastoDyn input file and places the values it reads in the InputFileData structure.
!!  It opens an echo file if requested and returns the (still-open) echo file to the calling routine.
!!  It also returns the names of the BldFile, FurlFile, and TrwFile for further reading of inputs.
SUBROUTINE ReadPrimaryFile( InputFile, InputFileData, BldFile, FurlFile, TwrFile, OutFileRoot, UnEc, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables
   INTEGER(IntKi),     INTENT(OUT)    :: UnEc                                !< I/O unit for echo file. If > 0, file is open for writing.
   INTEGER(IntKi),     INTENT(OUT)    :: ErrStat                             !< Error status

   CHARACTER(*),       INTENT(IN)     :: InputFile                           !< Name of the file containing the primary input data
   CHARACTER(*),       INTENT(OUT)    :: ErrMsg                              !< Error message
   CHARACTER(*),       INTENT(OUT)    :: TwrFile                             !< name of the file containing tower inputs
   CHARACTER(*),       INTENT(OUT)    :: FurlFile                            !< name of the file containing furling inputs
   CHARACTER(*),       INTENT(OUT)    :: BldFile(MaxBl)                      !< name of the files containing blade inputs
   CHARACTER(*),       INTENT(IN)     :: OutFileRoot                         !< The rootname of the echo file, possibly opened in this routine

   TYPE(ED_InputFile), INTENT(INOUT)  :: InputFileData                       !< All the data in the ElastoDyn input file

      ! Local variables:
   INTEGER(IntKi)               :: I                                         ! loop counter
!   INTEGER(IntKi)               :: NumOuts                                  ! Number of output channel names read from the file
   INTEGER(IntKi)               :: UnIn                                      ! Unit number for reading file
   INTEGER(IntKi)               :: IOS
   INTEGER(IntKi)               :: ErrStat2                                  ! Temporary Error status
   LOGICAL                      :: Echo                                      ! Determines if an echo file should be written
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                   ! Temporary Error message
   CHARACTER(ErrMsgLen)         :: ErrMsg_NoAllBldNdOuts                     ! Temporary Error message
   CHARACTER(*), PARAMETER      :: RoutineName = 'ReadPrimaryFile'
   CHARACTER(1024)              :: PriPath                                   ! Path name of the primary file
   CHARACTER(1024)              :: FTitle                                    ! "File Title": the 2nd line of the input file, which contains a description of its contents
   CHARACTER(200)               :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")

      ! Initialize some variables:
   ErrStat  =  ErrID_None
   ErrMsg   =  ""
   Echo = .FALSE.
   UnEc = -1                             ! Echo file not opened, yet
   CALL GetPath( InputFile, PriPath )    ! Input files will be relative to the path where the primary input file is located.


      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN


      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


      ! Allocate arrays for input, based on maximum allowed number of blades and outputs
   CALL AllocAry( InputFileData%BlPitch, MaxBl, 'BlPitch input array', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL AllocAry( InputFileData%PreCone, MaxBl, 'Precone input array', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL AllocAry( InputFileData%TipMass, MaxBl, 'TipMass input array', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   CALL AllocAry( InputFileData%OutList, MaxOutPts, "ElastoDyn Input File's Outlist", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Allocate array for holding the list of node outputs
   CALL AllocAry( InputFileData%BldNd_OutList, BldNd_MaxOutPts, "BldNd_Outlist", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file.
   ! If Echo is TRUE, rewind and write on the second try.

   I    = 1 ! the number of times we've read the file (used for the Echo variable)
   DO
   !-------------------------- HEADER ---------------------------------------------
      CALL ReadCom( UnIn, InputFile, 'File Header: Module Version (line 1)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF

      CALL ReadStr( UnIn, InputFile, FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF

   !---------------------- SIMULATION CONTROL --------------------------------------
      CALL ReadCom( UnIn, InputFile, 'Section Header: Simulation Control', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF

         ! Echo - Echo input to "<RootName>.ech".

      CALL ReadVar( UnIn, InputFile, Echo, 'Echo',   'Echo switch', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF


      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop

         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read

      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)

      CALL OpenEcho ( UnEc, TRIM(OutFileRoot)//'.ech', ErrStat2, ErrMsg2, ED_Ver )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF

      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(ED_Ver%Name)//' primary input file "'//TRIM( InputFile )//'":'

      REWIND( UnIn, IOSTAT=ErrStat2 )
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".',ErrStat,ErrMsg,RoutineName )
            IF ( ErrStat >= AbortErrLev ) RETURN
         END IF

   END DO

   IF (NWTC_VerboseLevel == NWTC_Verbose) THEN
      CALL WrScr( ' Heading of the '//TRIM(ED_Ver%Name)//' input file: ' )
      CALL WrScr( '   '//TRIM( FTitle ) )
   END IF

      ! Method - Integration method for loose coupling
   CALL ReadVar( UnIn, InputFile, InputFileData%method, "Method", "Requested integration method for ElastoDyn {1: RK4, 2: AB4, or 3: ABM4}", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! DT - Requested integration time for ElastoDyn (seconds):
   CALL ReadVar( UnIn, InputFile, Line, "DT", "Requested integration time for ElastoDyn (seconds)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "DEFAULT" ) /= 1 ) THEN ! If it's not "default", read this variable; otherwise use the value already stored in InputFileData%DT
         READ( Line, *, IOSTAT=IOS) InputFileData%DT
         IF ( IOS /= 0 ) THEN
            CALL CheckIOS ( IOS, InputFile, "DT", NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL Cleanup()
               RETURN
            END IF
            RETURN
         END IF
      END IF

   !---------------------- DEGREES OF FREEDOM --------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Feature Flags', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! FlapDOF1 - First flapwise blade mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%FlapDOF1, "FlapDOF1", "First flapwise blade mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! FlapDOF2 - Second flapwise blade mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%FlapDOF2, "FlapDOF2", "Second flapwise blade mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! EdgeDOF - Edgewise blade mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%EdgeDOF, "EdgeDOF", "Edgewise blade mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TeetDOF - Rotor-teeter DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetDOF, "TeetDOF", "Rotor-teeter DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! DrTrDOF - Drivetrain rotational-flexibility DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%DrTrDOF, "DrTrDOF", "Drivetrain rotational-flexibility DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! GenDOF - Generator DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%GenDOF, "GenDOF", "Generator DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! YawDOF - Nacelle-yaw DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%YawDOF, "YawDOF", "Nacelle-yaw DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TwFADOF1 - First tower fore-aft bending-mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwFADOF1, "TwFADOF1", "First tower fore-aft bending-mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TwFADOF2 - Second tower fore-aft bending-mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwFADOF2, "TwFADOF2", "Second tower fore-aft bending-mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TwSSDOF1 - First tower side-to-side bending-mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwSSDOF1, "TwSSDOF1", "First tower side-to-side bending-mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TwSSDOF2 - Second tower side-to-side bending-mode DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwSSDOF2, "TwSSDOF2", "Second tower side-to-side bending-mode DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmSgDOF - Platform horizontal surge translation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmSgDOF, "PtfmSgDOF", "Platform horizontal surge translation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmSwDOF - Platform horizontal sway translation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmSwDOF, "PtfmSwDOF", "Platform horizontal sway translation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmHvDOF - Platform vertical heave translation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmHvDOF, "PtfmHvDOF", "Platform vertical heave translation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmRDOF - Platform roll tilt rotation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmRDOF, "PtfmRDOF", "Platform roll tilt rotation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmPDOF - Platform pitch tilt rotation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmPDOF, "PtfmPDOF", "Platform pitch tilt rotation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmYDOF - Platform yaw rotation DOF (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmYDOF, "PtfmYDOF", "Platform yaw rotation DOF (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   !---------------------- INITIAL CONDITIONS --------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Initial Conditions', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! OoPDefl - Initial out-of-plane blade-tip displacement (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%OoPDefl, "OoPDefl", "Initial out-of-plane blade-tip displacement (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! IPDefl - Initial in-plane blade-tip deflection (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%IPDefl, "IPDefl", "Initial in-plane blade-tip deflection (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! BlPitch - Initial blade pitch angles (deg) (read from file in degrees and converted to radians here):
   CALL ReadAryLines( UnIn, InputFile, InputFileData%BlPitch, MaxBl, "BlPitch", "Initial blade pitch angles (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%BlPitch = InputFileData%BlPitch*D2R

      ! TeetDefl - Initial teeter angle (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetDefl, "TeetDefl", "Initial teeter angle (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%TeetDefl = InputFileData%TeetDefl*D2R

      ! Azimuth - Initial azimuth angle for blade 1 (degrees) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%Azimuth, "Azimuth", "Initial azimuth angle for blade 1 (degrees)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%Azimuth = InputFileData%Azimuth*D2R

      ! RotSpeed - Initial rotor speed (RPM) (read in RPM and converted to rad/sec here):
   CALL ReadVar( UnIn, InputFile, InputFileData%RotSpeed, "RotSpeed", "Initial rotor speed (RPM)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%RotSpeed = InputFileData%RotSpeed*RPM2RPS

      ! NacYaw - Initial nacelle-yaw angle (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacYaw, "NacYaw", "Initial nacelle-yaw angle (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%NacYaw = InputFileData%NacYaw*D2R

      ! TTDspFA - Initial fore-aft tower-top displacement (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TTDspFA, "TTDspFA", "Initial fore-aft tower-top displacement (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TTDspSS - Initial side-to-side tower-top displacement (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TTDspSS, "TTDspSS", "Initial side-to-side tower-top displacement (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmSurge - Initial horizontal surge translational displacement of platform (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmSurge, "PtfmSurge", "Initial horizontal surge translational displacement of platform (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmSway - Initial horizontal sway translational displacement of platform (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmSway, "PtfmSway", "Initial horizontal sway translational displacement of platform (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmHeave - Initial vertical heave translational displacement of platform (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmHeave, "PtfmHeave", "Initial vertical heave translational displacement of platform (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmRoll - Initial roll tilt rotational displacement of platform (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmRoll, "PtfmRoll", "Initial roll tilt rotational displacement of platform (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%PtfmRoll = InputFileData%PtfmRoll*D2R

      ! PtfmPitch - Initial pitch tilt rotational displacement of platform (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmPitch, "PtfmPitch", "Initial pitch tilt rotational displacement of platform (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%PtfmPitch = InputFileData%PtfmPitch*D2R

      ! PtfmYaw - Initial yaw rotational displacement of platform (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmYaw, "PtfmYaw", "Initial yaw rotational displacement of platform (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%PtfmYaw = InputFileData%PtfmYaw*D2R

   !---------------------- TURBINE CONFIGURATION -----------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Turbine Configuration', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NumBl - Number of blades (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NumBl, "NumBl", "Number of blades (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TipRad - Preconed blade-tip radius (distance from the rotor apex to the blade tip) (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TipRad, "TipRad", "Preconed blade-tip radius (distance from the rotor apex to the blade tip) (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! HubRad - Preconed hub radius (distance from the rotor apex to the blade root) (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubRad, "HubRad", "Preconed hub radius (distance from the rotor apex to the blade root) (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PreCone - Rotor precone angles (deg) (read from file in degrees and converted to radians here):
   CALL ReadAryLines( UnIn, InputFile, InputFileData%PreCone, MaxBl, "PreCone", "Rotor precone angles (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%PreCone = InputFileData%PreCone*D2R

      ! HubCM - Distance from rotor apex to hub mass (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubCM, "HubCM", "Distance from rotor apex to hub mass (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! UndSling - Undersling length (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%UndSling, "UndSling", "Undersling length (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Delta3 - Delta-3 angle for teetering rotors (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%Delta3, "Delta3", "Delta-3 angle for teetering rotors (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%Delta3 = InputFileData%Delta3*D2R

      ! AzimB1Up - Azimuth value to use for I/O when blade 1 points up (degrees) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%AzimB1Up, "AzimB1Up", "Azimuth value to use for I/O when blade 1 points up (degrees)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%AzimB1Up = InputFileData%AzimB1Up*D2R

      ! OverHang - Distance from yaw axis to rotor apex or teeter pin (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%OverHang, "OverHang", "Distance from yaw axis to rotor apex or teeter pin (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! ShftGagL - Distance from hub or teeter pin to shaft strain gages (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%ShftGagL, "ShftGagL", "Distance from hub or teeter pin to shaft strain gages (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! ShftTilt - Rotor shaft tilt angle (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%ShftTilt, "ShftTilt", "Rotor shaft tilt angle (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%ShftTilt = InputFileData%ShftTilt*D2R

      ! NacCMxn - Downwind distance from tower-top to nacelle CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacCMxn, "NacCMxn", "Downwind distance from tower-top to nacelle CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NacCMyn - Lateral distance from tower-top to nacelle CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacCMyn, "NacCMyn", "Lateral distance from tower-top to nacelle CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NacCMzn - Vertical distance from tower-top to nacelle CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacCMzn, "NacCMzn", "Vertical distance from tower-top to nacelle CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NcIMUxn - Downwind distance from the tower-top to the nacelle IMU (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NcIMUxn, "NcIMUxn", "Downwind distance from the tower-top to the nacelle IMU (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NcIMUyn - Lateral distance from the tower-top to the nacelle IMU (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NcIMUyn, "NcIMUyn", "Lateral distance from the tower-top to the nacelle IMU (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NcIMUzn - Vertical distance from the tower-top to the nacelle IMU (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%NcIMUzn, "NcIMUzn", "Vertical distance from the tower-top to the nacelle IMU (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! Twr2Shft - Vertical distance from the tower-top to the rotor shaft (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%Twr2Shft, "Twr2Shft", "Vertical distance from the tower-top to the rotor shaft (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TowerHt - Height of tower relative to ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TowerHt, "TowerHt", "Height of tower above ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TowerBsHt - Height of tower base relative to ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%TowerBsHt, "TowerBsHt", "Height of tower base above ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! PtfmCMxt - Downwind distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmCMxt, "PtfmCMxt", "Downwind distance from the ground [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! PtfmCMyt - Lateral distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmCMyt, "PtfmCMyt", "Lateral distance from the ground [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      ! PtfmCMzt - Vertical distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmCMzt, "PtfmCMzt", "Vertical distance from the ground [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmRefzt - Vertical distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform reference point (meters):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmRefzt, "PtfmRefzt", "Vertical distance from the ground [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform reference point (meters)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   !---------------------- MASS AND INERTIA ----------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Mass and Inertia', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TipMass - Tip-brake masses (kg):
   CALL ReadAryLines( UnIn, InputFile, InputFileData%TipMass, MaxBl, "TipMass", "Tip-brake masses (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! HubMass - Hub mass (kg):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubMass, "HubMass", "Hub mass (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! HubIner - Hub inertia about teeter axis (2-blader) or rotor axis (3-blader) (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubIner, "HubIner", "Hub inertia about teeter axis (2-blader) or rotor axis (3-blader) (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! GenIner - Generator inertia about HSS (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%GenIner, "GenIner", "Generator inertia about HSS (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NacMass - Nacelle mass (kg):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacMass, "NacMass", "Nacelle mass (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NacYIner - Nacelle yaw inertia (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%NacYIner, "NacYIner", "Nacelle yaw inertia (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! YawBrMass - Yaw bearing mass (kg):
   CALL ReadVar( UnIn, InputFile, InputFileData%YawBrMass, "YawBrMass", "Yaw bearing mass (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmMass - Platform mass (kg):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmMass, "PtfmMass", "Platform mass (kg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmRIner - Platform inertia for roll tilt rotation about the platform CM (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmRIner, "PtfmRIner", "Platform inertia for roll tilt rotation about the platform CM (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmPIner - Platform inertia for pitch tilt rotation about the platform CM (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmPIner, "PtfmPIner", "Platform inertia for pitch tilt rotation about the platform CM (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! PtfmYIner - Platform inertia for yaw rotation about the platform CM (kg m^2):
   CALL ReadVar( UnIn, InputFile, InputFileData%PtfmYIner, "PtfmYIner", "Platform inertia for yaw rotation about the platform CM (kg m^2)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   !---------------------- BLADE ---------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Blade', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! BldNodes - Number of blade nodes (per blade) used for analysis (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%BldNodes, "BldNodes", "Number of blade nodes (per blade) used for analysis (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      
      
      ! InpBl - Blade file Input data for individual blades (see BladeInputData type):
   DO I = 1,MaxBl
      CALL ReadVar ( UnIn, InputFile, BldFile(I), 'BldFile('//TRIM(Num2Lstr(I))//')', 'Name of the file containing properties for blade '//TRIM(Num2Lstr(I)), ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
      IF ( PathIsRelative( BldFile(I) ) ) BldFile(I) = TRIM(PriPath)//TRIM(BldFile(I))
   END DO

   !---------------------- ROTOR-TEETER --------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Rotor-Teeter', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TeetMod - Rotor-teeter spring/damper model switch (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetMod, "TeetMod", "Rotor-teeter spring/damper model switch (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TeetDmpP - Rotor-teeter damper position (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetDmpP, "TeetDmpP", "Rotor-teeter damper position (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%TeetDmpP = InputFileData%TeetDmpP*D2R

      ! TeetDmp - Rotor-teeter damping constant (N-m/(rad/s)):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetDmp, "TeetDmp", "Rotor-teeter damping constant (N-m/(rad/s))", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TeetCDmp - Rotor-teeter rate-independent Coulomb-damping (N-m):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetCDmp, "TeetCDmp", "Rotor-teeter rate-independent Coulomb-damping (N-m)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TeetSStP - Rotor-teeter soft-stop position (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetSStP, "TeetSStP", "Rotor-teeter soft-stop position (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%TeetSStP = InputFileData%TeetSStP*D2R

      ! TeetHStP - Rotor-teeter hard-stop position (deg) (read from file in degrees and converted to radians here):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetHStP, "TeetHStP", "Rotor-teeter hard-stop position (deg)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%TeetHStP = InputFileData%TeetHStP*D2R

      ! TeetSSSp - Rotor-teeter soft-stop linear-spring constant (N-m/rad):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetSSSp, "TeetSSSp", "Rotor-teeter soft-stop linear-spring constant (N-m/rad)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TeetHSSp - Rotor-teeter hard-stop linear-spring constant (N-m/rad):
   CALL ReadVar( UnIn, InputFile, InputFileData%TeetHSSp, "TeetHSSp", "Rotor-teeter hard-stop linear-spring constant (N-m/rad)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   !---------------------- DRIVETRAIN ----------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Drivetrain', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! GBoxEff - Gearbox efficiency (%) (read from file in % and converted to fraction here):
   CALL ReadVar( UnIn, InputFile, InputFileData%GBoxEff, "GBoxEff", "Gearbox efficiency (%)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   InputFileData%GBoxEff = InputFileData%GBoxEff*0.01_ReKi

      ! GBRatio - Gearbox ratio (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%GBRatio, "GBRatio", "Gearbox ratio (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! DTTorSpr - Drivetrain torsional spring (N-m/rad):
   CALL ReadVar( UnIn, InputFile, InputFileData%DTTorSpr, "DTTorSpr", "Drivetrain torsional spring (N-m/rad)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! DTTorDmp - Drivetrain torsional damper (N-m/(rad/s)):
   CALL ReadVar( UnIn, InputFile, InputFileData%DTTorDmp, "DTTorDmp", "Drivetrain torsional damper (N-m/(rad/s))", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   !---------------------- FURLING -------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Furling', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Furling - Use Additional Furling parameters? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%Furling, "Furling", "Use Additional Furling parameters? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! FurlFile - Name of the file containing furling parameters:
   CALL ReadVar ( UnIn, InputFile, FurlFile, 'FurlFile', 'Name of the file containing furling parameters', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   IF ( PathIsRelative( FurlFile ) ) FurlFile = TRIM(PriPath)//TRIM(FurlFile)

   !---------------------- TOWER ---------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Tower', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TwrNodes - Number of tower nodes used in the analysis (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrNodes, "TwrNodes", "Number of tower nodes used in the analysis (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! TwrFile - Name of the file containing tower properties:
   CALL ReadVar ( UnIn, InputFile, TwrFile, 'TwrFile', 'Name of the file containing tower properties', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
   IF ( PathIsRelative( TwrFile ) ) TwrFile = TRIM(PriPath)//TRIM(TwrFile)

   !---------------------- OUTPUT --------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Output', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! SumPrint - Print summary data to <RootName>.sum (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%SumPrint, "SumPrint", "Print summary data to <RootName>.sum (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! OutFile - Switch to determine where output will be placed: (1: in module output file only; 2: in glue code output file only; 3: both) (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%OutFile, "OutFile", "Switch to determine where output will be placed: (1: in module output file only; 2: in glue code output file only; 3: both) (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   !    OutFileFmt - Format for module tabular (time-marching) output: (1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both):
   !CALL ReadVar( UnIn, InputFile, InputFileData%OutFileFmt, "OutFileFmt", "Format for module tabular (time-marching) output: (1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both)", ErrStat2, ErrMsg2, UnEc)
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF ( ErrStat >= AbortErrLev ) RETURN

      ! TabDelim - Flag to cause tab-delimited text output (delimited by space otherwise) (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TabDelim, "TabDelim", "Flag to cause tab-delimited text output (delimited by space otherwise) (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! OutFmt - Format used for module's text tabular output (except time); resulting field should be 10 characters (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%OutFmt, "OutFmt", "Format used for module's text tabular output (except time); resulting field should be 10 characters (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! Tstart - Time to start module's tabular output (seconds):
   CALL ReadVar( UnIn, InputFile, InputFileData%Tstart, "Tstart", "Time to start module's tabular output (seconds)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! DecFact - Decimation factor for module's tabular output (1=output every step) (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%DecFact, "DecFact", "Decimation factor for module's tabular output (1=output every step) (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NTwGages - Number of tower strain gages (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NTwGages, "NTwGages", "Number of tower strain gages (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      IF ( InputFileData%NTwGages > SIZE(InputFileData%TwrGagNd) ) THEN
         CALL SetErrStat( ErrID_Warn, 'Number of tower strain gages exceeds '// &
                                      TRIM(Num2LStr(SIZE(InputFileData%TwrGagNd)))//'.', ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
         InputFileData%NTwGages = SIZE(InputFileData%TwrGagNd)
      END IF

      ! TwrGagNd - Nodes closest to the tower strain gages (-):
   CALL ReadAry( UnIn, InputFile, InputFileData%TwrGagNd, InputFileData%NTwGages, "TwrGagNd", "Nodes closest to the tower strain gages (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! NBlGages - Number of blade strain gages (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NBlGages, "NBlGages", "Number of blade strain gages (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      IF ( InputFileData%NBlGages > SIZE(InputFileData%BldGagNd) ) THEN
         CALL SetErrStat( ErrID_Warn, 'Number of blade strain gages exceeds '//&
                        TRIM(Num2LStr(SIZE(InputFileData%BldGagNd))) //'.', ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF         
         InputFileData%NBlGages = SIZE(InputFileData%BldGagNd)
      END IF

      ! BldGagNd - Nodes closest to the blade strain gages (-):
   CALL ReadAry( UnIn, InputFile, InputFileData%BldGagNd, InputFileData%NBlGages, "BldGagNd", "Nodes closest to the blade strain gages (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

   !---------------------- OUTLIST  --------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: OutList', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF

      ! OutList - List of user-requested output channels (-):
   CALL ReadOutputList ( UnIn, InputFile, InputFileData%OutList, InputFileData%NumOuts, 'OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc  )     ! Routine in NWTC Subroutine Library
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF


   !----------- OUTLIST  -----------------------------------------------------------
      ! In case there is something ill-formed in the additional nodal outputs section, we will simply ignore it and assume that this section does not exist.
   ErrMsg_NoAllBldNdOuts='Nodal outputs section of ElastoDyn input file not found or improperly formatted.'

   !----------- OUTLIST for BldNd -----------------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: OutList for Blade node channels', ErrStat2, ErrMsg2, UnEc )
   IF ( ErrStat2 >= AbortErrLev ) THEN
      InputFileData%BldNd_BladesOut = 0
      InputFileData%BldNd_NumOuts = 0
      call wrscr( trim(ErrMsg_NoAllBldNdOuts) )
      CALL Cleanup()
      RETURN
   ENDIF


     ! Number of blade nodes to output:  will modify this at some point for arrays
      ! TODO:  In a future release, allow this to be an array of N blade numbers (change BldNd_BladesOut to an array if we do that).
      !        Will likely require reading this line in as a string (BldNd_BladesOut_Str) and parsing it
   CALL ReadVar(  UnIn, InputFile, InputFileData%BldNd_BladesOut, 'BldNd_BladesOut', 'Which blades to output node data on.'//TRIM(Num2Lstr(I)), ErrStat2, ErrMsg2, UnEc )
   IF ( ErrStat2 >= AbortErrLev ) THEN
      InputFileData%BldNd_BladesOut = 0
      InputFileData%BldNd_NumOuts = 0
      call wrscr( trim(ErrMsg_NoAllBldNdOuts) )
      CALL Cleanup()
      RETURN
   ENDIF


      ! Which blades to output for:  will add this at some point
      ! TODO: Parse this string into an array of nodes to output at (one idea is to set an array of boolean to T/F for which nodes to output).  At present, we ignore it entirely.
   CALL ReadVar(  UnIn, InputFile, InputFileData%BldNd_BlOutNd_Str, 'BldNd_BlOutNd_Str', 'Which nodes to output node data on.'//TRIM(Num2Lstr(I)), ErrStat2, ErrMsg2, UnEc )
   IF ( ErrStat2 >= AbortErrLev ) THEN
      InputFileData%BldNd_BladesOut = 0
      InputFileData%BldNd_NumOuts = 0
      call wrscr( trim(ErrMsg_NoAllBldNdOuts) )
      CALL Cleanup()
      RETURN
   ENDIF


      ! Section header for outlist
   CALL ReadCom( UnIn, InputFile, 'Section Header: OutList', ErrStat2, ErrMsg2, UnEc )
   IF ( ErrStat2 >= AbortErrLev ) THEN
      InputFileData%BldNd_BladesOut = 0
      InputFileData%BldNd_NumOuts = 0
      call wrscr( trim(ErrMsg_NoAllBldNdOuts) )
      CALL Cleanup()
      RETURN
   ENDIF


      ! OutList - List of user-requested output channels at each node(-):
   CALL ReadOutputList ( UnIn, InputFile, InputFileData%BldNd_OutList, InputFileData%BldNd_NumOuts, 'BldNd_OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc  )     ! Routine in NWTC Subroutine Library
   IF ( ErrStat2 >= AbortErrLev ) THEN
      InputFileData%BldNd_BladesOut = 0
      InputFileData%BldNd_NumOuts = 0
      call wrscr( trim(ErrMsg_NoAllBldNdOuts) )
      CALL Cleanup()
      RETURN
   ENDIF
   !---------------------- END OF FILE -----------------------------------------
   
   call cleanup()
   RETURN


CONTAINS
   SUBROUTINE cleanup()
      CLOSE ( UnIn )
   END SUBROUTINE cleanup
END SUBROUTINE ReadPrimaryFile
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine allocates arrays for the blade mesh properties from the input file.
SUBROUTINE Alloc_BladeMeshInputProperties( BladeKInputFileMesh, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_BladeMeshInputData),   INTENT(INOUT)  :: BladeKInputFileMesh      !< Data for Blade K stored in the module's input file
   INTEGER(IntKi),                INTENT(OUT)    :: ErrStat                  !< Error status
   CHARACTER(*),                  INTENT(OUT)    :: ErrMsg                   !< Err msg


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
!> This routine allocates arrays for the blade properties from the input file.
SUBROUTINE Alloc_BladeInputProperties( BladeKInputFileData, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(BladeInputData),     INTENT(INOUT)  :: BladeKInputFileData      !< Data for Blade K stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                  !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                   !< Err message


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


      ! BJJ: note that these used to be allocated 2:PolyOrd  :

   CALL AllocAry  ( BladeKInputFileData%BldFl1Sh,  PolyOrd-1, 'BldFl1Sh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileData%BldFl2Sh,  PolyOrd-1, 'BldFl2Sh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN
   CALL AllocAry  ( BladeKInputFileData%BldEdgSh,  PolyOrd-1, 'BldEdgSh'  , ErrStat, ErrMsg )
   IF ( ErrStat /= ErrID_None ) RETURN


END SUBROUTINE Alloc_BladeInputProperties
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine allocates arrays for the tower properties from the input file.
SUBROUTINE Alloc_TowerInputProperties( InputFileData, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(ED_InputFile),       INTENT(INOUT)  :: InputFileData      !< All the data in the ElastoDyn input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat            !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg             !< Error message


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
!> This routine checks the blade file input data for errors.
SUBROUTINE ValidateBladeData ( BladeKInputFileData, ErrStat, ErrMsg )
!..................................................................................................................................
   TYPE(BladeInputData),     INTENT(IN)     :: BladeKInputFileData                 !< Data for Blade K stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              !< Error message

      ! local variables
   INTEGER                                  :: I                                   ! Loop counter
   INTEGER(IntKi)                           :: ErrStat2                            ! Error status
   CHARACTER(ErrMsgLen)                     :: ErrMsg2                             ! Temporary error message
   CHARACTER(*), PARAMETER                  :: RoutineName = 'ValidateBladeData'


   ErrStat = ErrID_None
   ErrMsg= ''


      ! Check that BlFract goes from 0.0 to 1.0 in increasing order:

   IF ( .NOT. EqualRealNos( BladeKInputFileData%BlFract(1), 0.0_ReKi ) ) THEN
      CALL SetErrStat( ErrID_Fatal,'BlFract(1) must be 0.0.',ErrStat,ErrMsg,RoutineName)
   END IF

   IF ( BladeKInputFileData%NBlInpSt /= 1_IntKi .AND. &
      .NOT. EqualRealNos( BladeKInputFileData%BlFract(BladeKInputFileData%NBlInpSt), 1.0_ReKi )  ) THEN
      CALL SetErrStat( ErrID_Fatal,'BlFract('//TRIM( Num2LStr( BladeKInputFileData%NBlInpSt ) )//') must be 1.0.',ErrStat,ErrMsg,RoutineName)
   END IF

   DO I = 2,BladeKInputFileData%NBlInpSt
      IF ( BladeKInputFileData%BlFract(I) <= BladeKInputFileData%BlFract(I-1) )  THEN
         CALL SetErrStat( ErrID_Fatal,'BlFract('//TRIM( Num2LStr( I ) )//') must be greater than BlFract('&
                                                      //TRIM( Num2LStr(I-1) )//').',ErrStat,ErrMsg,RoutineName)
      ENDIF
   END DO


   DO I = 1,BladeKInputFileData%NBlInpSt

         ! Check that PitchAx is contained in [0.0, 1.0]:
      IF ( ( BladeKInputFileData%PitchAx(I) ) < 0.0_ReKi .OR. ( BladeKInputFileData%PitchAx(I) > 1.0_ReKi ) )  THEN
         CALL SetErrStat( ErrID_Fatal,'PitchAx('//TRIM( Num2LStr( I ) )//') must be between 0 and 1 (inclusive).',ErrStat,ErrMsg,RoutineName)
      END IF

         ! Check that StrcTwst is contained in (-pi,pi] radians ( i.e., (-180.0, 180.0] degrees):
      IF ( ( BladeKInputFileData%StrcTwst(I) <= -pi ) .OR. ( BladeKInputFileData%StrcTwst(I) > pi ) )  THEN
         CALL SetErrStat( ErrID_Fatal,'StrcTwst('//TRIM( Num2LStr( I ) ) // &
                     ') must be greater than -180 and less than or equal to 180.',ErrStat,ErrMsg,RoutineName)
      END IF

         ! Check that BMassDen is contained in (0.0, inf):
      IF ( BladeKInputFileData%BMassDen(I) <= 0.0_ReKi )  THEN
         CALL SetErrStat( ErrID_Fatal,'BMassDen('//TRIM( Num2LStr( I ) )//') must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      END IF

         ! Check that FlpStff is contained in (0.0, inf):
      IF ( BladeKInputFileData%FlpStff (I) <= 0.0_ReKi )  THEN
         CALL SetErrStat( ErrID_Fatal,'FlpStff('//TRIM( Num2LStr( I ) )//') must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      END IF

         ! Check that EdgStff is contained in (0.0, inf):
      IF ( BladeKInputFileData%EdgStff (I) <= 0.0_ReKi )  THEN
         CALL SetErrStat( ErrID_Fatal,'EdgStff('//TRIM( Num2LStr( I ) )//') must be greater than zero.',ErrStat,ErrMsg,RoutineName)
      END IF

   END DO


      ! Check that the blade damping is not negative:

   IF ( ANY( BladeKInputFileData%BldFlDmp < 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal,'BldFlDmp must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( ANY( BladeKInputFileData%BldEdDmp < 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal,'BldEdDmp must not be negative.',ErrStat,ErrMsg,RoutineName)


      ! Check that the stiffness tuner isn't negative:

   IF ( ANY( BladeKInputFileData%FlStTunr <= 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal,'FlStTunr must be greater than zero.',ErrStat,ErrMsg,RoutineName)


      ! Check that the mode shape coefficients are valid:

   CALL ValidateModeShapeCoeffs(BladeKInputFileData%BldFl1Sh, 'blade flap mode 1', ErrStat2, ErrMsg2 )
   CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)

   CALL ValidateModeShapeCoeffs(BladeKInputFileData%BldFl2Sh, 'blade flap mode 2', ErrStat2, ErrMsg2 )
   CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)

   CALL ValidateModeShapeCoeffs(BladeKInputFileData%BldEdgSh, 'blade edge', ErrStat2, ErrMsg2 )
   CALL SetErrStat( ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName)

END SUBROUTINE ValidateBladeData
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine checks that the mode shape coefficients add to 1.0, within numerical tolerance (currently 0.0015).
SUBROUTINE ValidateModeShapeCoeffs( Coeffs, ShpDesc, ErrStat, ErrMsg )
!..................................................................................................................................
   REAL(ReKi),               INTENT(IN )    :: Coeffs(:)                           !< Mode shape coefficients
   CHARACTER(*),             INTENT(IN)     :: ShpDesc                             !< Description of the mode shape for the error message
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              !< Error message

      ! local variables
   REAL(ReKi)                               :: Displ                               ! Blade tip/tower top displacement for a mode shape


      ! Check that the mode shape coefficients add to 1.0:

   Displ = SUM( Coeffs )
! bjj this new check seems to be a bit too restrictive for the input data:
!   IF ( .NOT. EqualRealNos( Displ, 1.0_ReKi ) ) THEN
   IF ( ABS( Displ - 1.0_ReKi ) > 0.0015_ReKi ) THEN
      ErrStat = ErrID_Fatal
      ErrMsg  = ' Mode shape coefficients for '//TRIM(ShpDesc)//' must add to 1.0.'
   ELSE
      ErrStat = ErrID_None
      ErrMsg  = ''
   END IF


END SUBROUTINE ValidateModeShapeCoeffs
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine checks the tower file input data for errors
SUBROUTINE ValidateTowerData ( InputFileData, ErrStat, ErrMsg )
!..................................................................................................................................
   TYPE(ED_InputFile),     INTENT(IN   )    :: InputFileData                       !< Data stored in the module's input file
   INTEGER(IntKi),         INTENT(  OUT)    :: ErrStat                             !< Error status
   CHARACTER(*),           INTENT(  OUT)    :: ErrMsg                              !< Error message

      ! local variables
   INTEGER                                  :: I                                   ! Loop counter
   INTEGER(IntKi)                           :: ErrStat2                            ! Error status
   CHARACTER(ErrMsgLen)                     :: ErrMsg2                             ! Temporary error message
   CHARACTER(*), PARAMETER                  :: RoutineName = 'ValidateTowerData'

   ErrStat = ErrID_None
   ErrMsg= ''



      ! Check that HtFract goes from 0.0 to 1.0 in increasing order:

   IF ( .NOT. EqualRealNos( InputFileData%HtFract(1), 0.0_ReKi ) ) CALL SetErrStat( ErrID_Fatal, 'HtFract(1) must be 0.0.', ErrStat, ErrMsg, RoutineName)

   IF ( InputFileData%NTwInpSt /= 1 .AND. &
      .NOT. EqualRealNos( InputFileData%HtFract(InputFileData%NTwInpSt), 1.0_ReKi )  ) THEN
      CALL SetErrStat( ErrID_Fatal, 'HtFract('//TRIM( Num2LStr( InputFileData%NTwInpSt ) )//') must be 1.0.', ErrStat, ErrMsg, RoutineName)
   END IF

   DO I = 2,InputFileData%NTwInpSt
      IF ( InputFileData%HtFract(I) <= InputFileData%HtFract(I-1) )  THEN
         CALL SetErrStat( ErrID_Fatal, 'HtFract('//TRIM( Num2LStr( I ) )//') must be greater than HtFract('&
                                                      //TRIM( Num2LStr(I-1) )//').', ErrStat, ErrMsg, RoutineName)

      ENDIF
   END DO


      ! Check the input arrays:

   DO I = 1,InputFileData%NTwInpSt
      IF ( InputFileData%TMassDen(I) <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Fatal, 'TMassDen('//TRIM(Num2LStr( I ))//') must be greater than zero.', ErrStat, ErrMsg, RoutineName)
      END IF

      IF ( InputFileData%TwFAStif(I) <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Fatal, 'TwFAStif('//TRIM(Num2LStr( I ))//') must be greater than zero.', ErrStat, ErrMsg, RoutineName)
      END IF

      IF ( InputFileData%TwSSStif(I) <= 0.0_ReKi ) THEN
         CALL SetErrStat( ErrID_Fatal, 'TwSSStif('//TRIM(Num2LStr( I ))//') must be greater than zero.', ErrStat, ErrMsg, RoutineName)
      END IF
   END DO


      ! Check that the tower damping (TwrFADmp) is contained in the range [0, 100]:

   IF ( ANY( InputFileData%TwrFADmp < 0.0_ReKi ) .OR. ANY( InputFileData%TwrFADmp > 100.0_ReKi ) ) THEN
      CALL SetErrStat( ErrID_Fatal, 'TwrFADmp must be between 0 and 100 (inclusive).', ErrStat, ErrMsg, RoutineName)
   END IF

   IF ( ANY( InputFileData%TwrSSDmp < 0.0_ReKi ) .OR. ANY( InputFileData%TwrSSDmp > 100.0_ReKi ) ) THEN
      CALL SetErrStat( ErrID_Fatal, 'TwrSSDmp must be between 0 and 100 (inclusive).', ErrStat, ErrMsg, RoutineName)
   END IF


      ! Check that the tower tuners are positive numbers:

   IF ( ANY( InputFileData%FAStTunr <= 0.0_ReKi )  ) CALL SetErrStat( ErrID_Fatal, 'FAStTunr must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   IF ( ANY( InputFileData%SSStTunr <= 0.0_ReKi )  ) CALL SetErrStat( ErrID_Fatal, 'SSStTunr must be greater than zero.', ErrStat, ErrMsg, RoutineName )



      ! Validate the mode shape coefficients:

   CALL ValidateModeShapeCoeffs( InputFileData%TwFAM1Sh, 'tower fore-aft mode 1', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   CALL ValidateModeShapeCoeffs( InputFileData%TwFAM2Sh, 'tower fore-aft mode 2', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   CALL ValidateModeShapeCoeffs( InputFileData%TwSSM1Sh, 'tower side-to-side mode 1', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


   CALL ValidateModeShapeCoeffs( InputFileData%TwSSM2Sh, 'tower side-to-side mode 2', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

END SUBROUTINE ValidateTowerData
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine validates the furling inputs.
SUBROUTINE ValidateFurlData( InputFileData, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables:

   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                       !< All the data in the ElastoDyn input file

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              !< Error message

      ! Local variables:
!   CHARACTER(1024)                          :: TmpMsg                              ! a temporary message (so I don't have to keep typing the same error message)
   REAL(ReKi)                               :: SmallAngleLimit_Rad                 ! Largest input angle considered "small" (check in input file), radians
   CHARACTER(*), PARAMETER                  :: RoutineName = 'ValidateFurlData'

      ! Initialize error status and angle limit defined locally (in correct units)

   ErrStat = ErrID_None
   ErrMsg  = ''
   SmallAngleLimit_Rad = SmallAngleLimit_Deg*D2R


      ! note that all angles are assumed to be in radians here:

      ! Check that angles are in the range (-pi, pi] radians (i.e., (-180, 180] degrees ):
      ! NOTE: these are local subroutines, with ErrStat and ErrMsg INTENT(INOUT)

   CALL CheckAngle180Range( InputFileData%RotFurl,  'RotFurl'  )
   CALL CheckAngle180Range( InputFileData%TailFurl, 'TailFurl' )
   CALL CheckAngle180Range( InputFileData%RFrlSkew, 'RFrlSkew' )
   CALL CheckAngle180Range( InputFileData%TFrlSkew, 'TFrlSkew' )

   CALL CheckAngle180Range( InputFileData%RFrlUSSP, 'RFrlUSSP' )
   CALL CheckAngle180Range( InputFileData%TFrlUSSP, 'TFrlUSSP' )
   CALL CheckAngle180Range( InputFileData%RFrlUSDP, 'RFrlUSDP' )
   CALL CheckAngle180Range( InputFileData%TFrlUSDP, 'TFrlUSDP' )

   CALL CheckAngle180Range( InputFileData%RFrlDSSP, 'RFrlDSSP' )
   IF ( InputFileData%RFrlDSSP > InputFileData%RFrlUSSP ) THEN
      CALL SetErrStat( ErrID_Fatal,'RFrlDSSP must not be larger than RFrlUSSP.',ErrStat,ErrMsg,RoutineName)
   END IF

   CALL CheckAngle180Range( InputFileData%RFrlDSDP, 'RFrlDSDP' )
   IF ( InputFileData%RFrlDSDP > InputFileData%RFrlUSDP ) THEN
      CALL SetErrStat( ErrID_Fatal,'RFrlDSDP must not be larger than RFrlUSDP.',ErrStat,ErrMsg,RoutineName)
   END IF

   CALL CheckAngle180Range( InputFileData%TFrlDSSP, 'TFrlDSSP' )
   IF ( InputFileData%TFrlDSSP > InputFileData%TFrlUSSP ) THEN
      CALL SetErrStat( ErrID_Fatal,'TFrlDSSP must not be larger than TFrlUSSP.',ErrStat,ErrMsg,RoutineName)
   END IF

   CALL CheckAngle180Range( InputFileData%TFrlDSDP, 'TFrlDSDP' )
   IF ( InputFileData%TFrlDSDP > InputFileData%TFrlUSDP ) THEN
      CALL SetErrStat( ErrID_Fatal,'TFrlDSDP must not be larger than TFrlUSDP.',ErrStat,ErrMsg,RoutineName)
   END IF


      ! Check that tilt angles are in the range [-pi/2, pi/2] radians (i.e., [-90, 90] degrees ):
   IF ( ABS( InputFileData%RFrlTilt ) > PiBy2 )  THEN
      CALL SetErrStat(ErrID_Fatal,'RFrlTilt must be between -pi/2 and pi/2 radians (i.e., in the range [-90, 90] degrees).',ErrStat,ErrMsg,RoutineName)
   END IF
   IF ( ABS( InputFileData%TFrlTilt ) > PiBy2 )  THEN
      CALL SetErrStat(ErrID_Fatal,'TFrlTilt must be between -pi/2 and pi/2 radians (i.e., in the range [-90, 90] degrees).',ErrStat,ErrMsg,RoutineName)
   END IF

   

      ! Check for violations of the small-angle assumption (15-degree limit, using radians):

   IF ( ABS( InputFileData%ShftSkew ) > SmallAngleLimit_Rad )  THEN
      CALL SetErrStat( ErrID_Fatal,'ShftSkew should only be used to skew the shaft a few degrees away from the zero-yaw ' &
                //'position and must not be used as a replacement for the yaw angle. '&
                //'ShftSkew must be between -'//TRIM(Num2LStr(SmallAngleLimit_Rad))//' and ' &
                                              //TRIM(Num2LStr(SmallAngleLimit_Rad))//' radians.',ErrStat,ErrMsg,RoutineName )
   END IF


      ! Warn if tail is defined upwind of the tower:

   IF ( InputFileData%BoomCM_n(1) < 0.0_ReKi )  THEN   ! Print out warning when tail boom CM defined upwind of the tower.
      CALL SetErrStat( ErrID_Warn,'WARNING: Tail boom CM is defined upwind of the tower (BoomCM_n(1) < 0).',ErrStat,ErrMsg,RoutineName)
   ENDIF

   IF ( InputFileData%TFinCM_n(1) < 0.0_ReKi )  THEN   ! Print out warning when tail fin CM defined upwind of the tower.
      CALL SetErrStat( ErrID_Warn,'WARNING: Tail fin CM is defined upwind of the tower (TFinCM_n(1) < 0).',ErrStat,ErrMsg,RoutineName)
   ENDIF

      ! Check that mass, inertias, damping, etc. values aren't negative:

   IF (InputFileData%RFrlMass   < 0.0_ReKi) call SetErrStat(ErrID_Fatal, 'RFrlMass must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (InputFileData%BoomMass   < 0.0_ReKi) call SetErrStat(ErrID_Fatal, 'BoomMass must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (InputFileData%TFinMass   < 0.0_ReKi) call SetErrStat(ErrID_Fatal, 'TFinMass must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (InputFileData%RFrlIner   < 0.0_ReKi) call SetErrStat(ErrID_Fatal, 'RFrlIner must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (InputFileData%TFrlIner   < 0.0_ReKi) call SetErrStat(ErrID_Fatal, 'TFrlIner must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF (InputFileData%RFrlSpr    < 0.0_ReKi) call SetErrStat(ErrID_Fatal,  'RFrlSpr must not be negative.',ErrStat,ErrMsg,RoutineName)   
   IF ( InputFileData%RFrlDmp   < 0.0_ReKi) call SetErrStat(ErrID_Fatal,  'RFrlDmp must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%RFrlUSSpr < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'RFrlUSSpr must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%RFrlDSSpr < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'RFrlDSSpr must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%RFrlUSDmp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'RFrlUSDmp must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%RFrlDSDmp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'RFrlDSDmp must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%TFrlSpr   < 0.0_ReKi) call SetErrStat(ErrID_Fatal,  'TFrlSpr must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%TFrlDmp   < 0.0_ReKi) call SetErrStat(ErrID_Fatal,  'TFrlDmp must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%TFrlUSSpr < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TFrlUSSpr must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%TFrlDSSpr < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TFrlDSSpr must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%TFrlUSDmp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TFrlUSDmp must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%TFrlDSDmp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TFrlDSDmp must not be negative.',ErrStat,ErrMsg,RoutineName)


      ! Check that furling models are valid:

   IF ( InputFileData%TFrlMod < 0 .OR. InputFileData%TFrlMod > 2 )  THEN
      CALL SetErrStat( ErrID_Fatal, 'TFrlMod must be 0, 1, or 2.',ErrStat,ErrMsg,RoutineName)
   END IF

   IF ( InputFileData%RFrlMod < 0 .OR. InputFileData%RFrlMod > 2 )  THEN
      CALL SetErrStat( ErrID_Fatal, 'RFrlMod must be 0, 1, or 2.',ErrStat,ErrMsg,RoutineName )
   END IF


   !   ! bjj: THESE ARE checks for tail fin aerodynamics, which should be in aerodyn, in my opinion
   !IF ( InputFileData%TFinArea < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TFinArea must not be negative.',ErrStat,ErrMsg,RoutineName)

   !IF ( TFinMod < 0 .OR. TFinMod > 2 )  THEN
   !   CALL SetErrors( ErrID_Fatal,'TFinMod must be 0, 1, or 2.')
   !END IF


   RETURN

CONTAINS
   !-------------------------------------------------------------------------------------------------------------------------------
   !> This routine checks that an angle is in the range (-pi, pi] radians. If not, ErrStat = ErrID_Fatal
   !! Note that all values are assumed to be in radians, even if read in degrees (-180 deg, 180 deg]
   SUBROUTINE CheckAngle180Range( Var, VarDesc )
   !...............................................................................................................................
   REAL(ReKi),     INTENT(IN)    :: Var         ! Variable to check
   CHARACTER(*),   INTENT(IN)    :: VarDesc     ! Description of variable (used in error message)

      IF ( ( Var <= -pi ) .OR. ( Var > pi ) )  THEN
         CALL SetErrStat(ErrID_Fatal, TRIM(VarDesc)// &
            ' must be greater than -pi radians and less than or equal to pi radians '// &
            '(i.e., in the range (-180, 180] degrees).',ErrStat,ErrMsg,RoutineName)
      END IF

   END SUBROUTINE CheckAngle180Range
   !-------------------------------------------------------------------------------------------------------------------------------
 END SUBROUTINE ValidateFurlData
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine validates the inputs from the primary input file.
!! note that all angles are assumed to be in radians in this routine:
SUBROUTINE ValidatePrimaryData( InputFileData, BD4Blades, Linearize, MHK, ErrStat, ErrMsg )
!..................................................................................................................................

      ! Passed variables:

   TYPE(ED_InputFile),       INTENT(IN)     :: InputFileData                       !< All the data in the ElastoDyn input file
   LOGICAL,                  INTENT(IN)     :: BD4Blades                           !< Use BeamDyn for blades, thus ignore ElastoDyn blade info 
   LOGICAL,                  INTENT(IN)     :: Linearize                           !< Flag indicating glue code wants to linearize this module
   INTEGER(IntKi),           INTENT(IN)     :: MHK                                 !< MHK turbine type switch
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             !< Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              !< Error message

      ! Local variables:
   REAL(ReKi)                               :: SmallAngleLimit_Rad                 ! Largest input angle considered "small" (check in input file), radians
   INTEGER(IntKi)                           :: I                                   ! loop counter
   INTEGER(IntKi)                           :: K                                   ! blade number
   INTEGER(IntKi)                           :: FmtWidth                            ! width of the field returned by the specified OutFmt
   INTEGER(IntKi)                           :: ErrStat2                            ! Temporary error status
   CHARACTER(ErrMsgLen)                     :: ErrMsg2                             ! Temporary rror message
   CHARACTER(*), PARAMETER                  :: RoutineName = 'ValidatePrimaryData'

      ! Initialize error status and angle limit defined locally (in correct units)

   ErrStat = ErrID_None
   ErrMsg  = ''
   SmallAngleLimit_Rad = SmallAngleLimit_Deg*D2R

      ! Make sure the number of blades is valid:
   IF (BD4Blades) THEN
      IF ( ( InputFileData%NumBl < 1 ) .OR. ( InputFileData%NumBl > MaxBl ) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'NumBl must be 1, 2, or 3 for BeamDyn simulations.',ErrStat,ErrMsg,RoutineName)
      END IF
   ELSE
      IF ( ( InputFileData%NumBl < 1 ) .OR. ( InputFileData%NumBl > MaxBl ) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'NumBl must be 1, 2, or 3.',ErrStat,ErrMsg,RoutineName)
      END IF
   END IF

      ! Make sure the specified integration method makes sense:
   IF ( InputFileData%method .ne. Method_RK4) THEN
      IF ( InputFileData%method .ne. Method_AB4) THEN
         IF ( InputFileData%method .ne. Method_ABM4) THEN
            CALL SetErrStat( ErrID_Fatal, 'Integration method must be 1 (RK4), 2 (AB4), or 3 (ABM4)',ErrStat,ErrMsg,RoutineName )
         END IF
      END IF
   END IF


      ! make sure GBoxEff is 100% for now
   IF ( .NOT. EqualRealNos( InputFileData%GBoxEff, 1.0_ReKi ) .and. InputFileData%method == method_rk4 ) &
      CALL SetErrStat( ErrID_Fatal, 'GBoxEff must be 1 (i.e., 100%) when using RK4.',ErrStat,ErrMsg,RoutineName)


      ! Don't allow these parameters to be negative (i.e., they must be in the range (0,inf)):
   IF ( InputFileData%RotSpeed < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'RotSpeed must not be negative.',ErrStat,ErrMsg,RoutineName)
   
   IF ( InputFileData%TipRad < 0.0_ReKi) CALL SetErrStat(ErrID_Fatal,'TipRad must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%HubRad < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'HubRad must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%DTTorSpr < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'DTTorSpr must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%DTTorDmp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'DTTorDmp must not be negative.',ErrStat,ErrMsg,RoutineName)

   IF ( InputFileData%PtfmMass  < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'PtfmMass must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%PtfmRIner < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'PtfmRIner must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%PtfmPIner < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'PtfmPIner must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%PtfmYIner < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'PtfmYIner must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%YawBrMass < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'YawBrMass must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%NacMass   < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'NacMass must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%HubMass   < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'HubMass must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( MHK /= 2 ) THEN
      IF ( InputFileData%Twr2Shft  < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'Twr2Shft must not be negative.',ErrStat,ErrMsg,RoutineName)
   ELSEIF ( MHK == 2 ) THEN
      IF ( InputFileData%Twr2Shft  > 0.0_ReKi) call SetErrStat(ErrID_Fatal,'Twr2Shft must not be positive for a floating MHK turbine.',ErrStat,ErrMsg,RoutineName)
   ENDIF
      
   DO K=1,InputFileData%NumBl
      IF ( InputFileData%TipMass(K) < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TipMass('//TRIM( Num2LStr( K ) )//') must not be negative.',ErrStat,ErrMsg,RoutineName)
   ENDDO ! K

   IF ( InputFileData%NacYIner < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'NacYIner must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%GenIner  < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'GenIner must not be negative.',ErrStat,ErrMsg,RoutineName)
   IF ( InputFileData%HubIner  < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'HubIner must not be negative.',ErrStat,ErrMsg,RoutineName)

      ! Check that TowerHt is in the range [0,inf):
   IF ( MHK /= 2 ) THEN
      IF ( InputFileData%TowerHt <= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'TowerHt must be greater than zero.',ErrStat,ErrMsg,RoutineName )
   ELSEIF ( MHK == 2 ) THEN
      IF ( InputFileData%TowerHt >= 0.0_ReKi ) CALL SetErrStat( ErrID_Fatal, 'TowerHt must be less than zero for a floating MHK turbine.',ErrStat,ErrMsg,RoutineName )
   ENDIF

      ! Check that these integers are in appropriate ranges:
   IF ( InputFileData%TwrNodes < 1_IntKi ) CALL SetErrStat( ErrID_Fatal, 'TwrNodes must not be less than 1.',ErrStat,ErrMsg,RoutineName )

   IF ( InputFileData%BldNodes < 1_IntKi .AND. .NOT. BD4Blades) CALL SetErrStat( ErrID_Fatal, 'BldNodes must not be less than 1.',ErrStat,ErrMsg,RoutineName )

      ! Check that the gearbox efficiency is valid:
   IF ( ( InputFileData%GBoxEff <= 0.0_ReKi ) .OR. ( InputFileData%GBoxEff > 1.0_ReKi ) ) THEN
         CALL SetErrStat( ErrID_Fatal, 'GBoxEff must be in the range (0,1] (i.e., (0,100] percent).',ErrStat,ErrMsg,RoutineName )
      ENDIF

      ! warn if 2nd modes are enabled without their corresponding 1st modes

   IF ( InputFileData%FlapDOF2 .AND. ( .NOT. InputFileData%FlapDOF1 ) )  THEN  ! Print out warning when flap mode 1 is not enabled and flap mode 2 is enabled
      CALL SetErrStat( ErrID_Warn, '2nd blade flap mode is enabled without the 1st. '//&
                    ' This designation is recommended only for debugging purposes.',ErrStat,ErrMsg,RoutineName )
   ENDIF

   IF ( InputFileData%TwFADOF2 .AND. ( .NOT. InputFileData%TwFADOF1 ) )  THEN  ! Print out warning when tower fore-aft mode 1 is not enabled and fore-aft mode 2 is enabled
      CALL SetErrStat( ErrID_Warn, '2nd tower fore-aft mode is enabled without the 1st. '//&
                    ' This designation is recommended only for debugging purposes.',ErrStat,ErrMsg,RoutineName )
   ENDIF

   IF ( InputFileData%TwSSDOF2 .AND. ( .NOT. InputFileData%TwSSDOF1 ) )  THEN  ! Print out warning when tower side-to-side mode 1 is not enabled and side-to-side mode 2 is enabled
      CALL SetErrStat( ErrID_Warn, '2nd tower side-to-side mode is enabled without the 1st. '//&
                    ' This designation is recommended only for debugging purposes.',ErrStat,ErrMsg,RoutineName )
   ENDIF


      ! Check that turbine configuration makes sense:

   IF ( InputFileData%Furling ) THEN
      IF ( InputFileData%OverHang > 0.0_ReKi )  THEN   ! Print out warning when downwind turbine is modeled with furling.
         CALL SetErrStat( ErrID_Warn, 'Furling designation (Furling = True) specified for downwind rotor '// &
                    'configuration (OverHang > 0). Check for possible errors in the input file(s).',ErrStat,ErrMsg,RoutineName)
      END IF
   ENDIF

   IF ( MHK /= 2 ) THEN

      IF ( InputFileData%TowerBsHt >= InputFileData%TowerHt ) CALL SetErrStat( ErrID_Fatal, 'TowerBsHt must be less than TowerHt.',ErrStat,ErrMsg,RoutineName)

      IF ( InputFileData%PtfmCMzt  > InputFileData%TowerBsHt ) &
         CALL SetErrStat( ErrID_Fatal, 'PtfmCMzt must not be greater than TowerBsHt.',ErrStat,ErrMsg,RoutineName)
   
      IF ( InputFileData%PtfmRefzt  > InputFileData%TowerBsHt ) &
         CALL SetErrStat( ErrID_Fatal, 'PtfmRefzt must not be greater than TowerBsHt.',ErrStat,ErrMsg,RoutineName)

   ELSEIF ( MHK == 2 ) THEN

      IF ( InputFileData%TowerBsHt <= InputFileData%TowerHt ) CALL SetErrStat( ErrID_Fatal, 'TowerBsHt must be greater than TowerHt for a floating MHK turbine.',ErrStat,ErrMsg,RoutineName)
         
   ENDIF
     
   IF (.NOT. BD4Blades ) THEN
      IF (InputFileData%HubRad >= InputFileData%TipRad ) &
      CALL SetErrStat( ErrID_Fatal, 'HubRad must be less than TipRad.',ErrStat,ErrMsg,RoutineName)

      IF ( MHK /= 2 ) THEN
         IF ( InputFileData%TowerHt + InputFileData%Twr2Shft + InputFileData%OverHang*SIN(InputFileData%ShftTilt) &
                                    <= InputFileData%TipRad )  THEN
            CALL SetErrStat( ErrID_Fatal, 'TowerHt + Twr2Shft + OverHang*SIN(ShftTilt) must be greater than TipRad.',ErrStat,ErrMsg,RoutineName)
         END IF
      ELSEIF ( MHK == 2 ) THEN
         IF ( -InputFileData%TowerHt - InputFileData%Twr2Shft - InputFileData%OverHang*SIN(InputFileData%ShftTilt) &
                                    <= InputFileData%TipRad )  THEN
            CALL SetErrStat( ErrID_Fatal, 'TowerHt + Twr2Shft + OverHang*SIN(ShftTilt) must be greater than TipRad.',ErrStat,ErrMsg,RoutineName)
         END IF
      ENDIF
   END IF


   IF ( InputFileData%NumBl == 2_IntKi )  THEN
      IF ( ( InputFileData%TeetDefl <= -pi ) .OR. ( InputFileData%TeetDefl > pi ) )  &
         CALL SetErrStat( ErrID_Fatal, 'TeetDefl must be in the range (-pi, pi] radians (i.e., [-180,180] degrees).',ErrStat,ErrMsg,RoutineName)

      IF ( ABS( InputFileData%Delta3 ) >= PiBy2 )  &
         CALL SetErrStat( ErrID_Fatal, 'Delta3 must be in the range (pi/2, pi/2) radians (i.e., (-90, 90) degrees).',ErrStat,ErrMsg,RoutineName)

      IF ( ( InputFileData%TeetSStP < 0.0_ReKi ) .OR. ( InputFileData%TeetSStP > pi ) )  &
         CALL SetErrStat( ErrID_Fatal, 'TeetSStP must be in the range [0, pi] radians (i.e., [0,180] degrees).',ErrStat,ErrMsg,RoutineName)

      IF ( ( InputFileData%TeetDmpP < 0.0_ReKi ) .OR. ( InputFileData%TeetDmpP > pi ) )  &
         CALL SetErrStat( ErrID_Fatal, 'TeetDmpP must be in the range [0, pi] radians (i.e., [0,180] degrees).',ErrStat,ErrMsg,RoutineName)

      IF ( ( InputFileData%TeetHStP < InputFileData%TeetSStP ) .OR. ( InputFileData%TeetHStP > pi ) )  &
         CALL SetErrStat( ErrID_Fatal, 'TeetHStP must be in the range [TeetSStP, pi] radians (i.e., [TeetSStP, 180] degrees).',ErrStat,ErrMsg,RoutineName)

      IF ( ( InputFileData%TeetMod /= 0_IntKi ) .AND. ( InputFileData%TeetMod /= 1_IntKi ) .AND. &
           ( InputFileData%TeetMod /= 2_IntKi ) )  &
         CALL SetErrStat( ErrID_Fatal, 'TeetMod must be 0, 1, or 2.',ErrStat,ErrMsg,RoutineName)

      IF ( InputFileData%TeetDmp  < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TeetDmp must not be negative.',ErrStat,ErrMsg,RoutineName)
      IF ( InputFileData%TeetCDmp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TeetCDmp must not be negative.',ErrStat,ErrMsg,RoutineName)
      IF ( InputFileData%TeetSSSp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TeetSSSp must not be negative.',ErrStat,ErrMsg,RoutineName)
      IF ( InputFileData%TeetHSSp < 0.0_ReKi) call SetErrStat(ErrID_Fatal,'TeetHSSp must not be negative.',ErrStat,ErrMsg,RoutineName)
            
      ! TeetCDmp isn't allowed to be non-zero in this verison.
      IF (.NOT. EqualRealNos( InputFileData%TeetCDmp, 0.0_ReKi )) &
         CALL SetErrStat( ErrID_Fatal, 'TeetCDmp must be 0 in this version of ElastoDyn.',ErrStat,ErrMsg,RoutineName)
   ENDIF

      ! check these angles for appropriate ranges:
   IF ( ( InputFileData%NacYaw <= -pi ) .OR. ( InputFileData%NacYaw > pi ) ) &
      CALL SetErrStat( ErrID_Fatal, 'NacYaw must be in the range (-pi, pi] radians (i.e., (-180, 180] degrees).',ErrStat,ErrMsg,RoutineName)
   IF ( ( InputFileData%Azimuth  < 0.0_ReKi ) .OR. ( InputFileData%Azimuth  >= TwoPi ) ) &
      CALL SetErrStat( ErrID_Fatal, 'Azimuth must be in the range [0, 2pi) radians (i.e., [0, 360) degrees).',ErrStat,ErrMsg,RoutineName)
   IF ( ( InputFileData%AzimB1Up < 0.0_ReKi ) .OR. ( InputFileData%AzimB1Up >= TwoPi ) )  &
      CALL SetErrStat( ErrID_Fatal, 'AzimB1Up must be in the range [0, 2pi) radians (i.e., [0, 360) degrees).',ErrStat,ErrMsg,RoutineName)

   DO K=1,InputFileData%NumBl
      IF ( ( InputFileData%BlPitch(K) <= -pi ) .OR. ( InputFileData%BlPitch(K) > pi ) )  THEN
         CALL SetErrStat( ErrID_Fatal, 'BlPitch('//TRIM(Num2LStr(K))//')'//' must be greater than -pi radians and '// &
                                      'less than or equal to pi radians (i.e., in the range (-180, 180] degrees).',ErrStat,ErrMsg,RoutineName)
      END IF
      IF ( ABS( InputFileData%PreCone(K) ) >= PiBy2 )  &
         CALL SetErrStat( ErrID_Fatal, 'PreCone('//TRIM( Num2LStr( K ) )//') must be in the range (-pi/2, pi/2) '//&
                                      'radians (i.e., (-90, 90) degrees).',ErrStat,ErrMsg,RoutineName)
   END DO

      ! Check that these angles are in the range [-pi/2, pi/2] radians (i.e., [-90, 90] degrees ):
   IF ( ABS( InputFileData%ShftTilt ) > PiBy2 )  THEN
      CALL SetErrStat(ErrID_Fatal,'ShftTilt must be between -pi/2 and pi/2 radians (i.e., in the range [-90, 90] degrees).',ErrStat,ErrMsg,RoutineName)
   END IF


      ! Check for violations of the small-angle assumption (15-degree limit, using radians):
   IF ( ABS( InputFileData%PtfmRoll ) > SmallAngleLimit_Rad ) THEN
      CALL SetErrStat( ErrID_Fatal, 'PtfmRoll must be between -'//TRIM(Num2LStr(SmallAngleLimit_Rad))//' and ' &
                                    //TRIM(Num2LStr(SmallAngleLimit_Rad))//' radians.',ErrStat,ErrMsg,RoutineName)
   END IF

   IF ( ABS( InputFileData%PtfmPitch ) > SmallAngleLimit_Rad ) THEN
      CALL SetErrStat( ErrID_Fatal, 'PtfmPitch must be between -'//TRIM(Num2LStr(SmallAngleLimit_Rad))//' and ' &
                                    //TRIM(Num2LStr(SmallAngleLimit_Rad))//' radians.',ErrStat,ErrMsg,RoutineName)
   END IF

   IF ( ABS( InputFileData%PtfmYaw ) > SmallAngleLimit_Rad ) THEN
      CALL SetErrStat( ErrID_Fatal, 'PtfmYaw must be between -'//TRIM(Num2LStr(SmallAngleLimit_Rad))//' and ' &
                                    //TRIM(Num2LStr(SmallAngleLimit_Rad))//' radians.',ErrStat,ErrMsg,RoutineName)
   END IF

      ! Check the output parameters:
   IF ( InputFileData%DecFact < 1_IntKi )  CALL SetErrStat( ErrID_Fatal, 'DecFact must be greater than 0.',ErrStat,ErrMsg,RoutineName )

   IF ( ( InputFileData%NTwGages < 0_IntKi ) .OR. ( InputFileData%NTwGages > 9_IntKi ) )  THEN
      CALL SetErrStat( ErrID_Fatal, 'NTwGages must be between 0 and 9 (inclusive).',ErrStat,ErrMsg,RoutineName)
   ELSE
         ! Check to see if all TwrGagNd(:) analysis points are existing analysis points:

      DO I=1,InputFileData%NTwGages
         IF ( InputFileData%TwrGagNd(I) < 1_IntKi .OR. InputFileData%TwrGagNd(I) > InputFileData%TwrNodes ) THEN
            CALL SetErrStat( ErrID_Fatal, ' All TwrGagNd values must be between 1 and '//&
                           TRIM( Num2LStr( InputFileData%TwrNodes ) )//' (inclusive).',ErrStat,ErrMsg,RoutineName)
            EXIT ! stop checking this loop
         END IF
      END DO         
   
   END IF
         
         
   IF ( ( InputFileData%NBlGages < 0_IntKi ) .OR. ( InputFileData%NBlGages > 9_IntKi ) )  THEN
      CALL SetErrStat( ErrID_Fatal, 'NBlGages must be between 0 and 9 (inclusive).',ErrStat,ErrMsg,RoutineName)
   ELSE 

   ! Check to see if all BldGagNd(:) analysis points are existing analysis points:

      DO I=1,InputFileData%NBlGages
         IF ( InputFileData%BldGagNd(I) < 1_IntKi .OR. InputFileData%BldGagNd(I) > InputFileData%InpBlMesh(1)%BldNodes ) THEN
            CALL SetErrStat( ErrID_Fatal, ' All BldGagNd values must be between 1 and '//&
                 TRIM( Num2LStr( InputFileData%InpBlMesh(1)%BldNodes ) )//' (inclusive).',ErrStat,ErrMsg,RoutineName)
            EXIT ! stop checking this loop
         END IF
      END DO
      
   END IF

   !bjj: since ED doesn't actually use OutFmt at this point, I'm going to remove this check and warning message
   !!!!   ! Check that InputFileData%OutFmt is a valid format specifier and will fit over the column headings
   !!!!CALL ChkRealFmtStr( InputFileData%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
   !!!!CALL SetErrStat(ErrStat2, ErrMsg2,ErrStat,ErrMsg,RoutineName )
   !!!!IF ( FmtWidth /= ChanLen ) CALL SetErrStat(ErrID_Warn, 'OutFmt produces a column width of '//TRIM(Num2LStr(FmtWidth))//&
   !!!!                                        ' instead of '//TRIM(Num2LStr(ChanLen))//' characters.',ErrStat,ErrMsg,RoutineName )


   RETURN

END SUBROUTINE ValidatePrimaryData
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE ElastoDyn_IO   
!**********************************************************************************************************************************
