/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   map[dot]plus[dot]plus[dot]help[at]gmail                    *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


#ifndef _MAPSYS_H
#define _MAPSYS_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>


#if defined(_WIN32) || defined(_WIN64)
#  include <Windows.h>
#  include <tchar.h>
#else
#  include <unistd.h>
#endif



#if defined _WIN32 || defined __CYGWIN__ || defined _WINDOWS
#  ifdef BUILD_SHARED_LIC_LIBS
#    ifdef __GNUC__
#      define DLL_PUBLIC __attribute__ ((dllexport))
#    else
#      define DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
#    endif
#  else
#    define DLL_PUBLIC
#  endif
#else
#  if __GNUC__ >= 4
#    define DLL_PUBLIC __attribute__ ((visibility ("default")))
#  else
#    define DLL_PUBLIC
#  endif
#endif

#define MAP_EXTERNCALL DLL_PUBLIC

#if defined(_MSC_VER)
  typedef int bool;
  #define false 0
  #define true 1
//#  include "stdbool.h"
#  define map_snprintf _snprintf
#  define map_strcat(a,b,c) strcat_s(a,b,c)
#  define MAP_STRCPY(a,b,c) strcpy_s(a,b,c)
#else
#  include <stdbool.h>
#  define map_snprintf snprintf
#  define map_strcat(a,b,c) strncat(a,c,b)
#  define MAP_STRCPY(a,b,c) strcpy(a,c)
#endif


#ifndef BUILD_DEFS_H
#  define BUILD_DEFS_H
#  define BUILD_YEAR_CH0 (__DATE__[ 7])
#  define BUILD_YEAR_CH1 (__DATE__[ 8])
#  define BUILD_YEAR_CH2 (__DATE__[ 9])
#  define BUILD_YEAR_CH3 (__DATE__[10])
#  define BUILD_MONTH_IS_JAN (__DATE__[0] == 'J' && __DATE__[1] == 'a' && __DATE__[2] == 'n')
#  define BUILD_MONTH_IS_FEB (__DATE__[0] == 'F')
#  define BUILD_MONTH_IS_MAR (__DATE__[0] == 'M' && __DATE__[1] == 'a' && __DATE__[2] == 'r')
#  define BUILD_MONTH_IS_APR (__DATE__[0] == 'A' && __DATE__[1] == 'p')
#  define BUILD_MONTH_IS_MAY (__DATE__[0] == 'M' && __DATE__[1] == 'a' && __DATE__[2] == 'y')
#  define BUILD_MONTH_IS_JUN (__DATE__[0] == 'J' && __DATE__[1] == 'u' && __DATE__[2] == 'n')
#  define BUILD_MONTH_IS_JUL (__DATE__[0] == 'J' && __DATE__[1] == 'u' && __DATE__[2] == 'l')
#  define BUILD_MONTH_IS_AUG (__DATE__[0] == 'A' && __DATE__[1] == 'u')
#  define BUILD_MONTH_IS_SEP (__DATE__[0] == 'S')
#  define BUILD_MONTH_IS_OCT (__DATE__[0] == 'O')
#  define BUILD_MONTH_IS_NOV (__DATE__[0] == 'N')
#  define BUILD_MONTH_IS_DEC (__DATE__[0] == 'D')
#  define BUILD_MONTH_CH0 (__DATE__[ 0])
#  define BUILD_MONTH_CH1 (__DATE__[ 1])
#  define BUILD_MONTH_CH2 (__DATE__[ 2])
#  define BUILD_DAY_CH0 ((__DATE__[4] >= '0') ? (__DATE__[4]) : '0')
#  define BUILD_DAY_CH1 (__DATE__[ 5])
#endif // BUILD_DEFS_H


#ifdef DEBUG
#  define checkpoint() printf("Checkpoint: Line %d in file %s\n",__LINE__,__FILE__);
#else
#  define checkpoint() 
#endif // DEBUG


#define MAX_INIT_TYPE_STRING_LENGTH 255
#define TIME_BUFFER_SIZE 64
#define MAX_INIT_VERSION_STRING_LENGTH 99
#define MAX_INIT_COMPILING_DATA_STRING_LENGTH 25
#define MAP_ERROR_STRING_LENGTH 1024
#define MAP_HORIZONTAL_TOL 1E-2

#define PROGNAME "MAP++ (Mooring Analysis Program++)"
#define PROGVERSION "1.20.10"
#define CHECKERRQ(code) if(success!=MAP_SAFE) {set_universal_error(map_msg, ierr, code); break;} 
#define CHECKERRK(code) if(success!=MAP_SAFE) {set_universal_error(map_msg, ierr, code);} 
#define MAPFREE(obj) if(obj!=NULL) {free(obj); obj=NULL;} 
#define DEG2RAD 0.01745329251 /*  pi/180  */
#define RAD2DEG 57.2957795131 /*  180/pi  */
#define ARCSINH(x) log(x+sqrt(1+x*x))
#define SPACE_LENGTH 12
#define MACHINE_EPSILON 1e-16
#define MAP_BEGIN_ERROR_LOG do{ \
  ; 
#define MAP_END_ERROR_LOG } while(0);
#define MAP_RETURN_STATUS(x) \
  if (x==MAP_SAFE) {         \
    return MAP_SAFE;         \
  } else if (x==MAP_ERROR) { \
    return MAP_ERROR;        \
  } else {                   \
    return MAP_FATAL;        \
  };                             


/* Text Coloring (OS dependant)
 *
 * Not used any longer, but can color text at the terminal. If we are on a non-Unix OS, then:
 *   -- set the strings to "" (empty) so that garbage is not printed
 */
#ifdef __posix 
#  define MAP_COLOR_RED "\033[1;31m"
#  define MAP_COLOR_YELLOW "\033[1;33m"
#  define MAP_COLOR_BLUE "\033[1;34m"
#  define MAP_COLOR_END "\033[0m"
#elif __linux
#  define MAP_COLOR_RED "\033[1;31m"
#  define MAP_COLOR_YELLOW "\033[1;33m"
#  define MAP_COLOR_BLUE "\033[1;34m"
#  define MAP_COLOR_END "\033[0m"
#elif __unix
#  define MAP_COLOR_RED "\033[1;31m"
#  define MAP_COLOR_YELLOW "\033[1;33m"
#  define MAP_COLOR_BLUE "\033[1;34m"
#  define MAP_COLOR_END "\033[0m"
#elif __APPLE__
#  define MAP_COLOR_RED "\033[1;31m"
#  define MAP_COLOR_YELLOW "\033[1;33m"
#  define MAP_COLOR_BLUE "\033[1;34m"
#  define MAP_COLOR_END "\033[0m"
#else
#  define MAP_COLOR_RED ""
#  define MAP_COLOR_YELLOW ""
#  define MAP_COLOR_BLUE ""
#  define MAP_COLOR_END ""
#endif

void __get_machine_name( char* machineName );

#endif /* _MAPSYS_H */
