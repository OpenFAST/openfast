/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   map[dot]plus[dot]plus[dot]help[at]gmail                     *
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


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


#include "../src/map.h"


int main(int argc, char *argv[]) 
{
  void* none = NULL;
  int its = 0;
  int i = 0;
  char map_msg[1024] = "\0";
  MAP_ERROR_CODE ierr = MAP_SAFE;
  MAP_ERROR_CODE success = MAP_SAFE;

  double depth = 100;
  double g = 9.81;
  double rho = 1020;
  double dt = 0.5;
  double time = 0.0;
  double vessel_y_position = 0.0;
  char** header_array = NULL;
  char** unit_array = NULL;
  
  /* For clarity, all API functions are prefixed with 'map_'. The convention applied before, i.e., 'py_create_*_data()', seemed like
   * the function applied to python only. This isn't true; any function prototytped with MAP_EXTERNCALL is exposed to any lanugage, and 
   * is part of the API library. 
   */
  MAP_InitInputType_t* init_type = (MAP_InitInputType_t*)(uintptr_t)map_create_init_type(map_msg, &ierr); // @todo: check ierr==MAP_SAFE error for all
  MAP_InitOutputType_t* io_type = (MAP_InitOutputType_t*)(uintptr_t)map_create_initout_type(map_msg, &ierr);
  MAP_InputType_t* u_type = (MAP_InputType_t*)(uintptr_t)map_create_input_type(map_msg, &ierr);
  MAP_ParameterType_t* p_type = (MAP_ParameterType_t*)(uintptr_t)map_create_parameter_type(map_msg, &ierr);
  MAP_ConstraintStateType_t* z_type = (MAP_ConstraintStateType_t*)(uintptr_t)map_create_constraint_type(map_msg, &ierr);
  MAP_ContinuousStateType_t* x_type = (MAP_ContinuousStateType_t*)(uintptr_t)map_create_continuous_type(map_msg, &ierr);
  MAP_OutputType_t* y_type = (MAP_OutputType_t*)(uintptr_t)map_create_output_type(map_msg, &ierr);            
  MAP_OtherStateType_t* other_type = (MAP_OtherStateType_t*)(uintptr_t)map_create_other_type(map_msg, &ierr); 

  char prop_def[3][100];
  char node_def[5][100];
  char line_def[4][100];
  char option_def[12][100];

  /* The reason why the init_type needs to be set line-by-line for fortran (ISO_C_BINDING) 
   * legacy reasons. I haven't figured out a clean way to pass 2D char** arrays from fortran 
   * to C (because 2D arrays between C and fortran do no align in memory). 1D arrays do line 
   * up. A better, more generic, non-text driven way of setting model parameters should be use. 
   */
  strcpy(prop_def[0], "upper     0.07      160        600000000        1.0    1.0E8    0.6 -1.0    0.05 \0");
  strcpy(prop_def[1], "middle    0.08      180        700000000        1.0    1.0E8    0.6 -1.0    0.05 \0");
  strcpy(prop_def[2], "lower     0.085     200        800000000        1.0    1.0E8    0.6 -1.0    0.05 \0");
  strcpy(node_def[0], "1     fix       -300     0      depth   0     0      #      #       # \0");
  strcpy(node_def[1], "2     Connect  #-120    #0     #-70     0     0      0      0       0 \0");
  strcpy(node_def[2], "3     Connect  #-53     #0     #-20     0     0      0      0       0 \0");
  strcpy(node_def[3], "4     vessel    -25     15       0.0    0     0      #      #       # \0");
  strcpy(node_def[4], "5     Vessel    -25    -15       0.0    0     0      #      #       # \0");
  strcpy(line_def[0], "1        lower     122       1         2 omit_contact gx_pos \0");
  strcpy(line_def[1], "2        middle    122       2         3 omit_contact gx_pos \0");
  strcpy(line_def[2], "3        upper     50        3         4 omit_contact gx_pos \0");
  strcpy(line_def[3], "4        upper     50        3         5 omit_contact gx_pos \0");
  strcpy(option_def[0], " krylov_accelerator 0\0");
  strcpy(option_def[1], "inner_xtol 1e-9\0");
  strcpy(option_def[2], "outer_tol 1e-5\0");
  strcpy(option_def[3], " pg_cooked 1000 1\0");
  strcpy(option_def[4], "outer_fd\0");
  strcpy(option_def[5], "inner_max_its 500\0");
  strcpy(option_def[6], "outer_max_its 200\0");
  strcpy(option_def[7], " repeat 120 240\0");

  map_initialize_msqs_base(u_type, p_type, x_type, z_type, other_type, y_type, io_type);
  map_set_sea_depth(p_type, depth);
  map_set_gravity(p_type, g);
  map_set_sea_density(p_type, rho);
  
  /* set cable library data */
  strcpy(init_type->library_input_str, prop_def[0]); map_add_cable_library_input_text(init_type); 
  strcpy(init_type->library_input_str, prop_def[1]); map_add_cable_library_input_text(init_type); 
  strcpy(init_type->library_input_str, prop_def[2]); map_add_cable_library_input_text(init_type); //
   
  /* set node data */
  strcpy(init_type->node_input_str, node_def[0]); map_add_node_input_text(init_type);  
  strcpy(init_type->node_input_str, node_def[1]); map_add_node_input_text(init_type);
  strcpy(init_type->node_input_str, node_def[2]); map_add_node_input_text(init_type);//
  strcpy(init_type->node_input_str, node_def[3]); map_add_node_input_text(init_type);//
  strcpy(init_type->node_input_str, node_def[4]); map_add_node_input_text(init_type);//
  
  /* set line properties */
  strcpy(init_type->line_input_str, line_def[0]); map_add_line_input_text(init_type);  
  strcpy(init_type->line_input_str, line_def[1]); map_add_line_input_text(init_type); //
  strcpy(init_type->line_input_str, line_def[2]); map_add_line_input_text(init_type);  //
  strcpy(init_type->line_input_str, line_def[3]); map_add_line_input_text(init_type);  //
  
  /* set solver options */
  strcpy(init_type->option_input_str, option_def[0]); map_add_options_input_text(init_type);  
  strcpy(init_type->option_input_str, option_def[1]); map_add_options_input_text(init_type);  
  strcpy(init_type->option_input_str, option_def[2]); map_add_options_input_text(init_type);  
  strcpy(init_type->option_input_str, option_def[3]); map_add_options_input_text(init_type);  
  strcpy(init_type->option_input_str, option_def[4]); map_add_options_input_text(init_type);  
  strcpy(init_type->option_input_str, option_def[5]); map_add_options_input_text(init_type);  
  strcpy(init_type->option_input_str, option_def[6]); map_add_options_input_text(init_type);  
  strcpy(init_type->option_input_str, option_def[7]); map_add_options_input_text(init_type);  

  /* 1) read input file, set library input string, ect, in the MAP_InitInput_type_t structure 
   * 2) call map_init(). init_type->object is freed in map_init( ) at the end. 
   * 3) after initialization, delete init_type
   */


  strcpy(init_type->summary_file_name,"baseline.sum.map");  
  map_set_summary_file_name(init_type, map_msg, &ierr);

  map_init(init_type, u_type, p_type, x_type, NULL, z_type, other_type, y_type, io_type, &ierr, map_msg);
  if (ierr!=MAP_SAFE) {
    printf("%s\n",map_msg);
  };

//  /* OPTIONAL: if you want output headers    <-------------------------------------+   */
//  header_array = malloc(sizeof(char*)*(io_type->writeOutputHdr_Len));        //    |
//  unit_array = malloc(sizeof(char*)*(io_type->writeOutputUnt_Len));          //    |
//  for (i=0 ; i<io_type->writeOutputHdr_Len ; i++) {                          //    |
//    header_array[i] = malloc(sizeof(char)*17);                               //    |
//    unit_array[i] = malloc(sizeof(char)*17);                                 //    |
//  };                                                                         //    |
//  map_get_header_string(&io_type->writeOutputHdr_Len,header_array,other_type);//   |
//  map_get_unit_string(&io_type->writeOutputUnt_Len, unit_array, other_type); //    |
//  printf("\t");                                                              //    |
//  for (i=0 ; i<io_type->writeOutputHdr_Len ; i++) {                          //    |
//    printf("%s\t",header_array[i]);                                          //    |
//  };                                                                         //    |
//  printf("\n\t");                                                            //    |
//  for (i=0 ; i<io_type->writeOutputHdr_Len ; i++) {                          //    |
//    printf("%s\t",unit_array[i]);                                            //    |
//  };                                                                         //    |
//  printf("\n");                                                              //    |
//  for (i=0 ; i<io_type->writeOutputHdr_Len ; i++) {                          //    |
//    MAPFREE(header_array[i]);                                                //    |
//    MAPFREE(unit_array[i]);                                                  //    |
//  };                                                                         //    |
//  MAPFREE(header_array);                                                     //    |
//  MAPFREE(unit_array);                                                       //    |  
//  /* ------------------------------------------------------------------------------+   */

  MAPFREE(init_type); 
  MAPFREE(io_type); 

  do {
    /* 1) update the input states in u_type.x, u_type.y, u_type.z (fairlead displacements) 
     *     Alternatively, call py_offset_vessel to displace fairlead positions. 
     *  2) call map_update_states()       
     *  3) call map_calc_output() (optional is you don't want outlist to be updated)
     *  4) get outputs from y_type.fx, y_type.fy, y_type.fz (fairlead node sum-force). You have to calculate
     *     the moments manually, i.e., cross(r,f).
     */
    
    time = (double)its*dt;
//    vessel_y_position = its;
//    
//    u_type->y[0] = vessel_y_position;
//    /* vessel position:                   X    Y                  Z    phi  the  psi */
//    // map_offset_vessel(other_type, u_type, 0.0, vessel_y_position, 0.0, 0.0, 0.0, 0.0, map_msg, &ierr); // @todo: <----- inconsistent argument order, map_msg, ierr?
    map_update_states(time, its, u_type, p_type, x_type, NULL, z_type, other_type, &ierr, map_msg);    // @todo: <----- inconsistent argument order, ierr, map_msg?
    if (ierr!=MAP_SAFE) {
      printf("%s\n",map_msg);
    };

//    
//    printf("Time step %0.1f: ", time);
//    printf("Element 0 fairlead force: %0.2f  %0.2f  %0.2f\n", y_type->Fx[0], y_type->Fy[0], y_type->Fz[0]);
//
//    
//    /* OPTIONAL: if you want output headers    <-------------------------------------+   */
    map_calc_output(time, u_type, p_type, x_type, NULL, z_type,                //    |
                    other_type, y_type, &ierr, map_msg);                       //    |
    if (ierr!=MAP_SAFE) {
      printf("%s\n",map_msg);
    };

//    printf("\t");                                                              //    |
//    for (i=0 ; i<y_type->wrtOutput_Len ; i++) {                                //    |
//      printf("%0.2f\t",y_type->wrtOutput[i]);                                  //    |
//    };                                                                         //    |
//    printf("\n");                                                              //    |
//    /* ------------------------------------------------------------------------------+   */
//
//
    its++;
  } while (its<10); 
   
  /* first delete internal data type, then delete derived data types...  */
  map_end(u_type, p_type, x_type, NULL, z_type, other_type, y_type, &ierr, map_msg);
  success = map_free_types(u_type, p_type, x_type, z_type, other_type, y_type);  /* @todo: check success... why do I do this? nothing should return from API functions; 
                                                                                    arguments are passed by reference. Should pass ierr, map_msg as argument */

  /* 
     Need to make addtional deallocations in C and python; 
  */
  MAPFREE(other_type); 
  MAPFREE(y_type); 
  MAPFREE(u_type);
  MAPFREE(p_type);
  MAPFREE(z_type);
  MAPFREE(x_type);

  return 0;
}
