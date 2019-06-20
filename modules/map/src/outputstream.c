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


/**
 * @file 
 */


#include "outputstream.h"


// <<<<<<< HEAD
// /** @addtogroup FortranCall */
// /* @{ */
// MAP_EXTERNCALL void set_summary_file_name(MAP_InitInputType_t* init_type, char *map_msg, MAP_ERROR_CODE *ierr) {  
//   int length = 0;
//   char *temp = NULL;
//   InitializationData* init = init_type->object; 
//   
//   if (init->summaryFileName!=NULL) {
//     MAPFREE(init->summaryFileName);
//   };  
//   temp = strtok(init_type->summaryFileName, " \n&");
//   length = strlen(temp)+1;
//   init->summaryFileName = (char*)malloc(sizeof(char)*length);
//   strcpy(init->summaryFileName, temp);
//   init_type->summaryFileName[0] = 0;  
// };
// 
// 
// MAP_EXTERNCALL void get_header_string(int* n, char** str_array, MAP_OtherStateType_t* other_type)
// { 
//   int count = 0;    
//   ModelData* model_data = other_type->object;
//   VarTypePtr* vartype_ptr = NULL;
//   VarType* vartype = NULL;
// 
//   list_iterator_start(&model_data->yList->out_list_ptr);
//   while (list_iterator_hasnext(&model_data->yList->out_list_ptr)) { 
//     vartype_ptr = (VarTypePtr*)list_iterator_next(&model_data->yList->out_list_ptr);
//     strcpy(str_array[count],vartype_ptr->name);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list_ptr);     
// 
//   list_iterator_start(&model_data->yList->out_list);
//   while (list_iterator_hasnext(&model_data->yList->out_list)) { 
//     vartype = (VarType*)list_iterator_next(&model_data->yList->out_list);
//     strcpy(str_array[count],vartype->name);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list);     
//   /* @todo this should raise and error when count != n */
// };
// 
// 
// MAP_EXTERNCALL void get_unit_string(int* n, char** str_array, MAP_OtherStateType_t* other_type)
// { 
//   int count = 0;    
//   ModelData* model_data = other_type->object;
//   VarTypePtr* vartype_ptr = NULL;
//   VarType* vartype = NULL;
// 
//   list_iterator_start(&model_data->yList->out_list_ptr);
//   while (list_iterator_hasnext(&model_data->yList->out_list_ptr)) { 
//     vartype_ptr = (VarTypePtr*)list_iterator_next(&model_data->yList->out_list_ptr );
//     strcpy(str_array[count],vartype_ptr->units);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list_ptr);     
// 
//   list_iterator_start(&model_data->yList->out_list);
//   while (list_iterator_hasnext(&model_data->yList->out_list)) { 
//     vartype = (VarType*)list_iterator_next(&model_data->yList->out_list );
//     strcpy(str_array[count],vartype->units);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list);     
// };
// /* @} */
// 
// 
// =======
// >>>>>>> next
MAP_ERROR_CODE get_iteration_output_stream(MAP_OutputType_t *y_type, MAP_OtherStateType_t* other_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int count = 0;
    int i = 0;
  Domain* domain = other_type->object;
  VarTypePtr* vartype_ptr = NULL;
  VarType* vartype = NULL;
  
  /* this should only be run once at initialization */
  if (!y_type->wrtOutput) { /* if NULL, then wrtOutput is not initialized */
    /* first set size of out list */
    y_type->wrtOutput_Len = (int)list_size(&domain->y_list->out_list_ptr);
    y_type->wrtOutput_Len += (int)list_size(&domain->y_list->out_list);
    
    /* allocation of the array */
    y_type->wrtOutput = malloc(sizeof(double)*y_type->wrtOutput_Len);
  };

  list_iterator_start(&domain->y_list->out_list_ptr);
  while (list_iterator_hasnext(&domain->y_list->out_list_ptr)) { 
    vartype_ptr = (VarTypePtr*)list_iterator_next(&domain->y_list->out_list_ptr);
    y_type->wrtOutput[count] = *vartype_ptr->value;
    count++;
  };
  list_iterator_stop(&domain->y_list->out_list_ptr);     

  list_iterator_start(&domain->y_list->out_list);
  while (list_iterator_hasnext(&domain->y_list->out_list)) { 
    vartype = (VarType*)list_iterator_next(&domain->y_list->out_list);
    y_type->wrtOutput[count] = vartype->value;
    count++;    
  };
  list_iterator_stop(&domain->y_list->out_list);    

  //printf("C :\n");
  for(i=0 ; i<count ; i++)  {
	// printf("%f , ",y_type->wrtOutput[i]);
	 //printf("%d  %d\n",y_type->wrtOutput_Len,count);
  }//printf("\n");
  return MAP_SAFE;
};


MAP_ERROR_CODE write_summary_file(InitializationData* init, MAP_ParameterType_t* paramType, Domain* data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  FILE* file;
  MAP_ERROR_CODE err = MAP_SAFE;
  char time_buffer[TIME_BUFFER_SIZE] = "\0";
  time_t timer;

# if !defined(_MSC_VER)
  struct tm* tm_info;
  time(&timer);
  tm_info = localtime(&timer);
  strftime(time_buffer, TIME_BUFFER_SIZE, "%A %B %d-%Y at %H:%M:%S %p", tm_info);
# else
  struct tm tm_info;
  time(&timer);
  localtime_s(&tm_info, &timer);
  strftime(time_buffer, TIME_BUFFER_SIZE, "%A %B %d-%Y at %H:%M:%S %p", &tm_info);
# endif

  err = fopen_s(&file, (char*)init->summary_file_name->data, "w");
  if (err!=MAP_SAFE) {
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_95, "File name: <%s>", init->summary_file_name->data);
    return MAP_FATAL;
  };
  
  if (file==NULL) {
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_38, "File name: <%s>", init->summary_file_name->data);
    return MAP_FATAL;
  };
  
  //__get_machine_name(name);  

  /* May want to retag version number of git between major/minor releases. A trick to include git revision number
   * can be found archived here: http://stackoverflow.com/questions/1704907/how-can-i-get-my-c-code-to-automatically-print-out-its-git-version-hash
   *
   * $ git tag -a v1.00 -m "Release v1.00"
   */

  /* unboubtedly this will be useful on windows. Not sure how to track git version 
   * on windows with MSVC compiler, or for those who don't use version control? */
#  ifndef GITVERSION 
#    define GITVERSION PROGVERSION
#  endif

  fprintf(file,"Outputs were generated using MAP++ %s.\n\n", GITVERSION); 
  fprintf(file,"    %s\n\n", time_buffer); 
  fprintf(file,"    Gravity constant          [m/s^2]  : %1.2f\n", paramType->g); 
  fprintf(file,"    Sea density               [kg/m^3] : %1.2f\n", paramType->rho_sea);
  fprintf(file,"    Water depth               [m]      : %1.2f\n", paramType->depth);
  fprintf(file,"    Vessel reference position [m]      : %1.2f , %1.2f , %1.2f\n\n", data->vessel.ref_origin.x.value, data->vessel.ref_origin.y.value, data->vessel.ref_origin.z.value); 
  
  success = write_cable_library_information_to_summary_file(file, data); CHECKERRK(MAP_ERROR_6);
  success = write_node_information_to_summary_file(file, data, map_msg, ierr); CHECKERRK(MAP_ERROR_7);
  success = write_line_information_to_summary_file(file, data); CHECKERRK(MAP_ERROR_8);
  success = write_expanded_input_file_to_summary_file(file, init);
  fclose(file);  
  
  MAP_RETURN_STATUS(*ierr);
};


MAP_ERROR_CODE write_cable_library_information_to_summary_file(FILE* file, Domain* domain)
{
  CableLibrary* library_iter = NULL;  

  list_iterator_start(&domain->library);    
  while (list_iterator_hasnext(&domain->library)) { 
    library_iter = (CableLibrary*)list_iterator_next(&domain->library);    
    fprintf(file, "    Cable Type          : %s\n", library_iter->label->data);
    fprintf(file, "    Diameter     [m]    : %1.4f\n", library_iter->diam);
    fprintf(file, "    Mass Density [kg/m] : %1.2f\n", library_iter->mass_density);
    fprintf(file, "    EA           [N]    : %1.2f\n", library_iter->EA);
    fprintf(file, "    omega        [N/m]  : %1.2f\n", library_iter->omega);
    fprintf(file, "    CB                  : %1.2f\n\n", library_iter->cb);
  };
  list_iterator_stop(&domain->library);  
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_type_to_summary_file(const int num_col, const int count_to_four, const NodeType node_type, bstring line)
{
  int ret = 0;
  int j = 0;
  const int size = line->slen - 23*count_to_four;

  if (!num_col) { /* if this is the firs column, then line the right side with information */
    if (node_type==VESSEL) {
      ret = bcatcstr(line, "Type      |  VESSEL");
    } else if (node_type==FIX ) { 
      ret = bcatcstr(line, "Type      |  FIX");
    } else if (node_type==CONNECT ) {
      ret = bcatcstr(line, "Type      |  CONNECT");
    } else {
      ret = bcatcstr(line, "Type      |  UNKNOWN");
    };
  } else { /* otherwise, just print the information */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) { /* add white spaces */
      bconchar(line,' ');
    };

    if (node_type==VESSEL) {
      ret = bcatcstr(line, " VESSEL");
    } else if (node_type==FIX ) { 
      ret = bcatcstr(line, " FIX");
    } else if(node_type==CONNECT ) {
      ret = bcatcstr(line, " CONNECT");
    } else {
      ret = bcatcstr(line, " UNKNOWN");
    };
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_header_to_summary_file(const int num_col, const int count_to_four, const int node_num, bstring line)
{
  int ret = 0;
  int j = 0;
  const int size = line->slen - 23*count_to_four;
  bstring info = NULL;

  if (!num_col) { 
    info = bformat("          | Node %d Data", node_num);  
  } else {
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      bconchar(line,' ');
    };
    info = bformat("Node %d Data", node_num);
  };
  ret = bconcat(line, info);
  ret = bdestroy(info);

  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_x_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* x_pos, bstring line)
{
  int j = 0;
  const int size = line->slen-23*count_to_four; /* each node column is 23 character wide */ 
  double x = *x_pos->value;
  bstring info = NULL;

  if (!num_col) { 
    if (!x_pos->is_fixed) {
      if (x>=0.0) {
        info = bformat("X  [m]    | ( %1.3f)", x);
      } else {
        info = bformat("X  [m]    | (%1.3f)", x);
      };
    } else {
      if (x>0.0) {
        info = bformat("X  [m]    |   %1.3f", x);
      } else {
        info = bformat("X  [m]    |  %1.3f", x);
      };
    };
  } else { 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      bconchar(line,' ');
    };

    if (!x_pos->is_fixed) {
      if (x>=0.0) {
        info = bformat("( %1.3f)", x);
      } else {
        info = bformat("(%1.3f)", x);
      };
    } else {
      if (x>=0.0) {
        info = bformat("  %1.3f", x);
      } else {
        info = bformat(" %1.3f", x);
      };
    };
  };

  bconcat(line, info);
  bdestroy(info);
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_y_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* y_pos, bstring line)
{
  int j = 0;
  const int size = line->slen-23*count_to_four; /* each node column is 23 character wide */ 
  double y = *y_pos->value;
  bstring info = NULL;

  if (!num_col) { 
    if (!y_pos->is_fixed) {
      if (y>=0.0) {
        info = bformat("Y  [m]    | ( %1.3f)", y);
      } else {
        info = bformat("Y  [m]    | (%1.3f)", y);
      };
    } else {
      if (y>0.0) {
        info = bformat("Y  [m]    |   %1.3f", y);
      } else {
        info = bformat("Y  [m]    |  %1.3f", y);
      };
    };
  } else { 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      bconchar(line,' ');
    };

    if (!y_pos->is_fixed) {
      if (y>=0.0) {
        info = bformat("( %1.3f)", y);
      } else {
        info = bformat("(%1.3f)", y);
      };
    } else {
      if (y>=0.0) {
        info = bformat("  %1.3f", y);
      } else {
        info = bformat(" %1.3f", y);
      };
    };
  };

  bconcat(line, info);
  bdestroy(info);
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_z_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* z_pos, bstring line)
{
  int j = 0;
  const int size = line->slen-23*count_to_four; /* each node column is 23 character wide */ 
  double z = *z_pos->value;
  bstring info = NULL;

  if (!num_col) { 
    if (!z_pos->is_fixed) {
      if (z>=0.0) {
        info = bformat("Z  [m]    | ( %1.3f)", z);
      } else {
        info = bformat("Z  [m]    | (%1.3f)", z);
      };
    } else {
      if (z>0.0) {
        info = bformat("Z  [m]    |   %1.3f", z);
      } else {
        info = bformat("Z  [m]    |  %1.3f", z);
      };
    };
  } else { 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      bconchar(line,' ');
    };

    if (!z_pos->is_fixed) {
      if (z>=0.0) {
        info = bformat("( %1.3f)", z);
      } else {
        info = bformat("(%1.3f)", z);
      };
    } else {
      if (z>=0.0) {
        info = bformat("  %1.3f", z);
      } else {
        info = bformat(" %1.3f", z);
      };
    };
  };

  bconcat(line, info);
  bdestroy(info);
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_mass_information_to_summary_file(const int num_col, const int count_to_four, VarType* point_mass, bstring line)
{
  int j = 0;
  const int size = line->slen-23*count_to_four; /* each node column is 23 character wide */ 
  const double mass = point_mass->value;
  bstring info = NULL;

  if (!num_col) { 
    if (!point_mass->is_fixed) {
      if (mass>=0.0) {
        info = bformat("M  [kg]   | ( %1.3f)", mass);
      } else {
        info = bformat("M  [kg]   | (%1.3f)", mass);                       
      };
    } else {
      if (mass>0.0) {
        info = bformat("M  [kg]   |   %1.3f", mass);
      } else {
        info = bformat("M  [kg]   |  %1.3f", mass);
      };
    };
  } else { 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      bconchar(line,' ');
    };

    if (!point_mass->is_fixed) {
      if (mass>=0.0) {
        info = bformat("( %1.3f)", mass);
      } else {
        info = bformat("(%1.3f)", mass);
      };
    } else {
      if (mass>=0.0) {
        info = bformat("  %1.3f", mass);
      } else {
        info = bformat(" %1.3f", mass);
      };
    };
  };

  bconcat(line, info);
  bdestroy(info);
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_buoyancy_information_to_summary_file(const int num_col, const int count_to_four, VarType* point_buoy, bstring line)
{
  int j = 0;
  const int size = line->slen-23*count_to_four; /* each node column is 23 character wide */ 
  const double buoy = point_buoy->value;
  bstring info = NULL;
  
  if (!num_col) { 
    if (!point_buoy->is_fixed) {
      if (buoy>=0.0) {
        info = bformat("B  [m^3]  | ( %1.3f)", buoy);                       
      } else {
        info = bformat("B  [m^3]  | (%1.3f)", buoy);                       
      };
    } else {
      if (buoy>0.0) {
        info = bformat("B  [m^3]  |   %1.3f", buoy);
      } else {
        info = bformat("B  [m^3]  |  %1.3f", buoy);
      };
    };
  } else { 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      bconchar(line,' ');
    };

    if (!point_buoy->is_fixed) {
      if (buoy>=0.0) {
        info = bformat("( %1.3f)", buoy);
      } else {
        info = bformat("(%1.3f)", buoy);
      };
    } else {
      if (buoy>=0.0) {
        info = bformat("  %1.3f", buoy);
      } else {
        info = bformat(" %1.3f", buoy);
      };
    };
  };

  bconcat(line, info);
  bdestroy(info);
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_x_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* x_sum_force, bstring line)
{
  int j = 0;
  const int size = line->slen-23*count_to_four; /* each node column is 23 character wide */ 
  const double fx = *(x_sum_force->value);
  bstring info = NULL;

  if (!num_col) { 
    if (!x_sum_force->is_fixed) {
      if (fx>=0.0) {
        info = bformat("FX [N]    | ( %1.3f)", fx);
      } else {
        info = bformat("FX [N]    | (%1.3f)", fx);
      };
    } else {
      if (fx>0.0) {
        info = bformat("FX [N]    |   %1.3f", fx);
      } else {
        info = bformat("FX [N]    |  %1.3f", fx);
      };
    };
  } else { 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      bconchar(line,' ');
    };

    if (!x_sum_force->is_fixed) {
      if (fx>=0.0) {
        info = bformat("( %1.3f)", fx);
      } else {
        info = bformat("(%1.3f)", fx);
      };
    } else {
      if (fx>=0.0) {
        info = bformat("  %1.3f", fx);
      } else {
        info = bformat(" %1.3f", fx);
      };
    };
  };

  bconcat(line, info);
  bdestroy(info);
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_y_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* y_sum_force, bstring line)
{
  int j = 0;
  const int size = line->slen-23*count_to_four; /* each node column is 23 character wide */ 
  const double fy = *(y_sum_force->value);
  bstring info = NULL;

  if (!num_col) { 
    if (!y_sum_force->is_fixed) {
      if (fy>=0.0) {
        info = bformat("FY [N]    | ( %1.3f)", fy);
      } else {
        info = bformat("FY [N]    | (%1.3f)", fy);
      };
    } else {
      if (fy>0.0) {
        info = bformat("FY [N]    |   %1.3f", fy);
      } else {
        info = bformat("FY [N]    |  %1.3f", fy);
      };
    };
  } else { 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      bconchar(line,' ');
    };

    if (!y_sum_force->is_fixed) {
      if (fy>=0.0) {
        info = bformat("( %1.3f)", fy);
      } else {
        info = bformat("(%1.3f)", fy);
      };
    } else {
      if (fy>=0.0) {
        info = bformat("  %1.3f", fy);
      } else {
        info = bformat(" %1.3f", fy);
      };
    };
  };

  bconcat(line, info);
  bdestroy(info);
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_z_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* z_sum_force, bstring line)
{
  int j = 0;
  const int size = line->slen-23*count_to_four; /* each node column is 23 character wide */ 
  const double fz = *(z_sum_force->value);
  bstring info = NULL;

  if (!num_col) { 
    if (!z_sum_force->is_fixed) {
      if (fz>=0.0) {
        info = bformat("FZ [N]    | ( %1.3f)", fz);
      } else {
        info = bformat("FZ [N]    | (%1.3f)", fz);
      };
    } else {
      if (fz>0.0) {
        info = bformat("FZ [N]    |   %1.3f", fz);
      } else {
        info = bformat("FZ [N]    |  %1.3f", fz);
      };
    };
  } else { 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      bconchar(line,' ');
    };

    if (!z_sum_force->is_fixed) {
      if (fz>=0.0) {
        info = bformat("( %1.3f)", fz);
      } else {
        info = bformat("(%1.3f)", fz);
      };
    } else {
      if (fz>=0.0) {
        info = bformat("  %1.3f", fz);
      } else {
        info = bformat(" %1.3f", fz);
      };
    };
  };

  bconcat(line, info);
  bdestroy(info);
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_information_to_summary_file(FILE* file, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int num = 0;
  unsigned int i = 0;
  int col_cnt = 0; /* number of columns. The MAP summary file prints 4 columns of node information. After column 4, carriage return */
  
  bstring line0 = NULL; 
  bstring line1 = NULL; 
  bstring line2 = NULL; 
  bstring line3 = NULL; 
  bstring line4 = NULL; 
  bstring line5 = NULL;
  bstring line6 = NULL;
  bstring line7 = NULL;
  bstring line8 = NULL;
  bstring line9 = NULL;
  Node* node_iter = NULL;  
  const int FOUR = 4;
  const unsigned int num_nodes = list_size(&domain->node);  
  unsigned int col = 0;
  MAP_ERROR_CODE success = MAP_SAFE;

  MAP_BEGIN_ERROR_LOG;

  for (i=0 ; i<num_nodes ; i+=FOUR) {
    if (i+FOUR>num_nodes) {
      num = num_nodes-i;
    } else {
      num = FOUR;
    };    

    line0 = bformat("");
    line1 = bformat("");
    line2 = bformat("");
    line3 = bformat("");
    line4 = bformat("");
    line5 = bformat("");
    line6 = bformat("");
    line7 = bformat("");
    line8 = bformat("");
    line9 = bformat("");

    for (col=i ; col<i+num ; col++) {
      node_iter = (Node*)list_get_at(&domain->node, col);      
      success = write_node_header_to_summary_file(col-i, col_cnt, col+1, line0); CHECKERRQ(MAP_FATAL_70);
      success = write_node_type_to_summary_file(col-i, col_cnt, node_iter->type, line1); CHECKERRQ(MAP_FATAL_70);
      success = write_node_x_position_to_summary_file(col-i, col_cnt, &node_iter->position_ptr.x, line2); CHECKERRQ(MAP_FATAL_70);
      success = write_node_y_position_to_summary_file(col-i, col_cnt, &node_iter->position_ptr.y, line3); CHECKERRQ(MAP_FATAL_70);
      success = write_node_z_position_to_summary_file(col-i, col_cnt, &node_iter->position_ptr.z, line4); CHECKERRQ(MAP_FATAL_70);
      success = write_node_mass_information_to_summary_file(col-i, col_cnt, &node_iter->M_applied, line5); CHECKERRQ(MAP_FATAL_70);
      success = write_node_buoyancy_information_to_summary_file(col-i, col_cnt, &node_iter->B_applied, line6); CHECKERRQ(MAP_FATAL_70);
      success = write_node_x_sum_force_to_summary_file(col-i, col_cnt, &node_iter->sum_force_ptr.fx, line7); CHECKERRQ(MAP_FATAL_70);
      success = write_node_y_sum_force_to_summary_file(col-i, col_cnt, &node_iter->sum_force_ptr.fy, line8); CHECKERRQ(MAP_FATAL_70);
      success = write_node_z_sum_force_to_summary_file(col-i, col_cnt, &node_iter->sum_force_ptr.fz, line9); CHECKERRQ(MAP_FATAL_70);
      col_cnt++;
    };
    col_cnt = 0;
      
    fprintf(file, "%s\n",line0->data);
    fprintf(file, "          | -------------------------------------------------------------------------------------------\n");
    fprintf(file, "%s\n",line1->data);
    fprintf(file, "%s\n",line2->data);
    fprintf(file, "%s\n",line3->data);
    fprintf(file, "%s\n",line4->data);
    fprintf(file, "%s\n",line5->data);
    fprintf(file, "%s\n",line6->data);
    fprintf(file, "%s\n",line7->data);
    fprintf(file, "%s\n",line8->data);
    fprintf(file, "%s\n\n\n",line9->data);
      
    bdestroy(line0);
    bdestroy(line1);
    bdestroy(line2);
    bdestroy(line3);
    bdestroy(line4);
    bdestroy(line5);
    bdestroy(line6);
    bdestroy(line7);
    bdestroy(line8);
    bdestroy(line9);
  };

  MAP_END_ERROR_LOG;

  MAP_RETURN_STATUS(*ierr);
};


MAP_ERROR_CODE write_line_information_to_summary_file(FILE* file, Domain* domain)
{
  const unsigned int num_lines = list_size(&domain->line);
  Line* line_iter = NULL;  
  bstring line0 = NULL;
  bstring line1 = NULL;
  bstring line2 = NULL;
  bstring line3 = NULL;
  bstring line4 = NULL;
  bstring line5 = NULL;
  bstring line6 = NULL;
  bstring line7 = NULL;
  bstring line8 = NULL;
  bstring line9 = NULL;
  bstring line10 = NULL;
  unsigned int i = 0;  

  for (i=0 ; i<num_lines ; i++) {
    line_iter = (Line*)list_get_at(&domain->line, i);

    line0 = bformat("");
    line1 = bformat("");
    line2 = bformat("");
    line3 = bformat("");
    line4 = bformat("");
    line5 = bformat("");
    line6 = bformat("");
    line7 = bformat("");
    line8 = bformat("");
    line9 = bformat("");

    if (line_iter->Lu.value>0.0) {
      bconchar(line0,' ');      
    };
    if (line_iter->Lb>0.0) {
      bconchar(line1,' ');      
    };
    if (*(line_iter->H.value)>0.0) {
      bconchar(line2,' ');      
    };
    if (*(line_iter->V.value)>0.0) {
      bconchar(line3,' ');      
    };
    if (line_iter->T>0.0) {
      bconchar(line4,' ');      
    };
    if (line_iter->alpha>0.0) {
      bconchar(line5,' ');      
    };
    if (line_iter->H_at_anchor>0.0) {
      bconchar(line6,' ');      
    };
    if (line_iter->V_at_anchor>0.0) {
      bconchar(line7,' ');      
    };
    if (line_iter->T_at_anchor>0.0) {
      bconchar(line8,' ');      
    };
    if (line_iter->alpha_at_anchor>0.0) {
      bconchar(line9,' ');      
    };
        
    switch (line_iter->converge_reason) {
    case 0 :
      line10 = bformat("Improver input parameter for inner loop solver.");
      break;
    case 1 :
      line10 = bformat("SQRT( H^2 + V^2 ) <= INNER_FTOL.");
      break;
    case 2 :
      line10 = bformat("||x_i - x_{i-1}||  <= INNER_XTOL.");
      break;
    case 3 :
      line10 = bformat("Both SQRT( H^2 + V^2 )<=INNER_FTOL and ||x_i - x_{i-1}||  <= INNER_XTOL.");
      break;
    case 4 :
      line10 = bformat("The cosine of the angle between [H, V] and any column of the Jacobian is at most INNER_GTOL in absolute value.");
      break;
    case 5 :
      line10 = bformat("INNER_MAX_ITS reached. Try increasing the maximum allowable iteration count.");
      break;
    case 6 :
      line10 = bformat("INNER_FTOL si too small. SQRT( H^2 + V^2 ) cannot be reduced further.");
      break;
    case 7 :
      line10 = bformat("INNER_XTOL is too small. No further refinements can be made in [H, V].");
      break;
    case 8 :
      line10 = bformat("INNER_GTOL is too small. The iterated variable vector [H, V] is orthogonal to the Jacobian columns.");
      break;
    default :
      line10 = bformat("Uncaught error.");
      break;
    };

    fprintf(file, "                | Line %d\n", i+1);    
    fprintf(file, "                | ---------------------------------------\n");    
    fprintf(file, "Material        |  %s\n", line_iter->line_property->label->data);
    fprintf(file, "Lu        [m]   | %s%1.3f\n", line0->data, line_iter->Lu.value);
    fprintf(file, "Lb        [m]   | %s%1.3f\n", line1->data, line_iter->Lb); 
    fprintf(file, "H         [N]   | %s%1.3f\n", line2->data, *(line_iter->H.value));
    fprintf(file, "V         [N]   | %s%1.3f\n", line3->data, *(line_iter->V.value));
    fprintf(file, "T         [N]   | %s%1.3f\n", line4->data, line_iter->T);
    fprintf(file, "Alpha     [deg] | %s%1.3f\n", line5->data, line_iter->alpha*RAD2DEG); 
    fprintf(file, "HAnch     [N]   | %s%1.3f\n", line6->data, line_iter->H_at_anchor);
    fprintf(file, "VAnch     [N]   | %s%1.3f\n", line7->data, line_iter->V_at_anchor);
    fprintf(file, "TAnch     [N]   | %s%1.3f\n", line8->data, line_iter->T_at_anchor);
    fprintf(file, "AlphaAnch [deg] | %s%1.3f\n", line9->data, line_iter->alpha_at_anchor*RAD2DEG);
    fprintf(file, "L^2-Norm        |  %1.7g\n", line_iter->residual_norm);
    fprintf(file, "Function Evals  |  %d\n", line_iter->evals);
    fprintf(file, "Jacobian Evals  |  %d\n", line_iter->njac_evals);
    fprintf(file, "Term. criteria  |  %d : %s\n", line_iter->converge_reason, line10->data);
    fprintf(file, "\n\n");

    bdestroy(line0);
    bdestroy(line1);
    bdestroy(line2);
    bdestroy(line3);
    bdestroy(line4);
    bdestroy(line5);
    bdestroy(line6);
    bdestroy(line7);
    bdestroy(line8);
    bdestroy(line9);
    bdestroy(line10);
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_expanded_input_file_to_summary_file(FILE* file, InitializationData* init_data)
{
  const char* line = NULL;
  int len = 0;
  int i = 0;

  fprintf(file, "---------------------- LINE DICTIONARY ---------------------------------------\n");
  fprintf(file, "LineType  Diam      MassDenInAir   EA            CB\n");
  fprintf(file, "[-]       [m]       [kg/m]         [N]           [-]\n");  
  for (i=0 ; i<init_data->library_input_string->qty ; i++) {
    line = (char*)init_data->library_input_string->entry[i]->data;
    len = init_data->library_input_string->entry[i]->slen;
    fprintf(file, "%s", line);
    if (line[len-1]!='\n') {
      fprintf(file, "\n" );
    };
  };

  fprintf(file, "---------------------- NODE PROPERTIES ---------------------------------------\n");
  fprintf(file, "Node  Type       X       Y       Z      M     B     FX    FY    FZ\n");
  fprintf(file, "[-]   [-]       [m]     [m]     [m]    [kg]  [mˆ3]  [N]   [N]   [N]\n");
  for (i=0 ; i<init_data->expanded_node_input_string->qty ; i++) {    
    line = (char*)init_data->expanded_node_input_string->entry[i]->data;
    len = init_data->expanded_node_input_string->entry[i]->slen;
    fprintf(file, "%s", line);
    if (line[len-1]!='\n') {
      fprintf(file, "\n" );
    };
  };
  
  fprintf(file, "---------------------- LINE PROPERTIES ---------------------------------------\n");
  fprintf(file, "Line  LineType  UnstrLen  NodeAnch  NodeFair  Flags\n");
  fprintf(file, "[-]      [-]       [m]       [-]       [-]       [-]\n");
  for (i=0 ; i<init_data->expanded_line_input_string->qty ; i++) {    
    line = (char*)init_data->expanded_line_input_string->entry[i]->data;
    len = init_data->expanded_line_input_string->entry[i]->slen;
    fprintf(file, "%s", line);
    if (line[len-1]!='\n') {
      fprintf(file, "\n" );
    };
  };
  
  fprintf(file, "---------------------- SOLVER OPTIONS-----------------------------------------\n");
  fprintf(file, "Option \n");
  fprintf(file, "[-] \n");
  for (i=0 ; i<init_data->solver_options_string->qty ; i++) {    
    line = (char*)init_data->solver_options_string->entry[i]->data;
    len = init_data->solver_options_string->entry[i]->slen;
    fprintf(file, "%s", line);
    if (line[len-1]!='\n') {
      fprintf(file, "\n" );
    };
  };
  return MAP_SAFE;
};


#if !defined(_MSC_VER) && !defined(__MINGW32__)
MAP_ERROR_CODE fopen_s(FILE** f, const char* name, const char* mode)
{
  assert(f);
  *f = fopen(name, mode);
  if (!*f) {
    return MAP_FATAL;   
  }
  return MAP_SAFE;
};
#endif
