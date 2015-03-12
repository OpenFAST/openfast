/****************************************************************
 *   Copyright (C) 2014 Adrien Combourieu                       *
 *                                                              *
 ****************************************************************/


#include <iostream>

using namespace std;

//NOTE: including C functions into C++
extern "C" {

#include "../src/map.h"
#include "../src/mapinit.h"
#include "../src/outputstream.h"
#include "../src/mapapi.h"
#include "../src/freedata.h"

}


int main(int argc, char **argv)
{
   
  ///////// INITIALIZE with empty objects////////////
  MAP_InitInputType_t* _init_type = NULL;
  MAP_InitOutputType_t* _io_type = NULL;
  MAP_InputType_t* _u_type = NULL;
  MAP_ParameterType_t* _p_type = NULL;
  MAP_ConstraintStateType_t* _z_type = NULL;
  MAP_ContinuousStateType_t* _x_type = NULL;
  MAP_OutputType_t* _y_type = NULL;
  MAP_OtherStateType_t* _other_type = NULL;

  MAP_ERROR_CODE success = MAP_SAFE;  
  char map_msg[MAP_ERROR_STRING_LENGTH] = "\0";
  // char map_msg[255] = "\0";
  MAP_ERROR_CODE ierr = MAP_SAFE;
  
  ///////// ALLOCATE MEMORY ///////////////
  cout << "--------- Start of MEMORY ALLOCATION ----------" << endl;
  
  _init_type = (MAP_InitInputType_t*)(uintptr_t)map_create_init_type(map_msg, &ierr); // @todo: check ierr==MAP_SAFE error for all
  _io_type = (MAP_InitOutputType_t*)(uintptr_t)map_create_initout_type(map_msg, &ierr);
  _u_type = (MAP_InputType_t*)(uintptr_t)map_create_input_type(map_msg, &ierr);
  _p_type = (MAP_ParameterType_t*)(uintptr_t)map_create_parameter_type(map_msg, &ierr);
  _z_type = (MAP_ConstraintStateType_t*)(uintptr_t)map_create_constraint_type(map_msg, &ierr);
  _x_type = (MAP_ContinuousStateType_t*)(uintptr_t)map_create_continuous_type(map_msg, &ierr);
  _y_type = (MAP_OutputType_t*)(uintptr_t)map_create_output_type(map_msg, &ierr);            
  _other_type = (MAP_OtherStateType_t*)(uintptr_t)map_create_other_type(map_msg, &ierr); 
  
  cout << "--------- End of MEMORY ALLOCATION ----------" << endl;
  
  /////////// INIATILIZE THE OBJECTS WITH USER INPUT //////////
  cout << "--------- Start of INITIALIZATION ----------" << endl;
  ////// SET ENVIRONMENT PARAMETERS //////
  double depth = 100.0;
  double g = 9.81;
  double rho = 1025.0;

  map_initialize_msqs_base(_u_type, _p_type, _x_type, _z_type, _other_type, _y_type, _io_type);  
  map_set_sea_depth(_p_type, depth);
  map_set_gravity(_p_type, g);
  map_set_sea_density(_p_type, rho);
    
  
  ////// INPUT FILE PARSER (emulated here) /////////
  char line_def[2][255];
  char node_def[2][255];
  char element_def[1][255];
  char option_def[7][255];
 
  strcpy(line_def[0], "chainup   0.064929    138.73 545000000 1.0 1.0E8 0.6 -1.0 0.05\0");
  strcpy(line_def[1], "chainup2  0.064929    138.73 545000000 1.0 1.0E8 0.6 -1.0 0.05\0");
 
  strcpy(node_def[0], "1  fix    -50    0  -100  0 0 # # #\0");
  strcpy(node_def[1], "2  Vessel  0     0   0    0 0 # # #\0");
 
  //   strcpy(element_def[0], "1 chainup  100 1 3 gy_pos h_fair\0");
  //   strcpy(element_def[1], "2 chainup2  100 2 4 gy_pos h_fair\0");
  strcpy(element_def[0], "1 chainup  125 1 2 gy_pos h_fair\0");
  
  strcpy(option_def[0], "inner_ftol 1e-6\0");
  strcpy(option_def[1], "inner_gtol 1e-12\0");
  strcpy(option_def[2], "inner_xtol 1e-6\0");
  strcpy(option_def[3], "outer_tol 1e-4\0");
  strcpy(option_def[4], "inner_max_its 400\0");
  strcpy(option_def[5], "outer_max_its 400\0");
  strcpy(option_def[6], "repeat 12 280\0");
  
  /* set cable library data */
  strcpy(_init_type->library_input_str, line_def[0]); map_add_cable_library_input_text(_init_type); 
  strcpy(_init_type->library_input_str, line_def[1]); map_add_cable_library_input_text(_init_type); 
   
  /* set node data */
  strcpy(_init_type->node_input_str, node_def[0]); map_add_node_input_text(_init_type);  
  strcpy(_init_type->node_input_str, node_def[1]); map_add_node_input_text(_init_type);
  
  /* set line properties */
  strcpy(_init_type->line_input_str, element_def[0]); map_add_line_input_text(_init_type);  
  
  /* set solver options */
  strcpy(_init_type->option_input_str, option_def[0]); map_add_options_input_text(_init_type);  
  strcpy(_init_type->option_input_str, option_def[1]); map_add_options_input_text(_init_type);  
  strcpy(_init_type->option_input_str, option_def[2]); map_add_options_input_text(_init_type);  
  strcpy(_init_type->option_input_str, option_def[3]); map_add_options_input_text(_init_type);  
  strcpy(_init_type->option_input_str, option_def[4]); map_add_options_input_text(_init_type);  
  strcpy(_init_type->option_input_str, option_def[5]); map_add_options_input_text(_init_type);  
  //strcpy(_init_type->option_input_str, option_def[6]); map_add_options_input_text(_init_type);   
   
  strcpy(_init_type->summary_file_name,"baseline.sum.map\0"); map_set_summary_file_name(_init_type, map_msg, &ierr);
   
   
   
  /////////// CALL MAP_INIT //////////////////  
  //WARNING: adrien@marco : map_init can fail without raising warning or exception and program continues: is it expected ?
  map_init(_init_type, _u_type, _p_type, _x_type, NULL, _z_type, _other_type, _y_type, _io_type, &ierr, map_msg);
  if (ierr!=MAP_SAFE) {
    std::cout << map_msg << std::endl;
  };
 
  
   //////////// DEALLOCATE SOME MEMORY NEEDED ONLY FOR INIT HERE //////////////
   MAPFREE(_init_type); 
   MAPFREE(_io_type);
   
   cout << "--------- End of INITIALIZATION ----------" << endl;
   
   
   ///////// TIME SOLVER EMULATOR (forced displacements along X) /////////////
   
   cout << "--------- Start of SIMULATION ----------" << endl;
   
   int ntimesteps = 10;
   double dx = 1.0;
   Domain* data = (Domain*)_other_type->object;
   Vessel* vessel = &data->vessel;
   int N = _u_type->x_Len;
   
   
   for (int k = 0; k<ntimesteps; ++k) {     
     // map_reset_universal_error(map_msg, &ierr); /* should not need to do this. */     
     for (int i=0 ; i<N ; i++) {        
       _u_type->x[i] = k*dx;
       _u_type->y[i] = 0.0;
       _u_type->z[i] = 0.0;       
     };
     
     //NOTE: adrien @ marco : I put here only map_update_states and get rid of what seems to be unused input
     map_update_states(0, 0, _u_type, _p_type, _x_type, NULL, _z_type, _other_type, &ierr, map_msg);
     
     cout << "Time step: " << k << endl;
     cout << "Element 0 fairlead force: " << _y_type->Fx[0] << " " << _y_type->Fy[0] << " " << _y_type->Fz[0] << endl; 
   };   
   
   cout << "--------- End of SIMULATION ----------" << endl;
   
   //////////////////// MEMORY DEALLOCATION /////////////////
   cout << "--------- Start of DEALLOCATION ----------" << endl;   

   // MAP_ERROR_CODE success = MAP_SAFE;
   map_end(_u_type, _p_type, _x_type, NULL, _z_type, _other_type, _y_type, &ierr, map_msg);
   success = map_free_types(_u_type, _p_type, _x_type, _z_type, _other_type, _y_type); 
   MAPFREE(_other_type); 
   MAPFREE(_y_type); 
   MAPFREE(_u_type);
   MAPFREE(_p_type);
   MAPFREE(_z_type);
   MAPFREE(_x_type);
   
  cout << "--------- End of DEALLOCATION ----------" << endl;
  
  /////////////////////////////////////////////////////////////
  
  cout << "C++ driver test ends" << endl;
  
  return(0);
    
}
