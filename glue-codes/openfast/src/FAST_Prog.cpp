
#include "stdio.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include "FastLibAPI.h"

int main(int argc, char** argv) {
   if (argc != 2) {
      std::cerr << "Incorrect syntax. Expected syntax is `openfast_cpp input.fst`" << std::endl;
      return 1;
   }

   std::string input_file_name = argv[1];

   FastLibAPI fastlib = FastLibAPI(input_file_name);
   fastlib.fast_run();

   // // Get the hub position 
   // float absolute_position[3] = {};
   // float rotational_velocity[3] = {};
   // double orientation_dcm[9] = {};

   // fastlib.get_hub_position(absolute_position, rotational_velocity, orientation_dcm);

   // printf("%f %f %f\n", absolute_position[0], absolute_position[1], absolute_position[2]);
   // printf("%f %f %f\n", rotational_velocity[0], rotational_velocity[1], rotational_velocity[2]);
   // printf("%f %f %f\n", orientation_dcm[0], orientation_dcm[1], orientation_dcm[2]);
   // printf("%f %f %f\n", orientation_dcm[3], orientation_dcm[4], orientation_dcm[5]);
   // printf("%f %f %f\n", orientation_dcm[6], orientation_dcm[7], orientation_dcm[8]);

   return 0;
}
