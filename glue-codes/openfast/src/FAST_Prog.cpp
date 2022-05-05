
#include "stdio.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include "FastLibAPI.h"

int main(int argc, char** argv) {
   if (argc != 2) {
      std::cerr << "Incorrect syntax. Try: openfastcpp inputfile.yaml" << std::endl ;
      return 1;
   }

   std::string input_file_name = argv[1];

   FastLibAPI fastlib = FastLibAPI(input_file_name);
   fastlib.fast_run();
}
