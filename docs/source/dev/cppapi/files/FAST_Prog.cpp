#include "OpenFAST.H"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <mpi.h>

void readTurbineData(int iTurb, fast::fastInputs & fi, YAML::Node turbNode) {
  //Read turbine data for a given turbine using the YAML node
}

void readInputFile(fast::fastInputs & fi, std::string cInterfaceInputFile, double * tEnd) {
    //Read input data for a given turbine using the YAML node
}

int main() {
  int iErr;
  int nProcs;
  int rank;

  iErr = MPI_Init(NULL, NULL);
  iErr = MPI_Comm_size( MPI_COMM_WORLD, &nProcs);
  iErr = MPI_Comm_rank( MPI_COMM_WORLD, &rank);

  double tEnd ; // This doesn't belong in the OpenFAST - C++ API
  int ntEnd ; // This doesn't belong in the OpenFAST - C++ API

  std::string cDriverInputFile="cDriver.i";
  fast::OpenFAST FAST;
  fast::fastInputs fi ;
  readInputFile(fi, cDriverInputFile, &tEnd);
  ntEnd = tEnd/fi.dtFAST;  //Calculate the last time step

  FAST.setInputs(fi);
  // In a parallel simulation, multiple turbines have to be allocated to processors.
  // The C++ API can handle any allocation of turbines on an arbitrary number of processors
  FAST.allocateTurbinesToProcsSimple(); // Use this for a simple round robin allocation of turbines to processors.
  // Or allocate turbines to procs by calling "setTurbineProcNo(iTurbGlob, procId)" for each turbine.

  FAST.init();
  if (FAST.isTimeZero()) {
    FAST.solution0();
  }

  if( !FAST.isDryRun() ) {
    for (int nt = FAST.get_ntStart(); nt < ntEnd; nt++) {
      FAST.step();
    }
  }

  FAST.end() ;
  MPI_Finalize() ;

  return 0;
    
}
