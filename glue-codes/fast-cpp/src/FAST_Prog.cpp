#include "OpenFAST.H"
#include <iostream>
#include <mpi.h>

int main() {
  int iErr;
  int nProcs;
  int rank;
  double torque[] = {0.0,0.0,0.0};
  double thrust[] = {0.0,0.0,0.0};  

  iErr = MPI_Init(NULL, NULL);
  iErr = MPI_Comm_size( MPI_COMM_WORLD, &nProcs);
  iErr = MPI_Comm_rank( MPI_COMM_WORLD, &rank);

  std::string cDriverInputFile="cDriver.i";
  fast::OpenFAST FAST;
  try {
    FAST.readInputFile(cDriverInputFile);
  }
  catch( const std::runtime_error & ex) {
    std::cerr << ex.what() << std::endl ;
    std::cerr << "Program quitting now" << std::endl ;
    return 1;
  }

  FAST.allocateTurbinesToProcsSimple(); 
  // Or allocate turbines to procs by calling "setTurbineProcNo(iTurbGlob, procId)" for turbine.
  FAST.allocateInputData(); // Allocate memory for all inputs that are dependent on the number of turbines

  if ( !FAST.isDryRun() ) {
    FAST.init();
    // Set velocity at the aerodyn nodes here if necessary
    if (FAST.isTimeZero()) {
      FAST.solution0();
    }
  }

  if( !FAST.isDryRun() ) {
    for (int nt = FAST.get_ntStart(); nt <= FAST.get_ntEnd(); nt++) {
      FAST.step();
      if (FAST.isDebug()) {
	FAST.computeTorqueThrust(0,torque,thrust);
	std::cout.precision(16);
	std::cout << "Torque = " << torque[0] << " " << torque[1] << " " << torque[2] << std::endl ;
	std::cout << "Thrust = " << thrust[0] << " " << thrust[1] << " " << thrust[2] << std::endl ;      
      }
    }
  }

  FAST.end() ;
  MPI_Finalize() ;

  return 0;
    
}

