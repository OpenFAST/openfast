#include "OpenFAST.H"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <mpi.h>

inline bool checkFileExists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

void readTurbineData(int iTurb, fast::fastInputs & fi, YAML::Node turbNode) {

  //Read turbine data for a given turbine using the YAML node
  if (turbNode["turb_id"]) {
      fi.globTurbineData[iTurb].TurbID = turbNode["turb_id"].as<int>();
  } else {
      turbNode["turb_id"] = iTurb;
  }
  if (turbNode["FAST_input_filename"]) {
      fi.globTurbineData[iTurb].FASTInputFileName = turbNode["FAST_input_filename"].as<std::string>() ;  
  } else {
      fi.globTurbineData[iTurb].FASTInputFileName = "";
  }
  if (turbNode["restart_filename"]) {
      fi.globTurbineData[iTurb].FASTRestartFileName = turbNode["restart_filename"].as<std::string>() ;
  } else {
      fi.globTurbineData[iTurb].FASTRestartFileName = "";
  }
  if ( (fi.globTurbineData[iTurb].FASTRestartFileName == "") && (fi.globTurbineData[iTurb].FASTInputFileName == "") )
      throw std::runtime_error("Both FAST_input_filename and restart_filename are empty or not specified for Turbine " + std::to_string(iTurb));
  if (turbNode["turbine_base_pos"].IsSequence() ) {
      fi.globTurbineData[iTurb].TurbineBasePos = turbNode["turbine_base_pos"].as<std::vector<float> >() ;
  } else {
      fi.globTurbineData[iTurb].TurbineBasePos = std::vector<float>(3,0.0);
  }
  if (turbNode["turbine_hub_pos"].IsSequence() ) {
      fi.globTurbineData[iTurb].TurbineHubPos = turbNode["turbine_hub_pos"].as<std::vector<double> >() ;
  } else {
      fi.globTurbineData[iTurb].TurbineHubPos =  std::vector<double>(3,0.0);
  }
  if (turbNode["num_force_pts_blade"]) {
      fi.globTurbineData[iTurb].numForcePtsBlade = turbNode["num_force_pts_blade"].as<int>();
  } else {
      fi.globTurbineData[iTurb].numForcePtsBlade = 0;
  }
  if (turbNode["num_force_pts_tower"]) {
      fi.globTurbineData[iTurb].numForcePtsTwr = turbNode["num_force_pts_tower"].as<int>();
  } else {
      fi.globTurbineData[iTurb].numForcePtsTwr = 0;
  }
  if (turbNode["nacelle_cd"]) {
      fi.globTurbineData[iTurb].nacelle_cd = turbNode["nacelle_cd"].as<float>();
  } else {
      fi.globTurbineData[iTurb].nacelle_cd = 0.0;
  }
  if (turbNode["nacelle_area"]) {
      fi.globTurbineData[iTurb].nacelle_area = turbNode["nacelle_area"].as<float>();
  } else {
      fi.globTurbineData[iTurb].nacelle_area = 0.0;
  }
  if (turbNode["air_density"]) {
      fi.globTurbineData[iTurb].air_density = turbNode["air_density"].as<float>();
  } else {
      fi.globTurbineData[iTurb].air_density = 0.0;
  }
}

void readInputFile(fast::fastInputs & fi, std::string cInterfaceInputFile, double * tEnd, int * couplingMode, bool * setExpLawWind) {

  fi.comm = MPI_COMM_WORLD;

  // Check if the input file exists and read it
  if ( checkFileExists(cInterfaceInputFile) ) {

    YAML::Node cDriverInp = YAML::LoadFile(cInterfaceInputFile);

    fi.nTurbinesGlob = cDriverInp["n_turbines_glob"].as<int>();

    if (fi.nTurbinesGlob > 0) {
      
      if(cDriverInp["dry_run"]) {
	fi.dryRun = cDriverInp["dry_run"].as<bool>();
      } 
      
      if(cDriverInp["debug"]) {
	fi.debug = cDriverInp["debug"].as<bool>();
      }

      *couplingMode = 0; //CLASSIC is default
      if(cDriverInp["coupling_mode"]) {
          if ( cDriverInp["coupling_mode"].as<std::string>() == "strong" ) {
              *couplingMode = 1;
          } else if ( cDriverInp["coupling_mode"].as<std::string>() == "classic" ) {
              *couplingMode = 0;
          } else {
              throw std::runtime_error("coupling_mode is not well defined in the input file");
          }
      }
      
      if(cDriverInp["sim_start"]) {
	if (cDriverInp["sim_start"].as<std::string>() == "init") {
	  fi.simStart = fast::init;
	} else if(cDriverInp["sim_start"].as<std::string>() == "trueRestart") {
	  fi.simStart = fast::trueRestart;
	} else if(cDriverInp["sim_start"].as<std::string>() == "restartDriverInitFAST") {
	  fi.simStart = fast::restartDriverInitFAST;
	} else {
	  throw std::runtime_error("sim_start is not well defined in the input file");
	}
      }

      if(cDriverInp["t_start"]) { 
          fi.tStart = cDriverInp["t_start"].as<double>();
      } else {
          throw std::runtime_error("t_start is missing in the input file");
      }
      if(cDriverInp["t_end"]) {
          *tEnd = cDriverInp["t_end"].as<double>();
      } else {
          throw std::runtime_error("t_end is missing in the input file");
      }
      if(cDriverInp["n_checkpoint"]) {
          fi.nEveryCheckPoint = cDriverInp["n_checkpoint"].as<int>();
      } else {
          throw std::runtime_error("n_checkpoint is missing in the input file");
      }
      if (cDriverInp["dt_FAST"]) {
          fi.dtFAST = cDriverInp["dt_FAST"].as<double>();
      } else {
          throw std::runtime_error("dt_FAST is missing in the input file");
      }
      if (cDriverInp["n_substeps"]) {
          fi.nSubsteps = cDriverInp["n_substeps"].as<int>();
      } else {
          std::cout << "n_substeps is missing in the input file. Proceeding with assumption n_substeps=1." << std::endl ;
          fi.nSubsteps = 1;
      }
      if (cDriverInp["t_max"]) {
          fi.tMax = cDriverInp["t_max"].as<double>(); // t_max is the total duration to which you want to run FAST. This should be the same or greater than the max time given in the FAST fst file. Choose this carefully as FAST writes the output file only at this point if you choose the binary file output.
      } else {
          fi.tMax = 0;
      }

      if(cDriverInp["set_exp_law_wind"]) {
        *setExpLawWind = cDriverInp["set_exp_law_wind"].as<bool>();
      } else {
          *setExpLawWind = false;
      }
      
      if(cDriverInp["super_controller"]) {
	fi.scStatus = cDriverInp["super_controller"].as<bool>();
	fi.scLibFile = cDriverInp["sc_libfile"].as<std::string>();
	fi.numScInputs = cDriverInp["num_scinputs"].as<int>();
	fi.numScOutputs = cDriverInp["num_scoutputs"].as<int>();
      }
      
      fi.globTurbineData.resize(fi.nTurbinesGlob);
      for (int iTurb=0; iTurb < fi.nTurbinesGlob; iTurb++) {
	if (cDriverInp["Turbine" + std::to_string(iTurb)]) {
	  readTurbineData(iTurb, fi, cDriverInp["Turbine" + std::to_string(iTurb)] );
	} else {
	  throw std::runtime_error("Node for Turbine" + std::to_string(iTurb) + " not present in input file or I cannot read it");
	}
      }
      
    } else {
      throw std::runtime_error("Number of turbines <= 0 ");
    }
    
  } else {
    throw std::runtime_error("Input file " + cInterfaceInputFile + " does not exist or I cannot access it");
  }

}

int main() {
  int iErr;
  int nProcs;
  int rank;
  std::vector<double> torque (3, 0.0);
  std::vector<double> thrust (3, 0.0);  

  iErr = MPI_Init(NULL, NULL);
  iErr = MPI_Comm_size( MPI_COMM_WORLD, &nProcs);
  iErr = MPI_Comm_rank( MPI_COMM_WORLD, &rank);

  int couplingMode ; //CLASSIC (SOWFA style = 0) or STRONG (Conventional Serial Staggered - allow for outer iterations = 1)
  double tEnd ; // This doesn't belong in the FAST - C++ interface 
  int ntStart, ntEnd ; // This doesn't belong in the FAST - C++ interface
  bool setExpLawWind; // Set wind speed at Aerodyn nodes based on an exponential profile. Useful for testing the C++ API before running actuator line simulations.

  std::string cDriverInputFile="cDriver.i";
  fast::OpenFAST FAST;
  fast::fastInputs fi ;
  try {
      readInputFile(fi, cDriverInputFile, &tEnd, &couplingMode, &setExpLawWind);
  }
  catch( const std::runtime_error & ex) {
    std::cerr << ex.what() << std::endl ;
    std::cerr << "Program quitting now" << std::endl ;
    return 1;
  }

  FAST.setInputs(fi);
  FAST.allocateTurbinesToProcsSimple(); 
  // Or allocate turbines to procs by calling "setTurbineProcNo(iTurbGlob, procId)" for turbine.

  FAST.init();

  ntStart = fi.tStart/fi.dtFAST/fi.nSubsteps;  //Calculate the first time step
  ntEnd = tEnd/fi.dtFAST/fi.nSubsteps;  //Calculate the last time step

  if (setExpLawWind) {
      FAST.setExpLawWindSpeed();   
  }
  
  if (FAST.isTimeZero()) {
    FAST.solution0();
  }
  
  if( !FAST.isDryRun() ) {
    for (int nt = ntStart; nt < ntEnd; nt++) {
        if (couplingMode == 0) {
            // If running with a CFD solver, sample velocities at the actuator/velocity nodes here
            if (setExpLawWind)
                FAST.setExpLawWindSpeed();   
            for (int iSubstep=1; iSubstep < fi.nSubsteps+1; iSubstep++)
                FAST.step();
            // Get forces at actuator nodes and advance CFD solve by one time step here
        } else {
            for (int j=0; j < 2; j++) {
                // If running with a CFD solver, use 'FAST.predict_states()' to predict position and force at actuator nodes at the next time step on the first pass
                // Run a CFD time step as a 'predictor' to get velocity at the next time step
                // Sample and set velocity at the actuator/velocity nodes after the first cfd predictor
                if (setExpLawWind)
                    FAST.setExpLawWindSpeed();   
                FAST.update_states_driver_time_step();
            }
            // Call this after enough outer iterations have been done
            FAST.advance_to_next_driver_time_step();
        }
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

