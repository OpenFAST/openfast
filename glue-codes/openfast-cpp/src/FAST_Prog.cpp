#include "OpenFAST.H"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <mpi.h>

inline bool checkFileExists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

/// Optionally read in a value from a yaml node if present, else set it to a default value. Copied from github.com/Exawind/nalu-wind/include/NaluParsing.h
template<typename T>
void get_if_present(const YAML::Node & node, const std::string& key, T& result, const T& default_if_not_present = T())
{
    if (node[key]) {
        const YAML::Node value = node[key];
        result = value.as<T>();
    }
    else {
        int rank;
        int iErr = MPI_Comm_rank( MPI_COMM_WORLD, &rank);
        if(!rank)
            std::cout << key << " is missing in the input file. Proceeding with assumption " << key << " = " << default_if_not_present << std::endl ;
        result = default_if_not_present;
    }
}

/// Read a 'key' from a yaml node if it exists, else throw an error
template<typename T>
void get_required(const YAML::Node & node, const std::string& key, T& result)
{
    if (node[key]) {
        const YAML::Node value = node[key];
        result = value.as<T>();
    }
    else    {
        throw std::runtime_error("Error: parsing missing required key: " + key);
    }
}

void readTurbineData(int iTurb, fast::fastInputs & fi, YAML::Node turbNode) {

  //Read turbine data for a given turbine using the YAML node
  get_if_present(turbNode, "turb_id", fi.globTurbineData[iTurb].TurbID, iTurb);
  std::string emptyString = "";  
  get_if_present(turbNode, "FAST_input_filename", fi.globTurbineData[iTurb].FASTInputFileName);
  get_if_present(turbNode, "restart_filename", fi.globTurbineData[iTurb].FASTRestartFileName);
  if ( (fi.globTurbineData[iTurb].FASTRestartFileName == emptyString) && (fi.globTurbineData[iTurb].FASTInputFileName == emptyString) )
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
  get_if_present(turbNode, "num_force_pts_blade", fi.globTurbineData[iTurb].numForcePtsBlade, 0);
  get_if_present(turbNode, "num_force_pts_tower", fi.globTurbineData[iTurb].numForcePtsTwr, 0);
  float fZero = 0.0;
  get_if_present(turbNode, "nacelle_cd", fi.globTurbineData[iTurb].nacelle_cd, fZero);
  get_if_present(turbNode, "nacelle_area", fi.globTurbineData[iTurb].nacelle_area, fZero);
  get_if_present(turbNode, "air_density", fi.globTurbineData[iTurb].air_density, fZero);
}

void readInputFile(fast::fastInputs & fi, std::string cInterfaceInputFile, double * tEnd, int * couplingMode, bool * setExpLawWind) {

  fi.comm = MPI_COMM_WORLD;

  // Check if the input file exists and read it
  if ( checkFileExists(cInterfaceInputFile) ) {

    YAML::Node cDriverInp = YAML::LoadFile(cInterfaceInputFile);

    get_required(cDriverInp, "n_turbines_glob", fi.nTurbinesGlob);

    if (fi.nTurbinesGlob > 0) {
      
        get_if_present(cDriverInp, "dry_run", fi.dryRun, false);
        get_if_present(cDriverInp, "debug", fi.debug, false);

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
                fi.simStart = fast::INIT;
            } else if(cDriverInp["sim_start"].as<std::string>() == "trueRestart") {
                fi.simStart = fast::TRUERESTART;
            } else if(cDriverInp["sim_start"].as<std::string>() == "restartDriverInitFAST") {
                fi.simStart = fast::RESTARTDRIVERINITFAST;
            } else {
                throw std::runtime_error("sim_start is not well defined in the input file");
            }
        }
        
        get_required(cDriverInp, "t_start", fi.tStart);
        get_required(cDriverInp, "t_end", *tEnd);
        get_required(cDriverInp, "n_checkpoint", fi.nEveryCheckPoint);
        get_required(cDriverInp, "dt_FAST", fi.dtFAST);
        get_if_present(cDriverInp, "n_substeps", fi.nSubsteps, 1);
        get_required(cDriverInp, "t_max", fi.tMax); // t_max is the total duration to which you want to run FAST. This should be the same or greater than the max time given in the FAST fst file. 
        get_if_present(cDriverInp, "set_exp_law_wind", *setExpLawWind, false);
        
        get_if_present(cDriverInp, "super_controller", fi.scStatus, false);
        if(fi.scStatus) {
            get_required(cDriverInp, "sc_libfile", fi.scLibFile);
            get_required(cDriverInp, "num_scinputs", fi.numScInputs);
            get_required(cDriverInp, "num_scoutputs", fi.numScOutputs);
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
  readInputFile(fi, cDriverInputFile, &tEnd, &couplingMode, &setExpLawWind);

  FAST.setInputs(fi);
  FAST.allocateTurbinesToProcsSimple(); 
  // Or allocate turbines to procs by calling "setTurbineProcNo(iTurbGlob, procId)" for turbine.

  FAST.init();

  ntStart = fi.tStart/fi.dtFAST/fi.nSubsteps;  //Calculate the first time step
  ntEnd = tEnd/fi.dtFAST/fi.nSubsteps;  //Calculate the last time step

  if (setExpLawWind)
      FAST.setExpLawWindSpeed();   
  
  if (FAST.isTimeZero())
    FAST.solution0();
  
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

