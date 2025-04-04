#include "OpenFAST.H"
#include "yaml-cpp/yaml.h"
#include <iostream>
#include <cmath>
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
  std::string simType;
  get_if_present(turbNode, "sim_type", simType, std::string("ext-inflow"));
  if (simType == "ext-loads")
      fi.globTurbineData[iTurb].sType = fast::EXTLOADS;
  else
      fi.globTurbineData[iTurb].sType = fast::EXTINFLOW;

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
  fi.globTurbineData[iTurb].numForcePts =
      fi.globTurbineData[iTurb].numForcePtsBlade +
      fi.globTurbineData[iTurb].numForcePtsTwr;

  float fZero = 0.0;
  get_if_present(turbNode, "nacelle_cd", fi.globTurbineData[iTurb].nacelle_cd, fZero);
  get_if_present(turbNode, "nacelle_area", fi.globTurbineData[iTurb].nacelle_area, fZero);
  get_if_present(turbNode, "air_density", fi.globTurbineData[iTurb].air_density, fZero);

  if (simType == "ext-loads") {

      get_if_present(turbNode, "az_blend_mean", fi.globTurbineData[iTurb].azBlendMean, 20*360.0*M_PI/180.0); //20 revs
      get_if_present(turbNode, "az_blend_delta", fi.globTurbineData[iTurb].azBlendDelta, 3.0*360.0*M_PI/180.0);  // 3 rev

  }

}

void readInputFile(fast::fastInputs & fi, std::string cInterfaceInputFile, double *tStart, double * tEnd, int * couplingMode, bool * setExpLawWind, bool * setUniformXBladeForces, int * nIter, double *xBladeForce) {

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
        if (cDriverInp["n_iter"]) {
            *nIter = cDriverInp["n_iter"].as<int>();
            if (*nIter < 0) {
                *nIter = 1;
            }
        } else {
            *nIter = 1;
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

        get_required(cDriverInp, "t_start", *tStart);
        get_required(cDriverInp, "t_end", *tEnd);
        get_required(cDriverInp, "restart_freq", fi.restartFreq);
        get_if_present(cDriverInp, "output_freq", fi.outputFreq, 100);
        get_required(cDriverInp, "dt_driver", fi.dtDriver);
        get_required(cDriverInp, "t_max", fi.tMax); // t_max is the total duration to which you want to run FAST. This should be the same or greater than the max time given in the FAST fst file.
        get_if_present(cDriverInp, "set_exp_law_wind", *setExpLawWind, false);
        get_if_present(cDriverInp, "set_uniform_x_blade_forces", *setUniformXBladeForces, false);
        if (setUniformXBladeForces)
            get_if_present(cDriverInp, "x_blade_force", *xBladeForce, 0.0);

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

int main(int argc, char** argv) {

    if (argc != 2) {
        std::cerr << "Incorrect syntax. Try: openfastcpp inputfile.yaml" << std::endl ;
        return 1;
    }

    int iErr;
    int nProcs;
    int rank;
    std::vector<double> torque (3, 0.0);
    std::vector<double> thrust (3, 0.0);

    iErr = MPI_Init(NULL, NULL);
    iErr = MPI_Comm_size( MPI_COMM_WORLD, &nProcs);
    iErr = MPI_Comm_rank( MPI_COMM_WORLD, &rank);

    int couplingMode ; //CLASSIC (SOWFA style = 0) or STRONG (Conventional Serial Staggered - allow for outer iterations = 1)
    double tStart; // This doesn't belong in the C++ API
    double tEnd ; // This doesn't belong in the FAST - C++ interface
    int ntStart, ntEnd ; // This doesn't belong in the FAST - C++ interface
    int nSubsteps; //
    bool setExpLawWind; // Set wind speed at Aerodyn nodes based on an exponential profile. Useful for testing the C++ API before running actuator line simulations.
    bool setUniformXBladeForces; // Set uniform X blade forces on all blade nodes
    int nIter;
    double xBladeForce = 0.0;

    std::string cDriverInputFile=argv[1];
    fast::OpenFAST FAST;
    fast::fastInputs fi ;
    try {
        readInputFile(fi, cDriverInputFile, &tStart, &tEnd, &couplingMode, &setExpLawWind, &setUniformXBladeForces, &nIter, &xBladeForce);
    } catch( const std::runtime_error & ex) {
        std::cerr << ex.what() << std::endl ;
        std::cerr << "Program quitting now" << std::endl ;
        return 1;
    }

    FAST.setInputs(fi);
    FAST.allocateTurbinesToProcsSimple();
    // Or allocate turbines to procs by calling "setTurbineProcNo(iTurbGlob, procId)" for turbine.

    FAST.init();

    nSubsteps = fi.dtDriver/FAST.get_timestep();

    if ( FAST.isDryRun() ) {
        FAST.end() ;
        MPI_Finalize() ;
        return 0;
    }

    if (FAST.isTimeZero()) {
        if (setExpLawWind)
            FAST.setExpLawWindSpeed(0.0);

        FAST.solution0();
    }


    ntStart = tStart/fi.dtDriver;  //Calculate the first time step
    ntEnd = tEnd/fi.dtDriver;  //Calculate the last time step

    for (int nt = ntStart; nt < ntEnd; nt++) {
        if (couplingMode == 0) {
            // If running with a CFD solver, sample velocities at the actuator/velocity nodes here
            if (setExpLawWind)
                FAST.setExpLawWindSpeed( (nt+1)*fi.dtDriver );
            if (setUniformXBladeForces) {
                FAST.setUniformXBladeForces(xBladeForce);
            }

            for (int iSubstep=1; iSubstep < nSubsteps; iSubstep++) {
                FAST.step();
                std::cout << "iSubstep = " << iSubstep << std::endl ;
            }
            // Get forces at actuator nodes and advance CFD solve by one time step here
        } else {
            for (int j=0; j < nIter; j++) {
                // If running with a CFD solver, use 'FAST.predict_states()' to predict position and force at actuator nodes at the next time step on the first pass
                // Run a CFD time step as a 'predictor' to get velocity at the next time step
                // Sample and set velocity at the actuator/velocity nodes after the first cfd predictor
                if (setExpLawWind)
                    FAST.setExpLawWindSpeed( (nt+1)*fi.dtDriver );
                if (setUniformXBladeForces) {
                    FAST.setUniformXBladeForces(xBladeForce);
                }
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

    FAST.end() ;
    MPI_Finalize() ;

    return 0;
}
