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
    fi.globTurbineData[iTurb].TurbID = turbNode["turb_id"].as<int>();
    fi.globTurbineData[iTurb].FASTInputFileName = turbNode["FAST_input_filename"].as<std::string>();
    fi.globTurbineData[iTurb].FASTRestartFileName = turbNode["restart_filename"].as<std::string>();
    if (turbNode["turbine_base_pos"].IsSequence() ) {
        fi.globTurbineData[iTurb].TurbineBasePos = turbNode["turbine_base_pos"].as<std::vector<double> >();
    }
    if (turbNode["turbine_hub_pos"].IsSequence() ) {
        fi.globTurbineData[iTurb].TurbineHubPos = turbNode["turbine_hub_pos"].as<std::vector<double> >();
    }
    fi.globTurbineData[iTurb].numForcePtsBlade = turbNode["num_force_pts_blade"].as<int>();
    fi.globTurbineData[iTurb].numForcePtsTwr = turbNode["num_force_pts_tower"].as<int>();
    if (turbNode["nacelle_cd"]) fi.globTurbineData[iTurb].nacelle_cd = turbNode["nacelle_cd"].as<float>();
    if (turbNode["nacelle_area"]) fi.globTurbineData[iTurb].nacelle_area = turbNode["nacelle_area"].as<float>();
    if (turbNode["air_density"]) fi.globTurbineData[iTurb].air_density = turbNode["air_density"].as<float>();
}

void readInputFile(fast::fastInputs & fi, std::string cInterfaceInputFile, double * tEnd) {

    fi.comm = MPI_COMM_WORLD;

    // Check if the input file exists and read it
    if ( checkFileExists(cInterfaceInputFile) ) {

        YAML::Node cDriverInp = YAML::LoadFile(cInterfaceInputFile);

        fi.nTurbinesGlob = cDriverInp["nTurbinesGlob"].as<int>();

        if (fi.nTurbinesGlob > 0) {

            if(cDriverInp["dryRun"]) {
                fi.dryRun = cDriverInp["dryRun"].as<bool>();
            } 

            if(cDriverInp["debug"]) {
                fi.debug = cDriverInp["debug"].as<bool>();
            } 

            if(cDriverInp["simStart"]) {
                if (cDriverInp["simStart"].as<std::string>() == "init") {
                    fi.simStart = fast::init;
                } else if(cDriverInp["simStart"].as<std::string>() == "trueRestart") {
                    fi.simStart = fast::trueRestart;
                } else if(cDriverInp["simStart"].as<std::string>() == "restartDriverInitFAST") {
                    fi.simStart = fast::restartDriverInitFAST;
                } else {
                    throw std::runtime_error("simStart is not well defined in the input file");
                }
            }

            fi.tStart = cDriverInp["tStart"].as<double>();
            *tEnd = cDriverInp["tEnd"].as<double>();
            fi.nEveryCheckPoint = cDriverInp["nEveryCheckPoint"].as<int>();
            fi.dtFAST = cDriverInp["dtFAST"].as<double>();
            fi.tMax = cDriverInp["tMax"].as<double>(); // tMax is the total duration to which you want to run FAST. This should be the same or greater than the max time given in the FAST fst file. Choose this carefully as FAST writes the output file only at this point if you choose the binary file output.

            if(cDriverInp["superController"]) {
                fi.scStatus = cDriverInp["superController"].as<bool>();
                fi.scLibFile = cDriverInp["scLibFile"].as<std::string>();
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

    double tEnd ; // This doesn't belong in the FAST - C++ interface 
    int ntEnd ; // This doesn't belong in the FAST - C++ interface

    std::string cDriverInputFile=argv[1];
    fast::OpenFAST FAST;
    fast::fastInputs fi ;
    try {
        readInputFile(fi, cDriverInputFile, &tEnd);
    } catch( const std::runtime_error & ex) {
        std::cerr << ex.what() << std::endl ;
        std::cerr << "Program quitting now" << std::endl ;
        return 1;
    }

    // Calculate the last time step
    ntEnd = tEnd/fi.dtFAST;

    FAST.setInputs(fi);
    FAST.allocateTurbinesToProcsSimple(); 
    // Or allocate turbines to procs by calling "setTurbineProcNo(iTurbGlob, procId)" for turbine.

    FAST.init();

    if (FAST.isTimeZero()) FAST.solution0();

    if ( FAST.isDryRun() ) {
        FAST.end() ;
        MPI_Finalize() ;
        return 0;
    }

    for (int nt = FAST.get_ntStart(); nt < ntEnd; nt++) {
        FAST.step();
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
