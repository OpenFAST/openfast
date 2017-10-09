#ifndef OpenFAST_h
#define OpenFAST_h
#include "FAST_Library.h"
#include "sys/stat.h"
#include "math.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include "dlfcn.h"
#ifndef OMPI_SKIP_MPICXX
 #define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
 #define MPICH_SKIP_MPICXX
#endif
#include "mpi.h"
#include "SC.h"


namespace fast {

/**
 * A struct to hold all input data specific to each turbine
 */
struct globTurbineDataType {
  //!Integer id for every turbine
  int TurbID;
  //! The FAST Input file name. Typically a .fst file.
  std::string FASTInputFileName;
  //! The restart/checkpoint file name. 
  std::string FASTRestartFileName;
  //! The position of the base of the turbine in global coordinates
  std::vector<double> TurbineBasePos;
  //! The approximate position of the hub of the turbine in global coordinates
  std::vector<double> TurbineHubPos;
  //! Desired number of actuator points on each blade
  int numForcePtsBlade;
  //! Desired number of actuator points on the tower
  int numForcePtsTwr;
};

//! An id to indicate whether a particular actuator point is on the hub, node or tower    
enum ActuatorNodeType {
  HUB = 0,
  BLADE = 1,
  TOWER = 2,
  ActuatorNodeType_END
};

/** An id to indicate the start type of a simulation.
 * init - Start the simulation from scratch
 * trueRestart - Restart from a checkpoint file. Code expects checkpoint files for all parts of the simulation including the controller.
 * restartDriverInitFAST - Start all turbines from scratch and use the velocity data in 'velData.h5' file to run upto desired restart time, then continue the simulation like ''trueRestart'.
 */
enum simStartType {
  init = 0,
  trueRestart = 1,
  restartDriverInitFAST = 2,
  simStartType_END    
};

/**
 * A class to hold all input data for a simulation run through a OpenFAST C++ glue code
 */
class fastInputs {

 public: 

  //! MPI Communicator
  MPI_Comm comm;
  //! Total number of turbines in the simulation
  int nTurbinesGlob;
  //! The simulation will not run if dryRun is set to true. However, the simulation will read the input files, allocate turbines to processors and prepare to run the individual turbine instances. This flag is useful to test the setup of the simulation before running it.
  bool dryRun;
  //! Enable debug outputs if set to true
  bool debug;
  //! Start time of the simulation
  double tStart;
    //! Start type of the simulation: 'init', 'trueRestart' or 'restartDriverInitFAST'.
  simStartType simStart;
    //!Restart files will be written every so many time stneps
  int nEveryCheckPoint;
    //! Max time of the simulation
  double tMax;
    //! Time step for FAST. All turbines should have the same time step.
  double dtFAST;  

    //! Supercontroller status: True/False.
  bool scStatus;
    //! Name of the dynamic library containing the supercontroller implementation
  std::string scLibFile;
    //! Number of inputs and output to the supercontroller from/to each turbine
  int numScInputs, numScOutputs;

    //! Vector of turbine specific input data
  std::vector<globTurbineDataType>  globTurbineData;

  // Constructor 
  fastInputs() ;
 
  // Destructor
  ~fastInputs() {} ;

};


/**
 * A class to interface OpenFAST's fortran backend with a C++ driver program
 */    
class OpenFAST {

 private:

  //! MPI Communicator
  MPI_Comm mpiComm;
  //! The simulation will not run if dryRun is set to true. However, the simulation will read the input files, allocate turbines to processors and prepare to run the individual turbine instances. This flag is useful to test the setup of the simulation before running it.
  bool dryRun;        // If this is true, class will simply go through allocation and deallocation of turbine data
  //! Enable debug outputs if set to true
  bool debug;   // Write out extra information if this flags is turned on
  //! Vector of turbine specific input data
  std::vector<globTurbineDataType> globTurbineData;
  //! Number of turbines on this MPI rank
  int nTurbinesProc;
  //! Total number of turbines in the simulation
  int nTurbinesGlob;
  //! Start type of the simulation: 'init', 'trueRestart' or 'restartDriverInitFAST'.
  simStartType simStart;
  //! Is the time now zero: True/False
  bool timeZero;
  //! Time step for FAST. All turbines should have the same time step.
  double dtFAST;
  //! Max time of the simulation
  double tMax;
  //! The position of the base of the turbine in global coordinates for each turbine
  std::vector<std::vector<float> > TurbineBasePos;
  //! The approximate position of the hub of the turbine in global coordinates for each turbine
  std::vector<std::vector<float> > TurbineHubPos;
  //!Integer id for every turbine
  std::vector<int> TurbID;
  //! The FAST Input file name for every turbine. Typically a .fst file.
  std::vector<std::string > FASTInputFileName;
  //! The prefix for the checkpoint file name for each turbine. 
  std::vector<std::string > CheckpointFileRoot;
  //! Start time of the simulation
  double tStart;
  //! The current time step number 
  int nt_global;
  //! The starting time step number
  int ntStart;      // The time step to start the FAST simulation
  //! Restart files will be written every so many time stneps
  int nEveryCheckPoint;    // Check point files will be written every 'nEveryCheckPoint' time steps
  //! Number of blade for every turbine
  std::vector<int> numBlades;           // Number of blades
  //! Desired number of actuator points on each blade for every turbine
  std::vector<int> numForcePtsBlade;
  //! Desired number of actuator points on the tower for every turbine
  std::vector<int> numForcePtsTwr;
  //! Number of Aerodyn points on each blade for every turbine
  std::vector<int> numVelPtsBlade;
  //! Number of Aerodyn points on the tower for every turbine
  std::vector<int> numVelPtsTwr;
  //! Number of inputs and output to the supercontroller from/to each turbine    
  int numScOutputs;  // # outputs from the supercontroller == # inputs to the controller == NumSC2Ctrl
  int numScInputs;   // # inputs to the supercontroller == # outputs from the controller == NumCtrl2SC
  //! Number of inputs and outputs to the supercontroller aggregated across all turbines and processors
  std::vector<double> scOutputsGlob;  // # outputs from the supercontroller for all turbines
  std::vector<double> scInputsGlob;   // # inputs to the supercontroller for all turbines

  //! Velocity at force nodes - Store temporarily to interpolate to the velocity nodes
  std::vector<std::vector<std::vector<double> > > forceNodeVel;
  //! Position and velocity data at the velocity (aerodyn) nodes - (nTurbines, nTimesteps * nPoints * 6)
  std::vector<std::vector<double> > velNodeData;
  //! HDF-5 tag of file containing velocity (aerodyn) node data file
  hid_t velNodeDataFile; 

  //! Data structure to get forces and deflections from OpenFAST
  std::vector<OpFM_InputType_t> cDriver_Input_from_FAST;
  //! Data structure to send velocity information to OpenFAST    
  std::vector<OpFM_OutputType_t> cDriver_Output_to_FAST;

  //! Data structure to get inputs to supercontroller from turbines in OpenFAST
  std::vector<SC_InputType_t> cDriverSC_Input_from_FAST;
  //! Data structure to send outputs from supercontroller to turbines in OpenFAST    
  std::vector<SC_OutputType_t> cDriverSC_Output_to_FAST;

  // Turbine Number is DIFFERENT from TurbID. Turbine Number simply runs from 0:n-1 locally and globally.
  //! Mapping global turbine number to processor number
  std::map<int, int> turbineMapGlobToProc;
  //! Mapping local to global turbine number
  std::map<int, int> turbineMapProcToGlob;
  //! Reverse Mapping global turbine number to local turbine number
  std::map<int, int> reverseTurbineMapProcToGlob;
  //! Set of processors containing atleast one turbine 
  std::set<int> turbineSetProcs;
  //! Same as the turbineSetProcs, but as an integer array
  std::vector<int> turbineProcs; 

  // Supercontroller stuff
  bool scStatus;
  std::string scLibFile;
  // Dynamic load stuff copied from 'C++ dlopen mini HOWTO' on tldp.org
  void *scLibHandle ; 
  typedef SuperController* create_sc_t(); 
  create_sc_t * create_SuperController;
  typedef void destroy_sc_t(SuperController *); 
  destroy_sc_t * destroy_SuperController;
  SuperController * sc;

  //! Number of processors in a fastMPIGroup
  int fastMPIGroupSize;
  //! An MPI group created among all processors that simulate atleast one turbine
  MPI_Group fastMPIGroup;
  //! An MPI communicator for the MPI group created among all processors that simulate atleast one turbine
  MPI_Comm  fastMPIComm;
  //! MPI rank of processor on the fastMPIComm
  int fastMPIRank;

  //! Global MPI group
  MPI_Group worldMPIGroup;
  //! MPI rank of processor on global MPI Comm
  int worldMPIRank;

  //! Error status and Error message to communicate with OpenFAST
  int ErrStat;
  char ErrMsg[INTERFACE_STRING_LENGTH];  // make sure this is the same size as IntfStrLen in FAST_Library.f90

 public: 

  //! Constructor 
  OpenFAST() ;
  
  //! Destructor
  ~OpenFAST() {} ;

  //! Set inputs to OpenFAST through an object of the class fastInputs
  void setInputs(const fastInputs &);  

  //! Initialize the simulation - allocate memory for all data structures and initialize all turbines
  void init();
  //! Call FAST->solution0 for all turbines
  void solution0();
  //! Advance all turbines by 1 OpenFAST timestep. Also write velocity data at the Aerodyn nodes for all turbines to 'velData.h5'
  void step();
  //! Advance all turbines by 1 OpenFAST timestep. Will not write any velocity data to file
  void stepNoWrite();
  //! Call FAST->end for all turbines
  void end();
  //! Open 'velData.h5' file to in preparation to write velocity data for all turbines
  hid_t openVelocityDataFile(bool createFile);
  //! Read velocity data at the Aerodyn nodes for all turbines from 'velData.h5'
  void readVelocityData(int nTimesteps);
  //! Write velocity data at the Aerodyn nodes for all turbines to 'velData.h5'  
  void writeVelocityData(hid_t h5file, int iTurb, int iTimestep, OpFM_InputType_t iData, OpFM_OutputType_t oData);
  //! Close 'velData.h5' file 
  herr_t closeVelocityDataFile(int nt_global, hid_t velDataFile);

  //! Allocate turbine number 'iTurbGlob' to the processor with global MPI rank 'procNo'
  void setTurbineProcNo(int iTurbGlob, int procNo) { turbineMapGlobToProc[iTurbGlob] = procNo; }
  //! Allocate all turbines to processors in a round-robin fashion
  void allocateTurbinesToProcsSimple();
  //! Get the approximate hub position for turbine number 'iTurbGlob'. This is the value specified in the input to OpenFAST. 
  void getApproxHubPos(std::vector<double> & currentCoords, int iTurbGlob);
  //! Get the exact hub position for turbine number 'iTurbGlob'. This is avaiable only after OpenFAST has been initialized for a given turbine.
  void getHubPos(std::vector<double> & currentCoords, int iTurbGlob);
  //! Get a vector pointing downstream along the hub for turbine number 'iTurbGlob'
  void getHubShftDir(std::vector<double> & hubShftVec, int iTurbGlob);

  //! Get the node type (HUB, BLADE, TOWER) of velocity node number 'iNode' for turbine number 'iTurbGlob'
  ActuatorNodeType getVelNodeType(int iTurbGlob, int iNode);
  //! Get the coordinates of velocity node number 'iNode' for turbine number 'iTurbGlob'
  void getVelNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob);
  //! Set the velocity at velocity node 'iNode' for turbine number 'iTurbGlob'
  void setVelocity(std::vector<double> & velocity, int iNode, int iTurbGlob);
  //! Set the velocity at force node 'iNode' for turbine number 'iTurbGlob'
  void setVelocityForceNode(std::vector<double> & velocity, int iNode, int iTurbGlob);
  //! Map the velocity from the force nodes to the velocity nodes using linear interpolation along each blade and the tower.
  void interpolateVel_ForceToVelNodes();
  //! Get the node type (HUB, BLADE, TOWER) of force node number 'iNode' for turbine number 'iTurbGlob'
  ActuatorNodeType getForceNodeType(int iTurbGlob, int iNode);
  //! Get the coordinates of force node number 'iNode' for turbine number 'iTurbGlob'
  void getForceNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob);
  //! Get the tensor orientation of force node number 'iNode' for turbine number 'iTurbGlob'
  void getForceNodeOrientation(std::vector<double> & currentOrientation, int iNode, int iTurbGlob);
  //! Get the actuator force at force node 'iNode' for turbine number 'iTurbGlob'
  void getForce(std::vector<double> & force, int iNode, int iTurbGlob);
  //! Get the chord at force node 'iNode' for turbine number 'iTurbGlob'
  double getChord(int iNode, int iTurbGlob);

  //! Get the starting time step of the simulation
  int get_ntStart() { return ntStart; }
  //! Return a boolean flag whether the simulation is dryRun
  bool isDryRun() { return dryRun; }
  //! Return a boolean flag whether the simulation is debug
  bool isDebug() { return debug; }
  //! Get an enum of type 'simStartType' indicating the start type of the simulation
  simStartType get_simStartType() { return simStart; }
  //! Is the simulation time zero right now?
  bool isTimeZero() { return timeZero; }
  //! Get the global MPI rank of the processor containing turbine number 'iTurbGlob'
  int get_procNo(int iTurbGlob) { return turbineMapGlobToProc[iTurbGlob] ; }
  //! Get the local turbine number of the turbine number 'iTurbGlob'
  int get_localTurbNo(int iTurbGlob) { return reverseTurbineMapProcToGlob[iTurbGlob]; }
  //! Get the total number of turbines in the simulation
  int get_nTurbinesGlob() { return nTurbinesGlob; } 

  //! Get the number of blades in turbine number 'iTurbGlob'
  int get_numBlades(int iTurbGlob) { return get_numBladesLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Aerodyn/velocity nodes on each blade in turbine number 'iTurbGlob'
  int get_numVelPtsBlade(int iTurbGlob) { return get_numVelPtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Aerodyn/velocity nodes on the tower in turbine number 'iTurbGlob'
  int get_numVelPtsTwr(int iTurbGlob) { return get_numVelPtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the total number of Aerodyn/velocity nodes in turbine number 'iTurbGlob'
  int get_numVelPts(int iTurbGlob) { return get_numVelPtsLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Actuator/force nodes on each blade in turbine number 'iTurbGlob'
  int get_numForcePtsBlade(int iTurbGlob) { return get_numForcePtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Actuator/force nodes on the tower in turbine number 'iTurbGlob'
  int get_numForcePtsTwr(int iTurbGlob) { return get_numForcePtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the total number of Actuator/force nodes in turbine number 'iTurbGlob'
  int get_numForcePts(int iTurbGlob) { return get_numForcePtsLoc(get_localTurbNo(iTurbGlob)); }

  //! Compute the torque and thrust for turbine number 'iTurbGlob'
  void computeTorqueThrust(int iTurGlob, std::vector<double> &  torque, std::vector<double> &  thrust);

 private:

  //! Check whether the error status is ok. If not quit gracefully by printing the error message
  void checkError(const int ErrStat, const char * ErrMsg);
  //! Check whether a file with name "name" exists
  inline bool checkFileExists(const std::string& name);

  //! Allocate memory for data structures for all turbines on this processor
  void allocateMemory();

  //! Get the number of blades in local turbine number 'iTurbLoc'
  int get_numBladesLoc(int iTurbLoc) { return numBlades[iTurbLoc]; }
  //! Get the number of Aerodyn/velocity nodes on each blade in local turbine number 'iTurbLoc'
  int get_numVelPtsBladeLoc(int iTurbLoc) { return numVelPtsBlade[iTurbLoc]; }
  //! Get the number of Aerodyn/velocity nodes on the tower in local turbine number 'iTurbLoc'
  int get_numVelPtsTwrLoc(int iTurbLoc) { return numVelPtsTwr[iTurbLoc]; }
  //! Get the total number of Aerodyn/velocity nodes in local turbine number 'iTurbLoc'
  int get_numVelPtsLoc(int iTurbLoc) { return 1 + numBlades[iTurbLoc]*numVelPtsBlade[iTurbLoc] + numVelPtsTwr[iTurbLoc]; }
  //! Get the number of Actuator/force nodes on each blade in local turbine number 'iTurbLoc'
  int get_numForcePtsBladeLoc(int iTurbLoc) { return numForcePtsBlade[iTurbLoc]; }
  //! Get the number of Actuator/force nodes on the tower in local turbine number 'iTurbLoc'
  int get_numForcePtsTwrLoc(int iTurbLoc) { return numForcePtsTwr[iTurbLoc]; }
  //! Get the total number of Actuator/force nodes in local turbine number 'iTurbLoc'
  int get_numForcePtsLoc(int iTurbLoc) { return 1 + numBlades[iTurbLoc]*numForcePtsBlade[iTurbLoc] + numForcePtsTwr[iTurbLoc]; }

  void loadSuperController(const fastInputs & fi);
  void fillScInputsGlob() ;
  void fillScOutputsLoc() ;

  //! An example function to set velocities at the Aerodyn nodes    
  void setOutputsToFAST(OpFM_InputType_t cDriver_Input_from_FAST, OpFM_OutputType_t cDriver_Output_to_FAST) ;
  //! Apply the velocity data at the Aerodyn nodes in 'velData' to turbine number 'iTurb' at time step 'iPrestart' through the data structure 'cDriver_Output_to_FAST'
  void applyVelocityData(int iPrestart, int iTurb, OpFM_OutputType_t cDriver_Output_to_FAST, std::vector<double> & velData) ;

};

}

#endif
