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
#include "yaml-cpp/yaml.h"
#include "dlfcn.h"
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include "mpi.h"
#include "SC.h"


namespace fast {

struct globTurbineDataType {
  int TurbID;
  std::string FASTInputFileName;
  std::string FASTRestartFileName;
  std::vector<double> TurbineBasePos;
  std::vector<double> TurbineHubPos;
  int numForcePtsBlade;
  int numForcePtsTwr;
};

enum ActuatorNodeType {
  HUB = 0,
  BLADE = 1,
  TOWER = 2,
  ActuatorNodeType_END
};

class OpenFAST {

 private:

  bool dryRun;        // If this is true, class will simply go through allocation and deallocation of turbine data
  bool debug;   // Write out extra information if this flags is turned on
  bool timeZero; // Flag to determine whether the Solution0 function has to be called or not.
  globTurbineDataType *  globTurbineData;
  int nTurbinesProc;
  int nTurbinesGlob;
  bool restart;
  double dtFAST;
  double tMax;
  float ** TurbineBasePos;
  int * TurbID;
  char ** FASTInputFileName;
  char ** CheckpointFileRoot;
  double tStart, tEnd;
  int nt_global;           
  int ntStart, ntEnd;      // The time step to start and end the FAST simulation
  int nEveryCheckPoint;    // Check point files will be written every 'nEveryCheckPoint' time steps
  int * numBlades;           // Number of blades
  int * numForcePtsBlade;
  int * numForcePtsTwr;
  int * numVelPtsBlade;
  int * numVelPtsTwr;
  int numScOutputs;  // # outputs from the supercontroller == # inputs to the controller == NumSC2Ctrl
  int numScInputs;   // # inputs to the supercontroller == # outputs from the controller == NumCtrl2SC
  double ** scOutputsGlob;  // # outputs from the supercontroller for all turbines
  double ** scInputsGlob;   // # inputs to the supercontroller for all turbines
  
  OpFM_InputType_t * cDriver_Input_from_FAST;
  OpFM_OutputType_t * cDriver_Output_to_FAST;

  SC_InputType_t * cDriverSC_Input_from_FAST;
  SC_OutputType_t * cDriverSC_Output_to_FAST;

  // Turbine Number is DIFFERENT from TurbID. Turbine Number simply runs from 0:n-1 locally and globally.
  std::map<int, int> turbineMapGlobToProc; // Mapping global turbine number to processor number
  std::map<int, int> turbineMapProcToGlob; // Mapping local to global turbine number
  std::map<int, int> reverseTurbineMapProcToGlob; // Reverse Mapping global turbine number to local turbine number
  std::set<int> turbineSetProcs; // Set of processors containing atleast one turbine 
  int * turbineProcs; // Same as the turbineSetProcs, but as an integer array

  //Supercontroller stuff
  bool scStatus;
  std::string scLibFile;
  // Dynamic load stuff copied from 'C++ dlopen mini HOWTO' on tldp.org
  void *scLibHandle ; 
  typedef SuperController* create_sc_t(); 
  create_sc_t * create_SuperController;
  typedef void destroy_sc_t(SuperController *); 
  destroy_sc_t * destroy_SuperController;
  SuperController * sc;

  int fastMPIGroupSize;
  MPI_Group fastMPIGroup;
  MPI_Comm  fastMPIComm;
  int fastMPIRank;

  MPI_Group worldMPIGroup;
  int worldMPIRank;

  int ErrStat;
  char ErrMsg[INTERFACE_STRING_LENGTH];  // make sure this is the same size as IntfStrLen in FAST_Library.f90

 public: 

  // Constructor 
  OpenFAST() ;
  
  // Destructor
  ~OpenFAST() {} ;
  
  int readInputFile(std::string cInterfaceInputFile);  
  int readInputFile(const YAML::Node &);  
  void setRestart(const bool & isRestart);
  void setTstart(const double & cfdTstart);
  void setDt(const double & cfdDt);
  void setTend(const double & cfdTend);

  void allocateInputData();
  int init();
  int solution0();
  int step();
  void end();

  int setTurbineProcNo(int iTurbGlob, int procNo) { turbineMapGlobToProc[iTurbGlob] = procNo; }
  void allocateTurbinesToProcsSimple();
  void getHubPos(double *currentCoords, int iTurbGlob);

  ActuatorNodeType getVelNodeType(int iTurbGlob, int iNode);
  void getVelNodeCoordinates(double *currentCoords, int iNode, int iTurbGlob);
  void setVelocity(std::vector<double> & velocity, int iNode, int iTurbGlob);

  ActuatorNodeType getForceNodeType(int iTurbGlob, int iNode);
  void getForceNodeCoordinates(double *currentCoords, int iNode, int iTurbGlob);
  void getForceNodeOrientation(double *currentOrientation, int iNode, int iTurbGlob);
  void getForce(std::vector<double> & force, int iNode, int iTurbGlob);
  double getChord(int iNode, int iTurbGlob);

  int get_ntStart() { return ntStart; }
  int get_ntEnd() { return ntEnd; }
  bool isDryRun() { return dryRun; }
  bool isDebug() { return debug; }
  bool isTimeZero() { return timeZero; }
  int get_procNo(int iTurbGlob) { return turbineMapGlobToProc[iTurbGlob] ; } // Get processor number of a turbine with global id 'iTurbGlob'
  int get_localTurbNo(int iTurbGlob) { return reverseTurbineMapProcToGlob[iTurbGlob]; }
  int get_nTurbinesGlob() { return nTurbinesGlob; } 

  int get_numBlades(int iTurbGlob) { return get_numBladesLoc(get_localTurbNo(iTurbGlob)); }
  int get_numVelPtsBlade(int iTurbGlob) { return get_numVelPtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  int get_numVelPtsTwr(int iTurbGlob) { return get_numVelPtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  int get_numVelPts(int iTurbGlob) { return get_numVelPtsLoc(get_localTurbNo(iTurbGlob)); }
  int get_numForcePtsBlade(int iTurbGlob) { return get_numForcePtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  int get_numForcePtsTwr(int iTurbGlob) { return get_numForcePtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  int get_numForcePts(int iTurbGlob) { return get_numForcePtsLoc(get_localTurbNo(iTurbGlob)); }

  void computeTorqueThrust(int iTurbGlob, double * torque, double * thrust);

 private:

  void checkError(const int ErrStat, const char * ErrMsg);
  inline bool checkFileExists(const std::string& name);

  void readTurbineData(int iTurb, YAML::Node turbNode);

  int get_numBladesLoc(int iTurbLoc) { return numBlades[iTurbLoc]; }
  int get_numVelPtsBladeLoc(int iTurbLoc) { return numVelPtsBlade[iTurbLoc]; }
  int get_numVelPtsTwrLoc(int iTurbLoc) { return numVelPtsTwr[iTurbLoc]; }
  int get_numVelPtsLoc(int iTurbLoc) { return 1 + numBlades[iTurbLoc]*numVelPtsBlade[iTurbLoc] + numVelPtsTwr[iTurbLoc]; }
  int get_numForcePtsBladeLoc(int iTurbLoc) { return numForcePtsBlade[iTurbLoc]; }
  int get_numForcePtsTwrLoc(int iTurbLoc) { return numForcePtsTwr[iTurbLoc]; }
  int get_numForcePtsLoc(int iTurbLoc) { return 1 + numBlades[iTurbLoc]*numForcePtsBlade[iTurbLoc] + numForcePtsTwr[iTurbLoc]; }

  void loadSuperController(YAML::Node c);
  void fillScInputsGlob() ;
  void fillScOutputsLoc() ;

  void setOutputsToFAST(OpFM_InputType_t* cDriver_Input_from_FAST, OpFM_OutputType_t* cDriver_Output_to_FAST) ; // An example to set velocities at the Aerodyn nodes

};

}

#endif
