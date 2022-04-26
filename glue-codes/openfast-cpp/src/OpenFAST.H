#ifndef OpenFAST_h
#define OpenFAST_h
#include "FAST_Library.h"
#include "sys/stat.h"
#include <string>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>
#include "dlfcn.h"
//TODO: The skip MPICXX is put in place primarily to get around errors in OpenFOAM. This will cause problems if the driver program uses C++ API for MPI.
#ifndef OMPI_SKIP_MPICXX
 #define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
 #define MPICH_SKIP_MPICXX
#endif
#include "mpi.h"
#include "SC.h"


namespace fast {

struct globTurbineDataType {
  int TurbID;
  std::string FASTInputFileName;
  std::string FASTRestartFileName;
  std::vector<double> TurbineBasePos;
  std::vector<double> TurbineHubPos;
  std::string forcePtsBladeDistributionType;
  int numForcePtsBlade;
  int numForcePtsTwr;
  float nacelle_cd{0.0};
  float nacelle_area{0.0};
  float air_density{0.0};
};

enum ActuatorNodeType {
  HUB = 0,
  BLADE = 1,
  TOWER = 2,
  ActuatorNodeType_END
};

enum simStartType {
  init = 0,
  trueRestart = 1,
  restartDriverInitFAST = 2,
  simStartType_END
};


class fastInputs {

 public:

  MPI_Comm comm;
  int nTurbinesGlob;
  bool dryRun;
  bool debug;
  double tStart;
  simStartType simStart;
  int nEveryCheckPoint;
  double tMax;
  double dtFAST;

  bool scStatus;
  std::string scLibFile;
  std::vector<globTurbineDataType>  globTurbineData;

  // Constructor
  fastInputs() ;

  // Destructor
  ~fastInputs() {} ;

};


class OpenFAST {

 private:

  MPI_Comm mpiComm;
  bool dryRun;        // If this is true, class will simply go through allocation and deallocation of turbine data
  bool debug;   // Write out extra information if this flags is turned on
  std::vector<globTurbineDataType> globTurbineData;
  int nTurbinesProc;
  int nTurbinesGlob;
  simStartType simStart;
  bool timeZero;
  double dtFAST;
  double tMax;
  std::vector<std::vector<float> > TurbineBasePos;
  std::vector<std::vector<float> > TurbineHubPos;
  std::vector<int> TurbID;
  std::vector<std::string > FASTInputFileName;
  std::vector<std::string > CheckpointFileRoot;
  std::vector<float> nacelle_cd;
  std::vector<float> nacelle_area;
  std::vector<float> air_density;
  double tStart;
  int nt_global;
  int ntStart;      // The time step to start the FAST simulation
  int nEveryCheckPoint;    // Check point files will be written every 'nEveryCheckPoint' time steps
  std::vector<int> numBlades;           // Number of blades
  std::vector<std::string> forcePtsBladeDistributionType;
  std::vector<int> numForcePtsBlade;
  std::vector<int> numForcePtsTwr;
  std::vector<int> numVelPtsBlade;
  std::vector<int> numVelPtsTwr;

  std::vector<std::vector<std::vector<double> > > forceNodeVel; // Velocity at force nodes - Store temporarily to interpolate to the velocity nodes
  std::vector<std::vector<double> > velNodeData; // Position and velocity data at the velocity (aerodyn) nodes - (nTurbines, nTimesteps * nPoints * 6)
  hid_t velNodeDataFile; // HDF-5 tag of file containing velocity (aerodyn) node data file

  std::vector<OpFM_InputType_t> cDriver_Input_from_FAST;
  std::vector<OpFM_OutputType_t> cDriver_Output_to_FAST;

  // Turbine Number is DIFFERENT from TurbID. Turbine Number simply runs from 0:n-1 locally and globally.
  std::map<int, int> turbineMapGlobToProc; // Mapping global turbine number to processor number
  std::map<int, int> turbineMapProcToGlob; // Mapping local to global turbine number
  std::map<int, int> reverseTurbineMapProcToGlob; // Reverse Mapping global turbine number to local turbine number
  std::set<int> turbineSetProcs; // Set of processors containing at least one turbine
  std::vector<int> turbineProcs; // Same as the turbineSetProcs, but as an integer array

  //Supercontroller stuff
  bool scStatus;
  SuperController sc;
  scInitOutData scio;

  int fastMPIGroupSize;
  MPI_Group fastMPIGroup;
  MPI_Comm  fastMPIComm;
  int fastMPIRank;

  MPI_Group worldMPIGroup;
  int worldMPIRank;

  static int AbortErrLev;
  int ErrStat;
  char ErrMsg[INTERFACE_STRING_LENGTH];  // make sure this is the same size as IntfStrLen in FAST_Library.f90

 public:

  // Constructor
  OpenFAST() ;

  // Destructor
  ~OpenFAST() ;

  void setInputs(const fastInputs &);

  void init();
  void solution0();
  void step();
  void stepNoWrite();
  void end();

  // Compute the nacelle force
  void calc_nacelle_force(const float & u,
                          const float & v,
                          const float & w,
                          const float & cd,
                          const float & area,
                          const float & rho,
                          float & fx,
                          float & fy,
                          float & fz);

  hid_t openVelocityDataFile(bool createFile);
  void readVelocityData(int nTimesteps);
  void writeVelocityData(hid_t h5file, int iTurb, int iTimestep, OpFM_InputType_t iData, OpFM_OutputType_t oData);
  herr_t closeVelocityDataFile(int nt_global, hid_t velDataFile);
  void backupVelocityDataFile(int curTimeStep, hid_t & velDataFile);

  void setTurbineProcNo(int iTurbGlob, int procNo) { turbineMapGlobToProc[iTurbGlob] = procNo; }
  void allocateTurbinesToProcsSimple();
  void getApproxHubPos(double* currentCoords, int iTurbGlob, int nSize=3);
  void getHubPos(double* currentCoords, int iTurbGlob, int nSize=3);
  void getHubShftDir(double* hubShftVec, int iTurbGlob, int nSize=3);

  ActuatorNodeType getVelNodeType(int iTurbGlob, int iNode);
  void getVelNodeCoordinates(double* currentCoords, int iNode, int iTurbGlob, int nSize=3);
  void setVelocity(double* velocity, int iNode, int iTurbGlob, int nSize=3);
  void setVelocityForceNode(double* velocity, int iNode, int iTurbGlob, int nSize=3);
  void interpolateVel_ForceToVelNodes();
  ActuatorNodeType getForceNodeType(int iTurbGlob, int iNode);
  void getForceNodeCoordinates(double* currentCoords, int iNode, int iTurbGlob, int nSize=3);
  void getForceNodeOrientation(double* currentOrientation, int iNode, int iTurbGlob, int nSize=9);
  void getForce(double* force, int iNode, int iTurbGlob, int nSize=3);
  void getRelativeVelForceNode(double* vel, int iNode, int iTurbGlob, int nSize=3);
  double getChord(int iNode, int iTurbGlob);

  int get_ntStart() { return ntStart; }
  bool isDryRun() { return dryRun; }
  bool isDebug() { return debug; }
  simStartType get_simStartType() { return simStart; }
  bool isTimeZero() { return timeZero; }
  int get_procNo(int iTurbGlob) { return turbineMapGlobToProc[iTurbGlob] ; } // Get processor number of a turbine with global id 'iTurbGlob'
  int get_localTurbNo(int iTurbGlob) { return reverseTurbineMapProcToGlob[iTurbGlob]; }
  int get_nTurbinesGlob() { return nTurbinesGlob; }

  float get_nacelleArea(int iTurbGlob) { return get_nacelleAreaLoc(get_localTurbNo(iTurbGlob)); }
  float get_nacelleCd(int iTurbGlob) { return get_nacelleCdLoc(get_localTurbNo(iTurbGlob)); }
  float get_airDensity(int iTurbGlob) { return get_airDensityLoc(get_localTurbNo(iTurbGlob)); }
  int get_numBlades(int iTurbGlob) { return get_numBladesLoc(get_localTurbNo(iTurbGlob)); }
  int get_numVelPtsBlade(int iTurbGlob) { return get_numVelPtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  int get_numVelPtsTwr(int iTurbGlob) { return get_numVelPtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  int get_numVelPts(int iTurbGlob) { return get_numVelPtsLoc(get_localTurbNo(iTurbGlob)); }
  int get_numForcePtsBlade(int iTurbGlob) { return get_numForcePtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  int get_numForcePtsTwr(int iTurbGlob) { return get_numForcePtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  int get_numForcePts(int iTurbGlob) { return get_numForcePtsLoc(get_localTurbNo(iTurbGlob)); }

  void computeTorqueThrust(int iTurGlob, std::vector<double> &  torque, std::vector<double> &  thrust);


  inline
  void getHubPos(std::vector<double> & currentCoords, int iTurbGlob) {
    getHubPos(currentCoords.data(), iTurbGlob, currentCoords.size());
  }

  inline
  void getApproxHubPos(std::vector<double>& currentCoords, int iTurbGlob) {
    getApproxHubPos(currentCoords.data(), iTurbGlob, currentCoords.size());
  }

  inline
  void getHubShftDir(std::vector<double> & hubShftVec, int iTurbGlob) {
    getHubShftDir(hubShftVec.data(), iTurbGlob, hubShftVec.size());
  }

  inline
  void getVelNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob) {
    getVelNodeCoordinates(currentCoords.data(), iNode, iTurbGlob, currentCoords.size());
  }

  inline
  void setVelocity(std::vector<double> & currentVelocity, int iNode, int iTurbGlob) {
    setVelocity(currentVelocity.data(), iNode, iTurbGlob, currentVelocity.size());
  }

  inline
  void setVelocityForceNode(std::vector<double> & currentVelocity, int iNode, int iTurbGlob) {
    setVelocityForceNode(currentVelocity.data(), iNode, iTurbGlob, currentVelocity.size());
  }

  inline
  void getForceNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob) {
    getForceNodeCoordinates(currentCoords.data(), iNode, iTurbGlob, currentCoords.size());
  }

  inline
  void getForceNodeOrientation(std::vector<double> & currentOrientation, int iNode, int iTurbGlob) {
    getForceNodeOrientation(currentOrientation.data(), iNode, iTurbGlob, currentOrientation.size());
  }

  inline
  void getForce(std::vector<double> & currentForce, int iNode, int iTurbGlob) {
    getForce(currentForce.data(), iNode, iTurbGlob, currentForce.size());
  }

  inline
  void getRelativeVelForceNode(std::vector<double> & currentVelocity, int iNode, int iTurbGlob) {
    getRelativeVelForceNode(currentVelocity.data(), iNode, iTurbGlob, currentVelocity.size());
  }

 private:

  void checkError(const int ErrStat, const char * ErrMsg);
  inline bool checkFileExists(const std::string& name);

  void allocateMemory();

  float get_nacelleCdLoc(int iTurbLoc) { return nacelle_cd[iTurbLoc]; }
  float get_nacelleAreaLoc(int iTurbLoc) { return nacelle_area[iTurbLoc]; }
  float get_airDensityLoc(int iTurbLoc) { return air_density[iTurbLoc]; }
  int get_numBladesLoc(int iTurbLoc) { return numBlades[iTurbLoc]; }
  int get_numVelPtsBladeLoc(int iTurbLoc) { return numVelPtsBlade[iTurbLoc]; }
  int get_numVelPtsTwrLoc(int iTurbLoc) { return numVelPtsTwr[iTurbLoc]; }
  int get_numVelPtsLoc(int iTurbLoc) { return 1 + numBlades[iTurbLoc]*numVelPtsBlade[iTurbLoc] + numVelPtsTwr[iTurbLoc]; }
  int get_numForcePtsBladeLoc(int iTurbLoc) { return numForcePtsBlade[iTurbLoc]; }
  int get_numForcePtsTwrLoc(int iTurbLoc) { return numForcePtsTwr[iTurbLoc]; }
  int get_numForcePtsLoc(int iTurbLoc) { return 1 + numBlades[iTurbLoc]*numForcePtsBlade[iTurbLoc] + numForcePtsTwr[iTurbLoc]; }

  void loadSuperController(const fastInputs & fi);

  void setOutputsToFAST(OpFM_InputType_t cDriver_Input_from_FAST, OpFM_OutputType_t cDriver_Output_to_FAST) ; // An example to set velocities at the Aerodyn nodes
  void applyVelocityData(int iPrestart, int iTurb, OpFM_OutputType_t cDriver_Output_to_FAST, std::vector<double> & velData) ;

};

}

#endif
