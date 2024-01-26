#ifndef OpenFAST_h
#define OpenFAST_h
#include "FAST_Library.h"
#include "sys/stat.h"
#include <string>
#include <cstring>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <set>
#include <map>
#include <memory>
#include "netcdf.h"
#include "dlfcn.h"
//TODO: The skip MPICXX is put in place primarily to get around errors in ExternalInflow. This will cause problems if the driver program uses C++ API for MPI.
#ifndef OMPI_SKIP_MPICXX
 #define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
 #define MPICH_SKIP_MPICXX
#endif
#include "mpi.h"
#include "SC.h"


namespace fast {

//! An id to indicate the type of simulation for each turbine - Simple/Actuator with optional externally specified inflow or Blade-Resolved with externally specified loads
enum simType {
    EXTINFLOW = 0,
    EXTLOADS  = 1,
    simType_END
};

//! A data structure to hold all turbine related information
struct turbineDataType {
    //!Integer id for every turbine
    int TurbID;
    //! The FAST Input file name. Typically a .fst file.
    std::string FASTInputFileName;
    //! The restart/checkpoint file name.
    std::string FASTRestartFileName;
    //! Output file root
    std::string outFileRoot;
    //! The time step for OpenFAST for this turbine
    double dt;
    //! The position of the base of the turbine in global coordinates
    std::vector<float> TurbineBasePos;
    //! The approximate position of the hub of the turbine in global coordinates
    std::vector<double> TurbineHubPos;
    //! Simulation type
    simType sType;
    //! Number of blades
    int numBlades;
    //! Number of velocity nodes (AeroDyn) per blade
    int numVelPtsBlade;
    //! Number of velocity nodes (AeroDyn) on the tower
    int numVelPtsTwr;
    //! Total number of velocity nodes (AeroDyn)
    int numVelPts;
    //! Desired number of actuator points on each blade
    int numForcePtsBlade;
    //! Desired number of actuator points on the tower
    int numForcePtsTwr;
    //! Total number of actuator points
    int numForcePts;
    //! Node clustering type
    int nodeClusterType;    
    //! Inflow Type - 1 (InflowWind) or 2 (Externally specified)
    int inflowType;
    //! Drag coefficient of nacelle
    float nacelle_cd;
    //! Frontal area of the nacelle
    float nacelle_area;
    //! Air density around this turbine
    float air_density;
    //! Number of nodes at which the forces and deflections are computed for blade-resolved FSI on each blade
    std::vector<int> nBRfsiPtsBlade;
    //! Total number of BR fsi points on all blades combined
    int nTotBRfsiPtsBlade;
    //! Number of nodes at which the forces and deflections are computed for blade-resolved FSI on the tower
    int nBRfsiPtsTwr;
    //! The mean azimuth at which the loads are blended between AeroDyn and CFD
    double azBlendMean;
    //! The delta azimuth over which the the loads are blended between AeroDyn and CFD
    double azBlendDelta;
    //! Mean velocity at reference height
    double velMean;
    //! Compass angle of wind direction (in degrees)
    double windDir;
    //! Reference height for velocity profile
    double zRef;
    //! Shear exponent of velocity profile
    double shearExp;
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

//! A data structure to hold all velocity and force node information
struct turbVelForceNodeDataType {
    //! Blade location at velocity nodes
    std::vector<double> x_vel;
    //! Blade velocity at velocity nodes
    std::vector<double> xdot_vel;
    //! Sampled velocity at velocity nodes
    std::vector<double> vel_vel;
    //! Reference location at force nodes
    std::vector<double> xref_force;
    //! Blade location at force nodes
    std::vector<double> x_force;
    //! Blade velocity at force nodes
    std::vector<double> xdot_force;
    //! Blade orientation at force nodes
    std::vector<double> orient_force;
    //! Sampled velocity at force nodes
    std::vector<double> vel_force;
    //! Actuator force at force nodes
    std::vector<double> force;
    double x_vel_resid;
    double xdot_vel_resid;
    double vel_vel_resid;
    double x_force_resid;
    double xdot_force_resid;
    double orient_force_resid;
    double vel_force_resid;
    double force_resid;
};

//! An enum to keep track of information stored at different time steps
enum timeStep {
    STATE_NM2 = 0,
    STATE_NM1 = 1,
    STATE_N = 2,
    STATE_NP1 = 3,
    timeStep_END
};

//! A data structure to hold all loads and deflections information for blade-resolved FSI simulations
struct turbBRfsiDataType {
    //! Tower reference position
    std::vector<double> twr_ref_pos;
    //! Tower deflections
    std::vector<double> twr_def;
    //! Tower velocity
    std::vector<double> twr_vel;
    //! Blade radial location
    std::vector<double> bld_rloc;
    //! Blade chord
    std::vector<double> bld_chord;
    //! Blade reference position
    std::vector<double> bld_ref_pos;
    //! Blade deflections
    std::vector<double> bld_def;
    //! Blade velocity
    std::vector<double> bld_vel;
    //! Hub reference position
    std::vector<double> hub_ref_pos;
    //! Hub deflections
    std::vector<double> hub_def;
    //! Hub velocity
    std::vector<double> hub_vel;
    //! Nacelle reference position
    std::vector<double> nac_ref_pos;
    //! Nacelle deflections
    std::vector<double> nac_def;
    //! Nacelle velocity
    std::vector<double> nac_vel;
    //! Blade root reference position
    std::vector<double> bld_root_ref_pos;
    //! Blade root deformation
    std::vector<double> bld_root_def;
    //! Blade pitch
    std::vector<double> bld_pitch;

    //! Tower loads
    std::vector<double> twr_ld;
    //! Blade loads
    std::vector<double> bld_ld;
    double twr_def_resid;
    double twr_vel_resid;
    double bld_def_resid;
    double bld_vel_resid;
    double twr_ld_resid;
    double bld_ld_resid;
};

/**
 * A class to hold all input data for a simulation run through a OpenFAST C++ glue code
 */
class fastInputs {

 public:

  //! MPI Communicator
  MPI_Comm comm;
  //! Total number of turbines in the simulation
  int nTurbinesGlob{0};
  //! The simulation will not run if dryRun is set to true. However, the simulation will read the input files, allocate turbines to processors and prepare to run the individual turbine instances. This flag is useful to test the setup of the simulation before running it.
  bool dryRun{false};
  //! Enable debug outputs if set to true
  bool debug{false};
  //! Start time of the simulation
  double tStart{-1.0};
  //! Start type of the simulation: 'INIT', 'TRUERESTART' or 'RESTARTDRIVERINITFAST'.
  simStartType simStart;
  //!Restart files will be written every so many time stneps
  int restartFreq{-1};
  //!Output files will be written every so many time stneps
  int outputFreq{100};
  //! Max time of the simulation
  double tMax{0.0};
  //! Time step for driver.
  double dtDriver{0.0};
  //! Time step for openfast.
  double dtFAST{0.0};
  //! Supercontroller status: True/False.
  bool scStatus{false};
  //! Name of the dynamic library containing the supercontroller implementation
  std::string scLibFile{""};
  //! Number of inputs and output to the supercontroller from/to each turbine
  int numScInputs{0};
  int numScOutputs{0};

  //! Vector of turbine specific input data
  std::vector<turbineDataType>  globTurbineData;

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
  bool dryRun{false};        // If this is true, class will simply go through allocation and deallocation of turbine data
  //! Enable debug outputs if set to true
  bool debug{false};   // Write out extra information if this flags is turned on
  //! Number of turbines on this MPI rank
  int nTurbinesProc{0};
  //! Total number of turbines in the simulation
  int nTurbinesGlob{0};
  //! Start type of the simulation: 'INIT', 'TRUERESTART' or 'RESTARTDRIVERINITFAST'.
  simStartType simStart{fast::init};
  //! Offset between driver and openfast simulation time  - t_driver - t_openfast
  double driverOpenfastOffset_{0.0};
  //! Is the time now zero: True/False
  bool timeZero{false};
  //! Time step for FAST. All turbines on a given processor should have the same time step.
  double dtFAST{-1.0};
  //! Time step for Driver.
  double dtDriver{-1.0};
  //! Number of OpenFAST time steps per unit time step of the driver program
  int nSubsteps_{-1};
  //! Is this the first pass through a time step
  bool firstPass_{true};
  //! Max time of the simulation
  double tMax{-1.0};
  //! Start time of the simulation
  double tStart{-1.0};

  //! The current time step number
  int nt_global{0};
  //! The current nonlinear iteration
  int nlinIter_{0};
  //! The starting time step number
  int ntStart{0};
  //! Restart files will be written every so many time steps
  int restartFreq_{-1};
  //! Output files will be written every so many time steps
  int outputFreq_{100};

  //! Map of `{variableName : netCDF_ID}` obtained from the NetCDF C interface
  std::vector<std::string> ncOutVarNames_;
  std::unordered_map<std::string, int> ncOutVarIDs_;

  //! Map of `{dimName : netCDF_ID}` obtained from the NetCDF C interface
  std::vector<std::string> ncOutDimNames_;
  std::unordered_map<std::string, int> ncOutDimIDs_;

  //! Map of `{variableName : netCDF_ID}` obtained from the NetCDF C interface
  std::vector<std::string> ncRstVarNames_;
  std::unordered_map<std::string, int> ncRstVarIDs_;

  //! Map of `{dimName : netCDF_ID}` obtained from the NetCDF C interface
  std::vector<std::string> ncRstDimNames_;
  std::unordered_map<std::string, int> ncRstDimIDs_;

  std::vector<turbineDataType> globTurbineData; //All turbines
  std::vector<turbineDataType> turbineData;   // Only for turbines on the proc

  //! Velocity at force nodes - Store temporarily to interpolate to the velocity nodes
  std::vector<std::vector<std::vector<double> > > forceNodeVel; // Velocity at force nodes - Store temporarily to interpolate to the velocity nodes
  //! Position and velocity data at the velocity (aerodyn) nodes - (nTurbines, nTimesteps * nPoints * 6)
  std::vector<std::vector<double> > velNodeData; // Position and velocity data at the velocity (aerodyn) nodes - (nTurbines, nTimesteps * nPoints * 6)
  //! Array containing data at the velocity and force nodes
  std::vector<std::vector<turbVelForceNodeDataType>> velForceNodeData;
  //! Array containing forces and deflections data for blade-resolved FSI simulations.
  std::vector<std::vector<turbBRfsiDataType>> brFSIData;

  //! Data structure to get forces and deflections from ExternalInflow module in OpenFAST
  std::vector<ExtInfw_InputType_t> extinfw_i_f_FAST; // Input from OpenFAST
  //! Data structure to send velocity information to ExternalInflow module in OpenFAST
  std::vector<ExtInfw_OutputType_t> extinfw_o_t_FAST; // Output to OpenFAST

  //! Data structure to get deflections from ExternalLoads module in OpenFAST
  std::vector<ExtLdDX_InputType_t> extld_i_f_FAST; // Input from OpenFAST
  //! Data structure to get deflections from ExternalLoads module in OpenFAST
  std::vector<ExtLdDX_ParameterType_t> extld_p_f_FAST; // Parameter from OpenFAST
  //! Data structure to send force information to ExternalLoads module in OpenFAST
  std::vector<ExtLdDX_OutputType_t> extld_o_t_FAST; // Output to OpenFAST

  scInitOutData scio;

  // Mapping of local turbine number to global turbine and processor number
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
    bool scStatus{false};
  std::string scLibFile;
  // Dynamic load stuff copied from 'C++ dlopen mini HOWTO' on tldp.org
  void *scLibHandle ;
  typedef SuperController* create_sc_t();
  create_sc_t * create_SuperController;
  typedef void destroy_sc_t(SuperController *);
  destroy_sc_t * destroy_SuperController;
  std::unique_ptr<SuperController> sc;

  // MPI related book keeping for all processors containing turbines
  //! Number of processors in a fastMPIGroup
  int fastMPIGroupSize{-1};
  //! An MPI group created among all processors that simulate atleast one turbine
  MPI_Group fastMPIGroup;
  //! An MPI communicator for the MPI group created among all processors that simulate atleast one turbine
  MPI_Comm  fastMPIComm;
  //! MPI rank of processor on the fastMPIComm
  int fastMPIRank{-1};

  //! Global MPI group
  MPI_Group worldMPIGroup;
  //! MPI rank of processor on global MPI Comm
  int worldMPIRank{-1};

  //! Error status and Error message to communicate with OpenFAST
  int ErrStat{0};
  char ErrMsg[INTERFACE_STRING_LENGTH];  // make sure this is the same size as IntfStrLen in FAST_Library.f90
  static int AbortErrLev;

 public:

  //! Constructor
  OpenFAST() ;

  //! Destructor
  ~OpenFAST() {} ;

  //! Set inputs to OpenFAST through an object of the class fastInputs. Should be called on all MPI ranks.
  void setInputs(const fastInputs &);

  //! Check and set the number of sub-timesteps
  int checkAndSetSubsteps();

  //! Set driver time step and check point
  void setDriverTimeStep(double dt_driver);
  void setDriverCheckpoint(int nt_checkpoint_driver);

  //! Initialize the simulation - allocate memory for all data structures and initialize all turbines. Safe to call on all MPI ranks.
  void init();
  //! Call FAST->solution0 for all turbines. Safe to call on all MPI ranks.
  void solution0(bool writeFiles=true);
  //! Initialize velocity and force node data. Safe to call on all MPI ranks.
  void init_velForceNodeData();
  //! Set up before every OpenFAST time step. Safe to call on all MPI ranks.
  void prework();
  //! Update states to next time step by calling FAST_AdvanceStates_T and CalcOutputs_And_SolveForInputs. Safe to call on all MPI ranks.
  void update_states_driver_time_step(bool writeFiles=true);
  //! Copy the final predicted states from step t_global_next to actual states for that step. Safe to call on all MPI ranks.
  void advance_to_next_driver_time_step(bool writeFiles=true);
  //! Set external inputs for OpenFAST modules by interpolating to substep. Safe to call on all MPI ranks.
  void send_data_to_openfast(double ss_time);
  //! Set external inputs for OpenFAST modules at time step 't'. Safe to call on all MPI ranks.
  void send_data_to_openfast(fast::timeStep t);
  //! Get ouput data from OpenFAST modules. Safe to call on all MPI ranks.
  void get_data_from_openfast(fast::timeStep t);
  //! Extrapolate velocity and force node data to time step 'n+1' using data at 'n', 'n-1' and 'n-2'. Safe to call on all MPI ranks.
  void predict_states();
  //! Advance all turbines by 1 OpenFAST timestep. Safe to call on all MPI ranks.
  void step(bool writeFiles=true);
  //! Step function to be used with sub-stepping fast between time steps of the driver program. Safe to call on all MPI ranks.
  void step(double ss_time);
  //! Call FAST->end for all turbines. Safe to call on all MPI ranks.
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

  //! Allocate turbine number 'iTurbGlob' to the processor with global MPI rank 'procNo'. MUST be called from every MPI rank.
  void setTurbineProcNo(int iTurbGlob, int procNo) { turbineMapGlobToProc[iTurbGlob] = procNo; }
  //! Allocate all turbines to processors in a round-robin fashion. MUST be called from every MPI rank.
  void allocateTurbinesToProcsSimple();

  //! Get fast time step on this processor
  double get_timestep() { return dtFAST; }

  //! Get the approximate hub position for turbine number 'iTurbGlob'. This is the value specified in the input to OpenFAST. Must be called only from the processor containing the turbine.
  void getApproxHubPos(double* currentCoords, int iTurbGlob, int nSize=3);
  //! Get the exact hub position for turbine number 'iTurbGlob'. This is avaiable only after OpenFAST has been initialized for a given turbine. Must be called only from the processor containing the turbine.

  void getHubPos(double* currentCoords, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=3);
  //! Get a vector pointing downstream along the hub for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  void getHubShftDir(double* hubShftVec, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=3);

  //! Get the node type (HUB, BLADE, TOWER) of velocity node number 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  ActuatorNodeType getVelNodeType(int iTurbGlob, int iNode);
  //! Get the coordinates of velocity node number 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  void getVelNodeCoordinates(double* currentCoords, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=3);
  //! Set the velocity at velocity node 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  void setVelocity(double* velocity, int iNode, int iTurbGlob, int nSize=3);
  //! Set the velocity at force node 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  void setVelocityForceNode(double* velocity, int iNode, int iTurbGlob, int nSize=3);
  //! Map the velocity from the force nodes to the velocity nodes using linear interpolation along each blade and the tower. Safe to call from every MPI rank.
  void interpolateVel_ForceToVelNodes();
  //! Get the node type (HUB, BLADE, TOWER) of force node number 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  ActuatorNodeType getForceNodeType(int iTurbGlob, int iNode);
  //! Get the coordinates of force node number 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  void getForceNodeCoordinates(double* currentCoords, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=3);
  //! Get the tensor orientation of force node number 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  void getForceNodeOrientation(double* currentOrientation, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=3);
  //! Get the actuator force at force node 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  void getForce(double* force, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=3);

  void getRelativeVelForceNode(double* vel, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=3);

  //! Get the chord at force node 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  double getChord(int iNode, int iTurbGlob);
  //! Get the radial location/height along blade/tower at force node 'iNode' for turbine number 'iTurbGlob'. Must be called only from the processor containing the turbine.
  double getRHloc(int iNode, int iTurbGlob);

  //! Get processor containing turbine 'iTurbGlob'
  int getProc(int iTurbGlob) {return turbineMapGlobToProc[iTurbGlob];}

  //! Get the blade chord array 'bldRloc' of turbine number 'iTurbGlob'
  void getBladeChord(double * bldChord, int iTurbGlob);
  //! Get the blade node radial locations array 'bldRloc' of turbine number 'iTurbGlob'
  void getBladeRloc(double * bldRloc, int iTurbGlob);
  //! Get the blade reference positions array 'bldRefPos' of turbine number 'iTurbGlob'
  void getBladeRefPositions(double* bldRefPos, int iTurbGlob, int nSize=6);
  //! Get the blade root reference positions array 'bldRootRefPos' of turbine number 'iTurbGlob'
  void getBladeRootRefPositions(double* bldRootRefPos, int iTurbGlob, int nSize=6);
  //! Get the blade deflections array 'bldDefl' of turbine number 'iTurbGlob' at time step 't'
  void getBladeDisplacements(double* bldDefl, double* bldVel, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=6);
  //! Get the blade root deflections array 'bldRootDefl' of turbine number 'iTurbGlob' at time step 't'
  void getBladeRootDisplacements(double* bldRootDefl, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=6);
  //! Get the blade pitch 'bldPitch' of turbine number 'iTurbGlob'
  void getBladePitch(double* bldPitch, int iTurbGlob, int nSize=3);
  //! Get the tower reference positions array 'twrRefPos' of turbine number 'iTurbGlob'
  void getTowerRefPositions(double* twrRefPos, int iTurbGlob, int nSize=6);
  //! Get the tower deflections array 'twrDefl' of turbine number 'iTurbGlob' at time step 't'
  void getTowerDisplacements(double* twrDefl, double* twrVel, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=6);
  //! Get the hub reference position array 'hubRefPos' of turbine number 'iTurbGlob'
  void getHubRefPosition(double* hubRefPos, int iTurbGlob, int nSize=6);
  //! Get the hub deflections array 'hubDefl' of turbine number 'iTurbGlob' at time step 't'
  void getHubDisplacement(double* hubDefl, double* hubVel, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=6);
  //! Get the nacelle reference position array 'nacRefPos' of turbine number 'iTurbGlob'
  void getNacelleRefPosition(double* nacRefPos, int iTurbGlob, int nSize=6);
  //! Get the nacelle deflections array 'nacDefl' of turbine number 'iTurbGlob' at time step 't'
  void getNacelleDisplacement(double* nacDefl, double* nacVel, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=6);

  //! Set the blade forces array 'bldForce' for blade 'iBlade' of turbine number 'iTurbGlob' at time step 't'
  void setBladeForces(double* bldForce, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=6);
  //! Set the tower force array 'twrForce' of turbine number 'iTurbGlob' at time step 't'
  void setTowerForces(double* twrForce, int iTurbGlob, fast::timeStep t = fast::STATE_NP1, int nSize=6);


  //! Get all turbine parametric data
  void get_turbineParams(int iTurbGlob, turbineDataType & turbData);
  //! Get the starting time step of the simulation. Safe to call from every MPI rank.
  int get_ntStart() { return ntStart; }
  //! Return a boolean flag whether the simulation is dryRun. Safe to call from every MPI rank.
  bool isDryRun() { return dryRun; }
  //! Return a boolean flag whether the simulation is debug. Safe to call from every MPI rank.
  bool isDebug() { return debug; }
  //! Get an enum of type 'simStartType' indicating the start type of the simulation. Safe to call from every MPI rank.
  simStartType get_simStartType() { return simStart; }
  //! Is the simulation time zero right now? Safe to call from every MPI rank.
  bool isTimeZero() { return timeZero; }
  //! Get the global MPI rank of the processor containing turbine number 'iTurbGlob'. Safe to call from every MPI rank.
  int get_procNo(int iTurbGlob) { return turbineMapGlobToProc[iTurbGlob] ; }
  //! Get the local turbine number of the turbine number 'iTurbGlob'. Safe to call from every MPI rank.
  int get_localTurbNo(int iTurbGlob) { return reverseTurbineMapProcToGlob[iTurbGlob]; }
  //! Get the total number of turbines in the simulation. Safe to call from every MPI rank.
  int get_nTurbinesGlob() { return nTurbinesGlob; }

  //! Get the nacelle area of turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  float get_nacelleArea(int iTurbGlob) { return get_nacelleAreaLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the nacelle drag coefficient of turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  float get_nacelleCd(int iTurbGlob) { return get_nacelleCdLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the air density around turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  float get_airDensity(int iTurbGlob) { return get_airDensityLoc(get_localTurbNo(iTurbGlob)); }

  //! Get the number of blades in turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  int get_numBlades(int iTurbGlob) { return get_numBladesLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Aerodyn/velocity nodes on each blade in turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  int get_numVelPtsBlade(int iTurbGlob) { return get_numVelPtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Aerodyn/velocity nodes on the tower in turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  int get_numVelPtsTwr(int iTurbGlob) { return get_numVelPtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the total number of Aerodyn/velocity nodes in turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  int get_numVelPts(int iTurbGlob) { return get_numVelPtsLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Actuator/force nodes on each blade in turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  int get_numForcePtsBlade(int iTurbGlob) { return get_numForcePtsBladeLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the number of Actuator/force nodes on the tower in turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  int get_numForcePtsTwr(int iTurbGlob) { return get_numForcePtsTwrLoc(get_localTurbNo(iTurbGlob)); }
  //! Get the total number of Actuator/force nodes in turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  int get_numForcePts(int iTurbGlob) { return get_numForcePtsLoc(get_localTurbNo(iTurbGlob)); }

  //! Compute the torque and thrust for turbine number 'iTurbGlob'. Must be called only from processor containing the turbine.
  void computeTorqueThrust(int iTurGlob, double*  torque, double*  thrust, int nSize);

  inline
  void getApproxHubPos(std::vector<double>& currentCoords, int iTurbGlob) {
      getApproxHubPos(currentCoords.data(), iTurbGlob, currentCoords.size());
  }

  inline
  void getHubPos(std::vector<double>& currentCoords, int iTurbGlob, fast::timeStep t = fast::STATE_NP1) {
      getHubPos(currentCoords.data(), iTurbGlob, t, currentCoords.size());
  }


  inline
  void getHubShftDir(std::vector<double> & hubShftVec, int iTurbGlob, fast::timeStep t = fast::STATE_NP1) {
      getHubShftDir(hubShftVec.data(), iTurbGlob,  t, hubShftVec.size());
  }

  inline
  void getVelNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1) {
      getVelNodeCoordinates(currentCoords.data(), iNode, iTurbGlob, t, currentCoords.size());
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
  void getForceNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1) {
      getForceNodeCoordinates(currentCoords.data(), iNode, iTurbGlob,  t, currentCoords.size());
  }

  inline
  void getForceNodeOrientation(std::vector<double> & currentOrientation, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1) {
      getForceNodeOrientation(currentOrientation.data(), iNode, iTurbGlob,  t, currentOrientation.size());
  }

  inline
  void getForce(std::vector<double> & currentForce, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1) {
      getForce(currentForce.data(), iNode, iTurbGlob,  t, currentForce.size());
  }

  inline
  void getRelativeVelForceNode(std::vector<double> & currentVelocity, int iNode, int iTurbGlob, fast::timeStep t = fast::STATE_NP1) {
      getRelativeVelForceNode(currentVelocity.data(), iNode, iTurbGlob,  t, currentVelocity.size());
  }

 inline
 void getBladeRefPositions(std::vector<double> & bldRefPos, int iTurbGlob){
    getBladeRefPositions(bldRefPos.data(), nTurbinesGlob);
 }
 inline
 void getBladeDisplacements(std::vector<double> & bldDefl, std::vector<double> & bldVel, int iTurbGlob, fast::timeStep t = fast::STATE_NP1)
 {
   getBladeDisplacements(bldDefl.data(), bldVel.data(), iTurbGlob, t, bldDefl.size());
 }
 inline
 void getBladeRootRefPositions(std::vector<double> & bldRootRefPos, int iTurbGlob){
     getBladeRootRefPositions(bldRootRefPos.data(), iTurbGlob);
 }
 void getBladeRootDisplacements(std::vector<double> & bldRootDefl, int iTurbGlob, fast::timeStep t = fast::STATE_NP1)
 {
   getBladeRootDisplacements(bldRootDefl.data(), iTurbGlob, t, bldRootDefl.size());
 }
 inline
 void getBladePitch(std::vector<double> & bldPitch, int iTurbGlob)
 {
     getBladePitch(bldPitch.data(), iTurbGlob, bldPitch.size());
 }
 inline
 void getTowerRefPositions(std::vector<double> & twrRefPos, int iTurbGlob)
 {
     getTowerRefPositions(twrRefPos.data(), iTurbGlob, 6);
 }
 inline
 void getTowerDisplacements(std::vector<double> & twrDefl, std::vector<double> & twrVel, int iTurbGlob, fast::timeStep t = fast::STATE_NP1)
 {
   getTowerDisplacements(twrDefl.data(), twrVel.data(), iTurbGlob, t, twrDefl.size());
 }
 inline
 void getHubRefPosition(std::vector<double> & hubRefPos, int iTurbGlob)
 {
   getHubRefPosition(hubRefPos.data(), iTurbGlob, hubRefPos.size());
 }
 inline
 void getHubDisplacement(std::vector<double> & hubDefl, std::vector<double> & hubVel, int iTurbGlob, fast::timeStep t = fast::STATE_NP1)
 {
   getHubDisplacement(hubDefl.data(), hubVel.data(), iTurbGlob, t, hubDefl.size());
 }
 inline
 void getNacelleRefPosition(std::vector<double> & nacRefPos, int iTurbGlob)
 {
   getNacelleRefPosition(nacRefPos.data(), iTurbGlob, nacRefPos.size());
 }
 inline
 void getNacelleDisplacement(std::vector<double> & nacDefl, std::vector<double> & nacVel, int iTurbGlob, fast::timeStep t = fast::STATE_NP1)
 {
   getNacelleDisplacement(nacDefl.data(), nacVel.data(), iTurbGlob, t, nacDefl.size());
 }

 inline
 void setBladeForces(std::vector<double> & bldForce, int iTurbGlob, fast::timeStep t = fast::STATE_NP1)
 {
   setBladeForces(bldForce.data(), iTurbGlob, t, 6);
 }
 inline
 void setTowerForces(std::vector<double> & twrForce, int iTurbGlob, fast::timeStep t = fast::STATE_NP1)
 {
   setTowerForces(twrForce.data(), iTurbGlob, t, 6);
 }
 inline
 void computeTorqueThrust(int iTurbGlob, std::vector<double> & torque, std::vector<double> & thrust)
 {
   computeTorqueThrust(iTurbGlob, torque.data(), thrust.data(), torque.size());
 }

  //! An example function to set velocities at the Aerodyn nodes using a power law wind profile using an exponent of 0.2 and a reference wind speed of 10 m/s at 90 meters. Safe to call from every MPI rank.
    void setExpLawWindSpeed(double t) ; // An example to set velocities at the Aerodyn nodes

  //! An example function to set a uniform X force at all blade nodes. Safe to call from every MPI rank.
  void setUniformXBladeForces(double loadX);


private:

  //! Set state from another state
  void set_state_from_state(fast::timeStep fromState, fast::timeStep toState);

  //! Preprare the C+++ output file for a new OpenFAST simulation
  void prepareOutputFile(int iTurbLoc);
  //! Find the C++ output file for a restarted simulation
  void findOutputFile(int iTurbLoc);
  //! Write output data to file
  void writeOutputFile(int iTurbLoc, int n_t_global);

  //! Find the OpenFAST restart file from the C++ restart file for a restarted simulation
  void findRestartFile(int iTurbLoc);
  //! Preprare the C+++ restart file for a new OpenFAST simulation
  void prepareRestartFile(int iTurbLoc);

  //! Read velocity and force node data at time step 'n', 'n-1' and 'n-2' to allow for a clean restart
  void readRestartFile(int iTurbLoc, int n_t_global);
  //! Write velocity and force node data at time step 'n', 'n-1' and 'n-2' to allow for a clean restart
  void writeRestartFile(int iTurbLoc, int n_t_global);

  //! Create velocity data file in preparation to write velocity data
  void prepareVelocityDataFile(int iTurb);
  //! Open velocity data file to read velocity data
  int openVelocityDataFile(int iTurb);
  //! Read the number of nonlinear iterations for a given driver time step
  int read_nlin_iters(int iTurb, int iTimestep, int ncid);
  //! Read velocity data at the Aerodyn nodes from velocity data file
  void readVelocityData(int iTurb, int iTimestep, int iNlinIter, int ncid);
  //! Write velocity data at the Aerodyn nodes from velocity data file
  void writeVelocityData(int iTurb, int iTimestep, int nlinIter);

  //! Check whether the error status is ok. If not quit gracefully by printing the error message
  void checkError(const int ErrStat, const char * ErrMsg);
  //! Check whether a file with name "name" exists
  inline bool checkFileExists(const std::string& name);

  //! Allocate memory for data structures for all turbines on this processor
  void allocateMemory_preInit();
  //! Allocate more memory for each turbine after intialization/restart
  void allocateMemory_postInit(int iTurbLoc);

  //! Get the nacelle drag coefficient of local turbine number 'iTurbLoc'
  float get_nacelleCdLoc(int iTurbLoc) { return turbineData[iTurbLoc].nacelle_cd; }
  //! Get the nacelle area of local turbine number 'iTurbLoc'
  float get_nacelleAreaLoc(int iTurbLoc) { return turbineData[iTurbLoc].nacelle_area; }
  //! Get the air density around local turbine number 'iTurbLoc'
  float get_airDensityLoc(int iTurbLoc) { return turbineData[iTurbLoc].air_density; }

  //! Get the number of blades in local turbine number 'iTurbLoc'
  int get_numBladesLoc(int iTurbLoc) { return turbineData[iTurbLoc].numBlades; }
  //! Get the number of Aerodyn/velocity nodes on each blade in local turbine number 'iTurbLoc'
  int get_numVelPtsBladeLoc(int iTurbLoc) { return turbineData[iTurbLoc].numVelPtsBlade; }
  //! Get the number of Aerodyn/velocity nodes on the tower in local turbine number 'iTurbLoc'
  int get_numVelPtsTwrLoc(int iTurbLoc) { return turbineData[iTurbLoc].numVelPtsTwr; }
  //! Get the total number of Aerodyn/velocity nodes in local turbine number 'iTurbLoc'
  int get_numVelPtsLoc(int iTurbLoc) { return turbineData[iTurbLoc].numVelPts; }
  //! Get the number of Actuator/force nodes on each blade in local turbine number 'iTurbLoc'
  int get_numForcePtsBladeLoc(int iTurbLoc) { return turbineData[iTurbLoc].numForcePtsBlade; }
  //! Get the number of Actuator/force nodes on the tower in local turbine number 'iTurbLoc'
  int get_numForcePtsTwrLoc(int iTurbLoc) { return turbineData[iTurbLoc].numForcePtsTwr; }
  //! Get the total number of Actuator/force nodes in local turbine number 'iTurbLoc'
  int get_numForcePtsLoc(int iTurbLoc) { return turbineData[iTurbLoc].numForcePts; }

  //! Get reference positions of blade-resolved FSI nodes from OpenFAST
  void get_ref_positions_from_openfast(int iTurb);

  void loadSuperController(const fastInputs & fi);

  //! Apply the velocity data at the Aerodyn nodes in 'velData' to turbine number 'iTurb' at time step 'iPrestart' through the data structure 'cDriver_Output_to_FAST'
  void applyVelocityData(int iPrestart, int iTurb, ExtInfw_OutputType_t o_t_FAST, std::vector<double> & velData) ;

  //! Compute cross product a x b and store it into aCrossB
  void cross(double * a, double * b, double * aCrossB);
  //! Apply a Wiener-Milenkovic rotation 'wm' to a vector 'r' into 'rRot'. To optionally transpose the rotation, set 'tranpose=-1.0'.
  void applyWMrotation(double * wm, double * r, double *rRot, double transpose = 1.0);
  //! Apply a Direction Cosine Matrix rotation 'dcm' to a vector 'r' into 'rRot'. To optionally transpose the rotation, set 'tranpose=-1.0'.
  void applyDCMrotation(double * dcm, double * r, double *rRot, double transpose = 1.0);

};

}

#endif
