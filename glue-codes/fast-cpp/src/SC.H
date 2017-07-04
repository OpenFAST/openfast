#ifndef OMPI_SKIP_MPICXX
 #define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
 #define MPICH_SKIP_MPICXX
#endif

#include <sstream>
#include <iostream>
#include <string>
#include "yaml-cpp/yaml.h"
#include "hdf5.h"
#include "dlfcn.h"

class scInData {

public: 
    int nInputsTurbine;
    int nOutputsTurbine;
    int nInputsGlob;
    int nOutputsGlob;
    int nGlobStates;
    int nTurbineStates;
    std::string scLibFile;

};

class SuperController {

 private:
  
  int nTurbines;
  int nInputsTurbine;
  int nOutputsTurbine;
  int nInputsGlob;
  int nOutputsGlob;

  int nGlobStates; // Global states like time 
  std::vector<float> globStates;
  std::vector<float> globStates_np1;

  int nTurbineStates; // States for each turbine
  std::vector<float> turbineStates ;
  std::vector<float> turbineStates_np1 ;

  float d2R = 0.01745329251 ; //Degrees to Radians

  //Supercontroller stuff
  std::string scLibFile;
  // Dynamic load stuff copied from 'C++ dlopen mini HOWTO' on tldp.org
  void *scLibHandle ; 
  typedef void sc_updateStates_t(int nTurbines, int nInputsGlob, std::vector<float> & sc_inputsGlob, int nInputsTurbine, std::vector<float> & sc_inputsTurbine, int nGlobStates, std::vector<float> & globStates_n, std::vector<float> & globStates_np1, int nTurbineStates, std::vector<float> & turbineStates_n, std::vector<float> & turbineStates_np1); 
  sc_updateStates_t * sc_updateStates;
  typedef void sc_calcOutputs_t(int nTurbines, int nInputsGlob, std::vector<float> & sc_inputsGlob, int nInputsTurbine, std::vector<float> & sc_inputsTurbine, int nGlobStates, std::vector<float> & globStates, int nTurbineStates, std::vector<float> & turbineStates, int nOutputsGlob, std::vector<float> & sc_outputsGlob, int nOutputsTurbine, std::vector<float> & sc_outputsTurbine); 
  sc_calcOutputs_t * sc_calcOutputs;


 public:

  SuperController();

  ~SuperController() ;

  void init(int nTurbinesGlob);
  
  void load(scInData sci);

  void calcOutputs(std::vector<float> & sc_inputsGlob, std::vector<float> & sc_inputsTurbine, std::vector<float> & sc_outputsGlob, std::vector<float> & sc_outputsTurbine) ;

  void updateStates(std::vector<float> & sc_inputsGlob, std::vector<float> & sc_inputsTurbine) ; //Make a prediction for 'n+1' based on 'n'

  void advanceStates() ; //Advance states to time step 'n+1'

  int writeRestartFile(int n_t_global);

  int readRestartFile(int n_t_global);

  void end() {} ;
};

