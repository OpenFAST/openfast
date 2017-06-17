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
  std::vector<double> globStates;
  std::vector<double> globStates_np1;

  int nTurbineStates; // States for each turbine
  std::vector<double> turbineStates ;
  std::vector<double> turbineStates_np1 ;

  double d2R = 0.01745329251 ; //Degrees to Radians

  //Supercontroller stuff
  std::string scLibFile;
  // Dynamic load stuff copied from 'C++ dlopen mini HOWTO' on tldp.org
  void *scLibHandle ; 
  typedef void sc_updateStates_t(int nTurbines, int nInputsGlob, std::vector<double> & sc_inputsGlob, int nInputsTurbine, std::vector<double> & sc_inputsTurbine, int nGlobStates, std::vector<double> & globStates_n, std::vector<double> & globStates_np1, int nTurbineStates, std::vector<double> & turbineStates_n, std::vector<double> & turbineStates_np1); 
  sc_updateStates_t * sc_updateStates;
  typedef void sc_calcOutputs_t(int nTurbines, int nInputsGlob, std::vector<double> & sc_inputsGlob, int nInputsTurbine, std::vector<double> & sc_inputsTurbine, int nGlobStates, std::vector<double> & globStates, int nTurbineStates, std::vector<double> & turbineStates, int nOutputsGlob, std::vector<double> & sc_outputsGlob, int nOutputsTurbine, std::vector<double> & sc_outputsTurbine); 
  sc_calcOutputs_t * sc_calcOutputs;


 public:

  SuperController();

  ~SuperController() ;

  void init(int nTurbinesGlob);
  
  void load(scInData sci);

  void calcOutputs(std::vector<double> & sc_inputsGlob, std::vector<double> & sc_inputsTurbine, std::vector<double> & sc_outputsGlob, std::vector<double> & sc_outputsTurbine) ;

  void updateStates(std::vector<double> & sc_inputsGlob, std::vector<double> & sc_inputsTurbine) ; //Make a prediction for 'n+1' based on 'n'

  void advanceStates() ; //Advance states to time step 'n+1'

  int writeRestartFile(int n_t_global);

  int readRestartFile(int n_t_global);

  void end() {} ;
};

