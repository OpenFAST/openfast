
#include "SuperController_Types.h"
#include <sstream>
#include <iostream>
#include "hdf5.h"
#include <string>

class SuperController {

 private:
  
  int nTurbines;
  int nScInputs;
  int nScOutputs;

  int nGlobStates; // Global states like time 
  double * globStates;

  int nTurbineStates; // States for each turbine
  double ** turbineStates ;

  double d2R = 0.01745329251 ; //Degrees to Radians

 public:
  
  SuperController();

  virtual ~SuperController() ;

  virtual void init(int n, int numScInputs, int numScOutputs);
  
  virtual void calcOutputs(std::vector<double> & scOutputsGlob) ;

  virtual void updateStates(std::vector<double> & scInputsGlob) ;

  virtual int writeRestartFile(int n_t_global);

  virtual int readRestartFile(int n_t_global);
};
