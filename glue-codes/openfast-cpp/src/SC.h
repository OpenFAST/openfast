#ifndef OMPI_SKIP_MPICXX
 #define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
 #define MPICH_SKIP_MPICXX
#endif
#include "FAST_Library.h"
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "mpi.h"
#include "hdf5.h"
#include "dlfcn.h"

class scInitOutData {

public:
    int nInpGlobal;
    int nCtrl2SC;
    int nSC2CtrlGlob;
    int nSC2Ctrl;
    std::vector<float> from_SCglob;
    std::vector<std::vector<float>> from_SC;
};

class SuperController {

public:

    // Data structures to interface with OpenFAST per turbine
    // Unfortunately have to be public
    std::vector<SC_DX_InputType_t> ip_from_FAST; // At time step 'n+1'
    std::vector<SC_DX_OutputType_t> op_to_FAST;  // At time step 'n'

private:

    MPI_Comm  fastMPIComm;

    int nTurbinesGlob;
    int nTurbinesProc;
    std::map<int, int> turbineMapProcToGlob;

    int nCtrl2SC;
    int nSC2Ctrl;
    int nInpGlobal;
    int nSC2CtrlGlob;

    int nStatesGlobal; // Global states like time
    std::vector<float> globStates;
    std::vector<float> globStates_np1;

    int nStatesTurbine; // States for each turbine
    std::vector<float> turbineStates ;
    std::vector<float> turbineStates_np1 ;

    // Time 'n-1'
    std::vector<float> from_SC_nm1;  // # outputs from the supercontroller for turbines
    std::vector<float> to_SC_nm1;   // # inputs to the supercontroller from turbines
    std::vector<float> from_SCglob_nm1;  // # outputs from the supercontroller for glob
    std::vector<float> to_SCglob_nm1;   // # inputs to the supercontroller from glob
    // Time 'n'
    std::vector<float> from_SC_n;  // # outputs from the supercontroller for turbines
    std::vector<float> to_SC_n;   // # inputs to the supercontroller from turbines
    std::vector<float> from_SCglob_n;  // # outputs from the supercontroller for glob
    std::vector<float> to_SCglob_n;   // # inputs to the supercontroller from glob
    // Time 'n+1'
    std::vector<float> from_SC_np1;  // # outputs from the supercontroller for turbines
    std::vector<float> to_SC_np1;   // # inputs to the supercontroller from turbines
    std::vector<float> from_SCglob_np1;  // # outputs from the supercontroller for glob
    std::vector<float> to_SCglob_np1;   // # inputs to the supercontroller from glob

    int nParamGlobal;
    std::vector<float> paramGlobal;
    int nParamTurbine;
    std::vector<float> paramTurbine;

    int ErrStat;
    char ErrMsg[INTERFACE_STRING_LENGTH];  // make sure this is the same size as IntfStrLen in FAST_Library.f90

    float d2R = 0.01745329251 ; //Degrees to Radians

    //Supercontroller stuff
    std::string scLibFile;
    // Dynamic load stuff copied from 'C++ dlopen mini HOWTO' on tldp.org
    void *scLibHandle ;
    typedef void sc_init_t(int * nTurbinesGlob, int * nInpGlobal, int * nCtrl2SC, int * nParamGlobal, int * nParamTurbine, int * nStatesGlobal, int * nStatesTurbine, int * nSC2CtrlGlob, int * nSC2Ctrl, int *ErrStat, char * ErrMsg);
    sc_init_t * sc_init;
    bool sc_library_loaded = false;

    typedef void sc_getInitData_t(int * nTurbinesGlob, int * nParamGlobal, int * nParamTurbine, float * paramGlobal, float * paramTurbine, int * nSC2CtrlGlob, float * from_SCglob, int * nSC2Ctrl, float * from_SC, int * nStatesGlobal, float * globStates, int * nStatesTurbine, float * turbineStates, int *ErrStat, char * ErrMsg);
    sc_getInitData_t * sc_getInitData;

    typedef void sc_updateStates_t(double * t, int * nTurbinesGlob, int * nParamGlobal, float * paramGlobal, int * nParamTurbine, float * paramTurbine, int * nInpGlobal, float * to_SCglob, int * nCtrl2SC, float * to_SC, int * nStatesGlobal, float * statesGlob_n, float * statesGlob_np1, int * nStatesTurbine, float * statesTurbine_n, float * statesTurbine_np1, int * ErrStat, char * ErrMsg);
    sc_updateStates_t * sc_updateStates;

    typedef void sc_calcOutputs_t(double *  t, int * nTurbinesGlob, int * nParamGlobal, float * paramGlobal, int * nParamTurbine, float * paramTurbine, int * nInpGlobal, float * to_SCglob, int * nCtrl2SC, float * to_SC, int * nStatesGlobal, float * statesGlob, int * nStatesTurbine, float * statesTurbine, int * nSC2CtrlGlob, float * from_SCglob, int * nSC2Ctrl, float * from_SC, int * ErrStat, char * ErrMsg);
    sc_calcOutputs_t * sc_calcOutputs;


public:

    SuperController();

    ~SuperController() ;

    void init(scInitOutData & scio, int nTurbinesProc);
    void init_sc(scInitOutData & scio, int inNTurbinesProc, std::map<int, int> iTurbineMapProcToGlob, MPI_Comm inFastMPIComm);

    void load(int inNTurbinesGlob, std::string inScLibFile, scInitOutData & scio);

    void updateStates(double t) ; //Make a prediction for states at 'n+1' based on inputs and states at 'n'

    void calcOutputs_n(double t) ;
    void calcOutputs_np1(double t) ;

    void fastSCInputOutput() ; // Exchange input output information with OpenFAST turbines

    void advanceTime() ; //Advance states to time step 'n+1'

    int writeRestartFile(int n_t_global);

    int readRestartFile(int n_t_global);

    void end() {} ;
};
