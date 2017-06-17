#include "SC.h"

SuperController::SuperController():
nInputsTurbine(0),
nOutputsTurbine(0),
nInputsGlob(0),
nOutputsGlob(0),
nGlobStates(0),
nTurbineStates(0)
{
  
}

SuperController::~SuperController() {

    
  if(scLibHandle != NULL) {
      // close the library
      std::cout << "Closing SC library..." << std::endl;
      dlclose(scLibHandle);
  }
  
}

void SuperController::load(scInData sci)  {

    nInputsTurbine = sci.nInputsTurbine; 
    nOutputsTurbine = sci.nOutputsTurbine;
    nInputsGlob = sci.nInputsGlob;
    nOutputsGlob = sci.nOutputsGlob;
    nGlobStates = sci.nGlobStates;
    nTurbineStates = sci.nTurbineStates;
    scLibFile = sci.scLibFile;


    if (nInputsGlob < 0) 
        std::cerr << "Supercontroller: nInputsGlob has to be greater than zero" << std::endl ;
    
    if (nOutputsGlob < 0) 
        std::cerr << "Supercontroller: nOutputsGlob is less than zero." << std::endl ;

    if (nGlobStates < 0) 
        std::cerr << "Supercontroller: nGlobStates has to be greater than zero" << std::endl ;

    if (nTurbineStates < 0)
        std::cerr << "Supercontroller: nTurbineStates has to be greater than zero" << std::endl ;
    
    // open the library
    scLibHandle = dlopen(scLibFile.c_str(), RTLD_LAZY);
    if (!scLibHandle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';
    }
        
    sc_updateStates = (sc_updateStates_t*) dlsym(scLibHandle, "sc_updateStates");
    // reset errors
    const char *dlsym_error_us = dlerror();
    if (dlsym_error_us) {
        std::cerr << "Cannot load symbol 'sc_updateStates': " << dlsym_error_us << '\n';
        dlclose(scLibHandle);
    }
    
    sc_calcOutputs = (sc_calcOutputs_t*) dlsym(scLibHandle, "sc_calcOutputs");
    // reset errors
    const char *dlsym_error_co = dlerror();
    if (dlsym_error_co) {
        std::cerr << "Cannot load symbol 'sc_calcOutputs': " << dlsym_error_co << '\n';
        dlclose(scLibHandle);
    }

}

void SuperController::init(int nTurbinesGlob) {
  
  nTurbines = nTurbinesGlob;

  globStates.resize(nGlobStates); 
  globStates_np1.resize(nGlobStates); 

  turbineStates.resize(nTurbines*nTurbineStates);
  turbineStates_np1.resize(nTurbines*nTurbineStates);
  
  // Initialize the turbine states at time zero - Not sure how to do this. May be call calcOut?
 
}

void SuperController::calcOutputs(std::vector<double> & sc_inputsGlob, std::vector<double> & sc_inputsTurbine, std::vector<double> & sc_outputsGlob, std::vector<double> & sc_outputsTurbine) {

    sc_calcOutputs(nTurbines, nInputsGlob, sc_inputsGlob, nInputsTurbine, sc_inputsTurbine, nGlobStates, globStates, nTurbineStates, turbineStates, nOutputsGlob, sc_outputsGlob, nOutputsTurbine, sc_outputsTurbine);   

}

void SuperController::advanceStates() {

  for(int iState=0; iState<nGlobStates; iState++) {
    globStates[iState] = globStates_np1[iState];
  }
  
  for(int iTurb=0; iTurb < nTurbines; iTurb++) {
    for(int iState=0; iState < nTurbineStates; iState++) {
        turbineStates[iTurb*nTurbineStates + iState] = turbineStates_np1[iTurb*nTurbineStates + iState];
    }
  }

}

void SuperController::updateStates(std::vector<double> & sc_inputsGlob, std::vector<double> & sc_inputsTurbine) {


  // Meaning of scInputs
  // 0 - Time
  // 1 - GenTorque

  // Meaning of scOutputs
  // 0 - Minimum Blade pitch
  

  // Turbine 0
    /*  Vary PC_MinPit as a function of time: */
    /*  0-20s: 0 degrees */
    /*  20-40s: 1.5 degrees */
    /*  40-60s: 3 degrees */

  // Turbine 1
    /*  Vary PC_MinPit as a function of time: */
    /*  0-20s: 0.5 degrees */
    /*  20-40s: 1 degrees */
    /*  40-60s: 2.5 degrees */

    sc_updateStates(nTurbines, nInputsGlob, sc_inputsGlob, nInputsTurbine, sc_inputsTurbine, nGlobStates, globStates, globStates_np1, nTurbineStates, turbineStates, turbineStates_np1) ;

/*   //Copy inputs into states first */

    }

int SuperController::readRestartFile(int n_t_global) {

  hid_t restartFile = H5Fopen(("sc" + std::to_string(n_t_global) + ".chkp.h5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  
  {
    hid_t attr = H5Aopen(restartFile, "nTurbines", H5P_DEFAULT);
    herr_t ret = H5Aread(attr, H5T_NATIVE_INT, &nTurbines) ;
    H5Aclose(attr);

    attr = H5Aopen(restartFile, "nInputsTurbine", H5P_DEFAULT);
    ret = H5Aread(attr, H5T_NATIVE_INT, &nInputsTurbine) ;
    H5Aclose(attr);

    attr = H5Aopen(restartFile, "nOutputsTurbine", H5P_DEFAULT);
    ret = H5Aread(attr, H5T_NATIVE_INT, &nOutputsTurbine) ;
    H5Aclose(attr);

    attr = H5Aopen(restartFile, "nInputsGlob", H5P_DEFAULT);
    ret = H5Aread(attr, H5T_NATIVE_INT, &nInputsGlob) ;
    H5Aclose(attr);

    attr = H5Aopen(restartFile, "nOutputsGlob", H5P_DEFAULT);
    ret = H5Aread(attr, H5T_NATIVE_INT, &nOutputsGlob) ;
    H5Aclose(attr);

    attr = H5Aopen(restartFile, "nGlobStates", H5P_DEFAULT);
    ret = H5Aread(attr, H5T_NATIVE_INT, &nGlobStates) ;
    H5Aclose(attr);

    globStates.resize(nGlobStates); 
    globStates_np1.resize(nGlobStates); 

    attr = H5Aopen(restartFile, "nTurbineStates", H5P_DEFAULT);
    ret = H5Aread(attr, H5T_NATIVE_INT, &nTurbineStates) ;
    H5Aclose(attr);

    turbineStates.resize(nTurbines*nTurbineStates);
    turbineStates_np1.resize(nTurbines*nTurbineStates);

#ifdef DEBUG
    std::cout << "nTurbines = " << nTurbines << std::endl ;
    std::cout << "nInputsTurbine = " << nInputsTurbine << std::endl ;
    std::cout << "nOutputsTurbine = " << nOutputsTurbine << std::endl ;
    std::cout << "nInputsGlob = " << nInputsGlob << std::endl ;
    std::cout << "nOutputsGlob = " << nOutputsGlob << std::endl ;
    std::cout << "nGlobStates = " << nGlobStates << std::endl ;
    std::cout << "nTurbineStates = " << nTurbineStates << std::endl ;
#endif
   
  }

  if (nGlobStates > 0) {
    hid_t dataSet = H5Dopen2(restartFile, "/globStates", H5P_DEFAULT);
    herr_t status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, globStates.data());
    status = H5Dclose(dataSet);

    dataSet = H5Dopen2(restartFile, "/globStates_np1", H5P_DEFAULT);
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, globStates_np1.data());
    status = H5Dclose(dataSet);
  }
  
  if (nTurbineStates > 0) {
    hid_t dataSet = H5Dopen2(restartFile, "turbineStates", H5P_DEFAULT);
    herr_t status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, turbineStates.data());
    status = H5Dclose(dataSet);

    dataSet = H5Dopen2(restartFile, "turbineStates_np1", H5P_DEFAULT);
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, turbineStates_np1.data());
    status = H5Dclose(dataSet);
  }

#ifdef DEBUG
  for(int iTurb=0; iTurb < nTurbines; iTurb++) {
    for(int i=0; i < nTurbineStates; i++) {
      std::cout << "iTurb = " << iTurb << ", i = " << i << ",  " ;
      std::cout << turbineStates[iTurb*nTurbineStates + i] << std::endl ;
    }
  }
#endif
  herr_t status = H5Fclose(restartFile);
  
}


int SuperController::writeRestartFile(int n_t_global) {

  /* // HDF5 stuff to write states to restart file or read back from it */

  hid_t restartFile = H5Fcreate(("sc" + std::to_string(n_t_global) + ".chkp.h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  {
    hsize_t dims[1];
    dims[0] = 1;
    hid_t dataSpace = H5Screate_simple(1, dims, NULL);
    hid_t attr = H5Acreate2(restartFile, "nTurbines", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
    herr_t status = H5Awrite(attr, H5T_NATIVE_INT, &nTurbines);
    status = H5Aclose(attr);
    status = H5Sclose(dataSpace);

    dataSpace = H5Screate_simple(1, dims, NULL);
    attr = H5Acreate2(restartFile, "nInputsTurbine", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
    status = H5Awrite(attr, H5T_NATIVE_INT, &nInputsTurbine);
    status = H5Aclose(attr);
    status = H5Sclose(dataSpace);

    dataSpace = H5Screate_simple(1, dims, NULL);
    attr = H5Acreate2(restartFile, "nOutputsTurbine", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
    status = H5Awrite(attr, H5T_NATIVE_INT, &nOutputsTurbine);
    status = H5Aclose(attr);
    status = H5Sclose(dataSpace);

    dataSpace = H5Screate_simple(1, dims, NULL);
    attr = H5Acreate2(restartFile, "nInputsGlob", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
    status = H5Awrite(attr, H5T_NATIVE_INT, &nInputsGlob);
    status = H5Aclose(attr);
    status = H5Sclose(dataSpace);

    dataSpace = H5Screate_simple(1, dims, NULL);
    attr = H5Acreate2(restartFile, "nOutputsGlob", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
    status = H5Awrite(attr, H5T_NATIVE_INT, &nOutputsGlob);
    status = H5Aclose(attr);
    status = H5Sclose(dataSpace);

    dataSpace = H5Screate_simple(1, dims, NULL);
    attr = H5Acreate2(restartFile, "nGlobStates", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
    status = H5Awrite(attr, H5T_NATIVE_INT, &nGlobStates);
    status = H5Aclose(attr);
    status = H5Sclose(dataSpace);

    dataSpace = H5Screate_simple(1, dims, NULL);
    attr = H5Acreate2(restartFile, "nTurbineStates", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
    status = H5Awrite(attr, H5T_NATIVE_INT, &nTurbineStates);
    status = H5Aclose(attr);
    status = H5Sclose(dataSpace);
    
  }

  if (nGlobStates > 0) {
    hsize_t dims[1];
    dims[0] = nGlobStates;
    hid_t dataSpace = H5Screate_simple(1, dims, NULL);
    hid_t dataSet = H5Dcreate2(restartFile, "/globStates", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
    herr_t status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, globStates.data());
    
    status = H5Dclose(dataSet);
    status = H5Sclose(dataSpace);

    dataSpace = H5Screate_simple(1, dims, NULL);
    dataSet = H5Dcreate2(restartFile, "/globStates_np1", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, globStates_np1.data());
    
    status = H5Dclose(dataSet);
    status = H5Sclose(dataSpace);


  }
  
  if (nTurbineStates > 0) {

    hsize_t dims[2];
    dims[0] = nTurbines;
    dims[1] = nTurbineStates;

    hid_t dataSpace = H5Screate_simple(2, dims, NULL);
    hid_t dataSet = H5Dcreate2(restartFile, "turbineStates", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
    herr_t status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, turbineStates.data());
    
    status = H5Dclose(dataSet);
    status = H5Sclose(dataSpace);

    dataSpace = H5Screate_simple(2, dims, NULL);
    dataSet = H5Dcreate2(restartFile, "turbineStates_np1", H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);    
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, turbineStates_np1.data());
    
    status = H5Dclose(dataSet);
    status = H5Sclose(dataSpace);
  }

  herr_t status = H5Fclose(restartFile);

  return 0;

}


