#include "SC.h"

SuperController::SuperController():
nCtrl2SC(0),
nSC2Ctrl(0),
nInpGlobal(0),
nSC2CtrlGlob(0),
nStatesGlobal(0),
nStatesTurbine(0)
{

}

SuperController::~SuperController() {
    // close the library
    if (sc_library_loaded) {
        std::cout << "Closing SC library..." << std::endl;
        dlclose(scLibHandle);
    }
}

void SuperController::load(int inNTurbinesGlob, std::string inScLibFile, scInitOutData & scio)  {

    nTurbinesGlob = inNTurbinesGlob;
    scLibFile = inScLibFile;

    // open the library
    scLibHandle = dlopen(scLibFile.c_str(), RTLD_LAZY);
    if (!scLibHandle) {
        std::cerr << "Cannot open library: " << dlerror() << '\n';
    }
    sc_library_loaded = true;

    sc_init = (sc_init_t*) dlsym(scLibHandle, "sc_init");
    // reset errors
    const char *dlsym_error_i = dlerror();
    if (dlsym_error_i) {
        std::cerr << "Cannot load symbol 'sc_init': " << dlsym_error_i << '\n';
        dlclose(scLibHandle);
    }

    sc_getInitData = (sc_getInitData_t*) dlsym(scLibHandle, "sc_getInitData");
    // reset errors
    const char *dlsym_error_gid = dlerror();
    if (dlsym_error_gid) {
        std::cerr << "Cannot load symbol 'sc_getInitData': " << dlsym_error_gid << '\n';
        dlclose(scLibHandle);
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

    sc_init(&nTurbinesGlob, &nInpGlobal, &nCtrl2SC, &nParamGlobal, &nParamTurbine, &nStatesGlobal, &nStatesTurbine, &nSC2CtrlGlob, &nSC2Ctrl, &ErrStat, ErrMsg);

    if (nInpGlobal != 0)
        std::cerr << "Supercontroller: nInpGlobal has to be zero. Not implemented yet." << std::endl ;

    if (nCtrl2SC < 0)
        std::cerr << "Supercontroller: nCtrl2SC is less than zero." << std::endl ;

    if (nParamGlobal < 0)
        std::cerr << "Supercontroller: nParamGlobal is less than zero." << std::endl ;

    if (nParamTurbine < 0)
        std::cerr << "Supercontroller: nParamTurbine is less than zero." << std::endl ;

    if (nStatesGlobal < 0)
        std::cerr << "Supercontroller: nStatesGlobal is less than zero" << std::endl ;

    if (nStatesTurbine < 0)
        std::cerr << "Supercontroller: nStatesTurbine is less than zero" << std::endl ;

    if (nSC2CtrlGlob < 0)
        std::cerr << "Supercontroller: nSC2CtrlGlob is less than zero." << std::endl ;

    if (nSC2Ctrl < 0)
        std::cerr << "Supercontroller: nSC2Ctrl is less than zero." << std::endl ;

    scio.nInpGlobal = nInpGlobal;
    scio.nCtrl2SC = nCtrl2SC;
    scio.nSC2Ctrl = nSC2Ctrl;
    scio.nSC2CtrlGlob = nSC2CtrlGlob;

}

void SuperController::init(scInitOutData & scio, int nTurbinesProc) {
    ip_from_FAST.resize(nTurbinesProc) ;
    op_to_FAST.resize(nTurbinesProc) ;

    scio.nSC2CtrlGlob = 0;
    scio.nSC2Ctrl = 0;
    scio.nCtrl2SC = 0;

    scio.from_SCglob.resize(nSC2CtrlGlob);
    scio.from_SC.resize(nTurbinesProc);
    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        scio.from_SC[iTurb].resize(nSC2Ctrl);
    }
}

void SuperController::init_sc(scInitOutData & scio, int inNTurbinesProc, std::map<int, int> iTurbineMapProcToGlob, MPI_Comm inFastMPIComm) {

    fastMPIComm = inFastMPIComm;
    nTurbinesProc = inNTurbinesProc;
    turbineMapProcToGlob = iTurbineMapProcToGlob;

    if (nTurbinesProc > 0) {

        paramGlobal.resize(nParamGlobal);
        paramTurbine.resize(nTurbinesGlob*nParamTurbine);

        globStates.resize(nStatesGlobal);
        globStates_np1.resize(nStatesGlobal);

        turbineStates.resize(nTurbinesGlob*nStatesTurbine);
        turbineStates_np1.resize(nTurbinesGlob*nStatesTurbine);

        from_SC_nm1.resize(nTurbinesGlob*nSC2Ctrl);
        from_SC_n.resize(nTurbinesGlob*nSC2Ctrl);
        from_SC_np1.resize(nTurbinesGlob*nSC2Ctrl);

        to_SC_nm1.resize(nTurbinesGlob*nCtrl2SC);
        to_SC_n.resize(nTurbinesGlob*nCtrl2SC);
        to_SC_np1.resize(nTurbinesGlob*nCtrl2SC);

        from_SCglob_nm1.resize(nTurbinesGlob*nSC2CtrlGlob);
        from_SCglob_n.resize(nTurbinesGlob*nSC2CtrlGlob);
        from_SCglob_np1.resize(nTurbinesGlob*nSC2CtrlGlob);

        to_SCglob_nm1.resize(nTurbinesGlob*nInpGlobal);
        to_SCglob_n.resize(nTurbinesGlob*nInpGlobal);
        to_SCglob_np1.resize(nTurbinesGlob*nInpGlobal);

        sc_getInitData(&nTurbinesGlob, &nParamGlobal, &nParamTurbine, paramGlobal.data(), paramTurbine.data(), &nSC2CtrlGlob, from_SCglob_nm1.data(), &nSC2Ctrl, from_SC_nm1.data(), &nStatesGlobal, globStates.data(), &nStatesTurbine, turbineStates.data(), &ErrStat, ErrMsg);

        for(int i=0; i < nSC2CtrlGlob; i++) {
            scio.from_SCglob[i] = from_SCglob_nm1[i];
        }

        for (int iTurb = 0 ; iTurb < nTurbinesProc; iTurb++) {
            for(int i=0; i < nSC2Ctrl; i++) {
                scio.from_SC[iTurb][i] = from_SC_nm1[turbineMapProcToGlob[iTurb]*nSC2Ctrl + i];
            }
        }

    }

}

void SuperController::calcOutputs_n(double t) {

    if (nTurbinesProc > 0) {
        sc_calcOutputs(&t, &nTurbinesGlob, &nParamGlobal, paramGlobal.data(), &nParamTurbine, paramTurbine.data(), &nInpGlobal, to_SCglob_n.data(), &nCtrl2SC, to_SC_n.data(), &nStatesGlobal, globStates.data(), &nStatesTurbine, turbineStates.data(), &nSC2CtrlGlob, from_SCglob_n.data(), &nSC2Ctrl, from_SC_n.data(), &ErrStat, ErrMsg);
    }

}

void SuperController::calcOutputs_np1(double t) {

    if (nTurbinesProc > 0) {
        sc_calcOutputs(&t, &nTurbinesGlob, &nParamGlobal, paramGlobal.data(), &nParamTurbine, paramTurbine.data(), &nInpGlobal, to_SCglob_n.data(), &nCtrl2SC, to_SC_n.data(), &nStatesGlobal, globStates_np1.data(), &nStatesTurbine, turbineStates_np1.data(), &nSC2CtrlGlob, from_SCglob_np1.data(), &nSC2Ctrl, from_SC_np1.data(), &ErrStat, ErrMsg);
    }

}

void SuperController::updateStates(double t) {

    if (nTurbinesProc > 0) {
        sc_updateStates(&t, &nTurbinesGlob, &nParamGlobal, paramGlobal.data(), &nParamTurbine, paramTurbine.data(), &nInpGlobal, to_SCglob_n.data(), &nCtrl2SC, to_SC_n.data(), &nStatesGlobal, globStates.data(), globStates_np1.data(), &nStatesTurbine, turbineStates.data(), turbineStates_np1.data(), &ErrStat, ErrMsg);
    }

}

int SuperController::readRestartFile(int n_t_global) {

    if (nTurbinesProc > 0) {

        hid_t restartFile = H5Fopen(("sc" + std::to_string(n_t_global) + ".chkp.h5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

        {
            hid_t attr = H5Aopen(restartFile, "nTurbinesGlob", H5P_DEFAULT);
            herr_t ret = H5Aread(attr, H5T_NATIVE_INT, &nTurbinesGlob) ;
            H5Aclose(attr);

            attr = H5Aopen(restartFile, "nCtrl2SC", H5P_DEFAULT);
            ret = H5Aread(attr, H5T_NATIVE_INT, &nCtrl2SC) ;
            H5Aclose(attr);

            attr = H5Aopen(restartFile, "nSC2Ctrl", H5P_DEFAULT);
            ret = H5Aread(attr, H5T_NATIVE_INT, &nSC2Ctrl) ;
            H5Aclose(attr);

            attr = H5Aopen(restartFile, "nInpGlobal", H5P_DEFAULT);
            ret = H5Aread(attr, H5T_NATIVE_INT, &nInpGlobal) ;
            H5Aclose(attr);

            attr = H5Aopen(restartFile, "nSC2CtrlGlob", H5P_DEFAULT);
            ret = H5Aread(attr, H5T_NATIVE_INT, &nSC2CtrlGlob) ;
            H5Aclose(attr);

            attr = H5Aopen(restartFile, "nStatesGlobal", H5P_DEFAULT);
            ret = H5Aread(attr, H5T_NATIVE_INT, &nStatesGlobal) ;
            H5Aclose(attr);

            globStates.resize(nStatesGlobal);
            globStates_np1.resize(nStatesGlobal);

            attr = H5Aopen(restartFile, "nStatesTurbine", H5P_DEFAULT);
            ret = H5Aread(attr, H5T_NATIVE_INT, &nStatesTurbine) ;
            H5Aclose(attr);

            turbineStates.resize(nTurbinesGlob*nStatesTurbine);
            turbineStates_np1.resize(nTurbinesGlob*nStatesTurbine);

#ifdef DEBUG
            std::cout << "nTurbinesGlob = " << nTurbinesGlob << std::endl ;
            std::cout << "nCtrl2SC = " << nCtrl2SC << std::endl ;
            std::cout << "nSC2Ctrl = " << nSC2Ctrl << std::endl ;
            std::cout << "nInpGlobal = " << nInpGlobal << std::endl ;
            std::cout << "nSC2CtrlGlob = " << nSC2CtrlGlob << std::endl ;
            std::cout << "nStatesGlobal = " << nStatesGlobal << std::endl ;
            std::cout << "nStatesTurbine = " << nStatesTurbine << std::endl ;
#endif

        }

        if (nStatesGlobal > 0) {
            hid_t dataSet = H5Dopen2(restartFile, "/globStates", H5P_DEFAULT);
            herr_t status = H5Dread(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, globStates.data());
            status = H5Dclose(dataSet);

            dataSet = H5Dopen2(restartFile, "/globStates_np1", H5P_DEFAULT);
            status = H5Dread(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, globStates_np1.data());
            status = H5Dclose(dataSet);
        }

        if (nStatesTurbine > 0) {
            hid_t dataSet = H5Dopen2(restartFile, "turbineStates", H5P_DEFAULT);
            herr_t status = H5Dread(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, turbineStates.data());
            status = H5Dclose(dataSet);

            dataSet = H5Dopen2(restartFile, "turbineStates_np1", H5P_DEFAULT);
            status = H5Dread(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, turbineStates_np1.data());
            status = H5Dclose(dataSet);
        }

#ifdef DEBUG
        for(int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
            for(int i=0; i < nStatesTurbine; i++) {
                std::cout << "iTurb = " << iTurb << ", i = " << i << ",  " ;
                std::cout << turbineStates[iTurb*nStatesTurbine + i] << std::endl ;
            }
        }
#endif
        herr_t status = H5Fclose(restartFile);
    }

    return 0;
}


int SuperController::writeRestartFile(int n_t_global) {

  /* // HDF5 stuff to write states to restart file or read back from it */

    if (nTurbinesProc > 0) {

        hid_t restartFile = H5Fcreate(("sc" + std::to_string(n_t_global) + ".chkp.h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        {
            hsize_t dims[1];
            dims[0] = 1;
            hid_t dataSpace = H5Screate_simple(1, dims, NULL);
            hid_t attr = H5Acreate2(restartFile, "nTurbinesGlob", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            herr_t status = H5Awrite(attr, H5T_NATIVE_INT, &nTurbinesGlob);
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);

            dataSpace = H5Screate_simple(1, dims, NULL);
            attr = H5Acreate2(restartFile, "nCtrl2SC", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            status = H5Awrite(attr, H5T_NATIVE_INT, &nCtrl2SC);
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);

            dataSpace = H5Screate_simple(1, dims, NULL);
            attr = H5Acreate2(restartFile, "nSC2Ctrl", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            status = H5Awrite(attr, H5T_NATIVE_INT, &nSC2Ctrl);
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);

            dataSpace = H5Screate_simple(1, dims, NULL);
            attr = H5Acreate2(restartFile, "nInpGlobal", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            status = H5Awrite(attr, H5T_NATIVE_INT, &nInpGlobal);
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);

            dataSpace = H5Screate_simple(1, dims, NULL);
            attr = H5Acreate2(restartFile, "nSC2CtrlGlob", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            status = H5Awrite(attr, H5T_NATIVE_INT, &nSC2CtrlGlob);
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);

            dataSpace = H5Screate_simple(1, dims, NULL);
            attr = H5Acreate2(restartFile, "nStatesGlobal", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            status = H5Awrite(attr, H5T_NATIVE_INT, &nStatesGlobal);
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);

            dataSpace = H5Screate_simple(1, dims, NULL);
            attr = H5Acreate2(restartFile, "nStatesTurbine", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            status = H5Awrite(attr, H5T_NATIVE_INT, &nStatesTurbine);
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);

        }

        if (nStatesGlobal > 0) {
            hsize_t dims[1];
            dims[0] = nStatesGlobal;
            hid_t dataSpace = H5Screate_simple(1, dims, NULL);
            hid_t dataSet = H5Dcreate2(restartFile, "/globStates", H5T_NATIVE_FLOAT, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            herr_t status = H5Dwrite(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, globStates.data());

            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);

            dataSpace = H5Screate_simple(1, dims, NULL);
            dataSet = H5Dcreate2(restartFile, "/globStates_np1", H5T_NATIVE_FLOAT, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, globStates_np1.data());

            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);


        }

        if (nStatesTurbine > 0) {

            hsize_t dims[2];
            dims[0] = nTurbinesGlob;
            dims[1] = nStatesTurbine;

            hid_t dataSpace = H5Screate_simple(2, dims, NULL);
            hid_t dataSet = H5Dcreate2(restartFile, "turbineStates", H5T_NATIVE_FLOAT, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            herr_t status = H5Dwrite(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, turbineStates.data());

            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);

            dataSpace = H5Screate_simple(2, dims, NULL);
            dataSet = H5Dcreate2(restartFile, "turbineStates_np1", H5T_NATIVE_FLOAT, dataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataSet, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, turbineStates_np1.data());

            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
        }

        herr_t status = H5Fclose(restartFile);
    }

  return 0;

}

void SuperController::fastSCInputOutput() {

    // Transfers
    // to_SC_np1 <------ ip_from_FAST
    // op_to_FAST <------- from_SC_np1, from_SCglob_np1

    for(int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
        for(int iInput=0; iInput < nCtrl2SC; iInput++) {
            to_SC_np1[iTurb*nCtrl2SC + iInput] = 0.0; // Initialize to zero
        }
    }

    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        for(int iInput=0; iInput < nCtrl2SC; iInput++) {
            to_SC_np1[turbineMapProcToGlob[iTurb]*nCtrl2SC + iInput] = ip_from_FAST[iTurb].toSC[iInput] ;
        }
    }

    if (MPI_COMM_NULL != fastMPIComm) {
        MPI_Allreduce(MPI_IN_PLACE, to_SC_np1.data(), nCtrl2SC*nTurbinesGlob, MPI_FLOAT, MPI_SUM, fastMPIComm) ;
    }

    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        for(int iOutput=0; iOutput < nSC2Ctrl; iOutput++) {
            op_to_FAST[iTurb].fromSC[iOutput] = from_SC_np1[turbineMapProcToGlob[iTurb]*nSC2Ctrl + iOutput] ;
        }

        for(int iOutput=0; iOutput < nSC2CtrlGlob; iOutput++) {
            op_to_FAST[iTurb].fromSCglob[iOutput] = from_SCglob_np1[turbineMapProcToGlob[iTurb]*nSC2Ctrl + iOutput] ;
        }

    }

}


void SuperController::advanceTime() {

    if (nTurbinesProc > 0) {

        for(int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
            for(int iInput=0; iInput < nCtrl2SC; iInput++) {
                to_SC_nm1[iTurb*nCtrl2SC + iInput] = to_SC_n[iTurb*nCtrl2SC + iInput];
                to_SC_n[iTurb*nCtrl2SC + iInput] = to_SC_np1[iTurb*nCtrl2SC + iInput];
//            to_SC_np1[iTurb*nCtrl2SC + iInput] = Predictor?
            }
            for(int iOutput=0; iOutput < nSC2Ctrl; iOutput++) {
                from_SC_nm1[iTurb*nSC2Ctrl + iOutput] = from_SC_n[iTurb*nSC2Ctrl + iOutput];
                from_SC_n[iTurb*nSC2Ctrl + iOutput] = from_SC_np1[iTurb*nSC2Ctrl + iOutput];
            }
        }

        for(int iInput=0; iInput < nInpGlobal; iInput++) {
            to_SCglob_nm1[iInput] = to_SCglob_n[iInput];
            to_SCglob_n[iInput] = to_SCglob_np1[iInput];
            //to_SCglob_np1[iInput] = Predictor?
        }

        for(int iOutput=0; iOutput < nSC2CtrlGlob; iOutput++) {
            from_SCglob_nm1[iOutput] = from_SCglob_n[iOutput];
            from_SCglob_n[iOutput] = from_SCglob_np1[iOutput];
            //from_SCglob_np1[iOutput] = Predictor?
        }

        for(int iState=0; iState<nStatesGlobal; iState++) {
            globStates[iState] = globStates_np1[iState];
        }

        for(int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
            for(int iState=0; iState < nStatesTurbine; iState++) {
                turbineStates[iTurb*nStatesTurbine + iState] = turbineStates_np1[iTurb*nStatesTurbine + iState];
            }
        }
    }

}
