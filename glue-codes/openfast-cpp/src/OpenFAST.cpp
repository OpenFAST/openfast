#include "OpenFAST.H"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cassert>

int fast::OpenFAST::AbortErrLev = ErrID_Fatal; // abort error level; compare with NWTC Library

//Constructor
fast::fastInputs::fastInputs():
nTurbinesGlob(0),
dryRun(false),
debug(false),
tStart(-1.0),
nEveryCheckPoint(-1),
tMax(0.0),
dtFAST(0.0),
scStatus(false),
scLibFile("")
{
  //Nothing to do here
}






//Constructor
fast::OpenFAST::OpenFAST():
nTurbinesGlob(0),
nTurbinesProc(0),
scStatus(false),
simStart(fast::init),
timeZero(false)
{
}

fast::OpenFAST::~OpenFAST(){ }

inline bool fast::OpenFAST::checkFileExists(const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

void fast::OpenFAST::init() {
    // Temporary buffer to pass filenames to OpenFAST fortran subroutines
    char currentFileName[INTERFACE_STRING_LENGTH];

    allocateMemory();

    if (!dryRun) {
        switch (simStart) {

        case fast::trueRestart:

            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                /* note that this will set nt_global inside the FAST library */
                std::copy(
                    CheckpointFileRoot[iTurb].data(),
                    CheckpointFileRoot[iTurb].data() + (CheckpointFileRoot[iTurb].size() + 1),
                    currentFileName
                );
                FAST_OpFM_Restart(
                    &iTurb,
                    currentFileName,
                    &AbortErrLev,
                    &dtFAST,
                    &numBlades[iTurb],
                    &numVelPtsBlade[iTurb],
                    &ntStart,
                    &cDriver_Input_from_FAST[iTurb],
                    &cDriver_Output_to_FAST[iTurb],
                    &sc.ip_from_FAST[iTurb],
                    &sc.op_to_FAST[iTurb],
                    &ErrStat,
                    ErrMsg
                );
                checkError(ErrStat, ErrMsg);
                nt_global = ntStart;

                int nfpts = get_numForcePtsLoc(iTurb);
                forceNodeVel[iTurb].resize(nfpts);
                for (int k = 0; k < nfpts; k++) forceNodeVel[iTurb][k].resize(3) ;
            }

            if (nTurbinesProc > 0) velNodeDataFile = openVelocityDataFile(false);

            if(scStatus) {
                std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
                //sc.readRestartFile(nt_global);
            }

            break ;

        case fast::init:

            sc.init(scio, nTurbinesProc);
            if(scStatus) {
                std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
                // sc.init_sc(scio, nTurbinesProc, turbineMapProcToGlob, fastMPIComm);
                // sc.calcOutputs_n(0.0);
            }

            // this calls the Init() routines of each module

            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                int nodeClusterType = 0;
                if (forcePtsBladeDistributionType[iTurb] == "chordClustered")
                {
                    nodeClusterType = 1;
                }
                std::copy(
                    FASTInputFileName[iTurb].data(),
                    FASTInputFileName[iTurb].data() + (FASTInputFileName[iTurb].size() + 1),
                    currentFileName
                );
                FAST_OpFM_Init(
                    &iTurb,
                    &tMax,
                    currentFileName,
                    &TurbID[iTurb],
                    &scio.nSC2CtrlGlob, 
                    &scio.nSC2Ctrl, 
                    &scio.nCtrl2SC, 
                    scio.from_SCglob.data(), 
                    scio.from_SC[iTurb].data(),
                    &numForcePtsBlade[iTurb],
                    &numForcePtsTwr[iTurb],
                    TurbineBasePos[iTurb].data(),
                    &AbortErrLev,
                    &dtFAST,
                    &numBlades[iTurb],
                    &numVelPtsBlade[iTurb],
                    &nodeClusterType,
                    &cDriver_Input_from_FAST[iTurb],
                    &cDriver_Output_to_FAST[iTurb],
                    &sc.ip_from_FAST[iTurb],
                    &sc.op_to_FAST[iTurb],
                    &ErrStat,
                    ErrMsg
                );
                checkError(ErrStat, ErrMsg);

                timeZero = true;

                numVelPtsTwr[iTurb] = cDriver_Output_to_FAST[iTurb].u_Len - numBlades[iTurb]*numVelPtsBlade[iTurb] - 1;
                if(numVelPtsTwr[iTurb] == 0) {
                    numForcePtsTwr[iTurb] = 0;
                    std::cout << "Aerodyn doesn't want to calculate forces on the tower. All actuator points on the tower are turned off for turbine " << turbineMapProcToGlob[iTurb] << "." << std::endl ;
                }

                int nfpts = get_numForcePtsLoc(iTurb);
                forceNodeVel[iTurb].resize(nfpts);
                for (int k = 0; k < nfpts; k++) forceNodeVel[iTurb][k].resize(3) ;

                if ( isDebug() ) {
                    for (int iNode=0; iNode < get_numVelPtsLoc(iTurb); iNode++) {
                        std::cout << "Node " << iNode << " Position = " << cDriver_Input_from_FAST[iTurb].pxVel[iNode] << " " << cDriver_Input_from_FAST[iTurb].pyVel[iNode] << " " << cDriver_Input_from_FAST[iTurb].pzVel[iNode] << " " << std::endl ;
                    }
                }
            }

            if (nTurbinesProc > 0) velNodeDataFile = openVelocityDataFile(true);

            break ;

        case fast::restartDriverInitFAST:

            sc.init(scio, nTurbinesProc);
            if(scStatus) {
                std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
                // sc.init_sc(scio, nTurbinesProc, turbineMapProcToGlob, fastMPIComm);
                // sc.calcOutputs_n(0.0);
            }
            
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                int nodeClusterType = 0;
                if (forcePtsBladeDistributionType[iTurb] == "chordClustered")
                {
                    nodeClusterType = 1;
                }
                std::copy(
                    FASTInputFileName[iTurb].data(),
                    FASTInputFileName[iTurb].data() + (FASTInputFileName[iTurb].size() + 1),
                    currentFileName
                );
                FAST_OpFM_Init(
                    &iTurb,
                    &tMax,
                    currentFileName,
                    &TurbID[iTurb],
                    &scio.nSC2CtrlGlob, 
                    &scio.nSC2Ctrl, 
                    &scio.nCtrl2SC,
                    scio.from_SCglob.data(), 
                    scio.from_SC[iTurb].data(),
                    &numForcePtsBlade[iTurb],
                    &numForcePtsTwr[iTurb],
                    TurbineBasePos[iTurb].data(),
                    &AbortErrLev,
                    &dtFAST,
                    &numBlades[iTurb],
                    &numVelPtsBlade[iTurb],
                    &nodeClusterType,
                    &cDriver_Input_from_FAST[iTurb],
                    &cDriver_Output_to_FAST[iTurb],
                    &sc.ip_from_FAST[iTurb],
                    &sc.op_to_FAST[iTurb],
                    &ErrStat,
                    ErrMsg
                );
                checkError(ErrStat, ErrMsg);

                timeZero = true;

                numVelPtsTwr[iTurb] = cDriver_Output_to_FAST[iTurb].u_Len - numBlades[iTurb]*numVelPtsBlade[iTurb] - 1;
                if(numVelPtsTwr[iTurb] == 0) {
                    numForcePtsTwr[iTurb] = 0;
                    std::cout << "Aerodyn doesn't want to calculate forces on the tower. All actuator points on the tower are turned off for turbine " << turbineMapProcToGlob[iTurb] << "." << std::endl ;
                }

                int nfpts = get_numForcePtsLoc(iTurb);
                forceNodeVel[iTurb].resize(nfpts);
                for (int k = 0; k < nfpts; k++) forceNodeVel[iTurb][k].resize(3) ;

                if ( isDebug() ) {
                    for (int iNode=0; iNode < get_numVelPtsLoc(iTurb); iNode++) {
                        std::cout << "Node " << iNode << " Position = " << cDriver_Input_from_FAST[iTurb].pxVel[iNode] << " " << cDriver_Input_from_FAST[iTurb].pyVel[iNode] << " " << cDriver_Input_from_FAST[iTurb].pzVel[iNode] << " " << std::endl ;
                    }
                }
            }

            int nTimesteps;

            if (nTurbinesProc > 0) {
                readVelocityData(ntStart);
            }
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                applyVelocityData(0, iTurb, cDriver_Output_to_FAST[iTurb], velNodeData[iTurb]);
            }
            solution0() ;

            for (int iPrestart=0 ; iPrestart < ntStart; iPrestart++) {
                for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                    applyVelocityData(iPrestart, iTurb, cDriver_Output_to_FAST[iTurb], velNodeData[iTurb]);
                }
                stepNoWrite();
            }

            if (nTurbinesProc > 0) velNodeDataFile = openVelocityDataFile(false);

            break;

        case fast::simStartType_END:

            break;

        }
    }
}

void fast::OpenFAST::solution0() {

    if (!dryRun) {
        // set wind speeds at initial locations
        // for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        //     setOutputsToFAST(cDriver_Input_from_FAST[iTurb], cDriver_Output_to_FAST[iTurb]);
        // }

        if(scStatus) {
            std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
            // sc.fastSCInputOutput();
        }

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_OpFM_Solution0(&iTurb, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
        }

        timeZero = false;

        if (scStatus) {
            std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;            
            //  sc.calcOutputs_n(0.0);
            //  sc.fastSCInputOutput();
        }
    }
}

void fast::OpenFAST::step() {

    /* ******************************
    set inputs from this code and call FAST:
    ********************************* */

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        //  set wind speeds at original locations
        //     setOutputsToFAST(cDriver_Input_from_FAST[iTurb], cDriver_Output_to_FAST[iTurb]);

        // this advances the states, calls CalcOutput, and solves for next inputs. Predictor-corrector loop is imbeded here:
        // (note OpenFOAM could do subcycling around this step)

        writeVelocityData(velNodeDataFile, iTurb, nt_global, cDriver_Input_from_FAST[iTurb], cDriver_Output_to_FAST[iTurb]);

        if ( isDebug() ) {

            std::ofstream fastcpp_velocity_file;
            fastcpp_velocity_file.open("fastcpp_velocity.csv") ;
            fastcpp_velocity_file << "# x, y, z, Vx, Vy, Vz" << std::endl ;
            for (int iNode=0; iNode < get_numVelPtsLoc(iTurb); iNode++) {
                fastcpp_velocity_file << cDriver_Input_from_FAST[iTurb].pxVel[iNode] << ", " << cDriver_Input_from_FAST[iTurb].pyVel[iNode] << ", " << cDriver_Input_from_FAST[iTurb].pzVel[iNode] << ", " << cDriver_Output_to_FAST[iTurb].u[iNode] << ", " << cDriver_Output_to_FAST[iTurb].v[iNode] << ", " << cDriver_Output_to_FAST[iTurb].w[iNode] << " " << std::endl ;
            }
            fastcpp_velocity_file.close() ;
        }

        FAST_OpFM_Step(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);

        // Compute the force from the nacelle only if the drag coefficient is
        //   greater than zero
        if (nacelle_cd[iTurb]>0.) {
            calc_nacelle_force (
                cDriver_Output_to_FAST[iTurb].u[0],
                cDriver_Output_to_FAST[iTurb].v[0],
                cDriver_Output_to_FAST[iTurb].w[0],
                nacelle_cd[iTurb],
                nacelle_area[iTurb],
                air_density[iTurb],
                cDriver_Input_from_FAST[iTurb].fx[0],
                cDriver_Input_from_FAST[iTurb].fy[0],
                cDriver_Input_from_FAST[iTurb].fz[0]
            );
        }

        if ( isDebug() ) {
            std::ofstream actuatorForcesFile;
            actuatorForcesFile.open("actuator_forces.csv") ;
            actuatorForcesFile << "# x, y, z, fx, fy, fz" << std::endl ;
            for (int iNode=0; iNode < get_numForcePtsLoc(iTurb); iNode++) {
                actuatorForcesFile << cDriver_Input_from_FAST[iTurb].pxForce[iNode] << ", " << cDriver_Input_from_FAST[iTurb].pyForce[iNode] << ", " << cDriver_Input_from_FAST[iTurb].pzForce[iNode] << ", " << cDriver_Input_from_FAST[iTurb].fx[iNode] << ", " << cDriver_Input_from_FAST[iTurb].fy[iNode] << ", " << cDriver_Input_from_FAST[iTurb].fz[iNode] << " " << std::endl ;
            }
            actuatorForcesFile.close() ;
        }
    }

    if(scStatus) {
        std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
        // sc.updateStates(nt_global * dtFAST); // Predict state at 'n+1' based on inputs
        // sc.calcOutputs_np1( (nt_global + 1) * dtFAST);
        // sc.fastSCInputOutput();
    }

    nt_global = nt_global + 1;
    
    if(scStatus) {
        std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
        // sc.advanceTime(); // Advance states, inputs and outputs from 'n' to 'n+1'
    }

    if ( (((nt_global - ntStart) % nEveryCheckPoint) == 0 )  && (nt_global != ntStart) ) {
        // Use default FAST naming convention for checkpoint file
        // <RootName>.<nt_global>
        char dummyCheckPointRoot[INTERFACE_STRING_LENGTH] = " ";
        // Ensure that we have a null character
        dummyCheckPointRoot[1] = 0;

        if (nTurbinesProc > 0) backupVelocityDataFile(nt_global, velNodeDataFile);

        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_CreateCheckpoint(&iTurb, dummyCheckPointRoot, &ErrStat, ErrMsg);
            checkError(ErrStat, ErrMsg);
        }
        if(scStatus) {
            std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
            // if (fastMPIRank == 0) {
            // sc.writeRestartFile(nt_global);
            // }
        }
    }
}

void fast::OpenFAST::stepNoWrite() {

    /* ******************************
    set inputs from this code and call FAST:
    ********************************* */

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        //  set wind speeds at original locations
        //     setOutputsToFAST(cDriver_Input_from_FAST[iTurb], cDriver_Output_to_FAST[iTurb]);

        // this advances the states, calls CalcOutput, and solves for next inputs. Predictor-corrector loop is imbeded here:
        // (note OpenFOAM could do subcycling around this step)
        FAST_OpFM_Step(&iTurb, &ErrStat, ErrMsg);
        checkError(ErrStat, ErrMsg);

    }

    if(scStatus) {
        std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
        // sc.updateStates( nt_global * dtFAST); // Predict state at 'n+1' based on inputs
        // sc.calcOutputs_np1( (nt_global+1) * dtFAST);
        // sc.fastSCInputOutput();
    }

    nt_global = nt_global + 1;

    if(scStatus) {
        std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
        // sc.advanceTime(); // Advance states, inputs and outputs from 'n' to 'n+1'
    }
}

void fast::OpenFAST::calc_nacelle_force(const float & u, const float & v, const float & w, const float & cd, const float & area, const float & rho, float & fx, float & fy, float & fz) {
    // Calculate the force on the nacelle (fx,fy,fz) given the
    //   velocity sampled at the nacelle point (u,v,w),
    //   drag coefficient 'cd' and nacelle area 'area'

    // The velocity magnitude
    float Vmag = std::sqrt(u * u + v * v + w * w);

    // Velocity correction based on Martinez-Tossas PhD Thesis 2017
    // The correction samples the velocity at the center of the
    // Gaussian kernel and scales it to obtain the inflow velocity
    float epsilon_d = std::sqrt(2.0 / M_PI * cd * area);
    float correction = 1. / (1.0 - cd * area / (4.0 * M_PI * epsilon_d * epsilon_d));

    // Compute the force for each velocity component
    fx = rho * 1./2. * cd * area * Vmag * u * correction * correction;
    fy = rho * 1./2. * cd * area * Vmag * v * correction * correction;
    fz = rho * 1./2. * cd * area * Vmag * w * correction * correction;
}

void fast::OpenFAST::setInputs(const fast::fastInputs & fi ) {

    mpiComm = fi.comm;

    MPI_Comm_rank(mpiComm, &worldMPIRank);
    MPI_Comm_group(mpiComm, &worldMPIGroup);

    nTurbinesGlob = fi.nTurbinesGlob;

    if (nTurbinesGlob > 0) {

        dryRun = fi.dryRun;
        debug = fi.debug;

        tStart = fi.tStart;
        simStart = fi.simStart;
        nEveryCheckPoint = fi.nEveryCheckPoint;
        tMax = fi.tMax;
        loadSuperController(fi);
        dtFAST = fi.dtFAST;

        ntStart = int(tStart/dtFAST);

        if (simStart == fast::restartDriverInitFAST) {
            nt_global = 0;
        } else {
            nt_global = ntStart;
        }

        globTurbineData.resize(nTurbinesGlob);
        globTurbineData = fi.globTurbineData;

    } else {
        throw std::runtime_error("Number of turbines < 0 ");
    }
}

void fast::OpenFAST::checkError(const int ErrStat, const char * ErrMsg){
    if (ErrStat != ErrID_None){
        if (ErrStat >= AbortErrLev){
            throw std::runtime_error(ErrMsg);
        }
    }
}

void fast::OpenFAST::setOutputsToFAST(OpFM_InputType_t cDriver_Input_from_FAST, OpFM_OutputType_t cDriver_Output_to_FAST){

    // routine sets the u-v-w wind speeds used in FAST and the SuperController inputs

    for (int j = 0; j < cDriver_Output_to_FAST.u_Len; j++){
        cDriver_Output_to_FAST.u[j] = (float) 10.0*pow((cDriver_Input_from_FAST.pzVel[j] / 90.0), 0.2); // 0.2 power law wind profile using reference 10 m/s at 90 meters
        cDriver_Output_to_FAST.v[j] = 0.0;
        cDriver_Output_to_FAST.w[j] = 0.0;
    }

}

void fast::OpenFAST::getApproxHubPos(double* currentCoords, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Get hub position of Turbine 'iTurbGlob'
    for(int i =0; i<nSize; ++i){
        currentCoords[i] = globTurbineData[iTurbGlob].TurbineHubPos[i];
    }
}

void fast::OpenFAST::getHubPos(double* currentCoords, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Get hub position of Turbine 'iTurbGlob'
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    currentCoords[0] = cDriver_Input_from_FAST[iTurbLoc].pxVel[0] + TurbineBasePos[iTurbLoc][0] ;
    currentCoords[1] = cDriver_Input_from_FAST[iTurbLoc].pyVel[0] + TurbineBasePos[iTurbLoc][1] ;
    currentCoords[2] = cDriver_Input_from_FAST[iTurbLoc].pzVel[0] + TurbineBasePos[iTurbLoc][2] ;
}

void fast::OpenFAST::getHubShftDir(double* hubShftVec, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Get hub shaft direction of current turbine - pointing downwind
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int i=0; i<nSize; i++){
        hubShftVec[i] = cDriver_Input_from_FAST[iTurbLoc].pOrientation[i*3] ;
    }
}

void fast::OpenFAST::getVelNodeCoordinates(double* currentCoords, int iNode, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Set coordinates at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbLoc);
    currentCoords[0] = cDriver_Input_from_FAST[iTurbLoc].pxVel[iNode] + TurbineBasePos[iTurbLoc][0] ;
    currentCoords[1] = cDriver_Input_from_FAST[iTurbLoc].pyVel[iNode] + TurbineBasePos[iTurbLoc][1] ;
    currentCoords[2] = cDriver_Input_from_FAST[iTurbLoc].pzVel[iNode] + TurbineBasePos[iTurbLoc][2] ;
}

void fast::OpenFAST::getForceNodeCoordinates(double* currentCoords, int iNode, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Set coordinates at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    currentCoords[0] = cDriver_Input_from_FAST[iTurbLoc].pxForce[iNode] + TurbineBasePos[iTurbLoc][0] ;
    currentCoords[1] = cDriver_Input_from_FAST[iTurbLoc].pyForce[iNode] + TurbineBasePos[iTurbLoc][1] ;
    currentCoords[2] = cDriver_Input_from_FAST[iTurbLoc].pzForce[iNode] + TurbineBasePos[iTurbLoc][2] ;
}

void fast::OpenFAST::getForceNodeOrientation(double* currentOrientation, int iNode, int iTurbGlob, int nSize) {
    assert(nSize==9);
    // Set orientation at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    for(int i=0;i<9;i++) {
        currentOrientation[i] = cDriver_Input_from_FAST[iTurbLoc].pOrientation[iNode*9+i] ;
    }
}

void fast::OpenFAST::getRelativeVelForceNode(double* currentVelocity, int iNode, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Get relative velocity at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);

    currentVelocity[0] = forceNodeVel[iTurbLoc][iNode][0] - cDriver_Input_from_FAST[iTurbLoc].xdotForce[iNode];
    currentVelocity[1] = forceNodeVel[iTurbLoc][iNode][1] - cDriver_Input_from_FAST[iTurbLoc].ydotForce[iNode];
    currentVelocity[2] = forceNodeVel[iTurbLoc][iNode][2] - cDriver_Input_from_FAST[iTurbLoc].zdotForce[iNode];
}

void fast::OpenFAST::getForce(double* currentForce, int iNode, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Set forces at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    currentForce[0] = -cDriver_Input_from_FAST[iTurbLoc].fx[iNode] ;
    currentForce[1] = -cDriver_Input_from_FAST[iTurbLoc].fy[iNode] ;
    currentForce[2] = -cDriver_Input_from_FAST[iTurbLoc].fz[iNode] ;
}

double fast::OpenFAST::getChord(int iNode, int iTurbGlob) {
    // Return blade chord/tower diameter at current node of current turbine
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
    return cDriver_Input_from_FAST[iTurbLoc].forceNodesChord[iNode] ;
}

void fast::OpenFAST::setVelocity(double* currentVelocity, int iNode, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Set velocity at current node of current turbine -
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbLoc);
    cDriver_Output_to_FAST[iTurbLoc].u[iNode] = currentVelocity[0];
    cDriver_Output_to_FAST[iTurbLoc].v[iNode] = currentVelocity[1];
    cDriver_Output_to_FAST[iTurbLoc].w[iNode] = currentVelocity[2];
}

void fast::OpenFAST::setVelocityForceNode(double* currentVelocity, int iNode, int iTurbGlob, int nSize) {
    assert(nSize==3);
    // Set velocity at current node of current turbine -
    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);

    for(int i=0; i<nSize; ++i){
        forceNodeVel[iTurbLoc][iNode][i] = currentVelocity[i];
    }
}

void fast::OpenFAST::interpolateVel_ForceToVelNodes() {

    // Interpolates the velocity from the force nodes to the velocity nodes
    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        // Hub location
        cDriver_Output_to_FAST[iTurb].u[0] = forceNodeVel[iTurb][0][0];
        cDriver_Output_to_FAST[iTurb].v[0] = forceNodeVel[iTurb][0][1];
        cDriver_Output_to_FAST[iTurb].w[0] = forceNodeVel[iTurb][0][2];

        if ( isDebug() ) {
            std::ofstream actuatorVelFile;
            actuatorVelFile.open("actuator_velocity.csv") ;
            actuatorVelFile << "# x, y, z, Vx, Vy, Vz" << std::endl ;
            for (int iNode=0; iNode < get_numForcePtsLoc(iTurb); iNode++) {
                actuatorVelFile << cDriver_Input_from_FAST[iTurb].pxForce[iNode] << ", " << cDriver_Input_from_FAST[iTurb].pyForce[iNode] << ", " << cDriver_Input_from_FAST[iTurb].pzForce[iNode] << ", " << forceNodeVel[iTurb][iNode][0] << ", " << forceNodeVel[iTurb][iNode][1] << ", " << forceNodeVel[iTurb][iNode][2] << " " << std::endl ;
            }
            actuatorVelFile.close() ;
        }

        // Do the blades first
        int nBlades = get_numBladesLoc(iTurb);
        for(int iBlade=0; iBlade < nBlades; iBlade++) {
            // Create interpolating parameter - Distance from hub
            int nForcePtsBlade = get_numForcePtsBladeLoc(iTurb);
            std::vector<double> rDistForce(nForcePtsBlade) ;
            for(int j=0; j < nForcePtsBlade; j++) {
                int iNodeForce = 1 + iBlade * nForcePtsBlade + j ; //The number of actuator force points is always the same for all blades
                rDistForce[j] = std::sqrt(
                    (cDriver_Input_from_FAST[iTurb].pxForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pxForce[0])*(cDriver_Input_from_FAST[iTurb].pxForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pxForce[0])
                    + (cDriver_Input_from_FAST[iTurb].pyForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pyForce[0])*(cDriver_Input_from_FAST[iTurb].pyForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pyForce[0])
                    + (cDriver_Input_from_FAST[iTurb].pzForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pzForce[0])*(cDriver_Input_from_FAST[iTurb].pzForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pzForce[0])
                );
            }

            // Interpolate to the velocity nodes
            int nVelPtsBlade = get_numVelPtsBladeLoc(iTurb);
            for(int j=0; j < nVelPtsBlade; j++) {
                int iNodeVel = 1 + iBlade * nVelPtsBlade + j ; //Assumes the same number of velocity (Aerodyn) nodes for all blades
                double rDistVel = std::sqrt(
                    (cDriver_Input_from_FAST[iTurb].pxVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pxVel[0])*(cDriver_Input_from_FAST[iTurb].pxVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pxVel[0])
                    + (cDriver_Input_from_FAST[iTurb].pyVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pyVel[0])*(cDriver_Input_from_FAST[iTurb].pyVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pyVel[0])
                    + (cDriver_Input_from_FAST[iTurb].pzVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pzVel[0])*(cDriver_Input_from_FAST[iTurb].pzVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pzVel[0])
                );
                //Find nearest two force nodes
                int jForceLower = 0;
                while ( (rDistForce[jForceLower+1] < rDistVel) && ( jForceLower < (nForcePtsBlade-2)) ) {
                    jForceLower = jForceLower + 1;
                }
                int iNodeForceLower = 1 + iBlade * nForcePtsBlade + jForceLower ;
                double rInterp = (rDistVel - rDistForce[jForceLower])/(rDistForce[jForceLower+1]-rDistForce[jForceLower]);
                cDriver_Output_to_FAST[iTurb].u[iNodeVel] = forceNodeVel[iTurb][iNodeForceLower][0] + rInterp * (forceNodeVel[iTurb][iNodeForceLower+1][0] - forceNodeVel[iTurb][iNodeForceLower][0] );
                cDriver_Output_to_FAST[iTurb].v[iNodeVel] = forceNodeVel[iTurb][iNodeForceLower][1] + rInterp * (forceNodeVel[iTurb][iNodeForceLower+1][1] - forceNodeVel[iTurb][iNodeForceLower][1] );
                cDriver_Output_to_FAST[iTurb].w[iNodeVel] = forceNodeVel[iTurb][iNodeForceLower][2] + rInterp * (forceNodeVel[iTurb][iNodeForceLower+1][2] - forceNodeVel[iTurb][iNodeForceLower][2] );
            }
        }

        // Now the tower if present and used
        int nVelPtsTower = get_numVelPtsTwrLoc(iTurb);
        if ( nVelPtsTower > 0 ) {

            // Create interpolating parameter - Distance from first node from ground
            int nForcePtsTower = get_numForcePtsTwrLoc(iTurb);
            std::vector<double> hDistForce(nForcePtsTower) ;
            int iNodeBotTowerForce = 1 + nBlades * get_numForcePtsBladeLoc(iTurb); // The number of actuator force points is always the same for all blades
            for(int j=0; j < nForcePtsTower; j++) {
                int iNodeForce = iNodeBotTowerForce + j ;
                hDistForce[j] = std::sqrt(
                    (cDriver_Input_from_FAST[iTurb].pxForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pxForce[iNodeBotTowerForce])*(cDriver_Input_from_FAST[iTurb].pxForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pxForce[iNodeBotTowerForce])
                    + (cDriver_Input_from_FAST[iTurb].pyForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pyForce[iNodeBotTowerForce])*(cDriver_Input_from_FAST[iTurb].pyForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pyForce[iNodeBotTowerForce])
                    + (cDriver_Input_from_FAST[iTurb].pzForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pzForce[iNodeBotTowerForce])*(cDriver_Input_from_FAST[iTurb].pzForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pzForce[iNodeBotTowerForce])
                );
            }

            int iNodeBotTowerVel = 1 + nBlades * get_numVelPtsBladeLoc(iTurb); // Assumes the same number of velocity (Aerodyn) nodes for all blades
            for(int j=0; j < nVelPtsTower; j++) {
                int iNodeVel = iNodeBotTowerVel + j ;
                double hDistVel = std::sqrt(
                    (cDriver_Input_from_FAST[iTurb].pxVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pxVel[iNodeBotTowerVel])*(cDriver_Input_from_FAST[iTurb].pxVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pxVel[iNodeBotTowerVel])
                    + (cDriver_Input_from_FAST[iTurb].pyVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pyVel[iNodeBotTowerVel])*(cDriver_Input_from_FAST[iTurb].pyVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pyVel[iNodeBotTowerVel])
                    + (cDriver_Input_from_FAST[iTurb].pzVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pzVel[iNodeBotTowerVel])*(cDriver_Input_from_FAST[iTurb].pzVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pzVel[iNodeBotTowerVel])
                );
                //Find nearest two force nodes
                int jForceLower = 0;
                while ( (hDistForce[jForceLower+1] < hDistVel) && ( jForceLower < (nForcePtsTower-2)) )   {
                    jForceLower = jForceLower + 1;
                }
                int iNodeForceLower = iNodeBotTowerForce + jForceLower ;
                double rInterp = (hDistVel - hDistForce[jForceLower])/(hDistForce[jForceLower+1]-hDistForce[jForceLower]);
                cDriver_Output_to_FAST[iTurb].u[iNodeVel] = forceNodeVel[iTurb][iNodeForceLower][0] + rInterp * (forceNodeVel[iTurb][iNodeForceLower+1][0] - forceNodeVel[iTurb][iNodeForceLower][0] );
                cDriver_Output_to_FAST[iTurb].v[iNodeVel] = forceNodeVel[iTurb][iNodeForceLower][1] + rInterp * (forceNodeVel[iTurb][iNodeForceLower+1][1] - forceNodeVel[iTurb][iNodeForceLower][1] );
                cDriver_Output_to_FAST[iTurb].w[iNodeVel] = forceNodeVel[iTurb][iNodeForceLower][2] + rInterp * (forceNodeVel[iTurb][iNodeForceLower+1][2] - forceNodeVel[iTurb][iNodeForceLower][2] );
            }
        }
    }
}

void fast::OpenFAST::computeTorqueThrust(int iTurbGlob, std::vector<double> & torque, std::vector<double> & thrust) {

    //Compute the torque and thrust based on the forces at the actuator nodes
    std::vector<double> relLoc(3,0.0);
    std::vector<double> rPerpShft(3);
    thrust[0] = 0.0; thrust[1] = 0.0; thrust[2] = 0.0;
    torque[0] = 0.0; torque[1] = 0.0; torque[2] = 0.0;

    std::vector<double> hubShftVec(3);
    getHubShftDir(hubShftVec, iTurbGlob);

    int iTurbLoc = get_localTurbNo(iTurbGlob) ;
    for (int k=0; k < get_numBladesLoc(iTurbLoc); k++) {
        for (int j=0; j < numForcePtsBlade[iTurbLoc]; j++) {
            int iNode = 1 + numForcePtsBlade[iTurbLoc]*k + j ;

            thrust[0] = thrust[0] + cDriver_Input_from_FAST[iTurbLoc].fx[iNode] ;
            thrust[1] = thrust[1] + cDriver_Input_from_FAST[iTurbLoc].fy[iNode] ;
            thrust[2] = thrust[2] + cDriver_Input_from_FAST[iTurbLoc].fz[iNode] ;

            relLoc[0] = cDriver_Input_from_FAST[iTurbLoc].pxForce[iNode] - cDriver_Input_from_FAST[iTurbLoc].pxForce[0];
            relLoc[1] = cDriver_Input_from_FAST[iTurbLoc].pyForce[iNode] - cDriver_Input_from_FAST[iTurbLoc].pyForce[0];
            relLoc[2] = cDriver_Input_from_FAST[iTurbLoc].pzForce[iNode] - cDriver_Input_from_FAST[iTurbLoc].pzForce[0];

            double rDotHubShftVec = relLoc[0]*hubShftVec[0] + relLoc[1]*hubShftVec[1] + relLoc[2]*hubShftVec[2];
            for (int j=0; j < 3; j++)  rPerpShft[j] = relLoc[j] - rDotHubShftVec * hubShftVec[j];

            torque[0] = torque[0] + rPerpShft[1] * cDriver_Input_from_FAST[iTurbLoc].fz[iNode] - rPerpShft[2] * cDriver_Input_from_FAST[iTurbLoc].fy[iNode] + cDriver_Input_from_FAST[iTurbLoc].momentx[iNode] ;
            torque[1] = torque[1] + rPerpShft[2] * cDriver_Input_from_FAST[iTurbLoc].fx[iNode] - rPerpShft[0] * cDriver_Input_from_FAST[iTurbLoc].fz[iNode] + cDriver_Input_from_FAST[iTurbLoc].momenty[iNode] ;
            torque[2] = torque[2] + rPerpShft[0] * cDriver_Input_from_FAST[iTurbLoc].fy[iNode] - rPerpShft[1] * cDriver_Input_from_FAST[iTurbLoc].fx[iNode] + cDriver_Input_from_FAST[iTurbLoc].momentz[iNode] ;
        }
    }
}

fast::ActuatorNodeType fast::OpenFAST::getVelNodeType(int iTurbGlob, int iNode) {
    // Return the type of velocity node for the given node number. The node ordering (from FAST) is
    // Node 0 - Hub node
    // Blade 1 nodes
    // Blade 2 nodes
    // Blade 3 nodes
    // Tower nodes

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbGlob);
    if (iNode) {
        if ( (iNode + 1 - (get_numVelPts(iTurbLoc) - get_numVelPtsTwr(iTurbLoc)) ) > 0 ) {
            return TOWER;
        } else {
            return BLADE;
        }
    } else {
        return HUB;
    }
}

fast::ActuatorNodeType fast::OpenFAST::getForceNodeType(int iTurbGlob, int iNode) {
    // Return the type of actuator force node for the given node number. The node ordering (from FAST) is
    // Node 0 - Hub node
    // Blade 1 nodes
    // Blade 2 nodes
    // Blade 3 nodes
    // Tower nodes

    int iTurbLoc = get_localTurbNo(iTurbGlob);
    for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbGlob);
    if (iNode) {
        if ( (iNode + 1 - (get_numForcePts(iTurbLoc) - get_numForcePtsTwr(iTurbLoc)) ) > 0 ) {
            return TOWER;
        } else {
            return BLADE;
        }
    } else {
        return HUB;
    }
}

void fast::OpenFAST::allocateMemory() {

    for (int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
        if (dryRun) {
            if(worldMPIRank == 0) {
                std::cout << "iTurb = " << iTurb << " turbineMapGlobToProc[iTurb] = " << turbineMapGlobToProc[iTurb] << std::endl ;
            }
        }
        if(worldMPIRank == turbineMapGlobToProc[iTurb]) {
            turbineMapProcToGlob[nTurbinesProc] = iTurb;
            reverseTurbineMapProcToGlob[iTurb] = nTurbinesProc;
            nTurbinesProc++ ;
        }
        turbineSetProcs.insert(turbineMapGlobToProc[iTurb]);
    }

    int nProcsWithTurbines=0;
    turbineProcs.resize(turbineSetProcs.size());

    for (std::set<int>::const_iterator p = turbineSetProcs.begin(); p != turbineSetProcs.end(); p++) {
        turbineProcs[nProcsWithTurbines] = *p;
        nProcsWithTurbines++ ;
    }

    if (dryRun) {
        if (nTurbinesProc > 0) {
            std::ofstream turbineAllocFile;
            turbineAllocFile.open("turbineAlloc." + std::to_string(worldMPIRank) + ".txt") ;
            for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
                turbineAllocFile << "Proc " << worldMPIRank << " loc iTurb " << iTurb << " glob iTurb " << turbineMapProcToGlob[iTurb] << std::endl ;
            }
            turbineAllocFile.flush();
            turbineAllocFile.close() ;
        }
    }

    // Construct a group containing all procs running atleast 1 turbine in FAST
    MPI_Group_incl(worldMPIGroup, nProcsWithTurbines, &turbineProcs[0], &fastMPIGroup) ;
    int fastMPIcommTag = MPI_Comm_create(mpiComm, fastMPIGroup, &fastMPIComm);
    if (MPI_COMM_NULL != fastMPIComm) {
        MPI_Comm_rank(fastMPIComm, &fastMPIRank);
    }

    TurbID.resize(nTurbinesProc);
    TurbineBasePos.resize(nTurbinesProc);
    FASTInputFileName.resize(nTurbinesProc);
    CheckpointFileRoot.resize(nTurbinesProc);
    nacelle_cd.resize(nTurbinesProc);
    nacelle_area.resize(nTurbinesProc);
    air_density.resize(nTurbinesProc);
    numBlades.resize(nTurbinesProc);
    forcePtsBladeDistributionType.resize(nTurbinesProc);
    numForcePtsBlade.resize(nTurbinesProc);
    numForcePtsTwr.resize(nTurbinesProc);
    numVelPtsBlade.resize(nTurbinesProc);
    numVelPtsTwr.resize(nTurbinesProc);
    forceNodeVel.resize(nTurbinesProc);

    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

        TurbineBasePos[iTurb].resize(3);

        int globProc = turbineMapProcToGlob[iTurb];
        TurbID[iTurb] = globTurbineData[globProc].TurbID;
        FASTInputFileName[iTurb] = globTurbineData[globProc].FASTInputFileName ;
        CheckpointFileRoot[iTurb] = globTurbineData[globProc].FASTRestartFileName ;
        for(int i=0;i<3;i++) {
            TurbineBasePos[iTurb][i] = globTurbineData[globProc].TurbineBasePos[i];
        }
        forcePtsBladeDistributionType[iTurb] =  globTurbineData[globProc].forcePtsBladeDistributionType;
        numForcePtsBlade[iTurb] = globTurbineData[globProc].numForcePtsBlade;
        numForcePtsTwr[iTurb] = globTurbineData[globProc].numForcePtsTwr;
        nacelle_cd[iTurb] = globTurbineData[globProc].nacelle_cd;
        nacelle_area[iTurb] = globTurbineData[globProc].nacelle_area;
        air_density[iTurb] = globTurbineData[globProc].air_density;
    }

    // Allocate memory for Turbine datastructure for all turbines
    FAST_AllocateTurbines(&nTurbinesProc, &ErrStat, ErrMsg);

    // Allocate memory for OpFM Input types in FAST
    cDriver_Input_from_FAST.resize(nTurbinesProc) ;
    cDriver_Output_to_FAST.resize(nTurbinesProc) ;

    if(scStatus) {
        std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
        // scio.from_SC.resize(nTurbinesProc);
    }
}

void fast::OpenFAST::allocateTurbinesToProcsSimple() {
    // Allocate turbines to each processor - round robin fashion
    int nProcs ;
    MPI_Comm_size(mpiComm, &nProcs);
    for(int j = 0; j < nTurbinesGlob; j++)  turbineMapGlobToProc[j] = j % nProcs ;
}

void fast::OpenFAST::end() {
    // Deallocate types we allocated earlier

    if (nTurbinesProc > 0) closeVelocityDataFile(nt_global, velNodeDataFile);

    if ( !dryRun) {
        bool stopTheProgram = false;
        for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
            FAST_End(&iTurb, &stopTheProgram);
        }
        FAST_DeallocateTurbines(&ErrStat, ErrMsg);
    }

    MPI_Group_free(&fastMPIGroup);
    if (MPI_COMM_NULL != fastMPIComm) {
        MPI_Comm_free(&fastMPIComm);
    }
    MPI_Group_free(&worldMPIGroup);

    if(scStatus) {
        std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
        // sc.end();
    }
}

void fast::OpenFAST::readVelocityData(int nTimesteps) {

    int nTurbines;

    hid_t velDataFile = H5Fopen(("velDatafile." + std::to_string(worldMPIRank) + ".h5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    {
        hid_t attr = H5Aopen(velDataFile, "nTurbines", H5P_DEFAULT);
        herr_t ret = H5Aread(attr, H5T_NATIVE_INT, &nTurbines) ;
        H5Aclose(attr);
    }

    // Allocate memory and read the velocity data.
    velNodeData.resize(nTurbines);
    for (int iTurb=0; iTurb < nTurbines; iTurb++) {
        int nVelPts = get_numVelPtsLoc(iTurb) ;
        velNodeData[iTurb].resize(nTimesteps*nVelPts*6) ;
        hid_t dset_id = H5Dopen2(velDataFile, ("/turbine" + std::to_string(iTurb)).c_str(), H5P_DEFAULT);
        hid_t dspace_id = H5Dget_space(dset_id);

        hsize_t start[3]; start[1] = 0; start[2] = 0;
        hsize_t count[3]; count[0] = 1; count[1] = nVelPts; count[2] = 6;
        hid_t mspace_id = H5Screate_simple(3, count, NULL);

        for (int iStep=0; iStep < nTimesteps; iStep++) {
            start[0] = iStep;
            H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
            herr_t status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, mspace_id, dspace_id, H5P_DEFAULT, &velNodeData[iTurb][iStep*nVelPts*6] );
        }

        herr_t status = H5Dclose(dset_id);
    }
}

hid_t fast::OpenFAST::openVelocityDataFile(bool createFile) {

    hid_t velDataFile;
    if (createFile) {
        // Open the file in create mode
        velDataFile = H5Fcreate(("velDatafile." + std::to_string(worldMPIRank) + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        {
            hsize_t dims[1];
            dims[0] = 1;
            hid_t dataSpace = H5Screate_simple(1, dims, NULL);
            hid_t attr = H5Acreate2(velDataFile, "nTurbines", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            herr_t status = H5Awrite(attr, H5T_NATIVE_INT, &nTurbinesProc);
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);
            
            dataSpace = H5Screate_simple(1, dims, NULL);
            attr = H5Acreate2(velDataFile, "nTimesteps", H5T_NATIVE_INT, dataSpace, H5P_DEFAULT, H5P_DEFAULT) ;
            status = H5Aclose(attr);
            status = H5Sclose(dataSpace);
        }

        int ntMax = tMax/dtFAST ;

        for (int iTurb = 0; iTurb < nTurbinesProc; iTurb++) {
            int nVelPts = get_numVelPtsLoc(iTurb);
            hsize_t dims[3];
            dims[0] = ntMax; dims[1] = nVelPts; dims[2] = 6 ;

            hsize_t chunk_dims[3];
            chunk_dims[0] = 1; chunk_dims[1] = nVelPts; chunk_dims[2] = 6;
            hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(dcpl_id, 3, chunk_dims);

            hid_t dataSpace = H5Screate_simple(3, dims, NULL);
            hid_t dataSet = H5Dcreate(velDataFile, ("/turbine" + std::to_string(iTurb)).c_str(), H5T_NATIVE_DOUBLE, dataSpace, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

            herr_t status = H5Pclose(dcpl_id);
            status = H5Dclose(dataSet);
            status = H5Sclose(dataSpace);
        }

    } else {
        // Open the file in append mode
        velDataFile = H5Fopen(("velDatafile." + std::to_string(worldMPIRank) + ".h5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }

    return velDataFile;

}

herr_t fast::OpenFAST::closeVelocityDataFile(int nt_global, hid_t velDataFile) {
    herr_t status = H5Fclose(velDataFile) ;
    return status;
}

void fast::OpenFAST::backupVelocityDataFile(int curTimeStep, hid_t & velDataFile) {

    closeVelocityDataFile(curTimeStep, velDataFile);

    std::ifstream source("velDatafile." + std::to_string(worldMPIRank) + ".h5", std::ios::binary);
    std::ofstream dest("velDatafile." + std::to_string(worldMPIRank) + ".h5." + std::to_string(curTimeStep) + ".bak", std::ios::binary);

    dest << source.rdbuf();
    source.close();
    dest.close();

    velDataFile = openVelocityDataFile(false);
}

void fast::OpenFAST::writeVelocityData(hid_t h5File, int iTurb, int iTimestep, OpFM_InputType_t iData, OpFM_OutputType_t oData) {

    hsize_t start[3]; start[0] = iTimestep; start[1] = 0; start[2] = 0;
    int nVelPts = get_numVelPtsLoc(iTurb) ;
    hsize_t count[3]; count[0] = 1; count[1] = nVelPts; count[2] = 6;

    std::vector<double> tmpVelData;
    tmpVelData.resize(nVelPts * 6);

    for (int iNode=0 ; iNode < nVelPts; iNode++) {
        tmpVelData[iNode*6 + 0] = iData.pxVel[iNode];
        tmpVelData[iNode*6 + 1] = iData.pyVel[iNode];
        tmpVelData[iNode*6 + 2] = iData.pzVel[iNode];
        tmpVelData[iNode*6 + 3] = oData.u[iNode];
        tmpVelData[iNode*6 + 4] = oData.v[iNode];
        tmpVelData[iNode*6 + 5] = oData.w[iNode];
    }

    hid_t dset_id = H5Dopen2(h5File, ("/turbine" + std::to_string(iTurb)).c_str(), H5P_DEFAULT);
    hid_t dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, start, NULL, count, NULL);
    hid_t mspace_id = H5Screate_simple(3, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mspace_id, dspace_id, H5P_DEFAULT, tmpVelData.data());

    H5Dclose(dset_id);
    H5Sclose(dspace_id);
    H5Sclose(mspace_id);

    hid_t attr_id = H5Aopen_by_name(h5File, ".", "nTimesteps", H5P_DEFAULT, H5P_DEFAULT);
    herr_t status = H5Awrite(attr_id, H5T_NATIVE_INT, &iTimestep);
    status = H5Aclose(attr_id);

}

void fast::OpenFAST::applyVelocityData(int iPrestart, int iTurb, OpFM_OutputType_t cDriver_Output_to_FAST, std::vector<double> & velData) {
    int nVelPts = get_numVelPtsLoc(iTurb);
    for (int j = 0; j < nVelPts; j++){
        cDriver_Output_to_FAST.u[j] = velData[(iPrestart*nVelPts+j)*6 + 3];
        cDriver_Output_to_FAST.v[j] = velData[(iPrestart*nVelPts+j)*6 + 4];
        cDriver_Output_to_FAST.w[j] = velData[(iPrestart*nVelPts+j)*6 + 5];
    }
}

void fast::OpenFAST::loadSuperController(const fast::fastInputs & fi) {

    if(fi.scStatus) {
        std::cout << "Use of Supercontroller is not supported through the C++ API right now" << std::endl;
        // scStatus = fi.scStatus;
        // sc.load(fi.nTurbinesGlob, fi.scLibFile, scio);

    } else {
        scStatus = false;
    }
}
