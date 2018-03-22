#include "OpenFAST.H"

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
scLibFile(""),
numScInputsTurbine(0),
numScOutputsTurbine(0)
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

inline bool fast::OpenFAST::checkFileExists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

void fast::OpenFAST::init() {

  allocateMemory();

  if (!dryRun) {
    switch (simStart) {

    case fast::trueRestart:

     for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
       /* note that this will set nt_global inside the FAST library */
       FAST_OpFM_Restart(&iTurb, CheckpointFileRoot[iTurb].data(), &AbortErrLev, &dtFAST, &numBlades[iTurb], &numVelPtsBlade[iTurb], &ntStart, &cDriver_Input_from_FAST[iTurb], &cDriver_Output_to_FAST[iTurb], &scInputsTurbine_from_FAST[iTurb], &scOutputsTurbine_to_FAST[iTurb], &ErrStat, ErrMsg);
       checkError(ErrStat, ErrMsg);
       nt_global = ntStart;

       int nfpts = get_numForcePtsLoc(iTurb);
       forceNodeVel[iTurb].resize(nfpts);
       for (int k = 0; k < nfpts; k++) forceNodeVel[iTurb][k].resize(3) ;

     }

     if (nTurbinesProc > 0) velNodeDataFile = openVelocityDataFile(false);

     if(scStatus) {
	 sc.readRestartFile(nt_global);
     }
     
     break ;
   
    case fast::init:
     
     if(scStatus) {
        sc.init(nTurbinesGlob);
        sc.calcOutputs_n(0.0, scInputsGlob_n, scInputsTurbine_n, scOutputsGlob_n, scOutputsTurbine_n);
     }

     // this calls the Init() routines of each module
     for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
         FAST_OpFM_Init(&iTurb, &tMax, FASTInputFileName[iTurb].data(), &TurbID[iTurb], &numScOutputsGlob, &numScOutputsTurbine, &numScInputsTurbine, scOutputsGlob_n.data(), scOutputsTurbine_n.data(), &numForcePtsBlade[iTurb], &numForcePtsTwr[iTurb], TurbineBasePos[iTurb].data(), &AbortErrLev, &dtFAST, &numBlades[iTurb], &numVelPtsBlade[iTurb], &cDriver_Input_from_FAST[iTurb], &cDriver_Output_to_FAST[iTurb], &scInputsTurbine_from_FAST[iTurb], &scOutputsTurbine_to_FAST[iTurb], &ErrStat, ErrMsg);
       checkError(ErrStat, ErrMsg);
       
       timeZero = true;

       numVelPtsTwr[iTurb] = cDriver_Output_to_FAST[iTurb].u_Len - numBlades[iTurb]*numVelPtsBlade[iTurb] - 1;

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

     if(scStatus) {
         sc.init(nTurbinesGlob);
         sc.calcOutputs_n(0.0, scInputsGlob_n, scInputsTurbine_n, scOutputsGlob_n, scOutputsTurbine_n);
     }
     
     for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
         FAST_OpFM_Init(&iTurb, &tMax, FASTInputFileName[iTurb].data(), &TurbID[iTurb], &numScOutputsGlob, &numScOutputsTurbine, &numScInputsTurbine, scOutputsGlob_n.data(), scOutputsTurbine_n.data(), &numForcePtsBlade[iTurb], &numForcePtsTwr[iTurb], TurbineBasePos[iTurb].data(), &AbortErrLev, &dtFAST, &numBlades[iTurb], &numVelPtsBlade[iTurb], &cDriver_Input_from_FAST[iTurb], &cDriver_Output_to_FAST[iTurb], &scInputsTurbine_from_FAST[iTurb], &scOutputsTurbine_to_FAST[iTurb], &ErrStat, ErrMsg);
       checkError(ErrStat, ErrMsg);
       
       timeZero = true;

       numVelPtsTwr[iTurb] = cDriver_Output_to_FAST[iTurb].u_Len - numBlades[iTurb]*numVelPtsBlade[iTurb] - 1;

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
     
      
     if (scStatus) {
         fastScInputOutput(scInputsTurbine_from_FAST, scInputsTurbine_np1, scOutputsTurbine_np1, scOutputsTurbine_to_FAST);
     }
     
     for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {

       FAST_OpFM_Solution0(&iTurb, &ErrStat, ErrMsg);
       checkError(ErrStat, ErrMsg);

     }

     timeZero = false;

     if (scStatus) {
       sc.calcOutputs_n(0.0, scInputsGlob_n, scInputsTurbine_n, scOutputsGlob_n, scOutputsTurbine_n);
       fastScInputOutput(scInputsTurbine_from_FAST, scInputsTurbine_n, scOutputsTurbine_n, scOutputsTurbine_to_FAST);
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
       sc.updateStates(nt_global * dtFAST, scInputsGlob_n, scInputsTurbine_n); // Predict state at 'n+1' based on inputs
       sc.calcOutputs_np1( (nt_global + 1) * dtFAST, scInputsGlob_np1, scInputsTurbine_np1, scOutputsGlob_np1, scOutputsTurbine_np1);
       fastScInputOutput(scInputsTurbine_from_FAST, scInputsTurbine_np1, scOutputsTurbine_np1, scOutputsTurbine_to_FAST);
   }

   nt_global = nt_global + 1;
   
   if(scStatus) {
     scAdvanceTime(); // Advance states, inputs and outputs from 'n' to 'n+1'
   }
  
  if ( (((nt_global - ntStart) % nEveryCheckPoint) == 0 )  && (nt_global != ntStart) ) {
    //sprintf(CheckpointFileRoot, "../../CertTest/Test18.%d", nt_global);
    for (int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
      CheckpointFileRoot[iTurb] = " "; // if blank, it will use FAST convention <RootName>.nt_global
      FAST_CreateCheckpoint(&iTurb, CheckpointFileRoot[iTurb].data(), &ErrStat, ErrMsg);
      checkError(ErrStat, ErrMsg);
    }
    if(scStatus) {
      if (fastMPIRank == 0) {
          sc.writeRestartFile(nt_global);
      }
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
       sc.updateStates( nt_global * dtFAST, scInputsGlob_n, scInputsTurbine_n); // Predict state at 'n+1' based on inputs
       sc.calcOutputs_np1( (nt_global+1) * dtFAST, scInputsGlob_np1, scInputsTurbine_np1, scOutputsGlob_np1, scOutputsTurbine_np1);
       fastScInputOutput(scInputsTurbine_from_FAST, scInputsTurbine_np1, scOutputsTurbine_np1, scOutputsTurbine_to_FAST);
   }

   nt_global = nt_global + 1;
   
   if(scStatus) {
       scAdvanceTime(); // Advance states, inputs and outputs from 'n' to 'n+1'
   }
  
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

void fast::OpenFAST::getApproxHubPos(std::vector<double> & currentCoords, int iTurbGlob) {

  // Get hub position of Turbine 'iTurbGlob'
  currentCoords[0] = globTurbineData[iTurbGlob].TurbineHubPos[0];
  currentCoords[1] = globTurbineData[iTurbGlob].TurbineHubPos[1];
  currentCoords[2] = globTurbineData[iTurbGlob].TurbineHubPos[2];

}

void fast::OpenFAST::getHubPos(std::vector<double> & currentCoords, int iTurbGlob) {

  // Get hub position of Turbine 'iTurbGlob'
  int iTurbLoc = get_localTurbNo(iTurbGlob);
  currentCoords[0] = cDriver_Input_from_FAST[iTurbLoc].pxVel[0] + TurbineBasePos[iTurbLoc][0] ;
  currentCoords[1] = cDriver_Input_from_FAST[iTurbLoc].pyVel[0] + TurbineBasePos[iTurbLoc][1] ;
  currentCoords[2] = cDriver_Input_from_FAST[iTurbLoc].pzVel[0] + TurbineBasePos[iTurbLoc][2] ;
  
}

void fast::OpenFAST::getHubShftDir(std::vector<double> & hubShftVec, int iTurbGlob) {

  // Get hub shaft direction of current turbine - pointing downwind
  int iTurbLoc = get_localTurbNo(iTurbGlob);
  hubShftVec[0] = cDriver_Input_from_FAST[iTurbLoc].pOrientation[0] ;
  hubShftVec[1] = cDriver_Input_from_FAST[iTurbLoc].pOrientation[3] ;
  hubShftVec[2] = cDriver_Input_from_FAST[iTurbLoc].pOrientation[6] ;

}


void fast::OpenFAST::getVelNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob) {

  // Set coordinates at current node of current turbine 
  int iTurbLoc = get_localTurbNo(iTurbGlob);
  for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbLoc);
  currentCoords[0] = cDriver_Input_from_FAST[iTurbLoc].pxVel[iNode] + TurbineBasePos[iTurbLoc][0] ;
  currentCoords[1] = cDriver_Input_from_FAST[iTurbLoc].pyVel[iNode] + TurbineBasePos[iTurbLoc][1] ;
  currentCoords[2] = cDriver_Input_from_FAST[iTurbLoc].pzVel[iNode] + TurbineBasePos[iTurbLoc][2] ;
  
}

void fast::OpenFAST::getForceNodeCoordinates(std::vector<double> & currentCoords, int iNode, int iTurbGlob) {

  // Set coordinates at current node of current turbine 
  int iTurbLoc = get_localTurbNo(iTurbGlob);
  currentCoords[0] = cDriver_Input_from_FAST[iTurbLoc].pxForce[iNode] + TurbineBasePos[iTurbLoc][0] ;
  currentCoords[1] = cDriver_Input_from_FAST[iTurbLoc].pyForce[iNode] + TurbineBasePos[iTurbLoc][1] ;
  currentCoords[2] = cDriver_Input_from_FAST[iTurbLoc].pzForce[iNode] + TurbineBasePos[iTurbLoc][2] ;

}

void fast::OpenFAST::getForceNodeOrientation(std::vector<double> & currentOrientation, int iNode, int iTurbGlob) {

  // Set orientation at current node of current turbine 
  int iTurbLoc = get_localTurbNo(iTurbGlob);
  for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
  for(int i=0;i<9;i++) {
    currentOrientation[i] = cDriver_Input_from_FAST[iTurbLoc].pOrientation[iNode*9+i] ;
  }

}

void fast::OpenFAST::getForce(std::vector<double> & currentForce, int iNode, int iTurbGlob) {

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

void fast::OpenFAST::setVelocity(std::vector<double> & currentVelocity, int iNode, int iTurbGlob) {

  // Set velocity at current node of current turbine - 
  int iTurbLoc = get_localTurbNo(iTurbGlob);
  for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numVelPtsLoc(iTurbLoc);
  cDriver_Output_to_FAST[iTurbLoc].u[iNode] = currentVelocity[0];
  cDriver_Output_to_FAST[iTurbLoc].v[iNode] = currentVelocity[1];
  cDriver_Output_to_FAST[iTurbLoc].w[iNode] = currentVelocity[2];
}

void fast::OpenFAST::setVelocityForceNode(std::vector<double> & currentVelocity, int iNode, int iTurbGlob) {

  // Set velocity at current node of current turbine - 
  int iTurbLoc = get_localTurbNo(iTurbGlob);
  for(int j=0; j < iTurbLoc; j++) iNode = iNode - get_numForcePtsLoc(iTurbLoc);
  forceNodeVel[iTurbLoc][iNode][0] = currentVelocity[0];
  forceNodeVel[iTurbLoc][iNode][1] = currentVelocity[1];
  forceNodeVel[iTurbLoc][iNode][2] = currentVelocity[2];
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
	rDistForce[j] = sqrt( 
                             (cDriver_Input_from_FAST[iTurb].pxForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pxForce[0])*(cDriver_Input_from_FAST[iTurb].pxForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pxForce[0])  
		           + (cDriver_Input_from_FAST[iTurb].pyForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pyForce[0])*(cDriver_Input_from_FAST[iTurb].pyForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pyForce[0])  
		           + (cDriver_Input_from_FAST[iTurb].pzForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pzForce[0])*(cDriver_Input_from_FAST[iTurb].pzForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pzForce[0])  			
			    );
      }

      // Interpolate to the velocity nodes
      int nVelPtsBlade = get_numVelPtsBladeLoc(iTurb);
      for(int j=0; j < nVelPtsBlade; j++) {
	int iNodeVel = 1 + iBlade * nVelPtsBlade + j ; //Assumes the same number of velocity (Aerodyn) nodes for all blades
	double rDistVel = sqrt( 
			      (cDriver_Input_from_FAST[iTurb].pxVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pxVel[0])*(cDriver_Input_from_FAST[iTurb].pxVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pxVel[0])  
		            + (cDriver_Input_from_FAST[iTurb].pyVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pyVel[0])*(cDriver_Input_from_FAST[iTurb].pyVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pyVel[0])  
		            + (cDriver_Input_from_FAST[iTurb].pzVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pzVel[0])*(cDriver_Input_from_FAST[iTurb].pzVel[iNodeVel] - cDriver_Input_from_FAST[iTurb].pzVel[0])  			
			      );
	//Find nearest two force nodes
	int jForceLower = 0;
	while ( (rDistForce[jForceLower+1] < rDistVel) && ( jForceLower < (nForcePtsBlade-2)) )   {
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
	hDistForce[j] = sqrt( 
			     (cDriver_Input_from_FAST[iTurb].pxForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pxForce[iNodeBotTowerForce])*(cDriver_Input_from_FAST[iTurb].pxForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pxForce[iNodeBotTowerForce])  
                           + (cDriver_Input_from_FAST[iTurb].pyForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pyForce[iNodeBotTowerForce])*(cDriver_Input_from_FAST[iTurb].pyForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pyForce[iNodeBotTowerForce])  
			   + (cDriver_Input_from_FAST[iTurb].pzForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pzForce[iNodeBotTowerForce])*(cDriver_Input_from_FAST[iTurb].pzForce[iNodeForce] - cDriver_Input_from_FAST[iTurb].pzForce[iNodeBotTowerForce])	
			    );
      }
      
      
      int iNodeBotTowerVel = 1 + nBlades * get_numVelPtsBladeLoc(iTurb); // Assumes the same number of velocity (Aerodyn) nodes for all blades
      for(int j=0; j < nVelPtsTower; j++) {
	int iNodeVel = iNodeBotTowerVel + j ; 
	double hDistVel = sqrt( 
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

            relLoc[0] = cDriver_Input_from_FAST[iTurbLoc].pxForce[iNode] - cDriver_Input_from_FAST[iTurbLoc].pxForce[0] ;
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
    if ( (iNode + 1 - (get_numVelPts(iTurbLoc) - get_numVelPtsTwr(iTurbLoc)) ) > 0) {
      return TOWER; 
    }
    else {
      return BLADE;
    }
  }
  else {
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
    if ( (iNode + 1 - (get_numForcePts(iTurbLoc) - get_numForcePtsTwr(iTurbLoc)) ) > 0) {
      return TOWER; 
    }
    else {
      return BLADE;
    }
  }
  else {
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
  numBlades.resize(nTurbinesProc);
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
    numForcePtsBlade[iTurb] = globTurbineData[globProc].numForcePtsBlade;
    numForcePtsTwr[iTurb] = globTurbineData[globProc].numForcePtsTwr;

  }

  // Allocate memory for Turbine datastructure for all turbines
  FAST_AllocateTurbines(&nTurbinesProc, &ErrStat, ErrMsg);
  
  // Allocate memory for OpFM Input types in FAST
  cDriver_Input_from_FAST.resize(nTurbinesProc) ;
  cDriver_Output_to_FAST.resize(nTurbinesProc) ;
  
  scInputsTurbine_from_FAST.resize(nTurbinesProc) ;
  scOutputsTurbine_to_FAST.resize(nTurbinesProc) ;

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
  }
  
  MPI_Group_free(&fastMPIGroup);
  if (MPI_COMM_NULL != fastMPIComm) {
    MPI_Comm_free(&fastMPIComm);
  }
  MPI_Group_free(&worldMPIGroup);
  
  if(scStatus) {
      sc.end();
   
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

    scStatus = fi.scStatus;

    numScInputsTurbine = fi.numScInputsTurbine;
    numScOutputsTurbine = fi.numScOutputsTurbine;
    numScInputsGlob = fi.numScInputsGlob;
    numScOutputsGlob = fi.numScOutputsGlob;

    if ( numScOutputsTurbine > 0 ) {

        scOutputsTurbine_nm1.resize(nTurbinesGlob*numScOutputsTurbine) ;
        scOutputsTurbine_n.resize(nTurbinesGlob*numScOutputsTurbine) ;
        scOutputsTurbine_np1.resize(nTurbinesGlob*numScOutputsTurbine) ;
        for (int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
            for(int iOutput=0; iOutput < numScOutputsTurbine; iOutput++) {
                scOutputsTurbine_nm1[iTurb*numScOutputsTurbine + iOutput] = 0.0 ; // Initialize to zero
                scOutputsTurbine_n[iTurb*numScOutputsTurbine + iOutput] = 0.0 ; // Initialize to zero
                scOutputsTurbine_np1[iTurb*numScOutputsTurbine + iOutput] = 0.0 ; // Initialize to zero
            }
        }
    }

    if ( numScInputsTurbine > 0 ) {

        scInputsTurbine_nm1.resize(nTurbinesGlob*numScInputsTurbine) ;
        scInputsTurbine_n.resize(nTurbinesGlob*numScInputsTurbine) ;
        scInputsTurbine_np1.resize(nTurbinesGlob*numScInputsTurbine) ;
        for (int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
            for(int iInput=0; iInput < numScInputsTurbine; iInput++) {
                scInputsTurbine_nm1[iTurb*numScInputsTurbine + iInput] = 0.0 ; // Initialize to zero
                scInputsTurbine_n[iTurb*numScInputsTurbine + iInput] = 0.0 ; // Initialize to zero
                scInputsTurbine_np1[iTurb*numScInputsTurbine + iInput] = 0.0 ; // Initialize to zero
            }
        }
    }

    if ( numScOutputsGlob > 0 ) {
    
        scOutputsGlob_nm1.resize(numScOutputsGlob) ;
        scOutputsGlob_n.resize(numScOutputsGlob) ;
        scOutputsGlob_np1.resize(numScOutputsGlob) ;
        for(int iOutput=0; iOutput < numScOutputsGlob; iOutput++) {
            scOutputsGlob_nm1[iOutput] = 0.0 ; // Initialize to zero
            scOutputsGlob_n[iOutput] = 0.0 ; // Initialize to zero
            scOutputsGlob_np1[iOutput] = 0.0 ; // Initialize to zero
        }
    }

    if ( numScInputsGlob > 0 ) {

        scInputsGlob_nm1.resize(numScInputsGlob) ;
        scInputsGlob_n.resize(numScInputsGlob) ;
        scInputsGlob_np1.resize(numScInputsGlob) ;
        for(int iInput=0; iInput < numScInputsGlob; iInput++) {
            scInputsGlob_nm1[iInput] = 0.0 ; // Initialize to zero
            scInputsGlob_n[iInput] = 0.0 ; // Initialize to zero
            scInputsGlob_np1[iInput] = 0.0 ; // Initialize to zero
        }
    }

    scInData sci ;
    sci.nInputsTurbine = numScInputsTurbine ;
    sci.nOutputsTurbine = numScOutputsTurbine ;
    sci.nInputsGlob = numScInputsGlob ;
    sci.nOutputsGlob = numScOutputsGlob ;
    sci.nGlobStates = fi.numScGlobStates ;
    sci.nTurbineStates = fi.numScTurbineStates ;
    sci.scLibFile = fi.scLibFile ;
    sc.load(sci);

  } else {

    scStatus = false;
    numScInputsGlob = 0;
    numScOutputsGlob = 0;
    numScInputsTurbine = 0;
    numScOutputsTurbine = 0;
  }

}

void fast::OpenFAST::fastScInputOutput(std::vector<SC_InputType_t> & scIT_from_FAST, std::vector<float> & scIT, std::vector<float> & scOT, std::vector<SC_OutputType_t> & scOT_to_FAST) {
  
    // Transfers
    // scIT <------ scIT_from_FAST 
    // scOT_to_FAST <------- scOT
    
    for(int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
        for(int iInput=0; iInput < numScInputsTurbine; iInput++) {
            scIT[iTurb*numScInputsTurbine + iInput] = 0.0; // Initialize to zero 
        }
    }
    
    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        for(int iInput=0; iInput < numScInputsTurbine; iInput++) {
            scIT[turbineMapProcToGlob[iTurb]*numScInputsTurbine + iInput] = scIT_from_FAST[iTurb].toSC[iInput] ;
        }
    }
    
    if (MPI_COMM_NULL != fastMPIComm) {
        MPI_Allreduce(MPI_IN_PLACE, scIT.data(), numScInputsTurbine*nTurbinesGlob, MPI_DOUBLE, MPI_SUM, fastMPIComm) ;
    }
    
    
    for(int iTurb=0; iTurb < nTurbinesProc; iTurb++) {
        for(int iOutput=0; iOutput < numScOutputsTurbine; iOutput++) {
            scOT_to_FAST[iTurb].fromSC[iOutput] = scOT[turbineMapProcToGlob[iTurb]*numScOutputsTurbine + iOutput] ;
        }
    }
    
}


void fast::OpenFAST::scAdvanceTime() {

    for(int iTurb=0; iTurb < nTurbinesGlob; iTurb++) {
        for(int iInput=0; iInput < numScInputsTurbine; iInput++) {
            scInputsTurbine_nm1[iTurb*numScInputsTurbine + iInput] = scInputsTurbine_n[iTurb*numScInputsTurbine + iInput];
            scInputsTurbine_n[iTurb*numScInputsTurbine + iInput] = scInputsTurbine_np1[iTurb*numScInputsTurbine + iInput];
//            scInputsTurbine_np1[iTurb*numScInputsTurbine + iInput] = Predictor?
        }
        for(int iOutput=0; iOutput < numScOutputsTurbine; iOutput++) {
            scOutputsTurbine_nm1[iTurb*numScOutputsTurbine + iOutput] = scOutputsTurbine_n[iTurb*numScOutputsTurbine + iOutput];
            scOutputsTurbine_n[iTurb*numScOutputsTurbine + iOutput] = scOutputsTurbine_np1[iTurb*numScOutputsTurbine + iOutput];
        }
    }
    for(int iInput=0; iInput < numScInputsGlob; iInput++) {
        scInputsGlob_nm1[iInput] = scInputsGlob_n[iInput];
        scInputsGlob_n[iInput] = scInputsGlob_np1[iInput];
        //scInputsGlob_np1[iInput] = Predictor?
    }
    for(int iOutput=0; iOutput < numScOutputsGlob; iOutput++) {
        scOutputsGlob_nm1[iOutput] = scOutputsGlob_n[iOutput];
        scOutputsGlob_n[iOutput] = scOutputsGlob_np1[iOutput];
        //scOutputsGlob_np1[iOutput] = Predictor?
    }
    
    sc.advanceStates();
}

