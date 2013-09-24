/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        Numerics.cpp
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  	     
   Copyright Sept. 2013
   
   Author: Marco D. Masciola, 
   National Renewable Energy Laboratory, Golden, Colorado, USA
   
   This file is part of the Mooring Analysis Program (MAP).
   
   Licensed to the Apache Software Foundation (ASF) under one
   or more contributor license agreements.  See the NOTICE file
   distributed with this work for additional information
   regarding copyright ownership.  The ASF licenses this file
   to you under the Apache License, Version 2.0 (the
   "License"); you may not use this file except in compliance
   with the License.  You may obtain a copy of the License at
   
   http://www.apache.org/licenses/LICENSE-2.0
   
   Unless required by applicable law or agreed to in writing,
   software distributed under the License is distributed on an
   "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
   KIND, either express or implied.  See the License for the
   specific language governing permissions and limitations
   under the License.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
*/


#include "MAP_OtherStateType.h"

extern PetscErrorCode ResidualFunction( SNES , Vec , Vec , void* );
extern PetscErrorCode FormJacobian    ( SNES , Vec , Mat* , Mat* , MatStructure* , void* );



struct LoggerObject {
  std::string strLog;
};


std::ostream& 
operator<<( std::ostream       &ostr , 
            const LoggerObject &in   ) 
{
  ostr << in.strLog;
  return ostr;
}


void 
MPI_InitCallOnce( )
{
  static class CallOnce
  { 
  public: 
    CallOnce( ) { 
      int  num  = 0;
      char **arr = NULL;
      MPI_Init( &num, &arr );
    }
  } CallOnce_;
}

/**
 * Return a PETSc buffer to the MAP logger. Refer to PETSc documentation to sort out 
 * the details. The function LoggerObject and macro LOG_PETSC_INFO are designed into 
 * MAP (and are not part of the PETSc library). If DEBUG is decalred as a compile option,
 * then the logger will be written to std::cout. 
 *
 * @see  Numerics::InitializeSolver()
 * @see  MAP_SummaryLogger
 * @see  #define LOG_PETSC_INFO
 */
PetscErrorCode 
MAP_PETScPrintf( FILE       *fd      , 
                 const char format[] , 
                 va_list    argp     ) 
{
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  if (fd != stdout && fd != stderr) { 
    ierr = PetscVFPrintfDefault(fd,format,argp); CHKERRQ(ierr);
  } else {
    int BIG = 180;
    char buff[180];
    size_t length;
    
    ierr = PetscVSNPrintf(buff,BIG,format,&length,argp);CHKERRQ(ierr);   
    
    LoggerObject strPETSc = { std::string(buff) };
    LOG_PETSC_INFO( strPETSc );
  }
  PetscFunctionReturn(0);
}


/**
 * Recieve indidual words in the MAP input file options string
 *
 * @param  std::string  optionStr is a word related to one of the PETSc run-time options
 * @see    MAP_OtherStateType_class::SetSolverOptions( )
 * @see    MAP_Init( )
 */
void Numerics::
setNumericsOptionsString( const std::string &optionStr )
{
  options_string.push_back( optionStr );
};


/**
 * Initialize the numeric solver. Create data structures, allocate memory, determine which type
 * of Jacobian to use (analytical, finite different), ect...
 *
 * @param  other  MAP other states. This is used to access the UserData structure
 * @param  init   MAP initialization state. Not currently used at this time (and problably will never be needed)
 * @param  err
 * @param  msg
 */
int Numerics::
InitializeSolver( MAP_OtherStateType_class &other ,
                  MAP_InitInputType_class  &init  ,
                  MAP_ErrStat_class        &err   ,  
                  MAP_Message_class        &msg   )
{
  PetscScalar PETSc_number      = 0.0;	
  PetscReal   msqs_scaling      = 1.0;
  PetscBool   msqs_solver_types = PETSC_FALSE;
  PetscBool   permute           = PETSC_FALSE;
  PetscBool   msqs_fd_jacobian  = PETSC_FALSE;  
  int         num_eq            = other.user_data.sizeOfConstraint();
  int         arg_num           = options_string.size()+1;
  char        solve_option[256] = "none";
  char        fd_coloring[256]  = "none";
  char        ordering[256]     = MATORDERINGRCM;
  std::string help_str          = "";
  char        **arr             = new char*[ arg_num ];
  
  arr[0] = new char[ 1 ];
  strcpy( arr[0] , "" );
  
  for(unsigned int i=1 ; i<options_string.size()+1 ; i++ ){
    arr[i] = new char[ options_string[i-1].size() + 1 ];
    strcpy( arr[i] , options_string[i-1].c_str() );
  }

#ifndef DEBUG
  PetscVFPrintf = MAP_PETScPrintf;
#endif // if DEBUG isnt defined, the stdout is to the map summary file. If DEBUG is defined, 

  MPI_InitCallOnce( );
  PetscInitialize( &arg_num, &arr, (char *)0 ,"" ); 
  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);
  
  for( unsigned int i = 0 ; i<options_string.size()+1 ; i++ ) { delete [] arr[i]; }
  
  delete [] arr;
  
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"MAP MSQS options","SNES");CHKERRQ(ierr);
  {
    this->msqs_k = PETSC_FALSE;    

    // Are we using finite differencing or the analytical jacobian to solve the problem?
    help_str = "Compute the Jacobian using central finite-differencing. Coloring can be enabled with <-msqs_jacobian fdcolor>.";
    ierr = PetscOptionsList("-msqs_fd_jacobian", help_str.c_str(), "color",
                            MatOrderingList, fd_coloring, fd_coloring, sizeof(fd_coloring), &msqs_fd_jacobian );

    // if '-default' is in the MAP input file, then use the default solver settings
    help_str = "Set the convergence tolerance for the post check.";
    ierr = PetscOptionsReal( "-msqs_tol", help_str.c_str(), "", this->msqs_tol, &this->msqs_tol, NULL );CHKERRQ(ierr);

    // if '-default' is in the MAP input file, then use the default solver settings
    help_str = "Scaling factor for the Newton equation residuals. Chaning this value may help convergence";
    ierr = PetscOptionsReal( "-msqs_scaling", help_str.c_str(), "", msqs_scaling, &msqs_scaling, NULL );CHKERRQ(ierr);
    other.user_data.SetFunctionScaling( msqs_scaling );

    help_str = "Generate the linearized stiffness matrix to the MAP summary file.";
    ierr = PetscOptionsBool("-msqs_k", help_str.c_str(), 0, this->msqs_k, &this->msqs_k, 0 );CHKERRQ(ierr);

    help_str = "Premultiply the Jacobian by a permutation matrix to solve the problem under a new ordering. Changing the ordering may assist with ill-conditioned problems.";
    ierr = PetscOptionsList("-msqs_permute", help_str.c_str(), "natural, nd, 1wd, rcm, qmd, rowlength",
                            MatOrderingList, ordering, ordering, sizeof(ordering), &permute );

    help_str = "MSQS default solver options. The default options are line search (ls) and trust region (tr).";
    ierr = PetscOptionsList("-msqs_default", help_str.c_str(), "tr,ls",
                            MatOrderingList, solve_option, solve_option, sizeof(solve_option), &msqs_solver_types );

    // -help. If '-help' is in the input file, MAP won't solve
    ierr = PetscOptionsHasName( PETSC_NULL , "-help" , &this->help_flag );CHKERRQ( ierr );
  }
  PetscOptionsEnd();

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ( ierr );   // Create nonlinear solver context
  ierr = VecCreate( PETSC_COMM_WORLD, &x );CHKERRQ(ierr);       // Create matrix and vector data structures; set corresponding routines
  ierr = VecSetSizes( x, PETSC_DECIDE, num_eq );CHKERRQ(ierr);	// VecSetSizes( the vector , the local size , the global size )
  ierr = VecSetFromOptions( x );CHKERRQ(ierr);	                // Configures the vector from the options database.
  ierr = VecDuplicate( x, &r );CHKERRQ(ierr);                   // Creates a new vector of the same type as an existing vector.
	
  // Create Jacobian matrix data structure
  ierr = MatCreate( PETSC_COMM_WORLD, &J );CHKERRQ(ierr);
  ierr = MatSetSizes( J, PETSC_DECIDE, PETSC_DECIDE, num_eq, num_eq );CHKERRQ(ierr);
  ierr = MatSetFromOptions( J );CHKERRQ(ierr);
  ierr = MatSetUp( J );CHKERRQ(ierr);

  // Set function evaluation routine and vector.
  ierr = SNESSetFunction( snes, r, ResidualFunction, static_cast<void*>(&other.user_data) );CHKERRQ(ierr);
	
  // Customize nonlinear solver; set runtime options. Set linear solver defaults for this problem. By extracting the
  // KSP, KSP, and PC contexts from the SNES context, we can then directly call any KSP, KSP, and PC routines to set various options.
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

  // Default MAP solver options. Other (customized) solver options are possible as well.
  if (msqs_solver_types) { 
    if ( strcmp(solve_option,"tr")==0 ) {
      // if using default settings, use the trust region
      ierr = KSPSetType(ksp, KSPPREONLY );CHKERRQ(ierr); 
      ierr = PCSetType(pc, PCLU );CHKERRQ(ierr); 
      ierr = SNESSetType(snes, SNESNEWTONTR );CHKERRQ(ierr); 
      PCFactorReorderForNonzeroDiagonal( pc, 1e-10);
    } else if ( strcmp(solve_option,"ls")==0 ) {
      // if using default settings, use the trust region
      ierr = KSPSetType(ksp, KSPPREONLY );CHKERRQ(ierr); 
      ierr = PCSetType(pc, PCLU );CHKERRQ(ierr); 
      ierr = SNESSetType(snes, SNESNEWTONLS );CHKERRQ(ierr); 
      PCFactorReorderForNonzeroDiagonal( pc, 1e-10);
    } else {
      // if using default settings, use the trust region
      ierr = KSPSetType(ksp, KSPPREONLY );CHKERRQ(ierr); 
      ierr = PCSetType(pc, PCLU );CHKERRQ(ierr); 
      ierr = SNESSetType(snes, SNESNEWTONTR );CHKERRQ(ierr); 
      PCFactorReorderForNonzeroDiagonal( pc, 1e-10);
      std::string str = "";
      MAPSetUniversalErrorStat( MAP_WARNING_6, str, err, msg );
    }
  } else {
    // if using custom options
    ierr = KSPSetFromOptions( ksp );CHKERRQ(ierr);
    ierr = PCSetFromOptions( pc );CHKERRQ(ierr);
    ierr = SNESSetFromOptions( snes );CHKERRQ(ierr);
  } 

  // Reordering options are natural, nd, 1wd, rcm, qmd, rowlength"
  if (permute) { 
    if ( strcmp(ordering,"natural")==0 ) {
      PCFactorSetMatOrderingType( pc, MATORDERINGNATURAL );
    } else if ( strcmp(ordering,"nd")==0 ) {
      PCFactorSetMatOrderingType( pc, MATORDERINGND );
    } else if ( strcmp(ordering,"1wd")==0 ) {
      PCFactorSetMatOrderingType( pc, MATORDERING1WD );
    } else if ( strcmp(ordering,"rcm")==0 ) {
      PCFactorSetMatOrderingType( pc, MATORDERINGRCM );
    } else if ( strcmp(ordering,"qmd")==0 ) {
      PCFactorSetMatOrderingType( pc, MATORDERINGQMD );
    } else if ( strcmp(ordering,"rowlength")==0 ) {
      PCFactorSetMatOrderingType( pc, MATORDERINGROWLENGTH );
    } else {
      std::string str = "";
      PCFactorSetMatOrderingType( pc, MATORDERINGRCM );
      MAPSetUniversalErrorStat( MAP_WARNING_7, str, err, msg );    
    } 
  } else {
    PCFactorSetMatOrderingType( pc, ordering );
  }

  // Set Jacobian matrix data structure and Jacobian evaluation routine
  if ( msqs_fd_jacobian ) {
    if ( strcmp(fd_coloring,"color")==0 ) { 
      // the Jacobian is calculated using finite difference
      ISColoring    iscoloring;
      MatFDColoring matfdcoloring;
      MatStructure  matstruct;
      
      other.user_data.initializeJacobian( );
      FormJacobian( snes, x ,&J, &J, &matstruct, static_cast<void*>(&other.user_data) );
      
      MatGetColoring( J, MATCOLORINGSL, &iscoloring );
      MatFDColoringCreate( J, iscoloring, &matfdcoloring );
      ISColoringDestroy( &iscoloring );
      
      MatFDColoringSetFunction( matfdcoloring, (PetscErrorCode (*)(void))ResidualFunction, static_cast<void*>(&other.user_data) );
      MatFDColoringSetFromOptions( matfdcoloring );
      
      ierr = SNESSetJacobian( snes, J, J, SNESComputeJacobianDefaultColor, matfdcoloring );CHKERRQ(ierr);    
      MatFDColoringDestroy( &matfdcoloring );
    } else if ( strcmp(fd_coloring,"")==0 ) {
      ierr = SNESSetJacobian( snes, J, J, SNESComputeJacobianDefault, PETSC_NULL );CHKERRQ(ierr);    
    } else {
      std::string str = "";
      MAPSetUniversalErrorStat( MAP_WARNING_8, str, err, msg );
      ierr = SNESSetJacobian( snes, J, J, SNESComputeJacobianDefault, PETSC_NULL );CHKERRQ(ierr);    
    }
  } else { 
    /*
        Computing the Jacobian analytically
        Set the non-zero structure of the Jacobian
        
        J = [  A     B ]
            [ -B^T   C ]
  
        We are only setting the non-zero structure for the A and B blocks. The B block handles the 
        coupling between the force-balance equations and catenary equations
    */
    ierr = SNESSetJacobian( snes, J, J, FormJacobian, static_cast<void*>(&other.user_data) );CHKERRQ(ierr);    
    other.user_data.initializeJacobian( );
  }
 
  /*
      Pass initial guess into vector x
  */
  for ( int i=0 ; i<num_eq ; i++ ) {
    PETSc_number = boost::lexical_cast<PetscScalar> ( other.user_data.getConstraint(i) );
    VecSetValues( x, 1, &i, &PETSc_number, INSERT_VALUES );
  }

  return 0;
};


/** 
 * Solve the unknowns such that the residual is minimized The user should initialize the 
 * vector, x, with the initial guess for the nonlinear solver prior to calling 
 * SNESSolve().In particular, to employ an initial guess of zero, the user should 
 * explicitly set this vector to zero by calling VecSet().
 *
 * @param  other  this is used to access the UserData class (information about the problem
 *                such as the Jacobian, element connectivity, ect.
 * @param  err
 * @param  msg
 */
int Numerics::
PetscSolve( MAP_OtherStateType_class &other ,
            MAP_ErrStat_class              &err , 
            MAP_Message_class              &msg   )
{
  try {
    ierr = SNESSolve( snes, PETSC_NULL, x ); 
    if ( ierr != 0 ) throw MAP_ERROR_86; // throw error if the solver fails 
 
    ierr = SNESGetConvergedReason(snes , &reason ); CHKERRQ(ierr);
  } catch ( MAP_ERROR_CODE &code ) {    
    std::string str = "";
    MAPSetUniversalErrorStat( code, str, err, msg );
    CHKERRQ(ierr);
    return 1;
  }

  // now make sure the residuals meet our target tolerance -msqs_tol
  if( other.CheckResidualConvergence( err , msg ) != 0 ){
    std::string str = "";
    MAPSetUniversalErrorStat( MAP_ERROR_69, str, err, msg );
    return 1;
  }  
  
  this->PetscConvergeReason( err, msg );
  
  return 0;
};


/**
 * For the analytical Jacobian for the current iterates of x_in
 *
 * @see     UserData::getJacB
 * @see     UserData::getJacA
 * @param   snes_in  non-linear solver context
 * @param   x_in     input variables (the things being iterated). This is what 
 *                   changes the Jacobian at each iteration
 * @param   jac      Jacobian matrix
 * @param   b        is ussually the same as the jac
 * @param   flag     Same non-zero structure ALWAYS!
 * @param   ctx      User context. This accesses the Jacobian information in 
 *                   UserData::getJacA and UserData::getJacB
 * @return  PETSc error code if things go nuts unexpectedly
 */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian( SNES         snes_in , 
                             Vec          x_in    , 
                             Mat          *jac    ,
			     Mat          *b      , 
                             MatStructure *flag   , 
                             void         *ctx    ) 
{
  UserData    *data = static_cast<UserData*>(ctx);    
  const int    M = data->getNumNodeEqs();    
  PetscScalar *xx = NULL;
  PetscScalar  D[4];
  PetscInt     idAx[2] = {0,1};    
  double       K = data->GetFunctionScaling(); 

  //std::cout << "Forming Jacobian\n\n";

  //MatView( *jac , PETSC_VIEWER_STDOUT_SELF );
  MatZeroEntries( *jac );
  
  // Get pointer to vector data
  VecGetArray( x_in , &xx );    
  
//  data->UserJacobianEvaluations( xx );
//  std::cout << "In Jacobian : " << xx[0] << "\n"
//            << "              " << xx[1] << "\n"
//            << "              " << xx[2] << "\n"
//            << "              " << xx[3] << "\n"
//            << "              " << xx[4] << "\n"
//            << "              " << xx[5] << "\n\n";
  
  for ( int i=0 ; i<data->getNumJacAEntries() ; i++ ) {
    MatSetValue( *jac, data->Ai(i), data->Aj(i), K*data->getJacA(i), INSERT_VALUES );
  }
  
  // now fill the B (off-diagonal) block matrix and it's transpose
  for( int i=0 ; i<data->getNumJacBEntries() ; i++ ) {
    MatSetValue( *jac,   data->Bj(i), M+data->Bi(i), -K*data->getJacB(i), INSERT_VALUES ); // upper right J quadrant
    MatSetValue( *jac, M+data->Bi(i),   data->Bj(i),    data->getJacB(i), INSERT_VALUES ); // lower left J quadrant
  }

  // form jacobian for the catenary equations with respect to H, V, (and Lu, but this isn't support yet)
  for( int i=0 ; i<data->sizeOfElement() ; i++ ){
    D[0] = data->getdXdH( i );
    D[1] = data->getdXdV( i );
    D[2] = data->getdZdH( i );
    D[3] = data->getdZdV( i );
    
    idAx[0] = M+2*i;
    idAx[1] = M+2*i+1;
    
    // set D block diagonals
    MatSetValues( *jac, 2, idAx, 2, idAx, D, INSERT_VALUES );
  }
  
  *flag = SAME_NONZERO_PATTERN;
  
  // Restore vector
  VecRestoreArray( x_in , &xx );
  
  // Assemble matrix
  MatAssemblyBegin( *jac, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd  ( *jac, MAT_FINAL_ASSEMBLY );
  
//    MatView( *jac , PETSC_VIEWER_STDOUT_SELF );
//    std::cin.get();
    
  return(0);
};


/**
 * Calls the FunctionResidual metho in UserData to determine the function residual based on 
 * newly updated iterates.
 *
 * @see     UserData::FunctionEvaluations( )
 * @param   snes  PETSc non-linear solver context
 * @param   x     variables being minimized; x is the current iterate
 * @param   f     function residual
 * @param   ctx   user context. This is casted into UserData.
 * @return  PETSc error code if things fail
 */
#undef __FUNCT__
#define __FUNCT__ "ResidualFunction"
PetscErrorCode 
ResidualFunction( SNES snes , 
                  Vec  x    , 
                  Vec  f    , 
                  void *ctx )
{
  PetscErrorCode    ierr;
  const PetscScalar *xx;
  PetscScalar       *ff;
  MPI_Comm          comm;
  PetscMPIInt       size,rank;
  
  UserData *data = static_cast<UserData*>(ctx);
    
  ierr = PetscObjectGetComm((PetscObject)snes,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

  // Get pointers to vector data.
  //   - For default PETSc vectors, VecGetArray() returns a pointer to
  //     the data array.  Otherwise, the routine is implementation dependent.
  //   - You MUST call VecRestoreArray() when you no longer need access to
  //      the array.
  ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
  ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

  // this is the residual
  ff = data->UserFunctionEvaluations( ff, xx );
  
  // Restore vectors
  ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
  ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr); 
    
  return 0;
};


/**
 * Evaluate the function residuals and return them to the scalar FF
 *
 * @param  FF  array of scalars; the function residuals
 * @param  XX  the current state iterate; what is being solved
 * @todo   throw an exception for inconsistent number of equations/unknowns
 * @todo   Add in a stacling factor such that the sum force euquaion magnitudes can be decreased
 */
double *UserData::
UserFunctionEvaluations( PetscScalar       *FF , 
                         const PetscScalar *XX )
{
  //std::cout << "Evaluating functions\n\n";
  // The first step is to copy the constraint varaibles into XX (the state vector)
  for( int i=0 ; i<this->sizeOfConstraint() ; i++ ){
    this->setConstraint( i , XX[i] );
  }

  // set sum_FX, sum_FY and sum_FZ to zero
  for ( int i=0 ; i<this->sizeOfNode( ) ; i++ ){
    this->node[i]->SetSumForceToZero( );
  }
  
  // set fairlead and anchor nodes to zero. In other words, we are initializing the model. This calls
  // Node::fairlead->SetSumForceToZero();
  // Node::anchor->SetSumForceToZero();
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    this->element[i]->ResetNodes();
  }
  
  // update psi, l and h for each element being iterated
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    this->element[i]->UpdateElement( *errPtr , *msgPtr );
  }
  
  int cnt = 0;
  double K = this->GetFunctionScaling(); 
  
  for ( int i=0 ; i<this->sizeOfNode() ; i++ ) {	
    // solve X direction Newton equation
    if (this->node[i]->GetXNewtonEquationFlag()==true) {
      FF[cnt] = K*this->node[i]->f_x(); 
      //FF[cnt] = this->node[i]->f_x(); 
      cnt++;
    }
    
    // solve Y direction Newton equation
    if (this->node[i]->GetYNewtonEquationFlag()==true){
      //FF[cnt] = K*this->node[i]->f_y();
      FF[cnt] = this->node[i]->f_y();
      cnt++;
    }
    
    // solve Z direction Newton equation
    if (this->node[i]->GetZNewtonEquationFlag()==true){
      FF[cnt] = K*this->node[i]->f_z();
      //FF[cnt] = this->node[i]->f_z();
      cnt++;
    }
  }
  
  // Compute function :
  // For each element, the horizontal and vertical catenary equation is solved.  
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    FF[cnt] = this->element[i]->f_h(); // horizontal catenary eq. 
    cnt++;
    
    FF[cnt] = this->element[i]->f_v(); // vertical catenary eq. 
    cnt++;
  }
  
  // Make sure the number of equations is equal to the number of unknowns 
  assert( this->sizeOfConstraint() == cnt ); 
  
  return FF;
};


/** 
 * Converged:
 *   - SNES_CONVERGED_FNORM_ABS          :  2, ||F|| < atol 
 *   - SNES_CONVERGED_FNORM_RELATIVE     :  3, ||F|| < rtol*||F_initial|| 
 *   - SNES_CONVERGED_SNORM_RELATIVE     :  4, Newton computed step size small; || delta x || < stol 
 *   - SNES_CONVERGED_ITS                :  5, maximum iterations reached 
 *   - SNES_CONVERGED_TR_DELTA           :  7,
 * 
 * Diverged
 *   - SNES_DIVERGED_FUNCTION_DOMAIN     : -1, the new x location passed the function is not in the 
 *                                             domain of F
 *   - SNES_DIVERGED_FUNCTION_COUNT      : -2,
 *   - SNES_DIVERGED_LINEAR_SOLVE        : -3, the linear solve failed 
 *   - SNES_DIVERGED_FNORM_NAN           : -4,
 *   - SNES_DIVERGED_MAX_IT              : -5,
 *   - SNES_DIVERGED_LINE_SEARCH         : -6, the line search failed 
 *   - SNES_DIVERGED_INNER               : -7, inner solve failed 
 *   - SNES_DIVERGED_LOCAL_MIN           : -8, || J^T b || is small, implies converged to local 
 *                                             minimum of F()
 *   - SNES_CONVERGED_ITERATING          :  0 SNESConvergedReason;
 */
void Numerics::
PetscConvergeReason( MAP_ErrStat_class &err , 
                     MAP_Message_class &msg )
{
  try {
    switch( reason ){
    case 0 :                                                   
      msg.WriteConvergeReason("Converged (PETSc code 0).");            
      break;                                                   
    case 2 :                                                   
      msg.WriteConvergeReason("Converged (PETSc code 2: '||F|| < atol ').");
      break;                                                   
    case 3 :                                                   
      msg.WriteConvergeReason("Converged (PETSc code 3: '||F|| < rtol*||F_initial|| ').");
      break;                                                   
    case 4 :                                                   
      msg.WriteConvergeReason("Converged (PETSc code 4: 'Step size small; ||delta x|| < stol ').");
      break;                                                   
    case 5 :                                                   
      msg.WriteConvergeReason("Converged (PETSc code 5: 'Maximum iteration reached').");
      break;                                                                               
    case 7 :                                                                               
      msg.WriteConvergeReason("Converged (PETSc code 7).");
      break;                                                   
    case -1 :                                                                           
      throw MAP_ERROR_57;                                                               
      break;                                                                            
    case -2 :
      throw MAP_ERROR_58;                                                               
      break;                                                                            
    case -3 :
      throw MAP_ERROR_59;                                                               
      break;                                                                            
    case -4 :
      throw MAP_ERROR_60;                                                               
      break;                                                                            
    case -5 :
      throw MAP_ERROR_61;                                                               
      break;                                                                            
    case -6 :
      throw MAP_ERROR_62;                                                               
      break;                                                                            
    case -7 :
      // I am aware this case is contrived and doesn't really exist in practice         
      throw MAP_ERROR_63;                                                               
      break;                                                                            
    case -8 :
      throw MAP_ERROR_64;                                                               
      break;
    default :
      msg.WriteConvergeReason("MAP failed to converge.");
    }
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    MAPSetUniversalErrorStat( code, str, err, msg );
  }
};


/**
 * Clears all PETSc data structures and deallocate memory
 *
 * @param  MAP_ErrStat_class  err
 * @param  MAP_Message_class  msg
 */
int Numerics::
PetscEnd( MAP_ErrStat_class &err , 
          MAP_Message_class &msg )
{
    ierr = SNESGetIterationNumber( snes, &its );CHKERRQ(ierr);

    // Free work space.  PETSc variables/objects are destroyed when the program is done
    ierr = VecDestroy(&x);CHKERRQ(ierr); 
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = MatDestroy(&J);CHKERRQ(ierr);  
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    ierr = PetscFinalize();

    return 0;
}
