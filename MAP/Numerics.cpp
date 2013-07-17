/**
 * ====================================================================================================
 *                              Numerics.cpp
 * ====================================================================================================
 *	     
 * Copyright Sept. 2012
 * 
 * Author: Marco D. Masciola, 
 * National Renewable Energy Laboratory, Golden, Colorado, USA
 *
 * This file is part of the Mooring Analysis Program (MAP).
 *
 * MAP is free software: you can redistribute it and/or modify it under the terms 
 * of the GNU General Public License as published by the Free Software Foundation, 
 * either version 3 of the License, or (at your option) any later version.
 *
 * MAP is distributed in the hope that it will be useful, but WITHOUT ANY 
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along 
 * with MAP. If not, see:
 * 
 * < http://www.gnu.org/licenses/>
 * ====================================================================================================
 */


#include "MAP_OtherStateType.h"

extern PetscErrorCode ResidualFunction( SNES , Vec , Vec , void* );
extern PetscErrorCode FormJacobian    ( SNES , Vec , Mat* , Mat* , MatStructure* , void* );


/**
 * ====================================================================================================
 * setNumericsOptionsString
 * ====================================================================================================
 */
void Numerics::setNumericsOptionsString( const std::string &P ){
  options_string.push_back( P );
};


/**
 * ====================================================================================================
 * InitializeSolver
 * ====================================================================================================
 */
int Numerics::InitializeSolver( MAP_OtherStateType_class &T     ,
				MAP_InitInputType_class  &Init  ,
				MAP_ErrStat        &Error ,  
				MAP_Message        &Msg ){

  PetscScalar PETSc_number = 0.0;	
  int num_eq = T.user_data.sizeOfConstraint();

  int arg_num = options_string.size()+1;
  char **arr = new char*[ arg_num ];

  arr[0] = new char[ 1 ];
  strcpy( arr[0] , "" );
    
  for(unsigned int i=1 ; i<options_string.size()+1 ; i++ ){
    arr[i] = new char[ options_string[i-1].size() + 1 ];
    strcpy( arr[i] , options_string[i-1].c_str() );
  };//end FOR

  PetscInitialize( &arg_num, &arr, (char *)0 ,"" ); 

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);

  for( unsigned int i = 0 ; i<options_string.size()+1 ; i++ ){
    delete [] arr[i];
  };// END for

  delete [] arr;
    
  ierr = PetscOptionsHasName( PETSC_NULL , "-help" , &flg ); CHKERRQ( ierr );

//  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"MAP Coupled options","SNES");CHKERRQ(ierr);
//  {
//    this->FD_jacobian = PETSC_FALSE;
//    this->ptype = 1;
//    ierr = PetscOptionsBool("-finite_difference_jacobian","Compute the Jacobian using central finite-differencing.",0,this->FD_jacobian,&this->FD_jacobian,0);CHKERRQ(ierr);
//    ierr = PetscOptionsInt("-split_type","0: 1: Field splits is enabled.",0,this->ptype,&this->ptype,0);CHKERRQ(ierr);
//  }
//  PetscOptionsEnd();
    
  // Create nonlinear solver context
  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ( ierr );

  // Create matrix and vector data structures; set corresponding routines
  // Create vectors for solution and nonlinear function
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
	
  // Sets the local and global sizes, and checks to determine compatibility
  // VecSetSizes( the vector , the local size , the global size )
  ierr = VecSetSizes(x,PETSC_DECIDE,num_eq);CHKERRQ(ierr);
	
  // Configures the vector from the options database.
  // must be used after VecCreate but before the vector is used
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
	
  // Creates a new vector of the same type as an existing vector.
  ierr = VecDuplicate(x,&r);CHKERRQ(ierr);
	
  // Create Jacobian matrix data structure
  ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
  ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,num_eq,num_eq);CHKERRQ(ierr);
  ierr = MatSetFromOptions(J);CHKERRQ(ierr);
  ierr = MatSetUp(J);CHKERRQ(ierr);

  // Set function evaluation routine and vector.
  ierr = SNESSetFunction(snes, r, ResidualFunction, static_cast<void*>(&T.user_data) );CHKERRQ(ierr);

//  // Set Jacobian matrix data structure and Jacobian evaluation routine
//  if (this->FD_jacobian == PETSC_FALSE && this->ptype==1 ){
    ierr = SNESSetJacobian(snes , J , J , FormJacobian , static_cast<void*>(&T.user_data) );CHKERRQ(ierr);    
//  }
//  else if (this->FD_jacobian == PETSC_TRUE && this->ptype==0 ){
//    ierr = SNESSetJacobian(snes , J , J , SNESDefaultComputeJacobian , PETSC_NULL );CHKERRQ(ierr);    
//  } 
//  else if ( this->FD_jacobian==PETSC_TRUE && this->ptype==1 ){
//    throw MAP_ERROR_66;
//  }
//  else {
//    throw MAP_ERROR_67;
//  };
	
  // Customize nonlinear solver; set runtime options 
  // Set linear solver defaults for this problem. By extracting the
  // KSP, KSP, and PC contexts from the SNES context, we can then
  // directly call any KSP, KSP, and PC routines to set various options.
  ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr); 
	
  // Set SNES/KSP/KSP/PC runtime options, e.g.,
  //   -- snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
  // 
  // These options will override those specified above as long as
  // SNESSetFromOptions() is called _after_ any other customization
  // routines.	
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

  // Pass initial guess into vector x
  for ( int i=0 ; i<num_eq ; i++ ) {
    PETSc_number = boost::lexical_cast<PetscScalar> ( T.user_data.getConstraint(i) );
    VecSetValues( x, 1, &i, &PETSc_number, INSERT_VALUES );
  };//END for

  // set the non-zero structure of the Jacobian
  //
  //  J = [  A     B ]
  //      [ -B^T   C ]
  //  
  // We are only setting the non-zero structure for the A and B blocks. The B block handles the 
  // coupling between the force-balance equations and catenary equations
  T.user_data.initializeJacobian();

  return 0;
};


/**
 * ====================================================================================================
 * PetscSolve
 *
 *  @note : The user should initialize the vector, x, with the initial guess
 *          for the nonlinear solver prior to calling SNESSolve().  In particular,
 *          to employ an initial guess of zero, the user should explicitly set
 *           this vector to zero by calling VecSet().
 * ====================================================================================================
 */
int Numerics::PetscSolve( MAP_OtherStateType_class &other ,
                          MAP_ErrStat              &Error , 
                          MAP_Message              &Msg )
{
  ierr = SNESSolve( snes, PETSC_NULL, x ); CHKERRQ(ierr);
  ierr = SNESGetConvergedReason(snes , &reason ); CHKERRQ(ierr);
  
  if( other.CheckResidualConvergence( Error, Msg ) != 0 ){
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_69 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_69 );
    Error.set_error_key( MAP_ERROR );

    // Here is the exception message 
    Msg.RecordErrorToErrorList( str );

    str.erase(2,1);
    std::ostringstream S;
    S.str("");S.clear();
    S << ">>>> " +str +"    :    MAP error in file ";
    S << __FILE__;
    S << " on line ";
    S << __LINE__;
    S << " occurred.";
    Msg.WriteErrorToOutputFile( S.str() );

    this->ConvergeReason( Error, Msg );
  }  
  
  return 0;
};


/**
 * ====================================================================================================
 * FormJacobian
 * ====================================================================================================
 */
#undef __FUNCT__
#define __FUNCT__ "FormJacobian"
PetscErrorCode FormJacobian( SNES         snes_in , 
                             Vec          x_in    , 
                             Mat          *jac    ,
			     Mat          *b      , 
                             MatStructure *flag   , 
                             void         *ctx ) 
{
  UserData    *data = static_cast<UserData*>(ctx);    
  const int    M = data->getNumNodeEqs();    
  PetscScalar *xx = NULL;
  PetscScalar  D[4];
  PetscInt     idAx[2] = {0,1};    

  //MatView( *jac , PETSC_VIEWER_STDOUT_SELF );
  MatZeroEntries( *jac );

  // Get pointer to vector data
  VecGetArray( x_in , &xx );    

  for ( int i=0 ; i<data->getNumJacAEntries() ; i++ ) {
    MatSetValue( *jac , data->Ai(i) ,   data->Aj(i) ,  data->getJacA(i) , INSERT_VALUES );
  };//END for

    // form jacobian for the catenary equations with respect to H, V, and Lu
  for( int i=0 ; i<data->sizeOfElement() ; i++ ){
    D[0] = data->getdXdH( i );
    D[1] = data->getdXdV( i );
    D[2] = data->getdZdH( i );
    D[3] = data->getdZdV( i );

    idAx[0] = M+2*i;
    idAx[1] = M+2*i+1;
        
    // set D block diagonals
    MatSetValues( *jac , 2 , idAx , 2 , idAx , D , INSERT_VALUES );
  };

  // now fill the B (off-diagonal) block matrix and it's transpose
  for( int i=0 ; i<data->getNumJacBEntries() ; i++ ){
    MatSetValue( *jac , M+data->Bi(i) ,   data->Bj(i) ,  data->getJacB(i) , INSERT_VALUES );
    MatSetValue( *jac ,   data->Bj(i) , M+data->Bi(i) , -data->getJacB(i) , INSERT_VALUES );
  };

  *flag = SAME_NONZERO_PATTERN;
    
  // Restore vector
  VecRestoreArray( x_in , &xx );

  // Assemble matrix
  MatAssemblyBegin( *jac , MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd  ( *jac , MAT_FINAL_ASSEMBLY );
 
//    MatView( *jac , PETSC_VIEWER_STDOUT_SELF );
//    std::cin.get();
    
  return(0);
};


/**
 * ====================================================================================================
 * ResidualFunction
 * ====================================================================================================
 */
#undef __FUNCT__
#define __FUNCT__ "ResidualFunction"
PetscErrorCode ResidualFunction( SNES snes, Vec x, Vec f, void *ctx ){
    PetscErrorCode    ierr;
    const PetscScalar *xx;
    PetscScalar       *ff;
    MPI_Comm          comm;
    PetscMPIInt       size,rank;

    UserData *data = static_cast<UserData*>(ctx);
    
    ierr = PetscObjectGetComm((PetscObject)snes,&comm);CHKERRQ(ierr);
    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

    /**
     * Get pointers to vector data.
     *   -- For default PETSc vectors, VecGetArray() returns a pointer to
     *      the data array.  Otherwise, the routine is implementation dependent.
     *   -- You MUST call VecRestoreArray() when you no longer need access to
     *      the array.
     */
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);
    
    ff = data->FunctionEvaluations( ff, xx );
    
    // Restore vectors
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr); 
    
    return 0;
};


/**
 * ====================================================================================================
 * FunctionEvaluations
 * ====================================================================================================
 */
double *UserData::FunctionEvaluations( PetscScalar *FF , const PetscScalar *XX ){
  // The first step is to copy the constraint varaibles into XX (the state vector)
  for( int i=0 ; i<this->sizeOfConstraint() ; i++ ){
    this->setConstraint( i , XX[i] );
  }//END for

  // set sum_FX, sum_FY and sum_FZ to zero
  for ( int i=0 ; i<this->sizeOfNode( ) ; i++ ){
    this->node[i]->setSumForceToZero( );
  };

  // set fairlead and anchor nodes to zero. In other words, we are initializing the model. This calls
  // Node::fairlead->setSumForceToZero();
  // Node::anchor->setSumForceToZero();
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    this->element[i]->cleanNodes();
  };

  // update psi, l and h for each element being iterated
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    this->element[i]->updateElement( *error_status , *message );
  };

  int cnt = 0;
  double K = 1.0; 
  K = 1.0;//this->element[0]->getEA();

  for ( int i=0 ; i<this->sizeOfNode() ; i++ ) {	
    // solve X direction Newton equation
    if (this->node[i]->getXNewtonEquationFlag()==true){
      FF[cnt] = (1/K)*this->node[i]->f_x();
      cnt++;
    };// END if
	
    // solve Y direction Newton equation
    if (this->node[i]->getYNewtonEquationFlag()==true){
      FF[cnt] = (1/K)*this->node[i]->f_y();
      cnt++;
    };// END if

    // solve Z direction Newton equation
    if (this->node[i]->getZNewtonEquationFlag()==true){
      FF[cnt] = (1/K)*this->node[i]->f_z();
      cnt++;
    }
  };//END for

    // Compute function :
    // For each element, the X,Y,Z catenary equation is solved.  
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    FF[cnt] = this->element[i]->f_h(); // horizontal catenary eq. 
    cnt++;
    FF[cnt] = this->element[i]->f_v(); // horizontal catenary eq. 
    cnt++;
  };//END for
  
    // Make sure the number of equations is equal to the 
    // number of unknowns 
    //
    // @todo : restore the following line
  assert( this->sizeOfConstraint() == cnt ); 

  return FF;
};


/**
 * ====================================================================================================
 * ConvergeReason
 *
 * Converged:
 *   -- SNES_CONVERGED_FNORM_ABS          :  2, ||F|| < atol 
 *   -- SNES_CONVERGED_FNORM_RELATIVE     :  3, ||F|| < rtol*||F_initial|| 
 *   -- SNES_CONVERGED_SNORM_RELATIVE     :  4, Newton computed step size small; || delta x || < stol 
 *   -- SNES_CONVERGED_ITS                :  5, maximum iterations reached 
 *   -- SNES_CONVERGED_TR_DELTA           :  7,
 * 
 * Diverged
 *   -- SNES_DIVERGED_FUNCTION_DOMAIN     : -1, the new x location passed the function is not in the 
 *                                              domain of F
 *   -- SNES_DIVERGED_FUNCTION_COUNT      : -2,
 *   -- SNES_DIVERGED_LINEAR_SOLVE        : -3, the linear solve failed 
 *   -- SNES_DIVERGED_FNORM_NAN           : -4,
 *   -- SNES_DIVERGED_MAX_IT              : -5,
 *   -- SNES_DIVERGED_LINE_SEARCH         : -6, the line search failed 
 *   -- SNES_DIVERGED_INNER               : -7, inner solve failed 
 *   -- SNES_DIVERGED_LOCAL_MIN           : -8, || J^T b || is small, implies converged to local 
 *                                              minimum of F()
 *   -- SNES_CONVERGED_ITERATING          :  0 SNESConvergedReason;
 * ====================================================================================================
 */
void Numerics::ConvergeReason( MAP_ErrStat &Error, MAP_Message &Msg ){
  try {
    switch( reason ){
    /**
     * =======  PETSc Convergence codes  ======     <----------------------------------------------------------+
     */                                                                                            //          |
    case 0 :                                                                                       //          |
      Msg.WriteErrorToOutputFile("Converged (PETSc code 0).");                                                //          |
      break;                                                                                       //          |
    case 2 :                                                                                       //          |
      Msg.WriteErrorToOutputFile("Converged (PETSc code 2: '||F|| < atol ').");                               //          |
      break;                                                                                       //          |
    case 3 :                                                                                       //          |
      Msg.WriteErrorToOutputFile("Converged (PETSc code 3: '||F|| < rtol*||F_initial|| ').");                 //          |
      break;                                                                                       //          |
    case 4 :                                                                                       //          |
      Msg.WriteErrorToOutputFile("Converged (PETSc code 4: 'Step size small; || delta x || < stol ').");      //          |
      break;                                                                                       //          |
    case 5 :                                                                                       //          |
      Msg.WriteErrorToOutputFile("Converged (PETSc code 5: 'Maximum iteration reached').");                   //          |
      break;                                                                                       //          |
    case 7 :                                                                                       //          |
      Msg.WriteErrorToOutputFile("Converged (PETSc code 7).");                                                //          |
      break;                                                                                       //  --------+
      //============== <END> ===================================================================================
                
      /**
       * =======  PETSc Convergence codes  ======     <--------------------------------------------------------+
       */                                                                                          //          |
    case -1 :                                                                                      //          |
      throw MAP_ERROR_57;                                                                          //          |
      break;                                                                                       //          |
    case -2 :                                                                                      //          |
      throw MAP_ERROR_58;                                                                          //          |
      break;                                                                                       //          |
    case -3 :                                                                                      //          |
      throw MAP_ERROR_59;                                                                          //          |
      break;                                                                                       //          |
    case -4 :                                                                                      //          |
      throw MAP_ERROR_60;                                                                          //          |
      break;                                                                                       //          |
    case -5 :                                                                                      //          |
      throw MAP_ERROR_61;                                                                          //          |
      break;                                                                                       //          |
    case -6 :                                                                                      //          |
      throw MAP_ERROR_62;                                                                          //          |
      break;                                                                                       //          |
    case -7 :                                                                                      //          |
      // I am aware this case is contrived and doesn't really exist in practice                    //          |
      throw MAP_ERROR_63;                                                                          //          |
      break;                                                                                       //          |
    case -8 :                                                                                      //          |
      throw MAP_ERROR_64;                                                                          //          |
      break;                                                                                       //  --------+
      //============== <END> ===================================================================================
    };// END switch
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code );
    Error.set_error_key( MAP_ERROR );

    // Here is the exception message 
    Msg.RecordErrorToErrorList( str );

    str.erase(2,1);
    std::ostringstream S;
    S.str("");S.clear();
    S << ">>>> " +str +"    :    MAP error in file ";
    S << __FILE__;
    S << " on line ";
    S << __LINE__;
    S << " occurred.";
    Msg.WriteErrorToOutputFile( S.str() );
  };// END try
};


/**
 * ====================================================================================================
 * PetscEnd
 * ====================================================================================================
 */
int Numerics::End( MAP_ErrStat &Error, MAP_Message &Msg ){
    ierr = SNESGetIterationNumber(snes,&its);CHKERRQ(ierr);

    //    VecView(r,PETSC_VIEWER_STDOUT_WORLD);    
    //    VecView(x,PETSC_VIEWER_STDOUT_WORLD);    
    //    //    std::cout << "\n";
    //    MatView(J,PETSC_VIEWER_STDOUT_WORLD);
    //
    //    MatStructure  flag;
    //    SNESComputeJacobian(snes,x,&J,&J,&flag);
    //    MatView(J,PETSC_VIEWER_STDOUT_WORLD);

    // Free work space.  PETSc variables/objects are destroyed when the program is done
    //ierr = MatFDColoringDestroy(&matfdcoloring);CHKERRQ(ierr);
    ierr = VecDestroy(&x);CHKERRQ(ierr); 
    ierr = VecDestroy(&r);CHKERRQ(ierr);
    //ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
    //ierr = PCDestroy(&pc);CHKERRQ(ierr);
    ierr = MatDestroy(&J);CHKERRQ(ierr);  /* @todo : insert a check to see if J is allocated. Then delete it*/
    ierr = SNESDestroy(&snes);CHKERRQ(ierr);

    ierr = PetscFinalize();

    return 0;
};
