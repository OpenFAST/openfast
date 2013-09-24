/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        Numerics.h
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


#ifndef _NUMERICS_H
#define _NUMERICS_H


class MAP_OtherStateType_class;


/**
 * This is the guts of the numerics solver. It is dependent on the PETSc library. Current, MPI is not 
 * supported, and the PETSc uniprocessor case functions just fine for now (because the size of the 
 * unknowns is small).
 *
 * @property   SNES                      snes                  nonlinear solver context 
 * @property   KSP                       ksp;                  Krylov sup-space context
 * @property   PC                        pc;                   preconditioner context
 * @property   Vec                       x;                    solution vector
 * @property   Vec                       r;                    residual vectors
 * @property   Mat                       J;                    Jacobian matrix
 * @property   PetscErrorCode            ierr;                 PETSc error context
 * @property   PetscInt                  its;                  number of iterations to convergence
 * @property   PetscMPIInt               size;                 Not needed (for MPI)
 * @property   PetscMPIInit              rank;                 Not needed (for MPI)
 * @property   PetscScalar               *X;                         
 * @property   PetscBool                 help_flag;            Help flag. If '-help' is in the input file, 
 *                                                             this is true
 * @property   SNESConvergedReason       reason;               Reason the non-inear solver converged
 * @property   PetscBool                 msqs_fd_jacobian;     is the finite-differenced Jacobian used?
 * @property   PetscBool                 msqs_default_setting  are we using the default solver options?
 * @property   PetscReal                 msqs_tol;             set the msqs residual function tolerance
 * @property   std::vector<std::string>  options_string;       PETSc and MAP(MSQS) run-time options
 */
class 
Numerics 
{
private:
  SNES                     snes;         // nonlinear solver context 
  KSP                      ksp;          // linear solver context
  PC                       pc;           // preconditioner context
  Vec                      x;            // solution
  Vec                      r;            // residual vectors
  Mat                      J;            // Jacobian matrix
  PetscErrorCode           ierr;
  PetscInt                 its;
  PetscMPIInt              size;
  PetscMPIInt              rank;
  PetscScalar              *X;
  SNESConvergedReason      reason;

  // option built into MAP
  PetscBool                help_flag;
  PetscReal                msqs_tol;
  PetscBool                msqs_k;  

  // PETSc run-time options listed in the MAP input file
  std::vector<std::string> options_string;  

public:
  
Numerics( ) : 
  X           ( NULL       ) ,   // make sure X points to nothing
  help_flag   ( PETSC_TRUE ) ,   // default for the help flag. This is false if '-help' is defined in input 
  msqs_tol    (1e-2        ) ,
  msqs_k      ( PETSC_TRUE ) {}  // maximum value of all residuals to guarantee convergence
  
  ~Numerics( ) { /*ierr = PetscFinalize()*/ }
        
  int  InitializeSolver( MAP_OtherStateType_class &other ,
                         MAP_InitInputType_class  &init  ,
                         MAP_ErrStat_class              &err   ,  
                         MAP_Message_class              &msg   );
  void setNumericsOptionsString( const std::string &optionStr );

  int  PetscSolve           ( MAP_OtherStateType_class &other ,
                              MAP_ErrStat_class        &Error , 
                              MAP_Message_class        &Msg );
  void PetscConvergeReason ( MAP_ErrStat_class &Error, MAP_Message_class &Msg );
  int  PetscEnd       ( MAP_ErrStat_class &err, MAP_Message_class &msg );    

  PetscBool GetHelpFlag ( ) const { return help_flag;      }
  PetscReal GetMSQSTol  ( ) const { return this->msqs_tol; }
  PetscBool GetMsqsKFlag( ) const { return this->msqs_k;   } 
};

#endif // _NUMERICS_H
