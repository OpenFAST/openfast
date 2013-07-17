/**
 * ====================================================================================================
 *                              Numerics.h
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


#ifndef _NUMERICS_H
#define _NUMERICS_H


class MAP_OtherStateType_class;


/**
 * ====================================================================================================
 * Numerics
 *
 * 
 * ====================================================================================================
 */
class Numerics {
private:
  SNES                     snes;         // nonlinear solver context 
  KSP                      ksp;          // linear solver context
  PC                       pc;           // preconditioner context
  Vec                      x,r;          // solution, residual vectors
  Mat                      J;            // Jacobian matrix
  PetscErrorCode           ierr;
  PetscInt                 its;
  PetscMPIInt              size,rank;
  PetscScalar              pfive, *X;
  PetscBool                flg;
  SNESConvergedReason      reason;
  PetscBool                FD_jacobian;
 
//  // DM objects for field splitting
//  PetscInt ptype;
//  DM             da_elem , da_node , pack;
//  const PetscInt *lx_elem;
//  PetscInt       *lx_node, m, nprocs;
//  Vec            x_elem , X_node , F_node , F_elem;
//  PetscBool      view_draw , pass_dm;

  // PETSc run-time options listed in the MAP input file
  std::vector<std::string> options_string;  

public:
    
Numerics( ) : X ( NULL ) { }
  ~Numerics( ) { }
        
  int  InitializeSolver( MAP_OtherStateType_class &T     , 
                         MAP_InitInputType_class  &Init  , 
                         MAP_ErrStat              &Error , 
                         MAP_Message              &Msg );

  void setNumericsOptionsString( const std::string &P );

  int  PetscSolve           ( MAP_OtherStateType_class &other ,
                              MAP_ErrStat        &Error , 
                              MAP_Message        &Msg );
  void ConvergeReason       ( MAP_ErrStat &Error, MAP_Message &Msg );
  int  End                  ( MAP_ErrStat &Error, MAP_Message &Msg );
    
  PetscBool getFlg( ) { return flg; }
};

#endif // _NUMERICS_H
