/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   map[dot]plus[dot]plus[dot]help[at]gmail                    *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


#ifndef _NUMERIC_H
#define _NUMERIC_H

#include "map.h"


/**
 * @brief   Forms the LU matrix for inverting the Jacobian - step 1
 * @details Called by {@link node_solve_sequence} to form the 
 *          lower (L) and upper (U) triangular matrices to solve
 *          the system of equations. This is called in the outer
 *          (node) solve sequence. Back substitution through the
 *          function {@link lu_back_substitution} follows this 
 *          function call to complete the newton iteration.
 *          $\mathbf{A}=\mathbf{LU}$
 *          with 
 *          $\begin{bmatrix}
 *             a_{11} & a_{12} & a_{13} \\
 *             a_{21} & a_{22} & a_{23} \\
 *             a_{31} & a_{32} & a_{33} \\
 *           \end{bmatrix} =
 *           \begin{bmatrix}
 *             l_{11} & 0      & 0 \\
 *             l_{21} & l_{22} & 0 \\
 *             l_{31} & l_{32} & l_{33} \\
 *           \end{bmatrix}
 *           \begin{bmatrix}
 *             u_{11} & u_{12} & u_{13} \\
 *             0      & u_{22} & u_{23} \\
 *             0      & 0      & u_{33} \\
 *           \end{bmatrix}
 *           $
 * @param   ns, outer solve attributes; preserves l, jac, and u
 * @param   n, matrix size, nxn
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     node_solve_sequence()
 * @return  MAP error code
 */
MAP_ERROR_CODE root_finding_step(OuterSolveAttributes* ns, const int n, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, double* error, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Wrapper for cminpack's lmder function
 * @details Called by {@link solve_line} to solve the catenary
 *          equations. The minpack documentation describes the 
 *          the lmder function as: 
 *          <pre>
 *          The purpose of LMDER is to minimize the sum of the 
 *          squares of M nonlinear functions in N variables by 
 *          a modification of the Levenberg-Marquardt algorithm.  
 *          The user must provide a subroutine which calculates 
 *          the functions and the Jacobian.
 *          </pre>
 * @param   line, the line structure to be checked
 * @param   inner_opt,
 * @param   line_num, line number for error logging 
 * @param   time, current time
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     solve_line()
 * @return  MAP error code
 */
MAP_ERROR_CODE call_minpack_lmder(Line* line, InnerSolveAttributes* inner_opt, const int line_num, const float time, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Return the residual and Jacobian for each line
 * @details Passed as a function pointer to the cminpack lmder routine.
 *          Used in function {@link call_minpack_lmder} to solve the catenary
 *          equations. 
 * @param   line_ptr, the line structure
 * @param   m,
 * @param   n,
 * @param   x, current state vector
 * @param   fvec, residual vector
 * @param   fjac, jacobian matrix
 * @param   ldfjac, 
 * @param   iflag, solve status; 0 = error
 * @return  minpack error code
 */
int inner_function_evals(void* line_ptr, int m, int n, const __cminpack_real__* x, __cminpack_real__* fvec, __cminpack_real__* fjac, int ldfjac, int iflag);


/**
 * @brief   Forms the LU matrix for inverting the Jacobian - step 1
 * @details Called by {@link node_solve_sequence} to form the 
 *          lower (L) and upper (U) triangular matrices to solve
 *          the system of equations. This is called in the outer
 *          (node) solve sequence. Back substitution through the
 *          function {@link lu_back_substitution} follows this 
 *          function call to complete the newton iteration.
 *          $\mathbf{A}=\mathbf{LU}$
 *          with 
 *          $\begin{bmatrix}
 *             a_{11} & a_{12} & a_{13} \\
 *             a_{21} & a_{22} & a_{23} \\
 *             a_{31} & a_{32} & a_{33} \\
 *           \end{bmatrix} =
 *           \begin{bmatrix}
 *             l_{11} & 0      & 0 \\
 *             l_{21} & l_{22} & 0 \\
 *             l_{31} & l_{32} & l_{33} \\
 *           \end{bmatrix}
 *           \begin{bmatrix}
 *             u_{11} & u_{12} & u_{13} \\
 *             0      & u_{22} & u_{23} \\
 *             0      & 0      & u_{33} \\
 *           \end{bmatrix}
 *           $
 * @param   ns, outer solve attributes; preserves l, jac, and u
 * @param   n, matrix size, nxn
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     node_solve_sequence()
 * @return  MAP error code
 */
MAP_ERROR_CODE lu(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Solves $\mathbf{Ax}=\mathbf{b}$ through forward and backward substitution
 * @details Called by {@link node_solve_sequence} to solve: 
 *          $\mathbf{Ly}=\mathbf{b}$
 *          then:
 *          $\mathbf{Ux}=\mathbf{y}$
 *          The vector \mathbf{x} is used to increment the Newton step.
 *          Note that $\mathbf{x} = \mathbf{J}^{-1} (\textup{residual})$ with
 *          $\mathbf{x}=$ns->x
 * @param   ns, outer solve attributes; preserves l, u, x, and y
 * @param   n, matrix size, nxn
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     node_solve_sequence()
 * @return  MAP error code
 */
MAP_ERROR_CODE lu_back_substitution(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr);


#endif // _NUMERIC_H
