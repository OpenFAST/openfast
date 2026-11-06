/**********************************************************************************************************************************
 * LICENSING
 * Copyright (C) 2026 National Renewable Energy Laboratory
 *
 * This file is part of BeamDyn.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 **********************************************************************************************************************************/
/*
 * C interface to a standalone BeamDyn beam (see BeamDyn_C_Binding.f90).
 *
 * All routines follow the Fortran BIND(C) conventions used by the OpenFAST module C bindings: every argument is
 * passed by reference, strings are C_NULL_CHAR terminated char arrays, and matrices are flat arrays.  Direction
 * cosine matrices (DCMs) follow the OpenFAST convention (they map a vector from the global frame to the local
 * frame) and are stored row by row: dcm[3*i + j] is row i, column j.  Node arrays are stored node by node:
 * pos[3*n + k] is component k of node n.
 *
 * Precision: root kinematics (position, orientation, velocity, acceleration) and gravity are double; loads,
 * node motions, reaction loads, and channel values are float; orientation DCMs are always double.
 *
 * Error handling: each routine returns ErrStat_C (0 none, 1 info, 2 warning, 3 severe, 4 fatal) and a
 * C_NULL_CHAR terminated message in ErrMsg_C, which must be at least BD_C_ERRMSGLEN characters long.  After a
 * fatal error the module is left in a state where BD_C_End can still be called safely.
 *
 * Call sequence:
 *    BD_C_Init                 initialize the beam from an input file (returns node counts and output channel info)
 *    BD_C_GetRefPositions      (optional) reference positions of the output nodes and load nodes
 *    loop over time:
 *       BD_C_SetRootMotion     root motion at the time of the next BD_C_CalcOutput or BD_C_UpdateStates call
 *       BD_C_SetPointLoads     (optional) point loads at the finite element nodes
 *       BD_C_SetDistrLoads     (optional) distributed loads at the quadrature point nodes
 *       BD_C_CalcOutput        outputs at the current time (node motions, root reaction, output channels)
 *       BD_C_UpdateStates      advance the states from Time_C to TimeNext_C using the inputs set for TimeNext_C
 *    BD_C_PackStates / BD_C_UnpackStates   (optional) write/read a checkpoint file with the complete beam state
 *    BD_C_End
 */
#ifndef BEAMDYN_C_BINDING_H
#define BEAMDYN_C_BINDING_H

#ifdef __cplusplus
extern "C" {
#endif

/* Sizes used by the interface (must match the NWTC Library and BeamDyn_C_Binding.f90) */
#define BD_C_ERRMSGLEN      8197      /* ErrMsgLen + 1: minimum length of the ErrMsg_C buffer */
#define BD_C_INTFSTRLEN     1025      /* Maximum length of file name strings passed through the interface */
#define BD_C_CHANLEN        20        /* Length of each output channel name and unit string */
#define BD_C_MAXOUTPUTS     8000      /* Maximum number of output channels the interface can return */
#define BD_C_CHANNELBUFLEN  (BD_C_CHANLEN*BD_C_MAXOUTPUTS+1)   /* Minimum length of the channel name and unit buffers */

/* Error status values (ErrStat_C) */
#define BD_C_ERRID_NONE     0
#define BD_C_ERRID_INFO     1
#define BD_C_ERRID_WARN     2
#define BD_C_ERRID_SEVERE   3
#define BD_C_ERRID_FATAL    4

/*
 * Initialize a BeamDyn beam from a BeamDyn primary input file (the blade file it names is read from the file
 * system relative to the primary input file).  The root of the beam is placed at RootPos_C with orientation
 * RootOri_C, which is also the reference frame for the BeamDyn calculations.  On return, the node counts give the
 * sizes of the arrays used by the other routines.
 */
void BD_C_Init(
   const int    *InputFilePassed,         /* IN:  1: InputFileString_C holds the contents of the input file (lines separated by '\0');
                                                   0: InputFileString_C holds the path to the input file */
   const char  **InputFileString_C,       /* IN:  pointer to the input file contents or path string (see InputFilePassed) */
   const int    *InputFileStringLength_C, /* IN:  length of InputFileString_C */
   const char   *OutRootName_C,           /* IN:  root name for the summary (<root>.BD.sum), echo (<root>.BD.ech), and output (<root>.out) files */
   const double *RootPos_C,               /* IN:  [3] initial root position in the global frame (m) */
   const double *RootOri_C,               /* IN:  [9] initial root orientation DCM, global to root frame, row by row */
   const double *RootVel_C,               /* IN:  [6] initial root translational (0:2) and rotational (3:5) velocity, global frame (m/s, rad/s) */
   const double *Gravity_C,               /* IN:  [3] gravitational acceleration vector, global frame (m/s^2) */
   const double *DT_C,                    /* IN:  time step for BD_C_UpdateStates (s); the input file DTBeam must be DEFAULT or equal to this */
   const int    *InterpOrder_C,           /* IN:  input interpolation/extrapolation order: 1 (linear) or 2 (quadratic) */
   const int    *DynamicSolve_C,          /* IN:  1: dynamic solve; 0: static solve */
   const int    *WrOutputs_C,             /* IN:  1: write the output channels to <root>.out at each BD_C_CalcOutput call; 0: no file */
   int          *NumOutputNodes_C,        /* OUT: number of nodes on the blade motion output mesh */
   int          *NumPointLoadNodes_C,     /* OUT: number of nodes on the point load input mesh (finite element nodes) */
   int          *NumDistrLoadNodes_C,     /* OUT: number of nodes on the distributed load input mesh (quadrature point nodes) */
   int          *NumChannels_C,           /* OUT: number of output channels */
   char         *OutputChannelNames_C,    /* OUT: [BD_C_CHANNELBUFLEN] channel names, BD_C_CHANLEN characters each, '\0' terminated */
   char         *OutputChannelUnits_C,    /* OUT: [BD_C_CHANNELBUFLEN] channel units, BD_C_CHANLEN characters each, '\0' terminated */
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

/*
 * Reference (undeflected) positions of the output nodes and load input nodes.  Arrays are sized by the counts
 * returned from BD_C_Init.
 */
void BD_C_GetRefPositions(
   float        *OutputNodePos_C,         /* OUT: [3*NumOutputNodes] reference positions of the output nodes (m) */
   double       *OutputNodeOri_C,         /* OUT: [9*NumOutputNodes] reference orientation DCMs of the output nodes, global to node frame, row by row */
   float        *PointLoadNodePos_C,      /* OUT: [3*NumPointLoadNodes] reference positions of the point load nodes (m) */
   float        *DistrLoadNodePos_C,      /* OUT: [3*NumDistrLoadNodes] reference positions of the distributed load nodes (m) */
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

/*
 * Set the root motion used by the next BD_C_CalcOutput (motion at the current time) or BD_C_UpdateStates (motion
 * at the next time) call.
 */
void BD_C_SetRootMotion(
   const double *RootDisp_C,              /* IN:  [3] root translational displacement from the initial root position, global frame (m) */
   const double *RootOri_C,               /* IN:  [9] root orientation DCM, global to root frame, row by row */
   const double *RootVel_C,               /* IN:  [6] root translational (0:2) and rotational (3:5) velocity, global frame (m/s, rad/s) */
   const double *RootAcc_C,               /* IN:  [6] root translational (0:2) and rotational (3:5) acceleration, global frame (m/s^2, rad/s^2) */
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

/*
 * Set the point loads at the finite element nodes (global frame), used by the next BD_C_CalcOutput or
 * BD_C_UpdateStates call.
 */
void BD_C_SetPointLoads(
   const float  *PointLoads_C,            /* IN:  [6*NumPointLoadNodes] force (0:2) and moment (3:5) at each point load node (N, N-m) */
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

/*
 * Set the distributed loads per unit length at the quadrature point nodes (global frame), used by the next
 * BD_C_CalcOutput or BD_C_UpdateStates call.
 */
void BD_C_SetDistrLoads(
   const float  *DistrLoads_C,            /* IN:  [6*NumDistrLoadNodes] force (0:2) and moment (3:5) per unit length at each distributed load node (N/m, N-m/m) */
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

/*
 * Advance the states from Time_C to TimeNext_C (= Time_C + DT_C).  The inputs set by the BD_C_Set* routines are
 * taken as the inputs at TimeNext_C.  Calling again with the same Time_C repeats the step as a correction step with
 * the new inputs.
 */
void BD_C_UpdateStates(
   const double *Time_C,                  /* IN:  current time (s) */
   const double *TimeNext_C,              /* IN:  time to advance the states to (s) */
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

/*
 * Compute the outputs at Time_C from the current states and the inputs set by the BD_C_Set* routines.
 */
void BD_C_CalcOutput(
   const double *Time_C,                  /* IN:  current time (s) */
   float        *NodePos_C,               /* OUT: [3*NumOutputNodes] position of each output node, global frame (m) */
   double       *NodeOri_C,               /* OUT: [9*NumOutputNodes] orientation DCM of each output node, global to node frame, row by row */
   float        *NodeVel_C,               /* OUT: [6*NumOutputNodes] translational (0:2) and rotational (3:5) velocity of each output node, global frame (m/s, rad/s) */
   float        *NodeAcc_C,               /* OUT: [6*NumOutputNodes] translational (0:2) and rotational (3:5) acceleration of each output node, global frame (m/s^2, rad/s^2) */
   float        *RootReaction_C,          /* OUT: [6] reaction force (0:2) and moment (3:5) at the root, global frame (N, N-m) */
   float        *OutputChannelValues_C,   /* OUT: [NumChannels] output channel values */
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

/*
 * Write a checkpoint file <CheckpointRoot>.chkp holding the complete state of the beam (states at the current and
 * previous time, input history, inputs set through the interface, and time stepping information).  The file is
 * only valid for the same build and the same input file.
 */
void BD_C_PackStates(
   const char   *CheckpointRoot_C,        /* IN:  root name of the checkpoint file (".chkp" is appended) */
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

/*
 * Restore the state of the beam from a checkpoint file written by BD_C_PackStates.  BD_C_Init must have been
 * called with the same input file first.
 */
void BD_C_UnpackStates(
   const char   *CheckpointRoot_C,        /* IN:  root name of the checkpoint file (".chkp" is appended) */
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

/*
 * Free all memory held by the library and close any open files.  Safe to call after a fatal error or before
 * BD_C_Init.
 */
void BD_C_End(
   int          *ErrStat_C,               /* OUT: error status */
   char         *ErrMsg_C                 /* OUT: [BD_C_ERRMSGLEN] error message */
);

#ifdef __cplusplus
}
#endif

#endif /* BEAMDYN_C_BINDING_H */
