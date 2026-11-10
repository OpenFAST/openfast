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
 * Driver for the BeamDyn C interface.  It drives a standalone BeamDyn beam from C in the same way the BeamDyn
 * Fortran driver does (root rotating about the origin at a constant angular velocity, constant tip load, output file
 * with the requested channels) and provides checks of the interface against analytic results and against itself.
 *
 * Usage:  beamdyn_c_binding_driver <BeamDyn driver input file>
 *         beamdyn_c_binding_driver <mode> <BeamDyn primary input file> [key=value ...]
 *
 * With a single argument, the file is read as a BeamDyn driver input file (the input file of beamdyn_driver) and the
 * case is run through the C interface; <root>.out is written in the same format as beamdyn_driver writes it.  The
 * multi-point load table and GlbRotBladeT0 = False are not supported by this driver.
 *
 * Modes:
 *    run         Run the beam from t=0 to tmax and write <root>.out (same format as the BeamDyn driver output file).
 *    cantilever  Static tip deflection and first bending frequency of a uniform cantilever, compared with the
 *                Euler-Bernoulli values for the EI and mass per unit length given (EI=, mu=, tip=Fx).
 *    checkpoint  Run to tmax, then repeat the run with a checkpoint written at tchk and restored into a fresh
 *                instance; compare the tip motion and root reaction of the two runs.
 *
 * Options (key=value; vectors as comma separated values):
 *    root=NAME        root name for output files (default: run, cantilever, or checkpoint)
 *    dt=              time step (s)                                   [0.002]
 *    tmax=            end time (s)                                    [1.0]
 *    interp=          input interpolation order, 1 or 2               [1]
 *    dynamic=         1 dynamic solve, 0 static solve                 [1]
 *    wrout=           1 write <root>.out, 0 do not                    [1 for run, 0 otherwise]
 *    grav=gx,gy,gz    gravity vector (m/s^2)                          [0,0,0]
 *    pos=x,y,z        root position (m)                               [0,0,0]
 *    ori=r11,...,r33  root orientation DCM, global to root, row by row [identity]
 *    omega=wx,wy,wz   root angular velocity about the origin (rad/s)  [0,0,0]
 *    tip=Fx,Fy,Fz     point force at the tip node (N)                 [0,0,0]
 *    trelease=        time at which the tip force is removed (s)      [never]
 *    amp=  freq=      sinusoidal root displacement along x: amp*sin(2*pi*freq*t)  [0, 0]
 *    EI=  mu=         bending stiffness (N-m^2) and mass per unit length (kg/m) for the analytic comparisons
 *    tchk=            checkpoint time for the checkpoint mode (s)     [tmax/2]
 *    tol=             relative tolerance for pass/fail of the checks  [0.01]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "BeamDyn_C_Binding.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*----------------------------------------------------------------------------------------------------------------*/
/* Options                                                                                                          */
/*----------------------------------------------------------------------------------------------------------------*/
typedef struct {
   const char *input_file;
   char        root[BD_C_INTFSTRLEN];
   double      dt;
   double      tmax;
   int         interp;
   int         dynamic;
   int         wrout;
   double      grav[3];
   double      pos[3];
   double      ori[9];
   double      omega[3];
   double      tip[3];
   double      tipm[3];          /* point moment at the tip node (N-m) */
   double      distr[6];         /* constant distributed force and moment on every distributed load node (N/m, N-m/m) */
   int         use_distr;        /* set the distributed loads through the interface */
   double      t_release;
   double      amp;
   double      freq;
   double      EI;
   double      mu;
   double      t_chk;
   double      tol;
} Options;

static void set_defaults(Options *o, const char *mode)
{
   memset(o, 0, sizeof(*o));
   strncpy(o->root, mode, sizeof(o->root) - 1);
   o->dt        = 0.002;
   o->tmax      = 1.0;
   o->interp    = 1;
   o->dynamic   = 1;
   o->wrout     = (strcmp(mode, "run") == 0) ? 1 : 0;
   o->ori[0]    = 1.0;  o->ori[4] = 1.0;  o->ori[8] = 1.0;
   o->t_release = -1.0;
   o->t_chk     = -1.0;
   o->tol       = 0.01;
}

static int parse_vector(const char *s, double *v, int n)
{
   int   i;
   char *end;
   for (i = 0; i < n; i++) {
      v[i] = strtod(s, &end);
      if (end == s) return 1;
      s = end;
      if (*s == ',') s++;
   }
   return 0;
}

static int parse_options(int argc, char **argv, Options *o)
{
   int i, k;
   double v[9];
   for (i = 3; i < argc; i++) {
      const char *eq = strchr(argv[i], '=');
      if (eq == NULL) { fprintf(stderr, "Option '%s' is not of the form key=value\n", argv[i]); return 1; }
      size_t klen = (size_t)(eq - argv[i]);
      const char *val = eq + 1;
      int bad = 0;
#define KEY(name) (klen == strlen(name) && strncmp(argv[i], name, klen) == 0)
      if      (KEY("root"))     { strncpy(o->root, val, sizeof(o->root) - 1); }
      else if (KEY("dt"))       { o->dt = atof(val); }
      else if (KEY("tmax"))     { o->tmax = atof(val); }
      else if (KEY("interp"))   { o->interp = atoi(val); }
      else if (KEY("dynamic"))  { o->dynamic = atoi(val); }
      else if (KEY("wrout"))    { o->wrout = atoi(val); }
      else if (KEY("grav"))     { bad = parse_vector(val, v, 3); for (k = 0; k < 3; k++) o->grav[k] = v[k]; }
      else if (KEY("pos"))      { bad = parse_vector(val, v, 3); for (k = 0; k < 3; k++) o->pos[k]  = v[k]; }
      else if (KEY("ori"))      { bad = parse_vector(val, o->ori, 9); }
      else if (KEY("omega"))    { bad = parse_vector(val, o->omega, 3); }
      else if (KEY("tip"))      { bad = parse_vector(val, o->tip, 3); }
      else if (KEY("trelease")) { o->t_release = atof(val); }
      else if (KEY("amp"))      { o->amp = atof(val); }
      else if (KEY("freq"))     { o->freq = atof(val); }
      else if (KEY("EI"))       { o->EI = atof(val); }
      else if (KEY("mu"))       { o->mu = atof(val); }
      else if (KEY("tchk"))     { o->t_chk = atof(val); }
      else if (KEY("tol"))      { o->tol = atof(val); }
      else { fprintf(stderr, "Unknown option '%s'\n", argv[i]); return 1; }
#undef KEY
      if (bad) { fprintf(stderr, "Could not parse the value of option '%s'\n", argv[i]); return 1; }
   }
   return 0;
}

/*----------------------------------------------------------------------------------------------------------------*/
/* Error handling                                                                                                   */
/*----------------------------------------------------------------------------------------------------------------*/
static int  err_stat = 0;
static char err_msg[BD_C_ERRMSGLEN];

/* Print any message from the library; on a fatal error end the library and exit. */
static void check(const char *where)
{
   if (err_stat != BD_C_ERRID_NONE)
      printf("%s: ErrStat=%d\n%s\n", where, err_stat, err_msg);
   if (err_stat >= BD_C_ERRID_FATAL) {
      int  es;
      char em[BD_C_ERRMSGLEN];
      BD_C_End(&es, em);
      exit(EXIT_FAILURE);
   }
}

/*----------------------------------------------------------------------------------------------------------------*/
/* Root kinematics: rigid rotation about the origin at constant angular velocity (as in the BeamDyn driver), plus   */
/* an optional sinusoidal displacement along x.                                                                     */
/*----------------------------------------------------------------------------------------------------------------*/
static void cross(const double a[3], const double b[3], double c[3])
{
   c[0] = a[1]*b[2] - a[2]*b[1];
   c[1] = a[2]*b[0] - a[0]*b[2];
   c[2] = a[0]*b[1] - a[1]*b[0];
}

/* Rotation matrix (active rotation of a vector) for a rotation of angle about the unit vector k:  R = I + sin(a) K + (1-cos(a)) K^2 */
static void rodrigues(const double k[3], double angle, double R[9])
{
   double s = sin(angle), c = cos(angle), v = 1.0 - c;
   R[0] = c + k[0]*k[0]*v;        R[1] = k[0]*k[1]*v - k[2]*s;   R[2] = k[0]*k[2]*v + k[1]*s;
   R[3] = k[1]*k[0]*v + k[2]*s;   R[4] = c + k[1]*k[1]*v;        R[5] = k[1]*k[2]*v - k[0]*s;
   R[6] = k[2]*k[0]*v - k[1]*s;   R[7] = k[2]*k[1]*v + k[0]*s;   R[8] = c + k[2]*k[2]*v;
}

static void root_motion(const Options *o, double t, double disp[3], double ori[9], double vel[6], double acc[6])
{
   int    i, j, k;
   double w = sqrt(o->omega[0]*o->omega[0] + o->omega[1]*o->omega[1] + o->omega[2]*o->omega[2]);
   double R[9] = {1,0,0, 0,1,0, 0,0,1};
   double r0[3] = {o->pos[0], o->pos[1], o->pos[2]};
   double p[3], v[3], wxp[3], a[3];

   if (w > 0.0) {
      double khat[3] = {o->omega[0]/w, o->omega[1]/w, o->omega[2]/w};
      rodrigues(khat, w*t, R);
   }
   for (i = 0; i < 3; i++) p[i] = R[3*i]*r0[0] + R[3*i+1]*r0[1] + R[3*i+2]*r0[2];   /* rotated root position */
   cross(o->omega, p, v);                                                            /* velocity:  omega x p   */
   cross(o->omega, p, wxp);
   cross(o->omega, wxp, a);                                                          /* acceleration: omega x (omega x p) */

   for (i = 0; i < 3; i++) {
      disp[i]  = (p[i] - r0[i]);
      vel[i]   = v[i];
      vel[3+i] = o->omega[i];
      acc[i]   = a[i];
      acc[3+i] = 0.0;
   }
   /* The root frame rotates with the body: DCM(t) = DCM(0) * R^T */
   for (i = 0; i < 3; i++)
      for (j = 0; j < 3; j++) {
         ori[3*i+j] = 0.0;
         for (k = 0; k < 3; k++) ori[3*i+j] += o->ori[3*i+k] * R[3*j+k];
      }

   if (o->amp != 0.0) {
      double wf = 2.0*M_PI*o->freq;
      disp[0] += ( o->amp*sin(wf*t));
      vel[0]  += ( o->amp*wf*cos(wf*t));
      acc[0]  += (-o->amp*wf*wf*sin(wf*t));
   }
}

/*----------------------------------------------------------------------------------------------------------------*/
/* Time series recorded from the beam                                                                               */
/*----------------------------------------------------------------------------------------------------------------*/
typedef struct {
   int     n, nmax;
   double *t;
   double *tip;      /* [3*n] tip displacement from the reference tip position (m) */
   double *react;    /* [6*n] root reaction force and moment (N, N-m) */
   double  length;   /* distance from the root to the tip output node in the reference configuration (m) */
} Series;

static void series_alloc(Series *s, int nmax)
{
   s->n = 0; s->nmax = nmax;
   s->t     = (double*)calloc((size_t)nmax, sizeof(double));
   s->tip   = (double*)calloc((size_t)3*nmax, sizeof(double));
   s->react = (double*)calloc((size_t)6*nmax, sizeof(double));
   s->length = 0.0;
}

static void series_free(Series *s)
{
   free(s->t); free(s->tip); free(s->react);
   s->t = NULL; s->tip = NULL; s->react = NULL; s->n = 0;
}

/*----------------------------------------------------------------------------------------------------------------*/
/* Run the beam from t_start to t_end with the loads and root motion described by the options.                      */
/*   t_pack >= 0:      write a checkpoint when the simulation time reaches t_pack (and stop if stop_after_pack)      */
/*   unpack != 0:      restore the checkpoint right after initialization (the run then continues from t_start)       */
/*----------------------------------------------------------------------------------------------------------------*/
static void run_simulation(const Options *o, const char *root, double t_start, double t_end,
                           double t_pack, int stop_after_pack, const char *chkp_root, int unpack, Series *s)
{
   int    passed = 0, len = (int)strlen(o->input_file);
   const char *input_file = o->input_file;
   int    nOut = 0, nPL = 0, nDL = 0, nCh = 0;
   int    i, k, n, n0, iTip = 0, jTip = 0;
   char   rootbuf[BD_C_INTFSTRLEN], chkpbuf[BD_C_INTFSTRLEN];
   char  *names, *units;
   float *refpos, *plpos, *dlpos, *nodepos, *nodevel, *nodeacc, *ploads, *dloads, *chan;
   double *refori, *nodeori;
   float  react[6]; double disp[3], vel[6], acc[6], rootvel0[6];
   double ori[9], omega_x_r0[3], r0[3], d, dmax;
   double t, tn;

   memset(rootbuf, 0, sizeof(rootbuf));  strncpy(rootbuf, root, sizeof(rootbuf) - 1);
   memset(chkpbuf, 0, sizeof(chkpbuf));  if (chkp_root != NULL) strncpy(chkpbuf, chkp_root, sizeof(chkpbuf) - 1);

   /* Initial root velocity for a root rotating about the origin: omega x r0 */
   r0[0] = o->pos[0]; r0[1] = o->pos[1]; r0[2] = o->pos[2];
   cross(o->omega, r0, omega_x_r0);
   for (k = 0; k < 3; k++) { rootvel0[k] = omega_x_r0[k]; rootvel0[3+k] = o->omega[k]; }
   if (o->amp != 0.0) rootvel0[0] += (o->amp*2.0*M_PI*o->freq);

   names = (char*)malloc(BD_C_CHANNELBUFLEN);
   units = (char*)malloc(BD_C_CHANNELBUFLEN);

   BD_C_Init(&passed, &input_file, &len, rootbuf, o->pos, o->ori, rootvel0, o->grav, &o->dt, &o->interp,
             &o->dynamic, &o->wrout, &nOut, &nPL, &nDL, &nCh, names, units, &err_stat, err_msg);
   check("BD_C_Init");

   refpos  = (float*) malloc((size_t)3*nOut*sizeof(float));
   refori  = (double*)malloc((size_t)9*nOut*sizeof(double));
   plpos   = (float*) malloc((size_t)3*nPL *sizeof(float));
   dlpos   = (float*) malloc((size_t)3*nDL *sizeof(float));
   nodepos = (float*) malloc((size_t)3*nOut*sizeof(float));
   nodeori = (double*)malloc((size_t)9*nOut*sizeof(double));
   nodevel = (float*) malloc((size_t)6*nOut*sizeof(float));
   nodeacc = (float*) malloc((size_t)6*nOut*sizeof(float));
   ploads  = (float*) calloc((size_t)6*nPL, sizeof(float));
   dloads  = (float*) calloc((size_t)6*(nDL > 0 ? nDL : 1), sizeof(float));
   chan    = (float*) malloc((size_t)(nCh > 0 ? nCh : 1)*sizeof(float));

   BD_C_GetRefPositions(refpos, refori, plpos, dlpos, &err_stat, err_msg);
   check("BD_C_GetRefPositions");

   /* Tip: the output node and point load node farthest from the first output node */
   dmax = -1.0;
   for (i = 0; i < nOut; i++) {
      d = 0.0; for (k = 0; k < 3; k++) d += pow(refpos[3*i+k] - refpos[k], 2);
      if (d > dmax) { dmax = d; iTip = i; }
   }
   s->length = sqrt(dmax);
   dmax = -1.0;
   for (i = 0; i < nPL; i++) {
      d = 0.0; for (k = 0; k < 3; k++) d += pow(plpos[3*i+k] - refpos[k], 2);
      if (d > dmax) { dmax = d; jTip = i; }
   }

   printf("   Beam: %d output nodes, %d point load nodes, %d distributed load nodes, %d output channels, length %.6f m\n",
          nOut, nPL, nDL, nCh, s->length);

   if (unpack) {
      BD_C_UnpackStates(chkpbuf, &err_stat, err_msg);
      check("BD_C_UnpackStates");
   }

#define SET_INPUTS(time)                                                                             \
   do {                                                                                              \
      int kk;                                                                                        \
      root_motion(o, (time), disp, ori, vel, acc);                                                  \
      BD_C_SetRootMotion(disp, ori, vel, acc, &err_stat, err_msg);  check("BD_C_SetRootMotion");    \
      for (kk = 0; kk < 3; kk++) {                                                                   \
         ploads[6*jTip+kk]   = (o->t_release < 0.0 || (time) < o->t_release - 0.5*o->dt) ? (float)o->tip[kk]  : 0.0f; \
         ploads[6*jTip+3+kk] = (o->t_release < 0.0 || (time) < o->t_release - 0.5*o->dt) ? (float)o->tipm[kk] : 0.0f; \
      }                                                                                              \
      BD_C_SetPointLoads(ploads, &err_stat, err_msg);  check("BD_C_SetPointLoads");                 \
      if (o->use_distr) {                                                                            \
         int ii;                                                                                     \
         for (ii = 0; ii < nDL; ii++) for (kk = 0; kk < 6; kk++) dloads[6*ii+kk] = (float)o->distr[kk]; \
         BD_C_SetDistrLoads(dloads, &err_stat, err_msg);  check("BD_C_SetDistrLoads");              \
      }                                                                                              \
   } while (0)

#define RECORD(time)                                                                                 \
   do {                                                                                              \
      int kk;                                                                                        \
      if (s->n < s->nmax) {                                                                          \
         s->t[s->n] = (time);                                                                        \
         for (kk = 0; kk < 3; kk++) s->tip[3*s->n+kk] = (double)nodepos[3*iTip+kk] - (double)refpos[3*iTip+kk]; \
         for (kk = 0; kk < 6; kk++) s->react[6*s->n+kk] = (double)react[kk];                        \
         s->n++;                                                                                     \
      }                                                                                              \
   } while (0)

   n0 = (int)lround(t_start/o->dt);
   t  = n0*o->dt;

   SET_INPUTS(t);
   BD_C_CalcOutput(&t, nodepos, nodeori, nodevel, nodeacc, react, chan, &err_stat, err_msg);
   check("BD_C_CalcOutput");
   RECORD(t);

   for (n = n0; ; n++) {
      t  = n*o->dt;
      tn = (n+1)*o->dt;
      if (t_pack >= 0.0 && fabs(t - t_pack) < 0.5*o->dt) {
         BD_C_PackStates(chkpbuf, &err_stat, err_msg);
         check("BD_C_PackStates");
         printf("   Checkpoint written at t = %.6f s\n", t);
         if (stop_after_pack) break;
      }
      if (tn > t_end + 0.5*o->dt) break;

      SET_INPUTS(tn);
      BD_C_UpdateStates(&t, &tn, &err_stat, err_msg);
      check("BD_C_UpdateStates");
      BD_C_CalcOutput(&tn, nodepos, nodeori, nodevel, nodeacc, react, chan, &err_stat, err_msg);
      check("BD_C_CalcOutput");
      RECORD(tn);
   }
#undef SET_INPUTS
#undef RECORD

   BD_C_End(&err_stat, err_msg);
   check("BD_C_End");

   free(names); free(units); free(refpos); free(refori); free(plpos); free(dlpos);
   free(nodepos); free(nodeori); free(nodevel); free(nodeacc); free(ploads); free(dloads); free(chan);
}

/*----------------------------------------------------------------------------------------------------------------*/
/* Analysis helpers                                                                                                 */
/*----------------------------------------------------------------------------------------------------------------*/
/* Frequency and damping ratio of a decaying oscillation x(t) about zero for t >= t0, from the upward zero crossings
 * and the logarithmic decrement of the positive peaks.  Returns the number of cycles used (0 on failure). */
static int measure_oscillation(const Series *s, int comp, double t0, double *f_damped, double *zeta)
{
   int    i, ncross = 0, npeak = 0;
   double t_first = 0.0, t_last = 0.0, a_first = 0.0, a_last = 0.0;
   for (i = 1; i < s->n - 1; i++) {
      double x0 = s->tip[3*(i-1)+comp], x1 = s->tip[3*i+comp], x2 = s->tip[3*(i+1)+comp];
      if (s->t[i] < t0) continue;
      if (x0 < 0.0 && x1 >= 0.0) {                                   /* upward zero crossing, linear interpolation */
         double tc = s->t[i-1] + (0.0 - x0)/(x1 - x0)*(s->t[i] - s->t[i-1]);
         if (ncross == 0) t_first = tc;
         t_last = tc;
         ncross++;
      }
      if (x1 > 0.0 && x1 >= x0 && x1 > x2) {                         /* positive peak */
         if (npeak == 0) a_first = x1;
         a_last = x1;
         npeak++;
      }
   }
   if (ncross < 3 || npeak < 2) return 0;
   *f_damped = (ncross - 1)/(t_last - t_first);
   {
      double delta = log(a_first/a_last)/(npeak - 1);                 /* logarithmic decrement */
      *zeta = delta/sqrt(4.0*M_PI*M_PI + delta*delta);
   }
   return ncross - 1;
}

static int report(const char *what, double value, double reference, double tol)
{
   double err = (reference != 0.0) ? fabs(value - reference)/fabs(reference) : fabs(value - reference);
   int ok = err <= tol;
   printf("   %-44s %16.8e   reference %16.8e   rel. diff %10.3e   %s\n", what, value, reference, err, ok ? "PASS" : "FAIL");
   return ok ? 0 : 1;
}

/*----------------------------------------------------------------------------------------------------------------*/
/* Modes                                                                                                            */
/*----------------------------------------------------------------------------------------------------------------*/
static int mode_run(const Options *o)
{
   Series s;
   int    i, n = (int)lround(o->tmax/o->dt) + 2, fails = 0;
   series_alloc(&s, n);
   run_simulation(o, o->root, 0.0, o->tmax, -1.0, 0, NULL, 0, &s);

   printf("   Final time %.6f s: tip displacement [%.8e %.8e %.8e] m, root reaction force [%.8e %.8e %.8e] N\n",
          s.t[s.n-1], s.tip[3*(s.n-1)], s.tip[3*(s.n-1)+1], s.tip[3*(s.n-1)+2],
          s.react[6*(s.n-1)], s.react[6*(s.n-1)+1], s.react[6*(s.n-1)+2]);

   if (o->amp != 0.0 && o->freq > 0.0) {
      /* Prescribed sinusoidal root motion: the root reaction must follow the root acceleration.  Over the last
       * period, compare the amplitude of the x reaction force with the rigid-body inertial estimate
       * mu*L*amp*(2*pi*f)^2 (times the single-mode dynamic amplification if EI is given). */
      double T = 1.0/o->freq, fmax = 0.0, amax = 0.0;
      for (i = 0; i < s.n; i++) {
         if (s.t[i] < o->tmax - T) continue;
         if (fabs(s.react[6*i]) > fmax) fmax = fabs(s.react[6*i]);
         if (fabs(s.tip[3*i]) > amax) amax = fabs(s.tip[3*i]);
      }
      printf("   Root motion x = %.4e*sin(2*pi*%.4f*t): max |root Fx| over the last period %.8e N, max |tip x displacement| %.8e m\n",
             o->amp, o->freq, fmax, amax);
      if (fmax <= 0.0) { printf("   FAIL: the root reaction does not respond to the root motion\n"); fails++; }
      if (o->mu > 0.0) {
         double M = o->mu*s.length, a = o->amp*pow(2.0*M_PI*o->freq, 2), est = M*a;
         if (o->EI > 0.0) {
            double f1 = 1.8751040687*1.8751040687/(2.0*M_PI)*sqrt(o->EI/(o->mu*pow(s.length, 4)));
            est /= fabs(1.0 - pow(o->freq/f1, 2));
         }
         fails += report("max |root Fx| vs. inertial estimate", fmax, est, o->tol);
      }
   }
   series_free(&s);
   return fails;
}

static int mode_cantilever(const Options *o)
{
   Options os = *o, od = *o;
   Series  s;
   int     n, i, fails = 0, ncyc;
   double  L, F, d_static, d_settled = 0.0, d_eb, f_eb, f_d = 0.0, zeta = 0.0, f_n;
   char    root[BD_C_INTFSTRLEN];

   if (o->EI <= 0.0 || o->mu <= 0.0 || o->tip[0] == 0.0) {
      fprintf(stderr, "The cantilever mode needs EI=, mu=, and tip=Fx,0,0\n");
      return 1;
   }
   F = o->tip[0];

   /* (a) static solve: tip deflection under the tip force */
   printf("\n-- Static solve: tip force %.4e N\n", F);
   os.dynamic = 0; os.t_release = -1.0; os.amp = 0.0;
   snprintf(root, sizeof(root), "%s_static", o->root);
   series_alloc(&s, 20);
   run_simulation(&os, root, 0.0, 5.0*o->dt, -1.0, 0, NULL, 0, &s);
   L        = s.length;
   d_static = s.tip[3*(s.n-1)];
   series_free(&s);
   d_eb = F*L*L*L/(3.0*o->EI);
   fails += report("static tip deflection (static solve)", d_static, d_eb, o->tol);

   /* (b) dynamic solve: settle under the tip force, release, measure the first bending frequency from the tip motion */
   od.dynamic = 1; od.amp = 0.0;
   if (od.t_release < 0.0) od.t_release = 0.5*o->tmax;
   printf("\n-- Dynamic solve: tip force applied until t = %.3f s, free vibration until t = %.3f s\n", od.t_release, od.tmax);
   snprintf(root, sizeof(root), "%s_dynamic", o->root);
   n = (int)lround(o->tmax/o->dt) + 2;
   series_alloc(&s, n);
   run_simulation(&od, root, 0.0, od.tmax, -1.0, 0, NULL, 0, &s);
   for (i = 0; i < s.n; i++) if (s.t[i] <= od.t_release - 0.5*o->dt) d_settled = s.tip[3*i];
   fails += report("static tip deflection (settled dynamic solve)", d_settled, d_eb, o->tol);

   f_eb = 1.8751040687*1.8751040687/(2.0*M_PI)*sqrt(o->EI/(o->mu*pow(L, 4)));
   ncyc = measure_oscillation(&s, 0, od.t_release + 1.0/f_eb, &f_d, &zeta);
   if (ncyc == 0) {
      printf("   FAIL: could not measure the free vibration (not enough cycles after the release)\n");
      fails++;
   } else {
      f_n = f_d/sqrt(1.0 - zeta*zeta);
      printf("   Free vibration: %d cycles, damped frequency %.6f Hz, damping ratio %.5f\n", ncyc, f_d, zeta);
      fails += report("first bending frequency (undamped)", f_n, f_eb, o->tol);
   }
   series_free(&s);
   return fails;
}

static int mode_checkpoint(const Options *o)
{
   Options oc = *o;
   Series  s1, s2;
   int     n = (int)lround(o->tmax/o->dt) + 2, i, j, k, n0, fails = 0;
   double  dtip = 0.0, dreact = 0.0, tipmax = 0.0, reactmax = 0.0;
   char    root[BD_C_INTFSTRLEN], chkp[BD_C_INTFSTRLEN];

   if (oc.t_chk < 0.0) oc.t_chk = 0.5*o->tmax;
   oc.t_chk = lround(oc.t_chk/o->dt)*o->dt;

   printf("\n-- Straight run to t = %.3f s\n", o->tmax);
   series_alloc(&s1, n);
   snprintf(root, sizeof(root), "%s_straight", o->root);
   run_simulation(&oc, root, 0.0, o->tmax, -1.0, 0, NULL, 0, &s1);

   printf("\n-- Run to t = %.3f s and write a checkpoint\n", oc.t_chk);
   snprintf(root, sizeof(root), "%s_part1", o->root);
   snprintf(chkp, sizeof(chkp), "%s_checkpoint", o->root);
   series_alloc(&s2, n);
   run_simulation(&oc, root, 0.0, oc.t_chk, oc.t_chk, 1, chkp, 0, &s2);
   series_free(&s2);

   printf("\n-- Restore the checkpoint into a fresh instance and continue to t = %.3f s\n", o->tmax);
   series_alloc(&s2, n);
   snprintf(root, sizeof(root), "%s_part2", o->root);
   run_simulation(&oc, root, oc.t_chk, o->tmax, -1.0, 0, chkp, 1, &s2);

   /* Compare the restarted run with the straight run at the same times */
   n0 = (int)lround(oc.t_chk/o->dt);
   for (j = 0; j < s2.n; j++) {
      i = n0 + j;
      if (i >= s1.n) break;
      if (fabs(s1.t[i] - s2.t[j]) > 0.5*o->dt) { printf("   FAIL: time mismatch %f vs %f\n", s1.t[i], s2.t[j]); fails++; break; }
      for (k = 0; k < 3; k++) { dtip   = fmax(dtip,   fabs(s1.tip[3*i+k]   - s2.tip[3*j+k]));   tipmax   = fmax(tipmax,   fabs(s1.tip[3*i+k])); }
      for (k = 0; k < 6; k++) { dreact = fmax(dreact, fabs(s1.react[6*i+k] - s2.react[6*j+k])); reactmax = fmax(reactmax, fabs(s1.react[6*i+k])); }
   }
   printf("   Compared %d output times after the checkpoint\n", j);
   printf("   max |tip displacement difference|  %.3e m   (max |tip displacement| %.3e m)\n", dtip, tipmax);
   printf("   max |root reaction difference|     %.3e     (max |root reaction| %.3e)\n", dreact, reactmax);
   fails += report("restart vs. straight run, tip displacement diff", dtip, 0.0, 1e-6*(tipmax > 0.0 ? tipmax : 1.0));
   fails += report("restart vs. straight run, root reaction diff",   dreact, 0.0, 1e-6*(reactmax > 0.0 ? reactmax : 1.0));

   series_free(&s1);
   series_free(&s2);
   return fails;
}

/*----------------------------------------------------------------------------------------------------------------*/
/* BeamDyn driver input file (the input file of beamdyn_driver)                                                     */
/*----------------------------------------------------------------------------------------------------------------*/
static int dvr_line(FILE *fp, char *line, size_t n, const char *what)
{
   if (fgets(line, (int)n, fp) == NULL) { fprintf(stderr, "Driver input file: premature end of file while reading %s\n", what); return 1; }
   return 0;
}

static int dvr_skip(FILE *fp, int nlines, const char *what)
{
   char line[2048];
   int  i;
   for (i = 0; i < nlines; i++) if (dvr_line(fp, line, sizeof(line), what)) return 1;
   return 0;
}

/* Read nvals numbers from the start of the next line */
static int dvr_numbers(FILE *fp, double *v, int nvals, const char *what)
{
   char  line[2048], *s, *end;
   int   i;
   if (dvr_line(fp, line, sizeof(line), what)) return 1;
   s = line;
   for (i = 0; i < nvals; i++) {
      v[i] = strtod(s, &end);
      if (end == s) { fprintf(stderr, "Driver input file: invalid numerical input while reading %s\n", what); return 1; }
      s = end;
      while (*s == ',' ) s++;
   }
   return 0;
}

/* Read a Fortran logical (T, F, True, False, .true., .false.) from the start of the next line */
static int dvr_logical(FILE *fp, int *v, const char *what)
{
   char line[2048], *s;
   if (dvr_line(fp, line, sizeof(line), what)) return 1;
   s = line;
   while (*s == ' ' || *s == '\t' || *s == '.' || *s == '"') s++;
   if      (*s == 'T' || *s == 't') *v = 1;
   else if (*s == 'F' || *s == 'f') *v = 0;
   else { fprintf(stderr, "Driver input file: invalid logical input while reading %s\n", what); return 1; }
   return 0;
}

/* Read a BeamDyn driver input file into the options.  The primary input file name is returned in primary (resolved
 * relative to the directory of the driver input file) and the output root name is the driver file name without its
 * extension, as in beamdyn_driver. */
static int read_driver_file(const char *file, Options *o, char *primary, size_t nprimary)
{
   FILE  *fp = fopen(file, "r");
   char   line[2048], name[1024], *s, *e;
   double v[3], t_initial, t_final;
   int    i, k, flag, npl, nfinal;
   const char *slash;

   if (fp == NULL) { fprintf(stderr, "Cannot open the driver input file %s\n", file); return 1; }

#define DVR(call) do { if (call) { fclose(fp); return 1; } } while (0)
   DVR(dvr_skip(fp, 3, "the file header"));
   DVR(dvr_logical(fp, &o->dynamic, "DynamicSolve"));
   DVR(dvr_numbers(fp, &t_initial, 1, "t_initial"));
   DVR(dvr_numbers(fp, &t_final, 1, "t_final"));
   DVR(dvr_numbers(fp, &o->dt, 1, "dt"));
   DVR(dvr_skip(fp, 1, "the gravity section header"));
   for (k = 0; k < 3; k++) DVR(dvr_numbers(fp, &o->grav[k], 1, "the gravity vector"));
   DVR(dvr_skip(fp, 1, "the frame section header"));
   for (k = 0; k < 3; k++) DVR(dvr_numbers(fp, &o->pos[k], 1, "GlbPos"));
   DVR(dvr_skip(fp, 2, "the DCM comment lines"));
   for (i = 0; i < 3; i++) { DVR(dvr_numbers(fp, v, 3, "the root orientation DCM")); for (k = 0; k < 3; k++) o->ori[3*i+k] = v[k]; }
   DVR(dvr_logical(fp, &flag, "GlbRotBladeT0"));
   if (!flag) { fprintf(stderr, "Driver input file: GlbRotBladeT0 = False is not supported by this driver\n"); fclose(fp); return 1; }
   DVR(dvr_skip(fp, 1, "the root velocity section header"));
   for (k = 0; k < 3; k++) DVR(dvr_numbers(fp, &o->omega[k], 1, "RootVel"));
   DVR(dvr_skip(fp, 1, "the applied force section header"));
   for (k = 0; k < 6; k++) DVR(dvr_numbers(fp, &o->distr[k], 1, "DistrLoad"));
   for (k = 0; k < 3; k++) DVR(dvr_numbers(fp, &o->tip[k], 1, "TipLoad"));
   for (k = 0; k < 3; k++) DVR(dvr_numbers(fp, &o->tipm[k], 1, "TipLoad"));
   DVR(dvr_numbers(fp, v, 1, "NumPointLoads"));
   npl = (int)v[0];
   if (npl > 0) { fprintf(stderr, "Driver input file: the multi-point load table (NumPointLoads > 0) is not supported by this driver\n"); fclose(fp); return 1; }
   DVR(dvr_skip(fp, 3, "the point load table header and the primary input file section header"));
   DVR(dvr_line(fp, line, sizeof(line), "InputFile"));
#undef DVR
   fclose(fp);

   /* primary input file name: first token on the line, optionally quoted */
   s = line;
   while (*s == ' ' || *s == '\t') s++;
   if (*s == '"' || *s == '\'') { char q = *s++; e = strchr(s, q); } else { e = s; while (*e && *e != ' ' && *e != '\t' && *e != '\n' && *e != '\r') e++; }
   if (e == NULL || e == s) { fprintf(stderr, "Driver input file: cannot read the primary input file name\n"); return 1; }
   snprintf(name, sizeof(name), "%.*s", (int)(e - s), s);

   /* resolve relative to the directory of the driver input file; root name is the driver file without extension */
   slash = strrchr(file, '/');
   if (slash != NULL && name[0] != '/') snprintf(primary, nprimary, "%.*s/%s", (int)(slash - file), file, name);
   else                                 snprintf(primary, nprimary, "%s", name);
   snprintf(o->root, sizeof(o->root), "%s", file);
   s = strrchr(o->root, '.');
   if (s != NULL && (slash == NULL || s > o->root + (slash - file))) *s = '\0';

   if (t_initial != 0.0) { fprintf(stderr, "Driver input file: t_initial must be 0 for this driver\n"); return 1; }
   for (k = 0; k < 6; k++) if (o->distr[k] != 0.0) o->use_distr = 1;

   /* beamdyn_driver takes INT((t_final - t_initial)/dt) + 1 steps; a static solve stops after 9 load steps */
   nfinal = (int)(t_final/o->dt) + 1;
   if (!o->dynamic && nfinal > 9) nfinal = 9;
   o->tmax  = nfinal*o->dt;
   o->wrout = 1;
   return 0;
}

/*----------------------------------------------------------------------------------------------------------------*/
int main(int argc, char **argv)
{
   Options o;
   int     fails;
   char    primary[2048];

   if (argc == 2) {
      /* a BeamDyn driver input file */
      set_defaults(&o, "run");
      if (read_driver_file(argv[1], &o, primary, sizeof(primary))) return EXIT_FAILURE;
      o.input_file = primary;
      fails = mode_run(&o);
      printf("\n%s: %s\n", argv[1], fails == 0 ? "completed" : "FAILED");
      return fails == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
   }
   if (argc < 3) {
      fprintf(stderr, "Usage: %s <BeamDyn driver input file>\n       %s <run|cantilever|checkpoint> <BeamDyn primary input file> [key=value ...]\n", argv[0], argv[0]);
      return EXIT_FAILURE;
   }
   set_defaults(&o, argv[1]);
   o.input_file = argv[2];
   if (parse_options(argc, argv, &o)) return EXIT_FAILURE;

   if      (strcmp(argv[1], "run") == 0)        fails = mode_run(&o);
   else if (strcmp(argv[1], "cantilever") == 0) fails = mode_cantilever(&o);
   else if (strcmp(argv[1], "checkpoint") == 0) fails = mode_checkpoint(&o);
   else { fprintf(stderr, "Unknown mode '%s'\n", argv[1]); return EXIT_FAILURE; }

   printf("\n%s: %s\n", argv[1], fails == 0 ? "all checks passed" : "some checks FAILED");
   return fails == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
