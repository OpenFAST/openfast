#ifndef DATA_H
#include "registry.h"

typedef struct node_struct {

  int     node_kind ;
  int     type_type ;
  char          name[NAMELEN] ;
  char          mapsto[NAMELEN] ;
  char          nickname[NAMELEN] ;
  struct node_struct  * fields ;
  struct node_struct  * type ;
  int           ndims ;
  struct node_struct  * dims[MAXDIMS] ;
  int     proc_orient ;    /* ALL_[ZXY]_ON_PROC which dimension is all on processor */
  int           ntl ;
  int           subject_to_communication ;
  int           boundary_array ;
  int           boundary_array_4d ;
  char    use[NAMELEN] ;
  char    inival[NAMELEN] ;
  char    descrip[NAMELEN] ;
  char    units[NAMELEN] ;

/* Fields for 4D scalar arrays */
  int           scalar_array_member ;
  int           has_scalar_array_tendencies ;
  struct node_struct * members ;

/* I/O flags */
  int     restart ;
  int     boundary   ;
  int     namelist   ;
  char    namelistsection[NAMELEN] ;

/* Fields for Modname */
  struct node_struct * module_ddt_list ;

#ifdef FUTURE
  unsigned int     nest_mask ;
  struct node_struct * next4d ;
  char force_aux_fields[2048] ;
  char force_fcn_name[2048] ;
  char interpd_aux_fields[2048] ;
  char interpd_fcn_name[2048] ;
  char interpu_aux_fields[2048] ;
  char interpu_fcn_name[2048] ;
  char smoothu_fcn_name[2048] ;
  char smoothu_aux_fields[2048] ;
#endif

/* CTRL */
  int gen_wrapper ;
  struct node_struct * next ;

/* fields used by rconfig nodes */
  char nentries[NAMELEN] ;
  char howset[NAMELEN] ;
  char dflt[NAMELEN] ;

/* fields used by Dim nodes */

  char dim_name[32] ;
  char dim_data_name[NAMELEN] ;
  int  coord_axis ;   /* X, Y, Z, C */
                                 /* DOMAIN_STANDARD, NAMELIST, CONSTANT */
  int  len_defined_how ;  
  char assoc_nl_var_s[NAMELEN] ;  /* for NAMELIST */
  char assoc_nl_var_e[NAMELEN] ;  /* for NAMELIST */
  int  coord_start ;               /* for CONSTANT */
  int  coord_end ;                 /* for CONSTANT */
  int  dim_order ;                 /* order that dimensions are specified
                                      in framework */
  int  subgrid ;                  /* 1=subgrid dimension */
  int  deferred ;                 /* a deferred-shape dimension, that is, a colon */

/* fields used by Package nodes */
  char pkg_assoc[NAMELEN] ;
  char pkg_statevars[NAMELEN] ;
  char pkg_4dscalars[NAMELEN_LONG] ;

/* fields used by Comm (halo, period, xpose)  nodes */
  char comm_define[2*8192] ;

/* marker */
  int mark ;

} node_t ;

#ifndef DEFINE_GLOBALS
#  define EXTERN extern
#else
#  define EXTERN
#endif

EXTERN int sw_deref_kludge ;
EXTERN int sw_io_deref_kludge ;
EXTERN int sw_3dvar_iry_kludge ;
EXTERN int sw_distrib_io_layer ;
EXTERN int sw_limit_args ;
EXTERN int sw_dm_parallel  ;
EXTERN int sw_move  ;
EXTERN int sw_all_x_staggered ;
EXTERN int sw_all_y_staggered ;
EXTERN int sw_dm_serial_in_only ;
EXTERN int sw_fort_kludge ;
EXTERN char sw_commpath[NAMELEN] ;
EXTERN int sw_new_bdys ;  /* 20070207 JM support decomposed boundary arrays */
EXTERN int sw_unidir_shift_halo ;  /* 20100210 JM assume that halo to shift is same in both directions and only gen one of them */
EXTERN int sw_new_with_old_bdys ;  /* 20070207 JM for debugging interim phase, new comms w/ old data structs */

EXTERN node_t * Type ;
EXTERN node_t * Dim ;
EXTERN node_t * Packages ;
EXTERN node_t * Halos ;
EXTERN node_t * Periods ;
EXTERN node_t * Xposes ;
EXTERN node_t * FourD ;
EXTERN node_t * Swaps ;
EXTERN node_t * Cycles ;
EXTERN node_t * ModNames ;

EXTERN node_t Domain ;

EXTERN char t1[NAMELEN], t2[NAMELEN], t3[NAMELEN], t4[NAMELEN], t5[NAMELEN], t6[NAMELEN] ;
EXTERN char thiscom[4*NAMELEN] ;

EXTERN int max_time_level  ;  /* Maximum number of time levels of any state variable */

#define  P_XSB  1
#define  P_XEB  2
#define  P_YSB  3
#define  P_YEB  4


#define DATA_H
#endif
