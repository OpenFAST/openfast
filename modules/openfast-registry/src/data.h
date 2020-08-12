#ifndef DATA_H
#include "registry.h"

typedef struct node_struct {

  int     node_kind ;
  int     type_type ;
  char          name[NAMELEN] ;
  char          mapsto[NAMELEN] ;
  char          nickname[NAMELEN] ;
  struct node_struct  * fields ;
  struct node_struct  * params ;
  struct node_struct  * type ;
  struct node_struct  * module ;  /* type node pointer back to module node it is defined in */
  int    max_ndims;    // max number of dimensions (so we don't have hundreds of unused variables that produce warnings)
  int    containsPtr;  // if contains a pointer in type/subtype
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

/* I/O flags */
  int     restart ;
  int     boundary   ;
  int     namelist   ;
  char    namelistsection[NAMELEN] ;

/* Fields for Modname */
  struct node_struct * module_ddt_list ;


/* CTRL */
  int gen_periodic ;
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
  int  dim_param;                  /* for using PARAMETER dimension */
  char dim_param_name[NAMELEN];    /* for using PARAMETER dimension */

  int  dim_order ;                 /* order that dimensions are specified
                                      in framework */
  int  subgrid ;                  /* 1=subgrid dimension */
  int  deferred ;                 /* a deferred-shape dimension, that is, a colon */

  int  usefrom ;

/* fields used by Package nodes */
  char pkg_assoc[NAMELEN] ;
  char pkg_statevars[NAMELEN] ;
  char pkg_4dscalars[NAMELEN_LONG] ;

/* fields used by Comm (halo, period, xpose)  nodes */
  char comm_define[2*8192] ;

  int is_interface_type ;

/* marker */
  int mark ;

} node_t ;

#ifndef DEFINE_GLOBALS
#  define EXTERN extern
#else
#  define EXTERN
#endif

EXTERN int sw_output_template_force ;
EXTERN char sw_commpath[NAMELEN] ;
EXTERN char sw_modname_subst[NAMELEN] ;
EXTERN char sw_modnickname_subst[NAMELEN] ;
EXTERN int sw_new_bdys ;  /* 20070207 JM support decomposed boundary arrays */
EXTERN int sw_unidir_shift_halo ;  /* 20100210 JM assume that halo to shift is same in both directions and only gen one of them */
EXTERN int sw_new_with_old_bdys ;  /* 20070207 JM for debugging interim phase, new comms w/ old data structs */
EXTERN int sw_norealloc_lsh;  /* 20070207 addresses compilers like gfortran that do not /assume:realloc_lhs */
EXTERN int sw_ccode ;           /* 20130523 generate C code too */
EXTERN int sw_noextrap;
EXTERN char sw_shownodes ;

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
EXTERN char thiscom[NAMELEN] ;

EXTERN int max_time_level  ;  /* Maximum number of time levels of any state variable */

#define MAXINCLDIRS 50 
EXTERN int   nincldirs ;
EXTERN char IncludeDirs[MAXINCLDIRS][NAMELEN] ;
EXTERN char OutDir[NAMELEN];

#define  P_XSB  1
#define  P_XEB  2
#define  P_YSB  3
#define  P_YEB  4


#define DATA_H
#endif
