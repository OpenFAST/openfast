#ifndef REGISTRY_H
#define NAMELEN 512        
#define NAMELEN_LONG 12500 /*changed from 8192 to 12500 by PNNL on 12/22/2010*/
#define MAXDIMS 21
#define MAX_DYNCORES 50   /* ha ha, just kidding */
/* #define MAX_ARGLINE 175    WRF uses 128 by default, but the nested chem version hit the continuation line limit for efc so it had to be increased, wig 14-Oct-2004 */
#define MAX_ARGLINE 128   /* welp, 175 means lines longer than 130 chars, which is a Fortran no no */
#define MAX_TYPEDEFS 50   /* typedef history -ajb */
#define MAXTOKENS 100

/* defines of system commands */
#define UNIQSORT "/bin/sort -u"
#define CATCOMM  "/bin/cat"
#define RMCOMM   "/bin/rm"
#define MVCOMM   "/bin/mv"

#define DRIVER_LAYER     100
#define MEDIATION_LAYER  200

enum coord_axis      { COORD_X , COORD_Y , COORD_Z , COORD_C } ;
enum len_defined_how { DOMAIN_STANDARD , NAMELIST , CONSTANT } ;
enum type_type       { SIMPLE , DERIVED } ;
enum proc_orient     { ALL_Z_ON_PROC , ALL_X_ON_PROC , ALL_Y_ON_PROC } ;

/* wrapping options */
#define WRAP_HIDDEN_FIELD      2
#define WRAP_EXPOSED_FIELD     1
#define WRAP_NONE              0


/* node_kind  mask settings */
#define FIELD      1
#define PARAM      2
#define RCONFIG    4
#define FOURD      8
#define MEMBER    16
#define TYPE      32
#define DIM       64
#define MODNAME  128
#define HALO     256
#define PERIOD   512
#define SWAP    1024
#define CYCLE   2048
#define XPOSE   4096
#define FOURD1  8192
#define BDYONLY 16384

#define RESTART       0x02000000      /*   25 */
#define BOUNDARY      0x04000000      /*   26 */
#define INTERP_DOWN   0x08000000      /*   27 */
#define FORCE_DOWN    0x10000000      /*   28 */
#define INTERP_UP     0x20000000      /*   29 */
#define SMOOTH_UP     0x40000000      /*   20 */
#define METADATA      0x80000000      /*   31 */


#define REGISTRY_H
#endif

#ifdef WIN32
#define snprintf _snprintf
#endif

