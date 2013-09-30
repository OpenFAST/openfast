//!STARTOFREGISTRYGENERATEDFILE './MAP_Types.c'
//
//! WARNING This file is generated automatically by the FAST registry
//! Do not edit.  Your changes to this file will be lost.
//!
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MAP_TypesWin.h"
#include <windows.h >
//#include "MAP_FortranBindingWin.h"

#ifdef _WIN32
#include "stdbool.h"
#define CALL __declspec( dllexport ) //define something for Windows (64-bit)
#elif _WIN64
#include "stdbool.h"
#define CALL __declspec( dllexport ) //define something for Windows (64-bit)
#else
#include <stdbool.h>
#endif


//#define CALL __attribute__((dllexport) )
//MAP_InitInputType_t* CALL MAP_InitInput_Create() { return ((MAP_InitInputType_t*) malloc( sizeof(MAP_InitInputType_t()))) ; } ;
//void CALL MAP_InitInput_Delete(MAP_InitInputType_t *This) { free(This) ; } ;
//MAP_InitOutputType_t* CALL MAP_InitOutput_Create() { return ((MAP_InitOutputType_t*) malloc( sizeof(MAP_InitOutputType_t()))) ; } ;
//void CALL MAP_InitOutput_Delete(MAP_InitOutputType_t *This) { free(This) ; } ;
//MAP_ContinuousStateType_t* CALL MAP_ContState_Create() { return ((MAP_ContinuousStateType_t*) malloc( sizeof(MAP_ContinuousStateType_t()))) ; } ;
//void CALL MAP_ContState_Delete(MAP_ContinuousStateType_t *This) { free(This) ; } ;
//MAP_DiscreteStateType_t* CALL MAP_DiscState_Create() { return ((MAP_DiscreteStateType_t*) malloc( sizeof(MAP_DiscreteStateType_t()))) ; } ;
//void CALL MAP_DiscState_Delete(MAP_DiscreteStateType_t *This) { free(This) ; } ;
//MAP_OtherStateType_t* CALL MAP_OtherState_Create() { return ((MAP_OtherStateType_t*) malloc( sizeof(MAP_OtherStateType_t()))) ; } ;
//void CALL MAP_OtherState_Delete(MAP_OtherStateType_t *This) { free(This) ; } ;
//MAP_ConstraintStateType_t* CALL MAP_ConstrState_Create() { return ((MAP_ConstraintStateType_t*) malloc( sizeof(MAP_ConstraintStateType_t()))) ; } ;
//void CALL MAP_ConstrState_Delete(MAP_ConstraintStateType_t *This) { free(This) ; } ;
//MAP_ParameterType_t* CALL MAP_Param_Create() { return ((MAP_ParameterType_t*) malloc( sizeof(MAP_ParameterType_t()))) ; } ;
//void CALL MAP_Param_Delete(MAP_ParameterType_t *This) { free(This) ; } ;
//MAP_InputType_t* CALL MAP_Input_Create() { return ((MAP_InputType_t*) malloc( sizeof(MAP_InputType_t()))) ; } ;
//void CALL MAP_Input_Delete(MAP_InputType_t *This) { free(This) ; } ;
//MAP_OutputType_t* CALL MAP_Output_Create() { return ((MAP_OutputType_t*) malloc( sizeof(MAP_OutputType_t()))) ; } ;
//void CALL MAP_Output_Delete(MAP_OutputType_t *This) { free(This) ; } ;

int
C_MAP_PackInitInput( float * ReKiBuf,  int * Re_BufSz ,
                 double * DbKiBuf, int * Db_BufSz ,
                 int * IntKiBuf,   int * Int_BufSz ,
                 MAP_InitInputType_t *InData, char * ErrMsg, int *SizeOnly )
{
  int ErrStat ;
  int OnlySize ;
  int Re_BufSz2 ;
  int Db_BufSz2 ;
  int Int_BufSz2 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int one         = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes and subtypes, if any

  OnlySize = *SizeOnly ;

  *Re_BufSz = 0 ;
  *Db_BufSz = 0 ;
  *Int_BufSz = 0 ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  *Db_BufSz   += 1  ; // gravity
  *Db_BufSz   += 1  ; // sea_density
  *Db_BufSz   += 1  ; // depth
  if ( ! OnlySize ) {
    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
    DbKiBuf[Db_Xferred++] = InData->gravity ;
    DbKiBuf[Db_Xferred++] = InData->sea_density ;
    DbKiBuf[Db_Xferred++] = InData->depth ;
  }
  return(ErrStat) ;
}

int
C_MAP_UnpackInitInput( float * ReKiBuf,  
                 double * DbKiBuf, 
                 int * IntKiBuf,   
                 MAP_InitInputType_t *OutData, char * ErrMsg )
{
  int ErrStat ;
  int Re_BufSz2 = 0 ;
  int Db_BufSz2 = 0 ;
  int Int_BufSz2 = 0 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int Re_CurrSz = 0 ;
  int Db_CurrSz = 0 ;
  int Int_CurrSz = 0 ;
  int one        = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes, if any
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  OutData->gravity = DbKiBuf [ Db_Xferred ] ; 
  Db_Xferred   = Db_Xferred   + 1 ; 
  OutData->sea_density = DbKiBuf [ Db_Xferred ] ; 
  Db_Xferred   = Db_Xferred   + 1 ; 
  OutData->depth = DbKiBuf [ Db_Xferred ] ; 
  Db_Xferred   = Db_Xferred   + 1 ; 
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}

int
C_MAP_PackInitOutput( float * ReKiBuf,  int * Re_BufSz ,
                 double * DbKiBuf, int * Db_BufSz ,
                 int * IntKiBuf,   int * Int_BufSz ,
                 MAP_InitOutputType_t *InData, char * ErrMsg, int *SizeOnly )
{
  int ErrStat ;
  int OnlySize ;
  int Re_BufSz2 ;
  int Db_BufSz2 ;
  int Int_BufSz2 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int one         = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes and subtypes, if any

  OnlySize = *SizeOnly ;

  *Re_BufSz = 0 ;
  *Db_BufSz = 0 ;
  *Int_BufSz = 0 ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  if ( ! OnlySize ) {
    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
  }
  return(ErrStat) ;
}

int
C_MAP_UnpackInitOutput( float * ReKiBuf,  
                 double * DbKiBuf, 
                 int * IntKiBuf,   
                 MAP_InitOutputType_t *OutData, char * ErrMsg )
{
  int ErrStat ;
  int Re_BufSz2 = 0 ;
  int Db_BufSz2 = 0 ;
  int Int_BufSz2 = 0 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int Re_CurrSz = 0 ;
  int Db_CurrSz = 0 ;
  int Int_CurrSz = 0 ;
  int one        = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes, if any
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}

int
C_MAP_PackContState( float * ReKiBuf,  int * Re_BufSz ,
                 double * DbKiBuf, int * Db_BufSz ,
                 int * IntKiBuf,   int * Int_BufSz ,
                 MAP_ContinuousStateType_t *InData, char * ErrMsg, int *SizeOnly )
{
  int ErrStat ;
  int OnlySize ;
  int Re_BufSz2 ;
  int Db_BufSz2 ;
  int Int_BufSz2 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int one         = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes and subtypes, if any

  OnlySize = *SizeOnly ;

  *Re_BufSz = 0 ;
  *Db_BufSz = 0 ;
  *Int_BufSz = 0 ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  *Db_BufSz   += 1  ; // dummy
  if ( ! OnlySize ) {
    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
    DbKiBuf[Db_Xferred++] = InData->dummy ;
  }
  return(ErrStat) ;
}

int
C_MAP_UnpackContState( float * ReKiBuf,  
                 double * DbKiBuf, 
                 int * IntKiBuf,   
                 MAP_ContinuousStateType_t *OutData, char * ErrMsg )
{
  int ErrStat ;
  int Re_BufSz2 = 0 ;
  int Db_BufSz2 = 0 ;
  int Int_BufSz2 = 0 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int Re_CurrSz = 0 ;
  int Db_CurrSz = 0 ;
  int Int_CurrSz = 0 ;
  int one        = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes, if any
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  OutData->dummy = DbKiBuf [ Db_Xferred ] ; 
  Db_Xferred   = Db_Xferred   + 1 ; 
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}

int
C_MAP_PackDiscState( float * ReKiBuf,  int * Re_BufSz ,
                 double * DbKiBuf, int * Db_BufSz ,
                 int * IntKiBuf,   int * Int_BufSz ,
                 MAP_DiscreteStateType_t *InData, char * ErrMsg, int *SizeOnly )
{
  int ErrStat ;
  int OnlySize ;
  int Re_BufSz2 ;
  int Db_BufSz2 ;
  int Int_BufSz2 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int one         = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes and subtypes, if any

  OnlySize = *SizeOnly ;

  *Re_BufSz = 0 ;
  *Db_BufSz = 0 ;
  *Int_BufSz = 0 ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  *Db_BufSz   += 1  ; // dummy
  if ( ! OnlySize ) {
    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
    DbKiBuf[Db_Xferred++] = InData->dummy ;
  }
  return(ErrStat) ;
}

int
C_MAP_UnpackDiscState( float * ReKiBuf,  
                 double * DbKiBuf, 
                 int * IntKiBuf,   
                 MAP_DiscreteStateType_t *OutData, char * ErrMsg )
{
  int ErrStat ;
  int Re_BufSz2 = 0 ;
  int Db_BufSz2 = 0 ;
  int Int_BufSz2 = 0 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int Re_CurrSz = 0 ;
  int Db_CurrSz = 0 ;
  int Int_CurrSz = 0 ;
  int one        = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes, if any
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  OutData->dummy = DbKiBuf [ Db_Xferred ] ; 
  Db_Xferred   = Db_Xferred   + 1 ; 
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}

int
C_MAP_PackOtherState( float * ReKiBuf,  int * Re_BufSz ,
                 double * DbKiBuf, int * Db_BufSz ,
                 int * IntKiBuf,   int * Int_BufSz ,
                 MAP_OtherStateType_t *InData, char * ErrMsg, int *SizeOnly )
{
  int ErrStat ;
  int OnlySize ;
  int Re_BufSz2 ;
  int Db_BufSz2 ;
  int Int_BufSz2 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int one         = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes and subtypes, if any

  OnlySize = *SizeOnly ;

  *Re_BufSz = 0 ;
  *Db_BufSz = 0 ;
  *Int_BufSz = 0 ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  *Db_BufSz   += InData->FX_Len ; // FX 
  *Db_BufSz   += InData->FY_Len ; // FY 
  *Db_BufSz   += InData->FZ_Len ; // FZ 
  *Int_BufSz  += InData->u_index_Len ; // u_index 
  *Int_BufSz  += InData->p_index_Len ; // p_index 
  *Int_BufSz  += InData->x_index_Len ; // x_index 
  *Int_BufSz  += InData->xd_index_Len ; // xd_index 
  *Int_BufSz  += InData->z_index_Len ; // z_index 
  *Int_BufSz  += InData->y_index_Len ; // y_index 
  *Int_BufSz  += InData->o_index_Len ; // o_index 
  if ( ! OnlySize ) {
    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
    for ( i = 0 ; i < InData->FX_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->FX[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->FY_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->FY[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->FZ_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->FZ[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->u_index_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->u_index[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->p_index_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->p_index[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->x_index_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->x_index[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->xd_index_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->xd_index[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->z_index_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->z_index[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->y_index_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->y_index[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->o_index_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->o_index[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
  }
  return(ErrStat) ;
}

int
C_MAP_UnpackOtherState( float * ReKiBuf,  
                 double * DbKiBuf, 
                 int * IntKiBuf,   
                 MAP_OtherStateType_t *OutData, char * ErrMsg )
{
  int ErrStat ;
  int Re_BufSz2 = 0 ;
  int Db_BufSz2 = 0 ;
  int Int_BufSz2 = 0 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int Re_CurrSz = 0 ;
  int Db_CurrSz = 0 ;
  int Int_CurrSz = 0 ;
  int one        = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes, if any
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  if ( OutData->FX != NULL ) {
    memcpy( OutData->FX,&(DbKiBuf[ Db_Xferred ]),OutData->FX_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->FX_Len ; 
  }
  if ( OutData->FY != NULL ) {
    memcpy( OutData->FY,&(DbKiBuf[ Db_Xferred ]),OutData->FY_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->FY_Len ; 
  }
  if ( OutData->FZ != NULL ) {
    memcpy( OutData->FZ,&(DbKiBuf[ Db_Xferred ]),OutData->FZ_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->FZ_Len ; 
  }
  if ( OutData->u_index != NULL ) {
    memcpy( OutData->u_index,&(DbKiBuf[ Db_Xferred ]),OutData->u_index_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->u_index_Len ; 
  }
  if ( OutData->p_index != NULL ) {
    memcpy( OutData->p_index,&(DbKiBuf[ Db_Xferred ]),OutData->p_index_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->p_index_Len ; 
  }
  if ( OutData->x_index != NULL ) {
    memcpy( OutData->x_index,&(DbKiBuf[ Db_Xferred ]),OutData->x_index_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->x_index_Len ; 
  }
  if ( OutData->xd_index != NULL ) {
    memcpy( OutData->xd_index,&(DbKiBuf[ Db_Xferred ]),OutData->xd_index_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->xd_index_Len ; 
  }
  if ( OutData->z_index != NULL ) {
    memcpy( OutData->z_index,&(DbKiBuf[ Db_Xferred ]),OutData->z_index_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->z_index_Len ; 
  }
  if ( OutData->y_index != NULL ) {
    memcpy( OutData->y_index,&(DbKiBuf[ Db_Xferred ]),OutData->y_index_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->y_index_Len ; 
  }
  if ( OutData->o_index != NULL ) {
    memcpy( OutData->o_index,&(DbKiBuf[ Db_Xferred ]),OutData->o_index_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->o_index_Len ; 
  }
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}
CALL void __stdcall map_f2c_otherstate_fx ( double * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->FX != NULL ) free( dst->FX ) ;
  dst->FX = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->FX ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_fx_ ( MAP_OtherStateType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->FX ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
CALL void map_f2c_otherstate_fy ( double * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->FY != NULL ) free( dst->FY ) ;
  dst->FY = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->FY ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_fy_ ( MAP_OtherStateType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->FY ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
CALL void map_f2c_otherstate_fz_ ( double * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->FZ != NULL ) free( dst->FZ ) ;
  dst->FZ = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->FZ ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_fz_ ( MAP_OtherStateType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->FZ ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_flags_index_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->FLAGS_index != NULL ) free( dst->FLAGS_index ) ;
  dst->FLAGS_index = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->FLAGS_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_flags_index_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->FLAGS_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_plot_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->PLOT_flag != NULL ) free( dst->PLOT_flag ) ;
  dst->PLOT_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->PLOT_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_plot_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->PLOT_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_x_pos_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->X_POS_flag != NULL ) free( dst->X_POS_flag ) ;
  dst->X_POS_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->X_POS_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_x_pos_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->X_POS_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_y_pos_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->Y_POS_flag != NULL ) free( dst->Y_POS_flag ) ;
  dst->Y_POS_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->Y_POS_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_y_pos_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->Y_POS_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_z_pos_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->Z_POS_flag != NULL ) free( dst->Z_POS_flag ) ;
  dst->Z_POS_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->Z_POS_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_z_pos_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->Z_POS_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_x_force_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->X_FORCE_flag != NULL ) free( dst->X_FORCE_flag ) ;
  dst->X_FORCE_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->X_FORCE_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_x_force_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->X_FORCE_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_y_force_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->Y_FORCE_flag != NULL ) free( dst->Y_FORCE_flag ) ;
  dst->Y_FORCE_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->Y_FORCE_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_y_force_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->Y_FORCE_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_z_force_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->Z_FORCE_flag != NULL ) free( dst->Z_FORCE_flag ) ;
  dst->Z_FORCE_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->Z_FORCE_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_z_force_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->Z_FORCE_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_line_tension_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->LINE_TENSION_flag != NULL ) free( dst->LINE_TENSION_flag ) ;
  dst->LINE_TENSION_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->LINE_TENSION_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_line_tension_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->LINE_TENSION_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_omit_contact_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->OMIT_CONTACT_flag != NULL ) free( dst->OMIT_CONTACT_flag ) ;
  dst->OMIT_CONTACT_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->OMIT_CONTACT_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_omit_contact_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->OMIT_CONTACT_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_lay_length_flag_ ( bool * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  if ( dst->LAY_LENGTH_flag != NULL ) free( dst->LAY_LENGTH_flag ) ;
  dst->LAY_LENGTH_flag = (bool *)malloc(( *n1)*sizeof(bool)) ;
  p = (bool *)dst->LAY_LENGTH_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_lay_length_flag_ ( MAP_OtherStateType_t * src, bool * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  bool *p ;
  p = (bool *)src->LAY_LENGTH_flag ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_u_index_ ( int * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  if ( dst->u_index != NULL ) free( dst->u_index ) ;
  dst->u_index = (int *)malloc(( *n1)*sizeof(int)) ;
  p = (int *)dst->u_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_u_index_ ( MAP_OtherStateType_t * src, int * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  p = (int *)src->u_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_p_index_ ( int * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  if ( dst->p_index != NULL ) free( dst->p_index ) ;
  dst->p_index = (int *)malloc(( *n1)*sizeof(int)) ;
  p = (int *)dst->p_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_p_index_ ( MAP_OtherStateType_t * src, int * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  p = (int *)src->p_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_x_index_ ( int * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  if ( dst->x_index != NULL ) free( dst->x_index ) ;
  dst->x_index = (int *)malloc(( *n1)*sizeof(int)) ;
  p = (int *)dst->x_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_x_index_ ( MAP_OtherStateType_t * src, int * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  p = (int *)src->x_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_xd_index_ ( int * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  if ( dst->xd_index != NULL ) free( dst->xd_index ) ;
  dst->xd_index = (int *)malloc(( *n1)*sizeof(int)) ;
  p = (int *)dst->xd_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_xd_index_ ( MAP_OtherStateType_t * src, int * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  p = (int *)src->xd_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_z_index_ ( int * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  if ( dst->z_index != NULL ) free( dst->z_index ) ;
  dst->z_index = (int *)malloc(( *n1)*sizeof(int)) ;
  p = (int *)dst->z_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_z_index_ ( MAP_OtherStateType_t * src, int * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  p = (int *)src->z_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_y_index_ ( int * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  if ( dst->y_index != NULL ) free( dst->y_index ) ;
  dst->y_index = (int *)malloc(( *n1)*sizeof(int)) ;
  p = (int *)dst->y_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_y_index_ ( MAP_OtherStateType_t * src, int * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  p = (int *)src->y_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_otherstate_o_index_ ( int * src, MAP_OtherStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  if ( dst->o_index != NULL ) free( dst->o_index ) ;
  dst->o_index = (int *)malloc(( *n1)*sizeof(int)) ;
  p = (int *)dst->o_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_otherstate_o_index_ ( MAP_OtherStateType_t * src, int * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  int *p ;
  p = (int *)src->o_index ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}

int
C_MAP_PackConstrState( float * ReKiBuf,  int * Re_BufSz ,
                 double * DbKiBuf, int * Db_BufSz ,
                 int * IntKiBuf,   int * Int_BufSz ,
                 MAP_ConstraintStateType_t *InData, char * ErrMsg, int *SizeOnly )
{
  int ErrStat ;
  int OnlySize ;
  int Re_BufSz2 ;
  int Db_BufSz2 ;
  int Int_BufSz2 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int one         = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes and subtypes, if any

  OnlySize = *SizeOnly ;

  *Re_BufSz = 0 ;
  *Db_BufSz = 0 ;
  *Int_BufSz = 0 ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  *Db_BufSz   += InData->X_Len ; // X 
  *Db_BufSz   += InData->Y_Len ; // Y 
  *Db_BufSz   += InData->Z_Len ; // Z 
  *Db_BufSz   += InData->H_Len ; // H 
  *Db_BufSz   += InData->V_Len ; // V 
  if ( ! OnlySize ) {
    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
    for ( i = 0 ; i < InData->X_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->X[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->Y_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Y[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->Z_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Z[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->H_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->H[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->V_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->V[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
  }
  return(ErrStat) ;
}

int
C_MAP_UnpackConstrState( float * ReKiBuf,  
                 double * DbKiBuf, 
                 int * IntKiBuf,   
                 MAP_ConstraintStateType_t *OutData, char * ErrMsg )
{
  int ErrStat ;
  int Re_BufSz2 = 0 ;
  int Db_BufSz2 = 0 ;
  int Int_BufSz2 = 0 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int Re_CurrSz = 0 ;
  int Db_CurrSz = 0 ;
  int Int_CurrSz = 0 ;
  int one        = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes, if any
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  if ( OutData->X != NULL ) {
    memcpy( OutData->X,&(DbKiBuf[ Db_Xferred ]),OutData->X_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->X_Len ; 
  }
  if ( OutData->Y != NULL ) {
    memcpy( OutData->Y,&(DbKiBuf[ Db_Xferred ]),OutData->Y_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->Y_Len ; 
  }
  if ( OutData->Z != NULL ) {
    memcpy( OutData->Z,&(DbKiBuf[ Db_Xferred ]),OutData->Z_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->Z_Len ; 
  }
  if ( OutData->H != NULL ) {
    memcpy( OutData->H,&(DbKiBuf[ Db_Xferred ]),OutData->H_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->H_Len ; 
  }
  if ( OutData->V != NULL ) {
    memcpy( OutData->V,&(DbKiBuf[ Db_Xferred ]),OutData->V_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->V_Len ; 
  }
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}
void map_f2c_constrstate_x_ ( double * src, MAP_ConstraintStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->X != NULL ) free( dst->X ) ;
  dst->X = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->X ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_constrstate_x_ ( MAP_ConstraintStateType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->X ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_constrstate_y_ ( double * src, MAP_ConstraintStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->Y != NULL ) free( dst->Y ) ;
  dst->Y = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->Y ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_constrstate_y_ ( MAP_ConstraintStateType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->Y ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_constrstate_z_ ( double * src, MAP_ConstraintStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->Z != NULL ) free( dst->Z ) ;
  dst->Z = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->Z ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_constrstate_z_ ( MAP_ConstraintStateType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->Z ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_constrstate_h_ ( double * src, MAP_ConstraintStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->H != NULL ) free( dst->H ) ;
  dst->H = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->H ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_constrstate_h_ ( MAP_ConstraintStateType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->H ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_constrstate_v_ ( double * src, MAP_ConstraintStateType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->V != NULL ) free( dst->V ) ;
  dst->V = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->V ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_constrstate_v_ ( MAP_ConstraintStateType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->V ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}

int
C_MAP_PackParam( float * ReKiBuf,  int * Re_BufSz ,
                 double * DbKiBuf, int * Db_BufSz ,
                 int * IntKiBuf,   int * Int_BufSz ,
                 MAP_ParameterType_t *InData, char * ErrMsg, int *SizeOnly )
{
  int ErrStat ;
  int OnlySize ;
  int Re_BufSz2 ;
  int Db_BufSz2 ;
  int Int_BufSz2 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int one         = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes and subtypes, if any

  OnlySize = *SizeOnly ;

  *Re_BufSz = 0 ;
  *Db_BufSz = 0 ;
  *Int_BufSz = 0 ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  *Db_BufSz   += InData->Diam_Len ; // Diam 
  *Db_BufSz   += InData->MassDenInAir_Len ; // MassDenInAir 
  *Db_BufSz   += InData->EA_Len ; // EA 
  *Db_BufSz   += InData->CB_Len ; // CB 
  *Db_BufSz   += InData->Lu_Len ; // Lu 
  *Db_BufSz   += InData->X_Len ; // X 
  *Db_BufSz   += InData->Y_Len ; // Y 
  *Db_BufSz   += InData->Z_Len ; // Z 
  *Db_BufSz   += InData->FX_Len ; // FX 
  *Db_BufSz   += InData->FY_Len ; // FY 
  *Db_BufSz   += InData->FZ_Len ; // FZ 
  *Db_BufSz   += InData->M_Len ; // M 
  *Db_BufSz   += InData->B_Len ; // B 
  if ( ! OnlySize ) {
    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
    for ( i = 0 ; i < InData->Diam_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Diam[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->MassDenInAir_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->MassDenInAir[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->EA_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->EA[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->CB_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->CB[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->Lu_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Lu[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->X_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->X[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->Y_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Y[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->Z_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Z[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->FX_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->FX[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->FY_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->FY[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->FZ_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->FZ[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->M_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->M[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->B_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->B[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
  }
  return(ErrStat) ;
}

int
C_MAP_UnpackParam( float * ReKiBuf,  
                 double * DbKiBuf, 
                 int * IntKiBuf,   
                 MAP_ParameterType_t *OutData, char * ErrMsg )
{
  int ErrStat ;
  int Re_BufSz2 = 0 ;
  int Db_BufSz2 = 0 ;
  int Int_BufSz2 = 0 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int Re_CurrSz = 0 ;
  int Db_CurrSz = 0 ;
  int Int_CurrSz = 0 ;
  int one        = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes, if any
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  if ( OutData->Diam != NULL ) {
    memcpy( OutData->Diam,&(DbKiBuf[ Db_Xferred ]),OutData->Diam_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->Diam_Len ; 
  }
  if ( OutData->MassDenInAir != NULL ) {
    memcpy( OutData->MassDenInAir,&(DbKiBuf[ Db_Xferred ]),OutData->MassDenInAir_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->MassDenInAir_Len ; 
  }
  if ( OutData->EA != NULL ) {
    memcpy( OutData->EA,&(DbKiBuf[ Db_Xferred ]),OutData->EA_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->EA_Len ; 
  }
  if ( OutData->CB != NULL ) {
    memcpy( OutData->CB,&(DbKiBuf[ Db_Xferred ]),OutData->CB_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->CB_Len ; 
  }
  if ( OutData->Lu != NULL ) {
    memcpy( OutData->Lu,&(DbKiBuf[ Db_Xferred ]),OutData->Lu_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->Lu_Len ; 
  }
  if ( OutData->X != NULL ) {
    memcpy( OutData->X,&(DbKiBuf[ Db_Xferred ]),OutData->X_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->X_Len ; 
  }
  if ( OutData->Y != NULL ) {
    memcpy( OutData->Y,&(DbKiBuf[ Db_Xferred ]),OutData->Y_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->Y_Len ; 
  }
  if ( OutData->Z != NULL ) {
    memcpy( OutData->Z,&(DbKiBuf[ Db_Xferred ]),OutData->Z_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->Z_Len ; 
  }
  if ( OutData->FX != NULL ) {
    memcpy( OutData->FX,&(DbKiBuf[ Db_Xferred ]),OutData->FX_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->FX_Len ; 
  }
  if ( OutData->FY != NULL ) {
    memcpy( OutData->FY,&(DbKiBuf[ Db_Xferred ]),OutData->FY_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->FY_Len ; 
  }
  if ( OutData->FZ != NULL ) {
    memcpy( OutData->FZ,&(DbKiBuf[ Db_Xferred ]),OutData->FZ_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->FZ_Len ; 
  }
  if ( OutData->M != NULL ) {
    memcpy( OutData->M,&(DbKiBuf[ Db_Xferred ]),OutData->M_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->M_Len ; 
  }
  if ( OutData->B != NULL ) {
    memcpy( OutData->B,&(DbKiBuf[ Db_Xferred ]),OutData->B_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->B_Len ; 
  }
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}
void map_f2c_param_diam_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->Diam != NULL ) free( dst->Diam ) ;
  dst->Diam = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->Diam ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_diam_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->Diam ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_massdeninair_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->MassDenInAir != NULL ) free( dst->MassDenInAir ) ;
  dst->MassDenInAir = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->MassDenInAir ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_massdeninair_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->MassDenInAir ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_ea_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->EA != NULL ) free( dst->EA ) ;
  dst->EA = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->EA ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_ea_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->EA ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_cb_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->CB != NULL ) free( dst->CB ) ;
  dst->CB = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->CB ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_cb_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->CB ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_lu_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->Lu != NULL ) free( dst->Lu ) ;
  dst->Lu = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->Lu ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_lu_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->Lu ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_x_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->X != NULL ) free( dst->X ) ;
  dst->X = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->X ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_x_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->X ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_y_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->Y != NULL ) free( dst->Y ) ;
  dst->Y = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->Y ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_y_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->Y ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_z_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->Z != NULL ) free( dst->Z ) ;
  dst->Z = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->Z ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_z_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->Z ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_fx_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->FX != NULL ) free( dst->FX ) ;
  dst->FX = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->FX ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_fx_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->FX ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_fy_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->FY != NULL ) free( dst->FY ) ;
  dst->FY = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->FY ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_fy_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->FY ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_fz_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->FZ != NULL ) free( dst->FZ ) ;
  dst->FZ = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->FZ ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_fz_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->FZ ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_m_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->M != NULL ) free( dst->M ) ;
  dst->M = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->M ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_m_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->M ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_param_b_ ( double * src, MAP_ParameterType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->B != NULL ) free( dst->B ) ;
  dst->B = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->B ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_param_b_ ( MAP_ParameterType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->B ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}

int
C_MAP_PackInput( float * ReKiBuf,  int * Re_BufSz ,
                 double * DbKiBuf, int * Db_BufSz ,
                 int * IntKiBuf,   int * Int_BufSz ,
                 MAP_InputType_t *InData, char * ErrMsg, int *SizeOnly )
{
  int ErrStat ;
  int OnlySize ;
  int Re_BufSz2 ;
  int Db_BufSz2 ;
  int Int_BufSz2 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int one         = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes and subtypes, if any
  float   * Re_Position_Buf ;
  double  * Db_Position_Buf ;
  int     * Int_Position_Buf ;

  OnlySize = *SizeOnly ;

  *Re_BufSz = 0 ;
  *Db_BufSz = 0 ;
  *Int_BufSz = 0 ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  *Db_BufSz   += InData->X_Len ; // X 
  *Db_BufSz   += InData->Y_Len ; // Y 
  *Db_BufSz   += InData->Z_Len ; // Z 
  if ( ! OnlySize ) {
    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
    for ( i = 0 ; i < InData->X_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->X[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->Y_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Y[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->Z_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->Z[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
  }
  return(ErrStat) ;
}

int
C_MAP_UnpackInput( float * ReKiBuf,  
                 double * DbKiBuf, 
                 int * IntKiBuf,   
                 MAP_InputType_t *OutData, char * ErrMsg )
{
  int ErrStat ;
  int Re_BufSz2 = 0 ;
  int Db_BufSz2 = 0 ;
  int Int_BufSz2 = 0 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int Re_CurrSz = 0 ;
  int Db_CurrSz = 0 ;
  int Int_CurrSz = 0 ;
  int one        = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes, if any
  float   * Re_Position_Buf ;
  double  * Db_Position_Buf ;
  int     * Int_Position_Buf ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  if ( OutData->X != NULL ) {
    memcpy( OutData->X,&(DbKiBuf[ Db_Xferred ]),OutData->X_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->X_Len ; 
  }
  if ( OutData->Y != NULL ) {
    memcpy( OutData->Y,&(DbKiBuf[ Db_Xferred ]),OutData->Y_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->Y_Len ; 
  }
  if ( OutData->Z != NULL ) {
    memcpy( OutData->Z,&(DbKiBuf[ Db_Xferred ]),OutData->Z_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->Z_Len ; 
  }
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}
void map_f2c_input_x_ ( double * src, MAP_InputType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->X != NULL ) free( dst->X ) ;
  dst->X = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->X ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_input_x_ ( MAP_InputType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->X ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_input_y_ ( double * src, MAP_InputType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->Y != NULL ) free( dst->Y ) ;
  dst->Y = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->Y ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_input_y_ ( MAP_InputType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->Y ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_input_z_ ( double * src, MAP_InputType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->Z != NULL ) free( dst->Z ) ;
  dst->Z = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->Z ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_input_z_ ( MAP_InputType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->Z ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}

int
C_MAP_PackOutput( float * ReKiBuf,  int * Re_BufSz ,
                 double * DbKiBuf, int * Db_BufSz ,
                 int * IntKiBuf,   int * Int_BufSz ,
                 MAP_OutputType_t *InData, char * ErrMsg, int *SizeOnly )
{
  int ErrStat ;
  int OnlySize ;
  int Re_BufSz2 ;
  int Db_BufSz2 ;
  int Int_BufSz2 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int one         = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes and subtypes, if any
  float   * Re_Force_Buf ;
  double  * Db_Force_Buf ;
  int     * Int_Force_Buf ;

  OnlySize = *SizeOnly ;

  *Re_BufSz = 0 ;
  *Db_BufSz = 0 ;
  *Int_BufSz = 0 ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  *Db_BufSz   += InData->FX_Len ; // FX 
  *Db_BufSz   += InData->FY_Len ; // FY 
  *Db_BufSz   += InData->FZ_Len ; // FZ 
  if ( ! OnlySize ) {
    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;
    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;
    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;
    for ( i = 0 ; i < InData->FX_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->FX[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->FY_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->FY[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
    for ( i = 0 ; i < InData->FZ_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(DbKiBuf[Db_Xferred+i]), &(InData->FZ[i]), sizeof(double)) ;
      Db_Xferred++ ;
    }
  }
  return(ErrStat) ;
}

int
C_MAP_UnpackOutput( float * ReKiBuf,  
                 double * DbKiBuf, 
                 int * IntKiBuf,   
                 MAP_OutputType_t *OutData, char * ErrMsg )
{
  int ErrStat ;
  int Re_BufSz2 = 0 ;
  int Db_BufSz2 = 0 ;
  int Int_BufSz2 = 0 ;
  int Re_Xferred = 0 ;
  int Db_Xferred = 0 ;
  int Int_Xferred = 0 ;
  int Re_CurrSz = 0 ;
  int Db_CurrSz = 0 ;
  int Int_CurrSz = 0 ;
  int one        = 1 ;
  int i,i1,i2,i3,i4,i5 ;
 // buffers to store meshes, if any
  float   * Re_Force_Buf ;
  double  * Db_Force_Buf ;
  int     * Int_Force_Buf ;
  ReKiBuf = NULL ;
  DbKiBuf = NULL ;
  IntKiBuf = NULL ;
  if ( OutData->FX != NULL ) {
    memcpy( OutData->FX,&(DbKiBuf[ Db_Xferred ]),OutData->FX_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->FX_Len ; 
  }
  if ( OutData->FY != NULL ) {
    memcpy( OutData->FY,&(DbKiBuf[ Db_Xferred ]),OutData->FY_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->FY_Len ; 
  }
  if ( OutData->FZ != NULL ) {
    memcpy( OutData->FZ,&(DbKiBuf[ Db_Xferred ]),OutData->FZ_Len) ;
    Db_Xferred   = Db_Xferred   + OutData->FZ_Len ; 
  }
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}
void map_f2c_output_fx_ ( double * src, MAP_OutputType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->FX != NULL ) free( dst->FX ) ;
  dst->FX = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->FX ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_output_fx_ ( MAP_OutputType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->FX ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_output_fy_ ( double * src, MAP_OutputType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->FY != NULL ) free( dst->FY ) ;
  dst->FY = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->FY ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_output_fy_ ( MAP_OutputType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->FY ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
void map_f2c_output_fz_ ( double * src, MAP_OutputType_t * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  if ( dst->FZ != NULL ) free( dst->FZ ) ;
  dst->FZ = (double *)malloc(( *n1)*sizeof(double)) ;
  p = (double *)dst->FZ ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    *p++ = src[i1] ;
  }
}
void map_c2f_output_fz_ ( MAP_OutputType_t * src, double * dst , int *n1) {
  int i1,i2,i3,i4,i5 ;
  double *p ;
  p = (double *)src->FZ ;
  for (i1 = 0 ; i1 < *n1 ; i1++ )
  {
    dst[i1] = *p++ ;
  }
}
//!ENDOFREGISTRYGENERATEDFILE
