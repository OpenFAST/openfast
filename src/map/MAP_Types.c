//!STARTOFREGISTRYGENERATEDFILE './MAP_Types.c'
//
//! WARNING This file is generated automatically by the FAST registry
//! Do not edit.  Your changes to this file will be lost.
//!
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MAP_Types.h"


#ifdef _WIN32 //define something for Windows (32-bit)
  #include "stdbool.h"
  #define CALL __declspec( dllexport )
#elif _WIN64 //define something for Windows (64-bit)
  #include "stdbool.h"
  #define CALL __declspec( dllexport ) 
#else
  #include <stdbool.h>
  #define CALL 
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

CALL void MAP_F2C_InitOutput_WriteOutputHdr_C ( MAP_InitOutputType_t *type, char *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->WriteOutputHdr[i] = arr[i];
}

CALL void MAP_F2C_InitOutput_WriteOutputUnt_C ( MAP_InitOutputType_t *type, char *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->WriteOutputUnt[i] = arr[i];
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

CALL void MAP_F2C_OtherState_FX_C ( MAP_OtherStateType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->FX[i] = arr[i];
}

CALL void MAP_F2C_OtherState_FY_C ( MAP_OtherStateType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->FY[i] = arr[i];
}

CALL void MAP_F2C_OtherState_FZ_C ( MAP_OtherStateType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->FZ[i] = arr[i];
}

CALL void MAP_F2C_OtherState_u_index_C ( MAP_OtherStateType_t *type, int *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->u_index[i] = arr[i];
}

CALL void MAP_F2C_OtherState_p_index_C ( MAP_OtherStateType_t *type, int *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->p_index[i] = arr[i];
}

CALL void MAP_F2C_OtherState_x_index_C ( MAP_OtherStateType_t *type, int *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->x_index[i] = arr[i];
}

CALL void MAP_F2C_OtherState_xd_index_C ( MAP_OtherStateType_t *type, int *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->xd_index[i] = arr[i];
}

CALL void MAP_F2C_OtherState_z_index_C ( MAP_OtherStateType_t *type, int *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->z_index[i] = arr[i];
}

CALL void MAP_F2C_OtherState_y_index_C ( MAP_OtherStateType_t *type, int *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->y_index[i] = arr[i];
}

CALL void MAP_F2C_OtherState_o_index_C ( MAP_OtherStateType_t *type, int *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->o_index[i] = arr[i];
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

CALL void MAP_F2C_ConstrState_X_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->X[i] = arr[i];
}

CALL void MAP_F2C_ConstrState_Y_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->Y[i] = arr[i];
}

CALL void MAP_F2C_ConstrState_Z_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->Z[i] = arr[i];
}

CALL void MAP_F2C_ConstrState_H_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->H[i] = arr[i];
}

CALL void MAP_F2C_ConstrState_V_C ( MAP_ConstraintStateType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->V[i] = arr[i];
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

CALL void MAP_F2C_Param_Diam_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->Diam[i] = arr[i];
}

CALL void MAP_F2C_Param_MassDenInAir_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->MassDenInAir[i] = arr[i];
}

CALL void MAP_F2C_Param_EA_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->EA[i] = arr[i];
}

CALL void MAP_F2C_Param_CB_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->CB[i] = arr[i];
}

CALL void MAP_F2C_Param_Lu_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->Lu[i] = arr[i];
}

CALL void MAP_F2C_Param_X_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->X[i] = arr[i];
}

CALL void MAP_F2C_Param_Y_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->Y[i] = arr[i];
}

CALL void MAP_F2C_Param_Z_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->Z[i] = arr[i];
}

CALL void MAP_F2C_Param_FX_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->FX[i] = arr[i];
}

CALL void MAP_F2C_Param_FY_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->FY[i] = arr[i];
}

CALL void MAP_F2C_Param_FZ_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->FZ[i] = arr[i];
}

CALL void MAP_F2C_Param_M_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->M[i] = arr[i];
}

CALL void MAP_F2C_Param_B_C ( MAP_ParameterType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->B[i] = arr[i];
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
  float   * Re_PtFairleadDisplacement_Buf ;
  double  * Db_PtFairleadDisplacement_Buf ;
  int     * Int_PtFairleadDisplacement_Buf ;

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
  float   * Re_PtFairleadDisplacement_Buf ;
  double  * Db_PtFairleadDisplacement_Buf ;
  int     * Int_PtFairleadDisplacement_Buf ;
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

CALL void MAP_F2C_Input_X_C ( MAP_InputType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->X[i] = arr[i];
}

CALL void MAP_F2C_Input_Y_C ( MAP_InputType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->Y[i] = arr[i];
}

CALL void MAP_F2C_Input_Z_C ( MAP_InputType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->Z[i] = arr[i];
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
  float   * Re_PtFairleadLoad_Buf ;
  double  * Db_PtFairleadLoad_Buf ;
  int     * Int_PtFairleadLoad_Buf ;

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
  *Re_BufSz   += InData->writeOutput_Len ; // writeOutput 
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
    for ( i = 0 ; i < InData->writeOutput_Len ; i++ ) {
      if ( !OnlySize ) memcpy( &(ReKiBuf[Re_Xferred+i]), &(InData->writeOutput[i]), sizeof(float)) ;
      Re_Xferred++ ;
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
  float   * Re_PtFairleadLoad_Buf ;
  double  * Db_PtFairleadLoad_Buf ;
  int     * Int_PtFairleadLoad_Buf ;
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
  if ( OutData->writeOutput != NULL ) {
    memcpy( OutData->writeOutput,&(ReKiBuf[ Re_Xferred ]),OutData->writeOutput_Len) ;
    Re_Xferred   = Re_Xferred   + OutData->writeOutput_Len ; 
  }
  if ( ReKiBuf != NULL )  free(ReKiBuf) ;
  if ( DbKiBuf != NULL )  free(DbKiBuf) ;
  if ( IntKiBuf != NULL ) free(IntKiBuf) ;
  return(ErrStat) ;
}

CALL void MAP_F2C_Output_FX_C ( MAP_OutputType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->FX[i] = arr[i];
}

CALL void MAP_F2C_Output_FY_C ( MAP_OutputType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->FY[i] = arr[i];
}

CALL void MAP_F2C_Output_FZ_C ( MAP_OutputType_t *type, double *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->FZ[i] = arr[i];
}

CALL void MAP_F2C_Output_writeOutput_C ( MAP_OutputType_t *type, float *arr, int len )
{
  int i = 0;
  for( i=0 ; i<=len-1 ; i++ ) type->writeOutput[i] = arr[i];
}
//!ENDOFREGISTRYGENERATEDFILE
