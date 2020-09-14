#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifndef _WIN32
# include <strings.h>
#endif

#include "protos.h"
#include "registry.h"
#include "data.h"


#if 0
void
gen_c_unpack( FILE * fp, const node_t * ModName, char * inout, char * inoutlong )
{
  char tmp[NAMELEN], tmp2[NAMELEN], tmp3[NAMELEN], tmp4[NAMELEN], addnick[NAMELEN], nonick[NAMELEN] ;
  node_t *q, * r ;
  int d, idim, frst ;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
  sprintf(tmp,"%s",addnick) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Unpack%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
    return;//(1) ;
  }

fprintf(fp,"\nint\n") ;
fprintf(fp,"C_%s_Unpack%s( float * ReKiBuf,  \n",ModName->nickname,nonick) ;
fprintf(fp,"                 double * DbKiBuf, \n") ;
fprintf(fp,"                 int * IntKiBuf,   \n") ;
fprintf(fp,"                 %s_t *OutData, char * ErrMsg )\n", addnick) ;
fprintf(fp,"{\n") ;
fprintf(fp,"  int ErrStat = 0;\n") ;
fprintf(fp,"  int Re_BufSz2 = 0 ;\n") ;
fprintf(fp,"  int Db_BufSz2 = 0 ;\n") ;
fprintf(fp,"  int Int_BufSz2 = 0 ;\n") ;
fprintf(fp,"  int Re_Xferred = 0 ;\n") ;
fprintf(fp,"  int Db_Xferred = 0 ;\n") ;
fprintf(fp,"  int Int_Xferred = 0 ;\n") ;
fprintf(fp,"  int Re_CurrSz = 0 ;\n") ;
fprintf(fp,"  int Db_CurrSz = 0 ;\n") ;
fprintf(fp,"  int Int_CurrSz = 0 ;\n") ;
fprintf(fp,"  int one        = 1 ;\n") ;
fprintf(fp,"  int i,i1,i2,i3,i4,i5 ;\n") ;

  fprintf(fp," // buffers to store meshes, if any\n") ;
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( r->type == NULL ) {
      fprintf(stderr,"Registry warning generating %s_Unpack%s: %s has no type.\n",ModName->nickname,nonick,r->name) ;
      return ; // EARLY RETURN
    } else {
      if ( !strcmp( r->type->name, "meshtype" ) || (r->type->type_type == DERIVED && ! r->type->usefrom ) ) {
  fprintf(fp,"  float   * Re_%s_Buf ;\n",r->name) ;
  fprintf(fp,"  double  * Db_%s_Buf ;\n",r->name) ;
  fprintf(fp,"  int     * Int_%s_Buf ;\n",r->name) ;
      }
    }
  }
fprintf(fp,"  ReKiBuf = NULL ;\n") ;
fprintf(fp,"  DbKiBuf = NULL ;\n") ;
fprintf(fp,"  IntKiBuf = NULL ;\n") ;

   // Unpack data
  frst = 1 ;
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( r->type->type_type == DERIVED && ! r->type->usefrom && strcmp(make_lower_temp(r->type->mapsto),"meshtype") ) {
      char nonick2[NAMELEN] ;
      remove_nickname(ModName->nickname,r->type->name,nonick2) ;
  fprintf(fp," // first call %s_Pack%s to get correctly sized buffers for unpacking\n",
                        ModName->nickname,fast_interface_type_shortname(nonick2)) ;
  fprintf(fp,"  ErrStat = C_%s_Pack%s( Re_%s_Buf, &Re_BufSz2, Db_%s_Buf, &Db_BufSz2, Int_%s_Buf, &Int_BufSz2, &(OutData->%s%s), ErrMsg, &one ) ; // %s \n",
                    ModName->nickname,fast_interface_type_shortname(nonick2), r->name, r->name, r->name, r->name, dimstr_c(r->ndims),r->name ) ;

  fprintf(fp,"  if ( Re_%s_Buf != NULL ) {\n",r->name) ;
  fprintf(fp,"    memcpy( Re_%s_Buf, &(ReKiBuf[ Re_Xferred] ), Re_BufSz2 ) ;\n",r->name ) ;
  fprintf(fp,"    Re_Xferred += Re_BufSz2  ; // %s \n",r->name) ;
  fprintf(fp,"  }\n" ) ;
  fprintf(fp,"  if ( Db_%s_Buf != NULL ) {\n",r->name) ;
  fprintf(fp,"    memcpy( Db_%s_Buf, &(DbKiBuf[ Db_Xferred] ), Db_BufSz2 ) ;\n",r->name ) ;
  fprintf(fp,"    Db_Xferred += Db_BufSz2  ; // %s \n",r->name) ;
  fprintf(fp,"  }\n" ) ;
  fprintf(fp,"  if ( Int_%s_Buf != NULL ) {\n",r->name) ;
  fprintf(fp,"    memcpy( Int_%s_Buf, &(IntKiBuf[ Int_Xferred] ), Int_BufSz2 ) ;\n",r->name ) ;
  fprintf(fp,"    Int_Xferred += Int_BufSz2 ; // %s \n",r->name) ;
  fprintf(fp,"  }\n" ) ;
  fprintf(fp,"  ErrStat = C_%s_Unpack%s( Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, &(OutData->%s%s), ErrMsg ) ; // %s \n",
                  ModName->nickname,fast_interface_type_shortname(nonick2), r->name, r->name, r->name, r->name,
                  dimstr(r->ndims),
                  r->name ) ;
//  fprintf(fp,"  if ( Re_%s_Buf != NULL)  { free(Re_%s_Buf) ; Re_%s_Buf = NULL ;} \n",r->name, r->name, r->name) ;
//  fprintf(fp,"  if ( Db_%s_Buf != NULL)  { free(Db_%s_Buf) ; Db_%s_Buf = NULL  ;}\n",r->name, r->name, r->name) ;
//  fprintf(fp,"  if ( Int_%s_Buf != NULL) { free(Int_%s_Buf) ; Int_%s_Buf = NULL ;} \n",r->name, r->name, r->name) ;

    } else  {
      char * indent, * ty ;
      char arrayname[NAMELEN], tmp[NAMELEN], tmp2[NAMELEN] ;

      sprintf(arrayname,"OutData%%%s",r->name) ;
      sprintf(tmp2,"SIZE(OutData%%%s)",r->name) ;
      if      ( r->ndims==0 ) { strcpy(tmp3,"") ; }
      else if ( r->ndims==1 ) { strcpy(tmp3,"") ; }
      else if ( r->ndims==2 ) { sprintf(tmp3,"(1:(%s),1)",tmp2) ; }
      else if ( r->ndims==3 ) { sprintf(tmp3,"(1:(%s),1,1)",tmp2) ; }
      else if ( r->ndims==4 ) { sprintf(tmp3,"(1:(%s),1,1,1)",tmp2) ; }
      else if ( r->ndims==5 ) { sprintf(tmp3,"(1:(%s),1,1,1,1)",tmp2) ; }
      else                    { fprintf(stderr,"Registry WARNING: too many dimensions for %s\n",r->name) ; }
      indent = "" ;
      if ( !strcmp( r->type->mapsto, "REAL(ReKi)") ||
           !strcmp( r->type->mapsto, "REAL(DbKi)") ||
           !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
        if ( r->ndims > 0 && has_deferred_dim( r, 0 )) {
  fprintf(fp,"  if ( OutData->%s != NULL ) {\n", r->name ) ;
          indent = "  " ;
        }

        if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) ty = "Re" ;
        if      ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) ty = "Db" ;
        if      ( !strcmp( r->type->mapsto, "REAL(IntKi)")    ) ty = "Int" ;

        if ( r->ndims > 0 ) {
          if ( has_deferred_dim( r, 0 ) ) {
  fprintf(fp,"%s  memcpy( OutData->%s,&(%sKiBuf[ %s_Xferred ]),OutData->%s_Len) ;\n",indent,r->name,ty,ty,r->name) ;
  fprintf(fp,"%s  %s_Xferred   = %s_Xferred   + OutData->%s_Len ; \n",indent,ty,ty,r->name ) ;
          } else {
            int i ;
            strcpy(tmp2,"") ;
            for ( i = 0 ; i < r->ndims ; i++ )
            {
              sprintf(tmp,"((%d)-(%d)+1)",r->dims[i]->coord_end,r->dims[i]->coord_start) ;
              strcat(tmp2,tmp) ;
              if ( i < r->ndims-1 ) strcat(tmp2,"*") ;
            }
  fprintf(fp,"%s  memcpy( OutData->%s,&(%sKiBuf[ %s_Xferred ]),(%s)*sizeof(%s)) ;\n",
            indent,r->name,ty,ty,tmp2,C_type(r->type->mapsto)) ;
  fprintf(fp,"%s  %s_Xferred   = %s_Xferred   + (%s)*sizeof(%s) ; \n",
            indent,ty,ty,tmp2,C_type(r->type->mapsto) ) ;
          }
        } else {
  fprintf(fp,"%s  OutData->%s = %sKiBuf [ %s_Xferred ] ; \n",indent,r->name,ty,ty) ;
  fprintf(fp,"%s  %s_Xferred   = %s_Xferred   + 1 ; \n",indent,ty,ty ) ;
        }

        if ( r->ndims > 0 && has_deferred_dim( r, 0 )) {
  fprintf(fp,"  }\n" ) ;
        }

      }
    }
  }
  fprintf(fp,"  if ( ReKiBuf != NULL )  free(ReKiBuf) ;\n") ;
  fprintf(fp,"  if ( DbKiBuf != NULL )  free(DbKiBuf) ;\n") ;
  fprintf(fp,"  if ( IntKiBuf != NULL ) free(IntKiBuf) ;\n") ;
  fprintf(fp,"  return(ErrStat) ;\n") ;
  fprintf(fp,"}\n") ;
  return;//(0) ;
}

void
gen_c_pack( FILE * fp, const node_t * ModName, char * inout, char *inoutlong )
{
  char tmp[NAMELEN], tmp2[NAMELEN], tmp3[NAMELEN], addnick[NAMELEN], nonick[NAMELEN] ;
  node_t *q, * r ;
  int frst, d ;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
  sprintf(tmp,"%s",addnick) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Pack%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
    return;//(1) ;
  }
fprintf(fp,"\nint\n") ;
fprintf(fp,"C_%s_Pack%s( float * ReKiBuf,  int * Re_BufSz ,\n",ModName->nickname,nonick) ;
fprintf(fp,"                 double * DbKiBuf, int * Db_BufSz ,\n") ;
fprintf(fp,"                 int * IntKiBuf,   int * Int_BufSz ,\n") ;
fprintf(fp,"                 %s_t *InData, char * ErrMsg, int *SizeOnly )\n", addnick) ;
fprintf(fp,"{\n") ;
fprintf(fp,"  int ErrStat = 0;\n") ;
fprintf(fp,"  int OnlySize ;\n") ;
fprintf(fp,"  int Re_BufSz2 ;\n") ;
fprintf(fp,"  int Db_BufSz2 ;\n") ;
fprintf(fp,"  int Int_BufSz2 ;\n") ;
fprintf(fp,"  int Re_Xferred = 0 ;\n") ;
fprintf(fp,"  int Db_Xferred = 0 ;\n") ;
fprintf(fp,"  int Int_Xferred = 0 ;\n") ;
fprintf(fp,"  int one         = 1 ;\n") ;
fprintf(fp,"  int i,i1,i2,i3,i4,i5 ;\n") ;
fprintf(fp," // buffers to store meshes and subtypes, if any\n") ;

  for ( r = q->fields ; r ; r = r->next )
  {
    if ( r->type == NULL ) {
      fprintf(stderr,"Registry warning generating %s_Pack%s: %s has no type.\n",ModName->nickname,nonick,r->name) ;
      return ; // EARLY RETURN
    } else {
      if ( !strcmp( r->type->name, "meshtype" ) || (r->type->type_type == DERIVED && ! r->type->usefrom ) ) {
  fprintf(fp,"  float   * Re_%s_Buf ;\n",r->name) ;
  fprintf(fp,"  double  * Db_%s_Buf ;\n",r->name) ;
  fprintf(fp,"  int     * Int_%s_Buf ;\n",r->name) ;
      }
    }
  }

fprintf(fp,"\n") ;
fprintf(fp,"  OnlySize = *SizeOnly ;\n") ;
fprintf(fp,"\n") ;
fprintf(fp,"  *Re_BufSz = 0 ;\n") ;
fprintf(fp,"  *Db_BufSz = 0 ;\n") ;
fprintf(fp,"  *Int_BufSz = 0 ;\n") ;
fprintf(fp,"  ReKiBuf = NULL ;\n") ;
fprintf(fp,"  DbKiBuf = NULL ;\n") ;
fprintf(fp,"  IntKiBuf = NULL ;\n") ;
  frst = 1 ;
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( r->type->type_type == DERIVED && ! r->type->usefrom && strcmp(make_lower_temp(r->type->mapsto),"meshtype") ) {
      char nonick2[NAMELEN] ;
      remove_nickname(ModName->nickname,r->type->name,nonick2) ;
  fprintf(fp,"  ErrStat = C_%s_Pack%s( Re_%s_Buf,  &Re_BufSz2  ,\n",
                        ModName->nickname,fast_interface_type_shortname(nonick2), r->name) ;
  fprintf(fp,"                         Db_%s_Buf,  &Db_BufSz2  ,\n",r->name ) ;
  fprintf(fp,"                         Int_%s_Buf, &Int_BufSz2 , &(InData->%s%s), ErrMsg, &one ) ; // %s \n",
                         r->name, r->name, dimstr(r->ndims), r->name ) ;
  fprintf(fp,"  *Re_BufSz += Re_BufSz2  ; // %s\n",r->name ) ;
  fprintf(fp,"  *Db_BufSz += Db_BufSz2  ; // %s\n",r->name ) ;
  fprintf(fp,"  *Int_BufSz += Int_BufSz2  ; // %s\n",r->name ) ;
  fprintf(fp,"  if ( Re_%s_Buf != NULL)  { free(Re_%s_Buf) ; Re_%s_Buf = NULL ;} \n",r->name, r->name, r->name) ;
  fprintf(fp,"  if ( Db_%s_Buf != NULL)  { free(Db_%s_Buf) ; Db_%s_Buf = NULL  ;}\n",r->name, r->name, r->name) ;
  fprintf(fp,"  if ( Int_%s_Buf != NULL) { free(Int_%s_Buf) ; Int_%s_Buf = NULL ;} \n",r->name, r->name, r->name) ;
    } else if ( r->ndims == 0 ) {  // scalars
      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  *Re_BufSz   += 1  ; // %s\n",r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  *Db_BufSz   += 1  ; // %s\n",r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  *Int_BufSz  += 1  ; // %s\n",r->name ) ;
      }
    } else { // r->ndims > 0
      if ( r->dims[0]->deferred ) {
        if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  *Re_BufSz   += InData->%s_Len ; // %s \n", r->name , r->name ) ;
        }
        else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  *Db_BufSz   += InData->%s_Len ; // %s \n", r->name , r->name ) ;
        }
        else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  *Int_BufSz  += InData->%s_Len ; // %s \n", r->name , r->name ) ;
        }
      } else {
      }
    }
  }

  fprintf(fp,"  if ( ! OnlySize ) {\n") ;
   // Allocate buffers
  fprintf(fp,"    if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;\n") ;
  fprintf(fp,"    if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;\n") ;
  fprintf(fp,"    if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;\n") ;

   // Pack data
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( r->type->type_type == DERIVED && ! r->type->usefrom && strcmp(make_lower_temp(r->type->mapsto),"meshtype") ) {
      char nonick2[NAMELEN] ;
      remove_nickname(ModName->nickname,r->type->name,nonick2) ;
  fprintf(fp,"    ErrStat = C_%s_Pack%s( Re_%s_Buf,  &Re_BufSz2  ,\n",
                        ModName->nickname,fast_interface_type_shortname(nonick2), r->name) ;
  fprintf(fp,"                           Db_%s_Buf,  &Db_BufSz2  ,\n",r->name ) ;
  fprintf(fp,"                           Int_%s_Buf, &Int_BufSz2 , &(InData->%s%s), ErrMsg, &one ) ; // %s \n",
                            r->name, r->name, dimstr(r->ndims), r->name ) ;

  fprintf(fp,"    if ( Re_%s_Buf != NULL ) {\n",r->name) ;
  fprintf(fp,"      memcpy( &ReKiBuf[Re_Xferred], Re_%s_Buf, Re_BufSz2*sizeof(float) ) ;\n",r->name) ;
  fprintf(fp,"      Re_Xferred += Re_BufSz2 ;\n") ; 
  fprintf(fp,"    }\n" ) ;
  fprintf(fp,"    if ( Db_%s_Buf != NULL ) {\n",r->name) ;
  fprintf(fp,"      memcpy( &DbKiBuf[Db_Xferred], Db_%s_Buf, Db_BufSz2*sizeof(double) ) ;\n",r->name) ;
  fprintf(fp,"      Db_Xferred += Db_BufSz2 ;\n") ; 
  fprintf(fp,"    }\n" ) ;
  fprintf(fp,"    if ( Int_%s_Buf != NULL ) {\n",r->name) ;
  fprintf(fp,"      memcpy( &IntKiBuf[Int_Xferred], Int_%s_Buf, Int_BufSz2*sizeof(int) ) ;\n",r->name) ;
  fprintf(fp,"      Int_Xferred += Int_BufSz2 ;\n") ; 
  fprintf(fp,"    }\n" ) ;
  fprintf(fp,"    if ( Re_%s_Buf != NULL)  { free(Re_%s_Buf) ; Re_%s_Buf = NULL ;} \n",r->name, r->name, r->name) ;
  fprintf(fp,"    if ( Db_%s_Buf != NULL)  { free(Db_%s_Buf) ; Db_%s_Buf = NULL  ;}\n",r->name, r->name, r->name) ;
  fprintf(fp,"    if ( Int_%s_Buf != NULL) { free(Int_%s_Buf) ; Int_%s_Buf = NULL ;} \n",r->name, r->name, r->name) ;

    } else  {
      char * indent, *ty, *cty ;
      sprintf(tmp2,"InData->%s_Len)",r->name) ;
      if        ( r->ndims==0 ) {
        strcpy(tmp3,"") ;
      } else if ( r->ndims==1 ) {
        strcpy(tmp3,"") ;
      } else if ( r->ndims==2 ) {
        sprintf(tmp3,"(1:(%s),1)",tmp2) ;
      } else if ( r->ndims==3 ) {
        sprintf(tmp3,"(1:(%s),1,1)",tmp2) ;
      } else if ( r->ndims==4 ) {
        sprintf(tmp3,"(1:(%s),1,1,1)",tmp2) ;
      } else if ( r->ndims==5 ) {
        sprintf(tmp3,"(1:(%s),1,1,1,1)",tmp2) ;
      } else                    {
        fprintf(stderr,"Registry WARNING: too many dimensions for %s\n",r->name) ;
      }
      indent = "  " ;
      if ( !strcmp( r->type->mapsto, "REAL(ReKi)") ||
           !strcmp( r->type->mapsto, "REAL(DbKi)") ||
           !strcmp( r->type->mapsto, "INTEGER(IntKi)") )
      {
        if      ( !strcmp( r->type->mapsto, "REAL(ReKi)") )  {ty = "Re" ;  cty = "float"  ; }
        else if ( !strcmp( r->type->mapsto, "REAL(DbKi)") )  {ty = "Db" ;  cty = "double" ; }
        else if ( !strcmp( r->type->mapsto, "REAL(IntKi)") ) {ty = "Int" ; cty = "int"    ; }
        indent = "    " ;
        if ( r->ndims > 0 && has_deferred_dim( r, 0 )) {
  fprintf(fp,"%sfor ( i = 0 ; i < InData->%s_Len ; i++ ) {\n",indent, r->name ) ;
  fprintf(fp,"%s  if ( !OnlySize ) memcpy( &(%sKiBuf[%s_Xferred+i]), &(InData->%s[i]), sizeof(%s)) ;\n",
              indent,ty,ty,r->name,cty  ) ;
  fprintf(fp,"%s  %s_Xferred++ ;\n",indent,ty) ;
  fprintf(fp,"%s}\n",indent) ;
        } else if ( r->ndims == 0 ) {
  fprintf(fp,"    %sKiBuf[%s_Xferred++] = InData->%s ;\n",ty,ty,r->name) ;
        }
      }
    }
  }

fprintf(fp,"  }\n") ;
fprintf(fp,"  return(ErrStat) ;\n") ;
fprintf(fp,"}\n") ;
return;//(0) ;
}
#endif


void
gen_c_module( FILE * fph, node_t * ModName )
{
  node_t * q, * r ;
  int i ;
  char nonick[NAMELEN], star ;

  if ( strlen(ModName->nickname) > 0 ) {
// generate each derived data type
    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      if ( q->usefrom == 0 ) {
        if (*q->mapsto) remove_nickname(ModName->nickname, make_lower_temp(q->mapsto), nonick);
        fprintf(fph,  "  typedef struct %s {\n",q->mapsto) ;
        //if (!strcmp(make_lower_temp(nonick), "otherstatetype") !strcmp(make_lower_temp(nonick), "initinputtype")){
           fprintf(fph, "    void * object ;\n");
        //}
        for ( r = q->fields ; r ; r = r->next )
        {
          if ( r->type != NULL ) {
            star = ' ' ;
            if ( r->ndims > 0 ) {
               if ( has_deferred_dim(r, 0) ) star = '*';
            }
            if ( r->type->type_type == DERIVED ) {
              if ( strcmp(make_lower_temp(r->type->mapsto),"meshtype") ) { // do not output mesh types for C code,
                //fprintf(fph,"    struct %s %c%s",r->type->mapsto,star,r->name ) ;
              }
            } else {
              char tmp[NAMELEN] ; tmp[0] = '\0' ;
              if (*q->mapsto) remove_nickname( ModName->nickname, make_lower_temp(q->mapsto) , tmp ) ;
              if (r->ndims > 0 && has_deferred_dim(r, 0)) {
                fprintf(fph,"    %s * %s ; ",C_type( r->type->mapsto), r->name ) ;
                fprintf(fph,"    int %s_Len ;",r->name ) ;
              } else {
                char *p = r->type->mapsto;
                char buf[10];
// bjj: this assumes all character strings are defined with numeric lengths
// It should be modified to allow use of parameters, too. (and parameters defined in the registry should also be defined in the .h file)
                while (*p) { 
                  if (isdigit(*p)) { 
                    long val = strtol(p, &p, 10); 
                    snprintf(buf, 10, "%lu", val);
                  } else { 
                    p++;
                  }


                }    
                if (strcmp(C_type(r->type->mapsto), "char") == 0 ){ // if it's a char we need to add the array size
                   if (r->ndims == 0) 
                      fprintf(fph,"    %s %s[%s] ;",C_type( r->type->mapsto ),r->name,buf ) ;
                } else { // else, it's just a double or int value
                  fprintf(fph,"    %s %s ;",C_type( r->type->mapsto ),r->name ) ;
                }                
              }
            }
            for ( i = 0 ; i < r->ndims ; i++ )
            {
               if (!has_deferred_dim(r, 0) && (strcmp(C_type(r->type->mapsto), "char") || r->ndims == 0)) //  skip this for characters?
                fprintf(fph,"[%d] ;",r->dims[i]->coord_end - r->dims[i]->coord_start +1) ;
            }
            fprintf(fph, "\n");
          }
        }
        fprintf(fph,"  } %s_t ;\n", q->mapsto ) ;
      }
    }


    fprintf(fph,"  typedef struct %s_UserData {\n", ModName->nickname) ;
    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      remove_nickname(ModName->nickname,q->name,nonick) ;
      if ( is_a_fast_interface_type(nonick) ) {
        char temp[NAMELEN] ;
        sprintf(temp, "%s_t", q->mapsto ) ;
        fprintf(fph,"    %-30s %s_%s ;\n", temp, ModName->nickname, fast_interface_type_shortname(nonick) ) ;
      }
    }
    fprintf(fph,"  } %s_t ;\n", ModName->nickname ) ;

  }
}
