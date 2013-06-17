#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
# include <strings.h>
#endif

#include "protos.h"
#include "registry.h"
#include "data.h"

int
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
    return(1) ;
  }
fprintf(fp,"C_%s_Pack%s( float * ReKiBuf,  int * Re_BufSz ,\n",ModName->nickname,nonick) ;
fprintf(fp,"                 double * DbKiBuf, int * Db_BufSz ,\n") ;
fprintf(fp,"                 int * IntKiBuf,   int * Int_BufSz ,\n") ;
fprintf(fp,"                 struct MAP_InputType *InData, char * ErrMsg, int * SizeOnly )\n") ;
fprintf(fp,"{\n") ;
fprintf(fp,"  int OnlySize ;\n") ;
fprintf(fp,"  int Re_BufSz2 ;\n") ;
fprintf(fp,"  int Db_BufSz2 ;\n") ;
fprintf(fp,"  int Int_BufSz2 ;\n") ;
fprintf(fp,"  int Re_Xferred = 0 ;\n") ;
fprintf(fp,"  int Db_Xferred = 0 ;\n") ;
fprintf(fp,"  int Int_Xferred = 0 ;\n") ;
fprintf(fp,"  int i,i1,i2,i3,i4,i5 ;\n") ;
fprintf(fp," // buffers to store meshes and subtypes, if any\n") ;

  for ( r = q->fields ; r ; r = r->next )
  {
    if ( r->type == NULL ) {
      fprintf(stderr,"Registry warning generating %_Pack%s: %s has no type.\n",ModName->nickname,nonick,r->name) ;
      return ; // EARLY RETURN
    } else {
      if ( !strcmp( r->type->name, "meshtype" ) || (r->type->type_type == DERIVED && ! r->type->usefrom ) ) {
  fprintf(fp,"  float   * Re_%s_Buf\n",r->name) ;
  fprintf(fp,"  double  * Db_%s_Buf\n",r->name) ;
  fprintf(fp,"  int     * Int_%s_Buf\n",r->name) ;
      }
    }
  }

fprintf(fp,"\n") ;
fprintf(fp,"  OnlySize = *SizeOnly ;\n") ;
fprintf(fp,"\n") ;
fprintf(fp,"  *Re_BufSz = 0 ;\n") ;
fprintf(fp,"  *Db_BufSz = 0 ;\n") ;
fprintf(fp,"  *Int_BufSz = 0 ;\n") ;
  frst = 1 ;
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( r->type->type_type == DERIVED && ! r->type->usefrom ) {
      char nonick2[NAMELEN] ;
      remove_nickname(ModName->nickname,r->type->name,nonick2) ;
  fprintf(fp,"  ErrStat = C_%s_Pack%s( Re_%s_Buf,  int &Re_BufSz2  ,\n",
                        ModName->nickname,fast_interface_type_shortname(nonick2), r->name) ;
  fprintf(fp,"                         Db_%s_Buf,  int &Db_BufSz2  ,\n",r->name ) ;
  fprintf(fp,"                         Int_%s_Buf, int &Int_BufSz2 , InData%%%s%s, ErrMsg, 1 ) ! %s \n",
                         r->name, r->name, r->name, dimstr(r->ndims), r->name ) ;
  fprintf(fp,"  *Re_BufSz += Re_BufSz2  ) ! %s\n",r->name ) ;
  fprintf(fp,"  *Db_BufSz += Db_BufSz2  ) ! %s\n",r->name ) ;
  fprintf(fp,"  *Int_BufSz += Int_BufSz2  ) ! %s\n",r->name ) ;
  fprintf(fp,"  if ( Re_%s_Buf != NULL)  { free(Re_%s_Buf) ; Re_%s_Buf = NULL ;} \n",r->name, r->name, r->name) ;
  fprintf(fp,"  if ( Db_%s_Buf != NULL)  { free(Db_%s_Buf) ; Db_%s_Buf = NULL  ;}\n",r->name, r->name, r->name) ;
  fprintf(fp,"  if ( Int_%s_Buf != NULL) { free(Int_%s_Buf) ; Int_%s_Buf = NULL ;} \n",r->name, r->name, r->name) ;
    } else if ( r->ndims == 0 ) {  // scalars
      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  *Re_BufSz   += 1  ! %s\n",r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  *Db_BufSz   += 1  ! %s\n",r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  *Int_BufSz  += 1  ! %s\n",r->name ) ;
      }
    } else { // r->ndims > 0
      if ( r->dims[0]->deferred ) {
        if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  *Re_BufSz   += InData->%s.size() )  ! %s \n", r->name , r->name ) ;
        }
        else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  *Db_BufSz   += InData->%s.size() )  ! %s \n", r->name , r->name ) ;
        }
        else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  *Int_BufSz  += InData->%s.size() )  ! %s \n", r->name , r->name ) ;
        }
      } else {
      }
    }
  }

   // Allocate buffers
  fprintf(fp,"  if ( *Re_BufSz > 0  ) ReKiBuf  = (float  *)malloc(*Re_BufSz*sizeof(float) ) ;\n") ;
  fprintf(fp,"  if ( *Db_BufSz > 0  ) DbKiBuf  = (double *)malloc(*Db_BufSz*sizeof(double) ) ;\n") ;
  fprintf(fp,"  if ( *Int_BufSz > 0 ) IntKiBuf = (int    *)malloc(*Int_BufSz*sizeof(int) ) ;\n") ;


   // Pack data
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( r->type->type_type == DERIVED && ! r->type->usefrom ) {
      char nonick2[NAMELEN] ;
      remove_nickname(ModName->nickname,r->type->name,nonick2) ;
  fprintf(fp,"  ErrStat = C_%s_Pack%s( Re_%s_Buf,  int &Re_BufSz2  ,\n",
                        ModName->nickname,fast_interface_type_shortname(nonick2), r->name) ;
  fprintf(fp,"                         Db_%s_Buf,  int &Db_BufSz2  ,\n",r->name ) ;
  fprintf(fp,"                         Int_%s_Buf, int &Int_BufSz2 , InData%%%s%s, ErrMsg, 1 ) ! %s \n",
                            r->name, r->name, r->name, dimstr(r->ndims), r->name ) ;

  fprintf(fp,"  if ( Re_%s_Buf != NULL ) {\n",r->name) ;
  fprintf(fp,"    memcpy( &ReKiBuf[Re_Xferred], Re_%s_Buf, Re_BufSz2*sizeof(float) ) ;\n",r->name) ;
  fprintf(fp,"    Re_Xferred += Re_BufSz2 )\n") ; 
  fprintf(fp,"  }\n" ) ;
  fprintf(fp,"  if ( Db_%s_Buf != NULL ) {\n",r->name) ;
  fprintf(fp,"    memcpy( &DbKiBuf[Db_Xferred], Db_%s_Buf, Db_BufSz2*sizeof(double) ) ;\n",r->name) ;
  fprintf(fp,"    Db_Xferred += Db_BufSz2 )\n") ; 
  fprintf(fp,"  }\n" ) ;
  fprintf(fp,"  if ( Int_%s_Buf != NULL ) {\n",r->name) ;
  fprintf(fp,"    memcpy( &IntKiBuf[Int_Xferred], Int_%s_Buf, Int_BufSz2*sizeof(int) ) ;\n",r->name) ;
  fprintf(fp,"    Int_Xferred += Int_BufSz2 )\n") ; 
  fprintf(fp,"  }\n" ) ;


  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    IF ( .NOT. OnlySize ) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_%s_Buf)-1 ) = Db_%s_Buf\n",r->name,r->name) ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE(Db_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    IF ( .NOT. OnlySize ) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_%s_Buf)-1 ) = Int_%s_Buf\n",r->name,r->name) ;
  fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE(Int_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF( ALLOCATED(Re_%s_Buf) )  DEALLOCATE(Re_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ALLOCATED(Db_%s_Buf) )  DEALLOCATE(Db_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ALLOCATED(Int_%s_Buf) ) DEALLOCATE(Int_%s_Buf)\n",r->name, r->name) ;

    } else  {
      char * indent ;
//      sprintf(tmp2,"SIZE(InData%%%s)",r->name) ;
      sprintf(tmp2,"InData->%s.size())",r->name) ;
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
      indent = "" ;
      if ( !strcmp( r->type->mapsto, "REAL(ReKi)") ||
           !strcmp( r->type->mapsto, "REAL(DbKi)") ||
           !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
        if ( r->ndims > 0 && has_deferred_dim( r, 0 )) {
          indent = "  " ;
        }
        if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"%sfor ( i = 0 ; i < InData->%s.size() ; i++ ) {\n",indent, r->name ) ;
  fprintf(fp,"%s  if ( !OnlySize ) memcpy( &(ReKiBuf[Re_Xferred+i]), &(InData->%s[i]), sizeof(float)) ;\n",
              indent,r->name  ) ;
  fprintf(fp,"%s  Re_Xferred++ ;\n") ;
  fprintf(fp,"%s}\n") ;
        }
        if ( r->ndims > 0 && has_deferred_dim( r, 0 )) {
//  fprintf(fp,"  ENDIF\n") ;
        }
      }
    }
  }


fprintf(fp,"}\n") ;
}



int
gen_c_module( FILE * fp , node_t * ModName )
{
  node_t * p, * q, * r ;
  int i ;
  char nonick[NAMELEN], star ;

  if ( strlen(ModName->nickname) > 0 ) {
// generate each derived data type
    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      if ( q->usefrom == 0 ) {
        fprintf(fp,  "  struct %s_C {\n",q->mapsto) ;
        if ( sw_embed_class_ptr ) {
          fprintf(fp,"    %s * class ; // user must define a class named %s in C++ code\n",q->mapsto,q->mapsto) ;
          fprintf(fp,"    int *index ;\n") ;
          fprintf(fp,"    int indexLen ;\n") ;
        }
        for ( r = q->fields ; r ; r = r->next )
        {
          if ( r->type != NULL ) {
            star = ' ' ;
            if ( r->ndims > 0 ) {
              if ( r->dims[0]->deferred ) star = '*' ;
            }
            if ( r->type->type_type == DERIVED ) {
              fprintf(fp,"    struct %s %c%s",r->type->mapsto,star,r->name ) ;
            } else {
              char tmp[NAMELEN] ; tmp[0] = '\0' ;
              if ( q->mapsto) remove_nickname( ModName->nickname, make_lower_temp(q->mapsto) , tmp ) ;
              if ( r->ndims > 0 && r->dims[0]->deferred ) {
//                fprintf(fp,"    std::vector<%s> %s ",C_type( r->type->mapsto), r->name ) ;
                fprintf(fp,"    %s * %s ; ",C_type( r->type->mapsto), r->name ) ;
                fprintf(fp,"    int %sLen ",r->name ) ;
              } else {
                fprintf(fp,"    %s %s",C_type( r->type->mapsto ),r->name ) ;
              }
            }
            for ( i = 0 ; i < r->ndims ; i++ )
            {
              if ( ! r->dims[0]->deferred ) 
                fprintf(fp,"[((%d)-(%d)+1)]",r->dims[i]->coord_end,r->dims[i]->coord_start) ;
            }
            if ( r->ndims == 0 ) {
              if ( strlen(r->inival) > 0 ) {
                fprintf(fp," = %s ",  r->inival ) ;
              }
            }
            fprintf(fp,"; \n") ;
          }
        }
        fprintf(fp,"  } ;\n") ;
      }
    }


    fprintf(fp,"  struct %s_UserData {\n", ModName->nickname) ;
    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      remove_nickname(ModName->nickname,q->name,nonick) ;
fprintf(stderr,"%s %d\n",nonick,is_a_fast_interface_type(nonick)) ;
      if ( is_a_fast_interface_type(nonick) ) {
        fprintf(fp,"    struct %-30s %s_%s ;\n", q->mapsto, ModName->nickname, fast_interface_type_shortname(nonick) ) ;
      }
    }
    fprintf(fp,"  } ;\n") ;

    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      if ( q->usefrom == 0 ) {

        char * ddtname, * ddtnamelong, nonick[NAMELEN] ;
        ddtname = q->name ;

        remove_nickname(ModName->nickname,ddtname,nonick) ;

//fprintf(stderr,">> %s %s %s \n",ModName->name, ddtname, nonick) ;

        if ( is_a_fast_interface_type( nonick ) ) {
          ddtnamelong = nonick ;
          ddtname = fast_interface_type_shortname( nonick ) ;
        } else {
          ddtnamelong = ddtname ;
        }

  //      gen_copy( fp, ModName, ddtname, ddtnamelong ) ;
  //      gen_destroy( fp, ModName, ddtname, ddtnamelong ) ;
        gen_c_pack( fp, ModName, ddtname, ddtnamelong ) ;
  //      gen_unpack( fp, ModName, ddtname, ddtnamelong ) ;
      }
    }
  }
}
