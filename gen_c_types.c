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
gen_c_types( FILE * fp , node_t * ModName )
{
  node_t * p, * q, * r ;
  int i ;
  char nonick[NAMELEN] ;

  if ( strlen(ModName->nickname) > 0 ) {
// generate each derived data type
    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      if ( q->usefrom == 0 ) {
        fprintf(fp,"  struct %s {\n",q->mapsto) ;
        for ( r = q->fields ; r ; r = r->next )
        {
          if ( r->type != NULL ) {
            if ( r->type->type_type == DERIVED ) {
              fprintf(fp,"    struct %s ",r->type->mapsto ) ;
            } else {
              char tmp[NAMELEN] ; tmp[0] = '\0' ;
              if ( q->mapsto) remove_nickname( ModName->nickname, make_lower_temp(q->mapsto) , tmp ) ;
              fprintf(fp,"    %s ",C_type( r->type->mapsto) ) ;
            }
            if ( r->ndims > 0 )
            {
              if ( r->dims[0]->deferred )     // if one dim is deferred they all have to be; see check in type.c
              {
                fprintf(fp,"* %s ; \n", r->name) ;
#if 0
                for ( i = 0 ; i < r->ndims ; i++ )
                {
                  fprintf(fp,":") ;
                  if ( i < r->ndims-1 ) fprintf(fp,",") ;
                }
                fprintf(fp,"), ALLOCATABLE ") ;
#endif
              } else {
                fprintf(fp," %s") ;
                for ( i = 0 ; i < r->ndims ; i++ )
                {
                  fprintf(fp,"[((%d)-(%d)+1)]",r->dims[i]->coord_start,r->dims[i]->coord_end) ;
                }
                fprintf(fp,"; \n") ;
              }
            }
            if ( r->ndims == 0 ) {
              if ( strlen(r->inival) > 0 ) {
                fprintf(fp," %s = %s ;\n", r->name, r->inival ) ;
              } else {
                fprintf(fp," %s ;\n",r->name) ;
              }
            }
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

#if 0
    fprintf(fp,"CONTAINS\n") ;
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

        gen_copy( fp, ModName, ddtname, ddtnamelong ) ;
        gen_destroy( fp, ModName, ddtname, ddtnamelong ) ;
        gen_pack( fp, ModName, ddtname, ddtnamelong ) ;
        gen_unpack( fp, ModName, ddtname, ddtnamelong ) ;
      }
    }

    gen_modname_pack( fp, ModName ) ;
    gen_modname_unpack( fp, ModName ) ;
//    gen_rk4( fp, ModName ) ;
    gen_ExtrapInterp( fp, ModName, "Input", "inputtype" ) ;
    gen_ExtrapInterp( fp, ModName, "Output", "outputtype" ) ;

    fprintf(fp,"END MODULE %s_Types\n",ModName->name ) ;
#endif
  }
}
