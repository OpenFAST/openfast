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
              if ( r->ndims > 0 && r->dims[0]->deferred ) {
                fprintf(fp,"    std::vector<%s> %s ;\n",C_type( r->type->mapsto), r->name ) ;
              } else {
                fprintf(fp,"    %s ",C_type( r->type->mapsto ) ) ;
                for ( i = 0 ; i < r->ndims ; i++ )
                {
                  fprintf(fp,"[((%d)-(%d)+1)]",r->dims[i]->coord_start,r->dims[i]->coord_end) ;
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

  }
}
