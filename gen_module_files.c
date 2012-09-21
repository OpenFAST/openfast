#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
# include <strings.h>
#endif

#include "protos.h"
#include "registry.h"
#include "data.h"

#include "FAST_preamble.h"

int
gen_module( FILE * fp , const node_t * ModName )
{
  node_t * q, * r ;
  int i ;

// gen preamble
  {
    char ** p ;
    for ( p = FAST_preamble ; *p ; p++ ) { fprintf( fp, *p, ModName->name ) ; }
  }

// generate each derived data type
  for ( q = ModName->module_ddt_list ; q ; q = q->next )
  {
    fprintf(fp,"  TYPE, PUBLIC :: %s\n",q->mapsto) ;
    for ( r = q->fields ; r ; r = r->next )
    { 
      fprintf(fp,"    %s ",r->type->mapsto ) ;
      if ( r->ndims > 0 )
      {
        if ( r->dims[0]->deferred )     // if one dim is they all have to be; see check in type.c
        {
          fprintf(fp,", DIMENSION(") ;
          for ( i = 0 ; i < r->ndims ; i++ )
          {
            fprintf(fp,":") ;
            if ( i < r->ndims-1 ) fprintf(fp,",") ;
          }
          fprintf(fp,"), ALLOCATABLE ") ;
        } else {
          fprintf(fp,", DIMENSION(") ;
          for ( i = 0 ; i < r->ndims ; i++ )
          {
            fprintf(fp,"%d:%d",r->dims[i]->coord_start,r->dims[i]->coord_end) ;
            if ( i < r->ndims-1 ) fprintf(fp,",") ;
          }
          fprintf(fp,") ") ;
        }
      }
      fprintf(fp," :: %s \n",r->name) ;

    }
    fprintf(fp,"  END TYPE PUBLIC %s\n",q->mapsto) ;
  }
  fprintf(fp,"CONTAINS\n") ;
  fprintf(fp,"END MODULE %s_Types\n",ModName->name ) ;

}


int
gen_module_files ( char * dirname )
{
  FILE * fp ;
  char  fname[NAMELEN] ;
  char * fn ;

  node_t * p ;
  
  for ( p = ModNames ; p ; p = p->next )
  {
    if ( strlen(dirname) > 0 ) 
      { sprintf(fname,"%s/%s",dirname,p->name) ; }
    else                       
      { sprintf(fname,"%s",p->name) ; }
    if ((fp = fopen( fname , "w" )) == NULL ) return(1) ;
    print_warning(fp,fname) ;
    gen_module ( fp , p ) ;
    close_the_file( fp ) ;
  }
  return(0) ;
}



