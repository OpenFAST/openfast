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

// gen preamble
  {
    char ** p ;
    for ( p = FAST_preamble ; *p ; p++ ) { fprintf( fp, *p, ModName->name ) ; }
  }

// generate each derived data type
fprintf(stderr,">> %s \n",ModName->name ) ;
  for ( q = ModName->module_ddt_list ; q ; q = q->next )
  {
fprintf(stderr,">>   %s \n",q->name ) ;
    fprintf(fp,"q->name %s\n",q->name) ;
    for ( r = q->fields ; r ; r = r->next )
    { 
fprintf(stderr,">>     %s \n",r->name ) ;
    }
  }

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



