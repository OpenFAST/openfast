#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#define rindex(X,Y) strrchr(X,Y)
#define index(X,Y) strchr(X,Y)
#define bzero(X,Y) memset(X,0,Y)
#else
#  include <strings.h>
#endif

#include "registry.h"
#include "protos.h"
#include "data.h"

int
init_modname_table()
{
  ModNames = NULL ;
  return(0) ;
}

int
init_dim_table()
{
  Dim = NULL ;
  return(0) ;
}

node_t * 
new_node ( int kind )
{ node_t *p ; 
  p = (node_t *)malloc(sizeof(node_t)) ; 
  bzero(p,sizeof(node_t)); 
  p->node_kind = kind ; 

  p->fields          = NULL;
  p->params          = NULL;
  p->type            = NULL;
  p->module          = NULL;
  p->module_ddt_list = NULL;
  p->next            = NULL;
  //p->coord_end_param = NULL;
  strcpy(p->dim_param_name, "");
  p->dim_param = 0;
  p->type_type = 0;
  p->max_ndims = 0;
  p->containsPtr = 0;
  p->ndims = 0;
  p->deferred = 0;
  p->usefrom = 0;
  p->is_interface_type = 0;
  strcpy(p->name, "");
  strcpy(p->mapsto, "");
  strcpy(p->nickname, "");
  strcpy(p->descrip, "");
  strcpy(p->units, "");

  return (p) ; }

int
add_node_to_end ( node_t * node , node_t ** list )
{
  node_t * p ;
  if ( *list == NULL ) 
    { *list = node ; }
  else
  {
    for ( p = *list ; p->next != NULL ; p = p->next ) ;
    p->next = node ;
  }
  return(0) ;
}

int
add_node_to_beg ( node_t * node , node_t ** list )
{
  node_t * p ;
  if ( *list == NULL )
  { 
    *list = node ;
    (*list)->next = NULL ;
  }
  else
  {
//fprintf(stderr,"   add_node_to_beg: node %s to existing list. CH %s CN %08x\n", node->name,(*list)->name,(*list)->next) ;
//if ( (*list)->next ) fprintf(stderr,"    CN name %s\n",(*list)->next->name ) ;
    p = (*list) ;
    *list = node ;
    (*list)->next = p ;
  }
  return(0) ;
}


#if 0
int
add_node_to_end_4d ( node_t * node , node_t ** list )
{
  node_t * p ;
  if ( *list == NULL ) 
    { *list = node ; }
  else
  {
    for ( p = *list ; p->next4d != NULL ; p = p->next4d ) ;
    p->next4d = node ;
  }
  return(0) ;
}
#endif

#if 1

void
show_nodelist( node_t * p )
{
  show_nodelist1( p , 0 ) ;
}

void
show_nodelist1( node_t * p , int indent )
{
  if ( p == NULL ) return;
  show_node1( p, indent) ;
  show_nodelist1( p->next, indent ) ;
}

int
show_node( node_t * p )
{
  return(show_node1(p,0)) ;
}

int
show_node1( node_t * p, int indent )
{
  char spaces[] = "                           " ;
  char tmp[25] , t1[25] , t2[25] ;
  char * x, *ca, *ld, *ss, *se, *sg ;
  char *nodekind ;
  int nl ;
  int i ;

  if ( p == NULL ) return(1) ;
  strcpy(tmp, spaces) ;
  if ( indent >= 0 && indent < 20 ) tmp[indent] = '\0' ;

// this doesn't make much sense any more, ever since node_kind was 
// changed to a bit mask
  nodekind = "" ;
  if      ( p->node_kind & FIELD   ) nodekind = "FIELD" ;
  else if ( p->node_kind & MODNAME ) nodekind = "MODNAME" ;
  else if ( p->node_kind & TYPE    ) nodekind = "TYPE" ;

  switch ( p->node_kind )
  {
  case MODNAME :
    fprintf(stderr,"%s%s : %s nickname %s\n",tmp,nodekind,p->name,p->nickname) ;
    show_nodelist1(p->module_ddt_list, indent+1) ;
    break ;
  case FIELD   :
    fprintf(stderr,"%s%s : %10s ndims %1d\n",tmp,nodekind,p->name, p->ndims) ;
    for ( i = 0 ; i < p->ndims ; i++ )
    {
      sg = "" ;
      ca = "" ;
      switch ( p->dims[i]->coord_axis ) {
        case COORD_C : ca = "C" ; break ;
      }
      switch ( p->dims[i]->len_defined_how ) {
	case DOMAIN_STANDARD : ld = "STANDARD" ;  ss = "" ; se = "" ; break ;
	case CONSTANT        : ld = "CONSTANT" ;  sprintf(t1,"%d",p->dims[i]->coord_start) ; ss = t1 ;
						  sprintf(t2,"%d",p->dims[i]->coord_end  ) ; se = t2 ;
						  break ;
      }
      fprintf(stderr,"      dim %0d: {%s} %2s%s %10s %10s %10s\n",i,p->dims[i]->dim_name,ca,sg,ld,ss,se) ;
    }
    nl = 0 ;
    if ( strlen( p->use     ) > 0 ) {
       nl = 1 ; fprintf(stderr,"      use: %s",p->use) ;
    }
    if ( strlen( p->descrip ) > 0 ) { nl = 1 ; fprintf(stderr,"  descrip: %s",p->descrip) ;    }
    if ( nl == 1 ) fprintf(stderr,"\n") ;
    show_node1( p->type, indent+1 ) ;
    break ;
  case TYPE  :
    x = "derived" ;
    if ( p->type_type == SIMPLE ) x = "simple" ;
    fprintf(stderr,"%sTYPE : %10s %s ndims %1d\n",tmp,p->name,x, p->ndims) ;
    show_nodelist1( p->fields, indent+1 ) ;
    break ;
  case DIM   :
    break ;
  default :
    break ;
  }
  return(0) ;
}
#endif

int
set_mark ( int val , node_t * lst )
{
  node_t * p ;
  if ( lst == NULL ) return(0) ;
  for ( p = lst ; p != NULL ; p = p->next )
  {
    p->mark = val ;
    set_mark( val , p->fields ) ;
  }
  return(0) ;
}

#if 0
int
set_mark_4d ( int val , node_t * lst )
{
  node_t * p ;
  if ( lst == NULL ) return(0) ;
  for ( p = lst ; p != NULL ; p = p->next4d )
  {
    p->mark = val ;
    set_mark( val , p->fields ) ;
    set_mark( val , p->members ) ;
  }
  return(0) ;
}
#endif

