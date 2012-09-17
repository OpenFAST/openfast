#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# define rindex(X,Y) strrchr(X,Y)
# define index(X,Y) strchr(X,Y)
#else
# include <strings.h>
#endif

#include "registry.h"
#include "protos.h"
#include "data.h"
#include "sym.h"

/* read in the Registry file and build the internal representation of the registry */

#define MAXTOKENS 1000

/* fields for state entries (note, these get converted to field entries in the
   reg_parse routine; therefore, only TABLE needs to be looked at */
#define TABLE 0

/* fields for field entries (TABLE="typedef" and, with some munging,  TABLE="state") */
#define FIELD_MODNAME   1
#define FIELD_OF        2
#define FIELD_TYPE      3
#define FIELD_SYM       4
#define FIELD_DIMS      5
#define FIELD_CTRL      6
#define FIELD_DNAME     7
#define FIELD_DESCRIP   8
#define FIELD_UNITS     9

#define F_MODNAME  0
#define F_OF       1
#define F_TYPE     2
#define F_SYM      3
#define F_DIMS     4
#define F_CTRL     5
#define F_DNAME    6
#define F_DESCRIP  7
#define F_UNITS    8

/* fields for dimension entries (TABLE="dimspec") */
#define DIM_NAME       1
#define DIM_ORDER      2
#define DIM_SPEC       3

#define INLN_SIZE      8000
#define PARSELINE_SIZE 8000

int isNum( char c )
{
    if ( c < '0' || c > '9' ) return 0; 
    return 1 ;
}

int
pre_parse( char * dir, FILE * infile, FILE * outfile )
{
  /* Decreased size for SOA from 8192 to 8000--double check if necessary, Manish Shrivastava 2010 */
  char inln[INLN_SIZE], parseline[PARSELINE_SIZE], parseline_save[PARSELINE_SIZE] ;
  int found ; 
  char *p, *q ;
  char *tokens[MAXTOKENS], *toktmp[MAXTOKENS], newdims[NAMELEN_LONG], newdims4d[NAMELEN_LONG],newname[NAMELEN_LONG] ;
  int i, ii, len_of_tok ;
  char x, xstr[NAMELEN_LONG] ;
  int is4d, wantstend, wantsbdy ;
  int ifdef_stack_ptr = 0 ;
  int ifdef_stack[100] ;
  int inquote, retval ;

  ifdef_stack[0] = 1 ;
  retval = 0 ;

  parseline[0] = '\0' ;
  while ( fgets ( inln , INLN_SIZE , infile ) != NULL )
  {
/*** preprocessing directives ****/
    /* look for an include statement */
    for ( p = inln ; ( *p == ' ' || *p == '	' ) && *p != '\0' ; p++ ) ;
    if ( !strncmp( make_lower_temp(p) , "include", 7 ) &&  ! ( ifdef_stack_ptr >= 0 && ! ifdef_stack[ifdef_stack_ptr] ) ) {
      FILE *include_fp ;
      char include_file_name[128] ;
      p += 7 ; for ( ; ( *p == ' ' || *p == '	' ) && *p != '\0' ; p++ ) ;
      if ( strlen( p ) > 127 ) { fprintf(stderr,"Registry warning: invalid include file name: %s\n", p ) ; }
      else {
        sprintf( include_file_name , "%s/%s", dir , p ) ;
        if ( (p=index(include_file_name,'\n')) != NULL ) *p = '\0' ;
        fprintf(stderr,"opening %s\n",include_file_name) ;
        if (( include_fp = fopen( include_file_name , "r" )) != NULL ) {

          fprintf(stderr,"including %s\n",include_file_name ) ;
          pre_parse( dir , include_fp , outfile ) ;

          fclose( include_fp ) ;
        } else {
          fprintf(stderr,"Registry warning: cannot open %s. Ignoring.\n", include_file_name ) ;
        } 
      }
    }
    else if ( !strncmp( make_lower_temp(p) , "ifdef", 5 ) ) {
      char value[32] ;
      p += 5 ; for ( ; ( *p == ' ' || *p == '	' ) && *p != '\0' ; p++ ) ;
      strncpy(value, p, 31 ) ; value[31] = '\0' ;
      if ( (p=index(value,'\n')) != NULL ) *p = '\0' ;
      if ( (p=index(value,' ')) != NULL ) *p = '\0' ; if ( (p=index(value,'	')) != NULL ) *p = '\0' ; 
      ifdef_stack_ptr++ ;
      ifdef_stack[ifdef_stack_ptr] = ( sym_get(value) != NULL && ifdef_stack[ifdef_stack_ptr-1] ) ;
      if ( ifdef_stack_ptr >= 100 ) { fprintf(stderr,"Registry fatal: too many nested ifdefs\n") ; exit(1) ; }
      continue ;
    }
    else if ( !strncmp( make_lower_temp(p) , "ifndef", 6 ) ) {
      char value[32] ;
      p += 6 ; for ( ; ( *p == ' ' || *p == '	' ) && *p != '\0' ; p++ ) ;
      strncpy(value, p, 31 ) ; value[31] = '\0' ;
      if ( (p=index(value,'\n')) != NULL ) *p = '\0' ;
      if ( (p=index(value,' ')) != NULL ) *p = '\0' ; if ( (p=index(value,'	')) != NULL ) *p = '\0' ; 
      ifdef_stack_ptr++ ;
      ifdef_stack[ifdef_stack_ptr] = ( sym_get(value) == NULL && ifdef_stack[ifdef_stack_ptr-1] ) ;
      if ( ifdef_stack_ptr >= 100 ) { fprintf(stderr,"Registry fatal: too many nested ifdefs\n") ; exit(1) ; }
      continue ;
    }
    else if ( !strncmp( make_lower_temp(p) , "endif", 5 ) ) {
      ifdef_stack_ptr-- ; 
      if ( ifdef_stack_ptr < 0 ) { fprintf(stderr,"Registry fatal: unmatched endif\n") ; exit(1) ; }
      continue ;
    }
    else if ( !strncmp( make_lower_temp(p) , "define", 6 ) ) {
      char value[32] ;
      p += 6 ; for ( ; ( *p == ' ' || *p == '	' ) && *p != '\0' ; p++ ) ;
      strncpy(value, p, 31 ) ; value[31] = '\0' ;
      if ( (p=index(value,'\n')) != NULL ) *p = '\0' ;
      if ( (p=index(value,' ')) != NULL ) *p = '\0' ; if ( (p=index(value,'	')) != NULL ) *p = '\0' ; 
      sym_add( value ) ;
      continue ;
    }
    if ( ifdef_stack_ptr >= 0 && ! ifdef_stack[ifdef_stack_ptr] ) continue ;
/*** end of preprocessing directives ****/

    strcat( parseline , inln ) ;

    /* allow \ to continue the end of a line */
    if (( p = index( parseline,  '\\'  )) != NULL )
    {
      if ( *(p+1) == '\n' || *(p+1) == '\0' )
      {
        *p = '\0' ;
        continue ;  /* go get another line */
      }
    }
//    make_lower( parseline ) ;

    if (( p = index( parseline , '\n' )) != NULL  ) *p = '\0' ; /* discard newlines */

    /* check line and zap any # characters that are in double quotes */

    for ( p = parseline, inquote = 0 ; *p ; p++ ) {
      if      ( *p == '"' && inquote ) inquote = 0 ;
      else if ( *p == '"' && !inquote ) inquote = 1 ;
      else if ( *p == '#' && inquote ) *p = ' ' ;
      else if ( *p == '#' && !inquote ) { *p = '\0' ; break ; }
    }
    if ( inquote ) { retval=1 ; fprintf(stderr,"Registry error: unbalanced quotes in line:\n%s\n",parseline) ;}

    for ( i = 0 ; i < MAXTOKENS ; i++ ) tokens[i] = NULL ;
    i = 0 ;

    strcpy( parseline_save, parseline ) ;

    if ((tokens[i] = my_strtok(parseline)) != NULL ) i++ ;
    while (( tokens[i] = my_strtok(NULL) ) != NULL && i < MAXTOKENS ) i++ ;
    if ( i <= 0 ) continue ;

    for ( i = 0 ; i < MAXTOKENS ; i++ )
    {
      if ( tokens[i] == NULL ) tokens[i] = "-" ;
    }
/* remove quotes from quoted entries */
    for ( i = 0 ; i < MAXTOKENS ; i++ )
    {
      char * pp ;
      if ( tokens[i][0] == '"' ) tokens[i]++ ;
      if ((pp=rindex( tokens[i], '"' )) != NULL ) *pp = '\0' ;
    }
    if      ( !strcmp( tokens[ TABLE ] , "state" ) )
    {
        int inbrace = 0 ;
        strcpy( newdims, "" ) ;
        is4d = 0 ; wantstend = 0 ; wantsbdy = 0 ; 
        for ( i = 0 ; i < (len_of_tok = strlen(tokens[F_DIMS])) ; i++ )
        {
          x = tolower(tokens[F_DIMS][i]) ;
          if ( x == '{' ) { inbrace = 1 ; }
          if ( x == '}' ) { inbrace = 0 ; }
          sprintf(xstr,"%c",x) ;
          if ( x != 'b' || inbrace ) strcat ( newdims , xstr ) ;
        }
    }
normal:
    /* otherwise output the line as is */
    fprintf(outfile,"%s\n",parseline_save) ;
    parseline[0] = '\0' ;  /* reset parseline */
  }
  return(retval) ;
}

int
reg_parse( FILE * infile )
{
  /* Had to increase size for SOA from 4096 to 7000, Manish Shrivastava 2010 */
  char inln[INLN_SIZE], parseline[PARSELINE_SIZE] ;
  char *p, *q ;
  char *tokens[MAXTOKENS], *toktmp[MAXTOKENS] ; 
  int i, ii ;
  int defining_state_field, defining_rconfig_field, defining_i1_field ;

  parseline[0] = '\0' ;

  max_time_level = 1 ;

/* main parse loop over registry lines */
  while ( fgets ( inln , INLN_SIZE , infile ) != NULL )
  {
    strcat( parseline , inln ) ; 
    /* allow \ to continue the end of a line */
    if (( p = index( parseline,  '\\'  )) != NULL )
    {
      if ( *(p+1) == '\n' || *(p+1) == '\0' )
      {
	*p = '\0' ;
	continue ;  /* go get another line */
      }
    }

    make_lower( parseline ) ;
    if (( p = index( parseline , '#' ))  != NULL  ) *p = '\0' ; /* discard comments (dont worry about quotes for now) */
    if (( p = index( parseline , '\n' )) != NULL  ) *p = '\0' ; /* discard newlines */
    for ( i = 0 ; i < MAXTOKENS ; i++ ) tokens[i] = NULL ; 
    i = 0 ;

    if ((tokens[i] = my_strtok(parseline)) != NULL ) i++ ; 

    while (( tokens[i] = my_strtok(NULL) ) != NULL && i < MAXTOKENS ) i++ ;
    if ( i <= 0 ) continue ;

    for ( i = 0 ; i < MAXTOKENS ; i++ )
    {
      if ( tokens[i] == NULL ) tokens[i] = "-" ;
    }

/* remove quotes from quoted entries */
    for ( i = 0 ; i < MAXTOKENS ; i++ )
    {
      char * pp ;
      if ( tokens[i][0] == '"' ) tokens[i]++ ;
      if ((pp=rindex( tokens[i], '"' )) != NULL ) *pp = '\0' ;
    }

    defining_state_field = 0 ;
    defining_rconfig_field = 0 ;
    defining_i1_field = 0 ;

/* state entry */
    if      ( !strcmp( tokens[ TABLE ] , "state" ) )
    {
      /* turn a state entry into a typedef to define a field in the top-level built-in type domain */
      tokens[TABLE] = "typedef" ;
      for ( i = MAXTOKENS-1 ; i >= 2 ; i-- ) tokens[i] = tokens[i-1] ; /* shift the fields to the left */
      tokens[FIELD_OF] = "domain" ;
                 if ( !strcmp( tokens[FIELD_TYPE], "double" ) ) tokens[FIELD_TYPE] = "doubleprecision" ; 
      defining_state_field = 1 ;
    }

    /* NOTE: fall through */

/* typedef entry */
    if ( !strcmp( tokens[ TABLE ] , "typedef" ) )
    {
      node_t * field_struct ;
      node_t * type_struct ;
      node_t * modname_struct ;

// FAST registry construct a list of module nodes
      modname_struct = get_modname_entry( tokens[ FIELD_MODNAME ] ) ;
      if ( modname_struct == NULL ) 
      {
        modname_struct = new_node( MODNAME ) ;
        strcpy( modname_struct->name, tokens[FIELD_MODNAME] ) ;
         modname_struct->module_ddt_list = NULL ;
        add_node_to_end( modname_struct , &ModNames ) ;
      }
//

      if ( !defining_state_field && ! defining_i1_field && 
           !defining_rconfig_field && !strcmp(tokens[FIELD_OF],"domain") )
       { fprintf(stderr,"Registry warning: 'domain' is a reserved registry type name. Cannot 'typedef domain'\n") ; }


      type_struct = get_type_entry( tokens[ FIELD_OF ] ) ;
      if ( type_struct == NULL ) 
      {  
        type_struct = new_node( TYPE ) ;
        strcpy( type_struct->name, tokens[FIELD_OF] ) ;
        type_struct->type_type = DERIVED ;
fprintf(stderr,"calling add_node_to_end %d %s\n",__LINE__,tokens[FIELD_OF]) ;
        add_node_to_end( type_struct , &Type ) ;
      }

      field_struct = new_node( FIELD ) ;

      strcpy( field_struct->name, tokens[FIELD_SYM] ) ;

      if ( set_state_type( tokens[FIELD_TYPE], field_struct ) )
       { fprintf(stderr,"Registry warning: type %s used before defined \n",tokens[FIELD_TYPE] ) ; }

      if ( set_state_dims( tokens[FIELD_DIMS], field_struct ) )
       { fprintf(stderr,"Registry warning: some problem with dimstring %s\n", tokens[FIELD_DIMS] ) ; }

#ifdef FUTURE
      field_struct->restart  = 0 ; field_struct->boundary  = 0 ;
      for ( i = 0 ; i < MAX_STREAMS ; i++ ) { 
        reset_mask( field_struct->io_mask, i ) ;
      }
#endif

// process CTRL keys -- only 'h' (hidden) and 'e' (exposed).  Default is not to generate a wrapper,
// so something must be specified, either h or e
      {
	char prev = '\0' ;
	char x ;
        char tmp[NAMELEN], tmp1[NAMELEN], tmp2[NAMELEN] ;
	int len_of_tok ;
        char fcn_name[2048], aux_fields[2048] ;

        strcpy(tmp,tokens[FIELD_CTRL]) ;
        if (( p = index(tmp,'=') ) != NULL ) { *p = '\0' ; }
        for ( i = 0 ; i < strlen(tmp) ; i++ )
        {
	  x = tolower(tmp[i]) ;
          if        ( x == 'h' ) {
            field_struct->gen_wrapper = WRAP_HIDDEN_FIELD ;
          } else if ( x == 'e' ) {
            field_struct->gen_wrapper = WRAP_EXPOSED_FIELD ;
          } else {
            field_struct->gen_wrapper = WRAP_NONE ;  /* default */
          }
        }
      }

      field_struct->dname[0] = '\0' ;
      if ( strcmp( tokens[FIELD_DNAME], "-" ) ) /* that is, if not equal "-" */
        { strcpy( field_struct->dname , tokens[FIELD_DNAME] ) ; }
      strcpy(field_struct->descrip,"-") ;
      if ( strcmp( tokens[FIELD_DESCRIP], "-" ) ) /* that is, if not equal "-" */
        { strcpy( field_struct->descrip , tokens[FIELD_DESCRIP] ) ; }
      strcpy(field_struct->units,"-") ;
      if ( strcmp( tokens[FIELD_UNITS], "-" ) ) /* that is, if not equal "-" */
        { strcpy( field_struct->units , tokens[FIELD_UNITS] ) ; }

#ifdef OVERSTRICT
      if ( field_struct->type != NULL )
        if ( field_struct->type->type_type == DERIVED && field_struct->ndims > 0 )
          { fprintf(stderr,"Registry warning: type item %s of type %s can not be multi-dimensional ",
	  		   tokens[FIELD_SYM], tokens[FIELD_TYPE] ) ; }
#endif

      add_node_to_end( field_struct , &(type_struct->fields) ) ;
      add_node_to_end( field_struct, &(modname_struct->module_ddt_list) ) ;

    }

/* dimespec entry */
    else if ( !strcmp( tokens[ TABLE ] , "dimspec" ) )
    {
      node_t * dim_struct ;
      dim_struct = new_node( DIM ) ;
      if ( get_dim_entry ( tokens[DIM_NAME] ) != NULL )
        { fprintf(stderr,"Registry warning: dimspec (%s) already defined\n",tokens[DIM_NAME] ) ; }
      strcpy(dim_struct->dim_name,tokens[DIM_NAME]) ;
      if ( set_dim_len( tokens[DIM_SPEC], dim_struct ) )
        { fprintf(stderr,"Registry warning: problem with dimspec (%s)\n",tokens[DIM_SPEC] ) ; }

      add_node_to_end( dim_struct , &Dim ) ;
    }

#if 0
     fprintf(stderr,"vvvvvvvvvvvvvvvvvvvvvvvvvvv\n") ;
     show_nodelist( Type ) ;
     fprintf(stderr,"^^^^^^^^^^^^^^^^^^^^^^^^^^^\n") ;
#endif
     parseline[0] = '\0' ;  /* reset parseline */
  }

/* Domain is a type node with fields that are not part of any type. WRF "state" entries
   were these. They were simply fields of the data type for a domain (as opposed to
   fields within derived data types that were fields in a domain).  The FAST registry
   does not have the concept of a Domain.  Leave the following assignment here but 
   put a test around it so we do not segfault if there aren't any "state" entries. */
  if ( get_type_entry( "domain" ) ) {
    Domain = *(get_type_entry( "domain" )) ;
  } 

  return(0) ;

}

node_t *
get_dim_entry( char *s )
{
  node_t * p ;
  for ( p = Dim ; p != NULL ; p = p->next )
  {
    if ( !strcmp(p->dim_name, s ) ) {
      return( p ) ;
    }
  }
  /* not found, check if dimension is specified in line */
  if ( 1 ) {
    node_t * dim_struct ;
    dim_struct = new_node( DIM ) ;
    strncpy(dim_struct->dim_name,s,1) ;
    if ( set_dim_len( s, dim_struct ) )
    { 
      fprintf(stderr,"Registry warning: problem with dimspec (%s)\n",s ) ;
    }
    else
    {
      add_node_to_end( dim_struct , &Dim ) ;
      return( dim_struct ) ;
    }
  }
  return(NULL) ;
}

int
set_state_type( char * typename, node_t * state_entry )
{
  if ( typename == NULL ) return(1) ;
  return (( state_entry->type = get_type_entry( typename )) == NULL )  ;
}

int
set_dim_len ( char * dimspec , node_t * dim_entry )
{
  if      (!strcmp( dimspec , "standard_domain" ))
   { dim_entry->len_defined_how = DOMAIN_STANDARD ; }
  else if (!strncmp( dimspec, "constant=" , 9 ) || isNum(dimspec[0]) )
  {
    char *p, *colon, *paren ;
    p = isNum(dimspec[0])?dimspec:&(dimspec[9]) ;
    /* check for colon */
    if (( colon = index(p,':')) != NULL )
    {
      *colon = '\0' ;
      if (( paren = index(p,'(')) !=NULL )
      {
        dim_entry->coord_start = atoi(paren+1) ;
      }
      else
      {
        fprintf(stderr,"WARNING: illegal syntax (missing opening paren) for constant: %s\n",p) ;
      }
      dim_entry->coord_end   = atoi(colon+1) ;
    }
    else
    {
      dim_entry->coord_start = 1 ;
      dim_entry->coord_end   = atoi(p) ;
    }
    dim_entry->len_defined_how = CONSTANT ;
  }
  else if (!strncmp( dimspec, "namelist=", 9 ))
  {
    char *p, *colon ;

    p = &(dimspec[9]) ;
    /* check for colon */
    if (( colon = index(p,':')) != NULL )
    {
      *colon = '\0' ;
      strcpy( dim_entry->assoc_nl_var_s, p ) ;
      strcpy( dim_entry->assoc_nl_var_e, colon+1 ) ;
    }
    else
    {
      strcpy( dim_entry->assoc_nl_var_s, "1" ) ;
      strcpy( dim_entry->assoc_nl_var_e, p ) ;
    }
    dim_entry->len_defined_how = NAMELIST ;
  }
  else
  {
    return(1) ;
  }
  return(0) ;
}


/* integrity checking of dimension list */
int
check_dimspecs()
{
  return(0) ;
}

int
init_parser()
{
  return(0) ;
}
