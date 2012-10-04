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

/* fields for state entries (note, these get converted to field entries in the
   reg_parse routine; therefore, only TABLE needs to be looked at */
#define TABLE 0

/* fields for field entries (TABLE="typedef" and, with some munging,  TABLE="state") */
#define FIELD_MODNAME   1
#define FIELD_OF        2
#define FIELD_TYPE      3
#define FIELD_SYM       4
#define FIELD_DIMS      5
#define FIELD_INIVAL     6
#define FIELD_CTRL      7
#define FIELD_DESCRIP   8
#define FIELD_UNITS     9

#define F_MODNAME  0
#define F_OF       1
#define F_TYPE     2
#define F_SYM      3
#define F_DIMS     4
#define F_INIVAL    4
#define F_CTRL     6
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
  char *tokens[MAXTOKENS], *toktmp[MAXTOKENS],*ditto[MAXTOKENS] ; 
  int i, ii ;
  int defining_state_field, defining_rconfig_field, defining_i1_field ;

  parseline[0] = '\0' ;

  max_time_level = 1 ;

  for ( i = 0 ; i < MAXTOKENS ; i++ ) { ditto[i] = (char *)malloc(NAMELEN) ; strcpy(ditto[i],"-") ; }

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

    //make_lower( parseline ) ;
    if (( p = index( parseline , '#' ))  != NULL  ) *p = '\0' ; /* discard comments (dont worry about quotes for now) */
    if (( p = index( parseline , '\n' )) != NULL  ) *p = '\0' ; /* discard newlines */
    if (( p = index( parseline , '\r' )) != NULL  ) *p = '\0' ; /* discard carriage returns (happens on Windows)*/
    for ( i = 0 ; i < MAXTOKENS ; i++ ) tokens[i] = NULL ; 
    i = 0 ;

    if ((tokens[i] = my_strtok(parseline)) != NULL ) i++ ; 
    while (( tokens[i] = my_strtok(NULL) ) != NULL && i < MAXTOKENS ) i++ ;
    if ( i <= 0 ) continue ;


    for ( i = 0 ; i < MAXTOKENS ; i++ )
    {
      if ( tokens[i] == NULL ) tokens[i] = "-" ;
      if ( strcmp(tokens[i],"^") ) {   // that is, if *not* ^
        strcpy(ditto[i],tokens[i]) ;
      } else {                         // if is ^
        tokens[i] = ditto[i] ;
      }
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

/* typedef entry */
    if ( !strcmp( tokens[ TABLE ] , "typedef" ) )
    {
      node_t * field_struct ;
      node_t * type_struct ;
      node_t * modname_struct ;
      char tmpstr[NAMELEN] ;

// FAST registry, construct a list of module nodes
      strcpy(tmpstr, make_lower_temp(tokens[ FIELD_MODNAME ])) ;
      if ( (p = index(tmpstr,'/')) != NULL ) *p = '\0' ;
      modname_struct = get_modname_entry( tmpstr ) ;
      if ( modname_struct == NULL ) 
      {
        char *p ;
        modname_struct = new_node( MODNAME ) ;
        strcpy( modname_struct->name, tokens[FIELD_MODNAME] ) ;
        // if a shortname is indicated after a slash, record that, otherwise use full name for both 
        if ( (p = index(modname_struct->name,'/')) != NULL ) {
          *p = '\0' ;
          strcpy( modname_struct->nickname, p+1 ) ;
        } else {
          strcpy( modname_struct->nickname, modname_struct->name ) ;
        }
        
        modname_struct->module_ddt_list = NULL ;
        modname_struct->next            = NULL ;
        add_node_to_end( modname_struct , &ModNames ) ;
      }

// FAST registry, construct list of derived data types specified for the Module 
      if ( strcmp(modname_struct->nickname,"") ) {
        sprintf(tmpstr,"%s_%s",modname_struct->nickname,make_lower_temp( tokens[ FIELD_OF ])) ;
      } else {
        sprintf(tmpstr,"%s",make_lower_temp( tokens[ FIELD_OF ])) ; 
      }
      type_struct = get_entry( tmpstr, modname_struct->module_ddt_list ) ;
      if ( type_struct == NULL ) 
      {  
        type_struct = new_node( TYPE ) ;
        strcpy( type_struct->name, tmpstr ) ;
        if ( strcmp(modname_struct->nickname,"") ) {
          sprintf(type_struct->mapsto,"%s_%s",modname_struct->nickname, tokens[ FIELD_OF ]) ;
        } else {
          sprintf(type_struct->mapsto,"%s", tokens[ FIELD_OF ]) ; 
        }
        type_struct->type_type = DERIVED ;
        type_struct->next      = NULL ;
        add_node_to_end( type_struct, &(modname_struct->module_ddt_list) ) ;
      }

// FAST registry, construct the list of fields in the derived types in the Module
      field_struct = new_node( FIELD ) ;
      strcpy( field_struct->name, tokens[FIELD_SYM] ) ;
      if ( set_state_type( tokens[FIELD_TYPE], field_struct ) )
       { fprintf(stderr,"Registry warning: type %s used before defined for %s\n",tokens[FIELD_TYPE],tokens[FIELD_SYM] ) ; }
      if ( set_state_dims( tokens[FIELD_DIMS], field_struct ) )
       { fprintf(stderr,"Registry warning: some problem with dimstring %s for %s\n", tokens[FIELD_DIMS],tokens[FIELD_SYM] ) ; }
      if ( set_ctrl( tokens[FIELD_CTRL], field_struct ) )
       { fprintf(stderr,"Registry warning: some problem with ctrl %s for %s\n", tokens[FIELD_CTRL],tokens[FIELD_SYM] ) ; }

      field_struct->inival[0] = '\0' ;
      if ( strcmp( tokens[FIELD_INIVAL], "-" ) ) /* that is, if not equal "-" */
        { strcpy( field_struct->inival , tokens[FIELD_INIVAL] ) ; }
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

    }

/* dimspec entry */
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
      fprintf(stderr,"Registry warning: get_dim_entry: problem with dimspec (%s)\n",s ) ;
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
  node_t *p ;
  int retval ;
  
  if ( typename == NULL ) return(1) ;
  retval = 0 ;
  if ( ( state_entry->type = get_type_entry( make_lower_temp(typename) )) == NULL ) {
    if ( !strncmp(make_lower_temp(typename),"character",9) ) 
    {
      p = new_node( TYPE ) ;
      strcpy( p->name, make_lower_temp(typename) ) ;
      strcpy( p->mapsto, typename ) ;
      add_node_to_end( p , &(state_entry->type) ) ;
    } else {
      retval = 1 ;
    }
  }
  return(retval) ;
}

int
set_dim_len ( char * dimspec , node_t * dim_entry )
{
  dim_entry->deferred = 0 ;
  if      (!strcmp( dimspec , "standard_domain" ))
   { dim_entry->len_defined_how = DOMAIN_STANDARD ; }
  else if (!strncmp( dimspec, "constant=" , 9 ) || isNum(dimspec[0]) || dimspec[0] == ':' || dimspec[0] == '(' )
  {
    char *p, *colon, *paren ;
    p = (isNum(dimspec[0])||dimspec[0]==':'||dimspec[0]=='(')?dimspec:&(dimspec[9]) ;
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
        dim_entry->deferred = 1 ;
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

int
set_ctrl( char *ctrl , node_t * field_struct )
// process CTRL keys -- only 'h' (hidden) and 'e' (exposed).  Default is not to generate a wrapper,
// so something must be specified, either h or e
{
  char prev = '\0' ;
  char x ;
  char tmp[NAMELEN] ;
  char *p ;
  int i ;
  strcpy(tmp,ctrl) ;
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
