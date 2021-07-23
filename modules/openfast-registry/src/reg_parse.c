#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
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
#define FIELD_INIVAL    6
#define FIELD_CTRL      7
#define FIELD_DESCRIP   8
#define FIELD_UNITS     9

#define F_MODNAME  0
#define F_OF       1
#define F_TYPE     2
#define F_SYM      3
#define F_DIMS     4
#define F_INIVAL   5
#define F_CTRL     6
#define F_DESCRIP  7
#define F_UNITS    8

/* fields for dimension entries (TABLE="dimspec") */
#define DIM_NAME       1
//#define DIM_ORDER      2
#define DIM_SPEC       2

#define INLN_SIZE      8000
#define PARSELINE_SIZE 8000

int isNum( char c )
{
    if ( c < '0' || c > '9' ) return 0;
    return 1 ;
}

int
pre_parse( char * dir, FILE * infile, FILE * outfile, int usefrom_sw )
{
  /* Decreased size for SOA from 8192 to 8000--double check if necessary, Manish Shrivastava 2010 */
  char inln[INLN_SIZE], parseline[PARSELINE_SIZE], parseline_save[PARSELINE_SIZE] ;
  char *p, *q, *p1, *p2  ;
  char *tokens[MAXTOKENS] ;
  int i, ifile ;
  int ifdef_stack_ptr = 0 ;
  int ifdef_stack[100] ;
  int inquote, retval ;
  int foundit ;

  ifdef_stack[0] = 1 ;
  retval = 0 ;

  parseline[0] = '\0' ;
  while ( fgets ( inln , INLN_SIZE , infile ) != NULL )
  {
/*** preprocessing directives ****/
    /* look for an include statement */
    if (( p = index( inln , '\n' )) != NULL  ) *p = '\0' ; /* discard newlines */
    if (( p = index( inln , '\r' )) != NULL  ) *p = '\0' ; /* discard carriage returns (happens on Windows)*/
    for ( p = inln ; ( *p == ' ' || *p == '\t' ) && *p != '\0' ; p++ ) ;
    p1 = make_lower_temp(p) ;
    if ( (!strncmp( p1 , "include", 7 ) || !strncmp( p1, "usefrom", 7 ))  &&  ! ( ifdef_stack_ptr >= 0 && ! ifdef_stack[ifdef_stack_ptr] ) )
    {
      FILE *include_fp ;
      char include_file_name[NAMELEN] ;
      char include_file_name_tmp[NAMELEN] ;
      int checking_for_usefrom = !strncmp( p1, "usefrom", 7 ) ;
//fprintf(stderr,"checking_for_usefrom %d |%s|\n",checking_for_usefrom,p1) ;

      p += 7 ; for ( ; ( *p == ' ' ||  *p == '\t' ) && *p != '\0' ; p++ ) ;
      if ( strlen( p ) > 127 ) { fprintf(stderr,"Registry warning: invalid include file name: %s\n", p ) ; }
      else {
/* look in a few places for valid include files */
        foundit = 0 ;

       // See if it might be in the current directory
        sprintf( include_file_name , "%s", p ) ;            // first name in line from registry file, without the include or usefrom
        for ( p2 = include_file_name ; !( *p2 == ' ' || *p2 == '\t' || *p2 == '\n' ) && *p2 != '\0' ; p2++ ) {} 
        *p2 = '\0' ;     // drop tailing white space
        if ( (q=index(include_file_name,'\n')) != NULL ) *q = '\0' ;
        if (( include_fp = fopen( include_file_name , "r" )) != NULL )   { foundit = 1 ; goto gotit ; }

        // See if it might be in the directory specified (or whatever dir is). Don't remove spaces from the dir name though.
        sprintf( include_file_name , "%s", p ) ;            // first name in line from registry file, without the include or usefrom
        for ( p2 = include_file_name ; !( *p2 == ' ' || *p2 == '\t' || *p2 == '\n' ) && *p2 != '\0' ; p2++ ) {} 
        *p2 = '\0' ;     // drop tailing white space
        sprintf( include_file_name , "%s/%s", dir, p );              // set the dir + file
        if ( (q=index(include_file_name,'\n')) != NULL ) *q = '\0' ;
        if (( include_fp = fopen( include_file_name , "r" )) != NULL )   { foundit = 1 ; goto gotit ; }

        // Check in the list of include dirs
        for ( ifile = 0 ; ifile < nincldirs ; ifile++ ) {
          sprintf( include_file_name_tmp , "%s", p ) ;            // first name in line from registry file, without the include or usefrom
          for ( p2 = include_file_name_tmp ; !( *p2 == ' ' || *p2 == '\t' || *p2 == '\n' ) && *p2 != '\0' ; p2++ ) {}
          *p2 = '\0' ;     // drop tailing white space
          sprintf( include_file_name, "%s/%s", IncludeDirs[ifile] , include_file_name_tmp ) ;     // dir specified with -I
          if ( (q=index(include_file_name,'\n')) != NULL ) *q = '\0' ;
          if (( include_fp = fopen( include_file_name , "r" )) != NULL ) { foundit = 1 ; goto gotit ; }
        }

        // Cygwin specific -- assuming spaces in dir are ok.
        for ( ifile = 0 ; ifile < nincldirs ; ifile++ ) {
          int drive_specified = 0 ;
          sprintf( include_file_name_tmp , "%s", p ) ;            // first name in line from registry file, without the include or usefrom
          for ( p2 = include_file_name_tmp ; !( *p2 == ' ' || *p2 == '\t' || *p2 == '\n' ) && *p2 != '\0' ; p2++ ) {}  
          *p2 = '\0' ;
          sprintf( include_file_name , "%s/%s", IncludeDirs[ifile] , include_file_name_tmp ) ;     // dir munged for cigwin
          if ( include_file_name[0] == '/' ) {
            char tmp[NAMELEN], tmp2[NAMELEN], *dr ;
            strcpy( tmp2, include_file_name ) ;
            if ( !strncmp( tmp2, "/cygdrive/", 10 )) {
               strcpy(tmp,tmp2+11) ; // skip past /cygdrive/c
               strcpy(tmp2,tmp) ;
               drive_specified = 1 ;
            }
            for ( dr = "abcdefmy" ; *dr ; dr++ ) {
              sprintf(tmp,"%c:%s%s",*dr,(drive_specified)?"":"/cygwin",tmp2) ;
              strcpy( include_file_name, tmp ) ;
              for ( p2 = include_file_name ; !( *p2 == ' ' || *p2 == '\t' || *p2 == '\n' ) && *p2 != '\0' ; p2++ ) {}  
              *p2 = '\0' ;
              if ( (q=index(include_file_name,'\n')) != NULL ) *q = '\0' ;
              if (( include_fp = fopen( include_file_name , "r" )) != NULL ) { foundit = 1 ; goto gotit ; }
            }
          }
        }

gotit:
        if ( foundit ) {
          fprintf(stderr,"opening %s %s\n",include_file_name,
                                           (checking_for_usefrom || usefrom_sw)?"in usefrom mode":"" ) ;
          parseline[0] = '\0' ;
          pre_parse( dir , include_fp , outfile, ( checking_for_usefrom + usefrom_sw ) ) ;
          parseline[0] = '\0' ;
//          fprintf(stderr,"closing %s %s\n",include_file_name,
//                                           (checking_for_usefrom || usefrom_sw)?"in usefrom mode":"" ) ;
          fclose( include_fp ) ;
          continue ;
        } else {
          if ( ! checking_for_usefrom ) {
            fprintf(stderr,"Registry warning: cannot open %s . Ignoring.\n", include_file_name ) ;
          }
        }
      }
    }
    else if ( !strncmp( make_lower_temp(p) , "ifdef", 5 ) ) {
      char value[32] ;
      p += 5 ; for ( ; ( *p == ' ' || *p == '\t' ) && *p != '\0' ; p++ ) ;
      strncpy(value, p, 31 ) ; value[31] = '\0' ;
      if ( (p=index(value,'\n')) != NULL ) *p = '\0' ;
      if ( (p=index(value,' '))  != NULL ) *p = '\0' ;
      if ( (p=index(value,'\t')) != NULL ) *p = '\0' ;
      ifdef_stack_ptr++ ;
      ifdef_stack[ifdef_stack_ptr] = ( sym_get(value) != NULL && ifdef_stack[ifdef_stack_ptr-1] ) ;
      if ( ifdef_stack_ptr >= 100 ) { fprintf(stderr,"Registry fatal: too many nested ifdefs\n") ; exit(1) ; }
      continue ;
    }
    else if ( !strncmp( make_lower_temp(p) , "ifndef", 6 ) ) {
      char value[32] ;
      p += 6 ; for ( ; ( *p == ' ' || *p == '\t') && *p != '\0' ; p++ ) ;
      strncpy(value, p, 31 ) ; value[31] = '\0' ;
      if ( (p=index(value,'\n')) != NULL ) *p = '\0' ;
      if ( (p=index(value,' '))  != NULL ) *p = '\0' ;
      if ( (p=index(value,'\t')) != NULL ) *p = '\0' ;
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
      p += 6 ; for ( ; ( *p == ' ' || *p == '\t') && *p != '\0' ; p++ ) ;
      strncpy(value, p, 31 ) ; value[31] = '\0' ;
      if ( (p=index(value,'\n')) != NULL ) *p = '\0' ;
      if ( (p=index(value,' '))  != NULL ) *p = '\0' ;
      if ( (p=index(value,'\t')) != NULL ) *p = '\0' ;
      sym_add( value ) ;
      continue ;
    }
    if ( ifdef_stack_ptr >= 0 && ! ifdef_stack[ifdef_stack_ptr] ) continue ;
/*** end of preprocessing directives ****/
//fprintf(stderr,"parseline |%s|\n",parseline) ;
//fprintf(stderr,"inln |%s|\n",inln) ;

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

    // get parsline_save, the value written to the output file...
    //fprintf(stderr,"parseline_save |%s|\n",parseline_save) ;
    //strcpy(parseline_save, parseline);
    for (p = parseline; (*p == ' ' || *p == '\t') && *p != '\0'; p++);
    strcpy(parseline_save, p);  // get rid of leading spaces

    if (!strncmp(parseline_save, "typedef", 7))
    {
       char tmp[PARSELINE_SIZE], *x;
       strcpy(tmp, parseline_save);
       x = strpbrk(tmp, " \t"); // find the first space or tab
       if (usefrom_sw && x) {
          sprintf(parseline_save, "usefrom%i %s", usefrom_sw, x);
       }
    }

    // parse tokens from parseline
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



//normal:
    /* otherwise output the line as is */
    fprintf(outfile,"%s\n",parseline_save) ;
    parseline[0] = '\0' ;  /* reset parseline */
    parseline_save[0] = '\0' ;  /* reset parseline_save */
  }
  return(retval) ;
}

int
reg_parse( FILE * infile )
{
  /* Had to increase size for SOA from 4096 to 7000, Manish Shrivastava 2010 */
  char inln[INLN_SIZE], parseline[PARSELINE_SIZE] ;
  char *p ;
  char *tokens[MAXTOKENS],*ditto[MAXTOKENS] ;
  int i  ;
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

/* typedef, usefrom, and param entries */
//         || !strcmp(tokens[TABLE], "usefrom")
    if (  !strcmp( tokens[ TABLE ] , "typedef" )
        || !strncmp(tokens[TABLE], "usefrom", 7) 
        || !strcmp( tokens[ TABLE ] , "param"   ) )
    {
      node_t * param_struct ;
      node_t * field_struct ;
      node_t * type_struct ;
      node_t * modname_struct ;
      char tmpstr[NAMELEN], ddtname[NAMELEN] ;

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
      if (!strcmp(tokens[TABLE], "usefrom"))
      {
          modname_struct->usefrom = 1;
      } else if(!strncmp(tokens[TABLE], "usefrom", 7))
      {
          tokens[TABLE] += 7;
          if (!strcmp(tokens[TABLE], "1"))
          {
              modname_struct->usefrom = 1;
          }
          else
          {
              modname_struct->usefrom = 2;
          }
      }

      if ( !strcmp( tokens[ TABLE ] , "param" ) ) {
// FAST registry, construct list of params specified for the Module
        param_struct = new_node( PARAM ) ;
        sprintf(param_struct->name,"%s",tokens[ FIELD_SYM ]) ; // name of parameter
        if ( set_state_type( tokens[FIELD_TYPE], param_struct, Type, NULL  ) ) // Only search type list, not ddts for module
         { fprintf(stderr,"Registry warning: type %s used before defined for %s\n",tokens[FIELD_TYPE],tokens[FIELD_SYM] ) ; }
        if ( set_state_dims( tokens[FIELD_DIMS], param_struct ) )
         { fprintf(stderr,"Registry warning: some problem with dimstring %s for %s\n", tokens[FIELD_DIMS],tokens[FIELD_SYM] ) ; }
        param_struct->inival[0] = '\0' ;
        if ( strcmp( tokens[FIELD_INIVAL], "-" ) ) /* that is, if not equal "-" */
          { strcpy( param_struct->inival , tokens[FIELD_INIVAL] ) ; }
        strcpy(param_struct->descrip,"-") ;
        if ( strcmp( tokens[FIELD_DESCRIP], "-" ) ) /* that is, if not equal "-" */
          { strcpy( param_struct->descrip , tokens[FIELD_DESCRIP] ) ; }
        strcpy(param_struct->units,"-") ;
        if ( strcmp( tokens[FIELD_UNITS], "-" ) ) /* that is, if not equal "-" */
          { strcpy( param_struct->units , tokens[FIELD_UNITS] ) ; }

        add_node_to_end( param_struct , &(modname_struct->params) ) ;

      } else { // not param

// FAST registry, construct list of derived data types specified for the Module
//  Only the FAST interface defined types should have the Module's nickname prepended
        sprintf(ddtname,"%s",tokens[ FIELD_OF ]) ;
        modname_struct->is_interface_type = 0 ;
        if ( strcmp(modname_struct->nickname,"") ) {
          if ( is_a_fast_interface_type(tokens[FIELD_OF] ) ) {
            sprintf(ddtname,"%s_%s",modname_struct->nickname,tokens[ FIELD_OF ]) ;
            modname_struct->is_interface_type = 1 ;
          }
        }
        sprintf(tmpstr,"%s",make_lower_temp(ddtname)) ;
        type_struct = get_entry( tmpstr, modname_struct->module_ddt_list ) ;
        if ( type_struct == NULL && modname_struct->usefrom)
        {
            type_struct = get_entry( tmpstr, Type ) ;
        }

        if ( type_struct == NULL )
        {
          type_struct = new_node( TYPE ) ;
          strcpy( type_struct->name, tmpstr ) ;
          strcpy(type_struct->mapsto,ddtname) ;
          type_struct->type_type = DERIVED ;
          type_struct->next      = NULL ;
          type_struct->usefrom   = modname_struct->usefrom ;
          type_struct->module    = modname_struct ;
          add_node_to_end( type_struct,(type_struct->usefrom)? &Type : &(modname_struct->module_ddt_list ) ) ;
        }

// FAST registry, construct the list of fields in the derived types in the Module
        field_struct = new_node( FIELD ) ;
        strcpy( field_struct->name, tokens[FIELD_SYM] ) ;
        if ( set_state_type( tokens[FIELD_TYPE], field_struct, Type, modname_struct->module_ddt_list ) )
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
        field_struct->usefrom   = type_struct->usefrom ;

        add_node_to_end( field_struct , &(type_struct->fields) ) ;
      } // not param

    }

/* dimspec entry */
    else if ( !strcmp( tokens[ TABLE ] , "dimspec" ) )
    {
      node_t * dim_struct ;
      dim_struct = new_node( DIM ) ;
      if ( get_dim_entry ( tokens[DIM_NAME], 0 ) != NULL )
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
get_dim_entry( char *s, int sw ) // sw = 1 is used when checking an inline dimspec
{
  node_t * p ;
  for ( p = Dim ; p != NULL ; p = p->next )
  {
    if ( !strcmp(p->dim_name, s ) ) {
      return( p ) ;
    }
  }
  /* not found, check if dimension is specified in line */
  if ( 1  && sw ) {
    node_t * dim_struct ;
    dim_struct = new_node( DIM ) ;
    strcpy(dim_struct->dim_name,s) ;
//    strncpy(dim_struct->dim_name,s,1) ;
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
set_state_type( char * typename, node_t * state_entry, node_t * typelist, node_t * ddtlist )
{
  node_t *p ;
  int retval ;

  if ( typename == NULL ) return(1) ;
  retval = 0 ;
  if ( ( state_entry->type = get_entry( make_lower_temp(typename), ddtlist )) == NULL ) {
    if ( ( state_entry->type = get_entry( make_lower_temp(typename), typelist )) == NULL ) {
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
      else if ( isNum(*p) ) {
        dim_entry->coord_start = atoi(p) ;
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
  else /* if (param_dim != NULL) */ {
     dim_entry->coord_start = 1;
     dim_entry->len_defined_how = CONSTANT;
     strcpy(dim_entry->dim_param_name, dimspec);
     dim_entry->dim_param = 1;
  }
/*    else
  {
    return(1) ;
  }
*/
  return(0) ;
}

int
set_ctrl( char *ctrl , node_t * field_struct )
// process CTRL keys -- only '2pi' (interpolation of values with 2pi period).  Default is no special interpolation.
{
  char tmp[NAMELEN] ;
  char *p ;
  strcpy(tmp,ctrl) ;
  if (( p = index(tmp,'=') ) != NULL ) { *p = '\0' ; }
  if (!strcmp(make_lower_temp(tmp), "2pi")) {
      field_struct->gen_periodic = PERIOD_2PI;
  }
  else {
     field_struct->gen_periodic = PERIOD_NONE;
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

int
is_a_fast_interface_type( char *str )
{
   int retval ;

   retval =  (
     !strcmp(make_lower_temp(str), "initinputtype")       ||
     !strcmp(make_lower_temp(str), "initoutputtype")      ||
     !strcmp(make_lower_temp(str), "inputtype")           ||
     !strcmp(make_lower_temp(str), "outputtype")          ||
     !strcmp(make_lower_temp(str), "continuousstatetype") ||
     !strcmp(make_lower_temp(str), "discretestatetype")   ||
     !strcmp(make_lower_temp(str), "constraintstatetype") ||
     !strcmp(make_lower_temp(str), "otherstatetype")      ||
     !strcmp(make_lower_temp(str), "parametertype")       ||
     !strcmp(make_lower_temp(str), "miscvartype")         ||
     !strcmp(make_lower_temp(str), "partialoutputpinputtype") ||
     !strcmp(make_lower_temp(str), "partialcontstatepinputtype")         ||
     !strcmp(make_lower_temp(str), "partialdiscstatepinputtype")         ||
     !strcmp(make_lower_temp(str), "partialconstrstatepinputtype")       ||
            0 ) ;

   return(retval) ;
}

int
must_have_real_or_double( char *str )
{
   int retval ;

   retval =  (
     !strcmp(make_lower_temp(str), "inputtype")           ||
     !strcmp(make_lower_temp(str), "outputtype")          ||
     !strcmp(make_lower_temp(str), "continuousstatetype") ||
     !strcmp(make_lower_temp(str), "discretestatetype")   ||
     !strcmp(make_lower_temp(str), "constraintstatetype") ||
     !strcmp(make_lower_temp(str), "partialoutputpinputtype")            ||
     !strcmp(make_lower_temp(str), "partialcontstatepinputtype")         ||
     !strcmp(make_lower_temp(str), "partialdiscstatepinputtype")         ||
     !strcmp(make_lower_temp(str), "partialconstrstatepinputtype")       ||
            0 ) ;

   return(retval) ;
}

char *
fast_interface_type_shortname( char *str )
{
   char * retval, *str2;
   str2 = make_lower_temp(str);

   if        (  !strcmp(str2, "initinputtype") )  {
     retval = "InitInput" ;
   } else if (  !strcmp(str2, "initoutputtype") ) {
     retval = "InitOutput" ;
   } else if (  !strcmp(str2, "inputtype") ) {
     retval = "Input" ;
   } else if (  !strcmp(str2, "outputtype") ) {
     retval = "Output" ;
   } else if (  !strcmp(str2, "continuousstatetype") ) {
     retval = "ContState" ;
   } else if (  !strcmp(str2, "discretestatetype") )  {
     retval = "DiscState" ;
   } else if (  !strcmp(str2, "constraintstatetype") ) {
     retval = "ConstrState" ;
   } else if (  !strcmp(str2, "otherstatetype") ) {
     retval = "OtherState" ;
   } else if (  !strcmp(str2, "miscvartype") ) {
      retval = "Misc";
   } else if (  !strcmp(str2, "parametertype") ) {
     retval = "Param" ;
   } else if (  !strcmp(str2, "partialoutputpinputtype") ) {
     retval = "dYdu" ;
   } else if (  !strcmp(str2, "partialcontstatepinputtype") ) {
     retval = "dXdu" ;
   } else if (  !strcmp(str2, "partialdiscstatepinputtype") ) {
     retval = "dXddu" ;
   } else if (  !strcmp(str2, "partialconstrstatepinputtype") ) {
     retval = "dZdu" ;
   }
   else{ 
      retval = str; 
   }


   return(retval) ;
}

char *
std_case( char *str )    // returns the name in CamelBack case or just the name itself
{
   if      (  !strcmp(make_lower_temp(str), "initinputtype"))                {return("InitInputType");}
   else if (  !strcmp(make_lower_temp(str), "initoutputtype"))               {return("InitOutputType");}
   else if (  !strcmp(make_lower_temp(str), "inputtype"))                    {return("InputType");}
   else if (  !strcmp(make_lower_temp(str), "outputtype"))                   {return("OutputType");}
   else if (  !strcmp(make_lower_temp(str), "continuousstatetype"))          {return("ContinuousStateType");}
   else if (  !strcmp(make_lower_temp(str), "discretestatetype"))            {return("DiscreteStateType");}
   else if (  !strcmp(make_lower_temp(str), "constraintstatetype"))          {return("ConstraintStateType");}
   else if (  !strcmp(make_lower_temp(str), "otherstatetype"))               {return("OtherStateType");}
   else if (  !strcmp(make_lower_temp(str), "miscvartype"))                  {return("MiscVarType"); }
   else if (  !strcmp(make_lower_temp(str), "parametertype"))                {return("ParameterType"); }
   else if (  !strcmp(make_lower_temp(str), "partialoutputpinputtype"))      {return("PartialOutputPInputType");}
   else if (  !strcmp(make_lower_temp(str), "partialcontstatepinputtype"))   {return("PartialConstStatePInputType");}
   else if (  !strcmp(make_lower_temp(str), "partialdiscstatepinputtype"))   {return("PartialDiscStatePInputType");}
   else if (  !strcmp(make_lower_temp(str), "partialconstrstatepinputtype")) {return("PartialConstrStatePInputType");}
   else                                                                      {return(str);}
    // shouldn't happen
   return("") ;
}

