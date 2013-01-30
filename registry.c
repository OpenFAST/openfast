#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# include <io.h>
# define rindex(X,Y) strrchr(X,Y)
# define index(X,Y) strchr(X,Y)
#else
# include <sys/time.h>
# include <sys/resource.h>
# include <unistd.h>
# include <strings.h>
#endif

#define DEFINE_GLOBALS
#include "protos.h"
#include "registry.h"
#include "data.h"
#include "sym.h"

int
main( int argc, char *argv[], char *env[] )
{
  char fname_in[NAMELEN], dir[NAMELEN], fname_tmp[NAMELEN], command[NAMELEN] ;
  FILE * fp_in, *fp_tmp ;
  char * thisprog  ;
  int mypid ;
  int wrote_template ;
#ifndef _WIN32
  struct rlimit rlim ;
#endif

  mypid = (int) getpid() ;
  strcpy( thiscom, argv[0] ) ;
  argv++ ;

  sw_output_template_force = 0 ;
  strcpy( fname_in , "" ) ;

#ifndef _WIN32
  rlim.rlim_cur = RLIM_INFINITY ;
  rlim.rlim_max = RLIM_INFINITY ;
  setrlimit ( RLIMIT_STACK , &rlim ) ;
#endif


  fprintf(stderr,"----- FAST Registry  --------------\n") ;
  fprintf(stderr,"Revision $Rev$\n") ;
  fprintf(stderr,"Date $LastChangedDate$ \n" ) ;
  fprintf(stderr,"URL  $URL$\n" ) ;
  fprintf(stderr,"-----------------------------------\n") ;

  sym_forget() ;
//  thisprog = *argv ;
  thisprog = "registry.exe" ;
  strcpy(fname_in,"") ;
  wrote_template = 0 ;

  while (*argv) {
    if (*argv[0] == '-') {  /* an option */
      if (!strncmp(*argv,"-D",2)) {
        char * p ;
        p = *argv ;
        sym_add(p+2) ;
      }
      if (!strcmp(*argv,"-f")) {
        sw_output_template_force = 1 ;
      }
      if (!strcmp(*argv,"-template") || !strcmp(*argv,"-registry")) {
        char * arg ;
        arg = *argv ;
        argv++ ; if ( *argv ) { strcpy( sw_modname_subst,     *argv ) ; } else { goto usage ; }
        argv++ ; if ( *argv ) { strcpy( sw_modnickname_subst, *argv ) ; } else { goto usage ; }
        if (!strcmp(arg,"-template")) output_template(sw_modname_subst,sw_modnickname_subst,sw_output_template_force,0) ;
        if (!strcmp(arg,"-registry")) output_template(sw_modname_subst,sw_modnickname_subst,sw_output_template_force,1) ;
        wrote_template = 1 ;
      }
      if (!strncmp(*argv,"-h",2)) {
usage:
        fprintf(stderr,"Usage: %s [-D<MACRO>]  registryfile | [-f] [-template|-registry] ModuleName ModName \n",thisprog) ;
        exit(1) ;
      }
    }
    else  /* consider it an input file */
    {
      strcpy( fname_in , *argv ) ;
    }
    argv++ ;
  }
  if ( wrote_template ) exit(0) ;

  if ( !strcmp(fname_in,"") ) goto usage ;

#ifdef FUTURE
  gen_io_boilerplate() ;  /* 20091213 jm.  Generate the io_boilerplate_temporary.inc file */
#endif

  init_parser() ;
  init_type_table() ;
  init_dim_table() ;
  init_modname_table() ;

  if ( !strcmp(fname_in,"") ) fp_in = stdin ;
  else
    if (( fp_in = fopen( fname_in , "r" )) == NULL )
    {
      fprintf(stderr,"Registry program cannot open %s for reading. Ending.\n", fname_in ) ;
      exit(2) ;
    }
  
  sprintf( fname_tmp , "Registry_tmp.%d",mypid) ;
  if (( fp_tmp = fopen( fname_tmp  , "w" )) == NULL )
  {
    fprintf(stderr,"Registry program cannot open temporary %s for writing. Ending.\n", fname_tmp ) ;
    exit(2) ;
  }

  { char *e ;
    strcpy( dir , fname_in ) ;
    if ( ( e = rindex ( dir , '/' ) ) != NULL ) { *e = '\0' ; } else { strcpy( dir, "." ) ; } 
  }
  if ( pre_parse( dir, fp_in, fp_tmp ) ) {
    fprintf(stderr,"Problem with Registry File %s\n", fname_in ) ;
    goto cleanup ;
  }
  sym_forget() ;

  fclose(fp_in) ;
  fclose(fp_tmp) ;

  if (( fp_tmp = fopen( fname_tmp , "r" )) == NULL )
  {
    fprintf(stderr,"Registry program cannot open %s for reading. Ending.\n", fname_tmp ) ;
    goto cleanup ;
  }

  reg_parse(fp_tmp) ;

  fclose(fp_tmp) ;

  check_dimspecs() ;

  gen_module_files( "." ) ;


cleanup:
#ifdef _WIN32
   sprintf(command,"del /F /Q %s\n",fname_tmp );
#else
   sprintf(command,"/bin/rm -f %s\n",fname_tmp );
#endif
   system( command ) ;

   exit( 0 ) ;

}
#include "Template_data.c"
#include "Template_registry.c"

output_template( char * sw_modname_subst, char * sw_modnickname_subst, int * force, int sw  ) // sw = 0, template; 1 = registry
{
    char ** p ;
    FILE *fp ;
    char fname[NAMELEN] ;
    char tmp1[2096], tmp2[2096], tmp3[2096] ;
    if ( sw == 0 ) { sprintf(fname,"%s.f90",sw_modname_subst) ; }
    else           { sprintf(fname,"Registry_%s.txt",sw_modname_subst) ; }
    if ( ! force ) {
      if ( (fp = fopen( fname,"r" )) != NULL ) {
        fprintf(stderr,"Registry exiting. Attempt to overwrite file (%s) . Move out of the way or specify -f before -template option. \n", fname) ;
        exit(1) ;
      }
      fclose(fp) ;
    }
    if ( (fp = fopen( fname,"w" )) == NULL ) {
      fprintf(stderr,"Registry exiting. Failure opening %s.\n", fname ) ;
      exit(1) ;
    }
    if ( sw == 0 ) {
      for ( p = template_data ; *p ; p++ ) {
        strcpy(tmp1,*p) ;
        substitute(tmp1,"ModuleName",sw_modname_subst,tmp2) ;
        substitute(tmp2,"ModName",sw_modnickname_subst,tmp3) ;
        fprintf(fp,"%s\n",tmp3) ;
      }
    } else {
      for ( p = template_registry ; *p ; p++ ) {
        strcpy(tmp1,*p) ;
        substitute(tmp1,"ModuleName",sw_modname_subst,tmp2) ;
        substitute(tmp2,"ModName",sw_modnickname_subst,tmp3) ;
        fprintf(fp,"%s\n",tmp3) ;
      }
    }
    fclose(fp) ;
}



// would use regex for this but it does not seem to be uniformly or universally supported

int
substitute( char * str , char * match , char * replace, char * result )
{
   char * p, *q ;
   size_t n, m ;

   n = strlen( replace ) ;
   m = strlen( match ) ;
   for ( p = str , q = result ; *p ; )
   {
      if ( matches( p, match ) )
      {
        strncpy( q, replace, n ) ;
        q += n ; 
        p += m ;
      } else {
        *q = *p ;
        p++ ;
        q++ ;
      }
   }
   *q = '\0' ;
   strcpy( str, result ) ;
}

int 
matches( char * str , char * match )   // both must be null terminated
{
   char * p, * q ;
   int n, retval ;
   
   for ( n = 0, p = str, q = match ;  (*p && *q) ; p++, q++, n++ ) 
   {
     if ( *p != *q ) return(0) ;
   }
   if ( n != strlen(match) ) return(0) ;
   return(1) ;
}
