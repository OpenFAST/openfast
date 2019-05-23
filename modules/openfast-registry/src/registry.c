#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# include <io.h>
# define rindex(X,Y) strrchr(X,Y)
# define index(X,Y) strchr(X,Y)
# include <process.h>
# define getpid _getpid
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

void output_template( char * sw_modname_subst, char * sw_modnickname_subst, int force, int sw  );
int matches( char * str , char * match );

int
main( int argc, char *argv[], char *env[] )
{
  char fname_in[NAMELEN], dir[NAMELEN], fname_tmp[NAMELEN], command[NAMELEN] ;
  FILE * fp_in, *fp_tmp ;
  char * thisprog  ;
  char * thisprog_ver;
  int mypid ;
  int wrote_template ;
  int sw_keep = 0 ;
#ifndef _WIN32
  struct rlimit rlim ;
#endif

  mypid = (int) getpid() ;
  strcpy( thiscom, argv[0] ) ;
  argv++ ;

  sw_output_template_force = 0 ;
  sw_norealloc_lsh   = 1 ;
  sw_ccode           = 0 ;
  sw_noextrap        = 0 ;
  sw_shownodes       = 0 ;
  strcpy( fname_in , "" ) ;

#ifndef _WIN32
  rlim.rlim_cur = RLIM_INFINITY ;
  rlim.rlim_max = RLIM_INFINITY ;
  setrlimit ( RLIMIT_STACK , &rlim ) ;
#endif

   thisprog_ver = "FAST Registry";

  fprintf(stderr,"\n") ;
  fprintf(stderr,"----- %s --------------\n", thisprog_ver) ;
  fprintf(stderr,"----------------------------------------------------------\n") ;

  sym_forget() ;
 //thisprog = *argv ;
 // strcpy(thisprog, thiscom);
  thisprog = "registry.exe";
  strcpy(fname_in, "");
  strcpy(OutDir, "."); // if no OutDir is listed, use current directory
  wrote_template = 0;


  while (*argv) {

     if (!strncmp(*argv,"-D",2)) {
        char * p ;
        p = *argv ;
        sym_add(p+2) ;
      } else if (!strncmp(*argv,"/D=",3)) {
        char * p ;
        p = *argv ;
        sym_add(p+3) ;
      } else if (!strcmp(*argv,"-force") || !strcmp(*argv,"/force") ) {
        sw_output_template_force = 1 ;
      } else if (!strcmp(*argv,"-O") || !strcmp(*argv,"/O") ) {
        argv++ ; if ( *argv )  { strcpy( OutDir, *argv ) ; }
      } else if (!strcmp(*argv,"-I") || !strcmp(*argv,"/I") ) {
        argv++ ; if ( *argv ) { if( nincldirs < MAXINCLDIRS ) { strcpy( IncludeDirs[nincldirs++], *argv ) ; } }
      } else if (!strcmp(*argv, "-ccode") || !strcmp(*argv, "/ccode")) {
        sw_ccode = 1 ;
      } else if (!strcmp(*argv, "-noextrap") || !strcmp(*argv, "/noextrap")) {
          sw_noextrap = 1;
      } else if (!strncmp(*argv, "-shownodes", 4) || !strncmp(*argv, "/shownodes", 4)) {
        sw_shownodes = 1 ;
      } else if (!strcmp(*argv,"-template") || !strcmp(*argv,"-registry") ||
                 !strcmp(*argv,"/template") || !strcmp(*argv,"/registry")  ) {
        char * arg ;
        arg = *argv ;
        argv++ ; if ( *argv ) { strcpy( sw_modname_subst,     *argv ) ; } else { goto usage ; }
        argv++ ; if ( *argv ) { strcpy( sw_modnickname_subst, *argv ) ; } else { goto usage ; }
        if (!strcmp(arg+1,"template")) output_template(sw_modname_subst,sw_modnickname_subst,sw_output_template_force,0) ;
        if (!strcmp(arg+1,"registry")) output_template(sw_modname_subst,sw_modnickname_subst,sw_output_template_force,1) ;
        wrote_template = 1 ;
      } else if (!strcmp(*argv,"-h") || !strcmp(*argv,"/h")) {
usage:
//        fprintf(stderr,"Usage: %s [options] registryfile -or- \n",thisprog) ;
        fprintf(stderr, "Usage: %s registryfile [options] -or- \n",thiscom) ;
        fprintf(stderr, "          [-force] [-template|-registry] ModuleName ModName \n") ;
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "    -h                this summary\n");
        fprintf(stderr, "    -I <dir>          look for usefrom files in directory \"dir\"\n");
        fprintf(stderr, "    -O <dir>          generate types files in directory \"dir\"\n");
        fprintf(stderr, "    -noextrap         do not generate ModName_Input_ExtrapInterp or ModName_Output_ExtrapInterp routines\n");
        fprintf(stderr, "    -D<SYM>           define symbol for conditional evaluation inside registry file\n");
        fprintf(stderr, "    -ccode            generate additional code for interfacing with C/C++\n") ;
        fprintf(stderr, "    -keep             do not delete temporary files from registry program\n") ;
        fprintf(stderr, "    -shownodes        output a listing of the nodes in registry's AST\n") ;
        fprintf(stderr, "  === alternate usage for generating templates ===\n") ;
        fprintf(stderr, "    -template ModuleName ModName\n") ;
        fprintf(stderr, "                 Generate a template Module file none exists\n") ;
        fprintf(stderr, "    -registry ModuleName ModName\n") ;
        fprintf(stderr, "                 Generate a template registry file if none exists\n") ;
        fprintf(stderr, "    -force Force generating of template or registry file\n") ;
        fprintf(stderr, "  (the / character can be used in place of - when specifying options)\n") ;
        exit(1) ;
      } else if (!strcmp(*argv,"-keep") || !strcmp(*argv,"/keep") ) {
        sw_keep = 1 ;
      }
      else { /* consider it an input file */
        strcpy( fname_in , *argv ) ;
      }
      argv++ ;
  }
  if ( wrote_template ) exit(0) ;

  if ( !strcmp(fname_in,"") ) goto usage ;

#ifdef FUTURE
  gen_io_boilerplate() ;  /* 20091213 jm.  Generate the io_boilerplate_temporary.inc file */
#endif

  fprintf(stderr,"input file: %s\n",fname_in);

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
  if ( pre_parse( dir, fp_in, fp_tmp, 0 ) ) {
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

  if (sw_shownodes) {
    fprintf(stderr,"--- ModNames ---\n") ;
    show_nodelist(ModNames) ;
    fprintf(stderr,"--- Done ---\n") ;
  }

  gen_module_files( OutDir, thisprog_ver);

cleanup:
   if ( ! sw_keep ) {
#ifdef _WIN32
     sprintf(command,"del /F /Q %s\n",fname_tmp );
#else
     sprintf(command,"/bin/rm -f %s\n",fname_tmp );
#endif
     system( command ) ;
   }

   exit( 0 ) ;

}
#include "Template_data.c"
#include "Template_registry.c"

void
output_template( char * sw_modname_subst, char * sw_modnickname_subst, int force, int sw  ) // sw = 0, template; 1 = registry
{
    char ** p ;
    FILE *fp ;
    char fname[NAMELEN] ;
    char tmp1[2096], tmp2[2096], tmp3[2096] ;
    if ( sw == 0 ) { sprintf(fname,"%s.f90",sw_modname_subst) ; }
    else           { sprintf(fname,"%s_Registry.txt",sw_modname_subst) ; }

    if ( ! force ) { // check if file exists by trying to open file for reading. If the read is successful, exit program:
      if ( (fp = fopen( fname,"r" )) != NULL ) {
        fprintf(stderr,"Registry exiting. Attempt to overwrite file (%s) . Move out of the way or specify -force before -template option. \n", fname) ;
        exit(1) ;
      }
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

void
substitute( char * str , char * match , char * replace, char * result )
{
   char * p, *q ;
   char allup[NAMELEN], alllo[NAMELEN] ;
   size_t n, m ;
   int nmatch = 0 ;

   n = strlen( replace ) ;
   m = strlen( match ) ;
   strcpy(allup,replace) ; make_upper_case(allup) ;
   strcpy(alllo,replace) ; make_lower_case(alllo) ;
// watch for #defines, in which case first sub should be all upper, next all lower
   if ( str[0] == '#' ) {
     for ( p = str ; *p ; p++ ) {
       if ( matches( p, "define" ) ) nmatch = 2 ;
     }
   }

   for ( p = str , q = result ; *p ; )
   {
      if ( matches( p, match ) )
      {
        if        ( nmatch == 2 ) {
          strncpy( q, replace, n ) ;
          nmatch-- ;
        } else if ( nmatch == 1 ) {
          strncpy( q, alllo, n ) ;
          nmatch-- ;
        } else {
          strncpy( q, replace, n ) ;
        }
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
   int n  ;

   for ( n = 0, p = str, q = match ;  (*p && *q) ; p++, q++, n++ )
   {
     if ( *p != *q ) return(0) ;
   }
   if ( n != strlen(match) ) return(0) ;
   return(1) ;
}
