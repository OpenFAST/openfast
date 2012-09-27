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
gen_copy( FILE * fp, const node_t * ModName, char * inout ) 
{
  char tmp[NAMELEN] ;
  node_t *q, * r ;
  fprintf(fp," SUBROUTINE %s_Copy%s( Src%sData, Dst%sData, CtrlCode, ErrStat, ErrMsg )\n", ModName->nickname,inout,inout,inout ) ;
  fprintf(fp,"  TYPE(%s_%sType), INTENT(IN   ) :: Src%sData\n",ModName->nickname,inout,inout) ;
  fprintf(fp,"  TYPE(%s_%sType), INTENT(  OUT) :: Dst%sData\n",ModName->nickname,inout,inout) ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode\n") ; 
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ; 
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ; 
  fprintf(fp,"! \n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;

  sprintf(tmp,"%s_%sType",ModName->nickname,inout) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL ) 
  {
    fprintf(stderr,"Registry warning: generating %s_Copy%s: cannot find definition for %s\n",ModName->nickname,inout,tmp) ;
  } else {
    for ( r = q->fields ; r ; r = r->next )
    { 
      if ( !strcmp( r->type->name, "meshtype" ) ) {
        fprintf(fp,"  CALL MeshCopy( Src%sData%%%s, Dst%sData%%%s, CtrlCode, ErrStat, ErMsg )\n",inout,r->name,inout,r->name,inout) ;
      } else {
        fprintf(fp,"  Dst%sData%%%s = Src%sData%%%s\n",inout,r->name,inout,r->name,inout) ;
      }
    }
  }

  fprintf(fp," END SUBROUTINE %s_Copy%s\n\n", ModName->nickname,inout ) ;
  return(0) ;
}

int
gen_pack( FILE * fp, const node_t * ModName, char * inout )
{
  char tmp[NAMELEN], tmp2[NAMELEN] ;
  node_t *q, * r ;
  int frst ;

  sprintf(tmp,"%s_%sType",ModName->nickname,inout) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Copy%s: cannot find definition for %s\n",ModName->nickname,inout,tmp) ;
    return(1) ;
  }

  fprintf(fp," SUBROUTINE %s_Pack%s( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg )\n", ModName->nickname,inout ) ;
  fprintf(fp,"  REAL(ReKi),      ALLOCATABLE, INTENT(  OUT) :: ReKiBuf\n") ;
  fprintf(fp,"  REAL(DbKi),      ALLOCATABLE, INTENT(  OUT) :: DbKiBuf\n") ;
  fprintf(fp,"  INTEGER(IntKi)), ALLOCATABLE, INTENT(  OUT) :: IntKiBuf\n") ;
  fprintf(fp,"  TYPE(%s_%sType), INTENT(IN   ) :: InData\n",ModName->nickname,inout ) ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"    ! Local variables\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_Xfered\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_Xfered\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_Xfered\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_CurrSz\n") ;
  fprintf(fp," ! buffers to store meshes, if any\n") ;

  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  fprintf(fp,"  REAL(ReKi),    POINTER :: Re_%s_Buf\n",r->name) ;
  fprintf(fp,"  REAL(ReKi),    POINTER :: Db_%s_Buf\n",r->name) ;
  fprintf(fp,"  REAL(ReKi),    POINTER :: Int_%s_Buf\n",r->name) ;
    }
  }

  fprintf(fp,"    !\n") ;

  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp,"  Re_Xfered  = 0\n") ;
  fprintf(fp,"  Db_Xfered  = 0\n") ;
  fprintf(fp,"  Int_Xfered  = 0\n") ;

  frst = 1 ;
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  if ( frst == 1 ) { fprintf(fp," ! Allocate mesh buffers, if any (we'll also get sizes from these) \n") ; frst = 0 ;}
  fprintf(fp,"  CALL MeshPack( %sData%%%s, Re_%s_Buf,  ErrStat, ErrMess ) ! %s \n", inout,r->name,r->name,r->name ) ;
  fprintf(fp,"  CALL MeshPack( %sData%%%s, Db_%s_Buf,  ErrStat, ErrMess ) ! %s \n", inout,r->name,r->name,r->name ) ;
  fprintf(fp,"  CALL MeshPack( %sData%%%s, Int_%s_Buf, ErrStat, ErrMess ) ! %s \n", inout,r->name,r->name,r->name ) ;
    }
  }

  fprintf(fp," ! computing sizes\n") ;

  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  fprintf(fp,"  Re_BufSz  = Re_BufSz  + SIZE( Re_%s_Buf  ) ! %s\n",r->name,r->name ) ;
  fprintf(fp,"  Db_BufSz  = Db_BufSz  + SIZE( Db_%s_Buf  ) ! %s\n",r->name,r->name ) ;
  fprintf(fp,"  Int_BufSz = Int_BufSz + SIZE( Int_%s_Buf ) ! %s\n",r->name,r->name ) ;
    } else if ( r->ndims == 0 ) {  // scalars
      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  Re_BufSz   = Re_BufSz   + 1  ! %s\n",r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  Db_BufSz   = Db_BufSz   + 1  ! %s\n",r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  Int_BufSz  = Int_BufSz  + 1  ! %s\n",r->name ) ;
      }
    } else { // r->ndims > 0
      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  Re_BufSz    = Re_BufSz    + SIZE( %sData%%%s )  ! %s \n", inout, r->name , r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  Db_BufSz    = Db_BufSz    + SIZE( %sData%%%s )  ! %s \n", inout, r->name , r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  Int_BufSz   = Int_BufSz   + SIZE( %sData%%%s )  ! %s \n", inout, r->name , r->name ) ;
      }
    }

  }

   // Allocate buffers
  fprintf(fp,"  IF ( Re_BufSz  .GT. 0 ) ALLOCATE( ReKiBuf(  Re_BufSz  ) )\n") ;
  fprintf(fp,"  IF ( Db_BufSz  .GT. 0 ) ALLOCATE( DbKiBuf(  Db_BufSz  ) )\n") ;
  fprintf(fp,"  IF ( Int_BufSz .GT. 0 ) ALLOCATE( IntKiBuf( Int_BufSz ) )\n") ;

   // Pack data
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  fprintf(fp,"  IF(ALLOCATED(Re_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_%s_Buf)-1 ) = Re_%s_Buf\n",r->name,r->name,r->name) ;
  fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE(Re_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_%s_Buf)-1 ) = Db_%s_Buf\n",r->name,r->name) ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE(Db_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_%s_Buf)-1 ) = Int_%s_Buf\n",r->name,r->name) ;
  fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE(Int_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
    } else  {
      sprintf(tmp2,"SIZE(%sData%%%s)\n",inout,r->name) ;
      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  ReKiBuf ( Re_Xferred ) =  %sData%%%s\n",inout,r->name) ;
  fprintf(fp,"  Re_Xferred   = Re_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  DbKiBuf ( Db_Xferred ) =  %sData%%%s\n",inout,r->name) ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  IntKiBuf ( Int_Xferred ) =  %sData%%%s\n",inout,r->name) ;
  fprintf(fp,"  Int_Xferred   = Int_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
    }
  }

  
#if 0
    for ( r = q->fields ; r ; r = r->next )
    {
      if ( !strcmp( r->type->name, "meshtype" ) ) {
        fprintf(fp,"  CALL MeshPack( Src%sData%%%s, Dst%sData%%%s, CtrlCode, ErrStat, ErMsg )\n",inout,r->name,inout,r->name,inout) ;
      } else {
        fprintf(fp,"  Dst%sData%%%s = Src%sData%%%s\n",inout,r->name,inout,r->name,inout) ;
      }
    }
#endif

  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  fprintf(fp,"  IF( ASSOCIATED(Re_%s_Buf) DEALLOCATE(Re_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ASSOCIATED(Db_%s_Buf) DEALLOCATE(Db_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ASSOCIATED(Int_%s_Buf) DEALLOCATE(Int_%s_Buf)\n",r->name, r->name) ;
    }
  }

  fprintf(fp," END SUBROUTINE %s_Pack%s\n\n", ModName->nickname,inout ) ;
  return(0) ;
}


int
gen_destroy( FILE * fp, const node_t * ModName, char * inout )
{
  char tmp[NAMELEN] ;
  node_t *q, * r ;
  fprintf(fp," SUBROUTINE %s_Destroy%s( %sData, ErrStat, ErrMsg )\n",ModName->nickname,inout );
  fprintf(fp,"  TYPE(%s_%sType), INTENT(IN   ) :: %sData\n",ModName->nickname,inout,inout) ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"! \n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;

  sprintf(tmp,"%s_%sType",ModName->nickname,inout) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Copy%s: cannot find definition for %s\n",ModName->nickname,inout,tmp) ;
  } else {
    for ( r = q->fields ; r ; r = r->next )
    {
      if ( !strcmp( r->type->name, "meshtype" ) ) {
        fprintf(fp,"  CALL MeshDestroy( %sData%%%s, ErrStat, ErMsg )\n",inout,r->name,inout) ;
      } else if ( r->ndims > 0 ) {
        if ( r->dims[0]->deferred )     // if one dim is they all have to be; see check in type.c
        {
          fprintf(fp,"  DEALLOCATE(%s)\n",r->name) ;
        }
      }
    }
  }

  fprintf(fp," END SUBROUTINE %s_Destroy%s\n\n", ModName->nickname,inout ) ;
  return(0) ;
}



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

  gen_copy( fp, ModName, "Input" ) ;
  gen_destroy( fp, ModName, "Input" ) ;
  gen_pack( fp, ModName, "Input" ) ;
  gen_copy( fp, ModName, "Output" ) ;
  gen_destroy( fp, ModName, "Output" ) ;
  gen_pack( fp, ModName, "Output" ) ;

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
      { sprintf(fname,"%s/%s.f90",dirname,p->name) ; }
    else                       
      { sprintf(fname,"%s.f90",p->name) ; }
    if ((fp = fopen( fname , "w" )) == NULL ) return(1) ;
    print_warning(fp,fname) ;
    gen_module ( fp , p ) ;
    close_the_file( fp ) ;
  }
  return(0) ;
}



