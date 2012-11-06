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
gen_copy( FILE * fp, const node_t * ModName, char * inout, char * inoutlong ) 
{
  char tmp[NAMELEN] ;
  node_t *q, * r ;
  fprintf(fp," SUBROUTINE %s_Copy%s( Src%sData, Dst%sData, CtrlCode, ErrStat, ErrMsg )\n", ModName->nickname,inout,inout,inout ) ;
  fprintf(fp,"  TYPE(%s_%sType), INTENT(IN   ) :: Src%sData\n",ModName->nickname,inoutlong,inout) ;
  fprintf(fp,"  TYPE(%s_%sType), INTENT(  OUT) :: Dst%sData\n",ModName->nickname,inoutlong,inout) ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode\n") ; 
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ; 
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ; 
  fprintf(fp,"! \n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;

  sprintf(tmp,"%s_%sType",ModName->nickname,inoutlong) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL ) 
  {
    fprintf(stderr,"Registry warning: generating %s_Copy%s: cannot find definition for %s\n",ModName->nickname,inout,tmp) ;
  } else {
    for ( r = q->fields ; r ; r = r->next )
    { 
      if ( !strcmp( r->type->name, "meshtype" ) ) {
        fprintf(fp,"  CALL MeshCopy( Src%sData%%%s, Dst%sData%%%s, CtrlCode, ErrStat, ErrMsg )\n",inout,r->name,inout,r->name,inout) ;
      } else {
        fprintf(fp,"  Dst%sData%%%s = Src%sData%%%s\n",inout,r->name,inout,r->name,inout) ;
      }
    }
  }

  fprintf(fp," END SUBROUTINE %s_Copy%s\n\n", ModName->nickname,inout ) ;
  return(0) ;
}

int
gen_pack( FILE * fp, const node_t * ModName, char * inout, char *inoutlong )
{
  char tmp[NAMELEN], tmp2[NAMELEN], tmp3[NAMELEN] ;
  node_t *q, * r ;
  int frst ;

  sprintf(tmp,"%s_%sType",ModName->nickname,inoutlong) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Copy%s: cannot find definition for %s\n",ModName->nickname,inout,tmp) ;
    return(1) ;
  }

  fprintf(fp," SUBROUTINE %s_Pack%s( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )\n", ModName->nickname,inout ) ;
  fprintf(fp,"  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)\n") ;
  fprintf(fp,"  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)\n") ;
  fprintf(fp,"  TYPE(%s_%sType),  INTENT(IN   ) :: InData\n",ModName->nickname,inoutlong ) ;
  fprintf(fp,"  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),     INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly\n") ;
  fprintf(fp,"    ! Local variables\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: i     \n") ;
  fprintf(fp,"  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers\n") ;
  fprintf(fp," ! buffers to store meshes, if any\n") ;

  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  fprintf(fp,"  REAL(ReKi),     ALLOCATABLE :: Re_%s_Buf(:)\n",r->name) ;
  fprintf(fp,"  REAL(DbKi),     ALLOCATABLE :: Db_%s_Buf(:)\n",r->name) ;
  fprintf(fp,"  INTEGER(IntKi), ALLOCATABLE :: Int_%s_Buf(:)\n",r->name) ;
    }
  }
  fprintf(fp,"  OnlySize = .FALSE.\n") ;
  fprintf(fp,"  IF ( PRESENT(SizeOnly) ) THEN\n") ;
  fprintf(fp,"    OnlySize = SizeOnly\n") ;
  fprintf(fp,"  ENDIF\n") ;

  fprintf(fp,"    !\n") ;

  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;
  fprintf(fp,"  Re_BufSz  = 0\n") ;
  fprintf(fp,"  Db_BufSz  = 0\n") ;
  fprintf(fp,"  Int_BufSz  = 0\n") ;

  frst = 1 ;
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {

      if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(InData%%%s)\n",r->name  ) ;
      } else if ( r->ndims > 0 ) {
        fprintf(stderr,"Registry warning: %s, %s: Mesh elements of a structure can be only 0 or 1D\n",ModName->name,r->name) ;
      }
  if ( frst == 1 ) { fprintf(fp," ! Allocate mesh buffers, if any (we'll also get sizes from these) \n") ; frst = 0 ;}
  fprintf(fp,"  CALL MeshPack( InData%%%s%s, Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, ErrStat, ErrMsg, OnlySize ) ! %s \n", 
                                 r->name,(r->ndims==1)?"(i)":"",r->name,  r->name,    r->name,                              r->name ) ;
  fprintf(fp,"  IF(ALLOCATED(Re_%s_Buf)) Re_BufSz  = Re_BufSz  + SIZE( Re_%s_Buf  ) ! %s\n",r->name,r->name,r->name ) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf)) Db_BufSz  = Db_BufSz  + SIZE( Db_%s_Buf  ) ! %s\n",r->name,r->name,r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf))Int_BufSz = Int_BufSz + SIZE( Int_%s_Buf ) ! %s\n",r->name,r->name,r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Re_%s_Buf))  DEALLOCATE(Re_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf))  DEALLOCATE(Db_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf)) DEALLOCATE(Int_%s_Buf)\n",r->name, r->name) ;
      if ( r->ndims > 0 ) {
  fprintf(fp,"ENDDO\n") ;
      }
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
  fprintf(fp,"  Re_BufSz    = Re_BufSz    + SIZE( InData%%%s )  ! %s \n", r->name , r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  Db_BufSz    = Db_BufSz    + SIZE( InData%%%s )  ! %s \n", r->name , r->name ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  Int_BufSz   = Int_BufSz   + SIZE( InData%%%s )  ! %s \n", r->name , r->name ) ;
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
      if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(InData%%%s)\n",r->name  ) ;
      } else if ( r->ndims > 0 ) {
        fprintf(stderr,"Registry warning: %s, %s: Mesh elements of a structure can be only 0 or 1D\n",ModName->name,r->name) ;
      }
  if ( frst == 1 ) { fprintf(fp," ! Allocate mesh buffers, if any (we'll also get sizes from these) \n") ; frst = 0 ;}
  fprintf(fp,"  CALL MeshPack( InData%%%s%s, Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, ErrStat, ErrMsg, OnlySize ) ! %s \n",
                                 r->name,(r->ndims==1)?"(i)":"",r->name,  r->name,    r->name,                              r->name ) ;
  fprintf(fp,"  IF(ALLOCATED(Re_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    IF ( .NOT. OnlySize ) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_%s_Buf)-1 ) = Re_%s_Buf\n",r->name,r->name,r->name) ;
  fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE(Re_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    IF ( .NOT. OnlySize ) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_%s_Buf)-1 ) = Db_%s_Buf\n",r->name,r->name) ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE(Db_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    IF ( .NOT. OnlySize ) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_%s_Buf)-1 ) = Int_%s_Buf\n",r->name,r->name) ;
  fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE(Int_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF( ALLOCATED(Re_%s_Buf) )  DEALLOCATE(Re_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ALLOCATED(Db_%s_Buf) )  DEALLOCATE(Db_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ALLOCATED(Int_%s_Buf) ) DEALLOCATE(Int_%s_Buf)\n",r->name, r->name) ;
      if ( r->ndims > 0 ) {
  fprintf(fp,"ENDDO\n") ;
      }
    } else  {
      sprintf(tmp2,"SIZE(InData%%%s)",r->name) ;
      if      ( r->ndims==0 ) { strcpy(tmp3,"") ; }
      else if ( r->ndims==1 ) { strcpy(tmp3,"") ; }
      else if ( r->ndims==2 ) { sprintf(tmp3,"(1:(%s),1)",tmp2) ; }
      else if ( r->ndims==3 ) { sprintf(tmp3,"(1:(%s),1,1)",tmp2) ; }
      else if ( r->ndims==4 ) { sprintf(tmp3,"(1:(%s),1,1,1)",tmp2) ; }
      else if ( r->ndims==5 ) { sprintf(tmp3,"(1:(%s),1,1,1,1)",tmp2) ; }
      else                    { fprintf(stderr,"Registry WARNING: too man dimensions for %s\n",r->name) ; }
      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  IF ( .NOT. OnlySize ) ReKiBuf ( Re_Xferred:Re_Xferred+(%s)-1 ) =  InData%%%s%s\n",(r->ndims>0)?tmp2:"1",r->name,tmp3) ;
  fprintf(fp,"  Re_Xferred   = Re_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  IF ( .NOT. OnlySize ) DbKiBuf ( Db_Xferred:Db_Xferred+(%s)-1 ) =  InData%%%s\n",(r->ndims>0)?tmp2:"1",r->name) ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  IF ( .NOT. OnlySize ) IntKiBuf ( Int_Xferred:Int_Xferred+(%s)-1 ) =  InData%%%s\n",(r->ndims>0)?tmp2:"1",r->name) ;
  fprintf(fp,"  Int_Xferred   = Int_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
    }
  }

  fprintf(fp," END SUBROUTINE %s_Pack%s\n\n", ModName->nickname,inout ) ;
  return(0) ;
}

int
gen_unpack( FILE * fp, const node_t * ModName, char * inout, char * inoutlong )
{
  char tmp[NAMELEN], tmp2[NAMELEN] ;
  node_t *q, * r ;
  int frst ;

  sprintf(tmp,"%s_%sType",ModName->nickname,inoutlong) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Copy%s: cannot find definition for %s\n",ModName->nickname,inout,tmp) ;
    return(1) ;
  }

  fprintf(fp," SUBROUTINE %s_Unpack%s( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )\n", ModName->nickname,inout ) ;
  fprintf(fp,"  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)\n") ;
  fprintf(fp,"  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)\n") ;
  fprintf(fp,"  TYPE(%s_%sType), INTENT(  OUT) :: OutData\n",ModName->nickname,inoutlong ) ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"    ! Local variables\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: i\n") ;

  fprintf(fp," ! buffers to store meshes, if any\n") ;
  for ( r = q->fields ; r ; r = r->next ) 
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  fprintf(fp,"  REAL(ReKi),    ALLOCATABLE :: Re_%s_Buf(:)\n",r->name) ;
  fprintf(fp,"  REAL(DbKi),    ALLOCATABLE :: Db_%s_Buf(:)\n",r->name) ;
  fprintf(fp,"  INTEGER(IntKi),    ALLOCATABLE :: Int_%s_Buf(:)\n",r->name) ;
    }
  }
  fprintf(fp,"    !\n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;
  fprintf(fp,"  Re_BufSz  = 0\n") ;
  fprintf(fp,"  Db_BufSz  = 0\n") ;
  fprintf(fp,"  Int_BufSz  = 0\n") ;

   // Unpack data
  frst = 1 ;
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
      if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(OutData%%%s)\n",r->name  ) ;
      } else if ( r->ndims > 0 ) {
        fprintf(stderr,"Registry warning: %s, %s: Mesh elements of a structure can be only 0 or 1D\n",ModName->name,r->name) ;
      }
  if (frst == 1) {fprintf(fp," ! first call MeshPack to get correctly sized buffers for unpacking\n") ;frst=0;}
  fprintf(fp,"  CALL MeshPack( OutData%%%s%s, Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, ErrStat, ErrMsg , SizeOnly = .TRUE. ) ! %s \n",
                               r->name,(r->ndims==1)?"(i)":"",r->name,  r->name,    r->name,                     r->name ) ;
  fprintf(fp,"  IF(ALLOCATED(Re_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    Re_%s_Buf = ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_%s_Buf)-1 )\n",r->name,r->name,r->name) ;
  fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE(Re_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    Db_%s_Buf = DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_%s_Buf)-1 )\n",r->name,r->name) ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE(Db_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    Int_%s_Buf = IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_%s_Buf)-1 )\n",r->name,r->name) ;
  fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE(Int_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  CALL MeshUnPack( OutData%%%s%s, Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, ErrStat, ErrMsg ) ! %s \n",
                                 r->name,(r->ndims==1)?"(i)":"",r->name,  r->name,    r->name,                     r->name ) ;
  fprintf(fp,"  IF( ALLOCATED(Re_%s_Buf) )  DEALLOCATE(Re_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ALLOCATED(Db_%s_Buf) )  DEALLOCATE(Db_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ALLOCATED(Int_%s_Buf) ) DEALLOCATE(Int_%s_Buf)\n",r->name, r->name) ;
      if ( r->ndims > 0 ) {
  fprintf(fp,"ENDDO\n") ;
      }
    } else  {
      sprintf(tmp2,"SIZE(OutData%%%s)\n",r->name) ;
      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  OutData%%%s = ReKiBuf ( Re_Xferred )\n",r->name) ;
  fprintf(fp,"  Re_Xferred   = Re_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  OutData%%%s = DbKiBuf ( Db_Xferred )\n",r->name) ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  OutData%%%s = IntKiBuf ( Int_Xferred )\n",r->name) ;
  fprintf(fp,"  Int_Xferred   = Int_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
    }
  }
  fprintf(fp,"  Re_Xferred   = Re_Xferred-1\n") ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred-1\n") ;
  fprintf(fp,"  Int_Xferred  = Int_Xferred-1\n") ;
  fprintf(fp," END SUBROUTINE %s_Unpack%s\n\n", ModName->nickname,inout ) ;
  return(0) ;
}


int
gen_destroy( FILE * fp, const node_t * ModName, char * inout, char * inoutlong )
{
  char tmp[NAMELEN] ;
  node_t *q, * r ;
  fprintf(fp," SUBROUTINE %s_Destroy%s( %sData, ErrStat, ErrMsg )\n",ModName->nickname,inout,inout );
  fprintf(fp,"  TYPE(%s_%sType), INTENT(INOUT) :: %sData\n",ModName->nickname,inoutlong,inout) ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: i\n") ;
  fprintf(fp,"! \n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;

  sprintf(tmp,"%s_%sType",ModName->nickname,inoutlong) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Copy%s: cannot find definition for %s\n",ModName->nickname,inout,tmp) ;
  } else {
    for ( r = q->fields ; r ; r = r->next )
    {
      if ( !strcmp( r->type->name, "meshtype" ) ) {
        if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(%sData%%%s)\n",inout,r->name  ) ;
        } else if ( r->ndims > 0 ) {
          fprintf(stderr,"Registry warning: %s, %s: Mesh elements of a structure can be only 0 or 1D\n",ModName->name,r->name) ;
        }
        fprintf(fp,"  CALL MeshDestroy( %sData%%%s%s, ErrStat, ErrMsg )\n",inout,r->name,(r->ndims==1)?"(i)":"") ;
        if ( r->ndims > 0 ) {
  fprintf(fp,"ENDDO\n") ;
        }
      } else if ( r->ndims > 0 ) {
        if ( r->dims[0]->deferred )     // if one dim is they all have to be; see check in type.c
        {
          fprintf(fp,"  DEALLOCATE(%sData%%%s)\n",inout,r->name) ;
        }
      }
    }
  }

  fprintf(fp," END SUBROUTINE %s_Destroy%s\n\n", ModName->nickname,inout ) ;
  return(0) ;
}

static char *typenames[] = { "Input", "Param", "ContState", "DiscState", "ConstrState",
                             "OtherState", "Output", 0L } ;
static char **typename ;
static char *argtypenames[] = { "InData", "ParamData", "ContStateData", "DiscStateData", "ConstrStateData",
                                "OtherStateData", "OutData", 0L } ;
static char **argtypename ;

int
gen_modname_pack( FILE *fp , const node_t * ModName )
{
  char tmp[NAMELEN] ;

  node_t *q, * r ;
  fprintf(fp," SUBROUTINE %s_Pack( Re_RetAry, Db_RetAry, Int_RetAry, &\n",ModName->nickname) ;
  fprintf(fp,"                     InData, ParamData, ContStateData, DiscStateData, &\n") ;
  fprintf(fp,"                     ConstrStateData, OtherStateData, OutData, ErrStat, ErrMsg, &\n" ) ;
  fprintf(fp,"                     SizeOnly )\n") ;
  fprintf(fp,"  TYPE(%s_InputType),           INTENT(IN   ) :: InData\n",          ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ParameterType),       INTENT(IN   ) :: ParamData\n",       ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousStateType), INTENT(IN   ) :: ContStateData\n",   ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_DiscreteStateType),   INTENT(IN   ) :: DiscStateData\n",   ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ConstraintStateType), INTENT(IN   ) :: ConstrStateData\n", ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_OtherStateType),      INTENT(IN   ) :: OtherStateData\n",  ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_OutputType),          INTENT(IN   ) :: OutData\n",         ModName->nickname) ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,   INTENT(  OUT) :: Re_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,   INTENT(  OUT) :: Db_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,   INTENT(  OUT) :: Int_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),               INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),                 INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"  LOGICAL, OPTIONAL,            INTENT(IN   ) :: SizeOnly\n") ;
  fprintf(fp,"    ! Local variables\n" ) ;
  fprintf(fp,"  REAL(ReKi), ALLOCATABLE                :: Re_Ary(:)\n") ;
  fprintf(fp,"  REAL(DbKi), ALLOCATABLE                :: Db_Ary(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi), ALLOCATABLE            :: Int_Ary(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Re_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Re_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Re_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Db_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Db_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Db_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Int_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Int_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Int_CurrSz\n") ;

  fprintf(fp,"  INTEGER(IntKi)                         :: ErrStat2\n") ;
  fprintf(fp,"  CHARACTER(Len(ErrMsg))                 :: ErrMsg2\n" )  ;
  fprintf(fp,"  LOGICAL                                :: OnlySize ! if present and true, do not pack, just allocate buffers\n") ;
  fprintf(fp,"    ! Executable statements\n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp,"  OnlySize = .FALSE.\n") ;
  fprintf(fp,"  IF ( PRESENT(SizeOnly) ) THEN\n") ;
  fprintf(fp,"    OnlySize = SizeOnly\n") ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;

  for ( typename = typenames, argtypename = argtypenames ; *typename ; typename++ , argtypename++ ) {
  fprintf(fp,"    ! Pack %s\n",*typename) ;
  fprintf(fp,"  IF ( ALLOCATED( Re_Ary ) )  DEALLOCATE(Re_Ary)\n" ) ;
  fprintf(fp,"  IF ( ALLOCATED( Db_Ary ) )  DEALLOCATE(Db_Ary)\n" ) ;
  fprintf(fp,"  IF ( ALLOCATED( Int_Ary ) )  DEALLOCATE(Int_Ary)\n" ) ;
  fprintf(fp,"  CALL %s_Pack%s(Re_Ary,Db_Ary,Int_Ary,%s,ErrStat2,ErrMsg2,SizeOnly=.TRUE.)\n",ModName->nickname,*typename,*argtypename) ;
  fprintf(fp,"  IF ( ALLOCATED( Re_Ary ) ) THEN\n") ;
  fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE( Re_Ary )\n") ;
  fprintf(fp,"    DEALLOCATE(Re_Ary)\n" ) ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ALLOCATED( Db_Ary ) ) THEN\n") ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE( Db_Ary )\n") ;
  fprintf(fp,"    DEALLOCATE(Db_Ary)\n" ) ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ALLOCATED( Int_Ary ) ) THEN\n") ;
  fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE( Int_Ary )\n") ;
  fprintf(fp,"    DEALLOCATE(Int_Ary)\n" ) ;
  fprintf(fp,"  ENDIF\n") ;
  }
  fprintf(fp,"  Re_Xferred  = Re_Xferred - 1\n") ;
  fprintf(fp,"  Db_Xferred  = Db_Xferred - 1\n") ;
  fprintf(fp,"  Int_Xferred  = Int_Xferred - 1\n") ;
  fprintf(fp,"  IF ( ALLOCATED( Re_RetAry ) ) DEALLOCATE( Re_RetAry ) ;\n") ;
  fprintf(fp,"  IF ( Re_Xferred .GT. 0) ALLOCATE( Re_RetAry( Re_Xferred ) ) ;\n") ;
  fprintf(fp,"  IF ( ALLOCATED( Db_RetAry ) ) DEALLOCATE( Db_RetAry ) ;\n") ;
  fprintf(fp,"  IF ( Db_Xferred .GT. 0) ALLOCATE( Db_RetAry( Db_Xferred ) ) ;\n") ;
  fprintf(fp,"  IF ( ALLOCATED( Int_RetAry ) ) DEALLOCATE( Int_RetAry ) ;\n") ;
  fprintf(fp,"  IF ( Int_Xferred .GT. 0) ALLOCATE( Int_RetAry( Int_Xferred ) ) ;\n") ;

  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;

  for ( typename = typenames, argtypename = argtypenames ; *typename ; typename++ , argtypename++ ) {
    fprintf(fp,"    ! Pack %s\n",*typename) ;
    fprintf(fp,"  IF ( ALLOCATED( Re_Ary ) )  DEALLOCATE(Re_Ary)\n" ) ;
    fprintf(fp,"  IF ( ALLOCATED( Db_Ary ) )  DEALLOCATE(Db_Ary)\n" ) ;
    fprintf(fp,"  IF ( ALLOCATED( Int_Ary ) )  DEALLOCATE(Int_Ary)\n" ) ;
    fprintf(fp,"  CALL %s_Pack%s(Re_Ary,Db_Ary,Int_Ary,%s,ErrStat2,ErrMsg2)\n",ModName->nickname,*typename,*argtypename) ;
    fprintf(fp,"  IF ( ALLOCATED( Re_Ary ) ) THEN\n") ;
    fprintf(fp,"    IF ( .NOT. OnlySize ) Re_RetAry(Re_Xferred:Re_Xferred+SIZE(Re_Ary)-1)=Re_Ary\n") ;
    fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE( Re_Ary )\n") ;
    fprintf(fp,"    DEALLOCATE(Re_Ary)\n" ) ;
    fprintf(fp,"  ENDIF\n") ;
    fprintf(fp,"  IF ( ALLOCATED( Db_Ary ) ) THEN\n") ;
    fprintf(fp,"    IF ( .NOT. OnlySize ) Db_RetAry(Db_Xferred:Db_Xferred+SIZE(Db_Ary)-1)=Db_Ary\n") ;
    fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE( Db_Ary )\n") ;
    fprintf(fp,"    DEALLOCATE(Db_Ary)\n" ) ;
    fprintf(fp,"  ENDIF\n") ;
    fprintf(fp,"  IF ( ALLOCATED( Int_Ary ) ) THEN\n") ;
    fprintf(fp,"    IF ( .NOT. OnlySize ) Int_RetAry(Int_Xferred:Int_Xferred+SIZE(Int_Ary)-1)=Int_Ary\n") ;
    fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE( Int_Ary )\n") ;
    fprintf(fp,"    DEALLOCATE(Int_Ary)\n" ) ;
    fprintf(fp,"  ENDIF\n") ;
  }

  fprintf(fp,"  Re_Xferred   = Re_Xferred - 1\n") ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred - 1\n") ;
  fprintf(fp,"  Int_Xferred  = Int_Xferred - 1\n") ;
  fprintf(fp," END SUBROUTINE %s_Pack\n\n", ModName->nickname ) ;
}

int
gen_modname_unpack( FILE *fp , const node_t * ModName )
{
  char tmp[NAMELEN] ;

  node_t *q, * r ;
  fprintf(fp," SUBROUTINE %s_Unpack( Re_RetAry, Db_RetAry, Int_RetAry, &\n",ModName->nickname) ;
  fprintf(fp,"                     InData, ParamData, ContStateData, DiscStateData, &\n") ;
  fprintf(fp,"                     ConstrStateData, OtherStateData, OutData, ErrStat, ErrMsg )\n" ) ;
  fprintf(fp,"  TYPE(%s_InputType),           INTENT(  OUT) :: InData\n",          ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ParameterType),       INTENT(  OUT) :: ParamData\n",       ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousStateType), INTENT(  OUT) :: ContStateData\n",   ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_DiscreteStateType),   INTENT(  OUT) :: DiscStateData\n",   ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ConstraintStateType), INTENT(  OUT) :: ConstrStateData\n", ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_OtherStateType),      INTENT(  OUT) :: OtherStateData\n",  ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_OutputType),          INTENT(  OUT) :: OutData\n",         ModName->nickname) ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,   INTENT(IN   ) :: Re_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,   INTENT(IN   ) :: Db_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,   INTENT(IN   ) :: Int_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"    ! Local variables\n" ) ;
  fprintf(fp,"  REAL(ReKi), ALLOCATABLE                :: Re_Ary(:)\n") ;
  fprintf(fp,"  REAL(DbKi), ALLOCATABLE                :: Db_Ary(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi), ALLOCATABLE            :: Int_Ary(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Re_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Re_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Re_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Db_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Db_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Db_CurrSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Int_BufSz\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Int_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                         :: Int_CurrSz\n") ;

  fprintf(fp,"  INTEGER(IntKi)                         :: ErrStat2\n") ;
  fprintf(fp,"  CHARACTER(Len(ErrMsg))                 :: ErrMsg2\n" )  ;


  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;
  for ( typename = typenames, argtypename = argtypenames ; *typename ; typename++ , argtypename++ ) {
  fprintf(fp,"    ! Unpack %s\n",*typename) ;
  fprintf(fp,"  IF ( ALLOCATED( Re_Ary ) )  DEALLOCATE(Re_Ary)\n" ) ;
  fprintf(fp,"  IF ( ALLOCATED( Db_Ary ) )  DEALLOCATE(Db_Ary)\n" ) ;
  fprintf(fp,"  IF ( ALLOCATED( Int_Ary ) )  DEALLOCATE(Int_Ary)\n" ) ;
  fprintf(fp,"  CALL %s_Pack%s(Re_Ary,Db_Ary,Int_Ary,%s,ErrStat2,ErrMsg2,SizeOnly=.TRUE.)\n",ModName->nickname,*typename,*argtypename) ;
  fprintf(fp,"  IF ( ALLOCATED( Re_Ary ) ) THEN\n") ;
  fprintf(fp,"    Re_Ary = Re_RetAry(Re_Xferred:Re_Xferred+SIZE(Re_Ary)-1)\n") ;
  fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE( Re_Ary )\n") ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ALLOCATED( Db_Ary ) ) THEN\n") ;
  fprintf(fp,"    DB_Ary = Db_RetAry(Db_Xferred:Db_Xferred+SIZE(Db_Ary)-1)\n") ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE( Db_Ary )\n") ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ALLOCATED( Int_Ary ) ) THEN\n") ;
  fprintf(fp,"    Int_Ary = Int_RetAry(Int_Xferred:Int_Xferred+SIZE(Int_Ary)-1)\n") ;
  fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE( Int_Ary )\n") ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  CALL %s_Unpack%s(Re_Ary,Db_Ary,Int_Ary,%s,ErrStat2,ErrMsg2)\n",ModName->nickname,*typename,*argtypename) ;
  fprintf(fp,"  IF ( ALLOCATED( Re_Ary ) )  DEALLOCATE(Re_Ary)\n" ) ;
  fprintf(fp,"  IF ( ALLOCATED( Db_Ary ) )  DEALLOCATE(Db_Ary)\n" ) ;
  fprintf(fp,"  IF ( ALLOCATED( Int_Ary ) )  DEALLOCATE(Int_Ary)\n" ) ;
  }

  fprintf(fp,"  Re_Xferred   = Re_Xferred-1\n") ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred-1\n") ;
  fprintf(fp,"  Int_Xferred  = Int_Xferred-1\n") ;
  fprintf(fp," END SUBROUTINE %s_Unpack\n\n", ModName->nickname ) ;
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
        if ( r->dims[0]->deferred )     // if one dim is deferred they all have to be; see check in type.c
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
      if ( r->ndims == 0 && strlen(r->inival) > 0 ) {
        fprintf(fp," :: %s = %s \n", r->name, r->inival ) ;
      } else {
        fprintf(fp," :: %s \n",r->name) ;
      }

    }
    fprintf(fp,"  END TYPE %s\n",q->mapsto) ;
  }
  fprintf(fp,"CONTAINS\n") ;

  gen_copy( fp, ModName, "Input", "Input" ) ;
  gen_destroy( fp, ModName, "Input", "Input" ) ;
  gen_pack( fp, ModName, "Input", "Input" ) ;
  gen_unpack( fp, ModName, "Input", "Input" ) ;

  gen_copy( fp, ModName, "Output", "Output" ) ;
  gen_destroy( fp, ModName, "Output", "Output" ) ;
  gen_pack( fp, ModName, "Output", "Output" ) ;
  gen_unpack( fp, ModName, "Output", "Output" ) ;

  gen_copy( fp, ModName, "ContState", "ContinuousState" ) ;
  gen_destroy( fp, ModName, "ContState", "ContinuousState" ) ;
  gen_pack( fp, ModName, "ContState", "ContinuousState" ) ;
  gen_unpack( fp, ModName, "ContState", "ContinuousState" ) ;

  gen_copy( fp, ModName, "DiscState", "DiscreteState" ) ;
  gen_destroy( fp, ModName, "DiscState", "DiscreteState" ) ;
  gen_pack( fp, ModName, "DiscState", "DiscreteState" ) ;
  gen_unpack( fp, ModName, "DiscState", "DiscreteState" ) ;

  gen_copy( fp, ModName, "ConstrState", "ConstraintState" ) ;
  gen_destroy( fp, ModName, "ConstrState", "ConstraintState" ) ;
  gen_pack( fp, ModName, "ConstrState", "ConstraintState" ) ;
  gen_unpack( fp, ModName, "ConstrState", "ConstraintState" ) ;

  gen_copy( fp, ModName, "OtherState", "OtherState" ) ;
  gen_destroy( fp, ModName, "OtherState", "OtherState" ) ;
  gen_pack( fp, ModName, "OtherState", "OtherState" ) ;
  gen_unpack( fp, ModName, "OtherState", "OtherState" ) ;

  gen_copy( fp, ModName, "Param", "Parameter" ) ;
  gen_destroy( fp, ModName, "Param", "Parameter" ) ;
  gen_pack( fp, ModName, "Param", "Parameter" ) ;
  gen_unpack( fp, ModName, "Param", "Parameter" ) ;

  gen_modname_pack( fp, ModName ) ;
  gen_modname_unpack( fp, ModName ) ;

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
      { sprintf(fname,"%s/%s_Types.f90",dirname,p->name) ; }
    else                       
      { sprintf(fname,"%s_Types.f90",p->name) ; }
    if ((fp = fopen( fname , "w" )) == NULL ) return(1) ;
    print_warning(fp,fname) ;
    gen_module ( fp , p ) ;
    close_the_file( fp ) ;
  }
  return(0) ;
}



