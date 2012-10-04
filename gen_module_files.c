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

  fprintf(fp," SUBROUTINE %s_Pack%s( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )\n", ModName->nickname,inout ) ;
  fprintf(fp,"  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf\n") ;
  fprintf(fp,"  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf\n") ;
  fprintf(fp,"  INTEGER(IntKi)),  ALLOCATABLE, INTENT(  OUT) :: IntKiBuf\n") ;
  fprintf(fp,"  TYPE(%s_%sType),  INTENT(IN   ) :: InData\n",ModName->nickname,inout ) ;
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
  fprintf(fp,"  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers\n") ;
  fprintf(fp," ! buffers to store meshes, if any\n") ;

  fprintf(fp,"  OnlySize = .FALSE.\n") ;
  fprintf(fp,"  IF ( PRESENT(SizeOnly) ) THEN\n") ;
  fprintf(fp,"    OnlySize = SizeOnly\n") ;
  fprintf(fp,"  ENDIF\n") ;

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
  fprintf(fp,"  Re_Xferred  = 0\n") ;
  fprintf(fp,"  Db_Xferred  = 0\n") ;
  fprintf(fp,"  Int_Xferred  = 0\n") ;

  frst = 1 ;
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  if ( frst == 1 ) { fprintf(fp," ! Allocate mesh buffers, if any (we'll also get sizes from these) \n") ; frst = 0 ;}
  fprintf(fp,"  CALL MeshPack( %sData%%%s, Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, ErrStat, ErrMess, OnlySize ) ! %s \n", 
                               inout,  r->name,r->name,  r->name,    r->name,                              r->name ) ;
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
  fprintf(fp,"    IF ( OnlySize ) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_%s_Buf)-1 ) = Re_%s_Buf\n",r->name,r->name,r->name) ;
  fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE(Re_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    IF ( OnlySize ) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_%s_Buf)-1 ) = Db_%s_Buf\n",r->name,r->name) ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE(Db_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf)) THEN\n",r->name) ;
  fprintf(fp,"    IF ( OnlySize ) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_%s_Buf)-1 ) = Int_%s_Buf\n",r->name,r->name) ;
  fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE(Int_%s_Buf)\n",r->name) ;
  fprintf(fp,"  ENDIF\n" ) ;
    } else  {
      sprintf(tmp2,"SIZE(%sData%%%s)\n",inout,r->name) ;
      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  IF ( .NOT. OnlySize ) ReKiBuf ( Re_Xferred ) =  %sData%%%s\n",inout,r->name) ;
  fprintf(fp,"  Re_Xferred   = Re_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  IF ( .NOT. OnlySize ) DbKiBuf ( Db_Xferred ) =  %sData%%%s\n",inout,r->name) ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  IF ( .NOT. OnlySize ) IntKiBuf ( Int_Xferred ) =  %sData%%%s\n",inout,r->name) ;
  fprintf(fp,"  Int_Xferred   = Int_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
    }
  }

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
gen_unpack( FILE * fp, const node_t * ModName, char * inout )
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

  fprintf(fp," SUBROUTINE %s_Unpack%s( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg )\n", ModName->nickname,inout ) ;
  fprintf(fp,"  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf\n") ;
  fprintf(fp,"  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf\n") ;
  fprintf(fp,"  INTEGER(IntKi)), ALLOCATABLE, INTENT(IN   ) :: IntKiBuf\n") ;
  fprintf(fp,"  TYPE(%s_%sType), INTENT(IN   ) :: InData\n",ModName->nickname,inout ) ;
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
  fprintf(fp," ! buffers to store meshes, if any\n") ;
  for ( r = q->fields ; r ; r = r->next ) 
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  fprintf(fp,"  REAL(ReKi),    POINTER :: Re_%s_Buf\n",r->name) ;
  fprintf(fp,"  REAL(ReKi),    POINTER :: Db_%s_Buf\n",r->name) ;
  fprintf(fp,"  REAL(ReKi),    POINTER :: Int_%s_Buf\n",r->name) ;
    }
  }
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
  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;

   // Unpack data
  for ( r = q->fields ; r ; r = r->next )
  {
    if ( !strcmp( r->type->name, "meshtype" ) ) {
  frst = 1 ;
  if (frst == 1) {fprintf(fp," ! first call MeshPack to get correctly sized buffers for unpacking\n") ;frst=0;}
  fprintf(fp,"  CALL MeshPack( %sData%%%s, Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, ErrStat, ErrMess , SizeOnly = .TRUE. ) ! %s \n",
                               inout,  r->name,r->name,  r->name,    r->name,                     r->name ) ;
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
  fprintf(fp,"  CALL MeshUnPack( %sData%%%s, Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, ErrStat, ErrMess ) ! %s \n",
                                 inout,  r->name,r->name,  r->name,    r->name,                     r->name ) ;
  fprintf(fp,"  IF( ASSOCIATED(Re_%s_Buf) DEALLOCATE(Re_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ASSOCIATED(Db_%s_Buf) DEALLOCATE(Db_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF( ASSOCIATED(Int_%s_Buf) DEALLOCATE(Int_%s_Buf)\n",r->name, r->name) ;
    } else  {
      sprintf(tmp2,"SIZE(%sData%%%s)\n",inout,r->name) ;
      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  %sData%%%s = ReKiBuf ( Re_Xferred )\n",inout,r->name) ;
  fprintf(fp,"  Re_Xferred   = Re_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  %sData%%%s = DbKiBuf ( Db_Xferred )\n",inout,r->name) ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  %sData%%%s = IntKiBuf ( Int_Xferred )\n",inout,r->name) ;
  fprintf(fp,"  Int_Xferred   = Int_Xferred   + %s\n",(r->ndims>0)?tmp2:"1"  ) ;
      }
    }
  }
  fprintf(fp,"  Re_Xferred   = Re_Xferred-1\n") ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred-1\n") ;
  fprintf(fp,"  Int_Xferred  = Int_Xferred-1\n") ;
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
gen_modname_pack( FILE *fp , const node_t * ModName )
{
  char tmp[NAMELEN] ;
  char *typenames[] = { "Input", "Param", "ContState", "DiscState", "ConstrState",
                       "OtherState", "Output", 0L } ;
  char **typename ;

  node_t *q, * r ;
  fprintf(fp," SUBROUTINE %s_Pack( Re_RetAry, Db_RetAry, IntRetAry,\n",ModName->nickname) ;
  fprintf(fp,"                     InData, ParamData, ContStateData, DiscStateData, &\n") ;
  fprintf(fp,"                     ConstrStateData, OtherStateData, OutData, ErrStat, ErrMsg, &\n" ) ;
  fprintf(fp,"                     SizeOnly )\n") ;
  fprintf(fp,"  TYPE(%s_InputType),      INTENT(IN   ) :: InData\n",          ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ParameterType),  INTENT(IN   ) :: ParamData\n",       ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousType), INTENT(IN   ) :: ContStateData\n",   ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_DiscreteType),   INTENT(IN   ) :: DiscStateData\n",   ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ConstraintType), INTENT(IN   ) :: ConstrStateData\n", ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_OtherStateType), INTENT(IN   ) :: OtherStateData\n",  ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_OutputType),     INTENT(IN   ) :: OutData\n",         ModName->nickname) ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,  INTENT(  OUT) :: Re_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,  INTENT(  OUT) :: Db_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,  INTENT(  OUT) :: Int_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),          INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),            INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"  LOGICAL, OPTIONAL,       INTENT(IN   ) :: SizeOnly\n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp,"    ! Local variables\n" ) ;
  fprintf(fp,"  REAL(ReKi), POINTER                    :: Re_Ary(:)\n") ;
  fprintf(fp,"  REAL(DbKi), POINTER                    :: Db_Ary(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi), POINTER                :: Int_Ary(:)\n") ;
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
  fprintf(fp,"  OnlySize = .FALSE.\n") ;
  fprintf(fp,"  IF ( PRESENT(SizeOnly) ) THEN\n") ;
  fprintf(fp,"    OnlySize = SizeOnly\n") ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;

  for ( typename = typenames ; *typename ; typename++ ) {
  fprintf(fp,"    ! Pack %s\n",*typename) ;
  fprintf(fp,"  NULLIFY(Re_Ary) ; NULLIFY(Db_Ary) ; NULLIFY( Int_Ary )\n") ;
  fprintf(fp,"  CALL %s_Pack%s(Re_Ary,Db_Ary,Int_Ary,ErrStat2,ErrMsg2,SizeOnly=.TRUE.)\n",ModName->nickname,*typename) ;
  fprintf(fp,"  IF ( ASSOCIATED( Re_Ary ) ) THEN\n") ;
  fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE( Re_Ary )\n") ;
  fprintf(fp,"    DEALLOCATE(Re_Ary)\n" ) ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ASSOCIATED( Db_Ary ) ) THEN\n") ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE( Db_Ary )\n") ;
  fprintf(fp,"    DEALLOCATE(Db_Ary)\n" ) ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ASSOCIATED( Int_Ary ) ) THEN\n") ;
  fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE( Int_Ary )\n") ;
  fprintf(fp,"    DEALLOCATE(Int_Ary)\n" ) ;
  fprintf(fp,"  ENDIF\n") ;
  }
  fprintf(fp,"  Re_Xferred  = Re_Xferred - 1\n") ;
  fprintf(fp,"  Db_Xferred  = Db_Xferred - 1\n") ;
  fprintf(fp,"  Int_Xferred  = Int_Xferred - 1\n") ;
  fprintf(fp,"  IF ( ASSOCIATED( Re_RetAry ) DEALLOCATE( Re_RetAry ) ;\n") ;
  fprintf(fp,"  IF ( Re_Xferred .GT. 0) ALLOCATE( Re_RetAry( Re_Xferred ) ) ;\n") ;
  fprintf(fp,"  IF ( ASSOCIATED( Db_RetAry ) DEALLOCATE( Db_RetAry ) ;\n") ;
  fprintf(fp,"  IF ( Db_Xferred .GT. 0) ALLOCATE( Db_RetAry( Db_Xferred ) ) ;\n") ;
  fprintf(fp,"  IF ( ASSOCIATED( Int_RetAry ) DEALLOCATE( Int_RetAry ) ;\n") ;
  fprintf(fp,"  IF ( Int_Xferred .GT. 0) ALLOCATE( Int_RetAry( Int_Xferred ) ) ;\n") ;

  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;

  for ( typename = typenames ; *typename ; typename++ ) {
  fprintf(fp,"    ! Pack %s\n",*typename) ;
  fprintf(fp,"  NULLIFY(Re_Ary) ; NULLIFY(Db_Ary) ; NULLIFY( Int_Ary )\n") ;
  fprintf(fp,"  CALL %s_Pack%s(Re_Ary,Db_Ary,Int_Ary,ErrStat2,ErrMsg2)\n",ModName->nickname,*typename) ;
  fprintf(fp,"  IF ( ASSOCIATED( Re_Ary ) ) THEN\n") ;
  fprintf(fp,"    IF ( .NOT. OnlySize ) Re_RetAry(Re_Xferred:Re_Xferred+SIZE(Re_Ary)-1)=Re_Ary\n") ;
  fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE( Re_Ary )\n") ;
  fprintf(fp,"    DEALLOCATE(Re_Ary)\n" ) ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ASSOCIATED( Db_Ary ) ) THEN\n") ;
  fprintf(fp,"    IF ( .NOT. OnlySize ) Db_RetAry(Db_Xferred:Db_Xferred+SIZE(Db_Ary)-1)=Db_Ary\n") ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE( Db_Ary )\n") ;
  fprintf(fp,"    DEALLOCATE(Db_Ary)\n" ) ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ASSOCIATED( Int_Ary ) ) THEN\n") ;
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
  char *typenames[] = { "Input", "Param", "ContState", "DiscState", "ConstrState",
                       "OtherState", "Output", 0L } ;
  char **typename ;

  node_t *q, * r ;
  fprintf(fp," SUBROUTINE %s_Pack( Re_RetAry, Db_RetAry, IntRetAry,\n",ModName->nickname) ;
  fprintf(fp,"                     InData, ParamData, ContStateData, DiscStateData, &\n") ;
  fprintf(fp,"                     ConstrStateData, OtherStateData, OutData, ErrStat, ErrMsg\n" ) ;
  fprintf(fp,"  TYPE(%s_InputType),      INTENT(  OUT) :: InData\n",          ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ParameterType),  INTENT(  OUT) :: ParamData\n",       ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousType), INTENT(  OUT) :: ContStateData\n",   ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_DiscreteType),   INTENT(  OUT) :: DiscStateData\n",   ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ConstraintType), INTENT(  OUT) :: ConstrStateData\n", ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_OtherStateType), INTENT(  OUT) :: OtherStateData\n",  ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_OutputType),     INTENT(  OUT) :: OutData\n",         ModName->nickname) ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,  INTENT(IN   ) :: Re_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,  INTENT(IN   ) :: Db_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(B1Ki), ALLOCATABLE,  INTENT(IN   ) :: Int_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp,"    ! Local variables\n" ) ;
  fprintf(fp,"  REAL(ReKi), POINTER                    :: Re_Ary(:)\n") ;
  fprintf(fp,"  REAL(DbKi), POINTER                    :: Db_Ary(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi), POINTER                :: Int_Ary(:)\n") ;
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


  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;
  for ( typename = typenames ; *typename ; typename++ ) {
  fprintf(fp,"    ! Unpack %s\n",*typename) ;
  fprintf(fp,"  NULLIFY(Re_Ary) ; NULLIFY(Db_Ary) ; NULLIFY( Int_Ary )\n") ;
  fprintf(fp,"  CALL %s_Pack%s(Re_Ary,Db_Ary,Int_Ary,ErrStat2,ErrMsg2,SizeOnly=.TRUE.)\n",ModName->nickname,*typename) ;
  fprintf(fp,"  IF ( ASSOCIATED( Re_Ary ) ) THEN\n") ;
  fprintf(fp,"    Re_Ary = Re_RetAry(Re_Xferred:Re_Xferred+SIZE(Re_Ary)-1)\n") ;
  fprintf(fp,"    Re_Xferred = Re_Xferred + SIZE( Re_Ary )\n") ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ASSOCIATED( Db_Ary ) ) THEN\n") ;
  fprintf(fp,"    DB_Ary = Db_RetAry(Db_Xferred:Db_Xferred+SIZE(Db_Ary)-1)\n") ;
  fprintf(fp,"    Db_Xferred = Db_Xferred + SIZE( Db_Ary )\n") ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  IF ( ASSOCIATED( Int_Ary ) ) THEN\n") ;
  fprintf(fp,"    Int_Ary = Int_RetAry(Int_Xferred:Int_Xferred+SIZE(Int_Ary)-1)\n") ;
  fprintf(fp,"    Int_Xferred = Int_Xferred + SIZE( Int_Ary )\n") ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"  CALL %s_Unpack%s(Re_Ary,Db_Ary,Int_Ary,ErrStat2,ErrMsg2)\n",ModName->nickname,*typename) ;
  fprintf(fp,"  IF ( ASSOCIATED( Re_Ary ) )  DEALLOCATE(Re_Ary)\n" ) ;
  fprintf(fp,"  IF ( ASSOCIATED( Db_Ary ) )  DEALLOCATE(Db_Ary)\n" ) ;
  fprintf(fp,"  IF ( ASSOCIATED( Int_Ary ) )  DEALLOCATE(Int_Ary)\n" ) ;
  }

  fprintf(fp,"  Re_Xferred   = Re_Xferred-1\n") ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred-1\n") ;
  fprintf(fp,"  Int_Xferred  = Int_Xferred-1\n") ;
  fprintf(fp," END SUBROUTINE %s_Pack\n\n", ModName->nickname ) ;
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
    fprintf(fp,"  END TYPE PUBLIC %s\n",q->mapsto) ;
  }
  fprintf(fp,"CONTAINS\n") ;

  gen_copy( fp, ModName, "Input" ) ;
  gen_destroy( fp, ModName, "Input" ) ;
  gen_pack( fp, ModName, "Input" ) ;
  gen_unpack( fp, ModName, "Input" ) ;

  gen_copy( fp, ModName, "Output" ) ;
  gen_destroy( fp, ModName, "Output" ) ;
  gen_pack( fp, ModName, "Output" ) ;
  gen_unpack( fp, ModName, "Output" ) ;

  gen_copy( fp, ModName, "ContinuousState" ) ;
  gen_destroy( fp, ModName, "ContinuousState" ) ;
  gen_pack( fp, ModName, "ContinuousState" ) ;
  gen_unpack( fp, ModName, "ContinuousState" ) ;

  gen_copy( fp, ModName, "DiscreteState" ) ;
  gen_destroy( fp, ModName, "DiscreteState" ) ;
  gen_pack( fp, ModName, "DiscreteState" ) ;
  gen_unpack( fp, ModName, "DiscreteState" ) ;

  gen_copy( fp, ModName, "OtherState" ) ;
  gen_destroy( fp, ModName, "OtherState" ) ;
  gen_pack( fp, ModName, "OtherState" ) ;
  gen_unpack( fp, ModName, "OtherState" ) ;

  gen_copy( fp, ModName, "Parameter" ) ;
  gen_destroy( fp, ModName, "Parameter" ) ;
  gen_pack( fp, ModName, "Parameter" ) ;
  gen_unpack( fp, ModName, "Parameter" ) ;

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



