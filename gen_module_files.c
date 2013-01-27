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
  char tmp[NAMELEN], addnick[NAMELEN], nonick[NAMELEN] ;
  node_t *q, * r ;
  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
  fprintf(fp," SUBROUTINE %s_Copy%s( Src%sData, Dst%sData, CtrlCode, ErrStat, ErrMsg )\n",ModName->nickname,nonick,nonick,nonick ) ;
  fprintf(fp,"  TYPE(%s), INTENT(IN   ) :: Src%sData\n",addnick,nonick) ;
  fprintf(fp,"  TYPE(%s), INTENT(  OUT) :: Dst%sData\n",addnick,nonick) ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode\n") ; 
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ; 
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ; 
  fprintf(fp,"! Local \n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: i,j,k\n") ;
  fprintf(fp,"! \n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;

//  sprintf(tmp,"%s_%s",ModName->nickname,inoutlong) ;
//  sprintf(tmp,"%s",inoutlong) ;
  sprintf(tmp,"%s",addnick) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL ) 
  {
    fprintf(stderr,"Registry warning: generating %s_Copy%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
  } else {
    for ( r = q->fields ; r ; r = r->next )
    { 
      if ( r->type != NULL ) {
        if ( !strcmp( r->type->name, "meshtype" ) ) {
          if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(Src%sData%%%s)\n",nonick,r->name  ) ;
          } else if ( r->ndims > 0 ) {
            fprintf(stderr,"Registry warning: %s, %s: Mesh elements of a structure can be only 0 or 1D\n",ModName->name,r->name) ;
          }
          fprintf(fp,"  CALL MeshCopy( Src%sData%%%s%s, Dst%sData%%%s%s, CtrlCode, ErrStat, ErrMsg )\n",nonick,r->name,(r->ndims==1)?"(i)":"",nonick,r->name,(r->ndims==1)?"(i)":"") ;
          if ( r->ndims == 1 ) {
  fprintf(fp,"ENDDO\n") ;
          }

        } else if ( r->type->type_type == DERIVED && ! r->type->usefrom ) {
          char nonick2[NAMELEN] ;
          remove_nickname(ModName->nickname,r->type->name,nonick2) ;
          if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(Src%sData%%%s)\n",nonick,r->name  ) ;
          }
  fprintf(fp,"  CALL %s_Copy%s( Src%sData%%%s%s, Dst%sData%%%s%s, CtrlCode, ErrStat, ErrMsg )\n",
                                ModName->nickname,fast_interface_type_shortname(nonick2),nonick,r->name,(r->ndims==1)?"(i)":"",nonick,r->name,(r->ndims==1)?"(i)":"") ;
          if ( r->ndims == 1 ) {
  fprintf(fp,"ENDDO\n") ;
          }
        } else {
          fprintf(fp,"  Dst%sData%%%s = Src%sData%%%s\n",nonick,r->name,nonick,r->name) ;
        }
      } 
    }
  }

  fprintf(fp," END SUBROUTINE %s_Copy%s\n\n", ModName->nickname,nonick ) ;
  return(0) ;
}

int
gen_pack( FILE * fp, const node_t * ModName, char * inout, char *inoutlong )
{
  char tmp[NAMELEN], tmp2[NAMELEN], tmp3[NAMELEN], addnick[NAMELEN], nonick[NAMELEN] ;
  node_t *q, * r ;
  int frst ;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
//  sprintf(tmp,"%s_%s",ModName->nickname,inoutlong) ;
//  sprintf(tmp,"%s",inoutlong) ;
  sprintf(tmp,"%s",addnick) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Pack%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
    return(1) ;
  }

  fprintf(fp," SUBROUTINE %s_Pack%s( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )\n", ModName->nickname,nonick) ;
  fprintf(fp,"  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)\n") ;
  fprintf(fp,"  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)\n") ;
  fprintf(fp,"  TYPE(%s),  INTENT(IN   ) :: InData\n",addnick ) ;
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
    if ( r->type == NULL ) {
      fprintf(stderr,"Registry warning generating %_Pack%s: %s has no type.\n",ModName->nickname,nonick,r->name) ;
      return ; // EARLY RETURN
    } else {
      if ( !strcmp( r->type->name, "meshtype" ) || (r->type->type_type == DERIVED && ! r->type->usefrom ) ) {
  fprintf(fp,"  REAL(ReKi),     ALLOCATABLE :: Re_%s_Buf(:)\n",r->name) ;
  fprintf(fp,"  REAL(DbKi),     ALLOCATABLE :: Db_%s_Buf(:)\n",r->name) ;
  fprintf(fp,"  INTEGER(IntKi), ALLOCATABLE :: Int_%s_Buf(:)\n",r->name) ;
      }
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
  fprintf(fp,"  CALL MeshPack( InData%%%s%s, Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, ErrStat, ErrMsg, .TRUE. ) ! %s \n", 
                                 r->name,(r->ndims==1)?"(i)":"",r->name,  r->name,    r->name,                              r->name ) ;
  fprintf(fp,"  IF(ALLOCATED(Re_%s_Buf)) Re_BufSz  = Re_BufSz  + SIZE( Re_%s_Buf  ) ! %s\n",r->name,r->name,r->name ) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf)) Db_BufSz  = Db_BufSz  + SIZE( Db_%s_Buf  ) ! %s\n",r->name,r->name,r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf))Int_BufSz = Int_BufSz + SIZE( Int_%s_Buf ) ! %s\n",r->name,r->name,r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Re_%s_Buf))  DEALLOCATE(Re_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf))  DEALLOCATE(Db_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf)) DEALLOCATE(Int_%s_Buf)\n",r->name, r->name) ;
      if ( r->ndims == 1 ) {
  fprintf(fp,"ENDDO\n") ;
      }

    } else if ( r->type->type_type == DERIVED && ! r->type->usefrom ) {
      char nonick2[NAMELEN] ;
      remove_nickname(ModName->nickname,r->type->name,nonick2) ;
      if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(InData%%%s)\n",r->name  ) ;
      }
  fprintf(fp,"  CALL %s_Pack%s( Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, InData%%%s%s, ErrStat, ErrMsg, .TRUE. ) ! %s \n",
                        ModName->nickname,fast_interface_type_shortname(nonick2), r->name, r->name, r->name, r->name,
                        (r->ndims==1)?"(i)":"", r->name ) ;
  fprintf(fp,"  IF(ALLOCATED(Re_%s_Buf)) Re_BufSz  = Re_BufSz  + SIZE( Re_%s_Buf  ) ! %s\n",r->name,r->name,r->name ) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf)) Db_BufSz  = Db_BufSz  + SIZE( Db_%s_Buf  ) ! %s\n",r->name,r->name,r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf))Int_BufSz = Int_BufSz + SIZE( Int_%s_Buf ) ! %s\n",r->name,r->name,r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Re_%s_Buf))  DEALLOCATE(Re_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Db_%s_Buf))  DEALLOCATE(Db_%s_Buf)\n",r->name, r->name) ;
  fprintf(fp,"  IF(ALLOCATED(Int_%s_Buf)) DEALLOCATE(Int_%s_Buf)\n",r->name, r->name) ;
      if ( r->ndims == 1 ) {
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

      if ( r->ndims == 1 ) {
  fprintf(fp,"ENDDO\n") ;
      }

    } else if ( r->type->type_type == DERIVED && ! r->type->usefrom ) {
      char nonick2[NAMELEN] ;
      remove_nickname(ModName->nickname,r->type->name,nonick2) ;
      if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(InData%%%s)\n",r->name  ) ;
      }
  fprintf(fp,"  CALL %s_Pack%s( Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, InData%%%s%s, ErrStat, ErrMsg, OnlySize ) ! %s \n",
                        ModName->nickname,fast_interface_type_shortname(nonick2), r->name, r->name, r->name, r->name, 
                        (r->ndims==1)?"(i)":"",
                        r->name ) ;
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
      if ( r->ndims == 1 ) {
  fprintf(fp,"ENDDO\n") ;
      }

    } else  {
      char * indent ;
      sprintf(tmp2,"SIZE(InData%%%s)",r->name) ;
      if        ( r->ndims==0 ) { 
        strcpy(tmp3,"") ; 
      } else if ( r->ndims==1 ) { 
        strcpy(tmp3,"") ; 
      } else if ( r->ndims==2 ) { 
        sprintf(tmp3,"(1:(%s),1)",tmp2) ; 
      } else if ( r->ndims==3 ) { 
        sprintf(tmp3,"(1:(%s),1,1)",tmp2) ; 
      } else if ( r->ndims==4 ) { 
        sprintf(tmp3,"(1:(%s),1,1,1)",tmp2) ; 
      } else if ( r->ndims==5 ) { 
        sprintf(tmp3,"(1:(%s),1,1,1,1)",tmp2) ; 
      } else                    { 
        fprintf(stderr,"Registry WARNING: too many dimensions for %s\n",r->name) ; 
      }
      indent = "" ;
      if ( !strcmp( r->type->mapsto, "REAL(ReKi)") ||
           !strcmp( r->type->mapsto, "REAL(DbKi)") ||
           !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
        if ( r->ndims > 0 && has_deferred_dim( r, 0 )) {
  fprintf(fp,"  IF ( ALLOCATED(InData%%%s) ) THEN\n", r->name ) ;
          indent = "  " ;
        }
        if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
  fprintf(fp,"  %sIF ( .NOT. OnlySize ) ReKiBuf ( Re_Xferred:Re_Xferred+(%s)-1 ) =  %s(InData%%%s %s)\n",
             indent,(r->ndims>0)?tmp2:"1",(r->ndims>0)?"PACK":"",r->name,(r->ndims>0)?",.TRUE.":"") ;
  fprintf(fp,"  %sRe_Xferred   = Re_Xferred   + %s\n",indent,(r->ndims>0)?tmp2:"1"  ) ;
        }
        else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
  fprintf(fp,"  %sIF ( .NOT. OnlySize ) DbKiBuf ( Db_Xferred:Db_Xferred+(%s)-1 ) =  %s(InData%%%s %s)\n",
             indent,(r->ndims>0)?tmp2:"1",(r->ndims>0)?"PACK":"",r->name,(r->ndims>0)?",.TRUE.":"") ;
  fprintf(fp,"  %sDb_Xferred   = Db_Xferred   + %s\n",indent,(r->ndims>0)?tmp2:"1"  ) ;
        }
        else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
  fprintf(fp,"  %sIF ( .NOT. OnlySize ) IntKiBuf ( Int_Xferred:Int_Xferred+(%s)-1 ) = %s(InData%%%s %s)\n",
             indent,(r->ndims>0)?tmp2:"1",(r->ndims>0)?"PACK":"",r->name,(r->ndims>0)?",.TRUE.":"") ;
  fprintf(fp,"  %sInt_Xferred   = Int_Xferred   + %s\n",indent,(r->ndims>0)?tmp2:"1"  ) ;
        }
        if ( r->ndims > 0 && has_deferred_dim( r, 0 )) {
  fprintf(fp,"  ENDIF\n") ;
        }
      }
    }
  }

  fprintf(fp," END SUBROUTINE %s_Pack%s\n\n", ModName->nickname,nonick ) ;
  return(0) ;
}

int
gen_unpack( FILE * fp, const node_t * ModName, char * inout, char * inoutlong )
{
  char tmp[NAMELEN], tmp2[NAMELEN], tmp3[NAMELEN], tmp4[NAMELEN], addnick[NAMELEN], nonick[NAMELEN] ;
  node_t *q, * r ;
  int idim, frst ;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
//  sprintf(tmp,"%s_%s",ModName->nickname,inoutlong) ;
//  sprintf(tmp,"%s",inoutlong) ;
  sprintf(tmp,"%s",addnick) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Unpack%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
    return(1) ;
  }

  fprintf(fp," SUBROUTINE %s_Unpack%s( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )\n", ModName->nickname,nonick ) ;
  fprintf(fp,"  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)\n") ;
  fprintf(fp,"  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)\n") ;
  fprintf(fp,"  TYPE(%s), INTENT(  OUT) :: OutData\n",addnick ) ;
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
  fprintf(fp,"  LOGICAL, ALLOCATABLE           :: mask1(:)\n") ;
  fprintf(fp,"  LOGICAL, ALLOCATABLE           :: mask2(:,:)\n") ;
  fprintf(fp,"  LOGICAL, ALLOCATABLE           :: mask3(:,:,:)\n") ;
  fprintf(fp,"  LOGICAL, ALLOCATABLE           :: mask4(:,:,:,:)\n") ;
  fprintf(fp,"  LOGICAL, ALLOCATABLE           :: mask5(:,:,:,:,:)\n") ;

  fprintf(fp," ! buffers to store meshes, if any\n") ;
  for ( r = q->fields ; r ; r = r->next ) 
  {
    if ( r->type == NULL ) {
      fprintf(stderr,"Registry warning generating %_Unpack%s: %s has no type.\n",ModName->nickname,nonick,r->name) ;
      return ; // EARLY RETURN
    } else {
      if ( !strcmp( r->type->name, "meshtype" ) || (r->type->type_type == DERIVED && ! r->type->usefrom ) ) {
  fprintf(fp,"  REAL(ReKi),    ALLOCATABLE :: Re_%s_Buf(:)\n",r->name) ;
  fprintf(fp,"  REAL(DbKi),    ALLOCATABLE :: Db_%s_Buf(:)\n",r->name) ;
  fprintf(fp,"  INTEGER(IntKi),    ALLOCATABLE :: Int_%s_Buf(:)\n",r->name) ;
      }
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
  fprintf(fp,"  CALL MeshPack( OutData%%%s%s, Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, ErrStat, ErrMsg , .TRUE. ) ! %s \n",
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

    } else if ( r->type->type_type == DERIVED && ! r->type->usefrom ) {
      char nonick2[NAMELEN] ;
      remove_nickname(ModName->nickname,r->type->name,nonick2) ;
      if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(OutData%%%s)\n",r->name  ) ;
      }
  fprintf(fp," ! first call %s_Pack%s to get correctly sized buffers for unpacking\n",
                        ModName->nickname,fast_interface_type_shortname(nonick2)) ;
  fprintf(fp,"  CALL %s_Pack%s( Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, OutData%%%s%s, ErrStat, ErrMsg, .TRUE. ) ! %s \n",
                        ModName->nickname,fast_interface_type_shortname(nonick2), r->name, r->name, r->name, r->name, (r->ndims==1)?"(i)":"",r->name ) ;
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
  fprintf(fp,"  CALL %s_Unpack%s( Re_%s_Buf, Db_%s_Buf, Int_%s_Buf, OutData%%%s%s, ErrStat, ErrMsg ) ! %s \n",
                        ModName->nickname,fast_interface_type_shortname(nonick2), r->name, r->name, r->name, r->name, 
                        (r->ndims==1)?"(i)":"",
                        r->name ) ;
      if ( r->ndims == 1 ) {
  fprintf(fp,"ENDDO\n") ;
      }

    } else  {
      char * indent ;
      char arrayname[NAMELEN] ;
 
      sprintf(arrayname,"OutData%%%s",r->name) ;
      sprintf(tmp2,"SIZE(OutData%%%s)",r->name) ;
      if      ( r->ndims==0 ) { strcpy(tmp3,"") ; }
      else if ( r->ndims==1 ) { strcpy(tmp3,"") ; }
      else if ( r->ndims==2 ) { sprintf(tmp3,"(1:(%s),1)",tmp2) ; }
      else if ( r->ndims==3 ) { sprintf(tmp3,"(1:(%s),1,1)",tmp2) ; }
      else if ( r->ndims==4 ) { sprintf(tmp3,"(1:(%s),1,1,1)",tmp2) ; }
      else if ( r->ndims==5 ) { sprintf(tmp3,"(1:(%s),1,1,1,1)",tmp2) ; }
      else                    { fprintf(stderr,"Registry WARNING: too many dimensions for %s\n",r->name) ; }
      indent = "" ;
      if ( !strcmp( r->type->mapsto, "REAL(ReKi)") ||
           !strcmp( r->type->mapsto, "REAL(DbKi)") ||
           !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
        if ( r->ndims > 0 && has_deferred_dim( r, 0 )) {
  fprintf(fp,"  IF ( ALLOCATED(OutData%%%s) ) THEN\n", r->name ) ;
          indent = "  " ;
        }
        if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")     ) {
          if ( r->ndims > 0 ) { sprintf(tmp4,"Re_Xferred:Re_Xferred+(%s)-1",(r->ndims>0)?tmp2:"1") ; }
          else                { sprintf(tmp4,"Re_Xferred") ; }
          if ( r->ndims > 0 ) {
            gen_mask_alloc(fp, r->ndims, arrayname ) ;
  fprintf(fp,"  %sOutData%%%s = UNPACK(ReKiBuf( %s ),mask%d,OutData%%%s)\n",indent,r->name,tmp4,r->ndims,r->name) ;
  fprintf(fp,"  DEALLOCATE(mask%d)\n",r->ndims) ;
          } else {
  fprintf(fp,"  %sOutData%%%s%s = ReKiBuf ( %s )\n",indent,r->name,tmp3,tmp4) ;
          }
  fprintf(fp,"  %sRe_Xferred   = Re_Xferred   + %s\n",indent,(r->ndims>0)?tmp2:"1"  ) ;
        }
        else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")     ) {
          if ( r->ndims > 0 ) { sprintf(tmp4,"Db_Xferred:Re_Xferred+(%s)-1",(r->ndims>0)?tmp2:"1") ; }
          else                { sprintf(tmp4,"Db_Xferred") ; }
  fprintf(fp,"  %sOutData%%%s%s = DbKiBuf ( %s )\n",indent,r->name,tmp3,tmp4) ;
  fprintf(fp,"  %sDb_Xferred   = Db_Xferred   + %s\n",indent,(r->ndims>0)?tmp2:"1"  ) ;
        }
        else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ) {
          if ( r->ndims > 0 ) { sprintf(tmp4,"Int_Xferred:Re_Xferred+(%s)-1",(r->ndims>0)?tmp2:"1") ; }
          else                { sprintf(tmp4,"Int_Xferred") ; }
  fprintf(fp,"  %sOutData%%%s%s = IntKiBuf ( %s )\n",indent,r->name,tmp3,tmp4) ;
  fprintf(fp,"  %sInt_Xferred   = Int_Xferred   + %s\n",indent,(r->ndims>0)?tmp2:"1"  ) ;
        }
        if ( r->ndims > 0 && has_deferred_dim( r, 0 )) {
  fprintf(fp,"  ENDIF\n") ;
        }
      }
    }
  }
  fprintf(fp,"  Re_Xferred   = Re_Xferred-1\n") ;
  fprintf(fp,"  Db_Xferred   = Db_Xferred-1\n") ;
  fprintf(fp,"  Int_Xferred  = Int_Xferred-1\n") ;
  fprintf(fp," END SUBROUTINE %s_Unpack%s\n\n", ModName->nickname,nonick ) ;
  return(0) ;
}

gen_mask_alloc( FILE *fp, int ndims, char *tmp )
{
  if        ( ndims == 1 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1))); mask%d = .TRUE.\n",ndims,tmp,ndims) ;
  } else if ( ndims == 2 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1),SIZE(%s,2))); mask%d = .TRUE.\n",ndims,tmp,tmp,ndims) ;
  } else if ( ndims == 3 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1),SIZE(%s,2),SIZE(%s,3))); mask%d = .TRUE.\n",ndims,tmp,tmp,tmp,ndims) ;
  } else if ( ndims == 4 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1),SIZE(%s,2),SIZE(%s,3),SIZE(%s,4))); mask%d = .TRUE.\n",ndims,tmp,tmp,tmp,tmp,ndims) ;
  } else if ( ndims == 5 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1),SIZE(%s,2),SIZE(%s,3),SIZE(%s,4),SIZE(%s,5))); mask%d = .TRUE.\n",ndims,tmp,tmp,tmp,tmp,tmp,ndims) ;
  }
}



int
gen_destroy( FILE * fp, const node_t * ModName, char * inout, char * inoutlong )
{
  char tmp[NAMELEN], addnick[NAMELEN], nonick[NAMELEN] ;
  node_t *q, * r ;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
  fprintf(fp," SUBROUTINE %s_Destroy%s( %sData, ErrStat, ErrMsg )\n",ModName->nickname,nonick,nonick );
  fprintf(fp,"  TYPE(%s), INTENT(INOUT) :: %sData\n",addnick,nonick) ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: i\n") ;
  fprintf(fp,"! \n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;

//  sprintf(tmp,"%s_%s",ModName->nickname,inoutlong) ;
//  sprintf(tmp,"%s",inoutlong) ;
  sprintf(tmp,"%s",addnick) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Destroy%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
  } else {
    for ( r = q->fields ; r ; r = r->next )
    {
      if ( r->type == NULL ) {
        fprintf(stderr,"Registry warning generating %_Destroy%s: %s has no type.\n",ModName->nickname,nonick,r->name) ;
      } else {
        if ( !strcmp( r->type->name, "meshtype" ) ) {
          if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(%sData%%%s)\n",nonick,r->name  ) ;
          } else if ( r->ndims > 0 ) {
            fprintf(stderr,"Registry warning: %s, %s: Mesh elements of a structure can be only 0 or 1D\n",ModName->name,r->name) ;
          }
          fprintf(fp,"  CALL MeshDestroy( %sData%%%s%s, ErrStat, ErrMsg )\n",nonick,r->name,(r->ndims==1)?"(i)":"") ;
          if ( r->ndims > 0 ) {
  fprintf(fp,"ENDDO\n") ;
          }
        } else if ( r->type->type_type == DERIVED && ! r->type->usefrom ) {
          char nonick2[NAMELEN] ;
          remove_nickname(ModName->nickname,r->type->name,nonick2) ;
          if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(%sData%%%s)\n",nonick,r->name  ) ;
          }
          fprintf(fp,"  CALL %s_Destroy%s( %sData%%%s%s, ErrStat, ErrMsg )\n",
                          ModName->nickname,fast_interface_type_shortname(nonick2),nonick,r->name,(r->ndims==1)?"(i)":"") ;
          if ( r->ndims == 1 ) {
  fprintf(fp,"ENDDO\n") ;
          }
        } else if ( r->ndims > 0 ) {
          if ( r->dims[0]->deferred )     // if one dim is they all have to be; see check in type.c
          {
            if ( r->ndims == 1 ) {
  fprintf(fp,"DO i = 1, SIZE(%sData%%%s)\n",nonick,r->name  ) ;
            }
            fprintf(fp,"  DEALLOCATE(%sData%%%s)\n",nonick,r->name) ;
            if ( r->ndims == 1 ) {
  fprintf(fp,"ENDDO\n") ;
            }
          }
        }
      }
    }
  }

  fprintf(fp," END SUBROUTINE %s_Destroy%s\n\n", ModName->nickname,nonick ) ;
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
  fprintf(fp,"  REAL(ReKi), ALLOCATABLE,      INTENT(  OUT) :: Re_RetAry(:)\n") ;
  fprintf(fp,"  REAL(DbKi), ALLOCATABLE,      INTENT(  OUT) :: Db_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi), ALLOCATABLE,  INTENT(  OUT) :: Int_RetAry(:)\n") ;
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
  fprintf(fp,"  REAL(ReKi), ALLOCATABLE,      INTENT(IN   ) :: Re_RetAry(:)\n") ;
  fprintf(fp,"  REAL(DbKi), ALLOCATABLE,      INTENT(IN   ) :: Db_RetAry(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi), ALLOCATABLE,   INTENT(IN   ) :: Int_RetAry(:)\n") ;
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
gen_module( FILE * fp , node_t * ModName )
{
  node_t * q, * r ;
  int i ;

  if ( strlen(ModName->nickname) > 0 ) {
// gen preamble
    {
      char ** p ;
      for ( p = FAST_preamble ; *p ; p++ ) { fprintf( fp, *p, ModName->name ) ; }
    }

// generate each derived data type
    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      if ( q->usefrom == 0 ) {
        fprintf(fp,"  TYPE, PUBLIC :: %s\n",q->mapsto) ;
        for ( r = q->fields ; r ; r = r->next )
        { 
          if ( r->type != NULL ) {
            if ( r->type->type_type == DERIVED ) {
              fprintf(fp,"    TYPE(%s) ",r->type->mapsto ) ;
            } else {
              fprintf(fp,"    %s ",r->type->mapsto ) ;
            }
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
        }
        fprintf(fp,"  END TYPE %s\n",q->mapsto) ;
      }
    }
    fprintf(fp,"CONTAINS\n") ;
    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      if ( q->usefrom == 0 ) {

        char * ddtname, * ddtnamelong, nonick[NAMELEN] ;
        ddtname = q->name ;

        remove_nickname(ModName->nickname,ddtname,nonick) ;
        
//fprintf(stderr,">> %s %s %s \n",ModName->name, ddtname, nonick) ;

        if ( is_a_fast_interface_type( nonick ) ) {
          ddtnamelong = nonick ;
          ddtname = fast_interface_type_shortname( nonick ) ;
        } else {
          ddtnamelong = ddtname ;
        }

        gen_copy( fp, ModName, ddtname, ddtnamelong ) ;
        gen_destroy( fp, ModName, ddtname, ddtnamelong ) ;
        gen_pack( fp, ModName, ddtname, ddtnamelong ) ;
        gen_unpack( fp, ModName, ddtname, ddtnamelong ) ;
      }
    }

    gen_modname_pack( fp, ModName ) ;
    gen_modname_unpack( fp, ModName ) ;

    fprintf(fp,"END MODULE %s_Types\n",ModName->name ) ;
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
    if ( strlen( p->nickname ) > 0 ) {
      if ( strlen(dirname) > 0 ) 
        { sprintf(fname,"%s/%s_Types.f90",dirname,p->name) ; }
      else                       
        { sprintf(fname,"%s_Types.f90",p->name) ; }
      if ((fp = fopen( fname , "w" )) == NULL ) return(1) ;
      print_warning(fp,fname) ;
      gen_module ( fp , p ) ;
      close_the_file( fp ) ;
    }
  }
  return(0) ;
}

int
remove_nickname( char *nickname, char *src, char *dst )
{
  char tmp[NAMELEN];
  int n ;
  strcpy(tmp,make_lower_temp(nickname)) ;
  strcat(tmp,"_") ;
  n = strlen(tmp) ;
  if ( !strncmp(tmp,src,n) ) {
    strcpy(dst,&(src[n])) ;
  } else {
    strcpy(dst,src) ;
  }
}

int
append_nickname( char *nickname, char *src, char *dst )
{
  int n ;
  n = strlen(nickname) ;
  if ( n > 0 ) {
    sprintf(dst,"%s_%s",nickname,src) ;
  } else {
    strcpy(dst,src) ;
  }
}
