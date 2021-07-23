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

void gen_mask_alloc( FILE *fp, int ndims, char *tmp );

/**
 * ==============  Create the C2Farry Copy Subroutine in ModName_Types.f90 ======================
 *
 * In the C2F routines, we associate the pointer created in C with the variables in the
 * corresponding Fortran types.
 * ======================================================================================
 */
int
gen_copy_c2f( FILE         *fp        , // *.f90 file we are writting to
              const node_t *ModName   , // module name
              char         *inout     , // character string written out
              char         *inoutlong ) // not sure what this is used for
{
  node_t *q, *r ;
  char tmp[NAMELEN];
  char addnick[NAMELEN];
  char nonick[NAMELEN] ;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
  fprintf(fp," SUBROUTINE %s_C2Fary_Copy%s( %sData, ErrStat, ErrMsg, SkipPointers )\n", ModName->nickname, nonick,nonick );
  fprintf(fp,"    TYPE(%s), INTENT(INOUT) :: %sData\n"               , addnick, nonick                                      );
  fprintf(fp,"    INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n"                                                              );
  fprintf(fp,"    CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n"                                                               );
  fprintf(fp,"    LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers\n"                                                         );
  fprintf(fp,"    ! \n"                                                                                                     );
  fprintf(fp,"    LOGICAL                        :: SkipPointers_local\n");
  fprintf(fp,"    ErrStat = ErrID_None\n"                                                                                   );
  fprintf(fp,"    ErrMsg  = \"\"\n\n"                                                                                       );
  fprintf(fp,"    IF (PRESENT(SkipPointers)) THEN\n");
  fprintf(fp,"       SkipPointers_local = SkipPointers\n");
  fprintf(fp,"    ELSE\n");
  fprintf(fp,"       SkipPointers_local = .false.\n");
  fprintf(fp,"    END IF\n");

  sprintf(tmp,"%s",addnick) ;

  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_C2Fary_Copy%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
  } else {
    for ( r = q->fields ; r ; r = r->next )
    {
      if ( r->type != NULL ) {
        if ( r->type->type_type == DERIVED   ) { // && ! r->type->usefrom
          fprintf(stderr,"Registry WARNING: derived data type %s of type %s is not passed through C interface\n",r->name,r->type->name) ;
        } else {
            if ( is_pointer(r) ) {
                 fprintf(fp,"\n    ! -- %s %s Data fields\n",r->name,nonick) ;
                 fprintf(fp,"    IF ( .NOT. SkipPointers_local ) THEN\n");                 
                 fprintf(fp,"       IF ( .NOT. C_ASSOCIATED( %sData%%C_obj%%%s ) ) THEN\n",nonick,r->name) ;
                 fprintf(fp,"          NULLIFY( %sData%%%s )\n",nonick,r->name) ;
                 fprintf(fp,"       ELSE\n") ;
                 fprintf(fp,"          CALL C_F_POINTER(%sData%%C_obj%%%s, %sData%%%s, (/%sData%%C_obj%%%s_Len/))\n",nonick,r->name,nonick,r->name,nonick,r->name) ;
                 fprintf(fp,"       END IF\n") ;
                 fprintf(fp, "    END IF\n");
            }
            else if (!has_deferred_dim(r, 0)) {
               if (!strcmp(r->type->mapsto, "REAL(ReKi)") ||
                  !strcmp(r->type->mapsto, "REAL(SiKi)") ||
                  !strcmp(r->type->mapsto, "REAL(DbKi)") ||
                  !strcmp(r->type->mapsto, "REAL(R8Ki)") ||
                  !strcmp(r->type->mapsto, "INTEGER(IntKi)") ||
                  !strcmp(r->type->mapsto, "LOGICAL"))
               {
                  fprintf(fp, "    %sData%%%s = %sData%%C_obj%%%s\n", nonick, r->name, nonick, r->name);
               }
               else { // characters need to be copied differently
                  if (r->ndims == 0){
                     fprintf(fp, "    %sData%%%s = TRANSFER(%sData%%C_obj%%%s, %sData%%%s )\n", nonick, r->name, nonick, r->name, nonick, r->name);
                  }
               }
            }
        }
      }
    }
  }

  fprintf(fp," END SUBROUTINE %s_C2Fary_Copy%s\n\n", ModName->nickname,nonick ) ;
  return(0) ;
}

int
gen_copy_f2c(FILE         *fp, // *.f90 file we are writting to
    const node_t *ModName, // module name
    char         *inout, // character string written out
    char         *inoutlong) // not sure what this is used for
{
    node_t *q, *r;
    char tmp[NAMELEN];
    char addnick[NAMELEN];
    char nonick[NAMELEN];

    remove_nickname(ModName->nickname, inout, nonick);
    append_nickname((is_a_fast_interface_type(inoutlong)) ? ModName->nickname : "", inoutlong, addnick);
    fprintf(fp, " SUBROUTINE %s_F2C_Copy%s( %sData, ErrStat, ErrMsg, SkipPointers  )\n", ModName->nickname, nonick, nonick);
    fprintf(fp, "    TYPE(%s), INTENT(INOUT) :: %sData\n", addnick, nonick);
    fprintf(fp, "    INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n");
    fprintf(fp, "    CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n");
    fprintf(fp, "    LOGICAL,OPTIONAL,INTENT(IN   ) :: SkipPointers\n");
    fprintf(fp, "    ! \n");
    fprintf(fp, "    LOGICAL                        :: SkipPointers_local\n");
    fprintf(fp, "    ErrStat = ErrID_None\n");
    fprintf(fp, "    ErrMsg  = \"\"\n\n");
    fprintf(fp, "    IF (PRESENT(SkipPointers)) THEN\n");
    fprintf(fp, "       SkipPointers_local = SkipPointers\n");
    fprintf(fp, "    ELSE\n");
    fprintf(fp, "       SkipPointers_local = .false.\n");
    fprintf(fp, "    END IF\n");

    sprintf(tmp, "%s", addnick);

    if ((q = get_entry(make_lower_temp(tmp), ModName->module_ddt_list)) == NULL)
    {
        fprintf(stderr, "Registry warning: generating %s_F2C_Copy%s: cannot find definition for %s\n", ModName->nickname, nonick, tmp);
    }
    else {
        for (r = q->fields; r; r = r->next)
        {
            if (r->type != NULL) {
                if (r->type->type_type == DERIVED) { // && ! r->type->usefrom
                    fprintf(stderr, "Registry WARNING: derived data type %s of type %s is not passed through F-C interface\n", r->name, r->type->name);
                }
                else {
                    if (is_pointer(r)) {
                        fprintf(fp, "\n    ! -- %s %s Data fields\n", r->name, nonick);
                        fprintf(fp, "    IF ( .NOT. SkipPointers_local ) THEN\n");
                        fprintf(fp, "       IF ( .NOT. %s(%sData%%%s)) THEN \n", assoc_or_allocated(r), nonick, r->name);
                        fprintf(fp, "          %sData%%c_obj%%%s_Len = 0\n", nonick, r->name);
                        fprintf(fp, "          %sData%%c_obj%%%s = C_NULL_PTR\n", nonick, r->name);
                        fprintf(fp, "       ELSE\n");
                        fprintf(fp, "          %sData%%c_obj%%%s_Len = SIZE(%sData%%%s)\n", nonick, r->name, nonick, r->name);
                        fprintf(fp, "          IF (%sData%%c_obj%%%s_Len > 0) &\n", nonick, r->name);
                        fprintf(fp, "             %sData%%c_obj%%%s = C_LOC( %sData%%%s( LBOUND(%sData%%%s,1) ) ) \n", nonick, r->name, nonick, r->name, nonick, r->name );
                        fprintf(fp, "       END IF\n");
                        fprintf(fp, "    END IF\n");
                    }
                    else if (!has_deferred_dim(r, 0)) {
                        if (!strcmp(r->type->mapsto, "REAL(ReKi)") ||
                            !strcmp(r->type->mapsto, "REAL(SiKi)") ||
                            !strcmp(r->type->mapsto, "REAL(DbKi)") ||
                            !strcmp(r->type->mapsto, "REAL(R8Ki)") ||
                            !strcmp(r->type->mapsto, "INTEGER(IntKi)") ||
                            !strcmp(r->type->mapsto, "LOGICAL"))
                        {
                            fprintf(fp, "    %sData%%C_obj%%%s = %sData%%%s\n", nonick, r->name, nonick, r->name);
                        }
                        else { // characters need to be copied differently
                            if (r->ndims == 0) {
                                //fprintf(stderr, "Registry WARNING: character data type %s of type %s is not passed through F-C interface\n", r->name, r->type->name);
                                fprintf(fp, "    %sData%%C_obj%%%s = TRANSFER(%sData%%%s, %sData%%C_obj%%%s )\n", nonick, r->name, nonick, r->name, nonick, r->name);
                            }
                        }
                    }
                }
            }
        }
    }

    fprintf(fp, " END SUBROUTINE %s_F2C_Copy%s\n\n", ModName->nickname, nonick);
    return(0);
}


int
gen_copy( FILE * fp, const node_t * ModName, char * inout, char * inoutlong, const node_t * q_in )
{
  char tmp[NAMELEN], tmp2[NAMELEN], addnick[NAMELEN], nonick[NAMELEN] ;
  node_t *q, * r ;
  int d ;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
  fprintf(fp," SUBROUTINE %s_Copy%s( Src%sData, Dst%sData, CtrlCode, ErrStat, ErrMsg )\n",ModName->nickname,nonick,nonick,nonick ) ;
  fprintf(fp, "   TYPE(%s), INTENT(%s) :: Src%sData\n", addnick, (q_in->containsPtr == 1) ? "INOUT" : "IN", nonick);
//fprintf(fp, "   TYPE(%s), INTENT(INOUT) :: Src%sData\n", addnick, nonick);
  fprintf(fp,"   TYPE(%s), INTENT(INOUT) :: Dst%sData\n",addnick,nonick) ;
  fprintf(fp,"   INTEGER(IntKi),  INTENT(IN   ) :: CtrlCode\n") ;
  fprintf(fp,"   INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"   CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"! Local \n") ;
  fprintf(fp,"   INTEGER(IntKi)                 :: i,j,k\n") ;
  for (d = 1; d <= q_in->max_ndims; d++){
  fprintf(fp, "   INTEGER(IntKi)                 :: i%d, i%d_l, i%d_u  !  bounds (upper/lower) for an array dimension %d\n", d, d, d, d);
  }
  fprintf(fp,"   INTEGER(IntKi)                 :: ErrStat2\n") ;
  fprintf(fp,"   CHARACTER(ErrMsgLen)           :: ErrMsg2\n");
  fprintf(fp,"   CHARACTER(*), PARAMETER        :: RoutineName = '%s_Copy%s'\n", ModName->nickname, nonick);
  fprintf(fp, "! \n");
  fprintf(fp,"   ErrStat = ErrID_None\n") ;
  fprintf(fp,"   ErrMsg  = \"\"\n") ;

//  sprintf(tmp,"%s_%s",ModName->nickname,inoutlong) ;
//  sprintf(tmp,"%s",inoutlong) ;
  sprintf(tmp,"%s",addnick) ;

  sprintf(tmp2,"%s",make_lower_temp(tmp)) ;

  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Copy%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
  } else {
    for ( r = q->fields ; r ; r = r->next )
    {
      if ( r->type != NULL ) {

// check if this is an allocatable array:
        if ( r->ndims > 0 && has_deferred_dim(r,0) ) {
             fprintf(fp,"IF (%s(Src%sData%%%s)) THEN\n",assoc_or_allocated(r),nonick,r->name) ;
             strcpy(tmp,"") ;

             for (d = 1; d <= r->ndims; d++) {
                fprintf(fp, "  i%d_l = LBOUND(Src%sData%%%s,%d)\n", d, nonick, r->name, d);
                fprintf(fp, "  i%d_u = UBOUND(Src%sData%%%s,%d)\n", d, nonick, r->name, d);
                sprintf(tmp2, ",i%d_l:i%d_u", d, d);
                strcat(tmp, tmp2);
             }
//fprintf(fp," nonick=%s\n", nonick    );
  fprintf(fp,"  IF (.NOT. %s(Dst%sData%%%s)) THEN \n",assoc_or_allocated(r),nonick,r->name) ;
  fprintf(fp,"    ALLOCATE(Dst%sData%%%s(%s),STAT=ErrStat2)\n",nonick,r->name,(char*)&(tmp[1])) ;
  fprintf(fp,"    IF (ErrStat2 /= 0) THEN \n") ;
  fprintf(fp,"      CALL SetErrStat(ErrID_Fatal, 'Error allocating Dst%sData%%%s.', ErrStat, ErrMsg,RoutineName)\n",nonick,r->name);
  fprintf(fp,"      RETURN\n") ;
  fprintf(fp,"    END IF\n") ;

             if ( sw_ccode && is_pointer(r) ) { // bjj: this needs to be updated if we've got multiple dimension arrays
  fprintf(fp,"    Dst%sData%%c_obj%%%s_Len = SIZE(Dst%sData%%%s)\n",nonick,r->name,nonick,r->name) ;
  fprintf(fp,"    IF (Dst%sData%%c_obj%%%s_Len > 0) &\n",nonick,r->name) ;
  fprintf(fp,"      Dst%sData%%c_obj%%%s = C_LOC( Dst%sData%%%s(i1_l) ) \n",nonick,r->name, nonick,r->name ) ;
             }

  fprintf(fp,"  END IF\n") ; // end dest allocated/associated
        }

        if ( r->type->type_type == DERIVED  ) { // includes mesh and dll_type

            for (d = r->ndims; d >= 1; d--) {
  fprintf(fp,"    DO i%d = LBOUND(Src%sData%%%s,%d), UBOUND(Src%sData%%%s,%d)\n",d,nonick,r->name,d,nonick,r->name,d  ) ;
            }

            if (!strcmp(r->type->name, "meshtype")) {
  fprintf(fp,"      CALL MeshCopy( Src%sData%%%s%s, Dst%sData%%%s%s, CtrlCode, ErrStat2, ErrMsg2 )\n",nonick,r->name,dimstr(r->ndims),nonick,r->name,dimstr(r->ndims)) ;
  fprintf(fp,"         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)\n");
  fprintf(fp,"         IF (ErrStat>=AbortErrLev) RETURN\n");
            } else if ( !strcmp( r->type->name, "dll_type" ) ) {
  fprintf(fp,"      Dst%sData%%%s = Src%sData%%%s\n",nonick,r->name,nonick,r->name) ;
            }
            else { // && ! r->type->usefrom ) {
               char nonick2[NAMELEN];
               remove_nickname(r->type->module->nickname, r->type->name, nonick2);

               fprintf(fp, "      CALL %s_Copy%s( Src%sData%%%s%s, Dst%sData%%%s%s, CtrlCode, ErrStat2, ErrMsg2 )\n",
                  r->type->module->nickname, fast_interface_type_shortname(nonick2),
                  nonick, r->name, dimstr(r->ndims),
                  nonick, r->name, dimstr(r->ndims));
               fprintf(fp, "         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n");
               fprintf(fp, "         IF (ErrStat>=AbortErrLev) RETURN\n");

            }

            for ( d = r->ndims ; d >= 1 ; d-- ) {
  fprintf(fp,"    ENDDO\n") ;
            }
        } else { // not a derived type
  fprintf(fp, "    Dst%sData%%%s = Src%sData%%%s\n",nonick,r->name,nonick,r->name) ;
  if (sw_ccode && !is_pointer(r)){

     //if (!strcmp(r->type->mapsto, "REAL(ReKi)") ||
     //   !strcmp(r->type->mapsto, "REAL(SiKi)") ||
     //   !strcmp(r->type->mapsto, "REAL(DbKi)") ||
     //   !strcmp(r->type->mapsto, "REAL(R8Ki)") ||
     //   !strcmp(r->type->mapsto, "INTEGER(IntKi)") ||
     //   !strcmp(r->type->mapsto, "LOGICAL") ||
     //   r->ndims == 0)
             if ( r->ndims == 0 ) // scalar of any type OR a character array
             {
                //  fprintf(fp, "    Dst%sData%%C_obj%%%s = Dst%sData%%%s\n", nonick, r->name, nonick, r->name);
                   fprintf(fp, "    Dst%sData%%C_obj%%%s = Src%sData%%C_obj%%%s\n", nonick, r->name, nonick, r->name);
             }
          }
        }

// close IF (check on allocatable array)
        if ( r->ndims > 0 && has_deferred_dim(r,0) ) {
  fprintf(fp,"ENDIF\n") ;
        }

      } // if non-null field
    } // each field
  }

  fprintf(fp," END SUBROUTINE %s_Copy%s\n\n", ModName->nickname,nonick ) ;
  return(0) ;
}

void
gen_pack( FILE * fp, const node_t * ModName, char * inout, char *inoutlong )
{

  char tmp[NAMELEN], tmp2[NAMELEN], addnick[NAMELEN], nonick[NAMELEN] ;
  char nonick2[NAMELEN], indent[NAMELEN], mainIndent[6];
  node_t *q, * r ;
  int frst, d, i;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
//  sprintf(tmp,"%s_%s",ModName->nickname,inoutlong) ;
//  sprintf(tmp,"%s",inoutlong) ;
  sprintf(tmp,"%s",addnick) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_Pack%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
    return;//(1) ;
  }

  fprintf(fp, " SUBROUTINE %s_Pack%s( ReKiBuf, DbKiBuf, IntKiBuf, Indata, ErrStat, ErrMsg, SizeOnly )\n", ModName->nickname,nonick) ;
  fprintf(fp, "  REAL(ReKi),       ALLOCATABLE, INTENT(  OUT) :: ReKiBuf(:)\n") ;
  fprintf(fp, "  REAL(DbKi),       ALLOCATABLE, INTENT(  OUT) :: DbKiBuf(:)\n") ;
  fprintf(fp, "  INTEGER(IntKi),   ALLOCATABLE, INTENT(  OUT) :: IntKiBuf(:)\n") ;
  fprintf(fp, "  TYPE(%s),  INTENT(IN) :: InData\n",addnick ) ;
  fprintf(fp, "  INTEGER(IntKi),   INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp, "  CHARACTER(*),     INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp, "  LOGICAL,OPTIONAL, INTENT(IN   ) :: SizeOnly\n") ;
  fprintf(fp, "    ! Local variables\n") ;
  fprintf(fp, "  INTEGER(IntKi)                 :: Re_BufSz\n") ;
  fprintf(fp, "  INTEGER(IntKi)                 :: Re_Xferred\n") ;
  fprintf(fp, "  INTEGER(IntKi)                 :: Db_BufSz\n") ;
  fprintf(fp, "  INTEGER(IntKi)                 :: Db_Xferred\n") ;
  fprintf(fp, "  INTEGER(IntKi)                 :: Int_BufSz\n") ;
  fprintf(fp, "  INTEGER(IntKi)                 :: Int_Xferred\n") ;
  fprintf(fp, "  INTEGER(IntKi)                 :: i,i1,i2,i3,i4,i5\n") ;
  fprintf(fp, "  LOGICAL                        :: OnlySize ! if present and true, do not pack, just allocate buffers\n") ;
  fprintf(fp, "  INTEGER(IntKi)                 :: ErrStat2\n");
  fprintf(fp, "  CHARACTER(ErrMsgLen)           :: ErrMsg2\n");
  fprintf(fp, "  CHARACTER(*), PARAMETER        :: RoutineName = '%s_Pack%s'\n", ModName->nickname, nonick);

  fprintf(fp, " ! buffers to store subtypes, if any\n");
  fprintf(fp, "  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)\n");
  fprintf(fp, "  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)\n");
  fprintf(fp, "  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)\n\n");

  fprintf(fp,"  OnlySize = .FALSE.\n") ;
  fprintf(fp,"  IF ( PRESENT(SizeOnly) ) THEN\n") ;
  fprintf(fp,"    OnlySize = SizeOnly\n") ;
  fprintf(fp,"  ENDIF\n") ;
  fprintf(fp,"    !\n") ;

  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp,"  Re_BufSz  = 0\n") ;
  fprintf(fp,"  Db_BufSz  = 0\n") ;
  fprintf(fp,"  Int_BufSz  = 0\n") ;


  frst = 1;
  for (r = q->fields; r; r = r->next)
  {
     if (r->type == NULL) {
        fprintf(stderr, "Registry warning generating %s_Pack%s: %s has no type.\n", ModName->nickname, nonick, r->name);
        return; // EARLY RETURN
     }

    if (has_deferred_dim(r, 0)){
  //fprintf(fp, "\n");
  fprintf(fp, "  Int_BufSz   = Int_BufSz   + 1     ! %s allocated yes/no\n", r->name);

  fprintf(fp, "  IF ( %s(InData%%%s) ) THEN\n", assoc_or_allocated(r), r->name);
  fprintf(fp, "    Int_BufSz   = Int_BufSz   + 2*%d  ! %s upper/lower bounds for each dimension\n", r->ndims, r->name);
    }

    if (!strcmp(r->type->name, "meshtype") ||
        !strcmp(r->type->name, "dll_type") ||
        (r->type->type_type == DERIVED)   ) {  //  call individual routines to pack data from subtypes:

       if (frst == 1) {
  fprintf(fp, "   ! Allocate buffers for subtypes, if any (we'll get sizes from these) \n"); frst = 0;
       }

       for (d = r->ndims; d >= 1; d--) {
  fprintf(fp, "    DO i%d = LBOUND(InData%%%s,%d), UBOUND(InData%%%s,%d)\n", d, r->name, d, r->name, d);
       }
  fprintf(fp, "      Int_BufSz   = Int_BufSz + 3  ! %s: size of buffers for each call to pack subtype\n", r->name);

      if ( !strcmp( r->type->name, "meshtype" ) ) {
  fprintf(fp, "      CALL MeshPack( InData%%%s%s, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2, .TRUE. ) ! %s \n",
                                 r->name,dimstr(r->ndims),r->name ) ;
      } else if ( !strcmp( r->type->name, "dll_type" ) ) {
  fprintf(fp, "      CALL DLLTypePack( InData%%%s%s, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2, .TRUE. ) ! %s \n",
                                 r->name,dimstr(r->ndims), r->name ) ;
      } else if (r->type->type_type == DERIVED) { // && ! r->type->usefrom ) {
       remove_nickname(r->type->module->nickname, r->type->name, nonick2);
  fprintf(fp, "      CALL %s_Pack%s( Re_Buf, Db_Buf, Int_Buf, InData%%%s%s, ErrStat2, ErrMsg2, .TRUE. ) ! %s \n",
          r->type->module->nickname, fast_interface_type_shortname(nonick2), r->name,
          dimstr(r->ndims), r->name);
      }

  fprintf(fp, "        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)\n");
  fprintf(fp, "        IF (ErrStat >= AbortErrLev) RETURN\n\n");

  fprintf(fp, "      IF(ALLOCATED(Re_Buf)) THEN ! %s\n", r->name);
  fprintf(fp, "         Re_BufSz  = Re_BufSz  + SIZE( Re_Buf  )\n");
  fprintf(fp, "         DEALLOCATE(Re_Buf)\n");
  fprintf(fp, "      END IF\n");

  fprintf(fp, "      IF(ALLOCATED(Db_Buf)) THEN ! %s\n", r->name);
  fprintf(fp, "         Db_BufSz  = Db_BufSz  + SIZE( Db_Buf  )\n");
  fprintf(fp, "         DEALLOCATE(Db_Buf)\n");
  fprintf(fp, "      END IF\n");

  fprintf(fp, "      IF(ALLOCATED(Int_Buf)) THEN ! %s\n", r->name);
  fprintf(fp, "         Int_BufSz = Int_BufSz + SIZE( Int_Buf )\n");
  fprintf(fp, "         DEALLOCATE(Int_Buf)\n");
  fprintf(fp, "      END IF\n");

      for (d = r->ndims; d >= 1; d--) {
  fprintf(fp, "    END DO\n");
      }

    } else {  // intrinsic data types

       // do all dimensions of arrays (no need for loop over i%d)

      sprintf(tmp2, "SIZE(InData%%%s)", r->name);

      if      ( !strcmp( r->type->mapsto, "REAL(ReKi)")  ||
                !strcmp( r->type->mapsto, "REAL(SiKi)")   ) {
  fprintf(fp, "      Re_BufSz   = Re_BufSz   + %s  ! %s\n", (r->ndims>0) ? tmp2 : "1", r->name);
      }
      else if ( !strcmp( r->type->mapsto, "REAL(DbKi)")   ||
                !strcmp(r->type->mapsto, "REAL(R8Ki)")) {
  fprintf(fp, "      Db_BufSz   = Db_BufSz   + %s  ! %s\n", (r->ndims>0) ? tmp2 : "1", r->name);
      }
      else if ( !strcmp( r->type->mapsto, "INTEGER(IntKi)") ||
                !strcmp( r->type->mapsto, "LOGICAL" )       ) {
  fprintf(fp, "      Int_BufSz  = Int_BufSz  + %s  ! %s\n", (r->ndims>0) ? tmp2 : "1", r->name);
      }
      else /*if (!strcmp(r->type->mapsto, "CHARACTER")) */{
  fprintf(fp, "      Int_BufSz  = Int_BufSz  + %s*LEN(InData%%%s)  ! %s\n", (r->ndims>0) ? tmp2 : "1", r->name, r->name);
      }
      /*else
      {
  fprintf(fp,"!  missing buffer for %s\n",r->name ) ;
      }*/
    }

    if (has_deferred_dim(r, 0)){
  fprintf(fp, "  END IF\n");
    }
  // fprintf(fp, "\n"); // space between variables


  }

   // Allocate buffers
  fprintf(fp, "  IF ( Re_BufSz  .GT. 0 ) THEN \n");
  fprintf(fp, "     ALLOCATE( ReKiBuf(  Re_BufSz  ), STAT=ErrStat2 )\n");
  fprintf(fp, "     IF (ErrStat2 /= 0) THEN \n");
  fprintf(fp, "       CALL SetErrStat(ErrID_Fatal, 'Error allocating ReKiBuf.', ErrStat, ErrMsg,RoutineName)\n");
  fprintf(fp, "       RETURN\n");
  fprintf(fp, "     END IF\n");
  fprintf(fp, "  END IF\n");

  fprintf(fp, "  IF ( Db_BufSz  .GT. 0 ) THEN \n");
  fprintf(fp, "     ALLOCATE( DbKiBuf(  Db_BufSz  ), STAT=ErrStat2 )\n");
  fprintf(fp, "     IF (ErrStat2 /= 0) THEN \n");
  fprintf(fp, "       CALL SetErrStat(ErrID_Fatal, 'Error allocating DbKiBuf.', ErrStat, ErrMsg,RoutineName)\n");
  fprintf(fp, "       RETURN\n");
  fprintf(fp, "     END IF\n");
  fprintf(fp, "  END IF\n");

  fprintf(fp, "  IF ( Int_BufSz  .GT. 0 ) THEN \n");
  fprintf(fp, "     ALLOCATE( IntKiBuf(  Int_BufSz  ), STAT=ErrStat2 )\n");
  fprintf(fp, "     IF (ErrStat2 /= 0) THEN \n");
  fprintf(fp, "       CALL SetErrStat(ErrID_Fatal, 'Error allocating IntKiBuf.', ErrStat, ErrMsg,RoutineName)\n");
  fprintf(fp, "       RETURN\n");
  fprintf(fp, "     END IF\n");
  fprintf(fp, "  END IF\n");
  fprintf(fp, "  IF(OnlySize) RETURN ! return early if only trying to allocate buffers (not pack them)\n\n");

  if (sw_ccode) {
     fprintf(fp, "  IF (C_ASSOCIATED(InData%%C_obj%%object)) ");
     fprintf(fp, "CALL SetErrStat(ErrID_Severe,'C_obj%%object cannot be packed.',ErrStat,ErrMsg,RoutineName)\n\n");
  }


  fprintf(fp, "  Re_Xferred  = 1\n");
  fprintf(fp, "  Db_Xferred  = 1\n");
  fprintf(fp, "  Int_Xferred = 1\n\n");


   // Pack data
  for ( r = q->fields ; r ; r = r->next )
  {

     if (has_deferred_dim(r, 0)) {
        // store whether the data type is allocated and the bounds of each dimension
        fprintf(fp, "  IF ( .NOT. %s(InData%%%s) ) THEN\n", assoc_or_allocated(r), r->name);
        fprintf(fp, "    IntKiBuf( Int_Xferred ) = 0\n"); // not allocated
        fprintf(fp, "    Int_Xferred = Int_Xferred + 1\n");
        //fprintf(fp, "    IntKiBuf( Int_Xferred:Int_Xferred+2*%d-1 ) = 0\n", r->ndims, r->name);
        //fprintf(fp, "    Int_Xferred = Int_Xferred + 2*%d\n", r->ndims);
        fprintf(fp, "  ELSE\n");
        fprintf(fp, "    IntKiBuf( Int_Xferred ) = 1\n"); // allocated
        fprintf(fp, "    Int_Xferred = Int_Xferred + 1\n");
        for (d = 1; d <= r->ndims; d++) {
           fprintf(fp, "    IntKiBuf( Int_Xferred    ) = LBOUND(InData%%%s,%d)\n", r->name, d);
           fprintf(fp, "    IntKiBuf( Int_Xferred + 1) = UBOUND(InData%%%s,%d)\n", r->name, d);
           fprintf(fp, "    Int_Xferred = Int_Xferred + 2\n");
        }
        fprintf(fp, "\n");
        strcpy(mainIndent, "  ");
     } 
     else {
        strcpy(mainIndent, "");
    }


    if (!strcmp(r->type->name, "meshtype") ||
       !strcmp(r->type->name, "dll_type") ||
       (r->type->type_type == DERIVED)) {  //  call individual routines to pack data from subtypes:

       if (frst == 1) {
          fprintf(fp, "   ! Allocate buffers for subtypes, if any (we'll get sizes from these) \n"); frst = 0;
       }

       for (d = r->ndims; d >= 1; d--) {
          fprintf(fp, "    DO i%d = LBOUND(InData%%%s,%d), UBOUND(InData%%%s,%d)\n", d, r->name, d, r->name, d);
       }

       if (!strcmp(r->type->name, "meshtype")) {
          fprintf(fp, "      CALL MeshPack( InData%%%s%s, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2, OnlySize ) ! %s \n",
             r->name, dimstr(r->ndims), r->name);
       }
       else if (!strcmp(r->type->name, "dll_type")) {
          fprintf(fp, "      CALL DLLTypePack( InData%%%s%s, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2, OnlySize ) ! %s \n",
             r->name, dimstr(r->ndims), r->name);
       }
       else if (r->type->type_type == DERIVED) { // && ! r->type->usefrom ) {
          remove_nickname(r->type->module->nickname, r->type->name, nonick2);
          fprintf(fp, "      CALL %s_Pack%s( Re_Buf, Db_Buf, Int_Buf, InData%%%s%s, ErrStat2, ErrMsg2, OnlySize ) ! %s \n",
             r->type->module->nickname, fast_interface_type_shortname(nonick2), r->name,
             dimstr(r->ndims),r->name);
       }
       fprintf(fp, "        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)\n");
       fprintf(fp, "        IF (ErrStat >= AbortErrLev) RETURN\n\n");

       fprintf(fp, "      IF(ALLOCATED(Re_Buf)) THEN\n");
       fprintf(fp, "        IntKiBuf( Int_Xferred ) = SIZE(Re_Buf); Int_Xferred = Int_Xferred + 1\n");
       fprintf(fp, "        IF (SIZE(Re_Buf) > 0) ReKiBuf( Re_Xferred:Re_Xferred+SIZE(Re_Buf)-1 ) = Re_Buf\n");
       fprintf(fp, "        Re_Xferred = Re_Xferred + SIZE(Re_Buf)\n");
       fprintf(fp, "        DEALLOCATE(Re_Buf)\n");
       fprintf(fp, "      ELSE\n");
       fprintf(fp, "        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1\n");
       fprintf(fp, "      ENDIF\n");

       fprintf(fp, "      IF(ALLOCATED(Db_Buf)) THEN\n");
       fprintf(fp, "        IntKiBuf( Int_Xferred ) = SIZE(Db_Buf); Int_Xferred = Int_Xferred + 1\n");
       fprintf(fp, "        IF (SIZE(Db_Buf) > 0) DbKiBuf( Db_Xferred:Db_Xferred+SIZE(Db_Buf)-1 ) = Db_Buf\n");
       fprintf(fp, "        Db_Xferred = Db_Xferred + SIZE(Db_Buf)\n");
       fprintf(fp, "        DEALLOCATE(Db_Buf)\n");
       fprintf(fp, "      ELSE\n");
       fprintf(fp, "        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1\n");
       fprintf(fp, "      ENDIF\n");

       fprintf(fp, "      IF(ALLOCATED(Int_Buf)) THEN\n");
       fprintf(fp, "        IntKiBuf( Int_Xferred ) = SIZE(Int_Buf); Int_Xferred = Int_Xferred + 1\n");
       fprintf(fp, "        IF (SIZE(Int_Buf) > 0) IntKiBuf( Int_Xferred:Int_Xferred+SIZE(Int_Buf)-1 ) = Int_Buf\n");
       fprintf(fp, "        Int_Xferred = Int_Xferred + SIZE(Int_Buf)\n");
       fprintf(fp, "        DEALLOCATE(Int_Buf)\n");
       fprintf(fp, "      ELSE\n");
       fprintf(fp, "        IntKiBuf( Int_Xferred ) = 0; Int_Xferred = Int_Xferred + 1\n");
       fprintf(fp, "      ENDIF\n");

       for (d = r->ndims; d >= 1; d--) {
          fprintf(fp, "    END DO\n");
       }

    }
    else {  
       // intrinsic data types
       // do all dimensions of arrays (no need for loop over i%d)

       strcpy(indent, "  ");
       strcat(indent, mainIndent);
       for (d = r->ndims; d >= 1; d--) {
          fprintf(fp, "%s  DO i%d = LBOUND(InData%%%s,%d), UBOUND(InData%%%s,%d)\n", indent, d, r->name, d, r->name, d);
          strcat(indent, "  "); //create an indent
       }


       if (!strcmp(r->type->mapsto, "REAL(ReKi)") ||
          !strcmp(r->type->mapsto, "REAL(SiKi)")) {
          fprintf(fp, "%s  ReKiBuf(Re_Xferred) = InData%%%s%s\n", indent, r->name, dimstr(r->ndims));
          fprintf(fp, "%s  Re_Xferred = Re_Xferred + 1\n", indent);
       }
       else if (!strcmp(r->type->mapsto, "REAL(DbKi)") ||
                !strcmp(r->type->mapsto, "REAL(R8Ki)")) {
          fprintf(fp, "%s  DbKiBuf(Db_Xferred) = InData%%%s%s\n", indent, r->name, dimstr(r->ndims));
          fprintf(fp, "%s  Db_Xferred = Db_Xferred + 1\n", indent);
       }
       else if (!strcmp(r->type->mapsto, "INTEGER(IntKi)") ) {
          fprintf(fp, "%s  IntKiBuf(Int_Xferred) = InData%%%s%s\n", indent, r->name, dimstr(r->ndims));
          fprintf(fp, "%s  Int_Xferred = Int_Xferred + 1\n", indent);
       }
       else if (!strcmp(r->type->mapsto, "LOGICAL") ) {
          fprintf(fp, "%s  IntKiBuf(Int_Xferred) = TRANSFER(InData%%%s%s, IntKiBuf(1))\n", indent, r->name, dimstr(r->ndims));
          fprintf(fp, "%s  Int_Xferred = Int_Xferred + 1\n", indent);
       }

       else /*if (!strcmp(r->type->mapsto, "CHARACTER")) */{

          fprintf(fp, "%s  DO I = 1, LEN(InData%%%s)\n", indent, r->name);
          fprintf(fp, "%s    IntKiBuf(Int_Xferred) = ICHAR(InData%%%s%s(I:I), IntKi)\n", indent, r->name, dimstr(r->ndims));
          fprintf(fp, "%s    Int_Xferred = Int_Xferred + 1\n", indent);
          fprintf(fp, "%s  END DO ! I\n", indent);

       }

       for (d = r->ndims; d >= 1; d--) {
          strcpy(indent, "  ");
          strcat(indent, mainIndent);
          for (i = 1; i < d; i++) {
             strcat(indent, "  ");
          }
          fprintf(fp, "%s  END DO\n", indent);
       }

    }

    if (has_deferred_dim(r, 0)){
       fprintf(fp, "  END IF\n");
    }
  }

  fprintf(fp," END SUBROUTINE %s_Pack%s\n\n", ModName->nickname,nonick ) ;
  return;//(0) ;
}

void
gen_unpack( FILE * fp, const node_t * ModName, char * inout, char * inoutlong )
{
  char tmp[NAMELEN], tmp2[NAMELEN], indent[NAMELEN], addnick[NAMELEN], nonick[NAMELEN], nonick2[NAMELEN], mainIndent[6];
  node_t *q, * r ;
  int d, i ;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
//  sprintf(tmp,"%s_%s",ModName->nickname,inoutlong) ;
//  sprintf(tmp,"%s",inoutlong) ;
  sprintf(tmp,"%s",addnick) ;
  if (( q = get_entry( make_lower_temp(tmp),ModName->module_ddt_list ) ) == NULL )
  {
    fprintf(stderr,"Registry warning: generating %s_UnPack%s: cannot find definition for %s\n",ModName->nickname,nonick,tmp) ;
    return;//(1) ;
  }

  fprintf(fp," SUBROUTINE %s_UnPack%s( ReKiBuf, DbKiBuf, IntKiBuf, Outdata, ErrStat, ErrMsg )\n", ModName->nickname,nonick ) ;
  fprintf(fp,"  REAL(ReKi),      ALLOCATABLE, INTENT(IN   ) :: ReKiBuf(:)\n") ;
  fprintf(fp,"  REAL(DbKi),      ALLOCATABLE, INTENT(IN   ) :: DbKiBuf(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),  ALLOCATABLE, INTENT(IN   ) :: IntKiBuf(:)\n") ;
  fprintf(fp,"  TYPE(%s), INTENT(INOUT) :: OutData\n",addnick ) ;
  fprintf(fp,"  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"    ! Local variables\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Buf_size\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Re_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Db_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: Int_Xferred\n") ;
  fprintf(fp,"  INTEGER(IntKi)                 :: i\n") ;
  for (d = 1; d <= q->max_ndims; d++){
  fprintf(fp,"  INTEGER(IntKi)                 :: i%d, i%d_l, i%d_u  !  bounds (upper/lower) for an array dimension %d\n", d, d, d, d);
  }
  fprintf(fp, "  INTEGER(IntKi)                 :: ErrStat2\n");
  fprintf(fp, "  CHARACTER(ErrMsgLen)           :: ErrMsg2\n");
  fprintf(fp, "  CHARACTER(*), PARAMETER        :: RoutineName = '%s_UnPack%s'\n", ModName->nickname, nonick);

  fprintf(fp," ! buffers to store meshes, if any\n") ;
  fprintf(fp,"  REAL(ReKi),      ALLOCATABLE   :: Re_Buf(:)\n") ;
  fprintf(fp,"  REAL(DbKi),      ALLOCATABLE   :: Db_Buf(:)\n") ;
  fprintf(fp,"  INTEGER(IntKi),  ALLOCATABLE   :: Int_Buf(:)\n") ;
  fprintf(fp,"    !\n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp,"  Re_Xferred  = 1\n") ;
  fprintf(fp,"  Db_Xferred  = 1\n") ;
  fprintf(fp,"  Int_Xferred  = 1\n") ;


// BJJ: TODO:  if there are C types, we're going to have to associate with C data structures....

  // Unpack data
  for (r = q->fields; r; r = r->next)
  {

     strcpy(tmp, "");
     if (has_deferred_dim(r, 0)){
        // determine if the array was allocated when packed:
        fprintf(fp, "  IF ( IntKiBuf( Int_Xferred ) == 0 ) THEN  ! %s not allocated\n", r->name); // not allocated
        fprintf(fp, "    Int_Xferred = Int_Xferred + 1\n");
      //fprintf(fp, "    Int_Xferred = Int_Xferred + 2*%d\n", r->ndims);
        fprintf(fp, "  ELSE\n");
        fprintf(fp, "    Int_Xferred = Int_Xferred + 1\n");

        for (d = 1; d <= r->ndims; d++) {
           fprintf(fp, "    i%d_l = IntKiBuf( Int_Xferred    )\n", d); //fprintf(fp, "    IntKiBuf( Int_Xferred    ) = LBOUND(OutData%%%s,%d)\n", r->name, d);
           fprintf(fp, "    i%d_u = IntKiBuf( Int_Xferred + 1)\n", d); //fprintf(fp, "    IntKiBuf( Int_Xferred + 1) = UBOUND(OutData%%%s,%d)\n", r->name, d);
           fprintf(fp, "    Int_Xferred = Int_Xferred + 2\n");
           sprintf(tmp2, ",i%d_l:i%d_u", d, d);
           strcat(tmp, tmp2);
        }

        fprintf(fp, "    IF (%s(OutData%%%s)) DEALLOCATE(OutData%%%s)\n", assoc_or_allocated(r), r->name, r->name); // BJJ: need NULLIFY(), too?
        fprintf(fp, "    ALLOCATE(OutData%%%s(%s),STAT=ErrStat2)\n", r->name, (char*)&(tmp[1]));
        fprintf(fp, "    IF (ErrStat2 /= 0) THEN \n");
        fprintf(fp, "       CALL SetErrStat(ErrID_Fatal, 'Error allocating OutData%%%s.', ErrStat, ErrMsg,RoutineName)\n", r->name);
        fprintf(fp, "       RETURN\n");
        fprintf(fp, "    END IF\n");

        if (sw_ccode && is_pointer(r)) { // bjj: this needs to be updated if we've got multiple dimension arrays
          fprintf(fp, "    OutData%%c_obj%%%s_Len = SIZE(OutData%%%s)\n", r->name, r->name);
          fprintf(fp, "    IF (OutData%%c_obj%%%s_Len > 0) &\n", r->name);
          fprintf(fp, "       OutData%%c_obj%%%s = C_LOC( OutData%%%s(i1_l) ) \n", r->name, r->name);
        }
        strcpy(mainIndent, "  ");
     }
     else{
        for (d = 1; d <= r->ndims; d++) {
           fprintf(fp, "    i%d_l = LBOUND(OutData%%%s,%d)\n", d, r->name, d);
           fprintf(fp, "    i%d_u = UBOUND(OutData%%%s,%d)\n", d, r->name, d);
           sprintf(tmp2, ",i%d_l:i%d_u", d, d);
           strcat(tmp, tmp2);
        }
        strcpy(mainIndent, "");
     }

     if (!strcmp(r->type->name, "meshtype") ||
        !strcmp(r->type->name, "dll_type") ||
        (r->type->type_type == DERIVED)) {  //  call individual routines to pack data from subtypes:

        for (d = r->ndims; d >= 1; d--) {
           fprintf(fp, "    DO i%d = LBOUND(OutData%%%s,%d), UBOUND(OutData%%%s,%d)\n", d, r->name, d, r->name, d);
        }

        // initialize buffers to send to subtype-unpack routines:
        // reals:
        fprintf(fp, "      Buf_size=IntKiBuf( Int_Xferred )\n");
        fprintf(fp, "      Int_Xferred = Int_Xferred + 1\n");
        fprintf(fp, "      IF(Buf_size > 0) THEN\n");
        fprintf(fp, "        ALLOCATE(Re_Buf(Buf_size),STAT=ErrStat2)\n");
        fprintf(fp, "        IF (ErrStat2 /= 0) THEN \n");
        fprintf(fp, "           CALL SetErrStat(ErrID_Fatal, 'Error allocating Re_Buf.', ErrStat, ErrMsg,RoutineName)\n");
        fprintf(fp, "           RETURN\n");
        fprintf(fp, "        END IF\n");

        fprintf(fp, "        Re_Buf = ReKiBuf( Re_Xferred:Re_Xferred+Buf_size-1 )\n");
        fprintf(fp, "        Re_Xferred = Re_Xferred + Buf_size\n");
        fprintf(fp, "      END IF\n");

         // doubles:
        fprintf(fp, "      Buf_size=IntKiBuf( Int_Xferred )\n");
        fprintf(fp, "      Int_Xferred = Int_Xferred + 1\n");
        fprintf(fp, "      IF(Buf_size > 0) THEN\n");
        fprintf(fp, "        ALLOCATE(Db_Buf(Buf_size),STAT=ErrStat2)\n");
        fprintf(fp, "        IF (ErrStat2 /= 0) THEN \n");
        fprintf(fp, "           CALL SetErrStat(ErrID_Fatal, 'Error allocating Db_Buf.', ErrStat, ErrMsg,RoutineName)\n");
        fprintf(fp, "           RETURN\n");
        fprintf(fp, "        END IF\n");

        fprintf(fp, "        Db_Buf = DbKiBuf( Db_Xferred:Db_Xferred+Buf_size-1 )\n");
        fprintf(fp, "        Db_Xferred = Db_Xferred + Buf_size\n");
        fprintf(fp, "      END IF\n");

        // integers:
        fprintf(fp, "      Buf_size=IntKiBuf( Int_Xferred )\n");
        fprintf(fp, "      Int_Xferred = Int_Xferred + 1\n");
        fprintf(fp, "      IF(Buf_size > 0) THEN\n");
        fprintf(fp, "        ALLOCATE(Int_Buf(Buf_size),STAT=ErrStat2)\n");
        fprintf(fp, "        IF (ErrStat2 /= 0) THEN \n");
        fprintf(fp, "           CALL SetErrStat(ErrID_Fatal, 'Error allocating Int_Buf.', ErrStat, ErrMsg,RoutineName)\n");
        fprintf(fp, "           RETURN\n");
        fprintf(fp, "        END IF\n");

        fprintf(fp, "        Int_Buf = IntKiBuf( Int_Xferred:Int_Xferred+Buf_size-1 )\n");
        fprintf(fp, "        Int_Xferred = Int_Xferred + Buf_size\n");
        fprintf(fp, "      END IF\n");


        if (!strcmp(r->type->name, "meshtype")) {
           fprintf(fp, "      CALL MeshUnpack( OutData%%%s%s, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2 ) ! %s \n",
              r->name, dimstr(r->ndims), r->name);
        }
        else if (!strcmp(r->type->name, "dll_type")) {
           fprintf(fp, "      CALL DLLTypeUnpack( OutData%%%s%s, Re_Buf, Db_Buf, Int_Buf, ErrStat2, ErrMsg2 ) ! %s \n",
              r->name, dimstr(r->ndims), r->name);
        }
        else if (r->type->type_type == DERIVED) { // && ! r->type->usefrom ) {
           remove_nickname(r->type->module->nickname, r->type->name, nonick2);
           fprintf(fp, "      CALL %s_Unpack%s( Re_Buf, Db_Buf, Int_Buf, OutData%%%s%s, ErrStat2, ErrMsg2 ) ! %s \n",
              r->type->module->nickname, fast_interface_type_shortname(nonick2), r->name,
              dimstr(r->ndims), r->name);
        }
        fprintf(fp, "        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)\n");
        fprintf(fp, "        IF (ErrStat >= AbortErrLev) RETURN\n\n");

        fprintf(fp, "      IF(ALLOCATED(Re_Buf )) DEALLOCATE(Re_Buf )\n");
        fprintf(fp, "      IF(ALLOCATED(Db_Buf )) DEALLOCATE(Db_Buf )\n");
        fprintf(fp, "      IF(ALLOCATED(Int_Buf)) DEALLOCATE(Int_Buf)\n");

        for (d = r->ndims; d >= 1; d--) {
           fprintf(fp, "    END DO\n");
        }

     }
     else
     {
        strcpy(indent, "  ");
        strcat(indent, mainIndent);
        for (d = r->ndims; d >= 1; d--) {
           fprintf(fp, "%s  DO i%d = LBOUND(OutData%%%s,%d), UBOUND(OutData%%%s,%d)\n", indent, d, r->name, d, r->name, d);
           strcat(indent, "  "); //create an indent
        }


        if (!strcmp(r->type->mapsto, "REAL(ReKi)") ||
            !strcmp(r->type->mapsto, "REAL(SiKi)")) {
           if (sw_ccode && is_pointer(r)) {
              fprintf(fp, "%s  OutData%%%s%s = REAL(ReKiBuf(Re_Xferred), C_FLOAT)\n", indent, r->name, dimstr(r->ndims));
           }
           else if (!strcmp(r->type->mapsto, "REAL(SiKi)")) {
              fprintf(fp, "%s  OutData%%%s%s = REAL(ReKiBuf(Re_Xferred), SiKi)\n", indent, r->name, dimstr(r->ndims));
           }
           else {
              fprintf(fp, "%s  OutData%%%s%s = ReKiBuf(Re_Xferred)\n", indent, r->name, dimstr(r->ndims));
           }
           fprintf(fp, "%s  Re_Xferred = Re_Xferred + 1\n", indent);
        }
        else if (!strcmp(r->type->mapsto, "REAL(DbKi)") ||
                 !strcmp(r->type->mapsto, "REAL(R8Ki)")) {
           if (sw_ccode && is_pointer(r)) {
              fprintf(fp, "%s  OutData%%%s%s = REAL(DbKiBuf(Db_Xferred), C_DOUBLE)\n", indent, r->name, dimstr(r->ndims));
           }
           else if (!strcmp(r->type->mapsto, "REAL(R8Ki)")) {
              fprintf(fp, "%s  OutData%%%s%s = REAL(DbKiBuf(Db_Xferred), R8Ki)\n", indent, r->name, dimstr(r->ndims));
           }
           else {
              fprintf(fp, "%s  OutData%%%s%s = DbKiBuf(Db_Xferred)\n", indent, r->name, dimstr(r->ndims));
           }
           fprintf(fp, "%s  Db_Xferred = Db_Xferred + 1\n", indent);
        }
        else if (!strcmp(r->type->mapsto, "INTEGER(IntKi)")) {
           fprintf(fp, "%s  OutData%%%s%s = IntKiBuf(Int_Xferred)\n", indent, r->name, dimstr(r->ndims));
           fprintf(fp, "%s  Int_Xferred = Int_Xferred + 1\n", indent);
        }
        else if (!strcmp(r->type->mapsto, "LOGICAL")) {
           fprintf(fp, "%s  OutData%%%s%s = TRANSFER(IntKiBuf(Int_Xferred), OutData%%%s%s)\n", indent, r->name, dimstr(r->ndims), r->name, dimstr(r->ndims));
           fprintf(fp, "%s  Int_Xferred = Int_Xferred + 1\n", indent);
        }

        else /*if (!strcmp(r->type->mapsto, "CHARACTER")) */ {

           fprintf(fp, "%s  DO I = 1, LEN(OutData%%%s)\n", indent, r->name);
           fprintf(fp, "%s    OutData%%%s%s(I:I) = CHAR(IntKiBuf(Int_Xferred))\n", indent, r->name, dimstr(r->ndims));
           fprintf(fp, "%s    Int_Xferred = Int_Xferred + 1\n", indent);
           fprintf(fp, "%s  END DO ! I\n", indent);

        }

        for (d = r->ndims; d >= 1; d--) {
           strcpy(indent, "  ");
           strcat(indent, mainIndent);
           for (i = 1; i < d; i++) {
              strcat(indent, "  ");
           }
           fprintf(fp, "%s  END DO\n", indent);
        }

// need to move scalars and strings to the %c_obj% type, too!
// compare with copy routine

        if (sw_ccode && !is_pointer(r) && r->ndims == 0) {
              if (!strcmp(r->type->mapsto, "REAL(ReKi)") ||
                 !strcmp(r->type->mapsto, "REAL(SiKi)") ||
                 !strcmp(r->type->mapsto, "REAL(DbKi)") ||
                 !strcmp(r->type->mapsto, "REAL(R8Ki)") ||
                 !strcmp(r->type->mapsto, "INTEGER(IntKi)") ||
                 !strcmp(r->type->mapsto, "LOGICAL"))
              {
                 fprintf(fp, "      OutData%%C_obj%%%s = OutData%%%s\n", r->name, r->name);
              }
              else { // characters need to be copied differently
                 fprintf(fp, "      OutData%%C_obj%%%s = TRANSFER(OutData%%%s, OutData%%C_obj%%%s )\n", r->name, r->name, r->name);
              }
        }

     }

     if (has_deferred_dim(r, 0)){
        fprintf(fp, "  END IF\n");
     }
  }

  fprintf(fp," END SUBROUTINE %s_UnPack%s\n\n", ModName->nickname,nonick ) ;
  return;//(0) ;
}

void
gen_mask_alloc( FILE *fp, int ndims, char *tmp )
{
  if        ( ndims == 1 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1)))\n  mask%d = .TRUE.\n",ndims,tmp,ndims) ;
  } else if ( ndims == 2 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1),SIZE(%s,2)))\n  mask%d = .TRUE.\n",ndims,tmp,tmp,ndims) ;
  } else if ( ndims == 3 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1),SIZE(%s,2),SIZE(%s,3)))\n  mask%d = .TRUE.\n",ndims,tmp,tmp,tmp,ndims) ;
  } else if ( ndims == 4 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1),SIZE(%s,2),SIZE(%s,3),SIZE(%s,4)))\n  mask%d = .TRUE.\n",ndims,tmp,tmp,tmp,tmp,ndims) ;
  } else if ( ndims == 5 ) {
    fprintf(fp,"  ALLOCATE(mask%d(SIZE(%s,1),SIZE(%s,2),SIZE(%s,3),SIZE(%s,4),SIZE(%s,5)))\n  mask%d = .TRUE.\n",ndims,tmp,tmp,tmp,tmp,tmp,ndims) ;
  }
}



int
gen_destroy( FILE * fp, const node_t * ModName, char * inout, char * inoutlong )
{
  char tmp[NAMELEN], addnick[NAMELEN], nonick[NAMELEN] ;
  node_t *q, * r ;
  int d ;

  remove_nickname(ModName->nickname,inout,nonick) ;
  append_nickname((is_a_fast_interface_type(inoutlong))?ModName->nickname:"",inoutlong,addnick) ;
  fprintf(fp, " SUBROUTINE %s_Destroy%s( %sData, ErrStat, ErrMsg )\n",ModName->nickname,nonick,nonick );
  fprintf(fp, "  TYPE(%s), INTENT(INOUT) :: %sData\n",addnick,nonick) ;
  fprintf(fp, "  INTEGER(IntKi),  INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp, "  CHARACTER(*),    INTENT(  OUT) :: ErrMsg\n");
  fprintf(fp, "  CHARACTER(*),    PARAMETER :: RoutineName = '%s_Destroy%s'\n", ModName->nickname, nonick);
  fprintf(fp, "  INTEGER(IntKi)                 :: i, i1, i2, i3, i4, i5 \n");
  fprintf(fp,"! \n") ;
  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp, "  ErrMsg  = \"\"\n");

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
        fprintf(stderr,"Registry warning generating %s_Destroy%s: %s has no type.\n",ModName->nickname,nonick,r->name) ;
      } else {

     if ( r->ndims > 0 && has_deferred_dim(r,0) ) {
      fprintf(fp,"IF (%s(%sData%%%s)) THEN\n",assoc_or_allocated(r),nonick,r->name) ;
     }

     if (r->type->type_type == DERIVED){

        for (d = r->ndims; d >= 1; d--) {
           fprintf(fp, "DO i%d = LBOUND(%sData%%%s,%d), UBOUND(%sData%%%s,%d)\n", d, nonick, r->name, d, nonick, r->name, d);
        }

        if (!strcmp(r->type->name, "meshtype")) {
           fprintf(fp, "  CALL MeshDestroy( %sData%%%s%s, ErrStat, ErrMsg )\n", nonick, r->name, dimstr(r->ndims));
        }
        else if (!strcmp(r->type->name, "dll_type")) {
           fprintf(fp, "  CALL FreeDynamicLib( %sData%%%s%s, ErrStat, ErrMsg )\n", nonick, r->name, dimstr(r->ndims));
        }
        else { //if (r->type->type_type == DERIVED) { // && ! r->type->usefrom ) {
           char nonick2[NAMELEN];
           remove_nickname(r->type->module->nickname, r->type->name, nonick2);
           fprintf(fp, "  CALL %s_Destroy%s( %sData%%%s%s, ErrStat, ErrMsg )\n",
              r->type->module->nickname, fast_interface_type_shortname(nonick2), nonick, r->name, dimstr(r->ndims));
        }

        for (d = r->ndims; d >= 1; d--) {
           fprintf(fp, "ENDDO\n");
        }
     }
     if ( r->ndims > 0 && has_deferred_dim(r,0) ) {
         fprintf(fp,"  DEALLOCATE(%sData%%%s)\n",nonick,r->name) ;
         if ( is_pointer(r) ) {
            fprintf(fp, "  %sData%%%s => NULL()\n",nonick,r->name) ;
            if (sw_ccode){
               fprintf(fp, "  %sData%%C_obj%%%s = C_NULL_PTR\n", nonick, r->name);
               fprintf(fp, "  %sData%%C_obj%%%s_Len = 0\n", nonick, r->name);
            }
         }
         fprintf(fp,"ENDIF\n") ;
      }


      }
    }
  }

  fprintf(fp," END SUBROUTINE %s_Destroy%s\n\n", ModName->nickname,nonick ) ;
  return(0) ;
}


#define MAXRECURSE 9
// HERE
#if 0
void gen_extint_order(FILE *fp, const node_t *ModName, char * typnm, char * uy, const int order, node_t *r, char * deref, int recurselevel) {
   node_t *q, *r1 ;
   int j ;
   int mesh = 0 ;
   char derefrecurse[NAMELEN],tmp[NAMELEN] ;
   if ( recurselevel > MAXRECURSE ) {
     fprintf(stderr,"REGISTRY ERROR: too many levels of array subtypes\n") ;
     exit(9) ;
   }
   if ( r->type != NULL ) {

// check if this is an allocatable array:
     if ( r->ndims > 0 && has_deferred_dim(r,0) ) {
  fprintf(fp,"IF (%s(%s_out%s%%%s) .AND. %s(%s(1)%s%%%s)) THEN\n",assoc_or_allocated(r),uy,deref,r->name,
     assoc_or_allocated(r), uy, deref, r->name);
     }
     if ( r->type->type_type == DERIVED ) {
       if (( q = get_entry( make_lower_temp(r->type->name),ModName->module_ddt_list ) ) != NULL ) {
         for ( r1 = q->fields ; r1 ; r1 = r1->next )
         {
           sprintf(derefrecurse,"%s%%%s",deref,r->name) ;
           for ( j = r->ndims ; j > 0 ; j-- ) {

              fprintf(fp, "  DO i%d%d = LBOUND(%s_out%s,%d),UBOUND(%s_out%s,%d)\n", recurselevel, j, uy, derefrecurse, j, uy, derefrecurse, j);
             sprintf(derefrecurse,"%s%%%s(i%d%d)",deref,r->name,recurselevel,j) ;
           }
           gen_extint_order( fp, ModName, typnm, uy, order, r1, derefrecurse, recurselevel+1 ) ;
           for ( j = r->ndims ; j > 0 ; j-- ) {
  fprintf(fp,"  ENDDO\n") ;
           }
         }
       } else if ( !strcmp( r->type->mapsto, "MeshType" ) ) {
         for ( j = r->ndims ; j > 0 ; j-- ) {
            fprintf(fp, "  DO i%d%d = LBOUND(%s_out%s%%%s,%d),UBOUND(%s_out%s%%%s,%d)\n", 0, j, uy, deref, r->name, j, uy, deref, r->name, j);
         }

         if        ( order == 0 ) {
  fprintf(fp, "  CALL MeshCopy(%s(1)%s%%%s%s, %s_out%s%%%s%s, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )\n", uy, deref, r->name, dimstr(r->ndims)
               , uy, deref, r->name, dimstr(r->ndims));
         } else if ( order == 1 ) {
  fprintf(fp,"  CALL MeshExtrapInterp1(%s(1)%s%%%s%s, %s(2)%s%%%s%s, tin, %s_out%s%%%s%s, tin_out, ErrStat2, ErrMsg2 )\n"
     , uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
         } else if ( order == 2 ) {
  fprintf(fp,"  CALL MeshExtrapInterp2(%s(1)%s%%%s%s, %s(2)%s%%%s%s, %s(3)%s%%%s%s, tin, %s_out%s%%%s%s, tin_out, ErrStat2, ErrMsg2 )\n"
     , uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
         }
   fprintf(fp, "         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n");
   fprintf(fp, "         IF (ErrStat>=AbortErrLev) RETURN\n");

         for ( j = r->ndims ; j > 0 ; j-- ) {
  fprintf(fp,"  ENDDO\n") ;
         }
       } else {


          char nonick2[NAMELEN] ;
          remove_nickname(r->type->module->nickname,r->type->name,nonick2) ;
          strcpy(dimstr(r->ndims),"") ;
          for ( j = r->ndims ; j >= 1 ; j-- ) {
             fprintf(fp, "   DO i%d%d = LBOUND(%s_out%s%%%s,%d), UBOUND(%s_out%s%%%s,%d)\n", 0, j, uy, deref, r->name, j, uy, deref, r->name, j);
             if ( j == r->ndims ) strcat(dimstr(r->ndims),"(") ;
             sprintf(tmp,"i%d%d",0,j) ;
             if ( j == 1 ) strcat(tmp,")") ; else strcat(tmp,",") ;
             strcat(dimstr(r->ndims),tmp) ;
          }


  fprintf(fp,"      CALL %s_%s_ExtrapInterp( %s%s%%%s%s, tin, %s_out%s%%%s%s, tin_out, ErrStat2, ErrMsg2 )\n",
                                r->type->module->nickname,fast_interface_type_shortname(nonick2)
                                , uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
  fprintf(fp,"         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n");
  fprintf(fp,"         IF (ErrStat>=AbortErrLev) RETURN\n");


          for ( j = r->ndims ; j >= 1 ; j-- ) {
  fprintf(fp,"   ENDDO\n") ;
          }

       }
     } else if ( !strcmp( r->type->mapsto, "REAL(ReKi)") ||
                 !strcmp( r->type->mapsto, "REAL(SiKi)") ||
                 !strcmp( r->type->mapsto, "REAL(DbKi)")   ) {
       if        ( r->ndims==0 ) {
       } else if ( r->ndims==1 && order > 0 ) {
          fprintf(fp, "  ALLOCATE(b1(SIZE(%s_out%s%%%s,1)))\n", uy, deref, r->name);
          fprintf(fp, "  ALLOCATE(c1(SIZE(%s_out%s%%%s,1)))\n", uy, deref, r->name);
       } else if ( r->ndims==2 && order > 0 ) {
          fprintf(fp, "  ALLOCATE(b2(SIZE(%s_out%s%%%s,1),SIZE(%s_out%s%%%s,2) ))\n", uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "  ALLOCATE(c2(SIZE(%s_out%s%%%s,1),SIZE(%s_out%s%%%s,2) ))\n", uy, deref, r->name, uy, deref, r->name);
       } else if ( r->ndims==3 && order > 0 ) {
          fprintf(fp, "  ALLOCATE(b3(SIZE(%s_out%s%%%s,1),SIZE(%s_out%s%%%s,2), &\n", uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "              SIZE(%s_out%s%%%s,3)                     ))\n", uy, deref, r->name);
          fprintf(fp, "  ALLOCATE(c3(SIZE(%s_out%s%%%s,1),SIZE(%s_out%s%%%s,2), &\n", uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "              SIZE(%s_out%s%%%s,3)                     ))\n", uy, deref, r->name);
       } else if ( r->ndims==4 && order > 0 ) {
          fprintf(fp, "  ALLOCATE(b4(SIZE(%s_out%s%%%s,1),SIZE(%s_out%s%%%s,2), &\n", uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "              SIZE(%s_out%s%%%s,3),SIZE(%s_out%s%%%s,4) ))\n", uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "  ALLOCATE(c4(SIZE(%s_out%s%%%s,1),SIZE(%s_out%s%%%s,2), &\n", uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "              SIZE(%s_out%s%%%s,3),SIZE(%s_out%s%%%s,4) ))\n", uy, deref, r->name, uy, deref, r->name);
       } else if ( r->ndims==5 && order > 0 ) {
          fprintf(fp, "  ALLOCATE(b5(SIZE(%s_out%s%%%s,1),SIZE(%s_out%s%%%s,2), &\n", uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "              SIZE(%s_out%s%%%s,3),SIZE(%s_out%s%%%s,4), &\n", uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "              SIZE(%s_out%s%%%s,5)                         ))\n", uy, deref, r->name);
          fprintf(fp, "  ALLOCATE(c5(SIZE(%s_out%s%%%s,1),SIZE(%s_out%s%%%s,2), &\n", uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "              SIZE(%s_out%s%%%s,3),SIZE(%s_out%s%%%s,4), &\n", uy, deref, r->name, uy, deref, r->name);
  fprintf(fp, "              SIZE(%s_out%s%%%s,5)                     ))\n", uy, deref, r->name);
       } else                    {
          if (order > 0) fprintf(stderr, "Registry WARNING: too many dimensions for %s%%%s\n", deref, r->name);
       }

       if        ( order == 0 ) {
          fprintf(fp, "  %s_out%s%%%s = %s(1)%s%%%s\n", uy, deref, r->name, uy, deref, r->name);
       } else if ( order == 1 ) {
          fprintf(fp, "  b%d = -(%s(1)%s%%%s - %s(2)%s%%%s)/t(2)\n", r->ndims, uy, deref, r->name, uy, deref, r->name);
          fprintf(fp, "  %s_out%s%%%s = %s(1)%s%%%s + b%d * t_out\n", uy, deref, r->name, uy, deref, r->name, r->ndims);
       } else if ( order == 2 ) {
  fprintf(fp,"  b%d = (t(3)**2*(%s(1)%s%%%s - %s(2)%s%%%s) + t(2)**2*(-%s(1)%s%%%s + %s(3)%s%%%s))/(t(2)*t(3)*(t(2) - t(3)))\n",
     r->ndims, uy, deref, r->name, uy, deref, r->name, uy, deref, r->name, uy, deref, r->name);
  fprintf(fp,"  c%d = ( (t(2)-t(3))*%s(1)%s%%%s + t(3)*%s(2)%s%%%s - t(2)*%s(3)%s%%%s ) / (t(2)*t(3)*(t(2) - t(3)))\n",
     r->ndims, uy, deref, r->name, uy, deref, r->name, uy, deref, r->name);
  fprintf(fp,"  %s_out%s%%%s = %s(1)%s%%%s + b%d * t_out + c%d * t_out**2\n"
     , uy, deref, r->name, uy, deref, r->name, r->ndims, r->ndims);
       }
       if        ( r->ndims>=1 && order > 0 ) {
  fprintf(fp,"  DEALLOCATE(b%d)\n",r->ndims) ;
  fprintf(fp,"  DEALLOCATE(c%d)\n",r->ndims) ;
       }
     }
// check if this is an allocatable array:
     if ( r->ndims > 0 && has_deferred_dim(r,0) ) {
  fprintf(fp,"END IF ! check if allocated\n") ;
     }

   }
}
#endif
void gen_extint_order(FILE *fp, const node_t *ModName, char * typnm, char * uy, const int order, node_t *r, char * deref, int recurselevel) {
   node_t *q, *r1;
   int i, j;
   int mesh = 0;
   char derefrecurse[NAMELEN], indent[NAMELEN], tmp[NAMELEN];
   if (recurselevel > MAXRECURSE) {
      fprintf(stderr, "REGISTRY ERROR: too many levels of array subtypes\n");
      exit(9);
   }
   if (r->type != NULL) {

      // check if this is an allocatable array:
      if (r->ndims > 0 && has_deferred_dim(r, 0)) {
         fprintf(fp, "IF (%s(%s_out%s%%%s) .AND. %s(%s1%s%%%s)) THEN\n", assoc_or_allocated(r), uy, deref, r->name,
            assoc_or_allocated(r), uy, deref, r->name);
      }
      if (r->type->type_type == DERIVED) {

         if ((q = get_entry(make_lower_temp(r->type->name), ModName->module_ddt_list)) != NULL) {
            for (r1 = q->fields; r1; r1 = r1->next)
            {
               sprintf(derefrecurse, "%s%%%s", deref, r->name);

               for (j = r->ndims; j > 0; j--) {
                  fprintf(fp, "  DO i%d%d = LBOUND(%s_out%s,%d),UBOUND(%s_out%s,%d)\n", recurselevel, j, uy, derefrecurse, j, uy, derefrecurse, j);
               }


               if (r->ndims > 0) {
                   strcat(derefrecurse, "(");
                   for (j = 1; j <= r->ndims; j++) {
                       sprintf(tmp, "i%d%d", recurselevel, j);
                       strcat(derefrecurse, tmp);
                       if (j < r->ndims) {
                           strcat(derefrecurse, ",");
                       }
                   }
                   strcat(derefrecurse, ")");
               }

               gen_extint_order(fp, ModName, typnm, uy, order, r1, derefrecurse, recurselevel + 1);
               for (j = r->ndims; j > 0; j--) {
                  fprintf(fp, "  ENDDO\n");
               }
            }
         }

         else {

            for (j = r->ndims; j > 0; j--) {
               fprintf(fp, "  DO i%d = LBOUND(%s_out%s%%%s,%d),UBOUND(%s_out%s%%%s,%d)\n", j, uy, deref, r->name, j, uy, deref, r->name, j);
            }

            if (!strcmp(r->type->mapsto, "MeshType")) {
               if (order == 0) {
                  fprintf(fp, "      CALL MeshCopy(%s1%s%%%s%s, %s_out%s%%%s%s, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )\n", uy, deref, r->name, dimstr(r->ndims)
                     , uy, deref, r->name, dimstr(r->ndims));
               }
               else if (order == 1) {
                  fprintf(fp, "      CALL MeshExtrapInterp1(%s1%s%%%s%s, %s2%s%%%s%s, tin, %s_out%s%%%s%s, tin_out, ErrStat2, ErrMsg2 )\n"
                     , uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
               }
               else if (order == 2) {
                  fprintf(fp, "      CALL MeshExtrapInterp2(%s1%s%%%s%s, %s2%s%%%s%s, %s3%s%%%s%s, tin, %s_out%s%%%s%s, tin_out, ErrStat2, ErrMsg2 )\n"
                     , uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
               }
            }
            else {
               char nonick2[NAMELEN];
               remove_nickname(r->type->module->nickname, r->type->name, nonick2);

               if (order == 0) {
                  fprintf(fp, "      CALL %s_Copy%s(%s1%s%%%s%s, %s_out%s%%%s%s, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )\n", r->type->module->nickname, fast_interface_type_shortname(nonick2)
                     , uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
               }
               else if (order == 1) {
                  fprintf(fp, "      CALL %s_%s_ExtrapInterp1( %s1%s%%%s%s, %s2%s%%%s%s, tin, %s_out%s%%%s%s, tin_out, ErrStat2, ErrMsg2 )\n",
                     r->type->module->nickname, fast_interface_type_shortname(nonick2)
                     , uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
               }
               else if (order == 2) {
                  fprintf(fp, "      CALL %s_%s_ExtrapInterp2( %s1%s%%%s%s, %s2%s%%%s%s, %s3%s%%%s%s, tin, %s_out%s%%%s%s, tin_out, ErrStat2, ErrMsg2 )\n",
                     r->type->module->nickname, fast_interface_type_shortname(nonick2)
                     , uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
               }
            }

            fprintf(fp, "        CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n");
          //fprintf(fp, "        IF (ErrStat>=AbortErrLev) RETURN\n");
            for (j = r->ndims; j >= 1; j--) {
               fprintf(fp, "   ENDDO\n");
            }

         }
      }
      else if (!strcmp(r->type->mapsto, "REAL(ReKi)") ||
         !strcmp(r->type->mapsto, "REAL(SiKi)") ||
         !strcmp(r->type->mapsto, "REAL(R8Ki)") ||
         !strcmp(r->type->mapsto, "REAL(DbKi)")) {


         if (order == 0) {
            //bjj: this should probably have some "IF ALLOCATED" statements around it, but we're just calling
            // the copy routine
            fprintf(fp, "  %s_out%s%%%s = %s1%s%%%s\n", uy, deref, r->name, uy, deref, r->name);
         }
         else
            strcpy(indent, "");
            for (j = r->ndims; j > 0; j--) {
               fprintf(fp, "%s  DO i%d = LBOUND(%s_out%s%%%s,%d),UBOUND(%s_out%s%%%s,%d)\n", indent, j, uy, deref, r->name, j, uy, deref, r->name, j);
               strcat(indent, "  "); //create an indent
            }

            if (order == 1) {
               if (r->gen_periodic) {
                  fprintf(fp, "%s  CALL Angles_ExtrapInterp( %s1%s%%%s%s, %s2%s%%%s%s, tin, %s_out%s%%%s%s, tin_out )\n",
                     indent, uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
               }
               else {
                  fprintf(fp, "%s  b = -(%s1%s%%%s%s - %s2%s%%%s%s)\n", indent, uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
                  fprintf(fp, "%s  %s_out%s%%%s%s = %s1%s%%%s%s + b * ScaleFactor\n", indent, uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
               };
            }
            if (order == 2) {
               if (r->gen_periodic) {
                  fprintf(fp, "%s  CALL Angles_ExtrapInterp( %s1%s%%%s%s, %s2%s%%%s%s, %s3%s%%%s%s, tin, %s_out%s%%%s%s, tin_out )\n",
                     indent, uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
               }
               else {
                  fprintf(fp, "%s  b = (t(3)**2*(%s1%s%%%s%s - %s2%s%%%s%s) + t(2)**2*(-%s1%s%%%s%s + %s3%s%%%s%s))* scaleFactor\n",
                     indent, uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
                  fprintf(fp, "%s  c = ( (t(2)-t(3))*%s1%s%%%s%s + t(3)*%s2%s%%%s%s - t(2)*%s3%s%%%s%s ) * scaleFactor\n",
                     indent, uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
                  fprintf(fp, "%s  %s_out%s%%%s%s = %s1%s%%%s%s + b  + c * t_out\n",
                     indent, uy, deref, r->name, dimstr(r->ndims), uy, deref, r->name, dimstr(r->ndims));
               }
            }
            for (j = r->ndims; j >= 1; j--) {
               strcpy(indent, "");
               for (i = 1; i < j; i++) {
                  strcat(indent, "  ");
               }
               fprintf(fp, "%s  END DO\n", indent);
            }
         }
      // check if this is an allocatable array:
      if (r->ndims > 0 && has_deferred_dim(r, 0)) {
         fprintf(fp, "END IF ! check if allocated\n");
      }
   }

} // gen_extint_order

void calc_extint_order(FILE *fp, const node_t *ModName, node_t *r, int recurselevel, int *max_ndims, int *max_nrecurs, int *max_alloc_ndims) {
   node_t *q, *r1 ;
// bjj: make sure this is consistent with logic of gen_extint_order

   if ( r->type != NULL ) {
   //   if(r->ndims > *max_ndims  )* max_ndims = r->ndims;

      if (r->type->type_type == DERIVED) {
         if ((q = get_entry(make_lower_temp(r->type->name), ModName->module_ddt_list)) != NULL) {
            for (r1 = q->fields; r1; r1 = r1->next)
            {
               if (r->ndims > 0) {
                  if (recurselevel > *max_nrecurs) *max_nrecurs = recurselevel;
                  if (r->ndims     > *max_ndims  ) *max_ndims   = r->ndims;
               }
               calc_extint_order(fp, ModName, r1, recurselevel + 1, max_ndims, max_nrecurs, max_alloc_ndims);
            }
         }
         else if (!strcmp(r->type->mapsto, "MeshType")) {
            if (r->ndims > 0) {
               if (r->ndims > *max_ndims)* max_ndims = r->ndims;
            }
         }
         else {
            if (r->ndims >= 1) {
               if (r->ndims > *max_ndims)* max_ndims = r->ndims;
            }
         }

      }
      else if (!strcmp(r->type->mapsto, "REAL(ReKi)") ||
         !strcmp(r->type->mapsto, "REAL(SiKi)") ||
         !strcmp(r->type->mapsto, "REAL(R8Ki)") ||
         !strcmp(r->type->mapsto, "REAL(DbKi)")) {
         if (/*order > 0 &&*/ r->ndims > *max_alloc_ndims) *max_alloc_ndims = r->ndims;
         if (r->ndims > *max_ndims)* max_ndims = r->ndims;
      }


   }

   if ( recurselevel > MAXRECURSE ) {
     fprintf(stderr,"REGISTRY ERROR: too many levels of array subtypes\n") ;
     exit(9) ;
   }

}

#if 0
void
gen_ExtrapInterp( FILE *fp , const node_t * ModName, char * typnm, char * typnmlong )
{
  char nonick[NAMELEN] ;
  char *ddtname; char uy[2];
  node_t *q, * r ;
  int i, j, max_ndims, max_nrecurs, max_alloc_ndims;

  if (!strcmp(make_lower_temp(typnm), "output")){
     strcpy(uy,"y");
  }
  else{
     strcpy(uy, "u");
  }

  fprintf(fp,"\n") ;
  fprintf(fp," SUBROUTINE %s_%s_ExtrapInterp(%s, tin, %s_out, tin_out, ErrStat, ErrMsg )\n",ModName->nickname,typnm,uy,uy) ;
  fprintf(fp,"!\n") ;
  fprintf(fp, "! This subroutine calculates a extrapolated (or interpolated) %s %s_out at time t_out, from previous/future time\n", typnm, uy);
  fprintf(fp, "! values of %s (which has values associated with times in t).  Order of the interpolation is given by the size of %s\n", uy, uy);
  fprintf(fp,"!\n") ;
  fprintf(fp,"!  expressions below based on either\n") ;
  fprintf(fp,"!\n") ;
  fprintf(fp,"!  f(t) = a\n") ;
  fprintf(fp,"!  f(t) = a + b * t, or\n") ;
  fprintf(fp,"!  f(t) = a + b * t + c * t**2\n") ;
  fprintf(fp,"!\n") ;
  fprintf(fp,"!  where a, b and c are determined as the solution to\n") ;
  fprintf(fp, "!  f(t1) = %s1, f(t2) = %s2, f(t3) = %s3  (as appropriate)\n", uy, uy, uy);
  fprintf(fp,"!\n") ;
  fprintf(fp,"!..................................................................................................................................\n") ;
  fprintf(fp,"\n") ;


  fprintf(fp, " TYPE(%s_%s), INTENT(INOUT)  :: %s(:)      ! %s at t1 > t2 > t3\n", ModName->nickname, typnmlong, uy, typnm);
  fprintf(fp, " REAL(DbKi),         INTENT(IN   )  :: tin(:)      ! Times associated with the %ss\n", typnm);
//jm Modified from INTENT(  OUT) to INTENT(INOUT) to prevent ALLOCATABLE array arguments in the DDT
//jm from being maliciously deallocated through the call.See Sec. 5.1.2.7 of bonehead Fortran 2003 standard
  fprintf(fp, " TYPE(%s_%s), INTENT(INOUT)  :: %s_out     ! %s at tin_out\n", ModName->nickname, typnmlong, uy, typnm);
  fprintf(fp," REAL(DbKi),         INTENT(IN   )  :: tin_out     ! time to be extrap/interp'd to\n") ;
  fprintf(fp," INTEGER(IntKi),     INTENT(  OUT)  :: ErrStat   ! Error status of the operation\n") ;
  fprintf(fp," CHARACTER(*),       INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None\n") ;
  fprintf(fp,"   ! local variables\n") ;
  fprintf(fp, " REAL(DbKi) :: t(SIZE(tin))    ! Times associated with the %ss\n", typnm);
  fprintf(fp," REAL(DbKi) :: t_out           ! Time to which to be extrap/interpd\n") ;
  fprintf(fp," INTEGER(IntKi)                 :: order    ! order of polynomial fit (max 2)\n") ;
  fprintf(fp, " CHARACTER(*),    PARAMETER :: RoutineName = '%s_%s_ExtrapInterp'\n", ModName->nickname, typnm);

  max_ndims   = 0; // ModName->module_ddt_list->max_ndims; //bjj: this is max for module, not for typnmlong
  max_nrecurs = 0; // MAXRECURSE;
  max_alloc_ndims = 0;

  for (q = ModName->module_ddt_list; q; q = q->next)
  {
     if (q->usefrom == 0) {
        ddtname = q->name;
        remove_nickname(ModName->nickname, ddtname, nonick);
        if (!strcmp(nonick, make_lower_temp(typnmlong))) {
           for (r = q->fields; r; r = r->next)
           {
              // recursive
              calc_extint_order(fp, ModName, r, 0, &max_ndims, &max_nrecurs, &max_alloc_ndims);
           }
        }
     }
  }
  //fprintf(stderr, "ndims=%d nrecurs=%d %d\n\n", max_ndims, max_nrecurs, max_alloc_ndims);

  if (max_alloc_ndims >= 0){
  fprintf(fp," REAL(DbKi)                                 :: b0       ! temporary for extrapolation/interpolation\n") ;
  fprintf(fp," REAL(DbKi)                                 :: c0       ! temporary for extrapolation/interpolation\n") ;
  if (max_alloc_ndims >= 1){
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:)        :: b1       ! temporary for extrapolation/interpolation\n") ;
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:)        :: c1       ! temporary for extrapolation/interpolation\n") ;
  if (max_alloc_ndims >= 2){
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:,:)      :: b2       ! temporary for extrapolation/interpolation\n") ;
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:,:)      :: c2       ! temporary for extrapolation/interpolation\n") ;
  if (max_alloc_ndims >= 3){
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:,:,:)    :: b3       ! temporary for extrapolation/interpolation\n") ;
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:,:,:)    :: c3       ! temporary for extrapolation/interpolation\n") ;
  if (max_alloc_ndims >= 4){
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:,:,:,:)  :: b4       ! temporary for extrapolation/interpolation\n") ;
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:,:,:,:)  :: c4       ! temporary for extrapolation/interpolation\n") ;
  if (max_alloc_ndims >= 5){
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:,:,:,:,:):: b5       ! temporary for extrapolation/interpolation\n") ;
  fprintf(fp," REAL(DbKi),ALLOCATABLE,DIMENSION(:,:,:,:,:):: c5       ! temporary for extrapolation/interpolation\n") ;
   } // 5
   } // 4
   } // 3
   } // 2
   } // 1
  } // 0
  fprintf(fp," INTEGER(IntKi)                             :: ErrStat2 ! local errors\n");
  fprintf(fp," CHARACTER(ErrMsgLen)                       :: ErrMsg2  ! local errors\n");
  for ( j = 1 ; j <= max_ndims ; j++ ) {
    for ( i = 0 ; i <= max_nrecurs ; i++ ) {
  fprintf(fp," INTEGER                                    :: i%d%d    ! dim%d level %d counter variable for arrays of ddts\n",i,j,j,i) ;
    }
  }
  fprintf(fp,"    ! Initialize ErrStat\n") ;
  fprintf(fp," ErrStat = ErrID_None\n") ;
  fprintf(fp," ErrMsg  = \"\"\n") ;
  fprintf(fp,"    ! we'll subtract a constant from the times to resolve some \n") ;
  fprintf(fp,"    ! numerical issues when t gets large (and to simplify the equations)\n") ;
  fprintf(fp," t = tin - tin(1)\n") ;
  fprintf(fp," t_out = tin_out - tin(1)\n") ;
  fprintf(fp,"\n") ;
  fprintf(fp, " if ( size(t) .ne. size(%s)) then\n", uy);
  fprintf(fp,"    ErrStat = ErrID_Fatal\n") ;
  fprintf(fp, "    ErrMsg = ' Error in %s_%s_ExtrapInterp: size(t) must equal size(%s) '\n", ModName->nickname, typnm, uy);
  fprintf(fp,"    RETURN\n") ;
  fprintf(fp," endif\n") ;
  fprintf(fp, " if (size(%s) .gt. 3) then\n", uy);
  fprintf(fp,"    ErrStat = ErrID_Fatal\n") ;
  fprintf(fp, "    ErrMsg  = ' Error in %s_%s_ExtrapInterp: size(%s) must be less than 4 '\n", ModName->nickname, typnm, uy);
  fprintf(fp,"    RETURN\n") ;
  fprintf(fp," endif\n") ;

  fprintf(fp, " order = SIZE(%s) - 1\n", uy);

  fprintf(fp," IF ( order .eq. 0 ) THEN\n") ;
  for ( q = ModName->module_ddt_list ; q ; q = q->next )
  {
    if ( q->usefrom == 0 ) {
      ddtname = q->name ;
      remove_nickname(ModName->nickname,ddtname,nonick) ;
      if ( !strcmp( nonick, make_lower_temp(typnmlong) )) {
        for ( r = q->fields ; r ; r = r->next )
        {
          // recursive
          gen_extint_order( fp, ModName, typnm, uy, 0, r, "", 0 ) ;
        }
      }
    }
  }

  fprintf(fp," ELSE IF ( order .eq. 1 ) THEN\n") ;
fprintf(fp,"  IF ( EqualRealNos( t(1), t(2) ) ) THEN\n") ;
fprintf(fp,"    ErrStat = ErrID_Fatal\n") ;
fprintf(fp,"    ErrMsg  = ' Error in %s_%s_ExtrapInterp: t(1) must not equal t(2) to avoid a division-by-zero error.'\n",ModName->nickname,typnm) ;
fprintf(fp,"    RETURN\n") ;
fprintf(fp,"  END IF\n") ;
  for ( q = ModName->module_ddt_list ; q ; q = q->next )
  {

    if ( q->usefrom == 0 ) {
      ddtname = q->name ;
      remove_nickname(ModName->nickname,ddtname,nonick) ;
      if ( !strcmp( nonick, make_lower_temp(typnmlong) )) {
        for ( r = q->fields ; r ; r = r->next )
        {
          // recursive
          gen_extint_order( fp, ModName, typnm, uy, 1, r, "", 0 ) ;
        }
      }
    }
  }
  fprintf(fp," ELSE IF ( order .eq. 2 ) THEN\n") ;
fprintf(fp,"  IF ( EqualRealNos( t(1), t(2) ) ) THEN\n") ;
fprintf(fp,"    ErrStat = ErrID_Fatal\n") ;
fprintf(fp,"    ErrMsg  = ' Error in %s_%s_ExtrapInterp: t(1) must not equal t(2) to avoid a division-by-zero error.'\n",ModName->nickname,typnm) ;
fprintf(fp,"    RETURN\n") ;
fprintf(fp,"  END IF\n") ;
fprintf(fp,"  IF ( EqualRealNos( t(2), t(3) ) ) THEN\n") ;
fprintf(fp,"    ErrStat = ErrID_Fatal\n") ;
fprintf(fp,"    ErrMsg  = ' Error in %s_%s_ExtrapInterp: t(2) must not equal t(3) to avoid a division-by-zero error.'\n",ModName->nickname,typnm) ;
fprintf(fp,"    RETURN\n") ;
fprintf(fp,"  END IF\n") ;
fprintf(fp,"  IF ( EqualRealNos( t(1), t(3) ) ) THEN\n") ;
fprintf(fp,"    ErrStat = ErrID_Fatal\n") ;
fprintf(fp,"    ErrMsg  = ' Error in %s_%s_ExtrapInterp: t(1) must not equal t(3) to avoid a division-by-zero error.'\n",ModName->nickname,typnm) ;
fprintf(fp,"    RETURN\n") ;
fprintf(fp,"  END IF\n") ;

  for ( q = ModName->module_ddt_list ; q ; q = q->next )
  {
    if ( q->usefrom == 0 ) {
      ddtname = q->name ;
      remove_nickname(ModName->nickname,ddtname,nonick) ;
      if ( !strcmp( nonick, make_lower_temp(typnmlong) )) {
        for ( r = q->fields ; r ; r = r->next )
        {
          // recursive
          gen_extint_order( fp, ModName, typnm, uy, 2, r, "", 0 ) ;
        }
      }
    }
  }
  fprintf(fp," ELSE \n") ;
  fprintf(fp,"   ErrStat = ErrID_Fatal\n") ;
  fprintf(fp,"   ErrMsg = ' order must be less than 3 in %s_%s_ExtrapInterp '\n",ModName->nickname,typnm) ;
  fprintf(fp,"   RETURN\n") ;
  fprintf(fp," ENDIF \n") ;


  fprintf(fp," END SUBROUTINE %s_%s_ExtrapInterp\n",ModName->nickname,typnm) ;
  fprintf(fp,"\n") ;
}
#endif

void
gen_ExtrapInterp1(FILE *fp, const node_t * ModName, char * typnm, char * typnmlong, char * xtypnm, char * uy, const int max_ndims, const int max_nrecurs, const int max_alloc_ndims, const node_t *q)
{
   node_t *r;
   int i, j;

   fprintf(fp, "\n");
   fprintf(fp, " SUBROUTINE %s_%s_ExtrapInterp1(%s1, %s2, tin, %s_out, tin_out, ErrStat, ErrMsg )\n", ModName->nickname, typnm, uy, uy, uy);
   fprintf(fp, "!\n");
   fprintf(fp, "! This subroutine calculates a extrapolated (or interpolated) %s %s_out at time t_out, from previous/future time\n", typnm, uy);
   fprintf(fp, "! values of %s (which has values associated with times in t).  Order of the interpolation is 1.\n", uy);
   fprintf(fp, "!\n");
   fprintf(fp, "!  f(t) = a + b * t, or\n");
   fprintf(fp, "!\n");
   fprintf(fp, "!  where a and b are determined as the solution to\n");
   fprintf(fp, "!  f(t1) = %s1, f(t2) = %s2\n", uy, uy);
   fprintf(fp, "!\n");
   fprintf(fp, "!..................................................................................................................................\n");
   fprintf(fp, "\n");


   fprintf(fp, " TYPE(%s_%s), INTENT(%s)  :: %s1    ! %s at t1 > t2\n", ModName->nickname, typnmlong, (q->containsPtr == 1) ? "INOUT" : "IN", uy, typnm);
   fprintf(fp, " TYPE(%s_%s), INTENT(%s)  :: %s2    ! %s at t2 \n", ModName->nickname, typnmlong, (q->containsPtr == 1) ? "INOUT" : "IN", uy, typnm);
   fprintf(fp, " REAL(%s),         INTENT(IN   )          :: tin(2)   ! Times associated with the %ss\n", xtypnm, typnm);
   fprintf(fp, " TYPE(%s_%s), INTENT(INOUT)  :: %s_out ! %s at tin_out\n", ModName->nickname, typnmlong, uy, typnm);
   fprintf(fp, " REAL(%s),         INTENT(IN   )          :: tin_out  ! time to be extrap/interp'd to\n", xtypnm);
   fprintf(fp, " INTEGER(IntKi),     INTENT(  OUT)          :: ErrStat  ! Error status of the operation\n");
   fprintf(fp, " CHARACTER(*),       INTENT(  OUT)          :: ErrMsg   ! Error message if ErrStat /= ErrID_None\n");
   fprintf(fp, "   ! local variables\n");
   fprintf(fp, " REAL(%s)                                 :: t(2)     ! Times associated with the %ss\n", xtypnm, typnm);
   fprintf(fp, " REAL(%s)                                 :: t_out    ! Time to which to be extrap/interpd\n", xtypnm);
   fprintf(fp, " CHARACTER(*),                    PARAMETER :: RoutineName = '%s_%s_ExtrapInterp1'\n", ModName->nickname, typnm);


   fprintf(fp, " REAL(DbKi)                                 :: b        ! temporary for extrapolation/interpolation\n");
   fprintf(fp, " REAL(DbKi)                                 :: ScaleFactor ! temporary for extrapolation/interpolation\n");
   fprintf(fp, " INTEGER(IntKi)                             :: ErrStat2 ! local errors\n");
   fprintf(fp, " CHARACTER(ErrMsgLen)                       :: ErrMsg2  ! local errors\n");
   for (j = 1; j <= max_ndims; j++) {
      for (i = 0; i <= max_nrecurs; i++) {
         fprintf(fp, " INTEGER                                    :: i%d%d    ! dim%d level %d counter variable for arrays of ddts\n", i, j, j, i);
      }
   }
   for (j = 1; j <= max_ndims; j++) {
      fprintf(fp, " INTEGER                                    :: i%d    ! dim%d counter variable for arrays\n", j, j);
   }

   fprintf(fp, "    ! Initialize ErrStat\n");
   fprintf(fp, " ErrStat = ErrID_None\n");
   fprintf(fp, " ErrMsg  = \"\"\n");
   fprintf(fp, "    ! we'll subtract a constant from the times to resolve some \n");
   fprintf(fp, "    ! numerical issues when t gets large (and to simplify the equations)\n");
   fprintf(fp, " t = tin - tin(1)\n");
   fprintf(fp, " t_out = tin_out - tin(1)\n");
   fprintf(fp, "\n");

   fprintf(fp, "   IF ( EqualRealNos( t(1), t(2) ) ) THEN\n");
   fprintf(fp, "     CALL SetErrStat(ErrID_Fatal, 't(1) must not equal t(2) to avoid a division-by-zero error.', ErrStat, ErrMsg,RoutineName)\n");
   fprintf(fp, "     RETURN\n");
   fprintf(fp, "   END IF\n\n");

   fprintf(fp, "   ScaleFactor = t_out / t(2)\n");

   for (r = q->fields; r; r = r->next)
   {
      // recursive
      gen_extint_order(fp, ModName, typnm, uy, 1, r, "", 0);
   }


   fprintf(fp, " END SUBROUTINE %s_%s_ExtrapInterp1\n", ModName->nickname, typnm);
   fprintf(fp, "\n");
}

void
gen_ExtrapInterp2(FILE *fp, const node_t * ModName, char * typnm, char * typnmlong, char * xtypnm, char * uy, const int max_ndims, const int max_nrecurs, const int max_alloc_ndims, const node_t *q)
{
   node_t *r;
   int i, j;

   fprintf(fp, "\n");
   fprintf(fp, " SUBROUTINE %s_%s_ExtrapInterp2(%s1, %s2, %s3, tin, %s_out, tin_out, ErrStat, ErrMsg )\n", ModName->nickname, typnm, uy, uy, uy, uy);
   fprintf(fp, "!\n");
   fprintf(fp, "! This subroutine calculates a extrapolated (or interpolated) %s %s_out at time t_out, from previous/future time\n", typnm, uy);
   fprintf(fp, "! values of %s (which has values associated with times in t).  Order of the interpolation is 2.\n", uy);
   fprintf(fp, "!\n");
   fprintf(fp, "!  expressions below based on either\n");
   fprintf(fp, "!\n");
   fprintf(fp, "!  f(t) = a + b * t + c * t**2\n");
   fprintf(fp, "!\n");
   fprintf(fp, "!  where a, b and c are determined as the solution to\n");
   fprintf(fp, "!  f(t1) = %s1, f(t2) = %s2, f(t3) = %s3\n", uy, uy, uy);
   fprintf(fp, "!\n");
   fprintf(fp, "!..................................................................................................................................\n");
   fprintf(fp, "\n");

   fprintf(fp, " TYPE(%s_%s), INTENT(%s)  :: %s1      ! %s at t1 > t2 > t3\n", ModName->nickname, typnmlong, (q->containsPtr == 1) ? "INOUT" : "IN", uy, typnm);
   fprintf(fp, " TYPE(%s_%s), INTENT(%s)  :: %s2      ! %s at t2 > t3\n", ModName->nickname, typnmlong, (q->containsPtr == 1) ? "INOUT" : "IN", uy, typnm);
   fprintf(fp, " TYPE(%s_%s), INTENT(%s)  :: %s3      ! %s at t3\n", ModName->nickname, typnmlong, (q->containsPtr == 1) ? "INOUT" : "IN", uy, typnm);
   fprintf(fp, " REAL(%s),                 INTENT(IN   )  :: tin(3)    ! Times associated with the %ss\n", xtypnm, typnm);
   fprintf(fp, " TYPE(%s_%s), INTENT(INOUT)  :: %s_out     ! %s at tin_out\n", ModName->nickname, typnmlong, uy, typnm);
   fprintf(fp, " REAL(%s),                 INTENT(IN   )  :: tin_out   ! time to be extrap/interp'd to\n", xtypnm);

   fprintf(fp, " INTEGER(IntKi),             INTENT(  OUT)  :: ErrStat   ! Error status of the operation\n" );
   fprintf(fp, " CHARACTER(*),               INTENT(  OUT)  :: ErrMsg    ! Error message if ErrStat /= ErrID_None\n");
   fprintf(fp, "   ! local variables\n");
   fprintf(fp, " REAL(%s)                                 :: t(3)      ! Times associated with the %ss\n", xtypnm, typnm);
   fprintf(fp, " REAL(%s)                                 :: t_out     ! Time to which to be extrap/interpd\n", xtypnm);
   fprintf(fp, " INTEGER(IntKi)                             :: order     ! order of polynomial fit (max 2)\n");

   fprintf(fp, " REAL(DbKi)                                 :: b        ! temporary for extrapolation/interpolation\n");
   fprintf(fp, " REAL(DbKi)                                 :: c        ! temporary for extrapolation/interpolation\n");
   fprintf(fp, " REAL(DbKi)                                 :: ScaleFactor ! temporary for extrapolation/interpolation\n");
   fprintf(fp, " INTEGER(IntKi)                             :: ErrStat2 ! local errors\n");
   fprintf(fp, " CHARACTER(ErrMsgLen)                       :: ErrMsg2  ! local errors\n");
   fprintf(fp, " CHARACTER(*),            PARAMETER         :: RoutineName = '%s_%s_ExtrapInterp2'\n", ModName->nickname, typnm);
   for (j = 1; j <= max_ndims; j++) {
      for (i = 0; i <= max_nrecurs; i++) {
         fprintf(fp, " INTEGER                                    :: i%d%d    ! dim%d level %d counter variable for arrays of ddts\n", i, j, j, i);
      }
   }
   for (j = 1; j <= max_ndims; j++) {
      fprintf(fp, " INTEGER                                    :: i%d    ! dim%d counter variable for arrays\n", j, j);
   }
   fprintf(fp, "    ! Initialize ErrStat\n");
   fprintf(fp, " ErrStat = ErrID_None\n");
   fprintf(fp, " ErrMsg  = \"\"\n");
   fprintf(fp, "    ! we'll subtract a constant from the times to resolve some \n");
   fprintf(fp, "    ! numerical issues when t gets large (and to simplify the equations)\n");
   fprintf(fp, " t = tin - tin(1)\n");
   fprintf(fp, " t_out = tin_out - tin(1)\n");
   fprintf(fp, "\n");


   fprintf(fp, "   IF ( EqualRealNos( t(1), t(2) ) ) THEN\n");
   fprintf(fp, "     CALL SetErrStat(ErrID_Fatal, 't(1) must not equal t(2) to avoid a division-by-zero error.', ErrStat, ErrMsg,RoutineName)\n");
   fprintf(fp, "     RETURN\n");
   fprintf(fp, "   ELSE IF ( EqualRealNos( t(2), t(3) ) ) THEN\n");
   fprintf(fp, "     CALL SetErrStat(ErrID_Fatal, 't(2) must not equal t(3) to avoid a division-by-zero error.', ErrStat, ErrMsg,RoutineName)\n");
   fprintf(fp, "     RETURN\n");
   fprintf(fp, "   ELSE IF ( EqualRealNos( t(1), t(3) ) ) THEN\n");
   fprintf(fp, "     CALL SetErrStat(ErrID_Fatal, 't(1) must not equal t(3) to avoid a division-by-zero error.', ErrStat, ErrMsg,RoutineName)\n");
   fprintf(fp, "     RETURN\n");
   fprintf(fp, "   END IF\n\n");

   fprintf(fp, "   ScaleFactor = t_out / (t(2) * t(3) * (t(2) - t(3)))\n");

   

   for (r = q->fields; r; r = r->next)
   {
      // recursive
      gen_extint_order(fp, ModName, typnm, uy, 2, r, "", 0);
   }


   fprintf(fp, " END SUBROUTINE %s_%s_ExtrapInterp2\n", ModName->nickname, typnm);
   fprintf(fp, "\n");
}


void
gen_ExtrapInterp(FILE *fp, const node_t * ModName, char * typnm, char * typnmlong, char * xtypnm)
{
   char nonick[NAMELEN];
   char *ddtname; char uy[2];
   node_t *q, *r;
   int max_ndims, max_nrecurs, max_alloc_ndims;

   if (!strcmp(make_lower_temp(typnm), "output")){
      strcpy(uy, "y");
   }
   else{
      strcpy(uy, "u");
   }

   for (q = ModName->module_ddt_list; q; q = q->next)
   {
      if (q->usefrom == 0) {
         ddtname = q->name;
         remove_nickname(ModName->nickname, ddtname, nonick);
         if (!strcmp(nonick, make_lower_temp(typnmlong))) {

            fprintf(fp, "\n");
            fprintf(fp, " SUBROUTINE %s_%s_ExtrapInterp(%s, t, %s_out, t_out, ErrStat, ErrMsg )\n", ModName->nickname, typnm, uy, uy);
            fprintf(fp, "!\n");
            fprintf(fp, "! This subroutine calculates a extrapolated (or interpolated) %s %s_out at time t_out, from previous/future time\n", typnm, uy);
            fprintf(fp, "! values of %s (which has values associated with times in t).  Order of the interpolation is given by the size of %s\n", uy, uy);
            fprintf(fp, "!\n");
            fprintf(fp, "!  expressions below based on either\n");
            fprintf(fp, "!\n");
            fprintf(fp, "!  f(t) = a\n");
            fprintf(fp, "!  f(t) = a + b * t, or\n");
            fprintf(fp, "!  f(t) = a + b * t + c * t**2\n");
            fprintf(fp, "!\n");
            fprintf(fp, "!  where a, b and c are determined as the solution to\n");
            fprintf(fp, "!  f(t1) = %s1, f(t2) = %s2, f(t3) = %s3  (as appropriate)\n", uy, uy, uy);
            fprintf(fp, "!\n");
            fprintf(fp, "!..................................................................................................................................\n");
            fprintf(fp, "\n");


            fprintf(fp, " TYPE(%s_%s), INTENT(%s)  :: %s(:) ! %s at t1 > t2 > t3\n", ModName->nickname, typnmlong, (q->containsPtr == 1) ? "INOUT" : "IN", uy, typnm);
            fprintf(fp, " REAL(%s),                 INTENT(IN   )  :: t(:)           ! Times associated with the %ss\n", xtypnm, typnm);
            //jm Modified from INTENT(  OUT) to INTENT(INOUT) to prevent ALLOCATABLE array arguments in the DDT
            //jm from being maliciously deallocated through the call.See Sec. 5.1.2.7 of bonehead Fortran 2003 standard
            fprintf(fp, " TYPE(%s_%s), INTENT(INOUT)  :: %s_out ! %s at tin_out\n", ModName->nickname, typnmlong, uy, typnm);
            fprintf(fp, " REAL(%s),                 INTENT(IN   )  :: t_out           ! time to be extrap/interp'd to\n", xtypnm);
            fprintf(fp, " INTEGER(IntKi),             INTENT(  OUT)  :: ErrStat         ! Error status of the operation\n");
            fprintf(fp, " CHARACTER(*),               INTENT(  OUT)  :: ErrMsg          ! Error message if ErrStat /= ErrID_None\n");
            fprintf(fp, "   ! local variables\n");
            fprintf(fp, " INTEGER(IntKi)                             :: order           ! order of polynomial fit (max 2)\n");
            fprintf(fp, " INTEGER(IntKi)                             :: ErrStat2        ! local errors\n");
            fprintf(fp, " CHARACTER(ErrMsgLen)                       :: ErrMsg2         ! local errors\n");
            fprintf(fp, " CHARACTER(*),    PARAMETER                 :: RoutineName = '%s_%s_ExtrapInterp'\n", ModName->nickname, typnm);
            fprintf(fp, "    ! Initialize ErrStat\n");
            fprintf(fp, " ErrStat = ErrID_None\n");
            fprintf(fp, " ErrMsg  = \"\"\n");
            fprintf(fp, " if ( size(t) .ne. size(%s)) then\n", uy);
            fprintf(fp, "    CALL SetErrStat(ErrID_Fatal,'size(t) must equal size(%s)',ErrStat,ErrMsg,RoutineName)\n",uy);
            fprintf(fp, "    RETURN\n");
            fprintf(fp, " endif\n");

            fprintf(fp, " order = SIZE(%s) - 1\n", uy);

            fprintf(fp, " IF ( order .eq. 0 ) THEN\n");
            fprintf(fp, "   CALL %s_Copy%s(%s(1), %s_out, MESH_UPDATECOPY, ErrStat2, ErrMsg2 )\n", ModName->nickname, typnm, uy, uy);
            fprintf(fp, "     CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n");
            fprintf(fp, " ELSE IF ( order .eq. 1 ) THEN\n");
            fprintf(fp, "   CALL %s_%s_ExtrapInterp1(%s(1), %s(2), t, %s_out, t_out, ErrStat2, ErrMsg2 )\n", ModName->nickname, typnm, uy, uy, uy);
            fprintf(fp, "     CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n");
            fprintf(fp, " ELSE IF ( order .eq. 2 ) THEN\n");
            fprintf(fp, "   CALL %s_%s_ExtrapInterp2(%s(1), %s(2), %s(3), t, %s_out, t_out, ErrStat2, ErrMsg2 )\n", ModName->nickname, typnm, uy, uy, uy, uy);
            fprintf(fp, "     CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg,RoutineName)\n");
            fprintf(fp, " ELSE \n");
            fprintf(fp, "   CALL SetErrStat(ErrID_Fatal,'size(%s) must be less than 4 (order must be less than 3).',ErrStat,ErrMsg,RoutineName)\n", uy);
            fprintf(fp, "   RETURN\n");
            fprintf(fp, " ENDIF \n");

            fprintf(fp, " END SUBROUTINE %s_%s_ExtrapInterp\n", ModName->nickname, typnm);
            fprintf(fp, "\n");


            max_ndims = 0; // ModName->module_ddt_list->max_ndims; //bjj: this is max for module, not for typnmlong
            max_nrecurs = 0; // MAXRECURSE;
            max_alloc_ndims = 0;

            for (r = q->fields; r; r = r->next)
            {
               // recursive
               calc_extint_order(fp, ModName, r, 0, &max_ndims, &max_nrecurs, &max_alloc_ndims);
            }

            gen_ExtrapInterp1(fp, ModName, typnm, typnmlong, xtypnm, uy, max_ndims, max_nrecurs, max_alloc_ndims, q);
            gen_ExtrapInterp2(fp, ModName, typnm, typnmlong, xtypnm, uy, max_ndims, max_nrecurs, max_alloc_ndims, q);

         }
      }
   }



}







void
gen_rk4( FILE *fp , const node_t * ModName )
{
  char nonick[NAMELEN] ;
  char *ddtname ;
  node_t *q, * r ;
  int founddt, k ;

// make sure the user has dt in their parameter types
  founddt = 0 ;
  for ( q = ModName->module_ddt_list ; q ; q = q->next )
  {
    if ( q->usefrom == 0 ) {
      ddtname = q->name ;
      remove_nickname(ModName->nickname,ddtname,nonick) ;
      if ( !strcmp( nonick, "parametertype")) {
        for ( r = q->fields ; r ; r = r->next )
        {
          if ( !strcmp( r->type->mapsto, "REAL(ReKi)") ||
               !strcmp( r->type->mapsto, "REAL(SiKi)") ||
               !strcmp( r->type->mapsto, "REAL(R8Ki)") ||
               !strcmp( r->type->mapsto, "REAL(DbKi)"))
          {
            if ( !strcmp(make_lower_temp(r->name),"dt") ) {
              founddt = 1 ;
            }
          }
        }
      }
    }
  }
  if ( !founddt ) {
    fprintf(stderr,"Registry warning: cannot generate %s_RK4. Add dt to ParameterType for this module\n", ModName->nickname) ;
    return ;
  }


  fprintf(fp," SUBROUTINE %s_RK4(t, u, u_next, p, x, xd, z, OtherState, m, xdot, ErrStat, ErrMsg )\n",
                                                                                              ModName->nickname) ;
  fprintf(fp,"  REAL(DbKi),                   INTENT(IN   ) :: t           ! Current simulation time in seconds\n") ;
  fprintf(fp,"  TYPE(%s_InputType),           INTENT(IN   ) :: u           ! Inputs at t\n",  ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_InputType),           INTENT(IN   ) :: u_next      ! Inputs at t\n",  ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ParameterType),       INTENT(IN   ) :: p           ! Parameters\n",  ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousStateType), INTENT(INOUT) :: x           ! Continuous states at t on input at t + dt on output\n",
                                                                                              ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_DiscreteStateType),   INTENT(INOUT) :: xd          ! Discrete states at t\n",   ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ConstraintStateType), INTENT(IN   ) :: z           ! Constraint states at t (possibly a guess)\n",
                                                                                              ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_OtherStateType),      INTENT(INOUT) :: OtherState  ! Other states\n",  ModName->nickname) ;
  fprintf(fp, "  TYPE(%s_MiscVarType),         INTENT(INOUT) :: m           ! Misc/optimization variables\n", ModName->nickname);
  fprintf(fp, "  TYPE(%s_ContinuousStateType), INTENT(IN   ) :: xdot        ! Continuous states at t on input at t + dt on output\n",
                                                                                              ModName->nickname) ;
  fprintf(fp,"  INTEGER(IntKi),               INTENT(  OUT) :: ErrStat\n") ;
  fprintf(fp,"  CHARACTER(*),                 INTENT(  OUT) :: ErrMsg\n") ;
  fprintf(fp,"    ! Local variables\n" ) ;
  fprintf(fp,"  TYPE(%s_ContinuousStateType)                :: xdot_local     ! t derivatives of continuous states\n",
                                                                                             ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousStateType)                :: k1\n",
                                                                                             ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousStateType)                :: k2\n",
                                                                                             ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousStateType)                :: k3\n",
                                                                                             ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousStateType)                :: k4\n",
                                                                                             ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_ContinuousStateType)                :: x_tmp       ! Holds temporary modification to x\n",
                                                                                             ModName->nickname) ;
  fprintf(fp,"  TYPE(%s_InputType)                          :: u_interp\n",
                                                                                             ModName->nickname) ;
  fprintf(fp,"  REAL(ReKi)                                  :: alpha\n") ;

  fprintf(fp,"    ! Initialize ErrStat\n") ;

  fprintf(fp,"  ErrStat = ErrID_None\n") ;
  fprintf(fp,"  ErrMsg  = \"\"\n") ;
  fprintf(fp," !CALL %s_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, xdot_local, ErrStat, ErrMsg )\n",
                                                                                             ModName->nickname) ;
  fprintf(fp,"  alpha = 0.5\n") ;
  for ( k = 1 ; k <= 4 ; k++ )
  {
// generate statements for k1
  for ( q = ModName->module_ddt_list ; q ; q = q->next )
  {
    if ( q->usefrom == 0 ) {
      ddtname = q->name ;
      remove_nickname(ModName->nickname,ddtname,nonick) ;
      if ( !strcmp( nonick, "continuousstatetype")) {
        for ( r = q->fields ; r ; r = r->next )
        {
          if ( !strcmp( r->type->mapsto, "REAL(ReKi)") ||
             !strcmp(r->type->mapsto, "REAL(SiKi)") ||
             !strcmp(r->type->mapsto, "REAL(R8Ki)") ||
             !strcmp(r->type->mapsto, "REAL(DbKi)"))
          {
  fprintf(fp,"  k%d%%%s = p%%dt * xdot%s%%%s\n",k,r->name,(k<2)?"":"_local",r->name) ;
          }
        }
      }
    }
  }
// generate statements for x_tmp
  for ( q = ModName->module_ddt_list ; q ; q = q->next )
  {
    if ( q->usefrom == 0 ) {
      ddtname = q->name ;
      remove_nickname(ModName->nickname,ddtname,nonick) ;
      if ( !strcmp( nonick, "continuousstatetype")) {
        for ( r = q->fields ; r ; r = r->next )
        {
          if ( !strcmp( r->type->mapsto, "REAL(ReKi)") ||
             !strcmp(r->type->mapsto, "REAL(SiKi)") ||
             !strcmp(r->type->mapsto, "REAL(R8Ki)") ||
             !strcmp(r->type->mapsto, "REAL(DbKi)"))
          {
            if ( k < 4 ) {
  fprintf(fp,"  x_tmp%%%s = x%%%s + %s k%d%%%s\n",r->name,r->name,(k<3)?"0.5*":"",k,r->name) ;
            } else {
  fprintf(fp,"  x%%%s = x%%%s + ( k1%%%s + 2. * k2%%%s + 2. * k3%%%s  + k4%%%s ) / 6.\n",r->name,r->name,r->name,r->name,r->name,r->name) ;
            }
          }
        }
      }
    }
  }

  if (k == 1)  fprintf(fp,"  CALL %s_LinearInterpInput(u, u_next, u_interp, alpha, ErrStat, ErrMsg)\n",
                                                                                             ModName->nickname) ;
  if (k < 4 )fprintf(fp,"  CALL %s_CalcContStateDeriv( t+%sp%%dt, u_%s, p, x_tmp, xd, z, OtherState, m, xdot_local, ErrStat, ErrMsg )\n",
                                                                                             ModName->nickname,
                                                                                             (k<3)?"0.5*":"",
                                                                                             (k<3)?"interp":"next") ;
  fprintf(fp,"\n") ;
  }
  fprintf(fp," END SUBROUTINE %s_RK4\n",ModName->nickname) ;


}


void
gen_module( FILE * fp , node_t * ModName, char * prog_ver )
{
  node_t * p, * q, * r ;
  int i ;
  int ipass ;
  char nonick[NAMELEN] ;
  char tmp[NAMELEN] ;
  char ** p1;

  if ( strlen(ModName->nickname) > 0 ) {
// gen preamble
    {
      fprintf( fp, "! %s\n", prog_ver );

      for ( p1 = FAST_preamble ; *p1 ; p1++ ) { fprintf( fp, *p1, ModName->name ) ; }
    }
    for ( p = ModNames ; p ; p = p->next )
    {
      // Add use declarations for Modules that are included as "usefrom"
      if ( p->usefrom == 1 ) {
        if ( strcmp(make_lower_temp(p->name),"nwtc_library") ) {
          fprintf(fp,"USE %s_Types\n",p->name) ;
        }
      }
    }
    if ( sw_ccode ) {
// Generate a container object for the Fortran code to carry around a pointer to the CPP object(s)
      //fprintf(fp,"USE %s_C_Types\n",ModName->nickname) ;
      fprintf(fp,"!USE, INTRINSIC :: ISO_C_Binding\n") ; // this is inherited from NWTC_Library.f90, and older versions of gfortran complain about ambiguous data when we use this (it thinks it's declared twice; see http://gcc.gnu.org/ml/fortran/2013-04/msg00166.html )
    }

// if this is the NWTC Library, we're not going to print "USE NWTC_Library"
    if ( strcmp(make_lower_temp(ModName->name),"nwtc_library") == 0 ) {
      fprintf(fp,"USE SysSubs\n");
    } else {
      fprintf(fp,"USE NWTC_Library\n");
    }

    fprintf(fp,"IMPLICIT NONE\n") ;

#if 0
    if ( sw_ccode ) {
      fprintf(fp,"  TYPE MAP_In_C \n") ;
      fprintf(fp,"  ! This allows us to create an instance of a C++ \n") ;
      fprintf(fp,"  ! object in Fortran. From the perspective of \n") ;
      fprintf(fp,"  ! Fortran, this is seen as an address in memory\n") ;
      fprintf(fp,"     PRIVATE\n") ;
      fprintf(fp,"     TYPE(C_ptr) :: %s_UserData = C_NULL_ptr\n",ModName->nickname) ;
      fprintf(fp,"  END TYPE MAP_In_C \n") ;
    }
#endif

// generate parameters
    for ( q = ModName->params ; q ; q = q->next )
    {
      fprintf(fp,"    %s, PUBLIC, PARAMETER ",q->type->mapsto ) ;
      if ( q->ndims > 0 )
      {
        if ( q->dims[0]->deferred )
        {
          fprintf(stderr,"Registry warning: parameter %s can not have deferred type\n",q->name) ;
          fprintf(fp,"), ALLOCATABLE ") ;
        } else {
          fprintf(fp,", DIMENSION(") ;
          for ( i = 0 ; i < q->ndims ; i++ )
          {
            fprintf(fp,"%d:%d",q->dims[i]->coord_start,q->dims[i]->coord_end) ;
            if ( i < q->ndims-1 ) fprintf(fp,",") ;
          }
          fprintf(fp,") ") ;
        }
      }
      if ( strlen(q->inival) > 0 ) {
        if ( q->ndims > 0 ) {
          fprintf(fp," :: %s = (/%s/)", q->name, q->inival ) ;
        } else {
          fprintf(fp," :: %s = %s ", q->name, q->inival ) ;
        }
      } else {
        fprintf(fp," :: %s",q->name) ;
      }
      if ( strcmp( q->descrip, "-" ) || strcmp( q->units, "-" ) ) /* that is, if not equal "-" */ {
         fprintf(fp,"     ! %s [%s]", q->descrip, q->units) ;
      }
      fprintf(fp,"\n") ;
    }

// generate each derived data type
    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      if (*q->mapsto) remove_nickname( ModName->nickname, make_lower_temp(q->mapsto) , nonick ) ;
      fprintf(fp, "! =========  %s%s  =======\n", q->mapsto, (sw_ccode) ? "_C" : "");
    for ( ipass = (sw_ccode)?0:1 ; ipass < 2 ; ipass++ ) {   // 2 passes for C code, 1st pass generates bound ddt
      if ( q->usefrom == 0 ) {
        fprintf(fp,"  TYPE, %s :: %s%s\n",(ipass==0)?"BIND(C)":"PUBLIC",q->mapsto,(ipass==0)?"_C":"") ;
        if ( sw_ccode ) {
          if ( ipass == 0 ) {
//            q->containsPtr = 1;
             //if (!strcmp(make_lower_temp(nonick), "otherstatetype") || !strcmp(make_lower_temp(nonick), "initinputtype")){
                fprintf(fp, "   TYPE(C_PTR) :: object = C_NULL_PTR\n");
             //}
          } else {
            fprintf(fp,"    TYPE( %s_C ) :: C_obj\n",q->mapsto) ;
          }
        }
        for ( r = q->fields ; r ; r = r->next )
        {
          if ( r->type != NULL ) {
              // check max number of dimmensions
              // check if this type contains any pointers/meshes or types that have pointers/meshes
              if (r->ndims > q->max_ndims) q->max_ndims = r->ndims;
              if (r->ndims > ModName->module_ddt_list->max_ndims) ModName->module_ddt_list->max_ndims = r->ndims;
           if ( ipass == 0 ) {
              //r->containsPtr = 1;
              //q->containsPtr = 1;
              if        ( r->ndims == 0 && r->type->type_type != DERIVED ) {
                fprintf(fp,"    %s :: %s \n",c_types_binding( r->type->mapsto), r->name) ;
              } else if ( r->ndims >  0 && r->type->type_type != DERIVED ) {
                 if (r->dims[0]->deferred ) {
                  fprintf(fp,"    TYPE(C_ptr) :: %s = C_NULL_PTR \n", r->name) ;
                  fprintf(fp,"    INTEGER(C_int) :: %s_Len = 0 \n", r->name) ;
                 }
                 else {
                    if (strcmp(C_type(r->type->mapsto), "char")){
                  fprintf(fp,"    TYPE(C_PTR) :: %s(", r->name) ;
                  for ( i = 0 ; i < r->ndims ; i++ )
                  {
                    fprintf(fp,"%d",r->dims[i]->coord_end) ;
                    if ( i < r->ndims-1 ) fprintf(fp,",") ;
                  }
                  fprintf(fp,")\n") ;
                    }

                }
              }
           } else { // ipass /= 0
            if ( r->type->type_type == DERIVED ) {
               fprintf(fp,"    TYPE(%s) ",r->type->mapsto ) ;

               checkContainsMesh(r);
               if (r->containsPtr) q->containsPtr = 1;

               // bjj: we need to make sure these types map to reals, too
               tmp[0] = '\0' ;
               if (*q->mapsto ) remove_nickname( ModName->nickname, make_lower_temp(q->mapsto) , tmp ) ;
               if ( must_have_real_or_double(tmp) ) checkOnlyReals( q->mapsto, r );


            } else {
              tmp[0] = '\0' ;
              if (*q->mapsto ) remove_nickname( ModName->nickname, make_lower_temp(q->mapsto) , tmp ) ;
              if ( must_have_real_or_double(tmp) ) {
                if ( strncmp(r->type->mapsto,"REAL",4) ) {
                  fprintf(stderr,"Registry warning: %s contains a field (%s) whose type is not real or double: %s\n",
                   q->mapsto, r->name , r->type->mapsto ) ;
                }

              }
               if ( is_pointer(r) ) {
                  fprintf(fp,"    %s ",c_types_binding(r->type->mapsto) ) ;
               } else {
                  fprintf(fp,"    %s ",r->type->mapsto ) ;
               }
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
                  if ( is_pointer(r) ) {
                  fprintf(fp,"), POINTER ") ;
                  } else {
                  fprintf(fp,"), ALLOCATABLE ") ;
                  }

                } else {
                  fprintf(fp,", DIMENSION(") ;
                  for ( i = 0 ; i < r->ndims ; i++ )
                  {
                     if (r->dims[i]->dim_param == 0){
                        fprintf(fp, "%d:%d", r->dims[i]->coord_start, r->dims[i]->coord_end) ;
                     }
                     else {
                        //fprintf(stderr, "start, %s, %s, %s\n", dimspec, dim_entry->name, dim_entry->module);
                       // if (r->module != NULL) { node_t *param_dim = get_entry(r->dims[i]->dim_param_name, r->module->params); }

                        fprintf(fp, "%s", r->dims[i]->dim_param_name);
                     }
                    if ( i < r->ndims-1 ) fprintf(fp,",") ;
                  }
                  fprintf(fp,") ") ;
                }
            }


            if ( is_pointer( r ) ) {
              fprintf(fp," :: %s => NULL() ",r->name) ;
            } else if  ( r->ndims == 0 && strlen(r->inival) > 0 ) {
              fprintf(fp," :: %s = %s ", r->name, r->inival ) ;
            } else {
              fprintf(fp," :: %s ",r->name) ;
            }

            if ( strcmp( r->descrip, "-" ) || strcmp( r->units, "-" ) ) /* that is, if not equal "-" */ {
               fprintf(fp,"     !< %s [%s]", r->descrip, r->units) ;
            }
            fprintf(fp,"\n") ;
           } // ipass /= 0
          }
        }
        fprintf(fp,"  END TYPE %s%s\n",q->mapsto,(ipass==0)?"_C":"") ;
        //fprintf(stderr, "module %d type %d\n", ModName->module_ddt_list->max_ndims, q->max_ndims);

      }
  }
      fprintf(fp,"! =======================\n") ;
    }

    if ( sw_ccode ) {
      for ( q = ModName->module_ddt_list ; q ; q = q->next )
      {

         if ( q->usefrom == 0 ) {

            char * ddtname, * ddtnamelong, nonick[NAMELEN] ;
            ddtname = q->name ;

            remove_nickname(ModName->nickname,ddtname,nonick) ;

            if ( is_a_fast_interface_type( nonick ) ) {
               ddtnamelong = nonick ;
               ddtname = fast_interface_type_shortname( nonick ) ;
            } else {
               ddtnamelong = ddtname ;
            }

         }
      }
    } // sw_ccode


    fprintf(fp,"CONTAINS\n") ;
    for ( q = ModName->module_ddt_list ; q ; q = q->next )
    {
      if ( q->usefrom == 0 ) {

        char * ddtname, * ddtnamelong, nonick[NAMELEN] ;
        //ddtname = q->name ;
		ddtname = q->mapsto;

        remove_nickname(ModName->nickname,ddtname,nonick) ;

//fprintf(stderr,">> %s %s %s \n",ModName->name, ddtname, nonick) ;

        if ( is_a_fast_interface_type( nonick ) ) {
          ddtnamelong = nonick ;
          ddtname = fast_interface_type_shortname( nonick ) ;
        } else {
          ddtnamelong = ddtname ;
        }

        gen_copy( fp, ModName, ddtname, ddtnamelong , q) ;
        gen_destroy( fp, ModName, ddtname, ddtnamelong ) ;
        gen_pack( fp, ModName, ddtname, ddtnamelong ) ;
        gen_unpack( fp, ModName, ddtname, ddtnamelong ) ;
        if ( sw_ccode ) {
            gen_copy_c2f( fp, ModName, ddtname, ddtnamelong ) ;
            gen_copy_f2c(fp, ModName, ddtname, ddtnamelong);
        }

      }
    }
// bjj: removed gen_modname_pack and gen_modname_unpack because i don't see them being used any differently than the other pack/unpack routines 02/22/2014
//    gen_modname_pack( fp, ModName ) ;
//    gen_modname_unpack( fp, ModName ) ;
//    gen_rk4( fp, ModName ) ;

    if (strcmp(make_lower_temp(ModName->name), "airfoilinfo") == 0) { // make interpolation routines for AirfoilInfo module
        gen_ExtrapInterp(fp, ModName, "Output", "OutputType","ReKi");
        gen_ExtrapInterp(fp, ModName, "UA_BL_Type", "UA_BL_Type", "ReKi");
    } else if (!sw_noextrap) {
        if (strcmp(make_lower_temp(ModName->name), "dbemt") == 0) { // make interpolation routines for element-level DBEMT module
             
            gen_ExtrapInterp(fp, ModName, "ElementInputType", "ElementInputType", "DbKi");
        }

        gen_ExtrapInterp(fp, ModName, "Input", "InputType", "DbKi");
        gen_ExtrapInterp(fp, ModName, "Output", "OutputType", "DbKi");
    }

    fprintf(fp,"END MODULE %s_Types\n",ModName->name ) ;
  }

}


int
gen_module_files ( char * dirname, char * prog_ver )
{
  FILE * fp, *fph ;
  char  fname[NAMELEN], fname2[NAMELEN] ;

  node_t * p ;

  for ( p = ModNames ; p ; p = p->next )
  {
    if ( strlen( p->nickname ) > 0  && ! p->usefrom ) {
      fp = NULL ;

      if ( strlen(dirname) > 0 )
        { sprintf(fname,"%s/%s_Types.f90",dirname,p->name) ; }
      else
        { sprintf(fname,"%s_Types.f90",p->name) ; }
      sprintf(fname2, "%s_Types.f90", p->name);

      fprintf(stderr,"generating %s\n",fname) ;

      if ((fp = fopen( fname , "w" )) == NULL ) return(1) ;
      print_warning(fp,fname2, "");

      if ( sw_ccode == 1 ) {


        if ( strlen(dirname) > 0 )
          { sprintf(fname,"%s/%s_Types.h",dirname,p->name) ; }
        else
          { sprintf(fname, "%s_Types.h",p->name) ;}
        sprintf(fname2,"%s_Types.h",p->name) ;
        if ((fph = fopen( fname , "w" )) == NULL ) return(1) ;


        print_warning(fph,fname2, "//") ;

        fprintf(fph,"\n#ifndef _%s_TYPES_H\n",p->name);
        fprintf(fph,"#define _%s_TYPES_H\n\n",p->name);
        fprintf(fph,"\n#ifdef _WIN32 //define something for Windows (32-bit)\n");
        fprintf(fph,"#  include \"stdbool.h\"\n");
        fprintf(fph,"#  define CALL __declspec( dllexport )\n");
        fprintf(fph,"#elif _WIN64 //define something for Windows (64-bit)\n");
        fprintf(fph,"#  include \"stdbool.h\"\n");
        fprintf(fph,"#  define CALL __declspec( dllexport ) \n");
        fprintf(fph,"#else\n");
        fprintf(fph,"#  include <stdbool.h>\n");
        fprintf(fph,"#  define CALL \n");
        fprintf(fph,"#endif\n\n\n");
      }
      gen_module ( fp , p, prog_ver ) ;
      close_the_file( fp, "" ) ;
      if ( sw_ccode ) {
        gen_c_module ( fph , p ) ;

        fprintf(fph,"\n#endif // _%s_TYPES_H\n\n\n",p->name);
        close_the_file( fph,"//") ;

      }
    }
  }
  return(0) ;
}

void
remove_nickname( const char *nickname, char *src, char *dst )
{
  char tmp[NAMELEN];
  char srclo[NAMELEN];
  int n;
  strcpy(tmp,make_lower_temp(nickname)) ;
  strcpy(srclo, make_lower_temp(src));
  strcat(tmp,"_") ;
  n = strlen(tmp) ;
  if (!strncmp(tmp, srclo, n)) {
    strcpy(dst,&(src[n])) ;
  } else {
    strcpy(dst,src) ;
  }
}

void
append_nickname( const char *nickname, char *src, char *dst )
{
  int n ;
  n = strlen(nickname) ;
  if ( n > 0 ) {
    sprintf(dst,"%s_%s",nickname,src) ;
  } else {
    strcpy(dst,src) ;
  }
}

char * dimstr( int d )
{
  char * retval ;
  if        ( d == 0 ) {
    retval = "" ;
  } else if ( d == 1 ) {
    retval = "(i1)" ;
  } else if ( d == 2 ) {
    retval = "(i1,i2)" ;
  } else if ( d == 3 ) {
    retval = "(i1,i2,i3)" ;
  } else if ( d == 4 ) {
    retval = "(i1,i2,i3,i4)" ;
  } else if ( d == 5 ) {
    retval = "(i1,i2,i3,i4,i5)" ;
  } else {
    retval = " REGISTRY ERROR TOO MANY DIMS " ;
  }
  return(retval) ;

  //strcpy(dex, "");
  //strcat(dex, "(");
  //for (j = 1; j <= d; j++) {
  //   sprintf(tmp, "i%d%d", 0, j);
  //   strcat(dex, tmp);
  //   if (j == d) strcat(dex, ")"); else strcat(dex, ",");
  //}

}

char * dimstr_c( int d )
{
  char * retval ;
  if        ( d == 0 ) {
    retval = "" ;
  } else if ( d == 1 ) {
    retval = "[i1]" ;
  } else if ( d == 2 ) {
    retval = "[i2][i1]" ;
  } else if ( d == 3 ) {
    retval = "[i3][i2][i1]" ;
  } else if ( d == 4 ) {
    retval = "[i4][i3][i2][i1]" ;
  } else if ( d == 5 ) {
    retval = "[i5][i4][i3][i2][i1]" ;
  } else {
    retval = " REGISTRY ERROR TOO MANY DIMS " ;
  }
  return(retval) ;
}

void
checkOnlyReals( const char *q_mapsto, node_t * q) //, int recurselevel)
{
  node_t * r ;

  if ( q->type->type_type == DERIVED )
  {
     if ( strcmp( q->type->name, "meshtype" ) ) // skip meshes
     {
        for ( r = q->type->fields ; r ; r = r->next )
        {
           checkOnlyReals( q_mapsto, r);
        }
     }

  } else { // SIMPLE

     if ( strncmp(q->type->mapsto,"REAL",4) )
     {
         fprintf(stderr,"Registry warning: %s contains a field (%s) in a derived type whose type is not real or double: %s\n",
                q_mapsto, q->name , q->type->mapsto ) ;
     }

  }
  return;
}

void
checkContainsMesh( node_t * q) //, int recurselevel)
{
   node_t * r;

   if (q->type->type_type == DERIVED)
   {
      if (!strcmp(q->type->name, "meshtype") || !strcmp(q->type->name, "meshmaptype")){ // is a mesh or (a bad workaround for meshmaptype which contains meshtype in "usefrom" instead of "typedef")
         q->containsPtr = 1;
      }

      else {
         for (r = q->type->fields; r; r = r->next)
         {
            checkContainsMesh(r);
            if (r->containsPtr) q->containsPtr = 1;
         }
      }

   }

   return;
}
