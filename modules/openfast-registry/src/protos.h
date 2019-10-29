#ifndef PROTOS_H
#include "registry.h"
#include "data.h"

void substitute( char * str , char * match , char * replace, char * result );

int init_dim_table()   ;
char * make_lower( char * s1 ) ;
char * make_lower_temp( const char * s1 ) ;
int check_dimspecs();
int init_parser();
int is_a_fast_interface_type( char *str );
int pre_parse( char * dir, FILE * infile, FILE * outfile, int usefrom_sw );
int reg_parse( FILE * infile ) ;
int must_have_real_or_double( char *str );
int set_dim_len ( char * dimspec , node_t * dim_entry ) ;
int set_dim_order ( char * dimorder , node_t * dim_entry ) ;
int set_dim_orient ( char * dimorient , node_t * dim_entry ) ;
int add_node_to_end ( node_t * node , node_t ** list ) ;
int add_node_to_beg ( node_t * node , node_t ** list ) ;
int add_node_to_end_4d ( node_t * node , node_t ** list ) ;
int init_type_table() ;
int set_state_type ( char * , node_t *, node_t *, node_t * ) ;
int set_state_dims ( char * dims , node_t * node ) ;
int set_ctrl       ( char * ctrl , node_t * node ) ;
int gen_state_struct ( char * fname ) ;

#if 1
int show_node( node_t * p ) ;
int show_node1( node_t * p, int indent ) ;
void show_nodelist( node_t * p ) ;
void show_nodelist1( node_t * p , int indent ) ;
#endif

void gen_c_module( FILE * fph, node_t * ModName );

int gen_state_struct ( char * fname ) ;
int gen_decls ( FILE * fp ,  node_t * node , int sw_ranges, int sw_point , int mask , int layer ) ;
int gen_state_subtypes ( char * fname ) ;
int gen_state_subtypes1 ( FILE * fp , node_t * node , int sw_ranges, int sw_point , int mask ) ;
int print_warning( FILE * fp , char * fname, char * comment ) ;
void close_the_file( FILE * fp, char * comment  ) ;
int make_entries_uniq ( char * fname ) ;
int add_warning ( char * fname ) ;

int init_modname_table();
node_t * get_type_entry ( char * typename ) ;
node_t * get_modname_entry ( char * modname ) ;
node_t * get_rconfig_entry( char * name ) ;
node_t * get_entry ( char * name , node_t * node ) ;
node_t * get_entry_r ( char * name , char * use , node_t * node ) ;
node_t * get_dim_entry( char *s, int  ) ;
node_t * new_node ( int kind ) ;

node_t * get_4d_entry ( char * name ) ;
node_t * get_dimnode_for_coord ( node_t * node , int coord_axis ) ;
int      get_index_for_coord ( node_t * node , int coord_axis ) ;

char * my_strtok( char * s1 ) ;
char * strtok_rentr( char * s1 , char * s2, char ** tokpos ) ;

char * bdy_indicator( int bdy ) ;
char * make_upper_case ( char * str );
char * make_lower_case ( char * str );

char * field_name( char * tmp, node_t * p , int tag ) ;
char * field_name_bdy( char * tmp, node_t * p , int tag, int bdy  ) ;
char * dimension_with_colons( char * pre, char * tmp, node_t * p, char * post) ;
char * dimension_with_ones( char * pre, char * tmp, node_t * p, char * post) ;
char * dimension_with_ranges( char * ref , char * pre, int bdy , char * tmp, node_t * p, char * post, char * nlstructname ) ;
char * arrray_size_expression( char * refarg , char * pre , int bdy , char * tmp , node_t * p , char * post , char * nlstructname  ) ;
char * index_with_firstelem( char * pre , char * dref , int bdy , char * tmp , node_t * p , char * post ) ;

char * declare_array_as_pointer( char * tmp, node_t * p ) ;
char * field_type( char * tmp , node_t * p ) ;

/* For typedef history -ajb */
int init_typedef_history() ;
int add_typedef_name ( char * name ) ;
int get_num_typedefs() ;
char * get_typedef_name ( char * name ) ;
char * get_typename_i(int i) ;

int gen_alloc ( char * dirname ) ;
int gen_alloc1 ( char * dirname ) ;
int gen_alloc2 ( FILE * fp , char * structname , node_t * node, int *j, int *iguy, int *fraction, int numguys, int frac, int sw );

int gen_module_files ( char * dirname, char * prog_ver );
int gen_module_state_description ( char * dirname ) ;
int gen_module_state_description1 ( FILE * fp , node_t * node ) ;

void remove_nickname( const char *nickname, char *src, char *dst );
void append_nickname( const char *nickname, char *src, char *dst );
char * dimstr_c( int d );
void checkOnlyReals( const char *q_mapsto, node_t * q);
void checkContainsMesh(node_t * q);

int gen_scalar_indices ( char * dirname ) ;
int gen_scalar_indices1 ( FILE * fp, FILE ** fp2 ) ;

int gen_actual_args ( char * dirname ) ;
int gen_dummy_args ( char * dirname ) ;
int gen_dummy_decls ( char * dn ) ;
int gen_args ( char * dirname , int sw ) ;
int gen_args1 ( FILE * fp , char * outstr, char * structname , node_t * node , int *linelen , int sw , int deep ) ;

int gen_scalar_derefs ( char * dirname ) ;
int scalar_derefs ( char * dirname ) ;
int scalar_derefs1 ( FILE * fp , node_t * node, int direction ) ;

int set_mark ( int val , node_t * lst ) ;
int set_mark_4d ( int val , node_t * lst ) ;

int gen_i1_decls ( char * dn ) ;
int gen_get_nl_config ( char * dirname ) ;

int gen_config_assigns ( char * dirname ) ;
int gen_config_reads ( char * dirname ) ;

char * set_mem_order( node_t * node , char * str , int n  ) ;

int gen_wrf_io ( char * dirname ) ;
int set_dim_strs  ( node_t *node , char ddim[3][2][NAMELEN], char mdim[3][2][NAMELEN], char pdim[3][2][NAMELEN] , char * prepend, int sw_allow_stagger ) ;
int set_dim_strs2 ( node_t *node , char ddim[3][2][NAMELEN], char mdim[3][2][NAMELEN], char pdim[3][2][NAMELEN] , char * prepend, int sw_disregard_stag ) ;
int set_dim_strs3 ( node_t *node , char ddim[3][2][NAMELEN], char mdim[3][2][NAMELEN], char pdim[3][2][NAMELEN] , char * prepend, int sw_disregard_stag ) ;
int gen_wrf_io2 ( FILE * fp , char * fname , char * structname , char * fourdname , node_t * node , int sw_io ) ;

int gen_namelist_defines ( char * dirname , int sw_dimension ) ;
int gen_namelist_defaults ( char * dirname ) ;
int gen_namelist_script ( char * dirname ) ;

int gen_model_data_ord ( char * dirname ) ;

void get_elem ( char * structname , char * nlstructname , char * tx , int i , node_t * p , int first_last ) ;

int associated_with_4d_array( node_t * p ) ;

  
/* PGI Addition to resolve non-prototype function warnings  */
char * array_size_expression ( char *, char *, int, char *, node_t *, char * ,char * ); 
void range_of_dimension ( char *, char * , int, node_t *, char * );
void dimension_size_expression ( char *, char *, int, node_t *, char *);
int gen_alloc_count ( char *);
int gen_alloc_count1 ( char *);
int gen_ddt_write ( char * );
int gen_ddt_write1 ( FILE *, char *, node_t *);
int gen_dealloc ( char * );
int gen_dealloc1 ( char * );
int gen_dealloc2 ( FILE *, char *, node_t *);
int gen_scalar_tables ( FILE *);
int gen_scalar_tables_init ( FILE *);
int gen_scalar_indices_init ( FILE *);
int hash(char *);
int create_ht( char *** p ); 
int gen_nest_interp1 ( FILE *, node_t *, char *, int, int );
int gen_packs_halo ( FILE *fp , node_t *p, char *shw, int xy /* 0=y,1=x */ , int pu /* 0=pack,1=unpack */, char * packname, char * commname );
int gen_packs ( FILE *fp , node_t *p, int shw, int xy /* 0=y,1=x */ , int pu /* 0=pack,1=unpack */, char * packname, char * commname );
int gen_periods ( char * dirname , node_t * periods );
int gen_swaps ( char * dirname , node_t * swaps );
int gen_cycles ( char * dirname , node_t * cycles );
int gen_xposes ( char * dirname );
int gen_comm_descrips ( char * dirname );
int gen_shift (  char * dirname );
int gen_datacalls ( char * dirname );
int gen_nest_packing ( char * dirname );
int gen_nest_pack ( char * dirname );
int gen_nest_unpack ( char * dirname );
int gen_nest_packunpack ( FILE *fp , node_t * node , int dir, int down_path );
int count_fields ( node_t * node , int * d2 , int * d3 ,  char * fourd_names, int down_path );
int gen_debug (  char * dirname );

void reset_mask ( unsigned int * mask , int e ) ;
void set_mask ( unsigned int * mask , int e ) ;
int get_mask ( unsigned int * mask , int e ) ;

char * fast_interface_type_shortname (  char * ) ;
char * std_case( char * )  ;

char * dimstr( int ) ;

char * C_type ( char * ) ;
char * c_types_binding( char *s );
char * assoc_or_allocated( node_t  * r );
int is_pointer( node_t * r );
int has_deferred_dim( node_t * node, int noisy );

#define PROTOS_H
#endif

