/* H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */
/* H0 X                                                                            */
/* H0 X   libAtoms+QUIP: atomistic simulation library                              */
/* H0 X                                                                            */
/* H0 X   Portions of this code were written by                                    */
/* H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,      */
/* H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.          */
/* H0 X                                                                            */
/* H0 X   Copyright 2006-2010.                                                     */
/* H0 X                                                                            */
/* H0 X   These portions of the source code are released under the GNU General     */
/* H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html          */
/* H0 X                                                                            */
/* H0 X   If you would like to license the source code under different terms,      */
/* H0 X   please contact Gabor Csanyi, gabor@csanyi.net                            */
/* H0 X                                                                            */
/* H0 X   Portions of this code were written by Noam Bernstein as part of          */
/* H0 X   his employment for the U.S. Government, and are not subject              */
/* H0 X   to copyright in the USA.                                                 */
/* H0 X                                                                            */
/* H0 X                                                                            */
/* H0 X   When using this software, please cite the following reference:           */
/* H0 X                                                                            */
/* H0 X   http://www.libatoms.org                                                  */
/* H0 X                                                                            */
/* H0 X  Additional contributions by                                               */
/* H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras                 */
/* H0 X                                                                            */
/* H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */

#ifndef LIBATOMS_H
#define LIBATOMS_H

#define C_KEY_LEN    256

/* Dictionary types */

#define T_NONE       0
#define T_INTEGER    1
#define T_REAL       2
#define T_COMPLEX    3
#define T_LOGICAL    4
#define T_INTEGER_A  5
#define T_REAL_A     6
#define T_COMPLEX_A  7
#define T_LOGICAL_A  8
#define T_CHAR       9
#define T_CHAR_A     10
#define T_DATA       11
#define T_INTEGER_A2 12
#define T_REAL_A2    13

/* Property types -- legacy for backwards compatibility */

#define PROPERTY_INT 1
#define PROPERTY_REAL 2
#define PROPERTY_STR 3
#define PROPERTY_LOGICAL 4

/* Dictionary access macros */

#define INTEGER(loc) (*(int *)loc)
#define REAL(loc) (*(double *)loc)
#define LOGICAL(loc) (*(int *)loc)
#define COMPLEX_R(loc) (((double *)loc)[0])
#define COMPLEX_I(loc) (((double *)loc)[1])
#define CHAR(loc) ((char *)loc)
#define INTEGER_A(loc,i) (((int *)loc)[i])
#define REAL_A(loc,i) (((double *)loc)[i])
#define LOGICAL_A(loc,i) (((int *)loc)[i])
#define COMPLEX_A_R(loc,i)  (((double *)loc)[i*2+0])
#define COMPLEX_A_I(loc,i)  (((double *)loc)[i*2+1])
#define CHAR_A(loc, shape, i, j) (((char *)loc)[j*shape[0]+i])
#define INTEGER_A2(loc, shape, i, j) (((int *)loc)[j*shape[0]+i])
#define REAL_A2(loc, shape, i, j) (((double *)loc)[j*shape[0]+i])

/* error.f95 (via libAtoms_utils_no_module.f95) */

#ifdef CDEBUG
#define debug(fmt, ...) fprintf(stderr, fmt, ## __VA_ARGS__)
#else
#define debug(fmt, ...) 
#endif

#define ERROR_NONE                0
#define ERROR_UNSPECIFIED        -1
#define ERROR_IO                 -2
#define ERROR_IO_EOF             -3

#define INIT_ERROR if (error != NULL) *error = ERROR_NONE
#define RAISE_ERROR(info, ...) sprintf(error_h_info, info, ## __VA_ARGS__ ); error_h_line = __LINE__; error_h_kind = ERROR_UNSPECIFIED; c_push_error_with_info_(error_h_info, __FILE__, &error_h_line, &error_h_kind, strlen(error_h_info), strlen(__FILE__)); if (error != NULL) { *error = error_h_kind; return; } else c_error_abort_(error)
#define RAISE_ERROR_WITH_KIND(kind, info, ...) sprintf(error_h_info, info, ## __VA_ARGS__ ); error_h_line = __LINE__; error_h_kind = kind; c_push_error_with_info_(error_h_info, __FILE__, &error_h_line, &error_h_kind, strlen(error_h_info), strlen(__FILE__)); if (error != NULL) { *error = error_h_kind; return; } else c_error_abort_(error)
#define PASS_ERROR if (error != NULL && *error != ERROR_NONE) { error_h_line = __LINE__; c_push_error_(__FILE__, &error_h_line, error, strlen(__FILE__)); return; }
#define CLEAR_ERROR c_error_clear_stack_();

extern void c_push_error_with_info_(char*, char*, int*, int*, size_t, size_t);
extern void c_push_error_(char*, int*, int*, size_t);
extern void c_error_abort_(int *);
extern void c_error_clear_stack_(void);

extern char error_h_info[1000];
extern int error_h_line;
extern int error_h_kind;

/* quippy abort handler */

#include <setjmp.h>
extern jmp_buf environment_buffer;
extern char abort_message[1024];
void quippy_error_abort_(char *message, int len);
void quippy_error_abort_int_handler(int signum);

/* System.f95 (via libAtoms_utils_no_module.f95) */

extern void c_system_initialise_(int *verbosity);
#define system_initialise c_system_initialise_

/* Dictionary.f95 (via libAtoms_utils_no_module.f95) */

#ifndef SIZEOF_FORTRAN_T
#define SIZEOF_FORTRAN_T 12
#endif
typedef int fortran_t;

extern void c_dictionary_initialise_(fortran_t*);
#define dictionary_initialise c_dictionary_initialise_

extern void c_dictionary_finalise_(fortran_t*);
#define dictionary_finalise c_dictionary_finalise_

extern void c_dictionary_get_n_(fortran_t*, int*);
#define dictionary_get_n c_dictionary_get_n_

extern void c_dictionary_get_key_(fortran_t*, int*, char *, int*, int*, size_t);
#define dictionary_get_key c_dictionary_get_key_

extern void c_dictionary_query_key_(fortran_t*, char*, int*, int*, void*, int*, size_t);
#define dictionary_query_key c_dictionary_query_key_

extern void c_dictionary_query_index_(fortran_t*, int*, char*, int*, int*, void*, int*, size_t);
#define dictionary_query_index c_dictionary_query_index_

extern void c_dictionary_add_key_(fortran_t*, char*, int*, int*, void*, int*, size_t);
#define dictionary_add_key c_dictionary_add_key_

/* ExtendableStr.f95 (via libAtoms_utils_no_module.f95) */

extern void c_extendable_str_concat_(fortran_t*, char*, int *, int *, size_t);
#define extendable_str_concat c_extendable_str_concat_

/* Atoms.f95 (via libAtoms_utils_no_module.f95) */

void lattice_abc_to_xyz_(double cell_lengths[3], double cell_angles[3], double lattice[3][3]);
void lattice_xyz_to_abc_(double lattice[3][3], double cell_lengths[3], double cell_angles[3]);

#define LIBATOMS_DECLARE_CONFIG fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3], int *n_atom

/* xyz.c */

void read_xyz (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3], int *n_atom,
	       int compute_index, int frame, int *range, int string, int string_length, int *error);

void write_xyz (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3], int n_atom,
		int append, char *prefix, char *int_format, char *real_format, char *str_format, char *logical_format, 
		int string, fortran_t *estr, int update_index, int *error);

void query_xyz (char *filename, int compute_index, int frame, int *n_frame, int *n_atom, int *error);


/* netcdf.c */

void read_netcdf (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3], 
		  double cell_lengths[3], double cell_angles[3], int *cell_rotated,

		  int *n_atom, int frame, int zero, int *range, int irep, double rrep, int *error);
void write_netcdf (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3],
		   double cell_lengths[3], double cell_angles[3], int cell_rotated,
		   int n_atom, int n_label, int n_string, int frame, int netcdf4, int append,
		   int shuffle, int deflate, int deflate_level, int *error);
void query_netcdf (char *filename, int *n_frame, int *n_atom, int *n_label, int *n_string, int *error);


#endif
