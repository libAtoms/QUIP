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

/* C interface to Fortran error handling routines */

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
#define RAISE_ERROR(info, ...) sprintf(error_h_info, info, ## __VA_ARGS__ ); error_h_line = __LINE__; error_h_kind = ERROR_UNSPECIFIED; push_error_with_info(error_h_info, __FILE__, &error_h_line, &error_h_kind, strlen(error_h_info), strlen(__FILE__)); if (error != NULL) { *error = -1; return; } else error_abort(error)
#define PASS_ERROR if (error != NULL && *error != 0) { error_h_line = __LINE__; error_h_kind = ERROR_UNSPECIFIED; push_error(__FILE__, &error_h_line, &error_h_kind, strlen(__FILE__)); return; }
#define CLEAR_ERROR error_clear_stack();

/* Function pointers to Fortran push_error() and push_error_with_info() in error.f95 */
extern void (*push_error_with_info)(char*, char*, int*, int*, size_t, size_t);
extern void (*push_error)(char*, int*, int*, size_t);
extern void (*error_abort)(int *);
extern void (*error_clear_stack)();

extern char error_h_info[1000];
extern int error_h_line;
extern int error_h_kind;
