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

#define C_KEY_LEN    256

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

#define PROPERTY_INT 1
#define PROPERTY_REAL 2
#define PROPERTY_STR 3
#define PROPERTY_LOGICAL 4

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

extern void (*dictionary_get_n)(int*, int*);
extern void (*dictionary_get_key)(int*, int*, char *, int*, int*, size_t);
extern void (*dictionary_query_key)(int*, char*, int*, int*, void*, int*, size_t);
extern void (*dictionary_query_index)(int*, int*, char*, int*, int*, void*, int*, size_t);
extern void (*dictionary_add_key)(int*, char*, int*, int*, void*, int*, size_t);
