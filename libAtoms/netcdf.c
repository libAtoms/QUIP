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

/* Read and write Atoms objects in NetCDF 3 or 4 format */

#ifdef HAVE_NETCDF
#include <netcdf.h>
#endif

#include <stdio.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <float.h>

#include "libatoms.h"

#define absval(x)  ( (x) < 0 ? -(x) : (x) )
#define NETCDF_CHECK(s) if ((retval = (s))) { RAISE_ERROR_WITH_KIND(ERROR_IO, "%s",nc_strerror(retval)); }
#define DEG_TO_RAD (M_PI/180.0)
#define RAD_TO_DEG (180.0/M_PI)

// Utility functions

#ifdef HAVE_NETCDF
 
void replace_fill_values(fortran_t *params, fortran_t *properties, int irep, double rrep, int *error) {
  int i, j=0, k=0, n, d, type, shape[2];
  char key[C_KEY_LEN];
  void *loc;
  fortran_t *dictionaries[2];

  INIT_ERROR;
  dictionaries[0] = params;
  dictionaries[1] = properties;

  for (d=0; d<2; d++) {
    if (dictionaries[d] == NULL) continue;

    dictionary_get_n(dictionaries[d], &n);
    for (i=1; i <= n; i++) {
      dictionary_query_index(dictionaries[d], &i, key, &type, shape, &loc, error, C_KEY_LEN);
      PASS_ERROR;
      if (loc == NULL) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "replace_fill_values: NULL pointer for entry i=%d", i);
      }

      switch(type) {
      case(T_INTEGER):
      case(T_LOGICAL):
	if (INTEGER(loc) == NC_FILL_INT)
	  INTEGER(loc) = irep;
	break;

      case(T_REAL):
	if((REAL(loc) > 0) == (NC_FILL_DOUBLE > 0) && /* prevents potential overflow */
	 (absval(REAL(loc) - NC_FILL_DOUBLE) <= absval(DBL_EPSILON * NC_FILL_DOUBLE)))
	REAL(loc) = rrep;
      break;

      case(T_INTEGER_A):
      case(T_LOGICAL_A):
	for (j=0; j<shape[0]; j++)
	  if (INTEGER_A(loc,j) == NC_FILL_INT)
	    INTEGER_A(loc,j) = irep;
	break;

      case(T_REAL_A):
	for (j=0; j<shape[0]; j++)
	  if((REAL_A(loc,j) > 0) == (NC_FILL_DOUBLE > 0) && /* prevents potential overflow */
	     (absval(REAL_A(loc,j) - NC_FILL_DOUBLE) <= absval(DBL_EPSILON * NC_FILL_DOUBLE)))
	    REAL_A(loc,j) = rrep;
	break;

      case(T_INTEGER_A2):
	for (k=0; k<shape[1]; k++)
	  for (j=0; j<shape[0]; j++)
	    if (INTEGER_A2(loc,shape,j,k) == NC_FILL_INT)
	      INTEGER_A2(loc,shape,j,k) = irep;
	break;

      case(T_REAL_A2):
	for (k=0; k<shape[1]; k++)
	  for (j=0; j<shape[0]; j++)
	    if((REAL_A2(loc,shape,j,k) > 0) == (NC_FILL_DOUBLE > 0) && /* prevents potential overflow */
	       (absval(REAL_A2(loc,shape,j,k) - NC_FILL_DOUBLE) <= absval(DBL_EPSILON * NC_FILL_DOUBLE)))
	      REAL_A2(loc,shape,j,k) = rrep;
	break;
      }
    }
  }
}

void convert_from_netcdf_type(char *varname, nc_type vartype, int ndims, int *dimids, 
			      int frame_dim_id, int spatial_dim_id, int atom_dim_id, 
			      int label_dim_id, int string_dim_id, int spatial2_dim_id,
			      int n_spatial, int n_atom, int n_label, int n_string, int n_spatial2,
			      int *type, int *is_param, int *is_property, int *is_global, int *shape, int *error)
{
  int type_map[3][13];
  int i, j;

  INIT_ERROR;

  // Set up type_map. First index is ndims-1, second is vartype
  for (j=0; j<3; j++) 
    for (i=0; i<13; i++)
      type_map[j][i] = -1;

  type_map[0][NC_INT] = T_INTEGER;
  type_map[0][NC_FLOAT] = T_REAL;
  type_map[0][NC_DOUBLE] = T_REAL;

  type_map[1][NC_INT] = T_INTEGER_A;
  type_map[1][NC_FLOAT] = T_REAL_A;
  type_map[1][NC_DOUBLE] = T_REAL_A;
  type_map[1][NC_CHAR] = T_CHAR;

  type_map[2][NC_INT] = T_INTEGER_A2;
  type_map[2][NC_FLOAT] = T_REAL_A2;
  type_map[2][NC_DOUBLE] = T_REAL_A2;
  type_map[2][NC_CHAR] = T_CHAR_A;

  *is_param = 0;
  *is_property = 0;
  *is_global = 0;
  shape[0] = 0;
  shape[1] = 0;

  if (ndims == 0) 
    ; // one off variable, do nothing
  else if (ndims == 1) {
    if (dimids[0] == frame_dim_id) {
      // scalar per-frame parameter
      *is_param = 1;
      *type = type_map[ndims-1][vartype];
      if (*type == -1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: bad vartype %d for varname %s", vartype, varname);
      }
    }
    else if (dimids[0] == atom_dim_id) {
      // vector per-file parameter, needs to be read per frame
      *is_property = 1;
      *is_global = 1;
      *type = type_map[ndims][vartype];
      if (*type == -1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: bad vartype %d for varname %s", vartype, varname);
      }
      shape[0] = n_atom;
    }
    else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: unknown one-dimensional variable %s dimid=%d\n", varname, dimids[0]);
    }
  }
  else if (ndims == 2) {
    if (dimids[0] == frame_dim_id && dimids[1] == atom_dim_id) {
      // scalar per atom property
      *is_property = 1;
      *type = type_map[ndims-1][vartype];
      if (*type == -1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: bad vartype %d for varname %s", vartype, varname);
      }
      shape[0] = n_atom;
    }
    else if (dimids[0] == frame_dim_id && dimids[1] == spatial_dim_id) {
      // vector per-frame parameter
      *is_param = 1;
      *type = type_map[ndims-1][vartype];
      if (*type == -1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: bad vartype %d for varname %s", vartype, varname);
      }
      shape[0] = n_spatial;
    }
    else if (spatial2_dim_id != 0 && dimids[0] == frame_dim_id && dimids[1] == spatial2_dim_id) {
      // tensor per-frame parameter
      *is_param = 1;
      *type = type_map[ndims-1][vartype];
      if (*type == -1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: bad vartype %d for varname %s", vartype, varname);
      }
      shape[0] = n_spatial2;
    } else if (dimids[0] == frame_dim_id && dimids[1] == string_dim_id) {
      // string per-frame parameter
      *is_param = 1;
      *type = T_CHAR;
      shape[0] = n_string;
    }
    else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: unknown two dimensional variable %s\n", varname);
    }
  }
  else if (ndims == 3) {
    if (dimids[0] == frame_dim_id && dimids[1] == atom_dim_id) {
	  
      if (dimids[2] == label_dim_id) {
	// per atom string property
	*is_property = 1;
	*type = T_CHAR_A;
	shape[0] = n_label;
	shape[1] = n_atom;
      }
      else if (dimids[2] == spatial_dim_id) {
	// vector per atom property
	*is_property = 1;
	*type = type_map[ndims-1][vartype];
	if (*type == -1) {
	  RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: bad vartype %d for varname %s", vartype, varname);
	}
	shape[0] = n_spatial;
	shape[1] = n_atom;
      }
      else if (spatial2_dim_id != 0 && dimids[2] == spatial2_dim_id) {
	// tensor per atom property
	*is_property = 1;
	*type = type_map[ndims-1][vartype];
	if (*type == -1) {
	  RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: bad vartype %d for varname %s", vartype, varname);
	}
	shape[0] = n_spatial2;
	shape[1] = n_atom;
      }
    } else if (dimids[0] == frame_dim_id && dimids[1] == spatial_dim_id && dimids[2] == spatial_dim_id) {
      // per-frame integer or real 2D parameter array 
      *is_param = 1;
      *type = type_map[ndims-1][vartype];
      if (*type == -1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: bad vartype %d for varname %s", vartype, varname);
      }      
      shape[0] = n_spatial;
      shape[1] = n_spatial;
    }
    else {
      debug("spatial_dim_id = %d, frame_dim_id = %d\n", spatial_dim_id, frame_dim_id);
      RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: Unknown three dimensional variable %s dimids=[%d %d %d]", varname, dimids[0], dimids[1], dimids[2]);
    }
  }
  else {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_from_netcdf_type: Unknown variable %s, dimension %d", varname, ndims);
  }
}

void convert_to_netcdf_type(char *key, int type, int *shape, 
			    int frame_dim_id, int spatial_dim_id, int atom_dim_id, 
			    int label_dim_id, int string_dim_id, int spatial2_dim_id,
			    int n_spatial, int n_atom, int n_label, int n_string, int n_spatial2,
			    nc_type *vartype, int *ndims, int *dims, int *error)
{
  INIT_ERROR;

  *vartype = 0;
  *ndims = 0;
  dims[0] = 0;
  dims[1] = 0;
  dims[2] = 0;
  
  switch (type) {
  case(T_INTEGER):
  case(T_LOGICAL):
  case(T_INTEGER_A):
  case(T_LOGICAL_A):
  case(T_INTEGER_A2):
    *vartype = NC_INT;
    break;

  case(T_REAL):
  case(T_REAL_A):
  case(T_REAL_A2):
    *vartype = NC_DOUBLE;
    break;

  case(T_CHAR):
  case(T_CHAR_A):
    *vartype = NC_CHAR;
    break;

  default:
    RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_to_netcdf_type: bad type %d for key %s", type, key);
  }

  dims[0] = frame_dim_id;
  if (type == T_INTEGER || type == T_LOGICAL || type == T_REAL) {
    *ndims = 1;
  } else if (type == T_CHAR) {
    *ndims = 2;
    dims[1] = string_dim_id;
  } else if (type == T_INTEGER_A || type == T_LOGICAL_A || type == T_REAL_A) {
    *ndims = 2;
    if (shape[0] == n_spatial) {
      dims[1] = spatial_dim_id;
    } else if (shape[0] == n_atom) {
      dims[1] = atom_dim_id;
    } else if (spatial2_dim_id != 0 && shape[0] == n_spatial2) {
      dims[1] = spatial2_dim_id;
    } else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_to_netcdf_type: bad array length %d for key %s", shape[0], key);
    }
  } else if (type == T_INTEGER_A2 || type == T_REAL_A2) {
    *ndims = 3;
    if (shape[0] == n_spatial) {
      dims[2] = spatial_dim_id;
    } else if (spatial2_dim_id != 0 && shape[0] == n_spatial2) {
      dims[2] = spatial2_dim_id;
    } else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_to_netcdf_type: bad array length %d != n_spatial or n_spatial2 for key %s", shape[0], key);
    }
    if (shape[1] == n_spatial) {
      dims[1] = spatial_dim_id;
    } else if (shape[1] == n_atom) {
      dims[1] = atom_dim_id;
    } else if (spatial2_dim_id != 0 && shape[1] == n_spatial2) {
      dims[1] = spatial2_dim_id;
    } else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_to_netcdf_type: bad array length %d != n_spatial or n_spatial2 or n_atom for key %s", shape[1], key);
    }
  } else if (type == T_CHAR_A) {
    *ndims = 3;
    if (shape[0] == n_label) 
      dims[2] = label_dim_id;
    else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_to_netcdf_type: bad array length %d != n_label for key %s", shape[1], key);
    }
    if (shape[1] == n_atom)
      dims[1] = atom_dim_id;
    else {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "convert_to_netcdf_type: bad array length %d != n_atom for key %s", shape[2], key);      
    }
  }
}


void read_netcdf (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3], 
		  double cell_lengths[3], double cell_angles[3], int *cell_rotated,
		  int *n_atom, int frame, int zero, int *range, int irep, double rrep, int *error)
{
  int nc_id, i, n_selected;
  int retval, nvars, at_start, at_end;
  size_t start[3], count[3];
  char varname[NC_MAX_NAME+1];
  nc_type vartype, type_type;
  int ndims, dimids[NC_MAX_VAR_DIMS], natts, shape[2], type, type_att;
  size_t tmp_sizet;
  void *data, *tmp_data;
  int is_param, is_property, is_global;
  int frame_dim_id, spatial_dim_id, atom_dim_id, cell_spatial_dim_id,
    cell_angular_dim_id, label_dim_id, string_dim_id, spatial2_dim_id;
  int n_spatial, n_frame, n_label, n_string, n_spatial2;
  int tmp_type, tmp_shape[2], tmp_error;
  INIT_ERROR;

  // Open file
  NETCDF_CHECK(nc_open(filename, NC_NOWRITE, &nc_id));
  
  // Inquire dimensions
  NETCDF_CHECK(nc_inq_dimid(nc_id, "frame", &frame_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "spatial", &spatial_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "atom", &atom_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "cell_spatial", &cell_spatial_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "cell_angular", &cell_angular_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "label", &label_dim_id));
  if (nc_inq_dimid(nc_id, "string", &string_dim_id) != NC_NOERR) {
    // No strings in this file
    string_dim_id = 0;
  }
  if (nc_inq_dimid(nc_id, "spatial2", &spatial2_dim_id) != NC_NOERR) {
    spatial2_dim_id = 0;
  }

  // Get sizes of dimensions
  NETCDF_CHECK(nc_inq_dimlen(nc_id, spatial_dim_id, &tmp_sizet));
  n_spatial = (int)tmp_sizet;
  if (n_spatial != 3) {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "read_netcdf: number of spatial dimensions = %d != 3", n_spatial);
  }
  if (spatial2_dim_id != 0) {
    NETCDF_CHECK(nc_inq_dimlen(nc_id, spatial2_dim_id, &tmp_sizet));
    n_spatial2 = (int)tmp_sizet;
    if (n_spatial2 != 9) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "read_netcdf: number of spatial dimensions**2 = %d != 9", n_spatial2);
    }
  }
  NETCDF_CHECK(nc_inq_dimlen(nc_id, frame_dim_id, &tmp_sizet));
  n_frame = (int)tmp_sizet;
  NETCDF_CHECK(nc_inq_dimlen(nc_id, atom_dim_id, &tmp_sizet));
  *n_atom = (int)tmp_sizet;
  NETCDF_CHECK(nc_inq_dimlen(nc_id, label_dim_id, &tmp_sizet));
  n_label = (int)tmp_sizet;
  if (string_dim_id != 0) {
    NETCDF_CHECK(nc_inq_dimlen(nc_id, string_dim_id, &tmp_sizet));
    n_string = (int)tmp_sizet;
  }
  else {
    n_string = 0;
  }

  debug("read_netcdf: got %d frames, %d atoms\n", n_frame, *n_atom);

  if (frame < 0 || frame >= n_frame) {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "read_netcdf: frame %d out of range 0 <= frame < %d", frame, n_frame);
  }

  // Have we been asked to read only a specific range of atom indices?
  if (range[0] != 0 && range[1] != 0) {
    if (range[0] == -1 && range[1] == -1) {
      // special range  of [-1, -1] means don't read any atoms, only params and lattice
      *n_atom = 0;
      at_start = 0;
      at_end = -1;
    } else {
      if (range[0] < 1) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "read_netcdf: lower limit of range (%d) must be >= 1", range[0]);
      }
      if (range[1] > *n_atom) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "read_netcdf: upper limit of range (%d) must be <= %d", range[1], *n_atom);
      }
      if (range[1] <= range[0]) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "read_netcdf: upper limit of range (%d) must be > lower limit (%d)", range[1], range[0]);
      }
      *n_atom = range[1] - range[0] + 1;
      at_start = range[0]-1;
      at_end = range[1]-1;
    }
  }
  else {
    at_start = 0;
    at_end = *n_atom-1;
  }

  debug("read_netcdf: range=[%d %d]\n", range[0], range[1]);
  debug("read_netcdf: at_start=%d at_end=%d n_atom=%d\n", at_start, at_end, *n_atom);

  n_selected = 0;
  if (selected_properties != NULL) dictionary_get_n(selected_properties, &n_selected);

  // By default, assume cell is not rotated
  *cell_rotated = 0;
  
  // Loop over all variables in file
  NETCDF_CHECK(nc_inq_nvars(nc_id, &nvars));
  for (i=0; i<nvars; i++) {
    NETCDF_CHECK(nc_inq_var(nc_id, i, varname, &vartype, &ndims, dimids, &natts));
    
    if (strcmp(varname,"cell_lengths") == 0) {
      start[0] = frame;   start[1] = 0;
      count[0] = 1;       count[1] = 3;
      NETCDF_CHECK(nc_get_vara_double(nc_id, i, start, count, cell_lengths));
      continue;
    }
    
    if (strcmp(varname,"cell_angles") == 0) {
      start[0] = frame;   start[1] = 0;
      count[0] = 1;       count[1] = 3;
      NETCDF_CHECK(nc_get_vara_double(nc_id, i, start, count, cell_angles));
      continue;
    }

    if (strcmp(varname, "cell_rotated") == 0) {
      start[0] = frame;
      count[0] = 1;
      NETCDF_CHECK(nc_get_vara_int(nc_id, i, start, count, cell_rotated));
      continue;
    }

    if (strcmp(varname, "cell_lattice") == 0) {
      start[0] = frame;   start[1] = 0;   start[2] = 0;
      count[0] = 1;       count[1] = 3;   count[2] = 3;
      NETCDF_CHECK(nc_get_vara_double(nc_id, i, start, count, &(lattice[0][0])));
      continue;
    }
    
    if (strcmp(varname, "spatial") == 0 || strcmp(varname, "cell_spatial") == 0 ||
	strcmp(varname, "cell_angular") == 0)
      continue;

    // Translate names: coordinates -> pos and velocities -> velo
    if (strcmp(varname,"coordinates") == 0)
      strncpy(varname,"pos",NC_MAX_NAME);
    if (strcmp(varname,"velocities") == 0)
      strncpy(varname,"velo",NC_MAX_NAME);    
    convert_from_netcdf_type(varname, vartype, ndims, dimids, 
			     frame_dim_id, spatial_dim_id, atom_dim_id, label_dim_id, string_dim_id, spatial2_dim_id,
			     n_spatial, *n_atom, n_label, n_string, n_spatial2,
			     &type, &is_param, &is_property, &is_global, shape, error);
    PASS_ERROR;

    debug("read_netcdf: got variable \"%s\" type %d shape [%d %d] is_param=%d is_property=%d\n", varname, type, shape[0], shape[1], is_param, is_property);

    if (!is_param && !is_property) continue;

    // Are we filtering the properties we read?
    if (is_property && n_selected != 0) {
      dictionary_query_key(selected_properties, varname, &tmp_type, tmp_shape, &tmp_data, &tmp_error, strlen(varname));
      CLEAR_ERROR;
      if (tmp_error != ERROR_NONE) continue;
    }
     
    // Is there a "type" attribute for this variable? If so override type with it, 
    // for backwards compatibility. This is used to distinguish between Fortran integers
    // and logicals.
    if (nc_inq_att(nc_id, i, "type", &type_type, &tmp_sizet) == NC_NOERR && type_type == NC_INT) {
      nc_get_att_int(nc_id, i, "type", &type_att);
      if (is_property) {
	if (type_att < 1 || type_att > 4) {
	  RAISE_ERROR_WITH_KIND(ERROR_IO, "read_netcdf: type attribtue %d out of range 1 <= type_att <= 4", type_att);
	}
	switch (type_att) {
	case(PROPERTY_INT):
	  if (ndims == 2)
	    type = T_INTEGER_A;
	  else
	    type = T_INTEGER_A2;
	  break;
	case(PROPERTY_LOGICAL):
	  if (ndims == 2)
	    type = T_LOGICAL_A;
	  else {
	    RAISE_ERROR_WITH_KIND(ERROR_IO, "read_netcdf: logical property with ndims=%d no longer supported", ndims);
	  }
	  break;
	}
      }
      else
	type = type_att;
    }

    if (is_global) {
      start[0] = 0;
    }
    else {
      start[0] = frame;
    }
    start[1] = 0;
    start[2] = 0;

    if (is_property) {
      start[1] = at_start;
    }

    count[0] = 1;
    // If it's a 2D array we need to transpose shape
    if (type == T_INTEGER_A || type == T_REAL_A || type == T_LOGICAL_A || type == T_CHAR) {
      if (is_global) {
	count[0] = shape[0];
      }
      else {
	count[1] = shape[0];
      }
    } else {
      count[1] = shape[1];
      count[2] = shape[0];
    }

    // Create destination arrays in relevant Dictionary: params or properties
    if (is_param && params != NULL) {
      dictionary_add_key(params, varname, &type, shape, &data, error, strlen(varname));
      PASS_ERROR;
    }
    else if (is_property && properties != NULL) {
      dictionary_add_key(properties, varname, &type, shape, &data, error, strlen(varname));
      PASS_ERROR;
    } 

    // Do the reading.
    debug("read_netcdf: reading variable \"%s\" start=[%ld %ld %ld] count=[%ld %ld %ld]\n", varname, start[0], start[1], start[2], count[0], count[1], count[2]);
    switch(type) {
    case(T_INTEGER):
    case (T_LOGICAL):
      NETCDF_CHECK(nc_get_var1_int(nc_id, i, start, (int *)data));
      break;
    case(T_REAL):
      NETCDF_CHECK(nc_get_var1_double(nc_id, i, start, (double *)data));
      break;
    case(T_CHAR):
      memset(data, ' ', n_string); // pad with spaces for fortran
      NETCDF_CHECK(nc_get_vara_text(nc_id, i, start, count, (char *)data));
      if (strnlen(data, n_string) != n_string) { // overwrite '\0' with ' ' to fix _FillValue bugs
	memset(data+strnlen(data, n_string), ' ', n_string - strnlen(data, n_string));
      }
      break;
    case(T_INTEGER_A):
    case(T_LOGICAL_A):
      NETCDF_CHECK(nc_get_vara_int(nc_id, i, start, count, (int *)data));
      break;
    case(T_REAL_A):
      NETCDF_CHECK(nc_get_vara_double(nc_id, i, start, count, (double *)data));
      break;
    case(T_INTEGER_A2):
      NETCDF_CHECK(nc_get_vara_int(nc_id, i, start, count, (int *)data));
      break;
    case(T_REAL_A2):
      NETCDF_CHECK(nc_get_vara_double(nc_id, i, start, count, (double *)data));
      break;      
    case(T_CHAR_A):
      memset(data, ' ', (*n_atom)*n_label); // FIXME may need to change padding from '\0' to ' '
      NETCDF_CHECK(nc_get_vara_text(nc_id, i, start, count, (char *)data));
      break;      
    default:
      RAISE_ERROR_WITH_KIND(ERROR_IO, "read_netcdf: unknown variable \"%s\" type %d\n", varname, type);
    }
  }
  
  debug("read_netcdf: cell_lengths = [%f %f %f]\n", cell_lengths[0], cell_lengths[1], cell_lengths[2]);
  debug("read_netcdf: cell_angles = [%f %f %f]\n", cell_angles[0], cell_angles[1], cell_angles[2]);
  for (i=0; i<3; i++)
    cell_angles[i] *= DEG_TO_RAD;

  if (zero) {
    replace_fill_values(params, properties, irep, rrep, error);
    PASS_ERROR;
  }

  nc_close(nc_id);
} 

void write_netcdf (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3],
		   double cell_lengths[3], double cell_angles[3], int cell_rotated, 
		   int n_atom, int n_label, int n_string, int frame, int netcdf4, int append,
		   int shuffle, int deflate, int deflate_level, int *error)
{
  int nc_id, retval;
  int d, i, j;
  size_t start[3], count[3], tmp_sizet;
  int dims[NC_MAX_VAR_DIMS], check_dims[NC_MAX_VAR_DIMS];
  int ndims, check_ndims;
  char varname[NC_MAX_NAME+1];
  nc_type vartype, check_vartype;
  int natts, newfile;
  int nvars, ngatts, unlimdimid, newvar;
  int frame_dim_id, spatial_dim_id, atom_dim_id, cell_spatial_dim_id,
    cell_angular_dim_id, label_dim_id, string_dim_id, spatial2_dim_id;
  int spatial_var_id, cell_spatial_var_id, cell_angular_var_id, cell_lengths_var_id, cell_angles_var_id, 
    cell_rotated_var_id, cell_lattice_var_id;
  int n_spatial, n_frame, n_spatial2;
  int n, type, type_att, shape[2], var_id, tmp_type, tmp_shape[2], tmp_error;
  char key[C_KEY_LEN];
  void *data, *tmp_data;
  fortran_t *dictionaries[2];
  int add_cell_rotated, add_cell_lattice, rotate;
  double new_lattice[3][3];
  char fill[2] = " ";

  INIT_ERROR;
  dictionaries[0] = params;
  dictionaries[1] = properties;

  // Open the NetCDF file for writing
  if (append) {
    debug("write_netcdf: opening netcdf file %s for append\n", filename);
#ifdef NETCDF4
    NETCDF_CHECK(nc_open(filename, NC_WRITE, &nc_id));
#else
    NETCDF_CHECK(nc_open(filename, NC_64BIT_OFFSET | NC_WRITE, &nc_id));
#endif
  } else {
    debug("write_netcdf: creating netcdf file %s\n", filename);
#ifdef NETCDF4
    if (netcdf4) {
      NETCDF_CHECK(nc_set_default_format(NC_FORMAT_NETCDF4, NULL));
      NETCDF_CHECK(nc_create(filename, NC_NETCDF4 | NC_CLOBBER, &nc_id));
    } else {
      NETCDF_CHECK(nc_set_default_format(NC_FORMAT_64BIT, NULL));
      NETCDF_CHECK(nc_create(filename, NC_64BIT_OFFSET | NC_CLOBBER, &nc_id));
    }
#else
    NETCDF_CHECK(nc_create(filename, NC_64BIT_OFFSET | NC_CLOBBER, &nc_id));
#endif
  }
    
  NETCDF_CHECK(nc_inq(nc_id, &ndims, &nvars, &ngatts, &unlimdimid));
  newfile = ndims == 0;

  if (newfile) { // it's a new file
    // set dimension and variable information

    // Global attributes
    NETCDF_CHECK(nc_put_att_text(nc_id, NC_GLOBAL, "Conventions", strlen("AMBER"), "AMBER"));
    NETCDF_CHECK(nc_put_att_text(nc_id, NC_GLOBAL, "ConventionVersion", strlen("1.0"), "1.0"));
    NETCDF_CHECK(nc_put_att_text(nc_id, NC_GLOBAL, "application", strlen("libAtoms"), "libAtoms"));
    NETCDF_CHECK(nc_put_att_text(nc_id, NC_GLOBAL, "program", strlen("netcdf.c"), "netcdf.c"));
    NETCDF_CHECK(nc_put_att_text(nc_id, NC_GLOBAL, "programVersion", strlen(GIT_VERSION), GIT_VERSION));
    NETCDF_CHECK(nc_put_att_text(nc_id, NC_GLOBAL, "title", strlen("Atoms Object"), "Atoms Object"));
    
    // Dimensions
    n_spatial = 3;
    n_spatial2 = 9;
    NETCDF_CHECK(nc_def_dim(nc_id, "frame", NC_UNLIMITED, &frame_dim_id));
    NETCDF_CHECK(nc_def_dim(nc_id, "spatial", n_spatial, &spatial_dim_id));
    NETCDF_CHECK(nc_def_dim(nc_id, "atom", n_atom, &atom_dim_id));
    NETCDF_CHECK(nc_def_dim(nc_id, "cell_spatial", n_spatial, &cell_spatial_dim_id));
    NETCDF_CHECK(nc_def_dim(nc_id, "cell_angular", n_spatial, &cell_angular_dim_id));
    NETCDF_CHECK(nc_def_dim(nc_id, "label", n_label, &label_dim_id));
    NETCDF_CHECK(nc_def_dim(nc_id, "string", n_string, &string_dim_id));
    NETCDF_CHECK(nc_def_dim(nc_id, "spatial2", n_spatial2, &spatial2_dim_id));

    // Label variables
    dims[0] = spatial_dim_id;
    NETCDF_CHECK(nc_def_var(nc_id, "spatial", NC_CHAR, 1, dims, &spatial_var_id));

    dims[0] = cell_spatial_dim_id;
    NETCDF_CHECK(nc_def_var(nc_id, "cell_spatial", NC_CHAR, 1, dims, &cell_spatial_var_id));
      
    dims[0] = cell_angular_dim_id;
    dims[1] = label_dim_id;
    NETCDF_CHECK(nc_def_var(nc_id, "cell_angular", NC_CHAR, 2, dims, &cell_angular_var_id));
    
    dims[0] = frame_dim_id;
    dims[1] = spatial_dim_id;
    NETCDF_CHECK(nc_def_var(nc_id, "cell_lengths", NC_DOUBLE, 2, dims, &cell_lengths_var_id));
    NETCDF_CHECK(nc_def_var(nc_id, "cell_angles", NC_DOUBLE, 2, dims, &cell_angles_var_id));

    add_cell_rotated = 1;
    add_cell_lattice = 1;
    
#ifdef NETCDF4
    if (netcdf4) {
      NETCDF_CHECK(nc_def_var_deflate(nc_id, spatial_var_id, shuffle, deflate, deflate_level));
      NETCDF_CHECK(nc_def_var_deflate(nc_id, cell_spatial_var_id, shuffle, deflate, deflate_level));
      NETCDF_CHECK(nc_def_var_deflate(nc_id, cell_angular_var_id, shuffle, deflate, deflate_level));
      NETCDF_CHECK(nc_def_var_deflate(nc_id, cell_lengths_var_id, shuffle, deflate, deflate_level));
      NETCDF_CHECK(nc_def_var_deflate(nc_id, cell_angles_var_id, shuffle, deflate, deflate_level));
    }
#endif
  } else {
     // Inquire dimensions
    NETCDF_CHECK(nc_inq_dimid(nc_id, "frame", &frame_dim_id));
    NETCDF_CHECK(nc_inq_dimid(nc_id, "spatial", &spatial_dim_id));
    NETCDF_CHECK(nc_inq_dimid(nc_id, "atom", &atom_dim_id));
    NETCDF_CHECK(nc_inq_dimid(nc_id, "cell_spatial", &cell_spatial_dim_id));
    NETCDF_CHECK(nc_inq_dimid(nc_id, "cell_angular", &cell_angular_dim_id));
    NETCDF_CHECK(nc_inq_dimid(nc_id, "label", &label_dim_id));
    NETCDF_CHECK(nc_inq_dimid(nc_id, "string", &string_dim_id));
    NETCDF_CHECK(nc_inq_dimid(nc_id, "spatial2", &spatial2_dim_id));
    
    // Check sizes of dimensions
    NETCDF_CHECK(nc_inq_dimlen(nc_id, spatial_dim_id, &tmp_sizet));
    n_spatial = (int)tmp_sizet;
    if (tmp_sizet != 3) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: number of spatial dimensions = %d != 3", n_spatial);
    }

    NETCDF_CHECK(nc_inq_dimlen(nc_id, frame_dim_id, &tmp_sizet));
    n_frame = (int)tmp_sizet;

    NETCDF_CHECK(nc_inq_dimlen(nc_id, atom_dim_id, &tmp_sizet));
    if (tmp_sizet != n_atom) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: number of atoms defined in file = %d != %d", (int)tmp_sizet, n_atom);
    }

    NETCDF_CHECK(nc_inq_dimlen(nc_id, label_dim_id, &tmp_sizet));
    if (tmp_sizet != n_label) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: n_label defined in file = %d != %d", (int)tmp_sizet, n_label);
    }

    NETCDF_CHECK(nc_inq_dimlen(nc_id, string_dim_id, &tmp_sizet));
    if (tmp_sizet != n_string) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: n_string defined in file = %d != %d", (int)tmp_sizet, n_string);
    }

    if (spatial2_dim_id != 0) {
      NETCDF_CHECK(nc_inq_dimlen(nc_id, spatial2_dim_id, &tmp_sizet));
      n_spatial2 = (int)tmp_sizet;
      if (tmp_sizet != 9) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: n_spatial2 defined in file = %d != %d", (int)tmp_sizet, n_spatial2);
      }
    }

    // Inquire label variables
    NETCDF_CHECK(nc_inq_varid(nc_id, "spatial", &spatial_var_id));
    NETCDF_CHECK(nc_inq_var(nc_id, spatial_var_id, varname, &check_vartype, &check_ndims, check_dims, &natts));
    if (check_ndims != 1 || check_dims[0] != spatial_dim_id) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: \"spatial\" variable defined in variable incompatible ndims=%d check_dims=[%d %d %d]\n",
		  ndims, check_dims[0], check_dims[1], check_dims[2]);
    }

    NETCDF_CHECK(nc_inq_varid(nc_id, "cell_spatial", &cell_spatial_var_id));
    NETCDF_CHECK(nc_inq_var(nc_id, cell_spatial_var_id, varname, &check_vartype, &check_ndims, check_dims, &natts));
    if (check_ndims != 1 || check_dims[0] != cell_spatial_dim_id) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: \"cell_spatial\" variable defined in variable incompatible ndims=%d check_dims=[%d %d %d]\n",
		  ndims, check_dims[0], check_dims[1], check_dims[2]);
    }

    NETCDF_CHECK(nc_inq_varid(nc_id, "cell_angular", &cell_angular_var_id));
    NETCDF_CHECK(nc_inq_var(nc_id, cell_angular_var_id, varname, &check_vartype, &check_ndims, check_dims, &natts));
    if (check_ndims != 2 || check_dims[0] != cell_angular_dim_id || check_dims[1] != label_dim_id) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: \"cell_angular\" variable defined in variable incompatible ndims=%d check_dims=[%d %d %d]\n",
		  ndims, check_dims[0], check_dims[1], check_dims[2]);
    }

    NETCDF_CHECK(nc_inq_varid(nc_id, "cell_lengths", &cell_lengths_var_id));
    NETCDF_CHECK(nc_inq_var(nc_id, cell_lengths_var_id, varname, &check_vartype, &check_ndims, check_dims, &natts));
    if (check_ndims != 2 || check_dims[0] != frame_dim_id || check_dims[1] != spatial_dim_id) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: \"cell_lengths\" variable defined in variable incompatible ndims=%d check_dims=[%d %d %d]\n",
		  ndims, check_dims[0], check_dims[1], check_dims[2]);
    }

    NETCDF_CHECK(nc_inq_varid(nc_id, "cell_angles", &cell_angles_var_id));
    NETCDF_CHECK(nc_inq_var(nc_id, cell_angles_var_id, varname, &check_vartype, &check_ndims, check_dims, &natts));
    if (check_ndims != 2 || check_dims[0] != frame_dim_id || check_dims[1] != spatial_dim_id) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: \"cell_angles\" variable defined in variable incompatible ndims=%d check_dims=[%d %d %d]\n",
		  ndims, check_dims[0], check_dims[1], check_dims[2]);
    }

    add_cell_rotated = 0;
    if (nc_inq_varid(nc_id, "cell_rotated", &cell_rotated_var_id))
      add_cell_rotated = 1;
    else { 
      NETCDF_CHECK(nc_inq_var(nc_id, cell_rotated_var_id, varname, &check_vartype, &check_ndims, check_dims, &natts));
      if (check_ndims != 1 || check_dims[0] != frame_dim_id) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: \"cell_rotated\" variable defined in variable incompatible ndims=%d check_dims=[%d %d %d]\n",
			      ndims, check_dims[0], check_dims[1], check_dims[2]);
      }
    }

    add_cell_lattice = 0;
    if (nc_inq_varid(nc_id, "cell_lattice", &cell_lattice_var_id))
      add_cell_lattice = 1;
    else { 
      NETCDF_CHECK(nc_inq_var(nc_id, cell_lattice_var_id, varname, &check_vartype, &check_ndims, check_dims, &natts));
      if (check_ndims != 3 || check_dims[0] != frame_dim_id || check_dims[1] != spatial_dim_id || check_dims[2] != spatial_dim_id) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: \"cell_lattice\" variable defined in variable incompatible ndims=%d check_dims=[%d %d %d]\n",
			      ndims, check_dims[0], check_dims[1], check_dims[2]);      
      }
    }
  }

  nc_redef(nc_id);

  // Add variables cell_rotated and cell_lattice to describe orientation of cell
  // For backwards compatibility, it's not an error if these variables don't exist
  dims[0] = frame_dim_id;
  dims[1] = spatial_dim_id;
  dims[2] = spatial_dim_id;
  if (add_cell_rotated) NETCDF_CHECK(nc_def_var(nc_id, "cell_rotated", NC_INT, 1, dims, &cell_rotated_var_id));
  if (add_cell_lattice) NETCDF_CHECK(nc_def_var(nc_id, "cell_lattice", NC_DOUBLE, 3, dims, &cell_lattice_var_id));
#ifdef NETCDF4
  if (netcdf4) {
    if (add_cell_rotated) NETCDF_CHECK(nc_def_var_deflate(nc_id, cell_rotated_var_id, shuffle, deflate, deflate_level));
    if (add_cell_lattice) NETCDF_CHECK(nc_def_var_deflate(nc_id, cell_lattice_var_id, shuffle, deflate, deflate_level));
  }
#endif
  
  // Define variables for parameters (d=0) and properties (d=1)
  for (d=0; d<2; d++) {
    if (dictionaries[d] == NULL) continue;

    dictionary_get_n(dictionaries[d], &n);
    for (i=1; i <= n; i++) {
      dictionary_query_index(dictionaries[d], &i, key, &type, shape, &data, error, C_KEY_LEN);
      PASS_ERROR;
      if (data == NULL) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: NULL pointer for entry i=%d", i);
      }

      // null-terminate the Fortran string
      key[C_KEY_LEN-1] = '\0'; 
      if (strchr(key, ' ') == NULL) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: key %s not terminated with blank", key);
      }
      *strchr(key, ' ') = '\0';

      // See if we should be writing or skipping this property
      if (d ==1) {
	dictionary_query_key(selected_properties, key, &tmp_type, tmp_shape, &tmp_data, &tmp_error, strlen(key));
	CLEAR_ERROR;
	if (tmp_error != ERROR_NONE) continue;
      }

      if (strcasecmp(key, "Lattice") == 0 || strcasecmp(key, "Properties") == 0) continue;

      // Find equivalent NetCDF type and dimensions

      debug("%s n_spatial %d n_spatial2 %d spatial_dim_id %d spatial2_dim_id %d\n",
	    key, n_spatial, n_spatial2, spatial_dim_id, spatial2_dim_id);

      convert_to_netcdf_type(key, type, shape, frame_dim_id, spatial_dim_id, atom_dim_id,
			     label_dim_id, string_dim_id, spatial2_dim_id, 
			     n_spatial, n_atom, n_label, n_string, n_spatial2,
			     &vartype, &ndims, dims, error);
      PASS_ERROR;

      // Translate names: pos -> coordinates and velo -> velocities
      if (strcmp(key,"pos") == 0)
	strncpy(key,"coordinates",C_KEY_LEN);
      if (strcmp(key,"velo") == 0)
	strncpy(key,"velocities",C_KEY_LEN);    
      

      // Create variable if necessary
      newvar = 0;
      if (nc_inq_varid(nc_id, key, &var_id) != NC_NOERR) {
	debug("write_netcdf: defining variable %s type %d ndims %d dims=[%d %d %d]\n", key, vartype, ndims, dims[0], dims[1], dims[2]);
	NETCDF_CHECK(nc_def_var(nc_id, key, vartype, ndims, dims, &var_id));
	newvar = 1;

	// Set _FillValue for char variables to ' ' rather than default of '\0'
	if (type == T_CHAR || type == T_CHAR_A)  {
	  NETCDF_CHECK(nc_inq_var(nc_id, var_id, varname, &check_vartype, &check_ndims, check_dims, &natts));
	  NETCDF_CHECK(nc_put_att_text(nc_id, var_id, "_FillValue", strlen(fill), fill));
	}

      } else {
	NETCDF_CHECK(nc_inq_var(nc_id, var_id, varname, &check_vartype, &check_ndims, check_dims, &natts));
	if (vartype != check_vartype) {
	  RAISE_ERROR_WITH_KIND(ERROR_IO, "Mismatch in entry \"%s\" - type %d != %d. Perhaps duplicate variable name?\n", key, vartype, check_vartype);
	}
	if (ndims != check_ndims) {
	  RAISE_ERROR_WITH_KIND(ERROR_IO, "Mismatch in entry \"%s\" - ndims %d != %d. Perhaps duplicate variable name?\n", key, ndims, check_ndims);
	  
	}
	for (j=0; j<ndims; j++)
	  if (dims[j] != check_dims[j]) {
	    RAISE_ERROR_WITH_KIND(ERROR_IO, "Mismatch in entry %s - dimension %d ndims=%d dims=[%d %d %d] != check_dims=[%d %d %d]\n", key, j, ndims, dims[0], dims[1], dims[2], 
			check_dims[0], check_dims[1], check_dims[2]);
	  }
      }

      if (newfile || newvar) {
	// Backwards compatibility - add "type" attribute 
	type_att = type;
	if (dims[1] == atom_dim_id) {
	  switch(type) {
	  case(T_INTEGER_A):
	  case(T_INTEGER_A2):
	    type_att = PROPERTY_INT;
	    break;
	  case(T_REAL_A):
	  case(T_REAL_A2):
	    type_att = PROPERTY_REAL;
	    break;
	  case(T_LOGICAL_A):
	    type_att = PROPERTY_LOGICAL;
	    break;
	  case(T_CHAR_A):
	    type_att = PROPERTY_STR;
	    break;
	  default:
	    RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: cannot find property type for dictionary type %d", type);
	  }
	}
	nc_put_att_int(nc_id, var_id, "type", NC_INT, 1, &type_att);
      }
#ifdef NETCDF4
      if (netcdf4 && (newfile || newvar))
	NETCDF_CHECK(nc_def_var_deflate(nc_id, var_id, shuffle, deflate, deflate_level));
#endif
    }
  }

  NETCDF_CHECK(nc_enddef(nc_id));

  if (newfile) {
    // Put data for label variables
    NETCDF_CHECK(nc_put_var_text(nc_id, spatial_var_id, "xyz"));
    NETCDF_CHECK(nc_put_var_text(nc_id, cell_spatial_var_id, "abc"));
    
    start[0] = 0;    start[1] = 0;
    count[0] = 1;    count[1] = strlen("alpha");
    NETCDF_CHECK(nc_put_vara_text(nc_id, cell_angular_var_id, start, count, "alpha"));
    
    start[0] = 1;   start[1] = 0;
    count[0] = 1;   count[1] = strlen("beta");
    NETCDF_CHECK(nc_put_vara_text(nc_id, cell_angular_var_id, start, count, "beta"));
    
    start[0] = 2;   start[1] = 0;
    count[0] = 1;   count[1] = strlen("gamma");
    NETCDF_CHECK(nc_put_vara_text(nc_id, cell_angular_var_id, start, count, "gamma"));
  }

  // Put cell_lengths and cell_angles
  for (i=0; i<3; i++)
    cell_angles[i] *= RAD_TO_DEG;
  start[0] = frame;
  start[1] = 0;
  count[0] = 1;
  count[1] = 3;
  NETCDF_CHECK(nc_put_vara_double(nc_id, cell_lengths_var_id, start, count, cell_lengths));
  NETCDF_CHECK(nc_put_vara_double(nc_id, cell_angles_var_id, start, count, cell_angles));

  // Have atomic positions been rotated to align cell vector a with x axis?
  start[0] = frame;
  count[0] = 1;
  NETCDF_CHECK(nc_put_vara_int(nc_id, cell_rotated_var_id, start, count, &cell_rotated));
  
  // Also save original lattice so rotation can be undone when file is read
  start[2] = 0;
  count[2] = 3;
  NETCDF_CHECK(nc_put_vara_double(nc_id, cell_lattice_var_id, start, count, &(lattice[0][0])));

  // Put variables for parameters (d=0) and properties (d=1)
  for (d=0; d<2; d++) {
    if (dictionaries[d] == NULL) continue;

    dictionary_get_n(dictionaries[d], &n);
    for (i=1; i <= n; i++) {
      dictionary_query_index(dictionaries[d], &i, key, &type, shape, &data, error, C_KEY_LEN);
      PASS_ERROR;
      if (data == NULL) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: NULL pointer for entry i=%d", i);
      }

      // null-terminate the Fortran string
      key[C_KEY_LEN-1] = '\0'; 
      if (strchr(key, ' ') == NULL) {
	RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: key %s not terminated with blank", key);
      }
      *strchr(key, ' ') = '\0';

      // See if we should be writing or skipping this property
      if (d ==1) {
	dictionary_query_key(selected_properties, key, &tmp_type, tmp_shape, &tmp_data, &tmp_error, strlen(key));
	if (tmp_error != ERROR_NONE) continue;
      }
      
      if (strcasecmp(key, "Lattice") == 0 || strcasecmp(key, "Properties") == 0) continue;

      if (strcmp(key,"pos") == 0)
	strncpy(key,"coordinates",C_KEY_LEN);
      if (strcmp(key,"velo") == 0)
	strncpy(key,"velocities",C_KEY_LEN);    

      // Find variable ID
      NETCDF_CHECK(nc_inq_varid(nc_id, key, &var_id));

      //if (skip_variable(key,type,shape)) continue;

      start[0] = frame;
      start[1] = 0;
      start[2] = 0;

      count[0] = 1;
      // If it's a 2D array we need to transpose shape
      if (type == T_INTEGER_A || type == T_REAL_A || type == T_LOGICAL_A || type == T_CHAR) {
	count[1] = shape[0];
      } else {
	count[1] = shape[1];
	count[2] = shape[0];
      }

      debug("write_netcdf: writing variable \"%s\" start=[%ld %ld %ld] count=[%ld %ld %ld]\n", key, start[0], start[1], start[2], count[0], count[1], count[2]);

      // Do the writing.
      switch(type) {
      case(T_INTEGER):
      case (T_LOGICAL):
	NETCDF_CHECK(nc_put_vara_int(nc_id, var_id, start, count, (int *)data));
	break;
      case(T_REAL):
	NETCDF_CHECK(nc_put_vara_double(nc_id, var_id, start, count, (double *)data));
	break;
      case(T_CHAR):
	NETCDF_CHECK(nc_put_vara_text(nc_id, var_id, start, count, (char *)data));
	break;
      case(T_INTEGER_A):
      case(T_LOGICAL_A):
	NETCDF_CHECK(nc_put_vara_int(nc_id, var_id, start, count, (int *)data));
	break;
      case(T_REAL_A):
	NETCDF_CHECK(nc_put_vara_double(nc_id, var_id, start, count, (double *)data));
	break;
      case(T_INTEGER_A2):
	NETCDF_CHECK(nc_put_vara_int(nc_id, var_id, start, count, (int *)data));
	break;
      case(T_REAL_A2):
	NETCDF_CHECK(nc_put_vara_double(nc_id, var_id, start, count, (double *)data));
	break;      
      case(T_CHAR_A):
	NETCDF_CHECK(nc_put_vara_text(nc_id, var_id, start, count, (char *)data));
	break;      
      default:
	RAISE_ERROR_WITH_KIND(ERROR_IO, "write_netcdf: unknown variable %s type %d\n", key, type);
      }
    }
  }

  NETCDF_CHECK(nc_close(nc_id));
}

void query_netcdf (char *filename, int *n_frame, int *n_atom, int *n_label, int *n_string, int *error)
{
  int nc_id, retval;
  int frame_dim_id, spatial_dim_id, atom_dim_id, cell_spatial_dim_id,
    cell_angular_dim_id, label_dim_id, string_dim_id, n_spatial, spatial2_dim_id, n_spatial2;
  size_t tmp_sizet;

  INIT_ERROR;

  NETCDF_CHECK(nc_open(filename, NC_NOWRITE, &nc_id));

  // Inquire dimensions
  NETCDF_CHECK(nc_inq_dimid(nc_id, "frame", &frame_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "spatial", &spatial_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "atom", &atom_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "cell_spatial", &cell_spatial_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "cell_angular", &cell_angular_dim_id));
  NETCDF_CHECK(nc_inq_dimid(nc_id, "label", &label_dim_id));
  if (nc_inq_dimid(nc_id, "string", &string_dim_id) != NC_NOERR) {
    // No strings in this file
    string_dim_id = 0;
  }
  if (nc_inq_dimid(nc_id, "spatial2", &spatial2_dim_id) != NC_NOERR) {
    // No spatial2s in this file
    spatial2_dim_id = 0;
  }


  // Get sizes of dimensions
  NETCDF_CHECK(nc_inq_dimlen(nc_id, spatial_dim_id, &tmp_sizet));
  n_spatial = (int)tmp_sizet;
  if (n_spatial != 3) {
    RAISE_ERROR_WITH_KIND(ERROR_IO, "query_netcdf: number of spatial dimensions = %d != 3", n_spatial);
  }

  NETCDF_CHECK(nc_inq_dimlen(nc_id, frame_dim_id, &tmp_sizet));
  *n_frame = (int)tmp_sizet;
  NETCDF_CHECK(nc_inq_dimlen(nc_id, atom_dim_id, &tmp_sizet));
  *n_atom = (int)tmp_sizet;
  NETCDF_CHECK(nc_inq_dimlen(nc_id, label_dim_id, &tmp_sizet));
  *n_label = (int)tmp_sizet;
  if (string_dim_id != 0) {
    NETCDF_CHECK(nc_inq_dimlen(nc_id, string_dim_id, &tmp_sizet));
    *n_string = (int)tmp_sizet;
  }
  else {
    *n_string = 0;
  }
  if (spatial2_dim_id != 0) {
    NETCDF_CHECK(nc_inq_dimlen(nc_id, spatial2_dim_id, &tmp_sizet));
    n_spatial2 = (int)tmp_sizet;
    if (n_spatial2 != 9) {
      RAISE_ERROR_WITH_KIND(ERROR_IO, "query_netcdf: number of spatial dimensions**2 = %d != 9", n_spatial2);
    }
  }
  else {
    n_spatial2 = 0;
  }

  debug("query_netcdf: dimension information:\n");
  debug("  frame_dim_id = %d, n_frame = %d\n", frame_dim_id, *n_frame);
  debug("  spatial_dim_id = %d, n_spatial = %d\n", spatial_dim_id, n_spatial);
  debug("  atom_dim_id = %d, n_atom = %d\n", atom_dim_id, *n_atom);
  debug("  cell_spatial_dim_id = %d\n", cell_spatial_dim_id);
  debug("  cell_angular_dim_id = %d\n", cell_angular_dim_id);
  debug("  label_dim_id = %d, n_label = %d\n", label_dim_id, *n_label);
  debug("  string_dim_id = %d, n_string = %d\n", string_dim_id, *n_string);
  debug("  spatial2_dim_id = %d, n_spatial2 = %d\n\n", spatial2_dim_id, n_spatial2);

  NETCDF_CHECK(nc_close(nc_id));
}


#else
void read_netcdf (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3], 
		  double cell_lengths[3], double cell_angles[3], int *cell_rotated,
		  int *n_atom, int frame, int zero, int *range, int irep, double rrep, int *error)
{
  INIT_ERROR;
  RAISE_ERROR( "No NetCDF support compiled in. Recompile with HAVE_NETCDF=1. \n");
}

void write_netcdf (char *filename, fortran_t *params, fortran_t *properties, fortran_t *selected_properties, double lattice[3][3],
		   double cell_lengths[3], double cell_angles[3], int cell_rotated, 
		   int n_atom, int n_label, int n_string, int frame, int netcdf4, int append,
		   int shuffle, int deflate, int deflate_level, int *error)

{
  INIT_ERROR;
  RAISE_ERROR( "No NetCDF support compiled in. Recompile with HAVE_NETCDF=1. \n");
}

void query_netcdf (char *filename, int *n_frame, int *n_atom, int *n_label, int *n_string, int *error)
{
  INIT_ERROR;
  RAISE_ERROR( "No NetCDF support compiled in. Recompile with HAVE_NETCDF=1. \n");
}
#endif
