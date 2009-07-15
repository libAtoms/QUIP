/* Tools to convert between atomic trajectories in extended XYZ file format and
 * NetCDF format using the AMBER conventions. 
 *
 * To build main program compile with -DHAVE_NETCDF -DMAIN_PROGRAM and link with
 * NetCDF library version 3 or above. There's also a Fortran 95 API which calls
 * the cio* functions at the bottom of this file. This can be found in CInOuput.f95.
 *
 * If you're only interested in XYZ support, you can compile without -DHAVE_NETCDF.
 * For the Fortran API don't need to include the main program, so don't need
 * -DMAIN_PROGRAM.
 *
 * (c) James Kermode <jrk33@cam.ac.uk>
 * December 2008
 *
 */

#ifdef HAVE_NETCDF
#include <netcdf.h>
#endif

#ifdef DEBUG
#define debug(fmt, ...) fprintf(stderr, fmt, ## __VA_ARGS__)
#else
#define debug(fmt, ...) 
#endif

/* required for isblank() */
#ifndef __USE_ISOC99
#define __USE_ISOC99
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

#include "xyz_netcdf.h"

const char PROPERTY_STRINGS[5] = {'\0','I','R','S','L'};

const char *PROPERTY_TYPES[] = {"","int", "real", "string", "logical"};
const char *PARAM_TYPES[] = {"(none)", "int", "real", "complex", "logical",
			  "int array", "real array", "complex array",
			  "logical array", "string", "string array"};

#ifndef HAVE_ATOMEYE
/* pe() and vstrf() functions copied from AtomEye (c) Ju Li */

#define STRF_MAX_CHAR  1023

/* static memory driver of vsnprintf() */
char *vstrf (char *format, va_list ap)
{
    static char buffer[STRF_MAX_CHAR+1];
    vsnprintf(buffer, STRF_MAX_CHAR+1, format, ap);
    return (buffer);
} /* end strf() */

/* perror() driver using strf() */
void pe (char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    fprintf (stderr, "error: ");
    if (errno==0) vfprintf(stderr,format,ap);
    else perror(vstrf(format,ap));
    va_end(ap);
    //exit(1);
    abort();
    return;
} /* end pe() */
#endif

#ifdef HAVE_NETCDF
#ifdef MAIN_PROGRAM
#define netcdf_check(s) if ((retval = (s))) pe("NetCDF Error: %s %d %s\n", __FILE__, __LINE__, nc_strerror(retval));
#else
#define netcdf_check(s) if ((retval = (s))) { fprintf(stderr,"NetCDF Error: %s %d %s\n", __FILE__, __LINE__, nc_strerror(retval)); return 0; }
#endif
#endif

void atoms_init(Atoms *atoms) {
  atoms->int_data = NULL;
  atoms->real_data = NULL;
  atoms->str_data = NULL;
  atoms->logical_data = NULL;
  atoms->filter = NULL;
  atoms->initialised = 0;
  atoms->nc_in = -1;
  atoms->nc_out = -1;
}

void atoms_alloc(Atoms *atoms) {
  atoms->int_data = malloc(atoms->n_atom*atoms->n_int*sizeof(int));
  if (atoms->int_data == NULL) pe("Error allocating atoms->int_data");

  atoms->real_data = malloc(atoms->n_atom*atoms->n_real*sizeof(double));
  if (atoms->real_data == NULL) pe("Error allocating atoms->real_data");
  
  debug("allocating str_data %d*%d*%d\n", atoms->n_atom, atoms->n_str, PROPERTY_STRING_LENGTH);
  atoms->str_data = malloc(atoms->n_atom*atoms->n_str*PROPERTY_STRING_LENGTH*sizeof(char));
  if (atoms->str_data == NULL) pe("Error allocating atoms->str_data");

  atoms->logical_data = malloc(atoms->n_atom*atoms->n_logical*sizeof(int));
  if (atoms->logical_data == NULL) pe("Error allocating atoms->logical_data");
}

void atoms_free(Atoms *atoms) {
  if (atoms->int_data != NULL) {
    free(atoms->int_data);
    atoms->int_data     = NULL;
  }
  if (atoms->real_data != NULL) {
    free(atoms->real_data);
    atoms->real_data    = NULL;
  }
  if (atoms->str_data != NULL) {
    free(atoms->str_data);
    atoms->str_data     = NULL;
  }
  if (atoms->logical_data != NULL) {  
    free(atoms->logical_data);
    atoms->logical_data = NULL;
  }
}

void atoms_realloc(Atoms *atoms) {
  atoms_free(atoms);
  atoms_alloc(atoms);
}

void atoms_swap_properties(Atoms *atoms, int i, int j)
{
  int tmp_type, tmp_start, tmp_ncols, tmp_var_id1, tmp_var_id2, tmp_filter;
  char tmp_name[MAX_PROP_LENGTH];

  if (i != j) {
    tmp_type = atoms->property_type[j];
    tmp_start = atoms->property_start[j];
    tmp_ncols = atoms->property_ncols[j];
    tmp_var_id1 = atoms->property_var_id[j][NETCDF_IN];
    tmp_var_id2 = atoms->property_var_id[j][NETCDF_OUT];
    tmp_filter = atoms->property_filter[j];
    strcpy(tmp_name, atoms->property_name[j]);

    atoms->property_type[j] = atoms->property_type[i];
    atoms->property_start[j] = atoms->property_start[i];
    atoms->property_ncols[j] = atoms->property_ncols[i];
    atoms->property_var_id[j][NETCDF_IN] = atoms->property_var_id[i][NETCDF_IN];
    atoms->property_var_id[j][NETCDF_OUT] = atoms->property_var_id[i][NETCDF_OUT];
    atoms->property_filter[j] = atoms->property_filter[i];
    strcpy(atoms->property_name[j], atoms->property_name[i]);

    atoms->property_type[i] = tmp_type;
    atoms->property_start[i] = tmp_start;
    atoms->property_ncols[i] = tmp_ncols;
    atoms->property_var_id[i][NETCDF_IN] = tmp_var_id1;
    atoms->property_var_id[i][NETCDF_OUT] = tmp_var_id2;
    atoms->property_filter[i] = tmp_filter;
    strcpy(atoms->property_name[i], tmp_name);
  }
}

int atoms_add_param(Atoms *atoms, char *key, int type, int size, int var_id, int inout)
{
  int idx;
  
  idx = atoms->n_param;
  strcpy(atoms->param_key[idx],key);
  atoms->param_type[idx] = type;
  atoms->param_size[idx] = size;
  atoms->param_filter[idx] = 1;
  atoms->param_var_id[idx][inout] = var_id;
  atoms->n_param++;
  if (atoms->n_param == MAX_PARAM_COUNT) 
    pe("MAX_PARAM_COUNT (%d) exceeded", MAX_PARAM_COUNT);
 
  return idx;
}

int atoms_add_property(Atoms *atoms, char *key, int type, int ncols, int var_id, int inout)
{
  int idx;
  
  idx = atoms->n_property;
  strcpy(atoms->property_name[idx],key);
  atoms->property_type[idx] = type;
  atoms->property_ncols[idx] = ncols;
  atoms->property_filter[idx] = 1;
  atoms->property_var_id[idx][inout] = var_id;

  switch(atoms->property_type[idx]) {
  case(0):
    // do nothing
    break;
  case(PROPERTY_INT):
    atoms->property_start[idx] = atoms->n_int;
    atoms->n_int += ncols;
    break;
  case(PROPERTY_REAL):
    atoms->property_start[idx] = atoms->n_real;
    atoms->n_real += ncols;
    break;
  case(PROPERTY_STR):
    atoms->property_start[idx] = atoms->n_str;
    atoms->n_str += ncols;
    break;
  case(PROPERTY_LOGICAL):
    atoms->property_start[idx] = atoms->n_logical;
    atoms->n_logical += ncols;
    break;
  default:
    pe("Unknown property type %d", atoms->property_type[atoms->n_property]);
  }

  atoms->n_property++;
  if (atoms->n_property == MAX_ENTRY_COUNT) 
    pe("MAX_ENTRY_COUNT (%d) exceeded", MAX_ENTRY_COUNT);
 
  return idx;
}

int atoms_find_property(Atoms *atoms, char *key) 
{
  int i;

  for (i=0; i < atoms->n_property; i++)
    if (strcasecmp(atoms->property_name[i],key) == 0) return i;
  
  return -1;
}

int atoms_find_param(Atoms *atoms, char *key) 
{
  int i;

  for (i=0; i < atoms->n_param; i++)
    if (strcasecmp(atoms->param_key[i],key) == 0) return i;
  
  return -1;
}

void lattice_abc_to_xyz(double cell_lengths[3], double cell_angles[3],
			double lattice[3][3])
{
  double a, b, c, alpha, beta, gamma, cos_alpha, cos2_alpha, cos_beta, cos2_beta, cos_gamma, cos2_gamma,
    sin_gamma, sin2_gamma;

  a = cell_lengths[0];  b = cell_lengths[1]; c = cell_lengths[2];
  alpha = cell_angles[0]*M_PI/180.0; beta = cell_angles[1]*M_PI/180.0; gamma = cell_angles[2]*M_PI/180.0;
    
  cos_alpha = cos(alpha); cos2_alpha = cos_alpha*cos_alpha;
  cos_beta  = cos(beta);  cos2_beta  = cos_beta *cos_beta;
  cos_gamma = cos(gamma); cos2_gamma = cos_gamma*cos_gamma;
  sin_gamma = sin(gamma); sin2_gamma = sin_gamma*sin_gamma;

  lattice[0][0] = a;
  lattice[0][1] = 0.0;
  lattice[0][2] = 0.0;

  lattice[1][0] = b * cos_gamma;
  lattice[1][1] = b * sin_gamma;
  lattice[1][2] = 0.0;

  lattice[2][0] = c * cos_beta;
  lattice[2][1] = c * (cos_alpha - cos_beta*cos_gamma) / sin_gamma;
  lattice[2][2] = c * sqrt(1.0 - (cos2_alpha + cos2_beta - 2.0*cos_alpha*cos_beta*cos_gamma)/ sin2_gamma);
}

void lattice_xyz_to_abc(double lattice[3][3], double cell_lengths[3], double cell_angles[3])
{
  int i;

  for (i=0; i<3; i++)
    cell_lengths[i] = sqrt(lattice[i][0]*lattice[i][0] + 
			   lattice[i][1]*lattice[i][1] + 
			   lattice[i][2]*lattice[i][2]);

  cell_angles[0] = 180.0/M_PI*acos((lattice[1][0]*lattice[2][0] + 
				    lattice[1][1]*lattice[2][1] +
				    lattice[1][2]*lattice[2][2])/
				   (cell_lengths[1]*cell_lengths[2]));
		       
  cell_angles[1] = 180.0/M_PI*acos((lattice[2][0]*lattice[0][0] + 
				    lattice[2][1]*lattice[0][1] +
				    lattice[2][2]*lattice[0][2])/
				   (cell_lengths[2]*cell_lengths[0]));

  cell_angles[2] = 180.0/M_PI*acos((lattice[0][0]*lattice[1][0] + 
				    lattice[0][1]*lattice[1][1] +
				    lattice[0][2]*lattice[1][2])/
				   (cell_lengths[0]*cell_lengths[1]));
}


#ifdef HAVE_NETCDF
 
void replace_fill_values(Atoms *atoms, int irep, double rrep) {
  int i, j=0, n;

  // parameters
  for (i=0; i<atoms->n_param; i++) {
    switch(atoms->param_type[i]) {
    case(T_INTEGER):
      if (atoms->param_int[i] == NC_FILL_INT)
	atoms->param_int[i] = irep;
      break;
    case(T_REAL):
      if((atoms->param_real[i] > 0) == (NC_FILL_DOUBLE > 0) && /* prevents potential overflow */
	 (absval(atoms->param_real[i] - NC_FILL_DOUBLE) <= absval(DBL_EPSILON * NC_FILL_DOUBLE)))
	atoms->param_real[i] = rrep;
      break;
    case(T_INTEGER_A):
      for (j=0; j<3; j++)
	if (atoms->param_int_a[i][j] == NC_FILL_INT)
	  atoms->param_int_a[i][j] = irep;
      break;
    case(T_REAL_A):
      if((atoms->param_real_a[i][j] > 0) == (NC_FILL_DOUBLE > 0) && /* prevents potential overflow */
	 (absval(atoms->param_real_a[i][j] - NC_FILL_DOUBLE) <= absval(DBL_EPSILON * NC_FILL_DOUBLE)))
	atoms->param_real_a[i][j] = rrep;
      break;
    }
  }

  // properties
  for (i=0; i<atoms->n_property; i++) {
    switch(atoms->property_type[i]) {
    case(PROPERTY_INT):
      for (n=0; n<atoms->n_atom; n++)
	for (j=0; j<atoms->property_ncols[i]; j++)
	  if (property_int(atoms,i,j,n) == NC_FILL_INT)
	    property_int(atoms,i,j,n) = irep;
      break;

    case(PROPERTY_LOGICAL):
      for (n=0; n<atoms->n_atom; n++)
	for (j=0; j<atoms->property_ncols[i]; j++)
	  if (property_logical(atoms,i,j,n) == NC_FILL_INT)
	    property_logical(atoms,i,j,n) = irep;
      break;

    case(T_REAL):
      for (n=0; n<atoms->n_atom; n++)
	for (j=0; j<atoms->property_ncols[i]; j++)
	if((property_real(atoms,i,j,n) > 0) == (NC_FILL_DOUBLE > 0) && /* prevents potential overflow */
	   (absval(property_real(atoms,i,j,n) - NC_FILL_DOUBLE) <= absval(DBL_EPSILON * NC_FILL_DOUBLE))) {
	  property_real(atoms,i,j,n) = rrep;
	}
      break;
    }
  }
}


int read_netcdf (int nc_id, Atoms *atoms, int frame, int *atomlist, int natomlist, int query, 
		 int redefine, int realloc, int replacefill, int irep, double rrep)
{
  int i,j,n;
  int retval, nvars;
  size_t start1[1], start2[2], count2[2],start3[3],count3[3];
  double cell_lengths[3], cell_angles[3];
  char varname[NC_MAX_NAME+1];
  nc_type vartype;
  int ndims, dimids[NC_MAX_VAR_DIMS], natts;
  int *tmpint, *tmpint3, original_index=0, type;
  double *tmpreal, *tmpreal3;
  char *tmpchar;

  if (!atoms->initialised) {
    // get dimension and variable information

    debug("read_netcdf: nc_id = %d\n", nc_id);

    // Inquire dimensions
    netcdf_check(nc_inq_dimid(nc_id, "frame", &atoms->frame_dim_id[NETCDF_IN]));
    netcdf_check(nc_inq_dimid(nc_id, "spatial", &atoms->spatial_dim_id[NETCDF_IN]));
    netcdf_check(nc_inq_dimid(nc_id, "atom", &atoms->atom_dim_id[NETCDF_IN]));
    netcdf_check(nc_inq_dimid(nc_id, "cell_spatial", &atoms->cell_spatial_dim_id[NETCDF_IN]));
    netcdf_check(nc_inq_dimid(nc_id, "cell_angular", &atoms->cell_angular_dim_id[NETCDF_IN]));
    netcdf_check(nc_inq_dimid(nc_id, "label", &atoms->label_dim_id[NETCDF_IN]));
    netcdf_check(nc_inq_dimid(nc_id, "string", &atoms->string_dim_id[NETCDF_IN]));

    // Sizes of dimensions
    netcdf_check(nc_inq_dimlen(nc_id, atoms->frame_dim_id[NETCDF_IN], &atoms->n_frame));
    netcdf_check(nc_inq_dimlen(nc_id, atoms->atom_dim_id[NETCDF_IN], &atoms->n_atom_total));
    netcdf_check(nc_inq_dimlen(nc_id, atoms->label_dim_id[NETCDF_IN], &atoms->n_label));
    netcdf_check(nc_inq_dimlen(nc_id, atoms->string_dim_id[NETCDF_IN], &atoms->n_string));

    // Check string lengths (fixed at compile time for now...)
    if (atoms->n_label != PROPERTY_STRING_LENGTH) {
      fprintf(stderr,"atoms->n_label(%d) != PROPERTY_STRING_LENGTH(%d)\n", (int)atoms->n_label, PROPERTY_STRING_LENGTH);
      return 0;
    }
    if (atoms->n_string != PARAM_STRING_LENGTH) {
      fprintf(stderr,"atoms->n_string(%d) != PARAM_STRING_LENGTH(%d)\n", (int)atoms->n_string, PARAM_STRING_LENGTH);
      return 0;
    }

    atoms->n_property = 0;
    atoms->n_param = 0;
    atoms->n_int = 0;
    atoms->n_real = 0;
    atoms->n_str = 0;
    atoms->n_logical = 0;

    // Loop over all variables in file
    netcdf_check(nc_inq_nvars(nc_id, &nvars));
    for (i=0; i<nvars; i++) {
      netcdf_check(nc_inq_var(nc_id, i, varname, &vartype, &ndims, dimids, &natts));
      
      if (strcmp(varname,"cell_lengths") == 0) {
	atoms->cell_lengths_var_id[NETCDF_IN] = i;
	continue;
      }

      if (strcmp(varname,"cell_angles") == 0) {
	atoms->cell_angles_var_id[NETCDF_IN] = i;
	continue;
      }

      if (strcmp(varname, "spatial") == 0 || strcmp(varname, "cell_spatial") == 0 ||
	  strcmp(varname, "cell_angular") == 0)
	continue;

      netcdf_check(nc_get_att_int(nc_id, i, "type", &type));

      if (ndims == 0) 
	; // one off variable, do nothing
      else if (ndims == 1) {
	if (dimids[0] == atoms->frame_dim_id[NETCDF_IN]) {
	  // scalar per-frame parameter
	  atoms_add_param(atoms, varname, type, 1, i, NETCDF_IN);
	}
	else {
	  fprintf(stderr,"Unknown one dimensional variable %s\n", varname);
	  return 0;
	}
      }
      else if (ndims == 2) {
	if (dimids[0] == atoms->frame_dim_id[NETCDF_IN] && dimids[1] == atoms->atom_dim_id[NETCDF_IN]) {
	  // scalar per atom property
	  atoms_add_property(atoms, varname, type, 1, i, NETCDF_IN);
	}
	else if (dimids[0] == atoms->frame_dim_id[NETCDF_IN] && dimids[1] == atoms->spatial_dim_id[NETCDF_IN]) {
	  // vector per-frame paramter
	  atoms_add_param(atoms, varname, type, 3, i, NETCDF_IN);
	} else if (dimids[0] == atoms->frame_dim_id[NETCDF_IN] && dimids[1] == atoms->string_dim_id[NETCDF_IN]) {
	  // string per-frame parameter
	  atoms_add_param(atoms, varname, type, 1, i, NETCDF_IN);
	}
	else {
	  fprintf(stderr,"Unknown two dimensional variable %s\n", varname);
	  return 0;
	}
      }
      else if (ndims == 3) {
	if (dimids[0] == atoms->frame_dim_id[NETCDF_IN] && dimids[1] == atoms->atom_dim_id[NETCDF_IN]) {
	  
	  if (dimids[2] == atoms->label_dim_id[NETCDF_IN]) {
	    // per atom string property
	    atoms_add_property(atoms, varname, PROPERTY_STR, 1, i, NETCDF_IN);
	  }
	  else if (dimids[2] == atoms->spatial_dim_id[NETCDF_IN]) {
	    // vector per atom property
	    atoms_add_property(atoms, varname, type, 3, i, NETCDF_IN);
	  }
	}
	else {
	  fprintf(stderr,"Unknown three dimensional variable %s\n", varname);
	  return 0;
	}
      }
      else {
	fprintf(stderr,"Unknown variable %s, dimension %d\n", varname, ndims);
	return 0;
      }
    }

    // Translate names: coordinates -> pos and velocities -> velo
    for (i=0; i<atoms->n_property; i++) {
      if (strcmp(atoms->property_name[i],"coordinates") == 0)
	strcpy(atoms->property_name[i],"pos");
      if (strcmp(atoms->property_name[i],"velocities") == 0)
	strcpy(atoms->property_name[i],"velo");
    }

    debug("got %d frames, %d atoms %d props and %d params\n", atoms->n_frame, atoms->n_atom_total, 
	   atoms->n_property, atoms->n_param);

    if (atoms->filter != NULL) free(atoms->filter);
    atoms->filter = malloc(atoms->n_atom_total*sizeof(int));
    if (atoms->filter == NULL) {
      fprintf(stderr,"Error allocating atoms->filter\n");
      return 0;
    }

    if (atomlist) {
      for (i=0; i<atoms->n_atom_total; i++) atoms->filter[i] = 0;
      for (i=0; i<natomlist; i++) {
	if(atomlist[i] < 0 || atomlist[i] >= atoms->n_atom_total) {
	  fprintf(stderr,"filter atom %d out of range 0 < i < %d\n", atomlist[i], (int)atoms->n_atom_total);
	  return 0;
	}
	atoms->filter[atomlist[i]] = 1;
      }
      atoms->n_atom = natomlist;
      debug("filter applied. Selected %d/%d atoms.\n", natomlist, atoms->n_atom_total);
      // Add an int property for original index
      original_index = atoms_add_property(atoms, "original_index", PROPERTY_INT, 1, -1, NETCDF_IN);
    }
    else {
      for (i=0; i<atoms->n_atom_total; i++) atoms->filter[i] = 1;
      atoms->n_atom = atoms->n_atom_total;
    }

    debug("allocating %d atoms\n", atoms->n_atom);

    // reallocate to correct size
    if (realloc) atoms_realloc(atoms);
    atoms->initialised = 1;

    debug("done alloc\n");

    // Fill in original_index property
    if (atomlist) {
      n = 0;
      for (i=0; i<atoms->n_atom_total; i++) {
	if (!atoms->filter[i]) continue;
	//atoms->int_data[atoms->property_start[original_index]*atoms->n_atom + n] = i;
	property_int(atoms,original_index,0,n) = i;
	n++;
      }
    }
    
  } /* end if(nc_id == -1) */
  
  if (query) return atoms->n_frame;
  
  debug("allocating tmps\n");

  // Temporary arrays for netcdf reading when filtering
  //  if (atomlist) {
    tmpint = malloc(sizeof(int)*atoms->n_atom_total);
    tmpint3 = malloc(sizeof(int)*atoms->n_atom_total*3);
    tmpreal = malloc(sizeof(double)*atoms->n_atom_total);
    tmpreal3 = malloc(sizeof(double)*atoms->n_atom_total*3);
    tmpchar = malloc(sizeof(char)*atoms->n_atom_total*PROPERTY_STRING_LENGTH);

    if (tmpint == NULL || tmpreal == NULL || tmpchar == NULL) {
      fprintf(stderr,"Error allocating memory in read_netcdf\n");
      return 0;
    }
    //}

  debug("reading lattice\n");

  // Load lattice
  start2[0] = frame;
  start2[1] = 0;
  count2[0] = 1;
  count2[1] = 3;
  netcdf_check(nc_get_vara_double(nc_id, atoms->cell_lengths_var_id[NETCDF_IN], start2, count2, cell_lengths));
  netcdf_check(nc_get_vara_double(nc_id, atoms->cell_angles_var_id[NETCDF_IN], start2, count2, cell_angles));
  lattice_abc_to_xyz(cell_lengths, cell_angles, atoms->lattice);
  debug("reading params\n");

  // read parameters
  for (i=0; i<atoms->n_param; i++) {
    if (strcasecmp(atoms->param_key[i], "Lattice") == 0 ||
	strcasecmp(atoms->param_key[i], "Properties") == 0) continue;
    
    switch(atoms->param_type[i]) {
    case(T_INTEGER):
      start1[0] = frame;
      netcdf_check(nc_get_var1_int(nc_id, atoms->param_var_id[i][NETCDF_IN], start2, &atoms->param_int[i]));
      break;
    case(T_REAL):
      start1[0] = frame;
      netcdf_check(nc_get_var1_double(nc_id, atoms->param_var_id[i][NETCDF_IN], start2, &atoms->param_real[i]));
      break;
    case(T_CHAR):
      start2[0] = frame;
      start2[1] = 0;
      count2[0] = 1;
      count2[1] = PARAM_STRING_LENGTH;
      netcdf_check(nc_get_vara_text(nc_id, atoms->param_var_id[i][NETCDF_IN], start2, count2, atoms->param_value[i]));
      break;
    case(T_INTEGER_A):
      start2[0] = frame;
      start2[1] = 0;
      count2[0] = 1;
      count2[1] = 3;
      netcdf_check(nc_get_vara_int(nc_id, atoms->param_var_id[i][NETCDF_IN], start2, count2, &atoms->param_int_a[i][0]));
      break;
    case(T_REAL_A):
      start2[0] = frame;
      start2[1] = 0;
      count2[0] = 1;
      count2[1] = 3;
      netcdf_check(nc_get_vara_double(nc_id, atoms->param_var_id[i][NETCDF_IN], start2, count2, &atoms->param_real_a[i][0]));
      break;
    default:
      fprintf(stderr,"Unknown parameter %s type %d\n", atoms->param_key[i], atoms->param_type[i]);
      return 0;
    }
  }

  debug("reading properties\n");

  // Zero string data - this is needed as we're not null terminating our strings
  memset(atoms->str_data, ' ', atoms->n_atom*atoms->n_str*PROPERTY_STRING_LENGTH);

  // read properties
  for (i=0; i<atoms->n_property; i++) {
    
    if (atomlist && strcmp(atoms->property_name[i], "original_index") == 0)
      continue;

    debug("%d %s\n", i, atoms->property_name[i]);

    if (atoms->property_type[i] == PROPERTY_REAL) {
      if (atoms->property_ncols[i] == 1) {
	start2[0] = frame;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = atoms->n_atom_total;
	//if (atomlist) {
	  netcdf_check(nc_get_vara_double(nc_id, atoms->property_var_id[i][NETCDF_IN], start2, count2, tmpreal));
	  n = 0;
	  for (j=0; j<atoms->n_atom_total; j++) {
	    if (!atoms->filter[j]) continue;
	    //atoms->real_data[(atoms->property_start[i])*atoms->n_atom + n] = tmpreal[j];
	    property_real(atoms, i, 0, n) = tmpreal[j];
	    n++;
	  }
	  //}
	  //else
	  //netcdf_check(nc_get_vara_double(nc_id, atoms->property_var_id[i][NETCDF_IN], start2, count2, 
	  //				  &atoms->real_data[atoms->property_start[i]*atoms->n_atom]));
	  //netcdf_check(nc_get_vara_double(nc_id, atoms->property_var_id[i][NETCDF_IN], start2, count2, 
	  //				  &property_real(atoms,i,0,0)));
      } else {
/* 	for (k=0; k<3; k++) { */
/* 	  start3[0] = frame; */
/* 	  start3[1] = 0; */
/* 	  start3[2] = k; */
/* 	  count3[0] = 1; */
/* 	  count3[1] = atoms->n_atom_total; */
/* 	  count3[2] = 1; */
/* 	  //if (atomlist) { */
/* 	    netcdf_check(nc_get_vara_double(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3, tmpreal)); */
/* 	    n = 0; */
/* 	    for (j=0; j<atoms->n_atom_total; j++) { */
/* 	      if (!atoms->filter[j]) continue; */
/* 	      //atoms->real_data[(atoms->property_start[i]+k)*atoms->n_atom + n] = tmpreal[j]; */
/* 	      property_real(atoms, i, k, n) = tmpreal[j]; */
/* 	      /\*debug("setting %s %d %d = %f %f\n", atoms->property_name[i], atoms->property_start[i]+k, n, tmpreal[j], */
/* 		atoms->real_data[(atoms->property_start[i]+k)*atoms->n_atom + n]); *\/ */
/* 	      n++; */
/* 	    } */
/* 	    //} */
/* 	    //else */
/* 	    //netcdf_check(nc_get_vara_double(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3,  */
/* 	    //				    &atoms->real_data[(atoms->property_start[i] + k)*atoms->n_atom])); */
/* 	    //netcdf_check(nc_get_vara_double(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3,  */
/* 	    //				    &property_real(atoms, i, k, 0))); */

	start3[0] = frame;
	start3[1] = 0;
	start3[2] = 0;
	count3[0] = 1;
	count3[1] = atoms->n_atom_total;
	count3[2] = 3;
	netcdf_check(nc_get_vara_double(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3, tmpreal3));
	n = 0;
	for (j=0; j<atoms->n_atom_total; j++) {
	  if (!atoms->filter[j]) continue;
	  property_real(atoms, i, 0, n) = tmpreal3[3*j+0];
	  property_real(atoms, i, 1, n) = tmpreal3[3*j+1];
	  property_real(atoms, i, 2, n) = tmpreal3[3*j+2];
	  n++;
	}
      }
    } else if (atoms->property_type[i] == PROPERTY_INT) {
      if (atoms->property_ncols[i] == 1) {
	start2[0] = frame;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = atoms->n_atom_total;
	//if (atomlist) {
	  netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start2, count2, tmpint));
	  n = 0;
	  for (j=0; j<atoms->n_atom_total; j++) {
	    if (!atoms->filter[j]) continue;
	    //atoms->int_data[atoms->property_start[i]*atoms->n_atom + n] = tmpint[j];
	    property_int(atoms,i,0,n) = tmpint[j];
	    n++;
	  }
	  //} else 
	  //netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start2, count2,
	  //			      &atoms->int_data[(atoms->property_start[i])*atoms->n_atom]));
	  //netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start2, count2,
	  //			       &property_int(atoms,i,0,0)));
      } else {
/* 	for (k=0; k<3; k++) { */
/* 	  start3[0] = frame; */
/* 	  start3[1] = 0; */
/* 	  start3[2] = k; */
/* 	  count3[0] = 1; */
/* 	  count3[1] = atoms->n_atom_total; */
/* 	  count3[2] = 1; */
/* 	  //if (atomlist) { */
/* 	    netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3, tmpint)); */
/* 	    n = 0; */
/* 	    for (j=0; j<atoms->n_atom_total; j++) { */
/* 	      if (!atoms->filter[j]) continue; */
/* 	      //atoms->int_data[(atoms->property_start[i]+k)*atoms->n_atom + n] = tmpint[j]; */
/* 	      property_int(atoms,i,k,n) = tmpint[j]; */
/* 	      n++; */
/* 	    } */
/* 	    //} */
/* 	    //else */
/* 	    //netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3,  */
/* 	    //				 &atoms->int_data[(atoms->property_start[i] + k)*atoms->n_atom])); */
/* 	    //netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3,  */
/* 	    //				 &property_int(atoms,i,k,0))); */
/* 	} */
	start3[0] = frame;
	start3[1] = 0;
	start3[2] = 0;
	count3[0] = 1;
	count3[1] = atoms->n_atom_total;
	count3[2] = 3;
	netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3, tmpint3));
	n = 0;
	for (j=0; j<atoms->n_atom_total; j++) {
	  if (!atoms->filter[j]) continue;
	  property_int(atoms, i, 0, n) = tmpint3[3*j+0];
	  property_int(atoms, i, 1, n) = tmpint3[3*j+1];
	  property_int(atoms, i, 2, n) = tmpint3[3*j+2];
	  n++;
	}

      }
    }
    else if (atoms->property_type[i] == PROPERTY_LOGICAL) {
      if (atoms->property_ncols[i] == 1) {
	start2[0] = frame;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = atoms->n_atom_total;
	//if (atomlist) {
	  netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start2, count2, tmpint));
	  n = 0;
	  for (j=0; j<atoms->n_atom_total; j++) {
	    if (!atoms->filter[j]) continue;
	    //atoms->logical_data[atoms->property_start[i]*atoms->n_atom + n] = tmpint[j];
	    property_logical(atoms,i,0,n) = tmpint[j];
	    n++;
	  }
	  //} else 
	  //netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start2, count2,
	  //		      &atoms->logical_data[(atoms->property_start[i])*atoms->n_atom]));
	  //netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start2, count2,
	  //			       &property_logical(atoms,i,0,0)));

      } else {
/* 	for (k=0; k<3; k++) { */
/* 	  start3[0] = frame; */
/* 	  start3[1] = 0; */
/* 	  start3[2] = k; */
/* 	  count3[0] = 1; */
/* 	  count3[1] = atoms->n_atom_total; */
/* 	  count3[2] = 1; */
/* 	  //if (atomlist) { */
/* 	    netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3, tmpint)); */
/* 	    n = 0; */
/* 	    for (j=0; j<atoms->n_atom_total; j++) { */
/* 	      if (!atoms->filter[j]) continue; */
/* 	      //atoms->logical_data[(atoms->property_start[i]+k)*atoms->n_atom + n] = tmpint[j]; */
/* 	      property_logical(atoms,i,k,n) = tmpint[j]; */
/* 	      n++; */
/* 	    } */
/* 	    //} */
/* 	//else */
/* 	    //netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3,  */
/* 	    //				 &atoms->logical_data[(atoms->property_start[i] + k)*atoms->n_atom])); */
/* 	    //netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3,  */
/* 	//				 &property_logical(atoms,i,k,0))); */
	start3[0] = frame;
	start3[1] = 0;
	start3[2] = 0;
	count3[0] = 1;
	count3[1] = atoms->n_atom_total;
	count3[2] = 3;
	netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3, tmpint3));
	n = 0;
	for (j=0; j<atoms->n_atom_total; j++) {
	  if (!atoms->filter[j]) continue;
	  property_logical(atoms, i, 0, n) = tmpint3[3*j+0];
	  property_logical(atoms, i, 1, n) = tmpint3[3*j+1];
	  property_logical(atoms, i, 2, n) = tmpint3[3*j+2];
	  n++;
	}
      }
    }
    else if (atoms->property_type[i] == PROPERTY_STR) {
	start3[0] = frame;
	start3[1] = 0;
	start3[2] = 0;
	count3[0] = 1;
	count3[1] = atoms->n_atom_total;
	count3[2] = PROPERTY_STRING_LENGTH;
	//if (atomlist) {
	  netcdf_check(nc_get_vara_text(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3, tmpchar));
	  n = 0;
	  for (j=0; j<atoms->n_atom_total; j++) {
	    if (!atoms->filter[j]) continue;
	    //strncpy(&atoms->str_data[PROPERTY_STRING_LENGTH*(atoms->property_start[i]*atoms->n_atom + n)], 
	    //	    &tmpchar[PROPERTY_STRING_LENGTH*j], PROPERTY_STRING_LENGTH);
	    strncpy(&property_str(atoms,i,0,n), &tmpchar[PROPERTY_STRING_LENGTH*j], PROPERTY_STRING_LENGTH);
	    n++;
	  }
	  //} else 
	  //netcdf_check(nc_get_vara_text(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3,
	  //				&(atoms->str_data[PROPERTY_STRING_LENGTH*atoms->property_start[i]*atoms->n_atom])));
	  //netcdf_check(nc_get_vara_text(nc_id, atoms->property_var_id[i][NETCDF_IN], start3, count3,
	  //				&property_str(atoms,i,0,0)));

    }
  }

  if (replacefill) {
    debug("replacing fill values with %d,%f\n",irep,rrep);
    replace_fill_values(atoms, irep, rrep);
  }
  
  debug("done reading netcdf\n");

  //if (atomlist) {
    free(tmpint);
    free(tmpint3);
    free(tmpreal);
    free(tmpreal3);
    free(tmpchar);
    //}

  return atoms->n_frame;
} 

int write_netcdf(int nc_id, Atoms *atoms, int frame, int redefine, 
		 int shuffle, int deflate, int deflate_level)
{
  int retval;
  int i, j, ndims;
  char pname[MAX_PROP_LENGTH];
  size_t start1[1], count1[1], start2[2], count2[2], start3[3], count3[3];
  int dims1[1], dims2[2], dims3[3];
  double cell_lengths[3], cell_angles[3];
  char varname[NC_MAX_NAME+1];
  nc_type vartype;
  int dimids[NC_MAX_VAR_DIMS], natts, newfile;
  int *tmpint, *tmpint3;
  double *tmpreal, *tmpreal3;
  char *tmpchar;
  int nvars, ngatts, unlimdimid, newvar;

  netcdf_check(nc_inq(nc_id, &ndims, &nvars, &ngatts, &unlimdimid));
  newfile = ndims == 0;

  if (newfile) { // it's a new file
    // set dimension and variable information

    // Global attributes
    netcdf_check(nc_put_att_text(nc_id, NC_GLOBAL, "Conventions", strlen("AMBER"), "AMBER"));
    netcdf_check(nc_put_att_text(nc_id, NC_GLOBAL, "ConventionVersion", strlen("1.0"), "1.0"));
    netcdf_check(nc_put_att_text(nc_id, NC_GLOBAL, "application", strlen("libAtoms"), "libAtoms"));
    netcdf_check(nc_put_att_text(nc_id, NC_GLOBAL, "program", strlen("xyz2nc"), "xyz2nc"));
    netcdf_check(nc_put_att_text(nc_id, NC_GLOBAL, "programVersion", strlen(SVN_VERSION), SVN_VERSION));
    netcdf_check(nc_put_att_text(nc_id, NC_GLOBAL, "title", strlen("Atoms Object"), "Atoms Object"));
    
    // Dimensions
    netcdf_check(nc_def_dim(nc_id, "frame", NC_UNLIMITED, &atoms->frame_dim_id[NETCDF_OUT]));
    netcdf_check(nc_def_dim(nc_id, "spatial", 3, &atoms->spatial_dim_id[NETCDF_OUT]));
    netcdf_check(nc_def_dim(nc_id, "atom", atoms->n_atom, &atoms->atom_dim_id[NETCDF_OUT]));
    netcdf_check(nc_def_dim(nc_id, "cell_spatial", 3, &atoms->cell_spatial_dim_id[NETCDF_OUT]));
    netcdf_check(nc_def_dim(nc_id, "cell_angular", 3, &atoms->cell_angular_dim_id[NETCDF_OUT]));
    netcdf_check(nc_def_dim(nc_id, "label", PROPERTY_STRING_LENGTH, &atoms->label_dim_id[NETCDF_OUT]));
    netcdf_check(nc_def_dim(nc_id, "string", PARAM_STRING_LENGTH, &atoms->string_dim_id[NETCDF_OUT]));
      
    // Label variables
    dims1[0] = atoms->spatial_dim_id[NETCDF_OUT];
    netcdf_check(nc_def_var(nc_id, "spatial", NC_CHAR, 1, dims1, &atoms->spatial_var_id[NETCDF_OUT]));

    dims1[0] = atoms->cell_spatial_dim_id[NETCDF_OUT];
    netcdf_check(nc_def_var(nc_id, "cell_spatial", NC_CHAR, 1, dims1, &atoms->cell_spatial_var_id[NETCDF_OUT]));
      
    dims2[0] = atoms->cell_angular_dim_id[NETCDF_OUT];
    dims2[1] = atoms->label_dim_id[NETCDF_OUT];
    netcdf_check(nc_def_var(nc_id, "cell_angular", NC_CHAR, 2, dims2, &atoms->cell_angular_var_id[NETCDF_OUT]));
    
    dims2[0] = atoms->frame_dim_id[NETCDF_OUT];
    dims2[1] = atoms->spatial_dim_id[NETCDF_OUT];
    netcdf_check(nc_def_var(nc_id, "cell_lengths", NC_DOUBLE, 2, dims2, &atoms->cell_lengths_var_id[NETCDF_OUT]));
    netcdf_check(nc_def_var(nc_id, "cell_angles", NC_DOUBLE, 2, dims2, &atoms->cell_angles_var_id[NETCDF_OUT]));
    

#ifdef NETCDF4
    if (atoms->netcdf4) {
      netcdf_check(nc_def_var_deflate(nc_id, atoms->spatial_var_id[NETCDF_OUT], shuffle, deflate, deflate_level));
      netcdf_check(nc_def_var_deflate(nc_id, atoms->cell_spatial_var_id[NETCDF_OUT], shuffle, deflate, deflate_level));
      netcdf_check(nc_def_var_deflate(nc_id, atoms->cell_angular_var_id[NETCDF_OUT], shuffle, deflate, deflate_level));
      netcdf_check(nc_def_var_deflate(nc_id, atoms->cell_lengths_var_id[NETCDF_OUT], shuffle, deflate, deflate_level));
      netcdf_check(nc_def_var_deflate(nc_id, atoms->cell_angles_var_id[NETCDF_OUT], shuffle, deflate, deflate_level));
    }
#endif
  }

  if (newfile || redefine) {

    nc_redef(nc_id);
    
    // Define variables for per-frame parameters
    for (i=0; i<atoms->n_param; i++) {
      if (strcasecmp(atoms->param_key[i], "Lattice") == 0 ||
	  strcasecmp(atoms->param_key[i], "Properties") == 0) continue;
      if (!atoms->param_filter[i]) continue;

      newvar = 0;
      switch(atoms->param_type[i]) {
      case(T_INTEGER):
	dims1[0] = atoms->frame_dim_id[NETCDF_OUT];
	if (nc_inq_varid(nc_id, atoms->param_key[i], &atoms->param_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	  netcdf_check(nc_def_var(nc_id, atoms->param_key[i], NC_INT, 1, dims1, &atoms->param_var_id[i][NETCDF_OUT]));
	  newvar = 1;
	} else {
	  netcdf_check(nc_inq_var(nc_id, atoms->param_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	  if (ndims != 1) {
	    fprintf(stderr,"Mismatch in parameter %s - ndims(%d) != 1\n", atoms->param_key[i], ndims);
	    return 0;
	  }
	  if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT]) {
	    fprintf(stderr, "Mismatch in parameter %s - dimension 0 != frame\n",atoms->param_key[i]);
	    return 0;
	  }
	}
	break;
      case(T_REAL):
	dims1[0] = atoms->frame_dim_id[NETCDF_OUT];
	if (nc_inq_varid(nc_id, atoms->param_key[i], &atoms->param_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	  netcdf_check(nc_def_var(nc_id, atoms->param_key[i], NC_DOUBLE, 1, dims1, &atoms->param_var_id[i][NETCDF_OUT]));
	  newvar = 1;
	} else {
	  netcdf_check(nc_inq_var(nc_id, atoms->param_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	  if (ndims != 1) {
	    fprintf(stderr,"Mismatch in parameter %s - ndims(%d) != 1\n", atoms->param_key[i], ndims);
	    return 0;
	  }
	  if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT]) {
	    fprintf(stderr,"Mismatch in parameter %s - dimension 0 != frame\n",atoms->param_key[i]);
	    return 0;
	  }
	}
	break;
      case(T_CHAR):
	dims2[0] = atoms->frame_dim_id[NETCDF_OUT];
	dims2[1] = atoms->string_dim_id[NETCDF_OUT];
	if (nc_inq_varid(nc_id, atoms->param_key[i], &atoms->param_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	  netcdf_check(nc_def_var(nc_id, atoms->param_key[i], NC_CHAR, 2, dims2, &atoms->param_var_id[i][NETCDF_OUT]));
	  newvar = 1;
	} else {
	  netcdf_check(nc_inq_var(nc_id, atoms->param_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	  if (ndims != 2) {
	    fprintf(stderr,"Mismatch in parameter %s - ndims(%d) != 1\n", atoms->param_key[i], ndims);
	    return 0;
	  }
	  if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->string_dim_id[NETCDF_OUT]) {
	    fprintf(stderr,"Mismatch in parameter %s - dimension 0 != frame or dimension 1 != string\n",atoms->param_key[i]);
	    return 0;
	  }
	}
	break;
      case(T_INTEGER_A):
	dims2[0] = atoms->frame_dim_id[NETCDF_OUT];
	dims2[1] = atoms->spatial_dim_id[NETCDF_OUT];
	if (nc_inq_varid(nc_id, atoms->param_key[i], &atoms->param_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	  netcdf_check(nc_def_var(nc_id, atoms->param_key[i], NC_INT, 2, dims2, &atoms->param_var_id[i][NETCDF_OUT]));
	  newvar = 1;
	} else {
	  netcdf_check(nc_inq_var(nc_id, atoms->param_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	  if (ndims != 2) {
	    fprintf(stderr,"Mismatch in parameter %s - ndims(%d) != 1\n", atoms->param_key[i], ndims);
	    return 0;
	  }
	  if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->spatial_dim_id[NETCDF_OUT])  {
	    fprintf(stderr,"Mismatch in parameter %s - dimension 0 != frame or dimension 1 != spatial\n", atoms->param_key[i]);
	    return 0;
	  }
	}	
	break;
      case(T_REAL_A):
	dims2[0] = atoms->frame_dim_id[NETCDF_OUT];
	dims2[1] = atoms->spatial_dim_id[NETCDF_OUT];
	if (nc_inq_varid(nc_id, atoms->param_key[i], &atoms->param_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	  netcdf_check(nc_def_var(nc_id, atoms->param_key[i], NC_DOUBLE, 2, dims2, &atoms->param_var_id[i][NETCDF_OUT]));
	  newvar = 1;
	} else {
	  netcdf_check(nc_inq_var(nc_id, atoms->param_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	  if (ndims != 2) {
	    fprintf(stderr,"Mismatch in parameter %s - ndims(%d) != 1\n", atoms->param_key[i], ndims);
	    return 0;
	  }
	  if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->spatial_dim_id[NETCDF_OUT])  {
	    fprintf(stderr,"Mismatch in parameter %s - dimension 0 != frame or dimension 1 != spatial\n", atoms->param_key[i]);
	    return 0;
	  }
	}	
	break;
      default:
	fprintf(stderr,"Unknown parameter %s type %d when writing netcdf file\n", atoms->param_key[i],
		atoms->param_type[i]);
	return 0;
      }

      nc_put_att_int(nc_id, atoms->param_var_id[i][NETCDF_OUT], "type", NC_INT, 1, &(atoms->param_type[i]));
#ifdef NETCDF4
      if (atoms->netcdf4 && (newfile || newvar))
	netcdf_check(nc_def_var_deflate(nc_id, atoms->param_var_id[i][NETCDF_OUT], shuffle, deflate, deflate_level));
#endif
    }

    // Define variables for per-atom properties
    for (i=0; i<atoms->n_property; i++) {
      if (!atoms->property_filter[i]) continue;
      ndims = atoms->property_ncols[i];
      if (ndims != 1 && ndims!= 3) {
	fprintf(stderr,"Property %s ncols(%d) != 1 or 3\n", atoms->property_name[i], ndims);
	return 0;
      }

      newvar = 0;
      switch(atoms->property_type[i]) {
      case(PROPERTY_INT):
	if (ndims == 1) {
	  dims2[0] = atoms->frame_dim_id[NETCDF_OUT];
	  dims2[1] = atoms->atom_dim_id[NETCDF_OUT];
	  if (nc_inq_varid(nc_id, atoms->property_name[i], &atoms->property_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	    netcdf_check(nc_def_var(nc_id, atoms->property_name[i], NC_INT, 2, dims2, &atoms->property_var_id[i][NETCDF_OUT]));
	    newvar = 1;
	  } else {
	    netcdf_check(nc_inq_var(nc_id, atoms->property_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	    if (ndims != 2) {
	      fprintf(stderr,"Mismatch in property %s - ndims(%d) != 1\n", atoms->property_name[i], ndims);
	      return 0;
	    }
	    if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->atom_dim_id[NETCDF_OUT]) {
	      fprintf(stderr,"Mismatch in property %s - dimension 0 != frame or dimension 1 != atom\n", atoms->property_name[i]);
	      return 0;
	    }
	  }	
	} else {
	  dims3[0] = atoms->frame_dim_id[NETCDF_OUT];
	  dims3[1] = atoms->atom_dim_id[NETCDF_OUT];
	  dims3[2] = atoms->spatial_dim_id[NETCDF_OUT];
	  if (nc_inq_varid(nc_id, atoms->property_name[i], &atoms->property_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	    netcdf_check(nc_def_var(nc_id, atoms->property_name[i], NC_INT, 3, dims3, &atoms->property_var_id[i][NETCDF_OUT]));
	    newvar = 1;
	  } else {
	    netcdf_check(nc_inq_var(nc_id, atoms->property_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	    if (ndims != 3) {
	      fprintf(stderr,"Mismatch in property %s - ndims(%d) != 1\n", atoms->property_name[i], ndims);
	      return 0;
	    }
	    if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->atom_dim_id[NETCDF_OUT] || dimids[2] != atoms->spatial_dim_id[NETCDF_OUT]) {
	      fprintf(stderr,"Mismatch in property %s - dimension 0 != frame or dimension 1 != atom or dimension 2 != spatial\n", 
		      atoms->property_name[i]);
	      return 0;
	    }
	  }	
	}
	break;

      case(PROPERTY_REAL):
	
	// Name mangling for compatibility with AMBER/VMD
	if (strcmp(atoms->property_name[i], "pos") == 0) 
	  strcpy(pname, "coordinates");
	else if(strcmp(atoms->property_name[i], "velo") == 0)
	  strcpy(pname, "velocities");
	else
	  strcpy(pname, atoms->property_name[i]);
	
	if (ndims == 1) {
	  dims2[0] = atoms->frame_dim_id[NETCDF_OUT];
	  dims2[1] = atoms->atom_dim_id[NETCDF_OUT];
	  if (nc_inq_varid(nc_id, pname, &atoms->property_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	    netcdf_check(nc_def_var(nc_id, pname, NC_DOUBLE, 2, dims2, &atoms->property_var_id[i][NETCDF_OUT]));
	    newvar = 1;
	  } else {
	    netcdf_check(nc_inq_var(nc_id, atoms->property_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	    if (ndims != 2) {
	      fprintf(stderr,"Mismatch in property %s - ndims(%d) != 1\n", pname, ndims);
	      return 0;
	    }
	    if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->atom_dim_id[NETCDF_OUT])  {
	      pe("Mismatch in property %s - dimension 0 != frame or dimension 1 != atom\n", atoms->property_name[i]);
	      return 0;
	    }
	  }		  
	} else {
	  dims3[0] = atoms->frame_dim_id[NETCDF_OUT];
	  dims3[1] = atoms->atom_dim_id[NETCDF_OUT];
	  dims3[2] = atoms->spatial_dim_id[NETCDF_OUT];
	  if (nc_inq_varid(nc_id, pname, &atoms->property_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	    netcdf_check(nc_def_var(nc_id, pname, NC_DOUBLE, 3, dims3, &atoms->property_var_id[i][NETCDF_OUT]));
	    newvar = 1;
	  } else {
	    netcdf_check(nc_inq_var(nc_id, atoms->property_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	    if (ndims != 3) {
	      fprintf(stderr,"Mismatch in property %s - ndims(%d) != 1\n", atoms->property_name[i], ndims);
	      return 0;
	    }
	    if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->atom_dim_id[NETCDF_OUT] || 
		dimids[2] != atoms->spatial_dim_id[NETCDF_OUT]) {
	      pe("Mismatch in property %s - dimension 0 != frame or dimension 1 != atom or dimension 2 != spatial\n", 
		      pname);
	      return 0;
	    }
	  }	
	}
	break;

      case(PROPERTY_STR):
	if (ndims == 1) {
	  dims3[0] = atoms->frame_dim_id[NETCDF_OUT];
	  dims3[1] = atoms->atom_dim_id[NETCDF_OUT];
	  dims3[2] = atoms->label_dim_id[NETCDF_OUT];
	  if (nc_inq_varid(nc_id, atoms->property_name[i], &atoms->property_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	    netcdf_check(nc_def_var(nc_id, atoms->property_name[i], NC_CHAR, 3, dims3, &atoms->property_var_id[i][NETCDF_OUT]));
	    newvar = 1;
	  } else {
	    netcdf_check(nc_inq_var(nc_id, atoms->property_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	    if (ndims != 3) {
	      fprintf(stderr,"Mismatch in property %s - ndims(%d) != 1\n", atoms->property_name[i], ndims);
	      return 0;
	    }
	    if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->atom_dim_id[NETCDF_OUT] || dimids[2] != atoms->label_dim_id[NETCDF_OUT])  {
	      pe("Mismatch in property %s - dimension 0 != frame or dimension 1 != atom or dimension 2 != label\n", 
		      atoms->property_name[i]);
	      return 0;
	    }
	  }	
	}
	else {
	  fprintf(stderr,"String property %s has ndims == 3, not supported\n", atoms->property_name[i]);
	  return 0;
	}
	break;

      case(PROPERTY_LOGICAL):
	if (ndims == 1) {
	  dims2[0] = atoms->frame_dim_id[NETCDF_OUT];
	  dims2[1] = atoms->atom_dim_id[NETCDF_OUT];
	  if (nc_inq_varid(nc_id, atoms->property_name[i], &atoms->property_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	    netcdf_check(nc_def_var(nc_id, atoms->property_name[i], NC_INT, 2, dims2, &atoms->property_var_id[i][NETCDF_OUT]));
	    newvar = 1;
	  } else {

	    netcdf_check(nc_inq_var(nc_id, atoms->property_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	    if (ndims != 2) {
	      pe("Mismatch in property %s - ndims(%d) != 1\n", atoms->property_name[i], ndims);
	      return 0;
	    }
	    if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->atom_dim_id[NETCDF_OUT]) {
	      pe("Mismatch in property %s - dimension 0 != frame or dimension 1 != atom\n", atoms->property_name[i]);
	      return 0;
	    }
	  }	
	} else {
	  dims3[0] = atoms->frame_dim_id[NETCDF_OUT];
	  dims3[1] = atoms->atom_dim_id[NETCDF_OUT];
	  dims3[2] = atoms->spatial_dim_id[NETCDF_OUT];
	  if (nc_inq_varid(nc_id, atoms->property_name[i], &atoms->property_var_id[i][NETCDF_OUT]) != NC_NOERR) {
	    netcdf_check(nc_def_var(nc_id, atoms->property_name[i], NC_INT, 3, dims3, &atoms->property_var_id[i][NETCDF_OUT]));
	  } else {
	    netcdf_check(nc_inq_var(nc_id, atoms->property_var_id[i][NETCDF_OUT], varname, &vartype, &ndims, dimids, &natts));
	    if (ndims != 3) {
	      pe("Mismatch in property %s - ndims(%d) != 1\n", atoms->property_name[i], ndims);
	      return 0;
	    }
	    if (dimids[0] != atoms->frame_dim_id[NETCDF_OUT] || dimids[1] != atoms->atom_dim_id[NETCDF_OUT] || dimids[2] != atoms->spatial_dim_id[NETCDF_OUT])  {
	      pe("Mismatch in property %s - dimension 0 != frame or dimension 1 != atom or dimension 2 != spatial\n", 
		      atoms->property_name[i]);
	      return 0;
	    }
	  }	

	}

	break;

      default:
	  fprintf(stderr,"Unknown property type %d for property %s\n", atoms->property_type[i], 
		  pname);
	  return 0;
      }
      nc_put_att_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], "type", NC_INT, 1, &(atoms->property_type[i]));
#ifdef NETCDF4
      if (atoms->netcdf4 && (newfile || newvar))
	netcdf_check(nc_def_var_deflate(nc_id, atoms->property_var_id[i][NETCDF_OUT], shuffle, deflate, deflate_level));
#endif
    }
      
    netcdf_check(nc_enddef(nc_id));

    if (newfile) {
      // Put data for label variables
      netcdf_check(nc_put_var_text(nc_id, atoms->spatial_var_id[NETCDF_OUT], "xyz"));
      netcdf_check(nc_put_var_text(nc_id, atoms->cell_spatial_var_id[NETCDF_OUT], "abc"));
    
      start2[0] = 0;
      start2[1] = 0;
      count2[0] = 1;
      count2[1] = strlen("alpha");
      netcdf_check(nc_put_vara_text(nc_id, atoms->cell_angular_var_id[NETCDF_OUT], start2, count2, "alpha"));
      start2[0] = 1;
      start2[1] = 0;
      count2[0] = 1;
      count2[1] = strlen("beta");
      netcdf_check(nc_put_vara_text(nc_id, atoms->cell_angular_var_id[NETCDF_OUT], start2, count2, "beta"));
      start2[0] = 2;
      start2[1] = 0;
      count2[0] = 1;
      count2[1] = strlen("gamma");
      netcdf_check(nc_put_vara_text(nc_id, atoms->cell_angular_var_id[NETCDF_OUT], start2, count2, "gamma"));
    }

  } /* end if (newfile || redefine) */

  tmpint = malloc(sizeof(int)*atoms->n_atom);
  tmpint3 = malloc(sizeof(int)*atoms->n_atom*3);
  tmpreal = malloc(sizeof(double)*atoms->n_atom);
  tmpreal3 = malloc(sizeof(double)*atoms->n_atom*3);
  tmpchar = malloc(sizeof(char)*atoms->n_atom*PROPERTY_STRING_LENGTH);

  if (tmpint == NULL || tmpreal == NULL || tmpchar == NULL)  {
    fprintf(stderr,"Error allocating memory in write_netcdf\n");
    return 0;
  }

  debug("writing netcdf frame=%d\n", frame);

  // Put lattice
  lattice_xyz_to_abc(atoms->lattice, cell_lengths, cell_angles);
  start2[0] = frame;
  start2[1] = 0;
  count2[0] = 1;
  count2[1] = 3;
  netcdf_check(nc_put_vara_double(nc_id, atoms->cell_lengths_var_id[NETCDF_OUT], start2, count2, cell_lengths));
  netcdf_check(nc_put_vara_double(nc_id, atoms->cell_angles_var_id[NETCDF_OUT], start2, count2, cell_angles));

  // Write the data. First the parameters.
  debug("writing params\n");
  for (i=0; i<atoms->n_param; i++) {
      if (strcasecmp(atoms->param_key[i], "Lattice") == 0 ||
	  strcasecmp(atoms->param_key[i], "Properties") == 0) continue;
      if (!atoms->param_filter[i]) continue;

      debug("  writing %s\n", atoms->param_key[i]);

      switch(atoms->param_type[i]) {
      case(T_INTEGER):
	start1[0] = frame;
	count1[0] = 1;
	netcdf_check(nc_put_vara_int(nc_id, atoms->param_var_id[i][NETCDF_OUT], start1, count1, &atoms->param_int[i]));
	break;
      case(T_REAL):
	start1[0] = frame;
	count1[0] = 1;
	netcdf_check(nc_put_vara_double(nc_id, atoms->param_var_id[i][NETCDF_OUT], start1, count1, &atoms->param_real[i]));
	break;
      case(T_CHAR):
	start2[0] = frame;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = PARAM_STRING_LENGTH;
	netcdf_check(nc_put_vara_text(nc_id, atoms->param_var_id[i][NETCDF_OUT], start2, count2, atoms->param_value[i]));
	break;
      case(T_INTEGER_A):
	start2[0] = frame;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = 3;
	netcdf_check(nc_put_vara_int(nc_id, atoms->param_var_id[i][NETCDF_OUT], start2, count2, &atoms->param_int_a[i][0]));
	break;
      case(T_REAL_A):
	start2[0] = frame;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = 3;
	netcdf_check(nc_put_vara_double(nc_id, atoms->param_var_id[i][NETCDF_OUT], start2, count2, &atoms->param_real_a[i][0]));
	break;
      default:
	  fprintf(stderr,"Unknown parameter %s type %d when writing netcdf file\n", atoms->param_key[i],
		  atoms->param_type[i]);
	  return 0;
      }
  }
  
  debug("writing properties\n");
  // and now the properties
  for (i=0; i<atoms->n_property; i++) {
    debug("  writing %s\n", atoms->property_name[i]);
    if (!atoms->property_filter[i]) continue;
    ndims = atoms->property_ncols[i];
    switch(atoms->property_type[i]) {
    case(PROPERTY_INT):
      if (ndims == 1) {
	start2[0] = frame;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = atoms->n_atom;
	//netcdf_check(nc_put_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], start2, count2, 
	//			     &atoms->int_data[(atoms->property_start[i])*atoms->n_atom]));	
	for (j=0; j<atoms->n_atom; j++)
	  tmpint[j] = property_int(atoms, i, 0, j);
	netcdf_check(nc_put_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], start2, count2, 
				     tmpint));

      } else {
/* 	for (k=0; k<3; k++) { */
/* 	  start3[0] = frame; */
/* 	  start3[1] = 0; */
/* 	  start3[2] = k; */
/* 	  count3[0] = 1; */
/* 	  count3[1] = atoms->n_atom; */
/* 	  count3[2] = 1; */
/* 	  //netcdf_check(nc_put_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], start3, count3,  */
/* 	  //			       &atoms->int_data[(atoms->property_start[i] + k)*atoms->n_atom])); */
/* 	  for (j=0; j<atoms->n_atom; j++) */
/* 	    tmpint[j] = property_int(atoms, i, k, j); */
/* 	  netcdf_check(nc_put_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], start3, count3,  */
/* 				       tmpint)); */

/* 	} */
	start3[0] = frame;
	start3[1] = 0;
	start3[2] = 0;
	count3[0] = 1;
	count3[1] = atoms->n_atom;
	count3[2] = 3;
	for (j=0; j<atoms->n_atom; j++) {
	  tmpint3[3*j+0] = property_int(atoms, i, 0, j);
	  tmpint3[3*j+1] = property_int(atoms, i, 1, j);
	  tmpint3[3*j+2] = property_int(atoms, i, 2, j);
	}
	netcdf_check(nc_put_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], 
				     start3, count3, tmpint3));
      }
      break;

    case(PROPERTY_REAL):
      if (ndims == 1) {
	start2[0] = frame;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = atoms->n_atom;
	//netcdf_check(nc_put_vara_double(nc_id, atoms->property_var_id[i][NETCDF_OUT], start2, count2, 
	//				&atoms->real_data[atoms->property_start[i]*atoms->n_atom]));
	for (j=0; j<atoms->n_atom; j++)
	  tmpreal[j] = property_real(atoms, i, 0, j);
	netcdf_check(nc_put_vara_double(nc_id, atoms->property_var_id[i][NETCDF_OUT], start2, count2, 
					tmpreal));

      } else {
/* 	for (k=0; k<3; k++) { */
/* 	  start3[0] = frame; */
/* 	  start3[1] = 0; */
/* 	  start3[2] = k; */
/* 	  count3[0] = 1; */
/* 	  count3[1] = atoms->n_atom; */
/* 	  count3[2] = 1; */
/* 	  //netcdf_check(nc_put_vara_double(nc_id, atoms->property_var_id[i][NETCDF_OUT], start3, count3,  */
/* 	  //				  &atoms->real_data[(atoms->property_start[i] + k)*atoms->n_atom])); */
/* 	  for (j=0; j<atoms->n_atom; j++) */
/* 	    tmpreal[j] = property_real(atoms, i, k, j); */
/* 	  netcdf_check(nc_put_vara_double(nc_id, atoms->property_var_id[i][NETCDF_OUT], start3, count3,  */
/* 					  tmpreal)); */

/* 	} */
	start3[0] = frame;
	start3[1] = 0;
	start3[2] = 0;
	count3[0] = 1;
	count3[1] = atoms->n_atom;
	count3[2] = 3;
	for (j=0; j<atoms->n_atom; j++) {
	  tmpreal3[3*j+0] = property_real(atoms, i, 0, j);
	  tmpreal3[3*j+1] = property_real(atoms, i, 1, j);
	  tmpreal3[3*j+2] = property_real(atoms, i, 2, j);
	}
	netcdf_check(nc_put_vara_double(nc_id, atoms->property_var_id[i][NETCDF_OUT], 
				     start3, count3, tmpreal3));
      }
      break;

    case(PROPERTY_STR):
      start3[0] = frame;
      start3[1] = 0;
      start3[2] = 0;
      count3[0] = 1;
      count3[1] = atoms->n_atom;
      count3[2] = PROPERTY_STRING_LENGTH;
      //   netcdf_check(nc_put_vara_text(nc_id, atoms->property_var_id[i][NETCDF_OUT], start3, count3,
      //			    &(atoms->str_data[PROPERTY_STRING_LENGTH*atoms->property_start[i]*atoms->n_atom])));
      for (j=0; j<atoms->n_atom; j++)
	strncpy(tmpchar+PROPERTY_STRING_LENGTH*j, &property_str(atoms,i,0,j), PROPERTY_STRING_LENGTH);
      netcdf_check(nc_put_vara_text(nc_id, atoms->property_var_id[i][NETCDF_OUT], start3, count3,
				    tmpchar));
      break;

    case(PROPERTY_LOGICAL):
      if (ndims == 1) {
	start2[0] = frame;
	start2[1] = 0;
	count2[0] = 1;
	count2[1] = atoms->n_atom;
	//netcdf_check(nc_put_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], start2, count2,
	//			     &atoms->logical_data[(atoms->property_start[i])*atoms->n_atom]));
	for (j=0; j<atoms->n_atom; j++)
	  tmpint[j] = property_logical(atoms, i, 0, j);
	netcdf_check(nc_put_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], start2, count2,
				     tmpint));
      } else {
/* 	for (k=0; k<3; k++) { */
/* 	  start3[0] = frame; */
/* 	  start3[1] = 0; */
/* 	  start3[2] = k; */
/* 	  count3[0] = 1; */
/* 	  count3[1] = atoms->n_atom; */
/* 	  count3[2] = 1; */
/* 	  //netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], start3, count3,  */
/* 	  //			       &atoms->logical_data[(atoms->property_start[i] + k)*atoms->n_atom])); */
/* 	  for (j=0; j<atoms->n_atom; j++) */
/* 	    tmpint[j] = property_logical(atoms, i, k, j); */
/* 	  netcdf_check(nc_get_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], start3, count3,  */
/* 				       tmpint)); */

/* 	} */
	start3[0] = frame;
	start3[1] = 0;
	start3[2] = 0;
	count3[0] = 1;
	count3[1] = atoms->n_atom;
	count3[2] = 3;
	for (j=0; j<atoms->n_atom; j++) {
	  tmpint3[3*j+0] = property_logical(atoms, i, 0, j);
	  tmpint3[3*j+1] = property_logical(atoms, i, 1, j);
	  tmpint3[3*j+2] = property_logical(atoms, i, 2, j);
	}
	netcdf_check(nc_put_vara_int(nc_id, atoms->property_var_id[i][NETCDF_OUT], 
				     start3, count3, tmpint3));

      }
      break;

    default:
      fprintf(stderr,"Unknown property type %d for property %s\n", atoms->property_type[i], 
	      atoms->property_name[i]);
      return 0;
    }
  }

  free(tmpint);
  free(tmpint3);
  free(tmpreal);
  free(tmpreal3);
  free(tmpchar);

  debug("done writing netcdf\n");
  netcdf_check(nc_sync(nc_id));
  atoms->n_frame++;
  return 1;
}

#else
int read_netcdf (int ncid, Atoms *atoms, int frame, int *atomlist, int natomlist, int query, 
		 int redefine, int realloc, int replacefill, int irep, double rrep)
{
  pe("No NetCDF support compiled in.\n");
  return 0;
}

int write_netcdf(int ncid, Atoms *atoms, int frame, int redefine,
                 int shuffle, int deflate, int deflate_level)
{
  pe("No NetCDF support compiled in.\n");
  return 0;
}
#endif


/* Find starting positions of xyz frames within a file
 * Uses a disk cache to save recomputing if xyz
 * file hasn't been modified. Returns number of frames.
 *
 * File format: text
 * First line:  number of frames (int)
 * Subsequent nframes+1 lines: offset (long), natoms (int)
 * Last offset is end of final frame scanned.
 */

#include <libgen.h>

int xyz_find_frames(char *fname, long *frames, int *atoms) {
  FILE *in, *index;
  char *bname;
  char indexname[LINESIZE], linebuffer[LINESIZE], buf1[LINESIZE], buf2[LINESIZE];
  int natoms, i, nframes;
  int from_scratch, do_update; 
  struct stat xyz_stat, idx_stat;

  strcpy(indexname, fname);
  strcat(indexname, ".idx");

  if (stat(fname, &xyz_stat) != 0) {
    fprintf(stderr,"Cannot stat xyz file %s\n", fname);
    return 0;
  }
  
  from_scratch = stat(indexname, &idx_stat) != 0;

  if (from_scratch) {
    // Try to read from current dir instead
    strcpy(buf1, indexname);
    bname = basename(buf1);
    if (getcwd(buf2, LINESIZE) != NULL) {
      strcat(buf2, "/");
      strcat(buf2, bname);
      if (stat(buf2, &idx_stat) == 0) {
	fprintf(stderr,"Found index %s\n",buf2);
	strcpy(indexname,buf2);
	from_scratch = 0;
      }
    }
  }

  do_update = xyz_stat.st_mtime > idx_stat.st_mtime;
  
  if (!from_scratch) {
    debug("Reading XYZ index from file %s\n", indexname);
    index = fopen(indexname, "r");
    if (!fgets(linebuffer,LINESIZE,index)) {
      fprintf(stderr,"Index file %s is empty\n",indexname);
      return 0;
    }
    sscanf(linebuffer, "%d", &nframes);
    if (nframes+1 >= XYZ_MAX_FRAMES) {
      fprintf(stderr,"nframes(%d)+1 >= XYZ_MAX_FRAMES(%d)\n",nframes,XYZ_MAX_FRAMES);
      return 0;
    }
    for (i=0; i<=nframes; i++) {
      if (!fgets(linebuffer,LINESIZE,index)) {
	fprintf(stderr,"Premature end of indexfile %s\n",indexname);
	return 0;
      }
      sscanf(linebuffer, "%ld %d", &frames[i], &atoms[i]);
    }
    fclose(index);
  }

/*   for (i=0; i<nframes+1; i++) { */
/*     printf("%d: %ld\n", i, frames[i]); */
/*   } */

  if (from_scratch || do_update) {
    debug("Writing XYZ index to file %s\n", indexname);
    in = fopen(fname, "r");
    if (from_scratch) 
      nframes = 0;
    else {
      debug("Trying to update XYZ index... \n");

      // Try to seek past last frame, and check this looks
      // like the start of a new frame
      if (fseek(in,frames[nframes],SEEK_SET) != 0 ||
	  !fgets(linebuffer,LINESIZE,in) ||
	  sscanf(linebuffer, "%d", &natoms) != 1 ||
	  natoms != atoms[nframes]) {
	// Seek failed, we'll have to rebuild index from start
	fseek(in,0,SEEK_SET);
	nframes = 0;
	debug(" failed, rebuilding from scratch.\n");
      }
      else {
	// Rewind to start of frame
	fseek(in,frames[nframes],SEEK_SET);
      }

      // TODO - improve check - fails if number of atoms has changed
    }

    debug("starting to build index from file pos %ld nframes=%d\n", ftell(in), nframes);

    while (fgets(linebuffer,LINESIZE,in)) {
      frames[nframes] = ftell(in)-strlen(linebuffer);
      if (nframes+1 > XYZ_MAX_FRAMES) {
	fprintf(stderr,"nframes(%d)+1 > XYZ_MAX_FRAMES(%d)\n",nframes,XYZ_MAX_FRAMES);
	return 0;
      }
      if (sscanf(linebuffer, "%d", &natoms) != 1) {
	fprintf(stderr,"Malformed XYZ file %s at frame %d\n",fname,nframes);
	return 0;
      }

      atoms[nframes] = natoms;
      nframes++;

      // Skip the whole frame, as quickly as possible
      for (i=0; i<natoms+1; i++)
	if (!fgets(linebuffer,LINESIZE,in)) {
	  fseek(in, frames[nframes], SEEK_SET); // return file pointer to beginning of frame
	  nframes--; // incomplete last frame
	  goto XYZ_END;
	}
    }
  XYZ_END:
    frames[nframes] = ftell(in); // end of last frame in file
    atoms[nframes] = natoms;
    index = fopen(indexname, "w");
    if (index == NULL) {
      // Try to write in current dir instead
      strcpy(buf1, indexname);
      bname = basename(buf1);
      if (getcwd(buf2, LINESIZE) != NULL) {
	strcat(buf2, "/");
	strcat(buf2, bname);
	index = fopen(buf2, "w");
	debug("Writing index to %s\n", buf2);
      }
      if (index == NULL) {
	fprintf(stderr, "Cannot write index file.\n");
	fclose(in);
	return nframes;
      }
    } else
      debug("Writing index to %s\n", indexname);
    fprintf(index, "%d\n", nframes);
    for (i=0; i<=nframes; i++)
      fprintf(index, "%ld %d\n", frames[i], atoms[i]);
    fclose(in);
    fclose(index);
  }

/*   for (i=0; i<nframes+1; i++) { */
/*     printf("%d: %ld\n", i, frames[i]); */
/*   } */

  debug("xyz_find_frames %s: found %d complete frames\n", fname, nframes);

  return nframes;
}

int read_xyz (FILE *in, Atoms *atoms, int *atomlist, int natomlist, int frame, 
	      int query, int redefine, int realloc, int suppress, int override_lattice,
	      double lattice[3][3]) {
  int i,n, entry_count,j,k=0,ncols,m, atidx;
  char linebuffer[LINESIZE];
  char fields[MAX_ENTRY_COUNT][LINESIZE], subfields[MAX_ENTRY_COUNT][LINESIZE],
    finalfields[MAX_ENTRY_COUNT][LINESIZE];
  char *p, *p1, tmp_logical;
  int properties_idx, lattice_idx, nxyz, original_index=0, nfields=0, idx, offset, error;
  double tmpd;

  if (in != stdin && atoms->got_index)
    if (fseek(in, atoms->frames[frame], SEEK_SET) == -1) {
      fprintf(stderr,"cannot seek XYZ input file\n");
      return 0;
    }

  if (in != stdin || (in == stdin && query)) {
    if (!fgets(linebuffer,LINESIZE,in)) {
      if (!suppress) fprintf(stderr, "premature file ending - expecting number of atoms\n");
      return 0;
    }
    if (sscanf(linebuffer, "%d", &nxyz) != 1) {
      fprintf(stderr,"first line (%s) must be number of atoms\n", linebuffer);
      return 0;
    }

    // Read comment line, which should contain 'Lattice=' and 'Properties=' keys
    if (!fgets(linebuffer,LINESIZE,in)) {
      fprintf(stderr,"premature file ending - expecting comment line\n");
      return 0;
    }
    linebuffer[strlen(linebuffer)-1] = '\0';   // Remove trailing newline


    if (!strstr(linebuffer, "Lattice") && !strstr(linebuffer, "lattice")) {
      // It's not an extended XYZ file. Try to guess what's going on.
      // If comment line contains nine or more fields, assume last nine are
      // lattice in cartesian coordinates. 

      p = linebuffer;
      k = 0;
      while ((p1 = strsep(&p, " \t")) != NULL) {
	if (*p1 == '\0') continue;
	strcpy(fields[k++], p1);
      }
      
      if (k >= 9) {
	offset = k-9;
	error = 0;
	for (i=0; i<9; i++)
	  if (sscanf(fields[offset+i], "%lf", &tmpd) != 1) {
	    error = 1;
	    break;
	  }
	
	if ((p = strstr(linebuffer, "\n")) != NULL) *p = '\0';
	if (!error) {
	  sprintf(linebuffer, "%s Lattice=\"%s %s %s %s %s %s %s %s %s\"\n",
		  linebuffer, fields[offset+0], fields[offset+1], fields[offset+2], fields[offset+3],
		  fields[offset+4], fields[offset+5], fields[offset+6], fields[offset+7], fields[offset+8]);
	} else {
	  if (!override_lattice) {
	    fprintf(stderr,"Cannot extract lattice from line %s\n", linebuffer);
	    return 0;
	  }
	}
      } else {
	if (!override_lattice) {
	  fprintf(stderr,"Cannot extract lattice from line %s\n", linebuffer);
	  return 0;
	}
      }
    }
    
    if (!strstr(linebuffer, "Properties") && !strstr(linebuffer, "properties")) {
      // No Properties key. Add a default one.
      if ((p = strstr(linebuffer, "\n")) != NULL) *p = '\0';
      strcat(linebuffer, "Properties=species:S:1:pos:R:3\n");
    }

    // Parse parameters. First split on ", ', { or }
    p = linebuffer;
    k = 0;
    while ((p1 = strsep(&p, "\"'{}")) != NULL) {
      if (*p1 == '\0') continue;
      strcpy(fields[k++], p1);
    }
    
    // Now split things outside quotes on whitespace
    nfields = 0;
    for (i=0; i<k; i++) {
      if (i % 2 == 0) {
	p = fields[i];
	j = 0;
	while ((p1 = strsep(&p, " \t")) != NULL) {
	  if (*p1 == '\0') continue;
	  strcpy(subfields[j++], p1);
	}
	for (n=0; n<j; n++, nfields++) {
	  strcpy(finalfields[nfields],subfields[n]);
	}
	
      } else {
	strcat(finalfields[nfields-1],fields[i]);
      }
    }

    if (redefine || !atoms->initialised) {
      // Finally, split on '=' to get key/value pairs
      atoms->n_param = 0;
      for (i=0; i<nfields; i++) {
	strcpy(linebuffer, finalfields[i]);
	if ((p = strchr(linebuffer,'=')) == NULL) {
	  fprintf(stderr,"Badly formed key/value pair %s\n", linebuffer);
	  return 0;
	}
	
	*p = '\0';
	idx = atoms_add_param(atoms, linebuffer, 0, 0, -1, NETCDF_IN);
	strcpy(atoms->param_value[idx], p+1);
      }
      
      // Try to guess param types
      for (i=0; i<atoms->n_param; i++) {
	if (strcasecmp(atoms->param_key[i], "Lattice") == 0 ||
	    strcasecmp(atoms->param_key[i], "Properties") == 0) continue;
	
	strcpy(linebuffer, atoms->param_value[i]);
	k = 0;
	p = linebuffer;
	while ((p1 = strsep(&p, " ")) != NULL) {
	  if (*p1 == '\0') continue;
	  strcpy(fields[k++], p1);
	}
	if (k == 0) {
	  k = 1;
	  strcpy(fields[0], linebuffer);
	} 
	
	/*       if (k != 1 && k != 3) { */
	/* 	fprintf(stderr,"Parameter %s must have size 1 or 3\n", atoms->param_key[i]); */
	/* 	return 0; */
	/*       } */
	atoms->param_size[i] = k;
	
	for (j=0; j<k; j++) {
	  for (n=0; n<strlen(fields[j]); n++)
	    if (!isblank(fields[j][n]) && !isdigit(fields[j][n])) goto NOT_INT;
	}
	
	atoms->param_type[i] = (k == 1 ? T_INTEGER : T_INTEGER_A);
	continue;

      NOT_INT:
	for (j=0; j<k; j++)
	  if (strtod(fields[j], &p), strlen(p) != 0) goto NOT_REAL;
	
	atoms->param_type[i] = (k == 1 ? T_REAL : T_REAL_A);
	continue;
      
      NOT_REAL:
	// Fallback option: treat as a single string
	atoms->param_type[i] = T_CHAR;
	atoms->param_size[i] = 1;
	continue;
      }
      
      for (i=0; i<atoms->n_param; i++) {
	if ((atoms->param_type[i] == T_INTEGER_A || atoms->param_type[i] == T_REAL_A) &&
	    atoms->param_size[i] != 3) {
	  fprintf(stderr,"Parameter %s must have size 1 or 3, but got %d\n", atoms->param_key[i], atoms->param_size[i]);
	  return 0; 
	}
      }

      properties_idx = atoms_find_param(atoms, "Properties");
      if (properties_idx == -1) {
	fprintf(stderr,"missing properties parameter\n");
	return 0;
      }

      // Now parse properties
      strcpy(linebuffer, atoms->param_value[properties_idx]);
      p = linebuffer;
      k = 0;
      while ((p1 = strsep(&p, ":")) != NULL) {
	strcpy(fields[k++], p1);
      }
      
      atoms->n_property = 0;
      atoms->n_int = atoms->n_real = atoms->n_str = atoms->n_logical = 0;
      entry_count = 0;
      
      for (i=0; i<k/3; i++) {
	debug("got property %s:%s:%s\n", fields[3*i], fields[3*i+1], fields[3*i+2]);
	
	if (sscanf(fields[3*i+2], "%d", &ncols) != 1) {
	  fprintf(stderr,"Bad column count %s\n", fields[3*i+2]);
	  return 0;
	}

	entry_count += ncols;
	if (entry_count > MAX_ENTRY_COUNT) {
	  fprintf(stderr,"Maximum entry count(%d) exceeded\n", MAX_ENTRY_COUNT);
	  return 0;
	}

	if (strcmp(fields[3*i+1],"I") == 0) {
	  atoms_add_property(atoms, fields[3*i], PROPERTY_INT, ncols, -1, NETCDF_IN);
	} else if (strcmp(fields[3*i+1],"R") == 0) {
	  atoms_add_property(atoms, fields[3*i], PROPERTY_REAL, ncols, -1, NETCDF_IN);
	} else if (strcmp(fields[3*i+1],"S") == 0) {
	  atoms_add_property(atoms, fields[3*i], PROPERTY_STR, ncols, -1, NETCDF_IN);
	} else if (strcmp(fields[3*i+1],"L") == 0) {
	  atoms_add_property(atoms, fields[3*i], PROPERTY_LOGICAL, ncols, -1, NETCDF_IN);
	} else  {
	  fprintf(stderr,"Bad property type %s\n", fields[3*i+1]);
	  return 0;
	}
      }

      if (atoms->filter != NULL) free(atoms->filter);
      atoms->filter = malloc(nxyz*sizeof(int));
      if (atoms->filter == NULL) {
	fprintf(stderr,"Error allocating atoms->filter\n");
	return 0;
      }
      if (atomlist) {
	for (i=0; i<atoms->n_atom; i++) atoms->filter[i] = 0;
	for (i=0; i<natomlist; i++) {
	  if(atomlist[i] < 0 || atomlist[i] >= nxyz)  {
	    fprintf(stderr,"filter atom %d out of range 0 < i < %d\n", atomlist[i], nxyz);
	    return 0;
	  }
	  atoms->filter[atomlist[i]] = 1;
	}
	atoms->n_atom = natomlist;
	fprintf(stderr, "filter applied. Selected %d/%d atoms.\n", natomlist, nxyz);
	// Add an int property for original index
	original_index = atoms_add_property(atoms, "original_index", PROPERTY_INT, 1, -1, NETCDF_IN);
      }
      else {
	// Include all atoms
	atoms->n_atom = nxyz;
	for (i=0; i<atoms->n_atom; i++) atoms->filter[i] = 1;
      }

      if (realloc) atoms_realloc(atoms);
      atoms->initialised = 1;

    } /* end if (redefine) */
  } else {
    if (atomlist) {
      natomlist = atoms->n_atom;
    } else {
      nxyz = atoms->n_atom;
    }
  }

  if (query) return atoms->n_atom;

  // Check number of atoms
  if (atomlist) {
    if (natomlist != atoms->n_atom)  {
      fprintf(stderr,"Mismatch in number of atoms - expecting %d but got %d\n", (int)atoms->n_atom, natomlist);
      return 0;
    }
  } else {
    if (nxyz != atoms->n_atom) {
      fprintf(stderr,"Mismatch in number of atoms - expecting %d but got %d\n", (int)atoms->n_atom, nxyz);
      return 0;
    }
  }

  // Read parameters, according to definition
  for (i=0; i<nfields; i++) {
    if ((p = strchr(finalfields[i],'=')) == NULL)  {
      fprintf(stderr,"Badly formed key/value pair %s\n", finalfields[i]);
      return 0;
    }
    *p = '\0';
    
    // Look up key in definition
    j = atoms_find_param(atoms, finalfields[i]);
    if (j == -1) {
      fprintf(stderr,"Unknown parameter %s", finalfields[i]);
      return 0;
    }

    strcpy(atoms->param_value[j], p+1);
    if (strcasecmp(atoms->param_key[j], "Lattice") == 0 ||
	strcasecmp(atoms->param_key[j], "Properties") == 0) continue;

    strcpy(linebuffer, atoms->param_value[j]);
    if (atoms->param_type[j] == T_INTEGER_A || atoms->param_type[j] == T_REAL_A) {
      k = 0;
      p = linebuffer;
      while ((p1 = strsep(&p, " ")) != NULL) {
	if (*p1 == '\0') continue;
	strcpy(fields[k++], p1);
      }
      if (atoms->param_size[j] != k) {
	fprintf(stderr,"Mismatch in number of fields in parameter %s - expected %d but got %d\n", 
		atoms->param_key[j], atoms->param_size[j], k);
	return 0;
      }
    }

    switch(atoms->param_type[j]) {
    case(T_INTEGER):
      atoms->param_int[j] = strtol(atoms->param_value[j], &p, 10);
      break;
    case(T_REAL):
      atoms->param_real[j] = strtod(atoms->param_value[j], &p);
      break;
    case(T_INTEGER_A):
      for (m=0; m<k; m++)
	atoms->param_int_a[j][m] = strtol(fields[m], &p, 10);
      break;
    case(T_REAL_A):
      for (m=0; m<k; m++)
	atoms->param_real_a[j][m] = strtod(fields[m], &p);
      break;
    case(T_CHAR):
      break;
    default:
      fprintf(stderr,"Unkown param type %d\n", atoms->param_type[j]);
      return 0;
    }
  }

  // Read lattice
  if (override_lattice) {
    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
	atoms->lattice[i][j] = lattice[i][j];
  }
  else {
    lattice_idx = atoms_find_param(atoms, "Lattice");
    if (lattice_idx == -1)    {
      fprintf(stderr,"missing lattice parameter\n");
      return 0;
    }
    if (sscanf(atoms->param_value[lattice_idx], "%lf %lf %lf %lf %lf %lf %lf %lf %lf", 
	       &(atoms->lattice[0][0]), &(atoms->lattice[1][0]), &(atoms->lattice[2][0]), 
	       &(atoms->lattice[0][1]), &(atoms->lattice[1][1]), &(atoms->lattice[2][1]),
	       &(atoms->lattice[0][2]), &(atoms->lattice[1][2]), &(atoms->lattice[2][2])) != 9) {
      fprintf(stderr,"error reading lattice from string %s\n", atoms->param_value[lattice_idx]);
      return 0;
    }
  }

  // Check properties match definition
  properties_idx = atoms_find_param(atoms, "Properties");
  if (properties_idx == -1) {
    fprintf(stderr,"missing properties parameter\n");
    return 0;
  }
  
  strcpy(linebuffer, atoms->param_value[properties_idx]);
  p = linebuffer;
  k = 0;
  while ((p1 = strsep(&p, ":")) != NULL) {
    strcpy(fields[k++], p1);
  }
  
  if (k/3 + (atomlist ? 1 : 0) != atoms->n_property)  {
    fprintf(stderr,"Incorrect number of properties - expecting %d but got %d\n", atoms->n_property, k/3);
    return 0;
  }

  entry_count = 0;
  for (i=0; i<k/3; i++) {
    if(strcmp(atoms->property_name[i], fields[3*i]) != 0) {
      fprintf(stderr,"Mismatch in property %d - %s != %s\n", i, atoms->property_name[i], fields[3*i]);
      return 0;
    }

    if (sscanf(fields[3*i+2], "%d", &ncols) != 1) {
      fprintf(stderr,"Bad column count %s\n", fields[3*i+2]);
      return 0;
    }

    if (ncols != atoms->property_ncols[i]) {
      fprintf(stderr,"Mismatch in column count for property %s - expecting %d but got %d\n",
	      atoms->property_name[i], atoms->property_ncols[i], ncols);
      return 0;
    }

    entry_count += ncols;

    if (strcmp(fields[3*i+1],"I") == 0) {
      if (atoms->property_type[i] != PROPERTY_INT) {
	fprintf(stderr,"Mismatch in property type for property %s - expecting I but got %s\n",
		atoms->property_name[i], fields[3*i+1]);
	return 0;
      }
    } else if (strcmp(fields[3*i+1],"R") == 0) {
      if (atoms->property_type[i] != PROPERTY_REAL) {
	fprintf(stderr,"Mismatch in property type for property %s - expecting R but got %s\n",
		atoms->property_name[i], fields[3*i+1]);
	return 0;
      }
    } else if (strcmp(fields[3*i+1],"S") == 0) {
      if (atoms->property_type[i] != PROPERTY_STR) {
	fprintf(stderr,"Mismatch in property type for property %s - expecting S but got %s\n",
		atoms->property_name[i], fields[3*i+1]);
	return 0;
      }
    } else if (strcmp(fields[3*i+1],"L") == 0) {
      if (atoms->property_type[i] != PROPERTY_LOGICAL)	{
	fprintf(stderr,"Mismatch in property type for property %s - expecting L but got %s\n",
		atoms->property_name[i], fields[3*i+1]);
	return 0;
      }
    } else {
      fprintf(stderr,"Bad property type %s\n", fields[3*i+1]);
      return 0;
    }
  }

  // Zero string data - this is needed as we're not null terminating our strings
  memset(atoms->str_data, ' ', atoms->n_atom*atoms->n_str*PROPERTY_STRING_LENGTH);

  // Now it's just one line per atom
  n = 0;
  for (atidx=0; atidx < nxyz; atidx++) {
    if (!fgets(linebuffer,LINESIZE,in)) {
      fprintf(stderr, "premature file ending at atom %d\n",n);
      return 0;
    }
    if (!atoms->filter[atidx]) continue;

    k = 0;
    p = linebuffer;
    while ((p1 = strsep(&p, " \t\n")) != NULL) {
      if (*p1 == '\0') continue;
      strcpy(fields[k++], p1);
    }
    if (k != entry_count) {
      fprintf(stderr, "incomplete row, atom %d - got %d/%d entries\n", n, k, entry_count);
      for (i=0;i<k;i++) fprintf(stderr, "fields[%d] = %s, length %lu\n", i, fields[i], (unsigned long)strlen(fields[i]));
      return 0;
    }

    k = 0;
    for (i=0; i<atoms->n_property; i++) {
      if (atomlist && strcmp(atoms->property_name[i],"original_index") == 0) continue;
      switch(atoms->property_type[i]) {
      case(PROPERTY_INT): 
	
	for (j=0; j < atoms->property_ncols[i]; j++)
	  //if (sscanf(fields[k+j], "%d", &atoms->int_data[(atoms->property_start[i] + j)*(atoms->n_atom) + n]) != 1) 
	  if (sscanf(fields[k+j], "%d", &property_int(atoms, i, j, n)) != 1)  {
	    fprintf(stderr,"Can't convert int value %s\n", fields[k+j]);
	    return 0;
	  }
	k += atoms->property_ncols[i];
	break;

      case(PROPERTY_REAL): 
	for (j=0; j < atoms->property_ncols[i]; j++)
	  //if (sscanf(fields[k+j], "%lf", &atoms->real_data[(atoms->property_start[i] + j)*(atoms->n_atom) + n]) != 1) 
	  if (sscanf(fields[k+j], "%lf", &property_real(atoms, i, j, n)) != 1)  {
	    fprintf(stderr,"Can't convert real value %s\n", fields[k+j]);
	    return 0;
	  }
	k += atoms->property_ncols[i];
	break;

      case(PROPERTY_STR): 
	for (j=0; j < atoms->property_ncols[i]; j++) {
	  //if (sscanf(fields[k+j], "%10c", &atoms->str_data[PROPERTY_STRING_LENGTH*((atoms->property_start[i] + j)*(atoms->n_atom) + n)]) != 1) 
	  if (sscanf(fields[k+j], "%10c", &property_str(atoms, i, j, n)) != 1)  {
	    fprintf(stderr,"Can't convert str value %s\n", fields[k+j]);
	    return 0;
	  }
	}
	k += atoms->property_ncols[i];
	break;

      case(PROPERTY_LOGICAL): 
	for (j=0; j < atoms->property_ncols[i]; j++) {
	  if (sscanf(fields[k+j], "%c", &tmp_logical) != 1)  {
	    fprintf(stderr,"Can't convert logical value %s\n", fields[k+j]);
	    return 0;
	  }
	  if (tmp_logical == 'T') 
	    //atoms->logical_data[(atoms->property_start[i] + j)*(atoms->n_atom) + n] = 1;
	    property_logical(atoms,i,j,n) = 1;
	  else
	    //atoms->logical_data[(atoms->property_start[i] + j)*(atoms->n_atom) + n] = 0;
	    property_logical(atoms,i,j,n) = 0;
	}
	k += atoms->property_ncols[i];
	break;
      }
    }
    //if (atomlist) atoms->int_data[(atoms->property_start[original_index])*(atoms->n_atom) + n] = atidx;
    if (atomlist) property_int(atoms, original_index, 0, n) = atidx;
    n++;
  }

  return 1;
} 


int write_xyz(FILE *out, Atoms *atoms, char *int_format, char *real_format, char *str_format, char *logical_format, int swap) {
  char linebuffer[LINESIZE];
  char intf[LINESIZE], realf[LINESIZE], strf[LINESIZE], logf[LINESIZE];
  int i, j, n, lattice_idx, properties_idx;
  char *trimmed;
  int species_idx, pos_idx;

  if (out != stdout && fseek(out, 0, SEEK_END) != 0) {
    fprintf(stderr,"Cannot seek to end of file\n");
    return 0;
  }

  if (swap) {
    // Move species to first property
    species_idx = atoms_find_property(atoms, "species");
    if (species_idx == -1) {
      fprintf(stderr,"No species property, cannot write XYZ file\n");
      return 0;
    }
    atoms_swap_properties(atoms, species_idx, 0);

    // Move pos to second property
    pos_idx = atoms_find_property(atoms, "pos");
    if (pos_idx == -1) {
      fprintf(stderr,"No pos property, cannot write XYZ file\n");
      return 0;
    }
    atoms_swap_properties(atoms, pos_idx, 1);
  }

  lattice_idx = atoms_find_param(atoms, "Lattice");
  if (lattice_idx == -1) 
    lattice_idx = atoms_add_param(atoms, "Lattice", T_CHAR, 1, -1, NETCDF_IN);

  properties_idx = atoms_find_param(atoms, "Properties");
  if (properties_idx == -1) 
    properties_idx = atoms_add_param(atoms, "Properties", T_CHAR, 1, -1, NETCDF_IN);
  
  sprintf(atoms->param_value[lattice_idx], "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
	  atoms->lattice[0][0], atoms->lattice[1][0], atoms->lattice[2][0],
	  atoms->lattice[0][1], atoms->lattice[1][1], atoms->lattice[2][1],
	  atoms->lattice[0][2], atoms->lattice[1][2], atoms->lattice[2][2]);

  atoms->param_value[properties_idx][0] = '\0';
  for (i=0; i<atoms->n_property; i++) {
    if (!atoms->property_filter[i]) continue;
    sprintf(atoms->param_value[properties_idx], "%s%s:%c:%d:", atoms->param_value[properties_idx], 
	    atoms->property_name[i],
	    PROPERTY_STRINGS[atoms->property_type[i]],
	    atoms->property_ncols[i]);
  }
  atoms->param_value[properties_idx][strlen(atoms->param_value[properties_idx])-1] = '\0';

  sprintf(intf, "%%s%s", int_format);
  sprintf(realf, "%%s%s", real_format);
  sprintf(strf, "%%s%s", str_format);
  sprintf(logf, "%%s%s", logical_format);

  // Build parameter values
  linebuffer[0] = '\0';
  for (i=0; i<atoms->n_param; i++) {
    if (!atoms->param_filter[i]) continue;
    if (atoms->param_type[i] == T_INTEGER)
      sprintf(atoms->param_value[i], int_format, atoms->param_int[i]);
    else if (atoms->param_type[i] == T_REAL)
      sprintf(atoms->param_value[i], real_format, atoms->param_real[i]);
    else if (atoms->param_type[i] == T_INTEGER_A) {
      atoms->param_value[i][0]='\0';
      for (j=0; j<3; j++)
	sprintf(atoms->param_value[i], intf, atoms->param_value[i], atoms->param_int_a[i][j]);
    } else if (atoms->param_type[i] == T_REAL_A) {
      atoms->param_value[i][0]='\0';
      for (j=0; j<3; j++)
	sprintf(atoms->param_value[i], realf, atoms->param_value[i], atoms->param_real_a[i][j]);
    }
    trimmed = atoms->param_value[i];
    while (isblank(trimmed[0])) trimmed++;

    sprintf(linebuffer, "%s%s=%s%s%s ", linebuffer, atoms->param_key[i],
	    strchr(trimmed,' ') != NULL ? "\"" : "", 
	    trimmed,
	    strchr(trimmed,' ') != NULL ? "\"" : "");
  }

  fprintf(out, "%d\n", (int)atoms->n_atom);
  linebuffer[strlen(linebuffer)-1] = '\n';
  fputs(linebuffer, out);

  for (n=0; n<atoms->n_atom; n++) {
    linebuffer[0] = '\0';

    for (i=0; i<atoms->n_property; i++) {
      if (!atoms->property_filter[i]) continue;
      switch(atoms->property_type[i]) {
      case(PROPERTY_INT): 
	for (j=0; j < atoms->property_ncols[i]; j++)
	  //sprintf(linebuffer, intf, linebuffer, atoms->int_data[(atoms->property_start[i] + j)*(atoms->n_atom) + n]);
	  sprintf(linebuffer, intf, linebuffer, property_int(atoms, i, j, n));
	break;

      case(PROPERTY_REAL): 
	for (j=0; j < atoms->property_ncols[i]; j++)
	  //sprintf(linebuffer, realf, linebuffer, atoms->real_data[(atoms->property_start[i] + j)*(atoms->n_atom) + n]);
	  sprintf(linebuffer, realf, linebuffer, property_real(atoms, i, j, n));
	break;
      
      case(PROPERTY_STR):
	if (linebuffer[0] != '\0') strcat(linebuffer, " ");
	for (j=0; j < atoms->property_ncols[i]; j++)
	  //sprintf(linebuffer, strf, linebuffer, &atoms->str_data[PROPERTY_STRING_LENGTH*((atoms->property_start[i] + j)*(atoms->n_atom) + n)]);
	  sprintf(linebuffer, strf, linebuffer, &property_str(atoms, i, j, n));
	break;
      
      case(PROPERTY_LOGICAL): 
	for (j=0; j < atoms->property_ncols[i]; j++) {
	  //sprintf(linebuffer, logf, linebuffer, atoms->logical_data[(atoms->property_start[i] + j)*(atoms->n_atom) + n] ? 'T' : 'F');
	  sprintf(linebuffer, logf, linebuffer, property_logical(atoms, i, j, n));
	}
	break;
      }
    }
    strcat(linebuffer, "\n");
    fputs(linebuffer, out);
  }

  atoms->n_frame++;
  fflush(out);
  return 1;
}

int find_in_list(int *list, int n, int value)
{
  int i;

  for (i=0; i<n; i++) 
    if (list[i] == value) return 1;
  return 0;
}

void print_property(Atoms *at, char *name, char *intformat, char *realformat, char *strformat, char *logformat)
{
  int j, k, n;
  char linebuffer[LINESIZE], intf[10], realf[10], logf[10], strf[10];

  sprintf(intf, "%%s%s", intformat);
  sprintf(realf, "%%s%s", realformat);
  sprintf(strf, "%%s%s", strformat);
  sprintf(logf, "%%s%s", logformat);
  
  j = atoms_find_property(at, name);
  if (j == -1) return;

  for (n = 0; n < at->n_atom; n++) {
    linebuffer[0] ='\0';
    switch(at->property_type[j]) {
    case(PROPERTY_INT): 
      for (k=0; k < at->property_ncols[j]; k++)
	//sprintf(linebuffer, intf, linebuffer, at->int_data[(at->property_start[j] + k)*(at->n_atom) + n]);
	sprintf(linebuffer, intf, linebuffer, property_int(at, j, k, n));
      break;

    case(PROPERTY_REAL): 
      for (k=0; k < at->property_ncols[j]; k++)
	//sprintf(linebuffer, realf, linebuffer, at->real_data[(at->property_start[j] + k)*(at->n_atom) + n]);
      	sprintf(linebuffer, realf, linebuffer, property_real(at, j, k, n));
      break;
      
    case(PROPERTY_STR): 
      if (linebuffer[0] != '\0') strcat(linebuffer, " ");
      for (k=0; k < at->property_ncols[j]; k++)
	//sprintf(linebuffer, strf, linebuffer, &at->str_data[PROPERTY_STRING_LENGTH*((at->property_start[j] + k)*(at->n_atom) + n)]);
      	sprintf(linebuffer, strf, linebuffer, &property_str(at, j, k, n));
      break;
		
    case(PROPERTY_LOGICAL): 
      for (k=0; k < at->property_ncols[j]; k++) 
	//sprintf(linebuffer, logf, linebuffer, at->logical_data[(at->property_start[j] + k)*(at->n_atom) + n] ? 'T' : 'F');
	sprintf(linebuffer, logf, linebuffer, property_logical(at, j, k, n) ? 'T' : 'F');
      break;
    }
    strcat(linebuffer,"\n");
    printf(linebuffer);
  }
}

int sprint_param(char *linebuffer, Atoms *at, char *name, char *intformat, char *realformat)
{
  int j;
  char fmt[LINESIZE];

  j = atoms_find_param(at, name);
  if (j == -1) return 0;

  switch(at->param_type[j]) {
  case(T_INTEGER):
    sprintf(fmt, "%%s%s %s ", name, intformat);
    sprintf(linebuffer, fmt, linebuffer, at->param_int[j]);
    break;
  case(T_REAL):
    sprintf(fmt, "%%s%s %s ", name, realformat);
    sprintf(linebuffer, fmt, linebuffer, at->param_real[j]);
    break;
  case(T_INTEGER_A):
    sprintf(fmt, "%%s%s %s %s %s ", name, intformat, intformat, intformat);
    sprintf(linebuffer, fmt, linebuffer, at->param_int_a[j][0], at->param_int_a[j][1], at->param_int_a[j][2]);
    break;
  case(T_REAL_A):
    sprintf(fmt, "%%s%s %s %s %s ", name, realformat, realformat, realformat);
    sprintf(linebuffer, fmt, linebuffer, at->param_real_a[j][0], at->param_real_a[j][1], at->param_real_a[j][2]);
    break;
  case(T_CHAR):
    sprintf(linebuffer,"%s %s %s", linebuffer, name, at->param_value[j]);
    break;
  default:
    fprintf(stderr,"Unknown param type %d\n", at->param_type[j]);
    return 0;
  }
  return 1;
}


#ifdef MAIN_PROGRAM
void usage(int nc2xyz, int xyz2nc, int xyz2xyz, int nc2nc, int xyzstat, int ncstat)
{
  if (nc2xyz) {
    printf("Usage:\n  nc2xyz [OPTIONS] INFILE OUTFILE\n\n");
    printf("Converts INFILE in NetCDF format to OUTFILE in Extended XYZ format.\n");
    printf("Use - in place of OUTFILE to write the XYZ to stdout.\n\n");
  }
  else if (xyz2nc) {
    printf("Usage:\n  xyz2nc [OPTIONS] INFILE OUTFILE\n\n");
    printf("Converts INFILE in Extended XYZ format to OUTFILE in NetCDF format.\n");
    printf("Use - in place of INFILE to read the XYZ from stdin; in this case\n");
    printf("the -r, -c, -f, -n and -o options can't be used.\n\n");
  } else if (xyz2xyz) {
    printf("Usage:\n  xyz2xyz [OPTIONS] INFILE OUTFILE\n\n");
    printf("Read INFILE in Extended XYZ format and write to OUTFILE in XYZ format.\n");
    printf("Use \"-\" in place of INFILE to read the XYZ from stdin and/or use\n");
    printf("\"-\" in place of OUTFILE to write to stdout.\n\n");
  } else if (nc2nc) {
    printf("Usage:\n  nc2nc [OPTIONS] INFILE OUTFILE\n\n");
    printf("Read INFILE in NetCDF format and write to OUTFILE in NetCDF format.\n\n");
  } else if (xyzstat) {
    printf("Usage:\n  xyzstat [OPTIONS] INFILE\n\n");
    printf("Print statistics about the Extended XYZ file INFILE.\n");
    printf("With no options, xyzstat creates or updates an XYZ index file. The index for\n");
    printf("\"file.xyz\" is stored in text format in \"file.xyz.idx\".\n\n");
  } else if (ncstat) {
    printf("Usage:\n  ncstat [OPTIONS] INFILE\n\n");
    printf("Print statistics about the NetCDF file INFILE.\n");
  }

  printf("\n(Other conversion modes can be activated by naming the executable with one\n");
  printf("of the following variants: xyz2nc, nc2xyz, xyz2xyz, nc2nc, xyzstat, ncstat)\n\n");

  printf("General Options:\n");
  printf("  -h          Print this help.\n");
  printf("  -r RANGE    Range of frames to include. Should be either a single frame\n");
  printf("              number or a slice [start]:[stop][:step]. If -f is omitted,\n");
  printf("              default is all frames. With C indexing (default), slices are\n");
  printf("              exclusive of the \"stop\" frame: think of the loop\n");
  printf("                 for (i=start; i<stop; i+=step).\n"); 
  printf("              With Fortran indexing, slices include the \"stop\" frame.\n");
  printf("              Negative indices count backwards from the end of the file, with\n");
  printf("              -1 being the last frame for both indexing schemes. Here are some\n");
  printf("              examples of valid ranges:\n\n");
  printf("                C indexing    Fortran indexing    Description\n");
  printf("                ---------------------------------------------------------------\n");
  printf("                   0               1              first frame in file\n");
  printf("                  -1              -1              last frame in file\n");
  printf("                   0:              1:             all frames\n");
  printf("                   1:3             2:3            2nd and 3rd frames in file\n");
  printf("                   :-2             :-2            all but the last two frames\n");
  printf("                   ::2             ::2            entire file, every 3nd frame.\n\n");
  printf("  -a ATOMFILE  Restrict atoms selected to those indices listed in ATOMFILE.\n");
  printf("  -F           Use one-based  (Fortran-style) indexing for frames and atoms\n\n");
  printf("  -3           Use NetCDF 3 file format\n\n");

  printf("\nOptions specific to this variant:\n");

  if (ncstat || xyzstat) {
    printf("  -c           Print number of frames selected.\n");
    printf("  -f           Print indices of selected frames.\n");
    printf("  -p           List properties (per-atom variables).\n");
    printf("  -P           List parameters (per-frame variables).\n");
    if (xyzstat) {
      printf("  -n           Print number of selected atoms in selected frames (one per line).\n");
      printf("  -o           Print file offsets of start of selected frames (one per line).\n");
    }
    else {
      printf("  -n           Print number of selected atoms.\n");
    }
    printf("  -d VAR       Print value of variable for selected atoms and frames.\n");
  } else {
    printf("  -L LATTICE   Override lattice in all output frames with LATTICE, given in\n");
    printf("               Fortran ordering as \"R11 R21 R31 R12 R22 R32 R13 R23 R33\". If\n");
    printf("               three fields are given then the lattice is assumed to be cubic with\n");
    printf("               a=R11, b=R22, c=R33 and all other components zero. One field implies\n");
    printf("               a cubic cell with a=b=c=R11.\n");
    printf("  -p PROPS     Comma or colon separated list of properties (i.e. per-atom variables)\n");
    printf("               that should be included in output file, e.g.\n");
    printf("                 -p species,pos,velo\n");
    printf("               The default is to print all properties.\n");
    if (nc2xyz || xyz2xyz) {
      printf("               For a valid XYZ file, the first should be \"species\" and the second\n");
      printf("               should be something position-like. The order of properties is as\n");
      printf("               specified in this command.\n");
    }
    printf("  -P PARAMS    Comma or colon separated list of parameters (per-frame variables) to\n");
    printf("               include in output file. For example:\n");
    printf("                 -P Energy,MaxForce\n");
    printf("               The default is to output all parameters.\n");
    if (nc2xyz || xyz2xyz) {
      printf("               For XYZ output these parameters are in addition to the special\n"); 
      printf("               parameters \"Lattice\" and \"Properties\" which are always written.\n");
      printf("               The original order of parameters in the input file is preserved.\n");
    }
  }
  if (nc2xyz || xyz2xyz || ncstat || xyzstat) {
    printf("  -I INTFMT    printf(1) format string used to print integer variables.\n");
    printf("               The default is \"%%8d\".\n");
    printf("  -R REALFMT   printf(1) format string used to print real variables.\n");
    printf("               The default is \"%%16.8f\".\n");
  }
  if (ncstat || nc2nc || nc2xyz) {
    printf("  -z           Replace missing numerical values in NetCDF input with zeros.\n");
    printf("  -Z LEVEL     Set zlib deflate level of output file (0=No compression, 9=max).\n");
    printf("               Default is level 6.\n");
  }

  exit(0);
}


void parse_range(int rflag, int nframes, int offset, int gotcolon, 
		 char *start, char *stop, char *step, 
		 int *f_start, int *f_stop, int *f_step) 
{
  if (rflag) {
    // Calculate start:stop:step range, filling in defaults
    if (start == NULL || strcmp(start, "") == 0)
      *f_start = 0;
    else {
      if (sscanf(start, "%d", f_start) != 1) pe("Bad start frame %s\n", start);
      if (*f_start >= 0) *f_start -= offset;
    }
	  
    if (stop  == NULL || strcmp(stop, "") == 0) 
      if (gotcolon)
	*f_stop = nframes;
      else
	if (*f_start == -1)
	  *f_stop = nframes;
	else
	  *f_stop = *f_start + 1;
    else {
      if (sscanf(stop, "%d", f_stop) != 1) pe("Bad stop frame %s\n", stop);
    }

    if (step == NULL || strcmp(step, "") == 0)
      *f_step = 1;
    else
      if (sscanf(step, "%d", f_step) != 1) pe("Bad step %s\n", step);

    if (*f_start < 0) *f_start = nframes+*f_start;
    if (*f_stop < 0)  *f_stop = nframes+*f_stop;
    

    debug("start %d stop %d step %d\n", *f_start, *f_stop, *f_step);

    // Check for stupid ranges
    if (*f_step > 0 && *f_start >= *f_stop) pe("Start frame is after stop frame and step is positive\n");
    if (*f_step < 0 && *f_stop  >= *f_stop) pe("Stop frame is after start frame and step is negative\n");
    if (*f_stop > nframes) pe("Stop frame (%d) is after last frame (%d)\n",*f_stop,nframes);
    if (*f_step == 0) pe("Step is zero\n");
  }
  else {
    *f_start = 0;
    *f_stop  = nframes;
    *f_step  = 1;
  }
}

int main (int argc, char **argv) 
{
  Atoms at;
  int i, j, k, n, nframes, nsel, res;
  int ch, cflag, rflag, nflag, aflag, oflag, fflag, pflag, Pflag, Lflag, dflag, zflag, Zflag, nc3flag;
  int f_start, f_stop, f_step, gotcolon, offset;
  char *start, *stop, *step;
  char linebuffer[LINESIZE], exename[LINESIZE], optstr[LINESIZE];
  FILE *infile, *atomfile;
  FILE *outfile;
  int *atomlist, natomlist;
  char infilename[LINESIZE], outfilename[LINESIZE];
  char intformat[10], realformat[10], logformat[10], strformat[10];
  int xyz2xyz, xyz2nc, nc2xyz, nc2nc, xyzstat, ncstat;
#ifdef HAVE_NETCDF
  int retval;
#endif
  int printed_stats;
  char *p, *p1, nprops, nparams, params[MAX_ENTRY_COUNT][LINESIZE],
    props[MAX_ENTRY_COUNT][LINESIZE], dvar[LINESIZE];
  double lattice[3][3] = {{0.0,0.0,0.0},
			  {0.0,0.0,0.0},
			  {0.0,0.0,0.0}};
  int shuffle, deflate, deflate_level;
  int nc_in, nc_out;
  int allow_redefine = 0;

  // default deflation settings for now
  shuffle = 1; deflate = 1; deflate_level = 6;
  
  strcpy(exename, argv[0]);
  xyz2xyz  = strstr(exename,"xyz2xyz")   != NULL;  
  xyz2nc   = strstr(exename,"xyz2nc")    != NULL;   
  nc2xyz   = strstr(exename,"nc2xyz")    != NULL;   
  nc2nc    = strstr(exename,"nc2nc")     != NULL;   
  xyzstat  = strstr(exename,"xyzstat")   != NULL;
  ncstat   = strstr(exename,"ncstat")    != NULL;

  if (xyz2xyz + xyz2nc + nc2xyz + xyzstat + nc2nc + ncstat != 1)
    pe("Executable %s should be one of:\nxyz2nc, nc2xyz, xyz2xyz, nc2nc, xyzstat, ncstat\n", 
       exename);

  if(xyz2xyz || xyzstat) allow_redefine = 1;

  cflag = 0;
  rflag = 0;
  nflag = 0;
  aflag = 0;
  oflag = 0;
  fflag = 0;
  pflag = 0;
  Pflag = 0;
  Lflag = 0;
  dflag = 0;
  zflag = 0;
  Zflag = 0;
  offset = 0;
  nc3flag = 0;

  strcpy(intformat, "%8d");
  strcpy(realformat, "%16.8f");
  strcpy(logformat, "%5c");
  strcpy(strformat, "%.10s");

  // Build option string
  strcpy(optstr, "hr:a:F3");
  if (xyzstat || ncstat) 
    strcat(optstr, "cfnvpPd:");
  else
    strcat(optstr, "p:P:L:");
  if (xyzstat) 
    strcat(optstr, "o");
  if (nc2xyz || xyz2xyz || ncstat || xyzstat)
    strcat(optstr, "I:R:");
  if (nc2nc || nc2xyz || ncstat)
    strcat(optstr, "zZ:");

  while ((ch = getopt(argc, argv, optstr)) != -1) {
    switch (ch) {
      /* general options */
    case 'h':
      usage(nc2xyz, xyz2nc, xyz2xyz, nc2nc, xyzstat, ncstat);
      break;

    case 'r':
      rflag = 1;
      start = strsep(&optarg,":");
      gotcolon = optarg != NULL;
      stop = strsep(&optarg,":");
      step = strsep(&optarg, ":");
      break;

    case 'a':
      aflag = 1;
      if ((atomfile = fopen(optarg, "r")) == NULL) pe("Atom file %s cannot be opened", optarg);
      break;

    case 'F':
      offset=1;
      break;

    case '3':
      nc3flag=1;
      break;

      /* end of general options */

    case 'c':
      cflag = 1;
      break;

    case 'f':
      fflag = 1;
      break;

    case 'n':
      nflag = 1;
      break;

    case 'o':
      oflag = 1;
      break;

    case 'p':
      pflag = 1;
      if (!(xyzstat || ncstat)) {
	k = 0;
	p = optarg;
	while ((p1 = strsep(&p, ",:")) != NULL) {
	  if (*p1 == '\0') continue;
	  strcpy(props[k],p1);
	  if (++k == MAX_ENTRY_COUNT) pe("Too many properties specified on command line");
	}
	nprops = k;
#ifdef DEBUG
	debug("Got %d Properties:\n", nprops);
	for (i=0; i<nprops; i++)
	  printf("%s\n",props[i]);
#endif
      }
      break;
    case 'P':
      Pflag = 1;
      if (!(xyzstat || ncstat)) {
	k = 0;
	p = optarg;
	while ((p1 = strsep(&p, ",:")) != NULL) {
	  if (*p1 == '\0') continue;
	  strcpy(params[k],p1);
	  if (++k == MAX_ENTRY_COUNT) pe("Too many parameters specified on command line");
	}
	nparams = k;
#ifdef DEBUG
	debug("Got %d Parameters:\n", nparams);
	for (i=0; i<nparams; i++)
	  printf("%s\n",params[i]);
#endif
      }
      break;

    case 'd':
      dflag = 1;
      strcpy(dvar, optarg);
      debug("Got dump var %s\n" , dvar);
      break;

    case 'L':
      Lflag = 1;
      k = 0;
      strcpy(linebuffer, optarg);
      p = linebuffer;
      while ((p1 = strsep(&p, " \t")) != NULL) {
	if (*p1 == '\0') continue;
	k++;
      }
      if (k == 9) {
	if (sscanf(optarg, "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
		   &lattice[0][0], &lattice[1][0], &lattice[2][0],
		   &lattice[0][1], &lattice[1][1], &lattice[2][1],
		   &lattice[0][2], &lattice[1][2], &lattice[2][2]) != 9)
	  pe("Bad 9 component lattice %s", optarg);
      } else if (k == 3) {
	if (sscanf(optarg, "%lf %lf %lf", &lattice[0][0], &lattice[1][1], &lattice[2][2]) != 3)
	  pe("Bad 3 component lattice %s", optarg);
      } else if (k == 1) {
	if (sscanf(optarg, "%lf", &lattice[0][0]) != 1) 
	  pe("Bad 1 component lattice %s", optarg);
	lattice[1][1] = lattice[2][2] = lattice[0][0];
      } else
	pe("Bad lattice: must have 1, 3 or 9 components, but got %d", k);
      
      debug("got lattice =\n%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n",
	    lattice[0][0], lattice[1][0], lattice[2][0],
	    lattice[0][1], lattice[1][1], lattice[2][1],
	    lattice[0][2], lattice[1][2], lattice[2][2]);
      break;

    case 'I':
      strcpy(intformat, optarg);
      if (strchr(intformat, '%') == NULL) pe("Invalid int format string %s", intformat);
      debug("intformat %s\n", intformat);
      break;

    case 'R':
      strcpy(realformat, optarg);
      if (strchr(realformat, '%') == NULL) pe("Invalid real format string %s", realformat);
      debug("realformat %s\n", realformat);
      break;

    case 'z':
      zflag = 1;
      break;

    case 'Z':
      Zflag = 1;
      if (sscanf(optarg, "%d", &deflate_level) != 1) pe("Invalid deflate level %s\n", optarg);
      if (deflate_level < 0 || deflate_level > 9) pe("Deflate level (%d) must be in range 0..9\n", deflate_level);
      deflate = deflate_level != 0;
      break;

    case '?':
    default:
      usage(nc2xyz, xyz2nc, xyz2xyz, nc2nc, xyzstat, ncstat);
    }
  }
  argc -= optind;
  argv += optind;

  atomlist = NULL;
  natomlist = 0;
  if (aflag) {
    natomlist = 0;  // Count lines in atomlist file
    while (fgets(linebuffer,LINESIZE,atomfile)) natomlist++;
      
    atomlist = malloc(natomlist*sizeof(int));
    if (atomlist == NULL) pe("Error allocating atomlist");
    fseek(atomfile, 0, SEEK_SET);
    for (i=0; i<natomlist; i++) {
      if (!fgets(linebuffer,LINESIZE,atomfile)) pe("Premature end of atom file");
      if (sscanf(linebuffer, "%d", &atomlist[i]) != 1) pe("Error reading line %d of atom list: %s", i+1, linebuffer);
      atomlist[i] -= offset;
    }
    fclose(atomfile);
  }

  atoms_init(&at);

  if (xyz2xyz || xyz2nc || xyzstat) {

    if(xyz2xyz || xyz2nc) {
      if (argc != 2) usage(nc2xyz, xyz2nc, xyz2xyz, nc2nc, xyzstat, ncstat);
      strcpy(infilename, argv[0]);
      strcpy(outfilename, argv[1]);
    } else {
      if (argc != 1) usage(nc2xyz, xyz2nc, xyz2xyz, nc2nc, xyzstat, ncstat);
      strcpy(infilename, argv[0]);
    }

    if (strcmp(infilename, "-") == 0)  {
      infile = stdin;
      at.got_index = 0;

      if (xyzstat || rflag) 
	pe("Cannot create index file when reading XYZ from stdin\n");

    } else {

      // Only build index if one or more of -r, -o, -n, -f present
      // or if explicity asked to with "xyzstat" command
      at.got_index = 0;
      if (xyzstat || rflag) {
	nframes = xyz_find_frames(infilename, at.frames, at.atoms);
	if (nframes == 0) pe("Error building frame index");
	at.got_index = 1;
	parse_range(rflag, nframes, offset, gotcolon, start, stop, step, &f_start, &f_stop, &f_step);

	if (xyzstat && cflag) {
	  nsel = 0;
	  for (i=f_start; i<f_stop; i+=f_step) nsel++;
	  printf("%d frames\n",nsel);
	}
      }
      infile = fopen(infilename, "r");
    }

    if (xyz2xyz) {
      if (strcmp(outfilename, "-") == 0)
	outfile = stdout;
      else 
	outfile = fopen(outfilename, "w");
      if (outfile == NULL) pe("Cannot open XYZ file %s for writing", outfilename);
    }

    if (xyzstat && (pflag || Pflag)) {
      if (!read_xyz(infile, &at, atomlist, natomlist, 0, 1, 0, 1, 0, Lflag, lattice))
	pe("Error reading xyz header");
      if (Pflag) {
	debug("Parameters:\n");
	for (i=0; i<at.n_param; i++) {
	  if (strcasecmp(at.param_key[i], "Lattice") == 0 ||
	      strcasecmp(at.param_key[i], "Properties") == 0) continue;
	
	  printf("%s %s size %d\n", at.param_key[i], PARAM_TYPES[at.param_type[i]], at.param_size[i]);
	}
      }
      if (pflag) {
	debug("\nProperties:\n");
	for (i=0; i<at.n_property; i++) {
	  printf("%s %s size %d\n", at.property_name[i], PROPERTY_TYPES[at.property_type[i]], at.property_ncols[i]);
	}
      }
    }

    if (xyz2nc) {
      debug("Opening NetCDF file \"%s\" for writing\n", outfilename);
#ifdef HAVE_NETCDF
#ifdef NETCDF4
      if (nc3flag) {
	netcdf_check(nc_set_default_format(NC_FORMAT_64BIT, NULL));
	netcdf_check(nc_create(outfilename, NC_64BIT_OFFSET | NC_CLOBBER, &nc_out));
      } else {
	netcdf_check(nc_set_default_format(NC_FORMAT_NETCDF4, NULL));
	netcdf_check(nc_create(outfilename, NC_NETCDF4 | NC_CLOBBER, &nc_out));
      }
#else
      netcdf_check(nc_create(outfilename, NC_64BIT_OFFSET | NC_CLOBBER, &nc_out));
#endif
#else
      pe("No NetCDF support compiled in.\n");
#endif
    }

    n = 0;
    if (xyz2xyz || xyz2nc || (xyzstat && (dflag || oflag || nflag || fflag))) {
      if (infile == NULL) pe("Cannot open XYZ file for reading");
      if (at.got_index) {
	
	nframes = 0;
	for (i=f_start; i<f_stop; i+=f_step) nframes++;

	for (i=f_start; i < f_stop; i += f_step) {
	  if (!read_xyz(infile, &at, atomlist, natomlist, i, 0, allow_redefine, 1, 0, Lflag, lattice)) 
	    pe("Error reading frame %d", i);
	  
	  if (xyz2xyz || xyz2nc) {
	    if (pflag) {
	      for (j=0; j<at.n_property; j++)
		at.property_filter[j] = 0;
	      for (j=0; j<nprops; j++) {
		k = atoms_find_property(&at, props[j]);
		if (k == -1) pe("Property %s specified on command line not found", props[j]);
		at.property_filter[k] = 1;
		atoms_swap_properties(&at, k, j);
	      }
	    }
	    
	    if (Pflag) {
	      for (j=0; j<at.n_param; j++)
		at.param_filter[j] = 0;
	      for (j=0; j<nparams; j++) {
		k = atoms_find_param(&at, params[j]);
		if (k == -1) pe("Parameter %s specified on command line not found", params[j]);
		at.param_filter[k] = 1;
	      }
	    }
	    
	    if (xyz2xyz) {
	      if (!write_xyz(outfile, &at, intformat, realformat, strformat, logformat, !pflag))
		pe("Error writing xyz\n");
	    } else {
	      if (!write_netcdf(nc_out, &at, n++, 0, shuffle, deflate, deflate_level))
		pe("Error writing netcdf\n");
	    }
	    
	    fprintf(stderr, "Frame %d/%d                      %c", i+1, nframes, i + f_step >= f_stop ? '\n' : '\r');
	  } else {
	    // Print some things to stdout
	    printed_stats = 0;
	    linebuffer[0] = '\0';
	    if (fflag) sprintf(linebuffer, "frame %d ", i+offset);
	    if (nflag) sprintf(linebuffer, "%sn_atom %d ", linebuffer, (int)at.n_atom);
	    if (oflag) sprintf(linebuffer, "%soffset %ld ", linebuffer, at.frames[i]);
	  
	    if (dflag) {
	      j = atoms_find_param(&at, dvar);
	      if (j == -1) {
		j = atoms_find_property(&at, dvar);
		if (j == -1) pe("Variable %s not found", dvar);
		// It's a property. Print one atom per line
		if (fflag || nflag) {
		  // Flush stats line
		  strcat(linebuffer,"\n");
		  printf(linebuffer);
		  printed_stats = 1;
		}
		print_property(&at, dvar, intformat, realformat, strformat, logformat);
	      } else {
		// It's a param. Print one frame per line
		sprint_param(linebuffer, &at, dvar, intformat, realformat);
	      }
	    }
	    if (!printed_stats && (dflag || fflag || nflag || oflag)) {
	      strcat(linebuffer,"\n");
	      printf(linebuffer);
	    }
	  }
	}
      } else {
	while ((res = read_xyz(infile, &at, atomlist, natomlist, 0, 0, allow_redefine, 1, 1, Lflag, lattice))) {
	  debug("read frame %d\n", n);
	  
	  if (pflag) {
	    for (i=0; i<at.n_property; i++)
	      at.property_filter[i] = 0;
	    for (i=0; i<nprops; i++) {
	      j = atoms_find_property(&at, props[i]);
	      if (j == -1) pe("Property %s specified on command line not found", props[i]);
	      at.property_filter[j] = 1;
	    }
	  }

	  if (Pflag) {
	    for (i=0; i<at.n_param; i++)
	      at.param_filter[i] = 0;
	    for (i=0; i<nparams; i++) {
	      j = atoms_find_param(&at, params[i]);
	      if (j == -1) pe("Parameter %s specified on command line not found", params[i]);
	      at.param_filter[j] = 1;
	    }
	  }
	
	  if (xyz2xyz) {
	    if (!write_xyz(outfile, &at, intformat, realformat, strformat, logformat, !pflag))
	      pe("Error writing xyz\n");
	  } else {
	    if (!write_netcdf(nc_out, &at, n, 0, shuffle, deflate, deflate_level))
	      pe("Error writing netcdf\n");
	  }
	  
	  fprintf(stderr, "Frame %d                      \r", n);
	  n++;
	}
	fprintf(stderr, "\n");
	debug("read_xyz returned %d\n", res);
      }
 
      if (xyz2xyz && strcmp(outfilename, "-") != 0) fclose(outfile);
      
      if (strcmp(infilename, "-") != 0) fclose(infile);
#ifdef HAVE_NETCDF
      if (xyz2nc) netcdf_check(nc_close(nc_out));
#endif
    }
  } else if (nc2xyz || nc2nc || (ncstat && (dflag || fflag || nflag))) {

    if ((nc2xyz || nc2nc) && argc != 2)
      usage(nc2xyz, xyz2nc, xyz2xyz, nc2nc, xyzstat, ncstat);
    
    if (ncstat && argc != 1)
      usage(nc2xyz, xyz2nc, xyz2xyz, nc2nc, xyzstat, ncstat);
  
    strcpy(infilename, argv[0]);
    
    if (nc2xyz || nc2nc) strcpy(outfilename, argv[1]);
      
    if (nc2xyz) {
      if (strcmp(outfilename, "-") == 0)
	outfile = stdout;
      else 
	outfile = fopen(outfilename, "w");
      if (outfile == NULL) pe("Cannot open XYZ file %s for writing", outfilename);
    }
      
    debug("Opening NetCDF file \"%s\" for reading\n", infilename);
#ifdef HAVE_NETCDF
#ifdef NETCDF4
    if (nc3flag) {
      netcdf_check(nc_open(infilename, NC_64BIT_OFFSET | NC_NOWRITE, &nc_in));
    } else {
      netcdf_check(nc_open(infilename, NC_NOWRITE, &nc_in));
    }
#else
    netcdf_check(nc_open(infilename, NC_64BIT_OFFSET | NC_NOWRITE, &nc_in));
#endif
#else
      pe("No NetCDF support compiled in.\n");
#endif
      
    if (nc2nc) {
      debug("Opening NetCDF file \"%s\" for writing\n", outfilename);
#ifdef HAVE_NETCDF
#ifdef NETCDF4
      if (nc3flag) {
	netcdf_check(nc_set_default_format(NC_FORMAT_64BIT, NULL));
	netcdf_check(nc_create(outfilename, NC_64BIT_OFFSET | NC_CLOBBER, &nc_out));
      } else {
	netcdf_check(nc_set_default_format(NC_FORMAT_NETCDF4, NULL));
	netcdf_check(nc_create(outfilename, NC_NETCDF4 | NC_CLOBBER, &nc_out));
      }
#else
      netcdf_check(nc_create(outfilename, NC_64BIT_OFFSET | NC_CLOBBER, &nc_out));
#endif
#else
      pe("No NetCDF support compiled in.\n");
#endif
    }

    // Read NetCDF header
    nframes = read_netcdf(nc_in, &at, 0, atomlist, natomlist, 1, 0, 1, zflag, 0, 0.0);
    if (nframes == 0) pe("Error reading netcdf\n");
    parse_range(rflag, nframes, offset, gotcolon, start, stop, step, &f_start, &f_stop, &f_step);
    
    nframes = 0;
    for (i=f_start; i<f_stop; i+=f_step) nframes++;
    
    debug("start %d stop %d step %d\n", f_start, f_stop, f_step);
    
    n = 0;
    for (i = f_start; i < f_stop; i += f_step) {
      if (!read_netcdf(nc_in, &at, i, atomlist, natomlist, 0, 0, 1, zflag, 0, 0.0))
	pe("Error reading frame %d\n", i);
	
      if (!ncstat) {

	if (pflag) {
	  for (j=0; j<at.n_property; j++)
	    at.property_filter[j] = 0;
	  for (j=0; j<nprops; j++) {
	    k = atoms_find_property(&at, props[j]);
	    if (k == -1) pe("Property %s specified on command line not found", props[j]);
	    at.property_filter[k] = 1;
	    atoms_swap_properties(&at, k, j);
	  }
	}
	
	if (Pflag) {
	  for (j=0; j<at.n_param; j++)
	    at.param_filter[j] = 0;
	  for (j=0; j<nparams; j++) {
	    k = atoms_find_param(&at, params[j]);
	    if (k == -1) pe("Parameter %s specified on command line not found", params[j]);
	    at.param_filter[k] = 1;
	  }
	}

	if (nc2xyz) {
	  if (!write_xyz(outfile, &at, intformat, realformat, strformat, logformat, !pflag))
	    pe("Error writing xyz\n");
	} else {
	  if (!write_netcdf(nc_out, &at, n++, 0, shuffle, deflate, deflate_level))
	    pe("Error writing netcdf\n");
	}
	
	fprintf(stderr, "Frame %d/%d                      %c", i+1, nframes, i + f_step >= f_stop ? '\n' : '\r');
      } else {
	
	printed_stats = 0;
	linebuffer[0] = '\0';
	if (fflag) sprintf(linebuffer, "frame %d ", i+offset);
	if (nflag) sprintf(linebuffer, "%sn_atom %ld ", linebuffer, at.n_atom);
	  
	if (dflag) {
	  j = atoms_find_param(&at, dvar);
	  if (j == -1) {
	    j = atoms_find_property(&at, dvar);
	    if (j == -1) pe("Variable %s not found", dvar);
	    // It's a property. Print one atom per line
	    if (fflag || nflag) {
	      // Flush stats line
	      strcat(linebuffer,"\n");
	      printf(linebuffer);
	      printed_stats = 1;
	    }
	    print_property(&at, dvar, intformat, realformat, strformat, logformat);
	  } else {
	    // It's a param. Print one frame per line
	    sprint_param(linebuffer, &at, dvar, intformat, realformat);
	  }
	}
	if (!printed_stats && (!nflag || fflag || dflag)) {
	  strcat(linebuffer,"\n");
	  printf(linebuffer);
	}
	
      }
    }
      
#ifdef HAVE_NETCF
    netcdf_check(nc_close(nc_in));
    if (nc2xyz && strcmp(outfilename, "-") != 0) fclose(outfile);
    if (nc2nc) netcdf_check(nc_close(nc_out));
#endif
  } else {
#ifdef HAVE_NETCDF
    if (ncstat && argc != 1) 
      usage(nc2xyz, xyz2nc, xyz2xyz, nc2nc, xyzstat, ncstat);
    
    if (ncstat && (cflag || Pflag || pflag)) {
      strcpy(infilename, argv[0]);

      debug("Opening NetCDF file \"%s\" for reading\n", infilename);
#ifdef NETCDF4
      if (nc3flag) {
	netcdf_check(nc_open(infilename, NC_64BIT_OFFSET | NC_NOWRITE, &nc_in));
      } else {
	netcdf_check(nc_open(infilename, NC_NOWRITE, &nc_in));
      }
#else
      netcdf_check(nc_open(infilename, NC_64BIT_OFFSET | NC_NOWRITE, &nc_in));
#endif

      nframes = read_netcdf(nc_in, &at, 0, atomlist, natomlist, 1, 0, 1, zflag, 0, 0.0);
      if (nframes == 0) pe("Error reading netcdf");
      parse_range(rflag, nframes, offset, gotcolon, start, stop, step, &f_start, &f_stop, &f_step);
      
      if (cflag) {
	nframes = 0;
	for (i=f_start; i<f_stop; i+=f_step) nframes++;
	printf("%d frames\n",nframes);
      }
      
      if (Pflag) {
	debug("Parameters:\n");
	for (i=0; i<at.n_param; i++) {
	  if (strcasecmp(at.param_key[i], "Lattice") == 0 ||
	      strcasecmp(at.param_key[i], "Properties") == 0) continue;
	
	  printf("%s %s size %d\n", at.param_key[i], PARAM_TYPES[at.param_type[i]], at.param_size[i]);
	}
      }
      
      if (pflag) {
	debug("\nProperties:\n");
	for (i=0; i<at.n_property; i++) {
	  printf("%s %s size %d\n", at.property_name[i], PROPERTY_TYPES[at.property_type[i]], at.property_ncols[i]);
	}
      }

    netcdf_check(nc_close(nc_in));
#endif
    }
  }

  if (aflag) free(atomlist);
  atoms_free(&at);
  if (at.filter != NULL) {
    free(at.filter);
    at.filter = NULL;
  }
  
  return 0;
}


#endif /* MAIN_PROGRAM */


/********************************* FORTRAN API *********************************
 *                                                                             *
 * cioinit(), ciofree(), cioquery(), cioread() and ciowrite()	       *
 *									       *
 ********************************* FORTRAN API *********************************/

#define INPUT  0
#define OUTPUT 1
#define INOUT  2 

int cioquery(Atoms *at, int *frame) {
  if (at->format == XYZ_FORMAT) {
    if (at->xyz_in == NULL) return 0;
    return read_xyz((*at).xyz_in, at, NULL, 0, *frame, 1, 1, 0, !at->got_index, 0, NULL);
  } else if (at->format == NETCDF_FORMAT) {
    if (at->nc_in == 0) return 0;
    return read_netcdf(at->nc_in, at, *frame, NULL, 0, 1, 1, 0, 0, 0, 0.0);
  } else return 0;
}

int cioinit(Atoms **at, char *filename, int *action, int *append, int *netcdf4, 
	    int **n_frame, int **n_atom, int **n_int, int **n_real, int **n_str, int **n_logical,
	    int **n_param, int **n_property, char **property_name, int **property_type, int **property_ncols,
	    int **property_start, int **property_filter, char **param_name, int **param_type, int **param_size, char **param_value, 
	    int **param_int, double **param_real, int **param_int_a, double **param_real_a, int **param_filter, double **lattice,
	    int **got_index, int **pnetcdf4)
{
  char *p, *q;
  int i, z, *zp;
#ifdef HAVE_NETCDF
  int retval;
#endif

  *at = malloc(sizeof(Atoms));
  if (*at == NULL) return 0;
  atoms_init(*at);

  // Pass back some pointers into this Atoms structure
  if (property_name) *property_name = (char *)(**at).property_name;
  if (property_type) *property_type = (**at).property_type;
  if (property_ncols) *property_ncols = (**at).property_ncols;
  if (property_start) *property_start = (**at).property_start;
  if (property_filter) *property_filter = (**at).property_filter;
  
  if (param_name) *param_name = (char *)(**at).param_key;
  if (param_type) *param_type = (**at).param_type;
  if (param_size) *param_size = (**at).param_size;
  if (param_value) *param_value = (char *)(**at).param_value;
  if (param_int) *param_int = (**at).param_int;
  if (param_real) *param_real = (**at).param_real;
  if (param_int) *param_int_a = (int *)(**at).param_int_a;
  if (param_real) *param_real_a = (double *)(**at).param_real_a;
  if (param_filter) *param_filter = (**at).param_filter;

  if (lattice) *lattice = (double *)(**at).lattice;

  if (n_frame) *n_frame = (int *)&((**at).n_frame);
  if (n_atom) *n_atom = (int *)&((**at).n_atom);
  if (n_int) *n_int = &((**at).n_int);
  if (n_real) *n_real = &((**at).n_real);
  if (n_str) *n_str = &((**at).n_str);
  if (n_logical) *n_logical = &((**at).n_logical);
  if (n_param) *n_param = &((**at).n_param);
  if (n_property) *n_property = &((**at).n_property);
  if (got_index) *got_index = &((**at).got_index);
  if (pnetcdf4) *pnetcdf4 = &((**at).netcdf4);

  (**at).xyz_in = NULL;
  (**at).xyz_out = NULL;
  (**at).nc_in = -1;
  (**at).nc_out = -1;

  // Guess format from file extension
  if (filename[0] == '\0') {
    (**at).format = NO_FORMAT;
    return 1;
  } else {
    p = strrchr(filename, '/');
    if (p == NULL) p = filename;
    q = strchr(p, '.');
    if (q != NULL) {
      if (strstr(q, "nc") || strstr(q, "NC"))
	goto NETCDF;
      if (strstr(q, "xyz") || strstr(q, "XYZ"))
	goto XYZ;
    }
    if (strcmp(filename, "stdin") == 0) goto XYZ;
    if (strcmp(filename, "stdout") == 0) goto XYZ;
    fprintf(stderr, "unknown file format: %s\n", filename);
    return 0;
  }

  XYZ:
    if ((*action == INPUT || *action == INOUT) && strcmp(filename, "stdin") != 0) {
      **n_frame = xyz_find_frames(filename, (**at).frames, (**at).atoms);
      if (**n_frame == 0) {
	fprintf(stderr, "Error building XYZ index\n");
	return 0;
      }
      (**at).got_index = 1;
    }  else {
      **n_frame = 0;
      (**at).got_index = 0;
    }
    
    if (*action == INPUT) {
      if (strcmp(filename, "stdin") == 0) {
	(**at).xyz_in = stdin;
      } else if (strcmp(filename, "stdout") == 0) {
	fprintf(stderr, "cannot open stdout for action INPOUT");
	return 0;
      }
      else
	(**at).xyz_in = fopen(filename, "r");
      if ((**at).xyz_in == NULL) return 0;
      (**at).format = XYZ_FORMAT;
      if (strcmp(filename, "stdin") != 0) { 
	z = 0; zp = &z;  cioquery(*at, zp); 
      }
    } else if (*action == OUTPUT) {
      if (strcmp(filename, "stdin") == 0) {
	fprintf(stderr, "cannot open stdin for action OUTPUT");
	return 0;
      }
      if (strcmp(filename, "stdout") == 0)
	(**at).xyz_out = stdout;
      else {
	if (*append) 
	  (**at).xyz_out = fopen(filename, "a");
	else
	  (**at).xyz_out = fopen(filename, "w");
      }
      (**at).format = XYZ_FORMAT;
      if ((**at).xyz_out == NULL) return 0;
    } else if (*action == INOUT) {
      if (strcmp(filename, "stdin") == 0) {
	fprintf(stderr, "cannot open stdin for action INOUT");
	return 0;
      } else if (strcmp(filename, "stdout") == 0) {
	fprintf(stderr, "cannot open stdout for action INOUT");
	return 0;	
      }
      if (*append) {
	(**at).xyz_out = fopen(filename, "a+");
	(**at).xyz_in  = (**at).xyz_out;
      } else {
	(**at).xyz_out = fopen(filename, "w+");
	(**at).xyz_in  = (**at).xyz_out;
      }      
      if ((**at).xyz_out == NULL) return 0;
    (**at).format = XYZ_FORMAT;
      z = 0; zp = &z;  cioquery(*at, zp);
    } else {
      fprintf(stderr,"Bad action %d - should be one of INPUT(%d), OUTPUT(%d) or INOUT (%d)",
	      *action, INPUT, OUTPUT, INOUT);
      return 0;
    }
    return 1;
    
  NETCDF:
    (**at).got_index = 1;
    (**at).netcdf4 = *netcdf4;
#ifdef HAVE_NETCDF
    if (*action == INPUT) {

      debug("Opening NetCDF file \"%s\" for reading\n", filename);
#ifdef NETCDF4
      if ((**at).netcdf4) {
	netcdf_check(nc_open(filename, NC_NOWRITE, &((**at).nc_in)));
      } else {
	netcdf_check(nc_open(filename, NC_64BIT_OFFSET | NC_NOWRITE, &((**at).nc_in)));
      }
#else
      netcdf_check(nc_open(filename, NC_64BIT_OFFSET | NC_NOWRITE, &((**at).nc_in)));
#endif
      (**at).format = NETCDF_FORMAT;
      z = 0; zp = &z;  cioquery(*at, zp);
    } else if (*action == OUTPUT) {
      debug("Opening NetCDF file \"%s\" for writing\n", filename);
#ifdef HAVE_NETCDF
#ifdef NETCDF4
      if ((**at).netcdf4) {
	netcdf_check(nc_set_default_format(NC_FORMAT_NETCDF4, NULL));
	netcdf_check(nc_create(filename, NC_NETCDF4 | NC_CLOBBER, &((**at).nc_out)));
      } else {
	netcdf_check(nc_set_default_format(NC_FORMAT_64BIT, NULL));
	netcdf_check(nc_create(filename, NC_64BIT_OFFSET | NC_CLOBBER, &((**at).nc_out)));
      }
#else
      netcdf_check(nc_create(filename, NC_64BIT_OFFSET | NC_CLOBBER, &((**at).nc_out)));
#endif
      (**at).n_frame = 0;
      (**at).format = NETCDF_FORMAT;
#else
#endif
    } else if (*action == INOUT) {
#ifdef NETCDF4
      if ((**at).netcdf4) {
	netcdf_check(nc_open(filename, NC_WRITE, &((**at).nc_in)));
      } else {
	netcdf_check(nc_open(filename, NC_64BIT_OFFSET | NC_WRITE, &((**at).nc_in)));
      }
#else
      netcdf_check(nc_open(filename, NC_64BIT_OFFSET | NC_WRITE, &((**at).nc_in)));
#endif
      (**at).nc_out = (**at).nc_in;
      // do a query and copy all dimension and variable ids from IN to OUT
      
      (**at).format = NETCDF_FORMAT;
      z = 0; zp = &z; cioquery(*at, zp);

      for (i=0; i<(**at).n_property; i++)
	(**at).property_var_id[i][NETCDF_OUT] = (**at).property_var_id[i][NETCDF_IN];
      for (i=0; i<(**at).n_param; i++)
	(**at).param_var_id[i][NETCDF_OUT] = (**at).param_var_id[i][NETCDF_IN];

      (**at).frame_dim_id[NETCDF_OUT] = (**at).frame_dim_id[NETCDF_IN];
      (**at).spatial_dim_id[NETCDF_OUT] = (**at).spatial_dim_id[NETCDF_IN];
      (**at).atom_dim_id[NETCDF_OUT] = (**at).atom_dim_id[NETCDF_IN];
      (**at).cell_spatial_dim_id[NETCDF_OUT] = (**at).cell_spatial_dim_id[NETCDF_IN];
      (**at).cell_angular_dim_id[NETCDF_OUT] = (**at).cell_angular_dim_id[NETCDF_IN];
      (**at).label_dim_id[NETCDF_OUT] = (**at).label_dim_id[NETCDF_IN];
      (**at).string_dim_id[NETCDF_OUT] = (**at).string_dim_id[NETCDF_IN];
      (**at).cell_lengths_var_id[NETCDF_OUT] = (**at).cell_lengths_var_id[NETCDF_IN];
      (**at).cell_angles_var_id[NETCDF_OUT] = (**at).cell_angles_var_id[NETCDF_IN];
      (**at).spatial_var_id[NETCDF_OUT] = (**at).spatial_var_id[NETCDF_IN];
      (**at).cell_spatial_var_id[NETCDF_OUT] = (**at).cell_spatial_var_id[NETCDF_IN];
      (**at).cell_angular_var_id[NETCDF_OUT] = (**at).cell_angular_var_id[NETCDF_IN];
      
    } else {
      fprintf(stderr,"Bad action %d - should be one of INPUT(%d), OUTPUT(%d) or INOUT (%d)",
	      *action, INPUT, OUTPUT, INOUT);
      return 0;      
    }
    return 1;
#else
    pe("No NetCDF support compiled in");
    return 0;
#endif
}


void ciofree(Atoms *at) {
  // We don't call atoms_free as data belongs to fortran caller
  if (at->format == XYZ_FORMAT) {
    if (at->xyz_in  != NULL && at->xyz_in  != stdin)  fclose(at->xyz_in);
    if (at->xyz_in != at->xyz_out && at->xyz_out != NULL && at->xyz_out != stdout) fclose(at->xyz_out);
  } else if (at->format == NETCDF_FORMAT) {
#ifdef HAVE_NETCDF
    if (at->nc_in != -1) nc_close(at->nc_in);
    if (at->nc_out != -1 && at->nc_out != at->nc_in) nc_close(at->nc_out);
#else
    pe("No NetCDF support compiled in");
    return;
#endif
  }
  free(at->filter);
  free(at);
}

int cioread(Atoms *at, int *frame, int *int_data, double *real_data, char *str_data, 
	     int *logical_data, int *zero)
{
  int status;

  at->int_data = int_data;
  at->real_data = real_data;
  at->str_data = str_data;
  at->logical_data = logical_data;

  if (at->format == XYZ_FORMAT) {
    if (at->xyz_in == NULL) return 0;
    status = read_xyz(at->xyz_in, at, NULL, 0, *frame, 0, 0, 0, !at->got_index, 0, NULL);
    if (status == 0) return status;
    return status;
  } else if (at->format == NETCDF_FORMAT) {
    if (at->nc_in == 0) return 0;
    return read_netcdf(at->nc_in, at, *frame, NULL, 0, 0, 0, 0, *zero, 0, 0.0);
  } else return 0;
}

int ciowrite(Atoms *at, int *int_data, double *real_data, char *str_data, int *logical_data,
	      char *intformat, char *realformat, int *frame, int *shuffle, int *deflate,
	      int *deflate_level)
{
  at->int_data = int_data;
  at->real_data = real_data;
  at->str_data = str_data;
  at->logical_data = logical_data; 

  if (at->format == XYZ_FORMAT) {
    if (at->xyz_out == NULL) return 0;
    if (at->xyz_out != stdout && fseek(at->xyz_out, 0, SEEK_END) != 0) return 0;

    if (!write_xyz(at->xyz_out, at, intformat, realformat, "%.10s", "%5s", 1))
      return 0;
    if (at->xyz_out == stdout) fflush(at->xyz_out);
    return 1;
  } else if (at->format == NETCDF_FORMAT) {
    if (at->nc_out == 0) return 0;
    return write_netcdf(at->nc_out, at, *frame, 1, *shuffle, *deflate, *deflate_level);
  } else return 0;
}

int cioupdate(Atoms *at, int *int_data, double *real_data, char *str_data, int *logical_data)
{
  at->int_data = int_data;
  at->real_data = real_data;
  at->str_data = str_data;
  at->logical_data = logical_data; 

  return 1;
}
