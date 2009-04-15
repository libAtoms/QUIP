
#define XYZ_MAX_FRAMES 100000
#define LINESIZE 1024

#define MAX_ENTRY_COUNT 100
#define MAX_PROP_LENGTH 256
#define MAX_PARAM_COUNT 256

#define PROPERTY_STRING_LENGTH 10
#define PARAM_STRING_LENGTH 1024
#define PROPERTY_INT     1
#define PROPERTY_REAL    2
#define PROPERTY_STR     3
#define PROPERTY_LOGICAL 4

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

#define NETCDF_IN  0
#define NETCDF_OUT 1

#define XYZ_FORMAT 1
#define NETCDF_FORMAT 2
#define NO_FORMAT 3

typedef struct {
  int initialised;
  double lattice[3][3];
  int *int_data;
  double *real_data;
  int *logical_data;
  char *str_data;
  int *filter;
  int n_property;
  int entry_count;
  int property_type[MAX_ENTRY_COUNT];
  int property_start[MAX_ENTRY_COUNT];
  int property_ncols[MAX_ENTRY_COUNT];
  int property_var_id[MAX_ENTRY_COUNT][2];
  int property_filter[MAX_ENTRY_COUNT];
  char property_name[MAX_ENTRY_COUNT][MAX_PROP_LENGTH];
  int n_param; 
  char param_key[MAX_PARAM_COUNT][MAX_PROP_LENGTH];
  char param_value[MAX_PARAM_COUNT][PARAM_STRING_LENGTH];
  int param_type[MAX_PARAM_COUNT];
  int param_size[MAX_PARAM_COUNT];
  int param_var_id[MAX_PARAM_COUNT][2];
  int param_int[MAX_PARAM_COUNT];
  double param_real[MAX_PARAM_COUNT];
  int param_int_a[MAX_PARAM_COUNT][3];
  double param_real_a[MAX_PARAM_COUNT][3];
  int param_filter[MAX_ENTRY_COUNT];
  int n_int, n_real, n_str, n_logical;
  int frame_dim_id[2], spatial_dim_id[2], atom_dim_id[2], cell_spatial_dim_id[2],
    cell_angular_dim_id[2], label_dim_id[2], string_dim_id[2];
  int cell_lengths_var_id[2], cell_angles_var_id[2], spatial_var_id[2],
    cell_spatial_var_id[2], cell_angular_var_id[2];
  size_t n_frame, n_atom, n_label, n_string, n_atom_total;
  int nc_in, nc_out;
  FILE *xyz_in, *xyz_out;
  int format;
  long frames[XYZ_MAX_FRAMES];
  int atoms[XYZ_MAX_FRAMES];
  int got_index;
} Atoms;

// Localise in memory by property
#ifdef LOCALISE_BY_PROPERTY
#define property_int(atoms, i, j, n)     atoms->int_data[(atoms->property_start[i] + j)*atoms->n_atom + n]
#define property_real(atoms, i, j, n)    atoms->real_data[(atoms->property_start[i] + j)*atoms->n_atom + n]
#define property_str(atoms, i, j, n)     atoms->str_data[PROPERTY_STRING_LENGTH*((atoms->property_start[i] + j)*atoms->n_atom + n)]
#define property_logical(atoms, i, j, n) atoms->logical_data[(atoms->property_start[i] + j)*atoms->n_atom + n]
#else
// Localise in memory by atom (as in libAtoms). This is the default
#define property_int(atoms, i, j, n)     atoms->int_data[atoms->n_int*n + atoms->property_start[i] + j]
#define property_real(atoms, i, j, n)    atoms->real_data[atoms->n_real*n + atoms->property_start[i] + j]
#define property_str(atoms, i, j, n)     atoms->str_data[PROPERTY_STRING_LENGTH*(atoms->n_str*n + atoms->property_start[i] + j)]
#define property_logical(atoms, i, j, n) atoms->logical_data[atoms->n_logical*n + atoms->property_start[i] + j]
#endif 

#define absval(x)  ( (x) < 0 ? -(x) : (x) )


void atoms_init(Atoms *atoms);
void atoms_alloc(Atoms *atoms);
void atoms_free(Atoms *atoms);
void atoms_realloc(Atoms *atoms);
void atoms_swap_properties(Atoms *atoms, int i, int j);
int atoms_add_param(Atoms *atoms, char *key, int type, int size, int var_id, int inout);
int atoms_add_property(Atoms *atoms, char *key, int type, int ncols, int var_id, int inout);
int atoms_find_property(Atoms *atoms, char *key);
int atoms_find_param(Atoms *atoms, char *key);
void lattice_abc_to_xyz(double cell_lengths[3], double cell_angles[3], double lattice[3][3]);
void lattice_xyz_to_abc(double lattice[3][3], double cell_lengths[3], double cell_angles[3]);
int xyz_find_frames(char *fname, long *frames, int *atoms);
int read_netcdf (int ncid, Atoms *atoms, int frame, int *atomlist, int natomlist, int query, 
		 int redefine, int realloc, int replacefill, int irep, double rrep);
int write_netcdf(int ncid, Atoms *atoms, int frame, int redefine,
                 int shuffle, int deflate, int deflate_level);
int write_xyz(FILE *out, Atoms *atoms, char *int_format, char *real_format, char *str_format, char *logical_format, int swap);
int read_xyz (FILE *in, Atoms *atoms, int *atomlist, int natomlist, int frame, 
	      int query, int redefine, int realloc, int supress);
