/*******************************************/
/* Threaded Atomistic Configuration Viewer */
/*                                         */
/* Fri Mar 31 2000 Ju Li <liju99@mit.edu>  */
/*******************************************/

#ifndef _A_H
#define _A_H

#include <pthread.h>
/* default A/eV/amu unit system by not defining U_'s here */
#include <Atoms.h>
#include <AX.h>

#define XTERM_IDENTIFIER_SIZE     64
#define MAX_FILENAME_SIZE         512
/* completely loaded z-map in LII cache? */
#define DEFAULT_WIDTH             256
#define DEFAULT_HEIGHT            DEFAULT_WIDTH
#define DEFAULT_ATOM_R_RATIO      (1+GOLDEN_RATIO)
/* all bonds have the same radii: this looks best */
#define DEFAULT_BOND_RADIUS       0.25
#define DEFAULT_DBLECLICK_IN_MS   250
#define WIREFRAME_MODE_CONTRAST   0
#define WIREFRAME_MODE_NONE       1
#define WIREFRAME_MODE_RGBO       2
#define WIREFRAME_MODE_RGBK       3
#define WIREFRAME_MODE_RGB        4
#define DEFAULT_MGS_RATIO         0.9
#define COLOR_MODE_NORMAL         1
#define COLOR_MODE_COORD          2
#define COLOR_MODE_AUXILIARY      3
#define COLOR_MODE_SCRATCH        4
/* crystal shift needs to be slower for pattern to be traced */
#define XTAL_SHIFT_GEAR           0.2
#define RGB_URL JLCPL_URL"NetApp/JavaScript/Reader/rgb.html"
#define INITIAL_VIEW_RATIO        0.95
#define DEFAULT_CMAP_FNAME        "cmap.ps"
#define GEO_SHEAR_STRAIN          0
#define GEO_CENTRAL_SYMM          1
#define MAX_GEO_MEASURES          2
#define ATOM_STACK_SIZE           4

#define FILTER_PLANE_WIREFRAME_GRADES 4
typedef struct
{
    int wireframe_mode;  /* how to draw the filter plane wire-frames */
    double s0[3];        /* reference s[] on the filter plane */
    double dx_input[4];  /* original input of plane normal */
    double dx_cache[4];  /* plane normal in storage */
} Filter_Plane;

typedef struct
{
    double bgcolor[3];  /* background color r,g,b */
    int lx,ly; /* last pointer position */
    long last_button_press_time; /* for capturing double click */
    int anchor;  /* -1: hook[] as anchor, >=0: certain atom */
    int atom_stack[ATOM_STACK_SIZE];  /* keep track of selected atoms */
    double hook[3]; /* real space coordinates */
    int suppress_printout;  /* term printout */
    double delta;  /* rate of change (gear-box) */
    double mgs_radius; /* magic sphere for pointer rotation */
    double atom_r_ratio; /* compared to ATOM_RADIUS */
    int color_mode; /* NORMAL, COORD, AUXILIARY or SCRATCH */
    int last_color_mode;
    int wireframe_mode; /* draw H-box */
    int bond_mode; /* draw bonds */
    double bond_radius; /* in Angstrom */
    int bond_atom_color_need_update; /* delayed update due to atom color */
    int shell_viewer_mode; /* xv and ghostview */
    int xtal_mode; /* to allow for crystal translation and shear-encoding */
    double xtal_origin[3]; /* shifted xtal origin */
    int bond_xtal_origin_need_update; /* delayed bond update due to shift */
    int last_surface_id; /* for if the pointer is not pointing at H box */
    int mitosis; /* for spawning children thread */
    int auxiliary_idx; /* which auxiliary property should we render */
    int auxiliary_cmap; /* which colormap should we use */
    double auxiliary_threshold[CONFIG_MAX_AUXILIARY+MAX_GEO_MEASURES][2];
    int auxiliary_thresholds_saturation;
    /* 0:invisible if outside thresholds */
    int auxiliary_thresholds_rigid; /* 0: floating thresholds */
    int parallel_projection; /* parallel / perspective projections */
    int glob_advance; /* how fast the file list advances */
    Filter_Plane fp[AX_3D_MAX_FILTER_PLANE]; /* sidekick of AX_3D.fp */
    int just_activated_fp;  /* index of the fp under focus */
} Navigator;

/* parallel projection */
#define PARALLEL_AMP    1000.
#define BOND_U          0.55
/* local writable copy of MENDELEYEV */
#define ATOM_Color_R(Z) (Dmitri[Z].red)
#define ATOM_Color_G(Z) (Dmitri[Z].green)
#define ATOM_Color_B(Z) (Dmitri[Z].blue)
#define ATOM_Radius(Z)  (Dmitri[Z].charge_radius)
/* bonds are slightly tainted by what it is connected to. */
#define BOND_UPORTION   0.8
#define BONDMIX(u,c0,f0,c1,f1) ( BOND_UPORTION * (u) + \
  (1.-BOND_UPORTION) * ( (f0)*(c0) + (f1)*(c1) ) / ((f0) + (f1)) )
#define BOND_R(i,j) BONDMIX( BOND_U, B->BALL[i].r, \
  SQUARE(B->BALL[i].radius), B->BALL[j].r, SQUARE(B->BALL[j].radius) )
#define BOND_G(i,j) BONDMIX( BOND_U, B->BALL[i].g, \
  SQUARE(B->BALL[i].radius), B->BALL[j].g, SQUARE(B->BALL[j].radius) )
#define BOND_B(i,j) BONDMIX( BOND_U, B->BALL[i].b, \
  SQUARE(B->BALL[i].radius), B->BALL[j].b, SQUARE(B->BALL[j].radius) )
#define BONDCOLOR(i,j,k) AX_3D_AssignRGB( C->CYLINDER[k], \
  BOND_R(i,j), BOND_G(i,j), BOND_B(i,j) )
#define ATOM_RADIUS(Z) (MENDELEYEV[Z].charge_radius)
/* #define ATOM_RADIUS(Z) (MENDELEYEV[Z].empirical_radius) */

/* A.c: */
Aapp_Declare_Config;
extern double HI[3][3], volume, lengthscale, cm[3];
extern Chemtab ct[1];
extern Tp *tp;
extern Neighborlist N[1];
extern AX_3D_Balls B[1];
extern AX_3D_Cylinders C[1];
extern char fbasename[MAX_FILENAME_SIZE];
extern char config_fname[MAX_FILENAME_SIZE];
extern const double gearbox[10];
extern Navigator n[AX_MAXWIN];
extern char xterm_identifier[XTERM_IDENTIFIER_SIZE];
extern Window xterm_win;
extern Atom_coordination_color CoordColor[ATOM_COORDINATION_MAX+1];
extern struct Mendeleyev Dmitri[MENDELEYEV_MAX+1];

void paint_scene (int iw);
bool treatevent (int iw);
void thread_start(void *icopy);
void select_fbasename (char *raw);

/* primitives.c */
extern int temporary_disable_bond;
void Config_to_3D_Balls (double atom_r_ratio);
void Config_to_3D_Bonds (double bond_radius);
bool change_atom_r_ratio (double atom_r_ratio);
bool change_bond_radius (double bond_radius);
void hook_to_cylinder (int k, double *hook);
void reload_config (int iw, bool term_input_filename);
void script_animate (int iw);
bool set_glob_advance (int iw);
bool config_advance (int iw, int how_much);
bool config_advance_first (int iw);
bool config_advance_last (int iw);
bool load_atom_color_from_file (int iw);

/* viewport.c */
bool foo_advance (int iw, double delta);
bool pointer_advance (int iw, int to_x, int to_y);
bool translate (int iw, int i, double d);
bool pointer_translate (int iw, int to_x, int to_y);
void rotate (int iw, double R[3][3]);
void mgs (double x, double y, double r, double a[3]);
bool pointer_rotate (int iw, int to_x, int to_y);
bool axis_rotate (int iw, int i, double theta);

/* geo.c: */
typedef struct
{
    char *token;
    void (*fun) (double **);
    int should_evaluate;
    int has_evaluated;
} GeoList;
extern double *geo[MAX_GEO_MEASURES];
extern GeoList geolist[MAX_GEO_MEASURES];
extern double avg_M[3][3], avg_shear_strain, avg_coordination;
extern short *coordination, coordination_crystal;
extern int coordination_hist[ATOM_COORDINATION_MAX+1];
extern int shear_strain_subtract_mean, central_symm_neighbormax;
void evaluate_shear_strain (double **shear_strain);
bool change_shear_strain_subtract_mean (int iw);
void evaluate_central_symm (double **central_symm);
bool change_central_symm_neighbormax (int iw);
void geo_clear_has_evaluated_flags();
void geo_set_should_evaluate_flag (int i);
void evaluate_geo_measures();
void evaluate_geo_measures();

/* utils.c: */
AX_3D_Lines *plane_wireframe
(double dx[3], double d0, AX_Float r, AX_Float g, AX_Float b);
bool treat_numeral_event (int iw, int number);
bool shift_filter_plane (int iw, double delta);
bool capture_png (int iw);
bool capture_jpg (int iw);
bool capture_eps (int iw);
bool look_at_the_anchor (int iw);
bool observer_goto (int iw);
bool resize_window (int iw);
bool bgcolor_change (int iw);
void bond_atom_color_update (int iw);
bool assign_normal_color (int iw);
bool color_encode_strain (int iw, double amplification, double threshold);
bool color_encode_auxiliary (int iw);
bool load_auxiliary_from_file (int iw);
bool change_auxiliary_colormap (int iw);
void save_auxiliary_colormap (int iw, char *default_cmap_fname);
bool assign_coordination_color (int iw);
bool bond_color_change (int iw, int i, bool invisible);
bool atom_color_change (int iw, int i, bool invisible);
bool normal_color_change (int iw, int t, bool invisible);
bool assign_original_normal_color (int iw);
bool change_coordination_color (int iw, int i, bool invisible);
bool perspective_to_parallel (int iw);
bool parallel_to_perspective (int iw);
void xterm_get_focus (int iw);
void xterm_release_focus (int iw);
void reset_auxiliary_threshold (int iw, int i);
void atom_stack_insert (int iw, int atom);

/* xtal_shift.c: */
void atom_xtal_origin (double xtal_origin[3]);
void bond_xtal_origin_update (int iw);
bool xtal_shift (int iw, int i, double d);
bool pointer_grab_xtal_shift (int iw, int to_x, int to_y);
double Eyesight_Intersect_H_Surface
(double H[3][3], int surface_id,double x0[3],double V[3][3],
 double k0, double k1,double s[3], double x[3]);
int Eyesight_Intersect_H_Box (int iw, int to_x, int to_y, double xx[3]);
bool pointer_xtal_shift (int iw, int to_x, int to_y);
bool xtal_origin_goto (int iw);
bool xtal_origin_zero (int iw);

/* info.c: */
void evaluate_shear_strain();
void print_coordination_histogram();
bool print_atom (int iw, int i);
bool find_atom (int iw);
bool print_bond (int iw, int k);
bool print_status (int iw);
bool print_help (int iw);
void atom_pair (int j, int i, double dxji[4]);
void print_atom_pair_info (int iw, int j, int i);
double atom_triplet (int k, int j, int i, double dxkj[4], double dxij[4]);
void print_atom_triplet_info (int iw, int k, int j, int i);
void print_atom_quartet_info (int iw, int l, int k, int j, int i);
void save_atoms_in_monoclinic_filter (int iw);

/* rcut_patch.c: */
typedef struct
{
    int Zi,Zj;    /* atom chemical species */
    double rcut; /* nearest neighbor distance cutoff */
} Rcut_patch;
#define RCUT_PATCH_MAX  256
extern Rcut_patch rcut_patch [RCUT_PATCH_MAX];
extern int rcut_patching, rcut_patch_top, rcut_patch_item;
extern char rcut_patch_pairname[];
void start_rcut_patch(int iw);
bool apply_rcut_patch(int iw);
void finish_rcut_patch(int iw);

/* scratch.c: */
extern AX_Float *scratch_r, *scratch_g, *scratch_b;
extern int scratch_np, scratch_n[3];
bool scratch_color (int iw);
bool rescratch (int iw);
void scratch_free (int iw);

/* LeastSquareStrain.c: */
extern int ComputeLeastSquareStrain;
extern IsoAtomicReference ref[1];
extern char ref_fbasename[MAX_FILENAME_SIZE];
/* Free the auxiliary properties based on least-square strain */
void LeastSquareStrain_Free();
/* Append least-square strain as auxiliary properties */
void LeastSquareStrain_Append();

#endif
