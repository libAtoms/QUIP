/*************************************************/
/* Atoms: -llapack -lblas -lm                    */
/*        -lVecMat3 -lVecMat -lScalar -lIO       */
/*                                               */
/* Physical constants, macros, and configuration */
/* operators for atomistic simulations.          */
/*                                               */
/* Dec. 12, 1999  Ju Li <liju99@mit.edu>         */
/*************************************************/

#ifndef _Atoms_h
#define _Atoms_h

#include <IO.h>
#include <VecMat.h>
#include <VecMat3.h>
#include "NIST.h"

/* Config.c: */
/****************************************************/
/* Minimal Specification of atomistic configuration */
/****************************************************/

#define DIMENSION      3                /* pertaining to this world */
#define DIMENSION2    (DIMENSION*DIMENSION)
#define DIMENSION3     27               /* 3^DIMENSION (for enclosure) */
/* upper triangular representation of a symmetric DIMENSIONxDIMENSION matrix */
#define DIMSYM        (DIMENSION*(DIMENSION+1)/2)
/* symmetric 0..2 x 0..2 -> 0..5 */
#define SYMIDX(i,j)   ((5-MIN(i,j))*MIN(i,j)/2+MAX(i,j))
#define COMMENT        "# "             /* comment line escape sequence */
#define COMMENT_CHAR   2                /* escape sequence width */
#define COMMENT_SIZE  (COMMENT_CHAR+1)  /* escape sequence mem. size */
#define SYMBOL_CHAR    2                /* atomic symbol string width */
#define SYMBOL_SIZE   (SYMBOL_CHAR+1)   /* atomic symbol mem. size */
/* symbol string of the ith atom */
#define Symbol(a,i) ((char *)(a)+SYMBOL_SIZE*(i))
#define SYMBOL(i) ((*symbol)+SYMBOL_SIZE*(i))
#define SYM(i) Symbol(symbol,i)
#define SAFE_SYMBOL(to) { if (((to)[1]==' ')||((to)[1]==EOS)) \
  {(to)[1] = (to)[0]; (to)[0] = ' ';} (to)[2]=EOS; }
#define safe_symbol(from,to) { (to)[0]=(from)[0]; (to)[1]=(from)[1]; \
  SAFE_SYMBOL(to); }
#define COMPACT_SYMBOL(a) ((*(a)==' ')?((a)+1):(a))
#define DIFFERENT_SYMBOL(a,b) (((a)[0]!=(b)[0])||((a)[1]!=(b)[1]))
/* immobile atom support */
#define MOBILE(i)         ((*mass)[i] >  0)
#define IMMOBILE(i)       ((*mass)[i] <= 0)
#define MOBILE_ATOM(i)    ((mass)[i]  >  0)
#define IMMOBILE_ATOM(i)  ((mass)[i]  <= 0)

/************************************************************************/
/* Application (Aapp) is an external program where variables (np,H) and */
/* pointers (symbol,mass,s,s1) reside, and numerical values in whose    */
/* arrays are based on certain Application-Specific Unit System (ASUS). */
/* This library (Alib) serves as an interface between Aapp and the file */
/* system, converting the ASUS values in memory to and from files of    */
/* certain Universally Understood Unit Systems (UUUS). In order to do   */
/* that, Aapp must pass to Alib variables ulength_IN_A, umass_IN_AMU,   */
/* utime_IN_NS. On the Aapp side, those values are provided as macros   */
/* ULENGTH_IN_A, UMASS_IN_AMU, UTIME_IN_NS once the root ULENGTH_IN_M,  */
/* UMASS_IN_KG, UENERGY_IN_J are defined before including <Atoms.h>.    */
/************************************************************************/

/* Definition (namespace allocation) in application: */
#define Aapp_Define_Config int np=0; double H[DIMENSION][DIMENSION]; \
  char *symbol=NULL; double *mass=NULL; double *s=NULL; double *s1=NULL
/* For alignment or other reasons you may choose to not use the */
/* above macro as long as you define all components separately. */

/* To achieve entire-Aapp scope for multiple-source Aapp development */
#define Aapp_Declare_Config extern int np; \
  extern double H[DIMENSION][DIMENSION]; extern char *symbol; \
  extern double *mass; extern double *s; extern double *s1

/************************************************************************/
/* In addition to disk I/O, Alib also has the capabilities to carry out */
/* configuration analysis, transformation and memory allocation. It is  */
/* therefore necessary for Alib to be able to modify np, symbol, mass,  */
/* s, s1, their array contents, and the contents of H[][].              */
/************************************************************************/

/* Declaration in the argument list of Alib functions */
#define Alib_Declare_Config int *np, double H[DIMENSION][DIMENSION], \
  char **symbol, double **mass, double **s, double **s1, \
  double ulength_IN_A, double umass_IN_AMU, double utime_IN_NS
/* some Alib functions like Config_push/pop() modifies ulength_IN_A, etc. */
#define Alib_Declare_CONFIG int *np, double H[DIMENSION][DIMENSION], \
  char **symbol, double **mass, double **s, double **s1, \
  double *ulength_IN_A, double *umass_IN_AMU, double *utime_IN_NS

/* passing arguments from ordinary Alib function to ordinary Alib function */
#define Config_Alib_to_Alib np, H, symbol, mass, s, s1, \
  ulength_IN_A, umass_IN_AMU, utime_IN_NS
/* from ordinary Alib function to Config_push/pop() */
#define CONFIG_Alib_to_Alib np, H, symbol, mass, s, s1, \
  &ulength_IN_A, &umass_IN_AMU, &utime_IN_NS

/* passing arguments from Aapp to Alib function */
#define Config_Aapp_to_Alib &np, H, &symbol, &mass, &s, &s1, \
  ULENGTH_IN_A, UMASS_IN_AMU, UTIME_IN_NS

/* for use in Alib developments */
#define ulength_IN_M        (ulength_IN_A * A_IN_M)
#define uarea_IN_A2          SQUARE(ulength_IN_A)
#define A2_IN_uarea         (1./uarea_IN_A2)
#define uvolume_IN_A3        CUBE(ulength_IN_A)
#define uvolume_IN_CM3      (CUBE(ulength_IN_A) * A3_IN_CM3)
#define uvolume_IN_M3       (CUBE(ulength_IN_A) * A3_IN_M3)
#define umass_IN_G          (umass_IN_AMU * AMU_IN_G)
#define umass_IN_KG         (umass_IN_AMU * AMU_IN_KG)
#define utime_IN_S          (utime_IN_NS * NS_IN_S)
#define utime_IN_PS         (utime_IN_S * S_IN_PS)
#define utime_IN_FS         (utime_IN_S * S_IN_FS)
#define uenergy_IN_J        (umass_IN_KG * SQUARE(ulength_IN_M / utime_IN_S))
#define uenergy_IN_EV       (uenergy_IN_J * J_IN_EV)
#define uenergy_IN_KJ__MOL  (uenergy_IN_EV * EV_IN_KJ__MOL)
#define uforce_IN_N         (uenergy_IN_J / ulength_IN_M)
#define ustress_IN_PA       (uenergy_IN_J / uvolume_IN_M3)

/* fprintf to a stream but adding comment escape sequence at front */
void comment (FILE *out, char *fmt, ...);

/* Volume, number density and mass density  */
/* printouts; return mass density [g/cm^3]. */
double Config_analyze_density (Alib_Declare_Config, FILE *out);

typedef struct ConfigSpecies
{
    int first;
    int counter;
    struct ConfigSpecies *next;
} ConfigSpecies;
/* Chemical species and isotope analysis printouts; returns statistics */
ConfigSpecies *Config_analyze_species (Alib_Declare_Config, FILE *out);

/* find the nearest neighbor to atom w and return their bond length */
int Xtal_analyze_nearest_bond
(Alib_Declare_Config, int w, double *nearest_bond_IN_ulength, FILE *out);

/* Save configuration in line-based CFG ASCII file */
void Config_save (Alib_Declare_Config, bool PBC, char *fname);
#define Config_SAVE(Config_Alib_to_Alib,fname) \
  Config_save(Config_Alib_to_Alib, TRUE, fname)

/************************************************************/
/* CONFIG is more flexible than Config:                     */
/* 1. Atomic mass and chemical symbol (1st and 2nd entry in */
/*    Config) can be omitted if equal to the previous atom. */
/* 2. s1,s2,s3 (3rd,4th,5th entry in Config) can be stored  */
/*    in arbitrary precisions.                              */
/* 3. d(s1)/dt,d(s2)/dt,d(s3)/dt (6,7,8th entry in Config)  */
/*    may be omitted or stored in arbitrary precisions.     */
/* 4. New properties of arbitrary precision can be added.   */
/************************************************************/
/* #define CONFIG_MAX_AUXILIARY 16 */
#define CONFIG_MAX_AUXILIARY 48
void vCONFIG_save
(Alib_Declare_Config, bool PBC, char *fname, char *s_formats,
 char *velocity_formats, int num_auxiliary, va_list ap);
/* PBC enforced */
#define vCONFIG_SAVE(Config_Alib_to_Alib, fname, \
  s_formats, velocity_formats, num_auxiliary, ap) \
  vCONFIG_save(Config_Alib_to_Alib, TRUE, fname, s_formats, \
  velocity_formats, num_auxiliary, ap)
/* no velocity */
#define vCONFIG_NV_save(Config_Alib_to_Alib, PBC, fname, \
  s_formats, num_auxiliary, ap) \
  vCONFIG_save(Config_Alib_to_Alib, PBC, fname, s_formats, \
  NULL, num_auxiliary, ap)
/* no velocity & PBC enforced */
#define vCONFIG_NV_SAVE(Config_Alib_to_Alib, fname, \
  s_formats, num_auxiliary, ap) \
  vCONFIG_save(Config_Alib_to_Alib, TRUE, fname, s_formats, \
  NULL, num_auxiliary, ap)

void CONFIG_save
(Alib_Declare_Config, bool PBC, char *fname, char *s_formats,
 char *velocity_formats, int num_auxiliary, ...);
/* PBC enforced */
void CONFIG_SAVE
(Alib_Declare_Config, char *fname, char *s_formats,
 char *velocity_formats, int num_auxiliary, ...);
/* no velocity */
void CONFIG_NV_save
(Alib_Declare_Config, bool PBC, char *fname, char *s_formats,
 int num_auxiliary, ...);
/* no velocity & PBC enforced */
void CONFIG_NV_SAVE
(Alib_Declare_Config, char *fname, char *s_formats, int num_auxiliary, ...);

/* no auxiliary */
void CONFig_save(Alib_Declare_Config, bool PBC, char *fname,
                 char *s_formats, char *velocity_formats);
/* no auxiliary & PBC enforced */
#define CONFig_SAVE(Config_Alib_to_Alib, fname, s_formats, velocity_formats) \
  CONFig_save(Config_Alib_to_Alib, TRUE, fname, s_formats, velocity_formats)
/* no auxiliary & no velocity */
#define CONFig_NV_save(Config_Alib_to_Alib, PBC, fname, s_formats) \
  CONFig_save(Config_Alib_to_Alib, PBC, fname, s_formats, NULL)
/* no auxiliary & no velocity & PBC enforced  */
#define CONFig_NV_SAVE(Config_Alib_to_Alib, fname, s_formats) \
  CONFig_save(Config_Alib_to_Alib, TRUE, fname, s_formats, NULL)

/* Save configuration in Protein Data Bank format so it can be viewed  */
/* in external viewers: you can choose different origin for the saved  */
/* PDB configuration by nonzero new_origin_s0,s1,s2, and PBC specifies */
/* whether s[] should then be fold into [0,1)x[0,1)x[0,1). More info:  */
/* http://rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html.    */
void Config_save_as_pdb
(Alib_Declare_Config, double new_origin_s0, double new_origin_s1,
 double new_origin_s2, bool PBC, char *fname);
#define Config_save_as_PDB(Config_Alib_to_Alib,fname) \
  Config_save_as_pdb(Config_Alib_to_Alib,0,0,0,FALSE,fname)

typedef struct
{
    int NP;
    double h[DIMENSION][DIMENSION];
    char *symBOL;
    double *MASS;
    double *S;
    double *S1;
    double ULENGTH_in_A;
    double UMASS_in_AMU;
    double UTIME_in_NS;
} ConfigStack;

/* Save the current config, including allocation state, to a stack member */
void Config_push (Alib_Declare_Config, ConfigStack *cs);

/* Restore the config, including allocation state, from a stack member */
void Config_pop (ConfigStack *cs, Alib_Declare_Config);

/* Allocate new memory, then clone the configuration stored in "cs" */
void Config_clone (ConfigStack *cs, Alib_Declare_Config);

/* Save the current config, including allocation state, to a stack member */
void CONFIG_push(Alib_Declare_CONFIG, ConfigStack *cs);
#define CONFIG_PUSH(cs) CONFIG_push(CONFIG_Alib_to_Alib, cs)

/* Restore the config, including allocation state, from a stack member */
void CONFIG_pop(ConfigStack *cs, Alib_Declare_CONFIG);
#define CONFIG_POP(cs) CONFIG_pop(cs, CONFIG_Alib_to_Alib)

/* Retrieve atom i information in cs to atom j in current config */
#define Config_RETRIEVE(cs,i,Config_Alib_to_Alib,j,k) { \
  for (k=0; k<SYMBOL_SIZE; k++)  SYMBOL(j)[k] = Symbol(cs->symBOL,i)[k]; \
  (*mass)[j] = cs->MASS[i]; V3EQV( &(cs->S[DIMENSION*(i)]), \
  &((*s)[DIMENSION*(j)]) ); V3EQV( &(cs->S1[DIMENSION*(i)]), \
  &((*s1)[DIMENSION*(j)]) ); (j)++; }

/* Free memory pointers stored in a stack member */
void CONFIG_erase (ConfigStack *cs);

/* Allocate new memory pointers symbol, mass, s, s1 for current np */
void Config_alloc (Alib_Declare_Config);

/* NULL-safe release of memory pointers symbol, mass, s, s1 */
void Config_free (Alib_Declare_Config);

/* NULL-safe (re)allocation of symbol, mass, s, s1 for current np */
void Config_realloc (Alib_Declare_Config);


/** configuration transformations **/

/* particle-conservative transformations */

/* s0:=s0+s0_trans, s1:=s1+s1_trans, s2:=s2+s2_trans for all atoms */
void Config_translate
(double s0_trans, double s1_trans, double s2_trans, Alib_Declare_Config);
#define Config_Translate(s_trans,Config_Alib_to_Alib) \
  Config_translate(V3E(s_trans), Config_Alib_to_Alib)

/* PBC translate */
#define Config_TRANSLATE(s0_trans,s1_trans,s2_trans,Config_Alib_to_Alib) (\
  Config_translate(s0_trans,s1_trans,s2_trans,Config_Alib_to_Alib), \
  Config_fold_into_pbc(Config_Alib_to_Alib, 0.,0.,0.) )

/* Set the center of mass to S0,S1,S2; return the amount of shift in s[] */
double *Config_set_cm (double S0, double S1, double S2, Alib_Declare_Config);
#define Config_SET_CM(Config_Alib_to_Alib) \
  Config_set_cm(0.5,0.5,0.5,Config_Alib_to_Alib)

/* achieve x' = x*R by H' = H*R  */
void Config_rotate_via_H (Alib_Declare_Config, double R[3][3]);

/* Achieve x':=x*R by s':=s*H*R*H^-1 (you have to know what you are doing) */
void Config_rotate_via_s
(Alib_Declare_Config, double S0, double S1, double S2, double R[3][3]);

#define Config_rotate_VIA_s(Config_Alib_to_Alib,s0,R) \
  Config_rotate_via_s(Config_Alib_to_Alib, V3E(s0), R)

/* s,H -> snew,Hnew: (snew-s0) * Hnew := (s-s0) * H */
void Config_transplant
(Alib_Declare_Config, double S0, double S1, double S2, double Hnew[3][3]);

/* Fold all atoms' s0 into [S0,S0+1), s1 into [S1,S1+1), s2 into [S2,S2+1). */
void Config_fold_into_pbc (Alib_Declare_Config, double S0,double S1,double S2);
#define Config_fold_into_PBC(Config_Alib_to_Alib) \
  Config_fold_into_pbc(Config_Alib_to_Alib, 0.,0.,0.)
/* fold into [S0-0.5,S0+0.5),[S1-0.5,S1+0.5),[S2-0.5,S2+0.5) by PBC */
#define Config_centerfold(Config_Alib_to_Alib,S0,S1,S2) \
  Config_fold_into_pbc(Config_Alib_to_Alib,(S0)-0.5,(S1)-0.5,(S2)-0.5)
#define Config_CENTERFOLD(Config_Alib_to_Alib,s0) \
  Config_fold_into_pbc(Config_Alib_to_Alib,(s0)[0]-0.5,(s0)[1]-0.5,(s0)[2]-0.5)

/* Swap atoms i and j in index */
void Config_swapatom (Alib_Declare_Config, int i, int j);

/* Open up a rift in "direction" by expanding that cell */
/* edge by "ratio" while shrinking the s[] by as much.  */
void Config_openrift (Alib_Declare_Config, int direction, double ratio);
#define Config_openrifts(Config_Alib_to_Alib,ratio0,ratio1,ratio2) ( \
  Config_openrift(Config_Alib_to_Alib, 0, ratio0), \
  Config_openrift(Config_Alib_to_Alib, 1, ratio1), \
  Config_openrift(Config_Alib_to_Alib, 2, ratio2) )
#define Config_opencell(Config_Alib_to_Alib,ratio) ( \
  Config_openrift(Config_Alib_to_Alib, 0, ratio), \
  Config_openrift(Config_Alib_to_Alib, 1, ratio), \
  Config_openrift(Config_Alib_to_Alib, 2, ratio) )

/* Same as Config_openrift() but with s[]'s staying at the center */
void Config_Openrift (Alib_Declare_Config, int direction, double ratio);
#define Config_Openrifts(Config_Alib_to_Alib,ratio0,ratio1,ratio2) ( \
  Config_Openrift(Config_Alib_to_Alib, 0, ratio0), \
  Config_Openrift(Config_Alib_to_Alib, 1, ratio1), \
  Config_Openrift(Config_Alib_to_Alib, 2, ratio2) )
#define Config_Opencell(Config_Alib_to_Alib,ratio) ( \
  Config_Openrift(Config_Alib_to_Alib, 0, ratio), \
  Config_Openrift(Config_Alib_to_Alib, 1, ratio), \
  Config_Openrift(Config_Alib_to_Alib, 2, ratio) )

/* Perturb each atom in each H(i,:) direction by uniformly   */
/* distributed random (-s_amplitude[i]/2, s_amplitude[i]/2). */
void Config_perturb_s (Alib_Declare_Config, double s_amplitude[3]);

/* Same as Config_perturb_s() except center of mass does not move. */
void Config_Perturb_s (Alib_Declare_Config, double s_amplitude[3]);

/* Perturb every atom in each Cartesian direction by uniformly */
/* distributed random (-x_amplitude[i]/2, x_amplitude[i]/2).   */
void Config_perturb_x (Alib_Declare_Config, double x_amplitude[3]);

/* Same as Config_perturb_x() except center of mass does not move */
void Config_Perturb_x (Alib_Declare_Config, double x_amplitude[3]);

/* Convert tagged atoms to the kind specified   */
/* behind; return the number of converted atoms */
int Config_convert_taglist_to (char taglist[], Alib_Declare_Config,
                               char *atom_symbol, double atom_mass_IN_AMU);

#define Config_convert_taglist_TO(taglist,Config_Alib_to_Alib,Z) \
  Config_convert_taglist_to(taglist,Config_Alib_to_Alib, \
  ATOM_SYMBOL(Z),ATOM_MASS_IN_AMU(Z))

/* particle-reducing transformations */

/* Get rid of atoms with taglist[i]!=0; return the number of atoms lost */
int Config_rid_of_taglist (char taglist[], Alib_Declare_Config);
#define Config_rid_of_TAGLIST(taglist,Config_Alib_to_Alib,loss) ((loss)=\
  Config_rid_of_taglist(taglist,Config_Alib_to_Alib),free((void *)(taglist)),\
  (loss))

/* Delete atoms outside of s0[] * H + [0..1, 0..1, 0..1) * Hparallelepiped */
/* Return the number of atoms retained.                                    */
int Config_keep_atoms_in_parallelepiped
(V3 s0, M3 Hparallelepiped, Alib_Declare_Config);

/* Same as Config_keep_atoms_in_parallelepiped except */
/* H[][] is also changed to Hparallelepiped[][].      */
int Config_Keep_atoms_in_parallelepiped
(V3 s0, M3 Hparallelepiped, Alib_Declare_Config);

/* Get rid of those atoms on the blacklist; return the number of atoms lost */
int Config_rid_of_blacklist
(int how_many, int their_index[], Alib_Declare_Config);
#define Config_rid_of_BLACKLIST(how_many,their_index,Config_Alib_to_Alib,loss)\
  ((loss)=Config_rid_of_blacklist(how_many,their_index,Config_Alib_to_Alib),\
  free((void *)(their_index)),(loss))

/* Same as Config_rid_of_blacklist with a different input interface */
int Config_decimate (Alib_Declare_Config, int how_many, ...);
#define Config_DECIMATE(i,Config_Alib_to_Alib) \
  Config_decimate(Config_Alib_to_Alib,1,(int)(i))

/* Randomly get rid of prescribed amount of atoms; return "how_many" */
int Config_randomly_decimate (int how_many, Alib_Declare_Config);
#define Config_randomly_DECIMATE(proportion,Config_Alib_to_Alib) \
  Config_randomly_decimate(((int)(*np)*(proportion)),Config_Alib_to_Alib)

extern bool Config_punch_use_pbc;

/* punch a real space void |dx| <= r[A] centered at S0,S1,S2 */
int Config_punch_real_space_void
(double S0, double S1, double S2, double r_IN_A, Alib_Declare_Config);

/* punch a s-space void |ds| <= sr centered at S0,S1,S2 */
int Config_punch_s_space_void
(double S0, double S1, double S2, double sr, Alib_Declare_Config);

/* punch a real space ellipsoid dx * E[A^-2] * dx' <= 1 centered at S0,S1,S2 */
int Config_punch_real_space_ellipsoid
(double S0, double S1, double S2, double E_IN__A2[3][3], Alib_Declare_Config);

/*************************************************************************/
/* Same as Config_punch_real_space_ellipsoid() except E_IN__A2[3][3] is  */
/* deduced from one of the following inputs: 1,ax,ay,az,a_IN_A; 2,ax,ay, */
/* az,a_IN_A,bx,by,bz,b_in_A; 3,ax,ay,az,a_IN_A,bx,by,bz,b_in_A,c_in_A.  */
/* In which 1 means a fissure of halfwidth a_IN_A; 2 means an elliptical */
/* crack of radii a_IN_A, b_in_A; 3 means an true ellipsoid. If any of   */
/* the _IN_A's is <= 0, it is taken to be +infinity. And, ax,ay,az don't */
/* have to be unit vector, and bx,by,bz don't have to be perpendicular   */
/* to a - the normalization and orthogonalization are taken care of by   */
/* the subroutine. Notice that in order to produce a through fissure or  */
/* crack under PBC, one should use H[i][] as ax,ay,az, bx,by,bz, etc.    */
/*************************************************************************/
int Config_punch_real_space_Ellipsoid
(Alib_Declare_Config, double S0, double S1, double S2, int n_axis, ...);

/* Delete atoms more than halfway in real-space distance */
/* from point "in" to point "in+ds". No PBC is assumed.  */
int Config_vcut (double in_S0, double in_S1, double in_S2,
                 double ds0, double ds1, double ds2, Alib_Declare_Config);
#define Config_Vcut(s0,ds,Config_Alib_to_Alib) \
  Config_vcut(V3E(s0), V3E(ds), Config_Alib_to_Alib)

/* Cut configuration beyond (and including) (S0,S1,S2) in direction dx=ds*H. */
#define Config_VCUT(S0,S1,S2,ds0,ds1,ds2,Config_Alib_to_Alib) \
  Config_vcut((S0)-TINY*ds0,(S1)-TINY*ds1,(S2)-TINY*ds2,\
  TINY*ds0,TINY*ds1,TINY*ds2, Config_Alib_to_Alib)

/* Particle-increasing transformations */

/* Add an atom to config */
int Config_add_atom
(int index_after_insertion, char *atom_symbol, double atom_mass_IN_AMU,
 double atom_s0, double atom_s1, double atom_s2, double atom_s10_IN__NS,
 double atom_s11_IN__NS, double atom_s12_IN__NS, Alib_Declare_Config);

#define Config_add_ATOM(index_after_insertion,Z,atom_s0,atom_s1,atom_s2,\
  atom_s10_IN__NS,atom_s11_IN__NS,atom_s12_IN__NS,Config_Alib_to_Alib) \
  Config_add_atom(index_after_insertion,ATOM_SYMBOL(Z),ATOM_MASS_IN_AMU(Z),\
  atom_s0,atom_s1,atom_s2,atom_s10_IN__NS,atom_s11_IN__NS,atom_s12_IN__NS,\
  Config_Alib_to_Alib)
#define Config_ADD_ATOM(index_after_insertion,Z,atom_s0,atom_s1,atom_s2,\
  Config_Alib_to_Alib) \
  Config_add_atom(index_after_insertion,ATOM_SYMBOL(Z),ATOM_MASS_IN_AMU(Z),\
  atom_s0,atom_s1,atom_s2,0.,0.,0.,Config_Alib_to_Alib)

/* Concatenate "cs" configuration to the end of current configuration, */
/* keeping H[][] of the current configuration. "cs" will be erased.    */
int Config_cat (ConfigStack *cs, Alib_Declare_Config);

/* Make n0 x n1 x n2 stacked copies of the current configuration */
void Config_multiply (int n0, int n1, int n2, Alib_Declare_Config);


/**************************/
/* Configuration Builders */
/**************************/

/* crystal structure abstract description */
#define XTAL_MAX_NPA  24
typedef struct
{
    double npa;  /* number of atoms in a unit cell */
    double nkd;  /* number of different kinds of atoms */
    double alpha;  /* b/c angle in degrees */
    double beta;   /* a/c angle in degrees */
    double gamma;  /* a/b angle in degrees */
    double ba_ratio; /* >=0:no input; else -ba_ratio is the fallback ratio */
    double ca_ratio; /* >=0:no input; else -ca_ratio is the fallback ratio */
    /* above build H[][] - a convenience but not necessarily the unit cell */
    double h[3][3];  /* three edges of the unit cell in H[][] */
    double s[DIMENSION*XTAL_MAX_NPA];  /* atom coordinates in H[][] */
    int kind[XTAL_MAX_NPA];  /* 0 <= kind[i=0..npa-1] < nkd */
} Xtal_Abstract;

/* the following structure abstracts are provided as a public resource */
/* API: symbol0,mass0_IN_AMU, a_IN_A */
#define XTAL_fcc    0             /* 1-atom */
#define XTAL_FCC   (XTAL_fcc+1)   /* cubic 4-atom */
#define XTAL_Fcc   (XTAL_FCC+1)   /* orthorhombic 6-atom */
#define XTAL_bcc   (XTAL_Fcc+1)   /* 1-atom */
#define XTAL_BCC   (XTAL_bcc+1)   /* cubic 2-atom */
#define XTAL_Bcc   (XTAL_BCC+1)   /* orthorhombic 6-atom */
#define XTAL_bCC   (XTAL_Bcc+1)   /* another orthorhombic 6-atom */
#define XTAL_dia   (XTAL_bCC+1)   /* 2-atom */
#define XTAL_Dia   (XTAL_dia+1)   /* hexagonal 12-atom: triangle side in z */
#define XTAL_DiA   (XTAL_Dia+1)   /* hexagonal 12-atom: triangle side in x */
#define XTAL_DIA   (XTAL_DiA+1)   /* cubic 8-atom */
#define XTAL_sc    (XTAL_DIA+1)   /* 1-atom */
#define XTAL_SC    (XTAL_sc)
#define XTAL_rho   (XTAL_SC+1)    /* 1-atom rhombohedral */
#define XTAL_RHO   (XTAL_rho+1)   /* ?? orthorhombic converted */
/* API: symbol0,mass0_IN_AMU, a_IN_A,c_IN_A */
#define XTAL_hcp   (XTAL_RHO+1)   /* 2-atom */
#define XTAL_Hcp   (XTAL_hcp+1)   /* 4-atom orthorhombic */
#define XTAL_HCP   (XTAL_Hcp+1)   /* 4-atom orthorhombic converted */
#define XTAL_gra2H (XTAL_HCP+1)   /* 4-atom 2H hexagonal graphite */
#define XTAL_GRA2H (XTAL_gra2H+1) /* ?? orthorhombic converted */
#define XTAL_gra3R (XTAL_GRA2H+1) /* 6-atom 3R rhombohedral graphite */
#define XTAL_GRA3R (XTAL_gra3R+1) /* ?? orthorhombic converted */
/* http://www.phy.mtu.edu/faculty/info/jaszczak/structure.html */
#define XTAL_tet   (XTAL_GRA3R+1) /* 1-atom tetragonal */
#define XTAL_TET   XTAL_tet
/* API: symbol0,mass0_IN_AMU, a_IN_A,b_IN_A,c_IN_A */
#define XTAL_ort   (XTAL_TET+1)   /* 1-atom orthorhombic */
#define XTAL_ORT   XTAL_ort
/* API: symbol0,mass0_IN_AMU, symbol1,mass1_IN_AMU, a_IN_A */
#define XTAL_zns   (XTAL_ORT+1)   /* 2-atom Zinc-blende */
#define XTAL_ZNS   (XTAL_zns+1)   /* 8-atom cubic Zinc-blende */
#define XTAL_nacl  (XTAL_ZNS+1)   /* 2-atom */
#define XTAL_NACL  (XTAL_nacl+1)  /* 8-atom cubic */
#define XTAL_cscl  (XTAL_NACL+1)  /* 2-atom */
#define XTAL_CSCL  XTAL_cscl
/* API: symbol0,mass0_IN_AMU, symbol1,mass1_IN_AMU, a_IN_A,c_IN_A */
#define XTAL_wc    (XTAL_CSCL+1)  /* 2-atom */
#define XTAL_WC    (XTAL_wc+1)    /* 4-atom orthorhombic converted */
#define XTAL_ni3al (XTAL_WC+1)    /* cubic 4-atom */
#define XTAL_Ni3Al (XTAL_ni3al+1) /* orthorhombic 12-atom */
#define XTAL_ABSTRACTS_MAX  (XTAL_Ni3Al+1)
extern Xtal_Abstract Xtal_abstracts [XTAL_ABSTRACTS_MAX];

/******************************************************************/
/* Easy way to build a Xtal of n0 x n1 x n2 unit cells; each Xtal */
/* Abstract has its own argument list in the following fashion:   */
/* [(symbol_0,mass_0_IN_AMU),...], a_IN_A,(b_IN_A, c_IN_A); in    */
/* which only the symbols need to be verbatim. mass_i_IN_AMU, if  */
/* given a negative value, will be filled in from the periodic    */
/* table. a_IN_A, if given a negative value, will be deduced from */
/* the atomic radii. For symmetries with a=b=c, do not give b,c;  */
/* otherwise specify the free b or c or both in that order. If    */
/* given b_IN_A(c_IN_A) < 0, we will fill in the standard ratios. */
/******************************************************************/
void Config_rebuild_Xtal
(Alib_Declare_Config, int n0, int n1, int n2, Xtal_Abstract *st, ...);


/***************************************************************/
/* A general config like those loaded from PDB files has H[][] */
/* and s[], but s[] may not be in [0,1). Sometime that will    */
/* cause problem, like in constructing neighborlist. This      */
/* function redefines H[][],s[] to make s[] in [0,1).          */
/***************************************************************/
void Config_to_bounding_box_config (Alib_Declare_Config, FILE *info);

/*************************/
/* Configuration Loaders */
/*************************/

/*******************************************************************/
/* Load atomistic configuration from Protein Data Bank format. If  */
/* "info" is not NULL, we will report loading progress to "info".  */
/* If there is no CRYST1 tag in the file, the program will try to  */
/* assign H[][] as an upright bounding box to the atoms.           */
/*******************************************************************/
void Config_load_from_pdb (char *fname, FILE *info, Alib_Declare_Config);
#define Config_load_from_PDB(fname,Config_Alib_to_Alib) \
  Config_load_from_pdb(fname,stdout,Config_Alib_to_Alib)

/* possible auxiliary properties stored in CONFIG file: */
extern int CONFIG_num_auxiliary;
extern double *CONFIG_auxiliary[CONFIG_MAX_AUXILIARY];
extern char CONFIG_auxiliary_name[CONFIG_MAX_AUXILIARY][TERMSIZE];
extern char CONFIG_auxiliary_unit[CONFIG_MAX_AUXILIARY][TERMSIZE];

/* Free all auxiliary properties */
void Config_free_auxiliary();

/* Load atomistic configuration from Ju Li's CFG ASCII file */
void Config_load (char *fname, FILE *info, Alib_Declare_Config);
#define Config_LOAD(fname,Config_Alib_to_Alib) \
  Config_load(fname,stdout,Config_Alib_to_Alib)

#define CONFIG_PDB_LOADED  0
#define CONFIG_CFG_LOADED  1
/* Try to guess file format from filename suffix and   */
/* then load configuration using appropriate function. */
/* Return the file format loaded.                      */
int Config_Load (char *fname, FILE *info, Alib_Declare_Config);
#define CONFIG_LOAD(fname,Config_Alib_to_Alib) \
  Config_Load(fname,stdout,Config_Alib_to_Alib)

/* change s[] to its own image in [0,1)^3 */
void Config_TRIM (Alib_Declare_Config);

/* Create complete chemical disorder but maintaining the same kinetic energy */
void Config_chemical_randomize (Alib_Declare_Config);


/**************************************************************/
/* NIST.h can be included explicitly before including Atoms.h */
/* if one wants to use ULENGTH_IN_M,UENERGY_IN_J,UMASS_IN_KG  */
/* other than what are specified below as default (Angstrom   */
/* / eV / amu system). For example,                           */
/*                                                            */
/* code.c:                                                    */
/* ...                                                        */
/* #include <NIST.h>                                          */
/* #define ULENGTH_IN_M  CM_IN_M                              */
/* #define UMASS_IN_KG   G_IN_KG                              */
/* #define UENERGY_IN_J  G_IN_KG * CM_IN_M * CM_IN_M          */
/* #include <Atoms.h>                                         */
/* ...                                                        */
/*                                                            */
/* would in effect support CGS as the reduced unit system.    */
/**************************************************************/
/* Angstrom/eV/amu is often convenient for atomistic simulations, */
#ifndef ULENGTH_IN_M
#define ULENGTH_IN_M  A_IN_M
#endif
#ifndef UENERGY_IN_J
#define UENERGY_IN_J  EV_IN_J
#endif
#ifndef UMASS_IN_KG
#define UMASS_IN_KG  AMU_IN_KG
#endif
/* but one can override the defaults by defining them beforehand. */


/* ULENGTH_IN_M,UENERGY_IN_J,UMASS_IN_KG are set now; */
/* let us deduce the rest of the reduced unit system. */

#define M_IN_ULENGTH  (1./ULENGTH_IN_M)
#define J_IN_UENERGY  (1./UENERGY_IN_J)
#define KG_IN_UMASS   (1./UMASS_IN_KG)

#define ULENGTH_IN_A   (ULENGTH_IN_M * M_IN_A)
#define A_IN_ULENGTH    (1./ULENGTH_IN_A)
#define ULENGTH_IN_NM  (ULENGTH_IN_M * M_IN_NM)
#define NM_IN_ULENGTH   (1./ULENGTH_IN_NM)
#define ULENGTH_IN_UM  (ULENGTH_IN_M * M_IN_UM)
#define UM_IN_ULENGTH   (1./ULENGTH_IN_UM)
#define ULENGTH_IN_MM  (ULENGTH_IN_M * M_IN_MM)
#define MM_IN_ULENGTH   (1./ULENGTH_IN_MM)
#define ULENGTH_IN_CM  (ULENGTH_IN_M * M_IN_CM)
#define CM_IN_ULENGTH   (1./ULENGTH_IN_CM)

#define UAREA_IN_M2  (ULENGTH_IN_M * ULENGTH_IN_M)
#define M2_IN_UAREA   (1./UAREA_IN_M2)
#define UAREA_IN_CM2  (ULENGTH_IN_CM * ULENGTH_IN_CM)
#define CM2_IN_UAREA   (1./UAREA_IN_CM2)
#define UAREA_IN_A2  (ULENGTH_IN_A * ULENGTH_IN_A)
#define A2_IN_UAREA   (1./UAREA_IN_A2)
#define UAREA_IN_NM2  (ULENGTH_IN_NM * ULENGTH_IN_NM)
#define NM2_IN_UAREA   (1./UAREA_IN_NM2)
#define UAREA_IN_UM2  (ULENGTH_IN_UM * ULENGTH_IN_UM)
#define UM2_IN_UAREA   (1./UAREA_IN_UM2)

#define UVOLUME_IN_M3  (UAREA_IN_M2 * ULENGTH_IN_M)
#define M3_IN_UVOLUME   (1./UVOLUME_IN_M3)
#define UVOLUME_IN_A3  (UAREA_IN_A2 * ULENGTH_IN_A)
#define A3_IN_UVOLUME   (1./UVOLUME_IN_A3)
#define UVOLUME_IN_NM3  (UAREA_IN_NM2 * ULENGTH_IN_NM)
#define NM3_IN_UVOLUME   (1./UVOLUME_IN_NM3)
#define UVOLUME_IN_UM3  (UAREA_IN_UM2 * ULENGTH_IN_UM)
#define UM3_IN_UVOLUME   (1./UVOLUME_IN_UM3)
#define UVOLUME_IN_CM3  (UAREA_IN_CM2 * ULENGTH_IN_CM)
#define CM3_IN_UVOLUME   (1./UVOLUME_IN_CM3)

#define UNUMBERDENSITY_IN_MOL__M3   (1./AVO / UVOLUME_IN_M3)
#define MOL__M3_IN_UNUMBERDENSITY    (1./UNUMBERDENSITY_IN_MOL__M3)
#define UNUMBERDENSITY_IN_MOL__CM3  (1./AVO / UVOLUME_IN_CM3)
#define MOL__CM3_IN_UNUMBERDENSITY   (1./UNUMBERDENSITY_IN_MOL__CM3)

#define UENERGY_IN_EV  (UENERGY_IN_J * J_IN_EV)
#define EV_IN_UENERGY   (1./UENERGY_IN_EV)
#define UENERGY_IN_mEV  (UENERGY_IN_J * J_IN_mEV)
#define mEV_IN_UENERGY   (1./UENERGY_IN_mEV)
/* reaction heat */
#define UENERGY_IN_KJ__MOL  (UENERGY_IN_EV * EV_IN_KJ__MOL)
#define KJ__MOL_IN_UENERGY  (1./UENERGY_IN_KJ__MOL)

#define UFORCE_IN_N     (UENERGY_IN_J / ULENGTH_IN_M)
#define N_IN_UFORCE     (1./UFORCE_IN_N)
#define UFORCE_IN_PN    (UFORCE_IN_N * N_IN_PN)
#define PN_IN_UFORCE    (1./UFORCE_IN_PN)

#define UFORCE_IN_DYNE  (UFORCE_IN_N * N_IN_DYNE)
#define DYNE_IN_UFORCE  (1./UFORCE_IN_DYNE)

#define UTIME_IN_S   (ULENGTH_IN_M * sqrt(UMASS_IN_KG/UENERGY_IN_J))
#define S_IN_UTIME    (1./UTIME_IN_S)
#define UTIME_IN_FS  (UTIME_IN_S * S_IN_FS)
#define FS_IN_UTIME   (1./UTIME_IN_FS)
#define UTIME_IN_PS  (UTIME_IN_S * S_IN_PS)
#define PS_IN_UTIME   (1./UTIME_IN_PS)
#define UTIME_IN_NS  (UTIME_IN_S * S_IN_NS)
#define NS_IN_UTIME   (1./UTIME_IN_NS)
#define UTIME_IN_US  (UTIME_IN_S * S_IN_US)
#define US_IN_UTIME   (1./UTIME_IN_US)
#define UTIME_IN_MS  (UTIME_IN_S * S_IN_MS)
#define MS_IN_UTIME   (1./UTIME_IN_MS)

#define UFREQ_IN_HZ   (1./UTIME_IN_S)
#define HZ_IN_UFREQ    (1./UFREQ_IN_HZ)
#define UFREQ_IN_KHZ  (UFREQ_IN_HZ * HZ_IN_KHZ)
#define KHZ_IN_UFREQ   (1./UFREQ_IN_KHZ)
#define UFREQ_IN_MHZ  (UFREQ_IN_HZ * HZ_IN_MHZ)
#define MHZ_IN_UFREQ   (1./UFREQ_IN_MHZ)
#define UFREQ_IN_GHZ  (UFREQ_IN_HZ * HZ_IN_GHZ)
#define GHZ_IN_UFREQ   (1./UFREQ_IN_GHZ)
#define UFREQ_IN_THZ  (UFREQ_IN_HZ * HZ_IN_THZ)
#define THZ_IN_UFREQ   (1./UFREQ_IN_THZ)

#define UOMEGA_IN_R__S  (1./UTIME_IN_S)
#define R__S_IN_UOMEGA   (1./UOMEGA_IN_R__S)

#define UOMEGA_TO_HZ   (UOMEGA_IN_R__S * OMEGA_IN_R__S_TO_NU_IN_HZ)
#define HZ_TO_UOMEGA    (1./UOMEGA_TO_HZ)
#define UOMEGA_TO_MHZ  (UOMEGA_TO_HZ * HZ_IN_MHZ)
#define MHZ_TO_UOMEGA   (1./UOMEGA_TO_MHZ)
#define UOMEGA_TO_GHZ  (UOMEGA_TO_HZ * HZ_IN_GHZ)
#define GHZ_TO_UOMEGA   (1./UOMEGA_TO_GHZ)
#define UOMEGA_TO_THZ  (UOMEGA_TO_HZ * HZ_IN_THZ)
#define THZ_TO_UOMEGA   (1./UOMEGA_TO_THZ)

#define UVELOCITY_IN_M__S  (ULENGTH_IN_M / UTIME_IN_S)
#define M__S_IN_UVELOCITY  (1./UVELOCITY_IN_M__S)
#define UVELOCITY_IN_NM__PS  (ULENGTH_IN_NM / UTIME_IN_PS)
#define NM__PS_IN_UVELOCITY   (1./UVELOCITY_IN_NM__PS)

#define BOLZ_IN_UENERGY__K  (BOLZ_IN_J__K * J_IN_UENERGY)
#define MYBOLZ  BOLZ_IN_UENERGY__K
#define HBAR_IN_UENERGY_UTIME__R  (HBAR_IN_J_S__R * J_IN_UENERGY * S_IN_UTIME)
#define MYHBAR  HBAR_IN_UENERGY_UTIME__R

#define UMASS_IN_AMU  (UMASS_IN_KG * KG_IN_AMU)
#define AMU_IN_UMASS   (1./UMASS_IN_AMU)
#define UMASS_IN_G    (UMASS_IN_KG * KG_IN_G)
#define G_IN_UMASS     (1./UMASS_IN_G)

#define UMASSDENSITY_IN_KG__M3  (UMASS_IN_KG / UVOLUME_IN_M3)
#define KG__M3_IN_UMASSDENSITY   (1./UMASSDENSITY_IN_KG__M3)
#define UMASSDENSITY_IN_G__CM3  (UMASS_IN_G / UVOLUME_IN_CM3)
#define G__CM3_IN_UMASSDENSITY   (1./UMASSDENSITY_IN_G__CM3)

#define USTRESS_IN_PA   (UENERGY_IN_J / UVOLUME_IN_M3)
#define PA_IN_USTRESS    (1./USTRESS_IN_PA)
#define USTRESS_IN_MPA  (USTRESS_IN_PA * PA_IN_MPA)
#define MPA_IN_USTRESS   (1./USTRESS_IN_MPA)
#define USTRESS_IN_GPA  (USTRESS_IN_PA * PA_IN_GPA)
#define GPA_IN_USTRESS   (1./USTRESS_IN_GPA)
#define USTRESS_IN_ATM  (USTRESS_IN_PA * PA_IN_ATM)
#define ATM_IN_USTRESS   (1./USTRESS_IN_ATM)
#define USTRESS_IN_BAR  (USTRESS_IN_PA * PA_IN_BAR)
#define BAR_IN_USTRESS   (1./USTRESS_IN_BAR)
#define USTRESS_IN_MBAR (USTRESS_IN_PA * PA_IN_MBAR)
#define MBAR_IN_USTRESS  (1./USTRESS_IN_MBAR)

/* thermal conductivity: heat current = - kappa * grad(T) */
#define UKAPPA_IN_W__M_K  (UENERGY_IN_J / UTIME_IN_S / ULENGTH_IN_M)
#define W__M_K_IN_UKAPPA  (1./UKAPPA_IN_W__M_K)
#define UKAPPA_IN_W__CM_K  (UENERGY_IN_J / UTIME_IN_S / ULENGTH_IN_CM)
#define W__CM_K_IN_UKAPPA  (1./UKAPPA_IN_W__CM_K)

/* viscosity: stress = viscosity * strain rate */
#define UVISCOSITY_IN_PA_S  (USTRESS_IN_PA * UTIME_IN_S)
#define PA_S_IN_UVISCOSITY  (1./UVISCOSITY_IN_PA_S)

/* Periodic Table: */
/**********************************************************************/
/* http://www.t30.physik.tu-muenchen.de/lehrstuehle/T32/matpack/html/ */
/* Nuclear/Elements/properties.html, and http://www.webelements.com/. */
/**********************************************************************/

#define Z_H	1
#define Z_He	2
#define Z_Li	3
#define Z_Be	4
#define Z_B	5
#define Z_C	6
#define Z_N	7
#define Z_O	8
#define Z_F	9
#define Z_Ne	10
#define Z_Na	11
#define Z_Mg	12
#define Z_Al	13
#define Z_Si	14
#define Z_P	15
#define Z_S	16
#define Z_Cl	17
#define Z_Ar	18
#define Z_K	19
#define Z_Ca	20
#define Z_Sc	21
#define Z_Ti	22
#define Z_V	23
#define Z_Cr	24
#define Z_Mn	25
#define Z_Fe	26
#define Z_Co	27
#define Z_Ni	28
#define Z_Cu	29
#define Z_Zn	30
#define Z_Ga	31
#define Z_Ge	32
#define Z_As	33
#define Z_Se	34
#define Z_Br	35
#define Z_Kr	36
#define Z_Rb	37
#define Z_Sr	38
#define Z_Y	39
#define Z_Zr	40
#define Z_Nb	41
#define Z_Mo	42
#define Z_Tc	43
#define Z_Ru	44
#define Z_Rh	45
#define Z_Pd	46
#define Z_Ag	47
#define Z_Cd	48
#define Z_In	49
#define Z_Sn	50
#define Z_Sb	51
#define Z_Te	52
#define Z_I	53
#define Z_Xe	54
#define Z_Cs	55
#define Z_Ba	56
#define Z_La	57
#define Z_Ce	58
#define Z_Pr	59
#define Z_Nd	60
#define Z_Pm	61
#define Z_Sm	62
#define Z_Eu	63
#define Z_Gd	64
#define Z_Tb	65
#define Z_Dy	66
#define Z_Ho	67
#define Z_Er	68
#define Z_Tm	69
#define Z_Yb	70
#define Z_Lu	71
#define Z_Hf	72
#define Z_Ta	73
#define Z_W	74
#define Z_Re	75
#define Z_Os	76
#define Z_Ir	77
#define Z_Pt	78
#define Z_Au	79
#define Z_Hg	80
#define Z_Tl	81
#define Z_Pb	82
#define Z_Bi	83
#define Z_Po	84
#define Z_At	85
#define Z_Rn	86
#define Z_Fr	87
#define Z_Ra	88
#define Z_Ac	89
#define Z_Th	90
#define Z_Pa	91
#define Z_U	92
#define Z_Np	93
#define Z_Pu	94
#define Z_Am	95
#define Z_Cm	96
#define Z_Bk	97
#define Z_Cf	98
#define Z_Es	99
#define Z_Fm	100
#define Z_Md	101
#define Z_No	102
#define Z_Lr	103
#define Z_Rf	104
#define Z_Db	105
#define Z_Sg	106
#define Z_Bh	107
#define Z_Hs	108
#define Z_Mt	109
#define MENDELEYEV_MAX  Z_Mt

/* Atoms.c: */

/* standard atom coordination color encoding */
#define ATOM_COORDINATION_MAX 24
typedef struct
{
    char *name;
    double r,g,b;
} Atom_coordination_color;
extern const Atom_coordination_color ATOM_COORDINATION_COLOR
[ATOM_COORDINATION_MAX+1];

struct Mendeleyev
{
    char *symbol;      /* "Ar" */
    char *name;        /* "Argon" */
    
    int Z;             /* 18 */
    double A;          /* 39.948 [amu] */
    int period;        /* 3 (row) */
    int IUPAC;         /* group 18 (column) */

    char *old_IUPAC;        /* "VIIIB", or European group name */
    char *CAS;              /* "VIIIA", or American group name */
    char *classification;   /* "noble gas" */
    char *groupmates;       /* "He,Ne,Ar,Kr,Xe,Rn" */

    char *electron_config;             /* "3s2 3p6" */
    double first_ionization_energy;    /* 15.759 [eV] */

    double electronegativity_Allred;   /* NVL: not available */
    double electronegativity_Pauling;  /* NVL: not available */
    double electronegativity_Pearson;  /* 7.7 [eV] */

    char *aggregate_state;  /* "gas" at 20C and 1 atm */
    double mass_density;    /* 0.00166 [g/cm^3] */
    double melting_point;   /* 83.78 [K] */
    double boiling_point;   /* 87.29 [K] */

    char *space_group;        /* "Fm-3m" */
    int space_group_index;    /*  225 */
    char *lattice_structure;  /* "fcc" */
    double a;                 /* 5.256 [A] */
    double b;                 /* 5.256 [A] */
    double c;                 /* 5.256 [A] */
    double alpha;             /* 90 [degrees] */
    double beta;              /* 90 [degrees] */
    double gamma;             /* 90 [degrees] */

    double empirical_radius;  /* 1.8597 [A] */
    double charge_radius;     /* 1.9 [A] */
    double red;               /* 0.4 [out of 1] */
    double green;             /* 0.4 [out of 1] */
    double blue;              /* 0.4 [out of 1] */
    
};
/* empirical atomic radii: http://www.webelements.com/ */
/* charge radii: http://www.fhi-berlin.mpg.de/th/balsac/balm.47.html */
/* colors: http://www.chemicalgraphics.com/paul/Manual.html#CPK-Atom-Colors */
/* http://www.chemicalgraphics.com/paul/Distribution.html periodic.tab */

extern const struct Mendeleyev MENDELEYEV [MENDELEYEV_MAX+1];

/* conversion from atomic number Z to all other properties */
#define ATOM_SYMBOL(Z)           (MENDELEYEV[Z].symbol)
#define atom_symbol(Z)           blank_advance(MENDELEYEV[Z].symbol)
#define ATOM_NAME(Z)             (MENDELEYEV[Z].name)

#define ATOM_MASS_IN_AMU(Z)      (MENDELEYEV[Z].A)
#define ATOM_MASS_IN_KG(Z)       (ATOM_MASS_IN_AMU(Z)*AMU_IN_KG)
#define ATOM_MASS_IN_G(Z)        (ATOM_MASS_IN_AMU(Z)*AMU_IN_G)
#define ATOM_MASS_IN_UMASS(Z)    (ATOM_MASS_IN_AMU(Z)*AMU_IN_UMASS)
#define ATOM_PERIOD(Z)           (MENDELEYEV[Z].period)
#define ATOM_IUPAC(Z)            (MENDELEYEV[Z].IUPAC)

#define ATOM_OLD_IUPAC(Z)        (MENDELEYEV[Z].old_IUPAC)
#define ATOM_CAS(Z)              (MENDELEYEV[Z].CAS)
#define ATOM_CLASSIFICATION(Z)   (MENDELEYEV[Z].classification)
#define ATOM_GROUPMATES(Z)       (MENDELEYEV[Z].groupmates)

#define ATOM_ELECTRON_CONFIG(Z)  (MENDELEYEV[Z].electron_config)
#define ATOM_FIRST_IONIZATION_ENERGY_IN_EV(Z) \
  (MENDELEYEV[Z].first_ionization_energy)
#define ATOM_FIRST_IONIZATION_ENERGY_IN_UENERGY(Z) \
  (ATOM_FIRST_IONIZATION_ENERGY_IN_EV(Z) * EV_IN_UENERGY)

#define ATOM_ELECTRONEGATIVITY_ALLRED(Z) \
  (MENDELEYEV[Z].electronegativity_Allred)
#define ATOM_ELECTRONEGATIVITY_PAULING(Z) \
  (MENDELEYEV[Z].electronegativity_Pauling)
#define ATOM_ELECTRONEGATIVITY_PEARSON_IN_EV(Z) \
  (MENDELEYEV[Z].electronegativity_Pearson)

/* SSC stands for standard state condition: 20C (298K) and 1 atm */
#define ATOM_SSC_AGGREGATE_STATE(Z) (MENDELEYEV[Z].aggregate_state)
#define ATOM_SSC_MASS_DENSITY_IN_G__CM3(Z) (MENDELEYEV[Z].mass_density)
#define ATOM_SSC_MASS_DENSITY_IN_KG__M3(Z) \
  (ATOM_SSC_MASS_DENSITY_IN_G__CM3(Z) * G__CM3_IN_KG__M3)
#define ATOM_SSC_MASS_DENSITY_IN_UMASSDENSITY(Z) \
  (ATOM_SSC_MASS_DENSITY_IN_G__CM3(Z) * G__CM3_IN_UMASSDENSITY)

#define ATOM_SSC_ATOMIC_VOLUME_IN_A3(Z) ( \
  ATOM_MASS_IN_KG(Z) / ATOM_SSC_MASS_DENSITY_IN_KG__M3 * M3_IN_A3 )

#define ATOM_MELTING_POINT_IN_K(Z) (MENDELEYEV[Z].melting_point)
#define ATOM_BOILING_POINT_IN_K(Z) (MENDELEYEV[Z].boiling_point)

#define ATOM_EQLATTICE_SPACE_GROUP(Z)       (MENDELEYEV[Z].space_group)
#define ATOM_EQLATTICE_SPACE_GROUP_INDEX(Z) (MENDELEYEV[Z].space_group_index)
#define ATOM_EQLATTICE_STRUCTURE(Z)         (MENDELEYEV[Z].lattice_structure)
#define ATOM_EQLATTICE_a_IN_A(Z)      (MENDELEYEV[Z].a)
#define ATOM_EQLATTICE_b_IN_A(Z)      (MENDELEYEV[Z].b)
#define ATOM_EQLATTICE_c_IN_A(Z)      (MENDELEYEV[Z].c)
#define ATOM_EQLATTICE_alpha_IN_D(Z)  (MENDELEYEV[Z].alpha)
#define ATOM_EQLATTICE_beta_IN_D(Z)   (MENDELEYEV[Z].beta)
#define ATOM_EQLATTICE_gamma_IN_D(Z)  (MENDELEYEV[Z].gamma)
#define ATOM_EQLATTICE_a_IN_ULENGTH   (ATOM_EQLATTICE_a_IN_A(Z) * A_IN_ULENGTH)
#define ATOM_EQLATTICE_b_IN_ULENGTH   (ATOM_EQLATTICE_b_IN_A(Z) * A_IN_ULENGTH)
#define ATOM_EQLATTICE_c_IN_ULENGTH   (ATOM_EQLATTICE_c_IN_A(Z) * A_IN_ULENGTH)
#define ATOM_EQLATTICE_alpha_IN_R(Z)  (ATOM_EQLATTICE_alpha_IN_D(Z) * D_IN_R)
#define ATOM_EQLATTICE_beta_IN_R(Z)   (ATOM_EQLATTICE_beta_IN_D(Z)  * D_IN_R)
#define ATOM_EQLATTICE_gamma_IN_R(Z)  (ATOM_EQLATTICE_gamma_IN_D(Z) * D_IN_R)

#define ATOM_EMPIRICAL_RADIUS_IN_A(Z) (MENDELEYEV[Z].empirical_radius)
#define ATOM_CHARGE_RADIUS_IN_A(Z)    (MENDELEYEV[Z].charge_radius)
#define ATOM_COLOR_R(Z)               (MENDELEYEV[Z].red)
#define ATOM_COLOR_G(Z)               (MENDELEYEV[Z].green)
#define ATOM_COLOR_B(Z)               (MENDELEYEV[Z].blue)
#define ATOM_EMPIRICAL_RADIUS_IN_ULENGTH(Z) \
  (ATOM_EMPIRICAL_RADIUS_IN_A(Z) * A_IN_ULENGTH)
#define ATOM_CHARGE_RADIUS_IN_ULENGTH(Z) \
  (ATOM_CHARGE_RADIUS_IN_A(Z) * A_IN_ULENGTH)

/* From our periodic table, find the atom Z corresponding */
/* to the "symbol" string. If not found, return 0.        */
int search_atom_by_symbol (char *symbol);

/* From our periodic table, find the atom Z corresponding to */
/* the "symbol" string. If not found, print error and exit.  */
int Search_atom_by_symbol (char *symbol);

typedef struct
{
    int t;   /* number of different chemical species (symbols) */
    int Z [MENDELEYEV_MAX+1];  /* Z=0 if not known element */
    int count [MENDELEYEV_MAX+1]; /* number of occurrences */
    int first [MENDELEYEV_MAX+1]; /* index of first occurrence if count>0 */
} Chemtab;

typedef char Tp;

/************************************************************************/
/* Given an atomistic configuration, (re)bind sequential chemical index */
/* to each atom, which points to an entry of Z (and Zcount) in "ct".    */
/* When allocating the index, give priority to those that appear in     */
/* "specification" string, for instance "Si C" makes Si=0,C=1, " C O"   */
/* makes C=0,O=1, etc. If specification=NULL, treat it as "". Return    */
/* the total number of chemical species found and print a report.       */
/************************************************************************/
int rebind_ct (Alib_Declare_Config, char *specification,
               Chemtab *ct, Tp **tp, FILE *out);
#define rebind_CT(Config_Alib_to_Alib,specification,ct,tp) \
  rebind_ct(Config_Alib_to_Alib,specification,ct,tp,stdout)

/* Assign new symbol[] and mass[] according to new tp[] */
void rematch_ct (Chemtab *ct, Tp **tp, Alib_Declare_Config);
/* Assign new symbol[] and mass[] according to new tp[] */
/* and then rebind_ct() according to "specification".   */
#define REMATCH_ct(Config_Alib_to_Alib,specification,ct,tp,out) ( \
  rematch_ct(ct,tp,Config_Alib_to_Alib), \
  rebind_ct(Config_Alib_to_Alib,specification,ct,tp,out) )

/* Re-index the atoms according to their assigned chemical      */
/* index and mass. This is useful for efficient CONFIG storage. */
void Config_compact_index (Chemtab *ct, Tp **tp, Alib_Declare_Config);


/* Neighborlist.c: */

/*************************************************************/
/* bin-bin/bin-atom/atom-atom lists assuming initial [0,1)   */
/* s-bounds under H[][], which does not have to be PBC. If   */
/* max_strain_eigenvalue and max_atom_displacement are both  */
/* zero, it is taken to mean a disposable list and all extra */
/* memory is freed except a final compressed atom-atom list. */
/* Otherwise one registers anchor sa[] and H0[][], frees no  */
/* memory, and provides an uncompressed atom-atom list: at   */
/* each new configuration, the user should call Neighborlist */
/* Maintenance. Under PBC, the user himself should NOT fold  */
/* s[] into [0,1); Maintenance will do that only at the atom */
/* re_anchoring event, as the ANCHORS must be in [0,1). If a */
/* direction is not PBC, and a s[i] exceeds [0,1), it would  */
/* be stuck at the boundary and a warning is issued. The     */
/* atom-atom list is guaranteed to list all pairs <= RCUT    */
/* by making sure that at any atom's re_anchoring event,     */
/* all surrounding anchors <= RLIST are recorded.            */
/*                                                           */
/* Of all libAtoms modules, only this one require [0,1) s-   */
/* bounds. load_from_pdb(), for example, would not assign    */
/* bounding box H[][] if CRYST1 tag exists; yet a lot of PDB */
/* files with CRYST1 tag do exceed [0,1) bounds, and         */
/* load_from_pdb() does not trim them. Same for CFG format,  */
/* although it is recommended that a config which is not     */
/* meant to be non-PBC call Config_fold_into_PBC before it's */
/* saved. Therefore there are three options for this module  */
/* when it is given a non-[0,1) config: QUIT, FOLD_INTO_PBC, */
/* BOUNDING_BOX. There's no guarantee of the correct action, */
/* but a good guess would be that if the config is from PDB, */
/* BOUNDING_BOX should be used; otherwise FOLD_INTO_PBC.     */
/*************************************************************/

/* Under pairwise interaction, only atom i owns i-j bond */
#define own_pair(i,j) \
  ((((i)>(j))&&(((i)+(j))&1))||(((i)<(j))&&(!(((i)+(j))&1))))
/* for estimating nearest neighbor distances */
#define ATOM_RADIUS_IN_A(Z)  ATOM_EMPIRICAL_RADIUS_IN_A(Z)

typedef struct
{   /** form: **/
    int s_overflow_err_handler; /* QUIT, FOLD, or BOUNDING_BOX */
    int small_cell_err_handler; /* QUIT, or MULTIPLY */
    bool keepcounter;  /* if FALSE, then clear all counters */
    bool pairwise; /* whether to use 1/2 saving with own_pair(i,j) */
    bool pbc [DIMENSION];  /* whether to use PBC in that s-direction */
    bool track_dx;  /* such as for monitoring self-diffusion */
    double min_shrinkage;  /* (0,1]: for bin remeshing */
    double max_atom_displacement; /* >= 0:universal tether in reduced length */
    double *rcut;  /* cutoff radii in reduced length: t x t */
    int *neighbormax;
    /* likely maximum number of neighbors (BEFORE pairwise reduction): t */
    /** public **/
    int *idx, *list;  /* atom-atom list */
    double *dx, *dxa; /* counter: dx is atom drift, dxa is anchor drift */
    /** private **/
    double *rlist;  /* list radii in reduced length: t x t */
    double *rlist2;  /* list radii in reduced length: t x t */
    int nbin[DIMENSION+1];
    int *binbindx, *binbinlist;  /* bin-bin list */
    int *bindx, *binlist, *mybin;  /* bin-atom/atom-bin list */
    double H0[3][3];  /* at the start of a strain session */
    int remesh_counter;  /* counter: number of bin-remeshing */
    double *sa, *maxtether2;
    int stack_max, *stack;
    int reanchor_counter;  /* counter: number of re_anchoring events */
    int maintenance_counter;  /* counter: number of maintenance calls */
} Neighborlist;

/* s_overflow_err_handler field */
#define NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_QUIT             0
#define NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_FOLD_INTO_PBC    1
#define NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_BOUNDING_BOX     2
#define NEIGHBORLIST_DEF_S_OVERFLOW_ERR_HANDLER \
  NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_FOLD_INTO_PBC

/* small_cell_err_handler field */
#define NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_QUIT             0
#define NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_MULTIPLY         1
#define NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_NOCHECK          2
#define NEIGHBORLIST_DEF_SMALL_CELL_ERR_HANDLER \
  NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_QUIT

#define NEIGHBORLIST_DEF_PAIRWISE     TRUE
#define NEIGHBORLIST_DEF_PBC          TRUE
#define NEIGHBORLIST_DEF_TRACK_DX     FALSE
#define NEIGHBORLIST_DEF_MIN_SHRINKAGE          1.  /* not deformable */
#define NEIGHBORLIST_DEF_MAX_ATOM_DISPLACEMENT  0.  /* not movable */
/* Assuming the shell we keep and the shell we ignore has */
#define NEIGHBORLIST_RCUT_RATIO       ((1+SQRT2)/2)
/* distance ratio sqrt(2); for fcc this is 1st and 2nd    */
/* shell; for bcc this is 2nd and 3rd shell.              */
#define NEIGHBORLIST_BINATOM_RATIO    3.0   /* density inhomogeneity */
#define NEIGHBORLIST_STACK_MAX_RATIO  32    /* workspace for a single atom */
#define NEIGHBORLIST_BINDEX(N,i,j,k) \
  ( ((i) * (N)->nbin[1] + (j)) * (N)->nbin[2] + (k) )
/* only when s[] is strictly [0,1) */
#define NEIGHBORLIST_BINDEXS(N,s,i) NEIGHBORLIST_BINDEX( N, \
  INT((s)[DIMENSION*(i)]*(N)->nbin[0]), \
  INT((s)[DIMENSION*(i)+1]*(N)->nbin[1]), \
  INT((s)[DIMENSION*(i)+2]*(N)->nbin[2]) )
#define NEIGHBOR_TABLE(rlist,ct,i,j)  ((rlist)[(i) * ((ct)->t) + (j)])
#define NEIGHBORLIST_MAXBIN 50

/* Recreate form space and fill in defaults */
void Neighborlist_Recreate_Form
(Alib_Declare_Config, Chemtab *ct, Neighborlist *N);

/* Recreate neighborlist according to submitted form */
void Neighborlist_Recreate
(Alib_Declare_Config, FILE *info, Chemtab *ct, Tp **tp, Neighborlist *N);

#define Neighborlist_RECREATE(Config_Alib_to_Alib, ct, tp, N) \
  Neighborlist_Recreate(Config_Alib_to_Alib, stdout, ct, tp, N)

/* all-default neighborlist */
#define NEIGHBORLIST_recreate(Config_Alib_to_Alib, info, ct, tp, N) { \
  Neighborlist_Recreate_Form (Config_Alib_to_Alib, ct, N); \
  Neighborlist_Recreate (Config_Alib_to_Alib, info, ct, tp, N); }
#define NEIGHBORLIST_RECREATE(Config_Alib_to_Alib, ct, tp, N) { \
  Neighborlist_Recreate_Form (Config_Alib_to_Alib, ct, N); \
  Neighborlist_Recreate (Config_Alib_to_Alib, stdout, ct, tp, N); }

/* Warn and correct (to 0 or 1-TINY) if < 0 or >= 1 */
bool warn_and_correct_if_stuck_on_the_wall
(int np, double *s, int dim, FILE *info);

void Neighborlist_Maintenance
(Alib_Declare_Config, FILE *info, Chemtab *ct, Tp **tp, Neighborlist *N);

/* Self-check the neighborlist using brute force */
void Neighborlist_Check
(Alib_Declare_Config, FILE *info, Chemtab *ct, Tp **tp, Neighborlist *N);

/* Convert a neighborlist with pairwise saving */
/* to full (redundant) list image.             */
void Neighborlist_create_nonpairwise_image
(Alib_Declare_Config, Chemtab *ct, Tp **tp, Neighborlist *N, Neighborlist *M);

/* Convert a neighborlist with pairwise saving */
/* to full (redundant) & compressed list image */
void Neighborlist_create_nonpairwise_compressed_image
(Alib_Declare_Config, Chemtab *ct, Tp **tp, Neighborlist *N, Neighborlist *M);

/* Destroy created list image */
void Neighborlist_free_nonpairwise_image (Neighborlist *M);

/* Maximum number of records in the neighborlist */
int Neighborlist_max_records (Alib_Declare_Config, Neighborlist *N);

/* Maximum number of neighbors in the neighborlist */
int Neighborlist_max_neighbors (Alib_Declare_Config, Neighborlist *N);

/* Print appropriate allocation and occupation statistics */
void Neighborlist_print_statistics (int np, Neighborlist *N, FILE *info);
#define Neighborlist_PRINT_statistics(np,N) \
  Neighborlist_print_statistics(np,N,stdout)

/* Free all memory allocations and set NULL */
void Neighborlist_Free(Neighborlist *N);


/* Gr.c: */

/************************************************************/
/* Compute g(r)'s for system with multiple chemical species */
/************************************************************/

typedef struct
{   /** form: **/
    bool pbc[DIMENSION];    /* whether to use PBC in that s-direction */
    double *rcut;           /* cutoff radius in reduced length: t x t */
    int *mesh;              /* number of mesh points: t x t */
    double *accuracy;       /* approximate relative accuracy: t x t */
    char *fn_matlab;        /* MATLAB output filename */
    int autosaves;          /* number of automatic saves */
    /** kind of public **/
    FILE *info;             /* message channel */
    int save_freq;          /* should save per how many instances */
    int instances;          /* number of instances up to now */
    double *workspace;
    double **count;
    double *normalization;  /* how much each mesh bin SHOULD give by now */
    /** private **/
    double *rcut2;          /* cutoff radius (reduced) squared: t x t */
    double *s;              /* coordinates after PBC treatment */
    int nbin[DIMENSION+1];
    int *binbindx, *binbinlist;    /* bin-bin list */
    int *bindx, *binlist, *mybin;  /* bin-atom/atom-bin list */
} Gr;

/* form defaults: */
#define GR_DEF_PBC          TRUE
#define GR_DEF_RCUT_RATIO   4.0     /* nearest neighbor defaults */
#define GR_DEF_MESH         150
#define GR_DEF_ACCURACY     0.0025
#define GR_DEF_FN_MATLAB    "gr.m"
#define GR_DEF_AUTOSAVES    20
/* private parameters: */
#define GR_BINATOM_RATIO    3.0     /* density inhomogeneity */
#define GR_TABLE(a,ct,i,j)  ((a)[(i)*((ct)->t)+(j)])
#define GR_BINDEX(GR,i,j,k) \
  ( ((i) * (GR)->nbin[1] + (j)) * (GR)->nbin[2] + (k) )
/* only when s[] is strictly [0,1) */
#define GR_BINDEXS(GR,s,i) GR_BINDEX( GR, \
  INT((s)[DIMENSION*(i)]*(GR)->nbin[0]), \
  INT((s)[DIMENSION*(i)+1]*(GR)->nbin[1]), \
  INT((s)[DIMENSION*(i)+2]*(GR)->nbin[2]) )

/* Recreate form space and fill in defaults */
void Gr_Recreate_Form (Alib_Declare_Config, Chemtab *ct, Gr *GR);

/* Reallocate and clear counters: return the necessary number  */
/* of INDEPENDENT instances to achieve the desired accuracies. */
int Gr_Reset (Alib_Declare_Config, Chemtab *ct, Gr *GR, FILE *info);

/* Accumulate a configuration instance to counters */
void Gr_Accumulate (Alib_Declare_Config, Chemtab *ct, Tp *tp, Gr *GR);

/* Save result in counters to GR->fn_matlab */
void Gr_Save (Alib_Declare_Config, Chemtab *ct, Gr *GR, FILE *info);

/* Free all memory allocations and set NULL */
void Gr_Free (Gr *GR);


/***************************************************/
/* Frequently Used Macros in Atomistic Simulations */
/***************************************************/

/* user-supplied routine */
typedef double (*StaticPotential)
    (int, M3, double *, Tp *, Neighborlist *, double *, M3);
/* np, H, s, tp, N, f, stress; return total potential energy */

#define SYMMAT_ZERO(stress) ( (stress)[0][0]=0, (stress)[1][1]=0, \
  (stress)[2][2]=0, (stress)[0][1]=0, (stress)[0][2]=0, (stress)[1][2]=0 )

#define SYMMAT_accumulate(force,magnitude,direction,stress) ( \
  (stress)[0][0] += (force)[0] * (magnitude) * (direction)[0], \
  (stress)[1][1] += (force)[1] * (magnitude) * (direction)[1], \
  (stress)[2][2] += (force)[2] * (magnitude) * (direction)[2], \
  (stress)[0][1] += (force)[0] * (magnitude) * (direction)[1], \
  (stress)[0][2] += (force)[0] * (magnitude) * (direction)[2], \
  (stress)[1][2] += (force)[1] * (magnitude) * (direction)[2] )

#define SYMMAT_ACCUMULATE(force,displacement,stress) ( \
  (stress)[0][0] += (force)[0] * (displacement)[0], \
  (stress)[1][1] += (force)[1] * (displacement)[1], \
  (stress)[2][2] += (force)[2] * (displacement)[2], \
  (stress)[0][1] += (force)[0] * (displacement)[1], \
  (stress)[0][2] += (force)[0] * (displacement)[2], \
  (stress)[1][2] += (force)[1] * (displacement)[2] )

#define SYMMAT_accomplish(stress) ( (stress)[1][0] = (stress)[0][1], \
  (stress)[2][0] = (stress)[0][2], (stress)[2][1] = (stress)[1][2] )

#define SYMMAT_Accomplish(volume,stress) ( \
  SYMMAT_accomplish(stress), M3DividE(stress,-volume) )

#define SYMMAT_ACCOMPLISH(H,volume,stress) ( \
  (volume)=M3VOLUME(H), SYMMAT_Accomplish(volume,stress) )

#define SYMMAT_DOWNLOAD(x,S) ( S[0][0]=(x)[0], S[0][1]=S[1][0]=(x)[1], \
  S[0][2]=S[2][0]=(x)[2], S[1][1]=(x)[3], S[1][2]=S[2][1]=(x)[4], \
  S[2][2]=(x)[5] )

#define SYMMAT_UPLOAD(S,x) ( (x)[0]=S[0][0], (x)[1]=S[0][1], \
  (x)[2]=S[0][2], (x)[3]=S[1][1], (x)[4]=S[1][2], (x)[5]=S[2][2] )

#define SYMMAT_MULUPLOAD(multiplier,S,x) ( \
  (x)[0]=(multiplier)*S[0][0], (x)[1]=(multiplier)*S[0][1], \
  (x)[2]=(multiplier)*S[0][2], (x)[3]=(multiplier)*S[1][1], \
  (x)[4]=(multiplier)*S[1][2], (x)[5]=(multiplier)*S[2][2] )


/* Motion.c: */
/***************************************/
/* Mobility, Velocity & Kinetic Energy */
/***************************************/

typedef struct
{
    int nmobile;        /* number of mobile atoms */
    double mobilemass;  /* total mass of mobile atoms (reduced) */
    V3 s1_avg;          /* average s drift speed (reduced) */
    double kine;        /* total kinetic energy (reduced) */
    double T_in_K;      /* system averaged temperature */
} ConfigMotion;

/* Average drift, total kinetic energy, system T & printouts */
ConfigMotion *Config_analyze_motion (Alib_Declare_Config, FILE *out);

/* Set average mobile s1 drift to zero */
ConfigMotion *Config_motion_zero_drift (ConfigMotion *m, Alib_Declare_Config);

/* Multiply the current kinetic energy and temperature by "factor" */
ConfigMotion *Config_motion_scale_T
(double factor, ConfigMotion *m, Alib_Declare_Config);

/* Assign random thermal velocities such that temperature=T_in_K */
ConfigMotion *Config_motion_initialize (double T_in_K, Alib_Declare_Config);

/* "m" & the kinetic contribution to unsymmetrized/unnormalized stress */
void Config_motion_stress (Alib_Declare_Config, ConfigMotion *m, M3 stress);

/* Compute "m", "v", and kinetic contribution to "stress" & "ep" */
void Config_motion_stress_v_ep
(Alib_Declare_Config, ConfigMotion *m, M3 stress, double *v, double *ep);

/* TRUE if all atoms in the configuration are mobile */
bool Config_all_mobile (Alib_Declare_Config);


/* Voronoi.c: */
/*************************************************/
/* Build Grains by Voronoi Site-Rotation and Cut */
/*************************************************/

typedef struct
{
    V3 s;  /* site of rotation */
    M3 R;  /* rotation matrix in real space */
} VoronoiSite;

/* those within ratio * empirical bond length */
#define VORONOISPEC_DEF_MIN_BOND_RATIO  (1. / NEIGHBORLIST_RCUT_RATIO)

typedef struct
{
    int nv; /* number of Voronoi sites */
    VoronoiSite *v;
    bool cut_against_all_images; /* or only the nearest image */
    /* If two atoms get too close, a random one gets deleted: */
    double min_bond_ratio;   /* 0: default; < 0: no deletion */
    /* bond(i,j) = min_bond_ratio * (empirical radius i + radius j) */
} VoronoiSpec;

/************************************************************/
/* Build n0 x n1 x n2 random rotation Voronoi site lattice. */
/* If st == NULL, then a random position lattice is built.  */
/************************************************************/
VoronoiSpec *VoronoiSpec_rebuild
(int n0, int n1, int n2, Xtal_Abstract *st, VoronoiSpec *V);
#define VoronoiSpec_BCC(nc0,nc1,nc2,V) \
  VoronoiSpec_rebuild(nc0,nc1,nc2,Xtal_abstracts+XTAL_BCC,V)
#define BCC_VoronoiSpec(nc,V) VoronoiSpec_BCC((nc)[0],(nc)[1],(nc)[2],V)
#define VoronoiSpec_FCC(nc0,nc1,nc2,V) \
  VoronoiSpec_rebuild(nc0,nc1,nc2,Xtal_abstracts+XTAL_FCC,V)
#define FCC_VoronoiSpec(nc,V) VoronoiSpec_FCC((nc)[0],(nc)[1],(nc)[2],V)
#define Random_VoronoiSpec(n,V) VoronoiSpec_rebuild(n,1,1,NULL,V)
#define VoronoiSpec_Free(V) { \
  Free((V)->v); (V)->nv=(V)->cut_against_all_images=(V)->min_bond_ratio=0; }

/* Cut against "ss[]" and all images */
int Config_Voronoi_CUT (V3 s0, int nss, double *ss, Alib_Declare_Config);

/* Cut only against the nearest image of "ss[]" */
int Config_Voronoi_cut (V3 s0, int nss, double *ss, Alib_Declare_Config);

/* Generate grains based on specification "V" and the input configuration */
int Config_Voronoirize (VoronoiSpec *V, Alib_Declare_Config, FILE *info);


/* VASP.c: */
/*************************************************************/
/* Vienna Ab-initio Simulation Package configuration toolkit */
/*************************************************************/

/* Load atomistic configuration from POTCAR and POSCAR files */
void Config_load_from_VASP
(char *POTCAR_fname, char *POSCAR_fname, FILE *info, Alib_Declare_Config);


/* Gaussian.c: */
/**********************************/
/* Gaussian configuration toolkit */
/**********************************/

/* Load atomistic configuration from FChk file */
void Config_load_from_Gaussian_FChk
(char *FChk_fname, V3 padding, FILE *info, Alib_Declare_Config);
#define Config_LOAD_from_Gaussian_FChk(fname,padding,Config_Alib_to_Alib) \
  Config_load_from_Gaussian_FChk(fname,padding,stdout,Config_Alib_to_Alib)


/* LeastSquareStrain.c: */
/*****************/
/* Atomic Strain */
/*****************/

/* A unique integer that depends on the chemical symbol sequence */
unsigned int ConfigChecksum (Alib_Declare_Config);

/* IsoAtomic means two configurations with the */
/* same set of atoms and the same indexing.    */
typedef struct
{
    unsigned int checksum;
    double H[DIMENSION][DIMENSION];
    double *s;
} IsoAtomicReference;

/* Save a reference configuration. Remember to free it some time! */
void IsoAtomicReferenceReImprint
(Alib_Declare_Config, IsoAtomicReference *ref);

/* Free the reference configuration */
void IsoAtomicReferenceFree (IsoAtomicReference *ref);

/* Compute and save the least-square transformation */
/* matrix  dxij(new) \approx dxij(ref) J.           */
void ComputeLeastSquareDeformationGradient
(IsoAtomicReference *ref, Alib_Declare_Config, Neighborlist *N, M3 **J);


/* ChainofStates.c: */
/***************************************************************/
/* Toolkit for chain-of-states methods                         */
/*                                                             */
/* Henkelman, Jóhannesson and Jónsson,                         */
/* Methods for Finding Saddle Points and Minimum Energy Paths, */
/* in Progress on Theoretical Chemistry and Physics,           */
/* Ed. S. D. Schwartz (Kluwer Academic Publishers, 2000)       */
/* pp. 269-300.                                                */
/*                                                             */
/* Henkelman and Jónsson,                                      */
/* J Chem. Phys. 113 (2000) 9978-9985; 9901-9904.              */
/*                                                             */
/* The fundamental assumption is that np, symbol, mass, etc.   */
/* intrinsic properties do not change from chain node to chain */
/* node. Only H[][], s[], s1[] may change.                     */
/***************************************************************/

typedef struct
{
    int N;
    unsigned int checksum;
    double *H;
    double *s;
    double *s1;
    double alpha; /* the space is [X,Y]=[s*H0,alpha*H] */
    /* by default, alpha=np^(1/6) to equalize the stiffness */
} ChainofStates;

#define CHAIN_OF_STATES_NO_VELOCITY 1

/* Save current configuration to ChainofStates node i */
void ChainofStates_Push (Alib_Declare_Config, ChainofStates *c, int i);

/* Retrieve ChainofStates node i to current configuration */
void ChainofStates_Pop (ChainofStates *c, int i, Alib_Declare_Config);

/* Reload ChainofStates c based on a sorted UNIX glob file pattern */
void ChainofStates_Load (char *GlobPattern, FILE *info, int flags,
                         Alib_Declare_Config, ChainofStates *c);

/* Do not save velocity in the ChainofStates c */
#define ChainofStates_LoaD(GlobPattern,info,Config_Alib_to_Alib,c) \
  ChainofStates_Load(GlobPattern,info,CHAIN_OF_STATES_NO_VELOCITY, \
  Config_Alib_to_Alib,c)
#define ChainofStates_LOAD(GlobPattern,Config_Alib_to_Alib,c) \
  ChainofStates_Load(GlobPattern,stdout,CHAIN_OF_STATES_NO_VELOCITY,\
  Config_Alib_to_Alib,c)

/* NULL-safe freeing of memory H, s, s1 */
void ChainofStates_Free (ChainofStates *c);

/* Save chain of states c to a directory ("Cfg/") or pattern ("Cfg/#") */
void ChainofStates_Save
(Alib_Declare_Config, ChainofStates *c, FILE *info, char *pattern);

#define ChainofStates_SAVE(Config_Alib_to_Alib,c,pattern) \
  ChainofStates_Save(Config_Alib_to_Alib,c,stdout,pattern)


typedef struct
{
    double dH[9];
    double *ds;
    double *ds1;
    double dR2; /* according to [x, y] = [s*H0, alpha*H] measure */
    double dR;
} ChainofStates_Difference;

/* Compute the difference between node j and node i in hyperspace */
void ChainofStates_Diff (Alib_Declare_Config, ChainofStates *c, int j, int i,
                         ChainofStates_Difference *cd);

/* Total path length in [X,Y]=[s*H0,alpha*H] space */
double ChainofStates_TotalPathLength (int N, ChainofStates_Difference *cd);

/* NULL-safe freeing of memory */
void ChainofStates_Difference_Free(ChainofStates_Difference *cd);

/* Use linear interpolation in hyperspace to add nodes. plan[i] contains */
/* the number of new nodes to be added between original node i and i+1.  */
void ChainofStates_Add_Nodes
(int *plan, Alib_Declare_Config, ChainofStates *c);

typedef struct
{
    int size;
    Neighborlist *pool;
} ChainofStates_NeighborlistPool;

/* NULL-safe freeing of NeighborlistPool memory */
void ChainofStates_NeighborlistPool_Free(ChainofStates_NeighborlistPool *Np);

/* Create a neighborlist pool. cramped=0: a neighborlist for */
/* everyone;    cramped=1: everyone shares one neighborlist. */
void ChainofStates_NeighborlistPool_Recreate
(int cramped, ChainofStates *c, ChainofStates_NeighborlistPool *Np);

/* Swap out a member from the pool */
Neighborlist *ChainofStates_NeighborlistPool_Select
(ChainofStates_NeighborlistPool *Np, int i);

#endif  /* _Atoms_h */
