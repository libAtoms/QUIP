/***********************************************/
/* Atoms: -llapack -lblas -lm                  */
/*        -lVecMat3 -lVecMat -lScalar -lIO     */
/*                                             */
/* Frequently used Physical Constants, Macros, */
/* and Subroutines for Atomistic Simulation.   */
/*                                             */
/* Dec.12, 1999  Ju Li <liju99@mit.edu>        */
/***********************************************/

#ifndef _NIST_h
#define _NIST_h

#include <Scalar.h>

/* everything is measured against the SI system */

/* physical constants are from CODATA 98: */
/* http://physics.nist.gov/cuu/index.html */

/*******************************************************/
/* M for Meter, S for Second, KG for Kilogram          */
/* R for Rad, HZ for Hertz, N for Newton, J for Joule, */
/* W for Watt, PA for Pascal, _ for *, __ for /        */
/* K for Kelvin, C for Coulomb, T for Tesla            */
/*******************************************************/

#define M_IN_M      (1)
#define CM_IN_M     (1e-2)
#define M_IN_CM      (1./CM_IN_M)
#define CM2_IN_M2   (CM_IN_M * CM_IN_M)
#define M2_IN_CM2    (1./CM2_IN_M2)
/* barn is unit for cross-section */
#define BARN_IN_M2  (1e-24 * CM2_IN_M2)
#define M2_IN_BARN   (1./BARN_IN_M2)
#define CM3_IN_M3   (CM2_IN_M2 * CM_IN_M)
#define M3_IN_CM3    (1./CM3_IN_M3)
#define MM_IN_M     (1e-3)
#define M_IN_MM      (1./MM_IN_M)
#define UM_IN_M     (1e-6)
#define M_IN_UM      (1./UM_IN_M)
#define NM_IN_M     (1e-9)
#define M_IN_NM      (1./NM_IN_M)
#define A_IN_M      (1e-10)
#define M_IN_A       (1./A_IN_M)
#define A2_IN_M2    (A_IN_M * A_IN_M)
#define M2_IN_A2     (1./A2_IN_M2)
#define A3_IN_M3    (A2_IN_M2 * A_IN_M)
#define M3_IN_A3     (1./A3_IN_M3)
#define A_IN_CM     (A_IN_M * M_IN_CM)
#define CM_IN_A      (1./A_IN_CM)
#define A2_IN_CM2   (A_IN_CM * A_IN_CM)
#define CM2_IN_A2    (1./A2_IN_CM2)
#define A3_IN_CM3   (A2_IN_CM2 * A_IN_CM)
#define CM3_IN_A3    (1./A3_IN_CM3)
#define A_IN_MM     (A_IN_M * M_IN_MM)
#define MM_IN_A      (1./A_IN_MM)
#define A2_IN_MM2   (A_IN_MM * A_IN_MM)
#define MM2_IN_A2    (1./A2_IN_MM2)
#define A3_IN_MM3   (A2_IN_MM2 * A_IN_MM)
#define MM3_IN_A3    (1./A3_IN_MM3)
#define A_IN_UM     (A_IN_M * M_IN_UM)
#define UM_IN_A      (1./A_IN_UM)
#define A2_IN_UM2   (A_IN_UM * A_IN_UM)
#define UM2_IN_A2    (1./A2_IN_UM2)
#define A3_IN_UM3   (A2_IN_UM2 * A_IN_UM)
#define UM3_IN_A3    (1./A3_IN_UM3)
#define A_IN_NM     (A_IN_M * M_IN_NM)
#define NM_IN_A      (1./A_IN_NM)
#define A2_IN_NM2   (A_IN_NM * A_IN_NM)
#define NM2_IN_A2    (1./A2_IN_NM2)
#define A3_IN_NM3   (A2_IN_NM2 * A_IN_NM)
#define NM3_IN_A3    (1./A3_IN_NM3)

#define NM_IN_UM    (NM_IN_M * M_IN_UM)
#define UM_IN_NM     (1./NM_IN_UM)

#define FEET_IN_M   (3.048e-01)
#define M_IN_FEET    (1./FEET_IN_M)
#define INCH_IN_M   (FEET_IN_M / 12)
#define M_IN_INCH    (1./INCH_IN_M)

#define MILE_IN_M   (FEET_IN_M * 5280.)
#define M_IN_MILE   (1./MILE_IN_M)

#define YARD_IN_M   (FEET_IN_M * 3.)
#define M_IN_YARD   (1./YARD_IN_M)

#define S_IN_S     (1)
#define MS_IN_S    (1e-3)
#define S_IN_MS     (1./MS_IN_S)
#define US_IN_S    (1e-6)
#define S_IN_US     (1./US_IN_S)
#define NS_IN_S    (1e-9)
#define S_IN_NS     (1./NS_IN_S)
#define PS_IN_S    (1e-12)
#define S_IN_PS     (1./PS_IN_S)
#define FS_IN_S    (1e-15)
#define S_IN_FS     (1./FS_IN_S)
#define MIN_IN_S   (60.)
#define S_IN_MIN    (1./MIN_IN_S)
#define HOUR_IN_S  (24 * MIN_IN_S)
#define S_IN_HOUR   (1./HOUR_IN_S)
/* this is a solar year */
#define YEAR_IN_S  (31556925.974592)
#define S_IN_YEAR   (1./YEAR_IN_S)

#define KG_IN_KG   (1)
#define G_IN_KG    (1e-3)
#define KG_IN_G     (1./G_IN_KG)
#define TON_IN_KG  (1e3)
#define KG_IN_TON   (1./TON_IN_KG)

#define AVOIRDUPOIS_POUND_IN_KG (0.4536)
#define AVOIRDUPOIS_OUNCE_IN_KG (AVOIRDUPOIS_POUND_IN_KG / 16)
#define APOTHECARY_POUND_IN_KG  (0.37)
#define APOTHECARY_OUNCE_IN_KG  (APOTHECARY_POUND_IN_KG / 12)

#define N_IN_KG_M__S2   (1.)  /* Newton as force */
#define J_IN_KG_M2__S2  (1.)  /* Joule as energy */
#define W_IN_KG_M2__S3  (1.)  /* Watt as power   */
#define PA_IN_KG__M_S2  (1.)  /* Pascal as pressure */

#define MN_IN_N     (1e-3)
#define N_IN_MN      (1./MN_IN_N)
#define UN_IN_N     (1e-6)
#define N_IN_UN      (1./UN_IN_N)
#define NN_IN_N     (1e-9)
#define N_IN_NN      (1./NN_IN_N)
#define PN_IN_N     (1e-12)
#define N_IN_PN      (1./PN_IN_N)

/* dyne as force */
#define DYNE_IN_N  (0.00001)
#define N_IN_DYNE  (1./DYNE_IN_N)

#define ERG_IN_J   (1E-07)
#define J_IN_ERG   (1./ERG_IN_J)
#define ERG_IN_EV  (ERG_IN_J * J_IN_EV)
#define EV_IN_ERG  (1./ERG_IN_EV)

/* Rad versus Degree for angles */
#define R_IN_D     (180./PI)
#define D_IN_R     (PI/180.)

/* HZ is the unit of frequency (nu), one period per second */
#define HZ_IN__S   (1.)
#define KHZ_IN_HZ  (1e3)
#define HZ_IN_KHZ   (1./KHZ_IN_HZ)
#define MHZ_IN_HZ  (1e6)
#define HZ_IN_MHZ   (1./MHZ_IN_HZ)
#define GHZ_IN_HZ  (1e9)
#define HZ_IN_GHZ   (1./GHZ_IN_HZ)
#define THZ_IN_HZ  (1e12)
#define HZ_IN_THZ   (1./THZ_IN_HZ)

/* R__S is the unit of circular frequency (omega) */
#define NU_IN_HZ_TO_OMEGA_IN_R__S   (2 * PI)
#define OMEGA_IN_R__S_TO_NU_IN_HZ   (1./NU_IN_HZ_TO_OMEGA_IN_R__S)
#define OMEGA_IN_R__S_TO_NU_IN_MHZ  (OMEGA_IN_R__S_TO_NU_IN_HZ * HZ_IN_MHZ)
#define OMEGA_IN_R__S_TO_NU_IN_GHZ  (OMEGA_IN_R__S_TO_NU_IN_HZ * HZ_IN_GHZ)
#define OMEGA_IN_R__S_TO_NU_IN_THZ  (OMEGA_IN_R__S_TO_NU_IN_HZ * HZ_IN_THZ)
#define NU_IN_MHZ_TO_OMEGA_IN_R__S  (1./OMEGA_IN_R__S_TO_NU_IN_MHZ)
#define NU_IN_GHZ_TO_OMEGA_IN_R__S  (1./OMEGA_IN_R__S_TO_NU_IN_GHZ)
#define NU_IN_THZ_TO_OMEGA_IN_R__S  (1./OMEGA_IN_R__S_TO_NU_IN_THZ)
/* mnemonics */
#define HZ_IN_R__S    NU_IN_HZ_TO_OMEGA_IN_R__S
#define R__S_IN_HZ    OMEGA_IN_R__S_TO_NU_IN_HZ
#define OMEGA_TO_HZ   OMEGA_IN_R__S_TO_NU_IN_HZ
#define OMEGA_TO_MHZ  OMEGA_IN_R__S_TO_NU_IN_MHZ
#define OMEGA_TO_GHZ  OMEGA_IN_R__S_TO_NU_IN_GHZ
#define OMEGA_TO_THZ  OMEGA_IN_R__S_TO_NU_IN_THZ
#define HZ_TO_OMEGA   (1./OMEGA_TO_HZ)
#define MHZ_TO_OMEGA  (1./OMEGA_TO_MHZ)
#define GHZ_TO_OMEGA  (1./OMEGA_TO_GHZ)
#define THZ_TO_OMEGA  (1./OMEGA_TO_THZ)

#define J_IN_J     (1)
#define KJ_IN_J    (1000.)
#define J_IN_KJ    (1./KJ_IN_J)
/* calorie as heat & energy */
#define CAL_IN_J   (4.1868)
#define KCAL_IN_J  (CAL_IN_J * 1000.)

#define MPA_IN_PA   (1e6)
#define PA_IN_MPA    (1./MPA_IN_PA)
#define GPA_IN_PA   (1e9)
#define PA_IN_GPA    (1./GPA_IN_PA)
#define ATM_IN_PA   (1.013250e+05)
#define PA_IN_ATM    (1./ATM_IN_PA)
#define ATM_IN_MPA  (ATM_IN_PA * PA_IN_MPA)
#define MPA_IN_ATM   (1./ATM_IN_MPA)
#define ATM_IN_GPA  (ATM_IN_PA * PA_IN_GPA)
#define GPA_IN_ATM   (1./ATM_IN_GPA)
#define BAR_IN_PA   (1e+05)
#define PA_IN_BAR    (1./BAR_IN_PA)
#define BAR_IN_MPA  (BAR_IN_PA * PA_IN_MPA)
#define MPA_IN_BAR   (1./BAR_IN_MPA)
#define BAR_IN_GPA  (BAR_IN_PA * PA_IN_GPA)
#define GPA_IN_BAR    (1./BAR_IN_GPA)
#define MBAR_IN_BAR  (1e+06)
#define BAR_IN_MBAR   (1./MBAR_IN_BAR)
#define MBAR_IN_PA   (MBAR_IN_BAR * BAR_IN_PA)
#define PA_IN_MBAR    (1./MBAR_IN_PA)
#define MBAR_IN_MPA  (MBAR_IN_BAR * BAR_IN_MPA)
#define MPA_IN_MBAR   (1./MBAR_IN_MPA)
#define MBAR_IN_GPA  (MBAR_IN_BAR * BAR_IN_GPA)
#define GPA_IN_MBAR   (1./MBAR_IN_GPA)

/* viscosity */
#define POISEUILLE_IN_PA_S   (1.)
#define PA_S_IN_POISEUILLE   (1./POISEUILLE_IN_PA_S)
#define POISE_IN_POISEUILLE  (0.1)
#define POISEUILLE_IN_POISE  (1./POISE_IN_POISEUILLE)
#define CENTIPOISE_IN_POISE  (0.01)
#define POISE_IN_CENTIPOISE  (1./CENTIPOISE_IN_POISE)
#define CENTIPOISE_IN_POISEUILLE (CENTIPOISE_IN_POISE * POISE_IN_POISEUILLE)
#define POISEUILLE_IN_CENTIPOISE (1./CENTIPOISE_IN_POISEUILLE)
#define CENTIPOISE_IN_PA_S  (CENTIPOISE_IN_POISEUILLE * POISEUILLE_IN_PA_S)
#define PA_S_IN_CENTIPOISE  (1./CENTIPOISE_IN_PA_S)

/** this world's physical constants **/

/* Boltzmann constant kB */
#ifndef BOLZ_IN_J__K
#define BOLZ_IN_J__K  (1.380658e-23)
#endif
#define BOLZ           BOLZ_IN_J__K

/* Planck constant h */
#ifndef PLANCK_IN_J_S
#define PLANCK_IN_J_S   (6.6260755e-34)
#endif
#define PLANCK_IN_J__HZ  PLANCK_IN_J_S
#define PLANCK           PLANCK_IN_J_S
#define HBAR_IN_J_S__R  (PLANCK_IN_J_S / 2 / PI)
#define HBAR             HBAR_IN_J_S__R

/* Elementary charge e */
#ifndef ELEMENTARY_CHARGE_IN_C
#define ELEMENTARY_CHARGE_IN_C  (1.60217733e-19)
#endif
#define ELEMENTARY_CHARGE        ELEMENTARY_CHARGE_IN_C
#define EV_IN_J  ELEMENTARY_CHARGE_IN_C
#define J_IN_EV    (1./EV_IN_J)
#define MEV_IN_EV  1e+6
#define EV_IN_MEV   (1./MEV_IN_EV)
#define MEV_IN_J   (MEV_IN_EV * EV_IN_J)
#define J_IN_MEV    (1./MEV_IN_J)
#define mEV_IN_EV  1e-3
#define EV_IN_mEV   (1./mEV_IN_EV)
#define mEV_IN_J   (mEV_IN_EV * EV_IN_J)
#define J_IN_mEV    (1./mEV_IN_J)
#define KJ__MOL_IN_EV    (KJ_IN_J * J_IN_EV / AVO)
#define EV_IN_KJ__MOL     (1./KJ__MOL_IN_EV)
#define KCAL__MOL_IN_EV  (KCAL_IN_J * J_IN_EV / AVO)
#define EV_IN_KCAL__MOL   (1./KCAL__MOL_IN_EV)

/* (Electric) epsilon_0 */
#ifndef VACUUM_PERMITTIVITY_IN_C2__J_M  
#define VACUUM_PERMITTIVITY_IN_C2__J_M  (8.854187817e-12)
#endif
#define VACUUM_PERMITTIVITY             VACUUM_PERMITTIVITY_IN_C2__J_M

/* (Magnetic) mu_0 */
#ifndef VACUUM_PERMEABILITY_IN_T2_M3__J  
#define VACUUM_PERMEABILITY_IN_T2_M3__J  (4e-7 * PI)
#endif
#define VACUUM_PERMEABILITY              VACUUM_PERMEABILITY_IN_T2_M3__J

/* speed of light in vacuum */
#define LIGHT_SPEED_IN_M__S \
(1./sqrt(VACUUM_PERMITTIVITY * VACUUM_PERMEABILITY))
#define LIGHT_SPEED  LIGHT_SPEED_IN_M__S

/* wavenumber */
#define WAVENUMBER__CM_IN_HZ  (LIGHT_SPEED * M_IN_CM)
#define WAVENUMBER__CM_IN_J   (WAVENUMBER__CM_IN_HZ * PLANCK_IN_J__HZ)
#define WAVENUMBER__CM_IN_EV  (WAVENUMBER__CM_IN_J * J_IN_EV)
#define WAVENUMBER__CM_IN_mEV (WAVENUMBER__CM_IN_EV * EV_IN_mEV)

/* Avogadro constant */
#ifndef AVO
#define AVO        (6.02214199e+23)
#endif

#define AMU_IN_KG  (G_IN_KG / AVO)
#define KG_IN_AMU  (1./AMU_IN_KG)
#define AMU_IN_G   (1. / AVO)
#define G_IN_AMU   (1./AMU_IN_G)
    
#define M3__MOL_IN_A3__ATOM  (M3_IN_A3 / AVO)
#define A3__ATOM_IN_M3__MOL  (1./M3__MOL_IN_A3__ATOM)
    
#define MOL__M3_IN_ATOM__A3  A3__ATOM_IN_M3__MOL
#define ATOM__A3_IN_MOL__M3  M3__MOL_IN_A3__ATOM

#define G__CM3_IN_KG__M3   (G_IN_KG / CM3_IN_M3)
#define KG__M3_IN_G__CM3   (1./G__CM3_IN_KG__M3)
#define G__CM3_IN_AMU__A3  (G_IN_AMU / CM3_IN_A3)
#define AMU__A3_IN_G__CM3  (1./G__CM3_IN_AMU__A3)

#define M3__KG_IN_CM3__G   G__CM3_IN_KG__M3
#define CM3__G_IN_M3__KG   KG__M3_IN_G__CM3
#define A3__AMU_IN_CM3__G  G__CM3_IN_AMU__A3
#define CM3__G_IN_A3__AMU  AMU__A3_IN_G__CM3

/* electron's rest mass */
#ifndef ELECTRON_REST_MASS_IN_KG
#define ELECTRON_REST_MASS_IN_KG  (9.1093897e-31)
#endif
#define ELECTRON_REST_MASS         ELECTRON_REST_MASS_IN_KG

/* Newton's "universal" law */
#ifndef GRAVITATIONAL_CONSTANT_IN_J_M__KG2
#define GRAVITATIONAL_CONSTANT_IN_J_M__KG2  (6.6873e-11)
#endif
#define GRAVITATIONAL_CONSTANT              GRAVITATIONAL_CONSTANT_IN_J_M__KG2

/* hydrogen's electronic semi-classical orbit */
#ifndef BOHR_RADIUS_IN_M
#define BOHR_RADIUS_IN_M  (5.29177e-11)
#endif
#define BOHR_RADIUS_IN_A  (BOHR_RADIUS_IN_M * M_IN_A)
#define HARTREE_IN_J  (ELEMENTARY_CHARGE * ELEMENTARY_CHARGE \
		       / 4 / PI / VACUUM_PERMITTIVITY / BOHR_RADIUS_IN_M)
#define HARTREE_IN_EV  (HARTREE_IN_J * J_IN_EV)
#define HARTREE_IN_KJ__MOL  (HARTREE_IN_EV * EV_IN_KJ__MOL)
#define KJ__MOL_IN_HARTREE  (1./HARTREE_IN_KJ__MOL)
#define HARTREE_IN_KCAL__MOL  (HARTREE_IN_EV * EV_IN_KCAL__MOL)
#define KCAL__MOL_IN_HARTREE  (1./HARTREE_IN_KCAL__MOL)

#define RYDBERG_IN_EV  (HARTREE_IN_EV / 2.)
#define EV_IN_RYDBERG  (1./RYDBERG_IN_EV)

/* earth's average gravitational acceleration */
#ifndef EARTH_G_IN_M__S2
#define EARTH_G_IN_M__S2  (9.80665)
#endif
#define EARTH_G           EARTH_G_IN_M__S2

/* temperature measures */
#define ABSOLUTE_ZERO_SHIFT  273.15
#define CENTIGRADE_TO_K(C)  ((C)+ABSOLUTE_ZERO_SHIFT)
#define K_TO_CENTIGRADE(K)  ((K)-ABSOLUTE_ZERO_SHIFT)
#define FAHRENHEIT_TO_CENTIGRADE(F)  (((F)-32.)*5./9.)
#define CENTIGRADE_TO_FAHRENHEIT(C)  ((C)*9./5.+32.)
#define FAHRENHEIT_TO_K(F) (CENTIGRADE_TO_K(FAHRENHEIT_TO_CENTIGRADE(F)))
#define K_TO_FAHRENHEIT(K) (CENTIGRADE_TO_FAHRENHEIT(K_TO_CENTIGRADE(K)))

#endif  /* _NIST_h */
