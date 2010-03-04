//!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//!X
//!X     libAtoms: atomistic simulation library
//!X     
//!X     Copyright 2006-2007.
//!X
//!X     Authors: Gabor Csanyi, Steven Winfield, James Kermode
//!X     Contributors: Noam Bernstein, Alessio Comisso
//!X
//!X     The source code is released under the GNU General Public License,
//!X     version 2, http://www.gnu.org/copyleft/gpl.html
//!X
//!X     If you would like to license the source code under different terms,
//!X     please contact Gabor Csanyi, gabor@csanyi.net
//!X
//!X     When using this software, please cite the following reference:
//!X
//!X     http://www.libatoms.org
//!X
//!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//
//
// cutil.c. : C Utility Functions to do things that fortran 95 can't do
// 
//
//
// $Id: cutil.c,v 1.6 2008-02-11 18:47:38 nb326 Exp $
//
// $Log: not supported by cvs2svn $
// Revision 1.5  2007/10/02 21:04:42  nb326
// add external pointer_to
//
// Revision 1.4  2007/09/17 15:36:00  nb326
// Fix fisnan to receive pointer, and add c_increase_stack_()
//
// Revision 1.3  2007/08/30 14:32:59  gc121
//  now include math.h
//
// Revision 1.2  2007/08/10 16:53:14  gc121
// added fisnan
//
// Revision 1.1  2007/04/18 17:37:58  jrk33
// Group C utility functions in one file
//
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/resource.h>

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//X
//X  Constraint Pointers: 
//X
//X  Contains a struct used for accessing constraint subroutines 
//X  via pointers
//X
//X  The definition of the constraint subroutine in the struct   
//X  matches the fortran version:					     
//X  								     
//X   subroutine CONSTRAINT(pos, velo, t, data, C, dC_dr, dC_dt)	     
//X     real(dp), dimension(:),         intent(in)  :: pos, velo, data 
//X     real(dp),                       intent(in)  :: t
//X     real(dp),                       intent(out) :: C		     
//X     real(dp), dimension(size(pos)), intent(out) :: dC_dr	     
//X     real(dp),                       intent(out) :: dC_dt           
//X     ...
//X   end subroutine CONSTRAINT 
//X
//X  Remember fortran passes everything by reference, hence all the asterisks
//X  
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


typedef struct{
  void (*sub)(double*,double*,double*,double*,double*,double*,double*);
} VCONSTRAINTSUB_TABLE;

VCONSTRAINTSUB_TABLE constraintsub_table[20];
static int nconstraintsub = 0;

void register_constraint_sub_(void (*sub)(double*,double*,double*,double*,double*,double*,double*)){
  constraintsub_table[nconstraintsub++].sub = sub;
}

void call_constraint_sub_(int* i, double* pos, double* velo, double* t, 
			  double* data, double* C, double* dC_dr, double* dC_dt){
  constraintsub_table[*i].sub(pos, velo, t, data, C, dC_dr, dC_dt);
}

// some systems might not support isnan() from Fortran so wrap it here


int fisnan_(double *r)
{
  return isnan(*r);
}


// Abort and give a stack trace from fortran

void fabort_() {
  abort();
}

// Call system(3) from fortran

void system_command_(char* command, int* status, int len)
{
  char c_command[1025];
  int ret;

  strncpy(c_command, command, (len < 1024) ? len : 1024);
  c_command[ len < 1024 ? len : 1024] = 0;
  ret = system(c_command);
  fflush(stdout);
  if (status) *status = ret;
}

// increase stack from fortran
int c_increase_stack_(int *stack_size) {
  int stat;
  struct rlimit l;

  stat = 0;
  getrlimit(RLIMIT_STACK, &l);
  if (l.rlim_cur < *stack_size) { // Need to increase stack
    if (l.rlim_max >= *stack_size) { // New stack is below maximum
      l.rlim_cur = *stack_size;
      stat = setrlimit(RLIMIT_STACK, &l);
    } else { // New stack was more than max limit
      stat = l.rlim_max;
    }
  }

  return stat;
}

int pointer_to_(void *p) {
  return ((int) p);
}

typedef struct{
  void (*sub)(int*);
} VCALLBACKPOTSUB_TABLE;

VCALLBACKPOTSUB_TABLE callbackpotsub_table[20];
static int ncallbackpotsub = 0;

void register_callbackpot_sub_(void (*sub)(int*, int*, int*, int*, int*)){
  callbackpotsub_table[ncallbackpotsub++].sub = sub;
}

void call_callbackpot_sub_(int* i, int* at, int* calc_energy, int* calc_local_e, int* calc_force, int* calc_virial){
  callbackpotsub_table[*i].sub(at);
}
