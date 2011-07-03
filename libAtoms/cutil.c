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
// cutil.c. : C Utility Functions to do things that fortran 95 can't do
// 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef DARWIN
#include <sys/sysinfo.h>
#else
#include <sys/sysctl.h>
#endif


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
//X   subroutine CONSTRAINT(pos, velo, t, data, C, dC_dr, dC_dt, target_v)	     
//X     real(dp), dimension(:),         intent(in)  :: pos, velo, data 
//X     real(dp),                       intent(in)  :: t
//X     real(dp),                       intent(out) :: C		     
//X     real(dp), dimension(size(pos)), intent(out) :: dC_dr	     
//X     real(dp),                       intent(out) :: dC_dt           
//X     real(dp),                       intent(out) :: target_v
//X     ...
//X   end subroutine CONSTRAINT 
//X
//X  Remember fortran passes everything by reference, hence all the asterisks
//X  
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


typedef struct{
  void (*sub)(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*);
} VCONSTRAINTSUB_TABLE;

VCONSTRAINTSUB_TABLE constraintsub_table[20];
static int nconstraintsub = 0;

void register_constraint_sub_(void (*sub)(double*,double*,double*,double*,double*,double*,double*,double*,double*,double*,double*)){
  constraintsub_table[nconstraintsub++].sub = sub;
}

void call_constraint_sub_(int* i, double* pos, double* velo, double* mass, double* t, 
			  double* data, double* C, double* dC_dr, double* dC_dt, double* dcoll_dr, double* Z_coll, double* target_v){
  constraintsub_table[*i].sub(pos, velo, mass, t, data, C, dC_dr, dC_dt, dcoll_dr, Z_coll, target_v);
}

// some systems might not support isnan() from Fortran so wrap it here


int fisnan_(double *r)
{
  return isnan(*r);
}

int ffsize_(char *filename, int len) {
   struct stat st;
   char c_file[1024];

   strncpy(c_file, filename, (len < 1023 ? len : 1023));
   c_file[(len < 1023 ? len : 1023)] = 0;
   if (!stat(c_file, &st)) {
      return ((int)st.st_size);
   } else {
      return (-1);
   }
}


// Abort and give a stack trace from fortran

void fabort_() {
  abort();
}

void fusleep_(int *amt) {
   usleep((useconds_t)(*amt));
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

void mem_info_(double *total_mem, double *free_mem)
{
#ifndef DARWIN
   struct sysinfo s_info;
#else
   long int totalram;
   size_t size = sizeof(totalram);
#endif
   int error;

#ifndef DARWIN
   error = sysinfo(&s_info);
   *total_mem = s_info.totalram*s_info.mem_unit;
   *free_mem = s_info.freeram*s_info.mem_unit;
#else
   error = sysctlbyname("hw.memsize", &totalram, &size, NULL, 0);
   *total_mem = totalram;
   *free_mem = *total_mem; // What is the MacOS equivalent??
#endif
}

// Callback Potentials

typedef struct{
  void (*sub)(int*);
} VCALLBACKPOTSUB_TABLE;

VCALLBACKPOTSUB_TABLE callbackpotsub_table[20];
static int ncallbackpotsub = 0;

void register_callbackpot_sub_(void (*sub)(int*)){
  callbackpotsub_table[ncallbackpotsub++].sub = sub;
}

void call_callbackpot_sub_(int* i, int* at) {
  callbackpotsub_table[*i].sub(at);
}

// Temp data for error handling routines

char error_h_info[1000];
int error_h_line;
int error_h_kind;

// quippy error_abort() handler

#include <setjmp.h>
jmp_buf environment_buffer;
char abort_message[1024];

void quippy_error_abort_(char *message, int len)
{
  strncpy(abort_message, message, len);
  abort_message[len] = '\0';
  longjmp(environment_buffer,0);
}

void quippy_error_abort_int_handler(int signum)
{
  char message[] = "Interrupt occured";
  quippy_error_abort_(message, strlen(message));
}

