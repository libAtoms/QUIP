#include "stdio.h"

extern void quip_wrapper_simple_(int*,double*,int*,double*,double*,double*,double*);

main()
{
  int n=2;
  double lattice[3][3];
  int Z[2];
  double coord[2][3];
  double energy;
  double force[2][3];
  double virial[3][3];

  lattice[0][0]=lattice[0][1]=lattice[0][2]=0.0;
  lattice[1][0]=lattice[1][1]=lattice[1][2]=0.0;
  lattice[2][0]=lattice[2][1]=lattice[2][2]=0.0;
  lattice[0][0]=lattice[1][1]=lattice[2][2]=20.0;
  Z[0] = Z[1] = 29;
  coord[0][0] = -7.110371; coord[0][1] = -3.533572; coord[0][2] =  2.147261;
  coord[1][0] = -7.933029; coord[1][1] = -3.234956; coord[1][2] =  2.573383;

  quip_wrapper_simple_(&n, (double*)lattice,Z,(double*)coord,&energy,(double*)force,(double*)virial);

  printf("Energy = %e\n",energy);

}
