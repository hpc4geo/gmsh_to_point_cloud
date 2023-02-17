
/*

gcc -c -O2 -Wall -std=c99 parse.c
gcc -c -O2 -Wall -std=c99  point_in_tetra.c
gcc -O2 -Wall -std=c99 demo_parse.c -o demo parse.o point_in_tetra.o

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "parse.h"
#include "point_in_tetra.h"

void test_1(Mesh mesh, const double xp[],int method)
{
  long int np = 1,found;
  double xip[] = {0,0,0};
  long int ep[] = {-1};

  printf("------------------------------\n");
  switch (method) {
    case 0:
    PointLocation_BruteForce(mesh,np,xp,ep,xip,&found);
    break;
    case 1:
    PointLocation_PartitionedBoundingBox(mesh,np,(const double*)xp,ep,xip,&found);
    break;
  }
  printf("point x:[%+1.2e,%+1.2e,%+1.2e] -> e:%ld -> xi:[%+1.2e,%+1.2e,%+1.2e]\n",xp[0],xp[1],xp[2],ep[0],xip[0],xip[1],xip[2]);
}

void test_2_reg(Mesh mesh,int _nx,int _ny,int _nz,int method)
{
  clock_t start,diff;
  long int p,np,found;
  double *xp,*xip;
  long int *ep;
  long int i,j,k,nxp,nyp,nzp;
  double dx,dy,dz;

  printf("------------------------------\n");
  nxp = _nx;
  nyp = _ny;
  nzp = _nz;
  np = nxp * nyp * nzp;

  xp = (double*)malloc(np * 3 * sizeof(double));
  xip = (double*)malloc(np * 3 * sizeof(double));
  ep = (long int*)malloc(np * sizeof(long int));

  dx = 1.0e6 / (double)nxp;
  dy = 200000.0 / (double)nyp;
  dz = 1.0e6 / (double)nzp;

  p = 0;
  for (k=0; k<nzp; k++) {
    for (j=0; j<nyp; j++) {
      for (i=0; i<nxp; i++) {
        xp[3*p+0] = 0.0 + 0.5 * dx + i*dx;
        xp[3*p+1] = 0.0 + 0.5 * dy + j*dy;
        xp[3*p+2] = 0.0 + 0.5 * dz + k*dz;
        p++;
      }
    }
  }

  start = clock();
  switch (method) {
    case 0:
      PointLocation_BruteForce(mesh,np,(const double*)xp,ep,xip,&found);
    break;
    case 1:
      PointLocation_PartitionedBoundingBox(mesh,np,(const double*)xp,ep,xip,&found);
    break;
  }
  diff = clock() - start;
  double job_sec = (double)diff / (double)CLOCKS_PER_SEC;

  printf("npoints %ld\n",np);
  printf("cells   %d\n",mesh->ncell);
  printf("time %1.4e (sec)\n",job_sec);
  printf("time/point %1.4e (sec)\n",job_sec/(double)np);

  for (p=0; p<np; p++) {
    if (ep[p] == -1) {
      printf("p[%ld] x:[%+1.18e,%+1.18e,%+1.18e]\n",p,xp[3*p+0],xp[3*p+1],xp[3*p+2]);
    }
  }
  free(xp);
  free(xip);
  free(ep);
}

void test_3_rand(Mesh mesh,long int np,int method)
{
  clock_t start,diff;
  long int p,found;
  double *xp,*xip;
  long int *ep;

  printf("------------------------------\n");
  xp = (double*)malloc(np * 3 * sizeof(double));
  xip = (double*)malloc(np * 3 * sizeof(double));
  ep = (long int*)malloc(np * sizeof(long int));

  srand(0);
  for (p=0; p<np; p++) {
    double _x,_y,_z;
    _x = (double)(rand())/(double)(RAND_MAX + 1.0);
    _y = (double)(rand())/(double)(RAND_MAX + 1.0);
    _z = (double)(rand())/(double)(RAND_MAX + 1.0);
    xp[3*p+0] = 1.0e6 * _x;
    xp[3*p+1] = 200000.0 * _y;
    xp[3*p+2] = 1.0e6 * _z;
  }

  start = clock();
  switch (method) {
    case 0:
      PointLocation_BruteForce(mesh,np,(const double*)xp,ep,xip,&found);
    break;
    case 1:
      PointLocation_PartitionedBoundingBox(mesh,np,(const double*)xp,ep,xip,&found);
    break;
  }
  diff = clock() - start;
  double job_sec = (double)diff / (double)CLOCKS_PER_SEC;

  printf("npoints %ld\n",np);
  printf("cells   %d\n",mesh->ncell);
  printf("time %1.4e (sec)\n",job_sec);
  printf("time/point %1.4e (sec)\n",job_sec/(double)np);

  for (p=0; p<np; p++) {
    if (ep[p] == -1) {
      printf("p[%ld] x:[%+1.18e,%+1.18e,%+1.18e]\n",p,xp[3*p+0],xp[3*p+1],xp[3*p+2]);
    }
  }
  free(xp);
  free(xip);
  free(ep);
}

int main(int nargs,char *args[])
{
  Mesh mesh = NULL;

  parse_mesh("md.bin",&mesh);
  if (!mesh) { printf("mesh = NULL. Aborting.\n"); return(0); }
  {
    long int *data = NULL;

    parse_field(mesh,"region_cell.bin",'c',(void**)&data);
    if (!data) { printf("data = NULL. Aborting.\n"); return(0); }
    printf("cell regions %ld %ld %ld %ld %ld\n",data[0],data[1],data[2],data[3],data[4]);
    free(data);
  }
  {
    double *data = NULL;

    parse_field(mesh,"temperature_vertex.bin",'v',(void**)&data);
    if (!data) { printf("data = NULL. Aborting.\n"); return(0); }
    printf("vertex temperature %+1.3e %+1.3e %+1.3e %+1.3e %+1.3e\n",data[0],data[1],data[2],data[3],data[4]);
    free(data);
  }

  /*
  MeshView(mesh);
  */
  {
    const double xp[] = { +3.119925553910434246e+05,+1.317761797457933426e+05,+8.112700940109789371e+05 };
    test_1(mesh,xp,0);
    test_1(mesh,xp,1);
  }

  {
    const double xp[] = { 0,0,0, 1e6,0,0, 0,200000,0, 1e6,200000,0  ,  0,0,1e6, 1e6,0,1e6, 0,200000,1e6, 1e6,200000,1e6  };
    test_1(mesh,&xp[3*0],1);
    test_1(mesh,&xp[3*1],1);
    test_1(mesh,&xp[3*2],1);
    test_1(mesh,&xp[3*3],1);
    test_1(mesh,&xp[3*4],1);
    test_1(mesh,&xp[3*5],1);
    test_1(mesh,&xp[3*6],1);
    test_1(mesh,&xp[3*7],1);
  }

  test_2_reg(mesh,1000,1000,1,0);
  test_2_reg(mesh,1000,1000,1,1);

  test_3_rand(mesh,1000*1000*1,0);
  test_3_rand(mesh,1000*1000*1,1);

  MeshDestroy(&mesh);

  return(0);
}
