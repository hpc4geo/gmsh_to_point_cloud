
#ifndef __point_in_tetra_h__
#define __point_in_tetra_h__

#include "parse.h"

void PointLocation_BruteForce(
  Mesh dm,
  long int npoints,const double xp[],long int econtaining[],double xip[],
  long int *npoints_located);

void PointLocation_PartitionedBoundingBox(
  Mesh dm,
  long int npoints,const double xp[],long int econtaining[],double xip[],
  long int *npoints_located);

#endif
