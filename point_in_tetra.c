
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "parse.h"

#define _MAX_REAL DBL_MAX // 1.79769e+308
#define _MIN_REAL DBL_MIN // 2.22507e-308

#define _MIN(a,b) (((a)<(b))?(a):(b))
#define _MAX(a,b) (((a)>(b))?(a):(b))

#define PLOCATION_EPS 1.0e-8

typedef struct {
  double x,y,z;
} Point3d;

typedef struct _p_DM *DM;

void solve3x3(double A[3][3],double b[],double x[])
{
  double B[3][3];
  double t4, t6, t8, t10, t12, t14, t17;

  t4 = A[2][0] * A[0][1];
  t6 = A[2][0] * A[0][2];
  t8 = A[1][0] * A[0][1];
  t10 = A[1][0] * A[0][2];
  t12 = A[0][0] * A[1][1];
  t14 = A[0][0] * A[1][2];
  t17 = 0.1e1 / (t4 * A[1][2] - t6 * A[1][1] - t8 * A[2][2] + t10 * A[2][1] + t12 * A[2][2] - t14 * A[2][1]);

  B[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) * t17;
  B[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) * t17;
  B[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * t17;
  B[1][0] = -(-A[2][0] * A[1][2] + A[1][0] * A[2][2]) * t17;
  B[1][1] = (-t6 + A[0][0] * A[2][2]) * t17;
  B[1][2] = -(-t10 + t14) * t17;
  B[2][0] = (-A[2][0] * A[1][1] + A[1][0] * A[2][1]) * t17;
  B[2][1] = -(-t4 + A[0][0] * A[2][1]) * t17;
  B[2][2] = (-t8 + t12) * t17;

  x[0] = B[0][0]*b[0] + B[0][1]*b[1] + B[0][2]*b[2];
  x[1] = B[1][0]*b[0] + B[1][1]*b[1] + B[1][2]*b[2];
  x[2] = B[2][0]*b[0] + B[2][1]*b[1] + B[2][2]*b[2];
}

/* transform to reference cell [-1,1]^3 */
/*
N0 = 1 - 0.5(1+xi) - 0.5(1+eta) - 0.5(1+zeta)
N1 = 0.5(1+xi)
N2 = 0.5(1+eta)
N3 = 0.5(1+zeta)

x = N0(xi,eta,zeta).x0 + N1(xi,eta,zeta).x1 + N2(xi,eta,zeta).x2 + N3(xi,eta,zeta).x3
y
z

x  =   x0 - 0.5x0 - 0.5x0 - 0.5x0 - 0.5xi.x0 - 0.5eta.x0 - 0.5zeta.x0
0.5x1 + 0.5xi.x1
0.5x2 + 0.5eta.x2
0.5x3 + 0.5zeta.x3

x + 0.5x0 - 0.5x1 - 0.5x2 - 0.5x3 = xi(-0.5x0 + 0.5x1) + eta(-0.5x0 + 0.5x2) + zeta(-0.5x0 + 0.5x3)
*/
void DMTetComputeLocalCoordinateAffine3d(const double X[],const double coords[],double xip[])
{
  double A[3][3],b[3],v1[3],v2[3],v3[3],v4[3];
  double xp,yp,zp;
  int ii,d;

  xp = X[0];
  yp = X[1];
  zp = X[2];

  for (d=0; d<3; d++) {
    ii = 0;  v1[d] = coords[3*ii + d];
    ii = 1;  v2[d] = coords[3*ii + d];
    ii = 2;  v3[d] = coords[3*ii + d];
    ii = 3;  v4[d] = coords[3*ii + d];
  }

  b[0] = xp + 0.5*v1[0] - 0.5*v2[0] - 0.5*v3[0] - 0.5*v4[0];
  b[1] = yp + 0.5*v1[1] - 0.5*v2[1] - 0.5*v3[1] - 0.5*v4[1];
  b[2] = zp + 0.5*v1[2] - 0.5*v2[2] - 0.5*v3[2] - 0.5*v4[2];

  A[0][0] = -0.5*v1[0] + 0.5*v2[0];   A[0][1] = -0.5*v1[0] + 0.5*v3[0];   A[0][2] = -0.5*v1[0] + 0.5*v4[0];
  A[1][0] = -0.5*v1[1] + 0.5*v2[1];   A[1][1] = -0.5*v1[1] + 0.5*v3[1];   A[1][2] = -0.5*v1[1] + 0.5*v4[1];
  A[2][0] = -0.5*v1[2] + 0.5*v2[2];   A[2][1] = -0.5*v1[2] + 0.5*v3[2];   A[2][2] = -0.5*v1[2] + 0.5*v4[2];

  solve3x3(A,b,xip);
}

void PointLocation_BruteForce(
  Mesh dm,
  long int npoints,const double xp[],long int econtaining[],double xip[],long int *found)
{
  const int nsd = 3;
  int e,nelements_g=0,*element_g=NULL,npe_g=0;
  double *coords_g=NULL,elcoords_g[nsd*4],elcentroid[nsd],elcoords_normalized[nsd*4];
  long int p,i;
  bool point_located;
  Point3d point;
  double max_el_L;
  double tolerance;
  double point_normalized[nsd],xi[nsd];

  long int points_located = 0;


  tolerance = 1.0e-10;

  *found = 0;
  nelements_g = dm->ncell;
  npe_g = dm->points_per_cell;
  element_g = dm->cell;
  coords_g = dm->vert;

  for (p=0; p<npoints; p++) {

    point.x = xp[nsd*p+0];
    point.y = xp[nsd*p+1];
    point.z = xp[nsd*p+2];

    /* set defaults */
    econtaining[p] = -1;
    xip[nsd*p  ] = NAN;
    xip[nsd*p+1] = NAN;
    xip[nsd*p+2] = NAN;

    point_located = false;
    for (e=0; e<nelements_g; e++) {
      int *elnidx_g = NULL;
      double bbemin[nsd],bbemax[nsd],basis[4];

      /* get element -> node map for the triangle geometry */
      elnidx_g = &element_g[npe_g*e];

      bbemin[0] = bbemin[1] = bbemin[2] = _MAX_REAL;
      bbemax[0] = bbemax[1] = bbemax[2] = -_MAX_REAL;

      /* get element coordinates, compute the element bounding box and centroid */
      elcentroid[0] = 0.0;
      elcentroid[1] = 0.0;
      elcentroid[2] = 0.0;

      for (i=0; i<npe_g; i++) {
        int nidx = elnidx_g[i];

        elcoords_g[nsd*i+0] = coords_g[nsd*nidx  ];
        elcoords_g[nsd*i+1] = coords_g[nsd*nidx+1];
        elcoords_g[nsd*i+2] = coords_g[nsd*nidx+2];

        elcentroid[0] += 0.25 * elcoords_g[nsd*i+0];
        elcentroid[1] += 0.25 * elcoords_g[nsd*i+1];
        elcentroid[2] += 0.25 * elcoords_g[nsd*i+2];

        bbemin[0] = _MIN(bbemin[0],elcoords_g[nsd*i+0]);
        bbemin[1] = _MIN(bbemin[1],elcoords_g[nsd*i+1]);
        bbemin[2] = _MIN(bbemin[2],elcoords_g[nsd*i+2]);

        bbemax[0] = _MAX(bbemax[0],elcoords_g[nsd*i+0]);
        bbemax[1] = _MAX(bbemax[1],elcoords_g[nsd*i+1]);
        bbemax[2] = _MAX(bbemax[2],elcoords_g[nsd*i+2]);
      }
      /* if point outside the element bounding box - check next element */
      if (point.x < bbemin[0]) { continue; }
      if (point.x > bbemax[0]) { continue; }

      if (point.y < bbemin[1]) { continue; }
      if (point.y > bbemax[1]) { continue; }

      if (point.z < bbemin[2]) { continue; }
      if (point.z > bbemax[2]) { continue; }

      max_el_L = 0.0;
      max_el_L = _MAX(max_el_L,fabs(bbemax[0] - bbemin[0]));
      max_el_L = _MAX(max_el_L,fabs(bbemax[1] - bbemin[1]));
      max_el_L = _MAX(max_el_L,fabs(bbemax[2] - bbemin[2]));

      /* shift and normalize element coordinates */
      for (i=0; i<npe_g; i++) {
        elcoords_normalized[nsd*i+0] = 2.0*(elcoords_g[nsd*i+0] - elcentroid[0])/max_el_L;
        elcoords_normalized[nsd*i+1] = 2.0*(elcoords_g[nsd*i+1] - elcentroid[1])/max_el_L;
        elcoords_normalized[nsd*i+2] = 2.0*(elcoords_g[nsd*i+2] - elcentroid[2])/max_el_L;
      }

      /* shift and normalize point */
      point_normalized[0] = 2.0*(point.x - elcentroid[0])/max_el_L;
      point_normalized[1] = 2.0*(point.y - elcentroid[1])/max_el_L;
      point_normalized[2] = 2.0*(point.z - elcentroid[2])/max_el_L;

      /* perform an isoparametric inversion for xi */
      DMTetComputeLocalCoordinateAffine3d(point_normalized,elcoords_normalized,xi);

      /* check xi, eta, zeta are valid */
      basis[0] = 1.0 - 0.5*(1.0+xi[0]) - 0.5*(1.0+xi[1]) - 0.5*(1.0+xi[2]);
      basis[1] = 0.5*(1+xi[0]);
      basis[2] = 0.5*(1+xi[1]);
      basis[3] = 0.5*(1+xi[2]);

      if ((basis[0] < (0.0-PLOCATION_EPS)) || (basis[0] > (1.0+PLOCATION_EPS))) { continue; }
      if ((basis[1] < (0.0-PLOCATION_EPS)) || (basis[1] > (1.0+PLOCATION_EPS))) { continue; }
      if ((basis[2] < (0.0-PLOCATION_EPS)) || (basis[2] > (1.0+PLOCATION_EPS))) { continue; }
      if ((basis[3] < (0.0-PLOCATION_EPS)) || (basis[3] > (1.0+PLOCATION_EPS))) { continue; }

      //printf("  basis %+1.8e %+1.8e %+1.8e %+1.8e\n",basis[0],basis[1],basis[2],basis[3]);

      /* all the above passed - do an interpolation tests and check the interpolated coordinate is at least 0.01% of the original */
      {
        int d;
        double point_normalized_i[] = { 0.0, 0.0, 0.0 };
        bool diff_failed = false;
        double diff[3];

        for (d=0; d<nsd; d++) {
          for (i=0; i<4; i++) {
            point_normalized_i[d] += basis[i] * elcoords_normalized[nsd*i+d];
          }
        }
        //printf("  cell %d\n\t: exact  %+1.10e %+1.10e %+1.10e\n\t: interp %+1.10e %+1.10e %+1.10e\n",e,point_normalized[0],point_normalized[1],point_normalized[2],point_normalized_i[0],point_normalized_i[1],point_normalized_i[2]);
        for (d=0; d<nsd; d++) {

          diff[d] = fabs((point_normalized[d] - point_normalized_i[d])/1.0);
          if (diff[d] > tolerance) { diff_failed = true; }
        }
        if (diff_failed) {
          continue;
        }
      }

      point_located = true;
      points_located++;
      break;
    }

    if (point_located) {
      /* copy element index */
      econtaining[p] = e;

      /* copy xi */
      xip[nsd*p+0] = xi[0];
      xip[nsd*p+1] = xi[1];
      xip[nsd*p+2] = xi[2];
    }

  }

  /* rescale xi if required */
  for (p=0; p<npoints; p++) {

    /* only attempt to re-scale if the point was actually found  - otherwise xip is NAN */
    if (econtaining[p] == -1) { continue; }

    xip[nsd*p+0] = 0.5*(xip[nsd*p+0] + 1.0);
    xip[nsd*p+1] = 0.5*(xip[nsd*p+1] + 1.0);
    xip[nsd*p+2] = 0.5*(xip[nsd*p+2] + 1.0);
  }
  *found = points_located;

  printf("PointLocation_BruteForce:\n");
  printf("  npoints queried    : %ld\n",npoints);
  printf("  npoints located    : %ld\n",points_located);
  printf("  npoints not located: %ld\n",npoints-points_located);
}

void _initialize_points(long int npoints,long int econtaining[],double xip[])
{
  const int nsd = 3;
  long int p;

  for (p=0; p<npoints; p++) {
    /* set defaults */
    econtaining[p] = -1;
    xip[nsd*p  ] = NAN;
    xip[nsd*p+1] = NAN;
    xip[nsd*p+2] = NAN;
  }
}

void _finalize_points(long int npoints,long int econtaining[],double xip[])
{
  const int nsd = 3;
  long int p;

  for (p=0; p<npoints; p++) {
    /* only attempt to re-scale if the point was actually found  - otherwise xip is NAN */
    if (econtaining[p] == -1) { continue; }

    xip[nsd*p+0] = 0.5*(xip[nsd*p+0] + 1.0);
    xip[nsd*p+1] = 0.5*(xip[nsd*p+1] + 1.0);
    xip[nsd*p+2] = 0.5*(xip[nsd*p+2] + 1.0);
  }
}

bool _point_in_bounding_box(const double point[],const double cmin[],const double cmax[])
{
  if (point[0] < cmin[0]) { return(false); }
  if (point[0] > cmax[0]) { return(false); }

  if (point[1] < cmin[1]) { return(false); }
  if (point[1] > cmax[1]) { return(false); }

  if (point[2] < cmin[2]) { return(false); }
  if (point[2] > cmax[2]) { return(false); }
  return(true);
}

bool point_in_tet(const double point[],double xi[],int npe_g,const double elcoords_g[],double tolerance)
{
  const int nsd = 3;
  double bbemin[nsd],bbemax[nsd],elcentroid[nsd],elcoords_normalized[nsd*4];
  double point_normalized[nsd],point_normalized_i[nsd];
  double basis[4],max_el_L,diff;
  int i,d;
  bool point_in_bb = false;


  for (d=0; d<nsd; d++) {
    bbemin[d] = _MAX_REAL;
    bbemax[d] = -_MAX_REAL;
    elcentroid[d] = 0.0;
    point_normalized_i[d] = 0.0;
  }

  for (i=0; i<npe_g; i++) {
    for (d=0; d<nsd; d++) {
      elcentroid[d] += 0.25 * elcoords_g[nsd*i+d];
      bbemin[d] = _MIN(bbemin[d],elcoords_g[nsd*i+d]);
      bbemax[d] = _MAX(bbemax[d],elcoords_g[nsd*i+d]);
    }
  }

  point_in_bb = _point_in_bounding_box(point,(const double*)bbemin,(const double*)bbemax);
  if (point_in_bb == false) return(false);

  /* if point outside the element bounding box - check next element */
  max_el_L = 0.0;
  for (d=0; d<nsd; d++) {
    max_el_L = _MAX(max_el_L,fabs(bbemax[d] - bbemin[d]));
  }

  /* shift and normalize element coordinates */
  for (i=0; i<npe_g; i++) {
    for (d=0; d<nsd; d++) {
      elcoords_normalized[nsd*i+d] = 2.0*(elcoords_g[nsd*i+d] - elcentroid[d])/max_el_L;
    }
  }

  /* shift and normalize point */
  for (d=0; d<nsd; d++) {
    point_normalized[d] = 2.0*(point[d] - elcentroid[d])/max_el_L;
  }
  /* perform an isoparametric inversion for xi */
  DMTetComputeLocalCoordinateAffine3d(point_normalized,elcoords_normalized,xi);

  /* check xi, eta, zeta are valid */
  basis[0] = 1.0 - 0.5*(1.0+xi[0]) - 0.5*(1.0+xi[1]) - 0.5*(1.0+xi[2]);
  basis[1] = 0.5*(1+xi[0]);
  basis[2] = 0.5*(1+xi[1]);
  basis[3] = 0.5*(1+xi[2]);

  if ((basis[0] < (0.0-PLOCATION_EPS)) || (basis[0] > (1.0+PLOCATION_EPS))) { return(false); }
  if ((basis[1] < (0.0-PLOCATION_EPS)) || (basis[1] > (1.0+PLOCATION_EPS))) { return(false); }
  if ((basis[2] < (0.0-PLOCATION_EPS)) || (basis[2] > (1.0+PLOCATION_EPS))) { return(false); }
  if ((basis[3] < (0.0-PLOCATION_EPS)) || (basis[3] > (1.0+PLOCATION_EPS))) { return(false); }

  for (d=0; d<nsd; d++) {
    for (i=0; i<4; i++) {
      point_normalized_i[d] += basis[i] * elcoords_normalized[nsd*i+d];
    }
  }

  for (d=0; d<nsd; d++) {
    diff = fabs((point_normalized[d] - point_normalized_i[d])/1.0);
    if (diff > tolerance) { return(false); }
  }
  return(true);
}

void PointLocation_PartitionedBoundingBox(
  Mesh dm,
  long int npoints,const double xp[],long int econtaining[],double xip[],
  long int *_npoints_located)
{
  const int nsd = 3;
  const double tolerance = 1.0e-10;

  int d,e,*element_g=NULL,npe_g=0,pe,nelements_partition=0,part_idx=-1,npartitions;
  double *coords_g=NULL,elcoords_g[nsd*4];
  long int p,i;
  bool point_located;
  double xi[nsd];
  CellPartition part=NULL;
  long int npoints_located = 0;

  *_npoints_located = 0;
  npe_g = dm->points_per_cell;
  element_g = dm->cell;
  coords_g = dm->vert;
  npartitions = dm->npartition;

  _initialize_points(npoints,econtaining,xip);

  for (p=0; p<npoints; p++) {
    const double *_point = &xp[nsd*p];

    point_located = false;

    if (p > 0 && econtaining[p-1] >= 0) {
      int *elnidx_g = NULL;

      e = econtaining[p-1];

      elnidx_g = &element_g[npe_g*e];
      for (i=0; i<npe_g; i++) {
        int nidx = elnidx_g[i];
        for (d=0; d<nsd; d++) {
          elcoords_g[nsd*i+d] = coords_g[nsd*nidx+d];
        }
      }

      point_located = point_in_tet(_point,xi,npe_g,(const double*)elcoords_g,tolerance);
    }

    if (point_located) {
      /* copy element index */
      econtaining[p] = e;

      /* copy xi */
      xip[nsd*p+0] = xi[0];
      xip[nsd*p+1] = xi[1];
      xip[nsd*p+2] = xi[2];
      npoints_located++;
      continue;
    }

    /* if point outside bounding box - check next point */
    for (part_idx=0; part_idx<npartitions; part_idx++) {
      part = dm->partition[part_idx];
      const double *cmin = part->cmin;
      const double *cmax = part->cmax;
      bool point_in_bb = false;

      point_in_bb = _point_in_bounding_box(_point,cmin,cmax);
      if (point_in_bb == false) continue;

      nelements_partition = part->ncell;
      for (pe=0; pe<nelements_partition; pe++) {
        int *elnidx_g = NULL;

        /* get element -> node map for the triangle geometry */
        e = part->cell_list[pe];
        elnidx_g = &element_g[npe_g*e];

        for (i=0; i<npe_g; i++) {
          int nidx = elnidx_g[i];
          for (d=0; d<nsd; d++) {
            elcoords_g[nsd*i+d] = coords_g[nsd*nidx+d];
          }
        }

        point_located = point_in_tet(_point,xi,npe_g,(const double*)elcoords_g,tolerance);
        if (point_located) {
          break; // break from cell-partition sweep
        }
      }

      if (point_located) {
        break; // break from partition sweep
      }

    }

    if (point_located) {
      /* copy element index */
      econtaining[p] = e;

      /* copy xi */
      xip[nsd*p+0] = xi[0];
      xip[nsd*p+1] = xi[1];
      xip[nsd*p+2] = xi[2];

      npoints_located++;
    }

  }

  _finalize_points(npoints,econtaining,xip);
  *_npoints_located = npoints_located;

  printf("PointLocation_PartitionedBoundingBox:\n");
  printf("  npoints queried    : %ld\n",npoints);
  printf("  npoints located    : %ld\n",npoints_located);
  printf("  npoints not located: %ld\n",npoints-npoints_located);
}
