
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "parse.h"

void CellPartitionCreate(CellPartition *_c)
{
  CellPartition c = NULL;

  c = (struct _p_CellPartition*)malloc(sizeof(struct _p_CellPartition));
  memset(c,0,sizeof(struct _p_CellPartition));
  c->ncell = 0;
  c->cell_list = NULL;
  *_c = c;
}

void CellPartitionDestroy(CellPartition *_c)
{
  CellPartition c = NULL;

  if (_c) { c = *_c; }
  free(c->cell_list);
  free(c);
  c = NULL;
  *_c = c;
}

void MeshCreate(Mesh *_m)
{
  Mesh m = NULL;

  m = (struct _p_Mesh*)malloc(sizeof(struct _p_Mesh));
  memset(m,0,sizeof(struct _p_Mesh));
  m->nvert = 0;
  m->ncell = 0;
  m->vert = NULL;
  m->cell = NULL;
  m->npartition = 0;
  m->partition = NULL;
  *_m = m;
}

void MeshDestroy(Mesh *_m)
{
  Mesh m = NULL;
  int p;

  if (_m) { m = *_m; }
  free(m->vert);
  free(m->cell);
  for (p=0; p<m->npartition; p++) {
    CellPartitionDestroy(&m->partition[p]);
  }
  free(m->partition);
  m->partition = NULL;
  free(m);
  m = NULL;
  *_m = m;
}

void MeshView(Mesh m)
{
  printf("[Mesh]\n");
  printf("         nvert: %d\n",m->nvert);
  printf("      coor_dim: %d\n",m->coor_dim);
  printf("        ncells: %d\n",m->ncell);
  printf("verts-per-cell: %d",m->points_per_cell);
  switch (m->points_per_cell) {
    case 2:
    printf(" --> line\n");
    break;
    case 3:
    printf(" --> triangle\n");
    break;
    case 4:
    printf(" --> tetrahedral\n");
    break;
    default:
    printf("unknown cell type\n");
    break;
  }
  printf("    npartition: %d\n",m->npartition);
  {
    int p,d;
    for (p=0; p<m->npartition; p++) {
      printf("part %4d:\n",p);
      for (d=0; d<m->coor_dim; d++) {
        printf("  dir[%d] min %+1.6e max %+1.6e\n",d,m->partition[p]->cmin[d],m->partition[p]->cmax[d]);
      }
    }
  }
}

void parse_mesh(const char filename[],Mesh *m)
{
  FILE   *fp=NULL;
  int    nvert,coor_dim,ncell,pp_cell,nparts,p,nc;
  double *coor=NULL,cmin[3],cmax[3];
  int    *cell=NULL,*cell_list=NULL;
  Mesh   mesh = NULL;
  size_t bytes_read=0;

  *m = NULL;
  fp = fopen(filename, "rb");
  if (!fp) { printf("parse_mesh(): File %s was not found or read\n",filename); return; }

  bytes_read = fread(&nvert,sizeof(int),1,fp);
  bytes_read = fread(&coor_dim,sizeof(int),1,fp);
  //printf("nvert %d coor_dim %d\n",nvert,coor_dim);

  coor = (double*)malloc(sizeof(double)*nvert*coor_dim);
  memset(coor,0,sizeof(double)*nvert*coor_dim);
  bytes_read = fread(coor,sizeof(double),nvert*coor_dim,fp);

  bytes_read = fread(&ncell,sizeof(int),1,fp);
  bytes_read = fread(&pp_cell,sizeof(int),1,fp);
  //printf("ncell %d points-per-cell %d\n",ncell,pp_cell);

  cell = (int*)malloc(sizeof(int)*ncell*pp_cell);
  memset(cell,0,sizeof(int)*ncell*pp_cell);
  bytes_read = fread(cell,sizeof(int),ncell*pp_cell,fp);

  switch (pp_cell) {
    case 2:
    break;
    case 3:
    break;
    case 4:
    break;
    default:
    printf("parse_mesh(): unknown cell type - abort\n");
    free(coor);
    free(cell);
    fclose(fp);
    return;
  }

  MeshCreate(&mesh);
  mesh->nvert = nvert;
  mesh->coor_dim = coor_dim;
  mesh->vert = coor;
  mesh->ncell = ncell;
  mesh->points_per_cell = pp_cell;
  mesh->cell = cell;

  nparts = 0;
  bytes_read = fread(&nparts,sizeof(int),1,fp);

  mesh->npartition = nparts;
  mesh->partition = (CellPartition*)malloc(sizeof(CellPartition)*nparts);
  memset(mesh->partition,0,sizeof(CellPartition)*nparts);

  for (p=0; p<nparts; p++) {
    CellPartition cp;

    bytes_read = fread(cmin,sizeof(double),coor_dim,fp);
    bytes_read = fread(cmax,sizeof(double),coor_dim,fp);
    bytes_read = fread(&nc,sizeof(int),1,fp);

    //printf("part %d\n\txmin %+1.4e %+1.4e %+1.4e\n\txmax %+1.4e %+1.4e %+1.4e\n\tncells %d\n",p,cmin[0],cmin[1],cmin[2],cmax[0],cmax[1],cmax[2],nc);

    cell_list = (int*)malloc(sizeof(int)*nc);
    memset(cell_list,0,sizeof(int)*nc);
    bytes_read = fread(cell_list,sizeof(int),nc,fp);

    CellPartitionCreate(&cp);
    mesh->partition[p] = cp;
    cp->ncell = nc;
    cp->cell_list = cell_list;
    memcpy(cp->cmin,cmin,sizeof(double)*coor_dim);
    memcpy(cp->cmax,cmax,sizeof(double)*coor_dim);
  }
  fclose(fp);
  *m = mesh;
}

void parse_field(Mesh m,const char filename[],char ftypevoid,void **_data)
{
  FILE   *fp=NULL;
  int    len,dtype;
  void   *buffer;
  size_t bytes,bytes_item=0,bytes_read=0;
  bool   valid = true;

  *_data = NULL;
  fp = fopen(filename, "rb");
  if (!fp) { printf("parse_field(): File %s was not found or read\n",filename); return; }

  bytes_read = fread(&len,sizeof(int),1,fp);
  switch (ftypevoid) {
    case 'c':
      if (len != m->ncell) { valid = false; }
      break;
    case 'e':
      if (len != m->ncell) { valid = false; }
      break;
    case 'v':
      if (len != m->nvert) { valid = false; }
      break;
  }
  if (!valid) {
    printf("parse_field(): From %s -> unable to parse field. Value of ftypevoid is inconsistent with data.\n",filename);
    return;
  }

  bytes_read = fread(&dtype,sizeof(int),1,fp);
  switch (dtype) {
    case 10:
    bytes_item = sizeof(short);
    break;
    case 11:
    bytes_item = sizeof(int);
    break;
    case 12:
    bytes_item = sizeof(long int);
    break;

    case 20:
    bytes_item = sizeof(float);
    break;
    case 21:
    bytes_item = sizeof(double);
    break;

    default:
    bytes_item = 0;
    valid = false;
    break;
  }
  if (!valid) {
    printf("parse_field(): Unable to parse field data. Data type unrecognized\n");
    return;
  }

  bytes = len * bytes_item;

  buffer = (void*)malloc(bytes);
  memset(buffer,0,bytes);
  bytes_read = fread(buffer,bytes_item,len,fp);
  *_data = buffer;

  fclose(fp);

}

/*
int main(int nargs,char *args[])
{
  Mesh mesh = NULL;
  parse_mesh("md.bin",&mesh);
  MeshView(mesh);
  return(0);
}
*/
