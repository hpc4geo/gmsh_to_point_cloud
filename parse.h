
#ifndef __parse_h__
#define __parse_h__

typedef struct _p_CellPartition *CellPartition;
typedef struct _p_Mesh *Mesh;

struct _p_CellPartition {
  int ncell;
  int *cell_list;
  double cmin[3],cmax[3];
};

struct _p_Mesh {
  int nvert,ncell,coor_dim,points_per_cell;
  double *vert;
  int *cell;
  int npartition;
  CellPartition *partition;
};

void CellPartitionCreate(CellPartition *_c);
void CellPartitionDestroy(CellPartition *_c);
void MeshCreate(Mesh *_m);
void MeshDestroy(Mesh *_m);
void MeshView(Mesh m);
void parse_mesh(const char filename[],Mesh *m);
void parse_field(Mesh m,const char filename[],char ftypevoid,void **_data);


#endif
