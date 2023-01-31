
import numpy as np


# independent of cell type
def _adj_v2c(nvert, cells):
  ncell = cells.shape[0]

  v2c = [None]*nvert
  for kv in range(nvert):
    v2c[kv] = list()

  # init
  for kc in range(ncell):
    tri_v = cells[kc, :]
    for tri_vk in tri_v:
      v2c[tri_vk].append(kc)

  # make unique
  for kv in range(nvert):
    u = set(v2c[kv])
    v2c[kv] = list(u)

  return v2c


# independent of cell type
def _adj_v2v(cells, v2c):
  nvert = len(v2c)

  v2v = [None]*nvert
  for kv in range(nvert):
    v2v[kv] = list()

  for kv in range(nvert):
    celln = v2c[kv] # get cell neighbours
    for kc in celln:
      celln_kv = cells[kc, :] # get cells neighbour vertices

      # make sure v doesn't point to itself
      celln_kv = list(celln_kv)
      try:
        celln_kv.remove(kv)
      except:
        pass

      v2v[kv] += celln_kv

  # make unique
  for kv in range(nvert):
    u = set(v2v[kv])
    v2v[kv] = list(u)

  return v2v

def _edges_line(v_):
  point0 = v_[0]
  point1 = v_[1]
  return point0, point1

def _edges_triangle(v_):
  edge0 = set([v_[0], v_[1]])
  edge1 = set([v_[1], v_[2]])
  edge2 = set([v_[0], v_[2]])
  return edge0, edge1, edge2

def _edges_tetra(v_):
  face0 = set([v_[0], v_[1], v_[2]])
  face1 = set([v_[0], v_[1], v_[3]])
  face2 = set([v_[1], v_[2], v_[3]])
  face3 = set([v_[0], v_[2], v_[3]])
  return face0, face1, face2, face3


def _adj_c2c(celltype, cells, v2c):

  edge_selector = {"line":_edges_line , "triangle":_edges_triangle , "tetra":_edges_tetra}

  ncell = cells.shape[0]

  c2c = [None]*ncell
  for kc in range(ncell):
    c2c[kc] = list()

  for kc in range(ncell):
    v_ = cells[kc, :]
    #edge0 = set([v_[0], v_[1]])
    #edge1 = set([v_[1], v_[2]])
    #edge2 = set([v_[0], v_[2]])
    edges_ = edge_selector[celltype](v_)

    for vk in v_:
      neighbour_cells = v2c[vk]
      for nkc in neighbour_cells:
        if kc == nkc: continue # Don't add if the cell neighbour is this cell

        vn_ = cells[nkc, :]
        #nedge0 = set([vn_[0], vn_[1]])
        #nedge1 = set([vn_[1], vn_[2]])
        #nedge2 = set([vn_[0], vn_[2]])
        edgesn_ = edge_selector[celltype](vn_)

        for ee in edges_:

          match = False
          for een in edgesn_:
            if ee == een:
              match = True
              break
          if match == True:
            c2c[kc].append(nkc)

  # make unique
  for kc in range(ncell):
    u = set(c2c[kc])
    c2c[kc] = list(u)

  return c2c


class Mesh:

  def build_cell_adj(self):
    """
    Build and store the adjaceny tables.
    These include
      vertex to cell map (`v2c`)
      vertex to vertex map (`v2v`)
      cell to cell map (`c2c`)

               7 ----- 20
             /  \  B  /  \
            /    \   /    \
           /  A   \ /  C   \
          3 ----- 10 ------ 9
           \  D   /
            \    /
             \  /
              6

    * The vertex to cell map consists of all cells which share a given vertex.
    Example
      10 -> [A, B, C, D]
      3 -> [A, D]
      6 -> [D]
      9 -> [C]
      7 -> [A, B]
      20 -> [B, C]

    * The vertex to vertex map stores neighbour vertices.
    Neighbour vertices are defined as those vertices which are connected
    to any neighbour cell.
    The neighbour of a vertex cannot be itself.
    Example
      3 -> [6, 7, 10]
      6 -> [3, 10]
      7 -> [3, 10, 20]
      9 -> [10, 20]
      10 -> [3, 6, 7, 9, 10]
      20 -> [7, 9, 10]

    * The cell to cell map defines neighbour cells based on whether a cell
    shares an edge.
    In 1D (`line`) the edge consists of a single point.
    In 2D (`triangle`) the edge consists of a line segment consisting of two points.
    In 3D (`tetra`) the edge consists of a triangle defined by three points.
    Example
      A -> [B, D]
      B -> [A, C]
      C -> [B]
      D -> [A]
    """
    nvert = self.verts.shape[0]
    ncell = self.cells.shape[0]
    self.v2c = _adj_v2c(nvert, self.cells)
    self.v2v = _adj_v2v(self.cells, self.v2c)
    self.c2c = _adj_c2c(self.celltype, self.cells, self.v2c)


  def flatten_vertices(self):
    """
    Remove any unused vertices (i.e. vertices in the file but not referenced by the cells).
    Renumber the cells[] np.ndarray.
    Over-rides the values of self.verts.
    Stores the vertices selected (`self.selected`) and the global-to-local mapping (`self.mapping`).
    """
    nvert0 = self.verts.shape[0]
    ncell = self.cells.shape[0]
    tag = np.array( [False]*nvert0 )
    print('nvert-init', nvert0)
    for kc in range(ncell):
      v_ = self.cells[kc, :]
      tag[v_] = True

    selected_ = np.where(tag == True)
    self.selected = np.array(selected_[0], dtype=np.int64)
    nvert1 = self.selected.shape[0]
    print('nvert-final', nvert1)

    mapping = np.zeros(nvert0, dtype=np.int64)
    mapping[:] = -1

    # non-vectorized
    #for vk in range(nvert1):
    #  mapping[ self.selected[vk] ] = vk
    #
    # vectorized
    vk_array = np.arange(nvert1)
    mapping[ self.selected[:] ] = vk_array[:]

    # apply mapping to cells
    #print(self.cells[0,:])
    cells_ = self.cells.view()
    cells_.shape = ncell * self.points_per_entity
    # non-vectorized
    #for vk in range(cells_.shape[0]):
    #  cells_[vk] = mapping[ cells_[vk] ]
    # vectorized
    cells_[:] = mapping[ cells_[:] ]
    #print(self.cells[0,:])

    min_vidx = np.min(cells_)
    if min_vidx < 0:
      raise RuntimeError('re-mapping failed [min]')
    max_vidx = np.max(cells_)
    if max_vidx >= nvert1:
      raise RuntimeError('re-mapping failed [max]')

    self.mapping = mapping
    self.verts = self.verts[tag]


  def __init__(self, verts, cells, celltype):

    self.verts = verts

    self.coor_dim = verts.shape[1]

    if celltype == 'line':
      self.points_per_entity = 2
      self.cell_adj_max = 2
    elif celltype == 'triangle':
      self.points_per_entity = 3
      self.cell_adj_max = 3
    elif celltype == 'tetra':
      self.points_per_entity = 4
      self.cell_adj_max = 4
    else:
      raise RuntimeError('Cell type unsupported')

    self.celltype = celltype
    if cells.shape[1] != self.points_per_entity:
      raise RuntimeError('Points per cell is inconsistent')

    self.cells = cells

    self.flatten_vertices()

    self.build_cell_adj()
    self.partition = None


  def build_partitions(self, npart='auto'):
    """
    Partition mesh into `npart` subdomains.

    Input
    -----
      - npart: int (default = `'auto'`)
        Number of contiguous sub-domains the mesh will be partitioned into.
    """
    import metis as metis

    if npart == 'auto':
      ntarget_cells = 500.0
      npart = float(self.cells.shape[0]) / ntarget_cells
      npart = int(npart)
      if npart < 1: npart = 1

    graph = metis.adjlist_to_metis(self.c2c)
    partition = metis.part_graph(graph, nparts=npart, tpwgts=None, ubvec=None, recursive=False,
                     objtype='cut',
                     numbering=0,
                     contig=True,
                     seed=0)

    #print(partition)
    part_idx = np.array(partition[1], dtype=np.int32)

    self.partition = [None] * npart
    for p in range(npart):
      self.partition[p] = list()
    for k in range(part_idx.shape[0]):
      self.partition[ part_idx[k] ].append(k)
    # flatten and turn into numpy array
    for p in range(npart):
      cells = self.partition[p]
      self.partition[p] = np.array(cells, dtype=np.int32)
      #print(self.partition[p])


  def pyview(self):
    """
    Display line mesh or triangle mesh to the screnn.
    """
    import matplotlib.pyplot as plt

    ax = plt.figure(figsize=(8, 8))

    if self.celltype == 'line':
      for kc in range(self.cells.shape[0]):
        v_ = self.cells[kc, :]
        x_ = self.verts[v_, 0]
        y_ = self.verts[v_, 1]
        plt.plot(x_, y_, 'bo-')

    elif self.celltype == 'triangle':
      import matplotlib.tri as mtri
      coor2 = self.verts[:, 0:2]
      triang = mtri.Triangulation(coor2[:,0], coor2[:,1], self.cells)
      plt.triplot(triang, 'bo-')

    elif self.celltype == 'tetra':
      raise RuntimeError('Unsupported for 3d tetra')

    plt.show()
    #plt.savefig('mesh.pdf')
    plt.close()


  def __str__(self):
    m = 2
    msg = 'Mesh' + '\n'
    msg += ''*m + 'npoints:   ' + str(self.verts.shape[0])  + '\n'
    msg += ''*m + 'coord-dim: ' + str(self.coor_dim)  + '\n'
    msg += ''*m + 'cell type: ' + self.celltype  + '\n'
    msg += ''*m + 'ncells:    ' + str(self.cells.shape[0])
    return msg


  def vtu(self, fname, vdata=None, cdata=None):
    """
    Generate a VTU file of the mesh.

    Input
    -----
      - fname: str
        Name of file. Must end in .vtu
      - vdata: dict (optional)
        Vertex data. Dict key,value pairs define the field names and the np.ndarray storing the field.
        Example
          `
          temp_ic = np.zeros(mesh.verts.shape[0]
          point_data = {"temperature": temp_ic}
          mesh.vtu("demo.vtu", vdata=point_data)
          `
      - cdata: dict (optional)
        Cell data. Dict key,value pairs define the field names and a list of np.ndarray's for storing the field.
        Note the difference with `vdata` in that cell data fields need to be provided as an iterable.
        Example
          `
          region = np.zeros(mesh.cells.shape[0])
          cell_data = {"region": [regions]}
          mesh.vtu("demo.vtu", cdata=cell_data)
          `
    """
    import meshio as mio

    cells = {self.celltype: self.cells} #cells = {"triangle": tri}

    mesh = mio.Mesh( self.verts, cells, point_data=vdata)

    #if vdata is not None:
    #  mesh.point_data = vdata

    if cdata is not None:
      #mesh.cell_data = {"region": [regions], "partition": [part_idx]}
      mesh.cell_data = cdata

    if fname[-4:] != '.vtu':
      raise RuntimeError('Filename must end in .vtu')
    mesh.write(fname, file_format="vtu")


  def write(self, fname):
    """
    Write mesh to a binary file. The endianess used matches that of the machine
    writing the data - hence this is output file is not platform portable.

    Input
    -----
      - fname: str
        Name of output file.
    """
    import sys, struct

    endian = sys.byteorder
    print('write() -> endianness', endian)
    if endian == 'little':
      w_int = '<i' # i-int, l-long, L-unsigned long
    else:
      w_int = '>i'

    with open(fname, "wb") as fp:
      byt = struct.pack(w_int, np.int32(self.verts.shape[0])) # nverts
      fp.write(byt)
      byt = struct.pack(w_int, np.int32(self.verts.shape[1])) # coordinate_dimension
      fp.write(byt)
      self.verts.tofile(fp) # vertices
      print('  verts.dtype', self.verts.dtype)

      byt = struct.pack(w_int, np.int32(self.cells.shape[0])) # ncells
      fp.write(byt)
      byt = struct.pack(w_int, np.int32(self.cells.shape[1])) # points-per-cell
      fp.write(byt)
      cells_ = np.int32(self.cells)
      cells_.tofile(fp) # cells
      print('  cells.dtype', cells_.dtype)

      if self.partition is None:
        print('  -> no partitions found... generating default single partition')
        npart = 1
        byt = struct.pack(w_int, np.int32(npart)) # npartitions
        fp.write(byt)

        bb = np.min(self.verts, axis=0)
        bb.tofile(fp) # bbmin
        print('  bbox_min.dtype', bb.dtype)

        bb = np.max(self.verts, axis=0)
        bb.tofile(fp) # bbmax
        print('  bbox_max.dtype', bb.dtype)

        cp = np.arange(self.cells.shape[0], dtype=np.int32)
        byt = struct.pack(w_int, np.int32(cp.shape[0])) # ncells
        fp.write(byt)
        cp.tofile(fp)
        print('  part->celllist.dtype', cp.dtype) # cells
      else:
        npart = len(self.partition)
        print('  -> ' + str(npart) + ' partitions found')

        byt = struct.pack(w_int, np.int32(npart)) # npartitions
        fp.write(byt)

        for p in range(npart):
          cp = self.partition[p]
          cell_subset_ = self.cells[cp, :]
          cell_subset_v_ = cell_subset_.view()
          cell_subset_v_.shape = cell_subset_.shape[0] * cell_subset_.shape[1]

          verts_ = self.verts[cell_subset_v_, :]

          bb = np.min(verts_, axis=0)
          bb.tofile(fp)
          if p == 0: print('  bbox_min.dtype', bb.dtype)

          bb = np.max(verts_, axis=0)
          bb.tofile(fp)
          if p == 0: print('  bbox_max.dtype', bb.dtype)

          byt = struct.pack(w_int, np.int32(cp.shape[0])) # little endian
          fp.write(byt)

          cp.tofile(fp) # cell
          if p == 0: print('  part->celllist.dtype', cp.dtype)


  def write_fields(self, fname_prefix="./", cdata=None, vdata=None):
    """
    Write fields to a binary file.
    The order of the fields matches the order appearing in each dict.
    The endianess used matches that of the machine
    writing the data - hence this is output file is not platform portable.

    Input
    -----
      - fname: str
        Name of output file.
      - cell_data: dict
        Cell data defined on the mesh.
        The key/value pairs correspond to the fieldname and the data array.
      - vertex_data: dict
        Vertex data defined on the mesh.
        The key/value pairs correspond to the fieldname and the data array.
    """
    import sys, struct
    from collections import OrderedDict

    endian = sys.byteorder
    print('write_fields() -> endianness', endian)
    if endian == 'little':
      w_int = '<i' # i-int, l-long, L-unsigned long
    else:
      w_int = '>i'

    with open("numpy_dtypes.txt", "w") as fp:
      fp.write("np.short | np.int16 -> short -> 10" + "\n")
      fp.write("np.intc | np.int32 -> int -> 11" + "\n")
      fp.write("np.int_ | np.int64 -> long int -> 12" + "\n")

      fp.write("np.ushort | np.uint16 -> unsigned short -> 13" + "\n")
      fp.write("np.uintc | np.uint32 -> unsigned int -> 14" + "\n")
      fp.write("np.uint_ | np.uint64 -> unsigned long int -> 15" + "\n")

      fp.write("np.single | np.float32 -> float -> 20" + "\n")
      fp.write("np.single | np.float64 -> double -> 21" + "\n")
      fp.write("np.longdouble | np.float128 -> __float128 -> 22" + "\n")

    dtype_map = dict()
    dtype_map[np.short] = 10
    dtype_map[np.int16] = 10

    dtype_map[np.intc] = 11
    dtype_map[np.int32] = 11

    dtype_map[np.int_] = 12
    dtype_map[np.int64] = 12

    dtype_map[np.single] = 20
    dtype_map[np.float32] = 20

    dtype_map[np.double] = 21
    dtype_map[np.float64] = 21

    if cdata is not None:
      data_ = OrderedDict(cdata)
      keys_ = list(data_.keys())

      for k in keys_:
        fname = "".join([fname_prefix, str(k) , "_cell.bin"])
        print('  writing file', fname)
        with open(fname, "wb") as fp:
          byt = struct.pack(w_int, np.int32(self.cells.shape[0])) # ncells
          fp.write(byt)

          field_ = data_[k][0]
          byt = struct.pack(w_int, np.int32( dtype_map[field_.dtype.type] )) # data-type-key
          fp.write(byt)

          field_.tofile(fp) # cell field
          print('  cell-field_['+k+'].dtype', field_.dtype)


    if vdata is not None:
      data_ = OrderedDict(vdata)
      keys_ = list(data_.keys())

      for k in keys_:
        fname = "".join([fname_prefix, str(k) , "_vertex.bin"])
        print('  writing file', fname)
        with open(fname, "wb") as fp:
          byt = struct.pack(w_int, np.int32(self.verts.shape[0])) # nverts
          fp.write(byt)

          field_ = data_[k]
          byt = struct.pack(w_int, np.int32( dtype_map[field_.dtype.type] )) # data-type-key
          fp.write(byt)

          field_.tofile(fp)
          print('  vertex-field_['+k+'].dtype', field_.dtype)
