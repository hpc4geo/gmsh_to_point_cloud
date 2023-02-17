

import numpy as np
import discretedomain as dd

# Filter line segments (boundary) from the box example
def test_2d_line():
  import meshio as mio

  mesh = mio.read("data/criscros-regions-diri.msh", file_format="gmsh")
  verts = mesh.points
  cells = mesh.cells_dict['line'] # line, triangle, tetra

  bctag = mesh.get_cell_data(name='gmsh:physical', cell_type='line')

  dgeo = dd.Mesh(verts, cells, 'line')
  print(dgeo)
  #dgeo.pyview()

  # Collect cell field data into a dict
  cell_data = {"bctag": [bctag]}

  dgeo.vtu('md.vtu', cdata=cell_data)

  dgeo.write('md.bin')


# Filter triangles (from interior) from the box example
def test_2d_triangle():
  import meshio as mio

  mesh = mio.read("data/criscros-regions-diri.msh", file_format="gmsh")
  verts = mesh.points
  cells = mesh.cells_dict['triangle'] # line, triangle, tetra

  regions = mesh.get_cell_data(name='gmsh:physical', cell_type='triangle')

  dgeo = dd.Mesh(verts, cells, 'triangle')
  print(dgeo)
  #dgeo.pyview()

  # create a new field
  nvert = dgeo.verts.shape[0]

  temp_ic = np.zeros(nvert)

  verts = dgeo.verts
  x = verts[:, 0]
  y = verts[:,1]
  temp_ic[:] = x + y**2

  # Collect vertex field data into a dict
  point_data = {"temperature": temp_ic}

  cells = dgeo.cells
  index = np.arange(cells.shape[0])

  # Collect cell field data into a dict
  cell_data = {"region": [regions], "index": [index]}

  dgeo.vtu('md.vtu', vdata=point_data, cdata=cell_data)

  dgeo.write('md.bin')


# Filter tets (from interior) from the sphere example
def test_3d_tetra():
  import meshio as mio

  mesh = mio.read("data/sphere.msh", file_format="gmsh")
  verts = mesh.points
  cells = mesh.cells_dict['tetra'] # line, triangle, tetra

  regions = mesh.get_cell_data(name='gmsh:physical', cell_type='tetra')

  dgeo = dd.Mesh(verts, cells, 'tetra')
  print(dgeo)

  # create a new field
  nvert = dgeo.verts.shape[0]

  temp_ic = np.zeros(nvert)

  verts = dgeo.verts
  x = verts[:, 0]
  y = verts[:,1]
  temp_ic[:] = x + y**2

  # Collect vertex and cell field data into a dict
  point_data = {"temperature": temp_ic}
  cell_data = {"region": [regions]}

  dgeo.vtu('md.vtu', vdata=point_data, cdata=cell_data)

  dgeo.write('md.bin')

  dgeo.write_fields(cdata=cell_data, vdata=point_data)

# Filter surface triangles (from boundary) from the sphere example
def test_3d_triangle():
  import meshio as mio

  mesh = mio.read("data/sphere.msh", file_format="gmsh")
  verts = mesh.points
  cells = mesh.cells_dict['triangle'] # line, triangle, tetra

  regions = mesh.get_cell_data(name='gmsh:physical', cell_type='triangle')

  dgeo = dd.Mesh(verts, cells, 'triangle')
  print(dgeo)

  # create a new field
  nvert = dgeo.verts.shape[0]

  temp_ic = np.zeros(nvert)

  verts = dgeo.verts
  x = verts[:, 0]
  y = verts[:,1]
  temp_ic[:] = x + y**2

  # Collect vertex and cell field data into a dict
  point_data = {"temperature": temp_ic}
  cell_data = {"region": [regions]}

  #dgeo.vtu('md.vtu') # without any vertex/cell fields
  dgeo.vtu('md.vtu', vdata=point_data, cdata=cell_data)

  dgeo.write('md.bin')


def parse_i3elvis_tetra(filename):
  import meshio as mio

  mesh = mio.read(filename, file_format="gmsh")
  verts = mesh.points
  cells = mesh.cells_dict['tetra'] # line, triangle, tetra

  regions = mesh.get_cell_data(name='gmsh:physical', cell_type='tetra')

  dgeo = dd.Mesh(verts, cells, 'tetra')
  dgeo.build_partitions(npart=200)
  print(dgeo)

  # create a new field
  nvert = dgeo.verts.shape[0]

  temp_ic = np.zeros(nvert)

  verts = dgeo.verts
  x = verts[:, 0]
  y = verts[:,1]
  temp_ic[:] = x + y**2

  # Collect vertex and cell field data into a dict
  point_data = {"temperature": temp_ic}
  cell_data = {"region": [regions]}

  #dgeo.vtu('md.vtu', vdata=point_data, cdata=cell_data)

  dgeo.write('md.bin')

  dgeo.write_fields(cdata=cell_data, vdata=point_data)


if __name__=="__main__":
  #test_2d_line()
  #test_2d_triangle()
  #test_3d_triangle()
  test_3d_tetra()

  #fname = "/Users/dmay/Downloads/I3ELVIS_Gmsh_development/I3ELVIS.msh"
  #parse_i3elvis_tetra(fname)
