
import argparse
import numpy as np
import discretedomain as dd


def parse_regions_tetra(filename):
  import meshio as mio

  mesh = mio.read(filename, file_format="gmsh")
  verts = mesh.points
  cells = mesh.cells_dict['tetra'] # line, triangle, tetra

  regions = mesh.get_cell_data(name='gmsh:physical', cell_type='tetra')

  dgeo = dd.Mesh(verts, cells, 'tetra')
  dgeo.build_partitions(npart="auto")
  print(dgeo)

  # Collect vertex and cell field data into a dict
  cell_data = {"region": [regions]}

  dgeo.vtu('md.vtu', cdata=cell_data)

  dgeo.write('md.bin')

  dgeo.write_fields(cdata=cell_data)


if __name__=="__main__":

  fname = "data/I3ELVIS.msh"

  parser = argparse.ArgumentParser(description='Process .msh file.')
  parser.add_argument('-f', '--filename', help='GMSH generated .msh file', default=fname)

  args = parser.parse_args()
  print('Loading file', args.filename)

  parse_regions_tetra(args.filename)
