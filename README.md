# Workflow

### Preprocessing
* Create GMSH model.
* Mesh it, generate a `.msh` file.
* Parse the `.msh` file within Python using MeshIO.
* Push the mesh representation from MeshIO into a `Mesh` object (`discretedomain.Mesh`)
* Optionally partition the mesh into contiguous chunks - this will facilitate faster point location queries. See `discretedomain.Mesh.build_partitions()`.
* Create cell fields and vertex fields on `Mesh`, see `demo_mesh.py:test_3d_tetra()`.
* Emit mesh data to binary file using `discretedomain.Mesh.write()`.
* Emit field data to binary file using `discretedomain.Mesh.write_fields()`.
* Emit vtu file of `Mesh` and any cell or vertex fields using `Mesh.vtu()`.

### Simulation
* Load binary representation of `Mesh` using `parse.c:parse_mesh()`.
* Load binary representation of cell/vertex fields using `parse.c:parse_field()`.
* Perform point location queries using `point_in_tetra.c:PointLocation_PartitionedBoundingBox()`,



# Requirements

1. MeshIO  (https://pypi.org/project/meshio/1.2.0/)

    `pip install meshio`

   
2. METIS (http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)

    Installation instructions

    * Download this file `http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz`
    * `tar zxvf metis-xxx.tar.gz`
    * `cd metis-5.1.0`
    * `make config shared=1`
        * Note `shared=1` is required to ensure the build creates a shared / dynamic library which is essential for interoperability with Python.
    * `make all`
    * Lastly you will need to set the following environment variable `METIS_DLL` to indicate the full path to the metis shared/dynamic library. For example  
        `export METIS_DLL=${PWD}/build/Darwin-arm64/libmetis/libmetis.dylib`

3. METIS for Python (https://metis.readthedocs.io/en/latest/)

    `pip install metis`
    
4. Compile `parse.c`, `point_in_tetra.c` via the following
    * `gcc -c -O2 -Wall -std=c99 parse.c`
    * `gcc -c -O2 -Wall -std=c99  point_in_tetra.c`
    * A demo / testbed for the C routines are provided in `demo_parse.c`. This can be compiled via the following  
    `gcc -O2 -Wall -std=c99 demo_parse.c -o demo parse.o point_in_tetra.o`
    

# Examples
1. First execute the following `python parse_regions_from_gmsh.py`. This will load a default .msh file and output the mesh and mesh partitions (`md.bin`) and the region data associated with the cells (`region_cell.bin`).
2. Now execute `./demo`. This will load `md.bin` and `region_cell.bin` and perform point locations on particles with randomly generated coordinates.