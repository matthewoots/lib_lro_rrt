# Limited Range Octree RRT Library

### Description
lib_lro_rrt serves a sampling search server where the input is a **pointcloud** (using the PCL library), **current pose** and **destination** and finds a free path in around 1 millisecond, the library is compiled using CMake with a sample script.

### Setup
```bash
git clone https://github.com/matthewoots/lib_lro_rrt.git --recursive
cd lib_lro_rrt
mkdir build && cd build
cmake .. 
make
```
- Used in `PLASTO` package (https://github.com/matthewoots/plasto)

### Performance
Fast 3D RRT* search with some help from ZJU sampling based method
```bash
[lro] iterations(105)
[lro] Successful search complete after 1.15026ms
[lro] intermediate_nodes(78) iterations(105) path_size(2)
```

### Include in other projects:
To link this lib properly, add following in the `CMakeLists.txt`
```txt
add_subdirectory(lib_lro_rrt lro_rrt)

target_link_libraries(XXX
  ...
  lro_rrt
)
```
