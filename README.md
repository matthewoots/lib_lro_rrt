# Limited Range Octree RRT Library

### Description
lib_lro_rrt serves a sampling search server where the input is a **pointcloud** (using the PCL library), **current pose** and **destination** and finds a free path in sub millisecond, the library is compiled using CMake with a sample script.

### Setup
```bash
git clone https://github.com/matthewoots/lib_lro_rrt.git
cd lib_lro_rrt
mkdir build && cd build
cmake .. 
make
```

### Run Sample Executable
To run sample scripts go to `build` folder and launch `./lro_rrt_node`, the output in the console prints out until `ctrl-C` or an no valid path is found

Output of `./lro_rrt_node`
```bash
...
Iteration 715
start_position = -8.36553  16.3494  3.13021 end_position = -15.9612 -17.1414  1.42391 distance = 34.3837
local/original size = 4096 / 81920
transformed_start = 0 0 0 transformed_end =     34.3837 1.77636e-15 2.22045e-16
Sub-iterations taken = 1
Successful search complete after 3.083e-05s
Search complete with 3 nodes
Success! Size of 3 found
Duration = 1.33944ms
    node 0 = -8.36553  16.3494  3.13021
    node 1 = 16.3397 13.5121 5.37646
    node 2 = -15.9612 -17.1414  1.42391
...
```

### Include in other projects:
To link this lib properly, add following in the `CMakeLists.txt`
```
find_package(lro_rrt REQUIRED)
target_link_libraries(${PROJECT_NAME}_XX
  ...
  lro_rrt
)
```
