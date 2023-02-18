cmake_minimum_required(VERSION 2.8.3)
set(CMAKE_CXX_STANDARD 14)
project(lro_rrt)

# set default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wno-deprecated-declarations")

set(BUILD_SHARED_LIBS ON)

# file(GLOB INCLUDE_FILES "include/*.h" "kdtree/*.h")
# file(GLOB SRC_FILES "src/lro_rrt_server.cpp" "kdtree/*.c")

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
include_directories(
    ${EIGEN3_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
    include
    kdtree
)

# add_library(${PROJECT_NAME} STATIC ${SRC_FILES} ${INCLUDE_FILES})
add_library(${PROJECT_NAME} 
    src/server.cpp
    src/helper.cpp
    kdtree/kdtree.c
)
target_include_directories (${PROJECT_NAME} PUBLIC 
    ${PROJECT_SOURCE_DIR} 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/kdtree
)