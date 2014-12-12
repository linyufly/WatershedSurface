cmake_minimum_required(VERSION 2.8)

project(WatershedSurface)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(SurfaceExtractor xxx.cc)

target_link_libraries(SurfaceExtractor ${VTK_LIBRARIES})
