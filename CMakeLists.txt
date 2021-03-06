cmake_minimum_required(VERSION 2.8)

project(WatershedSurface)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(SphereGeneratorTest sphere_generator_test.cc
                                   sphere_generator.cc)
add_executable(SurfaceExtractorTest surface_extractor_test.cc
                                    surface_extractor.cc
                                    surface_smoother.cc
                                    surface_shader.cc
                                    util.cc
                                    lcsGeometry.cpp)
add_executable(SurfaceSmootherTest surface_smoother_test.cc
                                   surface_smoother.cc
                                   util.cc)
add_executable(SurfaceShaderTest surface_shader_test.cc
                                 surface_shader.cc
                                 util.cc
                                 lcsGeometry.cpp)

target_link_libraries(SphereGeneratorTest ${VTK_LIBRARIES})
target_link_libraries(SurfaceExtractorTest ${VTK_LIBRARIES})
target_link_libraries(SurfaceSmootherTest ${VTK_LIBRARIES})
target_link_libraries(SurfaceShaderTest ${VTK_LIBRARIES})
