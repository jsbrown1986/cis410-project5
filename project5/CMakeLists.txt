cmake_minimum_required(VERSION 3.10)

PROJECT(proj5)
SET(VTK_DIR c:/MyProjects/VTK-bin)
SET(CMAKE_VERBOSE_MAKEFILE ON)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE} ${VTK_DIR}/Rendering)

add_executable(proj5 proj5)


target_link_libraries(proj5 glu32)
target_link_libraries(proj5 opengl32)
if(VTK_LIBRARIES)
target_link_libraries(proj5 ${VTK_LIBRARIES})
else()
target_link_libraries(proj5 vtkHybrid)
endif()
